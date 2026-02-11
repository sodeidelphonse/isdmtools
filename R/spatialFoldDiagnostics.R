
#' Diagnostic for Spatial Folds Geometry
#'
#' @description
#' Internal function to evaluate the spatial structure of folds for blocked cross-validation.
#'
#' @param data_all An sf object containing pooled locations and fold IDs.
#' @param fold_col Character. Name of the fold ID column.
#' @param rho Numeric. Estimated spatial range (km).
#' @param plot Logical. If TRUE, generates a ggplot object.
#'
#' @importFrom units set_units
#' @keywords internal
#'
check_spatial_geometry <- function(data_all, fold_col = "folds_ids", rho = NULL, plot = TRUE) {

  data_all[[fold_col]] <- as.factor(data_all[[fold_col]])
  folds_ids <- levels(data_all[[fold_col]])

  # Folds centroids
  centroids <- data_all %>%
    dplyr::group_by(.data[[fold_col]]) %>%
    dplyr::summarise(geometry = sf::st_union(.data$geometry), .groups = "drop") %>%
    sf::st_centroid()

  coords_mat <- sf::st_coordinates(centroids)
  centroids_data <- sf::st_drop_geometry(centroids) %>%
    dplyr::mutate(x = coords_mat[, 1], y = coords_mat[, 2])

  # Internal distances (point-to-centroid)
  points_with_centroids <- data_all %>%
    dplyr::left_join(
      as.data.frame(centroids) %>%
        dplyr::select(!!rlang::sym(fold_col), centroid_geom = .data$geometry),
      by = fold_col
    )

  dist_raw <- sf::st_distance(
      sf::st_geometry(points_with_centroids),
      sf::st_geometry(points_with_centroids$centroid_geom),
      by_element = TRUE
    )
  points_with_centroids$dist_val <- as.numeric(units::set_units(dist_raw, "km"))

  # Inter-block gap
  min_gaps <- vapply(folds_ids, function(id) {
    block_pts <- data_all[data_all$folds_ids == id, ]
    other_pts <- data_all[data_all$folds_ids != id, ]
    if(nrow(other_pts) == 0) return(NA_real_)
    d_m <- min(sf::st_distance(block_pts, other_pts))
    as.numeric(units::set_units(d_m, "km"))
  }, FUN.VALUE = numeric(1))

  gap_df <- data.frame(
    folds_ids = names(min_gaps),
    min_gap_km = as.numeric(min_gaps),
    stringsAsFactors = FALSE
  )
  gap_df[[fold_col]] <- factor(gap_df[[fold_col]], levels = folds_ids)

  # Merging summary stats
  summary_stats <- points_with_centroids %>%
    sf::st_drop_geometry() %>%
    dplyr::group_by(.data[[fold_col]]) %>%
    dplyr::summarise(
      n_points = dplyr::n(),
      max_dist_km = max(.data$dist_val, na.rm = TRUE),
      mean_dist_km = mean(.data$dist_val, na.rm = TRUE),
      median_dist_km = stats::median(.data$dist_val, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::left_join(gap_df, by = fold_col) %>%
    dplyr::left_join(centroids_data, by = fold_col)

  # Independence logic
  if (!is.null(rho)) {
    summary_stats <- summary_stats %>%
      dplyr::mutate(
        gap_rho_ratio = .data$min_gap_km / rho,
        independence = dplyr::case_when(
          .data$min_gap_km == 0 ~ "Contiguous",
          .data$gap_rho_ratio < 1 ~ "Weakly Independent",
          .data$gap_rho_ratio >= 1 & .data$gap_rho_ratio < 2 ~ "Independent",
          TRUE ~ "Strongly Independent"
        )
      )
  }

  diag_plot <- NULL
  if (plot) {
    buffer_geoms <- sf::st_buffer(
      sf::st_geometry(centroids),
      dist = units::set_units(summary_stats$max_dist_km, "km")
      )

    buffers <- centroids %>%
      dplyr::mutate(geometry = buffer_geoms) %>%
      dplyr::left_join(as.data.frame(summary_stats), by = fold_col)

    rho_val <- if (is.null(rho)) "N/A" else rho
    sub_text <- bquote(paste("Estimated range (", rho, ") = ", .(rho_val), " km"))

    diag_plot <- ggplot2::ggplot() +
      ggplot2::geom_sf(data = buffers, ggplot2::aes(fill = .data[[fold_col]]),
                       alpha = 0.1, linetype = "dashed", color = "black") +
      ggplot2::geom_sf(data = data_all, ggplot2::aes(color = .data[[fold_col]]), size = 1) +
      ggplot2::geom_sf(data = centroids, ggplot2::aes(shape = "Fold Centroid"),
                       color = "black", size = 3, stroke = 1) +
      ggplot2::scale_shape_manual(name = "Reference", values = c("Fold Centroid" = 3)) +
      ggplot2::theme_minimal() +
      ggplot2::labs(title = "Spatial Fold Partitioning",
                    subtitle = sub_text,
                    fill = "Folds", color = "Folds"
                    )
  }

  return(list(summary = summary_stats, plot = diag_plot, rho = rho))
}


#' Check Spatial Folds Independence
#'
#' @description
#' Evaluates the geometric properties of spatial folds to ensure spatial independence
#' for cross-validation. This function verifies if block sizes and inter-block gaps
#' are sufficient relative to the prior or model's estimated spatial range (\eqn{\rho}).
#'
#' @param object A \code{DataFolds} object created by \code{create_folds()}.
#' @param rho Numeric. Optional. The spatial range (km) estimated from the exploratory
#' analysis (e.g., \link[blockCV]{cv_spatial_autocor}) and used as the block size or
#' the one estimated from the integrated model (e.g., the Matérn range parameter).
#' @param plot Logical. If \code{TRUE}, returns a diagnostic plot.
#' @param ... Additional arguments.
#'
#' @details
#' The function assesses independence based on the minimum gap between folds
#' compared to the spatial range (\eqn{\rho}):
#' \itemize{
#'   \item \strong{Contiguous}: Gap = 0. High risk of spatial leakage; observations in
#'   test folds are spatially correlated with training data.
#'   \item \strong{Weakly Independent}: 0 < Gap < \eqn{\rho}. A physical gap exists,
#'   but correlation remains above 0.1.
#'   \item \strong{Independent}: \eqn{\rho \le} Gap < \eqn{2\rho}. Spatial correlation
#'   is below 0.1 at the boundary; considered robust for most CV applications.
#'   \item \strong{Strongly Independent}: Gap \eqn{\ge 2\rho}. Spatial correlation
#'   is effectively zero, providing the most rigorous test of model extrapolation.
#' }
#'
#' @return An object of class \code{GeoDiagnostic}.
#' @export
#' @family blocks diagnostics
#'
#' @references
#' \itemize{
#'   \item Roberts DR, Bahn V, Ciuti S, Boyce MS, Elith J, Guillera-Arroita G, Hauenstein S, Lahoz-Monfort JJ, Schröder B, Thuiller W, et al. Cross-validation strategies for data with temporal, spatial, hierarchical, or phylogenetic structure. _Ecography_ (2017) 40:913–929. \doi{10.1111/ecog.02881}
#' }
#'
#' @examples
#' \dontrun{
#' library(sf)
#' library(terra)
#' library(ggplot2)
#' library(isdmtools)
#'
#' # Generate the data as a list of sf objects
#' set.seed(42)
#' presence_data <- data.frame(
#'  x = runif(100, 0, 4),
#'  y = runif(100, 6, 13),
#'  site = rbinom(100, 1, 0.6)
#' ) %>% st_as_sf(coords = c("x", "y"), crs = 4326)
#'
#' count_data <- data.frame(
#'   x = runif(50, 0, 4),
#'  y = runif(50, 6, 13),
#'  count = rpois(50, 5)
#' ) %>% st_as_sf(coords = c("x", "y"), crs = 4326)
#'
#' datasets_list <- list(Presence = presence_data, Count = count_data)
#'
#' ben_coords <- matrix(c(0, 6, 4, 6, 4, 13, 0, 13, 0, 6), ncol = 2, byrow = TRUE)
#' ben_sf <- st_sfc(st_polygon(list(ben_coords)), crs = 4326)
#' ben_sf <- st_sf(data.frame(name = "Benin"), ben_sf)
#'
#' # Create Folds using create_folds()
#' folds <- create_folds(
#'   datasets_list,
#'   region_polygon = ben_sf,
#'   k = 5,
#'   cv_method = "cluster"
#' )
#'
#' # Check Spatial Independence
#' # Assuming autocorrelation range (rho) is 150 km
#' spat_diag <- check_folds(folds, rho = 150, plot = TRUE)
#'
#' # View results
#' print(spat_diag)
#' plot(spat_diag)
#' }
check_folds <- function(object, ...) {
  UseMethod("check_folds")
}

#' @rdname check_folds
#' @method check_folds DataFolds
#' @export
check_folds.DataFolds <- function(object, rho = NULL, plot = TRUE, ...) {
  res <- check_spatial_geometry(
    data_all = object$data_all,
    fold_col = "folds_ids",
    rho = rho,
    plot = plot
  )
  class(res) <- "GeoDiagnostic"

  return(res)
}


#' Methods for GeoDiagnostic objects
#'
#' @param x A \code{GeoDiagnostic} object.
#' @param ... Additional arguments.
#'
#' @return
#' \itemize{
#'   \item \code{print}: Invisibly returns the original object.
#'   \item \code{plot}: Returns a \code{ggplot2} object.
#' }
#'
#' @export
#' @family blocks diagnostics
print.GeoDiagnostic <- function(x, ...) {
  cat("\n=== isdmtools: Spatial Fold Diagnostic ===\n")

  if (!is.null(x$rho)) {
    cat(paste0("Model Spatial Range (rho): ", x$rho, " km\n\n"))
    if ("independence" %in% names(x$summary)) {
      status_summary <- x$summary %>%
        dplyr::group_by(.data$independence) %>%
        dplyr::summarise(Count = dplyr::n(), .groups = "drop")
      print(as.data.frame(status_summary))
    }
  }

  cat("\nInternal Size (Max Distance to Fold Centroid):\n")
  print(summary(x$summary$max_dist_km))

  cat("\nInter-block Gap (Min Distance to Nearest Fold):\n")
  print(summary(x$summary$min_gap_km))

  if (any(x$summary$min_gap_km == 0, na.rm = TRUE)) {
    cat("\nWARNING: One or more folds are contiguous (Gap = 0). Potential spatial leakage.\n")
  }

  cat("==========================================\n")
  invisible(x)
}

#' @rdname print.GeoDiagnostic
#' @export
plot.GeoDiagnostic <- function(x, ...) {
  if (is.null(x$plot)) stop("No plot found.")

  rho_val <- if (is.null(x$rho)) "N/A" else x$rho

  p <- x$plot +
    ggplot2::labs(
      title = "Spatial Diagnostics: Size & Isolation",
      subtitle = bquote(paste("Estimated range (", rho, ") = ", .(rho_val), " km"))
    )
  return(p)
}


#--- Diagnosis in the environmental space -----

#' Check Environmental Balance of Folds
#'
#' @description Evaluates whether environmental covariates are well-balanced
#' across folds. It extracts values from a \code{SpatRaster} at point locations.
#' For large rasters, it uses a background sample to represent the study area.
#'
#' @param object A \code{DataFolds} object.
#' @param covariates A \code{SpatRaster} (terra) containing environmental layers.
#' It must have the same coordinate reference system (CRS) as the sf objects used for blocking.
#' @param plot_type Character. Either "density" (default) or "boxplot".
#' @param n_background Numeric. Number of background points to sample for environmental
#' space representation. Default 10,000.
#' @param ... Additional arguments passed to \code{sample_background}.
#'
#' @details
#' The function also calculates the environmental niche overlap using Schoener's D
#' metric (Schoener, 1968). This metric ranges from 0 (no overlap) to 1 (identical
#' niches).
#' \itemize{
#'   \item \strong{Schoener's D}: Quantifies how well each fold represents the
#'   available environmental space (the background). A median value is reported
#'   across all folds for each covariate.
#'   \item \strong{Interpretation}: Values > 0.6 generally indicate that the folds
#'   are representative of the study area's environmental conditions. Low values
#'   suggest that cross-validation results may be biased because the model is being
#'   tested on environmental conditions it rarely encountered during training.
#' }
#'
#' @return An object of class \code{EnvDiagnostic}.
#' @export
#' @family blocks diagnostics
#'
#' @examples
#' \dontrun{
#' library(sf)
#' library(terra)
#' library(ggplot2)
#' library(isdmtools)
#'
#' # Generate data as a list of sf objects
#' set.seed(42)
#' presence_data <- data.frame(
#'  x = runif(100, 0, 4),
#'  y = runif(100, 6, 13),
#'  site = rbinom(100, 1, 0.6)
#' ) %>% st_as_sf(coords = c("x", "y"), crs = 4326)
#'
#' count_data <- data.frame(
#'   x = runif(50, 0, 4),
#'  y = runif(50, 6, 13),
#'  count = rpois(50, 5)
#' ) %>% st_as_sf(coords = c("x", "y"), crs = 4326)
#'
#' datasets_list <- list(Presence = presence_data, Count = count_data)
#'
#' ben_coords <- matrix(c(0, 6, 4, 6, 4, 13, 0, 13, 0, 6), ncol = 2, byrow = TRUE)
#' ben_sf <- st_sfc(st_polygon(list(ben_coords)), crs = 4326)
#' ben_sf <- st_sf(data.frame(name = "Benin"), ben_sf)
#'
#' # a) Continuous covariates
#' r   <- rast(ben_sf, nrow = 100, ncol = 100, crs = 'epsg:4326')
#' r[] <- rnorm(ncell(r))
#' rtmp   <- r
#' rtmp[] <- runif(ncell(r), 5, 10)
#'
#' r_stk <- c(r, rtmp + r)
#' names(r_stk) <- c("cov1", "cov2")
#'
#' # Create Folds
#' folds <- create_folds(
#'   datasets_list,
#'   region_polygon = ben_sf,
#'   cv_method = "cluster"
#' )
#'
#' # Check Environmental Representation
#' env_diag <- suppressWarnings(check_env_balance(
#'   folds,
#'   covariates = r_stk,
#'   n_background = 5000)
#' )
#'
#' # View p-values in console
#' print(env_diag)
#'
#' # View density plots
#' plot(env_diag)
#'
#' # b) Mixture of continuous and categorical covariates
#' r_temp <- rast(extent = c(0, 4, 6, 13), res = 0.1, val = runif(2500, 15, 25))
#' r_land <- rast(extent = c(0, 4, 6, 13), res = 0.1, val = sample(1:3, 2500, TRUE))
#'
#' # Set up land cover as a factor
#' levels(r_land) <- data.frame(ID = 1:3, cover = c("Forest", "Grass", "Urban"))
#' env_stack <- c(r_temp, r_land)
#' names(env_stack) <- c("temperature", "land_use")
#'
#' # Run the diagnostic with 2000 cells over 2800 available
#' env_diag <- check_env_balance(
#'   folds,
#'   covariates = env_stack,
#'   n_background = 2000,
#'   plot_type = "boxplot"
#' )
#'
#' # Inspect results
#' # 'temperature' will show a p-value and Schoener's D
#' # 'land_use' will show a p-value (Chi-sq) and Schoener_D as NA
#' print(env_diag)
#' plot(env_diag)
#' }
check_env_balance <- function(object, ...) {
  UseMethod("check_env_balance")
}

#' @rdname check_env_balance
#' @method check_env_balance DataFolds
#' @export
check_env_balance.DataFolds <- function(object, covariates, plot_type = c("density", "boxplot"),
                                        n_background = 10000, ...) {

  plot_type <-match.arg(plot_type)

  # Prepare point data
  data_sf <- object$data_all %>%
    dplyr::filter(!is.na(.data$folds_ids))

  if (!inherits(covariates, "SpatRaster")) stop("covariates must be a SpatRaster.")
  if (!terra::same.crs(data_sf, covariates)) stop("CRS mismatch between points and raster.")

  fold_vals <- terra::extract(covariates, sf::st_coordinates(data_sf))
  cov_names <- names(fold_vals)

  raw_data <- data.frame(
    folds_ids = factor(data_sf$folds_ids),
    fold_vals,
    stringsAsFactors = FALSE
  )

  # Sample background points
  back_pts  <- sample_background(mask = covariates, n = n_background, values = FALSE, ...)
  back_vals <- terra::extract(covariates, back_pts$bg[, c("x", "y")])

  back_data <- data.frame(
    folds_ids = factor("Background"),
    back_vals,
    stringsAsFactors = FALSE
  )

  full_env_data <- dplyr::bind_rows(raw_data, back_data)

  # Statistical testing and overlap calculation
  stats_list <- lapply(cov_names, function(v) {
    if (is.numeric(raw_data[[v]])) {

      # Internal balance: are folds statistically different from each other?
      test <- stats::kruskal.test(raw_data[[v]] ~ raw_data$folds_ids)
      type <- "Continuous"

      # External balance: niche overlap between folds and background
      d_vals <- vapply(split(raw_data[[v]], raw_data$folds_ids),
                       function(x) calc_niche_overlap(x, back_data[[v]]), FUN.VALUE = numeric(1))
      overlap_val <- round(stats::median(d_vals, na.rm = TRUE), 3)
    } else {
      test <- stats::chisq.test(table(raw_data[[v]], raw_data$folds_ids))
      type <- "Categorical"
      overlap_val <- NA_real_
    }

    data.frame(
      Variable = v,
      Type = type,
      p_val = round(test$p.value, 4),
      Schoener_D = overlap_val,
      label = if (is.na(overlap_val)) {
        sprintf("%s\n(p = %.3f)", v, test$p.value)
      } else {
        sprintf("%s\n(p = %.3f, D = %.3f)", v, test$p.value, overlap_val)
      },
      stringsAsFactors = FALSE
    )
  })
  stats_df <- dplyr::bind_rows(stats_list)

  plot_list <- lapply(cov_names, function(v) {
    data.frame(
      Fold = full_env_data$folds_ids,
      Variable = v,
      Value = as.numeric(full_env_data[[v]]),
      stringsAsFactors = FALSE
    )
  })
  plot_df <- dplyr::bind_rows(plot_list)

  plot_df$Variable <- factor(plot_df$Variable,
                             levels = stats_df$Variable,
                             labels = stats_df$label)

  n_folds <- length(levels(raw_data$folds_ids))
  fold_colors <- .get_isdm_palette(n_folds)
  names(fold_colors) <- levels(raw_data$folds_ids)

  all_colors <- c("Background" = "#228B22", fold_colors)

  p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = .data$Value, fill = .data$Fold)) +
    ggplot2::facet_wrap(~ Variable, scales = "free") +
    ggplot2::scale_fill_manual(values = all_colors) +
    ggplot2::theme_bw() +
    ggplot2::labs(
      fill = "Group",
      title = "Environmental Diagnostics"
    )

  if (plot_type == "density") {
    p <- p + ggplot2::geom_density(alpha = 0.4)
  } else {
    p <- p + ggplot2::geom_boxplot(ggplot2::aes(y = .data$Fold))
  }

  res <- list(summary = stats_df, plot = p)
  class(res) <- "EnvDiagnostic"

  return(res)
}


#' Methods for EnvDiagnostic objects
#'
#' @param x An \code{EnvDiagnostic} object.
#' @param ... Additional arguments.
#'
#' @return
#' \itemize{
#'   \item \code{print}: Invisibly returns the original object.
#'   \item \code{plot}: Returns a \code{ggplot2} object for covariates' density plots.
#' }
#'
#' @export
#' @family blocks diagnostics
print.EnvDiagnostic <- function(x, ...) {
  cat("\n=== isdmtools: Environmental Balance Diagnostic ===\n")
  cat("Significance (p > 0.05 = Balanced)\n")
  cat("Overlap (D > 0.6 = Representative)\n\n")

  print(as.data.frame(x$summary[, c("Variable", "Type", "p_val", "Schoener_D")]))

  low_overlap <- x$summary$Variable[x$summary$Schoener_D < 0.5]
  if (length(low_overlap) > 0) {
    cat("\nWARNING: Poor environmental representation (D < 0.5) in:",
        paste(low_overlap, collapse = ", "), "\n")
  }

  cat("==================================================\n")
  invisible(x)
}

#' @rdname print.EnvDiagnostic
#' @export
plot.EnvDiagnostic <- function(x, ...) {
  if (is.null(x$plot)) {
    stop("No plot was found in the EnvDiagnostic object.")
  }
  return(x$plot)
}


#' Summarise Fold Diagnostics
#'
#' @description Combines geographic and environmental diagnostics into a
#' single unified report to evaluate the quality of a cross-validation scheme.
#'
#' @param geo_diag A \code{GeoDiagnostic} object.
#' @param env_diag An \code{EnvDiagnostic} object.
#'
#' @return An object of class \code{FoldsSummary}, which inherits from \code{data.frame}.
#' @export
#' @family blocks diagnostics
#'
#' @examples
#' \dontrun{
#' library(sf)
#' library(terra)
#' library(ggplot2)
#' library(isdmtools)
#'
#' # Generate points data
#' set.seed(42)
#' presence_data <- data.frame(
#'   x = runif(100, 0, 4),
#'   y = runif(100, 6, 13),
#'   site = rbinom(100, 1, 0.6)
#' ) %>% st_as_sf(coords = c("x", "y"), crs = 4326)
#'
#' count_data <- data.frame(
#'   x = runif(50, 0, 4),
#'   y = runif(50, 6, 13),
#'  count = rpois(50, 5)
#' ) %>% st_as_sf(coords = c("x", "y"), crs = 4326)
#'
#' datasets_list <- list(Presence = presence_data, Count = count_data)
#'
#' # Environmental data
#' set.seed(42)
#' r <- rast(extent = c(0, 4, 6, 13), nrow=100, ncol=100, crs='epsg:4326')
#' r[] <- rnorm(ncell(r))
#' rtmp   <- r
#' rtmp[] <- runif(ncell(r), 5, 10)
#'
#' r_stk <- c(r, rtmp + r)
#' names(r_stk) <- c("cov1", "cov2")
#'
#' # Create Folds
#' folds <- create_folds(datasets_list, cv_method = "cluster")
#'
#' # Spatial diagnostics
#' spat_diag <- check_folds(folds, plot = TRUE)
#'
#' # Environmental diagnostics
#' env_diag <- suppressWarnings(check_env_balance(
#'   folds,
#'   covariates = r_stk,
#'   n_background = 5000)
#' )
#'
#' # Combined diagnostics
#' summarise_fold_diagnostics(spat_diag, env_diag)
#' }
summarise_fold_diagnostics <- function(geo_diag, env_diag) {

  if (!inherits(geo_diag, "GeoDiagnostic") || !inherits(env_diag, "EnvDiagnostic")) {
    stop("Inputs must be 'GeoDiagnostic' and 'EnvDiagnostic' objects.")
  }

  # Captures spatial independence and fold size
  geo_df <- data.frame(
    Domain = "Geographic",
    Metric = c("Avg Internal Distance (km)", "Avg Inter-Fold Gap (km)"),
    Value = c(mean(geo_diag$summary$max_dist_km, na.rm = TRUE),
              mean(geo_diag$summary$min_gap_km, na.rm = TRUE)),
    Status = ifelse(any(geo_diag$summary$min_gap_km == 0), "Contiguous", "Separated"),
    stringsAsFactors = FALSE
  )

  # Captures niche representation and internal balance across folds
  env_summary <- env_diag$summary
  env_df <- data.frame(
    Domain = "Environmental",
    Metric = c("Median Overlap (D)", "Minimum p-value"),
    Value = c(stats::median(env_summary$Schoener_D, na.rm = TRUE),
              min(env_summary$p_val, na.rm = TRUE)),
    Status = ifelse(min(env_summary$p_val, na.rm = TRUE) < 0.05, "Biased", "Balanced"),
    stringsAsFactors = FALSE
  )

  res <- dplyr::bind_rows(geo_df, env_df)
  res$Value <- round(res$Value, 3)

  class(res) <- c("FoldsSummary", "data.frame")
  return(res)
}


#' Print FoldsSummary
#'
#' @param x An \code{EnvDiagnostic} object.
#' @param ... Additional arguments.
#' @return The \code{EnvDiagnostic} object invisibly.
#' @export
#' @family blocks diagnostics
#'
print.FoldsSummary <- function(x, ...) {
  cat("\n==========================================\n")
  cat("   isdmtools: INTEGRATED FOLD SUMMARY     \n")
  cat("==========================================\n\n")

  print.data.frame(x, row.names = FALSE)

  cat("\n------------------------------------------\n")

  # Schoener's D > 0.6 is a common threshold for good niche representativeness
  is_sep  <- all(x$Status[x$Domain == "Geographic"] == "Separated")
  is_bal  <- all(x$Status[x$Domain == "Environmental"] == "Balanced")
  avg_d   <- x$Value[x$Metric == "Median Overlap (D)"]

  cat("CONCLUSION: ")
  if (is_sep && is_bal && (is.na(avg_d) || avg_d > 0.6)) {
    cat("Folds are spatially independent \nand environmentally representative.\n")
  } else {
    cat("Potential issues detected. Review \nindividual diagnostic plots for bias.\n")
  }
  cat("==========================================\n")

  invisible(x)
}
