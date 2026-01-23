
#-----------------------------------------------------------------------
#--- Creating spatial cross-validation folds for mutisource datasets
#-----------------------------------------------------------------------

#' @title Create block cross-validation folds for multisource spatial datasets
#'
#' @description
#' A constructor function for the `DataFolds` S3 class. It binds multiple `sf` datasets
#' into a single object and generates spatially or environmentally-separated cross-validation folds.
#'
#' @param datasets A named list of `sf` objects. Each list element should be a
#'   spatial dataset with its name corresponding to the list element's name.
#' @param region.polygon An `sf` object representing the study area polygon.
#' @param k integer. It specifies the number of folds (k-fold cross-validation).
#' @param seed integer. It sets seed for reproducibility.
#' @param cv.method character. It specifies the spatial cross-validation method to use.
#' Options are `"cluster"` (default) or `"spatial"`, see \link[blockCV]{cv_cluster} or \link[blockCV]{cv_spatial} functions.
#' For `"block"`, `"buffer"`, `"location"`, or `"nndm"`, see the corresponding functions in the \code{spatialsample} package.
#' @param ... Additional arguments to be passed to the underlying blocking function in \code{blockCV} or \code{spatialsample}.
#'
#' @details
#' This function first binds all datasets into a single `sf` object. It then applies the chosen blocking method to create spatial folds.
#' The fold IDs are added to the combined data object, and the original datasets and other relevant information are
#' stored in the returned `DataFolds` object.
#'
#' The `"cluster"` method from the \code{blockCV} package, supports both spatial and environmental clustering.
#' Methods `"block"`, `"buffer"`, `"nndm"`, and `"location"`are useful for distance-based exclusion (buffering)
#' or leaving out specific groups/locations (e.g., using `group = "column_name"` with the \link[spatialsample]{spatial_leave_location_out_cv} method).
#' Use `"spatial"` for \code{blockCV} grid-blocking, or `"block"` for \code{spatialsample} grid-blocking.
#'
#' The behavior of \code{create_folds} depends on the \code{cv.method} chosen.
#' Several methods require specific arguments passed via the ellipsis (\code{...}):
#' \itemize{
#'   \item \bold{\code{location}}: Requires a \code{group} argument (character).
#'   This should be the name of a column in your datasets representing
#'   independent units like \code{"site_id"}, \code{"year"}, or \code{"observer"}.
#'   This tests model generalizability across these factors.
#'   \item \bold{\code{block} & \code{buffer}}: Accept \code{radius} and \code{buffer}
#'   arguments to define the size of the test/assessment areas and the width of
#'   the exclusion zones, respectively.
#'   \item \bold{\code{nndm}}: Requires \code{prediction_sites}, an \code{sf}
#'   object representing the area where the model will be projected.
#' }
#'
#' @return An S3 object of class `DataFolds` containing the combined data,
#' fold information, and the original datasets.
#' @export
#' @family spatial blocking methods
#'
#' @examples
#' \dontrun{
#' # Create some dummy sf data with different columns
#' library(sf)
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
#'   count = rpois(50, 5)
#' ) %>% st_as_sf(coords = c("x", "y"), crs = 4326)
#'
#' # Create a list of datasets
#' datasets_list <- list(Presence = presence_data, Count = count_data)
#'
#' # Create a dummy polygon for the region (e.g. Benin's minimum bounding rectangle)
#' ben_coords <- matrix(c(0, 6, 4, 6, 4, 13, 0, 13, 0, 6), ncol = 2, byrow = TRUE)
#' ben_sf <- st_sfc(st_polygon(list(ben_coords)), crs = 4326)
#' ben_sf <- st_sf(data.frame(name = "Benin"), ben_sf)
#'
#' # Create a DataFolds object using the default 'cluster' method
#' my_folds <- create_folds(datasets_list, ben_sf, k = 5)
#' print(my_folds)
#' }
#'
#' @references
#' Mahoney MJ, Johnson LK, Silge J, Frick H, Kuhn M, Beier CM. Assessing the performance of spatial cross-validation approaches for models of spatially structured data. _arXiv_ (2023) \doi{10.48550/arXiv.2303.07334}
#'
#' Roberts DR, Bahn V, Ciuti S, Boyce MS, Elith J, Guillera-Arroita G, Hauenstein S, Lahoz-Monfort JJ, Schröder B, Thuiller W, et al. Cross-validation strategies for data with temporal, spatial, hierarchical, or phylogenetic structure. _Ecography_ (2017) 40:913–929. \doi{10.1111/ecog.02881}
#'
#' Valavi R, Elith J, Lahoz-Monfort JJ, Guillera-Arroita G. blockCV: an R package for generating spatially or environmentally separated folds for k-fold cross-validation of species distribution models. _bioRxiv_ (2018). \doi{10.1101/357798}
#'
create_folds <- function(datasets, region.polygon = NULL, k = 5, seed = 23, cv.method = "cluster", ...) {

  xy_all <- bind_datasets(datasets)
  n_obs  <- nrow(xy_all)

  set.seed(seed)
  folds_cv <- switch(cv.method,
                     "cluster" = blockCV::cv_cluster(x = xy_all, k = k, ...),
                     "spatial" = blockCV::cv_spatial(x = xy_all, k = k, ...),
                     "nndm" = {
                       .check_suggests("spatialsample")
                       res <- spatialsample::spatial_nndm_cv(xy_all, ...)
                       list(folds_ids = .extract_spatialsample_ids(res, n_obs))
                     },
                     "buffer" = {
                       .check_suggests("spatialsample")
                       res <- spatialsample::spatial_buffer_vfold_cv(xy_all, v = k, ...)
                       list(folds_ids = .extract_spatialsample_ids(res, n_obs))
                     },
                     "location" = {
                       .check_suggests("spatialsample")
                       dots <- rlang::list2(...)
                       if (!"group" %in% names(dots)) {
                         stop("The 'group' argument (column name) is required for 'location out' CV.")
                       }
                       res <- spatialsample::spatial_leave_location_out_cv(xy_all, v = k, ...)
                       list(folds_ids = .extract_spatialsample_ids(res, n_obs))
                     },
                     "block" = {
                       .check_suggests("spatialsample")
                       res <- spatialsample::spatial_block_cv(xy_all, v = k, ...)
                       list(folds_ids = .extract_spatialsample_ids(res, n_obs))
                     },
                     stop("Invalid `cv.method`. Supported: 'cluster', 'spatial', 'block', 'buffer', 'location', 'nndm'.")
                     )

  xy_all$folds_ids <- folds_cv$folds_ids

  object <- list(
    data_all = xy_all,
    original_datasets = datasets,
    folds_info = folds_cv,
    dataset_names = names(datasets),
    k = k,
    region_polygon = region.polygon
  )
  class(object) <- "DataFolds"
  return(object)
}

#--- S3 method for DataFolds objects ----
#' @title Extract a specific fold from a DataFolds object
#'
#' @description
#' This method splits the data into training and testing sets for a given fold.
#' It uses the original datasets to ensure each returned `sf` object contains only
#' its original columns, avoiding `NA` values from the binding process.
#'
#' @param object A `DataFolds` S3 object.
#' @param fold An integer specifying the fold ID to be extracted as the test set.
#' @param ... Additional arguments (not used by this method).
#'
#' @return A list containing two named elements: `train` and `test`.
#' Each of these elements is a named list of `sf` objects, with names corresponding to the original datasets.
#' @method extract_fold DataFolds
#' @export
#' @family spatial blocking methods
#'
#' @examples
#' \dontrun{
#' # Create a list of some dummy sf data with different columns
#' library(sf)
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
#'   count = rpois(50, 5)
#' ) %>% st_as_sf(coords = c("x", "y"), crs = 4326)
#'
#' datasets_list <- list(Presence = presence_data, Count = count_data)
#'
#' # Create a dummy polygon for the region
#' ben_coords <- matrix(c(0, 6, 4, 6, 4, 13, 0, 13, 0, 6), ncol = 2, byrow = TRUE)
#' ben_sf <- st_sfc(st_polygon(list(ben_coords)), crs = 4326)
#' ben_sf <- st_sf(data.frame(name = "Benin"), ben_sf)
#'
#' # Create a DataFolds object with the default 'cluster' method
#' my_folds <- create_folds(datasets_list, ben_sf, k = 5)
#'
#' # Extract the desired fold (e.g. third fold of k-fold cross-validation)
#' splits_fold_3 <- extract_fold(my_folds, fold = 3)
#'
#' # Access the clean datasets for modeling
#' train_data_presence <- splits_fold_3$train$Presence
#' test_data_count <- splits_fold_3$test$Count
#'
#' # Check the columns to see that they are correct
#' names(train_data_presence)
#' names(test_data_count)
#' }
#'
extract_fold.DataFolds <- function(object, fold, ...) {
  if (!(fold %in% 1:object$k)) {
    stop(paste("Invalid fold number. Must be between 1 and", object$k))
  }

  train_folds_splits <- object$data_all %>% dplyr::filter(.data$folds_ids != fold & !is.na(.data$folds_ids))
  test_folds_splits <- object$data_all %>% dplyr::filter(.data$folds_ids == fold)

  train_splits_list <- split(train_folds_splits, train_folds_splits$datasetName)
  test_splits_list  <- split(test_folds_splits, test_folds_splits$datasetName)

  train_data <- purrr::map2(object$original_datasets, train_splits_list, sf::st_filter)
  test_data <- purrr::map2(object$original_datasets, test_splits_list, sf::st_filter)

  return(list(
    train = train_data,
    test = test_data
  ))
}

# --- Generic method ---
#' @title Extract a specific fold from a data partition
#'
#' @description
#' A generic method to extract a specific fold from an object that contains
#' data partitions, for use in cross-validation.
#'
#' @param object An object from which a spatial fold can be extracted (e.g. `DataFolds` object)
#' @param ... Additional arguments passed to specific methods.
#'
#' @details
#' This function is generic, meaning it provides a consistent interface for
#' different types of objects. The method dispatched depends on the class of the `object` argument.
#' The primary purpose is to abstract the process of accessing training and testing data for a given fold,
#' making it easier to write generic cross-validation loops.
#'
#' @return A list containing two named elements: `train` and `test`.
#' The structure of these elements depends on the specific method used.
#' For the `DataFolds` class, each element is a named list of `sf` objects,
#' with names corresponding to the original datasets.
#'
#' @export
#' @family spatial blocking methods
#'
extract_fold <- function(object, ...) {
  UseMethod("extract_fold")
}

# --- Default method ---
#' @method extract_fold default
extract_fold.default <- function(x) {
  cat("Default method for class", sQuote(class(x)), ".\n")
}

#--- Print method ---
#' @title Print a summary of a DataFolds object
#'
#' @description
#' A method to print a concise summary of the `DataFolds` object, including the
#' number of folds, the datasets included, and a breakdown of individuals per
#' dataset across all folds.
#'
#' @param x A `DataFolds` S3 object.
#' @param ... Additional arguments (not used by this method).
#'
#' @return The object invisibly.
#' @method print DataFolds
#' @export
#' @family spatial blocking methods
#'
print.DataFolds <- function(x, ...) {
  cat("A DataFolds S3 object with", x$k, "folds.\n")
  cat("Datasets included:", paste(x$dataset_names, collapse = ", "), "\n\n")
  cat("Summary of individuals per dataset:\n")

  summary_df <- x$data_all %>%
    dplyr::mutate(fold = ifelse(is.na(.data$folds_ids), "Excluded (Buffer)", as.character(.data$folds_ids))) %>%
    dplyr::group_by(.data$datasetName, .data$folds_ids) %>%
    dplyr::summarise(n = dplyr::n(), .groups = "drop")

 print(summary_df)

  invisible(x)
}

#--- Plot method ---
#' @title Plot block cross-validation folds for multisource datasets
#'
#' @description
#' A method to visualize the spatial blocks and the corresponding train/test
#' partitions for each fold of a `DataFolds` object.
#'
#' @param x A `DataFolds` S3 object.
#' @param nrow An integer specifying the number of rows needed for the panel plot.
#' The default is 1.
#' @param ... Additional arguments (not used by this method).
#'
#' @return A `ggplot2` object that can be printed or saved.
#' @method plot DataFolds
#' @export
#' @family spatial blocking methods
#'
#' @examples
#' \dontrun{
#' # Create some dummy sf data
#' library(sf)
#' library(dplyr)
#'
#' set.seed(42)
#' presence_data <- data.frame(
#'   x = runif(100, 0, 4),
#'   y = runif(100, 6, 13),
#'   site = rbinom(100, 1, 0.6)
#' ) %>%
#'   st_as_sf(coords = c("x", "y"), crs = 4326)
#'
#' count_data <- data.frame(
#'   x = runif(50, 0, 4),
#'   y = runif(50, 6, 13),
#'   count = rpois(50, 5)
#' ) %>%
#'   st_as_sf(coords = c("x", "y"), crs = 4326)
#'
#' # Create a list of datasets
#' datasets_list <- list(Presence = presence_data, Count = count_data)
#'
#' # Create a dummy polygon for the region (e.g., Benin)
#' ben_coords <- matrix(c(0, 6, 4, 6, 4, 13, 0, 13, 0, 6), ncol = 2, byrow = TRUE)
#' ben_sf <- st_sfc(st_polygon(list(ben_coords)), crs = 4326)
#' ben_sf <- st_sf(data.frame(name = "Benin"), ben_sf)
#'
#' # Create a DataFolds object for the example
#' my_folds <- create_folds(datasets_list, ben_sf, k = 5)
#'
#' # Now you can plot the object
#' plot_cv <- plot(my_folds)
#' print(plot_cv)
#'
#' # You can even customize the plot output (e.g. adjusting the axes breaks)
#' plot_cv <- plot_cv +
#'   ggplot2::scale_x_continuous(breaks = seq(0, 4, 1)) +
#'   ggplot2::scale_y_continuous(breaks = seq(6, 13, 2))
#' print(plot_cv)
#' }
#'
plot.DataFolds <- function(x, nrow = 1, ...) {

  num_datasets <- length(x$dataset_names)
  shapes <- c(16, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15)[1:num_datasets]
  names(shapes) <- x$dataset_names

  plot_data_list <- purrr::map(1:x$k, function(fold_id) {
    x$data_all %>%
      dplyr::mutate(
        Set = ifelse(.data$folds_ids == fold_id, "Test", "Train"),
        fold_panel = factor(fold_id)
      )
  })

  folds_xy_expanded <- dplyr::bind_rows(plot_data_list)

  plot_cv <- ggplot2::ggplot(folds_xy_expanded) +
    ggspatial::geom_sf(data = x$region_polygon, fill = NA, color = "grey20") +
    ggspatial::geom_sf(ggplot2::aes(color = .data$Set, shape = .data$datasetName), size = 1.2) +
    ggplot2::scale_color_manual(name = "Partition", values = c("Train" = "blue", "Test" = "orange")) +
    ggplot2::scale_shape_manual(name = "Dataset", values = shapes) +
    ggplot2::facet_wrap(~ fold_panel,
                        labeller = ggplot2::labeller(fold_panel = function(x) paste("Fold", x)),
                        nrow = nrow) +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(
      legend.position = "right",
      strip.text = ggplot2::element_text(size = 13),
      panel.grid.major = ggplot2::element_line(color = "grey80"),
      panel.grid.minor = ggplot2::element_line(color = "grey90")
    ) +
    ggplot2::labs(title = paste("Block Cross-Validation Folds"),
                  x = "Longitude",
                  y = "Latitude") +
    ggspatial::annotation_north_arrow(location = "tl", height = grid::unit(0.6, "cm"),
                                      width = grid::unit(0.3, "cm")) +
    ggspatial::annotation_scale(location = "br", bar_cols = c("grey60", "white"))

  return(plot_cv)
}

#-- Helper function to bind sf objects into a single object ----

#' @title Bind a list of spatial datasets for spatial blocking.
#' @description Helper function to create a single `sf` object with a `datasetName` column
#' useful for spatial blocking in cross-validation.
#'
#' @param datasets A named list of `sf` datasets having the same geographical projection and
#' coordinates reference system (CRS)
#'
#' @return An `sf` object containing all datasets along with their corresponding name.
#' @noRd
#'
bind_datasets <- function(datasets) {
  if (!is.list(datasets) || is.null(names(datasets))) {
    stop("Input must be a named list of 'sf' objects.", call. = FALSE)
  }

  for (ds_name in names(datasets)) {
    if (!inherits(datasets[[ds_name]], "sf")) {
      stop(sprintf("All elements in 'datasets' must be 'sf' objects. '%s' is not.", ds_name), call. = FALSE)
    }
    if (!inherits(sf::st_geometry(datasets[[ds_name]]), "sfc_POINT")) {
      warning(sprintf("Dataset '%s' in 'datasets' does not contain point geometries and could cause issue to values extraction.", ds_name), call. = FALSE)
    }
  }

  datasets_labeled <- purrr::map2(datasets, names(datasets), ~ .x %>% dplyr::mutate(datasetName = .y))
  bound_data <- dplyr::bind_rows(datasets_labeled)
  bound_data$datasetName <- factor(bound_data$datasetName, levels = names(datasets))

  return(bound_data)
}

#-- Other helper functions for 'spatialsample' blocking
.extract_spatialsample_ids <- function(rset, n) {
  ids <- rep(NA_integer_, n)
  for (i in seq_len(nrow(rset))) {
    assess_idx <- as.integer(rset$splits[[i]], data = "assessment")
    if (length(assess_idx) > 0) {
      ids[assess_idx] <- i
    }
  }
  return(ids)
}

.check_suggests <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(sprintf("Package '%s' is required for this CV method. Please install it.", pkg),
         call. = FALSE)
  }
}
