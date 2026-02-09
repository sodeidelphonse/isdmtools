
#-----------------------------------------------------------
#--- Function to compute an habitat suitability from ISDM
#-----------------------------------------------------------

#' @title Compute a unified suitability index from integrated spatial model predictions.
#'
#' @description
#' This function converts the linear predictor (`eta`) from a fitted integrated spatial model
#' into a unified suitability index, which can be interpreted as a probability of a species presence.
#' It also calculates the intensity(rate) or expected count for the count data, depending on whether the model has offset or not.
#'
#' @param x A `data.frame` containing `x` and `y` coordinates and the column(s) of the predicted linear predictor variables (e.g., mean, standard deviation and quantiles) or `SpatRaster`.
#' It can typically be a standardized grid-based output from a \link{prepare_predictions} call to various classes of spatial prediction on a linear scale, e.g. from the `PointedSDMs` or `inlabru` packages.
#' @param post_stat character. A vector specifying the column or layer name(s) to use for extracting the model predictions. Defaults to "mean".
#' @param output_format character. The desired output format and must be one of "prob" (probability-based suitability index), "response" (expected count or rate) or "linear" (linear predictor scale).
#' @param response_type character. The type of response data the model was fitted with.
#'   Must be one of "joint.po" (joint model including presence-only data), "count.pa" (joint model with count and presence-absence), "po" (single presence-only model), "count" (count model), or "pa" (presence-absence model).
#' @param has_offset logical. For `count.pa`, `count` and `pa` models, this should be `TRUE` if the linear predictor includes an explicit area offset.
#'   This argument is not used for "po" or "joint.po" models. Defaults to `FALSE`.
#' @param scale_independent logical. If `TRUE`, the scaling factor is set to 1, making the suitability index independent of the grid cell size. Defaults to `FALSE`.
#' @param projection character. The coordinate reference system (CRS) for the output raster. Defaults to `NULL`.
#' If `NULL` and `x` is a data.frame, an empty CRS is assigned to prevent errors.
#' @param ... Additional arguments passed to \link[terra]{rast}.
#'
#' @return A SpatRaster object representing the requested output format.
#'
#' @details
#' This function implements a unified framework for converting model predictions to a
#' comparable suitability index. The method relies on the *Inhomogeneous Poisson Process (IPP)*
#' theory, where the linear predictor `eta` is related to the probability of presence
#' via the inverse complementary log-log link as follows: \eqn{p(presence) = 1 - exp(-scaling \times exp(eta))}.
#'
#' The `scaling` factor is determined by the `response_type` and the `has_offset` arguments:
#' \itemize{
#'   \item For PO models (single or part of a joint model), `eta` is always a log-rate, and the `scaling` is set to the cell area.
#'   It is ignored if `scale_independent` is set to TRUE.
#'   \item For Count or PA models, if `has_offset = TRUE`, `eta` is a log-rate, and `scaling` is the cell area.
#'   \item In all other cases (`has_offset = FALSE`), `eta` is treated as a log of the expected count (count) or cloglog of probability (PA), and `scaling` is set to `1`.
#' }
#' If the raster is in a geographic coordinate system (longlat), the area is calculated in \eqn{km^2} using \link[terra]{cellSize}.
#' For projected systems, the area is the product of the resolutions (e.g., \eqn{km^2} if units are in km).
#'
#' @export
#' @family prediction analyses
#'
#' @examples
#' \dontrun{
#' library(terra)
#' # Simulate a sample data.frame with x, y, and a linear predictor
#' set.seed(42)
#' x <- expand.grid(x = seq(0, 50, 1), y = seq(0, 50, 1))
#' x$eta <- rnorm(nrow(x), mean = 0, sd = 1)
#'
#' # Create a simple raster object
#' x_rast <- rast(x)
#'
#' # Generate a suitability index for a Presence-Absence model
#' pa_probability <- suitability_index(
#'   x_rast,
#'   post_stat = "eta",
#'   response_type = "pa"
#' )
#' plot(pa_probability)
#'
#' # Create a binary map using a fixed threshold
#' binary_map <- app(pa_probability, function(x) ifelse(x < 0.5, 0, 1))
#' plot(binary_map)
#'
#' # Generate an expected mean (assume "eta" is from a count model)
#' count_expected_mean <- suitability_index(
#'   x_rast,
#'   post_stat = "eta",
#'   response_type = "count",
#'   output_format = "response"
#' )
#' plot(count_expected_mean)
#' }
#'
#' @references
#' Dorazio RM. Accounting for imperfect detection and survey bias in statistical analysis of presence-only data. _Global Ecology and Biogeography_ (2014) 23:1472–1484. \doi{10.1111/geb.12216}
#'
#' Fithian W, Elith J, Hastie T, Keith DA. Bias correction in species distribution models: pooling survey and collection data for multiple species. _Methods in Ecology and Evolution_ (2015) 6:424–438. \doi{10.1111/2041-210X.12242}
#'
#' Phillips SJ, Anderson RP, Dudík M, et al Opening the black box: an open‐source release of Maxent. _Ecography_ (2017) 40:887–893. \doi{10.1111/ecog.03049}
#'
suitability_index <- function(x,
                              post_stat = "mean",
                              output_format = c("prob", "response", "linear"),
                              response_type = c("joint.po", "count.pa", "po", "count", "pa"),
                              has_offset = FALSE,
                              scale_independent = FALSE,
                              projection = NULL, ...) {

  output_format <- match.arg(output_format)
  response_type <- match.arg(response_type)

  if (is.data.frame(x)) {
    names(x)[names(x) == "X"] <- "x"
    names(x)[names(x) == "Y"] <- "y"

    if (!all(c("x", "y") %in% names(x))) {
      stop("'x' must contain 'x' and 'y' coordinates when it is a data.frame.", call. = FALSE)
    }
    if (!all(post_stat %in% names(x))) {
      missing_stats <- post_stat[!post_stat %in% names(x)]
      stop(paste0("The following 'post_stat' columns were not found in the input data.frame: ",
                  paste(missing_stats, collapse = ", ")), call. = FALSE)
    }

    target_crs <- if(is.null(projection)) "" else projection
    eta_rast <- terra::rast(x[, c("x", "y", post_stat)], type = "xyz", crs = target_crs, ...)

  } else if (inherits(x, "SpatRaster")) {
    if (!all(post_stat %in% names(x))) {
      missing_stats <- post_stat[!post_stat %in% names(x)]
      stop(paste0("The following 'post_stat' layers were not found in the input raster: ",
                  paste(missing_stats, collapse = ", ")), call. = FALSE)
    }

    eta_rast <- x[[post_stat]]
    if (!is.null(projection)) terra::crs(eta_rast) <- projection

  } else {
    stop("Input 'x' must be a data.frame or a SpatRaster object.", call. = FALSE)
  }

  #--- Determine the correct scaling factor based on model type and has_offset
  if (scale_independent) {
    scaling_factor <- 1
  } else if (response_type %in% c("po", "joint.po") || has_offset) {
    if (isTRUE(terra::is.lonlat(eta_rast))) {
      scaling_factor <- terra::cellSize(eta_rast, unit = "km")
    } else {
      scaling_factor <- terra::res(eta_rast)[1] * terra::res(eta_rast)[2]
    }
  } else {
    scaling_factor <- 1
  }

  #--- The suitability index and the expected count/rate
  expected_response <- exp(eta_rast)
  prob_rast <- 1 - exp(-(scaling_factor * expected_response))

  out <- switch(output_format,
                prob = prob_rast,
                response = expected_response,
                linear = eta_rast
                )
  names(out) <- post_stat

  return(out)
}


#--- Generalized inverse cloglog function -------

#' Transform an intensity function into probability
#' @description Transform an intensity function into probability using the generalized inverse of `cloglog` function.
#'
#' @param eta The linear predictor from an integrated or individual model prediction
#' @param scaling A scaling factor which is 1 or equal to the area of grid cells of predicted response or intensity, depending on the model assumption.
#'
#' @return A `SpatRaster` object of predicted probabilities (suitability) or (relative/absolute) abundance depending on the response type.
#' @noRd
inv_cloglog <- function(eta, scaling = 1) {
  1 - exp(-scaling * exp(eta))
}


#--- Convert spatial predictions into a formatted output ---------------------

#' @title Obtain a formatted output from spatial predictions.
#' @description Function to transform a prediction data from various spatial models (e.g., `inlabru`, `PointedSDMs` or `GLMs` tools)
#' into an sf object if it is from points predictions (e.g., \link[fmesher]{fm_vertices}) or a data.frame with
#' the corresponding locations if it is from pixel grids (see \link[fmesher]{fm_pixels}). Spatial prediction data can also be obtained across a given region
#' using \code{expand.grid(x, y)} function, where 'x' and 'y' are geographical coordinates of the grids locations.
#'
#' @param prediction_data A model prediction which may be either an `sf` or a `data.frame` object or a raw prediction from the `inlabru-like` models.
#' The prediction can be on the response or linear predictor scale, depending on whether the output is for a model evaluation or visualization.
#' @param base_map An `sf` polygon having the same `crs` with the spatial locations used for predictions.
#'
#' @return A `data.frame` for grid-based predictions or `sf` object for point-based predictions.
#' @export
#' @family prediction analyses
#'
prepare_predictions <- function(prediction_data, base_map = NULL) {

  if(!is.null(base_map)) {
    if (!inherits(base_map, "sf")) {
      stop("'base_map' must be an sf object.", call. = FALSE)
    }
    if (!all(sf::st_geometry_type(base_map) %in% c("POLYGON", "MULTIPOLYGON"))) {
      stop("'base_map' must have polygon or multipolygon geometry.", call. = FALSE)
    }
  }

  data_to_prepare <- NULL
  if (inherits(prediction_data, "modISDM_predict")) {
    data_to_prepare <- prediction_data$predictions
  } else if (inherits(prediction_data, "bru_prediction")) {
    data_to_prepare <- prediction_data
  } else if (inherits(prediction_data, c("sf", "data.frame"))) {
    data_to_prepare <- prediction_data
  } else {
    stop("Unsupported prediction object class.", call. = FALSE)
  }

  if (any(c(".block", ".vertex") %in% names(data_to_prepare))) {
    if (!inherits(data_to_prepare, "sf")) {
      stop("Detected point-based predictions, but it is not an sf object.", call. = FALSE)
    }

    if (!is.null(base_map)) {
      pred_crs <- sf::st_crs(data_to_prepare)
      map_crs  <- sf::st_crs(base_map)

      if (is.na(pred_crs) || is.na(map_crs)) {
        warning("One of the inputs lacks a CRS. Spatial filtering may be inaccurate.", call. = FALSE)
      } else if (pred_crs != map_crs) {
        base_map <- sf::st_transform(base_map, pred_crs)
      }

      return(sf::st_filter(data_to_prepare, base_map))

    } else {
      return(data_to_prepare)
    }

  } else {
    if (inherits(data_to_prepare, "sf")) {
      coords <- sf::st_coordinates(data_to_prepare)

      df <- as.data.frame(data_to_prepare)
      geom_column <- attr(data_to_prepare, "sf_column")
      df[[geom_column]] <- NULL
      df <- df[, !(names(df) %in% c("x", "y", "X", "Y")), drop = FALSE]

      return(cbind(df, x = coords[, "X"], y = coords[, "Y"]))

    } else {
      df <- as.data.frame(data_to_prepare)
      names(df)[names(df) == "X"] <- "x"
      names(df)[names(df) == "Y"] <- "y"

      if (!all(c("x", "y") %in% names(df))) {
        stop("Input data lacks 'x' and 'y' coordinate columns for plotting.", call. = FALSE)
      }

      return(df)
    }
  }
}


