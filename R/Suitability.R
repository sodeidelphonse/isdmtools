
#-----------------------------------------------------------
#--- Function to compute an habitat suitability from ISDM
#-----------------------------------------------------------

#' @title Compute a unified suitability index from integrated spatial model predictions.
#'
#' @description
#' This function converts the linear predictor (`eta`) from a fitted integrated spatial model
#' into a unified suitability index, which can be interpreted as a probability of a species presence
#' or species *hotspots*. It also returns the intensity or expected count, depending on whether the model has offset or not.
#'
#' @param x A `data.frame` containing `X` and `Y` coordinates and the column(s) of the predicted linear predictor variables (e.g., mean, standard deviation and quantiles for Bayesian models).
#' It can typically be a standardised grid-based output from a \link{prepare_predictions} call to various classes of spatial prediction on a linear scale, e.g. from the `PointedSDMs` or `inlabru` packages.
#' It can also be a `SpatRaster` object containing the prediction variables as the layers' names.
#' @param post.stat A character vector specifying the column or layer name(s) to use for the predictions. Defaults to "mean".
#' @param output.format A character string indicating the desired output format. Must be one of "prob" (for a probability-based suitability index) or "response" (for the expected count or rate).
#' @param response.type A character string specifying the type of response data the model was fitted with.
#'   Must be one of "joint.po" (joint model including presence-only data), "count.pa" (joint model with count and presence-absence), "po" (single presence-only model), "count" (count model), or "pa" (presence-absence model).
#' @param has.offset A logical value. For `count.pa`, `count` and `pa` models, this should be `TRUE` if the linear predictor includes an explicit area offset.
#'   This argument is not used for "po" or "joint.po" models. Defaults to `FALSE` (sensible value).
#' @param projection A character string or object specifying the coordinate reference system (CRS) for the output raster. Defaults to `NULL`.
#' @param ... Additional arguments passed to the S4 method \link[terra]{rast} for signature 'data.frame'.
#'
#' @return A SpatRaster object (single or multi-layer) representing either the suitability
#'   index (probabilities, if `output.format = "prob"`) or the expected count/rate (if `output.format = "response"`).
#'
#' @details
#' This function implements a unified framework for converting model predictions to a
#' comparable suitability index. The method relies on the *Inhomogeneous Poisson Process (IPP)*
#' theory, where the linear predictor `eta` is related to the probability of presence
#' via the inverse complementary log-log link: `p(presence) = 1 - exp(-scaling * exp(eta))`.
#'
#' The `scaling` factor is determined by the `response.type` and the `has.offset` arguments:
#'
#' \itemize{
#'   \item For PO models (single or part of a joint model), `eta` is always a log-rate, and the `scaling` is automatically set to the cell area.
#'   \item For Count or PA models, if `has.offset = TRUE`, `eta` is a log-rate, and `scaling` is the cell area.
#'   \item In all other cases (`has.offset = FALSE`), `eta` is treated as a log of the expected count/probability, and `scaling` is set to `1`.
#' }
#'
#' @export
#' @family prediction analyses
#'
#' @examples
#' \dontrun{
#' # Simulate a sample data.frame with X, Y, and a linear predictor
#' set.seed(42)
#' x <- expand.grid(X = seq(0, 50, 1), Y = seq(0, 50, 1))
#' x$eta <- rnorm(nrow(x), mean = 0, sd = 1)
#'
#' # Generate a suitability index for a Presence-Absence model
#' pa_probability <- suitability_index(
#'   x,
#'   post.stat = "eta",
#'   response.type = "pa"
#' )
#' plot(pa_probability)
#'
#' # Generate an expected mean for a Count model
#' count_expected_mean <- suitability_index(
#'   x,
#'   post.stat = "eta",
#'   response.type = "count",
#'   output.format = "response"
#' )
#' plot(count_expected_mean)
#' }
#'
#' @references
#' Dorazio RM. Accounting for imperfect detection and survey bias in statistical analysis of presence-only data. _Global Ecology and Biogeography_ (2014) 23:1472–1484.
#'
#' Fithian W, Elith J, Hastie T, Keith DA. Bias correction in species distribution models: pooling survey and collection data for multiple species. _Methods in Ecology and Evolution_ (2015) 6:424–438. \doi{10.1111/2041-210X.12242}
#'
suitability_index <- function(x,
                              post.stat = "mean",
                              output.format = c("prob", "response"),
                              response.type = c("joint.po", "count.pa", "po", "count", "pa"),
                              has.offset = FALSE,
                              projection = NULL, ...) {

  output.format <- match.arg(output.format)
  response.type <- match.arg(response.type)

  if (is.data.frame(x)) {
    if (!all(c("X", "Y") %in% names(x))) {
      stop("'x' must contain 'X' and 'Y' coordinates when it is a data.frame.", call. = FALSE)
    }
    if (!all(post.stat %in% names(x))) {
      missing_stats <- post.stat[!post.stat %in% names(x)]
      stop(paste0("The following 'post.stat' columns were not found in the input data.frame: ",
                  paste(missing_stats, collapse = ", ")), call. = FALSE)
    }
    eta_rast <- terra::rast(x[, c("X", "Y", post.stat)], type = "xyz", crs = projection, ...)

  } else if (inherits(x, "SpatRaster")) {
    if (!all(post.stat %in% names(x))) {
      missing_stats <- post.stat[!post.stat %in% names(x)]
      stop(paste0("The following 'post.stat' layers were not found in the input raster: ",
                  paste(missing_stats, collapse = ", ")), call. = FALSE)
    }
    eta_rast <- x[[post.stat]]
  } else {
    stop("Input 'x' must be a data.frame or a SpatRaster object.", call. = FALSE)
  }

  #--- Determine the correct scaling factor based on model type and has.offset
  if (response.type %in% c("po", "joint.po")) {
    scaling_factor <- terra::res(eta_rast)[1] * terra::res(eta_rast)[2]
  } else if (has.offset) {
    scaling_factor <- terra::res(eta_rast)[1] * terra::res(eta_rast)[2]
  } else {
    scaling_factor <- 1
  }

  #--- The suitability index and the expected count/rate (response)
  prob_rast <- 1 - exp(-scaling_factor * exp(eta_rast))
  expected_response <- exp(eta_rast)

  names(prob_rast) <- post.stat
  names(expected_response) <- post.stat

  out <- if (output.format == "prob") {
    prob_rast
  } else {
    expected_response
  }

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
    #message("Detected point-based data by the presence of a .block column. Filtering to base map.")
    if (!inherits(data_to_prepare, "sf")) {
      stop("Detected point-based predictions, but it is not an sf object.", call. = FALSE)
    }

    if(!is.null(base_map)) {
      return(sf::st_filter(data_to_prepare, base_map))
    } else return(data_to_prepare)

  } else {
    #message("Detected grid-based data. Ensuring X, Y coordinates are present.")
    if (inherits(data_to_prepare, "sf")) {
      coords <- sf::st_coordinates(data_to_prepare)
      df <- as.data.frame(data_to_prepare)
      df$geometry <- NULL
      return(cbind(df, X = coords[, "X"], Y = coords[, "Y"]))

    } else {
      df <- as.data.frame(data_to_prepare)
      if (!all(c("X", "Y") %in% names(df))) {
        stop("Input data lacks 'X' and 'Y' coordinate columns for plotting.", call. = FALSE)
      }
      return(df)
    }
  }
}


