

#---------------------------------------------------------------------
#--- Function to generate background points without NA cells (10)
#---------------------------------------------------------------------

#' @title Generate background points
#' @description A constructor function to generate background points for computing evaluation scores for presence-only data.
#' It exclude NA cells from the sample, and eventually the observed locations if needed.
#'
#' Inputs:
#' @param mask A raster object (e.g. `spatRaster`) to be used as mask (preferably, the predicted intensity or habitat suitability).
#' @param points Spatial points (`data.frame`, `sf` or `SpatVector` objects) that would be excluded from the background sample.
#' @param method The sampling technique to select pixels from the raster mask (see \link[terra]{spatSample}), defaulted to `random`.
#' @param n An integer specifying the number of pseudo-absence points to sample for presence-only data. Default is 1000.
#' @param values A logical value indicating whether sampled cells values will be returned or not.
#' @param cells A logical value indicating whether sampled cells numbers will be returned.
#' @param xy A logical value indicating whether the locations of sampled cells will be returned.
#' @param as.points A logical value indicating whether spatial points object will be returned.
#' @param ... Additional arguments passed to \link[terra]{spatSample}.
#'
#' @return A list object containing a `spatRaster` object (with `points` excluded if not `NULL`) and the background points.
#' @export
#'
#' @examples
#' \dontrun{
#' library(terra)
#' set.seed(123)
#' r  <- rast(nrows = 100, ncols = 100, xmin = 0, xmax = 10, ymin = 0, ymax = 10)
#' terra::values(r) <- runif(ncell(r))
#' pts <- spatSample(r, size = 100, xy = TRUE, values = FALSE)
#'
#' # Requesting few points with their x and y coordinates
#' set.seed(235)
#' bg_sample1 <- sample_bg_points(r, points = pts, n = 500, xy = TRUE, cells = FALSE)
#' plot(bg_sample1)
#' print(bg_sample1)
#'
#' # Requesting points more than available non-NA cells
#' bg_sample2 <- sample_bg_points(r, points = pts, n = 10000, xy =TRUE, cells = FALSE)
#' dim(bg_sample2$bg)
#' plot(bg_sample2)
#' }
#'
#' @seealso \code{\link{plot.BackgroundPoints}}, \code{\link{print.BackgroundPoints}}
#'
sample_bg_points <- function(mask, points = NULL, n = 500, method = "random", values = FALSE,
                             cells = FALSE, xy = TRUE, as.points = FALSE, ...) {

  if (!any(c(cells, xy, as.points)) && values) {
    stop("At least one of 'cells', 'xy', or 'as.points' must be TRUE.")
  }

  if (!inherits(mask, c("SpatRaster", "RasterLayer"))) {
    stop(sprintf("'%s' must be a raster object from terra or raster package.", deparse(substitute(mask))), call. = FALSE)
  }

  if (inherits(mask, "RasterLayer")) {
    message(sprintf("Converting '%s' into a 'spatRaster' object.", deparse(substitute(mask))))
    mask <- terra::rast(mask)
  }

  if (!is.null(points)) {
    if (inherits(points, "SpatVector")) {
      if (terra::geomtype(points) != "points") {
        stop(sprintf("'%s' SpatVector must contain only point geometries.", deparse(substitute(points))), call. = FALSE)
      }
      pts <- terra::geom(points)[, c("x", "y")]
    } else if (inherits(points, "sf")) {
      if (!all(sf::st_geometry_type(points) %in% c("POINT", "MULTIPOINT"))) {
        stop(sprintf("'%s' sf object must contain only POINT geometries.", deparse(substitute(points))), call. = FALSE)
      }
      pts <- sf::st_coordinates(points)[, c("X", "Y")]
    } else if (inherits(points, c("SpatialPoints", "SpatialPointsDataFrame"))) {
      pts <- sp::coordinates(points)[, 1:2]
    } else if (is.matrix(points) || is.data.frame(points)) {
      if (ncol(points) < 2) {
        stop(sprintf("'%s' must have at least two columns (x and y).", deparse(substitute(points))), call. = FALSE)
      }
      pts <- as.matrix(points)[, 1:2]
    } else {
      stop("Unsupported 'points' format. Must be object from sf, sp, or terra, or a matrix/data.frame.", call. = FALSE)
    }

    excluded_cells <- terra::cellFromXY(mask, pts)
    mask_modified  <- mask
    mask_modified[excluded_cells] <- NA
  } else {
    warning("'points' argument is NULL. You can supply known 'points' to be excluded from the sample.", call. = FALSE)
    mask_modified <- mask
  }

  n_valid <- sum(!is.na(terra::values(mask_modified)))
  if (n > n_valid) {
    warning("Requested more background points than available non-NA pixels.", call. = FALSE)
  }

  result <- tryCatch({
    bg_points <- terra::spatSample(
      mask_modified,
      size = n,
      method = method,
      na.rm = TRUE,
      values = values,
      cells = cells,
      xy = xy,
      as.points = as.points,
      ...
    )
    structure(list(bg = bg_points, mask = mask_modified), class = "BackgroundPoints")
  }, error = function(e) {
    message(sprintf("⚠️ Sampling background points failed: %s", conditionMessage(e)))
    return(NULL)
  })

  return(result)
}


#-------------------------------------------
#--- Useful methods for BackgroundPoints
#-------------------------------------------

#--- Print method for BackgroundPoints class ---

#' @title Print method for the class `BackgroundPoints`
#'
#' @description
#' A method to print the Background points generated from `BackgroundPoints` object for model evaluation.
#'
#' @param x A `BackgroundPoints` S3 object.
#' @param ... Additional arguments (not used by this method).
#'
#' @return The object invisibly.
#' @method print BackgroundPoints
#' @export
#' @seealso \code{\link{plot.BackgroundPoints}}
#'
print.BackgroundPoints <- function(x, ...) {
  cat("Generated Background Points:\n")
  print_bg(x$bg)
  cat("Modified Mask available in `$mask`.\n")

  invisible(x)
}


#--- Plot method for BackgroundPoints objects ----

#' @title Plot method for the class `BackgroundPoints`.
#'
#' @description
#' A method to visualize the background points generated in `BackgroundPoints` object for model evaluation,
#' including the `points` locations that have been excluded from the sample (white color) if there are
#' provided to \code{sample_bg_points()} function.
#'
#' @param x A `BackgroundPoints` S3 object.
#' @param ... Additional arguments (not used by this method).
#'
#' @return The object invisibly.
#' @method plot BackgroundPoints
#'
#' @export
#' @importFrom graphics points
#'
#' @examples
#' \dontrun{
#' library(terra)
#' set.seed(123)
#' r  <- rast(nrows = 100, ncols = 100, xmin = 0, xmax = 10, ymin = 0, ymax = 10)
#' terra::values(r) <- runif(ncell(r))
#' pts <- spatSample(r, size = 100, xy = TRUE, values = FALSE)
#'
#' # Requesting few points with their x and y coordinates
#' set.seed(235)
#' bg_sample1 <- sample_bg_points(r, points = pts, n = 500, xy = TRUE, cells = FALSE)
#' plot(bg_sample1)
#' print(bg_sample1)
#'
#' # Requesting points more than available non-NA cells
#' bg_sample2 <- sample_bg_points(r, points = pts, n = 10000, xy =TRUE, cells = FALSE)
#' dim(bg_sample2$bg)
#' plot(bg_sample2)
#' }
#'
#' @seealso \code{\link{print.BackgroundPoints}}
#'
plot.BackgroundPoints <- function(x, ...) {

  extract_points <- function(bg, mask) {
    if (inherits(bg, "SpatVector")) return(terra::crds(bg))
    if (is.matrix(bg) || is.data.frame(bg)) {
      if (all(c("x", "y") %in% colnames(bg))) return(bg[, c("x", "y")])
      if (ncol(bg) == 1 || ncol(bg) == 2) return(terra::xyFromCell(mask, bg[, "cell"]))
    }
    stop("Unsupported format for background points.")
  }

  pts <- extract_points(bg = x$bg, mask = x$mask)
  plot(x$mask, ...)
  graphics::points(pts, col = "red", pch = 20, cex = 0.8)

  invisble(x)
}

#--- Helper function to print background points based on its type ----

#' @title Helper function to print background points
#'
#' @description Print background points based on its format.
#' @param A `BackgroundPoints` S3 object
#'
#' @return An object that can be printed.
#' @importFrom utils head tail
#' @noRd
#'
print_bg <- function(x) {
  if (inherits(x, "SpatVector")) {
    print(x)
  } else if (is.data.frame(x) || is.matrix(x)) {
    n <- nrow(x)
    print(utils::head(x, 10), row.names = FALSE)
    if (n > 10) {
      cat(rep("... ", ncol(x) + 1), "\n")
      tail_output <- utils::tail(x, 5)
      colnames(tail_output) <- rep("", ncol(tail_output))
      print(tail_output, row.names = FALSE)
    }
  } else {
    print(utils::head(x, 10))
  }
  cat("\n")
}


