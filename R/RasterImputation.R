

#---------------------------------------------------------------------------------------
#--- Imputing NA cells in SpatRaster with the nearest values using a moving window
#---------------------------------------------------------------------------------------

#' @title Fill NA cells with the nearest cells values
#'
#' @description
#' Function to impute a raster object (`SpatRaster` or `rasterLayer`) for spatial modeling tools that cannot handle NAs.
#' The function fill in NA cells with the nearest cells values using a moving window on missing cells.
#' @param x A raster file of the covariate to impute. For multiple `SpatRaster` or `rasterStack`, use `lapply()` and `lapply()` methods along this function.
#' @param boundary A spatial polygon object (`spatVector` or `sf`) to be used to mask cells outside the study region, as the output layer is extended. Defautls to `NULL`.
#' @param ... Additional arguments passed to \link[terra]{focal}.
#'
#' @return A `SpatRaster` object in which all NA cells are filled in by th nearest cells.
#' @export
#'
fill_na_near <- function(x, boundary = NULL, ...) {

  if(!inherits(x, c("SpatRaster", "RasterLayer"))) {
    stop(sprintf("'%' must be a 'SpatRaster' or 'RasterLayer' object.", deparse(substitute(x))), call. = FALSE)
  } else {
    if(inherits(x, "RasterLayer")) {
      message(sprintf("Converting '%s' into a 'spatRaster' object.", deparse(substitute(x))))
      filled <- terra::rast(x)
    } else {
      filled  <- x
    }
    w       <- 1
    to_fill <- terra::global(filled, function(x) any(is.na(x)))[,1]
    while(to_fill) {
      w       <- w + 2
      filled  <- terra::focal(filled, w = w, fun = mean, na.policy = "only", na.rm = TRUE, ...)
      to_fill <- terra::global(filled, function(x) any(is.na(x)))[,1]
    }

    if(!is.null(boundary)) {
      bndr <- boundary
      if (inherits(boundary, "SpatVector")) {
        if (terra::geomtype(boundary) != "polygons") {
          stop(sprintf("'%s' SpatVector must contain only polygons geometries.", boundary), call. = FALSE)
        }
      } else if (inherits(boundary, "sf")) {
        if (!all(sf::st_geometry_type(boundary) %in% c("POLYGON", "MULTIPOLYGON"))) {
          stop(sprintf("'%s' sf object must contains only POLYGON geometries.", boundary), call. = FALSE)
        }
      } else {
        stop(sprintf("Unsupported geometries. '%s' must be a polygon from sf or terra.", boundary), call. = FALSE)
      }
      filled  <- terra::mask(filled, bndr)
    }
  }
  names(filled) <- names(x)

  return(filled)
}
