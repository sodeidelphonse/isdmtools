
#----------------------------------------------------------
#--- Multi-panel plot for spatial data visualization
#----------------------------------------------------------

#' Generate Multi-panel Maps from Spatial Model Predictions
#'
#' This function creates a multi-panel map for visualizing multiple prediction
#' variables from a species distribution model or other spatial model. It is
#' designed to be flexible, handling both grid-based (data.frame) and
#' point-based (sf) spatial predictions.
#'
#' The function internally reshapes the data from a wide format (with a column
#' for each prediction variable) to a long format suitable for plotting with `ggplot2::facet_wrap()`.
#' It automatically selects the appropriate geometry (`geom_tile()` for grids and `geom_sf()` for points) and conditional scales.
#' Users can also add to the map other spatial vector layers or customize the plot using ggplot2 syntax if needed.
#'
#' @param data A data frame, `sf` or `SpatRaster` object containing the prediction data.
#'   For grid-based data frame, it must contain columns named "X" and "Y" representing pixels' coordinates.
#' @param var.names A character vector of column names in `data` to be
#'   plotted on separate panels.
#' @param base.map An `sf` object to be plotted as a base layer underneath the
#'   prediction data (e.g., a background simple polygon). Defaults to `NULL`.
#'   Users can add additional vector geometries if needed using a ggplot2 syntax.
#' @param color.gradient A vector of valid colors to be used in the fill/color gradient.
#'   Defaults to `map.pal("viridis", 100)`.
#' @param legend.title A character string for the title of the color legend.
#' @param panel.labels An optional character vector of labels for the facet panels.
#'   The order should correspond to `var.names`.
#' @param nrow The number of rows for `ggplot2::facet_wrap()`. Defaults to an
#'   optimal layout chosen by `ggplot2`.
#' @param xaxis.breaks A numeric vector specifying the breaks for the x-axis.
#' @param yaxis.breaks A numeric vector specifying the breaks for the y-axis.
#'
#' @return A `ggplot` object representing the multi-panel plot that can be customized by the user.
#' @export
#'
#' @examples
#' \dontrun{
#' # --- Example with grid-based data ---
#' # Simulate a data frame with coordinates and two prediction variables
#' grid_data <- expand.grid(X = 1:100, Y = 1:100)
#' grid_data$mean <- rnorm(10000, mean = grid_data$X / 100, sd = 0.1)
#' grid_data$sd <- rgamma(10000, shape = 2, scale = 0.2)
#'
#' # Simulate a boundary map (e.g., a simple polygon)
#' library(sf)
#' boundary <- st_sfc(st_polygon(list(cbind(c(0, 100, 100, 0, 0), c(0, 0, 100, 100, 0)))))
#' boundary_sf <- st_sf(data.frame(id = 1), geometry = boundary)
#'
#' # Generate the map
#' generate_maps(
#'   data = grid_data,
#'   var.names = c("mean", "sd"),
#'   base.map = boundary_sf,
#'   color.gradient = c("white", "skyblue", "navy"),
#'   legend.title = "Prediction Value",
#'   panel.labels = c("Mean", "StDev"),
#'   nrow = 1
#' )
#'
#' # --- Example with point-based data (sf) ---
#' # Simulate an sf object with point data
#' library(sf)
#' set.seed(123)
#' points_sf <- st_as_sf(grid_data[sample(1:10000, 1000), ],
#'                      coords = c("X", "Y"), crs = 32631) # UTM CRS
#' points_sf <- prepare_predictions(points_sf)
#'
#' # Generate the map
#' generate_maps(
#'   data = points_sf,
#'   var.names = c("mean", "sd"),
#'   base.map = boundary_sf,
#'   color.gradient = c("white", "orange", "red"),
#'   legend.title = "Prediction Value",
#'   panel.labels = c("Mean", "StDev"),
#'   nrow = 1
#' )
#' }
#'
#' @importFrom terra map.pal
#' @importFrom grid unit
#' @family prediction analyses
#'
generate_maps <- function(data,
                          var.names = c("mean", "sd"),
                          base.map = NULL,
                          color.gradient = map.pal("viridis", 100),
                          legend.title = NULL,
                          panel.labels = NULL,
                          nrow = NULL,
                          xaxis.breaks = NULL,
                          yaxis.breaks = NULL) {

  if (!inherits(data, c("data.frame", "sf", "SpatRaster"))) {
    stop("'data' must be a data.frame, an sf object, or a SpatRaster.", call. = FALSE)
  }

  if(inherits(data, "SpatRaster")) {
    data <- as.data.frame(data, xy = TRUE)
    data <- data %>% dplyr::rename(X = x, Y = y)
  }
  if (is.data.frame(data) && !inherits(data, "sf") && !all(c("X", "Y") %in% names(data))) {
    stop("'data' is a data.frame but is missing 'X' and 'Y' coordinate columns.", call. = FALSE)
  }

  if (!is.character(var.names) || length(var.names) == 0) {
    stop("'var.names' must be a non-empty character vector of columns names.", call. = FALSE)
  }

  if (!all(var.names %in% names(data))) {
    missing_vars <- setdiff(var.names, names(data))
    stop(paste("The following variables in 'var.names' were not found in 'data':",
               paste(missing_vars, collapse = ", ")), call. = FALSE)
  }

  if(!is.null(base.map)) {
    if (!inherits(base.map, "sf")) {
      stop("'base.map' must be an sf object.", call. = FALSE)
    }
    if (!all(sf::st_geometry_type(base.map) %in% c("POLYGON", "MULTIPOLYGON"))) {
      stop("'base.map' must have polygon or multipolygon geometry.", call. = FALSE)
    }
  }

  if (!is.character(color.gradient) || length(color.gradient) < 2) {
    stop("'color.gradient' must be a character vector of at least two colors.", call. = FALSE)
  }

  if (!is.null(legend.title) && is.character(legend.title) && length(legend.title) != 1) {
    stop("'legend.title' must be a single character string or NULL.", call. = FALSE)
  }

  if (!is.null(panel.labels)) {
    if (!is.character(panel.labels) || length(panel.labels) != length(var.names)) {
      stop("'panel.labels' must be a character vector with the same length as 'var.names'.", call. = FALSE)
    }
  }

  is_sf_data <- inherits(data, "sf")

  if (is_sf_data) {
    long_data <- reshape2::melt(
      data,
      measure.vars = var.names,
      variable.name = "prediction_var",
      value.name = "value"
    )
  } else {
    long_data <- reshape2::melt(
      data,
      id.vars = c("X", "Y"),
      measure.vars = var.names,
      variable.name = "prediction_var",
      value.name = "value"
    )
  }

  if (is.null(panel.labels)) {
    long_data$prediction_var <- factor(long_data$prediction_var, levels = var.names)
  } else {
    long_data$prediction_var <- factor(long_data$prediction_var, levels = var.names, labels = panel.labels)
  }

  if(!is.null(base.map)) {
    p <- ggplot2::ggplot(base.map) + ggplot2::geom_sf()
  } else{
    p <- ggplot2::ggplot()
  }

  if (is_sf_data) {
    p <- p +
      ggplot2::geom_sf(data = long_data, ggplot2::aes(color = value)) +
      ggplot2::scale_color_gradientn(colours = color.gradient, name = legend.title)
  } else {
    p <- p +
      ggplot2::geom_tile(data = long_data, ggplot2::aes(x = X, y = Y, fill = value)) +
      ggplot2::scale_fill_gradientn(colours = color.gradient, name = legend.title)
  }

  p <- p +
    ggplot2::labs(x = "Longitude", y = "Latitude") +
    ggplot2::scale_x_continuous(breaks = xaxis.breaks) +
    ggplot2::scale_y_continuous(breaks = yaxis.breaks) +
    ggplot2::facet_wrap(~ prediction_var,
                        strip.position = "top",
                        nrow = nrow) +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(
      axis.title.x = ggplot2::element_text(margin = ggplot2::margin(t = 10, b = 10)),
      axis.title.y = ggplot2::element_text(margin = ggplot2::margin(r = 10, l = 10)),
      panel.border = ggplot2::element_rect(color = "grey", fill = NA),
      strip.text = ggplot2::element_text(hjust = 0, vjust = 1)
    ) +
    ggspatial::annotation_north_arrow(location = "tl", height = grid::unit(0.6, "cm"), width = grid::unit(0.3, "cm")) +
    ggspatial::annotation_scale(location = "br", bar_cols = c("grey60", "white"))

  return(p)
}
