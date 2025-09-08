
#----------------------------------------------------------
#--- Multi-panel plot for spatial data visualization
#----------------------------------------------------------

#' Generate Multi-panel Maps from Spatial Model Predictions
#'
#' This function creates a multi-panel map for visualizing multiple prediction
#' variables from a species distribution model or other spatial model. It is
#' designed to be flexible, handling both grid-based (data.frame) and
#' point-based (sf) predictions.
#'
#' The function internally reshapes the data from a wide format (with a column
#' for each prediction variable) to a long format suitable for plotting with
#' `ggplot2::facet_wrap()`. It automatically selects the appropriate
#' geometry (`geom_tile` for grids, `geom_sf` for points) and conditional scales.
#'
#' @param data_to_plot A data frame, `sf` or `SpatRaster` object containing the prediction data.
#'   For grid-based data frame, it must contain columns named "X" and "Y" representing pixels' coordinates.
#' @param vars_to_plot A character vector of column names in `data_to_plot` to be
#'   plotted on separate panels.
#' @param base_map An `sf` object to be plotted as a base layer underneath the
#'   prediction data (e.g., a background simple polygon). Defaults to `NULL`.
#' @param layer_map An optional `sf` object representing the additional layer (e.g. line or point geometry).
#'   It is plotted on top of the prediction maps for context. Defaults to `NULL`.
#' @param color_gradient A vector of colors to be used in the fill/color gradient.
#'   Defaults to `rainbow(5)`.
#' @param legend_title A character string for the title of the color legend.
#' @param panel_labels An optional character vector of labels for the facet panels.
#'   The order should correspond to `vars_to_plot`.
#' @param nrow The number of rows for `ggplot2::facet_wrap()`. Defaults to an
#'   optimal layout chosen by `ggplot2`.
#' @param x_axis_breaks A numeric vector specifying the breaks for the x-axis.
#' @param y_axis_breaks A numeric vector specifying the breaks for the y-axis.
#'
#' @return A `ggplot` object representing the multi-panel map.
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
#'   data_to_plot = grid_data,
#'   vars_to_plot = c("mean", "sd"),
#'   base_map = boundary_sf,
#'   layer_map = boundary_sf,
#'   color_gradient = c("white", "skyblue", "navy"),
#'   legend_title = "Prediction Value",
#'   panel_labels = c("Mean", "Standard Deviation"),
#'   nrow = 1
#' )
#'
#' # --- Example with point-based data (sf) ---
#' # Simulate an sf object with point data
#' library(sf)
#' set.seed(123)
#' points_sf <- st_as_sf(grid_data[sample(1:10000, 1000), ],
#'                      coords = c("X", "Y"), crs = 32631) # UTM CRS
#'
#' # Generate the map
#' generate_maps(
#'   data_to_plot = points_sf,
#'   vars_to_plot = c("mean", "sd"),
#'   base_map = boundary_sf,
#'   layer_map = boundary_sf,
#'   color_gradient = c("white", "orange", "red"),
#'   legend_title = "Prediction Value",
#'   panel_labels = c("Mean", "StDev"),
#'   nrow = 2
#' )
#' }
#'
#' @importFrom ggplot2 ggplot geom_sf geom_tile aes scale_fill_gradientn scale_color_gradientn
#' @importFrom ggplot2 facet_wrap labs theme element_blank element_text element_rect
#' @importFrom reshape2 melt
#' @importFrom sf st_as_sf
#' @importFrom grDevices colorRampPalette rainbow terrain.colors heat.colors topo.colors
#' @importFrom grid unit
#'
generate_maps <- function(data_to_plot,
                          vars_to_plot = c("mean", "sd"),
                          base_map = NULL,
                          layer_map = NULL,
                          color_gradient = rainbow(5),
                          legend_title = NULL,
                          panel_labels = NULL,
                          nrow = NULL,
                          x_axis_breaks = NULL,
                          y_axis_breaks = NULL) {

  if (!is.data.frame(data_to_plot) && !inherits(data_to_plot, "sf") && !inherits(data_to_plot, "SpatRaster")) {
    stop("'data_to_plot' must be a data.frame, an sf object, or a SpatRaster.", call. = FALSE)
  }

  if(inherits(data_to_plot, "SpatRaster")) {
    data_to_plot <- as.data.frame(data_to_plot, xy = TRUE)
    data_to_plot <- data_to_plot %>% dplyr::rename(X = x, Y = y)
  }

  if (is.data.frame(data_to_plot) && !all(c("X", "Y") %in% names(data_to_plot))) {
    stop("'data_to_plot' is a data.frame but is missing 'X' and 'Y' coordinate columns.", call. = FALSE)
  }

  if (!is.character(vars_to_plot) || length(vars_to_plot) == 0) {
    stop("'vars_to_plot' must be a non-empty character vector of columns names.", call. = FALSE)
  }

  if (!all(vars_to_plot %in% names(data_to_plot))) {
    missing_vars <- setdiff(vars_to_plot, names(data_to_plot))
    stop(paste("The following variables in 'vars_to_plot' were not found in 'data_to_plot':",
               paste(missing_vars, collapse = ", ")), call. = FALSE)
  }

  if (!inherits(base_map, "sf")) {
    stop("'base_map' must be an sf object.", call. = FALSE)
  }
  if (!all(sf::st_geometry_type(base_map) %in% c("POLYGON", "MULTIPOLYGON"))) {
    stop("'base_map' must have polygon or multipolygon geometry.", call. = FALSE)
  }

  if (!is.null(layer_map) && !inherits(layer_map, "sf")) {
    stop("'layer_map' must be an sf object or NULL.", call. = FALSE)
  }

  if (!is.character(color_gradient) || length(color_gradient) < 2) {
    stop("'color_gradient' must be a character vector of at least two colors.", call. = FALSE)
  }

  if (!is.null(legend_title) && is.character(legend_title) && length(legend_title) != 1) {
    stop("'legend_title' must be a single character string or NULL.", call. = FALSE)
  }

  if (!is.null(panel_labels)) {
    if (!is.character(panel_labels) || length(panel_labels) != length(vars_to_plot)) {
      stop("'panel_labels' must be a character vector with the same length as 'vars_to_plot'.", call. = FALSE)
    }
  }

  # sf vs data.frame input
  is_sf_data <- inherits(data_to_plot, "sf")

  if (is_sf_data) {
    long_data <- reshape2::melt(
      data_to_plot,
      measure.vars = vars_to_plot,
      variable.name = "prediction_var",
      value.name = "value"
    )
  } else {
    long_data <- reshape2::melt(
      data_to_plot,
      id.vars = c("X", "Y"),
      measure.vars = vars_to_plot,
      variable.name = "prediction_var",
      value.name = "value"
    )
  }

  if (is.null(panel_labels)) {
    long_data$prediction_var <- factor(long_data$prediction_var, levels = vars_to_plot)
  } else {
    long_data$prediction_var <- factor(long_data$prediction_var, levels = vars_to_plot, labels = panel_labels)
  }

  # Add the main data layer and the conditional scale
  p <- ggplot2::ggplot(base_map) + ggplot2::geom_sf()

  if (is_sf_data) {
    p <- p +
      ggplot2::geom_sf(data = long_data, ggplot2::aes(color = value)) +
      ggplot2::scale_color_gradientn(colours = color_gradient, name = legend_title)
  } else {
    p <- p +
      ggplot2::geom_tile(data = long_data, ggplot2::aes(x = X, y = Y, fill = value)) +
      ggplot2::scale_fill_gradientn(colours = color_gradient, name = legend_title)
  }

  # Add the common layers for all panels
  geom_types <- sf::st_geometry_type(layer_map)

  if (all(geom_types %in% c("POINT", "MULTIPOINT"))) {
    p <- p + ggplot2::geom_sf(data = layer_map, fill = NA, color = "white", size = 1.5)
  } else if (all(geom_types %in% c("LINESTRING", "MULTILINESTRING"))) {
    p <- p + ggplot2::geom_sf(data = layer_map, fill = NA, color = "grey20", linewidth = 0.3)
  } else {
    warning(sprintf("`layer_map` has unhandled geometry type(s): %s. This layer is ignored.", paste(unique(geom_types), collapse = ", ")), call. = FALSE)
  }

  p <- p +
    ggplot2::labs(x = "Longitude", y = "Latitude") +
    ggplot2::scale_x_continuous(breaks = x_axis_breaks) +
    ggplot2::scale_y_continuous(breaks = y_axis_breaks) +
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
