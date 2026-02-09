
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
#'   For grid-based data frame, it must contain columns named "x" and "y" representing pixels' coordinates.
#' @param var_names character. The vector of column names in `data` to be
#'   plotted on separate panels.
#' @param base_map An `sf` object to be plotted as a base layer underneath the
#'   prediction data (e.g., a background simple polygon). Defaults to `NULL`.
#'   Users can add additional vector geometries if needed using a ggplot2 syntax.
#' @param color_gradient A vector of valid colors to be used in the fill/color gradient.
#'   Defaults to `map.pal("viridis", 100)`.
#' @param legend_title character. The title of the color legend.
#' @param panel_labels character. An optional vector of labels for the facet panels.
#'   The order should correspond to `var_names`.
#' @param nrow integer. The number of rows for `ggplot2::facet_wrap()`. Defaults to an
#'   optimal layout chosen by `ggplot2`.
#' @param xaxis_breaks A numeric vector specifying the breaks for the x-axis.
#' @param yaxis_breaks A numeric vector specifying the breaks for the y-axis.
#' @param annotate logical. If `TRUE`, add the north arrow and scale bar to the map. Defaults to `TRUE`
#'
#' @return A `ggplot` object representing the multi-panel plot that can be customized by the user.
#' @export
#' @examples
#' \dontrun{
#' # --- Example with grid-based data ---
#' # Simulate a data frame with coordinates and two prediction variables
#' grid_data <- expand.grid(x = 1:100, y = 1:100)
#' grid_data$mean <- rnorm(10000, mean = grid_data$x / 100, sd = 0.1)
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
#'   var_names = c("mean", "sd"),
#'   base_map = boundary_sf,
#'   color_gradient = c("white", "skyblue", "navy"),
#'   legend_title = "Prediction Value",
#'   panel_labels = c("Mean", "StDev"),
#'   nrow = 1
#' )
#'
#' # --- Example with point-based data (sf) ---
#' # Simulate an sf object with point data
#' library(sf)
#' set.seed(123)
#' points_sf <- st_as_sf(grid_data[sample(1:10000, 1000), ],
#'                      coords = c("x", "y"), crs = 32631) # UTM CRS
#' points_sf <- prepare_predictions(points_sf)
#'
#' # Generate the map
#' generate_maps(
#'   data = points_sf,
#'   var_names = c("mean", "sd"),
#'   base_map = boundary_sf,
#'   color_gradient = c("white", "orange", "red"),
#'   legend_title = "Prediction Value",
#'   panel_labels = c("Mean", "StDev"),
#'   nrow = 1
#' )
#' }
#'
#' @importFrom terra map.pal
#' @importFrom grid unit
#' @family prediction analyses
#'
generate_maps <- function(data,
                          var_names = c("mean", "sd"),
                          base_map = NULL,
                          color_gradient = map.pal("viridis", 100),
                          legend_title = NULL,
                          panel_labels = NULL,
                          nrow = NULL,
                          xaxis_breaks = NULL,
                          yaxis_breaks = NULL,
                          annotate = TRUE) {

  if (!inherits(data, c("data.frame", "sf", "SpatRaster"))) {
    stop("'data' must be a data.frame, an sf object, or a SpatRaster.", call. = FALSE)
  }

  if(inherits(data, "SpatRaster")) {
    data <- as.data.frame(data, xy = TRUE)
  }
  if (is.data.frame(data) && !inherits(data, "sf") && !all(c("x", "y") %in% names(data))) {
    stop("'data' is a data.frame but is missing 'x' and 'y' coordinate columns.", call. = FALSE)
  }

  if (!is.character(var_names) || length(var_names) == 0) {
    stop("'var_names' must be a non-empty character vector of columns names.", call. = FALSE)
  }

  if (!all(var_names %in% names(data))) {
    missing_vars <- setdiff(var_names, names(data))
    stop(paste("The following variables in 'var_names' were not found in 'data':",
               paste(missing_vars, collapse = ", ")), call. = FALSE)
  }

  if(!is.null(base_map)) {
    if (!inherits(base_map, "sf")) {
      stop("'base_map' must be an sf object.", call. = FALSE)
    }
    if (!all(sf::st_geometry_type(base_map) %in% c("POLYGON", "MULTIPOLYGON"))) {
      stop("'base_map' must have polygon or multipolygon geometry.", call. = FALSE)
    }
  }

  if (!is.character(color_gradient) || length(color_gradient) < 2) {
    stop("'color_gradient' must be a character vector of at least two colors.", call. = FALSE)
  }

  if (!is.null(legend_title) && is.character(legend_title) && length(legend_title) != 1) {
    stop("'legend_title' must be a single character string or NULL.", call. = FALSE)
  }

  if (!is.null(panel_labels)) {
    if (!is.character(panel_labels) || length(panel_labels) != length(var_names)) {
      stop("'panel_labels' must be a character vector with the same length as 'var_names'.", call. = FALSE)
    }
  }

  is_sf_data <- inherits(data, "sf")
  multi_panel <- length(var_names) > 1

  if(multi_panel) {
    long_data_list <- lapply(var_names, function(vn) {
      if (is_sf_data) {
        tmp <- data[, vn]
        names(tmp)[names(tmp) == vn] <- "value"
        tmp$prediction_var <- vn
        return(tmp)
      } else {
        tmp <- data[, c("x", "y", vn)]
        names(tmp)[names(tmp) == vn] <- "value"
        tmp$prediction_var <- vn
        return(tmp)
      }
    })

    long_data <- dplyr::bind_rows(long_data_list)

    if (is.null(panel_labels)) {
      long_data$prediction_var <- factor(long_data$prediction_var, levels = var_names)
    } else {
      long_data$prediction_var <- factor(long_data$prediction_var, levels = var_names, labels = panel_labels)
    }

  } else {
    long_data <- data
    fill_col  <- var_names
  }

  p <- ggplot2::ggplot(long_data)

  if (!is.null(base_map)) {
    if (is_sf_data && sf::st_crs(data) != sf::st_crs(base_map)) {
      base_map_to_plot <- sf::st_transform(base_map, sf::st_crs(data))
    } else {
      base_map_to_plot <- base_map
    }
    p <- p + ggplot2::geom_sf(data = base_map_to_plot, fill = NA, color = "grey20")
  }

  fill_col <- if(multi_panel) "value" else var_names

  if (is_sf_data) {
    p <- p +
      ggplot2::geom_sf(ggplot2::aes(color = .data[[fill_col]])) +
      ggplot2::scale_color_gradientn(colours = color_gradient, name = legend_title)
  } else {
    p <- p +
      ggplot2::geom_tile(ggplot2::aes(x = .data$x, y = .data$y, fill = .data[[fill_col]])) +
      ggplot2::scale_fill_gradientn(colours = color_gradient, name = legend_title)
  }

  p <- p +
    ggplot2::labs(x = "Longitude", y = "Latitude") +
    ggplot2::scale_x_continuous(breaks = xaxis_breaks) +
    ggplot2::scale_y_continuous(breaks = yaxis_breaks) +
    ggplot2::theme_bw(base_size = 12)

  # Conditional Faceting
  if (multi_panel) {
    p <- p + ggplot2::facet_wrap(~ prediction_var, strip.position = "top", nrow = nrow)
  } else {
    plot_title <- if(!is.null(panel_labels)) panel_labels else var_names
    p <- p + ggplot2::labs(title = plot_title)
  }

  p <- p + ggplot2::theme(
      axis.title.x = ggplot2::element_text(margin = ggplot2::margin(t = 10, b = 10)),
      axis.title.y = ggplot2::element_text(margin = ggplot2::margin(r = 10, l = 10)),
      panel.border = ggplot2::element_rect(color = "grey", fill = NA),
      strip.text = ggplot2::element_text(hjust = 0, vjust = 1),
      plot.title = ggplot2::element_text(size = 13, face = "bold", margin = ggplot2::margin(b = 10))
    )

  if (annotate && requireNamespace("ggspatial", quietly = TRUE)) {
    p <- p +
      ggspatial::annotation_north_arrow(location = "tl", height = grid::unit(0.6, "cm"), width = grid::unit(0.3, "cm")) +
      ggspatial::annotation_scale(location = "br", bar_cols = c("grey60", "white"))
  }

  return(p)
}
