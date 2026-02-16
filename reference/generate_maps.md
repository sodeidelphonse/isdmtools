# Generate Multi-panel Maps from Spatial Model Predictions

This function creates a multi-panel map for visualizing multiple
prediction variables from a species distribution model or other spatial
model. It is designed to be flexible, handling both grid-based
(data.frame) and point-based (sf) spatial predictions.

## Usage

``` r
generate_maps(
  data,
  var_names = c("mean", "sd"),
  base_map = NULL,
  color_gradient = map.pal("viridis", 100),
  legend_title = NULL,
  panel_labels = NULL,
  nrow = NULL,
  xaxis_breaks = NULL,
  yaxis_breaks = NULL,
  annotate = TRUE
)
```

## Arguments

- data:

  A data frame, `sf` or `SpatRaster` object containing the prediction
  data. For grid-based data frame, it must contain columns named "x" and
  "y" representing pixels' coordinates.

- var_names:

  character. The vector of column names in `data` to be plotted on
  separate panels.

- base_map:

  An `sf` object to be plotted as a base layer underneath the prediction
  data (e.g., a background simple polygon). Defaults to `NULL`. Users
  can add additional vector geometries if needed using a ggplot2 syntax.

- color_gradient:

  A vector of valid colors to be used in the fill/color gradient.
  Defaults to `map.pal("viridis", 100)`.

- legend_title:

  character. The title of the color legend.

- panel_labels:

  character. An optional vector of labels for the facet panels. The
  order should correspond to `var_names`.

- nrow:

  integer. The number of rows for
  [`ggplot2::facet_wrap()`](https://ggplot2.tidyverse.org/reference/facet_wrap.html).
  Defaults to an optimal layout chosen by `ggplot2`.

- xaxis_breaks:

  A numeric vector specifying the breaks for the x-axis.

- yaxis_breaks:

  A numeric vector specifying the breaks for the y-axis.

- annotate:

  logical. If `TRUE`, add the north arrow and scale bar to the map.
  Defaults to `TRUE`

## Value

A `ggplot` object representing the multi-panel plot that can be
customized by the user.

## Details

The function internally reshapes the data from a wide format (with a
column for each prediction variable) to a long format suitable for
plotting with
[`ggplot2::facet_wrap()`](https://ggplot2.tidyverse.org/reference/facet_wrap.html).
It automatically selects the appropriate geometry
([`geom_tile()`](https://ggplot2.tidyverse.org/reference/geom_tile.html)
for grids and
[`geom_sf()`](https://ggplot2.tidyverse.org/reference/ggsf.html) for
points) and conditional scales. Users can also add to the map other
spatial vector layers or customize the plot using ggplot2 syntax if
needed.

## See also

Other prediction analyses:
[`prepare_predictions()`](https://sodeidelphonse.github.io/isdmtools/reference/prepare_predictions.md),
[`suitability_index()`](https://sodeidelphonse.github.io/isdmtools/reference/suitability_index.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# --- Example with grid-based data ---
# Simulate a data frame with coordinates and two prediction variables
grid_data <- expand.grid(x = 1:100, y = 1:100)
grid_data$mean <- rnorm(10000, mean = grid_data$x / 100, sd = 0.1)
grid_data$sd <- rgamma(10000, shape = 2, scale = 0.2)

# Simulate a boundary map (e.g., a simple polygon)
library(sf)
boundary <- st_sfc(st_polygon(list(cbind(c(0, 100, 100, 0, 0), c(0, 0, 100, 100, 0)))))
boundary_sf <- st_sf(data.frame(id = 1), geometry = boundary)

# Generate the map
generate_maps(
  data = grid_data,
  var_names = c("mean", "sd"),
  base_map = boundary_sf,
  color_gradient = c("white", "skyblue", "navy"),
  legend_title = "Prediction Value",
  panel_labels = c("Mean", "StDev"),
  nrow = 1
)

# --- Example with point-based data (sf) ---
# Simulate an sf object with point data
library(sf)
set.seed(123)
points_sf <- st_as_sf(grid_data[sample(1:10000, 1000), ],
                     coords = c("x", "y"), crs = 32631) # UTM CRS
points_sf <- prepare_predictions(points_sf)

# Generate the map
generate_maps(
  data = points_sf,
  var_names = c("mean", "sd"),
  base_map = boundary_sf,
  color_gradient = c("white", "orange", "red"),
  legend_title = "Prediction Value",
  panel_labels = c("Mean", "StDev"),
  nrow = 1
)
} # }
```
