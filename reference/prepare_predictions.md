# Obtain a formatted output from spatial predictions.

Function to transform a prediction data from various spatial models
(e.g., `inlabru`, `PointedSDMs` or `GLMs` tools) into an sf object if it
is from points predictions (e.g.,
[fm_vertices](https://inlabru-org.github.io/fmesher/reference/fm_vertices.html))
or a data.frame with the corresponding locations if it is from pixel
grids (see
[fm_pixels](https://inlabru-org.github.io/fmesher/reference/fm_pixels.html)).
Spatial prediction data can also be obtained across a given region using
`expand.grid(x, y)` function, where 'x' and 'y' are geographical
coordinates of the grids locations.

## Usage

``` r
prepare_predictions(prediction_data, base_map = NULL)
```

## Arguments

- prediction_data:

  A model prediction which may be either an `sf` or a `data.frame`
  object or a raw prediction from the `inlabru-like` models. The
  prediction can be on the response or linear predictor scale, depending
  on whether the output is for a model evaluation or visualization.

- base_map:

  An `sf` polygon having the same `crs` with the spatial locations used
  for predictions.

## Value

A `data.frame` for grid-based predictions or `sf` object for point-based
predictions.

## See also

Other prediction analyses:
[`generate_maps()`](https://sodeidelphonse.github.io/isdmtools/reference/generate_maps.md),
[`suitability_index()`](https://sodeidelphonse.github.io/isdmtools/reference/suitability_index.md)

## Examples

``` r
if (FALSE) { # \dontrun{
library(dplyr)
library(sf)
set.seed(42)

# Simulate the prediction data
grid_df <- expand.grid(x = 0:50, y = 0:50)
grid_df <- grid_df %>%
  mutate(mu = (x + y) / 10 + rnorm(nrow(grid_df))) %>%
  mutate(sd = runif(nrow(grid_df), 0.5, 1.5))
head(grid_df)

# Convert the grid into a SpatRaster
grid_r <- terra::rast(grid_df, crs = "epsg:4326")

# a) A standard data.frame returns the same object
field_pp1 <- prepare_predictions(grid_df)
class(field_pp1)

# b) A grid-based object returns a data.frame
grid_sf <- st_as_sf(grid_df, coords = c("x", "y"), crs = "epsg:4326")
class(grid_sf) <- c("bru_prediction", "sf", "data.frame")

field_pp2 <- prepare_predictions(grid_sf)
print(class(field_pp2))

# c) A point-based prediction returns the original class
if(require("fmesher", quietly = TRUE)) {
  bnd  <- fm_nonconvex_hull(as.matrix(grid_df[, c("x", "y")]), convex = -0.10)
  mesh <- fm_mesh_2d(boundary = bnd, max.edge = c(3, 30), crs = "epsg:4326")
  vt <- fm_vertices(mesh, format = "sf")

# An inlabru-like prediction at mesh vertices (extended areas are imputed)
 sampled_vals <- terra::extract(grid_r, vt)
 sim_field <- vt %>%
   mutate(mean = dplyr::coalesce(sampled_vals$mu, mean(grid_df$mu, na.rm = TRUE)),
         sd    = dplyr::coalesce(sampled_vals$sd, mean(grid_df$sd, na.rm = TRUE)),
         q0.025 = mean - 1.96 * sd,
         q0.5   = mean,
         q0.975 = mean + 1.96 * sd,
         median = mean)
 class(sim_field) <- c("bru_prediction", "sf", "data.frame")

 field_pp3  <- prepare_predictions(sim_field)
 print(class(field_pp3))
 }
} # }
```
