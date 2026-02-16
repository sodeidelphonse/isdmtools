# Fill NA cells with the nearest cells values

Function to impute a raster object for spatial modeling tools that
cannot handle missing values in covariates. The function fill in NA
cells with the nearest cells values using a moving window on missing
cells. It is an iterative version of the
[focal](https://rspatial.github.io/terra/reference/focal.html) function
of `terra` to handle incomplete filling due to isolated NA cells
surrounded by other NAs.

## Usage

``` r
fill_na_near(
  x,
  boundary = NULL,
  fun = mean,
  na.policy = "only",
  na.rm = TRUE,
  start.window = 1,
  ...
)
```

## Arguments

- x:

  A raster layer (`SpatRaster` or `RasterLayer`) with missing values to
  fill in. For multiple rasters (`SpatRaster` or `RasterStack`), you can
  use the combination of
  [`lapply()`](https://rdrr.io/r/base/lapply.html) and
  [`rast()`](https://rspatial.github.io/terra/reference/rast.html)
  functions with this function.

- boundary:

  A optional spatial polygon object (`spatVector` or `sf`) to be used to
  mask extra cells outside the study region. It must have the same
  coordinates reference system (CRS) with the input raster. Defaults to
  `NULL`.

- fun:

  A function to compute a value for a cell based on the values of its
  neighbors. The default is `mean`. This function must take a vector of
  values and return a single value (e.g., mean, modal, min or max) or
  multiple values (e.g., quantile).

- na.policy:

  Character. Specifies which cells to fill. Must be one of "all"
  (compute for all cells), "only" (only for cells that are NA) or "omit"
  (skip cells that are NA). The default `"only"` ensures that only cells
  that are NA in the input raster are filled.

- na.rm:

  Logical. If `TRUE`, NA values in the neighborhood are ignored when
  computing the focal function. The default is `TRUE`.

- start.window:

  A positive odd integer specifying the starting size of the square
  focal window. The window will be `(2w+1) x (2w+1)`. The default is
  `1`, which corresponds to a 3x3 start window.

- ...:

  Additional arguments passed to the internal
  [focal](https://rspatial.github.io/terra/reference/focal.html)
  function. This includes arguments like `expand`, `silent`, `filename`,
  etc. Note that a custom `w` (e.g., a weights matrix) cannot be passed
  via `...` as this function's logic is built around an iterative square
  window.

## Value

A `SpatRaster` object that has all of its NA cells filled in.

## Examples

``` r
if (FALSE) { # \dontrun{
library(terra)
# Create a sample raster with some NA values
r <- terra::rast(nrows = 10, ncols = 10, res = 1, xmin = 0, ymin = 0)
r[] <- 1:100
r[c(10, 25, 50, 75, 90)] <- NA

# Fill the NAs using the default mean function with 3x3 window
r_filled <- fill_na_near(r)

# Start the filling with a larger 5x5 window
r_filled_large_start <- fill_na_near(r, start.window = 3)

plot(r, main = "Original raster")
plot(r_filled, main = "Filled with mean (3x3 start)")
plot(r_filled_large_start, main = "Filled with mean (5x5 start)")
} # }
```
