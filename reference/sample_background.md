# Generate background points

A constructor function to generate background points for various
purposes (e.g. computing evaluation scores for presence-only data). It
exclude NA cells from the sample, and eventually the observed locations
if needed.

## Usage

``` r
sample_background(
  mask,
  points = NULL,
  n = 1000,
  method = "random",
  cells = FALSE,
  xy = TRUE,
  as.points = FALSE,
  na.rm = TRUE,
  ...
)
```

## Arguments

- mask:

  `SpatRaster` object to be used as mask (preferably, the predicted
  intensity or habitat suitability).

- points:

  Spatial points (`data.frame`, `sf` or `SpatVector` objects) that would
  be excluded from the background sample.

- n:

  integer. The number of pseudo-absence points to sample for
  presence-only data. The default is 1000.

- method:

  character. The sampling technique to select pixels from the raster
  mask (see
  [spatSample](https://rspatial.github.io/terra/reference/sample.html)).
  It defaults to `random`.

- cells:

  logical. If `TRUE`, sampled cells numbers will be returned. The
  default is `FALSE`.

- xy:

  logical. If `TRUE`, the locations of sampled cells will be returned.
  The default is `TRUE`.

- as.points:

  logical. If `TRUE`, spatial points object will be returned. The
  default is `FALSE`.

- na.rm:

  logical. If `TRUE`, NA values will be excluded from the raster mask.
  It defaults to `TRUE`.

- ...:

  Additional arguments passed on to the internal
  [spatSample](https://rspatial.github.io/terra/reference/sample.html)
  function.

## Value

An S3 object with class `BackgroundPoints`, containing the modified
`SpatRaster` object and the generated background points.

## See also

Other BackgroundPoints methods:
[`BackgroundPoints-methods`](https://sodeidelphonse.github.io/isdmtools/reference/BackgroundPoints-methods.md)

## Examples

``` r
if (FALSE) { # \dontrun{
library(terra)
set.seed(123)
r  <- rast(nrows = 100, ncols = 100, xmin = 0, xmax = 10, ymin = 0, ymax = 10)
terra::values(r) <- runif(ncell(r))
pts <- spatSample(r, size = 100, xy = TRUE, values = FALSE)

# Requesting few points with their x and y coordinates
set.seed(235)
bg_sample1 <- sample_background(r, points = pts, n = 500, xy = TRUE, cells = FALSE)
plot(bg_sample1)
print(bg_sample1)

# Requesting points more than available non-NA cells
bg_sample2 <- sample_background(r, points = pts, n = 10000, xy =TRUE, cells = FALSE)
dim(bg_sample2$bg)
plot(bg_sample2)
} # }
```
