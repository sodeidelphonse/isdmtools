# Methods for `BackgroundPoints` objects

- `plot`: Visualizes the background points generated in
  `BackgroundPoints` object. The plot shows cells with NA values or
  those of the locations excluded from the sample (white color) if the
  `points` argument is provided to `sample_background`. The background
  points generated are colored red.

- `print`: Display few points (first and last ones) in the R session.

## Usage

``` r
# S3 method for class 'BackgroundPoints'
plot(x, ...)

# S3 method for class 'BackgroundPoints'
print(x, ...)
```

## Arguments

- x:

  A `BackgroundPoints` S3 object.

- ...:

  Additional arguments (not used by this method).

## Value

Invisibly returns the original object.

## See also

Other BackgroundPoints methods:
[`sample_background()`](https://sodeidelphonse.github.io/isdmtools/reference/sample_background.md)

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
