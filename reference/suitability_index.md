# Compute a unified suitability index from integrated spatial model predictions.

This function converts the linear predictor (`eta`) from a fitted
integrated spatial model into a unified suitability index, which can be
interpreted as a probability of a species presence. It also calculates
the intensity(rate) or expected count for the count data, depending on
whether the model has offset or not.

## Usage

``` r
suitability_index(
  x,
  post_stat = "mean",
  output_format = c("prob", "response", "linear"),
  response_type = c("joint.po", "count.pa", "po", "count", "pa"),
  has_offset = FALSE,
  scale_independent = FALSE,
  projection = NULL,
  ...
)
```

## Arguments

- x:

  A `data.frame` containing `x` and `y` coordinates and the column(s) of
  the predicted linear predictor variables (e.g., mean, standard
  deviation and quantiles) or `SpatRaster`. It can typically be a
  standardized grid-based output from a
  [prepare_predictions](https://sodeidelphonse.github.io/isdmtools/reference/prepare_predictions.md)
  call to various classes of spatial prediction on a linear scale, e.g.
  from the `PointedSDMs` or `inlabru` packages.

- post_stat:

  character. A vector specifying the column or layer name(s) to use for
  extracting the model predictions. Defaults to "mean".

- output_format:

  character. The desired output format and must be one of "prob"
  (probability-based suitability index), "response" (expected count or
  rate) or "linear" (linear predictor scale).

- response_type:

  character. The type of response data the model was fitted with. Must
  be one of "joint.po" (joint model including presence-only data),
  "count.pa" (joint model with count and presence-absence), "po" (single
  presence-only model), "count" (count model), or "pa" (presence-absence
  model).

- has_offset:

  logical. For `count.pa`, `count` and `pa` models, this should be
  `TRUE` if the linear predictor includes an explicit area offset. This
  argument is not used for "po" or "joint.po" models. Defaults to
  `FALSE`.

- scale_independent:

  logical. If `TRUE`, the scaling factor is set to 1, making the
  suitability index independent of the grid cell size. Defaults to
  `FALSE`.

- projection:

  character. The coordinate reference system (CRS) for the output
  raster. Defaults to `NULL`. If `NULL` and `x` is a data.frame, an
  empty CRS is assigned to prevent errors.

- ...:

  Additional arguments passed to
  [rast](https://rspatial.github.io/terra/reference/rast.html).

## Value

A SpatRaster object representing the requested output format.

## Details

This function implements a unified framework for converting model
predictions to a comparable suitability index. The method relies on the
*Inhomogeneous Poisson Process (IPP)* theory, where the linear predictor
`eta` is related to the probability of presence via the inverse
complementary log-log link as follows: \\p(presence) = 1 - exp(-scaling
\times exp(eta))\\.

The `scaling` factor is determined by the `response_type` and the
`has_offset` arguments:

- For PO models (single or part of a joint model), `eta` is always a
  log-rate, and the `scaling` is set to the cell area. It is ignored if
  `scale_independent` is set to TRUE.

- For Count or PA models, if `has_offset = TRUE`, `eta` is a log-rate,
  and `scaling` is the cell area.

- In all other cases (`has_offset = FALSE`), `eta` is treated as a log
  of the expected count (count) or cloglog of probability (PA), and
  `scaling` is set to `1`.

If the raster is in a geographic coordinate system (longlat), the area
is calculated in \\km^2\\ using
[cellSize](https://rspatial.github.io/terra/reference/cellSize.html).
For projected systems, the area is the product of the resolutions (e.g.,
\\km^2\\ if units are in km).

## References

Dorazio RM. Accounting for imperfect detection and survey bias in
statistical analysis of presence-only data. *Global Ecology and
Biogeography* (2014) 23:1472–1484.
[doi:10.1111/geb.12216](https://doi.org/10.1111/geb.12216)

Fithian W, Elith J, Hastie T, Keith DA. Bias correction in species
distribution models: pooling survey and collection data for multiple
species. *Methods in Ecology and Evolution* (2015) 6:424–438.
[doi:10.1111/2041-210X.12242](https://doi.org/10.1111/2041-210X.12242)

Phillips SJ, Anderson RP, Dudík M, et al Opening the black box: an
open‐source release of Maxent. *Ecography* (2017) 40:887–893.
[doi:10.1111/ecog.03049](https://doi.org/10.1111/ecog.03049)

## See also

Other prediction analyses:
[`generate_maps()`](https://sodeidelphonse.github.io/isdmtools/reference/generate_maps.md),
[`prepare_predictions()`](https://sodeidelphonse.github.io/isdmtools/reference/prepare_predictions.md)

## Examples

``` r
if (FALSE) { # \dontrun{
library(terra)
# Simulate a sample data.frame with x, y, and a linear predictor
set.seed(42)
x <- expand.grid(x = seq(0, 50, 1), y = seq(0, 50, 1))
x$eta <- rnorm(nrow(x), mean = 0, sd = 1)

# Create a simple raster object
x_rast <- rast(x)

# Generate a suitability index for a Presence-Absence model
pa_probability <- suitability_index(
  x_rast,
  post_stat = "eta",
  response_type = "pa"
)
plot(pa_probability)

# Create a binary map using a fixed threshold
binary_map <- app(pa_probability, function(x) ifelse(x < 0.5, 0, 1))
plot(binary_map)

# Generate an expected mean (assume "eta" is from a count model)
expected_mean <- suitability_index(
  x_rast,
  post_stat = "eta",
  response_type = "count",
  output_format = "response"
)
plot(expected_mean)
} # }
```
