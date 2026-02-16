# Calculate Niche Overlap (Schoener's D)

Computes the overlap between two distributions.

## Usage

``` r
calc_niche_overlap(x, y, n = 512)
```

## Arguments

- x:

  Numeric vector (e.g., predictions or spatial fold values).

- y:

  Numeric vector (e.g., observations or background values).

- n:

  Numeric. Number of points for density estimation. Default 512.

## Value

Numeric value between 0 and 1.

## Examples

``` r
if (FALSE) { # \dontrun{
library(isdmtools)

# Create two identical distributions (should have high overlap)
v1 <- rnorm(1000, mean = 10, sd = 1)
v2 <- rnorm(1000, mean = 10, sd = 1)
overlap_high <- calc_niche_overlap(v1, v2)

# Create two divergent distributions (should have low overlap)
v3 <- rnorm(1000, mean = 20, sd = 1)
overlap_low <- calc_niche_overlap(v1, v3)
} # }
```
