# Plot standardized Matern correlation or covariance

Visualization method for `std_matern_corr` objects using `ggplot2`. It
can display the Pairwise Correlation Function (PCF) for LGCP models or
the standard Matern covariance.

## Usage

``` r
# S3 method for class 'std_matern_corr'
plot(x, type = c("pair_cor", "covariance"), scaled = FALSE, ...)
```

## Arguments

- x:

  An object of class `std_matern_corr`.

- type:

  Character string. Either `"pair_cor"` (default for LGCP) or
  `"covariance"`.

- scaled:

  Logical. If `TRUE`, uses the min-max scaled PCF (`pair_cor_sc`) or the
  standardized covariance (\\C(r)/\sigma^2\\). Default is `FALSE`.

- ...:

  Additional arguments (currently ignored).

## Value

A `ggplot` object.

## See also

Other Matern covariance helpers:
[`solve_practical_range()`](https://sodeidelphonse.github.io/isdmtools/reference/solve_practical_range.md),
[`std_matern_corr()`](https://sodeidelphonse.github.io/isdmtools/reference/std_matern_corr.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Assuming 'cov_lgcp' is an object returned by std_matern_corr()
plot(cov_lgcp, type = "pair_cor", scaled = TRUE)

# To modify the plot
library(ggplot2)
p <- plot(cov_lgcp, type = "covariance")
p + theme_bw() + ggtitle("Spatial Decay")
} # }
```
