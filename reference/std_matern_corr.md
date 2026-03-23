# Standardise Matern covariance parameters and compute the pairwise correlation from LGCP models

Standardizes Matern covariance parameters across different modeling
engines (`INLA`, `spatstat`, `variofit`) and computes the pairwise
correlation function for log-Gaussian Cox process (LGCP) models.

## Usage

``` r
std_matern_corr(model, engine = c("kppm", "inla", "variofit"), r, nu = 1, ...)
```

## Arguments

- model:

  The fitted model object from `INLA`, `spatstat`, or `geoR`

- engine:

  Character string specifying the modeling engine utilised to estimate
  the covariance parameters: `"inla"`, `"kppm"`, or `"variofit"`.

- r:

  A numeric vector of distances at which to evaluate the covariance
  function.

- nu:

  The smoothing parameter of the process (it is \\\kappa\\ in `geoR`).

- ...:

  Additional arguments to pass on to internal function.

## Value

An object of class `std_matern_corr` (list) containing:

- `dist_vals`: The vector of distances (\\r\\) used for computation.

- `cov`: The computed Matern covariance: \\C(r)\\.

- `sigma2`: The standardized marginal variance (\\\sigma^2\\).

- `rho`: The standardized scale parameter (range).

- `nu`: The smoothness parameter.

- `engine`: The original modeling engine used.

For `"inla"` or `"kppm"` modeling engines, two additional elements are
included:

- `pair_cor`: The pairwise correlation function (PCF) for LGCP,
  calculated as \\g(r) = \exp(C(r))\\.

- `pair_cor_sc`: The standardized pairwise correlation function for
  cross-model comparison.

## Details

The function standardizes Matern covariance parameters and derives the
pairwise correlation function across different modeling engines for LGCP
models. The engines currently supported include those implemented in
`INLA/inlabru` (and related wrappers), `spatstat`, and `geoR`.

Note that the `variofit` class is included to derive the range and
`sigma2` parameters for comparison purposes; however, the exponential of
the resulting covariance has no interpretation in terms of pairwise
correlation. The same applies to covariance functions fitted to
geostatistical data in `INLA/inlabru`, which cannot be interpreted in
terms of pairwise correlation.

To facilitate visual comparison across different modeling engines, a
scaled version of the pairwise correlation function (`pair_cor_sc`) is
provided. This version applies a min-max normalization to the pairwise
correlation, mapping its values to a standardized range (typically \[0,
1\]). This ensures that the spatial decay and interaction strength can
be compared even when models have different marginal variances or
measurement scales.

## References

- Baddeley A, Rubak E, Turner R. Spatial point patterns: Methodology and
  applications with R. Boca Raton, FL: CHAPMAN & HALL CRC. (2015).

- Diggle PJ, Ribeiro PJ. Model-based Geostatistics. 1st ed. New York,
  NY: Springer. (2007).
  [doi:10.1007/978-0-387-48536-2](https://doi.org/10.1007/978-0-387-48536-2)

- Lindgren F, Rue H, Lindström J. An explicit link between Gaussian
  fields and Gaussian Markov random fields: the stochastic partial
  differential equation approach. *Journal of the Royal Statistical
  Society: Series B (Statistical Methodology)* (2011) 73:423–498.
  [doi:10.1111/j.1467-9868.2011.00777.x](https://doi.org/10.1111/j.1467-9868.2011.00777.x)

## See also

Other Matern covariance helpers:
[`plot.std_matern_corr()`](https://sodeidelphonse.github.io/isdmtools/reference/plot.std_matern_corr.md),
[`solve_practical_range()`](https://sodeidelphonse.github.io/isdmtools/reference/solve_practical_range.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Example with an INLA model
r_vec <- seq(0, 10, length.out = 100)
matern_results <- std_matern_corr(
  model = fit_inla,
  engine = "inla",
  r = r_vec,
  nu = 1
)

# Access the pair correlation function for LGCP
plot(matern_results$dist_vals, matern_results$cor,
  type = "l",
  ylab = "PCF", xlab = "Distance"
)
} # }
```
