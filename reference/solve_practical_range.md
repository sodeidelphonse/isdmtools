# Solve the practical range of a Matérn covariance function

Harmonizes spatial range parameters from different R packages (`INLA`,
`geoR`, `spatstat`) into a standardized "Practical Range." This is the
distance at which the spatial correlation drops to a specific threshold
(default is 0.10).

## Usage

``` r
solve_practical_range(
  param_val,
  nu,
  sigma_sq = 1,
  thresh = 0.1,
  package = c("inla", "geor", "spatstat")
)
```

## Arguments

- param_val:

  Numeric. The parameter value from the model (\\\rho\\ for INLA,
  \\\phi\\ for geoR, or \\\alpha\\ for spatstat). It must be positive.

- nu:

  Numeric. The smoothness parameter. It must be positive. For 2D SPDE
  models in INLA (where alpha = 2), the default is `nu = 1`. For an
  exponential covariance, use `nu = 0.5`.

- sigma_sq:

  Numeric. The partial sill (or marginal variance). Defaults to 1 for
  correlation focus.

- thresh:

  Numeric. The target correlation threshold. Defaults to 0.1 (10%).

- package:

  Character. One of `"inla"`, `"geor"`, or `"spatstat"`.

## Value

A numeric value representing the practical range in the same geographic
units as the input model parameter.

## Details

Different packages use different parameterisations for the Matérn
covariance:

- **INLA/inlabru:** Estimates a value close to the INLA range parameter
  (where correlation is ~ 0.139). If `thresh = 0.139`, the input
  `param_val` is returned almost as is. If a 5% threshold
  (`thresh = 0.05`) is desired, the function adjusts the INLA range
  accordingly.

- **geoR:** Uses a scale parameter \\\phi\\. The practical range is
  solved numerically based on \\\phi\\ and the smoothness \\\nu\\.

- **spatstat:** Uses a scale parameter \\\alpha\\. The function aligns
  this with the INLA-style practical range.

This harmonization ensures that the `rho` value used in `isdmtools`
diagnostic functions is consistent, regardless of the modeling engine
used for the exploratory analysis.

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

## Examples

``` r
# Estimated phi = 10 km with exponential covariance in `geoR`
solve_practical_range(param_val = 10, nu = 0.5, thresh = 0.1, package = "geor")
#> [1] 23.02585

# Estimated alpha = 10 km with Matérn covariance in `spatstat`
solve_practical_range(param_val = 13.10, nu = 1.5, thresh = 0.1, package = "spatstat")
#> [1] 29.41908
```
