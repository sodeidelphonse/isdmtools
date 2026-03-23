# Replicate response data from the posterior predictive distribution.

The function generates replicated response data from the posterior
samples of parameters for use in model diagnostics. It currently
supports Poisson, Negative Binomial, and Quasi-Poisson families.

## Usage

``` r
simulate_replicates(
  mu_matrix,
  family = c("poisson", "nbinomial", "quasi-poisson"),
  dispersion = NULL,
  ...
)
```

## Arguments

- mu_matrix:

  A matrix of posterior samples of the expected count rate (mu).

- family:

  Character string defining the likelihood family (e.g., "poisson",
  "nbinomial").

- dispersion:

  Optional parameter (e.g., size 'k' for the negative binomial family).

- ...:

  Additional arguments.

## Value

A matrix containing the simulated predictive samples

## See also

Other PPC diagnostics:
[`compute_ppc_stats()`](https://sodeidelphonse.github.io/isdmtools/reference/compute_ppc_stats.md),
[`ppc_stats-methods`](https://sodeidelphonse.github.io/isdmtools/reference/ppc_stats-methods.md)

## Examples

``` r
if (FALSE) { # \dontrun{
#--- Create a toy matrix of posterior samples for expected counts (mu)
# 5 sites (rows) and 3 posterior draws (columns)
mu_samples <- matrix(c(
  1.2, 1.5, 1.1,
  5.0, 4.8, 5.2,
  0.5, 0.2, 0.4,
  12.1, 11.5, 13,
  2.5, 2.2, 2.8
), nrow = 5, ncol = 3, byrow = TRUE)

#--- Generate Poisson replicates
# Returns a 5x3 matrix of random integer counts
sim_pois <- simulate_replicates(mu_samples, family = "poisson")
print(sim_pois)

#--- Generate Negative Binomial replicates with high variance
# Requires a dispersion (size) parameter
sim_nb <- simulate_replicates(mu_samples,
  family = "nbinomial",
  dispersion = 0.8
)
print(sim_nb)

#--- Integration with DHARMa (if installed)
if (requireNamespace("DHARMa", quietly = TRUE)) {
  # Define a mock observed response
  y_obs <- c(1, 5, 0, 12, 3)

  # Create DHARMa object
  res_dharma <- DHARMa::createDHARMa(
    simulatedResponse = sim_pois,
    observedResponse = y_obs,
    fittedPredictedResponse = rowMeans(mu_samples),
    integerResponse = TRUE
  )

  # Plot residuals
  par(mfrow = c(1, 2))
  plotQQunif(res_dharma)
  #    testDispersion(res_dharma)
}
} # }
```
