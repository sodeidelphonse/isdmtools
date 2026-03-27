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
# 10 sites (rows) and 100 posterior draws (columns)
set.seed(123)
n_sites <- 10
n_draws <- 100

# Generating mu based on a hypothetical spatial trend
mu_base <- c(1.2, 5.0, 0.5, 12.1, 2.5, 3.1, 0.8, 7.4, 1.1, 4.2)
mu_samples <- matrix(
  rlnorm(n_sites * n_draws, meanlog = log(mu_base), sdlog = 0.1),
  nrow = n_sites, ncol = n_draws
)

#--- Generate Poisson replicates
# Returns a 10x100 matrix of random integer counts
sim_pois <- simulate_replicates(mu_samples, family = "poisson")
print(sim_pois)

#--- Generate Negative Binomial replicates with high variance
# Dispersion (size) = 1.5 indicates moderate overdispersion
sim_nb <- simulate_replicates(mu_samples,
                              family = "nbinomial",
                              dispersion = 1.5)
print(sim_nb)

#--- Integration with DHARMa (if installed)
if (requireNamespace("DHARMa", quietly = TRUE)) {
  # Observed response (10 observations)
  y_obs <- rpois(n_sites, mu_base)

  # Create DHARMa object using the 50 replicates
  res_dharma <- DHARMa::createDHARMa(
    simulatedResponse = sim_pois,
    observedResponse = y_obs,
    fittedPredictedResponse = rowMeans(mu_samples),
    integerResponse = TRUE
  )

  # Plot residuals and test dispersion
  par(mfrow = c(1, 2))
  plotQQunif(res_dharma)
  testDispersion(res_dharma)
 }
} # }
```
