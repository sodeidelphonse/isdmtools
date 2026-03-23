# Compute Posterior Predictive Check (PPC) Statistics

This function performs a Bayesian Posterior Predictive Check for count
data models. It compares the observed Pearson Chi-squared statistic
against a distribution of statistics calculated from data replicated
from the posterior predictive distribution.

## Usage

``` r
compute_ppc_stats(
  samples_matrix,
  data_response,
  family = c("poisson", "nbinomial", "quasi-poisson"),
  dispersion = NULL
)
```

## Arguments

- samples_matrix:

  A matrix of posterior samples of the expected count rate (mu). Rows
  correspond to observations, and columns correspond to Monte Carlo
  samples.

- data_response:

  A numeric vector of the observed count data.

- family:

  Character string defining the likelihood family: "poisson",
  "nbinomial", or "quasi-poisson".

- dispersion:

  Optional numeric value for the dispersion parameter (e.g., size 'k'
  for Negative Binomial). Defaults to 1 if NULL for "nbinomial" or
  "quasi-poisson".

## Value

A list of class "ppc_stats" containing:

- T_obs:

  The observed Pearson Chi-squared statistic.

- T_rep:

  A vector of Pearson Chi-squared statistics from replicated data.

- p_value:

  The Bayesian p-value.

- family:

  The likelihood family used.

## Details

The function calculates the Pearson Chi-squared statistic for the
observed data \\T\_{obs}\\ and for each replicated dataset \\T\_{rep}\\
generated from the posterior samples.

The statistic is defined as: \$\$T(y, \mu) = \sum\_{i=1}^{n}
\frac{(y_i - \mu_i)^2}{Var(y_i \| \mu_i)}\$\$

**Interpretation:**

- **Bayesian p-value:** Represents the probability that the replicated
  data is "more extreme" than the observed data.

  - \\p\_{B} \approx 0.5\\: Indicates an excellent model fit.

  - \\p\_{B} \< 0.05\\ or \\p\_{B} \> 0.95\\: Indicates a significant
    lack of fit.

- **Overdispersion:** If \\T\_{obs}\\ is consistently larger than the
  \\T\_{rep}\\ distribution (\\p \approx 1\\).

- **Underdispersion:** If \\T\_{obs}\\ is consistently smaller than the
  \\T\_{rep}\\ distribution (\\p \approx 0\\).

## See also

Other PPC diagnostics:
[`ppc_stats-methods`](https://sodeidelphonse.github.io/isdmtools/reference/ppc_stats-methods.md),
[`simulate_replicates()`](https://sodeidelphonse.github.io/isdmtools/reference/simulate_replicates.md)

## Examples

``` r
if (FALSE) { # \dontrun{
#--- Setup mock data: 10 observations, 100 posterior samples
set.seed(123)
n_obs <- 10
n_sim <- 100

# Expected counts (mu) from a hypothetical model
mu_samples <- matrix(rlnorm(n_obs * n_sim, meanlog = 1),
  nrow = n_obs, ncol = n_sim
)

# Observed counts (y) - let's assume the model is a good fit
y_obs <- rpois(n_obs, rowMeans(mu_samples))

#--- Compute PPC statistics for a Poisson family
ppc_results <- compute_ppc_stats(
  samples_matrix = mu_samples,
  data_response = y_obs,
  family = "poisson"
)

#--- View the summary (triggers the `print.ppc_stats` method)
print(ppc_results)

#--- Visualize the fit (triggers the `plot.ppc_stats` method)
plot(ppc_results)

#--- Example with Negative Binomial and DHARMa integration
if (requireNamespace("DHARMa", quietly = TRUE)) {
  # Generate replicates of response for DHARMa
  nb_sims <- simulate_replicates(mu_samples, family = "nbinomial", dispersion = 2)

  res_dharma <- createDHARMa(
    simulatedResponse = nb_sims,
    observedResponse = y_obs,
    fittedPredictedResponse = rowMeans(mu_samples)
  )
  plot(res_dharma)
}
} # }
```
