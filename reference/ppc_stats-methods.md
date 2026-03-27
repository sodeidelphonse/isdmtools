# Posterior Predictive Check (PPC) Statistics

- `plot`: Visualizes the distribution of replicated T-statistics against
  the observed T-statistic using ggplot2.

- `print`: Display the Bayesian p-value and the dispersion status for
  the model fitted to the data.

## Usage

``` r
# S3 method for class 'ppc_stats'
print(x, ...)

# S3 method for class 'ppc_stats'
plot(x, ...)
```

## Arguments

- x:

  An object of class "ppc_stats".

- ...:

  Further arguments.

## Value

- `print`: Invisibly returns the original object.

- `plot`: Returns a `ggplot2` object.

## See also

Other PPC diagnostics:
[`compute_ppc_stats()`](https://sodeidelphonse.github.io/isdmtools/reference/compute_ppc_stats.md),
[`simulate_replicates()`](https://sodeidelphonse.github.io/isdmtools/reference/simulate_replicates.md)
