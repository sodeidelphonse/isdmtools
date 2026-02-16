# Methods for EnvDiagnostic objects

Methods for EnvDiagnostic objects

## Usage

``` r
# S3 method for class 'EnvDiagnostic'
print(x, ...)

# S3 method for class 'EnvDiagnostic'
plot(x, ...)
```

## Arguments

- x:

  An `EnvDiagnostic` object.

- ...:

  Additional arguments.

## Value

- `print`: Invisibly returns the original object.

- `plot`: Returns a `ggplot2` object for covariates' density plots.

## See also

Other diagnostic tools:
[`GeoDiagnostic-methods`](https://sodeidelphonse.github.io/isdmtools/reference/GeoDiagnostic-methods.md),
[`check_env_balance()`](https://sodeidelphonse.github.io/isdmtools/reference/check_env_balance.md),
[`check_folds()`](https://sodeidelphonse.github.io/isdmtools/reference/check_folds.md),
[`summarise_fold_diagnostics()`](https://sodeidelphonse.github.io/isdmtools/reference/summarise_fold_diagnostics.md)
