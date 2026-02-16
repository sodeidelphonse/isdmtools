# Methods for GeoDiagnostic objects

Methods for GeoDiagnostic objects

## Usage

``` r
# S3 method for class 'GeoDiagnostic'
print(x, ...)

# S3 method for class 'GeoDiagnostic'
plot(x, ...)
```

## Arguments

- x:

  A `GeoDiagnostic` object.

- ...:

  Additional arguments.

## Value

- `print`: Invisibly returns the original object.

- `plot`: Returns a `ggplot2` object.

## See also

Other diagnostic tools:
[`EnvDiagnostic-methods`](https://sodeidelphonse.github.io/isdmtools/reference/EnvDiagnostic-methods.md),
[`check_env_balance()`](https://sodeidelphonse.github.io/isdmtools/reference/check_env_balance.md),
[`check_folds()`](https://sodeidelphonse.github.io/isdmtools/reference/check_folds.md),
[`summarise_fold_diagnostics()`](https://sodeidelphonse.github.io/isdmtools/reference/summarise_fold_diagnostics.md)
