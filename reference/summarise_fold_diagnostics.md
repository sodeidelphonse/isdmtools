# Summarise Fold Diagnostics

Combines geographic and environmental diagnostics into a single unified
report to evaluate the quality of a cross-validation scheme.

## Usage

``` r
summarise_fold_diagnostics(geo_diag, env_diag)

# S3 method for class 'FoldsSummary'
print(x, ...)
```

## Arguments

- geo_diag:

  A `GeoDiagnostic` object.

- env_diag:

  An `EnvDiagnostic` object.

- x:

  A `FoldsSummary` object.

- ...:

  Additional arguments

## Value

- `summarise_fold_diagnostics`: An object of class `FoldsSummary`, which
  inherits from `data.frame`.

- `print`: The `FoldsSummary` object invisibly.

## See also

[`DataFolds-methods`](https://sodeidelphonse.github.io/isdmtools/reference/DataFolds-methods.md)
for interacting with `DataFolds` objects.

Other diagnostic tools:
[`EnvDiagnostic-methods`](https://sodeidelphonse.github.io/isdmtools/reference/EnvDiagnostic-methods.md),
[`GeoDiagnostic-methods`](https://sodeidelphonse.github.io/isdmtools/reference/GeoDiagnostic-methods.md),
[`check_env_balance()`](https://sodeidelphonse.github.io/isdmtools/reference/check_env_balance.md),
[`check_folds()`](https://sodeidelphonse.github.io/isdmtools/reference/check_folds.md)

## Examples

``` r
if (FALSE) { # \dontrun{
library(sf)
library(terra)
library(ggplot2)
library(isdmtools)

# Generate points data
set.seed(42)
presence_data <- data.frame(
  x = runif(100, 0, 4),
  y = runif(100, 6, 13),
  site = rbinom(100, 1, 0.6)
) %>% st_as_sf(coords = c("x", "y"), crs = 4326)

count_data <- data.frame(
  x = runif(50, 0, 4),
  y = runif(50, 6, 13),
 count = rpois(50, 5)
) %>% st_as_sf(coords = c("x", "y"), crs = 4326)

datasets_list <- list(Presence = presence_data, Count = count_data)

# Environmental data
set.seed(42)
r <- rast(extent = c(0, 4, 6, 13), nrow=100, ncol=100, crs='epsg:4326')
r[] <- rnorm(ncell(r))
rtmp   <- r
rtmp[] <- runif(ncell(r), 5, 10)

r_stk <- c(r, rtmp + r)
names(r_stk) <- c("cov1", "cov2")

# Create Folds
folds <- create_folds(datasets_list, cv_method = "cluster")

# Spatial diagnostics
spat_diag <- check_folds(folds, plot = TRUE)

# Environmental diagnostics
env_diag <- suppressWarnings(check_env_balance(
  folds,
  covariates = r_stk,
  n_background = 5000)
)

# Combined diagnostics
sum_diag <- summarise_fold_diagnostics(spat_diag, env_diag)
print(sum_diag)
} # }
```
