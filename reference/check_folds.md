# Check Spatial Folds Independence

Evaluates the geometric properties of spatial folds to ensure spatial
independence for cross-validation. This function verifies if block sizes
and inter-block gaps are sufficient relative to the prior or model's
estimated spatial range (\\\rho\\).

## Usage

``` r
check_folds(object, ...)

# S3 method for class 'DataFolds'
check_folds(object, rho = NULL, plot = TRUE, ...)
```

## Arguments

- object:

  A `DataFolds` object created by
  [`create_folds()`](https://sodeidelphonse.github.io/isdmtools/reference/create_folds.md).

- ...:

  Additional arguments.

- rho:

  Numeric. Optional. The spatial range (km) estimated from the
  exploratory analysis (e.g.,
  [cv_spatial_autocor](https://rdrr.io/pkg/blockCV/man/cv_spatial_autocor.html))
  and used as the block size or the one estimated from the integrated
  model (e.g., the Matérn range parameter).

- plot:

  Logical. If `TRUE`, returns a diagnostic plot.

## Value

An object of class `GeoDiagnostic`.

## Details

The function assesses independence based on the minimum gap between
folds compared to the spatial range (\\\rho\\):

- **Contiguous**: Gap = 0. High risk of spatial leakage; observations in
  test folds are spatially correlated with training data.

- **Weakly Independent**: 0 \< Gap \< \\\rho\\. A physical gap exists,
  but correlation remains above 0.1.

- **Independent**: \\\rho \le\\ Gap \< \\2\rho\\. Spatial correlation is
  below 0.1 at the boundary; considered robust for most CV applications.

- **Strongly Independent**: Gap \\\ge 2\rho\\. Spatial correlation is
  effectively zero, providing the most rigorous test of model
  extrapolation.

## References

- Roberts DR, Bahn V, Ciuti S, Boyce MS, Elith J, Guillera-Arroita G,
  Hauenstein S, Lahoz-Monfort JJ, Schröder B, Thuiller W, et al.
  Cross-validation strategies for data with temporal, spatial,
  hierarchical, or phylogenetic structure. *Ecography* (2017)
  40:913–929.
  [doi:10.1111/ecog.02881](https://doi.org/10.1111/ecog.02881)

## See also

[`DataFolds-methods`](https://sodeidelphonse.github.io/isdmtools/reference/DataFolds-methods.md)
for interacting with `DataFolds` objects.

Other diagnostic tools:
[`EnvDiagnostic-methods`](https://sodeidelphonse.github.io/isdmtools/reference/EnvDiagnostic-methods.md),
[`GeoDiagnostic-methods`](https://sodeidelphonse.github.io/isdmtools/reference/GeoDiagnostic-methods.md),
[`check_env_balance()`](https://sodeidelphonse.github.io/isdmtools/reference/check_env_balance.md),
[`summarise_fold_diagnostics()`](https://sodeidelphonse.github.io/isdmtools/reference/summarise_fold_diagnostics.md)

## Examples

``` r
if (FALSE) { # \dontrun{
library(sf)
library(terra)
library(ggplot2)
library(isdmtools)

# Generate the data as a list of sf objects
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

ben_coords <- matrix(c(0, 6, 4, 6, 4, 13, 0, 13, 0, 6), ncol = 2, byrow = TRUE)
ben_sf <- st_sfc(st_polygon(list(ben_coords)), crs = 4326)
ben_sf <- st_sf(data.frame(name = "Benin"), ben_sf)

# Create Folds using create_folds()
folds <- create_folds(
  datasets_list,
  region_polygon = ben_sf,
  k = 5,
  cv_method = "cluster"
)

# Check Spatial Independence
# Assuming autocorrelation range (rho) is 150 km
spat_diag <- check_folds(folds, rho = 150, plot = TRUE)

# View results
print(spat_diag)
plot(spat_diag)
} # }
```
