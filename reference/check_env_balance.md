# Check Environmental Balance of Folds

Evaluates whether environmental covariates are well-balanced across
folds. It extracts values from a `SpatRaster` at point locations. For
large rasters, it uses a background sample to represent the study area.

## Usage

``` r
check_env_balance(object, ...)

# S3 method for class 'DataFolds'
check_env_balance(
  object,
  covariates,
  plot_type = c("density", "boxplot"),
  n_background = 10000,
  ...
)
```

## Arguments

- object:

  A `DataFolds` object.

- ...:

  Additional arguments passed to `sample_background`.

- covariates:

  A `SpatRaster` (terra) containing environmental layers. It must have
  the same coordinate reference system (CRS) as the sf objects used for
  blocking.

- plot_type:

  Character. Either "density" (default) or "boxplot".

- n_background:

  Numeric. Number of background points to sample for environmental space
  representation. Default 10,000.

## Value

An object of class `EnvDiagnostic`.

## Details

The function also calculates the environmental niche overlap using
Schoener's D metric (Schoener, 1968). This metric ranges from 0 (no
overlap) to 1 (identical niches).

- **Schoener's D**: Quantifies how well each fold represents the
  available environmental space (the background). A median value is
  reported across all folds for each covariate.

- **Interpretation**: Values \> 0.6 generally indicate that the folds
  are representative of the study area's environmental conditions. Low
  values suggest that cross-validation results may be biased because the
  model is being tested on environmental conditions it rarely
  encountered during training.

## See also

[`DataFolds-methods`](https://sodeidelphonse.github.io/isdmtools/reference/DataFolds-methods.md)
for interacting with `DataFolds` objects.

Other diagnostic tools:
[`EnvDiagnostic-methods`](https://sodeidelphonse.github.io/isdmtools/reference/EnvDiagnostic-methods.md),
[`GeoDiagnostic-methods`](https://sodeidelphonse.github.io/isdmtools/reference/GeoDiagnostic-methods.md),
[`check_folds()`](https://sodeidelphonse.github.io/isdmtools/reference/check_folds.md),
[`summarise_fold_diagnostics()`](https://sodeidelphonse.github.io/isdmtools/reference/summarise_fold_diagnostics.md)

## Examples

``` r
if (FALSE) { # \dontrun{
library(sf)
library(terra)
library(ggplot2)
library(isdmtools)

# Generate data as a list of sf objects
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

# a) Continuous covariates
r   <- rast(ben_sf, nrow = 100, ncol = 100, crs = 'epsg:4326')
r[] <- rnorm(ncell(r))
rtmp   <- r
rtmp[] <- runif(ncell(r), 5, 10)

r_stk <- c(r, rtmp + r)
names(r_stk) <- c("cov1", "cov2")

# Create Folds
folds <- create_folds(
  datasets_list,
  region_polygon = ben_sf,
  cv_method = "cluster"
)

# Check Environmental Representation
env_diag <- suppressWarnings(check_env_balance(
  folds,
  covariates = r_stk,
  n_background = 5000)
)

# View p-values in console
print(env_diag)

# View density plots
plot(env_diag)

# b) Mixture of continuous and categorical covariates
set.seed(42)
r_temp <- rast(extent = c(0, 4, 6, 13), res = 0.1, val = runif(2500, 15, 25))
r_land <- rast(extent = c(0, 4, 6, 13), res = 0.1, val = sample(1:3, 2500, TRUE))

# Set up land cover as a factor
levels(r_land) <- data.frame(ID = 1:3, cover = c("Forest", "Grass", "Urban"))
env_stack <- c(r_temp, r_land)
names(env_stack) <- c("temperature", "land_use")

# Run the diagnostic with 2000 cells over 2800 available
env_diag <- check_env_balance(
  folds,
  covariates = env_stack,
  n_background = 2000,
  plot_type = "boxplot"
)

# Inspect results
# 'temperature' will show a p-value and Schoener's D
# 'land_use' will show a p-value (Chi-sq) and Schoener_D as NA
print(env_diag)
plot(env_diag)
} # }
```
