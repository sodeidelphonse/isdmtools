# Diagnostic for Spatial Folds Geometry

Internal function to evaluate the spatial structure of folds for blocked
cross-validation.

## Usage

``` r
check_spatial_geometry(
  data_all,
  fold_col = "folds_ids",
  rho = NULL,
  plot = TRUE
)
```

## Arguments

- data_all:

  An sf object containing pooled locations and fold IDs.

- fold_col:

  Character. Name of the fold ID column.

- rho:

  Numeric. Estimated spatial range (km).

- plot:

  Logical. If TRUE, generates a ggplot object.
