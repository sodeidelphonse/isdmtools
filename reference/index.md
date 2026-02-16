# Package index

## Folds Creation & Management

Core functions for spatial partitioning.

- [`create_folds()`](https://sodeidelphonse.github.io/isdmtools/reference/create_folds.md)
  [`bind_datasets()`](https://sodeidelphonse.github.io/isdmtools/reference/create_folds.md)
  : Create block cross-validation folds for multisource spatial datasets

- [`extract_fold()`](https://sodeidelphonse.github.io/isdmtools/reference/DataFolds-methods.md)
  [`plot(`*`<DataFolds>`*`)`](https://sodeidelphonse.github.io/isdmtools/reference/DataFolds-methods.md)
  [`print(`*`<DataFolds>`*`)`](https://sodeidelphonse.github.io/isdmtools/reference/DataFolds-methods.md)
  [`summary(`*`<DataFolds>`*`)`](https://sodeidelphonse.github.io/isdmtools/reference/DataFolds-methods.md)
  [`autoplot(`*`<DataFolds>`*`)`](https://sodeidelphonse.github.io/isdmtools/reference/DataFolds-methods.md)
  :

  Methods for manipulating data from `DataFolds` object.

## Spatial Folds Diagnostics

Evaluate the balance and quality of spatial folds.

- [`check_folds()`](https://sodeidelphonse.github.io/isdmtools/reference/check_folds.md)
  : Check Spatial Folds Independence
- [`check_env_balance()`](https://sodeidelphonse.github.io/isdmtools/reference/check_env_balance.md)
  : Check Environmental Balance of Folds
- [`summarise_fold_diagnostics()`](https://sodeidelphonse.github.io/isdmtools/reference/summarise_fold_diagnostics.md)
  [`print(`*`<FoldsSummary>`*`)`](https://sodeidelphonse.github.io/isdmtools/reference/summarise_fold_diagnostics.md)
  : Summarise Fold Diagnostics
- [`print(`*`<GeoDiagnostic>`*`)`](https://sodeidelphonse.github.io/isdmtools/reference/GeoDiagnostic-methods.md)
  [`plot(`*`<GeoDiagnostic>`*`)`](https://sodeidelphonse.github.io/isdmtools/reference/GeoDiagnostic-methods.md)
  : Methods for GeoDiagnostic objects
- [`print(`*`<EnvDiagnostic>`*`)`](https://sodeidelphonse.github.io/isdmtools/reference/EnvDiagnostic-methods.md)
  [`plot(`*`<EnvDiagnostic>`*`)`](https://sodeidelphonse.github.io/isdmtools/reference/EnvDiagnostic-methods.md)
  : Methods for EnvDiagnostic objects

## Model Evaluation

Background sampling and model performance evaluation.

- [`compute_metrics()`](https://sodeidelphonse.github.io/isdmtools/reference/compute_metrics.md)
  : Compute Evaluation Metrics for Integrated Spatial Models from
  Multisource Datasets

- [`sample_background()`](https://sodeidelphonse.github.io/isdmtools/reference/sample_background.md)
  : Generate background points

- [`plot(`*`<BackgroundPoints>`*`)`](https://sodeidelphonse.github.io/isdmtools/reference/BackgroundPoints-methods.md)
  [`print(`*`<BackgroundPoints>`*`)`](https://sodeidelphonse.github.io/isdmtools/reference/BackgroundPoints-methods.md)
  :

  Methods for `BackgroundPoints` objects

- [`print(`*`<ISDMmetrics>`*`)`](https://sodeidelphonse.github.io/isdmtools/reference/ISDMmetrics-methods.md)
  [`summary(`*`<ISDMmetrics>`*`)`](https://sodeidelphonse.github.io/isdmtools/reference/ISDMmetrics-methods.md)
  [`` `[`( ``*`<ISDMmetrics>`*`)`](https://sodeidelphonse.github.io/isdmtools/reference/ISDMmetrics-methods.md)
  [`plot(`*`<ISDMmetrics>`*`)`](https://sodeidelphonse.github.io/isdmtools/reference/ISDMmetrics-methods.md)
  [`as.data.frame(`*`<ISDMmetrics>`*`)`](https://sodeidelphonse.github.io/isdmtools/reference/ISDMmetrics-methods.md)
  [`get_background()`](https://sodeidelphonse.github.io/isdmtools/reference/ISDMmetrics-methods.md)
  : Methods for ISDMmetrics Objects

## Post-modelling analysis

Habitat suitability and niche overlap analysis.

- [`suitability_index()`](https://sodeidelphonse.github.io/isdmtools/reference/suitability_index.md)
  : Compute a unified suitability index from integrated spatial model
  predictions.
- [`prepare_predictions()`](https://sodeidelphonse.github.io/isdmtools/reference/prepare_predictions.md)
  : Obtain a formatted output from spatial predictions.
- [`calc_niche_overlap()`](https://sodeidelphonse.github.io/isdmtools/reference/calc_niche_overlap.md)
  : Calculate Niche Overlap (Schoener's D)

## Utilities & Helpers

Data cleaning and visualisation.

- [`fill_na_near()`](https://sodeidelphonse.github.io/isdmtools/reference/fill_na_near.md)
  : Fill NA cells with the nearest cells values
- [`generate_maps()`](https://sodeidelphonse.github.io/isdmtools/reference/generate_maps.md)
  : Generate Multi-panel Maps from Spatial Model Predictions
