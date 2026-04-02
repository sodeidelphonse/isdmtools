# Changelog

## isdmtools 0.4.0

### New Features:

- [`compute_ppc_stats()`](https://sodeidelphonse.github.io/isdmtools/reference/compute_ppc_stats.md):
  A wrapper to calculate Pearson Chi-squared statistics and Bayesian
  $p$-values to assess model fit.
- Added S3 methods: `print.ppc_stats` and `plot.ppc_stats` for immediate
  visual and console-based model diagnostics.
- [`simulate_replicates()`](https://sodeidelphonse.github.io/isdmtools/reference/simulate_replicates.md):
  A new vectorized simulation engine to generate $y_{rep}$ data from
  posterior samples. Supports `poisson`, `nbinomial`, and
  `quasi-poisson` families.
- Posterior Predictive Checks (PPC): Introduced a robust Bayesian
  diagnostic suite for count data models.
- `plot.std_matern_corr`: An S3 method for quick plot of the
  standardized outputs and pairwise correlation function.
- `std_matern_corr`: A wrapper to standardize Matern covariance
  parameters across estimation engines and compute the pairwise
  correlation from LGCP models.

### Improvements

- [`format_predictions()`](https://sodeidelphonse.github.io/isdmtools/reference/format_predictions.md):
  Updated to ensure extracted coordinates are consistently lowercase for
  internal processing while preserving existing X/Y data columns.
- Moved
  [`solve_practical_range()`](https://sodeidelphonse.github.io/isdmtools/reference/solve_practical_range.md)
  from `utils.R` to `matern_helpers.R` and optimized `uniroot` calls by
  removing the redundant `sigma2` argument.
- Added internal utility functions `matern_cov()` and `get_variance()`
  to standardize spatial covariance across engines and variance
  calculations for Poisson extensions.

## isdmtools 0.3.0

### Continuous Integration & Deployment

- Code styling: Established the `lintr` workflow for the automated
  detection of syntax errors and code styling issues.
- GitHub Actions Integration: Established a suite of automated workflows
  including R-CMD-check for cross-platform stability and `Codecov` for
  monitoring unit tests coverage.
- Automated Documentation: Implemented an automated pkgdown deployment
  pipeline to ensure the package website and vignettes are updated upon
  every push to the main branch.
- Website Launch: Deployed an official `pkgdown` site with a custom
  navigation bar and categorized tutorials.

### New Features

- Tutorials: Added a “Get Started” guide and an advanced “ISDM
  Evaluation Workflow” vignette with conditional evaluation for external
  dependencies.
- New S3 Class Architecture: Introduced `GeoDiagnostic`,
  `EnvDiagnostic`, and `FoldsSummary` classes to structure the spatial
  folds’ diagnostics results.
- Documentation Refactor: Implemented a grouped reference system to
  consolidate S3 methods for improving the package Reference index.

### Improvements

- R Compatibility: Full support for the native pipe `|>` across the
  package.
- Codecov coverage: Set up additional units tests to boost the Codecov
  test coverage.
- BugFixes: Refactoring `ISDMmetrics` methods to handle efficiently
  datasets’ names with underscores to avoid names collision during
  string splitting.
- Blocks Diagnostics: Refined
  [`check_folds()`](https://sodeidelphonse.github.io/isdmtools/reference/check_folds.md)
  and
  [`check_env_balance()`](https://sodeidelphonse.github.io/isdmtools/reference/check_env_balance.md)
  constructors for evaluating spatial independence and environmental
  balance of spatial folds.

## isdmtools 0.2.0

### New Features

- Added folds diagnostics tools: `check_folds` and `check_env_balance`
  are key methods operating on `DataFolds` objects to check the
  independence and representativeness of folds.

- Buffer-Aware Extraction: Updated
  [`extract_fold()`](https://sodeidelphonse.github.io/isdmtools/reference/DataFolds-methods.md)
  to automatically handle exclusion zones. Points falling within spatial
  buffers are now correctly identified as NA and excluded from both
  training and testing sets to prevent spatial autocorrelation bias.

- New CV Methods: Added support for four new `cv_method` options for
  advanced spatial cross-validation:

  - `"block"`: Grid-based blocking via “spatialsample” (similar to
    “spatial” in blockCV but with exclusion buffer).
  - `"nndm"`: Nearest Neighbor Distance Matching for matching prediction
    and validation environments.
  - `"buffer"`: Distance-based exclusion zones (Leave-One-Out with
    buffer).
  - `"location"`: Leave-location-out/Leave-group-out CV, enabling
    spatiotemporal validation (e.g., by year) and source-specific
    validation (e.g., by observers, sites or regions).

- Expanded Cross-Validation Engines: Integrated ‘spatialsample’ as a
  backend, providing a unified interface for advanced spatial resampling
  methods.

- Added `summary.DataFolds` method for providing clean aggregated
  statistics on data partitioning.

- Added relevant tests for `DataFolds` and `BackgroundPoints` classes
  and their constructors/methods.

### Improvements

- Refactoring: Remove unnecessary dependencies (`reshape2`, `purrr` and
  `ggspatial`) and ensure consistent argument naming (snake_case).
- Standardized S3 Outputs: The `DataFolds` object now standardizes
  internal indexing across different blocking engines, ensuring
  “folds_ids” are consistent regardless of the underlying package used.
- Improved Print Method: The `print.DataFolds` method now explicitly
  labels “Excluded” points, providing a clear summary of how many
  observations were buffered out of the validation process.
- Added robust unit tests for spatial buffer logic and multisource data
  fusion integrity.
- Unified CV Constructor:
  [`create_folds()`](https://sodeidelphonse.github.io/isdmtools/reference/create_folds.md)
  now acts as a high-level bridge between `blockCV` and `spatialsample`.

## isdmtools 0.1.0.9000

### Refactoring & Testing Framework

- Integrated unit testing framework using `testthat`.
- Standardized coordinates column names to lowercase `x` and `y` across
  all functions.
- Updated internal logic to use `.data` pronoun to resolve global
  variable notes in R-CMD-check.

### Improvements

- Implemented S3 methods for `ISDMmetrics`:
  [`print()`](https://rdrr.io/r/base/print.html),
  [`summary()`](https://rspatial.github.io/terra/reference/summary.html),
  [`plot()`](https://rspatial.github.io/terra/reference/plot.html),
  [`as.data.frame()`](https://rspatial.github.io/terra/reference/as.data.frame.html),
  and `[`.
- Added `ISDMmetrics` S3 class to handle outputs from
  [`compute_metrics()`](https://sodeidelphonse.github.io/isdmtools/reference/compute_metrics.md).
- Implemented `get_background` helper function for selected background
  points.
- Enhanced
  [`suitability_index()`](https://sodeidelphonse.github.io/isdmtools/reference/suitability_index.md):
  Added `"linear"` to `output.format` options.
- Updated `DESCRIPTION` with new minimum version requirements for
  package dependencies (`sf`, `terra`, `purrr`) to ensure compatibility
  with recent environment updates.

## isdmtools 0.1.0

- Initial release supporting “Integrating Presence-only and Abundance
  Data to Predict Baobab (*Adansonia digitata L.*) Distribution: A
  Bayesian Data Fusion Framework”. DOI: 10.21203/rs.3.rs-7871875/v1
