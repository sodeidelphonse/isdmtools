# Changelog

## isdmtools 0.3.0

### Continuous Integration & Quality Assurance

- **GitHub Actions Integration**: Established a suite of automated
  workflows including R-CMD-check for cross-platform stability and
  Codecov for monitoring unit test coverage.
- **Automated Documentation**: Implemented an automated pkgdown
  deployment pipeline to ensure the package website and vignettes are
  updated upon every push to the main branch.

### Major Changes

- **New S3 Class Architecture**: Introduced `GeoDiagnostic`,
  `EnvDiagnostic`, and `FoldsSummary` classes to structure the spatial
  folds’ diagnostics results.
- **Diagnostic Framework**: Refined
  [`check_folds()`](https://sodeidelphonse.github.io/isdmtools/reference/check_folds.md)
  and
  [`check_env_balance()`](https://sodeidelphonse.github.io/isdmtools/reference/check_env_balance.md)
  constructors to evaluate spatial independence and environmental
  balance of spatial folds.
- **Documentation Refactor**: Implemented a grouped reference system to
  consolidate S3 methods, significantly improving the package Reference
  index.

### Improvements

- **Website Launch**: Official `pkgdown` site deployment with a custom
  navigation bar and categorised tutorials.
- **Tutorials**: Added a “Get Started” guide and an advanced “ISDM
  Evaluation Workflow” vignette with conditional evaluation for external
  dependencies.

## isdmtools 0.2.0

### New Features

- Add folds diagnostics tools: `check_folds` and `check_env_balance` are
  key methods operating on `DataFolds` objects to check the independence
  and representativeness of generated folds.

- New CV Methods: Added support for four new `cv_method` options for
  advanced spatial cross-validation:

  - “block”: Grid-based blocking via “spatialsample” (similar to
    “spatial” in blockCV but with exclusion buffer).
  - “nndm”: Nearest Neighbor Distance Matching for matching prediction
    and validation environments.
  - “buffer”: Distance-based exclusion zones (Leave-One-Out with
    buffer).
  - “location”: Leave-location-out/Leave-group-out CV, enabling
    spatiotemporal validation (e.g., by year) and source-specific
    validation (e.g., by observers, sites or regions).

- Expanded Cross-Validation Engines: Integrated ‘spatialsample’ as a
  backend, providing a unified interface for advanced spatial resampling
  methods.

- Buffer-Aware Extraction: Updated
  [`extract_fold()`](https://sodeidelphonse.github.io/isdmtools/reference/DataFolds-methods.md)
  to automatically handle exclusion zones. Points falling within spatial
  buffers are now correctly identified as NA and excluded from both
  training and testing sets to prevent spatial autocorrelation bias.

- Added `summary.DataFolds` method for providing clean aggregated
  statistics on data partition.

- Added relevant tests for `DataFolds` and `BackgroundPoints` classes
  and their constructors/methods.

### Enhanced Features

- Refactoring: Remove unnecessary dependencies (`reshape2`, `purrr` and
  `ggspatial`) and ensure consistent argument naming (snake_case).
- Unified CV Constructor:
  [`create_folds()`](https://sodeidelphonse.github.io/isdmtools/reference/create_folds.md)
  now acts as a high-level bridge between `blockCV` and `spatialsample`.
- Standardized S3 Outputs: The `DataFolds` object now standardizes
  internal indexing across different blocking engines, ensuring
  “folds_ids” are consistent regardless of the underlying package used.
- Improved Print Method: The `print.DataFolds` method now explicitly
  labels “Excluded” points, providing a clear summary of how many
  observations were buffered out of the validation process.
- Added robust unit tests for spatial buffer logic and multisource data
  fusion integrity.

## isdmtools 0.1.0.9000

### Refactoring & internal checks

- Integrated unit testing framework using `testthat`.
- Standardized coordinate column names to lowercase `x` and `y` across
  all functions.
- Updated internal logic to use `.data` pronouns to resolve global
  variable notes.

### Enhanced Features

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
