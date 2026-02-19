# isdmtools 0.3.0

## Continuous Integration/Deployment & Quality Assurance
- GitHub Actions Integration: Established a suite of automated workflows including R-CMD-check for cross-platform stability and Codecov for monitoring unit test coverage.
- Automated Documentation: Implemented an automated pkgdown deployment pipeline to ensure the package website and vignettes are updated upon every push to the main branch.

## New Features
- Tutorials: Added a "Get Started" guide and an advanced "ISDM Evaluation Workflow" vignette with conditional evaluation for external dependencies.
- New S3 Class Architecture: Introduced `GeoDiagnostic`, `EnvDiagnostic`, and `FoldsSummary` classes to structure the spatial folds' diagnostics results.
- Documentation Refactor: Implemented a grouped reference system to consolidate S3 methods, significantly improving the package Reference index.

## Improvements
- Codecov coverage: Set up additional units tests to boost their coverage.
- Website Launch: Official `pkgdown` site deployment with a custom navigation bar and categorised tutorials.
- BugFixes: Refactoring `ISDMmetrics` methods to handle efficiently datasets' names with underscores to avoid names collision during string splitting.
- Blocks Diagnostics: Refined `check_folds()` and `check_env_balance()` constructors for evaluating spatial independence and environmental balance of spatial folds.

# isdmtools 0.2.0

## New Features
* Added folds diagnostics tools: `check_folds` and `check_env_balance` are key methods operating on `DataFolds` objects to check the independence and representativeness of generated folds.
* Buffer-Aware Extraction: Updated `extract_fold()` to automatically handle exclusion zones. 
Points falling within spatial buffers are now correctly identified as NA and excluded from both training and testing sets to prevent spatial autocorrelation bias.

* New CV Methods: Added support for four new `cv_method` options for advanced spatial cross-validation:
   * `"block"`: Grid-based blocking via "spatialsample" (similar to "spatial" in blockCV but with exclusion buffer).
   * `"nndm"`: Nearest Neighbor Distance Matching for matching prediction and validation environments.
   * `"buffer"`: Distance-based exclusion zones (Leave-One-Out with buffer).
   * `"location"`: Leave-location-out/Leave-group-out CV, enabling spatiotemporal validation (e.g., by year) and source-specific validation (e.g., by observers, sites or regions).
   
* Expanded Cross-Validation Engines: Integrated 'spatialsample' as a backend, providing a unified interface for advanced spatial resampling methods.
* Added `summary.DataFolds` method for providing clean aggregated statistics on data partition.
* Added relevant tests for `DataFolds` and `BackgroundPoints` classes and their constructors/methods. 

## Enhanced Features
* Refactoring: Remove unnecessary dependencies (`reshape2`, `purrr` and `ggspatial`) and ensure consistent argument naming (snake_case).
* Unified CV Constructor: `create_folds()` now acts as a high-level bridge between `blockCV` and `spatialsample`.
* Standardized S3 Outputs: The `DataFolds` object now standardizes internal indexing across different blocking engines, ensuring "folds_ids" are consistent regardless of the underlying package used.
* Improved Print Method: The `print.DataFolds` method now explicitly labels "Excluded" points, providing a clear summary of how many observations were buffered out of the validation process.
* Added robust unit tests for spatial buffer logic and multisource data fusion integrity.

# isdmtools 0.1.0.9000

## Refactoring & Testing Framework
* Integrated unit testing framework using `testthat`.
* Standardized coordinate column names to lowercase `x` and `y` across all functions.
* Updated internal logic to use `.data` pronouns to resolve global variable notes.

## Enhanced Features
* Implemented S3 methods for `ISDMmetrics`: `print()`, `summary()`, `plot()`, `as.data.frame()`, and `[`.
* Added `ISDMmetrics` S3 class to handle outputs from `compute_metrics()`.
* Implemented `get_background` helper function for selected background points.
* Enhanced `suitability_index()`: Added `"linear"` to `output.format` options.
* Updated `DESCRIPTION` with new minimum version requirements for package dependencies (`sf`, `terra`, `purrr`) to ensure compatibility with recent environment updates.

# isdmtools 0.1.0

* Initial release supporting "Integrating Presence-only and Abundance Data to Predict Baobab (*Adansonia digitata L.*) Distribution: A Bayesian Data Fusion Framework". DOI: 10.21203/rs.3.rs-7871875/v1
