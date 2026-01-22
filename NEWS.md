# isdmtools 0.1.0.9000

## New Features
* Added `spatialsample` blocking methods to `create_folds()` to support advanced spatial cross-validation.
* Integrated unit testing framework using `testthat`.

## Refactoring
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
