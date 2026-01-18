# isdmtools (development version)

* Added `ISDMmetrics` S3 class to handle outputs from `compute_metrics()`.
* Implemented S3 methods for `ISDMmetrics`: `print()`, `summary()`, `plot()`, `as.data.frame()`, and `[`.
* Enhanced `suitability_index()`: Added `"linear"` to `output.format` options.
* Updated `DESCRIPTION` with new minimum version requirements for package dependencies (`sf`, `terra`, `purrr`) to ensure compatibility with recent environment updates.


# isdmtools 0.1.0

* Initial release supporting "Integrating Presence-only and Abundance Data to Predict Baobab Distribution".
* DOI: 10.21203/rs.3.rs-7871875/v1
