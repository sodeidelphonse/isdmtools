# Compute Evaluation Metrics for Integrated Spatial Models from Multisource Datasets

This function computes a wide range of evaluation metrics for a
single-layer raster of model predictions against a list of one or more
point-based datasets. It is designed to handle different data types
(presence-only, presence-absence, and count data) and provides
individual metrics as well as dataset-weighted composite scores.

## Usage

``` r
compute_metrics(
  test_data,
  prob_raster = NULL,
  xy_excluded = NULL,
  expected_response = NULL,
  n_background = 1000,
  response_counts = "counts",
  response_pa = "present",
  threshold_method = c("best", "fixed"),
  best_method = c("youden", "closest.topleft"),
  fixed_threshold = NA_real_,
  best_threshold_policy = c("first", "last", "max.prec", "max.recall", "max.accu",
    "max.f1"),
  metrics = NULL,
  overall_roc_metrics = NULL,
  overall_error_metrics = NULL,
  is_pred_rate = FALSE,
  exposure = NULL,
  seed = 25,
  ...
)
```

## Arguments

- test_data:

  A named `list` of `sf` objects. Each `sf` object represents a
  different test dataset and must contain point geometries. The function
  will loop through each named dataset in the list. In particular,
  `test_data` can be a 'fold' from the
  [create_folds](https://sodeidelphonse.github.io/isdmtools/reference/create_folds.md)
  and
  [extract_fold](https://sodeidelphonse.github.io/isdmtools/reference/DataFolds-methods.md)
  outputs, if independent validation datasets are not available.

- prob_raster:

  A `SpatRaster` object with unique layer containing the model's
  predictions on a probability scale (0-1). It represents a suitability
  index, and its values are used to compute all ROC-based metrics (e.g.,
  AUC, TSS, F1 score). This argument is optional if only
  continuous-outcome metrics are requested for count data.

- xy_excluded:

  An optional `SpatVector` or `sf` object representing locations where
  pseudo-absence points should not be sampled, such as occupied areas or
  known background points. Only relevant for presence-only (PO) data.
  Default is `NULL`.

- expected_response:

  A `SpatRaster` object containing the model's predictions on a
  continuous scale (i.e. counts or rate if offset is used; see
  [suitability_index](https://sodeidelphonse.github.io/isdmtools/reference/suitability_index.md)).
  Its values are used to compute all continuous-outcome metrics (e.g.,
  RMSE, MAE, MAPE). This argument is required if a continuous-outcome
  metric is requested.

- n_background:

  integer. It specifies the number of pseudo-absence points to sample
  for presence-only data. Default is 1000 (see
  [sample_background](https://sodeidelphonse.github.io/isdmtools/reference/sample_background.md)).

- response_counts:

  character. The column name in the `sf` objects that contains observed
  counts. Default is 'counts' and must be standardized across all count
  data sets. Exceptionally, positive measurements (e.g. biomass) are
  supported by allowing exposure to its default value. In such cases,
  only continuous-outcome metrics can be requested.

- response_pa:

  character. The column name in the `sf` objects that contains
  presence-absence data (1 for presence, 0 for absence). Default is
  'present' and must be standardized across all PA data sets.

- threshold_method:

  character. The method to be used for selecting the threshold for
  converting probabilities to binary outcomes. Options are 'best' (using
  `best_method`) or 'fixed'. Default is "best".

- best_method:

  character. The method to be used for selecting the best threshold when
  `threshold_method` is 'best'. Options are 'youden' or
  'closest.topleft'. Default is "youden" criterion which maximizes the
  sensitivity and specificity.

- fixed_threshold:

  numeric. The value (0-1) to be used as the fixed threshold when
  `threshold_method` is 'fixed'. Default is `NA_real_`.

- best_threshold_policy:

  character. Specifies the policy for selecting a threshold when
  multiple thresholds yield the same 'best' value. Options are "first",
  "last", "max.prec" (max precision), "max.recall" (max recall),
  "max.accu" (max accuracy), or "max.f1" (max F1 score). Default is
  "first".

- metrics:

  character. A vector of the metric names to compute. If `NULL`, "auc"
  (area under the ROC curve), "tss" (true skill statistics), "accuracy",
  "F1" (F1 score), "precision", and "recall" are computed for ROC-based
  metrics while "rmse" (root mean squared error), "mae" (mean absolute
  error) and "r2" (pseudo R-squared) are computed for error-based
  metrics.

- overall_roc_metrics:

  character. A vector of a subset of ROC-based metrics to be used for
  the overall composite score (`TOT_ROC_SCORE`). Allowed options are
  "auc", "tss", "accuracy", and "F1". If `NULL`, the sensible default is
  "auc", "tss" and "accuracy". This metric is useful when the objective
  is to obtain a rapid overview of the rank of multiple candidate models
  fitted to datasets via blocked cross-validation using multi-criteria
  assessment.

- overall_error_metrics:

  character. A vector of a subset of continuous outcome metrics to be
  used for the overall composite score (`TOT_ERROR_SCORE`). Allowed
  options are "rmse", "mae", and "r2". If `NULL`, the default is "rmse"
  and "mae". In order to obtain an overall interpretable score, it is
  imperative to select metrics that have the same scale.

- is_pred_rate:

  logical. If `TRUE`, it indicates that the `expected_response` contains
  predictions at the intensity (per-unit-of-exposure) scale (typical for
  Bayesian models with offset from `inlabru`). If `FALSE`, it assumes
  predictions are at the original scale (e.g., counts). Default is
  `FALSE`.

- exposure:

  character. The column name in the `sf` objects that contains the
  exposure variable (offset). Only relevant for count (and sometimes
  presence-absence) data and must be standardized across all these types
  of datasets. If `is_pred_rate` is `TRUE`, observed counts are rescaled
  by this exposure variable. Default is `NULL`.

- seed:

  integer. It sets the seed for random number generation, used for
  pseudo-absence sampling to ensure reproducibility. Default is 25.

- ...:

  Additional arguments to be passed on to internal functions,
  particularly [coords](https://rdrr.io/pkg/pROC/man/coords.html)
  function.

## Value

An object of class `ISDMmetrics`. It is named `list` containing all
requested metrics. The names follow a consistent convention:

- `"<METRIC>_<DATASET_NAME>"`: Individual metric score for each dataset.

- `"<METRIC>_Comp"`: The sample-size-weighted composite score for a
  given metric across all valid datasets.

- `"TOT_ROC_SCORE"`: The overall ROC-based composite score, averaged
  across the selected `overall_roc_metrics`.

- `"TOT_ERROR_SCORE"`: The overall error-based composite score, averaged
  across `overall_error_metrics`.

The `ISDMmetrics` object is a list containing performance values for
individual datasets and composite scores. Use
[`as.data.frame()`](https://rspatial.github.io/terra/reference/as.data.frame.html)
to flatten these for cross-validation summaries.

## Details

The function handles three main data types and any combination thereof:

- **Presence-Absence (PA) Data:** The function uses the `response_pa`
  column and `prob_raster` to calculate all ROC-based metrics (see
  [coords](https://rdrr.io/pkg/pROC/man/coords.html), for more details
  on available metrics).

- **Count Data (or optionally measurements):** The function uses
  `expected_response` to calculate continuous-outcome metrics and can
  optionally use `prob_raster` to calculate ROC-based metrics for count
  data.

- **Presence-Only (PO) Data:** The function uses the presence points
  from the `sf` object (`xy_excluded`) and samples `n` pseudo-absence
  points from the study background (excluding `xy_excluded`) to create a
  presence-absence dataset for ROC-based metric calculations.

For models based on count data, if a user wants to compute both
continuous-outcome and ROC-based metrics, `expected_response` raster
must be supplied for the continuous metrics and `prob_raster` must also
be supplied for the ROC-based metrics. The `prob_raster` can be obtained
by converting the continuous-outcome prediction (e.g.,
`linear predictor`) to a suitability index using the
[suitability_index](https://sodeidelphonse.github.io/isdmtools/reference/suitability_index.md)
function.

The available continuous-outcome metrics are given as follows:

- **Root Mean Squared Error (RMSE)**: A measure of the average magnitude
  of the errors. It's the square root of the average of squared
  differences between prediction and actual observation. It gives higher
  weight to large errors. \$\$RMSE =
  \sqrt{\frac{1}{n}\sum\_{i=1}^{n}(\hat{y_i} - y_i)^2}\$\$.

- **Mean Absolute Error (MAE)**: A measure of the average magnitude of
  the errors without considering their direction. It is the average of
  the absolute differences between prediction and actual observation.
  \$\$MAE = \frac{1}{n}\sum\_{i=1}^{n}\|\hat{y_i} - y_i\|\$\$.

- **Mean Absolute Percentage Error (MAPE)**: A measure of prediction
  accuracy as a percentage. It is calculated as the average of the
  absolute percentage errors for each observation. It can be useful for
  comparing performance across different datasets or models. \$\$MAPE =
  \frac{100\\}{n}\sum\_{i=1}^{n}\|\frac{\hat{y_i} - y_i}{y_i}\|\$\$.

- **Pseudo R-squared (\\R^2\\)**: A measure of the proportion of
  variance in the observed data explained by the model's predictions.
  \$\$R^2 = 1 - \frac{SS\_{res}}{SS\_{tot}}\$\$ Where:

  - \\y_i\\ is the observed continuous value at location \\i\\.

  - \\\hat{y}\_i\\ is the predicted value from the model at location
    \\i\\ (e.g., the posterior mean of the predictions).

  - \\\bar{y}\\ is the mean of all observed values.

  - \\SS\_{res}\\ is the residual sum of squares, which measures the
    discrepancy between the observed and predicted values: \$\$SS\_{res}
    = \sum\_{i=1}^{n}(y_i - \hat{y}\_i)^2\$\$

  - \\SS\_{tot}\\ is the total sum of squares, which measures the total
    variance in the observed data: \$\$SS\_{tot} = \sum\_{i=1}^{n}(y_i -
    \bar{y})^2\$\$

A `weighted composite score` (`<METRIC>_Comp`) is computed for each
requested metric by taking the sample-size-weighted average across all
datasets where the metric was successfully calculated. A
`total composite score` (`TOT_ROC_SCORE` or `TOT_ERROR_SCORE`) is also
computed by averaging the selected metrics in the corresponding 'overall
metrics' character vector. It can be viewed as a quick *multi-criterion
decision metric* for multiple models comparison.

## See also

[`extract_fold`](https://sodeidelphonse.github.io/isdmtools/reference/DataFolds-methods.md),
[`sample_background`](https://sodeidelphonse.github.io/isdmtools/reference/sample_background.md)

Other ISDM evaluation methods:
[`ISDMmetrics-methods`](https://sodeidelphonse.github.io/isdmtools/reference/ISDMmetrics-methods.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Assuming you have dummy prediction rasters and a list of sf objects
# with 'counts' and 'present' columns for counts and PA data, respectively.

# Example 1: Compute metrics for a presence-absence model
# pa_metrics <- compute_metrics(
#   test_data = list(ds1 = my_pa_sf),
#   prob_raster = prob_raster,  # compulsory prob_raster
#   response_pa = "present"     # default labels column for all PA data
# )

# Example 2: Compute continuous-outcome metrics for a count-based model
# cont_metrics <- compute_metrics(
#   test_data = list(ds1 = my_count_sf),
#   expected_response = expected_raster, # prediction on count scale
#   response_count = "counts",           # default labels column for all counts
#   metrics = c("rmse", "mae", "mape")
# )

# Example 3: Compute both continuous and ROC-based metrics for a count model
# The user must first generate a suitability index (prob_raster & expected_response)
# from the linear scale prediction (pred_eta).

# expected_raster <- suitability_index(pred_eta,
#                    response_type = "count",
#                    output_format = "response")
# suitability_raster <- suitability_index(pred_eta,
#                       response_type = "count",
#                       output_format = "prob")
# full_metrics <- compute_metrics(
#   test_data = list(ds1 = my_count_sf),
#   prob_raster = suitability_raster,
#   expected_response = expected_raster,
#   metrics = c("rmse", "mae", "auc", "tss")
# )

# Example 4: Handle an inlabru-like model with an offset term
# The `expected_response` raster is at the intensity scale (rate).

# cont_metrics <- compute_metrics(
#   test_data = list(ds1 = my_count_sf),
#   prob_raster = suitability_raster,
#   expected_response = expected_raster,
#   metrics = c("rmse", "auc", "tss"),
#   is_pred_rate = TRUE,
#   exposure = "exposure_col" # Exposure column (e.g. sampling unit area)
# )

# Example 5: Compute dataset-specific and weighted composite metrics for a joint model
# expected_raster  <- suitability_index(pred_eta,
#                     response_type = "count.pa",
#                     output_format = "response")
# suitability_raster <- suitability_index(pred_eta,
#                       response_type = "count.pa",
#                       has_offset = FALSE)
# full_metrics <- compute_metrics(
#   test_data = list(ds1 = my_count_sf, ds2 = my_pa_sf),
#   prob_raster = suitability_raster,
#   expected_response = expected_raster,
#   metrics = c("rmse", "mae", "auc", "tss", "accuracy")
# )
} # }
```
