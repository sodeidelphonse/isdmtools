# Methods for ISDMmetrics Objects

Objects of class `ISDMmetrics` are returned by
[`compute_metrics`](https://sodeidelphonse.github.io/isdmtools/reference/compute_metrics.md).
These methods provide structured ways to view, summarize and manipulate
the evaluation results.

`get_background` is a helper function to extract the `BackgroundPoints`
object from an `ISDMmetrics` object.

## Usage

``` r
# S3 method for class 'ISDMmetrics'
print(x, ...)

# S3 method for class 'ISDMmetrics'
summary(object, ...)

# S3 method for class 'ISDMmetrics'
x[...]

# S3 method for class 'ISDMmetrics'
plot(x, include_composite = TRUE, ...)

# S3 method for class 'ISDMmetrics'
as.data.frame(x, ...)

get_background(x)
```

## Arguments

- x:

  An object of class `ISDMmetrics`.

- ...:

  Additional arguments passed on to the method.

- object:

  An object of class `ISDMmetrics`.

- include_composite:

  Logical. Should the weighted composite scores be included in the plot?
  Defaults to `TRUE`.

## Value

- `print`: Invisibly returns the original object.

- `summary`: Invisibly returns `NULL`.

- `plot`: Returns a `ggplot2` object.

- `[`: Returns a subset of `ISDMmetrics` object.

- `as.data.frame`: Returns a tidy `data.frame` in long format.

The `BackgroundPoints` object if present, otherwise `NULL`.

## See also

Other ISDM evaluation methods:
[`compute_metrics()`](https://sodeidelphonse.github.io/isdmtools/reference/compute_metrics.md)

## Examples

``` r
if (FALSE) { # \dontrun{
#--- Compute metrics for an Integrated SDM
# This object will contain metrics for e.g., Presence-only and Count data
eval_results <- compute_metrics(
  test_data = test_data,
  prob_raster = suitability_raster,
  expected_response = expected_raster,
  n_background = 1000,
  metrics = c("rmse", "mae", "auc", "tss"),
  is_pred_rate = TRUE, # model with offset
  exposure = "area"    # standardized exposure name across the counts data
)

#--- Quick view of the results
print(eval_results)

#--- Generate a full replication report
# summary.ISDMmetrics shows seeds, threshold logic, and prediction type
summary(eval_results)

#--- Visual comparison of metrics
plot(eval_results, include_composite = TRUE)

#--- Background visualization (for Presence-Only data)
bg_data <- get_background(eval_results)
if (!is.null(bg_data)) {
  plot(bg_data)
}

#--- Export to tabular format for external reports
results_df <- as.data.frame(eval_results)
head(results_df)

#--- Subset specific metrics for custom analysis
auc_only <- eval_results[grep("AUC", names(eval_results))]
} # }
```
