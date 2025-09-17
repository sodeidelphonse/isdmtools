
#--------------------------------------------------------------------------------------------------------
#--- Function to compute comprehensive evaluation metrics for integrated or standalone spatial models
#--------------------------------------------------------------------------------------------------------

#' @title Compute Evaluation Metrics for Integrated Spatial Models from Multisource Datasets
#'
#' @description This function computes a wide range of evaluation metrics for a single-layer raster of model predictions against a list of one or more point-based datasets. It is designed to handle different data types (presence-only, presence-absence, and count data) and provides individual metrics as well as dataset-weighted composite scores.
#'
#' @param test.data A named `list` of `sf` objects. Each `sf` object represents a different test dataset and must contain point geometries. The function will loop through each named dataset in the list. In particular, `test.data` can be a 'fold' of the 'test' set from `create_folds()` and `extract_fold()` outputs, if independent validation datasets are not available.
#' @param prob.raster A `SpatRaster` object with unique layer containing the model's predictions on a probability scale (0-1). It represents a suitability index, and its values are used to compute all ROC-based metrics (e.g., AUC, TSS, F1 score). This argument is optional if only continuous-outcome metrics are requested for count data.
#' @param expected.response A `SpatRaster` object containing the model's predictions on a continuous scale (i.e. counts or rate if offset used; see `suitability_index()`). Its values are used to compute all continuous-outcome metrics (e.g., RMSE, MAE, MAPE). This argument is required if a continuous-outcome metric is requested.
#' @param xy.excluded An optional `SpatVector` or `sf` object representing locations where pseudo-absence points should not be sampled, such as occupied areas or known background points. Only relevant for presence-only (PO) data. Default is `NULL`.
#' @param n.background An integer specifying the number of pseudo-absence points to sample for presence-only data. Default is 1000 (see \link{sample_bg_points}).
#' @param responseCounts A character string representing the column name in the `sf` objects that contains observed counts. Default is 'counts' and must be standardized across all count data sets.
#' @param responsePA A character string representing the column name in the `sf` objects that contains presence-absence data (1 for presence, 0 for absence). Default is 'present' and must be standardized across all PA data sets.
#' @param seed An integer for setting the seed for random number generation, used for pseudo-absence sampling to ensure reproducibility. Default is 25.
#' @param threshold.method A character string specifying how to select the threshold for converting probabilities to binary outcomes. Options are 'best' (using `best.method`) or 'fixed'. Default is "best".
#' @param best.method A character string specifying the method for selecting the best threshold when `threshold.method` is 'best'. Options are 'youden' or 'closest.topleft'. Default is "youden" criterion which maximizes the sensitivity and specificity.
#' @param fixed.threshold A numeric value (0-1) to use as the fixed threshold when `threshold.method` is 'fixed'. Default is `NA_real_`.
#' @param best.threshold.policy A character string specifying the policy for selecting a threshold when multiple thresholds yield the same 'best' value. Options are "first", "last", "max.prec" (max precision), "max.recall" (max recall), "max.accu" (max accuracy), or "max.f1" (max F1 score). Default is "first".
#' @param metrics A character vector of metric names to compute. If `NULL`, "auc" (area under the ROC curve), "tss" (true skill statistics), "accuracy", "F1" (F1 score), "precision", and "recall" are computed for ROC-based metrics while "rmse" and "mae" are computed for error-based metrics.
#' @param overall.roc.metrics A character vector specifying a subset of ROC-based metrics to be used for the overall composite score (`TOT_ROC_SCORE`). Allowed options are "auc", "tss", "accuracy", and "F1". If `NULL`, the sensible default is "auc", "tss" and "accuracy".
#' @param overall.error.metrics A character vector specifying a subset of continuous outcome metrics to be used for the overall composite score (`TOT_ERROR_SCORE`).
#' Allowed options are "rmse" (root mean squared error), "mae" (mean absolute error), "mape" (mean absolute percentage error) and "r2" (pseudo R-squared). If `NULL`, the default is "rmse" and "mae".
#' @param is.pred.rate A logical value. If `TRUE`, it indicates that the `expected.response` contains predictions at the intensity (per-unit-of-exposure) scale (typical for Bayesian models with offset from `inlabru`). If `FALSE`, it assumes predictions are at the original scale (e.g., counts). Default is `FALSE`.
#' @param exposure A character string representing the column name in the `sf` objects that contains the exposure variable (offset). Only relevant for count (and sometimes presence-absence) data and must be standardized across all these types of datasets.
#' If `is.pred.rate` is `TRUE`, observed counts are rescaled by this exposure variable. Default is `NULL`.
#' @param ... Additional arguments to be passed to internal functions, particularly \link[pROC]{coords} function.
#'
#' @details The function handles three main data types and any combination thereof:
#'
#' \itemize{
#'   \item \strong{Presence-Absence (PA) Data:} The function uses the `responsePA` column and `prob.raster` to calculate all ROC-based metrics (see \link[pROC]{coords}, for more details on available metrics).
#'   \item \strong{Count Data:} The function uses `expected.response` to calculate continuous-outcome metrics and can optionally use `prob.raster` to calculate ROC-based metrics.
#'   \item \strong{Presence-Only (PO) Data:} The function uses the presence points from the `sf` object (`xy.excluded`) and samples `n` pseudo-absence points from the study background (excluding `xy.excluded`) to create a presence-absence dataset for ROC-based metric calculations.
#' }
#'
#' For models based on count data, if a user wants to compute both continuous-outcome and ROC-based metrics, `expected.response` raster must be supplied for the continuous metrics and `prob.raster` must also be supplied for the ROC-based metrics.
#' The `prob.raster` can be obtained by converting the continuous-outcome prediction (e.g., `linear predictor`) to a suitability index using the \link{suitability_index} function.
#'
#' The available continuous-outcome metrics are given as follows:
#' \itemize{
#' \item **Root Mean Squared Error (RMSE)**: A measure of the average magnitude of the errors. It's the square root of the average of squared differences between prediction and actual observation. It gives higher weight to large errors.
#'    \deqn{RMSE = \sqrt{\frac{1}{n}\sum_{i=1}^{n}(\hat{y_i} - y_i)^2}}.
#' \item **Mean Absolute Error (MAE)**: A measure of the average magnitude of the errors without considering their direction. It is the average of the absolute differences between prediction and actual observation.
#'    \deqn{MAE = \frac{1}{n}\sum_{i=1}^{n}|\hat{y_i} - y_i|}.
#' \item **Mean Absolute Percentage Error (MAPE)**: A measure of prediction accuracy as a percentage. It is calculated as the average of the absolute percentage errors for each observation. It can be useful for comparing performance across different datasets or models.
#'    \deqn{MAPE = \frac{100\%}{n}\sum_{i=1}^{n}|\frac{\hat{y_i} - y_i}{y_i}|}.
#' \item **Pseudo R-squared (\eqn{R^2})**: A measure of the proportion of variance in the observed data explained by the model's predictions.
#'    \deqn{R^2 = 1 - \frac{SS_{res}}{SS_{tot}}}
#'    Where:
#'   \itemize{
#'   \item \eqn{y_i} is the observed continuous value at location \eqn{i}.
#'   \item \eqn{\hat{y}_i} is the predicted value from the model at location \eqn{i} (e.g., the posterior mean of the predictions).
#'   \item \eqn{\bar{y}} is the mean of all observed values.
#'   \item \eqn{SS_{res}} is the residual sum of squares, which measures the discrepancy between the observed and predicted values:
#'     \deqn{SS_{res} = \sum_{i=1}^{n}(y_i - \hat{y}_i)^2}
#'   \item \eqn{SS_{tot}} is the total sum of squares, which measures the total variance in the observed data:
#'     \deqn{SS_{tot} = \sum_{i=1}^{n}(y_i - \bar{y})^2}
#'   }
#' }
#'
#' A `weighted composite score` (`<METRIC>_Comp`) is computed for each requested metric by taking the sample-size-weighted average across all datasets where the metric was successfully calculated.
#' A `total composite score` (`TOT_ROC_SCORE` or `TOT_ERROR_SCORE`) is also computed by averaging the selected metrics in the corresponding `overall metrics` character vector.
#' It can be viewed as a quick *multi-criterion decision metric* for several models comparison.
#'
#' @return A named `list` containing all requested metrics. The names follow a consistent convention:
#' \itemize{
#'   \item `"<METRIC>_<DATASET_NAME>"`: Individual metric score for each dataset.
#'   \item `"<METRIC>_Comp"`: The sample-size-weighted composite score for a given metric across all valid datasets.
#'   \item `"TOT_ROC_SCORE"`: The overall ROC-based composite score, averaged across the selected `overall.roc.metrics`.
#'   \item `"TOT_ERROR_SCORE"`: The overall error-based composite score, averaged across `overall.error.metrics`.
#' }
#'
#' @examples
#' \dontrun{
#' # Assuming you have dummy rasters and a list of sf objects
#' # with 'counts' and 'present' columns.
#'
#' # Example 1: Compute metrics for a presence-absence model
#' # my_metrics <- compute_metrics(
#' #   test.data = list(ds1 = my_pa_sf),
#' #   prob.raster = my_prob_raster  # compulsory prob.raster
#' # )
#'
#' # Example 2: Compute continuous-outcome metrics for a count-based model
#' # cont_metrics <- compute_metrics(
#' #   test.data = list(ds1 = my_count_sf),
#' #   expected.response = expected_raster, # prediction raster on count scale
#' #   metrics = c("rmse", "mae", "mape")
#' # )
#'
#' # Example 3: Compute both continuous and ROC-based metrics for a count model
#' # The user must first generate a suitability index (prob_raster)
#' # from the linear scale prediction (pred_eta).
#' # expected_raster <- suitability_index(pred_eta, response.type = "count",
#' #                    output.format = "response")
#' # suitability_raster <- suitability_index(pred_eta, response.type = "count",
#' #                        output.format = "prob")
#' # full_metrics <- compute_metrics(
#' #   test.data = list(ds1 = my_count_sf),
#' #   prob.raster = suitability_raster,
#' #   expected.response = expected_raster,
#' #   metrics = c("rmse", "mae", "auc", "tss")
#' # )
#'
#' # Example 4: Handle an inlabru-like model with an offset term
#' # The `expected.response` raster is at the intensity scale (rate).
#' # cont_metrics <- compute_metrics(
#' #   test.data = list(ds1 = my_count_sf),
#' #   prob.raster = suitability_raster,
#' #   expected.response = expected_raster,
#' #   metrics = c("rmse", "auc", "tss"),
#' #   is.pred.rate = TRUE,
#' #   exposure = "exposure_col"
#' # )
#'
#' # Example 5: Compute dataset-specific and weighted composite metrics for a joint model
#' # expected_raster  <- suitability_index(pred_eta, response.type = "count.pa",
#' #                  output.format = "response")
#' # suitability_raster <- suitability_index(pred_eta,
#' #                    response.type = "count.pa", has.offset = FALSE)
#' # full_metrics <- compute_metrics(
#' #   test.data = list(ds1 = my_count_sf, ds2 = my_pa_sf),
#' #   prob.raster = suitability_raster,
#' #   expected.response = expected_raster,
#' #   metrics = c("rmse", "mae", "auc", "tss", "accuracy")
#' # )
#' }
#' @export
#' @seealso \code{\link{extract_fold}}, \code{\link{suitability_index}}
#'
compute_metrics <- function(test.data,
                            prob.raster = NULL,
                            xy.excluded = NULL,
                            expected.response = NULL,
                            n.background = 1000,
                            responseCounts = 'counts',
                            responsePA = 'present',
                            threshold.method = c("best", "fixed"),
                            best.method = c("youden", "closest.topleft"),
                            fixed.threshold = NA_real_,
                            best.threshold.policy = c("first", "last", "max.prec", "max.recall", "max.accu", "max.f1"),
                            metrics = NULL,
                            overall.roc.metrics = NULL,
                            overall.error.metrics = NULL,
                            is.pred.rate = FALSE,
                            exposure = NULL,
                            seed = 25, ...) {

  # Master list of allowed metrics
  error_metrics <- c("rmse", "mae", "mape", "r2")
  roc_metrics <- c("auc", "tss", "accuracy", "precision", "recall", "specificity",
                   "npv", "fpr", "fnr", "fdr", "F1", "threshold", "tpr", "tnr")
  master_allowed_metrics <- c(roc_metrics, error_metrics)

  # --- Critical validation for required inputs ----
  threshold.method <- match.arg(threshold.method)
  best.method      <- match.arg(best.method)
  best.threshold.policy <- match.arg(best.threshold.policy)

  has_roc_metrics <- any(roc_metrics %in% metrics)
  if (has_roc_metrics && is.null(prob.raster)) {
    stop("ROC-based metrics were requested but 'prob.raster' (suitability index) was not provided.", call. = FALSE)
  }

  has_error_metrics <- any(error_metrics %in% metrics)
  if (has_error_metrics && is.null(expected.response)) {
    stop("Continuous-outcome metrics were requested but 'expected.response' response raster was not provided.", call. = FALSE)
  }

  default_roc_metrics  <- c("auc", "tss", "accuracy", "F1", "recall", "precision")
  default_error_metrics <- c("rmse", "mae")

  if (is.null(metrics)) {
    if (!is.null(prob.raster) && !is.null(expected.response)) {
      metrics <- c(default_roc_metrics, default_error_metrics)
    } else if (!is.null(prob.raster)) {
      metrics <- default_roc_metrics
    } else if (!is.null(expected.response)) {
      metrics <- default_error_metrics
    }
  }

  # --- Validation of the main data ----
  if (!is.logical(is.pred.rate) || length(is.pred.rate) != 1) {
    stop("'is.pred.rate' must be a single logical value (TRUE or FALSE).", call. = FALSE)
  }

  if (threshold.method == "fixed") {
    if (!is.numeric(fixed.threshold) || length(fixed.threshold) != 1 || !is.finite(fixed.threshold) || fixed.threshold < 0 || fixed.threshold > 1) {
      stop("'fixed.threshold' must be a single finite numeric value between 0 and 1 when 'threshold.method' is 'fixed'.", call. = FALSE)
    }
  }

  if(!is.null(prob.raster) || !missing(prob.raster)) {
    if(!inherits(prob.raster, "SpatRaster")) {
      stop(sprintf("'%s' must be a SpatRaster object.", prob.raster), call. = FALSE)
    } else {
      vals <- terra::values(prob.raster)[,1]
      vals <- vals[!is.na(vals)]
      if (length(vals) == 0) {
        stop("Probabilities raster contains only NA values. Cannot compute metrics.", call. = FALSE)
      }

      if (any(vals < 0) || any(vals > 1)) {
        stop("Probabilities raster must be positive with range [0, 1]", call. = FALSE)
      }

      if (terra::nlyr(prob.raster) != 1) {
        stop(sprintf("The input '%s' must be a single-layer raster.", prob.raster), call. = FALSE)
      }
    }
  }

  if(!is.null(expected.response)) {
    if(!inherits(expected.response, "SpatRaster")) {
      stop(sprintf("'%s' must be a SpatRaster object.", deparse(substitute(expected.response))), call. = FALSE)
    }
    vals <- terra::values(expected.response)[,1]
    vals <- vals[!is.na(vals)]
    if(any(vals < 0)) {
      stop("'expected.response' must be positive values.", call. = FALSE)
    }
  }

  if(!inherits(test.data, "list") || is.null(names(test.data))) {
    stop(sprintf("'%s' must be a named list.", deparse(substitute(test.data))), call. = FALSE)
  }

  for (ds_name in names(test.data)) {
    if (!inherits(test.data[[ds_name]], "sf")) {
      stop(sprintf("All elements in 'test.data' must be 'sf' objects. '%s' is not.", ds_name), call. = FALSE)
    }
    if (!inherits(sf::st_geometry(test.data[[ds_name]]), "sfc_POINT")) {
      warning(sprintf("Dataset '%s' in 'test.data' does not contain point geometries and could cause issue to values extraction.", ds_name), call. = FALSE)
    }
  }

  # --- Validation of arguments for the overall composite scores ----
  if (!all(metrics %in% master_allowed_metrics)) {
    stop(sprintf("Invalid metric(s) requested. Allowed metrics are: %s", paste(master_allowed_metrics, collapse = ", ")), call. = FALSE)
  }

  if (is.null(overall.roc.metrics)) {
    roc_metrics_to_average <- c("auc", "tss", "accuracy")
    message(sprintf("No 'overall.roc.metrics' specified. Overall ROC score will be computed for available metrics among: %s.",
                    paste(toupper(roc_metrics_to_average), collapse = ", ")))
  } else {
    if (!is.character(overall.roc.metrics)) {
      stop("'overall.roc.metrics' must be a character vector of metric names.", call. = FALSE)
    }
    valid_options <- c("auc", "tss", "accuracy", "F1")
    if(!all(overall.roc.metrics %in% valid_options)) {
      stop(sprintf("Invalid metric(s) selected for 'overall.roc.metrics': %s. Allowed composite metrics are: %s.",
                   paste(setdiff(overall.roc.metrics, valid_options), collapse = ", "),
                   paste(valid_options, collapse = ", ")), call. = FALSE)
    }
    if (!all(overall.roc.metrics %in% metrics)) {
      stop(sprintf("Metric(s) selected for 'overall.roc.metrics' (%s) were not requested in the general 'metrics' argument. Please ensure all composite metrics are also in 'metrics'.",
                   paste(setdiff(overall.roc.metrics, metrics), collapse = ", ")), call. = FALSE)
    }
    if (length(overall.roc.metrics) < 1) {
      stop("At least one metric must be requested for 'overall.roc.metrics' to calculate an overall composite score.", call. = FALSE)
    }
    roc_metrics_to_average <- overall.roc.metrics
  }

  if (is.null(overall.error.metrics)) {
    error_metrics_to_average <- c("rmse", "mae")
    message(sprintf("No 'overall.error.metrics' specified. Overall continuous-outcome score will be computed for available metrics among: %s.",
                    paste(toupper(error_metrics_to_average), collapse = ", ")))
  } else {
    if (!is.character(overall.error.metrics)) {
      stop("'overall.error.metrics' must be a character vector of metric names.", call. = FALSE)
    }
    valid_options <- c("rmse", "mae", "mape")
    if(!all(overall.error.metrics %in% valid_options)) {
      stop(sprintf("Invalid metric(s) selected for 'overall.error.metrics': %s. Allowed composite metrics are: %s.",
                   paste(setdiff(overall.error.metrics, valid_options), collapse = ", "),
                   paste(valid_options, collapse = ", ")), call. = FALSE)
    }
    if (!all(overall.error.metrics %in% metrics)) {
      stop(sprintf("Metric(s) selected for 'overall.error.metrics' (%s) were not requested in the general 'metrics' argument. Please ensure all composite metrics are also in 'metrics'.",
                   paste(setdiff(overall.error.metrics, metrics), collapse = ", ")), call. = FALSE)
    }
    if (length(overall.error.metrics) < 1) {
      stop("At least one metric must be requested for 'overall.error.metrics' to calculate an overall composite score.", call. = FALSE)
    }
    error_metrics_to_average <- overall.error.metrics
  }

  # set seed for background sampling
  set.seed(seed)
  all_metrics <- list()

  # --- Loop through each dataset in test.data -----
  for (ds_name in names(test.data)) {
    current_data <- test.data[[ds_name]]
    current_sample_size <- nrow(current_data)
    loc <- sf::st_coordinates(current_data)[, c("X", "Y")]

    metrics_ds <- list(sample_size = current_sample_size)
    for (m in metrics) {
      metrics_ds[[m]] <- NA_real_
    }

    if (current_sample_size == 0) {
      message(sprintf("Skipping metrics for empty dataset: '%s'", ds_name))
      all_metrics[[ds_name]] <- metrics_ds
      next
    }

    is_count_data <- responseCounts %in% names(current_data)
    is_pa_data    <- responsePA %in% names(current_data)

    # --- Conditional continuous-outcome metrics calculation ---
    if (has_error_metrics && responseCounts %in% names(current_data)) {
      raw_expected <- terra::extract(expected.response, loc)[, 1]
      valid_idx    <- is.finite(raw_expected)

      if (any(valid_idx)) {
        if (isTRUE(is.pred.rate)) {
          if (is.null(exposure)) {
            warning("Predictions are marked as 'intensity', but no 'exposure' was provided. Assuming exposure is 1.", call. = FALSE)
            predicted <- raw_expected[valid_idx]
            observed  <- current_data[[responseCounts]][valid_idx]
          } else {
            predicted <- raw_expected[valid_idx]
            observed  <- current_data[[responseCounts]][valid_idx] / current_data[[exposure]][valid_idx]
          }
        } else {
          if (!is.null(exposure)) {
            warning("Predictions are not marked as 'intensity(rate)', but an 'exposure' was supplied. Assuming predictions are already on the count scale.", call. = FALSE)
          }
          predicted <- raw_expected[valid_idx]
          observed <- current_data[[responseCounts]][valid_idx]
        }

        if ("rmse" %in% metrics) metrics_ds$rmse <- rmse(observed, predicted)
        if ("mae" %in% metrics) metrics_ds$mae <- mae(observed, predicted)
        if ("mape" %in% metrics) metrics_ds$mape <- mape(observed, predicted)
        if ("r2" %in% metrics) metrics_ds$r2 <- r_squared(observed, predicted)
      } else {
        message(sprintf("Skipping continuous-outcome metrics for '%s' due to no valid predictions after filtering.", ds_name))
      }
    }

    # --- ROC-based metrics (AUC, TSS, Precision, Recall, Accuracy, etc.) ---
    if (has_roc_metrics) {
      raw_prob <- terra::extract(prob.raster, loc)[, 1]
      valid_idx <- is.finite(raw_prob)
      prob_roc <- raw_prob[valid_idx]

      eval_resp_roc <- NULL
      if (is_count_data) {
        eval_resp_roc <- ifelse(current_data[[responseCounts]] > 0, 1, 0)[valid_idx]
        metrics_ds$sample_size <- length(prob_roc)
      } else if (is_pa_data) {
        eval_resp_roc <- current_data[[responsePA]][valid_idx]
        metrics_ds$sample_size <- length(prob_roc)
      } else {
        bg_points_po <- sample_bg_points(mask = prob.raster, points = xy.excluded, n = n.background, xy = TRUE)
        if (!is.null(bg_points_po$bg)) {
          prob_bg_po <- terra::extract(prob.raster, bg_points_po$bg)[, 1]
          prob_bg_po <- prob_bg_po[is.finite(prob_bg_po)]

          if (length(prob_roc) > 0 && length(prob_bg_po) > 0) {
            eval_resp_roc <- c(rep(1, length(prob_roc)), rep(0, length(prob_bg_po)))
            prob_roc_po <- c(prob_roc, prob_bg_po)
            metrics_ds$sample_size <- length(prob_roc)
          } else {
            message(sprintf("Not enough valid presence (%s) or background (%s) points to compute ROC for '%s' (PO data).",
                            length(prob_roc), length(prob_bg_po), ds_name))
          }
        } else {
          message(sprintf("No background points sampled for '%s' (PO data). Cannot compute ROC metrics.", ds_name))
        }
      }

      # --- Check if we have valid data to compute ROC metrics ----
      if (!is.null(prob_roc) && length(unique(eval_resp_roc)) == 2 && length(prob_roc) > 0) {
        if (is_count_data || is_pa_data) {
          roc_obj <- pROC::roc(eval_resp_roc, prob_roc, quiet = TRUE)
        } else {
          roc_obj <- pROC::roc(eval_resp_roc, prob_roc_po, quiet = TRUE)
        }

        if ("auc" %in% metrics) {
          metrics_ds$auc <- as.numeric(roc_obj$auc)
        }
        if ("tss" %in% metrics) {
          metrics_ds$tss <- max(roc_obj$sensitivities + roc_obj$specificities - 1)
        }

        # List all metrics that coords() can return directly.
        metric_map_for_coords <- c(
          "threshold" = "threshold",
          "specificity" = "specificity",
          "recall" = "sensitivity", # sensitivity is recall
          "tpr" =  "sensitivity",   # sensisitivity is tpr
          "tnr" = "specificity",    # sepecificity is tnr
          "precision" = "ppv",      # ppv is precision
          "accuracy" = "accuracy",
          "npv" = "npv",
          "fpr" = "fpr",
          "fnr" = "fnr",
          "fdr" = "fdr"
        )
        coords_rets_to_pass <- unique(unname(metric_map_for_coords[names(metric_map_for_coords) %in% metrics]))
        coords_rets_to_pass <- coords_rets_to_pass[!is.na(coords_rets_to_pass) & !is.null(coords_rets_to_pass)]

        if (length(coords_rets_to_pass) > 0) {
          coords_x_val <- if (threshold.method == "best") "best" else fixed.threshold
          coords_results_df <- pROC::coords(roc = roc_obj,
                                            x = coords_x_val,
                                            ret = coords_rets_to_pass,
                                            best.method = best.method,
                                            transpose = FALSE, ...)

          temp_metrics_values <- list()
          for (coord_metric_name in names(metric_map_for_coords)) {
            pROC_col_name <- metric_map_for_coords[coord_metric_name]
            if (pROC_col_name %in% names(coords_results_df)) {
              temp_metrics_values[[coord_metric_name]] <- coords_results_df[[pROC_col_name]]
            } else {
              temp_metrics_values[[coord_metric_name]] <- rep(NA_real_, nrow(coords_results_df))
            }
          }

          if ("F1" %in% metrics) {
            temp_metrics_values$F1 <- 2 * (temp_metrics_values$precision * temp_metrics_values$recall) /
              (temp_metrics_values$precision + temp_metrics_values$recall)
            temp_metrics_values$F1[is.nan(temp_metrics_values$F1) | is.infinite(temp_metrics_values$F1)] <- -Inf
          }

          # Find row index for the chosen optimal threshold
          selected_idx <- 1
          if (threshold.method == "best" && nrow(coords_results_df) > 1) {
            if (best.threshold.policy == "first") {
              selected_idx <- 1
            } else if (best.threshold.policy == "last") {
              selected_idx <- nrow(coords_results_df)
            } else if (best.threshold.policy == "max.prec" && "precision" %in% names(temp_metrics_values)) {
              selected_idx <- which.max(temp_metrics_values$precision)
            } else if (best.threshold.policy == "max.recall" && "recall" %in% names(temp_metrics_values)) {
              selected_idx <- which.max(temp_metrics_values$recall)
            } else if (best.threshold.policy == "max.accu" && "accuracy" %in% names(temp_metrics_values)) {
              selected_idx <- which.max(temp_metrics_values$accuracy)
            } else if (best.threshold.policy == "max.f1" && "F1" %in% names(temp_metrics_values)) {
              selected_idx <- which.max(temp_metrics_values$F1)
            } else {
              message(sprintf("Policy '%s' cannot be applied for '%s' due to missing required metrics. Falling back to 'first' threshold.", best.threshold.policy, ds_name))
              selected_idx <- 1
            }
            if (length(selected_idx) > 1) {
              selected_idx <- selected_idx[1]
              message(sprintf("Multiple thresholds tied for '%s' with policy '%s'. Selecting the first one encountered.", ds_name, best.threshold.policy))
            }
          }

          for (metric_name_to_assign in names(temp_metrics_values)) {
            if (metric_name_to_assign %in% metrics) {
              metrics_ds[[metric_name_to_assign]] <- temp_metrics_values[[metric_name_to_assign]][selected_idx]
            }
          }
        }

      } else {
        message(sprintf("Skipping ROC-based metrics for '%s' due to single response class or no valid predictions.", ds_name))
      }
    }

    all_metrics[[ds_name]] <- metrics_ds
  }

  # --- Compute Weighted Composite Scores for requested metrics ----
  weighted_composite_scores <- list()

  for (metric_name in metrics) {
    if (metric_name == "threshold") {
      next
    }

    current_metric_values <- rep(NA_real_, length(test.data))
    current_sample_sizes <- rep(NA_real_, length(test.data))
    names(current_metric_values) <- names(test.data)
    names(current_sample_sizes) <- names(test.data)

    for (i in seq_along(test.data)) {
      ds_name <- names(test.data)[i]
      if (!is.null(all_metrics[[ds_name]]) && metric_name %in% names(all_metrics[[ds_name]])) {
        current_metric_values[ds_name] <- all_metrics[[ds_name]][[metric_name]]
        current_sample_sizes[ds_name] <- all_metrics[[ds_name]]$sample_size
      }
    }

    if (metric_name %in% error_metrics) {
      is_count_ds <- unlist(lapply(test.data, function(x) responseCounts %in% names(x)))
      names(is_count_ds) <- names(test.data)
      valid_indices <- !is.na(current_metric_values) &
        !is.na(current_sample_sizes) & current_sample_sizes > 0 & is_count_ds
    } else {
      valid_indices <- !is.na(current_metric_values) &
        !is.na(current_sample_sizes) & current_sample_sizes > 0
    }

    weighted_score <- NA_real_
    if (any(valid_indices)) {
      valid_metric_values <- current_metric_values[valid_indices]
      valid_weights <- current_sample_sizes[valid_indices]
      if (sum(valid_weights, na.rm = TRUE) > 0) {
        weighted_score <- sum(valid_metric_values * valid_weights) / sum(valid_weights, na.rm = TRUE)
      }
    }
    weighted_composite_scores[[paste0(toupper(metric_name), "_Comp")]] <- weighted_score
  }

  # Overall ROC-based Composite Score
  TOT_ROC_SCORE <- NA_real_
  relevant_roc_scores_names <- paste0(toupper(roc_metrics_to_average), "_Comp")
  available_relevant_roc_scores <- weighted_composite_scores[relevant_roc_scores_names]
  available_relevant_roc_scores <- unlist(available_relevant_roc_scores[!sapply(available_relevant_roc_scores, is.null)])
  available_relevant_roc_scores <- available_relevant_roc_scores[!is.na(available_relevant_roc_scores)]

  if (length(available_relevant_roc_scores) > 0) {
    TOT_ROC_SCORE <- sum(available_relevant_roc_scores) / length(available_relevant_roc_scores)
  } else {
    message("Cannot compute an overall ROC composite score as no relevant composite metrics were available or calculated.")
  }

  # Overall Continuous-outcome Composite Score
  TOT_ERROR_SCORE <- NA_real_
  relevant_cont_scores_names <- paste0(toupper(error_metrics_to_average), "_Comp")
  available_relevant_cont_scores <- weighted_composite_scores[relevant_cont_scores_names]
  available_relevant_cont_scores <- unlist(available_relevant_cont_scores[!sapply(available_relevant_cont_scores, is.null)])
  available_relevant_cont_scores <- available_relevant_cont_scores[!is.na(available_relevant_cont_scores)]

  if (length(available_relevant_cont_scores) > 0) {
    TOT_ERROR_SCORE <- sum(available_relevant_cont_scores) / length(available_relevant_cont_scores)
  } else {
    message("Cannot compute an overall continuous-outcome composite score as no relevant composite metrics were available or calculated.")
  }

  return_list <- list()
  for (ds_name in names(test.data)) {
    for (metric_name in metrics) {
      if (!is.null(all_metrics[[ds_name]]) && metric_name %in% names(all_metrics[[ds_name]])) {
        return_list[[paste0(toupper(metric_name), "_", ds_name)]] <- all_metrics[[ds_name]][[metric_name]]
      }
    }
  }

  return_list <- c(return_list, weighted_composite_scores)
  return_list$TOT_ROC_SCORE   <- TOT_ROC_SCORE
  return_list$TOT_ERROR_SCORE <- TOT_ERROR_SCORE

  return(return_list)
}


#--- 2) Functions to compute evaluation metrics for continuous-outcome responses ----

# RMSE (root mean square error),
rmse <- function(observed, predicted) {
  stopifnot(length(observed) == length(predicted))
  sqrt(mean((observed - predicted)^2, na.rm = TRUE))
}

# MAE (mean absolute error)
mae <- function(observed, predicted) {
  stopifnot(length(observed) == length(predicted))
  mean(abs(observed - predicted), na.rm = TRUE)
}

# MAPE (mean absolute percentage error) -- only valid for non-null response
mape <- function(observed, predicted) {
  stopifnot(length(observed) == length(predicted))
  100*mean(abs(observed - predicted)/observed, na.rm = TRUE)
}

# Pseudo R-squared
r_squared <- function(observed, predicted) {
  stopifnot(length(observed) == length(predicted))
  ss_res <- sum((observed - predicted)^2, na.rm = TRUE)
  ss_tot <- sum((observed - mean(observed, na.rm = TRUE))^2, na.rm = TRUE)
  1 - ss_res/ss_tot
}
