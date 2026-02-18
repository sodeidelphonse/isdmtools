
library(testthat)
library(terra)
library(sf)
library(isdmtools)

options(warn = -1)

test_that("compute_metrics handles Presence-Absence data", {

  # A small raster of probabilities (0 to 1)
  r_prob <- terra::rast(ncols = 10, nrows = 10, vals = runif(100), crs = "epsg:4326")
  names(r_prob) <- "prob"

  # One fold with presence (1) and absence (0)
  pts <- sf::st_as_sf(data.frame(
    x = c(0.1, 1, 0.5, 0.2),
    y = c(0.1, 1, 0.5, 0.8),
    present = c(1, 1, 0, 0)
  ), coords = c("x", "y"), crs = 4326)

  test_data_mock <- list(Data1 = pts)

  results <- compute_metrics(
    test_data = test_data_mock,
    prob_raster = r_prob,
    response_pa = "present",
    metrics = c("auc", "tss", "f1"),
    threshold_method = "best",
    best_method = "youden",
    n_background = 50
  )

  expect_type(results, "list")
  expect_s3_class(results, "ISDMmetrics")
  expect_output(print(results))
  expect_true(any(grepl("AUC_Data1|TSS_Data1|F1_Data1", names(results), ignore.case = TRUE)))

  # Check if background attribute is NULL
  bg_attr <- get_background(results)
  expect_type(bg_attr, "NULL")
})


test_that("compute_metrics handles Presence-Only data and background points", {

  r_prob <- terra::rast(ncols = 10, nrows = 10, vals = runif(100), crs = "epsg:4326")
  names(r_prob) <- "prob"

  # Function identifies PO and triggers background sampling.
  pts <- sf::st_as_sf(data.frame(
    x = seq(0.1, 0.5, length.out = 5),
    y = seq(0.1, 0.5, length.out = 5)
  ), coords = c("x", "y"), crs = 4326)

  test_data_mock <- list(PO_Data = pts)

  results <- compute_metrics(
    test_data = test_data_mock,
    prob_raster = r_prob,
    metrics = c("auc", "tss", "f1"),
    n_background = 50,
    seed = 123
  )
  bg_attr <- get_background(results)

  expect_s3_class(results, "ISDMmetrics")
  expect_s3_class(bg_attr, "BackgroundPoints")
  expect_equal(nrow(bg_attr$bg), 50)
  expect_output(print(results))

  expect_true(is.null(results$MAE))
  expect_true(any(grepl("AUC|TSS|F1", names(results))))
})


test_that("compute_metrics handles Count data and Error scores", {

  r_exp <- terra::rast(ncols = 10, nrows = 10, vals = seq(1, 5, length.out = 100),
                       crs = "epsg:4326")
  names(r_exp) <- "intensity"

  pts <- sf::st_as_sf(data.frame(
    x = c(0.2, 0.8, 1.2, 0.8),
    y = c(0.2, 0.8, 1.2, 1.5),
    counts = c(2, 5)
  ), coords = c("x", "y"), crs = 4326)
  test_data_mock <- list(CountData = pts)

  results <- compute_metrics(
    test_data = test_data_mock,
    expected_response = r_exp,
    metrics = c("rmse", "mae"),
    response_counts = "counts"
  )

  expect_false(is.na(results$TOT_ERROR_SCORE))
  expect_true(any(grepl("RMSE|MAE", names(results))))
  expect_true(is.na(results$TOT_ROC_SCORE))
  expect_true(is.null(results$TSS))
})


test_that("compute_metrics handles Presence-Absence with custom label name", {
  r_prob <- terra::rast(ncols = 10, nrows = 10, vals = runif(100),
                        crs = "epsg:4326")

  # Standardized PA data
  pts1 <- sf::st_as_sf(data.frame(
    x = c(0.1, 1, 0.5, 0.2),
    y = c(0.2, 1, 0.5, 0.8),
    status = c(1, 1, 0, 0)
  ), coords = c("x", "y"), crs = 4326)

  pts2 <- sf::st_as_sf(data.frame(
    x = c(0.1, 1.3, 4.5, 0.2),
    y = c(0.2, 1.8, 5.5, 0.8),
    status = c(1, 1, 1, 0)
  ), coords = c("x", "y"), crs = 4326)

  test_data = list(Data1_PA = pts1, Data2_PA = pts2)

  results <- compute_metrics(
    test_data = test_data,
    prob_raster = r_prob,
    response_pa = "status",
    metrics = c("auc", "precision")
    )

  # Should be NULL for PA data
  expect_type(get_background(results), "NULL")

  expect_true(any(grepl("AUC_Data1|PRECISION_Data1|AUC_Data2|PRECISION_Data2", names(results),
                        ignore.case = TRUE))
             )

  expect_true(any(grepl("AUC_Comp|PRECISION_Comp|TOT_ROC_SCORE|TOT_ERROR_SCORE", names(results),
                        ignore.case = TRUE))
             )

  expect_false("TSS_Data1" %in% names(results))
  expect_false("RMSE_Data1" %in% names(results))

  summ_out <- capture.output(print(results))
  expect_true(any(grepl("Datasets Evaluated: Data1_PA, Data2_PA", summ_out)))

  # We expect TOT_ROC_SCORE to be 0.500 as PREC is excluded
  expect_equal(results$TOT_ROC_SCORE, 0.500, tolerance = 0.001)
})

test_that("compute_metrics handles all valid overall ROC options", {
  r_prob <- terra::rast(ncols = 10, nrows = 10, vals = runif(100), crs = "epsg:4326")

  pts1 <- sf::st_as_sf(data.frame(
    x = c(0.1, 1, 0.5, 0.2),
    y = c(0.2, 1, 0.5, 0.8),
    status = c(1, 1, 0, 1)
  ), coords = c("x", "y"), crs = 4326)

  # Only request some robust ROC-based metrics for overall score
  valid_metrics <- c("auc", "accuracy", "f1")

  results <- compute_metrics(
    test_data = list(Data1 = pts1),
    prob_raster = r_prob,
    response_pa = "status",
    metrics = valid_metrics
  )

  metric_names <- names(results)
  expect_true(all(sapply(toupper(valid_metrics), function(m) any(grepl(m, metric_names)))))

  # TOT_ROC_SCORE should be a real number, not NA
  expect_false(is.na(results$TOT_ROC_SCORE))
})
