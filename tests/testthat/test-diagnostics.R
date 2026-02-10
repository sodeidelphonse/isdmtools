
library(testthat)
library(sf)
library(terra)
library(isdmtools)

options(warn = -1) # Suppresses warnings globally in the current session

test_that("GeoDiagnostic identifies distance categories correctly", {

  # Create dummy datasets
  set.seed(123)
  pts <- data.frame(v1 = c(runif(10, 0, 10), runif(10, 20, 30)),
                    v2 = c(runif(10, 0, 10), runif(10, 20, 30)),
                    fold_ids = rep(1:2, each = 5)
                   )
  pts_sf <- st_as_sf(pts, coords = c("v1", "v2"), crs = 32631)

  r <- rast(ext(0, 30, 0, 30), res = 1)
  values(r) <- runif(ncell(r))
  names(r) <- "env_var"

  # DataFolds object
  mock_folds <- list(data_all = pts_sf, k = 2)
  class(mock_folds) <- "DataFolds"

  # These points are > 10 units apart, so with rho=2, they should be Independent
  geo_res <- check_folds(mock_folds, rho = 2)

  expect_s3_class(geo_res, "GeoDiagnostic")
  expect_true("independence" %in% names(geo_res$summary))

  expect_true(all(c("x", "y") %in% names(geo_res$summary)))

  # Check if it detects a gap > 0
  expect_gt(min(geo_res$summary$min_gap_km, na.rm = TRUE), 0)
})

test_that("check_spatial_geometry works with realistic integer fold IDs", {

  pts <- data.frame(
    v1 = c(1, 1.1, 5, 5.1),
    v2 = c(1, 1.1, 5, 5.1),
    fold_ids = c(1, 1, 2, 2)
  )
  pts_sf <- sf::st_as_sf(pts, coords = c("v1", "v2"), crs = 32631)

  # The function should handle the internal factor conversion
  res <- check_spatial_geometry(pts_sf, fold_col = "fold_ids", rho = 2, plot = TRUE)

  expect_type(res, "list")
  expect_true("max_dist_km" %in% colnames(res$summary))
  expect_s3_class(res$plot, "ggplot")

  # Check that fold_ids in summary is still readable: factor or character
  expect_equal(nrow(res$summary), 2)
})


test_that("EnvDiagnostic extracts values and runs stats", {

  # Create dummy datasets
  set.seed(123)
  pts <- data.frame(v1 = c(runif(10, 0, 10), runif(10, 20, 30)),
                    v2 = c(runif(10, 0, 10), runif(10, 20, 30)),
                    fold_ids = rep(1:2, each = 5)
                   )
  pts_sf <- st_as_sf(pts, coords = c("v1", "v2"), crs = 32631)

  r <- rast(ext(0, 30, 0, 30), res = 1)
  values(r) <- runif(ncell(r))
  names(r) <- "env_var"

  mock_folds <- list(data_all = pts_sf, k = 2)
  class(mock_folds) <- "DataFolds"

  # Run the method
  env_res <- check_env_balance(mock_folds, covariates = r, n_background = 100)

  expect_s3_class(env_res, "EnvDiagnostic")
  expect_equal(nrow(env_res$summary), 1)
  expect_type(env_res$summary$p_val, "double")

  # Ensure background was included in the plot data
  plot_data <- ggplot2::ggplot_build(env_res$plot)$plot$data
  expect_true("Background" %in% plot_data$Fold)
})


test_that("EnvDiagnostic handles categorical variables correctly", {

  r_cat <- terra::rast(extent = c(0, 10, 0, 10), res = 1, crs = "EPSG:32631")

  # Create 3 habitat classes
  terra::values(r_cat) <- sample(c(1, 2, 3), terra::ncell(r_cat), replace = TRUE)
  levels(r_cat) <- data.frame(ID = 1:3, habitat = c("Forest", "Grassland", "Urban"))
  names(r_cat)  <- "land_cover"

  pts <- data.frame(
    x = runif(20, 0, 10),
    y = runif(20, 0, 10),
    fold_ids = rep(1:2, each = 10)
  )
  pts_sf <- sf::st_as_sf(pts, coords = c("x", "y"), crs = "EPSG:32631")

  mock_folds <- list(data_all = pts_sf, k = 2)
  class(mock_folds) <- "DataFolds"

  # Diagnostic
  env_diag <- suppressWarnings(
    check_env_balance(mock_folds, covariates = r_cat, n_background = 100)
  )

  # Expectations
  expect_s3_class(env_diag, "EnvDiagnostic")
  expect_equal(env_diag$summary$Type, "Categorical")

  # Schoener's D should be NA for categorical variables
  expect_true(is.na(env_diag$summary$Schoener_D))

  # P-value should still be calculated via Chi-squared statistic
  expect_type(env_diag$summary$p_val, "double")
})


test_that("Environmental overlap logic (Schoener's D) works correctly", {

  #--- Test the internal helper function
  # Create two identical distributions (should have high overlap)
  v1 <- rnorm(1000, mean = 10, sd = 1)
  v2 <- rnorm(1000, mean = 10, sd = 1)

  overlap_high <- isdmtools:::calc_niche_overlap(v1, v2)
  expect_gt(overlap_high, 0.8)

  # Create two divergent distributions (should have low overlap)
  v3 <- rnorm(1000, mean = 20, sd = 1)
  overlap_low <- isdmtools:::calc_niche_overlap(v1, v3)
  expect_lt(overlap_low, 0.2)

  #--- Test via check_env_balance method
  r <- terra::rast(extent = c(0, 10, 0, 10), res = 1, crs = "EPSG:32631")
  terra::values(r) <- seq(1, terra::ncell(r))
  names(r) <- "temp"

  # Create points that only exist in the 'low' values of the raster
  pts <- data.frame(x = runif(10, 0, 2), y = runif(10, 0, 2), fold_ids = rep(1:2, each = 5))
  pts_sf <- sf::st_as_sf(pts, coords = c("x", "y"), crs = "EPSG:32631")

  mock_folds <- list(data_all = pts_sf)
  class(mock_folds) <- "DataFolds"

  # We expect low Schoener_D because the points don't represent the full raster range
  env_diag <- check_env_balance(mock_folds, covariates = r, n_background = 100)

  expect_s3_class(env_diag, "EnvDiagnostic")
  expect_true("Schoener_D" %in% names(env_diag$summary))
  expect_type(env_diag$summary$Schoener_D, "double")

  # Points in a corner of the map should have lower overlap
  expect_lt(env_diag$summary$Schoener_D, 1.0)
})


test_that("EnvDiagnostic correctly integrates background in summary and plot", {

  r <- terra::rast(extent = c(0, 10, 0, 10), res = 1, val = 1:100)
  names(r) <- "env_var"

  # Points concentrated in the 'low' values of the raster (values 1-20)
  pts <- data.frame(x = runif(10, 0, 2), y = runif(10, 0, 2), fold_ids = rep(1:2, each = 5))
  pts_sf <- st_as_sf(pts, coords = c("x", "y"), crs = "EPSG:32631")

  mock_folds <- list(data_all = pts_sf)
  class(mock_folds) <- "DataFolds"

  # We use a small n_background for speed in testing
  env_diag  <- check_env_balance(mock_folds, covariates = r, n_background = 100)
  plot_data <- env_diag$plot$data

  # Ensure "Background" is present in the data used for the plot
  expect_true("Background" %in% unique(plot_data$Fold))

  # Since the points are in a corner, Schoener's D should be relatively low (< 1)
  expect_true("Schoener_D" %in% colnames(env_diag$summary))
  expect_lt(env_diag$summary$Schoener_D, 0.9)
})


test_that("summarise_fold_diagnostics is robust", {

  #--- GeoDiagnostic object: Folds are contiguous (gap = 0)
  geo_mock <- list(
    summary = data.frame(
      max_dist_km = c(10, 12),
      min_gap_km = c(0, 0),
      fold_ids = c(1, 2)
    )
  )
  class(geo_mock) <- "GeoDiagnostic"

  #-- EnvDiagnostic object: Folds are biased (p < 0.05) but have decent overlap
  env_mock <- list(
    summary = data.frame(
      Variable = "temp",
      Type = "Continuous",
      p_val = 0.01,
      Schoener_D = 0.75
    )
  )

  class(env_mock) <- "EnvDiagnostic"
  res <- summarise_fold_diagnostics(geo_mock, env_mock)

  # Assertions
  expect_s3_class(res, "FoldsSummary")
  expect_s3_class(res, "data.frame")

  # Did it catch the 'Contiguous' status from the Geo mock?
  geo_status <- res$Status[res$Domain == "Geographic"]
  expect_true(any(geo_status == "Contiguous"))

  # Did it catch the 'Biased' status from the Env mock?
  env_status <- res$Status[res$Metric == "Minimum p-value"]
  expect_equal(env_status, "Biased")

  # Ensure value rounding is working
  expect_equal(res$Value[res$Metric == "Median Overlap (D)"], 0.750)
})


test_that("summarise_fold_diagnostics handles categorical-only NAs gracefully", {

  # Mock GeoDiagnostic (spatially independent)
  geo_mock <- list(
    summary = data.frame(
      max_dist_km = c(50, 60),
      min_gap_km = c(10, 10),
      fold_ids = c(1, 2)
    )
  )
  class(geo_mock) <- "GeoDiagnostic"

  # Mock EnvDiagnostic (categorical only)
  env_mock <- list(
    summary = data.frame(
      Variable = "land_cover",
      Type = "Categorical",
      p_val = 0.08,         # Balanced
      Schoener_D = NA_real_ # Expected for categorical
    )
  )
  class(env_mock) <- "EnvDiagnostic"

  # Wrap in capture.output to verify the print logic doesn't crash
  output <- capture.output({
    summary_obj <- summarise_fold_diagnostics(geo_mock, env_mock)
    print(summary_obj)
  })

  #-- Assertions ---
  expect_s3_class(summary_obj, "FoldsSummary")

  # Check that Median Overlap is NA but didn't break the function
  expect_true(is.na(summary_obj$Value[summary_obj$Metric == "Median Overlap (D)"]))

  # Check that the conclusion was reached correctly despite NAs
  # Since p > 0.05 and gap > 0, it should be robust.
  expect_true(any(grepl("environmentally representative", output)))

  #-- Verify biased scenario (p < 0.05)
  env_mock_biased <- env_mock
  env_mock_biased$summary$p_val <- 0.01
  summary_biased <- summarise_fold_diagnostics(geo_mock, env_mock_biased)

  expect_equal(summary_biased$Status[summary_biased$Metric == "Minimum p-value"], "Biased")
})

test_that("summarise_fold_diagnostics provides a unified report", {

  r <- terra::rast(extent = c(0, 10, 0, 10), res = 1, crs = "EPSG:32631")
  terra::values(r) <- 1:100
  names(r) <- "env_var"

  # Create points with explicit names used in st_as_sf
  pts <- data.frame(
    x = c(1, 1.1, 9, 9.1),
    y = c(1, 1.1, 9, 9.1),
    fold_ids = factor(rep(1:2, each = 2))
  )
  pts_sf <- sf::st_as_sf(pts, coords = c("x", "y"), crs = "EPSG:32631")

  mock_folds <- list(data_all = pts_sf)
  class(mock_folds) <- "DataFolds"

  # Ensure check_folds is called on the DataFolds object
  geo_diag <- check_folds(mock_folds, rho = 0.5, plot = FALSE)
  env_diag <- check_env_balance(mock_folds, covariates = r, n_background = 100)

  # TEST Constructor
  report <- summarise_fold_diagnostics(geo_diag, env_diag)

  expect_s3_class(report, "FoldsSummary")
  expect_true(all(c("Geographic", "Environmental") %in% report$Domain))

  # Check that NA Schoener's D doesn't break print method
  env_diag$summary$Schoener_D[1] <- NA_real_
  expect_no_error(print(summarise_fold_diagnostics(geo_diag, env_diag)))
})
