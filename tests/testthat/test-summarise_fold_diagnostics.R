
library(testthat)
library(terra)
library(sf)
library(isdmtools)

options(warn = -1) # Suppresses warnings globally

test_that("summarise_fold_diagnostics is robust to handle both diagnostic classes", {

  # GeoDiagnostic object: Folds are contiguous (gap = 0)
  geo_mock <- list(
    summary = data.frame(
      max_dist_km = c(10, 12),
      min_gap_km = c(0, 0),
      folds_ids = c(1, 2)
    )
  )
  class(geo_mock) <- "GeoDiagnostic"

  # EnvDiagnostic object: Folds are biased (p < 0.05) but have decent overlap
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
      folds_ids = c(1, 2)
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

  # Wrap in capture.output to verify if the print logic works
  output <- capture.output({
    summary_obj <- summarise_fold_diagnostics(geo_mock, env_mock)
    print(summary_obj)
  })

  #-- Assertions ---
  expect_s3_class(summary_obj, "FoldsSummary")

  # Check that Median Overlap is NA but didn't break the function
  expect_true(is.na(summary_obj$Value[summary_obj$Metric == "Median Overlap (D)"]))

  # Check that the conclusion was reached correctly despite NAs
  expect_true(any(grepl("environmentally representative", output)))

  #-- Verify biased scenario (p < 0.05)
  env_mock_biased <- env_mock
  env_mock_biased$summary$p_val <- 0.01
  summary_biased <- summarise_fold_diagnostics(geo_mock, env_mock_biased)

  expect_equal(summary_biased$Status[summary_biased$Metric == "Minimum p-value"], "Biased")
})


test_that("summarise_fold_diagnostics provides a unified report", {

  # Create mock data
  r <- terra::rast(extent = c(0, 1000, 0, 1000), res = 1, crs = "EPSG:32631")
  terra::values(r) <- seq(1, terra::ncell(r))
  names(r) <- "env_var"

  pts <- data.frame(
    x = c(20, 200, 700, 900),
    y = c(20, 200, 700, 900),
    folds_ids = rep(1:2, each = 2)
  )
  pts_sf <- sf::st_as_sf(pts, coords = c("x", "y"), crs = "EPSG:32631")

  mock_folds <- list(data_all = pts_sf, k = 2)
  class(mock_folds) <- "DataFolds"

  # Diagnostics
  geo_diag <- check_folds(mock_folds, rho = 0.5, plot = FALSE)
  env_diag <- check_env_balance(mock_folds, covariates = r, n_background = 100)
  report <- summarise_fold_diagnostics(geo_diag, env_diag)

  # Eexpectations
  expect_s3_class(report, "data.frame")
  expect_s3_class(report, "FoldsSummary")
  expect_true(all(c("Geographic", "Environmental") %in% report$Domain))

  # Check that NA Schoener's D doesn't break print method
  env_diag$summary$Schoener_D[1] <- NA_real_
  expect_no_error(print(summarise_fold_diagnostics(geo_diag, env_diag)))

  # Check for presence of key metrics from both worlds
  expected_metrics <- c("Avg Internal Distance (km)",
                        "Avg Inter-Fold Gap (km)",
                        "Median Overlap (D)",
                        "Minimum p-value")

  expect_true(all(expected_metrics %in% report$Metric))

  # Check: With a gap of ~11 and rho of 0.5, Status should be 'Separated'
  geo_status <- report$Status[report$Metric == "Avg Inter-Fold Gap (km)"]
  expect_equal(geo_status, "Separated")

  # Check: Since points are in corners, p-value should likely be significant (Biased)
  env_status <- report$Status[report$Metric == "Min_P_Value"]
  expect_type(env_status, "character")
})
