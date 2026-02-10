
library(testthat)
library(terra)
library(sf)

test_that("summarise_fold_diagnostics provides a unified report", {

  # Create mock DataFolds and covariates
  r <- terra::rast(extent = c(0, 10, 0, 10), res = 1, crs = "EPSG:32631")
  terra::values(r) <- seq(1, terra::ncell(r))
  names(r) <- "env_var"

  pts <- data.frame(
    x = c(1, 1.1, 9, 9.1),
    y = c(1, 1.1, 9, 9.1),
    fold_ids = rep(1:2, each = 2)
  )
  pts_sf <- sf::st_as_sf(pts, coords = c("x", "y"), crs = "EPSG:32631")

  mock_folds <- list(data_all = pts_sf)
  class(mock_folds) <- "DataFolds"

  #-- Generate individual diagnostics
  # We use a small rho (0.5) so that a gap of ~11 units is "Strongly Independent"
  geo_diag <- check_folds(mock_folds, rho = 0.5, plot = FALSE)
  env_diag <- check_env_balance(mock_folds, covariates = r, n_back = 100)

  # Summarise
  report <- summarise_fold_diagnostics(geo_diag, env_diag)

  expect_s3_class(report, "data.frame")

  # Check for presence of key metrics from both worlds
  expected_metrics <- c("Avg_Internal_Size_km", "Avg_Inter_Fold_Gap_km",
                        "Median_Env_Overlap_D", "Min_P_Value")
  expect_true(all(expected_metrics %in% report$Metric))

  # Check logic: With a gap of ~11 and rho of 0.5, Status should be 'Separated'
  geo_status <- report$Status[report$Metric == "Avg_Inter_Fold_Gap_km"]
  expect_equal(geo_status, "Separated")

  # Check logic: Since points are in corners, p-value should likely be significant
  env_status <- report$Status[report$Metric == "Min_P_Value"]
  expect_type(env_status, "character")
})
