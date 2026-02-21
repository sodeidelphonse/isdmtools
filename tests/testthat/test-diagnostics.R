
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
                    folds_ids = rep(1:2, each = 5)
                   )
  pts_sf <- st_as_sf(pts, coords = c("v1", "v2"), crs = 32631)

  r <- rast(ext(0, 30, 0, 30), res = 1)
  values(r) <- runif(ncell(r))
  names(r) <- "env_var"
  crs(r) <- "EPSG:32631"

  # DataFolds object
  mock_folds <- list(data_all = pts_sf, k = 2)
  class(mock_folds) <- "DataFolds"

  # These points are > 10 units apart, so with rho=2, they should be Independent
  geo_res <- check_folds(mock_folds, rho = 2)

  expect_s3_class(geo_res, "GeoDiagnostic")

  expect_true("independence" %in% names(geo_res$summary))
  expect_true(all(c("x", "y") %in% names(geo_res$summary)))

  expect_gt(min(geo_res$summary$min_gap_km, na.rm = TRUE), 0)
})

test_that("check_spatial_geometry works with realistic integer fold IDs", {

  pts <- data.frame(
    v1 = c(1, 1.1, 5, 5.1),
    v2 = c(1, 1.1, 5, 5.1),
    folds_ids = c(1, 1, 2, 2)
  )
  pts_sf <- sf::st_as_sf(pts, coords = c("v1", "v2"), crs = 32631)

  # The function should handle the internal factor conversion
  res <- check_spatial_geometry(pts_sf, fold_col = "folds_ids", rho = 2, plot = TRUE)

  expect_type(res, "list")
  expect_true("max_dist_km" %in% colnames(res$summary))
  expect_s3_class(res$plot, "ggplot")

  # Check that folds_ids in summary is still readable: factor or character
  expect_equal(nrow(res$summary), 2)
})


test_that("EnvDiagnostic extracts values and runs stats", {

  # Create dummy datasets
  set.seed(123)
  pts <- data.frame(v1 = c(runif(10, 0, 10), runif(10, 20, 30)),
                    v2 = c(runif(10, 0, 10), runif(10, 20, 30)),
                    folds_ids = rep(1:2, each = 5)
                   )
  pts_sf <- st_as_sf(pts, coords = c("v1", "v2"), crs = "EPSG:32631")

  r <- rast(ext(0, 30, 0, 30), res = 1, crs = "EPSG:32631")
  values(r) <- runif(ncell(r))
  names(r) <- "env_var"

  mock_folds <- list(data_all = pts_sf, k = 2)
  class(mock_folds) <- "DataFolds"

  # Run the method
  env_res <- check_env_balance(mock_folds, covariates = r, n_background = 100)

  expect_s3_class(env_res, "EnvDiagnostic")
  expect_equal(nrow(env_res$summary), 1)
  expect_type(env_res$summary$p_val, "double")

  # Ensure background (BG) was included in the plot data
  plot_data <- ggplot2::ggplot_build(env_res$plot)$plot$data
  expect_true("BG" %in% plot_data$Fold)
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
    folds_ids = rep(1:2, each = 10)
  )
  pts_sf <- sf::st_as_sf(pts, coords = c("x", "y"), crs = "EPSG:32631")

  mock_folds <- list(data_all = pts_sf, k = 2)
  class(mock_folds) <- "DataFolds"

  # Diagnostic
  env_diag <- check_env_balance(mock_folds, covariates = r_cat, n_background = 100)

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
  pts <- data.frame(x = runif(10, 0, 2), y = runif(10, 0, 2), folds_ids = rep(1:2, each = 5))
  pts_sf <- sf::st_as_sf(pts, coords = c("x", "y"), crs = "EPSG:32631")

  mock_folds <- list(data_all = pts_sf)
  class(mock_folds) <- "DataFolds"

  # We expect low Schoener_D because the points don't represent the full raster range
  env_diag <- check_env_balance(mock_folds, covariates = r, n_background = 100)

  expect_s3_class(env_diag, "EnvDiagnostic")
  expect_true("Schoener_D" %in% names(env_diag$summary))
  expect_type(env_diag$summary$Schoener_D, "double")

  # Points in a corner of the map should have lower overlap
  expect_lt(env_diag$summary$Schoener_D, 0.9)
})


test_that("EnvDiagnostic correctly integrates background in summary and plot", {

  r <- terra::rast(extent = c(0, 10, 0, 10), res = 1, val = 1:100, crs = "EPSG:32631")
  names(r) <- "env_var"

  # Points concentrated in the 'low' values of the raster (values 1-20)
  pts <- data.frame(x = runif(10, 0, 2), y = runif(10, 0, 2), folds_ids = rep(1:2, each = 5))
  pts_sf <- st_as_sf(pts, coords = c("x", "y"), crs = "EPSG:32631")

  mock_folds <- list(data_all = pts_sf)
  class(mock_folds) <- "DataFolds"

  # We use a small n_background for speed in testing
  env_diag  <- check_env_balance(mock_folds, covariates = r, n_background = 100)
  plot_data <- env_diag$plot$data

  # Ensure background (BG) is present in the data used for the plot
  expect_true("BG" %in% unique(plot_data$Fold))

  # Since the points are in a corner, Schoener's D should be relatively low (< 1)
  expect_true("Schoener_D" %in% colnames(env_diag$summary))
  expect_lt(env_diag$summary$Schoener_D, 0.9)
})





