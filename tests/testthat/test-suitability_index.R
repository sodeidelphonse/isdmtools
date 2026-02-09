library(testthat)
library(terra)

test_that("suitability_index handles data.frame and SpatRaster inputs", {

  # Create dummy data
  df <- data.frame(x = c(0, 1), y = c(0, 1), mean = c(0, 1))
  projection <- "+proj=utm +zone=31 +units=km"

  # Test data.frame input
  res_df <- suitability_index(df, post_stat = "mean", projection = projection)
  expect_s4_class(res_df, "SpatRaster")

  # Test SpatRaster input
  r <- terra::rast(df, type = "xyz", crs = projection)
  res_rast <- suitability_index(r, post_stat = "mean")
  expect_s4_class(res_rast, "SpatRaster")

  expect_equal(terra::values(res_df), terra::values(res_rast))
})

test_that("scale_independent argument correctly modifies probability", {

  # Set up a raster with 10km2 cells (res = sqrt(10))
  side <- sqrt(10)
  r <- terra::rast(nrows=1, ncols=1, xmin=0, xmax=side, ymin=0, ymax=side)
  r[]      <- 0  # eta = 0, so exp(eta) = 1
  names(r) <- "mean"
  terra::crs(r) <- "+proj=utm +zone=31 +units=km"

  # Default (scale_independent = FALSE): p = 1 - exp(-10 * 1)
  p_scaled <- suitability_index(r, response_type = "po", scale_independent = FALSE)
  expect_equal(as.numeric(terra::values(p_scaled)), 1 - exp(-10), tolerance = 1e-6)

  # Scale Independent: p = 1 - exp(-1 * 1)
  p_indep <- suitability_index(r, response_type = "po", scale_independent = TRUE)
  expect_equal(as.numeric(terra::values(p_indep)), 1 - exp(-1), tolerance = 1e-6)
})



test_that("error handling triggers for missing columns/layers", {

  # Case 1: Data frame with wrong column name
  df <- data.frame(x = c(0, 1), y = c(0, 1), wrong_name = c(0, 1))
  expect_error(
    suitability_index(df, post_stat = "mean"),
    "The following 'post_stat' columns were not found"
  )

  # Case 2: Raster with wrong layer name
  r <- terra::rast(nrows = 5, ncols = 5)
  names(r) <- "wrong_name"
  expect_error(
    suitability_index(r, post_stat = "mean"),
    "The following 'post_stat' layers were not found"
  )
})

test_that("suitability_index handles lowercase x and y", {

  # Use a simple planar projection
  proj_planar <- "+proj=utm +zone=31 +units=km"

  # Testing lowercase standardization
  df <- data.frame(x = c(0, 1), y = c(0, 1), mean = c(0, 1))
  res <- suitability_index(df, post_stat = "mean", projection = proj_planar)
  expect_s4_class(res, "SpatRaster")

  # Testing uppercase X/Y auto-correction
  df_upper <- data.frame(X = c(0, 1), Y = c(0, 1), mean = c(0, 1))
  res_upper <- suitability_index(df_upper, post_stat = "mean", projection = proj_planar)
  expect_s4_class(res_upper, "SpatRaster")
  expect_equal(terra::values(res), terra::values(res_upper))
})

test_that("output_format returns correct scales", {
  r <- terra::rast(nrows = 5, ncols = 5, val = 2)
  names(r) <- "mean"

  # Linear should be exactly eta = 2
  expect_equal(unique(as.numeric(terra::values(suitability_index(r, output_format = "linear")))), 2)

  # Response should be exp(eta) regardless of scaling factor
  expect_equal(unique(as.numeric(terra::values(suitability_index(r, output_format = "response")))), exp(2))
})


test_that("longlat area calculation uses km2 for PO intensity scaling", {

  # Define a resolution roughly equivalent to 10km2 at 45 degrees North
  # 0.033 degrees is approx 3.6km. 3.6 * (3.6 * cos(45)) approx 10km2.
  res_deg <- 0.033
  r <- terra::rast(nrows = 10, ncols = 10,
                   xmin = 0, xmax = 10 * res_deg,
                   ymin = 45, ymax = 45 + 10 * res_deg)
  r[] <- -2  # constant log-intensity
  names(r) <- "mean"
  terra::crs(r) <- "+proj=longlat +datum=WGS84"

  # Calculate suitability
  res <- suitability_index(r, post_stat = "mean", response_type = "po", output_format = "prob")

  # Manual validation logic
  val <- as.numeric(terra::values(res)[1])
  area_km2 <- as.numeric(terra::cellSize(r, unit = "km")[1])

  # Ensure the area is in the expected 10km2 ballpark for this resolution
  expect_gt(area_km2, 8)
  expect_lt(area_km2, 12)

  # Verify the probability calculation matches the IPP formula
  expected_p <- 1 - exp(-(area_km2 * exp(-2)))
  expect_equal(val, expected_p, tolerance = 1e-6)
})

