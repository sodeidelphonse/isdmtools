library(testthat)
library(terra)
library(sf)

test_that("sample_background returns correct S3 object and sample size", {

  mask <- terra::rast(nrows=10, ncols=10, xmin=0, xmax=10, ymin=0, ymax=10, vals=1)
  n_sample <- 20
  bg_obj <-  suppressWarnings(sample_background(mask, n = n_sample, method = "random"))

  expect_s3_class(bg_obj, "BackgroundPoints")
  expect_named(bg_obj, c("bg", "mask"))
  expect_equal(nrow(bg_obj$bg), n_sample)
})


test_that("sample_background excludes the specific cell containing the presence", {

  mask <- terra::rast(nrows=5, ncols=5, xmin=0, xmax=5, ymin=0, ymax=5, vals=1)

  # Pick a point that is not at the center and record its cell
  pres_sf <- sf::st_as_sf(data.frame(x = 0.2, y = 4.8), coords = c("x", "y"))
  target_cell <- terra::cellFromXY(mask, matrix(c(0.2, 4.8), ncol=2))

  bg_obj <- suppressWarnings(sample_background(mask, points = pres_sf, n = 24))

  # Check if it actually returned the object and not NULL
  expect_s3_class(bg_obj, "BackgroundPoints")

  if (!is.null(bg_obj)) {
    bg_cells <- terra::cellFromXY(mask, as.matrix(bg_obj$bg[, c("x", "y")]))

    # The target_cell should NOT be in the background points cells (bg_cells)
    expect_false(target_cell %in% bg_cells,
                 info = paste("Cell", target_cell, "containing the presence was sampled as background!"))

    expect_true(is.na(bg_obj$mask[target_cell][1,1]))
  }
})


test_that("sample_background handles method arguments correctly", {
  mask <- terra::rast(nrows=10, ncols=10, vals=1)

  expect_no_error(
    bg_reg <- suppressWarnings(sample_background(mask, n = 10, method = "regular"))
  )
})
