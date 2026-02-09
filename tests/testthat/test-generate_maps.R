
library(testthat)
library(ggplot2)
library(sf)
library(terra)

test_that("generate_maps handles single panel correctly", {
  df <- data.frame(x = 1:5, y = 1:5, mean = 1)
  p <- generate_maps(df, var_names = "mean")

  expect_s3_class(p$data, "data.frame")
  expect_true("mean" %in% names(p$data))
  expect_s3_class(p$facet, "FacetNull")
})


test_that("generate_maps handles multi-panel map via melting", {
  df <- data.frame(x = 1:5, y = 1:5, mean = 1, sd = 0.5)
  p <- generate_maps(df, var_names = c("mean", "sd"))

  expect_true(all(c("prediction_var", "value") %in% names(p$data)))
  expect_equal(levels(p$data$prediction_var), c("mean", "sd"))
  expect_s3_class(p$facet, "FacetWrap")
})

test_that("generate_maps respects custom panel labels", {
  df <- data.frame(x = 1:5, y = 1:5, mean = 1, sd = 0.5)
  labels <- c("Suitability", "StDev")
  p <- generate_maps(df, var_names = c("mean", "sd"), panel_labels = labels)

  expect_equal(levels(p$data$prediction_var), labels)
})

test_that("generate_maps handles SpatRaster through data conversion", {
  r <- rast(nrows = 5, ncols = 5, vals = 1:25)
  names(r) <- "mean"

  p <- generate_maps(r, var_names = "mean")
  expect_true(all(c("x", "y", "mean") %in% names(p$data)))
  expect_no_error(ggplot2::ggplot_build(p))
})

test_that("generate_maps handles base.map and CRS matching", {
  pts <- st_as_sf(data.frame(x = 100, y = 100, mean = 1), coords = c("x", "y"), crs = 32631)
  poly <- st_as_sfc("POLYGON ((0 0, 1 0, 1 1, 0 1, 0 0))", crs = 4326)
  base_map <- st_sf(geometry = poly)

  expect_no_error(generate_maps(pts, var_names = "mean", base_map = base_map))
})
