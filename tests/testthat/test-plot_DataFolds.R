
library(testthat)
library(ggplot2)
library(sf)

test_that("plot.DataFolds creates the long data.frame correclty", {
  pts <- sf::st_as_sf(data.frame(x = 1:3, y = 1:3), coords = c("x", "y"), crs = 4326)
  folds <- create_folds(list(D1 = pts), k = 2, cv.method = "spatial", plot = FALSE)
  folds$data_all$folds_ids[1] <- c(1, 2, NA) # simulate a buffered point
  p <- plot(folds)

  expect_s3_class(p$data, "data.frame")
  expect_s3_class(p$data$fold_panel, "factor")
  expect_s3_class(p$data$Partition, "factor")
  expect_s3_class(p$facet, "FacetWrap")

  expect_no_error(ggplot2::ggplot_build(p))
  expect_true(all(c("Partition", "datasetName") %in% names(p$data)))
  expect_equal(levels(p$data$Partition), c("Train", "Test", "Excluded"))
})


test_that("plot.DataFolds uses correct color mapping", {
  pts <- sf::st_as_sf(data.frame(x = 1:3, y = 1:3), coords = c("x", "y"), crs = 4326)
  folds <- create_folds(list(D1 = pts), k = 2, cv.method = "spatial", plot = FALSE)
  folds$data_all$folds_ids[1:3] <- c(1, 2, NA)

  p_built <- ggplot2::ggplot_build(plot(folds))
  point_data <- p_built$data[[1]]
  actual_colors <- unique(point_data$colour)

  expect_equal(length(actual_colors), 3)
  expect_true("blue" %in% actual_colors)   # Train
  expect_true("orange" %in% actual_colors) # Test
  expect_true("grey" %in% actual_colors)   # Excluded
})
