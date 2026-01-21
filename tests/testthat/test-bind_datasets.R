
test_that("bind_datasets standardizes dataset Names and handles sf points", {
  p1 <- sf::st_as_sf(data.frame(x = 1, y = 1), coords = c("x","y"))
  p2 <- sf::st_as_sf(data.frame(x = 2, y = 2), coords = c("x","y"))
  datasets <- list(Presence = p1, Count = p2)

  result <- bind_datasets(datasets)

  expect_true("datasetName" %in% names(result))
  expect_s3_class(result$datasetName, "factor")
  expect_equal(levels(result$datasetName), c("Presence", "Count"))
})
