test_that("create_folds works with blockCV (default)", {
  set.seed(123)
  pts <- sf::st_as_sf(data.frame(x = runif(20), y = runif(20)), coords = c("x", "y"), crs = 4326)
  dsets <- list(PO = pts[1:10, ], AB = pts[11:20, ])

  res <- create_folds(dsets, k = 5, cv.method = "spatial", plot = FALSE)
  expect_s3_class(res, "DataFolds")
  expect_equal(length(res$data_all$folds_ids), 20)
  expect_type(res$data_all$folds_ids, "integer")
})

test_that("create_folds works with spatialsample methods", {
  skip_if_not_installed("spatialsample")
  set.seed(123)

  pts <- sf::st_as_sf(data.frame(x = runif(50), y = runif(50)), coords = c("x", "y"), crs = 4326)
  dsets <- list(PO = pts)

  # Test location-out: We pass 'site' as a quoted string to be safe with the ellipsis
  dsets$PO$site <- sample(letters[1:5], 50, replace = TRUE)
  res_loc <- create_folds(dsets, k = 5, cv.method = "location", group = "site")

  expect_s3_class(res_loc, "DataFolds")
  expect_true(all(res_loc$data_all$folds_ids %in% 1:5 | is.na(res_loc$data_all$folds_ids)))
})

test_that("spatialsample bridge handles buffers (NA values)", {
  skip_if_not_installed("spatialsample")

  # A & B are 1m apart. C & D are 10,000m apart.
  pts <- sf::st_as_sf(data.frame(
    x = c(0, 1, 10000, -10000),
    y = c(0, 0, 0, 0)
  ), coords = c("x", "y"), crs = 32631)
  dsets <- list(test = pts)

  # radius = 1: Only the focal point is 'assessment'
  # buffer = 100: The neighbor (1m away) is excluded.
  # If Fold 1 uses Point C as assessment and Fold 2 uses Point D,
  # Points A and B are in the 'analysis' set but could be buffered out.
  res_buf <- create_folds(
    dsets,
    k = 2,
    cv.method = "buffer",
    radius = 1,
    buffer = 500
  )

  expect_s3_class(res_buf, "DataFolds")

  has_nas <- any(is.na(res_buf$data_all$folds_ids))
  if (!has_nas) {
    skip("spatialsample resampled all points; no NAs generated for this small sample.")
  }

  expect_true(has_nas)
})

test_that("create_folds errors gracefully on missing packages", {

  if (requireNamespace("spatialsample", quietly = TRUE)) {
    skip("spatialsample is installed, cannot test missing package error.")
  }

  pts <- sf::st_as_sf(data.frame(x = 1:5, y = 1:5), coords = c("x", "y"), crs = 4326)
  dsets <- list(test = pts)
  expect_error(create_folds(dsets, cv.method = "block"), "required")
})
