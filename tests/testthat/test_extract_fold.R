
test_that("extract_fold correctly partitions data including NAs", {
  skip_if_not_installed("spatialsample")

  #-- Setup a controlled environment
  # Points A & B are close (one will be excluded), C & D are far.
  pts <- sf::st_as_sf(data.frame(
    x = c(0, 1, 1000, 2000),
    y = c(0, 0, 0, 0)
  ), coords = c("x", "y"), crs = 32631)
  dsets <- list(PO = pts)

  # Create folds with a buffer that forces an NA (k = 2 and a buffer of 500m)
  folds <- create_folds(dsets, k = 2, cv.method = "buffer", radius = 1, buffer = 500)

  # Identify how many were excluded globally
  n_total <- nrow(pts)
  n_excluded_global <- sum(is.na(folds$data_all$folds_ids))

  #-- Extract fold 1
  extracted <- extract_fold(folds, fold = 1)

  # Count observations in the resulting lists
  n_train <- nrow(extracted$train$PO)
  n_test  <- nrow(extracted$test$PO)

  # The points in the training set + the points in the test set
  # should equal the total points MINUS the points that were excluded.
  expect_equal(n_train + n_test, n_total - n_excluded_global)

  # Verify that the test set only contains points assigned to fold 1
  expect_equal(n_test, sum(folds$data_all$folds_ids == 1, na.rm = TRUE))
})

test_that("extract_fold handles multisource data fusion correctly", {
  pts1 <- sf::st_as_sf(data.frame(x = 1:5, y = 1:5), coords = c("x", "y"), crs = 32631)
  pts2 <- sf::st_as_sf(data.frame(x = 10:14, y = 10:14), coords = c("x", "y"), crs = 32631)
  dsets <- list(Presence = pts1, Abundance = pts2)

  folds <- create_folds(dsets, k = 2, cv.method = "spatial")
  res <- extract_fold(folds, fold = 1)

  # Verify structure: train and test should both be lists of 2
  expect_type(res$train, "list")
  expect_length(res$train, 2)
  expect_named(res$train, c("Presence", "Abundance"))
  expect_s3_class(res$train$Presence, "sf")
})
