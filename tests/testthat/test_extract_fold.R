
test_that("extract_fold correctly partitions data including NAs", {
  skip_if_not_installed("spatialsample")

  # Setup a controlled environment
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

  # Extract fold 1
  extracted <- extract_fold(folds, fold = 1)

  # Count observations in the resulting lists
  n_train <- nrow(extracted$train$PO)
  n_test  <- nrow(extracted$test$PO)

  # Training set points + testing set points = total points MINUS the points excluded.
  expect_equal(n_train + n_test, n_total - n_excluded_global)

  # Verify that the test set only contains points assigned to fold 1
  expect_equal(n_test, sum(folds$data_all$folds_ids == 1, na.rm = TRUE))
})

test_that("extract_fold handles multisource datasets correctly", {
  pts1 <- sf::st_as_sf(data.frame(x = 1:5, y = 1:5), coords = c("x", "y"), crs = 32631)
  pts2 <- sf::st_as_sf(data.frame(x = 10:14, y = 10:14), coords = c("x", "y"), crs = 32631)
  dsets <- list(Presence = pts1, Abundance = pts2)

  folds <- create_folds(dsets, k = 2, cv.method = "spatial")
  res <- extract_fold(folds, fold = 1)

  # Verify structure: train and test should both be lists of 2 elements
  expect_type(res$train, "list")
  expect_length(res$train, 2)
  expect_named(res$train, c("Presence", "Abundance"))
  expect_s3_class(res$train$Presence, "sf")
})


test_that("extract_fold preserves metadata and attributes", {

  # Create dummy datasets
  pts1 <- sf::st_as_sf(data.frame(x = 1:5, y = 1:5, env = runif(5)),
                       coords = c("x", "y"), crs = 32631)
  pts2 <- sf::st_as_sf(data.frame(x = 10:14, y = 10:14, count = rpois(5, 2)),
                       coords = c("x", "y"), crs = 32631)
  dsets <- list(Presence = pts1, Abundance = pts2)

  folds <- create_folds(dsets, k = 2, cv.method = "spatial")
  res   <- extract_fold(folds, fold = 1)

  # Check CRS preservation during filtering
  expect_equal(sf::st_crs(res$train$Presence), sf::st_crs(pts1))
  expect_equal(sf::st_crs(res$test$Abundance), sf::st_crs(pts2))

  # Check column preservation during filtering
  expect_true("env" %in% names(res$train$Presence))
  expect_true("count" %in% names(res$train$Abundance))

  # Check the list names (data fusion structure)
  expect_identical(names(res$train), names(dsets))
  expect_identical(names(res$test), names(dsets))
})

test_that("extract_fold handles empty subsets gracefully", {

  # Scenario: A dataset so small it only exists in one fold
  pts_large <- sf::st_as_sf(data.frame(x = 1:10, y = 1:10),
                            coords = c("x", "y"), crs = 32631)
  pts_small <- sf::st_as_sf(data.frame(x = 0.1, y = 0.1),
                            coords = c("x", "y"), crs = 32631)
  dsets <- list(Main = pts_large, Rare = pts_small)

  folds <- create_folds(dsets, k = 2, cv.method = "spatial")

  # Force the Rare point into Fold 1
  folds$data_all$folds_ids[folds$data_all$datasetName == "Rare"] <- 1

  # Extract Fold 2 (where 'Rare' should be absent from testing)
  res <- extract_fold(folds, fold = 2)

  # Rare should be in Training set for Fold 2, but empty in Testing for Fold 2
  expect_s3_class(res$test$Rare, "sf")
  expect_equal(nrow(res$test$Rare), 0)

  # Check that it still has the same geometry type even if empty
  expect_equal(sf::st_geometry_type(res$test$Rare)[0], sf::st_geometry_type(pts_small)[0])
})
