
#--- Testing prepare_predictions ----------

library(testthat)
library(sf)
library(terra)

test_that("prepare_predictions handles sf inputs and converts to clean x/y", {

  # Case 1: Test standard sf to data.frame conversion
  pts <- st_as_sf(data.frame(X = c(1, 2), Y = c(3, 4), val = c(10, 20)),
                  coords = c("X", "Y"), crs = 4326)

  res <- prepare_predictions(pts)

  expect_s3_class(res, "data.frame")
  expect_false(inherits(res, "sf"))
  expect_true(all(c("x", "y", "val") %in% names(res)))
  expect_equal(res$x, c(1, 2))

  # Case 2: Test prevention of duplicate x/y columns

  # If attributes already have an 'x', it should be replaced by the actual 'x'
  pts_dup <- st_as_sf(data.frame(X = 10, Y = 20, x = 999, mean = 0.5),
                      coords = c("X", "Y"), crs = 4326)
  res_dup <- prepare_predictions(pts_dup)

  expect_equal(sum(names(res_dup) == "x"), 1)
  expect_equal(res_dup$x, 10)  # Should be the coordinate, not the attribute 999
})

test_that("prepare_predictions handles non-standard geometry column names", {

  pts <- st_as_sf(data.frame(x = 1:2, y = 3:4, mean = 0.5), coords = c("x", "y"), crs = 4326)

  # Use the standard sf helper to rename the geometry column to 'geom'
  pts <- sf::st_set_geometry(pts, "geom")
  res <- prepare_predictions(pts)

  expect_s3_class(res, "data.frame")
  expect_null(res$geom) # confirms 'geom' was successfully identified and dropped
  expect_true("x" %in% names(res))
})

test_that("prepare_predictions handles model object classes", {

  # Test inlabru-style bru_prediction
  pts <- st_as_sf(data.frame(x = 1:2, y = 3:4, mean = c(0.1, 0.2)),
                  coords = c("x", "y"), crs = 4326)
  class(pts) <- c("bru_prediction", class(pts))

  res_bru <- prepare_predictions(pts)
  expect_true("mean" %in% names(res_bru))
  expect_s3_class(res_bru, "data.frame")

  # Test modISDM_predict internal class
  pred_sf  <- st_as_sf(data.frame(X = 1, Y = 2, mean = 0.5), coords = c("X", "Y"), crs = 4326)
  pred_list <- list(predictions = pred_sf)
  class(pred_list) <- "modISDM_predict"

  res_isdm <- prepare_predictions(pred_list)

  expect_s3_class(res_isdm, "data.frame")
  expect_true(all(c("x", "y") %in% names(res_isdm)))
  expect_equal(res_isdm$x, 1)
  expect_true("mean" %in% names(res_isdm))
})

test_that("prepare_predictions performs CRS-aware spatial filtering", {

  # Case 1: Point predictions in UTM (meters) with .vertex to trigger sf return
  pts <- st_as_sf(data.frame(x = c(500, 1500), y = c(500, 1500), .vertex = 1:2),
                  coords = c("x", "y"), crs = 32631)

  # Case 2: Base map in LongLat (degrees)
  poly <- st_as_sfc("POLYGON ((4.3 50.8, 4.5 50.8, 4.5 51.0, 4.3 51.0, 4.3 50.8))", crs = 4326)
  base_map <- st_sf(geometry = poly)

  # Test that it filters without error despite CRS mismatch
  expect_no_error(res <- prepare_predictions(pts, base_map = base_map))
  expect_s3_class(res, "sf")
})

test_that("prepare_predictions handles data.frame with uppercase X/Y", {
  df <- data.frame(X = 1:5, Y = 6:10, mean = 1:5)
  res <- prepare_predictions(df)

  expect_true(all(c("x", "y") %in% names(res)))
  expect_false("X" %in% names(res))
})

test_that("error handling for invalid inputs", {

  # Invalid base_map class
  expect_error(prepare_predictions(data.frame(x=1, y=1), base_map = "not_sf"),
               "'base_map' must be an sf object")

  # Invalid prediction data class
  expect_error(prepare_predictions(list(a = 1)),
               "Unsupported prediction object class")

  # Missing coordinates in data.frame
  expect_error(prepare_predictions(data.frame(val = 1:5)),
               "Input data lacks 'x' and 'y' coordinate columns")
})

if (require(fmesher, quietly = TRUE)) {

  test_that("prepare_predictions correctly handles and preserves bru_prediction structures", {
    set.seed(42)
    grid_df <- expand.grid(x = 0:20, y = 0:20)
    grid_df$mu <- (grid_df$x + grid_df$y) / 10 + rnorm(nrow(grid_df), 0, 0.1)
    grid_df$sd <- runif(nrow(grid_df), 0.1, 0.5)

    grid_r <- terra::rast(grid_df, crs = "epsg:4326")

    #--- Create Mesh and Vertices
    bnd  <- fm_nonconvex_hull(as.matrix(grid_df[, c("x", "y")]), convex = -0.10)
    mesh <- fm_mesh_2d(boundary = bnd, max.edge = c(5, 20), crs = "epsg:4326")
    vt   <- fm_vertices(mesh, format = "sf")

    #--- Extract outputs and Fill Buffer NAs with global values
    sampled_vals <- terra::extract(grid_r, vt)

    field_sim <- vt %>%
      dplyr::mutate(
        mean   = dplyr::coalesce(sampled_vals$mu, mean(grid_df$mu, na.rm = TRUE)),
        sd     = dplyr::coalesce(sampled_vals$sd, mean(grid_df$sd, na.rm = TRUE)),
        q0.025 = mean - 1.96 * sd,
        q0.975 = mean + 1.96 * sd,
        median = mean
      )

    # Mimic inlabru-like output
    class(field_sim) <- c("bru_prediction", "sf", "data.frame")

    # Assertions
    grid_pp3 <- prepare_predictions(field_sim)

    # Verify the class is maintained or transformed correctly
    expect_s3_class(grid_pp3, "bru_prediction")
    expect_s3_class(grid_pp3, "sf")

    # Verify crucial inlabru columns are present
    expect_true(".vertex" %in% names(grid_pp3))
    expect_true("mean" %in% names(grid_pp3))

    expect_false(any(is.na(grid_pp3$mean)))

    # Check that geometry is intact
    expect_true(inherits(sf::st_geometry(grid_pp3), "sfc_POINT"))
  })
}
