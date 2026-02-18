test_that("fill_na_near handles input validation and conversions", {

  r <- terra::rast(matrix(c(1, 2, NA, 4), nrow = 2))

  # Test input arguments
  expect_error(fill_na_near(data.frame(a = 1)), "must be a 'SpatRaster'")
  expect_error(fill_na_near(r, start.window = 2), "must be an odd number")
  expect_error(fill_na_near(r, w = 3), "argument cannot be set explicitly")
})

test_that("fill_na_near logic and boundary masking works", {

  # Create small raster with a central NA
  r <- terra::rast(ncols = 3, nrows = 3, vals = 1)
  r[5] <- NA

  # Test basic filling
  filled <- fill_na_near(r, start.window = 1)
  expect_false(any(is.na(terra::values(filled))))
  expect_equal(as.numeric(terra::values(filled)[5]), 1)

  expect_error(fill_na_near(r, boundary = "invalid"))

  # Test with a valid boundary (sf object)
  poly <- sf::st_as_sfc("POLYGON((0 0, 0 1, 1 1, 1 0, 0 0))")
  poly_df <- sf::st_sf(geometry = poly)

  terra::crs(r) <- "EPSG:4326"
  sf::st_crs(poly_df) <- "EPSG:4326"

  filled_bnd <- fill_na_near(r, boundary = poly_df)
  expect_s4_class(filled_bnd, "SpatRaster")
})
