
test_that("benin_shape is correctly integrated and spatially valid", {

  # Check if the object is accessible (lazy setup) with the right format
  expect_true(exists("benin_shape"))
  expect_s3_class(benin_shape, "sf")
  expect_s3_class(benin_shape, "data.frame")

  # Check the spatial integrity
  geom_type <- as.character(sf::st_geometry_type(benin_shape, by_geometry = FALSE))
  expect_match(geom_type, "POLYGON|MULTIPOLYGON")
  expect_equal(sf::st_crs(benin_shape)$epsg, 4326)

  # Benin is roughly between 0.7 and 3.9 longitude and 6.2 and 12.4 latitude
  bbox <- sf::st_bbox(benin_shape)
  expect_true(bbox$xmin < 1.5 && bbox$xmax > 2)
  expect_true(bbox$ymin < 7 && bbox$ymax > 11)
})

test_that("prepare_predictions works with internal data", {

  # Create a dummy point inside Benin
  test_point <- sf::st_sfc(sf::st_point(c(2.4, 9.3)), crs = 4326)

  # Check if the internal data can perform a spatial operation
  is_inside <- sf::st_intersects(test_point, benin_shape, sparse = FALSE)
  expect_true(is_inside[1,1])
})
