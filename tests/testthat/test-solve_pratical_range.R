
test_that("solve_practical_range harmonizes different package conventions correctly", {

  # Identity: INLA's internal range is defined at approx 0.139.
  # Using this threshold should return the input value nearly exactly.
  inla_val    <- 100
  res_inla_id <- solve_practical_range(param_val = inla_val, nu = 1.5,
                                       thresh = 0.139, package = "inla")
  expect_equal(res_inla_id, inla_val, tolerance = 0.01)

  # Dropping from 13.9% to 10% must result in a larger distance.
  res_10 <- solve_practical_range(param_val = inla_val, nu = 1.5,
                                  thresh = 0.10, package = "inla")
  expect_gt(res_10, inla_val)

  # A 5% practical range should always be larger than a 10% practical range
  res_05 <- solve_practical_range(param_val = inla_val, nu = 1.5,
                                  thresh = 0.05, package = "inla")
  expect_gt(res_05, res_10)

  print(res_10)
  print(res_05)

  # geoR: For Exponential correlation (nu = 0.5), corr(d) = exp(-d/phi).
  # Solve exp(-d/50) = 0.05 => -d/50 = ln(0.05) => d = -50 * -2.9957
  phi_val  <- 50
  res_geor <- solve_practical_range(param_val = phi_val, nu = 0.5,
                                    package = "geor", thresh = 0.05)
  expect_equal(res_geor, 149.78, tolerance = 0.01)  # ~ 3 * phi


  # Test Consistency
  # A 5% threshold must always yield a larger range than a 10% regardless of package and nu
  val_10 <- solve_practical_range(100, nu = 1.5, thresh = 0.1, package = "spatstat")
  val_05 <- solve_practical_range(100, nu = 1.5, thresh = 0.05, package = "spatstat")
  expect_gt(val_05, val_10)
})


test_that("solve_practical_range aligns with INLA default nu = 1", {

  # Using the INLA internal threshold (~0.139) should return the input.
  inla_val <- 100
  res_inla <- solve_practical_range(param_val = inla_val, nu = 1,
                                    thresh = 0.139, package = "inla")

  res_10 <- solve_practical_range(param_val = inla_val, nu = 1,
                                  thresh = 0.1, package = "inla")

  # Practical range (10%) must be slightly larger than the INLA range (13.9%)
  expect_gt(res_10, inla_val)

  # Validation of the 10% result (approx 1.13 x the INLA range for nu=1)
  expect_equal(res_10, 113.6, tolerance = 0.1)
})
