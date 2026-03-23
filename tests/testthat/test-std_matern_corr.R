test_that("std_matern_corr handles geoR (variofit) mock objects", {
  # Mock of variofit object
  mock_geor <- list(cov.pars = c(1.5, 0.5)) # sigma2 = 1.5, phi = 0.5
  class(mock_geor) <- "variofit"

  nu_test <- 1
  r_test <- seq(0, 1, length.out = 10)

  # Expected value : rho = phi * sqrt(8*nu)
  expected_rho <- 0.5 * sqrt(8 * nu_test)

  res <- std_matern_corr(mock_geor, engine = "variofit", r = r_test, nu = nu_test)

  expect_equal(res$sigma2, 1.5)
  expect_equal(res$rho, expected_rho)
  expect_null(res$pair_cor) # no PCF for variofit
})

test_that("std_matern_corr handles spatstat (kppm) mock objects", {
  # Mock of kppm object from spatstat
  mock_kppm <- list(par = c(0.8, 2.0)) # sigma2 = 0.8, alpha = 2.0
  class(mock_kppm) <- "kppm"

  # Expectation : rho = 2 * alpha
  expected_rho <- 2.0 * 2

  res <- std_matern_corr(mock_kppm, engine = "kppm", r = seq(0, 1), nu = 1)

  expect_equal(res$sigma2, 0.8)
  expect_equal(res$rho, expected_rho)
  expect_false(is.null(res$pair_cor))
})

test_that("std_matern_corr handles INLA mock objects", {
  mock_inla <- list(
    summary.hyperpar = data.frame(
      mean = c(3.5, 1.2), # range = 3.5, precision/stdev = 1.2
      row.names = c("Range for spatial", "Stdev for spatial")
    )
  )
  class(mock_inla) <- "inla"

  # Expected value : sigma2 = mean[2]^2, rho = mean[1]
  expected_sigma2 <- 1.2^2
  expected_rho <- 3.5

  res <- std_matern_corr(mock_inla, engine = "inla", r = seq(0, 1))

  expect_equal(res$sigma2, expected_sigma2)
  expect_equal(res$rho, expected_rho)
})

test_that("std_matern_corr throws error for unsupported classes", {
  unknown_obj <- list(a = 1)
  expect_error(
    std_matern_corr(unknown_obj, engine = "inla", r = 1),
    "must be an 'R-INLA', 'bru' or 'modISDM' object"
  )
})


test_that("Covariance at distance r=0 is exactly sigma2", {
  mock_inla <- list(
    summary.hyperpar = data.frame(mean = c(3, 1.5), row.names = c("range", "stdev"))
  )
  class(mock_inla) <- "inla"

  # sigma2 = 1.5^2 = 2.25
  res <- std_matern_corr(mock_inla, engine = "inla", r = c(0, 5, 10))

  # The first value(r=0) must be sigma2
  expect_equal(res$cov[1], res$sigma2)

  # The standardised Matern covariance rho(0) must be 1
  expect_equal(res$cov[1] / res$sigma2, 1)
})

test_that("Pairwise correlation at r=0 follows exp(sigma2)", {
  mock_kppm <- list(par = c(1.2, 0.5)) # sigma2 = 1.2
  class(mock_kppm) <- "kppm"

  res <- std_matern_corr(mock_kppm, engine = "kppm", r = 0)

  # g(0) = exp(C(0)) = exp(sigma2)
  expect_equal(res$pair_cor, exp(1.2))
})


test_that("matern_cov returns exact values at boundaries", {
  s2 <- 1.5
  rho <- 2.0
  nu_val <- 1

  # Test at exactly zero (must return sigma2)
  expect_equal(matern_cov(0, s2, rho, nu_val), s2)

  # Test for very small distances (should approach sigma2)
  expect_equal(matern_cov(1e-10, s2, rho, nu_val), s2, tolerance = 1e-6)

  # Test for very large distances (should approach zero)
  expect_lt(matern_cov(1000, s2, rho, nu_val), 1e-5)
})

test_that("matern_cov handles vector inputs correctly", {
  r_vec <- c(0, 1, 2, 5)
  s2 <- 1
  rho <- 1

  res <- matern_cov(r_vec, s2, rho)

  expect_length(res, length(r_vec))
  expect_equal(res[1], s2)

  # Covariance must be strictly decreasing
  expect_true(all(diff(res) < 0))
})


test_that("matern_cov handles nugget effect correctly", {
  s2 <- 2.0 # Partial sill
  rho <- 1.5 # Range
  nug <- 0.5 # Nugget

  # Total variance at r=0 must be sigma2 + nugget
  expect_equal(matern_cov(0, s2, rho, nugget = nug), s2 + nug)

  # At r > 0, the nugget should not be added (only spatial cov remains)
  # We compare with a version where nugget = 0
  r_dist <- 0.1
  expect_equal(
    matern_cov(r_dist, s2, rho, nugget = nug),
    matern_cov(r_dist, s2, rho, nugget = 0)
  )

  # If nugget is 0, it should fall back to standard behavior
  expect_equal(matern_cov(0, s2, rho, nugget = 0), s2)
})

test_that("matern_cov remains strictly decreasing with nugget", {
  r_seq <- c(0, 0.001, 1)
  res <- matern_cov(r_seq, sigma2 = 1, rho = 1, nugget = 0.2)

  # There should be a sharp drop from r=0 to r>0 due to the nugget
  expect_true(res[1] > res[2])
  expect_true(res[2] > res[3])
})
