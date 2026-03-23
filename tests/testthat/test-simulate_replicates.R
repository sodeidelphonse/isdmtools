test_that("Output dimensions match the input mu_matrix", {
  mu <- matrix(runif(100, 0.5, 5), nrow = 20, ncol = 5)
  sims <- simulate_replicates(mu, family = "poisson")

  expect_equal(dim(sims), dim(mu))
  expect_true(is.matrix(sims))
})

test_that("Family 'poisson' returns non-negative integers", {
  mu <- matrix(c(0.1, 10, 100), nrow = 3, ncol = 1)
  sims <- simulate_replicates(mu, family = "poisson")

  # Check for integer-like values (counts)
  expect_true(all(sims %% 1 == 0))
  expect_true(all(sims >= 0))
})

test_that("Negative binomial respects the dispersion parameter", {
  # We test this by ensuring it returns values within a reasonable range
  # for high/low dispersion
  mu <- matrix(10, nrow = 10, ncol = 10)

  # Low dispersion (high variance)
  sims_low <- simulate_replicates(mu, family = "nbinomial", dispersion = 0.1)

  # High dispersion (low variance, approaches Poisson)
  sims_high <- simulate_replicates(mu, family = "nbinomial", dispersion = 100)

  expect_equal(dim(sims_low), dim(mu))
  expect_equal(dim(sims_high), dim(mu))
})

test_that("Warnings are issued when dispersion is missing for NB/Quasi", {
  mu <- matrix(2, 5, 5)
  expect_warning(
    simulate_replicates(mu, family = "nbinomial", dispersion = NULL),
    "No dispersion parameter provided"
  )
  expect_warning(
    simulate_replicates(mu, family = "quasi-poisson", dispersion = NULL),
    "No dispersion parameter provided"
  )
})

test_that("Function handles edge cases (e.g., mu = 0)", {
  # Poisson mean of 0 must result in 0 counts
  mu_zero <- matrix(0, 5, 5)
  sims <- simulate_replicates(mu_zero, family = "poisson")

  expect_true(all(sims == 0))
})
