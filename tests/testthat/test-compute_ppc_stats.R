test_that("compute_ppc_stats returns correct structure", {
  mu <- matrix(rlnorm(50, 1, 0.2), 10, 5)
  y <- rpois(10, rowMeans(mu))

  res <- compute_ppc_stats(mu, y, family = "poisson")

  expect_s3_class(res, "ppc_stats")
  expect_length(res$T_rep, 5)
  expect_true(res$p_value >= 0 && res$p_value <= 1)
})

test_that("p-value is approximately 0.5 for perfectly matched data", {
  # Create a scenario where the model is 'true'
  set.seed(123)
  N <- 100
  S <- 500
  true_mu <- rlnorm(N, 1, 0.5)

  # Simulated response directly from the mean
  y_obs <- rpois(N, true_mu)

  # Samples matrix reflecting the true mean (ignoring parameter uncertainty for test)
  samples <- replicate(S, true_mu)

  res <- compute_ppc_stats(samples, y_obs, family = "poisson")

  # For a well-fitted model, p-value should not be at the extremes
  expect_gt(res$p_value, 0.05)
  expect_lt(res$p_value, 0.95)
})

test_that("Mismatched dimensions trigger an error", {
  # Create mismatched dummy data
  mu_bad <- matrix(runif(10), nrow = 5, ncol = 2) # 5 rows
  y_bad <- c(1, 2, 3) # 3 observations

  # Use a partial match of the actual error string
  expect_error(
    compute_ppc_stats(mu_bad, y_bad, family = "poisson"),
    regexp = "must match 'data_response' length"
  )
})

test_that("Negative Binomial handles dispersion correctly", {
  mu <- matrix(10, 10, 5)
  y <- rnbinom(10, mu = 10, size = 2)

  # If we provide the correct dispersion, it should pass
  expect_silent(compute_ppc_stats(mu, y, family = "nbinomial", dispersion = 2))
})


test_that("compute_ppc_stats accurately uses the variance engine", {
  # N=2, Sim=2
  mu <- matrix(c(2, 4, 2, 4), nrow = 2, ncol = 2)
  y_obs <- c(2, 4)

  # For Poisson, T_obs should be 0 because y_obs matches rowMeans(mu)
  res_pois <- compute_ppc_stats(mu, y_obs, family = "poisson")
  expect_equal(res_pois$T_obs, 0)

  # For NB with high dispersion (size=0.5), check if T_obs is calculated correctly
  # mu_E = c(2, 4). Var = mu + mu^2/0.5 -> c(2 + 4/0.5, 4 + 16/0.5) = c(10, 36)
  # But if y_obs matches mu_E, T_obs is still 0. Let's change y_obs to c(3, 5)
  y_new <- c(3, 5)

  # T_obs = (3-2)^2 / 10 + (5-4)^2 / 36 = 1/10 + 1/36 = 0.12777...
  res_nb <- compute_ppc_stats(mu, y_new, family = "nbinomial", dispersion = 0.5)
  expect_equal(res_nb$T_obs, (1 / 10 + 1 / 36))
})


test_that("compute_ppc_stats handles different likelihoods correctly", {
  # Setup: 2 observations, 100 samples
  set.seed(42)
  mu <- matrix(c(rep(2, 100), rep(10, 100)), nrow = 2, byrow = TRUE)
  y_obs <- c(2, 10) # Perfect fit: T_obs should be 0

  # Test Poisson
  res_pois <- compute_ppc_stats(mu, y_obs, family = "poisson")
  expect_equal(res_pois$T_obs, 0)
  expect_s3_class(res_pois, "ppc_stats")

  # Test Negative Binomial (T_obs should still be 0 if y == mu)
  res_nb <- compute_ppc_stats(mu, y_obs, family = "nbinomial", dispersion = 2)
  expect_equal(res_nb$T_obs, 0)

  # The p-value is the mean of (T_rep > T_obs). It must be between 0 and 1.
  expect_true(res_pois$p_value >= 0 && res_pois$p_value <= 1)
})

test_that("compute_ppc_stats throws error on dimension mismatch", {
  mu <- matrix(runif(10), nrow = 5) # 5 rows
  y_obs <- c(1, 2, 3) # 3 observations

  expect_error(
    compute_ppc_stats(mu, y_obs, family = "poisson"),
    "match 'data_response' length"
  )
})

test_that("simulate_replicates returns correct dimensions", {
  n_obs <- 5
  n_sim <- 10
  mu <- matrix(runif(n_obs * n_sim, 1, 5), nrow = n_obs)

  # Poisson
  sim_p <- simulate_replicates(mu, family = "poisson")
  expect_equal(dim(sim_p), c(n_obs, n_sim))
  expect_true(all(sim_p %% 1 == 0)) # must be integers

  # Negative binomial
  sim_nb <- simulate_replicates(mu, family = "nbinomial", dispersion = 1)
  expect_equal(dim(sim_nb), c(n_obs, n_sim))
})
