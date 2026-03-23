test_that("plot.std_matern_corr returns a ggplot object", {
  mock_res <- list(
    dist_vals = seq(0, 10, length.out = 50),
    cov = exp(-seq(0, 10, length.out = 50)), # simple declining
    pair_cor = exp(exp(-seq(0, 10, length.out = 50))),
    pair_cor_sc = seq(1, 0, length.out = 50), # Mock min-max
    sigma2 = 1,
    rho = 2,
    engine = "kppm"
  )
  class(mock_res) <- c("std_matern_corr", "list")

  # Test of the default plot (pair_cor)
  p1 <- plot(mock_res)
  expect_s3_class(p1, "ggplot")
  expect_equal(p1$labels$y, expression(g(r) == exp(C(r))))

  # Test of the plot with scaling (pair_cor_sc)
  p2 <- plot(mock_res, scaled = TRUE)
  expect_s3_class(p2, "ggplot")
  expect_match(p2$labels$y, "Min-Max")

  # Test of plot for covariance
  p3 <- plot(mock_res, type = "covariance")
  expect_s3_class(p3, "ggplot")
  expect_equal(p3$labels$y, "C(r)")
})

test_that("plot.std_matern_corr handles variofit objects correctly", {
  mock_vario <- list(
    dist_vals = seq(0, 5),
    cov = c(1, 0.8, 0.5, 0.3, 0.1, 0.05),
    sigma2 = 1,
    rho = 1,
    engine = "variofit"
  )
  class(mock_vario) <- c("std_matern_corr", "list")

  # Even if g(r) is requested, the covariance is plotted
  p_vario <- plot(mock_vario, type = "pair_cor")

  expect_s3_class(p_vario, "ggplot")
  expect_match(p_vario$labels$title, "Matern Covariance")
})


test_that("Plot y-axis start point is correct", {
  mock_res <- list(
    dist_vals = seq(0, 10),
    cov = 2 * exp(-seq(0, 10)),
    sigma2 = 2,
    engine = "variofit"
  )
  class(mock_res) <- c("std_matern_corr", "list")

  p <- plot(mock_res, type = "covariance", scaled = FALSE)

  # Check if the first point of the chart is sigma2
  expect_equal(p$data$y[1], 2)

  # If scaled = TRUE, the first point of the chart must be 1
  p_scaled <- plot(mock_res, type = "covariance", scaled = TRUE)
  expect_equal(p_scaled$data$y[1], 1)
})
