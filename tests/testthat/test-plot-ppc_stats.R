test_that("plot.ppc_stats returns a valid ggplot object", {
  res <- list(
    T_obs = 15,
    T_rep = rnorm(100, 15, 1),
    p_value = 0.5,
    family = "poisson"
  )
  class(res) <- "ppc_stats"

  p <- plot(res)
  expect_s3_class(p, "ggplot")

  # Check that the data inside the plot matches T_rep
  expect_equal(nrow(p$data), 100)
})

test_that("plot.ppc_stats generates correct title without tools dependency", {
  res <- list(T_obs = 10, T_rep = 8:12, p_value = 0.5, family = "poisson")
  class(res) <- "ppc_stats"

  p <- plot(res)

  # Check if "Poisson" is capitalized in the title
  expect_match(p$labels$title, "Poisson")
  expect_s3_class(p, "ggplot")
})
