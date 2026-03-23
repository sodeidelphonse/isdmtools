test_that("print.ppc_stats outputs correct strings and interpretations", {
  # Create a mock ppc_stats object for an 'Adequate Fit'
  res_good <- list(
    T_obs = 15.0,
    T_rep = seq(10, 20, length.out = 100),
    p_value = 0.50,
    family = "poisson"
  )
  class(res_good) <- "ppc_stats"

  out_good <- capture.output(print(res_good))

  #--- Check for key phrases in the output
  expect_true(any(grepl("Likelihood Family: poisson", out_good)))
  expect_true(any(grepl("Status:\\s+Adequate Fit", out_good)))
  expect_true(any(grepl("Variance is well-captured", out_good)))
  expect_match(paste(out_good, collapse = "\n"), "0\\.5")

  #--- Test an extreme p-value for 'Lack of Fit' (Overdispersion)
  res_bad <- res_good
  res_bad$p_value <- 0.98

  out_bad <- capture.output(print(res_bad))

  expect_true(any(grepl("Status:\\s+Potential Lack of Fit", out_bad)))
  expect_true(any(grepl("Possible Overdispersion", out_bad)))

  #--- Test an extreme p-value for 'Lack of Fit' (Underdispersion)
  res_under <- res_good
  res_under$p_value <- 0.02

  out_under <- capture.output(print(res_under))

  expect_true(any(grepl("Possible Underdispersion", out_under)))
})
