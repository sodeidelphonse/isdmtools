test_that("ISDMmetrics data structure handles NA and weights correctly", {

  # Setup mock data
  n_pres <- 300
  n_count <- 50

  auc_p <- 0.80
  auc_c <- 0.60
  rmse_p <- NA_real_ # Presence-only has no Error metric
  rmse_c <- 0.40

  # Composite Error ignores NA (just takes the Count value)
  auc_tot <- ((auc_p * n_pres) + (auc_c * n_count)) / (n_pres + n_count)

  m <- list(
    AUC_Presence  = auc_p,
    AUC_Count     = auc_c,
    RMSE_Presence = rmse_p,
    RMSE_Count    = rmse_c,
    TOT_ROC_SCORE = auc_tot,
    TOT_ERROR_SCORE = rmse_c
  )
  class(m) <- c("ISDMmetrics", "list")

  #--- Tests (as.data.frame) ---
  df <- as.data.frame(m)

  # Verify if the row actually exists before checking the value
  expect_true("RMSE_Presence" %in% df$Full_Name)

  val_rmse_p <- df$Value[df$Full_Name == "RMSE_Presence"]
  expect_true(is.na(val_rmse_p))

  # A. Check column existence and types
  expect_s3_class(df, "data.frame")
  expect_type(df$Value, "double")

  # B. Check NA handling for presence-only error metric
  val_rmse_p <- df$Value[df$Full_Name == "RMSE_Presence"]
  expect_true(is.na(val_rmse_p))

  # C. Check weighted calculation accuracy
  val_tot_roc <- df$Value[df$Full_Name == "TOT_ROC_SCORE"]
  expect_equal(val_tot_roc, 0.7714, tolerance = 1e-4)

  # D. Check source labeling logic

  # TOT_ names should be mapped to "Global"
  expect_equal(df$Source[df$Full_Name == "TOT_ERROR_SCORE"], "Global")

  # Individual names should be mapped to the Dataset name
  expect_equal(df$Source[df$Full_Name == "AUC_Presence"], "Presence")
})
