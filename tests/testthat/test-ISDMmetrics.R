test_that("ISDMmetrics data structure handles NA and weights correctly", {

  # Setup mock data
  n_pres <- 300
  n_count <- 50

  auc_p <- 0.80
  auc_c <- 0.70
  tss_p <- 0.70
  tss_c <- 0.60

  rmse_p <- NA_real_  # Presence-only has no Error metric
  rmse_c <- 0.40
  mae_p  <- NA_real_  # Presence-only has no Error metric
  mae_c  <- 0.60

  auc_comp <- ((auc_p * n_pres) + (auc_c * n_count)) / (n_pres + n_count)
  tss_comp <- ((tss_p * n_pres) + (tss_c * n_count)) / (n_pres + n_count)
  tot_roc <- (auc_comp + tss_comp)/2

  # Composite Error ignores NA (just takes the Count value)
  rmse_comp <- rmse_c
  mae_comp  <- mae_c
  tot_error <- (rmse_comp + mae_comp)/2

  m <- list(
    AUC_Presence  = auc_p,
    AUC_Count     = auc_c,
    TSS_Presence = tss_p,
    TSS_Count    = tss_c,
    RMSE_Presence = rmse_p,
    RMSE_Count    = rmse_c,
    MAE_Presence = mae_p,
    MAE_Count    = mae_c,
    AUC_Comp     = auc_comp,
    TSS_Comp     = tss_comp,
    RMSE_Comp    = rmse_comp,
    MAE_Comp     = mae_comp,
    TOT_ROC_SCORE = tot_roc,
    TOT_ERROR_SCORE = tot_error
  )

  class(m) <- c("ISDMmetrics", "list")

  #--- Tests (as.data.frame) ---
  df <- as.data.frame(m)

  # Verify if the row actually exists before checking the value
  expect_true("RMSE_Presence" %in% df$Full_Name)
  expect_true("MAE_Presence" %in% df$Full_Name)

  val_rmse_p <- df$Value[df$Full_Name == "RMSE_Presence"]
  expect_true(is.na(val_rmse_p))

  val_mae_p <- df$Value[df$Full_Name == "MAE_Presence"]
  expect_true(is.na(val_mae_p))

  # Check column existence and types
  expect_s3_class(df, "data.frame")
  expect_type(df$Value, "double")

  # Check NA handling for presence-only error metric
  val_rmse_p <- df$Value[df$Full_Name == "RMSE_Presence"]
  expect_true(is.na(val_rmse_p))

  val_mae_p <- df$Value[df$Full_Name == "MAE_Presence"]
  expect_true(is.na(val_mae_p))

  # Check weighted calculation accuracy
  val_auc_comp <- df$Value[df$Full_Name == "AUC_Comp"]
  expect_equal(val_auc_comp, 0.7857, tolerance = 1e-4)

  val_tss_comp <- df$Value[df$Full_Name == "TSS_Comp"]
  expect_equal(val_tss_comp, 0.6857, tolerance = 1e-4)

  val_rmse_comp <- df$Value[df$Full_Name == "RMSE_Comp"]
  expect_equal(val_rmse_comp, 0.4, tolerance = 1e-4)

  val_mae_comp <- df$Value[df$Full_Name == "MAE_Comp"]
  expect_equal(val_mae_comp, 0.6, tolerance = 1e-4)

  # Check overall calculation accuracy
  val_tot_roc <- df$Value[df$Full_Name == "TOT_ROC_SCORE"]
  expect_equal(val_tot_roc, 0.7357, tolerance = 1e-4)

  val_tot_error <- df$Value[df$Full_Name == "TOT_ERROR_SCORE"]
  expect_equal(val_tot_error, 0.5, tolerance = 1e-4)

  # Check data source labeling logic
  # "TOT_*" should be mapped to "Global" and "*_Comp" to "Weighted Composite"
  expect_equal(df$Source[df$Full_Name == "TOT_ERROR_SCORE"], "Global")
  expect_equal(df$Source[df$Full_Name == "TOT_ROC_SCORE"], "Global")

  expect_equal(df$Source[df$Full_Name == "AUC_Comp"], "Weighted Composite")
  expect_equal(df$Source[df$Full_Name == "TSS_Comp"], "Weighted Composite")

  # Individual names should be mapped to the Dataset name
  expect_equal(df$Source[df$Full_Name == "AUC_Presence"], "Presence")
  expect_equal(df$Source[df$Full_Name == "TSS_Presence"], "Presence")
})
