# --- Package documentation ---
#' @title A Toolkit for Integrated Species Distribution Models in R
#'
#' @description
#' `isdmtools` provides a set of tools for preparing, analyzing and visualizing spatial data for integrated species distribution models (ISDMs).
#' It is designed to help users prepare multisource spatial point datasets for block cross-validation, with a focus on Bayesian inference.
#' It analyses the habitat suitability from joint model predictions and maps the results. The package also provides a holistic view of model performance
#' by computing a comprehensive evaluation metrics, including ROC-based and continuous-outcome weighted composite scores.
#'
#' @details
#' For a complete guide on the workflow, see the package vignettes:
#' \itemize{
#'   \item \code{vignette("isdmtools", package = "isdmtools")}: Get started with the basics.
#'   \item \code{vignette("isdm-workflow", package = "isdmtools")}: Detailed evaluation workflow.
#' }
#'
#' @seealso
#' Useful links:
#' \itemize{
#'   \item Report bugs at \url{https://github.com/sodeidelphonse/isdmtools/issues}
#'   \item Official documentation: \url{https://sodeidelphonse.github.io/isdmtools/}
#' }
#'
#' @name isdmtools-package
#' @aliases isdmtools
#'
#' @importFrom dplyr %>%
#' @importFrom rlang .data
#' @keywords internal
"_PACKAGE"
