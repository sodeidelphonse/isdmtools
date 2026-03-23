#' @title Replicate response data from the posterior predictive distribution.
#'
#' @description
#' The function generates replicated response data from the posterior samples of
#' parameters for use in model diagnostics. It currently supports Poisson,
#' Negative Binomial, and Quasi-Poisson families.
#'
#' @param mu_matrix A matrix of posterior samples of the expected count rate (mu).
#' @param family Character string defining the likelihood family (e.g., "poisson", "nbinomial").
#' @param dispersion Optional parameter (e.g., size 'k' for the negative binomial family).
#' @param ... Additional arguments.
#'
#' @return A matrix containing the simulated predictive samples
#' @importFrom stats rpois rnbinom
#' @export
#' @family PPC diagnostics
#'
#' @examples
#' \dontrun{
#' #--- Create a toy matrix of posterior samples for expected counts (mu)
#' # 5 sites (rows) and 3 posterior draws (columns)
#' mu_samples <- matrix(c(
#'   1.2, 1.5, 1.1,
#'   5.0, 4.8, 5.2,
#'   0.5, 0.2, 0.4,
#'   12.1, 11.5, 13,
#'   2.5, 2.2, 2.8
#' ), nrow = 5, ncol = 3, byrow = TRUE)
#'
#' #--- Generate Poisson replicates
#' # Returns a 5x3 matrix of random integer counts
#' sim_pois <- simulate_replicates(mu_samples, family = "poisson")
#' print(sim_pois)
#'
#' #--- Generate Negative Binomial replicates with high variance
#' # Requires a dispersion (size) parameter
#' sim_nb <- simulate_replicates(mu_samples,
#'   family = "nbinomial",
#'   dispersion = 0.8
#' )
#' print(sim_nb)
#'
#' #--- Integration with DHARMa (if installed)
#' if (requireNamespace("DHARMa", quietly = TRUE)) {
#'   # Define a mock observed response
#'   y_obs <- c(1, 5, 0, 12, 3)
#'
#'   # Create DHARMa object
#'   res_dharma <- DHARMa::createDHARMa(
#'     simulatedResponse = sim_pois,
#'     observedResponse = y_obs,
#'     fittedPredictedResponse = rowMeans(mu_samples),
#'     integerResponse = TRUE
#'   )
#'
#'   # Plot residuals
#'   par(mfrow = c(1, 2))
#'   plotQQunif(res_dharma)
#'   #    testDispersion(res_dharma)
#' }
#' }
simulate_replicates <- function(mu_matrix,
                                family = c("poisson", "nbinomial", "quasi-poisson"),
                                dispersion = NULL, ...) {
  family <- match.arg(family)

  if (family %in% c("nbinomial", "quasi-poisson") && is.null(dispersion)) {
    warning(paste("No dispersion parameter provided for", family, "- defaulting to 1."))
    dispersion <- 1
  }
  N_obs <- nrow(mu_matrix)
  N_sim <- ncol(mu_matrix)

  if (family == "nbinomial") {
    return(matrix(stats::rnbinom(N_obs * N_sim, mu = mu_matrix, size = dispersion), N_obs, N_sim))
  } else {
    return(matrix(stats::rpois(N_obs * N_sim, lambda = mu_matrix), N_obs, N_sim))
  }
}


#' @title Compute Posterior Predictive Check (PPC) Statistics
#'
#' @description This function performs a Bayesian Posterior Predictive Check for count
#' data models. It compares the observed Pearson Chi-squared statistic against a
#' distribution of statistics calculated from data replicated from the posterior
#' predictive distribution.
#'
#' @param samples_matrix A matrix of posterior samples of the expected count rate (mu).
#'   Rows correspond to observations, and columns correspond to Monte Carlo samples.
#' @param data_response A numeric vector of the observed count data.
#' @param family Character string defining the likelihood family: "poisson",
#'   "nbinomial", or "quasi-poisson".
#' @param dispersion Optional numeric value for the dispersion parameter (e.g., size 'k'
#'   for Negative Binomial). Defaults to 1 if NULL for "nbinomial" or "quasi-poisson".
#'
#' @details
#' The function calculates the Pearson Chi-squared statistic for the observed data
#' \eqn{T_{obs}} and for each replicated dataset \eqn{T_{rep}} generated from
#' the posterior samples.
#'
#' The statistic is defined as:
#' \deqn{T(y, \mu) = \sum_{i=1}^{n} \frac{(y_i - \mu_i)^2}{Var(y_i | \mu_i)}}
#'
#' **Interpretation:**
#' \itemize{
#'   \item \strong{Bayesian p-value:} Represents the probability that the
#'   replicated data is "more extreme" than the observed data.
#'   \itemize{
#'     \item \eqn{p_{B} \approx 0.5}: Indicates an excellent model fit.
#'     \item \eqn{p_{B} < 0.05} or \eqn{p_{B} > 0.95}: Indicates a significant lack of fit.
#'   }
#'   \item \strong{Overdispersion:} If \eqn{T_{obs}} is consistently larger than
#'   the \eqn{T_{rep}} distribution (\eqn{p \approx 1}).
#'   \item \strong{Underdispersion:} If \eqn{T_{obs}} is consistently smaller
#'   than the \eqn{T_{rep}} distribution (\eqn{p \approx 0}).
#' }
#'
#' @return A list of class "ppc_stats" containing:
#' \item{T_obs}{The observed Pearson Chi-squared statistic.}
#' \item{T_rep}{A vector of Pearson Chi-squared statistics from replicated data.}
#' \item{p_value}{The Bayesian p-value.}
#' \item{family}{The likelihood family used.}
#'
#' @importFrom stats rlnorm
#' @export
#' @family PPC diagnostics
#'
#' @examples
#' \dontrun{
#' #--- Setup mock data: 10 observations, 100 posterior samples
#' set.seed(123)
#' n_obs <- 10
#' n_sim <- 100
#'
#' # Expected counts (mu) from a hypothetical model
#' mu_samples <- matrix(rlnorm(n_obs * n_sim, meanlog = 1),
#'   nrow = n_obs, ncol = n_sim
#' )
#'
#' # Observed counts (y) - let's assume the model is a good fit
#' y_obs <- rpois(n_obs, rowMeans(mu_samples))
#'
#' #--- Compute PPC statistics for a Poisson family
#' ppc_results <- compute_ppc_stats(
#'   samples_matrix = mu_samples,
#'   data_response = y_obs,
#'   family = "poisson"
#' )
#'
#' #--- View the summary (triggers the `print.ppc_stats` method)
#' print(ppc_results)
#'
#' #--- Visualize the fit (triggers the `plot.ppc_stats` method)
#' plot(ppc_results)
#'
#' #--- Example with Negative Binomial and DHARMa integration
#' if (requireNamespace("DHARMa", quietly = TRUE)) {
#'   # Generate replicates of response for DHARMa
#'   nb_sims <- simulate_replicates(mu_samples, family = "nbinomial", dispersion = 2)
#'
#'   res_dharma <- createDHARMa(
#'     simulatedResponse = nb_sims,
#'     observedResponse = y_obs,
#'     fittedPredictedResponse = rowMeans(mu_samples)
#'   )
#'   plot(res_dharma)
#' }
#' }
compute_ppc_stats <- function(samples_matrix, data_response,
                              family = c("poisson", "nbinomial", "quasi-poisson"),
                              dispersion = NULL) {
  family <- match.arg(family)
  N_obs <- nrow(samples_matrix)
  N_sim <- ncol(samples_matrix)

  if (N_obs != length(data_response)) {
    stop("The number of observations (rows) in 'samples_matrix' must match 'data_response' length.")
  }

  # Ensure dispersion is handled for the internal engines
  if (family %in% c("nbinomial", "quasi-poisson") && is.null(dispersion)) {
    warning(paste("No dispersion parameter provided for", family, "- defaulting to 1."))
    dispersion <- 1
  } else if (is.null(dispersion)) {
    dispersion <- 1 # default for Poisson (unused)
  }

  #--- Observed Statistic T_obs
  mu_E <- rowMeans(samples_matrix)
  Var_E <- get_variance(mu_E, family, dispersion)
  T_obs <- sum((data_response - mu_E)^2 / Var_E)

  #--- Replicated Statistics T_rep
  Y_rep_matrix <- simulate_replicates(samples_matrix, family, dispersion)

  # Calculate T_rep values across the simulation samples
  T_rep <- sapply(1:N_sim, function(i) {
    mu_i <- samples_matrix[, i]
    Var_i <- get_variance(mu_i, family, dispersion)
    sum((Y_rep_matrix[, i] - mu_i)^2 / Var_i)
  })

  #--- Return Results
  res <- list(
    T_obs   = T_obs,
    T_rep   = T_rep,
    p_value = mean(T_rep > T_obs),
    family  = family
  )
  class(res) <- "ppc_stats"

  return(res)
}


#' @title Posterior Predictive Check (PPC) Statistics
#'
#' @description
#' \itemize{
#' \item \code{plot}: Visualizes the distribution of replicated T-statistics against
#' the observed T-statistic using ggplot2.
#' \item \code{print}: Display the Bayesian p-value and the related status of the
#' model fit to the data.
#' }
#'
#' @param x An object of class "ppc_stats".
#' @param ... Further arguments.
#'
#' @return
#' \itemize{
#'   \item \code{print}: Invisibly returns the original object.
#'   \item \code{plot}: Returns a \code{ggplot2} object.
#' }
#'
#' @name ppc_stats-methods
#' @rdname ppc_stats-methods
#' @export
#' @family PPC diagnostics
print.ppc_stats <- function(x, ...) {
  cat("\n--- Bayesian Posterior Predictive Check ---\n")
  cat("Likelihood Family:", x$family, "\n")
  cat("Observed T-stat:  ", round(x$T_obs, 3), "\n")
  cat("Mean Rep T-stat:  ", round(mean(x$T_rep), 3), "\n")
  cat("Bayesian p-value: ", round(x$p_value, 4), "\n")
  cat("-------------------------------------------\n")

  # Automated Interpretation
  if (x$p_value > 0.95 || x$p_value < 0.05) {
    status <- "Potential Lack of Fit"
    diag <- if (x$p_value > 0.5) "Possible Overdispersion" else "Possible Underdispersion"
  } else {
    status <- "Adequate Fit"
    diag <- "Variance is well-captured by the model"
  }

  cat("Status:         ", status, "\n")
  cat("Diagnostic:     ", diag, "\n\n")
}


#' @rdname ppc_stats-methods
#' @export
plot.ppc_stats <- function(x, ...) {
  df <- data.frame(T_rep = x$T_rep)
  fam_cap <- paste0(toupper(substr(x$family, 1, 1)), substr(x$family, 2, nchar(x$family)))

  plot_title <- paste("Posterior Predictive Check:", fam_cap)
  p_label <- paste("Bayesian p-value =", round(x$p_value, 3))

  ggplot2::ggplot(df, ggplot2::aes(x = .data$T_rep)) +
    ggplot2::geom_density(fill = "skyblue", alpha = 0.4, color = "steelblue", linewidth = 0.8) +
    ggplot2::geom_vline(ggplot2::aes(xintercept = x$T_obs),
      color = "red", linetype = "dashed", linewidth = 1.2
    ) +
    ggplot2::annotate("text",
      x = Inf, y = Inf, label = p_label,
      hjust = 1.1, vjust = 2, fontface = "italic", size = 4
    ) +
    ggplot2::labs(
      title = plot_title,
      subtitle = "Red dashed line indicates observed T-statistic",
      x = "Pearson Chi-squared Statistic (T)",
      y = "Density"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold"),
      panel.grid.minor = ggplot2::element_blank()
    )
}
