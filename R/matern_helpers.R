#' @title Standardise Matern covariance parameters and compute the
#' pairwise correlation from LGCP models
#'
#' @description Standardizes Matern covariance parameters across different modeling
#' engines (\code{INLA}, \code{spatstat}, \code{variofit}) and computes the
#' pairwise correlation function for log-Gaussian Cox process (LGCP) models.
#'
#' @param model The fitted model object from \code{INLA}, \code{spatstat}, or \code{geoR}
#' @param engine Character string specifying the modeling engine utilised to estimate
#' the covariance parameters: \code{"inla"}, \code{"kppm"}, or \code{"variofit"}.
#' @param r A numeric vector of distances at which to evaluate the covariance function.
#' @param nu The smoothing parameter of the process (it is \eqn{\kappa} in `geoR`).
#' @param ... Additional arguments to pass on to internal function.
#'
#' @details
#' The function standardizes Matern covariance parameters and derives the pairwise
#' correlation function across different modeling engines for LGCP models. The engines
#' currently supported include those implemented in \code{INLA/inlabru} (and related wrappers),
#' \code{spatstat}, and \code{geoR}.
#'
#' Note that the \code{variofit} class is included to derive the range and \code{sigma2}
#' parameters for comparison purposes; however, the exponential of the resulting covariance
#' has no interpretation in terms of pairwise correlation. The same applies to covariance
#' functions fitted to geostatistical data in \code{INLA/inlabru}, which cannot be
#' interpreted in terms of pairwise correlation.
#'
#' To facilitate visual comparison across different modeling engines, a scaled version
#' of the pairwise correlation function (\code{pair_cor_sc}) is provided. This version
#' applies a min-max normalization to the pairwise correlation, mapping its values
#' to a standardized range (typically \[0, 1\]). This ensures that the spatial decay
#' and interaction strength can be compared even when models have different
#' marginal variances or measurement scales.
#'
#' @return An object of class \code{std_matern_corr} (list) containing:
#' \itemize{
#'   \item \code{dist_vals}: The vector of distances (\eqn{r}) used for computation.
#'   \item \code{cov}: The computed Matern covariance: \eqn{C(r)}.
#'   \item \code{sigma2}: The standardized marginal variance (\eqn{\sigma^2}).
#'   \item \code{rho}: The standardized scale parameter (range).
#'   \item \code{nu}: The smoothness parameter.
#'   \item \code{engine}: The original modeling engine used.
#' }
#'
#' For \code{"inla"} or \code{"kppm"} modeling engines, two additional elements are included:
#' \itemize{
#'   \item \code{pair_cor}: The pairwise correlation function (PCF) for LGCP,
#'   calculated as \eqn{g(r) = \exp(C(r))}.
#'   \item \code{pair_cor_sc}: The standardized pairwise correlation function for
#'   cross-model comparison.
#' }
#'
#' @export
#' @family Matern covariance helpers
#'
#' @examples
#' \dontrun{
#' # Example with an INLA model
#' r_vec <- seq(0, 10, length.out = 100)
#' matern_results <- std_matern_corr(
#'   model = fit_inla,
#'   engine = "inla",
#'   r = r_vec,
#'   nu = 1
#' )
#'
#' # Access the pair correlation function for LGCP
#' plot(matern_results$dist_vals, matern_results$cor,
#'   type = "l",
#'   ylab = "PCF", xlab = "Distance"
#' )
#' }
#'
#' @references
#' \itemize{
#' \item Baddeley A, Rubak E, Turner R. Spatial point patterns: Methodology and applications with R. Boca Raton, FL: CHAPMAN & HALL CRC. (2015).
#' \item Diggle PJ, Ribeiro PJ. Model-based Geostatistics. 1st ed. New York, NY: Springer. (2007). \doi{10.1007/978-0-387-48536-2}
#' \item Lindgren F, Rue H, Lindström J. An explicit link between Gaussian fields and Gaussian Markov random fields: the stochastic partial differential equation approach.
#' _Journal of the Royal Statistical Society: Series B (Statistical Methodology)_ (2011) 73:423–498. \doi{10.1111/j.1467-9868.2011.00777.x}
#' }
std_matern_corr <- function(model, engine = c("kppm", "inla", "variofit"), r, nu = 1, ...) {
  engine <- match.arg(engine)

  if (engine == "kppm") {
    if (!inherits(model, "kppm")) {
      stop(sprintf(
        "'%s' must be a 'kppm' model when 'engine' is 'kppm'.",
        model
      ), call. = FALSE)
    }
    sigma2 <- model$par[1]
    rho <- model$par[2] * 2 # rho = 2*alpha
  } else if (engine == "inla") {
    inla_classes <- c("modISDM", "bru", "iinla", "inla")
    if (!inherits(model, inla_classes)) {
      stop(sprintf(
        "'%s' must be an 'R-INLA', 'bru' or 'modISDM' object when 'engine' is 'inla'.",
        model
      ), call. = FALSE)
    }
    sigma2 <- (model$summary.hyperpar$mean[2])^2
    rho <- model$summary.hyperpar$mean[1]
  } else if (engine == "variofit") {
    if (!inherits(model, c("variomodel", "variofit"))) {
      stop(sprintf(
        "'%s' must be a 'variofit' model when 'engine' is 'variofit'.",
        model
      ), call. = FALSE)
    }
    sigma2 <- model$cov.pars[1]
    rho <- model$cov.pars[2] * sqrt(8 * nu) # rho = sqrt(8*nu)*phi
  }

  # Pairwise correlation function from an LGCP covariance
  c_r <- matern_cov(r, sigma2, rho, nu = nu, ...)
  pair_cor <- exp(c_r)
  pair_cor_sc <- (pair_cor - min(pair_cor, na.rm = TRUE)) / diff(range(pair_cor, na.rm = TRUE))

  common_elements <- list(
    dist_vals = r,
    cov = c_r,
    sigma2 = sigma2,
    rho = rho,
    nu = nu,
    engine = engine
  )

  if (engine == "variofit") {
    out <- common_elements
  } else {
    out <- c(common_elements, list(pair_cor = pair_cor, pair_cor_sc = pair_cor_sc))
  }

  class(out) <- c("std_matern_corr", "list")
  return(out)
}


#' Plot standardized Matern correlation or covariance
#'
#' @description Visualization method for \code{std_matern_corr} objects using \code{ggplot2}.
#' It can display the Pairwise Correlation Function (PCF) for LGCP models or the
#' standard Matern covariance.
#'
#' @param x An object of class \code{std_matern_corr}.
#' @param type Character string. Either \code{"pair_cor"} (default for LGCP) or
#' \code{"covariance"}.
#' @param scaled Logical. If \code{TRUE}, uses the min-max scaled PCF (\code{pair_cor_sc})
#' or the standardized covariance (\eqn{C(r)/\sigma^2}). Default is \code{FALSE}.
#' @param ... Additional arguments (currently ignored).
#'
#' @return A \code{ggplot} object.
#' @family Matern covariance helpers
#' @export
#'
#' @examples
#' \dontrun{
#' # Assuming 'cov_lgcp' is an object returned by std_matern_corr()
#' plot(cov_lgcp, type = "pair_cor", scaled = TRUE)
#'
#' # To modify the plot
#' library(ggplot2)
#' p <- plot(cov_lgcp, type = "covariance")
#' p + theme_bw() + ggtitle("Spatial Decay")
#' }
plot.std_matern_corr <- function(x, type = c("pair_cor", "covariance"), scaled = FALSE, ...) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for this plotting method.")
  }

  type <- match.arg(type)

  if (x$engine == "variofit") type <- "covariance"
  if (type == "pair_cor") {
    y_vals <- if (scaled) x$pair_cor_sc else x$pair_cor
    y_label <- if (scaled) "g(r) [Min-Max scaled]" else expression(g(r) == exp(C(r)))
    title_suffix <- "Pairwise Correlation (g(r))"
  } else {
    y_vals <- if (scaled) x$cov / x$sigma2 else x$cov
    y_label <- if (scaled) expression(rho(r) == C(r) / sigma^2) else "C(r)"
    title_suffix <- "Matern Covariance"
  }

  df_plot <- data.frame(r = x$dist_vals, y = y_vals)

  p <- ggplot2::ggplot(df_plot, ggplot2::aes(x = .data$r, y = .data$y)) +
    ggplot2::geom_line(color = "midnightblue", linewidth = 1) +
    ggplot2::labs(
      title = paste(x$engine, "-", title_suffix),
      subtitle = if (scaled) "Normalized values" else "Absolute values",
      x = "Distance (r)",
      y = y_label
    ) +
    ggplot2::theme_minimal()

  return(p)
}


#' Solve the practical range of a Matérn covariance function
#'
#' @description
#' Harmonizes spatial range parameters from different R packages (\code{INLA},
#' \code{geoR}, \code{spatstat}) into a standardized "Practical Range". This is
#' the distance at which the spatial correlation drops to a specific threshold
#' (default is 0.10).
#'
#' @details
#' Different packages use different parameterisations for the Matérn covariance:
#' \itemize{
#'   \item \strong{INLA/inlabru:} Estimates a value close to the INLA range parameter
#'   (where correlation is ~ 0.139). If \code{thresh = 0.139}, the input \code{param_val}
#'   is returned almost as is. If a 5% threshold (\code{thresh = 0.05}) is desired, the
#'   function adjusts the INLA range accordingly.
#'   \item \strong{geoR:} Uses a scale parameter \eqn{\phi}. The practical range
#'   is solved numerically based on \eqn{\phi} and the smoothness \eqn{\nu}.
#'   \item \strong{spatstat:} Uses a scale parameter \eqn{\alpha}. The function
#'   aligns this with the INLA-style practical range.
#' }
#'
#' This harmonization ensures that the \code{rho} value used in \code{isdmtools}
#' diagnostic functions is consistent, regardless of the modeling engine used
#' for the exploratory analysis.
#'
#' @param param_val Numeric. The parameter value from the model (\eqn{\rho} for INLA,
#' \eqn{\phi} for geoR, or \eqn{\alpha} for spatstat). It must be positive.
#' @param nu Numeric. The smoothness parameter. It must be positive.
#' For 2D SPDE models in INLA (where alpha = 2), the default is \code{nu = 1}.
#' For an exponential covariance, use \code{nu = 0.5}.
#' @param thresh Numeric. The target correlation threshold. Defaults to 0.1 (10%).
#' @param engine Character. One of \code{"inla"}, \code{"geor"}, or \code{"spatstat"}.
#'
#' @return A numeric value representing the practical range in the same
#' geographic units as the input model parameter.
#'
#' @importFrom stats uniroot
#' @export
#' @family Matern covariance helpers
#'
#' @examples
#' # Estimated phi = 10 km with exponential covariance in `geoR`
#' solve_practical_range(param_val = 10, nu = 0.5, thresh = 0.1, engine = "geor")
#'
#' # Estimated alpha = 13.10 km with Matérn covariance in `spatstat`
#' solve_practical_range(param_val = 13.10, nu = 1.5, thresh = 0.1, engine = "spatstat")
#'
#' @references
#' \itemize{
#' \item Baddeley A, Rubak E, Turner R. Spatial point patterns: Methodology and applications with R. Boca Raton, FL: CHAPMAN & HALL CRC. (2015).
#' \item Diggle PJ, Ribeiro PJ. Model-based Geostatistics. 1st ed. New York, NY: Springer. (2007). \doi{10.1007/978-0-387-48536-2}
#' \item Lindgren F, Rue H, Lindström J. An explicit link between Gaussian fields and Gaussian Markov random fields: the stochastic partial differential equation approach.
#' _Journal of the Royal Statistical Society: Series B (Statistical Methodology)_ (2011) 73:423–498. \doi{10.1111/j.1467-9868.2011.00777.x}
#' }
solve_practical_range <- function(param_val,
                                  nu,
                                  thresh = 0.1,
                                  engine = c("inla", "geor", "spatstat")) {
  engine <- match.arg(engine)

  # Corr(X * d) = thresh
  multiplier <- switch(engine,
    "geor"     = 1 / param_val,
    "spatstat" = sqrt(2 * nu) / param_val,
    "inla"     = sqrt(8 * nu) / param_val
  )

  matern_corr <- function(d, nu, X) {
    if (d == 0) {
      return(1)
    }
    res <- (X * d)^nu * besselK(X * d, nu) / (2^(nu - 1) * gamma(nu))
    return(res)
  }

  # Solve for d: MaternCorr(d) - (Target Corr) = 0
  target_corr <- thresh

  sol <- uniroot(function(d) matern_corr(d, nu, multiplier) - target_corr,
    interval = c(1e-8, param_val * 50)
  )

  return(sol$root)
}
