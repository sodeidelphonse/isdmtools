
#' Calculate Niche Overlap (Schoener's D)
#'
#' @description Computes the overlap between two distributions.
#' @param x Numeric vector (e.g., predictions or spatial fold values).
#' @param y Numeric vector (e.g., observations or background values).
#'
#' @param n Numeric. Number of points for density estimation. Default 512.
#' @return Numeric value between 0 and 1.
#' @export
#'
#' @examples
#' \dontrun{
#' library(isdmtools)
#'
#' # Create two identical distributions (should have high overlap)
#' v1 <- rnorm(1000, mean = 10, sd = 1)
#' v2 <- rnorm(1000, mean = 10, sd = 1)
#' overlap_high <- calc_niche_overlap(v1, v2)
#'
#' # Create two divergent distributions (should have low overlap)
#' v3 <- rnorm(1000, mean = 20, sd = 1)
#' overlap_low <- calc_niche_overlap(v1, v3)
#' }
#'
#' @references
#' Schoener TW. The Anolis Lizards of Bimini: Resource Partitioning in a Complex Fauna. _Ecology_(1968) 49:704–726. \doi{10.2307/1935534}.
#'
calc_niche_overlap <- function(x, y, n = 512) {
  if (!is.numeric(x) || !is.numeric(y)) return(NA_real_)
  if (length(stats::na.omit(x)) < 2 || length(stats::na.omit(y)) < 2) return(NA_real_)
  if (stats::sd(x, na.rm = TRUE) == 0 || stats::sd(y, na.rm = TRUE) == 0) return(NA_real_)

  # Estimate densities
  rng <- range(c(x, y), na.rm = TRUE)
  dx <- stats::density(x, from = rng[1], to = rng[2], n = n)$y
  dy <- stats::density(y, from = rng[1], to = rng[2], n = n)$y

  # Normalize to probability distributions
  px <- dx / sum(dx)
  py <- dy / sum(dy)

  return(1 - 0.5 * sum(abs(px - py)))
}


#' Solve the practical range of a Matérn covariance function
#'
#' @description
#' Harmonizes spatial range parameters from different R packages (\code{INLA},
#' \code{geoR}, \code{spatstat}) into a standardized "Practical Range." This is
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
#' @param sigma_sq Numeric. The partial sill (or marginal variance). Defaults to 1 for correlation focus.
#' @param thresh Numeric. The target correlation threshold. Defaults to 0.1 (10%).
#' @param package Character. One of \code{"inla"}, \code{"geor"}, or \code{"spatstat"}.
#'
#' @return A numeric value representing the practical range in the same
#' geographic units as the input model parameter.
#' @export
#' @importFrom stats uniroot
#'
#' @examples
#' # Estimated phi = 10 km with exponential covariance in `geoR`
#' solve_practical_range(param_val = 10, nu = 0.5, thresh = 0.1, package = "geor")
#'
#' # Estimated alpha = 10 km with Matérn covariance in `spatstat`
#' solve_practical_range(param_val = 13.10, nu = 1.5, thresh = 0.1, package = "spatstat")
#'
#' @references
#' \itemize{
#' \item Baddeley A, Rubak E, Turner R. Spatial point patterns: Methodology and applications with R. Boca Raton, FL: CHAPMAN & HALL CRC. (2015).
#' \item Diggle PJ, Ribeiro PJ. Model-based Geostatistics. 1st ed. New York, NY: Springer. (2007). \doi{10.1007/978-0-387-48536-2}
#' \item Lindgren F, Rue H, Lindström J. An explicit link between Gaussian fields and Gaussian Markov random fields: the stochastic partial differential equation approach.
#' _Journal of the Royal Statistical Society: Series B (Statistical Methodology)_ (2011) 73:423–498. \doi{10.1111/j.1467-9868.2011.00777.x}
#'}
solve_practical_range <- function(param_val,
                                  nu,
                                  sigma_sq = 1,
                                  thresh = 0.1,
                                  package = c("inla", "geor", "spatstat")
                                  ) {
  package <- match.arg(package)

  # Corr(X * d) = thresh/sigma_sq
  multiplier <- switch(package,
                       "geor"     = 1 / param_val,
                       "spatstat" = sqrt(2 * nu) / param_val,
                       "inla"     = sqrt(8 * nu) / param_val
                       )

  matern_corr <- function(d, nu, X) {
    if (d == 0) return(1)
    res <- (X * d)^nu * besselK(X * d, nu) / (2^(nu - 1) * gamma(nu))
    return(res)
  }

  # Solve for d: MaternCorr(d) - (Target Corr) = 0 (Corr = Cov/sigma_sq)
  target_corr <- thresh / sigma_sq

  sol <- uniroot(function(d) matern_corr(d, nu, multiplier) - target_corr,
                 interval = c(1e-8, param_val * 20))

  return(sol$root)
}


# Color Palette
.get_isdm_palette <- function(n) {
  mako_hex <- c(
    "#30123B", "#4454C4", "#4490FE", "#1EC4EB", "#2CF2B4",
    "#94FB68", "#F2E230", "#FEC029", "#F76E11", "#BB2001"
  )
  if (n <= length(mako_hex)) {
    return(mako_hex[seq_len(n)])
  } else { # Fallback for very high k
    return(grDevices::colorRampPalette(mako_hex)(n))
  }
}

# Check a package Namespace
.check_suggests <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(sprintf("Package '%s' is required for this function or option. Please install it.", pkg),
         call. = FALSE)
  }
}
