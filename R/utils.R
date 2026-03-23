
#' Calculate Niche Overlap (Schoener's D)
#'
#' @description Computes the overlap between two distributions.
#' @param x Numeric vector (e.g., predictions or spatial fold values).
#' @param y Numeric vector (e.g., observations or background values).
#' @param n Numeric. Number of points for density estimation. Default 512.
#'
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
  if (!is.numeric(x) || !is.numeric(y)) {
    return(NA_real_)
  }
  if (length(stats::na.omit(x)) < 2 || length(stats::na.omit(y)) < 2) {
    return(NA_real_)
  }
  if (stats::sd(x, na.rm = TRUE) == 0 || stats::sd(y, na.rm = TRUE) == 0) {
    return(NA_real_)
  }

  # Estimate densities
  rng <- range(c(x, y), na.rm = TRUE)
  dx <- stats::density(x, from = rng[1], to = rng[2], n = n)$y
  dy <- stats::density(y, from = rng[1], to = rng[2], n = n)$y

  # Normalize to probability distributions
  px <- dx / sum(dx)
  py <- dy / sum(dy)

  return(1 - 0.5 * sum(abs(px - py)))
}


#' Calculate the Matern covariance function (INLA-type)
#'
#' @param r The Euclidean distance
#' @param sigma2 The marginal variance of the spatial process
#' @param range The spatial range of the process
#' @param nu The smoothing parameter of the process (default is 1)
#' @param nugget The nugget effect (non-spatial variance, default is 0)
#' @param ... additional arguments
#' @noRd
matern_cov <- function(r, sigma2, rho, nu = 1, nugget = 0, ...) {
  kappa <- sqrt(8 * nu) / rho
  spat_cov <- ifelse(r > 0,
                     sigma2 * (2^(1 - nu) / gamma(nu)) *
                       (kappa * r)^nu * besselK(kappa * r, nu),
                     sigma2)

  # Add the nugget effect only at distance zero
  res <- ifelse(r == 0, spat_cov + nugget, spat_cov)
  return(res)
}


#' Get the variance of a count random variable
#' @param mu Numeric vector defining the expected count rate (mu).
#' @param family Character string defining the likelihood family (e.g., "poisson", "nbinomial").
#' @param dispersion Optional parameter (e.g., size 'k' for Negative binomial family)
#' @noRd
get_variance <- function(mu, family, dispersion = 1) {
  switch(family,
         "poisson" = mu,
         "nbinomial" = mu + (mu^2 / dispersion),
         "quasi-poisson" = dispersion * mu,
         stop("Unsupported family"))
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
      call. = FALSE
    )
  }
}
