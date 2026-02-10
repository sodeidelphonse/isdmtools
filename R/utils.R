
#' Calculate Niche Overlap (Schoener's D)
#'
#' @description Computes the overlap between two distributions.
#' @param x Numeric vector (e.g., predictions or fold values).
#' @param y Numeric vector (e.g., observations or background values).
#'
#' @param n Numeric. Number of points for density estimation. Default 512.
#' @return Numeric value between 0 and 1.
#' @noRd
#'
calc_niche_overlap <- function(x, y, n = 512) {
  if (!is.numeric(x) || !is.numeric(y)) return(NA_real_)
  if (length(na.omit(x)) < 2 || length(na.omit(y)) < 2) return(NA_real_)
  if (stats::sd(x, na.rm = TRUE) == 0 || stats::sd(y, na.rm = TRUE) == 0) return(NA_real_)

  # Define common range for both densities
  rng <- range(c(x, y), na.rm = TRUE)

  # Estimate densities
  dx <- stats::density(x, from = rng[1], to = rng[2], n = n)$y
  dy <- stats::density(y, from = rng[1], to = rng[2], n = n)$y

  # Normalize to probability distributions
  px <- dx / sum(dx)
  py <- dy / sum(dy)

  return(1 - 0.5 * sum(abs(px - py)))
}
