#' *stareal*: spatial-temporal areal estimation from survey data
#'
#' The model: binomial counts on area-time cells, logit prevalence = fixed
#' effects + a separable Gaussian field (Leroux in space, AR1 in time) + a
#' cell-level nugget. The empty cells are integrated out exactly, so the MCMC
#' runs on the observed cells only; the full table is reconstructed
#' afterwards by conditional simulation (Matheron's rule).
#'
#' @keywords internal
"_PACKAGE"

## usethis namespace: start
#' @importFrom stats dbeta dbinom dnorm plogis qlogis quantile rbinom rnorm runif toeplitz
## usethis namespace: end
NULL
