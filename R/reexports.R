# reexports.R
#
# Generics from other packages that stareal has methods for, exported so
# that they work without attaching those packages.
#
# ------------------------------------------------------------------

#' @importFrom loo loo
#' @export
loo::loo

#' @importFrom loo loo_compare
#' @export
loo::loo_compare
