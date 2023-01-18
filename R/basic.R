#' Basic quantile functions
#'
#' Basic quantile functions available for modification
#' `sqf_exp()`: Basic QF of exponential distribution
#' `sqf_unif()`: Basic QF of uniform distribution
#' `sqf_norm()`: Basic QF of normal distribution `qnorm()`.
#' `sqf_cauchy()`: Basic QF of Cauchy distribution`.
#' `sqf_halftriang()`: Basic QF of half-triangular distribution.
#' `sqf_halfcosine()`: Basic QF of half-cosine distribution.
#' `sqf_sech()`: Basic QF of hyperbolic secant distribution.
#' @param u depth
#' @param ... not used
#'
#' @return modified function
#' @rdname basicqf
#' @export
#' @examples
#' sqf_exp(0.6)
#' qexp(0.6, 1)
sqf_exp <- function(u, ...){
  -log(1-u)
}
#' @rdname basicqf
#' @export
sqf_unif <- function(u, ...){
  u
}
#' @rdname basicqf
#' @export
#' @importFrom stats qnorm
sqf_norm <- function(u, ...){
  qnorm(u)
}
#' @rdname basicqf
#' @export
sqf_cauchy <- function(u, ...){
  tan(pi*(u-0.5))
}
#' @rdname basicqf
#' @export
sqf_halftriang <- function(u, ...){
  -sqrt(1-u)
}
#' @rdname basicqf
#' @export
sqf_halfcosine <- function(u, ...){
  asin(u)
}
#' @rdname basicqf
#' @export
sqf_sech <- function(u,...){
  2/pi*log(tan(pi/2*u))
}
