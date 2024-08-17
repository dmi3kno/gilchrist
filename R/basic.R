#' Basic quantile functions
#'
#' Basic quantile functions available for modification
#' `s_exp()`: Basic QF of exponential distribution
#' `s_unif()`: Basic QF of uniform distribution
#' `s_norm()`: Basic QF of normal distribution `qnorm()`.
#' `s_cauchy()`: Basic QF of Cauchy distribution`.
#' `s_halftriang()`: Basic QF of half-triangular distribution.
#' `s_halfcosine()`: Basic QF of half-cosine distribution.
#' `s_sech()`: Basic QF of hyperbolic secant distribution.
#' @param u depth
#' @param ... not used
#'
#' @return modified function
#' @rdname basicqf
#' @export
#' @examples
#' s_exp(0.6)
#' qexp(0.6, 1)
s_exp <- function(u, ...){
  -log(1-u)
}
class(s_exp) <- c("function", "qf", "basic")
attr(s_exp, "expects") <- "U"
attr(s_exp, "returns") <- NA

#' @rdname basicqf
#' @export
s_unif <- function(u, ...){
  u
}
class(s_unif) <- c("function", "qf", "basic")
attr(s_unif, "expects") <- "U"
attr(s_unif, "returns") <- NA

#' @rdname basicqf
#' @export
#' @importFrom stats qnorm
s_norm <- function(u, ...){
  qnorm(u)
}
class(s_norm) <- c("function", "qf", "basic")
attr(s_norm, "expects") <- "U"
attr(s_norm, "returns") <- NA

#' @rdname basicqf
#' @export
s_cauchy <- function(u, ...){
  tan(pi*(u-0.5))
}
class(s_cauchy) <- c("function", "qf", "basic")
attr(s_cauchy, "expects") <- "U"
attr(s_cauchy, "returns") <- NA

#' @rdname basicqf
#' @export
s_halftriang <- function(u, ...){
  -sqrt(1-u)
}
class(s_halftriang) <- c("function", "qf", "basic")
attr(s_halftriang, "expects") <- "U"
attr(s_halftriang, "returns") <- NA

#' @rdname basicqf
#' @export
s_halfcosine <- function(u, ...){
  asin(u)
}
class(s_halfcosine) <- c("function", "qf", "basic")
attr(s_halfcosine, "expects") <- "U"
attr(s_halfcosine, "returns") <- NA

#' @rdname basicqf
#' @export
s_sech <- function(u,...){
  2/pi*log(tan(pi/2*u))
}
class(s_sech) <- c("function", "qf", "basic")
attr(s_sech, "expects") <- "U"
attr(s_sech, "returns") <- NA