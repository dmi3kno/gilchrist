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
class(s_exp) <- c("qf", "function")
attr(s_exp, "expects") <- "U"
attr(s_exp, "returns") <- NA
attr(s_exp, "math") <- r"--{-\ln(1- {&} )}--"

#' @rdname basicqf
#' @export
s_unif <- function(u, ...){
  u
}
class(s_unif) <- c("qf", "function")
attr(s_unif, "expects") <- "U"
attr(s_unif, "returns") <- NA
attr(s_unif, "math") <- r"--{ {&} }--"

#' @rdname basicqf
#' @export
#' @importFrom stats qnorm
s_norm <- function(u, ...){
  qnorm(u)
}
class(s_norm) <- c("qf", "function")
attr(s_norm, "expects") <- "U"
attr(s_norm, "returns") <- NA
attr(s_norm, "math") <- r"--{\Phi^{-1}\left( {&} \right)}--"

#' @rdname basicqf
#' @export
s_cauchy <- function(u, ...){
  tan(pi*(u-0.5))
}
class(s_cauchy) <- c("qf", "function")
attr(s_cauchy, "expects") <- "U"
attr(s_cauchy, "returns") <- NA
attr(s_cauchy, "math") <- r"--{\tan\left( \pi\left( {&}-\frac{1}{2} \right)  \right) }--"

#' @rdname basicqf
#' @export
s_halftriang <- function(u, ...){
  -sqrt(1-u)
}
class(s_halftriang) <- c("qf", "function")
attr(s_halftriang, "expects") <- "U"
attr(s_halftriang, "returns") <- NA
attr(s_halftriang, "math") <- r"--{-\sqrt{ 1- {&}} }--"


#' @rdname basicqf
#' @export
s_halfcosine <- function(u, ...){
  asin(u)
}
class(s_halfcosine) <- c("qf", "function")
attr(s_halfcosine, "expects") <- "U"
attr(s_halfcosine, "returns") <- NA
attr(s_halfcosine, "math") <- r"--{ \arcsin{ {&} } }--"

#' @rdname basicqf
#' @export
s_sech <- function(u,...){
  2/pi*log(tan(pi/2*u))
}
class(s_sech) <- c("qf", "function")
attr(s_sech, "expects") <- "U"
attr(s_sech, "returns") <- NA
attr(s_sech, "math") <- r"--{ \frac{2}{\pi}\ln \left( \tan\left(\frac{\pi}{2} {&} \right) \right) }--"