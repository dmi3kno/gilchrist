#' Q-transformations
#' @description
#' Some of the typical transformations of QFs, implementing a Q-transformation rule.
#'    `qtr_power()`: Raising of QF to a power. Returns \eqn{Q_1(u)^k}.
#'    `qtr_ipower()`: Raising of QF to an inverse power. Returns \eqn{Q_1(u)^{1/k}}.
#'    `qtr_exponentiate()`: Exponentiating the QF. Returns \eqn{k^Q_1(u)}.
#'    `qtr_fun()`: Q-transform with generic function without additional arguments. \eqn{.fun(Q_1(u))}.
#'
#' Note that today p-transformations can be performed by applying Q-transformations to standard uniform distribution
#' @param fun function
#' @param nm_pow character.  The name of the power parameter. The default name is `.pow`. The default value is 1
#' Should be a valid unique variable name other than "u"
#' @param .invert logical. Should the power parameter be inverted (1/.pow) before applying. Default FALSE
#' @return modified function
#' @rdname transformations
#' @export
#' @examples
#' qf_exp <- function(u)-log(1-u)
#' qf_weibull <- qtr_power(qf_exp, "k")
#' qf_weibull(0.5,k=1/5)
#' qweibull(0.5, shape = 5)
qtr_power <- function(fun, nm_pow=".pow", .invert=FALSE){
  f <- function(u, .pow=1, ...){
    if(.invert) .pow <- 1/.pow
    fun(u,...)^(.pow)
  }


  formals_ <- formals(f)
  body_ <- body(f)
  names(formals_)[names(formals_) == ".pow"] <- nm_pow
  body_ <- do.call(substitute, list(body_, list(.pow = as.symbol(nm_pow))))
  as.function(c(formals_, body_))
}

#' @param nm_base character.  The name of the base parameter. The default name is `.base`. The default value is `exp(1)` (Euler's constant).
# Should be a valid unique variable name other than "u"
#' @rdname transformations
#' @export
#' @examples
#' qf_norm <- qff_decorate(qnorm, nm_location="mu", nm_scale="sigma")
#' qf_lognorm <- qtr_exponentiate(qf_norm)
#' qf_lognorm(0.2, mu=2, sigma=0.1)
#' qlnorm(0.2, 2, 0.1)
qtr_exponentiate <- function(fun, nm_base=".base"){
  f <- function(u, .base=exp(1), ...)
    (.base)^fun(u,...)

  formals_ <- formals(f)
  body_ <- body(f)
  names(formals_)[names(formals_) == ".base"] <- nm_base
  body_ <- do.call(substitute, list(body_, list(.base = as.symbol(nm_base))))
  as.function(c(formals_, body_))
}


#' @param .fun function without arguments(or with all default arguments) to be applied as Q-transformation
#' @rdname transformations
#' @export
#' @examples
#' qtr_fun(sqf_exp,log1p)
qtr_fun <- function(fun, .fun){
  f <- function(u, ...)
    .fun(fun(u,...))
  f
}


#' @param x numeric. Fixed value to shift/scale/power the QF by
#' @rdname transformations
#' @export
qtr_shiftby <- function(fun, x=0){
  function(u, ...){
    x+fun(u, ...)
  }
}

#' @rdname transformations
#' @export
qtr_scaleby <- function(fun, x=1){
  function(u, ...){
    x*fun(u, ...)
  }
}

#' @rdname transformations
#' @export
qtr_powerby <- function(fun, x=1){
  function(u, ...){
    fun(u, ...)^x
  }
}

