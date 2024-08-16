#' Q-transformations
#' @description
#' Some of the typical transformations of QFs, implementing a Q-transformation rule.
#'    `qtr_lehmann1()`: Raising of QF to a power using Lehman Type I inverse exponentiation. Returns \eqn{Q_1(u)^{1/k}}.
#'    `qtr_exponentiate()`: Exponentiating the QF. Returns \eqn{k^{Q_1(u)}}. Default \eqn{k=e} Euler's constant
#'    `qtr_fun()`: Q-transform with generic function without additional arguments. \eqn{.fun(Q_1(u))}.
#'    `qtr_epsilon()`: unit-Q-transform using inverse epsilon function \eqn{\frac{(1+Q_1(u))^{1/\beta}-1}{(1+Q_1(u))^{1/\beta}+1}}.
#'    `qtr_shash()`: SHASH (sinh-asinh) q-transformation. \eqn{\text{sinh}(1/\delta(Q(u) - \epsilon)}
#'
#' Note that today p-transformations can be performed by applying Q-transformations to standard uniform distribution
#' @param fun function
#' @param pow numeric. Fixed value for the power parameter. Default is 1
#' @param nm_pow character.  The name of the power parameter. The default name is `.pow`. The default value is 1
#' Should be a valid unique variable name other than "u"
#' @param .invert logical. Should the power parameter be inverted (1/.pow) before applying. Default TRUE
#' @return modified function
#' @rdname qtransformations
#' @export
#' @examples
#' qf_exp <- function(u)-log(1-u)
#' qf_weibull <- qtr_lehmann1(qf_exp, "k")
#' qf_weibull(0.5, k = 1/5)
#' qweibull(0.5, shape = 5)
qtr_lehmann1 <- function(fun, nm_pow=".pow", pow=1, .invert=TRUE){
  f <- function(u, .pow=pow, ...){
    if(.invert) .pow <- 1/.pow
    fun(u,...)^(.pow)
  }

  formals_ <- formals(f)
  body_ <- body(f)
  names(formals_)[names(formals_) == ".pow"] <- nm_pow
  body_ <- do.call(substitute, list(body_, list(.pow = as.symbol(nm_pow))))
  as.function(c(formals_, body_))
}

#' @rdname qtransformations
#' @export
qtr_reflect_shift <- function(fun){
  function(u, ...){
    -fun(1-u, ...) + 1
  }
}

#' @rdname qtransformations
#' @export
qtr_shift_reciprocate <- function(fun){
  function(u, ...){
    1/(1+fun(1-u, ...))
  }
}

#' @rdname qtransformations
#' @export
qtr_odd <- function(fun){
  function(u, ...){
    fun(u, ...)/(1+fun(u, ...))
  }
}

#' @param pow numeric. Fixed value for the power parameter. Default is 1
#' @rdname qtransformations
#' @export
#' @examples
#' qtr_epsilon(qnorm)
qtr_epsilon <- function(fun, nm_pow=".pow", pow=1){
  f <- function(u, .pow=pow, ...){
    x <- fun(u,...)
    ((1+x)^(1/.pow)-1)/
      ((1+x)^(1/.pow)+1)
  }

  formals_ <- formals(f)
  body_ <- body(f)
  names(formals_)[names(formals_) == ".pow"] <- nm_pow
  body_ <- do.call(substitute, list(body_, list(.pow = as.symbol(nm_pow))))
  as.function(c(formals_, body_))
}

#' @param tail numeric. Fixed value for the tail parameter. Default is 1
#' @param nm_tail character.  The name of the tail thickness parameter. The default name is `.dlt`. 
#' The tail thickness parameter should be positive (default value is 1).
#' Should be a valid unique variable name other than "u"
#' The asymmetry parameter can be positive or negative (default value is 0).
#' @param asymm numeric. Default value for 
#' @param nm_asymm character.  The name of the asymmetry parameter. The default name is `.eps`. 
#' Should be a valid unique variable name other than "u"
#' @rdname qtransformations
#' @export
#' @examples
#' qtr_shash(qnorm)
qtr_shash <- function(fun, nm_tail=".dlt", nm_asymm=".eps", tail=1, asymm=0){
  f <- function(u, .eps=asymm, .dlt=tail, ...){
    sinh(1/.dlt*(fun(u,...)- .eps))
  }

  formals_ <- formals(f)
  body_ <- body(f)
  names(formals_)[names(formals_) == ".dlt"] <- nm_tail
  names(formals_)[names(formals_) == ".eps"] <- nm_asymm
  body_ <- do.call(substitute, list(body_, list(.dlt = as.symbol(nm_tail))))
  body_ <- do.call(substitute, list(body_, list(.eps = as.symbol(nm_asymm))))
  as.function(c(formals_, body_))
}

#' @param base numeric. Fixed value of the base parameter. The default value is exp(1) (Eulers constant).
#' @param nm_base character. The name of the base parameter.
#' Should be a valid unique variable name other than "u"
#' @rdname qtransformations
#' @export
#' @examples
#' qf_norm <- qtr_decorate(qnorm, nm_location="mu", nm_scale="sigma")
#' qf_lognorm <- qtr_exponentiate(qf_norm)
#' qf_lognorm(0.2, mu=2, sigma=0.1)
#' qlnorm(0.2, 2, 0.1)
qtr_exponentiate <- function(fun, nm_base=".base", base=exp(1), .invert=FALSE){
  f <- function(u, .base=base, ...){
    if(.invert) .base <- 1/.base
    (.base)^fun(u,...)
  }

  formals_ <- formals(f)
  body_ <- body(f)
  names(formals_)[names(formals_) == ".base"] <- nm_base
  body_ <- do.call(substitute, list(body_, list(.base = as.symbol(nm_base))))
  as.function(c(formals_, body_))
}


#' @param .fun function without arguments(or with all default arguments) to be applied as Q-transformation
#' @rdname qtransformations
#' @export
#' @examples
#' qtr_fun(sqf_exp, log1p)
qtr_fun <- function(fun, .fun){
  f <- function(u, ...)
    .fun(fun(u,...))
  f
}

#' @param .qf quantile function made with gilchrist (or regular function wrapped to safely accept optional arguments through ellipsis) to be applied as Q-transformation
#' @rdname qtransformations
#' @export
#' @examples
#' qtr_fun(sqf_exp,log1p)
qtr_qf <- function(fun, .qf){
  f <- function(u, ...)
    .qf(fun(u,...), ...)
  f
}

#' @rdname qtransformations
#' @export
qtr_oddITL <- function(fun){
  function(u, ...){
    sqrt(fun(u, ...)) / (1-sqrt(fun(u, ...)))
  }
}
