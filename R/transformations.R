#' Q-transformations
#' @description
#' Some of the typical transformations of QFs, implementing a Q-transformation rule.
#'    `qtr_power()`: Raising of QF to a power. Returns \eqn{Q_1(u)^k}.
#'    `qtr_exponentiate()`: Exponentiating the QF. Returns \eqn{k^Q_1(u)}.
#'    `qtr_fun()`: Q-transform with generic function without additional arguments. \eqn{.fun(Q_1(u))}.
#'    `qtr_epsilon()`: unit-Q-transform using inverse epsilon function \eqn{\frac{(1+Q_1(u))^{1/\beta}-1}{(1+Q_1(u))^{1/\beta}+1}}.
#'
#' Note that today p-transformations can be performed by applying Q-transformations to standard uniform distribution
#' @param fun function
#' @param nm_pow character.  The name of the power parameter. The default name is `.pow`. The default value is 1
#' Should be a valid unique variable name other than "u"
#' @param .invert logical. Should the power parameter be inverted (1/.pow) before applying. Default FALSE
#' @return modified function
#' @rdname qtransformations
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

#' @rdname qtransformations
#' @export
#' @examples
qtr_epsilon <- function(fun, nm_pow=".pow"){
  f <- function(u, .pow=1, ...){
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

#' p-transformations
#' @description
#' Some of the typical transformations of QFs, implementing a Q-transformation rule.
#'    `ptr_power()`: Raising of QF to a power. Returns \eqn{Q_1(u)^k}.
#'    `ptr_ipower()`: Raising of QF to an inverse power. Returns \eqn{Q_1(u)^{1/k}}.
#'    `ptr_exponentiate()`: Exponentiating the QF. Returns \eqn{k^Q_1(u)}.
#'    `ptr_fun()`: Q-transform with generic function without additional arguments. \eqn{.fun(Q_1(u))}.
#'    `ptr_KM()`: Kavya-Manoharan (KM) transformation \eqn{-\ln(1-u\frac{e-1}{e})}
#'
#' @param fun function
#' @param nm_pow character.  The name of the power parameter. The default name is `.pow`. The default value is 1
#' Should be a valid unique variable name other than "u"
#' @param .invert logical. Should the power parameter be inverted (1/.pow) before applying. Default FALSE
#' @return modified function
#' @rdname ptransformations
#' @export
#' @examples
#' qf_exp <- function(u)-log(1-u)
#' qf_weibull <- qtr_power(qf_exp, "k")
#' qf_weibull(0.5,k=1/5)
#' qweibull(0.5, shape = 5)
ptr_power <- function(fun, nm_pow=".pow", .invert=FALSE){
  f <- function(u, .pow=1, ...){
    if(.invert) .pow <- 1/.pow
    fun(u^(.pow),...)
  }

  formals_ <- formals(f)
  body_ <- body(f)
  names(formals_)[names(formals_) == ".pow"] <- nm_pow
  body_ <- do.call(substitute, list(body_, list(.pow = as.symbol(nm_pow))))
  as.function(c(formals_, body_))
}

#' @param nm_base character.  The name of the base parameter. The default name is `.base`. The default value is `exp(1)` (Euler's constant).
# Should be a valid unique variable name other than "u"
#' @rdname qtransformations
#' @export
#' @examples
#' qf_norm <- qff_decorate(qnorm, nm_location="mu", nm_scale="sigma")
#' qf_lognorm <- qtr_exponentiate(qf_norm)
#' qf_lognorm(0.2, mu=2, sigma=0.1)
#' qlnorm(0.2, 2, 0.1)
qtr_exponentiate <- function(fun, nm_base=".base", .invert=FALSE){
  f <- function(u, .base=exp(1), ...){
    if(.invert) .base <- 1/.base
    (.base)^fun(u,...)
  }

  formals_ <- formals(f)
  body_ <- body(f)
  names(formals_)[names(formals_) == ".base"] <- nm_base
  body_ <- do.call(substitute, list(body_, list(.base = as.symbol(nm_base))))
  as.function(c(formals_, body_))
}

#' @param nm_base character.  The name of the base parameter. The default name is `.base`. The default value is `exp(1)` (Euler's constant).
# Should be a valid unique variable name other than "u"
#' @rdname ptransformations
#' @export
ptr_exponentiate <- function(fun, nm_base=".base", .invert=FALSE){
  f <- function(u, .base=exp(1), ...){
   if(.invert) .base <- 1/.base
   fun((.base)^u,...)
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
#' qtr_fun(sqf_exp,log1p)
qtr_fun <- function(fun, .fun){
  f <- function(u, ...)
    .fun(fun(u,...))
  f
}

#' @param .fun function without arguments(or with all default arguments) to be applied as Q-transformation
#' @rdname ptransformations
#' @export
#' @examples
#' qtr_fun(sqf_exp,log1p)
ptr_fun <- function(fun, .fun){
  f <- function(u, ...)
    fun(.fun(u),...)
  f
}


#' @param x numeric. Fixed value to shift/scale/power the QF by
#' @rdname qtransformations
#' @export
qtr_shiftby <- function(fun, x=0){
  function(u, ...){
    x+fun(u, ...)
  }
}

#' @param x numeric. Fixed value to shift/scale/power the u by
#' @rdname ptransformations
#' @export
ptr_shiftby <- function(fun, x=0){
  function(u, ...){
    fun(x+u, ...)
  }
}

#' @rdname qtransformations
#' @export
qtr_scaleby <- function(fun, x=1){
  function(u, ...){
    x*fun(u, ...)
  }
}

#' @rdname ptransformations
#' @export
ptr_scaleby <- function(fun, x=1){
  function(u, ...){
    fun(x*u, ...)
  }
}

#' @rdname qtransformations
#' @export
qtr_powerby <- function(fun, x=1){
  function(u, ...){
    fun(u, ...)^x
  }
}

#' @rdname ptransformations
#' @export
ptr_powerby <- function(fun, x=1){
  function(u, ...){
    fun(u^x, ...)
  }
}

# Kavya-Manoharan (KM) p-transformation
#' @rdname ptransformations
#' @export
ptr_KM <- function(fun){
  function(u, ...){
    em1e <- expm1(1)/exp(1)
    fun(-log(1-u*em1e), ...)
  }
}
