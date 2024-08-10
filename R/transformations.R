#' Q-transformations
#' @description
#' Some of the typical transformations of QFs, implementing a Q-transformation rule.
#'    `qtr_lehmann1()`: Raising of QF to a power using Lehman Type I inverse exponentiation. Returns \eqn{Q_1(u)^{1/k}}.
#'    `qtr_exp()`: Exponentiating the QF. Returns \eqn{k^{Q_1(u)}}. Default \eqn{k=e} Euler's constant
#'    `qtr_fun()`: Q-transform with generic function without additional arguments. \eqn{.fun(Q_1(u))}.
#'    `qtr_epsilon()`: unit-Q-transform using inverse epsilon function \eqn{\frac{(1+Q_1(u))^{1/\beta}-1}{(1+Q_1(u))^{1/\beta}+1}}.
#'
#' Note that today p-transformations can be performed by applying Q-transformations to standard uniform distribution
#' @param fun function
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
qtr_lehmann1 <- function(fun, nm_pow=".pow", .invert=TRUE){
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
#' qtr_epsilon(qnorm)
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
#' Some of the typical transformations of QFs, implementing a p-transformation rule.
#'
#'    - `ptr_lehmann1()`: Lehman Type I inverse exponentiation. Returns \eqn{u^{1/k}}.
#'    - `ptr_lehmann2()`: Lehman Type II inverse exponentiation. Returns \eqn{1-(1-u)^{1/k}}.
#'    - `ptr_exp()`: Exponentiating the QF. Returns \eqn{k^Q_1(u)}. The value of k defaults to Euler's constant.
#'    - `ptr_fun()`: Q-transform with generic function without additional arguments. \eqn{.fun(Q_1(u))}.
#'    - `ptr_KM()`: Kavya-Manoharan (KM) transformation \eqn{-\ln(1-u\frac{e-1}{e})}
#'    - `ptr_DUS()`: Dinesh-Umesh-Sunjay (DUS) transformation \eqn{\ln(1-u+eu)}.
#'    - `ptr_cDUS()`: Complimentary (reflected and shifted) Dinesh-Umesh-Sunjay (DUS) transformation \eqn{1-\ln(u-eu+e)}.
#'    - `ptr_modi()`: Modi transformation \eqn{\frac{u\alpha^\beta}{1-u+\alpha^\beta})}
#'
#' @param fun function
#' @param nm_pow character.  The name of the power parameter. The default name is `.pow`. The default value is 1
#' Should be a valid unique variable name other than "u"
#' @param .invert logical. Should the power parameter be inverted (1/.pow) before applying. Default TRUE
#' @return modified function
#' @rdname ptransformations
#' @export
#' @examples
#' qf_exp <- function(u)-log(1-u)
#' qf_weibull <- qtr_lehmann1(qf_exp, "k")
#' qf_weibull(0.5,k=1/5)
#' qweibull(0.5, shape = 5)
ptr_lehmann1 <- function(fun, nm_pow=".pow", .invert=TRUE){
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


#' @rdname ptransformations
#' @export
ptr_lehmann2 <- function(fun, nm_pow=".pow", .invert=TRUE){
  f <- function(u, .pow=1, ...){
    if(.invert) .pow <- 1/.pow
    fun(1-(1-u)^(.pow),...)
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
#' qf_lognorm <- qtr_exp(qf_norm)
#' qf_lognorm(0.2, mu=2, sigma=0.1)
#' qlnorm(0.2, 2, 0.1)
qtr_exp <- function(fun, nm_base=".base", .invert=TRUE){
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
ptr_exp <- function(fun, nm_base=".base", .invert=TRUE){
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
# Dinesh-Umesh-Sanjay (DUS) p-transformation
#' @rdname ptransformations
#' @export
ptr_DUS <- function(fun){
  function(u, ...){
    fun(log(1-u+exp(1)*u), ...)
  }
}

# Complimentary Dinesh-Umesh-Sanjay (DUS) p-transformation (reflected and shifted by 1)
#' @rdname ptransformations
#' @export
ptr_cDUS <- function(fun){
  function(u, ...){
    fun(1-log(u-exp(1)*u+exp(1)), ...)
  }
}

# Modi p-transformation
#' @rdname ptransformations
#' @export
ptr_modi <- function(fun){
  function(u, alpha, beta, ...){
    fun(u*alpha^beta/(1-u+alpha^beta), ...)
  }
}

# SHASH (sinh-asinh) q-transformation
#' @param nm_tail character.  The name of the tail thickness parameter. The default name is `.dlt`. 
#' The tail thickness parameter should be positive (default value is 1).
#' Should be a valid unique variable name other than "u"
#' The asymmetry parameter can be positive or negative (default value is 0).
#' @param nm_asymm character.  The name of the asymmetry parameter. The default name is `.eps`. 
#' Should be a valid unique variable name other than "u"
#' @rdname qtransformations
#' @export
#' @examples
#' qtr_shash(qnorm)
qtr_shash <- function(fun, nm_tail=".dlt", nm_asymm=".eps"){
  f <- function(u, .eps=0, .dlt=1, ...){
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
