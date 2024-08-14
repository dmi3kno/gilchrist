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

#' @param x numeric. Fixed value to shift/scale/power the QF by
#' @rdname qtransformations
#' @export
qtr_shiftby <- function(fun, x=0){
  function(u, ...){
    x+fun(u, ...)
  }
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
qtr_scaleby <- function(fun, x=1){
  function(u, ...){
    x*fun(u, ...)
  }
}

#' @rdname qtransformations
#' @export
qtr_powerby <- function(fun, x=1){
  function(u, ...){
    fun(u, ...)^x
  }
}

#' @rdname qtransformations
#' @export
qtr_oddITL <- function(fun){
  function(u, ...){
    sqrt(fun(u, ...)) / (1-sqrt(fun(u, ...)))
  }
}

#' p-transformations
#' @description
#' Some of the typical transformations of QFs, implementing a p-transformation rule.
#'
#'    - `ptr_lehmann1()`: Lehman Type I inverse exponentiation. (U->U). Returns \eqn{u^{1/k}}.
#'    - `ptr_lehmann2()`: Lehman Type II inverse exponentiation. (U->U). Returns \eqn{1-(1-u)^{1/k}}.
#'    - `ptr_fun()`: p-transform with generic function without additional arguments. (U->U, U->R) \eqn{Q_1(.fun(u))}.
#'    - `ptr_half()`: p-transform into half-distribution. (U->U) Returns \eqn{Q_1((u+1)/2))}.
#'    - `ptr_oddITL()`: Odd Inverse Topp-Leone transformation \eqn{\frac{\sqrt{u}}{1-\sqrt{u}}=\frac{u+\sqrt{u}}{1-u}=\frac{\sqrt{u}(1+\sqrt{u})}{1-u}=(u^{-1/2}-1)^{-1} }
#'    - `ptr_DUS()`: Dinesh-Umesh-Sunjay (DUS) transformation. (U->U) Returns \eqn{\ln(1-u+eu)}.
#'    - `ptr_KM()`: Kavya-Manoharan (KM) transformation. Equal to reflected and shifted DUS trasnformation (U->U). Returns \eqn{-\ln(1-u\frac{e-1}{e})}
#'    - `ptr_modi1()`: Modi transformation \eqn{\frac{u\alpha^\beta}{1-u+\alpha^\beta}}. Only alpha is mandatory, while beta defaults to 1
#'    - `ptr_modi2()`: Modi transformation \eqn{\frac{u+u\alpha^\beta}{u+\alpha^\beta})}. Only alpha is mandatory, while beta defaults to 1
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

#' @param .qf quantile function made with gilchrist (or wrapped basic function able to accept optional arguments through ellipsis) to be applied as p-transformation
#' @rdname ptransformations
#' @export
#' @examples
#' qtr_fun(sqf_exp,log1p)
ptr_qf <- function(fun, .qf){
  f <- function(u, ...)
    fun(.qf(u, ...), ...)
  f
}

#' @rdname ptransformations
#' @export
ptr_half <- function(fun){
  function(u, ...)
    fun((u+1)/2,...)
}

#' @param x numeric. Fixed value to shift/scale/power the u by
#' @rdname ptransformations
#' @export
ptr_shiftby <- function(fun, x=0){
  function(u, ...){
    fun(x+u, ...)
  }
}

#' @rdname ptransformations
#' @export
ptr_scaleby <- function(fun, x=1){
  function(u, ...){
    fun(x*u, ...)
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

#' @rdname ptransformations
#' @export
ptr_oddITL <- function(fun){
  function(u, ...){
    fun( sqrt(u)/(1-sqrt(u)), ...)
  }
}

# Modi (Type I) p-transformation
#' @param nm_a character.  The name of the Modi parameter `alpha`. The default name is `.modialpha`.
#' @param nm_b character.  The name of the Modi parameter `beta`. The default name is `.modibeta`. Default value is 1. 
#' @rdname ptransformations
#' @export
ptr_modi1 <- function(fun, nm_a=".modialpha", nm_b=".modibeta"){
  f <- function(u, .modialpha, .modibeta=1, ...){
    fun(u * .modialpha ^ .modibeta/(1-u+ .modialpha ^ .modibeta), ...)
  }

  formals_ <- formals(f)
  body_ <- body(f)
  names(formals_)[names(formals_) == ".modialpha"] <- nm_a
  names(formals_)[names(formals_) == ".modibeta"] <- nm_b
  body_ <- do.call(substitute, list(body_, list(.modialpha = as.symbol(nm_a))))
  body_ <- do.call(substitute, list(body_, list(.modibeta = as.symbol(nm_b))))
  as.function(c(formals_, body_))
}

# Modi Type II p-transformation
#' @rdname ptransformations
#' @export
ptr_modi2 <- function(fun, nm_a=".modialpha", nm_b=".modibeta"){
  f <- function(u, .modialpha, .modibeta=1, ...){
    fun((u + u * .modialpha ^ .modibeta)/(u + .modialpha ^ .modibeta), ...)
  }

  formals_ <- formals(f)
  body_ <- body(f)
  names(formals_)[names(formals_) == ".modialpha"] <- nm_a
  names(formals_)[names(formals_) == ".modibeta"] <- nm_b
  body_ <- do.call(substitute, list(body_, list(.modialpha = as.symbol(nm_a))))
  body_ <- do.call(substitute, list(body_, list(.modibeta = as.symbol(nm_b))))
  as.function(c(formals_, body_))
}
