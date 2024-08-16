
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
