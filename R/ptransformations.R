
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
#' @param pow character.  The name of the power parameter. The default name is `.pow`. 
#' @param nm_pow character.  The name of the power parameter. The default name is `.pow`. The default value is 1
#' Should be a valid unique variable name other than "u"
#' @param .invert logical. Should the power parameter be inverted (1/.pow) before applying. Default TRUE
#' @return modified function
#' @rdname ptransformations
#' @export
#' @examples
#' qf_exp <- as_qf(function(u)-log(1-u))
#' qf_weibull <- qtr_lehmann1(qf_exp, "k")
#' qf_weibull(0.5,k=1/5)
#' qweibull(0.5, shape = 5)
ptr_lehmann1 <- function(fun, nm_pow=".pow", pow=1, .invert=TRUE){
  stopifnot("ptr_lehmann1() is expecting a quantile function"=inherits(fun, c("function", "qf")))
  f <- function(u, .pow=pow, ...){
    if(.invert) .pow <- 1/.pow
    fun(u^(.pow),...)
  }
  math_x <- math(fun)
  math_y <- paste0(r"--(&^{)--", prmtr(nm_pow, .invert) ,r"--(})--")
  math_f <-  fn_insert(math_x, math_y, br=FALSE)
  
  formals_ <- formals(f)
  body_ <- body(f)
  names(formals_)[names(formals_) == ".pow"] <- nm_pow
  body_ <- do.call(substitute, list(body_, list(.pow = as.symbol(nm_pow))))
  as_qf(as.function(c(formals_, body_)), math=math_f)
}

#' @rdname ptransformations
#' @export
ptr_lehmann2 <- function(fun, nm_pow=".pow", pow=1, .invert=TRUE){
  stopifnot("ptr_lehmann2() is expecting a quantile function"=inherits(fun, c("function", "qf")))
  f <- function(u, .pow=pow, ...){
    if(.invert) .pow <- 1/.pow
    fun(1-(1-u)^(.pow),...)
  }
  math_x <- math(fun)
  math_y <- paste0(r"--(1-(1-&)^{)--", prmtr(nm_pow, .invert) ,r"--(})--")
  math_f <-  fn_insert(math_x, math_y)

  formals_ <- formals(f)
  body_ <- body(f)
  names(formals_)[names(formals_) == ".pow"] <- nm_pow
  body_ <- do.call(substitute, list(body_, list(.pow = as.symbol(nm_pow))))
  as_qf(as.function(c(formals_, body_)), math=math_f)
}

#' @param .fun function without arguments(or with all default arguments) to be applied as Q-transformation
#' @rdname ptransformations
#' @export
#' @examples
#' qtr_fun(qf_exp,log1p)
ptr_fun <- function(fun, .fun){
  stopifnot("ptr_fun() is expecting a quantile function"=inherits(fun, c("function", "qf")))
  f <- function(u, ...)
    fun(.fun(u),...)
  
  math_x <- math(fun)
  math_y <- paste0(r"--(\text{)--", deparse(substitute(.fun)), "}(&)")
  math_f <-  fn_insert(math_x, math_y)

  as_qf(f, math=math_f)
}



#' @param .qf quantile function made with gilchrist (or wrapped basic function able to accept optional arguments through ellipsis) to be applied as p-transformation
#' @rdname ptransformations
#' @export
#' @examples
#' qtr_fun(qf_exp,log1p)
ptr_qf <- function(fun, .qf){
  stopifnot("ptr_qf() is expecting a quantile function"=inherits(fun, c("function", "qf")))
  stopifnot(".qf in ptr_qf() should be a quantile function"=inherits(.qf, c("function", "qf")))
  f <- function(u, ...)
    fun(.qf(u, ...), ...)

  math_x <- math(fun)
  math_y <- math(.qf)
  math_f <- fn_insert(math_x, math_y)

  as_qf(f, math=math_f)
}

#' @rdname ptransformations
#' @export
ptr_half <- function(fun){
  stopifnot("ptr_half() is expecting a quantile function"=inherits(fun, c("function", "qf")))
  f <- function(u, ...)
    fun((u+1)/2,...)

  math_x <- math(fun)
  math_y <- paste0(r"--(\frac{&+1}{2})--")
  math_f <- fn_insert(math_x, math_y, br=FALSE)

  as_qf(f, math = math_f)
}

#' Alpha-Power Type 1 p-transformation
#' 
#' @param nm_ap character. Name of AP parameter. Default is .ap
#' @param ap numeric. Fixed value for AP parameter. Default is exp(1), DUS-transformation.
#' @rdname ptransformations
#' @export
ptr_AP1 <- function(fun, nm_ap=".ap", ap=exp(1)){
  stopifnot("ptr_AP1() is expecting a quantile function"=inherits(fun, c("function", "qf")))
  f <- function(u, .ap=ap, ...){
    fun(log(1+ .ap*u -u)/log(.ap), ...)
  }
  math_x <- math(fun)
  math_y <- paste0(r"--(\frac{\ln \left( 1 +)--", prmtr(nm_ap),
                    r"--(& - & \right)}{\ln( )--", prmtr(nm_ap), r"--( )})--")
  math_f <- fn_insert(math_x, math_y)

  formals_ <- formals(f)
  body_ <- body(f)
  names(formals_)[names(formals_) == ".ap"] <- nm_ap
  body_ <- do.call(substitute, list(body_, list(.ap = as.symbol(nm_ap))))
  as_qf(as.function(c(formals_, body_)), math=math_f)
}

#' Alpha-Power Type 2 p-transformation
#' @rdname ptransformations
#' @export
ptr_AP2 <- function(fun, nm_ap=".ap", ap=exp(1)){
  stopifnot("ptr_AP2() is expecting a quantile function"=inherits(fun, c("function", "qf")))
  f <- function(u, .ap=ap, ...){
    fun(1-log(u- .ap*u +.ap)/log(.ap), ...)
  }
  math_x <- math(fun)
  math_y <- paste0(r"--(1-\frac{\ln \left( & -)--", prmtr(nm_ap),
                    r"--(& +)--", prmtr(nm_ap) ,r"--(\right)}{\ln( )--", 
                    prmtr(nm_ap), r"--( )})--")
  math_f <- fn_insert(math_x, math_y)

  formals_ <- formals(f)
  body_ <- body(f)
  names(formals_)[names(formals_) == ".ap"] <- nm_ap
  body_ <- do.call(substitute, list(body_, list(.ap = as.symbol(nm_ap))))
  as_qf(as.function(c(formals_, body_)), math=math_f)
}

# Marshall-Olkin (MO) p-transformation
#' @param mopar numeric. Default value for MO parameter. Default is 1.
#' @param nm_mopar character. Name of the Marchall-Olking parameter. Default `.mopar`
#' @rdname ptransformations
#' @export
ptr_MO <- function(fun, nm_mopar=".mopar", mopar=1) {
  stopifnot("ptr_MO() is expecting a quantile function"=inherits(fun, c("function", "qf")))
  f <- function(u, .mopar=mopar, ...){
    stopifnot("Marshall-Olkin parameter should be positive"=.mopar>0)
    fun(.mopar*u/( 1-(1-.mopar)*u ), ...)
  }

  math_x <- math(fun)
  math_y <- paste0(r"--(\frac{)--", prmtr(nm_mopar) ,r"--{ &}{1-\left(1-}--", 
                    prmtr(nm_mopar),  r"--{ \right)&} }--")
  math_f <- fn_insert(math_x, math_y)

  formals_ <- formals(f)
  body_ <- body(f)
  names(formals_)[names(formals_) == ".mopar"] <- nm_mopar
  body_ <- do.call(substitute, list(body_, list(.mopar = as.symbol(nm_mopar))))
  as_qf(as.function(c(formals_, body_)), math = math_f)
}

# Generalized Marshall-Olkin (Harris) p-transformation
#' @param mopar numeric. Default value for MO parameter. Default is 1.
#' @param nm_mopar character. Name of the Marchall-Olking parameter. Default `.mopar`. 
#' Should be a valid unique variable name other than "u"
#' @param pow character.  The name of the power parameter. The default name is `.pow`. 
#' @param nm_pow character.  The name of the power parameter. The default name is `.pow`. The default value is 1
#' Should be a valid unique variable name other than "u"
#' @rdname ptransformations
#' @export
ptr_GMO <- function(fun, nm_mopar=".mopar", nm_pow=".pow", mopar=1, pow=1) {
  stopifnot("ptr_GMO() is expecting a quantile function"=inherits(fun, c("function", "qf")))
  f <- function(u, .mopar=mopar, .pow=pow, ...){
    stopifnot("Generalized Marshall-Olkin parameter should be positive"=.mopar>0)
    stopifnot("Generalized Marshall-Olkin power parameter should be positive"=.pow>0)
    fun(1-(1-u)/( .mopar+(1-.mopar)*(1-u)^.pow )^(1/.pow), ...)
  }

  math_x <- math(fun)
  math_y <- paste0(r"--{1-\frac{1-&}{\left( }--", prmtr(nm_mopar) ,r"--{+\left(1-}--", 
                    prmtr(nm_mopar),  r"--( \right)\left(1-& \right)^)--", prmtr(nm_pow), 
                    r"--( \right)^\frac{1}{ )--", prmtr(nm_pow), r"--( }})--")
  math_f <- fn_insert(math_x, math_y)

  formals_ <- formals(f)
  body_ <- body(f)
  names(formals_)[names(formals_) == ".mopar"] <- nm_mopar
  names(formals_)[names(formals_) == ".pow"] <- nm_pow
  body_ <- do.call(substitute, list(body_, list(.mopar = as.symbol(nm_mopar))))
  body_ <- do.call(substitute, list(body_, list(.pow = as.symbol(nm_pow))))
  as_qf(as.function(c(formals_, body_)), math = math_f)
}

# Topp-Leone (TL) p-transformation
#' @param tlpar numeric. Default value for MO parameter. Default is 1.
#' @param nm_tlpar character. Name of the Marchall-Olking parameter. Default `.tlpar`
#' @rdname ptransformations
#' @export
ptr_TL <- function(fun, nm_tlpar=".tlpar", tlpar=1) {
  stopifnot("ptr_TL() is expecting a quantile function"=inherits(fun, c("function", "qf")))
  f <- function(u, .tlpar=tlpar, ...){
    stopifnot("Topp-Leone parameter should be positive"=.tlpar>0)
    fun(1-sqrt(1-u^(1/.tlpar)), ...)
  }

  math_x <- math(fun)
  math_y <- paste0(r"--(1-\sqrt{1-&^\frac{1}{ )--", prmtr(nm_tlpar), r"--( } })--")
  math_f <- fn_insert(math_x, math_y)

  formals_ <- formals(f)
  body_ <- body(f)
  names(formals_)[names(formals_) == ".tlpar"] <- nm_tlpar
  body_ <- do.call(substitute, list(body_, list(.tlpar = as.symbol(nm_tlpar))))
  as_qf(as.function(c(formals_, body_)), math=math_f)
}

#' @rdname ptransformations
#' @export
ptr_oddITL <- function(fun){
  stopifnot("ptr_oddITL() is expecting a quantile function"=inherits(fun, c("function", "qf")))
  f <- function(u, ...){
    fun( sqrt(u)/(1-sqrt(u)), ...)
  }

  math_x <- math(fun)
  math_y <- paste0(r"--( \frac{\sqrt{&} }{1-\sqrt{&} } )--")
  math_f <- fn_insert(math_x, math_y)

  as_qf(f, math=math_f)
}

# Modi (Type I) p-transformation
#' @param a numeric. Default value of the Modi parameter `alpha`. Default is 1
#' @param b numeric. Default value of the Modi parameter `beta`. Default is 2
#' @param nm_a character.  The name of the Modi parameter `alpha`. The default name is `.modialpha`.
#' @param nm_b character.  The name of the Modi parameter `beta`. The default name is `.modibeta`. Default value is 1. 
#' @rdname ptransformations
#' @export
ptr_modi1 <- function(fun, nm_a=".modialpha", nm_b=".modibeta", a=1, b=1){
  stopifnot("ptr_modi1() is expecting a quantile function"=inherits(fun, c("function", "qf")))
  f <- function(u, .modialpha=a, .modibeta=b, ...){
    fun(u * .modialpha ^ .modibeta/(1-u+ .modialpha ^ .modibeta), ...)
  }

  math_x <- math(fun)
  math_y <- paste0(r"--{\frac{& }--",prmtr(nm_a), r"--{^}--", prmtr(nm_b), r"--{ }{ }--",
                    r"--{1- & +}--", prmtr(nm_a), r"--{^}--", prmtr(nm_b), r"--(})--")
  math_f <- fn_insert(math_x, math_y)

  formals_ <- formals(f)
  body_ <- body(f)
  names(formals_)[names(formals_) == ".modialpha"] <- nm_a
  names(formals_)[names(formals_) == ".modibeta"] <- nm_b
  body_ <- do.call(substitute, list(body_, list(.modialpha = as.symbol(nm_a))))
  body_ <- do.call(substitute, list(body_, list(.modibeta = as.symbol(nm_b))))
  as_qf(as.function(c(formals_, body_)), math = math_f)
}

# Modi Type II p-transformation
#' @rdname ptransformations
#' @export
ptr_modi2 <- function(fun, nm_a=".modialpha", nm_b=".modibeta", a=1, b=1){
  stopifnot("ptr_modi2() is expecting a quantile function"=inherits(fun, c("function", "qf")))
  f <- function(u, .modialpha=a, .modibeta=b, ...){
    fun((u + u * .modialpha ^ .modibeta)/(u + .modialpha ^ .modibeta), ...)
  }

  math_x <- math(fun)
  math_y <- paste0(r"--{\frac{ & + & }--",prmtr(nm_a), r"--{^}--", prmtr(nm_b), r"--{ }{ }--",
                  r"--{& +}--",prmtr(nm_a), r"--{^}--", prmtr(nm_b), r"--(})--")
  math_f <- fn_insert(math_x, math_y)

  formals_ <- formals(f)
  body_ <- body(f)
  names(formals_)[names(formals_) == ".modialpha"] <- nm_a
  names(formals_)[names(formals_) == ".modibeta"] <- nm_b
  body_ <- do.call(substitute, list(body_, list(.modialpha = as.symbol(nm_a))))
  body_ <- do.call(substitute, list(body_, list(.modibeta = as.symbol(nm_b))))
  as_qf(as.function(c(formals_, body_)), math = math_f)
}
