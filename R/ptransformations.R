
#' p-transformations
#' @description
#' Some of the typical transformations of QFs, implementing a p-transformation rule.
#' \itemize{
#'    \item `ptr_lehmann1()`: Lehman Type I inverse exponentiation. (U->U). Returns \eqn{u^{1/k}}.
#'    \item `ptr_lehmann2()`: Lehman Type II inverse exponentiation. (U->U). Returns \eqn{1-(1-u)^{1/k}}.
#'    \item `ptr_fun()`: p-transform with generic function without additional arguments. (U->U, U->R) \eqn{Q_1(.fun(u))}.
#'    \item `ptr_half()`: p-transform into half-distribution. (U->U) Returns \eqn{Q_1((u+1)/2))}.
#'    \item `ptr_oddITL()`: Odd Inverse Topp-Leone transformation \eqn{\frac{\sqrt{u}}{1-\sqrt{u}}=\frac{u+\sqrt{u}}{1-u}=\frac{\sqrt{u}(1+\sqrt{u})}{1-u}=(u^{-1/2}-1)^{-1} }
#'    \item `ptr_DUS()`: Dinesh-Umesh-Sunjay (DUS) transformation. (U->U) Returns \eqn{\ln(1-u+eu)}.
#'    \item `ptr_KM()`: Kavya-Manoharan (KM) transformation. Equal to reflected and shifted DUS trasnformation (U->U). Returns \eqn{-\ln(1-u\frac{e-1}{e})}
#'    \item `ptr_modi1()`: Modi transformation \eqn{\frac{u\alpha^\beta}{1-u+\alpha^\beta}}. Only alpha is mandatory, while beta defaults to 1
#'    \item `ptr_modi2()`: Modi transformation \eqn{\frac{u+u\alpha^\beta}{u+\alpha^\beta})}. Only alpha is mandatory, while beta defaults to 1
#'    \item `ptr_arcsin()`: Arcsine transform \eqn{\frac{2}{\pi}\arcsin{u})}.
#'}
#' @param fun function
#' @param pow character.  The name of the power parameter. The default name is `.pow`.
#' @param nm_pow character.  The name of the power parameter. The default name is `.pow`. The default value is 1
#' Should be a valid unique variable name other than "u"
#' @param .invert logical. Should the power parameter be inverted (1/.pow) before applying. Default TRUE
#' @param pfn_pow function. Parameter transforming function for parameter pow
#' @return modified function
#' @md
#' @rdname ptransformations
#' @export
#' @examples
#' qf_exp <- as_qf(function(u)-log(1-u))
#' qf_weibull <- qtr_lehmann1(qf_exp, "k")
#' qf_weibull(0.5,k=1/5)
#' qweibull(0.5, shape = 5)
ptr_lehmann1 <- function(fun, nm_pow=".pow", pow=1, .invert=TRUE, pfn_pow=NULL){
  stopifnot("ptr_lehmann1() is expecting a quantile function"=inherits(fun, c("function", "qf")))
  f <- function(u, .pow=pow, ...){
    .pow <- prm_tr(.pow, .invert, pfn_pow)
    stopifnot("Power parameter should be non-negative"=.pow>=0)
    fun(u^(.pow),...)
  }
  math_x <- math(fun)
  math_y <- paste0(r"--(&^{)--", prm(nm_pow, .invert, pfn_pow) ,r"--(})--")
  math_f <-  fn_insert(math_x, math_y, br=FALSE)

  formals_ <- formals(f)
  body_ <- body(f)
  names(formals_)[names(formals_) == ".pow"] <- nm_pow
  body_ <- do.call(substitute, list(body_, list(.pow = as.symbol(nm_pow))))
  as_qf(as.function(c(formals_, body_)), math=math_f)
}

#' @rdname ptransformations
#' @export
ptr_lehmann2 <- function(fun, nm_pow=".pow", pow=1, .invert=TRUE, pfn_pow=NULL){
  stopifnot("ptr_lehmann2() is expecting a quantile function"=inherits(fun, c("function", "qf")))
  f <- function(u, .pow=pow, ...){
    .pow <- prm_tr(.pow, .invert, pfn_pow)
    stopifnot("Power parameter should be non-negative"=.pow>=0)
    fun(1-(1-u)^(.pow),...)
  }
  math_x <- math(fun)
  math_y <- paste0(r"--(1-(1-&)^{)--", prm(nm_pow, .invert, pfn_pow) ,r"--(})--")
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
#' @param ap numeric. Fixed value for AP parameter. Default is 1, DUS-transformation.
#' @param pfn_ap parameter transforming function. Default is exp()
#' @rdname ptransformations
#' @export
ptr_AP1 <- function(fun, nm_ap=".ap", ap=1, pfn_ap=exp){
  stopifnot("ptr_AP1() is expecting a quantile function"=inherits(fun, c("function", "qf")))
  f <- function(u, .ap=ap, ...){
    .ap <- prm_tr(.ap, FALSE, pfn_ap)
    stopifnot("Expecting positive ap not equal to 1"=(.ap != 1 && .ap>0))
    fun(log(1+ .ap*u -u)/log(.ap), ...)
  }
  math_x <- math(fun)
  math_y <- paste0(r"--(\frac{\ln \left( 1 +)--", prm(nm_ap, FALSE, pfn_ap),
                    r"--(& - & \right)}{\ln( )--", prm(nm_ap, FALSE, pfn_ap), r"--( )})--")
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
ptr_AP2 <- function(fun, nm_ap=".ap", ap=1, pfn_ap=exp){
  stopifnot("ptr_AP2() is expecting a quantile function"=inherits(fun, c("function", "qf")))
  f <- function(u, .ap=ap, ...){
    .ap <- prm_tr(.ap, FALSE, pfn_ap)
    stopifnot("Expecting positive ap not equal to 1"=(.ap != 1 && .ap>0))
    fun(1-log(u- .ap*u +.ap)/log(.ap), ...)
  }
  math_x <- math(fun)
  math_y <- paste0(r"--(1-\frac{\ln \left( & -)--", prm(nm_ap, FALSE, pfn_ap),
                    r"--(& +)--", prm(nm_ap, FALSE, pfn_ap) ,r"--(\right)}{\ln( )--",
                    prm(nm_ap, FALSE, pfn_ap), r"--( )})--")
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
#' @param pfn_mopar parameter transforming function. Default in none.
#' @rdname ptransformations
#' @export
ptr_MO <- function(fun, nm_mopar=".mopar", mopar=1, pfn_mopar=NULL) {
  stopifnot("ptr_MO() is expecting a quantile function"=inherits(fun, c("function", "qf")))
  f <- function(u, .mopar=mopar, ...){
    .mopar <- prm_tr(.mopar, FALSE, pfn_mopar)
    stopifnot("Marshall-Olkin parameter should be positive"=.mopar>0)
    fun(.mopar*u/( 1-(1-.mopar)*u ), ...)
  }

  math_x <- math(fun)
  math_y <- paste0(r"--(\frac{)--", prm(nm_mopar, FALSE, pfn_mopar) ,r"--{ &}{1-\left(1-}--",
                    prm(nm_mopar, FALSE, pfn_mopar),  r"--{ \right)&} }--")
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
ptr_GMO <- function(fun, nm_mopar=".mopar", nm_pow=".pow", mopar=1, pow=1, pfn_mopar=NULL, pfn_pow=NULL) {
  stopifnot("ptr_GMO() is expecting a quantile function"=inherits(fun, c("function", "qf")))
  f <- function(u, .mopar=mopar, .pow=pow, ...){
    .mopar <- prm_tr(.mopar, FALSE, pfn_mopar)
    .pow <- prm_tr(.pow, FALSE, pfn_pow)
    stopifnot("Generalized Marshall-Olkin parameter should be positive"=.mopar>0)
    stopifnot("Generalized Marshall-Olkin power parameter should be positive"=.pow>0)
    fun(1-(1-u)/( .mopar+(1-.mopar)*(1-u)^.pow )^(1/.pow), ...)
  }

  math_x <- math(fun)
  math_y <- paste0(r"--{1-\frac{1-&}{\left( }--", prm(nm_mopar, FALSE, pfn_mopar) ,r"--{+\left(1-}--",
                    prm(nm_mopar, FALSE, pfn_mopar),  r"--( \right)\left(1-& \right)^)--",
                   prm(nm_pow, FALSE, pfn_pow), r"--( \right)^\frac{1}{ )--",
                   prm(nm_pow, FALSE, pfn_pow), r"--( }})--")
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
#' @param pfn_tlpar parameter transforming function. Default is none.
#' @rdname ptransformations
#' @export
ptr_TL <- function(fun, nm_tlpar=".tlpar", tlpar=1, pfn_tlpar=NULL) {
  stopifnot("ptr_TL() is expecting a quantile function"=inherits(fun, c("function", "qf")))
  f <- function(u, .tlpar=tlpar, ...){
    .tlpar <- prm_tr(.tlpar, TRUE, pfn_tlpar)
    stopifnot("Topp-Leone parameter should be positive"=.tlpar>0)
    fun(1-sqrt(1-u^(.tlpar)), ...)
  }

  math_x <- math(fun)
  math_y <- paste0(r"--(1-\sqrt{1-&^)--", prm(nm_tlpar, TRUE, pfn_tlpar), r"--( })--")
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

#' @rdname ptransformations
#' @export
ptr_arcsin <- function(fun){
  stopifnot("ptr_arcsin() is expecting a quantile function"=inherits(fun, c("function", "qf")))
  f <- function(u, ...){
    fun( 2/pi*asin(u), ...)
  }

  math_x <- math(fun)
  math_y <- paste0(r"--( \frac{2}{\pi}\arcsin\left( & \right) )--")
  math_f <- fn_insert(math_x, math_y)

  as_qf(f, math=math_f)
}

# Modi (Type I) p-transformation
#' @param a numeric. Default value of the Modi parameter `alpha`. Default is 1
#' @param b numeric. Default value of the Modi parameter `beta`. Default is 2
#' @param nm_a character. The name of the Modi parameter `alpha`. The default name is `.modialpha`.
#' @param nm_b character. The name of the Modi parameter `beta`. The default name is `.modibeta`. Default value is 1.
#' @param pfn_a function. Parameter transforming function. Default is none.
#' @param pfn_b function. Parameter transforming function. Default is none.
#' @rdname ptransformations
#' @export
ptr_modi1 <- function(fun, nm_a=".modialpha", nm_b=".modibeta", a=1, b=1, pfn_a=NULL, pfn_b=NULL){
  stopifnot("ptr_modi1() is expecting a quantile function"=inherits(fun, c("function", "qf")))
  f <- function(u, .modialpha=a, .modibeta=b, ...){
    .modialpha <- prm_tr(.modialpha, FALSE, pfn_a)
    .modibeta <- prm_tr(.modibeta, FALSE, pfn_b)
    fun(u * .modialpha ^ .modibeta/(1-u+ .modialpha ^ .modibeta), ...)
  }

  math_x <- math(fun)
  math_y <- paste0(r"--{\frac{& }--",prm(nm_a, FALSE, pfn_a), r"--{^}--", prm(nm_b, FALSE, pfn_b), r"--{ }{ }--",
                    r"--{1- & +}--", prm(nm_a, FALSE, pfn_a), r"--{^}--", prm(nm_b, FALSE, pfn_b), r"--(})--")
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
ptr_modi2 <- function(fun, nm_a=".modialpha", nm_b=".modibeta", a=1, b=1, pfn_a=NULL, pfn_b=NULL){
  stopifnot("ptr_modi2() is expecting a quantile function"=inherits(fun, c("function", "qf")))
  f <- function(u, .modialpha=a, .modibeta=b, ...){
    .modialpha <- prm_tr(.modialpha, FALSE, pfn_a)
    .modibeta <- prm_tr(.modibeta, FALSE, pfn_b)
    fun((u + u * .modialpha ^ .modibeta)/(u + .modialpha ^ .modibeta), ...)
  }

  math_x <- math(fun)
  math_y <- paste0(r"--{\frac{ & + & }--",prm(nm_a, FALSE, pfn_a), r"--{^}--", prm(nm_b, FALSE, pfn_b), r"--{ }{ }--",
                  r"--{& +}--",prm(nm_a, FALSE, pfn_a), r"--{^}--", prm(nm_b, FALSE, pfn_b), r"--(})--")
  math_f <- fn_insert(math_x, math_y)

  formals_ <- formals(f)
  body_ <- body(f)
  names(formals_)[names(formals_) == ".modialpha"] <- nm_a
  names(formals_)[names(formals_) == ".modibeta"] <- nm_b
  body_ <- do.call(substitute, list(body_, list(.modialpha = as.symbol(nm_a))))
  body_ <- do.call(substitute, list(body_, list(.modibeta = as.symbol(nm_b))))
  as_qf(as.function(c(formals_, body_)), math = math_f)
}
