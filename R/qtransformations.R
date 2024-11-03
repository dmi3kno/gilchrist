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
#' @param pfn_pow function. Parameter transforming function for power argument. Default is none.
#' @return modified function
#' @rdname qtransformations
#' @export
#' @examples
#' qf_exp <- function(u)-log(1-u)
#' qf_weibull <- qtr_lehmann1(qf_exp, "k")
#' qf_weibull(0.5, k = 1/5)
#' qweibull(0.5, shape = 5)
qtr_lehmann1 <- function(fun, nm_pow=".pow", pow=1, .invert=TRUE, pfn_pow=NULL){
  stopifnot("qtr_lehmann1() is expecting a quantile function"=inherits(fun, c("function", "qf")))
  f <- function(u, .pow=pow, ...){
    .pow <- prm_tr(.pow, .invert, pfn_pow)
    stopifnot("Power parameter should be non-negative"=.pow>=0)
    fun(u,...)^(.pow)
  }
  math_y <- math(fun)
  math_x <- paste0(r"--(&^{)--", prm(nm_pow, .invert, pfn_pow) ,r"--(})--")
  math_f <-  fn_insert(math_x, math_y)

  formals_ <- formals(f)
  body_ <- body(f)
  names(formals_)[names(formals_) == ".pow"] <- nm_pow
  body_ <- do.call(substitute, list(body_, list(.pow = as.symbol(nm_pow))))
  as_qf(as.function(c(formals_, body_)), math = math_f)
}



#' @rdname qtransformations
#' @export
qtr_reflect_shift <- function(fun){
  stopifnot("qtr_reflect_shift() is expecting a quantile function"=inherits(fun, c("function", "qf")))
  f <- function(u, ...){
    -fun(1-u, ...) + 1
  }

  math_y <- math(fun)
  math_p <- paste0("1 - &")
  math_x <- paste0("- & +1")
  math_f <- fn_insert(fn_insert(math_x, math_y), math_p)

  as_qf(f, math = math_f)
}

#' Also know as Odd-Inverse transformation
#' @rdname qtransformations
#' @export
qtr_shift_reciprocate <- function(fun){
  stopifnot("qtr_shift_reciprocate() is expecting a quantile function"=inherits(fun, c("function", "qf")))
  f <- function(u, ...){
    1/(1+fun(1-u, ...))
  }

  math_y <- math(fun)
  math_p <- paste0("1 - &")
  math_x <- paste0(r"--(\frac{1}{& + 1})--")
  math_f <- fn_insert(fn_insert(math_x, math_y), math_p)
  as_qf(f, returns="U", math = math_f)
}

#' @param nm_offset character. Name of the offset argument. Default is `.offset`
#' @param offset numeric. Default value for odd offset. Defaults to 1
#' @param pfn_offset parameter transforming function for offset argument. Default is none.
#' @rdname qtransformations
#' @export
qtr_odd <- function(fun, nm_offset=".offset", offset=1, pfn_offset=NULL){
  stopifnot("qtr_odd() is expecting a quantile function"=inherits(fun, c("function", "qf")))
  f <- function(u, .offset=offset, ...){
    .offset <- prm_tr(.offset, FALSE, pfn_offset)
    fun(u, ...)/(.offset + fun(u, ...))
  }

  math_y <- math(fun)
  math_x <- paste0(r"--(\frac{&}{ )--", prm(nm_offset, FALSE, pfn_offset), r"--(+ & })--")
  math_f <-  fn_insert(math_x, math_y)

  formals_ <- formals(f)
  body_ <- body(f)
  names(formals_)[names(formals_) == ".offset"] <- nm_offset
  body_ <- do.call(substitute, list(body_, list(.offset = as.symbol(nm_offset))))
  as_qf(as.function(c(formals_, body_)), returns="U", math = math_f)
}

#' @param pow numeric. Fixed value for the power parameter. Default is 1
#' @rdname qtransformations
#' @export
#' @examples
#' qtr_epsilon(qnorm)
qtr_epsilon <- function(fun, nm_pow=".pow", pow=1, pfn_pow=NULL){
  stopifnot("qtr_epsilon() is expecting a quantile function"=inherits(fun, c("function", "qf")))
  f <- function(u, .pow=pow, ...){
    .pow <- prm_tr(.pow, TRUE, pfn_pow)
    x <- fun(u,...)
    ((1+x)^(1/.pow)-1)/
      ((1+x)^(1/.pow)+1)
  }

  math_y <- math(fun)
  math_x <- paste0(r"--(\frac{ \left(1+& \right)^)--", prm(nm_pow, TRUE, pfn_pow),
              r"--(-1}{ \left( 1+& \right)^)--", prm(nm_pow, TRUE, pfn_pow), r"--( +1})--")
  math_f <-  fn_insert(math_x, math_y)

  formals_ <- formals(f)
  body_ <- body(f)
  names(formals_)[names(formals_) == ".pow"] <- nm_pow
  body_ <- do.call(substitute, list(body_, list(.pow = as.symbol(nm_pow))))
  as_qf(as.function(c(formals_, body_)), returns="U", math=math_f)
}

#' @param tail numeric. Fixed value for the tail parameter. Default is 1
#' @param nm_tail character.  The name of the tail thickness parameter. The default name is `.dlt`.
#' The tail thickness parameter should be positive (default value is 1).
#' Should be a valid unique variable name other than "u"
#' The asymmetry parameter can be positive or negative (default value is 0).
#' @param asymm numeric. Default value for
#' @param nm_asymm character.  The name of the asymmetry parameter. The default name is `.eps`.
#' Should be a valid unique variable name other than "u"
#' @param pfn_tail,pfn_asymm functions. Parameter transforming functions for tail
#' and asymmetry parameters, respectively. Default is none.
#' @rdname qtransformations
#' @export
#' @examples
#' qtr_shash(qnorm)
qtr_shash <- function(fun, nm_tail=".dlt", nm_asymm=".eps", tail=1, asymm=0, pfn_tail=NULL, pfn_asymm=NULL){
  stopifnot("qtr_shash() is expecting a quantile function"=inherits(fun, c("function", "qf")))
  f <- function(u, .eps=asymm, .dlt=tail, ...){
    .dlt <- prm_tr(.dlt, TRUE, pfn_tail)
    .eps <- prm_tr(.eps, FALSE, pfn_asymm)
    sinh(.dlt*(fun(u,...)- .eps))
  }

  math_y <- math(fun)
  math_x <- paste0(r"--(\text{sinh}\left( )--", prm(nm_tail, TRUE, pfn_tail),
              r"--( \left(& - )--", prm(nm_asymm, FALSE, pfn_asymm),
              r"--( \right)\right)  )--")
  math_f <-  fn_insert(math_x, math_y)

  formals_ <- formals(f)
  body_ <- body(f)
  names(formals_)[names(formals_) == ".dlt"] <- nm_tail
  names(formals_)[names(formals_) == ".eps"] <- nm_asymm
  body_ <- do.call(substitute, list(body_, list(.dlt = as.symbol(nm_tail))))
  body_ <- do.call(substitute, list(body_, list(.eps = as.symbol(nm_asymm))))
  as_qf(as.function(c(formals_, body_)), math=math_f)
}

#' @param base numeric. Fixed value of the base parameter. The default value is 1 (Eulers constant).
#' @param nm_base character. The name of the base parameter.
#' Should be a valid unique variable name other than "u"
#' @param pfn_base function. Parameter transforming function for base argument. Default is exp()
#' @rdname qtransformations
#' @export
#' @examples
#' qf_norm <- qtr_decorate(qnorm, nm_location="mu", nm_scale="sigma")
#' qf_lognorm <- qtr_exponentiate(qf_norm)
#' qf_lognorm(0.2, mu=2, sigma=0.1)
#' qlnorm(0.2, 2, 0.1)
qtr_exponentiate <- function(fun, nm_base=".base", base=1, .invert=FALSE, pfn_base=exp){
  stopifnot("qtr_exponentiate() is expecting a quantile function"=inherits(fun, c("function", "qf")))
  f <- function(u, .base=base, ...){
    .base <- prm_tr(.base, .invert, pfn_base)
    (.base)^fun(u,...)
  }

  math_y <- math(fun)
  math_x <- paste0(prm(nm_base, .invert, pfn_base), r"--(^& )--")
  math_f <-  fn_insert(math_x, math_y)

  formals_ <- formals(f)
  body_ <- body(f)
  names(formals_)[names(formals_) == ".base"] <- nm_base
  body_ <- do.call(substitute, list(body_, list(.base = as.symbol(nm_base))))
  as_qf(as.function(c(formals_, body_)), math=math_f)
}


#' @param .fun function without arguments(or with all default arguments) to be applied as Q-transformation
#' @rdname qtransformations
#' @export
#' @examples
#' qtr_fun(s_exp, log1p)
qtr_fun <- function(fun, .fun){
  stopifnot("qtr_fun() is expecting a quantile function"=inherits(fun, c("function", "qf")))
  f <- function(u, ...)
    .fun(fun(u,...))

  math_y <- math(fun)
  math_x <- paste0(r"--(\text{)--", deparse(substitute(.fun)), "}(&)")
  math_f <- fn_insert(math_x, math_y)

  as_qf(f, math=math_f)
}

#' @param .qf quantile function made with gilchrist (or regular function wrapped to safely accept optional arguments through ellipsis) to be applied as Q-transformation
#' @rdname qtransformations
#' @export
qtr_qf <- function(fun, .qf){
  stopifnot("qtr_qf() is expecting a quantile function"=inherits(fun, c("function", "qf")))
  stopifnot(".qf in qtr_qf() should be a quantile function"=inherits(.qf, c("function", "qf")))
  f <- function(u, ...)
    .qf(fun(u,...), ...)

  math_y <- math(fun)
  math_x <- math(.qf)
  math_f <- fn_insert(math_x, math_y)

  as_qf(f, math=math_f)
}

#' @rdname qtransformations
#' @export
qtr_oddITL <- function(fun){
  stopifnot("qtr_oddITL() is expecting a quantile function"=inherits(fun, c("function", "qf")))
  f <- function(u, ...){
    sqrt(fun(u, ...)) / (1-sqrt(fun(u, ...)))
  }

  math_y <- math(fun)
  math_x <- paste0(r"--(\frac{ \sqrt{&} }{ 1 - \sqrt{&} })--")
  math_f <-  fn_insert(math_x, math_y)

  as_qf(f, math = math_f)
}
