#' Gilchrist rules
#'
#' @description
#' Gilchrist rules for constructing valid quantile functions implemented as function factories. All functions take \eqn{Q_1(u)} (the QF of \eqn{X}). Some also take \eqn{Q_2(u)} (the QF of \eqn{Y})
#' `qtr_reflect()`: Reflection rule. Returns \eqn{-Q_1(1-u)} (the QF of \eqn{-X}).
#' `qtr_reciprocate()`: Reciprocal rule. Returns \eqn{1/Q_1(1-u)} (the QF of \eqn{1/X}).
#' `qtr_add()`: Addition rule. Takes \eqn{Q_1(u), Q_2(u)}. Returns \eqn{Q_1(u)+Q_2(u)} - a valid QF of the sum of \eqn{X} and \eqn{Y}.
#' `qtr_mix()`: Linear combination rule with weight parameter `.wt`. Takes \eqn{Q_1(u), Q_2(u)}. Returns \eqn{aQ_1(u)+(1-a)Q_2(u)} - a weighted sum QF of \eqn{X} and \eqn{Y}.
#' `qtr_cmix()`: Complimentary linear combination rule with weight parameter `.wt`. Takes \eqn{Q_1(u), Q_2(u)}. Returns \eqn{(1-a)Q_1(u)+aQ_2(u)} - a weighted sum QF of \eqn{X} and \eqn{Y}.
#' `qtr_multiply()`:Multiplication rule for positive QFs.  Takes \eqn{Q_1(u), Q_2(u)>0} on all \eqn{[0,1]}. Returns \eqn{Q_1(u)Q_2(u)} - a product QF of \eqn{X} and \eqn{Y}.
#' `qtr_shift()`: Addition rule with a constant shift (adding location parameter). Takes \eqn{Q_1(u)}. Returns \eqn{a+Q_1(u)}.
#' `qtr_scale()`: Multiplication rule with a constant scale (multiplying by the scale parameter). Takes \eqn{Q_1(u)}. Returns \eqn{sQ_1(u)}.
#' @param fun,fun1,fun2 functions
#'
#' @return modified function
#' @rdname rules
#' @export
#'
#' @examples
#' qf_exp <- function(u)-log(1-u)
#' qf_logistic <- qtr_add(qf_exp, qtr_reflect(qf_exp))
#' qf_logistic(0.6)
#' qlogis(0.6)
qtr_reflect <- function(fun){
  f <- function(u, ...){
    -fun(1-u, ...)
  }

  math_y <- math(fun)
  math_x <- paste0(r"--(-&)--")
  math_p <- paste0("1-&")
  math_f <- fn_insert(fn_insert(math_x, math_y), math_p)

  as_qf(f, math = math_f)
}

#' @rdname rules
#' @export
qtr_reciprocate <- function(fun){
  f <- function(u,...){
    1/fun(1-u,...)
  }

  math_y <- math(fun)
  math_x <- paste0(r"--(\frac{1}{&})--")
  math_p <- paste0("1-&")
  math_f <- fn_insert(fn_insert(math_x, math_y), math_p)

  as_qf(f, math = math_f)
}

#' @rdname rules
#' @export
qtr_add <- function(fun1, fun2){
  f <- function(u, ...){
    fun1(u, ...) + fun2(u, ...)
  }
  math_y1 <- math(fun1)
  math_y2 <- math(fun2)
  math_f <- paste0(br(math_y1, left=r"--{\langle}--", right=r"--{\rangle}--"), "+", 
                    br(math_y2, left=r"--{\langle}--", right=r"--{\rangle}--"))
  as_qf(f, math = math_f)
}

#' @param wt numeric. Fixed value of weight parameter (for mixing). The default is 0.5. 
#' @param nm_wt character.  The name of the weight parameter (for mixing). The default name is `.wt`. The default value is 0.5
#' @rdname rules
#' @export
qtr_mix <- function(fun1, fun2, nm_wt=".wt", wt=0.5){
  f <- function(u, .wt=wt, ...){
    (.wt)*fun1(u, ...) + (1-.wt)*fun2(u, ...)
  }

  math_y1 <- math(fun1)
  math_y2 <- math(fun2)
  math_f <- paste0(prmtr(nm_wt, FALSE),  br(math_y1, left=r"--{\langle}--", right=r"--{\rangle}--"), 
        "+(1-",prmtr(nm_wt, FALSE), ")", br(math_y2, left=r"--{\langle}--", right=r"--{\rangle}--"))

  formals_ <- formals(f)
  body_ <- body(f)
  names(formals_)[names(formals_) == ".wt"] <- nm_wt
  body_ <- do.call(substitute, list(body_, list(.wt = as.symbol(nm_wt))))
  as_qf(as.function(c(formals_, body_)), math = math_f)
}

#' @rdname rules
#' @export
qtr_cmix <- function(fun1, fun2, nm_wt=".wt", wt=0.5){
  f <- function(u, .wt=wt, ...){
    (1-.wt)*fun1(u, ...) + (.wt)*fun2(u, ...)
  }

  math_y1 <- math(fun1)
  math_y2 <- math(fun2)
  math_f <- paste0("(1-", prmtr(nm_wt, FALSE), ")", br(math_y1, left=r"--{\langle}--", right=r"--{\rangle}--"), 
                          "+", prmtr(nm_wt, FALSE), br(math_y2, left=r"--{\langle}--", right=r"--{\rangle}--"))

  formals_ <- formals(f)
  body_ <- body(f)
  names(formals_)[names(formals_) == ".wt"] <- nm_wt
  body_ <- do.call(substitute, list(body_, list(.wt = as.symbol(nm_wt))))
  as_qf(as.function(c(formals_, body_)), math = math_f)
}

#' @rdname rules
#' @export
qtr_multiply <- function(fun1, fun2){
  f <- function(u,...){
    f1 <- fun1(u,...)
    f2 <- fun2(u,...)
    stopifnot(all(f1>=0), all(f2>=0))
    #if f1 and f2 are positive for all u \in [0,1]
    f1*f2
  }

  math_y1 <- math(fun1)
  math_y2 <- math(fun2)
  math_f <- paste0(br(math_y1, left=r"--{\langle}--", right=r"--{\rangle}--"), r"--(\times)--", 
                    br(math_y2, left=r"--{\langle }--", right=r"--{\rangle}--"))

  as_qf(f, math = math_f)
}

#' @param shift numeric. Fixed value for `.location`. The default value is 0
#' @param nm_shift character.  The name of the shift parameter. The default name is `.location`. The default value is 0
#' @rdname rules
#' @export
qtr_shift <- function(fun, nm_shift=".location", shift=0){
  f <- function(u, .location=shift, ...){
    (.location)+fun(u, ...)
  }
  math_y <- math(fun)
  math_f <- paste0(prmtr(nm_shift, FALSE), "+", br(math_y))

  formals_ <- formals(f)
  body_ <- body(f)
  names(formals_)[names(formals_) == ".location"] <- nm_shift
  body_ <- do.call(substitute, list(body_, list(.location = as.symbol(nm_shift))))
  as_qf(as.function(c(formals_, body_)), math = math_f)
}

#' @param scale numeric  Fixed value of the scale parameter. The default value is 1
#' @param nm_scale character.  The name of the scale parameter. The default name is `.scale`. The default value is 1
#' @param .invert logical. Should the scale parameter be inverted (1/.scale) before applying. Default FALSE
#' @rdname rules
#' @export
qtr_scale <- function(fun, nm_scale=".scale", scale=1, .invert=FALSE){
  f <- function(u, .scale=scale, ...){
    if(.invert) .scale <- 1/.scale
    (.scale)*fun(u, ...)
  }

  math_y <- math(fun)
  math_f <- paste0(prmtr(nm_scale, .invert), br(math_y))

  formals_ <- formals(f)
  body_ <- body(f)
  names(formals_)[names(formals_) == ".scale"] <- nm_scale
  body_ <- do.call(substitute, list(body_, list(.scale = as.symbol(nm_scale))))
  as_qf(as.function(c(formals_, body_)), math = math_f)
}

#' Decorate the basic QF with location and scale parameters
#'
#' @description Add location (`.location`) and scale (`.scale`) parameters (or the inverse scale for `qtr_idecorate()`)
#' Note! The parameter names can be changed by passing the names
#'
#' @param fun function to be decorated
#' @param nm_location character. The name of the location parameter. Default name `.location`. Default value is 0.
#' @param nm_scale character. The name of the scale parameter. Default name `.scale`. Default value is 1.
#' @param location numeric. Fixed value for location. Default is 0.
#' @param scale numeric. Fixed value for scale. Default is 1.
#' @param .invert logical. Should the scale parameter be inverted (1/.scale) before applying. Default FALSE
#'
#' @return modified function with location and scale parameter. The parameter names specified by the user.
#'
#' @examples
#' qtr_decorate(s_exp)

#' @rdname decorate
#' @export
qtr_decorate <- function(fun, nm_location=".location", nm_scale=".scale", location=0, scale=1, .invert=FALSE){
  f <- function(u, ...){
    s_fun <- qtr_scale(fun, nm_scale=nm_scale, scale=scale, .invert=.invert)
    ss_fun <- qtr_shift(s_fun, nm_shift=nm_location, shift=location)
    ss_fun(u,...)
  }

  math_y <- math(fun)
  math_f <- paste0(prmtr(nm_location, .invert), "+", prmtr(nm_scale, .invert), br(math_y))

  as_qf(f, math=math_f)
}