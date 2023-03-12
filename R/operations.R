#' Gilchrist rules
#'
#' @description
#' Gilchrist rules for constructing valid quantile functions implemented as function factories. All functions take \eqn{Q_1(u)} (the QF of \eqn{X}). Some also take \eqn{Q_2(u)} (the QF of \eqn{Y})
#' `qff_reflect()`: Reflection rule. Returns \eqn{-Q_1(1-u)} (the QF of \eqn{-X}).
#' `qff_reciprocate()`: Reciprocal rule. Returns \eqn{1/Q_1(1-u)} (the QF of \eqn{1/X}).
#' `qff_add()`: Addition rule. Takes \eqn{Q_1(u), Q_2(u)}. Returns \eqn{Q_1(u)+Q_2(u)} - a valid QF of the sum of \eqn{X} and \eqn{Y}.
#' `qff_mix()`: Linear combination rule with weight parameter `.wt`. Takes \eqn{Q_1(u), Q_2(u)}. Returns \eqn{aQ_1(u)+(1-a)Q_2(u)} - a weighted sum QF of \eqn{X} and \eqn{Y}.
#' `qff_cmix()`: Complimentary linear combination rule with weight parameter `.wt`. Takes \eqn{Q_1(u), Q_2(u)}. Returns \eqn{(1-a)Q_1(u)+aQ_2(u)} - a weighted sum QF of \eqn{X} and \eqn{Y}.
#' `qff_multiply()`:Multiplication rule for positive QFs.  Takes \eqn{Q_1(u), Q_2(u)>0} on all \eqn{[0,1]}. Returns \eqn{Q_1(u)Q_2(u)} - a product QF of \eqn{X} and \eqn{Y}.
#' `qff_shift()`: Addition rule with a constant shift (adding location parameter). Takes \eqn{Q_1(u)}. Returns \eqn{a+Q_1(u)}.
#' `qff_scale()`: Multiplication rule with a constant scale (multiplying by the scale parameter). Takes \eqn{Q_1(u)}. Returns \eqn{sQ_1(u)}.
#' @param fun,fun1,fun2 functions
#'
#' @return modified function
#' @rdname rules
#' @export
#'
#' @examples
#' qf_exp <- function(u)-log(1-u)
#' qf_logistic <- qff_add(qf_exp, qff_reflect(qf_exp))
#' qf_logistic(0.6)
#' qlogis(0.6)
qff_reflect <- function(fun){
  function(u, ...){
    -fun(1-u, ...)
  }
}

#' @rdname rules
#' @export
qff_reciprocate <- function(fun){
  function(u,...){
    1/fun(1-u,...)
  }
}

#' @rdname rules
#' @export
qff_add <- function(fun1, fun2){
  function(u, ...){
    fun1(u, ...) + fun2(u, ...)
  }
}

#' @param nm_wt character.  The name of the weight parameter (for mixing). The default name is `.wt`. The default value is 0.5
#' @rdname rules
#' @export
qff_mix <- function(fun1, fun2, nm_wt=".wt"){
   f <- function(u, .wt=0.5, ...){
    (.wt)*fun1(u, ...) + (1-.wt)*fun2(u, ...)
   }

  formals_ <- formals(f)
  body_ <- body(f)
  names(formals_)[names(formals_) == ".wt"] <- nm_wt
  body_ <- do.call(substitute, list(body_, list(.wt = as.symbol(nm_wt))))
  as.function(c(formals_, body_))
}

#' @rdname rules
#' @export
qff_cmix <- function(fun1, fun2, nm_wt=".wt"){
  f <- function(u, .wt=0.5, ...){
    (1-.wt)*fun1(u, ...) + (.wt)*fun2(u, ...)
  }

  formals_ <- formals(f)
  body_ <- body(f)
  names(formals_)[names(formals_) == ".wt"] <- nm_wt
  body_ <- do.call(substitute, list(body_, list(.wt = as.symbol(nm_wt))))
  as.function(c(formals_, body_))
}

#' @rdname rules
#' @export
qff_multiply <- function(fun1, fun2){
  function(u,...){
    f1 <- fun1(u,...)
    f2 <- fun2(u,...)
    stopifnot(all(f1>=0), all(f2>=0))
    #if f1 and f2 are positive for all u \in [0,1]
    f1*f2
  }
}

#' @param nm_shift character.  The name of the shift parameter. The default name is `.location`. The default value is 0
#' @rdname rules
#' @export
qff_shift <- function(fun, nm_shift=".location"){
  f <- function(u, .location=0, ...){
    (.location)+fun(u, ...)
  }

  formals_ <- formals(f)
  body_ <- body(f)
  names(formals_)[names(formals_) == ".location"] <- nm_shift
  body_ <- do.call(substitute, list(body_, list(.location = as.symbol(nm_shift))))
  as.function(c(formals_, body_))
}

#' @param nm_scale character.  The name of the scale parameter. The default name is `.scale`. The default value is 1
#' @param .invert logical. Should the scale parameter be inverted (1/.scale) before applying. Default FALSE
#' @rdname rules
#' @export
qff_scale <- function(fun, nm_scale=".scale", .invert=FALSE){
  f <- function(u, .scale=1, ...){
    if(.invert) .scale <- 1/.scale
    (.scale)*fun(u, ...)
  }

  formals_ <- formals(f)
  body_ <- body(f)
  names(formals_)[names(formals_) == ".scale"] <- nm_scale
  body_ <- do.call(substitute, list(body_, list(.scale = as.symbol(nm_scale))))
  as.function(c(formals_, body_))
}

#' Decorate the basic QF with location and scale parameters
#'
#' @description Add location (`.location`) and scale (`.scale`) parameters (or the inverse scale for `qff_idecorate()`)
#' Note! The parameter names can be changed by passing the names
#'
#' @param fun function to be decorated
#' @param nm_location character. The name of the location parameter. Default name `.location`. Default value is 0.
#' @param nm_scale character. The name of the scale parameter. Default name `.scale`. Default value is 1.
#' @param .invert logical. Should the scale parameter be inverted (1/.scale) before applying. Default FALSE
#'
#' @return modified function with location and scale parameter. The parameter names specified by the user.
#'
#' @examples
#' qff_decorate(s_exp)

#' @rdname decorate
#' @export
qff_decorate <- function(fun, nm_location=".location", nm_scale=".scale", .invert=FALSE){
  function(u, ...){
    s_fun <- qff_scale(fun, nm_scale, .invert=.invert)
    ss_fun <- qff_shift(s_fun, nm_location)
    ss_fun(u,...)
  }
}

