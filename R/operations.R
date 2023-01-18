#' Gilchrist rules
#'
#' @description
#' Gilchrist rules for constructing valid quantile functions implemented as function factories. All functions take \eqn{Q_1(u)} (the QF of \eqn{X}). Some also take \eqn{Q_2(u)} (the QF of \eqn{Y})
#' `qff_reflect()`: Reflection rule. Returns \eqn{-Q_1(1-u)} (the QF of \eqn{-X}).
#' `qff_reciprocate()`: Reciprocal rule. Returns \eqn{1/Q_1(1-u)} (the QF of \eqn{1/X}).
#' `qff_add()`: Addition rule. Takes \eqn{Q_1(u), Q_2(u)}. Returns \eqn{Q_1(u)+Q_2(u)} - a valid QF of the sum of \eqn{X} and \eqn{Y}.
#' `qff_mix()`: Linear combination rule with weight parameter `.wt`. Takes \eqn{Q_1(u), Q_2(u)}. Returns \eqn{aQ_1(u)+(1-a)Q_2(u)} - a weighted sum QF of \eqn{X} and \eqn{Y}.
#' `qff_imix()`: Complimentary linear combination rule with weight parameter `.wt`. Takes \eqn{Q_1(u), Q_2(u)}. Returns \eqn{(1-a)Q_1(u)+aQ_2(u)} - a weighted sum QF of \eqn{X} and \eqn{Y}.
#' `qff_multiply()`:Multiplication rule for positive QFs.  Takes \eqn{Q_1(u), Q_2(u)>0} on all \eqn{[0,1]}. Returns \eqn{Q_1(u)Q_2(u)} - a product QF of \eqn{X} and \eqn{Y}.
#'
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
qff_imix <- function(fun1, fun2, nm_wt=".wt"){
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
    stopifnot(all(f1>0), all(f2>0))
    #if f1 and f2 are positive for all u \in [0,1]
    f1*f2
  }
}

#' Decorate the basic QF with location and scale parameters
#'
#' @description Add location (`.location`) and scale (`.scale`) parameters
#' Note! The parameter names can be changed by passing the names
#'
#' @param fun function to be decorated
#' @param nm_location character. The name of the location parameter. Default name `.location`. Default value is 0.
#' @param nm_scale character. The name of the scale parameter. Default name `.scale`. Default value is 1.
#'
#' @return modified function with location and scale parameter. The parameter names specified by the user.
#'
#' @examples
#' qff_decorate(sf_exp)

#' @rdname decorate
#' @export
qff_decorate <- function(fun, nm_location=".location", nm_scale=".scale"){
  f <- function(u, .location=0, .scale=1, ...){
    .location+.scale*fun(u, ...)
  }

  formals_ <- formals(f)
  body_ <- body(f)
  names(formals_)[names(formals_) == ".location"] <- nm_location
  names(formals_)[names(formals_) == ".scale"] <- nm_scale
  body_ <- do.call(substitute, list(body_, list(.location = as.symbol(nm_location),
                                                .scale = as.symbol(nm_scale))))
  as.function(c(formals_, body_))
}
