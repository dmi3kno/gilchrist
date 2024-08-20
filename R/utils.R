#' Wrapper for producing classed functions. 
#'
#' Returns function of classes `c("function", "qf")`
#' Subclass can be added through optional `subclass` argument
#' @param fun function
#' @param subclass additional classes to be added to QF
#' @param expects domain flag for incoming function domain
#' @param returns domain flag for resulting function domain
#' @param math raw string for math describing quantile function
#' @export
#' @rdname qfclass
as_qf <- function(fun, subclass=NULL, expects=NA, returns=NA, math=NA){
  class(fun) <- unique(c("qf", subclass, "function"))
  attr(fun, "expects") <- expects
  attr(fun, "returns") <- returns
  attr(fun, "math") <- math
  fun
}

#' Accessors for function attributes
#' @export
#' @rdname qfclass
expects <- function(fun){
  attr(fun, "expects")
}

#' @export
#' @rdname qfclass
returns <- function(fun){
  attr(fun, "returns")
}

#' @export
#' @rdname qfclass
math <- function(fun){
  attr(fun, "math")
}

