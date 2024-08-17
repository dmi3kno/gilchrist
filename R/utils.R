#' Wrapper for producing classed functions. 
#'
#' Returns function of classes `c("function", "qf")`
#' Subclass can be added through optional `subclass` argument
#' @param fun function
#' @param subclass additional classes to be added to QF
#' @param expects domain flag for incoming function domain
#' @param returns domain flag for resulting function domain
#' @export
#' @rdname qfclass
as_qf <- function(fun, subclass=NULL, expects=NA, returns=NA){
  class(fun) <- unique(c("function", "qf", subclass))
  attr(fun, "expects") <- expects
  attr(fun, "returns") <- returns
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

