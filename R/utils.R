#' Wrapper for producing classed functions.
#'
#' Returns function of classes `c("function", "qf")`
#' Subclass can be added through optional `subclass` argument
#' @param fun function
#' @param subclass additional classes to be added to QF
#' @param expects domain flag for incoming function domain
#' @param returns domain flag for resulting function domain
#' @param math raw string for math describing quantile function
#' @param value string. Value of the attribute
#' @export
#' @rdname qfclass
as_qf <- function(fun, subclass=NULL, expects=NA, returns=NA, math=NA){
  class(fun) <- unique(c("qf", subclass, "function"))
  attr(fun, "expects") <- expects
  attr(fun, "returns") <- returns
  attr(fun, "math") <- math
  fun
}

#' @param x parameter value
#' @param .invert logical. Flag for whether parameter should be inverted
#' @param .fun function. Reference to parameter transforming function
#' @keywords internal
prm_tr <- function(x, .invert, .fun){
  if(!is.null(.fun)){
    stopifnot("Expecting a single-argument function as a parameter transform"=inherits(.fun, "function"))
    x=.fun(x)
  }
  if(.invert) return(1/x)
  x
}

#' Accessors and replacers for function attributes
#' @export
#' @rdname qfclass
expects <- function(fun){
  attr(fun, "expects")
}

#' @export
#' @rdname qfclass
`expects<-` <- function(fun, value){
  attr(fun, "expects") <- value
  fun
}

#' @export
#' @rdname qfclass
returns <- function(fun){
  attr(fun, "returns")
}

#' @export
#' @rdname qfclass
`returns<-` <- function(fun, value){
  attr(fun, "returns") <- value
  fun
}


#' @export
#' @rdname qfclass
math <- function(fun){
  attr(fun, "math")
}

#' @export
#' @rdname qfclass
`math<-` <- function(fun, value){
  attr(fun, "math") <- value
  fun
}
