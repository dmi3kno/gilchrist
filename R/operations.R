#' Gilchrist rules
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
qff_add <- function(fun1, fun2){
  function(u, ...){
    fun1(u, ...) + fun2(u, ...)
  }
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
#' @rdname rules
#' @export
qff_mix <- function(fun1, fun2){
  function(u, wt=0.5, ...){
    (wt)*fun1(u, ...) + (1-wt)*fun2(u, ...)
  }
}
#' @rdname rules
#' @export
qff_imix <- function(fun1, fun2){
  function(u, wt=0.5, ...){
    (1-wt)*fun1(u, ...) + (wt)*fun2(u, ...)
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
qff_decorate <- function(fun){
  function(u, location=0, scale=1, ...){
    location+scale*fun(u, ...)
  }
}
