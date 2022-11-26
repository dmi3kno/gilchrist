#' Q-transformations
#'
#' @param fun function
#'
#' @return modified function
#' @rdname transformations
#' @export
#' @examples
#' qf_exp <- function(u)-log(1-u)
#' qf_weibull <- qtr_power(qf_exp)
#' qf_weibull(0.5,pow=1/5)
#' qweibull(0.5, shape = 5)
qtr_power <- function(fun){
  function(u, pow=1, ...)
    fun(u,...)^pow
}
#' @rdname transformations
#' @export
qtr_exponentiate <- function(fun){
  function(u, base=1, ...)
    base^fun(u,...)
}


