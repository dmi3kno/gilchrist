#' Basic quantile functions
#'
#' @param u depth
#' @param ... not used
#'
#' @return modified function
#' @rdname basicqf
#' @export
#' @examples
#' sqf_exp(0.6)
#' qexp(0.6, 1)
sqf_exp <- function(u, ...){
  -log(1-u)
}
#' @rdname basicqf
#' @export
sqf_unif <- function(u, ...){
  u
}
