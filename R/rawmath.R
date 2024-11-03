#' Show raw math
#'
#' @param x character. Variable name.
#' @rdname rawmath
#' @keywords internal
vrbl <- function(x){
  grk_ltr <- c("alpha", "beta", "gamma", "Gamma",
    "delta", "Delta", "epsilon", "varepsilon", "zeta",
    "eta", "theta", "vartheta", "Theta", "iota", "kappa",
    "lambda", "Lambda", "mu", "nu", "xi", "Xi", "pi", "Pi", "varpi",
    "rho", "varrho", "sigma", "varsigma", "Sigma", "tau", "upsilon", "Upsilon",
    "phi", "varphi", "Phi", "chi", "psi", "Psi", "omega", "Omega")
  if(grepl(paste0(grk_ltr, collapse = "|"), x)) return(  paste0(r"--( {\)--", x ,r"--(})--") )
  paste0(r"--(\text{)--", x , r"--(})--")
}

#' @param x raw string. Math string to be bracketed.
#' @param left Left bracket
#' @param right Right bracket
#' @rdname rawmath
#' @keywords internal
br <- function(x, left="[", right="]"){
  paste0(r"--{\left}--",left, x, r"--{\right}--", right)
}

#' @param x raw string. Host math string.
#' @param y raw string. Math string to be inserted.
#' @param br logical. Should the inserted string be bracketed? Default is TRUE.
#' @param left raw string. Left bracket
#' @param right raw string. Right bracket
#' @param pl character. Single character indicating the placeholder to be replaced by inserted string
#' @rdname rawmath
#' @keywords internal
fn_insert <- function(x, y, br=TRUE, left="[", right="]", pl = "&"){
  if(br) y <- br(y, left=left, right=right)
  gsub(pattern=pl, replacement=y, x, fixed=TRUE)
}

#' @param x raw string. Parameter name. Will be passed `vrbl()`
#' @param .invert logical. Should the parameter be inverted (1/x). Default is FALSE
#' @param .fun parameter transforming function
#' @rdname rawmath
#' @keywords internal
prm <- function(x, .invert=FALSE, .fun=NULL){
  v <- vrbl(x)
  if(!is.null(.fun)) v <- paste0(r"--(\text{)--", deparse(substitute(.fun)),r"--(}\left( )--",v,r"--(\right) )--")
  if(.invert) return(paste0(r"--(\frac{1}{)--",v, r"--(})--"))
  paste0(r"--({)--",v, r"--( })--")
}

#' @param x raw string. Math string to be displayed
#' @param prfx raw string. Prefix string. Default is Q(u)=
#' @param pl character. Single character inidicating the placeholder.
#' @param u character. Variable name, which will be placed instead of placeholder.
#' @rdname rawmath
#' @export
#' @keywords internal
display <- function(.qf, prfx="Q(u)=", pl="&", u="u"){
    x <- gsub(pl, u, math(.qf))
    cat(paste0("$$", prfx,  x, "$$"))
}

#' @param x raw string. Math string to be displayed
#' @param prfx raw string. Prefix string. Default is Q(u)=
#' @param pl character. Single character inidicating the placeholder.
#' @param u character. Variable name, which will be placed instead of placeholder.
#' @rdname rawmath
#' @export
#' @keywords internal
inline <- function(.qf, prfx="Q(u)=", pl="&", u="u"){
  x <- gsub(pl, u, math(.qf))
  paste0("$", prfx,  x, "$")
}
