% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/transformations.R
\name{qtr_power}
\alias{qtr_power}
\alias{qtr_epsilon}
\alias{qtr_exponentiate}
\alias{qtr_fun}
\alias{qtr_shiftby}
\alias{qtr_scaleby}
\alias{qtr_powerby}
\title{Q-transformations}
\usage{
qtr_power(fun, nm_pow = ".pow", .invert = FALSE)

qtr_epsilon(fun, nm_pow = ".pow")

qtr_exponentiate(fun, nm_base = ".base", .invert = FALSE)

qtr_fun(fun, .fun)

qtr_shiftby(fun, x = 0)

qtr_scaleby(fun, x = 1)

qtr_powerby(fun, x = 1)
}
\arguments{
\item{fun}{function}

\item{nm_pow}{character.  The name of the power parameter. The default name is \code{.pow}. The default value is 1
Should be a valid unique variable name other than "u"}

\item{.invert}{logical. Should the power parameter be inverted (1/.pow) before applying. Default FALSE}

\item{nm_base}{character.  The name of the base parameter. The default name is \code{.base}. The default value is \code{exp(1)} (Euler's constant).}

\item{.fun}{function without arguments(or with all default arguments) to be applied as Q-transformation}

\item{x}{numeric. Fixed value to shift/scale/power the QF by}
}
\value{
modified function
}
\description{
Some of the typical transformations of QFs, implementing a Q-transformation rule.
\code{qtr_power()}: Raising of QF to a power. Returns \eqn{Q_1(u)^k}.
\code{qtr_exponentiate()}: Exponentiating the QF. Returns \eqn{k^Q_1(u)}.
\code{qtr_fun()}: Q-transform with generic function without additional arguments. \eqn{.fun(Q_1(u))}.
\code{qtr_epsilon()}: unit-Q-transform using inverse epsilon function \eqn{\frac{(1+Q_1(u))^{1/\beta}-1}{(1+Q_1(u))^{1/\beta}+1}}.

Note that today p-transformations can be performed by applying Q-transformations to standard uniform distribution
}
\examples{
qf_exp <- function(u)-log(1-u)
qf_weibull <- qtr_power(qf_exp, "k")
qf_weibull(0.5,k=1/5)
qweibull(0.5, shape = 5)
qf_norm <- qff_decorate(qnorm, nm_location="mu", nm_scale="sigma")
qf_lognorm <- qtr_exponentiate(qf_norm)
qf_lognorm(0.2, mu=2, sigma=0.1)
qlnorm(0.2, 2, 0.1)
qtr_fun(sqf_exp,log1p)
}
