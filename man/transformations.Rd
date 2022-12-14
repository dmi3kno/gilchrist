% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/transformations.R
\name{qtr_power}
\alias{qtr_power}
\alias{qtr_exponentiate}
\title{Q-transformations}
\usage{
qtr_power(fun)

qtr_exponentiate(fun)
}
\arguments{
\item{fun}{function}
}
\value{
modified function
}
\description{
Q-transformations
}
\examples{
qf_exp <- function(u)-log(1-u)
qf_weibull <- qtr_power(qf_exp)
qf_weibull(0.5,pow=1/5)
qweibull(0.5, shape = 5)
}
