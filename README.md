
<!-- README.md is generated from README.Rmd. Please edit that file -->

# gilchrist

<!-- badges: start -->

[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
<!-- badges: end -->

This package is made to honor the legacy and memory of Warren G.
Gilchrist (1932-2015). See Davies (2016) for a short biography of this
truly remarkable individual.

The goal of `{gilchrist}` is to implement Gilchrist QF transformation
rules (Gilchrist 2000) in R in the form of pipeable function factories.

## Installation

You can install the development version of gilchrist like so:

``` r
remotes::install_packages("dmi3kno/gilchrist")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(gilchrist)
library(magrittr)
## basic example code
```

In using and testing quantile function it is useful to have an
equi-spaced grid of probabilities

``` r
p_grd <- seq(0,1, by=0.2)
```

Start with standard exponential quantile function

$$Q(u)=-\ln(1-u)$$

and add location and scale parameters. We compare our hand-made
exponential quantile function to the standard function in R. Note that
`qexp` misses (defaults to 0) the location argument and reciprocates the
scale.

``` r
sqf_exp
#> function(u, ...){
#>   -log(1-u)
#> }
#> <bytecode: 0x558d49ac1780>
#> <environment: namespace:gilchrist>

qf_exp <- sqf_exp %>% 
  qff_decorate()

qf_exp(p_grd, location=0, scale=10)
#> [1]  0.000000  2.231436  5.108256  9.162907 16.094379       Inf
# compare to standard exponential quantile function. 
qexp(p_grd, 1/10)
#> [1]  0.000000  2.231436  5.108256  9.162907 16.094379       Inf
```

Let’s do a more challenging example. We will make a logistic
distribution in `{gilchrist}`. Logistic distribution consists of
exponential and reflected exponential distribution. This is how we do
it.

``` r
qf_logistic <- sqf_exp %>% 
  qff_add(
    sqf_exp %>% qff_reflect()
  ) %>% 
  qff_decorate()

qf_logistic(p_grd, location=4, scale=2)
#> [1]     -Inf 1.227411 3.189070 4.810930 6.772589      Inf
qlogis(p_grd, 4, 2)
#> [1]     -Inf 1.227411 3.189070 4.810930 6.772589      Inf
```

Can we add a little flatness to our newly made logistic distribution?
Lets make Flattened Skew-Logistic Distribution described in Sharma and
Chakrabarty (2020).

``` r
qf_fsld <- sqf_exp %>% 
  qff_mix(
    sqf_exp %>% qff_reflect()
  ) %>% 
  qff_add(sqf_unif) %>% 
  qff_decorate()

qf_fsld(p_grd, wt=0.21, location=4, scale=2)
#> [1]     -Inf 1.950808 3.566807 4.777738 5.923397      Inf
qpd::qfsld(p_grd, bt=2, k=1, dlt=0.21, a=4)
#> [1]     -Inf 1.950808 3.566807 4.777738 5.923397      Inf
```

Last one for this short tutorial. We make Weibull distribution. Weibull
distribution is a Q-transformed exponential distribution. The
transformation function is the the exponent $T(x)=x^a$.

``` r
qf_weibull <- sqf_exp %>% 
  qtr_power() %>% 
  qff_decorate()
qf_weibull(p_grd, location=0, scale=2, pow=1/4)
#> [1] 0.000000 1.374599 1.690823 1.956763 2.252675      Inf
qweibull(p_grd, shape=4, scale=2)
#> [1] 0.000000 1.374599 1.690823 1.956763 2.252675      Inf
```

Therefore, you can compose new quantile functions following Gilchrist
transformation rules.

## References

<div id="refs" class="references csl-bib-body hanging-indent">

<div id="ref-davies2016WarrenGilchrist19322015" class="csl-entry">

Davies, Neville. 2016. “Warren g. Gilchrist, 1932-2015.” *Journal of the
Royal Statistical Society. Series A (Statistics in Society)* 179 (3):
872–74. <https://doi.org/10.1111/rssa.12243>.

</div>

<div id="ref-gilchrist2000StatisticalModellingQuantile"
class="csl-entry">

Gilchrist, Warren. 2000. *Statistical Modelling with Quantile
Functions*. Boca Raton: Chapman & Hall/CRC.

</div>

<div id="ref-sharma2020QuantileBasedApproachSupervised"
class="csl-entry">

Sharma, Dreamlee, and Tapan Kumar Chakrabarty. 2020. “A Quantile-Based
Approach to Supervised Learning.” In *Applications of Machine Learning*,
edited by Prashant Johri, Jitendra Kumar Verma, and Sudip Paul, 321–40.
Singapore: Springer Singapore.
<https://doi.org/10.1007/978-981-15-3357-0_21>.

</div>

</div>
