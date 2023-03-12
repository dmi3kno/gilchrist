
<!-- README.md is generated from README.Rmd. Please edit that file -->

# gilchrist <a href='https://dmi3kno.github.io/gilchrist'><img src='man/figures/logo.png' align="right" height="200" /></a>

<!-- badges: start -->

[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
<!-- badges: end -->

This package is made to honor the legacy and memory of Warren G.
Gilchrist (1932-2015)[^1].

The goal of `{gilchrist}` is to implement Gilchrist QF transformation
rules (Gilchrist 2000) in R in the form of pipeable function factories.

## Installation

You can install the development version of gilchrist like so:

``` r
remotes::install_packages("dmi3kno/gilchrist")
```

## Gilchrist’s QF transformation rules

Gilchrist (2000) list the following rules for creating new quantile
functions out of existing ones.

| Original QF       | Rule                    | Resulting QF      | Resulting variable                             |
|-------------------|-------------------------|-------------------|------------------------------------------------|
| $Q_Y(u)$          | Relection rule          | $-Q(1-u)$         | QF of $-Y$                                     |
| $Q_Y(u)$          | Reciprocal rule         | $1/Q(1-u)$        | QF of $1/Y$                                    |
| $Q_1(u),Q_2(u)$   | Addition rule           | $Q_1(u)+Q_2(u)$   | Valid QF                                       |
| $Q_1(u),Q_2(u)$   | Linear combination rule | $aQ_1(u)+bQ_2(u)$ | Valid QF for a,b\>0                            |
| $Q_1(u),Q_2(u)>0$ | Multiplication rule     | $Q_1(u)Q_2(u)$    | Valid QF if $Q_1(u),Q_2(u)>0$                  |
| $Q_Y(u)$          | Q-transformation        | $T(Q_Y(u))$       | QF of $T(Y)$, for non-decreasing $T$           |
| $Q_Y(u)$          | p-transformation        | $Q_Y(H(u))$       | p-transformed $Q_Y(u)$, for non-decreasing $H$ |

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

`{gilchrist}` is a special package! It uses magrittr pipe to operate not
on data, but on functions. It is a function factory engine!

### Exponential

Start with standard exponential quantile function. Note that this “basic
quantile function” has no parameters.

$$S(u)=-\ln(1-u)$$

The built-in equivalent in `{gilchrist}` is `sqf_exp()`. It is a regular
R function, so we can inspect it.

``` r
s_exp
#> function(u, ...){
#>   -log(1-u)
#> }
#> <bytecode: 0x55e891d7ccf8>
#> <environment: namespace:gilchrist>
```

We will now add a scale parameter to our basic exponential QF to make it
like in
[Wikipedia](https://en.wikipedia.org/wiki/Exponential_distribution).

$$Q(u)=\frac{1}{\lambda}[-\ln(1-u)]$$

In order to remember that our scale parameter should be reciprocated we
call it “ilambda”.

``` r
qf_exp <- s_exp %>% 
  qff_scale(nm_scale="ilambda")
```

We compare our hand-made exponential quantile function to the standard
function in R. Note that `qexp` has a reciprocated scale.

``` r
qf_exp(p_grd, ilambda=10)
#> [1]  0.000000  2.231436  5.108256  9.162907 16.094379       Inf
# compare to standard exponential quantile function. 
qexp(p_grd, 1/10)
#> [1]  0.000000  2.231436  5.108256  9.162907 16.094379       Inf
```

`{gilchrist}` has several basic (parameterless) functions that you can
modify.

- `sqf_exp()`: Basic QF of exponential distribution
- `sqf_unif()`: Basic QF of uniform distribution
- `sqf_norm()`: Basic QF of normal distribution, a thinly wrapped
  `qnorm(u,0,1)`.
- `sqf_cauchy()`: Basic QF of Cauchy distribution.
- `sqf_halftriang()`⁠: Basic QF of half-triangular distribution.
- `sqf_halfcosine()`: Basic QF of half-cosine distribution
- `sqf_sech()`: Basic QF of hyperbolic secant distribution.

### Logistic

Let’s do a more challenging example. We will make a logistic
distribution in `{gilchrist}`. Logistic distribution consists of
exponential $-\ln(1-u)$ and reflected exponential $\ln(u)$
distributions.

$$Q(u)=\mu+s\ln\left(\frac{u}{1-u}\right)=\mu+s\left[\ln(u)-\ln(1-u)\right]$$
This is how we do it.

``` r
qf_logistic <- s_exp %>% 
  qff_add(
    s_exp %>% qff_reflect()
  ) %>% 
  qff_decorate("mu", "s")

qf_logistic(p_grd, mu=4, s=2)
#> [1]     -Inf 1.227411 3.189070 4.810930 6.772589      Inf
qlogis(p_grd, 4, 2)
#> [1]     -Inf 1.227411 3.189070 4.810930 6.772589      Inf
```

### Flattened Skew-Logistic

Can we add a little flatness to our newly made logistic distribution and
introduce the weights by the exponential components? Lets make Flattened
Skew-Logistic Distribution described in Sharma and Chakrabarty (2020).

$$Q(u)=\alpha+\beta[(1-\delta)\ln(u)-\delta\ln(1-u)+ku]$$

Note that the exponential distribution will gain a weight `delta` and
the reflected exponential will gain a weight `1-delta`, because this is
the order in which they are listed in `qff_mix`.

``` r
qf_fsld <- s_exp %>% 
  qff_mix(
    qff_reflect(s_exp),
    nm_wt="delta") %>% 
  qff_add(qff_scale(s_unif,"k")) %>% 
  qff_decorate(nm_location="alpha", nm_scale="beta")

qf_fsld(p_grd, delta=0.21, alpha=4, beta=2, k=1)
#> [1]     -Inf 1.950808 3.566807 4.777738 5.923397      Inf
qpd::qfsld(p_grd, bt=2, k=1, dlt=0.21, a=4)
#> [1]     -Inf 1.950808 3.566807 4.777738 5.923397      Inf
```

### Weibull

Last one for this short tutorial. We make Weibull distribution. Weibull
distribution is a Q-transformed exponential distribution. The
transformation function is the the exponent $T(x)=x^k$.

$$Q(u)=\lambda[-\ln(1-u)]^{1/k}$$ Again, to remember that the power
should be reciprocated, let’s call it “ik”.

``` r
qf_weibull <- s_exp %>% 
  qtr_power("ik") %>% 
  qff_scale("lambda")
qf_weibull(p_grd, lambda=2, ik=1/4)
#> [1] 0.000000 1.374599 1.690823 1.956763 2.252675      Inf
qweibull(p_grd, scale=2,  shape=4)
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

[^1]: See Davies (2016) for a short biography of this truly remarkable
    individual
