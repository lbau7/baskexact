
<!-- README.md is generated from README.Rmd. -->

# baskexact

<!-- badges: start -->

[![R-CMD-check](https://github.com/lbau7/baskexact/workflows/R-CMD-check/badge.svg)](https://github.com/lbau7/baskexact/actions)
[![Codecov test
coverage](https://codecov.io/gh/lbau7/baskexact/branch/main/graph/badge.svg)](https://codecov.io/gh/lbau7/baskexact?branch=main)
<!-- badges: end -->

## Overview

baskexact calculates the exact operating characteristics of a basket
trial using the design of Fujikawa et al. (
<https://doi.org/10.1002/bimj.201800404>) and related designs.

## Installation

Install the currenct CRAN version of baskexact:

``` r
install.packages("baskexact")
```

Or install the development version from GitHub:

``` r
# install.packages("devtools")
devtools::install_github("https://github.com/lbau7/baskexact")
```

## Usage

With baskexact you can calculate the exact operating characteristics
(type 1 error rate, power, expected number of correct decisions and
expected sample size) of single-stage and two-stage basket trials using
the design of Fujikawa et al.

At first, a design-object has to be created using either
`setupOneStageBasket` for a single-stage trial or `setupTwoStageBasket`
for a two-stage trial. For example:

``` r
library(blindrecalc) # the development version is used for the example
design <- setupOneStageBasket(k = 3, shape1 = 1, shape2 = 1, theta0 = 0.2)
```

`k` is the number of baskets, `shape1` and `shape2` are the two shape
parameters of the Beta-prior of the response probabilities of each
basket and `theta0` is the response probability under the null
hypothesis. Note that currently only common prior parameters and a
common null response probability are supported.

Use `toer` to calculate the type 1 error rate of a certain design:

``` r
toer(
  design = design,
  theta1 = NULL,
  n = 15,
  lambda = 0.99,
  weight_fun = weights_fujikawa,
  weight_params = list(epsilon = 1, tau = 0, logbase = 2),
  results = "group"
)

# $rejection_probabilities
# [1] 0.02533905 0.02533905 0.02533905
# 
# $fwer
# [1] 0.03845609
```

`theta1` refers to the true response probabilities under which the type
1 error rate is computed. Since `theta1 = NULL` is specified, the type 1
error rates under a global null hypothesis are calculated. `n` specifies
the sample size per basket. Currently only equal sample sizes are
supported. `lambda` is the probability cut-off for the posterior
probability cut-off to reject the null hypothesis. If the posterior
probability that the response probability of the basket is larger than
`theta0` is larger than `lambda`, then the null hypothesis is rejected.
`weight_fun` specifies which method should be used to calculate the
weights. With `weights_fujikawa` the weights are calculated based on the
Jensen-Shannon divergence between the individual posterior distributions
for the response probabilities. In `weight_params` a list of parameters
that further define the weights is given. See Fujikawa et al. (2020) for
details. `results` specifies whether only the family wise type 1 error
rate (option `fwer`) or also the basketwise type 1 error rates (option
`group`) are calculated.

To find the probability cut-off `lambda` such that a certain FWER is
maintained, use `adjust_lambda`, for example to find `lambda` such that
the FWER does not exceed 2.5% (note that all hypothesis are tested
one-sided):

``` r
adjust_lambda(
  design = design,
  alpha = 0.025,
  theta1 = NULL,
  n = 15,
  weight_fun = weights_fujikawa,
  weight_params = list(epsilon = 1, tau = 0, logbase = 2),
  prec_digits = 4
)

# $lambda
# [1] 0.9942
# 
# $toer
# [1] 0.02451947
```

`prec_digits` you can specify how many decimal places of `lambda` are
considered. Use `toer` with `lambda = 0.9941` to check that 0.9942 is
indeed the smallest probability cut-off with four decimals with a FWER
of at most 2.5%. Note that even when considering more decimal places the
actual FWER will generally below the nominal level (quite substantially
in some cases), since the outcome (number of responses) is discrete.

Use `pow` to calculate the power of the design:

``` r
pow(
  design = design,
  theta1 = c(0.5, 0.5, 0.5),
  n = 15,
  lambda = 0.9942,
  weight_fun = weights_fujikawa,
  weight_params = list(epsilon = 1, tau = 0, logbase = 2),
  results = "group"
)

# $rejection_probabilities
# [1] 0.9532226 0.9532226 0.9532226
# 
# $ewp
# [1] 0.9871411
```

`pow` has the same parameters as `toer`.
