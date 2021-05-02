
<!-- README.md is generated from README.Rmd. Please edit that file -->

# <img src="man/figures/rmm-hexagon.png" align="right" />

### Bayesian multiple membership multilevel regression models with endogenized weights

<!-- badges: start -->

[![R-CMD-check](https://github.com/benrosche/rmm/workflows/R-CMD-check/badge.svg)](https://github.com/benrosche/rmm/actions)
<!-- badges: end -->

The rmm package provides an interface to fit Bayesian multiple
membership multilevel models with endogenized weights using JAGS from
within R for a variety of outcomes (linear, logit, conditional logit,
Cox, Weibull).

Most multilevel analyses examine how lower-level units
(e.g. individuals) are affected by their embedding in
contextual/aggregate units at a higher level (e.g. neighborhoods) (=
macro-micro link). The rmm package uses the multiple membership
multilevel model to conceptually reverse this setup. It allows studying
how the effect of units at lower levels propagates to a higher level (=
micro-macro link).

Previous studies examining micro-macro links either aggregated or
disaggregated the data. Both approaches obstruct the inherent
aggregation problem, ignore dependencies among observations, which
induces excessive Type-I and Type-II error, and cannot separate micro-
from macro-level variance. The generalized multiple membership
multilevel model (MMMM) is able to overcome these problems by explicitly
modeling the aggregation from micro to macro level by including an
aggregation function in the regression model. It is a theoretically and
statistically sound solution to the study of micro-macro links with
regression analysis.

### With the **rmm** package, you can…

-   explicitly model how the effect of lower level units propagate to a
    higher level

-   examine the fit of different aggregation functions, such as the min,
    max, mean, or sum

-   uncover heterogeneity in the effect of lower-level units on
    higher-level entities

-   control dependencies that arise from a multiple membership
    (crisscrossing) data structure

-   estimate the (residual) variance at the lower and the higher level

The data structure that is being modeled looks like this:

…

This is how the `rmm()` function looks like:

``` r
data(coalgov)
rmm(Y ~ 1 + X_l2 + 
      X_l3 + hm(id=id_l3, type=RE) +
      mm(id(id_l2, id_l1), mmc(X_l1), mmw(w ~ 1/offset(n)^exp(-(X_w)))), 
    family="Gaussian", data=coalgov)
```

### Developers

I welcome contributions to the package! Feel free to submit changes for
review or contact me if you have any questions.

### Issues or Feature Requests

If you would like to log an issue or submit a feature request, please
create a new issue or comment on an existing issue on [GitHub
Issues](https://github.com/benrosche/rmm/issues) on this repo.

### Changelog

See [NEWS.md](https://github.com/benrosche/rmm/news/index.html) for the
package changelog.

More information can be found [here](http://benrosche.com/projects/rmm/)
