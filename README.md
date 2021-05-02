# rmm <img src="man/figures/rmm-hexagon.png" align="right" />

## An R package for fitting Bayesian multiple membership multilevel regression models 

The rmm package provides an interface to fit Bayesian multiple membership multilevel models with endogenized weights using JAGS from within R for a variety of outcomes (linear, logit, conditional logit, Cox, Weibull).

Most multilevel analyses examine how lower-level units (e.g. individuals) are affected by their embedding in contextual/aggregate units at a higher level (e.g. neighborhoods) (= macro-micro link). The rmm package uses the multiple membership multilevel model to conceptually reverse this setup. It allows studying how the effect of units at lower levels propagates to a higher level (= micro-macro link).

Previous studies examining micro-macro links either aggregated or disaggregated the data. Both approaches obstruct the inherent aggregation problem, ignore dependencies among observations, which induces excessive Type-I and Type-II error, and cannot separate micro- from macro-level variance. The generalized multiple membership multilevel model (MMMM) is able to overcome these problems by explicitly modeling the aggregation from micro to macro level by including an aggregation function in the regression model. It is a theoretically and statistically sound solution to the study of micro-macro links with regression analysis.

## With the **rmm** package, you can...

1. explicitly model how the effect of lower level units propagate to a higher level

2. examine the fit of different aggregation functions, such as the min/max, mean, sum function

3. uncover heterogeneity in the effect of lower-level units on higher-level entities

4. control dependencies that arise from crisscrossing data structures

5. estimate the (residual) variance at lower and higher level

Here is an example of **rmm** in action:

<img src="https://raw.githubusercontent.com/microsoft/wpa/main/man/figures/output2.gif" align="center" width=80% />

More information can be found [here](http://benrosche.com/projects/rmm/)

### Developers

I welcome contributions to the package! Feel free to submit changes for review or contact me if you have any questions.

### Issues or Feature Requests
If you would like to log an issue or submit a feature request, please create a new issue or comment on an existing issue on [GitHub Issues](https://github.com/benrosche/rmm/issues) on this repo.

### Changelog
See [NEWS.md](https://github.com/benrosche/rmm/news/index.html) for the package changelog.

