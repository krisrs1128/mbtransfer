
# Microbiome Intervention Analysis using `mbtransfer`

<!-- badges: start -->

[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/krisrs1128/mbtransfer_demo/HEAD?urlpath=rstudio)
<!-- badges: end -->

`mbtransfer` is an R package for modeling and inference of microbial
communities under dynamic environmental interventions. It supports
simulation of hypothetical community trajectories under user-specified
perturbations and can be used to select for taxa that are significantly
impacted by this change in a statistically rigorous sense. It
accomplishes this by learning a separate transfer function model for
each taxon, relating its changes in abundance to its immediate community
and environmental context. The fitted models can be used to measure the
significance of the intervention on individual taxa,

<center>
<img src=""/>
</center>

(larger mirror statistics correspond to stronger effects) and to
simulate trajectory differences between counterfactual treatments,

<center>
<img src=""/>
</center>

You can read more about the methodology in our preprint and can browse
complete examples in the [package
documentation](https://krisrs1128.github.io/mbtransfer).

## Installation

You can install `mbtransfer` from this repository using:

``` r
#install.packages("devtools")
devtools::install_github("krisrs1128/mbtransfer")
```

## Help

We welcome questions and comments about the package! You can reach us
either through [github](https://github.com/krisrs1128/mbtransfer/issues)
or [email](https://krisrs1128.github.io/LSLab/_includes/contact).
