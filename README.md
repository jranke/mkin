# mkin

[![](https://www.r-pkg.org/badges/version/mkin)](https://cran.r-project.org/package=mkin)
[![mkin status badge](https://jranke.r-universe.dev/badges/mkin)](https://jranke.r-universe.dev/ui#package:mkin)
[![Build Status](https://travis-ci.com/jranke/mkin.svg?branch=main)](https://app.travis-ci.com/github/jranke/mkin)
[![codecov](https://codecov.io/github/jranke/mkin/branch/main/graphs/badge.svg)](https://codecov.io/github/jranke/mkin)

The R package **mkin** provides calculation routines for the analysis of
chemical degradation data, including <b>m</b>ulticompartment <b>kin</b>etics as
needed for modelling the formation and decline of transformation products, or
if several degradation compartments are involved.

## Installation

You can install the latest released version from
[CRAN](https://cran.r-project.org/package=mkin) from within R:

```r
install.packages("mkin")
```

## Background

In the regulatory evaluation of chemical substances like plant protection
products (pesticides), biocides and other chemicals, degradation data play an
important role. For the evaluation of pesticide degradation experiments,
detailed guidance and helpful tools have been developed as detailed in
'Credits and historical remarks' below.

## Usage

For a start, have a look at the code examples provided for
[`plot.mkinfit`](https://pkgdown.jrwb.de/mkin/reference/plot.mkinfit.html)
and
[`plot.mmkin`](https://pkgdown.jrwb.de/mkin/reference/plot.mmkin.html), and
at the package vignettes
[`FOCUS L`](https://pkgdown.jrwb.de/mkin/articles/FOCUS_L.html) and
[`FOCUS D`](https://pkgdown.jrwb.de/mkin/articles/FOCUS_D.html).

## Documentation

The HTML documentation of the latest version released to CRAN is available at
[jrwb.de](https://pkgdown.jrwb.de/mkin/) and
[github](https://jranke.github.io/mkin/). Documentation of the development
version is found in the ['dev' subdirectory](https://pkgdown.jrwb.de/mkin/dev/).

## Features

### General

* Highly flexible model specification using
  [`mkinmod`](https://pkgdown.jrwb.de/mkin/reference/mkinmod.html),
  including equilibrium reactions and using the single first-order reversible
  binding (SFORB) model, which will automatically create two state variables
  for the observed variable.
* Model solution (forward modelling) in the function
  [`mkinpredict`](https://pkgdown.jrwb.de/mkin/reference/mkinpredict.html)
  is performed either using the analytical solution for the case of
  parent only degradation or some simple models involving a single transformation
  product, , an eigenvalue based solution if only simple first-order (SFO) or
  SFORB kinetics are used in the model, or using a numeric solver from the
  `deSolve` package (default is `lsoda`).
* The usual one-sided t-test for significant difference from zero is
  shown based on estimators for the untransformed parameters.
* Summary and plotting functions. The `summary` of an `mkinfit` object is in
  fact a full report that should give enough information to be able to
  approximately reproduce the fit with other tools.
* The chi-squared error level as defined in the FOCUS kinetics guidance
  (see below) is calculated for each observed variable.
* The 'variance by variable' error model which is often fitted using
  Iteratively Reweighted Least Squares (IRLS) can be specified as
  `error_model = "obs"`.

### Unique in mkin

* Three different error models can be selected using the argument `error_model`
  to the [`mkinfit`](https://pkgdown.jrwb.de/mkin/reference/mkinfit.html)
  function. A two-component error model similar to the one proposed by
  [Rocke and Lorenzato](https://pkgdown.jrwb.de/mkin/reference/sigma_twocomp.html)
  can be selected using the argument `error_model = "tc"`.
* Model comparisons using the Akaike Information Criterion (AIC) are supported
  which can also be used for non-constant variance. In such cases the FOCUS
  chi-squared error level is not meaningful.
* By default, kinetic rate constants and kinetic formation fractions are
  transformed internally using
  [`transform_odeparms`](https://pkgdown.jrwb.de/mkin/reference/transform_odeparms.html)
  so their estimators can more reasonably be expected to follow
  a normal distribution.
* When parameter estimates are backtransformed to match the model definition,
  confidence intervals calculated from standard errors are also backtransformed
  to the correct scale, and will not include meaningless values like negative
  rate constants or formation fractions adding up to more than 1, which cannot
  occur in a single experiment with a single defined radiolabel position.
* When a metabolite decline phase is not described well by SFO kinetics,
  SFORB kinetics can be used for the metabolite. Mathematically, the SFORB model
  is equivalent to the DFOP model. However, the SFORB model has the advantage
  that there is a mechanistic interpretation of the model parameters.
* Nonlinear mixed-effects models (hierarchical models) can be created from fits
  of the same degradation model to different datasets for the same compound by
  using the
  [nlme.mmkin](https://pkgdown.jrwb.de/mkin/reference/nlme.mmkin.html) and
  [saem.mmkin](https://pkgdown.jrwb.de/mkin/reference/saem.html)
  methods. Note that the convergence of the nlme fits depends on the quality of
  the data. Convergence is better for simple models and data for many groups
  (e.g. soils). The saem method uses the `saemix` package as a backend. Analytical
  solutions suitable for use with this package have been implemented for parent
  only models and the most important models including one metabolite (SFO-SFO
  and DFOP-SFO). Fitting other models with `saem.mmkin`, while it makes use
  of the compiled ODE models that mkin provides, has longer run times (from a couple
  of minutes to more than an hour).

### Performance

* Parallel fitting of several models to several datasets is supported, see for
  example
  [`plot.mmkin`](https://pkgdown.jrwb.de/mkin/reference/plot.mmkin.html).
* If a C compiler is installed, the kinetic models are compiled from automatically
  generated C code, see
  [vignette `compiled_models`](https://pkgdown.jrwb.de/mkin/articles/web_only/compiled_models.html).
  The autogeneration of C code was
  inspired by the [`ccSolve`](https://github.com/karlines/ccSolve) package. Thanks
  to Karline Soetaert for her work on that.
* Even if no compiler is installed, many degradation models still give
  [very good performance](https://pkgdown.jrwb.de/mkin/articles/web_only/benchmarks.html),
  as current versions of mkin also have [analytical solutions for some models
  with one metabolite](https://jrwb.de/performance-improvements-mkin/), and if
  SFO or SFORB are used for the parent compound, Eigenvalue based solutions of
  the degradation model are available.

## GUI

There is a graphical user interface that may be useful. Please
refer to its [documentation page](https://pkgdown.jrwb.de/gmkin/)
for installation instructions and a manual. It only supports
evaluations using (generalised) nonlinear regression, but
not simultaneous fits using nonlinear mixed-effects models.

## News

There is a list of changes for the latest [CRAN release](https://cran.r-project.org/package=mkin/news/news.html)
and one for each github branch, e.g. [the main branch](https://github.com/jranke/mkin/blob/main/NEWS.md).

## Credits and historical remarks

`mkin` would not be possible without the underlying software stack consisting of,
among others, R and the package [deSolve](https://cran.r-project.org/package=deSolve).
In previous version, `mkin` was also using the functionality of the
[FME](https://cran.r-project.org/package=FME) package. Please refer to the
[package page on CRAN](https://cran.r-project.org/package=mkin) for the full list
of imported and suggested R packages. Also, [Debian Linux](https://debian.org),
the vim editor and the [Nvim-R](https://github.com/jalvesaq/Nvim-R) plugin have
been invaluable in its development.

`mkin` could not have been written without me being introduced to regulatory fate
modelling of pesticides by Adrian Gurney during my time at Harlan Laboratories
Ltd (formerly RCC Ltd). `mkin` greatly profits from and largely follows
the work done by the
[FOCUS Degradation Kinetics Workgroup](http://esdac.jrc.ec.europa.eu/projects/degradation-kinetics),
as detailed in their guidance document from 2006, slightly updated in 2011 and
in 2014.

Also, it was inspired by the first version of KinGUI developed by
BayerCropScience, which is based on the MatLab runtime environment.

The companion package
[kinfit](http://kinfit.r-forge.r-project.org/kinfit_static/index.html) (now deprecated) was
[started in 2008](https://r-forge.r-project.org/scm/viewvc.php?view=rev&root=kinfit&revision=2) and
[first published](https://cran.r-project.org/src/contrib/Archive/kinfit/) on
CRAN on 01 May 2010.

The first `mkin` code was
[published on 11 May 2010](https://r-forge.r-project.org/scm/viewvc.php?view=rev&root=kinfit&revision=8) and the
[first CRAN version](https://cran.r-project.org/src/contrib/Archive/mkin/)
on 18 May 2010.

In 2011, Bayer Crop Science started to distribute an R based successor to KinGUI named
KinGUII whose R code is based on `mkin`, but which added, among other
refinements, a closed source graphical user interface (GUI), iteratively
reweighted least squares (IRLS) optimisation of the variance for each of the
observed variables, and Markov Chain Monte Carlo (MCMC) simulation
functionality, similar to what is available e.g. in the `FME` package.

Somewhat in parallel, Syngenta has sponsored the development of an `mkin` and
KinGUII based GUI application called CAKE, which also adds IRLS and MCMC, is
more limited in the model formulation, but puts more weight on usability.
CAKE is available for download from the [CAKE
website](https://cake-kinetics.org), where you can also
find a zip archive of the R scripts derived from `mkin`, published under the GPL
license.

Finally, there is
[KineticEval](https://github.com/zhenglei-gao/KineticEval), which contains
some further development of the scripts used for KinGUII.

Thanks to René Lehmann, formerly working at the Umweltbundesamt, for the nice
cooperation on parameter transformations, especially the isometric log-ratio
transformation that is now used for formation fractions in case there are more
than two transformation targets.

Many inspirations for improvements of mkin resulted from doing kinetic evaluations
of degradation data for my clients while working at Harlan Laboratories and
at Eurofins Regulatory AG, and now as an independent consultant.

Funding was received from the Umweltbundesamt in the course of the projects

- Project Number 27452 (Testing and validation of modelling software as an alternative
to ModelMaker 4.0, 2014-2015)
- Project Number 56703 (Optimization of gmkin for routine use in the Umweltbundesamt, 2015)
- Project Number 92570 (Update of Project Number 27452, 2017-2018)
- Project Number 112407 (Testing the feasibility of using an error model
  according to Rocke and Lorenzato for more realistic parameter estimates in
  the kinetic evaluation of degradation data, 2018-2019)
- Project Number 120667 (Development of objective criteria for the evaluation
  of the visual fit in the kinetic evaluation of degradation data, 2019-2020)
- Project Number 146839 (Checking the feasibility of using mixed-effects models for
  the derivation of kinetic modelling parameters from degradation studies, 2020-2021)
- Project Number 173340 (Application of nonlinear hierarchical models to the
  kinetic evaluation of chemical degradation data)

Thanks to everyone involved for collaboration and support!

Thanks are due also to Emmanuelle Comets, maintainer of the saemix package, for
her interest and support for using the SAEM algorithm and its implementation in
saemix for the evaluation of chemical degradation data.

## References

<table>
  <tr><td>Ranke J, Wöltjen J, Schmidt J, and Comets E (2021)
  Taking kinetic evaluations of degradation data to the next level with nonlinear mixed-effects models.
  <i>Environments</i>
  <b>8</b> (8) 71
  <a href='https://doi.org/10.3390/environments8080071'>doi:10.3390/environments8080071</a>
  </td></tr>

  <tr><td>Ranke J, Meinecke S (2019)
  Error Models for the Kinetic Evaluation of Chemical Degradation Data
  <i>Environments</i>
  <b>6</b> (12) 124
  <a href='https://doi.org/10.3390/environments6120124'>doi:10.3390/environments6120124</a>
  </td></tr>

  <tr><td>Ranke J, Wöltjen J, Meinecke S (2018)
  Comparison of software tools for kinetic evaluation of chemical degradation data
  <i>Environmental Sciences Europe</i>
  <b>30</b> 17
  <a href='https://doi.org/10.1186/s12302-018-0145-1'>doi:10.1186/s12302-018-0145-1</a>
  </td></tr>
</table>

## Development

Contributions are welcome!
