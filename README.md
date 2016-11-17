# mkin

[![](https://www.r-pkg.org/badges/version/mkin)](https://cran.r-project.org/package=mkin)

The R package **mkin** provides calculation routines for the analysis of
chemical degradation data, including <b>m</b>ulticompartment <b>kin</b>etics as
needed for modelling the formation and decline of transformation products, or
if several compartments are involved.

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

For a start, have a look a the code examples provided for
[`plot.mkinfit`](http://kinfit.r-forge.r-project.org/mkin_static/plot.mkinfit.html)
and 
[`plot.mmkin`](http://kinfit.r-forge.r-project.org/mkin_static/plot.mmkin.html), and 
at the package vignettes
[`FOCUS L`](http://kinfit.r-forge.r-project.org/mkin_static/vignettes/FOCUS_L.html) and
[`FOCUS D`](http://kinfit.r-forge.r-project.org/mkin_static/vignettes/FOCUS_D.html).

## Documentation

The [HTML documentation](http://jranke.github.io/mkin) is
maintained at the github project site.

## Features

* Highly flexible model specification using
  [`mkinmod`](http://kinfit.r-forge.r-project.org/mkin_static/mkinmod.html),
  including equilibrium reactions and using the single first-order 
  reversible binding (SFORB) model, which will automatically create
  two latent state variables for the observed variable.
* As of version 0.9-39, fitting of several models to several datasets, optionally in 
  parallel, is supported, see for example
  [`plot.mmkin`](http://kinfit.r-forge.r-project.org/mkin_static/plot.mmkin.html).
* Model solution (forward modelling) in the function
  [`mkinpredict`](http://kinfit.r-forge.r-project.org/mkin_static/mkinpredict.html) 
  is performed either using the analytical solution for the case of 
  parent only degradation, an eigenvalue based solution if only simple
  first-order (SFO) or SFORB kinetics are used in the model, or
  using a numeric solver from the `deSolve` package (default is `lsoda`).
* If a C compiler is installed, the kinetic models are compiled from automatically
  generated C code, see  
  [vignette `compiled_models`](http://kinfit.r-forge.r-project.org/mkin_static/vignettes/compiled_models.html).
  The autogeneration of C code was
  inspired by the [`ccSolve`](https://github.com/karlines/ccSolve) package. Thanks
  to Karline Soetaert for her work on that.
* By default, kinetic rate constants and kinetic formation fractions are
  transformed internally using
  [`transform_odeparms`](http://kinfit.r-forge.r-project.org/mkin_static/transform_odeparms.html)
  so their estimators can more reasonably be expected to follow
  a normal distribution. This has the side effect that no constraints
  are needed in the optimisation. Thanks to René Lehmann for the nice
  cooperation on this, especially the isometric logration transformation
  that is now used for the formation fractions.
* A side effect of this is that when parameter estimates are backtransformed
  to match the model definition, confidence intervals calculated from
  standard errors are also backtransformed to the correct scale, and will
  not include meaningless values like negative rate constants or 
  formation fractions adding up to more than 1, which can not occur in 
  a single experiment with a single defined radiolabel position.
* The usual one-sided t-test for significant difference from zero is nevertheless
  shown based on estimators for the untransformed parameters.
* Summary and plotting functions. The `summary` of an `mkinfit` object is in
  fact a full report that should give enough information to be able to
  approximately reproduce the fit with other tools.
* The chi-squared error level as defined in the FOCUS kinetics guidance
  (see below) is calculated for each observed variable.
* Iteratively reweighted least squares fitting is implemented in a similar way
  as in KinGUII and CAKE (see below). Simply add the argument
  `reweight = "obs"` to your call to `mkinfit` and a separate variance 
  componenent for each of the observed variables will be optimised
  in a second stage after the primary optimisation algorithm has converged.
* When a metabolite decline phase is not described well by SFO kinetics, 
  SFORB kinetics can be used for the metabolite.

## GUI

There is a graphical user interface that I consider useful for real work. Please
refer to its [documentation page](http://kinfit.r-forge.r-project.org/gmkin_static)
for installation instructions and a manual.
  
## News

Yes, there is a ChangeLog, for the latest [CRAN release](http://cran.r-project.org/web/packages/mkin/news.html)
and one for the [github master branch](https://github.com/jranke/mkin/blob/master/NEWS.md).

## Credits and historical remarks

`mkin` would not be possible without the underlying software stack consisting
of R and the packages [deSolve](http://cran.r-project.org/package=deSolve)
and [FME](http://cran.r-project.org/package=FME), to say the least.

It could not have been written without me being introduced to regulatory fate
modelling of pesticides by Adrian Gurney during my time at Harlan Laboratories
Ltd (formerly RCC Ltd). `mkin` greatly profits from and largely follows
the work done by the 
[FOCUS Degradation Kinetics Workgroup](http://focus.jrc.ec.europa.eu/dk),
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
[first CRAN version](http://cran.r-project.org/src/contrib/Archive/mkin)
on 18 May 2010.

In 2011, Bayer Crop Science started to distribute an R based successor to KinGUI named 
KinGUII whose R code is based on `mkin`, but which added, amongst other
refinements, a closed source graphical user interface (GUI), iteratively
reweighted least squares (IRLS) optimisation of the variance for each of the
observed variables, and Markov Chain Monte Carlo (MCMC) simulation
functionality, similar to what is available e.g. in the `FME` package.

Somewhat in parallel, Syngenta has sponsored the development of an `mkin` and
KinGUII based GUI application called CAKE, which also adds IRLS and MCMC, is
more limited in the model formulation, but puts more weight on usability.
CAKE is available for download from the [CAKE
website](http://showcase.tessella.com/products/cake), where you can also
find a zip archive of the R scripts derived from `mkin`, published under the GPL
license.

Finally, there is 
[KineticEval](http://github.com/zhenglei-gao/KineticEval), which contains 
a further development of the scripts used for KinGUII, so the different tools
will hopefully be able to learn from each other in the future as well.


## Development

Contributions are welcome! Your 
[mkin fork](https://help.github.com/articles/fork-a-repo) is just a mouse click
away... The master branch on github should always be in good shape, I implement 
new features in separate branches now. If you prefer subversion, project
members for the 
[r-forge project](http://r-forge.r-project.org/R/?group_id=615) are welcome as well.
Generally, the source code of the latest CRAN version should be available there.
You can also browse the source code at [cgit.jrwb.de/mkin](http://cgit.jrwb.de/mkin).
