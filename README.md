# mkin

The R package **mkin** provides calculation routines for the analysis of
chemical degradation data, including <b>m</b>ulticompartment <b>kin</b>etics as
needed for modelling the formation and decline of transformation products, or
if several compartments are involved.

## Installation

You can install the latest released version from 
[CRAN](http://cran.r-project.org/package=mkin) from within R:

```s
install.packages("mkin")
```

If looking for the latest features, you can install directly from 
[github](http://github.com/jranke/mkin), e.g.  using the `devtools` package.
Using `quick = TRUE` skips docs, multiple-architecture builds, demos, and
vignettes, to make installation as fast and painless as possible.

```s
require(devtools)
install_github("jranke/mkin", quick = TRUE)
```

## Background

In the regulatory evaluation of chemical substances like plant protection
products (pesticides), biocides and other chemicals, degradation data play an
important role. For the evaluation of pesticide degradation experiments, 
detailed guidance and helpful tools have been developed as detailed in
'Credits and historical remarks' below.

## Usage

The simplest usage example that I can think of, using model shorthand notation
(available since mkin 0.9-32) and a built-in dataset is

    library(mkin)
    fit <- mkinfit("SFO", FOCUS_2006_C)
    plot(fit, show_residuals = TRUE) 
    summary(fit)

A still very simple usage example including the definition of the same data in R
code would be

    example_data = data.frame(
      name = rep("parent", 9),
      time = c(0, 1, 3, 7, 14, 28, 63, 91, 119),
      value = c(85.1, 57.9, 29.9, 14.6, 9.7, 6.6, 4, 3.9, 0.6)
    )
    fit2 <- mkinfit("FOMC", example_data)
    plot(fit2, show_residuals = TRUE) 

A fairly complex usage example using another built-in dataset:

    data <- mkin_wide_to_long(schaefer07_complex_case, time = "time")

    model <- mkinmod(
      parent = list(type = "SFO", to = c("A1", "B1", "C1"), sink = FALSE),
      A1 = list(type = "SFO", to = "A2"),
      B1 = list(type = "SFO"),
      C1 = list(type = "SFO"),
      A2 = list(type = "SFO"), use_of_ff = "max")

    fit3 <- mkinfit(model, data, method.modFit = "Port")

    plot(fit3, show_residuals = TRUE) 
    summary(fit3)
    mkinparplot(fit3)

For more examples and to see results, have a look at the examples provided in the
[`mkinfit`](http://kinfit.r-forge.r-project.org/mkin_static/mkinfit.html)
documentation or the package vignettes referenced from the 
[mkin package documentation page](http://kinfit.r-forge.r-project.org/mkin_static/index.html)

## Features

* Highly flexible model specification using
  [`mkinmod`](http://kinfit.r-forge.r-project.org/mkin_static/mkinmod.html),
  including equilibrium reactions and using the single first-order 
  reversible binding (SFORB) model, which will automatically create
  two latent state variables for the observed variable.
* Model solution (forward modelling) in the function
  [`mkinpredict`](http://kinfit.r-forge.r-project.org/mkin_static/mkinpredict.html) 
  is performed either using the analytical solution for the case of 
  parent only degradation, an eigenvalue based solution if only simple
  first-order (SFO) or SFORB kinetics are used in the model, or
  using a numeric solver from the `deSolve` package (default is `lsoda`).
  These have decreasing efficiency, and are automatically chosen 
  by default. As of mkin 0.9-36, model solution for models with more than 
  one state observed variable (not only parent) is based on the
  [`ccSolve`](https://github.com/karlines/ccSolve) package, if installed.
  This is even faster than eigenvalue based solution, at least in the examples
  shown in the [vignette `compiled_models`](http://rawgit.com/jranke/mkin/master/vignettes/compiled_models.html)
* Model optimisation with 
  [`mkinfit`](http://kinfit.r-forge.r-project.org/mkin_static/mkinfit.html)
  internally using the `modFit` function from the `FME` package,
  but using the Port routine `nlminb` per default.
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
* Summary and plotting functions. The `summary` of an `mkinfit` object is in
  fact a full report that should give enough information to be able to
  approximately reproduce the fit with other tools.
* The chi-squared error level as defined in the FOCUS kinetics guidance
  (see below) is calculated for each observed variable.
* I recently added iteratively reweighted least squares in a similar way
  it is done in KinGUII and CAKE (see below). Simply add the argument
  `reweight = "obs"` to your call to `mkinfit` and a separate variance 
  componenent for each of the observed variables will be optimised
  in a second stage after the primary optimisation algorithm has converged.
* When a metabolite decline phase is not described well by SFO kinetics, 
  either IORE kinetics or SFORB kinetics can be used for the metabolite, 
  adding one respectively two parameters to the system.

## GUI

There is a graphical user interface that I consider useful for real work. Please
refer to its [documentation page](http://kinfit.r-forge.r-project.org/gmkin_static)
for installation instructions and a manual.
  
## News

Yes, there is a ![Changelog](NEWS.md).

## Credits and historical remarks

`mkin` would not be possible without the underlying software stack consisting
of R and the packages [deSolve](http://cran.r-project.org/package=deSolve),
[minpack.lm](http://cran.r-project.org/package=minpack.lm) and
[FME](http://cran.r-project.org/package=FME), to say the least.

It could not have been written without me being introduced to regulatory fate
modelling of pesticides by Adrian Gurney during my time at Harlan Laboratories
Ltd (formerly RCC Ltd). `mkin` greatly profits from and largely follows
the work done by the 
[FOCUS Degradation Kinetics Workgroup](http://focus.jrc.ec.europa.eu/dk),
as detailed in their guidance document from 2006, slightly updated in 2011.

Also, it was inspired by the first version of KinGUI developed by
BayerCropScience, which is based on the MatLab runtime environment.

The companion package 
[kinfit](http://kinfit.r-forge.r-project.org/kinfit_static/index.html) was 
[started in 2008](https://r-forge.r-project.org/scm/viewvc.php?view=rev&root=kinfit&revision=2) and 
[first published](http://cran.r-project.org/src/contrib/Archive/kinfit/) on
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
website](http://projects.tessella.com/cake), where you can also
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
