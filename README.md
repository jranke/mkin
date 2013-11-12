# mkin

The R package **mkin** provides calculation routines for the analysis of chemical
degradation data, including **m**ulticompartment **kin**etics as needed for modelling
the formation and decline of transformation products.

## Installation

You can install the latest released version from 
[CRAN](http://cran.r-project.org/package=mkin) from within R:

```s
install.packages('mkin')
```

A development version is usually available from [R-Forge](http://r-forge.r-project.org/R/?group_id=615):

```s
install.packages('mkin', repos = 'http://r-forge.r-project.org')
```

If R-Forge is lacking behind or if you prefer, you can install directly from
github using the `devtools` package:

```s
require(devtools)
install_github("mkin", "jranke")
```

## Usage

For a start, have a look at the examples provided in the 
[`mkinfit`](http://kinfit.r-forge.r-project.org/mkin_static/mkinfit.html)
Documentation 
or the package vignettes referenced from the 
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
  by default.
* Model optimisation with 
  [`mkinfit`](http://kinfit.r-forge.r-project.org/mkin_static/mkinfit.html)
  internally using the `modFit` function from the `FME` package,
  which uses the least-squares Levenberg-Marquardt algorithm from
  `minpack.lm` per default.
* Kinetic rate constants and kinetic formation fractions are transformed 
  internally using
  [`transform_odeparms`](http://kinfit.r-forge.r-project.org/mkin_static/transform_odeparms.html)
  so their estimators can more reasonably be expected to follow
  a normal distribution. This has the side effect that no constraints
  are needed in the optimisation. Thanks to René Lehmann for the nice
  cooperation on this, especially the isotropic logration transformation
  that is now used for the formation fractions.
* A side effect of this is that when parameter estimates are backtransformed
  to match the model definition, confidence intervals calculated from
  standard errors are also backtransformed to the correct scale, and will
  not include meaningless values (like negative rate constants or 
  formation fractions adding up to more than 1, which can not occur in 
  a single experiment with a single defined radiolabel position).
* Summary and plotting functions. The `summary` of an `mkinfit` object is in
  fact a full report that should give enough information to be able to
  approximately reproduce the fit with other tools.
* I recently added iteratively reweighted least squares in a similar way
  it is done in KinGUII and CAKE (see below). Simply add the argument
  `reweight = "obs"` to your call to `mkinfit` and a separate variance 
  componenent for each of the observed variables will be optimised
  in a second stage after the primary optimisation algorithm has converged.

## Credits

`mkin` would not be possible without the underlying software stack consisting
of R and the packages [deSolve](http://cran.r-project.org/package=deSolve),
[minpack.lm](http://cran.r-project.org/package=minpack.lm) and
[FME](http://cran.r-project.org/package=FME), to say the least.

Also, it was inspired by the first version of KinGUI developed by
BayerCropScience, which is based on the MatLab runtime environment.

Bayer has developed a successor named KinGUII whose R code is based on `mkin`, but which
added, amongst other refinements, a closed source graphical user interface
(GUI), iteratively reweighted least squares (IRLS) optimisation of the variance
for each of the observed variables, and Markov Chain Monte Carlo (MCMC)
simulation functionality, similar to what is available e.g. in the `FME`
package.

Syngenta has sponsored the development of an `mkin` (and KinGUII?) based GUI
application called CAKE, which adds IRLS and MCMC, is more limited in the model
formulation, but puts more weight on usability.  CAKE is available for download
from the [CAKE website](http://projects.tessella.com/cake), where you can also
find a zip archive of the R scripts derived from `mkin`, published under the GPL
license.

Finally, I just (2013-11-11) noticed the github repositories
[StudyKin](http://github.com/zhenglei-gao/StudyKin) and
[KineticEval](http://github.com/zhenglei-gao/KineticEval), the latter of which appears to be 
<<<<<<< HEAD
actively developed.
=======
actively developed, so the different tools will hopefully be able to learn
from each other in the future as well.
>>>>>>> master
