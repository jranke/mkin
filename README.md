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
[mkinfit Documentation](http://kinfit.r-forge.r-project.org/mkin_static/mkinfit.html)
or the package vignettes referenced from the 
[mkin package documentation page](http://kinfit.r-forge.r-project.org/mkin_static/mkinfit.html)

## Credits

`mkin` would not be possible without the underlying software stack consisting
of R and the packages [deSolve](http://cran.r-project.org/package=deSolve),
[minpack.lm](http://cran.r-project.org/package=minpack.lm) and
[FME](http://cran.r-project.org/package=FME), to say the least.

Also, it was inspired by the first version of KinGUI developed by BayerCropScience, which is based 
on the MatLab runtime environment.

Bayer has developed a successor named KinGUII which is based on mkin, but which added, amongst other 
refinements, a graphical user interface (GUI, which is not open source, iteratively
reweighted least squares (IRLS) optimisation of the variance for each of the observed variables,
and Markov Chain Monte Carlo (MCMC) simulation functionality, similar to what is available e.g. in the 
`FME` package.

Syngenta has sponsored the development of an mkin based GUI application called CAKE, which 
adds IRLS and MCMC, is more limited in the model formulation, but puts more weight on usability.
CAKE is available for download from the [CAKE website](http://projects.tessella.com/cake), where 
you can also find a zip archive of the R scripts derived from mkin, published under the GPL
license.

Finally, I just noticed the github repositories
[StudyKin](http://github.com/zhenglei-gao/StudyKin) and
[KineticEval](http://github.com/zhenglei-gao/KineticEval), the latter of which appears to be 
actively developed.
