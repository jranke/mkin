# CHANGES in mkin VERSION 0.9-32

## NEW FEATURES

- The number of degrees of freedom is difficult to define in the case of ilr transformation of formation fractions. Now for each source compartment the number of ilr parameters (=number of optimised parameters) is divided by the number of pathways to metabolites (=number of affected data series) which leads to fractional degrees of freedom in some cases.

- The default for the initial value for the first state value is now taken from the mean of the observations at time zero, if available.

- The kinetic model can alternatively be specified with a shorthand name for parent only degradation models, e.g. `SFO`, or `DFOP`.

- Optimisation method, number of model evaluations and time elapsed during optimisation are given in the summary of mkinfit objects.

- The maximum number of iterations in the optimisation algorithm can be specified using the argument `maxit.modFit` to the mkinfit function.

- mkinfit gives a warning when the fit does not converge (does not apply to SANN method). This warning is included in the summary.

## BUG FIXES

- Avoid plotting an artifical 0 residual at time zero in `mkinresplot`

- In the determination of the degrees of freedom in `mkinerrmin`, formation fractions were accounted for multiple times in the case of parallel formation of metabolites. See the new feature described above for the solution.

- `transform_rates=FALSE` in `mkinfit` now also works for FOMC and HS models.

- Initial values for formation fractions were not set in all cases.

- No warning was given when the fit did not converge when a method other than the default Levenberg-Marquardt method `Marq` was used.

## MINOR CHANGES

- Vignettes were rebuilt to reflect the changes in the summary method.

- Algorithm `Pseudo` was excluded because it needs user-defined parameter limits which are not supported.

- Algorithm `Newton` was excluded because of its different way to specify the maximum number of iterations and because it does not appear to provide additional benefits.

# CHANGES in mkin VERSION 0.9-31

## BUG FIXES

- The internal renaming of optimised parameters in Version 0.9-30 led to errors in the determination of the degrees of freedom for the chi2 error level calulations in `mkinerrmin()` used by the summary function.

# CHANGES in mkin VERSION 0.9-30 

## NEW FEATURES

- It is now possible to use formation fractions in combination with turning off the sink in `mkinmod()`.

## MAJOR CHANGES

- The original and the transformed parameters now have different names (e.g. `k_parent` and `log_k_parent`. They also differ in how many they are when we have formation fractions but no pathway to sink.

- The order of some of the information blocks in `print.summary.mkinfit.R()` has been ordered in a more logical way.

## MINOR CHANGES

- The vignette FOCUS_Z has been simplified to use formation fractions with turning off the sink, and slightly amended to use the new versions of DT50 values calculated since mkin 0.9-29.

- All vignettes have been rebuilt so they reflect all changes.

- The ChangeLog was renamed to NEWS.md and the entries were converted to markdown syntax compatible with the `tools::news()` function built into R.

- The test suite was overhauled. Tests of the DFOP and SFORB models with dataset FOCUS_2006_A were removed, as they were too much dependent on the optimisation algorithm and/or starting parameters, because the dataset is SFO (compare kinfit vignette).

- Also, the Schaefer complex case can now be fitted using formation fractions, and with the 'Port' optimisation method we also fit A2 in the same way as published in the Piacenza paper.

- Some more checks were introduced to `mkinfit()`, leading to warnings or stopping execution if unsupported combinations of methods and parameters are requested.

# CHANGES in mkin VERSION 0.9-29

- R/mkinresplot.R: Make it possible to specify `xlim`

- R/geometric_mean.R, man/geometric_mean.Rd: Add geometric mean function

- R/endpoints.R, man/endpoints.Rd: Calculate additional (pseudo)-DT50 values for FOMC, DFOP, HS and SFORB. Avoid calculation of formation fractions from rate constants when they are directly fitted

# CHANGES in mkin VERSION 0.9-28

- Do not backtransform confidence intervals for formation fractions if more than one compound is formed, as such parameters only define the pathways as a set

- Add historical remarks and some background to the main package vignette

- Correct 'isotropic' into 'isometric' for the ilr transformation

# CHANGES in mkin VERSION 0.9-27

- Fork the GUI into a separate package [gmkin](http://github.com/jranke/gmkin)

- DESCRIPTION, NAMESPACE, TODO: Adapt and add copyright information

- Remove files belonging to the GUI

- Possibility to fit without parameter transformations, using bounds as implemented in FME

- Add McCall 2,4,5-T dataset

- Enable selection of observed variables in plotting

- Add possibility to show residual plot in `plot.mkinfit`

- R/mkinparplot.R, man/mkinparplot.Rd: plot parameters with confidence intervals

# CHANGES in mkin VERSION 0.9-25

- Change vignette format from Sweave to knitr

- Split examples vignette to FOCUS_L and FOCUS_Z

- Remove warning about constant formation fractions in mkinmod as it was based on a misconception

- Restrict the unit test with the Schaefer data to parent and primary metabolites as formation fraction and DT50 for A2 are higly correlated and passing the test is platform dependent. For example, the test fails in 1 out of 14 platforms on CRAN as of today.

- Add Eurofins Regulatory AG copyright notices

- Import FME and deSolve instead of depending on them to have clean startup

- Add a starter function for the GUI: `gmkin()`

- Change the format of the workspace files of gmkin so they can be distributed and documented in the package

- Add gmkin workspace datasets FOCUS_2006_gmkin and FOCUS_2006_Z_gmkin

# CHANGES in mkin VERSION 0.9-24

- Bugfix re-enabling the fixing of any combination of initial values for state variables

- Default values for kinetic rate constants are not all 0.1 any more but are "salted" with a small increment to avoid numeric artefacts with the eigenvalue based solutions

# CHANGES in mkin VERSION 0.9-23

- Backtransform fixed ODE parameters for the summary

# CHANGES in mkin VERSION 0.9-22

- Get rid of the optimisation step in `mkinerrmin` - this was unnecessary. Thanks to KinGUII for the inspiration - actually this is equation 6-2 in FOCUS kinetics p. 91 that I had overlooked originally

- Fix `plot.mkinfit` as it passed graphical arguments like main to the solver

- Do not use `plot=TRUE` in `mkinfit()` example

- The first successful fits in the not so simple GUI

- Fix iteratively reweighted least squares for the case of many metabolites

- Unify naming of initial values of state variables

- Unify naming in dataframes of optimised and fixed parameters in the summary

- Show the weighting method for residuals in the summary

- Correct the output of the data in the case of manual weighting

- Implement IRLS assuming different variances for observed variables

- Do not use 0 values at time zero for chi2 error level calculations. This is the way it is done in KinGUII and it makes sense. It does impact the chi2 error levels in the output. Generally they seem to be lower for metabolites now, presumably because the mean of the observed values is higher

For a detailed list of changes to the mkin source please consult the commit history on http://github.com/jranke/mkin
