# mkin 0.9.49.11 (unreleased)

- Increase a test tolerance to make it pass on all CRAN check machines

# mkin 0.9.49.10 (2020-04-18)

- 'nlme.mmkin': An nlme method for mmkin row objects and an associated S3 class with print, plot, anova and endpoint methods

- 'mean_degparms, nlme_data, nlme_function': Three new functions to facilitate building nlme models from mmkin row objects

- 'endpoints': Don't return the SFORB list component if it's empty. This reduces distraction and complies with the documentation

- Article in compiled models: Add some platform specific code and suppress warnings about zero values being removed from the FOCUS D dataset

- 'plot.mmkin': Add the argument 'standardized' to avoid warnings that occurred when it was passed as part of the additional arguments captured by the dots (...)

- 'summary.mkinfit': Add AIC, BIC and log likelihood to the summary

# mkin 0.9.49.9 (2020-03-31)

- 'mkinmod': Use pkgbuild::has_compiler instead of Sys.which('gcc'), as the latter will often fail even if Rtools are installed

- 'mkinds': Use roxygen for documenting fields and methods of this R6 class

# mkin 0.9.49.8 (2020-01-09)

- 'aw': Generic function for calculating Akaike weights, methods for mkinfit objects and mmkin columns

- 'loftest': Add a lack-of-fit test

- 'plot_res', 'plot_sep' and 'mkinerrplot': Add the possibility to show standardized residuals and make it the default for fits with error models other than 'const'

- 'lrtest.mkinfit': Improve naming of the compared fits in the case of fixed parameters

- 'confint.mkinfit': Make the quadratic approximation the default, as the likelihood profiling takes a lot of time, especially if the fit has more than three parameters

# mkin 0.9.49.7 (2019-11-01)

- Fix a bug introduced in 0.9.49.6 that occurred if the direct optimisation yielded a higher likelihood than the three-step optimisation in the d_3 algorithm, which caused the fitted parameters of the three-step optimisation to be returned instead of the parameters of the direct optimisation

- Add a 'nobs' method for mkinfit objects, enabling the default 'BIC' method from the stats package. Also, add a 'BIC' method for mmkin column objects.

# mkin 0.9.49.6 (2019-10-31)

- Implement a likelihood ratio test as a method for 'lrtest' from the lmtest package

- Add an 'update' method for mkinfit objects which remembers fitted parameters if appropriate

- Add a 'residuals' method for mkinfit objects that supports scaling based on the error model

- Fix a bug in 'mkinfit' that prevented summaries of objects fitted with fixed parameters to be generated

- Add 'parms' and 'confint' methods for mkinfit objects. Confidence intervals based on the quadratic approximation as in the summary, and based on the profile likelihood

- Move long-running tests to tests/testthat/slow with a separate test log. They currently take around 7 minutes on my system

- 'mkinfit': Clean the code and return functions to calculate the log-likelihood and the sum of squared residuals

- Vignette 'twa.html': Add the maximum time weighted average formulas for the hockey stick model

- Support frameless plots ('frame = FALSE')

- Support to suppress the chi2 error level ('show_errmin = FALSE') in 'plot.mmkin'

- Update README and the introductory vignette

- Report 'OLS' as error_model_algorithm in the summary in the case that the default error_model ('const') is used

- Support summarizing 'mkinfit' objects generated with versions < 0.9.49.5

# mkin 0.9.49.5 (2019-07-04)

- Several algorithms for minimization of the negative log-likelihood for non-constant error models (two-component and variance by variable). In the case the error model is constant variance, least squares is used as this is more stable. The default algorithm 'd_3' tries direct minimization and a three-step procedure, and returns the model with the highest likelihood.

- The argument 'reweight.method' to mkinfit and mmkin is now obsolete, use 'error_model' and 'error_model_algorithm' instead

- Add a test that checks if we get the best known AIC for parent only fits to 12 test datasets. Add these test datasets for this purpose.

- New function 'mkinerrplot'. This function is also used for residual plots in 'plot.mmkin' if the argument 'resplot = "errmod"' is given, and in 'plot.mkinfit' if 'show_errplot' is set to TRUE.

- Remove dependency on FME, only use nlminb for optimisation ('Port' algorithm). I cannot remember cases where one of the other optimisation algorithms was preferable, except that I sometime used Levenberg-Marquardt for speed in cases where I did not expect to get trapped in a local minimum.

- Use the numDeriv package to calculate hessians. This results in slightly different confidence intervals, takes a bit longer, but is apparently more robust

- Add a simple benchmark vignette to document the impact on performance.

- The code for manual weighting was removed. This functionality might get added again at a later time. For the time being, please use an earlier version, e.g. 0.9.48.1 if you want to do manual weighting.

- The fitting time reported in the summary now includes the time used for calculation of the hessians

- Adapt tests

- Fix an error in the FOCUS chi2 error level calculations that occurred if parameters were specified in parms.ini that were not in the model. A warning was already issued, but when fitting in parallel via mmkin this could go unnoticed.

- Add example datasets obtained from risk assessment reports published by the European Food Safety Agency.

# mkin 0.9.48.1 (2019-03-04)

- Add the function 'logLik.mkinfit' which makes it possible to calculate an AIC for mkinfit objects

- Add the function 'AIC.mmkin' to make it easy to compare columns of mmkin objects

- 'add_err': Respect the argument giving the number of replicates in the synthetic dataset

- 'max_twa_parent': Support maximum time weighted average concentration calculations for the hockey stick (HS) model

- 'mkinpredict': Make the function generic and create a method for mkinfit objects

- 'mkinfit': Improve the correctness of the fitted two component error model by fitting the mean absolute deviance at each observation against the observed values, weighting with the current two-component error model

- 'tests/testthat/test_irls.R': Test if the components of the error model used to generate the data can be reproduced with moderate accuracy

- Add the function 'CAKE_export' to facilitate cross-checking of results

- Implement the logistic model (only tested for parent fits)

- 'nafta': Add evaluations according to the NAFTA guidance

# mkin 0.9.47.5 (2018-09-14)

- Make the two-component error model stop in cases where it is inadequate to avoid nls crashes on windows

- Move two vignettes to a location where they will not be built on CRAN (to avoid more NOTES from long execution times)

- Exclude more example code from testing on CRAN to avoid NOTES from long execution times

# mkin 0.9.47.3

- 'mkinfit': Improve fitting the error model for reweight.method = 'tc'. Add 'manual' to possible arguments for 'weight'

- Test that FOCUS_2006_C can be evaluated with DFOP and reweight.method = 'tc'

# mkin 0.9.47.2 (2018-07-19)

- 'sigma_twocomp': Rename 'sigma_rl' to 'sigma_twocomp' as the Rocke and Lorenzato model assumes lognormal distribution for large y. Correct references to the Rocke and Lorenzato model accordingly.

- 'mkinfit': Use 1.1 as starting value for N parameter of IORE models to obtain convergence in more difficult cases. Show parameter names when 'trace_parms' is 'TRUE'.

# mkin 0.9.47.1 (2018-02-06)

- Skip some tests on CRAN and winbuilder to avoid timeouts

- 'test_data_from_UBA_2014': Added this list of datasets containing experimental data used in the expertise from 2014

- 'mkinfit': Added the iterative reweighting method 'tc' using the two-component error model from Rocke and Lorenzato. NA values in the data are not returned any more.

- 'mkinfit': Work around a bug in the current FME version that prevented the convergence message to be returned in the case of non-convergence.

- 'summary.mkinfit': Improved output regarding weighting method. No predictions are returned for NA values in the model (see above).

- 'summary.mkinfit': Show versions of mkin and R used for fitting (not the ones used for the summary) if the fit was generated with mkin >= 0.9.47.1

# mkin 0.9.46.3 (2017-11-16)

- `README.md`, `vignettes/mkin.Rmd`: URLs were updated

- `synthetic_data_for_UBA`: Add the code used to generate the data in the interest of reproducibility

# mkin 0.9.46.2 (2017-10-10)

- Converted the vignette FOCUS_Z from tex/pdf to markdown/html

- `DESCRIPTION`: Add ORCID

# mkin 0.9.46.1 (2017-09-14)

- `plot.mkinfit`: Fix scaling of residual plots for the case of separate plots for each observed variable

- `plot.mkinfit`: Use all data points of the fitted curve for y axis scaling  for the case of separate plots for each observed variable

- Documentation updates

# mkin 0.9.46 (2017-07-24)

- Remove `test_FOMC_ill-defined.R` as it is too platform dependent

# mkin 0.9.45.2 (2017-07-24)

- Rename `twa` to `max_twa_parent` to avoid conflict with `twa` from my `pfm` package

- Update URLs in documentation

- Limit test code to one core to pass on windows

- Switch from `microbenchmark` to `rbenchmark` as the former is not supported on all platforms

# mkin 0.9.45.1 (2016-12-20)

## New features

- A `twa` function, calculating maximum time weighted average concentrations for the parent (SFO, FOMC and DFOP).

# mkin 0.9.45 (2016-12-08)

## Minor changes

- `plot.mkinfit` and `plot.mmkin`: If the plotting device is `tikz`, LaTeX markup is being used for the chi2 error in the graphs.

- Use `pkgdown`, the successor of `staticdocs` for generating static HTML documentation. Include example output and graphs also for `dontrun` sections.

- `plot.mkinfit`: Plotting does not fail any more when the compiled model is not available, e.g. because it was removed from the temporary directory. In this case, the uncompiled model is now used for plotting

# mkin 0.9.44 (2016-06-29)

## Bug fixes

- The test `test_FOMC_ill-defined` failed on several architectures, so the test is now skipped

# mkin 0.9.43 (2016-06-28)

## Major changes

- The title was changed to `Kinetic evaluations of chemical degradation data`

- `plot.mkinfit`: Add the possibility to show fits (and residual plots if requested) separately for the observed variables

- `plot.mkinfit`: Add the possibility to show the chi2 error levels in the plot, similar to the way they are shown in `plot.mmkin`

- `plot_sep`: Add this function as a convenience wrapper for plotting observed variables of mkinfit objects separately, with chi2 error values and residual plots.

- Vignettes: The main vignette `mkin` was converted to R markdown and updated. The other vignettes were also updated to show current improved functionality.

- The function `add_err` was added to the package, making it easy to generate simulated data using an error model based on the normal distribution

## Minor changes

- Remove an outdated reference to the inline package in the `compiled_models` vignette

- `mkinfit`: Do not error out in cases where the fit converges, but the Jacobian for the untransformed model cost can not be estimated. Give a warning instead and return NA for the t-test results.

- `summary.mkinfit`: Give a warning message when the covariance matrix can not be obtained.

- A test has been added to containing a corresponding edge case to check that the warnings are correctly issued and the fit does not terminate.

- `plot.mmkin`: Round the chi2 error value to three significant digits, instead of two decimal digits.

- `mkinfit`: Return the `err` values used on weighted fits as a column named `err`. Also include these inverse weights when the column `value` in the observed data is used, which is returned as `observed` in the data component of the mkinfit object.

## Bug fixes

- `endpoints`: When the name of a substance degrading to a metabolite (e.g. a parent compound) used in the model formulation ended in the letter `f`, some rate parameters could be listed as formation fractions with mixed up names. These would also appear in the summary.

- `mkinfit`: Check for all observed variables when checking if the user tried to fix formation fractions when fitting them using ilr transformation.

- `plot.mmkin`: Set the plot margins correctly, also in the case of a single fit to be plotted, so the main title is placed in a reasonable way.

- `plot.mkinfit`: Correct default values for `col_obs`, `pch_obs` and `lty_obs` for the case that `obs_vars` is specified.

# mkin 0.9.42 (2016-03-25)

## Major changes

- Add the argument `from_max_mean` to `mkinfit`, for fitting only the decline from the maximum observed value for models with a single observed variable

## Minor changes

- Add plots to `compiled_models` vignette

- Give an explanatory error message when mkinmod fails due to a missing definition of a target variable

- `print.mkinmod()`: Improve formatting when printing mkinmod model definitions

# mkin 0.9-41 (2015-11-09)

## Minor changes

- Add an R6 class `mkinds` representing datasets with a printing method

- Add a printing method for mkinmod objects

- Make it possible to specify arbitrary strings as names for the compounds in `mkinmod`, and show them in the plot

- Use an index.r file to group help topics in static documentation

## Bug fixes

- `print.summary.mkinfit()`: Avoid an error that occurred when printing summaries generated with mkin versions before 0.9-36

# mkin 0.9-40 (2015-07-21)

## Bug fixes

- `endpoints()`: For DFOP and SFORB models, where `optimize()` is used, make use of the fact that the DT50 must be between DT50_k1 and DT50_k2 (DFOP) or DT50_b1 and DT50_b2 (SFORB), as `optimize()` sometimes did not find the minimum. Likewise for finding DT90 values. Also fit on the log scale to make the function more efficient.

## Internal changes

- `DESCRIPTION`, `NAMESPACE`, `R/*.R`: Import (from) stats, graphics and methods packages, and qualify some function calls for non-base packages installed with R to avoid NOTES made by R CMD check --as-cran with upcoming R versions.

# mkin 0.9-39 (2015-06-26)

## Major changes

- New function `mmkin()`: This function takes a character vector of model shorthand names, or alternatively a list of mkinmod models, as well as a list of dataset as main arguments. It returns a matrix of mkinfit objects, with a row for each model and a column for each dataset. A subsetting method with single brackets is available. Fitting the models in parallel using the `parallel` package is supported.

- New function `plot.mmkin()`: Plots single-row or single-column `mmkin` objects including residual plots.

## Bug fixes

- `mkinparplot()`: Fix the x axis scaling for rate constants and formation fractions that got confused by the introduction of the t-values of transformed parameters.

# mkin 0.9-38 (2015-06-24)

## Minor changes

- `vignettes/compiled_models.html`: Show the performance improvement factor actually obtained when building the vignette, as well as mkin version, some system info and the CPU model used for building the vignette.

- `GNUMakefile`,`vignettes/*`: Clean up vignette generation and include table of contents in HTML vignettes.

## Bug fixes

- `mkinmod()`: When generating the C code for the derivatives, only declare the time variable when it is needed and remove the '-W-no-unused-variable' compiler flag as the C compiler used in the CRAN checks on Solaris does not know it.

# mkin 0.9-36 (2015-06-21)

## Major changes

- `summary.mkinfit()`: A one-sided t-test for significant difference of untransformed parameters from zero is now always shown, based on the assumption of normal distribution for estimators of all untransformed parameters. Use with caution, as this assumption is unrealistic e.g. for rate constants in these nonlinear kinetic models.

- If a compiler (gcc) is installed, use a version of the differential equation model compiled from C code, which is a huge performance boost for models where only the deSolve method works.

- `mkinmod()`: Create a list component $cf (of class CFuncList) in the list returned by mkinmod, if a version can be compiled from autogenerated C code (see above).

- `mkinfit()`: Set the default `solution_type` to `deSolve` when a compiled version of the model is present, except when an analytical solution is possible.

## Minor changes

- Added a simple showcase vignette with an evaluation of FOCUS example dataset D

# mkin 0.9-35 (2015-05-15)

## Major changes

- Switch from RUnit to testthat for testing

## Bug fixes

- `mkinparplot()`: Avoid warnings that occurred when not all confidence intervals were available in the summary of the fit

- `print.summary.mkinfit()`: Fix printing the summary for the case that the number of iterations is not available

- NAMESPACE: export S3 methods plot.mkinfit, summary.mkinfit and print.summary.mkinfit to satisfy R CMD check on R-devel

- `mkinparplot()`: Avoid warning in R CMD check about undeclared global variable `Lower`

## New features

- `mkinfit()`: Report successful termination when quiet = FALSE. This is helpful for more difficult problems fitted with reweight.method = obs, as no progress is often indicated during the reweighting.

- A first test using results established in the expertise written for the German Federal Environmental Agency (UBA) was added.

- Add synthetic datasets generated for expertise written for the German Federal Environmental Agency UBA

- Add tests based on these datasets

# mkin 0.9-34 (2014-11-22)

## New features

- Add the convenience function `mkinsub()` for creating the lists used in `mkinmod()`

- Add the possibility to fit indeterminate order rate equation (IORE) models using an analytical solution (parent only) or a numeric solution. Paths from IORE compounds to metabolites are supported when using formation fractions (use_of_ff = 'max'). Note that the numerical solution (method.ode = 'deSolve') of the IORE differential equations sometimes fails due to numerical problems.

- Switch to using the Port algorithm (using a model/trust region approach) per default. While needing more iterations than the Levenberg-Marquardt algorithm previously used per default, it is less sensitive to starting parameters.

## Minor changes

- The formatting of differential equations in the summary was further improved

- Always include 0 on y axis when plotting during the fit

# mkin 0.9-33 (2014-10-22)

## New features

- The initial value (state.ini) for the observed variable with the highest observed residue is set to 100 in case it has no time zero observation and `state.ini = "auto"`

- A basic unit test for `mkinerrmin()` was written

## Bug fixes

- `mkinfit()`: The internally fitted parameter for `g` was named `g_ilr` even when `transform_fractions=FALSE`

- `mkinfit()`: The initial value (state.ini) for the parent compound was not set when the parent was not the (only) variable with the highest value in the observed data.

- `mkinerrmin()`: When checking for degrees of freedom for metabolites, check if their time zero value is fixed instead of checking if the observed value is zero. This ensures correct calculation of degrees of freedom also in cases where the metabolite residue at time zero is greater zero.

- `plot.mkinfit()`: Avoid a warning message about only using the first component of ylim that occurred when ylim was specified explicitly

## Minor changes

- The formatting of differential equations in the summary was improved by wrapping overly long lines

- The FOCUS_Z vignette was rebuilt with the above improvement and using a width of 70 to avoid output outside of the grey area

- `print.summary.mkinfit()`: Avoid a warning that occurred when gmkin showed summaries ofinitial fits without iterations

- `mkinfit()`: Avoid a warning that occurred when summarising a fit that was performed with maxitmodFit = 0 as done in gmkin for configuring new fits.

# mkin 0.9-32 (2014-07-24)

## New features

- The number of degrees of freedom is difficult to define in the case of ilr transformation of formation fractions. Now for each source compartment the number of ilr parameters (=number of optimised parameters) is divided by the number of pathways to metabolites (=number of affected data series) which leads to fractional degrees of freedom in some cases.

- The default for the initial value for the first state value is now taken from the mean of the observations at time zero, if available.

- The kinetic model can alternatively be specified with a shorthand name for parent only degradation models, e.g. `SFO`, or `DFOP`.

- Optimisation method, number of model evaluations and time elapsed during optimisation are given in the summary of mkinfit objects.

- The maximum number of iterations in the optimisation algorithm can be specified using the argument `maxit.modFit` to the mkinfit function.

- mkinfit gives a warning when the fit does not converge (does not apply to SANN method). This warning is included in the summary.

## Bug fixes

- Avoid plotting an artifical 0 residual at time zero in `mkinresplot`

- In the determination of the degrees of freedom in `mkinerrmin`, formation fractions were accounted for multiple times in the case of parallel formation of metabolites. See the new feature described above for the solution.

- `transform_rates=FALSE` in `mkinfit` now also works for FOMC and HS models.

- Initial values for formation fractions were not set in all cases.

- No warning was given when the fit did not converge when a method other than the default Levenberg-Marquardt method `Marq` was used.

## Minor changes

- Vignettes were rebuilt to reflect the changes in the summary method.

- Algorithm `Pseudo` was excluded because it needs user-defined parameter limits which are not supported.

- Algorithm `Newton` was excluded because of its different way to specify the maximum number of iterations and because it does not appear to provide additional benefits.

# mkin 0.9-31 (2014-07-14)

## Bug fixes

- The internal renaming of optimised parameters in Version 0.9-30 led to errors in the determination of the degrees of freedom for the chi2 error level calulations in `mkinerrmin()` used by the summary function.

# mkin 0.9-30 (2014-07-11)

## New features

- It is now possible to use formation fractions in combination with turning off the sink in `mkinmod()`.

## Major changes

- The original and the transformed parameters now have different names (e.g. `k_parent` and `log_k_parent`. They also differ in how many they are when we have formation fractions but no pathway to sink.

- The order of some of the information blocks in `print.summary.mkinfit.R()` has been ordered in a more logical way.

## Minor changes

- The vignette FOCUS_Z has been simplified to use formation fractions with turning off the sink, and slightly amended to use the new versions of DT50 values calculated since mkin 0.9-29.

- All vignettes have been rebuilt so they reflect all changes.

- The ChangeLog was renamed to NEWS.md and the entries were converted to markdown syntax compatible with the `tools::news()` function built into R.

- The test suite was overhauled. Tests of the DFOP and SFORB models with dataset FOCUS_2006_A were removed, as they were too much dependent on the optimisation algorithm and/or starting parameters, because the dataset is SFO (compare kinfit vignette).

- Also, the Schaefer complex case can now be fitted using formation fractions, and with the 'Port' optimisation method we also fit A2 in the same way as published in the Piacenza paper.

- Some more checks were introduced to `mkinfit()`, leading to warnings or stopping execution if unsupported combinations of methods and parameters are requested.

# mkin 0.9-29 (2014-06-27)

- R/mkinresplot.R: Make it possible to specify `xlim`

- R/geometric_mean.R, man/geometric_mean.Rd: Add geometric mean function

- R/endpoints.R, man/endpoints.Rd: Calculate additional (pseudo)-DT50 values for FOMC, DFOP, HS and SFORB. Avoid calculation of formation fractions from rate constants when they are directly fitted

# mkin 0.9-28 (2014-05-20)

- Do not backtransform confidence intervals for formation fractions if more than one compound is formed, as such parameters only define the pathways as a set

- Add historical remarks and some background to the main package vignette

- Correct 'isotropic' into 'isometric' for the ilr transformation

# mkin 0.9-27 (2014-05-10)

- Fork the GUI into a separate package [gmkin](http://github.com/jranke/gmkin)

- DESCRIPTION, NAMESPACE, TODO: Adapt and add copyright information

- Remove files belonging to the GUI

- Possibility to fit without parameter transformations, using bounds as implemented in FME

- Add McCall 2,4,5-T dataset

- Enable selection of observed variables in plotting

- Add possibility to show residual plot in `plot.mkinfit`

- R/mkinparplot.R, man/mkinparplot.Rd: plot parameters with confidence intervals

- Change vignette format from Sweave to knitr

- Split examples vignette to FOCUS_L and FOCUS_Z

- Remove warning about constant formation fractions in mkinmod as it was based on a misconception

- Restrict the unit test with the Schaefer data to parent and primary metabolites as formation fraction and DT50 for A2 are higly correlated and passing the test is platform dependent. For example, the test fails in 1 out of 14 platforms on CRAN as of today.

- Add Eurofins Regulatory AG copyright notices

- Import FME and deSolve instead of depending on them to have clean startup

- Add a starter function for the GUI: `gmkin()`

- Change the format of the workspace files of gmkin so they can be distributed and documented in the package

- Add gmkin workspace datasets FOCUS_2006_gmkin and FOCUS_2006_Z_gmkin

# mkin 0.9-24 (2013-11-06)

- Bugfix re-enabling the fixing of any combination of initial values for state variables

- Default values for kinetic rate constants are not all 0.1 any more but are "salted" with a small increment to avoid numeric artefacts with the eigenvalue based solutions

- Backtransform fixed ODE parameters for the summary

# mkin 0.9-22 (2013-10-26)

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
