<!--
%\VignetteEngine{knitr::knitr}
%\VignetteIndexEntry{Example evaluation of FOCUS Laboratory Data L1 to L3}
-->

# Example evaluation of FOCUS Laboratory Data L1 to L3

## Laboratory Data L1

The following code defines example dataset L1 from the FOCUS kinetics
report, p. 284


```r
library("mkin")
```

```
## Loading required package: FME
## Loading required package: deSolve
## Loading required package: rootSolve
## Loading required package: minpack.lm
## Loading required package: MASS
## Loading required package: coda
## Loading required package: lattice
```

```r
FOCUS_2006_L1 = data.frame(t = rep(c(0, 1, 2, 3, 5, 7, 14, 21, 30), each = 2), 
    parent = c(88.3, 91.4, 85.6, 84.5, 78.9, 77.6, 72, 71.9, 50.3, 59.4, 47, 
        45.1, 27.7, 27.3, 10, 10.4, 2.9, 4))
FOCUS_2006_L1_mkin <- mkin_wide_to_long(FOCUS_2006_L1)
```


The next step is to set up the models used for the kinetic analysis. Note that
the model definitions contain the names of the observed variables in the data.
In this case, there is only one variable called `parent`.


```r
SFO <- mkinmod(parent = list(type = "SFO"))
FOMC <- mkinmod(parent = list(type = "FOMC"))
DFOP <- mkinmod(parent = list(type = "DFOP"))
```


The three models cover the first assumption of simple first order (SFO),
the case of declining rate constant over time (FOMC) and the case of two
different phases of the kinetics (DFOP). For a more detailed discussion
of the models, please see the FOCUS kinetics report.

The following two lines fit the model and produce the summary report
of the model fit. This covers the numerical analysis given in the 
FOCUS report.


```r
m.L1.SFO <- mkinfit(SFO, FOCUS_2006_L1_mkin, quiet = TRUE)
summary(m.L1.SFO)
```

```
## mkin version:    0.9.25 
## R version:       3.0.2 
## Date of fit:     Sun Nov 17 15:02:54 2013 
## Date of summary: Sun Nov 17 15:02:54 2013 
## 
## Equations:
## [1] d_parent = - k_parent_sink * parent
## 
## Method used for solution of differential equation system:
## analytical 
## 
## Weighting: none
## 
## Starting values for optimised parameters:
##               value   type transformed
## parent_0      100.0  state     100.000
## k_parent_sink   0.1 deparm      -2.303
## 
## Fixed parameter values:
## None
## 
## Optimised, transformed parameters:
##               Estimate Std. Error Lower Upper
## parent_0         92.50     1.3700 89.60 95.40
## k_parent_sink    -2.35     0.0406 -2.43 -2.26
## 
## Backtransformed parameters:
##               Estimate   Lower  Upper
## parent_0       92.5000 89.6000 95.400
## k_parent_sink   0.0956  0.0877  0.104
## 
## Residual standard error: 2.95 on 16 degrees of freedom
## 
## Chi2 error levels in percent:
##          err.min n.optim df
## All data    3.42       2  7
## parent      3.42       2  7
## 
## Estimated disappearance times:
##        DT50 DT90
## parent 7.25 24.1
## 
## Estimated formation fractions:
##             ff
## parent_sink  1
## 
## Parameter correlation:
##               parent_0 k_parent_sink
## parent_0         1.000         0.625
## k_parent_sink    0.625         1.000
## 
## Data:
##  time variable observed predicted residual
##     0   parent     88.3     92.47   -4.171
##     0   parent     91.4     92.47   -1.071
##     1   parent     85.6     84.04    1.561
##     1   parent     84.5     84.04    0.461
##     2   parent     78.9     76.38    2.524
##     2   parent     77.6     76.38    1.224
##     3   parent     72.0     69.41    2.588
##     3   parent     71.9     69.41    2.488
##     5   parent     50.3     57.33   -7.030
##     5   parent     59.4     57.33    2.070
##     7   parent     47.0     47.35   -0.352
##     7   parent     45.1     47.35   -2.252
##    14   parent     27.7     24.25    3.453
##    14   parent     27.3     24.25    3.053
##    21   parent     10.0     12.42   -2.416
##    21   parent     10.4     12.42   -2.016
##    30   parent      2.9      5.25   -2.351
##    30   parent      4.0      5.25   -1.251
```


A plot of the fit is obtained with the plot function for mkinfit objects.


```r
plot(m.L1.SFO)
```

![plot of chunk unnamed-chunk-4](figure/unnamed-chunk-4.png) 

The residual plot can be easily obtained by


```r
mkinresplot(m.L1.SFO, ylab = "Observed", xlab = "Time")
```

![plot of chunk unnamed-chunk-5](figure/unnamed-chunk-5.png) 


For comparison, the FOMC model is fitted as well, and the chi^2 error level
is checked.


```r
m.L1.FOMC <- mkinfit(FOMC, FOCUS_2006_L1_mkin, quiet = TRUE)
summary(m.L1.FOMC, data = FALSE)
```

```
## mkin version:    0.9.25 
## R version:       3.0.2 
## Date of fit:     Sun Nov 17 15:02:55 2013 
## Date of summary: Sun Nov 17 15:02:55 2013 
## 
## Equations:
## [1] d_parent = - (alpha/beta) * ((time/beta) + 1)^-1 * parent
## 
## Method used for solution of differential equation system:
## analytical 
## 
## Weighting: none
## 
## Starting values for optimised parameters:
##          value   type transformed
## parent_0   100  state     100.000
## alpha        1 deparm       0.000
## beta        10 deparm       2.303
## 
## Fixed parameter values:
## None
## 
## Optimised, transformed parameters:
##          Estimate Std. Error Lower Upper
## parent_0     92.5         NA    NA    NA
## alpha        25.6         NA    NA    NA
## beta         28.0         NA    NA    NA
## 
## Backtransformed parameters:
##          Estimate Lower Upper
## parent_0 9.25e+01    NA    NA
## alpha    1.35e+11    NA    NA
## beta     1.41e+12    NA    NA
## 
## Residual standard error: 3.05 on 15 degrees of freedom
## 
## Chi2 error levels in percent:
##          err.min n.optim df
## All data    3.62       3  6
## parent      3.62       3  6
## 
## Estimated disappearance times:
##        DT50 DT90
## parent 7.25 24.1
## 
## Estimated formation fractions:
##             ff
## parent_sink  1
## 
## Parameter correlation:
## Could not estimate covariance matrix; singular system:
```


Due to the higher number of parameters, and the lower number of degrees of
freedom of the fit, the chi^2 error level is actually higher for the FOMC
model (3.6%) than for the SFO model (3.4%). Additionally, the covariance
matrix can not be obtained, indicating overparameterisation of the model.
As a consequence, no standard errors for transformed parameters nor
confidence intervals for backtransformed parameters are available.

The chi^2 error levels reported in Appendix 3 and Appendix 7 to the FOCUS
kinetics report are rounded to integer percentages and partly deviate by one
percentage point from the results calculated by mkin. The reason for
this is not known. However, mkin gives the same chi^2 error levels
as the kinfit package.

Furthermore, the calculation routines of the kinfit package have been extensively
compared to the results obtained by the KinGUI software, as documented in the
kinfit package vignette. KinGUI is a widely used standard package in this field.
Therefore, the reason for the difference was not investigated further.

## Laboratory Data L2

The following code defines example dataset L2 from the FOCUS kinetics
report, p. 287


```r
FOCUS_2006_L2 = data.frame(t = rep(c(0, 1, 3, 7, 14, 28), each = 2), parent = c(96.1, 
    91.8, 41.4, 38.7, 19.3, 22.3, 4.6, 4.6, 2.6, 1.2, 0.3, 0.6))
FOCUS_2006_L2_mkin <- mkin_wide_to_long(FOCUS_2006_L2)
```


Again, the SFO model is fitted and a summary is obtained.


```r
m.L2.SFO <- mkinfit(SFO, FOCUS_2006_L2_mkin, quiet = TRUE)
summary(m.L2.SFO)
```

```
## mkin version:    0.9.25 
## R version:       3.0.2 
## Date of fit:     Sun Nov 17 15:02:55 2013 
## Date of summary: Sun Nov 17 15:02:55 2013 
## 
## Equations:
## [1] d_parent = - k_parent_sink * parent
## 
## Method used for solution of differential equation system:
## analytical 
## 
## Weighting: none
## 
## Starting values for optimised parameters:
##               value   type transformed
## parent_0      100.0  state     100.000
## k_parent_sink   0.1 deparm      -2.303
## 
## Fixed parameter values:
## None
## 
## Optimised, transformed parameters:
##               Estimate Std. Error  Lower  Upper
## parent_0        91.500      3.810 83.000 99.900
## k_parent_sink   -0.411      0.107 -0.651 -0.172
## 
## Backtransformed parameters:
##               Estimate  Lower  Upper
## parent_0        91.500 83.000 99.900
## k_parent_sink    0.663  0.522  0.842
## 
## Residual standard error: 5.51 on 10 degrees of freedom
## 
## Chi2 error levels in percent:
##          err.min n.optim df
## All data    14.4       2  4
## parent      14.4       2  4
## 
## Estimated disappearance times:
##        DT50 DT90
## parent 1.05 3.47
## 
## Estimated formation fractions:
##             ff
## parent_sink  1
## 
## Parameter correlation:
##               parent_0 k_parent_sink
## parent_0          1.00          0.43
## k_parent_sink     0.43          1.00
## 
## Data:
##  time variable observed predicted residual
##     0   parent     96.1  9.15e+01    4.634
##     0   parent     91.8  9.15e+01    0.334
##     1   parent     41.4  4.71e+01   -5.740
##     1   parent     38.7  4.71e+01   -8.440
##     3   parent     19.3  1.25e+01    6.779
##     3   parent     22.3  1.25e+01    9.779
##     7   parent      4.6  8.83e-01    3.717
##     7   parent      4.6  8.83e-01    3.717
##    14   parent      2.6  8.53e-03    2.591
##    14   parent      1.2  8.53e-03    1.191
##    28   parent      0.3  7.96e-07    0.300
##    28   parent      0.6  7.96e-07    0.600
```


The chi^2 error level of 14% suggests that the model does not fit very well.
This is also obvious from the plots of the fit and the residuals.


```r
par(mfrow = c(2, 1))
plot(m.L2.SFO)
mkinresplot(m.L2.SFO)
```

![plot of chunk unnamed-chunk-9](figure/unnamed-chunk-9.png) 


In the FOCUS kinetics report, it is stated that there is no apparent systematic
error observed from the residual plot up to the measured DT90 (approximately at
day 5), and there is an underestimation beyond that point.

We may add that it is difficult to judge the random nature of the residuals just 
from the three samplings at days 0, 1 and 3. Also, it is not clear _a
priori_ why a consistent underestimation after the approximate DT90 should be
irrelevant. However, this can be rationalised by the fact that the FOCUS fate
models generally only implement SFO kinetics.

For comparison, the FOMC model is fitted as well, and the chi^2 error level
is checked.


```r
m.L2.FOMC <- mkinfit(FOMC, FOCUS_2006_L2_mkin, quiet = TRUE)
par(mfrow = c(2, 1))
plot(m.L2.FOMC)
mkinresplot(m.L2.FOMC)
```

![plot of chunk unnamed-chunk-10](figure/unnamed-chunk-10.png) 

```r
summary(m.L2.FOMC, data = FALSE)
```

```
## mkin version:    0.9.25 
## R version:       3.0.2 
## Date of fit:     Sun Nov 17 15:02:56 2013 
## Date of summary: Sun Nov 17 15:02:56 2013 
## 
## Equations:
## [1] d_parent = - (alpha/beta) * ((time/beta) + 1)^-1 * parent
## 
## Method used for solution of differential equation system:
## analytical 
## 
## Weighting: none
## 
## Starting values for optimised parameters:
##          value   type transformed
## parent_0   100  state     100.000
## alpha        1 deparm       0.000
## beta        10 deparm       2.303
## 
## Fixed parameter values:
## None
## 
## Optimised, transformed parameters:
##          Estimate Std. Error  Lower  Upper
## parent_0   93.800      1.860 89.600 98.000
## alpha       0.318      0.187 -0.104  0.740
## beta        0.210      0.294 -0.456  0.876
## 
## Backtransformed parameters:
##          Estimate  Lower Upper
## parent_0    93.80 89.600  98.0
## alpha        1.37  0.901   2.1
## beta         1.23  0.634   2.4
## 
## Residual standard error: 2.63 on 9 degrees of freedom
## 
## Chi2 error levels in percent:
##          err.min n.optim df
## All data     6.2       3  3
## parent       6.2       3  3
## 
## Estimated disappearance times:
##         DT50 DT90
## parent 0.809 5.36
## 
## Estimated formation fractions:
##             ff
## parent_sink  1
## 
## Parameter correlation:
##          parent_0   alpha   beta
## parent_0   1.0000 -0.0955 -0.186
## alpha     -0.0955  1.0000  0.976
## beta      -0.1863  0.9757  1.000
```


The error level at which the chi^2 test passes is much lower in this case.
Therefore, the FOMC model provides a better description of the data, as less
experimental error has to be assumed in order to explain the data.

Fitting the four parameter DFOP model further reduces the chi^2 error level. 


```r
m.L2.DFOP <- mkinfit(DFOP, FOCUS_2006_L2_mkin, quiet = TRUE)
plot(m.L2.DFOP)
```

![plot of chunk unnamed-chunk-11](figure/unnamed-chunk-11.png) 


Here, the default starting parameters for the DFOP model obviously do not lead
to a reasonable solution. Therefore the fit is repeated with different starting
parameters.


```r
m.L2.DFOP <- mkinfit(DFOP, FOCUS_2006_L2_mkin, parms.ini = c(k1 = 1, k2 = 0.01, 
    g = 0.8), quiet = TRUE)
plot(m.L2.DFOP)
```

![plot of chunk unnamed-chunk-12](figure/unnamed-chunk-12.png) 

```r
summary(m.L2.DFOP, data = FALSE)
```

```
## mkin version:    0.9.25 
## R version:       3.0.2 
## Date of fit:     Sun Nov 17 15:02:57 2013 
## Date of summary: Sun Nov 17 15:02:57 2013 
## 
## Equations:
## [1] d_parent = - ((k1 * g * exp(-k1 * time) + k2 * (1 - g) * exp(-k2 * time)) / (g * exp(-k1 * time) + (1 - g) * exp(-k2 * time))) * parent
## 
## Method used for solution of differential equation system:
## analytical 
## 
## Weighting: none
## 
## Starting values for optimised parameters:
##          value   type transformed
## parent_0 1e+02  state    100.0000
## k1       1e+00 deparm      0.0000
## k2       1e-02 deparm     -4.6052
## g        8e-01 deparm      0.9803
## 
## Fixed parameter values:
## None
## 
## Optimised, transformed parameters:
##          Estimate Std. Error Lower Upper
## parent_0   93.900         NA    NA    NA
## k1          4.960         NA    NA    NA
## k2         -1.090         NA    NA    NA
## g          -0.282         NA    NA    NA
## 
## Backtransformed parameters:
##          Estimate Lower Upper
## parent_0   93.900    NA    NA
## k1        142.000    NA    NA
## k2          0.337    NA    NA
## g           0.402    NA    NA
## 
## Residual standard error: 1.73 on 8 degrees of freedom
## 
## Chi2 error levels in percent:
##          err.min n.optim df
## All data    2.53       4  2
## parent      2.53       4  2
## 
## Estimated disappearance times:
##        DT50 DT90
## parent   NA   NA
## 
## Estimated formation fractions:
##             ff
## parent_sink  1
## 
## Parameter correlation:
## Could not estimate covariance matrix; singular system:
```


Here, the DFOP model is clearly the best-fit model for dataset L2 based on the 
chi^2 error level criterion. However, the failure to calculate the covariance
matrix indicates that the parameter estimates correlate excessively. Therefore,
the FOMC model may be preferred for this dataset.

## Laboratory Data L3

The following code defines example dataset L3 from the FOCUS kinetics report,
p. 290.


```r
FOCUS_2006_L3 = data.frame(t = c(0, 3, 7, 14, 30, 60, 91, 120), parent = c(97.8, 
    60, 51, 43, 35, 22, 15, 12))
FOCUS_2006_L3_mkin <- mkin_wide_to_long(FOCUS_2006_L3)
```


SFO model, summary and plot:


```r
m.L3.SFO <- mkinfit(SFO, FOCUS_2006_L3_mkin, quiet = TRUE)
plot(m.L3.SFO)
```

![plot of chunk unnamed-chunk-14](figure/unnamed-chunk-14.png) 

```r
summary(m.L3.SFO)
```

```
## mkin version:    0.9.25 
## R version:       3.0.2 
## Date of fit:     Sun Nov 17 15:02:57 2013 
## Date of summary: Sun Nov 17 15:02:57 2013 
## 
## Equations:
## [1] d_parent = - k_parent_sink * parent
## 
## Method used for solution of differential equation system:
## analytical 
## 
## Weighting: none
## 
## Starting values for optimised parameters:
##               value   type transformed
## parent_0      100.0  state     100.000
## k_parent_sink   0.1 deparm      -2.303
## 
## Fixed parameter values:
## None
## 
## Optimised, transformed parameters:
##               Estimate Std. Error Lower Upper
## parent_0         74.90      8.460 54.20 95.60
## k_parent_sink    -3.68      0.326 -4.48 -2.88
## 
## Backtransformed parameters:
##               Estimate   Lower   Upper
## parent_0       74.9000 54.2000 95.6000
## k_parent_sink   0.0253  0.0114  0.0561
## 
## Residual standard error: 12.9 on 6 degrees of freedom
## 
## Chi2 error levels in percent:
##          err.min n.optim df
## All data    21.2       2  6
## parent      21.2       2  6
## 
## Estimated disappearance times:
##        DT50 DT90
## parent 27.4 91.1
## 
## Estimated formation fractions:
##             ff
## parent_sink  1
## 
## Parameter correlation:
##               parent_0 k_parent_sink
## parent_0         1.000         0.548
## k_parent_sink    0.548         1.000
## 
## Data:
##  time variable observed predicted residual
##     0   parent     97.8     74.87  22.9273
##     3   parent     60.0     69.41  -9.4065
##     7   parent     51.0     62.73 -11.7340
##    14   parent     43.0     52.56  -9.5634
##    30   parent     35.0     35.08  -0.0828
##    60   parent     22.0     16.44   5.5614
##    91   parent     15.0      7.51   7.4896
##   120   parent     12.0      3.61   8.3908
```


The chi^2 error level of 21% as well as the plot suggest that the model
does not fit very well. 

The FOMC model performs better:


```r
m.L3.FOMC <- mkinfit(FOMC, FOCUS_2006_L3_mkin, quiet = TRUE)
plot(m.L3.FOMC)
```

![plot of chunk unnamed-chunk-15](figure/unnamed-chunk-15.png) 

```r
summary(m.L3.FOMC, data = FALSE)
```

```
## mkin version:    0.9.25 
## R version:       3.0.2 
## Date of fit:     Sun Nov 17 15:02:57 2013 
## Date of summary: Sun Nov 17 15:02:57 2013 
## 
## Equations:
## [1] d_parent = - (alpha/beta) * ((time/beta) + 1)^-1 * parent
## 
## Method used for solution of differential equation system:
## analytical 
## 
## Weighting: none
## 
## Starting values for optimised parameters:
##          value   type transformed
## parent_0   100  state     100.000
## alpha        1 deparm       0.000
## beta        10 deparm       2.303
## 
## Fixed parameter values:
## None
## 
## Optimised, transformed parameters:
##          Estimate Std. Error Lower   Upper
## parent_0   97.000      4.550  85.3 109.000
## alpha      -0.862      0.170  -1.3  -0.424
## beta        0.619      0.474  -0.6   1.840
## 
## Backtransformed parameters:
##          Estimate  Lower   Upper
## parent_0   97.000 85.300 109.000
## alpha       0.422  0.273   0.655
## beta        1.860  0.549   6.290
## 
## Residual standard error: 4.57 on 5 degrees of freedom
## 
## Chi2 error levels in percent:
##          err.min n.optim df
## All data    7.32       3  5
## parent      7.32       3  5
## 
## Estimated disappearance times:
##        DT50 DT90
## parent 7.73  431
## 
## Estimated formation fractions:
##             ff
## parent_sink  1
## 
## Parameter correlation:
##          parent_0  alpha   beta
## parent_0    1.000 -0.151 -0.427
## alpha      -0.151  1.000  0.911
## beta       -0.427  0.911  1.000
```


The error level at which the chi^2 test passes is 7% in this case.

Fitting the four parameter DFOP model further reduces the chi^2 error level
considerably:


```r
m.L3.DFOP <- mkinfit(DFOP, FOCUS_2006_L3_mkin, quiet = TRUE)
plot(m.L3.DFOP)
```

![plot of chunk unnamed-chunk-16](figure/unnamed-chunk-16.png) 

```r
summary(m.L3.DFOP, data = FALSE)
```

```
## mkin version:    0.9.25 
## R version:       3.0.2 
## Date of fit:     Sun Nov 17 15:02:58 2013 
## Date of summary: Sun Nov 17 15:02:58 2013 
## 
## Equations:
## [1] d_parent = - ((k1 * g * exp(-k1 * time) + k2 * (1 - g) * exp(-k2 * time)) / (g * exp(-k1 * time) + (1 - g) * exp(-k2 * time))) * parent
## 
## Method used for solution of differential equation system:
## analytical 
## 
## Weighting: none
## 
## Starting values for optimised parameters:
##          value   type transformed
## parent_0 1e+02  state     100.000
## k1       1e-01 deparm      -2.303
## k2       1e-02 deparm      -4.605
## g        5e-01 deparm       0.000
## 
## Fixed parameter values:
## None
## 
## Optimised, transformed parameters:
##          Estimate Std. Error  Lower    Upper
## parent_0   97.700     1.4400 93.800 102.0000
## k1         -0.661     0.1330 -1.030  -0.2910
## k2         -4.290     0.0590 -4.450  -4.1200
## g          -0.123     0.0512 -0.265   0.0193
## 
## Backtransformed parameters:
##          Estimate   Lower    Upper
## parent_0  97.7000 93.8000 102.0000
## k1         0.5160  0.3560   0.7480
## k2         0.0138  0.0117   0.0162
## g          0.4570  0.4070   0.5070
## 
## Residual standard error: 1.44 on 4 degrees of freedom
## 
## Chi2 error levels in percent:
##          err.min n.optim df
## All data    2.23       4  4
## parent      2.23       4  4
## 
## Estimated disappearance times:
##        DT50 DT90
## parent 7.46  123
## 
## Estimated formation fractions:
##             ff
## parent_sink  1
## 
## Parameter correlation:
##          parent_0     k1      k2      g
## parent_0   1.0000  0.164  0.0131  0.425
## k1         0.1640  1.000  0.4648 -0.553
## k2         0.0131  0.465  1.0000 -0.663
## g          0.4253 -0.553 -0.6631  1.000
```


Here, a look to the model plot, the confidence intervals of the parameters 
and the correlation matrix suggest that the parameter estimates are reliable, and
the DFOP model can be used as the best-fit model based on the chi^2 error
level criterion for laboratory data L3.

## Laboratory Data L4

The following code defines example dataset L4 from the FOCUS kinetics
report, p. 293


```r
FOCUS_2006_L4 = data.frame(t = c(0, 3, 7, 14, 30, 60, 91, 120), parent = c(96.6, 
    96.3, 94.3, 88.8, 74.9, 59.9, 53.5, 49))
FOCUS_2006_L4_mkin <- mkin_wide_to_long(FOCUS_2006_L4)
```


SFO model, summary and plot:


```r
m.L4.SFO <- mkinfit(SFO, FOCUS_2006_L4_mkin, quiet = TRUE)
plot(m.L4.SFO)
```

![plot of chunk unnamed-chunk-18](figure/unnamed-chunk-18.png) 

```r
summary(m.L4.SFO, data = FALSE)
```

```
## mkin version:    0.9.25 
## R version:       3.0.2 
## Date of fit:     Sun Nov 17 15:02:58 2013 
## Date of summary: Sun Nov 17 15:02:58 2013 
## 
## Equations:
## [1] d_parent = - k_parent_sink * parent
## 
## Method used for solution of differential equation system:
## analytical 
## 
## Weighting: none
## 
## Starting values for optimised parameters:
##               value   type transformed
## parent_0      100.0  state     100.000
## k_parent_sink   0.1 deparm      -2.303
## 
## Fixed parameter values:
## None
## 
## Optimised, transformed parameters:
##               Estimate Std. Error Lower  Upper
## parent_0         96.40       1.95 91.70 101.00
## k_parent_sink    -5.03       0.08 -5.23  -4.83
## 
## Backtransformed parameters:
##               Estimate    Lower    Upper
## parent_0      96.40000 91.70000 1.01e+02
## k_parent_sink  0.00654  0.00538 7.95e-03
## 
## Residual standard error: 3.65 on 6 degrees of freedom
## 
## Chi2 error levels in percent:
##          err.min n.optim df
## All data    3.29       2  6
## parent      3.29       2  6
## 
## Estimated disappearance times:
##        DT50 DT90
## parent  106  352
## 
## Estimated formation fractions:
##             ff
## parent_sink  1
## 
## Parameter correlation:
##               parent_0 k_parent_sink
## parent_0         1.000         0.587
## k_parent_sink    0.587         1.000
```


The chi^2 error level of 3.3% as well as the plot suggest that the model
fits very well. 

The FOMC model for comparison


```r
m.L4.FOMC <- mkinfit(FOMC, FOCUS_2006_L4_mkin, quiet = TRUE)
plot(m.L4.FOMC)
```

![plot of chunk unnamed-chunk-19](figure/unnamed-chunk-19.png) 

```r
summary(m.L4.FOMC, data = FALSE)
```

```
## mkin version:    0.9.25 
## R version:       3.0.2 
## Date of fit:     Sun Nov 17 15:02:59 2013 
## Date of summary: Sun Nov 17 15:02:59 2013 
## 
## Equations:
## [1] d_parent = - (alpha/beta) * ((time/beta) + 1)^-1 * parent
## 
## Method used for solution of differential equation system:
## analytical 
## 
## Weighting: none
## 
## Starting values for optimised parameters:
##          value   type transformed
## parent_0   100  state     100.000
## alpha        1 deparm       0.000
## beta        10 deparm       2.303
## 
## Fixed parameter values:
## None
## 
## Optimised, transformed parameters:
##          Estimate Std. Error Lower   Upper
## parent_0   99.100      1.680 94.80 103.000
## alpha      -0.351      0.372 -1.31   0.607
## beta        4.170      0.564  2.73   5.620
## 
## Backtransformed parameters:
##          Estimate Lower  Upper
## parent_0   99.100 94.80 103.00
## alpha       0.704  0.27   1.83
## beta       65.000 15.30 277.00
## 
## Residual standard error: 2.31 on 5 degrees of freedom
## 
## Chi2 error levels in percent:
##          err.min n.optim df
## All data    2.03       3  5
## parent      2.03       3  5
## 
## Estimated disappearance times:
##        DT50 DT90
## parent  109 1644
## 
## Estimated formation fractions:
##             ff
## parent_sink  1
## 
## Parameter correlation:
##          parent_0  alpha   beta
## parent_0    1.000 -0.536 -0.608
## alpha      -0.536  1.000  0.991
## beta       -0.608  0.991  1.000
```


The error level at which the chi^2 test passes is slightly lower for the FOMC 
model. However, the difference appears negligible.

