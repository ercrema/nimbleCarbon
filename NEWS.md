# Version 0.2.5 (CRAN Version)
* Fixed bug in `dExponentialGrowth()` causing memory allocation problems.
* Added alternative parametrisation for logistic growth model (`dLogisticGrowth2()`)
* Updated plot functions to allow user defined x and y labels.
* Added option to specify line width in `modelPlot()`  

# Version 0.2.1
* NEW Distributions:
  * Trapezoidal Distribution.
  * Asymmetric Laplace Distribution
* Updated version of `modelPlot()` function.
* Fixed error in `compare.models()` following changes in WAIC storing after nimble version 0.12.x

# Version 0.1.2
* Fixes an installation error affecting Oracle Solaris 10, x86, 32 bit, R release, Oracle Developer Studio 12.6 

# Version 0.1.1
First release including the following key features:
* Custom population growth models (`dExponentialGrowth.R`, `dLogisticGrowth.R`,`dDoubleExponentialGrowth.R`, `dExponentialLogisticGrowth.R`, and `dLogisticExponentialGrowth.R`)
* A nimbleFunction for linear interpolation required for calibration (`interpLin.R`)
* Several diagnostic functions and posterior predictive checks, including agreement index (`agreementIndex.R`), correlation analyses of fitted and observed SPDs (`postPredCor.R`), simulation based comparison between observed SPD and fitted model (`postPredSPD.R`).
* Several plot functions for displaying posterior distributions (`postHPDplot.R`), growth models (`modelPlot.R`), and posterior predictive checks (`plot.spdppc.R`).
* Utility function for facilitating model comparison via WAIC (`compare.models.R`)
* IntCal20, ShCal20, and Marine20 calibration curves. 
* A vignette with sample workflows and a tutorial.
