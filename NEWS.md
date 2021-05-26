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
