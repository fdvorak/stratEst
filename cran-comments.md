---
title: "cran-comments"
author: "Fabian Dvorak"
date: "June 23 2020"
output: html_document
---

Thanks for your time. This is package version 1.0.0. This version introduces:

* S3 classes stratEst.model, stratEst.data, stratEst.strategy, and stratEst.check
* function stratEst.data() to reshape data
* function stratEst.simulate() to simulate data
* function stratEst.model() to fit models
* function stratEst.strategy() to generate strategies
* function stratEst.check() to check fitted models
* function stratEst.test() to test estimated model parameters

## Test environments
* windows 10 pro install, R 4.0.0 (local)
* windows server 2012 R2 x64, R 3.5.2 patched (via appveyor)
* windows, r-devel (via win-builder)
* windows, r-release (via win-builder)
* mac OS X	10.13.6, R 3.5.2 (via travis-ci)
* mac OS X	10.13.6, r-oldrel (via travis-ci)
* ubuntu 14.04.5 LTS , R 3.5.2 (via travis-ci)
* ubuntu 14.04.5 LTS , r-oldrel (via travis-ci)

## R CMD check results
windows 10 pro install, R 3.5.2 (local) 
windows server 2012 R2 x64, R 3.5.2 patched (via appveyor)  
windows, r-devel (via win-builder)  
windows, r-release (via win-builder)  
mac OS X	10.13.6, R 3.5.2 (via travis-ci)  
mac OS X	10.13.6, r-oldrel (via travis-ci)   
0 errors | 0 warnings | 0 notes

ubuntu 14.04.5 LTS , R 3.5.2 (via travis-ci)   
ubuntu 14.04.5 LTS , r-oldrel (via travis-ci)   
0 errors | 0 warnings | 1 note  
NOTE: installed size is  7.7Mb
      sub-directories of 1Mb or more:
      libs   7.0Mb
  
My understanding is that the inflation of the libs subdirectory is due to the use of Rcpp. The main functions of the stratEst package have been written in C++ using Rcpp and RcppArmadillo. These functions are essential for the functioning of the package. Without the speed up gained from those C++ functions, the package would become impractical.

## Downstream dependencies
I have also run R CMD check on downstream dependencies of the package.  

0 packages with problems

