---
title: "cran-comments"
author: "Fabian Dvorak"
date: "January 25 2019"
output: html_document
---

Thanks for your time. This is a patch release for version 0.1.1. This version:

* fixes a bug in the calculation of the score function for trembles
* fixes a bug in the selection procedure if the input object strategies is an integer
* adds the current working paper version of the project

## Test environments
* windows 10 pro install, R 3.5.2 (local)
* windows server 2012 R2 x64, R 3.5.2 patched (via appveyor)
* windows server 2008 R2 SP1, r-devel (via r hub)
* windows, r-devel (via win-builder)
* windows, r-release (via win-builder)
* mac OS X	10.13.6, R 3.5.2 (via travis-ci)
* debian linux, R-devel (via r hub)
* fedora linux, r-devel (via r hub)
* ubuntu 16.04 LTS, r-release (via r hub)
* ubuntu 14.04.5 LTS , R 3.5.2 (via travis-ci)

## R CMD check results
windows 10 pro install, R 3.5.2 (local)
windows server 2012 R2 x64, R 3.5.2 patched (via appveyor)
windows server 2008 R2 SP1, r-devel (via r hub)
windows, r-devel (via win-builder)
windows, r-release (via win-builder)
mac OS X	10.13.6, R 3.5.2 (via travis-ci)
0 errors | 0 warnings | 0 notes

debian linux, R-devel (via r hub)
fedora linux, r-devel (via r hub)
ubuntu 16.04 LTS, r-release (via r hub)
ubuntu 14.04.5 LTS , R 3.5.2 (via travis-ci)  
0 errors | 0 warnings | 1 note  
NOTE: installed size is  7.7Mb
      sub-directories of 1Mb or more:
      libs   7.0Mb
  
My understanding is that the inflation of the libs subdirectory is due to the use of Rcpp. The main functions of the stratEst package have been written in C++ using Rcpp and RcppArmadillo. These functions are essential for the functioning of the package. Without the speed up gained from those C++ functions, the package would become impractical.

## Downstream dependencies
I have also run R CMD check on downstream dependencies of the package.  

0 packages with problems

