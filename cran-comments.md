---
title: "cran-comments"
author: "Fabian Dvorak"
date: "December 16, 2018"
output: html_document
---

Thanks for your time. This is a patch release for version 0.1.0. This version:

* fixes the overload ambiguity of the log() function in C++ code (Error on Solaris)
* reports solutions with tremble parameters bigger than 0.5

## Test environments
* local windows 10 pro install, R 3.5.1
* mac OS X	10.13.6 (on travis-ci), R 3.5.0
* win-builder (r-devel and r-release)
* ubuntu 14.04.5 LTS (on travis-ci), 3.5.1

## R CMD check results
local windows 10 pro install  
0 errors | 0 warnings | 0 notes

mac OS X 10.13.6 (on travis-ci)  
0 errors | 0 warnings | 0 notes

win-builder (r-devel)   
0 errors | 0 warnings | 0 notes  

win-builder (r-release)     
0 errors | 0 warnings | 0 notes  

ubuntu 14.04.5 LTS (on travis-ci)  
0 errors | 0 warnings | 1 note  
NOTE: installed size is  7.1Mb
      sub-directories of 1Mb or more:
      libs   6.8Mb
  
My current understanding is that the inflation of the libs subdirectory is due to the use of Rcpp. The main functions of the stratEst package have been written in C++ using Rcpp and RcppArmadillo. These functions are essential for the functioning of the package. Without the speed up gained from those C++ functions, the package would become impractical.

## Downstream dependencies
I have also run R CMD check on downstream dependencies of the package.  

0 packages with problems

