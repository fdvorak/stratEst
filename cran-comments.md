---
title: "cran-comments"
author: "Fabian Dvorak"
date: "July 06 2020"
output: html_document
---

Thanks for your time. This is patch release 1.0.1. This version:

* fixes clang-UBSAN and gcc UBSAN errors
* makes sure RcppArmadillo package version >= 0.9.900.0.0 is installed
* fixes bug in strategy checking function of stratEst.simulate()
* corrects number of free parameters in the estimation of pure response probabilities
* corrects object shares.par with more than one sample

## Test environments
* windows 10 pro install, R 4.0.0 (local)
* windows x86_64-w64-mingw32, r-oldrel, 3.6.3 (with win-builder)
* windows x86_64-w64-mingw32, r-devel, 2020-06-19 r78718 (with win-builder)
* windows x86_64-w64-mingw32, r-release, 4.0.2 (with win-builder)
* mac OS X	10.13.6, r-oldrel, 3.6.3 (with travis-ci)
* mac OS X	10.13.6, r-release, 4.0.0 (with travis-ci)
* ubuntu 14.04.5 LTS , r-oldrel, 3.6.3 (with travis-ci)
* ubuntu 14.04.5 LTS , r-devel, 4.0.1 (with travis-ci)
* ubuntu 14.04.5 LTS , r-release, 4.0.0 (with travis-ci)

## R CMD check results
* windows x86_64-w64-mingw32, r-oldrel, 3.6.3 (with win-builder)
* windows x86_64-w64-mingw32, r-devel, 2020-06-19 r78718 (with win-builder)
* windows x86_64-w64-mingw32, r-release, 4.0.2 (with win-builder)
* mac OS X	10.13.6, r-oldrel, 3.6.3 (with travis-ci)
* mac OS X	10.13.6, r-release, 4.0.0 (with travis-ci)
* ubuntu 14.04.5 LTS , r-oldrel, 3.6.3 (with travis-ci)
* ubuntu 14.04.5 LTS , r-devel, 4.0.1 (with travis-ci)
* ubuntu 14.04.5 LTS , r-release, 4.0.0 (with travis-ci)

0 errors | 0 warnings | 0 notes

* windows 10 pro install, R 4.0.0 (local)

0 errors √ | 0 warnings √ | 1 note x

> checking for non-standard things in the check directory ... NOTE
  Found the following files/directories:
    'examples_i386' 'examples_x64' 'stratEst-Ex_i386.Rout'
    'stratEst-Ex_x64.Rout' 'tests_i386' 'tests_x64'
    
I couldn't fix this by adding these folders to Rbuildignore.

## Downstream dependencies
No downstream dependencies of the package.

