---
title: "cran-comments"
author: "Fabian Dvorak"
date: "September 15 2021"
output: html_document
---

Thanks for your time. I'm terribly sorry for this but I cannot reproduce the UBSAN errors on my windows machine (rhub::check_with_sanitizers gives me a preperror). I carefully went through the entire c++ code and tried to fix all potential problems. This is version 1.1.2:

* 2nd try to fix UBSAN problems of version 1.1.0

## Test environments
* windows 10 pro install, R 4.1.1                (local)
* windows server, x86_64, r-release, R 4.1.0     (with github actions)
* windows server, x86_64, r-oldrel, R 4.0.5      (with github actions)
* mac OS X	10.15.7, x86_64, r-release, R 4.1.0  (with github actions)
* mac OS X	10.15.7, x86_64, r-devel             (with github actions)
* mac OS X	10.15.7, x86_64, r-oldrel, R 4.0.5   (with github actions)
* ubuntu 18.04.5 LTS, x86_64, r-devel,           (with github actions)
* ubuntu 18.04.5 LTS, x86_64, r-oldrel, R 4.0.5  (with github actions)
* ubuntu 18.04.5 LTS, x86_64, r-release, R 4.1.0 (with github actions)

## R CMD check results
* windows 10 pro install, R 4.1.1                (local)
* windows server, x86_64, r-release, R 4.1.0     (with github actions)
* windows server, x86_64, r-oldrel, R 4.0.5      (with github actions)
* mac OS X	10.15.7, x86_64, r-release, R 4.1.0  (with github actions)
* mac OS X	10.15.7, x86_64, r-devel             (with github actions)
* mac OS X	10.15.7, x86_64, r-oldrel, R 4.0.5   (with github actions)

0 errors | 0 warnings | 0 notes

* ubuntu 18.04.5 LTS, x86_64, r-devel,           (with github actions)
* ubuntu 18.04.5 LTS, x86_64, r-oldrel, R 4.0.5  (with github actions)
* ubuntu 18.04.5 LTS, x86_64, r-release, R 4.1.0 (with github actions)

0 errors ✔ | 0 warnings ✔ | 1 note ✖

❯ checking installed package size ... NOTE
    installed size is 10.0Mb
    sub-directories of 1Mb or more:
      libs   8.7Mb

My understanding is that the inflation of the libs subdirectory is due to the use of Rcpp. The main functions of the stratEst package have been written in C++ using Rcpp and RcppArmadillo.


## Downstream dependencies
No downstream dependencies of the package according to revdepcheck::revdep_check.



