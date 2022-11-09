---
title: "cran-comments"
author: "Fabian Dvorak"
date: "November 09 2022"
output: html_document
---

Thanks for your time. This is version 1.1.5:

* fixes bug if input is NA in later periods

## Test environments
* windows 10 pro install, R 4.2.2                (local)
* windows server, x86_64, r-release, R 4.2.2     (with github actions)
* windows server, x86_64, r-oldrel, R 4.1.3      (with github actions)
* mac OS X	10.16, x86_64, r-release, R 4.2.2    (with github actions)
* mac OS X	10.16, x86_64, r-devel 2022-11-08    (with github actions)
* mac OS X	10.16, x86_64, r-oldrel, R 4.1.3     (with github actions)
* ubuntu 18.04.6 LTS, x86_64, r-devel 2022-11-08 (with github actions)
* ubuntu 18.04.6 LTS, x86_64, r-oldrel, R 4.1.3  (with github actions)
* ubuntu 18.04.6 LTS, x86_64, r-release, R 4.2.2 (with github actions)

## R CMD check results
* windows 10 pro install, R 4.2.2                (local)
* windows server, x86_64, r-release, R 4.2.2     (with github actions)
* windows server, x86_64, r-oldrel, R 4.1.3      (with github actions)
* mac OS X	10.16, x86_64, r-release, R 4.2.2    (with github actions)
* mac OS X	10.16, x86_64, r-devel 2022-11-08    (with github actions)
* mac OS X	10.16, x86_64, r-oldrel, R 4.1.3     (with github actions)

0 errors | 0 warnings | 0 notes

* ubuntu 18.04.6 LTS, x86_64, r-devel 2022-11-08 (with github actions)
* ubuntu 18.04.6 LTS, x86_64, r-oldrel, R 4.1.3  (with github actions)
* ubuntu 18.04.6 LTS, x86_64, r-release, R 4.2.2 (with github actions)

0 errors ✔ | 0 warnings ✔ | 1 note ✖

❯ checking installed package size ... NOTE
    installed size is 10.1Mb
    sub-directories of 1Mb or more:
      libs   8.7Mb

My understanding is that the inflation of the libs subdirectory is due to the use of Rcpp. The main functions of the stratEst package have been written in C++ using Rcpp and RcppArmadillo.


## Downstream dependencies
No downstream dependencies of the package according to revdepcheck::revdep_check.



