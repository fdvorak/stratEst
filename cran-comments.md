## Resubmission
This is a resubmission. This version:

* uses SHLIB_OPENMP_CFLAGS in both PKG_CXXFLAGS and PKG_LIBS

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
0 errors | 0 warnings | 1 note  
Note: Maintainer: 'Fabian Dvorak <fabian.dvorak@uni.kn>'     
New submission

win-builder (r-release)     
0 errors | 0 warnings | 1 note  
Note: Maintainer: 'Fabian Dvorak <fabian.dvorak@uni.kn>'      
New submission

ubuntu 14.04.5 LTS (on travis-ci)  
0 errors | 0 warnings | 1 note  
NOTE: installed size is  7.1Mb
      sub-directories of 1Mb or more:
      libs   6.8Mb
  
My current understanding is that the inflation of the libs subdirectory is due to the use of Rcpp. The main functions of the stratEst package have been written in C++ using Rcpp and RcppArmadillo. These functions are essential for the functioning of the package. Without the speed up gained from those C++ functions, the package would become impractical.

## Downstream dependencies
I have also run R CMD check on downstream dependencies of the package.  

No ERRORs or WARNINGs found.
