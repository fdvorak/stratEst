
<!-- README.md is generated from README.Rmd. Please edit that file -->

# stratEst

The stratEst package implements variants of the strategy estimation method (Dal Bo & Frechette, 2011; Breitmoser, 2015; Dvorak & Fehrler, 2018). Strategies are supplied by the user in the form of
deterministic finite-state automata. The package uses the EM algorithm
(Dempster, 1977) and the Newton-Raphson method to obtain
maximum-likelihood estimates of the population shares and choice
parameters of the strategies. The number and the complexity of
strategies can be restricted by the user or selected based on
information criteria. The package also features an extension of strategy
estimation in the spirit of latent class regression to assess the
effects of covariates on strategy use.

## Installation

To install the stratEst package from CRAN use:

``` r
install.packages("stratEst")
```

Ti install the development version of stratEst from github use:

``` r
install.packages("devtools")
devtools::install_github("fdvorak/stratEst")
```

## Example

This example shows how to replicate the strategy frequency estimation of Dal Bo and Frechette (2011). The following code reproduces the results of Table 7, page 424 of the paper.

``` r
library(stratEst)
model.DF2011 <- stratEst(data.DF2011,strategies.DF2011,sample.id="treatment" )
summary(model.DF2011)
```

## References

  - Breitmoser, Y. (2015): Cooperation, but no reciprocity: Individual
    strategies in the repeated prisoner’s dilemma, American Economic
    Review, 105, 2882-2910. \doi 10.1257/aer.20130675
  - Dal Bo, P. and G. R. Frechette (2011): The evolution of cooperation
    in infinitely repeated games: Experimental evidence, American
    Economic Review, 101, 411-429. \doi 10.1257/aer.101.1.411
  - Dempster, A., N. Laird, and D. B. Rubin (1977): Maximum likelihood
    from incomplete data via the EM algorithm," Journal of the Royal
    Statistical Society Series B, 39, 1-38.
