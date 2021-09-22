#' Round Method for stratEst.strategy
#' @param x An object of class \code{stratEst.strategy}.
#' @param digits Further arguments passed to or from other methods.
#' @return A \code{stratEst.strategy} object with rounded variable values. A data.frame with the following variables:
#' \item{prob.x}{the probability of choice \code{x}.}
#' \item{tremble}{the probability to observe a tremble.}
#' \item{tr(x)}{the deterministic state transitions of the strategy for input \code{x}.}
#' @export

round.stratEst.strategy <- function( x , digits = 0 ){
  for( i in 1:ncol(x)){
    x[,i] <- round(x[,i],digits )
  }
  return(x)
}
