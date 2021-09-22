#' Print Method for stratEst.strategy
#' @param x An object of class \code{stratEst.strategy}.
#' @param ... Further arguments passed to or from other methods.
#' @return No return value, prints a summary of the strategy to the console.
#' @export

print.stratEst.strategy <- function( x , ... ){
  x <- round.stratEst.strategy(x,digits=3)
  print.data.frame(x , ... )  #unlist("stratEst.strategy", NextMethod())

}
