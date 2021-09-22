#' Print Method for stratEst.check
#' @param x An object of class \code{stratEst.check}.
#' @param ... Further arguments passed to or from other methods.
#' @return Prints a \code{matrix} that contains the log-likelihood, the number of free model parameters, and the values of three information criteria in columns.
#' @export
print.stratEst.check <- function( x , ... ){

  c("stratEst.check", NextMethod())

}
