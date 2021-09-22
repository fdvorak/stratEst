#' Print Method for stratEst.data
#' @param x An object of class \code{stratEst.data}.
#' @param ... Further arguments passed to or from other methods.
#' @return Prints a \code{data.frame} object that contains the data.
#' @export

print.stratEst.data <- function( x , ... ){

  c("stratEst.data", NextMethod())

}
