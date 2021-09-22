#' Method dispatch for Generic Function summary
#' @param object An object to be summarized.
#' @param ... additional arguments affecting the result.
#' @return No return value, prints a summary of the datas to the console.
#' @export

summary.stratEst.data <- function( object , ...){

  c("stratEst.data", NextMethod())

}
