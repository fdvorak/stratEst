#' Class stratEst.data
#' @param x object to be tested.
#' @description Checks if an object is of class \code{stratEst.data}.
#' @details Objects of class \code{stratEst.data} are returned by the functions \code{stratEst.data()} and \code{stratEst.simulate()} of package \code{stratEst}.
#' @return \code{is.stratEst.data} returns \code{TRUE} if its argument is a \code{stratEst.data} object (that is, has "stratEst.data" amongst its classes) and \code{FALSE} otherwise.
#' @export
is.stratEst.data <- function( x ){
  inherits(x, "stratEst.data")
}
