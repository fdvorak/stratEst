#' Class stratEst.strategy
#' @param x object to be tested.
#' @description Checks if an object is of class \code{stratEst.strategy}.
#' @details Objects of class \code{stratEst.strategy} is returned by the function \code{stratEst.strategy()} of package \code{stratEst}.
#' @return \code{is.stratEst.strategy} returns \code{TRUE} if its argument is a \code{stratEst.strategy} object (that is, has "stratEst.strategy" among its classes) and \code{FALSE} otherwise.
#' @export
is.stratEst.strategy <- function( x ){
  inherits(x, "stratEst.strategy")
}
