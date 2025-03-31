#' Class stratEst.check
#' @param x object to be tested.
#' @description Checks if an object is of class \code{stratEst.check}.
#' @details Objects of class \code{stratEst.check} are returned by the function \code{stratEst.check()} of package \code{stratEst}.
#' @return \code{is.stratEst.check} returns \code{TRUE} if its argument is a \code{stratEst.check} object (that is, has "stratEst.check" among its classes) and \code{FALSE} otherwise.
#' @export
is.stratEst.check <- function( x ){
  inherits(x, "stratEst.check")
}
