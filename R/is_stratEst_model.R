#' Class stratEst.model
#' @param x object to be tested.
#' @description Checks if an object is of class \code{stratEst.model}.
#' @details Objects of class \code{stratEst.model} are returned by the estimation function \code{stratEst.model()} of package \code{stratEst}.
#' @return \code{is.stratEst.model} returns \code{TRUE} if its argument is a \code{stratEst.model} object (that is, has "stratEst.model" among its classes) and \code{FALSE} otherwise.
#' @export
is.stratEst.model <- function( x ){
  inherits(x, "stratEst.model")
}
