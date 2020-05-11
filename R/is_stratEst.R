#' Class stratEst
#' @param x object to be tested.
#' @description Checks if an object is of class \code{stratEst}.
#' @details Objects of class \code{stratEst} are returned by the estimation function \code{stratEst()} of package \code{stratEst}.
#' @export
is.stratEst <- function( x ){
  inherits(x, "stratEst")
}
