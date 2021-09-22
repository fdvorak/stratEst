#' Method dispatch for Generic Function Summary
#' @param object An object returned by the function\code{stratEst.check()}. An object of class \code{stratEst.check}.
#' @param ... additional arguments affecting the summary produced.
#' @return No return value, prints a summary of the model checks to the console.
#' @export

summary.stratEst.check <- function( object , ...){

  if( "stratEst.check" %in% class(object) == FALSE ){
    stop(paste("stratEst.summary error: The object ",as.character(object)," must be of class 'stratEst.check" ),sep="")
  }else{
    stratEst.check.return <- object
  }

  writeLines("model fit:")
  print(stratEst.check.return$fit)
  writeLines("")
  if( is.null(stratEst.check.return$chi.global) == FALSE ){
    writeLines("Chi^2 test of global model fit:")
    print(stratEst.check.return$chi.global)
    writeLines("")
  }
  if( is.null(stratEst.check.return$chi.local) == FALSE ){
    writeLines("Chi^2 test of local model fit:")
    print(stratEst.check.return$chi.local)
    writeLines("")
  }

}
