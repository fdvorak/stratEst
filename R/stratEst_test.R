#' Runs t-tests if model parameters differ from user defined values
#' @useDynLib stratEst,.registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @param model A fitted model returned by the estimation function\code{stratEst()}. An object of class \code{stratEst}.
#' @param par A character vector indicating the type parameters to be tested. The default is a vector with all possible elements: "shares", "probs", "trembles", and "coefficients".
#' @param values A vector of numeric values the parameter are compared to. Default is zero.
#' @param alternative A character string specifying the alternative hypothesis, must be one of "two.sided" (default), "greater" or "smaller". You can specify just the initial letter.
#' @param digits An integer which specifies the number of digits of the results.
#' @export
#' @return A data.frame with one row for each estimated parameter and 5 columns:
#' \item{estimate}{The parameter estimate.}
#' \item{diff}{The difference between the estimated parameter and the numeric value (if values unequal to zero are supplied).}
#' \item{std.error}{The standard error of the parameter estimate.}
#' \item{t-value}{The t-statistic.}
#' \item{res.degrees}{The residual degrees of freedom.}
#' \item{Pr(>|t|)}{The p-value of the test.}
#' @details The test function retrieves the parameter estimates and standard errors from the model and performs a two-sided t-test.
#' @examples
#' ## Nash equilibrium strategy of rock-paper-scissors
#'
#' @export
stratEst.test <- function( model, par = c("shares","probs","trembles","coefficients"), values = 0, alternative = "two.sided", digits = 4 ){

  # checks
  if( "stratEst.model" %in% class(model) == F ){
    stop("stratEst.test error: The object passed to the argument 'model' must be of class 'stratEst.model'.")
  }
  if( "character" %in% class(par) == F ){
    stop("stratEst.test error: The object passed to the argument 'par' must be of a character string or vector.")
  }else{
    for( i in 1:length(par) ){
      if( par[i] %in% c("shares","probs","trembles","coefficients") == F ){
        stop("stratEst.test error: The object passed to the argument 'par' can only contain the following chracters: 'shares','probs','trembles','coefficients'.")
      }
    }
  }
  if( "numeric" %in% class(values) == F ){
    stop("stratEst.test error: The object passed to the argument 'values' must be numeric.")
  }
  if( "numeric" %in% class(values) == F | length(digits) != 1 ){
    stop("stratEst.test error: The object passed to the argument 'digits' must be a positive integer.")
  }

  par_matrix <- NULL
  est <- NULL
  se <- NULL
  row_names <- NULL

  if( length(model$shares.par > 0) & is.null(model$coefficients) & "shares" %in% par ){
    est <- c(est,model$shares.par)
    ses <- model$shares.se
    ses[ ses == 0 ] = NA
    se <- c(se,ses)
    row_names <- c(row_names,paste("shares.par.",as.character(seq(1,length(model$shares.par),by = 1)),sep=""))
  }
  if( is.null(model$coefficients.par) == F  & "coefficients" %in% par ){
    est <- c(est,model$coefficients.par)
    ses <- model$coefficients.se
    ses[ ses == 0 ] = NA
    se <- c(se,ses)
    row_names <- c(row_names,paste("coefficients.par.",as.character(seq(1,length(model$coefficients.par),by = 1)),sep=""))
  }
  if( length(model$probs.par > 0) & length(model$probs.se > 0)  & "probs" %in% par ){
    est <- c(est,model$probs.par)
    ses <- model$probs.se
    ses[ ses == 0 ] = NA
    se <- c(se,ses)
    row_names <- c(row_names,paste("probs.par.",as.character(seq(1,length(model$probs.par),by = 1)),sep=""))
  }
  if( length(model$trembles.par > 0)  & "trembles" %in% par ){
    est <- c(est,model$trembles.par)
    ses <- model$trembles.se
    ses[ ses == 0 ] = NA
    se <- c(se,ses)
    row_names <- c(row_names,paste("trembles.par.",as.character(seq(1,length(model$trembles.par),by = 1)),sep=""))
  }
    diff <- est-values
    z <- diff/se
    # p-value
    if( alternative == "two.sided" ){
      p <- 2*stats::pt( abs(z) , model$res.degrees , lower = F )
    }else if( alternative == "greater" ){
      p <- stats::pt( z , model$res.degrees , lower = F )
    }else if( alternative == "smaller" ){
      p <- stats::pt( z , model$res.degrees , lower = T )
    }else{
      stop("stratEst error: The argument 'alternative' must be one of the following chracter strings: 'two.sided','greater' or 'smaller'.")
    }
    if( any( values != 0 ) ){
      test_matrix = cbind( est , diff , se , z , rep(model$res.degrees, length(est)) , p )
      colnames(test_matrix) <- c("estimate","diff","std.error","t-value","df","Pr(>|t|)")
    }else{
      test_matrix = cbind( est , se , z , rep(model$res.degrees, length(est)) , p )
      colnames(test_matrix) <- c("estimate","std.error","t-value","df","Pr(>|t|)")
    }
    rownames(test_matrix) <- row_names
    par_data <- data.frame(round(test_matrix,digits))

  # colnames
  if( any( values != 0 ) ){
    colnames(par_data) <- c("estimate","diff","std.error","t-value","df","Pr(>|t|)")
  }else{
    colnames(par_data) <- c("estimate","std.error","t-value","df","Pr(>|t|)")
  }

  return(par_data)

}

.onUnload <- function (libpath) { library.dynam.unload("stratEst", libpath)}
