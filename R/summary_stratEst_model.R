#' Dethod dispatch for Generic Function Summary
#' @param object An object returned by the estimation function\code{stratEst.model()}. An object of class \code{stratEst.model}.
#' @param ... additional arguments affecting the summary produced.
#' @export

summary.stratEst.model <- function( object , ...){

  if( "stratEst.model" %in% class(object) == F ){
    stop(paste("stratEst.summary error: The object ",as.character(object)," must be of class 'stratEst.model" ),sep="")
  }else{
    stratEst.return <- object
  }

  convergence_string <- ifelse( is.null(stratEst.return$convergence) , "no parameters estimated" , ifelse( max( stratEst.return$convergence[ is.na(stratEst.return$convergence) == F ] ) < 0.001 , "yes" , ifelse( any( stratEst.return$convergence[ is.na(stratEst.return$convergence) == F ] < 0.001 ) , "partial" , "no") ) )

  #stratEst.print <- function( stratEst.return  ){
  writeLines("==============================================================================================")
  writeLines("stratEst summary")
  writeLines("==============================================================================================")
  writeLines(paste("number of individuals: ",stratEst.return$num.ids,sep=""))
  writeLines(paste("number of observations: ",stratEst.return$num.obs,sep=""))
  writeLines(paste("number of model parameters: ",stratEst.return$num.par,sep=""))
  writeLines(paste("number of free model parameters: ",stratEst.return$free.par,sep=""))
  writeLines(paste("residual degrees of freedom: ", stratEst.return$res.degrees,sep=""))
  writeLines(paste("convergence: ",convergence_string,sep=""))
  writeLines("")
  writeLines("model fit:")
  writeLines(paste(rep("-",nchar("model fit:")),collapse = ""))
  writeLines(paste("log likelihood: ",round(stratEst.return$loglike,2),sep=""))
  writeLines(paste("model entropy: ",round(stratEst.return$entropy.model,2),sep=""))
  writeLines(paste("aic: ",round(stratEst.return$aic,2),sep=""))
  writeLines(paste("bic: ",round(stratEst.return$bic,2),sep=""))
  writeLines(paste("icl: ",round(stratEst.return$icl,2),sep=""))
  writeLines("")
  writeLines("shares:")
  writeLines(paste(rep("-",nchar("shares:")),collapse = ""))
  if( "list" %in% class(stratEst.return$shares) ){
    #writeLines(paste(rep("=",(sum(nchar(colnames(stratEst.return$shares[[1]])))+max(nchar(rownames(stratEst.return$shares[[1]])))+length(stratEst.return$shares[[1]]))), collapse = ""))
    print(round(do.call(rbind,stratEst.return$shares),2))
  }else{
    print(round(stratEst.return$shares,2))
  }
  writeLines("")
  if( is.null(stratEst.return$coefficients) == F ){
    writeLines("latent class coefficients:")
    writeLines(paste(rep("-",nchar("latent class coefficients:")),collapse = ""))
    print(round(stratEst.return$coefficients,3))
    writeLines("")
  }
  writeLines("strategies:")
  writeLines(paste(rep("-",nchar("strategies:")),collapse = ""))
  if( "list" %in% class(stratEst.return$strategies[[1]])  ){
    strategies_sample_list <- NULL
    strategies_print <- stratEst.return$strategies
    names_samples <- names(stratEst.return$strategies)
    for( i in 1:length(stratEst.return$strategies) ){
      strategies_sample <- do.call(rbind,strategies_print[[i]])
      row_names_strategies_sample <- rownames(strategies_sample)
      rownames(strategies_sample) =  paste( names_samples[i], "." , row_names_strategies_sample, sep="")
      strategies_sample_list <- rbind( strategies_sample_list , strategies_sample )
    }
    print(strategies_sample_list)
  }else{
#
#     for( i in 1:length(stratEst.return$strategies) ){
#       print(round(stratEst.return$strategies[[i]],3))
#     }
    print(do.call(rbind,stratEst.return$strategies))
    #print(round(do.call(rbind,stratEst.return$strategies),3))
  }
  writeLines("")
  writeLines("parameter estimates:")
  writeLines(paste(rep("-",nchar("parameter estimates:")),collapse = ""))

  par_matrix <- NULL

  if( length(stratEst.return$shares.par) > 1 & is.null(stratEst.return$coefficients)  ){
    par <- stratEst.return$shares.par
    se <- stratEst.return$shares.se
    se[ se == 0 ] = NA
    z <- abs(par/se)
    p <- stats::pt( z , stratEst.return$res.degrees , lower = F )
    share_matrix = cbind( par , stratEst.return$shares.quantiles , se , z , p )
    colnames(share_matrix) <- c("estimate",colnames(stratEst.return$shares.quantiles),"std.error","t-value","Pr(>|t|)")
    rownames(share_matrix) <- paste("shares.par.",as.character(seq(1,nrow(share_matrix),by = 1)),sep="")
    par_matrix = rbind(par_matrix,share_matrix)
  }
  if( is.null(stratEst.return$coefficients.par) == F ){
    par <- stratEst.return$coefficients.par
    se <- stratEst.return$coefficients.se
    se[ se == 0 ] = NA
    z <- abs(par/se)
    p <- stats::pt( z , stratEst.return$res.degrees , lower = F )
    coefficients_matrix = cbind( par ,  stratEst.return$coefficients.quantiles , se , z , p )
    colnames(coefficients_matrix) <- c("estimate",colnames(stratEst.return$coefficients.quantiles),"std. error","t-value","Pr(>|t|)")
    rownames(coefficients_matrix) <- paste("coefficients.par.",as.character(seq(1,nrow(coefficients_matrix),by = 1)),sep="")
    par_matrix = rbind(par_matrix,coefficients_matrix)

  }
  if( length(stratEst.return$probs.par > 0) & length(stratEst.return$probs.se > 0) ){
    par <- stratEst.return$probs.par
    se <- stratEst.return$probs.se
    se[ se == 0 ] = NA
    z <- abs(par/se)
    p <- stats::pt( z , stratEst.return$res.degrees , lower = F )
    response_matrix = cbind( par , stratEst.return$probs.quantiles , se , z , p )
    colnames(response_matrix) <- c("estimate",colnames(stratEst.return$probs.quantiles),"std.error","t-value","Pr(>|t|)")
    rownames(response_matrix) <- paste("probs.par.",as.character(seq(1,nrow(response_matrix),by = 1)),sep="")
    par_matrix = rbind(par_matrix,response_matrix)
  }
  if( length(stratEst.return$trembles.par > 0) ){
    par <- stratEst.return$trembles.par
    se <- stratEst.return$trembles.se
    se[ se == 0 ] = NA
    z <- abs(par/se)
    p <- stats::pt( z , stratEst.return$res.degrees , lower = F )
    tremble_matrix = cbind( par , stratEst.return$trembles.quantiles , se , z , p )
    colnames(tremble_matrix) <- c("estimate",colnames(stratEst.return$trembles.quantiles),"std.error","t-value","Pr(>|t|)")
    rownames(tremble_matrix) <- paste("trembles.par.",as.character(seq(1,nrow(tremble_matrix),by = 1)),sep="")
    par_matrix = rbind(par_matrix,tremble_matrix)
  }

  if( is.null(par_matrix) == F ){
    print(round(par_matrix,3))
    writeLines("")
  }

  writeLines("Please cite: Dvorak (2020). stratEst: strategy estimation in R.")
  writeLines("")


}