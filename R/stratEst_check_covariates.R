# checks the input object covariates

stratEst.check.covariates <- function( data , covariates ){

  covariate_mat <- NULL
  for( i in 1:length(covariates) ){
    if( covariates[i] %in% colnames(data) ){
      covariate_mat <- cbind(covariate_mat,data[,covariates[i]])
    }
    else{
      stop(paste("stratEst error: The data does not contain the variable '",covariates[i],"' "," specified as covariate.",sep=""))
    }
  }

  return( covariate_mat )

}
