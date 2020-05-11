# checks the input object shares

stratEst.simulate.check.shares <- function( shares , LCR , num_strats ){

  if( LCR ){
    shares = matrix( NA , num_strats , num_samples )
    warning("stratEst warning: Covariates specified. The input object 'shares' is ignored.");
  }else{
    if( class(shares) == "list" ){
      num_samples = length(shares)
      shares_matrix = matrix( NA , num_strats , num_samples )
      sample_names_shares = names(shares)
      null_names = is.null(sample_names_shares)
      for( i in 1:num_samples ){
        if( null_names ){
          sample_names_shares[i] = as.character(i)
          shares_vec = shares[[i]]
        }else{
          if( sample_names_shares[i] == "" ){
            sample_names_shares[i] = as.character(i)
            shares_vec = shares[[i]]
          }else{
            shares_vec = shares[[sample_names_shares[i]]]
          }
        }
        if( length(shares_vec) != num_strats ){
          stop("stratEst error: The elements of the input object 'shares' have to be numeric vectors with as many elements as there are strategies.");
        }
        if( any(is.na(shares_vec)) | class(shares_vec) != "numeric" ){
          stop("stratEst error: The elements of the input object 'shares' have to be numeric vectors. NA values are not allowed.");
        }
        shares_matrix[,i] = shares_vec
      }
      names(shares) = sample_names_shares
    }
    else if( class(shares) == "numeric" ){
      num_samples = 1
      sample_names_shares = "1"
      if( length(shares) != num_strats ){
        stop("stratEst error: The input object 'shares' has to be a numeric vector with as many elements as there are strategies.");
      }
      if( any(is.na(shares)) | class(shares) != "numeric"  ){
        stop("stratEst error: The elements of the input object 'shares' have to be numeric vectors. NA values are not allowed.");
      }
      shares_matrix = matrix( NA , num_strats , num_samples )
      shares_matrix[,1] = shares
    }else{
      stop("stratEst error: The input object 'shares' has to be a numeric vector or a list of numeric vectors with as many elements as there are strategies.");
    }
    shares = shares_matrix
  }

  stratEst.simulate.check.shares.return <- list("shares" = shares , "names_samples" = sample_names_shares , "num_samples" = num_samples )

  return( stratEst.simulate.check.shares.return )
}
