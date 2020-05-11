# checks other input objects

stratEst.check.other <- function( response , sample.specific , r.responses , r.trembles , select , min.strategies , crit , se , outer.runs , outer.tol , outer.max , inner.runs , inner.tol , inner.max , lcr.runs , lcr.tol , lcr.max , bs.samples , stepsize , penalty , verbose , quantiles ){

  # check response
  if ( response %in% c("mixed","pure") == F ){
    stop("stratEst error: The input object 'response' has to be one of the following: \"mixed\" or \"pure\". Default is \"mixed\".");
  }


  # check sample.specific
  specific_shares = F
  specific_responses = F
  specific_trembles = F
  specific_coefficients = F
  if( is.null(sample.specific) == F ){
    if( class( sample.specific ) != "character" ){
      stop("stratEst error: The input object 'sample.specific' has to be a character vector.");
    }
    for( i in 1:length( sample.specific ) ){
      if ( sample.specific[i] %in% c("shares","responses","trembles","coefficients") == F  ){
        stop("stratEst error: The input object 'sample.specific' should only contain the following characters: \"shares\", \"responses\", \"trembles\" or \"coefficients\".");
      }
    }
    specific_shares = ifelse( "shares" %in% sample.specific , T , F  )
    specific_responses = ifelse( "responses" %in% sample.specific , T , F  )
    specific_trembles = ifelse( "trembles" %in% sample.specific , T , F  )
    specific_coefficients = ifelse( "coefficients" %in% sample.specific , T , F  )
  }

  # check r.responses
  if ( r.responses %in% c("no","strategies","states","global") == F  ){
    stop("stratEst error: The input object 'r.responses' has to be one of the following: \"no\", \"strategies\", \"states\" or \"global\". Default is \"no\".");
  }

  # check r.trembles
  if ( r.trembles %in% c("no","strategies","states","global") == F  ){
    stop("stratEst error: The input object 'r.trembles' has to be one of the following: \"no\", \"strategies\", \"states\" or \"global\". Default is \"no\".");
  }

  # check select
  select_strategies = F
  select_responses = F
  select_trembles = F

  if( is.null(select) == F ){
    # check select
    if( class( select ) != "character" ){
      stop("stratEst error: The input object 'select' has to be a character vector.");
    }
    for( i in 1:length( select ) ){
      if ( select[i] %in% c("responses","trembles","strategies") == F  ){
        stop("stratEst error: The input object 'select' should only contain the following characters: \"strategies\", \"responses\" or \"trembles\".");
      }
      else{
        if( select[i] == "strategies" ){
          select_strategies = T
        }
        if( select[i] == "responses" ){
          select_responses = T
        }
        if( select[i] == "trembles" ){
          select_trembles = T
        }
      }
    }
  }

  # check min-strategies
  if ( min.strategies < 1 | min.strategies%%1 != 0 ){
    stop("stratEst error: The minimum number of strategies must be a positive integer. Default is 1.");
  }

  # check crit
  if ( crit %in% c("aic","bic","icl") == F ){
    stop("stratEst error: The input object 'crit' has to be one of the following: \"aic\", \"bic\", or \"icl\". Default is \"bic\".");
  }

  # check se
  if ( se %in% c("analytic","bootstrap") == F ){
    stop("stratEst error: The input object 'se' has to be one of the following: \"analytic\", or \"bootstrap\". Default is \"analytic\".");
  }

  # check outer.runs
  if ( outer.runs < 0 | outer.runs%%1 != 0 ){
    stop("stratEst error: The number of outer runs must be a positive integer. Default is 100.");
  }

  # check inner.runs
  if ( inner.runs < 0 | inner.runs%%1 != 0 ){
    stop("stratEst error: The number of inner runs must be a positive integer. Default is 100.");
  }

  # check lcr.runs
  if ( lcr.runs < 0 | lcr.runs%%1 != 0 ){
    stop("stratEst error: The number of lcr runs must be a positive integer. Default is 100.");
  }

  # check outer.max
  if ( outer.max < 0  | outer.max%%1 != 0){
    stop("stratEst error: The maximum of the number function evaluations of the outer runs must be a positive integer. Default is 1000.");
  }

  # check inner.max
  if ( inner.max < 0 | inner.max%%1 != 0 ){
    stop("stratEst error: The maximum of the number function evaluations of the inner runs must be a positive integer. Default is 100.");
  }

  # check lcr.max
  if ( lcr.max < 0 | lcr.max%%1 != 0 ){
    stop("stratEst error: The maximum of the number function evaluations of the lcr runs must be a positive integer. Default is 1000.");
  }

  # check outer.tol
  if ( outer.tol < 0 | outer.tol >=1 ){
    stop("stratEst error: The tolerance of the outer runs must be a small numeric value. Default is 0.");
  }

  # check inner.tol
  if ( inner.tol < 0 | inner.tol >=1 ){
    stop("stratEst error: The tolerance of the inner runs must be a small numeric value. Default is 0.");
  }

  # check lcr.tol
  if ( lcr.tol < 0 | lcr.tol >=1 ){
    stop("stratEst error: The tolerance of the lcr runs must be a small numeric value. Default is 0.");
  }


  # check bs.samples
  if ( bs.samples < 0  | bs.samples%%1 != 0){
    stop("stratEst error: The number of bootstrap samples specified by the argument 'bs.samples' must be a positive integer. Default is 1000.");
  }

  # check stepsize
  if ( stepsize < 0 ){
    stop("stratEst error: The newton stepsize specified by the argument 'newton.stepsize' must be a positive number. Default is 1.");
  }

  # check penalty
  if (  is.logical(penalty) == F){
    stop("stratEst error: The function argument 'penalty' must be boolean. Default is FALSE.");
  }

  # check verbose
  if (  class(verbose) != "logical"){
    stop("stratEst error: The input argument 'verbose' must be a logical value.");
  }
  else{
    print.messages = verbose[1]
    print.summary = F
  }

  # check print.summary
  if (  class(print.summary) != "logical"){
    stop("stratEst error: The input argument 'print.summary' must be a logical value.");
  }

  # check quantiles
  if (  class(quantiles) != "numeric"){
    stop("stratEst error: The input argument 'print.summary' must be a logical value.");
  }
  else{
    if( any(quantiles>1) | any(quantiles<0) ){
      stop("stratEst error: The elements of the input argument 'qunatiles' must numeric values between zero and one.");
    }
  }
  qunantile_vec <- quantiles

  stratEst.check.other.return = list( "select.strategies" = select_strategies , "select.responses" = select_responses , "select.trembles" = select_trembles, "specific.shares" = specific_shares , "specific.responses" = specific_responses , "specific.trembles" = specific_trembles, "specific.coefficients" = specific_coefficients ,  "quantile.vec" = qunantile_vec , "print.messages" = print.messages , "print.summary" = print.summary  )

  return(stratEst.check.other.return)

}
