# R function stratEst

#' @useDynLib stratEst
#' @importFrom Rcpp sourceCpp
#' @export

stratEst <- function( data, strategies, shares, covariates, cluster, response = "mixed", r.responses = "no", r.trembles = "global", select = "no", crit = "BIC", se = "yes", outer.runs = 10, outer.tol = 0, outer.max = 1000, inner.runs = 100, inner.tol = 0, inner.max = 10, lcr.runs = 1000, lcr.tol = 0, lcr.max = 1000, bs.samples = 1000, newton.stepsize = 1 , penalty = 0 ){
  if( missing(data) ) {
    stop("data not supplied")
  }
  if( ncol(data) != 5 ) {
    stop("data must have five columns (id,game,round,input,output)")
  }
  if( missing(strategies) ) {
    stop("strategies not supplied")
  }
  if( is.matrix(strategies) == F ){
    k = strategies
    n_inputs = length( unique( data[,4] ) )
    n_outputs = sum( unique( data[,5] ) != 0 )
    strat_states = rep(c(1:n_inputs),k)
    response_par = matrix(NA,n_inputs*k,n_outputs)
    transition_mat = rep(1,n_inputs*k) %*% t.default(c(1:n_inputs))
    strategies <- cbind(strat_states,response_par,transition_mat)
  }
  if( missing(shares) ) {
    k = sum( as.numeric( strategies[,1] == 1 ) )
    shares = rep( NA , k )
  }
  if( missing(covariates) ) {
    covariates = matrix(0,1,1)
  }
  if( missing(cluster) ) {
    cluster = matrix(0,1,1)
  }
  if ( outer.runs < 0 | outer.runs%%1 != 0 ){
    stop("Number of outer runs must be a positive integer. Default is 100.");
  }
  if ( inner.runs < 0 | inner.runs%%1 != 0 ){
    stop("Number of inner runs must be a positive integer. Default is 100.");
  }
  if ( outer.max < 0  | outer.max%%1 != 0){
    stop("Number of outer max evaluations must be a positive integer. Default is 1000.");
  }
  if ( inner.max < 0 | inner.max%%1 != 0 ){
    stop("Number of inner max evaluations must be a positive integer. Default is 100.");
  }
  if ( response != "mixed" & response != "pure" ){
    stop("response has to be one of the following: \"mixed\" or \"pure\". Default is \"mixed\".");
  }
  if ( r.responses != "no" & r.responses != "strategies" & r.responses != "states" & r.responses != "global"  ){
    stop("r_responses has to be one of the following: \"no\", \"strategies\", \"states\" or \"global\". Default is \"no\".");
  }
  if ( r.trembles != "no" & r.trembles != "strategies" & r.trembles != "states" & r.trembles != "global"  ){
    stop("r_trembles has to be one of the following: \"no\", \"strategies\", \"states\" or \"global\". Default is \"no\".");
  }
  if ( select != "no" & select != "strategies" & select != "responses"  & select != "trembles" & select != "both" & select != "all" ){
    stop("select has to be one of the following: \"no\", \"strategies\", \"responses\", \"trembles\", \"both\", or \"all\". Default is \"no\".");
  }
  if ( ( select == "strategies" | select == "all" ) && k == 1 ){
    stop("strategies cannot be selected if there is only one strategy.");
  }
  if ( crit != "AIC" & crit != "BIC" & crit != "ICL" ){
    stop("crit has to be one of the following: \"AIC\", \"BIC\", or \"ICL\". Default is \"BIC\".");
  }

  stratEst.output <- stratEst_cpp( data, strategies, shares, covariates, cluster, response, r.responses, r.trembles, select, crit, se, outer.runs, outer.tol, outer.max, inner.runs, inner.tol, inner.max, lcr.runs, lcr.tol, lcr.max, bs.samples, newton.stepsize, penalty )
  return(stratEst.output)
}
