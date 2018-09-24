#' Estimation function for strategy estimation
#' @useDynLib stratEst,.registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @param data A matrix with five columns which contains the data for the estimation in a long format. Each row in data represents one observation of one individual. The first column of the data matrix contains an ID variable which identifies the observations of the same participant across rows. The second column contains a variable which indicates the supergame number of the observation. If the game was only played once, all entries in the second column equal one. The third column contains a variable which indicates the round number of an observation. If only one round was played, all entries in the third column equal one. The fourth columns contain the inputs and the fifth column the outputs for the current round. The fourth and fifth columns may also contain zeros. Zeros have the following implications. In the input column a zero indicates that no input was observed in the current round and strategy responses are determined by the start state. A zero in the output column indicates that the output falls into some reference category. The response probabilities for the reference category are not estimated. Instead the remaining response probabilities are the remaining probabilities after subtracting the response probabilities of all other outputs from one.
#' @param strategies A matrix where each row corresponds to one state of a strategy, starting with the start state of an automaton. The first column enumerates the states of each strategy in ascending order. A value of one in the first column indicates the begin of a new strategy with its start state. The columns after the first column contain the collection of multinomial response vectors. The number of columns for the multinomial response vectors must correspond to the number of unique non-zero outputs in data. Without a reference output - which is labeled with a zero in the output column of data - the columns specify the complete multinomial response distribution for each unique value in the output column. In this case, the response probabilities in each row must sum to one. With a reference output, the response probability for the response labeled with zero is omitted and the response probabilities in each row must sum to a value smaller or equal to one. The remaining columns of the strategies matrix define the deterministic state transitions. The number of columns must equal the number of unique non-zero inputs in the data. The numbers in the first column indicate the next state of the automaton if the input is one. The numbers in the second column indicate the next state if the input is two and so on.
#' @param shares A column vector of strategy shares. The number of elements must correspond to the number of strategies defined in the strategies matrix. Elements which are NA are estimated from the data. If the object is not supplied, a share is estimated for every strategy defined in the strategies matrix.
#' @param covariates A matrix where each row corresponds to same row in data. Hence, the covariate matrix must have as many rows as the data matrix. Observations which have the same ID in data must also have the same vector of covariates. Missing value are not allowed. If covariates are supplied, a latent class regression model is estimated.
#' @param cluster A column vector which indicates the assignment of each row in data to cluter units.
#' @param response String which can be set to \code{"pure"} or \code{"mixed"}. If set to \code{"pure"} all response probabilities estimated from the data are pure responses. If set to \code{"mixed"} all response probabilities estimated from the data are mixed responses. The default is \code{"mixed"}.
#' @param r.responses A string which can be set to \code{"no"}, \code{"strategies"}, \code{"states"} or \code{"global"}. If set to \code{"strategies"}, the estimation function estimates strategies with one strategy specific vector of responses in every state of the strategy. If set to \code{"states"}, one state specific vector of responses is estimated for each state. If set to \code{"global"}, a single vector of responses is estimated which applies in every state of each strategy. Default is \code{"no"}.
#' @param r.trembles String which can be set to \code{"no"}, \code{"strategies"}, \code{"states"} or \code{"global"}. If set to \code{"strategies"}, the estimation unction estimates strategies with one strategy specific tremble probability. If set to  \code{"states"}, one state specific tremble probability is estimated for each state. If set to \code{"global"}, a single tremble is estimated which applies in every state of each strategy. Default is \code{"global"}.
#' @param select String which can be set to \code{"no"}, \code{"strategies"}, \code{"responses"}, \code{"trembles"}, \code{"both"}, and \code{"all"}. If set to \code{"strategies"}, \code{"responses"}, \code{"trembles"}, the number of strategies, responses, trembles respectively are selected based on the selection criterion specified in option \code{"crit"}. If set to \code{"both"}, the number of responses and trembles are selected. If set to \code{"all"}, the number of strategies, responses, and trembles are selected. Note that the selection of responses and trembles occurs within the scope of the restriction set to these parameters (E.g. if \code{r.responses} is set to \code{"strategies"}, \code{select = "responses"} will select responses within each strategy). Default is \code{"no"}.
#' @param min.strategies Integer which specifies the minimum number of strategies for strategy selection. The selection procedure stops if the minimum is reached.
#' @param crit String which can be set to \code{"bic"}, \code{"aic"} or \code{"icl"}. If set to \code{"bic"}, model selection based on the Bayesian Information criterion is performed. If set to \code{"aic"}, the Akaike Information criterion is used. If set to \code{"icl"} the Integrated Classification Likelihood criterion is used. Default is \code{"bic"}.
#' @param se String which can be set to \code{"no"}, \code{"yes"} or \code{"bs"}. If set to \code{"no"}, standard errors are not reported. If set to \code{"yes"}, analytic standard errors are reported. If set to \code{"bs"}, bootstrapped standard errors are reported for responses and trembles. Default is \code{"yes"}.
#' @param outer.runs Positive integer which stets the number of outer runs of the solver. Default is 10.
#' @param outer.tol Positive number which stets the tolerance of the continuation condition of the outer runs. The iterative algorithm stops if the relative decrease of the log-likelihood is smaller than \code{outer.tol}. Default is 0.
#' @param outer.max Positive integer which stets the maximum number of iterations of the outer runs of the solver. The iterative algorithm stops if it did not converge after \code{"outer.max"} iterations. Default is 1000.
#' @param inner.runs  Positive integer which stets the number of inner runs of the solver. Default is 100.
#' @param inner.tol Positive number which stets the tolerance of the continuation condition of the inner EM runs. The iterative algorithm stops if the relative decrease of the log-likelihood is smaller than \code{inner.tol}. Default is 0.
#' @param inner.max Positive integer which stets the maximum number of iterations of the inner EM runs. The iterative algorithm stops if it did not converge after \code{inner.max} iterations. Default is 100.
#' @param lcr.runs Positive integer which stets the number of estimation runs for latent class regression. Default is 100.
#' @param lcr.tol Positive number which stets the tolerance of the continuation condition of the Latent Class Regression runs. The iterative algorithm stops if the relative decrease of the log-likelihood is smaller than \code{lcr.tol}. Default is 0.
#' @param lcr.max Positive integer which stets the maximum number of iterations of the Latent Class Regression EM runs. The iterative algorithm stops if it did not converge after \code{lcr.max} iterations. Default is 1000.
#' @param bs.samples Positive integer which sets the number of bootstrap samples drawn with replacement.
#' @param print.messages Logical, if \code{TRUE} messages are printed which illustrate the status of the estimation process.
#' @note The strategy estimation method was introduced by (Dal Bo & Frechette 2011) to estimate the relative frequency of a fixed set of pure strategies in the indefinitely repeated prisoner's dilemma. Breitmoser (2015) extended the method to the estimation of behavior strategies. The \pkg{stratEst} package uses the EM algorithm (Dempster, Laird & Rubin 1977) and the Newton-Raphson method to obtain maximum-likelihood estimates for the population shares and response parameters of a set of candidate strategies. The package builds on other software contributions of the R community. To increase speed the estimation procedures, the package uses integration of C++ and R achieved by the Rcpp package (Eddelbuettel & Francois 2011) and the open source linear algebra library for the C++ language RppArmadillo (Sanderson & Curtin 2016).
#' @return The function returns a list with the following elements.
#' \item{shares}{Column vector which contains the estimates of population shares for the strategies. The first element corresponds to the first strategy defined in the strategy matrix, the second element to corresponds to the second strategy and to on. Can be used as input object of the estimation function.}
#' \item{strategies}{Matrix which contains the strategies of the model. Can be used as input object of the of the estimation function.}
#' \item{responses}{Column vector which contains the response probabilities of the strategies. The value -1 indicates that the corresponding response could not be estimated since data does not contain observations the model assigns to the corresponding state.}
#' \item{trembles}{Column vector which contains the tremble probabilities of the strategies. The value -1 indicates that the corresponding response could not be estimated since data does not contain observations the model assigns to the corresponding state.}
#' \item{coefficients}{Column vector which contains the Latent Class Regression coefficients for strategies.}
#' \item{response.mat}{Matrix which contains the estimates of the response probabilities for the columns of the strategy matrix which represent the response probabilities.}
#' \item{tremble.mat}{Matrix which contains the estimates of the tremble probabilities for the columns of the strategy matrix which represent the response probabilities.}
#' \item{coefficient.mat}{Matrix which contains the latent class regression coefficients of strategies in columns. Note that the coefficients of the first strategy are one by definition.}
#' \item{loglike}{The log-Likelihood value corresponding to the reported estimates. Bigger values indicate a better fit of the model to the data.}
#' \item{crit.val}{The value of the selection criterion defined under \code{crit}. Bigger values indicate a better fit of the model.}
#' \item{eval}{Number of iterations of the solver. The reported number is the sum of iterations performed in the inner and the outer run which led to the reported estimates.}
#' \item{tol.val}{The tolerance value in the last iteration.}
#' \item{entropy}{Entropy of the assignments.}
#' \item{assignments}{Matrix which contains the posterior probability assignments of individuals to strategies. The rows of the matrix correspond to the ID sorted in ascending order beginning with the individual with the lowest ID. The columns correspond to the strategies, starting with the first strategy defined in the strategy matrix in column one.}
#' \item{priors}{Matrix which contains the individual prior probabilities of individuals as predicted by the covariate vectors of the individuals. The rows correspond to the ID sorted in ascending order beginning with the individual with the lowest ID. The columns correspond to the strategies, starting with the first strategy defined in the strategy matrix.}
#' \item{shares.se}{Column vector which contains the standard errors of the estimated shares. The elements correspond to the vector of estimates.}
#' \item{responses.se}{Column vector which contains the standard errors of the reported responses. The elements correspond to the vector of estimates.}
#' \item{trembles.se}{Column vector which contains the standard errors of the reported trembles. The elements correspond to the vector of estimates.}
#' \item{coefficients.se}{Column vector which contains the standard errors of the reported coefficients. The elements correspond to the vector of estimates.}
#' \item{convergence}{Row vector which reports the maximum value of the score vector of the shares as the first element, responses as the second element, trembles as the third element, and LCR coefficients as the forth element. Small values indicate convergence of the algorithm to a (local) maximum.}
#' @description Performs variants of the strategy estimation method.
#' @details The estimation function \code{stratEst()} returns maximum-likelihood estimates for the population shares and response parameters of a set of candidate strategies given some data from an economic experiment. Candidate strategies can be supplied by the user in the form of deterministic finite-state automata. The number and the complexity of strategies can be restricted by the user or selected based on information criteria. stratEst also features latent class regression to assess the influence of covariates on strategy choice.
#' @references
#' Breitmoser, Y. (2015): Cooperation, but no reciprocity: Individual strategies in the repeated prisoner's dilemma, \emph{American Economic Review}, 105, 2882-2910.
#'
#' Dal Bo, P. and G. R. Frechette (2011): The evolution of cooperation in infinitely repeated games: Experimental evidence, \emph{American Economic Review}, 101, 411-429.
#'
#' Dempster, A., N. Laird, and D. B. Rubin (1977): Maximum likelihood from incomplete data via the EM algorithm," \emph{Journal of the Royal Statistical Society Series B}, 39, 1-38.
#'
#' Eddelbuettel, D. and R. Francois (2011): Rcpp: Seamless R and C++ Integration, \emph{Journal of Statistical Software}, 40, 1-18.
#'
#' Sanderson, C. and R. Curtin (2016): Armadillo: a template-based C++ library for linear algebra. \emph{Journal of Open Source Software}, 1-26.
#' @examples
#' ## Replication of Dal Bo and Frechette (2011), Table 7 on page 424
#' ## Results for the first treatment with delta = 1/2 and R = 32 (column 1 of Table 7)
#' data <- DF2011[DF2011$treatment == 1,]
#' strats <- rbind(ALLD,ALLC,GRIM,TFT,T2,WSLS)
#' stratEst(data,strats)
#'
#' ## Latent class regression with data from Dal Bo and Frechette (2011)
#' ## For the two treatments with R = 32, introduce a dummy which is one if delta = 3/4
#' data <- DF2011[DF2011$treatment == 1 | DF2011$treatment == 4,]
#' strats <- rbind(ALLD,ALLC,GRIM,TFT,T2,WSLS)
#' covar <- as.matrix(as.numeric(data$treatment == 4 ))
#' stratEst(data,strats,covariates = covar,select="strategies")
#' @export
stratEst <- function( data, strategies, shares, covariates, cluster, response = "mixed", r.responses = "no", r.trembles = "global", select = "no", min.strategies = 1, crit = "bic", se = "yes", outer.runs = 10, outer.tol = 0, outer.max = 1000, inner.runs = 100, inner.tol = 0, inner.max = 10, lcr.runs = 1000, lcr.tol = 0, lcr.max = 1000, bs.samples = 1000, print.messages = TRUE ){
  # crude argument checks
  # check data
  if( missing(data) ) {
    stop("data not supplied")
  }
  data_frame <- as.data.frame(data)
  id <- data_frame$id
  supergame <- data_frame$supergame
  period <- data_frame$period
  group <- data_frame$group
  cooperation <- data_frame$cooperation
  if( is.null(cooperation) ) {
    cooperation <- data_frame$coop
  }
  other_cooperation <- data_frame$other_cooperation
  if( is.null(other_cooperation) ) {
    other_cooperation <- data_frame$o_coop
  }
  input <- data_frame$input
  output <- data_frame$output

  if( is.null(id) ) {
    stop("data does not contain the variable: id (id)")
  }
  if( is.null(supergame) ) {
    stop("data does not contain the variable: supergame (match)")
  }
  if( is.null(period) ) {
    stop("data does not contain the variable: period (round)")
  }
  if( ( is.null(input) == F |  is.null(output) == F ) & ( is.null(group) == F |  is.null(cooperation) == F )  ){
      stop("make sure data either contains the variables group and cooperation or the variables input and output")
  }
  if( is.null(cooperation) == F & is.null(group) & is.null(other_cooperation) ) {
    stop("if data contains the variable cooperation, it must contain the variable group or other_cooperation")
  }
  if( is.null(group) == F & is.null(cooperation) ) {
    stop("if data contains the variable group, it must contain the variable cooperation")
  }
  if( is.null(output) == F & is.null(input)  ){
    stop("if data contains the variable output, it must contain the variable input")
  }
  if( is.null(input) == F & is.null(output)  ){
    stop("if data contains the variable input, it must contain the variable output")
  }
  if( missing(strategies) ) {
    stop("strategies not supplied")
  }
  # check covariates
  if( missing(covariates) ) {
    covariates = matrix(0,1,1)
    LCR = FALSE
  } else{
    LCR = TRUE
  }
  # check cluster
  if( missing(cluster) ) {
    cluster = matrix(0,1,1)
  }
  # check outer.runs
  if ( outer.runs < 0 | outer.runs%%1 != 0 ){
    stop("Number of outer runs must be a positive integer. Default is 100.");
  }
  # check inner.runs
  if ( inner.runs < 0 | inner.runs%%1 != 0 ){
    stop("Number of inner runs must be a positive integer. Default is 100.");
  }
  # check outer.max
  if ( outer.max < 0  | outer.max%%1 != 0){
    stop("Number of outer max evaluations must be a positive integer. Default is 1000.");
  }
  # check inner.max
  if ( inner.max < 0 | inner.max%%1 != 0 ){
    stop("Number of inner max evaluations must be a positive integer. Default is 100.");
  }
  # check response
  if ( response != "mixed" & response != "pure" ){
    stop("response has to be one of the following: \"mixed\" or \"pure\". Default is \"mixed\".");
  }
  # check r.responses
  if ( r.responses != "no" & r.responses != "strategies" & r.responses != "states" & r.responses != "global"  ){
    stop("r_responses has to be one of the following: \"no\", \"strategies\", \"states\" or \"global\". Default is \"no\".");
  }
  # check r.trembles
  if ( r.trembles != "no" & r.trembles != "strategies" & r.trembles != "states" & r.trembles != "global"  ){
    stop("r_trembles has to be one of the following: \"no\", \"strategies\", \"states\" or \"global\". Default is \"no\".");
  }
  # check select
  if ( select != "no" & select != "strategies" & select != "responses"  & select != "trembles" & select != "both" & select != "all" ){
    stop("select has to be one of the following: \"no\", \"strategies\", \"responses\", \"trembles\", \"both\", or \"all\". Default is \"no\".");
  }
  # check crit
  if ( crit != "aic" & crit != "bic" & crit != "icl" ){
    stop("crit has to be one of the following: \"aic\", \"bic\", or \"icl\". Default is \"bic\".");
  }
  # check bs.samples
  if ( bs.samples < 0  | bs.samples%%1 != 0){
    stop("Number of bootstrap samples must be a positive integer. Default is 1000.");
  }

  # for future use
  newton.stepsize = 1
  penalty = 0

  # transform PD data into input output data structure
  if( is.null(cooperation) == FALSE ){
    data <- transform_pd( data_frame )
    input <- data[,4]
    output <- data[,5]
  }else{
    # prepare data
    data <- cbind(id,supergame,period,input,output)
  }

  #check strategies
  if( is.matrix(strategies) == F ){
    n_strats = strategies
    n_inputs = length( unique( input ) )
    n_outputs = sum( unique( output ) != 0 )
    strat_states = rep(c(1:n_inputs),n_strats)
    response_par = matrix(NA,n_inputs*n_strats,n_outputs)
    transition_mat = rep(1,n_inputs*n_strats) %*% t.default(c(2:n_inputs))
    strategies <- cbind(strat_states,response_par,transition_mat)
  }
  else{
    n_strats = sum( as.numeric( strategies[,1] == 1 ) )
  }
  if ( ( select == "strategies" | select == "all" ) && n_strats == 1 ){
    stop("strategies cannot be selected if there is only one strategy.");
  }
  if ( n_strats <= min.strategies && ( select == "strategies" || select == "all" ) ){
    stop("the number of strategies supplied cannot be smaller or equal to the minimum number strategies when performing strategy selection.");
  }
  #check shares
  if( missing(shares) ) {
    shares = rep( NA , n_strats )
  }


  stratEst.output <- stratEst_cpp( data, strategies, shares, covariates, LCR, cluster, response, r.responses, r.trembles, select, min.strategies, crit, se, outer.runs, outer.tol, outer.max, inner.runs, inner.tol, inner.max, lcr.runs, lcr.tol, lcr.max, bs.samples, print.messages )
  return(stratEst.output)
}

.onUnload <- function (libpath) { library.dynam.unload("stratEst", libpath)}
