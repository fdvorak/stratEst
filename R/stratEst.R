#' Estimation function for strategy estimation
#' @useDynLib stratEst,.registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @param data Mandatory data frame object which contains the data for the estimation. Data has to be in the long format with variables in columns. Each row in \code{data} represents one observation. Three columns are mandatory: A column named \code{id} which identifies the observations of the same individual across the rows of the data frame. A column named \code{input} which indicates the type of information observed by the individual before giving a response. A column named \code{output} which contains the behavioral response of the individual after observing the input. If an individual plays the same game for more than one period with the same partner, \code{data} must contain a variable \code{period} which identifies the period within the game. If an individual plays the same game more than once with different partners, \code{data} must contain a variable \code{game} (or \code{supergame}) which identifies data from different games. For data from prisoner's dilemma experiments, two more data formats are possible. Instead of using the variables \code{input} and \code{output}, the data frame may also contain the variables \code{cooperation} and \code{other_cooperation}, or alternatively, the variables  \code{cooperation} and \code{group}. The variable \code{cooperation} should be a dummy which indicates if the participant cooperated in the current period. The variable \code{other_cooperation} should be a dummy which indicates if the other player cooperated in the current period. The variable \code{group} should be an identifier variable with a unique value for each unique match of two individuals.
#' @param strategies Mandatory input object. Can be either a positive integer or a matrix. If an integer is used, the estimation function will generate the respective number of memory-one strategies with as many states as there are unique input values in \code{data}. A matrix can be used to supply a set of customized strategies. In the matrix, each row corresponds to one state of a strategy, starting with the start state of an automaton. The first column enumerates the states of each strategy in ascending order. A value of one in the first column indicates the begin of a new strategy with its start state. The columns after the first column contain the collection of multinomial response vectors. The number of columns for the multinomial response vectors must correspond to the number of unique non-zero outputs in data. Without a reference output - which is labeled with a zero in the output column of data - the columns specify the complete multinomial response distribution for each unique value in the output column. In this case, the response probabilities in each row must sum to one. With a reference output, the response probability for the response labeled with zero is omitted and the response probabilities in each row must sum to a value smaller or equal to one. The remaining columns of the strategies matrix define the deterministic state transitions. The number of columns must equal the number of unique non-zero inputs in the data. The numbers in the first column indicate the next state of the automaton if the input is one. The numbers in the second column indicate the next state if the input is two and so on.
#' @param shares A column vector of strategy shares. The number of elements must correspond to the number of strategies defined in the strategies matrix. Elements which are NA are estimated from the data. If the object is not supplied, a share is estimated for every strategy defined in the strategies matrix.
#' @param sample.id A character indicating the variable in data which contains an id for samples. Individual observations must be nested in samples. The same must be true for clusters if specified. If more than one sample exists, shares are estimated for each sample. All other parameters are estimated for the data of all samples. If the object is not supplied, it is assumed that the data contains only one sample.
#' @param cluster.id A character indicating the variable in data which contains an id for clusters. Individual observations must be nested in clusters. for block-bootstrapped standard errors. Note that estimates will nevertheless be biased due to the non-linearity of the model.
#' @param covariates A character vector indicating the covariates in data for latent class regression. Rows with the same id must have the values of covariates. Missing value are not allowed. Whenever a character vector is supplied for the input object 'covariates', a latent class regression model is estimated.
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
#' @param stepsize Positive number which sets the stepsize of the Fisher scoring algorithm used to update the coefficients of the latent class regression model. Default is one. Values smaller than one slow down the convergence of the algorithm.
#' @param penalty A logical value indicating if the Firth penalty is used to estimate the coefficients of the latent class regression model. Default is FALSE. Independent of the logical value specified, the penalty is always used in the bootstrap of the standard errors of latent class regression coefficients.
#' @param print.messages Logical, if \code{TRUE} messages are printed which illustrate the status of the estimation process.
#' @note The strategy estimation method was introduced by (Dal Bo & Frechette 2011) to estimate the relative frequency of a fixed set of pure strategies in the indefinitely repeated prisoner's dilemma. Breitmoser (2015) extended the method to the estimation of behavior strategies. The \pkg{stratEst} package uses the EM algorithm (Dempster, Laird & Rubin 1977) and the Newton-Raphson method to obtain maximum-likelihood estimates for the population shares and response parameters of a set of candidate strategies. The package builds on other software contributions of the R community. To increase speed the estimation procedures, the package uses integration of C++ and R achieved by the Rcpp package (Eddelbuettel & Francois 2011) and the open source linear algebra library for the C++ language RppArmadillo (Sanderson & Curtin 2016).
#' @return The function returns a list with the following elements.
#' \item{shares}{Matrix which contains the estimates of population shares for the strategies. The order of rows corresponds to the order of strategies defined in the input object 'strategies'. Columns correspond to the samples defined in 'samples'. Can be used as input object of the estimation function.}
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
#' \item{state.obs}{A column vector with the number of weighted observations for each strategy state corresponding to the rows of \code{strategies}.}
#' \item{assignments}{Matrix which contains the posterior probability assignments of individuals to strategies. The rows of the matrix correspond to the ID sorted in ascending order beginning with the individual with the lowest ID. The columns correspond to the strategies, starting with the first strategy defined in the strategy matrix in column one.}
#' \item{priors}{Matrix which contains the individual prior probabilities of individuals as predicted by the covariate vectors of the individuals. The rows correspond to the ID sorted in ascending order beginning with the individual with the lowest ID. The columns correspond to the strategies, starting with the first strategy defined in the strategy matrix.}
#' \item{shares.se}{Matrix which contains the standard errors of the estimated shares. The elements correspond to the matrix of estimates.}
#' \item{responses.se}{Column vector which contains the standard errors of the reported responses. The elements correspond to the vector of estimates.}
#' \item{trembles.se}{Column vector which contains the standard errors of the reported trembles. The elements correspond to the vector of estimates.}
#' \item{coefficients.se}{Column vector which contains the standard errors of the reported coefficients. The elements correspond to the vector of estimates.}
#' \item{convergence}{Row vector which reports the maximum value of the score vector of shares, responses, trembles, and latent class regression coefficients. Small values indicate convergence of the algorithm to a (local) maximum of the negative log likelihood. Because of the sum to one constraint the score of the shares will not be close to zero if some of the shares are fixed.}
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
#' ## Fictitious data from a helping game
#' ## Participant 62 plays reciprocal strategy.
#' ## Participant 87 plays alternating strategy.
#' id <- c(62,62,62,62,87,87,87,87)
#' game <- c(4,4,4,4,4,4,4,4)
#' period <- c(1,2,3,4,1,2,3,4)
#' input <- c(0,1,2,3,0,1,3,2)
#' output <- c(2,2,1,2,2,1,2,1)
#' data <- as.data.frame(cbind(id,game,period,input,output))
#' strategies <- matrix(c(1,2,3,1,2,0.5,0,1,0.1,NA,0.5,1,0,0.9,NA,2,2,2,2,1,
#' 3,3,3,2,1,2,2,2,2,1,3,3,3,2,1),5,7)
#' model <- stratEst(data,strategies)
#'
#' ## Replication of Dal Bo and Frechette (2011), Table 7 on page 424
#' ## Results for the first treatment with delta = 1/2 and R = 32 (column 1 of Table 7)
#' data <- DF2011[DF2011$treatment == 1,]
#' strategies <- rbind(ALLD,ALLC,GRIM,TFT,WSLS,T2)
#' stratEst(data,strategies)
#'
#' ## Latent class regression with data from Dal Bo and Frechette (2011)
#' ## For the two treatments with R = 32, introduce a dummy which is one if delta = 3/4
#' dummy <- as.numeric(DF2011$treatment > 3 )
#' data <- as.data.frame(cbind(DF2011,dummy))
#' strats <- rbind(ALLD,TFT)
#' stratEst(data,strats,covariates = c("dummy"),lcr.runs = 500)
#' @export
stratEst <- function( data, strategies, shares , coefficients , sample.id , cluster.id , covariates, response = "mixed", r.responses = "no", r.trembles = "global", select = "no", min.strategies = 1, crit = "bic", se = "yes", outer.runs = 10, outer.tol = 0, outer.max = 1000, inner.runs = 100, inner.tol = 0, inner.max = 10, lcr.runs = 1000, lcr.tol = 0, lcr.max = 1000, bs.samples = 1000, stepsize = 1 , penalty = F , print.messages = TRUE ){
  # crude argument checks
  # check data
  if( missing(data) ) {
    stop("Mandatory input object data is missing.")
  }
  data_frame <- as.data.frame(data)
  id <- data_frame$id
  supergame <- data_frame$game
  if( is.null(supergame) ) {
    supergame <- data_frame$supergame
  }
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

  # generate variable sample if missing
  if( missing(sample.id) ){
    sample <- rep(1,length(id))
  }
  else{
    if( sample.id %in% colnames(data_frame) ){
      sample <- data_frame[,sample.id]
    }
    else{
      message_sample <- paste("The data does not contain the variable '",sample.id,"' "," specified as sample id.",sep="")
      stop(message_sample)
    }
  }

  data_frame$sample <- sample
  num_samples <- length(unique(sample))


  # generate variable cluster if missing
  if( missing(cluster.id) ){
    cluster <- matrix(0,1,1)
  }
  else{
    if( cluster.id %in% colnames(data_frame) ){
      cluster <- data_frame[,cluster.id]
    }
    else{
      message_cluster <- paste("The data does not contain the variable '",cluster.id,"' "," specified as cluster id.",sep="")
      stop(message_cluster)
    }
  }


  # generate variable covariates if missing
  if( missing(covariates) ){
    covariate_mat <- matrix(0,1,1)
    LCR = FALSE
  }
  else{
    covariate_mat <- NULL
    for( i in 1:length(covariates) ){
      if( covariates[i] %in% colnames(data_frame) ){
        covariate_mat <- cbind(covariate_mat,data_frame[,covariates[i]])
      }
      else{
        message_covariate <- paste("The data does not contain the variable '",covariates[i],"' "," specified as covariate.",sep="")
        stop(message_covariate)
      }
      LCR = TRUE
    }
  }

  # check mandatory data frame variables
  if( is.null(id) ) {
    stop("Data does not contain the variable: id")
  }
  if( ( is.null(input) == F |  is.null(output) == F ) & ( is.null(group) == F |  is.null(cooperation) == F )  ){
    stop("Make sure data contains the variables input and output. For data from the prisoner's dilemma the variables cooperation and group or cooperation and other_cooperation can also be used istead.")
  }
  if( is.null(cooperation) == F & is.null(group) & is.null(other_cooperation) ) {
    stop("If data contains the variable cooperation, it must contain the variable group or other_cooperation.")
  }
  if( is.null(group) == F & is.null(cooperation) ) {
    stop("If data contains the variable group, it must contain the variable cooperation.")
  }
  if( is.null(output) == F & is.null(input)  ){
    stop("If data contains the variable output, it must contain the variable input.")
  }
  if( is.null(input) == F & is.null(output)  ){
    stop("If data contains the variable input, it must contain the variable output.")
  }
  if( missing(strategies) ) {
    stop("Mandatory input object strategies is missing. Use either an integer or a strategy matrix. ")
  }

  # generate variable supergame if missing
  if( is.null(supergame) ) {
    data_frame$supergame <- rep(1,length(id))
    supergame <- data_frame$supergame
  }

  # generate variable period if missing
  if( is.null(period) ){
    data_frame$period <- rep(1,length(id))
    period <- data_frame$period
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
    stop("The input object response has to be one of the following: \"mixed\" or \"pure\". Default is \"mixed\".");
  }

  # check r.responses
  if ( r.responses != "no" & r.responses != "strategies" & r.responses != "states" & r.responses != "global"  ){
    stop("The input object r.responses has to be one of the following: \"no\", \"strategies\", \"states\" or \"global\". Default is \"no\".");
  }

  # check r.trembles
  if ( r.trembles != "no" & r.trembles != "strategies" & r.trembles != "states" & r.trembles != "global"  ){
    stop("The input object r.trembles has to be one of the following: \"no\", \"strategies\", \"states\" or \"global\". Default is \"no\".");
  }

  # check select
  if ( select != "no" & select != "strategies" & select != "responses"  & select != "trembles" & select != "both" & select != "all" ){
    stop("The input object select has to be one of the following: \"no\", \"strategies\", \"responses\", \"trembles\", \"both\", or \"all\". Default is \"no\".");
  }

  # check crit
  if ( crit != "aic" & crit != "bic" & crit != "icl" ){
    stop("The input object crit has to be one of the following: \"aic\", \"bic\", or \"icl\". Default is \"bic\".");
  }

  # check bs.samples
  if ( bs.samples < 0  | bs.samples%%1 != 0){
    stop("The number of bootstrap samples must be a positive integer. Default is 1000.");
  }

  # check stepsize
  if ( stepsize < 0 ){
    stop("The newton.stepsize must be a positive number. Default is 1.");
  }

  # check penalty
  if (  is.logical(penalty) == F){
    stop("The input argument 'penalty' must be boolean. Default is FALSE.");
  }

  # sort data
  # transform PD data into input output data structure
  if( is.null(cooperation) == FALSE ){
    # sort data frame
    data_frame <- data_frame[order(data_frame$period), ]
    data_frame <- data_frame[order(data_frame$supergame), ]
    data_frame <- data_frame[order(data_frame$id), ]
    # transform data
    data <- transform_pd( data_frame )
    input <- data[,4]
    output <- data[,5]
  }else{
    # prepare data
    data <- cbind(id,supergame,period,input,output,sample)
    # sort data
    data <- data[order(data[,3]), ]
    data <- data[order(data[,2]), ]
    data <- data[order(data[,1]), ]
  }

  #check strategies
  if( length(strategies) == 1 & sum( unique(strategies%%1) ) == 0 ){
    integer_strategies = T
    n_strats = strategies
    n_inputs = length( unique( input ) )
    n_outputs = sum( unique( output ) != 0 )
    strat_states = rep(c(1:n_inputs),n_strats)
    response_par = matrix(NA,n_inputs*n_strats,n_outputs)
    transition_mat = rep(1,n_inputs*n_strats) %*% t.default(c(2:n_inputs))
    strategies <- cbind(strat_states,response_par,transition_mat)
    tremble <- rep( NA , nrow(strategies) )
  }
  else if( is.data.frame(strategies) ){
    integer_strategies = F
    state <- strategies$state
    if( missing(state) ){
      stop("The input object 'strategies' does not contain the variable 'state'.")
    }
    # check and generate responses
    unique_outputs <- unique( output[output != 0] )
    sorted_unique_outputs <- sort( unique_outputs )
    num_unique_outputs <- length( sorted_unique_outputs )
    response_mat <- matrix(NA,nrow(strategies),num_unique_outputs)
    for( out in 1:num_unique_outputs ){
      r_string <- as.character(paste( "r",as.character(sorted_unique_outputs[out]), sep=""))
      if( r_string %in% colnames(strategies) ){
        response_mat[,out] <- strategies[,r_string]
      }
      else{
        message <- paste("stratEst error: There is an output with value ", sorted_unique_outputs[out] , " in the data but there is no column named '", r_string , "' in strategies.",sep="")
        stop(message)
      }
    }
    # check and generate transitions
    unique_inputs <- unique( input[input != 0] )
    order_inputs <- order( unique_inputs )
    sorted_unique_inputs <- sort( unique_inputs )
    num_unique_inputs <- length( sorted_unique_inputs )
    transition_mat <- matrix(NA,nrow(strategies),num_unique_inputs)
    for( ins in 1:num_unique_inputs ){
      t_string <- paste( "t",as.character(sorted_unique_inputs[ins]), sep="")
      if( t_string %in% colnames(strategies) ){
        transition_mat[,ins] <- strategies[,t_string]
      }
      else{
        message <- paste("stratEst error: There is an input with value ", sorted_unique_inputs[ins] , " in the data but there is no column named '", t_string , "' in strategies.",sep="")
        stop(message)
      }
      if( sum( transition_mat[,ins]%%1==0 ) < nrow(strategies) ){
        message <- paste("stratEst error: The transition columns in 'strategies' must be integers. Check the column named '", t_string , "'.",sep="")
        stop(message)
      }
      # check for tremble column
      if( "tremble" %in% colnames(strategies) ){
        tremble <- strategies[,"tremble"]
      }
      else{
        tremble <- rep( NA , nrow(strategies) )
      }
    }
    # generate strategy id if missing
    sid <- strategies$sid
    if( is.null(sid) ){
      sid <- rep(NA,nrow(strategies))
      n_strats <- 0
      for( i in 1:nrow(strategies) ){
        if( state[i] == 1 ){
          n_strats <- n_strats + 1
        }
        sid[i] <- n_strats
      }
    }
    else{
      n_strats <- length(unique(sid))
    }
    strategies_matrix = cbind(state,response_mat,transition_mat)
  }
  else{
    stop("The input object 'strategies' must either be an integer or a data frame.");
  }

  if ( ( select == "strategies" | select == "all" ) && n_strats == 1 ){
    stop("Strategies cannot be selected if there is only one strategy.");
  }


  #check shares
  if( missing(shares) ) {
    shares = matrix( NA , n_strats , num_samples )
  }

  #check coefficients
  if( missing(coefficients) ) {
   coefficients = matrix(0,1,1)
  }
  else{
    if( LCR = F){
      stop("There are no covariates specified for the coefficients. Use the input object 'covariates' to specify the names of the columns which contain the covariates in data.");
    }
   }

  # make coefficients input object and fixable
  cpp.output <- stratEst_cpp( data, strategies_matrix, shares , coefficients, covariate_mat, LCR, cluster, response, r.responses, r.trembles, select, min.strategies, crit, se, outer.runs, outer.tol, outer.max, inner.runs, inner.tol, inner.max, lcr.runs, lcr.tol, lcr.max, bs.samples, print.messages, integer_strategies, stepsize , penalty )
  # make data.frame out of strategies and skip responses, trembles
  stratEst.return <- list("shares" = cpp.output$shares, "strategies" = cpp.output$strategies, "responses" = cpp.output$responses, "trembles" = cpp.output$trembles,  "coefficients" = cpp.output$coefficients, "response.mat" = cpp.output$response.mat, "tremble.mat" = cpp.output$tremble.mat, "coefficient.mat" =  cpp.output$coefficient.mat, "loglike" = cpp.output$fit[1,1], "crit.val" = cpp.output$fit[1,2], "eval" = cpp.output$solver[1,1], "tol.val" = cpp.output$solver[1,2], "entropy" = cpp.output$fit[1,3], "state.obs" = cpp.output$state.obs, "assignments" = cpp.output$assignments, "priors" = cpp.output$priors, "shares.se" = cpp.output$shares.se, "responses.se" = cpp.output$responses.se, "trembles.se" = cpp.output$trembles.se, "coefficients.se" = cpp.output$coefficients.se, "shares.covar" = cpp.output$stats.list$shares.covar, "shares.score" =  cpp.output$stats.list$shares.score, "shares.fisher" = cpp.output$stats.list$shares.fisher, "responses.covar" = cpp.output$stats.list$responses.covar, "responses.score" = cpp.output$stats.list$responses.score, "responses.fisher" = cpp.output$stats.list$responses.fisher, "trembles.covar" = cpp.output$stats.list$trembles.covar, "trembles.score" = cpp.output$stats.list$trembles.score, "trembles.fisher" = cpp.output$stats.list$trembles.fisher, "coefficients.covar" = cpp.output$stats.list$coefficients.covar, "coefficients.score" = cpp.output$stats.list$coefficients.score, "coefficients.fisher" = cpp.output$stats.list$coefficients.fisher, "convergence" = cpp.output$convergence );
  # delete empty list entries
  stratEst.return <- stratEst.return[lapply(stratEst.return,length)>0]
  # reform strategies
  cbind()

  # return result
  return(stratEst.return)

}

.onUnload <- function (libpath) { library.dynam.unload("stratEst", libpath)}
