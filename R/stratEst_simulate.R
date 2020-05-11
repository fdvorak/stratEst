#' Simulation function for strategy estimation.
#' @useDynLib stratEst,.registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @param strategies Mandatory input object. Can be either a positive integer or a matrix. If an integer is used, the estimation function will generate the respective number of memory-one strategies with as many states as there are unique input values in \code{data}. A matrix can be used to supply a set of customized strategies. In the matrix, each row corresponds to one state of a strategy, starting with the start state of an automaton. The first column enumerates the states of each strategy in ascending order. A value of one in the first column indicates the begin of a new strategy with its start state. The columns after the first column contain the collection of multinomial response vectors. The number of columns for the multinomial response vectors must correspond to the number of unique non-zero outputs in data. Without a reference output - which is labeled with a zero in the output column of data - the columns specify the complete multinomial response distribution for each unique value in the output column. In this case, the response probabilities in each row must sum to one. With a reference output, the response probability for the response labeled with zero is omitted and the response probabilities in each row must sum to a value smaller or equal to one. The remaining columns of the strategies matrix define the deterministic state transitions. The number of columns must equal the number of unique non-zero inputs in the data. The numbers in the first column indicate the next state of the automaton if the input is one. The numbers in the second column indicate the next state if the input is two and so on.
#' @param shares A matrix of strategy shares. The number of columns must correspond to the number of samples in the data. The number of rows must correspond to the number of strategies defined in the strategies matrix. Elements which are NA are estimated from the data. If the object is not supplied, one share is estimated for every strategy defined in the strategies matrix for every sample.
#' @param coefficients Column vector which contains the latent class regression coefficients. The elements correspond to the vector of estimates.
#' @param covariates A character vector indicating the names of the variables in data that are the covariates of the latent class regression model. Rows with the same id must have the values of covariates. Missing value are not allowed. Whenever a character vector is supplied for the input object 'covariates', a latent class regression model is estimated.
#' @param num.ids Integer which specifies the number of individuals which are assigned to the stategies. Default is 100.
#' @param num.games Integer which specifies the number of games each individual plays. Default is 5 games.
#' @param num.periods A column vector with as many elments as there are games. The elements of the vector are integer vaules which specify the number of periods of each game. Default is 5 periods in each game.
#' @param fixed.assignment A logical indicating whether the assignment of individuals to strategies is fixed. If \code{FALSE} individuals are repeatedly assigned to strategies for each each game. If \code{TRUE}, individuals are only assigned once, and use the assigned strategy in each game. Default is \code{FALSE}.
#' @param input.na A logical indicating whether the input in the first period of a game is \code{NA}. If \code{FALSE}, and input is randomly selected for the first period. Default is \code{FALSE}.
#' @description Simulates data which can be used to test the strategy estimation function \code{stratEst()}.
#' @examples
#' ## Fictitious data from a helping game
#' ## Participant 62 plays reciprocal strategy.
#' ## Participant 87 plays alternating strategy.
#' @export
stratEst.simulate <- function( strategies , shares , coefficients , covariates , num.ids = 100 , num.games = 5 , num.periods  , fixed.assignment = FALSE , input.na = FALSE ){

  ##################################################################################################
  # CHECK INPUTS
  ##################################################################################################

  # check strategies
  if( missing(strategies) ) {
    stop("stratEst error: Mandatory input object 'strategies' is missing. Specify an integer or create a data.frame object.")
  }else{
    stratEst.simulate.check.strategies.return <- stratEst.simulate.check.strategies( strategies )
    strategies_matrix_list <- stratEst.simulate.check.strategies.return$strategies.matrix.list
    num_samples_strategies <- length(strategies_matrix_list)
    trembles_list <- stratEst.simulate.check.strategies.return$trembles.list
    num_strats <- stratEst.simulate.check.strategies.return$num_strats
    unique_inputs <- stratEst.simulate.check.strategies.return$unique.inputs
    unique_outputs <- stratEst.simulate.check.strategies.return$unique.outputs
    num_inputs <- stratEst.simulate.check.strategies.return$num.unique.inputs
    num_outputs <- stratEst.simulate.check.strategies.return$num.unique.outputs
    sid <- stratEst.simulate.check.strategies.return$sid
    names_strategies_lcr <- stratEst.simulate.check.strategies.return$names.strategies
  }

  if( missing( shares ) & ( missing( coefficients ) | missing( covariates ) ) ){
    stop("stratEst error: Either the input object 'shares' or the the inputs objects 'coefficients' and 'covariates' must be supplied.")
  }

  # check covariates
  if( missing(covariates) ){
    LCR = FALSE
  }
  else{
    stratEst.simulate.check.covariates.return <- stratEst.simulate.check.covariates(covariates)
    covariate_mat <- stratEst.simulate.check.covariates.return$covariate.mat
    num.ids <- stratEst.simulate.check.covariates.return$num.ids
    num.covariates <- stratEst.simulate.check.covariates.return$num.covariates
    names.covariates <- stratEst.simulate.check.covariates.return$names.covariates
    LCR = TRUE
  }

  #check shares
  if( missing(shares) == F ) {
    stratEst.simulate.check.shares.return <- stratEst.simulate.check.shares( shares , LCR , num_strats )
    shares = stratEst.simulate.check.shares.return$shares
    num_samples_shares = stratEst.simulate.check.shares.return$num_samples
    names_samples = stratEst.simulate.check.shares.return$names_samples
  }else{
    num_samples_shares = 1
    names_samples = "1"
  }

  #check coefficients
  if( missing(coefficients) ) {
    if( LCR ){
      stop("stratEst error: Input object 'coefficients' must be supplied.")
    }
  }else{
    coefficient_mat <- stratEst.simulate.check.coefficients( coefficients , covariate_mat , num_strats , names_strategies_lcr )
    if( LCR == F ){
      warning("stratEst warning: No covariates specified. The input object 'coefficients' is ignored.");
    }
  }

  # check num.ids
  num.ids <- as.integer(num.ids)
  if( class(num.ids) != "integer" ){
    stop("stratEst error: Input object 'num.ids' has to be of class 'integer'.")
  }else{
    if( num.ids < 0 ){
      stop("stratEst error: Input object 'num.ids' cannot be negative.")
    }
  }

  # check num.games
  if( class(num.games) != "numeric" ){
    stop("stratEst error: Input object 'num.games' has to be of class 'numeric'.")
  }else{
    num.games <- round(num.games)
    if( num.games < 0 ){
      stop("stratEst error: Input object 'num.games' cannot be negative.")
    }
  }

  # check num.periods
  if( missing( num.periods ) ){
    num.periods <- rep(5,num.games)
  }else{
    if( class(num.periods) != "numeric" ){
      stop("stratEst error: Input object 'num.periods' has to be of class 'numeric'.")
    }else{
      if( length(num.periods) != num.games ){
        stop("stratEst error: Input object 'num.periods' has to be a vector of length 'num.games'.")
      }
      num.periods <- round(num.periods)
      if( any(num.periods < 0 ) ){
        stop("stratEst error: Elements of the input object 'num.periods' cannot be negative.")
      }
    }
  }

  # check fixed.assignment
  if( class(fixed.assignment) != "logical" ){
    stop("stratEst error: Input object 'fixed.assignment' has to be of class 'logical'.")
  }

  # balance num sample information
  if( num_samples_strategies == num_samples_shares ){
    num_samples = num_samples_strategies
  }
  else if( num_samples_strategies == 1 & num_samples_shares > 1 ){
    num_samples = num_samples_shares
    for( sam in 2:num_samples ){
      strategies_matrix_list[[sam]] <- strategies_matrix_list[[1]]
      trembles_list[[sam]] <- trembles_list[[1]]
    }
  }
  else if( num_samples_strategies > 1 & num_samples_shares == 1 ){
    num_samples = num_samples_strategies
    shares <- matrix( shares , length(shares) , num_samples_strategies )
  }
  else{
    stop("stratEst error: The number of samples indicated in input object 'strategies' and 'shares' does not match.")
  }


  ##################################################################################################
  # GENERATE DATA
  ##################################################################################################

  data_matrix <- NULL
  max_id_sample <- 0

  for( sam in 1:num_samples ){

    strategies_matrix <- strategies_matrix_list[[sam]]
    trembles <- trembles_list[[sam]]
    trembles[is.na(trembles)] <- 0
    trembles <- trembles
    response_mat <- matrix(strategies_matrix[,2:(num_outputs+1)],nrow(strategies_matrix),num_outputs)
    tremble_mat <- matrix(rep(trembles,ncol(response_mat)),nrow(response_mat),ncol(response_mat))
    response_mat <- response_mat*(1-tremble_mat) + (1-response_mat)*tremble_mat/(ncol(response_mat)-1)
    transition_mat <- matrix(strategies_matrix[,(num_outputs+2):(num_inputs+1+num_outputs)],nrow(strategies_matrix),num_inputs)

    if( LCR ){
      priors_individuals = exp( covariate_mat %*% coefficient_mat)/( apply(exp(covariate_mat %*% coefficient_mat),1,sum) )
    }else{
      priors_individuals = t( matrix( shares[,sam] , nrow(shares) , num.ids ) )
    }

    # assign ids to strategies (strategies assignmen mat)
    strategy_assignment_mat <- matrix(NA,num.ids,num.games)
    if( fixed.assignment ){
      rand_vec = stats::runif(num.ids)
      for( i in 1:num.ids){
        assigned = F
        for( s in 1:num_strats ){
          if( assigned == F ){
            if( rand_vec[i] <= sum( priors_individuals[i,(1:s)]) ){
              strategy_assignment_mat[i,] <- rep( s , num.games )
              assigned = T
            }
          }
        }
      }
    }else{
      rand_mat = matrix(stats::runif(num.ids*num.games),num.ids,num.games)
      for( i in 1:num.ids){
        for( j in 1:num.games ){
          assigned = F
          for( s in 1:num_strats ){
            if( assigned == F ){
              if( rand_mat[i,j] <= sum( priors_individuals[i,(1:s)]) ){
                strategy_assignment_mat[i,j] <- s
                assigned = T
              }
            }
          }
        }
      }
    }

    # generate variables
    if( fixed.assignment ){
      id = max_id_sample + rep(c(1:num.ids),each=sum(num.periods))
    }else{
      id <- NULL
      for( i in 1:num.games ){
        id = max_id_sample + c(id,rep(c(((i-1)*num.ids + 1):((i-1)*num.ids + num.ids)),each=num.periods[i]))
      }
    }
    game = rep( NA,num.ids*sum(num.periods) )
    period = rep(NA,num.ids*sum(num.periods))
    input = rep(0,num.ids*sum(num.periods))
    output = rep(NA,num.ids*sum(num.periods))
    sample = rep(names_samples[sam],num.ids*sum(num.periods))
    strategy.id = rep(NA,num.ids*sum(num.periods))
    if( LCR ){
      covariate_vars = matrix(NA,num.ids*sum(num.periods),num.covariates)
    }

    # generate output
    counter <- 1
    for( i in 1:num.ids ){
      for( g in 1:num.games ){
        sid_sbj = strategy_assignment_mat[i,g]
        indices_strat <- sid == sid_sbj
        response_mat_sbj <- matrix(response_mat[ indices_strat , ],sum(indices_strat),ncol(response_mat))
        transition_mat_sbj <- matrix(transition_mat[ indices_strat , ],sum(indices_strat),ncol(transition_mat))
        state <- 1
        for( p in 1:num.periods[g] ){
          strategy.id[counter] <- sid_sbj
          game[counter] <- g
          period[counter] <- p
          if( p > 1 | input.na == FALSE ){
            input[counter] <- t(stats::rmultinom( 1 , 1 , rep(1,num_inputs)/num_inputs)) %*% c(1:num_inputs)
            state <- transition_mat_sbj[state,input[counter]]
          }
          output[counter] = t(stats::rmultinom( 1 , 1 , response_mat_sbj[state,] )) %*% c(1:num_outputs)
          if( LCR ){
              covariate_vars[counter,] <- covariate_mat[i,]
          }
          counter <- counter + 1
        }
      }
    }

    unique_output_values <- gsub("output.", "", unique_outputs)
    unique_input_values <- gsub("input.", "", unique_inputs)

    output <- unique_output_values[output]
    unique_input_values <- c(NA,unique_input_values)
    input <- unique_input_values[input+1]

    max_id_sample <- max(id)

    data_matrix_sample <- cbind(id,game,period,input,output,sample,strategy.id)

    if( LCR ){
      data_matrix_sample <- cbind(id,game,period,input,output,sample,strategy.id,covariate_vars)
    }

    data_matrix <- rbind(data_matrix,data_matrix_sample)

  }

  data <- as.data.frame( data_matrix )

  if( LCR ){
    colnames(data)[8:(8+num.covariates-1)] <- colnames(covariate_mat)
  }

  # make object of class stratEst.data
  attr(data, "class") <- c("stratEst.data","data.frame")

  return(data)
}
