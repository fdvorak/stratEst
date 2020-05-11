#' Creates a stratEst.strategy object.
#' @useDynLib stratEst,.registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @param inputs The names of the inputs, i.e. the levels of factor input in data.
#' @param outputs The names of the outputs, i.e. the levels of factor output in data. .
#' @param num.states An integer indicating the number states of the strategy.
#' @param responses A numeric vector which contains the response probabilities of the strategy.
#' @param trembles A numeric vector which contains the tremble probabilities of the strategy.
#' @param transitions  A vector of integers which contains the deterministic state transitions.
#' @return A stratEst.strategy object. A data.frame with the following variables:
#' \item{output.x}{The probability for output 'x'.}
#' \item{trembles}{The probability to observe a tremble.}
#' \item{input.x}{The deterministic state transition of the strategy after observing of input 'x'.}
#' @details The strategy generation function \code{stratEst.strategy()} creates a strategy as a data.frame.
#' @examples
#' ## Create the mixed nash strategy for the rock-paper-scissors game.
#' ins = c(NA,"rock","paper","scissors")
#' outs = c("rock","paper","scissors")
#' nash =  stratEst.strategy( inputs = ins , outputs = outs , num.states = 1 , responses = 1/3 )
#' @export
stratEst.strategy <- function( inputs , outputs , num.states , responses , trembles , transitions ){

  # check num.states
  if( missing( num.states ) == F ){
    if( class( num.states ) != "numeric" | length( num.states ) != 1 ){
      stop(paste("stratEst.strategy error: The input object 'num.states' must be a number.",sep=""))
    }
  }

  # check inputs
  if( missing( inputs ) ){
    stop(paste("stratEst.strategy error: The input object 'inputs' is missing.",sep=""))
  }
  else{
    if( class( inputs ) != "character"  ){
      stop(paste("stratEst.strategy error: The input object 'inputs' must be a character vector.",sep=""))
    }
  }

  # check outputs
  if( missing( outputs ) ){
    stop(paste("stratEst.strategy error: The input object 'outputs' is missing.",sep=""))
  }
  else{
    if( class( outputs ) != "character"  ){
      stop(paste("stratEst.strategy error: The input object 'outputs' must be a character vector.",sep=""))
    }
  }

  input_has_na <- as.numeric( any( is.na( inputs ) ) )
  num_inputs = length( inputs[ is.na( inputs ) == F ] )
  if( missing( num.states ) ){
    num.states = num_inputs + input_has_na
  }
  num_outputs = length( outputs )
  response_mat = matrix(NA,num.states,num_outputs)
  transition_mat = matrix(1,num_inputs,num.states)
  if( num.states > 1 ){
    transition_mat[] <- c((1+input_has_na):num.states)
  }
  transition_mat <- t(transition_mat)
  tremble_vec <- rep( NA , nrow(response_mat) )

  # check responses
  if( missing( responses ) == F ){
    if( all( is.na( responses ) ) == F ){
      if( class( responses ) != "numeric" ){
        stop(paste("stratEst.strategy error: Responses must be numeric.",sep=""))
      }
      if( any( is.na(responses) == F & ( responses < 0 | responses > 1 ) ) ){
        stop(paste("stratEst.strategy error: Responses must be values between zero and one.",sep=""))
      }
      response_mat <- t(response_mat)
      response_mat[] <- responses
      response_mat <- t(response_mat)
      response_mat_zeros <- response_mat
      response_mat_zeros[ is.na( response_mat_zeros )] = 0
      sums_responses <- apply( response_mat_zeros , 1 , sum )
      if( any( sums_responses > 1.0001 ) ){
        stop(paste("stratEst.strategy error: The column sum of responses cannot exceed one.",sep=""))
      }
    }
  }

  # check trembles
  if( missing( trembles ) == F ){
    if( all( is.na( trembles ) ) == F ){
      if( class( trembles ) != "numeric" ){
        stop(paste("stratEst.strategy error: Trembles must be numeric.",sep=""))
      }
      if( any( is.na(trembles) == F & ( trembles < 0 | trembles > 1 ) ) ){
        stop(paste("stratEst.strategy error: Trembles must be values between zero and one.",sep=""))
      }
      tremble_vec[] <- trembles
    }
  }

  # check transitions
  if( missing( transitions) == F ){
    if( class( transitions ) != "numeric" ){
      stop(paste("stratEst.strategy error: Transitions must be numeric.",sep=""))
    }
    if( any( transitions < 0 ) | any( transitions > num.states )  ){
      stop(paste("stratEst.strategy error: Transitions must be integers between one and the number of states (", as.character(num.states) , ").",sep=""))
    }
    transition_mat = t(transition_mat)
    transition_mat[] <- transitions
    transition_mat = t(transition_mat)
  }

  strategy <- as.data.frame(cbind(response_mat,tremble_vec,transition_mat))

  # column names
  output_names <- paste( "output." , outputs , sep ="" )
  input_names <- paste( "input." , inputs[ is.na( inputs ) == F ] , sep ="" )
  colnames(strategy) <- c(output_names,"tremble",input_names)

  # make object of class stratEst.strategy
  attr(strategy, "class") <- c("stratEst.strategy","data.frame")

  return(strategy)

}
