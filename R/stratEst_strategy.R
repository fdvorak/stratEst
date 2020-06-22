#' Creates a stratEst.strategy object.
#' @useDynLib stratEst,.registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @param choices A character vector indicating the levels of the variable 'output' in data.
#' @param inputs A character vector indicating the levels of the variable 'input' in data.
#' @param num.states An integer indicating the number states of the strategy.
#' @param prob.choices A numeric vector which contains the response probabilities of the strategy.
#' @param tr.inputs  A vector of integers which contains the deterministic state transitions.
#' @param trembles A numeric vector which contains the tremble probabilities of the strategy.
#' @return A stratEst.strategy object. A data.frame with the following variables:
#' \item{prob.x}{The probability of choice 'x'.}
#' \item{tremble}{The probability to observe a tremble.}
#' \item{tr(x)}{The deterministic state transition of the strategy after input 'x'.}
#' @details The strategy generation function \code{stratEst.strategy()} creates a strategy as a data.frame.
#' @examples
#' ## Nash equilibrium strategy of rock-paper-scissors
#' ins = c(NA,"rock","paper","scissors")
#' rps = c("rock","paper","scissors")
#' mixed = stratEst.strategy(choices = rps)
#' nash = stratEst.strategy(choices = rps, prob.choices = rep(1/3,3))
#' rock = stratEst.strategy(choices = rps, prob.choices = c(1,0,0))
#' @export
stratEst.strategy <- function( choices, inputs = NULL, prob.choices = NULL, tr.inputs = NULL, trembles = NULL, num.states = NULL ){

  # check num.states
  if( is.null(num.states) == F ){
    if( class( num.states ) != "numeric" | length( num.states ) != 1 ){
      stop(paste("stratEst.strategy error: The input object 'num.states' must be an integer.",sep=""))
    }
  }

  # check inputs
  if( is.null( inputs ) == F ){
    if( class( inputs ) != "character"  ){
      stop(paste("stratEst.strategy error: The input object 'inputs' must be a character vector.",sep=""))
    }
    input_has_na <- as.numeric( any( is.na( inputs ) ) )
    num_inputs = length( inputs[ is.na( inputs ) == F ] )
    if( is.null( num.states ) ){
      num.states = num_inputs + input_has_na
    }
    transition_mat = matrix(1,num_inputs,num.states)
    if( num.states > 1 ){
      transition_mat[] <- c((1+input_has_na):num.states)
    }
    transition_mat <- t(transition_mat)

    # check tr.inputs
    if( is.null( tr.inputs) == F ){
      if( class( tr.inputs ) != "numeric" ){
        stop(paste("stratEst.strategy error: tr.inputs must be numeric.",sep=""))
      }
      if( any( tr.inputs < 0 ) | any( tr.inputs > num.states )  ){
        stop(paste("stratEst.strategy error: tr.inputs must be integers between one and the number of states (", as.character(num.states) , ").",sep=""))
      }
      transition_mat = t(transition_mat)
      if( length(tr.inputs) > length(c(transition_mat))  ){
        stop(paste("stratEst.strategy error: There are more elements in tr.inputs than needed (", as.character(num.states*num_inputs) , ").",sep=""))
      }
      if( length(c(transition_mat)) %% length(tr.inputs) != 0  ){
        stop(paste("stratEst.strategy error: There number of elements in tr.inputs is not a multiple of the elements required (", as.character(num.states*num_inputs) , ").",sep=""))
      }
      transition_mat[] <- tr.inputs
      transition_mat = t(transition_mat)
    }
  }else{
    if( is.null(num.states) ){
      num.states = 1
    }
  }

  # check choices
  if( missing( choices ) ){
    stop(paste("stratEst.strategy error: The input object 'choices' is missing.",sep=""))
  }
  else{
    if( class( choices ) != "character"  ){
      stop(paste("stratEst.strategy error: The input object 'choices' must be a character vector.",sep=""))
    }
  }

  num_outputs = length( choices )
  response_mat = matrix(NA,num.states,num_outputs)
  tremble_vec <- rep( NA , nrow(response_mat) )

  # check prob.choices
  if( is.null( prob.choices ) == F ){
    if( all( is.na( prob.choices ) ) == F ){
      if( class( prob.choices ) != "numeric" ){
        stop(paste("stratEst.strategy error: prob.choices must be numeric.",sep=""))
      }
      if( any( is.na(prob.choices) == F & ( prob.choices < 0 | prob.choices > 1 ) ) ){
        stop(paste("stratEst.strategy error: prob.choices must be values between zero and one.",sep=""))
      }
      response_mat <- t(response_mat)
      if( length(prob.choices) > length(c(response_mat))  ){
        stop(paste("stratEst.strategy error: There are more elements in prob.choices than needed (", as.character(num.states*num_outputs) , ").",sep=""))
      }
      if( length(c(response_mat)) %% length(prob.choices) != 0  ){
        stop(paste("stratEst.strategy error: There number of elements in prob.choices is not a multiple of the elements required (", as.character(num.states*num_outputs) , ").",sep=""))
      }
      response_mat[] <- prob.choices
      response_mat <- t(response_mat)
      response_mat_zeros <- response_mat
      response_mat_zeros[ is.na( response_mat_zeros )] = 0
      sums_prob.choices <- apply( response_mat_zeros , 1 , sum )
      if( any( sums_prob.choices > 1.0001 ) ){
        stop(paste("stratEst.strategy error: The column sum of prob.choices cannot exceed one.",sep=""))
      }
    }
  }

  # check trembles
  if( is.null( trembles ) == F ){
    if( all( is.na( trembles ) ) == F ){
      if( class( trembles ) != "numeric" ){
        stop(paste("stratEst.strategy error: trembles must be numeric.",sep=""))
      }
      if( any( is.na(trembles) == F & ( trembles < 0 | trembles > 1 ) ) ){
        stop(paste("stratEst.strategy error: trembles must be values between zero and one.",sep=""))
      }
      if( length(trembles) > length(tremble_vec)  ){
        stop(paste("stratEst.strategy error: There are more elements in prob.choices than needed (", as.character(num.states) , ").",sep=""))
      }
      if( length(tremble_vec) %% length(trembles) != 0  ){
        stop(paste("stratEst.strategy error: There number of elements in trembles is not a multiple of the elements required (", as.character(num.states) , ").",sep=""))
      }
      tremble_vec[] <- trembles
    }
  }

  include_tremble = F
  if( is.null( prob.choices ) == F ){
    include_tremble = any( prob.choices[ is.na(prob.choices) == F ] == 0 )
  }

  if( is.null(inputs) == F & include_tremble == T ){
    strategy <- as.data.frame(cbind(response_mat,tremble_vec,transition_mat))
  }
  else if( is.null(inputs) == T & include_tremble == T ){
    strategy <- as.data.frame(cbind(response_mat,tremble_vec))
  }
  else if( is.null(inputs) == F & include_tremble == F ){
    strategy <- as.data.frame(cbind(response_mat,transition_mat))
  }
  else{
    strategy <- as.data.frame(response_mat)
  }

  # column names
  output_names <- paste( "prob." , choices , sep ="" )
  if( is.null(inputs) == F ){
    input_names <- paste( "tr(" , inputs[ is.na( inputs ) == F ] , ")" , sep ="" )
  }
  if( is.null(inputs) == F & include_tremble == T ){
    colnames(strategy) <- c(output_names,"tremble",input_names)
  }
  else if( is.null(inputs) == T & include_tremble == T ){
    colnames(strategy) <- c(output_names,"tremble")
  }
  else if( is.null(inputs) == F & include_tremble == F ){
    colnames(strategy) <- c(output_names,input_names)
  }
  else{
    colnames(strategy) <- c(output_names)
  }

  # make object of class stratEst.strategy
  attr(strategy, "class") <- c("stratEst.strategy","data.frame")

  return(strategy)

}

.onUnload <- function (libpath) { library.dynam.unload("stratEst", libpath)}
