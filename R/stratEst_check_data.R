# checks the input object data

stratEst.check.data <- function( data ){

  if( "data.frame" %in% class(data) == F ){
    stop("stratEst error: The input object 'data' must be an object of class 'stratEst.data'.")
  }

  # check mandatory data frame variables
  if( is.null(data$id) ) {
    stop("stratEst error: Data does not contain the variable 'id'.")
  }
  if( is.null(data$game) ) {
    stop("stratEst error: Data does not contain the variable 'game'.")
  }
  if( is.null(data$period) ) {
    stop("stratEst error: Data does not contain the variable 'period'.")
  }
  if( is.null(data$input) ) {
    stop("stratEst error: Data does not contain the variable 'input'.")
  }
  if( is.null(data$output) ) {
    stop("stratEst error: Data does not contain the variable 'output'.")
  }

  id <- data$id
  game <- data$game
  period <- data$period
  input <- data$input
  output <- data$output

  # id
  id_is_factor = FALSE
  if( class(id) == "factor" ){
    id_factor <- id
    id_is_factor = TRUE
  }
  id <- as.numeric(id)
  if( any( is.na( id ) ) ){
    stop("stratEst error: The variable 'id' in data cannot contain NA values.");
  }

  # input
  if( class(input) == "factor" ){
    input_factor <- input
  }
  else{
    input_factor <- as.factor( input )
  }
  input <- match(input_factor,sort(unique(input_factor)))
  input[is.na(input)] <- 0
  levels_input <- levels(input_factor)

  # output
  if( class(output) == "factor" ){
    output_factor <- output
  }
  else{
    output_factor <- as.factor( output )
  }
  output <- match(output_factor,sort(unique(output_factor)))
  levels_output <- levels( output_factor )
  if( any( is.na( output ) ) ){
    stop("stratEst error: The variable 'output' in data cannot contain NA values.");
  }

  stratEst.check.data.return <- list( "data" = data , "id" = id , "game" = game , "period" = period , "input" = input , "output" = output, "input.factor" = input_factor , "output.factor" = output_factor , "levels.input" = levels_input , "levels.output" = levels_output )

  return(stratEst.check.data.return)
}
