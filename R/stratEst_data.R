#' Creates a stratEst.data object.
#' @useDynLib stratEst,.registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @param data A data.frame object in the long format. The object \code{data} is mandatory.
#' @param choice A character indicating the variable in \code{data} which contains the discete choices. The default is "choice".
#' @param input A character vector indicating the names of the input generating variables in \code{data}. At least one input generating variable has to be specified.
#' @param input.lag A numeric vector indicating how many periods later each variable will affect the generation of the input. The vector must have as many elements as variables specified in the object \code{input}. The object \code{input.lag} is optional. The default is a vector of zeros.
#' @param input.sep An optional character used to seperate the values of the input generating variables. Default is no seperation character.
#' @param id A character indicating the variable in \code{data} which identifies observations of the same individual across the rows of the data frame. The default is "id".
#' @param game A character indicating the variable in \code{data} which identifies observations of the same individual in one game across the rows of the data frame. The default is "game".
#' @param period A character indicating the variable in \code{data} which identifies the period within a game. The object \code{period} is optional. The default is "period".
#' @param add An optional character vector indicating the names of variables in the global environment that should be added to the data frame.
#' @param drop An optional character vector indicating the names of variables in \code{data} that should be droped from the data frame.
#' @return A stratEst.data object, which is a data.frame in the long format with the following variables:
#' \item{data}{A data.frame object in the format required by the estimation function.}
#' @details The data generation function \code{stratEst.data()} transforms data into the format required by the estimation function \code{stratEst()}.
#' @references
#' Dal Bo, P. and G. R. Frechette (2011): The evolution of cooperation in infinitely repeated games: Experimental evidence, \emph{American Economic Review}, 101, 411-429.
#'
#' Fudenberg, D., Rand, D. G., and Dreber, A. (2012): Slow to Anger and Fast to Forgive: Cooperation in an Uncertain World, \emph{American Economic Review}, 102, 720-749.
#'
#' Wang, Z., Xu, B., and Zhou, H. (2014): Social cycling and conditional responses in the Rock-Paper-Scissors game. \emph{Scientific Reports} 4, 5830.
#'
#'#' @examples
#' ## Transform the rock-paper-scissors data of Wang, Xu, and Zhou (2014)
#' data.WXZ2014 <- stratEst.data(WXZ2014, input = c("choice"), choice = "choice", input.lag = 1)
#'
#' ## Transform the prisoner's dilemma data of Dal Bo and Frechette (2011).
#' data.DF2011 <- stratEst.data(DF2011, choice ="choice", input =c("choice","other.choice"), input.lag = 1)
#'
#' #' ## Transform the prisoner's dilemma data of Fudenberg, Rand, and Dreber (2012).
#' data.FRD2012 <- stratEst.data(data = FRD2012, choice ="choice", input =c("last.choice","last.other"))
#' @export
stratEst.data <- function( data, choice = "choice", id = "id", input = "input", input.lag = 0, input.sep = "", game = "game", period = "period", add = NULL, drop = NULL ){

  # CHECK INPUT ARGUMENTS
  # check data
  if( missing(data) ) {
    stop("stratEst.data error: Mandatory input object 'data' is missing.")
  }
  else{
    if( class(data) != "data.frame" ){
      stop("stratEst.data error: Input object 'data' must be a data frame.")
    }
  }

  # check id
  if( class(id) != "character" ){
    stop("stratEst.data error: Input object 'id' must be a character indicating the variable in 'data' which identifies the individuals.")
  }
  if( length(id) != 1 ){
    stop("stratEst.data error: Input object 'id' can only have one element.")
  }
  if( id %in% colnames(data) ) {
    id_var <- data[,id]
    norm_id <- match(id_var,sort(unique(id_var)))
    num_ids <- length(unique(norm_id))
  }else{
    stop(paste("stratEst.data error: The data does not contain the variable '",id,"' specified as variable 'id'.",sep=""))
  }

  # check choice
  if( choice %in% colnames(data) == F ) {
    stop(paste("stratEst.data error: The data does not contain the variable '",choice,"'.",sep=""))
  }
  else{
    output <- choice
    if( class(output) != "character" ){
      stop("stratEst.data error: Input object 'output' must be a character indicating the variable in 'data' which identifies the output.")
    }
    if( length(output) != 1 ){
      stop("stratEst.data error: Input object 'output' can only have one element.")
    }
    output_var <- data[,output]
  }

  # check input
  for( i in 1:length(input) ){
    if( input[i] %in% colnames(data) == F ) {
      stop(paste("stratEst.data error: The data does not contain the variable '",input[i],"'.",sep=""))
    }
  }
  if( class(input) != "character" ){
    stop("stratEst.data error: Input object 'input' must be a character vector indicating the variables in 'data' which generate the input.")
  }
  num_input_vars <- length(input)
  input_vars <- as.data.frame(matrix(NA,nrow(data),num_input_vars))
  norm_input_vars <- matrix(NA,nrow(data),num_input_vars)
  for( i in 1:num_input_vars ){
    if( input[i] %in% colnames(data) ) {
      input_vars[,i] <- as.factor(data[,input[i]])
    }else{
      stop(paste("stratEst.data error: The data does not contain the variable '",input[i],"' specified as variable to generate the input.",sep=""))
    }
    norm_input_vars[,i] <- match(input_vars[,i],sort(unique(input_vars[,i])))
  }


  # # check input.levels
  # if( class(input.levels) != "character" ){
  #   stop("stratEst.data error: Input object 'input.levels' must be a character vector.")
  # }

  # check input.lag
  if( class(input.lag) != "numeric" ){
    stop("stratEst.data error: Input object 'input.lag' must be numeric.")
  }
  if ( input.lag < 0 | input.lag%%1 != 0 ){
    stop("stratEst error: Input object 'input.lag' must be a positive integer. Default is zero.");
  }
  if( length(input.lag) != 1 ){
    stop("stratEst.data error: Input object 'input.lag' cannot have more than one element.")
  }
  input_lag <- as.integer(input.lag)


  # check game
  if( class(game) != "character" ){
    stop("stratEst.data error: Input object 'game' must be a character indicating the variable in 'data' which identifies the game.")
  }
  if( length(game) != 1 ){
    stop("stratEst.data error: Input object 'game' can only have one element.")
  }
  if( game %in% colnames(data) ) {
    game_var <- data[,game]
  }else{
    stop(paste("stratEst.data error: The data does not contain the variable '",game,"' specified as variable 'game'.",sep=""))
  }

  # check period
  if( class(period) != "character" ){
    stop("stratEst.data error: Input object 'period' must be a character indicating the variable in 'data' which identifies the period.")
  }
  if( length(period) != 1 ){
    stop("stratEst.data error: Input object 'period' can only have one element.")
  }
  if( period %in% colnames(data) ) {
    period_var <- data[,period]
  }else{
    stop(paste("stratEst.data error: The data does not contain the variable '",period,"' specified as variable 'period'.",sep=""))
  }

  # normalize and check game and period
  norm_game <- game_var
  norm_period <- period_var
  for( i in 1:max(norm_id)){
    norm_game[norm_id == i] <- match(game_var[norm_id == i],sort(unique(game_var[norm_id == i])))
    for( j in 1:max(norm_game[norm_id == i]) ){
      norm_period[norm_id == i & norm_game == j] <- match(period_var[norm_id == i & norm_game == j],sort(unique(period_var[norm_id == i & norm_game == j])))
      if( length(norm_period[norm_id == i & norm_game == j]) != max(norm_period[norm_id == i & norm_game == j]) ){
        stop("The same period cannot occur several times in the same game for the same id.")
      }
    }
  }

  # check drop and drop
  if( is.null(drop) == F ) {
    if( class(drop) != "character" ){
      stop("stratEst.data error: Input object 'drop' must be a character vector indicating the variables in 'data' that should be added to the data frame.")
    }
    num_vars_to_drop <- length(drop)
    for( i in 1:num_vars_to_drop ){
      if( drop[i] %in% colnames(data) ) {
        data <- subset(data, select = -which( colnames(data) == drop[i] ) )
      }else{
        stop(paste("stratEst.data error: The data does not contain the variable '", drop[i],"' that should be dropped from the data frame.",sep=""))
      }
    }
  }

  # check add and add
  if( is.null(add) == F ) {
    if( class(add) != "character" ){
      stop("stratEst.data error: Input object 'add' must be a character vector indicating the variables from the global environment that should be added to the data frame.")
    }
    num_vars_to_add <- length(add)
    for( i in 1:num_vars_to_add ){
      if( exists( add[i] ) ){
        old_names <- colnames(data)
        if( length( get( add[i] ) ) == nrow(data) ){
          data$add_name <- get(add[i])
          colnames(data) <- c(old_names,add[i])
        }
        else {
          stop(paste("stratEst.data error: The length variable '", add[i],"' that should be added to the data frame must correspond to the number of rows of 'data'.",sep=""))
        }
      }else{
        stop(paste("stratEst.data error: The variable '", add[i],"' that should be added to the data frame does not exist in the global environment.",sep=""))
      }
    }
  }
  ###############################################################################################
  # GENERATE INPUT VALUE MATRIX
  ###############################################################################################
  input_var <- rep(NA,nrow(data))
  input_NA_index <- rep(F,nrow(data))
  for( i in 1:num_input_vars ){
    if( i == 1 ){
      input_var <- input_vars[,i]
      input_NA_index[ is.na( input_var ) ] = T
      input_var <- as.character( input_var )
    }
    else{
      next_input_var <- input_vars[,i]
      input_NA_index[ is.na( next_input_var ) ] = T
      next_input_var <- as.character( next_input_var )
      input_var <- paste(input_var, next_input_var ,sep= input.sep )
    }
  }
  input_var[ input_NA_index ] = NA
  input_var <- as.factor(input_var)

  if( input.lag > 0 ){
    numeric_input_var <- as.numeric(input_var)

    #generate lag
    lagged_input_var <- rep(NA,length(numeric_input_var))
    lagged_input_var <- stratEst_data_cpp( norm_id , norm_game , norm_period , numeric_input_var , lagged_input_var, input_lag , num_ids )

    lagged_input_var <- as.factor(lagged_input_var)
    levels(lagged_input_var) <- levels(input_var)
    input_var <- lagged_input_var
  }

  # if( length(input.levels) != length(unique(input_var))){
  #   warning("stratEst.data warning: The number of supplied input.levels is not equal the number of levels of the generated input.")
  # }
  # input_var <- input.levels[input_var]

  stratEst.data.return <- data
  stratEst.data.return$id <- as.integer(id_var)
  stratEst.data.return$game <- as.integer(game_var)
  stratEst.data.return$period <- as.integer(period_var)
  stratEst.data.return$input <- as.factor(input_var)
  stratEst.data.return$choice <- as.factor(output_var)

  # make object of class stratEst.data
  attr(stratEst.data.return, "class") <- c("stratEst.data","data.frame")

  # return result
  return(stratEst.data.return)

}

.onUnload <- function (libpath) { library.dynam.unload("stratEst", libpath)}
