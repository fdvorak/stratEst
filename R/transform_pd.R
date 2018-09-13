# This function transforms data of teh Prisoner's Dilemma into the input output structure

transform_pd <- function( data ){
  id <- data$id
  supergame <- data$supergame
  period <- data$period
  group <- data$group
  output = data$cooperation
  if( is.null(output) ) {
    output <- data$coop
  }
  input = rep(NA,length(output))
  input[(period == 1)] = 0
  unique_ids = unique(id)
  if( is.null(data$other_cooperation) ) {
    input <- transform_cpp(id, supergame, period, group, output, input, unique_ids)
  }else{
    p_output = data$other_cooperation
    input <- transform_pd_cpp(id, supergame, period, output, p_output, input, unique_ids)
  }
  # prepare data
  data <- cbind(id,supergame,period,input,output)
  input <- data[,4]                                   # this is necessary for the column label input
  data <- cbind(id,supergame,period,input,output)     # do not delete
  return(data)
}
