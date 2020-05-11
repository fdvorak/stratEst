# checks the input object cluster.id

stratEst.check.cluster.id <- function( data , sample.id ){

  cluster.id <- NULL  # remove this to make the function work

  if( cluster.id %in% colnames(data) ){
    cluster <- data[,cluster.id]
  }
  else{
    stop(paste("stratEst error: The data does not contain the variable '",cluster.id,"' "," specified as cluster id.",sep=""))
  }
  if( class(cluster) != "factor" & class(cluster) != "integer" ){
    stop(paste("stratEst error: The variable '",cluster.id,"' ","has to be of class integer or factor.",sep=""))
  }
  cluster_factor = NULL
  if( class(cluster) == "factor" ){
    cluster_factor <- cluster
  }
  cluster <- as.numeric(cluster)

  stratEst.check.cluster.id.return = list( "cluster" = cluster , "cluster.factor" = cluster_factor )

  return( stratEst.check.cluster.id.return )

}
