
#' A function that generate simulation heterogeneous matrix 
#'
#' @param n the number of rows of each cluster
#' @param p the number of columns
#' @param group the number of clusters

make_simulation_matrix <- function( n, p, group){
  mat = list()
  
  for( i in 1:group ){
    u = rnorm(n = n, mean = i, sd = 1)
    v = matrix( data = rnorm(n = p, mean = i, sd = ) )
    lambda = i
    e = rnorm(n = n*p, mean = 0, sd = 1) # random noise
    
    mat[[i]] = matrix(u) %*% lambda %*% matrix(v, nrow = 1) + e
  }
  out = do.call(what = "rbind", args = mat)
  
  # row name & column name
  rownames(out) = paste0("r", 1:nrow(out))
  colnames(out) = paste0("c", 1:ncol(out))
  
  return( out )
}

#' A function that create missing value in matrix
#' 
#' @importFrom missForest prodNA
#' 
#' @param x the matrix
#' @param prob the ratio to make NA 

matrix_NA <- function( x, prob ){
  if( 0 > prob | prob > 1 ) stop( "'prob' must be between 0 and 1" )  
  out = missForest::prodNA(x = x, noNA = prob )
  return( out )
}