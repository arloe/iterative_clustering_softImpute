#' Training function
#' @param X the matrix with NAs
#' @param iteration maximum number of iterations
#' @param thresh convergence threshold

cluster_softImpute <- function( X, iteration = 500, thresh = 1e-5){
  
  pb <- txtProgressBar(min = 0, max = iteration, style = 3)
  for( i in 1:iteration ){
    
    # index of NA
    idx_NA = is.na(X)
    
    # initialize NA to 0
    if( i == 1 ){
      X_pre = X
      X_pre[is.na(X_pre)] = 0
    }
    # impute using cluster_softInpute
    X_post = model_fit(X = X_pre, X_org = X )
    
    # calculate Frobenius norm between pre-estimators and post-estimators
    e        = X_post - X_pre
    coverage = mean(x = e[idx_NA]^2)
    
    # if coverage is less than thresh, stop process
    # else update matrix
    if( coverage < thresh ){
      print("The coverage reached a value less than threshold")
      break
    }else{
      X_pre = X_post
    }
    setTxtProgressBar(pb = pb, value = i, title = "The ProgressBar")
  }
  return( X_post )
}

#' The function that calculate clustering and softImpute once
#' @importFrom softImpute softImpute complete
#' @importFrom wordspace dist.matrix
#' @importFrom dendextend find_k
#'
#' @param X
#' @param X_org the original 

model_fit <- function(X, X_org){

  ##### hierarchical clustering #####
  X_dist = wordspace::dist.matrix(M = X, as.dist = TRUE)
  fit_hc = hclust(d = X_dist, method = "ave")

  optimal_k = find_k(hc = fit_hc, dist_mat = X_dist )

  pred_cluster = cutree(tree = fit_hc, k = optimal_k)
  
  # split X matrix by predicted cluster
  split_X = lapply(X = 1:optimal_k, function(x) subset(x = X, subset = pred_cluster == x))
  split_X_org = lapply(X = 1:optimal_k, function(x) subset(x = X_org, subset = pred_cluster == x))
  
  # check the number of rows of each matrix
  # If the number of rows in the matrix is 1, the softimpute algorithm does not work
  K = 1:optimal_k
  K_one = which( x = sapply(X = split_X, FUN = nrow) == 1 )
  K_normal = K[ !K %in% K_one ]
  
  ##### softImpute #####
  fit_softimpute = lapply(X = split_X[K_normal], FUN = softImpute::softImpute )
  pred_X = sapply(X = K_normal, FUN = function(x) softImpute::complete(x = split_X_org[[x]], object = fit_softimpute[[x]] ) )
  pred_X = do.call(what = "rbind", args = pred_X)
  
  # concatenate with one row matrix 
  if( length( x = split_X[K_one] ) != 0 ) pred_X = rbind(pred_X, split_X[K_one])
  
  # sort
  pred_X = pred_X[ rownames(X), ]
  
  return( pred_X )
}


#' The function that find the most optimal number of cluster using average of silhouette width
#' 
#' @importFrom dendextend nleaves
#' @importFrom cluster silhouette
#' 
#' @param hierarchical clustering model
#' @param dist_mat matrix of distance

find_k <- function( hc, dist_mat ){
  if( class(x = hc) != "hclust" ) stop(" 'hc' must be 'hclust' class ")
  if( class(x = dist_mat) != "dist" ) stop(" 'dist_mat' must be 'dist' class ")
  # k range
  k_range = 2:min(10, (dendextend::nleaves(x = hc) - 1))
  
  # calculate average silhouette width
  out = data.frame( k = integer(), avg_sil = numeric() )
  for( i in 1:length(k_range) ){
    sil     = cluster::silhouette(x = cutree(tree = hc, k = k_range[i]), dist = dist_mat )
    avg_sil = summary(sil)$avg.width
    out[i, ] = c( k_range[i], avg_sil )
  }
  optimal_k = out$k[which.max(x = out$avg_sil)]
  
  return( optimal_k )
}


