#' Fast approximate clustering
#'
#' [TODO: find way to set default threads to max]
#' [TODO: how to estimat k? Sensitivity if wrong?]
#' [TODO: Description of algorithm; references]
#'  
#' @param data A numeric matrix of rows to cluster based on similarity of.
#' column features.
#' @param k Expected number of clusters.
#' @param num_threads Number of threads to use. Defaults to 1.
#'
#' @return Returns a list with two elements:
#' * `modularity` - [TODO]
#' * `communities` - An integer vector where the nth element is the nth
#'   row in the input matrix. Its value is the cluster that row has been assigned to.
#'
#' @examples
#' # [TODO]
#' 
#' @export
fastCluster <- function( data, k= 30, num_threads= 1 ) {
  init_nms <- nmslibR::NMSlib$new( input_data= data, space= 'l2', method= 'hnsw' )
  res <- init_nms$knn_Query_Batch( data, k= k, num_threads= num_threads )
  ind <- res$knn_idx
  
  links <- FastPG::rcpp_parallel_jce(ind)
  links <- links[ links[, 1] != 0 ]   
  links <- matrix( links, ncol= 3 )
  
  num_nodes <- length( union( links[, 1], links[, 2] ))
  
  FastPG::parallel_louvain( links, num_nodes ) 
}