#' Fast approximate clustering
#'
#' Clusters a numeric matrix assigning rows to clusters based on similarity of
#' features (columns). Follows the general outline of the PhenoGraph algorithm to do
#' this, but uses approximate parallel implementations when creating k
#' nearest-neighbor graphs (HNSW, hierarchical navigable small world) and when
#' clustering them (Grappolo, a parallel Louvain community detection
#' algorithm). This greatly increases clustering speed and hence the size of
#' data matrices that can be practically clustered, at the expense of
#' guaranteed optimal clustering.
#' 
#' Takes a variety of tuning parameters in addition to the main parameters of
#' `data`, `k`, and `num_threads`. These are all for optimizing the Grappolo
#' component of the clustering process and have reasonable defaults. More
#' information on the tuning parameters is given in the documentation for
#' \code{\link{parallel_louvain}}.
#' 
#' @section References: See the vignette
#'   \href{../doc/intro.html}{\code{vignette("intro.html", package =
#'   "FastPG")}}.
#' 
#' @param data A numeric matrix of rows to cluster based on the similarity of
#'   its column features.
#' @param k (3) How many nearest neighbors to choose for each point when
#'   creating the community graph.
#' @param num_threads (1) Number of threads to use.
#' 
#' @param coloring (1) Integer tuning flag between 0 and 3 that controls the
#'   type of distance-1 graph coloring. 0 = no coloring; 1 (default) =
#'   distance-1 graph coloring; 2= 1 with rebalancing; 3= Incomplete coloring
#'   with `numColors` colors.
#' @param minGraphSize (1,000) Integer tuning parameter. Change processing
#'   when graph size has reduced enough.
#' @param numColors (16) Integer tuning parameter between 1 and 1024. Limits
#'   graph coloring. Only used if `coloring=3` is set.
#' @param C_thresh (1e-6) Numeric tuning parameter > 0 and < 1. Change
#'   processing when modularity gain is too small. Not used if coloring is
#'   disabled.
#' @param threshold (1e-9) Numeric tuning parameter > 0 and < 1. Change
#'   processing when modularity gain is too small.
#' @param syncType (0) Integer tuning flag between 0 and 4 that controls
#'   synchronization, which retreat from parallelization back to more serial
#'   processing. 0 - (Default) No synchronization; 1 - Full synchronization;
#'   2 - Neighborhood synchronization, a hybrid between 0 and 1; 3 - Early
#'   termination; 4 - Full sync with early termination, a hybrid of 1 and 3.
#' @param basicOpt (1) Integer tuning flag of 0 or 1, controls internal data
#'   representation mode. 0 - A map/hash based structure; 1 - (Default) Use a
#'   vector/indexed structure.
#'
#' @return Returns a list with two elements:
#' * `modularity` - A measure of the connectedness of a clustered network.
#' When comparing different clusterings of the same network, the one with the
#' higher modularity is "better".
#' * `communities` - An integer vector where the i'th element is the i'th
#'   row in the input matrix. Its value is the cluster that row has been assigned to.
#'
#' @examples
#' \donttest{
#' # Note, this example requires the `flowCore` package, but only as part of
#' # obtaining a suitable large matrix to cluster.
#' 
#' # This downloads a 41.5 MB source binary file
#' url <- "https://github.com/lmweber/benchmark-data-Levine-32-dim/raw/master/data/Levine_32dim.fcs"
#' file <- "Levine_32dim.fcs"
#' download.file( url, file, mode="wb")
#'
#' # Convert to a suitably formatted data matrix to cluster,
#' # extracting only the data columns from the downloaded file.
#' dataColumns <- c( 5:36 )
#' data <-  flowCore::exprs( flowCore::read.FCS( file ))[ , dataColumns ]
#' 
#' # Cluster the matrix with fastCluster
#' clusters <- FastPG::fastCluster( data, k, num_threads )
#' }
#' 
#' @seealso
#'   * \code{\link{parallel_louvain}} (Grappolo)
#'   * \href{../../nmslibR/doc/the_nmslibR_package.html}{\code{vignette("the_nmslibR_package.html", package = "nmslibR")}}, (HNSW)
#'   * \href{../doc/intro.html}{\code{vignette("intro.html", package= "FastPG")}} (FastPG vignette)
#'   
#' @export
fastCluster <- function(
  data, k= 30, num_threads= 1,
  coloring= 1, minGraphSize= 1000, numColors= 16, C_thresh= 1e-6,
  threshold= 1e-9, syncType= 0, basicOpt= 1
) {
  init_nms <- nmslibR::NMSlib$new( input_data= data, space= 'l2', method= 'hnsw' )
  res <- init_nms$knn_Query_Batch( data, k= k, num_threads= num_threads )
  ind <- res$knn_idx
  
  links <- FastPG::rcpp_parallel_jce(ind)
  links <- dedup_links(links)
  
  FastPG::parallel_louvain(
    links, coloring= coloring, minGraphSize= minGraphSize, numColors= numColors,
    C_thresh= C_thresh, threshold= threshold, syncType= syncType,
    basicOpt= basicOpt
  )
}
