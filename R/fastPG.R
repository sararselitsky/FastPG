# FastPG.R - The fast phenograph implementation.

#' Load the data matrix from an FCS file
#'
#' @param data Can be a one element character vector giving the FCS file to
#'   load or URL to download from, or a connection to
#'   read the file from, or a flowCore::flowFrame object.
#' @param signalColumns The columns to keep.
#'   Either a vector of names or integers (positions). By default will keep
#'   all columns, which is probably not what you want. Any column not found
#'   causes an error, unless force = TRUE is set, in which case missing
#'   columns will be ignored with a warning.
#' @param drop Drop rows containing NA values.
#'   By default is `TRUE`; An `NA` in any column will drop that row from the
#'   matrix. Can be set to `FALSE` to keep all rows, or a vector of column
#'   names or positions to check. Applied before any columns are dropped by
#'   `signalColumns=`. Any column not found causes an error, unless force=
#'   TRUE is set, in which case missing columns will be ignored with a
#'   warning.
#' @param sample Select only this many rows of data.
#'   Default is all data rows. If `sample` exceeds the number of rows of data,
#'   will ignore with a warning. Applied after any rows are dropped due to
#'   `NA` values.
#' @param force Try and process anyway, despite errors.
#' @return A matrix of measurements for event rows (cells) x feature columns
#'   (signal channels).
#' @export
#'
#' @examples
fcsMatrix <- function( data, signalColumns= NULL, drop= TRUE, sample= FALSE, force= FALSE ) {
    
    if ( ! is(data, "flowFrame" )) {
        data <- read.FCS( data )
    }
    
    if (! file.exists(file) ) {
        stop( paste0( "Can't find FCS file: \"", file, "\"" ), call.= FALSE )
    }
    if (! isFCSfile(file)) {
        message <- paste0( "Appears not to be an FCS file: \"", file, "\"\n" )
        if (force) {
            message <- paste0( message, "\tAttempting to process anyway.")
            warning( message )
        } else {
            stop( message )
        }
    }
    
    fcsObj <- read.FCS( file )
}

### Quick and dirty get it running code, fix later
### TODO: Split to a couple of separate files and document

helpText <- "
USAGE:
./fastPG.Rscript [-g] <FCS_FILE>

-g          Switch to use parallel-louvain clustering with grappolo

<FCS-FILE>  Filename (path) to FCS format data file. May only contain
data columns.

NOTE:
Requires the fastPG.Rscript and the parallel_jc2.cpp file to be in the
same directory.
"


#' Parse command line
#'
#' @param args 
#'
#' @return
#' @export
#'
#' @examples
parseCli <- function( args= commandArgs( trailingOnly= TRUE ) ) {
    opts <- list(
        fcsFile= "",
        isGrappolo= FALSE,
        k= 30,
        nThread= 5,
        outFile= "",
        relFile= ""
    )
    if ( length(args) == 0 ) {
        stop( paste0( "Missing argument: FCS file name.\n", helpText), call.= FALSE )
    }
    else if (length(args) == 2 ) {
        if ( args[[1]] != "-g" ) {
            stop(  paste0( "Unknown option: '", args[[1]], "'.\n", helpText), call.= FALSE )
        }
        opts$isGrappolo <- TRUE
        opts$fcsFile <- args[[2]]
    }
    else {
        opts$fcsFile <- args[[1]]
    }
    
    if (! file.exists(opts$fcsFile) ) {
        stop( paste0( "Can't find FCS file: \"", opts$fcsFile, "\"" ), call.= FALSE )
    }
    
    if( opts$outFile == "") {
        opts$outFile= paste0(
            file.path(
                dirname( opts$fcsFile ),
                sub( "[.].*$","",  basename( opts$fcsFile ))
            ),
            ".clusterMembership"
        )
    }
    if( opts$relFile == "") {
        opts$relFile= paste0(
            file.path(
                dirname( opts$fcsFile ),
                sub( "[.].*$","",  basename( opts$fcsFile ))
            ),
            ".relations"
        )
    }
    
    opts   
}

#' Read in an FCS file
#'
#' @param fcsFile 
#'
#' @return
#' @export
#'
#' @examples
loadFCS <- function( fcsFile ) {
    # TODO column select - to be added as feature
    
    fcsDF <- as.data.frame( flowCore::exprs( flowCore::read.FCS( fcsFile )))
    

    as.matrix(fcsDF)
    
}

# 
#' Hierarchical navigable small world relations map
#'
#' @param dat 
#' @param k 
#' @param nThread 
#'
#' @return
#' @export
#'
#' @examples
hnsw <- function( dat, k= 30, nThread= 5 ) {
    
    
    ind <- NULL
    init_nms <- nmslibR::NMSlib$new( input_data= dat, space= 'l2', method= 'hnsw' )
    res <- init_nms$knn_Query_Batch( dat, k= k, num_threads= nThread )
    ind <- res$knn_idx
    
    links <- rcpp_parallel_jce( ind )
    
    relations <- as.data.frame( links )
    colnames( relations ) <- c( "from", "to", "weight" )
    badRowSelect <- relations$from == 0 | relations$to == 0 | relations$weight == 0
    relations[ ! badRowSelect, ]
}

#' Louvain clustering
#'
#' @param relations 
#'
#' @return
#' @export
#'
#' @examples
clusterLouvain <- function( relations ) {
    g <- igraph::graph.data.frame( relations, directed= FALSE )
    community <- igraph::cluster_louvain( g )
    
    # return cluster assignment
    community$membership
}

#' Output data to file for running Grapallo
#'
#' @param relations 
#' @param file 
#' @param k 
#' @param nThread 
#'
#' @return
#' @export
#'
#' @examples
outputForGrapallo <-function(
    relations, file, k,
    nThread= data.table::getDTthreads( FALSE )
) {
    nodes= nrow(relations)
    headerLine <- paste0( "# Nodes: ", nodes, " Edges: ", k * nodes )
    writeLines( headerLine, file )
    data.table::fwrite( relations, file, sep="\t", nThread= nThread, col.names=FALSE, append= TRUE )
}

#' Run Grapallo parallel louvain clustering
#'
#' @param file 
#' @param exec 
#' @param f 
#' @param c 
#' @param o 
#' @param v 
#'
#' @return
#' @export
#'
#' @examples
runGrapallo <- function(
    file,
    exec= "./bin/driverForGraphClustering",
    f= 8, c= 1, o= TRUE, v= TRUE
) {
    o <- ifelse(o, "-o", NULL)
    v <- ifelse(v, "-v", NULL)
    argsVec= c("-f", f, -c, c, o, v, file)
    system2(command = exec, argsVec)
}
