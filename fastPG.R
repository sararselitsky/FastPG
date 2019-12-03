#!/usr/bin/env Rscript

library(flowCore)
library(doMC)
registerDoMC(10)
library(data.table)
library(RANN)
library(Rcpp)
library(igraph)
library(nmslibR)

args = commandArgs( trailingOnly= TRUE )
if ( length(args) == 0 ) {
    stop( "Missing argument: FCS file name.", call.= FALSE )
}
fcsFile <- args[[1]]
if (! file.exists(fcsFile) ) {
    stop( paste0( "Can't find FCS file: \"", fcsFile, "\"" ), call.= FALSE )
}

l<-as.data.frame(exprs(read.FCS(fcsFile)))
l<-l[which(!is.na(l$label)),]

sourceCpp(file = "parallel_jc2.cpp")

# choose columns
dat<-as.matrix(l[,5:36])

# HNSW
ind<-NULL
init_nms <- NMSlib$new(input_data=dat, space='l2', method='hnsw')
# can set k
res <- init_nms$knn_Query_Batch(dat, k=30, num_threads=5)
ind<-res$knn_idx

links <- rcpp_parallel_jce(ind)
relations <- as.data.frame(links)
colnames(relations)<- c("from","to","weight")
relations <- relations[!(relations$from==0 | relations$to==0 | relations$weight==0),]

g <- graph.data.frame(relations, directed=FALSE)
community <- cluster_louvain(g)
# cluster assignment
nms<-community$membership

# output
out<-data.frame("cell"=rownames(dat),"cluster"=nms)
