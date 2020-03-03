
# FastPG

Fast phenograph-like clustering of millions of items with scores of features.

***This is still under development and is not necessarily returning correct results***

## Installation

This R package has dependencies outside of R, e.g. on a python package that wraps a library, `nmslibs`. If you have set up the dependencies, you can just install this from
GitHub in one of the normal ways. My favorite is:

```
# Requires the CRAN packages "remotes" and "BiocManager" to be installed.
BiocManager::install("sararselitsky/FastPG@rgrappolo")
```

A simpler way is to use a Docker container that is already set up and has the package pre-installed. You can just run it interactively, using the R command line inside it. To use the pre-built FastPG container on DockerHub:

```
docker pull jefferys/FastPG:latest
docker run -it --rm -v $PWD:$PWD -w $PWD jefferys/fastpg:latest
R
   # or
singularity pull docker://jefferys/FastPG:latest
singularity shell -B $PWD -C fastpg-latest.simg
R
```

You can also build a docker container from the `Dockerfile` included in the repository. The `build.sh` file in the `Docker/` directory automatically tags the container you build with the same names as used at DockerHub.

```
git clone --single-branch --branch rgrappolo https://github.com/sararselitsky/FastPG.git
cd FastPG/Docker
build.sh
```

You can just use this with docker as above, but you can't just use a local docker container with singularity versions less than 3.0. You have to push the container to some docker registry before you can run it. With 3.0+ you can pull a local container directly by using `docker-daemon://` instead of `docker://` to get a local container


## Running FastPG

Load the data to cluster into a matrix, rows will be matched by index to clusters.

[TODO: Include a function to pull a sample data file for use]

```
file <- "someFile.FCS"
dataColumns <- c( 5:35 ) # extract only the data columns, whatever they are
dat <-  flowCore::exprs( flowCore::read.FCS( file ))[ , dataColumns ]
```

Select the k for clustering, should be something like the number of clusters expected.
Sensitivity to this parameter is being explored. [TODO, advice]

```
k <- 30
```

Set the number of threads/cores to use: [TODO, defaults]

```
num_threads <- 4
```

Cluster the data:

```
clusters <- FastPG::fastCluster( data, k, num_threads )
```

This returns a list with two elements

* `$modularity` = [TODO...]
* `$communities` = An integer vector where the nth element is the nth element in the input data. Its value is the cluster that input element has been assigned to.

Calling `fastCluster()`  is equivalent to the following sequence of commands:

```
init_nms <- nmslibR::NMSlib$new( input_data= dat, space= 'l2', method= 'hnsw' )
res <- init_nms$knn_Query_Batch( dat, k= k, num_threads= num_threads )
ind <- res$knn_idx

links <- FastPG::rcpp_parallel_jce(ind)
links <- links[ links[, 1] != 0 ]
links <- matrix( links, ncol= 3 )

num_nodes <- length( union( links[, 1], links[, 2] ))

clusters <- FastPG::parallel_louvain( links, num_nodes )
```


