
# FastPG

Fast phenograph-like clustering of millions of items with scores of features.

***This is still under development and is not necessarily returning correct results***

## Installation

This R package has dependencies outside of R, e.g. on a python package that wraps a library, `nmslibs`. You can either set up the dependencies and install locally, pull down a prebuilt docker if you have a modern system, or build the docker container yourself locally.

### Local install

[TODO: Set up dependencies]

If you have the dependencies, you can just install this from
GitHub in one of the normal ways. My favorite is:

```
# Requires the CRAN packages "remotes" and "BiocManager" to be installed.
BiocManager::install("sararselitsky/FastPG@rgrappolo")
```

### Use a pre-built Docker container

A simpler way is to use a Docker container, `jefferys/fastpg:latest`, that is already set up and has the package pre-installed. You can just run it interactively, using the R command line inside it. To use the pre-built FastPG container on DockerHub you should run on a 64 bit intel/amd compatible machine with cpus that support the SSE3, SSE4.1, SSE4.2, and AVX instructions. This should be true for most systems. See "Building your own container" if if you have an older system, or if you want to take advantage of the AVX2 instructions.

```
docker pull jefferys/fastpg:latest
docker run -it --rm -v $PWD:$PWD -w $PWD jefferys/fastpg:latest
R
   # or
singularity pull docker://jefferys/fastpg:latest
singularity shell -B $PWD --pwd $PWD -C fastpg-latest.simg
R
```

Note that you should consider this container version like an "application" and not an environment. You will likely have problems installing additional packages into it. To do that, see "Extending the container".

### Building your own container

If you want to build your own container instead of pulling a pre-build one, you can use the `Dockerfile` included in the repository in the `Docker/` directory. The `build.sh` file in the same directory automatically tags the container you build with the names used at DockerHub, you can change these by editing the parameters in the build.sh file. This should automatically build a version that supports the cpu family of the machine built on.

```
git clone --single-branch --branch rgrappolo https://github.com/sararselitsky/FastPG.git
cd FastPG/Docker
# Edit docker build tags in build.sh for your use
build.sh --no-cache
```

You can just use this with docker as above, but you can't just use a local docker container with singularity versions less than 3.0. You have to push the container to some docker registry before you can run it. With 3.0+ you can pull a local container directly by using `docker-daemon://` instead of `docker://` to get a local container.

### Extending the container

If you want to add additional things to the container, you must build your own. You can extend the existing container either by editing the Dockerfile or by using the existing container as the FROM to base your own docker container on.

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
init_nms <- nmslibR::NMSlib$new( input_data= data, space= 'l2', method= 'hnsw' )
res <- init_nms$knn_Query_Batch( data, k= k, num_threads= num_threads )
ind <- res$knn_idx

links <- FastPG::rcpp_parallel_jce(ind)
links <- links[ links[, 1] != 0 ]
links <- matrix( links, ncol= 3 )

num_nodes <- length( union( links[, 1], links[, 2] ))

clusters <- FastPG::parallel_louvain( links, num_nodes )
```

## References
