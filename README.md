# Introduction to FastPG

Fast phenograph-like clustering of millions of items with scores of features.

## License

This package is licensed under the MIT license, except for the Grappolo C++ library by Mahantesh Halappanavar (hala@pnnl.gov) which is included and licensed under the BSD-3 license as described in the headings of the relevant *.cpp files. Some alterations to the original Grappolo source have been made to support installation and use within an R package. The original Grappolo C++ library is available at https://hpc.pnl.gov/people/hala/grappolo.html

## Installation

This R package has dependencies outside of R, e.g. on a python package that wraps a library, `nmslib` [@Boytsov-2013]. You can either set up the dependencies and install locally, pull down a pre-built [Docker](https://www.docker.com) container, or build the Docker container yourself locally.

### Local install

This R package depends on `nmslibR`[@Mouselimis-2018], which in turn depends on the python module `nmslib`. The instructions to install `nmslib` and its prerequisites are described in the `nmslibR` package [README](https://cran.r-project.org/web/packages/nmslibR/readme/README.html). Note that this setup is somewhat involved, especially if you are not familiar with setting up python modules.

Using python from R requires R to know where to find the python you want to run. In order to avoid specifying the python version every new session (as described in the `nmslibR` document linked above), you can set a shell startup environmental variable `RETICULATE_PYTHON` to the path to python, e.g.:

```{bash}
export ENV RETICULATE_PYTHON=/usr/bin/python3
```

Once you have set up the dependencies, you can just install the `FastPG` package from
GitHub in one of the normal ways. My favorite is:

```{r}
# Requires the CRAN packages "remotes" and "BiocManager" to be installed.
BiocManager::install("sararselitsky/FastPG")
```

### Use a pre-built Docker container

A simpler way is to use a Docker container, `jefferys/fastpg:latest`, that is already set up and has the package pre-installed. You can just run it interactively, using the R command line inside it. To use the pre-built [FastPG container from DockerHub](https://hub.docker.com/repository/docker/jefferys/fastpg) you should be running on a 64 bit Intel/AMD compatible machine.

```{bash}
docker pull jefferys/fastpg:latest
docker run -it --rm -v $PWD:$PWD -w $PWD jefferys/fastpg:latest
R
   # or (singularity < 3.0)
singularity pull docker://jefferys/fastpg:latest
singularity shell -B $PWD --pwd $PWD -C fastpg-latest.simg
R
   # or (singularity 3.0+)
singularity pull docker://jefferys/fastpg:latest
singularity shell -B $PWD --pwd $PWD -C fastpg_latest.sif
```

Note that you should consider this container version like an "application" and not an environment. You may have problems installing additional packages into it. To do that, see [Extending the Docker container](#extending-the-docker-container).

### Building your own Docker container

If you want to build your own container instead of pulling a pre-build one, you can use the `Dockerfile` included in the repository in the `Docker/` directory as a guide. The `build.sh` file in the same directory automatically builds and tags the container with the name and tags used by the pre-built container at DockerHub, you should change the tags by editing the parameters in the `build.sh` file, or by manually building it and tagging it yourself.

```{bash}
git clone --single-branch https://github.com/sararselitsky/FastPG.git
cd FastPG/Docker
# Edit Docker tags in build.sh for your use
./build.sh --no-cache
```

You can just use this with Docker as above, but you can't just use a local Docker container with singularity versions less than 3.0. You have to push the container to some Docker registry before you can run it. With 3.0+ you can pull a local container directly by using `docker-daemon://` instead of `docker://` to get a local container.

### Extending the Docker container

If you want to add additional things to the container, you should build your own. You can extend the existing container either by editing the provided  `Dockerfile` or by using the existing container as the `FROM` that your own Docker container is based on.

## Clustering with FastPG

Clustering is as simple as:

```{r}
clusters <- FastPG::fastCluster( data, k, num_threads )
```

The `fastCluster()` function takes a number of additional tuning parameters, but
those have reasonable defaults.

### The `data=` parameter

The main input is the numeric data to cluster as a matrix, where rows are elements to cluster and columns are the features that make elements similar or different. Any data matrix will do, for this example we extract a 265,627 x 32 numeric data matrix from the GitHub-published [clustering benchmark mass cytometry data set](https://github.com/lmweber/benchmark-data-Levine-32-dim), sourced from the phenograph paper [@Levine-2015].

```{r}
url <- "https://github.com/lmweber/benchmark-data-Levine-32-dim/raw/master/data/Levine_32dim.fcs"
file <- "Levine_32dim.fcs"
download.file( url, file, mode="wb") # This downloads a 41.5 MB binary file
dataColumns <- c( 5:36 ) # extract only the data columns
data <-  flowCore::exprs( flowCore::read.FCS( file ))[ , dataColumns ]
```

### The `k` parameter
To cluster a data set, a local neighborhood size needs to be specified as a parameter.

```{r}
k <- 30
```

### The `num_threads` parameter

The number of cpus to use should be specified, it defaults to 1. However, it will only limit the k nearest neighbors part of the clustering (see [Internal algorithm](internal-algorithm)). The rest of the clustering will use all available cpus regardless of what this is set to. 

```{r}
num_threads <- 4
```

### Results

`fastCluster()` returns a list with two elements

* `$modularity` = The modularity of the network created from the overlapping nearest neighbor graphs.
* `$communities` = An integer vector where the nth element is the nth element from the input data. The value is the cluster that each input element has been assigned to.

Caution: -1 indicates a point that was not clustered; each can be considered their own cluster, even though they are all labeled “-1”. It is probable that you will not have any singleton clusters, but they can occur.

## Internal algorithm

FastPG utilizes the same three main steps as the phenograph algorithm [@Levine-2015, @Chen-2016], but uses fast, parallel implementations.

* The k nearest neighbors determining step is implemented using hierarchical navigable small world graphs [@Malkov-2016] via the R library nmslibR [@Mouselimis-2018].
* The nearest-neighbor distances are generated using an included parallel Jaccard metric function.
* Clustering is implemented as community detection in the graph formed from the overlapping "k best friends" for each element. This is done using a parallel Louvain algorithm as implemented by Grappolo [@Lu-2015]. Code for Grappolo has been included within this package; the standalone application with additional functionality is available for download as described in the license section above.
  
Calling `fastCluster()`  is equivalent to and is essentially implemented as the following sequence of commands:

```{r}
# Approximate k nearest neighbours
init_nms <- nmslibR::NMSlib$new( input_data= data, space= 'l2', method= 'hnsw' )
res <- init_nms$knn_Query_Batch( data, k= k, num_threads= num_threads )
ind <- res$knn_idx

# Parallel Jaccard metric
links <- FastPG::rcpp_parallel_jce(ind)
links <- FastPG::dedup_links(links)


# Parallel Louvain clustering
clusters <- FastPG::parallel_louvain( links )
```

## References

Boytsov, Leonid, and Bilegsaikhan Naidan. 2013. “Engineering Efficient and Effective Non-Metric Space Library.” In *Similarity Search and Applications - 6th International Conference, SISAP 2013, A Coruña, Spain, October 2-4, 2013, Proceedings*, edited by Nieves R. Brisaboa, Oscar Pedreira, and Pavel Zezula, 8199:280–93. Lecture Notes in Computer Science. Springer. https://doi.org/10.1007/978-3-642-41062-8_28.

Chen, Hao, Mai Chan Lau, Michael Thomas Wong, Evan W Newell, Michael Poidinger, and Jinmiao Chen. 2016. *“Cytofkit: A Bioconductor Package for an Integrated Mass Cytometry Data Analysis Pipeline.”* PLoS Comput Biol 12 (9).

Levine, Jacob H, Erin F Simonds, Sean C Bendall, Kara L Davis, El-ad D Amir, Michelle D Tadmor, Oren Litvin, et al. 2015. “Data-Driven Phenotypic Dissection of Aml Reveals Progenitor-Like Cells That Correlate with Prognosis.” *Cell* 162 (1): 184–97. https://doi.org/10.1016/j.cell.2015.05.047.

Lu, Hao, Mahantesh Halappanavar, and Ananth Kalyanaraman. 2015. “Parallel Heuristics for Scalable Community Detection.” *Parallel Computing* 47: 19–37. https://doi.org/https://doi.org/10.1016/j.parco.2015.03.003.

Malkov, Yury A., and D. A. Yashunin. 2016. “Efficient and Robust Approximate Nearest Neighbor Search Using Hierarchical Navigable Small World Graphs.” *CoRR* abs/1603.09320. http://arxiv.org/abs/1603.09320.

Mouselimis, Lampros. 2018. *NmslibR: Non Metric Space (Approximate) Library.* https://CRAN.R-project.org/package=nmslibR.
