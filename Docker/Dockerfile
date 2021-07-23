FROM jefferys/bioconductor:3.13

RUN apt-get update \
 && apt-get install -y --no-install-recommends \
      build-essential \
      libtool \
      m4 \
      automake \
      autoconf \
      git

RUN R -e 'BiocManager::install( \
    c( "flowCore", "data.table", "RANN", "Rcpp", "RcppParallel", \
       "igraph", "RcppHNSW", "mclust", "MLmetrics", "remotes" ), \
    update= FALSE, \
    ask= FALSE, \
    version= "3.13" \
)'

RUN R -e 'BiocManager::install( \
    c( "sararselitsky/FastPG" ), \
    update= FALSE, \
    ask= FALSE, \
    version= "3.13" \
)'

# At end as args breaks all following caching. Timestamp changes every time...

ARG imageCreated
ARG toolVersion

ARG toolSource="https://github.com/sararselitsky/FastPG"
ARG tool="fastpg"
ARG brief="${tool} - Container for running the R package FastPG."
ARG toolLicense="MIT"

# From https://github.com/opencontainers/image-spec/blob/master/annotations.md
LABEL org.opencontainers.image.created="${imageCreated}"
LABEL org.opencontainers.image.authors="Stuart R. Jefferys <srj@unc.edu>"
LABEL org.opencontainers.image.url="https://hub.docker.com/r/jefferys/${tool}"
LABEL org.opencontainers.image.source="https://github.com/jefferys/Dockers/ToolDockers/${tool}"
LABEL org.opencontainers.image.version="${toolVersion}"
LABEL org.opencontainers.image.vendor="UNC - Lineberger"
LABEL org.opencontainers.image.licenses="${toolLicense}"
LABEL org.opencontainers.image.title="${brief}"
LABEL org.opencontainers.image.description="The fastpg R package for clustering large sample sets (millions) with intermediate features (order of 100). Includes source from Grappolo, BSD-3 licensed from Pacific Northwest National Laboratory"

# https://github.com/BioContainers/specs/blob/master/container-specs.md
LABEL base_image="jefferys/bioconductor:3.13_latest"
LABEL version="${toolVersion}"
LABEL software.version="${toolVersion}"
LABEL software="${tool}"
LABEL about.summary="${brief}"
LABEL about.home="${toolSource}"
LABEL about.license="${toolLicense}"
LABEL about.tags="debian, bioconductor, fastpg"
LABEL maintainer="Stuart R. Jefferys <srj@unc.edu>"

CMD ["bash"]
