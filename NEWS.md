# FastPG 0.0.7

* Replaced the Python library based HNSW KNN (nmslibR) with an Rcpp based one (RcppHNSW).
* Removed Python as a dependency.
* Upgraded the FastPG Docker to use R 4.1.0 and Bioconductor 3.13.
* Fixed examples using flowCore::read.FCS() to work with newer flowCore package.
* Updated README.md and intro.Rmd vignette to reflect changes.
* Added a `NEWS.md` file to track changes to the package.
