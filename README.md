# FastPG
Fast phenograph, CyTOF

[Under Construction]

Converting to a package; currently Just an R script and a cpp file; these must both be in the run directory.

Requires the following packages

* flowCore
* doMC
* data.table
* RANN
* Rcpp
* igraph
* nmslibR

### Usage:

From the directory that contains a copy of the fastPG.R script and the parallel_jc2.cpp source, run
```
./fastPG.R <file.fcs>
```