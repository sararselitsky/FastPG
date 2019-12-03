#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::depends(RcppParallel)]]
#include <RcppParallel.h>
#include <vector>
#include <algorithm>
#include <limits>
#include <cmath>

using namespace RcppParallel;

auto removeNANs = [](double number) -> bool
{
  return std::isnan(number);
};

struct Jce : public Worker {

  const RMatrix<double> mat;

  RMatrix<double> weights;

  const int ncols;

  Jce(const NumericMatrix mat, NumericMatrix weights)
    : mat(mat), weights(weights), ncols(mat.ncol()) {}

  void operator()(std::size_t begin, std::size_t end){
    for (std::size_t i = begin; i < end; i++) {
         for (std::size_t j = 0; j < ncols; j++) {

           RMatrix<double>::Row rowi = mat.row(i);
           if(mat(i,j)==mat(i,j)){
           int k = mat(i,j) - 1;
           RMatrix<double>::Row rowk = mat.row(k);

           std::vector<double> r1(rowi.length());
           std::vector<double> r2(rowk.length());
           std::vector<double> r3;
           std::copy(rowi.begin(),rowi.end(),r1.begin());
           std::copy(rowk.begin(),rowk.end(),r2.begin());

           std::vector<double> intersections;
           std::vector<double> u(r1.size()+r2.size());

           std::sort(r1.begin(),r1.end());
           std::sort(r2.begin(),r2.end());

           std::vector<double>::iterator r1end;
           std::vector<double>::iterator r2end;
           std::vector<double>::iterator uend;
           r1end = std::remove_if(r1.begin(),r1.end(),removeNANs);
           r2end = std::remove_if(r2.begin(),r2.end(),removeNANs);

           r1.resize(r1end-r1.begin());
           r2.resize(r2end-r2.begin());

           std::set_intersection(r1.begin(),r1.end(),
                                 r2.begin(),r2.end(),
                                 std::back_inserter(r3));

           weights((i*ncols)+j, 0) = i+1;
           weights((i*ncols)+j, 1) = k+1;
           weights((i*ncols)+j, 2) = ((double)r3.size()/(r1.size() + r2.size() - r3.size()))/2.0;

          }
         }
    }
  }
};
// [[Rcpp::export]]
NumericMatrix rcpp_parallel_jce(NumericMatrix mat) {

  // allocate the matrix we will return
  NumericMatrix rmat(mat.nrow()*mat.ncol(),3);

  // create the worker
  Jce jce(mat, rmat);

  // call it with parallelFor
  parallelFor(0, mat.nrow(), jce);

  return rmat;
}
