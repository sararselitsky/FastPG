// [[Rcpp::depends(BH)]]
#include <Rcpp.h>
#include <boost/container_hash/hash.hpp>
using namespace Rcpp;

struct Edge {
  std::tuple<long,long,double> edge;
  
  bool operator==(const Edge& e) const{
    long v1 = std::get<0>(this->edge);
    long v2 = std::get<1>(this->edge);
    long e_v1 = std::get<0>(e.edge);
    long e_v2 = std::get<1>(e.edge);
    return ((v1==e_v1)&&(v2==e_v2))
           ||
           ((v1==e_v2)&&(v2==e_v1));
  }
};

class edgeHashFn{
public:
  std::size_t operator()(const Edge& e) const{
    std::size_t seed = 0;
    if(std::get<0>(e.edge) < std::get<1>(e.edge)){
      boost::hash_combine(seed, std::get<0>(e.edge));
      boost::hash_combine(seed, std::get<1>(e.edge));
    }
    else{      
      boost::hash_combine(seed, std::get<1>(e.edge));
      boost::hash_combine(seed, std::get<0>(e.edge));
    }
    return seed;
  }
};

class edgeCompareFn{
public:
  std::size_t operator()(const Edge& e1,const Edge& e2) const{
    return (std::get<0>(e1.edge)==std::get<0>(e2.edge)
              &&
            std::get<1>(e1.edge)==std::get<1>(e2.edge))
           ||
           (std::get<0>(e1.edge)==std::get<1>(e2.edge)
              &&
            std::get<1>(e1.edge)==std::get<0>(e2.edge));
  }
};

//' Remove duplicate links
//'
//' @param links A numeric matrix of network edges
//' @return The matrix of edges with duplicates removed
//' @export
// [[Rcpp::export]]
NumericMatrix dedup_links(NumericMatrix links){
  std::unordered_set<Edge,edgeHashFn,edgeCompareFn> edgeSet;

  for(long i = 0; i<links.nrow();i++){
    if(links(i,0)!=0 && links(i,1)!=0 && links(i,2)!=0.0){
      Edge e = {std::make_tuple(links(i,0),links(i,1),links(i,2))};
      edgeSet.insert(e);
    }
  }

  NumericMatrix deduplicated_edges(edgeSet.size(),3);
  
  long row_index=0;
  for(Edge e:edgeSet){
    deduplicated_edges(row_index,0)=std::get<0>(e.edge);
    deduplicated_edges(row_index,1)=std::get<1>(e.edge);
    deduplicated_edges(row_index,2)=std::get<2>(e.edge);
    row_index++;
  }
  return deduplicated_edges;
}
