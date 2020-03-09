#include <Rcpp.h>
#include <omp.h>
#include<iostream>

#include "defs.h"
#include <unordered_map>
#include "basic_util.h"
#include "basic_comm.h"
#include "color_comm.h"
#include "sync_comm.h"

#include <set>
using namespace Rcpp;

void parse_SNAP(graph * G, NumericMatrix links, long num_vertices,std::unordered_map<long, long> &clusterLocalMap) {
 // printf("Parsing a SNAP formatted file as a general graph...\n");
 // printf("WARNING: Assumes that the graph is directed -- an edge is stored only once.\n");
 // printf("       : Graph will be stored as undirected, each edge appears twice.\n");
  int nthreads = 0;
#pragma omp parallel
{
  nthreads = omp_get_num_threads();
}
//printf("parse_SNAP: Number of threads: %d\n ", nthreads);

long   NV=0,  NE=0;
//string oneLine, myDelimiter(" "), myDelimiter2("\t"), oneWord; //Delimiter is a blank space
//char comment;

double time1, time2;
//ifstream fin;
//fin.open(fileName);
//if(!fin) {
//  cerr<<"Within Function: loadSNAPFileFormat() \n";
//  cerr<<"Could not open the file.. \n";
//  exit(1);
//}

/*do { //Parse the comment lines for problem size
  getline(fin, oneLine);
  cout<<"Read line: "<<oneLine<<endl;
  comment = oneLine[0];
  if (comment == '#') { //Check if this line has problem sizes
    StringTokenizer* ST = new StringTokenizer(oneLine, myDelimiter);
    if ( ST->HasMoreTokens() )
      oneWord = ST->GetNextToken(); //Ignore #
    if ( ST->HasMoreTokens() )
      oneWord = ST->GetNextToken(); //Ignore #
    if(oneWord == "Nodes:") {
      NV  = atol( ST->GetNextToken().c_str() ); //Number of Vertices
      oneWord = ST->GetNextToken(); //Ignore Edges:
      NE  = atol( ST->GetNextToken().c_str() ); //Number of Edges
    }
    delete ST;
  }
} while ( comment == '#');*/

NE = links.nrow();
NV = num_vertices;

//printf("|V|= %ld, |E|= %ld \n", NV, NE);
//printf("Weights will read from the file.\n");
//cout << oneLine <<endl;
/*---------------------------------------------------------------------*/
/* Read edge list: a U V W                                             */
/*---------------------------------------------------------------------*/
edge *tmpEdgeList = (edge *) malloc( NE * sizeof(edge)); //Every edge stored ONCE
assert( tmpEdgeList != NULL);
long Si, Ti;

//std::unordered_map<long, long> clusterLocalMap; //Renumber vertices contiguously from zero
//std::unordered_map<long, long>::iterator storedAlready;
std::unordered_map<long,long>::iterator storedAlready;
long numUniqueVertices = 0;

//Parse the first edge already read from the file and stored in oneLine
//long i=0;
double wt = 1.0;
//do {
for(long i = 0;i<NE;i++){
  //StringTokenizer* ST = new StringTokenizer(oneLine, myDelimiter2);
  //if ( ST->HasMoreTokens() )
  //  Si  = atol( ST->GetNextToken().c_str() );
  //if ( ST->HasMoreTokens() )
  //  Ti  = atol( ST->GetNextToken().c_str() );
  //if ( ST->HasMoreTokens() )
  //  wt  = atof( ST->GetNextToken().c_str() );
  //delete ST;
  Si = links(i,0);
  Ti = links(i,1);
  wt = links(i,2);
  
  storedAlready = clusterLocalMap.find(Si); //Check if it already exists
  if( storedAlready != clusterLocalMap.end() ) {	//Already exists
    Si = storedAlready->second; //Renumber the cluster id
  } else {
    clusterLocalMap[Si] = numUniqueVertices; //Does not exist, add to the map
    Si = numUniqueVertices; //Renumber the vertex id
    numUniqueVertices++; //Increment the number
  }
  
  storedAlready = clusterLocalMap.find(Ti); //Check if it already exists
  if( storedAlready != clusterLocalMap.end() ) {	//Already exists
    Ti = storedAlready->second; //Renumber the cluster id
  } else {
    clusterLocalMap[Ti] = numUniqueVertices; //Does not exist, add to the map
    Ti = numUniqueVertices; //Renumber the vertex id
    numUniqueVertices++; //Increment the number
  }
  tmpEdgeList[i].head   = Si;  //The S index
  tmpEdgeList[i].tail   = Ti;  //The T index: One-based indexing
  tmpEdgeList[i].weight = wt;     //default weight of one
  //cout<<" Adding edge ("<<Si<<", "<<Ti<<")\n";
}
//  i++;
  //Read-in the next line
//  getline(fin, oneLine);
//  if ((i % 99999) == 1) {
//    cout <<"Reading Line: "<<i<<endl;
//  }
//} while ( !fin.eof() );//End of while

//fin.close(); //Close the file
time2 = omp_get_wtime();
//printf("Done reading from file: NE= %ld. Time= %lf\n", NE, time2-time1);
//printf("Number of unique vertices: %ld \n", numUniqueVertices);

NV = numUniqueVertices;
///////////
time1 = omp_get_wtime();
long *edgeListPtr = (long *)  malloc((NV+1) * sizeof(long));
assert(edgeListPtr != NULL);
edge *edgeList = (edge *) malloc( 2*NE * sizeof(edge)); //Every edge stored twice
assert( edgeList != NULL);
time2 = omp_get_wtime();
//printf("Time for allocating memory for storing graph = %lf\n", time2 - time1);
#pragma omp parallel for
for (long i=0; i <= NV; i++)
  edgeListPtr[i] = 0; //For first touch purposes

//////Build the EdgeListPtr Array: Cumulative addition
time1 = omp_get_wtime();
#pragma omp parallel for
for(long i=0; i<NE; i++) {
  __sync_fetch_and_add(&edgeListPtr[tmpEdgeList[i].head+1], 1); //Leave 0th position intact
  __sync_fetch_and_add(&edgeListPtr[tmpEdgeList[i].tail+1], 1);
}
for (long i=0; i<NV; i++) {
  edgeListPtr[i+1] += edgeListPtr[i]; //Prefix Sum:
}
//The last element of Cumulative will hold the total number of characters
time2 = omp_get_wtime();
//printf("Done cumulative addition for edgeListPtrs:  %9.6lf sec.\n", time2 - time1);
//printf("Sanity Check: 2|E| = %ld, edgeListPtr[NV]= %ld\n", NE*2, edgeListPtr[NV]);
//printf("*********** (%ld)\n", NV);

//printf("About to build edgeList...\n");
time1 = omp_get_wtime();
//Keep track of how many edges have been added for a vertex:
long  *added  = (long *)  malloc( NV  * sizeof(long)); assert( added != NULL);
#pragma omp parallel for
for (long i = 0; i < NV; i++)
  added[i] = 0;
//printf("...\n");
//Build the edgeList from edgeListTmp:
#pragma omp parallel for
for(long i=0; i<NE; i++) {
  long head      = tmpEdgeList[i].head;
  long tail      = tmpEdgeList[i].tail;
  double weight  = tmpEdgeList[i].weight;
  
  long Where = edgeListPtr[head] + __sync_fetch_and_add(&added[head], 1);
  edgeList[Where].head = head;
  edgeList[Where].tail = tail;
  edgeList[Where].weight = weight;
  //Now add the counter-edge:
  Where = edgeListPtr[tail] + __sync_fetch_and_add(&added[tail], 1);
  edgeList[Where].head = tail;
  edgeList[Where].tail = head;
  edgeList[Where].weight = weight;    
}
time2 = omp_get_wtime();
//printf("Time for building edgeList = %lf\n", time2 - time1);

///////Store the vertex ids in a file////////
/*char filename2[256];
sprintf(filename2,"%s_vertexMap.txt", fileName);
//printf("Writing vertex map (new id -- old id) in file: %s\n", filename2);
FILE *fout;
fout = fopen(filename2, "w");
if (!fout) {
 // printf("Could not open the file \n");
  exit(1);
}*/
//Write the edges (lower triangle only):
//storedAlready = clusterLocalMap.begin();
//FILE *fout;
//fout = fopen("clustInfo.txt", "w");
//while (storedAlready != clusterLocalMap.end()) {
//  fprintf(fout, "%ld %ld\n", storedAlready->first, storedAlready->second);
//  storedAlready++;
//}
//fclose(fout);
//printf("Vertex map has been stored in file: %s\n",filename2);

G->sVertices    = NV;
G->numVertices  = NV;
G->numEdges     = NE;
G->edgeListPtrs = edgeListPtr;
G->edgeList     = edgeList;

//Clean up
free(tmpEdgeList);
free(added);
}//End of parse_SNAP()

double find_communities(graph * G, long* C_orig) {
  
  int nT = 1; //Default is one thread
#pragma omp parallel
{
  nT = omp_get_num_threads();
}
if (nT < 1) {
 // printf("The number of threads should be greater than one.\n");
  return 0;
}

// File Loading
double time1, time2;
//graph* G = (graph *) malloc (sizeof(graph));
double final_modularity = -1; 

//displayGraphCharacteristics(G);
//not sure, but looks like
// threadsOpt is not currently a
// real opt
/* orig code:
int threadsOpt = 0;
if(opts.threadsOpt)
threadsOpt =1;
threadsOpt =1;*/
int threadsOpt=1;

/* orig code
int replaceMap = 0;
if(  opts.basicOpt == 1 )
replaceMap = 1 
*/

int replaceMap = 1;

bool VF = true;

/* Vertex Following option */
//if( opts.VF ) {
if (VF){
 // printf("Vertex following is enabled.\n");
  time1 = omp_get_wtime();
  long numVtxToFix = 0; //Default zero
  long *C = (long *) malloc (G->numVertices * sizeof(long)); assert(C != 0);
  numVtxToFix = vertexFollowing(G,C); //Find vertices that follow other vertices
  if( numVtxToFix > 0) {  //Need to fix things: build a new graph
   // printf("Graph will be modified -- %ld vertices need to be fixed.\n", numVtxToFix);
    graph *Gnew = (graph *) malloc (sizeof(graph));
    long numClusters = renumberClustersContiguously(C, G->numVertices);
    buildNewGraphVF(G, Gnew, C, numClusters);
    //Get rid of the old graph and store the new graph
    free(G->edgeListPtrs);
    free(G->edgeList);
    free(G);
    G = Gnew;
  }
  free(C); //Free up memory
 // printf("Graph after modifications:\n");
 // displayGraphCharacteristics(G);
}//End of if( VF == 1 )


// Datastructures to store clustering information
long NV = G->numVertices;
//printf("NV (G->numVertices) = %ld",NV);
//long *C_orig = (long *) malloc (NV * sizeof(long)); assert(C_orig != 0);
//int *C_ints = (int *) malloc (NV * sizeof(int)); assert(C_ints != 0);

std::set<long> clustkeys;

graph* G_orig = (graph *) malloc (sizeof(graph)); //The original version of the graph
duplicateGivenGraph(G, G_orig);

//Call the clustering algorithm:
//Call the clustering algorithm:
bool strongScaling = false;
int coloring = 1;
int syncType = 0;
int numColors = 16;
double C_thresh = 0.01;
long minGraphSize = 100000;
double threshold = 0.000001;
int basicOpt = 1;
//if ( opts.strongScaling ) { //Strong scaling enabled
if(strongScaling){
  //Retain the original copy of the graph:
  graph* G_original = (graph *) malloc (sizeof(graph)); //The original version of the graph
  time1 = omp_get_wtime();
  duplicateGivenGraph(G, G_original);
  time2 = omp_get_wtime();
 // printf("Time to duplicate : %lf\n", time2-time1);
  
  //Run the algorithm in powers of two for the maximum number of threads available
  int curThread = 2; //Start with two threads
  while (curThread <= nT) {
   // printf("\n\n***************************************\n");
   // printf("Starting run with %d threads.\n", curThread);
   // printf("***************************************\n");
    //Call the clustering algorithm:
#pragma omp parallel for
    for (long i=0; i<G->numVertices; i++) {
      C_orig[i] = -1;
    }
    //if(opts.coloring != 0){
    if(coloring != 0) {
      runMultiPhaseColoring(G, C_orig, coloring, numColors, replaceMap, minGraphSize, threshold, C_thresh, curThread, threadsOpt);
    }else if(syncType != 0){
      runMultiPhaseSyncType(G, C_orig, syncType, minGraphSize, threshold, C_thresh, curThread, threadsOpt);
    }else{
      runMultiPhaseBasic(G, C_orig, basicOpt, minGraphSize, threshold, C_thresh, curThread,threadsOpt);
    }
    //Increment thread and revert back to original graph
    if (curThread < nT) {
      //Skip copying at the very end
      //Old graph is already destroyed in the above function
      G = (graph *) malloc (sizeof(graph)); //Allocate new space
      duplicateGivenGraph(G_original, G); //Copy the original graph to G
    }
    curThread = curThread*2; //Increment by powers of two
  }//End of while()
} else { //No strong scaling -- run once with max threads

 
#pragma omp parallel for
  for (long i=0; i<NV; i++) {
    C_orig[i] = -1;
  }
  //if(opts.coloring != 0){
  if(coloring != 0) {
    final_modularity = runMultiPhaseColoring(G, C_orig, coloring, numColors, replaceMap, minGraphSize, threshold, C_thresh, nT, threadsOpt);
    //}else if(opts.syncType != 0){
  }else if(syncType != 0){
    runMultiPhaseSyncType(G, C_orig, syncType, minGraphSize, threshold, C_thresh, nT,threadsOpt);
  }else{
    runMultiPhaseBasic(G, C_orig, basicOpt, minGraphSize, threshold, C_thresh, nT,threadsOpt);
  }
}


/* Tom - 1/17/2020
 Will this introduce a memory leak?
 //Cleanup:
 if(C_orig != 0) free(C_orig);
 //Do not free G here -- it will be done in another routine.
 */
//return C_orig;
//int *C_ints = (int *) malloc (NV * sizeof(int)); assert(C_ints != 0);

//C_ints = (int) C_orig;

//for(int i=0;i<NV;i++){
//  clustkeys.insert(C_orig[i]);
//}

//for(auto ck : clustkeys){
// // printf("%ld \n",ck);
//}
//FILE* out = fopen("c_orig.txt","w");
//for(long i = 0; i<NV;i++) {
//   fprintf(out,"%ld\n",C_orig[i]);
//}
//fclose(out);
//return NumericVector(C_orig,C_orig+(sizeof(C_orig)/sizeof(*C_orig)));
//return NumericVector(C_ints,C_ints+NV);
//return C_orig;
return final_modularity;
}//End of main()


//' Parallel Louvain clustering
//'
//' @param links A numeric matrix of network edges
//' @param num_vertices The number of nodes
//' @return A list with two elements
//' * `modularity` - 
//' * `communities` - A vector where the i'th value is the
//'   cluster number that the i'th node in the links matrix has been assigned to.
//' @export
// [[Rcpp::export]]
Rcpp::List parallel_louvain(NumericMatrix links, long num_vertices){
  NumericVector res(num_vertices);
  double modularity = -1;
  graph* G = (graph *) malloc (sizeof(graph));
  
  std::unordered_map<long,long> clusterLocalMap;
  
  parse_SNAP(G,links,num_vertices,clusterLocalMap);
  
  long *C_orig = (long *) malloc (num_vertices * sizeof(long)); assert(C_orig != 0);
  
  modularity = find_communities(G,C_orig);
  
  for(auto it = clusterLocalMap.begin();it != clusterLocalMap.end(); ++it){
    res[it->first-1]=(int)C_orig[it->second];
  }
  
  //printf("%d is the length of res 1\n",res.length());

  return Rcpp::List::create(Rcpp::Named("modularity")=modularity,
                            Rcpp::Named("communities")=res);
}
