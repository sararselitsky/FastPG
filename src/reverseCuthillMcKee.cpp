/*---------------------------------------------------------------------------*/
/*                                                                           */
/*                          Mahantesh Halappanavar                           */
/*                        High Performance Computing                         */
/*                Pacific Northwest National Lab, Richland, WA               */
/*                                                                           */
/*---------------------------------------------------------------------------*/
/*                                                                           */
/* Copyright (C) 2010 Mahantesh Halappanavar                                 */
/*                                                                           */
/*                                                                           */
/* This program is free software; you can redistribute it and/or             */
/* modify it under the terms of the GNU General Public License               */
/* as published by the Free Software Foundation; either version 2            */
/* of the License, or (at your option) any later version.                    */
/*                                                                           */
/* This program is distributed in the hope that it will be useful,           */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of            */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             */
/* GNU General Public License for more details.                              */
/*                                                                           */
/* You should have received a copy of the GNU General Public License         */
/* along with this program; if not, write to the Free Software               */
/* Foundation, Inc., 59 Temple Place-Suite 330,Boston,MA 02111-1307,USA.     */
/*                                                                           */
/*---------------------------------------------------------------------------*/

#include "defs.h"
#include "dataStructureHeap.h"

// Perform reverse Cuthill-McKee operation on the graph
// SSize indicates the size of Source
// Sets the pointer in isChordal from one direction ONLY
void algoReverseCuthillMcKee( graph *G, long *pOrder, int nThreads )
{
    printf("Within algoReverseCuthillMcKee() \n");
    if (nThreads < 1)
        omp_set_num_threads(1);
    else
        omp_set_num_threads(nThreads);
    int nT;
#pragma omp parallel
    {
        nT = omp_get_num_threads();
    }
    
    double time1=0, time2=0, total=0, totalTime=0;
    long    NV        = G->numVertices;
    long    NS        = G->sVertices;
    long    NT        = NV - NS;
    bool    isSym     = true;
    if(NT > 0)
        isSym = false; //A bipartite graph
    long    NE        = G->numEdges;
    long    *vtxPtr   = G->edgeListPtrs;
    edge    *vtxInd   = G->edgeList;
    printf("Vertices:%ld  Edges:%ld\n", NV, NE);
    
    //////STEP 1: Sort the vertices in order of their degree
    
    //Compute the degree of each vertex:
    time1 = omp_get_wtime();
    long *degree  = (long *) malloc (NV * sizeof(long)); assert(degree != 0);
#pragma omp parallel for
    for (long i=0; i<NV; i++) {
        degree[i] = vtxPtr[i+1] - vtxPtr[i];
    }
    //Add the vertices to a heap data structure:
    heap *myHeap = (heap*) malloc (sizeof(heap));
    heapInitializeToN(myHeap, NV);
    printf("Max heap size: %ld\n", myHeap->maxsize);
    term newTerm;
    
    long *visited  = (long *) malloc (NV * sizeof(long)); assert(visited != 0);
    long *R = (long *) malloc (NV * sizeof(long)); assert(R != 0);
    //Initialize the Vectors:
#pragma omp parallel for
    for (long i=0; i<NV; i++) {
        visited[i]= 0; //zero means not visited
        R[i] = -1; //Initialize the rank as -1
    }
    long howManyAdded = 0; //How many vertices have been added to the queue
    //Do not add isolated (degree=0) vertices:
    for (int i=0; i<NV; i++) {
        if(degree[i] > 0) {
            newTerm.id = i;
            newTerm.weight = degree[i];
            heapAdd(myHeap, newTerm); //Has to be added in serial
        } else {
            R[howManyAdded] = i; //Add them directly to the stack
            visited[i] = 1;
            howManyAdded++;
        }
    }
    printf("Heap size: %ld\n", myHeap->size);
    free(degree);
    time2 = omp_get_wtime();
    printf("Time for adding %ld elements to the heap: %lf\n",NV, time2-time1);
    totalTime += time2-time1;
    
    ////////STEP 2: Now perform the BFS
    
    
    //The Queue Data Structure for the Dominating Set:
    //The Queues are important for synchornizing the concurrency:
    //Have two queues - read from one, write into another
    // at the end, swap the two.
    long *Q    = (long *) malloc (NV * sizeof(long)); assert(Q != 0);
    long *Qtmp = (long *) malloc (NV * sizeof(long)); assert(Qtmp != 0);
    long *Qswap;
    
#pragma omp parallel for
    for (long i=0; i<NV; i++) {
        Q[i]= -1;
        Qtmp[i]= -1;
    }
    long QTail=0; //Tail of the queue (implicitly will represent the size)
    long QtmpTail=0; //Tail of the queue (implicitly will represent the size)
    
    //Get the smallest degree as the first source vertex
    term *data = myHeap->elements;
    heapRemoveMin(myHeap);     //Remove it from the heap
    printf("Source vertex: %ld\n", data[0].id);
    Q[0] = data[0].id; //Add the smallest vertex to the queue
    QTail = 1; //Increment the queue tail
    visited[data[0].id] = 1; //Mark the vertex as visited
    R[howManyAdded] = data[0].id; //Enter the vertex in the vector
    howManyAdded++;
    printf("R[0] = %ld\n", data[0].id);
    
    long nCC = 1;
    
    while (howManyAdded < NV) { //Process until all the vertices have been added to the queue
        printf("***********************************************\n");
        printf("Connected component       : %d     \n",    nCC);
        printf("***********************************************\n");
        //The size of Q1 is now QTail+1; the elements are contained in Q1[0] through Q1[Q1Tail]
        int nLoops=0; //Count number of iterations in the while loop
        while ( QTail > 0 ) {
            printf("Loop %d, QSize= %ld\n",nLoops,QTail);
            //KEY IDEA: Process all the members of the queue concurrently:
            time1 = omp_get_wtime();
#pragma omp parallel for
            for (long Qi=0; Qi<QTail; Qi++) {
                long v = Q[Qi];
                long adj1 = vtxPtr[v];
                long adj2 = vtxPtr[v+1];
                for(long k = adj1; k < adj2; k++) {
                    long x = vtxInd[k].tail;
                    //Has this neighbor been visited?
                    if ( __sync_fetch_and_add(&visited[x], 1) == 0 ) {
                        //Not visited: add it to the Queue
                        long whereInQ = __sync_fetch_and_add(&QtmpTail, 1);
                        Qtmp[whereInQ] = x;
                        long whereInR = __sync_fetch_and_add(&howManyAdded, 1);
                        //printf("R[%ld] = %ld\n", whereInR, x);
                        R[whereInR] = x; //Add it to the ranked list
                    }
                } //End of for loop on k: the neighborhood of v
            } //End of for(Qi)
            // Also end of the parallel region
            // Swap the two queues:
            Qswap = Q;
            Q = Qtmp; //Q now points to the second vector
            Qtmp = Qswap;
            QTail = QtmpTail; //Number of elements
            QtmpTail = 0; //Symbolic emptying of the second queue
            nLoops++;
            time2  = omp_get_wtime();
            total += time2-time1;
        } //end of while ( !Q.empty() )
        printf("***********************************************\n");
        printf("Number of iterations       : %d     \n", nLoops);
        printf("Time for Graph-traversal   : %lf sec\n", total);
        printf("Number of elements added   : %ld    \n", howManyAdded);
        printf("***********************************************\n");
        
        /////Now look for the next smallest id:
        bool found = 1;
        if(howManyAdded < NV) {
            printf("Looking for the next source vertex...(Heap size = %ld)\n", myHeap->size);
            term *data;
            do {
                data = myHeap->elements;
                //printf("Processing vertex: %ld Visited: %ld -- (Heap size = %ld)\n", data[0].id+1, visited[data[0].id], myHeap->size);
                heapRemoveMin(myHeap);
                if (myHeap->size == 0) { //No more elements to look at
                    found = 0;
                    printf("Quitting because there is nothing in the heap...(Heap size = %ld)\n", myHeap->size);
                    break;
                }
                //printf("Id: %ld  -- status: %ld\n", data[0].id, visited[data[0].id]);
            } while(visited[data[0].id] > 0); //Check if it has already been visited
            if (found == 0) {
                break; //break from the main loop
            } else { //Add the new source to the list
                printf("New source vertex: %ld\n", data[0].id);
                Q[0] = data[0].id; //Add the smallest vertex to the queue
                QTail = 1; //Increment the queue tail
                visited[data[0].id] = 1; //Mark the vertex as visited
                R[howManyAdded] = data[0].id; //Enter the vertex in the vector
                //printf("R[%ld] = %ld\n", howManyAdded, data[0].id);
                howManyAdded++; //Increment the #vertices that have been added to the stack
            }
        }//End of if()
        nCC++;
        totalTime += total;
        total = 0;
    }
    //Clean Up:
    free(Q);
    free(Qtmp);
    free(visited);
    free(myHeap->elements);
    free(myHeap);
    
    assert(howManyAdded == NV); //Sanity check before moving to next step
    //////STEP 3: Received a valid vector; reverse the order:
    if (isSym) { //A symmetric matrix
        for (long i=0; i<NV; i++) {
            //pOrder[i]= R[NV - i - 1];
            pOrder[R[i]]= NV - i - 1; //pOrder is a old2New index mapping
        }
    } else { //A bipartite graph:
        //STEP 3.1: Segregate the row and column vertices
        long rowCounter = 0;
        long colCounter = NS;
        long *Rprime    = (long *) malloc (NV * sizeof(long)); assert(Rprime != 0);
        for (long i=0; i<NV; i++) {
            Rprime[i]= -1;
        }
        for (long i=(NV-1); i>=0; i--) { //Go through the list in a reverse order
            if(R[i] < NS) { //A row vertex
                Rprime[rowCounter] = R[i];
                rowCounter++;
            } else { //A column vertex
                Rprime[colCounter] = R[i];
                colCounter++;
            }
        }//End of for(i)
        assert(rowCounter==NS); assert(colCounter==NV); //Sanity check
        //STEP 3.2: Now build the old2New map:
        for (long i=0; i<NV; i++) {
            pOrder[Rprime[i]] = i; //pOrder is a old2New index mapping
        }
        //Clean up:
        free(Rprime);
    }//End of else(bipartite graph)
    
    printf("***********************************************\n");
    printf("Number of connected components       : %d     \n", nCC);
    printf("Total Time                           : %lf sec\n", totalTime);
    printf("***********************************************\n");
    
    //Clean Up:
    free(R);
    
} //End of algoReverseCuthillMcKee

// Perform reverse Cuthill-McKee operation on the graph: Strict variant
// SSize indicates the size of Source
// Sets the pointer in isChordal from one direction ONLY
void algoReverseCuthillMcKeeStrict( graph *G, long *pOrder, int nThreads )
{
    printf("Within algoReverseCuthillMcKee() \n");
    if (nThreads < 1)
        omp_set_num_threads(1);
    else
        omp_set_num_threads(nThreads);
    int nT;
#pragma omp parallel
    {
        nT = omp_get_num_threads();
    }
    
    double time1=0, time2=0, total=0, totalTime=0;
    long    NV        = G->numVertices;
    long    NS        = G->sVertices;
    long    NT        = NV - NS;
    bool    isSym     = true;
    if(NT > 0)
        isSym = false; //A bipartite graph
    long    NE        = G->numEdges;
    long    *vtxPtr   = G->edgeListPtrs;
    edge    *vtxInd   = G->edgeList;
    printf("Vertices:%ld  Edges:%ld\n", NV, NE);
    
    //////STEP 1: Sort the vertices in order of their degree
    
    //Compute the degree of each vertex:
    time1 = omp_get_wtime();
    long *degree  = (long *) malloc (NV * sizeof(long)); assert(degree != 0);
#pragma omp parallel for
    for (long i=0; i<NV; i++) {
        degree[i] = vtxPtr[i+1] - vtxPtr[i];
    }
    //Add the vertices to a heap data structure:
    heap *myHeap = (heap*) malloc (sizeof(heap));
    heapInitializeToN(myHeap, NV);
    printf("Max heap size: %ld\n", myHeap->maxsize);
    term newTerm;
    
    long *visited  = (long *) malloc (NV * sizeof(long)); assert(visited != 0);
    long *R = (long *) malloc (NV * sizeof(long)); assert(R != 0);
    long *oneLevel = (long *) malloc (NV * sizeof(long)); assert(oneLevel != 0);
    //Initialize the Vectors:
#pragma omp parallel for
    for (long i=0; i<NV; i++) {
        visited[i]= 0; //zero means not visited
        R[i] = -1; //Initialize the rank as -1
        oneLevel[i] = -1; //Hold vertices at a given level
    }
    long howManyAdded = 0; //How many vertices have been added to the queue
    long howManyAddedLevel = 0;
    //Do not add isolated (degree=0) vertices:
    for (int i=0; i<NV; i++) {
        if(degree[i] > 0) {
            newTerm.id = i;
            newTerm.weight = degree[i];
            heapAdd(myHeap, newTerm); //Has to be added in serial
        } else {
            R[howManyAdded] = i; //Add them directly to the stack
            visited[i] = 1;
            howManyAdded++;
        }
    }
    printf("Heap size: %ld\n", myHeap->size);
    free(degree);
    time2 = omp_get_wtime();
    printf("Time for adding %ld elements to the heap: %lf\n",NV, time2-time1);
    totalTime += time2-time1;
    
    ////////STEP 2: Now perform the BFS
    
    
    //The Queue Data Structure for the Dominating Set:
    //The Queues are important for synchornizing the concurrency:
    //Have two queues - read from one, write into another
    // at the end, swap the two.
    long *Q    = (long *) malloc (NV * sizeof(long)); assert(Q != 0);
    long *Qtmp = (long *) malloc (NV * sizeof(long)); assert(Qtmp != 0);
    long *Qswap;
    
#pragma omp parallel for
    for (long i=0; i<NV; i++) {
        Q[i]= -1;
        Qtmp[i]= -1;
    }
    long QTail=0; //Tail of the queue (implicitly will represent the size)
    long QtmpTail=0; //Tail of the queue (implicitly will represent the size)
    
    //Get the smallest degree as the first source vertex
    term *data = myHeap->elements;
    heapRemoveMin(myHeap);     //Remove it from the heap
    printf("Source vertex: %ld\n", data[0].id);
    Q[0] = data[0].id; //Add the smallest vertex to the queue
    QTail = 1; //Increment the queue tail
    visited[data[0].id] = 1; //Mark the vertex as visited
    R[howManyAdded] = data[0].id; //Enter the vertex in the vector
    howManyAdded++;
    printf("R[0] = %ld\n", data[0].id);
    
    long nCC = 1;
    
    while (howManyAdded < NV) { //Process until all the vertices have been added to the queue
        printf("***********************************************\n");
        printf("Connected component       : %d     \n",    nCC);
        printf("***********************************************\n");
        //The size of Q1 is now QTail+1; the elements are contained in Q1[0] through Q1[Q1Tail]
        int nLoops=0; //Count number of iterations in the while loop
        while ( QTail > 0 ) {
            printf("Loop %d, QSize= %ld\n",nLoops,QTail);
            //KEY IDEA: Process all the members of the queue concurrently:
            time1 = omp_get_wtime();
            //#pragma omp parallel for
            for (long Qi=0; Qi<QTail; Qi++) {
                long v = Q[Qi];
                long adj1 = vtxPtr[v];
                long adj2 = vtxPtr[v+1];
                for(long k = adj1; k < adj2; k++) {
                    long x = vtxInd[k].tail;
                    //Has this neighbor been visited?
                    if ( __sync_fetch_and_add(&visited[x], 1) == 0 ) {
                        //Not visited: add it to the Queue
                        long whereInQ = __sync_fetch_and_add(&QtmpTail, 1);
                        Qtmp[whereInQ] = x;
                        //printf("R[%ld] = %ld\n", whereInR, x);
                    }
                } //End of for loop on k: the neighborhood of v
            } //End of for(Qi)
            // Also end of the parallel region
            //Now add all the vertices into R in a sorted order:
            heap *myTmpHeap = (heap*) malloc (sizeof(heap));
            heapInitializeToN(myTmpHeap, QtmpTail);
            for (int Ti=0; Ti<QtmpTail; Ti++) {
                newTerm.id = Qtmp[Ti];
                newTerm.weight = degree[Qtmp[Ti]];
                heapAdd(myTmpHeap, newTerm); //Has to be added in serial
            }
            //printf("Added %ld vertices to the queue\n", QtmpTail);
            term *tmpData;
            for (int Ti=0; Ti<QtmpTail; Ti++) {
                //printf("Removing %ld id.\n", Ti+1);
                tmpData = myTmpHeap->elements;
                heapRemoveMin(myTmpHeap);
                long whereInR = __sync_fetch_and_add(&howManyAdded, 1);
                R[whereInR] = tmpData[0].id; //Add it to the ranked list
            }
            free(myTmpHeap->elements);
            free(myTmpHeap);
            // Swap the two queues:
            Qswap = Q;
            Q = Qtmp; //Q now points to the second vector
            Qtmp = Qswap;
            QTail = QtmpTail; //Number of elements
            QtmpTail = 0; //Symbolic emptying of the second queue
            nLoops++;
            time2  = omp_get_wtime();
            total += time2-time1;
        } //end of while ( !Q.empty() )
        printf("***********************************************\n");
        printf("Number of iterations       : %d     \n", nLoops);
        printf("Time for Graph-traversal   : %lf sec\n", total);
        printf("Number of elements added   : %ld    \n", howManyAdded);
        printf("***********************************************\n");
        
        /////Now look for the next smallest id:
        bool found = 1;
        if(howManyAdded < NV) {
            printf("Looking for the next source vertex...(Heap size = %ld)\n", myHeap->size);
            term *data;
            do {
                data = myHeap->elements;
                //printf("Processing vertex: %ld Visited: %ld -- (Heap size = %ld)\n", data[0].id+1, visited[data[0].id], myHeap->size);
                heapRemoveMin(myHeap);
                if (myHeap->size == 0) { //No more elements to look at
                    found = 0;
                    printf("Quitting because there is nothing in the heap...(Heap size = %ld)\n", myHeap->size);
                    break;
                }
                //printf("Id: %ld  -- status: %ld\n", data[0].id, visited[data[0].id]);
            } while(visited[data[0].id] > 0); //Check if it has already been visited
            if (found == 0) {
                break; //break from the main loop
            } else { //Add the new source to the list
                printf("New source vertex: %ld\n", data[0].id);
                Q[0] = data[0].id; //Add the smallest vertex to the queue
                QTail = 1; //Increment the queue tail
                visited[data[0].id] = 1; //Mark the vertex as visited
                R[howManyAdded] = data[0].id; //Enter the vertex in the vector
                //printf("R[%ld] = %ld\n", howManyAdded, data[0].id);
                howManyAdded++; //Increment the #vertices that have been added to the stack
            }
        }//End of if()
        nCC++;
        totalTime += total;
        total = 0;
    }
    //Clean Up:
    free(Q);
    free(Qtmp);
    free(visited);
    free(myHeap->elements);
    free(myHeap);
    
    assert(howManyAdded == NV); //Sanity check before moving to next step
    //////STEP 3: Received a valid vector; reverse the order:
    if (isSym) { //A symmetric matrix
        for (long i=0; i<NV; i++) {
            //pOrder[i]= R[NV - i - 1];
            pOrder[R[i]]= NV - i - 1; //pOrder is a old2New index mapping
        }
    } else { //A bipartite graph:
        //STEP 3.1: Segregate the row and column vertices
        long rowCounter = 0;
        long colCounter = NS;
        long *Rprime    = (long *) malloc (NV * sizeof(long)); assert(Rprime != 0);
        for (long i=0; i<NV; i++) {
            Rprime[i]= -1;
        }
        for (long i=(NV-1); i>=0; i--) { //Go through the list in a reverse order
            if(R[i] < NS) { //A row vertex
                Rprime[rowCounter] = R[i];
                rowCounter++;
            } else { //A column vertex
                Rprime[colCounter] = R[i];
                colCounter++;
            }
        }//End of for(i)
        assert(rowCounter==NS); assert(colCounter==NV); //Sanity check
        //STEP 3.2: Now build the old2New map:
        for (long i=0; i<NV; i++) {
            pOrder[Rprime[i]] = i; //pOrder is a old2New index mapping
        }
        //Clean up:
        free(Rprime);
    }//End of else(bipartite graph)
    
    printf("***********************************************\n");
    printf("Number of connected components       : %d     \n", nCC);
    printf("Total Time                           : %lf sec\n", totalTime);
    printf("***********************************************\n");
    
    //Clean Up:
    free(R);
    
} //End of algoEdgeApproxDominatingEdgesLinearSearch

