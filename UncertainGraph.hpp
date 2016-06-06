//
//  UncertainGraph.hpp
//  testGraphCplus
//
//  Created by dongqingxiao on 4/11/16.
//  Copyright © 2016 dongqingxiao. All rights reserved.
//

#ifndef UncertainGraph_hpp
#define UncertainGraph_hpp

#include <stdio.h>
#include <igraph.h>
#include <iostream>
using namespace std;

class UncertainGraph{
public:
    igraph_t graph;
    long nv;
    long ne;
    igraph_vector_t pe;
protected:
    void degreeDistribution(igraph_vector_t * res, igraph_real_t maxDegree); // used for sigmaUniquess cal, it seems wrong thought
    void aggregateAK(igraph_vector_t * res, igraph_real_t maxDegree, igraph_vector_t * ak); // used for sigmaUniquess cal
    
    void aggregateReliablityDiff(igraph_vector_t * v_rep,igraph_vector_t * indicator,igraph_vector_t *res, double p);
    
    
    
public:
    UncertainGraph(long nv);      // empty constructor
    UncertainGraph(const UncertainGraph & obj); // copy const constructor
    UncertainGraph(UncertainGraph & obj); // copy constructor
    UncertainGraph& operator=(const UncertainGraph & obj); //
    ~UncertainGraph(); // destructor function
    
    void set_edges(igraph_vector_t * edges); // add edges to graph
    void set_edges_probs(igraph_vector_t * probs); // set edge probs to uncertain graph 
    
    void entropyReport(igraph_vector_t *res, igraph_real_t maxDegree);  // used only for uncertain graph, maxDegre=
    double testAgaist(igraph_vector_t * ak);
    void graphStastic(void); // basic graph statistics; linear function
    
    void degreeMetricRecord(string file);
    void graphDegreeVairance(double sigma_graph); //
    
    
    double diffconectPairAddEdge(double from, double to);
    void reliablityUtiitySubgraphCall(igraph_vector_t *ruv);
    void reliablityUtilitySubgraph(igraph_vector_t *ruv);
    void maxConnectedConponent(igraph_vector_t * re_edges, igraph_vector_t * re_pe, igraph_vector_t * startPos,long int * cnum);
    void reliablityUtiliy(igraph_vector_t * ruv); // used for uncertain graph to cal
    void reliablityUtiliyInit(igraph_vector_t * ruv); // read from file
    
    void reliablityUtiltyDiff(igraph_vector_t * ruv); // used for uncertail graph to cal
    
    
    void reliablity(igraph_vector_t * res); // used for uncertain for cal reliablity
    void reliablity(igraph_vector_t * res, string filePath); // used for storeing reliablity result into file
    void reliablity_record(long int sampleNum, string filePath);   // used for record relaiblity intermediate matrix =sample*nv
    void getDegrees(bool expected, igraph_vector_t *res);  // used for uncertain graph to extract adversary 
    
    
    
    
    
    void sigmaUniquess(igraph_vector_t * uv, igraph_vector_t ak, igraph_real_t maxDegree, igraph_real_t sigma);
    UncertainGraph generateObfuscation(igraph_real_t sigma, igraph_real_t * eps_res, igraph_vector_t * ak); // used for uncertain graph

    UncertainGraph randomGenerateObfuscation(igraph_real_t sigma, igraph_real_t * eps_res, igraph_vector_t * ak); // based on uniqueness val
    
    UncertainGraph greedyGenerateObfuscation(igraph_real_t sigma, igraph_real_t * eps_res, igraph_vector_t * ak); 
    
    
    
    UncertainGraph obfuscation(igraph_vector_t *ak, double * finalSigma); // used for uncertain graph obfuscation

   
    
    
    
    void aggregateReliablutyDiffEE(igraph_vector_t * ruv,igraph_vector_t *rue);
    UncertainGraph fixEdgeGraph(long int index,double val); // generate remove Graph via something subset
    UncertainGraph sampleGraph(igraph_vector_t * indicator); // generate sample graph with indicator
    UncertainGraph sampleGraph(void); // generate sample graph
    void print_graph(string filepath); // print graph to file
    
    void rawEstimate(igraph_vector_t * r_edge);
    void graphCheck();  // check input one is ok graph
};


UncertainGraph init_uncertain_from_file(string filepath); // used for init uncertain graph in edge format
UncertainGraph init_uncertain_OB_from_file(string filepath, long nv); // used for init uncertain graph in edge format with fixed number of vertices


#endif /* UncertainGraph_hpp */
