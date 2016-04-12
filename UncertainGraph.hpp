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
protected:
    igraph_vector_t pe;
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
    
    void reliablity(igraph_vector_t * res); // used for uncertain for cal reliablity
    void reliablity(igraph_vector_t * res, string filePath); // used for storeing reliablity result into file
    void getDegrees(bool expected, igraph_vector_t *res);  // used for uncertain graph to extract adversary 
    
    
    UncertainGraph sampleGraph(void); // generate sample graph
    void print_graph(string filepath);
    
    
};


UncertainGraph init_uncertain_from_file(string filepath); // used for init uncertain graph in edge format
UncertainGraph init_uncertain_OB_from_file(string filepath, long nv); // used for init uncertain graph in edge format with fixed number of vertices


#endif /* UncertainGraph_hpp */
