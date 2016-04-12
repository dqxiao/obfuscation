//
//  Graph.hpp
//  testGraphCplus
//
//  Created by dongqingxiao on 4/11/16.
//  Copyright © 2016 dongqingxiao. All rights reserved.
//

#ifndef Graph_hpp
#define Graph_hpp

#include <stdio.h>
#include <igraph.h>
#include <string>
#include <Configuration.hpp>
#include <UncertainGraph.hpp>
using namespace std;


class Graph{
protected:
    long nv;
    long ne;
    igraph_t graph;
    
public:
    Graph(long nv);      // empty constructor
    Graph(const Graph & obj); // copy const constructor
    Graph(Graph & obj); // copy constructor
    ~Graph(); // destructor function
    Graph(UncertainGraph & obj ); 
    void set_edges(igraph_vector_t * edges); // add edges to graph
    
    // graph property
    void graphStatstic(void);
    void metric(igraph_vector_t *res, double p); // used for uncertain graph for cal AD
    
    void reliablity(igraph_vector_t * res,double p); // used for uncertain/certain for cal reliablity
    void reliablity(igraph_vector_t * res, string filePath); // used for storeing reliablity result into file
    
    void entropyReport(igraph_vector_t *res, igraph_real_t maxDegree);  // used only for certain graph, maxDegre=
    double testAgaist(igraph_vector_t * ak, int k); // test against adversary knowledge
    void testDebug(igraph_vector_t * ak, int k);  // just for debuging the entropy function
    void selfTest(int k); // test against adversary knowledge=self.degrees
    
    // obfucation functions 
    void sigmaUniquess(igraph_vector_t * uv, igraph_vector_t ak, igraph_real_t maxDegree, igraph_real_t sigma);
    UncertainGraph generateObfuscation(igraph_real_t sigma, igraph_real_t * eps_res);
    UncertainGraph obfuscation();
    
    
    // obfucation against adversary knowledge
    UncertainGraph generateObfuscation(igraph_real_t sigma, igraph_real_t * eps_res, igraph_vector_t * ak);
    UncertainGraph obfuscation(igraph_vector_t *ak);
    
    
    
    // graph I/O
    void print_graph(string filepath);
};

// Graph I/O
Graph init_from_file(string filepath);
Graph init_from_Adj_File(string filepath);
//vector file I.O
void write_vector_file(igraph_vector_t *res, string filePath);
void init_vector_file(igraph_vector_t *res, string filePath);



#endif /* Graph_hpp */
