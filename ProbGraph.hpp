////
////  ProbGraph.hpp
////  testGraphCplus
////
////  Created by dongqingxiao on 3/19/16.
////  Copyright © 2016 dongqingxiao. All rights reserved.
////
//
//#ifndef ProbGraph_hpp
//#define ProbGraph_hpp
//
//#include <stdio.h>
//#include <igraph.h>
//#include <string>
//#include <fstream>
//#include <sstream>
//#include <vector>
//#include <boost/pending/disjoint_sets.hpp>
//
//using namespace std;
//
//class ProbGraph{
//private:
//    long int nv;
//    long int ne;
//    igraph_t graph;
//    void uncertain_entropyReport(igraph_vector_t * entropyReport,igraph_real_t maxDegree);
//    void certain_entropyReport(igraph_vector_t * entropyReport);
//    void ex_uncertain_entropyReport(igraph_vector_t * entropyReport,igraph_real_t maxDegree);
//    void sigmaUniquess(igraph_vector_t * uv, igraph_vector_t ak, igraph_real_t maxDegree, igraph_real_t sigma);
//    void prob_vertex(igraph_vector_t * prob_v,igraph_integer_t vid);
//    void print_uncertain_graph(string filepath);
//    void print_certain_graph(string filepath);
//public:
////    ProbGraph(const ProGraph &obj);  // copy constructor 
//    ProbGraph(igraph_integer_t n);  // empty constructor
//    ProbGraph(const ProbGraph & obj); // copy constuctor 
//    ~ProbGraph(){
//        printf("die\n");
//        igraph_destroy(&graph);
//    };
//    
//    long int getNV(void);
//    void maxDegrees(igraph_vector_t *deg);
//    igraph_real_t maxDegree(void);
//    void set_edges(const igraph_vector_t *edges);
//    void set_edges_prob(const igraph_vector_t *probs);
//    
//    void entropyReport(igraph_vector_t * entropyReport,igraph_real_t maxDegree);
//    igraph_real_t testAgainst(igraph_vector_t * ak);
//    igraph_real_t selfTest(void);
//    ProbGraph obfuscation(void);
//    void certainGraphStatstic(void);
//    void uncertainGraphStastic(void);
//    void certain_metrics(igraph_vector_t * res, igraph_real_t p);
//    
//    // suitable for small dataset later switch back to matrix
//    void uncertain_reliablityReport(igraph_spmatrix_t *res);
//    void certain_reliablityReport(igraph_spmatrix_t *res, igraph_real_t p);
//    void reliablityReport(void);
//
//    // suitable for medium dataset
//    
//    void certain_reliablityReport_file(string filePath, igraph_real_t p);
//    void uncertain_reliablityReport_file(string filePath);
//    
//    // suitable for large dataset estimation
//    void certain_reliablityReport_file_ex(string filePath, igraph_real_t p);
//    void uncertain_reliablityReport_file_ex(string filePath);
//    
//    // bounded reliablity report esitmation
//    void certain_reliablityReport_ex(igraph_vector_t *res,igraph_real_t p);
//    void certain_reliablityReport_ex_inc(igraph_vector_t *res,igraph_real_t p);
//    
//    void uncertain_reliablityReport_ex(igraph_vector_t *res);
//    void uncertain_reliablityReport_ex_store(igraph_vector_t *res, string filepath, string label);
//    void certain_reliablityReport_ex_store(igraph_vector_t *res, string filepath, string label);
//    
//    ProbGraph generateObfuscation(igraph_real_t sigma, igraph_real_t * eps_res);
//    ProbGraph sampleGraph(void);
//    void print_graph(string filepath); // all information will be store in graph
//    
//
//    
//    // just for debugging
//    
//    
//};
//
//
//ProbGraph certain_init_from_File(string path);
//ProbGraph uncertain_init_from_File(string path);
//ProbGraph init_from_File(string path);      //
//ProbGraph init_from_Adj_File(string path); // used for certain graph  only
//ProbGraph init_from_OB_File(string path,long int nv); // used for uncertain graph <-obfuscated one
//
//
//void init_vector_from_file(string path, igraph_vector_t *res);
//
//
////extern int k;
////extern int c;
////extern int attempt;
////extern double epsilon;
////extern double noise;
////extern char delimiter;
////extern bool uncertain;
////extern int sampleNum; // the number of sample  
//
//#endif /* ProbGraph_hpp */
