//
//  main.cpp
//  testGraphCplus
//
//  Created by dongqingxiao on 3/19/16.
//  Copyright © 2016 dongqingxiao. All rights reserved.
//

#include <iostream>
#include <Graph.hpp>
#include <UncertainGraph.hpp>
#include <map>
#include <DDCal.hpp>
using namespace std;
//using boost::math::normal;

int k;
int c;
int attempt;
double epsilon;
double noise;
char delimiter;
int sampleNum;

void graphCastTest(){
    delimiter='\t';
    string filepath="/Users/dongqingxiao/pythonEx/probGraph/exampleProbGraph.txt";
    UncertainGraph graph=init_uncertain_from_file(filepath);
    Graph g(graph);
    
    g.graphStatstic();
    
}

void graphTest(){
    // test baisc function about graph class
    k=100;
    c=2;
    attempt=1;
    epsilon=0.0001;
    noise=0.01;
    delimiter=' ';
    string filepath="/Users/dongqingxiao/Documents/uncetainGraphProject/graphs/dblp.txt";
    
    Graph g=init_from_file(filepath);
    
    g.graphStatstic();
    g.selfTest(100);
    igraph_real_t eps_res;
    UncertainGraph tpg=g.generateObfuscation(0.01, &eps_res);
}


void reliablityComparision(){
    
    string relfoler="/Users/dongqingxiao/Documents/uncetainGraphProject/allDataSet/relOutput/dblp/";
    string dataset="DBLP";
    string input_suffix="_InPG_rel.txt";
    string rep_suffix="rep_rel.txt";
    string ob_suffix="outpg_rel.txt";
    string method="GREEDY1-1";
    long int nv=824774;
    
    igraph_vector_t inRel,repRel,outRel;
    igraph_vector_init(&inRel,nv);
    igraph_vector_init(&repRel,nv);
    igraph_vector_init(&outRel,nv);
    
    
    init_vector_file(&inRel, relfoler+dataset+input_suffix);
    init_vector_file(&repRel, relfoler+method+"_"+rep_suffix);
    init_vector_file(&outRel, relfoler+method+"_"+ob_suffix);
    
    cout<<"relative error"<<endl;
    cout<<"inRel vs outRel : "<<cal_relative_error_vector(&inRel,&outRel)<<endl;
    
    cout<<"mean error"<<endl;
    cout<<"inRel vs outRel : "<<cal_mean_error_vector(&inRel, &outRel)<<endl;
    
    cout<<"mean error"<<endl;
    cout<<"inRel vs repRel : "<<cal_mean_error_vector(&inRel, &repRel)<<endl;
    
    cout<<"mean error"<<endl;
    cout<<"repRel vs outRel : "<<cal_mean_error_vector(&repRel, &outRel)<<endl;
    

}

int main(){
    graphTest();
    //graphCastTest();
    //reliablityComparision();
    return 0;
}