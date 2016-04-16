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
#include <Help.hpp>
#include <vector> 
#include <cmath>
#include <boost/tuple/tuple.hpp>

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
    UncertainGraph ug=init_uncertain_from_file(filepath);
    
    Graph g(ug);
    g.graphStatstic();
    ug.graphStastic();
    igraph_vector_t degs;
    igraph_vector_init(&degs,4);
    //ug.getDegrees(true, &degs);
    //print_vector(&degs, "expected degree");
    
    
    
}

void graphTest(){
    // test baisc function about certain class
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


void testObfuscation(){
    k=100;
    c=2;
    attempt=1;
    epsilon=0.0001;
    noise=0.01;
    delimiter='\t';
    igraph_vector_t ak;
    long int nv;
    string inputPath="/Users/dongqingxiao/Documents/uncetainGraphProject/allDataSet/input/";
    string obPath="/Users/dongqingxiao/Documents/uncetainGraphProject/allDataSet/obOutput/";
    string dataset="dblp";
    string method="ADR1_dblp_ob.txt";
    
    UncertainGraph ug=init_uncertain_from_file(inputPath+dataset+".txt");
    nv=ug.nv;
    igraph_vector_init(&ak,nv);
    ug.getDegrees(true, &ak);
    UncertainGraph uog=init_uncertain_OB_from_file(obPath+method, nv);
//    uog.testAgaist(&ak);
    ug.testAgaist(&ak);
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


void reliablityUtiltyTest(){
    // uncertain graph
    delimiter='\t';
    sampleNum=20;
    string filepath="/Users/dongqingxiao/Documents/uncetainGraphProject/allDataSet/input/dblp.txt";
    string testFilePath="/Users/dongqingxiao/Documents/uncetainGraphProject/allDataSet/progTest/testUncertainGraph.txt";
    string ruvPath="/Users/dongqingxiao/Documents/uncetainGraphProject/allDataSet/progTest/ruv_dblp.txt";
    UncertainGraph ug=init_uncertain_from_file(filepath);
    
    long int nv=ug.nv;
    
    igraph_vector_t ruv,rel;
    igraph_vector_init(&ruv,nv);
    
    igraph_vector_init(&rel,nv);
    

    ug.reliablityUtiliy(&ruv);

    vector_statstic(&ruv);
    
    write_vector_file(&ruv, ruvPath);
    
    igraph_vector_destroy(&ruv);
    
    
}



int main(){
   // graphTest();
   // graphCastTest();
    //reliablityComparision();
    //testObfuscation();
    reliablityUtiltyTest();
    return 0;
}