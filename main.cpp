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
double c;
int attempt;
double epsilon;
double noise;
char delimiter;
int sampleNum;
Option option;



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
    ug.getDegrees(false, &ak);
    UncertainGraph uog=init_uncertain_OB_from_file(obPath+method, nv);
    uog.testAgaist(&ak);
    //ug.testAgaist(&ak);
    
}

void reliablityComparision(){
    
    string relfoler="/Users/dongqingxiao/Documents/uncetainGraphProject/allDataSet/relOutput/dblp/";
    string dataset="DBLP";
    string input_suffix="_InPG_rel.txt";
    string rep_suffix="rep_rel.txt";
    string ob_suffix="outpg_rel.txt";
    string method="GREEDY1-1";
    string ourMethod="rand_dblp_ob_c1.1";
    
    long int nv=824774;
    
    igraph_vector_t inRel,repRel,outRel,randRel;
    igraph_vector_init(&inRel,nv);
    igraph_vector_init(&repRel,nv);
    igraph_vector_init(&outRel,nv);
    igraph_vector_init(&randRel,nv);
    
    
    
    init_vector_file(&inRel, relfoler+dataset+input_suffix);
    init_vector_file(&repRel, relfoler+method+"_"+"s2000_"+rep_suffix);
    init_vector_file(&outRel, relfoler+method+"_"+ob_suffix);
    init_vector_file(&randRel,relfoler+ourMethod+"_"+"rel.txt");
    
    
    
//    cout<<"mean error"<<endl;
//    cout<<"inRel vs outRel : "<<cal_mean_error_vector(&inRel, &outRel)<<endl;
//    
//    cout<<"mean error"<<endl;
//    cout<<"inRel vs repRel : "<<cal_mean_error_vector(&inRel, &repRel)<<endl;
//    
//    cout<<"mean error"<<endl;
//    cout<<"repRel vs outRel : "<<cal_mean_error_vector(&repRel, &outRel)<<endl;
//    
//    
//    cout<<"mean error"<<endl;
//    
//    cout<<"inRel vs randRel : "<<cal_mean_error_vector(&inRel, &randRel)<<endl;
    
    cout<<"sparse"<<endl;
    
    is_sparese(&inRel, nv);
 
    
    
    cout<<"mean relative Error"<<endl;
    
    cout<<"inRel vs randRel : "<<cal_relative_error_vector(&inRel, &randRel)<<endl;
    
    
    cout<<"mean error"<<endl;
    cout<<"inRel vs outRel : "<<cal_relative_error_vector(&inRel, &outRel)<<endl;
    
    cout<<"mean error"<<endl;
    cout<<"inRel vs repRel : "<<cal_relative_error_vector(&inRel, &repRel)<<endl;
    
    cout<<"mean error"<<endl;
    cout<<"repRel vs outRel : "<<cal_relative_error_vector(&repRel, &outRel)<<endl;
    
    
    cout<<"mean error"<<endl;
    
    cout<<"inRel vs randRel : "<<cal_mean_error_vector(&inRel, &randRel)<<endl;

}

void reliablityUtiltyTest(){
    // uncertain graph
    delimiter='\t';
    sampleNum=10000;
    string filepath="/Users/dongqingxiao/Documents/uncetainGraphProject/allDataSet/input/dblp.txt";
    string testFilePath="/Users/dongqingxiao/Documents/uncetainGraphProject/allDataSet/progTest/testUncertainGraph.txt";
    string ruvPath="/Users/dongqingxiao/Documents/uncetainGraphProject/allDataSet/progTest/ruv_dblp.txt";
    UncertainGraph ug=init_uncertain_from_file(filepath);
    
   // Graph g(ug);
    

    
    

    
    long int nv=ug.nv;
//
//
    
    igraph_vector_t ruv,rel;
    igraph_vector_init(&ruv,nv);
//
    igraph_vector_init(&rel,nv);
    
    
    ug.reliablityUtilitySubgraph(&ruv);
//    ug.reliablityUtiliy(&ruv);
////
////    
////
    vector_statstic(&ruv);
//
    write_vector_file(&ruv, ruvPath);
//
    igraph_vector_destroy(&ruv);
    
    
    
    
    
    // inverstigate degree
    
//    igraph_vector_t degres;
//    igraph_vector_init(&degres,nv);
//    
//    g.degrees(&degres);
//    string degfile="/Users/dongqingxiao/Documents/uncetainGraphProject/allDataSet/progTest/deg_dblp.txt";
//    vector_statstic(&degres);
//    write_vector_file(&degres,degfile);
//    
//    igraph_vector_destroy(&degres);
    
    
    
    
    
    
    
}

void randomPerturbationTest(){
    
    delimiter='\t';
    option=randPert; // set the random option
    c=1.3;
    attempt=1;
    epsilon=0.0001;
    noise=0.01;
    k=200;
    
    string filepath="/Users/dongqingxiao/Documents/uncetainGraphProject/allDataSet/input/dblp.txt";
    string testFilePath="/Users/dongqingxiao/Documents/uncetainGraphProject/allDataSet/progTest/testUncertainGraph.txt";
    string ruvPath="/Users/dongqingxiao/Documents/uncetainGraphProject/allDataSet/progTest/ruv_dblp.txt";
    string obFilePath="/Users/dongqingxiao/Documents/uncetainGraphProject/allDataSet/obOutput/rand_dblp_ob";
    
    string ssufix="_c"+to_string(c)+"_k"+to_string(k)+".txt" ;// setting suffix ;
    
    
    UncertainGraph ug=init_uncertain_from_file(filepath);
    
    
    
    long int nv=ug.nv;
    
    igraph_vector_t ak;
    igraph_real_t eps_res,sigma;
    igraph_vector_init(&ak,nv);
    ug.getDegrees(true, &ak);
    
    sigma=1;
//
//    
    UncertainGraph tpg=ug.generateObfuscation(sigma, &eps_res, &ak);
    
  // UncertainGraph tpg=ug.obfuscation(&ak);
    
   tpg.print_graph(obFilePath+ssufix);
    
//    ug.testAgaist(&ak);

    igraph_vector_destroy(&ak);
    

}

void generateReliablity(){
    
    delimiter='\t';
    string filepath="/Users/dongqingxiao/Documents/uncetainGraphProject/allDataSet/obOutput/rand_dblp_ob_c1.300000_k200.txt";
    string relPath="/Users/dongqingxiao/Documents/uncetainGraphProject/allDataSet/relOutput/dblp/rand_dblp_ob_c1.300000_k200_rel.txt";
    UncertainGraph uobg=init_uncertain_from_file(filepath);
    
    long int nv=uobg.nv;
    
    igraph_vector_t rv;
    igraph_vector_init(&rv, nv);
    sampleNum=2000;
    uobg.reliablity(&rv, relPath);
    
    igraph_vector_destroy(&rv);

}

void generateInReliablity(){
    
    delimiter='\t';
    string filepath="/Users/dongqingxiao/Documents/uncetainGraphProject/allDataSet/input/dblp.txt";
    string relPath="/Users/dongqingxiao/Documents/uncetainGraphProject/allDataSet/relInOutput/dblp/dblp_s2000.txt";
    UncertainGraph uobg=init_uncertain_from_file(filepath);
    
    long int cSampleNum=2000;
    uobg.reliablity_record(cSampleNum, relPath);
    
    
}

void basic_metric(){
    delimiter='\t';
    string filepath="/Users/dongqingxiao/Documents/uncetainGraphProject/allDataSet/input/";
    string dataset="dblp";
    
    UncertainGraph ug=init_uncertain_from_file(filepath+dataset+".txt");
    
    ug.graphStastic();
    
}

//
void exact_reliablityComparision(){
    long int nv=824774;
    
    for(long int i=0;i<nv;i++){
        for(long int j=i+1;j<nv;j++){
            for(int k=0;k<2;k++){
                //
            }
        }
        
        if(i%1000==0){
            cout<<"move"<<i<<endl;
        }
    }
}


int main(){
   // graphTest();
   // graphCastTest();
  //  reliablityComparision();
 //   testObfuscation();
//    reliablityUtiltyTest();
   // randomPerturbationTest();
    
    //generateReliablity();
    //generateInReliablity();
    
   // basic_metric();
   // exact_reliablityComparision();
    return 0;
}