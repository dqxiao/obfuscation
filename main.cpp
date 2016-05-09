//
//  main.cpp
//  testGraphCplus
//
//  Created by dongqingxiao on 3/19/16.
//  Copyright © 2016 dongqingxiao. All rights reserved.
//

#include <iostream>
#include "Graph.hpp"
#include "UncertainGraph.hpp"
#include <map>
#include "DDCal.hpp"
#include "Help.hpp"
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
FeatureCombineOption foption;



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
    
    
    
    string refdata="inpg_rel";
    long int nv=824774;
    
    vector<string> datasets;
    
//    datasets.push_back("rand_dblp_ob_c1.100000_k100_rel");
//    datasets.push_back("rand_dblp_ob_c1.100000_k60_rel");
//    datasets.push_back("rand_dblp_ob_c1.300000_k200_rel");

    
//    datasets.push_back("GREEDY5-1_dblpc3.000000_k300_ob_rel");
//    datasets.push_back("GREEDY1-1_dblpc3.000000_k200_ob_rel");
//    datasets.push_back("GREEDY5-1_dblpc2.000000_k100_ob_rel");
    
     // datasets.push_back("greedy_dblp_ob_c1.500000_k300_rel");
 //   datasets.push_back("rand_dblp_ob_c1.500000_k300_rel");
   // datasets.push_back("greedy_dblp_ob_c1.500000_k300_sigma1000.000000_rel");
  //  datasets.push_back("rand_dblp_ob_c1.700000_k300_rel");
    
    datasets.push_back("greedy_dblp_ob_c1.500000_k300_sigma1.000000EE_rel");
    
    igraph_vector_t ref;
    igraph_vector_init(&ref,nv);
    
    
    init_vector_file(&ref,relfoler+refdata+".txt");
    
    
    for(string dataset:datasets){
        
        igraph_vector_t test;
        
        igraph_vector_init(&test,nv);
        
        init_vector_file(&test, relfoler+dataset+".txt");
        
        cout<<"Original vs" <<dataset<<endl;
        
        cout<<"Mean Error"<<endl;
        
        cout<<cal_mean_error_vector(&ref,&test)<<endl;
        
        cout<<"Mean Relative Error"<<endl;
        
        cout<<cal_relative_error_vector(&ref, &test)<<endl;
        
        
        igraph_vector_destroy(&test);
    }
    
    igraph_vector_destroy(&ref);
    
    

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
//    option=greedPert;
//    foption=mutiply;
    
    c=1.5;
    attempt=5;
    epsilon=0.0001;
    noise=0.01;
    k=300;
    
    string filepath="/Users/dongqingxiao/Documents/uncetainGraphProject/allDataSet/input/dblp.txt";
    string obFilePath="/Users/dongqingxiao/Documents/uncetainGraphProject/allDataSet/obOutput/rand_dblp_ob";
    
    string ssufix="_c"+to_string(c)+"_k"+to_string(k)+".txt" ;// setting suffix ;
    
    
    UncertainGraph ug=init_uncertain_from_file(filepath);
    
    
    
    long int nv=ug.nv;
    
    igraph_vector_t ak;
    igraph_real_t eps_res,sigma;
    igraph_vector_init(&ak,nv);
    ug.getDegrees(true, &ak);
    
  
    sigma=1;
    
    
    UncertainGraph tpg=ug.generateObfuscation(sigma, &eps_res, &ak);
    
    
    
    tpg.print_graph(obFilePath+ssufix);
    
   

    igraph_vector_destroy(&ak);
    

}

void greedyPerturbationTest(){
    delimiter='\t';
    //option=randPert; // set the random option
    option=greedPert;
    foption=mutiply;
    
    c=1.5;
    attempt=2;
    epsilon=0.0001;
    noise=0.01;
    k=300;
    
    string filepath="/Users/dongqingxiao/Documents/uncetainGraphProject/allDataSet/input/dblp.txt";
    string obFilePath="/Users/dongqingxiao/Documents/uncetainGraphProject/allDataSet/obOutput/greedy_dblp_ob";
    igraph_real_t eps_res,sigma;
    
    
    string utiltySetting="EE_EX"; // existing edges p(e)
    
    string ssufix="_c"+to_string(c)+"_k"+to_string(k)+"search"+utiltySetting+".txt" ;// setting suffix ;
    
    
    UncertainGraph ug=init_uncertain_from_file(filepath);
    
    
    
    long int nv=ug.nv;
    
    igraph_vector_t ak;
    
    igraph_vector_init(&ak,nv);
    ug.getDegrees(true, &ak);
    // ug.testAgaist(&ak);
    
    
    
   // sigma=0.0625;
    //UncertainGraph tpg=ug.generateObfuscation(sigma, &eps_res, &ak);
   UncertainGraph tpg=ug.obfuscation(&ak);
    
    
    
    tpg.print_graph(obFilePath+ssufix);
    
    
    
    igraph_vector_destroy(&ak);

}

void generateReliablity(){
    
    delimiter='\t';
    string obfolder="/Users/dongqingxiao/Documents/uncetainGraphProject/allDataSet/obOutput/";
    string relfolder="/Users/dongqingxiao/Documents/uncetainGraphProject/allDataSet/relOutput/dblp/";
    
    long int nv=824774;
    vector<string> datasets;
    
//    datasets.push_back("GREEDY5-1_dblpc3.000000_k300_ob");
//    datasets.push_back("rand_dblp_ob_c1.100000_k60");
//    datasets.push_back("rand_dblp_ob_c1.300000_k200");
//    datasets.push_back("rand_dblp_ob_c1.500000_k300");
//    datasets.push_back("rand_dblp_ob_c1.100000_k100_rel");
   // datasets.push_back("greedy_dblp_ob_c1.500000_k300");
   // datasets.push_back("greedy_dblp_ob_c1.500000_k300_sigma1000.000000");
    //datasets.push_back("rand_dblp_ob_c1.700000_k300");
  //  datasets.push_back("greedy_dblp_ob_c1.500000_k300_sigma1.000000EE");
    
   // datasets.push_back("greedy_dblp_ob_c1.500000_k300searchEE");
   // datasets.push_back("greedy_dblp_ob_c1.500000_k200searchEE");
    
    //datasets.push_back("greedy_dblp_ob_c1.500000_k100searchEE");
    datasets.push_back("greedy_dblp_ob_c1.100000_k60searchEE");
    
    
    
    for(string dataset: datasets){
        UncertainGraph uobg=init_uncertain_OB_from_file(obfolder+dataset+".txt", nv);
        igraph_vector_t rv;
        igraph_vector_init(&rv, nv);
        sampleNum=2000;
        uobg.reliablity(&rv, relfolder+dataset+"_rel.txt");
        igraph_vector_destroy(&rv);
    }

}

void generateInReliablity(){
    
    delimiter='\t';
    string folder="/Users/dongqingxiao/Documents/uncetainGraphProject/allDataSet/obOutput/";
//    string dataset="GREEDY5-1_dblpc3.000000_k300_ob";
    string relFolder="/Users/dongqingxiao/Documents/uncetainGraphProject/allDataSet/relInOutput/dblp/";
    
    vector<string> datasets;
    
    datasets.push_back("");
    datasets.push_back("");
    for(string dataset: datasets){
        
        long int nv=824774;
        UncertainGraph uobg=init_uncertain_OB_from_file(folder+dataset+".txt", nv);
        
        long int cSampleNum=2000;
        uobg.reliablity_record(cSampleNum, relFolder+dataset+"_s2000_rel.txt");
    }
    
}

void basic_metric(){
    delimiter='\t';
    string filepath="/Users/dongqingxiao/Documents/uncetainGraphProject/allDataSet/input/";
    string dataset="dblp";
    
    UncertainGraph ug=init_uncertain_from_file(filepath+dataset+".txt");
    
    ug.graphStastic();
    
}


void exact_reliablityComparision_test(){
    int sampleNum=3;
    long int nv=3;
    boost::numeric::ublas::matrix<int> mref(nv,sampleNum);
    boost::numeric::ublas::matrix<int> mTest(nv,sampleNum);
    
    
    string folder="/Users/dongqingxiao/Documents/uncetainGraphProject/allDataSet/progTest/";
    
    string inDatast="relRef.txt";
    string randDatast="refTest.txt";
   
    
    init_boostMatrix_file(mref, folder+inDatast); // almost need half hour
    init_boostMatrix_file(mTest, folder+randDatast);
    
    cout<<"diff:"<<compare_In_Matrix(mref,mTest)<<endl;
    
    mref.clear();
    mTest.clear();
    
    
    
}

void exact_reliablityComparision(){
    int sampleNum=2000;
    long int nv=824774;
    boost::numeric::ublas::matrix<int> mref(nv,sampleNum);
    
    
    boost::numeric::ublas::matrix<int> mTest(nv,sampleNum);
        
    
    string folder="/Users/dongqingxiao/Documents/uncetainGraphProject/allDataSet/relInOutput/dblp/";
    
    string inDatast="dblp_s2000.txt";
    string ranDatast="rand_dblp_ob_c1.300000_k200_rel.txt";
    
    init_boostMatrix_file(mref, folder+inDatast); // almost need half hour
    init_boostMatrix_file(mTest, folder+ranDatast);
    
    
    
    
    cout<<"compare "<<inDatast<<" vs "<<ranDatast<<endl;
    double diff=sampleing_compare_In_Matrix(mref,mTest);
    cout<<"diff "<<diff<<endl;
    
    mref.clear();
    mTest.clear();
    
    
    
}

void timeestimation(){
    for(long int i=0;i<8247;i++){
        for(long int j=0;j<8247;j++){
            for(int k=0;k<2000;k++){
            }
            
        }
        if(i%1000==0){
            cout<<i<<endl;
        }
    }
}


void certainObfuscation(){
    string folder="/Users/dongqingxiao/Documents/uncetainGraphProject/allDataSet/output/";
    string file="GREEDY5-1_dblp";
    
    string obFolder="/Users/dongqingxiao/Documents/uncetainGraphProject/allDataSet/obOutput/";
    
    c=3;
    attempt=1;
    epsilon=0.0001;
    noise=0.01;
    k=300;
    delimiter=' ';
    
    Graph g=init_from_Adj_File(folder+file+".txt");
    
    string suffix="c"+to_string(c)+"_k"+to_string(k)+"_ob.txt";
    
    
    string filepath="/Users/dongqingxiao/Documents/uncetainGraphProject/allDataSet/input/dblp.txt";
    delimiter='\t';
    UncertainGraph ug=init_uncertain_from_file(filepath);
    
    igraph_vector_t ak;
    igraph_vector_init(&ak,ug.nv);
    ug.getDegrees(false, &ak);
    
    

    
   // g.selfTest(k);
    double sigma=32;
    double esp_res=1;
    //UncertainGraph tpg=g.generateObfuscation(sigma, &esp_res);
    UncertainGraph tpg=g.generateObfuscation(sigma, &esp_res, &ak);
    
    tpg.print_graph(obFolder+file+suffix);
    
}



void randomPerturbationTest_traffic(){
    
    delimiter='\t';
    option=randPert; // set the random option
    c=1.1;
    attempt=1;
    epsilon=0.0001;
    noise=0.01;
    k=10;
    
    string filepath="/Users/dongqingxiao/Documents/uncetainGraphProject/allDataSet/input/trafficUK.txt";
    string obFilePath="/Users/dongqingxiao/Documents/uncetainGraphProject/allDataSet/obOutput/rand_trafficUK_ob";
    
    string ssufix="_c"+to_string(c)+"_k"+to_string(k)+".txt" ;// setting suffix ;
    
    
    UncertainGraph ug=init_uncertain_from_file(filepath);
    
    
    
    long int nv=ug.nv;
    
    igraph_vector_t ak;
    igraph_real_t eps_res,sigma;
    igraph_vector_init(&ak,nv);
    ug.getDegrees(true, &ak);
    ug.testAgaist(&ak);
    
//    sigma=;
//    
//    
//    UncertainGraph tpg=ug.generateObfuscation(sigma, &eps_res, &ak);
//    
//    
//    
//    tpg.print_graph(obFilePath+ssufix);
    
    
    
    igraph_vector_destroy(&ak);
    
    
}


void nullEdgeReliablity(){
    delimiter='\t';
    
    string file="/Users/dongqingxiao/Documents/uncetainGraphProject/allDataSet/input/dblp.txt";
    
    string relOuput="/Users/dongqingxiao/Documents/uncetainGraphProject/allDataSet/relOutput/dblp/reliablityNode/nodeReNullEdege.txt";
    string relEdgeOutput="/Users/dongqingxiao/Documents/uncetainGraphProject/allDataSet/relOutput/dblp/reliablityNode/edgeDiffReliablity.txt";
    
//
//    sampleNum=2000;
    UncertainGraph ug=init_uncertain_from_file(file);
    igraph_vector_t ruv_n,ruv_e,ruv;
    
    igraph_vector_init(&ruv_n,ug.nv);
    igraph_vector_init(&ruv_e,ug.nv);
    igraph_vector_init(&ruv,ug.nv);
    igraph_vector_t rue;
    igraph_vector_init(&rue,ug.ne);
    
//
//    ug.reliablityUtiltyDiff(&ruv);
//
//    
//    //print_vector(&ruv, "reliablity by edges whose p(e)=0");
//    
//    write_vector_file(&ruv, relOuput);
//    

    
    init_vector_file(&ruv_n, relOuput);
    init_vector_file(&rue,relEdgeOutput);
    
    cout<<"reliablity utilty of each nodes from no existent one "<<endl;
    vector_statstic(&ruv_n);
    
    cout<<"reliablity utilty of each existing ones"<<endl;
    vector_statstic(&rue);
    ug.aggregateReliablutyDiffEE(&ruv_e, &rue);
    
   // igraph_vector_add(&ruv, &ruv_n);
    igraph_vector_add(&ruv, &ruv_e);
    cout<<"reliablity utilty of all"<<endl;
    vector_statstic(&ruv);
    
    
    write_vector_file(&ruv, "/Users/dongqingxiao/Documents/uncetainGraphProject/allDataSet/relOutput/dblp/reliablityNode/nodeRelDiff_E.txt");
    
    igraph_vector_destroy(&ruv_n);
    igraph_vector_destroy(&ruv_e);
    igraph_vector_destroy(&rue);
    
    
    
}


void testAgaist(){
    
    string refData="/Users/dongqingxiao/Documents/uncetainGraphProject/allDataSet/input/dblp.txt";
    string testFolder="/Users/dongqingxiao/Documents/uncetainGraphProject/allDataSet/obOutput/";
    
    vector<string> datasets;
    
    datasets.push_back("rand_dblp_ob_c1.500000_k300_sigma_1");
    datasets.push_back("greedy_dblp_ob_c1.500000_k300_sigma1.000000");
    delimiter='\t';
    k=300;
    
    UncertainGraph oug=init_uncertain_from_file(refData);
    igraph_vector_t ak;
    igraph_vector_init(&ak,oug.nv);
    
    oug.getDegrees(true, &ak);
    
    
    for(string datatset: datasets){
        UncertainGraph ug=init_uncertain_OB_from_file(testFolder+datatset+".txt", oug.nv);
        ug.testAgaist(&ak);
    }
    
    igraph_vector_destroy(&ak);
    
}

// change into script version later
int main(int argc, char *argv[]){
    
    
    
    
    

   // graphTest();
   // graphCastTest();
   // reliablityComparision();
 //   testObfuscation();
//    reliablityUtiltyTest();
    randomPerturbationTest();
   // greedyPerturbationTest();
  //  randomPerturbationTest_traffic();
    
    
    
    
   // generateReliablity();
   // generateInReliablity();
    
    
   // basic_metric();
    //exact_reliablityComparision();
    //timeestimation();
    //certainObfuscation();
    
    //nullEdgeReliablity();
   // testAgaist();
    return 0;
}