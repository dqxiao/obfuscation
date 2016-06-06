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
double fplusParmater;
bool debug;
WorkDataset w_dataset;



void graphTest(){
    // test baisc function about certain class
    k=100;
    c=2;
    attempt=5;
    epsilon=0.0001;
    noise=0.01;
    delimiter=' ';
    string filepath="/Users/dongqingxiao/Documents/uncetainGraphProject/graphs/dataset/dblp.txt";
    
    Graph g=init_from_file(filepath);
    
    g.graphStatstic();
    g.selfTest(100);
    igraph_real_t eps_res;
    UncertainGraph tpg=g.generateObfuscation(0.01, &eps_res);
    
    tpg.graphStastic();
    
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

void reliablityComparision(string dataset){
    
    string relfoler="/Users/dongqingxiao/Documents/uncetainGraphProject/allDataSet/relOutput/"+dataset+"/";
    
    
    long int nv=0;
    
    if(dataset=="hepth"){
        nv=14318;
    }
    if(dataset=="dblp"){
        nv=824774;
    }
    
    string refdata="ref_rel";
   
    
    vector<string> datasets;
    
    if(dataset=="hepth"){
//        datasets.push_back("rand_ob_c1.100000_k60_sigma0.000977");
//        datasets.push_back("rand_ob_c1.300000_k100_sigma0.250000");
//        datasets.push_back("rand_ob_c1.700000_k200_sigma17.125977");
//        datasets.push_back("rand_ob_c3.000000_k300_sigma1.875977");
//        
//        
//        datasets.push_back("greedy_ob_c1.100000_k60_sigma0.009766");
//        datasets.push_back("greedy_ob_c1.300000_k100_sigma0.091797");
//        datasets.push_back("greedy_ob_c2.000000_k200_sigma0.278320");
//        datasets.push_back("greedy_ob_c3.000000_k300_sigma0.419922");
//       
//        datasets.push_back("GREEDY5-1_hepthc2.000000_k100sigma0.031250_ob");
//        datasets.push_back("GREEDY5-1_hepthc2.000000_k60sigma0.001953_ob");
//        datasets.push_back("GREEDY5-1_hepthc3.000000_k200sigma0.062500_ob");
//        datasets.push_back("GREEDY5-1_hepthc3.000000_k300sigma1.000000_ob");

    }
    
    if(dataset=="dblp"){
        datasets.push_back("GREEDY5-1_dblpc2.000000_k60_ob");
        datasets.push_back("GREEDY5-1_dblpc2.000000_k100_ob");
        datasets.push_back("GREEDY1-1_dblpc3.000000_k200_ob");
        datasets.push_back("GREEDY5-1_dblpc3.000000_k300_ob");
        
    }
    
    igraph_vector_t ref;
    igraph_vector_init(&ref,nv);
    
    
    init_vector_file(&ref,relfoler+refdata+".txt");
    
    
    for(string dataset:datasets){
        
        igraph_vector_t test;
        
        igraph_vector_init(&test,nv);
        
        init_vector_file(&test, relfoler+dataset+"_rel.txt");
        
        cout<<"Original vs " <<dataset<<endl;
        
        cout<<"Mean Error"<<endl;
        
        cout<<cal_mean_error_vector(&ref,&test)<<endl;
        
        cout<<"Mean Relative Error"<<endl;
        
        cout<<cal_relative_error_vector(&ref, &test)<<endl;
        
        
        igraph_vector_destroy(&test);
    }
    
    igraph_vector_destroy(&ref);
    
    

}

void reliablityUtilty(string dataset){
    // uncertain graph
    delimiter='\t';
    sampleNum=10000;
    string filepath="/Users/dongqingxiao/Documents/uncetainGraphProject/allDataSet/input/"+dataset+".txt";
    string ruvPath="/Users/dongqingxiao/Documents/uncetainGraphProject/allDataSet/relOutput/"+dataset+"/reliablityNode/Node_EE.txt";
    string euvPath="/Users/dongqingxiao/Documents/uncetainGraphProject/allDataSet/relOutput/"+dataset+"/reliablityNode/Edge_EE.txt";
    UncertainGraph ug=init_uncertain_from_file(filepath);
    
    

    
    

    
    long int nv=ug.nv;
    long int ne=ug.ne;
    
    
    igraph_vector_t r_edge,ruv;
    
    igraph_vector_init(&r_edge,ne);
    igraph_vector_init(&ruv,nv);
    
    
    ug.rawEstimate(&r_edge);
    
    ug.aggregateReliablutyDiffEE(&ruv, &r_edge);
    
    
    write_vector_file(&r_edge, euvPath);
    write_vector_file(&ruv, ruvPath);
    
    igraph_vector_destroy(&r_edge);
    igraph_vector_destroy(&ruv);
    
    
}

void randomPerturbation(string dataset){
    
    delimiter='\t';
    option=randPert;
    
    
    
    string filepath="/Users/dongqingxiao/Documents/uncetainGraphProject/allDataSet/input/"+dataset+".txt";
    
    string obFilePath="/Users/dongqingxiao/Documents/uncetainGraphProject/allDataSet/obOutput/"+dataset+"/rand_ob";
    

    attempt=5;
    epsilon=0.001; // different epsilon setting
    noise=0.01;

    
    
    vector<int> ks;
    vector<double> cs;
    
//    ks.push_back(60);
//    ks.push_back(100);
//    ks.push_back(200);
//    ks.push_back(300);
    
//    cs.push_back(1.1);
//    cs.push_back(1.3);
//    cs.push_back(1.7);
//    cs.push_back(3);
    
    
    UncertainGraph ug=init_uncertain_from_file(filepath);
    long int nv=ug.nv;
    
    igraph_vector_t ak;
    igraph_vector_init(&ak,nv);
    ug.getDegrees(true, &ak);
    
    // self test
    
    
    for(int i=0;i<ks.size();i++){
    
        k=ks[i];
        c=cs[i];
       // ug.testAgaist(&ak);
        double finalSigma=0;
        UncertainGraph tpg=ug.obfuscation(&ak,&finalSigma);
        string fprefix="BI";
        string ssufix="_c"+to_string(c)+"_k"+to_string(k)+"_sigma"+to_string(finalSigma)+".txt" ;// setting suffix ;
      
        tpg.print_graph(obFilePath+fprefix+ssufix);
    }
    
    igraph_vector_destroy(&ak);
    

}

void greedyPerturbation(string dataset){
    delimiter='\t';
    option=greedPert;
    foption=mutiply;
    
    
    
    string filepath="/Users/dongqingxiao/Documents/uncetainGraphProject/allDataSet/input/"+dataset+".txt";
    
    string obFilePath="/Users/dongqingxiao/Documents/uncetainGraphProject/allDataSet/obOutput/"+dataset+"/greedy_ob";
    
    
    attempt=5;
    epsilon=0.001; // different epsilon setting
    noise=0.01;
    
    
    
    vector<int> ks;
    vector<double> cs;
    
    if(dataset=="hepth"){
    
        ks.push_back(60);
        ks.push_back(100);
        ks.push_back(200);
        ks.push_back(300);
    //
        cs.push_back(1.1);
        cs.push_back(1.3);
        cs.push_back(2);
        cs.push_back(3);
    }
    
    
   
    
    UncertainGraph ug=init_uncertain_from_file(filepath);
    long int nv=ug.nv;
    
    igraph_vector_t ak;
    igraph_vector_init(&ak,nv);
    ug.getDegrees(true, &ak);
    
    // self test
    
    
    for(int i=0;i<ks.size();i++){
        
        k=ks[i];
        c=cs[i];
        // ug.testAgaist(&ak);
        double finalSigma=0;
        UncertainGraph tpg=ug.obfuscation(&ak, &finalSigma);
        string fsufix="ER"; // ER: existing
        string ssufix="_c"+to_string(c)+"_k"+to_string(k)+"_sigma"+to_string(finalSigma)+".txt" ;// setting suffix ;
        
        tpg.print_graph(obFilePath+fsufix+ssufix);
    }
    
    igraph_vector_destroy(&ak);

}

void generateReliablityRef(string dataset){
    delimiter='\t';
    string infolder="/Users/dongqingxiao/Documents/uncetainGraphProject/allDataSet/input/";
    string relfolder="/Users/dongqingxiao/Documents/uncetainGraphProject/allDataSet/relOutput/"+dataset+"/";
    
    long int nv=0;
    if(dataset=="hepth"){
        nv=14318;
    }
    
    UncertainGraph ug=init_uncertain_from_file(infolder+dataset+".txt");
    igraph_vector_t rv;
    igraph_vector_init(&rv,nv);
    
    sampleNum=2000;
    
    ug.reliablity(&rv, relfolder+"ref"+"_rel.txt");
    igraph_vector_destroy(&rv);
}

void generateReliablity(string dataset){
    
    delimiter='\t';
    string obfolder="/Users/dongqingxiao/Documents/uncetainGraphProject/allDataSet/obOutput/"+dataset+"/";
    string relfolder="/Users/dongqingxiao/Documents/uncetainGraphProject/allDataSet/relOutput/"+dataset+"/";
    
//    long int nv=824774;
    
    long int nv=0;
    
    if(dataset=="hepth"){
        nv=14318;
    }
    if(dataset=="dblp"){
        nv=824774;
    }
    vector<string> datasets;
    
    
//    datasets.push_back("rand_ob_c1.100000_k60_sigma0.000977");
//    datasets.push_back("rand_ob_c1.300000_k100_sigma0.250000");
//    datasets.push_back("rand_ob_c1.700000_k200_sigma17.125977");
//    datasets.push_back("rand_ob_c3.000000_k300_sigma1.875977");
//    
//    datasets.push_back("greedy_ob_c1.100000_k60_sigma0.009766");
//    datasets.push_back("greedy_ob_c1.300000_k100_sigma0.091797");
//    datasets.push_back("greedy_ob_c2.000000_k200_sigma0.278320");
//    datasets.push_back("greedy_ob_c3.000000_k300_sigma0.419922");
    
//    datasets.push_back("GREEDY5-1_hepthc2.000000_k100sigma0.031250_ob");
//    datasets.push_back("GREEDY5-1_hepthc2.000000_k60sigma0.001953_ob");
//    datasets.push_back("GREEDY5-1_hepthc3.000000_k200sigma0.062500_ob");
//    datasets.push_back("GREEDY5-1_hepthc3.000000_k300sigma1.000000_ob");
    
    
   // datasets.push_back("greedy_obER_c3.000000_k300_sigma1.749023");
    datasets.push_back("rand_obBI_c3.000000_k300_sigma2.170898");
    
    
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

void basic_metric(string dataset){
    delimiter='\t';
    string filepath="/Users/dongqingxiao/Documents/uncetainGraphProject/allDataSet/input/";
    //string dataset="hepth";
    
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


void certainObfuscation(string dataset,string repDataset){
    string folder="/Users/dongqingxiao/Documents/uncetainGraphProject/allDataSet/output/";

    
    string obFolder="/Users/dongqingxiao/Documents/uncetainGraphProject/allDataSet/obOutput/"+dataset+"/"+repDataset;
    
    
    attempt=2;
    epsilon=0.001;
    noise=0.01;
 

    vector<int> ks;
    vector<double> cs;
    
    ks.push_back(60);
    ks.push_back(100);
    ks.push_back(200);
    ks.push_back(300);
    
    cs.push_back(2);
    cs.push_back(2);
    cs.push_back(3);
    cs.push_back(3);
    
    delimiter=' ';
    Graph g=init_from_Adj_File(folder+repDataset+".txt");
    string filepath="/Users/dongqingxiao/Documents/uncetainGraphProject/allDataSet/input/"+dataset+".txt";
    delimiter='\t';
    UncertainGraph ug=init_uncertain_from_file(filepath);
    
    
    
    igraph_vector_t ak;
    igraph_vector_init(&ak,ug.nv);
    ug.getDegrees(true, &ak);
    
    
    for(int i=0;i<ks.size();i++){
        k=ks[i];
        c=cs[i];
        
//        g.testAgaist(&ak, k);
        
        double finalSigma=0;
        UncertainGraph tpg=g.obfuscation(&ak,&finalSigma);
        string suffix="c"+to_string(c)+"_k"+to_string(k)+"sigma"+to_string(finalSigma)+"_ob.txt";
        
        tpg.print_graph(obFolder+suffix);
        
    }
    
    
    igraph_vector_destroy(&ak);
    
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



void generateDegreeMetrics(string dataset){
    
    delimiter='\t';
    string obfolder="/Users/dongqingxiao/Documents/uncetainGraphProject/allDataSet/obOutput/"+dataset+"/";
    string dMetricfolder="/Users/dongqingxiao/Documents/uncetainGraphProject/allDataSet/dMetricOutput/"+dataset+"/";
    
    string refFolder="/Users/dongqingxiao/Documents/uncetainGraphProject/allDataSet/input/";
    
    
    long int nv=0;
    sampleNum=2000;
    
    if(dataset=="hepth"){
        nv=14318;
    }
    if(dataset=="dblp"){
        nv=824774;
    }
    
    vector<string> datasets;
    
    if(dataset=="hepth"){
        
        datasets.push_back("rand_ob_c1.100000_k60_sigma0.000977");
        datasets.push_back("rand_ob_c1.300000_k100_sigma0.250000");
        datasets.push_back("rand_ob_c1.700000_k200_sigma17.125977");
        datasets.push_back("rand_ob_c3.000000_k300_sigma1.875977");


        datasets.push_back("greedy_ob_c1.100000_k60_sigma0.009766");
        datasets.push_back("greedy_ob_c1.300000_k100_sigma0.091797");
        datasets.push_back("greedy_ob_c2.000000_k200_sigma0.278320");
        datasets.push_back("greedy_ob_c3.000000_k300_sigma0.419922");

        datasets.push_back("GREEDY5-1_hepthc2.000000_k100sigma0.031250_ob");
        datasets.push_back("GREEDY5-1_hepthc2.000000_k60sigma0.001953_ob");
        datasets.push_back("GREEDY5-1_hepthc3.000000_k200sigma0.062500_ob");
        datasets.push_back("GREEDY5-1_hepthc3.000000_k300sigma1.000000_ob");
    }
    
    
    if(dataset=="dblp"){
//        datasets.push_back("GREEDY5-1_dblpc2.000000_k60_ob");
//        datasets.push_back("GREEDY5-1_dblpc2.000000_k100_ob");
//        datasets.push_back("GREEDY5-1_dblpc3.000000_k300_ob");
//        
//        datasets.push_back("rand_dblp_ob_c1.100000_k60");
//        datasets.push_back("rand_dblp_ob_c1.100000_k100");
//        datasets.push_back("rand_dblp_ob_c1.300000_k200");
//        datasets.push_back("rand_dblp_ob_c1.700000_k300");
//        
//        datasets.push_back("greedy_dblp_ob_c1.100000_k60searchEE");
//        datasets.push_back("greedy_dblp_ob_c1.500000_k100searchEE");
   //     datasets.push_back("greedy_dblp_ob_c1.500000_k200searchEE");
     //   datasets.push_back("greedy_dblp_ob_c1.500000_k300searchEE");
    }
    
    
    for(string data: datasets){
        
        UncertainGraph uobg=init_uncertain_OB_from_file(obfolder+data+".txt", nv);
        uobg.degreeMetricRecord(dMetricfolder+data+".txt");
        
    }
    
    // gen ref data
    
    UncertainGraph urefg=init_uncertain_from_file(refFolder+dataset+".txt");
    
    urefg.degreeMetricRecord(dMetricfolder+"ref"+".txt");
    
    
    
}



void degreMetricComparision(string dataset){
    
    string relfoler="/Users/dongqingxiao/Documents/uncetainGraphProject/allDataSet/dMetricOutput/"+dataset+"/";
    
    
    long int nv=0;
    
    if(dataset=="hepth"){
        nv=14318;
    }
    if(dataset=="dblp"){
        nv=824774;
    }
    
    string refdata="ref";
    
    
    vector<string> datasets;
    
    if(dataset=="hepth"){
        
                datasets.push_back("GREEDY5-1_hepthc2.000000_k60sigma0.001953_ob");
                datasets.push_back("GREEDY5-1_hepthc2.000000_k100sigma0.031250_ob");
                datasets.push_back("GREEDY5-1_hepthc3.000000_k200sigma0.062500_ob");
                datasets.push_back("GREEDY5-1_hepthc3.000000_k300sigma1.000000_ob");
        
        
                datasets.push_back("rand_ob_c1.100000_k60_sigma0.000977");
                datasets.push_back("rand_ob_c1.300000_k100_sigma0.250000");
                datasets.push_back("rand_ob_c1.700000_k200_sigma17.125977");
                datasets.push_back("rand_ob_c3.000000_k300_sigma1.875977");
        
        
                datasets.push_back("greedy_ob_c1.100000_k60_sigma0.009766");
                datasets.push_back("greedy_ob_c1.300000_k100_sigma0.091797");
                datasets.push_back("greedy_ob_c2.000000_k200_sigma0.278320");
                datasets.push_back("greedy_ob_c3.000000_k300_sigma0.419922");
        
        
        
        
    }
    
    if(dataset=="dblp"){
        datasets.push_back("GREEDY5-1_dblpc2.000000_k60_ob");
        datasets.push_back("GREEDY5-1_dblpc2.000000_k100_ob");
        datasets.push_back("GREEDY5-1_dblpc3.000000_k300_ob");
        
        datasets.push_back("rand_dblp_ob_c1.100000_k60");
        datasets.push_back("rand_dblp_ob_c1.100000_k100");
        datasets.push_back("rand_dblp_ob_c1.300000_k200");
        datasets.push_back("rand_dblp_ob_c1.700000_k300");
        
        datasets.push_back("greedy_dblp_ob_c1.100000_k60searchEE");
        datasets.push_back("greedy_dblp_ob_c1.500000_k100searchEE");
        datasets.push_back("greedy_dblp_ob_c1.500000_k200searchEE");
        datasets.push_back("greedy_dblp_ob_c1.500000_k300searchEE");
        
    }
    
    igraph_vector_t ref_d_metrics;
    igraph_vector_init(&ref_d_metrics,5);
    
    init_vector_file(&ref_d_metrics,relfoler+refdata+".txt");
    vector<string> metrics;
    
    metrics.push_back("NE");
    metrics.push_back("AD");
    metrics.push_back("MD");
    metrics.push_back("DV");
    metrics.push_back("PL");
    
    
    for(int i=0;i<5;i++){
        cout<<metrics[i]<<endl;
        
        for(string dataset: datasets){
            igraph_vector_t d_metrics;
            igraph_vector_init(&d_metrics,5);
            
            init_vector_file(&d_metrics,relfoler+dataset+".txt");
            
            cout<<dataset<<":"<<relative_error_metric(&ref_d_metrics,&d_metrics,i)<<endl;
            
            
            igraph_vector_destroy(&d_metrics);
        }
        
    }
    
    
    
    
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
    
    //string dataset="hepth";
//    string repDataset="GREEDY5-1_hepth";
   // w_dataset=hepth;
    
    //basic_metric(dataset);
    //randomPerturbation(dataset);
   // reliablityUtilty(dataset);
     // generateReliablityRef(dataset);
    
   // greedyPerturbation(dataset);
    //  certainObfuscation(dataset,repDataset);
    //generateReliablity(dataset);
    //reliablityComparision(dataset);
    
//    string dataset="dblp";
//    reliablityComparision(dataset);
    
    
  //  string dataset="dblp";
  // generateDegreeMetrics(dataset);
   // degreMetricComparision(dataset);
    
   graphTest();
  
    
    return 0;
}