//
//  main.cpp
//  testGraphCplus
//
//  Created by dongqingxiao on 3/19/16.
//  Copyright © 2016 dongqingxiao. All rights reserved.
//

#include <iostream>
#include <igraph.h>
#include <Help.hpp>
#include <ProbGraph.hpp>
#include <DDCal.hpp>
#include <time.h>
#include <random>
#include <map>
using namespace std;
//using boost::math::normal;

int k;
int c;
int attempt;
double epsilon;
double noise;
char delimiter;
bool uncertain;
int sampleNum;


//void remalloc_test(){
//    
////    igraph_vector_t probs;
////    igraph_real_t probs_t[]={0.7,0.9,0.8,0.8,0.1};
////    igraph_vector_init(&probs,5);
////    
//////    igraph_vector_view(&probs,probs_t,5);
////    print_vector(&probs, "prob sequence");
////    
////    igraph_vector_destroy(&probs);
//    
//    for(int i=0;i<10;i++){
//        ProbGraph pg(10000000);
//        
//        
//        
//    }
//}

void testEntropy(){
    double result=cal_entropy(0);
    double result_2=cal_entropy(0.000000001);
    printf("result:%f, result_2:%lf \n",result,result_2);
}

void testDiscreteDistribution(){
    random_device rd;
    std::mt19937 gen(rd());
    discrete_distribution<> d({4, 0, 1, 8000});
    map<int, int> m;
    for(int n=0; n<10000; ++n) {
        ++m[d(gen)];
    }
    for(auto p : m) {
        std::cout << p.first << " generated " << p.second << " times\n";
    }
}


void testIncComponent(){
    string ex_inc_path="/Users/dongqingxiao/pythonEx/probGraph/exampleGraph_test.txt";
    uncertain=false;
    delimiter='\t';
    ProbGraph pg=init_from_File(ex_inc_path);
    long int nv=pg.getNV();
    igraph_vector_t res;
    igraph_vector_init(&res,nv);
    
    pg.certain_reliablityReport_ex_inc(&res, 1);
    
    print_vector(&res,"res");
    
    
    
}




void graphStastic(ProbGraph pg){
    
    //pg.certainGraphStatstic();
//    pg.reliablityReport();
    igraph_vector_t entropyReport;
    igraph_vector_init(&entropyReport, 3);
    pg.entropyReport(&entropyReport, 2);
    print_vector(&entropyReport, "entropyReport");
}




void graphTest(ProbGraph pg){
//    igraph_vector_t deg;
//    igraph_vector_t probs;
////    igraph_vector_t entropyReport;
////    igraph_real_t maxDegree;
////    
////    igraph_vector_init(&deg, 0);
////    pg.maxDegrees(&deg);
////    print_vector(&deg, "degree");
////
//    igraph_vector_init(&probs, 0);
//    igraph_real_t probs_t[]={0.7,0.9,0.8,0.8,0.1};
//    igraph_vector_view(&probs,probs_t,5);
//    print_vector(&probs, "prob sequence");
//    
//    pg.set_edges_prob(&probs);
//    
//    pg.uncertainGraphStastic();
    
//
//    maxDegree=pg.maxDegree();
//    
//    igraph_vector_init(&entropyReport, maxDegree+1);
//    
//    pg.entropyReport(&entropyReport);
//    
//    print_vector(&entropyReport, "entropy Report");
    
   
//    igraph_vector_destroy(&deg);
//    igraph_vector_destroy(&probs);
//    igraph_vector_destroy(&entropyReport);
    
    
    
//    printf("selfTest epislon: %f \n",pg.selfTest());
    
    
//   
//    igraph_real_t epislon_t;
////////
////////    pg.certainGraphStatstic();
////////
////////    
//    ProbGraph tpg=pg.generateObfuscation(0.002, &epislon_t);
//////
//    printf("epislon:%f \n",epislon_t);
//    pg.selfTest();

   
    
    
//    tpg.uncertainGraphStastic();
    
}

void compare_test(string dataset, string inputDir, string repDir, string obDir, string repMethod, string repNum){
    
    string infilePath=inputDir+"/"+dataset+".txt";
    string obInfilePath=repDir+"/"+repMethod+repNum+"_"+dataset+".txt";
    string obOutPath=obDir+"/"+repMethod+repNum+"_"+dataset+"_ob"+".txt";
    string relOutputDir="/Users/dongqingxiao/Documents/uncetainGraphProject/allDataSet/relOutput";
    string relOutput=relOutputDir+"/"+dataset;
    long int nv=824774;
    
//    ProbGraph in_pg=init_from_File(infilePath);
//    long int cnv=in_pg.getNV();
//    cout<<"in_pg"<<endl;
//    igraph_vector_t in_rel;
//    igraph_vector_init(&in_rel,cnv);
   
//    in_pg.certainGraphStatstic();
//    in_pg.uncertainGraphStastic();
//    in_pg.uncertain_reliablityReport_ex_store(&in_rel, relOutput, "inpg");
//    
//    
//    ProbGraph rep_pg=init_from_Adj_File(obInfilePath);
//    cout<<"rep_pg"<<endl;
//    rep_pg.certainGraphStatstic();
//    igraph_vector_t rep_rel;
//    igraph_vector_init(&rep_rel,cnv);
    
//    rep_pg.certain_reliablityReport_ex_store(&rep_rel, relOutput, "reppg");
    
    ProbGraph out_pg=init_from_File(obOutPath);
    cout<<"ob_pg"<<endl;
    out_pg.uncertainGraphStastic();
    igraph_vector_t ob_rel;
    igraph_vector_init(&ob_rel,824774);
    out_pg.uncertain_reliablityReport_ex_store(&ob_rel, relOutput, dataset+"outpg");
    
    
    
//    igraph_vector_t inRel,repRel,obRel;
//    
//    string inRelPath,repRelPath,obRelPath;
//    
//    inRelPath=relOutput+"/"+"inpg"+"_rel.txt";
//    repRelPath=relOutput+"/"+"reppg"+"_rel.txt";
//    obRelPath=relOutput+"/"+"outpg"+"_rel.txt";
    
//    long int nv=824774;
//    
//    igraph_vector_init(&inRel, 824774);
//    igraph_vector_init(&repRel, 824774);
//    igraph_vector_init(&obRel, 824774);
//    
//    
//    init_vector_from_file(inRelPath, &inRel);
//    init_vector_from_file(repRelPath, &repRel);
//    init_vector_from_file(obRelPath, &obRel);
//    
//    cout<<"diff between inRel and obRel:"<<endl;
//    cout<<cal_distance_vector(&inRel, &obRel)<<endl;
//    
//    cout<<"diff between inRel and repRel:"<<endl;
//    cout<<cal_distance_vector(&inRel, &repRel)<<endl;
//    
//    cout<<"diff between repRel and obRel:"<<endl;
//    cout<<cal_distance_vector(&repRel, &obRel)<<endl;
//   
//    
//    cout<<"relative error between inRel and obRel" <<endl;
//    cout<<cal_relative_error_vector(&inRel, &obRel)<<endl;
//    cout<<"relative error between inRel and relRel" <<endl;
//    cout<<cal_relative_error_vector(&inRel, &repRel)<<endl;
//    cout<<"relative between repRel and obRel:"<<endl;
//    cout<<cal_relative_error_vector(&repRel, &obRel)<<endl;
    
    
    
    
    
}



/**
 * for generating obfucation 
 * input: certain graph
 */
void final_Test(string dataset,string inputDir, string outputDir){
    
    string inputpath=inputDir+"/"+dataset+".txt";
    string outputPath=outputDir+"/"+dataset+"_ob.txt";
    
    ProbGraph pg=init_from_Adj_File(inputpath);
   
    
    //ProbGraph pg=init_from_File(inputpath);
    pg.certainGraphStatstic();
    pg.selfTest();
    
    clock_t t;
    t=clock();
    t=clock()-t;
//    printf ("It clicks (%f seconds).\n",((float)t)/CLOCKS_PER_SEC);
    igraph_real_t eps_res;
    ProbGraph tpg=pg.generateObfuscation(0.01, &eps_res);
    printf("final eps_res:%f \n", eps_res);
    tpg.print_graph(outputPath);
    
    
}





int main(int argc, const char * argv[]) {
    
    //testProbGraph();
    k=100;
    c=2;
    attempt=1;
    epsilon=0.0001;
    noise=0.01;
    delimiter='\t';
    uncertain=true;
    sampleNum=50;
    
    
//
    igraph_i_set_attribute_table(&igraph_cattribute_table); // enable attribute handling
    string ex_path="/Users/dongqingxiao/pythonEx/probGraph/exampleGraph.txt";
    string ex_un_path="/Users/dongqingxiao/pythonEx/probGraph/exampleProbGraph.txt";
    string real_path="/Users/dongqingxiao/Documents/uncetainGraphProject/graphs/dblp.txt";
    string ex_inc_path="/Users/dongqingxiao/pythonEx/probGraph/exampleGraph_test.txt";
    
   //final_Test("GREEDY1-1_dblp","/Users/dongqingxiao/Documents/uncetainGraphProject/allDataSet/output","/Users/dongqingxiao/Documents/uncetainGraphProject/allDataSet/obOutput");
   // final_Test("dblp","/Users/dongqingxiao/Documents/uncetainGraphProject/graphs/","/Users/dongqingxiao/Documents/uncetainGraphProject/graphs/");
 compare_test("dblp","/Users/dongqingxiao/Documents/uncetainGraphProject/allDataSet/input","/Users/dongqingxiao/Documents/uncetainGraphProject/allDataSet/output",
         "/Users/dongqingxiao/Documents/uncetainGraphProject/allDataSet/obOutput","GREEDY1-","1");
    
    
   // testIncComponent();
    
    return 0;
}
