//
//  UncertainGraph.cpp
//  testGraphCplus
//
//  Created by dongqingxiao on 4/11/16.
//  Copyright © 2016 dongqingxiao. All rights reserved.
//

#include "UncertainGraph.hpp"
#include "Graph.hpp"
#include "DDCal.hpp"
#include "Help.hpp"
#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include <boost/math/distributions/normal.hpp>
#include <random>
#include <boost/math/special_functions/round.hpp>
using boost::math::normal;
using boost::math::iround;
using namespace std;



UncertainGraph::UncertainGraph(long n){
    igraph_empty(&graph, (igraph_real_t)n, IGRAPH_DIRECTED);// create empty graph
    nv=n;
    ne=0;
}

UncertainGraph::UncertainGraph(const UncertainGraph & obj){
    cout<<"const copy constructor"<<endl;
    igraph_copy(&graph, &obj.graph);
    nv=obj.nv;
    ne=obj.ne;
    //igraph_vector_init(&pe,ne);
    igraph_vector_copy(&pe, &obj.pe);
}

UncertainGraph::UncertainGraph(UncertainGraph &obj){
    cout<<"copy constructor"<<endl;
    igraph_copy(&graph, &obj.graph);
    nv=obj.nv;
    ne=obj.ne;
    //igraph_vector_init(&pe,ne);
    igraph_vector_copy(&pe, &obj.pe);
    
}

UncertainGraph& UncertainGraph::operator=(const UncertainGraph & obj){
    cout<<"assignment constructor" <<endl;
    igraph_copy(&graph, &obj.graph);
    nv=obj.nv;
    ne=obj.ne;
    //igraph_vector_init(&pe,ne);
    igraph_vector_copy(&pe, &obj.pe);
    return *this;
}




UncertainGraph::~UncertainGraph(){
    //cout<<"call uncertain deconstructor function"<<endl;
    igraph_destroy(&graph);
    igraph_vector_destroy(&pe);

}
void UncertainGraph::set_edges(igraph_vector_t *edges){
    igraph_add_edges(&graph, edges, 0);
    ne=igraph_vector_size(edges);
    ne/=2;
    igraph_vector_init(&pe,ne);
}
void UncertainGraph::set_edges_probs(igraph_vector_t *probs){
   // igraph_cattribute_EAN_setv(&graph, "prob", probs);
    igraph_vector_copy(&pe, probs);
    
}

void UncertainGraph::entropyReport(igraph_vector_t * entropyReport,igraph_real_t maxDegree){
    
    igraph_vector_t degs, edgeProbs,s, s_entropy;
    igraph_real_t degree;
    
    
    igraph_vector_init(&degs, nv);
    igraph_degree(&graph, &degs, igraph_vss_all(), IGRAPH_ALL,IGRAPH_LOOPS);
    igraph_vector_init(&edgeProbs, ne);
    igraph_vector_init(&s,maxDegree+1);
    igraph_vector_init(&s_entropy,maxDegree+1);
    
    
    igraph_vector_copy(&edgeProbs, &pe);
    
    // iterate over all vertices
    for(long int  i=0;i<nv;i++){
        igraph_vector_t prob_v,res_v;
        igraph_vector_t eids;
        
        igraph_vector_init(&eids,nv);
        igraph_incident(&graph, &eids, (igraph_real_t) i, IGRAPH_ALL);
        igraph_vector_init(&res_v, maxDegree+1);
        degree=VECTOR(degs)[i];
        igraph_vector_init(&prob_v, degree);
        
        if(degree!=0){
            for(long int j=0;j<degree;j++){
                long int eid=VECTOR(eids)[j];
                igraph_vector_set(&prob_v, j, VECTOR(edgeProbs)[eid]);
            }
            DDCal(&prob_v, &res_v,degree,maxDegree);
            for(long int j=0;j<min(degree,maxDegree)+1;j++){
                igraph_real_t val=VECTOR(res_v)[j];
                if(val<=numeric_limits<double>::epsilon()){
                    continue;
                }
                igraph_vector_set(&s,j,VECTOR(s)[j]+VECTOR(res_v)[j]);
                igraph_real_t addEn=cal_entropy(VECTOR(res_v)[j]);
                igraph_vector_set(&s_entropy,j,VECTOR(s_entropy)[j]+addEn);
            }
        }
        else{
            igraph_vector_set(&s,0,1+VECTOR(res_v)[0]);
        }
        
        if(i%10000==0){
            cout<<"move"<<i<<endl;
        }
        
        igraph_vector_destroy(&prob_v);
        igraph_vector_destroy(&eids);
        igraph_vector_destroy(&res_v);
    }
    
    // iterator over the range of degree
    for(long int i=0;i<maxDegree+1;i++){
        igraph_real_t sVal=VECTOR(s)[i];
        if(sVal!=0){
            igraph_real_t hVal=VECTOR(s_entropy)[i];
            igraph_vector_set(entropyReport, i, log2(sVal)+(hVal/sVal));
        }else{
            igraph_vector_set(entropyReport,i,-1);
        }
    }
    //
    igraph_vector_destroy(&degs);
    igraph_vector_destroy(&edgeProbs);
    igraph_vector_destroy(&s);
    igraph_vector_destroy(&s_entropy);
    
}

double UncertainGraph::testAgaist(igraph_vector_t *ak){
    int lessAn=0;
    igraph_real_t theshold;
    
    igraph_real_t maxDegree;
    igraph_vector_t ePResult;
    
    maxDegree=igraph_vector_max(ak);
    igraph_vector_init(&ePResult, maxDegree+1);
    cout<<"cal entropy "<<endl;
    entropyReport(&ePResult,maxDegree);
    cout<<"end entropy computation"<<endl;
    
    theshold=log2(k-2);
    for(long int i=0;i<nv;i++){
        igraph_real_t val=VECTOR(*ak)[i];
        double diff=theshold-VECTOR(ePResult)[(long int) val];
        if(diff>0.01){
            lessAn+=1;
        }
    }
    cout<<lessAn<<" vertices failed to obfuscated "<<endl;
    cout<<"epsilon:"<<(double)lessAn/nv<<endl;
    igraph_vector_destroy(&ePResult);
    return (double)lessAn/nv;
}

/**
 * Basic linear statstic
 */
void UncertainGraph::graphStastic(){
    double edgePSum=igraph_vector_sum(&pe);
    cout<<"uncertain statstic"<<endl;
    cout<<"|V|:"<<nv<<endl;
    cout<<"|E|:"<<ne<<endl;
    cout<<"|E|/|V|:"<<ne/nv<<endl;
    cout<<"edge probabiblity (mean):"<<edgePSum/ne<<endl;
    cout<<"exp degree (mean):"<<edgePSum*2/nv<<endl;
}

/**
 *
 */
void UncertainGraph::reliablityUtiliy(igraph_vector_t *ruv){
    igraph_vector_t e_edge;
    igraph_vector_t n_edge;
    igraph_vector_t r_edge;
    igraph_vector_t edges;
    igraph_vector_t peOb;
    vector<long int> sedges; // record ID of some edges which is very close to 1
    
    
    
    
    igraph_vector_init(&e_edge,ne);
    igraph_vector_init(&n_edge,ne);
    igraph_vector_init(&r_edge,ne);
    igraph_vector_init(&edges,2*ne);
    igraph_vector_init(&peOb, ne);
    
    double p=1.0/sampleNum;
    double nvp=0;
    
    cout<<"use sample method to cal utility "<<endl;
    
    
    for(long int i=0;i<sampleNum;i++){
        igraph_vector_t ind;
        igraph_vector_init(&ind,ne);
        
        cout<<i<<"th sample graph"<<endl;
        UncertainGraph sg=sampleGraph(&ind);
        
        
//        cout<<"ne:"<<sg.ne<<endl;
//        cout<<"indSum:"<<igraph_vector_sum(&ind)<<endl;
    

        Graph g(sg);
        
        if(g.getNE()!=igraph_vector_sum(&ind)){
            throw std::exception();
        }

        
        nvp=g.connectedVPairs();
        
        
        nvp=(double)nvp/nv;
        nvp=(double)nvp/(nv-1);
        nvp*=p;
        
        
        for(long int i=0;i<ne;i++){
            
            if(VECTOR(ind)[i]==1){
                igraph_vector_set(&e_edge,i,VECTOR(e_edge)[i]+nvp);
                igraph_vector_set(&peOb,i, VECTOR(peOb)[i]+p);
            }else{
                igraph_vector_set(&n_edge,i,VECTOR(n_edge)[i]+nvp);
            }
            
        }
        
        igraph_vector_destroy(&ind);
        
    }
    
    
   // print_vector(&peOb,"observed pe");
   // print_vector(&pe,"real pe");
    
    cout<<"mean diff"<<cal_mean_error_vector(&pe, &peOb)<<endl;
    
    for(long int i=0;i<ne;i++){
        double p_e=VECTOR(peOb)[i]; // the edge probability
        
        /*fix computation for such boundary edge */
        double eVal=VECTOR(e_edge)[i]/p_e;
        double nVal=VECTOR(n_edge)[i]/(1-p_e);
        double diff=0;
        if(p_e<p || p_e>1-p){
            sedges.push_back(i);
            
            if(p_e<p && eVal==0.0){
                igraph_vector_set(&r_edge,i,nVal);
                continue;
            }
            
            if(p_e>1-p && nVal==0.0){
                igraph_vector_set(&r_edge,i,eVal);
                continue;
            }
            
            
        }
        
        if(eVal>nVal){
            diff=eVal-nVal;
        }
        
        igraph_vector_set(&r_edge,i,diff);
    }
    
    // just for debugging
    
    cout<<"the number of edge ~0, ~1:" <<sedges.size()<<endl;
    
    
    sampleNum=20;
    long int i=0;
    for(auto item: sedges){
        cout<<i<<"special edge"<<endl;
        cout<<"p_e"<<VECTOR(peOb)[item]<<endl;
        i+=1;
        if(VECTOR(peOb)[item]>0.5){
            igraph_vector_t rel;
            double n_val,diff,e_val;
            igraph_vector_init(&rel,nv);
            UncertainGraph sg=fixEdgeGraph(item, 0.0);
        
            sg.reliablity(&rel);
            n_val=igraph_vector_sum(&rel);
            n_val/=nv;
            diff=0;
            e_val=VECTOR(r_edge)[item];
            
            if(e_val>n_val){
                diff=e_val-n_val;
            }
            igraph_vector_set(&r_edge,item,diff);
            
            igraph_vector_destroy(&rel);
        }else{
            igraph_vector_t rel;
            double n_val,diff,e_val;
            igraph_vector_init(&rel,nv);
            UncertainGraph sg=fixEdgeGraph(item, 1.0);
            sg.reliablity(&rel);
            e_val=igraph_vector_sum(&rel);
            e_val/=nv;
            diff=0;
            n_val=VECTOR(r_edge)[item];
            
            if(e_val>n_val){
                diff=e_val-n_val;
            }
            
            igraph_vector_set(&r_edge,item,diff);
            igraph_vector_destroy(&rel);
        }
    }
    
    sampleNum=200;
    
    
    
    
    string debugFile="/Users/dongqingxiao/Documents/uncetainGraphProject/allDataSet/progTest/edgeReliablityDiff.txt";
    cout<<"write the reliablity diff of each edge"<<endl;
    write_vector_file(&r_edge, debugFile);
    cout<<"statsitc about the reliablity of edges"<<endl;
    vector_statstic(&r_edge);
    long int pos=igraph_vector_which_max(&r_edge);
    cout<<"the max val is gained by this position"<<pos<<endl;
    cout<<"prob of edge:"<<VECTOR(peOb)[pos]<<endl;
    cout<<"prob of edge real:"<<VECTOR(pe)[pos]<<endl;
    
    
    
    
    
    
    igraph_get_edgelist(&graph, &edges, false);
    
    // aggregate by nodes
    for(long int i=0;i<ne;i++){
        long int from=VECTOR(edges)[2*i];
        long int to=VECTOR(edges)[2*i+1];
        
        double rdiff=VECTOR(r_edge)[i];
        rdiff*=VECTOR(peOb)[i];
        igraph_vector_set(ruv,from, VECTOR(*ruv)[from]+rdiff);
        igraph_vector_set(ruv, to, VECTOR(*ruv)[to]+rdiff);
    }
    
    //
    igraph_vector_destroy(&e_edge);
    igraph_vector_destroy(&n_edge);
    igraph_vector_destroy(&r_edge);
    igraph_vector_destroy(&edges);
    //done
}


void UncertainGraph::reliablity(igraph_vector_t * res){
    
    double p=1.0/sampleNum;
    
    for(long int i=0;i<sampleNum;i++){
        cout<<i<<"sample graph for reliablity"<<endl;
        UncertainGraph sg=sampleGraph();
        Graph g(sg); // cast to certain graph
        g.reliablity(res, p);
    }
    //done
}


void UncertainGraph::reliablity(igraph_vector_t * res, string filepath){
    reliablity(res);
    write_vector_file(res, filepath);
}

void UncertainGraph::getDegrees(bool expected, igraph_vector_t *res){
    igraph_vector_t degrees;
    igraph_vector_init(&degrees,nv);
    if(!expected){
        igraph_degree(&graph,&degrees,igraph_vss_all(), IGRAPH_ALL,IGRAPH_LOOPS);
        igraph_vector_copy(res,&degrees);
    }else{
        //need to for each node
        igraph_vector_t edges; // store edge in (u,v)
        igraph_vector_init(&edges, 2*ne);
        
        igraph_get_edgelist(&graph, &edges, false);
        
        for(long int i=0;i<ne;i++){
            long int from=VECTOR(edges)[2*i];
            long int to=VECTOR(edges)[2*i+1];
            double p=VECTOR(pe)[i];
            igraph_vector_set(&degrees,from, VECTOR(degrees)[from]+p);
            igraph_vector_set(&degrees,to, VECTOR(degrees)[to]+p);
        }
        
        for(long int i=0;i<nv;i++){
            //igraph_vector_copy(res,&degrees);
            int val=iround(VECTOR(degrees)[i]);
            igraph_vector_set(res,i,val);
        }
    }
    
    igraph_vector_destroy(&degrees);
    
    
}

/**
 * iterates over nodes
 */
void UncertainGraph::degreeDistribution(igraph_vector_t *res, igraph_real_t maxDegree){
    igraph_vector_t degreess;

    igraph_real_t degree;
    
    igraph_vector_init(&degreess,nv);
    igraph_degree(&graph,&degreess,igraph_vss_all(), IGRAPH_ALL,IGRAPH_LOOPS);
    
    
    
    for(long int i=0;i<nv;i++){
        igraph_vector_t prob_v,res_v;
        igraph_vector_t eids;
        
        igraph_vector_init(&eids,nv);
        igraph_incident(&graph, &eids, (igraph_real_t) i, IGRAPH_ALL);
        igraph_vector_init(&res_v, maxDegree+1);
        degree=VECTOR(degreess)[i];
        igraph_vector_init(&prob_v, degree);
        
        if(degree!=0){
            
            for(long int j=0;j<degree;j++){
                long int eid=VECTOR(eids)[j];
                igraph_vector_set(&prob_v, j, VECTOR(pe)[eid]);
            }
            
            DDCal(&prob_v, &res_v,degree,maxDegree);
            igraph_vector_add(res,&res_v);
        }
        
        
        igraph_vector_destroy(&prob_v);
        igraph_vector_destroy(&res_v);
        igraph_vector_destroy(&eids);

    }
    igraph_vector_destroy(&degreess);
    
    // done
}


void UncertainGraph::sigmaUniquess(igraph_vector_t * uv, igraph_vector_t ak, igraph_real_t maxDegree, igraph_real_t sigma){
    igraph_vector_t s_com,freqs;
    igraph_vector_t degrees;
    
    normal ns(0,sigma);
    
    igraph_vector_init(&freqs,maxDegree+1);
    igraph_vector_init(&s_com,maxDegree+1);
    igraph_vector_init(&degrees,nv);
    
    //calculate frequncy via uncertain semantic
    
    degreeDistribution(&freqs, maxDegree);

    
    // when sigma is very small
    // it approximate count the number of nodes with the same degree
    // iterate over all nodes
    for(long int i=0;i<maxDegree+1;i++){
        igraph_real_t val=0.0;
        for(long int j=0;j<maxDegree+1;j++){
            double d=cal_distance(i, j);
            double tmp=boost::math::pdf(ns,d);
            if(tmp==0 || d==1){
                tmp=numeric_limits<double>::epsilon();
            }
            tmp*=VECTOR(freqs)[j];
            val+=tmp;
        }
        
        igraph_vector_set(&s_com,i,1.0/val);
    }
    
    
    for(long int i=0;i<nv;i++){
        long int degree=VECTOR(ak)[i];
        igraph_vector_set(uv,i,VECTOR(s_com)[degree]);
    }
    
    
    igraph_vector_destroy(&s_com);
    igraph_vector_destroy(&freqs);
}


UncertainGraph UncertainGraph::generateObfuscation(igraph_real_t sigma, igraph_real_t * eps_res, igraph_vector_t * ak){
    UncertainGraph g(nv);
    
    
    // how to calculate ... 
    return g;
}

UncertainGraph UncertainGraph::obfuscation(igraph_vector_t *ak){
    
    igraph_real_t sigmaLow, sigmaUpper;
    igraph_real_t ep_res, tEpsilion;
    
    UncertainGraph tGraph((igraph_integer_t)nv);
    
    sigmaLow=0;
    sigmaUpper=1;
    
    
    while(true){
        cout<<"random search for sigmaUpper="<<sigmaUpper<<endl;
        ep_res=1;
        UncertainGraph pGraph=generateObfuscation(sigmaUpper, &ep_res, ak);
        
        if(ep_res>=epsilon){
            sigmaUpper*=2;
        }else{
            tGraph=pGraph;
            tEpsilion=ep_res;
            break;
        }
        
    }
    
    while(true){
        cout<<"random search for sigmaUpper="<<sigmaUpper<<endl;
        igraph_real_t sigma=(sigmaUpper+sigmaLow)/2;
        UncertainGraph pGraph=generateObfuscation(sigmaUpper, &ep_res, ak);
        if(ep_res>=epsilon){
            sigmaLow=sigma;
        }else{
            tGraph=pGraph;
            sigmaUpper=sigma;
            tEpsilion=ep_res;
        }
        
        if((sigmaUpper-sigmaLow)<0.0001){
            break;
        }
        
        
        
        
        
    }
    
    cout<<"inject perturbation sigma :"<<sigmaUpper<<endl;
    cout<<"tolerance level epislion:"<<tEpsilion<<endl;
    
    return tGraph;
}

/*used for testing*/
UncertainGraph UncertainGraph::fixEdgeGraph(long int index,double val){
    
    UncertainGraph g(*this);
    
    VECTOR(g.pe)[index]=val;
    
    return g;
}


UncertainGraph UncertainGraph::sampleGraph(igraph_vector_t * indicator){
    UncertainGraph g(nv);
    
    igraph_vector_t edges;
    vector<double> sampleEdges;
    
    random_device rand_dev;
    mt19937     generator(rand_dev());
    uniform_real_distribution<double> unDist(0.0,1.0);
    
    igraph_vector_init(&edges, 2*ne);
    igraph_get_edgelist(&graph, &edges, false);
    
    vector<int> vIDs;
    for(int i=0;i<ne;i++){
        vIDs.push_back(i);
    }
    
    
    random_shuffle(vIDs.begin(), vIDs.end());
    
    for(int j=0;j<ne;j++){
        
        int i=vIDs[j];
        igraph_real_t e_prob=VECTOR(pe)[i];
        
        float x=unDist(generator);
        
        if(e_prob>=x){
            igraph_vector_set(indicator,i,1);
            sampleEdges.push_back(VECTOR(edges)[2*i]);
            sampleEdges.push_back(VECTOR(edges)[2*i+1]);
        }else{
            igraph_vector_set(indicator,i,0);
        }
    }
    
    
    if(sampleEdges.size()!=2*igraph_vector_sum(indicator)){
        print_vector(indicator,"the stupid indicator \n");
        throw std::exception();
    }
    
    cout<<"sample Edge.size()"<<sampleEdges.size()<<endl;
    cout<<"edge sum in sample Graph:"<<igraph_vector_sum(indicator)<<endl;
    
    
    if(sampleEdges.size()!=0){
        igraph_vector_t sEdges;
        long int ssize=sampleEdges.size();
        igraph_vector_init(&sEdges,ssize);
        // stupid world for only accept double rather than int
        for(long int i=0;i<ssize;i++){
            igraph_vector_set(&sEdges, i, sampleEdges[i]);
        }
        
        g.set_edges(&sEdges);
        igraph_vector_destroy(&sEdges);
    }
    
    
    
    sampleEdges.clear();
    vIDs.clear();
    igraph_vector_destroy(&edges);
    
    return g;

}



UncertainGraph UncertainGraph::sampleGraph(){
    UncertainGraph g(nv);
    
    igraph_vector_t edges;
    vector<double> sampleEdges;
    random_device rand_dev;
    mt19937     generator(rand_dev());
    uniform_real_distribution<double> unDist(0.0,1.0);
    
    igraph_vector_init(&edges, 2*ne);
    igraph_get_edgelist(&graph, &edges, false);
    
    
    vector<int> vIDs;
    for(int i=0;i<ne;i++){
        vIDs.push_back(i);
    }
    
    
    random_shuffle(vIDs.begin(), vIDs.end());
    
    for(int j=0;j<ne;j++){
        
        int i=vIDs[j];
        igraph_real_t e_prob=VECTOR(pe)[i];
        
        float x=unDist(generator);
        
        if(e_prob>=x){
            sampleEdges.push_back(VECTOR(edges)[2*i]);
            sampleEdges.push_back(VECTOR(edges)[2*i+1]);
        }
    }
    
    
    if(sampleEdges.size()!=0){
        igraph_vector_t sEdges;
        long int ssize=sampleEdges.size();
        igraph_vector_init(&sEdges,ssize);
        // stupid world for only accept double rather than int
        for(long int i=0;i<ssize;i++){
            igraph_vector_set(&sEdges, i, sampleEdges[i]);
        }
        
        g.set_edges(&sEdges);
        igraph_vector_destroy(&sEdges);
    }
    
    
    sampleEdges.clear();
    vIDs.clear();
    
    igraph_vector_destroy(&edges);
   
    
    return g;
}



void UncertainGraph::print_graph(string filepath){
    igraph_vector_t eids;
    igraph_vector_t edgeProbs;
    
    igraph_vector_init(&eids,2*ne);
    igraph_vector_init(&edgeProbs, ne);
    igraph_get_edgelist(&graph, &eids, false);
    igraph_vector_copy(&edgeProbs,&pe);

    ofstream myfile (filepath);
    if(!myfile.is_open()){
        cout<<"sorry can't open "<<filepath<<endl;
    }
    for(long int i=0;i<ne;i++){
        igraph_real_t from, to, p;
        
        from=VECTOR(eids)[2*i];
        to=VECTOR(eids)[2*i+1];
        p=VECTOR(edgeProbs)[i];
        
        if(from>to){
            swap(from,to);
        }
        
        myfile<<(long int )from<<"\t"<<(long int )to<<"\t"<<(double)p<<endl;
    }
    
    myfile.close();
    igraph_vector_destroy(&eids);
    igraph_vector_destroy(&edgeProbs);
    cout<<"write uncertain graph into"<<filepath<<endl;
}

UncertainGraph init_uncertain_from_file(string filepath){
    ifstream graphFile(filepath);
    string line;
    string item;
    vector<double> v;
    igraph_vector_t v_graph,e_probs;
    vector<double> edge_prob;
    

    
    
    if(graphFile.is_open()){
        while (getline(graphFile,line)) {
            stringstream ss(line); // for split line
            int i=0;
            while(getline(ss,item,delimiter)){
                if(i<2){
                    v.push_back(stol(item));
                }else{
                    edge_prob.push_back(stof(item));
                }
                
                i+=1;
            }
        }
        
        
        graphFile.close();
    }else{
        cout<<"can't open"<<filepath<<endl;
    }
    
   // double * array=v.data();
    //igraph_vector_view(&v_graph, array, v.size());
    
    
    igraph_vector_init(&v_graph, v.size());
    
    for(long int i=0;i<v.size();i++){
        igraph_vector_set(&v_graph,i,v[i]);
    }
    igraph_real_t nv=igraph_vector_max(&v_graph);
    UncertainGraph pg(nv+1);
    pg.set_edges(&v_graph);
    
    //double * parray=edge_prob.data();
    
    cout<<"init edge probs:"<<edge_prob.size()<<endl;
    //igraph_vector_view(&e_probs, parray, edge_prob.size());
    igraph_vector_init(&e_probs,edge_prob.size());
    
    for(long int i=0;i<edge_prob.size();i++){
        igraph_vector_set(&e_probs,i,edge_prob[i]);
    }
    
    pg.set_edges_probs(&e_probs);
    
    cout<<"init uncertain graph from"<<filepath<<"in edges,p format"<<endl;
    
    igraph_vector_destroy(&v_graph);
    igraph_vector_destroy(&e_probs);
    return pg;

}

UncertainGraph init_uncertain_OB_from_file(string filepath, long nv){
    ifstream graphFile(filepath);
    string line;
    string item;
    vector<double> v;
    igraph_vector_t v_graph,e_probs;
    vector<double> edge_prob;
    if(graphFile.is_open()){
        while (getline(graphFile,line)) {
            stringstream ss(line); // for split line
            int i=0;
            while(getline(ss,item,delimiter)){
                if(i<2){
                    v.push_back(stol(item));
                }else{
                    edge_prob.push_back(stof(item));
                }

                i+=1;
            }
        }


        graphFile.close();
    }

    double * array=v.data();
    igraph_vector_view(&v_graph, array, v.size());
    UncertainGraph pg((igraph_real_t) nv);
    pg.set_edges(&v_graph);

    double * parray=edge_prob.data();
    cout<<"init edge probs:"<<edge_prob.size()<<endl;
    igraph_vector_view(&e_probs, parray, edge_prob.size());

    pg.set_edges_probs(&e_probs);

    cout<<"init uncertain graph from"<<filepath<<"in edges,p format"<<"OBfuscation output"<<endl;
    
    return pg;

}





