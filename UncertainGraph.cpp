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
#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include <boost/math/distributions/normal.hpp>
#include <random>

using namespace std;



UncertainGraph::UncertainGraph(long n){
    igraph_empty(&graph, (igraph_real_t)n, IGRAPH_DIRECTED);// create empty graph
    nv=n;
}

UncertainGraph::UncertainGraph(const UncertainGraph & obj){
    cout<<"const copy constructor"<<endl;
    igraph_copy(&graph, &obj.graph);
    nv=obj.nv;
    ne=obj.ne;
    igraph_vector_init(&pe,ne);
    igraph_vector_copy(&pe, &obj.pe);
}

UncertainGraph::UncertainGraph(UncertainGraph &obj){
    cout<<"copy constructor"<<endl;
    igraph_copy(&graph, &obj.graph);
    nv=obj.nv;
    ne=obj.ne;
    igraph_vector_init(&pe,ne);
    igraph_vector_copy(&pe, &obj.pe);
    
}

UncertainGraph& UncertainGraph::operator=(const UncertainGraph & obj){
    cout<<"assignment constructor" <<endl;
    igraph_copy(&graph, &obj.graph);
    nv=obj.nv;
    ne=obj.ne;
    igraph_vector_init(&pe,ne);
    igraph_vector_copy(&pe, &obj.pe);
    return *this;
}


UncertainGraph::~UncertainGraph(){
    igraph_destroy(&graph);
    igraph_vector_destroy(&pe);
}
void UncertainGraph::set_edges(igraph_vector_t *edges){
    igraph_add_edges(&graph, edges, 0);
    ne=igraph_vector_size(edges);
    ne/=2;
}
void UncertainGraph::set_edges_probs(igraph_vector_t *probs){
   // igraph_cattribute_EAN_setv(&graph, "prob", probs);
    igraph_vector_init(&pe,ne);
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
    cout<<lessAn<<"vertices failed to obfuscated "<<endl;
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

void UncertainGraph::reliablity(igraph_vector_t * res){
    
    double p=1.0/sampleNum;
    
    for(long int i=0;i<sampleNum;i++){
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
        
        igraph_vector_copy(res,&degrees);
    }
    
    igraph_vector_destroy(&degrees);
    
    
}


UncertainGraph UncertainGraph::sampleGraph(){
    UncertainGraph g(nv);
    
    igraph_vector_t edges, edgeProbs,sEdges;
    vector<double> sampleEdges;
    random_device rand_dev;
    mt19937     generator(rand_dev());
    uniform_real_distribution<double> unDist(0.0,1.0);
    
    igraph_vector_init(&edges, 2*ne);
    igraph_vector_init(&edgeProbs,ne);
    igraph_get_edgelist(&graph, &edges, false);
    igraph_vector_copy(&edgeProbs,&pe);
    
    vector<int> vIDs;
    for(int i=0;i<ne;i++){
        vIDs.push_back(i);
    }
    
    
    random_shuffle(vIDs.begin(), vIDs.end());
    
    for(int j=0;j<ne;j++){
        
        int i=vIDs[j];
        igraph_real_t e_prob=VECTOR(edgeProbs)[i];
        
        float x=unDist(generator);
        
        if(e_prob>=x){
            sampleEdges.push_back(VECTOR(edges)[2*i]);
            sampleEdges.push_back(VECTOR(edges)[2*i+1]);
        }
    }
    
    
    if(sampleEdges.size()!=0){
        igraph_vector_init(&sEdges,sampleEdges.size());
        // stupid world for only accept double rather than int
        double * array=sampleEdges.data();
        igraph_vector_view(&sEdges, array, sampleEdges.size());
        
        g.set_edges(&sEdges);
    }
    
    
    
    
    igraph_vector_destroy(&edges);
    igraph_vector_destroy(&edgeProbs);
    // igraph_vector_destroy(&sEdges);
    
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
    
    double * array=v.data();
    igraph_vector_view(&v_graph, array, v.size());
    igraph_real_t nv=igraph_vector_max(&v_graph);
    
    UncertainGraph pg(nv+1);
    pg.set_edges(&v_graph);
    
    double * parray=edge_prob.data();
    
    cout<<"init edge probs:"<<edge_prob.size()<<endl;
    igraph_vector_view(&e_probs, parray, edge_prob.size());
    
    pg.set_edges_probs(&e_probs);
    
    cout<<"init uncertain graph from"<<filepath<<"in edges,p format"<<endl;
    
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


