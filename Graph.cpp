//
//  Graph.cpp
//  testGraphCplus
//
//  Created by dongqingxiao on 4/11/16.
//  Copyright © 2016 dongqingxiao. All rights reserved.
//

#include "Graph.hpp"
#include <vector>
#include <fstream>
#include <sstream>
#include <iostream>
#include <boost/pending/disjoint_sets.hpp>
#include <boost/math/distributions/normal.hpp>
#include <math.h>
#include "Help.hpp"
#include "DDCal.hpp"
#include <random>
#include <map>
#include "truncated_normal.hpp"
using boost::math::normal;


Graph::Graph(long n){
    igraph_empty(&graph, (igraph_real_t)n, IGRAPH_DIRECTED);// create empty graph
    nv=n;
    //cout<<"init nv"<<nv<<endl;
}

Graph::Graph(const Graph & obj){
    igraph_copy(&graph, &obj.graph);
    nv=obj.nv;
    ne=obj.ne;
}

Graph::Graph(Graph &obj){
    igraph_copy(&graph, &obj.graph);
    nv=obj.nv;
    ne=obj.ne;
}

Graph::Graph(UncertainGraph & obj){
    igraph_copy(&graph,&obj.graph);
    nv=obj.nv;
    ne=obj.ne;
}

Graph::~Graph(){
    //cout<<"call graph deconstuction function"<<endl;
    igraph_destroy(&graph);
}
void Graph::set_edges(igraph_vector_t *edges){
    igraph_add_edges(&graph, edges, 0);
    ne=igraph_vector_size(edges);
    ne/=2;
}

long int Graph::getNE(){
    return ne;
}

void Graph::graphStatstic(){
    igraph_vector_t degrees;
    igraph_real_t nne,ad,md;
    
    igraph_vector_init(&degrees, nv);
    igraph_degree(&graph,&degrees,igraph_vss_all(), IGRAPH_ALL, true);
    
    nne=igraph_vector_sum(&degrees);
    nne/=2;
    md=igraph_vector_max(&degrees);
    ad=2*nne/nv;
    
    printf("certain stastic \n");
    printf("the number of vertice: %li \n",nv);
    printf("NE\tAD\tMD\n");
    printf("%2f\t%2f\t%2f\n", nne,ad,md);
    igraph_vector_destroy(&degrees);
}

void Graph::metric(igraph_vector_t *res,double p){
    igraph_vector_t degrees;
    igraph_real_t md;
    int pos_md=0;
    
    igraph_vector_init(&degrees, nv);
    igraph_degree(&graph,&degrees,igraph_vss_all(), IGRAPH_ALL, true);
    md=igraph_vector_max(&degrees);
    //update
    igraph_vector_set(res,pos_md,VECTOR(*res)[pos_md]+md*p);
    
    igraph_vector_destroy(&degrees);
}

void Graph::degrees(igraph_vector_t * degRes){
    igraph_degree(&graph, degRes, igraph_vss_all(), IGRAPH_ALL, true);
    // done 
}


igraph_real_t Graph::connectedVPairs(void){
    igraph_vector_t edges;
    std::map<long int, long int> c_component;
    igraph_real_t nVpairs;
    
    
    
    nVpairs=0;
    igraph_vector_init(&edges,2*ne);
   
    igraph_get_edgelist(&graph, &edges, false);
    vector<long int> rank (nv);
    vector<long int> parent(nv);
    
    for(int i=0;i<nv;i++){
        rank[i]=i;
        parent[i]=i;
    }
    
    boost::disjoint_sets<long int *, long int * > ds(&rank[0],&parent[0]);
    
    // link via each edge
    for(long int i=0;i<ne;i++){
        long int from=VECTOR(edges)[2*i];
        long int to=VECTOR(edges)[2*i+1];
        ds.union_set(from,to);
    }
    
    
    for(long int i=0;i<nv;i++){
        long int rep=ds.find_set(i);
        ++c_component[rep];
    }
    //
//    cout<<"connected componet:"<<c_component.size()<<endl;
    
    
    for(auto entry: c_component){
        nVpairs+=permuateCal(entry.second);
    }

    
    igraph_vector_destroy(&edges);
    return nVpairs;
}




double Graph::diffconectPairAddEdge(double nFrom, double nTo){
    
    igraph_vector_t edges;
    std::map<long int, long int> c_component;
    
    double diff=0;
    
    
    
    igraph_vector_init(&edges,2*ne);
    
    igraph_get_edgelist(&graph, &edges, false);
    vector<long int> rank (nv);
    vector<long int> parent(nv);
    
    for(int i=0;i<nv;i++){
        rank[i]=i;
        parent[i]=i;
    }
    
    boost::disjoint_sets<long int *, long int * > ds(&rank[0],&parent[0]);
    
    // link via each edge
    for(long int i=0;i<ne;i++){
        long int from=VECTOR(edges)[2*i];
        long int to=VECTOR(edges)[2*i+1];
        ds.union_set(from,to);
    }
    
   
    
    long int nFromRep=ds.find_set(nFrom);
    long int nToRep=ds.find_set(nTo);
    
    if(nFromRep==nToRep){
        return 0;
    }else{
        
        // diff cluster
            for(long int i=0;i<nv;i++){
                long int rep=ds.find_set(i);
                ++c_component[rep];
            }
        
        
        long int nFromNum=c_component[nFromRep];
        long int nToNum=c_component[nToRep];
        
        
        //return 2*nFromNum*nToNum;
        
        diff=2*nFromNum*nToNum;

    }
    
    
    
    
    

    
    
    
    return diff;
 
}

void Graph::reliablity_record(igraph_vector_t *res){
    
    igraph_vector_t edges;
  
    
    igraph_vector_init(&edges,2*ne);
    vector<long int> rank (nv);
    vector<long int> parent(nv);
 
    
    igraph_get_edgelist(&graph, &edges, false);
    for(int i=0;i<nv;i++){
        rank[i]=i;
        parent[i]=i;
    }
    
    boost::disjoint_sets<long int *, long int * > ds(&rank[0],&parent[0]);
    
    // link via each edge
    for(long int i=0;i<ne;i++){
        long int from=VECTOR(edges)[2*i];
        long int to=VECTOR(edges)[2*i+1];
        ds.union_set(from,to);
    }
    
    //
    for(long int i=0;i<nv;i++){
        int rep=(int)ds.find_set(i);
        igraph_vector_set(res,i, rep);
    }
    
    rank.clear();
    parent.clear();
    
    igraph_vector_destroy(&edges);
    

}



void Graph::reliablity(igraph_vector_t *res, double p){
    
    igraph_vector_t edges;
    igraph_vector_t v_component; // record its represetnative
    igraph_vector_t c_component; // record the number of count
    
    igraph_vector_init(&edges,2*ne);
    vector<long int> rank (nv);
    vector<long int> parent(nv);
    igraph_vector_init(&v_component,nv);
    igraph_vector_init(&c_component,nv);
    
    igraph_get_edgelist(&graph, &edges, false);
    for(int i=0;i<nv;i++){
        rank[i]=i;
        parent[i]=i;
    }
    
    boost::disjoint_sets<long int *, long int * > ds(&rank[0],&parent[0]);
    
    // link via each edge
    for(long int i=0;i<ne;i++){
        long int from=VECTOR(edges)[2*i];
        long int to=VECTOR(edges)[2*i+1];
        ds.union_set(from,to);
    }
    
    //
    for(long int i=0;i<nv;i++){
        long int rep=ds.find_set(i);
        igraph_vector_set(&v_component,i,rep);
        igraph_vector_set(&c_component,rep,VECTOR(c_component)[rep]+1);
    }
    
    //record val
    for(long int i=0;i<nv;i++){
        double inc=VECTOR(c_component)[(int)(VECTOR(v_component)[i])]-1;
        inc*=p;
        igraph_vector_set(res,i,VECTOR(*res)[i]+inc);
    }

    igraph_vector_destroy(&v_component);
    igraph_vector_destroy(&c_component);
    igraph_vector_destroy(&edges);
    
    
}

void Graph::reliablity(igraph_vector_t *res, string path){
    reliablity(res, 1.0);
    write_vector_file(res, path);
}



void Graph::entropyReport(igraph_vector_t *entropyReport, igraph_real_t maxDegree){
    igraph_vector_t degs;
    igraph_vector_init(&degs, 0);
    igraph_degree(&graph, &degs, igraph_vss_all(), IGRAPH_ALL,IGRAPH_LOOPS);
 
    // iterate over all vertices
    for(long int i=0;i<nv;i++){
        long int degree=VECTOR(degs)[i];
        igraph_vector_set(entropyReport,degree,VECTOR(*entropyReport)[degree]+1);
    }
    
    
    for(long int i=0;i<maxDegree+1;i++){
        long int val=VECTOR(*entropyReport)[i];
        if(val!=0) {
            igraph_vector_set(entropyReport,i,log2(val));
        }else{
            igraph_vector_set(entropyReport,i,0);
        }
    }
    
    igraph_vector_destroy(&degs);
}

double Graph::testAgaist(igraph_vector_t *ak, int k){
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
    return lessAn/nv;
}


void Graph::testDebug(igraph_vector_t *ak, int k){
    
    igraph_vector_t counter;
    igraph_real_t lessAn;
    igraph_real_t maxDegree;
    
    lessAn=0;
    maxDegree=igraph_vector_max(ak);
    
    igraph_vector_init(&counter,maxDegree+1);
    
    for(long int i=0;i<nv;i++){
        long int degree=VECTOR(*ak)[i];
        igraph_vector_set(&counter, degree, VECTOR(counter)[degree]+1);
    }

    
    for(long int i=0;i<maxDegree+1;i++){
        if(VECTOR(counter)[i]<k){
            lessAn+=VECTOR(counter)[i];
        }
    }
    
    cout<<"number of vertices failed to obfuscated"<<lessAn<<endl;
    cout<<"epsilon"<< lessAn/nv<<endl;
}


void Graph::selfTest(int k){
    igraph_vector_t expectDeg;
    igraph_real_t  eps_res;
    igraph_vector_init(&expectDeg,nv);
    igraph_degree(&graph,&expectDeg,igraph_vss_all(), IGRAPH_ALL, IGRAPH_LOOPS);
    eps_res=testAgaist(&expectDeg, k);
    
    igraph_vector_destroy(&expectDeg);
}

/**
 *
 */
void Graph::sigmaUniquess(igraph_vector_t * uv, igraph_vector_t ak, igraph_real_t maxDegree, igraph_real_t sigma){
    igraph_vector_t s_com,freqs,degrees;
    
    normal ns(0,sigma);
    
    igraph_vector_init(&freqs,maxDegree+1);
    igraph_vector_init(&s_com,maxDegree+1);
    igraph_vector_init(&degrees,nv);
    
    igraph_degree(&graph, &degrees, igraph_vss_all(), IGRAPH_ALL,IGRAPH_LOOPS);
    
    for(long int i=0;i<nv;i++){
        long int degree=VECTOR(degrees)[i];
        igraph_vector_set(&freqs, degree, VECTOR(freqs)[degree]+1);
    }
    
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
    
    // map use adversary knowledge
    for(long int i=0;i<nv;i++){
        long int degree=VECTOR(ak)[i];
        igraph_vector_set(uv,i,VECTOR(s_com)[degree]);
    }
    
    
    igraph_vector_destroy(&s_com);
    igraph_vector_destroy(&freqs);
}


UncertainGraph Graph::generateObfuscation(igraph_real_t sigma, igraph_real_t * eps_res, igraph_vector_t * ak){
    
    UncertainGraph pg((igraph_integer_t)nv);

    igraph_real_t epsilonStart,maxDegree;
    igraph_vector_t degrees,uv,eids;
    
    igraph_matrix_t  EIndicator;
    vector<double> lowNodes;
    vector<int> edgeShuffle;
    vector<Node_UN> nodeUNs;
    
    
    
    epsilonStart=1;
    igraph_vector_init(&degrees,nv);
    igraph_degree(&graph,&degrees,igraph_vss_all(), IGRAPH_ALL,IGRAPH_LOOPS);
    maxDegree=igraph_vector_max(ak);
    igraph_vector_init(&uv,nv);
    
    
    printf("start the computation of the sigma-uniqueness of all V \n");
    sigmaUniquess(&uv, *ak, maxDegree, sigma);
    printf("finish the computation \n");
    
    
    
    int skip_count=(int) lround(epsilon*nv/2)+1;
    
    
    for(long int i=0;i<nv;i++){
        nodeUNs.push_back(Node_UN(i,VECTOR(degrees)[i],VECTOR(uv)[i]));
        lowNodes.push_back(VECTOR(uv)[i]);
    }
    
    
    sort(nodeUNs.begin(),nodeUNs.end());
    
    for(int i=0;i<skip_count;i++){
        Node_UN  nu=nodeUNs[i];
        long int pos=nu.nodeID;
        lowNodes[pos]=0;
    }
    
    discrete_distribution<long int> distribution(lowNodes.begin(),lowNodes.end());
    uniform_real_distribution<double> unDist(0.0,1.0);
    uniform_real_distribution<double> unGenDist(0.0,1.0);
    default_random_engine sgen; // random engine
    
    
    
    
    printf("number of vertices:%li \n",nv);
    
    
    igraph_matrix_init(&EIndicator,nv,nv);
    igraph_vector_init(&eids, 0);
    igraph_get_edgelist(&graph, &eids, false);
    
    
    
    /*init E indicator matrix*/
    for(long int i=0;i<2*ne;i+=2){
        long int from=VECTOR(eids)[i];
        long int to=VECTOR(eids)[i+1];
        if(from>to){
            swap(from,to);
        }
        igraph_matrix_set(&EIndicator, from, to, 1);
    }
    
    
    
    
    /*randomized geneation candiates*/
    for(int tn=0;tn<attempt;tn++){
        
        long int ce=c*ne;
        igraph_vector_t EC,ue,pe;
        long int k=0;
        long count=ne;
        igraph_matrix_t ECIndicator;
        igraph_real_t ueSum,sigma_Sum;
        ueSum=0;
        sigma_Sum=0;
        
        random_device rd;
        
        std::mt19937 gen(rd());
        
        
        // initiate EC
        igraph_vector_init(&EC, 2*ce);
        igraph_matrix_init(&ECIndicator, nv, nv);
        
        
        
        for(long int i=0;i<2*ne;i+=2){
            long int from=VECTOR(eids)[i];
            long int to=VECTOR(eids)[i+1];
            if(from>to){
                swap(from,to);
            }
            igraph_matrix_set(&ECIndicator,from,to,1);
        }
        
        
        
        printf("iteration %d\n",tn);
        printf("select %li edges from %li nodes \n", ce,nv-skip_count);
        
        // select edge
        while(true){
            
            long int u=distribution(gen);
            long int v=distribution(gen);
            
            
            if(u==v){
                continue;
            }
            
            if(u>v){
                swap(u,v);
            }
            
            if(MATRIX(EIndicator,u,v)==1){
                // if edge exist in G
                if(MATRIX(ECIndicator, u, v)==1){
                    count-=1;
                    igraph_matrix_set(&ECIndicator,u,v,0);
                }else{
                    
                }
                
            }else{
                
                if(MATRIX(ECIndicator, u, v)==0){
                    igraph_vector_set(&EC, k++,u);
                    igraph_vector_set(&EC, k++,v);
                    
                    count+=1;
                    igraph_matrix_set(&ECIndicator,u,v,1);
                }
                
            }
            
            if(count==ce){
                break;
            }
            
        }
        
        long int exist=0;
        
        for(long int i=0;i<2*ne;i+=2){
            long int from=VECTOR(eids)[i];
            long int to=VECTOR(eids)[i+1];
            if(from>to){
                swap(from,to);
            }
            
            if(MATRIX(ECIndicator,from,to)==1){
                igraph_vector_set(&EC, k++, from);
                igraph_vector_set(&EC, k++, to);
                exist+=1;
            }
            
        }
        printf("finish edge selection:%li edges \n", k);
        printf("select edge from existing edge :%li\n",exist);
        
        
        if(k!=2*ce){
            throw std::exception();
        }
        
        
        
        igraph_vector_init(&ue,ce);
        igraph_vector_init(&pe,ce);
        
        
        
        
        
        for(long int i=0;i<2*ce;i+=2){
            igraph_real_t eVal=0;
            eVal+=VECTOR(uv)[int(VECTOR(EC)[i])];
            eVal+=VECTOR(uv)[int(VECTOR(EC)[i+1])];
            eVal/=2;
            ueSum+=eVal;
            igraph_vector_set(&ue,i/2,eVal);
        }
        
        
        igraph_vector_scale(&ue,1.0/ueSum);
        
        int seed=1246789091;
        int unCount=0;
        double reSum=0;
        printf("start inject uncertainty \n");
        
        // add one random shuffle to help
        
        
        for(long int i=0;i<ce;i++){
            
            
            igraph_real_t sigma_e=sigma*ce*VECTOR(ue)[i];
            
            double w=unDist(gen);
            sigma_Sum+=sigma_e;
            igraph_real_t re=0;
            
            if(w<noise){
                re=unGenDist(sgen);
                unCount+=1;
            }else{
                re=truncated_normal_ab_sample(0.0, sqrt((double) sigma_e), 0.0, 1.0, seed);
                // re=truncated_normal_ab_sample(0.0,  sigma_e, 0.0, 1.0, seed);
            }
            
            long int from=VECTOR(EC)[2*i];
            long int to=VECTOR(EC)[2*i+1];
        
            
            if(from>to){
                swap(from,to);
            }
            
            if(MATRIX(EIndicator,from,to)==1){
                igraph_vector_set(&pe, i, 1-re);
            }else{
                igraph_vector_set(&pe, i, re);
            }
            //
            reSum+=re;
        }
        
        
        printf("finish inject uncertainty \n");
        
        printf("inject uncertainty %f \n", reSum/ce);
        
        printf("inject uncertainty over white noise %f \n", (double)unCount/ce);
        
        printf("sum of edge :%f \n", igraph_vector_sum(&pe));
        
        printf("average of sigme %f \n",sigma_Sum/ce);
        
        printf("init edge and probs \n");
        UncertainGraph pGraph((igraph_integer_t)nv);
        pGraph.set_edges(&EC);
        pGraph.set_edges_probs(&pe);
        
        igraph_vector_destroy(&ue);
        igraph_vector_destroy(&pe);
        igraph_vector_destroy(&EC);
        igraph_matrix_destroy(&ECIndicator);
        double epsilon_G=pGraph.testAgaist(ak);
        
        if(epsilon_G<epsilonStart){
            pg=pGraph;
            epsilonStart=epsilon_G;
            cout<<"find better obfuscation one :"<<epsilon_G<<endl;
        }
        
        cout<<"finish "<<tn<<" attempt" <<endl;
    }
    
    
    
    //igraph_vector_destroy(&cuv);
    igraph_vector_destroy(&degrees);
    igraph_vector_destroy(&uv);
    igraph_matrix_destroy(&EIndicator);
    
    
    *eps_res=epsilonStart;
    
    return pg;
}

UncertainGraph Graph::generateObfuscation(igraph_real_t sigma, igraph_real_t * eps_res){
    igraph_vector_t degrees;
    igraph_vector_init(&degrees, nv);
    igraph_degree(&graph,&degrees,igraph_vss_all(), IGRAPH_ALL,IGRAPH_LOOPS);
    UncertainGraph tpg=generateObfuscation(sigma, eps_res, &degrees);
    igraph_vector_destroy(&degrees);
    return tpg;
}

UncertainGraph Graph::obfuscation(){
    igraph_vector_t degrees;
    igraph_vector_init(&degrees, nv);
    igraph_degree(&graph,&degrees,igraph_vss_all(), IGRAPH_ALL,IGRAPH_LOOPS);
    UncertainGraph tpg=obfuscation(&degrees);
    igraph_vector_destroy(&degrees);
    return tpg;
}


UncertainGraph Graph::obfuscation(igraph_vector_t * ak){
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


void Graph::print_graph(string filepath){
    igraph_vector_t eids;
    igraph_vector_init(&eids,2*ne);
    
    igraph_get_edgelist(&graph, &eids, false);
    
    ofstream myfile (filepath);
    if(!myfile.is_open()){
        cout<<"sorry can't open"<<filepath<<endl;
    }
    for(long int i=0;i<ne;i++){
        igraph_real_t from, to;
        from=VECTOR(eids)[2*i];
        to=VECTOR(eids)[2*i+1];
        
        if(from>to){
            swap(from,to);
        }
        myfile<<from<<"\t"<<to<<endl;
    }
    
    myfile.close();
    igraph_vector_destroy(&eids);
    cout<<"write graph into"<<filepath<<"in edges format"<<endl;
}



/**
 * init  deterministic graph from file
 */
Graph init_from_file(string path){
    
    ifstream graphFile(path);
    string line;
    string item;
    vector<double> v;
    igraph_vector_t v_graph;
    if(graphFile.is_open()){
        while (getline(graphFile,line)) {
            stringstream ss(line); // for split line
            while(getline(ss,item,delimiter)){
                v.push_back(stol(item));
            }
        }
        
        
        graphFile.close();
    }else{
        cout<<"can't open"<<path<<endl;
    }
    
    double * array=v.data();
    igraph_vector_view(&v_graph, array, v.size());
    igraph_real_t nv=igraph_vector_max(&v_graph);
    
    Graph g(nv+1);
    g.set_edges(&v_graph);
    cout<<"init graph from"<<path<<" in edges format"<<endl;
    
    return g;
}

Graph init_from_Adj_File(string filepath){
    ifstream graphFile(filepath);
    string line;
    string item;
    vector<double> v;
    
    igraph_vector_t v_graph;
    long int cv=0;
    if(graphFile.is_open()){
        while(getline(graphFile,line)){
            stringstream ss(line);
            while(getline(ss,item,' ')){
                long int temp=cv;
                long int nv=stol(item);
                if(temp<nv){
                    v.push_back(temp);
                    v.push_back(nv);
                }
            }
            cv+=1;
        }
        graphFile.close();
    }
    
    double * array=v.data();
    igraph_vector_view(&v_graph, array, v.size());
    
    Graph pg((igraph_real_t)cv);
    pg.set_edges(&v_graph);
    
    cout<<"init graph from"<<filepath<<" in adj format"<<endl;
    
    return pg;

}



