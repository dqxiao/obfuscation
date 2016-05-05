
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
#include <boost/pending/disjoint_sets.hpp>
#include <random>
#include <map>
#include <boost/math/special_functions/round.hpp>
#include "truncated_normal.hpp"
#include <unordered_map>
#include <boost/numeric/ublas/matrix.hpp>
using boost::math::normal;
using boost::math::iround;
//using namespace boost::numeric::ublas;
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
    
    igraph_vector_t degs,s, s_entropy;
    igraph_real_t degree;
    
    
    igraph_vector_init(&degs, nv);
    igraph_degree(&graph, &degs, igraph_vss_all(), IGRAPH_ALL,IGRAPH_LOOPS);
    igraph_vector_init(&s,maxDegree+1);
    igraph_vector_init(&s_entropy,maxDegree+1);
    
    
    //igraph_vector_copy(&edgeProbs, &pe);
    
    // iterate over all vertices
    for(long int i=0;i<nv;i++){
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
                igraph_vector_set(&prob_v, j, VECTOR(pe)[eid]);
            }
            DDCal(&prob_v, &res_v,degree,maxDegree);
            for(long int j=0;j<min(degree,maxDegree)+1;j++){
                igraph_real_t val=VECTOR(res_v)[j];
//                if(val<=numeric_limits<double>::epsilon()){
//                    continue;
//                }
                if(val==0){
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
        
        if(i%400000==0){
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
    cout<<"cal entropy, it takes a while ~10 mins "<<endl;
    entropyReport(&ePResult,maxDegree);
    cout<<"end entropy computation"<<endl;
    
    theshold=log2(k); //
    for(long int i=0;i<nv;i++){
        igraph_real_t val=VECTOR(*ak)[i];
        double diff=theshold-VECTOR(ePResult)[(long int) val];
        if(diff>0.01){
            lessAn+=1;
            //cout<<"node degree: "<< val<<endl;
            //cout<<"diff:"<<diff<<endl;
            //cout<<"log(Y): "<<VECTOR(ePResult)[(long int) val]<<endl;
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
    cout<<"|E|/|V|:"<<2*(double)ne/nv<<endl;
    cout<<"edge probabiblity (mean):"<<edgePSum/ne<<endl;
    cout<<"exp degree (mean):"<<edgePSum*2/nv<<endl;
}


void UncertainGraph::rawEstimate(igraph_vector_t *r_edge){
    
    igraph_vector_t e_edge;
    igraph_vector_t n_edge;
    igraph_vector_t peOb;
    
    igraph_vector_init(&e_edge,ne);
    igraph_vector_init(&n_edge,ne);
    igraph_vector_init(&peOb, ne);
    
    double p=1.0/sampleNum;
    double nvp=0;
    
    cout<<"use sample method to cal utility "<<endl;
    
    
    for(long int i=0;i<sampleNum;i++){
        igraph_vector_t ind;
        igraph_vector_init(&ind,ne);
        
        cout<<i<<"th sample graph"<<endl;
        UncertainGraph sg=sampleGraph(&ind);
        
        
        Graph g(sg);
        
        if(g.getNE()!=igraph_vector_sum(&ind)){
            throw std::exception();
        }
        
        
        nvp=g.connectedVPairs();
        
        
        nvp=(double)nvp/nv;
        //nvp=(double)nvp/(nv-1);
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
        
        
        
        if(p_e<10*p){
            //sedges.push_back(i);
            igraph_vector_set(r_edge,i,nVal);
            continue;
        }
        
        if(p_e>1-10*p){
            //sedges.push_back(i);
            igraph_vector_set(r_edge,i,eVal);
            continue;
        }
        
        
        
        
        if(eVal>nVal && nVal!=0){
            diff=eVal-nVal;
            igraph_vector_set(r_edge,i,diff);
        }else{
            igraph_vector_set(r_edge,i,0);
        }
        
        
        
    }

    
//    string debugFile="/Users/dongqingxiao/Documents/uncetainGraphProject/allDataSet/progTest/edgeReliablityDiff.txt";
//    cout<<"write the reliablity diff of each edge"<<endl;
//    write_vector_file(r_edge, debugFile);
//    cout<<"statsitc about the reliablity of edges"<<endl;
//    vector_statstic(r_edge);
//    long int pos=igraph_vector_which_max(r_edge);
//    cout<<"the max val is gained by this position"<<pos<<endl;
//    cout<<"prob of edge:"<<VECTOR(peOb)[pos]<<endl;
//    cout<<"prob of edge real:"<<VECTOR(pe)[pos]<<endl;
    
    igraph_vector_destroy(&e_edge);
    igraph_vector_destroy(&n_edge);
    
}


void UncertainGraph::maxConnectedConponent(igraph_vector_t * re_edges, igraph_vector_t * re_pe, igraph_vector_t * startPos, long int * cnums){
    
    igraph_vector_t edges, v_rep;
    
    igraph_vector_init(&edges,2*ne);
    igraph_get_edgelist(&graph, &edges, false);
    igraph_vector_init(&v_rep,nv);
    
    vector<long int> rank (nv);
    vector<long int> parent(nv);
    
    vector<ProbEdge> probEdges;
    
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
        long int vrep=ds.find_set(i);
        igraph_vector_set(&v_rep,i,vrep);
    }
    
    for(long int i=0;i<ne;i++){
        long int from=VECTOR(edges)[2*i];
        long int to=VECTOR(edges)[2*i+1];
        double   p=VECTOR(pe)[i];
        
        long int from_rep=VECTOR(v_rep)[from];
        long int to_rep=VECTOR(v_rep)[to];
        
        if(from_rep!=to_rep){
            throw std::exception();
        }
        //sorted by rep, from, to
        probEdges.push_back(ProbEdge(from,to,from_rep,p));
    }
    
    sort(probEdges.begin(),probEdges.end());
    

    long int repCount=0;
    long c_rep=-1;
    
    for(long int i=0;i<ne;i++){
        ProbEdge pEdge=probEdges[i];
        
        if(pEdge.rep!=c_rep){
            igraph_vector_set(startPos,repCount,i);
            c_rep=pEdge.rep;
            repCount+=1;
        }
        
        
        igraph_vector_set(re_edges,2*i,pEdge.from);
        igraph_vector_set(re_edges,2*i,pEdge.to);
        igraph_vector_set(re_pe,i,pEdge.pe);
        
       
    }
    
    igraph_vector_set(startPos, repCount, ne);
    
    *cnums=repCount;
    
    cout<<"repCount:"<<repCount<<endl;
    cout<<"ne"<<ne<<endl;
    cout<<"last startPos:"<<VECTOR(*startPos)[repCount-1]<<endl;
   
   
    
    
    

    rank.clear();
    parent.clear();
    probEdges.clear();
    igraph_vector_destroy(&edges);
    igraph_vector_destroy(&v_rep);
    
    
    
}

double UncertainGraph::diffconectPairAddEdge(double from, double to){
    
    double diff=0.0;
    double p=1.0/sampleNum;
    for(long int i=0;i<sampleNum;i++){
        UncertainGraph sg=sampleGraph();
        Graph g(sg);
        
        if(i%10==0){
            cout<<i<<" diff connect pair add edge"<<endl;
        }
        
        diff=diff+(g.diffconectPairAddEdge(from, to)*p);
        
    }
 
    return diff;
}



void UncertainGraph::reliablityUtiitySubgraphCall(igraph_vector_t * res){
    igraph_vector_t e_edge,n_edge,pe_ob ;
    
    
  
    igraph_vector_init(&e_edge,ne);
    igraph_vector_init(&n_edge,ne);
    igraph_vector_init(&pe_ob,ne);
   
    
    double p=1.0/sampleNum;
    for(long int i=0;i<sampleNum;i++){
        igraph_vector_t indicator;
        igraph_vector_init(&indicator,ne);
        
        if(i%10==0){
            cout<<i<<"utiltiy via subgraph call"<<endl;
        }
        
        UncertainGraph usg=sampleGraph(&indicator);
        Graph g(usg);
        
        double nvp=(double)g.connectedVPairs();
        nvp*=p;
        
        for(long int i=0;i<ne;i++){
            if(VECTOR(indicator)[i]==1){
                igraph_vector_set(&e_edge,i,VECTOR(e_edge)[i]+nvp);
                igraph_vector_set(&pe_ob,i,VECTOR(pe_ob)[i]+p);
            }else{
                igraph_vector_set(&n_edge,i,VECTOR(n_edge)[i]+nvp);
            }
        }
        
        igraph_vector_destroy(&indicator);
        
    }
    
    
    for(long int i=0;i<ne;i++){
        double spe=VECTOR(pe_ob)[i];
        
        if(spe>1-10*p){
            continue;
        }
        
        if(spe<10*p){
            continue;
        }
        
        double eval=VECTOR(e_edge)[i]/spe;
        double nval=VECTOR(n_edge)[i]/(1-spe);
        
        
        double diff=eval-nval;
        
        igraph_vector_set(res,i,diff);
        
        
        
        
    }
    
    
    igraph_vector_destroy(&e_edge);
    igraph_vector_destroy(&n_edge);
    igraph_vector_destroy(&pe_ob);
    
}





void UncertainGraph::reliablityUtilitySubgraph(igraph_vector_t *ruv){
    
    igraph_vector_t re_edge,re_pe,start_pos; // reoriginzed edges, pe, and start pos
    igraph_vector_t relDiff_edge; // relDiff edges
    long int cnums;
    cnums=0;
    
    igraph_vector_init(&re_edge,2*ne);
    igraph_vector_init(&re_pe,ne);
    igraph_vector_init(&start_pos,nv);
    maxConnectedConponent(&re_edge,&re_pe,&start_pos,&cnums);
    
    
    igraph_vector_init(&relDiff_edge,ne);
    
    
    cout<<"there are "<<cnums<<" components"<<endl;
    
    // before debuging
    igraph_vector_t numComponents;
    igraph_vector_init(&numComponents,cnums);
    
    for(long int i=0;i<cnums;i++){
        long int start=VECTOR(start_pos)[i];
        long int end=VECTOR(start_pos)[i+1];
        igraph_vector_set(&numComponents,i,end-start);
    }
    
    write_vector_file(&numComponents, "/Users/dongqingxiao/Documents/uncetainGraphProject/allDataSet/progTest/numEdge.txt");
    
    vector_statstic(&numComponents);
    
    igraph_vector_destroy(&numComponents);
    
    
    
    for(long int i=0;i<cnums;i++){
        long int start=VECTOR(start_pos)[i];
        long int end=VECTOR(start_pos)[i+1];
        
        // only one edge inside
        cout<<i<<"th components"<<endl;
        if(end-start<100){
            
            for(long int pos=start;pos<end;pos++){
                igraph_vector_set(&relDiff_edge,pos,0);
            }

        }else{
            
            
            long int size=end-start;
            igraph_vector_t subedges, subpes;
            igraph_vector_t subrdiff;
            long int sedges; // some edges whose probability ~0, ~1
            
            igraph_vector_init(&subedges,2*size);
            igraph_vector_init(&subpes,size);
            igraph_vector_init(&subrdiff,size);
            
            if(size<500){
            
                sampleNum=100;
            }else{
                sampleNum=1000;
                if(size/1000>5){
                    sampleNum=10000;
                }
            }
            double p=1.0/sampleNum;
            sedges=0;
            
            
            cout<<size<<"edges"<<endl;
            
            for(long int i=0;i<size;i++){
                long relPos=i+start;
                
                igraph_vector_set(&subedges,2*i,VECTOR(re_edge)[2*relPos]);
                igraph_vector_set(&subedges,2*i+1,VECTOR(re_edge)[2*relPos+1]);
                igraph_vector_set(&subpes,i,VECTOR(re_pe)[relPos]);
                
                if(VECTOR(re_pe)[relPos]>0.89 || VECTOR(re_pe)[relPos]<0.11){
                    sedges+=1;
                }
            }
        
            
            
            UncertainGraph usubg(nv);
            
            usubg.set_edges(&subedges);
            usubg.set_edges_probs(&subpes);
            
            
            
            
           
            
            
            usubg.reliablityUtiitySubgraphCall(&subrdiff);
            
            
            //refine work local pos
            
            if(size<100){
                sampleNum=10;
            }else{
                sampleNum=100;
            }
            
            
            cout<<"it take a while for each iteration for refined"<<endl;
            cout<<sedges<<"~0 ~1 edges in"<<endl;
            
            
            if(sedges==0){
                continue;
            }
            for(long int i=0;i<size;i++){
                
                if(VECTOR(subpes)[i]>0.89 || VECTOR(subpes)[i]<0.11){
//                    UncertainGraph ufsg=usubg.fixEdgeGraph(i, 0.0);
//                    
//                    long int from=VECTOR(subedges)[2*i];
//                    long int to=VECTOR(subedges)[2*i+1];
//                    
//                    double subedif=ufsg.diffconectPairAddEdge(from, to);
                    
                    igraph_vector_set(&subrdiff,i,0);
                    
                    continue;
                }
                
                
                
            }
            
            cout<<"merge result into real edges"<<endl;
            
            for(long int i=0;i<size;i++){
                long relPos=i+start;
                igraph_vector_set(&relDiff_edge,relPos,VECTOR(subrdiff)[i]);
            }
            
        
            
            igraph_vector_destroy(&subedges);
            igraph_vector_destroy(&subpes);
            igraph_vector_destroy(&subrdiff);
        
        }
        
        // finish=0
        
        
        
        
        
      
        
    }
    
    string debugRefFile="/Users/dongqingxiao/Documents/uncetainGraphProject/allDataSet/progTest/edgeReliablityDiffRef.txt";
    cout<<"write the reliablity diff of each edge"<<endl;
    write_vector_file(&relDiff_edge, debugRefFile);
    cout<<"statsitc about the reliablity of edges"<<endl;
    vector_statstic(&relDiff_edge);
    long int posRef=igraph_vector_which_max(&relDiff_edge);
    cout<<"the max val is gained by this position"<<posRef<<endl;
    cout<<"prob of edge real:"<<VECTOR(pe)[posRef]<<endl;
    
    
    
    igraph_vector_destroy(&re_edge);
    igraph_vector_destroy(&re_pe);
    igraph_vector_destroy(&start_pos);
    igraph_vector_destroy(&relDiff_edge);
    
    
}

// just for current use
//

void UncertainGraph::reliablityUtiliyInit(igraph_vector_t *ruv){
    
    // file path
    string filepath="/Users/dongqingxiao/Documents/uncetainGraphProject/allDataSet/relOutput/dblp/reliablityNode/nodeRelDiff_E.txt";
    
    init_vector_file(ruv, filepath);
    reverseVector(ruv);
    cout<<"reverse reliablity utility vector"<<endl;
    vector_statstic(ruv);
    
    // done
}


void UncertainGraph::aggregateReliablutyDiffEE(igraph_vector_t * ruv, igraph_vector_t * rue){
    
    igraph_vector_t edges;
    igraph_vector_init(&edges,2*ne);
    
    igraph_get_edgelist(&graph, &edges, false);
    
    for(long int i=0;i<ne;i++){
        int from=VECTOR(edges)[2*i];
        int to=VECTOR(edges)[2*i+1];
        double val=VECTOR(*rue)[i];
        igraph_vector_set(ruv,from, VECTOR(*ruv)[from]+val);
        igraph_vector_set(ruv,to,VECTOR(*ruv)[to]+val);
    }
    
    igraph_vector_scale(ruv, 1.0/(nv-1)); // scale to 
    igraph_vector_destroy(&edges);
}


void UncertainGraph::aggregateReliablityDiff(igraph_vector_t *v_rep, igraph_vector_t *indicator, igraph_vector_t *res, double p){
    
    
    igraph_vector_t c_res, edges;
    igraph_vector_init(&c_res,nv);
    igraph_vector_init(&edges,2*ne);
    
    igraph_get_edgelist(&graph, &edges, false);
    
    
    // thanks god
    
    std::map<int, int> c_componet; // component: # of v in inside
    long int allPairs=0;
    // aggregate
    for(int i=0;i<nv;i++){
        int rep=(int)VECTOR(*v_rep)[i];
        ++c_componet[rep];
    }
    
    
    for(auto item:c_componet){
        allPairs+=item.second*item.second;
    }
    
    //
    
    for(int i=0;i<nv;i++){
        int rep=(int)VECTOR(*v_rep)[i];
        int repNum=c_componet[rep];
        
        long int diffVal=allPairs-(repNum*repNum);
        diffVal*=2*repNum;
        
        igraph_vector_set(&c_res,i,diffVal);
    }
    
    
    // fixed the problem for existing edges
    for(int i=0;i<ne;i++){
        
        int from=VECTOR(edges)[2*i];
        int to=VECTOR(edges)[2*i+1];
        
        int repFrom=VECTOR(*v_rep)[from];
        int repTo=VECTOR(*v_rep)[to];
        
        if(repFrom==repTo){
            // do no thing
        }else{
            int repNumFrom=c_componet[repFrom];
            int repNumTo=c_componet[repTo];
            long int refDiff=2*repNumFrom*repNumTo;
            igraph_vector_set(&c_res,from, VECTOR(c_res)[from]-refDiff);
            igraph_vector_set(&c_res,to, VECTOR(c_res)[to]-refDiff);
        }
    }
    
    
    // done
    
    
    
    
    
    igraph_vector_scale(&c_res, p);
    igraph_vector_add(res,&c_res);
    igraph_vector_destroy(&c_res);
    igraph_vector_destroy(&edges);
    c_componet.clear();
    
}



void UncertainGraph::reliablityUtiltyDiff(igraph_vector_t *ruv){
    
    double p=1.0/sampleNum;
    
    for(int i=0;i<sampleNum;i++){
        igraph_vector_t ind, vRep;
        igraph_vector_init(&ind,ne);
        igraph_vector_init(&vRep,nv);
       
        
        cout<<i<<"sample for reliablity Utility difference"<<endl;
        UncertainGraph sg=sampleGraph(&ind);
        Graph g(sg); // cast to certain graph
        
        
        g.connectedComponent(&vRep);
        
        aggregateReliablityDiff(&vRep, &ind, ruv, p);
        
        
        igraph_vector_destroy(&ind);
        igraph_vector_destroy(&vRep);
    }
    
    double scale_nv_p=1.0/(nv*(nv-1));
    
    igraph_vector_scale(ruv, scale_nv_p);
}


void UncertainGraph::reliablityUtiliy(igraph_vector_t *ruv){
   
    igraph_vector_t r_edge;
    igraph_vector_t edges;
    
    vector<long int> sedges; // record ID of some edges which is very close to 1
    

  
    igraph_vector_init(&r_edge,ne);
    igraph_vector_init(&edges,2*ne);
    
    
    
    // just for debugging
    cout<<"the number of edge ~0, ~1:" <<sedges.size()<<endl;
    
    init_vector_file(&r_edge,"/Users/dongqingxiao/Documents/uncetainGraphProject/allDataSet/progTest/edgeReliablityDiff.txt");
    
    
    for(long int i=0;i<ne;i++){
        double p_v=VECTOR(pe)[i];
        if(p_v>1-10.0/sampleNum){
            sedges.push_back(i);
        }
    }
    

    
    
    
    sampleNum=10;
    long int i=0;
    for(auto item: sedges){
        
        
        cout<<i<<"special edge"<<endl;
        cout<<"p_e"<<VECTOR(pe)[item]<<endl;
        i+=1;
        if(VECTOR(pe)[item]>0.5){
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
    
    
    
    
    string debugRefFile="/Users/dongqingxiao/Documents/uncetainGraphProject/allDataSet/progTest/edgeReliablityDiffRef.txt";
    cout<<"write the reliablity diff of each edge"<<endl;
    write_vector_file(&r_edge, debugRefFile);
    cout<<"statsitc about the reliablity of edges"<<endl;
    vector_statstic(&r_edge);
    long int posRef=igraph_vector_which_max(&r_edge);
    cout<<"the max val is gained by this position"<<posRef<<endl;
   // cout<<"prob of edge:"<<VECTOR(peOb)[posRef]<<endl;
    cout<<"prob of edge real:"<<VECTOR(pe)[posRef]<<endl;
    
    
    
    
    
    
    igraph_get_edgelist(&graph, &edges, false);
    
    // aggregate by nodes
    for(long int i=0;i<ne;i++){
        long int from=VECTOR(edges)[2*i];
        long int to=VECTOR(edges)[2*i+1];
        
        double rdiff=VECTOR(r_edge)[i];
        rdiff*=VECTOR(pe)[i];
        igraph_vector_set(ruv,from, VECTOR(*ruv)[from]+rdiff);
        igraph_vector_set(ruv, to, VECTOR(*ruv)[to]+rdiff);
    }
    
    //
    
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

void UncertainGraph::reliablity_record(long int csampleNum, string filepath){
    boost::numeric::ublas::matrix<int> m (nv,csampleNum); //
    for(long int i=0;i<csampleNum;i++){
        cout<<i<<"sample Graph for reliablity record" <<endl;
        
        UncertainGraph sg=sampleGraph();
        Graph g(sg);
        igraph_vector_t sc;
        igraph_vector_init(&sc,nv);
        g.reliablity_record(&sc);
        
        for(long int index=0;index<nv;index++){
            m(index,i)=VECTOR(sc)[index];
        }
        igraph_vector_destroy(&sc);
        
    }
    
    write_boostMatrix_file(m, filepath); //
    
    m.clear(); //clear
  
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
        
        igraph_vector_destroy(&edges);
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


void UncertainGraph::aggregateAK(igraph_vector_t *res, igraph_real_t maxDegree, igraph_vector_t *ak){
    
    //
    for(long int i=0;i<nv;i++){
        long int val=VECTOR(*ak)[i];
        igraph_vector_set(res,val,VECTOR(*res)[val]+1);
    }
    //done
}





void UncertainGraph::sigmaUniquess(igraph_vector_t * uv, igraph_vector_t ak, igraph_real_t maxDegree, igraph_real_t sigma){
    igraph_vector_t s_com,freqs;
    //igraph_vector_t degrees;
    
    
    
    normal ns(0,sigma);
    
    igraph_vector_init(&freqs,maxDegree+1);
    igraph_vector_init(&s_com,maxDegree+1);
   // igraph_vector_init(&degrees,nv);
    

    
    aggregateAK(&freqs, maxDegree, &ak);
    
    // when sigma is very small
    // it approximate count the number of nodes with the same degree
    // iterate over all nodes
    for(long int i=0;i<maxDegree+1;i++){
        igraph_real_t val=0.0;
        for(long int j=0;j<maxDegree+1;j++){
            double d=cal_distance(i, j);
            double tmp=boost::math::pdf(ns,d);
//            if(tmp==0 || d==1){
//                tmp=numeric_limits<double>::epsilon();
//            }
            tmp*=VECTOR(freqs)[j];
            val+=tmp;
        }
        
        igraph_vector_set(&s_com,i,1.0/val);
    }
    
    

    
    
    for(long int i=0;i<nv;i++){
        long int degree=VECTOR(ak)[i];
        igraph_vector_set(uv,i,VECTOR(s_com)[degree]);
    }
    
    
    //just for debugging purpose
    
    //write_vector_file(uv, "/Users/dongqingxiao/Documents/uncetainGraphProject/allDataSet/progTest/uniqness.txt");
    //write_vector_file(&ak, "/Users/dongqingxiao/Documents/uncetainGraphProject/allDataSet/progTest/adversary.txt");
    
    
    
    
    igraph_vector_destroy(&s_com);
    igraph_vector_destroy(&freqs);

}


UncertainGraph UncertainGraph::randomGenerateObfuscation(igraph_real_t sigma, igraph_real_t *eps_res, igraph_vector_t *ak){
    
    UncertainGraph pg(nv);
    
    igraph_real_t epsilonStart,maxDegree;
    igraph_vector_t degrees,uv,eids;
    
   
    vector<double> lowNodes;
    vector<int> edgeShuffle;
    vector<Node_UN> nodeUNs;
    
    
    
    epsilonStart=1;
    igraph_vector_init(&degrees,nv);
    igraph_degree(&graph,&degrees,igraph_vss_all(), IGRAPH_ALL,IGRAPH_LOOPS);
    maxDegree=igraph_vector_max(ak);
    igraph_vector_init(&uv,nv);
    
    
    printf("start the computation of the sigma-uniqueness of all V \n");
    
    double unique_sigma=0.4;
    sigmaUniquess(&uv, *ak, maxDegree, unique_sigma);
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
    
    
  
    igraph_vector_init(&eids, 0);
    igraph_get_edgelist(&graph, &eids, false);
    
    
    unordered_map<Edge, double>EIndicator;
    
    for(long int i=0;i<ne;i++){
        long int from=VECTOR(eids)[2*i];
        long int to=VECTOR(eids)[2*i+1];
        double p_e=VECTOR(pe)[i];
        
        if(from>to){
            swap(from,to);
        }
        // use this one as matrix
        EIndicator[Edge(from,to)]=p_e;
    }

    
    
    
    
    
    /*randomized geneation candiates*/
    for(int tn=0;tn<attempt;tn++){
        
        long int ce=lround(c*ne);
        igraph_vector_t EC,ue,tpe; // just for temporary uncertain graph
        long int k=0;
        long count=ne;
       
        
        
        igraph_real_t ueSum,sigma_Sum;
        ueSum=0;
        sigma_Sum=0;
        
        random_device rd;
        
        std::mt19937 gen(rd());
        

        igraph_vector_init(&EC, 2*ce);
        
        unordered_map<Edge, double> ECIndicator;
        for(long int i=0;i<ne;i++){
            long int from=VECTOR(eids)[2*i];
            long int to=VECTOR(eids)[2*i+1];
        
            if(from>to){
                swap(from,to);
            }
            ECIndicator[Edge(from,to)]=1;
        }
//
        
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
            
            double pickP=unDist(gen);

        
            double p_uv=0;
            
            auto search=EIndicator.find(Edge(u,v));
            if(search!=EIndicator.end()){
                p_uv=search->second;
            }
            
            
            int existEC=0;
            auto ecSearch=ECIndicator.find(Edge(u,v));
            // exist or not
            if(ecSearch!=ECIndicator.end()){
                existEC=ecSearch->second;
            }
            
            if(p_uv>pickP){
                
                
                // remove edge
                if(existEC==1){
                    count-=1;
                    ecSearch->second=0;
                }else{
                }
            }else{
                // add one edge
                
                if(existEC==0 && p_uv==0){
                    igraph_vector_set(&EC, k++,u);
                    igraph_vector_set(&EC, k++,v);

                    count+=1;
                    
                    if(ecSearch!=ECIndicator.end()){
                        ecSearch->second=1;
                    }else{
                        ECIndicator[Edge(u,v)]=1;
                    }
                }

            }
            
            
            if(count==ce){
                break;
            }
            
        }
        
        long int exist=0;
        
        for(long int i=0;i<ne;i++){
            long int from=VECTOR(eids)[2*i];
            long int to=VECTOR(eids)[2*i+1];
            if(from>to){
                swap(from,to);
            }
            
            int existEC=0;
            auto ecSearch=ECIndicator.find(Edge(from,to));
            
            if(ecSearch!=ECIndicator.end()){
                existEC=ecSearch->second;
            }

            
            if(existEC==1){
                igraph_vector_set(&EC, k++, from);
                igraph_vector_set(&EC, k++, to);
                exist+=1;
            }
            
        }
        
       
        printf("finish edge selection:%li edges \n", k);
        printf("select edge from existing edge :%li\n",exist);
        
        printf("size of EC : %li\n", igraph_vector_size(&EC));
        printf("cout:%li, ce:%li, k:%li \n",count,ce,k);
        
        
        
        if(k!=2*ce){
            throw std::exception();
        }
        
        
        
        igraph_vector_init(&ue,ce);
        igraph_vector_init(&tpe,ce);
        
        
        
        
        
        for(long int i=0;i<ce;i++){
            igraph_real_t eVal=0;
            eVal+=VECTOR(uv)[int(VECTOR(EC)[2*i])];
            eVal+=VECTOR(uv)[int(VECTOR(EC)[2*i+1])];
            eVal/=2;
            ueSum+=eVal;
            igraph_vector_set(&ue,i,eVal);
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
            
            
            long int from=VECTOR(EC)[2*i];
            long int to=VECTOR(EC)[2*i+1];
            
            
            
            
            if(from>to){
                swap(from,to);
            }
            
            //
            double old_p_e=0;
            
            auto search=EIndicator.find(Edge(from,to));
            if(search!=EIndicator.end()){
                old_p_e=search->second;
            }
            
            
            if(w<noise){
                re=unGenDist(sgen);
                unCount+=1;
                re*=(2*(0.5-old_p_e));
            }else{
                re=truncated_normal_ab_sample(0.0, sigma_e, 0.0, 1.0, seed);
                re*=(2*(0.5-old_p_e));
                
            }
            
            // generate the solution
            old_p_e+=re;
            
            if(old_p_e<0 || old_p_e>1){
                throw std::exception();
            }
            igraph_vector_set(&tpe,i,old_p_e);
            
    
            reSum+=re;
        }
        
        
        printf("finish inject uncertainty \n");
        
        printf("inject uncertainty %f \n", reSum/ce);
        
        printf("inject uncertainty over white noise %f \n", (double)unCount/ce);
        
        printf("sum of edge :%f \n", igraph_vector_sum(&tpe));
        
        printf("average of sigme %f \n",sigma_Sum/ce);
        
        printf("init edge and probs \n");
        UncertainGraph pGraph((igraph_integer_t)nv);
        pGraph.set_edges(&EC);
        pGraph.set_edges_probs(&tpe);
        
        igraph_vector_destroy(&ue);
        igraph_vector_destroy(&tpe);
        igraph_vector_destroy(&EC);
        ECIndicator.clear();
        
        // check this graph
        pGraph.graphCheck();
        
        double epsilon_G=pGraph.testAgaist(ak);
        
        if(epsilon_G<epsilonStart){
            pg=pGraph;
            epsilonStart=epsilon_G;
            cout<<"find better obfuscation one :"<<epsilon_G<<endl;
        }
        
        cout<<"finish "<<tn<<" attempt" <<endl;
        
        
    }
    
    
    
    //igraph_vector_destroy(&cuv);
    EIndicator.clear();
    igraph_vector_destroy(&degrees);
    igraph_vector_destroy(&uv);
    igraph_vector_destroy(&eids);
    lowNodes.clear();
    edgeShuffle.clear();
    nodeUNs.clear();
    
    
    *eps_res=epsilonStart;
    
    
    
    
    
    return pg;
    
}


UncertainGraph UncertainGraph::greedyGenerateObfuscation(igraph_real_t sigma, igraph_real_t *eps_res, igraph_vector_t *ak){
    
    UncertainGraph pg(nv);
    
    igraph_real_t epsilonStart,maxDegree;
    igraph_vector_t degrees,uv,eids;
    
    
    vector<double> lowNodes;
    vector<int> edgeShuffle;
    vector<Node_UN> nodeUNs;
    
    
    
    epsilonStart=1;
    igraph_vector_init(&degrees,nv);
    igraph_degree(&graph,&degrees,igraph_vss_all(), IGRAPH_ALL,IGRAPH_LOOPS);
    maxDegree=igraph_vector_max(ak);
    igraph_vector_init(&uv,nv);
    
    
    printf("start the computation of the sigma-uniqueness of all V \n");
    
    double unique_sigma=0.4;
    sigmaUniquess(&uv, *ak, maxDegree, unique_sigma);
    printf("finish the computation \n");
    
 
    
    
    igraph_vector_t ruv;
    igraph_vector_init(&ruv,nv);
    reliablityUtiliyInit(&ruv); // done
    
    
    
    
    
    
    int skip_count=(int) lround(epsilon*nv/2)+1;
    
    // just for debugging
    igraph_vector_t  fuv;
    igraph_vector_init(&fuv,nv);
    
    for(long int i=0;i<nv;i++){
        
       double val=featureCombine(VECTOR(uv)[i], VECTOR(ruv)[i]);
        
        nodeUNs.push_back(Node_UN(i,VECTOR(degrees)[i],val));
        lowNodes.push_back(val);
    }
    
    // for debugging
   
    
    
    sort(nodeUNs.begin(),nodeUNs.end());
    
    for(int i=0;i<skip_count;i++){
        Node_UN  nu=nodeUNs[i];
        long int pos=nu.nodeID;
       // cout<<"give up pos"<<pos<<endl;
        lowNodes[pos]=0;
    }
    
    
    //just for curiousity
    
//    for(int i=0;i<nv;i++){
//        if(lowNodes[i]!=0){
//            double val=lowNodes[i];
//            val/=VECTOR(ruv)[i];
//           // val*=1-VECTOR(ruv)[i];
//            lowNodes[i]=val;
//        }
//    }
    
    
  //  write_vector_file(&fuv, "/Users/dongqingxiao/Documents/uncetainGraphProject/allDataSet/relOutput/dblp/reliablityNode/nodechoic.txt");
    
    
    discrete_distribution<long int> distribution(lowNodes.begin(),lowNodes.end());
    uniform_real_distribution<double> unDist(0.0,1.0);
    uniform_real_distribution<double> unGenDist(0.0,1.0);
    default_random_engine sgen; // random engine
    
    
    
    
    printf("number of vertices:%li \n",nv);
    
    
    
    igraph_vector_init(&eids, 0);
    igraph_get_edgelist(&graph, &eids, false);
    
    
    unordered_map<Edge, double>EIndicator;
    
    for(long int i=0;i<ne;i++){
        long int from=VECTOR(eids)[2*i];
        long int to=VECTOR(eids)[2*i+1];
        double p_e=VECTOR(pe)[i];
        
        if(from>to){
            swap(from,to);
        }
        // use this one as matrix
        EIndicator[Edge(from,to)]=p_e;
    }
    
    
    
    
    
    
    /*randomized geneation candiates*/
    for(int tn=0;tn<attempt;tn++){
        
        long int ce=lround(c*ne);
        igraph_vector_t EC,ue,tpe; // just for temporary uncertain graph
        long int k=0;
        long count=ne;
        
        
        
        igraph_real_t ueSum,sigma_Sum;
        ueSum=0;
        sigma_Sum=0;
        
        random_device rd;
        
        std::mt19937 gen(rd());
        
        
        igraph_vector_init(&EC, 2*ce);
        
        unordered_map<Edge, double> ECIndicator;
        for(long int i=0;i<ne;i++){
            long int from=VECTOR(eids)[2*i];
            long int to=VECTOR(eids)[2*i+1];
            
            if(from>to){
                swap(from,to);
            }
            ECIndicator[Edge(from,to)]=1;
        }
        //
        
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
            
            double pickP=unDist(gen);
            
            
            double p_uv=0;
            
            auto search=EIndicator.find(Edge(u,v));
            if(search!=EIndicator.end()){
                p_uv=search->second;
            }
            
            
            int existEC=0;
            auto ecSearch=ECIndicator.find(Edge(u,v));
            // exist or not
            if(ecSearch!=ECIndicator.end()){
                existEC=ecSearch->second;
            }
            
            if(p_uv>pickP){
                
                
                // remove edge
                if(existEC==1){
                    count-=1;
                    ecSearch->second=0;
                }else{
                }
            }else{
                // add one edge
                
                if(existEC==0 && p_uv==0){
                    igraph_vector_set(&EC, k++,u);
                    igraph_vector_set(&EC, k++,v);
                    
                    count+=1;
                    
                    if(ecSearch!=ECIndicator.end()){
                        ecSearch->second=1;
                    }else{
                        ECIndicator[Edge(u,v)]=1;
                    }
                }
                
            }
            
            
            if(count==ce){
                break;
            }
            
        }
        
        long int exist=0;
        
        for(long int i=0;i<ne;i++){
            long int from=VECTOR(eids)[2*i];
            long int to=VECTOR(eids)[2*i+1];
            if(from>to){
                swap(from,to);
            }
            
            int existEC=0;
            auto ecSearch=ECIndicator.find(Edge(from,to));
            
            if(ecSearch!=ECIndicator.end()){
                existEC=ecSearch->second;
            }
            
            
            if(existEC==1){
                igraph_vector_set(&EC, k++, from);
                igraph_vector_set(&EC, k++, to);
                exist+=1;
            }
            
        }
        
        
        printf("finish edge selection:%li edges \n", k);
        printf("select edge from existing edge :%li\n",exist);
        printf("size of EC : %li\n", igraph_vector_size(&EC));
        printf("cout:%li, ce:%li, k:%li \n",count,ce,k);
        
        
        
        if(k!=2*ce){
            throw std::exception();
        }
        
        
        
        igraph_vector_init(&ue,ce);
        igraph_vector_init(&tpe,ce);
        
        
        
        
        
        for(long int i=0;i<ce;i++){
            igraph_real_t eVal=0;
            eVal+=VECTOR(uv)[int(VECTOR(EC)[2*i])];
            eVal+=VECTOR(uv)[int(VECTOR(EC)[2*i+1])];
            eVal/=2;
            ueSum+=eVal;
            igraph_vector_set(&ue,i,eVal);
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
            
            
            long int from=VECTOR(EC)[2*i];
            long int to=VECTOR(EC)[2*i+1];
            
            
            
            
            if(from>to){
                swap(from,to);
            }
            
            //
            double old_p_e=0;
            
            auto search=EIndicator.find(Edge(from,to));
            if(search!=EIndicator.end()){
                old_p_e=search->second;
            }
            
            
            if(w<noise){
                re=unGenDist(sgen);
                unCount+=1;
                re*=(2*(0.5-old_p_e));
            }else{
                re=truncated_normal_ab_sample(0.0, sigma_e, 0.0, 1.0, seed);
                re*=(2*(0.5-old_p_e));
                
            }
            
            // generate the solution
            old_p_e+=re;
            
            if(old_p_e<0 || old_p_e>1){
                throw std::exception();
            }
            igraph_vector_set(&tpe,i,old_p_e);
            
            
            reSum+=re;
        }
        
        
        printf("finish inject uncertainty \n");
        
        printf("inject uncertainty %f \n", reSum/ce);
        
        printf("inject uncertainty over white noise %f \n", (double)unCount/ce);
        
        printf("sum of edge :%f \n", igraph_vector_sum(&tpe));
        
        printf("average of sigme %f \n",sigma_Sum/ce);
        
        printf("init edge and probs \n");
        UncertainGraph pGraph((igraph_integer_t)nv);
        pGraph.set_edges(&EC);
        pGraph.set_edges_probs(&tpe);
        
        igraph_vector_destroy(&ue);
        igraph_vector_destroy(&tpe);
        igraph_vector_destroy(&EC);
        ECIndicator.clear();
        
        // check this graph
        pGraph.graphCheck();
        
        double epsilon_G=pGraph.testAgaist(ak);
        
        if(epsilon_G<epsilonStart){
            pg=pGraph;
            epsilonStart=epsilon_G;
            cout<<"find better obfuscation one :"<<epsilon_G<<endl;
        }
        
        cout<<"finish "<<tn<<" attempt" <<endl;
        
        
    }
    
    
    
    //igraph_vector_destroy(&cuv);
    EIndicator.clear();
    igraph_vector_destroy(&degrees);
    igraph_vector_destroy(&uv);
    igraph_vector_destroy(&eids);
    lowNodes.clear();
    edgeShuffle.clear();
    nodeUNs.clear();
    
    
    *eps_res=epsilonStart;
    
    
    
    
    
    return pg;
    
}

UncertainGraph UncertainGraph::generateObfuscation(igraph_real_t sigma, igraph_real_t * eps_res, igraph_vector_t * ak){

    switch (option) {
        case randPert:
            return randomGenerateObfuscation(sigma, eps_res, ak);
            break;
            
        case greedPert:
            return greedyGenerateObfuscation(sigma, eps_res, ak);
            break;
    }
}

UncertainGraph UncertainGraph::obfuscation(igraph_vector_t *ak){
    
    igraph_real_t sigmaLow, sigmaUpper, final_sigma;
    igraph_real_t ep_res, tEpsilion;
    
    UncertainGraph tGraph((igraph_integer_t)nv);
    
    sigmaLow=0;
    sigmaUpper=0.25; // this is one can be claimed by coders for speeding up the execution
    
    
    while(true){
        cout<<"random search for sigmaUpper="<<sigmaUpper<<endl;
        ep_res=1;
        UncertainGraph pGraph=generateObfuscation(sigmaUpper, &ep_res, ak);
        
        if(ep_res>=epsilon){
            sigmaLow=sigmaUpper;
            sigmaUpper*=2;
        }else{
            tGraph=pGraph;
            tEpsilion=ep_res;
            break;
        }
        
    }
    // sigmaupper pass , sigmaUpper
    
    final_sigma=sigmaUpper;
    
    while((sigmaUpper-sigmaLow)>0.00001 ){
        
        igraph_real_t sigma_mid=(sigmaUpper+sigmaLow)/2;
        cout<<"random search for sigmaUpper="<<sigma_mid<<endl;
        ep_res=1;
        UncertainGraph pGraph=generateObfuscation(sigma_mid, &ep_res, ak);
        if(ep_res>=epsilon){
            sigmaLow=sigma_mid;
        }else{
            tGraph=pGraph;
            if(sigma_mid<final_sigma){
                final_sigma=sigma_mid;
            }
            sigmaUpper=sigma_mid;
            tEpsilion=ep_res;
        }
        
        
        
        
        
        
    }
    
    cout<<"inject perturbation sigma :"<<final_sigma<<endl;
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
        
       ;
        
        
        //cout<<"clear i:"<<i<<endl;
        if(e_prob>=x){
           // cout<<"set i:"<<i<<"to 1"<<endl;
            igraph_vector_set(indicator,i,1);
            sampleEdges.push_back(VECTOR(edges)[2*i]);
            sampleEdges.push_back(VECTOR(edges)[2*i+1]);
        }else{
           
            igraph_vector_set(indicator,i,0);
        }
    }
    
    
    if(sampleEdges.size()!=2*igraph_vector_sum(indicator)){
        print_vector(indicator,"the stupid indicator \n");
        cout<<"sample Edge.size()"<<sampleEdges.size()<<endl;
        throw std::exception();
    }
    
   
   // cout<<"edge sum in sample Graph:"<<igraph_vector_sum(indicator)<<endl;
    
    
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


// just for debugging

void UncertainGraph::graphCheck(){
    igraph_vector_t deg;
    
    igraph_vector_init(&deg, nv);
    getDegrees(false, &deg);
    
    int edgeNum=igraph_vector_sum(&deg)/2;
    if(edgeNum!=ne){
        throw std::exception();
    }
}


