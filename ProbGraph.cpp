//
//  ProbGraph.cpp
//  testGraphCplus
//
//  Created by dongqingxiao on 3/19/16.
//  Copyright © 2016 dongqingxiao. All rights reserved.
//

#include "ProbGraph.hpp"
#include "Help.hpp"
#include "DDCal.hpp"
#include "truncated_normal.hpp"
#include <math.h>
#include <limits.h>
#include <random>
#include <boost/math/distributions/normal.hpp>
#include <map>
using boost::math::normal;




ProbGraph::ProbGraph(igraph_integer_t n){
    igraph_empty(&graph, n, IGRAPH_DIRECTED);// create empty graph
    nv=n;
}

long int ProbGraph::getNV(void){
    return nv;
}



ProbGraph::ProbGraph(const ProbGraph & obj){
    printf("copy constructor\n");
    igraph_copy(&graph, &obj.graph);
    nv=obj.nv;
    ne=obj.ne;
}


void ProbGraph::set_edges(const igraph_vector_t *edges){
    igraph_add_edges(&graph, edges, 0);
    ne=igraph_vector_size(edges);
    ne/=2;
}

void ProbGraph::set_edges_prob(const igraph_vector_t *probs){
    igraph_cattribute_EAN_setv(&graph, "prob", probs);
}



void ProbGraph::maxDegrees(igraph_vector_t *deg){
    igraph_degree(&graph, deg, igraph_vss_all(), IGRAPH_ALL,IGRAPH_LOOPS);
}

igraph_real_t ProbGraph::maxDegree(void){
    igraph_vector_t degs;
    igraph_vector_init(&degs, 0);
    maxDegrees(&degs);
    return igraph_vector_max(&degs);
}


void ProbGraph::ex_uncertain_entropyReport(igraph_vector_t * entropyReport,igraph_real_t maxDegree){
    igraph_vector_t degs, edgeProbs;
    igraph_matrix_t deMatrix;
    
    igraph_real_t degree;
    
    
    igraph_vector_init(&degs, nv);
    igraph_degree(&graph, &degs, igraph_vss_all(), IGRAPH_ALL,IGRAPH_LOOPS);
    igraph_vector_init(&edgeProbs, ne);
    
    igraph_matrix_init(&deMatrix, nv, maxDegree+1);
    igraph_matrix_null(&deMatrix);
    
    
    
    EANV(&graph,"prob",&edgeProbs);
    // iterate over all vertices
    for(long int  i=0;i<nv;i++){
        igraph_vector_t prob_v,res_v, eids;
        
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
        }else{
            // set the prob equals=0
            igraph_vector_set(&res_v,0,1);
        }
        
        igraph_matrix_set_row(&deMatrix,&res_v,i);
        
        
        igraph_vector_destroy(&prob_v);
        igraph_vector_destroy(&eids);
        igraph_vector_destroy(&res_v);
    }
    
    // iterator over the range of degree
    for(long int i=0;i<maxDegree+1;i++){
        igraph_vector_t s_column;
        igraph_real_t sc_sum=0.0;
        igraph_real_t entropy=0.0;
        igraph_vector_init(&s_column,nv);
        igraph_matrix_get_col(&deMatrix, &s_column, i);
        sc_sum=igraph_vector_sum(&s_column);
        if(sc_sum<=nv*numeric_limits<double>::epsilon()){
            
            continue;
        }
        
      
        
        igraph_vector_scale(&s_column, 1.0/sc_sum);
        
        
        
        
        for(int j=0;j<nv;j++){
            double p= VECTOR(s_column)[j];
            double val=cal_entropy(p);
            entropy+=val;
        }
        
    
        igraph_vector_set(entropyReport,i,entropy);
        igraph_vector_destroy(&s_column);
    }
    
    
    
    
    igraph_vector_destroy(&degs);
    igraph_vector_destroy(&edgeProbs);
    igraph_matrix_destroy(&deMatrix);
 
    
}



void ProbGraph::uncertain_entropyReport(igraph_vector_t * entropyReport,igraph_real_t maxDegree){
    
    igraph_vector_t degs, edgeProbs, s,s_entropy;
    
    igraph_real_t degree;
    
    
    igraph_vector_init(&degs, nv);
    igraph_degree(&graph, &degs, igraph_vss_all(), IGRAPH_ALL,IGRAPH_LOOPS);
    igraph_vector_init(&edgeProbs, ne);
    
    
    igraph_vector_init(&s,maxDegree+1);
    igraph_vector_init(&s_entropy,maxDegree+1);
   
    
    EANV(&graph,"prob",&edgeProbs);
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
        
//        if(i%10000==0){
//            cout<<"at least move"<<i<<endl;
//        }
   

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

// rewrite the interface
void ProbGraph::certain_entropyReport(igraph_vector_t * entropyReport){
    igraph_vector_t degs;
    igraph_real_t maxDegree;
    igraph_vector_init(&degs, 0);
    igraph_degree(&graph, &degs, igraph_vss_all(), IGRAPH_ALL,IGRAPH_LOOPS);
    
    maxDegree=igraph_vector_max(&degs);
    for(long int i=0;i<nv;i++){
        long int degree=VECTOR(degs)[i];
        igraph_vector_set(entropyReport,degree,VECTOR(*entropyReport)[degree]+1);
    }
    
    
    
    for(long int i=0;i<maxDegree+1;i++){
        long int val=VECTOR(*entropyReport)[i];
        if(val!=0) {
            igraph_vector_set(entropyReport,i,log2(val));
        }else{
            igraph_vector_set(entropyReport,i,-1);
        }
    }
    
    igraph_vector_destroy(&degs);
}





void ProbGraph::entropyReport(igraph_vector_t * entropyReport,igraph_real_t maxDegree){
    igraph_bool_t uncertain;
    
    uncertain=igraph_cattribute_has_attr(&graph,IGRAPH_ATTRIBUTE_EDGE,"prob");
    
    if(uncertain){
        uncertain_entropyReport(entropyReport,maxDegree);
    }else{
        certain_entropyReport(entropyReport);
    }
    
}


/*
 * iterate over vertex and cal its prob with degree=val 
 * val in [0, maxDegree]
 */
void ProbGraph::prob_vertex(igraph_vector_t * prob_v,igraph_integer_t vid){
    igraph_vector_t eids;
    igraph_real_t  degree;
    igraph_vector_t edgeProbs;
    
    
    igraph_vector_init(&eids,nv);
    igraph_vector_init(&edgeProbs,ne);
    EANV(&graph,"prob",&edgeProbs);
    igraph_incident(&graph, &eids, vid, IGRAPH_ALL);
    degree=igraph_vector_size(&eids);

    
    for(long int i=0;i<degree;i++){
        long int eid=VECTOR(eids)[i];
        igraph_vector_set(prob_v, i, VECTOR(edgeProbs)[eid]);
    }

    igraph_vector_destroy(&eids);
    igraph_vector_destroy(&edgeProbs);
}


igraph_real_t ProbGraph::testAgainst(igraph_vector_t * ak){
    
    igraph_real_t eResult=0;
    igraph_real_t theshold;
    
    igraph_real_t maxDegree;
    igraph_vector_t ePResult;
    
    maxDegree=igraph_vector_max(ak);
    igraph_vector_init(&ePResult, maxDegree+1);
    printf("init entropy report \n");
    entropyReport(&ePResult,maxDegree);
    printf("finsh init entropy Report \n");
    
    theshold=log2(k-2);
    for(long int i=0;i<nv;i++){
        igraph_real_t val=VECTOR(*ak)[i];
        double diff=theshold-VECTOR(ePResult)[(long int) val];
        if(diff>0.01){
            eResult+=1;
            printf("less deg=%f, diff=%f \n",val, diff);
        }
    }
    
    printf("lessAn:%f \n",eResult);
    igraph_vector_destroy(&ePResult);
    return eResult/nv;
}

igraph_real_t testBug(igraph_vector_t * degrees,long int nv){
    
    igraph_real_t maxDegree=igraph_vector_max(degrees);
    
    printf("maxDegre:%f \n",maxDegree);
    
    igraph_vector_t counter;
    igraph_real_t lessAn;

    
    lessAn=0;
    
    igraph_vector_init(&counter,maxDegree+1);
    
    for(long int i=0;i<nv;i++){
        long int degree=VECTOR(*degrees)[i];
        igraph_vector_set(&counter, degree, VECTOR(counter)[degree]+1);
    }
    
//    printf("the samllest number of nodes with same degree :%li \n", VECTOR(counter)[(long int)maxDegree]);
    
//    printf("deg=0 %f \n", VECTOR(counter)[0]);
    for(long int i=0;i<maxDegree+1;i++){
        if(VECTOR(counter)[i]<k){
            lessAn+=VECTOR(counter)[i];
            if(VECTOR(counter)[i]!=0){
               // printf("deg=%li \n",i);
            }
        }
    }
    
    printf("lessAn count:%f \n",lessAn);
    printf("epislon: %f \n",lessAn/nv);
    return lessAn/nv;
}


/*
 *self Test
 */
igraph_real_t ProbGraph::selfTest(void){
    igraph_vector_t expectDeg;
    igraph_vector_init(&expectDeg,nv);
    igraph_degree(&graph,&expectDeg,igraph_vss_all(), IGRAPH_ALL, IGRAPH_LOOPS);
    return testBug(&expectDeg, nv);
}


void ProbGraph::sigmaUniquess(igraph_vector_t * uv, igraph_vector_t ak, igraph_real_t maxDegree, igraph_real_t sigma){
    igraph_vector_t s_com,freqs;
    
    normal ns(0,sigma);

    igraph_vector_init(&freqs,maxDegree+1);
    igraph_vector_init(&s_com,maxDegree+1);
    
    for(long int i=0;i<nv;i++){
        long int degree=VECTOR(ak)[i];
        igraph_vector_set(&freqs, degree, VECTOR(freqs)[degree]+1);
    }
    
    
    printf("deg=0: %f\n", VECTOR(freqs)[0]);
//    printf("deg=135: %f\n", VECTOR(freqs)[135]);
    
    

    // when sigma is very small
    // it approximate count the number of nodes with the same degree
    for(long int i=0;i<maxDegree+1;i++){
        igraph_real_t val=0.0;
        // iterate over all nodes
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



ProbGraph ProbGraph::generateObfuscation(igraph_real_t sigma, igraph_real_t * eps_res){
    
    ProbGraph pg((igraph_integer_t)nv);
    
    
    igraph_real_t epsilonStart,maxDegree;
    igraph_vector_t degrees,uv,eids;
    
    igraph_matrix_t  EIndicator;
    vector<double> lowNodes;
    vector<int> edgeShuffle;
    vector<Node_UN> nodeUNs;
    
    
    
    epsilonStart=1;
    igraph_vector_init(&degrees,0);
    igraph_degree(&graph,&degrees,igraph_vss_all(), IGRAPH_ALL,IGRAPH_LOOPS);
    maxDegree=igraph_vector_max(&degrees);
    igraph_vector_init(&uv,nv);
    
    
    printf("start the computation of the sigma-uniqueness of all V \n");
    sigmaUniquess(&uv, degrees, maxDegree, sigma);
    printf("finish the computation \n");

    
    
    int skip_count=(int) lround(epsilon*nv/2)+1;


    for(long int i=0;i<nv;i++){
        nodeUNs.push_back(Node_UN(i,VECTOR(degrees)[i],VECTOR(uv)[i]));
        lowNodes.push_back(VECTOR(uv)[i]);
    }
    
    
    sort(nodeUNs.begin(),nodeUNs.end());
    
    for(int i=0;i<skip_count;i++){
        Node_UN  nu=nodeUNs[i];
       // cout<<i<<"node information"<<nu<<endl;
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
            
            
//            if(re==0){
//                printf("**I hate the world \n");
//            }
            
//            if(re<0.1){
//                re*=5;
//            }
//
            
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
        
        ProbGraph pGraph((igraph_integer_t)nv);
        pGraph.set_edges(&EC);
        pGraph.set_edges_prob(&pe);
        printf("init edge and probs \n");
        
        igraph_vector_destroy(&ue);
        igraph_vector_destroy(&pe);
        igraph_vector_destroy(&EC);
        igraph_matrix_destroy(&ECIndicator);
        
        
        
        pGraph.uncertainGraphStastic();
        
        
        
        igraph_real_t epsilon_G=pGraph.testAgainst(&degrees);
        printf("epsilon_G:%f \n",epsilon_G);
        
        
        if(epsilon_G<epsilonStart){
            pg=pGraph;
            epsilonStart=epsilon_G;
            printf("find better obfuscation one \n");
        }
            
        
        
        
        
        printf("end one iteration\n");
        
    }

    
    
    //igraph_vector_destroy(&cuv);
    igraph_vector_destroy(&degrees);
    igraph_vector_destroy(&uv);
    igraph_matrix_destroy(&EIndicator);
    
  
    *eps_res=epsilonStart;
    
    return pg;
}





ProbGraph ProbGraph::obfuscation(){
    igraph_real_t sigmaLow, sigmaUpper;
    igraph_real_t ep_res, tEpsilion;
    
    ProbGraph tGraph((igraph_integer_t)nv);
    
    sigmaLow=0;
    sigmaUpper=1;
    
    
    while(true){
        printf("search for sigmaUpper:%f \n", sigmaUpper);
        ep_res=1;
        ProbGraph pGraph=generateObfuscation(sigmaUpper, &ep_res);
        
        if(ep_res>=epsilon){
            sigmaUpper*=2;
        }else{
            tGraph=pGraph;
            tEpsilion=ep_res;
            break;
        }
        
    }
    
    while(true){
        printf("2 search for sigmaUpper:%f \n", sigmaUpper);
        igraph_real_t sigma=(sigmaUpper+sigmaLow)/2;
        ProbGraph pGraph=generateObfuscation(sigmaUpper, &ep_res);
        if(ep_res>=epsilon){
            sigmaLow=sigma;
            printf("up \n");
        }else{
            printf("down \n");
            tGraph=pGraph;
            sigmaUpper=sigma;
            tEpsilion=ep_res;
        }
        
        if((sigmaUpper-sigmaLow)<0.0001){
            break;
        }
        
       
        
       
    
    }
    
    
    
    printf("inject perturbation sigma:%f \n",sigmaUpper);
    printf("tolerance level epislion:%f \n",tEpsilion);
    
    return tGraph;
}


void ProbGraph::certainGraphStatstic(void){
    
    igraph_vector_t degrees;
    igraph_real_t nne,ad,md;
    
    igraph_vector_init(&degrees, nv);
    

    igraph_degree(&graph,&degrees,igraph_vss_all(), IGRAPH_ALL, true);
    
    
    
    nne=igraph_vector_sum(&degrees);
    nne/=2;
    
    md=igraph_vector_max(&degrees);
    
    
    ad=2*nne/nv;
    
    printf("certain stastic \n");
    printf("max degree gained from node :%li \n",igraph_vector_which_max(&degrees));
    
    printf("the number of vertice: %li \n",nv);
    printf("NE\tAD\tMD\n");
    printf("%2f\t%2f\t%2f\n", nne,ad,md);
    
    igraph_vector_destroy(&degrees);

}
/** 
 * AD: maximal degree of certain graph 
 *
 */
void ProbGraph::certain_metrics(igraph_vector_t * res, igraph_real_t p){
    igraph_vector_t degrees;
    igraph_real_t md;
    igraph_vector_init(&degrees, nv);
    igraph_degree(&graph,&degrees,igraph_vss_all(), IGRAPH_ALL, true);
    md=igraph_vector_max(&degrees);
    // set the maximal degree value;
    igraph_vector_set(res,0,VECTOR(*res)[0]+md*p);
    igraph_vector_destroy(&degrees);
}


void ProbGraph::uncertainGraphStastic(void){
    igraph_vector_t  expected_Degree,edge_probs, eids;
    igraph_vector_t sampleRes;
    
    igraph_real_t nne,ad,md,epMean;
    
    igraph_vector_init(&expected_Degree,nv);
    igraph_vector_init(&edge_probs,ne);
    igraph_vector_init(&eids,2*ne);
    igraph_vector_init(&sampleRes,3); // this case we only cal AD and ..
    
    
    EANV(&graph,"prob",&edge_probs);
    igraph_get_edgelist(&graph, &eids, false);

    printf("call uncertain Graph statstic \n");
    
    for(long int i=0;i<ne;i++){
        long int from=VECTOR(eids)[2*i];
        long int to=VECTOR(eids)[2*i+1];
        igraph_real_t p=VECTOR(edge_probs)[i];
        igraph_vector_set(&expected_Degree, from,VECTOR(expected_Degree)[from]+p);
        igraph_vector_set(&expected_Degree, to,VECTOR(expected_Degree)[to]+p);
    }
    
    nne=igraph_vector_sum(&expected_Degree);
    nne/=2;
    
    ad=2*nne/nv;

    
    epMean=igraph_vector_sum(&edge_probs)/igraph_vector_size(&edge_probs);

    
    double p=1.0/sampleNum;
    
    for(int i=0;i<sampleNum;i++){
        ProbGraph sg=sampleGraph();
        sg.certain_metrics(&sampleRes, p);
    }
    
    md=VECTOR(sampleRes)[0];
    
    printf("uncentrain stastic \n");
    printf("nv:%ld \n",nv);
    printf("edge probability(mean): %2f \n", epMean);
    printf("NE\tAD\tMD\n");
    printf("%2f\t%2f\t%2f\n", nne,ad,md);

    
    igraph_vector_destroy(&expected_Degree);
    igraph_vector_destroy(&edge_probs);
    igraph_vector_destroy(&eids);
    igraph_vector_destroy(&sampleRes);
    
    
}


ProbGraph init_from_Adj_File(string path){
    ifstream graphFile(path);
    string line;
    string item;
    vector<double> v;
    
    igraph_vector_t v_graph;
    long int cv=0;
    if(graphFile.is_open()){
        while(getline(graphFile,line)){
            stringstream ss(line);
            //printf("cv:%ld \n",cv);
            while(getline(ss,item,' ')){
                long int temp=cv;
                long int nv=stol(item);
            
                //printf("temp:%ld,nv:%ld \n",temp,nv);
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
    
    ProbGraph pg((igraph_real_t)cv);
    pg.set_edges(&v_graph);
    printf("initiation graph structure from file with edge\n");
    
    return pg;
    
}


ProbGraph init_from_OB_File(string path,long int nv){
    
    ifstream graphFile(path);
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
    ProbGraph pg((igraph_real_t) nv);
    pg.set_edges(&v_graph);
    
    double * parray=edge_prob.data();
    printf("edge prob size:%li \n",edge_prob.size());
    
    igraph_vector_view(&e_probs, parray, edge_prob.size());

   // print_vector(&e_probs, "e_probs");
    pg.set_edges_prob(&e_probs);
    
    printf("initiation graph structure from file with edge prob \n");
    
    return pg;

}




ProbGraph uncertain_init_from_File(string path){
    ifstream graphFile(path);
    string line;
    string item;
    vector<double> v;
    igraph_vector_t v_graph,e_probs;
    vector<double> edge_prob;

    cout<<"openfile:"<<path<<endl;
    cout<<"delimiter:"<<delimiter<<endl;
    
    
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
        cout<<"can't open file"<<endl;
    }
    
    double * array=v.data();
    igraph_vector_view(&v_graph, array, v.size());
    igraph_real_t nv=igraph_vector_max(&v_graph);
    
    ProbGraph pg(nv+1);
    pg.set_edges(&v_graph);
    
    double * parray=edge_prob.data();
    printf("edge prob size:%li \n",edge_prob.size());

    igraph_vector_view(&e_probs, parray, edge_prob.size());
    
    //print_vector(&e_probs, "e_probs");
    pg.set_edges_prob(&e_probs);
    
    printf("initiation graph structure from file with edge prob \n");
    
    return pg;
}


ProbGraph certain_init_from_File(string path){
    
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
    }
    
    double * array=v.data();
    igraph_vector_view(&v_graph, array, v.size());
    igraph_real_t nv=igraph_vector_max(&v_graph);
    ProbGraph pg(nv+1);
    pg.set_edges(&v_graph);
    
    printf("initiation graph structure from file \n");
    
    return pg;
}



ProbGraph init_from_File(string path){
    
    if(uncertain){
        return uncertain_init_from_File(path);
    }else{
        return certain_init_from_File(path);
    }
}




void ProbGraph::certain_reliablityReport(igraph_spmatrix_t *res, igraph_real_t p){
    

    igraph_vector_t nv_commponet;
    igraph_vector_init(&nv_commponet,nv);
    int k=0;
    
 
    
    for(long int i=0;i<nv;i++){
        
        if(i%1==0){
            cout<<"at least move 1"<<endl;
        }
        if(VECTOR(nv_commponet)[i]==0){
        
            igraph_vector_t cres;
            //igraph_vector_t cres_t;
            igraph_vector_init(&cres,nv);
        
            igraph_subcomponent(&graph, &cres, i, IGRAPH_ALL);
            igraph_real_t size=igraph_vector_size(&cres);
            k+=1;
            igraph_vector_set(&nv_commponet,i,k);
            
            igraph_vector_sort(&cres);
            
            for(long int j=0;j<size;j++){
                long int n_j=VECTOR(cres)[j];
                igraph_vector_set(&nv_commponet,n_j,k);
                
                igraph_spmatrix_set(res, i, n_j, p);
               
                
                for(long int k=j+1;k<size;k++){
                    long int n_k=VECTOR(cres)[k];
                    igraph_spmatrix_set(res, n_j, n_k, p);
                    
                }
            }
            
            
            
            
            
        
   
            igraph_vector_destroy(&cres);
           
        }
    }
    
    
    igraph_vector_destroy(&nv_commponet);
  
    
}


void ProbGraph::certain_reliablityReport_file(string filePath, igraph_real_t p){
    igraph_vector_t nv_commponet,cstart,crecord;
    
    igraph_vector_init(&nv_commponet,nv+1);
    igraph_vector_init(&cstart,nv);
    igraph_vector_init(&crecord,nv);
    
    int k=0;
    int c_start=0;
    
    
    for(int i=0;i<nv;i++){
        //        if(i%10000==0){
        //            cout<<"move"<<i<<endl;
        //        }
        
        if(VECTOR(nv_commponet)[i]==0){
            
            
            igraph_vector_t cres;
            igraph_vector_init(&cres,nv);
            
            int code=igraph_subcomponent(&graph, &cres, i, IGRAPH_ALL);
            
            if(code==IGRAPH_EINVMODE||code==IGRAPH_EINVVID|| code==IGRAPH_ENOMEM ){
                printf("serach error\n");
            }
            
            
            long int size=igraph_vector_size(&cres);
            
            k+=1;
            
            igraph_vector_set(&cstart,k,c_start);
            
            for(long int j=0;j<size;j++){
                
                long int n_j=VECTOR(cres)[j];
                igraph_vector_set(&crecord,c_start++,n_j);
                igraph_vector_set(&nv_commponet,n_j,k);
            }
            
            
            igraph_vector_destroy(&cres);
            
        }
        
    }
    
    igraph_vector_set(&cstart,k+1,nv);
    
    cout<<"sub component done"<<endl;
    
    
    ofstream myfile (filePath);
    
    if(!myfile.is_open()){
        cout<<"no open file"<<endl;
        cout<<"file path:"<<filePath<<endl;
    }else{
        cout<<"file path:"<<filePath<<endl;
    }
    
    for(long int i=0;i<nv;i++){
        long int c=VECTOR(nv_commponet)[i];
        long int start=VECTOR(cstart)[c];
        long int end=VECTOR(cstart)[c+1];
        
        for(long int j=start;j<end;j++){
            if(j!=end-1){
                myfile<<VECTOR(crecord)[j]<<"\t"<<p<<"\t";
            }else{
                myfile<<VECTOR(crecord)[j]<<"\t"<<p;
            }
        }
        
        myfile<<endl;
    }
    
    myfile.close();
    cout<<"writing finish"<<endl;
    
    igraph_vector_destroy(&nv_commponet);
    igraph_vector_destroy(&cstart);
    igraph_vector_destroy(&crecord);
    
    
}



void ProbGraph::uncertain_reliablityReport_ex(igraph_vector_t *res){
    
    double p=1.0/sampleNum;
    
    for(int i=0;i<sampleNum;i++){
        ProbGraph sg=sampleGraph();
        cout<<i<<"iteartion"<<endl;
       // sg.certain_reliablityReport_ex(res, p);
        sg.certain_reliablityReport_ex_inc(res, p);
    }
}

void ProbGraph::uncertain_reliablityReport_ex_store(igraph_vector_t *res, string filepath, string label){
    // init res
    uncertain_reliablityReport_ex(res);
    //
    string path=filepath+"/"+label+"_rel.txt";
    ofstream myfile(path);
    if(!myfile.is_open()){cout<<"can not open file"<<path<<endl;}
    for(long int i=0;i<nv;i++){
        myfile<<VECTOR(*res)[i]<<endl;
    }
    myfile.close();
    cout<<"done"<<endl;
}

void ProbGraph::certain_reliablityReport_ex_store(igraph_vector_t *res, string filepath, string label){
    
    certain_reliablityReport_ex(res, 1);
    string path=filepath+"/"+label+"_rel.txt";
    ofstream myfile(path);
    if(!myfile.is_open()){cout<<"can not open file"<<path<<endl;}
    for(long int i=0;i<nv;i++){
        myfile<<VECTOR(*res)[i]<<endl;
    }
    myfile.close();
    cout<<"done"<<endl;
}

/**
 * compute connected component incrementally
 */
void ProbGraph::certain_reliablityReport_ex_inc(igraph_vector_t *res,igraph_real_t p){
    
    
    igraph_vector_t edges;
    igraph_vector_init(&edges,2*ne);
    igraph_vector_t v_component; // record its represetnative
    igraph_vector_t c_component; // record the number of count
    
    
    
    
    vector<long int> rank (nv);
    vector<long int> parent(nv);
    igraph_get_edgelist(&graph, &edges, false);
    
    
    for(int i=0;i<nv;i++){
        rank[i]=i;
        parent[i]=i;
    }
    
    
    igraph_vector_init(&v_component,nv);
    igraph_vector_init(&c_component,nv);
    
    
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
    
    // for debugging
    
    
    //record val
    for(long int i=0;i<nv;i++){
        double inc=VECTOR(c_component)[(int)(VECTOR(v_component)[i])]-1;
        inc*=p;
        igraph_vector_set(res,i,VECTOR(*res)[i]+inc);
    }

    
    // done 
}




/** cal reliablity report via call
 * igraph_subcomponent
 */
void ProbGraph::certain_reliablityReport_ex(igraph_vector_t *res, igraph_real_t p){
    igraph_vector_t component_count;
    
    igraph_vector_init(&component_count,nv);
    
    cout<<"start subcomponent finding "<<endl;
    for(int i=0;i<nv;i++){
        
        if(i%10000==0){
            cout<<"move"<<i<<endl;
        }
        
        if(VECTOR(component_count)[i]==0){
            
            
            igraph_vector_t cres;
            igraph_vector_init(&cres,nv);
            
            int code=igraph_subcomponent(&graph, &cres, i, IGRAPH_ALL);
            
            if(code==IGRAPH_EINVMODE||code==IGRAPH_EINVVID|| code==IGRAPH_ENOMEM ){
                printf("serach error\n");
            }
            
            
            long int size=igraph_vector_size(&cres);
            
            
            
            for(long int j=0;j<size;j++){
                
                long int n_j=VECTOR(cres)[j];
                igraph_vector_set(&component_count, n_j, size-1);
                
            }
            
            
            igraph_vector_destroy(&cres);
            
        }
        
    }

    
    for(long int i=0;i<nv;i++){
        igraph_vector_set(res,i,VECTOR(*res)[i]+VECTOR(component_count)[i]*p);
    }
    
    cout<<"one iteration done"<<endl;
    
}


void ProbGraph::certain_reliablityReport_file_ex(string filePath, igraph_real_t p){
    
    igraph_vector_t component_count;
    
    igraph_vector_init(&component_count,nv);

    cout<<"start subcomponent finding "<<endl;
    for(int i=0;i<nv;i++){
        
        if(i%10000==0){
            cout<<"move"<<i<<endl;
        }
        
        if(VECTOR(component_count)[i]==0){
            
            
            igraph_vector_t cres;
            igraph_vector_init(&cres,nv);
            
            int code=igraph_subcomponent(&graph, &cres, i, IGRAPH_ALL);
            
            if(code==IGRAPH_EINVMODE||code==IGRAPH_EINVVID|| code==IGRAPH_ENOMEM ){
                printf("serach error\n");
            }
            
            
            long int size=igraph_vector_size(&cres);
        
          
            
            for(long int j=0;j<size;j++){
                
                long int n_j=VECTOR(cres)[j];
                igraph_vector_set(&component_count, n_j, size-1);
                
            }
            
            
            igraph_vector_destroy(&cres);
            
        }
        
    }
    
   
    
    cout<<"sub component done"<<endl;
    
    
    ofstream myfile (filePath);
    
    if(!myfile.is_open()){
        cout<<"no open file"<<endl;
        cout<<"file path:"<<filePath<<endl;
    }else{
        cout<<"file path:"<<filePath<<endl;
    }
    
    for(long int i=0;i<nv;i++){
        long int size=VECTOR(component_count)[i];
        myfile<<((double)size)*p<<endl;
    }
    
    myfile.close();
    
    
}


void ProbGraph::uncertain_reliablityReport(igraph_spmatrix_t *res){
    
    
    double p=1.0/sampleNum;
    
    for(int i=0;i<sampleNum;i++){
        ProbGraph sg=sampleGraph();
        sg.certain_reliablityReport(res,p);
    }
    
    
}





string filePathGen(string filePath,long int i){
    string tempath=filePath+"/"+"temp"+"_"+to_string(i)+".txt";
    return tempath;
}

string mergePathGen(string filePath, long int i){
    string tempath=filePath+"/"+"merge"+"_"+to_string(i)+".txt";
    return tempath;
}


string mergeFile_ex(string tempPath, string mergePath, long int i,string filepath, long int nv){
    if(i==0){
        return tempPath;
    }else{
        string mergePath=mergePathGen(filepath,i);
        ofstream nMerge(mergePath);
        ifstream tempfile(tempPath);
        ifstream mergefile(mergePath);
        string online;
        string seline;
        for(long int i=0;i<nv;i++){
            //hash map
            getline(mergefile,online);
            getline(tempfile,seline);
            double oldval=stod(online);
            double newVal=stod(seline);
            
            nMerge<<oldval+newVal<<endl;
        }
        
        return mergePath;
    }

}


// it doesn't work since the file would be too huge
// may be we can try using bit represetnation
string mergeFile(string tempPath, string mergePath, long int i,string filepath, long int nv){
    if(i==0){
        return tempPath;
    }else{
        string mergePath=mergePathGen(filepath,i);
        ofstream nMerge(mergePath);
        ifstream tempfile(tempPath);
        ifstream mergefile(mergePath);
        string online;
        string seline;
        for(long int i=0;i<nv;i++){
            //hash map
            map<long int, double> valMap;
            getline(tempfile,online);
            
            stringstream ss(online);
            string item;
            long int j=0;
            long int pos;
            double val;
            while(getline(ss,item,'\t')){
                if(j%2==0){
                    pos=stol(item);
                }else{
                    val=stod(item);
                    //
                    auto it=valMap.find(pos);
                    if(it==valMap.end()){
                        valMap.insert(pair<long int, double>(pos,val));
                    }else{
                        it->second=it->second+val;
                    }
                }
                j+=1;
            }
            
            getline(mergefile,online);
            
            while(getline(ss,item,'\t')){
                if(j%2==0){
                    pos=stol(item);
                }else{
                    val=stod(item);
                    //
                    auto it=valMap.find(pos);
                    if(it==valMap.end()){
                        valMap.insert(pair<long int, double>(pos,val));
                    }else{
                        it->second=it->second+val;
                    }
                }
                j+=1;
            }
            
            
            long int size=valMap.size();
            int  count=0;
            for(auto vit:valMap){
                if(count!=size-1){
                    nMerge<<vit.first<<"\t"<<vit.second<<"\t";
                }else{
                    nMerge<<vit.first<<"\t"<<vit.second;
                }
            }
            
            nMerge<<endl;
        }
        
        /**/
        
        
        
        return mergePath;
    }
    
}

void ProbGraph::uncertain_reliablityReport_file(string filePath){
    
    
    double p=1.0/sampleNum;
    
    string mergepath="empty";
    
    for(long int i=0;i<sampleNum;i++){
        ProbGraph sg=sampleGraph();
        string tempath=filePathGen(filePath,i);
        sg.certain_reliablityReport_file_ex(tempath, p);
        //mergepath=mergeFile_ex(tempath,mergepath,i, filePath,nv);
        cout<<"finish "<<i <<"th iteration"<<endl;
    }

}








void ProbGraph::reliablityReport(){
    


    igraph_spmatrix_t res;
    igraph_spmatrix_init(&res, nv, nv);
    
    if(uncertain){
        uncertain_reliablityReport(&res);
    }else{
        certain_reliablityReport(&res,1);
    }
    

    
}


ProbGraph ProbGraph::sampleGraph(void){
   
    ProbGraph g(nv);
    
    igraph_vector_t edges, edgeProbs,sEdges;
    vector<double> sampleEdges;
    random_device rand_dev;
    mt19937     generator(rand_dev());
    uniform_real_distribution<double> unDist(0.0,1.0);
    
    igraph_vector_init(&edges, 2*ne);
    igraph_vector_init(&edgeProbs,ne);
    igraph_get_edgelist(&graph, &edges, false);
    EANV(&graph,"prob",&edgeProbs);
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




void ProbGraph::print_uncertain_graph(string filepath){
    //
    igraph_vector_t eids;
    igraph_vector_t edgeProbs;
    
    igraph_vector_init(&eids,2*ne);
    igraph_vector_init(&edgeProbs, ne);
    igraph_get_edgelist(&graph, &eids, false);
    EANV(&graph,"prob",&edgeProbs);
    
    //
    
    ofstream myfile (filepath);
    if(!myfile.is_open()){
        cout<<"sorry can't open this file "<<endl;
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
}


void ProbGraph::print_certain_graph(string filepath){
    //
    igraph_vector_t eids;
    igraph_vector_init(&eids,2*ne);
    
    igraph_get_edgelist(&graph, &eids, false);
    
    ofstream myfile (filepath);
    if(!myfile.is_open()){
        cout<<"sorry can't open this file "<<endl;
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
}


void ProbGraph::print_graph(string filepath){
    
    igraph_bool_t c_uncertain;
    c_uncertain=igraph_cattribute_has_attr(&graph,IGRAPH_ATTRIBUTE_EDGE,"prob");
    
    if(c_uncertain){
        print_uncertain_graph(filepath);
    }else{
        print_certain_graph(filepath);
    }
}



void init_vector_from_file(string path, igraph_vector_t *res){
    
    ifstream myfile(path);
    long int size;
    string line;
    
    size=igraph_vector_size(res);
    
    if(!myfile.is_open()){
        cout<<"can not open "<<path<<endl;
        exit(0);
    }
    
    
    for(long int i=0;i<size;i++){
        getline(myfile,line);
        igraph_vector_set(res,i, stod(line));
    }

    // done 
    
}




