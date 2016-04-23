//
//  Help.hpp
//  testGraphCplus
//
//  Created by dongqingxiao on 3/19/16.
//  Copyright © 2016 dongqingxiao. All rights reserved.
//

#ifndef Help_hpp
#define Help_hpp

#include <stdio.h>
#include <igraph.h>
#include <string>
#include <iostream>
#include <unordered_map>
#include <vector>
#include <map>
#include <boost/numeric/ublas/matrix.hpp>
using namespace std;

// just for debugging 
void print_vector(igraph_vector_t *v, std::string exp);
void print_matrix(const igraph_matrix_t *m,std::string exp);
void print_sp_matrix(const igraph_spmatrix_t *m, std::string exp);

void print_real_t_array(igraph_real_t * v, long int length,std::string exp);


class Node_UN{

    public :
    
    long int nodeID;
    double node_comm;
    long int nodeDegree;
    
    Node_UN(long int id, long int degree, double comm): nodeID(id), node_comm(comm),nodeDegree(degree){};

   
    
    bool operator < (const Node_UN & other) const {
        if(node_comm > other.node_comm){return true;}
        
        if(node_comm==other.node_comm){
            if(nodeDegree > other.nodeDegree){
                return true;
            }else{
                if(nodeDegree==other.nodeDegree){
                    return nodeID >other.nodeID;
                }
            }
            
        }
        
        return false;
    };
    
   friend  std::ostream & operator<<(std::ostream & strm, const Node_UN &a){
        
        strm << "ID:"<<a.nodeID <<"Degree:"<<a.nodeDegree<<"Un:"<<a.node_comm;
       return strm;

    }
   

};


class NodeRep{
public:
    long int nodeID;
    long int repID;
    NodeRep(long int id, long int rep): nodeID(id),repID(rep){};
    
    bool operator <(const NodeRep & other) const {
        if(repID<other.repID){return true;}
        if(repID==other.repID){
            return nodeID<other.nodeID;
        }
        return false;
    };
};



class ProbEdge{
public:
    long int from;
    long int to;
    long int rep;
    double pe;
    
    ProbEdge(long int f, long int t, long int r, double p): from(f),to(t),rep(r),pe(p){};
    
    bool operator<(const ProbEdge & other) const{
        if(rep<other.rep){return true;}
        if(rep>other.rep){return false;}
        
        if(from<other.from){return true;}
        if(from>other.from){return false;}
        
        return to<other.to;
    };
};

typedef std::pair<long int, long int>  Edge;

namespace std {
    template <>
    struct hash<Edge> {
        size_t operator()(const Edge &m) const {
            // hash method here.
            return hash<int>()(m.first) ^ hash<int> ()(m.second);
        }
    };
}








//vector file I.O
void write_vector_file(igraph_vector_t *res, string filePath);
void write_boostMatrix_file(boost::numeric::ublas::matrix<int> &m, string filepath);

void init_vector_file(igraph_vector_t *res, string filePath);
void init_boostMatrix_file(boost::numeric::ublas::matrix<int> &m, string filepath); // load everything into memory guess how much it needed..
void vector_statstic(igraph_vector_t *input);

#endif /* Help_hpp */
