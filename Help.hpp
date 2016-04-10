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

#endif /* Help_hpp */
