//
//  DDCal.hpp
//  testGraphCplus
//
//  Created by dongqingxiao on 3/21/16.
//  Copyright © 2016 dongqingxiao. All rights reserved.
//

#ifndef DDCal_hpp
#define DDCal_hpp

#include <stdio.h>
#include <igraph.h>

void DDCal(const igraph_vector_t * input, igraph_vector_t * res, igraph_real_t degree, igraph_real_t maxDegree);
igraph_real_t cal_entropy(igraph_real_t p);

double cal_distance(long int i, long int j);

double cal_distance_matrix(igraph_matrix_t * one, igraph_matrix_t * second, igraph_real_t nv);


double cal_distance_vector(igraph_vector_t * one, igraph_vector_t * second);

/**
 * one: reference 
 * second: test
 */
double cal_relative_error_vector(igraph_vector_t * one, igraph_vector_t * second);
double cal_mean_error_vector(igraph_vector_t * ref, igraph_vector_t * test);

double relative_error_metric(igraph_vector_t * one, igraph_vector_t * second, int i);



long int permuateCal(long int number); 


void is_sparese(igraph_vector_t * input, long int nv);

void reDistribute(igraph_vector_t * uv);
void reverseVector(igraph_vector_t * ruv);

double featureCombine(double uniq, double robust);

double standard_vaiance_vector(igraph_vector_t * ruv);



#endif /* DDCal_hpp */
