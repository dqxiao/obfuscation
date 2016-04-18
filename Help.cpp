//
//  Help.cpp
//  testGraphCplus
//
//  Created by dongqingxiao on 3/19/16.
//  Copyright © 2016 dongqingxiao. All rights reserved.
//

#include "Help.hpp"
#include <fstream>

void call_print_vector(igraph_vector_t *v, FILE *f) {
    long int i;
    for (i=0; i<igraph_vector_size(v); i++) {
        fprintf(f, " %f",  VECTOR(*v)[i]);
    }
    fprintf(f, "\n");
}

void print_vector(igraph_vector_t *v, std::string exp){
    
    fprintf(stdout, "%s\n",exp.c_str());
    call_print_vector(v, stdout);
}

void call_print_real_t_array(igraph_real_t * v, long int length,FILE *f){
    long int i;
    
    for(i=0;i<length;i++){
         fprintf(f, " %f",  *(v+i));
    }
    fprintf(f,"\n");
    
}


void print_real_t_array(igraph_real_t * v, long int length,std::string exp){
    fprintf(stdout,"%s\n",exp.c_str());
    call_print_real_t_array(v, length, stdout);
}


void cal_print_matrix(const igraph_matrix_t *m) {
    long int nrow=igraph_matrix_nrow(m);
    long int ncol=igraph_matrix_ncol(m);
    long int i, j;
    igraph_real_t val;
    
    for (i=0; i<nrow; i++) {
        printf("%li:", i);
        for (j=0; j<ncol; j++) {
            val = MATRIX(*m, i, j);
            if (igraph_is_inf(val)) {
                if (val < 0) {
                    printf("-inf");
                } else {
                    printf(" inf");
                }
            } else {
                printf(" %3.3f", val);
            }
        }
        printf("\n");
    }
}


void print_matrix(const igraph_matrix_t *m,std::string exp){
    fprintf(stdout, "%s\n",exp.c_str());
    cal_print_matrix(m);
}


void print_sp_matrix(const igraph_spmatrix_t *m, std::string exp){
     fprintf(stdout, "%s\n",exp.c_str());
     igraph_spmatrix_fprint(m, stdout);
}


void write_vector_file(igraph_vector_t *res, string filePath){
    ofstream myfile(filePath);
    long int nv=igraph_vector_size(res);
    if(!myfile.is_open()){cout<<"can not open file"<<filePath<<endl;}
    for(long int i=0;i<nv;i++){
        myfile<<VECTOR(*res)[i]<<endl;
    }
    myfile.close();
    cout<<"write"<<nv<<"line"<<endl;
}

void init_vector_file(igraph_vector_t *res, string filePath){
    
    ifstream myfile(filePath);
    long int size;
    string line;
    
    size=igraph_vector_size(res);
    
    if(!myfile.is_open()){
        cout<<"can not open "<<filePath<<endl;
        exit(0);
    }
    
    for(long int i=0;i<size;i++){
        getline(myfile,line);
        igraph_vector_set(res,i, stod(line));
    }
}

// just for debugging purpose
void vector_statstic(igraph_vector_t *input){
    
    double max,min,mean,sum;
    
    long int size=igraph_vector_size(input);
    max=igraph_vector_max(input);
    min=igraph_vector_min(input);
    sum=igraph_vector_sum(input);
    mean=(double)sum/size;
    
    cout<<"basic stastic for this vector"<<endl;
    cout<<"max"<<"\t"<<"min"<<"\t"<<"mean"<<endl;
    
    cout<<max<<"\t"<<min<<"\t"<<mean<<endl;
    
    cout<<"min get by:"<<igraph_vector_which_min(input);
    
    cout<<"suggested interval: "<<(max-min)/200<<endl;
    
}

