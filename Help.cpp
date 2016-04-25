//
//  Help.cpp
//  testGraphCplus
//
//  Created by dongqingxiao on 3/19/16.
//  Copyright © 2016 dongqingxiao. All rights reserved.
//

#include "Help.hpp"
#include <random>
#include <fstream>
#include <set>

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


void print_map(map<int,vector<int> > m){
    for(auto mitem:m){
        cout<<mitem.first<<":";
        for(auto p:mitem.second){
            cout<<p<<",";
        }
        cout<<endl;
    }
    //
}

void print_boost_matrix(boost::numeric::ublas::matrix<int> &m){
    int nrow=m.size1();
    int ncol=m.size2();
    
    for(int i=0;i<nrow;i++){
        for(int j=0;j<ncol;j++){
            if(j==0){
                cout<<m(i,j);
            }else{
                cout<<"\t"<<m(i,j);
            }
        }
        cout<<endl;
    }
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

// pass by reference 
void write_boostMatrix_file(boost::numeric::ublas::matrix<int>&m, string filePath){
    ofstream myfile(filePath);
   // long int nv=igraph_vector_size(res);
    if(!myfile.is_open()){cout<<"can not open file"<<filePath<<endl;}
    
    int nrow=(int)m.size1();
    int ncol=(int)m.size2();
    
    for(int i=0;i<nrow;i++){
        for(int j=0;j<ncol;j++){
            if(j==0){
                myfile<<m(i,j);
            }else{
                myfile<<'\t'<<m(i,j);
            }
            
        }
        myfile<<endl;
    }
    myfile.close();
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

// pass by reference
void init_boostMatrix_file(boost::numeric::ublas::matrix<int> &m, string filepath){
    ifstream myfile(filepath);
    string line;
    string item;
    
    
    if(!myfile.is_open()){
        cout<<"can not open "<<filepath<<endl;
        exit(0);
    }
    
    // reverse the matrix
    int nrow=0;
    int ncol=0;
    while(getline(myfile,line)){
        ncol=0;
        stringstream ss(line);
        while(getline(ss,item,'\t')){
            m(nrow,ncol)=stoi(item);
            ncol+=1;
        }
        nrow+=1;
    }
    // matrix: each row is one ...
    
    
    myfile.close();
   
    cout<<"finish load file"<<filepath<<endl;
    
    // done 
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
    
    cout<<"min get by:"<<igraph_vector_which_min(input)<<endl;
    cout<<"max get by:"<<igraph_vector_which_max(input)<<endl;
    
    cout<<"suggested interval: "<<(max-min)/200<<endl;
    
}


void reorginsedByrow(boost::numeric::ublas::matrix<int> &m, map<int, vector<int> > & rep, map<int, vector<int> > & repPos){
    // input: matrix
    // nrow=nv
    // ncol=sampleNum
    
    int nrow=(int)m.size1();
    int ncol=(int)m.size2();
    
    
    //for each row
    for(int i=0;i<nrow;i++){
       //used for search 
        vector<SampleRep> temp;
        
        for(int j=0;j<ncol;j++){
            int val=m(i,j);
            temp.push_back(SampleRep(j,val));
        }
        
        sort(temp.begin(),temp.end());
        
        int c_rep=-1;
        int pos=0;
        for(auto sr: temp){
            if(c_rep!=sr.repNum){
                rep[i].push_back(sr.repNum);
                repPos[i].push_back(pos);
                c_rep=sr.repNum;
            }
            m(i,pos)=sr.sampleNum;
            pos++;
        }
    
        
        temp.clear();
    
    }
    
    cout<<"reorgianzed by row aggregation"<<endl;
    
}


void intersection_result(vector<int> rep_f, vector<int> rep_s,vector<int> &result){
    
    int i=0;
    int j=0;
    
    int m=(int)rep_f.size();
    int n=(int)rep_s.size();
    
    while(i<m && j<n){
        if(rep_f[i]<rep_s[j]){
            i++;
        }else if (rep_f[i]> rep_s[j]){
            j++;
        }else{
            // equal
            // push current index into result
            result.push_back(i);
            result.push_back(j);//
            i++;
            j++;
        }
        
    }
    
}


int intersection_cout(vector<int> rep_f, vector<int> rep_s, int fstart, int fend, int sstart, int ssend){
    
    int i=fstart;
    int j=fend;
    
    int m=fend;
    int n=ssend;
    
    int count=0;
    while(i<m && j<n){
        if(rep_f[i]<rep_s[j]){
            i++;
        }else if (rep_f[i]> rep_s[j]){
            j++;
        }else{
            // equal
            count+=1;
            i++;
            j++;
        }
        
    }
    
    return count;
}


int intersection_matrix_cout(boost::numeric::ublas::matrix<int> &matrix, int fstart, int fend, int sstart, int ssend, int first, int second){
    
    int i=fstart;
    int j=sstart;
    
    int m=fend;
    int n=ssend;
    
    int count=0;
    
    while(i<m && j<n){
        if(matrix(first,i)<matrix(second,j)){
            i++;
        }else if (matrix(first,i)>matrix(second,j)){
            j++;
        }else{
            // equal
            count+=1;
            i++;
            j++;
        }
        
    }
    
    return count;
}


double compareRowDirect(boost::numeric::ublas::matrix<int> &m, int first, int second){
    
    double result=0;
    
    int ncol=(int)m.size2();
    
    for(int j=0;j<ncol;j++){
        if(m(first,j)==m(second,j)){
            result+=1;
        }
        
    }
    return result;
}


int compareRow(boost::numeric::ublas::matrix<int> &m, map<int, vector<int> > & rep, map<int, vector<int> > & repPos,int first,int second){

    double result=0;
    
    int lastElement=(int)m.size2(); // sampleNum
    vector<int> rep_f=rep[first];
    vector<int> rep_s=rep[second];
    vector<int> rep_p_f=repPos[first];
    vector<int> rep_p_s=repPos[second];
    
    
    vector<int> rep_inter;
    
    
    intersection_result(rep_f,rep_s,rep_inter);
    
    if(rep_inter.size()==0){
        return 0; // no common one
    }
    
    for(int i=0;i<rep_inter.size();i+=2){
       
        int fIndex=rep_inter[i];
        int sIndex=rep_inter[i+1];
        
        int fstart=rep_p_f[fIndex];
        int fend=lastElement;
        if(fIndex!=rep_p_f.size()-1){
            fend=rep_p_f[fIndex+1];
        }
        
        int sstart=rep_p_s[sIndex];
        int send=lastElement;
        if(sIndex!=rep_p_s.size()-1){
            send=rep_p_s[sIndex+1];
        }
        
        
        int c_count=intersection_matrix_cout(m,fstart,fend,sstart,send, first,second);
        
        result+=c_count;
    }
    
    rep_f.clear();
    rep_s.clear();
    rep_p_f.clear();
    rep_p_s.clear();
    
    return result;
}


double compare_In_Matrix(boost::numeric::ublas::matrix<int> &ref, boost::numeric::ublas::matrix<int> &test){
    map<int,vector<int> > reps_mref; // repstart pos;
    map<int,vector<int> > reppos_mref;
    map<int,vector<int> > reps_mTest; // repstart pos;
    map<int,vector<int> > reppos_mTest;
    
   
    
    
    reorginsedByrow(ref, reps_mref,reppos_mref);
    reorginsedByrow(test, reps_mTest, reppos_mTest);
    
    
    
    
    
    
    
    int nv=(int)ref.size1(); // get nv
    int c_sample=(int) ref.size2(); // get the number of sample
    
    double diff=0;
    
    for(int first=0;first<nv;first++){
        // never try same pair
        for(int second=first+1;second<nv;second++){
            int refVal=compareRow(ref, reps_mref, reppos_mref, first,second);
            int testVal=compareRow(test, reps_mTest, reppos_mTest, first,second);
            int cdiff=std::abs(refVal-testVal);
            diff+=(double) cdiff/c_sample; // pair
        }
        
    }
    diff*=2;
   

    reps_mref.clear();
    reppos_mref.clear();
    reps_mTest.clear();
    reppos_mTest.clear();
    
    return diff;
}




double sampleing_compare_In_Matrix(boost::numeric::ublas::matrix<int> &ref, boost::numeric::ublas::matrix<int> &test){
    
    double diff=0.0;
    
    int nv=(int)ref.size1(); // nv
    int vsampleNum=2000;
    std::set<int> sampleVertices;
    
    if(nv==0){
        throw std::exception();
    }
//    uniform_real_distribution<int> unDist(0,nv-1);
//    

//    
//    
//    random_device rd;
//    
//    std::mt19937 gen(rd());
    
    while(sampleVertices.size()!=vsampleNum){
        int newVal=rand()%nv;
        sampleVertices.insert(newVal);
        
    }
    
    vector<int> svec(sampleVertices.begin(),sampleVertices.end());
    
    
    sort(svec.begin(),svec.end());
    int c_sample=(int)ref.size2();
    for(int i=0;i<svec.size();i++){
        for(int j=i+1;j<svec.size();j++){
            int first=svec[i];
            int second=svec[j];
            
            if(first>nv || second>nv){
                throw std::exception();
            }
            
            int refVal=compareRowDirect(ref, first, second);
            int testVal=compareRowDirect(test, first, second);
            int cdiff=std::abs(refVal-testVal);
            diff+=(double) cdiff/c_sample; // pair
        }
    }
    
    diff*=2;
    diff/=vsampleNum;
    diff/=vsampleNum-1;
    
    return diff;
}

