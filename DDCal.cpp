//
//  DDCal.cpp
//  testGraphCplus
//
//  Created by dongqingxiao on 3/21/16.
//  Copyright © 2016 dongqingxiao. All rights reserved.
//

#include "DDCal.hpp"
#include "Configuration.hpp"
#include <math.h>
#include <boost/math/distributions/normal.hpp>
#include <algorithm>
#include "Help.hpp"
//#include <limits.h>
using boost::math::normal;

/*
 * O(n^2) where n is the size of input
 */
void exact_cal(const igraph_vector_t * input, igraph_vector_t *res,long int degree,long int maxDegree){
    
    igraph_vector_t temp;
    igraph_vector_init(&temp,degree+1);
    
    igraph_vector_set(&temp, 0, 1);
    
    for(long int i=0;i<degree;i++){
        igraph_real_t p=VECTOR(*input)[i];        
        for(long j=i+1;j>=1;j--){
            igraph_vector_set(&temp,j,VECTOR(temp)[j]*(1-p)+VECTOR(temp)[j-1]*p);
        }
        
        igraph_vector_set(&temp,0,VECTOR(temp)[0]*(1-p));
        
    }
    
    // copy to real result 
    
    for(long int i=0;i<std::min(degree,maxDegree)+1;i++){
        igraph_vector_set(res,i, VECTOR(temp)[i]);
    }

    igraph_vector_destroy(&temp);
    // done 
    
}
/*
 * O(n) where n is the size of input
 */
void normal_cal(const igraph_vector_t * input, igraph_vector_t *res, long int degree,long int maxDegree){
    
    igraph_real_t mu=igraph_vector_sum(input);
    igraph_real_t sigma=0;
    igraph_real_t need;
    for(long int i=0;i<degree;i++){
        igraph_real_t p=VECTOR(*input)[i];
        if(p>1 || p<0){
            cout<<"p"<<p<<endl;
            throw std::exception();
        }
        
        sigma+=p*(1-p);
    }
    
    
    sigma=sqrt(sigma);
    
    if(isnan(sigma) || sigma<=0){
        //vector_statstic(input);
        
        throw std::exception();
    }
    if(isnan(mu) || mu<=0){
         throw std::exception();
    }
    
    normal ns(mu,sigma);
    double lval,hval;
    double realVal;
    lval=boost::math::cdf(ns,-0.5);
    need=std::min(degree, maxDegree);
    for(long int i=0;i<need+1;i++){
        hval=boost::math::cdf(ns,i+0.5);
//
        realVal=hval-lval;
        if(realVal<100*std::numeric_limits<double>::epsilon()){
            realVal=boost::math::pdf(ns,i-1)+boost::math::pdf(ns,i);
            realVal/=2;
        }
        igraph_vector_set(res, i, realVal);
        lval=hval;
    }
    
    
}

void DDCal(const igraph_vector_t * input,  igraph_vector_t * res, igraph_real_t degree, igraph_real_t maxDegree){
    
    //
    if(degree<20){
        exact_cal(input,res,degree,maxDegree);
    }else{
        normal_cal(input, res,degree,maxDegree);
    }

}


igraph_real_t cal_entropy(igraph_real_t p){
    
    if(isnan(p)){
        
    }
    if(p<=0){
        if(p!=0){
        }
        return 0;
    }
    return -1*p*log2(p);
}


double cal_distance(long int i, long int j){
    double t=std::abs(i-j);
    return t;
    
}


double cal_distance_matrix(igraph_matrix_t * one, igraph_matrix_t * second, igraph_real_t nv){
    
    double dis_sum=0.0;
    
    
    for(long int i=0;i<nv;i++){
        for(long int j=i+1;j<nv;i++){
            double oVal=MATRIX(*one, i,j);
            double sVal=MATRIX(*second, i,j);
            dis_sum+=std::abs(oVal-sVal);
        }
    }
    
    return dis_sum;
}


double cal_distance_vector(igraph_vector_t * one, igraph_vector_t * second){
    
    double dis_sum=0.0;
    long int oSize=igraph_vector_size(one);
    long int sSize=igraph_vector_size(second);
    if(oSize!=sSize){
        throw std::exception();
    }
    
    
    for(long int i=0;i<oSize;i++){
        dis_sum+=std::abs((double)VECTOR(*one)[i]-(double)VECTOR(*second)[i]);
    }
    
    return dis_sum;
}

double cal_mean_error_vector(igraph_vector_t * ref, igraph_vector_t * test){
    
    double dis_sum=0;
    
    long int refLen=igraph_vector_size(ref);
    long int testLen=igraph_vector_size(test);
    
    if(refLen!=testLen){
        throw std::exception();
    }
    
    for(long int i=0;i<refLen;i++){
        double mean_error=std::abs((double)VECTOR(*ref)[i]-(double)VECTOR(*test)[i]);
        dis_sum+=mean_error;
    }
    
    dis_sum/=refLen-1;
    dis_sum/=refLen;
    
    return dis_sum;
}


double cal_relative_error_vector(igraph_vector_t * ref, igraph_vector_t * test){
//    double dis_sum=cal_distance_vector(one, second);
//    
//    double val=igraph_vector_sum(one);
//    
//    
//    return dis_sum/val;
    
    double relError=0.0;
    
    long int refLen=igraph_vector_size(ref);
    long int testLen=igraph_vector_size(test);
    
    if(refLen!=testLen){
        throw std::exception();
    }
    
    for(long int i=0;i<refLen;i++){
        double mean_error=std::abs((double)VECTOR(*ref)[i]-(double)VECTOR(*test)[i]);
        
        double maxVal=std::max((double)VECTOR(*ref)[i], (double)VECTOR(*test)[i]);
        if(maxVal!=0){
            double rel_error=mean_error/maxVal;
            relError+=rel_error;
        }
        
        
        
    }
    
    return relError/refLen;
    
}



long int permuateCal(long int number){
    return number*(number-1);
}


void is_sparese(igraph_vector_t * input, long int nv){
    
    long int sum=0;
    for(long int i=0;i<nv;i++){
        sum+=VECTOR(*input)[i];
    }
    
    //
    printf("average:%f \n", (double)sum/(nv*(nv-1)) );
    
    
}


void reDistribute(igraph_vector_t * uv){
    
    double valSum=igraph_vector_sum(uv);
    igraph_vector_scale(uv, 1.0/valSum);
    
}


void reverseVector(igraph_vector_t * ruv){
    
    double minVal=igraph_vector_min(ruv);
    double maxVal=igraph_vector_max(ruv);

    double diff=maxVal-minVal;
    int  size=(int)igraph_vector_size(ruv);
    
    double sigma=standard_vaiance_vector(ruv);
    
    normal ns(0,0.4);
    
    for(int i=0;i<size;i++){

        
        double val=VECTOR(*ruv)[i];
        
        val=(val-minVal)/diff; // into 0,1
        //val=boost::math::pdf(ns,val);
        
        igraph_vector_set(ruv,i,val);
    }
    
}




double featureCombine(double uniq, double robust){
    //option 1: mutilple
    if(foption==mutiply){
        return uniq*robust;
    }
    if(foption==fplus){
        return uniq+fplusParmater*robust;
    }
    return uniq;
}


double standard_vaiance_vector(igraph_vector_t * ruv){
    // E(x^2)-E(x)^2
    
    double exSum=0;
    double ex2Sum=0;
    double sigma=0.1;
    int size=(int)igraph_vector_size(ruv);
    
    for(long int i=0;i<size;i++){
        double val=VECTOR(*ruv)[i];
       
        exSum+=val;
        ex2Sum+=(val*val);
    }
    
    
    double varX=(ex2Sum/size)-(exSum/size)*(exSum/size);
    
    sigma=sqrt(varX);
    if(sigma==0){
        return 0.1;
    }
    
    return sigma;
}

double relative_error_metric(igraph_vector_t * one, igraph_vector_t * second, int i){
    
    double result=0;
    
    result=std::abs(VECTOR(*one)[i]-VECTOR(*second)[i]);
    
    result/=std::abs(VECTOR(*one)[i]);
    
    return result;
}
