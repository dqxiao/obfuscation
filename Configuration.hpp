//
//  Configuration.hpp
//  testGraphCplus
//
//  Created by dongqingxiao on 4/11/16.
//  Copyright © 2016 dongqingxiao. All rights reserved.
//

#ifndef Configuration_hpp
#define Configuration_hpp

#include <stdio.h>



extern int k;
extern double c;
extern int attempt;
extern double epsilon;
extern double noise;
extern char delimiter;
extern int sampleNum; // the number of sample
extern bool uncertain;


enum Option {randPert=1,greedPert=2};
enum WorkDataset {dblp=1, hepth=2};


enum FeatureCombineOption {mutiply=1,fplus=2};
extern double fplusParmater;

extern Option option;
extern FeatureCombineOption foption;
extern WorkDataset w_dataset;

#endif /* Configuration_hpp */
