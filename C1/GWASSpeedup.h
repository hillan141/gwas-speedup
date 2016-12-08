#ifndef GWAS_H
#define GWAS_H

#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <iostream>
#include <algorithm>
#include <cassert>
#include <string>
#include <cmath>
#include <vector>
#include <sstream>
using namespace std;

class GWASSpeedup
{
public:
  template <class T>
      void buildMatrix(vector<string> &src,vector< vector<T> > &tar);	
    vector< vector<bool> > ind2phe;
    vector< vector<double> > ind2marker;
    vector< vector<double> > ind2cov;
    vector<double> computeAssociations(const vector<string> &phenotype, const vector<string> &genotype, const vector<double> &covariate);
};

#endif // GWAS_H
