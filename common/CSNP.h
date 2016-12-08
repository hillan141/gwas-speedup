#ifndef __CSNP_H__
#define __CSNP_H__

#include <vector>

using namespace std;

// Main genotype storage, ordered by SNP
class CSNP
{
 public:

  vector<bool> one; // SNP-major mode genotypes
  vector<bool> two; 
  
};

#endif

