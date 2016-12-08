#ifndef __FAMILY_H__
#define __FAMILY_H__

#include "Individual.h"

class Family;
class Individual;

using namespace std;

class Family
{
 public:

  Family() 
    { 
      include = false;
      parents = false;
      discordant_parents = false;
      singleton = false;
      sibship = false;
      TDT = false;
      pat = mat = NULL;
      kid.clear();
    }
  
  void copy(const Family & rhs)
    {
      include = rhs.include;
      pat = rhs.pat;
      mat = rhs.mat;
      kid.clear();
      for (unsigned int c=0; c<rhs.kid.size(); c++)
	kid.push_back(rhs.kid[c]);
    }

  Family & operator = (const Family & rhs)
  {
    copy(rhs);
    return *this;
  }
  
  bool include;
  bool parents;
  bool sibship;
  bool discordant_parents;
  bool singleton;
  bool TDT;

  Individual * pat;
  Individual * mat;
  vector<Individual *> kid;

  // Between-family genotypic score
  double B;
};

#endif

