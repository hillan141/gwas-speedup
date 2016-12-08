#ifndef __INDIV_H__
#define __INDIV_H__

#include <string>
#include <vector>

#include "gvar.h"
#include "Family.h"
#include "WMLocus.h"

typedef vector<double> vector_t;

class Individual;
class WMLocus;

using namespace std;

class Individual { 
 public:

  string fid;
  string iid;

  Individual() { 
    fid=iid=pat=mat=""; 
    ip=im=-1;
    sex=false; phenotype=-9; 
    sexcode="";
    aff=false;
    covar=-9;
    bcovar=false;
    clist.resize(0);
    clistMissing.resize(0);
    plist.resize(0);
    plistMissing.resize(0);
    missing=false;
    missing2=false;
    flag=true;
    one.resize(0);
    two.resize(0);
    sol=0;
    founder=true;
    pp=pm=NULL;
    family=NULL;
    kids.resize(0);
    pperson=this;
    T=W=B=0;
    gvar.resize(0);
  }


  // Parental codes
  string pat;      
  string mat;

  // Pointers to parents
  Individual * pp;
  Individual * pm;

  // Parent slot number
  int ip;
  int im;

  // Relatedness functions
  int countMeioses(Individual*);

  // Permuted self
  Individual * pperson;

  // Children (pointers, slot numbers)
  vector<Individual*> kids;
  vector<int> ikids;

  bool sex;
  string sexcode;
  double phenotype;
  bool aff;
  double covar;
  bool bcovar;

  vector_t clist; // multiple covariates
  vector<bool> clistMissing;
  
  vector_t plist; // multiple phenotypes
  vector<bool> plistMissing;
  
  bool missing;
  bool missing2;
  bool flag; 

  int sol;
  bool founder;
  Family * family;

  // SNP data
  vector<bool> one; // Person-major mode genotypes
  vector<bool> two;

  vector<bool>::iterator i1;
  vector<bool>::iterator i2;

  // Generic variant data
  vector<GVariant*> gvar;

  

  // Weighted, multi-allelic single marker
  
  WMLocus wmlocus;

  // For QFAM, within and total scores (temporary variables)
  double T;
  double B;
  double W;
};

#endif

