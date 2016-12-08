#include <iostream>
#include <time.h>
#include "Individual.h"
#include "Family.h"
#include "CSNP.h"
#include "PlinkIO.h"
#include "CArgs.h"
#include "Z.h"
#include "options.h"
#include "helper.h"
#include "read_bed.h"

using namespace std;

template <typename T>
void delete_pointed_to(T* const ptr)
{
  // cerr << "delete_ptr_to" << endl;
  delete ptr;
}

void dumpStdout( vector<string> &genotypes, vector<string> &phenotypes, vector<double> &covariates) {
  int N, P, M, C;
  N = genotypes.size();
  P = 1;
  M = genotypes[0].size();
  C = (int)covariates.size() / N;

  cout << N << endl;
  cout << P << endl;
  for (int i=0;i<N;i++) {
    cout << phenotypes[i] << endl;
  }
  
  cout << N << endl;
  cout << M << endl;
  for (int i=0;i<N;i++) {
    cout <<  genotypes[i] << endl;
  }

  cout<< N << endl;
  cout << C << endl;
  int k=0;
  for (int i=0;i<N;i++) {
    for (int j=0;j<C;j++) {
      cout << covariates[k++] << endl;
    }
  }
} // dumpStdout


void readPlinkData( string fileStub, vector<string> &genotypes, vector<string> &phenotypes, vector<double> &covariates) {
  // read Plink files (.fam, .bim, .bed, .cov) and load variables genotypes, phenotypes, covariates
  
  clock_t begin, end;
  double time_spent;

  vector<CSNP*> SNP;
  vector<Individual*> sample; 
  vector<Locus*> locus; 
  PlinkIO *pio = new PlinkIO();
  // cout << "Starting test_read_bin.cpp" << endl;
  // cout << "par::verbose : " << par::verbose << endl;

  string myBedFile = fileStub + ".bed"; 
  string myFamFile = fileStub + ".fam";
  string myMapFile = fileStub + ".bim";
  string myCovFile = fileStub + ".cov";

  pio->checkFileExists(myFamFile);
  pio->checkFileExists(myBedFile);
  pio->checkFileExists(myMapFile);
  pio->checkFileExists(myCovFile);
  
  par::silent = 1;
  par::famfile = myFamFile;
  par::bitfilename_map = myMapFile;
  par::bitfilename = myBedFile;
  par::covar_filename = myCovFile;
  par::clist_filename = myCovFile;
  
  begin = clock();
  pio->readBinData( myFamFile, myMapFile, myBedFile);
  end = clock();
  time_spent = (double)(end-begin)/CLOCKS_PER_SEC;
  cerr << "Read Plink fam/bim/bed files: " << time_spent << " sec." << endl;

  begin = clock();
  pio->readCovListFile();
  end = clock();
  time_spent = (double)(end-begin)/CLOCKS_PER_SEC;
  cerr << "Read covariate list file: " << time_spent << " sec." << endl;
  //dimensions of the problem  
  int N, P, M, C;
  sample = pio->sample;
  N = sample.size();
  locus = pio->locus;
  M = locus.size();
  P = 1;
  SNP = pio->SNP;
  C = sample[0]->clist.size();
  
  int numvalid = 0;
  vector<int> ivalid;
  ivalid.resize(0);
  for (int i=0;i<N;i++) {
    // cout << sample[i]->iid << " " << sample[i]->missing << endl;
    if (!sample[i]->missing) {
      ivalid.push_back(i);
    }
  }
  cerr << "total non-missing individuals = " << ivalid.size() << endl;

  // load the output variables
  
  //phenotypes
  begin = clock();
  phenotypes.resize(ivalid.size());
  for (int i=0; i<ivalid.size(); i++) {
    double pheno = sample[ivalid[i]]->phenotype;
    // contestant code wants 0/1/9 coding, not 1/2/-9 PLINK-style coding
    if (pheno==-9) pheno = 9;
    if (pheno==2 || pheno==1) pheno -=1;
    std::ostringstream strs;
    strs << pheno;
    phenotypes[i] = strs.str();
  }
  end = clock();
  time_spent = (double)(end-begin)/CLOCKS_PER_SEC;
  cerr << "Load phenotypes: " << time_spent << " sec." << endl;

  // marker dosages
  begin = clock();
  genotypes.resize(ivalid.size());
  vector<CSNP*>::iterator iSNP;
  int temp;
  // char buf[2];
  // traverse dosage in order (donor, marker)
  //outer loop: donor
  for (int i=0; i<ivalid.size(); i++) {
    iSNP = SNP.begin();
    genotypes[i] = "";
    // genotypes[i].reserve(SNP.size());
    std::ostringstream strs;
    while (iSNP != SNP.end()) {
      // cerr << sample[ivalid[i]]->iid << " " << (*iSNP)->one[ivalid[i]] << ":" <<  (*iSNP)->two[ivalid[i]] << endl;
      temp = (2-((*iSNP)->one[ivalid[i]]+(*iSNP)->two[ivalid[i]]));
      // sprintf(buf, "%d", temp);
      // genotypes[i].append(buf, 1);
      strs << temp;
      iSNP++;
    } // end loci for donor i
    genotypes[i] = strs.str();
  } // end donor i
  end = clock();
  time_spent = (double)(end-begin)/CLOCKS_PER_SEC;
  cerr << "Load genotypes: " << time_spent << " sec." << endl;

  // covariates
  begin = clock();
  covariates.clear();
  covariates.resize(0);
  for (int i=0; i<ivalid.size(); i++) {
    for (int j=0; j<C; j++) {
      covariates.push_back( sample[ivalid[i]]->clist[j]);
    }
  }
  end = clock();
  time_spent = (double)(end-begin)/CLOCKS_PER_SEC;
  cerr << "Load covariates: " << time_spent << " sec." << endl;

  // clear the pio pointed-to-objects - will not be freed on deletion of pio

  // free memory
  begin = clock();
  for_each(pio->locus.begin(), pio->locus.end(), delete_pointed_to<Locus>);
  for_each(pio->SNP.begin(), pio->SNP.end(), delete_pointed_to<CSNP>);
  for_each(pio->sample.begin(), pio->sample.end(), delete_pointed_to<Individual>);

  delete pio; // large IO object no longer needed
  end = clock();
  time_spent = (double)(end-begin)/CLOCKS_PER_SEC;
  cerr << "Free memory: " << time_spent << " sec." << endl;

} // readPlinkData
