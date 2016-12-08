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
  bool is_transposed_genotype = true;

  if (!is_transposed_genotype) {
    N = genotypes.size();
    M = genotypes[0].size();
  } else {
    N = genotypes[0].size();
    M = genotypes.size();
  }
  P = 1;
  C = (int)covariates.size() / N;

  cout << N << endl;
  cout << P << endl;
  for (int i=0;i<N;i++) {
    cout << phenotypes[i] << endl;
  }
  

  cout << N << endl;
  cout << M << endl;

  
  if (!is_transposed_genotype) { // original non-transposed genotypes
    for (int i=0;i<N;i++) {
      cout <<  genotypes[i] << endl;
    }
  } else {
    for (int i=0;i<M;i++) {
      cout <<  genotypes[i] << endl;
    }
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


void readPlinkData( string fileStub, vector<string> &genotypes, vector<string> &phenotypes, vector<double> &covariates, vector<string> &locus_name, int firstMarker, int lastMarker) {
  // read Plink files (.fam, .bim, .bed, .cov) and load variables genotypes, phenotypes, covariates
  clock_t begin, end;
  double time_spent;
  PlinkIO *pio = new PlinkIO();
  // cout << "Starting test_read_bin.cpp" << endl;
  // cout << "par::verbose : " << par::verbose << endl;

  string myBedFile = fileStub + ".bed"; // "../../../UC_N6000_M7000_P1_C5/plink.bed"; // test-data/CG10kNHWhg19clean_maf1pct.bed";
  string myFamFile = fileStub + ".fam"; // "../../../UC_N6000_M7000_P1_C5/plink.fam"; // test-data/CG10kNHWhg19clean_maf1pct.fam";
  string myMapFile = fileStub + ".bim"; // "../../../UC_N6000_M7000_P1_C5/plink.bim"; // test-data/CG10kNHWhg19clean_maf1pct.bim";
  string myCovFile = fileStub + ".cov"; // "../../../UC_N6000_M7000_P1_C5/UC_N6000_M7000_P1_C5.cov";
  
  par::silent = 1;
  par::famfile = myFamFile;
  par::bitfilename_map = myMapFile;
  par::bitfilename = myBedFile;
  par::covar_filename = myCovFile;
  par::clist_filename = myCovFile;
  
  begin = clock();
  //~ pio->readBinData( myFamFile, myMapFile, myBedFile);
  int locus_size = pio->readBinDataInterval(myFamFile, myMapFile, myBedFile, genotypes, locus_name, firstMarker, lastMarker);
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
  vector<Individual*> &sample = pio->sample;
  N = sample.size();
  vector<Locus*> &locus = pio->locus;
  M = locus.size();
  M = locus_size;
  P = 1;
  vector<CSNP*> &SNP = pio->SNP;
  C = sample[0]->clist.size();

  /* dump the locus_name's
  for (int i=0; i<locus_name.size(); i++) {
    cerr << "i=" << i << "name = " << locus_name[i] << endl; 
  }
  */

  int numvalid = 0;
  vector<int> ivalid;
  ivalid.resize(0);
  for (int i=0;i<N;i++) {
    if (!sample[i]->missing) {
      ivalid.push_back(i);
    }
  }
  cerr << "total non-missing individuals = " << ivalid.size() << endl;

  // load the output variables
  
  //phenotypes
  begin = clock();
  phenotypes.resize(ivalid.size());
  for (int i=0, vsize=ivalid.size(); i<vsize; i++) {
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
  //~ genotypes.clear();
  //~ genotypes.resize(SNP.size());
  //~ for (int i=0, vsize=SNP.size(); i < vsize; ++ i) {
    //~ SNP[i]->getString(genotypes[i]);
  //~ }
  
  end = clock();
  time_spent = (double)(end-begin)/CLOCKS_PER_SEC;
  cerr << "Load genotypes: " << time_spent << " sec." << endl;

  // covariates
  begin = clock();
  {
	  int vsize = ivalid.size();
	  covariates.clear();
	  covariates.resize(vsize*C);
	  int k=0;
	  for (int i=0; i<vsize; i++) {
		for (int j=0; j<C; j++) {
		  covariates[k++] = sample[ivalid[i]]->clist[j];
		}
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
