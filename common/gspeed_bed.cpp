#include <unistd.h>
#include <time.h>
#include <iostream>

#include "GWASSpeedup.h"
#include "read_bed.h"

using namespace std;

/* This is a simple wrapper or test-harness around the GWASSpeedup 
 contestant code.  It reads input data in .bed/.fam/.map/.cov format 
 from STDIN and passes it into the contestant code public method
 computeAssociations().  The primary output goes to STDOUT, and consists of:

 - the chi-squared stat for each marker in the input
 - the elapsed time in seconds (in the last line)

 In addition, files 'doudouille_estimate.txt' and 'doudouille_final.txt' may
 be saved to the working directory, containing the first-iteration estimate
 and final chi-squared statistics, respectively.

 Important note: the contestant code computes chi-squared statistics,
 unlike PLINK 1.07, which reports Wald statistics in the "STAT" column of
 the PLINK output file.  The chi-squared statistic is the square of the PLINK
 Wald statistic.
*/

int main(int argc, char **argv)
{
  clock_t begin, end;
  double time_spent;
  int N, M, P, C;
  double x;
  vector<string> genotypes, phenotypes;
  vector<double> covariates;

  int c;
  string filestub;
  bool gotFilestub = false;
  
  while ((c = getopt (argc, argv, "b:")) != -1)
         switch (c)
           {
	   case 'b':
	     filestub = string(optarg);
	     gotFilestub = true;
	     break;
	   default:
	     cerr << "There was a problem processing your command line arguments" << endl;
	     exit(-1);
	   }

  if (! gotFilestub) {
    cerr << "You must specify the prefix of a set of input .bed/.fam/.bim/.cov files using the '-b' option" << endl;
    cerr << "For example, if you specify '-b foobar', then the files foobar.bed, foobar.fam, foobar.bim, " << endl;
    cerr << "and foobar.cov must all exist." << endl;
    exit(-1);
  }

  readPlinkData(filestub, genotypes, phenotypes, covariates);
  M = genotypes[0].size();
  P = 1;

  GWASSpeedup* filter = new GWASSpeedup();
  begin = clock();
  vector<double> result = filter->computeAssociations(phenotypes,genotypes,covariates);
  end = clock();
  time_spent = (double)(end-begin)/CLOCKS_PER_SEC;
  assert(result.size() == M * P);
  for (int i=0;i<result.size();++i){
    printf("%.6f\n",result[i]);
  } 
  printf("computeAssociations: %0.3f sec\n", time_spent);
  fflush(stdout);
  return 0;
}

