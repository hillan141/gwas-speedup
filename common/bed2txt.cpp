#include <unistd.h>
#include <time.h>
#include "read_bed.h"
#include <iostream>
using namespace std;

/* Convert PLINK .bed/.bam/format to the simple text format
 needed by GWASSpeedup contest code
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
    cerr << "you must specify an input .bed file using the '-b' option" << endl;
    exit(-1);
  }
  readPlinkData(filestub, genotypes, phenotypes, covariates);
  dumpStdout(genotypes, phenotypes, covariates);
  return 0;
}

