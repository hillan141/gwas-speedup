#include <unistd.h>
#include <time.h>

#include "read_bed.h"

int main(int argc, char **argv)
{
  clock_t begin, end;
  double time_spent;
  int N, M, P, C;
  double x;
  vector<string> genotypes, phenotypes, locus_name;
  vector<double> covariates;

  int c;
  string filestub;
  int headflag = 0;
  while ((c = getopt (argc, argv, "hb:")) != -1)
         switch (c)
           {
	   case 'b':
	     filestub = string(optarg);
	     break;
	   case 'h':
	     headflag = 1;
	   default:
	     abort();
	   }

  readPlinkData(filestub, genotypes, phenotypes, covariates, locus_name, 0, 1);
  dumpStdout(genotypes, phenotypes, covariates);
  return 0;
}

