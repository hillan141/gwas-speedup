#include <unistd.h>
#include <time.h>

#include "GWASSpeedup.h"
#include "read_bed.h"

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
  while ((c = getopt (argc, argv, "b:")) != -1)
         switch (c)
           {
	   case 'b':
	     filestub = string(optarg);
	     break;
	   default:
	     abort();
	   }

  readPlinkData(filestub, genotypes, phenotypes, covariates);
  M = genotypes.size();
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
  printf("%0.3f\n", time_spent);
  fflush(stdout);
  return 0;
}

