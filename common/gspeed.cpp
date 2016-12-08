#include "GWASSpeedup.h"
#include <time.h>

/* This is a simple wrapper or test-harness around the GWASSpeedup 
 contestant code.  It reads input data (phenotypes, genotypes, covariates) 
 from STDIN and passes it into the contestant code public method
 computeAssociations().  The output goes to STDOUT, and consists of:

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

int main()
{
  clock_t begin, end;
  double time_spent;
  int N, M, P, C;
  double x;
  vector<string> genotypes, phenotypes;
  vector<double> covariates;

  //phenotype matrix
  cin >> N >> P;
  phenotypes.resize(N);
  for (int i=0;i<N;++i){
    cin >> phenotypes[i];
  }
  //genotype matrix
  cin >> N >> M;
  genotypes.resize(N);
  for (int i=0;i<N;++i){
    cin >> genotypes[i];
  }
  //covariate matrix
  cin >> N >> C;
  covariates.clear();
  for (int i=0;i<N;++i){
    for (int j=0;j<C;++j){
      cin >> x;
      covariates.push_back(x);
    }
  }
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

