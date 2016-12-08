#include "GWASSpeedup.h"
#include <time.h>

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

