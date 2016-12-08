#ifndef __LOGISTIC_H__
#define __LOGISTIC_H__

class LogisticModel {
      public:

          size_t nind;       // number of individuals: nind = N - (# missing values in marker col)
          size_t np;         // number of regressors:  np = C+2

          matrix_t Xt; // size(Xt) = np * nind

          bool all_valid; // flag to record errors
          bool S0_made;
          
          matrix_t S0;

          // logistic model: Y ~ logistic(X * coef)
          vector_t Y; // size(Y) = nind, whether a individual has the phenotype

          // - all-1s vector (const term)
          // - marker vector
          // - C covariate vectors
          vector_t coef;  // size(coef) = np

          // variables used internally by fitLM()
          matrix_t S; // size(S) = np*np, covariance between regressors
          vector_t p;     // probabilities: logistic(X * coef)
          vector_t V;     // diagonal p(1-p)

          vector_t t1, t2, s1, ncoef; // temporary vectors
          matrix_t ST; // temporary matrix
          
      public:
          //=================from model.cpp============================
        
          /**
           * Check for multicollinearity (near-dependence of regressors, i.e., 
      columns of X)
           * Set all_valid = false (and subsequently Z=0) in the rare case of 
      multicollinearity
      #if 0
          bool checkVIF()
          {
              TIMER(checkVIF);
              // Calculate correlation matrix for X
              // Skip intercept (all-1s vector stored in first col of X)

              int p = X.size();
              if (p<2) return false;

              int q = X[0].size() - 1;     
              if ( q < 2 ) return true;

              vector_t m(q);
              matrix_t c, ct;
              sizeMatrix(c,q,q);

              for (int i=0; i<p; i++)
                  for (int j=0; j<q; j++)
                      m[j] += X[i][j+1];

              for (int j=0; j<q; j++)
                  m[j] /= (double)p;
            
              for (int i=0; i<p; i++)
                  for (int j1=0; j1<q; j1++)
                      for (int j2=j1; j2<q; j2++)
                          c[j1][j2] += ( X[i][j1+1] - m[j1] ) * ( X[i][j2+1] - 
      m[j2] );
            
              for (int j1=0; j1<q; j1++)
                  for (int j2=j1; j2<q; j2++)
                      c[j1][j2] /= (double)(p-1);
            
              for (int j1=0; j1<q; j1++)
                  for (int j2=j1+1; j2<q; j2++) {
                      c[j1][j2] /= sqrt( c[j1][j1] * c[j2][j2] );
                      c[j2][j1] = c[j1][j2];
          
                      if ( c[j2][j1] > 0.999 ) {
                          return false;
                      }
                  }        

              // Any item with zero variance?
              for (int j=0; j<q; j++) {
                  if ( c[j][j] == 0 || ! realnum( c[j][j] ) )
                      return false;
                  c[j][j] = 1;
              }
              
              // Get inverse
              if ( !invertMatr(c, ct) ) {
                  return false;
              }

              // Calculate VIFs
              //double maxVIF = 0;
              for (int j=0;j<q;j++) {
                  // r^2 = 1 - 1/x where x is diagonal element of inverted
                  // correlation matrix
                  // As VIF = 1 / ( 1 - r^2 ) , implies VIF = x

                  if ( c[j][j] > 50 ) //par::vif_threshold = 50 in options.cpp, 
      default
                      return false;
              }
              return true;
          }
      #endif
	  */

      #define not_set ~0u
          size_t cur_gen, cur_phe;
          typedef vector<size_t> XGen;
          XGen xgen;
          vector_t XtXKM;
          struct CacheXGen {
              matrix_t Xt;
          };
          struct CachePhe {
              bool all_valid;
              vector_t coef;
              matrix_t S;
      #if processing_order == by_estimate_start1 || processing_order == by_estimate_start2
      # define SAVE_Y
              vector_t Y;
      #endif
      #if processing_order == by_estimate_start1 || processing_order == by_estimate_start2
      # define SAVE_P
              vector_t p;
      #endif
          };
          struct CacheGen {
              bool all_valid;
              bool S0_made;
              XGen xgen;
              vector_t XtXKM;
              matrix_t S0;
              vector<CachePhe> cache_Phe;
          };
          map<XGen, CacheXGen> cache_XGen;
          vector<CacheGen> cache_Gen;

          vector< vector<char> > phenotypesT, genotypesT;
          matrix_t ind2cov;

	  // Constructor
	  LogisticModel(const vector<string>& phenotypes, const vector<string>& genotypes, const vector<double>& covariates);

          CacheGen* get_cache_gen(size_t i, bool alloc)
          {
              if ( cache_Gen.empty() ) {
                  if ( !alloc ) return 0;
                  cache_Gen.resize(genotypesT.size());
                  forIter ( i, cache_Gen ) {
                      i->all_valid = 1;
                      i->S0_made = 0;
                  }
              }
              return &cache_Gen[i];
          }

          CachePhe* get_cache_phe(size_t i, size_t j, bool alloc)
          {
              CacheGen* c = get_cache_gen(i, alloc);
              if ( !c ) return 0;
              if ( c->cache_Phe.empty() ) {
                  if ( !alloc ) return 0;
                  c->cache_Phe.resize(phenotypesT.size());
                  forIter ( i, c->cache_Phe ) {
                      i->all_valid = true;
                  }
              }
              return &c->cache_Phe[j];
          }

	  void unselect_phenotype(bool save = 1);
	  void select_phenotype(size_t phe);
	  void unselect_genotype(bool save = 1);
	  void select_genotype(size_t gen);
	  void update_S0_25();
	  void update_S_V();
          // fit the logistic model
          void fitLM_genotype();
         
      #define DELTA_EPS 1e-3
      #define START_STEPS 2
      #define DELTAS_EPS 1e-10
          //#define NCOEF_EPS 1e-8
	  double getResult(size_t k);

          double fitLM_phenotype_estimate(size_t k);
          // fit the logistic model
	  // this is the default "full fit"
          double fitLM_phenotype();
        
          // fit the logistic model
          double fitLM_phenotype_start(size_t k);

          // fit the logistic model
	  double fitLM_phenotype_step(size_t k, size_t it);
          // fit the logistic model
	  double fitLM_phenotype_continue(size_t it0);
          // return variance of each estimated coefficient
          vector_t getVar();
          // return variance of each estimated coefficient
	  double getVar(int testParameter);
          // return association statistic (chi-square from Wald test) computed by logistic model
          double getAssociationStat(int testParameter);
      };

#endif // __LOGISTIC_H__
