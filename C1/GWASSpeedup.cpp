//Coder: venco 
//Submission: 8 

//#define CLOCK
//#define FTRACE

#include "GWASSpeedup.h"
#include "commons.h"
#include "logistic.h"

#include "../common/read_bed.h"
#include "../common/dcdflib.h"
#include "../common/ipmpar.h"

using namespace std;

#define sequential 0
#define by_estimate 1
#define by_start 2
#define by_estimate_start1 3
#define by_estimate_start2 4
#define processing_order by_estimate

bool realnum(double d)
{
  double zero = 0;
  if (d != d || d == 1/zero || d == -1/zero) 
    return false;
  else
    return true;
}

double chiprobP(double x, double df)
{

  if ( ! realnum(x) ) return -9;

  double p, q;
  int st = 0; // error variable
  int w = 1; // function variable
  double bnd = 1; // boundary function

  // NCP is set to 0
  cdfchi(&w,&p,&q,&x,&df,&st,&bnd);

  // Check status
  if (st != 0 ) return -9;

  // Return p-value
  return q;
}

     
//=================from logistic.cpp=========================

const int mAB = 1;
const int mBB = 2;

vector<double> GWASSpeedup::computeAssociations(const vector<string>& phenotypes,
						  const vector<string>& genotypes,
						  const vector<double>& covariates)
{
  ofstream fout;
  clock_t begin, end;
  double time_spent;
 
 size_t N = genotypes.size();
 size_t M = genotypes[0].size();
 size_t P = phenotypes[0].size();
 size_t T = M*P;
 
 vector<double> ret(T);
 
 begin = clock();
 LogisticModel model(phenotypes, genotypes, covariates);
 end = clock();
 time_spent = (double)(end-begin)/CLOCKS_PER_SEC;
 cerr << "init_data " << time_spent << endl;

 size_t np = model.np;
 double seq_time = 100; // force estimate-ordering 5./225e6*T*(N+np*np+50)*np;
 // TR(seq_time);
 
 cout << "seq_time = " << seq_time << " TIME_LIMIT = " << TIME_LIMIT  << endl;

 if ( seq_time < TIME_LIMIT || processing_order == sequential ) {
   cout << "doing sequential processing, no estimate" << endl;
   forN ( i, M ) {
     model.select_genotype(i);
     forN ( j, P ) {
       model.select_phenotype(j);
       double t = model.fitLM_phenotype();;
       ret[i*P+j] = t;
       model.unselect_phenotype(0);
     }
     model.unselect_genotype(0);
   }
 }
 else {
   const size_t TIME_BLOCK = 32;
#if processing_order == by_estimate
   cout << "doing by_estimate ordering" << endl;
   begin = clock();
   typedef pair<double, pair<size_t, size_t> > SortValue;
   vector<SortValue> sv; sv.reserve(T);
   SortValue v;
   forN ( i, M ) {
     model.select_genotype(i);
     forN ( j, P ) {
       model.select_phenotype(j);
       double t = model.fitLM_phenotype_estimate(i*P+j);
       v.first = -t;
       v.second.first = i;
       v.second.second = j;
       sv.push_back(v);
       ret[i*P+j] = t;
     }
   }
   end = clock();
   time_spent = (double)(end-begin)/CLOCKS_PER_SEC;
   cerr << "by_estimate_ordering " << time_spent << endl;   

   // first estimate for stat now complete for all M*P problems
   #ifdef ESTIMATE
   fout.open("venco_estimate.txt");
   if (!fout) {
     cout << "Failed to open file\n";
     exit(-1);
   }
   for (int i=0;i<M;i++) {
     for (int j=0;j<P;j++) {
       fout << ret[i*P+j] << endl;
     }
   }
   fout.close();
   #endif

   // fTR(get_run_time());
   {
     TIMER(sort sv);
     sort(ALL(sv));
   }
   begin = clock();
   forN ( k, 1000 ) {
     SortValue& v = sv[k];
     size_t i = v.second.first;
     size_t j = v.second.second;
     model.select_genotype(i);
     model.select_phenotype(j);
     model.fitLM_phenotype_estimate(i*P+j);
   }
   forN ( k, T ) {
     const SortValue& v = sv[k];
     size_t i = v.second.first;
     size_t j = v.second.second;
     model.select_genotype(i);
     model.select_phenotype(j);
     ret[i*P+j] = model.fitLM_phenotype();
   }
   end = clock();
   time_spent = (double)(end-begin)/CLOCKS_PER_SEC;
   cerr << "done final_fits " << time_spent << endl;   
   #ifdef ESTIMATE
   fout.open("venco_final.txt");
   if (!fout) {
     cout << "Failed to open file\n";
     exit(-1);
   }
   for (int i=0;i<M;i++) {
     for (int j=0;j<P;j++) {
       fout << ret[i*P+j] << endl;
     }
   }
   fout.close();
   #endif
#endif // processing_order == by_estimate (the default)
   // MARK START OF ALT ESTIMATES
#if processing_order == by_start
   typedef pair<double, pair<size_t, size_t> > SortValue;
   vector<SortValue> sv; sv.reserve(T);
   SortValue v;
   forN ( i, M ) {
     model.select_genotype(i);
     forN ( j, P ) {
       model.select_phenotype(j);
       double t = model.fitLM_phenotype_start(i*P+j);
       v.first = -t;
       v.second.first = i;
       v.second.second = j;
       sv.push_back(v);
       ret[i*P+j] = t;
     }
   }
   // TR(get_run_time());
   {
     TIMER(sort sv);
     sort(ALL(sv));
   }
   forN ( k, T ) {
     const SortValue& v = sv[k];
     size_t i = v.second.first;
     size_t j = v.second.second;
     model.select_genotype(i);
     model.select_phenotype(j);
     ret[i*P+j] = model.fitLM_phenotype_continue(1);
   }
#endif
#if processing_order == by_estimate_start1
   typedef pair<double, pair<size_t, size_t> > SortValue;
   vector<SortValue> sv; sv.reserve(T);
   SortValue v;
   forN ( i, M ) {
     model.select_genotype(i);
     forN ( j, P ) {
       model.select_phenotype(j);
       double t = model.fitLM_phenotype_estimate(i*P+j);
       v.first = -t;
       v.second.first = i;
       v.second.second = j;
       sv.push_back(v);
       ret[i*P+j] = t;
     }
   }
   // TR(get_run_time());
   {
     TIMER(sort sv);
     sort(ALL(sv));
   }
   multiset<SortValue> sq1;
   size_t k1 = 0;
   forN ( k, T ) {
     for ( size_t T1 = min(T, (k+1)*2); k1 < T1; ++k1 ) {
       SortValue& v = sv[k1];
       size_t i = v.second.first;
       size_t j = v.second.second;
       model.select_genotype(i);
       model.select_phenotype(j);
       double t = model.fitLM_phenotype_start(i*P+j);
       ret[i*P+j] = t;
       v.first = -t;
       sq1.insert(v);
     }
     const SortValue& v = *sq1.begin();
     size_t i = v.second.first;
     size_t j = v.second.second;
     sq1.erase(sq1.begin());
     model.select_genotype(i);
     model.select_phenotype(j);
     //model.getResult(i*P+j);
     ret[i*P+j] = model.fitLM_phenotype_continue(2);
     //ret[i*P+j] = model.getResult(i*P+j);
   }
#endif
#if processing_order == by_estimate_start2
   typedef pair<double, pair<size_t, size_t> > SortValue;
   vector<SortValue> sv; sv.reserve(T);
   SortValue v;
   forN ( i, M ) {
     model.select_genotype(i);
     forN ( j, P ) {
       model.select_phenotype(j);
       double t = model.fitLM_phenotype_estimate(i*P+j);
       v.first = -t;
       v.second.first = i;
       v.second.second = j;
       sv.push_back(v);
       ret[i*P+j] = t;
     }
   }
   // TR(get_run_time());
   {
     TIMER(sort sv);
     sort(ALL(sv));
   }
   multiset<SortValue> sq1, sq2;
   size_t k1 = 0;
   forN ( k, T ) {
     for ( size_t T1 = min(T, (k+1)*2); k1 < T1; ++k1 ) {
       SortValue& v = sv[k1];
       size_t i = v.second.first;
       size_t j = v.second.second;
       model.select_genotype(i);
       model.select_phenotype(j);
       double t = model.fitLM_phenotype_start(i*P+j);
       ret[i*P+j] = t;
       v.first = -t;
       sq1.insert(v);
     }
     while ( !sq1.empty() && sq2.size() < (k+1) ) {
       SortValue v = *sq1.begin(); sq1.erase(sq1.begin());
       size_t i = v.second.first;
       size_t j = v.second.second;
       model.select_genotype(i);
       model.select_phenotype(j);
       //model.getResult(i*P+j);
       double t = model.fitLM_phenotype_step(i*P+j, 1);
       ret[i*P+j] = t;
       v.first = -t;
       sq2.insert(v);
     }
     const SortValue& v = *sq2.begin();
     size_t i = v.second.first;
     size_t j = v.second.second;
     sq2.erase(sq2.begin());
     model.select_genotype(i);
     model.select_phenotype(j);
     //model.getResult(i*P+j);
     ret[i*P+j] = model.fitLM_phenotype_continue(2);
     //ret[i*P+j] = model.getResult(i*P+j);
   } // forN( k, T)
#endif
   // MARK END OF ALT ESTIMATES
 }

 fout.open("venco_final.txt");
   if (!fout) {
     cout << "Failed to open file\n";
     exit(-1);
   }
   for (int i=0;i<M;i++) {
     for (int j=0;j<P;j++) {
       fout << ret[i*P+j] << endl;
     }
   }
   fout.close();

 return ret;
} // GWASSpeedup::computeAssociations

      //==========================================================
      #ifdef HOME_RUN
      inline bool read_string(FILE* f, bool bin, string& s, size_t L)
      {
          if ( !bin ) {
              int c; 
              do {
                  c = getc(f);
              } while ( c < ' ' );
              ungetc(c, f);
          }
          s.resize(L);
          if ( fread(&s[0], L, 1, f) != 1 ) return 0;
          return 1;
      }

      inline void write_string(FILE* f, bool bin, const string& s, size_t L)
      {
          fwrite(s.data(), L, 1, f);
          if ( !bin ) putc('\n', f);
      }

      void read_strings(FILE* f, bool bin, vector<string>& ss, const char* name)
      {
          size_t N, L;
          if ( bin ) {
              if ( fread(&N, sizeof(N), 1, f) != 1 ) error(name);
              if ( fread(&L, sizeof(L), 1, f) != 1 ) error(name);
          }
          else {
              if ( fscanf(f, "%u%u", &N, &L) != 2 ) error(name);
          }
          ss.resize(N);
          forN ( i, N ) {
              if ( !read_string(f, bin, ss[i], L) ) error(name);
          }
      }

      void write_strings(FILE* f, bool bin, const vector<string>& ss)
      {
          size_t N = ss.size(), L = ss[0].size();
          if ( bin ) {
              fwrite(&N, sizeof(N), 1, f);
              fwrite(&L, sizeof(L), 1, f);
          }
          else {
              fprintf(f, "%u %u\n", N, L);
          }
          forN ( i, N ) {
              write_string(f, bin, ss[i], L);
          }
      }

      void read_arr(FILE* f, bool bin, vector<double>& vv, const char* name)
      {
          size_t N, L;
          if ( bin ) {
              if ( fread(&N, sizeof(N), 1, f) != 1 ) error(name);
              if ( fread(&L, sizeof(L), 1, f) != 1 ) error(name);
          }
          else {
              if ( fscanf(f, "%u%u", &N, &L) != 2 ) error(name);
          }
          size_t S = N*L;
          vv.reserve(S);
          if ( bin ) {
              vv.resize(S);
              if ( fread(&vv[0], S*sizeof(vv[0]), 1, f) != 1 ) error(name);
          }
          else {
              forN ( i, S ) {
                  double v;
                  if ( fscanf(f, "%lf", &v) != 1 ) error(name);
                  vv.push_back(v);
              }
          }
      }

      void read_arr(FILE* f, vector<float>& vv, const char* name, size_t S)
      {
          vv.resize(S);
          if ( fread(&vv[0], S*sizeof(vv[0]), 1, f) != 1 ) error(name);
      }

      void write_arr(FILE* f, bool bin, const vector<double>& vv, size_t N)
      {
          size_t L = vv.size()/N;
          if ( bin ) {
              fwrite(&N, sizeof(N), 1, f);
              fwrite(&L, sizeof(L), 1, f);
          }
          else {
              fprintf(f, "%u %u\n", N, L);
          }
          if ( bin ) {
              fwrite(&vv[0], vv.size()*sizeof(vv[0]), 1, f);
          }
          else {
              size_t j = 0;
              forIter ( i, vv ) {
                  fprintf(f, "%.15lg", *i);
                  char c = ' ';
                  if ( ++j == L ) {
                      j = 0;
                      c = '\n';
                  }
                  putc(c, f);
              }
          }
      }

      void write_arr(FILE* f, bool bin, const vector<double>& vv)
      {
          if ( bin ) {
              fwrite(&vv[0], vv.size()*sizeof(vv[0]), 1, f);
          }
          else {
              forIter ( i, vv ) {
                  fprintf(f, "%.6lg\n", *i);
              }
          }
      }

//AH: in PFE environ, this #if evaluates to false, causing missing definition
// later at line 3238.  Appears to be a back-hack for Topcoder's 4.03 gcc?
// force this def by commenting out #if, beware of implications
/* #if (__GNUC__*100+__GNUC_MINOR__*10+__GNUC_PATCHLEVEL__) <= 403 */
      inline uint32_t __builtin_bswap32(uint32_t v)
      {
          return (v>>24)+((v>>8)&0xff00)+((v<<8)&0xff0000)+(v<<24);
      }
/*      #endif */

#ifdef MAIN 
      int main(int argc, char** argv)
      {
	FILE* out = stdout; 
	bool out_bin = 0;
	clock_t begin, end;
	double time_spent;
	/*
          FILE* in = stdin; bool in_bin = 0;
          FILE* ref = 0;

          FILE* save = 0; bool save_bin = 0;
          forNF ( i, 1, size_t(argc) ) {
              if ( strcmp(argv[i], "-s") == 0 && i+1 < size_t(argc) ) {
                  save = fopen(argv[++i], "wb");
                  save_bin = 0;
                  continue;
              }
              if ( strcmp(argv[i], "-sb") == 0 && i+1 < size_t(argc) ) {
                  save = fopen(argv[++i], "wb");
                  save_bin = 1;
                  continue;
              }
              if ( strcmp(argv[i], "-i") == 0 && i+1 < size_t(argc) ) {
                  in = fopen(argv[++i], "rb");
                  in_bin = 0;
                  continue;
              }
              if ( strcmp(argv[i], "-ib") == 0 && i+1 < size_t(argc) ) {
                  in = fopen(argv[++i], "rb"); 
                  in_bin = 1;
                  continue;
              }
              if ( strcmp(argv[i], "-rb") == 0 && i+1 < size_t(argc) ) {
                  //ref = fopen(argv[++i], "rb"); 
                  continue;
              }
              if ( strcmp(argv[i], "-o") == 0 && i+1 < size_t(argc) ) {
                  out = fopen(argv[++i], "wb");
                  out_bin = 0;
                  continue;
              }
              if ( strcmp(argv[i], "-ob") == 0 && i+1 < size_t(argc) ) {
                  out = fopen(argv[++i], "wb");
                  out_bin = 1;
                  continue;
              }
          }
          if ( argc == 1 ) {
              in = fopen("examples/seed1in.bin", "rb");
              in_bin = 1;
          }
	*/
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
	
	GWASSpeedup* obj = 0;
	vector<string> genotypes, phenotypes;
	vector<double> covariates, result;
	time_t timev;
	time_t starttime = time(&timev);
	cerr << "===== Start GWASSpeedup run: " << asctime(localtime(&starttime));
	{
	    /*              {
                  //TIMER(read);
                  read_strings(in, in_bin, phenotypes, "phenotypes");
                  read_strings(in, in_bin, genotypes, "genotypes");
                  read_arr(in, in_bin, covariates, "covariates");
              }
              if ( ref ) {
                  size_t G = genotypes[0].size();
                  size_t P = phenotypes[0].size();
                  read_arr(ref, reference, "reference", G*P);
                  forN ( i, G*P ) {
                      uint32_t* r = (uint32_t*)&reference[i];
                      *r = __builtin_bswap32(*r);
                  }
                  fclose(ref);
              }
              if ( save ) {
                  //TIMER(save);
                  write_strings(save, save_bin, phenotypes);
                  write_strings(save, save_bin, genotypes);
                  write_arr(save, save_bin, covariates, phenotypes.size());
                  fclose(save);
              }
	    */
	  readPlinkData(filestub, genotypes, phenotypes, covariates);
	  time_t end_read_data = time(&timev);
	  cerr << " readPlinkData time : " << end_read_data - starttime << " sec" << endl;
	  obj = new GWASSpeedup();
	  begin = clock();   
	  obj->computeAssociations(phenotypes,genotypes,covariates).swap(result);
	  end = clock();
	  time_spent = (double)(end-begin)/CLOCKS_PER_SEC;
	  cerr << "computeAssociations: " << time_spent << " sec" << endl;
	  assert(result.size() == genotypes[0].size()*phenotypes[0].size());
	  for (int i = 0; i < result.size(); ++ i) {
	    result[i] = chiprobP(result[i], 1.0);
	  }
	  if ( out_bin ) {
	    //TIMER(write);
	    write_arr(out, out_bin, result);
	    fclose(out);
	    out = 0;
	  } 
	} // crazy block
          if ( out ) {
              write_arr(out, out_bin, result);
              fflush(out);
          }
          fflush(stdout);
	  time_t end_main_time = time(&timev);
	  cerr << " end-to-end main() run time : " << end_main_time - starttime << " sec" << endl;
	  return 0;
      } // main
#endif
# endif
