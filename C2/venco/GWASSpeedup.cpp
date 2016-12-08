/**
 * exp_ps : exponential using single precision SSE from fast math library for float.
 */
/**
	@brief fast math library for float
	@author herumi
	@url http://homepage1.nifty.com/herumi/
	@note modified new BSD license
	http://opensource.org/licenses/BSD-3-Clause

	cl /Ox /Ob2 /arch:SSE2 /fp:fast bench.cpp -I../xbyak /EHsc /DNOMINMAX
	g++ -O3 -fomit-frame-pointer -fno-operator-names -march=core2 -mssse3 -mfpmath=sse -ffast-math -fexcess-precision=fast
*/

#include "../common/dcdflib.h"
#include "../common/ipmpar.h"

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

#include "pmmintrin.h"
#include <cstdlib>
#include <cmath>

//{{{
namespace fmath {
  const size_t EXP_TABLE_SIZE = 10;
  union fi {
    float f;
    unsigned int i;
  };
  inline unsigned int mask(int x)
  {
    return (1U << x) - 1;
  }
  template<typename T>
    inline const T* cast_to(const void *p)
    {
      return reinterpret_cast<const T*>(p);
    }
  template<size_t N = EXP_TABLE_SIZE>
    struct ExpVar {
      enum {
        s = N,
        n = 1 << s,
        f88 = 0x42b00000
      };
      float minX[4];
      float maxX[4];
      float a[4];
      float b[4];
      float f1[4];
      unsigned int i127s[4];
      unsigned int mask_s[4];
      unsigned int i7fffffff[4];
      unsigned int tbl[n];
      ExpVar()
      {
        float log_2 = log(2.0);
        for (int i = 0; i < 4; i++) {
          maxX[i] = 88;
          minX[i] = -88;
          a[i] = n / log_2;
          b[i] = log_2 / n;
          f1[i] = 1.0f;
          i127s[i] = 127 << s;
          i7fffffff[i] = 0x7fffffff;
          mask_s[i] = mask(s);
        }

        for (int i = 0; i < n; i++) {
          float y = pow(2.0f, (float)i / n);
          fi fi;
          fi.f = y;
          tbl[i] = fi.i & mask(23);
        }
      }
    };

  template<size_t EXP_N = EXP_TABLE_SIZE>
    struct C {
      static const ExpVar<EXP_N> expVar;
    };

  template<size_t EXP_N>
    __attribute__((aligned(16))) const ExpVar<EXP_N> C<EXP_N>::expVar;

  inline __m128 exp_ps(__m128 x)
  {
    const ExpVar<>& expVar = C<>::expVar;

    __m128i limit = _mm_castps_si128(_mm_and_ps(x, *cast_to<__m128>(expVar.i7fffffff)));
    int over = _mm_movemask_epi8(_mm_cmpgt_epi32(limit, *cast_to<__m128i>(expVar.maxX)));
    if (over) {
      x = _mm_min_ps(x, _mm_load_ps(expVar.maxX));
      x = _mm_max_ps(x, _mm_load_ps(expVar.minX));
    }

    __m128i r = _mm_cvtps_epi32(_mm_mul_ps(x, *cast_to<__m128>(expVar.a)));
    __m128 t = _mm_sub_ps(x, _mm_mul_ps(_mm_cvtepi32_ps(r), *cast_to<__m128>(expVar.b)));
    t = _mm_add_ps(t, *cast_to<__m128>(expVar.f1));

    __m128i v4 = _mm_and_si128(r, *cast_to<__m128i>(expVar.mask_s));
    __m128i u4 = _mm_add_epi32(r, *cast_to<__m128i>(expVar.i127s));
    u4 = _mm_srli_epi32(u4, expVar.s);
    u4 = _mm_slli_epi32(u4, 23);

    unsigned int v0, v1, v2, v3;
    v0 = _mm_cvtsi128_si32(v4);
    v1 = ((int) (unsigned short) __builtin_ia32_vec_ext_v8hi ((__v8hi)(__m128i)(v4), (int)(2)));
    v2 = ((int) (unsigned short) __builtin_ia32_vec_ext_v8hi ((__v8hi)(__m128i)(v4), (int)(4)));
    v3 = ((int) (unsigned short) __builtin_ia32_vec_ext_v8hi ((__v8hi)(__m128i)(v4), (int)(6)));

    __m128 t0, t1, t2, t3;

    t0 = _mm_set_ss(*(const float*)&expVar.tbl[v0]);
    t1 = _mm_set_ss(*(const float*)&expVar.tbl[v1]);
    t2 = _mm_set_ss(*(const float*)&expVar.tbl[v2]);
    t3 = _mm_set_ss(*(const float*)&expVar.tbl[v3]);

    t1 = _mm_movelh_ps(t1, t3);
    t1 = _mm_castsi128_ps(_mm_slli_epi64(_mm_castps_si128(t1), 32));
    t0 = _mm_movelh_ps(t0, t2);
    t0 = _mm_or_ps(t0, t1);

    t0 = _mm_or_ps(t0, _mm_castsi128_ps(u4));
    t = _mm_mul_ps(t, t0);
    return t;
  }

}

//}}}

#include <cstdio>
#include <iostream>
#include <algorithm>
#include <cassert>
#include <string>
#include <vector>
#include <sstream>
#include <sys/time.h>
#include <cstring>
#include <cmath>
#include <iostream>
#include <sstream>
#include "GWASSpeedup.h"
#include "../common/read_bed.h"

using namespace std;

double getTime()
{
    timeval tv;
    gettimeofday(&tv, 0);
    return tv.tv_sec + tv.tv_usec * 1e-6;
}

//#define PROFILING

//{{{

//-------------------------
// Code imported from venco
//
#define forIter(I,C) for(typeof((C).end()) I=(C).begin(); I!=(C).end(); ++I) 
      typedef float Double;

      template<class V>struct aligned_alloc1{
		  typedef V value_type;
          typedef size_t size_type;
          typedef ptrdiff_t difference_type;
          typedef V*pointer;
          typedef const V*const_pointer;
          typedef V&reference;
          typedef const V&const_reference;
          template<class U>struct rebind{typedef aligned_alloc1<U>other;
          };
          aligned_alloc1() throw(){}
          aligned_alloc1(const aligned_alloc1&) throw(){}
          template<class U>aligned_alloc1(const aligned_alloc1<U>&) throw(){}
          ~aligned_alloc1()throw(){}
          V*allocate(size_t c,const void* =0){
			  return(V*)_mm_malloc(sizeof(V)*c,16);
			}
          void deallocate(V* p,size_t ){_mm_free(p);}
          size_t max_size()const throw(){return ~size_t(0)/sizeof(V);}
          void construct(V* p,const V&v){::new((void *)p)V(v);}
          void destroy(V* p){p->~V();}
      };
      template<class V,class U>inline
      bool operator==(const aligned_alloc1<V>&,const aligned_alloc1<U>&){
		  return true;
		}
      template<typename V,typename U>inline
      bool operator!=(const aligned_alloc1<V>&,const aligned_alloc1<U>&){
		  return false;
		}

      typedef vector<Double, aligned_alloc1<Double> > vector_t;
      typedef vector<vector_t> matrix_t;
      
      //resize the matrix to row x col, and clear to zero.
      inline void sizeMatrix(matrix_t &matrix,size_t row,size_t col)
      {
          matrix.resize(row);
          forIter ( i, matrix ) i->resize(col);
      }

//
// end of venco's code
// 

#ifdef PROFILING 
__attribute__((noinline))
#endif
inline void logistic_sse(float *vect, int n) {
  // Equivalent C++ code
  //  for(int i = 0; i < size; i++) vect[i] = 1. / (1 + exp(-vect[i]));
  __m128 zero = _mm_setzero_ps();
  __m128 one = _mm_set1_ps(1.);
  for(int i = 0; i < n; i+= 4) {
    __m128 a = _mm_load_ps(vect + i);
    a = _mm_sub_ps(zero, a);
    a = fmath::exp_ps(a);
    a = _mm_add_ps(a, one);
    a = _mm_div_ps(one, a);
    _mm_store_ps(vect+i, a);
  }
}



#ifdef PROFILING 
__attribute__((noinline))
#endif
inline void compute_V_and_P_minus_Y(float* p, float* v, const float* y, int n) {
  //  equivalent C++ code
  /* for(int i = 0; i < n; i++) {
       v[i] = p[i]*(1-p[i]);
       p[i] -= y[i];
     }*/

  __m128 _one = _mm_set1_ps(1.);
  for(int i = 0; i < n; i += 4) {
    __m128 _p = _mm_load_ps(p + i);
    __m128 _one_minus_p = _mm_sub_ps(_one, _p);
    _mm_store_ps(v+i, _mm_mul_ps(_p, _one_minus_p));
    __m128 _y = _mm_load_ps(y + i);
    _mm_store_ps(p+i, _mm_sub_ps(_p, _y));
  }
}

//~ template<int d>
#ifdef PROFILING 
__attribute__((noinline))
#endif
//~ inline void mult_tmatrix_nxd_vect_d(float** tm, const float* vect, float* dest, int n) {
inline void mult_tmatrix_nxd_vect_d(int d, float** tm, const float* vect, float* dest, int n) {
// equivalent C++ code
/*for(int i = 0; i < n; i++) {
    dest[i] = 0.;
    for(int j = 0; j < d; j++) 
      dest[i] += tm[j][i] * vect[j];
  }*/

  __m128 w1 = _mm_load1_ps(vect);
  __m128 w2 = _mm_load1_ps(vect+1);
  __m128 w3 = _mm_load1_ps(vect+2);
  __m128 w4 = _mm_load1_ps(vect+3);
  for(int i = 0; i < n; i += 4) {
    __m128 r1 = _mm_load_ps(tm[0] + i);
    __m128 r2 = _mm_load_ps(tm[1] + i);
    __m128 r3 = _mm_load_ps(tm[2] + i);
    __m128 r4 = _mm_load_ps(tm[3] + i);
    r1 = _mm_mul_ps(r1, w1);
    r2 = _mm_mul_ps(r2, w2);
    r3 = _mm_mul_ps(r3, w3);
    r4 = _mm_mul_ps(r4, w4);
    r1 = _mm_add_ps(r1, r2);
    r3 = _mm_add_ps(r3, r4);
    r1 = _mm_add_ps(r1, r3);
    _mm_store_ps(dest + i, r1);
  }

  int k = 4;
  for(; k + 3 < d; k += 4) {
    w1 = _mm_load1_ps(vect+k);
    w2 = _mm_load1_ps(vect+k+1);
    w3 = _mm_load1_ps(vect+k+2);
    w4 = _mm_load1_ps(vect+k+3);

    for(int i = 0; i < n; i += 4) {
      __m128 r1 = _mm_load_ps(tm[k] + i);
      __m128 r2 = _mm_load_ps(tm[k+1] + i);
      __m128 r3 = _mm_load_ps(tm[k+2] + i);
      __m128 r4 = _mm_load_ps(tm[k+3] + i);

      r1 = _mm_mul_ps(r1, w1);
      r2 = _mm_mul_ps(r2, w2);
      r3 = _mm_mul_ps(r3, w3);
      r4 = _mm_mul_ps(r4, w4);
      r1 = _mm_add_ps(r1, r2);
      r3 = _mm_add_ps(r3, r4);
      r1 = _mm_add_ps(r1, r3);
      r1 = _mm_add_ps(r1, _mm_load_ps(dest + i));
      _mm_store_ps(dest + i, r1);
    }
  }

  switch(d%4) {
    case 3:
      w1 = _mm_load1_ps(vect+k);
      w2 = _mm_load1_ps(vect+k+1);
      w3 = _mm_load1_ps(vect+k+2);

      for(int i = 0; i < n; i += 4) {
        __m128 r1 = _mm_load_ps(tm[k] + i);
        __m128 r2 = _mm_load_ps(tm[k+1] + i);
        __m128 r3 = _mm_load_ps(tm[k+2] + i);

        r1 = _mm_mul_ps(r1, w1);
        r2 = _mm_mul_ps(r2, w2);
        r3 = _mm_mul_ps(r3, w3);
        r1 = _mm_add_ps(r1, r2);
        r3 = _mm_add_ps(r3, _mm_load_ps(dest + i));
        r1 = _mm_add_ps(r1, r3);
        _mm_store_ps(dest + i, r1);
      }
      break;
    case 2:
      w1 = _mm_load1_ps(vect+k);
      w2 = _mm_load1_ps(vect+k+1);    
      for(int i = 0; i < n; i += 4) {
        __m128 r1 = _mm_load_ps(tm[k] + i);
        __m128 r2 = _mm_load_ps(tm[k+1] + i);
        r1 = _mm_mul_ps(r1, w1);
        r2 = _mm_mul_ps(r2, w2);
        r1 = _mm_add_ps(r1, r2);
        r1 = _mm_add_ps(r1, _mm_load_ps(dest + i));
        _mm_store_ps(dest + i, r1);
      }
      break;
    case 1:
      w1 = _mm_load1_ps(vect+k);
      for(int i = 0; i < n; i += 4) {
        __m128 r1 = _mm_load_ps(tm[k] + i);
        r1 = _mm_mul_ps(r1, w1);
        r1 = _mm_add_ps(r1, _mm_load_ps(dest + i));
        _mm_store_ps(dest + i, r1);
      }
  }
}




//~ template<int d>
#ifdef PROFILING 
__attribute__((noinline))
#endif
//~ inline void mult_matrix_dxn_vect_n(float** m, const float* vect, float dest[12], int n) {
//~ inline void mult_matrix_dxn_vect_n(int d, float** m, const float* vect, float dest[12], int n) {
inline void mult_matrix_dxn_vect_n(int d, float** m, const float* vect, float *dest, int n) {
  // equivalent C++ code
  /* for(int i = 0; i < n; i++) 
       for(int j = 0; j < d; j++) 
         dest[j] = m[j][i] * vect[i]; */

  int k = 0;
  for(; k+3 < d; k+=4) {
    __m128 s1 = _mm_setzero_ps();
    __m128 s2 = _mm_setzero_ps();
    __m128 s3 = _mm_setzero_ps();
    __m128 s4 = _mm_setzero_ps();
    for(int i = 0; i < n; i += 4) {
      __m128 v = _mm_load_ps(vect + i);
      __m128 a1 = _mm_load_ps(m[k] + i);
      __m128 a2 = _mm_load_ps(m[k+1] + i);
      __m128 a3 = _mm_load_ps(m[k+2] + i);
      __m128 a4 = _mm_load_ps(m[k+3] + i);
      a1 = _mm_mul_ps(a1, v);
      a2 = _mm_mul_ps(a2, v);
      a3 = _mm_mul_ps(a3, v);
      a4 = _mm_mul_ps(a4, v);
      s1 = _mm_add_ps(s1, a1);
      s2 = _mm_add_ps(s2, a2);
      s3 = _mm_add_ps(s3, a3);
      s4 = _mm_add_ps(s4, a4);
    }
    s1 = _mm_hadd_ps(s1, s2);
    s3 = _mm_hadd_ps(s3, s4);
    s1 = _mm_hadd_ps(s1, s3);
    _mm_store_ps(dest + k, s1);
  }

  __m128 s1 = _mm_setzero_ps();
  __m128 s2 = _mm_setzero_ps();
  __m128 s3 = _mm_setzero_ps();

  switch(d%4) {
    case 3:
      for(int i = 0; i < n; i += 4) {
        __m128 v = _mm_load_ps(vect + i);
        __m128 a1 = _mm_load_ps(m[k] + i);
        __m128 a2 = _mm_load_ps(m[k+1] + i);
        __m128 a3 = _mm_load_ps(m[k+2] + i);
        a1 = _mm_mul_ps(a1, v);
        a2 = _mm_mul_ps(a2, v);
        a3 = _mm_mul_ps(a3, v);
        s1 = _mm_add_ps(s1, a1);
        s2 = _mm_add_ps(s2, a2);
        s3 = _mm_add_ps(s3, a3);
      }
      s1 = _mm_hadd_ps(s1, s2);
      s3 = _mm_hadd_ps(s3, s3);
      s1 = _mm_hadd_ps(s1, s3);
      _mm_store_ps(dest + k, s1);
      break;
    case 2:
      for(int i = 0; i < n; i += 4) {
        __m128 v = _mm_load_ps(vect + i);
        __m128 a1 = _mm_load_ps(m[k] + i);
        __m128 a2 = _mm_load_ps(m[k+1] + i);
        a1 = _mm_mul_ps(a1, v);
        a2 = _mm_mul_ps(a2, v);
        s1 = _mm_add_ps(s1, a1);
        s2 = _mm_add_ps(s2, a2);
      }
      s1 = _mm_hadd_ps(s1, s2);
      s1 = _mm_hadd_ps(s1, s1);
      _mm_store_ps(dest + k, s1);
      break;
    case 1:
      for(int i = 0; i < n; i += 4) {
        __m128 v = _mm_load_ps(vect + i);
        __m128 a1 = _mm_load_ps(m[k] + i);
        a1 = _mm_mul_ps(a1, v);
        s1 = _mm_add_ps(s1, a1);
      }
      s1 = _mm_hadd_ps(s1, s1);
      s1 = _mm_hadd_ps(s1, s1);
      _mm_store_ps(dest + k, s1);
      break;
  }
}

#ifdef PROFILING 
__attribute__((noinline))
#endif
inline float dot_product(float* v1, float* v2, int n) {
  // Equivalent C++ code
  /* float r = 0;
     for(int i = 0; i < n; i++) r += v1[i]*v2[i];
     return r;*/
  __m128 sum = _mm_setzero_ps();
  for(int i = 0; i < n; i += 4) {
    __m128 a = _mm_load_ps(v1 + i);
    __m128 b = _mm_load_ps(v2 + i);
    sum = _mm_add_ps(sum, _mm_mul_ps(a, b));
  }
  sum = _mm_hadd_ps(sum, sum);
  sum = _mm_hadd_ps(sum, sum);
  float ret;
  _mm_store_ss(&ret, sum);
  return ret;
}


#ifdef PROFILING 
__attribute__((noinline))
#endif
inline  float triple_product(const float* v1, const float* v2, const float* v3, int n) {
  // Equivalent C++ code
  /* float r = 0;
     for(int i = 0; i < n; i++) r += v1[i]*v2[i]*v3[i];
     return r;*/

  __m128 sum = _mm_setzero_ps();
  for(int i = 0; i < n; i += 4) {
    __m128 a = _mm_load_ps(v1 + i);
    __m128 b = _mm_load_ps(v2 + i);
    __m128 c = _mm_load_ps(v3 + i);
    sum = _mm_add_ps(sum, _mm_mul_ps(_mm_mul_ps(a, b), c));
  }
  sum = _mm_hadd_ps(sum, sum);
  sum = _mm_hadd_ps(sum, sum);
  float ret;
  _mm_store_ss(&ret, sum);

  return ret;
}

#ifdef PROFILING 
__attribute__((noinline))
#endif
inline void compute_two_diag_triple_product(const float* a, const float* b, const float* v, float& raa, float& rab, float& rbb, int n) {
  /* equivalent C++ code
     raa = rab = rbb = 0.;
     for(int i = 0; i < n; i++) {
     raa += a[i]*a[i]*v[i];
     rab += a[i]*b[i]*v[i];
     rbb += b[i]*b[i]*v[i];
  }*/

  __m128 saa = _mm_setzero_ps();
  __m128 sab = _mm_setzero_ps();
  __m128 sbb = _mm_setzero_ps();
  for(int i = 0; i < n; i += 4) {
    __m128 _v = _mm_load_ps(v + i);
    __m128 _a = _mm_load_ps(a + i);
    __m128 _b = _mm_load_ps(b + i);
    __m128 _av = _mm_mul_ps(_a, _v);
    __m128 _bv = _mm_mul_ps(_b, _v);
    saa = _mm_add_ps(saa, _mm_mul_ps(_a, _av));
    sab = _mm_add_ps(sab, _mm_mul_ps(_a, _bv));
    sbb = _mm_add_ps(sbb, _mm_mul_ps(_b, _bv));
  }
  saa = _mm_hadd_ps(saa, sab);
  sbb = _mm_hadd_ps(sbb, sbb);
  saa = _mm_hadd_ps(saa, sbb);
  float  __attribute__ ((aligned (16))) ret[4];
  _mm_store_ps(ret, saa);
  raa = ret[0];
  rab = ret[1];
  rbb = ret[2];
}

#ifdef PROFILING 
__attribute__((noinline))
#endif
inline void compute_three_triple_product(const float* b, const float* a1, const float* a2, const float* a3, const float* v, float& r1, float& r2, float& r3, int n) {
  // equivalent C++ code
  /*r1 = r2 = r3 = 0.;
    for(int i = 0; i < n; i++) {
    r1 += a1[i]*b[i]*v[i];
    r2 += a2[i]*b[i]*v[i];
    r3 += a3[i]*b[i]*v[i];
    }*/
  __m128 s1 = _mm_setzero_ps();
  __m128 s2 = _mm_setzero_ps();
  __m128 s3 = _mm_setzero_ps();
  for(int i = 0; i < n; i += 4) {
    __m128 _a1 = _mm_load_ps(a1 + i);
    __m128 _a2 = _mm_load_ps(a2 + i);
    __m128 _a3 = _mm_load_ps(a3 + i);
    __m128 _v = _mm_load_ps(v + i);
    __m128 _b = _mm_load_ps(b + i);
    _b = _mm_mul_ps(_b, _v);
    s1 = _mm_add_ps(s1, _mm_mul_ps(_a1, _b));
    s2 = _mm_add_ps(s2, _mm_mul_ps(_a2, _b));
    s3 = _mm_add_ps(s3, _mm_mul_ps(_a3, _b));
  }
  s1 = _mm_hadd_ps(s1, s2);
  s3 = _mm_hadd_ps(s3, s3);
  s1 = _mm_hadd_ps(s1, s3);
  float  __attribute__ ((aligned (16))) ret[4];
  _mm_store_ps(ret, s1);
  r1 = ret[0];
  r2 = ret[1];
  r3 = ret[2];
}

#ifdef PROFILING 
__attribute__((noinline))
#endif
inline void compute_two_plus_one_triple_product(const float* b, const float* a1, const float* a2, const float* v, float& r1, float& r2, float& r3, int n) {
  // equivalent C++ code
  /*  r1 = r2 = r3 = 0.;
      for(int i = 0; i < n; i++) {
      r1 += b[i]*b[i]*v[i];
      r2 += a1[i]*b[i]*v[i];
      r3 += a2[i]*b[i]*v[i];
      }*/
  __m128 s1 = _mm_setzero_ps();
  __m128 s2 = _mm_setzero_ps();
  __m128 s3 = _mm_setzero_ps();
  for(int i = 0; i < n; i += 4) {
    __m128 _a1 = _mm_load_ps(a1 + i);
    __m128 _a2 = _mm_load_ps(a2 + i);
    __m128 _b = _mm_load_ps(b + i);
    __m128 _v = _mm_load_ps(v + i);
    __m128 _bv = _mm_mul_ps(_b, _v);
    s1 = _mm_add_ps(s1, _mm_mul_ps(_b, _bv));
    s2 = _mm_add_ps(s2, _mm_mul_ps(_a1, _bv));
    s3 = _mm_add_ps(s3, _mm_mul_ps(_a2, _bv));
  }
  s1 = _mm_hadd_ps(s1, s2);
  s3 = _mm_hadd_ps(s3, s3);
  s1 = _mm_hadd_ps(s1, s3);
  float  __attribute__ ((aligned (16))) ret[4];
  _mm_store_ps(ret, s1);
  r1 = ret[0];
  r2 = ret[1];
  r3 = ret[2];
}


inline void compute_hessian(int d, float* const * m, const float* v, matrix_t &dest, int n) {
  // compute all the elements of the hessian 3 by 3 (reading only 3 or 4 vector of size n at each time.
  for(int j = 0; j+3 < d; j += 3) {
    compute_two_diag_triple_product(m[j], m[j+1], v, dest[j][j], dest[j+1][j], dest[j+1][j+1], n);
    compute_two_plus_one_triple_product(m[j+2], m[j+1], m[j], v, dest[j+2][j+2], dest[j+2][j+1], dest[j+2][j], n);
    for(int i = j+3; i < d; i++) compute_three_triple_product(m[i], m[j], m[j+1], m[j+2], v, dest[i][j], dest[i][j+1], dest[i][j+2], n);
  }
  switch(d%3) {
    case 0: compute_two_plus_one_triple_product(m[d-3], m[d-2], m[d-1], v, dest[d-3][d-3], dest[d-2][d-3], dest[d-1][d-3], n);
    case 2: compute_two_diag_triple_product(m[d-2], m[d-1], v, dest[d-2][d-2], dest[d-1][d-2], dest[d-1][d-1], n); break;
    case 1: dest[d-1][d-1] = triple_product(m[d-1], m[d-1], v, n);
  }

}  

inline void compute_last_row_hessian(int d, float* const * m, const float* v, matrix_t &dest, int n) {
  // only compute the last row (an column) of the hessian
  compute_two_plus_one_triple_product(m[d-1], m[d-2], m[d-3], v, dest[d-1][d-1], dest[d-1][d-2], dest[d-1][d-3], n);
  int j = d-4;
  for(; j - 2 >= 0; j-=3) compute_three_triple_product(m[d-1], m[j], m[j-1], m[j-2], v,  dest[d-1][j], dest[d-1][j-1], dest[d-1][j-2], n);
  for(; j >= 0; j--) dest[d-1][j] = triple_product(m[d-1], m[j], v, n);
}

//}}}


const double TIME_LIMIT = 10.;
timeval BEGIN;
inline void start() {gettimeofday(&BEGIN, NULL);}
inline double time() {
  timeval NOW; gettimeofday(&NOW, NULL);
  return NOW.tv_sec - BEGIN.tv_sec + (NOW.tv_usec - BEGIN.tv_usec)/1000000.;
}

vector<double> RET;

int N; // number of individuals
int M; // number of markers
int P; // number of phenotypes
int C; // number of covariate columns
int D; // regression space dimension D = C+2

struct TestCase {
  //~ float  __attribute__ ((aligned (16))) coef[12];
  vector_t coef;
  float** X;
  float* Y;
  float W;
  int n;
  int iteration;
  bool restart;
  float delta_coef;
  float min_delta_coef;
  TestCase() : X(0), Y(0), W(0), n(N), iteration(0), restart(false), delta_coef(0), min_delta_coef(1e6) {}
  inline bool bad() {return (delta_coef != delta_coef || (iteration > 3 && delta_coef > 100 && delta_coef > 2*min_delta_coef));}
};
  inline bool operator<(const TestCase& T1, const TestCase& T2) { return T1.W < T2.W;}


matrix_t _H;
matrix_t _L;
vector_t _LD;
vector_t _grad;
vector_t _ncoef;

vector_t _P;
vector_t _V;
matrix_t __Y;
vector_t __one;
vector_t __zero;
matrix_t __cov;
matrix_t __marker;

vector<matrix_t> start_H;
matrix_t start_P;
matrix_t start_V;
matrix_t start_grad;
matrix_t start_coef;

vector<int> missing_marker;
float* ___X;
float** __X;

inline void LDL_cholesky_decompsition(int n, const matrix_t &A,  matrix_t &L, float *LD)
{
  // Compute a Cholesky factorization A=LL* (Cholesky–Banachiewicz algorithm)
  // Input A, Output L, dimension n, A and L can be the same array
  // see = http://en.wikipedia.org/wiki/Cholesky_decomposition
  for(int j=0; j<n; ++j) {
    float diag_val = A[j][j];
    for(int k=0; k<j; ++k)
      diag_val -= L[j][k]*L[j][k]*LD[k];
      
    L[j][j] = 1.0f;
    if(diag_val<0) diag_val=1e-6;
    LD[j] = diag_val;
    diag_val = 1 / diag_val;
    
    for(int i=j+1; i<n; ++i) {
      float val = A[i][j];
      for(int k=0; k<j; ++k) val -= L[i][k]*L[j][k]*LD[k];
      L[i][j] = val * diag_val;
    }
  }
}

inline void LDL_solve_linear_system(int n, const  matrix_t &L, const float *LD, 
	float *Y, float *X)
	//~ float Y[12], float X[12])
{
  // Solve Mult(LL*,X)=Y with L lower triangular cholesky factorization
  // dimension n
  // Solving Mult(L, Xtmp1)=Y
  for(int i=0; i<n; ++i) {
	float s = Y[i];
    for(int j=0; j<i; ++j) s -= L[i][j]*X[j];
    X[i] = s;
  }
  // Solving Mult(DL*, X)=Xtmp1
  for(int i=n-1;i>=0;--i) {
	float s = X[i];  
    s /= LD[i];
    for(int j=i+1; j<n; ++j) s -= L[j][i]*X[j];
    X[i] = s;
  }
}

inline double LDL_compute_wald(int n, const  matrix_t &L, const float *LD)
{
  // Solve Mult(LL*,X)=Y with L lower triangular cholesky factorization
  // dimension n
  // Solving Mult(L, Xtmp1)=Y
	float *X = new float[n]; 

  for(int i=0; i<n-1; ++i) {
	float s = 0;
    for(int j=0; j<i; ++j) s -= L[i][j]*X[j];
    X[i] = s;
  }

	float res = 1;
	for(int j=0; j<n-1; ++j) res -= L[n-1][j]*X[j];

	delete X;
  return ( res / LD[n-1] );
}

void update_coef(int d, TestCase &T) {
  T.iteration++;

  LDL_cholesky_decompsition(d, _H, _L, &_LD[0]);

  memset(&_ncoef[0], 0, d*sizeof(float));
   
  LDL_solve_linear_system(d, _L, &_LD[0], &_grad[0], &_ncoef[0]);

  // Update coefficients, and check for 
  // convergence
  T.delta_coef = 0;

  for (int j=0; j< d; j++) T.delta_coef += abs(_ncoef[j]);

  float alpha = 1.;

  if(T.iteration > 7 && T.delta_coef > 10) alpha = 0.5;
  for (int j=0; j< d; j++) T.coef[j] -= alpha*_ncoef[j];
  if(T.delta_coef < T.min_delta_coef) T.min_delta_coef = T.delta_coef;
  
  float nW = LDL_compute_wald(d,_L, &_LD[0]);

  if(nW <= 0 || nW != nW || T.coef[d-1] != T.coef[d-1]) nW = 0;
  else nW = T.coef[d-1]*T.coef[d-1] / nW;
  T.W = nW;

  if(T.bad()) {
    if(T.restart) T.W = 0.;
    else T.W = 100.;
    memset(&T.coef[0], 0, d*sizeof(float));
    T.restart = true;
    T.iteration = 0;
    T.delta_coef = 1.;
    T.min_delta_coef = 1e6;
    return;
  }
}

void _newton_iteration(int d, TestCase &T) {
  // P[i] = SUM_j coef[j] * _X[i][j]
  mult_tmatrix_nxd_vect_d(d,T.X, &T.coef[0], &_P[0], T.n);

  // P[i] = 1 / (1 + exp(-P[i]))
  logistic_sse(&_P[0], T.n);

  //for(int i = 0; i < T.n; i++) {
  //   _V[i] = _P[i]*(1-_P[i]);
  //   _P[i] -= Y[i];
  // }
  compute_V_and_P_minus_Y(&_P[0], &_V[0], T.Y, T.n);

  compute_hessian(d,T.X, &_V[0], _H, T.n);

  mult_matrix_dxn_vect_n(d,T.X, &_P[0], &_grad[0], T.n);

  update_coef(d,T);
}


void _newton_first_iteration(int d, TestCase &T, const float* P, const float* V, 
	const matrix_t &H, const float *grad, const vector<int> &missing)  {
  memcpy(&_P[0], P, T.n*sizeof(float));
  memcpy(&_V[0], V, T.n*sizeof(float));
  for(int j=0; j<d; ++j)
	memcpy(&_H[j][0], &H[j][0], d*sizeof(float));
  memcpy(&_grad[0], grad, d*sizeof(float));

  float* X[d-1];
  X[0] = &__one[0];
  for(int i = 0; i < C; i++) X[1+i] = &__cov[i][0];

	int ms = missing.size();
  for(int k = 0; k < ms; k++) { // remove the contribution of samples that where removed in this regression but used in the global initial regression
    int m = missing[k];
    for(int i = 0; i < d-1; i++) 
      for(int j = 0; j <= i; j++) 
        _H[i][j] -= X[i][m]*X[j][m]*V[m];

    for(int i = 0; i < d-1; i++) 
      _grad[i] -= X[i][m]*_P[m];
  }

  compute_last_row_hessian(d,T.X, &_V[0], _H, T.n);  // update hessian with new column
  _grad[d-1] = dot_product(T.X[d-1], &_P[0], T.n);      // update gradient with new column

  update_coef(d,T);
}



float* get_Y(int phenotype) {
  return &__Y[phenotype][0];
}

void set_X(int marker) {
	
	int L = 4*((N+3)/4);
  memcpy(__X[C+1], &__marker[marker][0], L*sizeof(float));

	missing_marker.clear();
  for(int i = 0; i < N; i++) 
    if(__marker[marker][i] < 0) {
      missing_marker.push_back(i);
      for(int j = 0; j < D; j++)
        __X[j][i] = 0.;
    }
}

void reset_X(int marker) {
	forIter( m, missing_marker) {
		__X[0][*m] = 1.0f;
		for(int i = 0; i < C; i++) 
			__X[i+1][*m] = __cov[i][*m];
	}
}


void init_data(const vector<string> &phenotype, const vector<string> &genotype, const vector<double> &covariate) {
  N = phenotype.size();
  C = covariate.size() / N;
  P = phenotype[0].size();
  M = genotype.size();
  D = C+2;

	int L = 4*((N+3)/4);

  _V.resize(L);
  _P.resize(L);

  RET.resize(M*P);
	__Y.resize(P);
  for(int j = 0; j < P; j++) {
    __Y[j].resize(L, float(0.0f - '0') );
    for(int i = 0; i < N; i++)
      __Y[j][i] += phenotype[i][j];
     for(int i = N; i < L; i++)
		__Y[j][i] = 0.0f;
  }
  
  int t0 = getTime();
  __marker.resize(M);
  float float2 = float('2');
  for(int j = 0; j < M; j++) {
    __marker[j].resize(L, float2);
    for(int i = 0; i < N; i++)
      __marker[j][i] -= genotype[j][i];
    for(int i = N; i < L; i++)
      __marker[j][i] = 0.0f;
  }	
  cerr << "	genotype load time: " << getTime() - t0 << endl;
  
	__cov.resize(C);
	for(int j = 0; j < C; j++)
		__cov[j].resize(L);
	
	int k=0;
	for(int i = 0; i < N; i++)
		for(int j = 0; j < C; j++)
		__cov[j][i] = covariate[k++];

	for(int j = 0; j < C; j++)
		for(int i = N; i < L; i++)
		__cov[j][i] = 0.0f;
		
	__one.resize(L, 1.f);
	for(int i = N; i < L; i++)
		__one[i] = 0.f;

	__zero.resize(L, 0.);
    
  __X = new float*[D];
  ___X = (float*) _mm_malloc(D*L*sizeof(float), 16);
  for(int i = 0; i < D; i++)
    __X[i] = ___X + L*i;

  memcpy(__X[0], &__one[0], L*sizeof(float));
  for(int i = 0; i < C; i++) 
    memcpy(__X[1+i], &__cov[i][0], L*sizeof(float));


	start_P.resize(P);
	start_V.resize(P);
	start_coef.resize(P);
	start_grad.resize(P);
	start_H.resize(P);	
		
	_ncoef.assign(D,0.0f);
	_grad.assign(4*((D+3)/4), 0.0f);
	_H.assign(D, vector_t(D,0.0f));
	_L.assign(D, vector_t(D,0.0f));
	_LD.assign(D, 0.0f);
	
}

int cmp(const void *T1, const void *T2) {
  if (((TestCase*)T1)->W < ((TestCase*)T2)->W) return -1;
  else if (((TestCase*)T1)->W > ((TestCase*)T2)->W) return 1;
  else return 0;
}

// compute regression for each phenotype with covariate but without marker
// The solution will be the initial step for the regression with the same phenotype + one the additionnal column.
void init_global_case() {
  for(int i = 0; i < P; i++) {
    TestCase T;
    T.coef.assign(D,0);
    T.X = new float*[D];
    T.X[0] = &__one[0];
    for(int j = 0; j < C; j++) T.X[1+j] = &__cov[j][0];
    T.Y = get_Y(i);
    for(int j = 0; j < 10; j++) {
      //~ newton_iteration(T, D-1);
      _newton_iteration(D-1, T);
      if(T.delta_coef < 1e-6) break; // TODO PARAM
    }

	start_coef[i] = T.coef;
    // save intermediate computation that can be reused in the first newton step of each marker regression for the same phenotype
    _newton_iteration(D-1, T);

	start_H[i] = _H;		
    start_P[i] = _P;
    start_V[i] = _V;
    start_grad[i] = _grad;
  }
}

//inline bool cmp(const int& i1, const int& i2) {return TC[i1] < TC[i2];}

void compute(bool approximate) {
	int index = 0;	
  for (int marker = 0; marker < M; marker++) {
	  set_X(marker);		
    for (int phenotype = 0; phenotype < P; phenotype++) {		
      TestCase T;
      T.X = __X;
      T.Y = get_Y(phenotype);
      T.coef = start_coef[phenotype];
      
      _newton_first_iteration(D, T, &start_P[phenotype][0], &start_V[phenotype][0], 
		start_H[phenotype], &start_grad[phenotype][0], missing_marker);
      while(!approximate) {
		  float old_W = T.W;
		  _newton_iteration(D, T);

	//         cerr << T.index << " " << T.iteration << " " << T.W << " " << T.coef[D-1] << " " << T.delta_coef << endl;

		  if( T.iteration >= 15 || ( T.iteration >= 8 && abs(1. - T.delta_coef) < 1e-3) ||
			  (T.delta_coef < 1e-3 && abs(T.W - old_W) < 0.1)  ||
			  (T.W < 1e-3 - 1e-6 && T.iteration > 2 && T.delta_coef < 10.) || 
			  (T.W == 0 && T.restart)
			) {
	//          cerr << T.index << " OK" << endl;
			break;
		  }
		}
		if((T.iteration >= 8 && abs(1. - T.delta_coef) < 1e-3) || T.W != T.W || 
			T.W < 1e-3 - 1e-6 || T.W > 1e8) T.W = 1e-3 - 1e-6;
			RET[index] = T.W;
		index++;	      
    }
  reset_X(marker);
  }
}



vector<double> GWASSpeedup::computeAssociations(const vector<string> &phenotype, 
	const vector<string> &genotype, const vector<double> &covariate){
  int t0 = getTime();
  start();
  init_data(phenotype, genotype, covariate);
  cerr << "	Sub GWASSpeedup data load time : " << getTime() - t0 << endl;
  init_global_case();
  compute(__approximate);
  cerr << "Computed statistics ";
  if (__approximate) {
    cerr << "(using approximation)" << endl;
  } else {
    cerr << "(using complete solution)" << endl;
  }
  cerr << "total time : " << getTime() - t0 << endl;
  return RET;
}

//==========================================================
#ifdef HOME_RUN
#ifdef MAIN

int main(int argc, char** argv)
{
	int c;
	string filestub;
	int firstMarker = 0;
	int lastMarker = 99;
	bool approximate = false;
	while ((c = getopt (argc, argv, "ab:f:l:")) != -1)
	  switch (c)
	    {
	    case 'a':
	      approximate = true;
	      break;
	    case 'b':
	      filestub = string(optarg);
	      break;
	    case 'f':
	      firstMarker = atoi(optarg);
	      break;
	    case 'l':
	      lastMarker = atoi(optarg);
	      break;
	    default:
	      abort();
	    }
	
	GWASSpeedup* obj = 0;
	vector<string> genotypes, phenotypes;
	vector<double> covariates, result;
	vector<string> locus_name;
	time_t timev;
	time_t starttime = time(&timev);
	cerr << "===== Start GWASSpeedup run: " << asctime(localtime(&starttime));
	{
	  readPlinkData(filestub, genotypes, phenotypes, covariates, locus_name, firstMarker, lastMarker);
	  time_t end_read_data = time(&timev);
	  cerr << " readPlinkData time : " << end_read_data - starttime << " sec" << endl;
	  obj = new GWASSpeedup(approximate);
	  time_t start_compute_associations = time(&timev);
	  obj->computeAssociations(phenotypes,genotypes,covariates).swap(result);
	  time_t end_compute_associations = time(&timev);
	  cerr << "computeAssociations: " << (end_compute_associations-start_compute_associations) << " sec" << endl;
	  for (int i = 0; i < result.size(); ++ i) {
	    result[i] = chiprobP(result[i], 1.0);
	  }
	  assert(result.size() == genotypes.size()*phenotypes[0].size());
	  assert(result.size() == locus_name.size());

	  stringstream buf;
	  buf << std::scientific;
	  buf.precision(4);
	  int ms = buf.str().max_size()/100;
	  int rs = result.size();
	  for(int i=0;i<rs;++i) {
	    buf << locus_name[i] << "\t" << result[i] << "\n";
	    if(i==ms || i==rs-1) {
	      cout << buf.str();
	      buf.str("");
	    }
	  }
	}
	time_t end_main_time = time(&timev);
	cerr << " end-to-end main() run time : " << end_main_time - starttime << " sec" << endl;
	return 0;
} // main
#endif
#endif

