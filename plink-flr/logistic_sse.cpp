

//////////////////////////////////////////////////////////////////
//                                                              //
//           PLINK (c) 2005-2009 Shaun Purcell                  //
//                                                              //
// This file is distributed under the GNU General Public        //
// License, Version 2.  Please see the file COPYING for more    //
// details                                                      //
//                                                              //
//////////////////////////////////////////////////////////////////

#ifdef LOGISTIC_SSE

#include <iostream>
#include <iomanip>
#include <cmath>

#include "logistic.h"
#include "plink.h"
#include "helper.h"
#include "options.h"
#include "stats.h"


#include "pmmintrin.h"


/**
  exp_ps : exponential using single precision SSE from fast math library for float.
	@brief fast math library for float
	@author herumi
	@url http://homepage1.nifty.com/herumi/
	@note modified new BSD license
	http://opensource.org/licenses/BSD-3-Clause
*/
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



#ifdef PROFILING 
__attribute__((noinline))
#endif
inline void logistic_sse(float *vect, int n) {
  // Equivalent C++ code
  //  for(int i = 0; i < n; i++) vect[i] = 1. / (1 + exp(-vect[i]));
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

#ifdef PROFILING 
__attribute__((noinline))
#endif
inline void mult_tmatrix_nxd_vect_d(float** tm, const float* vect, float* dest, int n, int d) {
// equivalent C++ code
/*for(int i = 0; i < n; i++) {
    dest[i] = 0.;
    for(int j = 0; j < d; j++) 
      dest[i] += tm[j][i] * vect[j];
  }*/

  __m128 w1, w2, w3, w4;
  int k = 0;

  if(d < 4) memset(dest, 0, n*sizeof(float));
  else {
    w1 = _mm_load1_ps(vect);
    w2 = _mm_load1_ps(vect+1);
    w3 = _mm_load1_ps(vect+2);
    w4 = _mm_load1_ps(vect+3);
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

    for(k = 4; k + 3 < d; k += 4) {
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




#ifdef PROFILING 
__attribute__((noinline))
#endif
inline void mult_matrix_dxn_vect_n(float** m, const float* vect, float* dest, int n, int d) {
  // equivalent C++ code
  /* for(int j = 0; j < d; j++) dest[j] = 0.;
     for(int i = 0; i < n; i++) 
       for(int j = 0; j < d; j++) 
         dest[j] += m[j][i] * vect[i]; */

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

#ifdef PROFILING 
__attribute__((noinline))
#endif
inline void compute_hessian(float* const * m, const float* v, float** dest, int n,  int d) {
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


// solve linear system using cholesky decomposition
#ifdef PROFILING 
__attribute__((noinline))
#endif
void  solve_linear_system(const float* const* L, const float* Y, float* X, int d) {
  for(int i = 0; i < d; i++) {
    float s = Y[i];
    for(int j = 0; j < i; j++) s -= L[i][j]*X[j];
    X[i] = s / L[i][i];
  }
  for(int i = d-1; i >= 0; i--) {
    float s = X[i];
    for(int j = d-1; j > i; j--) s -= L[j][i]*X[j];
    X[i] = s / L[i][i];
  }
}

// we can get a H^-1 column by solving HX = (0,0,...,0,1) 
// as we only need the the last coef of X, we can do only half of the resolution using only one triangular matrix instead of two.
#ifdef PROFILING 
__attribute__((noinline))
#endif
float compute_wald(const float* const* L, int d) {
  float Y[d]; 
  float X[d];
  memset(Y, 0, d*sizeof(float));
  Y[d-1] = 1.;

  // only first triangular resolution of solve_linear_system
  for(int i = 0; i < d; i++) {
    float s = Y[i];
    for(int j = 0; j < i; j++) s -= L[i][j]*X[j];
    X[i] = s / L[i][i];
  }
  return X[d-1] / L[d-1][d-1];
}



// see http://en.wikipedia.org/wiki/Cholesky_decomposition
#ifdef PROFILING 
__attribute__((noinline))
#endif
void  cholesky_decompsition(const float* const * A, float** L, int d) {
  for(int j = 0; j < d; j++) {
    float s = A[j][j];
    for(int k = 0; k < j; k++) s -= L[j][k]*L[j][k];
    if(s < 0) L[j][j] = 1e-6;
    else L[j][j] = sqrt(s);

    for(int i = j+1; i < d; i++) {
      s = A[i][j];
      for(int k = 0; k < j; k++) s -= L[i][k]*L[j][k];
      L[i][j] = s / L[j][j];
    }
  }
}


class FastLM {
  private:
    float** X;
    float* Y;
    float* grad; 

    float* _H;
    float** H;
    float* _L;

    int np;
    int nind;
    float min_delta_coef;

    FastLM(const FastLM&) {
      // make copy constructor private to forbid copy at compilation time.
      // as the class contains a lot of data, it would be inefficient to copy it.
    }   
  public:
    float* P;
    float* V;

    float* coef;
    int iteration;
    float delta_coef;
    float** L;

    FastLM(const matrix_t &_X, float* _Y, float* _P, float* _V) {
      Y = _Y;
      iteration = 0;
      nind = _X.size();
      np = _X[0].size();
      min_delta_coef = 1e9;

      int sz = 4*((nind+3)/4)*sizeof(float);
      X = new float*[np];

      for(int i = 0; i < np; i++) {
        X[i] = (float*) _mm_malloc(sz, 16);
        for(int j = 0; j < nind; j++) X[i][j] = _X[j][i];
        for(int j = nind; j < 4*((nind+3)/4); j++) X[i][j] = 0.;
      }


      V = _V;
      P = _P;
      memset(P, 0, sz);

      sz = 4*((np+3)/4)*sizeof(float);  
      grad = (float*) _mm_malloc(sz, 16);
      coef = (float*) _mm_malloc(sz, 16);
      _H = (float*) _mm_malloc(sz*np, 16);
      _L = (float*) _mm_malloc(sz*np, 16);

      memset(_H, 0, sz*np);
      memset(_L, 0, sz*np);
      memset(grad, 0, sz);
      memset(coef, 0, sz);

      H = new float*[np];
      L = new float*[np];
      for(int i = 0; i < np; i++) {
        H[i] = _H + i*4*((np+3)/4);
        L[i] = _L + i*4*((np+3)/4);
      }
    }

    ~FastLM() {
      for(int i = 0; i < np; i++) _mm_free(X[i]);
      delete[] X;
      delete[] H;
      delete[] L;
      _mm_free(_H);
      _mm_free(_L);
      _mm_free(grad);
      _mm_free(coef);
    }



    bool make_iteration() {
      iteration++;

      // P[i] = SUM_j coef[j] * _X[i][j]
      mult_tmatrix_nxd_vect_d(X, coef, P, nind, np);

      // P[i] = 1 / (1 + exp(-P[i]))
      logistic_sse(P, nind);

      //for(int i = 0; i < T.n; i++) {
      //   _V[i] = _P[i]*(1-_P[i]);
      //   _P[i] -= Y[i];
      // }
      compute_V_and_P_minus_Y(P, V, Y, nind);

      compute_hessian(X, V, H, nind, np);

      mult_matrix_dxn_vect_n(X, P, grad, nind, np);

      cholesky_decompsition(H, L, np);

      float dcoef[np];
      memset(dcoef, 0, sizeof(dcoef));
      solve_linear_system(L, grad, dcoef, np);

      delta_coef = 0.;
      for (int i = 0; i < np; i++) {
        delta_coef += abs(dcoef[i]);
        coef[i] -= dcoef[i];
      }
      if(delta_coef < min_delta_coef) min_delta_coef = delta_coef;

      if(delta_coef != delta_coef || (iteration > 4 && delta_coef > 20 && delta_coef > 2*min_delta_coef)) {
        memset(coef, 0, np*sizeof(float));
        return false;
      }

      return true;
    }
};


LogisticModel::LogisticModel(Plink * p_)
{
  P = p_;
  nc = 0;
  cluster = false;
  Y = 0;
  p = 0;
  V = 0;
}

LogisticModel::~LogisticModel()
{
  if(Y) _mm_free(Y);
  if(p) _mm_free(p);
  if(V) _mm_free(V);
}

void LogisticModel::setDependent()
{
  if(nind == 0 || nind != X.size()) // check that X and nind are coorectly initialized;
  error("Internal error: bad call to LogisticModel::setDependent()");

  int sz = 4*((nind+3)/4);
  if(Y) _mm_free(Y);
  Y = (float*) _mm_malloc(sz*sizeof(float), 16);

  for (int i = 0, j = 0; i < miss.size() && j < nind; i++) if(!miss[i]) Y[j++] = P->sample[i]->pperson->aff;
  for (int j = nind; j < sz; j++) Y[j] = 0.;

  if(p) _mm_free(p);
  p = (float*) _mm_malloc(sz*sizeof(float), 16);
  if(V) _mm_free(V);
  V = (float*) _mm_malloc(sz*sizeof(float), 16);
}


// INPUT : X, Y, nind, np
// OUTPUT : coef, S, (p, V)
void LogisticModel::fitLM() 
{
  coef.resize(np);
  sizeMatrix(S,np,np);

  if (np==0 || nind==0 || ! all_valid )
    return;

  if (par::verbose)
  {
    for (int i=0; i<nind; i++)
    {
      cout << i << "\t"
        << Y[i] << "\t";
      for (int j=0; j<np; j++)
        cout << X[i][j] << "\t";
      cout << "\n";
    }
  }



  FastLM lm(X, Y, p, V);

  while(true) {

    if(!lm.make_iteration()) {
      all_valid = false;
      return;
    }

    if(lm.iteration >= 8 && abs(1. - lm.delta_coef) < 1e-3) {
      all_valid = false;
      return;
    }

    if(lm.delta_coef < 1e-4
        || lm.iteration >= 15
        //          || (lm.delta_coef < 1e-3 && abs(T.W - old_W) < 0.1)  
        //          || (T.W < 1e-3 - 1e-6 && T.iteration > 2 && T.delta_coef < 10.) || 
        //              (T.W == 0 && T.restart)
      ) break;
  }


  // Compute S
  for(int i = 0; i < np; i++) {
    float y[np], x[np];
    memset(y, 0, sizeof(y));
    y[i] = 1.;
    solve_linear_system(lm.L, y, x, np);
    for(int j = 0; j < np; j++) S[i][j] = x[j];
  }

  for(int i = 0; i < np; i++) coef[i] = lm.coef[i];

  if ( cluster )
    HuberWhite();

  if (par::verbose)
  {
    cout << "beta\n";
    display(coef);
    cout << "Sigma\n";
    display(S);
    cout << "\n";
  }
}

vector_t LogisticModel::getCoefs()
{
  return coef;
}

vector_t LogisticModel::getVar()
{  
  vector_t var(np);
  for (int i=0; i<np; i++) 
    var[i] = S[i][i];
  return var;
}

vector_t LogisticModel::getSE()
{  
  vector_t var(np);
  for (int i=0; i<np; i++) 
    var[i] = sqrt(S[i][i]);
  return var;
}

void LogisticModel::reset()
{
  np=0;
  nind=0;
  coef.clear();
  S.clear();  
  if(Y) _mm_free(Y);
  Y = 0;
  if(p) _mm_free(p);
  p = 0;
  if(V) _mm_free(V);
  V = 0;
  X.clear();
  miss.clear();
}


void LogisticModel::displayResults(ofstream & OUT, Locus * loc)
{

  vector_t var;
  
  if ( all_valid ) 
    var = getVar();
  else
    {
      var.clear();
      var.resize(np,0);
    }
  

  for (int p=1; p<np; p++) // Skip intercept
    {
 
      bool okay = var[p] < 1e-20 || !realnum(var[p]) ? false : all_valid;
      
      double se = 0; 
      double Z = 0;
      double pvalue = 1;
      
      if (okay)
	{
	  se = sqrt(var[p]);
	  Z = coef[p] / se;
	  //	  pvalue = pT(Z,Y.size()-np);
	  pvalue = chiprobP(Z*Z,1);
	}
      
      // If filtering p-values
      if ( (!par::pfilter) || pvalue <= par::pfvalue ) 
	{	 

	  // Skip covariates?
	  if ( par::no_show_covar && p != testParameter )
	    continue;

	  OUT << setw(4) << loc->chr  << " " 
	      << setw(par::pp_maxsnp) << loc->name << " " 
	      << setw(10) << loc->bp << " " 
	      << setw(4) << loc->allele1 << " "
	      << setw(10) << label[p] << " "
	      << setw(8) << nind << " ";
	  
	  if (okay)
	    {
	      if ( par::return_beta )
		OUT << setw(10) << coef[p] << " ";
	      else
		OUT << setw(10) << exp(coef[p]) << " ";
	      
	      if (par::display_ci)
		{
		  OUT << setw(8) << se << " ";

		  if ( par::return_beta )
		    OUT << setw(8) << coef[p] - par::ci_zt * se << " "
			<< setw(8) << coef[p] + par::ci_zt * se << " ";	    
		  else
		    OUT << setw(8) << exp(coef[p] - par::ci_zt * se) << " "
			<< setw(8) << exp(coef[p] + par::ci_zt * se) << " ";	    

		}	      

	      OUT << setw(12) << Z << " "
		  << setw(12) << pvalue;
	    }
	  else
	    {
	      OUT << setw(10) << "NA" << " ";
	      
	      if (par::display_ci)
		OUT << setw(8) << "NA" << " "  
		    << setw(8) << "NA" << " "
		    << setw(8) << "NA" << " ";	    
	      
	      OUT << setw(12) << "NA" << " "
		  << setw(12) << "NA";
	    }
	  
	  OUT << "\n";	
	}
    }
}

double LogisticModel::getPValue()
{  
  vector_t var = getVar();
  bool okay = var[testParameter] < 1e-20 || !realnum(var[testParameter]) ? false : all_valid;

  if (okay)
  {
    double se = sqrt(var[testParameter]);
    double Z = coef[testParameter] / se;	  
    return chiprobP(Z*Z,1);
  }
  else return 1;
}

vector_t LogisticModel::getPVals()
{
  int tmp = testParameter;
  vector_t res;
  for ( testParameter = 1; testParameter < np; testParameter++)
    res.push_back( getPValue() );
  testParameter = tmp;
  return res;
}


double LogisticModel::getLnLk()
{
  // Return -2 * sample log-likelihood
  // We assume the model is fit, and all Y's are either 0 or 1

  double lnlk = 0;

  for (int i=0; i<nind; i++)
  {
    double t = 0;
    for (int j=0; j<np; j++)
      t += coef[j] * X[i][j];	    
    lnlk += Y[i] == 1 ? log( 1/(1+exp(-t))) : log(1 - (1/(1+exp(-t))) );
  }

  return -2 * lnlk;

}


void LogisticModel::HuberWhite()
{

  // Calculate sandwich variance estimators, potentially allowing for
  // clustered data

  // Works to update the S matrix, variance/covariance matrix

  // Originally, S will contain this, uncorrected

  // Calcuate S = (XtX)^-1


  matrix_t S0 = S; 

  vector<vector_t> sc(nc);
  for (int i=0; i<nc; i++)
    sc[i].resize(np,0);

  for (int i=0; i<nind; i++)
  {
    double err = Y[i] - p[i];      
    for (int j=0; j<np; j++)
      sc[clst[i]][j] += err * X[i][j];
  }

  matrix_t meat;
  sizeMatrix(meat, np, np);
  for (int k=0; k<nc; k++)
  {      
    for (int i=0; i<np; i++)
      for (int j=0; j<np; j++)
        meat[i][j] += sc[k][i] * sc[k][j];

  }

  matrix_t tmp1;
  multMatrix( S0 , meat, tmp1);
  multMatrix( tmp1 , S0, S);

}

bool LogisticModel::need_checkVIF() {
  return false;
}


#endif //#ifdef LOGISTIC_SSE

