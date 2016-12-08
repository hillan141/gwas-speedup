#include "commons.h"
#include "logistic.h"    

// utils (non-class members)

double start_time, end_time;

//report errors
void error(const char* s)
{
  cerr << s << endl;
  abort();
}

      void vector_index_error(const void* v, size_t index, size_t size)
      {
	// fTR(v|index|size);
          abort();
      }
      const size_t XKM = 0;
      const size_t XK1 = 1;
      vector<float> reference;

      inline size_t unaligned(const void*p) {
          return size_t(p)&15;
      }

      inline bool aligned(const void*p) {
          return !unaligned(p);
      }

      ostream& operator<<(ostream& out, const vector_t& m)
      {
          out << "{";
          forIter ( j, m )
              out << "  " << *j;
          out << " }";
          return out;
      }

      ostream& operator<<(ostream& out, const matrix_t& m)
      {
          out << "{\n";
          forIter ( i, m ) {
              forIter ( j, *i )
                  out << "  " << *j;
              out << "\n";
          }
          out << "}";
          return out;
      }
      //===============from stats.cpp====================

      #ifdef __SSE3__
      inline __m128d _mm_hadd_d(__m128d a)
      {
          return _mm_hadd_pd(a, a);
      }
      inline __m128 _mm_hadd_s(__m128 a)
      {
          a = _mm_hadd_ps(a, a);
          a = _mm_hadd_ps(a, a);
          return a;
      }
      #else
      inline __m128d _mm_hadd_d(__m128d a)
      {
          return _mm_add_sd(a, _mm_shuffle_pd(a, a, 3));
      }
      inline __m128 _mm_hadd_s(__m128 a)
      {
          a = _mm_add_ps(a, _mm_shuffle_ps(a, a, 0xee));
          a = _mm_add_ps(a, _mm_shuffle_ps(a, a, 0x55));
          return a;
      }
      #endif

      static __m128 row_mask;

      inline void init_row_mask(size_t N)
      {
          static const uint32_t mm[4][4] ALIGN16 = {
              { ~0, ~0, ~0, ~0 },
              { ~0, 0, 0, 0 },
              { ~0, ~0, 0, 0 },
              { ~0, ~0, ~0, 0 }
          };
          row_mask = _mm_load_ps((const float*)mm[N%4]);
      }

      inline void add_tail_zeros(vector_t& v)
      {
          size_t S = v.size();
          size_t S4 = (S+3)&-4;
          v.resize(S4);
          v.resize(S);
      }

      inline void mul_row(size_t N, double* ap, double b)
      {
          if ( unaligned(ap) ) {
              *ap *= b;
              ++ap;
              --N;
          }
          __m128d bb = _mm_set1_pd(b);
          const double* ae = ap+N;
          for ( const double* ae4 = ap+(N&-4); ap < ae4; ap += 4 ) {
              __m128d a0 = _mm_load_pd(ap  );
              __m128d a1 = _mm_load_pd(ap+2);
              a0 = _mm_mul_pd(a0, bb);
              a1 = _mm_mul_pd(a1, bb);
              _mm_store_pd(ap  , a0);
              _mm_store_pd(ap+2, a1);
          }
          if ( ap+2 <= ae ) {
              __m128d a0 = _mm_load_pd(ap  );
              a0 = _mm_mul_pd(a0, bb);
              _mm_store_pd(ap  , a0);
              ap += 2;
          }
          if ( ap < ae ) {
              *ap *= b;
          }
      }

      inline void xmul_row(size_t N, float* ap, float b)
      {
          assert(aligned(ap));
          __m128 bb = _mm_set1_ps(b);
          const float* ae = ap+N;
          for ( ; ap < ae-4; ap += 8 ) {
              __m128 a0 = _mm_load_ps(ap  );
              __m128 a1 = _mm_load_ps(ap+4);
              a0 = _mm_mul_ps(a0, bb);
              a1 = _mm_mul_ps(a1, bb);
              _mm_store_ps(ap  , a0);
              _mm_store_ps(ap+4, a1);
          }
          if ( ap < ae ) {
              __m128 a0 = _mm_load_ps(ap  );
              a0 = _mm_mul_ps(a0, bb);
              _mm_store_ps(ap  , a0);
          }
      }

      inline void xsub_mul_row(size_t N, const float* ap, float b, float* rp)
      {
          assert(aligned(ap));
          __m128 bb = _mm_set1_ps(b);
          const float* ae = ap+N;
          for ( ; ap < ae-4; ) {
              __m128 a0 = _mm_load_ps(ap  );
              __m128 a1 = _mm_load_ps(ap+4);
              ap += 8;
              __m128 r0 = _mm_load_ps(rp  );
              __m128 r1 = _mm_load_ps(rp+4);
              a0 = _mm_mul_ps(a0, bb);
              a1 = _mm_mul_ps(a1, bb);
              r0 = _mm_sub_ps(r0, a0);
              r1 = _mm_sub_ps(r1, a1);
              _mm_store_ps(rp  , r0);
              _mm_store_ps(rp+4, r1);
              rp += 8;
          }
          if ( ap < ae ) {
              __m128 a0 = _mm_load_ps(ap  );
              __m128 r0 = _mm_load_ps(rp  );
              a0 = _mm_mul_ps(a0, bb);
              r0 = _mm_sub_ps(r0, a0);
              _mm_store_ps(rp  , r0);
          }
      }

      inline void sub_mul_row(size_t N, const double* ap, double b, double* rp)
      {
          if ( unaligned(ap) ) {
              *rp -= *ap * b;
              ++ap;
              ++rp;
              --N;
          }
          __m128d bb = _mm_set1_pd(b);
          const double* ae = ap+N;
          for ( ; ap <= ae-4; ) {
              __m128d a0 = _mm_load_pd(ap  );
              __m128d a1 = _mm_load_pd(ap+2);
              ap += 4;
              __m128d r0 = _mm_load_pd(rp  );
              __m128d r1 = _mm_load_pd(rp+2);
              a0 = _mm_mul_pd(a0, bb);
              a1 = _mm_mul_pd(a1, bb);
              r0 = _mm_sub_pd(r0, a0);
              r1 = _mm_sub_pd(r1, a1);
              _mm_store_pd(rp  , r0);
              _mm_store_pd(rp+2, r1);
              rp += 4;
          }
          if ( ap <= ae-2 ) {
              __m128d a0 = _mm_load_pd(ap  );
              ap += 2;
              __m128d r0 = _mm_load_pd(rp  );
              a0 = _mm_mul_pd(a0, bb);
              r0 = _mm_sub_pd(r0, a0);
              _mm_store_pd(rp  , r0);
              rp += 2;
          }
          if ( ap < ae ) {
              *rp -= *ap * b;
          }
      }

      inline double& sum_mul(size_t N, const double* ap, const double* bp, 
      double& r)
      {
          assert(aligned(ap));
          const double* ae =ap+N;
          __m128d t0 = _mm_setzero_pd();
          __m128d t1 = _mm_setzero_pd();
          for ( ; ap <= ae-4; ) {
              __m128d a0 = _mm_load_pd(ap);
              __m128d b0 = _mm_load_pd(bp);
              __m128d a1 = _mm_load_pd(ap+2);
              __m128d b1 = _mm_load_pd(bp+2);
              ap += 4;
              bp += 4;
              t0 = _mm_add_pd(t0, _mm_mul_pd(a0, b0));
              t1 = _mm_add_pd(t1, _mm_mul_pd(a1, b1));
          }
          for ( ; ap <= ae-2; ) {
              __m128d a0 = _mm_load_pd(ap); ap += 2;
              __m128d b0 = _mm_load_pd(bp); bp += 2;
              t0 = _mm_add_pd(t0, _mm_mul_pd(a0, b0));
          }
          if ( ap < ae ) {
              __m128d a0 = _mm_load_sd(ap);
              __m128d b0 = _mm_load_sd(bp);
              t0 = _mm_add_pd(t0, _mm_mul_sd(a0, b0));
          }
          t0 = _mm_add_pd(t0, t1);
          t0 = _mm_hadd_d(t0);
          _mm_store_sd(&r, t0);
          return r;
      }

      inline void sum_mul2(size_t N,
                           const double* ap,
                           const double* bp,
                           const double* cp,
                           double* r)
      {
          const double* ae = ap+N;
          assert(aligned(ap));
          __m128d t0 = _mm_setzero_pd();
          __m128d s0 = _mm_setzero_pd();
          for ( ; ap <= ae-4; ) {
              __m128d a0 = _mm_load_pd(ap);
              __m128d b0 = _mm_load_pd(bp);
              __m128d c0 = _mm_load_pd(cp);
              __m128d a1 = _mm_load_pd(ap+2);
              __m128d b1 = _mm_load_pd(bp+2);
              __m128d c1 = _mm_load_pd(cp+2);
              ap += 4;
              bp += 4;
              cp += 4;
              b0 = _mm_mul_pd(b0, a0);
              c0 = _mm_mul_pd(c0, a0);
              t0 = _mm_add_pd(t0, b0);
              s0 = _mm_add_pd(s0, c0);
              b1 = _mm_mul_pd(b1, a1);
              c1 = _mm_mul_pd(c1, a1);
              t0 = _mm_add_pd(t0, b1);
              s0 = _mm_add_pd(s0, c1);
          }
          for ( ; ap <= ae-2; ) {
              __m128d a0 = _mm_load_pd(ap); ap += 2;
              __m128d b0 = _mm_load_pd(bp); bp += 2;
              __m128d c0 = _mm_load_pd(cp); cp += 2;
              b0 = _mm_mul_pd(b0, a0);
              c0 = _mm_mul_pd(c0, a0);
              t0 = _mm_add_pd(t0, b0);
              s0 = _mm_add_pd(s0, c0);
          }
          if ( ap < ae ) {
              __m128d a0 = _mm_load_sd(ap);
              __m128d b0 = _mm_load_sd(bp);
              __m128d c0 = _mm_load_sd(cp);
              b0 = _mm_mul_pd(b0, a0);
              c0 = _mm_mul_pd(c0, a0);
              t0 = _mm_add_pd(t0, b0);
              s0 = _mm_add_pd(s0, c0);
          }
          t0 = _mm_hadd_d(t0);
          s0 = _mm_hadd_d(s0);
          _mm_store_sd(r  , t0);
          _mm_store_sd(r+1, s0);
      }

      inline float& sum_mul(size_t N, const float* ap, const float* bp, float& 
r)
      {
          assert(aligned(ap));
          const float* ae =ap+N;
          __m128 t0 = _mm_setzero_ps();
          for ( ; ap <= ae-8; ) {
              __m128 a0 = _mm_load_ps(ap);
              __m128 b0 = _mm_load_ps(bp);
              __m128 a1 = _mm_load_ps(ap+4);
              __m128 b1 = _mm_load_ps(bp+4);
              ap += 8;
              bp += 8;
              t0 = _mm_add_ps(t0, _mm_mul_ps(a0, b0));
              t0 = _mm_add_ps(t0, _mm_mul_ps(a1, b1));
          }
          if ( ap <= ae-4 ) {
              __m128 a0 = _mm_load_ps(ap);
              __m128 b0 = _mm_load_ps(bp);
              ap += 4;
              bp += 4;
              t0 = _mm_add_ps(t0, _mm_mul_ps(a0, b0));
          }
          while ( ap < ae ) {
              __m128 a0 = _mm_load_ss(ap);
              ++ap;
              __m128 b0 = _mm_load_ss(bp);
              ++bp;
              t0 = _mm_add_ps(t0, _mm_mul_ss(a0, b0));
          }
          t0 = _mm_hadd_s(t0);
          _mm_store_ss(&r, t0);
          return r;
      }

      inline void sum_mul12(size_t N,
                            const float* ap,
                            const float* bp,
                            const float* cp,
                            float* r)
      {
          const float* ae = ap+N;
          assert(aligned(ap));
          __m128 t0 = _mm_setzero_ps();
          __m128 s0 = _mm_setzero_ps();
          for ( ; ap <= ae-8; ) {
              __m128 a0 = _mm_load_ps(ap);
              __m128 b0 = _mm_load_ps(bp);
              __m128 c0 = _mm_load_ps(cp);
              __m128 a1 = _mm_load_ps(ap+4);
              __m128 b1 = _mm_load_ps(bp+4);
              __m128 c1 = _mm_load_ps(cp+4);
              ap += 8;
              bp += 8;
              cp += 8;
              b0 = _mm_mul_ps(b0, a0);
              c0 = _mm_mul_ps(c0, a0);
              t0 = _mm_add_ps(t0, b0);
              s0 = _mm_add_ps(s0, c0);
              b1 = _mm_mul_ps(b1, a1);
              c1 = _mm_mul_ps(c1, a1);
              t0 = _mm_add_ps(t0, b1);
              s0 = _mm_add_ps(s0, c1);
          }
          if ( ap <= ae-4 ) {
              __m128 a0 = _mm_load_ps(ap);
              __m128 b0 = _mm_load_ps(bp);
              __m128 c0 = _mm_load_ps(cp);
              ap += 4;
              bp += 4;
              cp += 4;
              b0 = _mm_mul_ps(b0, a0);
              c0 = _mm_mul_ps(c0, a0);
              t0 = _mm_add_ps(t0, b0);
              s0 = _mm_add_ps(s0, c0);
          }
          while ( ap < ae ) {
              __m128 a0 = _mm_load_ss(ap);
              ++ap;
              __m128 b0 = _mm_load_ss(bp);
              ++bp;
              __m128 c0 = _mm_load_ss(cp);
              ++cp;
              b0 = _mm_mul_ss(b0, a0);
              c0 = _mm_mul_ss(c0, a0);
              t0 = _mm_add_ss(t0, b0);
              s0 = _mm_add_ss(s0, c0);
          }
          t0 = _mm_hadd_s(t0);
          s0 = _mm_hadd_s(s0);
          _mm_store_ss(r  , t0);
          _mm_store_ss(r+1, s0);
      }

      inline __m128 _mm_slli_ps(__m128 v, int c)
      {
          return (__m128)_mm_slli_si128((__m128i)v, c*4);
      }

      inline float& xsum_mul(size_t N, const float* ap, const float* bp, float& 
      r)
      {
          assert(aligned(ap));
          const float* ae =ap+N;
          __m128 t0 = _mm_setzero_ps();
          __m128 t1 = _mm_setzero_ps();
          for ( ; ap < ae-4; ) {
              __m128 a0 = _mm_load_ps(ap);
              __m128 b0 = _mm_load_ps(bp);
              __m128 a1 = _mm_load_ps(ap+4);
              __m128 b1 = _mm_load_ps(bp+4);
              ap += 8;
              bp += 8;
              t0 = _mm_add_ps(t0, _mm_mul_ps(a0, b0));
              t1 = _mm_add_ps(t1, _mm_mul_ps(a1, b1));
          }
          if ( ap < ae ) {
              __m128 a0 = _mm_load_ps(ap);
              __m128 b0 = _mm_load_ps(bp);
              t0 = _mm_add_ps(t0, _mm_mul_ps(a0, b0));
          }
          t0 = _mm_add_ps(t0, t1);
          t0 = _mm_hadd_s(t0);
          _mm_store_ss(&r, t0);
          return r;
      }

      inline void xsum_mul12(size_t N,
                             const float* ap,
                             const float* bp,
                             const float* cp,
                             float* r)
      {
          const float* ae = ap+N;
          assert(aligned(ap));
          __m128 t0 = _mm_setzero_ps();
          __m128 s0 = _mm_setzero_ps();
          for ( ; ap < ae-4; ) {
              __m128 a0 = _mm_load_ps(ap);
              __m128 b0 = _mm_load_ps(bp);
              __m128 c0 = _mm_load_ps(cp);
              __m128 a1 = _mm_load_ps(ap+4);
              __m128 b1 = _mm_load_ps(bp+4);
              __m128 c1 = _mm_load_ps(cp+4);
              ap += 8;
              bp += 8;
              cp += 8;
              b0 = _mm_mul_ps(b0, a0);
              c0 = _mm_mul_ps(c0, a0);
              t0 = _mm_add_ps(t0, b0);
              s0 = _mm_add_ps(s0, c0);
              b1 = _mm_mul_ps(b1, a1);
              c1 = _mm_mul_ps(c1, a1);
              t0 = _mm_add_ps(t0, b1);
              s0 = _mm_add_ps(s0, c1);
          }
          if ( ap < ae ) {
              __m128 a0 = _mm_load_ps(ap);
              __m128 b0 = _mm_load_ps(bp);
              __m128 c0 = _mm_load_ps(cp);
              b0 = _mm_mul_ps(b0, a0);
              c0 = _mm_mul_ps(c0, a0);
              t0 = _mm_add_ps(t0, b0);
              s0 = _mm_add_ps(s0, c0);
          }
          t0 = _mm_hadd_s(t0);
          s0 = _mm_hadd_s(s0);
          _mm_store_ss(r  , t0);
          _mm_store_ss(r+1, s0);
      }

      inline void xsum_mul21(size_t N,
                             const float* a1p, const float* a2p, const float* 
      b1p,
                             float& r1, float& r2)
      {
          assert(aligned(a1p));
          __m128 s11 = _mm_setzero_ps();
          __m128 s21 = _mm_setzero_ps();
          size_t a2p_off = (const char*)a2p-(const char*)a1p;
          size_t b1p_off = (const char*)b1p-(const char*)a1p;
      #define a2p (const float*)((const char*)a1p+a2p_off)
      #define b1p (const float*)((const char*)a1p+b1p_off)
          const float* a1e = a1p+N;
          do {
              __m128 b1 = _mm_load_ps(b1p);
              __m128 a1 = _mm_load_ps(a1p);
              __m128 a2 = _mm_load_ps(a2p);
              a1p += 4;
              __m128 p11 = _mm_mul_ps(a1, b1);
              __m128 p21 = _mm_mul_ps(a2, b1);
              s11 = _mm_add_ps(s11, p11);
              s21 = _mm_add_ps(s21, p21);
          } while ( a1p < a1e );
      #undef a2p
      #undef b1p
          __m128 s = _mm_hadd_ps(s11, s21);
          s = _mm_hadd_ps(s, s);
          _mm_store_ss(&r1, s);
          s = _mm_shuffle_ps(s, s, 0x55);
          _mm_store_ss(&r2, s);
      }

      inline void xsum_mul22(size_t N,
                             const float* a1p, const float* a2p,
                             const float* b1p, const float* b2p,
                             float* r1, float* r2)
      {
          assert(aligned(a1p));
          __m128 s11 = _mm_setzero_ps();
          __m128 s12 = _mm_setzero_ps();
          __m128 s21 = _mm_setzero_ps();
          __m128 s22 = _mm_setzero_ps();
          size_t a2p_off = (const char*)a2p-(const char*)a1p;
          size_t b1p_off = (const char*)b1p-(const char*)a1p;
          size_t b2p_off = (const char*)b2p-(const char*)a1p;
      #define a2p (const float*)((const char*)a1p+a2p_off)
      #define b1p (const float*)((const char*)a1p+b1p_off)
      #define b2p (const float*)((const char*)a1p+b2p_off)
          const float* a1e = a1p+N;
          do {
              __m128 a1 = _mm_load_ps(a1p);
              __m128 b1 = _mm_load_ps(b1p);
              __m128 p11 = _mm_mul_ps(a1, b1);
              s11 = _mm_add_ps(s11, p11);
              __m128 b2 = _mm_load_ps(b2p);
              __m128 a2 = _mm_load_ps(a2p); a1p += 4;
              __m128 p12 = _mm_mul_ps(a1, b2);
              s12 = _mm_add_ps(s12, p12);
              __m128 p21 = _mm_mul_ps(b1, a2);
              __m128 p22 = _mm_mul_ps(b2, a2);
              s21 = _mm_add_ps(s21, p21);
              s22 = _mm_add_ps(s22, p22);
          } while ( a1p < a1e );
      #undef a2p
      #undef b1p
      #undef b2p
          __m128 s1 = _mm_hadd_ps(s11, s12);
          __m128 s2 = _mm_hadd_ps(s21, s22);
          __m128 s = _mm_hadd_ps(s1, s2);
          _mm_store_sd((double*)r1, (__m128d)s);
          s = _mm_shuffle_ps(s, s, 0xee);
          _mm_store_sd((double*)r2, (__m128d)s);
      }

      inline void mul_row(size_t N, const double* ap, const double* bp, double* 
      rp)
      {
          const double* ae = ap+N;
          for ( ; ap <= ae-4; ) {
              __m128d a0 = _mm_load_pd(ap);
              __m128d a1 = _mm_load_pd(ap+2);
              ap += 4;
              __m128d b0 = _mm_load_pd(bp);
              __m128d b1 = _mm_load_pd(bp+2);
              bp += 4;
              a0 = _mm_mul_pd(a0, b0);
              a1 = _mm_mul_pd(a1, b1);
              _mm_store_pd(rp, a0);
              _mm_store_pd(rp+2, a1);
              rp += 4;
          }
          if ( ap <= ae-2 ) {
              __m128d a0 = _mm_load_pd(ap);
              ap += 2;
              __m128d b0 = _mm_load_pd(bp);
              bp += 2;
              a0 = _mm_mul_pd(a0, b0);
              _mm_store_pd(rp, a0);
              rp += 2;
          }
          if ( ap < ae ) {
              __m128d a0 = _mm_load_sd(ap);
              __m128d b0 = _mm_load_sd(bp);
              a0 = _mm_mul_sd(a0, b0);
              _mm_store_sd(rp, a0);
          }
      }

      inline void mul_row(size_t N, const float* ap, const float* bp, float* rp)
      {
          const float* ae = ap+N;
          for ( ; ap <= ae-8; ) {
              __m128 a0 = _mm_load_ps(ap);
              __m128 a1 = _mm_load_ps(ap+4);
              __m128 b0 = _mm_load_ps(bp);
              __m128 b1 = _mm_load_ps(bp+4);
              ap += 8;
              bp += 8;
              a0 = _mm_mul_ps(a0, b0);
              a1 = _mm_mul_ps(a1, b1);
              _mm_store_ps(rp  , a0);
              _mm_store_ps(rp+4, a1);
              rp += 8;
          }
          if ( ap <= ae-4 ) {
              __m128 a0 = _mm_load_ps(ap);
              ap += 4;
              __m128 b0 = _mm_load_ps(bp);
              bp += 4;
              a0 = _mm_mul_ps(a0, b0);
              _mm_store_ps(rp, a0);
              rp += 4;
          }
          while ( ap < ae ) {
              __m128 a0 = _mm_load_ss(ap);
              ap += 1;
              __m128 b0 = _mm_load_ss(bp);
              bp += 1;
              a0 = _mm_mul_ss(a0, b0);
              _mm_store_ss(rp, a0);
              rp += 1;
          }
      }

      inline void xmul_row(size_t N, const float* ap, const float* bp, float* 
rp)
      {
          assert(aligned(ap));
          const float* ae = ap+N;
          for ( ; ap < ae-4; ) {
              __m128 a0 = _mm_load_ps(ap);
              __m128 a1 = _mm_load_ps(ap+4);
              __m128 b0 = _mm_load_ps(bp);
              __m128 b1 = _mm_load_ps(bp+4);
              ap += 8;
              bp += 8;
              a0 = _mm_mul_ps(a0, b0);
              a1 = _mm_mul_ps(a1, b1);
              _mm_store_ps(rp  , a0);
              _mm_store_ps(rp+4, a1);
              rp += 8;
          }
          if ( ap < ae ) {
              __m128 a0 = _mm_load_ps(ap);
              __m128 b0 = _mm_load_ps(bp);
              a0 = _mm_mul_ps(a0, b0);
              _mm_store_ps(rp, a0);
          }
      }


     //resize the matrix to row x col, and clear to zero.
      inline void sizeMatrix(matrix_t &matrix,size_t row,size_t col)
      {
          matrix.resize(row);
          forIter ( i, matrix ) i->resize(col);
      }

 inline void clearVector(vector_t& v)
{
  void** pp = (void**)&v;
  pp[1] = pp[0];
}

inline void clearVector(vector<size_t>& v)
{
  void** pp = (void**)&v;
  pp[1] = pp[0];
}

inline bool realnum(double d)
{
  double zero = 0;
  if (d != d || d == 1/zero || d == -1/zero){
    return false;
  } else {
    return true;
  }
}


inline void xsizeVector(vector_t& v, size_t s)
{
  size_t ps = v.size();
  if ( s != ps ) {
    size_t s4 = (s+3)&-4;
    v.resize(s4);
    for ( size_t i = s; i < s4; ++i ) v[i] = 0;
    v.resize(s);
  }
}

inline void sub_row(size_t N, const double* ap, const double* bp, double* rp)
{
  const double* ae = ap+N;
  for ( ; ap <= ae-4; ) {
    __m128d a0 = _mm_load_pd(ap);
    __m128d a1 = _mm_load_pd(ap+2);
    ap += 4;
    __m128d b0 = _mm_load_pd(bp);
    __m128d b1 = _mm_load_pd(bp+2);
    bp += 4;
    a0 = _mm_sub_pd(a0, b0);
    a1 = _mm_sub_pd(a1, b1);
    _mm_store_pd(rp, a0);
    _mm_store_pd(rp+2, a1);
    rp += 4;
  }
  if ( ap <= ae-2 ) {
    __m128d a0 = _mm_load_pd(ap);
              ap += 2;
              __m128d b0 = _mm_load_pd(bp);
              bp += 2;
              a0 = _mm_sub_pd(a0, b0);
              _mm_store_pd(rp, a0);
              rp += 2;
  }
  if ( ap < ae ) {
    __m128d a0 = _mm_load_sd(ap);
    __m128d b0 = _mm_load_sd(bp);
    a0 = _mm_sub_sd(a0, b0);
    _mm_store_sd(rp, a0);
  }
}

inline void sub_row(size_t N, const float* ap, const float* bp, float* rp)
{
  const float* ae = ap+N;
  for ( ; ap <= ae-8; ) {
    __m128 a0 = _mm_load_ps(ap);
    __m128 b0 = _mm_load_ps(bp);
    __m128 a1 = _mm_load_ps(ap+4);
    __m128 b1 = _mm_load_ps(bp+4);
    ap += 8;
    bp += 8;
    a0 = _mm_sub_ps(a0, b0);
    a1 = _mm_sub_ps(a1, b1);
    _mm_store_ps(rp  , a0);
    _mm_store_ps(rp+4, a1);
    rp += 8;
  }
  if ( ap <= ae-4 ) {
    __m128 a0 = _mm_load_ps(ap);
    ap += 4;
    __m128 b0 = _mm_load_ps(bp);
    bp += 4;
    a0 = _mm_sub_ps(a0, b0);
    _mm_store_ps(rp, a0);
    rp += 4;
  }
  while ( ap < ae ) {
    __m128 a0 = _mm_load_ss(ap);
    ap += 1;
    __m128 b0 = _mm_load_ss(bp);
    bp += 1;
    a0 = _mm_sub_ss(a0, b0);
    _mm_store_ss(rp, a0);
    rp += 1;
  }
}

inline void xsub_row(size_t N, const float* ap, const float* bp, float* 
		     rp)
{
  assert(aligned(ap));
  const float* ae = ap+N;
  for ( ; ap < ae-4; ) {
    __m128 a0 = _mm_load_ps(ap);
    __m128 b0 = _mm_load_ps(bp);
    __m128 a1 = _mm_load_ps(ap+4);
    __m128 b1 = _mm_load_ps(bp+4);
    ap += 8;
    bp += 8;
    a0 = _mm_sub_ps(a0, b0);
    a1 = _mm_sub_ps(a1, b1);
    _mm_store_ps(rp  , a0);
    _mm_store_ps(rp+4, a1);
    rp += 8;
  }
  if ( ap < ae ) {
    __m128 a0 = _mm_load_ps(ap);
    __m128 b0 = _mm_load_ps(bp);
    a0 = _mm_sub_ps(a0, b0);
    _mm_store_ps(rp, a0);
  }
}

inline void sub_row(size_t N, const double* ap, double b, double* rp)
{
  const double* ae = ap+N;
  __m128d bb = _mm_set1_pd(b);
  for ( ; ap <= ae-4; ) {
    __m128d a0 = _mm_load_pd(ap);
    __m128d a1 = _mm_load_pd(ap+2);
    ap += 4;
    a0 = _mm_sub_pd(a0, bb);
    a1 = _mm_sub_pd(a1, bb);
    _mm_store_pd(rp  , a0);
    _mm_store_pd(rp+2, a1);
    rp += 4;
  }
  if ( ap <= ae-2 ) {
    __m128d a0 = _mm_load_pd(ap);
    ap += 2;
    a0 = _mm_sub_pd(a0, bb);
    _mm_store_pd(rp, a0);
    rp += 2;
  }
  if ( ap < ae ) {
    __m128d a0 = _mm_load_sd(ap);
    a0 = _mm_sub_sd(a0, bb);
    _mm_store_sd(rp, a0);
          }
}

inline void sub_row(size_t N, const float* ap, float b, float* rp)
{
  const float* ae = ap+N;
  __m128 bb = _mm_set1_ps(b);
  for ( ; ap <= ae-8; ) {
    __m128 a0 = _mm_load_ps(ap  );
    __m128 a1 = _mm_load_ps(ap+4);
    ap += 8;
              a0 = _mm_sub_ps(a0, bb);
              a1 = _mm_sub_ps(a1, bb);
              _mm_store_ps(rp  , a0);
              _mm_store_ps(rp+4, a1);
              rp += 8;
  }
  if ( ap <= ae-4 ) {
    __m128 a0 = _mm_load_ps(ap);
    ap += 4;
    a0 = _mm_sub_ps(a0, bb);
    _mm_store_ps(rp, a0);
    rp += 4;
  }
  while ( ap < ae ) {
    __m128 a0 = _mm_load_ss(ap);
    ap += 1;
    a0 = _mm_sub_ss(a0, bb);
    _mm_store_ss(rp, a0);
    rp += 1;
  }
}

inline void add_mul_row(size_t N, const double* ap, double b, double* rp)
{
  __m128d bb = _mm_set1_pd(b);
  const double* ae = ap+N;
  for ( ; ap <= ae-4; ) {
    __m128d a0 = _mm_load_pd(ap);
    __m128d a1 = _mm_load_pd(ap+2);
    ap += 4;
    __m128d r0 = _mm_load_pd(rp);
    __m128d r1 = _mm_load_pd(rp+2);
    a0 = _mm_mul_pd(a0, bb);
    a1 = _mm_mul_pd(a1, bb);
    r0 = _mm_add_pd(r0, a0);
    r1 = _mm_add_pd(r1, a1);
    _mm_store_pd(rp, r0);
    _mm_store_pd(rp+2, r1);
    rp += 4;
  }
  if ( ap+2 <= ae ) {
    __m128d a0 = _mm_load_pd(ap);
    ap += 2;
    __m128d r0 = _mm_load_pd(rp);
    a0 = _mm_mul_pd(a0, bb);
    r0 = _mm_add_pd(r0, a0);
    _mm_store_pd(rp, r0);
    rp += 2;
  }
  if ( ap < ae ) {
    __m128d a0 = _mm_load_sd(ap);
    __m128d r0 = _mm_load_sd(rp);
    a0 = _mm_mul_sd(a0, bb);
    r0 = _mm_add_sd(r0, a0);
    _mm_store_sd(rp, r0);
  }
}

inline void add_mul_row(size_t N, const float* ap, float b, float* rp)
{
  __m128 bb = _mm_set1_ps(b);
  const float* ae = ap+N;
  for ( ; ap <= ae-8; ) {
    __m128 a0 = _mm_load_ps(ap  );
    __m128 a1 = _mm_load_ps(ap+4);
    ap += 8;
    __m128 r0 = _mm_load_ps(rp  );
    __m128 r1 = _mm_load_ps(rp+4);
    a0 = _mm_mul_ps(a0, bb);
    a1 = _mm_mul_ps(a1, bb);
    r0 = _mm_add_ps(r0, a0);
    r1 = _mm_add_ps(r1, a1);
    _mm_store_ps(rp  , r0);
    _mm_store_ps(rp+4, r1);
    rp += 8;
  }
  if ( ap <= ae-4 ) {
    __m128 a0 = _mm_load_ps(ap);
    ap += 4;
    __m128 r0 = _mm_load_ps(rp);
    a0 = _mm_mul_ps(a0, bb);
    r0 = _mm_add_ps(r0, a0);
    _mm_store_ps(rp, r0);
    rp += 4;
  }
  while ( ap < ae ) {
    __m128 a0 = _mm_load_ss(ap);
    ap += 1;
    __m128 r0 = _mm_load_ss(rp);
    a0 = _mm_mul_ss(a0, bb);
    r0 = _mm_add_ss(r0, a0);
    _mm_store_ss(rp, r0);
    rp += 1;
  }
}


      inline void subVector(const vector_t& a, const vector_t& b, vector_t& r)
      {
          size_t N = a.size();
          r.resize(N);
          sub_row(N, &a[0], &b[0], &r[0]);
      }

      inline void subVector(const vector_t& a, Double b, vector_t& r)
      {
          size_t N = a.size();
          r.resize(N);
          sub_row(N, &a[0], b, &r[0]);
      }

      inline Double& dot(const vector_t& a, const vector_t& b, Double& r)
      {
          return sum_mul(a.size(), &a[0], &b[0], r);
      }

      inline Double& xdot(const vector_t& a, const vector_t& b, Double& r)
      {
          return xsum_mul(a.size(), &a[0], &b[0], r);
      }

      inline void xdot12(const vector_t& a, const vector_t& b, const vector_t& 
c,
                         Double* r)
      {
          xsum_mul12(a.size(), &a[0], &b[0], &c[0], r);
      }

      inline void xdot21(const vector_t& a1, const vector_t& a2,
                         const vector_t& b,
                         Double& r1, Double& r2)
      {
          return xsum_mul21(a1.size(), &a1[0], &a2[0], &b[0], r1, r2);
      }

      inline void xdot22(const vector_t& a1, const vector_t& a2,
                         const vector_t& b1, const vector_t& b2,
                         Double* r1, Double* r2)
      {
          xsum_mul22(a1.size(), &a1[0], &a2[0], &b1[0], &b2[0], r1, r2);
      }


inline void xdot2(const vector_t& a1, const vector_t& a2,
		  size_t i0, const matrix_t& b,
		  vector_t& r1, vector_t& r2)
{
  for ( size_t i = i0, np = b.size(); i < np; ++i ) {
    if ( i < np-1 ) {
      xdot22(a1, a2, b[i], b[i+1], &r1[i], &r2[i]);
      ++i;
    }
    else {
      xdot21(a1, a2, b[i], r1[i], r2[i]);
      break;
    }
  }
}

inline void dot12(const vector_t& a, const vector_t& b, const vector_t& c,
		  Double* r)
{
  sum_mul12(a.size(), &a[0], &b[0], &c[0], r);
}

inline void dot(const vector_t& a, size_t i0, const matrix_t& b, vector_t& r)
{
  for ( size_t i = i0, np = b.size(); i < np; ++i ) {
    if ( i < np-1 ) {
      dot12(a, b[i], b[i+1], &r[i]);
      ++i;
    }
    else {
      dot(a, b[i], r[i]);
      break;
    }
  }
}

inline void xdot(const vector_t& a, size_t i0, const matrix_t& b, 
		 vector_t& r)
{
  for ( size_t i = i0, np = b.size(); i < np; ++i ) {
    if ( i < np-1 ) {
      xdot12(a, b[i], b[i+1], &r[i]);
      ++i;
    }
    else {
      xdot(a, b[i], r[i]);
      break;
    }
  }
}

inline void dot(const vector_t& a, const matrix_t& b, vector_t& r)
{
  dot(a, 0, b, r);
}

inline void xdot(const vector_t& a, const matrix_t& b, vector_t& r)
{
  xdot(a, 0, b, r);
}

inline void calc_p_v(size_t N, double* pp, double* vp)
{
  if ( 1 ) {
    const double* pe = pp+N;
    for ( const double* pe2 = pp+(N&-2); pp < pe2; pp += 2, vp += 2 ) 
      {
	__m128d x = _mm_load_pd(pp);
	__m128d s = _mm_mul_pd(x, x);
	__m128d a = _mm_add_pd(s, _mm_set1_pd(3960.));
	__m128d b = _mm_add_pd(s, _mm_set1_pd(1232.));
	a = _mm_add_pd(_mm_mul_pd(a, s), _mm_set1_pd(2162160.));
	b = _mm_add_pd(_mm_mul_pd(b, s), _mm_set1_pd(336336.));
	a = _mm_add_pd(_mm_mul_pd(a, s), _mm_set1_pd(302702400.));
	b = _mm_add_pd(_mm_mul_pd(b, s), _mm_set1_pd(23063040.));
	a = _mm_add_pd(_mm_mul_pd(a, s), _mm_set1_pd(8821612800.));
	b = _mm_add_pd(_mm_mul_pd(b, s), _mm_set1_pd(196035840.));
	x = _mm_mul_pd(x, _mm_set1_pd(1/180.));
	__m128d d = _mm_div_pd(_mm_mul_pd(x, a), b);
	d = _mm_min_pd(d, _mm_set1_pd(.5));
	d = _mm_max_pd(d, _mm_set1_pd(-.5));
	__m128d p = _mm_sub_pd(_mm_set1_pd(.5), d);
	__m128d v = _mm_sub_pd(_mm_set1_pd(.25), _mm_mul_pd(d, d));
	_mm_store_pd(pp, p);
	_mm_store_pd(vp, v);
      }
    if ( pp < pe ) {
      double x = *pp;
      double s = x*x;
      double a = 
	(((s+3960.)*s+2162160.)*s+302702400.)*s+8821612800.;
      double b = (((s+1232.)*s+336336.)*s+23063040.)*s+196035840.;
      double d = (1/180.)*x*a/b;
      d = min(.5, max(-.5, d));
      *pp=.5-d;
      *vp=.25-d*d;
    }
  }
  else {
    forN ( i, N ) {
      double x = pp[i];
      double v = 1/(1+__builtin_exp(x));
      pp[i]=v;
      vp[i]=v*(1-v);
    }
  }
}

      inline void calc_p_v(size_t N, float* pp, float* vp)
      {
          if ( 1 ) {
              const float* pe = pp+N;
              for ( ; pp <= pe-4; pp += 4, vp += 4 ) {
                  __m128 x = _mm_load_ps(pp);
                  __m128 s = _mm_mul_ps(x, x);
                  __m128 a = _mm_add_ps(s, _mm_set1_ps(3960.));
                  __m128 b = _mm_add_ps(s, _mm_set1_ps(1232.));
                  a = _mm_add_ps(_mm_mul_ps(a, s), _mm_set1_ps(2162160.));
                  b = _mm_add_ps(_mm_mul_ps(b, s), _mm_set1_ps(336336.));
                  a = _mm_add_ps(_mm_mul_ps(a, s), _mm_set1_ps(302702400.));
                  b = _mm_add_ps(_mm_mul_ps(b, s), _mm_set1_ps(23063040.));
                  a = _mm_add_ps(_mm_mul_ps(a, s), _mm_set1_ps(8821612800.));
                  b = _mm_add_ps(_mm_mul_ps(b, s), _mm_set1_ps(196035840.));
                  x = _mm_mul_ps(x, _mm_set1_ps(1/180.));
                  __m128 d = _mm_div_ps(_mm_mul_ps(x, a), b);
                  d = _mm_min_ps(d, _mm_set1_ps(.5));
                  d = _mm_max_ps(d, _mm_set1_ps(-.5));
                  __m128 p = _mm_sub_ps(_mm_set1_ps(.5), d);
                  __m128 v = _mm_sub_ps(_mm_set1_ps(.25), _mm_mul_ps(d, d));
                  _mm_store_ps(pp, p);
                  _mm_store_ps(vp, v);
              }
              while ( pp < pe ) {
                  float x = *pp;
                  float s = x*x;
                  float a = (((s+3960.)*s+2162160.)*s+302702400.)*s+8821612800.;
                  float b = (((s+1232.)*s+336336.)*s+23063040.)*s+196035840.;
                  float d = (1/180.)*x*a/b;
                  d = min(.5f, max(-.5f, d));
                  *pp=.5-d;
                  *vp=.25-d*d;
                  ++pp;
                  ++vp;
              }
          }
          else {
              forN ( i, N ) {
                  double x = pp[i];
                  double v = 1/(1+__builtin_exp(x));
                  pp[i]=v;
                  vp[i]=v*(1-v);
              }
          }
      }

      inline void calc_p_v(vector_t& p, vector_t& V)
      {
          calc_p_v(p.size(), &p[0], &V[0]);
      }

inline void transpose(const matrix_t& a,
		      matrix_t& at)
{
  //TIMER(transpose);
  size_t R = a.size(), C = a[0].size();
  sizeMatrix(at,C,R);
  forN ( i, R ) forN ( j, C ) at[j][i] = a[i][j];
}

      inline void transposein(matrix_t& a)
      {
          //TIMER(transpose);
          size_t N = a.size();
          forN ( i, N ) forNF ( j, i+1, N ) SWAP(a[j][i], a[i][j]);
      }

      inline void multMatrixNT(const matrix_t& a,
                               const matrix_t& bt,
                               matrix_t& c)
      {
          TIMER(multMatrixNT);
          size_t cr = a.size(), cc = bt.size();
          assert(a[0].size() == bt[0].size());
          sizeMatrix(c,cr,cc);
          if ( 1 || cr >= cc ) {
              forN ( i, cr ) forN ( j, cc )
                  dot(a[i], bt[j], c[i][j]);
          }
          else {
              forN ( j, cc ) forN ( i, cr )
                  dot(bt[j], a[i], c[i][j]);
          }
      }

      void multMatrix(const matrix_t& a,
                      const matrix_t& b,
                      matrix_t& c)
      {
          TIMER(multMatrix);

          size_t ar = a.size();
          size_t br = b.size();
          if (ar == 0 || br == 0)
              error("Internal error: multiplying 0-sized matrices");
        
          size_t ac = a[0].size();
          size_t bc = b[0].size();
          if ( ac != br )
              error("Internal error: non-conformable matrices in multMatrix()"); 

        
          size_t cr = ar;
          size_t cc = bc;

          sizeMatrix(c,cr,cc);
        
          forN ( i, ar ) forN ( j, bc ) forN ( k, ac )
              c[i][j] += a[i][k] * b[k][j];
      }

      inline double pythag(double a, double b)
      {
          return __builtin_sqrt(a*a+b*b);
          a = abs(a);
          b = abs(b);
          if (a > b) return a*sqrt(1+SQR(b/a));
          else return (b == 0.0 ? 0.0 : b*sqrt(1+SQR(a/b)));
      }

      static vector_t rv1;

      // compute singular value decomposition
      bool svdcmp(matrix_t& a, 
                  vector_t& w, 
                  matrix_t& v)
      {
          TIMER(svdcmp);
          bool flag;
          int i,its,j,jj,k,l=0,nm;
          double anorm,c,f,g,h,s,scale,x,y,z;
          double volatile temp;

          int m=a.size();
          int n=a[0].size();

          rv1.resize(n);
          g=scale=anorm=0.0;
          {
              TIMER(svdcmp 1);
              for (i=0;i<n;i++) {
                  l=i+2;
                  rv1[i]=scale*g;
                  g=s=scale=0.0;
                  if (i < m) {
                      for (k=i;k<m;k++) scale += fabs(a[k][i]);
                      if (scale != 0.0) {
                          double inv_scale = 1/scale;
                          for (k=i;k<m;k++) {
                              a[k][i] *= inv_scale;
                              s += a[k][i]*a[k][i];
                          }
                          f=a[i][i];
                          g = -SIGN(sqrt(s),f);
                          h=f*g-s;
                          a[i][i]=f-g;
                          double inv_h = 1/h;
                          for (j=l-1;j<n;j++) {
                              for (s=0.0,k=i;k<m;k++) s += a[k][i]*a[k][j];
                              f=s*inv_h;
                              for (k=i;k<m;k++) a[k][j] += f*a[k][i];
                          }
                          for (k=i;k<m;k++) a[k][i] *= scale;
                      }
                  }
                  w[i]=scale *g;
                  g=s=scale=0.0;
                  if (i+1 <= m && i+1 != n) {
                      for (k=l-1;k<n;k++) scale += fabs(a[i][k]);
                      if (scale != 0.0) {
                          double inv_scale = 1/scale;
                          for (k=l-1;k<n;k++) {
                              a[i][k] *= inv_scale;
                              s += a[i][k]*a[i][k];
                          }
                          f=a[i][l-1];
                          g = -SIGN(sqrt(s),f);
                          h=f*g-s;
                          a[i][l-1]=f-g;
                          double inv_h = 1/h;
                          for (k=l-1;k<n;k++) rv1[k]=a[i][k]*inv_h;
                          for (j=l-1;j<m;j++) {
                              for (s=0.0,k=l-1;k<n;k++) s += a[j][k]*a[i][k];
                              for (k=l-1;k<n;k++) a[j][k] += s*rv1[k];
                          }
                          for (k=l-1;k<n;k++) a[i][k] *= scale;
                      }
                  }
                  anorm=MAX(anorm,double(fabs(w[i])+fabs(rv1[i])));
              }
          }
          {
              TIMER(svdcmp 2);
              for (i=n-1;i>=0;i--) {
                  if (i < n-1) {
                      if (g != 0.0) {
                          double mul = 1/(a[i][l]*g);
                          for (j=l;j<n;j++)
                              v[j][i]=a[i][j]*mul;
                          for (j=l;j<n;j++) {
                              for (s=0.0,k=l;k<n;k++) s += a[i][k]*v[k][j];
                              for (k=l;k<n;k++) v[k][j] += s*v[k][i];
                          }
                      }
                      for (j=l;j<n;j++) v[i][j]=v[j][i]=0.0;
                  }
                  v[i][i]=1.0;
                  g=rv1[i];
                  l=i;
              }
          }
          {
              TIMER(svdcmp 3);
              for (i=MIN(m,n)-1;i>=0;i--) {
                  l=i+1;
                  g=w[i];
                  for (j=l;j<n;j++) a[i][j]=0.0;
                  if (g != 0.0) {
                      g=1.0/g;
                      double mul = g/a[i][i];
                      for (j=l;j<n;j++) {
                          for (s=0.0,k=l;k<m;k++) s += a[k][i]*a[k][j];
                          f=s*mul;
                          for (k=i;k<m;k++) a[k][j] += f*a[k][i];
                      }
                      for (j=i;j<m;j++) a[j][i] *= g;
                  } else for (j=i;j<m;j++) a[j][i]=0.0;
                  ++a[i][i];
              }
          }
          {
              TIMER(svdcmp 4);
              for (k=n-1;k>=0;k--) {
                  for (its=0;its<30;its++) {
                      flag=true;
                      for (l=k;l>=0;l--) {
                          nm=l-1;
                          temp=fabs(rv1[l])+anorm;
                          if (temp == anorm) {
                              flag=false;
                              break;
                          }
                          temp=fabs(w[nm])+anorm;
                          if (temp == anorm) break;
                      }
                      if (flag) {
                          c=0.0;
                          s=1.0;
                          for (i=l;i<k+1;i++) {
                              f=s*rv1[i];
                              rv1[i]=c*rv1[i];
                              temp = fabs(f)+anorm;
                              if (temp == anorm) break;
                              g=w[i];
                              h=pythag(f,g);
                              w[i]=h;
                              h=1.0/h;
                              c=g*h;
                              s = -f*h;
                              for (j=0;j<m;j++) {
                                  y=a[j][nm];
                                  z=a[j][i];
                                  a[j][nm]=y*c+z*s;
                                  a[j][i]=z*c-y*s;
                              }
                          }
                      }
                      z=w[k];
                      if (l == k) {
                          if (z < 0.0) {
                              w[k] = -z;
                              for (j=0;j<n;j++) v[j][k] = -v[j][k];
                          }
                          break;
                      }
                      if (its == 29) 
                          return false; // cannot converge: multi-collinearity?
                      x=w[l];
                      nm=k-1;
                      y=w[nm];
                      g=rv1[nm];
                      h=rv1[k];
                      f=((y-z)*(y+z)+(g-h)*(g+h))/(2.0*h*y);
                      g=pythag(f,1.0);
                      f=((x-z)*(x+z)+h*((y/(f+SIGN(g,f)))-h))/x;
                      c=s=1.0;
                      for (j=l;j<=nm;j++) {
                          i=j+1;
                          g=rv1[i];
                          y=w[i];
                          h=s*g;
                          g=c*g;
                          z=pythag(f,h);
                          rv1[j]=z;
                          c=f/z;
                          s=h/z;
                          f=x*c+g*s;
                          g=g*c-x*s;
                          h=y*s;
                          y *= c;
                          for (jj=0;jj<n;jj++) {
                              x=v[jj][j];
                              z=v[jj][i];
                              v[jj][j]=x*c+z*s;
                              v[jj][i]=z*c-x*s;
                          }
                          z=pythag(f,h);
                          w[j]=z;
                          if (z) {
                              z=1.0/z;
                              c=f*z;
                              s=h*z;
                          }
                          f=c*g+s*y;
                          x=c*y-s*g;
                          for (jj=0;jj<m;jj++) {
                              y=a[jj][j];
                              z=a[jj][i];
                              a[jj][j]=y*c+z*s;
                              a[jj][i]=z*c-y*s;
                          }
                      }
                      rv1[l]=0.0;
                      rv1[k]=f;
                      w[k]=x;
                  }
              }
          }
          return true;
      }

      inline void rotate(size_t n, double* a, double* b, double c, double s)
      {
          if ( 1 ) {
              const double* ae = a+n;
              __m128d cc = _mm_set1_pd(c);
              __m128d ss = _mm_set1_pd(s);
              for ( const double* ae2 = a+(n&-2); a != ae2; a += 2, b += 2 ) {
                  __m128d xx = _mm_load_pd(a);
                  __m128d zz = _mm_load_pd(b);
                  __m128d xc = _mm_mul_pd(xx, cc);
                  __m128d xs = _mm_mul_pd(xx, ss);
                  __m128d zc = _mm_mul_pd(zz, cc);
                  __m128d zs = _mm_mul_pd(zz, ss);
                  __m128d aa = _mm_add_pd(xc, zs);
                  __m128d bb = _mm_sub_pd(zc, xs);
                  _mm_store_pd(a, aa);
                  _mm_store_pd(b, bb);
              }
              if ( a < ae ) {
                  double x=*a;
                  double z=*b;
                  *a=z*s+x*c;
                  *b=z*c-x*s;
              }
          }
          else {
              forN ( i, n ) {
                  double x=a[i];
                  double z=b[i];
                  a[i]=x*c+z*s;
                  b[i]=z*c-x*s;
              }
          }
      }

      inline void rotate(size_t n, float* a, float* b, float c, float s)
      {
          forN ( i, n ) {
              double x=a[i];
              double z=b[i];
              a[i]=x*c+z*s;
              b[i]=z*c-x*s;
          }
      }

      // compute singular value decomposition
      bool svdcmpt(matrix_t& a, 
                   vector_t& w, 
                   matrix_t& v)
      {
          TIMER(svdcmp);
          bool flag;
          int i,its,j,k,l=0,nm;
          double anorm,c,f,g,h,s,scale,x,y,z;
          double volatile temp;

          int m=a.size();
          int n=a[0].size();

          rv1.resize(n);
          g=scale=anorm=0.0;
          {
              TIMER(svdcmpt 1);
              for (i=0;i<n;i++) {
                  l=i+2;
                  rv1[i]=scale*g;
                  g=s=scale=0.0;
                  if (i < m) {
                      for (k=i;k<m;k++) scale += fabs(a[i][k]);
                      if (scale != 0.0) {
                          double inv_scale = 1/scale;
                          for (k=i;k<m;k++) {
                              s += SQR(a[i][k] *= inv_scale);
                          }
                          f=a[i][i];
                          g = -SIGN(sqrt(s),f);
                          h=f*g-s;
                          a[i][i]=f-g;
                          double inv_h = 1/h;
                          for (j=l-1;j<n;j++) {
                              for (s=0.0,k=i;k<m;k++) s += a[i][k]*a[j][k];
                              f=s*inv_h;
                              for (k=i;k<m;k++) a[j][k] += f*a[i][k];
                          }
                          for (k=i;k<m;k++) a[i][k] *= scale;
                      }
                  }
                  w[i]=scale *g;
                  g=s=scale=0.0;
                  if (i+1 <= m && i+1 != n) {
                      for (k=l-1;k<n;k++) scale += fabs(a[k][i]);
                      if (scale != 0.0) {
                          double inv_scale = 1/scale;
                          for (k=l-1;k<n;k++) {
                              s += SQR(a[k][i] *= inv_scale);
                          }
                          f=a[l-1][i];
                          g = -SIGN(sqrt(s),f);
                          h=f*g-s;
                          a[l-1][i]=f-g;
                          double inv_h = 1/h;
                          for (k=l-1;k<n;k++) rv1[k]=a[k][i]*inv_h;
                          for (j=l-1;j<m;j++) {
                              for (s=0.0,k=l-1;k<n;k++) s += a[k][j]*a[k][i];
                              for (k=l-1;k<n;k++) a[k][j] += s*rv1[k];
                          }
                          for (k=l-1;k<n;k++) a[k][i] *= scale;
                      }
                  }
                  anorm=MAX(anorm,double(fabs(w[i])+fabs(rv1[i])));
              }
          }
          {
              TIMER(svdcmpt 2);
              for (i=n-1;i>=0;i--) {
                  if (i < n-1) {
                      if (g != 0.0) {
                          double mul = 1/(a[l][i]*g);
                          for (j=l;j<n;j++)
                              v[i][j]=a[j][i]*mul;
                          for (j=l;j<n;j++) {
                              for (s=0.0,k=l;k<n;k++) s += a[k][i]*v[j][k];
                              for (k=l;k<n;k++) v[j][k] += s*v[i][k];
                          }
                      }
                      for (j=l;j<n;j++) v[j][i]=v[i][j]=0.0;
                  }
                  v[i][i]=1.0;
                  g=rv1[i];
                  l=i;
              }
          }
          {
              TIMER(svdcmpt 3);
              for (i=MIN(m,n)-1;i>=0;i--) {
                  l=i+1;
                  g=w[i];
                  for (j=l;j<n;j++) a[j][i]=0.0;
                  if (g != 0.0) {
                      g=1.0/g;
                      double mul = g/a[i][i];
                      for (j=l;j<n;j++) {
                          for (s=0.0,k=l;k<m;k++) s += a[i][k]*a[j][k];
                          f=s*mul;
                          for (k=i;k<m;k++) a[j][k] += f*a[i][k];
                      }
                      for (j=i;j<m;j++) a[i][j] *= g;
                  } else for (j=i;j<m;j++) a[i][j]=0.0;
                  ++a[i][i];
              }
          }
          {
              TIMER(svdcmpt 4);
              for (k=n-1;k>=0;k--) {
                  for (its=0;its<30;its++) {
                      flag=true;
                      for (l=k;l>=0;l--) {
                          nm=l-1;
                          temp=fabs(rv1[l])+anorm;
                          if (temp == anorm) {
                              flag=false;
                              break;
                          }
                          temp=fabs(w[nm])+anorm;
                          if (temp == anorm) break;
                      }
                      if (flag) {
                          c=0.0;
                          s=1.0;
                          for (i=l;i<k+1;i++) {
                              f=s*rv1[i];
                              rv1[i]=c*rv1[i];
                              temp = fabs(f)+anorm;
                              if (temp == anorm) break;
                              g=w[i];
                              h=pythag(f,g);
                              w[i]=h;
                              h=1.0/h;
                              c=g*h;
                              s = -f*h;
                              rotate(m, &a[nm][0], &a[i][0], c, s);
                          }
                      }
                      z=w[k];
                      if (l == k) {
                          if (z < 0.0) {
                              w[k] = -z;
                              for (j=0;j<n;j++) v[k][j] = -v[k][j];
                          }
                          break;
                      }
                      if (its == 29) 
                          return false; // cannot converge: multi-collinearity?
                      x=w[l];
                      nm=k-1;
                      y=w[nm];
                      g=rv1[nm];
                      h=rv1[k];
                      f=((y-z)*(y+z)+(g-h)*(g+h))/(2.0*h*y);
                      g=pythag(f,1.0);
                      f=((x-z)*(x+z)+h*((y/(f+SIGN(g,f)))-h))/x;
                      c=s=1.0;
                      for (j=l;j<=nm;j++) {
                          i=j+1;
                          g=rv1[i];
                          y=w[i];
                          h=s*g;
                          g=c*g;
                          z=pythag(f,h);
                          rv1[j]=z;
                          double inv_z = 1/z;
                          c=f*inv_z;
                          s=h*inv_z;
                          f=x*c+g*s;
                          g=g*c-x*s;
                          h=y*s;
                          y *= c;
                          rotate(n, &v[j][0], &v[i][0], c, s);
                          z=pythag(f,h);
                          w[j]=z;
                          if (z) {
                              z=1.0/z;
                              c=f*z;
                              s=h*z;
                          }
                          f=c*g+s*y;
                          x=c*y-s*g;
                          rotate(n, &a[j][0], &a[i][0], c, s);
                      }
                      rv1[l]=0.0;
                      rv1[k]=f;
                      w[k]=x;
                  }
              }
          }
          return true;
      }

inline void multVector(const vector_t& a, const vector_t& b, vector_t& r)
{
  size_t N = a.size();
  r.resize(N);
  mul_row(N, &a[0], &b[0], &r[0]);
}

      static vector_t w, t;

      // invert matrix using SVD
      bool svd_inverse(matrix_t& u, matrix_t& v)
      {
          TIMER(svd_inverse);
        
          const double eps = 1e-24; 
        
          size_t n = u.size();
        
          sizeMatrix(v, n, n);
          w.resize(n);
          t.resize(n);
          if ( !svdcmp(u,w,v) ) return 0;
        
          // Look for singular values
          double wmax = *max_element(ALL(w));
          double wmin = wmax * eps;
          forN ( i, n ) w[i] = w[i] < wmin ? 0 : 1/w[i];
          
          // results matrix
          forN ( i, n ) {
              multVector(u[i], w, t);
              dot(t, v, u[i]);
          }
          return 1;
      }

      // invert matrix using SVD
      bool svd_inverse_sym(matrix_t& u, matrix_t& v)
      {
          TIMER(svd_inverse_sym);
        
          const double eps = 1e-24; 
        
          size_t n = u.size();
        
          sizeMatrix(v, n, n);
          w.resize(n);
          t.resize(n);
          if ( !svdcmpt(u,w,v) ) return 0;
          transposein(u);
          transposein(v);
        
          // Look for singular values
          double wmax = *max_element(ALL(w));
          double wmin = wmax * eps;
          forN ( i, n ) w[i] = w[i] < wmin ? 0 : 1/w[i];
          
          // results matrix
          forN ( i, n ) {
              multVector(u[i], w, t);
              dot(t, i, v, u[i]);
          }
          forN ( i, n ) {
              forNF ( j, i+1, n ) {
                  u[j][i] = u[i][j];
              }
          }
          return 1;
      }

      inline bool invertGauss(matrix_t& m)
      {
          const size_t MAX_N = 16;
          static Double t[MAX_N][MAX_N*2];
          TIMER(invertGauss);
          size_t n = m.size(), S = (2*n+3)&-4;
          assert(n <= MAX_N);
          assert(S <= MAX_N*2);
          const Double EPS = 1e-12;
          {
              //TIMER(invertGauss input);
              forN ( i, n ) {
                  const vector_t& s = m[i];
                  Double* d = t[i];
                  copy(ALL(s), d);
                  fill(d+n, d+S, 0.);
                  d[n+i] = 1;
              }
          }
          static bool used[MAX_N];
          fill_n(used, n, false);
          static size_t order[MAX_N];
          forN ( k, n ) {
              size_t s = 0;
              Double max_v = 0;
              forN ( i, n ) {
                  if ( used[i] ) continue;
                  Double v = abs(t[i][k]);
                  if ( v > max_v ) {
                      max_v = v;
                      s = i;
                  }
              }
              if ( max_v <= EPS ) {
                  return 0;
              }
              used[s] = true;
              order[k] = s;
              Double m = 1/t[s][k];
              size_t k0 = (k+1)&-4;
              //TIMER(invertGauss work);
              xmul_row(S-k0, t[s]+k0, m);
              forN ( i, n ) {
                  if ( i == s ) continue;
                  xsub_mul_row(S-k0, t[s]+k0, t[i][k], &t[i][k0]);
              }
          }
          {
              //TIMER(invertGauss output);
              forN ( i, n ) {
                  const Double* s = t[order[i]];
                  vector_t& d = m[i];
                  d.assign(s+n, s+2*n);
              }
          }
          return 1;
      }

      inline bool invertMatr(matrix_t& m, matrix_t& t)
      {
          //matrix_t m0 = m, a = m;
          if ( invertGauss(m) ) {
              return 1;
          }
          return 0;
          return svd_inverse_sym(m, t);
      }

      /*
      inline void dot2(const vector_t& a1, const vector_t& a2,
                       size_t i0, const matrix_t& b,
                       vector_t& r1, vector_t& r2)
      {
          for ( size_t i = i0, np = b.size(); i < np; ++i ) {
              if ( i < np-1 ) {
                  dot22(a1, a2, b[i], b[i+1], &r1[i], &r2[i]);
                  ++i;
              }
              else {
                  dot21(a1, a2, b[i], r1[i], r2[i]);
                  break;
              }
          }
      }
      */

inline void addScaledVector(const vector_t& a, Double b, vector_t& r)
      {
	add_mul_row(a.size(), &a[0], b, &r[0]);
      }


//
// class members
//

LogisticModel::LogisticModel(const vector<string>& phenotypes,
			     const vector<string>& genotypes,
			     const vector<double>& covariates)
{
  TIMER(LogisticModel init);
  size_t N = phenotypes.size();
  size_t C = covariates.size()/N;
  size_t P = phenotypes[0].size();
  size_t G = genotypes[0].size();
  // TR(N|G|P|C);
  
  {
    //TIMER(LogisticModel init gen/phe);
    phenotypesT.resize(P, vector<char>(N));
    genotypesT.resize(G, vector<char>(N));
    forN ( i, N ) {
      forN ( j, P ) {
	phenotypesT[j][i] = phenotypes[i][j]-'0';
      }
      forN ( j, G ) {
	genotypesT[j][i] = '2'-genotypes[i][j];
      }
    }
  }
  {
    //TIMER(LogisticModel init ind2cov);
    assert(N*C == covariates.size());
    ind2cov.resize(N);
    forN ( i, N ) {
      ind2cov[i].resize(C);
      forN ( j, C ) {
	ind2cov[i][j]=(covariates[i*C+j]);
      }
    }
  }
  cur_gen = not_set;
  cur_phe = not_set;
  
  np = C+2;
  nind = 0;
  
  t1.reserve(N+4);
  t2.reserve(N+4);
  s1.reserve(np);
  
  all_valid = true;
} // LogisticModel constructor

double LogisticModel::getAssociationStat(int testParameter)
{  
  if ( !all_valid ) return 0;
  double var = getVar(testParameter);
  if ( var < 1e-20 || !realnum(var) ) return 0;
  return SQR(coef[testParameter])/var;
}

double LogisticModel::getVar(int testParameter)
{  
  return S[testParameter][testParameter];
}

vector_t LogisticModel::getVar()
{  
  vector_t var(np);
  forN ( i, np )
    var[i] = S[i][i];
  return var;
}


double LogisticModel::fitLM_phenotype_continue(size_t it0){
              TIMER(fitLM_phenotype_continue);

              xsizeVector(t1, nind);
              s1.resize(np);
              ncoef.resize(np);

              //////////////////////////////////////
              // Newton-Raphson to fit logistic model
              forNF ( it, it0, 20 ) {
                  {
                      TIMER(fitLM_phenotype_continue ncoef);
                      subVector(Y, p, t1);
                      dot(t1, Xt, s1);
                      dot(s1, S, ncoef);
                  }
                  // Update coefficients, and check for 
                  // convergence
                  double delta = 0;
                  {
                      TIMER(fitLM_phenotype_continue coef);
                      forN ( j, np ) {
                          delta += abs(ncoef[j]);
                          coef[j] += ncoef[j];
                      }
                  }
                  if ( delta < DELTA_EPS ) {
                      break;
                  }
      #ifdef NCOEF_EPS
                  if ( it > 1 && abs(ncoef[XKM]/coef[XKM]) < NCOEF_EPS ) {
                      break;
                  }
      #endif
                  // Determine p and V
                  {
                      TIMER(fitLM_phenotype_continue p);
                      p.assign(nind, 0);
                      forN ( i, np ) {
                          addScaledVector(Xt[i], -coef[i], p);
                      }
                  }
                  {
                      TIMER(fitLM_phenotype_continue v);
                      calc_p_v(p, V);
                  }
                  // Update coefficients (equivalent R code in comment)
                  // b <- b + solve( t(X) %*% V %*% X ) %*% t(X) %*% ( y - p ) 
                  {
                      TIMER(fitLM_phenotype_continue S);
                      update_S_V();
                  }
                  {
                      TIMER(fitLM_phenotype_continue S svd);
                      if ( !invertMatr(S, ST) ) {
                          return 0;
                      }
                  }
              }
              return getAssociationStat(XKM);
          }
        


double LogisticModel::fitLM_phenotype_step(size_t k, size_t it){
	    TIMER(fitLM_phenotype_continue);

              xsizeVector(t1, nind);
              s1.resize(np);
              ncoef.resize(np);

              //////////////////////////////////////
              // Newton-Raphson to fit logistic model
              do {
                  {
                      TIMER(fitLM_phenotype_step ncoef);
                      subVector(Y, p, t1);
                      dot(t1, Xt, s1);
                      dot(s1, S, ncoef);
                  }
                  // Update coefficients, and check for 
                  // convergence
                  double delta = 0;
                  {
                      TIMER(fitLM_phenotype_step coef);
                      forN ( j, np ) {
                          delta += abs(ncoef[j]);
                          coef[j] += ncoef[j];
                      }
                  }
                  if ( delta < DELTA_EPS ) {
                      break;
                  }
                  if ( it > 1 ) {
      #ifdef NCOEF_EPS
                      if ( abs(ncoef[XKM]/coef[XKM]) < NCOEF_EPS ) break;
      #endif
                  }
                  // Determine p and V
                  {
                      TIMER(fitLM_phenotype_step p);
                      p.assign(nind, 0);
                      forN ( i, np ) {
                          addScaledVector(Xt[i], -coef[i], p);
                      }
                  }
                  {
                      TIMER(fitLM_phenotype_step v);
                      calc_p_v(p, V);
                  }
                  // Update coefficients (equivalent R code in comment)
                  // b <- b + solve( t(X) %*% V %*% X ) %*% t(X) %*% ( y - p ) 
                  {
                      TIMER(fitLM_phenotype_step S);
                      update_S_V();
                  }
                  {
                      TIMER(fitLM_phenotype_step S svd);
                      if ( !invertMatr(S, ST) ) {
                          return 0;
                      }
                  }
              } while (0);
              return getResult(k);
          }
        

double LogisticModel::fitLM_phenotype_start(size_t k){
              TIMER(fitLM_phenotype_start);

              assert(cur_gen != not_set);
              assert(cur_phe != not_set);
              assert(!Xt.empty());
              assert(!Y.empty());
              assert(!S0.empty());
              S = S0;

              coef.resize(np);
              V.assign(nind, .25);

              xsizeVector(t1, nind);
              s1.resize(np);
              ncoef.resize(np);

              //////////////////////////////////////
              // Newton-Raphson to fit logistic model
              forN ( it, START_STEPS ) {
                  {
                      TIMER(fitLM_phenotype_start ncoef);
                      subVector(Y, .5, t1);
                      xdot(t1, Xt, s1);
                      dot(s1, S, coef);
                  }
                  double delta = 0;
                  forN ( i, np ) delta += abs(coef[i]);
                  if ( delta < DELTAS_EPS ) {
                      break;
                  }
                  // Determine p and V
                  {
                      TIMER(fitLM_phenotype_start p);
                      p.assign(nind, 0);
                      forN ( i, np ) {
                          addScaledVector(Xt[i], -coef[i], p);
                      }
                  }
                  {
                      TIMER(fitLM_phenotype_start v);
                      calc_p_v(p, V);
                  }
                  // Update coefficients (equivalent R code in comment)
                  // b <- b + solve( t(X) %*% V %*% X ) %*% t(X) %*% ( y - p ) 
                  {
                      TIMER(fitLM_phenotype_start S);
                      update_S_V();
                  }
                  {
                      TIMER(fitLM_phenotype_start S svd);
                      if ( !invertMatr(S, ST) ) {
                          return -1e9;
                      }
                  }
              }
              return getResult(k);
              return getAssociationStat(XKM);
          }


double LogisticModel::fitLM_phenotype(){
  TIMER(fitLM phenotype);
  
  assert(cur_gen != not_set);
  assert(cur_phe != not_set);
  assert(!Xt.empty());
  assert(!Y.empty());
  assert(!S0.empty());
  S = S0;
  xsizeVector(t1, nind);
  s1.resize(np);
  ncoef.resize(np);
  
  coef.assign(np, 0);
  p.assign(nind, .5);
  V.assign(nind, .25);
  
  #ifdef DEBUG
  cout << Xt << endl;
  /* for (int i=0; i<nind; i++)
     {
     cout << i << "\t"
     << Y[i] << "\t";
     for (int j=0; j<np; j++)
     cout << Xt[j][i] << "\t";
     cout << "\n";
     }
  */
  cerr << "model size(Y) = " << Y.size() << endl;
  cerr << "model size(Xt) = " << Xt.size() << endl;
  #endif
  //////////////////////////////////////
  // Newton-Raphson to fit logistic model
  {
    forN ( it, 20 ) {
      {
	TIMER(fitLM_phenotype ncoef);
	xsub_row(nind, &Y[0], &p[0], &t1[0]);
	xdot(t1, Xt, s1);
	dot(s1, S, ncoef);
      }
      // Update coefficients, and check for 
      // convergence
      double delta = 0;
      {
	TIMER(fitLM_phenotype coef);
	forN ( j, np ) {
	  delta += abs(ncoef[j]);
	  coef[j] += ncoef[j];
	}
      }
      if ( delta < DELTA_EPS ) {
	break;
      }
#ifdef NCOEF_EPS
      if ( it > 1 && abs(ncoef[XKM]/coef[XKM]) < NCOEF_EPS ) {
	break;
      }
#endif
      // Determine p and V
      {
	TIMER(fitLM_phenotype p);
	p.assign(nind, 0);
	forN ( i, np ) {
	  addScaledVector(Xt[i], -coef[i], p);
	}
      }
      {
	TIMER(fitLM_phenotype v);
	calc_p_v(p, V);
      }
      // Update coefficients (equivalent R code in comment)
      // b <- b + solve( t(X) %*% V %*% X ) %*% t(X) %*% ( y - p )
      {
	TIMER(fitLM_phenotype S);
	update_S_V();
      }
      {
	TIMER(fitLM_phenotype S svd);
	if ( !invertMatr(S, ST) ) {
	  return 0;
	}
      }
    }
  }
  #ifdef DEBUG
  cout << coef << endl;
  cout << getVar() << endl;
  cout << S << endl;
  #endif
  return getAssociationStat(XKM);
}


double LogisticModel::fitLM_phenotype_estimate(size_t k)
	  // this is the default first approximation for problem k=(i*P + j) of M*P regressions
	  // Is this a marker only fit, Y ~ M ?
          {
              TIMER(fitLM_phenotype_estimate);
              assert(cur_gen != not_set);
              assert(cur_phe != not_set);
              assert(!Xt.empty());
              assert(!Y.empty());
              assert(!S0.empty());

              xsizeVector(t1, nind);
              s1.resize(np);

              TIMER(fitLM_phenotype_estimate coef);
              Double nc;
              sub_row(nind, &Y[0], .5, &t1[0]);
              s1.resize(np);
	      // set s1 = t1 . Xt
              xdot(t1, Xt, s1);
	      // set nc = s1 . S0[0]
              dot(s1, S0[XKM], nc);
	      // set s = S[0][0], r = ...  // preparing single column fit Y ~ marker only?
              double s = S0[XKM][XKM], r = SQR(nc)/s;
              // this is never done!
	      if ( 0 ) {
                  Double ss = 0;
		  // for t=0..np accumulate squared S0[0][t] / SQR(s)
                  forN ( t, np ) ss += SQR(S0[XKM][t]);
                  ss /= SQR(s);
                  Double ref = reference.empty()? r: reference[k], err = ref/r-1;
                  Double s1s = accumulate(1+ALL(s1), .0);
                  // TR(k|setprecision(6)|r|nc|s|ss|ref|err|s1[0]|s1s);
              }
	      // debug dump:
#ifdef DEBUG 
	      cout << "==== Start _phenotype_estimate head  for k=" << k << endl;
	      for (int i=0; i<13; i++)
		{
		  cout << i << "\t"
		       << Y[i] << "\t";
		  for (int j=0; j<np; j++)
		    cout << Xt[j][i] << "\t";
		  cout << "\n";
		}
	      cout << "S0 matrix:" << endl;
	      for (int i=0;i<np;i++) {
		for (int j=0;j<np;j++) {
		  cout << S0[i][j] << "\t";
		}
		cout << "\n";
	      }
	      cout << "end S0 matrix" << endl;		
	      cout << "XKM = " << XKM << endl;
	      cout << "==== End _phenotype_estimate head for k=" << k << endl;
#endif // DEBUG
              return r;
          }
          


          double LogisticModel::getResult(size_t k)
          {
              if ( 1 ) {
                  TIMER(getResult);
                  Double nc;
                  xsub_row(nind, &Y[0], &p[0], &t1[0]);
                  dot(t1, Xt, s1);
                  dot(s1, S[XKM], nc);
                  Double c = coef[XKM], s = S[XKM][XKM];
                  Double r = SQR(c)/s;
                  if ( 0 ) {
                      Double ss = 0;
                      forN ( t, np ) ss += SQR(S[XKM][t]);
                      ss /= SQR(s);
                      Double ref = reference.empty()? r: reference[k], err =ref/r-1, est = nc/c;
                      // TR(k|setprecision(6)|r|c|nc|s|ss|ref|err|est);
                  }
                  return r;
              }
              //return SQR(coef[XKM])/S0[XKM][XKM];
              return getAssociationStat(XKM);
          }

          void LogisticModel::fitLM_genotype() {
              TIMER(fitLM genotype);
              // check for multicollinearity; quit if columns of X are near-dependent
              all_valid = np && nind;// initially all valid
              //all_valid = checkVIF();
              if ( !all_valid ) return;
              {
                  TIMER(fitLM_genotype S);
                  sizeMatrix(S0,np,np);
                  update_S0_25();
              }
              {
                  TIMER(fitLM_genotype S svd);
                  if ( !invertMatr(S0, ST) ) {
                      all_valid=false;
                      return;
                  }
              }
              S0_made = 1;
          }
 
          void LogisticModel::update_S_V()
          {
              forN ( i, np ) {
                  xmul_row(nind, &Xt[i][0], &V[0], &t1[0]);
                  if ( i+1 < np ) {
                      xmul_row(nind, &Xt[i+1][0], &V[0], &t2[0]);
                      xdot(t1, Xt[i], S[i][i]);
                      xdot2(t1, t2, i+1, Xt, S[i], S[i+1]);
                      ++i;
                  }
                  else {
                      xdot(t1, i, Xt, S[i]);
                  }
              }
              forN ( i, np ) forNF ( j, i, np ) {
                  S[j][i]=S[i][j];
              }
          }

void LogisticModel::update_S0_25()
{
  forN ( i, np ) {
    dot(Xt[i], i, Xt, S0[i]);
  }
  forN ( i, np ) forNF ( j, i, np ) {
    S0[j][i]=S0[i][j]*=.25;
  }
}

          void LogisticModel::select_genotype(size_t gen)
          {
              if ( gen == cur_gen ) return;
              unselect_genotype();

              TIMER(LogisticModel set_genotype);
              cur_gen = gen;
              const vector<char>& genotypes = genotypesT[cur_gen];
              size_t N = genotypes.size();
              if ( CacheGen* c = get_cache_gen(cur_gen, 0) ) {
                  all_valid = c->all_valid;
                  assert(xgen.empty());
                  swap(c->xgen, xgen);
                  assert(XtXKM.empty());
                  swap(c->XtXKM, XtXKM);
                  assert(S0.empty());
                  swap(c->S0, S0);
                  S0_made = c->S0_made;
              }
              else {
                  all_valid = 1;
                  S0_made = 0;
                  clearVector(xgen);
                  clearVector(XtXKM);
              }
              if ( XtXKM.empty() ) {
                  TIMER(set_genotype init);
                  XtXKM.reserve(N+4);
                  forN ( i, N ) {
                      int markerValue = genotypes[i];
                      if ( markerValue < 0 ) {
                          xgen.push_back(i);
                      }
                      else {
                          XtXKM.push_back(markerValue);
                      }
                  }
                  add_tail_zeros(XtXKM);
              }
              CacheXGen& cx = cache_XGen[xgen];
              assert(Xt.empty());
              swap(cx.Xt, Xt);
              if ( Xt.empty() ) {
                  TIMER(set_genotype Xt);
                  Xt.assign(np, vector_t());
                  forN ( i, np ) Xt[i].reserve(N+4);
                  size_t C = ind2cov[0].size();
                  forN ( i, N ) {
                      int markerValue = genotypes[i];
                      if ( markerValue < 0 ) continue;
                      forN ( j, C ) {
                          Xt[j+2].push_back(ind2cov[i][j]);
                      }
                  }
                  nind=Xt[2].size();
                  forN ( j, C ) add_tail_zeros(Xt[j+2]);
                  Xt[XK1].assign(nind, 1);
                  add_tail_zeros(Xt[XK1]);
              }
              else {
                  nind=Xt[2].size();
              }
              assert(Xt[XKM].empty());
              assert(!XtXKM.empty());
              swap(Xt[XKM], XtXKM);
              if ( !S0_made ) {
                  fitLM_genotype();
              }
          }


void LogisticModel::unselect_genotype(bool save)
          {
              if ( cur_gen == not_set ) return;
              unselect_phenotype(save);
              {
                  CacheXGen& cx = cache_XGen[xgen];
                  assert(!Xt.empty());
                  assert(XtXKM.empty());
                  swap(Xt[XKM], XtXKM);
                  assert(cx.Xt.empty());
                  swap(cx.Xt, Xt);
              }
              if ( save ) {
                  CacheGen& c = *get_cache_gen(cur_gen, 1);
                  assert(Xt.empty());
                  assert(c.xgen.empty());
                  swap(c.xgen, xgen);
                  assert(c.XtXKM.empty());
                  swap(c.XtXKM, XtXKM);
                  c.all_valid = all_valid;
                  assert(c.S0.empty());
                  swap(c.S0, S0);
                  c.S0_made = S0_made;
              }
              else {
                  clearVector(xgen);
                  clearVector(XtXKM);
                  all_valid = 1;
                  S0_made = 0;
              }
              cur_gen = not_set;
          }

         void LogisticModel::select_phenotype(size_t phe)
          {
              if ( phe == cur_phe ) return;
              unselect_phenotype();
              cur_phe = phe;

              TIMER(LogisticModel set_phenotype);
              if ( CachePhe* c = get_cache_phe(cur_gen, cur_phe, 0) ) {
                  all_valid = c->all_valid;
                  assert(coef.empty());
                  swap(c->coef, coef);
                  assert(S.empty());
                  swap(c->S, S);
      #ifdef SAVE_Y
                  assert(Y.empty());
                  swap(c->Y, Y);
      #else
                  clearVector(Y);
      #endif
      #ifdef SAVE_P
                  assert(p.empty());
                  swap(c->p, p);
      #else
                  clearVector(p);
      #endif
              }
              else {
                  clearVector(coef);
                  clearVector(Y);
                  clearVector(p);
              }
              if ( Y.empty() ) {
                  TIMER(set_phenotype Y);
                  const vector<char>& phenotypes = phenotypesT[cur_phe];
                  if ( xgen.empty() ) {
                      Y.assign(ALL(phenotypes));
                  }
                  else {
                      size_t i1 = 0;
                      forIter ( k, xgen ) {
                          size_t i = *k;
                          Y.insert(Y.end(), &phenotypes[i1], &phenotypes[i]);
                          i1 = i+1;
                      }
                      Y.insert(Y.end(), &phenotypes[i1], &*phenotypes.end());
                  }
                  //Y.resize(nind);
                  add_tail_zeros(Y);
              }
              if ( p.empty() && !coef.empty() ) {
                  TIMER(set_phenotype V);
                  assert(!S.empty());
                  p.assign(nind, 0);
                  forN ( i, np ) {
                      addScaledVector(Xt[i], -coef[i], p);
                  }
                  //V.resize(nind);
                  add_tail_zeros(V);
                  calc_p_v(p, V);
              }
          }
          
 
void LogisticModel::unselect_phenotype(bool save)
{
  if ( cur_phe == not_set ) return;
  if ( save ) {
    CachePhe& c = *get_cache_phe(cur_gen, cur_phe, 1);
    c.all_valid = all_valid;
    assert(c.coef.empty());
    swap(c.coef, coef);
    assert(c.S.empty());
    swap(c.S, S);
#ifdef SAVE_Y
    assert(c.Y.empty());
    swap(c.Y, Y);
#endif
#ifdef SAVE_P
    assert(c.p.empty());
    swap(c.p, p);
#endif
  }
  else {
    clearVector(coef);
    clearVector(Y);
    clearVector(p);
  }
  cur_phe = not_set;
}

