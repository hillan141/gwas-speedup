#ifndef __COMMON_H__
#define __COMMON_H__

using namespace std;

#include <cstdlib>
#include <vector>

#define TIME_LIMIT 10

#include <stdint.h>
#include <cstdio>
#include <cstring>
#include <unistd.h>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <algorithm>
#include <cassert>
#include <string>
#include <cmath>
#include <map>
#include <set>
#include <numeric>
#include <strstream>
#include <sstream>
#include <xmmintrin.h>
#ifdef __SSE3__
# include <pmmintrin.h>
#endif
#include <sys/time.h>

# include <pmmintrin.h>

#define sequential 0
#define by_estimate 1
#define by_start 2
#define by_estimate_start1 3
#define by_estimate_start2 4
#define processing_order by_estimate


#define NOP() do{}while(0)

#ifndef CPU_FREQ
# define CPU_FREQ 3.4e9
#endif
#ifndef SPEED
# define SPEED 1
#endif
#define STRINGx(a) #a
#define STRING(a) STRINGx(a)
#define NAME2x(a,b) a##b
#define NAME2(a,b) NAME2x(a,b)
#define ALL(C) (C).begin(), (C).end()
#define forIter(I,C) for(typeof((C).end()) I=(C).begin(); I!=(C).end(); ++I)
#define forNF(I,F,C) for(size_t I=(F); I<(C); ++I)
#define forN(I,C) forNF(I,0,C)
#define ALIGN16 __attribute__((aligned(16)))

#ifndef HOME_RUN
# define NTRACE
# ifdef assert
#  undef assert
# endif
# define assert(m) NOP()
#endif

// #define NTR xNTR
// int NTR = 200000;
struct c_UU { 
  mutable std::ostrstream str;
  c_UU(const char*f, int l, const char*v) { 
    str << f << ':' << l << ' ';
    while ( char c = *v++ ) 
      str.put(c == '|'? ' ': c); 
  }
  ~c_UU() { 
    str.put('\n'); 
    size_t c = str.pcount();
    const char* s = str.str(); 
    str.freeze(0); 
    ::write(2, s, c); 
  }
  template<class X> const c_UU& operator|(const X&x) const
  { str << ' ' << x; return*this; }
};

/*
#ifdef FTRACE
# define fTR(v) (c_UU(__FILE__,__LINE__,#v)|'='|v)
#else
# define fTR(v) NOP()
#endif
#define xTR(v) do{if(--::NTR>=0)fTR(v);else ++::NTR;}while(0)
#ifdef NTRACE
# define TR(v) NOP()
#else
# define TR(v) xTR(v)
#endif
*/

typedef unsigned long long u64;
typedef u64 cpu_clock_t; 
struct clock_t2 { 
  cpu_clock_t clock; 
  unsigned call; 
};

/* extern clock_t2 start_clock; 
extern unsigned clock_calls;
      inline cpu_clock_t get_clock2(){
          cpu_clock_t c; asm volatile("rdtsc":"=A"(c):); return c;
      }
      inline clock_t2 get_clock(){
          clock_t2 c;
      #ifdef CLOCK
          c.clock = get_clock2(); c.call = ++clock_calls;
      #else
          c.clock = 0; c.call = 0;
      #endif
          return c;
      }
      cpu_clock_t operator-(const clock_t2& t1, const clock_t2& t2) {
          cpu_clock_t c = t1.clock-t2.clock, w = (t1.call-t2.call)*11ULL;
          return c > w? c-w: 0;
      }
      void init_clock(){get_clock();clock_calls=0;start_clock=get_clock();}
      struct Clock { cpu_clock_t c; Clock():c(get_clock()-start_clock) {} };
      ostream& operator<<(ostream& out,const Clock& c) {
          return out<<fixed<<setprecision(6)<<c.c*(1./CPU_FREQ);
      }
*/
      /////////////////////////////////////////////////////////////////////////////

      /////////////////////////////////////////////////////////////////////////////

/*
      struct TimerS{cpu_clock_t clc;unsigned cnt;const char*f;int l;const 
      char*n;};
      ostream& operator<<(ostream& out,TimerS& s){
          cpu_clock_t t=get_clock()-start_clock; out<<s.f<<':'<<s.l<<':';
          int l=s.l, w = 30+(l<=9)+(l<=99)+(l<=999)+(l<=9999);
          out<<setw(w)<<s.n<<':'<<fixed<<setw(7)<<setprecision(2)<<1e2*s.clc/t;
          out<<"% "<<s.clc*(SPEED/CPU_FREQ)<<" "<<s.cnt;
          s.clc=s.cnt=0; return out;
      }
*/      

template<class V>struct aligned_alloc{typedef V value_type;
typedef size_t size_type;typedef ptrdiff_t difference_type;
typedef V*pointer;typedef const V*const_pointer;
typedef V&reference;typedef const V&const_reference;
template<class U>struct rebind{typedef aligned_alloc<U>other;};
aligned_alloc() throw(){}aligned_alloc(const aligned_alloc&) throw(){}
template<class U>aligned_alloc(const aligned_alloc<U>&) throw(){}
~aligned_alloc()throw(){}
V*allocate(size_t c,const void* =0){return(V*)_mm_malloc(sizeof(V)*c, 
							   16);}
void deallocate(V* p,size_t ){_mm_free(p);}
size_t max_size()const throw(){return ~size_t(0)/sizeof(V);}
void construct(V* p,const V&v){::new((void *)p)V(v);}
void destroy(V* p){p->~V();}
}; // aligned_alloc


template<class V,class U>inline
bool operator==(const aligned_alloc<V>&,const aligned_alloc<U>&){return 
true;}
template<typename V,typename U>inline
bool operator!=(const aligned_alloc<V>&,const aligned_alloc<U>&){return 
false;}

#if defined(CLOCK)
struct Timer{Timer(TimerS&s);~Timer();TimerS&s;clock_t2 t;};
TimerS* clocks[99];size_t clocks_cnt;Timer::~Timer(){ 
  s.clc+=get_clock()-t; }
Timer::Timer(TimerS&s):s(s){if(!s.cnt++)clocks[clocks_cnt++]=&s;t=get_clock();} // Timer::Timer
struct out_clocks { out_clocks() { init_clock(); } ~out_clocks() {
    forN(i,clocks_cnt)if(clocks[i]->cnt)cerr<<*clocks[i]<<'\n';clocks_cnt=0;
}}; // out_clocks
# define TIMER(n) static TimerS NAME2(clock_,__LINE__)=		 \
  {0,0,__FILE__,__LINE__,#n};					 \
  Timer NAME2(timer_,__LINE__)(NAME2(clock_,__LINE__))
#define TOTAL_TIMER(name) out_clocks NAME2(clocks_,__LINE__); TIMER(name)
#else
#define TIMER(name) NOP()
#define TOTAL_TIMER(name) NOP()
#endif  // CLOCK

typedef float Double;      
typedef vector<Double, aligned_alloc<Double> > vector_t;
typedef vector<vector_t> matrix_t;


/* inline double get_time() {
  timeval tv; gettimeofday(&tv, 0); return tv.tv_sec+tv.tv_usec*1e-6;
}
inline double get_run_time() {
  double t = get_time()-start_time;
#ifdef HOME_RUN
  t *= SPEED;
#endif
  return t;
}
*/

extern void error(const char* s);

      //==============from helper.h======================
      template<class T>
      inline const T SQR(const T a)
      {
          return a*a;
      }

      template<class T>
      inline const T MAX(const T &a, const T &b)
      {
          return b > a ? (b) : (a);
      }

      template<class T>
      inline const T MIN(const T &a, const T &b)
      {
          return b < a ? (b) : (a);
      }

      template<class T>
      inline const T SIGN(const T &a, const T &b)
      {
          return b >= 0 ? (a >= 0 ? a : -a) : (a >= 0 ? -a : a);
      }

      template<class T>
      inline void SWAP(T &a, T &b)
      {
          T dum=a; a=b; b=dum;
      }


#endif // __COMMON_H__
