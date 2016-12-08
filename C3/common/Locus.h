#ifndef __LOCUS_H__
#define __LOCUS_H__

#include <string>

using namespace std;

class Locus {
 public:
  Locus() { chr=0; name=""; allele1=""; allele2=""; freq=0; pos=0; bp=0; nm=0; }
  ~Locus() { 
    // cerr << "Destroying locus: " << name << endl;
  }
  int chr;
  string name;
  string allele1;
  string allele2;

  double freq;     // of allele1
  double pos;      // cM map positions
  int bp;          // base-pair position
  int nm;          // number of non-missing alleles

  // Copy constructor
  Locus(const Locus& h1) { copy(h1); }
  Locus & operator= (const Locus & h1) { copy(h1); return *this; }
  
  void copy(const Locus &h1)
    {
      chr = h1.chr;
      name = h1.name;
      allele1 = h1.allele1;
      allele2 = h1.allele2;
      freq = h1.freq;
      pos = h1.pos;
      bp = h1.bp;
      nm = h1.nm;
    }
  
  bool operator< (const Locus & p2) const
    {
      return (chr < p2.chr || (chr == p2.chr && bp < p2.bp) );
    }
  
  bool operator== (const Locus & p2) const
    {
      return ( name == p2.name );
    }


};


#endif

