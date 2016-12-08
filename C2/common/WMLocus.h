#ifndef __WMLOCUS__H
#define __WMLOUCS_H

#include <string>
#include <vector>

using namespace std;

class WMLocus {
 public:
  WMLocus() { chr=0; name=""; }

  void reset()
    {
      allele.clear();
      weight.clear();
    }

  int chr;
  string name;
  
  vector<string> allele;
  vector<double> weight;
};

#endif
