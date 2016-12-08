#ifndef __CSNP_H__
#define __CSNP_H__

#include <vector>

using namespace std;

// Main genotype storage, ordered by SNP
class CSNP
{
 public:

  vector<bool> one; // SNP-major mode genotypes
  vector<bool> two; 
  
  char getChar(int i) {
    return '2' - (one[i] + two[i]);
  }
  
  string getString() {
    string ret = "";
    for (int i = 0; i < one.size(); ++ i) {
        ret += '2' - (one[i] + two[i]);
    }
    return ret;
  }
  
  void getString(string &ret) {
	int vsize = one.size();   
    ret.assign(vsize,'2');
    for (int i = 0; i < vsize; ++ i) {
        ret[i] -=  one[i] + two[i];
    }    
  }
  
  void getString(char ret[]) {
    for (int i = 0; i < one.size(); ++ i) {
        ret[i] = '2' - (one[i] + two[i]);
    }
  }
  
  void fill(char* &buffer, int ss) {
    vector<bool>::iterator iter_1 = one.begin();
    vector<bool>::iterator iter_2 = two.begin();
    vector<bool>::iterator iter_3 = one.end();
    vector<bool>::iterator iter_4 = iter_1 +4*(one.size()/4);
    for (;iter_1 != iter_4;) {
        char ch = *buffer;
        ++ buffer;
	*(iter_1 ++) = (ch & 1);
	*(iter_2 ++) = (ch & 2);
	*(iter_1 ++) = (ch & 4);
	*(iter_2 ++) = (ch & 8);
	*(iter_1 ++) = (ch & 16);
	*(iter_2 ++) = (ch & 32);
	*(iter_1 ++) = (ch & 64);
	*(iter_2 ++) = (ch & 128);
    }
    if(iter_1 != iter_3) {
	char ch = *buffer;
	++ buffer;
	*(iter_1 ++) = (ch & 1);
	*(iter_2 ++) = (ch & 2);	
	if (iter_1 != iter_3) {
	    *(iter_1 ++) = (ch & 4);
	    *(iter_2 ++) = (ch & 8);
		if (iter_1 != iter_3) {
		    *(iter_1 ++) = (ch & 16);
		    *(iter_2 ++) = (ch & 32);
		}
	}
    }    	   
  }
};

inline bool operator <(const CSNP &a, const CSNP &b)
{
    return a.one < b.one || a.one == b.one && a.two < b.two;
}

#endif

