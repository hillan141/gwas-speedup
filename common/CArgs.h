#ifndef __CARGS_H__
#define __CARGS_H__

#include <string>
#include <vector>

using namespace std;

class CArgs
{
 public:
  CArgs(int,char**);

  int count()
    { return n; }
  
  bool any()
    { return n > 1 ? true : false; }
  
  void fromScript(string);
  void fromPriorLog(string);
  bool find(string);
  string value(string);
  int value_int(string);
  double value_double(string);
  long unsigned int value_lui(string);
  // void check_unused_options(Plink &);
  bool parseOptions(string,string);

  vector<string> value(string,int);
  vector<string> varValue(string);
  
  vector<string> a;
  
 private:
  int n;
  vector<bool> parsed;
  vector<bool> option;
  vector<string> root_command;
  vector<string> original;
  map<string,string> optionLabel;
};

#endif
