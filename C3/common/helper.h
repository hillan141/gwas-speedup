#ifndef __HELPER_H__
#define __HELPER_H__

#include <string>
#include <iostream>
#include <vector>
#include "PlinkIO.h"

using namespace std;

template <class T>
bool from_string(T& t, 
                 const std::string& s, 
                 std::ios_base& (*f)(std::ios_base&))
{
  std::istringstream iss(s);
  return !(iss >> f >> t).fail();
}

string int2str(int n);
void shutdown();
vector<int> getChromosomeMarkerRange(PlinkIO & P, int c);
string leftWindowEdge(PlinkIO & P, int bp, int chr);
string rightWindowEdge(PlinkIO & P, int bp, int chr);
int getMarkerNumber(PlinkIO &,string);
vector<int> getWindowRange(PlinkIO &P, int);
int getMarkerChromosome(PlinkIO &,string);
void error(string msg);
string bool_as_text(bool b);

#endif
