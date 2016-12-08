#ifndef __PLINKIO_H__
#define __PLINKIO_H__

#include <string>
#include <vector>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <map>
#include <algorithm>
#include <bitset>
#include "Locus.h"
#include "Individual.h"
#include "CSNP.h"

using namespace std;

class PlinkIO {
  vector<Locus> ordered;
  // checkFileExists(string fileName);
 public:
  vector<Locus*> locus;
  vector<CSNP*> SNP;
  vector<Individual*> sample;
  vector<string> clistname;
  int n;
  int nl;
 PlinkIO() : ordered(), locus(), SNP(), sample(), clistname(), n(0), nl(0) {
    cerr << "Constructing PlinkIO object..." << endl;
    ordered.resize(0);
    locus.resize(0);
    SNP.resize(0);
    sample.resize(0);
    clistname.resize(0);
    n=nl=0;
  }

  ~PlinkIO() {
    cerr <<  "Destroying PlinkIO object..." << endl;
  }

  bool openBinaryFile(string s, ifstream & BIT);
  void readBinData( string famfile, string bitfilename_map, string bitfilename);
  void readFamFile(string famfile);
  void checkFileExists(string f);
  void checkFileExists(vector<string> f);
  void printLOG(string s);
  void error(string msg);
  void setMarkerRange();
  bool readCovariateFile();
  bool readCovListFile();
};

#endif

