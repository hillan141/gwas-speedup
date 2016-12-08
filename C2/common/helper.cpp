#include <iostream>
#include <ostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <cmath>
#include <algorithm>
#include <cerrno>
#include <cstdio>
#include <ctime>

#include "helper.h"
#include "options.h"

using namespace std;

string bool_as_text(bool b)
{
    std::stringstream converter;
    converter << b;
    return converter.str();
}

string int2str(int n)
{
  std::ostringstream s2( std::stringstream::out );
  s2 << n;
  return s2.str();
}

void shutdown()
{

  time_t curr=time(0);
  string tdstamp = ctime(&curr);

  if (!par::silent) cout << "\nAnalysis finished: " + tdstamp +"\n";
  cout << "\nAnalysis finished: " + tdstamp +"\n";
  exit(0);
}

vector<int> getChromosomeMarkerRange(PlinkIO & P, int c)
{
  vector<int> m(2);
  m[0] = -1; // first
  m[1] = -1; // last
  
  for (int i=0;i<P.locus.size();i++)
    {
      if (P.locus[i]->chr==c) 
        {
          if (i<=m[0] || m[0]==-1) m[0]=i;
          if (i>=m[1] || m[1]==-1) m[1]=i;        
        }      
    }
  return m;  
}

string leftWindowEdge(PlinkIO & P, int chr, int bp)
{
  // Get nearest SNP 
  Locus * marker = NULL;
  int distance = -1;
    
  vector<Locus*>::iterator loc = P.locus.begin();
  while ( loc != P.locus.end() )
    {
      if ( (*loc)->chr == chr) 
	if ( (*loc)->bp >= bp )
	  if ( (*loc)->bp - bp < distance || ! marker )
	    {
	      distance =  (*loc)->bp - bp;
	      marker = *loc;
	    }
      loc++;
    }

  if (!marker) 
    error("Could not place marker for left window edge");

  return marker->name;  

}


string rightWindowEdge(PlinkIO & P, int chr, int bp)
{

  // Get nearest SNP 
  Locus * marker = NULL;
  int distance = -1;
  
  vector<Locus*>::iterator loc = P.locus.begin();
  while ( loc != P.locus.end() )
    {
      if ( (*loc)->chr == chr) 
	if ( (*loc)->bp <= bp )
	  if ( bp - (*loc)->bp < distance || ! marker )
	    {
	      distance =  bp - (*loc)->bp;
	      marker = *loc;
	    }
      loc++;
    }

  if (!marker) 
    error("Could not place marker for right window edge");

  return marker->name;  
}

int getMarkerNumber(PlinkIO & P, string m)
{
  for (int i=0;i<P.locus.size();i++)
    if (P.locus[i]->name==m) return i;
  return -1;  
}

vector<int> getWindowRange(PlinkIO &P, int s)
{
  vector<int> m(2);
  m[0] = s; // first SNP
  m[1] = s; // last SNP
  
  // move backwards
  int x=s;
  int chr=P.locus[s]->chr;
  int bp=P.locus[s]->bp;
  int win = (int)(par::window * 1000); // half window size in bases
  int nl = P.locus.size() - 1;

  while ( 1 ) 
    {
      if ( x==0 ) break;
      // Move one position, until on different chromosome, or outside window
      x--;
      if ( P.locus[x]->chr != chr ) { x++; break; }
      if ( bp - P.locus[x]->bp > win ) { x++; break; }
    }
  m[0]=x;
  x=s;
  while ( 1 ) 
    {
      if ( x== nl ) break;
      x++;
      if ( P.locus[x]->chr != chr ) { x--; break; }
      if ( P.locus[x]->bp - bp > win ) { x--; break; }
    }
  m[1]=x;
  return m;
}

int getMarkerChromosome(PlinkIO & P, string m)
{
  for (int i=0;i<P.locus.size();i++)
    if (P.locus[i]->name==m) return P.locus[i]->chr;
  return -1;
}

void error(string msg)
{
  cerr << "\nERROR: " << msg << "\n";
  cerr << "\nERROR: " << msg << "\n";  
  exit(1);
}
