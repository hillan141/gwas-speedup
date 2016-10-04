

//////////////////////////////////////////////////////////////////
//                                                              //
//           PLINK (c) 2005-2008 Shaun Purcell                  //
//                                                              //
// This file is distributed under the GNU General Public        //
// License, Version 2.  Please see the file COPYING for more    //
// details                                                      //
//                                                              //
//////////////////////////////////////////////////////////////////


#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <map>
#include <algorithm>
#include <bitset>

#include "plink.h"
#include "options.h"
#include "helper.h"

extern ofstream LOG;

void Plink::readBinData()
{

  if ( par::do_not_load_snps )
    {
      printLOG("Skipping SNP and genotype information...\n");
      checkFileExists(par::famfile);
      readFamFile(par::famfile);
      return;
    }

  // We cannot assume the file will be in order, as it might have been 
  // previusly created by a --merge/--bmerge command

  // Check files exist

  checkFileExists(par::famfile);
  checkFileExists(par::bitfilename_map);
  checkFileExists(par::bitfilename);

  // extract selected SNPs
  checkFileExists(par::extract_file);
  ifstream INFILE(par::extract_file.c_str(),ios::in);
  INFILE.clear();
  printLOG("Reading list of SNPs ");
  printLOG("to extract [ "
		 + par::extract_file + " ] ... ");

  map<string, int> selectedLocusIndex;
  int cs = 0;
  while (!INFILE.eof())
    {
      string m;
      INFILE >> m;
      if (m=="") continue;
      selectedLocusIndex.insert(make_pair(m, -1));
      cs++;
    }
  INFILE.close();
  printLOG(int2str(cs)+" read\n");

// end extract selected SNPs

  

printLOG("Reading map (extended format) from [ " 
	   + par::bitfilename_map + " ] \n");

  vector<Locus> ordered;

  ifstream MAP(par::bitfilename_map.c_str(), ios::in);
  MAP.clear();
  
  int c=0;
  map<string, int>::iterator sli;
  while(!MAP.eof())
    {
           
      Locus * loc = new Locus;
      
      MAP >> loc->chr   // will automatically by numeric
	  >> loc->name 
	  >> loc->pos   
	  >> loc->bp     
	  >> loc->allele1
	  >> loc->allele2;
      
      if ( MAP.eof() )
	{
	  delete loc; 
	  continue;
	}
      else if ( MAP.fail() )
	{
	  delete loc;
	  error("Problem reading BIM file, line " + int2str(c+1) + "\n");
	}

      // Use the frequency slot temporarily to 
      // store order information
      loc->freq = c++;

      // Check that cM/M specification looks correct, if 
      // we want to perform a plink-based analysis
      if (par::plink && (!par::cm_map) && (loc->pos > 50) )
	error("Looks like you need to specify --cm ??");
      
      // Convert cM to M map distances
      if (par::cm_map) loc->pos /= 100;
      
      
      // Always included, but not always in correct order

      if (loc->name!="") 
	{
	  sli = selectedLocusIndex.find(loc->name);
	  if (sli != selectedLocusIndex.end()) {
	    sli->second = (int)loc->freq;
	    locus.push_back(loc);
	    ordered.push_back(*loc);
	  }
	}  
      else
	delete loc;
    }
  
  /* for (sli=selectedLocusIndex.begin(); sli != selectedLocusIndex.end(); sli++) {
    cout << "selectedLocus: " << sli->first << " -> " << sli->second << endl;
  }
  */

  printLOG( int2str(locus.size()) 
	    + " markers to be included from [ " +
	    par::bitfilename_map + " ]\n");
  
  MAP.clear();
  MAP.close();

  
  if ( locus.size() == 0 ) 
    shutdown();
  
  ///////////////////////////////////////////////
  // Build ordered table, so that genotypes can 
  // be inserted in correct order; then swap locus 
  // file over

  // Sort vector of pointers to Locus
  stable_sort(locus.begin(),locus.end(),less<Locus*>());
  
  // Sort normal vector Locus
  stable_sort(ordered.begin(),ordered.end());
  
  c=0;
  for (int i=0; i<locus.size(); i++)
    {
      // swap file order into locus position
      // and make all same chromosome
      ordered[i].bp = (int)ordered[i].freq;
      ordered[i].chr = 1;
      // keep track of genetic order, but only 
      ordered[i].freq = c++;
    }

  // resort to get lookup table
  stable_sort(ordered.begin(),ordered.end());  


  ///////////////////////////////////////
  // Do we want to look at all the data?
  
  vector<int> include(0);
  int nl_actual = locus.size();

  if ( (!par::plink) && (!par::run_chr==0) ) 
    {

      // Get range
      setMarkerRange();						
      
      // And set to not import all markers outside range
      // 0..nl_all scale: par::run_start..par::run_end

      nl_actual = 0;
      
      for (int j=0; j<ordered.size(); j++)
	{
	  int fp = (int)ordered[j].freq;
	  
	  if ( fp < par::run_start || fp > par::run_end ) 
	    include.push_back(-1);
	  else
	 {
	   include.push_back(fp);
	   nl_actual++;
	 }
	}
      

      //              0  1  2  3  4  5  6   7  8  9
      // We now have -1 -1 -1  3  4  5  6  -1 -1 -1
      // but we want -1 -1 -1  0  1  2  3  -1 -1 -1
      

      for (int j=0; j<ordered.size(); j++)
	{
	  if ( include[j] > -1 ) 
	    include[j] -= par::run_start ;
	}
      
    }
  else
    {
      // If we do want to look at all the data
      for (int j=0; j<ordered.size(); j++)
	include.push_back((int)ordered[j].freq);
    }
 
  

  //////////////////////////////
  // Read individual information

  readFamFile(par::famfile);

  
  printLOG("Reading genotype bitfile from [ " 
	   + par::bitfilename + " ] \n");
  
  ifstream BIT;
  bool bfile_SNP_major = openBinaryFile(par::bitfilename, BIT);
      
      
      /////////////////////////////////////////
      // Read entire genotype data into a temp
      // array -- not do not use as default method
      // as it does not speed things up that much
      // but uses twice the memory
      
      //   vector<char> memblock;
      
      //   if ( par::fast_binary )
      //     {
      //       ifstream::pos_type fbegin = BIT.tellg();
      //       BIT.seekg(0, ios::end);
      //       ifstream::pos_type fend = BIT.tellg();
      //       ifstream::pos_type size = fend-fbegin+1;
      //       memblock.resize(size);
      //       BIT.seekg(fbegin);
      //       BIT.read(&memblock[0], size);
      //       BIT.close();            
      //     }
      
      
      //////////////////////////////
      // Allocate space for SNPs
      
  printLOG("Start allocation of SNP (resize)\n");

      if (bfile_SNP_major)
	{
	  for (int i=0; i<nl_actual; i++)
	    {
	      CSNP * newlocus = new CSNP;
	      newlocus->one.resize( sample.size() );
	      newlocus->two.resize( sample.size() );
	      SNP.push_back(newlocus);
	    }     
	}
      else
	{
	  vector<Individual*>::iterator person = sample.begin();
	  while ( person != sample.end() )
	    {
	      (*person)->one.resize(nl_actual);
	      (*person)->two.resize(nl_actual);
	      person++;
	    }      
	}
      printLOG("Done allocation of SNP (resize)\n");
      
      ///////////////////////////
      // SNP-major mode

      printLOG("Using collapsed IO\n");
      
      if (bfile_SNP_major)
	{
	  
	  CSNP * snp;
	  int ss = sample.size();
	  const long each = (sample.size() * 2 + 7) / 8; // number of bytes to store one 2-bit genotype for all N individuals
	  const long length = each * 1; // total bytes to store all N genotypes for 1 locus
	  char* buffer = new char[length + 1];
	  char *buffer_base = buffer;

	  // Outer loop for SNPs
	  int s=0;
	  while (s<locus.size()) // for all SNPs
	    {    
	      // Do we want to include this SNP?
	      if ( include[s] > -1 ) {
		snp = SNP[ include[s] ];	  
		// seek and load snp with data
		sli = selectedLocusIndex.find(ordered[s].name);
		if (sli != selectedLocusIndex.end()) {
		  int snp_index = sli->second;
		  BIT.seekg( (each * snp_index) + 3, ios::beg);
		  // cerr << "For " << ordered[s].name << " index= " << snp_index << " seeking to pos = " << each*snp_index + 3 << endl;
		  // Inner loop for individuals within the  outer locus      
		  int indx = 0;
		  buffer = buffer_base;
		  BIT.read(buffer, length);  // all N individuals for 1 locus
		  if (!BIT) {
		    cerr << "error: only " << BIT.gcount() << " chars could be read" << endl;
		    error("Problem with the BED file...has the FAM/BIM file been changed?\n");
		  }
		  // cerr << "Read s=" << s << " length=" << length << endl;
		  while ( indx < ss )
		    {
		      //bitset<8> b;
		      //b = *buffer;
		      char ch = *buffer;
		      buffer++; 
		  int c=0;
		  while (c<7 && indx < ss ) 
		    {
		      if (snp)
			{
			  snp->one[indx] = ch & 1; // b[c++];
			  ch >>= 1;
			  snp->two[indx] = ch & 1; // b[c++];
			  ch >>= 1;
			  c+=2;
			}
		      else
			{
			  ch >>= 2;
			  c+=2;
			}
		      ++indx;
		    }
		    } // end individual loop	  
		} else {
		  error("Cannot find the current SNP in the selectedLocusIndex\n");
		}
	      }  else {
		// snp is not to be included, do nothing
		snp = NULL;
	      }
	      // next SNP
	      s++;
	    }
	  // Set file mode
	  par::SNP_major = true;
	}
      
      
      ////////////////////////////////////
      // Individual-major mode
      
      else
	{
	  
	  // Outer loop for individuals
	  vector<Individual*>::iterator person = sample.begin();
	  while ( person != sample.end() )
	    {
	      
	      // Inner loop for SNPs
	      int s=0;
	      while (s<locus.size()) // for all SNPs
		{
		  
		  char ch[1];
		  BIT.read(ch,1);
		  if (!BIT) error("Problem with the BED file... has the FAM/BIM file been changed?\n");
		  
		  bitset<8> b;
		  b = ch[0];	  
		  
		  int c=0;
		  
		  while (c<7 && s<locus.size())
		    {
		      if ( include[s] > -1 )
			{
			  (*person)->one[ include[s] ] = b[c++];
			  (*person)->two[ include[s] ] = b[c++];	      
			}
		      else
			{
			  c+=2;
			}
		      s++;
		    }	 	  
		}
	      
	      person++;
	    }
	  
	  // Set file mode
	  par::SNP_major = false;
	}
      
      // Check that we got what we expected
      /*
      char ch[1];
      BIT.read(ch,1);
      if (BIT) 
	error("Problem with the BED file... has the FAM/BIM file been changed?\n");
      */   
      BIT.clear();
      BIT.close();
      
      
    
  // Free any buffer memory used
  //   if ( par::fast_binary )
  //     memblock.clear();
  
  
  ////////////////////////////////////////
  // If need be, now prune the MAP file 
  // i.e. if --chr or --from/--to were used
  
  if ( (!par::plink) && (!par::run_chr==0) )
    {
      
      vector<Locus*> l0(0);
      for(int l=0; l < locus.size(); l++)
	{    
	  if ( !( l < par::run_start || l > par::run_end ) ) 
	    l0.push_back(locus[l]);
	  else
	    delete locus[l];
	}
      locus.clear();
      locus = l0;
    }
  
}




bool Plink::openBinaryFile(string s, ifstream & BIT)
{

  BIT.open(s.c_str(), ios::in | ios::binary);

  // 1) Check for magic number
  // 2) else check for 0.99 SNP/Ind coding
  // 3) else print warning that file is too old
  
  char ch[1];
  BIT.read(ch,1);
  bitset<8> b;
  b = ch[0];	  
  
  bool bfile_SNP_major = false;
  bool v1_bfile = true;

  // If v1.00 file format
  // Magic numbers for .bed file: 00110110 11011000 = v1.00 bed file
  if (   ( b[2] && b[3] && b[5] && b[6] ) && 
       ! ( b[0] || b[1] || b[4] || b[7] )    )
    {

     // Next number
     BIT.read(ch,1);
     b = ch[0];	  
     if (   ( b[0] && b[1] && b[3] && b[4] ) && 
          ! ( b[2] || b[5] || b[6] || b[7] )    )
      {
        // Read SNP/Ind major coding
        BIT.read(ch,1);
        b = ch[0];	  
        if ( b[0] ) bfile_SNP_major = true;
        else bfile_SNP_major = false;

        if (bfile_SNP_major) 
  	  printLOG("Detected that binary PED file is v1.00 SNP-major mode\n");
        else
	  printLOG("Detected that binary PED file is v1.00 individual-major mode\n");

      } else v1_bfile = false;
      
    } else v1_bfile = false;


  // Reset file if < v1
  if ( ! v1_bfile ) 
   {
    printLOG("Warning, old BED file <v1.00 : will try to recover...\n");
    printLOG("  but you should --make-bed from PED )\n");
    BIT.close();
    BIT.clear();
    BIT.open(s.c_str(), ios::in | ios::binary);
    BIT.read(ch,1);
    b = ch[0];	  
  }

  // If 0.99 file format
  if ( (!v1_bfile) && ( b[1] || b[2] || b[3] || b[4] || b[5] || b[6] || b[7] ) )
    {
      printLOG("\n *** Possible problem: guessing that BED is < v0.99      *** \n");
      printLOG(" *** High chance of data corruption, spurious results    *** \n");
      printLOG(" *** Unles you are _sure_ this really is an old BED file *** \n");
      printLOG(" *** you should recreate PED -> BED                      *** \n\n");

      bfile_SNP_major = false;
      BIT.close();
      BIT.clear();
      BIT.open(s.c_str(), ios::in | ios::binary);
    }
  else if ( ! v1_bfile ) 
    {
      if ( b[0] ) bfile_SNP_major = true;
      else bfile_SNP_major = false;

      printLOG("Binary PED file is v0.99\n");
      
      if (bfile_SNP_major) 
	printLOG("Detected that binary PED file is in SNP-major mode\n");
      else
	printLOG("Detected that binary PED file is in individual-major mode\n");
    }

 return bfile_SNP_major;

}
