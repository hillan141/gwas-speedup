#include <string>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <map>
#include <algorithm>
#include <bitset>

#include "Locus.h"
#include "PlinkIO.h"
#include "options.h"
#include "helper.h"
#include "CSNP.h"
#include "nlist.h"

using namespace std;

void PlinkIO::error( string msg) {
  cerr << "\nERROR: " << msg << "\n";
  exit(1);
}

void PlinkIO::readBinData( string famfile, string bitfilename_map, string bitfilename) {
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

  printLOG("Reading map (extended format) from [ " 
	   + par::bitfilename_map + " ] \n");
  
  vector<Locus> ordered;

  ifstream MAP(par::bitfilename_map.c_str(), ios::in);
  MAP.clear();
  
  int c=0;
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
	  locus.push_back(loc);
	  ordered.push_back(*loc);
	}  
      else
	delete loc;
    }
  
  
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
      
      
      ///////////////////////////
      // SNP-major mode
      
      if (bfile_SNP_major)
	{
	  
	  CSNP * snp;
	  
	  // Outer loop for SNPs
	  int s=0;
	  while (s<locus.size()) // for all SNPs
	    {
	      
	      // Do we want to include this SNP?
	      if ( include[s] > -1 )
		snp = SNP[ include[s] ];	  
	      else
		snp = NULL;
	      
	      // Inner loop for individuals
	      
	      //	  vector<Individual*>::iterator person = sample.begin();
	      int indx = 0;
	      int ss = sample.size();
	      
	      while ( indx < ss )
		{
		  
		  bitset<8> b;
		  //            if ( par::fast_binary )
		  // 		{
		  // 		  b = memblock[indx++];
		  // 		}
		  //	      else
		  //		{
		  char ch[1];
		  BIT.read(ch,1);
		  if (!BIT) 
		    error("Problem with the BED file...has the FAM/BIM file been changed?\n");		
		  b = ch[0];	  
		  //		}
		  
		  int c=0;
		  
		  while (c<7 && indx < ss ) 
		    {
		      if (snp)
			{
			  //  		      snp->one.push_back( b[c++] );
			  //  		      snp->two.push_back( b[c++] );
			  
			  snp->one[indx] = b[c++];
			  snp->two[indx] = b[c++];
			  
			}
		      else
			{
			  c+=2;
			}
		      // 		  ++person;		  
		      ++indx;
		    }
		  
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
      
      char ch[1];
      BIT.read(ch,1);
      if (BIT) 
	error("Problem with the BED file... has the FAM/BIM file been changed?\n");
            
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
  
  // initialize n, used later by readCovListFile
  n = sample.size();
  nl = locus.size();
}

void PlinkIO::readFamFile(string filename) {

  //////////////////////////////
  // Read pedigree information

  // Initially, assume a binary trait
  par::qt = false;
  par::bt = true; 

  printLOG("Reading pedigree information from [ " 
	   + filename + " ] \n");

  checkFileExists( filename );

  ifstream PED;
  PED.open(filename.c_str());
  PED.clear();

  vector<Individual*> ambiguous;
  
  int nmale = 0;
  int nfemale = 0;
  int nambig = 0;

  int c=0;
  while(!PED.eof())
    {

      Individual * person = new Individual;

      // First 6 obligatory fields
      string phenotype;

      PED >> person->fid 
	  >> person->iid 
	  >> person->pat
	  >> person->mat
	  >> person->sexcode
	  >> phenotype;
      
      // Are we using 0/1 coding?
      if (par::coding01) 
	{
	  if ( phenotype == "1" ) 
	    phenotype = "2";      
	  else if ( phenotype == "0" )
	    phenotype = "1";
	  else 
	    phenotype = "0";
	}

      // Skip last empty line that gets read
      if (person->fid=="") 
	{
	  delete person;
	  break;
	}

      // Check for reserved family ID code
      if ( person->fid=="FID" )
	error("FID is a reserved ID... please select a different family ID");

      // Check sex
      if (person->sexcode=="1")
	{
	  person->sex = true; // male
	  nmale++;
	}
      else if (person->sexcode=="2")
	{
	  person->sex = false;  // female (default)
	  nfemale++;
	}
      else 
	{
	  ambiguous.push_back(person);
	  nambig++;
	  if (!par::ignore_missing_sex)
	    person->missing = true;
	}

      
      ///////////////
      // A non-founder?
      
      person->founder = (person->pat == "0" && person->mat == "0") ? true : false;
      
      //////////////////////////////
      // Test for quantitative trait
      
      if (phenotype == par::missing_phenotype)
	person->missing = true;
      else
	{
	  if ( ! from_string<double>( person->phenotype, phenotype, std::dec ) )
	    person->missing = true;
	  else	  
	    if (phenotype != "0" && 
		phenotype != "1" && 
		phenotype != "2" ) 
	      {
		par::qt = true;
		par::bt = false; 
	      }
	}
      
      
      // Increase person counter
      c++;
      
      // Add individual to list
      sample.push_back(person);
      
    } // while !PED.eof()
    
  PED.clear();
  PED.close();
  

  // If a binary trait, now make 0 missing also
  // i.e. if we never saw other than missing, 0, 1 or 2
  
  if (par::bt)
    for (int i=0; i<sample.size(); i++)
      if ( sample[i]->phenotype == 0 )
	sample[i]->missing = true;

  
  printLOG(int2str(c)+" individuals read from [ " 
	   + filename + " ] \n");
  int nm=0;
  for (int i=0;i<sample.size();i++)
    if(!sample[i]->missing) nm++;
  printLOG(int2str(nm) 
	   + " individuals with nonmissing phenotypes\n");
  
  if (par::bt) 
    {

      if (par::coding01) 
	printLOG("Assuming a disease phenotype (0=unaff, 1=aff, other=miss)\n");
      else
	{
	  printLOG("Assuming a disease phenotype (1=unaff, 2=aff, 0=miss)\n");
	  if (par::missing_phenotype!="0")
	    printLOG("Missing phenotype value is also " 
		     + par::missing_phenotype + "\n");
	}

      int ncase = 0;
      int ncontrol = 0;
      int nmissing = 0;

      for (int i=0; i<sample.size(); i++)
	if ( sample[i]->missing )
	  nmissing++;
	else if ( sample[i]->phenotype == 1 )
	  ncontrol++;
	else if ( sample[i]->phenotype == 2 )
	  ncase++;
      printLOG(int2str(ncase)+" cases, "
	       +int2str(ncontrol)+" controls and "
	       +int2str(nmissing)+" missing\n");
      
    }
  else 
    {
      printLOG("Assuming a quantitative trait\n");
      printLOG("Missing phenotype value is " + par::missing_phenotype + "\n");
    }

  // Display sex counts
  printLOG(int2str(nmale)+" males, "+int2str(nfemale)
	   +" females, and "+int2str(nambig)+" of unspecified sex\n");


  // Display list of ambiguously-sexed individuals?
  if (ambiguous.size()>0)
    {
      printLOG("Warning, found "+int2str(ambiguous.size())
	       +" individuals with ambiguous sex codes\n");
      if (!par::ignore_missing_sex)
	printLOG("These individuals will be set to missing ( or use --allow-no-sex )\n");      
      string f = par::output_file_name + ".nosex";
      printLOG("Writing list of these individuals to [ "+f+" ]\n");
      ofstream AMB;
      AMB.open(f.c_str(), ifstream::out);
      for (int i=0; i<ambiguous.size(); i++)
	AMB << ambiguous[i]->fid << "\t" << ambiguous[i]->iid << "\n";
      AMB.close();      
      ambiguous.clear();
    }


}

void PlinkIO::checkFileExists(string f)
{

  ifstream inp;
  
  inp.open(f.c_str(), ifstream::in);
  if(inp.fail())
    {
      inp.clear(ios::failbit);
      inp.close();
      string msg = "No file [ " + f + " ] exists.";
      error(msg);
    }
  inp.close();
  return;

}

void PlinkIO::checkFileExists(vector<string> f)
{
  for (int k=0; k<f.size(); k++)
    checkFileExists(f[k]);
  return;
}

void PlinkIO::printLOG(string s)
{
  cerr << s;
  cerr.flush();
  
  if (!par::silent)
    {
      cout << s;
      cout.flush();
    }
}

bool PlinkIO::openBinaryFile(string s, ifstream & BIT)
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


void PlinkIO::setMarkerRange()
{

  // If chromosome code >0, implies a specific chromosome
  if (par::run_chr>0 && !par::position_window)
    {
      // Get first and last markers on this chromosome
      vector<int> m = getChromosomeMarkerRange((*this),par::run_chr);
      if(m[0]==-1 || m[1]==-1) error("--chr {chromosome} not found:"+int2str(par::run_chr));
      par::run_start = m[0];
      par::run_end = m[1];
    }
  else if (par::position_window)
    {
      // Physical position specified (chromosome and range)
      par::m1 = leftWindowEdge(*this, par::run_chr, par::from_window);
      par::m2 = rightWindowEdge(*this, par::run_chr, par::to_window);
      par::run_start = getMarkerNumber( (*this), par::m1 );
      par::run_end = getMarkerNumber( (*this), par::m2 );       
    }
  else
    {
      // Two SNPs specified (or a SNP and a range)

      // If a specific range on one chromosome is specified
      par::run_start = getMarkerNumber((*this),par::m1);
      par::run_end = getMarkerNumber((*this),par::m2);
      
      if (par::run_start==-1) error("--from {marker} not found");
      if (par::run_end==-1) error("--to {marker} not found");    

      // Do we require a window around a specific SNP?
      if ( par::run_start == par::run_end && par::window > 0 )
        {
          vector<int> m = getWindowRange( *this , par::run_start ); 
          par::run_start = m[0];
          par::run_end = m[1];
        }

  
      if (getMarkerChromosome((*this),par::m1) != getMarkerChromosome((*this),par::m2))
        {
          string msg = "--from {marker} and --to {marker} must lie on same chromosome";
          msg += "\nwhereas these lie on chromosomes "+int2str(getMarkerChromosome((*this),par::m1));   
          msg += " and "+int2str(getMarkerChromosome((*this),par::m2));
          error(msg);
        }
    }

  // Get order right
  if (par::run_start > par::run_end)
    {
      int tmp = par::run_start;
      par::run_start = par::run_end;
      par::run_end = tmp;
    }
  int ccode = locus[par::run_start]->chr;
  printLOG("Scan region on chromosome " + int2str(ccode)
           + " from [ " + locus[par::run_start]->name 
           + " ] to [ " + locus[par::run_end]->name 
           + " ]\n");
  
}

bool PlinkIO::readCovariateFile()
{

  // This will set individuals as missing 
  // if they have a missing value for the 
  // covariate, or do not appear in the file
  
  checkFileExists(par::covar_filename);
  ifstream COV(par::covar_filename.c_str(), ios::in);
  
  map<string,Individual*> uid;
  map<string,Individual*>::iterator ii;
  set<Individual*> hasCovariate;

  for (int i=0; i<sample.size(); i++)
    {
      uid.insert(make_pair(sample[i]->fid+"_"+sample[i]->iid,sample[i]));      
    }
  
  int nvalid=0;
  while (!COV.eof())
    {
      string pfid, piid, cov;
      
      char cline[par::MAX_LINE_LENGTH];
      COV.getline(cline,par::MAX_LINE_LENGTH,'\n');
      
      // convert to string
      string sline = cline;
      if (sline=="") continue;

      string buf; 
      stringstream ss(sline); 
      vector<string> tokens; 
      while (ss >> buf)
	tokens.push_back(buf);
      
      // Trying to read past last column 
      if (tokens.size() < 2+par::mult_covar) 
      {

	if (! par::qmatch )
	  {
	    for (int i0=0; i0<tokens.size(); i0++)
	      printLOG(tokens[i0]+" ");
	    printLOG("\n");
	  }

	COV.close();
        return false;
      } 
      
      pfid = tokens[0];
      piid = tokens[1];
      cov = tokens[1 + par::mult_covar];

      
      ii = uid.find(pfid+"_"+piid);
      if (ii != uid.end() )
	{
	  Individual * person = ii->second;


	  // Set covariate value
	  bool badValue = false;
	  if ( ! from_string<double>( person->covar, cov, std::dec ) )
	    badValue = true;
	  
	  // Note that we've seen a covariate for this individual
	  hasCovariate.insert(person);
	  
	  // Was this missing? 
	  if (cov == par::missing_phenotype || badValue )
	    {
	      person->missing = true;
	    }
	  else
	    nvalid++;
	  
	  
	}
    }
  COV.close();


  // Set to missing any individuals for who we did not see the covariate

  vector<Individual*>::iterator person = sample.begin();
  while ( person != sample.end() )
    {
      if ( hasCovariate.find( *person ) == hasCovariate.end() )
	(*person)->missing = true;
      person++;
    }
  
  
  printLOG("Reading covariate from [ " + par::covar_filename + " ] with ");
  printLOG("nonmissing values for "+int2str(nvalid)+" individuals\n");
  
  return true;
}

bool PlinkIO::readCovListFile()
{

  // This will set individuals as missing if they have a missing value
  // for the covariate, or do not appear in the file
  
  // if number of individuals not initialized, before, do it now
  if (n==0) {
    n = sample.size();
    cerr << "warning: number of individuals not initialized before calling readCovListFile" << endl;
  }

  ifstream COV(par::clist_filename.c_str(), ios::in);
  
  
  map<string,Individual*> uid;
  map<string,Individual*>::iterator ii;
  set<Individual*> hasCovariate;

  for (int i=0; i<sample.size(); i++)
    {
      uid.insert(make_pair(sample[i]->fid+"_"+sample[i]->iid, sample[i]));
    }

  // If need (for later selection) keep explicit track of what is missing
  map<Individual*, vector<bool> > isMissing;
  map<Individual*,bool> originalPersonMissingStatus;

  par::clist_number = -1;

  int nvalid=0;
  while (!COV.eof())
    {

      string pfid, piid;
      vector_t clist;

      char cline[par::MAX_LINE_LENGTH];
      COV.getline(cline,par::MAX_LINE_LENGTH,'\n');
      
      // convert to string
      string sline = cline;
      if (sline=="") continue;
      
      string buf; 
      stringstream ss(sline); 
      vector<string> tokens; 
      while (ss >> buf)
	tokens.push_back(buf);
      
      if ( par::clist_number < 0 ) 
	{
	  par::clist_number = tokens.size() - 2;
	  // Assign default headers (can be overwritten)
	  clistname.resize(par::clist_number);
	  for (int c=0; c<par::clist_number; c++)
	    clistname[c] = "COV"+int2str(c+1);	  
	}
      else if (tokens.size() != par::clist_number + 2 ) 
	{
	  printLOG("Line:\n"+sline+"\n");
	  COV.close();
	  return false;
	} 
      
      pfid = tokens[0];
      piid = tokens[1];
      
      ii = uid.find(pfid+"_"+piid);
      if (ii != uid.end() )
	{
	  Individual * person = ii->second;
	  
	  // Track individual covariate missing status
	  person->clistMissing.resize(par::clist_number);
	  
	  // Store original missingness status for this person
	  originalPersonMissingStatus.insert(make_pair( person, person->missing ));

	  vector<bool> missing_status;
	  
	  // Were any missing/bad values?
	  bool okay = true;

	  // Add covariate values to clist
	  person->clist.clear();
	  for (int c=2; c<par::clist_number+2; c++)
	    {
	      double t = 0;
	      if ( ! from_string<double>( t, tokens[c], std::dec ) )
		okay = false;
	      person->clist.push_back( t );
	    }
	  
	  // Note that we've seen a covariate for this individual
	  hasCovariate.insert(person);
	  	  
	  for (int c=0; c<par::clist_number; c++)
	    {
	      if ( tokens[c+2] == par::missing_phenotype )
		{
		  okay = false;
		  missing_status.push_back(false);
		  person->clistMissing[c] = true;
		}
	      else
		{
		  missing_status.push_back(true);
		  person->clistMissing[c] = false;
		}
	    }
	  
	  if (!okay)
	    person->missing = true;
	  else
	    nvalid++;

	  // Record, if we will use this below
	  if ( par::clist_selection )
	    isMissing.insert(make_pair( person, missing_status ));
	  
	}
      else if ( pfid == "FID" && piid == "IID" )
	{
	  // This is a header row -- read in covariate names
	  for (int c=0; c<par::clist_number; c++)
	    clistname[c] = tokens[c+2];
	}
    }
  COV.close();


  
  
  // Set to missing any individuals for who we did not see the covariate
  // But also fill their covariate list with missing values

  vector<Individual*>::iterator person = sample.begin();
  vector<bool> dummy_missing_status(n,false);
  while ( person != sample.end() )
    {
      if ( hasCovariate.find( *person ) == hasCovariate.end() )
	{
	  (*person)->missing = true;
	  (*person)->clist.clear();
	  (*person)->clist.resize(par::clist_number, -9 );
	  (*person)->clistMissing.clear();
	  (*person)->clistMissing.resize(par::clist_number, true );
	  
	  if ( par::clist_selection )
	    isMissing.insert(make_pair( (*person), dummy_missing_status ));
	}
      person++;
    }
  

  printLOG("Reading " 
	   + int2str(par::clist_number) 
	   + " covariates from [ " 
	   + par::clist_filename 
	   + " ] with ");
  printLOG("nonmissing values for "
	   +int2str(nvalid)
	   +" individuals\n");
  

  /////////////////////////////////////////////////////////
  // Do we actually want to keep all these covariates?
  
  if ( par::clist_selection_number || par::clist_selection_name )
    {
      vector<int> covlist;
      if ( par::clist_selection_number )
	{
	  NList nl(par::clist_number);
	  covlist = nl.deparseNumberList(par::clist_selection_string);
	}
      else
	{
	  map<string,int> mapping;
	  for (int c=0; c<par::clist_number; c++)
	    mapping.insert(make_pair( clistname[c],c));
	  NList nl(par::clist_number);
	  covlist = nl.deparseStringList(par::clist_selection_string,&mapping);
	}
      
      int nvalid = 0;
      for (int i=0; i<n; i++)
	{
	  Individual * person = sample[i];

	  vector_t tmp = person->clist;
	  vector<bool> tmpMissing = person->clistMissing;

	  person->clist.clear();
	  person->clistMissing.clear();

	  // Reset per-person missing code
	  person->missing = originalPersonMissingStatus.find( person )->second;
	  
	  vector<bool> missing_status = isMissing.find( person )->second;

	  bool okay = true;
	  for (int c=0; c<covlist.size(); c++)
	    {
	      person->clist.push_back( tmp[ covlist[c] ] );
	      person->clistMissing.push_back( tmpMissing[ covlist[c] ] );

	      if ( ! missing_status[covlist[c]] )
		{ 
		  person->missing = true;
		  okay = false; 
		}
	    }
	  
	  if ( okay ) nvalid++;
	}

      // Reset sample-wide values (names, number)

      vector<string> tmp = clistname;
      clistname.clear();
      for (int c=0; c<covlist.size(); c++)
	clistname.push_back( tmp[ covlist[c] ] );

      printLOG("Selected subset of " + int2str(covlist.size()) + " from " 
	       + int2str(par::clist_number) + " covariates\n");
      
      par::clist_number = covlist.size();
      
      printLOG("For these, nonmissing covariate values for "
	       +int2str(nvalid)+" individuals\n");
  

    }
  return true;
} // readCovListFile
