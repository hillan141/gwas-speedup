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
  
	locus.reserve(1<<20);

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
	  locus.push_back(loc); 
      else
	delete loc;
    }
    
    vector<Locus> ordered(locus.size());
    for(int i=0, vsize=ordered.size(); i<vsize; ++i)
		ordered[i] = *locus[i];
  
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


  for (int i=0, vsize=locus.size(); i<vsize; i++)
    {
      // swap file order into locus position
      // and make all same chromosome
      Locus * loc = &ordered[i];
      loc->bp = (int)loc->freq;
      loc->chr = 1;
      // keep track of genetic order, but only 
      loc->freq = i;
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
      int vsize = ordered.size();
      include.resize(vsize);
      for (int j=0; j<vsize; j++)
	{
	  int fp = (int)ordered[j].freq;
	  
	  if ( fp < par::run_start || fp > par::run_end ) 
	    include[j] = -1;
	  else
	 {
	   include[j] = fp;
	   nl_actual++;
	 }
	}
      

      //              0  1  2  3  4  5  6   7  8  9
      // We now have -1 -1 -1  3  4  5  6  -1 -1 -1
      // but we want -1 -1 -1  0  1  2  3  -1 -1 -1
      

      for (int j=0, vsize=ordered.size(); j<vsize; j++)
	{
	  if ( include[j] > -1 ) 
	    include[j] -= par::run_start ;
	}
      
    }
  else
    {
		int vsize = ordered.size();
      include.resize(vsize);
      // If we do want to look at all the data
      for (int j=0; j<vsize; j++)
		include[j] = (int)ordered[j].freq;
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
		int vsize = sample.size();
		SNP.resize(nl_actual);
	  for (int i=0; i<nl_actual; i++)
	    {
	      CSNP * newlocus = new CSNP;
	      newlocus->one.resize( vsize, false );
	      newlocus->two.resize( vsize, false );
	      SNP[i] = newlocus;
	    }
	}
      else
	{
	  vector<Individual*>::iterator person = sample.begin();
	  vector<Individual*>::iterator iend = sample.end();
	  while ( person != iend )
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
	  
	  // Outer loop for SNPs
	  const int each = (sample.size() * 2 + 7) / 8;
	  const int length = each * locus.size();
	  char* buffer = new char[length + 1];
	  const char* backup = buffer;
	  BIT.read(buffer, length);

	  if (!BIT) 
		    error("Problem with the BED file...has the FAM/BIM file been changed?\n");
	  for (int s = 0, vsize=locus.size(), ssize=sample.size(); s < vsize; ++ s) // for all SNPs
	    {
	      // Do we want to include this SNP?
	      int id = include[s];
	      if ( id == -1 ) {
            buffer += each;
	        continue;
	      }
	      SNP[id]->fill(buffer, ssize);
	    }
	  delete [] backup;
	  // Set file mode
	  par::SNP_major = true;
	}
      
      
      ////////////////////////////////////
      // Individual-major mode
      
      else
	{
	  
	  // Outer loop for individuals
	  vector<Individual*>::iterator person = sample.begin();
	  vector<Individual*>::iterator iend = sample.end();
	  int vsize = locus.size();
	  while ( person != iend )
	    {
	      
	      // Inner loop for SNPs
	      int s=0;
	      while (s<vsize) // for all SNPs
		{
		  
		  char ch[1];
		  BIT.read(ch,1);
		  if (!BIT) error("Problem with the BED file... has the FAM/BIM file been changed?\n");
		  
		  bitset<8> b;
		  b = ch[0];	  
		  
		  int c=0;
		  
		  while (c<7 && s<vsize)
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
      for(int l=0, vsize=locus.size(); l < vsize; l++)
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

#include <time.h>
int PlinkIO::readBinData2( string famfile, string bitfilename_map,
			   string bitfilename, vector<string> &genotypes ) {
  
  // We cannot assume the file will be in order, as it might have been 
  // previusly created by a --merge/--bmerge command

  // Check files exist
  clock_t begin, end;
  double time_spent;  
  begin = clock();

  
  checkFileExists(par::famfile);
  checkFileExists(par::bitfilename_map);
  checkFileExists(par::bitfilename);


    // Open the bit file now to know
    ifstream BIT;
    bool bfile_SNP_major = openBinaryFile(par::bitfilename, BIT);
    
  end = clock();
  time_spent = (double)(end-begin)/CLOCKS_PER_SEC;  
  cerr << "	Sub init time: " << time_spent << endl;
  begin = clock();
    
    vector<Locus> ordered;
    if( !bfile_SNP_major || ((!par::plink) && (!par::run_chr==0)) || par::do_not_load_snps ) {
      error("wrong set of options");
    }

  printLOG("Reading map (extended format) from [ " 
	   + par::bitfilename_map + " ] \n");
  
  ifstream MAP(par::bitfilename_map.c_str(), ios::in);
  MAP.clear();

  int locus_size=0;
  {
    int loc_int;
    string loc_name;
    string loc_string;   
    double loc_pos;               
    while(true)
      {
	MAP >> loc_int   // will automatically by numeric
	    >> loc_name 
	    >> loc_pos   
	    >> loc_int     
	    >> loc_string
	    >> loc_string;
	
	if ( MAP.eof() )
	    break;
	else if ( MAP.fail() )
	    error("Problem reading BIM file, line " + 
	      int2str(locus_size+1) + "\n");
	else if(loc_name=="" )
	    continue;

	// Check that cM/M specification looks correct, if 
	// we want to perform a plink-based analysis
	if (par::plink && (!par::cm_map) && (loc_pos > 50) )
	  error("Looks like you need to specify --cm ??");
	
	locus_size++;
      } // for each locus
      //~ string loc_line;
      //~ while (std::getline(MAP, loc_line))
        //~ locus_size++;
  
  } 

  printLOG( int2str(locus_size) 
	    + " markers to be included from [ " +
	    par::bitfilename_map + " ]\n");
  
  MAP.clear();
  MAP.close();

  if ( locus_size == 0 ) 
    shutdown();

  end = clock();
  time_spent = (double)(end-begin)/CLOCKS_PER_SEC;  
    cerr << "	Sub map parsing time: " << time_spent << endl;
  begin = clock();  
 
  //////////////////////////////
  // Read individual information
  readFamFile(par::famfile);

  end = clock();
  time_spent = (double)(end-begin)/CLOCKS_PER_SEC;  
  cerr << "	Sub fam parsing time: " << time_spent << endl;

  begin = clock();  
  
  printLOG("Reading genotype bitfile from [ " 
	   + par::bitfilename + " ] \n");
  
  //~ ifstream BIT;
  //~ bool bfile_SNP_major = openBinaryFile(par::bitfilename, BIT);
      
      /////////////////////////////////////////
      // Read entire genotype data into a temp
      // array -- not do not use as default method
      // as it does not speed things up that much
      // but uses twice the memory
      
         //~ vector<char> memblock;
      //~ 
      //~ 
           //~ {
             //~ ifstream::pos_type fbegin = BIT.tellg();
             //~ BIT.seekg(0, ios::end);
             //~ ifstream::pos_type fend = BIT.tellg();
             //~ ifstream::pos_type size = fend-fbegin+1;
             //~ memblock.resize(size);
             //~ BIT.seekg(fbegin);
             //~ BIT.read(&memblock[0], size);
             //~ ////~ BIT.close();            
           //~ }
      
      //~ ///////////////////////////
      //~ // SNP-major mode
      //~ // Outer loop for SNPs
  const int each = (sample.size() * 2 + 7) / 8; //number of bytes to store one 2-bit genotype for all N individuals
  const int length = each * locus_size; // total bytes to store all N*M genotypes
  char* buffer = new char[length + 1];
  const char* backup = buffer;
  BIT.read(buffer, length); // read the entire bed file in SNP-major order
  if (!BIT) 
    error("Problem with the BED file...has the FAM/BIM file been changed?\n");

      //~ char* buffer = &memblock[0];
      genotypes.clear();
      genotypes.resize(locus_size);
      for (int s = 0, ssize=sample.size(); s < locus_size; ++ s) // for all SNPs
      {
	genotypes[s].assign(ssize,'2'); // sample.size assigments of '2' to initialize all individuals at s-th genotype
	string& cur_genotype = genotypes[s];
	int ssize4 = 4*(ssize/4);   // rounded number of individuals
	int i=0;  
	for (;i<ssize4;) {
	  char ch = *buffer;
	  ++ buffer;
	  cur_genotype[i++] -= (ch & 1) + ((ch & 2)>>1);
	  cur_genotype[i++] -= ((ch & 4)>>2) + ((ch & 8)>>3);
	  cur_genotype[i++] -= ((ch & 16)>>4) + ((ch & 32)>>5);
	  cur_genotype[i++] -= ((ch & 64)>>6) + ((ch & 128)>>7);
	}
	if(i != ssize) { // rem(individuals, 4), need to tack these 1-3 individuals on the end
	      char ch = *buffer;
	      ++ buffer;
	      cur_genotype[i++] -= (ch & 1) + ((ch & 2)>>1);
	      if(i != ssize) {
		  cur_genotype[i++] -= ((ch & 4)>>2) + ((ch & 8)>>3);
		  if(i != ssize) {
		      cur_genotype[i++] -= ((ch & 64)>>6) + ((ch & 128)>>7);
		}
	    }
	  }    	   
      }  
      //~ delete [] backup;
      // Set file mode
      par::SNP_major = true;
      
      // Check that we got what we expected
      
      char ch[1];
      BIT.read(ch,1);
      if (BIT) 
	error("Problem with the BED file... has the FAM/BIM file been changed?\n");
            
      BIT.clear();
      BIT.close();
      
  end = clock();
  time_spent = (double)(end-begin)/CLOCKS_PER_SEC;
  cerr << "	Sub bit parsing time: " << time_spent << endl;
  begin = clock();  
  // Free any buffer memory used
     //~ if ( par::fast_binary )
       //~ memblock.clear();
  
  // initialize n, used later by readCovListFile
  n = sample.size();
  nl = locus_size;
  return locus_size;
} // readBinData2

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
    for (int i=0, vsize=sample.size(); i<vsize; i++)
      if ( sample[i]->phenotype == 0 )
	sample[i]->missing = true;

  
  printLOG(int2str(c)+" individuals read from [ " 
	   + filename + " ] \n");
  int nm=0;
  for (int i=0, vsize=sample.size();i<vsize;i++)
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

      for (int i=0, vsize=sample.size(); i<vsize; i++)
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

  for (int i=0, vsize=sample.size(); i<vsize; i++)
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


// Convert dataset from Individual-major to SNP-major format

void PlinkIO::Ind2SNP()
{
  printLOG("Converting data to SNP-major format\n");

  SNP.clear();
  
  vector<Individual*>::iterator person = sample.begin();
  
  // Initialise SNP positions per person
  while ( person != sample.end() )
    {
      (*person)->i1 = (*person)->one.end()-1;
      (*person)->i2 = (*person)->two.end()-1;
      person++;
    }
  
  // Copy, per SNP
  int l = 0;
  int nl_all = locus.size();
  while ( l < nl_all  )
    {
 
      CSNP * newlocus = new CSNP;
      
      person = sample.begin();
       
      while ( person != sample.end() )
	{
	  // Add genotype to SNP-major storage
	  newlocus->one.push_back( *((*person)->i1) );
	  newlocus->two.push_back( *((*person)->i2) );
	  
	  // Shift one SNP back
	  (*person)->i1--;
	  (*person)->i2--;

	  // Remove individual-major storage
	  (*person)->one.pop_back();
	  (*person)->two.pop_back();
	  
	  // Advance to next person
	  person++;
	}
            
      // And add this new SNP to the main list
      SNP.push_back(newlocus);
      
      // Next SNP
      l++;
    }
 
  // We finally need to reverse the order of these 
  reverse(SNP.begin(), SNP.end());

  par::SNP_major = true;
}

// Convert dataset from SNP-major format to Individual-major
void PlinkIO::SNP2Ind()
{
  
  printLOG("Converting data to Individual-major format\n");
    
  vector<Individual*>::iterator person = sample.begin();
  
  // Make sure these containers are empty
  while ( person != sample.end() )
    {
      (*person)->one.clear();
      (*person)->two.clear();
      person++;
    }
  
  ///////////////////////////////
  // Iterate over SNPs
  
  vector<CSNP*>::iterator s = SNP.begin();
 
  while ( s != SNP.end() )
    {
            
      /////////////////////////////
      // Iterate over individuals
      
      vector<bool>::iterator i1 = (*s)->one.begin();
      vector<bool>::iterator i2 = (*s)->two.begin();
      vector<Individual*>::iterator person = sample.begin();
      
      while ( person != sample.end() )
	{
	  
	  // Add SNP alleles
	  (*person)->one.push_back(*i1);
	  (*person)->two.push_back(*i2);
	  
 	  // Shift one SNP back
 	  i1++;
 	  i2++;
	  	  
 	  // Advance to next person
 	  person++;
 	}
      
      // For this SNP, remove SNP-major storage completely
      delete (*s);
      
      // Next SNP
      s++;      
    }
  
  SNP.clear();
  par::SNP_major = false;
}

int PlinkIO::readBinDataInterval( string famfile, string bitfilename_map,
				  string bitfilename, vector<string> &genotypes, vector<string> &locus_name, int firstMarker, int lastMarker ) {
  
  if (firstMarker < 0 || lastMarker < 0 || lastMarker < firstMarker) {
    error("invalid marker range in readBinDataInterval");
  }
  int numMarkers = lastMarker - firstMarker + 1;

  // We cannot assume the file will be in order, as it might have been 
  // previusly created by a --merge/--bmerge command

  // Check files exist
  clock_t begin, end;
  double time_spent;  
  begin = clock();

  
  checkFileExists(par::famfile);
  checkFileExists(par::bitfilename_map);
  checkFileExists(par::bitfilename);


    // Open the bit file now to know
    ifstream BIT;
    bool bfile_SNP_major = openBinaryFile(par::bitfilename, BIT);
    
  end = clock();
  time_spent = (double)(end-begin)/CLOCKS_PER_SEC;  
  cerr << "	Sub init time: " << time_spent << endl;
  begin = clock();
    
    vector<Locus> ordered;
    if( !bfile_SNP_major || ((!par::plink) && (!par::run_chr==0)) || par::do_not_load_snps ) {
      error("wrong set of options");
    }

  printLOG("Reading map (extended format) from [ " 
	   + par::bitfilename_map + " ] \n");
  
  ifstream MAP(par::bitfilename_map.c_str(), ios::in);
  MAP.clear();

  int locus_size=0;
  int loc_curr_index = 0;
  {
    int loc_int;
    string loc_name;
    string loc_string;   
    double loc_pos;               

    locus_name.resize(0);
    while(true)
      {
	MAP >> loc_int   // will automatically by numeric
	    >> loc_name 
	    >> loc_pos   
	    >> loc_int     
	    >> loc_string
	    >> loc_string;
	
	if ( MAP.eof() )
	    break;
	else if ( MAP.fail() )
	    error("Problem reading BIM file, line " + 
	      int2str(locus_size+1) + "\n");
	else if(loc_name=="" )
	    continue;

	// Check that cM/M specification looks correct, if 
	// we want to perform a plink-based analysis
	if (par::plink && (!par::cm_map) && (loc_pos > 50) )
	  error("Looks like you need to specify --cm ??");
	if (loc_curr_index >= firstMarker && loc_curr_index <= lastMarker) {
	  locus_size++;
	  locus_name.push_back(loc_name);
	}
	loc_curr_index++;
      } // for each locus
      //~ string loc_line;
      //~ while (std::getline(MAP, loc_line))
        //~ locus_size++;
  
  } 

  printLOG( int2str(locus_size) 
	    + " markers on range [" + int2str(firstMarker) + "-" + int2str(lastMarker) + "] to be included from [ " +
	    par::bitfilename_map + " ]\n");
  
  MAP.clear();
  MAP.close();

  if ( locus_size == 0 ) 
    shutdown();

  end = clock();
  time_spent = (double)(end-begin)/CLOCKS_PER_SEC;  
    cerr << "	Sub map parsing time: " << time_spent << endl;
  begin = clock();  
 
  //////////////////////////////
  // Read individual information
  readFamFile(par::famfile);

  end = clock();
  time_spent = (double)(end-begin)/CLOCKS_PER_SEC;  
  cerr << "	Sub fam parsing time: " << time_spent << endl;

  begin = clock();  
  
  printLOG("Reading genotype bitfile from [ " 
	   + par::bitfilename + " ] \n");
  
  //~ ifstream BIT;
  //~ bool bfile_SNP_major = openBinaryFile(par::bitfilename, BIT);
      
      /////////////////////////////////////////
      // Read entire genotype data into a temp
      // array -- not do not use as default method
      // as it does not speed things up that much
      // but uses twice the memory
      
         //~ vector<char> memblock;
      //~ 
      //~ 
           //~ {
             //~ ifstream::pos_type fbegin = BIT.tellg();
             //~ BIT.seekg(0, ios::end);
             //~ ifstream::pos_type fend = BIT.tellg();
             //~ ifstream::pos_type size = fend-fbegin+1;
             //~ memblock.resize(size);
             //~ BIT.seekg(fbegin);
             //~ BIT.read(&memblock[0], size);
             //~ ////~ BIT.close();            
           //~ }
      
      //~ ///////////////////////////
      //~ // SNP-major mode
      //~ // Outer loop for SNPs
  const long each = (sample.size() * 2 + 7) / 8; //number of bytes to store one 2-bit genotype for all N individuals
  const long length = each * locus_size; // total bytes to store all N*M genotypes
  char* buffer = new char[length + 1];
  const char* backup = buffer;
  cerr << "locus_size= " << locus_size << " sample.size= " << sample.size() << " each =" << each << " length= " << length << endl;
  // need to offset to firstMarker, then read to lastMarker
  BIT.seekg( each * firstMarker, ios::cur);
  if (!BIT) {
    error("failed to seek to offset position = " + each*firstMarker);
  }
  cerr << "seek-ed to position = " << each*firstMarker << endl;
  BIT.read(buffer, length); // read the entire bed file in SNP-major order
  if (BIT) {
    cerr << "all characters read successfully." << endl;
  } else {
    cerr << "error: only " << BIT.gcount() << " chars could be read" << endl;
    error("Problem with the BED file...has the FAM/BIM file been changed?\n");
  }
      //~ char* buffer = &memblock[0];
      genotypes.clear();
      genotypes.resize(locus_size);
      for (int s = 0, ssize=sample.size(); s < locus_size; ++ s) // for all SNPs
      {
	genotypes[s].assign(ssize,'2'); // sample.size assigments of '2' to initialize all individuals at s-th genotype
	string& cur_genotype = genotypes[s];
	int ssize4 = 4*(ssize/4);   // rounded number of individuals
	int i=0;  
	for (;i<ssize4;) {
	  char ch = *buffer;
	  ++ buffer;
	  cur_genotype[i++] -= (ch & 1) + ((ch & 2)>>1);
	  cur_genotype[i++] -= ((ch & 4)>>2) + ((ch & 8)>>3);
	  cur_genotype[i++] -= ((ch & 16)>>4) + ((ch & 32)>>5);
	  cur_genotype[i++] -= ((ch & 64)>>6) + ((ch & 128)>>7);
	}
	if(i != ssize) { // rem(individuals, 4), need to tack these 1-3 individuals on the end
	      char ch = *buffer;
	      ++ buffer;
	      cur_genotype[i++] -= (ch & 1) + ((ch & 2)>>1);
	      if(i != ssize) {
		  cur_genotype[i++] -= ((ch & 4)>>2) + ((ch & 8)>>3);
		  if(i != ssize) {
		      cur_genotype[i++] -= ((ch & 64)>>6) + ((ch & 128)>>7);
		}
	    }
	  }    	   
      }  
      //~ delete [] backup;
      // Set file mode
      par::SNP_major = true;
      
      // Check that we got what we expected
      /*
      char ch[1];
      BIT.read(ch,1);
      if (BIT) 
	error("Problem with the BED file... has the FAM/BIM file been changed?\n");
      */
      
      BIT.clear();
      BIT.close();
      
  end = clock();
  time_spent = (double)(end-begin)/CLOCKS_PER_SEC;
  cerr << "	Sub bit parsing time: " << time_spent << endl;
  begin = clock();  
  // Free any buffer memory used
     //~ if ( par::fast_binary )
       //~ memblock.clear();
  
  // initialize n, used later by readCovListFile
  n = sample.size();
  nl = locus_size;
  return locus_size;
} // readBinData2
