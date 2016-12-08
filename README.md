# gwas-speedup
## Code for GWAS logistic regression project

Date: 6 October 2016

This directory includes source code associated with the publication
"Stepwise Distributed Open Innovation Contests for Software Development - Acceleration of Genome-Wide Association Analysis".

__LICENSE__

All source code here is released under the GNU General Public License Version 2 (GPL2).  See COPYING.txt for a copy of the license.

__PREREQUISITES__

Codes have been tested on Red Hat Enterprise Linux 6 with gcc 4.4.7, 32 GB RAM.  We cannot guarantee compilation or execution in other environments. 

__CONTENTS__

The codes are contained in subdirectories.  See the manuscript for more detailed descriptions of each code.

C1: A top-scoring TopCoder contest submission, amended to read PLINK binary filesets
C2: C1, amended to accelerate the initial reading of PLINK binary data
C3: C2, with additional multithreading
plink-flr: PLINK 1.07 source code, amended to include the contest-winning logistic regression code

An example PLINK fileset is included in the t/ subdirectory, containing synthetic data for N=6000 individuals, M=7000 markers, P=1 phenotype and C=5 covariates.

__RUNNING__

A 'test' makefile target is provided that will execute the codes on the example PLINK fileset.

To run the codes:

For C1:
cd C1 
make test 

For C2: 
cd C2/venco 
make test 

For C3: 
cd C3/venco 
make test 

For plink-flr: 
cd plink-flr 
make test  


