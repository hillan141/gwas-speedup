CURRDIR = ../venco
export CXXFLAGS = -I$(CURRDIR) -O2 -fno-operator-names -march=core2 -mssse3 -mfpmath=sse -ffast-math -Wno-deprecated -DMAIN -DHOME_RUN# -fomit-frame-pointer
CC = $(CXX)
VPATH = ../common

SMALL_DATADIR = ../../data/small
SMALL_TESTBED = $(SMALL_DATADIR)/UC_N100_M5000_P1_C7
SMALL_PFILE = UC_N100_M5000_P1_C7.assoc.logistic.pval

MID_DATADIR = ../../data/mid
MID_TESTBED = $(MID_DATADIR)/UC_N100_M50000_P1_C7
MID_PFILE = UC_N100_M50000_P1_C7.assoc.logistic.pval
LARGE_TESTBED = $(LARGE_DATADIR)/UC_N100_M700000_P1_C7
LARGE_PFILE = UC_N100_M700000_P1_C7.assoc.logistic.pval

VLARGE_DATADIR = ../../../../t # ../../data/vlarge
# TESTBED = ../../../../t/UC_N6000_M7000_P1_C5/plink
TESTBED = ../../t/UC_N6000_M7000_P1_C5/plink
VLARGE_TESTBED = $(VLARGE_DATADIR)/UC_N6678_M645863_P1_C7
VLARGE_PFILE = UC_N6678_M645863_P1_C7.assoc.logistic.pval
VLARGE_URL = https://docs.google.com/uc?export=download&id=0B0Z89c1A1-hTMjNERXduOGU1Y1E

NEWIO_COMMON = ../../../alegro3x/src/common

OBJS = GWASSpeedup.o ../common/read_bed.o ../common/PlinkIO.o ../common/options.o ../common/nlist.o ../common/helper.o ../common/dcdflib.o

OBJS_NEWIO = GWASSpeedup.o $(NEWIO_COMMON)/read_bed.o $(NEWIO_COMMON)/PlinkIO.o $(NEWIO_COMMON)/options.o $(NEWIO_COMMON)/nlist.o $(NEWIO_COMMON)/helper.o $(NEWIO_COMMON)/dcdflib.o

GWASSpeedup.o : GWASSpeedup.cpp GWASSpeedup.h

gspeed_bed_newio : $(OBJS_NEWIO)
	$(CXX) $(CXXFLAGS) -o gspeed_bed_newio $(OBJS_NEWIO)

gspeed_bed : $(OBJS)
	$(CXX) $(CXXFLAGS) -o gspeed_bed $(OBJS)

test : gspeed_bed
	./gspeed_bed -f 0 -l 6999 -b $(TESTBED)

smalltest : gspeed_bed
	./gspeed_bed -f 0 -l 4999 -b $(SMALL_TESTBED) > small.out
	python checkDiff.py small.out $(SMALL_DATADIR)/$(SMALL_PFILE)

midtest : gspeed_bed
	./gspeed_bed -f 0 -l 49999 -b $(MID_DATADIR)/$(MID_TESTBED) > mid.out
	python checkDiff.py mid.out $(MID_DATADIR)/$(MID_PFILE)

largetest : gspeed_bed
	./gspeed_bed -f 0 -l 699999 -b $(LARGE_DATADIR)/$(LARGE_TESTBED) > large.out
	python checkDiff.py large.out $(LARGE_DATADIR)/$(LARGE_PFILE)

vlargetest : gspeed_bed
	@echo "download the v. large test set from '$(VLARGE_URL)'"
	@echo "and extract the contents to $(VLARGE_DATADIR)"
	@echo "then run the commands:"
	@echo "gspeed_bed -f 0 -l 645862 -b $(VLARGE_DATADIR)/$(VLARGE_TESTBED) > vlarge.out"
	@echo "and to check result:"
	@echo "python checkDiff.py vlarge.out $(VLARGE_DATADIR)/$(VLARGE_PFILE)"

test_ms_io : gspeed_bed
	for rep in 1 2 3 4 5 ; do \
		bsub -q repair -K -app large -o test_io.out -e test_io.err ./gspeed_bed -f 0 -l 645862 -b $(TESTBED) && sleep 5 && grep 'Read Plink fam' test_io.err >> io_log.log; \
	done

test_ms_mt : gspeed_bed
	for rep in 1 2 3 4 5 ; do \
		bsub -q repair -K -app large -o test_mt.out -e test_mt.err ./gspeed_bed -f 0 -l 645862 -b $(TESTBED) && sleep 5 && grep 'total time' test_mt.err >> mt_log.log; \
	done

test_ms_io_new : gspeed_bed_newio
	bsub -K -app large -o test_io_new.out -e test_io_new.err ./gspeed_bed_newio -f 0 -l 645862 -b $(TESTBED) 

clean-test :
	rm -f small.out mid.out large.out
clean : 
	rm -f *.o *~ ../common/*.o
