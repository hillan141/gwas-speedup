CURRDIR = ../C1
export CXXFLAGS = -I$(CURRDIR) -O2 -fno-operator-names -march=core2 -mssse3 -mfpmath=sse -ffast-math -Wno-deprecated -DMAIN -DHOME_RUN 
CC = $(CXX)
VPATH = ../common
TESTBED = ../t/UC_N6000_M7000_P1_C5/plink
OBJS = GWASSpeedup.o logistic.o ../common/read_bed.o ../common/PlinkIO.o ../common/options.o ../common/nlist.o ../common/helper.o ../common/dcdflib.o

GWASSpeedup.o : GWASSpeedup.cpp GWASSpeedup.h

gspeed.o : gspeed.cpp GWASSpeedup.h

gspeed : gspeed.o GWASSpeedup.o

gspeed_bed : $(OBJS)
	$(CXX) $(CXXFLAGS) -o gspeed_bed $(OBJS)

test : gspeed_bed
		./gspeed_bed -b $(TESTBED)

clean : 
	rm -f *.o *~
