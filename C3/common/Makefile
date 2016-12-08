CURRDIR = ../common
export CPPFLAGS = -I$(CURRDIR) -O3 -fomit-frame-pointer -fno-operator-names -march=core2 -mssse3 -mfpmath=sse -ffast-math -DDEBUG
CC = $(CXX)
# VPATH = ../../../plink

OBJ = options.o PlinkIO.o helper.o nlist.o read_bed.o 

read_bed.o : read_bed.h read_bed.cpp Individual.h Family.h gvar.h WMLocus.h helper.h helper.cpp PlinkIO.h nlist.h nlist.cpp

PlinkIO.o : PlinkIO.cpp

read_bed : $(OBJ) 

bed2txt.o : bed2txt.cpp 

bed2txt : $(OBJ)

test_bed2txt : bed2txt 
	./bed2txt -b ../../../../../t/UC_N6000_M7000_P1_C5/plink > test.bed2txt
	diff test.bed2txt test.bed2txt.ref

clean : 
	rm -f *.o *~
