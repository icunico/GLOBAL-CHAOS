include compiler.inc
OPT = #-no-wrap-margin

EXE = adj.xx

SOURCE = bioptimod_memory.o fm34.o SLAE.o adj_new.o main_adj.o
SOURCE_LIB = bioptimod_memory.o fm34.o SLAE.o adj_new.o 

all: $(SOURCE)
	$(FC) $(OPT) $(SOURCE) -o $(EXE)

%.o: %.f90
	$(FC) $(OPT) -c $< -o $*.o

%.o: %.F90
	$(FC) $(OPT) -c $< -o $*.o

libadj.a: $(SOURCE_LIB)
	ar rcs libadj.a $(SOURCE_LIB)
	cp libadj.a ../lib/
	cp *.mod ../include
clear:
	rm *.o
	rm ../include/*.mod
	rm ../lib/*
	rm libadj.a 

clean: clear

SAVE = *.f90 *.txt *.nml Makefile

tar: 
	tar -zcf backup_`date +%F`.tar.gz $(SAVE)
	echo 'backup_'`date +%F`'.tar.gz created'
lib: libadj.a

