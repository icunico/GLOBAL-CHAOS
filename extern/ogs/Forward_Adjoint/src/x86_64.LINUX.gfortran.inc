FC      = gfortran $(OPENMP_FLAG)
CC      = gcc
F90     = gfortran $(OPENMP_FLAG)
F77     = gfortran
CPP     = gfortran -cpp -E
LD      = $(FC) $(OPENMP_FLAG)

