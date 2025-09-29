FC      = mpiifort $(OPENMP_FLAG)
CC      = mpiicc
F90     = mpiifort $(OPENMP_FLAG)
F77     = mpiifort
CPP     = mpiifort -cpp -E
LD      = $(FC) $(OPENMP_FLAG)

