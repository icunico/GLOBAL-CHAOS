source $HOME/sequence.sh

ifort -no-wrap-margin bioptimod_memory.f90 fm34.f90 SLAE.F90 adj_new.f90 main_adj.f90 -o adj.xx
