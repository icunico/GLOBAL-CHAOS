program  bio_optical_adjoint
use bioptimod_memory 
use adj_3stream, only: compute_3stream_adjoint
IMPLICIT NONE
!local
integer           :: i
CHARACTER(len=32) :: arg
double precision  :: Q

! Init

Q = 4.0D0
sunz = 0.0D0

CALL getarg(1, arg)
READ(arg,fmt=*) nw

call allocate_mem()

    open (unit=15, file="surfdata.txt", status='old',    &
          access='sequential', form='formatted', action='read' )

    do i=1, nw

       read(15,*) wavelength(i), Rrs0p_sat(i), Ed0mOASIM(i), Es0mOASIM(i), sunz

       write(*,*) 'Wavelenght ', wavelength(i)
       write(*,*) 'Rrs',  Rrs0p_sat(i)
       write(*,*) 'Ed0m', Ed0mOASIM(i)
       write(*,*) 'Es0m', Es0mOASIM(i)
       write(*,*) 'sunz', sunz
       write(*,*) "________________"

    end do

    close(unit=15)

call init_coeff()


Ed0m(:) = Ed0mOASIM(:)
Es0m(:) = Es0mOASIM(:)

call Rrs0p_to_Eu0m(Q)

!Eu0m(:) = Eu0m_sat(:)

! compute
 call compute_3stream_adjoint()

!finalize
end program bio_optical_adjoint

