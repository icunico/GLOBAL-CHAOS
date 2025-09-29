program  bio_optical_adjoint
use bioptimod_memory,  only: nw, wavelength, Rrs0p_to_Eu0m, Ed0m, Es0m, Eu0m
use adj_3stream, only: compute_3stream_adjoint
!local
integer :: i
double precision :: Q
double precision :: Ed0mOASIM(nw), Es0mOASIM(nw)
double precision :: Rrs0p_sat(nw), Eu0m_sat(nw)

! Init

Q = 4.0D0

open (unit=15, file="surfdata.txt", status='old',    &
      access='sequential', form='formatted', action='read' )

do i=1, nw

   read(15,*) wavelength(i), Rrs0p_sat(i), Ed0mOASIM(i), Es0mOASIM(i)

   write(*,*) wavelength(i)
   write(*,*) Rrs0p_sat(i)
   write(*,*) Ed0mOASIM(i)
   write(*,*) Es0mOASIM(i)
   write(*,*) "________________"

end do

close(unit=15)

Ed0m(:) = Ed0mOASIM(:)
Es0m(:) = Es0mOASIM(:)

call Rrs0p_to_Eu0m(Rrs0p_sat, Ed0mOASIM, Es0mOASIM, Q, Eu0m_sat)

write(*,*) "Eu0m_sat", Eu0m_sat

Eu0m(:) = Eu0m_sat(:)

! compute
call compute_3stream_adjoint()

!finalize
end program bio_optical_adjoint

