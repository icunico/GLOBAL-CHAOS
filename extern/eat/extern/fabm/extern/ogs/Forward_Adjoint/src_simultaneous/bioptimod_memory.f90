MODULE bioptimod_memory

!   PUBLIC:: test_3stream_adjoint, get_derivs

!   PRIVATE
    integer                     :: nw      ! Number of walength solved 
    integer         , parameter :: no  = 7 ! Number of bioptical constituents
! downward planar irradiance [direct and diffuse] just below the sea surface
! provided by 
    double precision,allocatable :: wavelength(:)
    double precision,allocatable :: Ed0m(:), Es0m(:), Eu0m(:)
    double precision,allocatable :: Ed0mOASIM(:), Es0mOASIM(:)
    double precision,allocatable :: Rrs0p_sat(:), Eu0m_sat(:)
! Remote Sensing reflectance just above the Sea surface
    double precision,allocatable :: Rrs0p(:)
    double precision             :: sunz ! solar zenith angle
    double precision             :: a_init, b_init, bb_init 
    double precision,parameter   :: rd_init = 1.0D0
    double precision,parameter   :: rs_init = 1.5D0
    double precision,parameter   :: ru_init = 3.0D0
    double precision,parameter   :: vs_init = 0.83D0
    double precision,parameter   :: vu_init = 0.4D0
    double precision,allocatable :: a_coeff(:,:), b_coeff(:,:), bb_coeff(:,:)
    character(80)                :: opt_const_name(7)

    CONTAINS

subroutine allocate_mem()
   IMPLICIT NONE

   allocate(wavelength(nw), Ed0m(nw), Es0m(nw), Eu0m(nw), Rrs0p(nw))
   allocate(Ed0mOASIM(nw), Es0mOASIM(nw),Rrs0p_sat(nw), Eu0m_sat(nw))
   write(*,*) "nw = ", nw
   write(*,*) 'inRrs0p',size(Rrs0p)

end subroutine allocate_mem

subroutine  Rrs0p_to_Eu0m(Q)
 implicit none
  double precision              :: Q
  double precision, parameter   :: T=0.52D0,GammaQ=1.7D0
! local
  double precision              :: Rrs0m(nw)
! this subroutine applies two corrections 
! 1) derive Rrs0m using the correction by Lee et al. 2002
!write(*,*) 'Rrs0m',size(Rrs0m)
  Rrs0m(:) = Rrs0p_sat(:)/(T+GammaQ*Rrs0p_sat(:))
! 2) correct for the Raman Scattering at surface

  Eu0m(:)  = Rrs0m(:) * (Ed0m(:) + Es0m(:)) * Q

end subroutine Rrs0p_to_Eu0m

subroutine init_coeff()
implicit none 
!local
  integer          :: i,ii,j
  integer,parameter :: n_opt=7
  character(80)    :: header(n_opt+1)
  double precision :: a(9,n_opt+1), b(9,n_opt+1), bb(9,n_opt+1)

  allocate(a_coeff(nw,no), b_coeff(nw,no), bb_coeff(nw,no))


! absorption
  open(11,file='./abs_coeff.txt',status='old')


  read(11,22) (header(j),j=1,n_opt+1)


  do j=1, n_opt
     opt_const_name(j)=header(j+1)
  enddo

  write(*,*) 'Header ', (TRIM( header(j) ),j=1,n_opt+1)
  write(*,*)(TRIM( opt_const_name(j) ),j=1,n_opt)

  do i = 1, 9 
    read(11,23) (a(i,j),j=1,n_opt+1)
    write(*,*) a(i,:)
  end do
  close(11) 

    do ii =1,nw 
       do i = 1, 9 
           if(wavelength(ii) .EQ. a(i,1)) then
               do j=1,no
                   a_coeff(ii,j) = a(i,j+1)
               enddo
           endif
       enddo
    enddo
!   scattering
  open(11,file='./scat_coeff.txt',status='old')

  read(11,22) (header(j),j=1,n_opt+1)
  write(*,*) 'Header ', (TRIM( header(j) ),j=1,n_opt+1)
! write(*,*) header

  do i = 1, 9
    read(11,23) (b(i,j),j=1,n_opt+1)
  end do
    do ii =1,nw
       do i = 1, 9
       if(wavelength(ii) .EQ. b(i,1)) then
           do j=1,no
               b_coeff(ii,j) = b(i,j+1)
           enddo
       endif
       enddo
    enddo

  close(11) 

! back scattering
  open(11,file='./back_scat_coeff.txt',status='old')
  read(11,22) (header(j),j=1,n_opt+1)
  write(*,*) 'Header ', (TRIM( header(j) ),j=1,n_opt+1)

  do i = 1, 9
    read(11,23) (bb(i,j),j=1,n_opt+1)
  end do
    do ii =1,nw
       do i = 1, 9
       if(wavelength(ii) .EQ. bb(i,1)) then
           do j=1,no
               bb_coeff(ii,j) = bb(i,j+1)
           enddo
       endif
       enddo
    enddo

  close(11)

22    FORMAT(A6,A8,A8,A8,A8,A8,A8,A8)
23    FORMAT(f6.2,f8.5,f8.5,f8.5,f8.5,f8.5,f8.5,f8.5)

  write(*,*) 'Absorption parameters'
  do i =1,nw
     write(*,*) wavelength(i), (a_coeff(i,j),j=1,no)
  enddo
  write(*,*) 'Scattering  parameters'
  do i =1,nw
     write(*,*) wavelength(i), (b_coeff(i,j),j=1,no)
  enddo
  write(*,*) 'Back scattering  parameters'
  do i =1,nw
     write(*,*) wavelength(i), (bb_coeff(i,j),j=1,no)
  enddo
end subroutine init_coeff

END MODULE bioptimod_memory


