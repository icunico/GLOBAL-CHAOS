MODULE bioptimod_memory

!   PRIVATE
    integer         , parameter :: nw  = 1  ! Number of walength solved  for inversion
    integer         , parameter :: nlt = 33 ! Number of walength considered
    character(1024) :: command_line ! parser argument variable
    integer :: nlev
    integer :: nphy
    double precision :: rd, rs, ru, vs, vu
! downltard planar irradiance [direct and diffuse] just below the sea surface
! provided by 
    integer                       :: wavelength(nlt)
    double precision              :: Ed0m(nw), Es0m(nw), Eu0m(nw)
! Remote Sensing reflectance just above the Sea surface
    double precision              :: Rrs0p(nw)
    double precision, allocatable :: lam(:), lam1(:), lam2(:), aw(:), bw(:), bbw(:), acdom(:), anap(:), bnap(:), bbnap(:) 
    double precision, allocatable :: ac(:,:), ac_ps(:,:), bc(:,:), bbc(:,:)


    CONTAINS

subroutine  Rrs0p_to_Eu0m(inRrs0p, inEd0m, inEs0m, inQ, outEu0m)
 implicit none
  double precision, intent(IN)  :: inRrs0p(nw), inEd0m(nw), inEs0m(nw), inQ  
  double precision, intent(OUT) :: outEu0m(nw)  

! local
  double precision              :: Rrs0m(nw)
! this subroutine applies two corrections 
! 1) derive Rrs0m using the correction by Lee et al. 2002
  Rrs0m(:) = inRrs0p(:)
! 2) correct for the Raman Scattering at surface

  outEu0m(:)  = Rrs0m(:) * (inEd0m(:) + inEs0m(:)) * inQ

end subroutine Rrs0p_to_Eu0m

subroutine compute_total_a_b_bb(nlt, nphy, nlev, chl, C, cdom, nap, a, b, bb)
        implicit none
        integer, intent (IN) :: nlt, nphy, nlev 
        double precision, intent (IN) :: chl(nlev,nphy), C(nlev,nphy), cdom(nlev), nap(nlev)
        double precision, intent(OUT) :: a(nlev,nlt), b(nlev,nlt), bb(nlev,nlt)
        ! local variables
        integer i,wl,p

        do i=1,nlev
           do wl=1,nlt
              a(i,wl)  = aw(wl)  + DOT_PRODUCT(ac(:,wl),chl(i,:))    +acdom(wl)*cdom(i) +  anap(wl)*nap(i)
              write(*,*) 'OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO'
              write(*,*) 'phyto-abs', DOT_PRODUCT(ac(:,wl),chl(i,:))
              write(*,*) 'phyto-ac', ac(i,:)
              write(*,*) 'phyto-chl', chl(i,:)

              b(i,wl)  = bw(wl)  + DOT_PRODUCT(bc(:,wl),  C(i,:))                       +  bnap(wl)*nap(i)
              write(*,*) 'phyto-scatter', DOT_PRODUCT(bc(:,wl),C(i,:))

              bb(i,wl) = bbw(wl) + DOT_PRODUCT(bbc(:,wl)*bc(:,wl), C(i,:))              + bbnap(wl)*nap(i)
              write(*,*) 'phyto-back-scatter', DOT_PRODUCT(bbc(:,wl)*bc(:,wl),C(i,:))
              write(*,*) 'OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO'
           enddo
        enddo
        
end subroutine compute_total_a_b_bb

subroutine  allocate_bio_optical_parameters(nlt, nphy, nlev)
  implicit none
  integer, intent(IN)  :: nlt, nphy, nlev

  allocate(lam(nlt),lam1(nlt),lam2(nlt),aw(nlt),bw(nlt),bbw(nlt))
  allocate(acdom(nlt),anap(nlt),bnap(nlt),bbnap(nlt))
  allocate(ac(nphy,nlt),ac_ps(nphy,nlt),bc(nphy,nlt),bbc(nphy,nlt))

end subroutine allocate_bio_optical_parameters

     subroutine read_command_line
       integer :: exenamelength
       integer :: io, io2

       command_line = ""
       call get_command(command = command_line,status = io)
       if (io==0) then
         call get_command_argument(0,length = exenamelength,status = io2)
         if (io2==0) then
           command_line = "&cmd "//adjustl(trim(command_line(exenamelength+1:)))//" /"
         else
           command_line = "&cmd "//adjustl(trim(command_line))//" /"
         end if
       else
         write(*,*) io,"Error getting command line."
       end if
     end subroutine

     subroutine parse_command_line
       character(256) :: msg
       namelist /cmd/  nlev, nphy
       integer :: io

       if (len_trim(command_line)>0) then
         msg = ''
         read(command_line,nml = cmd,iostat = io,iomsg = msg)
         if (io/=0) then
           write(*,*) "Error parsing the command line or cmd.conf " // msg
           error stop
         end if
       end if
     end subroutine

     subroutine read_coefficients()
     ! namelist example
     ! &coeff
     ! rd=1.0
     ! rs=1.5
     ! ru=3.0
     ! vs=0.83
     ! vu=0.4
     ! /
         implicit none
         namelist /coeff/ rd, rs, ru, vs, vu
         open (unit=20,file="coeff.nml")
         read (20,nml=coeff)
         write (*,*) ' Reading normalization coefficients '
         write (*,nml=coeff)
         write (*,*) ' ++++++++++++++++++++++++++++++++++ '
         close (unit=20)
     end subroutine

     subroutine read_1d_ascii(fname, nlev, var1d)
     ! eg. Read vertical grid mesh
     implicit none
     character(*) :: fname
     integer, intent(in) :: nlev
     double precision, intent(out) :: var1d(nlev)
     ! local variable
     integer :: i
     
     open (unit=15, file=TRIM(fname), status='old',    &
     access='sequential', form='formatted', action='read' )

     do i=1, nlev

         read(15,*) var1d(i)
         write(*,*) var1d(i)
         write(*,*) "________________"

     end do

     close(unit=15)
     end subroutine

     subroutine read_2d_ascii(fname, nlev, nphy, var2d)
     ! eg Read chlorophyll concentration
     implicit none
     character(*) :: fname
     integer, intent(in) ::nlev, nphy
     double precision, intent(out) :: var2d(nlev,nphy)
     ! local variable
     integer :: i,col
     open (unit=16, file=TRIM(fname), status='old',    &
      access='sequential', form='formatted', action='read' )

     do i = 1,nlev
         read(16,*) (var2d(i,col),col=1,nphy)
         write(*,*) (var2d(i,col),col=1,nphy)
         write(*,*) "________________"
     end do

     close(unit=16)
     end subroutine

     subroutine write_2d_ascii(fname, nlev, ncol, var2d)
     ! eg Read chlorophyll concentration
     implicit none
     character(*) :: fname
     integer, intent(in) ::nlev, ncol
     double precision, intent(in) :: var2d(nlev,ncol)
     ! local variable
     integer :: i,col
     open (unit=17, file=TRIM(fname), status='new',    &
      access='sequential', form='formatted', action='write' )

     do i = 1,nlev
         write(17,"(33F8.2)") (var2d(i,col),col=1,ncol)
     end do

     close(unit=17)
     end subroutine



END MODULE bioptimod_memory


