!************************************************************************************!
!***   The module that introduces the structure to work with the 3-stream light   ***!
!*** TYPE(PAR) is the new type for that. It keeps reflectance, scattering, and    ***!
!*** backscattering, together with cosines and three more parameters rd, rs, ru   ***!
!*** as a single onject with subroutines to work with.                            ***!
!***   These objects can be added, substracted, divided, multiplied, also can be  ***!
!*** multiplied on real numbers.                                                  ***!
!***   Any declared object MUST BE INITIALIZED. This can be done by the function  ***!
!*** par(), which can either take two integers nw and nl (numbers of wavelengths) ***!
!*** and vertical layers, or can clone an existing object. Also, you can call     ***!
!*** a method: call p%par(nw,nl) or call p%init_par(nw,nl)                        ***!
!***   When an object is destroyed, all inner strcutres are freed internally.     ***!
!***   Scalar parameters vd, vs, vu and rd, rs, ru are set by default,            ***!
!*** but can be modified using p%set('vd', value), etc; also, the value can be    ***!
!*** requested by p%get1 (read this as 'get one').                                ***!
!***   2D arrays of parameters a, b, bb, can be set or requested by p%set, p%get  ***!
!*** either in the form p%set('a',values), or by indicating the wavelength number ***!
!*** of depth layer number, or both.                                              ***!
!***   It is possible to set a parameter to a constant value:                     ***!
!*** call p%set_const('a', 1)                                                     ***!
!*** OR reset to zeros: call p%reset                                              ***!
!***   An object can print itself: p%print. Without arguments, this just prints,  ***!
!*** or 3-symbol prefix can be given to be added before e.g. a: 123.              ***!
!***   This is to print stuff like d/da: 123.                                     ***!
!*** It is possible to use a subroutine for printing: call print_par(a+b)         ***!
!************************************************************************************!

MODULE par_mod
implicit none

PUBLIC:: par, operator(+),  operator(-), operator(*), operator(/), print_par,calc_total_parameters

PRIVATE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!These are the constructors. Use them as follows:
!par:: p
!call p%par(nw, nl)
! OR
!p = par(nw,nl)
!OTHER constructors can be added, e.g., to convert something to a parameter object
INTERFACE PAR
    module procedure:: init_par
    module procedure:: return_default_par
    module procedure:: clone_par
END INTERFACE PAR
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!These are overloaded operators
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
INTERFACE OPERATOR(+)              !add two objects, like p1+p2
    module procedure:: add_pars
END INTERFACE OPERATOR(+)

INTERFACE OPERATOR(-)              !substract two objects, like p1-p2
    module procedure:: subs_pars
END INTERFACE OPERATOR(-)

INTERFACE OPERATOR(*)              !multiply an object and a number, like 5+p or p+5. Only for real and dble.
    module procedure:: par_times_x
    module procedure:: par_times_x_real
    module procedure:: x_times_par
    module procedure:: x_times_par_real
    module procedure:: par_times_par
END INTERFACE OPERATOR(*)

INTERFACE OPERATOR(/)             !divide an object on a number. Like p/8. Real or dble. 
    module procedure:: par_div_x
    module procedure:: par_div_x_real
    module procedure:: par_div_par
END INTERFACE OPERATOR(/)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!! The structure !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!This is the structure to contain parameters, together with subs to work with them
TYPE PAR
    integer:: n_wavelen   !number of wavelength 
    integer:: n_layers    !number of layers
    integer:: n_opt       !number of optical constituents 
    double precision, dimension(:,:,:), allocatable:: a, b, bb         !depth-dependent parameters
    double precision, dimension(:,:), allocatable:: Ta, Tb, Tbb        !Total coefficients
    double precision, dimension(:,:), allocatable:: opt_const      ! optical active constituents e.g. water, chlorophyll etc etc 
    double precision, dimension(:,:), allocatable:: a_coeff,b_coeff,bb_coeff ! optical active constituents coefficients e.g. water, chlorophyll etc etc 
    double precision:: vd, vs, vu, rd, rs, ru                        !scalar parameters

    contains                                                         !inner subroutines

    procedure, public:: init_par         !a constructor
    final:: del_par                      !the destructor
    generic, public:: set => set_2d,set_3d, set_for, set_for1, set1, set_2d_real,set_3d_real, set_for_real, set_for1_real, set1_real !different set mehods
    generic, public:: get => get_2d, get_3d, get_for       !several get methods
    procedure, public :: get1                      !separate get-one-value function, to be used as a+p%get1('a',1,2) 
    generic, public :: print => print_me, add_text !print the contents, simply or with a prefix
    procedure, public :: reset                     !put zeros to all parameters
    procedure, public :: set_const                 !set one array to contain constant value
    !private methods, all within generic public ones declared above
    procedure, private :: set_2d
    procedure, private :: set_3d
    procedure, private :: set_for
    procedure, private :: set_for1
    procedure, private :: set1
    procedure, private :: set_2d_real
    procedure, private :: set_3d_real
    procedure, private :: set_for_real
    procedure, private :: set_for1_real
    procedure, private :: set1_real
    procedure, private :: get_2d
    procedure, private :: get_3d
    procedure, private :: get_for
    procedure, private :: print_me
    procedure, private :: add_text
END TYPE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

CONTAINS

!Initilaizer. Must be called prior to using the structure
subroutine init_par(this, nw, nl, nopt)
use bioptimod_memory, ONLY:a_coeff, b_coeff, bb_coeff
    class(par), intent(inout):: this
    integer, intent(in):: nw, nl, nopt
    integer :: j,k
    character(80)    :: header(nopt)
    this%n_wavelen = nw
    this%n_layers  = nl
    this%n_opt     = nopt

    allocate(this%a(nl,nw,nopt))
    allocate(this%b(nl,nw,nopt))
    allocate(this%bb(nl,nw,nopt))

    allocate(this%opt_const(nl,nopt))

    allocate(this%Ta(nl,nw))
    allocate(this%Tb(nl,nw))
    allocate(this%Tbb(nl,nw))

    allocate(this%a_coeff(nw,nopt))
    allocate(this%b_coeff(nw,nopt))
    allocate(this%bb_coeff(nw,nopt))
    !default values
    this%rd = 1.0
    this%rs = 1.5
    this%ru = 3.0
    this%vs = 0.83
    this%vu = 0.4
    this%vd = 0.625 !some value from getmud 
    ! from configuration files
    this%a_coeff=a_coeff
    this%b_coeff=b_coeff
    this%bb_coeff=bb_coeff

    open(11,file='./init_opt_const.txt',status='old')
    read(11,22) (header(j),j=1,nopt)

    do k=1,nl
         read(11,23) (this%opt_const(k,j),j=1,nopt)
    enddo

    close(11)
    write(*,*) "opt_const", this%opt_const

22   FORMAT(A8,A8,A8,A8,A8,A8,A8)
23   FORMAT(f8.5,f8.5,f8.5,f8.5,f8.5,f8.5,f8.5)

end subroutine init_par

!Function as the initializer. Use this as p=par(nw,nl).
function return_default_par(nw, nl, nopt)
    type(par):: return_default_par
    integer, intent(in):: nw, nl, nopt
    call return_default_par%init_par(nw, nl, nopt)
end function return_default_par

!Function as the initializer. Clones an object
function clone_par(p)
    type(par), intent(in):: p
    type(par):: clone_par
    call clone_par%init_par(p%n_wavelen, p%n_layers, p%n_opt)
    clone_par%a = p%a
    clone_par%b = p%b
    clone_par%bb = p%bb
    clone_par%opt_const = p%opt_const
    clone_par%Ta = p%Ta
    clone_par%Tb = p%Tb
    clone_par%Tbb = p%Tbb
    clone_par%vd = p%vd
    clone_par%vs = p%vs
    clone_par%vu = p%vu
    clone_par%rd = p%rd
    clone_par%rs = p%rs
    clone_par%ru = p%ru
end function clone_par

!destructor. Is called automatically when an object is destroyed
subroutine del_par(this)
    type(par), intent(inout):: this
    if (allocated(this%a)) then 
       deallocate(this%a)
    else
       write(*,*) 'Trying to deallocate unallocated this%a'
    endif
    if (allocated(this%b)) then 
       deallocate(this%b)
    else
       write(*,*) 'Trying to deallocate unallocated this%b'
    endif
    if (allocated(this%bb)) then 
       deallocate(this%bb)
    else
       write(*,*) 'Trying to deallocate unallocated this%bb'
    endif
    if (allocated(this%opt_const)) then 
       deallocate(this%opt_const)
    else
       write(*,*) 'Trying to deallocate unallocated this%opt_const'
    endif
        if (allocated(this%Ta)) then
       deallocate(this%Ta)
    else
       write(*,*) 'Trying to deallocate unallocated this%Ta'
    endif
    if (allocated(this%Tb)) then
       deallocate(this%Tb)
    else
       write(*,*) 'Trying to deallocate unallocated this%Tb'
    endif
    if (allocated(this%Tbb)) then
       deallocate(this%Tbb)
    else
       write(*,*) 'Trying to deallocate unallocated this%Tbb'
    endif
end subroutine del_par

!to set a constant value of a parameter
subroutine set_const(this, what, value)
    class(par), intent(inout):: this
    character(len=*), intent(in):: what
    double precision:: value
    select case(trim(what))
    case('a')
        this%a = value
    case('b')
        this%b = value
    case('bb')
        this%bb = value
    case('opt_const')
        this%opt_const = value
    case('Ta')
        this%Ta = value
    case('Tb')
        this%Tb = value
    case('Tbb')
        this%Tbb = value
    case('vd')
        this%vd = value
    case('vs')
        this%vs = value
    case('vu')
        this%vu = value
    case('rd')
        this%rd = value
    case('rs')
        this%rs = value
    case('ru')
        this%ru = value
    case default
        print*, 'Error: wrong "what" parameter ', what, ' in set_const'
    end select
end subroutine set_const

subroutine calc_total_parameters(this)
class(par), intent(inout):: this
!local
integer i,j,k

do i=1,this%n_opt
   do j=1,this%n_wavelen
      do k=1,this%n_layers
         this%a(k,j,i) = this%a_coeff(j,i)*this%opt_const(k,i)
         this%b(k,j,i) = this%b_coeff(j,i)*this%opt_const(k,i)
         this%bb(k,j,i) = this%bb_coeff(j,i)*this%opt_const(k,i)
      enddo
   enddo
enddo

this%Ta=0.0D0
this%Tb=0.0D0
this%Tbb=0.0D0

do i=1,this%n_opt
   do j=1,this%n_wavelen
      do k=1,this%n_layers
         this%Ta(k,j) = this%Ta(k,j) + this%a(k,j,i)
         this%Tb(k,j) = this%Tb(k,j) + this%b(k,j,i)
         this%Tbb(k,j) = this%Tbb(k,j) + this%bb(k,j,i)
      enddo
   enddo
enddo

end subroutine calc_total_parameters

!partial case of the 'set' routine. Sets the 2D wavelength-layers array to one of the parameters a, b, bb, vd
subroutine set_2d(this, what, values)
    class(par), intent(inout):: this
    character(len=*), intent(in):: what
    double precision, dimension(:,:):: values
    integer, dimension(2):: theshape
    theshape = shape(values)
    select case(trim(what))
    case('Ta')
        if(theshape(1).ne.this%n_layers) print*, 'Error: wrong dim 1 of "values'
        if(theshape(2).ne.this%n_wavelen) print*, 'Error: wrong dim 2 of "values'
        this%Ta = values
    case('Tb')
        if(theshape(1).ne.this%n_layers) print*, 'Error: wrong dim 1 of "values'
        if(theshape(2).ne.this%n_wavelen) print*, 'Error: wrong dim 2 of "values'
        this%Tb = values
    case('Tbb')
        if(theshape(1).ne.this%n_layers) print*, 'Error: wrong dim 1 of "values'
        if(theshape(2).ne.this%n_wavelen) print*, 'Error: wrong dim 2 of "values'
        this%Tbb = values
    case('opt_const')
        if(theshape(1).ne.this%n_layers) print*, 'Error: wrong dim 1 of "values'
        if(theshape(2).ne.this%n_opt) print*, 'Error: wrong dim 2 of "values'
        this%Tbb = values
    case default
        print*, 'Error: wrong "what" parameter ', what, ' in set_2d'
    end select
end subroutine set_2d

!partial case of the 'set' routine. Sets the 2D wavelength-layers array to one of the parameters a, b, bb, vd
subroutine set_3d(this, what, values)
    class(par), intent(inout):: this
    character(len=*), intent(in):: what
    double precision, dimension(:,:,:):: values
    integer, dimension(3):: theshape
    theshape = shape(values)
    select case(trim(what))
    case('a')
        if(theshape(1).ne.this%n_layers) print*, 'Error: wrong dim 1 of "values'
        if(theshape(2).ne.this%n_wavelen) print*, 'Error: wrong dim 2 of "values'
        if(theshape(3).ne.this%n_opt) print*, 'Error: wrong dim 3 of "values'
        this%a = values
    case('b')
        if(theshape(1).ne.this%n_layers) print*, 'Error: wrong dim 1 of "values'
        if(theshape(2).ne.this%n_wavelen) print*, 'Error: wrong dim 2 of "values'
        if(theshape(3).ne.this%n_opt) print*, 'Error: wrong dim 3 of "values'
        this%b = values
    case('bb')
        if(theshape(1).ne.this%n_layers) print*, 'Error: wrong dim 1 of "values'
        if(theshape(2).ne.this%n_wavelen) print*, 'Error: wrong dim 2 of "values'
        if(theshape(3).ne.this%n_opt) print*, 'Error: wrong dim 3 of "values'
        this%bb = values
    case default
        print*, 'Error: wrong "what" parameter ', what, ' in set_2d'
    end select
end subroutine set_3d

subroutine set_2d_real(this, what, values)
    class(par), intent(inout):: this
    character(len=*), intent(in):: what
    real, dimension(:,:):: values
  call this%set_2d(what, dble(values))
end subroutine set_2d_real

subroutine set_3d_real(this, what, values)
    class(par), intent(inout):: this
    character(len=*), intent(in):: what
    real, dimension(:,:,:):: values
  call this%set_3d(what, dble(values))
end subroutine set_3d_real


!partial case of the 'set' routine. Sets either the 1D for layers at fixed wavelength, or vice versa, or a single value for wv and layer,
!or one of vs, vu, rd, rs, ru for all wavelength
subroutine set_for(this, what, values, wln, layern, optn)
    class(par), intent(inout):: this
    character(len=*), intent(in):: what
    double precision, dimension(:), intent(in):: values
    integer, intent(in), optional:: wln, layern, optn
    if (present(optn)) then
        if(present(wln).and. .not. present(layern)) then
            select case(trim(what))
            case('a')
                this%a(:,wln,optn) = values
                if(size(values).ne.this%n_layers) print*, 'Error: wrong dim of "values"'
            case('b')
                this%b(:,wln,optn) = values
                if(size(values).ne.this%n_layers) print*, 'Error: wrong dim of "values"'
            case('bb')
                this%bb(:,wln,optn) = values
                if(size(values).ne.this%n_layers) print*, 'Error: wrong dim of "values"'
            case default
                print*, 'Error: wrong "what" parameter ', what, ' in set_for #1'
            end select
        endif
        if( .not. present(wln).and. present(layern)) then
            if(size(values).ne.this%n_wavelen) print*, 'Error: wrong dim of "values"'
            select case(trim(what))
            case('a')
                this%a(layern,:,optn) = values
            case('b')
                this%b(layern,:,optn) = values
            case('bb')
                this%bb(layern,:,optn) = values
            case default
                print*, 'Error: wrong "what" parameter ', what, ' in set_for #2'
            end select
        endif
        if(present(wln).and. present(layern)) then
            select case(trim(what))
            case('a')
                this%a(layern,wln,optn) = values(1)
            case('b')
                this%b(layern,wln,optn) = values(1)
            case('bb')
                this%bb(layern,wln,optn) = values(1)
            case default
                print*, 'Error: wrong "what" parameter ', what, ' in set_for #3'
            end select
        endif
    else
        if(present(wln).and. .not. present(layern)) then
            select case(trim(what))
            case('Ta')
                this%Ta(:,wln) = values
                if(size(values).ne.this%n_layers) print*, 'Error: wrong dim of "values"'
            case('Tb')
                this%Tb(:,wln) = values
                if(size(values).ne.this%n_layers) print*, 'Error: wrong dim of "values"'
            case('Tbb')
                this%Tbb(:,wln) = values
                if(size(values).ne.this%n_layers) print*, 'Error: wrong dim of "values"'
            case default
                print*, 'Error: wrong "what" parameter ', what, ' in set_for #1'
            end select
        endif
        if( .not. present(wln).and. present(layern)) then
            if(size(values).ne.this%n_wavelen) print*, 'Error: wrong dim of "values"'
            select case(trim(what))
            case('Ta')
                this%Ta(layern,:) = values
            case('Tb')
                this%Tb(layern,:) = values
            case('Tbb')
                this%Tbb(layern,:) = values
            case default
                print*, 'Error: wrong "what" parameter ', what, ' in set_for #2'
            end select
        endif
        if(present(wln).and. present(layern)) then
            select case(trim(what))
            case('Ta')
                this%Ta(layern,wln) = values(1)
            case('Tb')
                this%Tb(layern,wln) = values(1)
            case('Tbb')
                this%Tbb(layern,wln) = values(1)
            case default
                print*, 'Error: wrong "what" parameter ', what, ' in set_for #3'
            end select
        endif
        if(.not.(present(wln) .or. present(layern))) then
            if(size(values).ne.this%n_wavelen) print*, 'Error: wrong dim of "values"'
            select case(trim(what))
            case('vd')
                this%vd= values(1)
            case('vs')
                this%vs= values(1)
            case('vu')
                this%vu= values(1)
            case('rs')
                this%rs= values(1)
            case('ru')
                this%ru= values(1)
            case('rd')
                this%rd= values(1)
            case default
                print*, 'Error: wrong "what" parameter ', what, ' in set_for #4'
            end select
        endif
    endif
end subroutine set_for

!partial case of the 'set' routine. Sets either the 1D for layers at fixed wavelength, or vice versa, or a single value for wv and layer,
!or one of vs, vu, rd, rs, ru for all wavelength
subroutine set_for1(this, what, value)
    class(par), intent(inout):: this
    character(len=*), intent(in):: what
    double precision, intent(in):: value
    select case(trim(what))
    case('vd')
        this%vd= value
    case('vs')
        this%vs= value
    case('vu')
        this%vu= value
    case('rs')
        this%rs= value
    case('ru')
        this%ru= value
    case('rd')
        this%rd= value
    case default
        print*, 'Error: wrong "what" parameter ', what, ' in set_for1'
    end select
end subroutine set_for1

!partial case of the 'set' routine. Sets either the 1D for layers at fixed wavelength, or vice versa, or a single value for wv and layer,
!or one of vs, vu, rd, rs, ru for all wavelength
subroutine set_for1_real(this, what, value)
    class(par), intent(inout):: this
    character(len=*), intent(in):: what
    real, intent(in):: value
    select case(trim(what))
    case('vd')
        this%vd= value
    case('vs')
        this%vs= value
    case('vu')
        this%vu= value
    case('rs')
        this%rs= value
    case('ru')
        this%ru= value
    case('rd')
        this%rd= value
    case default
        print*, 'Error: wrong "what" parameter ', what, ' in set_for1_real'
    end select
end subroutine set_for1_real

subroutine set_for_real(this, what, values, wln, layern,optn)
    class(par), intent(inout):: this
    character(len=*), intent(in):: what
    real, dimension(:), intent(in):: values
    integer, intent(in), optional:: wln, layern,optn
  call this%set_for(what, dble(values), wln, layern,optn)
end subroutine set_for_real

!partial case of the 'set' routine. Sets a single value for wv and layer. Does almost the same as set1d.
subroutine set1(this, what, value, wln, layern, optn)
    class(par), intent(inout):: this
    character(len=*), intent(in):: what
    double precision, intent(in):: value
    integer, intent(in):: wln, layern
    integer, intent(in),optional:: optn
        select case(trim(what))
        case('a')
            this%a(layern,wln,optn) = value
        case('b')
            this%b(layern,wln,optn) = value
        case('bb')
            this%bb(layern,wln,optn) = value
        case('Ta')
            this%Ta(layern,wln) = value
        case('Tb')
            this%Tb(layern,wln) = value
        case('Tbb')
            this%Tbb(layern,wln) = value
        case default
            print*, 'Error: wrong "what" parameter ', what, ' in set1'
        end select
end subroutine set1

subroutine set1_real(this, what, value, wln, layern,optn)
    class(par), intent(inout):: this
    character(len=*), intent(in):: what
    real, intent(in):: value
    integer, intent(in):: wln, layern, optn
  call this%set1(what, dble(value), wln, layern,optn)
end subroutine set1_real

!partial case of the 'get' routine. Gets the 2D wavelength-layers array of one of the parameters a, b, bb, vd
subroutine get_2d(this, what, values)
    class(par), intent(inout):: this
    character(len=*), intent(in):: what
    double precision, dimension(:,:), intent(out):: values
    integer, dimension(2):: theshape
    theshape = shape(values)
    if(theshape(1).ne.this%n_layers) print*, 'Error: wrong dim 1 of "values'
    if(theshape(2).ne.this%n_wavelen) print*, 'Error: wrong dim 2 of "values'
    select case(trim(what))
    case('Ta')
        values = this%Ta
    case('Tb')
        values = this%Tb
    case('Tbb')
        values = this%Tbb
    case default
        print*, 'Error: wrong "what" parameter'
    end select
end subroutine get_2d

!partial case of the 'get' routine. Gets the 2D wavelength-layers array of one of the parameters a, b, bb, vd
subroutine get_3d(this, what, values)
    class(par), intent(inout):: this
    character(len=*), intent(in):: what
    double precision, dimension(:,:,:), intent(out):: values
    integer, dimension(3):: theshape
    theshape = shape(values)
    if(theshape(1).ne.this%n_layers) print*, 'Error: wrong dim 1 of "values'
    if(theshape(2).ne.this%n_wavelen) print*, 'Error: wrong dim 2 of "values'
    if(theshape(3).ne.this%n_opt) print*, 'Error: wrong dim 3 of "values'
    select case(trim(what))
    case('a')
        values = this%a
    case('b')
        values = this%b
    case('bb')
        values = this%bb
    case('vd')
        values = this%vd
    case default
        print*, 'Error: wrong "what" parameter'
    end select
end subroutine get_3d

!partial case of the 'get' routine. Gets either the 1D for layers at fixed wavelength, or vice versa, or a single value for wv and layer,
!or one of vs, vu, rd, rs, ru for all wavelength
subroutine get_for(this, what, values, wln, layern, optn)
    class(par), intent(inout):: this
    character(len=*), intent(in):: what
    double precision, dimension(:), intent(out):: values
    integer, intent(in), optional:: wln, layern, optn
    if (present(optn)) then
       if(present(wln).and. .not. present(layern)) then
           select case(trim(what))
           case('a')
               values = this%a(:,wln,optn)
               if(size(values).ne.this%n_layers) print*, 'Error: wrong dim of "values"'
           case('b')
               values = this%b(:,wln,optn)
               if(size(values).ne.this%n_layers) print*, 'Error: wrong dim of "values"'
           case('bb')
               values = this%bb(:,wln,optn)
               if(size(values).ne.this%n_layers) print*, 'Error: wrong dim of "values"'
           case default
               print*, 'Error: wrong "what" parameter'
           end select
       endif
       if(.not. present(wln).and. present(layern)) then
           if(size(values).ne.this%n_wavelen) print*, 'Error: wrong dim of "values"'
           select case(trim(what))
           case('a')
               values = this%a(layern,:,optn)
           case('b')
               values = this%b(layern,:,optn)
           case('bb')
               values = this%bb(layern,:,optn)
           case default
               print*, 'Error: wrong "what" parameter'
           end select
       endif
       if(present(wln).and. present(layern)) then
           select case(trim(what))
           case('a')
               values(1) = this%a(layern,wln,optn)
           case('b')
               values(1) = this%b(layern,wln,optn)
           case('bb')
               values(1) = this%bb(layern,wln,optn)
           case default
               print*, 'Error: wrong "what" parameter'
           end select
       endif
    else
       if(present(wln).and. .not. present(layern)) then
           select case(trim(what))
           case('Ta')
               values = this%Ta(:,wln)
               if(size(values).ne.this%n_layers) print*, 'Error: wrong dim of "values"'
           case('Tb')
               values = this%Tb(:,wln)
               if(size(values).ne.this%n_layers) print*, 'Error: wrong dim of "values"'
           case('Tbb')
               values = this%Tbb(:,wln)
               if(size(values).ne.this%n_layers) print*, 'Error: wrong dim of "values"'
           case default
               print*, 'Error: wrong "what" parameter'
           end select
       endif
       if(.not. present(wln).and. present(layern)) then
           if(size(values).ne.this%n_wavelen) print*, 'Error: wrong dim of "values"'
           select case(trim(what))
           case('Ta')
               values = this%Ta(layern,:)
           case('Tb')
               values = this%Tb(layern,:)
           case('Tbb')
               values = this%Tbb(layern,:)
           case default
               print*, 'Error: wrong "what" parameter'
           end select
       endif
       if(present(wln).and. present(layern)) then
           select case(trim(what))
           case('Ta')
               values(1) = this%Ta(layern,wln)
           case('Tb')
               values(1) = this%Tb(layern,wln)
           case('Tbb')
               values(1) = this%Tbb(layern,wln)
           case default
               print*, 'Error: wrong "what" parameter'
           end select
       endif
       if(.not.(present(wln) .or. present(layern))) then
           select case(trim(what))
           case('vd')
               values(1) = this%vd
           case('vs')
               values(1) = this%vs
           case('vu')
               values(1) = this%vu
           case('rs')
               values(1) = this%rs
           case('ru')
               values(1) = this%ru
           case('rd')
               values(1) = this%rd
           case default
               print*, 'Error: wrong "what" parameter'
           end select
       endif
    endif
end subroutine get_for

!a function to get a single chosen value for a given wavlength number and layer number
function get1(this, what, wln, layern,optn)
    class(par), intent(inout):: this
    character(len=*), intent(in):: what
    integer, intent(in), optional:: wln
    integer, intent(in), optional:: layern
    integer, intent(in), optional:: optn
    double precision:: get1
    if(present(wln).and.present(layern) .and. present(optn) ) then !return a wavelength- and depth-dependent value
        select case(trim(what))
        case('a')
            get1 = this%a(layern,wln,optn)
        case('b')
            get1 = this%b(layern,wln,optn)
        case('bb')
            get1 = this%bb(layern,wln,optn)
        case('vd')
            get1 = this%vd
        case('vs')
            get1 = this%vs
        case('vu')
            get1 = this%vu
        case('rd')
            get1 = this%rd
        case('rs')
            get1 = this%rs
        case('ru')
            get1 = this%ru
        case default
            print*, 'Error: wrong "what" parameter ', what, 'get1 #2'
        end select
    elseif (present(wln).and.present(layern) ) then
        select case(trim(what))
        case('Ta')
            get1 = this%Ta(layern,wln)
        case('Tb')
            get1 = this%Tb(layern,wln)
        case('Tbb')
            get1 = this%Tbb(layern,wln)
        case('vd')
            get1 = this%vd
        case('vs')
            get1 = this%vs
        case('vu')
            get1 = this%vu
        case('rd')
            get1 = this%rd
        case('rs')
            get1 = this%rs
        case('ru')
            get1 = this%ru
        case default
            print*, 'Error: wrong "what" parameter ', what, 'get1 #2'
        end select
    else !if at least one is absent, then a constant parameter is required
        select case(trim(what))
        case('vd')
            get1 = this%vd
        case('vs')
            get1 = this%vs
        case('vu')
            get1 = this%vu
        case('rd')
            get1 = this%rd
        case('rs')
            get1 = this%rs
        case('ru')
            get1 = this%ru
        case default
            print*, 'Error: wrong "what" parameter ', what, 'get1 #3'
        end select
    endif
end function get1

function add_pars(p1,p2)
 type(par), intent(in):: p1, p2
 type(par):: add_pars
    add_pars=par(p1)
    add_pars%a  = p1%a  + p2%a 
    add_pars%b  = p1%b  + p2%b 
    add_pars%bb = p1%bb + p2%bb
    add_pars%opt_const = p1%opt_const + p2%opt_const
    add_pars%Ta  = p1%Ta  + p2%Ta 
    add_pars%Tb  = p1%Tb  + p2%Tb 
    add_pars%Tbb = p1%Tbb + p2%Tbb
end function add_pars

function subs_pars(p1,p2)
 type(par), intent(in):: p1, p2
 type(par):: subs_pars
    subs_pars=par(p1)
    subs_pars%a  = p1%a  - p2%a 
    subs_pars%b  = p1%b  - p2%b 
    subs_pars%bb = p1%bb - p2%bb
    subs_pars%opt_const = p1%opt_const - p2%opt_const
    subs_pars%Ta  = p1%Ta  - p2%Ta 
    subs_pars%Tb  = p1%Tb  - p2%Tb 
    subs_pars%Tbb = p1%Tbb - p2%Tbb
end function subs_pars

function par_times_x(p,x)
 type(par), intent(in):: p
 double precision, intent(in):: x
 type(par):: par_times_x 
    par_times_x = par(p)
    par_times_x%a  = p%a  * x
    par_times_x%b  = p%b  * x
    par_times_x%bb = p%bb * x
    par_times_x%opt_const = p%opt_const * x
    par_times_x%Ta  = p%Ta  * x
    par_times_x%Tb  = p%Tb  * x
    par_times_x%Tbb = p%Tbb * x
end function par_times_x

function x_times_par(x,p)
 double precision, intent(in):: x
 type(par), intent(in):: p
 type(par):: x_times_par 
 x_times_par = par_times_x(p,x)
end function x_times_par

function par_times_x_real(p,x)
 type(par), intent(in):: p
 real, intent(in):: x
 type(par):: par_times_x_real 
    par_times_x_real = p * dble(x)
end function par_times_x_real

function x_times_par_real(x,p)
 real, intent(in):: x
 type(par), intent(in):: p
 type(par):: x_times_par_real 
 x_times_par_real = par_times_x(p,dble(x))
end function x_times_par_real

function par_times_par(p1,p2)
 type(par), intent(in):: p1,p2
 type(par):: par_times_par 
    par_times_par = par(p1)
    par_times_par%a  = p1%a  * p2%a 
    par_times_par%b  = p1%b  * p2%b 
    par_times_par%bb = p1%bb * p2%bb
    par_times_par%opt_const = p1%opt_const * p2%opt_const
    par_times_par%Ta  = p1%Ta  * p2%Ta 
    par_times_par%Tb  = p1%Tb  * p2%Tb 
    par_times_par%Tbb = p1%Tbb * p2%Tbb
end function par_times_par

function par_div_x(p,x)
 type(par), intent(in):: p
 double precision, intent(in):: x
 type(par):: par_div_x 
    par_div_x = par(p)
    par_div_x%a  = p%a  / x
    par_div_x%b  = p%b  / x
    par_div_x%bb = p%bb / x
    par_div_x%opt_const  = p%opt_const  / x
    par_div_x%Ta  = p%Ta  / x
    par_div_x%Tb  = p%Tb  / x
    par_div_x%Tbb = p%Tbb / x
end function par_div_x

function par_div_x_real(p,x)
 type(par), intent(in):: p
 real, intent(in):: x
 type(par):: par_div_x_real 
    par_div_x_real = p / dble(x)
end function par_div_x_real

function par_div_par(p1,p2)
 type(par), intent(in):: p1,p2
 type(par):: par_div_par 
    par_div_par = par(p1)
    par_div_par%a  = p1%a  / p2%a
    par_div_par%b  = p1%b  / p2%b
    par_div_par%bb = p1%bb / p2%bb
    par_div_par%opt_const  = p1%opt_const  / p2%opt_const
    par_div_par%Ta  = p1%Ta  / p2%Ta
    par_div_par%Tb  = p1%Tb  / p2%Tb
    par_div_par%Tbb = p1%Tbb / p2%Tbb
end function par_div_par

subroutine print_me(this)
    class(par), intent(inout):: this
    integer:: wv,op
    do op=1,this%n_opt
        print*, 'Opt const #', op, ':'
        do wv=1,this%n_wavelen
            print*, 'Wavelength #', wv, ':'
            print 123, 'a ', this%a(:,wv,op)
            print 123, 'b ', this%b(:,wv,op)
            print 123, 'bb', this%bb(:,wv,op)
        enddo
    enddo
    do wv=1,this%n_wavelen
        print*, 'Wavelength #', wv, ':'
        print 123, 'Ta ', this%Ta(:,wv)
        print 123, 'Tb ', this%Tb(:,wv)
        print 123, 'Tbb', this%Tbb(:,wv)
    enddo
    print 124, 'vd,vs,vu', this%vd, this%vs, this%vu
    print 124, 'rd,rs,ru', this%rd, this%rs, this%ru
123 format (2X,A2,":",1X,666(F9.6,1X))       
124 format (2X,A8,":",1X,666(F9.6,1X))       
end subroutine print_me

subroutine add_text(this, text)
    class(par), intent(inout):: this
    character(len=3), intent(in):: text
    integer:: wv,op
    do op=1,this%n_opt
        print*, 'Opt const #', op, ':'
        do wv=1,this%n_wavelen
            print*, 'Wavelength #', wv, ':'
            print 125, text//'a ', this%a(:,wv,op)
            print 125, text//'b ', this%b(:,wv,op)
            print 125, text//'bb', this%bb(:,wv,op)
        enddo
    enddo
    do wv=1,this%n_wavelen
        print*, 'Wavelength #', wv, ':'
        print 125, text//'Ta ', this%Ta(:,wv)
        print 125, text//'Tb ', this%Tb(:,wv)
        print 125, text//'Tbb', this%Tbb(:,wv)
    enddo
    print 126, text//'vd,vs,vu', this%vd, this%vs, this%vu
    print 126, text//'rd,rs,ru', this%rd, this%rs, this%ru
125 format (2X,A5,":",1X,666(F9.6,1X))       
126 format (2X,A11,":",1X,666(F9.6,1X))       
end subroutine add_text

subroutine print_par(p, text)
type(par), intent(in):: p
character(len=3), optional:: text
type(par):: tmp
tmp=par(p)
if(present(text)) then
    call tmp%add_text(text)
else
    call tmp%print
endif
end subroutine print_par

subroutine reset(this)
    class(par), intent(inout):: this
    this%a=0.0d0
    this%b=0.0d0
    this%bb=0.0d0
    this%opt_const=0.0d0
    this%Ta=0.0d0
    this%Tb=0.0d0
    this%Tbb=0.0d0
end subroutine reset

END MODULE Par_mod

SUBROUTINE TEST_par
use par_mod
type(par):: p, p2
p=par(1,1,1)
call p%set('a',[1.0],wln=1)
call p%set('b',[2.0],wln=1)
call p%set('bb',[3.0],wln=1)
call p%set('vd',4.0d0)
call p%set('vs',5.0)
call p%set('vu',[6.0])
call p%set('rd',[7.0])
call p%set('rs',[8.0])
call p%set('ru',[9.0])
print*, 'initial p:'
call p%print
p2=par(p)
p2=p2*5.0d0
print*, 'p*5:'
call p2%print
p2=1./5.*p2;
print*, 'p2=p/5:'
call p2%print
p2=p+0.7*p2
print*, 'p+0.7p2:'
call p2%print
p2 = p-p
print*, 'null:'
call p2%print
p2 = p/8.
print*, 'p/8:'
call p2%print
print*, '================'
print*, 'Here is p:'
call p%print
print*, '///'
call p2%print
print*, 'test ='
p2=p
print*, 'p2==p?'
call p2%print
p2=p2*1.1
print*, 'p2\=p?'
call p2%print
call p%print
print*, '=========='
print*, 'a=',p%get1('a',1,1)
print*, 'vd=',p%get1('vd')
print*, 'RESET:'
call p2%reset
call p2%print
call p%set_const('a', 42.0d0)
call p%print('>  ')
print*, '>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
call print_par(p/p)
call print_par((p+0.2*p)/p)
END SUBROUTINE TEST_par

!PROGRAM TEST
!call TEST_par
!END PROGRAM TEST
