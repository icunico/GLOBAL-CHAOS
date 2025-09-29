module ogs_bfm_shared

   use fabm_types
   use fabm_standard_variables

   implicit none

   public

   integer,parameter :: LightPeriodFlag=1

   real(rk),parameter :: CMass         = 12.011_rk
   real(rk),parameter :: ONE           = 1._rk
   real(rk),parameter :: ZERO          = 0._rk
   real(rk),parameter :: BASETEMP      = 10._rk
   real(rk),parameter :: p_small       = 1.0E-20_rk
   real(rk),parameter :: qnRPIcX       = 1.26E-02_rk
   real(rk),parameter :: qpRPIcX       = 7.86E-04_rk
   real(rk),parameter :: qsRPIcX       = 15._rk/106._rk/CMass
   real(rk),parameter :: ZeroX         = 1e-8_rk
   real(rk),parameter :: pi            = acos(-1._rk)
   real(rk),parameter :: deg2rad       = pi/180._rk
   real(rk),parameter :: WtoQuanta     = 4.57_rk
   real(rk),parameter :: SEC_PER_DAY   = 86400.0_rk
   real(rk),parameter :: SUNQ          = 24.0_rk
   real(rk),parameter :: HOURS_PER_DAY = 24.0_rk
   real(rk),parameter :: MW_C          = 12.0_rk
   real(rk),parameter :: C2ALK         = 2.0_rk/MW_C   ! Conversion factor between inorganic carbon and alkalinity
   real(rk),parameter :: p_atm0         = 1013.25_rk    !reference sea level pressure
   real(rk),parameter :: ZERO_KELVIN   = -273.15_rk;
   real(rk),parameter :: h_planck      = 6.6256E-34   !Plancks constant J sec
   real(rk),parameter :: c_light       = 2.998E8      !speed of light m/sec
   real(rk),parameter :: oavo          = 1.0D0/6.023E23   ! 1/Avogadros number

   real(rk)           :: flPTN6r    ! total rate of formation of reduction equivalent [mmolHS/m3/d] computed in PelBac and used in PelChem   
   real(rk)           :: qccPPY     ! PIC:POC ration in P2: compputed in Phyto and used in MicroZoo and MesoZoo only for prey P2 (crapy solution)

#ifdef IRON
   logical,parameter :: use_iron = .true.
#else
   logical,parameter :: use_iron = .false.
#endif

   ! Aggregate diagnostics for e.g., carbon budgets.
   type (type_bulk_standard_variable),parameter :: total_chlorophyll = type_bulk_standard_variable(name='total_chlorophyll',units='mg/m^3',aggregate_variable=.true.)
   type (type_bulk_standard_variable),parameter :: total_oxygen = type_bulk_standard_variable(name='total_oxygen',units='mmolO2/m^3',aggregate_variable=.true.)
   type (type_bulk_standard_variable),parameter :: total_reduction_equivalent = type_bulk_standard_variable(name='total_reduction_equivalent',units='mmolEq/m^3',aggregate_variable=.true.)
   type (type_bulk_standard_variable),parameter :: alkalinity = type_bulk_standard_variable(name='alkalinity',units='mmolEq/m^3',aggregate_variable=.true.)
   type (type_bulk_standard_variable),parameter :: total_calcite_in_biota = type_bulk_standard_variable(name='total_calcite_in_biota',units='mg C/m^3',aggregate_variable=.true.)
   type (type_bulk_standard_variable),parameter :: total_silicate = type_bulk_standard_variable(name='total_silicate',units='mmolSi/m^3',aggregate_variable=.true.)
   type (type_bulk_standard_variable),parameter :: secchi_depth = type_bulk_standard_variable(name='secchi_depth',units='m')

   ! Aggregate variables for benthic bioturbation and bioirrigation (summed over all fauna).
   type (type_bulk_standard_variable),parameter :: total_bioturbation_activity = type_bulk_standard_variable(name='total_bioturbation_activity',units='mg C/m^2/d',aggregate_variable=.true.)
   type (type_bulk_standard_variable),parameter :: total_bioirrigation_activity = type_bulk_standard_variable(name='total_bioirrigation_activity',units='mg C/m^2/d',aggregate_variable=.true.)

   ! Spectral light variables
   real(rk), allocatable, dimension(:)                 :: lam,lam1,lam2,aw,bw,bbw,apoc,bpoc,bbpoc,WtoQ,acdom_min
   real(rk), allocatable, dimension(:)                 :: Ed_0,Es_0
   real(rk), allocatable, dimension(:,:)               :: ac,ac_ps,bc,bbc,acdom
!  real(rk), allocatable, dimension(:,:)               :: a_array, b_array, bb_array
!  type (type_surface_standard_variable),parameter     :: surf_direct_downward_irradiance_250_nm = type_surface_standard_variable(name='surf_direct_downward_irradiance_250_nm',units='W/m2')
!  type (type_surface_standard_variable),parameter     :: surf_diffuse_downward_irradiance_250_nm = type_surface_standard_variable(name='surf_diffuse_downward_irradiance_250_nm',units='W/m2')
!  type (type_surface_standard_variable),parameter     :: surf_diffuse_direct_irradiance_325_nm = type_surface_standard_variable(name='surf_direct_downward_irradiance_325_nm',units='W/m2')
!  type (type_surface_standard_variable),parameter     :: surf_diffuse_downward_irradiance_325_nm = type_surface_standard_variable(name='surf_diffuse_downward_irradiance_325_nm',units='W/m2')
    type (type_bulk_standard_variable),parameter       :: PAR_tot = type_bulk_standard_variable(name='PAR_tot',units='<UNITS>',aggregate_variable=.true.)
    type (type_bulk_standard_variable),parameter       :: PAR_P1 = type_bulk_standard_variable(name='PAR_P1',units='<UNITS>',aggregate_variable=.true.)
    type (type_bulk_standard_variable),parameter       :: PAR_P2 = type_bulk_standard_variable(name='PAR_P2',units='<UNITS>',aggregate_variable=.true.)
    type (type_bulk_standard_variable),parameter       :: PAR_P3 = type_bulk_standard_variable(name='PAR_P3',units='<UNITS>',aggregate_variable=.true.)
    type (type_bulk_standard_variable),parameter       :: PAR_P4 = type_bulk_standard_variable(name='PAR_P4',units='<UNITS>',aggregate_variable=.true.)
    type (type_bulk_standard_variable),parameter       :: PAR_P5 = type_bulk_standard_variable(name='PAR_P5',units='<UNITS>',aggregate_variable=.true.)
    type (type_bulk_standard_variable),parameter       :: PAR_P6 = type_bulk_standard_variable(name='PAR_P6',units='<UNITS>',aggregate_variable=.true.)
    type (type_bulk_standard_variable),parameter       :: PAR_P7 = type_bulk_standard_variable(name='PAR_P7',units='<UNITS>',aggregate_variable=.true.)
    type (type_bulk_standard_variable),parameter       :: PAR_P8 = type_bulk_standard_variable(name='PAR_P8',units='<UNITS>',aggregate_variable=.true.)
    type (type_bulk_standard_variable),parameter       :: PAR_P9 = type_bulk_standard_variable(name='PAR_P9',units='<UNITS>',aggregate_variable=.true.)

   ! Aggregate chlorophyll per optical type:
   type (type_bulk_standard_variable),parameter :: chlorophyll_P1 = type_bulk_standard_variable(name='chlorophyll_P1',units='mg m-3',aggregate_variable=.true.)
   type (type_bulk_standard_variable),parameter :: chlorophyll_P2 = type_bulk_standard_variable(name='chlorophyll_P2',units='mg m-3',aggregate_variable=.true.)
   type (type_bulk_standard_variable),parameter :: chlorophyll_P3 = type_bulk_standard_variable(name='chlorophyll_P3',units='mg m-3',aggregate_variable=.true.)
   type (type_bulk_standard_variable),parameter :: chlorophyll_P4 = type_bulk_standard_variable(name='chlorophyll_P4',units='mg m-3',aggregate_variable=.true.)
   type (type_bulk_standard_variable),parameter :: chlorophyll_P5 = type_bulk_standard_variable(name='chlorophyll_P5',units='mg m-3',aggregate_variable=.true.)
   type (type_bulk_standard_variable),parameter :: chlorophyll_P6 = type_bulk_standard_variable(name='chlorophyll_P6',units='mg m-3',aggregate_variable=.true.)
   type (type_bulk_standard_variable),parameter :: chlorophyll_P7 = type_bulk_standard_variable(name='chlorophyll_P7',units='mg m-3',aggregate_variable=.true.)
   type (type_bulk_standard_variable),parameter :: chlorophyll_P8 = type_bulk_standard_variable(name='chlorophyll_P8',units='mg m-3',aggregate_variable=.true.)
   type (type_bulk_standard_variable),parameter :: chlorophyll_P9 = type_bulk_standard_variable(name='chlorophyll_P9',units='mg m-3',aggregate_variable=.true.)

   ! Aggregate carbon per optical type:
   type (type_bulk_standard_variable),parameter :: carbon_P1 = type_bulk_standard_variable(name='carbon_P1',units='mg m-3',aggregate_variable=.true.)
   type (type_bulk_standard_variable),parameter :: carbon_P2 = type_bulk_standard_variable(name='carbon_P2',units='mg m-3',aggregate_variable=.true.)
   type (type_bulk_standard_variable),parameter :: carbon_P3 = type_bulk_standard_variable(name='carbon_P3',units='mg m-3',aggregate_variable=.true.)
   type (type_bulk_standard_variable),parameter :: carbon_P4 = type_bulk_standard_variable(name='carbon_P4',units='mg m-3',aggregate_variable=.true.)
   type (type_bulk_standard_variable),parameter :: carbon_P5 = type_bulk_standard_variable(name='carbon_P5',units='mg m-3',aggregate_variable=.true.)
   type (type_bulk_standard_variable),parameter :: carbon_P6 = type_bulk_standard_variable(name='carbon_P6',units='mg m-3',aggregate_variable=.true.)
   type (type_bulk_standard_variable),parameter :: carbon_P7 = type_bulk_standard_variable(name='carbon_P7',units='mg m-3',aggregate_variable=.true.)
   type (type_bulk_standard_variable),parameter :: carbon_P8 = type_bulk_standard_variable(name='carbon_P8',units='mg m-3',aggregate_variable=.true.)
   type (type_bulk_standard_variable),parameter :: carbon_P9 = type_bulk_standard_variable(name='carbon_P9',units='mg m-3',aggregate_variable=.true.)
   
   ! Standard benthic variables used to make implicit based on matching standard names coupling possible.
   type (type_horizontal_standard_variable),parameter :: depth_of_sediment_column = type_horizontal_standard_variable(name='depth_of_sediment_column',units='m')
   type (type_horizontal_standard_variable),parameter :: diffusivity_in_sediment_layer_1 = type_horizontal_standard_variable(name='diffusivity_in_sediment_layer_1',units='m^2/d')
   type (type_horizontal_standard_variable),parameter :: diffusivity_in_sediment_layer_2 = type_horizontal_standard_variable(name='diffusivity_in_sediment_layer_2',units='m^2/d')
   type (type_horizontal_standard_variable),parameter :: diffusivity_in_sediment_layer_3 = type_horizontal_standard_variable(name='diffusivity_in_sediment_layer_3',units='m^2/d')
   type (type_horizontal_standard_variable),parameter :: particulate_diffusivity_due_to_bioturbation = type_horizontal_standard_variable(name='particulate_diffusivity_due_to_bioturbation',units='m^2/d')
   type (type_horizontal_standard_variable),parameter :: bioturbation_depth = type_horizontal_standard_variable(name='bioturbation_depth',units='m')
   type (type_horizontal_standard_variable),parameter :: sediment_porosity = type_horizontal_standard_variable(name='sediment_porosity',units='-')
   type (type_horizontal_standard_variable),parameter :: depth_of_bottom_interface_of_layer_1 = type_horizontal_standard_variable(name='depth_of_bottom_interface_of_layer_1',units='m')
   type (type_horizontal_standard_variable),parameter :: depth_of_bottom_interface_of_layer_2 = type_horizontal_standard_variable(name='depth_of_bottom_interface_of_layer_2',units='m')
   type (type_horizontal_standard_variable),parameter :: pelagic_benthic_transfer_constant = type_horizontal_standard_variable(name='pelagic_benthic_transfer_constant',units='d/m')
   type (type_horizontal_standard_variable),parameter :: sediment_erosion = type_horizontal_standard_variable(name='sediment_erosion',units='m/d')

   ! Aggregate absorption and backscatter.
   type (type_bulk_standard_variable),parameter :: particulate_organic_absorption_coefficient = type_bulk_standard_variable(name='particulate_organic_absorption_coefficient',units='1/m',aggregate_variable=.true.)
   type (type_bulk_standard_variable),parameter :: particulate_organic_backscatter_coefficient = type_bulk_standard_variable(name='particulate_organic_backscatter_coefficient',units='1/m',aggregate_variable=.true.)

   ! Gelbstoff absorption.
   type (type_horizontal_standard_variable),parameter :: gelbstoff_absorption_from_satellite = type_horizontal_standard_variable(name='gelbstoff_absorption_from_satellite',units='1/m')

   ! Suspended Particle Matter.
   type (type_horizontal_standard_variable),parameter :: spm_from_satellite = type_horizontal_standard_variable(name='spm_from_satellite',units='g/m3')   
   
   ! Zenith angle.
   type (type_horizontal_standard_variable),parameter :: zenith_angle = type_horizontal_standard_variable(name='zenith_angle',units='degrees')

   contains

    ! temperature dependency for Q10 function
    elemental function eTq(t, q10)

        IMPLICIT NONE
        real(rk),intent(IN) :: t, q10
        real(rk)            :: eTq

        eTq = exp( log(q10) * (t-BASETEMP) / BASETEMP)

    end function eTq

    ! Michaelis-Menten saturation curve
    elemental FUNCTION MM(x, m)

        IMPLICIT NONE
        real(rk),intent(IN) :: x, m
        real(rk)            :: MM

        MM = x / (x + m + p_small)

    end function MM

    ! Convert values in 0 or 1 according to input field
    elemental FUNCTION INSW(x)

        IMPLICIT NONE
        real(rk),intent(IN) :: x
        real(rk)            :: INSW
 
        INSW = ZERO
        if (x > ZERO ) INSW=ONE 

    end function INSW

    ! Michaelis-Menten saturation curve at power
    elemental FUNCTION MM_POWER(x, m, p)

        IMPLICIT NONE
        real(rk),intent(IN) :: x, m
        integer   ,intent(IN) :: p
        real(rk)            :: MM_POWER

        MM_POWER = x**p / ( x**p+ m**p + p_small)

    end function MM_POWER

    subroutine linear_regression(x, y, n, a, b)
      
        IMPLICIT NONE
        real(rk), intent(IN) :: x(:), y(:)
        integer,  intent(IN) :: n
        real(rk), intent(OUT) :: a, b
        real(rk) :: s1,s2,s3,s4
        integer :: i
        do i=1,n
           s1=s1+x(i)
           s2=s2+x(i)**2
           s3=s3+y(i)
           s4=s4+x(i)*y(i)
        enddo
        b=((n*s4)-(s1*s3))/((n*s2)-(s1**2))
        a=(s3-(s1*b))/n
        
    end subroutine linear_regression
       subroutine calculate_integral_weights(xl, xr, n, x, w)
      integer,  intent(in)  :: n
      real(rk), intent(in)  :: x(n), xl, xr
      real(rk), intent(out) :: w(n)

      integer  :: i
      real(rk) :: deltax, f

      w = 0._rk

      if (x(1) > xr) then
         ! All points to right of desired range.
         w(1) = xr - xl
         return
      elseif (x(n) < xl) then
         ! All points to left of desired range.
         w(n) =  xr - xl
         return
      end if

      do i = 2, n
         deltax = x(i) - x(i-1)
         if (x(i) >= xl .and. x(i) <= xr) then
            ! Right-hand point in desired range.
            if (x(i-1) >= xl) then
               ! Whole segment in range
               ! Integral: 0.5*(y1+y2)*(x2-x1)
               w(i-1:i) = w(i-1:i) + 0.5_rk * deltax
            else
               ! Segment crosses left boundary
               ! y at boundary: yb = y1 + (y2-y1)/(x2-x1)*(xl-x1) = y1*(1-(xl-x1)/(x2-x1)) + y2*(xl-x1)/(x2-x1)
               ! integral = 0.5*(yb+y2)*(x2-xl)
               ! yb+y2 = y1*(1-(xl-x1)/(x2-x1)) + y2*(1+(xl-x1)/(x2-x1))
               f = (xl - x(i-1)) / deltax
               w(i-1) = w(i-1) + (1.0_rk - f) * (x(i) - xl) * 0.5_rk
               w(i  ) = w(i  ) + (1.0_rk + f) * (x(i) - xl) * 0.5_rk
            end if
         elseif (x(i) > xr .and. x(i-1) < xr) then
            ! Right-hand point beyond desired range, left-hand point before right range boundary.
            if (x(i-1) >= xl) then
               ! Segment crosses right boundary
               ! y at boundary: yb = y1 + (y2-y1)/(x2-x1)*(xr-x1) = y1*(1-(xr-x1)/(x2-x1)) + y2*(xr-x1)/(x2-x1)
               ! integral = 0.5*(y1+yb)*(xr-x1)
               ! y1+yb = y1*(2-(xr-x1)/(x2-x1)) + y2*(xr-x1)/(x2-x1)
               f = (xr - x(i-1)) / deltax
               w(i-1) = w(i-1) + (2.0_rk - f) * (xr - x(i-1)) * 0.5_rk
               w(i  ) = w(i  ) + (         f) * (xr - x(i-1)) * 0.5_rk
            else
               ! Segment crosses both boundaries
               ! y at centre: yc = y1 + (y2-y1)/(x2-x1)*(0.5*(xl+xr)-x1) = y1*(1-(0.5*(xl+xr)-x1)/(x2-x1)) + y2*(0.5*(xl+xr)-x1)/(x2-x1)
               ! integral: (xr-xl)*yc
               f = (0.5_rk*(xl+xr)-x(i-1))/deltax
               w(i-1) = w(i-1) + (1._rk-f)*(xr-xl)
               w(i  ) = w(i  ) + (      f)*(xr-xl)
            end if
         end if
      end do
      if (x(1)>xl) w(1) = w(1) + (x(1)-xl)
      if (x(n)<xr) w(n) = w(n) + (xr-x(n))
   end subroutine


end module
