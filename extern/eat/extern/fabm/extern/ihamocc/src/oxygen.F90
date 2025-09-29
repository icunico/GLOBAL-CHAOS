#include "fabm_driver.h"

module ihamocc_oxygen

   use fabm_types
   use ihamocc_shared

   implicit none

   private

   type, extends(type_base_model), public :: type_ihamocc_oxygen
      type (type_dependency_id) :: id_psao, id_ptho
      type (type_surface_dependency_id) :: id_psicomo, id_pfu10, id_ppao
      type (type_state_variable_id) :: id_oxygen
      type (type_surface_diagnostic_variable_id) ::  id_oxflux
      type (type_diagnostic_variable_id) :: id_satoxy
   contains
      procedure :: initialize
      procedure :: do_surface
      procedure :: do
   end type type_ihamocc_oxygen

      real(rk), parameter :: OX0  = -173.4292_rk !VOLUMETRIC SOLUBILITY CONSTANTS FOR O2 IN ML/L from moist air at one atm total pressure. 
      real(rk), parameter :: OX1  = 249.6339_rk  !Table 2 in WEISS, R. F. (1970) THE SOLUBILITY OF NITROGEN OXYGEN AND ARGON IN WATER AND SEAWATER. 
      real(rk), parameter :: OX2  = 143.3483_rk  !DEEP-SEA RESEARCH, VOL. 17, 721-735.
      real(rk), parameter :: OX3  = -21.8492_rk  !
      real(rk), parameter :: OX4  = -0.033096_rk !
      real(rk), parameter :: OX5  = 0.014259_rk  !
      real(rk), parameter :: OX6  = -0.0017_rk   !
      real(rk), parameter :: ato2 = 196800._rk   ! atmospheric oxygen concentration in PPM        

contains

   subroutine initialize(self, configunit)
      class (type_ihamocc_oxygen), intent(inout), target :: self
      integer,                  intent(in)            :: configunit
      
      call self%register_state_variable(self%id_oxygen, 'oxygen', 'kmol m-3', 'Dissolved oxygen', minimum=0.0_rk)

      call self%register_dependency(self%id_psao,    standard_variables%practical_salinity)
      call self%register_dependency(self%id_ptho,    standard_variables%temperature)
      call self%register_dependency(self%id_psicomo, standard_variables%ice_area_fraction)
      call self%register_dependency(self%id_pfu10,   standard_variables%wind_speed)
      call self%register_dependency(self%id_ppao,    standard_variables%surface_air_pressure) ! surface air pressure in pascal
            
      call self%register_diagnostic_variable(self%id_oxflux, 'oxflux', 'kmol m-2 s-1', 'oxygen surface flux')
      call self%register_diagnostic_variable(self%id_satoxy, 'satoxy', 'kmol m-3',     'oxygen solubility')
   end subroutine
   
   subroutine do_surface(self, _ARGUMENTS_DO_SURFACE_)
      class (type_ihamocc_oxygen), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_SURFACE_

      real(rk) :: oxy, t, t2, t3, t4, tk, tk100, s, psao, ptho, oxflux, sco2, satoxy, rpp0, ppao, oxygen, psicomo
      real(rk) :: pfu10, kwo2
            
      _SURFACE_LOOP_BEGIN_
         _GET_(self%id_oxygen, oxygen)
         _GET_(self%id_ptho, ptho)
         _GET_(self%id_psao, psao)
         _GET_SURFACE_(self%id_psicomo, psicomo)
         _GET_SURFACE_(self%id_ppao, ppao)
         _GET_SURFACE_(self%id_pfu10, pfu10)
         
         t = min(40._rk,max(-3._rk,ptho))
         t2   = t**2._rk
         t3   = t**3._rk
         t4   = t**4._rk
         tk = t + tzero
         tk100 = tk/100.0_rk
         s = min(40._rk,max( 25._rk,psao))
         sco2  = 1920.4_rk - 135.6_rk *t + 5.2122_rk*t2 - 0.10939_rk *t3 + 0.00093777_rk*t4 ! Schmidt numbers according to Wanninkhof (2014), Table 1
         
         ! solubility of O2 (Weiss, R.F. 1970, Deep-Sea Res., 17, 721-735) for moist air at 1 atm; multiplication with oxyco converts to kmol/m^3/atm
         oxy = ox0+ox1/tk100+ox2*log(tk100)+ox3*tk100+s*(ox4+ox5*tk100+ox6*tk100**2._rk)
		 satoxy = exp(oxy)*oxyco
         
         kwo2  = (1._rk-psicomo) * Xconvxa * pfu10**2._rk*(660._rk/sco2)**0.5_rk 
         rpp0 = ppao/atm2pa
         
         ! Surface flux of oxygen
		 oxflux=kwo2*(oxygen-satoxy*(ato2/196800._rk)*rpp0) ! originally multiplied by dtbgc (ts in s) to get absolute change. Removed as FABM rates-of-change has units s-1
         _SET_SURFACE_DIAGNOSTIC_(self%id_oxflux, oxflux)
         _ADD_SURFACE_FLUX_(self%id_oxygen, -oxflux) ! NIC: positive flux indicates air -> water exchange; negative indicates water -> air exchange
      _SURFACE_LOOP_END_
   end subroutine do_surface  
   
   subroutine do(self, _ARGUMENTS_DO_)
      class (type_ihamocc_oxygen), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_

      real(rk) :: ptho, psao, t, tk, tk100, s, oxy, satoxy
      
      _LOOP_BEGIN_
         _GET_(self%id_ptho, ptho)
         _GET_(self%id_psao, psao)
         
         t = min(40._rk,max(-3._rk,ptho))
         tk = t + tzero
         tk100 = tk/100.0_rk
         s = min(40._rk,max( 25._rk,psao))

         ! solubility of O2 (Weiss, R.F. 1970, Deep-Sea Res., 17, 721-735) for moist air at 1 atm; multiplication with oxyco converts to kmol/m^3/atm
         oxy = ox0+ox1/tk100+ox2*log(tk100)+ox3*tk100+s*(ox4+ox5*tk100+ox6*tk100**2._rk)
		 satoxy = exp(oxy)*oxyco                  
         
         _SET_DIAGNOSTIC_(self%id_satoxy,satoxy)
      _LOOP_END_
   end subroutine do
end module ihamocc_oxygen
