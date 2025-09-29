#include "fabm_driver.h"

module ihamocc_cfc

   use fabm_types
   use ihamocc_shared

   implicit none

   private

   type, extends(type_base_model), public :: type_ihamocc_cfc
      type (type_dependency_id) :: id_psao, id_ptho
      type (type_surface_dependency_id) :: id_psicomo, id_pfu10, id_ppao, id_atm_cfc11, id_atm_cfc12, id_atm_sf6
      type (type_state_variable_id) :: id_cfc11, id_cfc12, id_sf6
      type (type_surface_diagnostic_variable_id) ::  id_atmf11, id_atmf12, id_atmsf6
   contains
      procedure :: initialize
      procedure :: do_surface
   end type type_ihamocc_cfc

contains

   subroutine initialize(self, configunit)
      class (type_ihamocc_cfc), intent(inout), target :: self
      integer,                  intent(in)            :: configunit
      
      call self%register_state_variable(self%id_cfc11, 'cfc11', 'kmol m-3', 'Dissolved cfc 11 gas', minimum=0.0_rk)
      call self%register_state_variable(self%id_cfc12, 'cfc12', 'kmol m-3', 'Dissolved cfc 12 gas', minimum=0.0_rk)
      call self%register_state_variable(self%id_sf6,   'sf6',   'kmol m-3', 'Dissolved sf 6 gas', minimum=0.0_rk)

      call self%register_dependency(self%id_ppao,      standard_variables%surface_air_pressure) ! surface air pressure in pascal
      call self%register_dependency(self%id_psao,      standard_variables%practical_salinity)
      call self%register_dependency(self%id_ptho,      standard_variables%temperature)
      call self%register_dependency(self%id_psicomo,   standard_variables%ice_area_fraction)
      call self%register_dependency(self%id_pfu10,     standard_variables%wind_speed)
      call self%register_dependency(self%id_atm_cfc11, 'atm_cfc11', 'ppt', 'atmospheric cfc11 concentration')
      call self%register_dependency(self%id_atm_cfc12, 'atm_cfc12', 'ppt', 'atmospheric cfc12 concentration')
      call self%register_dependency(self%id_atm_sf6,   'atm_sf6',   'ppt', 'atmospheric sf6 concentration') 
      
      call self%register_diagnostic_variable(self%id_atmf11, 'atmf11', 'kmol m-2 s-1', 'cfc 11 surface flux')
      call self%register_diagnostic_variable(self%id_atmf12, 'atmf12', 'kmol m-2 s-1', 'cfc 12 surface flux')
      call self%register_diagnostic_variable(self%id_atmsf6, 'atmsf6', 'kmol m-2 s-1', 'sf 6 surface flux')
   end subroutine
   
   subroutine do_surface(self, _ARGUMENTS_DO_SURFACE_)
      class (type_ihamocc_cfc), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_SURFACE_

      real(rk) :: t, t2, t3, t4, tk, tk100, s, ptho, psao, psicomo, pfu10, sch_11, sch_12, sch_sf, atm_cfc11, atm_cfc12, atm_sf6
      real(rk) :: ppao, a_11, a_12, a_sf, kw_11, kw_12, kw_sf, cfc11, cfc12, sf6, flx11, flx12, flxsf
            
      _SURFACE_LOOP_BEGIN_
         _GET_SURFACE_(self%id_atm_cfc11, atm_cfc11)
         _GET_SURFACE_(self%id_atm_cfc12, atm_cfc12)
         _GET_SURFACE_(self%id_atm_sf6, atm_sf6)
         _GET_(self%id_cfc11, cfc11)
         _GET_(self%id_cfc12, cfc12)
         _GET_(self%id_sf6, sf6)

         _GET_(self%id_ptho, ptho)
         _GET_(self%id_psao, psao)
         _GET_SURFACE_(self%id_psicomo, psicomo)
         _GET_SURFACE_(self%id_ppao, ppao)
         _GET_SURFACE_(self%id_pfu10, pfu10)
         
         t = min(40._rk,max(-3._rk,ptho))
         t2   = t**2
         t3   = t**3
         t4   = t**4
         tk = t + tzero
         tk100 = tk/100.0_rk
         s = min(40._rk,max( 25._rk,psao))
         
         sch_11= 3579.2_rk - 222.63_rk*t + 7.5749_rk*t2 - 0.14595_rk *t3 + 0.0011874_rk *t4
         sch_12= 3828.1_rk - 249.86_rk*t + 8.7603_rk*t2 - 0.1716_rk  *t3 + 0.001408_rk  *t4
         sch_sf= 3177.5_rk - 200.57_rk*t + 6.8865_rk*t2 - 0.13335_rk *t3 + 0.0010877_rk *t4
         
         a_11 = exp(-229.9261_rk + 319.6552_rk*(100._rk/tk) + 119.4471_rk*log(tk100)  & ! solubility of cfc11,12 (mol/(l*atm)) (Warner and Weiss 1985) and sf6 from eq. 6 of Bullister et al. (2002) These are the alpha in (1b) of the ocmpic2 howto
         &         -1.39165_rk*(tk100)**2 + s*(-0.142382_rk + 0.091459_rk*(tk100)  &
         &         -0.0157274_rk*(tk100)**2)) 
         a_12 = exp(-218.0971_rk + 298.9702_rk*(100._rk/tk) + 113.8049_rk*log(tk100)  &
         &         -1.39165_rk*(tk100)**2 + s*(-0.143566_rk + 0.091015_rk*(tk100)  &
         &         -0.0153924_rk*(tk100)**2)) 
         a_sf = exp(-80.0343_rk  + 117.232_rk *(100._rk/tk) +  29.5817_rk*log(tk100)  &
         &         +s*(0.033518_rk-0.0373942_rk*(tk100)+0.00774862_rk*(tk100)**2)) 

         a_11 = 1e-12_rk * a_11 ! conversion from mol/(l * atm) to kmol/(m3 * pptv) 
         a_12 = 1e-12_rk * a_12
         a_sf = 1e-12_rk * a_sf
      
         kw_11 = (1._rk-psicomo) * Xconvxa * pfu10**2*(660._rk/sch_11)**0.5_rk
         kw_12 = (1._rk-psicomo) * Xconvxa * pfu10**2*(660._rk/sch_12)**0.5_rk
         kw_sf = (1._rk-psicomo) * Xconvxa * pfu10**2*(660._rk/sch_sf)**0.5_rk

         ! Use conversion of 9.86923e-6 [std atm / Pascal]
         flx11 = kw_11*(a_11*atm_cfc11*ppao*9.86923_rk*1.e-6_rk-cfc11) ! Surface flux of cfc11
         flx12 = kw_12*(a_12*atm_cfc12*ppao*9.86923_rk*1.e-6_rk-cfc12) ! Surface flux of cfc12
         flxsf = kw_sf*(a_sf*atm_sf6*ppao*9.86923_rk*1.e-6_rk-sf6) ! Surface flux of sf6
         
         _ADD_SURFACE_FLUX_(self%id_cfc11, flx11) ! NIC: positive flux indicates air -> water exchange; negative indicates water -> air exchange
         _ADD_SURFACE_FLUX_(self%id_cfc12, flx12) ! NIC: positive flux indicates air -> water exchange; negative indicates water -> air exchange
         _ADD_SURFACE_FLUX_(self%id_sf6, flxsf) ! NIC: positive flux indicates air -> water exchange; negative indicates water -> air exchange
         _SET_SURFACE_DIAGNOSTIC_(self%id_atmf11, flx11)
         _SET_SURFACE_DIAGNOSTIC_(self%id_atmf12, flx12)
         _SET_SURFACE_DIAGNOSTIC_(self%id_atmsf6, flxsf)
      _SURFACE_LOOP_END_
   end subroutine do_surface
end module ihamocc_cfc
