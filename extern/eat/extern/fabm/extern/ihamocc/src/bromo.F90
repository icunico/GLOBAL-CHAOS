#include "fabm_driver.h"

module ihamocc_bromo

   use fabm_types
   use ihamocc_shared

   implicit none

   private

   type, extends(type_base_model), public :: type_ihamocc_bromo
      type (type_dependency_id) :: id_psao, id_ptho, id_hi, id_Kw, id_phosy, id_bkopal, id_depth, id_uv, id_swa_clim
      type (type_surface_dependency_id) :: id_psicomo, id_pfu10, id_atmbromo, id_ppao
      type (type_state_variable_id) :: id_oxygen, id_silica, id_bromo
      type (type_surface_diagnostic_variable_id) ::  id_bromoflx
      real(rk) :: rbro, fbro1, fbro2
   contains
      procedure :: initialize
      procedure :: do_surface
      procedure :: do
   end type type_ihamocc_bromo

contains

   subroutine initialize(self, configunit)
      class (type_ihamocc_bromo), intent(inout), target :: self
      integer,                  intent(in)            :: configunit
      
      call self%get_parameter(self%rbro, 'rbro',  '-', 'rbro',  default=2.4e-6_rk*rnit)
      call self%get_parameter(self%fbro1,'fbro1', '-', 'fbro1', default=1._rk*rnit)
      call self%get_parameter(self%fbro2,'fbro2', '-', 'fbro2', default=1._rk*rnit)
      
      call self%register_state_variable(self%id_bromo, 'bromoform', 'kmol m-3', 'Dissolved bromoform', minimum=0.0_rk)

      call self%register_dependency(self%id_depth,    standard_variables%depth)
      call self%register_dependency(self%id_ptho,     standard_variables%temperature)
      call self%register_dependency(self%id_psicomo,  standard_variables%ice_area_fraction)
      call self%register_dependency(self%id_pfu10,    standard_variables%wind_speed)
      call self%register_dependency(self%id_ppao,     standard_variables%surface_air_pressure) ! surface air pressure in pascal
      call self%register_dependency(self%id_atmbromo, 'atmbromo', '-',            'surface air bromoform mixing ratio') 
      call self%register_dependency(self%id_uv,       'uv',       'W m-2',        'remaining uv light not absorbed above')
      call self%register_dependency(self%id_hi,       'hi',       'mol kg-1',     'Hydrogen ion concentration')
      call self%register_dependency(self%id_Kw,       'kW',       'mol kg-1',     'Water dissociation product')
      call self%register_dependency(self%id_phosy,    'phosy',    'kmol m-3 s-1', 'photosynthetic production')
      call self%register_dependency(self%id_bkopal,   'bkopal',   'kmol Si m-3',  'half sat. constant for opal')
      call self%register_dependency(self%id_swa_clim, 'swa_clim', '??',           'swa climatology field')
      call self%register_dependency(self%id_phosy,    'phosy',    'kmol m-3 s-1', 'photosynthetic rate')

      call self%register_state_dependency(self%id_silica, 'silica', 'kmol m-3', 'Silicid acid (Si(OH)4)')
      
      call self%register_diagnostic_variable(self%id_bromoflx, 'bromoflx', 'kmol m-2 s-1', 'Bromoform surface flux')

   end subroutine
   
   subroutine do_surface(self, _ARGUMENTS_DO_SURFACE_)
      class (type_ihamocc_bromo), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_SURFACE_

      real(rk) :: t, t2, t3, tk, ptho, ppao, atmbrf, bromo, psicomo, pfu10, flx_bromo, kw_bromo, a_bromo, sch_bromo, atbrf
            
      _SURFACE_LOOP_BEGIN_
         _GET_(self%id_bromo, bromo)
         _GET_(self%id_ptho, ptho)
         _GET_SURFACE_(self%id_psicomo, psicomo)
         _GET_SURFACE_(self%id_atmbromo, atbrf)
         _GET_SURFACE_(self%id_ppao, ppao)
         _GET_SURFACE_(self%id_pfu10, pfu10)
         
         t = min(40._rk,max(-3._rk,ptho))
         t2   = t**2
         t3   = t**3
         tk = t + tzero
         
         sch_bromo= 4662.8_rk - 319.45_rk*t + 9.9012_rk*t2 - 0.1159_rk*t3 ! Stemmler et al. (2015; Biogeosciences) Eq. (9); Quack and Wallace (2003; GBC)
        
         a_bromo = exp(13.16_rk - 4973._rk*(1._rk/tk)) !Henry's law constant [dimensionless] for Bromoform from Quack and Wallace (2003; GBC)

         kw_bromo=(1._rk-psicomo) * 1.e-2_rk/3600._rk*(0.222_rk*pfu10**2+0.33_rk*pfu10)*(660._rk/sch_bromo)**0.5_rk ! Stemmler et al. (2015; Biogeosciences) Eq. (8) 1.e-2/3600 = conversion from [cm hr-1]/[m s-1]^2 to [ms-1]/[m s-1]^2

         
        flx_bromo=kw_bromo*(atbrf/a_bromo*1.e-12_rk*ppao*1.e-5_rk/(tk*0.083_rk) - bromo) ! Quack and Wallace (2003) eq. 1: flux = kw*(Cw - Ca/H) ; kw[m s-1]; Cw[kmol m-3]; 
                                                                                         ! Convert Ca(atbrf) from 
                                                                                         !  [pptv]    to [ppp]      by multiplying with 1e-12 (ppp = parts per part, dimensionless)
                                                                                         !  [ppp]     to [mol L-1]  by multiplying with pressure[bar]/(SST[K]*R[L bar K-1 mol-1]); R=0,083
                                                                                         !  [mol L-1] to [kmol m-3] by multiplying with 1 
        _ADD_SURFACE_FLUX_(self%id_bromo, -flx_bromo) ! NIC: positive flux indicates air -> water exchange; negative indicates water -> air exchange
        _SET_SURFACE_DIAGNOSTIC_(self%id_bromoflx, -flx_bromo)
      _SURFACE_LOOP_END_
   end subroutine do_surface
   
   subroutine do(self, _ARGUMENTS_DO_)
      class (type_ihamocc_bromo), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_

      real(rk) :: t, tk, hi, ah1, Kb1, rocbromo, lsub, bromo, Kw, phosy, bkopal, uv, swa_clim, bro_beta, bro_uv, depth, ptho, silica, avsil  
      
      _LOOP_BEGIN_
         _GET_(self%id_bromo, bromo)
         _GET_(self%id_hi, hi)
         _GET_(self%id_Kw,Kw)
         _GET_(self%id_phosy,phosy)
         _GET_(self%id_bkopal,bkopal)
         _GET_(self%id_uv,uv)
         _GET_(self%id_swa_clim,swa_clim)
         _GET_(self%id_depth, depth)             

         ! Carbon chemistry: Calculate equilibrium constants and solve for [H+] and
         ! carbonate alkalinity (ac)
         t    = min(40._rk,max(-3._rk,ptho))
         tk   = t + tzero
   
         ah1  = hi
         
         Kb1=2.05e12_rk*exp(-1.073e5_rk/(8.314_rk*tk))*dtbgc ! Degradation to hydrolysis (Eq. 2-4 of Stemmler et al., 2015) A1=1.23e17 mol min-1 => 2.05e12 kmol sec-1 
         
         lsub=7.33e-10_rk*exp(1.250713e4_rk*(1._rk/298._rk-1._rk/tk))*dtbgc ! Degradation to halogen substitution (Eq. 5-6 of Stemmler et al., 2015)

         rocbromo = -(Kb1*Kw/ah1 + lsub)*bromo !NIC: This expression is bugged (unit mismatch in Kb1*Kw/ah1) and should be rewritten, pending a decision from Jerry

         if (depth<=100_rk) then
             _GET_(self%id_silica, silica)
             _GET_(self%id_phosy, phosy)
             avsil = max(0.0_rk,silica)
             bro_beta = self%rbro*(self%fbro1*avsil/(avsil+bkopal)+self%fbro2*bkopal/(avsil+bkopal))
             
             if (swa_clim > 0._rk) then
                  bro_uv = 0.0333_rk*0.3_rk*uv/swa_clim*bromo
             else
                  bro_uv = 0.0_rk
             endif
             
             rocbromo = rocbromo + bro_beta*phosy - bro_uv
         endif
         _ADD_SOURCE_(self%id_bromo,rocbromo/dtbgc)
      _LOOP_END_
   end subroutine do
end module ihamocc_bromo
