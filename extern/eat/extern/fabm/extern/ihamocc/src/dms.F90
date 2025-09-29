#include "fabm_driver.h"

module ihamocc_dms

   use fabm_types
   use ihamocc_shared

   implicit none

   private

   type, extends(type_base_model), public :: type_ihamocc_dms
      type (type_dependency_id) :: id_ptho, id_hi, id_delcar, id_delsil, id_depth, id_light, id_pi_ph
      type (type_surface_dependency_id) :: id_psicomo, id_pfu10
      type (type_state_variable_id) :: id_dms
      type (type_surface_diagnostic_variable_id) ::  id_atmdms
      real(rk):: dmsp1, dmsp2, dmsp3, dmsp4, dmsp5, dmsp6, dms_gamma
      logical :: with_dmsph
   contains
      procedure :: initialize
      procedure :: do_surface
      procedure :: do
   end type type_ihamocc_dms

contains

   subroutine initialize(self, configunit)
      class (type_ihamocc_dms), intent(inout), target :: self
      integer,                  intent(in)            :: configunit
      
      call self%get_parameter(self%with_dmsph, 'with_dmsph', '-','with_dmsph',                  default=.false.)
      call self%get_parameter(self%dms_gamma,  'dms_gamma',  '-','dms_ph scaling factor',       default=0.87_rk)
      call self%get_parameter(self%dmsp6,      'dmsp6',      '-','0 half saturation microbial', default=1.e-8_rk) !Parameter are a result from kettle optimisation 02.03.04
      call self%get_parameter(self%dmsp5,      'dmsp5',      '-','production with delsil',      default=0.025_rk) !Parameter are a result from kettle optimisation 02.03.04. Following Kloster et al., 06 Table 1, but increased by a factor of ~2
      call self%get_parameter(self%dmsp4,      'dmsp4',      '-','production with delcar',      default=0.125_rk) !Parameter are a result from kettle optimisation 02.03.04. Following Kloster et al., 06 Table 1, but reduced by ~7%
      call self%get_parameter(self%dmsp3,      'dmsp3',      '-','dms parameter 3',             default=0.0864_rk) !Following Kloster et al., 06 Table 1 with 50% reduction to reduce bacterial removal and increase dms emissions
      call self%get_parameter(self%dmsp2,      'dmsp2',      '-','dms parameter 2',             default=0.0011_rk) !Following Kloster et al., 06 Table 1
      call self%get_parameter(self%dmsp1,      'dmsp1',      '-','dms parameter 1',             default=10._rk) !2*5. production with temp
                  
      call self%register_state_variable(self%id_dms, 'dms', 'kmol m-3', 'dimethyl sulfide concentration', minimum=0.0_rk)
 
      call self%register_dependency(self%id_ptho,    standard_variables%temperature)
      call self%register_dependency(self%id_psicomo, standard_variables%ice_area_fraction)
      call self%register_dependency(self%id_pfu10,   standard_variables%wind_speed)
      call self%register_dependency(self%id_depth,   standard_variables%depth)
      call self%register_dependency(self%id_light,   'light',  'W m-2',        'remaining PAR light not absorbed above')
      if (self%with_dmsph) then
          call self%register_dependency(self%id_pi_ph,   'pi_pi',  'mol kg-1',       'PI pH') ! NIC: this appears to be just the ph value of seawater(?)
          call self%register_dependency(self%id_hi,      'hi',     'mol kg-1',       'Hydrogen ion concentration')
      endif
      
      call self%register_dependency(self%id_delcar,  'delcar', 'kmol m-3 d-1', 'delcar')
      call self%register_dependency(self%id_delsil,  'delsil', 'kmol m-3 d-1', 'delsil')

      call self%register_diagnostic_variable(self%id_atmdms, 'atmdms', 'kmol m2 s-1', 'dimethyl sulfide surface flux')
   end subroutine
   
   subroutine do_surface(self, _ARGUMENTS_DO_SURFACE_)
      class (type_ihamocc_dms), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_SURFACE_

      real(rk) :: t, t2, t3, t4, ptho, psicomo, pfu10, schdms, kwdms, dmsflux, dms, scdms
            
      _SURFACE_LOOP_BEGIN_
         _GET_(self%id_dms, dms)

         _GET_(self%id_ptho, ptho)
         _GET_SURFACE_(self%id_psicomo, psicomo)
         _GET_SURFACE_(self%id_pfu10, pfu10)
         
         t = min(40._rk,max(-3._rk,ptho))
         t2   = t**2
         t3   = t**3
         t4   = t**4

         scdms = 2855.7_rk - 177.63_rk*t + 6.0438_rk*t2 - 0.11645_rk *t3 + 0.00094743_rk*t4 
         kwdms = (1._rk-psicomo) * Xconvxa * pfu10**2*(660._rk/scdms)**0.5_rk 
         
         dmsflux = kwdms*dms ! Surface flux of dms
         
         _ADD_SURFACE_FLUX_(self%id_dms, -dmsflux) ! NIC: positive flux indicates air -> water exchange; negative indicates water -> air exchange
         _SET_SURFACE_DIAGNOSTIC_(self%id_atmdms, dmsflux)
      _SURFACE_LOOP_END_
   end subroutine do_surface
   
   subroutine do(self, _ARGUMENTS_DO_)
      class (type_ihamocc_dms), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_

      real(rk) :: hi, pi_ph, temp, ptho, dms, depth, dms_bac, light, dms_ph, delsil, delcar, dmsprod, dms_uv
      
      _LOOP_BEGIN_
         _GET_(self%id_ptho, ptho)
         _GET_(self%id_dms, dms)
         _GET_(self%id_depth, depth)           
         _GET_(self%id_delsil, delsil)
         _GET_(self%id_delcar, delcar)
         temp = min(40._rk,max(-3._rk,ptho))
         dms_bac = self%dmsp3*abs(temp+3._rk)*dms*(dms/(self%dmsp6+dms))
         if (depth<=100_rk) then
             ! DMS sources/sinks
             _GET_(self%id_light, light)  
             if (self%with_dmsph) then
                 _GET_(self%id_hi, hi)
                 _GET_(self%id_pi_ph, pi_ph)
                 dms_ph  = 1._rk + (-log10(hi) - pi_ph)*self%dms_gamma
             else
                 dms_ph  = 1._rk
             endif
             dmsprod = (self%dmsp5*delsil+self%dmsp4*delcar)*(1._rk+1._rk/(temp+self%dmsp1)**2._rk)*dms_ph
             dms_uv  = self%dmsp2*light*dms
         endif
         _ADD_SOURCE_(self%id_dms, (dmsprod-dms_bac-dms_uv)/dtbgc)
      _LOOP_END_
   end subroutine do
end module ihamocc_dms
