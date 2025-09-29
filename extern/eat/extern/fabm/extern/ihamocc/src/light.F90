#include "fabm_driver.h"

module ihamocc_light

   use fabm_types
   use ihamocc_shared

   implicit none

   private

   type, extends(type_base_model), public :: type_ihamocc_light
      type (type_state_variable_id)      :: id_phy 
      type (type_dependency_id)          :: id_dz
      type (type_surface_dependency_id)  :: id_strahl
      type (type_diagnostic_variable_id) :: id_uv, id_light, id_chl
      real(rk) :: atten_w, atten_uv, ctochl
   contains
      procedure :: initialize
      procedure :: do_column
   end type type_ihamocc_light

contains

   subroutine initialize(self, configunit)
      class (type_ihamocc_light), intent(inout), target :: self
      integer,                  intent(in)            :: configunit
      
      call self%get_parameter(self%atten_w,  'atten_w',  'm-1', 'yellow substances attenuation coefficient',                   default=0.04_rk)
      call self%get_parameter(self%atten_uv, 'atten_uv', 'm-1', 'uv attenuation coefficient',                                  default=0.33_rk) 
      call self%get_parameter(self%ctochl,   'ctochl',   '-',   'Carbon to chlorophyl ratio',                                  default=60._rk)

      call self%register_state_dependency(self%id_phy, 'phy', 'kmol m-3', 'phytoplankton')
      
      call self%register_dependency(self%id_strahl,     standard_variables%surface_downwelling_shortwave_flux)
      call self%register_dependency(self%id_dz,         standard_variables%cell_thickness)
      
      call self%register_diagnostic_variable(self%id_uv,     'uv',     'W m-2', 'remaining uv light not absorbed above',source=source_do_column)
      call self%register_diagnostic_variable(self%id_light,  'light',  'W m-2', 'remaining PAR light not absorbed above',source=source_do_column)
      call self%register_diagnostic_variable(self%id_chl,    'chl',    'mg m-3', 'chlorophyl concentration',source=source_do_column)
   end subroutine
   
   subroutine do_column(self, _ARGUMENTS_DO_COLUMN_)
      class (type_ihamocc_light), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_COLUMN_

      real(rk) :: strahl, dz, atten_phyt, atten, abs_bgc, light, uv, avphy, absorption, absorption_uv, abs_uv, chl

      _GET_SURFACE_(self%id_strahl,strahl)
      absorption = 1._rk
      absorption_uv = 1._rk
      _DOWNWARD_LOOP_BEGIN_
         _GET_(self%id_dz,dz)
         _GET_(self%id_phy,avphy)
         
         ! Average light intensity in layer k
         chl = rcar*(12._rk/self%ctochl)*1.e6_rk * max(0._rk,avphy)
         atten_phyt = 0.03_rk*chl  ! phytoplankton attenuation in 1/m 
         atten = self%atten_w + atten_phyt
         abs_bgc = ((absorption/atten)*(1._rk-exp(-atten*dz)))/dz
         abs_uv  = ((absorption_uv/self%atten_uv)*(1._rk-exp(-self%atten_uv*dz)))/dz
                  
         ! Radiation intensity I_0 at the top of next layer
         absorption    = absorption    * exp(-atten*dz)
         absorption_uv = absorption_uv * exp(-self%atten_uv*dz)
         
         uv    = abs_uv*strahl
         light = abs_bgc*strahl
         
         _SET_DIAGNOSTIC_(self%id_chl,chl)
         _SET_DIAGNOSTIC_(self%id_uv,uv)
         _SET_DIAGNOSTIC_(self%id_light,light)
      _DOWNWARD_LOOP_END_
   end subroutine do_column
end module ihamocc_light
