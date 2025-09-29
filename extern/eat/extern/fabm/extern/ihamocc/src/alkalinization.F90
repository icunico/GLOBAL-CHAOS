#include "fabm_driver.h"

module ihamocc_alkalinization

   use fabm_types
   use ihamocc_shared

   implicit none

   private

   type, extends(type_base_model), public :: type_ihamocc_alkalinization
      type (type_surface_dependency_id) :: id_oafx
      type (type_state_variable_id) :: id_alkali
   contains
      procedure :: initialize
      procedure :: do_surface
   end type type_ihamocc_alkalinization

contains

   subroutine initialize(self, configunit)
      class (type_ihamocc_alkalinization), intent(inout), target :: self
      integer,                  intent(in)            :: configunit
      
      call self%register_state_dependency(self%id_alkali, 'alkali', 'kmol m-3', 'Alkalinity')

      call self%register_dependency(self%id_oafx, 'oafx', 'kmol m-2 s-1', 'ocean alkalinization flux')
   end subroutine
   
   subroutine do_surface(self, _ARGUMENTS_DO_SURFACE_)
      class (type_ihamocc_alkalinization), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_SURFACE_

      real(rk) :: oafx
      
      _SURFACE_LOOP_BEGIN_
         _GET_SURFACE_(self%id_oafx, oafx)
         
         _ADD_SURFACE_FLUX_(self%id_alkali, oafx) ! NIC: positive flux indicates air -> water exchange; negative indicates water -> air exchange
      _SURFACE_LOOP_END_
   end subroutine do_surface
end module ihamocc_alkalinization
