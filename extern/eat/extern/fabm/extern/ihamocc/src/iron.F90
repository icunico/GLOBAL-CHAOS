#include "fabm_driver.h"

module ihamocc_iron

   use fabm_types
   use ihamocc_shared

   implicit none

   private

   type, extends(type_base_model), public :: type_ihamocc_iron
      type (type_surface_dependency_id) :: id_dust
      type (type_state_variable_id) :: id_iron, id_fdust
   contains
      procedure :: initialize
      procedure :: do_surface
      procedure :: do
   end type type_ihamocc_iron

contains

   subroutine initialize(self, configunit)
      class (type_ihamocc_iron), intent(inout), target :: self
      integer,                  intent(in)            :: configunit
      
      call self%register_state_variable(self%id_iron,  'iron',  'kmol m-3', 'dissolved iron', minimum=0.0_rk)
      call self%add_to_aggregate_variable(standard_variables%total_iron,       self%id_iron, scale_factor=1e9_rk)
      
      call self%register_state_variable(self%id_fdust, 'fdust', 'kg m-3',   'non-aggregated dust deposition', minimum=0.0_rk)

      call self%register_dependency(self%id_dust, 'dust', 'kg m-2 s-1', 'dust deposition rate')
   end subroutine
   
   subroutine do_surface(self, _ARGUMENTS_DO_SURFACE_)
      class (type_ihamocc_iron), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_SURFACE_

      real(rk) :: dust, perc_diron, roc_iron
      
      _SURFACE_LOOP_BEGIN_
         _GET_SURFACE_(self%id_dust, dust)
         
         perc_diron = 0.6_rk * 0.035_rk * 0.01_rk / 55.85_rk
         roc_iron = perc_diron*dust
         
         _ADD_SURFACE_FLUX_(self%id_fdust, dust) ! NIC: positive flux indicates air -> water exchange; negative indicates water -> air exchange
         _ADD_SURFACE_FLUX_(self%id_iron, roc_iron) ! NIC: positive flux indicates air -> water exchange; negative indicates water -> air exchange
      _SURFACE_LOOP_END_
   end subroutine do_surface

   subroutine do(self, _ARGUMENTS_DO_)
      class (type_ihamocc_iron), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_

      real(rk) :: fesoly, relaxfe, iron
      
      _LOOP_BEGIN_
         _GET_(self%id_iron, iron)
         
         fesoly=0.5_rk*1.e-9_rk            ! max. diss. iron concentration in deep water 
         relaxfe = 0.05_rk/365._rk/dtbgc      
         _ADD_SOURCE_(self%id_iron,- relaxfe*MAX(iron-fesoly,0._rk))
      _LOOP_END_
   end subroutine do
end module ihamocc_iron
