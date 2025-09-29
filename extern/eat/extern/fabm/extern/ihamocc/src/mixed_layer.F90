#include "fabm_driver.h"

module ihamocc_mixed_layer !this is basically a copy of the pisces implementation of mixed layer depth

   use fabm_types
   use ihamocc_shared

   implicit none

   private

   type, extends(type_base_model), public :: type_ihamocc_mixed_layer
      type (type_surface_diagnostic_variable_id) :: id_kmle
      type (type_dependency_id)                  :: id_depth, id_avt
      type (type_bottom_dependency_id)           :: id_bdepth
      real(rk) :: avt_c, minh
   contains
      procedure :: initialize
      procedure :: do_column
   end type type_ihamocc_mixed_layer

contains

   subroutine initialize(self, configunit)
      class (type_ihamocc_mixed_layer), intent(inout), target :: self
      integer,                  intent(in)            :: configunit
      
      call self%get_parameter(self%avt_c, 'avt_c', 'm2 s-1', 'critical vertical diffusivity', default=5.e-4_rk)
      call self%get_parameter(self%minh,  'minh',  'm',      'minimum thickness',             default=10._rk)
      
      call self%register_dependency(self%id_bdepth, standard_variables%bottom_depth)
      call self%register_dependency(self%id_depth,  standard_variables%depth)
      call self%register_dependency(self%id_avt,    type_interior_standard_variable(name='vertical_tracer_diffusivity', units='m2 s-1'))
      
      call self%register_diagnostic_variable(self%id_kmle, 'kmle', 'm', 'mixed layer depth',source=source_do_column)
   end subroutine
   
   subroutine do_column(self, _ARGUMENTS_DO_COLUMN_)
      class (type_ihamocc_mixed_layer), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_COLUMN_

      real(rk) :: h, avt, depth
      
      _GET_BOTTOM_(self%id_bdepth, h)
      _UPWARD_LOOP_BEGIN_
         _GET_(self%id_avt, avt)
         _GET_(self%id_depth, depth)
         if (avt < self%avt_c .and. depth > self%minh) h = depth
      _UPWARD_LOOP_END_
      _SET_SURFACE_DIAGNOSTIC_(self%id_kmle, h)
   end subroutine do_column
end module ihamocc_mixed_layer
