#include "fabm_driver.h"

module ihamocc_tracer
   use fabm_types
   use fabm_particle

   use ihamocc_shared

   implicit none

   private

   type, extends(type_particle_model), public :: type_ihamocc_tracer
      type (type_state_variable_id) :: id_c, id_fe, id_si
   contains
      procedure :: initialize
   end type

contains

   subroutine initialize(self, configunit)
      class (type_ihamocc_tracer), intent(inout), target :: self
      integer,                    intent(in)            :: configunit

      logical :: has_carbon, has_nitrogen, has_phosphorus, has_silicon, has_iron

      call self%get_parameter(has_carbon,     'has_carbon',     '', 'tracer contains carbon',     default=.false.)
      call self%get_parameter(has_nitrogen,   'has_nitrogen',   '', 'tracer contains nitrogen',   default=.false.)
      call self%get_parameter(has_phosphorus, 'has_phosphorus', '', 'tracer contains phosphorus', default=.false.)
      call self%get_parameter(has_silicon,    'has_silicon',    '', 'tracer contains silicon',    default=.false.)
      call self%get_parameter(has_iron,       'has_iron',       '', 'tracer contains iron',       default=.false.)

      if (has_phosphorus) then
          call self%register_state_variable(self%id_c, 'c', 'kmol P m-3', 'concentration in P units', minimum=0.0_rk)
          call self%add_to_aggregate_variable(standard_variables%total_phosphorus, self%id_c, scale_factor=1e6_rk)
          if (has_carbon)     call self%add_to_aggregate_variable(standard_variables%total_carbon,     self%id_c, scale_factor=rcar * 1e6_rk)
          if (has_nitrogen)   call self%add_to_aggregate_variable(standard_variables%total_nitrogen,   self%id_c, scale_factor=rnit * 1e6_rk)
          if (has_iron)       call self%add_to_aggregate_variable(standard_variables%total_iron,       self%id_c, scale_factor=riron * 1e9_rk)
      endif
      
      if (has_silicon) then
         call self%register_state_variable(self%id_si, 'si', 'kmol Si m-3', 'concentration in silicon units', minimum=0.0_rk)
         call self%add_to_aggregate_variable(standard_variables%total_silicate, self%id_si, scale_factor=1e6_rk)
      end if
   end subroutine initialize

end module