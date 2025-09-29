#include "fabm_driver.h"

module ihamocc_sediment_bypass

   use fabm_types
   use ihamocc_shared

   implicit none

   private

   type, extends(type_base_model), public :: type_ihamocc_sediment_bypass
      type (type_bottom_diagnostic_variable_id), allocatable :: id_pool(:)
      type (type_diagnostic_variable_id), allocatable ::  id_flux_out(:)
      type (type_bottom_dependency_id), allocatable          :: id_flux(:)
      type (type_state_variable_id), allocatable             :: id_tracer(:)
      type (type_bottom_dependency_id)                       :: id_bdepth
      integer :: ntracers
   contains
      procedure :: initialize
      procedure :: do
   end type type_ihamocc_sediment_bypass

contains
   subroutine initialize(self, configunit)
      class (type_ihamocc_sediment_bypass), intent(inout), target :: self
      integer,                  intent(in)            :: configunit
      integer           :: i
      character(len=64) :: index
      
      ! Register number of sinking tracers to bypass sediment
      call self%get_parameter(self%ntracers, 'ntracers', '-','number of sinking tracer couplings to setup', default=0)
      
      ! Register dependencies
      call self%register_dependency(self%id_bdepth, standard_variables%bottom_depth)
      
      allocate(self%id_flux(self%ntracers))
      allocate(self%id_pool(self%ntracers))
      allocate(self%id_tracer(self%ntracers))
      allocate(self%id_flux_out(self%ntracers))
      do i=1, self%ntracers
          write(index,'(i0)') i
          call self%register_dependency(self%id_flux(i),  'flux'//trim(index),   'kmol m-2 s-1',  'bottom flux '//trim(index))
      
          ! Register diagnostics (fake state vars to receive sinking material)
          call self%register_diagnostic_variable(self%id_pool(i), 'pool'//trim(index),'kmol m-2', 'target bottom pool for pelagic tracer '//trim(index),source=source_constant, output=output_none, act_as_state_variable=.true.)
          call self%register_diagnostic_variable(self%id_flux_out(i), 'pool'//trim(index)//'flux','kmol m-2', 'flux to target bottom pool for pelagic tracer '//trim(index),source=source_do, output=output_instantaneous)
          
          ! couple to access the incoming flux
          call self%request_coupling(self%id_flux(i),'./pool'//trim(index)//'_sms_tot')
      
          ! target tracer to recieve flux evenly distributed back into the water column 
          call self%register_state_dependency(self%id_tracer(i), 'tracer'//trim(index), 'kmol m-3', 'target tracer to receive reintroduced bottom flux '//trim(index))
      end do
   end subroutine
   
   subroutine do(self, _ARGUMENTS_DO_)
      class (type_ihamocc_sediment_bypass), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_

      integer                             :: i
      real(rk)                            :: bdepth
      real(rk), dimension(self%ntracers)  :: flux, source
      _LOOP_BEGIN_
         _GET_BOTTOM_(self%id_bdepth, bdepth)
         do i=1,self%ntracers
            _GET_BOTTOM_(self%id_flux(i), flux(i)) ! get tracer settlement rate in the lowest layer
         end do
         source = flux/bdepth ! calculate water column tracer concentration increase
         do i=1,self%ntracers
            _SET_DIAGNOSTIC_(self%id_flux_out(i), flux(i))
            _ADD_SOURCE_(self%id_tracer(i),source(i)) ! add tracer to final recipient
         end do
      _LOOP_END_
   end subroutine do
end module ihamocc_sediment_bypass
