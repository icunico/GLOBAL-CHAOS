#include "fabm_driver.h"

! Fennel & Neumann 1996 NPZD model - phytoplankton component

module fussmann_Z
   use fabm_types

   implicit none

   private

   type, extends(type_base_model), public :: type_fussmann_Z
      ! Variable identifiers
      type (type_state_variable_id)      :: id_Z1
      type (type_state_variable_id)      :: id_Y1
      type (type_diagnostic_variable_id) :: id_lvZ1
      type (type_diagnostic_variable_id) :: id_grazingY1

      ! Model parameters
      real(rk) :: aY1Z1,bY1Z1,dZ1
   contains
      procedure :: initialize
      procedure :: do
   end type

contains

   subroutine initialize(self, configunit)
      class (type_fussmann_Z), intent(inout), target :: self
      integer,                        intent(in)            :: configunit

      real(rk), parameter :: d_per_s = 1.0_rk/86400.0_rk
      real(rk)            :: w_p

      self%dt = 86400._rk

      ! Store parameter values in our own derived type
      ! NB: all rates must be provided in values per day and are converted here to values per second.
      call self%get_parameter(self%aY1Z1,    'aY1Z1',    'd',             'growth',                         default=1.0_rk)
      call self%get_parameter(self%bY1Z1,    'bY1Z1',    'd',             'inverse half saturation',        default=5.0_rk)
      call self%get_parameter(self%dZ1,      'dZ1',      'd-1',           'mortality',                      default=1.0_rk)

      ! Register state variables
      call self%register_state_variable(self%id_Z1, 'DWZ1', 'd-1', 'dryweight', 0.0_rk, minimum=0.0_rk)
      call self%register_state_dependency(self%id_Y1, 'DWY1', 'd-1', 'dryweight')

      ! Register diagnostic variables
      call self%register_diagnostic_variable(self%id_lvZ1,    'growth1', 'd-2', 'growth rate')
      call self%register_diagnostic_variable(self%id_grazingY1,  'grazingY1', 'd-2', 'grazing1')

      ! Register environmental dependencies

      ! Contribute to light attentuation
   end subroutine initialize

   subroutine do(self, _ARGUMENTS_DO_)
      class (type_fussmann_Z), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_

      real(rk)            :: Z1,Y1
      real(rk)            :: grazingY1,lvZ1
      real(rk), parameter :: secs_pr_day = 86400.0_rk

      ! Enter spatial loops (if any)
      _LOOP_BEGIN_

         ! Retrieve current (local) state variable values.
         _GET_(self%id_Y1,Y1)           ! predatorY1
         _GET_(self%id_Z1,Z1)           ! predatorZ1

         ! Retrieve current environmental conditions.

         ! Light acclimation formulation based on surface light intensity.

         ! Loss rate of phytoplankton to detritus depends on local light intensity.

         ! Define some intermediate quantities that will be reused multiple times.
         lvZ1 = Z1 * ( (self%aY1Z1*Y1) / (1.0+self%bY1Z1*Y1) - self%dZ1 )
         grazingY1 = - Z1 * ( (self%aY1Z1*Y1) / (1.0+self%bY1Z1*Y1) )  
         _SET_ODE_(self%id_Z1, lvZ1)
         _SET_ODE_(self%id_Y1, grazingY1)
         _SET_DIAGNOSTIC_(self%id_lvZ1, lvZ1)
         _SET_DIAGNOSTIC_(self%id_grazingY1, grazingY1)
         
         ! Set temporal derivatives

         ! If an externally maintained ...

         ! Export diagnostic variables

      ! Leave spatial loops (if any)
   _LOOP_END_
   end subroutine do



end module fussmann_Z

!-----------------------------------------------------------------------
! Copyright Bolding & Bruggeman ApS - GNU Public License - www.gnu.org
!-----------------------------------------------------------------------
