#include "fabm_driver.h"

! Fennel & Neumann 1996 NPZD model - phytoplankton component

module fussmann_D
   use fabm_types

   implicit none

   private

   type, extends(type_base_model), public :: type_fussmann_D
      ! Variable identifiers
      type (type_state_variable_id)      :: id_D1
      type (type_state_variable_id)      :: id_C1, id_C2, id_X1, id_Y1, id_Z1
      type (type_diagnostic_variable_id) :: id_death, id_lvD1

      ! Model parameters
      real(rk) :: dC1, dC2, dX1, dY1, dZ1
   contains
      procedure :: initialize
      procedure :: do
   end type

contains

   subroutine initialize(self, configunit)
      class (type_fussmann_D), intent(inout), target :: self
      integer,                        intent(in)            :: configunit

      real(rk), parameter :: d_per_s = 1.0_rk/86400.0_rk
      real(rk)            :: w_p

      self%dt = 86400._rk

      ! Store parameter values in our own derived type
      ! NB: all rates must be provided in values per day and are converted here to values per second.
      call self%get_parameter(self%dC1,      'dC1',      'd',         'mortality',                      default=1.0_rk)
      call self%get_parameter(self%dC2,      'dC2',      'd',         'mortality',                      default=1.0_rk)
      call self%get_parameter(self%dX1,      'dX1',      'd',         'mortality',                      default=1.0_rk)
      call self%get_parameter(self%dY1,      'dY1',      'd',         'mortality',                      default=1.0_rk)
      call self%get_parameter(self%dZ1,      'dZ1',      'd',         'mortality',                      default=1.0_rk)

      ! Register state variables
      call self%register_state_variable(self%id_D1, 'DWD1', 'd-1', 'detritus', 0.0_rk, minimum=0.0_rk)
      call self%register_state_dependency(self%id_C1, 'DWC1', 'd-1', 'dryweight')
      call self%register_state_dependency(self%id_C2, 'DWC2', 'd-1', 'dryweight')
      call self%register_state_dependency(self%id_X1, 'DWX1', 'd-1', 'dryweight')
      call self%register_state_dependency(self%id_Y1, 'DWY1', 'd-1', 'dryweight')
      call self%register_state_dependency(self%id_Z1, 'DWZ1', 'd-1', 'dryweight')

      ! Register diagnostic variables
      call self%register_diagnostic_variable(self%id_death,    'death', 'd-2', 'death')

      ! Register environmental dependencies

      ! Contribute to light attentuation
   end subroutine initialize

   subroutine do(self, _ARGUMENTS_DO_)
      class (type_fussmann_D), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_

      real(rk)            :: D1,C1,C2,X1,Y1,Z1
      real(rk)            :: death
      real(rk), parameter :: secs_pr_day = 86400.0_rk

      ! Enter spatial loops (if any)
      _LOOP_BEGIN_

         ! Retrieve current (local) state variable values.
         _GET_(self%id_D1,D1)           ! detritus
         _GET_(self%id_C1,C1)           ! predatorC1
         _GET_(self%id_C2,C2)           ! predatorC2
         _GET_(self%id_X1,X1)           ! predatorX1
         _GET_(self%id_Y1,Y1)           ! predatorY1
         _GET_(self%id_Z1,Z1)           ! predatorZ1

         ! Retrieve current environmental conditions.

         ! Light acclimation formulation based on surface light intensity.

         ! Loss rate of phytoplankton to detritus depends on local light intensity.

         ! Define some intermediate quantities that will be reused multiple times.
         death = C1 * self%dC1+ C2 * self%dC2+ X1 * self%dX1+Y1 * self%dY1+Z1 * self%dZ1

         _SET_ODE_(self%id_D1, death)
         _SET_ODE_(self%id_C1, -C1 * self%dC1 )
         _SET_ODE_(self%id_C2, -C2 * self%dC2 )
         _SET_ODE_(self%id_X1, -X1 * self%dX1)
         _SET_ODE_(self%id_Y1, -Y1 * self%dY1)
         _SET_ODE_(self%id_Z1, -Z1 * self%dZ1)
         _SET_DIAGNOSTIC_(self%id_death, death)
         ! Set temporal derivatives

         ! If an externally maintained ...

         ! Export diagnostic variables

      ! Leave spatial loops (if any)
   _LOOP_END_
   end subroutine do



end module fussmann_D

!-----------------------------------------------------------------------
! Copyright Bolding & Bruggeman ApS - GNU Public License - www.gnu.org
!-----------------------------------------------------------------------
