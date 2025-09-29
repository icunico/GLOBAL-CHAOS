#include "fabm_driver.h"

! Fennel & Neumann 1996 NPZD model - phytoplankton component

module fussmann_N
   use fabm_types

   implicit none

   private

   type, extends(type_base_model), public :: type_fussmann_N
      ! Variable identifiers
      type (type_state_variable_id)      :: id_N1
      type (type_state_variable_id)      :: id_D1
      type (type_diagnostic_variable_id) :: id_lvN1

      ! Model parameters
      real(rk) :: tau 
   contains
      procedure :: initialize
      procedure :: do
   end type

contains

   subroutine initialize(self, configunit)
      class (type_fussmann_N), intent(inout), target :: self
      integer,                        intent(in)            :: configunit

      real(rk), parameter :: d_per_s = 1.0_rk/86400.0_rk
      real(rk)            :: w_p

      self%dt = 86400._rk

      ! Store parameter values in our own derived type
      call self%get_parameter(self%tau,    'tau',    'd-1',        'time nutrients',            default=2.5_rk) 
     ! NB: all rates must be provided in values per day and are converted here to values per second.

      ! Register state variables
      call self%register_state_variable(self%id_N1, 'DWN1', 'd-1', 'detritus', 0.0_rk, minimum=0.0_rk)
      call self%register_state_dependency(self%id_D1, 'DWD1', 'd-1', 'dryweight')

      ! Register diagnostic variables
      call self%register_diagnostic_variable(self%id_lvN1,    'detritus', 'd-2', 'death')


      ! Register environmental dependencies

      ! Contribute to light attentuation
   end subroutine initialize

   subroutine do(self, _ARGUMENTS_DO_)
      class (type_fussmann_N), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_

      real(rk)            :: N1,D1
      real(rk)            :: lvN1
      real(rk), parameter :: secs_pr_day = 86400.0_rk

      ! Enter spatial loops (if any)
      _LOOP_BEGIN_

         ! Retrieve current (local) state variable values.
         _GET_(self%id_N1,N1)           ! nutrients
         _GET_(self%id_D1,D1)           ! detritus

         ! Retrieve current environmental conditions.

         ! Light acclimation formulation based on surface light intensity.

         ! Loss rate of phytoplankton to detritus depends on local light intensity.

         ! Define some intermediate quantities that will be reused multiple times.
         
         lvN1=(1/self%tau)*D1
         _SET_ODE_(self%id_N1, lvN1)
         _SET_ODE_(self%id_D1, -lvN1)
         _SET_DIAGNOSTIC_(self%id_lvN1, lvN1)

         ! Set temporal derivatives

         ! If an externally maintained ...

         ! Export diagnostic variables

      ! Leave spatial loops (if any)
   _LOOP_END_
   end subroutine do



end module fussmann_N

!-----------------------------------------------------------------------
! Copyright Bolding & Bruggeman ApS - GNU Public License - www.gnu.org
!-----------------------------------------------------------------------
