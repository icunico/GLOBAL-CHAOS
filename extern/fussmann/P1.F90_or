#include "fabm_driver.h"

! Fennel & Neumann 1996 NPZD model - phytoplankton component

module fussmann_P1
   use fabm_types

   implicit none

   private

   type, extends(type_base_model), public :: type_fussmann_P1
      ! Variable identifiers
      type (type_state_variable_id)      :: id_P1
      type (type_state_variable_id)      :: id_P2
      type (type_state_variable_id)      :: id_P3
      type (type_diagnostic_variable_id) :: id_lvP1
      type (type_diagnostic_variable_id) :: id_lvP2
      type (type_diagnostic_variable_id) :: id_lvP3

      ! Model parameters
      real(rk) :: r1,K1,r2,K2,r3,K3 
   contains
      procedure :: initialize
      procedure :: do
   end type

contains

   subroutine initialize(self, configunit)
      class (type_fussmann_P1), intent(inout), target :: self
      integer,                        intent(in)            :: configunit

      real(rk), parameter :: d_per_s = 1.0_rk/86400.0_rk
      real(rk)            :: w_p

      self%dt = 86400._rk

      ! Store parameter values in our own derived type
      ! NB: all rates must be provided in values per day and are converted here to values per second.
      call self%get_parameter(self%r1,    'r1',    'd-1',        'growth rate ',            default=2.5_rk)
      call self%get_parameter(self%r2,    'r2',    'd-1',        'growth rate ',            default=2.5_rk)
      call self%get_parameter(self%r3,    'r3',    'd-1',        'growth rate ',            default=2.5_rk)
      call self%get_parameter(self%K1,    'K1',    'mg DW m-1',  'carrying capacity',       default=1.0_rk)
      call self%get_parameter(self%K2,    'K2',    'mg DW m-1',  'carrying capacity',       default=1.0_rk)
      call self%get_parameter(self%K3,    'K3',    'mg DW m-1',  'carrying capacity',       default=1.0_rk)

      ! Register state variables
      call self%register_state_variable(self%id_P1, 'DW1', 'd-1', 'dryweight', 0.0_rk, minimum=0.0_rk)
      call self%register_state_variable(self%id_P2, 'DW2', 'd-1', 'dryweight', 0.0_rk, minimum=0.0_rk)
      call self%register_state_variable(self%id_P3, 'DW3', 'd-1', 'dryweight', 0.0_rk, minimum=0.0_rk)
!      call self%register_state_variable(self%id_prey2, 'DW2', 'd-1', 'dryweight', 0.0_rk, minimum=0.0_rk)


      ! Register diagnostic variables
      call self%register_diagnostic_variable(self%id_lvP1,  'growth1', 'd-2', 'growth rate')
      call self%register_diagnostic_variable(self%id_lvP2,  'growth2', 'd-2', 'growth rate')
      call self%register_diagnostic_variable(self%id_lvP3,  'growth3', 'd-2', 'growth rate')
!      call self%register_diagnostic_variable(self%id_lvprey2,  'growth2', 'd-2', 'growth rate')


      ! Register environmental dependencies

      ! Contribute to light attentuation
   end subroutine initialize

   subroutine do(self, _ARGUMENTS_DO_)
      class (type_fussmann_P1), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_

      real(rk)            :: P1,P2,P3
      real(rk)            :: lvP1,lvP2,lvP3
      real(rk), parameter :: secs_pr_day = 86400.0_rk

      ! Enter spatial loops (if any)
      _LOOP_BEGIN_

         ! Retrieve current (local) state variable values.
         _GET_(self%id_P1,P1)         ! prey 1
         _GET_(self%id_P2,P2)         ! prey 2
         _GET_(self%id_P3,P3)         ! prey 3

         ! Retrieve current environmental conditions.

         ! Light acclimation formulation based on surface light intensity.

         ! Loss rate of phytoplankton to detritus depends on local light intensity.

         ! Define some intermediate quantities that will be reused multiple times.
         lvP1 = P1*self%r1*(1.0-P1/self%K1)
         lvP2 = P2*self%r2*(1.0-P2/self%K2)
         lvP3 = P3*self%r3*(1.0-P3/self%K3)
         _SET_ODE_(self%id_P1, lvP1)
         _SET_ODE_(self%id_P2, lvP2)
         _SET_ODE_(self%id_P3, lvP3)
         _SET_DIAGNOSTIC_(self%id_lvP1, lvP1)
         _SET_DIAGNOSTIC_(self%id_lvP2, lvP2)
         _SET_DIAGNOSTIC_(self%id_lvP3, lvP3)

         ! Set temporal derivatives

         ! If an externally maintained ...

         ! Export diagnostic variables

      ! Leave spatial loops (if any)
   _LOOP_END_
   end subroutine do


end module fussmann_P1

!-----------------------------------------------------------------------
! Copyright Bolding & Bruggeman ApS - GNU Public License - www.gnu.org
!-----------------------------------------------------------------------

