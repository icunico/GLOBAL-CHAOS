#include "fabm_driver.h"

! Fennel & Neumann 1996 NPZD model - phytoplankton component

module fussmann_X
   use fabm_types

   implicit none

   private

   type, extends(type_base_model), public :: type_fussmann_X
      ! Variable identifiers
      type (type_state_variable_id)      :: id_X1
      type (type_state_variable_id)      :: id_C1
      type (type_state_variable_id)      :: id_C2
      type (type_state_variable_id)      :: id_P1
      type (type_state_variable_id)      :: id_P2
      type (type_state_variable_id)      :: id_P3
      type (type_diagnostic_variable_id) :: id_lvX1
      type (type_diagnostic_variable_id) :: id_grazingP1
      type (type_diagnostic_variable_id) :: id_grazingP2
      type (type_diagnostic_variable_id) :: id_grazingP3
      type (type_diagnostic_variable_id) :: id_grazingC1
      type (type_diagnostic_variable_id) :: id_grazingC2

      ! Model parameters
      real(rk) :: aP1X1,bP1X1, aP2X1,bP2X1, aP3X1,bP3X1, aC1X1,bC1X1,aC2X1,bC2X1,dX1
   contains
      procedure :: initialize
      procedure :: do
   end type

contains

   subroutine initialize(self, configunit)
      class (type_fussmann_X), intent(inout), target :: self
      integer,                        intent(in)            :: configunit

      real(rk), parameter :: d_per_s = 1.0_rk/86400.0_rk
      real(rk)            :: w_p

      self%dt = 86400._rk

      ! Store parameter values in our own derived type
      ! NB: all rates must be provided in values per day and are converted here to values per second.
      call self%get_parameter(self%aP1X1,    'aP1X1',    'd',           'growth',                         default=1.0_rk)
      call self%get_parameter(self%aP2X1,    'aP2X1',    'd',           'growth',                         default=1.0_rk)
      call self%get_parameter(self%aP3X1,    'aP3X1',    'd',           'growth',                         default=1.0_rk)
      call self%get_parameter(self%aC1X1,    'aC1X1',    'd',           'growth',                         default=1.0_rk)
      call self%get_parameter(self%aC2X1,    'aC2X1',    'd',           'growth',                         default=1.0_rk)
      call self%get_parameter(self%bP1X1,    'bP1X1',    'd',           'inverse half saturation',        default=5.0_rk)
      call self%get_parameter(self%bP2X1,    'bP2X1',    'd',           'inverse half saturation',        default=5.0_rk)
      call self%get_parameter(self%bP3X1,    'bP3X1',    'd',           'inverse half saturation',        default=5.0_rk)
      call self%get_parameter(self%bC1X1,    'bC1X1',    'd',           'inverse half saturation',        default=5.0_rk)
      call self%get_parameter(self%bC2X1,    'bC2X1',    'd',           'inverse half saturation',        default=5.0_rk)
      call self%get_parameter(self%dX1,      'dX1',      'd-1',           'mortality',                      default=1.0_rk)

      ! Register state variables
      call self%register_state_variable(self%id_X1, 'DWX1', 'd-1', 'dryweight', 0.0_rk, minimum=0.0_rk)
      call self%register_state_variable(self%id_C1, 'DWC1', 'd-1', 'dryweight', 0.0_rk, minimum=0.0_rk)
      call self%register_state_variable(self%id_C2, 'DWC2', 'd-1', 'dryweight', 0.0_rk, minimum=0.0_rk)
      call self%register_state_dependency(self%id_P1, 'DWP1', 'd-1', 'dryweight')
      call self%register_state_dependency(self%id_P2, 'DWP2', 'd-1', 'dryweight')
      call self%register_state_dependency(self%id_P3, 'DWP3', 'd-1', 'dryweight')

      ! Register diagnostic variables
      call self%register_diagnostic_variable(self%id_lvX1,    'growth1', 'd-2', 'growth rate')
      call self%register_diagnostic_variable(self%id_grazingC1,  'grazingC1', 'd-2', 'grazing1')
      call self%register_diagnostic_variable(self%id_grazingC2,  'grazingC2', 'd-2', 'grazing1')
      call self%register_diagnostic_variable(self%id_grazingP1,  'grazingP1', 'd-2', 'grazing1')
      call self%register_diagnostic_variable(self%id_grazingP2,  'grazingP2', 'd-2', 'grazing1')
      call self%register_diagnostic_variable(self%id_grazingP3,  'grazingP3', 'd-2', 'grazing1')

      ! Register environmental dependencies

      ! Contribute to light attentuation
   end subroutine initialize

   subroutine do(self, _ARGUMENTS_DO_)
      class (type_fussmann_X), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_

      real(rk)            :: X1,C1,C2,P1,P2,P3
      real(rk)            :: grazingP1,grazingP2,grazingP3,grazingC1,grazingC2,lvX1
      real(rk), parameter :: secs_pr_day = 86400.0_rk

      ! Enter spatial loops (if any)
      _LOOP_BEGIN_

         ! Retrieve current (local) state variable values.
         _GET_(self%id_X1,X1)           ! predatorX1
         _GET_(self%id_C1,C1)           ! predatoriC1
         _GET_(self%id_C2,C1)           ! predatoriC2
         _GET_(self%id_P1,P1)           ! prey1
         _GET_(self%id_P2,P1)           ! prey2
         _GET_(self%id_P3,P1)           ! prey3

         ! Retrieve current environmental conditions.

         ! Light acclimation formulation based on surface light intensity.

         ! Loss rate of phytoplankton to detritus depends on local light intensity.

         ! Define some intermediate quantities that will be reused multiple times.
         lvX1 = X1 * ( (self%aP1X1*P1+self%aP2X1*P2+self%aP3X1*P3+self%aC1X1*C1+self%aC2X1*C2) / (1.0+self%bP1X1*P1+self%bP2X1*P2+self%bP3X1*P3+self%bC1X1*C1+self%bC2X1*C2) - self%dX1 )
         grazingP1 = - X1 * ( self%aP1X1*P1 / (1.0+self%bP1X1*P1+self%bP2X1*P2+self%bP3X1*P3+self%bC1X1*C1+self%bC2X1*C2) )  
         grazingP2 = - X1 * ( self%aP2X1*P2 / (1.0+self%bP1X1*P1+self%bP2X1*P2+self%bP3X1*P3+self%bC1X1*C1+self%bC2X1*C2) )  
         grazingP3 = - X1 * ( self%aP3X1*P3 / (1.0+self%bP1X1*P1+self%bP2X1*P2+self%bP3X1*P3+self%bC1X1*C1+self%bC2X1*C2) )  
         grazingC1 = - X1 * ( self%aC1X1*C1 / (1.0+self%bP1X1*P1+self%bP2X1*P2+self%bP3X1*P3+self%bC1X1*C1+self%bC2X1*C2) )  
         grazingC2 = - X1 * ( self%aC2X1*C2 / (1.0+self%bP1X1*P1+self%bP2X1*P2+self%bP3X1*P3+self%bC1X1*C1+self%bC2X1*C2) )  
         _SET_ODE_(self%id_X1, lvX1)
         _SET_ODE_(self%id_C1, grazingC1)
         _SET_ODE_(self%id_C2, grazingC2)
         _SET_ODE_(self%id_P1, grazingP1)
         _SET_ODE_(self%id_P2, grazingP2)
         _SET_ODE_(self%id_P3, grazingP3)
         _SET_DIAGNOSTIC_(self%id_lvX1, lvX1)
         _SET_DIAGNOSTIC_(self%id_grazingC1, grazingC1)
         _SET_DIAGNOSTIC_(self%id_grazingC2, grazingC2)
         _SET_DIAGNOSTIC_(self%id_grazingP1, grazingP1)
         _SET_DIAGNOSTIC_(self%id_grazingP2, grazingP2)
         _SET_DIAGNOSTIC_(self%id_grazingP3, grazingP3)
         
         ! Set temporal derivatives

         ! If an externally maintained ...

         ! Export diagnostic variables

      ! Leave spatial loops (if any)
   _LOOP_END_
   end subroutine do



end module fussmann_X

!-----------------------------------------------------------------------
! Copyright Bolding & Bruggeman ApS - GNU Public License - www.gnu.org
!-----------------------------------------------------------------------
