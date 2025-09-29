#include "fabm_driver.h"

! Fennel & Neumann 1996 NPZD model - phytoplankton component

module fussmann_C1
   use fabm_types

   implicit none

   private

   type, extends(type_base_model), public :: type_fussmann_C1
      ! Variable identifiers
      type (type_state_variable_id)      :: id_C1
      type (type_state_variable_id)      :: id_C2
      type (type_state_variable_id)      :: id_P1
      type (type_state_variable_id)      :: id_P2
      type (type_state_variable_id)      :: id_P3
      type (type_diagnostic_variable_id) :: id_lvC1
      type (type_diagnostic_variable_id) :: id_lvC2
      type (type_diagnostic_variable_id) :: id_grazingP1
      type (type_diagnostic_variable_id) :: id_grazingP2
      type (type_diagnostic_variable_id) :: id_grazingP3

      ! Model parameters
      real(rk) :: aP1C1,bP1C1,aP2C1,bP2C1,aP3C1,bP3C1,dC1,aP1C2,bP1C2,aP2C2,bP2C2,aP3C2,bP3C2,dC2
   contains
      procedure :: initialize
      procedure :: do
   end type

contains

   subroutine initialize(self, configunit)
      class (type_fussmann_C1), intent(inout), target :: self
      integer,                        intent(in)            :: configunit

      real(rk), parameter :: d_per_s = 1.0_rk/86400.0_rk
      real(rk)            :: w_p

      self%dt = 86400._rk

      ! Store parameter values in our own derived type
      ! NB: all rates must be provided in values per day and are converted here to values per second.
      call self%get_parameter(self%aP1C1,    'aP1C1',    'd',           'growth',                         default=1.0_rk)
      call self%get_parameter(self%aP2C1,    'aP2C1',    'd',           'growth',                         default=1.0_rk)
      call self%get_parameter(self%aP3C1,    'aP3C1',    'd',           'growth',                         default=1.0_rk)
      call self%get_parameter(self%BP1C1,    'bP1C1',    'd',           'inverse half saturation',        default=5.0_rk)
      call self%get_parameter(self%BP2C1,    'bP2C1',    'd',           'inverse half saturation',        default=5.0_rk)
      call self%get_parameter(self%BP3C1,    'bP3C1',    'd',           'inverse half saturation',        default=5.0_rk)
      call self%get_parameter(self%dC1,      'dC1',      'd',           'mortality',                      default=1.0_rk)
      call self%get_parameter(self%aP1C2,    'aP1C2',    'd',           'growth',                         default=1.0_rk)
      call self%get_parameter(self%aP2C2,    'aP2C2',    'd',           'growth',                         default=1.0_rk)
      call self%get_parameter(self%aP3C2,    'aP3C2',    'd',           'growth',                         default=1.0_rk)
      call self%get_parameter(self%BP1C2,    'bP1C2',    'd',           'inverse half saturation',        default=5.0_rk)
      call self%get_parameter(self%BP2C2,    'bP2C2',    'd',           'inverse half saturation',        default=5.0_rk)
      call self%get_parameter(self%BP3C2,    'bP3C2',    'd',           'inverse half saturation',        default=5.0_rk)
      call self%get_parameter(self%dC2,      'dC2',      'd-1',           'mortality',                      default=1.0_rk)

      ! Register state variables
      call self%register_state_variable(self%id_C1, 'DWC1', 'd-1', 'dryweight', 0.0_rk, minimum=0.0_rk)
      call self%register_state_variable(self%id_C2, 'DWC2', 'd-1', 'dryweight', 0.0_rk, minimum=0.0_rk)
      call self%register_state_dependency(self%id_P1, 'DWP1', 'd-1', 'dryweight')
      call self%register_state_dependency(self%id_P2, 'DWP2', 'd-1', 'dryweight')
      call self%register_state_dependency(self%id_P3, 'DWP3', 'd-1', 'dryweight')

      ! Register diagnostic variables
      call self%register_diagnostic_variable(self%id_lvC1,    'growth1', 'd-2', 'growth rate')
      call self%register_diagnostic_variable(self%id_lvC2,    'growth2', 'd-2', 'growth rate')
      call self%register_diagnostic_variable(self%id_grazingP1,  'grazingP1', 'd-2', 'grazing1')
      call self%register_diagnostic_variable(self%id_grazingP2,  'grazingP2', 'd-2', 'grazing1')
      call self%register_diagnostic_variable(self%id_grazingP3,  'grazingP3', 'd-2', 'grazing1')

      ! Register environmental dependencies

      ! Contribute to light attentuation
   end subroutine initialize

   subroutine do(self, _ARGUMENTS_DO_)
      class (type_fussmann_C1), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_

      real(rk)            :: C1,C2,P1,P2,P3
      real(rk)            :: grazingP1,grazingP2,grazingP3,lvC1,lvC2
      real(rk), parameter :: secs_pr_day = 86400.0_rk

      ! Enter spatial loops (if any)
      _LOOP_BEGIN_

         ! Retrieve current (local) state variable values.
         _GET_(self%id_C1,C1)           ! predator1
         _GET_(self%id_C2,C2)           ! predator2
         _GET_(self%id_P1,P1)           ! prey1
         _GET_(self%id_P2,P2)           ! prey2
         _GET_(self%id_P3,P3)           ! prey3

         ! Retrieve current environmental conditions.

         ! Light acclimation formulation based on surface light intensity.

         ! Loss rate of phytoplankton to detritus depends on local light intensity.

         ! Define some intermediate quantities that will be reused multiple times.
         lvC1 = C1 * ( (self%aP1C1*P1+self%aP2C1*P2+self%aP3C1*P3) / (1.0+self%bP1C1*P1+self%bP2C1*P2+self%bP3C1*P3) - self%dC1 )
         lvC2 = C2 * ( (self%aP1C2*P1+self%aP2C2*P2+self%aP3C2*P3) / (1.0+self%bP1C2*P1+self%bP2C2*P2+self%bP3C2*P3) - self%dC2 )
         grazingP1 = - C1 * ( self%aP1C1*P1 / (1.0+self%bP1C1*P1+self%bP2C1*P2+self%bP3C1*P3) ) - C2 * ( self%aP1C2*P1 / (1.0+self%bP1C2*P1+self%bP2C2*P2+self%bP3C2*P3) )  
         grazingP2 = - C1 * ( self%aP2C1*P2 / (1.0+self%bP1C1*P1+self%bP2C1*P2+self%bP3C1*P3) ) - C2 * ( self%aP2C2*P2 / (1.0+self%bP1C2*P1+self%bP2C2*P2+self%bP3C2*P3) )  
         grazingP3 = - C1 * ( self%aP3C1*P3 / (1.0+self%bP1C1*P1+self%bP2C1*P2+self%bP3C1*P3) ) - C2 * ( self%aP3C2*P3 / (1.0+self%bP1C2*P1+self%bP2C2*P2+self%bP3C2*P3) )  
         _SET_ODE_(self%id_C1, lvC1)
         _SET_ODE_(self%id_C2, lvC2)
         _SET_ODE_(self%id_P1, grazingP1)
         _SET_ODE_(self%id_P2, grazingP2)
         _SET_ODE_(self%id_P3, grazingP3)
         _SET_DIAGNOSTIC_(self%id_lvC1, lvC1)
         _SET_DIAGNOSTIC_(self%id_lvC2, lvC2)
         _SET_DIAGNOSTIC_(self%id_grazingP1, grazingP1)
         _SET_DIAGNOSTIC_(self%id_grazingP2, grazingP2)
         _SET_DIAGNOSTIC_(self%id_grazingP3, grazingP3)
         
         ! Set temporal derivatives

         ! If an externally maintained ...

         ! Export diagnostic variables

      ! Leave spatial loops (if any)
   _LOOP_END_
   end subroutine do



end module fussmann_C1

!-----------------------------------------------------------------------
! Copyright Bolding & Bruggeman ApS - GNU Public License - www.gnu.org
!-----------------------------------------------------------------------
