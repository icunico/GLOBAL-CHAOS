#include "fabm_driver.h"

module ihamocc_preformed_tracer

   use fabm_types
   use ihamocc_shared

   implicit none

   private

   type, extends(type_base_model), public :: type_ihamocc_preformed_tracer
      type (type_state_variable_id)      :: id_oxygen, id_phosph, id_alkali, id_sco212
      type (type_diagnostic_variable_id) :: id_prefo2, id_prefpo4, id_prefalk, id_prefdic
   contains
      procedure :: initialize
      procedure :: do
   end type type_ihamocc_preformed_tracer

contains

   subroutine initialize(self, configunit)
      class (type_ihamocc_preformed_tracer), intent(inout), target :: self
      integer,                  intent(in)            :: configunit
      
      call self%register_state_dependency(self%id_oxygen, 'oxygen', 'kmol m-3', 'Dissolved oxygen')
      call self%register_state_dependency(self%id_alkali, 'alkali', 'kmol m-3', 'Alkalinity')
      call self%register_state_dependency(self%id_sco212, 'sco212', 'kmol m-3', 'Dissolved co2')
      call self%register_state_dependency(self%id_phosph, 'phosph', 'kmol m-3', 'dissolved phosphate')
      
      call self%register_diagnostic_variable(self%id_prefo2,     'prefo2',     'kmol m-3', 'preformed Dissolved oxygen')
      call self%register_diagnostic_variable(self%id_prefpo4,    'prefpo4',    'kmol m-3', 'preformed dissolved phosphate')
      call self%register_diagnostic_variable(self%id_prefalk,    'prefalk',    'kmol m-3', 'preformed Alkalinity')
      call self%register_diagnostic_variable(self%id_prefdic,    'prefdic',    'kmol m-3', 'preformed Dissolved co2')
   end subroutine
   
   subroutine do(self, _ARGUMENTS_DO_)
      class (type_ihamocc_preformed_tracer), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_

      real(rk) :: oxygen, alkali, sco212, phosph
      
      _LOOP_BEGIN_
         _GET_(self%id_oxygen, oxygen)
         _GET_(self%id_alkali, alkali)
         _GET_(self%id_sco212, sco212)
         _GET_(self%id_phosph, phosph)
         
         _SET_DIAGNOSTIC_(self%id_prefo2,oxygen)
         _SET_DIAGNOSTIC_(self%id_prefalk,alkali)
         _SET_DIAGNOSTIC_(self%id_prefdic,sco212)
         _SET_DIAGNOSTIC_(self%id_prefpo4,phosph)
      _LOOP_END_
   end subroutine do
end module ihamocc_preformed_tracer
