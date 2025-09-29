#include "fabm_driver.h"

!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model 
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!BOP
!
! !ROUTINE: PelChem
!
! DESCRIPTION
!       This process describes the additional dynamics of dissolved
!       compounds in the watercolumn. Parameterized processes are:
!       - nitrification
!       - denitrification
!       - reoxidation of reduction equivalents
!       - dissolution of biogenic silica
!       This function also calls the carbonate system dynamics
!       (INCLUDE_PELCO2) and iron dynamics (INCLUDE_PELFE)
!       if activated
!
! !INTERFACE
 module bfm_BenthicLayer

   use fabm_types
   use ogs_bfm_shared
   use ogs_bfm_pelagic_base
!  
!
! !AUTHORS
!   Original version by P. Ruardij and M. Vichi
!
!
!
! !REVISION_HISTORY
!
! COPYING
!   
!   Copyright (C) 2015 BFM System Team (bfm_st@lists.cmcc.it)
!   Copyright (C) 2006 P. Ruardij, M. Vichi
!   (rua@nioz.nl, vichi@bo.ingv.it)
!
!   This program is free software; you can redistribute it and/or modify
!   it under the terms of the GNU General Public License as published by
!   the Free Software Foundation;
!   This program is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTEABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!   GNU General Public License for more details.
!
!EOP
!-------------------------------------------------------------------------!
!BOC
!
!
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Implicit typing is never allowed
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  IMPLICIT NONE

  private

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Local Variables
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  type, extends(type_ogs_bfm_pelagic_base), public :: type_BenthicLayer
      type (type_dependency_id) :: id_gdept_n
      type (type_diagnostic_variable_id) :: id_isBen
      real(rk) :: p_BenDepth
   contains
      procedure :: initialize
      procedure :: do_column 
  end type


contains

   subroutine initialize(self, configunit)
   class (type_BenthicLayer), intent(inout), target :: self
      integer,                      intent(in)            :: configunit

      call self%register_diagnostic_variable(self%id_isBen, 'isBen', '', 'boolean for pelagic', source=source_do_column)
! Parameters (described in subroutine initialize, below)
      call self%register_dependency(self%id_gdept_n, standard_variables%depth)
      call self%get_parameter(self%p_BenDepth, 'p_BenDepth',   '[m]',             'depth of benthic process activations',default=12000.0D0)
   end subroutine

   subroutine do_column(self, _ARGUMENTS_DO_COLUMN_)
      class (type_BenthicLayer), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_COLUMN_

      real(rk) :: gdept_n
      real(rk) :: isBen

      _DOWNWARD_LOOP_BEGIN_
         _GET_(self%id_gdept_n, gdept_n)
         IF( gdept_n < self%p_BenDepth ) THEN

            isBen = 0.0D0
         ELSE
            isBen = 1.0D0
         ENDIF
            _SET_DIAGNOSTIC_(self%id_isBen, isBen)
      _DOWNWARD_LOOP_END_
   end subroutine

end module

