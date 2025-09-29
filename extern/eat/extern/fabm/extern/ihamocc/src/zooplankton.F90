#include "fabm_driver.h"

module ihamocc_zooplankton

   use fabm_types
   use ihamocc_shared

   implicit none

   private

   type, extends(type_base_model), public :: type_ihamocc_zooplankton
      type (type_dependency_id) :: id_ptho, id_phytomi, id_depth
      type (type_state_variable_id) :: id_zoo, id_phy, id_silica, id_sco212, id_phosph, id_det, id_doc
      type (type_diagnostic_variable_id) :: id_gratpoc, id_pommor, id_dimmor, id_graton, id_grawa, id_domex
      real(rk) :: grami, grazra, bkzoo, epsher, zinges, spemor, gammaz, ecan
   contains
      procedure :: initialize
      procedure :: do
   end type type_ihamocc_zooplankton

contains

   subroutine initialize(self, configunit)
      class (type_ihamocc_zooplankton), intent(inout), target :: self
      integer,                  intent(in)            :: configunit
      
      call self%get_parameter(self%grami,  'grami',  'kmol P/m3','minimum concentration of zooplankton', default=1.0e-10_rk) 
      call self%get_parameter(self%grazra, 'grazra', '1/d','zooplankton grazing rate',                   default=1.2_rk) 
      call self%get_parameter(self%bkzoo,  'bkzoo',  'kmol P/m3','zooplankton half sat. constant',       default=8.e-8_rk) 
      call self%get_parameter(self%epsher, 'epsher', '-','ingestion fraction of grazing',                default=0.8_rk) 
      call self%get_parameter(self%zinges, 'zinges', '-','assimilation efficiency',                      default=0.6_rk)
      call self%get_parameter(self%spemor, 'spemor', '1/d','quadratic mortality constant',               default=3.e6_rk)
      call self%get_parameter(self%gammaz, 'gammaz', '1/d','excretion rate',                             default=0.06_rk)
      call self%get_parameter(self%ecan,   'ecan',   'fraction of mortality as PO_4',                    default=0.95_rk) 

      call self%register_state_variable(self%id_zoo, 'zoo', 'kmol/m^3', 'zooplankton', initial_value=5.0e-9_rk , minimum=self%grami)
      call self%add_to_aggregate_variable(standard_variables%total_carbon,     self%id_zoo, scale_factor=rcar * 1e6_rk)
      call self%add_to_aggregate_variable(standard_variables%total_nitrogen,   self%id_zoo, scale_factor=rnit * 1e6_rk)
      call self%add_to_aggregate_variable(standard_variables%total_phosphorus, self%id_zoo, scale_factor=1e6_rk)
      call self%add_to_aggregate_variable(standard_variables%total_iron,       self%id_zoo, scale_factor=riron * 1e9_rk)

      call self%register_dependency(self%id_depth,        standard_variables%depth)
      call self%register_dependency(self%id_ptho,      standard_variables%temperature)
      call self%register_dependency(self%id_phytomi,   'phytomi', 'kmol P/m3', 'minimum concentration of phytoplankton')
      
      call self%register_state_dependency(self%id_phy, 'phy',     'kmol/m^3',  'phytoplankton')
      call self%register_state_dependency(self%id_doc, 'doc',     'kmol/m^3',  'Dissolved organic carbon')

      call self%register_diagnostic_variable(self%id_pommor,  'pommor',  'kmol/m3/d', 'zooplankton particulate export from mortality')
      call self%register_diagnostic_variable(self%id_dimmor,  'dimmor',  'kmol/m3/d', 'zooplankton dissolved inorganic export from mortality')
      call self%register_diagnostic_variable(self%id_domex,   'domex',   'kmol/m3/d', 'zooplankton dissolved organic export')
      call self%register_diagnostic_variable(self%id_graton,  'graton',  'kmol/m3/d', 'zooplankton sloppy feeding inorganic release rate')
      call self%register_diagnostic_variable(self%id_grawa,   'grawa',   'kmol/m3/d', 'zooplankton assimilation rate')
      call self%register_diagnostic_variable(self%id_gratpoc, 'gratpoc', 'kmol/m3/d', 'zooplankton sloppy feeding particulate release rate')
   end subroutine
   
   subroutine do(self, _ARGUMENTS_DO_SURFACE_)
      class (type_ihamocc_zooplankton), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_
      
      real(rk) :: ptho, avphy, phytomi, avgra, grazing, graton, gratpoc, grawa, zoothresh, zoomor, excdoc, pommor, dimmor, domex, depth
      
      _LOOP_BEGIN_
         _GET_(self%id_ptho, ptho)
         _GET_(self%id_phy, avphy)
         _GET_(self%id_phytomi, phytomi)
         _GET_(self%id_zoo, avgra)
         _GET_(self%id_depth, depth)
         
         if (depth<=100_rk) then
             grazing = MAX(0.0_rk,avgra*self%grazra*(avphy-phytomi)/(avphy+self%bkzoo)) ! NIC: Changed from BLOM-iHAMOCC, now identical to formulation in Six and Maier-Reimer (1996)
             graton = self%epsher*(1._rk-self%zinges)*grazing
             gratpoc = (1._rk-self%epsher)*grazing
             grawa = self%epsher*self%zinges*grazing
         else
             grazing = 0.0_rk
             graton = 0.0_rk
             gratpoc = 0.0_rk
             grawa = 0.0_rk
         endif
         zoothresh = MAX(0._rk,avgra-2._rk*self%grami)
         zoomor = self%spemor*zoothresh*zoothresh           ! *10 compared to linear in tropics (tinka)
         excdoc = self%gammaz*zoothresh                     ! excretion of doc by zooplankton
         pommor = zoomor*(1._rk-self%ecan)
         dimmor = zoomor*self%ecan
         domex = excdoc
         
         _ADD_SOURCE_(self%id_zoo, (grawa-excdoc-zoomor)/dtbgc)
         _ADD_SOURCE_(self%id_doc, excdoc/dtbgc)
         _ADD_SOURCE_(self%id_phy, -grazing/dtbgc)
         
         _SET_DIAGNOSTIC_(self%id_pommor,  pommor)
         _SET_DIAGNOSTIC_(self%id_gratpoc, gratpoc)
         _SET_DIAGNOSTIC_(self%id_dimmor,  dimmor)
         _SET_DIAGNOSTIC_(self%id_domex,   domex)
         _SET_DIAGNOSTIC_(self%id_graton,  graton)
         _SET_DIAGNOSTIC_(self%id_grawa,   grawa)
      _LOOP_END_
   end subroutine do
end module ihamocc_zooplankton
