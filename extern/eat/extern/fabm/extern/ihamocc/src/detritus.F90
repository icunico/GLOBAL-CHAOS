#include "fabm_driver.h"

module ihamocc_detritus

   use fabm_types
   use ihamocc_shared

   implicit none

   private

   type, extends(type_base_model), public :: type_ihamocc_detritus
      type (type_dependency_id) :: id_depth, id_ptho, id_dimmor, id_pommor, id_phymor, id_phosy, id_exud, id_graton, id_grawa, id_gratpoc, id_hi
      type (type_dependency_id) :: id_phyrem, id_satoxy, id_wdust_in, id_wopal_in, id_wcal_in, id_wpoc_in, id_wnos_in
      type (type_surface_dependency_id) :: id_kmle
      type (type_state_variable_id) :: id_phy, id_silica, id_nos, id_oxygen, id_sco212, id_phosph, id_ano3, id_alkali, id_calc, id_opal, id_iron, id_an2o
      type (type_state_variable_id) :: id_fdust, id_adust, id_det, id_doc, id_gasnit
      type (type_diagnostic_variable_id) :: id_wopal, id_wcal, id_wpoc, id_wdust, id_bkopal, id_rem, id_remin2o, id_delcar, id_rdnit1, id_remin, id_dustagg
      type (type_diagnostic_variable_id) :: id_wnos, id_aggregate, id_pocrem, id_docrem, id_delcar_part, id_delsil, id_dnit
      type (type_bottom_state_variable_id) :: id_det_bot, id_calc_bot, id_opal_bot, id_fdust_bot, id_alkali_bot
      type (type_bottom_diagnostic_variable_id) :: id_flux_opal
      logical  :: AGG, WLIN
      real(rk) :: remido, bkopal, ropal, rcalc, calmax, drempoc, dremopal, dremn2o, dremsul, relaxfe, nmldmin
      real(rk) :: claydens, safe, SinkExp, cellsink, alow1, alar1, Stick, wmin, wmax, wline, wpoc, wcal, wopal, FractDim, cellmass, shear
   contains
      procedure :: initialize
      procedure :: do
      procedure :: do_bottom
      procedure :: get_vertical_movement
      procedure :: check_state
   end type type_ihamocc_detritus
   
      real(rk), parameter :: rdnit0 = 0.8_rk*ro2ut ! moles nitrate lost for remineralisation of 1 mole P. Paulmier et al. 2009, Table 1 and equation 18. Note that their R_0=ro2ut-2*rnit.
      real(rk), parameter :: rdnit1 = 0.8_rk*ro2ut-rnit !Paulmier et al. 2009, Table 1 and equation 18. Note that their R_0=ro2ut-2*rnit.
      real(rk), parameter :: rdnit2 = 0.4_rk*ro2ut !Paulmier et al. 2009, Table 1 and equation 18. Note that their R_0=ro2ut-2*rnit.      
      real(rk), parameter :: rdn2o1 = 2._rk*ro2ut-2.5_rk*rnit !Paulmier et al. 2009, Table 1 and equation 18. Note that their R_0=ro2ut-2*rnit.      
      real(rk), parameter :: rdn2o2 = 2._rk*ro2ut-2._rk*rnit !Paulmier et al. 2009, Table 1 and equation 18. Note that their R_0=ro2ut-2*rnit.      
      
contains

   subroutine initialize(self, configunit)
      class (type_ihamocc_detritus), intent(inout), target :: self
      integer,                  intent(in)            :: configunit
    
      call self%get_parameter(self%remido,  'remido',   'd-1',         'DOM remineralization rate',                                 default=0.004_rk)
      call self%get_parameter(self%drempoc, 'drempoc',  'd-1',         'deep sea poc remineralisation rate',                        default=0.025_rk)
      call self%get_parameter(self%dremopal,'dremopal', 'd-1',         'deep sea opal remineralisation rate',                       default=0.003_rk)
      call self%get_parameter(self%dremn2o, 'dremn2o',  'd-1',         'deep sea n2o remineralisation rate',                        default=0.01_rk)
      call self%get_parameter(self%dremsul, 'dremsul',  'd-1',         'deep sea sulphate remineralisation rate',                   default=0.005_rk)
      call self%get_parameter(self%bkopal,  'bkopal',   'kmol Si m-3', 'half sat. constant for opal',                               default=5.e-6_rk) !i.e. 0.04 mmol P/m3
      call self%get_parameter(self%ropal,   'ropal',    '-',           'opal to organic phosphorous production ratio',              default=10.5_rk)
      call self%get_parameter(self%rcalc,   'rcalc',    '-',           'calcium carbonate to organic phosphorous production ratio', default=14._rk)
      call self%get_parameter(self%claydens,'claydens', 'kg m-3',      'clay (quartz) density',                                     default=2600._rk)
      call self%get_parameter(self%relaxfe, 'relaxfe',  '-',           'relaxfe',                                                   default=1.3699e-4_rk)
      call self%get_parameter(self%wpoc,    'wpoc',     'm d-1',       'poc sinking speed',                                         default=5._rk)
      call self%get_parameter(self%wcal,    'wcal',     'm d-1',       'calcium carbonate sinking speed',                           default=30._rk)
      call self%get_parameter(self%wopal,   'wopal',    'm d-1',       'opal sinking speed',                                        default=30._rk)
      
      call self%get_parameter(self%AGG,     'AGG',      '-',            'turn on aggregations',                                      default=.FALSE.)
      if (self%AGG) then
        call self%get_parameter(self%calmax,   'calmax',   '-',     'calmax',                                                                  default=0.2_rk)
        call self%get_parameter(self%FractDim, 'FractDim', '-',     'FractDim',                                                                default=1.62_rk)
        call self%get_parameter(self%cellmass, 'cellmass', 'nmol P','cellmass',                                                                default=0.012_rk/rnit)
        call self%get_parameter(self%nmldmin,  'nmldmin',  '-',     'minimum particle number in mixed layer',                                  default=0.1_rk)
        call self%get_parameter(self%safe,     'safe',     '-',     'safe',                                                                    default=1.e-6_rk)
        call self%get_parameter(self%SinkExp,  'SinkExp',  '-',     'SinkExp',                                                                 default=0.62_rk)
        call self%get_parameter(self%cellsink, 'cellsink', '',      'cellsink',                                                                default=1.4_rk)
        call self%get_parameter(self%alow1,    'alow1',    'cm',    'diameter of smallest particle',                                           default=0.002_rk)
        call self%get_parameter(self%alar1,    'alar1',    'cm',    'diameter of largest particle for size dependent aggregation and sinking', default=0.5_rk)
        call self%get_parameter(self%Stick,    'Stick',    '-',     'Stick',                                                                   default=0.15_rk)
        call self%get_parameter(self%shear,    'shear',    'd-1',   'shear in the mixed layer',                                                default=43200._rk)
        
        call self%register_state_variable(self%id_nos,   'nos',   'g-1', 'marine snow aggregates per g sea water', minimum=0.0_rk)
        call self%register_state_variable(self%id_adust, 'adust', 'g-1', 'dust aggregates per g sea water', minimum=0.0_rk)
        
        call self%register_diagnostic_variable(self%id_aggregate, 'aggregate', 'g-1 d-1',    'marine snow aggregation')
        call self%register_diagnostic_variable(self%id_dustagg,   'dustagg',   'kg m-3 d-1', 'dust particle aggregation')
        call self%register_diagnostic_variable(self%id_wnos,      'wnos',      'm d-1',      'sinking speed of particles')

        call self%register_dependency(self%id_wnos_in, 'wnos', 'm d-1', 'sinking speed of particles')
      endif
   
      call self%get_parameter(self%WLIN, 'WLIN', '-','second sinking scheme', default=.FALSE.)
      if (self%WLIN .and. .not. self%AGG) then
        call self%get_parameter(self%wmin,  'wmin',  'm d-1', 'minimum sinking speed',                default=1._rk)
        call self%get_parameter(self%wmax,  'wmax',  'm d-1', 'maximum sinking speed',                default=60._rk)
        call self%get_parameter(self%wline, 'wline', 'm d-1', 'constant describing incr. with depth', default=0.025_rk)
      endif
    
      call self%register_state_variable(self%id_det,    'det',    'kmol P m-3', 'detritus', minimum=0.0_rk)
      call self%add_to_aggregate_variable(standard_variables%total_carbon,     self%id_det, scale_factor=rcar * 1e6_rk)
      call self%add_to_aggregate_variable(standard_variables%total_nitrogen,   self%id_det, scale_factor=rnit * 1e6_rk)
      call self%add_to_aggregate_variable(standard_variables%total_phosphorus, self%id_det, scale_factor=1e6_rk)
      call self%add_to_aggregate_variable(standard_variables%total_iron,       self%id_det, scale_factor=riron * 1e9_rk)
      
      call self%register_diagnostic_variable(self%id_wpoc,        'wpoc',        'm d-1',            'poc sinking speed', missing_value=0.0_rk, source=source_do)
      call self%register_diagnostic_variable(self%id_wcal,        'wcal',        'm d-1',            'calcium carbonate sinking speed')
      call self%register_diagnostic_variable(self%id_wopal,       'wopal',       'm d-1',            'opal sinking speed')
      call self%register_diagnostic_variable(self%id_wdust,       'wdust',       'm d-1',            'dust sinking speed')
      call self%register_diagnostic_variable(self%id_bkopal,      'bkopal',      'kmol Si m-3',      'half sat. constant for opal') !for bromo.f90
      call self%register_diagnostic_variable(self%id_pocrem,      'pocrem',      'kmol P m-3 d-1',   'deep remineralization of POC') !for cisonew.f90 & natdic.f90
      call self%register_diagnostic_variable(self%id_docrem,      'docrem',      'kmol P m-3 d-1',   'deep remineralization of DOC') !for cisonew.f90 & natdic.f90
      call self%register_diagnostic_variable(self%id_remin2o,     'remin2o',     'kmol P m-3 d-1',   'remin2o') !for cisonew.f90 & natdic.f90
      call self%register_diagnostic_variable(self%id_remin,       'remin',       'kmol P m-3 d-1',   'remin') !for cisonew.f90 & natdic.f90
      call self%register_diagnostic_variable(self%id_delcar_part, 'delcar_part', '-',                'calcium fraction of export') !for cisonew.f90 
      call self%register_diagnostic_variable(self%id_delcar,      'delcar',      'kmol P m-3 d-1',   'delcar') !for natdic.f90
      call self%register_diagnostic_variable(self%id_delsil,      'delsil',      'kmol P m-3 d-1',   'delsil') !for natdic.f90
      call self%register_diagnostic_variable(self%id_rdnit1,      'rdnit1',      '-',                'rdnit1') !for natdic.f90
      call self%register_diagnostic_variable(self%id_dnit,        'dnit',        'kmol m-3 s-1',     'denitrification rate')

      call self%register_diagnostic_variable(self%id_flux_opal,      'flux_opal', 'kmol Si m-2 d-1', 'bottom flux opal', source=source_do_bottom, output=output_instantaneous)
      
      call self%register_state_dependency(self%id_doc,        'doc',        'kmol P m-3',  'dissolvecd organic carbon')
      call self%register_state_dependency(self%id_silica,     'silica',     'kmol Si m-3', 'Silicid acid (Si(OH)4)')
      call self%register_state_dependency(self%id_phosph,     'phosph',     'kmol P m-3',  'phosphate')
      call self%register_state_dependency(self%id_opal,       'opal',       'kmol Si m-3', 'Biogenic silica')
      call self%register_state_dependency(self%id_det_bot,    'det_bot',    'kmol P m-2',  'target bottom pool for detritus')
      call self%register_state_dependency(self%id_calc_bot,   'calc_bot',   'kmol C m-2',  'target bottom pool for calcium')
      call self%register_state_dependency(self%id_alkali_bot, 'alkali_bot', 'kmol m-2',    'target bottom pool for alkali')
      call self%register_state_dependency(self%id_opal_bot,   'opal_bot',   'kmol Si m-2', 'target bottom pool for opal')
      call self%register_state_dependency(self%id_fdust_bot,  'fdust_bot',  'kmol Fe m-2', 'target bottom pool for dust')
      
      call self%register_state_dependency(self%id_gasnit, 'gasnit', 'kmol N m-3', 'Gaseous nitrogen (N2)')
      call self%register_state_dependency(self%id_phy,    'phy',    'kmol P m-3', 'phytoplankton')
      call self%register_state_dependency(self%id_oxygen, 'oxygen', 'kmol O m-3', 'Dissolved oxygen')
      call self%register_state_dependency(self%id_sco212, 'sco212', 'kmol C m-3', 'Dissolved co2')
      call self%register_state_dependency(self%id_ano3,   'ano3',   'kmol N m-3', 'Dissolved nitrate')
      call self%register_state_dependency(self%id_an2o,   'an2o',   'kmol N m-3', 'laughing gas')
      call self%register_state_dependency(self%id_alkali, 'alkali', 'kmol m-3',   'Alkalinity')
      call self%register_state_dependency(self%id_calc,   'calc',   'kmol C m-3', 'Calcium carbonate')
      call self%register_state_dependency(self%id_iron,   'iron',   'kmol Fe m-3','dissolved iron')
      call self%register_state_dependency(self%id_fdust,  'fdust',  'kg m-3',     'non-aggregated dust deposition')
      
      call self%register_dependency(self%id_kmle,   'kmle',   'm',              'Mixed layer depth')
      call self%register_dependency(self%id_phymor, 'phymor', 'kmol P m-3 d-1', 'photosynthetic mortality rate')
      call self%register_dependency(self%id_phyrem, 'phyrem', 'kmol P m-3 d-1', 'photosynthetic remineralization rate')
      call self%register_dependency(self%id_pommor, 'pommor', 'kmol P m-3 d-1', 'zooplankton particulate export from mortality')
      call self%register_dependency(self%id_dimmor, 'dimmor', 'kmol P m-3 d-1', 'zooplankton dissolved inorganic export from mortality')
      call self%register_dependency(self%id_phosy,  'phosy',  'kmol P m-3 d-1', 'photosynthetic rate')
      call self%register_dependency(self%id_exud,   'exud',   'kmol P m-3 d-1', 'phytoplankton exudation rate')
      call self%register_dependency(self%id_graton, 'graton', 'kmol P m-3 d-1', 'zooplankton sloppy feeding inorganic release rate')
      call self%register_dependency(self%id_grawa,  'grawa',  'kmol P m-3 d-1', 'zooplankton assimilation rate')
      call self%register_dependency(self%id_gratpoc,'gratpoc','kmol P m-3 d-1', 'zooplankton sloppy feeding particulate release rate')
      call self%register_dependency(self%id_satoxy, 'satoxy', 'kmol O m-3',     'oxygen solubility')
      
      call self%register_dependency(self%id_depth, standard_variables%depth)
      call self%register_dependency(self%id_ptho, standard_variables%temperature)
      
      call self%register_dependency(self%id_wdust_in, 'wdust', 'm d-1', 'dust sinking speed')
      call self%register_dependency(self%id_wopal_in, 'wopal', 'm d-1', 'opal sinking speed')
      call self%register_dependency(self%id_wcal_in,  'wcal',  'm d-1', 'calcium carbonate sinking speed')
      call self%register_dependency(self%id_wpoc_in,  'wpoc',  'm d-1', 'poc sinking speed')
   end subroutine
   
   subroutine do(self, _ARGUMENTS_DO_)
      class (type_ihamocc_detritus), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_

      real(rk) :: phy, oxygen, det, doc, phymor, pomex, phosy, avsil, silica, bacfra, export, avmass, delsil, delcar, avnos, delcar_part
      real(rk) :: anosloss, nos, exud, graton, grawa, nos_roc, zmornos, pocrem, docrem, phyrem, pommor, dimmor, an2o_roc, dtr, zdis
      real(rk) :: dissopal, opal, iron, doc_roc, phosph_roc, ano3_roc, det_roc, ptho, depth, gratpoc, satoxy, ano3, an2o, temp
      real(rk) :: sco212_roc, alkali_roc, oxygen_roc, calc_roc, silica_roc, opal_roc, iron_roc, fdust, remin2o, gasnit_roc, remin_out
      real(rk) :: sterzo, remin, opalrem, aou, refra, dustd1, dustd2, dustd3, dustsink, kmle, pupper, plower, vsmall, snow, fshear, fsh
      real(rk) :: fse, eps, e1, e2, e3, e4, es1, es3, TopF, TMFac, TSFac, alar2, alar3, TopM, wmass, wnumb, alow2, alow3, sagg1, sagg2, sagg4
      real(rk) :: shear_agg, sett_agg, effsti, aggregate, dfirst, dshagg, dsett, dustagg, wpoc, wcal, wopal, wnos, dnit
      
      _LOOP_BEGIN_                                    
         _GET_(self%id_ptho, ptho)               
         _GET_(self%id_depth, depth)             
         _GET_(self%id_phy, phy)                 
         _GET_(self%id_oxygen, oxygen)           
         _GET_(self%id_det, det)                 
         _GET_(self%id_doc, doc)                 
         _GET_(self%id_phymor, phymor)           
         _GET_(self%id_phyrem, phyrem)           
         _GET_(self%id_pommor, pommor)           
         _GET_(self%id_dimmor, dimmor)
         _GET_(self%id_phosy, phosy)
         _GET_(self%id_silica, silica)
         _GET_(self%id_opal, opal)
         _GET_(self%id_gratpoc, gratpoc)
         _GET_(self%id_graton, graton)
         _GET_(self%id_iron, iron)
         _GET_(self%id_satoxy, satoxy)
         _GET_(self%id_ano3, ano3)
         _GET_(self%id_an2o, an2o)
         _GET_(self%id_fdust, fdust)
         
         _SET_DIAGNOSTIC_(self%id_bkopal, self%bkopal) !pass bkopal parameter to bromo module.
         _SET_DIAGNOSTIC_(self%id_rdnit1,  rdnit1) !for natdic.f90

         temp = min(40._rk,max(-3._rk,ptho))

         doc_roc     = 0.0_rk
         phosph_roc  = 0.0_rk
         ano3_roc    = 0.0_rk
         det_roc     = 0.0_rk
         sco212_roc  = 0.0_rk
         alkali_roc  = 0.0_rk
         oxygen_roc  = 0.0_rk
         calc_roc    = 0.0_rk
         silica_roc  = 0.0_rk
         opal_roc    = 0.0_rk
         iron_roc    = 0.0_rk
         gasnit_roc  = 0.0_rk
         an2o_roc    = 0.0_rk
         
         dnit = 0.0_rk
         delcar_part = 0.0_rk
         delcar = 0.0_rk
         delsil = 0.0_rk
         zdis = 0.01_rk / ((self%FractDim + 0.01_rk)*self%cellmass)
         if (depth<=100_rk) then ! in photic zone             
             bacfra = self%remido*doc
             export = pommor + gratpoc + phymor
             avsil = max(0.0_rk,silica)
             
             delcar_part = self%rcalc * self%bkopal/(avsil+self%bkopal)
             if (self%AGG) then ! if aggregations are turned on
                 delsil = MIN(self%ropal*phosy*avsil/(avsil+self%bkopal),0.5_rk*avsil)
                 delcar = self%rcalc*MIN(self%calmax*phosy,(phosy-delsil/self%ropal))
             else
                 delsil = MIN(self%ropal*export*avsil/(avsil+self%bkopal),0.5_rk*avsil)
                 delcar = delcar_part * export
             endif
             ! net remineralization rate
             dtr = bacfra+graton+dimmor
             
             dissopal = self%dremopal*opal
             
             if (self%AGG) then! if aggregations are turned on
                 _GET_(self%id_nos, nos)
                 _GET_(self%id_exud, exud)
                 _GET_(self%id_graton, graton)
                 _GET_(self%id_grawa, grawa)
                 avmass = det + phy
                 
                 nos_roc  = 0.0_rk
                 if (avmass > 0._rk) then
                     avnos = nos
                     anosloss = (phosy-exud-graton-grawa)*avnos/avmass
                     nos_roc = anosloss
                 endif
                 
                 zmornos = pommor * zdis * 1.e+6_rk
                 nos_roc = nos_roc + zmornos
             endif
             
             doc_roc    = doc_roc    - bacfra
             phosph_roc = phosph_roc + dtr
             ano3_roc   = ano3_roc   + dtr*rnit
             det_roc    = det_roc    + export
             sco212_roc = sco212_roc + rcar*dtr - delcar
             alkali_roc = alkali_roc - 2._rk*delcar - (rnit + 1._rk)*dtr
             oxygen_roc = oxygen_roc - dtr*ro2ut
             calc_roc   = calc_roc   + delcar
             silica_roc = silica_roc + dissopal - delsil
             opal_roc   = opal_roc   + delsil - dissopal
             iron_roc   = iron_roc   + dtr*riron
             
             pocrem = 0.0_rk
             remin2o = 0.0_rk
             remin_out = 0.0_rk
             docrem = bacfra
         else ! below photic zone
             if (oxygen > 5.0e-8_rk) then
                 pocrem = MIN(self%drempoc*det,0.33_rk*oxygen/ro2ut)
                 docrem = MIN(self%remido*doc,0.33_rk*oxygen/ro2ut)
             else
                 pocrem = 0._rk
                 docrem = 0._rk
             endif
             sterzo = pommor + dimmor
             remin = pocrem + docrem + phyrem
             remin2o = 0.0_rk
             opalrem = self%dremopal*0.1_rk*(temp+3._rk)*opal
             
             aou = satoxy-oxygen
             refra = 1._rk+3._rk*(0.5_rk+sign(0.5_rk,aou-1.97e-4_rk))
             
             if (self%AGG) then! if aggregations are turned on
                 _GET_(self%id_nos, nos)
                 avmass = det + phy
                 if (avmass > 0._rk) then
                     avnos = nos
                     nos_roc = -remin*avnos/avmass
                 endif
                 zmornos = sterzo * zdis * 1.e+6_rk
                 nos_roc = nos_roc + zmornos
             endif

             gasnit_roc = gasnit_roc - remin*1.e-4_rk*ro2ut*refra
             an2o_roc   = an2o_roc   + remin*1.e-4_rk*ro2ut*refra
             det_roc    = det_roc    - pocrem + phymor + sterzo + gratpoc
             doc_roc    = doc_roc    - docrem
             phosph_roc = phosph_roc + remin
             ano3_roc   = ano3_roc   + remin*rnit
             sco212_roc = sco212_roc + rcar*remin
             alkali_roc = alkali_roc - (rnit + 1._rk)*remin
             oxygen_roc = oxygen_roc - ro2ut*remin - remin*1.e-4_rk*ro2ut*refra*0.5_rk
             iron_roc   = iron_roc   + remin*riron
             opal_roc   = opal_roc   - opalrem
             silica_roc = silica_roc + opalrem
             
             if (oxygen < 5.0e-7_rk) then
                 remin   = 0.05_rk * self%drempoc * MIN(det,0.5_rk * ano3 / rdnit1)
                 remin_out = remin
                 remin2o = self%dremn2o * MIN(det,0.003_rk * an2o / rdn2o1)
                 pocrem = pocrem + remin !+ remin2o
                 if (self%AGG) then! if aggregations are turned on
                     _GET_(self%id_nos, nos)
                     avmass = det + phy
                     if (avmass > 0.) then
                         avnos = nos
                         nos_roc = nos_roc - (remin + remin2o)*avnos/avmass
                     endif
                 endif
                 alkali_roc = alkali_roc + (rdnit1-1._rk)*remin-remin2o
                 sco212_roc = sco212_roc + rcar*(remin+remin2o)
                 det_roc    = det_roc    - (remin+remin2o)
                 phosph_roc = phosph_roc + (remin+remin2o)
                 ano3_roc   = ano3_roc   - rdnit1*remin
                 gasnit_roc = gasnit_roc + rdnit2*remin+rdn2o2*remin2o
                 an2o_roc   = an2o_roc   - rdn2o1*remin2o
                 iron_roc   = iron_roc   + riron*(remin+remin2o)
                 dnit = rdnit0*remin
                 if (ano3 < 3.0e-6_rk) then
                     remin = self%dremsul*det
                     pocrem = pocrem + remin
                     if (self%AGG) then! if aggregations are turned on
                         _GET_(self%id_nos, nos)
                         avmass = det + phy
                         if (avmass > 0.) then
                             avnos = nos
                             nos_roc = nos_roc - remin*avnos/avmass
                         endif
                     endif
                     det_roc    = det_roc    - remin
                     alkali_roc = alkali_roc - (rnit+1._rk)*remin
                     sco212_roc = sco212_roc + rcar*remin
                     phosph_roc = phosph_roc + remin
                     ano3_roc   = ano3_roc   + rnit*remin
                     iron_roc   = iron_roc   + riron*remin
                 endif
             endif
         endif
         dustd1 = 0.0001_rk !cm = 1 um, boundary between clay and silt
         dustd2=dustd1*dustd1
         dustd3=dustd2*dustd1
         dustsink = (9.81_rk * 86400._rk / 18._rk * (self%claydens - 1025._rk) / 1.567_rk * 1000._rk * dustd2 * 1.e-4_rk) ! g * sec per day / 18.  !excess density / dyn. visc.          
         if (self%AGG) then! if aggregations are turned on
             _GET_(self%id_nos, nos)
             _GET_SURFACE_(self%id_kmle, kmle)
             
             pupper = self%safe/((self%FractDim+self%safe)*self%cellmass)
             plower = 1._rk/(1.1_rk*self%cellmass)
             vsmall = 1.e-9_rk
             
             avmass = det + phy
             snow = avmass*1.e6_rk
             if (avmass > 0._rk) then
                 fsh = 0.163_rk * self%shear
                 fse = 0.125_rk * 3.1415927_rk * self%cellsink * 100._rk
                 if (depth <= kmle) then
                     avnos = max(self%nmldmin,nos)
                     fshear = fsh
                 else
                     fshear = 0.0_rk
                 endif
                 avnos = max(snow*pupper,avnos)
                 avnos = min(snow*plower,avnos)
                 eps   = ((1._rk+ self%FractDim)*snow-avnos*self%cellmass)/(snow-avnos*self%cellmass)

                 if (abs(eps-3._rk) < 1.e-15_rk)                            eps = 3._rk + vsmall
                 if (abs(eps-4._rk) < 1.e-15_rk)                            eps = 4._rk + vsmall
                 if (abs(eps-3._rk-self%SinkExp) < 1.e-15_rk)               eps = 3._rk + vsmall + self%SinkExp
                 if (abs(eps-1._rk-self%SinkExp-self%FractDim) < 1.e-15_rk) eps = 1._rk + vsmall + self%SinkExp + self%FractDim
                 
                 e1 = 1._rk - eps
                 e2 = 2._rk - eps
                 e3 = 3._rk - eps
                 e4 = 4._rk - eps
                 es1 = e1 + self%SinkExp
                 es3 = e3 + self%SinkExp
                 TopF = (self%alar1/self%alow1)**e1
                 TMFac = (self%alar1/self%alow1)**self%FractDim
                 TSFac = (self%alar1/self%alow1)**self%SinkExp
                 alar2 = self%alar1 * self%alar1
                 alar3 = alar2 * self%alar1
                 alow2 = self%alow1 * self%alow1
                 alow3 = alow2 * self%alow1
                 TopM = TopF * TMFac
                 
                 wmass = self%cellsink * ( (self%FractDim+e1)/ (self%FractDim+es1) + TopM * TSFac * self%SinkExp / (self%FractDim+es1)) ! SINKING SPEED FOR THIS LAYER
                 wnumb = self%cellsink * (e1/es1 + TopF*TSFac*self%SinkExp/es1)
                 
                 ! AGGREGATION shear kernel:
                 sagg1 = (TopF-1._rk) * (TopF*alar3-alow3) * e1 / e4 + 3._rk * (TopF*self%alar1-self%alow1) * (TopF*alar2-alow2) * e1 * e1 / (e2*e3)
                 sagg2 = TopF*((alar3 + 3._rk * (alar2*self%alow1*e1/e2 + self%alar1*alow2*e1/e3) + alow3*e1/e4) - TopF*alar3*(1._rk+3._rk*(e1/e2 + e1/e3) + e1/e4))
                 sagg4 = TopF * TopF * 4._rk * alar3
                 shear_agg = (sagg1+sagg2+sagg4) * fshear
                 
                 ! AGGREGATION settlement kernel:
                 sagg1 = (TopF * TopF * alar2 * TSFac - alow2) * self%SinkExp / (es3 * e3 * (es3 + e1)) + alow2 * ((1._rk - TopF * TSFac) / (e3 * es1) - (1._rk - TopF) / (es3*e1))
                 sagg2 = TopF * e1 * (TSFac * ( alow2 - TopF * alar2) / e3 - (alow2 - TopF * alar2 * TSFac) / es3)
                 sett_agg =  (e1*e1*sagg1+sagg2) * fse

                 effsti = self%Stick * opal*1.e+6_rk/self%ropal/ ((opal * 1.e+6_rk / self%ropal) + snow)

                 aggregate = (shear_agg+sett_agg) * effsti * avnos * avnos
                                  
                 ! DUST shear kernel:
                 dfirst = dustd3 + 3._rk * dustd2 * self%alar1 + 3._rk * dustd1 * alar2 + alar3
                 dshagg = e1 * fsh * (dfirst * TopF / e1 - ((TopF-1._rk)/e1*dustd3 + 3._rk*(TopF*self%alar1-self%alow1)/e2*dustd2 + 3._rk*(TopF*alar2-alow2)/e3*dustd1 + (TopF*alar3-alow3)/e4))

                 ! DUST settlement kernel:
                 dsett = fse * dustd2 * ((e1+self%SinkExp*TopF*TSFac)/es1-dustsink/self%cellsink)

                 dustagg = effsti * avnos * fdust * (dshagg+dsett)
             else
                 wmass = self%cellsink
                 wnumb = 0._rk
                 aggregate = 0._rk
                 dustagg = 0._rk
             endif
             wpoc  = wmass
             wcal  = wmass
             wopal = wmass
             wnos  = wnumb
             nos_roc   = nos_roc - aggregate
             
             _SET_DIAGNOSTIC_(self%id_aggregate, aggregate)
             _SET_DIAGNOSTIC_(self%id_dustagg, dustagg)
             _SET_DIAGNOSTIC_(self%id_wnos, wnos)
             
             _ADD_SOURCE_(self%id_nos, nos_roc/dtbgc)
             _ADD_SOURCE_(self%id_adust, dustagg/dtbgc)
             _ADD_SOURCE_(self%id_fdust, -dustagg/dtbgc)
         else if (self%WLIN) then
             wpoc = min(self%wmin+self%wline*depth, self%wmax)
             wcal = self%wcal
             wopal= self%wopal
         else
             wpoc = self%wpoc
             wcal = self%wcal
             wopal= self%wopal
         endif
         _SET_DIAGNOSTIC_(self%id_wdust,       dustsink)
         _SET_DIAGNOSTIC_(self%id_wpoc,        wpoc)
         _SET_DIAGNOSTIC_(self%id_wopal,       wopal)
         _SET_DIAGNOSTIC_(self%id_wcal,        wcal)
         _SET_DIAGNOSTIC_(self%id_pocrem,      pocrem) !for cisonew.f90
         _SET_DIAGNOSTIC_(self%id_docrem,      docrem) !for cisonew.f90
         _SET_DIAGNOSTIC_(self%id_remin2o,     remin2o)!for cisonew.f90
         _SET_DIAGNOSTIC_(self%id_remin,       remin_out)!for natdic.f90
         _SET_DIAGNOSTIC_(self%id_delcar_part, delcar_part) !for cisonew.f90
         _SET_DIAGNOSTIC_(self%id_delcar,      delcar) !for natdic.f90
         _SET_DIAGNOSTIC_(self%id_delsil,      delsil)
         _SET_DIAGNOSTIC_(self%id_dnit,        dnit)
         
         _ADD_SOURCE_(self%id_calc,   calc_roc  /dtbgc)
         _ADD_SOURCE_(self%id_doc,    doc_roc   /dtbgc)
         _ADD_SOURCE_(self%id_oxygen, oxygen_roc/dtbgc)
         _ADD_SOURCE_(self%id_opal,   opal_roc  /dtbgc)
         _ADD_SOURCE_(self%id_silica, silica_roc/dtbgc)
         _ADD_SOURCE_(self%id_gasnit, gasnit_roc/dtbgc)
         _ADD_SOURCE_(self%id_an2o,   an2o_roc  /dtbgc)
         _ADD_SOURCE_(self%id_det,    det_roc   /dtbgc)
         _ADD_SOURCE_(self%id_alkali, alkali_roc/dtbgc)
         _ADD_SOURCE_(self%id_sco212, sco212_roc/dtbgc)
         _ADD_SOURCE_(self%id_phosph, phosph_roc/dtbgc)
         _ADD_SOURCE_(self%id_ano3,   ano3_roc  /dtbgc)
         _ADD_SOURCE_(self%id_iron,   iron_roc  /dtbgc)
      _LOOP_END_
    end subroutine do

    subroutine do_bottom(self, _ARGUMENTS_DO_BOTTOM_)
      class (type_ihamocc_detritus), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_BOTTOM_
      
      real(rk) :: wdust, wpoc, wopal, wcal, wnos, det, calc, opal, fdust, phy, nos, adust, det_bot_roc, fdust_bot_roc
      
      _BOTTOM_LOOP_BEGIN_
         _GET_(self%id_wdust_in,wdust)
         _GET_(self%id_wpoc_in, wpoc)
         _GET_(self%id_wopal_in,wopal)
         _GET_(self%id_wcal_in, wcal)
         
         _GET_(self%id_det,   det)
         _GET_(self%id_calc,  calc)
         _GET_(self%id_opal,  opal)
         _GET_(self%id_fdust, fdust)
         
         _ADD_BOTTOM_FLUX_(self%id_det,   -wpoc*det/dtbgc)
         _ADD_BOTTOM_FLUX_(self%id_calc,  -wcal*calc/dtbgc)
         _ADD_BOTTOM_FLUX_(self%id_opal,  -wopal*opal/dtbgc)
         _ADD_BOTTOM_FLUX_(self%id_fdust, -wdust*fdust/dtbgc)
         
         det_bot_roc = wpoc*det
         fdust_bot_roc = wdust*fdust
         
         if (self%AGG) then
             _GET_(self%id_wnos_in, wnos)
             _GET_(self%id_phy,     phy)
             _GET_(self%id_nos,     nos)
             _GET_(self%id_adust,   adust)
             
             det_bot_roc = det_bot_roc + wpoc*phy
             fdust_bot_roc = fdust_bot_roc + wpoc*adust
             
             _ADD_BOTTOM_FLUX_(self%id_phy,   -wpoc*phy/dtbgc)
             _ADD_BOTTOM_FLUX_(self%id_nos,   -wnos*nos/dtbgc)
             _ADD_BOTTOM_FLUX_(self%id_adust, -wpoc*adust/dtbgc)
         endif
         _ADD_BOTTOM_SOURCE_(self%id_calc_bot, wcal*calc/dtbgc)
         _ADD_BOTTOM_SOURCE_(self%id_alkali_bot, wcal*2._rk*calc/dtbgc)
         _ADD_BOTTOM_SOURCE_(self%id_opal_bot, wopal*opal/dtbgc)
         _ADD_BOTTOM_SOURCE_(self%id_det_bot,  det_bot_roc/dtbgc)
         _ADD_BOTTOM_SOURCE_(self%id_fdust_bot,fdust_bot_roc/dtbgc)
         _SET_BOTTOM_DIAGNOSTIC_(self%id_flux_opal, wopal*opal/dtbgc)
      _BOTTOM_LOOP_END_
   end subroutine do_bottom
    
   subroutine get_vertical_movement(self, _ARGUMENTS_GET_VERTICAL_MOVEMENT_)
      class (type_ihamocc_detritus), intent(in) :: self
      _DECLARE_ARGUMENTS_GET_VERTICAL_MOVEMENT_
      
      real(rk) :: wdust, wpoc, wopal, wcal, wnos
      
      _LOOP_BEGIN_
         _GET_(self%id_wdust_in,wdust)
         _GET_(self%id_wpoc_in, wpoc)
         _GET_(self%id_wopal_in,wopal)
         _GET_(self%id_wcal_in, wcal)
         
         _ADD_VERTICAL_VELOCITY_(self%id_det,  -wpoc/dtbgc)
         _ADD_VERTICAL_VELOCITY_(self%id_calc, -wcal/dtbgc)
         _ADD_VERTICAL_VELOCITY_(self%id_opal, -wopal/dtbgc)
         _ADD_VERTICAL_VELOCITY_(self%id_fdust,-wdust/dtbgc)
         if (self%AGG) then
             _GET_(self%id_wnos_in, wnos)
             _ADD_VERTICAL_VELOCITY_(self%id_phy,  -wpoc/dtbgc)
             _ADD_VERTICAL_VELOCITY_(self%id_nos,  -wnos/dtbgc)
             _ADD_VERTICAL_VELOCITY_(self%id_adust,-wpoc/dtbgc)
         endif
      _LOOP_END_
   end subroutine get_vertical_movement
    
   subroutine check_state(self, _ARGUMENTS_CHECK_STATE_)
      class (type_ihamocc_detritus), intent(in) :: self
      _DECLARE_ARGUMENTS_CHECK_STATE_

      real(rk) :: phy, det, nos, kmle, pupper, plower, avmass, snow, nos_new, depth
      
      _LOOP_BEGIN_
         if (self%AGG) then
             _GET_(self%id_nos, nos)
             _GET_SURFACE_(self%id_kmle, kmle)
             _GET_(self%id_phy, phy)                 
             _GET_(self%id_det, det)  
             _GET_(self%id_depth, depth)
             
             pupper = self%safe/((self%FractDim+self%safe)*self%cellmass)
             plower = 1._rk/(1.1_rk*self%cellmass)
             
             avmass = det + phy
             snow = avmass*1.e6_rk
             if (avmass > 0._rk) then
                 if (depth <= kmle) then
                     nos = max(self%nmldmin,nos)
                 endif  
                 nos = max(snow*pupper,nos)
                 nos = min(snow*plower,nos)
             else
                 nos = 0.0_rk
             endif
             _SET_(self%id_nos,nos)
         endif
      _LOOP_END_
    end subroutine check_state
    
end module ihamocc_detritus
