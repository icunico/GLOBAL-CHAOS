#include "fabm_driver.h"

module ihamocc_natdic

   use fabm_types
   use ihamocc_shared

   implicit none

   private

   type, extends(type_base_model), public :: type_ihamocc_natdic
      type (type_dependency_id) :: id_psao, id_ptho, id_prho, id_prb, id_pddpo, id_phosy, id_phyrem, id_dimmor
      type (type_dependency_id) :: id_graton, id_pocrem, id_docrem, id_rdnit1, id_dano3, id_nathi_in, id_remin2o
      type (type_dependency_id) :: id_remin, id_delcar, id_depth, id_wcal
      type (type_surface_dependency_id) :: id_pfu10, id_psicomo, id_ppao
      type (type_state_variable_id) :: id_natcalc, id_natsco212, id_natalkali, id_phosph, id_silica, id_oxygen, id_ano3
      type (type_diagnostic_variable_id) :: id_nathi, id_natco3, id_natomegaA, id_natomegaC
      type (type_surface_diagnostic_variable_id) :: id_natco2fxd, id_natco2fxu, id_natpco2d
      real(rk) :: atco2_nat
   contains
      procedure :: initialize
      procedure :: do_surface
      procedure :: do
      procedure :: get_vertical_movement
   end type type_ihamocc_natdic

contains

   subroutine initialize(self, configunit)
      class (type_ihamocc_natdic), intent(inout), target :: self
      integer,                  intent(in)            :: configunit
      
      call self%get_parameter(self%atco2_nat, 'atco2_nat', 'ppm','CMIP6 pre-industrial reference CO2 atm concentration', default=284.32_rk)
      
      ! Register state variables
      call self%register_state_variable(self%id_natcalc,   'natcalc',   'kmol m-3', 'Natural Calcium carbonate', minimum=0.0_rk)
      call self%register_state_variable(self%id_natsco212, 'natsco212', 'kmol m-3', 'Dissolved natural co2', minimum=0.0_rk)
      call self%register_state_variable(self%id_natalkali, 'natalkali', 'kmol m-3', 'Natural alkalinity', minimum=0.0_rk)
      
      ! Register diagnostic variables
      call self%register_diagnostic_variable(self%id_nathi,     'nathi',     'mol kg-1',     'Natural Hydrogen ion concentration', missing_value=1.e-20_rk)
      call self%register_diagnostic_variable(self%id_natco3,    'natco3',    'kmol m-3',     'Natural Dissolved carbonate (CO3)')
      call self%register_diagnostic_variable(self%id_natco2fxd, 'natco2fxd', 'kmol m-2 s-1', 'Natural Downwards co2 surface flux')
      call self%register_diagnostic_variable(self%id_natco2fxu, 'natco2fxu', 'kmol m-2 s-1', 'Natural Downwards co2 surface flux')
      call self%register_diagnostic_variable(self%id_natpco2d,  'natpco2d',  'microatm',     'Natural Dry air co2 pressure')
      call self%register_diagnostic_variable(self%id_natomegaA, 'natomegaA', '-',            'natomegaA')
      call self%register_diagnostic_variable(self%id_natomegaC, 'natomegaC', '-',            'natomegaC')

      call self%register_dependency(self%id_psao,      standard_variables%practical_salinity)
      call self%register_dependency(self%id_ptho,      standard_variables%temperature)
      call self%register_dependency(self%id_prho,      standard_variables%density)
      call self%register_dependency(self%id_prb,       standard_variables%pressure)
      call self%register_dependency(self%id_pddpo,     standard_variables%cell_thickness)
      call self%register_dependency(self%id_pfu10,     standard_variables%wind_speed)
      call self%register_dependency(self%id_psicomo,   standard_variables%ice_area_fraction)
      call self%register_dependency(self%id_ppao,      standard_variables%surface_air_pressure) ! surface air pressure in pascal
      call self%register_dependency(self%id_depth,     standard_variables%depth)
      call self%register_dependency(self%id_nathi_in,  'nathi',     'mol kg-1',     'Natural Hydrogen ion concentration')
      call self%register_dependency(self%id_remin2o,   'remin2o',   'kmol m-3 d-1', 'remin2o')
      call self%register_dependency(self%id_remin,     'remin',     'kmol m-3 d-1', 'remin')
      call self%register_dependency(self%id_phyrem,    'phyrem',    'kmol m-3 d-1', 'photosynthetic remineralization rate')
      call self%register_dependency(self%id_dimmor,    'dimmor',    'kmol m-3 d-1', 'zooplankton dissolved inorganic export from mortality')
      call self%register_dependency(self%id_phosy,     'phosy',     'kmol m-3 d-1', 'photosynthetic rate')
      call self%register_dependency(self%id_graton,    'graton',    'kmol m-3 d-1', 'zooplankton sloppy feeding inorganic release rate')
      call self%register_dependency(self%id_pocrem,    'pocrem',    'kmol m-3 d-1', 'deep remineralization of POC') 
      call self%register_dependency(self%id_docrem,    'docrem',    'kmol m-3 d-1', 'deep remineralization of DOC')
      call self%register_dependency(self%id_delcar,    'delcar',    'kmol m-3 d-1', 'delcar')
      call self%register_dependency(self%id_rdnit1,    'rdnit1',    '-',            'rdnit1') !for natdic.f90
      call self%register_dependency(self%id_dano3,     'dano3',     'kmol m-3 d-1', 'dano3') 
      call self%register_dependency(self%id_wcal,      'wcal',      'm d-1',        'calcium carbonate sinking speed')
      
      ! Register environmental dependencies
      call self%register_state_dependency(self%id_silica, 'silica', 'kmol m-3', 'Silicid acid (Si(OH)4)')
      call self%register_state_dependency(self%id_phosph, 'phosph', 'kmol m-3', 'Dissolved hosphate')
      call self%register_state_dependency(self%id_oxygen, 'oxygen', 'kmol m-3',  'Dissolved oxygen')
      call self%register_state_dependency(self%id_ano3,   'ano3',   'kmol m-3', 'Dissolved nitrate')
   end subroutine
   
   subroutine do_surface(self, _ARGUMENTS_DO_SURFACE_)
      class (type_ihamocc_natdic), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_SURFACE_

      real(rk) :: t, t2, t3, t4, tk, tk100, s, psao, ptho, prho, prb, natsco212, natalkali, silica, phosph, hi, cu, ac, K1, K2, pco2, scco2, ppao
      real(rk) :: pfu10, kwco2, rpp0, fluxu, fluxd, ta, nathi, psicomo, pddpo, rrho, tc, sit, pt, ah1, Ksi, Kw, Kf, Khd, Ks1, Kh, Kb, K2p, K3p
      real(rk) :: Kspa, Kspc, K1p, natcu, natcb, natcc, natco3, natpco2, natfluxd, natfluxu
      integer  :: niter
      
      _SURFACE_LOOP_BEGIN_
         _GET_(self%id_natsco212, natsco212)
         _GET_(self%id_natalkali, natalkali)
         _GET_(self%id_ptho, ptho)
         _GET_(self%id_psao, psao)
         _GET_(self%id_nathi_in, nathi)
         _GET_(self%id_silica, silica)
         _GET_SURFACE_(self%id_psicomo, psicomo)
         _GET_SURFACE_(self%id_pfu10, pfu10)
         _GET_SURFACE_(self%id_ppao, ppao)
         _GET_(self%id_prb, prb)
         _GET_(self%id_pddpo, pddpo)
         _GET_(self%id_prho, prho)
         _GET_(self%id_phosph, phosph)

         ! Carbon chemistry: Calculate equilibrium constants and solve for [H+] and
         ! carbonate alkalinity (ac)
         t    = min(40._rk,max(-3._rk,ptho))
         t2   = t**2
         t3   = t**3
         t4   = t**4
         tk   = t + tzero
         tk100= tk/100.0_rk
         s    = min(40._rk,max( 25._rk,psao))
         rrho = prho/1000.0_rk                ! seawater density [kg/m3]->[g/cm3]
         prb  = prb*10._rk  !convert from dbar to bar. ORIGINAL: ptiestu(i,j,k)*98060*1.027e-6_rk ! pressure in unit bars, 98060 = onem
   
         tc   = natsco212 / rrho  ! convert to mol/kg
         ta   = natalkali / rrho
         sit  = silica / rrho
         pt   = phosph / rrho
         ah1  = nathi
         niter = 20
   
         CALL CARCHM_KEQUI(t,s,prb,Kh,Khd,K1,K2,Kb,Kw,Ks1,Kf,Ksi,             &
                           K1p,K2p,K3p,Kspc,Kspa)
   
         CALL CARCHM_SOLVE(s,tc,ta,sit,pt,K1,K2,Kb,Kw,Ks1,Kf,Ksi,K1p,K2p,K3p, &
                           ah1,ac,niter)
            
         ! Determine natural CO2*, HCO3- and CO3-- concentrations (in mol/kg soln)
         natcu = ( 2._rk * tc - ac ) / ( 2._rk + K1 / ah1 )
      
         natpco2 = natcu * 1.e6_rk / Kh ! Determine CO2 pressure and fugacity (in micoatm)   NOTE: equation below for pCO2 needs requires CO2 in mol/kg

         scco2 = 2116.8_rk - 136.25_rk*t + 4.7353_rk*t2 - 0.092307_rk*t3 + 0.0007555_rk*t4 ! Schmidt numbers according to Wanninkhof (2014), Table 1

         kwco2 = (1._rk-psicomo) * Xconvxa * pfu10**2._rk*(660._rk/scco2)**0.5_rk    ! Transfer (piston) velocity kw according to Wanninkhof (2014), in units of ms-1 

         ! Ratio P/P_0, where P is the local SLP and P_0 is standard pressure (1 atm). This is
         ! used in all surface flux calculations where atmospheric concentration is given as a
         ! mixing ratio (i.e. partial presure = mixing ratio*SLP/P_0 [atm])
         rpp0 = ppao/atm2pa

         natfluxd=self%atco2_nat*rpp0*kwco2*dtbgc*Kh*1.e-6_rk*rrho ! Kh is in mol/kg/atm. Multiply by rrho (g/cm^3) to get fluxes in kmol/m^2   NOTE: originally multiplied by dtbgc (86400s/d). Removed as FABM rates-of-change has units s-1
         natfluxu=natpco2      *kwco2*dtbgc*Kh*1.e-6_rk*rrho
         natfluxu=min(natfluxu,natfluxd-(1.e-5_rk - natsco212)*pddpo) !JT set limit for CO2 outgassing to avoid negative DIC concentration, set minimum DIC concentration to 1e-5 kmol/m3 
         
         _ADD_SURFACE_FLUX_(self%id_natsco212, (natfluxd-natfluxu)/dtbgc)  !Nic: divided by the time step to get instantaneous rate of change
         
         _SET_SURFACE_DIAGNOSTIC_(self%id_natco2fxd, natfluxd/dtbgc) ! Save up- and downward components of carbon fluxes for output
         _SET_SURFACE_DIAGNOSTIC_(self%id_natco2fxu, natfluxu/dtbgc)
         _SET_SURFACE_DIAGNOSTIC_(self%id_natpco2d, natcu * 1.e6 / Khd) ! Save pco2 w.r.t. dry air for output
      _SURFACE_LOOP_END_
   end subroutine do_surface   
   
   subroutine do(self, _ARGUMENTS_DO_SURFACE_)
      class (type_ihamocc_natdic), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_
      
      real(rk) :: ptho, psao, prho, prb, natsco212, natalkali, nathi, silica, phosph, natcalc, t, s, rrho, tc, alkali, ta, sit, pt, ah1
      real(rk) :: Kw, Kb, K1, Ks1, Ksi, K2, Khd, Kf, Kh, K3p, Kspc, Kspa, K1p, K2p, ac, natcu, natcb, natcc, cc, natco3, natomega, natOmegaA
      real(rk) :: natOmegaC, natsupsat, natundsa, natdissol, natcalc_roc, natalkali_roc, natsco212_roc, delcar, docrem, phosy, graton
      real(rk) :: dimmor, phyrem, pocrem, remin2o, remin, rdnit1, dano3, oxygen, ano3, depth
      integer  :: niter
      
      _LOOP_BEGIN_
         _GET_(self%id_ptho, ptho)
         _GET_(self%id_psao, psao)
         _GET_(self%id_prho, prho)
         _GET_(self%id_prb, prb)
         _GET_(self%id_natsco212, natsco212)
         _GET_(self%id_natalkali, natalkali)
         _GET_(self%id_nathi_in, nathi)
         _GET_(self%id_silica, silica)
         _GET_(self%id_phosph, phosph)
         _GET_(self%id_natcalc, natcalc)

         ! Carbon chemistry: Calculate equilibrium constants and solve for [H+] and
         ! carbonate alkalinity (ac)
         t    = min(40._rk,max(-3._rk,ptho))
         s    = min(40._rk,max( 25._rk,psao))
         rrho = prho/1000.0_rk                ! seawater density [kg/m3]->[g/cm3]
         prb  = prb/10._rk  !convert from dbar to bar. ORIGINAL: ptiestu(i,j,k)*98060*1.027e-6_rk ! pressure in unit bars, 98060 = onem
         !
         tc   = natsco212 / rrho  ! convert to mol/kg
         ta   = natalkali / rrho
         sit  = silica / rrho
         pt   = phosph / rrho
         ah1  = nathi
         niter = 20
   
         CALL CARCHM_KEQUI(t,s,prb,Kh,Khd,K1,K2,Kb,Kw,Ks1,Kf,Ksi,             &
                           K1p,K2p,K3p,Kspc,Kspa)
   
         CALL CARCHM_SOLVE(s,tc,ta,sit,pt,K1,K2,Kb,Kw,Ks1,Kf,Ksi,K1p,K2p,K3p, &
                           ah1,ac,niter)
   
         if(ah1.gt.0._rk) then
            _SET_DIAGNOSTIC_(self%id_nathi, max(1.e-20_rk,ah1))
         else
            _SET_DIAGNOSTIC_(self%id_nathi, nathi)
         endif

         ! Determine CO2*, HCO3- and CO3-- concentrations (in mol/kg soln)
         natcu = ( 2._rk * tc - ac ) / ( 2._rk + K1 / ah1 )
         natcb = K1 * natcu / ah1
         natcc = K2 * natcb / ah1
         !co2star=cu
   
         ! Carbonate ion concentration, convert from mol/kg to kmol/m^3 
         natco3  = natcc * rrho 
   
         ! Deep ocean processes
         natomega = ( calcon * s / 35._rk ) * natcc           ! Determine Omega Calcite/Aragonite and dissolution of caco3 based on OmegaC:
         natOmegaA = natomega / Kspa                          !   omegaC=([CO3]*[Ca])/([CO3]sat*[Ca]sat)
         natOmegaC = natomega / Kspc                          !   Following Sarmiento and Gruber book, assumed that [Ca]=[Ca]sat
         natsupsat=natco3-natco3/natOmegaC     !   Thus, [CO3]sat=[CO3]/OmegaC. 
         natundsa=MAX(0._rk,-natsupsat)
         natdissol=MIN(natundsa,0.05_rk*natcalc)

         _SET_DIAGNOSTIC_(self%id_natomegaA, natOmegaA)
         _SET_DIAGNOSTIC_(self%id_natomegaC, natOmegaC)
         _SET_DIAGNOSTIC_(self%id_natco3, natco3*rrho)
         
         natcalc_roc   = - natdissol
         natalkali_roc = 2._rk*natdissol
         natsco212_roc = natdissol
         
         !bgc sources/sinks
         _GET_(self%id_delcar,delcar)
         _GET_(self%id_docrem,docrem)
         _GET_(self%id_phosy,phosy)
         _GET_(self%id_graton,graton)
         _GET_(self%id_dimmor,dimmor)
         _GET_(self%id_phyrem,phyrem)
         _GET_(self%id_pocrem,pocrem)
         _GET_(self%id_remin2o,remin2o)
         _GET_(self%id_remin,remin)
         _GET_(self%id_rdnit1,rdnit1)
         _GET_(self%id_dano3,dano3)
         _GET_(self%id_oxygen,oxygen)
         _GET_(self%id_ano3,ano3)
         _GET_(self%id_depth,depth)

         !rocs
         natsco212_roc = natsco212_roc - delcar + rcar*(pocrem - phosy + graton + dimmor + docrem + phyrem + remin2o)
         natcalc_roc   = natcalc_roc   + delcar
         natalkali_roc = natalkali_roc - 2._rk*delcar - (rnit+1._rk)*((pocrem-remin) - phosy + graton + dimmor + docrem + phyrem) - dano3
         
         if (depth>100_rk .and. oxygen < 5.0e-7_rk .and. ano3 >= 3.0e-6_rk) then
             natalkali_roc = natalkali_roc + (rdnit1-1._rk)*remin - remin2o
         endif
         
         _ADD_SOURCE_(self%id_natcalc,  natcalc_roc/dtbgc)
         _ADD_SOURCE_(self%id_natalkali,natalkali_roc/dtbgc)
         _ADD_SOURCE_(self%id_natsco212,natsco212_roc/dtbgc)
      _LOOP_END_
   end subroutine do

   subroutine get_vertical_movement(self, _ARGUMENTS_GET_VERTICAL_MOVEMENT_)
      class (type_ihamocc_natdic), intent(in) :: self
      _DECLARE_ARGUMENTS_GET_VERTICAL_MOVEMENT_
      
      real(rk) :: wcal
      
      _LOOP_BEGIN_
         _GET_(self%id_wcal,wcal)
         
         _ADD_VERTICAL_VELOCITY_(self%id_natcalc, -wcal/dtbgc)
      _LOOP_END_
   end subroutine get_vertical_movement
end module ihamocc_natdic
