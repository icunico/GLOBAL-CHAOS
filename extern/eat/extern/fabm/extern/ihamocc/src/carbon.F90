#include "fabm_driver.h"

module ihamocc_carbon

   use fabm_types
   use ihamocc_shared

   implicit none

   private

   type, extends(type_base_model), public :: type_ihamocc_carbon
      type (type_dependency_id) :: id_psao, id_ptho, id_prho, id_prb, id_hi_in, id_pddpo
      type (type_surface_dependency_id) :: id_atco2, id_pfu10, id_psicomo, id_ppao
      type (type_state_variable_id) :: id_sco212, id_alkali, id_calc, id_phosph, id_silica
      type (type_diagnostic_variable_id) :: id_hi, id_co2star, id_co3, id_omegaA, id_omegaC, id_Kw, id_dissol
      type (type_surface_diagnostic_variable_id) ::  id_dicsat, id_co2fxd, id_co2fxu, id_pco2d, id_pco2m, id_kwco2sol, id_kwco2d, id_co2sold, id_co2solm
   contains
      procedure :: initialize
      procedure :: do_surface
      procedure :: do
   end type type_ihamocc_carbon
contains

   subroutine initialize(self, configunit)
      class (type_ihamocc_carbon), intent(inout), target :: self
      integer,                  intent(in)            :: configunit
      
      call self%register_state_variable(self%id_sco212, 'sco212', 'kmol m-3', 'Dissolved co2', minimum=0.0_rk)
      call self%add_to_aggregate_variable(standard_variables%total_carbon,     self%id_sco212, scale_factor=1e6_rk)
      call self%register_state_variable(self%id_alkali, 'alkali', 'kmol m-3', 'Alkalinity', minimum=0.0_rk)
      call self%register_state_variable(self%id_calc,   'calc',   'kmol m-3', 'Calcium carbonate', minimum=0.0_rk)
      call self%add_to_aggregate_variable(standard_variables%total_carbon,     self%id_calc, scale_factor=1e6_rk)
      
      call self%register_diagnostic_variable(self%id_Kw,       'Kw',       'mol kg-1',                  'Water dissociation product')
      call self%register_diagnostic_variable(self%id_hi,       'hi',       'mol kg-1',                  'Hydrogen ion concentration',missing_value=1.e-20_rk)
      call self%register_diagnostic_variable(self%id_co2star,  'co2star',  'mol kg-1',                  'Dissolved CO2 (CO2*)')
      call self%register_diagnostic_variable(self%id_co3,      'co3',      'kmol m-3',                  'Dissolved carbonate (CO3)')
      call self%register_diagnostic_variable(self%id_dicsat,   'dicsat',   'kmol m-3',                  'Saturated dic')
      call self%register_diagnostic_variable(self%id_omegaA,   'omegaA',   '-',                         'omegaA')
      call self%register_diagnostic_variable(self%id_omegaC,   'omegaC',   '-',                         'omegaC')
      call self%register_diagnostic_variable(self%id_co2fxd,   'co2fxd',   'kmol m-2 s-1',              'Downwards co2 surface flux')
      call self%register_diagnostic_variable(self%id_co2fxu,   'co2fxu',   'kmol m-2 s-1',              'Upwards co2 surface flux')
      call self%register_diagnostic_variable(self%id_pco2d,    'pco2d',    'microatm',                  'Dry air co2 pressure')
      call self%register_diagnostic_variable(self%id_pco2m,    'pco2m',    'microatm',                  'Moist air co2 pressure')
      call self%register_diagnostic_variable(self%id_kwco2sol, 'kwco2sol', 'm mol s-1 kg-1 microatm-1', 'kwco2sol')
      call self%register_diagnostic_variable(self%id_kwco2d,   'kwco2d',   'm s-1',                     'kwco2d')
      call self%register_diagnostic_variable(self%id_co2sold,  'co2sold',  'mol kg-1 atm-1',            'co2sold')
      call self%register_diagnostic_variable(self%id_co2solm,  'co2solm',  'mol kg-1 atm-1',            'co2solm')
      call self%register_diagnostic_variable(self%id_dissol,   'dissol',   'kmol m-3',                  'dissol')
      
      call self%register_dependency(self%id_ppao,         standard_variables%surface_air_pressure) ! surface air pressure in pascal
      call self%register_dependency(self%id_psao,         standard_variables%practical_salinity)
      call self%register_dependency(self%id_ptho,         standard_variables%temperature)
      call self%register_dependency(self%id_prho,         standard_variables%density)
      call self%register_dependency(self%id_prb,          standard_variables%pressure)
      call self%register_dependency(self%id_pddpo,        standard_variables%cell_thickness)
      call self%register_dependency(self%id_pfu10,        standard_variables%wind_speed)
      call self%register_dependency(self%id_psicomo,      standard_variables%ice_area_fraction)
      call self%register_dependency(self%id_atco2,        'atco2',  '-',        'surface air carbon dioxide mixing ratio') ! atmospheric co2 mixing ratio (i.e. partial presure = mixing ratio*SLP/P_0 [atm]) 
      call self%register_dependency(self%id_hi_in,        'hi',     'mol kg-1', 'Hydrogen ion concentration')
      call self%register_state_dependency(self%id_silica, 'silica', 'kmol m-3', 'Silicid acid (Si(OH)4)')
      call self%register_state_dependency(self%id_phosph, 'phosph', 'kmol m-3', 'Dissolved phosphate')
   end subroutine
   
   subroutine do_surface(self, _ARGUMENTS_DO_SURFACE_)
      class (type_ihamocc_carbon), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_SURFACE_

      real(rk) :: t, t2, t3, t4, tk, tk100, s, psao, ptho, prho, prb, sco212, alkali, silica, phosph, hi, cu, ac, K1, K2, pco2
      real(rk) :: pddpo, scco2, ppao, pfu10, kwco2, rpp0, fluxu, fluxd, ta, atco2, psicomo, rrho, tc, sit, pt, ah1, Ks1, Kw, Kf
      real(rk) :: Kh, Kb, Ksi, Khd, K3p, K1p, Kspc, K2p, Kspa, tc_sat, dicsat
      integer ::  niter
      
      _SURFACE_LOOP_BEGIN_
         _GET_(self%id_ptho, ptho)
         _GET_(self%id_psao, psao)
         _GET_(self%id_pddpo, pddpo)
         _GET_(self%id_prho, prho)
         _GET_(self%id_prb, prb)
         _GET_(self%id_sco212, sco212)
         _GET_(self%id_alkali, alkali)
         _GET_(self%id_hi_in, hi)
         _GET_(self%id_silica, silica)
         _GET_(self%id_phosph, phosph)
         _GET_SURFACE_(self%id_atco2, atco2)
         _GET_SURFACE_(self%id_ppao, ppao)
         _GET_SURFACE_(self%id_pfu10, pfu10)
         _GET_SURFACE_(self%id_psicomo, psicomo)

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
   
         tc   = sco212 / rrho  ! convert to mol/kg
         ta   = alkali / rrho
         sit  = silica / rrho
         pt   = phosph / rrho
         ah1  = hi
         niter = 20
         
         CALL CARCHM_KEQUI(t,s,prb,Kh,Khd,K1,K2,Kb,Kw,Ks1,Kf,Ksi,             &
                           K1p,K2p,K3p,Kspc,Kspa)
   
         CALL CARCHM_SOLVE(s,tc,ta,sit,pt,K1,K2,Kb,Kw,Ks1,Kf,Ksi,K1p,K2p,K3p, &
                           ah1,ac,niter)
   
         ! Determine CO2*, HCO3- and CO3-- concentrations (in mol/kg soln)
         cu = ( 2._rk * tc - ac ) / ( 2._rk + K1 / ah1 )
  
         pco2 = cu * 1.e6_rk / Kh ! Determine CO2 pressure and fugacity (in micoatm)   NOTE: equation below for pCO2 needs requires CO2 in mol/kg

         scco2 = 2116.8_rk - 136.25_rk*t + 4.7353_rk*t2 - 0.092307_rk*t3 + 0.0007555_rk *t4 ! Schmidt numbers according to Wanninkhof (2014), Table 1

         kwco2 = (1._rk-psicomo) * Xconvxa * pfu10**2._rk*(660._rk/scco2)**0.5_rk    ! Transfer (piston) velocity kw according to Wanninkhof (2014), in units of ms-1 

         ! Ratio P/P_0, where P is the local SLP and P_0 is standard pressure (1 atm). This is
         ! used in all surface flux calculations where atmospheric concentration is given as a
         ! mixing ratio (i.e. partial presure = mixing ratio*SLP/P_0 [atm])
         rpp0 = ppao/atm2pa

         fluxd=atco2*rpp0*kwco2*dtbgc*Kh*1.e-6_rk*rrho ! Kh is in mol/kg/atm. Multiply by rrho (g/cm^3) to get fluxes in kmol/m^2   
         fluxu=pco2      *kwco2*dtbgc*Kh*1.e-6_rk*rrho
         fluxu=min(fluxu,fluxd-(1.e-5_rk - sco212)*pddpo) !JT set limit for CO2 outgassing to avoid negative DIC concentration, set minimum DIC concentration to 1e-5 kmol/m3 

         ! Calculate saturation DIC concentration in mixed layer
         ta = alkali / rrho
         CALL carchm_solve_DICsat(s,atco2*rpp0,ta,sit,pt,Kh,K1,K2,Kb,Kw,Ks1,Kf, &
                                 Ksi,K1p,K2p,K3p,tc_sat,niter)
         dicsat = tc_sat * rrho ! convert mol/kg to kmlo/m^3         
         
         _ADD_SURFACE_FLUX_(self%id_sco212, (fluxd-fluxu)/dtbgc)  !Nic: divided by the time step to get instantaneous rate of change
         
         _SET_SURFACE_DIAGNOSTIC_(self%id_dicsat, dicsat) !NOTE: Nic: Implemented as surface diagnostic. If required further down the water column, a subroutine with a vertical loop will be implemented.
         _SET_SURFACE_DIAGNOSTIC_(self%id_co2fxd, fluxd/dtbgc) ! Save up- and downward components of carbon fluxes for output
         _SET_SURFACE_DIAGNOSTIC_(self%id_co2fxu, fluxu/dtbgc)
         _SET_SURFACE_DIAGNOSTIC_(self%id_pco2d, cu * 1.e6 / Khd) ! Save pco2 w.r.t. dry air for output
         _SET_SURFACE_DIAGNOSTIC_(self%id_pco2m, pco2) !pCO2 wrt moist air
         _SET_SURFACE_DIAGNOSTIC_(self%id_kwco2sol, kwco2*Kh*1e-6) ! Save product of piston velocity and solubility for output
         _SET_SURFACE_DIAGNOSTIC_(self%id_kwco2d, kwco2)
         _SET_SURFACE_DIAGNOSTIC_(self%id_co2sold, Khd)
         _SET_SURFACE_DIAGNOSTIC_(self%id_co2solm, Kh)
      _SURFACE_LOOP_END_
   end subroutine do_surface   
   
   subroutine do(self, _ARGUMENTS_DO_)
      class (type_ihamocc_carbon), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_

      real(rk) :: t, t2, t3, t4, tk, tk100, s, psao, ptho, prho, prb, sco212, silica, phosph, alkali, hi, cu, ac, K1, K2, cb, cc, co2star
      real(rk) :: tc, ta, sit, pt, ah1, c03, omega, OmegaA, OmegaC, calc, rrho, Ks1, Kw, co3, Kf, Kb, Ksi, K3p, Kspa, K1p, Kspc, Khd, Kh
      real(rk) :: supsat, undsa, K2p, dissol
      integer  :: niter
      
      _LOOP_BEGIN_
         _GET_(self%id_ptho, ptho)
         _GET_(self%id_psao, psao)
         _GET_(self%id_prho, prho)
         _GET_(self%id_prb, prb)
         _GET_(self%id_sco212, sco212)
         _GET_(self%id_alkali, alkali)
         _GET_(self%id_hi_in, hi)
         _GET_(self%id_silica, silica)
         _GET_(self%id_phosph, phosph)
         _GET_(self%id_calc, calc)

         ! Carbon chemistry: Calculate equilibrium constants and solve for [H+] and
         ! carbonate alkalinity (ac)
         t    = min(40._rk,max(-3._rk,ptho))
         s    = min(40._rk,max( 25._rk,psao))
         rrho = prho/1000.0_rk                ! seawater density [kg/m3]->[g/cm3]
         prb  = prb/10._rk  !convert from dbar to bar. ORIGINAL: ptiestu(i,j,k)*98060*1.027e-6_rk ! pressure in unit bars, 98060 = one
   
         tc   = sco212 / rrho  ! convert to mol/kg
         ta   = alkali / rrho
         sit  = silica / rrho
         pt   = phosph / rrho
         ah1  = hi
         niter = 20
   
         CALL CARCHM_KEQUI(t,s,prb,Kh,Khd,K1,K2,Kb,Kw,Ks1,Kf,Ksi,             &
                           K1p,K2p,K3p,Kspc,Kspa)
   
         CALL CARCHM_SOLVE(s,tc,ta,sit,pt,K1,K2,Kb,Kw,Ks1,Kf,Ksi,K1p,K2p,K3p, &
                           ah1,ac,niter)
   
         if(ah1 .gt. 0._rk) then
            _SET_DIAGNOSTIC_(self%id_hi, max(1.e-20_rk,ah1))
         else
            _SET_DIAGNOSTIC_(self%id_hi, hi)
         endif
         
         ! Determine CO2*, HCO3- and CO3-- concentrations (in mol/kg soln)
         cu = ( 2._rk * tc - ac ) / ( 2._rk + K1 / ah1 )
         cb = K1 * cu / ah1
         cc = K2 * cb / ah1
         co2star=cu
   
         ! Carbonate ion concentration, convert from mol/kg to kmol/m^3 
         co3  = cc * rrho 
   
         ! -----------------------------------------------------------------
         ! Deep ocean processes
         omega = ( calcon * s / 35._rk ) * cc           ! Determine Omega Calcite/Aragonite and dissolution of caco3 based on OmegaC:
         OmegaA = omega / Kspa                          !   omegaC=([CO3]*[Ca])/([CO3]sat*[Ca]sat)
         OmegaC = omega / Kspc                          !   Following Sarmiento and Gruber book, assumed that [Ca]=[Ca]sat
         supsat=co3-co3/OmegaC     !   Thus, [CO3]sat=[CO3]/OmegaC. 
         undsa=MAX(0._rk,-supsat)
         dissol=MIN(undsa,0.05_rk*calc)

         _SET_DIAGNOSTIC_(self%id_co2star, co2star)
         _SET_DIAGNOSTIC_(self%id_co3, co3)
         _SET_DIAGNOSTIC_(self%id_omegaA, OmegaA)
         _SET_DIAGNOSTIC_(self%id_omegaC, OmegaC)
         _SET_DIAGNOSTIC_(self%id_Kw, Kw)
         _SET_DIAGNOSTIC_(self%id_dissol, dissol)
         _ADD_SOURCE_(self%id_calc,-dissol/dtbgc)
         _ADD_SOURCE_(self%id_alkali,(2._rk*dissol)/dtbgc)
         _ADD_SOURCE_(self%id_sco212,dissol/dtbgc)
      _LOOP_END_
   end subroutine do
end module ihamocc_carbon
