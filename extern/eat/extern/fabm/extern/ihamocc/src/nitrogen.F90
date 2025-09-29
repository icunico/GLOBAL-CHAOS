#include "fabm_driver.h"

module ihamocc_nitrogen

   use fabm_types
   use ihamocc_shared

   implicit none

   private

   type, extends(type_base_model), public :: type_ihamocc_nitrogen
      type (type_surface_dependency_id)          :: id_ndep_in, id_psicomo, id_pfu10, id_ppao, id_atn2, id_kmle
      type (type_state_variable_id)              :: id_alkali, id_ano3, id_phosph, id_gasnit, id_an2o, id_oxygen 
      type (type_dependency_id)                  :: id_depth, id_ptho, id_psao
      type (type_surface_diagnostic_variable_id) :: id_ndep
      type (type_diagnostic_variable_id)         :: id_dano3
      real(rk) :: tf0, tf1, tf2, tff, bluefix, an0, an1, an2, an3, an4, an5, an6, al1, al2, al3, al4, bl1, bl2, bl3
   contains                                                                       
      procedure :: initialize                                                     
      procedure :: do_surface                                                    
      procedure :: do                                                             
   end type type_ihamocc_nitrogen                                                 
                                     
   real(rk), parameter :: an0 = -172.4965_rk
   real(rk), parameter :: an1 = 248.4262_rk
   real(rk), parameter :: an2 = 143.0738_rk
   real(rk), parameter :: an3 = -21.7120_rk
   real(rk), parameter :: an4 = -0.049781_rk
   real(rk), parameter :: an5 = 0.025018_rk
   real(rk), parameter :: an6 = -0.0034861_rk
   real(rk), parameter :: al1 = -165.8806_rk
   real(rk), parameter :: al2 = 222.8743_rk
   real(rk), parameter :: al3 = 92.0792_rk
   real(rk), parameter :: al4 = -1.48425_rk
   real(rk), parameter :: bl1 = -0.056235_rk
   real(rk), parameter :: bl2 = 0.031619_rk
   real(rk), parameter :: bl3 = -0.0048472_rk
   real(rk), parameter :: atn2 = 802000._rk   ! atmospheric nitrogen concentration in PPM
   real(rk), parameter :: atn2o = 3.e-7_rk   ! atmospheric laughing gas concentration in PPM
   
contains                                                                          

   subroutine initialize(self, configunit)
      class (type_ihamocc_nitrogen), intent(inout), target :: self
      integer,                  intent(in)            :: configunit
      
      call self%get_parameter(self%tf2,     'tf2',     '-',   'N-fix T dependency parameter 2',              default=-0.0042_rk)
      call self%get_parameter(self%tf1,     'tf1',     '-',   'N-fix T dependency parameter 1',              default=0.2253_rk)
      call self%get_parameter(self%tf0,     'tf0',     '-',   'N-fix T dependency parameter 0',              default=-2.7819_rk)
      call self%get_parameter(self%tff,     'tff',     '-',   'Trichodesmium max growth rate',               default=0.2395_rk)
      call self%get_parameter(self%bluefix, 'bluefix', 'd-1', 'nitrogen fixation rate',                      default=0.005_rk)
      
      ! Register state variables
      call self%register_state_variable(self%id_gasnit, 'gasnit', 'kmol m-3', 'Gaseous nitrogen (N2)', minimum=0.0_rk)
      call self%register_state_variable(self%id_an2o,   'an2o',   'kmol m-3', 'laughing gas', minimum=0.0_rk)
      call self%register_state_variable(self%id_ano3,   'ano3',   'kmol m-3', 'dissolved nitrate', minimum=0.0_rk)
      call self%add_to_aggregate_variable(standard_variables%total_nitrogen,   self%id_gasnit, scale_factor=1e6_rk)
      call self%add_to_aggregate_variable(standard_variables%total_nitrogen,   self%id_an2o, scale_factor=1e6_rk)
      call self%add_to_aggregate_variable(standard_variables%total_nitrogen,   self%id_ano3, scale_factor=1e6_rk)

      ! Register state dependencies
      call self%register_state_dependency(self%id_alkali, 'alkali', 'kmol m-3', 'Alkalinity')
      call self%register_state_dependency(self%id_phosph, 'phosph', 'kmol m-3', 'Dissolved hosphate')
      call self%register_state_dependency(self%id_oxygen, 'oxygen', 'kmol m-3', 'Dissolved oxygen')
      
      ! Register environmental dependencies
      call self%register_dependency(self%id_ndep_in, 'ndep_in', 'kmol m-2 s-1', 'nitrogen deposition flux')
      call self%register_dependency(self%id_kmle,    'kmle',    'm',            'Mixed layer depth')
      
      call self%register_dependency(self%id_psao,    standard_variables%practical_salinity)
      call self%register_dependency(self%id_psicomo, standard_variables%ice_area_fraction)
      call self%register_dependency(self%id_pfu10,   standard_variables%wind_speed)
      call self%register_dependency(self%id_ppao,    standard_variables%surface_air_pressure) ! surface air pressure in pascal
      call self%register_dependency(self%id_depth,   standard_variables%depth)
      call self%register_dependency(self%id_ptho,    standard_variables%temperature)
      
      ! Register diagnostics
      call self%register_diagnostic_variable(self%id_ndep,  'ndep',  'kmol m-2 s-1', 'nitrogen deposition flux')
      call self%register_diagnostic_variable(self%id_dano3, 'dano3', 'kmol m-3 d-1', 'nitrogen fixation rate')
   end subroutine
   
   subroutine do_surface(self, _ARGUMENTS_DO_SURFACE_)
      class (type_ihamocc_nitrogen), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_SURFACE_

      real(rk) :: ndep, ptho, psao, psicomo, ppao, pfu10, t, t2, t3, t4, tk, tk100, s, scn2, scn2o, ani, anisa, rs, satn2o, kwn2
      real(rk) :: gasnit, kwn2o, niflux, n2oflux, rpp0, an2o
      
      _SURFACE_LOOP_BEGIN_
         _GET_SURFACE_(self%id_ndep_in, ndep)
         _GET_(self%id_ptho, ptho)
         _GET_(self%id_psao, psao)
         _GET_SURFACE_(self%id_psicomo, psicomo)
         _GET_SURFACE_(self%id_ppao, ppao)
         _GET_SURFACE_(self%id_pfu10, pfu10)
         _GET_(self%id_gasnit, gasnit)
         _GET_(self%id_an2o, an2o)

         t = min(40._rk,max(-3._rk,ptho))
         t2   = t**2
         t3   = t**3
         t4   = t**4
         tk = t + tzero
         tk100 = tk/100.0_rk
         s = min(40._rk,max( 25._rk,psao))

         scn2  = 2304.8 - 162.75*t + 6.2557*t2 - 0.13129 *t3 + 0.0011255 *t4
         scn2o = 2356.2 - 166.38*t + 6.3952*t2 - 0.13422 *t3 + 0.0011506 *t4
         
         ! solubility of N2 (Weiss, R.F. 1970, Deep-Sea Res., 17, 721-735) for moist air at 1 atm; multiplication with oxyco converts to kmol/m^3/atm
         ani=an0+an1/tk100+an2*log(tk100)+an3*tk100+s*(an4+an5*tk100+an6*tk100**2._rk)
         anisa=exp(ani)*oxyco

         ! solubility of laughing gas  (Weiss and Price 1980, Marine Chemistry, 8, 347-359) 
         ! for moist air at 1 atm in kmol/m^3/atm
         rs=al1+al2/tk100+al3*log(tk100)+al4*tk100**2._rk+s*(bl1+bl2*tk100+bl3*tk100**2._rk)
         satn2o=exp(rs)
         
         kwn2  = (1._rk-psicomo) * Xconvxa * pfu10**2._rk*(660._rk/scn2)**0.5_rk
         kwn2o = (1._rk-psicomo) * Xconvxa * pfu10**2._rk*(660._rk/scn2o)**0.5_rk 
         rpp0 = ppao/atm2pa

         ! Surface flux of gaseous nitrogen (same piston velocity as for O2)
         niflux=kwn2*(gasnit-anisa*(atn2/802000._rk)*rpp0) 

         ! Surface flux of laughing gas (same piston velocity as for O2 and N2)
         n2oflux=kwn2o*(an2o-satn2o*atn2o*rpp0) 
         
         _ADD_SURFACE_FLUX_(self%id_an2o,       -n2oflux) ! NIC: positive flux indicates air -> water exchange; negative indicates water -> air exchange
         _ADD_SURFACE_FLUX_(self%id_gasnit,     -niflux)  ! NIC: positive flux indicates air -> water exchange; negative indicates water -> air exchange
         _ADD_SURFACE_FLUX_(self%id_ano3,       ndep)    ! NIC: positive flux indicates air -> water exchange; negative indicates water -> air exchange
         _ADD_SURFACE_FLUX_(self%id_alkali,     -ndep)   ! NIC: positive flux indicates air -> water exchange; negative indicates water -> air exchange
         _SET_SURFACE_DIAGNOSTIC_(self%id_ndep, ndep)
      _SURFACE_LOOP_END_
   end subroutine do_surface
   
   subroutine do(self, _ARGUMENTS_DO_)
      class (type_ihamocc_nitrogen), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_

      real(rk) :: ndep, kmle, depth, ano3, phosph, gasnit_roc, oxygen_roc, ano3_roc, alkali_roc, ptho, ttemp, nfixtfac

      _LOOP_BEGIN_
         _GET_SURFACE_(self%id_kmle,  kmle)
         _GET_(self%id_depth, depth)
         _GET_(self%id_ano3,  ano3)
         _GET_(self%id_phosph,phosph)
         gasnit_roc = 0.0_rk
         oxygen_roc = 0.0_rk 
         alkali_roc = 0.0_rk
         ano3_roc   = 0.0_rk
         if (depth <= kmle .and. ano3 < rnit*phosph) then
             _GET_(self%id_ptho, ptho)
             ttemp = min(40._rk,max(-3._rk,ptho))
                
             nfixtfac = MAX(0._rk,self%tf2*ttemp*ttemp + self%tf1*ttemp + self%tf0)/self%tff ! Temperature dependence of nitrogen fixation, Kriest and Oschlies 2015.

             ano3_roc   = ano3*(-self%bluefix*nfixtfac) + self%bluefix*nfixtfac*rnit*phosph
             gasnit_roc = -ano3_roc/2._rk
             oxygen_roc = -ano3_roc*1.25_rk ! Note: to fix one mole N2 requires: N2+H2O+y*O2 = 2* HNO3 <-> y=2.5 mole O2. I.e., to release one mole HNO3 = H+ + NO3- requires 1.25 mole O2
             alkali_roc = -ano3_roc ! Nitrogen fixation followed by remineralisation and nitrification decreases alkalinity by 1 mole per mole nitrogen fixed (Wolf-Gladrow et al. 2007)
         endif
         
         _ADD_SOURCE_(self%id_ano3,   ano3_roc/dtbgc) 
         _ADD_SOURCE_(self%id_gasnit, gasnit_roc/dtbgc)
         _ADD_SOURCE_(self%id_oxygen, oxygen_roc/dtbgc)
         _ADD_SOURCE_(self%id_alkali, alkali_roc/dtbgc)

         _SET_DIAGNOSTIC_(self%id_dano3, ano3_roc)
      _LOOP_END_
   end subroutine do
end module ihamocc_nitrogen
