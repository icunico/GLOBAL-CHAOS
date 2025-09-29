#include "fabm_driver.h"
    
module ihamocc_shared
   use fabm_types
   implicit none
   
   public 
   
   real(rk), parameter :: Xconvxa = 6.97e-07_rk!oxygen.f90 !carbon.f90 !cfc.f90 !dms.f90 !natdic.f90 !nitrogen.f90      ! Wanninkhof's a=0.251 converted from [cm hr-1]/[m s-1]^2 to [ms-1]/[m s-1]^2       NIC: from carchm.f90
   real(rk), parameter :: atm2pa = 101325.0_rk !oxygen.f90 !carbon.f90 !nitrogen.f90                        ! conversion factor from atmospheres to pascal
   real(rk), parameter :: tzero = 273.15_rk    !oxygen.f90 !carbon.f90 !bromo.f90 !cfc.f90                     ! absolute min temperature (*C) 
   
   ! mo_control_bgc parameters
   real(rk), parameter :: dtbgc = 86400.0_rk   ! time step length [sec].
   real(rk), parameter :: dtb = 1.0_rk         ! time step length [days].
   
   ! mo_chemcon parameters
   real(rk), parameter :: BOR1=0.000232_rk     !BORON CONCENTRATION IN SEA WATER IN G/KG PER O/OO CL (RILEY AND SKIRROW, 1965, P.250)
   real(rk), parameter :: BOR2=1./10.811_rk    !INVERSE OF ATOMIC WEIGHT OF BORON [G**-1] (USED TO CONVERT SPECIFIC TOTAL BORAT INTO CONCENTRATIONS)
   real(rk), parameter :: SALCHL=1./1.80655_rk !CONVERSION FACTOR SALINITY -> CHLORINITY (AFTER WOOSTER ET AL., 1969)
   real(rk), parameter :: CALCON=0.01028_rk    !carbon.f90 & natdic.f90    !SET MEAN TOTAL [CA++] IN SEAWATER (MOLES/KG) (SEE BROECKER A. PENG, 1982, P. 26; [CA++](MOLES/KG)=1.028E-2*(S/35.); Value taken from Sarmiento and Gruber, 2006, p. 365
   real(rk), parameter :: OXYCO=1./22414.4_rk  !oxygen.f90 & nitrogen.f90  !INVERS OF NORMAL MOLAL VOLUME OF AN IDEAL GAS [mol/ml] at 0C

   ! beleg_parm parameters
   !---------------------------------------------------------------
   ! extended redfield ratio declaration
   ! Note: stoichiometric ratios are based on Takahashi etal. (1985)
   ! P:N:C:-O2 + 1:16:122:172
    real(rk), parameter :: ro2ut=172._rk !detritus.f90 & phytoplankton
    real(rk), parameter :: rcar=122._rk  !phytoplankton.f90, detritus.f90, natdic.f90 & cisonew.f90
    real(rk), parameter :: rnit=16._rk   !detritus.f90, natdic.f90, nitrogen.f90 & bromo.f90
    real(rk), parameter :: riron= 5._rk*rcar*1.e-6_rk !phytoplankton.f90 & detritus.f90 ! fe to P ratio in organic matter
    
contains
    
    subroutine carchm_kequi(temp,saln,prb,Kh,Khd,K1,K2,Kb,Kw,Ks1,Kf,Ksi,K1p,K2p,K3p,Kspc,Kspa) !Calculate equilibrium constant for the carbonate system
!     *REAL*    *temp*    - potential temperature [degr C].
!     *REAL*    *saln*    - salinity [psu].
!     *REAL*    *prb*     - pressure [bar].
!     *REAL*    *Kh*      - equilibrium constant Kh  =  [CO2]/pCO2, moist air.
!     *REAL*    *Khd*     - equilibrium constant Kh  =  [CO2]/pCO2, dry air.
!     *REAL*    *K1*      - equilibrium constant K1  = [H][HCO3]/[H2CO3].
!     *REAL*    *K2*      - equilibrium constant K2  = [H][CO3]/[HCO3].
!     *REAL*    *Kb*      - equilibrium constant Kb  = [H][BO2]/[HBO2].
!     *REAL*    *Kw*      - equilibrium constant Kw  = [H][OH].
!     *REAL*    *Ks1*     - equilibrium constant Ks1 = [H][SO4]/[HSO4].
!     *REAL*    *Kf*      - equilibrium constant Kf  = [H][F]/[HF].
!     *REAL*    *Ksi*     - equilibrium constant Ksi = [H][SiO(OH)3]/[Si(OH)4].
!     *REAL*    *K1p*     - equilibrium constant K1p = [H][H2PO4]/[H3PO4].
!     *REAL*    *K2p*     - equilibrium constant K2p = [H][HPO4]/[H2PO4].
!     *REAL*    *K3p*     - equilibrium constant K3p = [H][PO4]/[HPO4].
!     *REAL*    *Kspc*    - equilibrium constant Kspc= [Ca2+]T [CO3]T.
!     *REAL*    *Kspa*    - equilibrium constant Kspa= [Ca2+]T [CO3]T.
      REAL(rk),    INTENT(IN)    :: temp,saln,prb
      REAL(rk),    INTENT(OUT)   :: Kh,Khd,K1,K2,Kb,Kw,Ks1,Kf,Ksi,K1p,K2p,K3p,Kspc,Kspa

      ! Local varibles
      INTEGER                    :: js
      REAL(rk)                   :: tk,tk100,invtk,dlogtk
      REAL(rk)                   :: s,is,is2,sqrtis,s15,s2,sqrts,scl
      REAL(rk)                   :: nKhwe74,deltav,deltak,zprb,zprb2
      REAL(rk)                   :: lnkpok0(11)
      
      ! Local parameters
      real(rk), parameter :: ac1= -162.8301_rk!Constants for CO2 solubility in mol/kg/atm from moist air at one atm total pressure. 
      real(rk), parameter :: ac2= 218.2968_rk !Table 6 in WEISS, R.F., NITROUS OXIDE SOLUBILITY IN WATER AND SEAWATER, Marine Chemistry, 8, 347-359, 1980
      real(rk), parameter :: ac3= 90.9241_rk  !
      real(rk), parameter :: ac4= -1.47696_rk !
      real(rk), parameter :: bc1= 0.025695_rk !
      real(rk), parameter :: bc2= -0.025225_rk!
      real(rk), parameter :: bc3= 0.0049867_rk!
      real(rk), parameter :: ad1= -60.2409_rk !Constants for CO2 solubility in mol/kg/atm for dry air at one atm total pressure. 
      real(rk), parameter :: ad2= 93.4517_rk  !Table 1 in WEISS, R.F., CARBON DIOXIDE IN WATER AND SEAWATER: THE SOLUBILITY OF A NON - IDEAL GAS, Marine Chemistry, 2, 203-215, 1974
      real(rk), parameter :: ad3= 23.3585_rk  !
      real(rk), parameter :: bd1= 0.023517_rk !
      real(rk), parameter :: bd2= -0.023656_rk!
      real(rk), parameter :: bd3= 0.0047036_rk!
      real(rk), parameter :: rgas = 83.131_rk !Gas constant, value as used by Millero (1995)
      REAL(rk), DIMENSION(11) :: a0, a1, a2, b0, b1, b2                                               !Constants needed for pressure correction of equilibrium constants. 
      DATA a0 /-25.5_rk, -15.82_rk, -29.48_rk, -25.60_rk, -18.03_rk, -9.78_rk, -48.76_rk, &           !F. Millero, Thermodynamics of the carbon dioxide system in the oceans, 
               -46._rk, -14.51_rk, -23.12_rk, -26.57_rk/                                              !Geochimica et Cosmochimica Acta, Vol. 59, No. 4, pp. 661-677, 1995
      DATA a1 /0.1271_rk, -0.0219_rk, 0.1622_rk, 0.2324_rk, 0.0466_rk, -0.0090_rk,     &              !
               0.5304_rk, 0.5304_rk, 0.1211_rk, 0.1758_rk, 0.2020_rk/                                 !
      DATA a2 /0.0_rk, 0.0_rk, 2.608e-3_rk, -3.6246e-3_rk, 0.316e-3_rk,             &                 !
              -0.942e-3_rk, 0.0_rk, 0.0_rk, -0.321e-3_rk, -2.647e-3_rk, -3.042e-3_rk/                 !
      DATA b0 /-3.08e-3_rk, 1.13e-3_rk, -2.84e-3_rk, -5.13e-3_rk, -4.53e-3_rk,      &                 !
               -3.91e-3_rk, -11.76e-3_rk, -11.76e-3_rk, -2.67e-3_rk, -5.15e-3_rk,   &                 !
               -4.08e-3_rk/                                                                           !
      DATA b1 /0.0877e-3_rk, -0.1475e-3_rk, 0.0_rk, 0.0794e-3_rk, 0.09e-3_rk,       &                 !
               0.054e-3_rk, 0.3692e-3_rk, 0.3692e-3_rk, 0.0427e-3_rk,            &                    !
               0.09e-3_rk, 0.0714e-3_rk/                                                              !
      DATA b2 /0.0_rk, 0.0_rk, 0.0_rk, 0.0_rk, 0.0_rk, 0.0_rk, 0.0_rk, 0.0_rk, 0.0_rk, 0.0_rk, 0.0_rk/!

      s = MAX(25._rk,saln)
      tk = temp + tzero
      tk100 = tk/100.0_rk
      invtk = 1.0_rk / tk
      dlogtk = log(tk)
      is = 19.924_rk * s / ( 1000._rk - 1.005_rk * s )
      is2 = is * is
      sqrtis = SQRT(is)
      s15    = s**1.5_rk
      s2     = s * s
      sqrts  = SQRT(s)
      scl    = s * salchl
      
      ! Kh = [CO2]/ p CO2
      ! Weiss (1974), refitted for moist air Weiss and Price (1980) [mol/kg/atm]
      nKhwe74 = ac1+ac2/tk100+ac3*log(tk100)+ac4*tk100**2+s*(bc1+bc2*tk100+bc3*tk100**2)
      Kh      = exp( nKhwe74 )
      ! Khd = [CO2]/ p CO2
      ! Weiss (1974) for dry air [mol/kg/atm]
      nKhwe74 = ad1+ad2/tk100+ad3*log(tk100)+s*(bd1+bd2*tk100+bd3*tk100**2)
      Khd     = exp( nKhwe74 )
      ! K1 = [H][HCO3]/[H2CO3]   ; K2 = [H][CO3]/[HCO3]
      ! Millero p.664 (1995) using Mehrbach et al. data on seawater scale
      K1 = 10**( -1.0_rk * ( 3670.7_rk * invtk - 62.008_rk + 9.7944_rk * dlogtk - 0.0118_rk * s + 0.000116_rk * s2 ) )
      K2 = 10**( -1.0_rk * ( 1394.7_rk * invtk + 4.777_rk - 0.0184_rk * s + 0.000118_rk * s2 ) )
      ! Kb = [H][BO2]/[HBO2] !
      ! Millero p.669 (1995) using DATA from Dickson (1990)
      Kb = exp( ( -8966.90_rk - 2890.53_rk  * sqrts - 77.942_rk  * s + 1.728_rk * s15 - 0.0996_rk * s2 ) * invtk +    &
                ( 148.0248_rk + 137.1942_rk * sqrts + 1.62142_rk * s ) +                                        &
                ( -24.4344_rk - 25.085_rk   * sqrts - 0.2474_rk  * s ) * dlogtk + 0.053105_rk * sqrts * tk )
      ! K1p = [H][H2PO4]/[H3PO4] ; K2p = [H][HPO4]/[H2PO4] ; K3p = [H][PO4]/[HPO4]
      ! DOE(1994) eq 7.2.20 with footnote using data from Millero (1974)
      K1p = exp( -4576.752_rk * invtk + 115.525_rk - 18.453_rk * dlogtk + ( -106.736_rk * invtk + 0.69171_rk ) *      &
                 sqrts + ( -0.65643_rk * invtk - 0.01844_rk ) * s )
      K2p = exp( -8814.715_rk * invtk + 172.0883_rk - 27.927_rk * dlogtk + ( -160.340_rk * invtk + 1.3566_rk ) *      &
                 sqrts + ( 0.37335_rk * invtk - 0.05778_rk ) *s );
      K3p = exp( -3070.75_rk * invtk - 18.141_rk + ( 17.27039_rk * invtk + 2.81197_rk ) * sqrts + ( -44.99486_rk *    &
                 invtk - 0.09984_rk ) * s );
      ! Ksi = [H][SiO(OH)3]/[Si(OH)4]
      ! Millero p.671 (1995) using data from Yao and Millero (1995)
      Ksi = exp( -8904.2_rk * invtk + 117.385_rk - 19.334_rk * dlogtk + ( -458.79_rk * invtk + 3.5913_rk ) * sqrtis   & 
             + ( 188.74_rk * invtk - 1.5998_rk) * is + ( -12.1652_rk * invtk + 0.07871_rk) * is2 +                 &
                 log(1.0_rk-0.001005_rk*s))
      ! Kw = [H][OH] 
      ! Millero p.670 (1995) using composite data
      Kw = exp( -13847.26_rk * invtk + 148.9652_rk - 23.6521_rk * dlogtk + ( 118.67_rk * invtk - 5.977_rk + 1.0495_rk *  &
                dlogtk ) * sqrts - 0.01615_rk * s)
      ! Ks = [H][SO4]/[HSO4]
      ! Dickson (1990, J. chem. Thermodynamics 22, 113)
      Ks1 = exp( -4276.1_rk * invtk + 141.328_rk - 23.093_rk * dlogtk + ( -13856._rk * invtk + 324.57_rk - 47.986_rk *   &
                 dlogtk ) * sqrtis + ( 35474._rk * invtk - 771.54_rk + 114.723_rk * dlogtk ) * is - 2698._rk *     &
                 invtk * is**1.5_rk + 1776._rk * invtk * is2 + log(1.0_rk - 0.001005_rk * s ) )
      ! Kf = [H][F]/[HF]
      ! Dickson and Riley (1979) -- change pH scale to total
      Kf = exp( 1590.2_rk * invtk - 12.641_rk + 1.525_rk * sqrtis + log( 1.0_rk - 0.001005_rk * s ) + log( 1.0_rk + (    &
                0.1400_rk / 96.062_rk ) * scl / Ks1 ) )
      ! Kspc (calcite)
      ! apparent solubility product of calcite : Kspc = [Ca2+]T [CO32-]T
      ! where $[]_T$ refers to the equilibrium total (free + complexed) ion concentration.
      !          Mucci 1983 mol/kg-soln
      Kspc = 10._rk**( -171.9065_rk - 0.077993_rk * tk + 2839.319_rk / tk + 71.595_rk * log10( tk ) + ( - 0.77712_rk +    &
                   0.0028426_rk * tk + 178.34_rk / tk ) * sqrts - 0.07711_rk * s + 0.0041249_rk * s15 );
      ! Kspa (aragonite)
      ! apparent solubility product of aragonite : Kspa = [Ca2+]T [CO32-]T
      ! where $[]_T$ refers to the equilibrium total (free + complexed) ion concentration.
      !          Mucci 1983 mol/kg-soln
      Kspa = 10._rk**( -171.945_rk - 0.077993_rk * tk + 2903.293_rk / tk  + 71.595_rk * log10( tk ) + ( -0.068393_rk +    &
                   0.0017276_rk * tk + 88.135_rk / tk ) * sqrts - 0.10018_rk * s + 0.0059415_rk * s15 );
      
      
      !---------------------- Pressure effect on Ks (Millero, 95) --------------------
      ! index: K1 1, K2 2, Kb 3, Kw 4, Ks 5, Kf 6, Kspc 7, Kspa 8, K1p 9, K2p 10, K3p 11
      DO js = 1,11
         deltav      = a0(js) + a1(js) * temp + a2(js) * temp * temp
         deltak      = b0(js) + b1(js) * temp + b2(js) * temp * temp
         zprb        = prb / ( rgas * tk )
         zprb2       = prb * zprb
         lnkpok0(js) = - ( deltav * zprb + 0.5_rk * deltak * zprb2 )
      ENDDO
      
      K1   = K1   * exp( lnkpok0(1)  )
      K2   = K2   * exp( lnkpok0(2)  )
      Kb   = Kb   * exp( lnkpok0(3)  )
      Kw   = Kw   * exp( lnkpok0(4)  )
      Ks1  = Ks1  * exp( lnkpok0(5)  )
      Kf   = Kf   * exp( lnkpok0(6)  )
      Kspc = Kspc * exp( lnkpok0(7)  )
      Kspa = Kspa * exp( lnkpok0(8)  )
      K1p  = K1p  * exp( lnkpok0(9)  )
      K2p  = K2p  * exp( lnkpok0(10) )
      K3p  = K3p  * exp( lnkpok0(11) )
   
   end subroutine carchm_kequi

   subroutine carchm_solve(saln,tc,ta,sit,pt,K1,K2,Kb,Kw,Ks1,Kf,Ksi,K1p,K2p,K3p,ah1,ac,niter) !Solve carbon chemistry.
!     *REAL*    *saln*    - salinity [psu].
!     *REAL*    *tc*      - total DIC concentraion [mol/kg].
!     *REAL*    *ta*      - total alkalinity [eq/kg].
!     *REAL*    *sit*     - silicate concentration [mol/kg].
!     *REAL*    *pt*      - phosphate concentration [mol/kg].
!     *REAL*    *K1*      - equilibrium constant K1  = [H][HCO3]/[H2CO3].
!     *REAL*    *K2*      - equilibrium constant K2  = [H][CO3]/[HCO3].
!     *REAL*    *Kb*      - equilibrium constant Kb  = [H][BO2]/[HBO2].
!     *REAL*    *Kw*      - equilibrium constant Kw  = [H][OH].
!     *REAL*    *Ks1*     - equilibrium constant Ks1 = [H][SO4]/[HSO4].
!     *REAL*    *Kf*      - equilibrium constant Kf  = [H][F]/[HF].
!     *REAL*    *Ksi*     - equilibrium constant Ksi = [H][SiO(OH)3]/[Si(OH)4].
!     *REAL*    *K1p*     - equilibrium constant K1p = [H][H2PO4]/[H3PO4].
!     *REAL*    *K2p*     - equilibrium constant K2p = [H][HPO4]/[H2PO4].
!     *REAL*    *K3p*     - equilibrium constant K3p = [H][PO4]/[HPO4].
!     *REAL*    *ah1*     - hydrogen ion concentration.
!     *REAL*    *ac*      - carbonate alkalinity.
!     *INTEGER* *niter*   - maximum number of iteration
      REAL(rk),    INTENT(IN)    :: saln,tc,ta,sit,pt
      REAL(rk),    INTENT(IN)    :: K1,K2,Kb,Kw,Ks1,Kf,Ksi,K1p,K2p,K3p
      REAL(rk),    INTENT(INOUT) :: ah1
      REAL(rk),    INTENT(OUT)   :: ac
      INTEGER, INTENT(IN)        :: niter
        
      ! Parameters to set accuracy of iteration 
      REAL(rk),    PARAMETER     :: eps=5.e-5_rk
      
      ! Local varibles
      INTEGER                    :: jit
      REAL(rk)                   :: s,scl,borat,sti,ft
      REAL(rk)                   :: hso4,hf,hsi,hpo4,ab,aw,ah2o,ah2,erel

      ! Calculate concentrations for borate, sulfate, and fluoride; see Dickson, A.G.,
      ! Sabine, C.L. and Christian, J.R. (Eds.) 2007. Guide to best practices 
      ! for ocean CO2 measurements. PICES Special Publication 3, chapter 5 p. 10
      s = MAX(25._rk,saln)
      scl = s * salchl
      borat = bor1 * scl * bor2           ! Uppstrom (1974)
      sti = 0.14_rk * scl / 96.062_rk     ! Morris & Riley (1966)
      ft = 0.000067_rk * scl / 18.9984_rk ! Riley (1965)


      iflag: DO jit = 1,niter
         hso4 = sti / ( 1._rk + Ks1 / ( ah1 / ( 1._rk + sti / Ks1 ) ) )
         hf   = 1._rk / ( 1._rk + Kf / ah1 )
         hsi  = 1._rk/ ( 1._rk + ah1 / Ksi )
         hpo4 = ( K1p * K2p * ( ah1 + 2._rk * K3p ) - ah1**3._rk ) /    & 
                ( ah1**3._rk + K1p * ah1**2._rk + K1p * K2p * ah1 + K1p * K2p * K3p )
         ab   = borat / ( 1._rk + ah1 / Kb )
         aw   = Kw / ah1 - ah1 / ( 1._rk + sti / Ks1 )
         ac   = ta + hso4 - sit * hsi - ab - aw + ft * hf - pt * hpo4
         ah2o = SQRT( ( tc - ac )**2._rk + 4._rk * ( ac * K2 / K1 ) * ( 2._rk * tc - ac ) )
         ah2  = 0.5_rk * K1 / ac *( ( tc - ac ) + ah2o )
         erel = ( ah2 - ah1 ) / ah2
         if (abs( erel ).ge.eps) then
            ah1 = ah2
         else
            exit iflag
         endif
      ENDDO iflag

   end subroutine carchm_solve

   subroutine carchm_solve_DICsat(saln,pco2,ta,sit,pt,Kh,K1,K2,Kb,Kw,Ks1,Kf,Ksi,K1p,K2p,K3p,tc_sat,niter) !Solve DICsat from TALK and pCO2.
!     *REAL*    *saln*    - salinity [psu].
!     *REAL*    *pco2*    - partial pressure of CO2 [ppm].
!     *REAL*    *ta*      - total alkalinity [eq/kg].
!     *REAL*    *sit*     - silicate concentration [mol/kg].
!     *REAL*    *pt*      - phosphate concentration [mol/kg].
!     *REAL*    *Kh*      - equilibrium constant K0  = [H2CO3]/pCO2.
!     *REAL*    *K1*      - equilibrium constant K1  = [H][HCO3]/[H2CO3].
!     *REAL*    *K2*      - equilibrium constant K2  = [H][CO3]/[HCO3].
!     *REAL*    *Kb*      - equilibrium constant Kb  = [H][BO2]/[HBO2].
!     *REAL*    *Kw*      - equilibrium constant Kw  = [H][OH].
!     *REAL*    *Ks1*     - equilibrium constant Ks1 = [H][SO4]/[HSO4].
!     *REAL*    *Kf*      - equilibrium constant Kf  = [H][F]/[HF].
!     *REAL*    *Ksi*     - equilibrium constant Ksi = [H][SiO(OH)3]/[Si(OH)4].
!     *REAL*    *K1p*     - equilibrium constant K1p = [H][H2PO4]/[H3PO4].
!     *REAL*    *K2p*     - equilibrium constant K2p = [H][HPO4]/[H2PO4].
!     *REAL*    *K3p*     - equilibrium constant K3p = [H][PO4]/[HPO4].
!     *REAL*    *tc_sat*  - saturated total DIC concentration [mol/kg].
!     *INTEGER* *niter*   - maximum number of iteration
      REAL(rk),    INTENT(IN)    :: saln,pco2,ta,sit,pt
      REAL(rk),    INTENT(IN)    :: Kh,K1,K2,Kb,Kw,Ks1,Kf,Ksi,K1p,K2p,K3p
      REAL(rk),    INTENT(OUT)   :: tc_sat
      INTEGER, INTENT(IN)        :: niter
     
      ! Parameters to set accuracy of iteration 
      REAL(rk),    PARAMETER     :: eps=5.e-5_rk
   
      ! Local varibles
      INTEGER                    :: jit
      REAL(rk)                   :: s,scl,borat,sti,ft
      REAL(rk)                   :: hso4,hf,hsi,hpo4,ab,aw,ah2o,ah2,erel
      REAL(rk)                   :: dic_h2co3,dic_hco3,dic_co3,ah1,ac
   
      ! Calculate concentrations for borate, sulfate, and fluoride; see Dickson, A.G.,
      ! Sabine, C.L. and Christian, J.R. (Eds.) 2007. Guide to best practices 
      ! for ocean CO2 measurements. PICES Special Publication 3, chapter 5 p. 10
      s = MAX(25._rk,saln)
      scl = s * salchl
      borat = bor1 * scl * bor2            ! Uppstrom (1974)
      sti = 0.14_rk * scl / 96.062_rk      ! Morris & Riley (1966)
      ft = 0.000067_rk * scl / 18.9984_rk  ! Riley (1965)
      ah1=1.e-8_rk
      dic_h2co3 = Kh * pco2 * 1.e-6_rk 
      
      iflag: DO jit = 1,niter
         hso4 = sti / ( 1._rk + Ks1 / ( ah1 / ( 1._rk + sti / Ks1 ) ) )
         hf   = 1._rk / ( 1._rk + Kf / ah1 )
         hsi  = 1._rk/ ( 1._rk + ah1 / Ksi )
         hpo4 = ( K1p * K2p * ( ah1 + 2._rk * K3p ) - ah1**3._rk ) /    & 
                ( ah1**3._rk + K1p * ah1**2._rk + K1p * K2p * ah1 + K1p * K2p * K3p )
         ab   = borat / ( 1._rk + ah1 / Kb )
         aw   = Kw / ah1 - ah1 / ( 1._rk + sti / Ks1 )
         ac   = ta + hso4 - sit * hsi - ab - aw + ft * hf - pt * hpo4
         ah2o = SQRT((K1*dic_h2co3)**2._rk + 4._rk*ac*2._rk*K1*k2*dic_h2co3) 
         ah2  = (K1*dic_h2co3 + ah2o)/(2._rk*ac)
         erel = ( ah2 - ah1 ) / ah2
         if (abs( erel ).ge.eps) then
            ah1 = ah2
         else
            exit iflag
         endif
      ENDDO iflag
      
      dic_hco3  = Kh * K1 *      pco2 * 1.e-6_rk / ah1
      dic_co3   = Kh * K1 * K2 * pco2 * 1.e-6_rk / ah1**2._rk
      tc_sat    = dic_h2co3 + dic_hco3 + dic_co3 
   end subroutine carchm_solve_DICsat
end module