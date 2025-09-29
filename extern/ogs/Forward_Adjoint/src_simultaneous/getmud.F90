      subroutine getmud(sunz,mud)
      IMPLICIT NONE 
!  Computes average cosine for direct irradiance just below
!  sea surface
      real(8), intent(in)  :: sunz
      real(8), intent(out) :: mud
      real(8), parameter   :: refrac_idx = 1.341D0
      real(8)              :: rad 
      real(8)              :: rsza, sinszaw
      real(8)              :: szaw, rmudl
 
!  Compute average cosine for direct irradiance in the water 
!  column given solar zenith angle (in degrees) at surface.
      rad    = 180.0D0/dacos(-1.0D0) ! radians
      rsza = sunz/rad
      sinszaw = sin(rsza)/refrac_idx
      szaw = asin(sinszaw)
      mud  = cos(szaw)
!     rmudl = 1.0D0/cos(szaw)   !avg cosine direct (1 over)
!     rmud = min(rmudl,1.5D0)
!     rmud = max(rmud,0.0D0)
      
 
      return
      end
