MODULE adj_3stream

USE par_mod

PUBLIC:: test_3stream_adjoint, compute_3stream_adjoint, get_derivs, solve_direct

PRIVATE

integer, parameter:: d=1,s=2,u=3
double precision, parameter:: one = 1.0d0, zero = 0.0d0

CONTAINS

!solves both problems (direct and adjoint) and get the derivatives wrt all parameters
!INTENT(IN):
!n: the number of layers, so than there are n+1 horizons, from surface to bottom.
!z: depths of the horizons. z(1)=0 (is supposed to be), z(n+1)=H (bottom depth).
!m: the vertical grid size. This can be any points, normally >= 1 per layer
!zz: the vertical grid, of size m
!b, vd: parameters for Ed component; constant within a layer
!rd: a constant used for Ed; rd~1
!a, bb: parameters for all components; constant within a layer
!vs, vu: parameters for Es, Eu components; constant within a layer
!rs, ru: constants used for Es, Eu.
!W: the derivative of the functional G at the surface (must depend on surface data only!) wrt Eu.
!E: the (3,m) array that contains the solution to the direct problem: Ed, Es, Eu at the zz grid.
!INTENT(OUT): 
!derivs is the 4n+5 vector that contains derivatives of the functional wrt to: 
!  all a (1:n); all b (n+1:2n); all bb (2n+1:3n); all vd (3n+1:4n); vs,vu,rd,rs,ru
!sol_adj is the optinal (3,m) array to return the solution to the adjoint problem
!See the Overleaf document https://www.overleaf.com/project/5cbefd4a70921e1466457de2
subroutine get_derivs(m, zz, n, z, p, W, E, derivs, sol_adj)
use bioptimod_memory, only: nw,no
 implicit none
 integer, intent(in):: n, m
 double precision, dimension(m), intent(in):: zz !must be some grid from z(1)=0 to z(n+1)=bottom
 type(par), intent(inout):: p
 double precision, dimension(n+1), intent(in):: z
 double precision, intent(in):: W(nw)
 double precision, dimension(3,m,nw), intent(in):: E
 type(par), intent(out):: derivs !wrt a(n), b(n), b(n), vd(n), also wrt rd, rs, ru, vs, vu
 double precision, dimension(3,m,nw), intent(out), optional:: sol_adj
 double precision, dimension(3,m,nw):: lambda
 double precision, dimension(3):: lambdav,Ev
 double precision, dimension(3,3):: dA_dopt ! derivative of array over optical parameter
 double precision, dimension(m):: integrand
 double precision:: vd, vs, vu, rd, rs, ru
 double precision:: ac,bc,bbc
 integer:: L, k, iopt,wl
 vd=p%get1('vd'); vs=p%get1('vs'); vu=p%get1('vu')
 rd=p%get1('rd'); rs=p%get1('rs'); ru=p%get1('ru')
 derivs = par(p)   !clone the structure, the values to be replaced
 call derivs%reset !nullify the values
 call solve_adj(m, zz, n, z, p, W, lambda)
 if(present(sol_adj)) sol_adj = lambda
 !first let us get the derivatives wrt optical contituents
 do iopt=1,no
 do k=1,n
     integrand = 0.0d0
     do L=1,m                                            !loop over the grid
         if(zz(L).ge.z(k) .and. zz(L).lt.z(k+1)) then    !skip points above or below the layer
            do wl=1,nw
               ac  = p%a_coeff(wl,iopt)
               bc  = p%b_coeff(wl,iopt)
               bbc = p%bb_coeff(wl,iopt)

               dA_dopt(1,:) = [- (ac +bc)/vd     ,             0.0D0,           0.0D0]
               dA_dopt(2,:) = [ (bc - rd*bbc)/vd , - (ac +rs*bbc)/vs,       ru*bbc/vu]
               dA_dopt(3,:) = [- rd*bbc/vd       , -rs*bbc/vs       , (ac +ru*bbc)/vu]

               Ev      = [     E(d,L,wl),      E(s,L,wl),      E(u,L,wl)]
               lambdav = [lambda(d,L,wl), lambda(s,L,wl), lambda(u,L,wl)]
               integrand(L) = integrand(L) +  DOT_PRODUCT(lambdav,(MATMUL(dA_dopt,Ev))) !the expression under the integral
             enddo
         endif
     enddo !L: grid step within a layer
     derivs%opt_const(k,iopt) =  -trapz(zz, integrand)                     !integrating
!    call derivs%set('a', -trapz(zz, integrand), 1, k)                     !integrating
 enddo !k: the layer

 enddo !iopt: the optical constituent
end subroutine get_derivs

!solve the adjoint problem: build the matrix, solve the SLAE, get the coefficients and finally the solution lambda
!INTENT(IN):
!n: the number of layers, so than there are n+1 horizons, from surface to bottom.
!z: depths of the horizons. z(1)=0 (is supposed to be), z(n+1)=H (bottom depth).
!m: the vertical grid size. This can be any points, normally >= 1 per layer
!zz: the vertical grid, of size m
!b, vd: parameters for Ed component; constant within a layer
!rd: a constant used for Ed; rd~1
!a, bb: parameters for all components; constant within a layer
!vs, vu: parameters for Es, Eu components; constant within a layer
!rs, ru: constants used for Es, Eu.
!W: the derivative of the functional G at the surface (must depend on surface data only!) wrt Eu.
!INTENT(OUT): 
!lambda: the (3,m) array that contains the solution to the adjoint problem: lambda_d, lambda_s, lambda_u at the zz grid.
!See the Overleaf document https://www.overleaf.com/project/5cbefd4a70921e1466457de2
subroutine solve_adj(m, zz, n, z, p, W, lambda)
use bioptimod_memory, only: nw
use Tridiagonal, only: SLAE3diag
 implicit none
 integer, intent(in):: n, m
 double precision, dimension(m), intent(in):: zz !an argument, must be between z(1) and z(n+1)
 type(par), intent(inout):: p
 double precision, dimension(n+1), intent(in):: z
 double precision, dimension(3,m,nw),intent(out):: lambda
 double precision, intent(in):: W(nw)
 double precision, dimension(0:2*n-1):: diag, diagu, diagd, RHS
 double precision, dimension(n):: Ad, km, kp
 double precision, dimension(n,3):: Advec, kmvec, kpvec
 double precision, dimension(2*n):: psipm
 integer                         :: wl
 do wl = 1,nw
        !prepare the 3-diagonal matrix:
           call matrix_adj_3stream(n, z, wl, p, W(wl), &  !input
                   & diag, diagu, diagd, RHS,  &  !output: 3-diag matrix and RHS
                   & Ad, kp, km,               &  !eigenvalues
                   & Advec,kmvec,kpvec         &  !eigenvectors
                   & )                                
           call SLAE3diag(2*n,1,diagd, diag, diagu, RHS, psipm)                  !solve the SLAE for psi
           call get_lambda(m, zz, n,z, Ad, kp, km, Advec,kpvec,kmvec, psipm, lambda(:,:,wl)) !build the solution using the obtained psi
 end do
end subroutine solve_adj

!Prepare the 3-diagonal matrix for the adjoint problem
!INTENT(IN):
!n: the number of layers, so than there are n+1 horizons, from surface to the bottom.
!z: depths of the horizons. z(1)=0 (is supposed to be), z(n+1)=H (bottom depth).
!b, vd: parameters for Ed component; constant within a layer
!rd: a constant used for Ed; rd~1
!a, bb: parameters for all components; constant within a layer
!vs, vu: parameters for Es, Eu components; constant within a layer
!rs, ru: constants used for Es, Eu.
!W: the derivative of the functional G at the surface (must depend on surface data only!) wrt Eu.
!OUTPUT: 
!diag is the main diagonal, has the size 2n
!diagu is the diagonal above the main one, has the size 2n, but the last element is never used
!diagd is the diagonal below the main one, has the size 2n, but the first element is never used
!RHS is the vector in the right-hand side of the system. Is zero except for the first component, which is W
!See the Overleaf document https://www.overleaf.com/project/5cbefd4a70921e1466457de2
subroutine matrix_adj_3stream(n, z, wl, p, W,   &  !input
        & diag, diagu, diagd, RHS,          &  !output: 3-diag matrix and RHS
        & Ad, kp, km,                       &  !eigenvalues
        & Advec,kmvec,kpvec                 &  !eigenvectors
        & )                                
 implicit none
 integer, intent(in):: n,wl
 type(par), intent(inout):: p
 double precision, dimension(n+1), intent(in):: z
 double precision, intent(in):: W
 double precision, dimension(0:2*n-1), intent(out):: diag, diagu, diagd, RHS
 double precision, dimension(n), intent(out):: Ad, km, kp
 double precision, dimension(n,3), intent(out):: Advec, kmvec, kpvec
 double precision, dimension(n):: rkp,rkm, RRkp, RRkm
 double precision:: Fd, Bd, Bs, Bu, Cs, Cu, DD, dz
 double precision:: a,b,bb,vd,vs,vu,rd,rs,ru
 integer:: k
 vd=p%get1('vd'); vs=p%get1('vs'); vu=p%get1('vu')
 rd=p%get1('rd'); rs=p%get1('rs'); ru=p%get1('ru')
 do k=1,n
     a=p%get1('Ta',wl,k)
     b=p%get1('Tb',wl,k)
     bb=p%get1('Tbb',wl,k)
     !matrix
     Ad(k) = (a+b)/vd
     Fd = (b-rd*bb)/vd
     Bd = rd*bb/vd
     Cs = (a+rs*bb)/vs
     Cu = (a+ru*bb)/vu
     Bs = rs*bb/vs
     Bu = ru*bb/vu
     !eigenvalues:  Ad = Ad,
     DD = 0.5*(Cs+Cu+sqrt((Cs+Cu)*(Cs+Cu)-4*Bs*Bu))
     kp(k) = DD-Cu
     km(k) = DD-Cs
     !eigenvectors
     rkp(k) = Bs/DD
     rkm(k) = Bu/DD
     Advec(k,:) = [one, zero, zero]
     RRkp(k) = -(Fd*rkp(k) + Bd) / (Ad(k)-Cs+DD)
     RRkm(k) =  (Fd + Bd*rkm(k)) / (Ad(k)+Cu-DD)
!    write(*,*) 'NORM1', Ad(k)-Cs+DD
!    write(*,*) 'NORM2', Ad(k)+Cu-DD
     kmvec(k,:) = [RRkp(k), -rkp(k), one]
     kpvec(k,:) = [RRkm(k), one, -rkm(k)]
 enddo
 !the matrix, 3 diagonals
 diagd(0) = zero !not used
 diag(0) = -rkm(1)*exp(-kp(1)*z(1)) !xi0
 diagu(0) = exp(km(1)*z(1)) !nu0
 do k=1,n-1
         dz = z(k+1)-z(k)
         diag (2*k) = -(rkm(k)-rkm(k+1)) !xi 
         diagu(2*k) = 1-rkp(k+1)*rkm(k) !nu
         diagd(2*k) = -(rkm(k)-rkp(k))*exp(-km(k)*dz) !eps !bug corrected: minus sign in the exp
         diag (2*k-1) = -(rkp(k+1)-rkp(k))*exp(-km(k)*dz) !beta 
         diagu(2*k-1) = 1-rkm(k+1)*rkp(k+1) !gamma
         diagd(2*k-1) = -(1.-rkm(k)*rkm(k+1))*exp(kp(k)*dz) !alpha   
 enddo
 dz = z(n+1)-z(n)
 diagu(2*n-1) = zero !not used
 diagd(2*n-1) = exp(kp(n)*dz) !alpha
 diag (2*n-1) = -rkp(n)*exp(-km(n)*dz) !beta
! print*, 'dd', diagd
! print*, 'd ', diag
! print*, 'du', diagu
 !right-hand side
 RHS = 0.
 RHS(0) = W
 return
end subroutine matrix_adj_3stream

!get the solution lambda of the adjoint problem, given the coefficients psi^+, psi^- and the eigenvalues and eigenvectors
!INTENT(IN):
!n: the number of layers, so than there are n+1 horizons, from surface to the bottom.
!z: depths of the horizons. z(1)=0 (is supposed to be), z(n+1)=H (bottom depth).
!m: the vertical grid size. This can be any points, normally >= 1 per layer
!zz: the vertical grid, of size m
!Ad, kp, km: eigenvalues of the 3x3 matrix of the adjoint problem (A_d, k^-, k^+) for every layer
!Advec, kmvec, kpvec: the corresponding eigenvectors, of dimension 3, for every layer
!psipm: the coefficients \psi^+, \psi^- obtained by solving the 3-diag SLAE for the adjoint problem.
!OUTPUT: 
!lambda: three lambda functions (the solution to the adjoint problem) at the zz grid
subroutine get_lambda(m, zz, n, z,              &  !grid size, grid, depth level size, depth levels
           & Ad, kp, km,                        &  !eigenvalues
           & Advec,kpvec,kmvec,                 &  !eigenvectors
           & psipm,                             &  !the solution of the SLAE: both psi^+ and psi^- as odd and even items of the vector
           & lambda)                               !the solution of the adjoint problem (output)
 implicit none
 double precision, dimension(m), intent(in):: zz   !the vertical grid, must be between z(1) and z(n+1)
 integer, intent(in):: n, m
 double precision, dimension(n+1), intent(in):: z  !the layers, z(1)=0, z(n+1)=bottom
 double precision, dimension(n), intent(in):: Ad, km, kp             !the eigenvalues of the 3-matrix
 double precision, dimension(n,3), intent(in):: Advec, kmvec, kpvec  !the eigenvectors
 double precision, dimension(2*n), intent(in):: psipm                !the coefficients \psi obtained by solving a 3-diag SLAE
 double precision, dimension(3,m), intent(out):: lambda              !three lambda functions: \lambda_d, \lambda_s, \lambda_u
 double precision, dimension(n,4):: cPsi                             !aux: coefficients for the capital Psi
 double precision, dimension(n):: Psi                                !the capital Psi
 integer:: k, item
 call coef_Psi(n,z, Ad, kp, km, kmvec,kpvec, cPsi) !get the coefficients to calculate \Psi
 call get_captial_Psi(n,cPsi,psipm,Psi)            !calculate \Psi using the obtained coefficients 
 do k=1,n                                          !for each layer
         do item=1,m                               !find grid nodes that ara within the layer
                 if(zz(item)>=z(k) .and. zz(item)<=z(k+1)) then
                         lambda(:,item) =  Psi(k)*Advec(k,:)*exp(Ad(k) *(zz(item)-z(k))) + &
                                 & psipm(2*k)  *kmvec(k,:)*exp(-km(k)*(zz(item)-z(k))) + &
                                 & psipm(2*k-1)*kpvec(k,:)*exp(kp(k) *(zz(item)-z(k)))
                 endif
         enddo
 enddo
end subroutine get_lambda

!coefficients to get the captial Psi given the psi^- and psi^+
subroutine coef_Psi(n,z,                        &
           & Ad, kp, km,                        &  !eigenvalues
           & kmvec,kpvec,                       &  !eigenvectors
           & cPsi                               &  !the result: four coefficients for each layer
           )
 implicit none
 integer, intent(in):: n
 double precision, dimension(n+1), intent(in):: z
 double precision, dimension(n), intent(in):: Ad, km, kp
 double precision, dimension(n,3), intent(in):: kmvec, kpvec
 double precision, dimension(n,4), intent(out):: cPsi
 double precision:: dz, RRkp, RRkm
 integer:: k
 RRkp = kmvec(n,1)
 RRkm = kpvec(n,1)
 dz = z(n+1)-z(n)
 cPsi(n,1) = -RRkp*exp(-km(n)*dz)/exp(Ad(n)*dz) !coef for \psi_n^-
 cPsi(n,2) = -RRkm*exp( kp(n)*dz)/exp(Ad(n)*dz) !coef for \psi_n^+
 cPsi(n,3) = 0. !formal coef for \psi_{n+1}^+ : not used
 cPsi(n,4) = 0. !formal coef for \psi_{n+1}^+ : not used
 do k=n-1,1,-1
         dz = z(k+1)-z(k)
         RRkp = kmvec(k,1)
         RRkm = kpvec(k,1)
         cPsi(k,1) = -RRkp*exp(-km(k)*dz)/exp(Ad(k)*dz) !coef for \psi_k^-
         cPsi(k,2) = -RRkm*exp( kp(k)*dz)/exp(Ad(k)*dz) !coef for \psi_k^+
         RRkp = kmvec(k+1,1)
         RRkm = kpvec(k+1,1)
         cPsi(k,3) =  RRkp/exp(Ad(k)*dz) !coef for \psi_{k+1}^-
         cPsi(k,4) =  RRkm/exp(Ad(k)*dz) !coef for \psi_{k+1}^+
 enddo
end subroutine coef_Psi

!get the captial Psi given the psi^- and psi^+ and the coefficients
subroutine get_captial_Psi(n,cPsi,psipm,Psi)
 implicit none
 integer, intent(in):: n
 double precision, dimension(n,4), intent(in):: cPsi
 double precision, dimension(2*n), intent(in):: psipm !1,3,...=c_k^+, 2,4,6,...=c_k^-, k=1..n
 double precision, dimension(n), intent(out):: Psi
 integer:: k
 Psi(n) = psipm(2*n)*cPsi(n,1) + psipm(2*n-1)*cPsi(n,2)
 do k=n-1,1,-1
         Psi(k) = Psi(k+1) + psipm(2*k)*cPsi(k,1) + psipm(2*k-1)*cPsi(k,2) + psipm(2*k+2)*cPsi(k,3) + psipm(2*k+1)*cPsi(k,4)  !psipm has psi^+ as odd and psi_- as even, so psi^-(k)=psipm(2k), psi^+(k)=psipm(2k-1)
 enddo
end subroutine get_captial_Psi

!solve the direct problem
 subroutine solve_direct(m, zz, n, z, nlt, p, EdOASIM, EsOASIM, E)
 use Tridiagonal, only: SLAE3diag
 implicit none
 integer, intent(in):: n,m                                             !n of layers and size of the vertical grid
 integer, intent(in):: nlt                                             !n of wavelenghts to be considered
 double precision, dimension(m), intent(in):: zz                       !must be some grid from z(1)=0 to z(n+1)=bottom
 double precision, dimension(n+1), intent(in):: z                      !layer boundaries (depth levels); z(1)=0 (must be), z(n+1) = bottom
 type(par), intent(inout):: p                                             !input depth-dependent data
 double precision, dimension(nlt),intent(in) :: EdOASIM, EsOASIM       !boundary values
 double precision, dimension(3,m,nlt), intent(out):: E                 !the 3-stream solution on the zz grid
 double precision, dimension(n,nlt):: cd, x, y, kp, km, rkp, rkm, ekp, ekm !stuff used in the formulae, see the Overleaf doc https://www.overleaf.com/project/5cbefd4a70921e1466457de2
 double precision, dimension(n,2):: kmvec, kpvec                   !eigenvectors
 double precision, dimension(n+1,nlt):: Ed                             !solution of the first equation (which is independent of the other two) on the z grid
 double precision, dimension(nlt):: NORM
 double precision :: dz, dz1
 double precision, dimension(nlt) :: Fd, Bd, Cs, Cu, Bs, Bu, DD        !aux variables
 double precision, dimension(0:2*n-2,nlt):: diag, diagd, diagu, RHS    !the matrix, the right-hand side
 double precision, dimension(0:2*n-1,nlt):: cpm                        !the solution c_+ (0,2,4,...), c_- (1,3,5,...), up to cpm(2n-1)=c^-_n=0
 double precision, dimension(nlt):: ak,bk,bbk
 double precision:: vd,vs,vu,rd,rs,ru
 integer:: k, item,wl                                                  !counters
 logical :: isSingular(n,nlt)

 isSingular(:,:) = .False.

 vd=p%vd; vs=p%vs; vu=p%vu
 rd=p%rd; rs=p%rs; ru=p%ru
 Ed(1,:) = EdOASIM(:)                                                  !boundary value: Ed at the surface is just given
 do wl=1,nlt 
 do k=1,n                                                              !for each layer
! Ta(layern,wln)
     ak=p%Ta(k,wl)
     bk=p%Tb(k,wl)
     bbk=p%Tbb(k,wl)
!    write(*,*) 'ak', ak
!    write(*,*) 'bk', bk
!    write(*,*) 'bbk', bbk

!     call p%get('Ta',ak,layern=k)
!     call p%get('Tb',bk,layern=k)
!     call p%get('Tbb',bbk,layern=k)
         !matrix
         dz = z(k+1)-z(k)                                              !layer thickness
         cd(k,wl) = (ak(wl)+bk(wl))/vd                             !coefficient of the separate equation
         Ed(k+1,wl) = Ed(k,wl)*exp(-cd(k,wl)*dz)                          !solution of the separate equation
         Fd(wl) = (bk(wl)-rd*bbk(wl))/vd                           !components of the auxiliary matrix for the partial solution
         Bd(wl) = rd*bbk(wl)/vd
         Cs(wl) = (ak(wl)+rs*bbk(wl))/vs
         Cu(wl) = (ak(wl)+ru*bbk(wl))/vu
         Bs(wl) = rs*bbk(wl)/vs
         Bu(wl) = ru*bbk(wl)/vu
         !aux variables x,y
         NORM(wl) = ((cd(k,wl)-Cs(wl))*(cd(k,wl)+Cu(wl))+Bs(wl)*Bu(wl))
!        write(*,*) 'NORM', NORM(wl)
         if (NORM(wl) .LT. 1.0D-2) then
                 WRITE(*,*) 'SINGULIRTY'
                 isSingular(k,wl) = .TRUE.
         endif  

         if (isSingular(k,wl)) then
             x(k,wl) = Fd(wl)      !explicit solution of the aux SLAE
             y(k,wl) = -Bd(wl)
         else
             x(k,wl) = (-Fd(wl)*(Cu(wl)+cd(k,wl)) -Bd(wl)*Bu(wl))  / NORM(wl)       !explicit solution of the aux SLAE
             y(k,wl) = (-Fd(wl)*Bs(wl) +Bd(wl)*(-Cs(wl)+cd(k,wl))) / NORM(wl)
         endif 
         !eigenvalues
         DD(wl) = 0.5*(Cs(wl)+Cu(wl)+sqrt((Cs(wl)+Cu(wl))*(Cs(wl)+Cu(wl))-4*Bs(wl)*Bu(wl)))
!        write(*,*) 'DD', DD(wl)
!        write(*,*) 'SQRT', (Cs(wl)+Cu(wl))*(Cs(wl)+Cu(wl))-4*Bs(wl)*Bu(wl)
!        write(*,*) 'Cs+Cu', (Cs(wl)+Cu(wl))
!        write(*,*) 'Bs(wl)*Bu(wl)', 4*Bs(wl)*Bu(wl)
         kp(k,wl) = DD(wl)-Cu(wl)                                                !k^+
         km(k,wl) = DD(wl)-Cs(wl)                                                !k^-
         !eigenvectors
         rkp(k,wl) = Bs(wl)/DD(wl)                                               !r_k^+
         rkm(k,wl) = Bu(wl)/DD(wl)                                               !r_k^-
         kmvec(k,:) = [ rkm(k,wl), one]
         kpvec(k,:) = [ one, rkp(k,wl)]
         !exponents
         ekp(k,wl) = exp(-kp(k,wl)*dz)                                      !e_k^+
         ekm(k,wl) = exp(-km(k,wl)*dz)                                      !e_k^-
 enddo
 !the matrix, 3 diagonals. The matrix is (2n-1)x(2n-1), so are the diagonals: from 0 to 2n-2
 diagd(0,wl) = zero                           !not used
 diag(0,wl) = one                             !zeta0
 diagu(0,wl) = rkm(1,wl)*exp(-km(1,wl)*(z(2)-z(1))) !eta0
 RHS(0,wl) = EsOASIM(wl) - x(1,wl)*EdOASIM(wl)           !theta0
 do k=1,n-1
         diag (2*k-1,wl) =  rkm(k,wl)-rkm(k+1,wl)             !beta_k 
         diagu(2*k-1,wl) = -(1.-rkp(k+1,wl)*rkm(k+1,wl))       !gamma_k
         diagd(2*k-1,wl) =  (1.-rkp(k,wl)*rkm(k+1,wl))*ekp(k,wl) !alpha_k
         diag (2*k,wl) = -(rkp(k+1,wl)-rkp(k,wl))             !zeta_k 
         diagu(2*k,wl) = -(1.-rkm(k+1,wl)*rkp(k,wl))*ekm(k+1,wl) !nu_k
         diagd(2*k,wl) = 1-rkp(k,wl)*rkm(k,wl)                !eps_k
!  k   Singular
!  -
!  k+1 Singular
         if (( isSingular(k,wl)) .AND. ( isSingular(k+1,wl))) then
             RHS(2*k-1,wl) = (x(k+1,wl)-x(k,wl) - (y(k+1,wl)-y(k,wl))*rkm(k+1,wl))*z(k)*Ed(k+1,wl) !delta_k
             RHS(2*k,wl) = (y(k+1,wl)-y(k,wl) - (x(k+1,wl)-x(k,wl))*rkp(k,wl))*z(k)*Ed(k+1,wl) !theta_k
         endif

!  k   Non Singular
!  -
!  k+1 Singular
         if (( .NOT. isSingular(k,wl)) .AND. ( isSingular(k+1,wl))) then
             RHS(2*k-1,wl) = (z(k)*x(k+1,wl)-x(k,wl) - (z(k)*y(k+1,wl)-y(k,wl))*rkm(k+1,wl))*Ed(k+1,wl) !delta_k
             RHS(2*k,wl) = (z(k)*y(k+1,wl)-y(k,wl) - (z(k)*x(k+1,wl)-x(k,wl))*rkp(k,wl))*Ed(k+1,wl) !theta_k
         endif
!  k   Singular
!  -
!  k+1 Non Singular
         if (( isSingular(k,wl)) .AND. ( .NOT. isSingular(k+1,wl))) then
             RHS(2*k-1,wl) = (x(k+1,wl)-z(k)*x(k,wl) - (y(k+1,wl)-z(k)*y(k,wl))*rkm(k+1,wl))*Ed(k+1,wl) !delta_k
             RHS(2*k,wl) = (y(k+1,wl)-z(k)*y(k,wl) - (x(k+1,wl)-z(k)*x(k,wl))*rkp(k,wl))*Ed(k+1,wl) !theta_k
         endif
!  k   Non Singular
!  -
!  k+1 Non Singular
         if (( .NOT. isSingular(k,wl)) .AND. ( .NOT. isSingular(k+1,wl))) then
             RHS(2*k-1,wl) = (x(k+1,wl)-x(k,wl) - (y(k+1,wl)-y(k,wl))*rkm(k+1,wl))*Ed(k+1,wl) !delta_k
             RHS(2*k,wl) = (y(k+1,wl)-y(k,wl) - (x(k+1,wl)-x(k,wl))*rkp(k,wl))*Ed(k+1,wl) !theta_k
         endif
 enddo
 diagu(2*n-2,wl)=0.                           !not used
 cpm(2*n-1,wl)=0.                             !c^-_n=0 is just known, not a part of the solution
 if(n>1) then
        call SLAE3diag(2*n-1,nlt,diagd, diag, diagu, RHS, cpm(0:2*n-2,wl)) !solve the SLAE     
 else
        cpm(0,wl)=RHS(0,wl)
 endif
 do k=1,n                                                                        !loop over the layers
         do item=1,m                                                             !loop over the zz grid
                 if(zz(item)>=z(k) .and. zz(item)<=z(k+1)) then                  !find the horizons within the current layer
                         dz = zz(item) - z(k)                                    !from the depth to the top layer boundary
                         dz1= z(k+1) - zz(item)                                  !from the depth to the bottom layer boundary
                         E(d,item,wl) = Ed(k,wl)*exp(-cd(k,wl)*dz)                        !solution of the independent equation
                                                                                 !two other components: (bug corrected: minus sign in the 2nd exp)
                             if ( isSingular(k,wl)) then
                                 E(s:u,item,wl) = cpm(2*k-2,wl)*kpvec(k,:)*exp(-kp(k,wl)*dz) + & 
                                     & cpm(2*k-1,wl)*kmvec(k,:)*exp(-km(k,wl)*dz1)+ &
                                     & [x(k,wl), y(k,wl)]*zz(item)*E(d,item,wl)
                              else
                                 E(s:u,item,wl) = cpm(2*k-2,wl)*kpvec(k,:)*exp(-kp(k,wl)*dz) + & 
                                     & cpm(2*k-1,wl)*kmvec(k,:)*exp(-km(k,wl)*dz1)+ &
                                     & [x(k,wl), y(k,wl)]*E(d,item,wl)
                              endif
                 else if(zz(item).ge.z(k+1)) then                                !skip the rest if deeper than the current layer
                         exit
                 endif
         enddo
 enddo
 enddo !nlt
 return
 end subroutine solve_direct

!the trapezoid numerical integration
      PURE FUNCTION trapz(x,y) ! trapz integration: int{ydx}
      double precision trapz
      double precision, intent(in), dimension(:):: x, y
        integer:: n, i
        n = size(x)
        trapz = 0.
        do i = 2,n
         trapz = trapz + (y(i) + y(i-1))/2.0*(x(i)-x(i-1))
        end do
       END FUNCTION trapz
      
!the trapezoid numerical cumulative integration
      PURE FUNCTION cumtrapz(x,y) ! trapz integration
        double precision, intent(in), dimension(:):: x, y
        integer:: n
        double precision, dimension(size(x)):: cumtrapz
        integer:: i
        n = size(x)
        cumtrapz(1) = 0.
        do i = 2,n
         cumtrapz(i) = cumtrapz(i-1) + (y(i) + y(i-1))/2*(x(i)-x(i-1))
        end do
       END FUNCTION cumtrapz 

subroutine compute_3stream_adjoint
 use bioptimod_memory, only: nw,no, wavelength, Ed0m, Es0m, Eu0m, sunz,&
                             a_init, b_init, bb_init, &
                             rd_init, rs_init, ru_init, &
                             vs_init, vu_init ,opt_const_name
 use par_mod,  only: par, calc_total_parameters
 use grad_descent
 implicit none
 integer, parameter:: n=1, m=10
!integer, parameter:: n=4, m=5
 integer:: i,k,j,inw,ino
 double precision, dimension(n+1), parameter :: z=[0., 9.] !layers
 !double precision, dimension(n+1), parameter:: z=[0., 5., 10., 15., 20.] !layers
 !double precision, dimension(n+1), parameter:: z=[0., 1.25, 2.5, 3.75, 5.] !layers
 double precision, dimension(m), parameter :: zz = [((i-1)/dble(m-1)*z(n+1),i=1,m)]    !regular grid from z(1) to z(n+1)
 type(par):: p, p_n, savepar, derivs
 double precision, dimension(n,nw)   :: aa
 double precision                    :: mud
 double precision, dimension(nw)     :: rd, rs, ru, vs, vu, W, EdOASIM, EsOASIM
 double precision, dimension(3,m,nw) :: lambda, E                  !solutions to two problems
 double precision, dimension(m,nw)   :: integrand                    !grid function for numerical integration
! double precision, parameter:: fac = 1.e2, maxstep = 1.e-1        !Reasonable
 double precision, parameter         :: fac = 1.e2, maxstep = 1.e-1        !aux default step = 1.e-1
 double precision  :: step        !aux default fac = 1.e2, step = 1.e-1
 double precision  :: FF, FF_n        !aux default fac = 1.e2, step = 1.e-1
 double precision  :: WTOT
 double precision  :: X(no),low(no),Eps
 integer           :: L, iter, lvl, err
 integer           :: descent
 integer, parameter:: niter=30                                  !number of iterations
 logical, parameter:: solution=.false.                         !print the solution?
 logical, parameter:: perturbe_all=.false.                     !all params or just one or two?
 logical           :: individual_perturbation = .false. !change a, b, bb only
 logical           :: perturb_vector(9)                 !only valid if individual_perturbation .eq. true: then any parameter can be changed
 logical           :: step_back
 character*18      :: fileout
 character*7       :: w_string
 integer           :: error

 step = maxstep
 Eps  = 5.0E-5

 call getmud(sunz,mud)

 perturb_vector(:)= .FALSE.
 perturb_vector(1)= .FALSE. ! water properties
 perturb_vector(2)= .TRUE.
 perturb_vector(3)= .TRUE.
 perturb_vector(4)= .TRUE.
 perturb_vector(5)= .FALSE.
 perturb_vector(6)= .TRUE.
 perturb_vector(7)= .TRUE.

 p = par(nw,n,no)
 call p%set_const('vd',mud)
 call p%set_const('vs',vs_init)
 call p%set_const('vu',vu_init)
 call p%set_const('rd',rd_init)
 call p%set_const('rs',rs_init)
 call p%set_const('ru',ru_init)
 ! Initialize optical constituents
 !p%opt_const(:,1) =1.0D0 ! Water concentration is nomilal 1
!p%opt_const(:,2) =0.5D0 ! Diatoms [mg Chla/m3]
! p%opt_const(:,3) =0.5D0 ! Nano    [mg Chla/m3]
! p%opt_const(:,4) =0.5D0 ! Pico    [mg Chla/m3]
! p%opt_const(:,5) =0.5D0 ! Dino    [mg Chla/m3]
! p%opt_const(:,6) =0.5D0 ! CDOM    [mg C/m3]
! p%opt_const(:,7) =0.05D0 ! NAP    [mg C/m3]
 call calc_total_parameters(p)
 EdOASIM   = Ed0m(:)
 EsOASIM   = Es0m(:)

 fileout   = 'solver.txt'
 open (unit=15, file=fileout, status='unknown',    &
       access='sequential', form='formatted', action='readwrite' )

!
!    write (unit=w_string,fmt='(F6.2)') wavelength(inw)

!    write(*,*) ' EdOASIM', EdOASIM(inw) 
!    write(*,*) ' EsOASIM', EsOASIM(inw) 

     !solve the direct problem for the perturbed parameters
     call solve_direct(m, zz, n, z, nw, p, EdOASIM(:), EsOASIM(:), E(:,:,:))

 do inw=1,nw
     if(solution) then
             print*, 'The perturbed solution E:'
             print*, E(d,:,1)
             print*, E(s,:,1)
             print*, E(u,:,1)
     endif
     !the derivative of the functional wrt to the value to be compared with the
     !observed
     W(inw) = fac*(E(u,1,inw) - Eu0m(inw))
 enddo

     !solve the adjoint problem and get the derivatives
     call get_derivs(m, zz, n, z, p, W, E, derivs, lambda)
     if(solution) then
             print*, 'The solution lambda:'
             print*, lambda(d,:,1)
             print*, lambda(s,:,1)
             print*, lambda(u,:,1)
     endif

!    print*, 'The derivatives:'
!    call derivs%print('d/d')
     !try to improve the set of parameters
!    p_n = par(p)
  descent=1
  if (descent .EQ. 0 ) then  
     iter=0
     do iter=1,niter
             print*,'====================', iter
             !step backward along the grad
             do ino =2,no
                 if (perturb_vector(ino))  p%opt_const(1,ino) =p%opt_const(1,ino) -derivs%opt_const(1,ino) *step
             enddo
             call calc_total_parameters(p)

             !solve the direct problem for the perturbed parameters
             do inw=1,nw
             step_back = .FALSE.
             if (p%a(1,inw,2) .LT. 0.) then
                     step_back = .TRUE.
             endif
             enddo

             if (step_back) then
                     do ino =2,no
                        p%opt_const(1,ino) =p%opt_const(1,ino) + derivs%opt_const(1,ino) *step
                     enddo
                     step =step/1.5
                     call calc_total_parameters(p)
             endif
             if (step .LT. maxstep/500.) then
                     write(*,*) "Failed to converge"
                     STOP
             endif
                call solve_direct(m, zz, n, z, nw, p, EdOASIM, EsOASIM, E(:,:,:))
             WTOT =0.0D0
             do inw=1,nw
                 W = fac*(E(u,1,inw) - Eu0m(inw))
                 write(*,*)  'Wavelenght: ', wavelength(inw), ',Functional  score: ',  fac*0.5*(E(u,1,inw) - Eu0m(inw))**2 !the functional itself
                 write(*,*)   wavelength(inw), ', Eu(model): ',  E(u,1,inw), ', Eu(sat): ', Eu0m(inw)
                 WTOT = WTOT + fac*0.5*(E(u,1,inw) - Eu0m(inw))**2
             enddo
             write(*,*)  'total functional score ', WTOT
             do ino =1,no
                 write(*,*)  'Optical constituent: ',  TRIM(opt_const_name(ino)), &
                    ' concentration  ',  p%opt_const(1,ino)
             enddo

             write(*,*) '============'

             if(solution) then
                     print*, 'The improved solution E:'
                     print*, E(d,:,inw)
                     print*, E(s,:,inw)
                     print*, E(u,:,inw)
             endif
             !the derivative of the functional wrt to the value to be compared with
             !the observed

             !solve the adjoint problem and get the derivatives
             call get_derivs(m, zz, n, z, p, W, E, derivs, lambda)
             if(solution) then
                     print*, 'The solution lambda:'
                     print*, lambda(d,:,1)
                     print*, lambda(s,:,1)
                     print*, lambda(u,:,1)
             endif
!            print*, 'The derivatives:'
!            call derivs%print('d/d')
!            call p%print('   ')
!            call print_par(p/savepar)
     enddo
     close(unit=15)

    endif

    if (descent .EQ. 1 ) then  
       do j=1,no
            X(j) = p%opt_const(1,j)
            low(j) = 0.01D0
       enddo
       call FM34(p, functional, gradient, X, Eps, error,verbose=.true.,low=low)
    endif

    open(15,file='./res_bio.txt',status='new')
    write(15,22) "LAYER   ", "W       ", "P1      ","P2      ","P3      ","P4      ","CDOM    ","NAP     "

    do k=1,n
         write(15,23) k, (p%opt_const(k,j),j=1,no)
    enddo

    close(15)


    do inw=1,nw
    write (unit=w_string,fmt='(F6.2)') wavelength(inw)
    fileout   = 'res_opt' // TRIM(w_string) // '.txt'
    open(15,file=fileout,status='new')
    write(15,22) "LAYER   ", "Z_UP     ", "Z_DOWN  ","ZEN     ","EdOAS   ","EsOAS   ","EuSAT   ","EdMOD   ","EsMOD   ", & 
                 "EuMOD   ", "EuMODtop"

    do k=1,n
         write(15,23) k, z(k), z(k+1), sunz, EdOASIM(inw), EsOASIM(inw), Eu0m(inw), E(d,m,inw), E(s,m,inw), & 
                      E(u,m,inw), E(u,1,inw)
    enddo

    close(15)
    enddo

22   FORMAT(A8,A8,A8,A9,A8,A8,A8,A8,A8,A8,A8)
23   FORMAT(I8,f8.5,f8.5,f9.5,f8.5,f8.5,f8.5,f8.5,f8.5,f8.5,f8.5)

     return

 contains

     function functional(mypar,X, const)
         double precision:: functional
         type(par),intent(inout) :: mypar
         double precision, dimension(:), intent(in), optional:: const     !any extra parameters
         double precision, dimension(:), intent(inout) :: X     !any extra parameters
         do ino =2,no
             if (perturb_vector(ino))  then
                     mypar%opt_const(1,ino) = X(ino)
             endif
         enddo
         call calc_total_parameters(mypar)
         call solve_direct(m, zz, n, z, nw, mypar, EdOASIM(:), EsOASIM(:), E(:,:,:))
         functional = 0.0D0
         do inw=1,nw
             functional = functional + fac*0.5*(E(u,1,inw) - Eu0m(inw))**2
         enddo
     end function functional

     function gradient(mypar,X,grad,const)
         implicit none
         double precision:: gradient
         type(par), intent(inout) :: mypar
         double precision, dimension(:), intent(in), optional:: const     !any extra parameters
         double precision, dimension(:), intent(inout) :: X
         double precision, dimension(:), intent(out) :: grad    !any extra parameters
         double precision :: grad2 
         double precision :: myW(nw)
         type(par)        :: derivs
         do ino =2,no
             if (perturb_vector(ino))  then
                     mypar%opt_const(1,ino) = X(ino)
             endif
         enddo
         call calc_total_parameters(mypar)
         call solve_direct(m, zz, n, z, nw, mypar, EdOASIM(:), EsOASIM(:), E(:,:,:))
         do inw=1,nw
             myW(inw) = fac*(E(u,1,inw) - Eu0m(inw))
         enddo
         call get_derivs(m, zz, n, z, mypar, myW, E, derivs)

         grad(:) = 0.0D0
         grad2 = 0.0D0

         do ino =2,no
             if (perturb_vector(ino))  then
                     grad(ino) = derivs%opt_const(1,ino)
                     grad2     = derivs%opt_const(1,ino) * derivs%opt_const(1,ino)
             endif
         enddo
         gradient = sqrt(grad2)
     end function gradient


end subroutine compute_3stream_adjoint

END MODULE adj_3stream

!program test
!use adj_3stream, only: test_3stream_adjoint
!call test_3stream_adjoint

!end program test

