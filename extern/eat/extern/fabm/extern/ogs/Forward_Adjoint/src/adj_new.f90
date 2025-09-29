MODULE adj_3stream

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
subroutine get_derivs(m, zz, n, z, a, b, bb, rd, rs, ru, vd, vs, vu, W, E, derivs, sol_adj)
 implicit none
 integer, intent(in):: n, m
 double precision, dimension(m), intent(in):: zz !must be some grid from z(1)=0 to z(n+1)=bottom
 double precision, dimension(n), intent(in):: a, b, bb, vd
 double precision, dimension(n+1), intent(in):: z
 double precision, intent(in):: rd, rs, ru, vs, vu, W
 double precision, dimension(3,m), intent(in):: E
 double precision, dimension(4*n+5), intent(out):: derivs !wrt a(n), b(n), b(n), vd(n), also wrt rd, rs, ru, vs, vu
 double precision, dimension(3,m), intent(out), optional:: sol_adj
 double precision, dimension(3,m):: lambda
 double precision, dimension(m):: integrand
 integer:: L, k
 call solve_adj(m, zz, n, z, a, b, bb, rd, rs, ru, vd, vs, vu, W, lambda)
 if(present(sol_adj)) sol_adj = lambda
 !first let us get the derivatives wrt a(n)
 do k=1,n
     integrand = 0.
     do L=1,m                                            !loop over the grid
         if(zz(L).ge.z(k) .and. zz(L).lt.z(k+1)) then    !skip points above or below the layer
             integrand(L) = -lambda(d,L)*E(d,L)/vd(k) - lambda(s,L)*E(s,L)/vs + lambda(u,L)*E(u,L)/vu !the expression under the integral
         endif
     enddo !L: grid step within a layer
     derivs(k)= -trapz(zz, integrand)                    !integrating
 enddo !k: the layer
 !now, let us get the derivatives wrt to b(n)
 do k=1,n
     integrand = 0.
     do L=1,m
         if(zz(L).ge.z(k) .and. zz(L).lt.z(k+1)) then
             integrand(L) = ( -lambda(d,L) + lambda(s,L) ) / vd(k) * E(d,L)
         endif
     enddo !L: grid step within a layer
     derivs(n+k)= -trapz(zz, integrand)
 enddo !k: the layer
 !now, let us get the derivatives wrt to bb(n)
 do k=1,n
     integrand = 0.
     do L=1,m
         if(zz(L).ge.z(k) .and. zz(L).lt.z(k+1)) then
             integrand(L) = -(lambda(u,L)+lambda(s,L)) * (rd/vd(k)*E(d,L)+rs/vs*E(s,L)-ru/vu*E(u,L))
         endif
     enddo !L: grid step within a layer
     derivs(2*n+k)= -trapz(zz, integrand)
 enddo !k: the layer
 !now, let us get the derivatives wrt to vd(n)
 do k=1,n
     integrand = 0.
     do L=1,m
         if(zz(L).ge.z(k) .and. zz(L).lt.z(k+1)) then
             integrand(L) = ( lambda(d,L)*(a(k)+b(k)) - lambda(s,L)*(b(k)-rd*bb(k)) + lambda(u,L)*rd*bb(k) )*E(d,L)/vd(k)/vd(k)
         endif
     enddo !L: grid step within a layer
     derivs(3*n+k)= -trapz(zz, integrand)
 enddo !k: the layer
 !now, let us get the derivatives wrt to vs
 integrand = 0.
 do k=1,n
     do L=1,m
         if(zz(L).ge.z(k) .and. zz(L).lt.z(k+1)) then
             integrand(L) = ( lambda(s,L)*(a(k)+rs*bb(k)) + lambda(u,L)*rs*bb(k) )*E(s,L)/vs/vs
         endif
     enddo !L: grid step within a layer
 enddo !k: the layer
 derivs(4*n+1)= -trapz(zz, integrand)
 !now, let us get the derivatives wrt to vu
 integrand = 0.
 do k=1,n
     do L=1,m
         if(zz(L).ge.z(k) .and. zz(L).lt.z(k+1)) then
             integrand(L) = -( lambda(u,L)*(a(k)+ru*bb(k)) + lambda(s,L)*ru*bb(k) )*E(u,L)/vu/vu
         endif
     enddo !L: grid step within a layer
 enddo !k: the layer
 derivs(4*n+2)= -trapz(zz, integrand)
 !now, let us get the derivatives wrt to rd
 integrand = 0.
 do k=1,n
     do L=1,m
         if(zz(L).ge.z(k) .and. zz(L).lt.z(k+1)) then
             integrand(L) = -( lambda(u,L) + lambda(s,L) )*E(d,L)*bb(k)/vd(k)
         endif
     enddo !L: grid step within a layer
 enddo !k: the layer
 derivs(4*n+3)= -trapz(zz, integrand)
 !now, let us get the derivatives wrt to rs
 integrand = 0.
 do k=1,n
     do L=1,m
         if(zz(L).ge.z(k) .and. zz(L).lt.z(k+1)) then
             integrand(L) = -( lambda(u,L) + lambda(s,L) )*E(s,L)*bb(k)/vs
         endif
     enddo !L: grid step within a layer
 enddo !k: the layer
 derivs(4*n+4)= -trapz(zz, integrand)
 !now, let us get the derivatives wrt to ru
 integrand = 0.
 do k=1,n
     do L=1,m
         if(zz(L).ge.z(k) .and. zz(L).lt.z(k+1)) then
             integrand(L) = ( lambda(u,L) + lambda(s,L) )*E(u,L)*bb(k)/vu
         endif
     enddo !L: grid step within a layer
 enddo !k: the layer
 derivs(4*n+5)= -trapz(zz, integrand)
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
subroutine solve_adj(m, zz, n, z, a, b, bb, rd, rs, ru, vd, vs, vu, W, lambda)
use Tridiagonal, only: SLAE3diag
 implicit none
 integer, intent(in):: n, m
 double precision, dimension(m), intent(in):: zz !an argument, must be between z(1) and z(n+1)
 double precision, dimension(n), intent(in):: a, b, bb, vd
 double precision, dimension(n+1), intent(in):: z
 double precision, dimension(3,m),intent(out):: lambda
 double precision, intent(in):: rd, rs, ru, vs, vu, W
 double precision, dimension(0:2*n-1):: diag, diagu, diagd, RHS
 double precision, dimension(n):: Ad, km, kp
 double precision, dimension(n,3):: Advec, kmvec, kpvec
 double precision, dimension(2*n):: psipm
        !prepare the 3-diagonal matrix:
           call matrix_adj_3stream(n, z, a, b, bb, rd, rs, ru, vd, vs, vu, W, &  !input
                   & diag, diagu, diagd, RHS,                                 &  !output: 3-diag matrix and RHS
                   & Ad, kp, km,                                              &  !eigenvalues
                   & Advec,kmvec,kpvec                                        &  !eigenvectors
                   & )                                
!               print*, 'eigen kp:', kp
!               print*, 'eigen km:', km
           call SLAE3diag(2*n,1,diagd, diag, diagu, RHS, psipm)                  !solve the SLAE for psi
           call get_lambda(m, zz, n,z, Ad, kp, km, Advec,kpvec,kmvec, psipm, lambda) !build the solution using the obtained psi
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
subroutine matrix_adj_3stream(n, z, a, b, bb, rd, rs, ru, vd, vs, vu, W,    &  !input
        & diag, diagu, diagd, RHS,                                          &  !output: 3-diag matrix and RHS
        & Ad, kp, km,                                                       &  !eigenvalues
        & Advec,kmvec,kpvec                                                 &  !eigenvectors
        & )                                
 implicit none
 integer, intent(in):: n
 double precision, dimension(n), intent(in):: a, b, bb, vd
 double precision, dimension(n+1), intent(in):: z
 double precision, intent(in):: rd, rs, ru, vs, vu, W
 double precision, dimension(0:2*n-1), intent(out):: diag, diagu, diagd, RHS
 double precision, dimension(n), intent(out):: Ad, km, kp
 double precision, dimension(n,3), intent(out):: Advec, kmvec, kpvec
 double precision, dimension(n):: rkp,rkm, RRkp, RRkm
 double precision:: Fd, Bd, Bs, Bu, Cs, Cu, DD, dz
 integer:: k
 do k=1,n
         !matrix
         Ad(k) = (a(k)+b(k))/vd(k)
         Fd = (b(k)-rd*bb(k))/vd(k)
         Bd = rd*bb(k)/vd(k)
         Cs = (a(k)+rs*bb(k))/vs
         Cu = (a(k)+ru*bb(k))/vu
         Bs = rs*bb(k)/vs
         Bu = ru*bb(k)/vu
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
 subroutine solve_direct(m, zz, n, z, nlt, a, b, bb, rd, rs, ru, vd, vs, vu, EdOASIM, EsOASIM, E, E_ave)
 use Tridiagonal, only: SLAE3diag
 implicit none
 integer, intent(in):: n,m                                             !n of layers and size of the vertical grid
 integer, intent(in):: nlt                                             !n of wavelenghts to be considered
 double precision, dimension(m), intent(in):: zz                       !must be some grid from z(1)=0 to z(n+1)=bottom
 double precision, dimension(n+1), intent(in):: z                      !layer boundaries (depth levels); z(1)=0 (must be), z(n+1) = bottom
 double precision, dimension(n,nlt), intent(in):: a, b, bb, vd         !input depth-dependent data
 double precision, dimension(n,nlt)  :: rvd 
 double precision, intent(in):: rd, rs, ru, vs, vu                     !input parameters
 double precision, dimension(nlt),intent(in) :: EdOASIM, EsOASIM                       !boundary values
 double precision, dimension(3,m,nlt), intent(out):: E                     !the 3-stream solution on the zz grid
 double precision, dimension(3,n,nlt), intent(out):: E_ave                 !the 3-stream solutionaveraged over the discrete levels
 double precision, dimension(n,nlt):: cd, x, y, kp, km, rkp, rkm, ekp, ekm !stuff used in the formulae, see the Overleaf doc https://www.overleaf.com/project/5cbefd4a70921e1466457de2
 double precision, dimension(n,2,nlt):: kmvec, kpvec                       !eigenvectors
 double precision, dimension(n+1,nlt):: Ed                                 !solution of the first equation (which is independent of the other two) on the z grid
 double precision :: dz, dz1
 double precision, dimension(nlt) :: Fd, Bd, Cs, Cu, Bs, Bu, DD                !aux variables
 double precision, dimension(0:2*n-2,nlt):: diag, diagd, diagu, RHS        !the matrix, the right-hand side
 double precision, dimension(0:2*n-1,nlt):: cpm                            !the solution c_+ (0,2,4,...), c_- (1,3,5,...), up to cpm(2n-1)=c^-_n=0
 double precision  :: aux_ave1, aux_ave2
 integer:: k, item,wl                                                     !counters
 Ed(1,:) = EdOASIM(:)                                                    !boundary value: Ed at the surface is just given
 rvd=min(1.0D0/vd(:,:),1.5D0)
 rvd=max(rvd,0.0D0)

 do k=1,n                                                             !for each layer
         !matrix
         dz = z(k+1)-z(k)                                             !layer thickness
         cd(k,:) = (a(k,:)+b(k,:))*rvd(k,:)                                    !coefficient of the separate equation
         Ed(k+1,:) = Ed(k,:)*exp(-cd(k,:)*dz)                               !solution of the separate equation
         Fd(:) = (b(k,:)-rd*bb(k,:))*rvd(k,:)                                   !components of the auxiliary matrix for the partial solution
         Bd(:) = rd*bb(k,:)*rvd(k,:)
         Cs(:) = (a(k,:)+rs*bb(k,:))/vs
         Cu(:) = (a(k,:)+ru*bb(k,:))/vu
         Bs(:) = rs*bb(k,:)/vs
         Bu(:) = ru*bb(k,:)/vu
         !aux variables x,y
         x(k,:) = (-Fd(:)*(Cu(:)+cd(k,:)) -Bd(:)*Bu(:))  / ((cd(k,:)-Cs(:))*(cd(k,:)+Cu(:))+Bs(:)*Bu(:))       !explicit solution of the aux SLAE
         y(k,:) = (-Fd(:)*Bs(:) +Bd(:)*(-Cs(:)+cd(k,:))) / ((cd(k,:)-Cs(:))*(cd(k,:)+Cu(:))+Bs(:)*Bu(:))
         !eigenvalues
         DD(:) = 0.5*(Cs(:)+Cu(:)+sqrt((Cs(:)+Cu(:))*(Cs(:)+Cu(:))-4*Bs(:)*Bu(:)))
         kp(k,:) = DD(:)-Cu(:)                                                !k^+
         km(k,:) = DD(:)-Cs(:)                                                !k^-
         !eigenvectors
         rkp(k,:) = Bs(:)/DD(:)                                               !r_k^+
         rkm(k,:) = Bu(:)/DD(:)                                               !r_k^-
         do wl = 1,nlt ! lopp on wavelenghts
             kmvec(k,:,wl) = [ rkm(k,wl), one]
             kpvec(k,:,wl) = [ one, rkp(k,wl)]
         enddo
         !exponents
         ekp(k,:) = exp(-kp(k,:)*dz)                                      !e_k^+
         ekm(k,:) = exp(-km(k,:)*dz)                                      !e_k^-
 enddo
 !the matrix, 3 diagonals. The matrix is (2n-1)x(2n-1), so are the diagonals: from 0 to 2n-2
 diagd(0,:) = zero                           !not used
 diag(0,:) = one                             !zeta0
 diagu(0,:) = rkm(1,:)*exp(-km(1,:)*(z(2)-z(1))) !eta0
 RHS(0,:) = EsOASIM(:) - x(1,:)*EdOASIM(:)           !theta0
 do k=1,n-1
         diag (2*k-1,:) =  rkm(k,:)-rkm(k+1,:)             !beta_k 
         diagu(2*k-1,:) = -(1.-rkp(k+1,:)*rkm(k+1,:))       !gamma_k
         diagd(2*k-1,:) =  (1.-rkp(k,:)*rkm(k+1,:))*ekp(k,:) !alpha_k
         RHS(2*k-1,:) = (x(k+1,:)-x(k,:) - (y(k+1,:)-y(k,:))*rkm(k+1,:))*Ed(k+1,:) !delta_k
         diag (2*k,:) = -(rkp(k+1,:)-rkp(k,:))             !zeta_k 
         diagu(2*k,:) = -(1.-rkm(k+1,:)*rkp(k,:))*ekm(k+1,:) !nu_k
         diagd(2*k,:) = 1-rkp(k,:)*rkm(k,:)                !eps_k
         RHS(2*k,:) = (y(k+1,:)-y(k,:) - (x(k+1,:)-x(k,:))*rkp(k,:))*Ed(k+1,:) !theta_k
 enddo
 diagu(2*n-2,:)=0.                           !not used
 cpm(2*n-1,:)=0.                             !c^-_n=0 is just known, not a part of the solution
 if(n>1) then
        call SLAE3diag(2*n-1,nlt,diagd, diag, diagu, RHS, cpm(0:2*n-2,:)) !solve the SLAE     
 else
        cpm(0,:)=RHS(0,:)
 endif
 do k=1,n                                                                        !loop over the layers
         do item=1,m                                                             !loop over the zz grid
                 if(zz(item)>=z(k) .and. zz(item)<=z(k+1)) then                  !find the horizons within the current layer
                         dz = zz(item) - z(k)                                    !from the depth to the top layer boundary
                         dz1= z(k+1) - zz(item)                                  !from the depth to the bottom layer boundary
                         E(d,item,:) = Ed(k,:)*exp(-cd(k,:)*dz)                        !solution of the independent equation
                                                                                 !two other components: (bug corrected: minus sign in the 2nd exp)
                         do wl=1,nlt                                                         
                             E(s:u,item,wl) = cpm(2*k-2,wl)*kpvec(k,:,wl)*exp(-kp(k,wl)*dz) + & 
                                     & cpm(2*k-1,wl)*kmvec(k,:,wl)*exp(-km(k,wl)*dz1)+ &
                                     & [x(k,wl), y(k,wl)]*E(d,item,wl)
                         enddo
                 else if(zz(item).ge.z(k+1)) then                                !skip the rest if deeper than the current layer
                         exit
                 endif
         enddo
 enddo
! Averages over layers
 do k=1,n                                                                        !loop over the layers
                         dz = z(k+1) - z(k)                                    !from the depth to the top layer boundary
                         E_ave(d,k,:) = Ed(k,:)*(1.0-exp(-cd(k,:)*dz))/(cd(k,:)*dz)              !solution of the independent equation
                                                                                 !two other components: (bug corrected: minus sign in the 2nd exp)
                         do wl=1,nlt                                                         
                             aux_ave1 = (1.0-exp(-kp(k,wl)*dz))/(kp(k,wl)*dz)
                             aux_ave2 = (1.0-exp(-km(k,wl)*dz))/(km(k,wl)*dz)
                             E_ave(s:u,k,wl) = cpm(2*k-2,wl)*kpvec(k,:,wl)*aux_ave1 + & 
                                     & cpm(2*k-1,wl)*kmvec(k,:,wl)*aux_ave2 + &
                                     & [x(k,wl), y(k,wl)]*E_ave(d,k,wl)
                         enddo
 enddo
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

subroutine test_3stream_adjoint
 use bioptimod_memory, only: nw       
 use grad_descent, only: FM34 
 implicit none
 integer, parameter:: n=4, m=5
 integer:: i
 double precision, dimension(n+1), parameter:: z=[0., 1.25, 2.5, 3.75, 5.]                             !layers
 double precision, dimension(m), parameter:: zz = [((i-1)/dble(m-1)*z(n+1),i=1,m)]    !regular grid from z(1) to z(n+1)
 double precision, dimension(n,nw):: a, b, bb, vd
 double precision, dimension(n,nw):: aa
 double precision, dimension(nw) :: rd, rs, ru, vs, vu, W
 double precision,dimension(nw) :: EdOASIM, EsOASIM
 double precision, dimension(4*n+5,nw):: savepar, pp 
 double precision, dimension(4*n+5,nw):: derivs                   !wrt a(n), b(n), b(n), vd(n), also wrt vs, vu, rd, rs, ru
 double precision, dimension(3,m,nw):: lambda, E                  !solutions to two problems
 double precision, dimension(3,n,nw):: E_ave                      !
 double precision, dimension(m):: integrand                    !grid function for numerical integration
 double precision, parameter:: fac = 1.e2, step = 1.e-1        !aux
 double precision:: Eu0_given = 0.577999                       !the "observed" value
 integer:: L, iter ,err
 integer, parameter:: niter=7                                 !number of iterations
 logical, parameter:: solution=.not..false.                    !print the solution?
 logical, parameter:: perturbe_all=.false.                     !all params or just one or two?
 !some values, just for testing
 a = 1.0D0; b=1.0D0; bb=1.0D0; vd=0.42D0;
 rd=1.0D0; rs=1.5D0; ru=3.0D0; vs=0.83D0; vu=0.4D0; 
 !save the values for future comparing
 savepar(1:n,:) = a; savepar(n+1:2*n,:) = b; savepar(2*n+1:3*n,:) = bb; savepar(3*n+1:4*n,:) = vd; 
 savepar(4*n+1,:) = vs;  savepar(4*n+2,:) = vu;  savepar(4*n+3,:) = rd;  savepar(4*n+4,:) = rs;  savepar(4*n+5,:) = ru; 
 savepar = savepar/100. !to be in %
 EdOASIM(1) = 2.0
 EsOASIM(1) = 1.5
 !solve the direct problem
 call solve_direct(m, zz, n, z, nw, a, b, bb, rd(1), rs(1), ru(1), vd, vs(1), vu(1), EdOASIM, EsOASIM, E, E_ave)
 if(solution) then
         print*, 'The solution E:'
         print*, E(d,:,1)
         print*, E(s,:,1)
         print*, E(u,:,1)
 endif
 !let the obtained value be the observed value
 Eu0_given = E(u,1,1)
 !change the parameters
 a=a*1.1;
 if(perturbe_all) then
     b=b*1.01; 
     bb=bb*1.01; 
     vd=vd*1.01;
     vs=vs*1.01; 
     vu=vu*1.01; 
     rd=rd*1.01; 
     rs=rs*1.01; 
     ru=ru*1.01; 
 endif
 !solve the direct problem for the perturbed parameters
 call solve_direct(m, zz, n, z, nw, a, b, bb, rd(1), rs(1), ru(1), vd, vs(1), vu(1), EdOASIM, EsOASIM, E, E_ave)
 if(solution) then
         print*, 'The perturbed solution E:'
         print*, E(d,:,1)
         print*, E(s,:,1)
         print*, E(u,:,1)
 endif
 print*, 'Functional: ', fac*0.5*(E(u,1,1) - Eu0_given)**2 !the functional itself
 !the derivative of the functional wrt to the value to be compared with the observed
 W = fac*(E(u,1,1) - Eu0_given)
 !solve the adjoint problem and get the derivatives
 call get_derivs(m, zz, n, z, a, b, bb, rd(1), rs(1), ru(1), vd, vs(1), vu(1), W(1), E, derivs, lambda)
 if(solution) then
         print*, 'The solution lambda:'
         print*, lambda(d,:,1)
         print*, lambda(s,:,1)
         print*, lambda(u,:,1)
 endif
 print*, 'W=',W
 print*, 'The derivatives:'
 print*, 'd/da:', derivs(1:n,1)
 print*, 'd/db:', derivs(n+1:2*n,1)
 print*, 'd/dbb:', derivs(2*n+1:3*n,1)
 print*, 'd/dvd:', derivs(3*n+1:4*n,1)
 print*, 'd/d(vs, vu, rd, rs, ru):', derivs(4*n+1:4*n+5,1)
 !try to improve the set of parameters
 do iter=1,niter
         print*,'===================='
         !step backward along the grad
         a=a-derivs(1:n,:)*step/iter
         if(perturbe_all) then
             b=b-derivs(n+1:2*n,:)*step
             bb=bb-derivs(2*n+1:3*n,:)*step
             vd=vd-derivs(3*n+1:4*n,:)*step*0.1 !this one is quite sensible
             vs=vs-derivs(4*n+1,:)*step
             vu=vu-derivs(4*n+2,:)*step
             rd=rd-derivs(4*n+3,:)*step
             rs=rs-derivs(4*n+4,:)*step
             ru=ru-derivs(4*n+5,:)*step
         endif
         !solve the direct problem for the perturbed parameters
         call solve_direct(m, zz, n, z, 1, a, b, bb, rd(1), rs(1), ru(1), vd, vs(1), vu(1), EdOASIM, EsOASIM, E, E_ave)
         if(solution) then
                 print*, 'The improved solution E:'
                 print*, E(d,:,1)
                 print*, E(s,:,1)
                 print*, E(u,:,1)
         endif
         print*, 'Functional: ', fac*0.5*(E(u,1,1) - Eu0_given)**2 !the functional value
         !the derivative of the functional wrt to the value to be compared with the observed
         W = fac*(E(u,1,1) - Eu0_given)
         !solve the adjoint problem and get the derivatives
         call get_derivs(m, zz, n, z, a, b, bb, rd(1), rs(1), ru(1), vd, vs(1), vu(1), W(1), E, derivs, lambda)
         if(solution) then
                 print*, 'The solution lambda:'
                 print*, lambda(d,:,1)
                 print*, lambda(s,:,1)
                 print*, lambda(u,:,1)
         endif
         print*, 'The derivatives:'
         print*, 'd/da:', derivs(1:n,1)
         print*, 'd/db:', derivs(n+1:2*n,1)
         print*, 'd/dbb:', derivs(2*n+1:3*n,1)
         print*, 'd/dvd:', derivs(3*n+1:4*n,1)
         print*, 'd/d(vs, vu, rd, rs, ru):', derivs(4*n+1:4*n+5,1)
         print*, 'parameters:'
         print*, '    a         ', a
         print*, '    b         ', b
         print*, '    bb        ', bb
         print*, '    vd        ', vd
         print*, 'vs,vu,rd,rs,ru', vs,vu,rd,rs,ru
         print*, 'parameters wrt the original values:'
         print*, '    a         ', a/savepar(1:n,:)
         print*, '    b         ', b/savepar(n+1:2*n,:)
         print*, '    bb        ', bb/savepar(2*n+1:3*n,:)
         print*, '    vd        ', vd/savepar(3*n+1:4*n,:)
         print*, 'vs,vu,rd,rs,ru', [vs,vu,rd,rs,ru]/savepar(4*n+1:4*n+5,1)
 enddo

 print*, ''
 print*, '================================================='
 print*, '============ DESCENT ============================'
 print*, '================================================='
 pp = (100.*savepar)*1.10
 print*, 'parameters wrt the original values, %:', pp/savepar
 print*, 'func(par)=', functional(pp)
 print*, '================================================='
! call FM34(functional, Gradient, pp, 1.0d-3, err)
 err = 0
 select case(err)
 case(0)
     print*, 'OK!!!'
 case(66)
     print*, 'Error', err, "'Number of iterations kmax exceeded'"
 case(67)
     print*, 'Error', err, "'Number of iterations cmax exceeded'"
 case(65)
     print*, 'Error', err, "'Possible instability'"
 case default
     print*, 'Unknown error', err
 end select
 print*, '================================================='
 print*, 'parameters wrt the original values, %:', pp/savepar
 print*, 'func(par)=', functional(pp)
 print*, '|grad(par)|=', gradient(pp,derivs)
 print*, 'grad(par)=', derivs
 print*, '================================================='

 contains

     function functional(par, const) 
         double precision:: functional
         double precision,dimension(:,:),intent(in):: par
         double precision, dimension(:,:), intent(in), optional:: const     !any extra parameters
         a=par(1:n,:)
         b=par(n+1:2*n,:)
         bb=par(2*n+1:3*n,:)
         vd=par(3*n+1:4*n,:)
         vs=par(4*n+1,:)
         vu=par(4*n+2,:)
         rd=par(4*n+3,:)
         rs=par(4*n+4,:)
         ru=par(4*n+5,:)
         call solve_direct(m, zz, n, z, 1, a, b, bb, rd(1), rs(1), ru(1), vd, vs(1), vu(1), EdOASIM, EsOASIM, E, E_ave)
         functional =  fac*0.5*(E(u,1,1) - Eu0_given)**2 !the functional value
     end function functional

     function gradient(par,grad,const) 
         double precision:: gradient
         double precision,dimension(:,:),intent(in):: par
         double precision,dimension(:,:),intent(out):: grad
         double precision, dimension(:), intent(in), optional:: const     !any extra parameters
         a=par(1:n,:)
         b=par(n+1:2*n,:)
         bb=par(2*n+1:3*n,:)
         vd=par(3*n+1:4*n,:)
         vs=par(4*n+1,:)
         vu=par(4*n+2,:)
         rd=par(4*n+3,:)
         rs=par(4*n+4,:)
         ru=par(4*n+5,:)
         call solve_direct(m, zz, n, z, 1, a, b, bb, rd(1), rs(1), ru(1), vd, vs(1), vu(1), EdOASIM, EsOASIM, E, E_ave)
         W = fac*(E(u,1,1) - Eu0_given)
         call get_derivs(m, zz, n, z, a, b, bb, rd(1), rs(1), ru(1), vd, vs(1), vu(1), W(1), E, grad, lambda)
         gradient = sqrt(sum(grad*grad))
     end function gradient

end subroutine test_3stream_adjoint

subroutine compute_3stream_adjoint
 use bioptimod_memory, only: nw, wavelength, Ed0m, Es0m, Eu0m
 use grad_descent, only: FM34 
 implicit none
 integer, parameter:: n=4, m=5
 integer:: i,inw
 double precision, dimension(n+1), parameter:: z=[0., 5., 10., 15., 20.] !layers
 !double precision, dimension(n+1), parameter:: z=[0., 1.25, 2.5, 3.75, 5.] !layers
 double precision, dimension(m), parameter:: zz = [((i-1)/dble(m-1)*z(n+1),i=1,m)]    !regular grid from z(1) to z(n+1)
 double precision, dimension(n,nw):: a, b, bb, vd
 double precision, dimension(n,nw):: aa
 double precision, dimension(nw) :: rd, rs, ru, vs, vu, W, EdOASIM, EsOASIM
 double precision, dimension(4*n+5,nw):: savepar, pp 
 double precision, dimension(4*n+5,nw):: derivs                   !wrt a(n), b(n), b(n), vd(n), also wrt vs, vu, rd, rs, ru
 double precision, dimension(3,m,nw):: lambda, E                  !solutions to two problems
 double precision, dimension(3,n,nw):: E_ave                      !
 double precision, dimension(m,nw):: integrand                    !grid function for numerical integration
 double precision, parameter:: fac = 1.e2, step = 1.e-1        !aux default step = 1.e-1
 double precision:: Eu0_given = 0.577999                       !the "observed" value
 integer:: L, iter, lvl, err
 integer, parameter:: niter=7                                  !number of iterations
 logical, parameter:: solution=.false.                         !print the solution?
 logical, parameter:: perturbe_all=.false.                     !all params or just one or two?
 logical           :: perturb_vector(9)  
 character*12      :: fileout
 character*4       :: w_string

 perturb_vector(:)= .FALSE.
 perturb_vector(1)= .TRUE.
 perturb_vector(2)= .TRUE.
 perturb_vector(3)= .TRUE.
 do inw=1,nw
 !some values, just for testing
     a = 0.1D0; b=0.1D0; bb=0.05D0; vd=0.42D0;
     rd=1.0D0; rs=1.5D0; ru=3.0D0; vs=0.83D0; vu=0.4D0; 
 !save the values for future comparing
     savepar(1:n,:) = a; savepar(n+1:2*n,:) = b; savepar(2*n+1:3*n,:) = bb;
     savepar(3*n+1:4*n,:) = vd; 
     savepar(4*n+1,:) = vs;  savepar(4*n+2,:) = vu;  savepar(4*n+3,:) = rd;
     savepar(4*n+4,:) = rs;  savepar(4*n+5,:) = ru; 
     savepar = savepar/100. !to be in %
!
     write (unit=w_string,fmt='(I0.4)') wavelength(inw)
     fileout   = 'sol_' // TRIM(w_string) // '.txt'
     open (unit=15, file=fileout, status='unknown',    &
           access='sequential', form='formatted', action='readwrite' )

     EdOASIM   = Ed0m(:)
     EsOASIM   = Es0m(:)
     Eu0_given = Eu0m(1)
!    Eu0_given = E(u,1)
     !solve the direct problem for the perturbed parameters
     call solve_direct(m, zz, n, z, 1, a, b, bb, rd(1), rs(1), ru(1), vd, vs(1), vu(1), EdOASIM, EsOASIM, E, E_ave)
     if(solution) then
             print*, 'The perturbed solution E:'
             print*, E(d,:,1)
             print*, E(s,:,1)
             print*, E(u,:,1)
     endif
     print*, 'Functional: ', fac*0.5*(E(u,1,1) - Eu0_given)**2 !the functional itself
     !the derivative of the functional wrt to the value to be compared with the
     !observed
     W(1) = fac*(E(u,1,1) - Eu0_given)
     !solve the adjoint problem and get the derivatives
     call get_derivs(m, zz, n, z, a, b, bb, rd(1), rs(1), ru(1), vd, vs(1), vu(1), W(1), E, derivs, lambda)
     if(solution) then
             print*, 'The solution lambda:'
             print*, lambda(d,:,1)
             print*, lambda(s,:,1)
             print*, lambda(u,:,1)
     endif

     print*, 'W=',W
     print*, 'The derivatives:'
     print*, 'd/da:', derivs(1:n,1)
     print*, 'd/db:', derivs(n+1:2*n,1)
     print*, 'd/dbb:', derivs(2*n+1:3*n,1)
     print*, 'd/dvd:', derivs(3*n+1:4*n,1)
     print*, 'd/d(vs, vu, rd, rs, ru):', derivs(4*n+1:4*n+5,1)
     !try to improve the set of parameters
     do iter=1,niter
             print*,'===================='
             !step backward along the grad
             if (perturb_vector(1))  a=a-derivs(1:n,:)*step
             if (perturb_vector(2))  b=b-derivs(n+1:2*n,:)*step
             if (perturb_vector(3))  bb=bb-derivs(2*n+1:3*n,:)*step
             if (perturb_vector(4))  vd=vd-derivs(3*n+1:4*n,:)*step !this one is quite sensible
             if (perturb_vector(5))  vs=vs-derivs(4*n+1,:)*step
             if (perturb_vector(6))  vu=vu-derivs(4*n+2,:)*step
             if (perturb_vector(7))  rd=rd-derivs(4*n+3,:)*step
             if (perturb_vector(8))  rs=rs-derivs(4*n+4,:)*step
             if (perturb_vector(9))  ru=ru-derivs(4*n+5,:)*step
             !solve the direct problem for the perturbed parameters
             call solve_direct(m, zz, n, z, 1, a, b, bb, rd(1), rs(1), ru(1), vd, vs(1), vu(1), EdOASIM, EsOASIM, E, E_ave)
             if(solution) then
                     print*, 'The improved solution E:'
                     print*, E(d,:,1)
                     print*, E(s,:,1)
                     print*, E(u,:,1)
             endif
             print*, 'Functional: ', fac*0.5*(E(u,1,1) - Eu0_given)**2 !the functional value
             !the derivative of the functional wrt to the value to be compared with
             !the observed
             W = fac*(E(u,1,1) - Eu0_given)
             !solve the adjoint problem and get the derivatives
             call get_derivs(m, zz, n, z, a, b, bb, rd(1), rs(1), ru(1), vd, vs(1), vu(1), W(1), E, derivs, lambda)
             if(solution) then
                     print*, 'The solution lambda:'
                     print*, lambda(d,:,1)
                     print*, lambda(s,:,1)
                     print*, lambda(u,:,1)
             endif
             print*, 'The derivatives:'
             print*, 'd/da:', derivs(1:n,:)
             print*, 'd/db:', derivs(n+1:2*n,:)
             print*, 'd/dbb:', derivs(2*n+1:3*n,:)
             print*, 'd/dvd:', derivs(3*n+1:4*n,:)
             print*, 'd/d(vs, vu, rd, rs, ru):', derivs(4*n+1:4*n+5,:)
             print*, 'parameters:'
             print*, '    a         ', a
             print*, '    b         ', b
             print*, '    bb        ', bb
             print*, '    vd        ', vd
             print*, 'vs,vu,rd,rs,ru', vs,vu,rd,rs,ru
             print*, 'parameters wrt the original values:'
             print*, '    a         ', a/savepar(1:n,:)
             print*, '    b         ', b/savepar(n+1:2*n,:)
             print*, '    bb        ', bb/savepar(2*n+1:3*n,:)
             print*, '    vd        ', vd/savepar(3*n+1:4*n,:)
             print*, 'vs,vu,rd,rs,ru', [vs,vu,rd,rs,ru]/savepar(4*n+1:4*n+5,1)
             do lvl=1,n ! loop on levels
                 write(15,*) iter, z(lvl), z(lvl+1), a(lvl,1), b(lvl,1), bb(lvl,1), vd(lvl,1), vs, vu, rd, rs, ru
             enddo
             write(15,*) '####'
     enddo
     close(unit=15)
    
     print*, ''
     print*, '================================================='
     print*, '============ DESCENT ============================'
     print*, '================================================='
     pp = (100.*savepar)*1.10
     print*, 'parameters wrt the original values, %:', pp/savepar
     print*, 'func(par)=', functional(pp)
     print*, '================================================='
!    call FM34(functional, Gradient, pp, 1.0d-3, err)
     err = 0
     select case(err)
     case(0)
         print*, 'OK!!!'
     case(66)
         print*, 'Error', err, "'Number of iterations kmax exceeded'"
     case(67)
         print*, 'Error', err, "'Number of iterations cmax exceeded'"
     case(65)
         print*, 'Error', err, "'Possible instability'"
     case default
         print*, 'Unknown error', err
     end select
     print*, '================================================='
     print*, 'parameters wrt the original values, %:', pp/savepar
     print*, 'func(par)=', functional(pp)
     print*, '|grad(par)|=', gradient(pp,derivs)
     print*, 'grad(par)=', derivs
     print*, '================================================='

 enddo
 contains

     function functional(par, const) 
         double precision:: functional
         double precision,dimension(:,:),intent(in):: par
         double precision, dimension(:), intent(in), optional:: const     !any extra parameters
         a=par(1:n,:)
         b=par(n+1:2*n,:)
         bb=par(2*n+1:3*n,:)
         vd=par(3*n+1:4*n,:)
         vs=par(4*n+1,:)
         vu=par(4*n+2,:)
         rd=par(4*n+3,:)
         rs=par(4*n+4,:)
         ru=par(4*n+5,:)
         call solve_direct(m, zz, n, z, 1, a, b, bb, rd(1), rs(1), ru(1), vd, vs(1), vu(1), EdOASIM, EsOASIM, E, E_ave)
         functional =  fac*0.5*(E(u,1,1) - Eu0_given)**2 !the functional value
     end function functional

     function gradient(par,grad,const) 
         double precision:: gradient
         double precision,dimension(:,:),intent(in):: par
         double precision,dimension(:,:),intent(out):: grad
         double precision, dimension(:), intent(in), optional:: const     !any extra parameters
         a=par(1:n,:)
         b=par(n+1:2*n,:)
         bb=par(2*n+1:3*n,:)
         vd=par(3*n+1:4*n,:)
         vs=par(4*n+1,:)
         vu=par(4*n+2,:)
         rd=par(4*n+3,:)
         rs=par(4*n+4,:)
         ru=par(4*n+5,:)
         call solve_direct(m, zz, n, z, 1, a, b, bb, rd(1), rs(1), ru(1), vd, vs(1), vu(1), EdOASIM, EsOASIM, E, E_ave)
         W = fac*(E(u,1,1) - Eu0_given)
         call get_derivs(m, zz, n, z, a, b, bb, rd(1), rs(1), ru(1), vd, vs(1), vu(1), W(1), E, grad, lambda)
         gradient = sqrt(sum(grad*grad))
     end function gradient
end subroutine compute_3stream_adjoint

END MODULE adj_3stream

!program test
!use adj_3stream, only: test_3stream_adjoint
!call test_3stream_adjoint

!end program test

