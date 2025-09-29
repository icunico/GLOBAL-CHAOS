module grad_descent

use par_mod

    implicit none

public:: FM34, Cmax, Kmax

integer:: Cmax=1000, Kmax=400

private

interface
    double precision function minimized_fun(mypar, x, const)                 !minimized functional must look like that
    use par_mod
    type(par),intent(inout) :: mypar
    double precision, intent(inout):: x(:)
    double precision, dimension(:), intent(in), optional:: const     !any extra parameters
    end function minimized_fun

    double precision function Grad_of_minimized(mypar, x, g, const)         !user-calculated gradient must look like that
    use par_mod
    type(par),intent(inout) :: mypar
    double precision, intent(inout):: x(:)
    double precision, intent(out) :: g(:)
    double precision, dimension(:), intent(in), optional:: const     !any extra parameters
    end function Grad_of_minimized
end interface

contains

subroutine FM34(mypar, functional, gradient, X, Eps, error, low, high, verbose, const)
use bioptimod_memory, only: no
use par_mod
!Minimizing a function f(x1,x2,...,xn)
!by the Fletcher-Reeves method (conjugate gradients)

double precision, intent(in):: Eps     !precision of the x
type(par),intent(inout) :: mypar
double precision, intent(inout):: X(:) !initial assumption and the result, its size n is the reference for all arrays
double precision, intent(in), optional:: low(:), high(:) !lower and higher bounds for X, of size n
integer, intent(out):: error           !error code: 0=ok, 66=try higher Kmax, 65=possible instability, 67 - try higher Cmax
procedure(minimized_fun):: functional         !the function to minimize, to be dble(1:n)->dble; see the interface above
procedure(Grad_of_minimized):: gradient    !the function to get the gradient (to be returned in the 2nd arg), returns the norm of the grad; see the interface.
logical, intent(in),optional:: verbose !if present and true, step info is printed
double precision, dimension(:), intent(in), optional:: const     !any extra parameters, any number of them
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
logical:: say
double precision, allocatable,dimension(:):: P, Q, D ,G
double precision:: Z, FP, FQ, FR, GP, GQ, GR, G0, G1, G2, G3, GK
double precision:: QX, DD, HH, W, WW, ZZ, AK, EE
integer i, k, n, cc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
say = present(verbose).and.verbose     !is true, then comment what is going on
!begin
!initial values
EE=Eps
if (EE<Eps) EE=Eps 
n=size(X); cc=0; Error=0
allocate(P(n), Q(n), D(n), G(n))
if(present(low)) X=max(X,low)
if(present(high)) X=min(X,high)
if(say) then
    print *, '  The current values:'
    print 333, (i, X(i), i=1, N)
endif
P(:)=X(:)
Z=functional(mypar,X,const); FP=Z
if(say) print *, 'Iteration',cc,' Value of function',Z
G0=gradient(mypar,X, G, const); GK=G0
!Take the antigradient as the first direction
D(:)=-G(:)
write(*,*) D
k=0
L600: do
    GP=sum(G(:)*D(:))
    I1100: if (.not.(abs(GP)<Eps)) then 
        !Determine the first step
        QX=abs(2.0*FP/GP); if (QX>1.0) QX=1.0
        if (GP>0.0) then
            !Change the descent direction to the opposite one
            do i=1, N
                X(i)=P(i)-QX*D(i); P(i)=X(i)
            enddo
            if(present(low)) X=max(X,low)
            if(present(high)) X=min(X,high)
            Z=functional(mypar,X,const); FP=Z !; print *, 'Possible instability?',' Value of function',Z
            G0=gradient(mypar,X, G, const); G1=G0
            k=k+1
            if(k<2) then 
                cycle
            else !Abort
                Error=65
                exit L600 !possible instability
            endif
        endif
        HH=QX
        L700:do
            !Find the next point Q
            Q(:)=P(:)+HH*D(:)
            X(:)=Q(:)
            if(present(low)) X=max(X,low)
            if(present(high)) X=min(X,high)
            Z=functional(mypar,X,const); FQ=Z
            G0=gradient(mypar,X, G, const); G2=G0
            GQ=sum(G(:)*D(:))
            if ((FP>FQ).and.(GQ<0.0)) then        
                if(say) print*,'    HH=',2.0*HH,' GP=',GP,' GQ=',GQ
                !Double the step to catch the min point to the uncertainty interval
                !do i=1,n; P(i)=Q(i); enddo
                !FP=FQ; GP=GQ; G1=G2
                HH=2.0*HH
            else
                exit L700
            endif
        enddo L700
        L1100:do k=0,kmax
            !Cubic interpolation
            ZZ=3.0*(FP-FQ)/HH+GP+GQ
            WW=ZZ*ZZ-GP*GQ
            if (WW<0.0) WW=0.0
            W=sqrt(WW)
            DD=HH*(1.0-(GQ+W-ZZ)/(GQ-GP+2.0*W))
            X(:)=P(:)+DD*D(:)
            if(present(low)) X=max(X,low)
            if(present(high)) X=min(X,high)
            Z=functional(mypar,X,const); FR=Z
            G0=gradient(mypar,X, G, const); G3=G0
            if (abs(G3)<EE) exit L600
            !Evaluate the gradient in the new point
            GR=sum(G(:)*D(:))
            if(say) then
                print 333, (i, X(i), i=1, N)
                print *, 'K=',k,' h=',HH,' fun=',FR,' grad=',G3
                print *, ' FP=',FP,' FQ=',FQ,' Z=FR=',FR
                print *, ' GP=',GP,' GQ=',GQ,' GR=',GR
            endif
            if ((Z<=FP).and.(Z<=FQ)) exit L1100 !k !was goto1100
            !Repeat the cubic interpolation
            !The new interval would be either HH-DD or DD
            if (GR<0.0) then
                HH=HH-DD
                P(:)=X(:)
                FP=Z; GP=GR; G1=G0
            else
                HH=DD
                Q(:)=X(:)
                FQ=Z; GQ=GR; G2=G0
            endif
        enddo L1100
        if(k>kmax) then !Abort if number of iterations exceeded
            Error=66
            exit L600 !66 means 'number of iterations exceeded'
        endif
    endif I1100
    !1100 was here
    !Check the stop criteria
    do i=1,n
        if(abs(X(i)+HH)==abs(X(i))) exit L600
    enddo
    cc=cc+1
    if (cc==Cmax) then
        error = 67 !67 means 'number of iterations exceeded'
        exit L600
    endif
    if (mod(cc,n)==0) then
        !print *, '  The current values'
        !print 333, (i, X(i), i=1, N)
        P(:)=X(:)
        Z=functional(mypar,X,const); FP=Z
        !print *, 'Iteration',cc,' Value of function',Z
        G0=gradient(mypar,X, G, const); GK=G0
        !Take the antigradient as the first direction
        D(:)=-G(:)
        k=0
        cycle
    endif
    !Find the adjoint direction
    AK=G3*G3/(GK*GK)
    D(:)=-G(:)+AK*D(:)
    P(:)=X(:)
    FP=FR; G1=G0; GK=G0; 
enddo L600
!The minimum found
!print *, 'The minimum'
!print 333, (i, X(i), i=1,n)
!print *, 'The number of iterations =',cc,', value of function =',Z
!print *, 'Error',Error
deallocate(P,Q,D)
333 format((2X,'X[',I1,']=',F9.6))
end subroutine FM34

end module grad_descent

