module grad_descent

    public:: FM34, Cmax, Kmax

    integer:: Cmax=1000, Kmax=40

    private

interface
  double precision function minimized_fun(x,const)
      double precision, intent(in):: X(:)
      double precision, dimension(:), intent(in), optional:: const     !any extra parameters
  end function minimized_fun
  double precision function Grad_of_minimized(x, g, const)
      double precision, intent(in):: x(:)
      double precision, intent(out):: g(:)
      double precision, dimension(:), intent(in), optional:: const     !any extra parameters
  end function Grad_of_minimized
end interface

contains

subroutine FM34(Fun, Grad, X, Eps, error, low, high, verbose, const)
!Minimizing a function f(x1,x2,...,xn)
!by the Fletcher-Reeves method (conjugate gradients)
double precision, intent(in):: Eps     !precision of the x
double precision, intent(inout):: X(:) !initial assumption and the result, its size n is the reference for all arrays
double precision, intent(in), optional:: low(:), high(:) !lower and higher bounds for X, of size n
integer, intent(out):: error           !error code: 0=ok, 66=try higher Kmax, 65=possible instability, 67 - try higher Cmax
procedure(minimized_fun):: Fun         !the function to minimize, to be dble(1:n)->dble; see the interface above
procedure(Grad_of_minimized):: Grad    !the function to get the gradient (to be returned in the 2nd arg), returns the norm of the grad; see the interface.
logical, intent(in),optional:: verbose !if present and true, step info is printed
double precision, dimension(:), intent(in), optional:: const     !any extra parameters, any number of them
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
logical:: say
double precision, allocatable,dimension(:):: P, Q, D, G
double precision:: Z, FP, FQ, FR, GP, GQ, GR, G0, G1, G2, G3, GK
double precision:: QX, DD, HH, W, WW, ZZ, AK, EE
integer i, k, n, cc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
say = present(verbose).and.verbose
!begin
    !initial values
    EE=Eps; if (EE<Eps) EE=Eps 
    n=size(X); cc=0; Error=0
    allocate(P(n), Q(n), D(n), G(n))
    if(present(low)) X=max(X,low)
    if(present(high)) X=min(X,high)
400 continue
    !print *, '  The current values'
    !print 333, (i, X(i), i=1, N)
    do i=1,n; P(i)=X(i); enddo
    Z=Fun(X,const); FP=Z
    !print *, 'Iteration',cc,' Value of function',Z
    G0=Grad(X, G, const); GK=G0
    !Take the antigradient as the first direction
    do i=1,n; D(i)=-G(i); enddo
    k=0
600 continue
    GP=0.0
    do i=1,n; GP=GP+G(i)*D(i); enddo
    if (abs(GP)<Eps) goto 1100
    !Determine the first step
    QX=abs(2.0*FP/GP); if (QX>1.0) QX=1.0
    if (GP>0.0) then
        !Change the descent direction to the opposite one
        do i=1, N
            X(i)=P(i)-QX*D(i); P(i)=X(i)
        enddo
        if(present(low)) X=max(X,low)
        if(present(high)) X=min(X,high)
        Z=Fun(X,const); FP=Z !; print *, 'Possible instability?',' Value of function',Z
        G0=Grad(X, G, const); G1=G0
        k=k+1; if (k<2) goto 600
        !Abort
        Error=65; goto 1300  !possible instability
    endif
    HH=QX
700 continue
    !Find the next point Q
    do i=1, N
        Q(i)=P(i)+HH*D(i); X(i)=Q(i)
    enddo
    if(present(low)) X=max(X,low)
    if(present(high)) X=min(X,high)
    Z=Fun(X,const); FQ=Z
    G0=Grad(X, G, const); G2=G0
    GQ=0.0
    do i=1, N; GQ=GQ+G(i)*D(i); enddo
    if ((FP>FQ).and.(GQ<0.0)) then        
        if(say) print*,'    HH=',2.0*HH,' GP=',GP,' GQ=',GQ
        !Double the step to catch the min point to the uncertainty interval
        !do i=1,n; P(i)=Q(i); enddo
        !FP=FQ; GP=GQ; G1=G2
        HH=2.0*HH; goto 700
    endif
    k=0;
860 continue
    !Cubic interpolation
    ZZ=3.0*(FP-FQ)/HH+GP+GQ
    WW=ZZ*ZZ-GP*GQ; if (WW<0.0) WW=0.0
    W=sqrt(WW)
    DD=HH*(1.0-(GQ+W-ZZ)/(GQ-GP+2.0*W))
    do i=1,n; X(i)=P(i)+DD*D(i); enddo
    if(present(low)) X=max(X,low)
    if(present(high)) X=min(X,high)
    Z=Fun(X,const); FR=Z
    G0=Grad(X, G, const); G3=G0
    if (abs(G3)<EE) goto 1300
    !Evaluate the gradient in the new point
    GR=0.0
    do i=1,n; GR=GR+G(i)*D(i); enddo;
        if(say) then
            print 333, (i, X(i), i=1, N)
            print *, 'K=',k,' h=',HH,' fun=',FR,' grad=',G3
            print *, ' FP=',FP,' FQ=',FQ,' Z=FR=',FR
            print *, ' GP=',GP,' GQ=',GQ,' GR=',GR
        endif
    if ((Z<=FP).and.(Z<=FQ)) goto 1100
    !Repeat the cubic interpolation
    !The new interval would be either HH-DD or DD
    if (GR<0.0) then
        HH=HH-DD
        do i=1,n; P(i)=X(i); enddo
        FP=Z; GP=GR; G1=G0
        k=k+1; if (k<Kmax) goto 860
    else
        HH=DD
        do i=1,n; Q(i)=X(i); enddo
        FQ=Z; GQ=GR; G2=G0
        k=k+1; if (k<Kmax) goto 860
    endif
    !Abort
    Error=66; goto 1300 !66 means 'number of iterations exceeded'
1100 continue
    !Check the stop criteria
    do i=1,n
        if(abs(X(i)+HH)==abs(X(i))) goto 1300
    enddo
    cc=cc+1
    if (cc==Cmax) then
        error = 67 !67 means 'number of iterations exceeded'
        goto 1300
    endif
    if (mod(cc,n)==0) goto 400
    !Find the adjoint direction
    AK=G3*G3/(GK*GK)
    do i=1,n
        D(i)=-G(i)+AK*D(i); P(i)=X(i)
    enddo
    FP=FR; G1=G0; GK=G0; goto 600
1300 continue
    !The minimum found
    !print *, 'The minimum'
    !print 333, (i, X(i), i=1,n)
    !print *, 'The number of iterations =',cc,', value of function =',Z
    !print *, 'Error',Error
    deallocate(P,Q,D,G)
333 format((2X,'X[',I1,']=',F9.6))
end subroutine FM34

end module grad_descent

