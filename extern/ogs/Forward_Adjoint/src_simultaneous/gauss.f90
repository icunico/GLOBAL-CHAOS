
MODULE GAUSS
IMPLICIT NONE

PUBLIC:: SolveSLAEGauss, SolveOverdet

PRIVATE

double precision:: eps

CONTAINS

    !Solving a system of linear algebraic equations using the Gauss method of elimination.
    !Needs just det .ne. 0
    SUBROUTINE SolveSLAEGauss(n, m, A, x)
        integer, intent(in):: n, m                          !matrix is nxn, m right-hand sides
        double precision, dimension(n,n), intent(in):: A    !the matrix
        double precision, dimension(n,m), intent(inout):: x !in: the right-hand sides; out: the solution
        double precision, dimension(n,n):: B,D              !A=BD, B is down-triangle, D is up-triangle. To be computed.
        integer:: i,j,k
        !calculate B and D such that A=BD
        B = 0.0d0
        D = 0.0d0
        do i=1,n
            do k=1,i
                B(i,k) = A(i,k)
                do j=1,k-1
                    B(i,k)=B(i,k) - B(i,j)*D(j,k)
                enddo
            enddo
            D(i,i)=1.
            do k=i+1,n
                D(i,k) = A(i,k)
                do j=1,i-1
                    D(i,k)=D(i,k) - B(i,j)*D(j,k)
                enddo
                D(i,k) = D(i,k) / B(i,i)
            enddo
        enddo
        !Gauss elimination: having BDX=b, first solve By=b to get y=Dx
        do i=1,n
            x(i,:) = x(i,:) / B(i,i)
            B(i,:)=B(i,:) / B(i,i)
            do k=i+1,n
                x(k,:)= x(k,:)-B(k,i)*x(i,:)
                B(k,:)= B(k,:)-B(i,:)*B(k,i)
            enddo
        enddo
        !then solve Dx=y
        do i=n,1,-1
            x(i,:) = x(i,:) / D(i,i)
            D(i,:)=D(i,:) / D(i,i)
            do k=i-1,1,-1
                x(k,:)= x(k,:)-D(k,i)*x(i,:)
                D(k,:)= D(k,:)-D(i,:)*D(k,i)
            enddo
        enddo
        !here we have the solution in x. If anything is wrong, then det.eq.0
    END SUBROUTINE SolveSLAEGauss

    !Solves an overdetermined SLAE using the LSQ method
    SUBROUTINE SolveOverdet(ne, nx, A, b, x, misfit)
        !ne is number of equations (can be many)
        !mx is the number of unknowns
        !A is the matrix of size ne x mx,
        !b is the right-hand side of size ne
        !x is the solution of size nx
        integer, intent(in):: ne, nx                      !the matrix is ne x nx
        double precision, dimension(ne,nx),intent(in):: A !the matrix
        double precision, dimension(ne),intent(in):: b    !the right-hand side
        double precision, dimension(nx),intent(out):: x   !the solution to be got
        double precision, intent(out):: misfit            !the LSQ error. If large, then the equations contradict each other too much.
        double precision, dimension(nx,nx):: AA           !A transposed multiplied on A (A'*A)
        double precision, dimension(ne):: Ax              !A multiplied on x for comparing with b and evaluating the error
        integer:: i,j,k
        !multiply tr(A) on A and on x
        do i=1,nx
            do j=1,nx
                AA(i,j)=sum(A(:,i)*A(:,j))                !AA=A'*A
            enddo
            x(i)=sum(A(:,i)*b)                            !x=A'*b
        enddo
        call SolveSLAEGauss(nx,1,AA,x)                    !Solve the system A'*A*x=A'*b, which is the solution to the LSQ minimization problem
        !check
        Ax = 0.                                           !get Ax=A*x
        do i=1,ne
            Ax(i)=sum(A(i,:)*x)
        enddo
        misfit = norm2(Ax-b)                              !get the LSQ norm of Ax-b
    END SUBROUTINE SolveOverdet

END MODULE GAUSS

program test
    use gauss
    double precision, dimension(6,3):: A
    double precision, dimension(6,1):: b,bx
    double precision, dimension(3,1):: X
    double precision:: err
    A(1,:) = [1,2,3]
    A(2,:) = [0,1,2]
    A(3,:) = [-1,0,0]
    A(4,:) = [1,2,3]   *1.00001
    A(5,:) = [0,1,2]   *1.00002
    A(6,:) = [-1,0,0]  *0.99999
    b(1,:) = [6]       *1.00001
    b(2,:) = [3]       *1.00002
    b(3,:) = [-1]      *0.99999
    b(4,:) = [6]       /1.00001
    b(5,:) = [3]       /1.00002
    b(6,:) = [-1]      /0.99999
    bx=b                                                  !remember the RHS
    call SolveSLAEGauss(3,1,A(1:3,1:3),bx(1:3,1))         !test the Gauss solver, taking just three first equations
    print*, 'The solution is ', bx(1:3,1)
    err=0.
    do j=1,1
        do i=1,3
            err=err+abs(sum(A(i,:)*bx(1:3,j)) - B(i,j))
        enddo
    enddo
    print*, 'Error=',err
    call SolveOverdet(6,3,A,b,x,err)                      !Solve the full overditermined system
    print*, 'LSQ Sol:', x
    print*, 'LSQ error=',err
    err=0.
    do i=1,6
        err=err+abs(sum(A(i,:)*x(:,1)) - B(i,1))          !get the L1-norm error, abs(x-y) instead of (x-y)^2
    enddo
    print*, 'L1 error=',err

end program

