   program  prova_mult
   IMPLICIT NONE

   double precision :: A(3,3)
   double precision :: Ev(3),B(3),C
   double precision :: lambdav(3)

   A(1,:)= [1.0 , 0.0 , 0.0]
   A(2,:)= [1.0 , 1.0 , 0.0]
   A(3,:)= [1.0 , 1.0 , 1.0]

   Ev(:) = [1.0, 2.0 , 3.0]
   lambdav(:) = [0.0, 2.0 , 1.0]

   B=MATMUL(A,Ev)
   write(*,*) B
   C=DOT_PRODUCT(lambdav,B)
   write(*,*) C

end program

  
