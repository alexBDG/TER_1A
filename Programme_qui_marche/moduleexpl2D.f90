MODULE expl2D
  IMPLICIT NONE
CONTAINS
  SUBROUTINE Alexandre(dt, dx, dy, alpha, L, C, T, U)
    IMPLICIT NONE
    REAL, INTENT(IN)::dt, dx, dy, alpha
    REAL, DIMENSION(:,:,:), INTENT(OUT):: U
    REAL :: A, B, D
    INTEGER :: i,j,k, L, C, T
    B= aplha*(dt/(dx**2))
    D=alpha*(dt/(dy**2))
    A= 1-2*B-2*D
    L=size(U,1)
    C=size(U,2)
    U=size(U,3)
    DO i=1,L
       DO j=1,C
          U(i,j,1)=f(i*dx,j*dy)
       END DO
    END DO

    DO k=1,T-1
       U(1,1,k+1)=A*U(1,1,k)+B*U(2,1,k)+B*U(L,1,k)+D*U(1,2,k)+D*U(1,C,k)
       U(1,C,k+1)=A*U(1,C,k)+B*U(2,C,k)+B*U(L,C,k)+D*U(1,1,k)+D*U(1,C-1,k)
       DO i=2,L-1
          U(i,1,k+1)=A*U(i,1,k)+B*U(i+1,1,k)+B*U(i-1,1,k)+D*U(i,2,k)+D*U(i,C,k)
          U(i,C,k+1)=A*U(i,C,k)+B*U(i+1,C,k)+B*U(i-1,C,k)+D*U(i,1,k)+D*U(i,C-1,k)
       END DO
       U(L,1,k+1)=A*U(L,1,k)+B*U(1,1,k)+B*U(L-1,1,k)+D*U(L,2,k)+D*U(L,C,k)
       U(L,C,k+1)=A*U(L,C,k)+B*U(1,C,k)+B*U(L-1,C,k)+D*U(L,1,k)+D*U(L,C-1,k)
       DO j=2,C-1
          U(1,j,k+1)=A*U(1,j,k)+B*U(2,j,k)+B*U(L,j,k)+D*U(1,j+1,k)+D*U(1,j-1,k)
          DO i=2,L-1
             U(i,j,k+1)=A*U(i,j,k)+B*U(i+1,j,k)+B*U(i-1,j,k)+D*U(i,j+1,k)+D*U(i,j-1,k)
          END DO
          U(L,j,k+1)=A*U(L,j,k)+B*U(1,j,k)+B*U(L-1,j,k)+D*U(L,j+1,k)+D*U(L,j-1,k)
       END DO
    END DO

  END SUBROUTINE Alexandre

  FUNCTION f(x,y)
    IMPLICIT NONE
    REAL, INTENT(IN) :: x,y
    REAL :: f
    f=-1
  END FUNCTION f


END MODULE expl2D







