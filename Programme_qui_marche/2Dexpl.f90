PROGRAM 2Dexpl
USE expl2D
IMPLICIT NONE

    REAL :: dt, dx, dy, alpha, cfl
    REAL, DIMENSION(:,:,:), ALLOCATABLE :: U
    INTEGER :: i,j,k, L, C, T

PRINT*,'L,C,cfl,T,alpha=?'
READ*, L, C, cfl, T, alpha
dx = 1/L
dy=1/C
dt=cfl*dx

CALL Alexandre(dt, dx, dy, alpha, L, C, T, U)

  OPEN(unit=1,file='file.dat',status='NEW')
  DO i=1,L
     Do j=1,C
        WRITE(1,fmt=*)U(i,j),U(i,j)
     END DO
  END DO
