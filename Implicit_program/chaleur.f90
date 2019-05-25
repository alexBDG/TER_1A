PROGRAM chaleur
  USE module
  IMPLICIT NONE

  INTEGER :: i, j
  REAL, PARAMETER :: kappa=1, lambda=100, cfl=0.9, tmax=0.1, w=0.045, h=10000
  INTEGER, PARAMETER :: imax=200, istop=1000

!  CALL Euler_implicit(kappa,lambda,imax,cfl,tmax,w,h,istop)

  CALL Euler_implicit_dynamique(kappa,lambda,imax,cfl,tmax,w,h,istop)

END PROGRAM chaleur
