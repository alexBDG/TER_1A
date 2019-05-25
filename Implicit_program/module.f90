MODULE module
  IMPLICIT NONE

CONTAINS


!FONCTIONS UTILES POUR CHOLESKY
  SUBROUTINE matA(kappa,Dx,Dt,c,l)
    IMPLICIT NONE
    REAL, INTENT(IN) :: kappa, Dx, Dt
    REAL, DIMENSION(:), INTENT(OUT) :: c, l
    INTEGER :: i, j
    l = 1 + 2*kappa*Dt/(Dx**2)
    c = -kappa*Dt/(Dx**2)
  END SUBROUTINE matA

  SUBROUTINE cholesky(c,l)
    IMPLICIT NONE
    REAL, DIMENSION(:), INTENT(INOUT) :: c, l
    INTEGER :: N, i, k
    REAL :: s
    N=size(l)
    l(1)=SQRT(l(1))
    DO i=1,N-1
       c(i)=c(i)/l(i)
       l(i+1)=SQRT(l(i+1)-c(i)**2)
    END DO
  END SUBROUTINE cholesky

  SUBROUTINE solve_cholesky(c,l,b)
    IMPLICIT NONE
    REAL, DIMENSION(:), INTENT(INOUT) :: c, l, b
    INTEGER :: i, N
    N=size(l)
    b(1)=b(1)/l(1)
    DO i=2,N
       b(i)=(b(i)-c(i-1)*b(i-1))/l(i)
    END DO
    b(N)=b(N)/l(N)
    DO i=1,N-1
       b(N-i)=(b(N-i)-c(N-i)*b(N-i+1))/l(N-i)
    END DO
  END SUBROUTINE solve_cholesky

  SUBROUTINE second_membre(b,imax,Dt,lambda,cfl,i,istop,w,Dx,h)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: imax, istop, i
    REAL, INTENT(IN) :: Dt, lambda, cfl, w, h, Dx
    REAL, DIMENSION(imax), INTENT(OUT) :: b
    INTEGER :: j
    DO j=1,imax
       b(j)=b(j)+4*lambda*Dt*b(j)*(1-b(j)**2)+Dt*GaussDouble(istop,i,w,Dx,h)
    END DO
    b(1)=b(1)-cfl/2.
    b(imax)=b(imax)-cfl/2.
  END SUBROUTINE second_membre


! FONCTIONS INITIALES
  SUBROUTINE b0_1(b)
    IMPLICIT NONE
    REAL, DIMENSION(:), INTENT(OUT) :: b
    INTEGER :: i
    b=1.
    DO i=1,INT(size(b)/2)
       b(i)=-1
    END DO
  END SUBROUTINE b0_1

  SUBROUTINE b0_2(b)
    IMPLICIT NONE
    REAL, DIMENSION(:), INTENT(OUT) :: b
    INTEGER :: i
    b=-1.
  END SUBROUTINE b0_2


!FONCTIONS GAUSSIENNES
  FUNCTION Gauss(istop,i,w,Dx,h)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: i, istop
    REAL, INTENT(IN) :: Dx, w, h
    REAL :: Gauss
    IF (i<istop) THEN
       Gauss=h*exp(-(i*Dx-0.5)**2/(2*w**2))
    ELSE
       Gauss=0.
    END IF
  END FUNCTION Gauss

  FUNCTION GaussDouble(istop,i,w,Dx,h)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: i, istop
    REAL, INTENT(IN) :: Dx, w, h
    REAL :: GaussDouble
    IF (i<istop) THEN
       GaussDouble=h*exp(-(i*Dx-0.25)**2/(2*w**2))+h*exp(-(i*Dx-0.25)**2/(2*w**2))
    ELSE
       GaussDouble=0.
    END IF
  END FUNCTION GaussDouble


!NERF DU PROGRAMME : EULER IMPLICIT
  SUBROUTINE Euler_implicit(kappa,lambda,imax,cfl,tmax,w,h,istop)
    IMPLICIT NONE
    REAL, INTENT(IN) :: kappa, lambda, cfl, tmax, w, h
    INTEGER, INTENT(IN) :: imax, istop
    REAL, DIMENSION(:), ALLOCATABLE :: b, c, l
    REAL :: Dt, Dx
    INTEGER :: i

    Dx=1./(imax+1)
    Dt=cfl*kappa*Dx**2/(2*kappa)

    ALLOCATE(c(imax-1),l(imax),b(imax))

    CALL b0_2(b)

    DO i=1,10!INT(tmax/Dt)
       CALL matA(kappa,Dx,Dt,c,l)
       CALL cholesky(c,l)
       CALL second_membre(b,imax,Dt,lambda,cfl,i,istop,w,Dx,h)
       CALL solve_cholesky(c,l,b)
    END DO

    OPEN(unit=2,file='resultat_impl')
    WRITE(2,fmt=*)0,-1,tanh(-SQRT(2*lambda/kappa)*0.5)
    DO i=1,imax
       WRITE(2,'(F8.7,X,F8.5,X,F8.5)')i*Dx,b(i),tanh(SQRT(2*lambda/kappa)*(i*Dx-0.5))
    END DO
    WRITE(2,fmt=*)1,1,tanh(SQRT(2*lambda/kappa)*0.5)

    DEALLOCATE(c,l,b)

  END SUBROUTINE Euler_implicit


!ECRITURE DYNAMIQUE
  SUBROUTINE WriteGnuplot(x,y,name)
    REAL, DIMENSION(:), INTENT (INOUT) :: x,y
    CHARACTER(len=10), INTENT(IN) :: name
    INTEGER :: i
    OPEN(unit=1,file=name)
print*,name
    DO i=1,size(x)
       WRITE(1,fmt=*)x(i),y(i)
    END DO
    CLOSE(unit=1)
  END SUBROUTINE WriteGnuplot

function EntierToString(n) result(chaine)
  implicit none
  integer :: n, n0
  character(len=4) :: chaine
  n0 = n 
  chaine(4:4) = achar(48+mod(n0,10))
  n0 = n0/10
  chaine(3:3) = achar(48+mod(n0,10))
  n0 = n0/10
  chaine(2:2) = achar(48+mod(n0,10))
  n0 = n0/10
  chaine(1:1) = achar(48+mod(n0,10))
end function EntierToString


  !NERF DU PROGRAMME : EULER IMPLICIT VERSION DYNAMIQUE
  SUBROUTINE Euler_implicit_dynamique(kappa,lambda,imax,cfl,tmax,w,h,istop)
    IMPLICIT NONE
    REAL, INTENT(IN) :: kappa, lambda, cfl, tmax, w, h
    INTEGER, INTENT(IN) :: imax, istop
    REAL, DIMENSION(:), ALLOCATABLE :: b, c, l
    REAL, DIMENSION(imax) :: x
    REAL :: Dt, Dx
    INTEGER :: i

    integer :: Ndisplay
    character(len=6) entier
    Ndisplay = 1

    Dx=1./(imax+1)
    Dt=cfl*kappa*Dx**2/(2*kappa)

    x=(/(i*1./(imax+1),i=1,imax)/)

    ALLOCATE(c(imax-1),l(imax),b(imax))

    CALL b0_2(b)

    DO i=1,INT(tmax/Dt)
       CALL matA(kappa,Dx,Dt,c,l)
       CALL cholesky(c,l)
       CALL second_membre(b,imax,Dt,lambda,cfl,i,istop,w,Dx,h)
       CALL solve_cholesky(c,l,b)

       if (mod(i, Ndisplay) .eq. 0) then
!print*,i,"Un" // EntierToString(i/Ndisplay) // ".dat"
          call WriteGnuPlot(x, b, "Un" // EntierToString(i/Ndisplay) // ".dat")
       end if

    END DO

    DEALLOCATE(c,l,b)

  END SUBROUTINE Euler_implicit_dynamique

END MODULE module
