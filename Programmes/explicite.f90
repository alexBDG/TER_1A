PROGRAM solution
  USE resolutionlu
  IMPLICIT NONE 
  INTEGER::i,j
  INTEGER, PARAMETER :: n = 100000, p = 100
  REAL::dt,dx
  REAL,DIMENSION(0:p,0:n)::U
  REAL,DIMENSION(0:p)::X
  REAL,DIMENSION(0:n)::T
  dt=1./n
  dx=1./p
  X=(/(i*1./p,i=0,p)/)
  T=(/(j*1./n,j=0,n)/)

#Fonction de dirac
!  DO j=1,p/2-1
!     U(j,0)=-1
!     U(j+p/2,0)=1
!  END DO

#Fonction sinus carré
!  DO j=1,p
!    U(j,0)=sin(3.1415*real(j))
!  END DO

#Fonction nulle
  DO j=1,p
     U(j,0)=0
  END DO

#Conditions aux limites
  Do i=0,n
     U(0,i)=0
     U(p,i)=0
  END DO

#Calcul des temps
  DO j=0,n-1
     Do i=1,p-1
         U(i,j+1)=U(i,j)+(U(i+1,j)-2*U(i,j)+U(i-1,j))*0.1-j*dt*f(U(i,j))+j*dt*h(U(i,j))
     END DO
     U(0,j)=0
     U(p,j)=0
  END DO 

#Affichage du temps voulu
  DO i=0,p-1
     PRINT*,X(i),U(i,10)
  END DO 

END PROGRAM solution
