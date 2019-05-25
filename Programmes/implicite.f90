PROGRAM  Eulerimplicit
  USE resolutionlu
  IMPLICIT NONE
  INTEGER,parameter::N=1000 
  REAL,parameter::Tmax=0.00005
  REAl::u,cfl,dx,dt
  INTEGER::i,m
  REAL,DIMENSION(0:n+1)::X
  REAL,DIMENSION(:),ALLOCATABLE::T
  REAl,DIMENSION(:,:),ALLOCATABLE::S
  dx=1./(N+1)
  X=(/ ( i*dx,i=0,N+1) /)
  u=1.
  cfl=0.9
  dt=cfl*(dx**2)/(2*u)
 
  CAll euler(u,cfl,Tmax,N,S,m)
  ALLOCATE(T(0:m))
  Do i=1,N
      PRINT*,X(i),S(i,m)
  END DO 

CONTAINS 
  SUBROUTINE euler(u,cfl,Tmax,N,S,m)
    IMPLICIT NONE 
    REAL,intent(in)::u,cfl,Tmax
    INTEGER,intent(in)::N
    INTEGER,INTENT(out)::m
    REAl,DIMENSION(:,:),ALLOCATABLE,intent(out)::S
    REAL,DIMENSION(0:N+1)::X,Y,b
    REAL,DIMENSION(:),ALLOCATABLE::T
    REAL,DIMENSION(0:N+1,0:N+1)::A,L,G
    INTEGER::i,j
    REAL::dt,dx
    dx=1./(N+1)
    dt=cfl*(dx**2)/(2*u)
    m=Tmax/dt
    X=(/(i*dx,i=0,n+1)/)
    ALLOCATE(T(0:m),S(0:N+1,0:m))
    T=(/(i*dt,i=0,m)/)
    DO i=0,N+1
       Do j=0,N+1
          IF ((i==j+1).or.(i==j-1)) THEN 
             A(i,j)=-cfl/2
          ELSE IF (i==j) THEN 
             A(i,j)=1+cfl
          ELSE 
             A(i,j)=0
          END IF
       END DO
    END DO
    Do j=0,N+1
       S(j,0)=f(real(X(j)))
    END DO 

    CALL decompo(A,L,G)
    Do i=0,m-1
       Do j=0,N+1
          b(j)=S(j,i)
       END DO
       X=resolution_inf(L,b)
       Y=resolution_sup(G,x)
       DO j=0,N+1
          S(j,i+1)=y(j)
       END DO
    END DO
  END SUBROUTINE euler
END PROGRAM Eulerimplicit
