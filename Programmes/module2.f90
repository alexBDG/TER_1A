MODULE resolutionlu 
  IMPLICIT NONE 
CONTAINS 
  
  SUBROUTINE decompo(A,L,U)
    IMPLICIT NONE 
    REAL,DIMENSION(:,:),INTENT(in)::A
    REAL,DIMENSION(size(A(:,1)),size(A(1,:))),INTENT(out)::L,U
    REAL,DIMENSION(size(A(:,1)),size(A(1,:)))::M
    INTEGER::i,j,n,k
    n=size(A(1,:))
    DO i=1,n
       DO j=1,n
          M(i,j)=A(i,j)
       END DO
    END DO
    DO j=1,n-1
       DO i=j+1,n
          M(i,j)=M(i,j)/M(j,j)
       END DO
       DO i=j+1,n
          DO k=j+1,n
             M(i,k)=M(i,k)-M(i,j)*M(j,k)
          END DO
       END DO
    END DO
    DO i=1,n
       DO j=1,n
          IF (i<j) THEN 
             U(i,j)=M(i,j)
             L(i,j)=0
          ELSE IF   (i==j)  THEN
             U(i,j)=M(i,j)
             L(i,j)=1
          ELSE
             L(i,j)=M(i,j)
             U(i,j)=0
          END IF
       END DO
    END DO
  END SUBROUTINE decompo
  

  FUNCTION resolution_inf(L,b)result(x)
    IMPLICIT NONE
    REAL,DIMENSION(:,:),INTENT(in)::L
    REAL,DIMENSION(:),INTENT(in)::b
    REAL,DIMENSION(size(b))::x
    INTEGER::i,j,n
    REAL::m
    n=size(b)
    x(1)=b(1)
    DO i=2,n
       m=0
       Do j=1,i-1
          m=m+L(i,j)*x(j)
       END DO
       x(i)=b(i)-m
    END DO
  END FUNCTION resolution_inf
  

  FUNCTION resolution_sup(U,x)result(y)
    IMPLICIT NONE 
    REAL,DIMENSION(:,:),INTENT(in)::U
    REAL,DIMENSION(:),INTENT(in)::x
    REAL,DIMENSION(size(x))::y
    INTEGER::i,j,n
    REAL::m 
    n=size(x)
    y(n)=x(n)/U(n,n)
    DO i=n-1,1,-1
       m=0
       DO j=n,i+1,-1
          m=m+U(i,j)*y(j)
       END DO
       y(i)=(x(i)-m)/U(i,i)
    END DO
  END FUNCTION resolution_sup
  

  FUNCTION verif(A,x,z)result(d)
    IMPLICIT NONE 
    REAL,DIMENSION(:,:),INTENT(in)::A
    REAL,DIMENSION(:),INTENT(in)::x,z
    REAL,DIMENSION(size(z))::M
    LOGICAL::d
    INTEGER::i,j,n
    n=size(x)
    d=.true.
    DO i=1,n
       DO j=1,n
          M(i)=M(i)+A(i,j)*x(j)
       END DO
    END DO
    DO i=1,n
        If (M(i)/=z(i)) THEN 
           d=.false.
        END IF 
    END DO 
  END FUNCTION verif
  

  FUNCTION f(x)
    IMPLICIT NONE 
    REAL,INTENT(in)::x
    REAL::f
    f=6*exp(-4096*((x-0.5)**2))
  END FUNCTION f
  
 FUNCTION Q(x)
    IMPLICIT NONE 
    REAL,INTENT(in)::x
    REAL::Q
    Q=4*x*(1-x**2)
 END FUNCTION q

END MODULE  resolutionlu

