Program Str
 Implicit none
 INTEGER, PARAMETER:: IO = 12 ! input-output unit
 INTEGER NI, NJ, NITER, SMAX
 INTEGER I,J
 REAL L,H,U0,MU,Nu,R0,P0,h_s
 REAL dx,dy,CFL,EPS
 REAL,ALLOCATABLE :: X_Node(:,:),Y_Node(:,:)
 REAL,ALLOCATABLE :: U_n(:,:),V_n(:,:),P_n(:,:)
 REAL,ALLOCATABLE :: Re_x(:), cf_thi(:), Cf(:), Umax(:)
 write(*,*) 'Read input file'
 open(IO,FILE='InputPR.txt')
 read(IO,*) L
 read(IO,*) H
 read(IO,*) h_s
 read(IO,*) NI
 read(IO,*) NJ
 read(IO,*) EPS
 read(IO,*) SMAX
 read(IO,*) U0
 read(IO,*) MU
 read(IO,*) R0
 read(IO,*) P0
 CLOSE(IO)
 allocate(X_Node(NI,NJ)) ! mesh nodes X-coordinates
 allocate(Y_Node(NI,NJ)) ! mesh nodes Y-coordinates
!----------------- Node variables -----------------------------
 allocate(U_n(NI,NJ)) ! Velocity U
 allocate(V_n(NI,NJ)) ! Velocity V
 allocate(P_n(NI,NJ)) ! Pressure
 allocate(Re_x(NI))
 allocate(cf_thi(NI))
 allocate(Cf(NI))
 allocate(Umax(NI))
!----------------- Coordinate of nodes ------------------------
 dx=L/(NI-1)
 dy=H/(NJ-1)
 DO I=1,NI
     DO J=1,NJ
         X_Node(I,J)=(I-1)*dx
         Y_Node(I,J)=(J-1)*dy
     END DO
 END DO
!----------------- Parameters ------------------------
 NU=MU/R0
 write(*,*)'L= ', L, 'NI= ', NI, 'dx= ', dx
 write(*,*)'H= ', H, 'NJ= ', NJ, 'dy= ', dy
 write(*,*)'ReL= ', U0*L/NU
 !pause
!----------------- Initial fields -----------------------------
 DO I=1,NI
     DO J=1,NJ
         U_n(I,J)=U0
         V_n(I,J)=0.0
         P_n(I,J)=P0
     ENDDO
 ENDDO
!---------------- Solve Prandtl equations ---------------------
 write(*,*) 'Solve Prandtl equations'
 call Prandtl(NI,NJ,U_N,V_N,P_N,dx,dy,NU,U0,SMAX,EPS,h_s)
!----------------- Output data ------------------------------
 write(*,*) 'Output data'
 Open(IO,FILE='Resgid.plt')
 Call Output_Fields(IO,NI,NJ,X_Node,Y_Node,U_n,V_n,P_n)
 Close(IO)

 Open(2,FILE='tw.plt')
 write(2,*) 'VARIABLES = "Re_x", "Cf","cf_b","x"'
 do i=1, Ni
 !tw(i) = NU*(-3*U_n(I,1)+4*U_n(I,2)-U_n(I,3))/(2*dy)!??????????
 cf(i)=-(MU*((U_n(i,1) - U_n(i,2))/dy))/(0.5*R0*U0**2)
 Re_x(i) = U0*dx*(i-1)/NU
 cf_thi(i) = 0.664/sqrt(Re_X(i))
 IF (I>=2) then
 write(2,*) Re_x(i), cf(i),cf_thi(i), x_node(i,1)
 end if
 enddo
 Close(2)

OPEN(3,FILE='Umax.plt')
WRITE(3,*) 'VARIABLES = "x", "Umax"'
DO i=1, Ni
Umax(i)=MAXVAL(U_n(I,1:NJ))
WRITE(3,*) x_node(i,1), Umax(i)
ENDDO
CLOSE(3)


 END PROGRAM



!*****************************************************************************
 SUBROUTINE Prandtl(NI,NJ,U,V,P,dx,dy,NU,U0,SMAX,EPS,h_s)
 IMPLICIT NONE
!input
 INTEGER NI,NJ,NITTER,SMAX
 REAL dx,dy,NU,U0,EPS,h_s
 REAL,DIMENSION(NI,NJ)::U,V,P
!LOCAL
 INTEGER I,J,S,y_s
 REAL NORM_U,NORM_V
 REAL,DIMENSION(NJ)::UI0,UI1,VI0,VI1
 REAL,DIMENSION(NJ)::A,B,C,D
!INITIAL CONDITIONS FOR INLET BOUNADRY

 do j=1,NJ
     if (dy*j .lt. h_s) then
        U(1,j)=U0
     else
        U(1,j)=0.0
     endif
     V(1,J)=0.0
 END DO

 UI0=U(1,:)
 VI0=0.0
 UI1=UI0
 VI1=VI0

!SOLVING
 DO I=2,NI
     S=0
     NORM_U=10.0
     NORM_V=10.0
     UI0=U(I-1,:)
     VI0=V(I-1,:)
     UI1=UI0
     VI1=VI0
     DO WHILE (S.LE.SMAX .AND. MAX(NORM_U,NORM_V).GE.EPS)
         !U-VELOCITY
         A(1)=0.0
         B(1)=1.0
         C(1)=0.0
         D(1)=0.0

         DO J=2, NJ-1
            ! A(j)=-VI0(j-1)/(2.0*dy)-Nu/(dy**2.0)
            ! B(j)=UI0(j)/(dx)+2.0*Nu/(dy**2.0)
            ! C(j)=VI0(j+1)/(2.0*dy)-Nu/(dy**2.0)
            ! D(j)=(U(i-1,j)*U(i-1,j))/(dx)
              A(J) = (-2.0*NU*dx - VI0(J-1)*dx*dy)
              B(J) = (UI0(J)*2.0*dy*dy+4.0*NU*dx)
              C(J) = (VI0(J+1)*dx*dy - 2.0*NU*dx)
              D(J) = ((U(i-1,J)*U(I-1,J))*2.0*dy*dy)
         END DO

         A(NJ) = 0.0
         B(NJ) = 1.0
         C(NJ) = 0.0
         D(NJ) = 0.0

         CALL Progonka(NJ, A, B, C, D, UI1)
        !PRINT*,UI1(10),'U1',I,J

         !V-VELOCITY

         VI1(1)=0.0
         DO J=2,NJ !(UI1(J)-U(I-1,J)+UI1(J-1)-U(I-1,J-1))
            VI1(J)=VI1(j-1)-(0.5*DY/DX)*(UI1(J)-U(I-1,J)+UI1(J-1)-U(I-1,J-1))
         END DO
        !PRINT*,VI1(10),'V1'
         !NORMS
         NORM_U=MAXVAL(ABS(UI1(2:NJ)-UI0(2:NJ)))/MAXVAL(ABS(UI1(1:NJ)))
         NORM_V=MAXVAL(ABS(VI1(2:NJ)-VI0(2:NJ)))/MAXVAL(ABS(VI1(1:NJ)))
         UI0=UI1
         VI0=VI1
         S=S+1
     END DO
     WRITE(*,*) 'S=',S,'NORM_U',NORM_U,'NORM_V=',NORM_V,'U',UI1(17),'V',VI1(17),'Remain',NI-i
     U(I,:)=UI1(:)
     V(I,:)=VI1(:)
 END DO

 END SUBROUTINE

 Subroutine Progonka(n, a, b, c, d, u)
 integer :: i, n
 real(4) :: u(n), a(n), b(n), c(n), d(n), alf(n), bet(n)
 alf(1) = - c(1)/b(1)
 bet(1) = d(1)/b(1)
 do i = 2, n-1
 alf(i) = - c(i)/(b(i) + a(i) * alf(i-1))
 bet(i) = (d(i) - a(i) * bet(i-1)) / (b(i) + a(i) * alf(i-1))
 end do
 u(n)=(d(n) - a(n) * bet(n-1)) / (b(n) + a(n) * alf(n-1))
 do i = n - 1, 1 ,-1
 u(i) = alf(i) * u(i + 1) + bet(i)
 end do
 END SUBROUTINE
!*****************************************************************************
 SUBROUTINE Output_Fields(IO,NI,NJ,X,Y,U,V,P)
 IMPLICIT NONE
 INTEGER NI,NJ,IO
 REAL,DIMENSION(NI,NJ):: X,Y
 REAL,DIMENSION(NI,NJ):: U,V,P
 Write(IO,*) 'VARIABLES = "X", "Y", "U", "V", "P"'
 Write(IO,*) 'ZONE I=',NI,', J=',NJ, ', DATAPACKING=BLOCK'
 Write(IO,'(100E25.16)') X(1:NI,1:NJ)
 Write(IO,'(100E25.16)') Y(1:NI,1:NJ)
 Write(IO,'(100E25.16)') U(1:NI,1:NJ)
 Write(IO,'(100E25.16)') V(1:NI,1:NJ)
 Write(IO,'(100E25.16)') P(1:NI,1:NJ)
 END SUBROUTINE
!****************************************************************************
