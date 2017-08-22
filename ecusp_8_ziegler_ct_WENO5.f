      program main
      implicit none
      integer NX,NY,MP
      PARAMETER (NX=401,NY=401,MP=401)
      double precision:: ql(-5:MP+5,8), qr(-5:MP+5,8)

      double precision:: q(5,-5:NX+5,-5:NY+5,8) !Runge-Kutta 5,N=800,7 unkonwn
      double precision:: rho(-5:NX+5,-5:NY+5),
     $                   u(-5:NX+5,-5:NY+5),
     $                   v(-5:NX+5,-5:NY+5),
     $                   w(-5:NX+5,-5:NY+5),
     $                   B_x(-5:NX+5,-5:NY+5),
     $                   B_y(-5:NX+5,-5:NY+5),
     $                   B_z(-5:NX+5,-5:NY+5),
     $                   p(-5:NX+5,-5:NY+5)
      double precision:: fflux(-5:NX+5,-5:NY+5,8),
     $                   gflux(-5:NX+5,-5:NY+5,8),
     $                   source(-5:NX+5,-5:NY+5,8)

c---------Ziegler:
      double precision:: BxA(5,-5:NX+5,-5:NY+5),
     $                   ByA(5,-5:NX+5,-5:NY+5),
     $                   BzA(5,-5:NX+5,-5:NY+5),
     $                   Exx(-5:NX+5,-5:NY+5),
     $                   Exy(-5:NX+5,-5:NY+5),
     $                   Exz(-5:NX+5,-5:NY+5),
     $                   Eyx(-5:NX+5,-5:NY+5),
     $                   Eyy(-5:NX+5,-5:NY+5),
     $                   Eyz(-5:NX+5,-5:NY+5),
     $                   Exf(-5:NX+5,-5:NY+5),
     $                   Eyf(-5:NX+5,-5:NY+5),
     $                   Ezf(-5:NX+5,-5:NY+5),
     $                   B1L(-5:NX+5,-5:NY+5),
     $                   B1R(-5:NX+5,-5:NY+5),
     $                   B2L(-5:NX+5,-5:NY+5),
     $                   B2R(-5:NX+5,-5:NY+5),
     $                   ExL,ExR,EyL,EyR,EzL,EzR,ap,an,x05,y05

c---------END of Ziegler:
      double precision B2,b_2,v2,a2
      double precision gamma ,delta
      integer blen,IRK,NRK,Scheme,SCHEME1,cont,EXAMPLE,CT

c-----local variable
      integer i,j,i1,i2,k,isave
      double precision x,y,dx,dy,dt,CFL,temp,time,finaltime,pstar,
     $                ra,ra2,pi,Fr,Radius,xstart,ystart
      double precision du(8)
      double precision rhol,ul,vl,wl,bxl,byl,bzl,rhoel,pl,hl,
     $                 rhor,ur,vr,wr,bxr,byr,bzr,rhoer,pr,hr,
     $                 x_max,y_max,a,c_f,ua,va,wa
      DOUBLE PRECISION U05,V05,D05,P05,C05,
     $                 U15,V15,D15,P15,C15,
     $                 al,ar,c_al,c_ar,c_fl,c_fr,c_sl,c_sr,
     $                 C_MID,QML,QMR,SIGN_L,SIGN_R,
     $                 AFA_LP,AFA_RN,BTA_L,BTA_R,
     $                 DETAP,DETAN,QMLP,QMLN,
     $                 QMRP,QMRN,QM_MID,PHI,
     $                 QM_MID_P,QM_MID_N,C_P,C_N,
     $                 PLP,PRN,DLP,DRN,SL_P,SR_N
      double precision output(4)
      double precision temp1,temp2


      isave=0
      pi=2.0d0*asin(1.0d0)
      q=0.0

      output(1)=0.15
      output(2)=0.3
      output(3)=2.0
      output(4)=3.0
      blen = 4
      delta=1e-12

      example=3
      write(*,*)'Choose Scheme'
      write(*,*)'1------Orszag-Tang'
      write(*,*)'2------Kelvin-Hemhotz'
      write(*,*)'3------Rotor'
c      read(*,*)example

       CT=1    !    ''''''!£¿£¿£¿               
      write(*,*)'0-Normal ethod'
      write(*,*)'1-CT method'
c      read(*,*)CT

      scheme=5
      write(*,*)'Choose Scheme'
      write(*,*)'0------Ziegler TVD'
      write(*,*)'1------1-Upwind'
      write(*,*)'3------3-WENO'
      write(*,*)'5------5-WENO'
c      read(*,*)scheme

      scheme1=33
      write(*,*)'Choose Riemman Solver'
      write(*,*)'0  ------Ziegler TVD'
      write(*,*)'3  ------Zha-3'
      write(*,*)'33 ------Van Leer (RUN)'
      write(*,*)'333------different speed(fixed)'
      write(*,*)'334------different speed(testing)'
c      read(*,*)scheme1

      CFL=0.05
      finaltime=0.15 ! (Orszag=3, Rotor=0.15)
      NRK=3
      write(*,*)'INPUT CFL,T,NRK'
c      read(*,*)CFL,finaltime,NRK


      cont=0
      write(*,*)'New-0/Continue-1'
c      read(*,*)cont

!      open(13,file='datain.txt')
!      read(13,*)example
!      read(13,*)CT
!      read(13,*)scheme
!      read(13,*)scheme1
!      read(13,*)CFL,finaltime,NRK
!      read(13,*)cont
      

      open(10,file='input.bak')
      write(10,*)'Example:',example
      write(10,*)'CT method:',CT
      write(10,*)'Scheme:',scheme
      write(10,*)'Solver:',scheme1
      write(10,*)'CFL,finaltime,NRK=',CFL,finaltime,NRK
      write(10,*)'Continue',cont
      close(10)


c------------------------------------------
      SELECT CASE (EXAMPLE)

      CASE (1)
      gamma = 5.0d0/3.0d0

      dx=2.0*pi/(NX-1.0)
      dy=2.0*pi/(NY-1.0)
      xstart=0.0
      ystart=0.0

c--------------2D Rotor
      CASE (3)
      gamma = 1.4

      dx=1.0/(NX-1.0)
      dy=1.0/(NY-1.0)
      xstart=0.5
      ystart=0.5

      END SELECT
c------------------------------------------

      if (cont.eq.1) then
      open(10,file='final_t.dat')
      read(10,*)time
      close(10)

      open(10,file='rho.dat')
      read(10,*)
      read(10,*)
      read(10,*)
      do i=1,NX
         do j=1,NY
            read(10,*)x,y,rho(i,j),u(i,j),v(i,j),w(i,j),
     $                 b_x(i,j),b_y(i,j),b_z(i,j),p(i,j)
c-------------
            q(1,i,j,1)=rho(i,j)
            q(1,i,j,2)=rho(i,j)*u(i,j)
            q(1,i,j,3)=rho(i,j)*v(i,j)
            q(1,i,j,4)=rho(i,j)*w(i,j)
            q(1,i,j,5)=B_x(i,j)
            q(1,i,j,6)=B_y(i,j)
            q(1,i,j,7)=B_z(i,j)
            q(1,i,j,8)=0.5*rho(i,j)*(u(i,j)**2+v(i,j)**2+w(i,j)**2)+
     $            0.5*(B_x(i,j)**2+B_y(i,j)**2+B_z(i,j)**2)+
     $            p(i,j)/(gamma-1.)
         end do
      end do
      close(10)

c-------------
      else
      time=0.0
c-----initial condition
      do i=1-blen,NX+blen
         x=(i-1.0)*dx-xstart
         x05=(i-1.0-0.5)*dx-xstart
         do j=1-blen,NY+blen
            y=(j-1.0)*dy-ystart
            y05=(j-1.0-0.5)*dy-ystart

            
            SELECT CASE (EXAMPLE)
            CASE(1)
c-------------Orszag-Tang
            rho(i,j)=gamma**2
            u(i,j)=-sin(y)
            v(i,j)=sin(x)
            w(i,j)=0.0
            B_x(i,j)=-sin(y)
            B_y(i,j)=sin(2.0*x)
            B_z(i,j)=0.0
            p(i,j)=gamma
c-------------2D Rotor    Toth-2000_JCP
            CASE(3)

            radius=sqrt(x*x+y*y)

            if (radius.lt.0.1) then
               fr=1.0
               radius=0.1
            elseif (radius.gt.0.115) then
               fr=0.0
            else
               fr=200.0/3.0*(0.115-radius)
            end if                        

            rho(i,j)=1.0+9.*fr
            p(i,j)=1.0
            u(i,j)=-2.*fr*y/radius
            v(i,j)=2.*fr*x/radius


            w(i,j)=0.0
            B_x(i,j)=5.0/sqrt(4.0*pi)  !!!!!!!!!!!!!!!!!!!!!!!
            B_y(i,j)=0.0
            B_z(i,j)=0.0

c-----BEGIN---:http://flash.uchicago.edu/website/codesupport
c            fr=(0.115-radius)/(radius-0.1) 
c            if (radius.lt.0.1) then
c                fr=1.0 !!!!!!!!!!!!!!!!!!
c               rho(i,j)=10.
c               u(i,j)=-fr*y/0.1
c               v(i,j)= fr*x/0.1
c            elseif (radius.gt.0.115) then
c               rho(i,j)=1.0
c               u(i,j)=0.0
c               v(i,j)=0.0
c            else
c               rho(i,j)=1.0+9.0*fr
c               u(i,j)=-fr*y/radius
c               v(i,j)= fr*x/radius
c            end if
c            p(i,j)=1.0
c            w(i,j)=0.0
c            B_x(i,j)=5.0/sqrt(4.0*pi)   !!!!!!!!!!!!!!!!!!!!!!!
c            B_y(i,j)=0.0
c            B_z(i,j)=0.0

c-------END---:http://flash.uchicago.edu/website/codesupport

c-----BEGIN---:Mishra & Tadmor 2009 ETH (downloaded)
c            if (radius.lt.0.1) then
c               fr=1.0
c            elseif (radius.gt.0.115) then
c               fr=0.0
c            else
c               fr=200.0/3.0*(0.115-radius)
c            end if
c            rho(i,j)=1.0+9.*fr
cc            p(i,j)= 1.0
cc            u(i,j)=-2.*fr*y/radius
cc            v(i,j)= 2.*fr*x/radius
c            w(i,j)=0.0
c            B_x(i,j)=5.0/sqrt(4.0*pi)  !!!!!!!!!!!!!!!!!!!!!!!
c            B_y(i,j)=0.0
c            B_z(i,j)=0.0
c            if (radius.lt.0.1) then
c               u(i,j)= (10.*y-5.0)
c               v(i,j)=-(10.*x-5.0)
c            elseif (radius.gt.0.115) then
c               u(i,j)=0.0
c               v(i,j)=0.0
c            else
c               u(i,j)= (10.*y-5.0)*fr
c               v(i,j)=-(10.*x-5.0)*fr
c            end if
c            p(i,j)=0.5
c-----END---:Mishra & Tadmor 2009 ETH 

c-------------initial values for 2D Rotor Problem

            END SELECT


            q(1,i,j,1)=rho(i,j)
            q(1,i,j,2)=rho(i,j)*u(i,j)
            q(1,i,j,3)=rho(i,j)*v(i,j)
            q(1,i,j,4)=rho(i,j)*w(i,j)
            q(1,i,j,5)=B_x(i,j)
            q(1,i,j,6)=B_y(i,j)
            q(1,i,j,7)=B_z(i,j)
            q(1,i,j,8)=0.5*rho(i,j)*(u(i,j)**2+v(i,j)**2+w(i,j)**2)+
     $            0.5*(B_x(i,j)**2+B_y(i,j)**2+B_z(i,j)**2)+
     $            p(i,j)/(gamma-1.)

         end do
      end do
      end if
      
      
      
c-----initial condition FOR CONSTRAINED TRANSPORT METHOD-------------------
c=============== B is at half point===================================
      if (CT.EQ.1) then
      do i=1-blen,NX+blen
         do j=1-blen,NY+blen
            
            SELECT CASE (EXAMPLE)
c-------------2D Rotor    Toth-2000_JCP
            CASE(3)

            BxA(1,i,j)=5.0/sqrt(4.0*pi)  !!!!!!!!!!!!!!!!!!!!!!!
            ByA(1,i,j)=0.0
            BzA(1,i,j)=0.0


c-------------initial values for 2D Rotor Problem

            END SELECT

            q(1,i,j,5)=BxA(1,i,j)
            q(1,i,j,6)=ByA(1,i,j)
            q(1,i,j,7)=BzA(1,i,j)
            q(1,i,j,8)=0.5*rho(i,j)*(u(i,j)**2+v(i,j)**2+w(i,j)**2)+
     $            0.5*(B_x(i,j)**2+B_y(i,j)**2+B_z(i,j)**2)+
     $            p(i,j)/(gamma-1.)

         end do
      end do
      end if !if (CT.EQ.1)
            
c--------end of initial condition
c--------boundary for Orszag-Tang-----------------------------------
            SELECT CASE (EXAMPLE)
            CASE(1)
      do j=1-blen,NY+blen
         do k=1-blen,0
            do i1=1,8
               q(1,k,j,i1)=q(1,NX+k-1,j,i1)
               q(1,NX-k+1,j,i1)=q(1,2-k,j,i1)
            end do
          end do
       end do

      do i=1-blen,NX+blen
         do k=1-blen,0
            do i1=1,8
               q(1,i,k,i1)=q(1,i,NY+k-1,i1)
               q(1,i,NY-k+1,i1)=q(1,i,2-k,i1)
            end do
          end do
       end do

c--------end of Boundary


c--------boundary for 2D Rotor Problem---------------------------------
c--------Boundary: free outflow 
            CASE(3)
      do j=1-blen,NY+blen
         do k=1-blen,0
            do i1=1,8
               q(1,k,j,i1)=q(1,1,j,i1)
               q(1,NX-k+1,j,i1)=q(1,NX,j,i1)
            end do
          end do
       end do
      do i=1-blen,NX+blen
         do k=1-blen,0
            do i1=1,8
               q(1,i,k,i1)=q(1,i,1,i1) 
               q(1,i,NY-k+1,i1)=q(1,i,NY,i1) 
            end do
          end do
       end do

            END SELECT
c--------------------------------------------------------------
c------------Due to BxA,ByA,BzA are constants------------------
      DO I=1-blen,NX+blen
         DO J=1-blen,NY+blen
            BxA(1,I,J)=q(1,I,J,5)
            ByA(1,I,J)=q(1,I,J,6)
            BzA(1,I,J)=q(1,I,J,7)
         END DO
      END DO

c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

c
c---------main loop---------------------------------
 100  x_max=0.0
      y_max=0.0
      do i=1,NX
         do j=1,NY
            rho(i,j)=q(1,i,j,1)
            u(i,j)=q(1,i,j,2)/rho(i,j)
            v(i,j)=q(1,i,j,3)/rho(i,j)
            w(i,j)=q(1,i,j,4)/rho(i,j)
            B_x(i,j)=q(1,i,j,5)
            B_y(i,j)=q(1,i,j,6)
            B_z(i,j)=q(1,i,j,7)
            B2=B_x(i,j)**2+B_y(i,j)**2+B_z(i,j)**2
            v2=u(i,j)**2+v(i,j)**2+w(i,j)**2
            p(i,j)=(gamma-1.)*(q(1,i,j,8)-0.5*(rho(i,j)*v2+B2))

c-----------------
           if (rho(i,j).le.0.0) then
              write(*,*)'time_rho',i,j,rho(i,j)
!              read(*,*)k
!              rho(i,j)=0.01
               rho(i,j)=-rho(i,j)     
            rho(i,j)=rho(i-1,j-1)
            q(1,i,j,8)=0.5*rho(i,j)*(u(i,j)**2+v(i,j)**2+w(i,j)**2)+
     $            0.5*(B_x(i,j)**2+B_y(i,j)**2+B_z(i,j)**2)+
     $            p(i,j)/(gamma-1.)

           end if
           if (p(i,j).le.0.0) then
              write(*,*)'time_p',i,j,p(i,j)
!              read(*,*)k
!              p(i,j)=0.01
                p(i,j)=-p(i,j)
!                p(i,j)=p(i-1,j-1)
            q(1,i,j,8)=0.5*rho(i,j)*(u(i,j)**2+v(i,j)**2+w(i,j)**2)+
     $            0.5*(B_x(i,j)**2+B_y(i,j)**2+B_z(i,j)**2)+
     $            p(i,j)/(gamma-1.)

           end if
           if (B2.lt.0.0) then
              write(*,*)'time_B',i,j,B_x(i,j),B_y(i,j),B_z(i,j)
!              read(*,*)k
           end if
           if (v2.lt.0.0) then
              write(*,*)'time_v',i,j,u(i,j),v(i,j),w(i,j)
!              read(*,*)k
           end if
c-----------------

            a=sqrt(gamma*p(i,j)/rho(i,j))


            b_2=B2/rho(i,j)
            a2=a**2
            temp=(a2+b_2)**2-4.0*a2*B_x(i,j)**2/rho(i,j)
            if (temp.gt.0.0) then
               c_f=sqrt(0.5*(a2+b_2+sqrt(temp)))
            else
               write(*,*)'time1,temp=',temp
            !  read(*,*)k
               if (a2+b_2.le.0.0) then
                   write(*,*)'time2,a2+b_2=',a2+b_2
              !     read(*,*)k
               end if
               c_f=sqrt(0.5*(a2+b_2))
            end if
c            c_f=sqrt(0.5*(a2+b_2+
c     $          sqrt((a2+b_2)**2-4.0*a2*B_x(i,j)**2/rho(i,j))))
            x_max=max(x_max,abs(u(i,j))+c_f)
            temp=(a2+b_2)**2-4.0*a2*B_y(i,j)**2/rho(i,j)
            if (temp.gt.0.0) then
               c_f=sqrt(0.5*(a2+b_2+sqrt(temp)))
            else
           !   read(*,*)k
               write(*,*)'time_y1,temp=',temp
               if (a2+b_2.le.0.0) then
            !  read(*,*)k
                   write(*,*)'time_y2,a2+b_2=',a2+b_2
               end if
               c_f=sqrt(0.5*(a2+b_2))
            end if
c            c_f=sqrt(0.5*(a2+b_2+
c     $          sqrt((a2+b_2)**2-4.0*a2*B_y(i,j)**2/rho(i,j))))
            y_max=max(y_max,abs(v(i,j))+c_f)
         end do
      end do
      write(*,*)'ttttime',x_max,y_max
      dt=CFL/(x_max/dx+y_max/dy)
c      dt=CFL*min(dx/x_max,dy/y_max)
      time=time+dt
      write(*,*)time      
c-----------Runge-Kutta Mathod----------------------
      DO IRK=1,NRK

c      do i=1-blen,NX+blen
c         do j=1-blen,NY+blen
c            rho(i,j)=q(1,i,j,1)
c            u(i,j)=q(1,i,j,2)/rho(i,j)
c            v(i,j)=q(1,i,j,3)/rho(i,j)
c            w(i,j)=q(1,i,j,4)/rho(i,j)
c            B_x(i,j)=q(1,i,j,5)
c            B_y(i,j)=q(1,i,j,6)
c            B_z(i,j)=q(1,i,j,7)
c            B2=B_x(i,j)**2+B_y(i,j)**2+B_z(i,j)**2
c            v2=u(i,j)**2+v(i,j)**2+w(i,j)**2
c            p(i,j)=(gamma-1.)*(q(1,i,j,8)-0.5*(rho(i,j)*v2+B2))
c         end do
c      end do



c---------------X-DIRECTION-------------------------------------
      IF (CT.eq.1) then
c      DO J=0,NY+1
      DO J=-1,NY+1
         k=6  !By
         call reconstruct(scheme,NX,ByA(1,:,j),B1L(:,j),B1R(:,j))
         k=7  !Bz
         call reconstruct(scheme,NX,BzA(1,:,j),B2L(:,j),B2R(:,j))
      END DO
      END IF

c      do j=1,NY
      do j=-1,NY+1

c--------reconstruction of Magnetic field, different from Y-DIRECTION,
c        B_x^E,B_y^E,B_z^E [E===>L]
c        B_x^W,B_y^W,B_z^W [W===>R(i-1)]
c        BxA--->Bx average at (i+1/2,j,k)
c        ByA--->By average at (i,j+1/2,k)
c        BzA--->Bz average at (i,j,k+1/2)

c         DO I=0,NX+1
c            QL(i,5)=BxA(i,j)                  ! B^E(i,j,k)
c            QL(i,6)=0.5*(ByA(i,j)+ByA(i,j-1)+...)
c            QL(i,7)=0.5*(BzA(i,j)+BzA(i,j-1)+...)

c            QR(i,5)=BxA(i,j)              ! B^W(i,j,k)
c            QR(i,6)=0.5*(ByA(i+1,j)+ByA(i+1,j-1)-...)
c            QR(i,7)=0.5*(BzA(i+1,j)+BzA(i+1,j-1)-...)
c         END DO


      IF (CT.eq.1) then
      do k=1,4
         call reconstruct(scheme,NX,q(1,:,j,k),ql(:,k),qr(:,k))
      end do

         k=5  !Bx
         DO I=-1,NX+1
            QL(i,5)=BxA(1,i,j)                  ! B^E(i,j,k)
            QR(i,5)=BxA(1,i,j)              ! B^W(i,j,k)
         END DO
         !k=6,7
c         DO I=0,NX+1
         DO I=-1,NX+1
            QL(i,6)=0.5*(B1L(i,j)+B1L(i,j-1))
            QR(i,6)=0.5*(B1R(i,j)+B1R(i,j-1))
            QL(i,7)=0.5*(B2L(i,j)+B2L(i,j-1))
            QR(i,7)=0.5*(B2R(i,j)+B2R(i,j-1))

         END DO

c---------------new CT
        do k=5,7
           call reconstruct(scheme,NX,q(1,:,j,k),ql(:,k),qr(:,k))
        end do
c---------------new CT
         
         k=8
         call reconstruct(scheme,NX,q(1,:,j,k),ql(:,k),qr(:,k))
        ELSE ! IF (CT.eq.1) then

c---------reconstruction--------------
      do k=1,8 
         call reconstruct(scheme,NX,q(1,:,j,k),ql(:,k),qr(:,k))
      end do
       END IF ! IF (CT.eq.1) then
c--------------calculate eigenvectors---------------
c      do i=0,NX+1
      do i=-1,NX+1

         rhol=ql(i,1)
         ul=ql(i,2)/rhol
         vl=ql(i,3)/rhol
         wl=ql(i,4)/rhol
         bxl=ql(i,5)
         byl=ql(i,6)
         bzl=ql(i,7)
         rhoel=ql(i,8)
         B2=Bxl**2+Byl**2+Bzl**2
         v2=ul**2+vl**2+wl**2
         pl=(gamma-1.)*(rhoel-0.5*(rhol*v2+B2))
         hl=(rhoel+pl+0.5*B2)/rhol
         P05=pl+B2/2.0

         b_2=B2/rhol
         al=sqrt(gamma*pl/rhol)         
         a2=al**2
         c_al=abs(Bxl/sqrt(rhol))
         temp1=(a2+b_2)**2-4.0*a2*Bxl**2/rhol

c-----------------
           if (pl.le.0. .or.rhol.le.0.0.or.temp1.le.0.) then
c              write(*,*)'reconstruct_xL',i,j,pl,rhol,temp1
c              write(*,*)temp1

c-----------------
         rhol=q(1,i,j,1)
         ul=q(1,i,j,2)/rhol
         vl=q(1,i,j,3)/rhol
         wl=q(1,i,j,4)/rhol
         bxl=q(1,i,j,5)
         byl=0.5*(q(1,i,j,6)+q(1,i,j-1,6))
         bzl=0.5*(q(1,i,j,7)+q(1,i,j-1,7))

c         byl=0.5*(q(1,i,j,6)+q(1,i,j+1,6))
c         bzl=0.5*(q(1,i,j,7)+q(1,i,j+1,7))

         rhoel=q(1,i,j,8)
         B2=Bxl**2+Byl**2+Bzl**2
         v2=ul**2+vl**2+wl**2
         pl=(gamma-1.)*(rhoel-0.5*(rhol*v2+B2))
         hl=(rhoel+pl+0.5*B2)/rhol
         P05=pl+B2/2.0

         b_2=B2/rhol
         al=sqrt(gamma*pl/rhol)         
         a2=al**2
         c_al=abs(Bxl/sqrt(rhol))
         temp1=(a2+b_2)**2-4.0*a2*Bxl**2/rhol
c         write(*,*)'----------XL',temp1,a2,b_2,Bxl,rhol,pl
c              read(*,*)k
c----------------
           end if
c-----------------

c         al=sqrt(gamma*pl/rhol)         
c         a2=al**2
c         c_al=abs(Bxl/sqrt(rhol))
c         temp1=(a2+b_2)**2-4.0*a2*Bxl**2/rhol
c         temp2=0.5*(a2+b_2-sqrt((a2+b_2)**2-4.0*a2*Bxl**2/rhol))
c         if (temp1.le.0.) then
c            write(*,*)'temp1-xl',i,j,temp1
c         end if
c         if (temp2.le.0.) then
c            write(*,*)'temp2-xl',i,j,temp2
c         end if
         c_fl=sqrt(0.5*(a2+b_2+sqrt((a2+b_2)**2-4.0*a2*Bxl**2/rhol)))
         c_sl=sqrt(0.5*(a2+b_2-sqrt((a2+b_2)**2-4.0*a2*Bxl**2/rhol)))

         
         rhor=qr(i,1)
         ur=qr(i,2)/rhor
         vr=qr(i,3)/rhor
         wr=qr(i,4)/rhor
         bxr=qr(i,5)
         byr=qr(i,6)
         bzr=qr(i,7)
         rhoer=qr(i,8)
         B2=Bxr**2+Byr**2+Bzr**2
         v2=ur**2+vr**2+wr**2
         pr=(gamma-1.)*(rhoer-0.5*(rhor*v2+B2))
         hr=(rhoer+pr+0.5*B2)/rhor
         P15=pr+B2/2.0

         b_2=B2/rhor
         ar=sqrt(gamma*pr/rhor)     
         a2=ar**2
         c_ar=abs(Bxr/sqrt(rhor))
         temp1=(a2+b_2)**2-4.0*a2*Bxr**2/rhor

c-----------------
           if (pr.le.0. .or.rhor.le.0.0.or.temp1.le.0.0) then
c              write(*,*)'reconstruct_xR',i,j,pr,rhor,temp1
c              write(*,*)temp1
c------------------
         rhor=q(1,i+1,j,1)
         ur=q(1,i+1,j,2)/rhor
         vr=q(1,i+1,j,3)/rhor
         wr=q(1,i+1,j,4)/rhor
         bxr=q(1,i,j,5)
         byr=0.5*(q(1,i+1,j,6)+q(1,i+1,j-1,6))
         bzr=0.5*(q(1,i+1,j,7)+q(1,i+1,j-1,7))

c         byr=0.5*(q(1,i+1,j,6)+q(1,i+1,j+1,6))
c         bzr=0.5*(q(1,i+1,j,7)+q(1,i+1,j+1,7))


         rhoer=q(1,i+1,j,8)
         B2=Bxr**2+Byr**2+Bzr**2
         v2=ur**2+vr**2+wr**2
         pr=(gamma-1.)*(rhoer-0.5*(rhor*v2+B2))
         hr=(rhoer+pr+0.5*B2)/rhor
         P15=pr+B2/2.0

         b_2=B2/rhor
         ar=sqrt(gamma*pr/rhor)     
         a2=ar**2
         c_ar=abs(Bxr/sqrt(rhor))
         temp1=(a2+b_2)**2-4.0*a2*Bxr**2/rhor
c         write(*,*)'----------XR',temp1,a2,b_2,Bxr,rhor,pr
c         write(*,*)'----------1',rhoer,0.5*(rhor*v2+B2)
c         write(*,*)'??????????1',ur,vr,wr,bxr,byr,bzr

c------------------              
c              read(*,*)k
           end if
c-----------------
    
c         ar=sqrt(gamma*pr/rhor)     
c         a2=ar**2
c         c_ar=abs(Bxr/sqrt(rhor))
c         temp1=(a2+b_2)**2-4.0*a2*Bxr**2/rhor
c         temp2=0.5*(a2+b_2-sqrt((a2+b_2)**2-4.0*a2*Bxr**2/rhor))
c         if (temp1.le.0.) then
c            write(*,*)'temp1-xr',i,j,temp1
c         end if
c         if (temp2.le.0.) then
c            write(*,*)'temp2-xr',i,j,temp2
c         end if
         c_fr=sqrt(0.5*(a2+b_2+sqrt((a2+b_2)**2-4.0*a2*Bxr**2/rhor)))
         c_sr=sqrt(0.5*(a2+b_2-sqrt((a2+b_2)**2-4.0*a2*Bxr**2/rhor)))



         d05=rhol
         d15=rhor
         u05=ul
         u15=ur

c--------Ziegler----------------------
c   ap: A^(+)(i+1/2,j,k)
c   an: A^(-)(i+1/2,j,k)
c   Electric field     Exx(i,j): E^x_x(i+1/2,j)
c   Electric field     Exy(i,j): E^x_y(i+1/2,j)
c   Electric field     Exz(i,j): E^x_z(i+1/2,j)
         ap=dmax1(ur+c_fr,ul+c_fl,0.0d0)   !u_R: u^W(i+1,j); u_L: u^E(i,j)
         an=dmin1(ur-c_fr,ul-c_fl,0.0d0)
         temp=ap-an
c--------
c         ExL=wl*byl-vl*bzl
c         ExR=wr*byr-vr*bzr
         EyL=ul*bzl-wl*bxl
         EyR=ur*bzr-wr*bxr
         EzL=vl*bxl-ul*byl
         EzR=vr*bxr-ur*byr
c--------
         Exx(i,j)=(                    ap*an*(BxR-BxL))/temp
         Exy(i,j)=(ap*(-EzL)-an*(-EzR)+ap*an*(ByR-ByL))/temp
         Exz(i,j)=(ap*EyL   -an*EyR   +ap*an*(BzR-BzL))/temp



c--------End of Ziegler---------------
c--------NOTE: p05 and p15 are p_star

c----------------- CASE (ZHA3)
         SELECT CASE (SCHEME1)

         CASE (0)  ! Ziegler TVD
             FFLUX(I,J,1)=(ap*ql(i,2)-an*qr(i,2)
     $                    +ap*an*(qr(i,1)-ql(i,1)))/temp
             FFLUX(I,J,2)=( ap*(rhol*ul**2+p05-bxl*bxl)
     $                     -an*(rhor*ur**2+p15-bxr*bxr)
     $                    +ap*an*(qr(i,2)-ql(i,2)))/temp
             FFLUX(I,J,3)=( ap*(rhol*ul*vl-bxl*byl)
     $                     -an*(rhor*ur*vr-bxr*byr)
     $                    +ap*an*(qr(i,3)-ql(i,3)))/temp
             FFLUX(I,J,4)=( ap*(rhol*ul*wl-bxl*bzl)
     $                     -an*(rhor*ur*wr-bxr*bzr)
     $                    +ap*an*(qr(i,4)-ql(i,4)))/temp
             FFLUX(I,J,8)=(ap*(rhol*ul*hl-bxl*(bxl*ul+byl*vl+bzl*wl))
     $                    -an*(rhor*ur*hr-bxr*(bxr*ur+byr*vr+bzr*wr))
     $                    +ap*an*(qr(i,8)-ql(i,8)))/temp

             FFLUX(I,J,5)=0.0
             FFLUX(I,J,6)=(ap*(-EzL)-an*(-EzR)
     $                    +ap*an*(qr(i,6)-ql(i,6)))/temp
             FFLUX(I,J,7)=(ap*EyL-an*EyR
     $                    +ap*an*(qr(i,7)-ql(i,7)))/temp

         CASE (3)

c------------convective terms (Fast wave speed)

c--------------------------------------------------------------

         CASE (33)

            c05=c_fl
            c15=c_fr

            c05=c_fl+al
            c15=c_fr+ar

            C_MID = 0.5*(C05 + C15)

            QML = U05/C_MID
            QMR = U15/C_MID

            if(QML.ge.0.)then
               sign_l=1.0
            else
               sign_l=-1.0
            end if

            if(QMR.ge.0.)then 
               sign_r=1.0
            else
               sign_r=-1.0
            end if

            AFA_LP = 0.5*( 1.0 + sign_l ) 
            AFA_RN = 0.5*( 1.0 - sign_r ) 

            BTA_L = -MAX(0.0, (1.0-INT(ABS(QML))) )
            BTA_R = -MAX(0.0, (1.0-INT(ABS(QMR))) )

            C_P = AFA_LP*(1.0 +BTA_L)*QML-0.25*BTA_L*(QML+1.0)**2
            C_N = AFA_RN*(1.0 +BTA_R)*QMR+0.25*BTA_R*(QMR-1.0)**2

            DLP = AFA_LP*(1. + BTA_L) - 0.5*BTA_L*(1.0+QML)
            DRN = AFA_RN*(1. + BTA_R) - 0.5*BTA_R*(1.0-QMR)

            FFLUX(I,J,1)=C_MID*(D05*C_P     + D15*C_N )

	    FFLUX(I,J,2)=C_MID*(D05*C_P*UL + D15*C_N*UR)

	    FFLUX(I,J,3)=C_MID*(D05*C_P*VL + D15*C_N*VR)

	    FFLUX(I,J,4)=C_MID*(D05*C_P*WL + D15*C_N*WR)

	    FFLUX(I,J,5)=0.0

            FFLUX(I,J,6)=C_MID*(C_P*BYL + C_N*BYR)

            FFLUX(I,J,7)=C_MID*(C_P*BZL + C_N*BZR)

            FFLUX(I,J,8)=C_MID*(C_P*RHOEL + C_N*RHOER ) 
c---------   momentum eq.(2)
            FFLUX(I,J,2) = FFLUX(I,J,2) + 
     $                  (DLP*(P05-BXL**2) + DRN*(P15-BXR**2))
c---------   energy eq. (7)
            SL_P = AFA_LP*(1. + BTA_L)*QML -0.25*BTA_L*(1.0+QML)**2
     $                                              *(2.-QML)*QML
            SR_N = AFA_RN*(1. + BTA_R)*QMR +0.25*BTA_R*(1.0-QMR)**2
     $                                              *(2.+QMR)*QMR
            FFLUX(I,J,8)=FFLUX(I,J,8)+C_MID*(C_P+C_N)*
     $                           (P05*DLP + P15*DRN) !(1) BETTER THAN (2)
c---------

            FFLUX(I,J,6) = FFLUX(I,J,6) - (DLP*BXL*VL + DRN*BXR*VR)
            FFLUX(I,J,7) = FFLUX(I,J,7) - (DLP*BXL*WL + DRN*BXR*WR)

c----------
            FFLUX(I,J,8)=FFLUX(I,J,8)-
     $                (DLP*BXL*(UL*BXL+VL*BYL+WL*BZL)
     $               + DRN*BXR*(UR*BXR+VR*BYR+WR*BZR))
            FFLUX(I,J,3) = FFLUX(I,J,3) - (DLP*BXL*BYL + DRN*BXR*BYR)
            FFLUX(I,J,4) = FFLUX(I,J,4) - (DLP*BXL*BZL + DRN*BXR*BZR)

c--------------------------------------------------------------

         CASE (333)
            c05=al
            c15=ar


            C_MID = 0.5*(C05 + C15)

            QML = U05/C_MID
            QMR = U15/C_MID

            if(QML.ge.0.)then
               sign_l=1.0

            else
               sign_l=-1.0
            end if

            if(QMR.ge.0.)then 

               sign_r=1.0
            else
               sign_r=-1.0
            end if

            AFA_LP = 0.5*( 1.0 + sign_l ) 
            AFA_RN = 0.5*( 1.0 - sign_r ) 

            BTA_L = -MAX(0.0, (1.0-INT(ABS(QML))) )
            BTA_R = -MAX(0.0, (1.0-INT(ABS(QMR))) )

            C_P = AFA_LP*(1.0 +BTA_L)*QML-0.25*BTA_L*(QML+1.0)**2
            C_N = AFA_RN*(1.0 +BTA_R)*QMR+0.25*BTA_R*(QMR-1.0)**2

            DLP = AFA_LP*(1. + BTA_L) - 0.5*BTA_L*(1.0+QML)
            DRN = AFA_RN*(1. + BTA_R) - 0.5*BTA_R*(1.0-QMR)

            FFLUX(I,J,1)=C_MID*(D05*C_P     + D15*C_N )

	    FFLUX(I,J,2)=C_MID*(D05*C_P*UL + D15*C_N*UR)

	    FFLUX(I,J,3)=C_MID*(D05*C_P*VL + D15*C_N*VR)

	    FFLUX(I,J,4)=C_MID*(D05*C_P*WL + D15*C_N*WR)

	    FFLUX(I,J,5)=0.0

            FFLUX(I,J,6)=C_MID*(C_P*BYL + C_N*BYR)

            FFLUX(I,J,7)=C_MID*(C_P*BZL + C_N*BZR)

            FFLUX(I,J,8)=C_MID*(C_P*RHOEL + C_N*RHOER ) 
c---------   momentum eq.(2)
            FFLUX(I,J,2) = FFLUX(I,J,2) + 
     $                  (DLP*(P05-BXL**2) + DRN*(P15-BXR**2))
c---------   energy eq. (8)
            FFLUX(I,J,8)=FFLUX(I,J,8)+C_MID*(C_P+C_N)*
     $                           (P05*DLP + P15*DRN) !(1) BETTER THAN (2)
c---------

            FFLUX(I,J,6) = FFLUX(I,J,6) - (DLP*BXL*VL + DRN*BXR*VR)
            FFLUX(I,J,7) = FFLUX(I,J,7) - (DLP*BXL*WL + DRN*BXR*WR)


c--------------using different Mach number

c            FFLUX(I,J,3) = FFLUX(I,J,3) - (DLP*BXL*BYL + DRN*BXR*BYR)
c            FFLUX(I,J,4) = FFLUX(I,J,4) - (DLP*BXL*BZL + DRN*BXR*BZR)

c            FFLUX(I,J,8)=FFLUX(I,J,8)-
c     $                (DLP*BXL*(UL*BXL+VL*BYL+WL*BZL)
c     $               + DRN*BXR*(UR*BXR+VR*BYR+WR*BZR))

c----------------------
            c05=c_fl
            c15=c_fr


            C_MID = 0.5*(C05 + C15)

            QML = U05/C_MID
            QMR = U15/C_MID

            if(QML.ge.0.)then
               sign_l=1.0
            else
               sign_l=-1.0
            end if

            if(QMR.ge.0.)then 
               sign_r=1.0
            else
               sign_r=-1.0
            end if
            AFA_LP = 0.5*( 1.0 + sign_l ) 
            AFA_RN = 0.5*( 1.0 - sign_r ) 

            BTA_L = -MAX(0.0, (1.0-INT(ABS(QML))) )
            BTA_R = -MAX(0.0, (1.0-INT(ABS(QMR))) )

            AFA_LP = 0.5*( 1.0 + sign_l ) 
            AFA_RN = 0.5*( 1.0 - sign_r ) 

            BTA_L = -MAX(0.0, (1.0-INT(ABS(QML))) )
            BTA_R = -MAX(0.0, (1.0-INT(ABS(QMR))) )

            C_P = AFA_LP*(1.0 +BTA_L)*QML-0.25*BTA_L*(QML+1.0)**2
            C_N = AFA_RN*(1.0 +BTA_R)*QMR+0.25*BTA_R*(QMR-1.0)**2

            DLP = AFA_LP*(1. + BTA_L) - 0.5*BTA_L*(1.0+QML)
            DRN = AFA_RN*(1. + BTA_R) - 0.5*BTA_R*(1.0-QMR)

            FFLUX(I,J,3) = FFLUX(I,J,3) - (DLP*BXL*BYL + DRN*BXR*BYR)
            FFLUX(I,J,4) = FFLUX(I,J,4) - (DLP*BXL*BZL + DRN*BXR*BZR)

            FFLUX(I,J,8)=FFLUX(I,J,8)-
     $                (DLP*BXL*(UL*BXL+VL*BYL+WL*BZL)
     $               + DRN*BXR*(UR*BXR+VR*BYR+WR*BZR))

c----------------------------------------------Optimized

        END SELECT



         Exx(i,j)=0.0
         Exy(i,j)=FFLUX(I,J,6)
         Exz(i,j)=FFLUX(I,J,7)



       end do


      end do !        j=0,NY+1


c_______________________________________________________________
c_______________________________________________________________




c---------------Y-DIRECTION-------------------------------------
      IF (CT.eq.1) then
c      DO j=0,NX+1
      DO j=-1,NX+1
         k=5  !Bx
         call reconstruct(scheme,NY,BxA(1,j,:),B1L(j,:),B1R(j,:))
         k=7  !Bz
         call reconstruct(scheme,NY,BzA(1,j,:),B2L(j,:),B2R(j,:))
      END DO
      END IF ! IF (CT.eq.1) then

c      do j=1,NX
      do j=-1,NX+1

      IF (CT.eq.1) then
      do k=1,4
         call reconstruct(scheme,NY,q(1,j,:,k),ql(:,k),qr(:,k))
      end do
         k=6  !By
c         DO I=0,NY+1
         DO I=-1,NY+1
            QL(i,6)=ByA(1,j,i)                  ! B^E(i,j,k)
            QR(i,6)=ByA(1,j,i)                  ! B^W(i,j,k)
         END DO
         !k=5,7
c         DO I=0,NX+1
         DO I=-1,NX+1
            QL(i,5)=0.5*(B1L(j,i)+B1L(j-1,i))
            QR(i,5)=0.5*(B1R(j,i)+B1R(j-1,i))
            QL(i,7)=0.5*(B2L(j,i)+B2L(j-1,i))
            QR(i,7)=0.5*(B2R(j,i)+B2R(j-1,i))

c            QL(i,5)=0.5*(B1L(j,i)+B1L(j+1,i))
c            QR(i,5)=0.5*(B1R(j,i)+B1R(j+1,i))
c            QL(i,7)=0.5*(B2L(j,i)+B2L(j+1,i))
c            QR(i,7)=0.5*(B2R(j,i)+B2R(j+1,i))
         END DO

c---------------new CTM
        do k=5,7
           call reconstruct(scheme,NX,q(1,j,:,k),ql(:,k),qr(:,k))
        end do
c---------------new CTM

         k=8
         call reconstruct(scheme,NY,q(1,j,:,k),ql(:,k),qr(:,k))
      ELSE ! IF (CT.eq.1) then

c---------reconstruction--------------
      do k=1,8
         call reconstruct(scheme,NY,q(1,j,:,k),ql(:,k),qr(:,k))
      end do
      END IF ! IF (CT.eq.1) then

c--------------calculate eigenvectors---------------
c      do i=0,NY+1
      do i=-1,NY+1

         rhol=ql(i,1)
         ul=ql(i,2)/rhol
         vl=ql(i,3)/rhol
         wl=ql(i,4)/rhol
         bxl=ql(i,5)
         byl=ql(i,6)
         bzl=ql(i,7)
         rhoel=ql(i,8)
         B2=Bxl**2+Byl**2+Bzl**2
         v2=ul**2+vl**2+wl**2
         pl=(gamma-1.)*(rhoel-0.5*(rhol*v2+B2))
         hl=(rhoel+pl+0.5*B2)/rhol
         P05=pl+B2/2.0

         b_2=B2/rhol
         al=sqrt(gamma*pl/rhol)         
         a2=al**2
         c_al=abs(Byl/sqrt(rhol))
         temp1=(a2+b_2)**2-4.0*a2*Bxl**2/rhol
c-----------------
           if (pl.le.0. .or.rhol.le.0.0.or.temp1.le.0.0) then
c              write(*,*)'reconstruct_yL',j,i,pl,rhol,temp1
c              write(*,*)temp1



         rhol=q(1,j,i,1)
         ul=q(1,j,i,2)/rhol
         vl=q(1,j,i,3)/rhol
         wl=q(1,j,i,4)/rhol
         bxl=0.5*(q(1,j,i,5)+q(1,j-1,i,5))
         byl=q(1,j,i,6)
         bzl=0.5*(q(1,j,i,7)+q(1,j-1,i,7))


c         bxl=0.5*(q(1,j,i,5)+q(1,j+1,i,5))
c         bzl=0.5*(q(1,j,i,7)+q(1,j+1,i,7))

         rhoel=q(1,j,i,8)
         B2=Bxl**2+Byl**2+Bzl**2
         v2=ul**2+vl**2+wl**2
         pl=(gamma-1.)*(rhoel-0.5*(rhol*v2+B2))
         hl=(rhoel+pl+0.5*B2)/rhol
         P05=pl+B2/2.0

         b_2=B2/rhol
         al=sqrt(gamma*pl/rhol)         
         a2=al**2
         c_al=abs(Byl/sqrt(rhol))
         temp1=(a2+b_2)**2-4.0*a2*Byl**2/rhol
c         write(*,*)'----------YL',temp1,a2,b_2,Byl,rhol,pl
c              read(*,*)k


           end if
c-----------------
c         al=sqrt(gamma*pl/rhol)         
c         a2=al**2
c         c_al=abs(Byl/sqrt(rhol))
c         temp1=(a2+b_2)**2-4.0*a2*Bxl**2/rhol
c         temp2=0.5*(a2+b_2-sqrt((a2+b_2)**2-4.0*a2*Bxl**2/rhol))
c         if (temp1.le.0.) then
c            write(*,*)'temp1-yl',j,i,temp1
c         end if
c         if (temp2.le.0.) then
c            write(*,*)'temp2-yl',j,i,temp2
c         end if
         c_fl=sqrt(0.5*(a2+b_2+sqrt((a2+b_2)**2-4.0*a2*Byl**2/rhol)))
         c_sl=sqrt(0.5*(a2+b_2-sqrt((a2+b_2)**2-4.0*a2*Byl**2/rhol)))

         rhor=qr(i,1)
         ur=qr(i,2)/rhor
         vr=qr(i,3)/rhor
         wr=qr(i,4)/rhor
         bxr=qr(i,5)
         byr=qr(i,6)
         bzr=qr(i,7)
         rhoer=qr(i,8)
         B2=Bxr**2+Byr**2+Bzr**2
         v2=ur**2+vr**2+wr**2
         pr=(gamma-1.)*(rhoer-0.5*(rhor*v2+B2))
         hr=(rhoer+pr+0.5*B2)/rhor
         P15=pr+B2/2.0

         b_2=B2/rhor
         ar=sqrt(gamma*pr/rhor)         
         a2=ar**2
         c_ar=abs(Byr/sqrt(rhor))
         temp1=(a2+b_2)**2-4.0*a2*Byr**2/rhor
c-----------------
           if (pr.le.0. .or.rhor.le.0.0.or.temp1.le.0.0) then
c              write(*,*)'reconstruct_yR',j,i,pr,rhor,temp1
c              write(*,*)temp1
c              read(*,*)k


         rhor=q(1,j,i+1,1)
         ur=q(1,j,i+1,2)/rhor
         vr=q(1,j,i+1,3)/rhor
         wr=q(1,j,i+1,4)/rhor
         bxr=0.5*(q(1,j,i+1,5)+q(1,j-1,i+1,5))
         byr=q(1,j,i,6)
         bzr=0.5*(q(1,j,i+1,7)+q(1,j-1,i+1,7))

c         bxr=0.5*(q(1,j,i+1,5)+q(1,j+1,i+1,5))
c         bzr=0.5*(q(1,j,i+1,7)+q(1,j+1,i+1,7))


         rhoer=q(1,j,i+1,8)
         B2=Bxr**2+Byr**2+Bzr**2
         v2=ur**2+vr**2+wr**2
         pr=(gamma-1.)*(rhoer-0.5*(rhor*v2+B2))
         hr=(rhoer+pr+0.5*B2)/rhor
         P15=pr+B2/2.0

         b_2=B2/rhor
         ar=sqrt(gamma*pr/rhor)         
         a2=ar**2
         c_ar=abs(Byr/sqrt(rhor))

         temp1=(a2+b_2)**2-4.0*a2*Byr**2/rhor
c         write(*,*)'----------YR',temp1,a2,b_2,Byr,rhor,pr
c              read(*,*)k
c         write(*,*)'----------2',rhoer,0.5*(rhor*v2+B2)
         




           end if
c-----------------
c         ar=sqrt(gamma*pr/rhor)         
c         a2=ar**2
c         c_ar=abs(Byr/sqrt(rhor))
c         temp1=(a2+b_2)**2-4.0*a2*Bxr**2/rhor
c         temp2=0.5*(a2+b_2-sqrt((a2+b_2)**2-4.0*a2*Bxr**2/rhor))
c         if (temp1.le.0.) then
c            write(*,*)'temp1-yr',j,i,temp1
c         end if
c         if (temp2.le.0.) then
c            write(*,*)'temp2-yr',j,i,temp2
c         end if
         c_fr=sqrt(0.5*(a2+b_2+sqrt((a2+b_2)**2-4.0*a2*Byr**2/rhor)))
         c_sr=sqrt(0.5*(a2+b_2-sqrt((a2+b_2)**2-4.0*a2*Byr**2/rhor)))



         d05=rhol
         d15=rhor
         V05=vl
         V15=vr

c--------Ziegler----------------------
c   ap: b^(+)(i,j+1/2,k)

c   an: b^(-)(i,j+1/2,k)
c   Electric field     Eyx(i,j): E^y_x(i+1/2,j)
c   Electric field     Eyy(i,j): E^y_y(i+1/2,j)
c   Electric field     Eyz(i,j): E^y_z(i+1/2,j)
         ap=dmax1(vr+c_fr,vl+c_fl,0.0d0)   !v_R: v^S(i+1,j); v_L: v^N(i,j)
         an=dmin1(vr-c_fr,vl-c_fl,0.0d0)
         temp=ap-an
c--------
         ExL=wl*byl-vl*bzl
         ExR=wr*byr-vr*bzr
c         EyL=ul*bzl-wl*bxl
c         EyR=ur*bzr-wr*bxr
         EzL=vl*bxl-ul*byl
         EzR=vr*bxr-ur*byr
c--------
         Eyx(j,i)=(ap*EzL  -  an*EzR+  ap*an*(BxR-BxL))/temp
         Eyy(j,i)=(                    ap*an*(ByR-ByL))/temp
         Eyz(j,i)=(ap*(-ExL)-an*(-ExR)+ap*an*(BzR-BzL))/temp



c--------End of Ziegler---------------

c***********Z-Direction**************************************
cc--------Ziegler----------------------
cc   ap: c^(+)(i,j,k+1/2)
cc   an: c^(-)(i,j,k+1/2)
cc   Electric field     Ezx(i,j): E^z_x(i,j,k+1/2)
cc   Electric field     Ezy(i,j): E^z_y(i,j,k+1/2)
cc   Electric field     Ezz(i,j): E^z_z(i,j,k+1/2)
c         ap=dmax1(wr+c_fr,wl+c_fl,0.0)   !w_R: v^B(i+1,j); v_L: v^T(i,j)
c         an=dmin1(wr-c_fr,wl-c_fl,0.0)
c         temp=ap-an
cc--------
c         ExL=wl*byl-vl*bzl
c         ExR=wr*byr-vr*bzr
c         EyL=ul*bzl-wl*bxl
c         EyR=ur*bzr-wr*bxr
c         EzL=vl*bxl-ul*byl
c         EzR=vr*bxr-ur*byr
cc--------
c         Ezx(i,j)=(ap*(-EyL)-an*(-EyR)+ap*an*(BxR-BxL))/temp
c         Ezy(i,j)=(ap*ExL   -an*ExR   +ap*an*(ByR-ByL))/temp
c         Ezz(i,j)=(                    ap*an*(BzR-BzL))/temp
c--------End of Ziegler---------------
c***********************************************************



c---------Fluxes for Magnetic field at surface[i-1=====>i-1/2]
c         d(Bx)/dT_(i-1/2)=-(Ezf(i-1,j+1,k)-Ezf(i-1,j,k))/Dy
c                          +(Eyf(i-1,j,k+1)-Eyf(i-1,j,k))/Dz
c         d(By)/dT_(j-1/2)= (Ezf(i,j-1,k)-Ezf(i-1,j-1,k))/Dx
c                          -(Exf(i,j-1,k)-Exf(i-1,j,k-1))/Dz
c         d(Bz)/dT_(k-1/2)=-(Eyf(i,j,k-1)-Eyf(i-1,j,k-1))/Dx
c                          +(Exf(i,j,k-1)-Exf(i,j-1,k-1))/Dy
c
c

c--------NOTE: p05 and p15 are p_star

c----------------- CASE (ZHA3)
         SELECT CASE (SCHEME1)
         CASE (0)  ! Ziegler TVD
             GFLUX(J,I,1)=(ap*ql(i,3)-an*qr(i,3)
     $                    +ap*an*(qr(i,1)-ql(i,1)))/temp
             GFLUX(J,I,2)=( ap*(rhol*ul*vl-byl*bxl)
     $                     -an*(rhor*ur*vr-byr*bxr)
     $                    +ap*an*(qr(i,2)-ql(i,2)))/temp
             GFLUX(J,I,3)=( ap*(rhol*vl**2+p05-byl*byl)
     $                     -an*(rhor*vr**2+p15-byr*byr)
     $                    +ap*an*(qr(i,3)-ql(i,3)))/temp
             GFLUX(J,I,4)=( ap*(rhol*vl*wl-byl*bzl)
     $                     -an*(rhor*vr*wr-byr*bzr)
     $                    +ap*an*(qr(i,4)-ql(i,4)))/temp
             GFLUX(J,I,8)=(ap*(rhol*vl*hl-byl*(bxl*ul+byl*vl+bzl*wl))
     $                    -an*(rhor*vr*hr-byr*(bxr*ur+byr*vr+bzr*wr))
     $                    +ap*an*(qr(i,8)-ql(i,8)))/temp

             GFLUX(J,I,5)=(ap*EzL-an*EzR
     $                    +ap*an*(qr(i,5)-ql(i,5)))/temp
             GFLUX(J,I,6)=0.0
             GFLUX(J,I,7)=(ap*(-ExL)-an*(-ExR)
     $                    +ap*an*(qr(i,7)-ql(i,7)))/temp

            CASE (3)

c------------convective terms (Fast wave speed)

c--------------------------------------------------------------

         CASE (33)

            c05=c_fl
            c15=c_fr

c            c05=c_fl+al
c            c15=c_fr+ar

            C_MID = 0.5*(C05 + C15)

            QML = V05/C_MID
            QMR = V15/C_MID

            if(QML.ge.0.)then
               sign_l=1.0
            else
               sign_l=-1.0
            end if

            if(QMR.ge.0.)then 
               sign_r=1.0
            else
               sign_r=-1.0
            end if

            AFA_LP = 0.5*( 1.0 + sign_l ) 
            AFA_RN = 0.5*( 1.0 - sign_r ) 

            BTA_L = -MAX(0.0, (1.0-INT(ABS(QML))) )
            BTA_R = -MAX(0.0, (1.0-INT(ABS(QMR))) )

            C_P = AFA_LP*(1.0 +BTA_L)*QML-0.25*BTA_L*(QML+1.0)**2
            C_N = AFA_RN*(1.0 +BTA_R)*QMR+0.25*BTA_R*(QMR-1.0)**2

            DLP = AFA_LP*(1. + BTA_L) - 0.5*BTA_L*(1.0+QML)
            DRN = AFA_RN*(1. + BTA_R) - 0.5*BTA_R*(1.0-QMR)

            GFLUX(J,I,1)=C_MID*(D05*C_P     + D15*C_N )

	    GFLUX(J,I,2)=C_MID*(D05*C_P*UL + D15*C_N*UR)

	    GFLUX(J,I,3)=C_MID*(D05*C_P*VL + D15*C_N*VR)

	    GFLUX(J,I,4)=C_MID*(D05*C_P*WL + D15*C_N*WR)

            GFLUX(J,I,5)=C_MID*(C_P*BXL + C_N*BXR)

	    GFLUX(J,I,6)=0.0

            GFLUX(J,I,7)=C_MID*(C_P*BZL + C_N*BZR)

            GFLUX(J,I,8)=C_MID*(C_P*RHOEL + C_N*RHOER ) 
c---------   momentum eq.(3)
            GFLUX(J,I,3) = GFLUX(J,I,3) + 
     $                  (DLP*(P05-BYL**2) + DRN*(P15-BYR**2))
c---------   energy eq. (7)
            SL_P = AFA_LP*(1. + BTA_L)*QML -0.25*BTA_L*(1.0+QML)**2
     $                                              *(2.-QML)*QML
            SR_N = AFA_RN*(1. + BTA_R)*QMR +0.25*BTA_R*(1.0-QMR)**2
     $                                              *(2.+QMR)*QMR
            GFLUX(J,I,8)=GFLUX(J,I,8)+C_MID*(C_P+C_N)*
     $                           (P05*DLP + P15*DRN) !(1) BETTER THAN (2)
c---------

            GFLUX(J,I,5) = GFLUX(J,I,5) - (DLP*BYL*UL + DRN*BYR*UR)
            GFLUX(J,I,7) = GFLUX(J,I,7) - (DLP*BYL*WL + DRN*BYR*WR)

c----------
            GFLUX(J,I,8)=GFLUX(J,I,8)-
     $                (DLP*BYL*(UL*BXL+VL*BYL+WL*BZL)
     $               + DRN*BYR*(UR*BXR+VR*BYR+WR*BZR))
            GFLUX(J,I,2) = GFLUX(J,I,2) - (DLP*BXL*BYL + DRN*BXR*BYR)
            GFLUX(J,I,4) = GFLUX(J,I,4) - (DLP*BYL*BZL + DRN*BYR*BZR)

c--------------------------------------------------------------

         CASE (333)
            c05=al
            c15=ar


            C_MID = 0.5*(C05 + C15)

            QML = V05/C_MID
            QMR = V15/C_MID

            if(QML.ge.0.)then
               sign_l=1.0
            else
               sign_l=-1.0
            end if

            if(QMR.ge.0.)then 
               sign_r=1.0
            else
               sign_r=-1.0
            end if

            AFA_LP = 0.5*( 1.0 + sign_l ) 
            AFA_RN = 0.5*( 1.0 - sign_r ) 

            BTA_L = -MAX(0.0, (1.0-INT(ABS(QML))) )
            BTA_R = -MAX(0.0, (1.0-INT(ABS(QMR))) )

            C_P = AFA_LP*(1.0 +BTA_L)*QML-0.25*BTA_L*(QML+1.0)**2
            C_N = AFA_RN*(1.0 +BTA_R)*QMR+0.25*BTA_R*(QMR-1.0)**2

            DLP = AFA_LP*(1. + BTA_L) - 0.5*BTA_L*(1.0+QML)
            DRN = AFA_RN*(1. + BTA_R) - 0.5*BTA_R*(1.0-QMR)

            GFLUX(J,I,1)=C_MID*(D05*C_P     + D15*C_N )

	    GFLUX(J,I,2)=C_MID*(D05*C_P*UL + D15*C_N*UR)

	    GFLUX(J,I,3)=C_MID*(D05*C_P*VL + D15*C_N*VR)

	    GFLUX(J,I,4)=C_MID*(D05*C_P*WL + D15*C_N*WR)

            GFLUX(J,I,5)=C_MID*(C_P*BXL + C_N*BXR)

	    GFLUX(J,I,6)=0.0

            GFLUX(J,I,7)=C_MID*(C_P*BZL + C_N*BZR)

            GFLUX(J,I,8)=C_MID*(C_P*RHOEL + C_N*RHOER ) 
c---------   momentum eq.(3)
            GFLUX(J,I,3) = GFLUX(J,I,3) + 
     $                  (DLP*(P05-BYL**2) + DRN*(P15-BYR**2))
c---------   energy eq. (8)
            GFLUX(J,I,8)=GFLUX(J,I,8)+C_MID*(C_P+C_N)*
     $                           (P05*DLP + P15*DRN) !(1) BETTER THAN (2)
c---------

            GFLUX(J,I,5) = GFLUX(J,I,5) - (DLP*BYL*UL + DRN*BYL*UR)
            GFLUX(J,I,7) = GFLUX(J,I,7) - (DLP*BYL*WL + DRN*BYR*WR)



c----------------------
            c05=c_fl
            c15=c_fr


            C_MID = 0.5*(C05 + C15)

            QML = V05/C_MID
            QMR = V15/C_MID

            if(QML.ge.0.)then
               sign_l=1.0
            else
               sign_l=-1.0
            end if

            if(QMR.ge.0.)then 
               sign_r=1.0
            else
               sign_r=-1.0
            end if
            AFA_LP = 0.5*( 1.0 + sign_l ) 
            AFA_RN = 0.5*( 1.0 - sign_r ) 

            BTA_L = -MAX(0.0, (1.0-INT(ABS(QML))) )
            BTA_R = -MAX(0.0, (1.0-INT(ABS(QMR))) )

            AFA_LP = 0.5*( 1.0 + sign_l ) 
            AFA_RN = 0.5*( 1.0 - sign_r ) 

            BTA_L = -MAX(0.0, (1.0-INT(ABS(QML))) )
            BTA_R = -MAX(0.0, (1.0-INT(ABS(QMR))) )

            C_P = AFA_LP*(1.0 +BTA_L)*QML-0.25*BTA_L*(QML+1.0)**2
            C_N = AFA_RN*(1.0 +BTA_R)*QMR+0.25*BTA_R*(QMR-1.0)**2

            DLP = AFA_LP*(1. + BTA_L) - 0.5*BTA_L*(1.0+QML)
            DRN = AFA_RN*(1. + BTA_R) - 0.5*BTA_R*(1.0-QMR)

            GFLUX(J,I,2) = GFLUX(J,I,2) - (DLP*BXL*BYL + DRN*BXR*BYR)
            GFLUX(J,I,4) = GFLUX(J,I,4) - (DLP*BYL*BZL + DRN*BYR*BZR)

            GFLUX(J,I,8)=GFLUX(J,I,8)-
     $                (DLP*BYL*(UL*BXL+VL*BYL+WL*BZL)
     $               + DRN*BYR*(UR*BXR+VR*BYR+WR*BZR))

c----------------------------------------------Optimized

        END SELECT
         Eyx(j,i)=GFLUX(J,I,5)
         Eyy(j,i)=0.0
         Eyz(j,i)=GFLUX(J,I,7)

      end do
      end do !        j=0,NX+1

c--------------Ziegler for 3D--------------------------
C      FFLUX(I,J,5)= 0.0
C      FFLUX(I,J,6)=-Ezf
C      FFLUX(I,J,7)= Eyf
c------
C      GFLUX(I,J,5)= Ezf
C      GFLUX(I,J,6)= 0.0
C      GFLUX(I,J,7)=-Exf
c------
C      HFLUX(I,J,5)=-Eyf
C      HFLUX(I,J,6)= Exf
C      HFLUX(I,J,7)= 0.0

      IF (CT.eq.1) then

c      DO I=0,NX+1
c         DO J=0,NY+1
c            Exf(I,J)=-Eyz(I,J)
c            Eyf(I,J)= Exz(I,J)
c            Ezf(I,J)=(-Exy(I,J)-Exy(I,J+1)+Eyx(I,J)+Eyx(I+1,J))/4.0

c            Exf(I,J)=-Eyz(I,J)
c            Eyf(I,J)= Exz(I,J)
c            Ezf(I,J)=(-Exy(I,J)-Exy(I,J+1)+Eyx(I,J)+Eyx(I+1,J))/4.0
c         END DO
c      END DO

      DO I=0,NX+1
         DO J=0,NY+1
            Exf(I,J)=-Eyz(I,J)
            Eyf(I,J)= Exz(I,J)
c            Ezf(I,J)=(-Exy(I,J)-Exy(I,J-1)+Eyx(I,J)+Eyx(I-1,J))/4.0
            Ezf(I,J)=(-Exy(I,J)-Exy(I,J+1)+Eyx(I,J)+Eyx(I+1,J))/4.0
c------------
c            FFLUX(I,J,5)= 0.0
c            FFLUX(I,J,6)=-Ezf(I,J)
c            FFLUX(I,J,7)= Eyf(I,J)
c------
c            GFLUX(I,J,5)= Ezf(I,J)
c            GFLUX(I,J,6)= 0.0
c            GFLUX(I,J,7)=-Exf(I,J)

c------------

         END DO
      END DO
      END IF !IF (CT.eq.1) then
c_______________________________________________________________
c_______________________________________________________________


c--------source term
      do i=1,NX
         do j=1,NY
            temp= (q(1,i+1,j,5)-q(1,i-1,j,5))/(2.*dx)
     $           +(q(1,i,j+1,6)-q(1,i,j-1,6))/(2.*dy) 
cc     $           +(q(1,i,j,k+1,7)-q(1,i,j,k-1,7))/(2.*dz) 
              
c              temp=source(i,j,1)

              source(i,j,1)=0.0
              source(i,j,2)=temp*q(1,i,j,5)
              source(i,j,3)=temp*q(1,i,j,6)
              source(i,j,4)=temp*q(1,i,j,7)
              ua=q(1,i,j,2)/q(1,i,j,1)
              va=q(1,i,j,3)/q(1,i,j,1)
              wa=q(1,i,j,4)/q(1,i,j,1)
              source(i,j,5)=temp*ua
              source(i,j,6)=temp*va
              source(i,j,7)=temp*wa
              source(i,j,8)=temp*
     $                     (ua*q(1,i,j,5)+va*q(1,i,j,6)+wa*q(1,i,j,7))
          end do
       end do


       source = 0.0





















c--------Runge-Kutta method for time marching--------------
c    (1)  the first-order Runge-Kutta scheme
      IF (NRK.EQ.1) THEN
         do i=1,NX
            do j=1,NY
              do i1=1,8
               q(1,i,j,i1)=q(1,i,j,i1)
     $                        -dt*(fflux(i,j,i1)-fflux(i-1,j,i1))/dx
     $                        -dt*(gflux(i,j,i1)-gflux(i,j-1,i1))/dy
     $                        -dt*source(i,j,i1)
              end do
            end do
         end do
      END IF
c    (2)  the second-order Runge-Kutta scheme
      IF (NRK.EQ.2) THEN
          if (IRK.eq.1) then
             q(5,:,:,:)=q(1,:,:,:)
             do i=1,NX
               do j=1,NY
                 do i1=1,8
                   q(1,i,j,i1)=q(1,i,j,i1)
     $                  -dt*(fflux(i,j,i1)-fflux(i-1,j,i1))/dx
     $                  -dt*(gflux(i,j,i1)-gflux(i,j-1,i1))/dy
     $                  -dt*source(i,j,i1)
                 end do
               end do
             end do
          end if
          if (IRK.eq.2) then
             do i=1,NX
               do j=1,NY
                 do i1=1,8
                   q(1,i,j,i1)=0.5*(q(5,i,j,i1)+q(1,i,j,i1)
     $                      -dt*(fflux(i,j,i1)-fflux(i-1,j,i1))/dx
     $                      -dt*(gflux(i,j,i1)-gflux(i,j-1,i1))/dy
     $                      -dt*source(i,j,i1))
                 end do
               end do
             end do
          end if
      END iF
c    (3)  the third-order Runge-Kutta scheme
      IF (NRK.EQ.3) THEN
          if (IRK.eq.1) then
             q(5,:,:,:)=q(1,:,:,:)
             do i=1,NX
               do j=1,NY
                 do i1=1,8
                  q(1,i,j,i1)=q(1,i,j,i1)
     $                  -dt*(fflux(i,j,i1)-fflux(i-1,j,i1))/dx
     $                  -dt*(gflux(i,j,i1)-gflux(i,j-1,i1))/dy
     $                  -dt*source(i,j,i1)
                  if (q(1,i,j,5)**2+q(1,i,j,6)**2
     $                        +q(1,i,j,7)**2.lt.0.0) then
                      write(*,*)'RK2',i,j,q(1,i,j,5),q(1,i,j,6),
     $                           q(1,i,j,7)                           
                   end if
                 end do
               end do
             end do
             
c--------------new CTM
             if (CT.EQ.1) then
             BxA(5,:,:)=BxA(1,:,:)
             ByA(5,:,:)=ByA(1,:,:)
             BzA(5,:,:)=BzA(1,:,:)
             do i=1,NX
               do j=1,NY
c                  BxA(1,i,j)=BxA(1,i,j)
c     $                  -dt*(fflux(i,j,5)-fflux(i-1,j,5))/dx
c     $                  -dt*(gflux(i,j,5)-gflux(i,j-1,5))/dy
c     $                  -dt*source(i,j,5)
c                  ByA(1,i,j)=ByA(1,i,j)
c     $                  -dt*(fflux(i,j,6)-fflux(i-1,j,6))/dx
c     $                  -dt*(gflux(i,j,6)-gflux(i,j-1,6))/dy
c     $                  -dt*source(i,j,6)
c                  BzA(1,i,j)=BzA(1,i,j)
c     $                  -dt*(fflux(i,j,7)-fflux(i-1,j,7))/dx
c     $                  -dt*(gflux(i,j,7)-gflux(i,j-1,7))/dy
c     $                  -dt*source(i,j,7)     
     
                  BxA(1,i,j)=BxA(1,i,j)
     $                  -dt*(Ezf(i,j)-Ezf(i,j-1))/dy
     $                  -dt*source(i,j,5)
                  ByA(1,i,j)=ByA(1,i,j)
     $                  +dt*(Ezf(i,j)-Ezf(i-1,j))/dx
     $                  -dt*source(i,j,6)
                  BzA(1,i,j)=BzA(1,i,j)
     $                  -dt*(Eyf(i,j)-Eyf(i-1,j))/dx
     $                  +dt*(Exf(i,j)-Exf(i,j-1))/dy
     $                  -dt*source(i,j,7)     
     
               end do
             end do
             
             end if
c--------------new CTM             
          end if

          IF (IRK.EQ.2) THEN
             q(4,:,:,:)=q(1,:,:,:)
             do i=1,NX
               do j=1,NY
                 do i1=1,8
                   q(1,i,j,i1)=0.25*(3.*q(5,i,j,i1)+q(4,i,j,i1)
     $                    -dt*(fflux(i,j,i1)-fflux(i-1,j,i1))/dx
     $                    -dt*(gflux(i,j,i1)-gflux(i,j-1,i1))/dy
     $                    -dt*source(i,j,i1))
                  if (q(1,i,j,5)**2+q(1,i,j,6)**2
     $                        +q(1,i,j,7)**2.lt.0.0) then
                      write(*,*)'RK2',i,j,q(1,i,j,5),q(1,i,j,6),
     $                           q(1,i,j,7)
                   end if
                 end do
               end do
             end do
c--------------new CTM
             if (CT.EQ.1) then
             BxA(4,:,:)=BxA(1,:,:)
             ByA(4,:,:)=ByA(1,:,:)
             BzA(4,:,:)=BzA(1,:,:)
             do i=1,NX
               do j=1,NY
C                   BxA(1,i,j)=0.25*(3.*BxA(5,i,j)+BxA(4,i,j)
c     $                    -dt*(fflux(i,j,5)-fflux(i-1,j,5))/dx
c     $                    -dt*(gflux(i,j,5)-gflux(i,j-1,5))/dy
c     $                    -dt*source(i,j,5))
c                   ByA(1,i,j)=0.25*(3.*ByA(5,i,j)+ByA(4,i,j)
c     $                    -dt*(fflux(i,j,6)-fflux(i-1,j,6))/dx
c     $                    -dt*(gflux(i,j,6)-gflux(i,j-1,6))/dy
c     $                    -dt*source(i,j,6))
c                   BzA(1,i,j)=0.25*(3.*BzA(5,i,j)+BzA(4,i,j)
c     $                    -dt*(fflux(i,j,7)-fflux(i-1,j,7))/dx
c     $                    -dt*(gflux(i,j,7)-gflux(i,j-1,7))/dy
c     $                    -dt*source(i,j,7))

                   BxA(1,i,j)=0.25*(3.*BxA(5,i,j)+BxA(4,i,j)
     $                    -dt*(Ezf(i,j)-Ezf(i,j-1))/dy
     $                    -dt*source(i,j,5))
                   ByA(1,i,j)=0.25*(3.*ByA(5,i,j)+ByA(4,i,j)
     $                    +dt*(Ezf(i,j)-Ezf(i-1,j))/dx
     $                    -dt*source(i,j,6))
                   BzA(1,i,j)=0.25*(3.*BzA(5,i,j)+BzA(4,i,j)
     $                    -dt*(Eyf(i,j)-Eyf(i-1,j))/dx
     $                    +dt*(Exf(i,j)-Exf(i,j-1))/dy
     $                    -dt*source(i,j,7))
     
               end do
             end do
             
             end if
c--------------new CTM                    
          end if
          if (IRK.eq.3) then
             do i=1,NX
               do j=1,NY
                 do i1=1,8
                   q(1,i,j,i1)=(q(5,i,j,i1)+2.*q(1,i,j,i1)
     $               -2.0*dt*(fflux(i,j,i1)-fflux(i-1,j,i1))/dx
     $               -2.0*dt*(gflux(i,j,i1)-gflux(i,j-1,i1))/dy
     $               -2.0*dt*source(i,j,i1))/3.0

                   if (q(1,i,j,5)**2+q(1,i,j,6)**2
     $                        +q(1,i,j,7)**2.lt.0.0) then
                      write(*,*)'RK3',i,j,q(1,i,j,5),q(1,i,j,6),
     $                           q(1,i,j,7)
                   end if
                 end do
               end do
             end do
c--------------new CTM
             if (CT.EQ.1) then
             do i=1,NX
               do j=1,NY
c                   BxA(1,i,j)=(BxA(5,i,j)+2.*BxA(1,i,j)
c     $               -2.0*dt*(fflux(i,j,5)-fflux(i-1,j,5))/dx
c     $               -2.0*dt*(gflux(i,j,5)-gflux(i,j-1,5))/dy
c     $               -2.0*dt*source(i,j,5))/3.0
c                   ByA(1,i,j)=(ByA(5,i,j)+2.*ByA(1,i,j)
c     $               -2.0*dt*(fflux(i,j,6)-fflux(i-1,j,6))/dx
c     $               -2.0*dt*(gflux(i,j,6)-gflux(i,j-1,6))/dy
c     $               -2.0*dt*source(i,j,6))/3.0
c                   BzA(1,i,j)=(BzA(5,i,j)+2.*BzA(1,i,j)
c     $               -2.0*dt*(fflux(i,j,7)-fflux(i-1,j,7))/dx
c     $               -2.0*dt*(gflux(i,j,7)-gflux(i,j-1,7))/dy
c     $               -2.0*dt*source(i,j,7))/3.0

                   BxA(1,i,j)=(BxA(5,i,j)+2.*BxA(1,i,j)
     $               -2.0*dt*(Ezf(i,j)-Ezf(i,j-1))/dy
     $               -2.0*dt*source(i,j,5))/3.0
                   ByA(1,i,j)=(ByA(5,i,j)+2.*ByA(1,i,j)
     $               +2.0*dt*(Ezf(i,j)-Ezf(i-1,j))/dx
     $               -2.0*dt*source(i,j,6))/3.0
                   BzA(1,i,j)=(BzA(5,i,j)+2.*BzA(1,i,j)
     $               -2.0*dt*(Eyf(i,j)-Eyf(i-1,j))/dx
     $               +2.0*dt*(Exf(i,j)-Exf(i,j-1))/dy
     $               -2.0*dt*source(i,j,7))/3.0
     
               end do
             end do
             
             end if
c--------------new CTM                     
          end if
      ENDIF

ccccttttttttttttttttttttttttttttttttttt
      IF (CT.EQ.1) then
          SELECT CASE (EXAMPLE)
          CASE(1)

c--------boundary for Orszag-Tang-----------------------------------
          do j=1-blen,NY+blen
             do k=1-blen,0
               BxA(1,k,j)=BxA(1,NX+k-1,j)
               BxA(1,NX-k+1,j)=BxA(1,2-k,j)

               ByA(1,k,j)=ByA(1,NX+k-1,j)
               ByA(1,NX-k+1,j)=ByA(1,2-k,j)

               BzA(1,k,j)=BzA(1,NX+k-1,j)
               BzA(1,NX-k+1,j)=BzA(1,2-k,j)
               
            end do
          end do

          do i=1-blen,NX+blen
            do k=1-blen,0
               BxA(1,i,k)=BxA(1,i,NY+k-1)
               BxA(1,i,NY-k+1)=BxA(1,i,2-k)

               ByA(1,i,k)=ByA(1,i,NY+k-1)
               ByA(1,i,NY-k+1)=ByA(1,i,2-k)

               BzA(1,i,k)=BzA(1,i,NY+k-1)
               BzA(1,i,NY-k+1)=BzA(1,i,2-k)
               
            end do
          end do
c--------End of Boundary for Orszag-Tang-----------------------------------

c--------boundary for 2D Rotor Problem-------------------------------------
          CASE(3)
          do j=1-blen,NY+blen
             do k=1-blen,0
               BxA(1,k,j)=BxA(1,1,j)
               BxA(1,NX-k+1,j)=BxA(1,NX,j)

               ByA(1,k,j)=ByA(1,1,j)
               ByA(1,NX-k+1,j)=ByA(1,NX,j)

               BzA(1,k,j)=BzA(1,1,j)
               BzA(1,NX-k+1,j)=BzA(1,NX,j)
               
             end do             
          end do

          do i=1-blen,NX+blen
             do k=1-blen,0
               BxA(1,i,k)=BxA(1,i,1)
               BxA(1,i,NY-k+1)=BxA(1,i,NY)

               ByA(1,i,k)=ByA(1,i,1)
               ByA(1,i,NY-k+1)=ByA(1,i,NY)

               BzA(1,i,k)=BzA(1,i,1)
               BzA(1,i,NY-k+1)=BzA(1,i,NY)
               
             end do
          end do
c--------End of boundary for 2D Rotor Problem------------------------------
         END SELECT
c==========================================================================      
         do i=1,NX
            do j=1,NY
               temp1=0.5*(BxA(1,i-1,j)+BxA(1,i,j))
               temp2=0.5*(ByA(1,i,j-1)+ByA(1,i,j))
               q(1,i,j,8)=q(1,i,j,8)+
     $                    0.5*(temp1**2+temp2**2+BzA(1,i,j)**2
     $                       -q(1,i,j,5)**2-q(1,i,j,6)**2-q(1,i,j,7)**2)
               q(1,i,j,5)=temp1
               q(1,i,j,6)=temp2
               q(1,i,j,7)=BzA(1,i,j)
            end do
         end do
           
      END IF !IF (CT.EQ.1) then
ccccttttttttttttttttttttttttttttttttttt




            SELECT CASE (EXAMPLE)
            CASE(1)

c--------boundary for Orszag-Tang-----------------------------------
      do j=1-blen,NY+blen
         do k=1-blen,0
            do i1=1,8
               q(1,k,j,i1)=q(1,NX+k-1,j,i1)
               q(1,NX-k+1,j,i1)=q(1,2-k,j,i1)
            end do
          end do
       end do

      do i=1-blen,NX+blen
         do k=1-blen,0
            do i1=1,8
               q(1,i,k,i1)=q(1,i,NY+k-1,i1)
               q(1,i,NY-k+1,i1)=q(1,i,2-k,i1)
            end do
          end do
       end do
c--------End of Boundary for Orszag-Tang-----------------------------------

c--------boundary for 2D Rotor Problem-------------------------------------
            CASE(3)
      do j=1-blen,NY+blen
         do k=1-blen,0
            do i1=1,8
               q(1,k,j,i1)=q(1,1,j,i1)
               q(1,NX-k+1,j,i1)=q(1,NX,j,i1)
            end do
          end do
       end do

      do i=1-blen,NX+blen
         do k=1-blen,0
            do i1=1,8
               q(1,i,k,i1)=q(1,i,1,i1)
               q(1,i,NY-k+1,i1)=q(1,i,NY,i1)
            end do
          end do
       end do
c--------End of boundary for 2D Rotor Problem------------------------------
            END SELECT

c-----------Runge-Kutta Mathod----------------------
      END DO !        IRK=1,NRK


c----------output--------------------------
      if (time-dt.lt.output(1).and.time.ge.output(1)) then
         open(10,file='rho1.dat')
         isave=1
      end if
      if (time-dt.lt.output(2).and.time.ge.output(2)) then
         open(10,file='rho2.dat')
         isave=1
      end if
      if (time-dt.lt.output(3).and.time.ge.output(3)) then
         open(10,file='rho3.dat')
         isave=1
      end if
      if (time-dt.lt.output(4).and.time.ge.output(4)) then
         open(10,file='rho4.dat')
         isave=1
      end if
      if (isave.eq.1) then
c-------output---------------
      WRITE(10,*)'TITLE="FRONT STEP"'
      WRITE(10,*)'VARIABLES=X,Y,RHO,u,v,w,Bx,By,Bz,p'
      WRITE(10,*)'ZONE I=',NY,'J=',NX
      do i=1,NX
         x=(i-1.0)*dx-xstart
         do j=1,NY
            y=(j-1.0)*dy-ystart
            rho(i,j)=q(1,i,j,1)
            u(i,j)=q(1,i,j,2)/rho(i,j)
            v(i,j)=q(1,i,j,3)/rho(i,j)
            w(i,j)=q(1,i,j,4)/rho(i,j)
            B_x(i,j)=q(1,i,j,5)
            B_y(i,j)=q(1,i,j,6)
            B_z(i,j)=q(1,i,j,7)
            B2=B_x(i,j)**2+B_y(i,j)**2+B_z(i,j)**2
            v2=u(i,j)**2+v(i,j)**2+w(i,j)**2
            p(i,j)=(gamma-1.)*(q(1,i,j,8)-0.5*(rho(i,j)*v2+B2))
            write(10,*)x,y,rho(i,j),u(i,j),v(i,j),w(i,j),
     $                 b_x(i,j),b_y(i,j),b_z(i,j),p(i,j)
         end do
      end do
      close(10)
      isave=0
      end if
c-------end of output---------------
         
      i=72
      j=94
      write(*,*)'72,94',q(1,i,j,5),q(1,i,j,6),q(1,i,j,7)



c---------------------------------
      if (time.lt.finaltime) goto 100

      open(10,file='final_t.dat')
      WRITE(10,*)time
      close(10)

c-------output---------------
      open(10,file='rho.dat')
      WRITE(10,*)'TITLE="FRONT STEP"'
      WRITE(10,*)'VARIABLES=X,Y,RHO,u,v,w,Bx,By,Bz,p'
      WRITE(10,*)'ZONE I=',NY,'J=',NX
      do i=1,NX
         x=(i-1.0)*dx
         do j=1,NY
            y=(j-1.0)*dy
            rho(i,j)=q(1,i,j,1)
            u(i,j)=q(1,i,j,2)/rho(i,j)
            v(i,j)=q(1,i,j,3)/rho(i,j)
            w(i,j)=q(1,i,j,4)/rho(i,j)
            B_x(i,j)=q(1,i,j,5)
            B_y(i,j)=q(1,i,j,6)
            B_z(i,j)=q(1,i,j,7)
            B2=B_x(i,j)**2+B_y(i,j)**2+B_z(i,j)**2
            v2=u(i,j)**2+v(i,j)**2+w(i,j)**2
            p(i,j)=(gamma-1.)*(q(1,i,j,8)-0.5*(rho(i,j)*v2+B2))
            write(10,*)x,y,rho(i,j),u(i,j),v(i,j),w(i,j),
     $                 b_x(i,j),b_y(i,j),b_z(i,j),p(i,j)
         end do
      end do
      close(10)
     
      open(11,file='rholine.dat')
       WRITE(11,*)'TITLE="FRONT STEP"'
      WRITE(11,*)'VARIABLES=X,RHO,u,v,w,Bx,By,Bz,p'
      WRITE(11,*)'ZONE I=',NY
      do i=1,NX
         x=(i-1.0)*dx
            j=i
            y=(j-1.0)*dy
            rho(i,j)=q(1,i,j,1)
            u(i,j)=q(1,i,j,2)/rho(i,j)
            v(i,j)=q(1,i,j,3)/rho(i,j)
            w(i,j)=q(1,i,j,4)/rho(i,j)
            B_x(i,j)=q(1,i,j,5)
            B_y(i,j)=q(1,i,j,6)
            B_z(i,j)=q(1,i,j,7)
            B2=B_x(i,j)**2+B_y(i,j)**2+B_z(i,j)**2
            v2=u(i,j)**2+v(i,j)**2+w(i,j)**2
            p(i,j)=(gamma-1.)*(q(1,i,j,8)-0.5*(rho(i,j)*v2+B2))
            write(11,*)x,rho(i,j),u(i,j),v(i,j),w(i,j),
     $                 b_x(i,j),b_y(i,j),b_z(i,j),p(i,j)
      end do
      close(11)
      end
        subroutine reconstruct(scheme,N,q,ql,qr)
      implicit none
      integer, intent(in):: scheme, N
      double precision, intent(in):: q(-5:N+5) 
      double precision, intent(out):: ql(-5:N+5),qr(-5:N+5)
      double precision is0,is1,is2,
     $                  q0,q1,q2,
     $                  a0,a1,a2,
     $                  w0,w1,w2
      double precision eps,ss
      double precision tao5,tao0,tao2,dmin,tao1
      integer i,j,k,m
      INTEGER START_P(1000),END_P(1000)
      double precision::A(1000),B(1000),C(1000),
     &D(1000),XH(1000),H(1000),HR(1000)
      integer ntc
      eps=1e-2
      ss=1e-2
      select case(scheme)
c----------- first order
      case (1)
        do i=0,N
           ql(i)=q(i)
           qr(i)=q(i+1)
        end do
c----------- second order TVD
      case (2)
c----------- third order WENO
      case (3)
        do i=-1,N+1
           q0=(-q(i-1)+3.*q(i))/2.
           q1=0.5*(q(i)+q(i+1))

           is0=(q(i-1)-q(i))**2
           is1=(q(i)-q(i+1))**2

           a0=1./3./(eps+is0)**2
           a1=2./3./(eps+is1)**2

           w0=a0/(a0+a1)
           w1=a1/(a0+a1)
           
           ql(i)=w0*q0+w1*q1

c------
           q0=(-q(i+2)+3.*q(i+1))/2.
           q1=0.5*(q(i+1)+q(i))

           is0=(q(i+2)-q(i+1))**2
           is1=(q(i+1)-q(i))**2

           a0=1./3./(eps+is0)**2
           a1=2./3./(eps+is1)**2

           w0=a0/(a0+a1)
           w1=a1/(a0+a1)
           
           qr(i)=w0*q0+w1*q1
        end do

c----------- fifth order WENO
      case (5)
        do i=-1,N+1
          is0=13./12.*(q(i-2)-2.*q(i-1)+q(i))**2+
     $                0.25*(q(i-2)-4.*q(i-1)+3.*q(i))**2
           is1=13./12.*(q(i-1)-2.*q(i)+q(i+1))**2+
     $                0.25*(q(i-1)-q(i+1))**2
           is2=13./12.*(q(i)-2.*q(i+1)+q(i+2))**2+
     $                0.25*(3.*q(i)-4.*q(i+1)+q(i+2))**2
           
           tao5=abs(is0-is2)
           
           q0=(2.*q(i-2)-7.*q(i-1)+11.*q(i))/6.0
           q1=(-q(i-1)+5.*q(i)+2.*q(i+1))/6.0
           q2=(2.*q(i)+5.*q(i+1)-q(i+2))/6.0
          
                a0=0.1*(1.0d0+(tao5/(ss+IS0)))
                a1=0.6*(1.0d0+(tao5/(ss+IS1)))
                a2=0.3*(1.0d0+(tao5/(ss+IS2)))          
                w0=a0/(a0+a1+a2)
                w1=a1/(a0+a1+a2)
                w2=a2/(a0+a1+a2)
     
           ql(i)=w0*q0+w1*q1+w2*q2

c---
           q0=(2.*q(i+3)-7.*q(i+2)+11.*q(i+1))/6.0
           q1=(-q(i+2)+5.*q(i+1)+2.*q(i))/6.0
           q2=(2.*q(i+1)+5.*q(i)-q(i-1))/6.0

           is0=13./12.*(q(i+3)-2.*q(i+2)+q(i+1))**2+
     $                0.25*(q(i+3)-4.*q(i+2)+3.*q(i+1))**2
           is1=13./12.*(q(i+2)-2.*q(i+1)+q(i))**2+
     $                0.25*(q(i+2)-q(i))**2
           is2=13./12.*(q(i+1)-2.*q(i)+q(i-1))**2+
     $                0.25*(3.*q(i+1)-4.*q(i)+q(i-1))**2
             tao5=abs(is0-is2)
          
                a0=0.1*(1.0d0+(tao5/(ss+IS0)))
                a1=0.6*(1.0d0+(tao5/(ss+IS1)))
                a2=0.3*(1.0d0+(tao5/(ss+IS2)))          
                w0=a0/(a0+a1+a2)
                w1=a1/(a0+a1+a2)
                w2=a2/(a0+a1+a2)
           qr(i)=w0*q0+w1*q1+w2*q2
        end do
        
c----------- fifth order upstream scheme (fixed points)
      case (6)
    
        do i=-1,N+1
           ql(i)=(2.*q(i-2)-13.*q(i-1)+47.*q(i)
     $           +27.*q(i+1)-3.*q(i+2))/60.
           qr(i)=(2.*q(i+3)-13.*q(i+2)+47.*q(i+1)
     $           +27.*q(i)-3.*q(i-1))/60.           
        end do
      case(56)
      
      do j=1,2
      if (j.eq.1) i=-1
      if (j.eq.2) i=N+1
           is0=13./12.*(q(i-2)-2.*q(i-1)+q(i))**2+
     $                0.25*(q(i-2)-4.*q(i-1)+3.*q(i))**2
           is1=13./12.*(q(i-1)-2.*q(i)+q(i+1))**2+
     $                0.25*(q(i-1)-q(i+1))**2
           is2=13./12.*(q(i)-2.*q(i+1)+q(i+2))**2+
     $                0.25*(3.*q(i)-4.*q(i+1)+q(i+2))**2
           
           tao5=abs(is0-is2)
           
           q0=(2.*q(i-2)-7.*q(i-1)+11.*q(i))/6.0
           q1=(-q(i-1)+5.*q(i)+2.*q(i+1))/6.0
           q2=(2.*q(i)+5.*q(i+1)-q(i+2))/6.0
          
                a0=0.1*(1.0d0+(tao5/(ss+IS0)))
                a1=0.6*(1.0d0+(tao5/(ss+IS1)))
                a2=0.3*(1.0d0+(tao5/(ss+IS2)))          
                w0=a0/(a0+a1+a2)
                w1=a1/(a0+a1+a2)
                w2=a2/(a0+a1+a2)
!                ql(i)=w0*Q0+w1*Q31+w2*Q32
     
           ql(i)=w0*q0+w1*q1+w2*q2
      enddo
!      i=0
!      ql(i)=(-2.*q(i-1)+22.0*q(i)+57.*q(i+1)-
!     $         23.*q(i+2)+7.0*q(i+3)-q(i+4))/60.0d0
!      i=n+1
!       ql(i)=(-q(i-3)+7.0*q(i-2)-23.*q(i-1)+
!     $         57.*q(i)+22.0*q(i+1)-2.*q(i+2))/60.0d0
     
      start_p=N+4
       NTC=1
       start_p(NTC)=0
       do I=0,N
          IS0=13.d0*(q(I-2)-2.*q(I-1)+q(I))**2/12.d0
     $                     +(q(I-2)-4.*q(I-1)+3.*q(I))**2/4.d0
          IS1=13.d0*(q(I-1)-2.*q(I)+q(I+1))**2/12.d0
     $                     +(q(I-1)-q(I+1))**2/4.d0
          IS2=13.d0*(q(I)-2.*q(I+1)+q(I+2))**2/12.d0
     $                     +(3.*q(I)-4.*q(I+1)+q(I+2))**2/4.d0
          tao5=abs(is0-is2)
          tao0=abs(is0-is1)
          tao1=abs(is2-is1)
          dmin=dmin1(is0,is1,is2)
          if (tao5.gt.dmin) then
!            if(tao5.lt.0)then         
             end_p(NTC)=i
             NTC=NTC+1
             start_p(NTC)=i+1

             q0=(2.*q(i-2)-7.*q(i-1)+11.*q(i))/6.0
             q1=(-q(i-1)+5.*q(i)+2.*q(i+1))/6.0
             q2=(2.*q(i)+5.*q(i+1)-q(i+2))/6.0
             
              a0=0.1*(1.0d0+(tao5/(ss+IS0)))
              a1=0.6*(1.0d0+(tao5/(ss+IS1)))
              a2=0.3*(1.0d0+(tao5/(ss+IS2))) 
             
                w0=a0/(a0+a1+a2)
                w1=a1/(a0+a1+a2)
                w2=a2/(a0+a1+a2)
                ql(i)=w0*Q0+w1*Q1+w2*Q2
           
                    else 
!            ql(i)=1.d0/60*q(i+3)-2.d0/15*q(i+2)
!     &          +37.d0/60*q(i+1)+1.d0/60*q(i-2)
!     &          -2.d0/15*q(i-1)+37.d0/60*q(i)

!                 qr(i)=w0*Q0+w1*Q1+w2*Q2
!          end if
       end if
       end do            ! end of (do i=5,n-5)
       end_p(NTC)=N+1
      do k=1,ntc
          m=end_p(k)-start_p(k)

          if (m.ge.1) then
             do j=1,m
                i=start_p(k)+j-1
                a(j)=1.0d0/3.0d0
                b(j)=1.0d0
                c(j)=a(j)
                d(j)=(q(i+2)+29.d0*q(i+1)+29.d0*q(i)+q(i-1))
     $            /36.d0
             end do
             d(1)=d(1)-a(1)*ql(start_p(k)-1)
             d(m)=d(m)-c(m)*ql(end_p(k))
             call tridiagsolve(m,a,b,c,d,xh)
             do j=1,m
                ql(start_p(k)+j-1)=xh(j)
             end do
           
          end if
        enddo            
c****************************************************************
c****************************************************************
c********************Negative PART*******************************
       do j=1,2
c       do j=4,n-4
          if (j.eq.1) i=-1
          if (j.eq.2) i=N+1
c          i=j
          IS0=13.*(q(I+3)-2.*q(I+2)+q(I+1))**2/12.
     $                     +(q(I+3)-4.*q(I+2)+3.*q(I+1))**2/4.
          IS1=13.*(q(I+2)-2.*q(I+1)+q(I))**2/12.
     $                     +(q(I+2)-q(I))**2/4.
          IS2=13.*(q(I+1)-2.*q(I)+q(I-1))**2/12.
     $                     +(3.*q(I+1)-4.*q(I)+q(I-1))**2/4.

c---------------
          tao5=abs(is0-is2)

           q0=(2.*q(i+3)-7.*q(i+2)+11.*q(i+1))/6.0
           q1=(-q(i+2)+5.*q(i+1)+2.*q(i))/6.0
           q2=(2.*q(i+1)+5.*q(i)-q(i-1))/6.0
      

            a0=0.1*(1.0d0+(tao5/(ss+IS0)))
            a1=0.6*(1.0d0+(tao5/(ss+IS1)))
            a2=0.3*(1.0d0+(tao5/(ss+IS2)))          

            w0=a0/(a0+a1+a2)
            w1=a1/(a0+a1+a2)
            w2=a2/(a0+a1+a2)

            qr(i)=w0*Q0+w1*Q1+w2*Q2
         end do
!
!
!
c-----------------6th-order biased  h(4) and h(n-4)--------
!       i=0
!       qr(i)=(-2.*q(i-1)+22.0*q(i)+57.*q(i+1)-
!     $         23.*q(i+2)+7.0*q(i+3)-q(i+4))/60.0d0
!       i=n
!       qr(i)=(-q(i-3)+7.0*q(i-2)-23.*q(i-1)+
!     $         57.*q(i)+22.0*q(i+1)-2.*q(i+2))/60.0d0

c       goto 2345

       start_p=N+4
       NTC=1
       start_p(NTC)=0
       do I=0,N
          IS0=13.*(q(I+3)-2.*q(I+2)+q(I+1))**2/12.
     $                     +(q(I+3)-4.*q(I+2)+3.*q(I+1))**2/4.
          IS1=13.*(q(I+2)-2.*q(I+1)+q(I))**2/12.
     $                     +(q(I+2)-q(I))**2/4.
          IS2=13.*(q(I+1)-2.*q(I)+q(I-1))**2/12.
     $                     +(3.*q(I+1)-4.*q(I)+q(I-1))**2/4.
          tao5=abs(is0-is2)
          tao0=abs(is0-is1)
          tao1=abs(is2-is1)
          dmin=dmin1(is0,is1,is2)
          if (tao5.gt.dmin) then 
!            if(tao5.lt.0)then        
             end_p(NTC)=i
             NTC=NTC+1
             start_p(NTC)=i+1

           q0=(2.*q(i+3)-7.*q(i+2)+11.*q(i+1))/6.0
           q1=(-q(i+2)+5.*q(i+1)+2.*q(i))/6.0
           q2=(2.*q(i+1)+5.*q(i)-q(i-1))/6.0
      

            a0=0.1*(1.0d0+(tao5/(ss+IS0)))
            a1=0.6*(1.0d0+(tao5/(ss+IS1)))
            a2=0.3*(1.0d0+(tao5/(ss+IS2)))          

            w0=a0/(a0+a1+a2)
            w1=a1/(a0+a1+a2)
            w2=a2/(a0+a1+a2)

            qr(i)=w0*Q0+w1*Q1+w2*Q2
            else 
!            qr(i)=1.d0/60*q(i+3)-2.d0/15*q(i+2)
!     &          +37.d0/60*q(i+1)+1.d0/60*q(i-2)
!     &          -2.d0/15*q(i-1)+37.d0/60*q(i)
!                 qr(i)=w0*Q0+w1*Q1+w2*Q2
          end if
       end do            ! end of (do i=5,n-5)
       end_p(NTC)=N+1

       do k=1,NTC

          m=end_p(k)-start_p(k)

          if (m.ge.1) then
             do j=1,m
                i=start_p(k)+j-1
                a(j)=1.0d0/3.0d0
                b(j)=1.0d0
                c(j)=a(j)
                d(j)=(q(i+2)+29.d0*q(i+1)+29.d0*q(i)+q(i-1))
     $            /36.d0
             end do
             d(1)=d(1)-a(1)*qr(start_p(k)-1)
             d(m)=d(m)-c(m)*qr(end_p(k))
             call tridiagsolve(m,a,b,c,d,xh)
             do j=1,m
                qr(start_p(k)+j-1)=xh(j)
             end do
             
          end if
       end do                 


      end select
     
      return
        end
        
        
      subroutine tridiagsolve(n,a,b,c,d,x)
c     Forward sweep and back substitution algorithm for
c     tridiagonal matrix
        implicit none
        integer, intent(in)::
     $        n                  ! dimension
        double precision, dimension(n)::
     $        a,b,c,d            !a(i)*x(i-1)+b(i)*x(i)+c(i)*x(i+1)=d(i)
        double precision, dimension(n),intent(out)::
     $        x 
c       local variable
        integer::
     $        i       
c             start
c             Modify the coefficients
        c(1)=c(1)/b(1)
        d(1)=d(1)/b(1)
        do i=2,n
           c(i)=c(i)/(b(i)-c(i-1)*a(i))
           d(i)=(d(i)-d(i-1)*a(i))/(b(i)-c(i-1)*a(i))
        end do
c             Back substitute
        x(n)=d(n)
        do i=n-1,1,-1
           x(i)=d(i)-c(i)*x(i+1)
        end do

        return

        end

         
!      subroutine reconstruct(scheme,N,q,ql,qr)
!      implicit DOUBLE PRECISION(a-h,l,o-z)
!      PARAMETER(NN=1689)
!       PARAMETER (NX=201,NY=201,MP=201)
!      integer, intent(in):: scheme, N
!      double precision, intent(in):: q(-5:N+5) 
!      double precision, intent(out):: ql(-5:N+5),qr(-5:N+5)
!      double precision is0,is1,is2,
!     $                  q0,q1,q2,
!     $                  a0,a1,a2,
!     $                  w0,w1,w2
!      double precision temp
!      double precision eps
!      integer i
!      INTEGER START_P(1000),END_P(1000)
!       double precision::A(NN),B(NN),C(NN),D(NN),XH(NN),H(NN),HR(NN)
!      
!      eps=1e-10
!      C40=1.D0/35.D0
!	C41=12.D0/35.D0
!	C42=18.D0/35.D0
!        C43=4.D0/35.D0
!        eps=1e-2
!         SS=1E-20
!	   SSS=SS
!         a300=1.d0/3.d0
!         a301=-7.d0/6.d0
!         a302=11.d0/6.d0
!         a310=-1.d0/6.d0
!         a311=5.d0/6.d0
!         a312=1.d0/3.d0
!         a320=1.d0/3.d0
!         a321=5.d0/6.d0
!         a322=-1.d0/6.d0
!         c30=1.d0/10.d0
!         c31=6.d0/10.d0
!         c32=3.d0/10.d0
!       dx=2.0*pi/(NX-1.0)
!       dy=2.0*pi/(NY-1.0)
!
!      select case(scheme)
!c----------- Ziegler TVD--------------------
!!        NOTE: ql(i)=====>q^E(i)
!!              qr(i)=====>q^W(i+1)
!      case (0)
!        do i=-1,N+1
!           temp=q(i+1)-q(i-1)
!           if (temp.eq.0) then
!               ql(i)=q(i)
!           else
!               ql(i)=q(i)+max((q(i+1)-q(i))*(q(i)-q(i-1)),
!     $               0.0d0)/temp
!           end if
!           temp=q(i+2)-q(i)
!           if (temp.eq.0) then
!               qr(i)=q(i+1)
!           else
!               qr(i)=q(i+1)-max((q(i+2)-q(i+1))*(q(i+1)-q(i)),
!     $               0.0d0)/temp
!           end if
!        end do
!
!
!
!c----------- first order
!      case (1)
!        do i=-1,N+1
!           ql(i)=q(i)
!           qr(i)=q(i+1)
!        end do
!c----------- second order TVD
!      case (2)
!c----------- third order WENO
!      case (3)
!        do i=-1,N+1
!           q0=(-q(i-1)+3.*q(i))/2.
!           q1=0.5*(q(i)+q(i+1))
!
!           is0=(q(i-1)-q(i))**2
!           is1=(q(i)-q(i+1))**2
!
!           a0=1./3./(eps+is0)**2
!           a1=2./3./(eps+is1)**2
!
!           w0=a0/(a0+a1)
!           w1=a1/(a0+a1)
!           
!           ql(i)=w0*q0+w1*q1
!
!c------
!           q0=(-q(i+2)+3.*q(i+1))/2.
!           q1=0.5*(q(i+1)+q(i))
!
!           is0=(q(i+2)-q(i+1))**2
!           is1=(q(i+1)-q(i))**2
!
!           a0=1./3./(eps+is0)**2
!           a1=2./3./(eps+is1)**2
!
!           w0=a0/(a0+a1)
!           w1=a1/(a0+a1)
!           
!           qr(i)=w0*q0+w1*q1
!        end do
!
!c----------- fifth order WENO
!      case (5)
!        do i=-1,N+1
!           q0=(2.*q(i-2)-7.*q(i-1)+11.*q(i))/6.0
!           q1=(-q(i-1)+5.*q(i)+2.*q(i+1))/6.0
!           q2=(2.*q(i)+5.*q(i+1)-q(i+2))/6.0
!
!           is0=13./12.*(q(i-2)-2.*q(i-1)+q(i))**2+
!     $                0.25*(q(i-2)-4.*q(i-1)+3.*q(i))**2
!           is1=13./12.*(q(i-1)-2.*q(i)+q(i+1))**2+
!     $                0.25*(q(i-1)-q(i+1))**2
!           is2=13./12.*(q(i)-2.*q(i+1)+q(i+2))**2+
!     $                0.25*(3.*q(i)-4.*q(i+1)+q(i+2))**2
!
!           a0=0.1/(eps+is0)**2
!           a1=0.6/(eps+is1)**2
!           a2=0.3/(eps+is2)**2
!
!           w0=a0/(a0+a1+a2)
!           w1=a1/(a0+a1+a2)
!           w2=a2/(a0+a1+a2)
!     
!           ql(i)=w0*q0+w1*q1+w2*q2
!
!c---
!           q0=(2.*q(i+3)-7.*q(i+2)+11.*q(i+1))/6.0
!           q1=(-q(i+2)+5.*q(i+1)+2.*q(i))/6.0
!           q2=(2.*q(i+1)+5.*q(i)-q(i-1))/6.0
!
!           is0=13./12.*(q(i+3)-2.*q(i+2)+q(i+1))**2+
!     $                0.25*(q(i+3)-4.*q(i+2)+3.*q(i+1))**2
!           is1=13./12.*(q(i+2)-2.*q(i+1)+q(i))**2+
!     $                0.25*(q(i+2)-q(i))**2
!           is2=13./12.*(q(i+1)-2.*q(i)+q(i-1))**2+
!     $                0.25*(3.*q(i+1)-4.*q(i)+q(i-1))**2
!     
!          
!          tao=abs(is0-is2)
!
!           a0=0.1/(eps+is0)**2
!           a1=0.6/(eps+is1)**2
!           a2=0.3/(eps+is2)**2
!           
!           
!!             a0=C30*(1.0d0+(tao/(ss+IS0)))
!!            a1=C31*(1.0d0+(tao/(ss+IS1)))
!!            a2=C32*(1.0d0+(tao/(ss+IS2)))   
!
!           w0=a0/(a0+a1+a2)
!           w1=a1/(a0+a1+a2)
!           w2=a2/(a0+a1+a2)
!     
!           qr(i)=w0*q0+w1*q1+w2*q2
!        end do
!
!
!
!      case(55)
!      
!      
!      
!      do j=1,2
!      if (j.eq.1) i=-1
!      if (j.eq.2) i=N+1
!           is0=13./12.*(q(i-2)-2.*q(i-1)+q(i))**2+
!     $                0.25*(q(i-2)-4.*q(i-1)+3.*q(i))**2
!           is1=13./12.*(q(i-1)-2.*q(i)+q(i+1))**2+
!     $                0.25*(q(i-1)-q(i+1))**2
!           is2=13./12.*(q(i)-2.*q(i+1)+q(i+2))**2+
!     $                0.25*(3.*q(i)-4.*q(i+1)+q(i+2))**2
!           
!           tao5=abs(is0-is2)
!           
!           q0=(2.*q(i-2)-7.*q(i-1)+11.*q(i))/6.0
!           q1=(-q(i-1)+5.*q(i)+2.*q(i+1))/6.0
!           q2=(2.*q(i)+5.*q(i+1)-q(i+2))/6.0
!          
!           a0=0.1/(eps+is0)**2
!           a1=0.6/(eps+is1)**2
!           a2=0.3/(eps+is2)**2
!
!           w0=a0/(a0+a1+a2)
!           w1=a1/(a0+a1+a2)
!           w2=a2/(a0+a1+a2)
!     
!           ql(i)=w0*q0+w1*q1+w2*q2
!      enddo
!      i=-1
!      ql(i)=(-2.*q(i-1)+22.0*q(i)+57.*q(i+1)-
!     $         23.*q(i+2)+7.0*q(i+3)-q(i+4))/60.0d0
!      i=n+1
!       ql(i)=(-q(i-3)+7.0*q(i-2)-23.*q(i-1)+
!     $         57.*q(i)+22.0*q(i+1)-2.*q(i+2))/60.0d0
!     
!      start_p=N+5
!       NTC=1
!       start_p(NTC)=0
!       do I=0,N
!          IS0=13.d0*(q(I-2)-2.*q(I-1)+q(I))**2/12.d0
!     $                     +(q(I-2)-4.*q(I-1)+3.*q(I))**2/4.d0
!          IS1=13.d0*(q(I-1)-2.*q(I)+q(I+1))**2/12.d0
!     $                     +(q(I-1)-q(I+1))**2/4.d0
!          IS2=13.d0*(q(I)-2.*q(I+1)+q(I+2))**2/12.d0
!     $                     +(3.*q(I)-4.*q(I+1)+q(I+2))**2/4.d0
!          tao5=abs(is0-is2)
!          tao0=abs(is0-is1)
!          tao1=abs(is2-is1)
!          dmin=dmin1(is0,is1,is2)
!          
!!          tao5=abs(q(i+1)-q(i))
!!          dmin=dmax1(dx,dy)
!          
!          if (tao5.gt.dmin) then         
!             end_p(NTC)=i
!             NTC=NTC+1
!             start_p(NTC)=i+1
!
!             Q30=A300*q(i-2)+A301*q(i-1)+A302*q(i)
!             Q31=A310*q(i-1)+A311*q(i)+A312*q(i+1)
!             Q32=A320*q(i)+A321*q(i+1)+A322*q(i+2)
!
!
!                aa0=C30*(1.0d0+(tao5/(ss+IS0)))
!                aa1=C31*(1.0d0+(tao5/(ss+IS1)))
!                aa2=C32*(1.0d0+(tao5/(ss+IS2)))          
!                w0=aa0/(aa0+aa1+aa2)
!                w1=aa1/(aa0+aa1+aa2)
!                w2=aa2/(aa0+aa1+aa2)
!                ql(i)=w0*Q30+w1*Q31+w2*Q32
!           
!        
!       end if
!       end do            ! end of (do i=5,n-5)
!       end_p(NTC)=N+1
!      do k=1,NTC
!          m=end_p(k)-start_p(k)
!
!          if (m.ge.1) then
!             do j=1,m
!                i=start_p(k)+j-1
!                a(j)=1.0d0/3.0d0
!                b(j)=1.0d0
!                c(j)=a(j)
!                d(j)=(q(i+2)+29.d0*q(i+1)+29.d0*q(i)+q(i-1))
!     $            /36.d0
!             end do
!             d(1)=d(1)-a(1)*ql(start_p(k)-1)
!             d(m)=d(m)-c(m)*ql(end_p(k))
!             call tridiagsolve(m,a,b,c,d,xh)
!             do j=1,m
!                ql(start_p(k)+j-1)=xh(j)
!             end do
!           
!          end if
!        enddo            
!c****************************************************************
!c****************************************************************
!c********************Negative PART*******************************
!       do j=1,2
!c       do j=4,n-4
!          if (j.eq.1) i=-1
!          if (j.eq.2) i=N+1
!c          i=j
!          IS0=13.*(q(I+3)-2.*q(I+2)+q(I+1))**2/12.
!     $                     +(q(I+3)-4.*q(I+2)+3.*q(I+1))**2/4.
!          IS1=13.*(q(I+2)-2.*q(I+1)+q(I))**2/12.
!     $                     +(q(I+2)-q(I))**2/4.
!          IS2=13.*(q(I+1)-2.*q(I)+q(I-1))**2/12.
!     $                     +(3.*q(I+1)-4.*q(I)+q(I-1))**2/4.
!
!c---------------
!          tao5=abs(is0-is2)
!
!          Q30=A300*q(i+3)+A301*q(i+2)+A302*q(i+1)
!          Q31=A310*q(i+2)+A311*q(i+1)+A312*q(i)
!          Q32=A320*q(i+1)+A321*q(i)+A322*q(i-1)
!      
!
!            aa0=C30*(1.0d0+(tao5/(ss+IS0)))
!            aa1=C31*(1.0d0+(tao5/(ss+IS1)))
!            aa2=C32*(1.0d0+(tao5/(ss+IS2)))          
!
!            w0=aa0/(aa0+aa1+aa2)
!            w1=aa1/(aa0+aa1+aa2)
!            w2=aa2/(aa0+aa1+aa2)
!
!            qr(i)=w0*Q30+w1*Q31+w2*Q32
!         end do
!
!
!
!c-----------------6th-order biased  h(4) and h(n-4)--------
!       i=-1
!       qr(i)=(-2.*q(i-1)+22.0*q(i)+57.*q(i+1)-
!     $         23.*q(i+2)+7.0*q(i+3)-q(i+4))/60.0d0
!       i=n+1
!       qr(i)=(-q(i-3)+7.0*q(i-2)-23.*q(i-1)+
!     $         57.*q(i)+22.0*q(i+1)-2.*q(i+2))/60.0d0
!
!c       goto 2345
!
!       start_p=N+5
!       NTC=1
!       start_p(NTC)=0
!       do I=0,N
!          IS0=13.*(q(I+3)-2.*q(I+2)+q(I+1))**2/12.
!     $                     +(q(I+3)-4.*q(I+2)+3.*q(I+1))**2/4.
!          IS1=13.*(q(I+2)-2.*q(I+1)+q(I))**2/12.
!     $                     +(q(I+2)-q(I))**2/4.
!          IS2=13.*(q(I+1)-2.*q(I)+q(I-1))**2/12.
!     $                     +(3.*q(I+1)-4.*q(I)+q(I-1))**2/4.
!          tao5=abs(is0-is2)
!          tao0=abs(is0-is1)
!          tao1=abs(is2-is1)
!          dmin=dmin1(is0,is1,is2)
!!          
!!          tao5=abs(q(i+1)-q(i))
!!          dmin=dmax1(dx,dy)
!          
!          if (tao5.gt.dmin) then         
!             end_p(NTC)=i
!             NTC=NTC+1
!             start_p(NTC)=i+1
!
!             Q30=A300*q(i+3)+A301*q(i+2)+A302*q(i+1)
!             Q31=A310*q(i+2)+A311*q(i+1)+A312*q(i)
!             Q32=A320*q(i+1)+A321*q(i)+A322*q(i-1)
!       
!
!                 aa0=C30*(1.0d0+(tao5/(ss+IS0)))
!                 aa1=C31*(1.0d0+(tao5/(ss+IS1)))
!                 aa2=C32*(1.0d0+(tao5/(ss+IS2)))          
!
!                 w0=aa0/(aa0+aa1+aa2)
!                 w1=aa1/(aa0+aa1+aa2)
!                 w2=aa2/(aa0+aa1+aa2)
!
!                 qr(i)=w0*Q30+w1*Q31+w2*Q32
!          end if
!       end do            ! end of (do i=5,n-5)
!       end_p(NTC)=N+1
!
!       do k=1,NTC
!
!          m=end_p(k)-start_p(k)
!
!          if (m.ge.1) then
!             do j=1,m
!                i=start_p(k)+j-1
!                a(j)=1.0d0/3.0d0
!                b(j)=1.0d0
!                c(j)=a(j)
!                d(j)=(q(i+2)+29.d0*q(i+1)+29.d0*q(i)+q(i-1))
!     $            /36.d0
!             end do
!             d(1)=d(1)-a(1)*qr(start_p(k)-1)
!             d(m)=d(m)-c(m)*qr(end_p(k))
!             call tridiagsolve(m,a,b,c,d,xh)
!             do j=1,m
!                qr(start_p(k)+j-1)=xh(j)
!             end do
!          end if
!       end do                 
!           
!      
!      case(56)
!      
!      
!      
!      do j=1,2
!      if (j.eq.1) i=-1
!      if (j.eq.2) i=N+1
!           is0=13./12.*(q(i-2)-2.*q(i-1)+q(i))**2+
!     $                0.25*(q(i-2)-4.*q(i-1)+3.*q(i))**2
!           is1=13./12.*(q(i-1)-2.*q(i)+q(i+1))**2+
!     $                0.25*(q(i-1)-q(i+1))**2
!           is2=13./12.*(q(i)-2.*q(i+1)+q(i+2))**2+
!     $                0.25*(3.*q(i)-4.*q(i+1)+q(i+2))**2
!           
!           tao5=abs(is0-is2)
!           
!           q0=(2.*q(i-2)-7.*q(i-1)+11.*q(i))/6.0
!           q1=(-q(i-1)+5.*q(i)+2.*q(i+1))/6.0
!           q2=(2.*q(i)+5.*q(i+1)-q(i+2))/6.0
!          
!           a0=0.1/(eps+is0)**2
!           a1=0.6/(eps+is1)**2
!           a2=0.3/(eps+is2)**2
!
!           w0=a0/(a0+a1+a2)
!           w1=a1/(a0+a1+a2)
!           w2=a2/(a0+a1+a2)
!     
!           ql(i)=w0*q0+w1*q1+w2*q2
!      enddo
!      i=-1
!      ql(i)=(-2.*q(i-1)+22.0*q(i)+57.*q(i+1)-
!     $         23.*q(i+2)+7.0*q(i+3)-q(i+4))/60.0d0
!      i=n+1
!       ql(i)=(-q(i-3)+7.0*q(i-2)-23.*q(i-1)+
!     $         57.*q(i)+22.0*q(i+1)-2.*q(i+2))/60.0d0
!     
!      start_p=N+5
!       NTC=1
!       start_p(NTC)=0
!       do I=0,N
!          IS0=13.d0*(q(I-2)-2.*q(I-1)+q(I))**2/12.d0
!     $                     +(q(I-2)-4.*q(I-1)+3.*q(I))**2/4.d0
!          IS1=13.d0*(q(I-1)-2.*q(I)+q(I+1))**2/12.d0
!     $                     +(q(I-1)-q(I+1))**2/4.d0
!          IS2=13.d0*(q(I)-2.*q(I+1)+q(I+2))**2/12.d0
!     $                     +(3.*q(I)-4.*q(I+1)+q(I+2))**2/4.d0
!          tao5=abs(is0-is2)
!          tao0=abs(is0-is1)
!          tao1=abs(is2-is1)
!          dmin=dmin1(is0,is1,is2)
!          if (tao5.gt.dmin) then         
!             end_p(NTC)=i
!             NTC=NTC+1
!             start_p(NTC)=i+1
!
!             Q30=A300*q(i-2)+A301*q(i-1)+A302*q(i)
!             Q31=A310*q(i-1)+A311*q(i)+A312*q(i+1)
!             Q32=A320*q(i)+A321*q(i+1)+A322*q(i+2)
!
!
!                aa0=C30*(1.0d0+(tao5/(ss+IS0)))
!                aa1=C31*(1.0d0+(tao5/(ss+IS1)))
!                aa2=C32*(1.0d0+(tao5/(ss+IS2)))          
!                w0=aa0/(aa0+aa1+aa2)
!                w1=aa1/(aa0+aa1+aa2)
!                w2=aa2/(aa0+aa1+aa2)
!                ql(i)=w0*Q30+w1*Q31+w2*Q32
!           
!        
!       else
!           
!             Q30=A300*q(i-2)+A301*q(i-1)+A302*q(i)
!             Q31=A310*q(i-1)+A311*q(i)+A312*q(i+1)
!             Q32=A320*q(i)+A321*q(i+1)+A322*q(i+2)
!       
!                ql(i)=C30*Q30+C31*Q31+C32*Q32      
!                
!                ql(i)=(2.0d0*q(i-2)-13.0d0*q(i-1)+47.0d0*q(i)
!     $                +27.0d0*q(i+1)-3.d0*q(i+2))/60.d0     
!            
!           end if          
!                     ! end of (do i=5,n-5)
!
!       end do            ! end of (do k=1,NTC)     
!c****************************************************************
!c****************************************************************
!c********************Negative PART*******************************
!       do j=1,2
!c       do j=4,n-4
!          if (j.eq.1) i=-1
!          if (j.eq.2) i=N+1
!c          i=j
!          IS0=13.*(q(I+3)-2.*q(I+2)+q(I+1))**2/12.
!     $                     +(q(I+3)-4.*q(I+2)+3.*q(I+1))**2/4.
!          IS1=13.*(q(I+2)-2.*q(I+1)+q(I))**2/12.
!     $                     +(q(I+2)-q(I))**2/4.
!          IS2=13.*(q(I+1)-2.*q(I)+q(I-1))**2/12.
!     $                     +(3.*q(I+1)-4.*q(I)+q(I-1))**2/4.
!
!c---------------
!          tao5=abs(is0-is2)
!
!          Q30=A300*q(i+3)+A301*q(i+2)+A302*q(i+1)
!          Q31=A310*q(i+2)+A311*q(i+1)+A312*q(i)
!          Q32=A320*q(i+1)+A321*q(i)+A322*q(i-1)
!      
!
!            aa0=C30*(1.0d0+(tao5/(ss+IS0)))
!            aa1=C31*(1.0d0+(tao5/(ss+IS1)))
!            aa2=C32*(1.0d0+(tao5/(ss+IS2)))          
!
!            w0=aa0/(aa0+aa1+aa2)
!            w1=aa1/(aa0+aa1+aa2)
!            w2=aa2/(aa0+aa1+aa2)
!
!            qr(i)=w0*Q30+w1*Q31+w2*Q32
!         end do
!
!
!
!c-----------------6th-order biased  h(4) and h(n-4)--------
!       i=-1
!       qr(i)=(-2.*q(i-1)+22.0*q(i)+57.*q(i+1)-
!     $         23.*q(i+2)+7.0*q(i+3)-q(i+4))/60.0d0
!       i=n+1
!       qr(i)=(-q(i-3)+7.0*q(i-2)-23.*q(i-1)+
!     $         57.*q(i)+22.0*q(i+1)-2.*q(i+2))/60.0d0
!
!c       goto 2345
!
!       start_p=N+5
!       NTC=1
!       start_p(NTC)=0
!       do I=0,N
!          IS0=13.*(q(I+3)-2.*q(I+2)+q(I+1))**2/12.
!     $                     +(q(I+3)-4.*q(I+2)+3.*q(I+1))**2/4.
!          IS1=13.*(q(I+2)-2.*q(I+1)+q(I))**2/12.
!     $                     +(q(I+2)-q(I))**2/4.
!          IS2=13.*(q(I+1)-2.*q(I)+q(I-1))**2/12.
!     $                     +(3.*q(I+1)-4.*q(I)+q(I-1))**2/4.
!          tao5=abs(is0-is2)
!          tao0=abs(is0-is1)
!          tao1=abs(is2-is1)
!          dmin=dmin1(is0,is1,is2)
!          if (tao5.gt.dmin) then         
!             end_p(NTC)=i
!             NTC=NTC+1
!             start_p(NTC)=i+1
!
!             Q30=A300*q(i+3)+A301*q(i+2)+A302*q(i+1)
!             Q31=A310*q(i+2)+A311*q(i+1)+A312*q(i)
!             Q32=A320*q(i+1)+A321*q(i)+A322*q(i-1)
!       
!
!                 aa0=C30*(1.0d0+(tao5/(ss+IS0)))
!                 aa1=C31*(1.0d0+(tao5/(ss+IS1)))
!                 aa2=C32*(1.0d0+(tao5/(ss+IS2)))          
!
!                 w0=aa0/(aa0+aa1+aa2)
!                 w1=aa1/(aa0+aa1+aa2)
!                 w2=aa2/(aa0+aa1+aa2)
!
!                 qr(i)=w0*Q30+w1*Q31+w2*Q32
!          else
!           
!             Q30=A300*q(i-2)+A301*q(i-1)+A302*q(i)
!             Q31=A310*q(i-1)+A311*q(i)+A312*q(i+1)
!             Q32=A320*q(i)+A321*q(i+1)+A322*q(i+2)
!       
!                qr(i)=C30*Q30+C31*Q31+C32*Q32      
!                
!                qr(i)=(2.0d0*q(i-2)-13.0d0*q(i-1)+47.0d0*q(i)
!     $                +27.0d0*q(i+1)-3.d0*q(i+2))/60.d0     
!            
!           end if          
!       end do            ! end of (do k=1,NTC)        
!           
!
!      end select
!     
!      return
!        end
!        
!        
!      subroutine tridiagsolve(n,a,b,c,d,x)
!c     Forward sweep and back substitution algorithm for
!c     tridiagonal matrix
!        implicit none
!        integer, intent(in)::
!     $        n                  ! dimension
!        double precision, dimension(n)::
!     $        a,b,c,d            !a(i)*x(i-1)+b(i)*x(i)+c(i)*x(i+1)=d(i)
!        double precision, dimension(n),intent(out)::
!     $        x 
!c       local variable
!        integer::
!     $        i       
!c             start
!c             Modify the coefficients
!        c(1)=c(1)/b(1)
!        d(1)=d(1)/b(1)
!        do i=2,n
!           c(i)=c(i)/(b(i)-c(i-1)*a(i))
!           d(i)=(d(i)-d(i-1)*a(i))/(b(i)-c(i-1)*a(i))
!        end do
!c             Back substitute
!        x(n)=d(n)
!        do i=n-1,1,-1
!           x(i)=d(i)-c(i)*x(i+1)
!        end do
!
!        return
!
!        end
