      program myfirst
      implicit none
      real*8 :: tau,source,qvpp,qvs,PP,esat,ks,temp,dp,time,temp0
      real*8 :: rhoa,rhoa0,sp,sp2,vol,avgconc,cql,h,thetapp,cdp,cd
      real*8 :: exner,racp,edr,p0,p1,rm,sumrp,delt
      real*8 :: aa11,aa22
      real*8 :: thetap,qvp,deltaqp,dtheta,gamma,gamma1
      character*1 :: name
      integer :: iter,ntmic,ntot
      integer :: ndrop,dt
      real, parameter :: diffvnd = 2.55d-5              ! Coefficient of diffusion of water vapour in air [m**2/s]
      real, parameter :: ka = 2.48d-2                   ! Thermal conductivity of air [J/msK]
      real,parameter :: pi=3.1415
      real,parameter ::  KK = 8.54d-11
      real,parameter :: grav=9.8
      real,parameter :: visc = 0.16d-4
      real,parameter :: lat = 2.477d6
      real,parameter :: up=10.
      real,parameter ::  ra=287.0
      real,parameter :: cp=1005.0
      real,parameter :: rv= 461.5
      real,parameter :: rhow=1000.0
      real,parameter :: eps=ra/rv
      real,parameter :: latovercp=lat/cp
      do 160 dt=1,1
      delt= 1.E-02*10**dt
      write(name,'(I0)') dt
      ntot=300/delt
      edr = 0.05
      OPEN(UNIT=16,FILE='test'//name//'.out',ACCESS='APPEND')
      print*,'ntot=',ntot,'delt=',delt,'file=test'//name//'.out'
      h = PI/1.3*(visc**3./EDR)**.25
      vol = h**3
      rm = 8.0d-6
      avgconc=2.d8 !ndrop/(64**3*vol)
      ndrop=avgconc*64**3*vol
      cd=rm**3*4./3.*pi*rhow*avgconc
      temp=280.15
      p1=80000.0
      pp=p1
      p0=1.d5
      racp=ra/cp
      exner=(pp/p0)**racp
      thetapp=temp/exner
      rhoa=pp/(ra*temp)
      rhoa0=rhoa
      print*,rhoa
      cql=4.0d0*pi*rhow/(rhoa*vol)
      esat = 2.53d11*exp(-5.42d3/temp)
      ks=1.0/(lat**2*eps*rhow/(Ka*Ra*Temp**2)+Ra*Temp*rhow/(eps*diffvnd*esat)) !kk
      qvs = eps*esat/(PP-esat)
      qvpp = qvs
!      sp = 0.
!      sp2= 0.2
!      do 6 iter=1,100
!      if(abs((sp-sp2)/sp2) .lt. 0.01) goto 6 
!      tau = qvs/(cql*rm*ks*(1.0d0 + qvpp*lat**2*(1.0d0 + qvs/eps)&
!                      /(Rv*Temp**2*cp)))
!      source = grav*UP*(lat/(Rv*TEMP**2*Cp) &
!                     -1/(Ra*TEMP))*(PP/(PP-esat))
!      sp2=sp
!      sp = tau*source/(vol*avgconc)
!      qvpp=(sp+1.0d0)*qvs
! 6    continue      
      sp = 0.
      qvpp=(sp+1.0d0)*qvs
      do 100 ntmic=1,ntot
         dp = rhoa*grav*Up*delt
         time = ntmic*delt
         PP = PP-DP
         exner = (PP/P0)**RACP         
         sumrp = rm*ndrop
         cql=4.0d0*pi*rhow/(rhoa*vol)
         Cdp=delt*cql*KS*SUMRP/(64**3)
!         cdp=0.d0
         temp0=temp
         do 16 iter=1,2
            temp=thetapp*exner+latovercp*Cdp*Sp
!            temp=temp0-grav/cp*delt*up+latovercp*cdp*sp
	    esat=2.53d11*exp(-5.42d3/temp)
 	    qvs = eps*esat/(PP-Esat)
            aa11=qvpp-Cdp*Sp-(1.0d0+Sp)*Qvs
            aa22 = eps*pp/(pp-esat)**2 &
              * esat*lat/(Rv*temp**2) &
              * latovercp*Cdp
	    aa22 = -CDP-(1.D0+SP)*AA22-QVS
            sp =  sp-aa11/aa22
   16     continue
          temp=thetapp*exner+latovercp*Cdp*Sp
         ks=1.0/(lat**2*eps*rhow/(Ka*Ra*Temp**2)+Ra*Temp*rhow/(eps*diffvnd*esat)) !kk
!          temp=temp0-grav/cp*delt*up+latovercp*cdp*sp
         !print*,(thetapp*exner-temp0)/delt*up
          rhoa=pp/(Ra*temp)
          exner = ((pp+dp)/p0)**racp
  	  deltaqp = cdp*sp
          cd=cd+deltaqp
          thetap=thetapp
          qvp=qvpp
          thetapp=THETAPP+LATOVERCP*DELTAQP/EXNER
          dtheta=latovercp*deltaqp/exner/delt
          QVPP    = QVPP - DELTAQP
          rm=sqrt(rm**2+2*ks*delt*sp) 
            if(mod(ntmic,ntot/1000) .eq. 0 .or. ntot .le. 1000) then
              write(16,*) time,ks,pp,up*delt*ntmic,temp,sp,&
!6 terms
                    dp,cdp,rm,thetapp,qvpp,cd
            endif
  100 enddo
      gamma=(temp-283.15)/(up*delt*ntmic)
      gamma1=-grav/cp +latovercp*cd/(up*time)
      print*,'gamma=',gamma,'gamma1=',gamma1,rhoa,cd
      close(unit=16)
  160 enddo
      end program myfirst
