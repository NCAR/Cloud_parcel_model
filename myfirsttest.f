      program myfirst
      implicit none
      real*8 :: tau,source,qvpp,qvs,PP,esat,ks,temp,dp,time
      real*8 :: rhoa,sp,sp2,vol,avgconc,cql,h,thetapp,cdp
      real*8 :: exner,racp,edr,p0,p1,rm,sumrp,delt
      real*8 :: aa11,aa22
      real*8 :: thetap,qvp,deltaqp,dtheta
      character*2 :: name
      integer :: iter,ntmic,ntot
      integer :: ndrop,dt
      real,parameter :: pi=3.1415
      real,parameter ::  KK = 8.54d-11
      real,parameter :: grav=9.8
      real,parameter :: visc = 0.16d-4
      real,parameter :: lat = 2.477d6
      real,parameter :: up=2.5
      real,parameter ::  ra=287.0
      real,parameter :: cp=1005.0
      real,parameter :: rv= 461.5
      real,parameter :: rhow=1000.0
      real,parameter :: eps=ra/rv
      real,parameter :: latovercp=lat/cp
      do 160 dt=1,10
      delt= 1.E-05*dt*5
      write(name,'(I0)') dt
      ntot=300/delt
      edr = 0.05
      OPEN(UNIT=16,FILE='test'//name//'.out',ACCESS='APPEND')
      h = PI/1.3*(visc**3./EDR)**.25
      vol = h**3
      rm = 5.0d-6
      ndrop=1280000
      avgconc=ndrop/(64**3*vol)
      temp=283.15
      p1=90000.0
      pp=p1
      p0=1.d5
      racp=ra/cp
      exner=(pp/p0)**racp
      thetapp=temp/exner
      rhoa=pp/(ra*temp)
      ks=kk
      cql=4.0d0*pi*rhow/(rhoa*vol)
      esat = 2.53d11*exp(-5.42d3/temp)
      qvs = eps*esat/(PP-esat)
      qvpp = qvs
      sp = 1.
      sp2= 2.
      do 6 iter=1,100
      if(abs((sp-sp2)/sp2) .lt. 0.01) goto 6 
c      print*,(sp-sp2)/sp2
      tau = qvs/(cql*rm*ks*(1.0d0 + qvpp*lat**2*(1.0d0 + qvs/eps)
     +                 /(Rv*Temp**2*cp)))
      source = grav*UP*(lat/(Rv*TEMP**2*Cp)
     +                -1/(Ra*TEMP))*(PP/(PP-esat))
      sp2=sp
      sp = tau*source/(vol*avgconc)
      qvpp=(sp+1.0d0)*qvs
c          write(*,*) 'iter=,',iter,'esat=',esat,'cql=',cql
c          write(*,*) 'sp=',sp,'tau=,',tau,'source=',source
c          write(*,*) 'qvpp=',qvpp,'qvs=',qvs,'ks=',ks
 6    continue      
      dp = rhoa*grav*Up*delt
      print*,rhoa*grav*Up*1.d-3,rhoa*grav*Up*1.d-4,rhoa*grav*Up*1.d-2
      print*,'dp'
      do 100 ntmic=1,ntot
         time = ntmic*delt
         PP = P1-DP*ntmic
         exner = (PP/P0)**RACP         
         sumrp = rm*ndrop
         cql=4.0d0*pi*rhow/(rhoa*vol)
         Cdp=delt*cql*KS*SUMRP/(64**3)
         do 16 iter=1,2
            temp=thetapp*exner+latovercp*Cdp*Sp
	    esat=2.53d11*exp(-5.42d3/temp)
 	    qvs = eps*esat/(PP-Esat)
            aa11=qvpp-Cdp*Sp-(1.0d0+Sp)*Qvs
            aa22 = 
     1       eps*pp/(pp-esat)**2 
     1         * esat*lat/(Rv*temp**2) 
     1         * latovercp*Cdp
	    aa22 = -CDP-(1.D0+SP)*AA22-QVS
            sp =  sp-aa11/aa22
   16     continue
          temp=thetapp*exner+latovercp*Cdp*Sp
          rhoa=pp/(Ra*temp)
          exner = ((pp+dp)/p0)**racp
  	  deltaqp = cdp*sp
          thetap=thetapp
          qvp=qvpp
          thetapp=THETAPP+LATOVERCP*DELTAQP/EXNER
          dtheta=latovercp*deltaqp/exner/delt
           QVPP    = QVPP - DELTAQP
         ! rm=sqrt(rm**2+2*ks*delt*sp) 
            if(mod(time,0.5) .eq. 0.d0) then
              print*,'time=',time,'dtheta=',dtheta
              write(16,*) time,pp,up*delt*ntmic,temp,sp,
     +               latovercp*deltaqp/exner,thetapp,dtheta
            endif
  100 enddo
      close(unit=16)
  160 enddo
      end program myfirst
