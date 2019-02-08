      program myfirst
      implicit none
      real*8 :: tau,source,qvpp,qvs,PP,esat,ks,temp,dp,time,temp0
      real*8 :: rhoa,rhoa0,sp,sp2,vol,avgconc,cql,h,thetapp,cdp,cd
      real*8 :: exner,racp,p0,p1,sumrp,edr
      real*8   :: rm,rm1,rm0,curv,solu
      real*8 :: aa11,aa22,vtemp,vtemp1,delt
      real*8 :: thetap,qvp,deltaqp,dtheta,gamma,gamma1
      character*1 :: name
      real(8), allocatable, dimension(:) :: rad,rad1,dsd,mass
      integer :: nbins,iinit,ifinal
      integer :: disp
      integer :: iter,ntmic,ntot,i
      real(8) :: ndrop
      integer, parameter :: NN =64
      real, parameter :: diffvnd = 2.55d-5              ! Coefficient of diffusion of water vapour in air [m**2/s]
      real, parameter :: ka = 2.48d-2                   ! Thermal conductivity of air [J/msK]
      real,parameter :: pi=3.1415
      real,parameter ::  KK = 8.54d-11
      real,parameter :: grav=9.8
      real,parameter :: visc = 0.16d-4
      real,parameter :: lat = 2.477d6
      real,parameter :: up=5.d0
      real,parameter ::  ra=287.0
      real,parameter :: cp=1005.0
      real,parameter :: rv= 461.5
      real,parameter :: rhow=1000.0
      real,parameter :: eps=ra/rv
      real,parameter :: latovercp=lat/cp
      real,parameter :: rho_ccn = 1726.d0 !kg/m**3
196   format(1x,6(f16.8,2x))
145   format(1x,51(e16.8,2x))

      
      delt= 1.d-01
      ntot=300./delt
      print*,ntot
      OPEN(UNIT=16,FILE='test.out',ACCESS='APPEND')
      OPEN(UNIT=50,FILE='test.dsd',ACCESS='APPEND')
      OPEN(UNIT=51,FILE='test.rad',ACCESS='APPEND')
      !initialize aerosols
      disp = 30 !35 urban, 30 marine
      nbins = 40
      allocate (dsd(nbins),rad(nbins),rad1(nbins),mass(nbins))
      rm =1.d-5
      rad=rm
      dsd=5.d0
!      ndrop=2.d2 !cm^-3
     call iaerosol(disp,rad,mass,dsd,nbins,ndrop,iinit,ifinal,rm)
      write(50,145) 0.,(dsd(i), i=1,nbins)
      write(51,145) 0.,(rad(i), i=1,nbins)
     print*,'ndrop',ndrop,'disp',disp
      rad1=rad
      h = 1.d-2 !.01m=1cm
      vol = h**3
      avgconc= ndrop/vol!ndrop/(64**3*vol) !m^-3
      Cd=rm**3*4./3.*pi*rhow*avgconc
      temp=293.15
      p1=95000.0
      pp=p1
      p0=1.d5
      racp= ra/cp
      exner=(pp/p0)**racp
      thetapp=temp/exner
      rhoa=pp/(ra*temp)
      ndrop=avgconc*vol
      rhoa0=rhoa
      cql=4.0d0*pi*rhow/(rhoa*vol)
      esat = 2.53d11*exp(-5.42d3/temp)
!      ks=1.0/(lat**2*eps*rhow/(Ka*Ra*Temp**2)+Ra*Temp*rhow/(eps*diffvnd*esat)) !kk
      ks=kk
      qvs = eps*esat/(PP-esat)
      qvpp = qvs
      sp = 0.d0
      qvpp=qvs
      do 100 ntmic=1,ntot
         time = ntmic*delt
         dp = rhoa*grav*Up*delt
         PP = PP-DP
         exner = (PP/P0)**RACP         
         sumrp = rm*ndrop
         cql=4.0d0*pi*rhow/(rhoa*vol)
         Cdp=delt*cql*KS*SUMRP
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
           !print*,'sp',sp,'a11,aa22',aa11,aa22,'temp',temp,'esat',esat,'qvs',qvs
   16     continue
          !sp=0.003
          deltaqp=cdp*sp
          temp=thetapp*exner+latovercp*deltaqp
!         ks=1.0/(lat**2*eps*rhow/(Ka*Ra*Temp**2)+Ra*Temp*rhow/(eps*diffvnd*esat)) !kk
          ks=kk
!          temp=temp0-grav/cp*delt*up+latovercp*cdp*sp
          rhoa=pp/(Ra*temp)
          avgconc=ndrop/(vol**3*rhoa)
!          exner = ((pp+dp)/p0)**racp
!  	  deltaqp = rm**3*4./3.*pi*rhow*avgconc-cd!cdp*sp
!          cd=cd+deltaqp
!          cd=rm**3*4./3.*pi*rhow*avgconc
          thetap=thetapp
          qvp=qvpp
          thetapp=THETAPP+LATOVERCP*DELTAQP/EXNER
!          dtheta=latovercp*deltaqp/exner/delt
         !print*,'cd',cd,'deltaqp',deltaqp,'dtheta',dtheta,'thetapp',thetapp
          QVPP    = QVPP - DELTAQP
!          vtemp=2.0d0*ks*delt
!         vtemp1=3.0d0*ks*delt
          rm=0.d0

!          rm=sqrt(rm**2+2*ks*delt*sp)
          do i = 1,nbins
             curv=3.3d-7/temp !curvature effect coefficient !in m
             solu=4.3d0*2.d0*mass(i)/132.14d0*1.d-6 !solute effect coefficient ms=132.14 for ammonium sulfate !m^3
            ! print*,i,'rad(i)**2',rad(i)**2,vtemp*sp+2.d0*delt*ks*(-curv/rad(i)+solu/rad(i)**3)
	      rm0=3.0d0*delt*ks*(sp*rad(i)-curv)+rad(i)**3!*fv(i)
!             rm0=rad(i)**2+vtemp*sp +2.d0*delt*ks*(-curv/rad(i)+solu/rad(i)**3)
!             rm0=3.d0*delt*ks*(sp*rad(i)-curv+solu/rad(i)**2)+rad(i)**3
!             if (rm0 .gt. 0) rad(i)=sqrt(rm0)
             if(rm0 .gt. 0) rad(i)=(rm0)**(1.d0/3.d0)
!             if(rad(i) .lt. rad1(i)) rad(i)=rad1(i)
             rm=rm+rad(i)*dsd(i)

          enddo 
             rm=rm/ndrop
            if(mod(ntmic,ntot/1000) .eq. 0 .or. ntot .le. 1000) then
              write(16,*) time,up*delt*ntmic,Sp,rm,maxval(rad),rad(8),pp,temp,thetapp,qvpp,sumrp,rhoa
              write(50,145) time,(dsd(i), i=1,nbins)
              write(51,145) time,(rad(i), i=1,nbins)
            endif
  100 enddo
!      gamma=(temp-283.15)/(up*delt*ntmic)
!      gamma1=-grav/cp +latovercp*cd/(up*time)
      ! print*,'gamma=',gamma,'gamma1=',gamma1,rhoa,cd

      close(unit=16)
      deallocate (dsd, rad, rad1,mass)
      close(unit=50)


      end program myfirst

    SUBROUTINE IAEROSOL(disp,rad,mass,nrad,nbins,ndrop,iinit,ifinal,rm)
! This subroutine determines the initial position and size of all droplets
!! for droplets locations & ID# & random # generator
  implicit none

  ! --- argument ---
  integer :: nbins
  real(8), dimension(nbins) :: rad,nrad,mass
  integer :: disp,iinit,ifinal
  real(8) :: rm,ndrop

  ! --- local --

  integer :: i
  real(8), allocatable, dimension(:) :: wid,dNdlogr,dNdr
  real(8) :: r1,n1,logsig,rmin,rmax,rho_ccn
  real(8) :: logrmin,logrmax,rad_power,bin_factor
  real,parameter :: pi=3.14159265

 111 format(a20,i3)

! Set everything to zero.
  rad = 0.0
  ndrop=0.0 
  rm=0.d0
allocate(wid(nbins),dNdlogr(nbins),dNdr(nbins))
dNdlogr = 0.0
dNdr = 0.0d0
wid =0.d0
mass=0.d0
nrad=0.
!size dispersion
  if (disp .eq. 30) then !lulin maritime case
     rmin = 6.d-9
     rad(1)=rmin
     rho_ccn=1726 !kg/m**3 for Ammonium sulfate
     mass(1)=4.d0/3.d0*pi*rmin**3*rho_ccn 
     bin_factor=2.d0
     wid(1)=rad(1)*(bin_factor**(1.d0/3.d0)-1)
     do i=2,disp
        mass(i)=mass(1)*bin_factor**i
        rad_power=real(i)/3.d0
        rad(i)=rad(1)*bin_factor**rad_power
        wid(i)=rad(i)-rad(i-1)
        n1=133.d0
        r1=0.0039d-6
        logsig=.657d0
        logsig=log(10.d0)*logsig
        dNdlogr(i) = n1/(sqrt(2.0d0*pi) *logsig) * exp(-((log(rad(i))-log(r1))/(sqrt(2.0d0)*logsig))**2)
        n1=66.6d0
        r1=.133d-6
        logsig=.21d0
        logsig=log(10.d0)*logsig
        dNdlogr(i)= dNdlogr(i)+n1/(sqrt(2.0d0*pi) *logsig) * exp(-((log(rad(i))-log(r1))/(sqrt(2.0d0)*logsig))**2)
        n1=3.06d0
        r1=.29d-6
        logsig=.396d0
        logsig=log(10.d0)*logsig
        dNdlogr(i)= dNdlogr(i)+n1/(sqrt(2.0d0*pi) *logsig) * exp(-((log(rad(i))-log(r1))/(sqrt(2.0d0)*logsig))**2)
     enddo
     do i=1,disp
        rm=rm+rad(i)**3*dNdlogr(i)*wid(i)/rad(i)
        ndrop=ndrop+dNdlogr(i)*wid(i)/rad(i)      
        dNdr(i)=dNdlogr(i)/rad(i)
        nrad(i)=dNdr(i)*wid(i)
     enddo
  elseif (disp .eq. 35) then !Lulin rural
     rmin = 6.d-9
     rad(1)=rmin
     rho_ccn=1726.d0
     mass(1)=4.d0/3.d0*pi*rmin**3*rho_ccn
     bin_factor=2.d0 !mass increment
     wid(1)=rad(1)*(bin_factor**(1.d0/3.0d0)-1)
     do i=2,disp
        mass(1)=mass(1)*bin_factor**i
        rad_power=real(i)/3.0
        rad(i)=rad(1)*bin_factor**rad_power
        wid(i)=rad(i)-rad(i-1)
        n1=6650.d0
        r1= 0.00739d-6
        logsig= .225d0
        logsig=log(10.d0)*logsig
        dNdlogr(i) = n1/(sqrt(2.0d0*pi) *logsig) * exp(-((log(rad(i))-log(r1))/(sqrt(2.0d0)*logsig))**2)
        n1= 147.d0    
        r1= .0269d-6
        logsig= .557
        logsig=log(10.d0)*logsig
        dNdlogr(i) = dNdlogr(i)+n1/(sqrt(2.0d0*pi) *logsig) * exp(-((log(rad(i))-log(r1))/(sqrt(2.0d0)*logsig))**2)
        n1= 1990.d0
        r1=.0419d-6
        logsig= .266
        logsig=log(10.d0)*logsig
        dNdlogr(i) = dNdlogr(i)+n1/(sqrt(2.0d0*pi) *logsig) * exp(-((log(rad(i))-log(r1))/(sqrt(2.0d0)*logsig))**2)
     enddo
     do i=1,disp
        rm=rm+rad(i)**3*dNdlogr(i)*wid(i)/rad(i)
        ndrop=ndrop+dNdlogr(i)*wid(i)/rad(i)
        dNdr(i)=dNdlogr(i)/rad(i)
        nrad(i)=dNdr(i)*wid(i)
     enddo
 
  elseif(disp .eq. 39 ) then !Jaenicke1988
     rmin = 6.d-9
     rmax = 5.d-6
     logrmin=10.d0**floor(log10(rmin))
     logrmax=10.d0**floor(log10(rmax))
     iinit=Nint(rmin/logrmin)
     ifinal=(floor(log10(rmax))-floor(log10(rmin)))*9+floor(rmax/logrmax)
     r1=0.
     n1=0.
     rm=0.
     logsig=0.
     do i=iinit,ifinal
          rad(i) = real(mod(i,9))*10.d0**(-7+i/9)
          if (mod(i,9) .eq. 0) then
             rad(i) = 9.0d0*10.d0**(-7+i/9-1)
          endif
          n1=133.d0
          r1=0.0039d-6
          logsig=.657d0
          dNdlogr(i) = n1/(sqrt(2.0d0*pi) *logsig) * exp(-((log10(rad(i))-log10(r1))/(sqrt(2.0d0)*logsig))**2)
          n1=66.6d0
          r1=.133d-6
          logsig=.21d0
          dNdlogr(i)= dNdlogr(i)+n1/(sqrt(2.0d0*pi) *logsig) * exp(-((log10(rad(i))-log10(r1))/(sqrt(2.0d0)*logsig))**2)
          n1=3.06d0
          r1=.29d-6
          logsig=.396d0
          dNdlogr(i)= dNdlogr(i)+n1/(sqrt(2.0d0*pi) *logsig) * exp(-((log10(rad(i))-log10(r1))/(sqrt(2.0d0)*logsig))**2)
          rm=rm+rad(i)**3*dNdlogr(i)*10.0d0**(floor(log10(rad(i))))/rad(i)/log(10.d0)
          ndrop=ndrop+dNdlogr(i)*10.0d0**(floor(log10(rad(i))))/rad(i)/log(10.d0)
     enddo

  elseif (disp .eq. 32) then !Jaenicke1988 maritime case !6d-7~5d-4cm
     !r<1d-5
     rmin = 6.d-7
     rmax = 5.d-4
     logrmin=10.d0**floor(log10(rmin))
     logrmax=10.d0**floor(log10(rmax))
     iinit=Nint(rmin/logrmin)
     ifinal=(floor(log10(rmax))-floor(log10(rmin)))*9+floor(rmax/logrmax)
     r1=0.
     n1=0.
     rm=0.
     logsig=0.
     do i=iinit,ifinal
          rad(i) = real(mod(i,9))*10.d0**(-7+i/9)
          if (mod(i,9) .eq. 0) then
             rad(i) = 9.0d0*10.d0**(-7+i/9-1)
          endif
          if(rad(i) .lt. 1.d-5) then
             n1 = 1.33d2
             r1 = 3.9d-7
             logsig = .657d0
             dNdlogr(i) = n1/(sqrt(2.0d0*pi) *logsig) * exp(-((log10(rad(i))-log10(r1))/(sqrt(2.0d0)*logsig))**2)
          elseif(rad(i) .lt. 1.d-4) then
             n1 = 6.66d1
             r1 = 1.33d-5
             logsig = .21d0
             dNdlogr(i) = n1/(sqrt(2.0*pi) *logsig) * exp(-((log10(rad(i))-log10(r1))/(sqrt(2.0d0)*logsig))**2)
             r1 = 2.9d-5
             logsig = .396d0
             dNdlogr(i) = dNdlogr(i) +n1/(sqrt(2.0*pi) *logsig) * exp(-((log10(rad(i))-log10(r1))/(sqrt(2.0d0)*logsig))**2)
          else !r>1.d-4
             n1=2.84d-1
             r1=2.0d-4
             logsig = .5d0
             dNdlogr(i) = n1/(sqrt(2.0*pi) *logsig) * exp(-((log10(rad(i))-log10(r1))/(sqrt(2.0d0)*logsig))**2)
          endif
          rm=rm+rad(i)**3*dNdlogr(i)*10.0d0**(floor(log10(rad(i))))/rad(i)/log(10.d0)
          ndrop=ndrop+dNdlogr(i)*10.0d0**(floor(log10(rad(i))))/rad(i)/log(10.d0)

     enddo

  endif ! disp
     rm = (rm/ndrop)**(1.d0/3.d0)

     print*,"rm", rm,ndrop,10.d0**(floor(log10(rad(7))))
deallocate(wid,dNdlogr,dNdr)

 299   format(1x,(e12.6),3(f8.6))
 199   format(1x,3(f8.6))


  end SUBROUTINE IAEROSOL

real function log2(x)
  implicit none
  real(8), intent(in) :: x

  log2 = log(x) / log(2.d0)
end function
