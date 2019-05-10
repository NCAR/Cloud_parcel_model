program Parcel
        !aerosol activation for parcel of 1cm^3 in volumeume
        !///////////////////////////////////////////////////////////////////////////////!!
        ! list of variables 
        ! qvpp, qvs are the water vapor mixing ratio & its repective saturated value
        ! esat water vapor pressure
        ! temperature, temp_old, thetapp are new & old temperature & potential temperature
        ! delt: time step 
        ! time_prep: spin up time
        ! rhoa, rhow, rho_ccn, density of air, water & solute
        ! sp, seq are supersaturation wrt water & equilibrium value over solute particle
        ! h & volume are width & volume of parcel box
        ! lwc: liquid water content
        ! dsd: droplet size distribution with nbinsout width
        ! ndrop: total droplet number
        ! rad, rad_ccn: wet & dry radius of the particle
        ! kappa : hygroscopicity, see petters & Kredenweis 2007
        ! solu : solute term calculated either from kappa (isolu=1) or from classical format,
        !see Jensen & Nugent 2017, eqn (1)
        ! disp: DSD dispersion flags
        ! 
        !///////////////////////////////////////////////////////////////////////////////!!
      implicit none
      real*8            :: qvpp,qvs,esat,ks,temperature,temp_old,time,delt,time_prep
      real*8            :: rhoa,sp,volume,h,thetapp,seq,seq1,seq2
      real*8            :: exner,racp,p1,PP,sumrp
      real*8            :: rm,rm0,curv,solu
      real*8            :: deltaqp,lwc,cql,cdp
      character*6       :: name
      real(8), allocatable, dimension(:) :: rad,rad_ccn,dsd,dr3,rad_wet,kappa
      integer           :: nbins,nbinsout,iinit,ifinal, isolu
      integer           :: disp,GCCN
      integer           :: iter,ntmic,ntot,i
      real*8            :: ndrop
      real*8            :: diffvnd1,diffvnd2,ka1,ka2
      real, parameter :: p0=1.0e5 !reference pressure
      real, parameter :: diffvnd = 2.55d-5              ! Coefficient of diffusion of water vapour in air [m**2/s]
      real, parameter :: ka = 2.48d-2                   ! Thermal conductivity of air [J/msK]
      real,parameter :: pi=3.1415
      real,parameter ::  KK = 8.54d-11
      real,parameter :: grav=9.8
      real,parameter :: visc = 1.78e05!0.16d-4!1.78e-5
      real,parameter :: lat = 2.5e-6!2.477d6!2.5e-6
      real,parameter :: ra=287.0 
      real,parameter :: cp=1004.0!1005.0 
      real,parameter :: rv= 467!461.5!467
      real,parameter :: m_w=18.d-3 !molecular weight of water
      real,parameter :: rhow=1000.0 !density of water
      real,parameter :: eps=ra/rv
      real,parameter :: latovercp=lat/cp
      real*8,parameter :: sigma_sa=7.61d-2
      real,parameter :: upp=0.0 !spin-up updraft
      real :: up,vh !updraft velo and van Hoff factor 
      real :: m_s !molecular weight of solute; ammonium sulfate=132.14d-3; NaCl = 58.44d-3
      real :: rho_ccn!kg/m**3 for ammonium sulfate =1726.d0 for NaCl=2160.d0
196   format(1x,6(f16.8,2x))
145   format(1x,100(e16.8,2x))

      time_prep=.0
      delt= 1.d-04 
      ntot=300./delt
      name = 'simple'
      !------------------------------setup output files-------------------------------!!
      OPEN(UNIT=16,FILE=name//'.out',ACCESS='APPEND')!parcel mean variables
      OPEN(UNIT=50,FILE=name//'.dsd',ACCESS='APPEND')!number concentration of each bin
      OPEN(UNIT=51,FILE=name//'.rad',ACCESS='APPEND')!droplet size of each bin
      OPEN(UNIT=60,FILE=name//'.test',ACCESS='APPEND')!output the tested variables
      !--------------------------------------initialize aerosols----------------------!!
      !------------------------------initial size distribution------------------------!!
      !////GCCN=1 insert Jensen & Nugent 2017 Giant CCN           ////////////////////!!
      !//////   =2 monotonic seeding r=1micron;                   ////////////////////!!
      !//////=3 add 3 mode lognormal seeding distribution by Cooper et al. 1997///////!!
      !///disp=35 Xue10 urban, 30 Xue10 marine, 31 JN17 polluted, 32 NJ17 pristine////!!
      !////disp=1 monotonic initial sizes ////////////////////////////////////////////!!
      !-------------------------------------------------------------------------------!!
      disp = 1
      GCCN = 0 
      isolu = 1            !=1 kappa form solute term; =2 classical solute term
      up = 2.0             !  updraft velocity
      sp = -14.39d-2       !supersaturation
      temperature=284.3d0  !initial temperature
      p1=93850.0d0         !initial pressure
      nbins = 100
      allocate (dsd(nbins),rad(nbins),rad_ccn(nbins),dr3(nbins),rad_wet(nbins),kappa(nbins))
      if (disp .le. 2) then !bi-disperse
	      rho_ccn=1726.d0!ammonium sulfate
	      m_s=132.14d-3
	      vh = 3.!2.
      elseif (disp .eq. 30 .or. disp .eq. 35 ) then !Xue10  case
	      rho_ccn=1726.d0!2160.d0!1726.d0 !ammonium sulfate
	      m_s=132.14d-3!58.44d-3!132.14d-3 
	      vh = 3.!2.
      elseif (disp .eq. 31 .or. disp .eq. 32) then !NJ17
	      rho_ccn=2160.d0 !NaCl
	      m_s=58.44d-3
	      vh = 2.
	      GCCN=1
      endif
      !-------------------------------------------------------!
      rm =0.d0
      rad=rm
      dsd=0.d0
      dr3=0.d0
      call iaerosol(disp,rad_ccn,dsd,nbins,ndrop,rho_ccn,rm,nbinsout,GCCN)
      ! to get dry radius rad_ccn
      print*,'maximum binsize is',nbinsout
      if (GCCN .ne. 0) print*,'GCCN is on, value = ',GCCN
      write(50,145) 0.,(dsd(i), i=1,nbinsout)
      write(51,145) 0.,(rad(i), i=1,nbinsout)
      print*,'number of drops ',ndrop
      print*,'dispersion type ',disp
!-------------initialize variables------------
      if(isolu .eq. 1) then
         kappa(1:nbinsout)=vh*m_w/m_s*rho_ccn/rhow
         print*, 'kappa = ', kappa(1)
      endif
      h = 1.d-2 !.01m=1cm
      volume = h**3
      pp=p1
      sumrp=0.d0
      racp= ra/cp
      exner=(pp/p0)**racp
      thetapp=temperature/exner
      cql=4.0d0*pi*rhow/volume
      esat = 2.53d11*exp(-5.42d3/temperature)
      ks = 1.d0/(rhow*Rv*temperature/(esat*diffvnd)+rhow*Lat/(Ka*temperature)*(Lat/(Rv*temperature)-1))
      qvs = eps*esat/(PP-esat)
      qvpp= (sp+1.d0)*qvs !kg/m^3
      rhoa=pp/(Ra*(1+18.d0/29.d0*qvpp)*temperature)
!--------------first guess the radius of wet aerosol--------!
!--------------start with r_wet=1.5*r_d-----------------------!
      rad_wet=rad_ccn*1.5d0
      diffvnd1=1.d-5*(0.015*temperature-1.9)
      ka1=1.5d-11*temperature**3-4.8d-8*temperature**2+1.d-4*temperature-3.9d-4
      do i=1,nbinsout
         seq=sp+.01d0
	      seq1=seq+.01d0
	      seq2=seq+.01d0
      do while(abs(seq-sp) .gt. 1.d-7 .and. seq2 .ne. seq) 
	      seq2=seq1
	      seq1=seq
         diffvnd2=diffvnd1*1.d0/(rad_wet(i)/(rad_wet(i)+0.104d-6)+diffvnd1/(rad_wet(i)*0.036)*sqrt(2.d0*pi/(Ra*temperature)))
         ka2=ka1*1.d0/(rad_wet(i)/(rad_wet(i)+.216d-6)+ka1/(rad_wet(i)*.7*rhoa*cp)*sqrt(2.d0*pi/(Ra*temperature)))
         ks = 1.d0/(rhow*Rv*temperature/(esat*diffvnd2)+rhow*Lat/(Ka2*temperature)*(Lat/(Rv*temperature)-1))
         curv=2.d0*sigma_sa/(Rv*rhow*temperature*rad_wet(i))  

         if(isolu .eq. 1) then !kappa
            solu=(rad_wet(i)**3-rad_ccn(i)**3)/(rad_wet(i)**3-(1-kappa(i))*rad_ccn(i)**3)
         elseif (isolu .eq. 2) then !classical solute term
            solu=exp(-vh*m_w/m_s*rho_ccn/rhow*rad_ccn(i)**3/(rad_wet(i)**3-rad_ccn(i)**3)) !solute effect coefficient ms=132.14 for ammonium sulfate !m^3
         endif !isolu

         seq = solu * exp(curv)-1.d0           
         if(seq .gt. sp .and. rad_wet(i) .gt. rad_ccn(i)) then
	         rad_wet(i)=rad_wet(i)-rad_ccn(i)*1.d-6
	      elseif(seq .lt. sp) then 
	         rad_wet(i)=rad_wet(i)+rad_ccn(i)*1.d-6
	      endif
      enddo
      enddo
      rad=rad_wet !droplet radius
!--------------spin-up------------------
if(time_prep .ne. 0) then
      do 200 ntmic=1,int(time_prep/delt*2.d0)
         time = ntmic*delt/2.d0-time_prep
         pp=rhoa*Ra*(1.d0+18.d0/29.d0*qvpp)*temperature
         exner = (PP/P0)**RACP
         sumrp=sum(dr3(1:nbinsout)*dsd(1:nbinsout))/3.d0
         deltaqp=cql*sumrp!new cdp condensation
         temp_old=temperature
         temperature=temp_old-grav/cp*delt/2.0d0*upp+latovercp*deltaqp!mark new
         esat=2.53d11*exp(-5.42d3/temperature)
         qvs = eps*esat/(PP-esat)
         rhoa= rhoa*(-grav*upp/(Ra*temperature)*delt/2.d0-(temperature-temp_old)/temperature)+rhoa
          thetapp=thetapp+latovercp*deltaqp/exner
          qvpp    = qvpp - deltaqp
          sp = qvpp/qvs-1.d0 !mark new sp
          rm=0.d0
          diffvnd1=1.d-5*(0.015*temperature-1.9)
          ka1=1.5d-11*temperature**3-4.8d-8*temperature**2+1.d-4*temperature-3.9d-4
          do i = 1,nbinsout
             curv=2.d0*sigma_sa/(Rv*rhow*temperature*rad(i)) !curvature effect coefficient !unit in m
            if(isolu .eq. 1) then !kappa
               solu=(rad(i)**3-rad_ccn(i)**3)/(rad(i)**3-(1-kappa(i))*rad_ccn(i)**3)
            elseif (isolu .eq. 2) then !classical solute term
               solu=exp(-vh*m_w/m_s*rho_ccn*rad_ccn(i)**3/rhow/(rad(i)**3-rad_ccn(i)**3))
            endif !isolu
            
            diffvnd2=diffvnd1*1.d0/(rad(i)/(rad(i)+0.104d-6)+diffvnd1/(rad(i)*0.036)*sqrt(2.d0*pi/(Ra*temperature)))     
            ka2=ka1*1.d0/(rad(i)/(rad(i)+.216d-6)+ka1/(rad(i)*.7*rhoa*cp)*sqrt(2.d0*pi/(Ra*temperature)))
            ks =1.d0/(rhow*Rv*temperature/(esat*diffvnd2)+rhow*Lat/(Ka2*temperature)*(Lat/(Rv*temperature)-1))
            !!caculate equillibrium supersat.
             seq = solu * exp(curv)-1.d0
            print*,'solu=',solu, 'seq=',seq, 'supersat=',sp-seq!mark
             rm0=rad(i)*3.0d0*delt*ks*(sp-seq)+rad(i)**3 !!r^3 scheme
             if(rm0 .gt. rad_ccn(i)**3 ) then
                dr3(i)=rm0-rad(i)**3
                rad(i)=rm0**(1.d0/3.d0)
             else
                dr3(i)=0.d0
                rad(i)=rad_ccn(i)
             endif
          enddo
          rm=(sum(rad(1:nbinsout)**3*dsd(1:nbinsout))/ndrop)**(1.d0/3.0d0)
          if(mod(ntmic,int(1.d0/delt)) .eq. 0 .or. int(time_prep/delt*2.d0) .le. 1000) then
             write(16,*) time,0,Sp,0.,rad(15),rad(40),pp,temperature,281.5-grav/cp*delt*up*ntmic,&
               thetapp,qvpp,qvs,rm,rhoa,LATOVERCP*DELTAQP/EXNER,deltaqp
             write(50,145) time,(dsd(i), i=1,nbinsout)
             write(51,145) time,(rad(i), i=1,nbinsout)
	  endif
  200 enddo
endif ! spin_up
lwc=sum(4.d0/3.d0*pi*rad(1:nbinsout)**3*rhow*dsd(1:nbinsout))
write(16,*) 0.d0-time_prep,up*delt*ntmic-287.6,Sp,wc,rad(15),rad(40),pp,temperature,281.5-grav/cp*delt*up,& !8
            thetapp,qvpp,qvs,rm,rhoa,LATOVERCP*DELTAQP/EXNER,deltaqp
!--------------evolumeution------------------
      do 100 ntmic=1,ntot
         time = ntmic*delt
         pp=rhoa*Ra*(1+m_w/29.d-3*qvpp)*temperature
         exner = (PP/P0)**RACP         
         sumrp=sum(dr3(1:nbinsout)*dsd(1:nbinsout))/3.d0
	      deltaqp=cql*sumrp!cdp condensation
         temp_old=temperature
         temperature=temp_old-grav/cp*delt*up+latovercp*deltaqp!mark new
         esat=2.53d11*exp(-5.42d3/temperature)
         qvs = eps*esat/(PP-esat)
         rhoa= rhoa*(-grav*up/(Ra*temperature)*delt-(temperature-temp_old)/temperature)+rhoa
         thetapp=thetapp+latovercp*deltaqp/exner
         qvpp    = qvpp - deltaqp
	      sp = qvpp/qvs-1.d0 !mark new sp
         rm=0.d0
         diffvnd1=1.d-5*(0.015*temperature-1.9)
         ka1=1.5d-11*temperature**3-4.8d-8*temperature**2+1.d-4*temperature-3.9d-4
          do i = 1,nbinsout
             curv=2.d0*sigma_sa/(Rv*rhow*temperature*rad(i)) !curvature effect coefficient !in m
            if(isolu .eq. 1) then !kappa
               solu=(rad(i)**3-rad_ccn(i)**3)/(rad(i)**3-(1-kappa(i))*rad_ccn(i)**3)
            elseif (isolu .eq. 2) then !classical solute term
               solu=exp(-vh*m_w/m_s*rho_ccn*rad_ccn(i)**3/rhow/(rad(i)**3-rad_ccn(i)**3)) !solute effect coefficient ms=132.14 for ammonium sulfate !m^3
            endif !isolu
            diffvnd2=diffvnd1*1.d0/(rad(i)/(rad(i)+0.104d-6)+diffvnd1/(rad(i)*0.036)*sqrt(2.d0*pi/(Ra*temperature)))
            ka2=ka1*1.d0/(rad(i)/(rad(i)+.216d-6)+ka1/(rad(i)*.7*rhoa*cp)*sqrt(2.d0*pi/(Ra*temperature)))
            ks =1.d0/(rhow*Rv*temperature/(esat*diffvnd2)+rhow*Lat/(Ka2*temperature)*(Lat/(Rv*temperature)-1))
!!--------------caculate equillibrium supersat.
            seq=solu*exp(curv)-1.d0
            rm0=rad(i)*3.0d0*delt*ks*(sp-seq)+rad(i)**3 !!r^3 
	         if(rm0 .gt. rad_ccn(i)**3) then
	            dr3(i)=rm0-rad(i)**3
	            rad(i)=rm0**(1.d0/3.d0)
	         else
		         dr3(i)=0.d0
               rad(i)=rad_ccn(i)
	         endif
          enddo 
            rm=(sum(rad(1:nbinsout)**3*dsd(1:nbinsout))/ndrop)**(1.d0/3.0d0)
            if(mod(ntmic,int(1./delt)) .eq. 0 .or. ntot .le. 1000) then
               lwc=sum(4.d0/3.d0*pi*rad(1:nbinsout)**3*rhow*dsd(1:nbinsout))
              write(16,*) time,up*delt*ntmic-287.6,Sp,lwc,rad(15),rad(40),pp,temperature,281.5-grav/cp*delt*up*ntmic,&
                thetapp,qvpp,qvs,rm,rhoa,LATOVERCP*DELTAQP/EXNER,deltaqp
              write(50,145) time,(dsd(i), i=1,nbinsout)
              write(51,145) time,(rad(i), i=1,nbinsout)
            endif
  100 enddo

      close(unit=16)
      deallocate (dsd, rad, rad_ccn,dr3,rad_wet)
      close(unit=50)


      end program parcel

    SUBROUTINE IAEROSOL(disp,rad,nrad,nbins,ndrop,rho_ccn,rm,nbinsout,GCCN)
!---- This subroutine determines the initial position and size of all droplets
  implicit none

  ! --- argument ---
  integer :: nbins,nbinsout,GCCN
  real(8), dimension(nbins) :: rad,nrad
  integer :: disp
  real(8) :: rm,ndrop

  ! --- local --

  integer :: i,iinit,ifinal
  real(8), allocatable, dimension(:) :: wid,dNdlogr,dNdr
  real(8) :: r1,n1,logsig,rmin,rmax
  real :: rho_ccn
  real(8) :: logrmin,logrmax,rad_power,bin_factor
  real,parameter :: pi=3.14159265
 145   format(1x,100(e16.8,2x))
 111 format(a20,i3)

! Set everything to zero.
rad = 0.0
ndrop=0.0 
rm=0.d0
nrad=0.
allocate(wid(nbins),dNdlogr(nbins),dNdr(nbins))
dNdlogr = 0.0
dNdr = 0.0d0
wid =0.d0
!size dispersion
  if (disp .eq. 1) then !mono disperse
      nbinsout=1
      rad(nbinsout)=1.d-7
      ndrop=100
      nrad=ndrop
  elseif (disp .eq. 30) then !lulin 2010 maritime case
     nbinsout=39
     rmin = 6.d-9
     rad(1)=rmin
     bin_factor=2.d0
     wid(1)=rad(1)*(bin_factor**(1.d0/3.d0)-1.d0)
     do i=2,nbinsout
        rad_power=real(i)/3.d0
        rad(i)=rad(1)*bin_factor**rad_power
        wid(i)=rad(i)-rad(i-1)
     enddo
     do i=1,nbinsout
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
	      if (GCCN .eq. 3) then !add 3 mode lognormal seeding distribution by Cooper et al. 1997
	      !mode 1
	         n1=100.d0
	         r1=.15d-6
	         logsig=.2d0
	         logsig=log(10.d0)*logsig
            dNdlogr(i)= dNdlogr(i)+n1/(sqrt(2.0d0*pi) *logsig) * exp(-((log(rad(i))-log(r1))/(sqrt(2.0d0)*logsig))**2)
	      !mode 2
            n1=100.d0*1.7d-4
            r1=.5d-6
            logsig=.4d0
            logsig=log(10.d0)*logsig
            dNdlogr(i)= dNdlogr(i)+n1/(sqrt(2.0d0*pi) *logsig) * exp(-((log(rad(i))-log(r1))/(sqrt(2.0d0)*logsig))**2)
	      !mode 3
            n1=100.d0*3.d-7
            r1=5.d-6
            logsig=.6d0
            logsig=log(10.d0)*logsig
            dNdlogr(i)= dNdlogr(i)+n1/(sqrt(2.0d0*pi) *logsig) * exp(-((log(rad(i))-log(r1))/(sqrt(2.0d0)*logsig))**2)
	      endif!GCCN=3
     enddo
     do i=1,nbinsout 
        dNdr(i)=dNdlogr(i)/rad(i)
        nrad(i)=dNdr(i)*wid(i)
     enddo
  elseif (disp .eq. 31) then !Jensen&Nugent 2017 modified polluted case
     nbinsout=30
     rmin = 1.d-8
     rad(1)=rmin
     bin_factor=2.0d0 !mass increment
     wid(1)=rad(1)*(bin_factor**(1.d0/3.d0)-1.d0)
     do i=2,nbinsout
        rad_power=real(i)/3.d0
        rad(i)=rad(1)*bin_factor**rad_power
        wid(i)=rad(i)-rad(i-1)
     enddo
     do i=1,nbinsout
        n1=48.d0!160!48.d0
        r1=.029d-6
        logsig=1.36d0
        logsig=log(logsig)
        dNdlogr(i)= dNdlogr(i)+n1/(sqrt(2.0d0*pi) *logsig) *exp(-((log(rad(i))-log(r1))/(sqrt(2.0d0)*logsig))**2)
        n1=125.d0!380!125.d0
        r1=.071d-6
        logsig=1.57d0
        logsig=log(logsig)
        dNdlogr(i)= dNdlogr(i)+n1/(sqrt(2.0d0*pi) *logsig)*exp(-((log(rad(i))-log(r1))/(sqrt(2.0d0)*logsig))**2)
     enddo
     dNdr(1:nbinsout)=dNdlogr(1:nbinsout)/rad(1:nbinsout)
     nrad(1:nbinsout)=dNdr(1:nbinsout)*wid(1:nbinsout)


  elseif (disp .eq. 32) then !Jensen&Nugent 2017 pristine case
     nbinsout=39
     rmin = 1.d-8
     rad(1)=rmin
     bin_factor=2.0d0 !mass increment
     wid(1)=rad(1)*(bin_factor**(1.d0/3.d0)-1.d0)
     do i=2,nbinsout
        rad_power=real(i)/3.d0
        rad(i)=rad(1)*bin_factor**rad_power
        wid(i)=rad(i)-rad(i-1)
     enddo
     do i=1,nbinsout
        n1=125.d0
        r1=.011d-6
        logsig=1.2d0
        logsig=log(logsig)
        dNdlogr(i)= dNdlogr(i)+n1/(sqrt(2.0d0*pi) *logsig)*exp(-((log(rad(i))-log(r1))/(sqrt(2.0d0)*logsig))**2)
        n1=65.d0
        r1=.06d-6
        logsig=1.7d0
        logsig=log(logsig)
        dNdlogr(i)= dNdlogr(i)+n1/(sqrt(2.0d0*pi) *logsig)*exp(-((log(rad(i))-log(r1))/(sqrt(2.0d0)*logsig))**2)
     enddo
        dNdr(1:nbinsout) = dNdlogr(1:nbinsout)/rad(1:nbinsout)
        nrad(1:nbinsout) = dNdr(1:nbinsout)*wid(1:nbinsout)
  elseif (disp .eq. 35) then !Lulin 2010 rural
     rmin = 6.d-9
     rad(1)=rmin
     bin_factor=2.d0 !mass increment
     wid(1)=rad(1)*(bin_factor**(1.d0/3.0d0)-1)
     do i=2,35
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
     do i=1,35
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
     enddo
  endif ! disp
  if(GCCN .eq. 1) then !JN17 GCCN
     do i=1,42
        rad(nbinsout+i)=.8d-6+0.2d-6*real(i-1)
     enddo
        nrad(nbinsout+1)=0.1118d0
        nrad(nbinsout+2)=.06849d0!1micron
        nrad(nbinsout+3)=.0384d0
        nrad(nbinsout+4)=.02182d0
        nrad(nbinsout+5)=.0133d0
        nrad(nbinsout+6)=.8496d-2
        nrad(nbinsout+7)=.5486d-2!2micron
        nrad(nbinsout+8)=.3805d-2
        nrad(nbinsout+9)=.2593d-2
        nrad(nbinsout+10)=.1919d-2
        nrad(nbinsout+11)=.1278d-2
        nrad(nbinsout+12)=.9884d-3!3micron
        nrad(nbinsout+13)=.7779d-3
        nrad(nbinsout+14)=.5195d-3
        nrad(nbinsout+15)=.4005d-3
        nrad(nbinsout+16)=.3769d-3
        nrad(nbinsout+17)=.2653d-3!4micron
        nrad(nbinsout+18)=.2124d-3
        nrad(nbinsout+19)=.1378d-3
        nrad(nbinsout+20)=.1214d-3
        nrad(nbinsout+21)=.1009d-3
        nrad(nbinsout+22)=.1222d-3!5micron
        nrad(nbinsout+23)=.5064d-4
        nrad(nbinsout+24)=.383d-4
        nrad(nbinsout+25)=.5547d-4
        nrad(nbinsout+26)=.2145d-4
        nrad(nbinsout+27)=.1295d-4!6micron
        nrad(nbinsout+28)=.4323d-4
        nrad(nbinsout+29)=.2626d-4
        nrad(nbinsout+30)=.305d-4
        nrad(nbinsout+31)=.4385d-5
        nrad(nbinsout+32)=.4372d-5!7micron
        nrad(nbinsout+33)=.4465d-5
        nrad(nbinsout+34)=.4395d-5
        nrad(nbinsout+35)=.4427d-5
        nrad(nbinsout+36)=.4411d-5
        nrad(nbinsout+37)=.0d0    !8micron
        nrad(nbinsout+38)=.0d0
        nrad(nbinsout+39)=.0d0
        nrad(nbinsout+40)=.4522d-5
        nrad(nbinsout+41)=.0d0
        nrad(nbinsout+42)=.4542d-5!9micron
        nbinsout=nbinsout+42
  elseif (GCCN==2) then !some simple one size GCCN
        rad(nbinsout+1)=1.d-6
        nrad(nbinsout+1)=10
        nbinsout=nbinsout+1
  endif!GCCN
  ndrop=sum(nrad(1:nbinsout)) !total number
  write(50,145) 0.,(dNdlogr(i), i=1,nbinsout) !initial dry size distribution
  rm = (sum(rad(1:nbinsout)**3*nrad(1:nbinsout))/ndrop)**(1.d0/3.d0)
  print*,"rm", rm, 'nbins',nbins,'nbinsout',nbinsout,'ndrop',ndrop
deallocate(wid,dNdlogr,dNdr)

 299   format(1x,(e12.6),3(f8.6))
 199   format(1x,3(f8.6))


  end SUBROUTINE IAEROSOL

real function log2(x)
  implicit none
  real(8), intent(in) :: x

  log2 = log(x) / log(2.d0)
end function
