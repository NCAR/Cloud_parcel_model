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
      use parameter_mod
      implicit none
      real(dp)            :: qvpp,qvs,esat,temp_old,time,delt,time_prep
      real(dp)            :: sp,volume,h,thetapp,seq,seq1,seq2
      real(dp)            :: exner,racp,p1,PP,sumrp
      real(dp)            :: rm,rm0,curv,solu
      real(dp)            :: deltaqp,lwc,cql,cdp
      character*4       :: name
      real(dp), allocatable, dimension(:) :: rad,rad_ccn,dsd,dr3,rad_wet,kappa
      integer           :: iinit,ifinal, isolu
      integer           :: GCCN
      integer           :: iter,ntmic,ntot,i
      real(dp)            :: ndrop       
      real,parameter :: upp=0.0 !spin-up updraft
      real :: up !updraft 
196   format(1x,6(f16.8,2x))  
145   format(1x,100(e16.8,2x))

      time_prep=.0
      delt= 1.d-04 
      ntot=3./delt
      name = 'IUGG'
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
      idebug = 0           !=1 switch on debug mode: print out variables
      disp = 20 
      GCCN = 0 
      isolu = 1            !=1 kappa form solute term; =2 classical solute term
      up = 2.0             !  updraft velocity
      sp = -14.39d-2       !supersaturation
      temperature=284.3d0  !initial temperature
      p1=93850.0d0         !initial pressure
      nbins = 100
      allocate (dsd(nbins),rad(nbins),rad_ccn(nbins),dr3(nbins),rad_wet(nbins),kappa(nbins))
      if (disp .eq. 2) then ! mono or bi-disperse
	      rho_ccn=1726.d0!ammonium sulfate
	      m_s=132.14d-3
	      vh = 3.!2.
         kappa(1:nbinsout)=vh*m_w/m_s*rho_ccn/rhow
      elseif (disp .eq. 30 .or. disp .eq. 35 ) then !Xue10  case
	      rho_ccn=1726.d0!2160.d0!1726.d0 !ammonium sulfate
	      m_s=132.14d-3!58.44d-3!132.14d-3 
	      vh = 3.!2.
         kappa(1:nbinsout)=vh*m_w/m_s*rho_ccn/rhow
      elseif (disp .eq. 31 .or. disp .eq. 32) then !NJ17
	      rho_ccn=2160.d0 !NaCl
	      m_s=58.44d-3
	      vh = 2.
	      GCCN=1
         kappa(1:nbinsout)=vh*m_w/m_s*rho_ccn/rhow
      elseif (disp .eq. 20) then !IUGG case mono backgound + GCCN
         kappa(1) =  0.2d0
         kappa(2) =  1.2d0
         GCCN  =  2  !monodisperse GCCN
      endif
      !-------------------------------------------------------!
      rm =0.d0
      rad=rm
      dsd=0.d0
      dr3=0.d0
      call iaerosol(rad_ccn,dsd,ndrop,rm,GCCN)
      print*,'radius',rad_ccn(1:2)
      ! to get dry radius rad_ccn
      print*,'maximum binsize is',nbinsout
      if (GCCN .ne. 0) print*,'GCCN is on, value = ',GCCN
      write(50,145) 0.,(dsd(i), i=1,nbinsout)
      write(51,145) 0.,(rad(i), i=1,nbinsout)
      print*,'number of drops ',ndrop
      print*,'dispersion type ',disp
!-------------initialize variables------------
!      if(isolu .eq. 1) then
!         print*, 'kappa = ', kappa(1)
!      endif
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
   if (idebug .eq. 1) print*,'dry radius',rad_ccn(1:2)
   call wetradius(isolu,sp,dsd,rad_ccn,rad_wet,kappa)

      rad=rad_wet !droplet radius
      if (idebug .eq. 1) print*,'dry & wet radius',rad_ccn(1:2),rad(1:2)
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
            if (idebug .eq. 1) print*,'id',i,'solu=',solu, 'seq=',seq, 'supersat=',sp-seq!mark
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
             write(16,*) time,0,Sp,0.,rad(1),rad(2),pp,temperature,&
               thetapp,qvpp,qvs,rm,rhoa,LATOVERCP*DELTAQP/EXNER,deltaqp
             write(50,145) time,(dsd(i), i=1,nbinsout)
             write(51,145) time,(rad(i), i=1,nbinsout)
	  endif
  200 enddo
endif ! spin_up
lwc=sum(4.d0/3.d0*pi*rad(1:nbinsout)**3*rhow*dsd(1:nbinsout))*1e6
write(16,*) 0.d0-time_prep,0.d0,Sp,lwc,rad(1),rad(2),pp,temperature, & !8
            thetapp,qvpp,qvs,rm,rhoa,LATOVERCP*DELTAQP/EXNER,deltaqp
!--------------evolution------------------
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
            elseif (isolu .eq. 0) then !no solute effect
               solu=1.d0
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
            if (idebug .eq. 1) print*,'id',i,'rad',rad(i),'solu',solu,'seq',seq,'supersat=',sp-seq,'s_env',sp
          enddo 
            rm=(sum(rad(1:nbinsout)**3*dsd(1:nbinsout))/ndrop)**(1.d0/3.0d0)
            if(mod(ntmic,int(1./delt)) .eq. 0 .or. ntot .le. 1000) then
               lwc=sum(4.d0/3.d0*pi*rad(1:nbinsout)**3*rhow*dsd(1:nbinsout))*1e6
              write(16,*) time,up*delt*ntmic,Sp,lwc,rad(1),rad(2),pp,temperature, &
                thetapp,qvpp,qvs,rm,rhoa,LATOVERCP*DELTAQP/EXNER,deltaqp
              write(50,145) time,(dsd(i), i=1,nbinsout)
              write(51,145) time,(rad(i), i=1,nbinsout)
            endif
  100 enddo

      close(unit=16)
      deallocate (dsd, rad, rad_ccn,dr3,rad_wet)
      close(unit=50)


      end program parcel


real function log2(x)
  implicit none
  real(8), intent(in) :: x

  log2 = log(x) / log(2.d0)
end function