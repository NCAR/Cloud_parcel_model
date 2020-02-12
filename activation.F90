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
        ! rad_wet, rad_ccn: wet & dry radius of the particle
        ! kappa : hygroscopicity, see petters & Kredenweis 2007
        ! solu : solute term calculated either from kappa (isolu=1) or from classical format,
        !see Jensen & Nugent 2017, eqn (1)
        ! disp: DSD dispersion flags
        ! GCCN: natural GCCN & seeding particle size flags
        ! 
        !///////////////////////////////////////////////////////////////////////////////!!
      use parameter_mod
      implicit none
      real*8            :: qvpp,qvs,esat,temp_old,time,delt,time_prep
      real*8            :: sp,volume,h,thetapp,seq,seq1,seq2
      real*8            :: exner,racp,p1,PP,sumrp
      real*8            :: rm,rm0,curv,solu
      real*8            :: deltaqp,lwc,cql,cdp
      character*50        :: name
      !real*8, allocatable, dimension(:) :: rad_wet,rad_ccn,nrad,dr3,kappa
      integer           :: GCCN,iseed,isolu
      integer           :: iter,ntmic,ntot,i,irad,imax
      real*8            :: ndrop
      real              :: time_tot
      real,parameter :: upp=0.0 !spin-up updraft
      real :: up !updraft 
196   format(1x,6(f16.8,2x))  
145   format(1x,101(e16.8,2x))
      !----------read parameters from file-------------------------!
      open(90,file='parameter.dat',status='unknown')
      read(90,*) 
      read(90,*) 
      read(90,*) name
      read(90,*)  
      read(90,*) time_tot,time_prep
      read(90,*)  
      read(90,*) disp,GCCN,iseed
      read(90,*)  
      read(90,*) up
      read(90,*)  
      read(90,*) sp
      read(90,*)  
      read(90,*) isolu
      read(90,*)  
      read(90,*) temperature
      read(90,*)  
      read(90,*) p1
      read(90,*)  
      read(90,*) idebug
      close(unit=90)
      !!----------------------------------------------------------!!
      !------------------------------setup output files-------------------------------!!
      OPEN(UNIT=16,FILE=trim(name)//'.out',ACCESS='APPEND')!parcel mean variables
      OPEN(UNIT=50,FILE=trim(name)//'.nrad',ACCESS='APPEND')!number concentration of each bin
      OPEN(UNIT=51,FILE=trim(name)//'.rad',ACCESS='APPEND')!droplet size of each bin
      OPEN(UNIT=60,FILE=trim(name)//'.info',ACCESS='APPEND')!output the tested variables
      OPEN(UNIT=52,FILE=trim(name)//'.dndr',ACCESS='APPEND')!output the remapped dsd (dr=1micron)
      !--------------------------------------initialize aerosols----------------------!!
      !------------------------------initial size distribution------------------------!!
      !////GCCN=1 insert Jensen & Nugent 2017 Giant CCN           ////////////////////!!
      !//////   =2 monotonic seeding r=1micron;                   ////////////////////!!
      !//////=3 add 3 mode lognormal seeding distribution by Cooper et al. 1997///////!!
      !///disp=35 Xue10 urban, 30 Xue10 marine, 31 JN17 polluted, 32 NJ17 pristine////!!
      !////disp=1 monotonic initial sizes ////////////////////////////////////////////!!
      !-------------------------------------------------------------------------------!!
!      time_prep=.0
      delt= 1.d-04 
      ntot=time_tot/delt
      nbins = 100
      nbins2 = 50*iseed !seeded particle bins
      nbinsout2 =0 !actual seeded bins initialize
      allocate(nrad(nbins+nbins2),rad_ccn(nbins+nbins2),dr3(nbins+nbins2))
      allocate(rad_wet(nbins+nbins2),kappa(nbins+nbins2))
      allocate(wid(nbins+nbins2),dNdlogr(nbins+nbins2),dNdr(nbins+nbins2))
      rm =0.d0
      rad_ccn=rm
      rad_wet=1.5*rad_ccn
      kappa=0.d0
      nrad=0.d0
      dr3=0.d0
      ndrop=0.d0
      ! to get dry radius rad_ccn DSD & kappa
      call iaerosol(ndrop,rm,GCCN,iseed)
      print*,'maximum binsize is',nbinsout+nbinsout2
      !-------------------------------------------------------!
      write(50,145) 0.,(nrad(i), i=1,nbinsout+nbinsout2)
      write(51,145) 0.,(rad_ccn(i), i=1,nbinsout+nbinsout2)
      print*,'number of drops ',ndrop
      print*,'dispersion type ',disp
!-------------initialize variables------------
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
      write(60,*) '-------------initial conditions--------------'
      write(60,*) 'initial temperature = ',temperature
      write(60,*) 'initial RH = ',sp+1.d0
      write(60,*) 'initial pressure=',p1
      write(60,*) 'updraft velocity=',up
      write(60,*) 'debug mode flag=', idebug
      write(60,*) 'seeding = ',iseed,'GCCN=',GCCN
      write(60,*) 'kappa=',kappa(1:nbinsout+nbinsout2)
      write(60,*) '----------------------------------------------'
!--------------first guess the radius of wet aerosol--------!
!--------------start with r_wet=1.5*r_d-----------------------!
   if (idebug .eq. 1) print*,'dry radius',rad_ccn(1:2)
   call wetradius(isolu,sp)!,nrad,rad_ccn,rad_wet,kappa)
   if (idebug .eq. 1) then
      print*,'dry & wet radius',rad_ccn(1:2),rad_wet(1:2)
      print*,'nbinsout=',nbinsout,'nbinsout2=',nbinsout2
      print*,'nbins=',nbins,'nbins2=',nbins2
      print*, 'get wet radius','radwet',rad_wet
   endif
!--------------spin-up------------------
if(time_prep .ne. 0) then
      do 200 ntmic=1,int(time_prep/delt*2.d0)
         time = ntmic*delt/2.d0-time_prep
         pp=rhoa*Ra*(1.d0+18.d0/29.d0*qvpp)*temperature
         exner = (PP/P0)**RACP
         sumrp=sum(dr3(1:nbinsout)*nrad(1:nbinsout))/3.d0
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
          do i = 1,nbinsout+nbinsout2
             curv=2.d0*sigma_sa/(Rv*rhow*temperature*rad_wet(i)) !curvature effect coefficient !unit in m
            if(isolu .eq. 1) then !kappa
               solu=(rad_wet(i)**3-rad_ccn(i)**3)/(rad_wet(i)**3-(1-kappa(i))*rad_ccn(i)**3)
            elseif (isolu .eq. 2) then !classical solute term
               solu=exp(-vh*m_w/m_s*rho_ccn*rad_ccn(i)**3/rhow/(rad_wet(i)**3-rad_ccn(i)**3))
            endif !isolu
            
            diffvnd2=diffvnd1*1.d0/(rad_wet(i)/(rad_wet(i)+0.104d-6)+diffvnd1/(rad_wet(i)*0.036)*sqrt(2.d0*pi/(Ra*temperature)))     
            ka2=ka1*1.d0/(rad_wet(i)/(rad_wet(i)+.216d-6)+ka1/(rad_wet(i)*.7*rhoa*cp)*sqrt(2.d0*pi/(Ra*temperature)))
            ks =1.d0/(rhow*Rv*temperature/(esat*diffvnd2)+rhow*Lat/(Ka2*temperature)*(Lat/(Rv*temperature)-1))
            !!caculate equillibrium supersat.
             seq = solu * exp(curv)-1.d0
            !if (idebug .eq. 1) print*,'id',i,'solu=',solu, 'seq=',seq, 'supersat=',sp-seq!mark
             rm0=rad_wet(i)*3.0d0*delt*ks*(sp-seq)+rad_wet(i)**3 !!r^3 scheme
             if(rm0 .gt. rad_ccn(i)**3 ) then
                dr3(i)=rm0-rad_wet(i)**3
                rad_wet(i)=rm0**(1.d0/3.d0)
             else
                dr3(i)=0.d0
                rad_wet(i)=rad_ccn(i)
             endif
          enddo
          rm=(sum(rad_wet(1:nbinsout+nbinsout2)**3*nrad(1:nbinsout+nbinsout2))/ndrop)**(1.d0/3.0d0)
          if(mod(ntmic,int(1.d0/delt)) .eq. 0 .or. int(time_prep/delt*2.d0) .le. 1000) then
             write(16,*) time,0,Sp,0.,ndrop,pp,temperature,&
               thetapp,qvpp,qvs,rm,rhoa,LATOVERCP*DELTAQP/EXNER,deltaqp
             write(50,145) time,(nrad(i), i=1,nbinsout+nbinsout2)
             write(51,145) time,(rad_wet(i), i=1,nbinsout+nbinsout2)
	  endif
  200 enddo
endif ! spin_up
lwc=sum(4.d0/3.d0*pi*rad_wet(1:nbinsout+nbinsout2)**3*rhow*nrad(1:nbinsout+nbinsout2))*1e6
write(16,*) 0.d0-time_prep,0.d0,Sp,lwc,ndrop,pp,temperature, & !8
            thetapp,qvpp,qvs,rm,rhoa,LATOVERCP*DELTAQP/EXNER,deltaqp
!--------------evolution------------------
print*,'evolution of droplet size starting with rm=',rm
      do 100 ntmic=1,ntot
         time = ntmic*delt
         pp=rhoa*Ra*(1.d0+m_w/m_d*qvpp)*temperature!pp-rhoa*grav*delt*up!
         exner = (PP/P0)**RACP         
         sumrp=sum(dr3(1:nbinsout+nbinsout2)*nrad(1:nbinsout+nbinsout2))/3.d0
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
          do i = 1,nbinsout+nbinsout2
             curv=2.d0*sigma_sa/(Rv*rhow*temperature*rad_wet(i)) !curvature effect coefficient !in m
            if(isolu .eq. 1) then !kappa
               solu=(rad_wet(i)**3-rad_ccn(i)**3)/(rad_wet(i)**3-(1-kappa(i))*rad_ccn(i)**3)
            elseif (isolu .eq. 2) then !classical solute term
               solu=exp(-vh*m_w/m_s*rho_ccn*rad_ccn(i)**3/rhow/(rad_wet(i)**3-rad_ccn(i)**3)) !solute effect coefficient ms=132.14 for ammonium sulfate !m^3
            elseif (isolu .eq. 0) then !no solute effect
               solu=1.d0
            endif !isolu
            diffvnd2=diffvnd1*1.d0/(rad_wet(i)/(rad_wet(i)+0.104d-6)+diffvnd1/(rad_wet(i)*0.036)*sqrt(2.d0*pi/(Ra*temperature)))
            ka2=ka1*1.d0/(rad_wet(i)/(rad_wet(i)+.216d-6)+ka1/(rad_wet(i)*.7*rhoa*cp)*sqrt(2.d0*pi/(Ra*temperature)))
            ks =1.d0/(rhow*Rv*temperature/(esat*diffvnd2)+rhow*Lat/(Ka2*temperature)*(Lat/(Rv*temperature)-1))
!!--------------caculate equillibrium supersat.
            seq=solu*exp(curv)-1.d0
            rm0=rad_wet(i)*3.0d0*delt*ks*(sp-seq)+rad_wet(i)**3 !!r^3 
	         if(rm0 .gt. rad_ccn(i)**3) then
	            dr3(i)=rm0-rad_wet(i)**3
	            rad_wet(i)=rm0**(1.d0/3.d0)
	         else
		         dr3(i)=0.d0
               rad_wet(i)=rad_ccn(i)
	         endif
            !if (idebug .eq. 1) print*,'id',i,'rad_wet',rad_wet(i),'solu',solu,'seq',seq,'supersat=',sp-seq,'s_env',sp
          enddo 
            rm=(sum(rad_wet(1:nbinsout+nbinsout2)**3*nrad(1:nbinsout+nbinsout2))/ndrop)**(1.d0/3.0d0)
            if(mod(ntmic,int(1./delt)) .eq. 0 .or. ntot .le. 1000) then
               lwc=sum(4.d0/3.d0*pi*rad_wet(1:nbinsout+nbinsout2)**3*rhow*nrad(1:nbinsout+nbinsout2))*1e6
               dNdr=0.d0! initialize dNdr
               imax=1
               ndrop=0
               do i = 1,nbinsout+nbinsout2
                  irad=ceiling(rad_wet(i)/1.d-6) !map to the bin index in dNdr
                  imax=max(imax,irad)
                  if (imax .gt. nbins) then
                     print*,'index out of the maximum dimension of bins'
                     imax=nbins
                     irad=imax
                  endif
                  dNdr(irad)=dNdr(irad)+nrad(i)
               enddo
               ndrop=sum(dNdr(1:imax))
               !rm=0.d0
               !do i=2,imax
               !   rm=rm+(real(i)*1.d-6)**3 * dndr(i)
               !enddo
               !rm=(rm/sum(dndr(2:imax)))**(1.d0/3.d0)
               !ndrop=sum(nrad(1:nbinsout+nbinsout2))
               write(52,145) time, (dNdr(i),i=1,nbins)
               write(16,*) time,up*delt*ntmic,Sp,lwc,ndrop,pp,temperature, &
                  thetapp,qvpp,qvs,rm,rhoa,LATOVERCP*DELTAQP/EXNER,deltaqp
               write(50,145) time,(nrad(i), i=1,nbinsout+nbinsout2)
               write(51,145) time,(rad_wet(i), i=1,nbinsout+nbinsout2)
                          !map rad_wet into dNdr (dr=1micron)               
            endif

  100 enddo

      close(unit=16)
      deallocate (nrad, rad_wet, rad_ccn, dr3, kappa)
      deallocate(wid,dNdlogr,dNdr)
      close(unit=50)


      end program parcel


real function log2(x)
  implicit none
  real(8), intent(in) :: x

  log2 = log(x) / log(2.d0)
end function
