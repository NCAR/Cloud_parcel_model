    SUBROUTINE IAEROSOL(rad,nrad,ndrop,rm,GCCN)
!---- This subroutine determines the initial position and size of all droplets
     use parameter_mod
     implicit none
  ! --- argument ---
     integer :: GCCN
     real*8, dimension(nbins) :: rad,nrad
     real*8 :: rm,ndrop,ccn_ndrop,ccn_rad,gccn_ndrop,gccn_rad
  ! --- local --
     integer :: i,iinit,ifinal
     real*8, allocatable, dimension(:) :: wid,dNdlogr,dNdr
     real*8 :: r1,n1,logsig,rmin,rmax
     real*8 :: logrmin,logrmax,rad_power,bin_factor
 145 format(1x,100(e16.8,2x))
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
open(80,file='parameter.dat',status='unknown')
do i = 1,20
   read(80,*)
enddo
read(80,*) ccn_ndrop,ccn_rad
read(80,*)
read(80,*) gccn_ndrop,gccn_rad
!size dispersion
  if (disp .eq. 1) then !mono disperse
      nbinsout=1
      rad(nbinsout)=1.d-7
      ndrop=11.2!100 ! number concentration percc
      nrad=ndrop
  elseif (disp .eq. 20) then !IUGG monodisperse case
     nbinsout = 1

     nrad(nbinsout)=ccn_ndrop
     rad(nbinsout)=ccn_rad!number per cc

  elseif (disp .eq. 30) then !lulin 2010 maritime case
     nbinsout=50
     rmin = 1.d-9
     rad(1)=rmin
     bin_factor=2.d0
     wid(1)=rad(1)*(bin_factor**(1.d0/3.d0)-1.d0)
     do i=2,nbinsout
        rad_power=real(i-1)/3.d0
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
        rad_power=real(i-1)/3.d0
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
        rad_power=real(i-1)/3.d0
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
  elseif (disp .eq. 34) then !modified from Seinfeld and Pandis, 2006. p.343, urban polluted
     nbinsout=39
     rmin = 1.d-8
     rad(1)=rmin
     bin_factor=2.0d0 !mass increment
     wid(1)=rad(1)*(bin_factor**(1.d0/3.d0)-1.d0)
     do i=2,nbinsout
        rad_power=real(i-1)/3.d0
        rad(i)=rad(1)*bin_factor**rad_power
        wid(i)=rad(i)-rad(i-1)
     enddo
      do i=1,nbinsout
        n1= 7100.d0/6.d0
        r1= 0.00585d-6
        logsig= .232d0
        logsig=log(10.d0)*logsig
        dNdlogr(i) = n1/(sqrt(2.0d0*pi) *logsig) * exp(-((log(rad(i))-log(r1))/(sqrt(2.0d0)*logsig))**2)
        n1= 6320.d0/6.d0  
        r1= .0187d-6
        logsig= .250d0
        logsig=log(10.d0)*logsig
        dNdlogr(i) = dNdlogr(i)+n1/(sqrt(2.0d0*pi) *logsig) * exp(-((log(rad(i))-log(r1))/(sqrt(2.0d0)*logsig))**2)
        n1= 960.d0/6.d0
        r1=.0755d-6
        logsig= .204
        logsig=log(10.d0)*logsig
        dNdlogr(i) = dNdlogr(i)+n1/(sqrt(2.0d0*pi) *logsig) * exp(-((log(rad(i))-log(r1))/(sqrt(2.0d0)*logsig))**2)
      enddo
     do i=1,nbinsout
        dNdr(i)=dNdlogr(i)/rad(i)
        nrad(i)=dNdr(i)*wid(i)
     enddo
 
  elseif (disp .eq. 35) then !Lulin 2010 rural
     rmin = 6.d-9
     nbinsout=35
     rad(1)=rmin
     bin_factor=2.d0 !mass increment
     wid(1)=rad(1)*(bin_factor**(1.d0/3.0d0)-1)
      do i=2,nbinsout
        rad_power=real(i-1)/3.0
        rad(i)=rad(1)*bin_factor**rad_power
        wid(i)=rad(i)-rad(i-1)
      enddo
      do i=1,nbinsout
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

     nrad(nbinsout+1)=gccn_ndrop
     rad(nbinsout+1)=gccn_rad!number per cc
     nbinsout=nbinsout+1
  elseif (GCCN .eq. 3) then !add 3 mode lognormal seeding distribution by Cooper et al. 1997
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

      dNdr(1:nbinsout)=dNdr(1:nbinsout)+dNdlogr(1:nbinsout)/rad(1:nbinsout)
      nrad(1:nbinsout)=nrad(1:nbinsout)+dNdr(1:nbinsout)*wid(1:nbinsout)
  endif!GCCN
  ndrop=sum(nrad(1:nbinsout)) !total number 
  rm = (sum(rad(1:nbinsout)**3*nrad(1:nbinsout))/ndrop)**(1.d0/3.d0)
  if (idebug .eq. 1) print*,"rm", rm, 'nbins',nbins,'nbinsout',nbinsout,'ndrop',ndrop
deallocate(wid,dNdlogr,dNdr)
close(unit=80)
 299   format(1x,(e12.6),3(f8.6))
 199   format(1x,3(f8.6))


     end SUBROUTINE IAEROSOL
