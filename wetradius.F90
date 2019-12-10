SUBROUTINE wetradius(isolu,sp0,nrad,rad_ccn,rad_wet,kappa)
   ! this subroutine calculate the equilibrium radius at current thermodynamic environment
   use parameter_mod
   implicit none
   integer :: id,isolu
   real(dp) :: esat,seq, seq1,seq2
   real(dp) :: curv, solu, sp_c,sp0
   real, parameter :: sp_init=-1.d-2

   real(dp), dimension(nbins) :: nrad,rad_ccn,rad_wet,kappa
!/////////calculate necessary parameters& coefficients
   diffvnd1=1.d-5*(0.015*temperature-1.9)
   ka1=1.5d-11*temperature**3-4.8d-8*temperature**2+1.d-4*temperature-3.9d-4
   esat = 2.53d11*exp(-5.42d3/temperature)
   !--------------first guess the radius of wet aerosol--------!
   !--------------start with r_wet=1.5*r_d-----------------------!
   rad_wet=rad_ccn*1.5d0
   !--------------iterate until get equilibrium saturation status around each droplet surface----!
   sp_c=min(sp0,sp_init)
   print*, 'initial supersaturation for computing equilibrium wet radius=',sp_c
   do 100 id=1,nbinsout
      seq=sp_c+.01d0
	   seq1=seq+.01d0
	   seq2=seq+.01d0
      if (rad_ccn(id) .le. 5.d-6) then
      do 200 while(abs(seq-sp_c) .gt. 1.d-7 .and. seq2 .ne. seq) 
         seq2=seq1
         seq1=seq
         diffvnd2=diffvnd1*1.d0/(rad_wet(id)/(rad_wet(id)+0.104d-6)+diffvnd1/(rad_wet(id)*0.036)*sqrt(2.d0*pi/(Ra*temperature)))
         ka2=ka1*1.d0/(rad_wet(id)/(rad_wet(id)+.216d-6)+ka1/(rad_wet(id)*.7*rhoa*cp)*sqrt(2.d0*pi/(Ra*temperature)))
         ks = 1.d0/(rhow*Rv*temperature/(esat*diffvnd2)+rhow*Lat/(Ka2*temperature)*(Lat/(Rv*temperature)-1))
         curv=2.d0*sigma_sa/(Rv*rhow*temperature*rad_wet(id)) 
         if (isolu .eq. 1) then !kappa
            solu=(rad_wet(id)**3-rad_ccn(id)**3)/(rad_wet(id)**3-(1-kappa(id))*rad_ccn(id)**3)
         elseif (isolu .eq. 2) then !classical solute term
            solu=exp(-vh*m_w/m_s*rho_ccn/rhow*rad_ccn(id)**3/(rad_wet(id)**3-rad_ccn(id)**3)) !solute effect coefficient ms=132.14 for ammonium sulfate !m^3
         endif !isolu
         seq = solu * exp(curv)-1.d0           
         if (seq .gt. sp_c .and. rad_wet(id) .gt. rad_ccn(id)) then
	         rad_wet(id)=rad_wet(id)-rad_ccn(id)*1.d-6
	      elseif(seq .lt. sp_c) then 
	         rad_wet(id)=rad_wet(id)+rad_ccn(id)*1.d-6
	      endif
200   enddo !dowhile
      else !big r
      rad_wet(id)=1.5*rad_ccn(id)
      endif

100   enddo !drop id
end SUBROUTINE wetradius
