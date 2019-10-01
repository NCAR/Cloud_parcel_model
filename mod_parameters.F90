module parameter_mod
    ! for mic block
  implicit none
      INTEGER, PARAMETER :: DP = SELECTED_REAL_KIND(14)
      real, parameter :: p0=1.0e5 !reference pressure
      real, parameter :: diffvnd = 2.55d-5              ! Coefficient of diffusion of water vapour in air [m**2/s]
      real, parameter :: ka = 2.48d-2                   ! Thermal conductivity of air [J/msK]
      real,parameter :: pi=3.1415926
      real,parameter ::  KK = 8.54d-11
      real,parameter :: grav=9.8
      real,parameter :: visc = 1.78e05!0.16d-4!1.78e-5
      real,parameter :: lat = 2.477d6!2.5e6
      real,parameter :: ra=287.0 
      real,parameter :: cp=1004.0!1005.0 
      real,parameter :: rv= 467!461.5!467
      real,parameter :: m_w=18.d-3 !molecular weight of water
      real,parameter :: m_d=29.0d-3 !molecular weight of dry air
      real,parameter :: rhow=1000.0 !density of water
      real,parameter :: eps=ra/rv
      real,parameter :: latovercp=lat/cp
      real :: m_s !molecular weight of solute; ammonium sulfate=132.14d-3; NaCl = 58.44d-3
      real :: rho_ccn,vh !van Hoff factor and density of ccn !kg/m**3 for ammonium sulfate =1726.d0 for NaCl=2160.d0
      real*8,parameter :: sigma_sa=7.61d-2
      integer :: disp,nbins,nbinsout
      integer :: idebug
      real*8  ::  temperature,ks,rhoa
      real*8  :: diffvnd1,diffvnd2,ka1,ka2
end module parameter_mod
