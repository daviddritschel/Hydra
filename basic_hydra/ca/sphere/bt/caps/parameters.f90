module parameters

! This module contains all the modifiable parameters for 
! the suite of spe f90 files.

!***Numerical parameters:***
integer,parameter:: ng=N_G, nt=2*ng
integer,parameter:: nperiod=N_PER
double precision,parameter:: dtmax=DT_MAX,tsave=T_SAV,tsim=T_SIM
double precision,parameter:: dq=PV_JUMP,cdamp=C_DAMP
! ng, nt : number of grid boxes in latitude & longitude (inversion grid)
! nperiod: total number of periods to run (each of length tsim)
! dtmax  : maximum allowed timestep
! tsave  : approximate time interval between data saves (of fields)
! tsim   : simulation duration (one "period")
! dq     : PV jump across all contours
! cdamp  : this times zz_rms is the damping rate on wavenumber ng;
!          cdamp*zz_rms*(m/(ng*cos(lat)))^4 is the damping applied to
!          wavenumber m; latitude damping at the same rate is also done

!***Physical parameters:***
double precision,parameter:: asp=ASPECT,omega=OMEGA
double precision,parameter:: rekman=R_EKMAN,esr=ESR
integer,parameter:: ksr=KSR,iseed=ISEED
! asp    : the height:width aspect ratio of the ellipsoidal surface
! omega  : the planetary rotation rate
! rekman : Ekman friction rate (the relative vorticity zeta
!          is relaxed back to zero at the rate rekman)
! esr    : Enstrophy injection rate per unit of time
! ksr    : Wavenumber centroid of enstrophy injection; 
!          the enstrophy spectrum is proportional to 
!          k^5*exp(-2k^2/k0^2) --- see ranspec
! iseed  : seed for initialising stochastic forcing
!----------------------------------------------------------------

end module
