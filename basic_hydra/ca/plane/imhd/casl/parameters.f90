module parameters

! This module contains all the modifiable parameters (except nz & pi below)
! for the suite of casl f90 files.

 !The number pi (*** non-modifiable ***):
double precision,parameter:: pi=3.141592653589793238462643383279502884197169399375105820974944592307816d0

!***Numerical parameters:***
integer,parameter:: nx=N_X,ny=N_Y

 !Number of contours used for representing the PV variation:
integer,parameter:: ncontq=N_CONTQ
! ncontq : used to compute the PV contour interval from 
!          dq = (qq_max-qq_min)/ncontq

 !Simulation time length etc..
double precision,parameter:: tsim=T_SIM
double precision,parameter:: tgsave=T_GSAVE,tcsave=T_CSAVE
! tsim   : total duration of the simulation
! tgsave : grid data save time increment 
! tcsave : contour data save time increment (approximate)

!***Physical parameters:***
double precision,parameter:: ellx=L_X,elly=L_Y
double precision,parameter:: kd=K_D,beta=PV_GRAD
double precision,parameter:: b0=B_0,eta=MAG_DIFF
double precision,parameter:: rtherm=R_THERM,rekman=R_EKMAN
double precision,parameter:: esr=E_SR,vorvor=VOR_VOR
integer,parameter:: ivor=I_VOR,iseed=I_SEED
! ellx   : domain width in x (periodic, centred at 0)
! elly   : domain width in y (periodic, centred at 0)
! kd     : Rossby deformation wavenumber associated with the baroclinic mode
! beta   : planetary vorticity gradient
! b0     : mean zonal magnetic field in x direction
! eta    : magnetic diffusivity
! rtherm : thermal damping rate
! rekman : Ekman   damping rate
! esr    : enstrophy input rate (via pairs of point vortices which
!          are converted to gridded vorticity and added to qd)
! vorvor : the mean vorticity associated with the point vortex 
!          upon placement on the grid
! ivor   : use 1 for arbitrarily spaced pairs (effectively monopoles)
!          use 2 for dipoles concentrated at a point, and 0 otherwise
! iseed  : seed for initialising point vortex forcing
!----------------------------------------------------------------

end module
