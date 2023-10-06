module parameters

! This module contains all the modifiable parameters (except pi below)
! for the suite of casl f90 files.

 !Initial time loop from which to begin the simulation:
integer,parameter:: loopinit=0
 !loopinit > 0 may be used for continuing a simulation (see replace below)

 !Logical to indicate replacing existing data with the output of a new
 !simulation or to append to existing data:
logical,parameter:: replace=.false.

! Number of grid boxes in the x & y directions (inversion grid):
integer,parameter:: nx=N_X
integer,parameter:: ny=N_Y

 !Number of contours used for representing PV:
integer,parameter:: ncontq=N_CONTQ
! ncontq : used to compute the PV contour interval from 
!          dq = (qq_max-qq_min)/ncontq

! Simulation time length etc..
double precision,parameter:: tsim=T_SIM
double precision,parameter:: tgsave=T_GSAVE
double precision,parameter:: tcsave=T_CSAVE
! tsim   : total duration of the simulation
! tgsave : grid data save time increment 
! tcsave : contour data save time increment (approximate)
! ***NOTE*** tcsave should always be an integer multiple of tgsave

 !The number pi (*** non-modifiable ***):
double precision,parameter:: pi=3.141592653589793238462643383279502884197169399375105820974944592307816d0

!***Physical parameters:***
double precision,parameter:: kd=K_D
double precision,parameter:: ellx=L_X
double precision,parameter:: elly=L_Y
double precision,parameter:: beta=PV_GRAD
double precision,parameter:: ubar=U_BAR
! kd     : Rossby deformation wavenumber
! ellx   : domain width in x (periodic,  centred at 0)
! elly   : domain width in y (free slip, centred at 0)
! beta   : planetary vorticity gradient
! ubar   : a uniform x velocity added to that found after inversion
!----------------------------------------------------------------

end module
