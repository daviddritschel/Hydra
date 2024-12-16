module parameters

! This module contains all the modifiable parameters (except nz & pi below)
! for the suite of casl f90 files.

 !The number pi (*** non-modifiable ***):
double precision,parameter:: pi=3.141592653589793238462643383279502884197169399375105820974944592307816d0

!***Numerical parameters:***
integer,parameter:: nx=N_X,ny=N_Y

 !Number of contours used for representing the buoyancy variation:
integer,parameter:: ncontq=N_CONTQ
! ncontq : used to compute the buoyancy contour interval from 
!          dq = (qq_max-qq_min)/ncontq

 !Simulation time length etc..
double precision,parameter:: tsim=T_SIM
double precision,parameter:: tgsave=T_GSAVE,tcsave=T_CSAVE
! tsim   : total duration of the simulation
! tgsave : grid data save time increment 
! tcsave : contour data save time increment (approximate)

!***Physical parameters:***
double precision,parameter:: ellx=L_X,elly=L_Y
double precision,parameter:: depth=L_Z
double precision,parameter:: cdamp=C_DAMP,nnu=POW_HYPER
double precision,parameter:: kappa=T_DIFF,dcdx=DC_DX,dcdy=DC_DY
! ellx   : domain width in x (periodic, centred at 0)
! elly   : domain width in y (periodic, centred at 0)
! depth  : N*H/f, the scaled depth
! cdamp  : this times |zeta|_rms*(k/k_max)^(2*nnu) is the hyperdiffusivity 
!          coefficient, where nnu is the power specified above
! kappa  : if > 0, allows diffusive evolution of a tracer field;
!          then kappa is the diffusivity
! dcdx   : Mean value of tracer gradient dc/dx
! dcdy   : Mean value of tracer gradient dc/dy
!----------------------------------------------------------------

end module
