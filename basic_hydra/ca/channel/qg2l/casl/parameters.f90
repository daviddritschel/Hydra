module parameters

! This module contains all the modifiable parameters (except nz & pi below)
! for the suite of casl f90 files.

 !Initial time loop from which to begin the simulation:
integer,parameter:: loopinit=0
 !loopinit > 0 may be used for continuing a simulation (see replace below)

 !Logical to indicate replacing existing data with the output of a new
 !simulation or to append to existing data:
logical,parameter:: replace=.false.

 !Number of fluid layers (here always 2):
integer,parameter:: nz=2

! Number of grid boxes in the x & y directions (inversion grid):
integer,parameter:: nx=N_X
integer,parameter:: ny=N_Y

 !Number of contours used for representing PV in either layer:
integer,parameter:: ncontq=N_CONTQ
! ncontq : used to compute the PV contour interval from 
!          dq = (qq_max-qq_min)/ncontq (in each layer)

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
double precision,parameter:: alpha=DENSITY_RATIO
double precision,parameter:: h1=H_1
double precision,parameter:: kdbar=K_DBAR
double precision,parameter:: ellx=L_X
double precision,parameter:: elly=L_Y
double precision,parameter:: ymin=Y_M
double precision,parameter:: beta=PV_GRAD
double precision,parameter:: utopo=U_TOPO
double precision,parameter:: ftopo=F_TOPO
double precision,parameter:: epsilon=FPV_GRAD1
double precision,parameter:: rtherm1=R_THERM1
double precision,parameter:: rtherm2=R_THERM2
double precision,parameter:: rekman=R_EKMAN
double precision,parameter:: eirate=EI_RATE
double precision,parameter:: vorvor=VOR_VOR
double precision,parameter:: rheton=R_HETON
integer,parameter:: iseed=I_SEED
logical,parameter:: topogr=TOP_FLAG
! h1     : fractional thickness of the lower layer
! kdbar  : f*sqrt{(H1+H2)/(g*H1*H2*(1-alpha))} where H1 & H2 are the
!          layer depths, f is the Coriolis frequency and g is gravity
! ellx   : domain width in x (periodic, centred at 0)
! elly   : domain width in y
! ymin   : minimum value of y (ymax = ymin + elly; see constants.f90)
! beta   : planetary vorticity gradient
! utopo  : a uniform x velocity added e.g. when topography is present
! ftopo  : oscillation frequency of utopo; if non-zero, the mean barotropic
!          zonal velocity is utopo*sin(ftopo*t)
! epsilon: the mean PV gradient is eps*beta in the lower layer initially
!          (this results in a uniform velocity in each layer, and this is
!          used in the inversion to keep the boundary velocities fixed)
! rtherm1: thermal damping rate (1/tau_1) of lower layer thickness
! rtherm2: thermal damping rate (1/tau_2) of upper layer thickness
! rekman : Ekman   damping rate (1/tau_E)
! eirate : enstrophy input rate (via heton pairs of point vortices 
!          which are converted to gridded vorticity and added to qd)
! vorvor : the mean vorticity magnitude associated with the point 
!          vortex upon placement on the grid
! rheton : the radius of each of the pair of vortices added as a heton
! iseed  : seed for initialising point vortex forcing
! topogr : logical variable to indicate presence of topography
!----------------------------------------------------------------

end module
