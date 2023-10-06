module parameters

! This module contains all the modifiable parameters for 
! the suite of casl f90 files.

 !Domain grid dimensions:
integer,parameter:: nx=N_X,ny=N_Y

 !Domain limits:
double precision,parameter:: xmin=X_MIN,xmax=X_MAX
double precision,parameter:: ymin=Y_MIN,ymax=Y_MAX

 !Inclination of the bottom of the box relative to x-axis:
double precision,parameter:: themindeg=TH_MIN,themaxdeg=TH_MAX
double precision,parameter:: theinideg=TH_INI,omegag=OMEGA_G
! themin/max : Min/max tilt angles of the box (in degrees)
! theini     : Initial tilt angle of the box (in degrees)
! omegag     : Tilt oscillation angular frequency (in rad/s)

 !Particle settling velocity:
double precision,parameter:: vs=V_S

 !Number of contours used for representing buoyancy and vorticity:
integer,parameter:: ncontb=N_CONTB,ncontz=N_CONTZ
! ncontb : used to compute the buoyancy contour interval from 
!          dbb = (bb_max-bb_min)/ncontb 
! ncontz : used to compute the vorticity contour interval from 
!          dzz = (zz_2/zz_1)/ncontz  where zz_n is the L_n norm

 !Simulation time length etc..
double precision,parameter:: tsim=T_SIM
double precision,parameter:: tgsave=T_GSAVE,tcsave=T_CSAVE
! tsim   : total duration of the simulation
! tgsave : grid data save time increment 
! tcsave : contour data save time increment (approximate)

 !(Hyper-)Viscosity parameters (for residual vorticity only):
integer,parameter:: nnu=N_NU
double precision,parameter:: prediss=PRE_DISS
! nnu      : power of damping term (1 -> Laplacian)
! prediss  : if nnu=1 --> nu = prediss * (bbmax-bbmin)/kmax^(3/2) 
!          : where kmax is the max wavenumber magnitude 
!          : if nnu>1 --> nu = prediss/kmax^(2*nnu)   
!          : (Note for nnu=3, prediss=80 is recommended)

end module
