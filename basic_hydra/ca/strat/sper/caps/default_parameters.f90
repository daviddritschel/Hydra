module parameters

! This module contains all the modifiable parameters for 
! the suite of casl f90 files.

 !Domain grid dimensions:
integer,parameter:: nx=1024,ny=128

 !Domain width in x and limits in y:
double precision,parameter:: ellx=51200d0
double precision,parameter:: ymin=0.d0,ymax=6400d0

 !Simulation duration and data save interval:
double precision,parameter:: tsim=900.d0,tgsave=10.d0,tcsave=100.d0

 !Number of contours used for representing buoyancy and vorticity:
integer,parameter:: ncontb=100,ncontz=20
! ncontb : used to compute the buoyancy contour interval from 
!          dbb = (bb_max-bb_min)/ncontb 
! ncontz : used to compute the vorticity contour interval from 
!          dzz = (zz_2/zz_1)/ncontz  where zz_n is the L_n norm
!          over those points with vorticity satisfying |zz| > zz_2

 !Reference translational velocity (added to u):
double precision,parameter:: uref=0.d0

 !(Hyper-)viscosity parameters (for residual vorticity only):
integer,parameter:: nnu=3
double precision,parameter:: prediss=10.d0
! nnu      : power of damping term (1 -> Laplacian)
! prediss  : the damping rate at wavenumber k is
!          : zz_char * prediss * (k/kmax)^(2*nnu)
!          : where kmax is the maximum x or y wavenumber
!          : nnu = 3 and prediss = 10.d0 are recommended

end module
