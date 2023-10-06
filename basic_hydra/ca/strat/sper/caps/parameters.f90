module parameters

! This module contains all the modifiable parameters for 
! the suite of caps f90 files.

 !Domain grid dimensions:
integer,parameter:: nx=N_X,ny=N_Y

 !Domain width in x and limits in y:
double precision,parameter:: ellx=L_X
double precision,parameter:: ymin=Y_MIN,ymax=Y_MAX

 !Simulation duration and data save interval:
double precision,parameter:: tsim=T_SIM,tgsave=T_GSAVE,tcsave=T_CSAVE

 !Number of contours used for representing buoyancy and vorticity:
integer,parameter:: ncontb=N_CONTB,ncontz=N_CONTZ
! ncontb : used to compute the buoyancy contour interval from 
!          dbb = (bb_max-bb_min)/ncontb 
! ncontz : used to compute the vorticity contour interval from 
!          dzz = (zz_2/zz_1)/ncontz  where zz_n is the L_n norm
!          over those points with vorticity satisfying |zz| > zz_2
 
 !Reference translational velocity (added to u):
double precision,parameter:: uref=U_REF

 !(Hyper-)viscosity parameters (for residual vorticity only):
integer,parameter:: nnu=N_NU
double precision,parameter:: prediss=PRE_DISS
! nnu      : power of damping term (1 -> Laplacian)
! prediss  : the damping rate at wavenumber k is
!          : zz_char * prediss * (k/kmax)^(2*nnu)
!          : where kmax is the maximum x or y wavenumber
!          : nnu = 3 and prediss = 10.d0 are recommended

end module
