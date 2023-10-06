module parameters

! This module contains all the modifiable parameters for 
! the suite of casl f90 files.

 !Domain grid dimensions:
integer,parameter:: nx=N_X,ny=N_Y

 !Domain width in x (physical and conformal), and limits in y:
double precision,parameter:: ellxp=L_X,ellx=LT_X
double precision,parameter:: ymin=Y_MIN,ymax=Y_MAX

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
 
 !Reference translational velocity (added to u):
double precision,parameter:: uref=U_REF

 !Tide parameters:
double precision,parameter:: utidemax=UT_MAX,ttide=T_TIDE

 !(Hyper-)Viscosity parameters (for residual vorticity only):
integer,parameter:: nnu=N_NU
double precision,parameter:: prediss=PRE_DISS
! nnu      : power of damping term (1 -> Laplacian)
! prediss  : if nnu=1 --> nu = prediss * (bbmax-bbmin)/kmax^(3/2) 
!          : where kmax is the max wavenumber magnitude 
!          : if nnu>1 --> nu = prediss/kmax^(2*nnu)   
!          : (Note for nnu=3, prediss=80 is recommended)

end module
