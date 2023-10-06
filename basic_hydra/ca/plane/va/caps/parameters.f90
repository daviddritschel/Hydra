module parameters

! Module containing all the modifiable parameters (except pi below).

! The number pi (*** non-modifiable ***):
double precision,parameter:: pi=3.141592653589793238462643383279502884197169399375105820974944592307816d0

! ==> Numerical parameters <==
integer,parameter:: ng=N_G
integer,parameter:: ncont=N_CONTQ,nsamp=N_SAMP
double precision,parameter:: dt=T_STEP,tsim=T_SIM
double precision,parameter:: tgsave=T_GSAVE,tcsave=T_CSAVE
! ng     : inversion grid resolution in both x and y
!          (Note: the domain is a 2*pi periodic box.)
! ncont  : Number of PV jumps (contours) spanning q_min to q_max
! nsamp  : nsamp x nsamp equally spaced sample points are used
!          to record the divergence each time step for a
!          frequency spectrum analysis (post processing)
! dt     : time step (fixed)
! tsim   : total duration of the simulation
! tgsave : grid data save time increment
! tcsave : contour data save time increment

! ==> Physical parameters <==
double precision,parameter:: cof=COR_FREQ,cgw=C_GW,hbar=H_BAR
double precision,parameter:: cdamp=C_DAMP,nnu=POW_HYPER
! cof    : Constant Coriolis frequency
! cgw    : Short-scale gravity wave speed
! hbar   : Mean fluid height H
! cdamp  : This times f is the damping rate on wavenumber ng/2
!----------------------------------------------------------------

end module parameters
