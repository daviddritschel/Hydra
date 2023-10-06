module parameters

! Module containing all the modifiable parameters (except pi below).

! The number pi (*** non-modifiable ***):
double precision,parameter:: pi=3.141592653589793238462643383279502884197169399375105820974944592307816d0

! ==> Numerical parameters <==
integer,parameter:: ng=N_G,ncont=N_CONTQ
double precision,parameter:: tsim=T_SIM,tgsave=T_GSAVE,tcsave=T_CSAVE
! ng     : inversion grid resolution in both x and y
!          (Note: the domain is a 2*pi periodic box.)
! ncont  : Number of PV jumps (contours) spanning q_min to q_max
! tsim   : total duration of the simulation
! tgsave : grid data save time increment (approximate)
! tcsave : contour data save time increment (approximate)

! ==> Physical parameters <==
double precision,parameter:: cof=COR_FREQ,cgw=C_GW,eta=MAG_DIFF
double precision,parameter:: cdamp=C_DAMP,nnu=POW_HYPER
! cof    : Constant Coriolis frequency f
! cgw    : Short-scale gravity wave speed c
! eta    : magnetic diffusivity
! cdamp  : cdamp*(f+|zeta|_rms)*(k/k_max)^(2*nnu) is the hyperdiffusivity
!          coefficient, where nnu is the power specified above
!----------------------------------------------------------------

end module parameters
