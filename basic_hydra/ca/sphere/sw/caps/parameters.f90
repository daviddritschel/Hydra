module parameters

! Module containing all the modifiable parameters (except pi below).

! The number pi (*** non-modifiable ***):
double precision,parameter:: pi=3.141592653589793238462643383279502884197169399375105820974944592307816d0

! ==> Numerical parameters <==
integer,parameter:: ng=N_G, nt=2*ng
integer,parameter:: ncont=N_CONTQ
double precision,parameter:: dt=T_STEP,tsim=T_SIM
double precision,parameter:: tgsave=T_GSAVE,tcsave=T_CSAVE
! ng     : number of  latitudinal grid points (using a half grid)
! nt     : number of longitudinal grid points (always equal to 2*ng)
! ncont  : number of PV jumps used to represent f = 2*Omega*sin(lat)
! dt     : time step (fixed)
! tsim   : total duration of the simulation
! tgsave : grid data save time increment
! tcsave : contour data save time increment

! ==> Physical parameters <==
double precision,parameter:: fpole=COR_FREQ,cgw=C_GW
double precision,parameter:: Rocp=R_CP
double precision,parameter:: drate=D_RATE,nnu=POW_HYPER
double precision,parameter:: rth=R_THERM,rek=R_EKMAN
double precision,parameter:: brms=B_RMS,tb=T_B
integer,parameter:: nbeg=N_BEG, nend=N_END
integer,parameter:: iseed=I_SEED
! fpole      : Coriolis frequency f at the north pole (same as 2*Omega)
! cgw        : short-scale gravity wave speed c (= fpole/k_d = fpole*L_d)
! Rocp       : gas constant / specific heat at constant pressure (R/c_p)
! drate      : drate*fpole is the hyperviscous damping rate on wavenumber ng
! nnu        : power of hyperviscosity
! rth        : thermal damping rate
! rek        : Ekman damping rate
! brms       : r.m.s. topographic amplitude
! tb         : Markovian de-correlation time of the topographic forcing
! nbeg, nend : range in order n of spherical harmonics Y_n^m used
! iseed      : a random number seed used for topographic forcing

end module parameters
