module parameters

! This module contains all the modifiable parameters (except pi below)
! for the suite of imhd caps f90 codes.

 !The number pi (*** non-modifiable ***):
double precision,parameter:: pi=3.141592653589793238462643383279502884197169399375105820974944592307816d0

! ==> Numerical parameters <==
integer,parameter:: nx=N_X,ny=N_Y
integer,parameter:: ncontq=N_CONTQ,npm=NPT_MAX
double precision,parameter:: tsim=T_SIM,tgsave=T_GSAVE,tcsave=T_CSAVE
! nx & ny: grid resolution in x & y
! ncontq : number of PV contour levels over the range Q = <q^2>/<|q|>
!          where |q| > q_rms; the contour interval dq = Q/ncontq
! npm    : maximum number of contour nodes
! tsim   : total duration of the simulation
! tgsave : grid data save time increment 
! tcsave : contour data save time increment (approximate)

! ==> Physical parameters <==
double precision,parameter:: ellx=L_X,elly=L_Y
double precision,parameter:: kd=K_D,beta=PV_GRAD
double precision,parameter:: b0=B_0,eta=MAG_DIFF
double precision,parameter:: rtherm=R_THERM,rekman=R_EKMAN
double precision,parameter:: cdamp=C_DAMP,nnu=POW_HYPER
double precision,parameter:: esrz=E_SRZ,powz=POW_Z,tcz=TC_Z
double precision,parameter:: esra=E_SRA,powa=POW_A,tca=TC_A
integer,parameter:: k0z=K0_Z,k0a=K0_A,iseed=I_SEED
logical,parameter:: tracer=T_RACER
! ellx   : domain width in x (periodic, centred at 0)
! elly   : domain width in y (periodic, centred at 0)
! kd     : Rossby deformation wavenumber associated with the baroclinic mode
! beta   : planetary vorticity gradient
! b0     : mean zonal magnetic field in x direction
! eta    : magnetic diffusivity
! rtherm : thermal damping rate
! rekman : Ekman   damping rate
! cdamp  : this times |zeta|_rms*(k/k_max)^(2*nnu) is the hyperdiffusivity
!          coefficient for qd, where nnu is the power specified above
! esrz   : enstrophy input rate (vorticity)
! powz   : spectral exponent (vorticity)
! tcz    : correlation time (vorticity)
! k0z    : peak wavenumber (vorticity)
! esra   : enstrophy input rate (magnetic potential)
! powa   : spectral exponent (magnetic potential)
! tca    : correlation time (magnetic potential)
! k0a    : peak wavenumber (magnetic potential)
!        - a spectrum proportional to k^{2p-1}*exp(-2(k/k_0)^2)
!        - is used, where p = powz or powa and k_0 = k0z or k0a
! iseed  : seed for initialising forcing
! tracer : logical used to monitor circulation (if true)
!----------------------------------------------------------------

end module
