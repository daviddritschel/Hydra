module parameters

! This module contains all the modifiable parameters for the suite of 
! ellipsoidal casl f90 files.

!***Numerical parameters:***
integer,parameter:: ng=N_G, nt=2*ng
integer,parameter:: nperiod=N_PER, nsteps=N_STEP
double precision,parameter:: dt=T_STEP, dq=PV_JUMP
! ng, nt : number of grid boxes in latitude & longitude (inversion grid)
! nperiod: total number of periods to run (each of length dt*nsteps)
! nsteps : number of time steps per period (and per data save)
! dt     : fixed time step (no forcing in this code)
! dq     : PV jump across all contours

!***Physical parameters:***
double precision,parameter:: asp=ASPECT,omega=OMEGA
! asp    : the height:width aspect ratio of the ellipsoidal surface
! omega  : the planetary rotation rate

end module
