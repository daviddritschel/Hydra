module common

 !Import contants, parameters and common arrays needed for inversion etc:
use constants
use contours
use spectral
use generic

 !Velocity field, vorticity & buoyancy (physical):
double precision:: uu(0:ny,0:nxm1),vv(0:ny,0:nxm1)
double precision:: zz(0:ny,0:nxm1),bb(0:ny,0:nxm1)

 !Prognostic fields (spectral):
double precision:: zs(0:nxm1,0:ny),zd(0:nxm1,0:ny)

 !Time, time step, fractions thereof and "twist":
double precision:: t,dt,dt2,dt4,twist

 !Background potential energy calculation flag:
integer:: iene

 !Used for saving data:
integer:: igrids,iconts
logical:: gsave,csave

end module
