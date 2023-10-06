module common

 !Import contants, parameters and common arrays needed for inversion etc:
use constants
use variables
use contours
use spectral
use generic

 !Define quantities which need to be preserved between recontouring and evolution:

 !Gridded buoyancy & vorticity fields:
double precision:: bb(0:ny,0:nx),zs(0:ny,0:nx),zd(0:ny,0:nx)

 !Time stepping parameters:
double precision:: dtmax,tfin,tgrid

 !Background potential energy calculation flag:
integer:: iene

end module
