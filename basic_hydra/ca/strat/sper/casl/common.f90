module common

 !Import contants, parameters and common arrays needed for inversion etc:
use constants
use variables
use contours
use spectral
use generic

 !Define quantities which need to be preserved between recontouring and evolution:

 !Gridded buoyancy & vorticity fields:
double precision:: bb(0:ny,0:nxm1),zs(0:ny,0:nxm1),zd(0:ny,0:nxm1)

 !For semi-Lagrangian advection: 
double precision:: xig(0:nxm1),yig(0:ny)
double precision,parameter:: xigmax=dble(nx),yigmax=dble(ny) 

 !Time stepping parameters:
double precision:: dtmax,tfin,tgrid

 !Background potential energy calculation flag:
integer:: iene

end module
