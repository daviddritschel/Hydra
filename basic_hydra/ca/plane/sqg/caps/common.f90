module common

 !Import contants, parameters and common arrays needed for inversion etc:
use constants
use variables
use contours
use spectral

 !Define quantities to be preserved between recontouring and evolution:

 !buoyancy residue:
double precision:: qr(ny,nx)

 !Spectral prognostic field (note array order):
double precision:: qs(nx,ny),qspre(nx,ny)

 !tracer (if present):
double precision,allocatable,dimension(:,:):: cs,cspre

 !Time stepping parameters:
double precision:: dtmax,tfin,tgrid

end module
