module common

 !Import contants, parameters and common arrays needed for inversion etc:
use constants
use variables
use contours
use spectral
use generic

 !Define quantities to be preserved between recontouring and evolution:

 !buoyancy residue:
double precision:: qr(ny,nx)

 !Spectral prognostic fields (note array order):
double precision:: qs(nx,ny)

 !Time stepping parameters:
double precision:: dtmax,tfin,tgrid

end module
