module common

 !Import contants, parameters and common arrays needed for inversion etc:
use constants
use variables
use spectral
use contours

 !Define quantities to be preserved between recontouring and evolution:

 !PV fields:
double precision:: qs(ng,nt),qd(ng,nt)

end module
