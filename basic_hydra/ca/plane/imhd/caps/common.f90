module common

 !Import contants, parameters and common arrays needed for inversion etc:
use constants
use contours
use spectral

 !Define quantities to be preserved between recontouring and evolution:

 !PV residue:
double precision:: qr(ny,nx)

 !Spectral prognostic fields (note array order):
double precision:: aa(nx,ny),qs(nx,ny)

 !Rates of change of vorticity and magnetic potential due to forcing (if used):
double precision:: dzdt(ny,nx),dadt(ny,nx)

 !Time stepping parameters:
double precision:: t,dt,dt2,dt3,dt6,tgrid,twist
integer:: igrids

end module
