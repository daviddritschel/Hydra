module common

 !Module containing all global common areas

 !Import contants, parameters and common arrays:
use constants
use contours
use spectral

!-----------------------------------------------------------------------
 !Define quantities to be preserved between recontouring and evolution:
!-----------------------------------------------------------------------

 !PV residual for recontouring):
double precision:: qr(ng,ng)

 !Velocity, dimensionless height anomaly and vorticity in physical space:
double precision:: uu(ng,ng),vv(ng,ng),hh(ng,ng),zz(ng,ng)

 !Spectral prognostic fields:
double precision:: qs(ng,ng),ds(ng,ng),gs(ng,ng)
double precision:: bxs(ng,ng),bys(ng,ng),bzs(ng,ng)

 !Time stepping parameters:
double precision:: t,dt,dt2,dt4,dt2i,dt4i,twist
integer:: igrids,iconts

end module
