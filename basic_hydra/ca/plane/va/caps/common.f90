module common

 !Module containing all global common areas

 !Import contants, parameters and common arrays:
use constants
use contours
use spectral

!-----------------------------------------------------------------------
 !Define quantities to be preserved between recontouring and evolution:
!-----------------------------------------------------------------------

 !PV anomaly (full and residual for recontouring) and relative vorticity:
double precision:: qq(ng,ng),qr(ng,ng),zz(ng,ng)

 !Velocity field, dimensionless height anomaly & non-hydrostatic pressure:
double precision:: uu(ng,ng),vv(ng,ng),hh(ng,ng),ppn(ng,ng)

 !Spectral prognostic fields:
double precision:: qs(ng,ng),ds(ng,ng),gs(ng,ng)

 !Time and twist parameter (used for contour surgery):
double precision:: t,twist

 !Number of time steps between grid and contour saves:
double precision:: ngsave,ncsave

end module
