module variables

 !Include all variables which may change during the course of a simulation.

 !Time & time step:
double precision:: t,dt,hfdt,qudt
 !Domain tilt angle and its cosine and sine:
double precision:: theta,ctheta,stheta
integer:: igrids

end module
