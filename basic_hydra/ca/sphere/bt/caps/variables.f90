module variables

 !Include all variables which may change during the course of a simulation
 !and need to be preserved through recontouring

implicit none

 !Time-stepping variables:
double precision:: t,dt,dt2,dt3,dt6
 !Time to integrate to over the next period:
double precision:: tmax
 !Indices for data writes:
integer:: iref,idump

end module
