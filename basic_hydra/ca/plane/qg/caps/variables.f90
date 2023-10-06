module variables

 !Include all variables which may change during the course of a simulation.

 !Time, time step & fractional time steps:
double precision:: t,dt,dt2,dt3,dt6

 !Twist parameter for timing of surgery and (when forcing) the total number
 !of vortices added:
double precision:: twist,totnvor

!Counters for writing direct-access data (gridded and contours):
integer:: igrec,icrec

end module
