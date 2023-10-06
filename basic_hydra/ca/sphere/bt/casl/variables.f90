module variables

use parameters

 !Include all variables which may change during the course of a simulation:

implicit none

 !Relative vorticity field:
double precision:: qq(ng,nt)

 !Gridded velocity field (including poles):
double precision:: uu(0:ng+1,nt),vv(0:ng+1,nt)

 !The current time:
double precision:: t

 !Index for direct access data writes:
integer:: idump

end module
