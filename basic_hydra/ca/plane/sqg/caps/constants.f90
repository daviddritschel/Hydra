module constants

 !Include all modifiable parameters for use below:
use parameters

! Contains all the non-modifiable parameters as well as all 
! quantities which never change throughout a simulation
! for the suite of casl f90 codes.

 !Integer kind parameters used in contours.f90 and 
 !congen.f90: dbleint represents up to 10^18
 !            halfint represents up to 127
integer,parameter:: dbleint=selected_int_kind(16)
integer,parameter:: halfint=selected_int_kind(-1)

integer,parameter:: nxm1=nx-1,nym1=ny-1

 !Grid dimensions used in write statements:
integer,parameter:: ngridp=nx*ny,nbytes=4*(ngridp+1)

 !Fine grid used normally in contour -> grid conversion: 
integer,parameter:: mgf=4,nxf=mgf*nx,nyf=mgf*ny
 !mgf:  fine grid/coarse grid ratio (4 is required by subroutine 
 !      coarsen in contours.f90)

 !Ultra-fine grid used in contouring: 
integer,parameter:: mgu=16,nxu=mgu*nx,nyu=mgu*ny
 !mgu:  ultra-fine grid/coarse grid ratio (16 is the default)

 !Maximum number of contour levels (used in surgery and congen):
integer,parameter:: nlevm=2000
 !nlevm: up to 2*nlevm contour levels are allowed

 !Maximum number of contour nodes:
integer,parameter:: npm=400*nx*ny
 !Maximum number of contours:
integer,parameter:: nm=npm/20+npm/200
 !Maximum number of nodes on any single contour:
integer,parameter:: nprm=npm/10
 !Maximum number of nodes in any contour level:
integer,parameter:: nplm=npm/4

 !Generic double precision numerical constants: 
double precision,parameter:: zero=0.d0,one=1.d0,two=2.d0,three=3.d0
double precision,parameter:: four=4.d0,six=6.d0,f12=one/two,f13=one/three
double precision,parameter:: f14=one/four,f16=one/six
double precision,parameter:: twopi=two*pi
double precision,parameter:: small=1.d-12,small3=small*small*small
double precision,parameter:: oms=one-small

 !Domain lengths and inverses:
double precision,parameter:: hlx=f12*ellx,hlxi=one/hlx,xmin=-hlx,xmax=hlx
double precision,parameter:: hly=f12*elly,hlyi=one/hly,ymin=-hly,ymax=hly

 !Basic constants:
double precision,parameter:: domarea=ellx*elly,aspect=ellx/elly
double precision,parameter:: glx=ellx/dble(nx),glxi=dble(nx)/ellx
double precision,parameter:: gly=elly/dble(ny),glyi=dble(ny)/elly
double precision,parameter:: garea=glx*gly,dsumi=one/dble(nx*ny)

 !Control parameters:
logical,parameter:: tracer=(kappa > zero)

end module
