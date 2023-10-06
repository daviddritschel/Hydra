module constants

 !Include all modifiable parameters for use below:
use parameters

! Contains all the non-modifiable parameters as well as all 
! quantities which never change throughout a simulation
! for the suite of casl f90 codes.

 !Grid dimensions +/- 1 & 2:
integer,parameter:: nxp1=nx+1,nxm1=nx-1,nxm2=nx-2
integer,parameter:: nyp1=ny+1,nym1=ny-1

 !Grid dimensions used in write statements:                                                       
integer,parameter:: ngridp=nx*(ny+1),nbytes=4*(ngridp+1)

 !Fine grid used normally in contour -> grid conversion: 
integer,parameter:: mgf=4,nxf=mgf*nx,nyf=mgf*ny
 !mgf:  fine grid/coarse grid ratio (4 is required by subroutine 
 !      coarsen in contours.f90)
 !Fine grid dimensions +/- 1 & 2:
integer,parameter:: nxfm1=nxf-1,nxfm2=nxf-2
integer,parameter:: nyfp1=nyf+1,nyfm1=nyf-1

 !Ultra-fine grid used in contouring: 
integer,parameter:: mgu=16,nxu=mgu*nx,nyu=mgu*ny
 !mgu:  ultra-fine grid/coarse grid ratio (16 is the default)
 !Ultra-fine grid dimensions +/- 1 & 2:
integer,parameter:: nxum1=nxu-1,nxum2=nxu-2
integer,parameter:: nyup1=nyu+1,nyum1=nyu-1

 !Number of divisions of a contour segment for computing db/dx:
integer,parameter:: ndiv=mgu

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
double precision,parameter:: zero=0.d0,one=1.d0 ,two=2.d0,three=3.d0
double precision,parameter:: four=4.d0,five=5.d0,six=6.d0,twelve=12.d0
double precision,parameter:: f12=one/two,    f13=one/three,f23=two/three
double precision,parameter:: f14=one/four,   f32=three/two,f43=four/three
double precision,parameter:: f53=five/three, f56=five/six, f74=7.d0/four
double precision,parameter:: f112=one/twelve,f1112=11.d0/twelve,f16=one/six
double precision,parameter:: f18=one/8.d0,f116=one/16.d0
double precision,parameter:: pi=3.14159265358979d0,twopi=two*pi
double precision,parameter:: small=1.d-12,small3=small*small*small
double precision,parameter:: oms=one-small

 !Domain lengths and inverses:
double precision,parameter:: hlx=f12*ellx,hlxi=one/hlx,xmin=-hlx,xmax=hlx
double precision,parameter:: ycen=(ymin+ymax)/two,elly=ymax-ymin
double precision,parameter:: hly=oms*f12*elly,ybeg=ycen-hly

 !Set number of timesteps between surgery & reset (see Fontane & D. 2009 JCP)
double precision,parameter:: dtfac=pi/25.d0 
integer,parameter:: nstep=nint(8.d0*pi/(10.d0*dtfac)) 

 !Basic constants:
double precision,parameter:: aspect=ellx/elly
double precision,parameter:: glx=ellx/dble(nx),glxi=dble(nx)/ellx
double precision,parameter:: gly=elly/dble(ny),glyi=dble(ny)/elly
double precision,parameter:: garea=glx*gly,dsumi=one/dble(nx*ny)
double precision,parameter:: ftide=twopi/ttide

end module
