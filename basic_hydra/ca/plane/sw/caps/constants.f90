module constants

 !Module containing all non-modifiable parameters.

use parameters

 !Integer kind parameters used in contours.f90 and congen.f90:
 !  =>  dbleint represents up to 10^18
 !  =>  halfint represents up to 127
integer,parameter:: dbleint=selected_int_kind(16)
integer,parameter:: halfint=selected_int_kind(-1)

 !Size of record used in unformatted writes of real*8 (dble prec) data:
integer,parameter:: nbytes=8*(ng*ng+1)

 !Fine grid used normally in contour -> grid conversion: 
integer,parameter:: mgf=4,ngf=mgf*ng
 !mgf:  fine grid/coarse grid ratio (4 is the standard value)

 !Ultra-fine grid used in contouring: 
integer,parameter:: mgu=16,ngu=mgu*ng
 !mgu:  ultra-fine grid/coarse grid ratio (16 is the standard value)

 !Maximum number of contour levels (used in surgery and congen):
integer,parameter:: nlevm=2000
 !nlevm: up to 2*nlevm contour levels are allowed (very generous)

 !Maximum number of contour nodes:
integer,parameter:: npm=800*ng*ng
 !Maximum number of contours:
integer,parameter:: nm=npm/20+npm/200
 !Maximum number of nodes on any single contour:
integer,parameter:: nprm=npm/10
 !Maximum number of nodes in any contour level:
integer,parameter:: nplm=npm/4

 !For recording divergence at selected grid points (for a frequency spectrum):
integer,parameter:: nsamp2=nsamp*nsamp

 !Generic double precision numerical constants:
double precision,parameter:: zero=0.d0,one=1.d0,two=2.d0,three=3.d0
double precision,parameter:: four=4.d0,six=6.d0,f12=one/two,f13=one/three
double precision,parameter:: f14=one/four,f16=one/six,f32=three/two
double precision,parameter:: twopi=two*pi,thrpi=three*pi,pinv=one/pi
double precision,parameter:: small=1.d-12,small3=small*small*small
double precision,parameter:: oms=one-small

 !Grid constants:
double precision,parameter:: domarea=twopi*twopi
double precision,parameter:: gl=twopi/dble(ng),gli=dble(ng)/twopi
double precision,parameter:: garea=gl*gl,dsumi=one/dble(ng*ng)

! Time step related parameters:
double precision,parameter:: dt2=dt*f12,dt4=dt*f14
double precision,parameter:: dt2i=one/dt2,dt4i=one/dt4

! Squared gravity wave speed:
double precision,parameter:: csq=cgw**2

end module
