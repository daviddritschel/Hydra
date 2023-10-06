module constants

 !Include all modifiable parameters for use below:
use parameters

! Contains all the non-modifiable parameters as well as all 
! quantities which never change throughout a simulation
! for the suite of casl f90 codes.

 !Grid dimensions +/- 1 & 2:
integer,parameter:: nxp1=nx+1,nxm1=nx-1,nxm2=nx-2
integer,parameter:: nyp1=ny+1,nym1=ny-1

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

 !For reading & writing direct access data:
integer,parameter:: ngridp=nx*nyp1,nbytes=4*(ngridp+1)

 !Maximum number of contour levels (used in surgery and congen):
integer,parameter:: nlevm=2000
 !nlevm: up to 2*nlevm contour levels are allowed

 !Maximum number of contour nodes:
integer,parameter:: npm=625*nx*ny
 !Maximum number of contours:
integer,parameter:: nm=npm/20+npm/200
 !Maximum number of nodes on any single contour:
integer,parameter:: nprm=npm/10
 !Maximum number of nodes in any contour level:
integer,parameter:: nplm=npm/2

 !Generic double precision numerical constants: 
double precision,parameter:: zero=0.d0,one=1.d0 ,two=2.d0,three=3.d0
double precision,parameter:: four=4.d0,five=5.d0,six=6.d0,twelve=12.d0
double precision,parameter:: f12=one/two,    f13=one/three,f23=two/three
double precision,parameter:: f14=one/four,   f32=three/two,f43=four/three
double precision,parameter:: f53=five/three, f56=five/six, f74=7.d0/four
double precision,parameter:: f112=one/twelve,f1112=11.d0/twelve,f16=one/six
double precision,parameter:: f18=one/8.d0,f116=one/16.d0
double precision,parameter:: twopi=two*pi
double precision,parameter:: small=1.d-12,small3=small*small*small
double precision,parameter:: oms=one-small

 !Domain lengths and inverses:
double precision,parameter:: hlx=f12*ellx,hlxi=one/hlx,xmin=-hlx,xmax=hlx
double precision,parameter:: ymax=ymin+elly,ycen=(ymin+ymax)/two
double precision,parameter:: hly=oms*f12*elly,ybeg=ycen-hly

 !Fractional thickness of the upper layer (note: h1 + h2 = 1):
double precision,parameter:: h2=one-h1,h1h2=h1*h2
double precision,parameter:: h1inv=one/h1,h2inv=one/h2

 !Constants used to define relative vorticity and interface displacements:
double precision,parameter:: alphac=one-alpha,kdbarsq=kdbar*kdbar
double precision,parameter:: h1kdbarsq=h1*kdbarsq
double precision,parameter:: h1h2kdbarsq=h1h2*kdbarsq
double precision,parameter:: h1h2ackdbarsq=h1h2*alphac*kdbarsq
double precision,parameter:: h1m=h1/(h1+alpha*h2),h2m=alpha*h2/(h1+alpha*h2)

 !Define Rossby deformation wavenumbers (for each mode):
double precision,parameter:: gamma1=f12-sqrt(f14-alphac*h1*h2)
double precision,parameter:: gamma2=f12+sqrt(f14-alphac*h1*h2)
double precision,parameter:: kd1sq=gamma1*kdbarsq, kd1=sqrt(kd1sq)
double precision,parameter:: kd2sq=gamma2*kdbarsq, kd2=sqrt(kd2sq)

 !Define layer -> mode transformation coefficients:
double precision,parameter:: vec11=h1/(one-gamma1)          !mode 1, layer 1
double precision,parameter:: vec12=(h2-gamma1)/(one-gamma1) !mode 1, layer 2
double precision,parameter:: vec21=h1/(h2-gamma2)           !mode 2, layer 1
double precision,parameter:: vec22=one                      !mode 2, layer 2

 !Define mode -> layer transformation coefficients:
double precision,parameter:: determinant=vec11*vec22-vec12*vec21
double precision,parameter:: vect11=vec22/determinant       !layer 1, mode 1
double precision,parameter:: vect12=-vec12/determinant      !layer 1, mode 2
double precision,parameter:: vect21=-vec21/determinant      !layer 2, mode 1
double precision,parameter:: vect22=vec11/determinant       !layer 2, mode 2

 !Initialise coefficients for decomposing energy & enstrophy into modes:
double precision,parameter:: coeff1=h1*vect11**2+alpha*h2*vect21**2
double precision,parameter:: coeff2=h1*vect12**2+alpha*h2*vect22**2

 !Constant used to determine edge velocities below:
double precision,parameter:: momcon=-(h1*vect12+alpha*h2*vect22)/ & 
                                     (h1*vect11+alpha*h2*vect21)

 !Fraction mean PV gradient in upper layer (used only in initialisation):
double precision,parameter:: lambda=one+(one-epsilon)* &
                                    (vec11*gamma2-momcon*vec21*gamma1)/ & 
                                    (vec12*gamma2-momcon*vec22*gamma1)

 !Mean modal velocities (used as boundary conditions on mean zonal velocity):
double precision,parameter:: u2hmean=(one-epsilon)*beta*determinant/ &
                             (kdbarsq*(vec12*gamma2-momcon*vec22*gamma1))
double precision,parameter:: u1hmean=momcon*u2hmean

 !Mean layer velocities (used only for initialisation and diagnostic purposes):
double precision,parameter:: u1lmean=vect11*u1hmean+vect12*u2hmean
double precision,parameter:: u2lmean=vect21*u1hmean+vect22*u2hmean

 !Thermal damping rates:
double precision,parameter:: therm1=rtherm1*h1inv,therm2=rtherm2*h2inv

 !Maximum time step:
double precision,parameter:: dtmax=0.1d0*min(tgsave,tcsave)

 !Set maximum Rossby wave frequency (used in adapt in evolution.f90):
double precision,parameter:: srwfm=beta/((twopi/ellx)**2+kd1sq)

 !Basic constants:
double precision,parameter:: domarea=ellx*elly,aspect=ellx/elly
double precision,parameter:: glx=ellx/dble(nx),glxi=dble(nx)/ellx
double precision,parameter:: gly=elly/dble(ny),glyi=dble(ny)/elly
double precision,parameter:: garea=glx*gly,dsumi=one/dble(nx*ny)
double precision,parameter:: dnxi=one/dble(nx),dnyi=one/dble(ny),hgly=gly/two

 !Logical control variables:
logical,parameter:: heating=((rtherm1 .gt. zero) .or. (rtherm2 .gt. zero))
logical,parameter:: friction=(rekman .gt. zero)
logical,parameter:: damping=(heating .or. friction)
logical,parameter:: barot=(alpha .eq. one)
logical,parameter:: stoch=(eirate .gt. zero)
logical,parameter:: tide=(ftopo .ne. zero)

end module
