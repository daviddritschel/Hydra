program slug

use constants

! This routine sets up initial buoyancy and vorticity fields for
! casl in a rectangular aperiodic domain.
! The initial condition is a buoyancy anomaly, with vorticity. 
! Adapted by Gabriel Rooney (UKMO) on 8 Feb 2023 from slug.f90

implicit double precision (a-h,o-z)

 !Local vorticity and buoyancy arrays:
double precision:: zz(0:ny,0:nxm1),bb(0:ny,0:nxm1)

! local constants
double precision:: smallh, sigma
double precision:: smallh_alt_g, smallh_alt_omega
double precision:: sigma_v2
double precision:: x_v1, x_v2
double precision:: theta
integer:: expt

! index of vortex experiments
! with suggested domain sizes (WxH)
!
write(*,*) ' Choose one of the following experiments:'
write(*,*) ' 1 ( 6x3,  512x256) t= 40 : single balanced vortex'
write(*,*) ' 2 ( 6x3,  512x256) t= 40 : single balanced vortex but without sign changes'
write(*,*) ' 3 ( 6x3,  512x256) t= 40 : vortex with b and zeta fields out of balance'
write(*,*) ' 4 (12x3, 1024x256) t= 60 : colliding balanced vortices, m != 0'
write(*,*) ' 5 (12x3, 1024x256) t= 60 : colliding balanced vortices, m = 0'
write(*,*) ' 6 ( 6x3,  512x256) t= 40 : single "shielded" vortex (h > J1 root)'
write(*,*) ' 7 (12x3, 1024x256) t=120 : chasing balanced vortices, m != 0'
read(*,*) expt

write(*,*)
write(*,'(a,2(f6.2))') ' We consider a domain of height and length: ',elly,ellx
write(*,'(a,2(i6))') ' The y and x grid resolution is: ',ny,nx
gridrat=dble(nx)*elly/(dble(ny)*ellx)
write(*,'(a,2(f6.2))') ' Note, the y:x grid length ratio is: ',gridrat
write(*,*)

! domain x coord seems to be centred on x=0
xc = 0.0  ! centre in x of the domain
w  = 1.0  ! width 1
h  = 1.0  ! height 1
d  = 0.0  ! not used

! NOTE that the hydra fields of buoyancy and vorticity
! are the negative of those in my equations.

! CHOOSE WHICH EXPERIMENT
! hardwired here
expt = 1

! VALUES OF PARAMETERS

! default (balanced vortex)
smallh = 2.0
sigma = 1.0

! intentionally imbalanced vortex
! because the values of h are different
! for buoyancy and vorticity
smallh_alt_g     = 2.5
smallh_alt_omega = 1.5

! sigma value of 2nd vortex in collision
! note opposite sign to first vortex
sigma_v2 = -0.8

! sigma value of 2nd vortex in chase / catch-up
! note same sign as first vortex
sigma_v3 = 0.5

! initial x positions of colliding vortices
!                 and of chasing vortices
x_v1 = -3.0
x_v2 =  3.0

! value of h bigger than the first root of Ji(h)
smallh_big = 5.0

 !Set up buoyancy and vorticity distributions:
! gl stands for grid length
! for instance glx = ellx / dble(nx)

if (expt == 1) then

! EXPT1 : the steady solution
   
   do ix=0,nxm1
      x=xmin+glx*dble(ix)-xc ! x position relative to xc
      do iy=0,ny
         y=gly*dble(iy) ! y position

         r=sqrt(x**2+y**2)

         if (r <= 1) then
            theta = atan2( y, x )

            zz(iy,ix) = - sigma &
                 * ( bessel_jn(1, smallh * r) / bessel_jn(1, smallh) ) &
                 * sin( theta )
            bb(iy,ix) = ( sigma / smallh )**2 &
                 * ( bessel_jn(1, smallh * r) / bessel_jn(1, smallh) - r ) &
                 * sin( theta )
            ! NOW change sign compared to my conventions
            zz(iy,ix) = -1.0 * zz(iy,ix)
            bb(iy,ix) = -1.0 * bb(iy,ix)
         else
            zz(iy,ix) = 0.0
            bb(iy,ix) = 0.0
         endif

      enddo
   enddo

else if (expt == 2) then

! EXPT2 : the steady solution but without the sign changes
!         so the vortex should be positively buoyant;
!         what happens then?
   
   do ix=0,nxm1
      x=xmin+glx*dble(ix)-xc ! x position relative to xc
      do iy=0,ny
         y=gly*dble(iy) ! y position

         r=sqrt(x**2+y**2)

         if (r <= 1) then
            theta = atan2( y, x )

            zz(iy,ix) = - sigma &
                 * ( bessel_jn(1, smallh * r) / bessel_jn(1, smallh) ) &
                 * sin( theta )
            bb(iy,ix) = ( sigma / smallh )**2 &
                 * ( bessel_jn(1, smallh * r) / bessel_jn(1, smallh) - r ) &
                 * sin( theta )
            ! DON'T change sign compared to my conventions
         else
            zz(iy,ix) = 0.0
            bb(iy,ix) = 0.0
         endif

      enddo
   enddo

else if (expt == 3) then
   
! EXPT3 : the steady solution but using different values of h
!         for buoyancy and vorticity, so it shouldn't be properly balanced

   do ix=0,nxm1
      x=xmin+glx*dble(ix)-xc ! x position relative to xc
      do iy=0,ny
         y=gly*dble(iy) ! y position

         r=sqrt(x**2+y**2)

         if (r <= 1) then
            theta = atan2( y, x )

            ! use different h values to unbalance the fields
            zz(iy,ix) = - sigma &
                 * ( bessel_jn(1, smallh_alt_omega * r) / bessel_jn(1, smallh_alt_omega) ) &
                 * sin( theta )
            bb(iy,ix) = ( sigma / smallh_alt_g )**2 &
                 * ( bessel_jn(1, smallh_alt_g * r) / bessel_jn(1, smallh_alt_g) - r ) &
                 * sin( theta )
            ! NOW change sign compared to my conventions
            zz(iy,ix) = -1.0 * zz(iy,ix)
            bb(iy,ix) = -1.0 * bb(iy,ix)
         else
            zz(iy,ix) = 0.0
            bb(iy,ix) = 0.0
         endif

      enddo
   enddo

else if (expt == 4) then
   
! EXPT4 : two balanced m!=0 vortices with opposite signed sigmas
!         so they travel towards each other and collide

   ! initialise first vortex
   do ix=0,nxm1
      x=xmin+glx*dble(ix)-x_v1 ! x position relative to x_v1
      do iy=0,ny
         y=gly*dble(iy) ! y position

         r=sqrt(x**2+y**2)

         if (r <= 1) then
            theta = atan2( y, x )

            zz(iy,ix) = - sigma &
                 * ( bessel_jn(1, smallh * r) / bessel_jn(1, smallh) ) &
                 * sin( theta )
            bb(iy,ix) = ( sigma / smallh )**2 &
                 * ( bessel_jn(1, smallh * r) / bessel_jn(1, smallh) - r ) &
                 * sin( theta )
            ! NOW change sign compared to my conventions
            zz(iy,ix) = -1.0 * zz(iy,ix)
            bb(iy,ix) = -1.0 * bb(iy,ix)
         else
            ! zero the background for the first vortex
            zz(iy,ix) = 0.0
            bb(iy,ix) = 0.0
         endif

      enddo
   enddo

   ! initialise second vortex
   do ix=0,nxm1
      x=xmin+glx*dble(ix)-x_v2 ! x position relative to x_v2
      do iy=0,ny
         y=gly*dble(iy) ! y position

         r=sqrt(x**2+y**2)

         if (r <= 1) then
            theta = atan2( y, x )

            zz(iy,ix) = - sigma_v2 &
                 * ( bessel_jn(1, smallh * r) / bessel_jn(1, smallh) ) &
                 * sin( theta )
            bb(iy,ix) = ( sigma_v2 / smallh )**2 &
                 * ( bessel_jn(1, smallh * r) / bessel_jn(1, smallh) - r ) &
                 * sin( theta )
            ! NOW change sign compared to my conventions
            zz(iy,ix) = -1.0 * zz(iy,ix)
            bb(iy,ix) = -1.0 * bb(iy,ix)
         endif
         ! DON'T zero the background for the second vortex

      enddo
   enddo

else if (expt == 5) then
   
! EXPT5 : two balanced m=0 vortices with opposite signed sigmas
!         so they travel towards each other and collide

   ! initialise first vortex
   do ix=0,nxm1
      x=xmin+glx*dble(ix)-x_v1 ! x position relative to x_v1
      do iy=0,ny
         y=gly*dble(iy) ! y position

         r=sqrt(x**2+y**2)

         if (r <= 1) then
            theta = atan2( y, x )

            zz(iy,ix) = - sigma * r * sin( theta )
            bb(iy,ix) = ( sigma**2 / 8.0 ) &
                         * ( r - r**3 ) * sin( theta )
            ! NOW change sign compared to my conventions
            zz(iy,ix) = -1.0 * zz(iy,ix)
            bb(iy,ix) = -1.0 * bb(iy,ix)
         else
            ! zero the background for the first vortex
            zz(iy,ix) = 0.0
            bb(iy,ix) = 0.0
         endif

      enddo
   enddo

   ! initialise second vortex
   do ix=0,nxm1
      x=xmin+glx*dble(ix)-x_v2 ! x position relative to x_v2
      do iy=0,ny
         y=gly*dble(iy) ! y position

         r=sqrt(x**2+y**2)

         if (r <= 1) then
            theta = atan2( y, x )

            zz(iy,ix) = - sigma_v2 * r * sin( theta )
            bb(iy,ix) = ( sigma_v2**2 / 8.0 ) &
                         * ( r - r**3 ) * sin( theta )
            ! NOW change sign compared to my conventions
            zz(iy,ix) = -1.0 * zz(iy,ix)
            bb(iy,ix) = -1.0 * bb(iy,ix)
         endif
         ! DON'T zero the background for the second vortex

      enddo
   enddo

else if (expt == 6) then

! EXPT6 : as expt 1 but witl a bigger value of h
   
   do ix=0,nxm1
      x=xmin+glx*dble(ix)-xc ! x position relative to xc
      do iy=0,ny
         y=gly*dble(iy) ! y position

         r=sqrt(x**2+y**2)

         if (r <= 1) then
            theta = atan2( y, x )

            zz(iy,ix) = - sigma &
                 * ( bessel_jn(1, smallh_big * r) / bessel_jn(1, smallh_big) ) &
                 * sin( theta )
            bb(iy,ix) = ( sigma / smallh_big )**2 &
                 * ( bessel_jn(1, smallh_big * r) / bessel_jn(1, smallh_big) - r ) &
                 * sin( theta )
            ! NOW change sign compared to my conventions
            zz(iy,ix) = -1.0 * zz(iy,ix)
            bb(iy,ix) = -1.0 * bb(iy,ix)
         else
            zz(iy,ix) = 0.0
            bb(iy,ix) = 0.0
         endif

      enddo
   enddo


else if (expt == 7) then
   
! EXPT7 : two balanced m!=0 vortices with same signed sigmas
!         but the "rear" one has a bigger value, hence faster
!         so it should catch the other one up

   ! initialise first vortex
   do ix=0,nxm1
      x=xmin+glx*dble(ix)-x_v1 ! x position relative to x_v1
      do iy=0,ny
         y=gly*dble(iy) ! y position

         r=sqrt(x**2+y**2)

         if (r <= 1) then
            theta = atan2( y, x )

            zz(iy,ix) = - sigma &
                 * ( bessel_jn(1, smallh * r) / bessel_jn(1, smallh) ) &
                 * sin( theta )
            bb(iy,ix) = ( sigma / smallh )**2 &
                 * ( bessel_jn(1, smallh * r) / bessel_jn(1, smallh) - r ) &
                 * sin( theta )
            ! NOW change sign compared to my conventions
            zz(iy,ix) = -1.0 * zz(iy,ix)
            bb(iy,ix) = -1.0 * bb(iy,ix)
         else
            ! zero the background for the first vortex
            zz(iy,ix) = 0.0
            bb(iy,ix) = 0.0
         endif

      enddo
   enddo

   ! initialise second vortex
   do ix=0,nxm1
      x=xmin+glx*dble(ix)-x_v2 ! x position relative to x_v2
      do iy=0,ny
         y=gly*dble(iy) ! y position

         r=sqrt(x**2+y**2)

         if (r <= 1) then
            theta = atan2( y, x )

            zz(iy,ix) = - sigma_v3 &
                 * ( bessel_jn(1, smallh * r) / bessel_jn(1, smallh) ) &
                 * sin( theta )
            bb(iy,ix) = ( sigma_v3 / smallh )**2 &
                 * ( bessel_jn(1, smallh * r) / bessel_jn(1, smallh) - r ) &
                 * sin( theta )
            ! NOW change sign compared to my conventions
            zz(iy,ix) = -1.0 * zz(iy,ix)
            bb(iy,ix) = -1.0 * bb(iy,ix)
         endif
         ! DON'T zero the background for the second vortex

      enddo
   enddo

endif

 !Write vorticity distribution to file:
open(20,file='zz_init.r8',form='unformatted', &
    & access='direct',status='replace',recl=2*nbytes)
write(20,rec=1) zero,zz
close(20)

 !Write buoyancy distribution to file:
open(20,file='bb_init.r8',form='unformatted', &
    & access='direct',status='replace',recl=2*nbytes)
write(20,rec=1) zero,bb
close(20)

 !Write information file:
open(12,file='input_for_modon',status='replace')
write(12,*) expt, ' !Experiment number --- see modon.f90'
close(12)
      
end program
