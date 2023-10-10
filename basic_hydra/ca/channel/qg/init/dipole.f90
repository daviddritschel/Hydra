program dipole
! Initialises a dipolar vortex in beta (from Larichev & Reznik 1976)

use contours

implicit double precision(a-h,o-z)

double precision:: qq(0:ny,0:nxm1)
double precision:: qa(0:nyfp1,0:nxfm1)
double precision:: ygf(0:nyf),bety(0:ny)
double precision:: c,k,ksq

write(*,*) ' Propagation speed, c? '
read(*,*) c

!-----------------------------------------------------------------
 !Initialise constants and arrays for contour advection:
call init_contours

!-----------------------------------------------------------------
psq=beta/c+kdsq
p=sqrt(psq)
a=p*bessk1(p)/bessk(p,2)

k=4.d0
dk=1.d0
do while (abs(dk).gt. 1.d-12)
  rj0=bessj0(k)
  rj1=bessj1(k)
  rj2=2.d0*rj1/k-rj0
  dk=-(k*rj1+a*rj2)/(k*rj0+a*(rj1-2.d0*rj2/k))
  k=k+dk
enddo
ksq=k**2

!-----------------------------------------------------------------
ai=ksq+kdsq
bi=ksq+psq
di=bessj1(k)
hi=psq/ksq
ae=(kdsq-psq)/bessk1(p)

do iy=0,nyf
  ygf(iy)=ymin+glyf*dble(iy)
enddo

do ix=0,nxfm1
  do iy=0,nyf
    r=sqrt(xgf(ix)**2+ygf(iy)**2)+1.d-15
    if (r .ge. one) then
      f=ae*bessk1(p*r)/r
    else
      f=ai*(one+hi*(one-bessj1(k*r)/(r*di)))-bi
    endif
    qa(iy,ix)=c*ygf(iy)*f
  enddo
enddo

 !Ensure average PV is zero:
qavg=zero

 !Average the fine-grid PV field in qa to the coarser grid (ny,nx):
call coarsen(qa,qq)

 !Add beta*y:
do iy=0,ny
  bety(iy)=beta*(ymin+gly*dble(iy))
enddo

do ix=0,nxm1
  do iy=0,ny
    qq(iy,ix)=qq(iy,ix)+bety(iy)
  enddo
enddo

 !Write data:
open(11,file='qq_init.r8',form='unformatted', &
    & access='direct',status='replace',recl=2*nbytes)
write(11,rec=1) zero,qq
close(11)

contains

!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

double precision function bessj0(x)
! Evaluates the Bessel function J_0(x) for any real x

implicit double precision (a-h,o-z)

data p1,p2,p3,p4,p5/1.d0,-0.1098628627d-2,0.2734510407d-4, &
                    -0.2073370639d-5,0.2093887211d-6/
data q1,q2,q3,q4,q5/-0.1562499995d-1,0.1430488765d-3, &
                    -0.6911147651d-5,0.7621095161d-6,-0.934945152d-7/
data r1,r2,r3,r4,r5,r6/57568490574.d0,-13362590354.d0,651619640.7d0, &
                       -11214424.18d0,77392.33017d0,-184.9052456d0/
data s1,s2,s3,s4,s5,s6/57568490411.d0,1029532985.d0, &
                       9494680.718d0,59272.64853d0,267.8532712d0,1.d0/

if (abs(x) .lt. 8.d0) then
  y=x**2
  bessj0=(r1+y*(r2+y*(r3+y*(r4+y*(r5+y*r6))))) &
        /(s1+y*(s2+y*(s3+y*(s4+y*(s5+y*s6)))))
else
  ax=abs(x)
  z=8.d0/ax
  y=z**2
  xx=ax-0.785398164d0
  bessj0=sqrt(0.636619772d0/ax)* &
           (cos(xx)*(p1+y*(p2+y*(p3+y*(p4+y*p5)))) &
         -z*sin(xx)*(q1+y*(q2+y*(q3+y*(q4+y*q5)))))
endif

return
end function bessj0

!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

double precision function bessj1(x)
! Evaluates the Bessel function J_1(x) for any real x

implicit double precision (a-h,o-z)

data r1,r2,r3,r4,r5,r6/72362614232.d0,-7895059235.d0,242396853.1d0, &
                       -2972611.439d0,15704.48260d0,-30.16036606d0/
data s1,s2,s3,s4,s5,s6/144725228442.d0,2300535178.d0, &
                       18583304.74d0,99447.43394d0,376.9991397d0,1.d0/
data p1,p2,p3,p4,p5/1.d0,0.183105d-2,-0.3516396496d-4, &
                    0.2457520174d-5,-0.240337019d-6/ 
data q1,q2,q3,q4,q5/0.04687499995d0,-0.2002690873d-3, &
                    0.8449199096d-5,-0.88228987d-6,0.105787412d-6/

if (abs(x) .lt. 8.d0) then
  y=x**2
  bessj1=x*(r1+y*(r2+y*(r3+y*(r4+y*(r5+y*r6))))) &
          /(s1+y*(s2+y*(s3+y*(s4+y*(s5+y*s6)))))
else
  ax=abs(x)
  z=8.d0/ax
  y=z**2
  xx=ax-2.356194491d0
  bessj1=sqrt(0.636619772d0/ax)* &
           (cos(xx)*(p1+y*(p2+y*(p3+y*(p4+y*p5)))) &
         -z*sin(xx)*(q1+y*(q2+y*(q3+y*(q4+y*q5)))))*sign(1.d0,x)
endif

return
end function bessj1

!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

double precision function bessj2(x)
! Evaluates the Bessel function J_2(x) for any real x

implicit double precision (a-h,o-z)

if (abs(x) .gt. 0.d0) then
  bessj2=2.d0*bessj1(x)/x-bessj0(x)
else
  bessj2=0.d0
endif

return
end function bessj2

!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

double precision function bessi(x,m)
! Evaluates the modified Bessel functions I_m(x) for any real x.  

implicit double precision (a-h,o-z)
parameter (acc=40.d0)
parameter (bigno=1.d10,bigni=1.d-10)

if (m .eq. 0) then
  bessi=bessi0(x)
else if (m .eq. 1) then
  bessi=bessi1(x)
else
  if (x .eq. 0.d0) then
    bessi = 0.d0
  else
    cc = 2.d0/x
    bip = 0.d0
    ans = 0.d0
    bi = 1.d0
    do j = 2*(m + int(sqrt(acc*dble(m)))),1,-1
      bim = bip + cc*dble(j)*bi
      bip = bi
      bi  = bim
      if (bi .gt. bigno) then
        ans = ans*bigni
        bi  =  bi*bigni
        bip = bip*bigni
      endif
      if (j .eq. m) ans = bip
    enddo
    bessi = ans*bessi0(x)/bi
  endif
endif

return
end function bessi

!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

double precision function bessi0(x)
! Evaluates the modified Bessel function I_0(x) for any real x

implicit double precision (a-h,o-z)

data p1,p2,p3,p4,p5,p6,p7/1.0d0,3.5156229d0,3.0899424d0, &
                          1.2067492d0,0.2659732d0,0.360768d-1,0.45813d-2/
data q1,q2,q3,q4,q5,q6,q7,q8,q9/0.39894228d0,0.1328592d-1, &
                    0.225319d-2,-0.157565d-2,0.916281d-2,-0.2057706d-1, &
                   0.2635537d-1,-0.1647633d-1,0.392377d-2/

if (abs(x) .lt. 3.75d0) then
  y = (x/3.75d0)**2
  bessi0 = p1+y*(p2+y*(p3+y*(p4+y*(p5+y*(p6+y*p7)))))
else
  ax = abs(x)
  y = 3.75d0/ax
  bessi0 = (exp(ax)/sqrt(ax))* &
           (q1+y*(q2+y*(q3+y*(q4+y*(q5+y*(q6+y*(q7+y*(q8+y*q9))))))))
endif

return
end function bessi0

!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

double precision function bessi1(x)
! Evaluates the modified Bessel function I_1(x) for any real x

implicit double precision (a-h,o-z)

!save p1,p2,p3,p4,p5,p6,p7,q1,q2,q3,q4,q5,q6,q7,q8,q9

data p1,p2,p3,p4,p5,p6,p7/0.5d0,0.87890594d0,0.51498869d0, &
             0.15084934d0,0.2658733d-1,0.301532d-2,0.32411d-3/
data q1,q2,q3,q4,q5,q6,q7,q8,q9/0.39894228d0,-0.3988024d-1, &
                   -0.362018d-2,0.163801d-2,-0.1031555d-1,0.2282967d-1, &
                  -0.2895312d-1,0.1787654d-1,-0.420059d-2/

if (abs(x) .lt. 3.75d0) then
  y = (x/3.75d0)**2
  bessi1 = x*(p1+y*(p2+y*(p3+y*(p4+y*(p5+y*(p6+y*p7))))))
else
  ax = abs(x)
  y = 3.75/ax
  bessi1 = (exp(ax)/sqrt(ax))* &
           (q1+y*(q2+y*(q3+y*(q4+y*(q5+y*(q6+y*(q7+y*(q8+y*q9))))))))
endif

return
end function bessi1

!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

double precision function bessk(x,m)
! Evaluates the modified Bessel functions K_m(x) for any real x > 0.  

implicit double precision (a-h,o-z)

if (m .eq. 0) then
  bessk = bessk0(x)
else if (m .eq. 1) then
  bessk = bessk1(x)
else
  cc = 2.d0/x
  bkm = bessk0(x)
  bk  = bessk1(x)
  do j = 1, m-1
    bkp = bkm + cc*dble(j)*bk
    bkm = bk
    bk  = bkp
  enddo
  bessk = bk
endif

return
end function bessk

!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

double precision function bessk0(x)
! Evaluates the modified Bessel function K_0(x) for real x > 0
! (makes use of bessi0.F)

implicit double precision (a-h,o-z)

data p1,p2,p3,p4,p5,p6,p7/-0.57721566d0,0.42278420d0,0.23069756d0, &
                          0.3488590d-1,0.262698d-2,0.10750d-3,0.74d-5/
data q1,q2,q3,q4,q5,q6,q7/1.25331414d0,-0.7832358d-1,0.2189568d-1, &
                          -0.1062446d-1,0.587872d-2,-0.251540d-2,0.53208d-3/

if (x .le. 2.d0) then
  y = x*x/4.d0
  bessk0 = -log(x/2.d0)*bessi0(x) &
           +(p1+y*(p2+y*(p3+y*(p4+y*(p5+y*(p6+y*p7))))))
else
  y = 2.d0/x
  bessk0 = (exp(-x)/sqrt(x))* &
           (q1+y*(q2+y*(q3+y*(q4+y*(q5+y*(q6+y*q7))))))
endif

return
end function bessk0

!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

double precision function bessk1(x)
! Evaluates the modified Bessel function K_1(x) for real x > 0
! (makes use of bessi1.F)

implicit double precision (a-h,o-z)

!save p1,p2,p3,p4,p5,p6,p7,q1,q2,q3,q4,q5,q6,q7

data p1,p2,p3,p4,p5,p6,p7/1.0d0,       0.15443144d0,-0.67278579d0, &
            -0.18156897d0,-0.1919402d-1,-0.110404d-2,  -0.4686d-4/
data q1,q2,q3,q4,q5,q6,q7/1.25331414d0,0.23498619d0,-0.3655620d-1, &
            0.1504268d-1, -0.780353d-2, 0.325614d-2, -0.68245d-3/

if (x .le. 2.d0) then
  y = x*x/4.d0
  bessk1 = (log(x/2.d0)*bessi1(x)) + (1.d0/x)* &
           (p1+y*(p2+y*(p3+y*(p4+y*(p5+y*(p6+y*p7))))))
else
  y = 2.d0/x
  bessk1 = (exp(-x)/sqrt(x))* &
           (q1+y*(q2+y*(q3+y*(q4+y*(q5+y*(q6+y*q7))))))
endif

return
end function bessk1

end program dipole
