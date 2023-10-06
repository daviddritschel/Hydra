program stair
!  ---------------------------------------------------------------------------
!  |  Sorts the PV in each layer as a decreasing or an increasing function   |
!  |  depending on the mean PV gradient at the initial time.  This is used   |
!  |  to find the equivalent latitude y_e(q) and its inverse q_e(y).  The    |
!  |  gradient dq_e/dy, appropriately normalised (see below), is used to     |
!  |  identify the number of jets and measure the degree of "staircasing".   |
!  |                                                                         |
!  |     Output files (in diagnostics directory):                            |
!  |  jets.asc     - fraction PV variation due to the jets in each layer     |
!  |  qvar.asc     - (adjusted) PV variation / (beta*L_y) in each layer      |
!  |  ipdf1.asc    - log_10 PDFs of homogenisation, PV gradient and          |
!  |                 jet widths for layer 1                                  |
!  |  ipdf2.asc    - log_10 PDFs of homogenisation, PV gradient and          |
!  |                 jet widths for layer 2                                  |
!  |  hov1.r4      - q_e(y,t) for layer 1                                    |
!  |  hov2.r4      - q_e(y,t) for layer 2                                    |
!  |  hov-plotting - information for plotting hov1 & hov2.r4                 |
!  |  yeq.r4       - y_e(q,t) in both layers; plot using "yeqv"              |
!  |  pdf.r4       - scaled dq_e/dy vs scaled y in both layers;              |
!  |                 plot using "pdfv"                                       |
!  |                                                                         |
!  |  Completed 30 January 2015 by D G Dritschel @ St Andrews                |
!  ---------------------------------------------------------------------------

 !Import constants:
use constants

implicit none
 !Number of bins used to form homogenising and jet intensity PDFs:
integer,parameter:: ndiv=100
real,parameter:: divi=float(ndiv), div=1./divi, ds=1./float(ny)

 !Used to remove high boundary PV from the PDF p(q):
real,parameter:: cmin=0.5

real:: qq(0:ny,0:nxm1,nz), t
real:: pdfpvg(ndiv,nz), pdfds(ndiv,nz), pdfhom(ndiv,nz)
real:: pdf(ny,nz), s(ny), q(0:ny+1)
real:: yeq(0:ny+1,nz), y(0:ny), qe(0:ny)

real:: yfac, pfac, tbeg
real:: qmin, qmax, qbot, qtop, dir, dbin, dbini
real:: p, sum, ww, fhom, ps0, ps1, ps2, sbar, dels, pdfmax
real:: fpvg(nz), delq(nz), nhom(nz), npvg(nz), fnhom(nz), fnpvg(nz)

integer:: loop, iread, ix, iy, iz, k, m, nsamp

logical sample, peak

!---------------------------------------------------------------------
 !Open files containing the PV in each layer:
open(21,file='evolution/qq1.r4',form='unformatted', &
      access='direct',status='old',recl=nbytes)
open(22,file='evolution/qq2.r4',form='unformatted', &
      access='direct',status='old',recl=nbytes)
 
 !Open output files:
open(31,file='diagnostics/jets.asc',status='replace')
open(32,file='diagnostics/qvar.asc',status='replace')
open(41,file='diagnostics/yeq.r4',form='unformatted',access='direct', &
                    status='replace',recl=4*(1+4*ny))
open(42,file='diagnostics/pdf.r4',form='unformatted',access='direct', &
                    status='replace',recl=4*(1+4*ny))
open(61,file='diagnostics/hov1.r4',form='unformatted',access='stream', &
                    status='replace')
open(62,file='diagnostics/hov2.r4',form='unformatted',access='stream', &
                    status='replace')
t=0.
write(61) t
write(62) t

write(*,*) ' Time to start accumulating intensity PDFs?'
read(*,*) tbeg

 !Initialise intensity PDFs:
do iz=1,nz
  nhom(iz)=0.
  npvg(iz)=0.
  do m=1,ndiv
    pdfhom(m,iz)=0.
    pdfpvg(m,iz)=0.
    pdfds(m,iz)=0.
  enddo
enddo

 !y grid lines (used for making the hovmuller diagrams):
do iy=0,ny
  y(iy)=ymin+gly*float(iy)
enddo

 !Normalised PV values ranging from 0 to 1:
do iy=1,ny
  s(iy)=float(iy-1)/float(ny-1)
enddo

 !Used for converting PDFs into y_e(q):
yfac=garea/ellx

!---------------------------------------------------------------
 !Read data and process:
loop=0
nsamp=0
do
  loop=loop+1
  iread=0
  read(21,rec=loop,iostat=iread) t,qq(:,:,1)
  if (iread .ne. 0) exit 
  read(22,rec=loop,iostat=iread) t,qq(:,:,2)

  write(*,'(a,f13.5)') ' Processing t = ',t

   !Logical used to decide when to accumulate intensity PDFs:
  sample=(t .ge. tbeg)
  if (sample) nsamp=nsamp+1

   !Loop over layers:
  do iz=1,nz

     !Compute min & max PV:
    qmin=qq(0,0,iz)
    qmax=qq(0,0,iz)

    do ix=0,nxm1
      do iy=0,ny
        qmin=min(qmin,qq(iy,ix,iz))
        qmax=max(qmax,qq(iy,ix,iz))
      enddo
    enddo

     !Increment in PV used for binning PV:
    dbin=(qmax-qmin)/float(ny-1)
    dbini=1./dbin

     !Initialise probability:
    do k=1,ny
      pdf(k,iz)=0.
    enddo

     !Edge values at iy = 0 & ny count half:
    do ix=0,nxm1
      k=nint((qq(0 ,ix,iz)-qmin)*dbini)+1
      pdf(k,iz)=pdf(k,iz)+0.5
      k=nint((qq(ny,ix,iz)-qmin)*dbini)+1
      pdf(k,iz)=pdf(k,iz)+0.5
    enddo

     !Interior values count 1:
    do ix=0,nxm1
      do iy=1,ny-1
        k=nint((qq(iy,ix,iz)-qmin)*dbini)+1
        pdf(k,iz)=pdf(k,iz)+1.
      enddo
    enddo

     !Normalise probabilities so they integrate to 1 over a rescaled
     !q range going from 0 to 1 (without loss of generality):
    sum=0.
    do k=1,ny
      sum=sum+pdf(k,iz)
    enddo
    sum=float(ny)/sum
    do k=1,ny
      pdf(k,iz)=pdf(k,iz)*sum
    enddo

     !Remove edges where PDF < cmin:
    k=1
    do while (pdf(k,iz) .lt. cmin)
      k=k+1
    enddo
    qmin=qmin+dbin*float(k-1)

    k=ny
    do while (pdf(k,iz) .lt. cmin)
      k=k-1
    enddo
    qmax=qmax-dbin*float(ny-k)

     !For output:
    delq(iz)=(qmax-qmin)/(beta*elly)

     !Increment in PV used for binning PV:
    dbin=(qmax-qmin)/float(ny)
    dbini=1./dbin

     !Initialise probability:
    do k=1,ny
      pdf(k,iz)=0.
    enddo

     !Edge values at iy = 0 & ny count half:
    do ix=0,nxm1
      k=nint((qq(0 ,ix,iz)-qmin)*dbini)+1
      if (k*(nyp1-k) .gt. 0) pdf(k,iz)=pdf(k,iz)+0.5
      k=nint((qq(ny,ix,iz)-qmin)*dbini)+1
      if (k*(nyp1-k) .gt. 0) pdf(k,iz)=pdf(k,iz)+0.5
    enddo

     !Interior values count 1:
    do ix=0,nxm1
      do iy=1,ny-1
        k=nint((qq(iy,ix,iz)-qmin)*dbini)+1
        if (k*(nyp1-k) .gt. 0) pdf(k,iz)=pdf(k,iz)+1.
      enddo
    enddo

     !Normalise probabilities so they integrate to 1 over a rescaled
     !q range going from 0 to 1 (without loss of generality):
    sum=0.
    do k=1,ny
      sum=sum+pdf(k,iz)
    enddo
    sum=float(ny)/sum
    do k=1,ny
      pdf(k,iz)=pdf(k,iz)*sum
    enddo

    if (sample) then
       !Diagnose fractional degree of homogenisation:
      peak=.false.
      do k=1,ny
        if (peak) then
           !We're in a region where the PDF > 1
          if (pdf(k,iz) .gt. 1.) then
             !We're still in this region; accumulate PDF and keep max value:
            fhom=fhom+pdf(k,iz)*ds
            pdfmax=max(pdfmax,pdf(k,iz))
            if (k .eq. ny) then
               !Reached upper edge of domain; accumulate intensity PDF:
              m=int(fhom*divi)+1
               !Weight the peak by a function of its maximum value:
              ww=(pdfmax-1.)/pdfmax
              pdfhom(m,iz)=pdfhom(m,iz)+ww
              nhom(iz)=nhom(iz)+ww
            endif
          else
             !The PDF < 1 here; record peak:
            m=int(fhom*divi)+1
             !Weight the peak by a function of its maximum value:
            ww=(pdfmax-1.)/pdfmax
            pdfhom(m,iz)=pdfhom(m,iz)+ww
            nhom(iz)=nhom(iz)+ww
            peak=.false.
          endif
        else
           !No peak started yet
          if (pdf(k,iz) .gt. 1.) then
             !Start a new peak:
            fhom=pdf(k,iz)*ds
            pdfmax=pdf(k,iz)
            peak=.true.
          endif
        endif
      enddo
    endif

     !---------------------------------------------------------
     !Next form y_e(q):
    qbot=0.
    qtop=0.
    do ix=0,nxm1
      qbot=qbot+qq(0 ,ix,iz)
      qtop=qtop+qq(ny,ix,iz)
    enddo
    qbot=qbot/dble(nx)
    qtop=qtop/dble(nx)

     !Sign of PV gradient in layer 1:
    dir=sign(1.,qtop-qbot)

     !Scaling factor to convert PDF into y_eq:
    pfac=yfac/sum

    if (dir .gt. 0.) then
       !Positive mean gradient:
      q(0)=qmin
      yeq(0,iz)=ymin
      q(1)=qmin+0.5*dbin
      yeq(1,iz)=pdf(1,iz)*pfac+ymin
      do k=2,ny
        q(k)=q(k-1)+dbin
        yeq(k,iz)=yeq(k-1,iz)+pdf(k,iz)*pfac
      enddo
      q(ny+1)=qmax
      yeq(ny+1,iz)=ymax
       !Invert equivalent latitude to get q_e:
      qe(0)=qmin
      k=0
      do iy=1,ny-1
        do while (y(iy) .gt. yeq(k+1,iz))
          k=k+1
        enddo
        p=(y(iy)-yeq(k,iz))/(yeq(k+1,iz)-yeq(k,iz))
        qe(iy)=q(k)+p*(q(k+1)-q(k))
      enddo
      qe(ny)=qmax
      write(60+iz) qe
    else
       !Negative mean gradient:
      q(0)=qmin
      yeq(0,iz)=ymax
      q(1)=qmin+0.5*dbin
      yeq(1,iz)=ymax-pdf(1,iz)*pfac
      do k=2,ny
        q(k)=q(k-1)+dbin
        yeq(k,iz)=yeq(k-1,iz)-pdf(k,iz)*pfac
      enddo
      q(ny+1)=qmax
      yeq(ny+1,iz)=ymin
       !Invert equivalent latitude to get q_e:
      qe(0)=qmax
      k=ny+1
      do iy=1,ny-1
        do while (y(iy) .gt. yeq(k-1,iz))
          k=k-1
        enddo
        p=(y(iy)-yeq(k,iz))/(yeq(k-1,iz)-yeq(k,iz))
        qe(iy)=q(k)+p*(q(k-1)-q(k))
      enddo
      qe(ny)=qmin
      write(60+iz) qe
    endif

     !Create a new PDF proportional to dqe/dy and with unit integral
     !over s = (y-ymin)/L_y:
    do iy=1,ny
      pdf(iy,iz)=qe(iy)-qe(iy-1)
    enddo
    sum=0.
    do k=1,ny
      sum=sum+pdf(k,iz)
    enddo
    sum=float(ny)/sum
    do k=1,ny
      pdf(k,iz)=pdf(k,iz)*sum
    enddo

    if (sample) then
       !Diagnose fractional degree that the PV is made up of jumps:
      fpvg(iz)=0.
      peak=.false.
      do k=1,ny
        if (peak) then
           !We're in a region where the PDF > 1
          if (pdf(k,iz) .gt. 1.) then
             !We're still in this region; accumulate PDF and keep max value:
            ps0=ps0+pdf(k,iz)*ds
            ps1=ps1+pdf(k,iz)*s(k)*ds
            ps2=ps2+pdf(k,iz)*s(k)**2*ds
            pdfmax=max(pdfmax,pdf(k,iz))
            if (k .eq. ny) then
               !Reached upper edge of domain; accumulate intensity PDF:
              m=int(ps0*divi)+1
               !Weight the peak by a function of its maximum value:
              ww=(pdfmax-1.)/pdfmax
              pdfpvg(m,iz)=pdfpvg(m,iz)+ww
              npvg(iz)=npvg(iz)+ww
              fpvg(iz)=fpvg(iz)+ww*ps0

               !Mean location in scaled y:
              sbar=ps1/ps0
               !Width of peak in scaled y:
              dels=2.*sqrt(abs(ps2/ps0-sbar**2))
              m=int(dels*divi)+1
              pdfds(m,iz)=pdfds(m,iz)+ww
            endif
          else
             !The PDF < 1 here; record peak:
            m=int(ps0*divi)+1
             !Weight the peak by a function of its maximum value:
            ww=(pdfmax-1.)/pdfmax
            pdfpvg(m,iz)=pdfpvg(m,iz)+ww
            npvg(iz)=npvg(iz)+ww
            fpvg(iz)=fpvg(iz)+ww*ps0

             !Mean location in scaled y:
            sbar=ps1/ps0
             !Width of peak in scaled y:
            dels=2.*sqrt(abs(ps2/ps0-sbar**2))
            m=int(dels*divi)+1
            pdfds(m,iz)=pdfds(m,iz)+ww
            peak=.false.
          endif
        else
           !No peak started yet
          if (pdf(k,iz) .gt. 1.) then
             !Start a new peak:
            ps0=pdf(k,iz)*ds
            ps1=pdf(k,iz)*s(k)*ds
            ps2=pdf(k,iz)*s(k)**2*ds
            pdfmax=pdf(k,iz)
            peak=.true.
          endif
        endif
      enddo
    endif

     !End of loop over layers:
  enddo

   !Write diagnostic data for this time:
  if (sample) write(31,'(f10.2,2(1x,f9.7))') t,fpvg(1),fpvg(2)
  write(32,'(f10.2,2(1x,f9.6))')             t,delq(1),delq(2)
  write(41,rec=loop) t,s(1:ny),yeq(1:ny,1),s(1:ny),yeq(1:ny,2)
  write(42,rec=loop) t,s(1:ny),pdf(1:ny,1),s(1:ny),pdf(1:ny,2)

   !End of loop over time:
enddo

 !Close files:
close(21)
close(22)
close(31)
close(41)
close(42)
close(61)
close(62)

 !------------------------------------------------------------------
 !Finalise intensity PDFs (ensure unit integral):
do iz=1,nz
  sum=0.
  do m=1,ndiv
    sum=sum+pdfhom(m,iz)
  enddo
  sum=float(ndiv)/sum
  do m=1,ndiv
    pdfhom(m,iz)=pdfhom(m,iz)*sum
  enddo

  sum=0.
  do m=1,ndiv
    sum=sum+pdfpvg(m,iz)
  enddo
  sum=float(ndiv)/sum
  do m=1,ndiv
    pdfpvg(m,iz)=pdfpvg(m,iz)*sum
  enddo

  sum=0.
  do m=1,ndiv
    sum=sum+pdfds(m,iz)
  enddo
  sum=float(ndiv)/sum
  do m=1,ndiv
    pdfds(m,iz)=pdfds(m,iz)*sum
  enddo
enddo

 !Write intensity PDFs:
open(20,file='diagnostics/ipdf1.asc',status='replace')
do m=1,ndiv
  write(20,'(1x,f9.7,3(1x,f9.5))') div*(float(m)-0.5), &
          log10(pdfhom(m,1)+1.e-10),log10(pdfpvg(m,1)+1.e-10), &
          log10(pdfds(m,1)+1.e-10)
enddo
close(20)

open(20,file='diagnostics/ipdf2.asc',status='replace')
do m=1,ndiv
  write(20,'(1x,f9.7,3(1x,f9.5))') div*(float(m)-0.5), &
          log10(pdfhom(m,2)+1.e-10),log10(pdfpvg(m,2)+1.e-10), &
          log10(pdfds(m,2)+1.e-10)
enddo
close(20)

 !Get time mean number of homogeneous regions and jets:
do iz=1,nz
  fnhom(iz)=nhom(iz)/float(nsamp)
  fnpvg(iz)=npvg(iz)/float(nsamp)
enddo

write(*,*)
write(*,'(a,f9.5)') ' Time mean number of homogeneous regions in layer 1 = ', &
     fnhom(1)
write(*,'(a,f9.5)') '                Time mean number of jets in layer 1 = ', &
     fnpvg(1)

write(*,*)
write(*,'(a,f9.5)') ' Time mean number of homogeneous regions in layer 2 = ', &
     fnhom(2)
write(*,'(a,f9.5)') '                Time mean number of jets in layer 2 = ', &
     fnpvg(2)

write(*,*)
write(*,*) ' Intensity PDFs for layers 1 & 2 are in ipdf1.asc & ipdf2.asc'
write(*,*)
write(*,*) ' See yeqv for y_e(q,t) vs t and pdfv for p_q(y) vs t'
write(*,*)
write(*,*) ' See yeqv for y_e(q,t) vs t and pdfv for p_q(y) vs t'
write(*,*)
write(*,*) '  The fractional PV variation in the jets is in jets.asc'
write(*,*) ' while the (adjusted) PV variation / beta is in qvar.asc'
write(*,*)

loop=loop-1
open(66,file='diagnostics/hov-plotting',status='replace')
write(66,*) &
   ' To plot Hovmuller diagrams of q_e(t,y) for each layer, type:'
write(* ,*) &
   ' To plot Hovmuller diagrams of q_e(t,y) for each layer, type:'
write(66,*)
write(* ,*)
write(66,'(a,i5,1x,i4)') ' dataview diagnostics/hov1.r4 -ndim ',loop,ny+1
write(* ,'(a,i5,1x,i4)') ' dataview diagnostics/hov1.r4 -ndim ',loop,ny+1
write(66,*)
write(* ,*)
write(66,'(a,i5,1x,i4)') ' dataview diagnostics/hov2.r4 -ndim ',loop,ny+1
write(* ,'(a,i5,1x,i4)') ' dataview diagnostics/hov2.r4 -ndim ',loop,ny+1
write(66,*)
write(* ,*)
close(66)

 !End main program
end program stair
!=======================================================================
