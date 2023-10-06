program hov

! Creates hovmuller diagrams of the data created by running zonal.f90  

! Written 4 September 2021 by D G Dritschel @ St Andrews

use parameters

implicit none

double precision:: zuu(ng),zhh(ng),zek(ng),zqq(ng),zvq(ng),zqg(0:ng)
double precision:: t,dli,ombar
integer:: iread, loop, j

!---------------------------------------------------------------------
 !Open all input data files:
open(51,file='evolution/zu.asc',status='old')
open(52,file='evolution/zh.asc',status='old')
open(53,file='evolution/zk.asc',status='old')
open(54,file='evolution/zq.asc',status='old')
open(55,file='evolution/zf.asc',status='old')

open(61,file='hovu.r8',form='unformatted',access='stream',status='replace')
open(62,file='hovh.r8',form='unformatted',access='stream',status='replace')
open(63,file='hovk.r8',form='unformatted',access='stream',status='replace')
open(64,file='hovq.r8',form='unformatted',access='stream',status='replace')
open(65,file='hovf.r8',form='unformatted',access='stream',status='replace')
open(66,file='hovg.r8',form='unformatted',access='stream',status='replace')

write(61) 0.d0
write(62) 0.d0
write(63) 0.d0
write(64) 0.d0
write(65) 0.d0
write(66) 0.d0

loop=0
do
  iread=0
  read(51,*,iostat=iread) t
  if (iread .ne. 0) exit 
  read(52,*) t
  read(53,*) t
  read(54,*) t
  read(55,*) t
  do j=1,ng
    read(51,*) zuu(j)
    read(52,*) zhh(j)
    read(53,*) zek(j)
    read(54,*) zqq(j)
    read(55,*) zvq(j)
  enddo

   !Remove mean angular velocity from zuu:
  ombar=(f1112*(zuu(1)+zuu(ng))+sum(zuu(2:ngm1)))*rsumi
  zuu=zuu-ombar*clat

  write(61) zuu
  write(62) zhh
  write(63) zek
  write(64) zqq
  write(65) zvq
   !Compute PV gradient at full latitudes (use zero polar values):
  zqg(0)=0.d0
  zqg(ng)=0.d0
  dli=dble(ng)/pi
  do j=1,ng-1
    zqg(j)=dli*(zqq(j+1)-zqq(j))
  enddo
  write(66) zqg

  loop=loop+1
   !End of loop over time:
enddo

 !Close files:
close(51)
close(52)
close(53)
close(54)
close(55)

close(61)
close(62)
close(63)
close(64)
close(65)
close(66)

open(88,file='hov-plotting',status='replace')
write(88,*) ' To plot Hovmuller diagrams of u_bar(t,phi), etc, type:'
write(* ,*) ' To plot Hovmuller diagrams of u_bar(t,phi), etc, type:'
write(* ,*)
write(88,'(a,i5,1x,i4)') ' dataview hovu.r8 -ndim ',loop,ng
write(88,'(a,i5,1x,i4)') ' dataview hovh.r8 -ndim ',loop,ng
write(88,'(a,i5,1x,i4)') ' dataview hovk.r8 -ndim ',loop,ng
write(88,'(a,i5,1x,i4)') ' dataview hovq.r8 -ndim ',loop,ng
write(88,'(a,i5,1x,i4)') ' dataview hovf.r8 -ndim ',loop,ng
write(88,'(a,i5,1x,i4)') ' dataview hovg.r8 -ndim ',loop,ng+1
write(* ,'(a,i5,1x,i4)') ' dataview hovu.r8 -ndim ',loop,ng
write(* ,'(a,i5,1x,i4)') ' dataview hovh.r8 -ndim ',loop,ng
write(* ,'(a,i5,1x,i4)') ' dataview hovk.r8 -ndim ',loop,ng
write(* ,'(a,i5,1x,i4)') ' dataview hovq.r8 -ndim ',loop,ng
write(* ,'(a,i5,1x,i4)') ' dataview hovf.r8 -ndim ',loop,ng
write(* ,'(a,i5,1x,i4)') ' dataview hovg.r8 -ndim ',loop,ng+1
write(88,*)
write(* ,*)
close(88)

 !End main program
end program hov
