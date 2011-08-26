attributes(global) subroutine convol_part1_kernel(ekinvec,n1,n2,n3,s)
       real*8,device :: ekinvec(0:(n1+1)*(n2+1)*(n3+1)-1)
       integer, value :: n1,n2,n3,s
       integer :: i, j, kb, k, tx, ty, tz

! Start execution, first get thread indices

       tx = threadidx%x-1

! thread indices over grid

       i = (blockidx%x-1) * blockDim%x + tx

if (i .le. (n1+1)*(n2+1)*(n3+1)-1-s) then
ekinvec(s+i) = 0.d0
end if

    end subroutine convol_part1_kernel

   attributes(global) subroutine convol_part2_kernel(psig,y,ekinvec,filter, ibyz_c,n1,n2,n3,lowfil,lupfil)
       real*8,device :: psig(0:n1,1:2,0:n2,1:2,0:n3,1:2),y(0:n1,1:2,0:n2,1:2,0:n3,1:2),ekinvec(0:n1,0:n2,0:n3),filter(lowfil:lupfil,4)
       integer*4,device :: ibyz_c(2,0:n2,0:n3)
       integer*4, value :: n1,n2,n3,lowfil,lupfil,l,llow,lhigh
       real*8 :: add111,ekin,t211,add112
       integer*4 :: i, j, kb, k, bx, tx, ty, tz, ka

! submatrices are declared to be in CUDA shared memory

       integer*4, shared :: ibyz_csub(0:1,0:3,1:2),lupfilsub
       real*8, shared :: filtersub(-14:14,1:4),psigsub(0:46,1:2,0:1,1:1,0:3,1:1)

! Start execution, first get thread indices

       tx = threadidx%x-1
       ty = threadidx%y-1
       tz = threadidx%z-1

! thread indices over block

       bx = tx+(ty*blockDim%x)+(tz*blockDim%x*blockDim%y)

! thread indices over grid

       i = (blockidx%x-1) * blockDim%x + tx
       j = (blockidx%y-1) * blockDim%y + ty

!Filter coefficients into shared memory 29*4*8 B=0.464 KB
if (-bx .le. 0 .and. -bx .ge. lowfil .and. bx .le. lupfil) then
filtersub(-bx,1) = filter(-bx,1)
filtersub(-bx,2) = filter(-bx,2)
filtersub(bx,1) = filter(bx,1)
filtersub(bx,2) = filter(bx,2)
endif

! (1/2) d^2/dx^2
if (i .le. n2 .and. j .le. n3) then
!ibyz_c into shared memory 4*4*2*4 B= 0.128 KB
!To avoid Non-stride-1 accesses we change the arrangement of the array.
ibyz_csub(tx,ty,1) = ibyz_c(1,i,j)
ibyz_csub(tx,ty,2) = ibyz_c(2,i,j)
call syncthreads()
do kb=ibyz_csub(tx,ty,1),ibyz_csub(tx,ty,2),32
   add111=0.d0
!part of psig into shared memory (32+15)*2*4*2*8 B=6.016 KB
k = ibyz_csub(tx,ty,2)-kb
if (tz .le. 13 .and. tz .le. k .and. tz+kb-14 .ge. 0) then
psigsub(tz,1,tx,1,ty,1) = psig(kb+tz-14,1,i,1,j,1)
psigsub(tz,2,tx,1,ty,1) = psig(kb+tz-14,2,i,1,j,1)
end if
call syncthreads()
if (k .le. 13 .and. tz .le. 13) then
psigsub(tz+k+1,1,tx,1,ty,1) = psig(kb+tz-14+k+1,1,i,1,j,1)
psigsub(tz+k+1,2,tx,1,ty,1) = psig(kb+tz-14+k+1,2,i,1,j,1)
else
psigsub(tz+14,1,tx,1,ty,1) = psig(kb+tz,1,i,1,j,1)
psigsub(tz+14,2,tx,1,ty,1) = psig(kb+tz,2,i,1,j,1)
end if
call syncthreads()

!we only use the llow do loop when needed

if (tz .le. k) then
llow = -(tz+kb)+14

if (-kb+14 .ge. 1) then
 if (llow .ge. 1) then
 do l=llow,13
 add111=add111+psigsub(tz+l,1,tx,1,ty,1)*filtersub(-14+l,1)+psigsub(tz+l,2,tx,1,ty,1)*filtersub(-14+l,2)
 end do
 else
 do l=0,13
  add111=add111+psigsub(tz+l,1,tx,1,ty,1)*filtersub(-14+l,1)+psigsub(tz+l,2,tx,1,ty,1)*filtersub(-14+l,2)
  enddo
 end if
else
  do l=0,13
  add111=add111+psigsub(tz+l,1,tx,1,ty,1)*filtersub(-14+l,1)+psigsub(tz+l,2,tx,1,ty,1)*filtersub(-14+l,2)
  enddo
end if

end if

! Copy data from the end of the array to the beginning
call syncthreads()
if (tz .le. 13 .and. tz .le. k) then
psigsub(tz,1,tx,1,ty,1)=psigsub(tz+14,1,tx,1,ty,1)
psigsub(tz,2,tx,1,ty,1)=psigsub(tz+14,2,tx,1,ty,1)
end if
if (tz .le. 13 .and. k .ge. 14) then
psigsub(tz+14,1,tx,1,ty,1)=psigsub(tz+28,1,tx,1,ty,1)
psigsub(tz+14,2,tx,1,ty,1)=psigsub(tz+28,2,tx,1,ty,1)
end if
call syncthreads()
if (tz .ge. 28 .and. k .ge. 28) then 
psigsub(tz,1,tx,1,ty,1)=psigsub(tz+14,1,tx,1,ty,1) 
psigsub(tz,2,tx,1,ty,1)=psigsub(tz+14,2,tx,1,ty,1) 
end if 
call syncthreads()


! Fill array with new data
if (tz .le. 13 .and. tz .le. k-1) then
psigsub(tz+32,1,tx,1,ty,1)=psig(tz+kb+32,1,i,1,j,1)
psigsub(tz+32,2,tx,1,ty,1)=psig(tz+kb+32,2,i,1,j,1)
end if
call syncthreads()
if (tz .le. 13 .and. k .le. 31) then
psigsub(tz+k+1,1,tx,1,ty,1)=psig(tz+kb+k+1,1,i,1,j,1)
psigsub(tz+k+1,2,tx,1,ty,1)=psig(tz+kb+k+1,2,i,1,j,1)
end if
call syncthreads()

lupfilsub = 14
!we only use the lhigh do loop when needed
if (tz .le. k) then
lhigh = n1-(tz+kb)
 
if (n1-(kb+blockDim%z) .le. lupfilsub-1) then
 if (lhigh .le. lupfilsub-1) then
 do l=0,lhigh
 add111=add111+psigsub(tz+l,1,tx,1,ty,1)*filtersub(l,1)+psigsub(tz+l,2,tx,1,ty,1)*filtersub(l,2)
 end do
 else
 do l=0,lupfilsub
 add111=add111+psigsub(tz+l,1,tx,1,ty,1)*filtersub(l,1)+psigsub(tz+l,2,tx,1,ty,1)*filtersub(l,2)
 end do
 end if
else
 do l=0,lupfilsub
 add111=add111+psigsub(tz+l,1,tx,1,ty,1)*filtersub(l,1)+psigsub(tz+l,2,tx,1,ty,1)*filtersub(l,2)
 end do 
end if

call syncthreads()
y(kb+tz,1,i,1,j,1)=y(kb+tz,1,i,1,j,1)+add111
ekinvec(kb+tz,i,j)=add111*psigsub(tz,1,tx,1,ty,1)
end if

   enddo
end if

    end subroutine convol_part2_kernel


        attributes(global) subroutine convol_part3_kernel(psig,y,ekinvec,filter,ibyz_f,n1,n2,n3,lowfil,lupfil)
       real*8,device :: psig(0:n1,1:2,0:n2,1:2,0:n3,1:2),y(0:n1,1:2,0:n2,1:2,0:n3,1:2),ekinvec(0:n1,0:n2,0:n3),filter(lowfil:lupfil,4)
       integer*4,device :: ibyz_f(2,0:n2,0:n3)
       integer*4, value :: n1,n2,n3,lowfil,lupfil,l,llow,lhigh
       real*8 :: add111,ekin,t211,add112
       integer*4 :: i, j, kb, k, bx, tx, ty, tz, ka

! submatrices are declared to be in CUDA shared memory

       integer*4, shared :: ibyz_fsub(0:3,0:0,1:2),lupfilsub
       real*8, shared :: filtersub(-14:14,1:4),psigsub(0:60,1:1,0:3,1:1,0:0,1:1)

! Start execution, first get thread indices

       tx = threadidx%x-1
       ty = threadidx%y-1
       tz = threadidx%z-1

! thread indices over block

       bx = tx+(ty*blockDim%x)+(tz*blockDim%x*blockDim%y)

! thread indices over grid

       i = (blockidx%x-1) * blockDim%x + tx
       j = (blockidx%y-1) * blockDim%y + ty

!Filter coefficients into shared memory 29*8 B=0.232 KB
if (-bx .le. 0 .and. -bx .ge. lowfil .and. bx .le. lupfil) then
filtersub(-bx,3) = filter(-bx,3)
filtersub(bx,3) = filter(bx,3)
endif
call syncthreads()

! (1/2) d^2/dx^2
if (i .le. n2 .and. j .le. n3) then
!ibyz_f into shared memory 4*4*2*4 B= 0.128 KB
!To avoid Non-stride-1 accesses we change the arrangement of the array.
ibyz_fsub(tx,ty,1) = ibyz_f(1,i,j)
ibyz_fsub(tx,ty,2) = ibyz_f(2,i,j)
call syncthreads()
do kb=ibyz_fsub(tx,ty,1),ibyz_fsub(tx,ty,2),32
   add111=0.d0
!part of psig into shared memory (32+29)*4*1*8 B=1.920 KB
k = ibyz_fsub(tx,ty,2)-kb
if (tz .le. 28 .and. tz .le. k .and. tz+kb-14 .ge. 0) then
psigsub(tz,1,tx,1,ty,1) = psig(kb+tz-14,1,i,1,j,1)
end if
call syncthreads()
if (k .le. 28 .and. tz .le. 28) then
psigsub(tz+k+1,1,tx,1,ty,1) = psig(kb+tz-14+k+1,1,i,1,j,1)
else
psigsub(tz+29,1,tx,1,ty,1) = psig(kb+tz+15,1,i,1,j,1)
end if
call syncthreads()

!we only use the llow do loop when needed

if (tz .le. k) then
llow = -(tz+kb)+14

if (-kb+14 .ge. 1) then
 if (llow .ge. 1) then
 do l=llow,13
 add111=add111+psigsub(tz+l,1,tx,1,ty,1)*filtersub(-14+l,3)
 end do
 else
 do l=0,13
  add111=add111+psigsub(tz+l,1,tx,1,ty,1)*filtersub(-14+l,3)
  enddo
 end if
else
  do l=0,13
  add111=add111+psigsub(tz+l,1,tx,1,ty,1)*filtersub(-14+l,3)
  enddo
end if

end if

lupfilsub = 14
!we only use the lhigh do loop when needed
if (tz .le. k) then
lhigh = n1-(tz+kb)
 
if (n1-(kb+blockDim%z) .le. lupfilsub-1) then
 if (lhigh .le. lupfilsub-1) then
 do l=0,lhigh
 add111=add111+psigsub(tz+l+14,1,tx,1,ty,1)*filtersub(l,3)
 end do
 else
 do l=0,lupfilsub
 add111=add111+psigsub(tz+l+14,1,tx,1,ty,1)*filtersub(l,3)
 end do
 end if
else
 do l=0,lupfilsub
 add111=add111+psigsub(tz+l+14,1,tx,1,ty,1)*filtersub(l,3)
 end do 
end if

y(kb+tz,2,i,1,j,1)=y(kb+tz,2,i,1,j,1)+add111
ekinvec(kb+tz,i,j)=ekinvec(kb+tz,i,j)+add111*psig(kb+tz,2,i,1,j,1)
call syncthreads()
end if

   enddo
end if

    end subroutine convol_part3_kernel


        attributes(global) subroutine convol_part4_kernel(psig,y,ekinvec,filter,ibxz_c,n1,n2,n3,lowfil,lupfil)
       real*8,device :: psig(0:n1,1:2,0:n2,1:2,0:n3,1:2),y(0:n1,1:2,0:n2,1:2,0:n3,1:2),ekinvec(0:n1,0:n2,0:n3),filter(lowfil:lupfil,4)
       integer*4,device :: ibxz_c(2,0:n1,0:n3)
       integer*4, value :: n1,n2,n3,lowfil,lupfil,l,llow,lhigh
       real*8 :: add111,ekin,t211,add112
       integer*4 :: i, j, kb, k, bx, tx, ty, tz, ka

! submatrices are declared to be in CUDA shared memory

       integer*4, shared :: ibxz_csub(0:3,0:1,1:2),lupfilsub
       real*8, shared :: filtersub(-14:14,1:4),psigsub(0:3,1:1,0:46,1:2,0:1,1:1)

! Start execution, first get thread indices

       tx = threadidx%x-1
       ty = threadidx%y-1
       tz = threadidx%z-1

! thread indices over block

       bx = tx+(ty*blockDim%x)+(tz*blockDim%x*blockDim%y)

! thread indices over grid

       i = (blockidx%x-1) * blockDim%x + tx
       j = (blockidx%y-1) * blockDim%y + ty

!Filter coefficients into shared memory 29*4*8 B=0.464 KB
if (-bx .le. 0 .and. -bx .ge. lowfil .and. bx .le. lupfil) then
filtersub(-bx,1) = filter(-bx,1)
filtersub(-bx,2) = filter(-bx,2)
filtersub(bx,1) = filter(bx,1)
filtersub(bx,2) = filter(bx,2)
endif

! (1/2) d^2/dy^2
if (i .le. n1 .and. j .le. n3) then
!ibxz_c into shared memory 4*4*2*4 B= 0.128 KB
!To avoid Non-stride-1 accesses we change the arrangement of the array.
ibxz_csub(tx,ty,1) = ibxz_c(1,i,j)
ibxz_csub(tx,ty,2) = ibxz_c(2,i,j)
call syncthreads()
do kb=ibxz_csub(tx,ty,1),ibxz_csub(tx,ty,2),32
   add111=0.d0
!part of psig into shared memory (32+15)*4*2*2*8 B=6.016 KB
k = ibxz_csub(tx,ty,2)-kb
if (tz .le. 13 .and. tz .le. k .and. tz+kb-14 .ge. 0) then
psigsub(tx,1,tz,1,ty,1) = psig(i,1,kb+tz-14,1,j,1)
psigsub(tx,1,tz,2,ty,1) = psig(i,1,kb+tz-14,2,j,1)
end if
call syncthreads()
if (k .le. 13 .and. tz .le. 13) then
psigsub(tx,1,tz+k+1,1,ty,1) = psig(i,1,kb+tz-14+k+1,1,j,1)
psigsub(tx,1,tz+k+1,2,ty,1) = psig(i,1,kb+tz-14+k+1,2,j,1)
else
psigsub(tx,1,tz+14,1,ty,1) = psig(i,1,kb+tz,1,j,1)
psigsub(tx,1,tz+14,2,ty,1) = psig(i,1,kb+tz,2,j,1)
end if
call syncthreads()

!we only use the llow do loop when needed

if (tz .le. k) then
llow = -(tz+kb)+14

if (-kb+14 .ge. 1) then
 if (llow .ge. 1) then
 do l=llow,13
 add111=add111+psigsub(tx,1,tz+l,1,ty,1)*filtersub(-14+l,1)+psigsub(tx,1,tz+l,2,ty,1)*filtersub(-14+l,2)
 end do
 else
 do l=0,13
  add111=add111+psigsub(tx,1,tz+l,1,ty,1)*filtersub(-14+l,1)+psigsub(tx,1,tz+l,2,ty,1)*filtersub(-14+l,2)
  enddo
 end if
else
  do l=0,13
  add111=add111+psigsub(tx,1,tz+l,1,ty,1)*filtersub(-14+l,1)+psigsub(tx,1,tz+l,2,ty,1)*filtersub(-14+l,2)
  enddo
end if

end if

! Copy data from the end of the array to the beginning
call syncthreads()
if (tz .le. 13 .and. tz .le. k) then
psigsub(tx,1,tz,1,ty,1)=psigsub(tx,1,tz+14,1,ty,1)
psigsub(tx,1,tz,2,ty,1)=psigsub(tx,1,tz+14,2,ty,1)
end if
if (tz .le. 13 .and. k .ge. 14) then
psigsub(tx,1,tz+14,1,ty,1)=psigsub(tx,1,tz+28,1,ty,1)
psigsub(tx,1,tz+14,2,ty,1)=psigsub(tx,1,tz+28,2,ty,1)
end if
call syncthreads()
if (tz .ge. 28 .and. k .ge. 28) then 
psigsub(tx,1,tz,1,ty,1)=psigsub(tx,1,tz+14,1,ty,1) 
psigsub(tx,1,tz,2,ty,1)=psigsub(tx,1,tz+14,2,ty,1) 
end if 
call syncthreads()

! Fill array with new data
if (tz .le. 13 .and. tz .le. k-1) then
psigsub(tx,1,tz+32,1,ty,1)=psig(i,1,tz+kb+32,1,j,1)
psigsub(tx,1,tz+32,2,ty,1)=psig(i,1,tz+kb+32,2,j,1)
end if
call syncthreads()
if (tz .le. 13 .and. k .le. 31) then
psigsub(tx,1,tz+k+1,1,ty,1)=psig(i,1,tz+kb+k+1,1,j,1)
psigsub(tx,1,tz+k+1,2,ty,1)=psig(i,1,tz+kb+k+1,2,j,1)
end if
call syncthreads()

!instead of the min. in the cpu version we use if clauses, this is faster as:

lupfilsub = 14
!we only use the lhigh do loop when needed
if (tz .le. k) then
lhigh = n2-(tz+kb)
 
if (n2-(kb+blockDim%z) .le. lupfilsub-1) then
 if (lhigh .le. lupfilsub-1) then
 do l=0,lhigh
 add111=add111+psigsub(tx,1,tz+l,1,ty,1)*filtersub(l,1)+psigsub(tx,1,tz+l,2,ty,1)*filtersub(l,2)
 end do
 else
 do l=0,lupfilsub
 add111=add111+psigsub(tx,1,tz+l,1,ty,1)*filtersub(l,1)+psigsub(tx,1,tz+l,2,ty,1)*filtersub(l,2)
 end do
 end if
else
 do l=0,lupfilsub
 add111=add111+psigsub(tx,1,tz+l,1,ty,1)*filtersub(l,1)+psigsub(tx,1,tz+l,2,ty,1)*filtersub(l,2)
 end do 
end if

y(i,1,kb+tz,1,j,1)=y(i,1,kb+tz,1,j,1)+add111
ekinvec(i,kb+tz,j)=ekinvec(i,kb+tz,j)+add111*psigsub(tx,1,tz,1,ty,1)
call syncthreads()
end if

   enddo
end if

    end subroutine convol_part4_kernel


 attributes(global) subroutine convol_part5_kernel(psig,y,ekinvec,filter,ibxz_f,n1,n2,n3,lowfil,lupfil)
       real*8,device :: psig(0:n1,1:2,0:n2,1:2,0:n3,1:2),y(0:n1,1:2,0:n2,1:2,0:n3,1:2),ekinvec(0:n1,0:n2,0:n3),filter(lowfil:lupfil,4)
       integer*4,device :: ibxz_f(2,0:n1,0:n3)
       integer*4, value :: n1,n2,n3,lowfil,lupfil,l,llow,lhigh
       real*8 :: add111,ekin,t211,add112
       integer*4 :: i, j, kb, k, bx, tx, ty, tz, ka

! submatrices are declared to be in CUDA shared memory

       integer*4, shared :: ibxz_fsub(0:3,0:0,1:2),lupfilsub
       real*8, shared :: filtersub(-14:14,1:4),psigsub(0:3,1:1,0:60,1:1,0:0,1:1)

! Start execution, first get thread indices

       tx = threadidx%x-1
       ty = threadidx%y-1
       tz = threadidx%z-1

! thread indices over block

       bx = tx+(ty*blockDim%x)+(tz*blockDim%x*blockDim%y)

! thread indices over grid

       i = (blockidx%x-1) * blockDim%x + tx
       j = (blockidx%y-1) * blockDim%y + ty

!Filter coefficients into shared memory 29*8 B=0.232 KB
if (-bx .le. 0 .and. -bx .ge. lowfil .and. bx .le. lupfil) then
filtersub(-bx,3) = filter(-bx,3)
filtersub(bx,3) = filter(bx,3)
endif
call syncthreads()

! (1/2) d^2/dy^2
if (i .le. n1 .and. j .le. n3) then
!ibyz_f into shared memory 4*4*2*4 B= 0.128 KB
!To avoid Non-stride-1 accesses we change the arrangement of the array.
ibxz_fsub(tx,ty,1) = ibxz_f(1,i,j)
ibxz_fsub(tx,ty,2) = ibxz_f(2,i,j)
call syncthreads()
do kb=ibxz_fsub(tx,ty,1),ibxz_fsub(tx,ty,2),32
   add111=0.d0
!part of psig into shared memory (32+29)*4*1*8 B=1.920 KB
k = ibxz_fsub(tx,ty,2)-kb

 if (tz .le. 28 .and. tz .le. k .and. tz+kb-14 .ge. 0) then
 psigsub(tx,1,tz,1,ty,1) = psig(i,1,kb+tz-14,1,j,1)
 end if
 if (k .le. 28 .and. tz .le. 28) then
 psigsub(tx,1,tz+k+1,1,ty,1) = psig(i,1,kb+tz-14+k+1,1,j,1)
 else
 psigsub(tx,1,tz+29,1,ty,1) = psig(i,1,kb+tz+15,1,j,1)
 end if
 call syncthreads()

!we only use the llow do loop when needed

if (tz .le. k) then
llow = -(tz+kb)+14

if (-kb+14 .ge. 1) then
 if (llow .ge. 1) then
 do l=llow,13
 add111=add111+psigsub(tx,1,tz+l,1,ty,1)*filtersub(-14+l,3)
 end do
 else
 do l=0,13
  add111=add111+psigsub(tx,1,tz+l,1,ty,1)*filtersub(-14+l,3)
  enddo
 end if
else
  do l=0,13
  add111=add111+psigsub(tx,1,tz+l,1,ty,1)*filtersub(-14+l,3)
  enddo
end if

end if


lupfilsub = 14
!we only use the lhigh do loop when needed
if (tz .le. k) then
lhigh = n2-(tz+kb)
 
 if (lhigh .le. lupfilsub-1) then
 do l=0,lhigh
 add111=add111+psigsub(tx,1,tz+l+14,1,ty,1)*filtersub(l,3)
 end do
 else
 do l=0,lupfilsub
 add111=add111+psigsub(tx,1,tz+l+14,1,ty,1)*filtersub(l,3)
 end do
 end if

ekinvec(i,kb+tz,j)=ekinvec(i,kb+tz,j)+add111*psig(i,1,kb+tz,2,j,1)
y(i,1,kb+tz,2,j,1)=y(i,1,kb+tz,2,j,1)+add111
call syncthreads()
end if

   enddo
end if

    end subroutine convol_part5_kernel

        attributes(global) subroutine convol_part6_kernel(psig,y,ekinvec,filter,ibxy_c,n1,n2,n3,lowfil,lupfil)
       real*8,device :: psig(0:n1,1:2,0:n2,1:2,0:n3,1:2),y(0:n1,1:2,0:n2,1:2,0:n3,1:2),ekinvec(0:n1,0:n2,0:n3),filter(lowfil:lupfil,4)
       integer*4,device :: ibxy_c(2,0:n1,0:n2)
       integer*4, value :: n1,n2,n3,lowfil,lupfil,l,llow,lhigh
       real*8 :: add111,ekin,t211,add112
       integer*4 :: i, j, kb, k, bx, tx, ty, tz, ka

! submatrices are declared to be in CUDA shared memory

       integer*4, shared :: ibxy_csub(0:7,0:0,1:2),lupfilsub
       real*8, shared :: filtersub(-14:14,1:4),psigsub(0:7,1:1,0:0,1:1,0:46,1:2)

! Start execution, first get thread indices

       tx = threadidx%x-1
       ty = threadidx%y-1
       tz = threadidx%z-1

! thread indices over block

       bx = tx+(ty*blockDim%x)+(tz*blockDim%x*blockDim%y)

! thread indices over grid

       i = (blockidx%x-1) * blockDim%x + tx
       j = (blockidx%y-1) * blockDim%y + ty

!Filter coefficients into shared memory 29*4*8 B=0.464 KB
if (-bx .le. 0 .and. -bx .ge. lowfil .and. bx .le. lupfil) then
filtersub(-bx,1) = filter(-bx,1)
filtersub(-bx,2) = filter(-bx,2)
filtersub(bx,1) = filter(bx,1)
filtersub(bx,2) = filter(bx,2)
endif

! (1/2) d^2/dz^2

if (i .le. n1 .and. j .le. n2) then
!ibyz_c into shared memory 4*4*2*4 B= 0.128 KB
!To avoid Non-stride-1 accesses we change the arrangement of the array.
ibxy_csub(tx,ty,1) = ibxy_c(1,i,j)
ibxy_csub(tx,ty,2) = ibxy_c(2,i,j)
call syncthreads()
do kb=ibxy_csub(tx,ty,1),ibxy_csub(tx,ty,2),32
   add111=0.d0
!part of psig into shared memory (32+15)*4*1*2*8 B=3.008 KB
k = ibxy_csub(tx,ty,2)-kb

if (tz .le. 13 .and. tz .le. k .and. tz+kb-14 .ge. 0) then
psigsub(tx,1,ty,1,tz,1) = psig(i,1,j,1,kb+tz-14,1)
psigsub(tx,1,ty,1,tz,2) = psig(i,1,j,1,kb+tz-14,2)
end if
call syncthreads()
if (k .le. 13 .and. tz .le. 13) then
psigsub(tx,1,ty,1,tz+k+1,1) = psig(i,1,j,1,kb+tz-14+k+1,1)
psigsub(tx,1,ty,1,tz+k+1,2) = psig(i,1,j,1,kb+tz-14+k+1,2)
else
psigsub(tx,1,ty,1,tz+14,1) = psig(i,1,j,1,kb+tz,1)
psigsub(tx,1,ty,1,tz+14,2) = psig(i,1,j,1,kb+tz,2)
end if
call syncthreads()

!we only use the llow do loop when needed

if (tz .le. k) then
llow = -(tz+kb)+14

if (-kb+14 .ge. 1) then
 if (llow .ge. 1) then
 do l=llow,13
 add111=add111+psigsub(tx,1,ty,1,tz+l,1)*filtersub(-14+l,1)+psigsub(tx,1,ty,1,tz+l,2)*filtersub(-14+l,2)
 end do
 else
 do l=0,13
  add111=add111+psigsub(tx,1,ty,1,tz+l,1)*filtersub(-14+l,1)+psigsub(tx,1,ty,1,tz+l,2)*filtersub(-14+l,2)
  enddo
 end if
else
  do l=0,13
  add111=add111+psigsub(tx,1,ty,1,tz+l,1)*filtersub(-14+l,1)+psigsub(tx,1,ty,1,tz+l,2)*filtersub(-14+l,2)
  enddo
end if

end if

! Copy data from the end of the array to the beginning
call syncthreads()
if (tz .le. 13 .and. tz .le. k) then
psigsub(tx,1,ty,1,tz,1)=psigsub(tx,1,ty,1,tz+14,1)
psigsub(tx,1,ty,1,tz,2)=psigsub(tx,1,ty,1,tz+14,2)
end if
if (tz .le. 13 .and. k .ge. 14) then
psigsub(tx,1,ty,1,tz+14,1)=psigsub(tx,1,ty,1,tz+28,1)
psigsub(tx,1,ty,1,tz+14,2)=psigsub(tx,1,ty,1,tz+28,2)
end if
call syncthreads()
if (tz .ge. 28 .and. k .ge. 28) then 
psigsub(tx,1,ty,1,tz,1)=psigsub(tx,1,ty,1,tz+14,1) 
psigsub(tx,1,ty,1,tz,2)=psigsub(tx,1,ty,1,tz+14,2) 
end if 
call syncthreads()

! Fill array with new data
if (tz .le. 13 .and. tz .le. k-1) then
psigsub(tx,1,ty,1,tz+32,1)=psig(i,1,j,1,tz+kb+32,1)
psigsub(tx,1,ty,1,tz+32,2)=psig(i,1,j,1,tz+kb+32,2)
end if
call syncthreads()
if (tz .le. 13 .and. k .le. 31) then
psigsub(tx,1,ty,1,tz+k+1,1)=psig(i,1,j,1,tz+kb+k+1,1)
psigsub(tx,1,ty,1,tz+k+1,2)=psig(i,1,j,1,tz+kb+k+1,2)
end if
call syncthreads()


!instead of the min. in the cpu version we use if clauses, this is faster as:

lupfilsub = 14
!we only use the lhigh do loop when needed
if (tz .le. k) then
lhigh = n3-(tz+kb)
 
if (n3-(kb+blockDim%z) .le. lupfilsub-1) then
 if (lhigh .le. lupfilsub-1) then
 do l=0,lhigh
 add111=add111+psigsub(tx,1,ty,1,tz+l,1)*filtersub(l,1)+psigsub(tx,1,ty,1,tz+l,2)*filtersub(l,2)
 end do
 else
 do l=0,lupfilsub
 add111=add111+psigsub(tx,1,ty,1,tz+l,1)*filtersub(l,1)+psigsub(tx,1,ty,1,tz+l,2)*filtersub(l,2)
 end do
 end if
else
 do l=0,lupfilsub
 add111=add111+psigsub(tx,1,ty,1,tz+l,1)*filtersub(l,1)+psigsub(tx,1,ty,1,tz+l,2)*filtersub(l,2)
 end do 
end if

ekinvec(i,j,kb+tz)=ekinvec(i,j,kb+tz)+add111*psigsub(tx,1,ty,1,tz,1)
y(i,1,j,1,kb+tz,1)=y(i,1,j,1,kb+tz,1)+add111
call syncthreads()
end if

   enddo
end if

    end subroutine convol_part6_kernel

attributes(global) subroutine convol_part7_kernel(psig,y,ekinvec,filter, ibxy_f,n1,n2,n3,lowfil,lupfil)
       real*8,device :: psig(0:n1,1:2,0:n2,1:2,0:n3,1:2),y(0:n1,1:2,0:n2,1:2,0:n3,1:2),ekinvec(0:n1,0:n2,0:n3),filter(lowfil:lupfil,4)
       integer*4,device :: ibxy_f(2,0:n1,0:n2)
       integer*4, value :: n1,n2,n3,lowfil,lupfil,l,llow,lhigh
       real*8 :: add111,ekin,t211,add112
       integer*4 :: i, j, kb, k, bx, tx, ty, tz, ka

! submatrices are declared to be in CUDA shared memory

       integer*4, shared :: ibxy_fsub(0:3,0:0,1:2),lupfilsub
       real*8, shared :: filtersub(-14:14,1:4),psigsub(0:3,1:1,0:0,1:1,0:60,1:1)

! Start execution, first get thread indices

       tx = threadidx%x-1
       ty = threadidx%y-1
       tz = threadidx%z-1

! thread indices over block

       bx = tx+(ty*blockDim%x)+(tz*blockDim%x*blockDim%y)

! thread indices over grid

       i = (blockidx%x-1) * blockDim%x + tx
       j = (blockidx%y-1) * blockDim%y + ty

!Filter coefficients into shared memory 29*8 B=0.232 KB
if (-bx .le. 0 .and. -bx .ge. lowfil .and. bx .le. lupfil) then
filtersub(-bx,3) = filter(-bx,3)
filtersub(bx,3) = filter(bx,3)
endif
call syncthreads()

! (1/2) d^2/dz^2
if (i .le. n1 .and. j .le. n2) then
!ibyz_f into shared memory 4*4*2*4 B= 0.128 KB
!To avoid Non-stride-1 accesses we change the arrangement of the array.
ibxy_fsub(tx,ty,1) = ibxy_f(1,i,j)
ibxy_fsub(tx,ty,2) = ibxy_f(2,i,j)
call syncthreads()
do kb=ibxy_fsub(tx,ty,1),ibxy_fsub(tx,ty,2),32
   add111=0.d0
!part of psig into shared memory (32+29)*4*1*8 B=1.920 KB
k = ibxy_fsub(tx,ty,2)-kb
if (tz .le. 28 .and. tz .le. k .and. tz+kb-14 .ge. 0) then
psigsub(tx,1,ty,1,tz,1) = psig(i,1,j,1,kb+tz-14,1)
end if
call syncthreads()
if (k .le. 28 .and. tz .le. 28) then
psigsub(tx,1,ty,1,tz+k+1,1) = psig(i,1,j,1,kb+tz-14+k+1,1)
else
psigsub(tx,1,ty,1,tz+29,1) = psig(i,1,j,1,kb+tz+15,1)
end if
call syncthreads()

!we only use the llow do loop when needed

if (tz .le. k) then
llow = -(tz+kb)+14

if (-kb+14 .ge. 1) then
 if (llow .ge. 1) then
 do l=llow,13
 add111=add111+psigsub(tx,1,ty,1,tz+l,1)*filtersub(-14+l,3)
 end do
 else
 do l=0,13
  add111=add111+psigsub(tx,1,ty,1,tz+l,1)*filtersub(-14+l,3)
  enddo
 end if
else
  do l=0,13
  add111=add111+psigsub(tx,1,ty,1,tz+l,1)*filtersub(-14+l,3)
  enddo
end if

end if

lupfilsub = 14
!we only use the lhigh do loop when needed
if (tz .le. k) then
lhigh = n3-(tz+kb)
 
if (n3-(kb+blockDim%z) .le. lupfilsub-1) then
 if (lhigh .le. lupfilsub-1) then
 do l=0,lhigh
 add111=add111+psigsub(tx,1,ty,1,tz+l+14,1)*filtersub(l,3)
 end do
 else
 do l=0,lupfilsub
 add111=add111+psigsub(tx,1,ty,1,tz+l+14,1)*filtersub(l,3)
 end do
 end if
else
 do l=0,lupfilsub
 add111=add111+psigsub(tx,1,ty,1,tz+l+14,1)*filtersub(l,3)
 end do 
end if

ekinvec(i,j,kb+tz)=ekinvec(i,j,kb+tz)+add111*psig(i,1,j,1,kb+tz,2)
y(i,1,j,1,kb+tz,2)=y(i,1,j,1,kb+tz,2)+add111
call syncthreads()
end if

   enddo
end if

    end subroutine convol_part7_kernel


attributes(global) subroutine convol_part8_kernel(psig,y,ekinvec,filter,ibyz_f,nfl1,nfl2,nfl3,nfu1,nfu2,nfu3,n1,n2,n3,lowfil,lupfil)
       real*8,device :: psig(0:n1,1:2,0:n2,1:2,0:n3,1:2),y(0:n1,1:2,0:n2,1:2,0:n3,1:2),ekinvec(0:n1,0:n2,0:n3),filter(lowfil:lupfil,4)
       integer*4,device :: ibyz_f(2,0:n2,0:n3)
       integer*4, value :: nfl1,nfl2,nfl3,nfu1,nfu2,nfu3,n1,n2,n3,lowfil,lupfil,l,llow,lhigh
       real*8 :: ekin,t111,t112,t121,t211,t122,t212,t221,t222
       integer*4 :: i, j, kb, k, bx, tx, ty, tz, ka

! submatrices are declared to be in CUDA shared memory

       integer*4, shared :: ibyz_fsub(0:0,0:7,1:2),lupfilsub
       real*8, shared :: filtersub(-14:14,1:4),psigsub(0:46,1:2,0:0,1:2,0:7,2:2)

! Start execution, first get thread indices

       tx = threadidx%x-1
       ty = threadidx%y-1
       tz = threadidx%z-1

! thread indices over block

       bx = tx+(ty*blockDim%x)+(tz*blockDim%x*blockDim%y)

! thread indices over grid

       i = (blockidx%x-1) * blockDim%x + tx
       j = (blockidx%y-1) * blockDim%y + ty

!Filter coefficients into shared memory 29*4*8 B=0.928 KB
if (-bx .le. 0 .and. -bx .ge. lowfil .and. bx .le. lupfil) then
filtersub(-bx,1) = filter(-bx,1)
filtersub(-bx,2) = filter(-bx,2)
filtersub(-bx,3) = filter(-bx,3)
filtersub(-bx,4) = filter(-bx,4)
filtersub(bx,1) = filter(bx,1)
filtersub(bx,2) = filter(bx,2)
filtersub(bx,3) = filter(bx,3)
filtersub(bx,4) = filter(bx,4)
endif

! (1/2) d^2/dx^2
if (i .ge. nfl2 .and. i .le. nfu2 .and. j .ge. nfl3 .and. j .le. nfu3) then
!ibyz_c into shared memory 4*4*2*4 B= 0.128 KB
!To avoid Non-stride-1 accesses we change the arrangement of the array.
ibyz_fsub(tx,ty,1) = ibyz_f(1,i,j)
ibyz_fsub(tx,ty,2) = ibyz_f(2,i,j)
call syncthreads()
do kb=ibyz_fsub(tx,ty,1),ibyz_fsub(tx,ty,2),32
t112=0.d0;t121=0.d0;t122=0.d0;t212=0.d0;t221=0.d0;t222=0.d0;t211=0.d0
!part of psig into shared memory (32+15)*1*8*4*8 B=12.032 KB
k = ibyz_fsub(tx,ty,2)-kb
if (tz .le. 13 .and. tz .le. k .and. tz+kb-14 .ge. 0) then
psigsub(tz,1,tx,1,ty,2) = psig(kb+tz-14,1,i,1,j,2)
psigsub(tz,2,tx,1,ty,2) = psig(kb+tz-14,2,i,1,j,2)
psigsub(tz,1,tx,2,ty,2) = psig(kb+tz-14,1,i,2,j,2) 
psigsub(tz,2,tx,2,ty,2) = psig(kb+tz-14,2,i,2,j,2) 
end if
call syncthreads()
if (k .le. 13 .and. tz .le. 13) then
psigsub(tz+k+1,1,tx,1,ty,2) = psig(kb+tz-14+k+1,1,i,1,j,2)
psigsub(tz+k+1,2,tx,1,ty,2) = psig(kb+tz-14+k+1,2,i,1,j,2)
psigsub(tz+k+1,1,tx,2,ty,2) = psig(kb+tz-14+k+1,1,i,2,j,2)
psigsub(tz+k+1,2,tx,2,ty,2) = psig(kb+tz-14+k+1,2,i,2,j,2)
else
psigsub(tz+14,1,tx,1,ty,2) = psig(kb+tz,1,i,1,j,2)
psigsub(tz+14,2,tx,1,ty,2) = psig(kb+tz,2,i,1,j,2)
psigsub(tz+14,1,tx,2,ty,2) = psig(kb+tz,1,i,2,j,2) 
psigsub(tz+14,2,tx,2,ty,2) = psig(kb+tz,2,i,2,j,2)
end if
call syncthreads()

!we only use the llow do loop when needed

if (tz .le. k) then
llow = nfl1-(tz+kb)+14

if (nfl1-kb+14 .ge. 1) then
 if (llow .ge. 1) then
 do l=llow,13
  t112=t112+psigsub(tz+l,1,tx,1,ty,2)*filtersub(-14+l,1)+psigsub(tz+l,2,tx,1,ty,2)*filtersub(-14+l,2)
  t122=t122+psigsub(tz+l,1,tx,2,ty,2)*filtersub(-14+l,1)+psigsub(tz+l,2,tx,2,ty,2)*filtersub(-14+l,2)
  t212=t212+psigsub(tz+l,1,tx,1,ty,2)*filtersub(-14+l,3)+psigsub(tz+l,2,tx,1,ty,2)*filtersub(-14+l,4)
  t222=t222+psigsub(tz+l,1,tx,2,ty,2)*filtersub(-14+l,3)+psigsub(tz+l,2,tx,2,ty,2)*filtersub(-14+l,4)
 end do
 else
 do l=0,13
  t112=t112+psigsub(tz+l,1,tx,1,ty,2)*filtersub(-14+l,1)+psigsub(tz+l,2,tx,1,ty,2)*filtersub(-14+l,2)
  t122=t122+psigsub(tz+l,1,tx,2,ty,2)*filtersub(-14+l,1)+psigsub(tz+l,2,tx,2,ty,2)*filtersub(-14+l,2)
  t212=t212+psigsub(tz+l,1,tx,1,ty,2)*filtersub(-14+l,3)+psigsub(tz+l,2,tx,1,ty,2)*filtersub(-14+l,4)
  t222=t222+psigsub(tz+l,1,tx,2,ty,2)*filtersub(-14+l,3)+psigsub(tz+l,2,tx,2,ty,2)*filtersub(-14+l,4)
 enddo
 end if
else
 do l=0,13
  t112=t112+psigsub(tz+l,1,tx,1,ty,2)*filtersub(-14+l,1)+psigsub(tz+l,2,tx,1,ty,2)*filtersub(-14+l,2)
  t122=t122+psigsub(tz+l,1,tx,2,ty,2)*filtersub(-14+l,1)+psigsub(tz+l,2,tx,2,ty,2)*filtersub(-14+l,2)
  t212=t212+psigsub(tz+l,1,tx,1,ty,2)*filtersub(-14+l,3)+psigsub(tz+l,2,tx,1,ty,2)*filtersub(-14+l,4)
  t222=t222+psigsub(tz+l,1,tx,2,ty,2)*filtersub(-14+l,3)+psigsub(tz+l,2,tx,2,ty,2)*filtersub(-14+l,4)
 enddo
end if

end if

! Copy data from the end of the array to the beginning
call syncthreads()
if (tz .le. 13 .and. tz .le. k) then
psigsub(tz,1,tx,1,ty,2)=psigsub(tz+14,1,tx,1,ty,2)
psigsub(tz,2,tx,1,ty,2)=psigsub(tz+14,2,tx,1,ty,2)
psigsub(tz,1,tx,2,ty,2)=psigsub(tz+14,1,tx,2,ty,2)
psigsub(tz,2,tx,2,ty,2)=psigsub(tz+14,2,tx,2,ty,2)
end if
if (tz .le. 13 .and. k .ge. 14) then
psigsub(tz+14,1,tx,1,ty,2)=psigsub(tz+28,1,tx,1,ty,2)
psigsub(tz+14,2,tx,1,ty,2)=psigsub(tz+28,2,tx,1,ty,2)
psigsub(tz+14,1,tx,2,ty,2)=psigsub(tz+28,1,tx,2,ty,2)
psigsub(tz+14,2,tx,2,ty,2)=psigsub(tz+28,2,tx,2,ty,2)
end if
call syncthreads()
if (tz .ge. 28 .and. k .ge. 28) then 
psigsub(tz,1,tx,1,ty,2)=psigsub(tz+14,1,tx,1,ty,2) 
psigsub(tz,2,tx,1,ty,2)=psigsub(tz+14,2,tx,1,ty,2)
psigsub(tz,1,tx,2,ty,2)=psigsub(tz+14,1,tx,2,ty,2) 
psigsub(tz,2,tx,2,ty,2)=psigsub(tz+14,2,tx,2,ty,2) 
end if 
call syncthreads()


! Fill array with new data
if (tz .le. 13 .and. tz .le. k-1) then
psigsub(tz+32,1,tx,1,ty,2)=psig(tz+kb+32,1,i,1,j,2)
psigsub(tz+32,2,tx,1,ty,2)=psig(tz+kb+32,2,i,1,j,2)
psigsub(tz+32,1,tx,2,ty,2)=psig(tz+kb+32,1,i,2,j,2) 
psigsub(tz+32,2,tx,2,ty,2)=psig(tz+kb+32,2,i,2,j,2)
end if
call syncthreads()
if (tz .le. 13 .and. k .le. 31) then
psigsub(tz+k+1,1,tx,1,ty,2)=psig(tz+kb+k+1,1,i,1,j,2)
psigsub(tz+k+1,2,tx,1,ty,2)=psig(tz+kb+k+1,2,i,1,j,2)
psigsub(tz+k+1,1,tx,2,ty,2)=psig(tz+kb+k+1,1,i,2,j,2)
psigsub(tz+k+1,2,tx,2,ty,2)=psig(tz+kb+k+1,2,i,2,j,2)
end if
call syncthreads()


lupfilsub = 14
!we only use the lhigh do loop when needed
if (tz .le. k) then
lhigh = nfu1-(tz+kb)
 
if (nfu1-(kb+blockDim%z) .le. lupfilsub-1) then
 if (lhigh .le. lupfilsub-1) then
 do l=0,lhigh
 t112=t112+psigsub(tz+l,1,tx,1,ty,2)*filtersub(l,1)+psigsub(tz+l,2,tx,1,ty,2)*filtersub(l,2)
 t122=t122+psigsub(tz+l,1,tx,2,ty,2)*filtersub(l,1)+psigsub(tz+l,2,tx,2,ty,2)*filtersub(l,2)
 t212=t212+psigsub(tz+l,1,tx,1,ty,2)*filtersub(l,3)+psigsub(tz+l,2,tx,1,ty,2)*filtersub(l,4)
 t222=t222+psigsub(tz+l,1,tx,2,ty,2)*filtersub(l,3)+psigsub(tz+l,2,tx,2,ty,2)*filtersub(l,4)
 end do
 else
 do l=0,lupfilsub
 t112=t112+psigsub(tz+l,1,tx,1,ty,2)*filtersub(l,1)+psigsub(tz+l,2,tx,1,ty,2)*filtersub(l,2)
 t122=t122+psigsub(tz+l,1,tx,2,ty,2)*filtersub(l,1)+psigsub(tz+l,2,tx,2,ty,2)*filtersub(l,2)
 t212=t212+psigsub(tz+l,1,tx,1,ty,2)*filtersub(l,3)+psigsub(tz+l,2,tx,1,ty,2)*filtersub(l,4)
 t222=t222+psigsub(tz+l,1,tx,2,ty,2)*filtersub(l,3)+psigsub(tz+l,2,tx,2,ty,2)*filtersub(l,4)
 end do
 end if
else
 do l=0,lupfilsub
 t112=t112+psigsub(tz+l,1,tx,1,ty,2)*filtersub(l,1)+psigsub(tz+l,2,tx,1,ty,2)*filtersub(l,2)
 t122=t122+psigsub(tz+l,1,tx,2,ty,2)*filtersub(l,1)+psigsub(tz+l,2,tx,2,ty,2)*filtersub(l,2)
 t212=t212+psigsub(tz+l,1,tx,1,ty,2)*filtersub(l,3)+psigsub(tz+l,2,tx,1,ty,2)*filtersub(l,4)
 t222=t222+psigsub(tz+l,1,tx,2,ty,2)*filtersub(l,3)+psigsub(tz+l,2,tx,2,ty,2)*filtersub(l,4)
 end do 
end if

ekinvec(kb+tz,i,j)=ekinvec(kb+tz,i,j)+t112*psigsub(tz,1,tx,1,ty,2)+t212*psigsub(tz,2,tx,1,ty,2)+t122*psigsub(tz,1,tx,2,ty,2)+t222*psigsub(tz,2,tx,2,ty,2)

y(kb+tz,1,i,1,j,2)=y(kb+tz,1,i,1,j,2)+t112
y(kb+tz,2,i,1,j,2)=y(kb+tz,2,i,1,j,2)+t212
y(kb+tz,1,i,2,j,2)=y(kb+tz,1,i,2,j,2)+t122
y(kb+tz,2,i,2,j,2)=y(kb+tz,2,i,2,j,2)+t222
call syncthreads()

end if

   enddo
end if

    end subroutine convol_part8_kernel

attributes(global) subroutine convol_part9_kernel(psig,y,ekinvec,filter,ibyz_f,nfl1,nfl2,nfl3,nfu1,nfu2,nfu3,n1,n2,n3,lowfil,lupfil)
       real*8,device :: psig(0:n1,1:2,0:n2,1:2,0:n3,1:2),y(0:n1,1:2,0:n2,1:2,0:n3,1:2),ekinvec(0:n1,0:n2,0:n3),filter(lowfil:lupfil,4)
       integer*4,device :: ibyz_f(2,0:n2,0:n3)
       integer*4, value :: nfl1,nfl2,nfl3,nfu1,nfu2,nfu3,n1,n2,n3,lowfil,lupfil,l,llow,lhigh
       real*8 :: ekin,t111,t112,t121,t211,t122,t212,t221,t222
       integer*4 :: i, j, kb, k, bx, tx, ty, tz, ka

! submatrices are declared to be in CUDA shared memory

       integer*4, shared :: ibyz_fsub(0:0,0:7,1:2),lupfilsub
       real*8, shared :: filtersub(-14:14,1:4),psigsub(0:46,1:2,0:0,2:2,0:7,1:1)

! Start execution, first get thread indices

       tx = threadidx%x-1
       ty = threadidx%y-1
       tz = threadidx%z-1

! thread indices over block

       bx = tx+(ty*blockDim%x)+(tz*blockDim%x*blockDim%y)

! thread indices over grid

       i = (blockidx%x-1) * blockDim%x + tx
       j = (blockidx%y-1) * blockDim%y + ty

!Filter coefficients into shared memory 29*4*8 B=0.928 KB
if (-bx .le. 0 .and. -bx .ge. lowfil .and. bx .le. lupfil) then
filtersub(-bx,1) = filter(-bx,1)
filtersub(-bx,2) = filter(-bx,2)
filtersub(-bx,3) = filter(-bx,3)
filtersub(-bx,4) = filter(-bx,4)
filtersub(bx,1) = filter(bx,1)
filtersub(bx,2) = filter(bx,2)
filtersub(bx,3) = filter(bx,3)
filtersub(bx,4) = filter(bx,4)
endif

! (1/2) d^2/dx^2
! (1/2) d^2/dx^2
if (i .ge. nfl2 .and. i .le. nfu2 .and. j .ge. nfl3 .and. j .le. nfu3) then
!ibyz_c into shared memory 4*4*2*4 B= 0.128 KB
!To avoid Non-stride-1 accesses we change the arrangement of the array.
ibyz_fsub(tx,ty,1) = ibyz_f(1,i,j)
ibyz_fsub(tx,ty,2) = ibyz_f(2,i,j)
call syncthreads()
do kb=ibyz_fsub(tx,ty,1),ibyz_fsub(tx,ty,2),32
t112=0.d0;t121=0.d0;t122=0.d0;t212=0.d0;t221=0.d0;t222=0.d0;t211=0.d0
!part of psig into shared memory (32+15)*4*4*2*8 B=12.032 KB
k = ibyz_fsub(tx,ty,2)-kb
if (tz .le. 13 .and. tz .le. k .and. tz+kb-14 .ge. 0) then
psigsub(tz,1,tx,2,ty,1) = psig(kb+tz-14,1,i,2,j,1)
psigsub(tz,2,tx,2,ty,1) = psig(kb+tz-14,2,i,2,j,1)
end if
call syncthreads()
if (k .le. 13 .and. tz .le. 13) then
psigsub(tz+k+1,1,tx,2,ty,1) = psig(kb+tz-14+k+1,1,i,2,j,1)
psigsub(tz+k+1,2,tx,2,ty,1) = psig(kb+tz-14+k+1,2,i,2,j,1)
else
psigsub(tz+14,1,tx,2,ty,1) = psig(kb+tz,1,i,2,j,1)
psigsub(tz+14,2,tx,2,ty,1) = psig(kb+tz,2,i,2,j,1)
end if
call syncthreads()

!we only use the llow do loop when needed

if (tz .le. k) then
llow = nfl1-(tz+kb)+14

if (nfl1-kb+14 .ge. 1) then
 if (llow .ge. 1) then
 do l=llow,13
  t121=t121+psigsub(tz+l,1,tx,2,ty,1)*filtersub(-14+l,1)+psigsub(tz+l,2,tx,2,ty,1)*filtersub(-14+l,2)
  t221=t221+psigsub(tz+l,1,tx,2,ty,1)*filtersub(-14+l,3)+psigsub(tz+l,2,tx,2,ty,1)*filtersub(-14+l,4)
 end do
 else
 do l=0,13
  t121=t121+psigsub(tz+l,1,tx,2,ty,1)*filtersub(-14+l,1)+psigsub(tz+l,2,tx,2,ty,1)*filtersub(-14+l,2)
  t221=t221+psigsub(tz+l,1,tx,2,ty,1)*filtersub(-14+l,3)+psigsub(tz+l,2,tx,2,ty,1)*filtersub(-14+l,4) 
 enddo
 end if
else
 do l=0,13
  t121=t121+psigsub(tz+l,1,tx,2,ty,1)*filtersub(-14+l,1)+psigsub(tz+l,2,tx,2,ty,1)*filtersub(-14+l,2)
  t221=t221+psigsub(tz+l,1,tx,2,ty,1)*filtersub(-14+l,3)+psigsub(tz+l,2,tx,2,ty,1)*filtersub(-14+l,4)
 enddo
end if

end if

! Copy data from the end of the array to the beginning
call syncthreads()
if (tz .le. 13 .and. tz .le. k) then
psigsub(tz,1,tx,2,ty,1)=psigsub(tz+14,1,tx,2,ty,1)
psigsub(tz,2,tx,2,ty,1)=psigsub(tz+14,2,tx,2,ty,1)
end if
if (tz .le. 13 .and. k .ge. 14) then
psigsub(tz+14,1,tx,2,ty,1)=psigsub(tz+28,1,tx,2,ty,1)
psigsub(tz+14,2,tx,2,ty,1)=psigsub(tz+28,2,tx,2,ty,1)
end if
call syncthreads()
if (tz .ge. 28 .and. k .ge. 28) then 
psigsub(tz,1,tx,2,ty,1)=psigsub(tz+14,1,tx,2,ty,1) 
psigsub(tz,2,tx,2,ty,1)=psigsub(tz+14,2,tx,2,ty,1) 
end if 
call syncthreads()


! Fill array with new data
if (tz .le. 13 .and. tz .le. k-1) then
psigsub(tz+32,1,tx,2,ty,1)=psig(tz+kb+32,1,i,2,j,1)
psigsub(tz+32,2,tx,2,ty,1)=psig(tz+kb+32,2,i,2,j,1)
end if
call syncthreads()
if (tz .le. 13 .and. k .le. 31) then
psigsub(tz+k+1,1,tx,2,ty,1)=psig(tz+kb+k+1,1,i,2,j,1)
psigsub(tz+k+1,2,tx,2,ty,1)=psig(tz+kb+k+1,2,i,2,j,1)
end if
call syncthreads()


lupfilsub = 14
!we only use the lhigh do loop when needed
if (tz .le. k) then
lhigh = nfu1-(tz+kb)
 
if (nfu1-(kb+blockDim%z) .le. lupfilsub-1) then
 if (lhigh .le. lupfilsub-1) then
 do l=0,lhigh
 t121=t121+psigsub(tz+l,1,tx,2,ty,1)*filtersub(l,1)+psigsub(tz+l,2,tx,2,ty,1)*filtersub(l,2)
 t221=t221+psigsub(tz+l,1,tx,2,ty,1)*filtersub(l,3)+psigsub(tz+l,2,tx,2,ty,1)*filtersub(l,4)
 end do
 else
 do l=0,lupfilsub
 t121=t121+psigsub(tz+l,1,tx,2,ty,1)*filtersub(l,1)+psigsub(tz+l,2,tx,2,ty,1)*filtersub(l,2)
 t221=t221+psigsub(tz+l,1,tx,2,ty,1)*filtersub(l,3)+psigsub(tz+l,2,tx,2,ty,1)*filtersub(l,4)
 end do
 end if
else
 do l=0,lupfilsub
 t121=t121+psigsub(tz+l,1,tx,2,ty,1)*filtersub(l,1)+psigsub(tz+l,2,tx,2,ty,1)*filtersub(l,2)
 t221=t221+psigsub(tz+l,1,tx,2,ty,1)*filtersub(l,3)+psigsub(tz+l,2,tx,2,ty,1)*filtersub(l,4)
 end do 
end if

ekinvec(kb+tz,i,j)=ekinvec(kb+tz,i,j)+t121*psigsub(tz,1,tx,2,ty,1)+t221*psigsub(tz,2,tx,2,ty,1)
y(kb+tz,1,i,2,j,1)=y(kb+tz,1,i,2,j,1)+t121
y(kb+tz,2,i,2,j,1)=y(kb+tz,2,i,2,j,1)+t221
call syncthreads()

end if

   enddo
end if

    end subroutine convol_part9_kernel


attributes(global) subroutine convol_part11_kernel(psig,y,ekinvec,filter,ibyz_f,nfl1,nfl2,nfl3,nfu1,nfu2,nfu3,n1,n2,n3,lowfil,lupfil)
       real*8,device :: psig(0:n1,1:2,0:n2,1:2,0:n3,1:2),y(0:n1,1:2,0:n2,1:2,0:n3,1:2),ekinvec(0:n1,0:n2,0:n3),filter(lowfil:lupfil,4)
       integer*4,device :: ibyz_f(2,0:n2,0:n3)
       integer*4, value :: nfl1,nfl2,nfl3,nfu1,nfu2,nfu3,n1,n2,n3,lowfil,lupfil,l,llow,lhigh
       real*8 :: ekin,t111,t112,t121,t211,t122,t212,t221,t222
       integer*4 :: i, j, kb, k, bx, tx, ty, tz, ka

! submatrices are declared to be in CUDA shared memory

       integer*4, shared :: ibyz_fsub(0:0,0:7,1:2),lupfilsub
       real*8, shared :: filtersub(-14:14,4:4),psigsub(0:46,2:2,0:0,1:1,0:7,1:1)

! Start execution, first get thread indices

       tx = threadidx%x-1
       ty = threadidx%y-1
       tz = threadidx%z-1

! thread indices over block

       bx = tx+(ty*blockDim%x)+(tz*blockDim%x*blockDim%y)

! thread indices over grid

       i = (blockidx%x-1) * blockDim%x + tx
       j = (blockidx%y-1) * blockDim%y + ty

!Filter coefficients into shared memory 29*4*8 B=0.232 KB
if (-bx .le. 0 .and. -bx .ge. lowfil .and. bx .le. lupfil) then
filtersub(-bx,4) = filter(-bx,4)
filtersub(bx,4) = filter(bx,4)
endif

! (1/2) d^2/dx^2
if (i .ge. nfl2 .and. i .le. nfu2 .and. j .ge. nfl3 .and. j .le. nfu3) then
!ibyz_c into shared memory 4*4*2*4 B= 0.128 KB
!To avoid Non-stride-1 accesses we change the arrangement of the array.
ibyz_fsub(tx,ty,1) = ibyz_f(1,i,j)
ibyz_fsub(tx,ty,2) = ibyz_f(2,i,j)
call syncthreads()
do kb=ibyz_fsub(tx,ty,1),ibyz_fsub(tx,ty,2),32
t112=0.d0;t121=0.d0;t122=0.d0;t212=0.d0;t221=0.d0;t222=0.d0;t211=0.d0
!part of psig into shared memory (32+15)*4*4*2*8 B=12.032 KB
k = ibyz_fsub(tx,ty,2)-kb
if (tz .le. 13 .and. tz .le. k .and. tz+kb-14 .ge. 0) then
psigsub(tz,2,tx,1,ty,1) = psig(kb+tz-14,2,i,1,j,1)
end if
call syncthreads()
if (k .le. 13 .and. tz .le. 13) then
psigsub(tz+k+1,2,tx,1,ty,1) = psig(kb+tz-14+k+1,2,i,1,j,1)
else
psigsub(tz+14,2,tx,1,ty,1) = psig(kb+tz,2,i,1,j,1)
end if
call syncthreads()

!we only use the llow do loop when needed

if (tz .le. k) then
llow = nfl1-(tz+kb)+14

if (nfl1-kb+14 .ge. 1) then
 if (llow .ge. 1) then
 do l=llow,13
  t211=t211+psigsub(tz+l,2,tx,1,ty,1)*filtersub(-14+l,4)
 end do
 else
 do l=0,13
  t211=t211+psigsub(tz+l,2,tx,1,ty,1)*filtersub(-14+l,4)
 enddo
 end if
else
 do l=0,13
  t211=t211+psigsub(tz+l,2,tx,1,ty,1)*filtersub(-14+l,4)
 enddo
end if

end if

! Copy data from the end of the array to the beginning
call syncthreads()
if (tz .le. 13 .and. tz .le. k) then
psigsub(tz,2,tx,1,ty,1)=psigsub(tz+14,2,tx,1,ty,1)
end if
if (tz .le. 13 .and. k .ge. 14) then
psigsub(tz+14,2,tx,1,ty,1)=psigsub(tz+28,2,tx,1,ty,1)
end if
call syncthreads()
if (tz .ge. 28 .and. k .ge. 28) then 
psigsub(tz,2,tx,1,ty,1)=psigsub(tz+14,2,tx,1,ty,1) 
end if 
call syncthreads()


! Fill array with new data
if (tz .le. 13 .and. tz .le. k-1) then
psigsub(tz+32,2,tx,1,ty,1)=psig(tz+kb+32,2,i,1,j,1)
end if
call syncthreads()
if (tz .le. 13 .and. k .le. 31) then
psigsub(tz+k+1,2,tx,1,ty,1)=psig(tz+kb+k+1,2,i,1,j,1)
end if
call syncthreads()


lupfilsub = 14
!we only use the lhigh do loop when needed
if (tz .le. k) then
lhigh = nfu1-(tz+kb)
 
if (nfu1-(kb+blockDim%z) .le. lupfilsub-1) then
 if (lhigh .le. lupfilsub-1) then
 do l=0,lhigh
 t211=t211+psigsub(tz+l,2,tx,1,ty,1)*filtersub(l,4)
 end do
 else
 do l=0,lupfilsub
 t211=t211+psigsub(tz+l,2,tx,1,ty,1)*filtersub(l,4)
 end do
 end if
else
 do l=0,lupfilsub
 t211=t211+psigsub(tz+l,2,tx,1,ty,1)*filtersub(l,4)
 end do 
end if

ekinvec(kb+tz,i,j)=ekinvec(kb+tz,i,j)+t211*psigsub(tz,2,tx,1,ty,1)
y(kb+tz,2,i,1,j,1)=y(kb+tz,2,i,1,j,1)+t211
call syncthreads()

end if

   enddo
end if

    end subroutine convol_part11_kernel

attributes(global) subroutine convol_part12_kernel(psig,y,ekinvec,filter, ibxz_f,nfl1,nfl2,nfl3,nfu1,nfu2,nfu3,n1,n2,n3,lowfil,lupfil)
       real*8,device :: psig(0:n1,1:2,0:n2,1:2,0:n3,1:2),y(0:n1,1:2,0:n2,1:2,0:n3,1:2),ekinvec(0:n1,0:n2,0:n3),filter(lowfil:lupfil,4)
       integer*4,device :: ibxz_f(2,0:n1,0:n3)
       integer*4, value :: nfl1,nfl2,nfl3,nfu1,nfu2,nfu3,n1,n2,n3,lowfil,lupfil,l,llow,lhigh
       real*8 :: ekin,t111,t112,t121,t211,t122,t212,t221,t222
       integer*4 :: i, j, kb, k, bx, tx, ty, tz, ka

! submatrices are declared to be in CUDA shared memory

       integer*4, shared :: ibxz_fsub(0:0,0:7,1:2),lupfilsub
       real*8, shared :: filtersub(-14:14,1:4),psigsub(0:0,1:2,0:46,1:2,0:7,2:2)

! Start execution, first get thread indices

       tx = threadidx%x-1
       ty = threadidx%y-1
       tz = threadidx%z-1

! thread indices over block

       bx = tx+(ty*blockDim%x)+(tz*blockDim%x*blockDim%y)

! thread indices over grid

       i = (blockidx%x-1) * blockDim%x + tx
       j = (blockidx%y-1) * blockDim%y + ty

!Filter coefficients into shared memory 29*4*8 B=0.928 KB
if (-bx .le. 0 .and. -bx .ge. lowfil .and. bx .le. lupfil) then
filtersub(-bx,1) = filter(-bx,1)
filtersub(-bx,2) = filter(-bx,2)
filtersub(-bx,3) = filter(-bx,3)
filtersub(-bx,4) = filter(-bx,4)
filtersub(bx,1) = filter(bx,1)
filtersub(bx,2) = filter(bx,2)
filtersub(bx,3) = filter(bx,3)
filtersub(bx,4) = filter(bx,4)
endif

! (1/2) d^2/dx^2
if (i .ge. nfl1 .and. i .le. nfu1 .and. j .ge. nfl3 .and. j .le. nfu3) then
!ibxz_f into shared memory 4*4*2*4 B= 0.128 KB
!To avoid Non-stride-1 accesses we change the arrangement of the array.
ibxz_fsub(tx,ty,1) = ibxz_f(1,i,j)
ibxz_fsub(tx,ty,2) = ibxz_f(2,i,j)
call syncthreads()
do kb=ibxz_fsub(tx,ty,1),ibxz_fsub(tx,ty,2),32
t112=0.d0;t121=0.d0;t122=0.d0;t212=0.d0;t221=0.d0;t222=0.d0;t211=0.d0
!part of psig into shared memory (32+15)*1*8*4*8 B=12.032 KB
k = ibxz_fsub(tx,ty,2)-kb
if (tz .le. 13 .and. tz .le. k .and. tz+kb-14 .ge. 0) then
psigsub(tx,1,tz,1,ty,2) = psig(i,1,kb+tz-14,1,j,2)
psigsub(tx,2,tz,1,ty,2) = psig(i,2,kb+tz-14,1,j,2)
psigsub(tx,1,tz,2,ty,2) = psig(i,1,kb+tz-14,2,j,2) 
psigsub(tx,2,tz,2,ty,2) = psig(i,2,kb+tz-14,2,j,2) 
end if
call syncthreads()
if (k .le. 13 .and. tz .le. 13) then
psigsub(tx,1,tz+k+1,1,ty,2) = psig(i,1,kb+tz-14+k+1,1,j,2)
psigsub(tx,2,tz+k+1,1,ty,2) = psig(i,2,kb+tz-14+k+1,1,j,2)
psigsub(tx,1,tz+k+1,2,ty,2) = psig(i,1,kb+tz-14+k+1,2,j,2)
psigsub(tx,2,tz+k+1,2,ty,2) = psig(i,2,kb+tz-14+k+1,2,j,2)
else
psigsub(tx,1,tz+14,1,ty,2) = psig(i,1,kb+tz,1,j,2)
psigsub(tx,2,tz+14,1,ty,2) = psig(i,2,kb+tz,1,j,2)
psigsub(tx,1,tz+14,2,ty,2) = psig(i,1,kb+tz,2,j,2) 
psigsub(tx,2,tz+14,2,ty,2) = psig(i,2,kb+tz,2,j,2)
end if
call syncthreads()

!we only use the llow do loop when needed

if (tz .le. k) then
llow = nfl2-(tz+kb)+14

if (nfl2-kb+14 .ge. 1) then
 if (llow .ge. 1) then
 do l=llow,13
  t112=t112+psigsub(tx,1,tz+l,1,ty,2)*filtersub(-14+l,1)+psigsub(tx,1,tz+l,2,ty,2)*filtersub(-14+l,2)
  t122=t122+psigsub(tx,1,tz+l,1,ty,2)*filtersub(-14+l,3)+psigsub(tx,1,tz+l,2,ty,2)*filtersub(-14+l,4)
  t212=t212+psigsub(tx,2,tz+l,1,ty,2)*filtersub(-14+l,1)+psigsub(tx,2,tz+l,2,ty,2)*filtersub(-14+l,2)
  t222=t222+psigsub(tx,2,tz+l,1,ty,2)*filtersub(-14+l,3)+psigsub(tx,2,tz+l,2,ty,2)*filtersub(-14+l,4)
 end do
 else
 do l=0,13
  t112=t112+psigsub(tx,1,tz+l,1,ty,2)*filtersub(-14+l,1)+psigsub(tx,1,tz+l,2,ty,2)*filtersub(-14+l,2)
  t122=t122+psigsub(tx,1,tz+l,1,ty,2)*filtersub(-14+l,3)+psigsub(tx,1,tz+l,2,ty,2)*filtersub(-14+l,4)
  t212=t212+psigsub(tx,2,tz+l,1,ty,2)*filtersub(-14+l,1)+psigsub(tx,2,tz+l,2,ty,2)*filtersub(-14+l,2)
  t222=t222+psigsub(tx,2,tz+l,1,ty,2)*filtersub(-14+l,3)+psigsub(tx,2,tz+l,2,ty,2)*filtersub(-14+l,4)
 enddo
 end if
else
 do l=0,13
  t112=t112+psigsub(tx,1,tz+l,1,ty,2)*filtersub(-14+l,1)+psigsub(tx,1,tz+l,2,ty,2)*filtersub(-14+l,2)
  t122=t122+psigsub(tx,1,tz+l,1,ty,2)*filtersub(-14+l,3)+psigsub(tx,1,tz+l,2,ty,2)*filtersub(-14+l,4)
  t212=t212+psigsub(tx,2,tz+l,1,ty,2)*filtersub(-14+l,1)+psigsub(tx,2,tz+l,2,ty,2)*filtersub(-14+l,2)
  t222=t222+psigsub(tx,2,tz+l,1,ty,2)*filtersub(-14+l,3)+psigsub(tx,2,tz+l,2,ty,2)*filtersub(-14+l,4)
 enddo
end if

end if

! Copy data from the end of the array to the beginning
call syncthreads()
if (tz .le. 13 .and. tz .le. k) then
psigsub(tx,1,tz,1,ty,2)=psigsub(tx,1,tz+14,1,ty,2)
psigsub(tx,2,tz,1,ty,2)=psigsub(tx,2,tz+14,1,ty,2)
psigsub(tx,1,tz,2,ty,2)=psigsub(tx,1,tz+14,2,ty,2)
psigsub(tx,2,tz,2,ty,2)=psigsub(tx,2,tz+14,2,ty,2)
end if
if (tz .le. 13 .and. k .ge. 14) then
psigsub(tx,1,tz+14,1,ty,2)=psigsub(tx,1,tz+28,1,ty,2)
psigsub(tx,2,tz+14,1,ty,2)=psigsub(tx,2,tz+28,1,ty,2)
psigsub(tx,1,tz+14,2,ty,2)=psigsub(tx,1,tz+28,2,ty,2)
psigsub(tx,2,tz+14,2,ty,2)=psigsub(tx,2,tz+28,2,ty,2)
end if
call syncthreads()
if (tz .ge. 28 .and. k .ge. 28) then 
psigsub(tx,1,tz,1,ty,2)=psigsub(tx,1,tz+14,1,ty,2) 
psigsub(tx,2,tz,1,ty,2)=psigsub(tx,2,tz+14,1,ty,2)
psigsub(tx,1,tz,2,ty,2)=psigsub(tx,1,tz+14,2,ty,2) 
psigsub(tx,2,tz,2,ty,2)=psigsub(tx,2,tz+14,2,ty,2) 
end if 
call syncthreads()


! Fill array with new data
if (tz .le. 13 .and. tz .le. k-1) then
psigsub(tx,1,tz+32,1,ty,2)=psig(i,1,tz+kb+32,1,j,2)
psigsub(tx,2,tz+32,1,ty,2)=psig(i,2,tz+kb+32,1,j,2)
psigsub(tx,1,tz+32,2,ty,2)=psig(i,1,tz+kb+32,2,j,2) 
psigsub(tx,2,tz+32,2,ty,2)=psig(i,2,tz+kb+32,2,j,2)
end if
call syncthreads()
if (tz .le. 13 .and. k .le. 31) then
psigsub(tx,1,tz+k+1,1,ty,2)=psig(i,1,tz+kb+k+1,1,j,2)
psigsub(tx,2,tz+k+1,1,ty,2)=psig(i,2,tz+kb+k+1,1,j,2)
psigsub(tx,1,tz+k+1,2,ty,2)=psig(i,1,tz+kb+k+1,2,j,2)
psigsub(tx,2,tz+k+1,2,ty,2)=psig(i,2,tz+kb+k+1,2,j,2)
end if
call syncthreads()


lupfilsub = 14
!we only use the lhigh do loop when needed
if (tz .le. k) then
lhigh = nfu2-(tz+kb)
 
if (nfu2-(kb+blockDim%z) .le. lupfilsub-1) then
 if (lhigh .le. lupfilsub-1) then
 do l=0,lhigh
  t112=t112+psigsub(tx,1,tz+l,1,ty,2)*filtersub(l,1)+psigsub(tx,1,tz+l,2,ty,2)*filtersub(l,2)
  t122=t122+psigsub(tx,1,tz+l,1,ty,2)*filtersub(l,3)+psigsub(tx,1,tz+l,2,ty,2)*filtersub(l,4)
  t212=t212+psigsub(tx,2,tz+l,1,ty,2)*filtersub(l,1)+psigsub(tx,2,tz+l,2,ty,2)*filtersub(l,2)
  t222=t222+psigsub(tx,2,tz+l,1,ty,2)*filtersub(l,3)+psigsub(tx,2,tz+l,2,ty,2)*filtersub(l,4)
 end do
 else
 do l=0,lupfilsub
  t112=t112+psigsub(tx,1,tz+l,1,ty,2)*filtersub(l,1)+psigsub(tx,1,tz+l,2,ty,2)*filtersub(l,2)
  t122=t122+psigsub(tx,1,tz+l,1,ty,2)*filtersub(l,3)+psigsub(tx,1,tz+l,2,ty,2)*filtersub(l,4)
  t212=t212+psigsub(tx,2,tz+l,1,ty,2)*filtersub(l,1)+psigsub(tx,2,tz+l,2,ty,2)*filtersub(l,2)
  t222=t222+psigsub(tx,2,tz+l,1,ty,2)*filtersub(l,3)+psigsub(tx,2,tz+l,2,ty,2)*filtersub(l,4)
 end do
 end if
else
 do l=0,lupfilsub
  t112=t112+psigsub(tx,1,tz+l,1,ty,2)*filtersub(l,1)+psigsub(tx,1,tz+l,2,ty,2)*filtersub(l,2)
  t122=t122+psigsub(tx,1,tz+l,1,ty,2)*filtersub(l,3)+psigsub(tx,1,tz+l,2,ty,2)*filtersub(l,4)
  t212=t212+psigsub(tx,2,tz+l,1,ty,2)*filtersub(l,1)+psigsub(tx,2,tz+l,2,ty,2)*filtersub(l,2)
  t222=t222+psigsub(tx,2,tz+l,1,ty,2)*filtersub(l,3)+psigsub(tx,2,tz+l,2,ty,2)*filtersub(l,4)
 end do 
end if

ekinvec(i,kb+tz,j)=ekinvec(i,kb+tz,j)+t112*psigsub(tx,1,tz,1,ty,2)+t212*psigsub(tx,2,tz,1,ty,2)+t122*psigsub(tx,1,tz,2,ty,2)+t222*psigsub(tx,2,tz,2,ty,2)

y(i,1,kb+tz,1,j,2)=y(i,1,kb+tz,1,j,2)+t112
y(i,2,kb+tz,1,j,2)=y(i,2,kb+tz,1,j,2)+t212
y(i,1,kb+tz,2,j,2)=y(i,1,kb+tz,2,j,2)+t122
y(i,2,kb+tz,2,j,2)=y(i,2,kb+tz,2,j,2)+t222
call syncthreads()

end if

   enddo
end if

    end subroutine convol_part12_kernel


attributes(global) subroutine convol_part14_kernel(psig,y,ekinvec,filter,ibxz_f,nfl1,nfl2,nfl3,nfu1,nfu2,nfu3,n1,n2,n3,lowfil,lupfil)
       real*8,device :: psig(0:n1,1:2,0:n2,1:2,0:n3,1:2),y(0:n1,1:2,0:n2,1:2,0:n3,1:2),ekinvec(0:n1,0:n2,0:n3),filter(lowfil:lupfil,4)
       integer*4,device :: ibxz_f(2,0:n1,0:n3)
       integer*4, value :: nfl1,nfl2,nfl3,nfu1,nfu2,nfu3,n1,n2,n3,lowfil,lupfil,l,llow,lhigh
       real*8 :: ekin,t111,t112,t121,t211,t122,t212,t221,t222
       integer*4 :: i, j, kb, k, bx, tx, ty, tz, ka

! submatrices are declared to be in CUDA shared memory

       integer*4, shared :: ibxz_fsub(0:0,0:7,1:2),lupfilsub
       real*8, shared :: filtersub(-14:14,1:4),psigsub(0:0,2:2,0:46,1:2,0:7,1:1)

! Start execution, first get thread indices

       tx = threadidx%x-1
       ty = threadidx%y-1
       tz = threadidx%z-1

! thread indices over block

       bx = tx+(ty*blockDim%x)+(tz*blockDim%x*blockDim%y)

! thread indices over grid

       i = (blockidx%x-1) * blockDim%x + tx
       j = (blockidx%y-1) * blockDim%y + ty

!Filter coefficients into shared memory 29*4*8 B=0.928 KB
if (-bx .le. 0 .and. -bx .ge. lowfil .and. bx .le. lupfil) then
filtersub(-bx,1) = filter(-bx,1)
filtersub(-bx,2) = filter(-bx,2)
filtersub(-bx,3) = filter(-bx,3)
filtersub(-bx,4) = filter(-bx,4)
filtersub(bx,1) = filter(bx,1)
filtersub(bx,2) = filter(bx,2)
filtersub(bx,3) = filter(bx,3)
filtersub(bx,4) = filter(bx,4)
endif

! (1/2) d^2/dx^2
if (i .ge. nfl1 .and. i .le. nfu1 .and. j .ge. nfl3 .and. j .le. nfu3) then
!ibxz_f into shared memory 4*4*2*4 B= 0.128 KB
!To avoid Non-stride-1 accesses we change the arrangement of the array.
ibxz_fsub(tx,ty,1) = ibxz_f(1,i,j)
ibxz_fsub(tx,ty,2) = ibxz_f(2,i,j)
call syncthreads()
do kb=ibxz_fsub(tx,ty,1),ibxz_fsub(tx,ty,2),32
t112=0.d0;t121=0.d0;t122=0.d0;t212=0.d0;t221=0.d0;t222=0.d0;t211=0.d0
!part of psig into shared memory (32+15)*4*4*2*8 B=12.032 KB
k = ibxz_fsub(tx,ty,2)-kb
if (tz .le. 13 .and. tz .le. k .and. tz+kb-14 .ge. 0) then
psigsub(tx,2,tz,1,ty,1) = psig(i,2,kb+tz-14,1,j,1)
psigsub(tx,2,tz,2,ty,1) = psig(i,2,kb+tz-14,2,j,1)
end if
call syncthreads()
if (k .le. 13 .and. tz .le. 13) then
psigsub(tx,2,tz+k+1,1,ty,1) = psig(i,2,kb+tz-14+k+1,1,j,1)
psigsub(tx,2,tz+k+1,2,ty,1) = psig(i,2,kb+tz-14+k+1,2,j,1)
else
psigsub(tx,2,tz+14,1,ty,1) = psig(i,2,kb+tz,1,j,1)
psigsub(tx,2,tz+14,2,ty,1) = psig(i,2,kb+tz,2,j,1)
end if
call syncthreads()

!we only use the llow do loop when needed

if (tz .le. k) then
llow = nfl2-(tz+kb)+14

if (nfl2-kb+14 .ge. 1) then
 if (llow .ge. 1) then
 do l=llow,13
  t211=t211+psigsub(tx,2,tz+l,1,ty,1)*filtersub(-14+l,1)+psigsub(tx,2,tz+l,2,ty,1)*filtersub(-14+l,2)
  t221=t221+psigsub(tx,2,tz+l,1,ty,1)*filtersub(-14+l,3)+psigsub(tx,2,tz+l,2,ty,1)*filtersub(-14+l,4)
 end do
 else
 do l=0,13
  t211=t211+psigsub(tx,2,tz+l,1,ty,1)*filtersub(-14+l,1)+psigsub(tx,2,tz+l,2,ty,1)*filtersub(-14+l,2) 
  t221=t221+psigsub(tx,2,tz+l,1,ty,1)*filtersub(-14+l,3)+psigsub(tx,2,tz+l,2,ty,1)*filtersub(-14+l,4)
 enddo
 end if
else
 do l=0,13
  t211=t211+psigsub(tx,2,tz+l,1,ty,1)*filtersub(-14+l,1)+psigsub(tx,2,tz+l,2,ty,1)*filtersub(-14+l,2)  
  t221=t221+psigsub(tx,2,tz+l,1,ty,1)*filtersub(-14+l,3)+psigsub(tx,2,tz+l,2,ty,1)*filtersub(-14+l,4)
 enddo
end if

end if

! Copy data from the end of the array to the beginning
call syncthreads()
if (tz .le. 13 .and. tz .le. k) then
psigsub(tx,2,tz,1,ty,1) = psigsub(tx,2,tz+14,1,ty,1)
psigsub(tx,2,tz,2,ty,1) = psigsub(tx,2,tz+14,2,ty,1)
end if
if (tz .le. 13 .and. k .ge. 14) then
psigsub(tx,2,tz+14,1,ty,1) = psigsub(tx,2,tz+28,1,ty,1)
psigsub(tx,2,tz+14,2,ty,1) = psigsub(tx,2,tz+28,2,ty,1)
end if
call syncthreads()
if (tz .ge. 28 .and. k .ge. 28) then 
psigsub(tx,2,tz,1,ty,1) = psigsub(tx,2,tz+14,1,ty,1)
psigsub(tx,2,tz,2,ty,1) = psigsub(tx,2,tz+14,2,ty,1)
end if 
call syncthreads()


! Fill array with new data
if (tz .le. 13 .and. tz .le. k-1) then
psigsub(tx,2,tz+32,1,ty,1) = psig(i,2,kb+tz+32,1,j,1)
psigsub(tx,2,tz+32,2,ty,1) = psig(i,2,kb+tz+32,2,j,1)
end if
call syncthreads()
if (tz .le. 13 .and. k .le. 31) then
psigsub(tx,2,tz+k+1,1,ty,1) = psig(i,2,tz+kb+k+1,1,j,1)
psigsub(tx,2,tz+k+1,2,ty,1) = psig(i,2,tz+kb+k+1,2,j,1)
end if
call syncthreads()


lupfilsub = 14
!we only use the lhigh do loop when needed
if (tz .le. k) then
lhigh = nfu2-(tz+kb)
 
if (nfu2-(kb+blockDim%z) .le. lupfilsub-1) then
 if (lhigh .le. lupfilsub-1) then
 do l=0,lhigh
  t211=t211+psigsub(tx,2,tz+l,1,ty,1)*filtersub(l,1)+psigsub(tx,2,tz+l,2,ty,1)*filtersub(l,2)
  t221=t221+psigsub(tx,2,tz+l,1,ty,1)*filtersub(l,3)+psigsub(tx,2,tz+l,2,ty,1)*filtersub(l,4)
 end do
 else
 do l=0,lupfilsub
  t211=t211+psigsub(tx,2,tz+l,1,ty,1)*filtersub(l,1)+psigsub(tx,2,tz+l,2,ty,1)*filtersub(l,2)
  t221=t221+psigsub(tx,2,tz+l,1,ty,1)*filtersub(l,3)+psigsub(tx,2,tz+l,2,ty,1)*filtersub(l,4)
 end do
 end if
else
 do l=0,lupfilsub
  t211=t211+psigsub(tx,2,tz+l,1,ty,1)*filtersub(l,1)+psigsub(tx,2,tz+l,2,ty,1)*filtersub(l,2)
  t221=t221+psigsub(tx,2,tz+l,1,ty,1)*filtersub(l,3)+psigsub(tx,2,tz+l,2,ty,1)*filtersub(l,4)
 end do 
end if

ekinvec(i,kb+tz,j)=ekinvec(i,kb+tz,j)+t211*psigsub(tx,2,tz,1,ty,1)+t221*psigsub(tx,2,tz,2,ty,1)
y(i,2,kb+tz,1,j,1)=y(i,2,kb+tz,1,j,1)+t211
y(i,2,kb+tz,2,j,1)=y(i,2,kb+tz,2,j,1)+t221
call syncthreads()

end if

   enddo
end if

    end subroutine convol_part14_kernel


attributes(global) subroutine convol_part15_kernel(psig,y,ekinvec,filter,ibxz_f,nfl1,nfl2,nfl3,nfu1,nfu2,nfu3,n1,n2,n3,lowfil,lupfil)
       real*8,device :: psig(0:n1,1:2,0:n2,1:2,0:n3,1:2),y(0:n1,1:2,0:n2,1:2,0:n3,1:2),ekinvec(0:n1,0:n2,0:n3),filter(lowfil:lupfil,4)
       integer*4,device :: ibxz_f(2,0:n1,0:n3)
       integer*4, value :: nfl1,nfl2,nfl3,nfu1,nfu2,nfu3,n1,n2,n3,lowfil,lupfil,l,llow,lhigh
       real*8 :: ekin,t111,t112,t121,t211,t122,t212,t221,t222
       integer*4 :: i, j, kb, k, bx, tx, ty, tz, ka

! submatrices are declared to be in CUDA shared memory

       integer*4, shared :: ibxz_fsub(0:0,0:7,1:2),lupfilsub
       real*8, shared :: filtersub(-14:14,4:4),psigsub(0:0,1:1,0:46,2:2,0:7,1:1)

! Start execution, first get thread indices

       tx = threadidx%x-1
       ty = threadidx%y-1
       tz = threadidx%z-1

! thread indices over block

       bx = tx+(ty*blockDim%x)+(tz*blockDim%x*blockDim%y)

! thread indices over grid

       i = (blockidx%x-1) * blockDim%x + tx
       j = (blockidx%y-1) * blockDim%y + ty

!Filter coefficients into shared memory 29*4*8 B=0.232 KB
if (-bx .le. 0 .and. -bx .ge. lowfil .and. bx .le. lupfil) then
filtersub(-bx,4) = filter(-bx,4)
filtersub(bx,4) = filter(bx,4)
endif

! (1/2) d^2/dx^2
if (i .ge. nfl1 .and. i .le. nfu1 .and. j .ge. nfl3 .and. j .le. nfu3) then
!ibyz_c into shared memory 4*4*2*4 B= 0.128 KB
!To avoid Non-stride-1 accesses we change the arrangement of the array.
ibxz_fsub(tx,ty,1) = ibxz_f(1,i,j)
ibxz_fsub(tx,ty,2) = ibxz_f(2,i,j)
call syncthreads()
do kb=ibxz_fsub(tx,ty,1),ibxz_fsub(tx,ty,2),32
t112=0.d0;t121=0.d0;t122=0.d0;t212=0.d0;t221=0.d0;t222=0.d0;t211=0.d0
!part of psig into shared memory (32+15)*4*4*2*8 B=12.032 KB
k = ibxz_fsub(tx,ty,2)-kb
if (tz .le. 13 .and. tz .le. k .and. tz+kb-14 .ge. 0) then
psigsub(tx,1,tz,2,ty,1) = psig(i,1,kb+tz-14,2,j,1)
end if
call syncthreads()
if (k .le. 13 .and. tz .le. 13) then
psigsub(tx,1,tz+k+1,2,ty,1) = psig(i,1,kb+tz-14+k+1,2,j,1)
else
psigsub(tx,1,tz+14,2,ty,1) = psig(i,1,kb+tz,2,j,1)
end if
call syncthreads()

!we only use the llow do loop when needed

if (tz .le. k) then
llow = nfl2-(tz+kb)+14

if (nfl2-kb+14 .ge. 1) then
 if (llow .ge. 1) then
 do l=llow,13
  t121=t121+psigsub(tx,1,tz+l,2,ty,1)*filtersub(-14+l,4)
 end do
 else
 do l=0,13
  t121=t121+psigsub(tx,1,tz+l,2,ty,1)*filtersub(-14+l,4)
 enddo
 end if
else
 do l=0,13
  t121=t121+psigsub(tx,1,tz+l,2,ty,1)*filtersub(-14+l,4)
 enddo
end if

end if

! Copy data from the end of the array to the beginning
call syncthreads()
if (tz .le. 13 .and. tz .le. k) then
psigsub(tx,1,tz,2,ty,1)=psigsub(tx,1,tz+14,2,ty,1)
end if
if (tz .le. 13 .and. k .ge. 14) then
psigsub(tx,1,tz+14,2,ty,1)=psigsub(tx,1,tz+28,2,ty,1)
end if
call syncthreads()
if (tz .ge. 28 .and. k .ge. 28) then 
psigsub(tx,1,tz,2,ty,1)=psigsub(tx,1,tz+14,2,ty,1) 
end if 
call syncthreads()


! Fill array with new data
if (tz .le. 13 .and. tz .le. k-1) then
psigsub(tx,1,tz+32,2,ty,1)=psig(i,1,tz+kb+32,2,j,1)
end if
call syncthreads()
if (tz .le. 13 .and. k .le. 31) then
psigsub(tx,1,tz+k+1,2,ty,1)=psig(i,1,tz+kb+k+1,2,j,1)
end if
call syncthreads()


lupfilsub = 14
!we only use the lhigh do loop when needed
if (tz .le. k) then
lhigh = nfu2-(tz+kb)
 
if (nfu2-(kb+blockDim%z) .le. lupfilsub-1) then
 if (lhigh .le. lupfilsub-1) then
 do l=0,lhigh
 t121=t121+psigsub(tx,1,tz+l,2,ty,1)*filtersub(l,4)
 end do
 else
 do l=0,lupfilsub
 t121=t121+psigsub(tx,1,tz+l,2,ty,1)*filtersub(l,4)
 end do
 end if
else
 do l=0,lupfilsub
 t121=t121+psigsub(tx,1,tz+l,2,ty,1)*filtersub(l,4)
 end do 
end if

ekinvec(i,kb+tz,j)=ekinvec(i,kb+tz,j)+t121*psigsub(tx,1,tz,2,ty,1)
y(i,1,kb+tz,2,j,1)=y(i,1,kb+tz,2,j,1)+t121
call syncthreads()

end if

   enddo
end if

    end subroutine convol_part15_kernel

attributes(global) subroutine convol_part16_kernel(psig,y,ekinvec,filter, ibxy_f,nfl1,nfl2,nfl3,nfu1,nfu2,nfu3,n1,n2,n3,lowfil,lupfil)
       real*8,device :: psig(0:n1,1:2,0:n2,1:2,0:n3,1:2),y(0:n1,1:2,0:n2,1:2,0:n3,1:2),ekinvec(0:n1,0:n2,0:n3),filter(lowfil:lupfil,4)
       integer*4,device :: ibxy_f(2,0:n1,0:n2)
       integer*4, value :: nfl1,nfl2,nfl3,nfu1,nfu2,nfu3,n1,n2,n3,lowfil,lupfil,l,llow,lhigh
       real*8 :: ekin,t111,t112,t121,t211,t122,t212,t221,t222
       integer*4 :: i, j, kb, k, bx, tx, ty, tz, ka

! submatrices are declared to be in CUDA shared memory

       integer*4, shared :: ibxy_fsub(0:0,0:7,1:2),lupfilsub
       real*8, shared :: filtersub(-14:14,1:4),psigsub(0:0,2:2,0:7,1:2,0:46,1:2)

! Start execution, first get thread indices

       tx = threadidx%x-1
       ty = threadidx%y-1
       tz = threadidx%z-1

! thread indices over block

       bx = tx+(ty*blockDim%x)+(tz*blockDim%x*blockDim%y)

! thread indices over grid

       i = (blockidx%x-1) * blockDim%x + tx
       j = (blockidx%y-1) * blockDim%y + ty

!Filter coefficients into shared memory 29*4*8 B=0.928 KB
if (-bx .le. 0 .and. -bx .ge. lowfil .and. bx .le. lupfil) then
filtersub(-bx,1) = filter(-bx,1)
filtersub(-bx,2) = filter(-bx,2)
filtersub(-bx,3) = filter(-bx,3)
filtersub(-bx,4) = filter(-bx,4)
filtersub(bx,1) = filter(bx,1)
filtersub(bx,2) = filter(bx,2)
filtersub(bx,3) = filter(bx,3)
filtersub(bx,4) = filter(bx,4)
endif

! (1/2) d^2/dx^2
if (i .ge. nfl1 .and. i .le. nfu1 .and. j .ge. nfl2 .and. j .le. nfu2) then
!ibxz_f into shared memory 4*4*2*4 B= 0.128 KB
!To avoid Non-stride-1 accesses we change the arrangement of the array.
ibxy_fsub(tx,ty,1) = ibxy_f(1,i,j)
ibxy_fsub(tx,ty,2) = ibxy_f(2,i,j)
call syncthreads()
do kb=ibxy_fsub(tx,ty,1),ibxy_fsub(tx,ty,2),32
t112=0.d0;t121=0.d0;t122=0.d0;t212=0.d0;t221=0.d0;t222=0.d0;t211=0.d0
!part of psig into shared memory (32+15)*1*8*4*8 B=12.032 KB
k = ibxy_fsub(tx,ty,2)-kb
if (tz .le. 13 .and. tz .le. k .and. tz+kb-14 .ge. 0) then
psigsub(tx,2,ty,1,tz,1) = psig(i,2,j,1,kb+tz-14,1)
psigsub(tx,2,ty,1,tz,2) = psig(i,2,j,1,kb+tz-14,2)
psigsub(tx,2,ty,2,tz,1) = psig(i,2,j,2,kb+tz-14,1) 
psigsub(tx,2,ty,2,tz,2) = psig(i,2,j,2,kb+tz-14,2) 
end if
call syncthreads()
if (k .le. 13 .and. tz .le. 13) then
psigsub(tx,2,ty,1,tz+k+1,1) = psig(i,2,j,1,kb+tz-14+k+1,1)
psigsub(tx,2,ty,1,tz+k+1,2) = psig(i,2,j,1,kb+tz-14+k+1,2)
psigsub(tx,2,ty,2,tz+k+1,1) = psig(i,2,j,2,kb+tz-14+k+1,1)
psigsub(tx,2,ty,2,tz+k+1,2) = psig(i,2,j,2,kb+tz-14+k+1,2)
else
psigsub(tx,2,ty,1,tz+14,1) = psig(i,2,j,1,kb+tz,1)
psigsub(tx,2,ty,1,tz+14,2) = psig(i,2,j,1,kb+tz,2)
psigsub(tx,2,ty,2,tz+14,1) = psig(i,2,j,2,kb+tz,1) 
psigsub(tx,2,ty,2,tz+14,2) = psig(i,2,j,2,kb+tz,2)
end if
call syncthreads()

!we only use the llow do loop when needed

if (tz .le. k) then
llow = nfl3-(tz+kb)+14

if (nfl3-kb+14 .ge. 1) then
 if (llow .ge. 1) then
 do l=llow,13
  t211=t211+psigsub(tx,2,ty,1,tz+l,1)*filtersub(-14+l,1)+psigsub(tx,2,ty,1,tz+l,2)*filtersub(-14+l,2)
  t212=t212+psigsub(tx,2,ty,1,tz+l,1)*filtersub(-14+l,3)+psigsub(tx,2,ty,1,tz+l,2)*filtersub(-14+l,4)
  t221=t221+psigsub(tx,2,ty,2,tz+l,1)*filtersub(-14+l,1)+psigsub(tx,2,ty,2,tz+l,2)*filtersub(-14+l,2)
  t222=t222+psigsub(tx,2,ty,2,tz+l,1)*filtersub(-14+l,3)+psigsub(tx,2,ty,2,tz+l,2)*filtersub(-14+l,4)
 end do
 else
 do l=0,13
  t211=t211+psigsub(tx,2,ty,1,tz+l,1)*filtersub(-14+l,1)+psigsub(tx,2,ty,1,tz+l,2)*filtersub(-14+l,2)
  t212=t212+psigsub(tx,2,ty,1,tz+l,1)*filtersub(-14+l,3)+psigsub(tx,2,ty,1,tz+l,2)*filtersub(-14+l,4)
  t221=t221+psigsub(tx,2,ty,2,tz+l,1)*filtersub(-14+l,1)+psigsub(tx,2,ty,2,tz+l,2)*filtersub(-14+l,2)
  t222=t222+psigsub(tx,2,ty,2,tz+l,1)*filtersub(-14+l,3)+psigsub(tx,2,ty,2,tz+l,2)*filtersub(-14+l,4)
 enddo
 end if
else
 do l=0,13
  t211=t211+psigsub(tx,2,ty,1,tz+l,1)*filtersub(-14+l,1)+psigsub(tx,2,ty,1,tz+l,2)*filtersub(-14+l,2)
  t212=t212+psigsub(tx,2,ty,1,tz+l,1)*filtersub(-14+l,3)+psigsub(tx,2,ty,1,tz+l,2)*filtersub(-14+l,4)
  t221=t221+psigsub(tx,2,ty,2,tz+l,1)*filtersub(-14+l,1)+psigsub(tx,2,ty,2,tz+l,2)*filtersub(-14+l,2)
  t222=t222+psigsub(tx,2,ty,2,tz+l,1)*filtersub(-14+l,3)+psigsub(tx,2,ty,2,tz+l,2)*filtersub(-14+l,4)
 enddo
end if

end if

! Copy data from the end of the array to the beginning
call syncthreads()
if (tz .le. 13 .and. tz .le. k) then
psigsub(tx,2,ty,1,tz,1)=psigsub(tx,2,ty,1,tz+14,1)
psigsub(tx,2,ty,1,tz,2)=psigsub(tx,2,ty,1,tz+14,2)
psigsub(tx,2,ty,2,tz,1)=psigsub(tx,2,ty,2,tz+14,1)
psigsub(tx,2,ty,2,tz,2)=psigsub(tx,2,ty,2,tz+14,2)
end if
if (tz .le. 13 .and. k .ge. 14) then
psigsub(tx,2,ty,1,tz+14,1)=psigsub(tx,2,ty,1,tz+28,1)
psigsub(tx,2,ty,1,tz+14,2)=psigsub(tx,2,ty,1,tz+28,2)
psigsub(tx,2,ty,2,tz+14,1)=psigsub(tx,2,ty,2,tz+28,1)
psigsub(tx,2,ty,2,tz+14,2)=psigsub(tx,2,ty,2,tz+28,2)
end if
call syncthreads()
if (tz .ge. 28 .and. k .ge. 28) then 
psigsub(tx,2,ty,1,tz,1)=psigsub(tx,2,ty,1,tz+14,1) 
psigsub(tx,2,ty,1,tz,2)=psigsub(tx,2,ty,1,tz+14,2)
psigsub(tx,2,ty,2,tz,1)=psigsub(tx,2,ty,2,tz+14,1) 
psigsub(tx,2,ty,2,tz,2)=psigsub(tx,2,ty,2,tz+14,2) 
end if 
call syncthreads()


! Fill array with new data
if (tz .le. 13 .and. tz .le. k-1) then
psigsub(tx,2,ty,1,tz+32,1)=psig(i,2,j,1,tz+kb+32,1)
psigsub(tx,2,ty,1,tz+32,2)=psig(i,2,j,1,tz+kb+32,2)
psigsub(tx,2,ty,2,tz+32,1)=psig(i,2,j,2,tz+kb+32,1) 
psigsub(tx,2,ty,2,tz+32,2)=psig(i,2,j,2,tz+kb+32,2)
end if
call syncthreads()
if (tz .le. 13 .and. k .le. 31) then
psigsub(tx,2,ty,1,tz+k+1,1)=psig(i,2,j,1,tz+kb+k+1,1)
psigsub(tx,2,ty,1,tz+k+1,2)=psig(i,2,j,1,tz+kb+k+1,2)
psigsub(tx,2,ty,2,tz+k+1,1)=psig(i,2,j,2,tz+kb+k+1,1)
psigsub(tx,2,ty,2,tz+k+1,2)=psig(i,2,j,2,tz+kb+k+1,2)
end if
call syncthreads()


lupfilsub = 14
!we only use the lhigh do loop when needed
if (tz .le. k) then
lhigh = nfu3-(tz+kb)
 
if (nfu3-(kb+blockDim%z) .le. lupfilsub-1) then
 if (lhigh .le. lupfilsub-1) then
 do l=0,lhigh
  t211=t211+psigsub(tx,2,ty,1,tz+l,1)*filtersub(l,1)+psigsub(tx,2,ty,1,tz+l,2)*filtersub(l,2)
  t212=t212+psigsub(tx,2,ty,1,tz+l,1)*filtersub(l,3)+psigsub(tx,2,ty,1,tz+l,2)*filtersub(l,4)
  t221=t221+psigsub(tx,2,ty,2,tz+l,1)*filtersub(l,1)+psigsub(tx,2,ty,2,tz+l,2)*filtersub(l,2)
  t222=t222+psigsub(tx,2,ty,2,tz+l,1)*filtersub(l,3)+psigsub(tx,2,ty,2,tz+l,2)*filtersub(l,4)
 end do
 else
 do l=0,lupfilsub
  t211=t211+psigsub(tx,2,ty,1,tz+l,1)*filtersub(l,1)+psigsub(tx,2,ty,1,tz+l,2)*filtersub(l,2)
  t212=t212+psigsub(tx,2,ty,1,tz+l,1)*filtersub(l,3)+psigsub(tx,2,ty,1,tz+l,2)*filtersub(l,4)
  t221=t221+psigsub(tx,2,ty,2,tz+l,1)*filtersub(l,1)+psigsub(tx,2,ty,2,tz+l,2)*filtersub(l,2)
  t222=t222+psigsub(tx,2,ty,2,tz+l,1)*filtersub(l,3)+psigsub(tx,2,ty,2,tz+l,2)*filtersub(l,4)
 end do
 end if
else
 do l=0,lupfilsub
  t211=t211+psigsub(tx,2,ty,1,tz+l,1)*filtersub(l,1)+psigsub(tx,2,ty,1,tz+l,2)*filtersub(l,2)
  t212=t212+psigsub(tx,2,ty,1,tz+l,1)*filtersub(l,3)+psigsub(tx,2,ty,1,tz+l,2)*filtersub(l,4)
  t221=t221+psigsub(tx,2,ty,2,tz+l,1)*filtersub(l,1)+psigsub(tx,2,ty,2,tz+l,2)*filtersub(l,2)
  t222=t222+psigsub(tx,2,ty,2,tz+l,1)*filtersub(l,3)+psigsub(tx,2,ty,2,tz+l,2)*filtersub(l,4)
 end do 
end if

ekinvec(i,j,kb+tz)=ekinvec(i,j,kb+tz)+t211*psigsub(tx,2,ty,1,tz,1)+t212*psigsub(tx,2,ty,1,tz,2)+t221*psigsub(tx,2,ty,2,tz,1)+t222*psigsub(tx,2,ty,2,tz,2)

y(i,2,j,1,kb+tz,1)=y(i,2,j,1,kb+tz,1)+t211
y(i,2,j,1,kb+tz,2)=y(i,2,j,1,kb+tz,2)+t212
y(i,2,j,2,kb+tz,1)=y(i,2,j,2,kb+tz,1)+t221
y(i,2,j,2,kb+tz,2)=y(i,2,j,2,kb+tz,2)+t222
call syncthreads()

end if

   enddo
end if

    end subroutine convol_part16_kernel


attributes(global) subroutine convol_part18_kernel(psig,y,ekinvec,filter,ibxy_f,nfl1,nfl2,nfl3,nfu1,nfu2,nfu3,n1,n2,n3,lowfil,lupfil)
       real*8,device :: psig(0:n1,1:2,0:n2,1:2,0:n3,1:2),y(0:n1,1:2,0:n2,1:2,0:n3,1:2),ekinvec(0:n1,0:n2,0:n3),filter(lowfil:lupfil,4)
       integer*4,device :: ibxy_f(2,0:n1,0:n2)
       integer*4, value :: nfl1,nfl2,nfl3,nfu1,nfu2,nfu3,n1,n2,n3,lowfil,lupfil,l,llow,lhigh
       real*8 :: ekin,t111,t112,t121,t211,t122,t212,t221,t222
       integer*4 :: i, j, kb, k, bx, tx, ty, tz, ka

! submatrices are declared to be in CUDA shared memory

       integer*4, shared :: ibxy_fsub(0:0,0:7,1:2),lupfilsub
       real*8, shared :: filtersub(-14:14,1:4),psigsub(0:0,1:1,0:7,2:2,0:46,1:2)

! Start execution, first get thread indices

       tx = threadidx%x-1
       ty = threadidx%y-1
       tz = threadidx%z-1

! thread indices over block

       bx = tx+(ty*blockDim%x)+(tz*blockDim%x*blockDim%y)

! thread indices over grid

       i = (blockidx%x-1) * blockDim%x + tx
       j = (blockidx%y-1) * blockDim%y + ty

!Filter coefficients into shared memory 29*4*8 B=0.928 KB
if (-bx .le. 0 .and. -bx .ge. lowfil .and. bx .le. lupfil) then
filtersub(-bx,1) = filter(-bx,1)
filtersub(-bx,2) = filter(-bx,2)
filtersub(-bx,3) = filter(-bx,3)
filtersub(-bx,4) = filter(-bx,4)
filtersub(bx,1) = filter(bx,1)
filtersub(bx,2) = filter(bx,2)
filtersub(bx,3) = filter(bx,3)
filtersub(bx,4) = filter(bx,4)
endif

! (1/2) d^2/dx^2
if (i .ge. nfl1 .and. i .le. nfu1 .and. j .ge. nfl2 .and. j .le. nfu2) then
!ibxz_f into shared memory 4*4*2*4 B= 0.128 KB
!To avoid Non-stride-1 accesses we change the arrangement of the array.
ibxy_fsub(tx,ty,1) = ibxy_f(1,i,j)
ibxy_fsub(tx,ty,2) = ibxy_f(2,i,j)
call syncthreads()
do kb=ibxy_fsub(tx,ty,1),ibxy_fsub(tx,ty,2),32
t112=0.d0;t121=0.d0;t122=0.d0;t212=0.d0;t221=0.d0;t222=0.d0;t211=0.d0
!part of psig into shared memory (32+15)*4*4*2*8 B=12.032 KB
k = ibxy_fsub(tx,ty,2)-kb
if (tz .le. 13 .and. tz .le. k .and. tz+kb-14 .ge. 0) then
psigsub(tx,1,ty,2,tz,1) = psig(i,1,j,2,kb+tz-14,1)
psigsub(tx,1,ty,2,tz,2) = psig(i,1,j,2,kb+tz-14,2)
end if
call syncthreads()
if (k .le. 13 .and. tz .le. 13) then
psigsub(tx,1,ty,2,tz+k+1,1) = psig(i,1,j,2,kb+tz-14+k+1,1)
psigsub(tx,1,ty,2,tz+k+1,2) = psig(i,1,j,2,kb+tz-14+k+1,2)
else
psigsub(tx,1,ty,2,tz+14,1) = psig(i,1,j,2,kb+tz,1)
psigsub(tx,1,ty,2,tz+14,2) = psig(i,1,j,2,kb+tz,2)
end if
call syncthreads()

!we only use the llow do loop when needed

if (tz .le. k) then
llow = nfl3-(tz+kb)+14

if (nfl3-kb+14 .ge. 1) then
 if (llow .ge. 1) then
 do l=llow,13
  t121=t121+psigsub(tx,1,ty,2,tz+l,1)*filtersub(-14+l,1)+psigsub(tx,1,ty,2,tz+l,2)*filtersub(-14+l,2)
  t122=t122+psigsub(tx,1,ty,2,tz+l,1)*filtersub(-14+l,3)+psigsub(tx,1,ty,2,tz+l,2)*filtersub(-14+l,4)
 end do
 else
 do l=0,13
  t121=t121+psigsub(tx,1,ty,2,tz+l,1)*filtersub(-14+l,1)+psigsub(tx,1,ty,2,tz+l,2)*filtersub(-14+l,2)
  t122=t122+psigsub(tx,1,ty,2,tz+l,1)*filtersub(-14+l,3)+psigsub(tx,1,ty,2,tz+l,2)*filtersub(-14+l,4)
 enddo
 end if
else
 do l=0,13
  t121=t121+psigsub(tx,1,ty,2,tz+l,1)*filtersub(-14+l,1)+psigsub(tx,1,ty,2,tz+l,2)*filtersub(-14+l,2)
  t122=t122+psigsub(tx,1,ty,2,tz+l,1)*filtersub(-14+l,3)+psigsub(tx,1,ty,2,tz+l,2)*filtersub(-14+l,4)
 enddo
end if

end if

! Copy data from the end of the array to the beginning
call syncthreads()
if (tz .le. 13 .and. tz .le. k) then
psigsub(tx,1,ty,2,tz,1)=psigsub(tx,1,ty,2,tz+14,1)
psigsub(tx,1,ty,2,tz,2)=psigsub(tx,1,ty,2,tz+14,2)
end if
if (tz .le. 13 .and. k .ge. 14) then
psigsub(tx,1,ty,2,tz+14,1)=psigsub(tx,1,ty,2,tz+28,1)
psigsub(tx,1,ty,2,tz+14,2)=psigsub(tx,1,ty,2,tz+28,2)
end if
call syncthreads()
if (tz .ge. 28 .and. k .ge. 28) then 
psigsub(tx,1,ty,2,tz,1)=psigsub(tx,1,ty,2,tz+14,1) 
psigsub(tx,1,ty,2,tz,2)=psigsub(tx,1,ty,2,tz+14,2)
end if 
call syncthreads()


! Fill array with new data
if (tz .le. 13 .and. tz .le. k-1) then
psigsub(tx,1,ty,2,tz+32,1)=psig(i,1,j,2,tz+kb+32,1)
psigsub(tx,1,ty,2,tz+32,2)=psig(i,1,j,2,tz+kb+32,2)
end if
call syncthreads()
if (tz .le. 13 .and. k .le. 31) then
psigsub(tx,1,ty,2,tz+k+1,1)=psig(i,1,j,2,tz+kb+k+1,1)
psigsub(tx,1,ty,2,tz+k+1,2)=psig(i,1,j,2,tz+kb+k+1,2)
end if
call syncthreads()


lupfilsub = 14
!we only use the lhigh do loop when needed
if (tz .le. k) then
lhigh = nfu3-(tz+kb)
 
if (nfu3-(kb+blockDim%z) .le. lupfilsub-1) then
 if (lhigh .le. lupfilsub-1) then
 do l=0,lhigh
  t121=t121+psigsub(tx,1,ty,2,tz+l,1)*filtersub(l,1)+psigsub(tx,1,ty,2,tz+l,2)*filtersub(l,2)
  t122=t122+psigsub(tx,1,ty,2,tz+l,1)*filtersub(l,3)+psigsub(tx,1,ty,2,tz+l,2)*filtersub(l,4)
 end do
 else
 do l=0,lupfilsub
  t121=t121+psigsub(tx,1,ty,2,tz+l,1)*filtersub(l,1)+psigsub(tx,1,ty,2,tz+l,2)*filtersub(l,2)
  t122=t122+psigsub(tx,1,ty,2,tz+l,1)*filtersub(l,3)+psigsub(tx,1,ty,2,tz+l,2)*filtersub(l,4)
  end do
 end if
else
 do l=0,lupfilsub
  t121=t121+psigsub(tx,1,ty,2,tz+l,1)*filtersub(l,1)+psigsub(tx,1,ty,2,tz+l,2)*filtersub(l,2)
  t122=t122+psigsub(tx,1,ty,2,tz+l,1)*filtersub(l,3)+psigsub(tx,1,ty,2,tz+l,2)*filtersub(l,4)
 end do 
end if

ekinvec(i,j,kb+tz)=ekinvec(i,j,kb+tz)+t121*psigsub(tx,1,ty,2,tz,1)+t122*psigsub(tx,1,ty,2,tz,2)

y(i,1,j,2,kb+tz,1)=y(i,1,j,2,kb+tz,1)+t121
y(i,1,j,2,kb+tz,2)=y(i,1,j,2,kb+tz,2)+t122
call syncthreads()

end if

   enddo
end if

    end subroutine convol_part18_kernel

attributes(global) subroutine convol_part19_kernel(psig,y,ekinvec,filter,ibxy_f,nfl1,nfl2,nfl3,nfu1,nfu2,nfu3,n1,n2,n3,lowfil,lupfil)
       real*8,device :: psig(0:n1,1:2,0:n2,1:2,0:n3,1:2),y(0:n1,1:2,0:n2,1:2,0:n3,1:2),ekinvec(0:n1,0:n2,0:n3),filter(lowfil:lupfil,4)
       integer*4,device :: ibxy_f(2,0:n1,0:n2)
       integer*4, value :: nfl1,nfl2,nfl3,nfu1,nfu2,nfu3,n1,n2,n3,lowfil,lupfil,l,llow,lhigh
       real*8 :: ekin,t111,t112,t121,t211,t122,t212,t221,t222
       integer*4 :: i, j, kb, k, bx, tx, ty, tz, ka

! submatrices are declared to be in CUDA shared memory

       integer*4, shared :: ibxy_fsub(0:0,0:7,1:2),lupfilsub
       real*8, shared :: filtersub(-14:14,4:4),psigsub(0:0,1:1,0:7,1:1,0:46,1:2)

! Start execution, first get thread indices

       tx = threadidx%x-1
       ty = threadidx%y-1
       tz = threadidx%z-1

! thread indices over block

       bx = tx+(ty*blockDim%x)+(tz*blockDim%x*blockDim%y)

! thread indices over grid

       i = (blockidx%x-1) * blockDim%x + tx
       j = (blockidx%y-1) * blockDim%y + ty

!Filter coefficients into shared memory 29*4*8 B=0.232 KB
if (-bx .le. 0 .and. -bx .ge. lowfil .and. bx .le. lupfil) then
filtersub(-bx,4) = filter(-bx,4)
filtersub(bx,4) = filter(bx,4)
endif

! (1/2) d^2/dx^2
if (i .ge. nfl1 .and. i .le. nfu1 .and. j .ge. nfl2 .and. j .le. nfu2) then
!ibxz_f into shared memory 4*4*2*4 B= 0.128 KB
!To avoid Non-stride-1 accesses we change the arrangement of the array.
ibxy_fsub(tx,ty,1) = ibxy_f(1,i,j)
ibxy_fsub(tx,ty,2) = ibxy_f(2,i,j)
call syncthreads()
do kb=ibxy_fsub(tx,ty,1),ibxy_fsub(tx,ty,2),32
t112=0.d0;t121=0.d0;t122=0.d0;t212=0.d0;t221=0.d0;t222=0.d0;t211=0.d0
!part of psig into shared memory (32+15)*4*4*2*8 B=12.032 KB
k = ibxy_fsub(tx,ty,2)-kb
if (tz .le. 13 .and. tz .le. k .and. tz+kb-14 .ge. 0) then
psigsub(tx,1,ty,1,tz,2) = psig(i,1,j,1,kb+tz-14,2)
end if
call syncthreads()
if (k .le. 13 .and. tz .le. 13) then
psigsub(tx,1,ty,1,tz+k+1,2) = psig(i,1,j,1,kb+tz-14+k+1,2)
else
psigsub(tx,1,ty,1,tz+14,2) = psig(i,1,j,1,kb+tz,2)
end if
call syncthreads()

!we only use the llow do loop when needed

if (tz .le. k) then
llow = nfl3-(tz+kb)+14

if (nfl3-kb+14 .ge. 1) then
 if (llow .ge. 1) then
 do l=llow,13
  t112=t112+psigsub(tx,1,ty,1,tz+l,2)*filtersub(-14+l,4)
 end do
 else
 do l=0,13
  t112=t112+psigsub(tx,1,ty,1,tz+l,2)*filtersub(-14+l,4)
 enddo
 end if
else
 do l=0,13
  t112=t112+psigsub(tx,1,ty,1,tz+l,2)*filtersub(-14+l,4)
 enddo
end if

end if

! Copy data from the end of the array to the beginning
call syncthreads()
if (tz .le. 13 .and. tz .le. k) then
psigsub(tx,1,ty,1,tz,2)=psigsub(tx,1,ty,1,tz+14,2)
end if
if (tz .le. 13 .and. k .ge. 14) then
psigsub(tx,1,ty,1,tz+14,2)=psigsub(tx,1,ty,1,tz+28,2)
end if
call syncthreads()
if (tz .ge. 28 .and. k .ge. 28) then 
psigsub(tx,1,ty,1,tz,2)=psigsub(tx,1,ty,1,tz+14,2) 
end if 
call syncthreads()


! Fill array with new data
if (tz .le. 13 .and. tz .le. k-1) then
psigsub(tx,1,ty,1,tz+32,2)=psig(i,1,j,1,tz+kb+32,2)
end if
call syncthreads()
if (tz .le. 13 .and. k .le. 31) then
psigsub(tx,1,ty,1,tz+k+1,2)=psig(i,1,j,1,tz+kb+k+1,2)
end if
call syncthreads()


lupfilsub = 14
!we only use the lhigh do loop when needed
if (tz .le. k) then
lhigh = nfu3-(tz+kb)
 
if (nfu3-(kb+blockDim%z) .le. lupfilsub-1) then
 if (lhigh .le. lupfilsub-1) then
 do l=0,lhigh
  t112=t112+psigsub(tx,1,ty,1,tz+l,2)*filtersub(l,4)
 end do
 else
 do l=0,lupfilsub
  t112=t112+psigsub(tx,1,ty,1,tz+l,2)*filtersub(l,4)
  end do
 end if
else
 do l=0,lupfilsub
  t112=t112+psigsub(tx,1,ty,1,tz+l,2)*filtersub(l,4)
 end do 
end if

ekinvec(i,j,kb+tz)=ekinvec(i,j,kb+tz)+t112*psigsub(tx,1,ty,1,tz,2)

y(i,1,j,1,kb+tz,2)=y(i,1,j,1,kb+tz,2)+t112
call syncthreads()

end if

   enddo
end if

    end subroutine convol_part19_kernel

     attributes(global) subroutine ones_kernel(onesvec,n1,n2,n3)
       real*8,device :: onesvec(0:(n1+1)*(n2+1)*(n3+1)-1)
       integer*4, value :: n1,n2,n3
       integer :: i, j, kb, k, bx, bz, tx, ty, tz

! Start execution, first get thread indices
 
       tx = threadidx%x-1

! thread indices over grid

       i = (blockidx%x-1) * blockDim%x + tx

if (i .le. (n1+1)*(n2+1)*(n3+1)-1) then
onesvec(i) = 1.d0
end if

end subroutine ones_kernel

