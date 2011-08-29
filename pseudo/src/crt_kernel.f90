attributes(global) subroutine zeros_kernel(psig,n1,n2,n3,s)
       real*8,device :: psig(0:(n1+1)*(n2+1)*(n3+1)*8-1)
       integer, value :: n1,n2,n3,s
       integer :: i, j, kb, k, tx, ty, tz

! Start execution, first get thread indices

       tx = threadidx%x-1

! thread indices over grid

       i = (blockidx%x-1) * blockDim%x + tx

!Set the psig vector to 0

if (i .le. (n1+1)*(n2+1)*(n3+1)*8-1-s) then
psig(s+i) = 0.d0
end if

end subroutine zeros_kernel

attributes(global) subroutine crt_part1_kernel(allwproj, psig,nl1_c,nl2_c,nl3_c,nu1_c,nu2_c,nu3_c, nl1_f, nl2_f, nl3_f, nu1_f, nu2_f, nu3_f,n1,n2,n3)
       real*8,device :: psig(0:n1,1:2,0:n2,1:2,0:n3,1:2),allwproj(0:n1+n2+n3+2,2) 
       integer, value :: nu1_c,nu2_c,nu3_c,nl3_f,nu3_f,nl2_f,nu2_f,nl1_f,nu1_f,nl1_c,nl2_c,nl3_c,n1,n2,n3

       integer :: i, j, kb, k, tx, ty, tz

! submatrices are declared to be in CUDA shared memory

       real*8, shared :: wprojxsub(0:15,1:2), wprojysub(0:15,1:2), wprojzsub(0:1,1:2)

! Start execution, first get my thread indices

       tx = threadidx%x-1
       ty = threadidx%y-1
       tz = threadidx%z-1

! thread indices over grid

       i = (blockidx%x-1) * blockDim%x + tx
       j = (blockidx%y-1) * blockDim%y + ty

!The two loops using shared memory
if (i .le. nu1_c .and. j .le. nu2_c) then
 wprojxsub(tx,1) = allwproj(i,1)
 wprojysub(ty,1) = allwproj(j+n1+1,1)
 do kb=0,nu3_c-tz,2
 wprojzsub(tz,1) = allwproj(tz+kb+n1+n2+2,1)
 call syncthreads()
  psig(i,1,j,1,kb+tz,1)=wprojxsub(tx,1)*wprojysub(ty,1)*wprojzsub(tz,1)
 call syncthreads()
 end do
end if

!Second loop
if (i .le. nu1_f .and. j .le. nu2_f .and. i .ge. nl1_f .and. j .ge. nl2_f) then
 wprojxsub(tx,2) = allwproj(i,2)
 wprojysub(ty,2) = allwproj(j+n1+1,2)
do kb=nl3_f,nu3_f-tz,2
wprojzsub(tz,1) = allwproj(tz+kb+n1+n2+2,1)
wprojzsub(tz,2) = allwproj(tz+kb+n1+n2+2,2)
 call syncthreads()
 psig(i,2,j,1,kb+tz,1)=wprojxsub(tx,2)*wprojysub(ty,1)*wprojzsub(tz,1)
 psig(i,1,j,2,kb+tz,1)=wprojxsub(tx,1)*wprojysub(ty,2)*wprojzsub(tz,1)
 psig(i,1,j,1,kb+tz,2)=wprojxsub(tx,1)*wprojysub(ty,1)*wprojzsub(tz,2)
 psig(i,1,j,2,kb+tz,2)=wprojxsub(tx,1)*wprojysub(ty,2)*wprojzsub(tz,2)
 psig(i,2,j,2,kb+tz,1)=wprojxsub(tx,2)*wprojysub(ty,2)*wprojzsub(tz,1)
 psig(i,2,j,1,kb+tz,2)=wprojxsub(tx,2)*wprojysub(ty,1)*wprojzsub(tz,2)
 psig(i,2,j,2,kb+tz,2)=wprojxsub(tx,2)*wprojysub(ty,2)*wprojzsub(tz,2)
 call syncthreads()
end do
end if

    end subroutine crt_part1_kernel



attributes(global) subroutine crt_part2_kernel( allwproj,psig,nl1_c,nl2_c,nl3_c,nu1_c,nu2_c,nu3_c, nl1_f, nl2_f, nl3_f, nu1_f, nu2_f, nu3_f,n1,n2,n3)
       real*8,device :: psig(0:n1,1:2,0:n2,1:2,0:n3,1:2),allwproj(0:n1+n2+n3+2,2)
       integer, value :: nu1_c,nu2_c,nu3_c,nl3_f,nu3_f,nl2_f,nu2_f,nl1_f,nu1_f,nl1_c,nl2_c,nl3_c,n1,n2,n3

       integer :: i, j, kb, k, tx, ty, tz

! submatrices are declared to be in CUDA shared memory

       real*8, shared :: wprojxsub(0:15,1:2), wprojysub(0:15,1:2), wprojzsub(0:1,1:2)

! Start execution, first get thread indices

       tx = threadidx%x-1
       ty = threadidx%y-1
       tz = threadidx%z-1

! thread indices over grid

       i = (blockidx%x-1) * blockDim%x + tx
       j = (blockidx%y-1) * blockDim%y + ty

!The two loops using shared memory
if (i .le. nu1_c .and. j .le. nu2_c) then
 wprojxsub(tx,1) = allwproj(i,1)
 wprojysub(ty,1) = allwproj(j+n1+1,1)
 do kb=0,nu3_c-tz,2
 wprojzsub(tz,1) = allwproj(tz+kb+n1+n2+2,1)
 call syncthreads()
  psig(i,1,j,1,kb+tz,1)=psig(i,1,j,1,kb+tz,1)+wprojxsub(tx,1)*wprojysub(ty,1)*wprojzsub(tz,1)
 call syncthreads()
 end do
end if

!Second loop
if (i .le. nu1_f .and. j .le. nu2_f .and. i .ge. nl1_f .and. j .ge. nl2_f) then
 wprojxsub(tx,2) = allwproj(i,2)
 wprojysub(ty,2) = allwproj(j+n1+1,2)
do kb=nl3_f,nu3_f-tz,2
wprojzsub(tz,1) = allwproj(tz+kb+n1+n2+2,1)
wprojzsub(tz,2) = allwproj(tz+kb+n1+n2+2,2)
 call syncthreads()
 psig(i,2,j,1,kb+tz,1)=psig(i,2,j,1,kb+tz,1)+wprojxsub(tx,2)*wprojysub(ty,1)*wprojzsub(tz,1)
 psig(i,1,j,2,kb+tz,1)=psig(i,1,j,2,kb+tz,1)+wprojxsub(tx,1)*wprojysub(ty,2)*wprojzsub(tz,1)
 psig(i,1,j,1,kb+tz,2)=psig(i,1,j,1,kb+tz,2)+wprojxsub(tx,1)*wprojysub(ty,1)*wprojzsub(tz,2)
 psig(i,1,j,2,kb+tz,2)=psig(i,1,j,2,kb+tz,2)+wprojxsub(tx,1)*wprojysub(ty,2)*wprojzsub(tz,2)
 psig(i,2,j,2,kb+tz,1)=psig(i,2,j,2,kb+tz,1)+wprojxsub(tx,2)*wprojysub(ty,2)*wprojzsub(tz,1)
 psig(i,2,j,1,kb+tz,2)=psig(i,2,j,1,kb+tz,2)+wprojxsub(tx,2)*wprojysub(ty,1)*wprojzsub(tz,2)
 psig(i,2,j,2,kb+tz,2)=psig(i,2,j,2,kb+tz,2)+wprojxsub(tx,2)*wprojysub(ty,2)*wprojzsub(tz,2)
call syncthreads()
end do
end if

    end subroutine crt_part2_kernel



