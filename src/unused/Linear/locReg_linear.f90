!> Count for each orbital and each process the number of overlapping orbitals.
subroutine determine_overlap_from_descriptors(iproc, nproc, orbs, orbsig, lzd, lzdig, op, comon)
use module_base
use module_types
use module_interfaces, except_this_one => determine_overlap_from_descriptors
implicit none

! Calling arguments
integer,intent(in) :: iproc, nproc
type(orbitals_data),intent(in) :: orbs, orbsig
type(local_zone_descriptors),intent(in) :: lzd, lzdig
type(overlapParameters),intent(inout) :: op
type(p2pComms),intent(inout) :: comon
! Local variables
integer :: jproc, iorb, jorb, ioverlapMPI, ioverlaporb, ilr, jlr, ilrold
!integer :: is1, ie1, is2, ie2, is3, ie3
!integer :: js1, je1, js2, je2, js3, je3
integer :: iiorb, istat, iall, noverlaps, ierr, noverlapsmaxp
!logical :: ovrlpx, ovrlpy, ovrlpz
logical :: isoverlap
integer :: i1, i2, ii
integer :: onseg
logical,dimension(:,:),allocatable :: overlapMatrix
integer,dimension(:),allocatable :: noverlapsarr, displs, recvcnts, overlaps_comon
integer,dimension(:,:),allocatable :: overlaps_op
integer,dimension(:,:,:),allocatable :: overlaps_nseg
!integer,dimension(:,:,:),allocatable :: iseglist, jseglist
character(len=*),parameter :: subname='determine_overlap_from_descriptors'

allocate(overlapMatrix(orbsig%norb,maxval(orbs%norb_par(:,0))), stat=istat)
call memocc(istat, overlapMatrix, 'overlapMatrix', subname)
allocate(noverlapsarr(orbs%norbp), stat=istat)
call memocc(istat, noverlapsarr, 'noverlapsarr', subname)
allocate(displs(0:nproc-1), stat=istat)
call memocc(istat, displs, 'displs', subname)
allocate(recvcnts(0:nproc-1), stat=istat)
call memocc(istat, recvcnts, 'recvcnts', subname)
allocate(overlaps_nseg(orbsig%norb,orbs%norbp,2), stat=istat)
call memocc(istat, overlaps_nseg, 'overlaps_nseg', subname)

    overlapMatrix=.false.
    overlaps_nseg = 0
    ioverlapMPI=0 ! counts the overlaps for the given MPI process.
    ilrold=-1
    do iorb=1,orbs%norbp
        ioverlaporb=0 ! counts the overlaps for the given orbital.
        iiorb=orbs%isorb+iorb
        ilr=orbs%inWhichLocreg(iiorb)
!        call get_start_and_end_indices(lzd%llr(ilr), is1, ie1, is2, ie2, is3, ie3)
        do jorb=1,orbsig%norb
            jlr=orbsig%inWhichLocreg(jorb)
            call check_overlap_cubic_periodic(lzd%Glr,lzd%llr(ilr),lzdig%llr(jlr),isoverlap)
!            call get_start_and_end_indices(lzd%llr(jlr), js1, je1, js2, je2, js3, je3)
!            ovrlpx = ( is1<=je1 .and. ie1>=js1 )
!            ovrlpy = ( is2<=je2 .and. ie2>=js2 )
!            ovrlpz = ( is3<=je3 .and. ie3>=js3 )
!            if(ovrlpx .and. ovrlpy .and. ovrlpz) then
             !write(*,*) 'second: iproc, ilr, jlr, isoverlap', iproc, ilr, jlr, isoverlap
             if(isoverlap) then
                ! From the viewpoint of the box boundaries, an overlap between ilr and jlr is possible.
                ! Now explicitely check whether there is an overlap by using the descriptors.
                !!call overlapbox_from_descriptors(lzd%llr(ilr)%d%n1, lzd%llr(ilr)%d%n2, lzd%llr(ilr)%d%n3, &
                !!     lzd%llr(jlr)%d%n1, lzd%llr(jlr)%d%n2, lzd%llr(jlr)%d%n3, &
                !!     lzd%glr%d%n1, lzd%glr%d%n2, lzd%glr%d%n3, &
                !!     lzd%llr(ilr)%ns1, lzd%llr(ilr)%ns2, lzd%llr(ilr)%ns3, &
                !!     lzd%llr(jlr)%ns1, lzd%llr(jlr)%ns2, lzd%llr(jlr)%ns3, &
                !!     lzd%glr%ns1, lzd%glr%ns2, lzd%glr%ns3, &
                !!     lzd%llr(ilr)%wfd%nseg_c, lzd%llr(jlr)%wfd%nseg_c, &
                !!     lzd%llr(ilr)%wfd%keygloc, lzd%llr(ilr)%wfd%keyvloc, lzd%llr(jlr)%wfd%keygloc, lzd%llr(jlr)%wfd%keyvloc, &
                !!     n1_ovrlp, n2_ovrlp, n3_ovrlp, ns1_ovrlp, ns2_ovrlp, ns3_ovrlp, nseg_ovrlp)
                call check_overlap_from_descriptors_periodic(lzd%llr(ilr)%wfd%nseg_c, lzdig%llr(jlr)%wfd%nseg_c,&
                     lzd%llr(ilr)%wfd%keyglob, lzdig%llr(jlr)%wfd%keyglob, &
                     isoverlap, onseg)
                if(isoverlap) then
                    ! There is really an overlap
                    overlapMatrix(jorb,iorb)=.true.
                    ioverlaporb=ioverlaporb+1
                    overlaps_nseg(ioverlaporb,iorb,1)=onseg
                    if(ilr/=ilrold) then
                        ! if ilr==ilrold, we are in the same localization region, so the MPI prosess
                        ! would get the same orbitals again. Therefore the counter is not increased
                        ! in that case.
                        ioverlapMPI=ioverlapMPI+1
                    end if
                else
                    overlapMatrix(jorb,iorb)=.false.
                end if
             else
                overlapMatrix(jorb,iorb)=.false.
             end if
        end do
        noverlapsarr(iorb)=ioverlaporb
        ilrold=ilr
    end do
    !comon%noverlaps(jproc)=ioverlapMPI
    noverlaps=ioverlapMPI

!call mpi_allreduce(overlapMatrix, orbs%norb*maxval(orbs%norb_par(:,0))*nproc, mpi_sum bigdft_mpi%mpi_comm, ierr)

! Communicate op%noverlaps and comon%noverlaps
    if (nproc > 1) then
       call mpi_allgatherv(noverlapsarr, orbs%norbp, mpi_integer, op%noverlaps, orbs%norb_par, &
            orbs%isorb_par, mpi_integer, bigdft_mpi%mpi_comm, ierr)
    else
       call vcopy(orbs%norb,noverlapsarr(1),1,op%noverlaps(1),1)
    end if
    do jproc=0,nproc-1
       recvcnts(jproc)=1
       displs(jproc)=jproc
    end do
    if (nproc > 1) then
       call mpi_allgatherv(noverlaps, 1, mpi_integer, comon%noverlaps, recvcnts, &
            displs, mpi_integer, bigdft_mpi%mpi_comm, ierr)
    else
       comon%noverlaps=noverlaps
    end if


allocate(op%overlaps(maxval(op%noverlaps),orbs%norb), stat=istat)
call memocc(istat, op%overlaps, 'op%overlaps', subname)
!!allocate(comon%overlaps(maxval(comon%noverlaps),0:nproc-1), stat=istat)
!!call memocc(istat, comon%overlaps, 'comon%overlaps', subname)
allocate(overlaps_op(maxval(op%noverlaps),orbs%norbp), stat=istat)
call memocc(istat, overlaps_op, 'overlaps_op', subname)
if(orbs%norbp>0) call to_zero(maxval(op%noverlaps)*orbs%norbp,overlaps_op(1,1))
allocate(overlaps_comon(comon%noverlaps(iproc)), stat=istat)
call memocc(istat, overlaps_comon, 'overlaps_comon', subname)

! Now we know how many overlaps have to be calculated, so determine which orbital overlaps
! with which one. This is essentially the same loop as above, but we use the array 'overlapMatrix'
! which indicates the overlaps.

! Initialize to some value which will never be used.
op%overlaps=-1
!!comon%overlaps=-1

iiorb=0
ioverlapMPI=0 ! counts the overlaps for the given MPI process.
ilrold=-1
do iorb=1,orbs%norbp
    ioverlaporb=0 ! counts the overlaps for the given orbital.
    iiorb=orbs%isorb+iorb
    ilr=orbs%inWhichLocreg(iiorb)
    do jorb=1,orbsig%norb
        jlr=orbsig%inWhichLocreg(jorb)
        if(overlapMatrix(jorb,iorb)) then
            ioverlaporb=ioverlaporb+1
            ! Determine the number of segments of the fine grid in the overlap
            call check_overlap_from_descriptors_periodic(lzd%llr(ilr)%wfd%nseg_c, lzdig%llr(jlr)%wfd%nseg_f,&
                 lzd%llr(ilr)%wfd%keyglob(1,1), lzdig%llr(jlr)%wfd%keyglob(1,1+lzdig%llr(jlr)%wfd%nseg_c), &
                 isoverlap, overlaps_nseg(ioverlaporb,iorb,2))
            !op%overlaps(ioverlaporb,iiorb)=jorb
            overlaps_op(ioverlaporb,iorb)=jorb
            if(ilr/=ilrold) then
                ! if ilr==ilrold, we are in th same localization region, so the MPI prosess
                ! would get the same orbitals again. Therefore the counter is not increased
                ! in that case.
                ioverlapMPI=ioverlapMPI+1
                !comon%overlaps(ioverlapMPI,iproc)=jorb
                overlaps_comon(ioverlapMPI)=jorb
            end if
        end if
    end do 
    ilrold=ilr
end do

! Allocate the overlap wavefunctions_descriptors
! and copy the nseg_c
noverlapsmaxp=maxval(op%noverlaps(orbs%isorb+1:orbs%isorb+orbs%norbp))
allocate(op%wfd_overlap(noverlapsmaxp,orbs%norbp), stat=istat)
do i2=1,orbs%norbp
    do i1=1,noverlapsmaxp
        call nullify_wavefunctions_descriptors(op%wfd_overlap(i1,i2))
        op%wfd_overlap(i1,i2)%nseg_c = overlaps_nseg(i1,i2,1)
        op%wfd_overlap(i1,i2)%nseg_f = overlaps_nseg(i1,i2,2)
        call allocate_wfd(op%wfd_overlap(i1,i2),subname)
    end do
end do

!Now redo the loop for the keygs
iiorb=0
do iorb=1,orbs%norbp
    ioverlaporb=0 ! counts the overlaps for the given orbital.
    iiorb=orbs%isorb+iorb
    ilr=orbs%inWhichLocreg(iiorb)
    do jorb=1,orbsig%norb
        jlr=orbsig%inWhichLocreg(jorb)
        if(overlapMatrix(jorb,iorb)) then
            ioverlaporb=ioverlaporb+1
           ! Determine the keyglob, keyvglob, nvctr of the coarse grid
           call get_overlap_from_descriptors_periodic(lzd%llr(ilr)%wfd%nseg_c, lzdig%llr(jlr)%wfd%nseg_c, &
                lzd%llr(ilr)%wfd%keyglob(1,1), lzdig%llr(jlr)%wfd%keyglob(1,1),  &
                .true.,op%wfd_overlap(ioverlaporb,iorb)%nseg_c, op%wfd_overlap(ioverlaporb,iorb)%nvctr_c,&
                op%wfd_overlap(ioverlaporb,iorb)%keyglob(1,1), op%wfd_overlap(ioverlaporb,iorb)%keyvglob(1))
           ! Determine the keyglob, keyvglob, nvctr of the fine grid
           if(op%wfd_overlap(ioverlaporb,iorb)%nseg_f > 0) then
              call get_overlap_from_descriptors_periodic(lzd%llr(ilr)%wfd%nseg_c, lzdig%llr(jlr)%wfd%nseg_f, &
                   lzd%llr(ilr)%wfd%keyglob(1,1), lzdig%llr(jlr)%wfd%keyglob(1,1+lzdig%llr(jlr)%wfd%nseg_c),  &
                   .true.,op%wfd_overlap(ioverlaporb,iorb)%nseg_f, op%wfd_overlap(ioverlaporb,iorb)%nvctr_f,&
                   op%wfd_overlap(ioverlaporb,iorb)%keyglob(1,op%wfd_overlap(ioverlaporb,iorb)%nseg_c+1), &
                   op%wfd_overlap(ioverlaporb,iorb)%keyvglob(op%wfd_overlap(ioverlaporb,iorb)%nseg_c+1))
           else
              op%wfd_overlap(ioverlaporb,iorb)%nvctr_f = 0
           end if
        end if
    end do 
end do


displs(0)=0
recvcnts(0)=comon%noverlaps(0)
do jproc=1,nproc-1
    recvcnts(jproc)=comon%noverlaps(jproc)
    displs(jproc)=displs(jproc-1)+recvcnts(jproc-1)
end do
!!if (nproc > 1) then
!!   call mpi_allgatherv(overlaps_comon, comon%noverlaps(iproc), mpi_integer, comon%overlaps, recvcnts, &
!!        displs, mpi_integer, bigdft_mpi%mpi_comm, ierr)
!!else
!!   call vcopy(comon%noverlaps(iproc),overlaps_comon(1),1,comon%overlaps(1,0),1)
!!end if
ii=maxval(op%noverlaps)
displs(0)=0
recvcnts(0)=ii*orbs%norb_par(0,0)
do jproc=1,nproc-1
    recvcnts(jproc)=ii*orbs%norb_par(jproc,0)
    displs(jproc)=displs(jproc-1)+recvcnts(jproc-1)
end do
if (nproc > 1) then
   call mpi_allgatherv(overlaps_op, ii*orbs%norbp, mpi_integer, op%overlaps, recvcnts, &
        displs, mpi_integer, bigdft_mpi%mpi_comm, ierr)
else
   call vcopy(ii*orbs%norbp,overlaps_op(1,1),1,op%overlaps(1,1),1)
end if


iall=-product(shape(overlapMatrix))*kind(overlapMatrix)
deallocate(overlapMatrix, stat=istat)
call memocc(istat, iall, 'overlapMatrix', subname)

iall=-product(shape(noverlapsarr))*kind(noverlapsarr)
deallocate(noverlapsarr, stat=istat)
call memocc(istat, iall, 'noverlapsarr', subname)

iall=-product(shape(overlaps_op))*kind(overlaps_op)
deallocate(overlaps_op, stat=istat)
call memocc(istat, iall, 'overlaps_op', subname)

iall=-product(shape(overlaps_comon))*kind(overlaps_comon)
deallocate(overlaps_comon, stat=istat)
call memocc(istat, iall, 'overlaps_comon', subname)

iall=-product(shape(displs))*kind(displs)
deallocate(displs, stat=istat)
call memocc(istat, iall, 'displs', subname)

iall=-product(shape(recvcnts))*kind(recvcnts)
deallocate(recvcnts, stat=istat)
call memocc(istat, iall, 'recvcnts', subname)

iall=-product(shape(overlaps_nseg))*kind(overlaps_nseg)
deallocate(overlaps_nseg, stat=istat)
call memocc(istat, iall, 'overlaps_nseg', subname)

end subroutine determine_overlap_from_descriptors
