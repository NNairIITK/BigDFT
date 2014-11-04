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



subroutine determine_overlapdescriptors_from_descriptors(llr_i, llr_j, glr, olr)
use module_base
use module_types
implicit none

! Calling arguments
type(locreg_descriptors),intent(in) :: llr_i, llr_j, glr
type(locreg_descriptors),intent(out) :: olr

! Local variables
integer :: n1_ovrlp, n2_ovrlp, n3_ovrlp, ns1_ovrlp, ns2_ovrlp, ns3_ovrlp
character(len=*),parameter :: subname='determine_overlapdescriptors_from_descriptors'



! Determine the values describing the localization region.
! First coarse region.
call overlapbox_from_descriptors(llr_i%d%n1, llr_i%d%n2, llr_i%d%n3, &
     llr_j%d%n1, llr_j%d%n2, llr_j%d%n3, &
     glr%d%n1, glr%d%n2, glr%d%n3, &
     llr_i%ns1, llr_i%ns2, llr_i%ns3, &
     llr_j%ns1, llr_j%ns2, llr_j%ns3, &
     glr%ns1, glr%ns2, glr%ns3, &
     llr_i%wfd%nseg_c, llr_j%wfd%nseg_c, &
     llr_i%wfd%keygloc, llr_i%wfd%keyvloc, llr_j%wfd%keygloc, llr_j%wfd%keyvloc, &
     olr%d%n1, olr%d%n2, olr%d%n3, olr%ns1, olr%ns2, olr%ns3, olr%wfd%nseg_c)

! Now the fine region.
call overlapbox_from_descriptors(llr_i%d%n1, llr_i%d%n2, llr_i%d%n3, &
     llr_j%d%n1, llr_j%d%n2, llr_j%d%n3, &
     glr%d%n1, glr%d%n2, glr%d%n3, &
     llr_i%ns1, llr_i%ns2, llr_i%ns3, &
     llr_j%ns1, llr_j%ns2, llr_j%ns3, &
     glr%ns1, glr%ns2, glr%ns3, &
     llr_i%wfd%nseg_f, llr_j%wfd%nseg_f, &
     llr_i%wfd%keygloc(1,llr_i%wfd%nseg_c+min(1,llr_i%wfd%nseg_f)), llr_i%wfd%keyvloc(llr_i%wfd%nseg_c+min(1,llr_i%wfd%nseg_f)), &
     llr_j%wfd%keygloc(1,llr_j%wfd%nseg_c+min(1,llr_j%wfd%nseg_f)), llr_j%wfd%keyvloc(llr_j%wfd%nseg_c+min(1,llr_j%wfd%nseg_f)), &
     n1_ovrlp, n2_ovrlp, n3_ovrlp, ns1_ovrlp, ns2_ovrlp, ns3_ovrlp, olr%wfd%nseg_f)
 
     ! Determine the boundary for the fine part.
     ! ns1_ovrlp etc is in global coordinates, but olr%d%nfl1 etc is in local coordinates, so correct this.
     olr%d%nfl1=ns1_ovrlp-olr%ns1
     olr%d%nfu1=olr%d%nfl1+n1_ovrlp
     olr%d%nfl2=ns2_ovrlp-olr%ns2
     olr%d%nfu2=olr%d%nfl2+n2_ovrlp
     olr%d%nfl3=ns2_ovrlp-olr%ns2
     olr%d%nfu3=olr%d%nfl3+n3_ovrlp

if(llr_i%geocode/=llr_j%geocode) then
    write(*,*) 'ERROR: llr_i%geocode/=llr_j%geocode'
    stop
end if
olr%geocode=llr_i%geocode

! Dimensions for interpolating scaling function grid
olr%d%n1i=2*olr%d%n1+31
olr%d%n2i=2*olr%d%n2+31
olr%d%n3i=2*olr%d%n3+31



! Allocate the descriptor structures
call allocate_wfd(olr%wfd,subname)


! some checks
call check_overlapregion(glr, llr_i, llr_j, olr)



! Fill the descriptors for the coarse part.
call overlapdescriptors_from_descriptors(llr_i%d%n1, llr_i%d%n2, llr_i%d%n3, &
     llr_j%d%n1, llr_j%d%n2, llr_j%d%n3, &
     glr%d%n1, glr%d%n2, glr%d%n3, &
     olr%d%n1, olr%d%n2, olr%d%n3, &
     llr_i%ns1, llr_i%ns2, llr_i%ns3, &
     llr_j%ns1, llr_j%ns2, llr_j%ns3, &
     glr%ns1, glr%ns2, glr%ns3, &
     olr%ns1, olr%ns2, olr%ns3, &
     llr_i%wfd%nseg_c, llr_j%wfd%nseg_c, olr%wfd%nseg_c, &
     llr_i%wfd%keygloc, llr_i%wfd%keyvloc, llr_j%wfd%keygloc, llr_j%wfd%keyvloc, &
     olr%wfd%keygloc, olr%wfd%keyvloc, &
     olr%wfd%nvctr_c)

! Fill the descriptors for the fine part.
call overlapdescriptors_from_descriptors(llr_i%d%n1, llr_i%d%n2, llr_i%d%n3, &
     llr_j%d%n1, llr_j%d%n2, llr_j%d%n3, &
     glr%d%n1, glr%d%n2, glr%d%n3, &
     olr%d%n1, olr%d%n2, olr%d%n3, &
     llr_i%ns1, llr_i%ns2, llr_i%ns3, &
     llr_j%ns1, llr_j%ns2, llr_j%ns3, &
     glr%ns1, glr%ns2, glr%ns3, &
     olr%ns1, olr%ns2, olr%ns3, &
     llr_i%wfd%nseg_f, llr_j%wfd%nseg_f, olr%wfd%nseg_f, &
     llr_i%wfd%keygloc(1,llr_i%wfd%nseg_c+min(1,llr_i%wfd%nseg_f)), llr_i%wfd%keyvloc(llr_i%wfd%nseg_c+min(1,llr_i%wfd%nseg_f)), &
     llr_j%wfd%keygloc(1,llr_j%wfd%nseg_c+min(1,llr_j%wfd%nseg_f)), llr_j%wfd%keyvloc(llr_j%wfd%nseg_c+min(1,llr_j%wfd%nseg_f)), &
     olr%wfd%keygloc(1,olr%wfd%nseg_c+min(1,olr%wfd%nseg_f)), olr%wfd%keyvloc(olr%wfd%nseg_c+min(1,olr%wfd%nseg_f)), &
     olr%wfd%nvctr_f)



end subroutine determine_overlapdescriptors_from_descriptors



!> Creates the overlap descriptor from two input overlap descriptors. The dimension of the overlap box (i.e.
!! ns1_k,ns2_k,ns3_k (starting indices) and n1_k,n2_k,n3_k (length) have to be determined earlier (using
!! the subroutine overlapbox_from_descriptors)).
!! Calling arguments: *_i refers to overlap region i (input)
!!                    *_j refers to overlap region j (input)
!!                    *_g refers to the global region (input)
!!                    *_k refers to the overlap region (input/output)
!> SHOULD NOW CHANGE THIS BECAUSE GLOBAL COORDINATES ARE NOW ALWAYS AVAILABLE
subroutine overlapdescriptors_from_descriptors(n1_i, n2_i, n3_i, n1_j, n2_j, n3_j, n1_g, n2_g, n3_g, n1_k, n2_k, n3_k, &
           ns1_i, ns2_i, ns3_i, ns1_j, ns2_j, ns3_j, ns1_g, ns2_g, ns3_g, ns1_k, ns2_k, ns3_k, &
           nseg_i, nseg_j, nseg_k, &
           keyg_i, keyv_i, keyg_j, keyv_j, keyg_k, keyv_k, nvctr_k)
use module_base
use module_types
implicit none

! Calling arguments
integer,intent(in) :: n1_i, n2_i, n3_i, n1_j, n2_j, n3_j, n1_g, n2_g, n3_g, n1_k, n2_k, n3_k
integer,intent(in) :: ns1_i, ns2_i, ns3_i, ns1_j, ns2_j, ns3_j, ns1_g, ns2_g, ns3_g, ns1_k, ns2_k, ns3_k
integer,intent(in) :: nseg_i, nseg_j, nseg_k
integer,dimension(2,nseg_i),intent(in) :: keyg_i
integer,dimension(nseg_i),intent(in) :: keyv_i
integer,dimension(2,nseg_j),intent(in) :: keyg_j
integer,dimension(nseg_j),intent(in) :: keyv_j
integer,dimension(2,nseg_k),intent(out) :: keyg_k
integer,dimension(nseg_k),intent(out) :: keyv_k
integer,intent(out) :: nvctr_k

! Local variables
integer :: iseg, jseg, kseg, knvctr, istart, jstart, kstart, istartg, jstartg, kstartg
integer :: iend, jend, kend, iendg, jendg, kendg, transform_index
character(len=1) :: increase

! Initialize some counters
iseg=min(1,nseg_i)
jseg=min(1,nseg_j)
kseg=0
knvctr=0
nvctr_k=0

! Quick return if possible
if(nseg_i==0 .or. nseg_j==0) return

segment_loop: do

    ! Starting point in local coordinates
    istart=keyg_i(1,iseg)
    jstart=keyg_j(1,jseg)

    ! Get the global counterparts
    istartg=transform_index(istart, n1_i, n2_i, n3_i, n1_g, n2_g, n3_g, ns1_i-ns1_g, ns2_i-ns2_g, ns3_i-ns3_g)
    jstartg=transform_index(jstart, n1_j, n2_j, n3_j, n1_g, n2_g, n3_g, ns1_j-ns1_g, ns2_j-ns2_g, ns3_j-ns3_g)

    ! Ending point in local coordinates
    iend=keyg_i(2,iseg)
    jend=keyg_j(2,jseg)

    ! Get the global counterparts
    iendg=transform_index(iend, n1_i, n2_i, n3_i, n1_g, n2_g, n3_g, ns1_i-ns1_g, ns2_i-ns2_g, ns3_i-ns3_g)
    jendg=transform_index(jend, n1_j, n2_j, n3_j, n1_g, n2_g, n3_g, ns1_j-ns1_g, ns2_j-ns2_g, ns3_j-ns3_g)

    ! Determine starting and ending point of the common segment in global coordinates.
    kstartg=max(istartg,jstartg)
    kendg=min(iendg,jendg)
    if((iendg<=jendg .and. iseg<nseg_i) .or. jseg==nseg_j) then
        increase='i'
    else if(jseg<nseg_j) then
        increase='j'
    end if

    ! Check whether this common segment has a non-zero length
    if(kendg-kstartg+1>0) then
        kseg=kseg+1

        ! Transform the starting and ending point to the overlap localization region.
        kstart=transform_index(kstartg, n1_g, n2_g, n3_g, n1_k, n2_k, n3_k, ns1_g-ns1_k, ns2_g-ns2_k, ns3_g-ns3_k)
        kend=transform_index(kendg, n1_g, n2_g, n3_g, n1_k, n2_k, n3_k, ns1_g-ns1_k, ns2_g-ns2_k, ns3_g-ns3_k)

        ! Assign the values to the descriptors
        keyg_k(1,kseg)=kstart
        keyg_k(2,kseg)=kend
        keyv_k(kseg)=knvctr+1

        knvctr=knvctr+kendg-kstartg+1
    end if

    ! Check whether all segments of both localization regions have been processed
    if(iseg>=nseg_i .and. jseg>=nseg_j) exit segment_loop

    ! Increase the segment index
    if(increase=='i') then
        iseg=iseg+1
    else if(increase=='j') then
        jseg=jseg+1
    end if

end do segment_loop

nvctr_k=knvctr


end subroutine overlapdescriptors_from_descriptors

subroutine check_overlapregion(glr, llr_i, llr_j, olr)

  use module_base
  use module_types
  implicit none

  ! Calling arguments
  type(locreg_descriptors),intent(in) :: glr, llr_i, llr_j, olr

  if(olr%ns1<glr%ns1) then
      write(*,'(a,2(i0,a))') 'ERROR: olr%ns1 = ', olr%ns1, ' < ', glr%ns1, '= glr%ns1'
      stop
  end if
  if(olr%ns2<glr%ns2) then
      write(*,'(a,2(i0,a))') 'ERROR: olr%ns2 = ', olr%ns2, ' < ', glr%ns2, '= glr%ns2'
      stop
  end if
  if(olr%ns3<glr%ns3) then
      write(*,'(a,2(i0,a))') 'ERROR: olr%ns3 = ', olr%ns3, ' < ', glr%ns3, '= glr%ns3'
      stop
  end if
  if(olr%ns1+olr%d%n1>glr%ns1+glr%d%n1) then
      write(*,'(a,2(i0,a))') 'ERROR: olr%ns1+olr%d%n1 = ', olr%ns1+olr%d%n1, ' < ', glr%ns1+glr%d%n1, '= glr%ns1+glr%d%n1'
      stop
  end if
  if(olr%ns2+olr%d%n2>glr%ns2+glr%d%n2) then
      write(*,'(a,2(i0,a))') 'ERROR: olr%ns2+olr%d%n2 = ', olr%ns2+olr%d%n2, ' < ', glr%ns2+glr%d%n2, '= glr%ns2+glr%d%n2'
      stop
  end if
  if(olr%ns3+olr%d%n3>glr%ns3+glr%d%n3) then
      write(*,'(a,2(i0,a))') 'ERROR: olr%ns3+olr%d%n3 = ', olr%ns3+olr%d%n3, ' < ', glr%ns3+glr%d%n3, '= glr%ns3+glr%d%n3'
      stop
  end if
  
  if(olr%ns1<llr_i%ns1) then
      write(*,'(a,2(i0,a))') 'ERROR: olr%ns1 = ', olr%ns1, ' < ', llr_i%ns1, '= llr_i%ns1'
      stop
  end if
  if(olr%ns2<llr_i%ns2) then
      write(*,'(a,2(i0,a))') 'ERROR: olr%ns2 = ', olr%ns2, ' < ', llr_i%ns2, '= llr_i%ns2'
      stop
  end if
  if(olr%ns3<llr_i%ns3) then
      write(*,'(a,2(i0,a))') 'ERROR: olr%ns3 = ', olr%ns3, ' < ', llr_i%ns3, '= llr_i%ns3'
      stop
  end if
  if(olr%ns1+olr%d%n1>llr_i%ns1+llr_i%d%n1) then
      write(*,'(a,2(i0,a))') 'ERROR: olr%ns1+olr%d%n1 = ', olr%ns1+olr%d%n1, ' < ', llr_i%ns1+llr_i%d%n1, '= llr_i%ns1+llr_i%d%n1'
      stop
  end if
  if(olr%ns2+olr%d%n2>llr_i%ns2+llr_i%d%n2) then
      write(*,'(a,2(i0,a))') 'ERROR: olr%ns2+olr%d%n2 = ', olr%ns2+olr%d%n2, ' < ', llr_i%ns2+llr_i%d%n2, '= llr_i%ns2+llr_i%d%n2'
      stop
  end if
  if(olr%ns3+olr%d%n3>llr_i%ns3+llr_i%d%n3) then
      write(*,'(a,2(i0,a))') 'ERROR: olr%ns3+olr%d%n3 = ', olr%ns3+olr%d%n3, ' < ', llr_i%ns3+llr_i%d%n3, '= llr_i%ns3+llr_i%d%n3'
      stop
  end if
  
  if(olr%ns1<llr_j%ns1) then
      write(*,'(a,2(i0,a))') 'ERROR: olr%ns1 = ', olr%ns1, ' < ', llr_j%ns1, '= llr_j%ns1'
      stop
  end if
  if(olr%ns2<llr_j%ns2) then
      write(*,'(a,2(i0,a))') 'ERROR: olr%ns2 = ', olr%ns2, ' < ', llr_j%ns2, '= llr_j%ns2'
      stop
  end if
  if(olr%ns3<llr_j%ns3) then
      write(*,'(a,2(i0,a))') 'ERROR: olr%ns3 = ', olr%ns3, ' < ', llr_j%ns3, '= llr_j%ns3'
      stop
  end if
  if(olr%ns1+olr%d%n1>llr_j%ns1+llr_j%d%n1) then
      write(*,'(a,2(i0,a))') 'ERROR: olr%ns1+olr%d%n1 = ', olr%ns1+olr%d%n1, ' < ', llr_j%ns1+llr_j%d%n1, '= llr_j%ns1+llr_j%d%n1'
      stop
  end if
  if(olr%ns2+olr%d%n2>llr_j%ns2+llr_j%d%n2) then
      write(*,'(a,2(i0,a))') 'ERROR: olr%ns2+olr%d%n2 = ', olr%ns2+olr%d%n2, ' < ', llr_j%ns2+llr_j%d%n2, '= llr_j%ns2+llr_j%d%n2'
      stop
  end if
  if(olr%ns3+olr%d%n3>llr_j%ns3+llr_j%d%n3) then
      write(*,'(a,2(i0,a))') 'ERROR: olr%ns3+olr%d%n3 = ', olr%ns3+olr%d%n3, ' < ', llr_j%ns3+llr_j%d%n3, '= llr_j%ns3+llr_j%d%n3'
      stop
  end if

end subroutine check_overlapregion



!> Gives the dimensions of the overlap box resulting from the overlap of two wavefunction descriptors and
!> the number of segments of the resulting overlap descriptor.
!> Calling arguments: *_i refers to overlap region i (input)
!>                    *_j refers to overlap region j (input)
!>                    *_g refers to the global region (input)
!>                    *_k refers to the overlap region (output)
subroutine overlapbox_from_descriptors(n1_i, n2_i, n3_i, n1_j, n2_j, n3_j, n1_g, n2_g, n3_g, &
           ns1_i, ns2_i, ns3_i, ns1_j, ns2_j, ns3_j, ns1_g, ns2_g, ns3_g, &
           nseg_i, nseg_j, &
           keyg_i, keyv_i, keyg_j, keyv_j, &
           n1_k, n2_k, n3_k, ns1_k, ns2_k, ns3_k, nseg_k)
  use module_base
  use module_types
  implicit none

  ! Calling arguments
  integer,intent(in) :: n1_i, n2_i, n3_i, n1_j, n2_j, n3_j, n1_g, n2_g, n3_g
  integer,intent(in) :: ns1_i, ns2_i, ns3_i, ns1_j, ns2_j, ns3_j, ns1_g, ns2_g, ns3_g
  integer :: nseg_i, nseg_j
  integer,dimension(2,nseg_i),intent(in) :: keyg_i
  integer,dimension(nseg_i),intent(in) :: keyv_i
  integer,dimension(2,nseg_j),intent(in) :: keyg_j
  integer,dimension(nseg_j),intent(in) :: keyv_j
  integer,intent(out) :: n1_k, n2_k, n3_k, ns1_k, ns2_k, ns3_k, nseg_k

  ! Local variables
  integer :: iseg, jseg, istart, jstart, istartg, jstartg, kstartg
  integer :: iend, jend, iendg, jendg, kendg, transform_index
  integer :: kxs, kys, kzs, kxe, kye, kze, kxemax, kyemax, kzemax
  character(len=1) :: increase


  ! Initialize the return values such that they represent a box with no volume
  ns1_k=ns1_g+n1_g+1
  ns2_k=ns2_g+n2_g+1
  ns3_k=ns3_g+n3_g+1
  n1_k=0
  n2_k=0
  n3_k=0
  nseg_k=0

  ! Quick return if possible
  if(nseg_i==0 .or. nseg_j==0) return

  ! Initialize some counters
  iseg=min(1,nseg_i)
  jseg=min(1,nseg_j)
  kxemax=0
  kyemax=0
  kzemax=0


  segment_loop: do

      ! Starting point in local coordinates
      istart=keyg_i(1,iseg)
      jstart=keyg_j(1,jseg)

      ! Get the global counterparts
      istartg=transform_index(istart, n1_i, n2_i, n3_i, n1_g, n2_g, n3_g, ns1_i-ns1_g, ns2_i-ns2_g, ns3_i-ns3_g)
      jstartg=transform_index(jstart, n1_j, n2_j, n3_j, n1_g, n2_g, n3_g, ns1_j-ns1_g, ns2_j-ns2_g, ns3_j-ns3_g)

      ! Ending point in local coordinates
      iend=keyg_i(2,iseg)
      jend=keyg_j(2,jseg)

      ! Get the global counterparts
      iendg=transform_index(iend, n1_i, n2_i, n3_i, n1_g, n2_g, n3_g, ns1_i-ns1_g, ns2_i-ns2_g, ns3_i-ns3_g)
      jendg=transform_index(jend, n1_j, n2_j, n3_j, n1_g, n2_g, n3_g, ns1_j-ns1_g, ns2_j-ns2_g, ns3_j-ns3_g)

      ! Determine starting and ending point of the common segment in global coordinates.
      kstartg=max(istartg,jstartg)
      kendg=min(iendg,jendg)

      ! Determine which segment counter should be increased.
      if((iendg<=jendg .and. iseg<nseg_i) .or. jseg==nseg_j) then
          increase='i'
      else if(jseg<nseg_j) then
          increase='j'
      end if

      ! Check whether this common segment has a non-zero length
      if(kendg-kstartg+1>0) then
          nseg_k=nseg_k+1

          ! Get the global coordinates of this segment
          call get_coordinates(kstartg, n1_g, n2_g, n3_g, kxs, kys, kzs)
          call get_coordinates(kendg, n1_g, n2_g, n3_g, kxe, kye, kze)

          ! Check whether this segment enlarges the overlap box
          if(kxs<ns1_k) ns1_k=kxs
          if(kys<ns2_k) ns2_k=kys
          if(kzs<ns3_k) ns3_k=kzs
          if(kxe>kxemax) kxemax=kxe
          if(kye>kyemax) kyemax=kye
          if(kze>kzemax) kzemax=kze

      end if

      ! Check whether all segments of both localization regions have been processed.
      if(iseg>=nseg_i .and. jseg>=nseg_j) exit segment_loop

      ! Increase the segment index
      if(increase=='i') then
          iseg=iseg+1
      else if(increase=='j') then
          jseg=jseg+1
      end if

  end do segment_loop

  ! n1_k etc is the length of the segment, but kxemax etc is the end position of the segment, 
  ! therefore subtract the starting position
  n1_k=kxemax-ns1_k
  n2_k=kyemax-ns2_k
  n3_k=kzemax-ns3_k

end subroutine overlapbox_from_descriptors



subroutine get_overlap_from_descriptors_periodic(nseg_i, nseg_j, keyg_i, keyg_j,  &
           isoverlap,onseg, onvctr, keyglob, keyvglob)
  use module_base
  use module_types
  implicit none
  ! Calling arguments
  integer :: nseg_i, nseg_j
  integer,dimension(2,nseg_i),intent(in) :: keyg_i
  integer,dimension(2,nseg_j),intent(in) :: keyg_j
  logical, intent(in) :: isoverlap
  integer, intent(in) :: onseg
  integer, intent(out) :: onvctr
  integer, dimension(2,max(onseg,1)),intent(out) :: keyglob
  integer, dimension(max(onseg,1)), intent(out) :: keyvglob
  ! Local variables
  integer :: iseg, jseg, knvctr, istart, jstart, kstartg
  integer :: iend, jend, kendg, nseg_k

  if(.not. isoverlap) then
  !initialize the variable
  onvctr = 0
  return
  end if

  nseg_k = 0
  knvctr = 1
  iseg = 1
  jseg= 1
  segment_loop: do

      ! Starting point already in global coordinates
      istart=keyg_i(1,iseg)
      jstart=keyg_j(1,jseg)

      ! Ending point already in global coordinates
      iend=keyg_i(2,iseg)
      jend=keyg_j(2,jseg)

      ! Determine starting and ending point of the common segment in global coordinates.
      kstartg=max(istart,jstart)
      kendg=min(iend,jend)

      ! Check whether this common segment has a non-zero length
      ! and assign the keys
      if(kendg-kstartg+1>0) then
          nseg_k=nseg_k+1
          keyglob(1,nseg_k) = kstartg
          keyglob(2,nseg_k) = kendg
          keyvglob(nseg_k) =  knvctr
          knvctr = knvctr + kendg-kstartg+1
      end if

      ! Check whether all segments of both localization regions have been processed.
      if(iseg>=nseg_i .and. jseg>=nseg_j) exit segment_loop

      ! Increase the segment index
      if((iend<=jend .and. iseg<nseg_i) .or. jseg==nseg_j) then
          iseg=iseg+1
      else if(jseg<nseg_j) then
          jseg=jseg+1
      end if

  end do segment_loop

  onvctr = knvctr-1
  !check if everything matches
  if(onseg .ne. nseg_k) then
    print *,'onseg ',onseg,' nseg_k ',nseg_k
    stop 'get_overlap_from_descriptors_periodic: number of segments not right'
  end if

end subroutine get_overlap_from_descriptors_periodic


!> Transform an index from localization region A to localization region B.
!! Calling arguments:
!!   ist: index to transform (with respcet to coordinate system of locreg A)
!!   n1a, n2a, n3a:             box size of locreg A
!!   n1b, n2b, n3b:             box size of locreg B
!!   nshift1, nshift2, nshift3: nsa-nsb, where nsa,nsb are the starting points of the boxes of A,B (for all 3 dimensions)
function transform_index(ist, n1a, n2a, n3a, n1b, n2b, n3b, nshift1, nshift2, nshift3)
  implicit none

  ! Calling arguments
  integer,intent(in) :: ist, n1a, n2a, n3a, n1b, n2b, n3b, nshift1, nshift2, nshift3
  integer :: transform_index

  ! Local variables
  integer :: ii, ix, iy, iz, ixg, iyg, izg, istg

    ! Get the coordinates with respect to localization region A
    ii = ist - 1
    iz = ii / ((n1a+1) * (n2a+1))
    ii = ii - iz * ((n1a+1) * (n2a+1))
    iy = ii / (n1a+1)
    ix = ii - iy * (n1a+1)

    ! Transform ix, iy, iz to the coordinates with respect to localization region B
    ! WARNING: MODULO ONLY CORRECT WHEN COMPARING TO GLR
    ! FOR NOW THIS FUNCTION ALWAYS DOES THIS... SHOULD DISAPEAR NOW THE KEYGLOB IS DEFINED
  !!  izg = modulo(iz + nshift3, n3b)
  !!  iyg = modulo(iy + nshift2, n2b)
  !!  ixg = modulo(ix + nshift1, n1b)
    izg = iz + nshift3
    iyg = iy + nshift2
    ixg = ix + nshift1

    ! Transform ist to its counterpart in the coordinate system of B
    istg = izg*(n1b+1)*(n2b+1) + iyg*(n1b+1) + ixg + 1
    
    ! Assign istg to the value that is passed back.
    transform_index=istg

end function transform_index


!> Get the coordinates of ist with respect to its localization region
!! Calling arguments:
!!  ist          index for which coordinates shall be calculated
!!  n1, n2, n3   box sizes
!!  ix, iy, iz   coordinates of ist
subroutine get_coordinates(ist, n1, n2, n3, ix, iy, iz)
  implicit none

  ! Calling arguments
  integer,intent(in) :: ist, n1, n2, n3
  integer,intent(out) :: ix, iy, iz

  ! Local variable
  integer :: ii

    ! Get the coordinates ix, iy, iz
    ii = ist - 1
    !ii = ist
    iz = ii / ((n1+1) * (n2+1))
    ii = ii - iz * ((n1+1) * (n2+1))
    iy = ii / (n1+1)
    ix = ii - iy * (n1+1)

  end subroutine get_coordinates

subroutine transform_ISFcoordinates(direction,iin1,iin2,iin3,Glr,Llr,iout1,iout2,iout3,ishift1, ishift2, ishift3)
  use module_types
  implicit none
  integer, intent(in) :: direction                      !< integer specifying the direction of the transformation (0 from local to global, 1 from global to local)
  integer, intent(in) :: iin1, iin2, iin3               !< the input coordinate
  type(locreg_descriptors), intent(in) :: Glr           !< global region descriptors
  type(locreg_descriptors), intent(in) :: Llr           !< Localization regiondescriptors
  integer, intent(out) :: iout1, iout2, iout3           !< ouput coordinates
  integer, intent(out) :: ishift1, ishift2, ishift3     !< shift between coordinates: in_coord + shift = out_coord

  !This routines supposes that the specified coordinate is part of both regions.
  if(direction==0) then
  !from local to global
  iout1 = iin1 -(Glr%nsi1 - Llr%nsi1)
  ishift1 = Llr%nsi1 - Glr%nsi1
  iout2 = iin2 -(Glr%nsi2 - Llr%nsi2)
  ishift2 = Llr%nsi2 - Glr%nsi2
  iout3 = iin3 -(Glr%nsi3 - Llr%nsi3)
  ishift3 = Llr%nsi3 - Glr%nsi3
  if(iout1 > Glr%nsi1 + Glr%d%n1i .and. Glr%geocode/='F') then
     iout1 = modulo(iout1,Glr%d%n1i+1)+Glr%nsi1
     ishift1 = ishift1 - Glr%d%n1i
  end if
  if(iout2 > Glr%nsi2 + Glr%d%n2i .and. Glr%geocode=='P') then
     iout2 = modulo(iout2,Glr%d%n2i+1)+Glr%nsi2
     ishift2 = ishift2 - Glr%d%n2i
  end if
  if(iout3 > Glr%nsi3 + Glr%d%n3i .and. Glr%geocode/='F') then
     iout3 = modulo(iout3,Glr%d%n3i+1)+Glr%nsi3
     ishift3 = ishift3 - Glr%d%n3i
  end if

  else if(direction==1) then
  !from global to local
  iout1 = iin1 -(Llr%nsi1 - Glr%nsi1)
  ishift1 = Glr%nsi1 - Llr%nsi1
  iout2 = iin2 -(Llr%nsi2 - Glr%nsi2)
  ishift2 = Glr%nsi2 - Llr%nsi2
  iout3 = iin3 -(Llr%nsi3 - Glr%nsi3)
  ishift3 = Glr%nsi3 - Llr%nsi3
  if(iout1 < 0 .and. Glr%geocode/='F') then
     iout1 = iin1 + Glr%d%n1i - Llr%nsi1
     ishift1 = ishift1 + Glr%d%n1i 
  end if
  if(iout2 < 0 .and. Glr%geocode=='P') then
     iout2 = iin2 + Glr%d%n2i - Llr%nsi2
     ishift2 = ishift2 + Glr%d%n2i 
  end if
  if(iout3 < 0 .and. Glr%geocode/='F') then
     iout3 = iin3 + Glr%d%n3i - Llr%nsi3
     ishift3 = ishift3 + Glr%d%n3i 
  end if

  else
  STOP 'transform_ISFcoordinates: wrong descriptor, should be 0 or 1'
  end if

end subroutine transform_ISFcoordinates


subroutine fracture_periodic_zone_ISF(nzones,Glr,Llr,outofzone,astart,aend)

  use module_base
  use module_types

  implicit none

  ! Subroutine Scalar Arguments
  integer,intent(in) :: nzones
  type(locreg_descriptors),intent(in) :: Glr  ! Global grid descriptor
  type(locreg_descriptors),intent(in) :: Llr  ! Localization grid descriptors 
  !Subroutine Array Arguments
  integer,dimension(3),intent(in) :: outofzone  ! array indicating the directions in which the locreg exceeds the Glr
  integer,dimension(3,nzones),intent(out) :: astart !
  integer,dimension(3,nzones),intent(out) :: aend !
  !local variables
  integer :: ii,index,jj
  integer,dimension(3) :: alrs,alre,Gend,Gstart,period

! Start and end of Global region
  Gstart(1) = Glr%nsi1
  Gstart(2) = Glr%nsi2
  Gstart(3) = Glr%nsi3
  Gend(1) = Glr%nsi1 + Glr%d%n1i
  Gend(2) = Glr%nsi2 + Glr%d%n2i
  Gend(3) = Glr%nsi3 + Glr%d%n3i

! Periodicity of the system
  period(1) = Glr%d%n1i + 1
  period(2) = Glr%d%n2i + 1
  period(3) = Glr%d%n3i + 1

! Start and end of local region
  alrs(1) = Llr%nsi1
  alrs(2) = Llr%nsi2
  alrs(3) = Llr%nsi3
  alre(1) = Llr%nsi1 + Llr%d%n1i
  alre(2) = Llr%nsi2 + Llr%d%n2i
  alre(3) = Llr%nsi3 + Llr%d%n3i

  if(outofzone(1) <= 0 .and. outofzone(2) <= 0 .and. outofzone(3) <= 0 ) then !Nothing to do
     astart(1,1) = alrs(1)
     astart(2,1) = alrs(2)
     astart(3,1) = alrs(3)
     aend(1,1) = alre(1)
     aend(2,1) = alre(2)
     aend(3,1) = alre(3)
     return
  end if

!assign the first zone (necessarily without shift) and initiliaze the others
  do ii=1,3
     astart(ii,:) = alrs(ii)
     aend(ii,:) = min(Gend(ii),alre(ii))
  end do

!assign the other zones
  index = 2
  do ii=1,2
     if(outofzone(ii) > 0) then    !Translation: X,Y
        astart(ii,index) =  Gstart(ii)
        aend(ii,index) = modulo(alre(ii),period(ii))
        index = index + 1
     end if
     do jj=ii+1,2
        if(outofzone(ii) > 0 .and. outofzone(jj) > 0) then  !Translation: X+Y
           astart(ii,index) = Gstart(ii)
           astart(jj,index) = Gstart(jj)
           aend(ii,index) = modulo(alre(ii),period(ii))
           aend(jj,index) = modulo(alre(jj),period(jj))
           index = index + 1
        end if
     end do
  end do

END SUBROUTINE fracture_periodic_zone_ISF

