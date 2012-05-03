subroutine plotOrbitals(iproc, orbs, Glr, phi, nat, rxyz, hxh, hyh, hzh, it)
!
! Plots the orbitals
!
use module_base
use module_types
implicit none

! Calling arguments
integer:: iproc
type(orbitals_data), intent(inout) :: orbs
type(locreg_descriptors), intent(in) :: Glr
real(8),dimension((Glr%wfd%nvctr_c+7*Glr%wfd%nvctr_f)*orbs%nspinor*orbs%norbp):: phi
integer:: nat
real(8),dimension(3,nat):: rxyz
real(8):: hxh, hyh, hzh
integer:: it

integer:: ix, iy, iz, ix0, iy0, iz0, iiAt, jj, iorb, i1, i2, i3, istart, ii, istat
integer:: unit1, unit2, unit3
real(8),dimension(:),allocatable:: phir
type(workarr_sumrho) :: w
character(len=10):: c1, c2, c3
character(len=50):: file1, file2, file3

allocate(phir(Glr%d%n1i*Glr%d%n2i*Glr%d%n3i), stat=istat)

call initialize_work_arrays_sumrho(Glr,w)

istart=0

unit1=10*iproc+7
unit2=10*iproc+8
unit3=10*iproc+9

!write(*,*) 'write, orbs%nbasisp', orbs%norbp
    orbLoop: do iorb=1,orbs%norbp
        phir=0.d0
        call daub_to_isf(Glr,w,phi(istart+1),phir(1))
        iiAt=orbs%inwhichlocreg(orbs%isorb+iorb)
        ix0=nint(rxyz(1,iiAt)/hxh)
        iy0=nint(rxyz(2,iiAt)/hyh)
        iz0=nint(rxyz(3,iiAt)/hzh)

        jj=0
        write(c1,'(i0)') iproc
        write(c2,'(i0)') iorb
        write(c3,'(i0)') it
        file1='orbs_'//trim(c1)//'_'//trim(c2)//'_'//trim(c3)//'_x'
        file2='orbs_'//trim(c1)//'_'//trim(c2)//'_'//trim(c3)//'_y'
        file3='orbs_'//trim(c1)//'_'//trim(c2)//'_'//trim(c3)//'_z'
        open(unit=unit1, file=trim(file1))
        open(unit=unit2, file=trim(file2))
        open(unit=unit3, file=trim(file3))
        do i3=1,Glr%d%n3i
            do i2=1,Glr%d%n2i
                do i1=1,Glr%d%n1i
                   jj=jj+1
                   ! z component of point jj
                   iz=jj/(Glr%d%n2i*Glr%d%n1i)
                   ! Subtract the 'lower' xy layers
                   ii=jj-iz*(Glr%d%n2i*Glr%d%n1i)
                   ! y component of point jj
                   iy=ii/Glr%d%n1i
                   ! Subtract the 'lower' y rows
                   ii=ii-iy*Glr%d%n1i
                   ! x component
                   ix=ii
!if(phir(jj)>1.d0) write(*,'(a,3i7,es15.6)') 'WARNING: ix, iy, iz, phir(jj)', ix, iy, iz, phir(jj)
                   if(iy==ix0 .and. iz==iz0) write(unit1,*) ix, phir(jj)
                   ! Write along y-axis
                   if(ix==ix0 .and. iz==iz0) write(unit2,*) iy, phir(jj)
                   ! Write along z-axis
                   if(ix==ix0 .and. iy==iy0) write(unit3,*) iz, phir(jj)


                end do
            end do
        end do
        close(unit=unit1)
        close(unit=unit2)
        close(unit=unit3)

        istart=istart+(Glr%wfd%nvctr_c+7*Glr%wfd%nvctr_f)*orbs%nspinor

    end do orbLoop

call deallocate_work_arrays_sumrho(w)
deallocate(phir, stat=istat)


end subroutine plotOrbitals



subroutine compressMatrix(norb, mad, mat, lmat)
  use module_base
  use module_types
  implicit none
  
  ! Calling arguments
  integer,intent(in):: norb
  type(matrixDescriptors),intent(in):: mad
  real(8),dimension(norb**2),intent(in):: mat
  real(8),dimension(mad%nvctr),intent(out):: lmat
  
  ! Local variables
  integer:: iseg, jj, jorb, iiorb, jjorb
  
  
  jj=0
  do iseg=1,mad%nseg
      do jorb=mad%keyg(1,iseg),mad%keyg(2,iseg)
          jj=jj+1
          lmat(jj)=mat(jorb)
      end do
  end do
  if(jj/=mad%nvctr) then
      write(*,'(a,2(2x,i0))') 'ERROR in compressMatrix: jj/=mad%nvctr',jj,mad%nvctr
      stop
  end if
  
end subroutine compressMatrix



subroutine compressMatrix2(iproc, nproc, orbs, mad, mat, lmat, sendcounts, displs)
  use module_base
  use module_types
  implicit none
  
  ! Calling arguments
  integer,intent(in):: iproc, nproc
  type(orbitals_data),intent(in):: orbs
  type(matrixDescriptors),intent(in):: mad
  real(8),dimension(orbs%norb**2),intent(in):: mat
  real(8),dimension(mad%nvctr),intent(out):: lmat
  integer,dimension(0:nproc-1),intent(out):: sendcounts, displs
  
  ! Local variables
  integer:: iseg, jj, jorb, iiorb, jjorb, jjproc, jjprocold, ncount
  
  sendcounts=0
  displs=0
  
  jj=0
  ncount=0
  jjprocold=0
  displs(0)=0
  do iseg=1,mad%nseg
      do jorb=mad%keyg(1,iseg),mad%keyg(2,iseg)
          jj=jj+1
          lmat(jj)=mat(jorb)
          
          ncount=ncount+1
          jjorb=(jorb-1)/orbs%norb+1
          jjproc=orbs%onWhichMPI(jjorb)
          if(jjproc>jjprocold) then
              ! This part of the matrix is calculated by a new MPI process.
              sendcounts(jjproc-1)=ncount-1
              displs(jjproc)=displs(jjproc-1)+sendcounts(jjproc-1)
              ncount=1
              jjprocold=jjproc
          end if
      end do
  end do
  !sendcounts(nproc-1)=ncount
  sendcounts(jjproc)=ncount !last process
  if(jj/=mad%nvctr) then
      write(*,'(a,2(2x,i0))') 'ERROR in compressMatrix: jj/=mad%nvctr',jj,mad%nvctr
      stop
  end if

  if(sum(sendcounts)/=mad%nvctr) then
      write(*,'(a,2(2x,i0))') 'ERROR in compressMatrix2: sum(sendcounts)/=mad%nvctr',sum(sendcounts),mad%nvctr
      stop
  end if

  !if(iproc==0) then
  !    do jjproc=0,nproc-1
  !        write(*,'(a,3i8)') 'jjproc, displs(jjproc), sendcounts(jjproc)', jjproc, displs(jjproc), sendcounts(jjproc)
  !    end do
  !end if
  
end subroutine compressMatrix2



subroutine compressMatrixPerProcess(iproc, nproc, orbs, mad, mat, size_lmat, lmat)
  use module_base
  use module_types
  implicit none
  
  ! Calling arguments
  integer,intent(in):: iproc, nproc, size_lmat
  type(orbitals_data),intent(in):: orbs
  type(matrixDescriptors),intent(in):: mad
  real(8),dimension(orbs%norb**2),intent(in):: mat
  real(8),dimension(size_lmat),intent(out):: lmat
  
  ! Local variables
  integer:: iseg, jj, jorb, jjorb, jjproc
  
  
  jj=0
  do iseg=1,mad%nseg
      do jorb=mad%keyg(1,iseg),mad%keyg(2,iseg)
          jjorb=(jorb-1)/orbs%norb+1
          jjproc=orbs%onWhichMPI(jjorb)
          if(iproc==jjproc) then
              jj=jj+1
              lmat(jj)=mat(jorb)
          end if
      end do
  end do
  
end subroutine compressMatrixPerProcess




subroutine uncompressMatrix(norb, mad, lmat, mat)
  use module_base
  use module_types
  implicit none
  
  ! Calling arguments
  integer,intent(in):: norb
  type(matrixDescriptors),intent(in):: mad
  real(8),dimension(mad%nvctr),intent(in):: lmat
  real(8),dimension(norb**2),intent(out):: mat
  
  ! Local variables
  integer:: iseg, jj, jorb, iiorb, jjorb
  
  mat=0.d0
  
  jj=0
  do iseg=1,mad%nseg
      do jorb=mad%keyg(1,iseg),mad%keyg(2,iseg)
          jj=jj+1
          mat(jorb)=lmat(jj)
      end do
  end do
  if(jj/=mad%nvctr) then
      write(*,'(a,2(2x,i0))') 'ERROR in uncompressMatrix: jj/=mad%nvctr',jj,mad%nvctr
      stop
  end if
  
end subroutine uncompressMatrix





subroutine dgemm_compressed2(iproc, nproc, norb, nsegline, nseglinemax, keygline, nsegmatmul, keygmatmul, a, b, c)
!! ATTENTION: A MUST BE SYMMETRIC
use module_base
use module_types
implicit none

! Calling arguments
integer,intent(in):: iproc, nproc, norb, nseglinemax, nsegmatmul
integer,dimension(2,nsegmatmul),intent(in):: keygmatmul
integer,dimension(norb):: nsegline
!integer,dimension(2,maxval(nsegline),norb):: keygline
integer,dimension(2,nseglinemax,norb):: keygline
real(8),dimension(norb,norb),intent(in):: a, b
real(8),dimension(norb,norb),intent(out):: c

! Local variables
integer:: iseg, i, irow, icolumn, k, iorb, jorb, korb, jseg, j, jrow, jcolumn, ii
integer:: ierr, istart, iend, iiseg, jjseg, ncount
real(8):: tt, ddot
logical:: iistop, jjstop




c=0.d0
ii=0
do iseg=1,nsegmatmul
    do i=keygmatmul(1,iseg),keygmatmul(2,iseg)
        ii=ii+1
        ! Get the row and column index
        irow=(i-1)/norb+1
        icolumn=i-(irow-1)*norb
        !c(irow,icolumn)=ddot(norb, a(1,irow), 1, b(1,icolumn), 1)
        iiseg=1
        jjseg=1
        iistop=.false.
        jjstop=.false.
        do
            istart=max(keygline(1,iiseg,irow),keygline(1,jjseg,icolumn))
            iend=min(keygline(2,iiseg,irow),keygline(2,jjseg,icolumn))
            ncount=iend-istart+1

            if(ncount>0) then
                tt=ddot(ncount, a(istart,irow), 1, b(istart,icolumn), 1)
            else
                tt=0.d0
            end if
            c(irow,icolumn) = c(irow,icolumn) + tt
            if(iiseg==nsegline(irow)) iistop=.true.
            if(jjseg==nsegline(icolumn)) jjstop=.true.
            if(iistop .and. jjstop) exit
            if((keygline(1,iiseg,irow)<=keygline(1,jjseg,icolumn) .or. jjstop) .and. .not.iistop) then
                iiseg=iiseg+1
            else
                jjseg=jjseg+1
            end if
        end do
        !if(iproc==0) write(*,'(3(a,i0),a,es15.6)') 'process ',iproc,': c(',irow,',',icolumn,')=',c(irow,icolumn)
    end do
end do
!write(*,*) 'ii, norb**2', ii, norb**2
!!do icolumn=1,norb
!!    do irow=1,norb
!!        if(iproc==0) write(200,*) icolumn, irow, c(irow,icolumn)
!!    end do
!!end do



end subroutine dgemm_compressed2



subroutine dgemm_compressed_parallel(iproc, nproc, norb, nsegline, nseglinemax, keygline, &
           nsegmatmul, keygmatmul, norb_par, isorb_par, norbp, a, b, c)
!! ATTENTION: A MUST BE SYMMETRIC
use module_base
use module_types
implicit none

! Calling arguments
integer,intent(in):: iproc, nproc, norb, norbp, nseglinemax, nsegmatmul
integer,dimension(2,nsegmatmul),intent(in):: keygmatmul
integer,dimension(norb):: nsegline
integer,dimension(2,nseglinemax,norb):: keygline
integer,dimension(0:nproc-1),intent(in):: norb_par, isorb_par
real(8),dimension(norb,norb),intent(in):: a, b
real(8),dimension(norb,norb),intent(out):: c

! Local variables
integer:: iseg, i, irow, icolumn, k, iorb, jorb, korb, jseg, j, jrow, jcolumn, ii
integer:: ierr, istart, iend, iiseg, jjseg, ncount, jproc, istat, iall, iirow, iicolumn
real(8):: tt, ddot
logical:: iistop, jjstop
integer,dimension(:),allocatable:: sendcounts, displs
real(8),dimension(:,:),allocatable:: c_loc
character(len=*),parameter:: subname='dgemm_compressed_parallel'


allocate(c_loc(norb,norbp), stat=istat)
call memocc(istat, c_loc, 'c_loc', subname)

!c=0.d0
c_loc=0.d0
ii=0
do iseg=1,nsegmatmul
    do i=keygmatmul(1,iseg),keygmatmul(2,iseg)
        ii=ii+1
        ! Get the row and column index
        !irow=(i-1)/norb+1
        !icolumn=i-(irow-1)*norb
        icolumn=(i-1)/norb+1
        irow=i-(icolumn-1)*norb
        !if(irow>isorb_par(iproc) .and. irow<=isorb_par(min(iproc+1,nproc-1))) then
        if((icolumn>isorb_par(iproc) .and. icolumn<=isorb_par(min(iproc+1,nproc-1)))&
             .or. (iproc==nproc-1 .and. icolumn>isorb_par(iproc))) then
            !iirow=irow-isorb_par(iproc)
            iicolumn=icolumn-isorb_par(iproc)
            ! This process handles this entry of the matrix
            !c(irow,icolumn)=ddot(norb, a(1,irow), 1, b(1,icolumn), 1)
            iiseg=1
            jjseg=1
            iistop=.false.
            jjstop=.false.
            !write(*,'(a,3(i0,a))') 'process ',iproc,' calculates entry (',irow,',',iicolumn,')'
            do
                istart=max(keygline(1,iiseg,irow),keygline(1,jjseg,icolumn))
                iend=min(keygline(2,iiseg,irow),keygline(2,jjseg,icolumn))
                ncount=iend-istart+1

                if(ncount>0) then
                    tt=ddot(ncount, a(istart,irow), 1, b(istart,icolumn), 1)
                    !tt=ddot(ncount, a(istart,icolumn), 1, b(istart,irow), 1)
                else
                    tt=0.d0
                end if
                !c(irow,icolumn) = c(irow,icolumn) + tt
                !c_loc(icolumn,iirow) = c_loc(icolumn,iirow) + tt
                c_loc(irow,iicolumn) = c_loc(irow,iicolumn) + tt
                if(iiseg==nsegline(irow)) iistop=.true.
                if(jjseg==nsegline(icolumn)) jjstop=.true.
                if(iistop .and. jjstop) exit
                if((keygline(1,iiseg,irow)<=keygline(1,jjseg,icolumn) .or. jjstop) .and. .not.iistop) then
                    iiseg=iiseg+1
                else
                    jjseg=jjseg+1
                end if
            end do
            !write(*,'(5(a,i0),a,es15.6)') 'process ',iproc,': c_loc(',irow,',',iicolumn,')=c(',irow,',',icolumn,')=',c_loc(irow,iicolumn)
        end if
    end do
end do
!write(*,*) 'ii, norb**2', ii, norb**2

! Communicate the matrix.
allocate(sendcounts(0:nproc-1), stat=istat)
call memocc(istat, sendcounts, 'sendcounts', subname)
allocate(displs(0:nproc-1), stat=istat)
call memocc(istat, displs, 'displs', subname)

displs(0)=0
do jproc=0,nproc-1
    sendcounts(jproc)=norb*norb_par(jproc)
    if(jproc>0) displs(jproc)=displs(jproc-1)+sendcounts(jproc-1)
end do
if (nproc > 1) then
   call mpi_allgatherv(c_loc(1,1), sendcounts(iproc), mpi_double_precision, c(1,1), sendcounts, displs, &
        mpi_double_precision, mpi_comm_world, ierr)
else
   call vcopy(sendcounts(iproc),c_loc(1,1),1,c(1,1),1)
end if

iall=-product(shape(sendcounts))*kind(sendcounts)
deallocate(sendcounts, stat=istat)
call memocc(istat, iall, 'sendcounts', subname)
iall=-product(shape(displs))*kind(displs)
deallocate(displs, stat=istat)
call memocc(istat, iall, 'displs', subname)
iall=-product(shape(c_loc))*kind(c_loc)
deallocate(c_loc, stat=istat)
call memocc(istat, iall, 'c_loc', subname)

!!do icolumn=1,norb
!!    do irow=1,norb
!!        if(iproc==0) write(201,*) icolumn, irow, c(irow,icolumn)
!!    end do
!!end do


end subroutine dgemm_compressed_parallel






subroutine plotGrid(iproc, nproc, norb, nspinor, nspin, orbitalNumber, llr, glr, atoms, rxyz, hx, hy, hz)
  use module_base
  use module_types
  implicit none
  
  ! Calling arguments
  integer,intent(in):: iproc, nproc, norb, nspinor, nspin, orbitalNumber
  type(locreg_descriptors),intent(in):: llr, glr
  type(atoms_data),intent(in)::atoms
  real(8),dimension(3,atoms%nat),intent(in):: rxyz
  real(8),intent(in):: hx, hy, hz
  
  ! Local variables
  integer:: iseg, jj, j0, j1, ii, i3, i2, i0, i1, i, ishift, iat, ldim, gdim, jjj, istat
  character(len=10):: num
  character(len=20):: filename
  real(8),dimension(:),allocatable:: lphi, phi


    ldim=llr%wfd%nvctr_c+7*llr%wfd%nvctr_f
    gdim=glr%wfd%nvctr_c+7*glr%wfd%nvctr_f
    allocate(lphi(ldim), stat=istat)
    allocate(phi(gdim), stat=istat)
    lphi=1.d0
    phi=0.d0
    call Lpsi_to_global2(iproc, nproc, ldim, gdim, norb, nspinor, nspin, glr, llr, lphi, phi)
  
    write(num,'(i0)') orbitalNumber
    filename='orbital_'//trim(num)
  
    open(unit=2000+iproc,file=trim(filename)//'.xyz',status='unknown')
    !write(2000+iproc,*) llr%wfd%nvctr_c+llr%wfd%nvctr_f+atoms%nat,' atomic'
    write(2000+iproc,*) glr%wfd%nvctr_c+glr%wfd%nvctr_f+llr%wfd%nvctr_c+llr%wfd%nvctr_f+atoms%nat,' atomic'
    if (atoms%geocode=='F') then
       write(2000+iproc,*)'complete simulation grid with low and high resolution points'
    else if (atoms%geocode =='S') then
       write(2000+iproc,'(a,2x,3(1x,1pe24.17))')'surface',atoms%alat1,atoms%alat2,atoms%alat3
    else if (atoms%geocode =='P') then
       write(2000+iproc,'(a,2x,3(1x,1pe24.17))')'periodic',atoms%alat1,atoms%alat2,atoms%alat3
    end if

   do iat=1,atoms%nat
      write(2000+iproc,'(a6,2x,3(1x,e12.5),3x)') trim(atoms%atomnames(atoms%iatype(iat))),rxyz(1,iat),rxyz(2,iat),rxyz(3,iat)
   end do

  
    jjj=0
    do iseg=1,glr%wfd%nseg_c
       jj=glr%wfd%keyvloc(iseg)
       j0=glr%wfd%keygloc(1,iseg)
       j1=glr%wfd%keygloc(2,iseg)
       ii=j0-1
       i3=ii/((glr%d%n1+1)*(glr%d%n2+1))
       ii=ii-i3*(glr%d%n1+1)*(glr%d%n2+1)
       i2=ii/(glr%d%n1+1)
       i0=ii-i2*(glr%d%n1+1)
       i1=i0+j1-j0
       do i=i0,i1
           jjj=jjj+1
           if(phi(jjj)==1.d0) write(2000+iproc,'(a4,2x,3(1x,e10.3))') '  lg ',&
                real(i,kind=8)*hx,real(i2,kind=8)*hy,real(i3,kind=8)*hz
           write(2000+iproc,'(a4,2x,3(1x,e10.3))') '  g ',real(i,kind=8)*hx,&
                real(i2,kind=8)*hy,real(i3,kind=8)*hz
       enddo
    enddo

    ishift=glr%wfd%nseg_c  
    ! fine part
    do iseg=1,glr%wfd%nseg_f
       jj=glr%wfd%keyvloc(ishift+iseg)
       j0=glr%wfd%keygloc(1,ishift+iseg)
       j1=glr%wfd%keygloc(2,ishift+iseg)
       ii=j0-1
       i3=ii/((glr%d%n1+1)*(glr%d%n2+1))
       ii=ii-i3*(glr%d%n1+1)*(glr%d%n2+1)
       i2=ii/(glr%d%n1+1)
       i0=ii-i2*(glr%d%n1+1)
       i1=i0+j1-j0
       do i=i0,i1
          jjj=jjj+1
          if(phi(jjj)==1.d0) write(2000+iproc,'(a4,2x,3(1x,e10.3))') '  lG ',real(i,kind=8)*hx,real(i2,kind=8)*hy,real(i3,kind=8)*hz
          write(2000+iproc,'(a4,2x,3(1x,e10.3))') '  G ',real(i,kind=8)*hx,real(i2,kind=8)*hy,real(i3,kind=8)*hz
          jjj=jjj+6
       enddo
    enddo
  
    close(unit=2000+iproc)

end subroutine plotGrid



subroutine local_potential_dimensions(Lzd,orbs,ndimfirstproc)
  use module_base
  use module_types
  use module_xc
  implicit none
  integer, intent(in) :: ndimfirstproc
  type(local_zone_descriptors), intent(inout) :: Lzd
  type(orbitals_data), intent(inout) :: orbs
  !local variables
  character(len=*), parameter :: subname='local_potential_dimensions'
  logical :: newvalue
  integer :: i_all,i_stat,ii,iilr,ilr,iorb,iorb2,nilr,ispin
  integer,dimension(:,:),allocatable:: ilrtable
  
  if(Lzd%nlr > 1) then
     allocate(ilrtable(orbs%norbp,2),stat=i_stat)
     call memocc(i_stat,ilrtable,'ilrtable',subname)
     !call to_zero(orbs%norbp*2,ilrtable(1,1))
     ilrtable=0
     ii=0
     do iorb=1,orbs%norbp
        newvalue=.true.
        !localization region to which the orbital belongs
        ilr = orbs%inwhichlocreg(iorb+orbs%isorb)
        !spin state of the orbital
        if (orbs%spinsgn(orbs%isorb+iorb) > 0.0_gp) then
           ispin = 1       
        else
           ispin=2
        end if
        !check if the orbitals already visited have the same conditions
        loop_iorb2: do iorb2=1,orbs%norbp
           if(ilrtable(iorb2,1) == ilr .and. ilrtable(iorb2,2)==ispin) then
              newvalue=.false.
              exit loop_iorb2
           end if
        end do loop_iorb2
        if (newvalue) then
           ii = ii + 1
           ilrtable(ii,1)=ilr
           ilrtable(ii,2)=ispin    !SOMETHING IS NOT WORKING IN THE CONCEPT HERE... ispin is not a property of the locregs, but of the orbitals
        end if
     end do
     !number of inequivalent potential regions
     nilr = ii

     !calculate the dimension of the potential in the gathered form
     lzd%ndimpotisf=0
     do iilr=1,nilr
        ilr=ilrtable(iilr,1)
        do iorb=1,orbs%norbp
           !put the starting point
           if (orbs%inWhichLocreg(iorb+orbs%isorb) == ilr) then
              !assignment of ispot array to the value of the starting address of inequivalent
              orbs%ispot(iorb)=lzd%ndimpotisf + 1
              if(orbs%spinsgn(orbs%isorb+iorb) <= 0.0_gp) then
                 orbs%ispot(iorb)=lzd%ndimpotisf + &
                      1 + lzd%llr(ilr)%d%n1i*lzd%llr(ilr)%d%n2i*lzd%llr(ilr)%d%n3i
              end if
           end if
        end do
        lzd%ndimpotisf = lzd%ndimpotisf + &
             lzd%llr(ilr)%d%n1i*lzd%llr(ilr)%d%n2i*lzd%llr(ilr)%d%n3i*orbs%nspin
     end do
     !part which refers to exact exchange (only meaningful for one region)
     if (xc_exctXfac() /= 0.0_gp) then
        lzd%ndimpotisf = lzd%ndimpotisf + &
             max(max(lzd%llr(ilr)%d%n1i*lzd%llr(ilr)%d%n2i*lzd%llr(ilr)%d%n3i*orbs%norbp,ndimfirstproc*orbs%norb),1)
     end if

  else 
     allocate(ilrtable(1,2),stat=i_stat)
     call memocc(i_stat,ilrtable,'ilrtable',subname)
     nilr = 1
     ilrtable=1

     !calculate the dimension of the potential in the gathered form
     lzd%ndimpotisf=0
     do iorb=1,orbs%norbp
        !assignment of ispot array to the value of the starting address of inequivalent
        orbs%ispot(iorb)=lzd%ndimpotisf + 1
        if(orbs%spinsgn(orbs%isorb+iorb) <= 0.0_gp) then
           orbs%ispot(iorb)=lzd%ndimpotisf + &
                1 + lzd%Glr%d%n1i*lzd%Glr%d%n2i*lzd%Glr%d%n3i
        end if
     end do
     lzd%ndimpotisf = lzd%ndimpotisf + &
          lzd%Glr%d%n1i*lzd%Glr%d%n2i*lzd%Glr%d%n3i*orbs%nspin
          
     !part which refers to exact exchange (only meaningful for one region)
     if (xc_exctXfac() /= 0.0_gp) then
        lzd%ndimpotisf = lzd%ndimpotisf + &
             max(max(lzd%Glr%d%n1i*lzd%Glr%d%n2i*lzd%Glr%d%n3i*orbs%norbp,ndimfirstproc*orbs%norb),1)
     end if


  end if


  i_all=-product(shape(ilrtable))*kind(ilrtable)
  deallocate(ilrtable,stat=i_stat)
  call memocc(i_stat,i_all,'ilrtable',subname)

end subroutine local_potential_dimensions



subroutine print_orbital_distribution(iproc, nproc, orbs, derorbs)
use module_base
use module_types
implicit none

integer,intent(in):: iproc, nproc
type(orbitals_data),intent(in):: orbs, derorbs

! Local variables
integer:: jproc, len1, len2, space1, space2
logical:: written

write(*,'(1x,a)') '------------------------------------------------------------------------------------'
written=.false.
write(*,'(1x,a)') '>>>> Partition of the basis functions among the processes.'
do jproc=1,nproc-1
    if(orbs%norb_par(jproc,0)<orbs%norb_par(jproc-1,0)) then
        len1=1+ceiling(log10(dble(jproc-1)+1.d-5))+ceiling(log10(dble(orbs%norb_par(jproc-1,0)+1.d-5)))
        len2=ceiling(log10(dble(jproc)+1.d-5))+ceiling(log10(dble(nproc-1)+1.d-5))+&
             ceiling(log10(dble(orbs%norb_par(jproc,0)+1.d-5)))
        if(len1>=len2) then
            space1=1
            space2=1+len1-len2
        else
            space1=1+len2-len1
            space2=1
        end if
        write(*,'(4x,a,2(i0,a),a,a)') '| Processes from 0 to ',jproc-1,' treat ',&
            orbs%norb_par(jproc-1,0), ' orbitals,', repeat(' ', space1), '|'
        write(*,'(4x,a,3(i0,a),a,a)')  '| processes from ',jproc,' to ',nproc-1,' treat ', &
            orbs%norb_par(jproc,0),' orbitals.', repeat(' ', space2), '|'
        written=.true.
        exit
    end if
end do
if(.not.written) then
    write(*,'(4x,a,2(i0,a),a,a)') '| Processes from 0 to ',nproc-1, &
        ' treat ',orbs%norbp,' orbitals. |'!, &
end if
write(*,'(1x,a)') '-----------------------------------------------'

written=.false.
write(*,'(1x,a)') '>>>> Partition of the basis functions including the derivatives among the processes.'
do jproc=1,nproc-1
    if(derorbs%norb_par(jproc,0)<derorbs%norb_par(jproc-1,0)) then
        len1=1+ceiling(log10(dble(jproc-1)+1.d-5))+ceiling(log10(dble(derorbs%norb_par(jproc-1,0)+1.d-5)))
        len2=ceiling(log10(dble(jproc)+1.d-5))+ceiling(log10(dble(nproc-1)+1.d-5))+&
             ceiling(log10(dble(derorbs%norb_par(jproc,0)+1.d-5)))
        if(len1>=len2) then
            space1=1
            space2=1+len1-len2
        else
            space1=1+len2-len1
            space2=1
        end if
        write(*,'(4x,a,2(i0,a),a,a)') '| Processes from 0 to ',jproc-1,' treat ',&
            derorbs%norb_par(jproc-1,0), ' orbitals,', repeat(' ', space1), '|'
        write(*,'(4x,a,3(i0,a),a,a)')  '| processes from ',jproc,' to ',nproc-1,' treat ', &
            derorbs%norb_par(jproc,0),' orbitals.', repeat(' ', space2), '|'
        written=.true.
        exit
    end if
end do
if(.not.written) then
    write(*,'(4x,a,2(i0,a),a,a)') '| Processes from 0 to ',nproc-1, &
        ' treat ',derorbs%norbp,' orbitals. |'
end if
write(*,'(1x,a)') '------------------------------------------------------------------------------------'


end subroutine print_orbital_distribution
