!> @file 
!!   Miscellaneous routines for linear toolbox
!! @author
!!   Copyright (C) 2011-2012 BigDFT group 
!!   This file is distributed under the terms of the
!!   GNU General Public License, see ~/COPYING file
!!   or http://www.gnu.org/copyleft/gpl.txt .
!!   For the list of contributors, see ~/AUTHORS 
 

!> Plots the orbitals
subroutine plotOrbitals(iproc, orbs, Glr, phi, nat, rxyz, hxh, hyh, hzh, it)
use module_base
use module_types
implicit none

! Calling arguments
integer :: iproc
type(orbitals_data), intent(inout) :: orbs
type(locreg_descriptors), intent(in) :: Glr
real(kind=8), dimension((Glr%wfd%nvctr_c+7*Glr%wfd%nvctr_f)*orbs%nspinor*orbs%norbp) :: phi
integer :: nat
real(kind=8), dimension(3,nat) :: rxyz
real(kind=8) :: hxh, hyh, hzh
integer :: it

integer :: ix, iy, iz, ix0, iy0, iz0, iiAt, jj, iorb, i1, i2, i3, istart, ii, istat
integer :: unit1, unit2, unit3
real(kind=8), dimension(:), allocatable :: phir
type(workarr_sumrho) :: w
character(len=10) :: c1, c2, c3
character(len=50) :: file1, file2, file3

allocate(phir(Glr%d%n1i*Glr%d%n2i*Glr%d%n3i), stat=istat)

call initialize_work_arrays_sumrho(Glr,w)

istart=0

unit1=10*iproc+7
unit2=10*iproc+8
unit3=10*iproc+9

!write(*,*) 'write, orbs%nbasisp', orbs%norbp
    orbLoop: do iorb=1,orbs%norbp
        !!phir=0.d0
        call to_zero(Glr%d%n1i*Glr%d%n2i*Glr%d%n3i, phir(1))
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




subroutine compressMatrix2(iproc, nproc, orbs, mad, mat, lmat, sendcounts, displs)
  use module_base
  use module_types
  implicit none
  
  ! Calling arguments
  integer, intent(in) :: iproc, nproc
  type(orbitals_data), intent(in) :: orbs
  type(matrixDescriptors), intent(in) :: mad
  real(kind=8), dimension(orbs%norb**2), intent(in) :: mat
  real(kind=8), dimension(mad%nvctr), intent(out) :: lmat
  integer, dimension(0:nproc-1), intent(out) :: sendcounts, displs
  
  ! Local variables
  integer :: iseg, jj, jorb, jjorb, jjproc, jjprocold, ncount
  
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




subroutine uncompressMatrix(norb, mad, lmat, mat)
  use module_base
  use module_types
  implicit none
  
  ! Calling arguments
  integer, intent(in) :: norb
  type(matrixDescriptors), intent(in) :: mad
  real(kind=8), dimension(mad%nvctr), intent(in) :: lmat
  real(kind=8), dimension(norb**2), intent(out) :: mat
  
  ! Local variables
  integer :: iseg, ii, jorb

  
  call to_zero(norb**2, mat(1))
  
  !$omp parallel do default(private) shared(mad,lmat,mat)

  do iseg=1,mad%nseg
      ii=0
      do jorb=mad%keyg(1,iseg),mad%keyg(2,iseg)
          mat(jorb)=lmat(mad%keyv(iseg)+ii)
          ii=ii+1
      end do
  end do

  !$omp end parallel do
  
end subroutine uncompressMatrix




subroutine plotGrid(iproc, norb, nspinor, nspin, orbitalNumber, llr, glr, atoms, rxyz, hx, hy, hz)
  use module_base
  use module_types
  implicit none
  
  ! Calling arguments
  integer, intent(in) :: iproc, norb, nspinor, nspin, orbitalNumber
  type(locreg_descriptors), intent(in) :: llr, glr
  type(atoms_data), intent(in) ::atoms
  real(kind=8), dimension(3,atoms%nat), intent(in) :: rxyz
  real(kind=8), intent(in) :: hx, hy, hz
  
  ! Local variables
  integer :: iseg, jj, j0, j1, ii, i3, i2, i0, i1, i, ishift, iat, ldim, gdim, jjj, istat
  character(len=10) :: num
  character(len=20) :: filename
  real(kind=8), dimension(:), allocatable :: lphi, phi


    ldim=llr%wfd%nvctr_c+7*llr%wfd%nvctr_f
    gdim=glr%wfd%nvctr_c+7*glr%wfd%nvctr_f
    allocate(lphi(ldim), stat=istat)
    allocate(phi(gdim), stat=istat)
    lphi=1.d0
    !!phi=0.d0
    call to_zero(gdim, phi(1))
    call Lpsi_to_global2(iproc, ldim, gdim, norb, nspinor, nspin, glr, llr, lphi, phi)
  
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
  integer, dimension(:,:), allocatable :: ilrtable
  
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



subroutine print_orbital_distribution(iproc, nproc, orbs)
use module_base
use module_types
implicit none

integer, intent(in) :: iproc, nproc
type(orbitals_data), intent(in) :: orbs

! Local variables
integer :: jproc, len1, len2, space1, space2
logical :: written

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

!!written=.false.
!!write(*,'(1x,a)') '>>>> Partition of the basis functions including the derivatives among the processes.'
!!do jproc=1,nproc-1
!!    if(derorbs%norb_par(jproc,0)<derorbs%norb_par(jproc-1,0)) then
!!        len1=1+ceiling(log10(dble(jproc-1)+1.d-5))+ceiling(log10(dble(derorbs%norb_par(jproc-1,0)+1.d-5)))
!!        len2=ceiling(log10(dble(jproc)+1.d-5))+ceiling(log10(dble(nproc-1)+1.d-5))+&
!!             ceiling(log10(dble(derorbs%norb_par(jproc,0)+1.d-5)))
!!        if(len1>=len2) then
!!            space1=1
!!            space2=1+len1-len2
!!        else
!!            space1=1+len2-len1
!!            space2=1
!!        end if
!!        write(*,'(4x,a,2(i0,a),a,a)') '| Processes from 0 to ',jproc-1,' treat ',&
!!            derorbs%norb_par(jproc-1,0), ' orbitals,', repeat(' ', space1), '|'
!!        write(*,'(4x,a,3(i0,a),a,a)')  '| processes from ',jproc,' to ',nproc-1,' treat ', &
!!            derorbs%norb_par(jproc,0),' orbitals.', repeat(' ', space2), '|'
!!        written=.true.
!!        exit
!!    end if
!!end do
!!if(.not.written) then
!!    write(*,'(4x,a,2(i0,a),a,a)') '| Processes from 0 to ',nproc-1, &
!!        ' treat ',derorbs%norbp,' orbitals. |'
!!end if
!!write(*,'(1x,a)') '------------------------------------------------------------------------------------'


end subroutine print_orbital_distribution

