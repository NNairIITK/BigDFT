!> @file
!!  Routines to do restart
!! @author
!!    Copyright (C) 2007-2011 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 


!>  Copy old wavefunctions from psi to psi_old
subroutine copy_old_wavefunctions(nproc,orbs,n1,n2,n3,wfd,psi,&
     n1_old,n2_old,n3_old,wfd_old,psi_old)
  use module_base
  use module_types
  implicit none
  integer, intent(in) :: nproc,n1,n2,n3
  type(orbitals_data), intent(in) :: orbs
  type(wavefunctions_descriptors), intent(inout) :: wfd,wfd_old
  integer, intent(out) :: n1_old,n2_old,n3_old
  real(wp), dimension(:), pointer :: psi,psi_old
  !Local variables
  character(len=*), parameter :: subname='copy_old_wavefunctions'
  real(kind=8), parameter :: eps_mach=1.d-12
  integer :: iseg,j,ind1,iorb,i_all,i_stat,oidx,sidx !n(c) nvctrp_old
  real(kind=8) :: tt

  wfd_old%nvctr_c = wfd%nvctr_c
  wfd_old%nvctr_f = wfd%nvctr_f
  wfd_old%nseg_c  = wfd%nseg_c
  wfd_old%nseg_f  = wfd%nseg_f

  !allocations
  call allocate_wfd(wfd_old,subname)

  do iseg=1,wfd_old%nseg_c+wfd_old%nseg_f
     wfd_old%keyglob(1,iseg)    = wfd%keyglob(1,iseg) 
     wfd_old%keyglob(2,iseg)    = wfd%keyglob(2,iseg)
     wfd_old%keygloc(1,iseg)    = wfd%keygloc(1,iseg)
     wfd_old%keygloc(2,iseg)    = wfd%keygloc(2,iseg)
     wfd_old%keyvloc(iseg)      = wfd%keyvloc(iseg)
     wfd_old%keyvglob(iseg)      = wfd%keyvglob(iseg)
  enddo
  !deallocation
  call deallocate_wfd(wfd,subname)

  n1_old = n1
  n2_old = n2
  n3_old = n3

  !add the number of distributed point for the compressed wavefunction
  tt=dble(wfd_old%nvctr_c+7*wfd_old%nvctr_f)/dble(nproc)
  !n(c) nvctrp_old=int((1.d0-eps_mach*tt) + tt)

  allocate(psi_old((wfd_old%nvctr_c+7*wfd_old%nvctr_f)*orbs%norbp*orbs%nspinor+ndebug),&
       stat=i_stat)
  call memocc(i_stat,psi_old,'psi_old',subname)

  do iorb=1,orbs%norbp
     tt=0.d0
     oidx=(iorb-1)*orbs%nspinor+1
     do sidx=oidx,oidx+orbs%nspinor-1
        do j=1,wfd_old%nvctr_c+7*wfd_old%nvctr_f
           ind1=j+(wfd_old%nvctr_c+7*wfd_old%nvctr_f)*(sidx-1)
           psi_old(ind1)= psi(ind1)
           tt=tt+real(psi(ind1),kind=8)**2
        enddo
     end do

     tt=sqrt(tt)
     if (abs(tt-1.d0) > 1.d-8) then
        write(*,*)'wrong psi_old',iorb,tt
        stop 
     end if
  enddo
  !deallocation
  i_all=-product(shape(psi))*kind(psi)
  deallocate(psi,stat=i_stat)
  call memocc(i_stat,i_all,'psi',subname)

END SUBROUTINE copy_old_wavefunctions


!>   Reformat wavefunctions if the mesh have changed (in a restart)
subroutine reformatmywaves(iproc,orbs,at,&
     hx_old,hy_old,hz_old,n1_old,n2_old,n3_old,rxyz_old,wfd_old,psi_old,&
     hx,hy,hz,n1,n2,n3,rxyz,wfd,psi)
  use module_base
  use module_types
  implicit none
  integer, intent(in) :: iproc,n1_old,n2_old,n3_old,n1,n2,n3
  real(gp), intent(in) :: hx_old,hy_old,hz_old,hx,hy,hz
  type(wavefunctions_descriptors), intent(in) :: wfd,wfd_old
  type(atoms_data), intent(in) :: at
  type(orbitals_data), intent(in) :: orbs
  real(gp), dimension(3,at%nat), intent(in) :: rxyz,rxyz_old
  real(wp), dimension(wfd_old%nvctr_c+7*wfd_old%nvctr_f,orbs%nspinor*orbs%norbp), intent(in) :: psi_old
  real(wp), dimension(wfd%nvctr_c+7*wfd%nvctr_f,orbs%nspinor*orbs%norbp), intent(out) :: psi
  !Local variables
  character(len=*), parameter :: subname='reformatmywaves'
  logical :: reformat,perx,pery,perz
  integer :: iat,iorb,j,i_stat,i_all,jj,j0,j1,ii,i0,i1,i2,i3,i,iseg,nb1,nb2,nb3
  real(gp) :: tx,ty,tz,displ,mindist
  real(wp), dimension(:,:,:), allocatable :: psifscf
  real(wp), dimension(:,:,:,:,:,:), allocatable :: psigold

  !conditions for periodicity in the three directions
  perx=(at%geocode /= 'F')
  pery=(at%geocode == 'P')
  perz=(at%geocode /= 'F')

  !buffers realted to periodicity
  !WARNING: the boundary conditions are not assumed to change between new and old
  call ext_buffers_coarse(perx,nb1)
  call ext_buffers_coarse(pery,nb2)
  call ext_buffers_coarse(perz,nb3)


  allocate(psifscf(-nb1:2*n1+1+nb1,-nb2:2*n2+1+nb2,-nb3:2*n3+1+nb3+ndebug),stat=i_stat)
  call memocc(i_stat,psifscf,'psifscf',subname)

  tx=0.0_gp 
  ty=0.0_gp
  tz=0.0_gp

  do iat=1,at%nat
     tx=tx+mindist(perx,at%alat1,rxyz(1,iat),rxyz_old(1,iat))**2
     ty=ty+mindist(pery,at%alat2,rxyz(2,iat),rxyz_old(2,iat))**2
     tz=tz+mindist(perz,at%alat3,rxyz(3,iat),rxyz_old(3,iat))**2
  enddo
  displ=sqrt(tx+ty+tz)
!  write(100+iproc,*) 'displacement',dis
!  write(100+iproc,*) 'rxyz ',rxyz
!  write(100+iproc,*) 'rxyz_old ',rxyz_old

  !reformatting criterion
  if (hx == hx_old .and. hy == hy_old .and. hz == hz_old .and. &
       wfd_old%nvctr_c  == wfd%nvctr_c .and. wfd_old%nvctr_f == wfd%nvctr_f .and.&
       n1_old  == n1  .and. n2_old == n2 .and. n3_old == n3  .and.  displ <  1.d-3  ) then
     reformat=.false.
     if (iproc==0) then
        write(*,'(1x,a)',advance='NO')&
         'The wavefunctions do not need reformatting and can be imported directly...   '
       !  '-------------------------------------------------------------- Wavefunctions Restart'
     end if
  else
     reformat=.true.
     if (iproc==0) then
        write(*,'(1x,a)')&
         'The wavefunctions need reformatting because:                                 '
        if (hx /= hx_old .or. hy /= hy_old .or. hz /= hz_old) then 
           write(*,"(4x,a,6(1pe20.12))") &
                '  hgrid_old /= hgrid  ',hx_old,hy_old,hz_old,hx,hy,hz
        else if (wfd_old%nvctr_c /= wfd%nvctr_c) then
           write(*,"(4x,a,2i8)") &
                'nvctr_c_old /= nvctr_c',wfd_old%nvctr_c,wfd%nvctr_c
        else if (wfd_old%nvctr_f /= wfd%nvctr_f)  then
           write(*,"(4x,a,2i8)") &
                'nvctr_f_old /= nvctr_f',wfd_old%nvctr_f,wfd%nvctr_f
        else if (n1_old /= n1  .or. n2_old /= n2 .or. n3_old /= n3 )  then  
           write(*,"(4x,a,6i5)") &
                'cell size has changed ',n1_old,n1  , n2_old,n2 , n3_old,n3
        else
           write(*,"(4x,a,3(1pe19.12))") &
                'molecule was shifted  ' , tx,ty,tz
        endif
           write(*,"(1x,a)",advance='NO')& 
                'Reformatting...'
     end if
     !calculate the new grid values
     
!check
!        write(100+iproc,'(1x,a)')&
!         'The wavefunctions need reformatting because:                                 '
!        if (hgrid_old.ne.hgrid) then 
!           write(100+iproc,"(4x,a,1pe20.12)") &
!                '  hgrid_old /= hgrid  ',hgrid_old, hgrid
!        else if (wfd_old%nvctr_c.ne.wfd%nvctr_c) then
!           write(100+iproc,"(4x,a,2i8)") &
!                'nvctr_c_old /= nvctr_c',wfd_old%nvctr_c,wfd%nvctr_c
!        else if (wfd_old%nvctr_f.ne.wfd%nvctr_f)  then
!           write(100+iproc,"(4x,a,2i8)") &
!                'nvctr_f_old /= nvctr_f',wfd_old%nvctr_f,wfd%nvctr_f
!        else if (n1_old.ne.n1  .or. n2_old.ne.n2 .or. n3_old.ne.n3 )  then  
!           write(100+iproc,"(4x,a,6i5)") &
!                'cell size has changed ',n1_old,n1  , n2_old,n2 , n3_old,n3
!        else
!           write(100+iproc,"(4x,a,3(1pe19.12))") &
!                'molecule was shifted  ' , tx,ty,tz
!        endif
!checkend
  end if

  do iorb=1,orbs%norbp*orbs%nspinor

     if (.not. reformat) then
!write(100+iproc,*) 'no reformatting' 

        do j=1,wfd_old%nvctr_c
           psi(j,iorb)=psi_old(j, iorb)
        enddo
        do j=1,7*wfd_old%nvctr_f-6,7
           psi(wfd%nvctr_c+j+0,iorb)=psi_old(wfd%nvctr_c+j+0,iorb)
           psi(wfd%nvctr_c+j+1,iorb)=psi_old(wfd%nvctr_c+j+1,iorb)
           psi(wfd%nvctr_c+j+2,iorb)=psi_old(wfd%nvctr_c+j+2,iorb)
           psi(wfd%nvctr_c+j+3,iorb)=psi_old(wfd%nvctr_c+j+3,iorb)
           psi(wfd%nvctr_c+j+4,iorb)=psi_old(wfd%nvctr_c+j+4,iorb)
           psi(wfd%nvctr_c+j+5,iorb)=psi_old(wfd%nvctr_c+j+5,iorb)
           psi(wfd%nvctr_c+j+6,iorb)=psi_old(wfd%nvctr_c+j+6,iorb)
        enddo

     else

        allocate(psigold(0:n1_old,2,0:n2_old,2,0:n3_old,2+ndebug),stat=i_stat)
        call memocc(i_stat,psigold,'psigold',subname)

        call razero(8*(n1_old+1)*(n2_old+1)*(n3_old+1),psigold)

        ! coarse part
        do iseg=1,wfd_old%nseg_c
           jj=wfd_old%keyvloc(iseg)
           j0=wfd_old%keygloc(1,iseg)
           j1=wfd_old%keygloc(2,iseg)
           ii=j0-1
           i3=ii/((n1_old+1)*(n2_old+1))
           ii=ii-i3*(n1_old+1)*(n2_old+1)
           i2=ii/(n1_old+1)
           i0=ii-i2*(n1_old+1)
           i1=i0+j1-j0
           do i=i0,i1
              psigold(i,1,i2,1,i3,1) = psi_old(i-i0+jj,iorb)
           enddo
        enddo

        ! fine part
        do iseg=1,wfd_old%nseg_f
           jj=wfd_old%keyvloc(wfd_old%nseg_c + iseg)
           j0=wfd_old%keygloc(1,wfd_old%nseg_c + iseg)
           j1=wfd_old%keygloc(2,wfd_old%nseg_c + iseg)
           ii=j0-1
           i3=ii/((n1_old+1)*(n2_old+1))
           ii=ii-i3*(n1_old+1)*(n2_old+1)
           i2=ii/(n1_old+1)
           i0=ii-i2*(n1_old+1)
           i1=i0+j1-j0
           do i=i0,i1
              psigold(i,2,i2,1,i3,1)=psi_old(wfd_old%nvctr_c+1+7*(i-i0+jj-1), iorb)
              psigold(i,1,i2,2,i3,1)=psi_old(wfd_old%nvctr_c+2+7*(i-i0+jj-1), iorb)
              psigold(i,2,i2,2,i3,1)=psi_old(wfd_old%nvctr_c+3+7*(i-i0+jj-1), iorb)
              psigold(i,1,i2,1,i3,2)=psi_old(wfd_old%nvctr_c+4+7*(i-i0+jj-1), iorb)
              psigold(i,2,i2,1,i3,2)=psi_old(wfd_old%nvctr_c+5+7*(i-i0+jj-1), iorb)
              psigold(i,1,i2,2,i3,2)=psi_old(wfd_old%nvctr_c+6+7*(i-i0+jj-1), iorb)
              psigold(i,2,i2,2,i3,2)=psi_old(wfd_old%nvctr_c+7+7*(i-i0+jj-1), iorb)
           enddo
        enddo

!write(100+iproc,*) 'norm psigold ',dnrm2(8*(n1_old+1)*(n2_old+1)*(n3_old+1),psigold,1)

        call reformatonewave(displ,wfd,at,hx_old,hy_old,hz_old, & !n(m)
             n1_old,n2_old,n3_old,rxyz_old,psigold,hx,hy,hz,&
             n1,n2,n3,rxyz,psifscf,psi(1,iorb))

        i_all=-product(shape(psigold))*kind(psigold)
        deallocate(psigold,stat=i_stat)
        call memocc(i_stat,i_all,'psigold',subname)
     end if
  end do

  i_all=-product(shape(psifscf))*kind(psifscf)
  deallocate(psifscf,stat=i_stat)
  call memocc(i_stat,i_all,'psifscf',subname)

  if (iproc==0) write(*,"(1x,a)")'done.'

END SUBROUTINE reformatmywaves

integer function wave_format_from_filename(iproc, filename)
  use module_types
  implicit none
  integer, intent(in) :: iproc
  character(len=*), intent(in) :: filename

  integer :: isuffix

  wave_format_from_filename = WF_FORMAT_NONE

  isuffix = index(filename, ".etsf", back = .true.)
  if (isuffix > 0) then
     wave_format_from_filename = WF_FORMAT_ETSF
     if (iproc ==0) write(*,*) "Reading wavefunctions in ETSF file format."
  else
     isuffix = index(filename, ".bin", back = .true.)
     if (isuffix > 0) then
        wave_format_from_filename = WF_FORMAT_BINARY
        if (iproc ==0) write(*,*) "Reading wavefunctions in BigDFT binary file format."
     else
        wave_format_from_filename = WF_FORMAT_PLAIN
        if (iproc ==0) write(*,*) "Reading wavefunctions in plain text file format."
     end if
  end if
end function wave_format_from_filename

!>  Reads wavefunction from file and transforms it properly if hgrid or size of simulation cell
!!  have changed
subroutine readmywaves(iproc,filename,iformat,orbs,n1,n2,n3,hx,hy,hz,at,rxyz_old,rxyz,  & 
     wfd,psi,orblist)
  use module_base
  use module_types
  use module_interfaces, except_this_one => readmywaves
  implicit none
  integer, intent(in) :: iproc,n1,n2,n3, iformat
  real(gp), intent(in) :: hx,hy,hz
  type(wavefunctions_descriptors), intent(in) :: wfd
  type(orbitals_data), intent(inout) :: orbs
  type(atoms_data), intent(in) :: at
  real(gp), dimension(3,at%nat), intent(in) :: rxyz
  real(gp), dimension(3,at%nat), intent(out) :: rxyz_old
  real(wp), dimension(wfd%nvctr_c+7*wfd%nvctr_f,orbs%nspinor,orbs%norbp), intent(out) :: psi
  character(len=*), intent(in) :: filename
  integer, dimension(orbs%norb), optional :: orblist
  !Local variables
  character(len=*), parameter :: subname='readmywaves'
  logical :: perx,pery,perz
  integer :: ncount1,ncount_rate,ncount_max,iorb,i_stat,i_all,ncount2,nb1,nb2,nb3,iorb_out,ispinor
  real(kind=4) :: tr0,tr1
  real(kind=8) :: tel
  real(wp), dimension(:,:,:), allocatable :: psifscf
  !integer, dimension(orbs%norb) :: orblist2

  call cpu_time(tr0)
  call system_clock(ncount1,ncount_rate,ncount_max)

  if (iformat == WF_FORMAT_ETSF) then
     !construct the orblist or use the one in argument
     !do nb1 = 1, orbs%norb
     !orblist2(nb1) = nb1
     !if(present(orblist)) orblist2(nb1) = orblist(nb1) 
     !end do

     call read_waves_etsf(iproc,filename // ".etsf",orbs,n1,n2,n3,hx,hy,hz,at,rxyz_old,rxyz,  & 
          wfd,psi)
  else if (iformat == WF_FORMAT_BINARY .or. iformat == WF_FORMAT_PLAIN) then
     !conditions for periodicity in the three directions
     perx=(at%geocode /= 'F')
     pery=(at%geocode == 'P')
     perz=(at%geocode /= 'F')

     !buffers related to periodicity
     !WARNING: the boundary conditions are not assumed to change between new and old
     call ext_buffers_coarse(perx,nb1)
     call ext_buffers_coarse(pery,nb2)
     call ext_buffers_coarse(perz,nb3)

     allocate(psifscf(-nb1:2*n1+1+nb1,-nb2:2*n2+1+nb2,-nb3:2*n3+1+nb3+ndebug),stat=i_stat)
     call memocc(i_stat,psifscf,'psifscf',subname)

     do iorb=1,orbs%norbp!*orbs%nspinor

!!$        write(f4,'(i4.4)') iorb+orbs%isorb*orbs%nspinor
!!$        if (exists) then
!!$           filename_ = filename//".bin."//f4
!!$           open(unit=99,file=filename_,status='unknown',form="unformatted")
!!$        else
!!$           filename_ = trim(filename)//"."//f4
!!$           open(unit=99,file=trim(filename_),status='unknown')
!!$        end if
!!$           call readonewave(99, .not.exists,iorb+orbs%isorb*orbs%nspinor,iproc,n1,n2,n3, &
!!$                & hx,hy,hz,at,wfd,rxyz_old,rxyz,&
!!$                psi(1,iorb),orbs%eval((iorb-1)/orbs%nspinor+1+orbs%isorb),psifscf)

        do ispinor=1,orbs%nspinor
           if(present(orblist)) then
              call open_filename_of_iorb(99,(iformat == WF_FORMAT_BINARY),filename, &
                   & orbs,iorb,ispinor,iorb_out, orblist(iorb+orbs%isorb))
           else
              call open_filename_of_iorb(99,(iformat == WF_FORMAT_BINARY),filename, &
                   & orbs,iorb,ispinor,iorb_out)
           end if           
           call readonewave(99, (iformat == WF_FORMAT_PLAIN),iorb_out,iproc,n1,n2,n3, &
                & hx,hy,hz,at,wfd,rxyz_old,rxyz,&
                psi(1,ispinor,iorb),orbs%eval(orbs%isorb+iorb),psifscf)
           close(99)
        end do

     end do

     i_all=-product(shape(psifscf))*kind(psifscf)
     deallocate(psifscf,stat=i_stat)
     call memocc(i_stat,i_all,'psifscf',subname)

  else
     write(0,*) "Unknown wavefunction file format from filename."
     stop
  end if

  call cpu_time(tr1)
  call system_clock(ncount2,ncount_rate,ncount_max)
  tel=dble(ncount2-ncount1)/dble(ncount_rate)
  write(*,'(a,i4,2(1x,1pe10.3))') '- READING WAVES TIME',iproc,tr1-tr0,tel
END SUBROUTINE readmywaves

!> Verify the presence of a given file
subroutine verify_file_presence(filerad,orbs,iformat,nproc)
  use module_base
  use module_types
  use module_interfaces
  implicit none
  integer, intent(in) :: nproc
  character(len=*), intent(in) :: filerad
  type(orbitals_data), intent(in) :: orbs
  integer, intent(out) :: iformat
  !local variables
  character(len=500) :: filename
  logical :: onefile,allfiles
  integer :: iorb,ispinor,iorb_out,ierr
  
  allfiles=.true.

  !first try with plain files
  loop_plain: do iorb=1,orbs%norbp
     do ispinor=1,orbs%nspinor
        call filename_of_iorb(.false.,trim(filerad),orbs,iorb,ispinor,filename,iorb_out)
        inquire(file=filename,exist=onefile)
        allfiles=allfiles .and. onefile
        if (.not. allfiles) then
           exit loop_plain
        end if
     end do
  end do loop_plain
  !reduce the result among the other processors
  if (nproc > 1) call mpiallred(allfiles,1,MPI_LAND,bigdft_mpi%mpi_comm,ierr)
 
  if (allfiles) then
     iformat=WF_FORMAT_PLAIN
     return
  end if

  !Otherwise  test binary files.
  allfiles = .true.
  loop_binary: do iorb=1,orbs%norbp
     do ispinor=1,orbs%nspinor
        call filename_of_iorb(.true.,trim(filerad),orbs,iorb,ispinor,filename,iorb_out)

        inquire(file=filename,exist=onefile)
        allfiles=allfiles .and. onefile
        if (.not. allfiles) then
           exit loop_binary
        end if

     end do
  end do loop_binary
  !reduce the result among the other processors
  if (nproc > 1) call mpiallred(allfiles,1,MPI_LAND,bigdft_mpi%mpi_comm,ierr)

  if (allfiles) then
     iformat=WF_FORMAT_BINARY
     return
  end if

  !otherwise, switch to normal input guess
  iformat=WF_FORMAT_NONE

end subroutine verify_file_presence

!> Associate to the absolute value of orbital a filename which depends of the k-point and 
!! of the spin sign
subroutine filename_of_iorb(lbin,filename,orbs,iorb,ispinor,filename_out,iorb_out,iiorb)
  use module_base
  use module_types
  implicit none
  character(len=*), intent(in) :: filename
  logical, intent(in) :: lbin
  integer, intent(in) :: iorb,ispinor
  type(orbitals_data), intent(in) :: orbs
  character(len=*) :: filename_out
  integer, intent(out) :: iorb_out
  integer, intent(in), optional :: iiorb
  !local variables
  character(len=1) :: spintype,realimag
  character(len=4) :: f3
  character(len=5) :: f4
  character(len=8) :: completename
  integer :: ikpt
  real(gp) :: spins

  !calculate k-point
  ikpt=orbs%iokpt(iorb)
  write(f3,'(a1,i3.3)') "k", ikpt !not more than 999 kpts

  !see if the wavefunction is real or imaginary
  if(modulo(ispinor,2)==0) then
     realimag='I'
  else
     realimag='R'
  end if
  !calculate the spin sector
  spins=orbs%spinsgn(orbs%isorb+iorb)
  if(orbs%nspinor == 4) then
     if (ispinor <=2) then
        spintype='A'
     else
        spintype='B'
     end if
  else
     if (spins==1.0_gp) then
        spintype='U'
     else
        spintype='D'
     end if
  end if
  !no spin polarization if nspin=1
  if (orbs%nspin==1) spintype='N'

  !calculate the actual orbital value
  iorb_out=iorb+orbs%isorb-(ikpt-1)*orbs%norb
  if(present(iiorb)) iorb_out = iiorb
  !purge the value from the spin sign
  if (spins==-1.0_gp) iorb_out=iorb_out-orbs%norbu

  !value of the orbital 
  write(f4,'(a1,i4.4)') "b", iorb_out

  !complete the information in the name of the orbital
  completename='-'//f3//'-'//spintype//realimag
  if (lbin) then
     filename_out = trim(filename)//completename//".bin."//f4
  else
     filename_out = trim(filename)//completename//"."//f4
  end if

  !print *,'filename: ',filename_out
end subroutine filename_of_iorb


!> Associate to the absolute value of orbital a filename which depends of the k-point and 
!! of the spin sign
subroutine open_filename_of_iorb(unitfile,lbin,filename,orbs,iorb,ispinor,iorb_out,iiorb)
  use module_base
  use module_types
  use module_interfaces, except_this_one => open_filename_of_iorb
  implicit none
  character(len=*), intent(in) :: filename
  logical, intent(in) :: lbin
  integer, intent(in) :: iorb,ispinor,unitfile
  type(orbitals_data), intent(in) :: orbs
  integer, intent(out) :: iorb_out
  integer, intent(in), optional :: iiorb
  !local variables
  character(len=500) :: filename_out

  if(present(iiorb)) then   
     call filename_of_iorb(lbin,filename,orbs,iorb,ispinor,filename_out,iorb_out,iiorb) 
  else
     call filename_of_iorb(lbin,filename,orbs,iorb,ispinor,filename_out,iorb_out)
  end if
  if (lbin) then
     open(unit=unitfile,file=trim(filename_out),status='unknown',form="unformatted")
  else
     open(unit=unitfile,file=trim(filename_out),status='unknown')
  end if

end subroutine open_filename_of_iorb

!>   Write all my wavefunctions in files by calling writeonewave
subroutine writemywaves(iproc,filename,iformat,orbs,n1,n2,n3,hx,hy,hz,at,rxyz,wfd,psi)
  use module_types
  use module_base
  use module_interfaces, except_this_one => writeonewave
  implicit none
  integer, intent(in) :: iproc,n1,n2,n3,iformat
  real(gp), intent(in) :: hx,hy,hz
  type(atoms_data), intent(in) :: at
  type(orbitals_data), intent(in) :: orbs
  type(wavefunctions_descriptors), intent(in) :: wfd
  real(gp), dimension(3,at%nat), intent(in) :: rxyz
  real(wp), dimension(wfd%nvctr_c+7*wfd%nvctr_f,orbs%nspinor,orbs%norbp), intent(in) :: psi
  character(len=*), intent(in) :: filename
  !Local variables
  integer :: ncount1,ncount_rate,ncount_max,iorb,ncount2,iorb_out,ispinor
  real(kind=4) :: tr0,tr1
  real(kind=8) :: tel

  if (iproc == 0) write(*,"(1x,A,A,a)") "Write wavefunctions to file: ", trim(filename),'.*'
  if (iformat == WF_FORMAT_ETSF) then
     call write_waves_etsf(iproc,filename,orbs,n1,n2,n3,hx,hy,hz,at,rxyz,wfd,psi)
  else
     call cpu_time(tr0)
     call system_clock(ncount1,ncount_rate,ncount_max)

     ! Plain BigDFT files.
     do iorb=1,orbs%norbp
        do ispinor=1,orbs%nspinor
           call open_filename_of_iorb(99,(iformat == WF_FORMAT_BINARY),filename, &
                & orbs,iorb,ispinor,iorb_out)           
           call writeonewave(99,(iformat == WF_FORMAT_PLAIN),iorb_out,n1,n2,n3,hx,hy,hz, &
                at%nat,rxyz,wfd%nseg_c,wfd%nvctr_c,wfd%keygloc(1,1),wfd%keyvloc(1),  & 
                wfd%nseg_f,wfd%nvctr_f,wfd%keygloc(1,wfd%nseg_c+1),wfd%keyvloc(wfd%nseg_c+1), & 
                psi(1,ispinor,iorb),psi(wfd%nvctr_c+1,ispinor,iorb), &
                orbs%eval(iorb+orbs%isorb))
           close(99)
        end do
     enddo

     call cpu_time(tr1)
     call system_clock(ncount2,ncount_rate,ncount_max)
     tel=dble(ncount2-ncount1)/dble(ncount_rate)
     write(*,'(a,i4,2(1x,1pe10.3))') '- WRITE WAVES TIME',iproc,tr1-tr0,tel
     !write(*,'(a,1x,i0,a)') '- iproc',iproc,' finished writing waves'
  end if

END SUBROUTINE writemywaves

subroutine read_wave_to_isf(lstat, filename, ln, iorbp, hx, hy, hz, &
     & n1, n2, n3, nspinor, psiscf)
  use module_base
  use module_types
  use module_interfaces, except_this_one => read_wave_to_isf

  implicit none

  integer, intent(in) :: ln
  character, intent(in) :: filename(ln)
  integer, intent(in) :: iorbp
  integer, intent(out) :: n1, n2, n3, nspinor
  real(gp), intent(out) :: hx, hy, hz
  real(wp), dimension(:,:,:,:), pointer :: psiscf
  logical, intent(out) :: lstat

  character(len = 1024) :: filename_
  integer :: wave_format_from_filename, iformat, i
  
  write(filename_, "(A)") " "
  do i = 1, ln, 1
     filename_(i:i) = filename(i)
  end do

  ! Find format from name.
  iformat = wave_format_from_filename(1, trim(filename_))

  ! Call proper wraping routine.
  if (iformat == WF_FORMAT_ETSF) then
     call readwavetoisf_etsf(lstat, trim(filename_), iorbp, hx, hy, hz, &
          & n1, n2, n3, nspinor, psiscf)
  else
     call readwavetoisf(lstat, trim(filename_), (iformat == WF_FORMAT_PLAIN), hx, hy, hz, &
          & n1, n2, n3, nspinor, psiscf)
  end if
END SUBROUTINE read_wave_to_isf

subroutine free_wave_to_isf(psiscf)
  use module_base
  implicit none
  real(wp), dimension(:,:,:,:), pointer :: psiscf

  integer :: i_all, i_stat

  i_all=-product(shape(psiscf))*kind(psiscf)
  deallocate(psiscf,stat=i_stat)
  call memocc(i_stat,i_all,'psiscf',"free_wave_to_isf_etsf")
END SUBROUTINE free_wave_to_isf

subroutine read_wave_descr(lstat, filename, ln, &
     & norbu, norbd, iorbs, ispins, nkpt, ikpts, nspinor, ispinor)
  use module_types
  implicit none
  integer, intent(in) :: ln
  character, intent(in) :: filename(ln)
  integer, intent(out) :: norbu, norbd, nkpt, nspinor
  integer, intent(out) :: iorbs, ispins, ikpts, ispinor
  logical, intent(out) :: lstat

  character(len = 1024) :: filename_
  integer :: wave_format_from_filename, iformat, i
  character(len = 1024) :: testf
  
  write(filename_, "(A)") " "
  do i = 1, ln, 1
     filename_(i:i) = filename(i)
  end do

  ! Find format from name.
  iformat = wave_format_from_filename(1, trim(filename_))

  ! Call proper wraping routine.
  if (iformat == WF_FORMAT_ETSF) then
     call readwavedescr_etsf(lstat, trim(filename_), norbu, norbd, nkpt, nspinor)
     iorbs = 1
     ispins = 1
     ikpts = 1
     ispinor = 1
  else
     call readwavedescr(lstat, trim(filename_), iorbs, ispins, ikpts, ispinor, nspinor, testf)
     norbu = 0
     norbd = 0
     nkpt = 0
  end if
END SUBROUTINE read_wave_descr


subroutine writeonewave_linear(unitwf,useFormattedOutput,iorb,n1,n2,n3,hx,hy,hz,locregCenter,&
     locrad,confPotOrder,confPotprefac,nat,rxyz, nseg_c,nvctr_c,keyg_c,keyv_c,  &
     nseg_f,nvctr_f,keyg_f,keyv_f, &
     psi_c,psi_f,eval,onwhichatom)
  use module_base
  implicit none
  logical, intent(in) :: useFormattedOutput
  integer, intent(in) :: unitwf,iorb,n1,n2,n3,nat,nseg_c,nvctr_c,nseg_f,nvctr_f,confPotOrder
  real(gp), intent(in) :: hx,hy,hz,locrad,confPotprefac
  real(wp), intent(in) :: eval
  integer, dimension(nseg_c), intent(in) :: keyv_c
  integer, dimension(nseg_f), intent(in) :: keyv_f
  integer, dimension(2,nseg_c), intent(in) :: keyg_c
  integer, dimension(2,nseg_f), intent(in) :: keyg_f
  real(wp), dimension(nvctr_c), intent(in) :: psi_c
  real(wp), dimension(7,nvctr_f), intent(in) :: psi_f
  real(gp), dimension(3,nat), intent(in) :: rxyz
  real(gp), dimension(3), intent(in) :: locregCenter
  integer, intent(in) :: onwhichatom
  !local variables
  integer :: iat,jj,j0,j1,ii,i0,i1,i2,i3,i,iseg,j
  real(wp) :: tt,t1,t2,t3,t4,t5,t6,t7

  if (useFormattedOutput) then
     write(unitwf,*) iorb,eval
     write(unitwf,*) hx,hy,hz
     write(unitwf,*) n1,n2,n3
     write(unitwf,*) locregCenter(1),locregCenter(2),locregCenter(3),onwhichatom,locrad,&
          confPotOrder,confPotprefac
     write(unitwf,*) nat
     do iat=1,nat
     write(unitwf,'(3(1x,e24.17))') (rxyz(j,iat),j=1,3)
     enddo
     write(unitwf,*) nvctr_c, nvctr_f
  else
     write(unitwf) iorb,eval
     write(unitwf) hx,hy,hz
     write(unitwf) n1,n2,n3
     write(unitwf) locregCenter(1),locregCenter(2),locregCenter(3),onwhichatom,locrad,&
          confPotOrder,confPotprefac
     write(unitwf) nat
     do iat=1,nat
     write(unitwf) (rxyz(j,iat),j=1,3)
     enddo
     write(unitwf) nvctr_c, nvctr_f
  end if

  ! coarse part
  do iseg=1,nseg_c
     jj=keyv_c(iseg)
     j0=keyg_c(1,iseg)
     j1=keyg_c(2,iseg)
     ii=j0-1
     i3=ii/((n1+1)*(n2+1))
     ii=ii-i3*(n1+1)*(n2+1)
     i2=ii/(n1+1)
     i0=ii-i2*(n1+1)
     i1=i0+j1-j0
     do i=i0,i1
        tt=psi_c(i-i0+jj)
        if (useFormattedOutput) then
           write(unitwf,'(3(i4),1x,e19.12)') i,i2,i3,tt
        else
           write(unitwf) i,i2,i3,tt
        end if
     enddo
  enddo

  ! fine part
  do iseg=1,nseg_f
     jj=keyv_f(iseg)
     j0=keyg_f(1,iseg)
     j1=keyg_f(2,iseg)
     ii=j0-1
     i3=ii/((n1+1)*(n2+1))
     ii=ii-i3*(n1+1)*(n2+1)
     i2=ii/(n1+1)
     i0=ii-i2*(n1+1)
     i1=i0+j1-j0
     do i=i0,i1
        t1=psi_f(1,i-i0+jj)
        t2=psi_f(2,i-i0+jj)
        t3=psi_f(3,i-i0+jj)
        t4=psi_f(4,i-i0+jj)
        t5=psi_f(5,i-i0+jj)
        t6=psi_f(6,i-i0+jj)
        t7=psi_f(7,i-i0+jj)
        if (useFormattedOutput) then
           write(unitwf,'(3(i4),7(1x,e17.10))') i,i2,i3,t1,t2,t3,t4,t5,t6,t7
        else
           write(unitwf) i,i2,i3,t1,t2,t3,t4,t5,t6,t7
        end if
     enddo
  enddo

  if (verbose >= 2) write(*,'(1x,i0,a)') iorb,'th wavefunction written'

END SUBROUTINE writeonewave_linear


subroutine writeLinearCoefficients(unitwf,useFormattedOutput,n1,n2,n3,hx,hy,hz,nat,rxyz,&
           norb,ntmb,nvctr_c,nvctr_f,coeff,eval)
  use module_base
  implicit none
  logical, intent(in) :: useFormattedOutput
  integer, intent(in) :: unitwf,norb,n1,n2,n3,nat,ntmb,nvctr_c,nvctr_f
  real(gp), intent(in) :: hx,hy,hz
  real(wp), dimension(ntmb,norb), intent(in) :: coeff
  real(wp), dimension(norb), intent(in) :: eval
  real(gp), dimension(3,nat), intent(in) :: rxyz
  !local variables
  integer :: iat,i,j,iorb
  real(wp) :: tt

  ! Write the Header
  if (useFormattedOutput) then
     write(unitwf,*) norb,ntmb
     write(unitwf,*) hx,hy,hz
     write(unitwf,*) n1,n2,n3
     write(unitwf,*) nat
     do iat=1,nat
     write(unitwf,'(3(1x,e24.17))') (rxyz(j,iat),j=1,3)
     enddo
     write(unitwf,*) nvctr_c, nvctr_f
     do iorb=1,norb
     write(unitwf,*) iorb,eval(iorb)
     enddo
  else
     write(unitwf) norb, ntmb
     write(unitwf) hx,hy,hz
     write(unitwf) n1,n2,n3
     write(unitwf) nat
     do iat=1,nat
     write(unitwf) (rxyz(j,iat),j=1,3)
     enddo
     write(unitwf) nvctr_c, nvctr_f
     do iorb=1,norb
     write(unitwf) iorb,eval(iorb)
     enddo
  end if

  ! Now write the coefficients
  do i = 1, norb
     do j = 1, ntmb
          tt = coeff(j,i)
          if (useFormattedOutput) then
             write(unitwf,'(2(i6,1x),e19.12)') i,j,tt
          else
             write(unitwf) i,j,tt
          end if
     end do
  end do  
  if (verbose >= 2) write(*,'(1x,a)') 'Wavefunction coefficients written'

END SUBROUTINE writeLinearCoefficients


!>   Write all my wavefunctions in files by calling writeonewave                                                                                            
subroutine writemywaves_linear(iproc,filename,iformat,Lzd,orbs,norb,hx,hy,hz,at,rxyz,psi,coeff,eval)
  use module_types
  use module_base
  use module_interfaces, except_this_one => writeonewave
  implicit none
  integer, intent(in) :: iproc,iformat
  integer, intent(in) :: norb   !< number of orbitals, not basis functions
  real(gp), intent(in) :: hx,hy,hz
  type(atoms_data), intent(in) :: at
  type(orbitals_data), intent(in) :: orbs         !< orbs describing the basis functions
  type(local_zone_descriptors), intent(in) :: Lzd
  real(gp), dimension(3,at%nat), intent(in) :: rxyz
  real(wp), dimension(max(orbs%npsidim_orbs,orbs%npsidim_comp)), intent(in) :: psi  ! Should be the real linear dimension and not the global
  real(wp), dimension(orbs%norb,norb), intent(in) :: coeff
  real(wp), dimension(norb), intent(in) :: eval
  character(len=*), intent(in) :: filename
  !Local variables
  integer :: ncount1,ncount_rate,ncount_max,iorb,ncount2,iorb_out,ispinor,ilr,shift,ii,iat
  integer :: jorb,jlr
  real(kind=4) :: tr0,tr1
  real(kind=8) :: tel

  if (iproc == 0) write(*,"(1x,A,A,a)") "Write wavefunctions to file: ", trim(filename),'.*'

  if (iformat == WF_FORMAT_ETSF) then
      stop 'Linear scaling with ETSF writing not implemented yet'
!     call write_waves_etsf(iproc,filename,orbs,n1,n2,n3,hx,hy,hz,at,rxyz,wfd,psi)
  else
     call cpu_time(tr0)
     call system_clock(ncount1,ncount_rate,ncount_max)

     ! Write the TMBs in the Plain BigDFT files.
     ! Use same ordering as posinp and llr generation
     ii = 0
     do iat = 1, at%nat
        do iorb=1,orbs%norbp
           if(iat == orbs%onwhichatom(iorb+orbs%isorb)) then
              shift = 1
              do jorb = 1, iorb-1 
                 jlr = orbs%inwhichlocreg(jorb+orbs%isorb)
                 shift = shift + Lzd%Llr(jlr)%wfd%nvctr_c+7*Lzd%Llr(jlr)%wfd%nvctr_f
              end do
              ii = ii + 1
              ilr = orbs%inwhichlocreg(iorb+orbs%isorb)
              do ispinor=1,orbs%nspinor
                 call open_filename_of_iorb(99,(iformat == WF_FORMAT_BINARY),filename, &
                    & orbs,iorb,ispinor,iorb_out)
                 call writeonewave_linear(99,(iformat == WF_FORMAT_PLAIN),iorb_out,Lzd%Glr%d%n1,Lzd%Glr%d%n2,Lzd%Glr%d%n3,&
                    & Lzd%hgrids(1),Lzd%hgrids(2),Lzd%hgrids(3), &
                    & Lzd%Llr(ilr)%locregCenter,Lzd%Llr(ilr)%locrad, 4, 0.0d0, &  !put here the real potentialPrefac and Order
                    & at%nat,rxyz,Lzd%Llr(ilr)%wfd%nseg_c,Lzd%Llr(ilr)%wfd%nvctr_c,&
                    & Lzd%Llr(ilr)%wfd%keyglob,Lzd%Llr(ilr)%wfd%keyvglob, &
                    & Lzd%Llr(ilr)%wfd%nseg_f,Lzd%Llr(ilr)%wfd%nvctr_f,&
                    & Lzd%Llr(ilr)%wfd%keyglob(1,Lzd%Llr(ilr)%wfd%nseg_c+1), &
                    & Lzd%Llr(ilr)%wfd%keyvglob(Lzd%Llr(ilr)%wfd%nseg_c+1), &
                    & psi(shift),psi(Lzd%Llr(ilr)%wfd%nvctr_c+shift),orbs%eval(iorb+orbs%isorb),&
                    & orbs%onwhichatom(iorb+orbs%isorb))
                 close(99)
              end do
           end if
        enddo
     end do

    ! Now write the coefficients to file
    ! Must be careful, the orbs%norb is the number of basis functions
    ! while the norb is the number of orbitals.
    if(iproc == 0) then
      if(iformat == WF_FORMAT_PLAIN) then
         open(99, file=filename//'_coeff.bin', status='unknown',form='formatted')
      else
         open(99, file=filename//'_coeff.bin', status='unknown',form='unformatted')
      end if
      call writeLinearCoefficients(99,(iformat == WF_FORMAT_PLAIN),Lzd%Glr%d%n1,Lzd%Glr%d%n2,Lzd%Glr%d%n3,&
           Lzd%hgrids(1),Lzd%hgrids(2),Lzd%hgrids(3),at%nat,rxyz,norb,orbs%norb,Lzd%Glr%wfd%nvctr_c,Lzd%Glr%wfd%nvctr_f,&
           coeff,eval)
      close(99)
    end if
     call cpu_time(tr1)
     call system_clock(ncount2,ncount_rate,ncount_max)
     tel=dble(ncount2-ncount1)/dble(ncount_rate)
     write(*,'(a,i4,2(1x,1pe10.3))') '- WRITE WAVES TIME',iproc,tr1-tr0,tel
     !write(*,'(a,1x,i0,a)') '- iproc',iproc,' finished writing waves'
  end if

END SUBROUTINE writemywaves_linear


subroutine readonewave_linear(unitwf,useFormattedInput,iorb,iproc,n1,n2,n3,&
     & hx,hy,hz,at,wfd,rxyz_old,rxyz,locrad,locregCenter,confPotOrder,&
     & confPotprefac,psi,eval,psifscf,onwhichatom,lr)
  use module_base
  use module_types
  use internal_io
  use module_interfaces
  implicit none
  logical, intent(in) :: useFormattedInput
  integer, intent(in) :: unitwf,iorb,iproc,n1,n2,n3
  type(wavefunctions_descriptors), intent(in) :: wfd
  type(atoms_data), intent(in) :: at
  real(gp), intent(in) :: hx,hy,hz
  real(gp), dimension(3,at%nat), intent(in) :: rxyz
  integer, intent(out) :: confPotOrder
  real(gp), intent(out) :: locrad, confPotprefac
  real(wp), intent(out) :: eval
  real(gp), dimension(3), intent(out) :: locregCenter
  real(gp), dimension(3,at%nat), intent(out) :: rxyz_old
  real(wp), dimension(wfd%nvctr_c+7*wfd%nvctr_f), intent(out) :: psi
  real(wp), dimension(*), intent(out) :: psifscf !this supports different BC
  integer, dimension(*), intent(in) :: onwhichatom
  type(locreg_descriptors), intent(in) :: lr

  !local variables
  character(len=*), parameter :: subname='readonewave_linear'
  character(len = 256) :: error
  logical :: perx,pery,perz,lstat
  integer :: iorb_old,n1_old,n2_old,n3_old,iat,nvctr_c_old,nvctr_f_old,i_all,iiat
  integer :: i1,i2,i3,iel,i_stat,iall,onwhichatom_tmp
  real(gp) :: tx,ty,tz,displ,hx_old,hy_old,hz_old,mindist
  real(gp) :: tt,t1,t2,t3,t4,t5,t6,t7
  real(wp), dimension(:,:,:,:,:,:), allocatable :: psigold
  character(len=12) :: orbname
  !write(*,*) 'INSIDE readonewave'

  call io_read_descr_linear(unitwf, useFormattedInput, iorb_old, eval, n1_old, n2_old, n3_old, &
       & hx_old, hy_old, hz_old, lstat, error, nvctr_c_old, nvctr_f_old, rxyz_old, at%nat,&
       & locrad, locregCenter, confPotOrder, confPotprefac, onwhichatom_tmp)

  if (.not. lstat) call io_error(trim(error))
  if (iorb_old /= iorb) stop 'readonewave_linear'

  !conditions for periodicity in the three directions
  perx=(at%geocode /= 'F')
  pery=(at%geocode == 'P')
  perz=(at%geocode /= 'F')

  tx=0.0_gp
  ty=0.0_gp
  tz=0.0_gp
  do iat=1,at%nat
     tx=tx+mindist(perx,at%alat1,rxyz(1,iat),rxyz_old(1,iat))**2
     ty=ty+mindist(pery,at%alat2,rxyz(2,iat),rxyz_old(2,iat))**2
     tz=tz+mindist(perz,at%alat3,rxyz(3,iat),rxyz_old(3,iat))**2
  enddo
  displ=sqrt(tx+ty+tz)

  if (hx_old == hx .and. hy_old == hy .and. hz_old == hz .and.&
       n1_old == n1  .and. n2_old == n2 .and. n3_old == n3 .and. displ <= 1.d-3) then

     if (iproc == 0) write(*,*) 'wavefunctions need NO reformatting'
     call read_psi_compress(unitwf, useFormattedInput, nvctr_c_old, nvctr_f_old, psi, lstat, error)
     if (.not. lstat) call io_error(trim(error))

  else

     if (iproc == 0 .and. iorb == 1) then
        write(*,*) 'wavefunctions need reformatting'
        if (hx_old /= hx .or. hy_old /= hy .or. hz_old /= hz) write(*,"(1x,A,6F14.10)") &
             'because hgrid_old /= hgrid',hx_old,hy_old,hz_old,hx,hy,hz
        if (n1_old /= n1  .or. n2_old /= n2 .or. n3_old /= n3 ) &
             write(*,*) 'because cell size has changed',n1_old,n1,n2_old,n2,n3_old,n3
        if (displ > 1.d-3 ) write(*,*) 'large displacement of molecule',displ
     end if

!     ! NOT SURE YET WHAT SHOULD BE DONE FOR LINEAR CASE, so just stop
!     if(iproc==0) write(*,*) 'This is forbidden for now in linear case!'
!     call mpi_finalize(i_all)
!     stop 
!needs fixing below
     ! also need to add derivative functions, which needs orbs and lzd
     allocate(psigold(0:n1_old,2,0:n2_old,2,0:n3_old,2+ndebug),stat=i_stat)
     call memocc(i_stat,psigold,'psigold',subname)

     call razero(8*(n1_old+1)*(n2_old+1)*(n3_old+1),psigold)
     do iel=1,nvctr_c_old
        if (useFormattedInput) then
           read(unitwf,*) i1,i2,i3,tt
        else
           read(unitwf) i1,i2,i3,tt
        end if
        psigold(i1,1,i2,1,i3,1)=tt
     enddo
     do iel=1,nvctr_f_old
        if (useFormattedInput) then
           read(unitwf,*) i1,i2,i3,t1,t2,t3,t4,t5,t6,t7
        else
           read(unitwf) i1,i2,i3,t1,t2,t3,t4,t5,t6,t7
        end if
        psigold(i1,2,i2,1,i3,1)=t1
        psigold(i1,1,i2,2,i3,1)=t2
        psigold(i1,2,i2,2,i3,1)=t3
        psigold(i1,1,i2,1,i3,2)=t4
        psigold(i1,2,i2,1,i3,2)=t5
        psigold(i1,1,i2,2,i3,2)=t6
        psigold(i1,2,i2,2,i3,2)=t7
     enddo

     ! onwhichatom should be replaced with that read from file for consistent reordering -
     ! ordering can change if positions have moved - already done previously
     !print*,onwhichatom(iorb), onwhichatom_tmp
     !onwhichatom(iorb) = onwhichatom_tmp
     iiat=onwhichatom(iorb)

     !call reformat_one_supportfunction here to be consistent with cubic
     call reformat_one_supportfunction(iiat,displ,wfd,at,hx_old,hy_old,hz_old, & !n(m)
          n1_old,n2_old,n3_old,rxyz_old,psigold,hx,hy,hz,&
          n1,n2,n3,rxyz,psifscf,psi)

     i_all=-product(shape(psigold))*kind(psigold)
     deallocate(psigold,stat=i_stat)
     call memocc(i_stat,i_all,'psigold',subname)

  endif

  ! DEBUG
  !write(orbname,*) iorb
  !call plot_wf(trim(adjustl(orbname)),1,at,1.0_dp,lr,hx,hy,hz,rxyz,psi)
  ! END DEBUG 

END SUBROUTINE readonewave_linear                                                     

subroutine io_read_descr_linear(unitwf, formatted, iorb_old, eval, n1_old, n2_old, n3_old, &
       & hx_old, hy_old, hz_old, lstat, error, nvctr_c_old, nvctr_f_old, rxyz_old, nat, &
       & locrad, locregCenter, confPotOrder, confPotprefac,onwhichatom)
    use module_base
    use module_types
    use internal_io
    implicit none

    integer, intent(in) :: unitwf
    logical, intent(in) :: formatted
    integer, intent(out) :: iorb_old
    integer, intent(out) :: n1_old, n2_old, n3_old
    real(gp), intent(out) :: hx_old, hy_old, hz_old
    logical, intent(out) :: lstat
    real(wp), intent(out) :: eval
    integer, intent(out) :: confPotOrder
    real(gp), intent(out) :: locrad, confPotprefac
    real(gp), dimension(3), intent(out) :: locregCenter
    character(len =256), intent(out) :: error
    ! Optional arguments
    integer, intent(out), optional :: nvctr_c_old, nvctr_f_old
    integer, intent(in), optional :: nat
    real(gp), dimension(:,:), intent(out), optional :: rxyz_old
    integer, intent(out) :: onwhichatom

    character(len = *), parameter :: subname = "io_read_descr_linear"
    integer :: i, iat, i_stat, nat_
    real(gp) :: rxyz(3)

    lstat = .false.
    write(error, "(A)") "cannot read psi description."
    if (formatted) then
       read(unitwf,*,iostat=i_stat) iorb_old,eval
       if (i_stat /= 0) return

       read(unitwf,*,iostat=i_stat) hx_old,hy_old,hz_old
       if (i_stat /= 0) return

       read(unitwf,*,iostat=i_stat) n1_old,n2_old,n3_old
       if (i_stat /= 0) return

       read(unitwf,*,iostat=i_stat) (locregCenter(i),i=1,3),onwhichatom,&
            locrad,confPotOrder, confPotprefac
       if (i_stat /= 0) return
       !write(*,*) 'reading ',nat,' atomic positions' !*

       if (present(nat) .And. present(rxyz_old)) then
          read(unitwf,*,iostat=i_stat) nat_
          if (i_stat /= 0) return
          ! Sanity check
          if (size(rxyz_old, 2) /= nat) stop "Mismatch in coordinate array size."
          if (nat_ /= nat) stop "Mismatch in coordinate array size."
          do iat=1,nat
             read(unitwf,*,iostat=i_stat) (rxyz_old(i,iat),i=1,3)
             if (i_stat /= 0) return

          enddo
       else
          read(unitwf,*,iostat=i_stat) nat_
          if (i_stat /= 0) return
          do iat=1,nat_
             read(unitwf,*,iostat=i_stat)
             if (i_stat /= 0) return
          enddo
       end if
       if (present(nvctr_c_old) .and. present(nvctr_f_old)) then
          read(unitwf,*,iostat=i_stat) nvctr_c_old, nvctr_f_old
          if (i_stat /= 0) return
       else
          read(unitwf,*,iostat=i_stat) i, iat
          if (i_stat /= 0) return
       end if
    else
       read(unitwf,iostat=i_stat) iorb_old,eval
       if (i_stat /= 0) return

       read(unitwf,iostat=i_stat) hx_old,hy_old,hz_old
       if (i_stat /= 0) return
       read(unitwf,iostat=i_stat) n1_old,n2_old,n3_old
       if (i_stat /= 0) return
       read(unitwf,iostat=i_stat) (locregCenter(i),i=1,3),onwhichatom,&
            locrad,confPotOrder, confPotprefac
       if (i_stat /= 0) return
       if (present(nat) .And. present(rxyz_old)) then
          read(unitwf,iostat=i_stat) nat_
          if (i_stat /= 0) return
          ! Sanity check
          if (size(rxyz_old, 2) /= nat) stop "Mismatch in coordinate array size." 
          if (nat_ /= nat) stop "Mismatch in coordinate array size."
          do iat=1,nat
             read(unitwf,iostat=i_stat)(rxyz_old(i,iat),i=1,3)
             if (i_stat /= 0) return
          enddo
       else
          read(unitwf,iostat=i_stat) nat_
          if (i_stat /= 0) return
          do iat=1,nat_
             read(unitwf,iostat=i_stat) rxyz
             if (i_stat /= 0) return
          enddo
       end if
       if (present(nvctr_c_old) .and. present(nvctr_f_old)) then
          read(unitwf,iostat=i_stat) nvctr_c_old, nvctr_f_old
          if (i_stat /= 0) return
       else
          read(unitwf,iostat=i_stat) i, iat
          if (i_stat /= 0) return
       end if
    end if
    lstat = .true.

END SUBROUTINE io_read_descr_linear

subroutine io_read_descr_coeff(unitwf, formatted, norb_old, ntmb_old, n1_old, n2_old, n3_old, &
       & hx_old, hy_old, hz_old, lstat, error, nvctr_c_old, nvctr_f_old, rxyz_old, nat)
    use module_base
    use module_types
    use internal_io
    implicit none
    integer, intent(in) :: unitwf
    logical, intent(in) :: formatted
    integer, intent(out) :: norb_old, ntmb_old
    integer, intent(out) :: n1_old, n2_old, n3_old
    real(gp), intent(out) :: hx_old, hy_old, hz_old
    logical, intent(out) :: lstat
    character(len =256), intent(out) :: error
    ! Optional arguments
    integer, intent(out), optional :: nvctr_c_old, nvctr_f_old
    integer, intent(in), optional :: nat
    real(gp), dimension(:,:), intent(out), optional :: rxyz_old

    character(len = *), parameter :: subname = "io_read_descr_linear"
    integer :: i, iat, i_stat, nat_
    real(gp) :: rxyz(3)

    lstat = .false.
    write(error, "(A)") "cannot read psi description."
    if (formatted) then
       read(unitwf,*,iostat=i_stat) norb_old, ntmb_old
       if (i_stat /= 0) return
       read(unitwf,*,iostat=i_stat) hx_old,hy_old,hz_old
       if (i_stat /= 0) return
       read(unitwf,*,iostat=i_stat) n1_old,n2_old,n3_old
       if (i_stat /= 0) return
       !write(*,*) 'reading ',nat,' atomic positions'
       if (present(nat) .And. present(rxyz_old)) then
          read(unitwf,*,iostat=i_stat) nat_
          if (i_stat /= 0) return
          ! Sanity check
          if (size(rxyz_old, 2) /= nat) stop "Mismatch in coordinate array size."
          if (nat_ /= nat) stop "Mismatch in coordinate array size."
          do iat=1,nat
             read(unitwf,*,iostat=i_stat) (rxyz_old(i,iat),i=1,3)
             if (i_stat /= 0) return
          enddo
       else
          read(unitwf,*,iostat=i_stat) nat_
          if (i_stat /= 0) return
          do iat=1,nat_
             read(unitwf,*,iostat=i_stat)
             if (i_stat /= 0) return
          enddo
       end if
       if (present(nvctr_c_old) .and. present(nvctr_f_old)) then
          read(unitwf,*,iostat=i_stat) nvctr_c_old, nvctr_f_old
          if (i_stat /= 0) return
       else
          read(unitwf,*,iostat=i_stat) i, iat
          if (i_stat /= 0) return
       end if
    else
       read(unitwf,iostat=i_stat) norb_old, ntmb_old
       if (i_stat /= 0) return
       read(unitwf,iostat=i_stat) hx_old,hy_old,hz_old
       if (i_stat /= 0) return
       read(unitwf,iostat=i_stat) n1_old,n2_old,n3_old
       if (i_stat /= 0) return
       if (present(nat) .And. present(rxyz_old)) then
          read(unitwf,iostat=i_stat) nat_
          if (i_stat /= 0) return
          ! Sanity check
          if (size(rxyz_old, 2) /= nat) stop "Mismatch in coordinate array size." 
          if (nat_ /= nat) stop "Mismatch in coordinate array size."
          do iat=1,nat
             read(unitwf,iostat=i_stat)(rxyz_old(i,iat),i=1,3)
             if (i_stat /= 0) return
          enddo
       else
          read(unitwf,iostat=i_stat) nat_
          if (i_stat /= 0) return
          do iat=1,nat_
             read(unitwf,iostat=i_stat) rxyz
             if (i_stat /= 0) return
          enddo
       end if
       if (present(nvctr_c_old) .and. present(nvctr_f_old)) then
          read(unitwf,iostat=i_stat) nvctr_c_old, nvctr_f_old
          if (i_stat /= 0) return
       else
          read(unitwf,iostat=i_stat) i, iat
          if (i_stat /= 0) return
       end if
    end if
    lstat = .true.
END SUBROUTINE io_read_descr_coeff


subroutine read_coeff_minbasis(unitwf,useFormattedInput,iproc,n1,n2,n3,norb,ntmb,&
     & hx,hy,hz,at,rxyz_old,rxyz,coeff,eval,norb_change)
  use module_base
  use module_types
  use internal_io
  use module_interfaces
  implicit none
  logical, intent(in) :: useFormattedInput
  integer, intent(in) :: unitwf,iproc,n1,n2,n3,norb,ntmb
  type(atoms_data), intent(in) :: at
  real(gp), intent(in) :: hx,hy,hz
  real(gp), dimension(3,at%nat), intent(in) :: rxyz
  real(gp), dimension(3,at%nat), intent(out) :: rxyz_old
  real(wp), dimension(ntmb,norb), intent(out) :: coeff
  real(wp), dimension(norb), intent(out) :: eval
  logical, intent(out) :: norb_change
  !local variables
  character(len=*), parameter :: subname='readonewave_linear'
  character(len = 256) :: error
  logical :: perx,pery,perz,lstat
  integer :: norb_old,n1_old,n2_old,n3_old,iat,nvctr_c_old,nvctr_f_old,i_stat,i_all
  integer :: ntmb_old, i1, i2,i,j,iorb,iorb_old
  real(wp) :: tt
  real(gp) :: tx,ty,tz,displ,hx_old,hy_old,hz_old,mindist

  norb_change = .false.

  !write(*,*) 'INSIDE readonewave'
  call io_read_descr_coeff(unitwf, useFormattedInput, norb_old, ntmb_old, n1_old, n2_old, n3_old, &
       & hx_old, hy_old, hz_old, lstat, error, nvctr_c_old, nvctr_f_old, rxyz_old, at%nat)
  if (.not. lstat) call io_error(trim(error))

  !conditions for periodicity in the three directions
  perx=(at%geocode /= 'F')
  pery=(at%geocode == 'P')
  perz=(at%geocode /= 'F')

  tx=0.0_gp
  ty=0.0_gp
  tz=0.0_gp
  do iat=1,at%nat
     tx=tx+mindist(perx,at%alat1,rxyz(1,iat),rxyz_old(1,iat))**2
     ty=ty+mindist(pery,at%alat2,rxyz(2,iat),rxyz_old(2,iat))**2
     tz=tz+mindist(perz,at%alat3,rxyz(3,iat),rxyz_old(3,iat))**2
  enddo
  displ=sqrt(tx+ty+tz)

  if (norb == norb_old) then

     ! read the eigenvalues
     if (useFormattedInput) then
        do iorb=1,norb
           read(unitwf,*,iostat=i_stat) iorb_old,eval(iorb)
           if (iorb_old /= iorb) stop 'read_coeff_minbasis'
        enddo
     else 
        do iorb=1,norb
           read(unitwf,iostat=i_stat) iorb_old,eval(iorb)
           if (iorb_old /= iorb) stop 'read_coeff_minbasis'
        enddo
        if (i_stat /= 0) stop 'Problem reading the coefficients'
     end if

     if (iproc == 0) write(*,*) 'coefficients need NO reformatting'

     ! Now read the coefficients
     do i = 1, norb
        do j = 1, ntmb
           if (useFormattedInput) then
              read(unitwf,*,iostat=i_stat) i1,i2,tt
           else
              read(unitwf,iostat=i_stat) i1,i2,tt
           end if
           if (i_stat /= 0) stop 'Problem reading the coefficients'
           coeff(j,i) = tt  
        end do
     end do
     if (verbose >= 2) write(*,'(1x,a)') 'Wavefunction coefficients read'
  else
     ! tmbs themselves should be ok, but need to recalculate the coefficients
     if (norb < norb_old) then ! for now if we have too many, just eliminate highest
        ! read the eigenvalues
        if (useFormattedInput) then
           do iorb=1,norb
              read(unitwf,*,iostat=i_stat) iorb_old,eval(iorb)
              if (iorb_old /= iorb) stop 'read_coeff_minbasis'
           enddo
           do iorb=norb+1,norb_old
              read(unitwf,*,iostat=i_stat) iorb_old,tt
           end do
        else 
           do iorb=1,norb
              read(unitwf,iostat=i_stat) iorb_old,eval(iorb)
              if (iorb_old /= iorb) stop 'read_coeff_minbasis'
           enddo
           do iorb=norb+1,norb_old
              read(unitwf,iostat=i_stat) iorb_old,tt
           end do
           if (i_stat /= 0) stop 'Problem reading the coefficients'
        end if

        if (iproc == 0) write(*,*) 'Eliminating coefficients for highest',norb_old-norb,'states'

        do i = 1, norb_old
           do j = 1, ntmb
              if (useFormattedInput) then
                 read(unitwf,*,iostat=i_stat) i1,i2,tt
              else
                 read(unitwf,iostat=i_stat) i1,i2,tt
              end if
              if (i_stat /= 0) stop 'Problem reading the coefficients'
              if (i <= norb) coeff(j,i) = tt  
           end do
        end do
        if (verbose >= 2) write(*,'(1x,a)') 'Wavefunction coefficients read'    
     else
        if (iproc == 0) write(*,*) 'Not enough orbitals in coefficients, resetting them to the identity'
        ! read the eigenvalues
        if (useFormattedInput) then
           do iorb=1,norb_old
              read(unitwf,*,iostat=i_stat) iorb_old,eval(iorb)
              if (iorb_old /= iorb) stop 'read_coeff_minbasis'
           enddo
           do iorb=norb_old+1,norb
              eval(iorb) = 0.0_dp ! set to zero for lack of better idea
           end do
        else 
           do iorb=1,norb_old
              read(unitwf,iostat=i_stat) iorb_old,eval(iorb)
              if (iorb_old /= iorb) stop 'read_coeff_minbasis'
           enddo
           do iorb=norb_old+1,norb
              eval(iorb) = 0.0_dp ! set to zero for lack of better idea
           end do
           if (i_stat /= 0) stop 'Problem reading the coefficients'
        end if

        norb_change = .true.
        do i = 1, norb_old
           do j = 1, ntmb
              if (useFormattedInput) then
                 read(unitwf,*,iostat=i_stat) i1,i2,tt
              else
                 read(unitwf,iostat=i_stat) i1,i2,tt
              end if
              if (i_stat /= 0) stop 'Problem reading the coefficients'
              coeff(j,i) = tt  
           end do
        end do
        do i=norb_old+1,norb
           do j = 1, ntmb
              if (j==i) then
                 coeff(j,i)=1.0_dp  
              else 
                 coeff(j,i)=0.0_dp  
              end if
           end do
        end do
     end if
  end if



END SUBROUTINE read_coeff_minbasis


!>  Reads wavefunction from file and transforms it properly if hgrid or size of simulation cell                                                                                                                                                                                                                                                                                                                                   
!!  have changed
subroutine readmywaves_linear(iproc,filename,iformat,norb,Lzd,orbs,at,rxyz_old,rxyz,  & 
    psi,coeff,eval,norb_change,orblist)
  use module_base
  use module_types
  use module_interfaces, except_this_one => readmywaves_linear
  implicit none
  integer, intent(in) :: iproc, iformat,norb
  type(orbitals_data), intent(inout) :: orbs  ! orbs related to the basis functions
  type(local_zone_descriptors), intent(in) :: Lzd
  type(atoms_data), intent(in) :: at
  real(gp), dimension(3,at%nat), intent(in) :: rxyz
  real(gp), dimension(3,at%nat), intent(out) :: rxyz_old
  real(wp), dimension(orbs%npsidim_orbs), intent(out) :: psi  
  real(gp), dimension(norb,orbs%norb),intent(out) :: coeff
  real(gp), dimension(norb),intent(out) :: eval
  character(len=*), intent(in) :: filename
  logical, intent(out) :: norb_change
  integer, dimension(orbs%norb), optional :: orblist
  !Local variables
  character(len=*), parameter :: subname='readmywaves_linear'
  integer :: ncount1,ncount_rate,ncount_max,iorb,i_stat,i_all,ncount2
  integer :: iorb_out,ispinor,ilr,ind,nb1,nb2,nb3,n1,n2,n3
  integer :: confPotOrder
  real(gp) :: locrad, confPotprefac
  real(gp), dimension(3) :: locregCenter
  real(kind=4) :: tr0,tr1
  real(kind=8) :: tel
  real(wp), dimension(:,:,:), allocatable :: psifscf
  logical :: perx, pery, perz
  !integer, dimension(orbs%norb) :: orblist2

  call cpu_time(tr0)
  call system_clock(ncount1,ncount_rate,ncount_max)

  if (iformat == WF_FORMAT_ETSF) then
     stop 'Linear scaling with ETSF writing not implemented yet'
     !construct the orblist or use the one in argument
     !do nb1 = 1, orbs%norb
     !orblist2(nb1) = nb1
     !if(present(orblist)) orblist2(nb1) = orblist(nb1) 
     !end do

     !call read_waves_etsf(iproc,filename // ".etsf",orbs,n1,n2,n3,hx,hy,hz,at,rxyz_old,rxyz,  & 
     !     wfd,psi)
  else if (iformat == WF_FORMAT_BINARY .or. iformat == WF_FORMAT_PLAIN) then
     !conditions for periodicity in the three directions
     perx=(at%geocode /= 'F')
     pery=(at%geocode == 'P')
     perz=(at%geocode /= 'F')

     !buffers related to periodicity
     !WARNING: the boundary conditions are not assumed to change between new and old
     call ext_buffers_coarse(perx,nb1)
     call ext_buffers_coarse(pery,nb2)
     call ext_buffers_coarse(perz,nb3)
     n1 = Lzd%Glr%d%n1
     n2 = Lzd%Glr%d%n2
     n3 = Lzd%Glr%d%n3
     allocate(psifscf(-nb1:2*n1+1+nb1,-nb2:2*n2+1+nb2,-nb3:2*n3+1+nb3+ndebug),stat=i_stat)
     call memocc(i_stat,psifscf,'psifscf',subname)
     !allocate(psifscf(1,1,1+ndebug),stat=i_stat)
     !call memocc(i_stat,psifscf,'psifscf',subname)
     ind = 1
     do iorb=1,orbs%norbp!*orbs%nspinor
        ilr = orbs%inwhichlocreg(iorb+orbs%isorb)
        do ispinor=1,orbs%nspinor
           if(present(orblist)) then
              call open_filename_of_iorb(99,(iformat == WF_FORMAT_BINARY),filename, &
                   & orbs,iorb,ispinor,iorb_out, orblist(iorb+orbs%isorb))
           else
              call open_filename_of_iorb(99,(iformat == WF_FORMAT_BINARY),filename, &
                   & orbs,iorb,ispinor,iorb_out)
           end if  
         
           call readonewave_linear(99, (iformat == WF_FORMAT_PLAIN),iorb_out,iproc,&
                Lzd%Glr%d%n1,Lzd%Glr%d%n2,Lzd%Glr%d%n3,Lzd%hgrids(1),Lzd%hgrids(2),&
                Lzd%hgrids(3),at,Lzd%Llr(ilr)%wfd,rxyz_old,rxyz,locrad,locregCenter,&
                confPotOrder,confPotPrefac,psi(ind),orbs%eval(orbs%isorb+iorb),psifscf,&
                orbs%onwhichatom,Lzd%Llr(ilr))

           close(99)
           ind = ind + Lzd%Llr(ilr)%wfd%nvctr_c+7*Lzd%Llr(ilr)%wfd%nvctr_f
        end do

     end do

     i_all=-product(shape(psifscf))*kind(psifscf)
     deallocate(psifscf,stat=i_stat)
     call memocc(i_stat,i_all,'psifscf',subname)

     !Open the coefficient file 
     if(iformat == WF_FORMAT_PLAIN) then
        open(99,file=filename//'_coeff.bin',status='unknown',form='formatted')
     else if(iformat == WF_FORMAT_BINARY) then
        open(99,file=filename//'_coeff.bin',status='unknown',form='unformatted')
     else
        stop 'Coefficient format not implemented'
     end if
     call read_coeff_minbasis(99,(iformat == WF_FORMAT_PLAIN),iproc,Lzd%Glr%d%n1,Lzd%Glr%d%n2,Lzd%Glr%d%n3,norb,orbs%norb,&
     & Lzd%hgrids(1),Lzd%hgrids(2),Lzd%hgrids(3),at,rxyz_old,rxyz,coeff,eval,norb_change)
     close(99)
  else
     write(0,*) "Unknown wavefunction file format from filename."
     stop
  end if

  call cpu_time(tr1)
  call system_clock(ncount2,ncount_rate,ncount_max)
  tel=dble(ncount2-ncount1)/dble(ncount_rate)
  write(*,'(a,i4,2(1x,1pe10.3))') '- READING WAVES TIME',iproc,tr1-tr0,tel
END SUBROUTINE readmywaves_linear


subroutine initialize_linear_from_file(iproc,nproc,filename,iformat,Lzd,orbs,at,rxyz,orblist)
  use module_base
  use module_types
  use module_defs
  use module_interfaces, except_this_one => initialize_linear_from_file
  implicit none
  integer, intent(in) :: iproc, nproc, iformat
  type(orbitals_data), intent(inout) :: orbs  !< orbs related to the basis functions, inwhichlocreg generated in this routine
  type(atoms_data), intent(in) :: at
  real(gp), dimension(3,at%nat), intent(in) :: rxyz
  character(len=*), intent(in) :: filename
  type(local_zone_descriptors), intent(inout) :: Lzd !< must already contain Glr and hgrids
  integer, dimension(orbs%norb), optional :: orblist
  !Local variables
  character(len=*), parameter :: subname='initialize_linear_from_file'
  character(len =256) :: error
  logical :: lstat, consistent, perx, pery, perz
  integer :: ilr, ierr, iorb_old, iorb, ispinor, iorb_out, n1_old, n2_old, n3_old
  integer :: nlr, i_stat, i_all,confPotOrder, confPotOrder_old, onwhichatom_tmp, iat
! integer :: jorb
  real(gp) :: hx_old, hy_old, hz_old
! real(gp) :: mindist
  real(kind=8) :: eval
!  real(kind=8) :: dx,dy,dz,dist,eval
  real(gp), dimension(orbs%norb):: locrad, confPotprefac
  real(gp), dimension(3,at%nat) :: rxyz_old
  real(gp), dimension(3,orbs%norb) :: locregCenter
  integer, dimension(orbs%norb) :: nvctr_c, nvctr_f
  integer, dimension(:), allocatable :: norbsperlocreg
  real(kind=8), dimension(:), allocatable :: lrad
  real(gp), dimension(:,:), allocatable :: cxyz
  logical, dimension(:), allocatable :: calcbounds

  ! NOTES:
  ! The orbs%norb family must be all constructed before this routine
  ! This can be done from the input.lin since the number of basis functions should be fixed.

  call to_zero(3*orbs%norb,locregCenter(1,1))
  call to_zero(orbs%norb,locrad(1))
  call to_zero(orbs%norb,confPotprefac(1))
  consistent = .true.

  ! First read the headers (reading is distributed) and then the information is communicated to all procs.
  ! Then each proc generates a group of lrs that are communicated to all others.
  if (iformat == WF_FORMAT_ETSF) then
     stop 'Linear scaling with ETSF writing not implemented yet'
  else if (iformat == WF_FORMAT_BINARY .or. iformat == WF_FORMAT_PLAIN) then
     orbs%onwhichatom = 0
     loop_iorb: do iorb=1,orbs%norbp!*orbs%nspinor
        do ispinor=1,orbs%nspinor
           if(present(orblist)) then
              call open_filename_of_iorb(99,(iformat == WF_FORMAT_BINARY),filename, &
                   & orbs,iorb,ispinor,iorb_out, orblist(iorb+orbs%isorb))
           else
              call open_filename_of_iorb(99,(iformat == WF_FORMAT_BINARY),filename, &
                   & orbs,iorb,ispinor,iorb_out)
           end if    

           call io_read_descr_linear(99,(iformat == WF_FORMAT_PLAIN), iorb_old, eval, n1_old, n2_old, n3_old, &
                & hx_old, hy_old, hz_old, lstat, error, nvctr_c(iorb+orbs%isorb), nvctr_f(iorb+orbs%isorb),&
                & rxyz_old, at%nat, locrad(iorb+orbs%isorb), locregCenter(1,iorb+orbs%isorb), confPotOrder,&
                & confPotprefac(iorb+orbs%isorb), onwhichatom_tmp)

           ! get locregcenters from new atomic positions
           orbs%onwhichatom(iorb+orbs%isorb) = onwhichatom_tmp
           locregcenter(:,iorb+orbs%isorb) = rxyz(:,onwhichatom_tmp)
           ! DEBUG: print*,iproc,iorb,iorb+orbs%isorb,iorb_old,iorb_out
           if (.not. lstat) then ; write(*,*) trim(error) ; stop; end if
           if (iorb_old /= iorb_out) stop 'initialize_linear_from_file'
           close(99)
           !TO DO: confPotOrder_old should be read from input.lin
           if(iorb==1) confPotOrder_old = confPotOrder
           confPotOrder_old = confPotOrder
        end do
     end do loop_iorb
  else
     write(0,*) "Unknown wavefunction file format from filename."
     stop
  end if

  Lzd%nlr = orbs%norb
  ! Communication of the quantities
  if (nproc > 1)  call mpiallred(orbs%onwhichatom(1),orbs%norb,MPI_SUM,bigdft_mpi%mpi_comm,ierr)
  if (nproc > 1)  call mpiallred(locregCenter(1,1),3*orbs%norb,MPI_SUM,bigdft_mpi%mpi_comm,ierr)
  if (nproc > 1)  call mpiallred(locrad(1),orbs%norb,MPI_SUM,bigdft_mpi%mpi_comm,ierr)
  if (nproc > 1)  call mpiallred(confPotprefac(1),orbs%norb,MPI_SUM,bigdft_mpi%mpi_comm,ierr)


  allocate(Lzd%Llr(Lzd%nlr),stat=i_stat)
  do ilr=1,lzd%nlr
     lzd%Llr(ilr)=default_locreg()
  end do

  allocate(cxyz(3,Lzd%nlr),stat=i_stat)
  call memocc(i_stat,cxyz,'cxyz',subname)
  allocate(calcbounds(Lzd%nlr),stat=i_stat)
  call memocc(i_stat,calcbounds,'calcbounds',subname)
  allocate(lrad(Lzd%nlr), stat=i_stat)
  call memocc(i_stat, lrad, 'lrad', subname)
  allocate(norbsperlocreg(Lzd%nlr), stat=i_stat)
  call memocc(i_stat, norbsPerLocreg, 'norbsPerLocreg', subname)
  norbsPerLocreg=1 !should be norbsPerLocreg  

  ! Put the llr in posinp order
  ilr = 0
  do iat = 1, at%nat
     do iorb=1,orbs%norb
        if(iat == orbs%onwhichatom(iorb)) then
           ilr = ilr + 1
           cxyz(1,ilr) = locregCenter(1,iorb)
           cxyz(2,ilr) = locregCenter(2,iorb)
           cxyz(3,ilr) = locregCenter(3,iorb)
           lrad(ilr) = locrad(iorb)
        end if
     end do
  end do
  

  i_all = -product(shape(orbs%inwhichlocreg))*kind(orbs%inwhichlocreg)
  deallocate(orbs%inwhichlocreg,stat=i_stat)
  call memocc(i_stat,i_all,'orbs%inwhichlocreg',subname)
  call assignToLocreg2(iproc, nproc, orbs%norb, orbs%norb_par, at%nat, Lzd%nlr, orbs%nspin, &
       norbsperlocreg, cxyz, orbs%inwhichlocreg) 

  ! Set calcbounds correctly
  calcbounds = .false.
  do iorb = 1, orbs%norbp
     ilr = orbs%inwhichlocreg(iorb+orbs%isorb)
     calcbounds(ilr) = .true.
  end do

!TO DO: CUBIC LOCREGS
  call determine_locregSphere_parallel(iproc,nproc,Lzd%nlr,cxyz,lrad,Lzd%hgrids(1),&
       Lzd%hgrids(2),Lzd%hgrids(3),at,orbs,Lzd%Glr,Lzd%Llr,calcbounds)

  i_all = -product(shape(cxyz))*kind(cxyz)
  deallocate(cxyz,stat=i_stat)
  call memocc(i_stat,i_all,'cxyz',subname)
  i_all = -product(shape(calcbounds))*kind(calcbounds)
  deallocate(calcbounds,stat=i_stat)
  call memocc(i_stat,i_all,'calcbounds',subname)
  i_all = -product(shape(lrad))*kind(lrad)
  deallocate(lrad,stat=i_stat)
  call memocc(i_stat,i_all,'lrad',subname)
  i_all = -product(shape(norbsperlocreg))*kind(norbsperlocreg)
  deallocate(norbsperlocreg,stat=i_stat)
  call memocc(i_stat,i_all,'norbsPerLocreg',subname)

END SUBROUTINE initialize_linear_from_file

subroutine check_consistency(Lzd, at, hx_old, hy_old, hz_old, n1_old, n2_old, n3_old, &
           rxyz_old,rxyz,confPotOrder,confPotOrder_old,consistent)
  use module_base
  use module_types
  implicit none
  integer, intent(in) :: confPotOrder,confPotOrder_old, n1_old, n2_old, n3_old
  type(atoms_data), intent(in) :: at
  real(gp), intent(in) :: hx_old, hy_old, hz_old
  real(gp), dimension(3,at%nat), intent(in) :: rxyz, rxyz_old
  type(local_zone_descriptors), intent(in) :: Lzd !< must already contain Glr and hgrids
  logical, intent(out) :: consistent
  ! Local variables
  logical :: perx, pery, perz
  integer :: iat
  real(gp):: tx, ty, tz, displ, mindist  

  !conditions for periodicity in the three directions
  perx=(at%geocode /= 'F')
  pery=(at%geocode == 'P')
  perz=(at%geocode /= 'F')

  tx=0.0_gp
  ty=0.0_gp
  tz=0.0_gp
  do iat=1,at%nat
     tx=tx+mindist(perx,at%alat1,rxyz(1,iat),rxyz_old(1,iat))**2
     ty=ty+mindist(pery,at%alat2,rxyz(2,iat),rxyz_old(2,iat))**2
     tz=tz+mindist(perz,at%alat3,rxyz(3,iat),rxyz_old(3,iat))**2
  enddo
  displ=sqrt(tx+ty+tz)
  consistent = .true.
  if(hx_old /= Lzd%hgrids(1) .or. hy_old /= Lzd%hgrids(2) .or. hz_old /= Lzd%hgrids(3)) then
    write(*,"(1x,A,6F14.10)") 'Stopping because hgrid_old /= hgrid',hx_old,hy_old,hz_old,&
         Lzd%hgrids(1),Lzd%hgrids(2),Lzd%hgrids(3)
    consistent = .false.
  else if (n1_old /= Lzd%Glr%d%n1  .or. n2_old /= Lzd%Glr%d%n2 .or. n3_old /= Lzd%Glr%d%n3 ) then
    write(*,"(1x,A,6I14)") 'Stopping because global cell size',&
    n1_old,Lzd%Glr%d%n1,n2_old,Lzd%Glr%d%n2,n3_old,Lzd%Glr%d%n3
    consistent = .false.
  else if(displ > 1.d-3 ) then
    write(*,*) 'Stopping because of large displacement of molecule',displ
    consistent = .false.
  else if(confpotOrder /= confPotOrder_old) then
    write(*,*) 'Stopping because of inconsistent confPotOrder',confPotOrder,confPotOrder_old 
    consistent = .false.
  end if

END SUBROUTINE check_consistency


!> Copy old support functions from phi to phi_old
subroutine copy_old_supportfunctions(orbs,lzd,phi,lzd_old,phi_old)
  use module_base
  use module_types
  implicit none
  type(orbitals_data), intent(in) :: orbs
  type(local_zone_descriptors), intent(inout) :: lzd,lzd_old
  real(wp), dimension(:), pointer :: phi,phi_old
  !Local variables
  character(len=*), parameter :: subname='copy_old_supportfunctions'
  integer :: iseg,j,ind1,iorb,i_stat,ii,iiorb,ilr
  real(kind=8) :: tt
! integer :: i_all

  ! First copy global quantities
  call nullify_locreg_descriptors(lzd_old%glr)

  lzd_old%glr%wfd%nvctr_c = lzd%glr%wfd%nvctr_c
  lzd_old%glr%wfd%nvctr_f = lzd%glr%wfd%nvctr_f
  lzd_old%glr%wfd%nseg_c  = lzd%glr%wfd%nseg_c
  lzd_old%glr%wfd%nseg_f  = lzd%glr%wfd%nseg_f

  !allocations
  call allocate_wfd(lzd_old%glr%wfd,subname)

  do iseg=1,lzd_old%glr%wfd%nseg_c+lzd_old%glr%wfd%nseg_f
     lzd_old%glr%wfd%keyglob(1,iseg)    = lzd%glr%wfd%keyglob(1,iseg) 
     lzd_old%glr%wfd%keyglob(2,iseg)    = lzd%glr%wfd%keyglob(2,iseg)
     lzd_old%glr%wfd%keygloc(1,iseg)    = lzd%glr%wfd%keygloc(1,iseg)
     lzd_old%glr%wfd%keygloc(2,iseg)    = lzd%glr%wfd%keygloc(2,iseg)
     lzd_old%glr%wfd%keyvloc(iseg)      = lzd%glr%wfd%keyvloc(iseg)
     lzd_old%glr%wfd%keyvglob(iseg)     = lzd%glr%wfd%keyvglob(iseg)
  enddo
  !!!deallocation
  !!call deallocate_wfd(lzd%glr%wfd,subname)

  !!lzd_old%glr%d%n1 = lzd%glr%d%n1
  !!lzd_old%glr%d%n2 = lzd%glr%d%n2
  !!lzd_old%glr%d%n3 = lzd%glr%d%n3
  call copy_grid_dimensions(lzd%glr%d, lzd_old%glr%d)


  lzd_old%nlr=lzd%nlr
  nullify(lzd_old%llr)
  nullify(lzd_old%doHamAppl)
  allocate(lzd_old%llr(lzd_old%nlr))
  do ilr=1,lzd_old%nlr
      call nullify_locreg_descriptors(lzd_old%llr(ilr))
  end do

  lzd_old%hgrids(1)=lzd%hgrids(1)
  lzd_old%hgrids(2)=lzd%hgrids(2)
  lzd_old%hgrids(3)=lzd%hgrids(3)
 
  !!ii=0
  !!do ilr=1,lzd_old%nlr

  !!    ! Now copy local quantities

  !!    lzd_old%llr(ilr)%wfd%nvctr_c = lzd%llr(ilr)%wfd%nvctr_c
  !!    lzd_old%llr(ilr)%wfd%nvctr_f = lzd%llr(ilr)%wfd%nvctr_f
  !!    lzd_old%llr(ilr)%wfd%nseg_c  = lzd%llr(ilr)%wfd%nseg_c
  !!    lzd_old%llr(ilr)%wfd%nseg_f  = lzd%llr(ilr)%wfd%nseg_f

  !!    !allocations
  !!    call allocate_wfd(lzd_old%llr(ilr)%wfd,subname)

  !!    do iseg=1,lzd_old%llr(ilr)%wfd%nseg_c+lzd_old%llr(ilr)%wfd%nseg_f
  !!       lzd_old%llr(ilr)%wfd%keyglob(1,iseg)    = lzd%llr(ilr)%wfd%keyglob(1,iseg) 
  !!       lzd_old%llr(ilr)%wfd%keyglob(2,iseg)    = lzd%llr(ilr)%wfd%keyglob(2,iseg)
  !!       lzd_old%llr(ilr)%wfd%keygloc(1,iseg)    = lzd%llr(ilr)%wfd%keygloc(1,iseg)
  !!       lzd_old%llr(ilr)%wfd%keygloc(2,iseg)    = lzd%llr(ilr)%wfd%keygloc(2,iseg)
  !!       lzd_old%llr(ilr)%wfd%keyvloc(iseg)      = lzd%llr(ilr)%wfd%keyvloc(iseg)
  !!       lzd_old%llr(ilr)%wfd%keyvglob(iseg)     = lzd%llr(ilr)%wfd%keyvglob(iseg)
  !!    enddo
  !!    !!!deallocation
  !!    !!call deallocate_wfd(lzd%llr(ilr)%wfd,subname)

  !!    !!lzd_old%llr(ilr)%d%n1 = lzd%llr(ilr)%d%n1
  !!    !!lzd_old%llr(ilr)%d%n2 = lzd%llr(ilr)%d%n2
  !!    !!lzd_old%llr(ilr)%d%n3 = lzd%llr(ilr)%d%n3
  !!    call copy_grid_dimensions(lzd%llr(ilr)%d, lzd_old%llr(ilr)%d)

  !!    ii = ii + lzd_old%llr(ilr)%wfd%nvctr_c + 7*lzd_old%llr(ilr)%wfd%nvctr_f

  !!end do

  ii=0
  do iorb=1,orbs%norbp
      iiorb=orbs%isorb+iorb
      ilr=orbs%inwhichlocreg(iiorb)
      call copy_locreg_descriptors(lzd%llr(ilr), lzd_old%llr(ilr), subname)
      ii = ii + lzd_old%llr(ilr)%wfd%nvctr_c + 7*lzd_old%llr(ilr)%wfd%nvctr_f
  end do
  allocate(phi_old(ii+ndebug),stat=i_stat)
  call memocc(i_stat,phi_old,'phi_old',subname)

  ! Now copy the suport functions
  ind1=0
  do iorb=1,orbs%norbp
      tt=0.d0
      iiorb=orbs%isorb+iorb
      ilr=orbs%inwhichlocreg(iiorb)
      do j=1,lzd_old%llr(ilr)%wfd%nvctr_c+7*lzd_old%llr(ilr)%wfd%nvctr_f
          ind1=ind1+1
          phi_old(ind1)=phi(ind1)
          tt=tt+real(phi(ind1),kind=8)**2
      end do
      tt=sqrt(tt)
      if (abs(tt-1.d0) > 1.d-8) then
         write(*,*)'wrong phi_old',iiorb,tt
         stop 
      end if
  end do

  !!!deallocation
  !!i_all=-product(shape(phi))*kind(phi)
  !!deallocate(phi,stat=i_stat)
  !!call memocc(i_stat,i_all,'phi',subname)

END SUBROUTINE copy_old_supportfunctions


subroutine copy_old_coefficients(norb_KS, norb_tmb, coeff, coeff_old)
  use module_base
  implicit none

  ! Calling arguments
  integer,intent(in):: norb_KS, norb_tmb
  real(8),dimension(:,:),pointer:: coeff, coeff_old

  ! Local variables
  character(len=*),parameter:: subname='copy_old_coefficients'
  integer:: istat
!  integer:: iall

  allocate(coeff_old(norb_tmb,norb_KS),stat=istat)
  call memocc(istat,coeff_old,'coeff_old',subname)

  call vcopy(norb_KS*norb_tmb, coeff(1,1), 1, coeff_old(1,1), 1)

  !!iall=-product(shape(coeff))*kind(coeff)
  !!deallocate(coeff,stat=istat)
  !!call memocc(istat,iall,'coeff',subname)

END SUBROUTINE copy_old_coefficients


subroutine copy_old_inwhichlocreg(norb_tmb, inwhichlocreg, inwhichlocreg_old, onwhichatom, onwhichatom_old)
  use module_base
  implicit none

  ! Calling arguments
  integer,intent(in):: norb_tmb
  integer,dimension(:),pointer:: inwhichlocreg, inwhichlocreg_old, onwhichatom, onwhichatom_old

  ! Local variables
  character(len=*),parameter:: subname='copy_old_inwhichlocreg'
  integer :: istat
!  integer:: iall

  allocate(inwhichlocreg_old(norb_tmb),stat=istat)
  call memocc(istat,inwhichlocreg_old,'inwhichlocreg_old',subname)
  call vcopy(norb_tmb, inwhichlocreg(1), 1, inwhichlocreg_old(1), 1)
  !!iall=-product(shape(inwhichlocreg))*kind(inwhichlocreg)
  !!deallocate(inwhichlocreg,stat=istat)
  !!call memocc(istat,iall,'inwhichlocreg',subname)


  allocate(onwhichatom_old(norb_tmb),stat=istat)
  call memocc(istat,onwhichatom_old,'onwhichatom_old',subname)
  call vcopy(norb_tmb, onwhichatom(1), 1, onwhichatom_old(1), 1)
  !!iall=-product(shape(onwhichatom))*kind(onwhichatom)
  !!deallocate(onwhichatom,stat=istat)
  !!call memocc(istat,iall,'onwhichatom',subname)

END SUBROUTINE copy_old_inwhichlocreg


!> Reformat wavefunctions if the mesh have changed (in a restart)
subroutine reformat_supportfunctions(iproc,orbs,at,lzd_old,&
           rxyz_old,ndim_old,phi_old,lzd,rxyz,ndim,phi)
  use module_base
  use module_types
  implicit none
  integer, intent(in) :: iproc,ndim_old,ndim
  type(orbitals_data), intent(in) :: orbs
  type(local_zone_descriptors), intent(in) :: lzd_old,lzd
  type(atoms_data), intent(in) :: at
  real(gp), dimension(3,at%nat), intent(in) :: rxyz,rxyz_old
  real(wp), dimension(ndim_old), intent(in) :: phi_old
  real(wp), dimension(ndim), intent(out) :: phi
  !Local variables
  character(len=*), parameter :: subname='reformatmywaves'
  logical :: reformat,perx,pery,perz
  integer :: iorb,j,i_stat,i_all,jj,j0,j1,ii,i0,i1,i2,i3,i,iseg,nb1,nb2,nb3,jstart,jstart_old,iiorb,ilr,iiat
  integer:: n1_old,n2_old,n3_old,n1,n2,n3,ierr,idir,jstart_old_der,ncount
  real(gp) :: tx,ty,tz,displ,mindist,tt
  real(wp), dimension(:,:,:), allocatable :: phifscf
  real(wp), dimension(:,:,:,:,:,:), allocatable :: phigold
  real(wp), dimension(:), allocatable :: phi_old_der
  integer, dimension(0:5) :: reformat_reason
!  real(gp) :: dnrm2
!  integer :: iat

  reformat_reason=0


  !!do ilr=1,lzd%nlr
  !!    write(*,*) 'iproc, assoc(new)',iproc, associated(lzd%llr(ilr)%wfd%keyvloc)
  !!    write(*,*) 'iproc, assoc(old)',iproc, associated(lzd_old%llr(ilr)%wfd%keyvloc)
  !!end do

  !!do i_stat=1,ndim_old
  !!    write(800+iproc,*) i_stat,phi_old(i_stat)
  !!end do

  !conditions for periodicity in the three directions
  perx=(at%geocode /= 'F')
  pery=(at%geocode == 'P')
  perz=(at%geocode /= 'F')

  !buffers realted to periodicity
  !WARNING: the boundary conditions are not assumed to change between new and old
  call ext_buffers_coarse(perx,nb1)
  call ext_buffers_coarse(pery,nb2)
  call ext_buffers_coarse(perz,nb3)

  ! Calculate the average shift
  !!tx=0.0_gp 
  !!ty=0.0_gp
  !!tz=0.0_gp
  !!do iat=1,at%nat
  !!   tx=tx+mindist(perx,at%alat1,rxyz(1,iat),rxyz_old(1,iat))**2
  !!   ty=ty+mindist(pery,at%alat2,rxyz(2,iat),rxyz_old(2,iat))**2
  !!   tz=tz+mindist(perz,at%alat3,rxyz(3,iat),rxyz_old(3,iat))**2
  !!   !!if (iproc==0) write(333,'(i6,3es15.6)') iat, tx, ty, tz
  !!enddo
  !!!!if (iproc==0) write(333,*) '========================================'
  !!displ=sqrt(tx+ty+tz)/sqrt(dble(at%nat))
  !!if (iproc==0) write(*,*) 'mean shift of the atoms',displ

  !!if(displ<1.d-2) then
  !!    restart_method=LINEAR_HIGHACCURACY
  !!    if(iproc==0) write(*,'(1x,a)') 'Method after restart: high accuracy'
  !!else
  !!    restart_method=LINEAR_LOWACCURACY
  !!    if(iproc==0) write(*,'(1x,a)') 'Method after restart: low accuracy'
  !!end if

  allocate(phi_old_der(3*ndim_old),stat=i_stat)
  call memocc(i_stat,phi_old_der,'phi_old_der',subname)

  ! Get the derivatives of the support functions
  call get_derivative_supportfunctions(ndim_old, lzd_old%hgrids(1), lzd_old, orbs, phi_old, phi_old_der)

  jstart_old=1
  jstart_old_der=1
  jstart=1
  do iorb=1,orbs%norbp
      iiorb=orbs%isorb+iorb
      ilr=orbs%inwhichlocreg(iiorb)
      iiat=orbs%onwhichatom(iiorb)

      tx=mindist(perx,at%alat1,rxyz(1,iiat),rxyz_old(1,iiat))**2
      ty=mindist(pery,at%alat2,rxyz(2,iiat),rxyz_old(2,iiat))**2
      tz=mindist(perz,at%alat3,rxyz(3,iiat),rxyz_old(3,iiat))**2

      displ=sqrt(tx+ty+tz)

      n1_old=lzd_old%Glr%d%n1
      n2_old=lzd_old%Glr%d%n2
      n3_old=lzd_old%Glr%d%n3
      n1=lzd%Glr%d%n1
      n2=lzd%Glr%d%n2
      n3=lzd%Glr%d%n3

      !reformatting criterion
      if (lzd%hgrids(1) == lzd_old%hgrids(1) .and. lzd%hgrids(2) == lzd_old%hgrids(2) &
            .and. lzd%hgrids(3) == lzd_old%hgrids(3) .and. &
            lzd_old%llr(ilr)%wfd%nvctr_c  == lzd%llr(ilr)%wfd%nvctr_c .and. &
            lzd_old%llr(ilr)%wfd%nvctr_f == lzd%llr(ilr)%wfd%nvctr_f .and.&
            n1_old  == n1  .and. n2_old == n2 .and. n3_old == n3  .and.  displ <  1.d-3  ) then
          reformat_reason(0) = reformat_reason(0) + 1
          reformat=.false.
          !if (iproc==0) then
          !   write(*,'(1x,a)',advance='NO')&
          !    'The wavefunctions do not need reformatting and can be imported directly...   '
          !  !  '-------------------------------------------------------------- Wavefunctions Restart'
          !end if
      else
          reformat=.true.
          !if (iproc==0) then
              !write(*,'(1x,a)')&
              ! 'The wavefunctions need reformatting because:                                 '
              if (lzd%hgrids(1) /= lzd_old%hgrids(1) .or. lzd%hgrids(2) /= lzd_old%hgrids(2) &
                  .or. lzd%hgrids(3) /= lzd_old%hgrids(3)) then 
                 reformat_reason(1) = reformat_reason(1) + 1
                 !!write(*,"(4x,a,6(1pe20.12))") &
                 !!     '  hgrid_old /= hgrid  ',lzd_old%hgrids(1),lzd_old%hgrids(2),lzd_old%hgrids(3),&
                 !!     lzd%hgrids(1),lzd%hgrids(2),lzd%hgrids(3)
              end if
              if (lzd_old%llr(ilr)%wfd%nvctr_c /= lzd%llr(ilr)%wfd%nvctr_c) then
                 reformat_reason(2) = reformat_reason(2) + 1
                 !!write(*,"(4x,a,2i8)") &
                 !!     'nvctr_c_old /= nvctr_c',lzd_old%llr(ilr)%wfd%nvctr_c,lzd%llr(ilr)%wfd%nvctr_c
              end if
              if (lzd_old%llr(ilr)%wfd%nvctr_f /= lzd%llr(ilr)%wfd%nvctr_f)  then
                 reformat_reason(3) = reformat_reason(3) + 1
                 !!write(*,"(4x,a,2i8)") &
                 !!     'nvctr_f_old /= nvctr_f',lzd_old%llr(ilr)%wfd%nvctr_f,lzd%llr(ilr)%wfd%nvctr_f
              end if
              if (n1_old /= n1  .or. n2_old /= n2 .or. n3_old /= n3 )  then  
                 !!reformat_reason(4) = reformat_reason(4) + 1
                 !!write(*,"(4x,a,6i5)") &
                 !!     'cell size has changed ',n1_old,n1  , n2_old,n2 , n3_old,n3
              end if
              if (displ >=  1.d-3) then
                 reformat_reason(5) = reformat_reason(5) + 1
                 !!write(*,"(4x,a,3(1pe19.12))") &
                 !!     'molecule was shifted  ' , tx,ty,tz
              endif
                 !!write(*,"(1x,a)",advance='NO')& 
                 !!     'Reformatting...'
          !end if
         !calculate the new grid values
         
    !check
    !        write(100+iproc,'(1x,a)')&
    !         'The wavefunctions need reformatting because:                                 '
    !        if (hgrid_old.ne.hgrid) then 
    !           write(100+iproc,"(4x,a,1pe20.12)") &
    !                '  hgrid_old /= hgrid  ',hgrid_old, hgrid
    !        else if (wfd_old%nvctr_c.ne.wfd%nvctr_c) then
    !           write(100+iproc,"(4x,a,2i8)") &
    !                'nvctr_c_old /= nvctr_c',wfd_old%nvctr_c,wfd%nvctr_c
    !        else if (wfd_old%nvctr_f.ne.wfd%nvctr_f)  then
    !           write(100+iproc,"(4x,a,2i8)") &
    !                'nvctr_f_old /= nvctr_f',wfd_old%nvctr_f,wfd%nvctr_f
    !        else if (n1_old.ne.n1  .or. n2_old.ne.n2 .or. n3_old.ne.n3 )  then  
    !           write(100+iproc,"(4x,a,6i5)") &
    !                'cell size has changed ',n1_old,n1  , n2_old,n2 , n3_old,n3
    !        else
    !           write(100+iproc,"(4x,a,3(1pe19.12))") &
    !                'molecule was shifted  ' , tx,ty,tz
    !        endif
    !checkend
      end if
   
   
      if (.not. reformat) then
          !write(100+iproc,*) 'no reformatting' 
   
          do j=1,lzd_old%llr(ilr)%wfd%nvctr_c
              phi(jstart)=phi_old(jstart_old)
              jstart=jstart+1
              jstart_old=jstart_old+1
          end do
          do j=1,7*lzd_old%llr(ilr)%wfd%nvctr_f-6,7
              phi(jstart+0)=phi_old(jstart_old+0)
              phi(jstart+1)=phi_old(jstart_old+1)
              phi(jstart+2)=phi_old(jstart_old+2)
              phi(jstart+3)=phi_old(jstart_old+3)
              phi(jstart+4)=phi_old(jstart_old+4)
              phi(jstart+5)=phi_old(jstart_old+5)
              phi(jstart+6)=phi_old(jstart_old+6)
              jstart=jstart+7
              jstart_old=jstart_old+7
          end do
   
      else
   
          allocate(phifscf(-nb1:2*n1+1+nb1,-nb2:2*n2+1+nb2,-nb3:2*n3+1+nb3+ndebug),stat=i_stat)
          call memocc(i_stat,phifscf,'phifscf',subname)

          allocate(phigold(0:n1_old,2,0:n2_old,2,0:n3_old,2+ndebug),stat=i_stat)
          call memocc(i_stat,phigold,'phigold',subname)
   
          call razero(8*(n1_old+1)*(n2_old+1)*(n3_old+1),phigold(0,1,0,1,0,1))

          ! Add the derivatives to the basis functions
          do idir=1,3
              tt=rxyz(idir,iiat)-rxyz_old(idir,iiat)
              ncount = lzd_old%llr(ilr)%wfd%nvctr_c+7*lzd_old%llr(ilr)%wfd%nvctr_f
              call daxpy(ncount, tt, phi_old_der(jstart_old_der), 1, phi_old(jstart_old), 1)
              jstart_old_der = jstart_old_der + ncount
          end do
 
          ! coarse part
          do iseg=1,lzd_old%llr(ilr)%wfd%nseg_c
             jj=lzd_old%llr(ilr)%wfd%keyvglob(iseg)
             j0=lzd_old%llr(ilr)%wfd%keyglob(1,iseg)
             j1=lzd_old%llr(ilr)%wfd%keyglob(2,iseg)
             ii=j0-1
             i3=ii/((n1_old+1)*(n2_old+1))
             ii=ii-i3*(n1_old+1)*(n2_old+1)
             i2=ii/(n1_old+1)
             i0=ii-i2*(n1_old+1)
             i1=i0+j1-j0
             do i=i0,i1
                phigold(i,1,i2,1,i3,1) = phi_old(jstart_old)
                jstart_old=jstart_old+1
             end do
          end do
   
          ! fine part
          do iseg=1,lzd_old%llr(ilr)%wfd%nseg_f
             jj=lzd_old%llr(ilr)%wfd%keyvglob(lzd_old%llr(ilr)%wfd%nseg_c + iseg)
             j0=lzd_old%llr(ilr)%wfd%keyglob(1,lzd_old%llr(ilr)%wfd%nseg_c + iseg)
             j1=lzd_old%llr(ilr)%wfd%keyglob(2,lzd_old%llr(ilr)%wfd%nseg_c + iseg)
             ii=j0-1
             i3=ii/((n1_old+1)*(n2_old+1))
             ii=ii-i3*(n1_old+1)*(n2_old+1)
             i2=ii/(n1_old+1)
             i0=ii-i2*(n1_old+1)
             i1=i0+j1-j0
             do i=i0,i1
                   phigold(i,2,i2,1,i3,1)=phi_old(jstart_old+0)
                   phigold(i,1,i2,2,i3,1)=phi_old(jstart_old+1)
                   phigold(i,2,i2,2,i3,1)=phi_old(jstart_old+2)
                   phigold(i,1,i2,1,i3,2)=phi_old(jstart_old+3)
                   phigold(i,2,i2,1,i3,2)=phi_old(jstart_old+4)
                   phigold(i,1,i2,2,i3,2)=phi_old(jstart_old+5)
                   phigold(i,2,i2,2,i3,2)=phi_old(jstart_old+6)
                jstart_old=jstart_old+7
             end do
          end do
   
          !write(100+iproc,*) 'norm phigold ',dnrm2(8*(n1_old+1)*(n2_old+1)*(n3_old+1),phigold,1)
          !write(*,*) 'iproc,norm phigold ',iproc,dnrm2(8*(n1_old+1)*(n2_old+1)*(n3_old+1),phigold,1)
   
          !!call reformatonewave(displ,lzd%llr(ilr)%wfd,at,lzd_old%hgrids(1),lzd_old%hgrids(2),lzd_old%hgrids(3), & !n(m)
          !!     n1_old,n2_old,n3_old,rxyz_old,phigold,lzd%hgrids(1),lzd%hgrids(2),lzd%hgrids(3),&
          !!     n1,n2,n3,rxyz,phifscf,phi(jstart))
          call reformat_one_supportfunction(iiat,displ,lzd%llr(ilr)%wfd,at,lzd_old%hgrids(1),lzd_old%hgrids(2),lzd_old%hgrids(3), & !n(m)
               n1_old,n2_old,n3_old,rxyz_old,phigold,lzd%hgrids(1),lzd%hgrids(2),lzd%hgrids(3),&
               n1,n2,n3,rxyz,phifscf,phi(jstart))

          jstart=jstart+lzd%llr(ilr)%wfd%nvctr_c+7*lzd%llr(ilr)%wfd%nvctr_f
   
          i_all=-product(shape(phifscf))*kind(phifscf)
          deallocate(phifscf,stat=i_stat)
          call memocc(i_stat,i_all,'phifscf',subname)
   
          i_all=-product(shape(phigold))*kind(phigold)
          deallocate(phigold,stat=i_stat)
          call memocc(i_stat,i_all,'phigold',subname)

      end if

      if (iproc==0) write(*,"(1x,a)")'done.'

  end do

  i_all=-product(shape(phi_old_der))*kind(phi_old_der)
  deallocate(phi_old_der,stat=i_stat)
  call memocc(i_stat,i_all,'phi_old_der',subname)


  call mpiallred(reformat_reason(0), 6, mpi_sum, mpi_comm_world, ierr)
  if (iproc==0) then
        write(*,'(1x,a)') 'Overview of the reformatting (several categories may apply):'
        write(*,'(3x,a,i0)') '- No reformating required: ', reformat_reason(0)
        write(*,'(3x,a,i0)') '- Grid spacing has changed: ', reformat_reason(1)
        write(*,'(3x,a,i0)') '- number of coarse grid points has changed: ', reformat_reason(2)
        write(*,'(3x,a,i0)') '- number of fine grid points has changed: ', reformat_reason(3)
        write(*,'(3x,a,i0)') '- box size has changed: ', reformat_reason(4)
        write(*,'(3x,a,i0)') '- molecule was shifted: ', reformat_reason(5)
  end if

END SUBROUTINE reformat_supportfunctions




