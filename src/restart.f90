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
  integer :: iseg,nvctrp_old,j,ind1,iorb,i_all,i_stat,oidx,sidx
  real(kind=8) :: tt

  wfd_old%nvctr_c = wfd%nvctr_c
  wfd_old%nvctr_f = wfd%nvctr_f
  wfd_old%nseg_c  = wfd%nseg_c
  wfd_old%nseg_f  = wfd%nseg_f

  !allocations
  call allocate_wfd(wfd_old,subname)

  do iseg=1,wfd_old%nseg_c+wfd_old%nseg_f
     wfd_old%keyg(1,iseg)    = wfd%keyg(1,iseg)
     wfd_old%keyg(2,iseg)    = wfd%keyg(2,iseg)
     wfd_old%keyv(iseg)      = wfd%keyv(iseg)
  enddo
  !deallocation
  call deallocate_wfd(wfd,subname)

  n1_old = n1
  n2_old = n2
  n3_old = n3

  !add the number of distributed point for the compressed wavefunction
  tt=dble(wfd_old%nvctr_c+7*wfd_old%nvctr_f)/dble(nproc)
  nvctrp_old=int((1.d0-eps_mach*tt) + tt)

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
           jj=wfd_old%keyv(iseg)
           j0=wfd_old%keyg(1,iseg)
           j1=wfd_old%keyg(2,iseg)
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
           jj=wfd_old%keyv(wfd_old%nseg_c + iseg)
           j0=wfd_old%keyg(1,wfd_old%nseg_c + iseg)
           j1=wfd_old%keyg(2,wfd_old%nseg_c + iseg)
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

        call reformatonewave(iproc,displ,wfd,at,hx_old,hy_old,hz_old, &
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


!>  Reads wavefunction from file and transforms it properly if hgrid or size of simulation cell
!!  have changed
subroutine readmywaves(iproc,filename,orbs,n1,n2,n3,hx,hy,hz,at,rxyz_old,rxyz,  & 
     wfd,psi)
  use module_base
  use module_types
  implicit none
  integer, intent(in) :: iproc,n1,n2,n3
  real(gp), intent(in) :: hx,hy,hz
  type(wavefunctions_descriptors), intent(in) :: wfd
  type(orbitals_data), intent(inout) :: orbs
  type(atoms_data), intent(in) :: at
  real(gp), dimension(3,at%nat), intent(in) :: rxyz
  real(gp), dimension(3,at%nat), intent(out) :: rxyz_old
  real(wp), dimension(wfd%nvctr_c+7*wfd%nvctr_f,orbs%nspinor,orbs%norbp), intent(out) :: psi
  character(len=*), intent(in) :: filename
  !Local variables
  character(len=*), parameter :: subname='readmywaves'
  logical :: perx,pery,perz,exists
  integer :: ncount1,ncount_rate,ncount_max,iorb,i_stat,i_all,ncount2,nb1,nb2,nb3,iorb_out,ispinor
  real(kind=4) :: tr0,tr1
  real(kind=8) :: tel
  real(wp), dimension(:,:,:), allocatable :: psifscf

  inquire(file=trim(filename)//".etsf",exist=exists)
  if (exists) then
     if (iproc ==0) write(*,*) "Reading wavefunctions in ETSF file format."
     call read_waves_etsf(iproc,filename//".etsf",orbs,n1,n2,n3,hx,hy,hz,at,rxyz_old,rxyz,  & 
     wfd,psi)
  else
     call cpu_time(tr0)
     call system_clock(ncount1,ncount_rate,ncount_max)

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

     inquire(file=trim(filename)//".bin.0001",exist=exists)
     if (exists) then
        if (iproc ==0) write(*,*) "Reading wavefunctions in BigDFT binary file format."
     else
        if (iproc ==0) write(*,*) "Reading wavefunctions in plain text file format."
     end if

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
           call open_filename_of_iorb(99,exists,filename,orbs,iorb,ispinor,iorb_out)           
           call readonewave(99, .not.exists,iorb_out,iproc,n1,n2,n3, &
                & hx,hy,hz,at,wfd,rxyz_old,rxyz,&
                psi(1,ispinor,iorb),orbs%eval(orbs%isorb+iorb),psifscf)
           close(99)
        end do

     end do

     i_all=-product(shape(psifscf))*kind(psifscf)
     deallocate(psifscf,stat=i_stat)
     call memocc(i_stat,i_all,'psifscf',subname)

     call cpu_time(tr1)
     call system_clock(ncount2,ncount_rate,ncount_max)
     tel=dble(ncount2-ncount1)/dble(ncount_rate)
     write(*,'(a,i4,2(1x,1pe10.3))') '- READING WAVES TIME',iproc,tr1-tr0,tel
  end if
END SUBROUTINE readmywaves

subroutine verify_file_presence(orbs,allfiles)
  use module_base
  use module_types
  implicit none
  type(orbitals_data), intent(in) :: orbs
  logical, intent(out) :: allfiles
  !local variables
  character(len=50) :: filename
  logical :: onefile
  integer :: iorb,ispinor,iorb_out,ierr
  
  allfiles=.true.

  !first try with plain files
  loop_plain: do iorb=1,orbs%norbp
     do ispinor=1,orbs%nspinor
        call filename_of_iorb(.false.,"wavefunction",orbs,iorb,ispinor,filename,iorb_out)
        inquire(file=filename,exist=onefile)
        allfiles=allfiles .and. onefile
        if (.not. allfiles) then
           exit loop_plain
        end if
     end do
  end do loop_plain
  !reduce the result among the other processors
  call mpiallred(allfiles,1,MPI_LAND,MPI_COMM_WORLD,ierr)
 
  !Otherwise  test binary files.
  if (.not. allfiles) then           
     allfiles = .true.
     loop_binary: do iorb=1,orbs%norbp
        do ispinor=1,orbs%nspinor
           call filename_of_iorb(.true.,"wavefunction",orbs,iorb,ispinor,filename,iorb_out)

           inquire(file=filename,exist=onefile)
           allfiles=allfiles .and. onefile
           if (.not. allfiles) then
              exit loop_binary
           end if

        end do
     end do loop_binary
     !reduce the result among the other processors
     call mpiallred(allfiles,1,MPI_LAND,MPI_COMM_WORLD,ierr)
  end if

end subroutine verify_file_presence

!> Associate to the absolute value of orbital a filename which depends of the k-point and 
!! of the spin sign
subroutine filename_of_iorb(lbin,filename,orbs,iorb,ispinor,filename_out,iorb_out)
  use module_base
  use module_types
  implicit none
  character(len=*), intent(in) :: filename
  logical, intent(in) :: lbin
  integer, intent(in) :: iorb,ispinor
  type(orbitals_data), intent(in) :: orbs
  character(len=*) :: filename_out
  integer, intent(out) :: iorb_out
  !local variables
  character(len=1) :: spintype,realimag
  character(len=3) :: f3
  character(len=4) :: f4
  character(len=7) :: completename
  integer :: ikpt
  real(gp) :: spins

  !calculate k-point
  ikpt=orbs%iokpt(iorb)
  write(f3,'(i3.3)') ikpt !not more than 999 kpts
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
  !purge the value from the spin sign
  if (spins==-1.0_gp) iorb_out=iorb_out-orbs%norbu

  !value of the orbital 
  write(f4,'(i4.4)') iorb_out

  !complete the information in the name of the orbital
  completename='-'//f3//'-'//spintype//realimag

  if (lbin) then
     filename_out = trim(filename)//completename//".bin."//f4
  else
     filename_out = trim(filename)//completename//"."//f4
  end if

end subroutine filename_of_iorb


!> Associate to the absolute value of orbital a filename which depends of the k-point and 
!! of the spin sign
subroutine open_filename_of_iorb(unitfile,lbin,filename,orbs,iorb,ispinor,iorb_out)
  use module_base
  use module_types
  implicit none
  character(len=*), intent(in) :: filename
  logical, intent(in) :: lbin
  integer, intent(in) :: iorb,ispinor,unitfile
  type(orbitals_data), intent(in) :: orbs
  integer, intent(out) :: iorb_out
  !local variables
  character(len=50) :: filename_out

  call filename_of_iorb(lbin,filename,orbs,iorb,ispinor,filename_out,iorb_out)

  if (lbin) then
     open(unit=unitfile,file=trim(filename_out),status='unknown',form="unformatted")
  else
     open(unit=unitfile,file=trim(filename_out),status='unknown')
  end if

end subroutine open_filename_of_iorb

!>   Write all my wavefunctions in files by calling writeonewave
subroutine writemywaves(iproc,filename,orbs,n1,n2,n3,hx,hy,hz,at,rxyz,wfd,psi)
  use module_types
  use module_base
  implicit none
  integer, intent(in) :: iproc,n1,n2,n3
  real(gp), intent(in) :: hx,hy,hz
  type(atoms_data), intent(in) :: at
  type(orbitals_data), intent(in) :: orbs
  type(wavefunctions_descriptors), intent(in) :: wfd
  real(gp), dimension(3,at%nat), intent(in) :: rxyz
  real(wp), dimension(wfd%nvctr_c+7*wfd%nvctr_f,orbs%nspinor,orbs%norbp), intent(in) :: psi
  character(len=*), intent(in) :: filename
  !Local variables
  integer :: ncount1,ncount_rate,ncount_max,iorb,ncount2,isuffix,iorb_out,ispinor
  real(kind=4) :: tr0,tr1
  real(kind=8) :: tel

  if (iproc == 0) write(*,"(1x,A,A)") "Write wavefunctions to file: ", trim(filename)
  isuffix = index(filename, ".etsf", back = .true.)
  if (isuffix <= 0) isuffix = index(filename, ".etsf.nc", back = .true.)
  if (isuffix > 0) then
     call write_waves_etsf(iproc,filename,orbs,n1,n2,n3,hx,hy,hz,at,rxyz,wfd,psi)
  else
     call cpu_time(tr0)
     call system_clock(ncount1,ncount_rate,ncount_max)

     isuffix = index(filename, ".bin", back = .true.)

     ! Plain BigDFT files.
     do iorb=1,orbs%norbp!*orbs%nspinor

!!$        write(f4,'(i4.4)')  iorb+orbs%isorb*orbs%nspinor
!!$        filename_ = trim(filename)//"."//f4
!!$        if (verbose > 2) write(*,*) 'opening ',filename_
!!$        if (isuffix <= 0) then
!!$           open(unit=99,file=filename_,status='unknown')
!!$        else
!!$           open(unit=99,file=trim(filename_),status='unknown',form="unformatted")
!!$        end if
!!$        call writeonewave(99,(isuffix <= 0),iorb+orbs%isorb*orbs%nspinor,n1,n2,n3,hx,hy,hz,at%nat,rxyz,  & 
!!$             wfd%nseg_c,wfd%nvctr_c,wfd%keyg(1,1),wfd%keyv(1),  & 
!!$             wfd%nseg_f,wfd%nvctr_f,wfd%keyg(1,wfd%nseg_c+1),wfd%keyv(wfd%nseg_c+1), & 
!!$             psi(1,iorb),psi(wfd%nvctr_c+1,iorb), &
!!$             orbs%eval((iorb-1)/orbs%nspinor+1+orbs%isorb))
!!$        close(99)

        do ispinor=1,orbs%nspinor
           call open_filename_of_iorb(99,(isuffix > 0),filename,orbs,iorb,ispinor,iorb_out)           
           call writeonewave(99,(isuffix <= 0),iorb_out,n1,n2,n3,hx,hy,hz,at%nat,rxyz,  & 
                wfd%nseg_c,wfd%nvctr_c,wfd%keyg(1,1),wfd%keyv(1),  & 
                wfd%nseg_f,wfd%nvctr_f,wfd%keyg(1,wfd%nseg_c+1),wfd%keyv(wfd%nseg_c+1), & 
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
