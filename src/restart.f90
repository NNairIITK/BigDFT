!> @file
!!  Routines to do restart
!! @author
!!    Copyright (C) 2007-2011 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 


!> Copy old wavefunctions from psi to psi_old
subroutine copy_old_wavefunctions(nproc,orbs,n1,n2,n3,wfd,psi,&
     n1_old,n2_old,n3_old,wfd_old,psi_old)
  use module_base
  use module_types
  use yaml_output
  implicit none
  integer, intent(in) :: nproc,n1,n2,n3
  type(orbitals_data), intent(in) :: orbs
  type(wavefunctions_descriptors), intent(inout) :: wfd,wfd_old
  integer, intent(out) :: n1_old,n2_old,n3_old
  real(wp), dimension(:), pointer :: psi,psi_old
  !Local variables
  character(len=*), parameter :: subname='copy_old_wavefunctions'
  !real(kind=8), parameter :: eps_mach=1.d-12
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
        call yaml_warning('wrong psi_old' // trim(yaml_toa(iorb)) // trim(yaml_toa(tt)))
        !write(*,*)'wrong psi_old',iorb,tt
        stop 
     end if
  enddo
  !deallocation
  i_all=-product(shape(psi))*kind(psi)
  deallocate(psi,stat=i_stat)
  call memocc(i_stat,i_all,'psi',subname)

END SUBROUTINE copy_old_wavefunctions


!> Reformat wavefunctions if the mesh have changed (in a restart)
subroutine reformatmywaves(iproc,orbs,at,&
     hx_old,hy_old,hz_old,n1_old,n2_old,n3_old,rxyz_old,wfd_old,psi_old,&
     hx,hy,hz,n1,n2,n3,rxyz,wfd,psi)
  use module_base
  use module_types
  use yaml_output
  implicit none
  integer, intent(in) :: iproc,n1_old,n2_old,n3_old,n1,n2,n3
  real(gp), intent(in) :: hx_old,hy_old,hz_old,hx,hy,hz
  type(wavefunctions_descriptors), intent(in) :: wfd,wfd_old
  type(atoms_data), intent(in) :: at
  type(orbitals_data), intent(in) :: orbs
  real(gp), dimension(3,at%astruct%nat), intent(in) :: rxyz,rxyz_old
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
  perx=(at%astruct%geocode /= 'F')
  pery=(at%astruct%geocode == 'P')
  perz=(at%astruct%geocode /= 'F')

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

  do iat=1,at%astruct%nat
     tx=tx+mindist(perx,at%astruct%cell_dim(1),rxyz(1,iat),rxyz_old(1,iat))**2
     ty=ty+mindist(pery,at%astruct%cell_dim(2),rxyz(2,iat),rxyz_old(2,iat))**2
     tz=tz+mindist(perz,at%astruct%cell_dim(3),rxyz(3,iat),rxyz_old(3,iat))**2
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
        call yaml_map('Reformating Wavefunctions',.false.)
      !write(*,'(1x,a)',advance='NO') 'The wavefunctions do not need reformatting and can be imported directly...'
      !  '-------------------------------------------------------------- Wavefunctions Restart'
     end if
  else
     reformat=.true.
     if (iproc==0) then
        call yaml_map('Reformating wavefunctions',.true.)
        call yaml_open_map('Reformatting for')
        !write(*,'(1x,a)') 'The wavefunctions need reformatting because:'
        if (hx /= hx_old .or. hy /= hy_old .or. hz /= hz_old) then 
           call yaml_open_map('hgrid modified',flow=.true.)
              call yaml_map('hgrid_old', (/ hx_old,hy_old,hz_old /),fmt='(1pe20.12)')
              call yaml_map('hgrid', (/ hx,hy,hz /), fmt='(1pe20.12)')
           call yaml_close_map()
           !write(*,"(4x,a,6(1pe20.12))") 'hgrid_old /= hgrid  ',hx_old,hy_old,hz_old,hx,hy,hz
        else if (wfd_old%nvctr_c /= wfd%nvctr_c) then
           call yaml_map('nvctr_c modified', (/ wfd_old%nvctr_c,wfd%nvctr_c /))
           !write(*,"(4x,a,2i8)") 'nvctr_c_old /= nvctr_c',wfd_old%nvctr_c,wfd%nvctr_c
        else if (wfd_old%nvctr_f /= wfd%nvctr_f)  then
           call yaml_map('nvctr_f modified', (/ wfd_old%nvctr_f,wfd%nvctr_f /))
           !write(*,"(4x,a,2i8)") 'nvctr_f_old /= nvctr_f',wfd_old%nvctr_f,wfd%nvctr_f
        else if (n1_old /= n1  .or. n2_old /= n2 .or. n3_old /= n3 )  then  
           call yaml_map('Cell size has changed ', (/ n1_old,n1  , n2_old,n2 , n3_old,n3 /))
           !write(*,"(4x,a,6i5)") 'cell size has changed ',n1_old,n1  , n2_old,n2 , n3_old,n3
        else
           call yaml_map('Molecule was shifted' ,  (/ tx,ty,tz /), fmt='(1pe19.12)')
           !write(*,"(4x,a,3(1pe19.12))") 'molecule was shifted  ' , tx,ty,tz
        endif
        !write(*,"(1x,a)",advance='NO') 'Reformatting...'
        call yaml_close_map()
     end if
     !calculate the new grid values
     
!check
!        write(100+iproc,'(1x,a)') 'The wavefunctions need reformatting because:'
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

  !if (iproc==0) write(*,"(1x,a)")'done.'

END SUBROUTINE reformatmywaves


integer function wave_format_from_filename(iproc, filename)
  use module_types
  use yaml_output
  implicit none
  integer, intent(in) :: iproc
  character(len=*), intent(in) :: filename

  integer :: isuffix

  wave_format_from_filename = WF_FORMAT_NONE

  isuffix = index(filename, ".etsf", back = .true.)
  if (isuffix > 0) then
     wave_format_from_filename = WF_FORMAT_ETSF
     if (iproc ==0) call yaml_comment('Reading wavefunctions in ETSF file format.')
     !if (iproc ==0) write(*,*) "Reading wavefunctions in ETSF file format."
  else
     isuffix = index(filename, ".bin", back = .true.)
     if (isuffix > 0) then
        wave_format_from_filename = WF_FORMAT_BINARY
        if (iproc ==0) call yaml_comment('Reading wavefunctions in BigDFT binary file format.')
        !if (iproc ==0) write(*,*) "Reading wavefunctions in BigDFT binary file format."
     else
        wave_format_from_filename = WF_FORMAT_PLAIN
        if (iproc ==0) call yaml_comment('Reading wavefunctions in plain text file format.')
        !if (iproc ==0) write(*,*) "Reading wavefunctions in plain text file format."
     end if
  end if
end function wave_format_from_filename


!> Reads wavefunction from file and transforms it properly if hgrid or size of simulation cell
!!  have changed
subroutine readmywaves(iproc,filename,iformat,orbs,n1,n2,n3,hx,hy,hz,at,rxyz_old,rxyz,  & 
     wfd,psi,orblist)
  use module_base
  use module_types
  use yaml_output
  use module_interfaces, except_this_one => readmywaves
  implicit none
  integer, intent(in) :: iproc,n1,n2,n3, iformat
  real(gp), intent(in) :: hx,hy,hz
  type(wavefunctions_descriptors), intent(in) :: wfd
  type(orbitals_data), intent(inout) :: orbs
  type(atoms_data), intent(in) :: at
  real(gp), dimension(3,at%astruct%nat), intent(in) :: rxyz
  real(gp), dimension(3,at%astruct%nat), intent(out) :: rxyz_old
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
     perx=(at%astruct%geocode /= 'F')
     pery=(at%astruct%geocode == 'P')
     perz=(at%astruct%geocode /= 'F')

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
     call yaml_warning('Unknown wavefunction file format from filename.')
     stop
  end if

  call cpu_time(tr1)
  call system_clock(ncount2,ncount_rate,ncount_max)
  tel=dble(ncount2-ncount1)/dble(ncount_rate)

  if (iproc == 0) then 
     call yaml_open_sequence('Reading Waves Time')
     call yaml_sequence(advance='no')
     call yaml_open_map(flow=.true.)
     call yaml_map('Process',iproc)
     call yaml_map('Timing',(/ real(tr1-tr0,kind=8),tel /),fmt='(1pe10.3)')
     call yaml_close_map()
     call yaml_close_sequence()
  end if
  !write(*,'(a,i4,2(1x,1pe10.3))') '- READING WAVES TIME',iproc,tr1-tr0,tel
END SUBROUTINE readmywaves


!> Verify the presence of a given file
subroutine verify_file_presence(filerad,orbs,iformat,nproc,nforb)
  use module_base
  use module_types
  use module_interfaces, except_this_one=>verify_file_presence
  implicit none
  integer, intent(in) :: nproc
  character(len=*), intent(in) :: filerad
  type(orbitals_data), intent(in) :: orbs
  integer, intent(out) :: iformat
  integer, optional, intent(in) :: nforb
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
        ! for fragment calculations, this number could be less than the number of orbs - find a better way of doing this
        if (.not. present(nforb)) then
           allfiles=allfiles .and. onefile
        else if (iorb_out<=nforb) then
           allfiles=allfiles .and. onefile
        end if
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
        ! for fragment calculations, this number could be less than the number of orbs - find a better way of doing this
        if (.not. present(nforb)) then
           allfiles=allfiles .and. onefile
        else if (iorb_out<=nforb) then
           allfiles=allfiles .and. onefile
        end if
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


!> Write all my wavefunctions in files by calling writeonewave
subroutine writemywaves(iproc,filename,iformat,orbs,n1,n2,n3,hx,hy,hz,at,rxyz,wfd,psi)
  use module_types
  use module_base
  use yaml_output
  use module_interfaces, except_this_one => writeonewave
  implicit none
  integer, intent(in) :: iproc,n1,n2,n3,iformat
  real(gp), intent(in) :: hx,hy,hz
  type(atoms_data), intent(in) :: at
  type(orbitals_data), intent(in) :: orbs
  type(wavefunctions_descriptors), intent(in) :: wfd
  real(gp), dimension(3,at%astruct%nat), intent(in) :: rxyz
  real(wp), dimension(wfd%nvctr_c+7*wfd%nvctr_f,orbs%nspinor,orbs%norbp), intent(in) :: psi
  character(len=*), intent(in) :: filename
  !Local variables
  integer :: ncount1,ncount_rate,ncount_max,iorb,ncount2,iorb_out,ispinor
  real(kind=4) :: tr0,tr1
  real(kind=8) :: tel

  if (iproc == 0) call yaml_map('Write wavefunctions to file', trim(filename) // '.*')
  !if (iproc == 0) write(*,"(1x,A,A,a)") "Write wavefunctions to file: ", trim(filename),'.*'
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
                at%astruct%nat,rxyz,wfd%nseg_c,wfd%nvctr_c,wfd%keygloc(1,1),wfd%keyvloc(1),  & 
                wfd%nseg_f,wfd%nvctr_f,wfd%keygloc(1,wfd%nseg_c+1),wfd%keyvloc(wfd%nseg_c+1), & 
                psi(1,ispinor,iorb),psi(wfd%nvctr_c+1,ispinor,iorb), &
                orbs%eval(iorb+orbs%isorb))
           close(99)
        end do
     enddo

     call cpu_time(tr1)
     call system_clock(ncount2,ncount_rate,ncount_max)
     tel=dble(ncount2-ncount1)/dble(ncount_rate)
     if (iproc == 0) then
        call yaml_open_sequence('Write Waves Time')
        call yaml_sequence(advance='no')
        call yaml_open_map(flow=.true.)
        call yaml_map('Process',iproc)
        call yaml_map('Timing',(/ real(tr1-tr0,kind=8),tel /),fmt='(1pe10.3)')
        call yaml_close_map()
        call yaml_close_sequence()
     end if
     !write(*,'(a,i4,2(1x,1pe10.3))') '- WRITE WAVES TIME',iproc,tr1-tr0,tel
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


subroutine writeonewave_linear(unitwf,useFormattedOutput,iorb,n1,n2,n3,ns1,ns2,ns3,hx,hy,hz,locregCenter,&
     locrad,confPotOrder,confPotprefac,nat,rxyz, nseg_c,nvctr_c,keyg_c,keyv_c,  &
     nseg_f,nvctr_f,keyg_f,keyv_f, &
     psi_c,psi_f,eval,onwhichatom)
  use module_base
  use yaml_output
  implicit none
  logical, intent(in) :: useFormattedOutput
  integer, intent(in) :: unitwf,iorb,n1,n2,n3,ns1,ns2,ns3,nat,nseg_c,nvctr_c,nseg_f,nvctr_f,confPotOrder
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
     write(unitwf,*) ns1,ns2,ns3
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
     write(unitwf) ns1,ns2,ns3
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

  if (verbose >= 2 .and. bigdft_mpi%iproc==0) call yaml_map('Wavefunction written No.',iorb)
  !if (verbose >= 2) write(*,'(1x,i0,a)') iorb,'th wavefunction written'

END SUBROUTINE writeonewave_linear


subroutine writeLinearCoefficients(unitwf,useFormattedOutput,nat,rxyz,&
           ntmb,norb,coeff,eval)
  use module_base
  use yaml_output
  implicit none
  logical, intent(in) :: useFormattedOutput
  integer, intent(in) :: unitwf,nat,ntmb,norb
  real(wp), dimension(ntmb,ntmb), intent(in) :: coeff
  real(wp), dimension(ntmb), intent(in) :: eval
  real(gp), dimension(3,nat), intent(in) :: rxyz
  !local variables
  integer :: iat,i,j,iorb
  real(wp) :: tt

  ! Write the Header
  if (useFormattedOutput) then
     write(unitwf,*) ntmb,norb
     write(unitwf,*) nat
     do iat=1,nat
     write(unitwf,'(3(1x,e24.17))') (rxyz(j,iat),j=1,3)
     enddo
     do iorb=1,ntmb
     write(unitwf,*) iorb,eval(iorb)
     enddo
  else
     write(unitwf) ntmb, norb
     write(unitwf) nat
     do iat=1,nat
     write(unitwf) (rxyz(j,iat),j=1,3)
     enddo
     do iorb=1,ntmb
     write(unitwf) iorb,eval(iorb)
     enddo
  end if

  ! Now write the coefficients
  do i = 1, ntmb
     do j = 1, ntmb
          tt = coeff(j,i)
          if (useFormattedOutput) then
             write(unitwf,'(2(i6,1x),e19.12)') i,j,tt
          else
             write(unitwf) i,j,tt
          end if
     end do
  end do  
  if (verbose >= 2 .and. bigdft_mpi%iproc==0) call yaml_map('Wavefunction coefficients written',.true.)

END SUBROUTINE writeLinearCoefficients


!write Hamiltonian, overlap and kernel matrices in tmb basis
subroutine write_linear_matrices(iproc,nproc,filename,iformat,tmb,at,rxyz)
  use module_types
  use module_base
  use yaml_output
  use module_interfaces, except_this_one => writeonewave
  implicit none
  integer, intent(in) :: iproc,nproc,iformat
  character(len=*), intent(in) :: filename 
  type(DFT_wavefunction), intent(inout) :: tmb
  type(atoms_data), intent(inout) :: at
  real(gp),dimension(3,at%astruct%nat),intent(in) :: rxyz

  integer :: iorb, jorb, i_stat, i_all, iat, jat
  character(len=*),parameter :: subname='write_linear_matrices'

  if (iproc==0) then
     if(iformat == WF_FORMAT_PLAIN) then
        open(99, file=filename//'hamiltonian.bin', status='unknown',form='formatted')
     else
        open(99, file=filename//'hamiltonian.bin', status='unknown',form='unformatted')
     end if

     allocate(tmb%linmat%ham%matrix(tmb%linmat%ham%full_dim1,tmb%linmat%ham%full_dim1), stat=i_stat)
     call memocc(i_stat, tmb%linmat%ham%matrix, 'tmb%linmat%ham%matrix', subname)

     call uncompressMatrix(iproc,tmb%linmat%ham)

     do iorb=1,tmb%linmat%ham%full_dim1
        iat=tmb%orbs%onwhichatom(iorb)
        do jorb=1,tmb%linmat%ham%full_dim1
           jat=tmb%orbs%onwhichatom(jorb)
           if (iformat == WF_FORMAT_PLAIN) then
              write(99,'(2(i6,1x),e19.12,2(1x,i6))') iorb,jorb,tmb%linmat%ham%matrix(iorb,jorb),iat,jat
           else
              write(99) iorb,jorb,tmb%linmat%ham%matrix(iorb,jorb),iat,jat
           end if
        end do
     end do

     i_all = -product(shape(tmb%linmat%ham%matrix))*kind(tmb%linmat%ham%matrix)
     deallocate(tmb%linmat%ham%matrix,stat=i_stat)
     call memocc(i_stat,i_all,'tmb%linmat%ham%matrix',subname)

     close(99)

     if(iformat == WF_FORMAT_PLAIN) then
        open(99, file=filename//'overlap.bin', status='unknown',form='formatted')
     else
        open(99, file=filename//'overlap.bin', status='unknown',form='unformatted')
     end if

     allocate(tmb%linmat%ovrlp%matrix(tmb%linmat%ovrlp%full_dim1,tmb%linmat%ovrlp%full_dim1), stat=i_stat)
     call memocc(i_stat, tmb%linmat%ovrlp%matrix, 'tmb%linmat%ovrlp%matrix', subname)

     call uncompressMatrix(iproc,tmb%linmat%ovrlp)

     do iorb=1,tmb%linmat%ovrlp%full_dim1
        iat=tmb%orbs%onwhichatom(iorb)
        do jorb=1,tmb%linmat%ovrlp%full_dim1
           jat=tmb%orbs%onwhichatom(jorb)
           if (iformat == WF_FORMAT_PLAIN) then
              write(99,'(2(i6,1x),e19.12,2(1x,i6))') iorb,jorb,tmb%linmat%ovrlp%matrix(iorb,jorb),iat,jat
           else
              write(99) iorb,jorb,tmb%linmat%ovrlp%matrix(iorb,jorb),iat,jat
           end if
        end do
     end do

     i_all = -product(shape(tmb%linmat%ovrlp%matrix))*kind(tmb%linmat%ovrlp%matrix)
     deallocate(tmb%linmat%ovrlp%matrix,stat=i_stat)
     call memocc(i_stat,i_all,'tmb%linmat%ovrlp%matrix',subname)

     close(99)

     if(iformat == WF_FORMAT_PLAIN) then
        open(99, file=filename//'density_kernel.bin', status='unknown',form='formatted')
     else
        open(99, file=filename//'density_kernel.bin', status='unknown',form='unformatted')
     end if

     allocate(tmb%linmat%denskern%matrix(tmb%linmat%denskern%full_dim1,tmb%linmat%denskern%full_dim1), stat=i_stat)
     call memocc(i_stat, tmb%linmat%denskern%matrix, 'tmb%linmat%denskern%matrix', subname)

     call uncompressMatrix(iproc,tmb%linmat%denskern)

     do iorb=1,tmb%linmat%denskern%full_dim1
        iat=tmb%orbs%onwhichatom(iorb)
        do jorb=1,tmb%linmat%denskern%full_dim1
           jat=tmb%orbs%onwhichatom(jorb)
           if (iformat == WF_FORMAT_PLAIN) then
              write(99,'(2(i6,1x),e19.12,2(1x,i6))') iorb,jorb,tmb%linmat%denskern%matrix(iorb,jorb),iat,jat
           else
              write(99) iorb,jorb,tmb%linmat%denskern%matrix(iorb,jorb),iat,jat
           end if
        end do
     end do

     i_all = -product(shape(tmb%linmat%denskern%matrix))*kind(tmb%linmat%denskern%matrix)
     deallocate(tmb%linmat%denskern%matrix,stat=i_stat)
     call memocc(i_stat,i_all,'tmb%linmat%denskern%matrix',subname)

     close(99)

  end if

  ! calculate 'onsite' overlap matrix as well - needs double checking

  allocate(tmb%linmat%ovrlp%matrix(tmb%linmat%ovrlp%full_dim1,tmb%linmat%ovrlp%full_dim1), stat=i_stat)
  call memocc(i_stat, tmb%linmat%ovrlp%matrix, 'tmb%linmat%ovrlp%matrix', subname)

  call tmb_overlap_onsite(iproc, nproc, at, tmb, rxyz)

  if (iproc==0) then
     if(iformat == WF_FORMAT_PLAIN) then
        open(99, file=filename//'overlap_onsite.bin', status='unknown',form='formatted')
     else
        open(99, file=filename//'overlap_onsite.bin', status='unknown',form='unformatted')
     end if

     do iorb=1,tmb%linmat%denskern%full_dim1
        iat=tmb%orbs%onwhichatom(iorb)
        do jorb=1,tmb%linmat%denskern%full_dim1
           jat=tmb%orbs%onwhichatom(jorb)
           if (iformat == WF_FORMAT_PLAIN) then
              write(99,'(2(i6,1x),e19.12,2(1x,i6))') iorb,jorb,tmb%linmat%ovrlp%matrix(iorb,jorb),iat,jat
           else
              write(99) iorb,jorb,tmb%linmat%ovrlp%matrix(iorb,jorb),iat,jat
           end if
        end do
     end do

  end if

  i_all = -product(shape(tmb%linmat%ovrlp%matrix))*kind(tmb%linmat%ovrlp%matrix)
  deallocate(tmb%linmat%ovrlp%matrix,stat=i_stat)
  call memocc(i_stat,i_all,'tmb%linmat%ovrlp%matrix',subname)

end subroutine write_linear_matrices


subroutine tmb_overlap_onsite(iproc, nproc, at, tmb, rxyz)

  use module_base
  use module_types
  use module_interfaces
  use module_fragments
  implicit none

  ! Calling arguments
  integer,intent(in) :: iproc, nproc
  type(atoms_data), intent(inout) :: at
  type(DFT_wavefunction),intent(in):: tmb
  real(gp),dimension(3,at%astruct%nat),intent(in) :: rxyz

  ! Local variables
  logical :: reformat
  integer :: iorb,i_stat,i_all,jstart,jstart_tmp
  integer :: iiorb,ilr,iiat,j,iis1,iie1,i1
  integer :: ilr_tmp,iiat_tmp,ndim_tmp,ndim,norb_tmp
  integer, dimension(3) :: ns,ns_tmp,n,n_tmp
  real(gp), dimension(3) :: centre_old_box, centre_new_box, da
  real(wp), dimension(:,:,:,:,:,:), allocatable :: phigold
  real(wp), dimension(:), pointer :: psi_tmp, psit_c_tmp, psit_f_tmp, norm
  integer, dimension(0:7) :: reformat_reason
  type(collective_comms) :: collcom_tmp
  type(local_zone_descriptors) :: lzd_tmp
  real(gp) :: tol
  character(len=*),parameter:: subname='tmb_overlap_onsite'
  type(fragment_transformation) :: frag_trans

  ! move all psi into psi_tmp all centred in the same place and calculate overlap matrix
  tol=1.d-3
  reformat_reason=0

  !arbitrarily pick the middle one as assuming it'll be near the centre of structure
  !and therefore have large fine grid
  norb_tmp=tmb%orbs%norb/2
  ilr_tmp=tmb%orbs%inwhichlocreg(norb_tmp) 
  iiat_tmp=tmb%orbs%onwhichatom(norb_tmp)

  ! find biggest instead
  !do ilr=1,tmb%lzr%nlr
  !  if (tmb%lzd%llr(ilr)%wfd%nvctr_c
  !end do

  ! Determine size of phi_old and phi
  ndim_tmp=0
  ndim=0
  do iorb=1,tmb%orbs%norbp
      iiorb=tmb%orbs%isorb+iorb
      ilr=tmb%orbs%inwhichlocreg(iiorb)
      ndim=ndim+tmb%lzd%llr(ilr)%wfd%nvctr_c+7*tmb%lzd%llr(ilr)%wfd%nvctr_f
      ndim_tmp=ndim_tmp+tmb%lzd%llr(ilr_tmp)%wfd%nvctr_c+7*tmb%lzd%llr(ilr_tmp)%wfd%nvctr_f
  end do

  ! should integrate bettwer with existing reformat routines, but restart needs tidying anyway
  allocate(psi_tmp(ndim_tmp),stat=i_stat)
  call memocc(i_stat,psi_tmp,'psi_tmp',subname)

  jstart=1
  jstart_tmp=1
  do iorb=1,tmb%orbs%norbp
      iiorb=tmb%orbs%isorb+iorb
      ilr=tmb%orbs%inwhichlocreg(iiorb)
      iiat=tmb%orbs%onwhichatom(iiorb)

      n(1)=tmb%lzd%Llr(ilr)%d%n1
      n(2)=tmb%lzd%Llr(ilr)%d%n2
      n(3)=tmb%lzd%Llr(ilr)%d%n3
      n_tmp(1)=tmb%lzd%Llr(ilr_tmp)%d%n1
      n_tmp(2)=tmb%lzd%Llr(ilr_tmp)%d%n2
      n_tmp(3)=tmb%lzd%Llr(ilr_tmp)%d%n3
      ns(1)=tmb%lzd%Llr(ilr)%ns1
      ns(2)=tmb%lzd%Llr(ilr)%ns2
      ns(3)=tmb%lzd%Llr(ilr)%ns3
      ns_tmp(1)=tmb%lzd%Llr(ilr_tmp)%ns1
      ns_tmp(2)=tmb%lzd%Llr(ilr_tmp)%ns2
      ns_tmp(3)=tmb%lzd%Llr(ilr_tmp)%ns3

      !theta=0.d0*(4.0_gp*atan(1.d0)/180.0_gp)
      !newz=(/1.0_gp,0.0_gp,0.0_gp/)
      !centre_old(:)=rxyz(:,iiat)
      !centre_new(:)=rxyz(:,iiat_tmp)
      !shift(:)=centre_new(:)-centre_old(:)

      allocate(frag_trans%discrete_operations(0),stat=i_stat)
      call memocc(i_stat,frag_trans%discrete_operations,'frag_trans%discrete_operations',subname)

      frag_trans%theta=0.0d0*(4.0_gp*atan(1.d0)/180.0_gp)
      frag_trans%rot_axis=(/1.0_gp,0.0_gp,0.0_gp/)
      frag_trans%rot_center(:)=rxyz(:,iiat)
      frag_trans%rot_center_new(:)=rxyz(:,iiat_tmp)

      call reformat_check(reformat,reformat_reason,tol,at,tmb%lzd%hgrids,tmb%lzd%hgrids,&
           tmb%lzd%llr(ilr)%wfd%nvctr_c,tmb%lzd%llr(ilr)%wfd%nvctr_c,&
           tmb%lzd%llr(ilr_tmp)%wfd%nvctr_c,tmb%lzd%llr(ilr_tmp)%wfd%nvctr_c,&
           n,n_tmp,ns,ns_tmp,frag_trans,centre_old_box,centre_new_box,da)  

      if (.not. reformat) then ! copy psi into psi_tmp
          do j=1,tmb%lzd%llr(ilr_tmp)%wfd%nvctr_c
              psi_tmp(jstart_tmp)=tmb%psi(jstart)
              jstart=jstart+1
              jstart_tmp=jstart_tmp+1
          end do
          do j=1,7*tmb%lzd%llr(ilr)%wfd%nvctr_f-6,7
              psi_tmp(jstart_tmp+0)=tmb%psi(jstart+0)
              psi_tmp(jstart_tmp+1)=tmb%psi(jstart+1)
              psi_tmp(jstart_tmp+2)=tmb%psi(jstart+2)
              psi_tmp(jstart_tmp+3)=tmb%psi(jstart+3)
              psi_tmp(jstart_tmp+4)=tmb%psi(jstart+4)
              psi_tmp(jstart_tmp+5)=tmb%psi(jstart+5)
              psi_tmp(jstart_tmp+6)=tmb%psi(jstart+6)
              jstart=jstart+7
              jstart_tmp=jstart_tmp+7
          end do
   
      else
          allocate(phigold(0:n(1),2,0:n(2),2,0:n(3),2+ndebug),stat=i_stat)
          call memocc(i_stat,phigold,'phigold',subname)

          call psi_to_psig(n,tmb%lzd%llr(ilr)%wfd%nvctr_c,tmb%lzd%llr(ilr)%wfd%nvctr_f,&
               tmb%lzd%llr(ilr)%wfd%nseg_c,tmb%lzd%llr(ilr)%wfd%nseg_f,&
               tmb%lzd%llr(ilr)%wfd%keyvloc,tmb%lzd%llr(ilr)%wfd%keygloc,jstart,tmb%psi(jstart),phigold)

          call reformat_one_supportfunction(tmb%lzd%llr(ilr_tmp)%wfd,tmb%lzd%llr(ilr_tmp)%geocode,&
               tmb%lzd%hgrids,n,phigold,tmb%lzd%hgrids,n_tmp,centre_old_box,centre_new_box,da,&
               frag_trans,psi_tmp(jstart_tmp))


          i_all = -product(shape(frag_trans%discrete_operations))*kind(frag_trans%discrete_operations)
          deallocate(frag_trans%discrete_operations,stat=i_stat)
          call memocc(i_stat,i_all,'frag_trans%discrete_operations',subname)

          jstart_tmp=jstart_tmp+tmb%lzd%llr(ilr_tmp)%wfd%nvctr_c+7*tmb%lzd%llr(ilr_tmp)%wfd%nvctr_f
   
          i_all=-product(shape(phigold))*kind(phigold)
          deallocate(phigold,stat=i_stat)
          call memocc(i_stat,i_all,'phigold',subname)

      end if

  end do

  call print_reformat_summary(iproc,reformat_reason)

  ! now that they are all in one lr, need to calculate overlap matrix
  ! make lzd_tmp contain all identical lrs
  lzd_tmp%linear=tmb%lzd%linear
  lzd_tmp%nlr=tmb%lzd%nlr
  lzd_tmp%lintyp=tmb%lzd%lintyp
  lzd_tmp%ndimpotisf=tmb%lzd%ndimpotisf
  lzd_tmp%hgrids(:)=tmb%lzd%hgrids(:)

  call nullify_locreg_descriptors(lzd_tmp%glr)
  call copy_locreg_descriptors(tmb%lzd%glr, lzd_tmp%glr, subname)

  iis1=lbound(tmb%lzd%llr,1)
  iie1=ubound(tmb%lzd%llr,1)
  allocate(lzd_tmp%llr(iis1:iie1), stat=i_stat)
  do i1=iis1,iie1
     call nullify_locreg_descriptors(lzd_tmp%llr(i1))
     call copy_locreg_descriptors(tmb%lzd%llr(ilr_tmp), lzd_tmp%llr(i1), subname)
  end do

  call nullify_collective_comms(collcom_tmp)
  call init_collective_comms(iproc, nproc, ndim_tmp, tmb%orbs, lzd_tmp, collcom_tmp)

  allocate(psit_c_tmp(sum(collcom_tmp%nrecvcounts_c)), stat=i_stat)
  call memocc(i_stat, psit_c_tmp, 'psit_c_tmp', subname)

  allocate(psit_f_tmp(7*sum(collcom_tmp%nrecvcounts_f)), stat=i_stat)
  call memocc(i_stat, psit_f_tmp, 'psit_f_tmp', subname)

  call transpose_localized(iproc, nproc, ndim_tmp, tmb%orbs, collcom_tmp, &
       psi_tmp, psit_c_tmp, psit_f_tmp, lzd_tmp)

  ! normalize psi
  allocate(norm(tmb%orbs%norb), stat=i_stat)
  call memocc(i_stat, norm, 'norm', subname)
  call normalize_transposed(iproc, nproc, tmb%orbs, collcom_tmp, psit_c_tmp, psit_f_tmp, norm)
  i_all = -product(shape(norm))*kind(norm)
  deallocate(norm,stat=i_stat)
  call memocc(i_stat,i_all,'norm',subname)

  call calculate_pulay_overlap(iproc, nproc, tmb%orbs, tmb%orbs, collcom_tmp, collcom_tmp, &
       psit_c_tmp, psit_c_tmp, psit_f_tmp, psit_f_tmp, tmb%linmat%ovrlp%matrix)

  call deallocate_collective_comms(collcom_tmp, subname)
  call deallocate_local_zone_descriptors(lzd_tmp, subname)

  i_all = -product(shape(psit_c_tmp))*kind(psit_c_tmp)
  deallocate(psit_c_tmp,stat=i_stat)
  call memocc(i_stat,i_all,'psit_c_tmp',subname)

  i_all = -product(shape(psit_f_tmp))*kind(psit_f_tmp)
  deallocate(psit_f_tmp,stat=i_stat)
  call memocc(i_stat,i_all,'psit_f_tmp',subname)

  i_all = -product(shape(psi_tmp))*kind(psi_tmp)
  deallocate(psi_tmp,stat=i_stat)
  call memocc(i_stat,i_all,'psi_tmp',subname)

END SUBROUTINE tmb_overlap_onsite


!> Write all my wavefunctions in files by calling writeonewave
subroutine writemywaves_linear(iproc,filename,iformat,npsidim,Lzd,orbs,nksorb,at,rxyz,psi,coeff)
  use module_types
  use module_base
  use yaml_output
  use module_interfaces, except_this_one => writeonewave
  implicit none
  integer, intent(in) :: iproc,iformat,npsidim,nksorb
  !integer, intent(in) :: norb   !< number of orbitals, not basis functions
  type(atoms_data), intent(in) :: at
  type(orbitals_data), intent(in) :: orbs         !< orbs describing the basis functions
  type(local_zone_descriptors), intent(in) :: Lzd
  real(gp), dimension(3,at%astruct%nat), intent(in) :: rxyz
  real(wp), dimension(npsidim), intent(in) :: psi  ! Should be the real linear dimension and not the global
  real(wp), dimension(orbs%norb,orbs%norb), intent(in) :: coeff
  character(len=*), intent(in) :: filename
  !Local variables
  integer :: ncount1,ncount_rate,ncount_max,iorb,ncount2,iorb_out,ispinor,ilr,shift,ii,iat
  integer :: jorb,jlr
  real(kind=4) :: tr0,tr1
  real(kind=8) :: tel

  if (iproc == 0) call yaml_map('Write wavefunctions to file', trim(filename)//'.*')
  !if (iproc == 0) write(*,"(1x,A,A,a)") "Write wavefunctions to file: ", trim(filename),'.*'

  if (iformat == WF_FORMAT_ETSF) then
      stop 'Linear scaling with ETSF writing not implemented yet'
!     call write_waves_etsf(iproc,filename,orbs,n1,n2,n3,hx,hy,hz,at,rxyz,wfd,psi)
  else
     call cpu_time(tr0)
     call system_clock(ncount1,ncount_rate,ncount_max)

     ! Write the TMBs in the Plain BigDFT files.
     ! Use same ordering as posinp and llr generation
     ii = 0
     do iat = 1, at%astruct%nat
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
                 call writeonewave_linear(99,(iformat == WF_FORMAT_PLAIN),iorb_out,&
                    & Lzd%Llr(ilr)%d%n1,Lzd%Llr(ilr)%d%n2,Lzd%Llr(ilr)%d%n3,&
                    & Lzd%Llr(ilr)%ns1,Lzd%Llr(ilr)%ns2,Lzd%Llr(ilr)%ns3,& 
                    & Lzd%hgrids(1),Lzd%hgrids(2),Lzd%hgrids(3), &
                    & Lzd%Llr(ilr)%locregCenter,Lzd%Llr(ilr)%locrad, 4, 0.0d0, &  !put here the real potentialPrefac and Order
                    & at%astruct%nat,rxyz,Lzd%Llr(ilr)%wfd%nseg_c,Lzd%Llr(ilr)%wfd%nvctr_c,&
                    & Lzd%Llr(ilr)%wfd%keygloc,Lzd%Llr(ilr)%wfd%keyvloc, &
                    & Lzd%Llr(ilr)%wfd%nseg_f,Lzd%Llr(ilr)%wfd%nvctr_f,&
                    & Lzd%Llr(ilr)%wfd%keygloc(1,Lzd%Llr(ilr)%wfd%nseg_c+1), &
                    & Lzd%Llr(ilr)%wfd%keyvloc(Lzd%Llr(ilr)%wfd%nseg_c+1), &
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
      call writeLinearCoefficients(99,(iformat == WF_FORMAT_PLAIN),at%astruct%nat,rxyz,orbs%norb,nksorb,coeff,orbs%eval)
      close(99)
    end if
     call cpu_time(tr1)
     call system_clock(ncount2,ncount_rate,ncount_max)
     tel=dble(ncount2-ncount1)/dble(ncount_rate)
     if (iproc == 0) then
        call yaml_open_sequence('Write Waves Time')
        call yaml_sequence(advance='no')
        call yaml_open_map(flow=.true.)
        call yaml_map('Process',iproc)
        call yaml_map('Timing',(/ real(tr1-tr0,kind=8),tel /),fmt='(1pe10.3)')
        call yaml_close_map()
        call yaml_close_sequence()
     end if
     !write(*,'(a,i4,2(1x,1pe10.3))') '- WRITE WAVES TIME',iproc,tr1-tr0,tel
     !write(*,'(a,1x,i0,a)') '- iproc',iproc,' finished writing waves'
  end if

END SUBROUTINE writemywaves_linear


subroutine readonewave_linear(unitwf,useFormattedInput,iorb,iproc,n,ns,&
     & hgrids,at,wfd,rxyz_old,rxyz,locrad,locregCenter,confPotOrder,&
     & confPotprefac,psi,eval,onwhichatom,lr,glr,reformat_reason)
  use module_base
  use module_types
  use internal_io
  use module_interfaces
  use yaml_output
  use module_fragments
  implicit none
  logical, intent(in) :: useFormattedInput
  integer, intent(in) :: unitwf,iorb,iproc
  integer, dimension(3), intent(in) :: n,ns
  type(wavefunctions_descriptors), intent(in) :: wfd
  type(atoms_data), intent(in) :: at
  real(gp), dimension(3), intent(in) :: hgrids
  real(gp), dimension(3,at%astruct%nat), intent(in) :: rxyz
  integer, intent(out) :: confPotOrder
  real(gp), intent(out) :: locrad, confPotprefac
  real(wp), intent(out) :: eval
  real(gp), dimension(3), intent(out) :: locregCenter
  real(gp), dimension(3,at%astruct%nat), intent(out) :: rxyz_old
  real(wp), dimension(:), pointer :: psi
  integer, dimension(*), intent(in) :: onwhichatom
  type(locreg_descriptors), intent(in) :: lr, glr
  integer, dimension(0:7), intent(out) :: reformat_reason

  !local variables
  character(len=*), parameter :: subname='readonewave_linear'
  character(len = 256) :: error
  logical :: lstat,reformat
  integer :: iorb_old,nvctr_c_old,nvctr_f_old,i_all,iiat
  integer :: i1,i2,i3,iel,i_stat,onwhichatom_tmp
  integer, dimension(3) :: ns_old,n_old
  real(gp) :: tol
  real(gp), dimension(3) :: hgrids_old, centre_old_box, centre_new_box, da
  real(gp) :: tt,t1,t2,t3,t4,t5,t6,t7
  real(wp), dimension(:,:,:,:,:,:), allocatable :: psigold
  type(fragment_transformation) :: frag_trans
  ! DEBUG
  character(len=12) :: orbname
  real(wp), dimension(:), allocatable :: gpsi


  call io_read_descr_linear(unitwf, useFormattedInput, iorb_old, eval, n_old(1), n_old(2), n_old(3), &
       ns_old(1), ns_old(2), ns_old(3), hgrids_old, lstat, error, onwhichatom_tmp, locrad, &
       locregCenter, confPotOrder, confPotprefac, nvctr_c_old, nvctr_f_old, at%astruct%nat, rxyz_old)

  if (.not. lstat) call io_error(trim(error))
  if (iorb_old /= iorb) stop 'readonewave_linear'

  iiat=onwhichatom(iorb)
  tol=1.d-3

  frag_trans%theta=20.0d0*(4.0_gp*atan(1.d0)/180.0_gp)
  frag_trans%rot_axis=(/1.0_gp,0.0_gp,0.0_gp/)
  frag_trans%rot_center(:)=(/7.8d0,11.8d0,11.6d0/)
  frag_trans%rot_center_new(:)=(/7.8d0,11.2d0,11.8d0/)

  allocate(frag_trans%discrete_operations(0),stat=i_stat)
  call memocc(i_stat,frag_trans%discrete_operations,'frag_trans%discrete_operations',subname)

  call reformat_check(reformat,reformat_reason,tol,at,hgrids,hgrids_old,&
       nvctr_c_old,nvctr_f_old,wfd%nvctr_c,wfd%nvctr_f,&
       n_old,n,ns_old,ns,frag_trans,centre_old_box,centre_new_box,da)  

  if (.not. reformat) then
     call read_psi_compress(unitwf, useFormattedInput, nvctr_c_old, nvctr_f_old, psi, lstat, error)
     if (.not. lstat) call io_error(trim(error))
  else
     ! add derivative functions at a later date? (needs orbs and lzd)
     allocate(psigold(0:n_old(1),2,0:n_old(2),2,0:n_old(3),2+ndebug),stat=i_stat)
     call memocc(i_stat,psigold,'psigold',subname)

     call razero(8*(n_old(1)+1)*(n_old(2)+1)*(n_old(3)+1),psigold)
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

     ! NB assuming here geocode is the same in glr and llr
     call reformat_one_supportfunction(wfd,at%astruct%geocode,hgrids_old,n_old,psigold,hgrids,n, &
         centre_old_box,centre_new_box,da,frag_trans,psi)

     i_all = -product(shape(frag_trans%discrete_operations))*kind(frag_trans%discrete_operations)
     deallocate(frag_trans%discrete_operations,stat=i_stat)
     call memocc(i_stat,i_all,'frag_trans%discrete_operations',subname)

     i_all=-product(shape(psigold))*kind(psigold)
     deallocate(psigold,stat=i_stat)
     call memocc(i_stat,i_all,'psigold',subname)

  endif

  ! DEBUG - plot in global box - CHECK WITH REFORMAT ETC IN LRs
  allocate (gpsi(glr%wfd%nvctr_c+7*glr%wfd%nvctr_f),stat=i_stat)
  call memocc(i_stat,gpsi,'gpsi',subname)

  call to_zero(glr%wfd%nvctr_c+7*glr%wfd%nvctr_f,gpsi)
  call Lpsi_to_global2(iproc, lr%wfd%nvctr_c+7*lr%wfd%nvctr_f, glr%wfd%nvctr_c+7*glr%wfd%nvctr_f, &
       1, 1, 1, glr, lr, psi, gpsi)

  write(orbname,*) iorb
  call plot_wf(trim(adjustl(orbname)),1,at,1.0_dp,glr,hgrids(1),hgrids(2),hgrids(3),rxyz,gpsi)
  !call plot_wf(trim(adjustl(orbname)),1,at,1.0_dp,lr,hx,hy,hz,rxyz,psi)

  i_all=-product(shape(gpsi))*kind(gpsi)
  deallocate(gpsi,stat=i_stat)
  call memocc(i_stat,i_all,'gpsi',subname)
  ! END DEBUG 


END SUBROUTINE readonewave_linear


subroutine io_read_descr_linear(unitwf, formatted, iorb_old, eval, n_old1, n_old2, n_old3, &
       & ns_old1, ns_old2, ns_old3, hgrids_old, lstat, error, onwhichatom, locrad, locregCenter, &
       & confPotOrder, confPotprefac, nvctr_c_old, nvctr_f_old, nat, rxyz_old)
    use module_base
    use module_types
    use internal_io
    use yaml_output
    implicit none

    integer, intent(in) :: unitwf
    logical, intent(in) :: formatted
    integer, intent(out) :: iorb_old
    integer, intent(out) :: n_old1, n_old2, n_old3, ns_old1, ns_old2, ns_old3
    real(gp), dimension(3), intent(out) :: hgrids_old
    logical, intent(out) :: lstat
    real(wp), intent(out) :: eval
    real(gp), intent(out) :: locrad
    real(gp), dimension(3), intent(out) :: locregCenter
    character(len =256), intent(out) :: error
    integer, intent(out) :: onwhichatom
    integer, intent(out) :: confPotOrder
    real(gp), intent(out) :: confPotprefac
    ! Optional arguments
    integer, intent(out), optional :: nvctr_c_old, nvctr_f_old
    integer, intent(in), optional :: nat
    real(gp), dimension(:,:), intent(out), optional :: rxyz_old

    integer :: i, iat, i_stat, nat_
    real(gp) :: rxyz(3)

    lstat = .false.
    write(error, "(A)") "cannot read psi description."
    if (formatted) then
       read(unitwf,*,iostat=i_stat) iorb_old,eval
       if (i_stat /= 0) return
       read(unitwf,*,iostat=i_stat) hgrids_old(1),hgrids_old(2),hgrids_old(3)
       if (i_stat /= 0) return
       read(unitwf,*,iostat=i_stat) n_old1,n_old2,n_old3
       if (i_stat /= 0) return
       read(unitwf,*,iostat=i_stat) ns_old1,ns_old2,ns_old3
       if (i_stat /= 0) return
       read(unitwf,*,iostat=i_stat) (locregCenter(i),i=1,3),onwhichatom,&
            locrad,confPotOrder, confPotprefac
       if (i_stat /= 0) return
       !call yaml_map('Reading atomic positions',nat)
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

       read(unitwf,iostat=i_stat) hgrids_old(1),hgrids_old(2),hgrids_old(3)
       if (i_stat /= 0) return
       read(unitwf,iostat=i_stat) n_old1,n_old2,n_old3
       if (i_stat /= 0) return
       read(unitwf,iostat=i_stat) ns_old1,ns_old2,ns_old3
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

subroutine io_read_descr_coeff(unitwf, formatted, norb_old, ntmb_old, &
       & lstat, error, nat, rxyz_old)
    use module_base
    use module_types
    use internal_io
    implicit none
    integer, intent(in) :: unitwf
    logical, intent(in) :: formatted
    integer, intent(out) :: norb_old, ntmb_old
    logical, intent(out) :: lstat
    character(len =256), intent(out) :: error
    ! Optional arguments
    integer, intent(in), optional :: nat
    real(gp), dimension(:,:), intent(out), optional :: rxyz_old

    integer :: i, iat, i_stat, nat_
    real(gp) :: rxyz(3)

    lstat = .false.
    write(error, "(A)") "cannot read coeff description."
    if (formatted) then
       read(unitwf,*,iostat=i_stat) ntmb_old, norb_old
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
       !read(unitwf,*,iostat=i_stat) i, iat
       !if (i_stat /= 0) return
    else
       read(unitwf,iostat=i_stat) ntmb_old, norb_old
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
       !read(unitwf,iostat=i_stat) i, iat
       !if (i_stat /= 0) return
    end if
    lstat = .true.
END SUBROUTINE io_read_descr_coeff


subroutine read_coeff_minbasis(unitwf,useFormattedInput,iproc,ntmb,norb_old,coeff,eval,nat,rxyz_old)
  use module_base
  use module_types
  use internal_io
  use module_interfaces, except_this_one => read_coeff_minbasis
  use yaml_output
  implicit none
  logical, intent(in) :: useFormattedInput
  integer, intent(in) :: unitwf,iproc,ntmb
  integer, intent(out) :: norb_old
  real(wp), dimension(ntmb,ntmb), intent(out) :: coeff
  real(wp), dimension(ntmb), intent(out) :: eval
  integer, optional, intent(in) :: nat
  real(gp), dimension(:,:), optional, intent(out) :: rxyz_old

  !local variables
  character(len = 256) :: error
  logical :: lstat
  integer :: i_stat
  integer :: ntmb_old, i1, i2,i,j,iorb,iorb_old
  real(wp) :: tt

  call io_read_descr_coeff(unitwf, useFormattedInput, norb_old, ntmb_old, &
       & lstat, error, nat, rxyz_old)
  if (.not. lstat) call io_error(trim(error))

  if (ntmb_old /= ntmb_old) then
     if (iproc == 0) write(error,"(A)") 'error in read coeffs, ntmb_old/=ntmb'
     call io_error(trim(error))
  end if

  ! read the eigenvalues
  if (useFormattedInput) then
     do iorb=1,ntmb
        read(unitwf,*,iostat=i_stat) iorb_old,eval(iorb)
        if (iorb_old /= iorb) stop 'read_coeff_minbasis'
     enddo
  else 
     do iorb=1,ntmb
        read(unitwf,iostat=i_stat) iorb_old,eval(iorb)
        if (iorb_old /= iorb) stop 'read_coeff_minbasis'
     enddo
     if (i_stat /= 0) stop 'Problem reading the eigenvalues'
  end if

  !if (iproc == 0) write(*,*) 'coefficients need NO reformatting'

  ! Now read the coefficients
  do i = 1, ntmb
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


END SUBROUTINE read_coeff_minbasis



!> Reads wavefunction from file and transforms it properly if hgrid or size of simulation cell                                                                                                                                                                                                                                                                                                                                   
!!  have changed
subroutine readmywaves_linear_new(iproc,dir_output,filename,iformat,at,tmb,rxyz_old,rxyz,&
       ref_frags,input_frag,frag_calc,orblist)
  use module_base
  use module_types
  use yaml_output
  use module_fragments
  use internal_io
  use module_interfaces, except_this_one => readmywaves_linear_new
  implicit none
  integer, intent(in) :: iproc, iformat
  type(atoms_data), intent(in) :: at
  type(DFT_wavefunction), intent(inout) :: tmb
  real(gp), dimension(3,at%astruct%nat), intent(in) :: rxyz
  real(gp), dimension(3,at%astruct%nat), intent(out) :: rxyz_old
  character(len=*), intent(in) :: dir_output, filename
  type(fragmentInputParameters), intent(in) :: input_frag
  type(system_fragment), dimension(input_frag%nfrag_ref), intent(inout) :: ref_frags
  logical, intent(in) :: frag_calc
  integer, dimension(tmb%orbs%norb), intent(in), optional :: orblist
  !Local variables
  integer :: ncount1,ncount_rate,ncount_max,ncount2
  integer :: iorb_out,ispinor,ilr,ind,i_all,i_stat,iorb_old
  integer :: confPotOrder,onwhichatom_tmp,unitwf
  real(gp) :: locrad, confPotprefac
  real(gp), dimension(3) :: locregCenter, mol_centre, mol_centre_new
  real(kind=4) :: tr0,tr1
  real(kind=8) :: tel,eval
  character(len=256) :: error, full_filename
  logical :: lstat
  character(len=*), parameter :: subname='readmywaves_linear_new'
  ! to eventually be part of the fragment structure?
  integer :: ndim_old, iiorb, ifrag, norb, isorb, ifrag_ref, isfat, iorbp, iforb, isforb, iiat, iat, itmb, jtmb, jsforb, ierr
  type(local_zone_descriptors) :: lzd_old
  real(wp), dimension(:), pointer :: psi_old
  type(phi_array), dimension(:), pointer :: phi_array_old
  type(fragment_transformation), dimension(:), pointer :: frag_trans_orb, frag_trans_frag
  real(gp), dimension(:,:), allocatable :: rxyz_ref, rxyz_new, rxyz4_ref, rxyz4_new
  real(gp), dimension(:), allocatable :: dist
  integer, dimension(:), allocatable :: ipiv
  logical :: skip

  ! DEBUG
  character(len=12) :: orbname
  real(wp), dimension(:), allocatable :: gpsi

open(16)
  call cpu_time(tr0)
  call system_clock(ncount1,ncount_rate,ncount_max)

  ! check file format
  if (iformat == WF_FORMAT_ETSF) then
     stop 'Linear scaling with ETSF writing not implemented yet'
  else if (iformat /= WF_FORMAT_BINARY .and. iformat /= WF_FORMAT_PLAIN) then
     call yaml_warning('Unknown wavefunction file format from filename.')
     stop
  end if

  ! to be fixed
  if (present(orblist)) then
     stop 'orblist no longer functional in initialize_linear_from_file due to addition of fragment calculation'
  end if

  ! lzd_old => ref_frags(onwhichfrag)%frag_basis%lzd
  ! orbs_old -> ref_frags(onwhichfrag)%frag_basis%forbs ! <- BUT problem with it not being same type
  ! phi_array_old => ref_frags(onwhichfrag)%frag_basis%phi

  ! change parallelization later, for now all procs read the same number of tmbs as before
  ! initialize fragment lzd and phi_array_old for fragment, then allocate lzd_old which points to appropriate fragment entry
  ! for now directly using lzd_old etc - less efficient if fragments are used multiple times

  ! use above information to generate lzds
  call nullify_local_zone_descriptors(lzd_old)
  call nullify_locreg_descriptors(lzd_old%glr)
  lzd_old%nlr=tmb%orbs%norb
  allocate(lzd_old%Llr(lzd_old%nlr),stat=i_stat)
  do ilr=1,lzd_old%nlr
     call nullify_locreg_descriptors(lzd_old%llr(ilr))
  end do

  ! has size of new orbs, will possibly point towards the same tmb multiple times
  allocate(phi_array_old(tmb%orbs%norbp), stat=i_stat)
  do iorbp=1,tmb%orbs%norbp
     nullify(phi_array_old(iorbp)%psig)
  end do

  !allocate(frag_trans_orb(tmb%orbs%norbp))

  unitwf=99
  isforb=0
  isfat=0
  do ifrag=1,input_frag%nfrag
     ! find reference fragment this corresponds to
     ifrag_ref=input_frag%frag_index(ifrag)
     ! loop over orbitals of this fragment
     loop_iforb: do iforb=1,ref_frags(ifrag_ref)%fbasis%forbs%norb
        loop_iorb: do iorbp=1,tmb%orbs%norbp
           iiorb=iorbp+tmb%orbs%isorb
           !ilr = ref_frags(ifrag)%fbasis%forbs%inwhichlocreg(iiorb)
           ilr=tmb%orbs%inwhichlocreg(iiorb)
           iiat=tmb%orbs%onwhichatom(iiorb)

           ! check if this ref frag orbital corresponds to the orbital we want
           if (iiorb/=iforb+isforb) cycle
           do ispinor=1,tmb%orbs%nspinor
              ! if this is a fragment calculation frag%dirname will contain fragment directory, otherwise it will be empty
              ! bit of a hack to use orbs here not forbs, but different structures so this is necessary - to clean somehow
              full_filename=trim(dir_output)//trim(input_frag%dirname(ifrag_ref))//trim(filename)

              call open_filename_of_iorb(unitwf,(iformat == WF_FORMAT_BINARY),full_filename, &
                   & tmb%orbs,iorbp,ispinor,iorb_out,iforb)
                   !& ref_frags(ifrag_ref)%fbasis%forbs,iforb,ispinor,iorb_out)
  
              ! read headers, reading lzd info directly into lzd_old, which is otherwise nullified
              call io_read_descr_linear(unitwf, (iformat == WF_FORMAT_PLAIN), iorb_old, eval, &
                   Lzd_old%Llr(ilr)%d%n1,Lzd_old%Llr(ilr)%d%n2,Lzd_old%Llr(ilr)%d%n3, &
                   Lzd_old%Llr(ilr)%ns1,Lzd_old%Llr(ilr)%ns2,Lzd_old%Llr(ilr)%ns3, lzd_old%hgrids, &
                   lstat, error, onwhichatom_tmp, Lzd_old%Llr(ilr)%locrad, Lzd_old%Llr(ilr)%locregCenter, &
                   confPotOrder, confPotprefac, Lzd_old%Llr(ilr)%wfd%nvctr_c, Lzd_old%Llr(ilr)%wfd%nvctr_f, &
                   ref_frags(ifrag_ref)%astruct_frg%nat, rxyz_old(:,isfat+1:isfat+ref_frags(ifrag_ref)%astruct_frg%nat))
                   !ref_frags(ifrag_ref)%astruct_frg%nat, rxyz_old(1,isfat+1))

              ! in general this might point to a different tmb
              allocate(phi_array_old(iorbp)%psig(0:Lzd_old%Llr(ilr)%d%n1,2,0:Lzd_old%Llr(ilr)%d%n2,2,&
                   0:Lzd_old%Llr(ilr)%d%n3,2), stat=i_stat)
              call memocc(i_stat, phi_array_old(iorbp)%psig, 'phi_array_old(iorb)%psig', subname)

              !read phig directly
              call read_psig(unitwf, (iformat == WF_FORMAT_PLAIN), Lzd_old%Llr(ilr)%wfd%nvctr_c, Lzd_old%Llr(ilr)%wfd%nvctr_f, &
                   Lzd_old%Llr(ilr)%d%n1, Lzd_old%Llr(ilr)%d%n2, Lzd_old%Llr(ilr)%d%n3, phi_array_old(iorbp)%psig, lstat, error)
              if (.not. lstat) call io_error(trim(error))

              ! DEBUG: print*,iproc,iorb,iorb+orbs%isorb,iorb_old,iorb_out

              !! define fragment transformation - should eventually be done automatically...
              !! first fragment is shifted only, hack here that second fragment should be rotated
              !if (ifrag==1) then
              !   frag_trans_orb(iorbp)%theta=0.0d0*(4.0_gp*atan(1.d0)/180.0_gp)
              !   frag_trans_orb(iorbp)%rot_axis=(/1.0_gp,0.0_gp,0.0_gp/)
              !   frag_trans_orb(iorbp)%rot_center(:)=rxyz_old(:,iiat)
              !   frag_trans_orb(iorbp)%rot_center_new(:)=rxyz(:,iiat)
              !else
              !   ! unnecessary recalculation here, to be tidied later
              !   mol_centre=0.0d0
              !   mol_centre_new=0.0d0
              !   do iat=1,ref_frags(ifrag_ref)%astruct_frg%nat
              !      mol_centre(:)=mol_centre(:)+rxyz_old(:,isfat+iat)
              !      mol_centre_new(:)=mol_centre_new(:)+rxyz(:,isfat+iat)
              !   end do
              !   mol_centre=mol_centre/real(ref_frags(ifrag_ref)%astruct_frg%nat,gp)
              !   mol_centre_new=mol_centre_new/real(ref_frags(ifrag_ref)%astruct_frg%nat,gp)
              !   frag_trans_orb(iorbp)%theta=30.0d0*(4.0_gp*atan(1.d0)/180.0_gp)
              !   frag_trans_orb(iorbp)%rot_axis=(/1.0_gp,0.0_gp,0.0_gp/)
              !   frag_trans_orb(iorbp)%rot_center(:)=mol_centre(:) ! take as average for now
              !   frag_trans_orb(iorbp)%rot_center_new(:)=mol_centre_new(:) ! to get shift, mol is rigidly shifted so could take any, rather than centre
              !end if

              !write(*,'(a,x,2(i2,x),4(f5.2,x),6(f7.3,x))'),'trans',ifrag,iiorb,frag_trans_orb(iorbp)%theta,&
              !     frag_trans_orb(iorbp)%rot_axis, &
              !     frag_trans_orb(iorbp)%rot_center,frag_trans_orb(iorbp)%rot_center_new

              if (.not. lstat) then
                 call yaml_warning(trim(error))
                 stop
              end if
              if (iorb_old /= iorb_out) then
                 call yaml_warning('Initialize_linear_from_file')
                 stop
              end if
              close(unitwf)
              
           end do
        end do loop_iorb
     end do loop_iforb
     isforb=isforb+ref_frags(ifrag_ref)%fbasis%forbs%norb
     isfat=isfat+ref_frags(ifrag_ref)%astruct_frg%nat        
  end do

  ! reformat fragments
  nullify(psi_old)

!if several fragments do this, otherwise find 3 nearest neighbours (use sort in time.f90) and send rxyz arrays with 4 atoms

  if (input_frag%nfrag>1) then
     ! Find fragment transformations for each fragment, then put in frag_trans array for each orb
     allocate(frag_trans_frag(input_frag%nfrag))
     do ifrag=1,input_frag%nfrag
        nullify(frag_trans_frag(ifrag)%discrete_operations)
     end do

     isfat=0
     isforb=0
     do ifrag=1,input_frag%nfrag
        ! find reference fragment this corresponds to
        ifrag_ref=input_frag%frag_index(ifrag)

        ! check if we need this fragment transformation on this proc
        skip=.true.
        do iforb=1,ref_frags(ifrag_ref)%fbasis%forbs%norb
           do iorbp=1,tmb%orbs%norbp
              iiorb=iorbp+tmb%orbs%isorb
              ! check if this ref frag orbital corresponds to the orbital we want
              if (iiorb==iforb+isforb) then
                 skip=.false.
                 exit
              end if
           end do
        end do

        if (skip) then
           isfat=isfat+ref_frags(ifrag_ref)%astruct_frg%nat     
           isforb=isforb+ref_frags(ifrag_ref)%fbasis%forbs%norb
           cycle
        end if

        allocate(rxyz_ref(3,ref_frags(ifrag_ref)%astruct_frg%nat), stat=i_stat)
        call memocc(i_stat, rxyz_ref, 'rxyz_ref', subname)
        allocate(rxyz_new(3,ref_frags(ifrag_ref)%astruct_frg%nat), stat=i_stat)
        call memocc(i_stat, rxyz_new, 'rxyz_ref', subname)

        do iat=1,ref_frags(ifrag_ref)%astruct_frg%nat
           rxyz_new(:,iat)=rxyz(:,isfat+iat)
           rxyz_ref(:,iat)=rxyz_old(:,isfat+iat)
        end do

        ! use center of fragment for now, could later change to center of symmetry
        frag_trans_frag(ifrag)%rot_center=frag_center(ref_frags(ifrag_ref)%astruct_frg%nat,rxyz_ref)
        frag_trans_frag(ifrag)%rot_center_new=frag_center(ref_frags(ifrag_ref)%astruct_frg%nat,rxyz_new)

        ! shift rxyz wrt center of rotation
        do iat=1,ref_frags(ifrag_ref)%astruct_frg%nat
           rxyz_ref(:,iat)=rxyz_ref(:,iat)-frag_trans_frag(ifrag)%rot_center
           rxyz_new(:,iat)=rxyz_new(:,iat)-frag_trans_frag(ifrag)%rot_center_new
        end do

        call find_frag_trans(ref_frags(ifrag_ref)%astruct_frg%nat,rxyz_ref,rxyz_new,frag_trans_frag(ifrag))

        i_all = -product(shape(rxyz_ref))*kind(rxyz_ref)
        deallocate(rxyz_ref,stat=i_stat)
        call memocc(i_stat,i_all,'rxyz_ref',subname)
        i_all = -product(shape(rxyz_new))*kind(rxyz_new)
        deallocate(rxyz_new,stat=i_stat)
        call memocc(i_stat,i_all,'rxyz_new',subname)

        write(*,'(A,I3,1x,I3,1x,3(F12.6,1x),F12.6)') 'ifrag,ifrag_ref,rot_axis,theta',&
             ifrag,ifrag_ref,frag_trans_frag(ifrag)%rot_axis,frag_trans_frag(ifrag)%theta/(4.0_gp*atan(1.d0)/180.0_gp)

        isfat=isfat+ref_frags(ifrag_ref)%astruct_frg%nat     
        isforb=isforb+ref_frags(ifrag_ref)%fbasis%forbs%norb
     end do

     allocate(frag_trans_orb(tmb%orbs%norbp))
     isforb=0
     isfat=0
     do ifrag=1,input_frag%nfrag
        ! find reference fragment this corresponds to
        ifrag_ref=input_frag%frag_index(ifrag)
        ! loop over orbitals of this fragment
        do iforb=1,ref_frags(ifrag_ref)%fbasis%forbs%norb
           do iorbp=1,tmb%orbs%norbp
              iiorb=iorbp+tmb%orbs%isorb
              ! check if this ref frag orbital corresponds to the orbital we want
              if (iiorb/=iforb+isforb) cycle

              allocate(frag_trans_orb(iorbp)%discrete_operations(size(frag_trans_frag(ifrag)%discrete_operations)),stat=i_stat)
              call memocc(i_stat,frag_trans_orb(iorbp)%discrete_operations,'frag_trans_orb(iorbp)%discrete_operations',subname)

              frag_trans_orb(iorbp)%rot_center=frag_trans_frag(ifrag)%rot_center
              frag_trans_orb(iorbp)%rot_center_new=frag_trans_frag(ifrag)%rot_center_new
              frag_trans_orb(iorbp)%rot_axis=(frag_trans_frag(ifrag)%rot_axis)
              frag_trans_orb(iorbp)%theta=frag_trans_frag(ifrag)%theta
              call dcopy(size(frag_trans_frag(ifrag)%discrete_operations),frag_trans_frag(ifrag)%discrete_operations,1,&
                   frag_trans_orb(iorbp)%discrete_operations,1)

              !write(*,'(a,x,2(i2,x),4(f5.2,x),6(f7.3,x))'),'trans2',ifrag,iiorb,frag_trans_orb(iorbp)%theta,&
              !     frag_trans_orb(iorbp)%rot_axis, &
              !     frag_trans_orb(iorbp)%rot_center,frag_trans_orb(iorbp)%rot_center_new
           end do
        end do
        isforb=isforb+ref_frags(ifrag_ref)%fbasis%forbs%norb
        isfat=isfat+ref_frags(ifrag_ref)%astruct_frg%nat     

        ! not associated on every proc
        if (associated(frag_trans_frag(ifrag)%discrete_operations)) then
           i_all = -product(shape(frag_trans_frag(ifrag)%discrete_operations))*kind(frag_trans_frag(ifrag)%discrete_operations)
           deallocate(frag_trans_frag(ifrag)%discrete_operations,stat=i_stat)
           call memocc(i_stat,i_all,'frag_trans_frag(ifrag)%discrete_operations',subname)
        end if
     end do

     deallocate(frag_trans_frag)
  else
     ! only 1 'fragment', calculate rotation/shift atom wise, using nearest neighbours
     allocate(frag_trans_orb(tmb%orbs%norbp))

     allocate(rxyz4_ref(3,min(4,ref_frags(ifrag_ref)%astruct_frg%nat)), stat=i_stat)
     call memocc(i_stat, rxyz4_ref, 'rxyz4_ref', subname)
     allocate(rxyz4_new(3,min(4,ref_frags(ifrag_ref)%astruct_frg%nat)), stat=i_stat)
     call memocc(i_stat, rxyz4_new, 'rxyz4_ref', subname)

     isforb=0
     isfat=0
     do ifrag=1,input_frag%nfrag
        ! find reference fragment this corresponds to
        ifrag_ref=input_frag%frag_index(ifrag)

        allocate(rxyz_ref(3,ref_frags(ifrag_ref)%astruct_frg%nat), stat=i_stat)
        call memocc(i_stat, rxyz_ref, 'rxyz_ref', subname)
        allocate(rxyz_new(3,ref_frags(ifrag_ref)%astruct_frg%nat), stat=i_stat)
        call memocc(i_stat, rxyz_new, 'rxyz_ref', subname)
        allocate(dist(ref_frags(ifrag_ref)%astruct_frg%nat), stat=i_stat)
        call memocc(i_stat, dist, 'dist', subname)
        allocate(ipiv(ref_frags(ifrag_ref)%astruct_frg%nat), stat=i_stat)
        call memocc(i_stat, ipiv, 'ipiv', subname)

        ! loop over orbitals of this fragment
        do iforb=1,ref_frags(ifrag_ref)%fbasis%forbs%norb
           do iorbp=1,tmb%orbs%norbp
              iiorb=iorbp+tmb%orbs%isorb
              if (iiorb/=iforb+isforb) cycle

              do iat=1,ref_frags(ifrag_ref)%astruct_frg%nat
                 rxyz_new(:,iat)=rxyz(:,isfat+iat)
                 rxyz_ref(:,iat)=rxyz_old(:,isfat+iat)
              end do

              iiat=tmb%orbs%onwhichatom(iiorb)

              ! use atom position
              frag_trans_orb(iorbp)%rot_center=rxyz_old(:,iiat)
              frag_trans_orb(iorbp)%rot_center_new=rxyz(:,iiat)

              ! shift rxyz wrt center of rotation
              do iat=1,ref_frags(ifrag_ref)%astruct_frg%nat
                 rxyz_ref(:,iat)=rxyz_ref(:,iat)-frag_trans_orb(iorbp)%rot_center
                 rxyz_new(:,iat)=rxyz_new(:,iat)-frag_trans_orb(iorbp)%rot_center_new
              end do

              ! find distances from this atom
              do iat=1,ref_frags(ifrag_ref)%astruct_frg%nat
                   dist(iat)=-dsqrt(rxyz_ref(1,iat)**2+rxyz_ref(2,iat)**2+rxyz_ref(3,iat)**2)
              end do             

              ! sort atoms into neighbour order
              call sort_positions(ref_frags(ifrag_ref)%astruct_frg%nat,dist,ipiv)

              ! take atom and 3 nearest neighbours
              do iat=1,min(4,ref_frags(ifrag_ref)%astruct_frg%nat)
                 rxyz4_ref(:,iat)=rxyz_ref(:,ipiv(iat))
                 rxyz4_new(:,iat)=rxyz_new(:,ipiv(iat))
              end do

              call find_frag_trans(min(4,ref_frags(ifrag_ref)%astruct_frg%nat),rxyz4_ref,rxyz4_new,frag_trans_orb(iorbp))

              write(*,'(A,I3,1x,I3,1x,3(F12.6,1x),F12.6)') 'ifrag,iorb,rot_axis,theta',&
                   ifrag,iiorb,frag_trans_orb(iorbp)%rot_axis,frag_trans_orb(iorbp)%theta/(4.0_gp*atan(1.d0)/180.0_gp)

           end do
        end do
        isforb=isforb+ref_frags(ifrag_ref)%fbasis%forbs%norb
        isfat=isfat+ref_frags(ifrag_ref)%astruct_frg%nat     

        i_all = -product(shape(ipiv))*kind(ipiv)
        deallocate(ipiv,stat=i_stat)
        call memocc(i_stat,i_all,'ipiv',subname)
        i_all = -product(shape(dist))*kind(dist)
        deallocate(dist,stat=i_stat)
        call memocc(i_stat,i_all,'dist',subname)
        i_all = -product(shape(rxyz_ref))*kind(rxyz_ref)
        deallocate(rxyz_ref,stat=i_stat)
        call memocc(i_stat,i_all,'rxyz_ref',subname)
        i_all = -product(shape(rxyz_new))*kind(rxyz_new)
        deallocate(rxyz_new,stat=i_stat)
        call memocc(i_stat,i_all,'rxyz_new',subname)

     end do

     i_all = -product(shape(rxyz4_ref))*kind(rxyz4_ref)
     deallocate(rxyz4_ref,stat=i_stat)
     call memocc(i_stat,i_all,'rxyz4_ref',subname)
     i_all = -product(shape(rxyz4_new))*kind(rxyz4_new)
     deallocate(rxyz4_new,stat=i_stat)
     call memocc(i_stat,i_all,'rxyz4_new',subname)

  end if
close(16)

  call reformat_supportfunctions(iproc,at,rxyz_old,rxyz,.false.,tmb,ndim_old,lzd_old,frag_trans_orb,psi_old,phi_array_old)

  deallocate(frag_trans_orb)

  do iorbp=1,tmb%orbs%norbp
     !nullify/deallocate here as appropriate, in future may keep
     i_all = -product(shape(phi_array_old(iorbp)%psig))*kind(phi_array_old(iorbp)%psig)
     deallocate(phi_array_old(iorbp)%psig,stat=i_stat)
     call memocc(i_stat,i_all,'phi_array_old(iorbp)%psig',subname)
  end do

  deallocate(phi_array_old,stat=i_stat)
  call deallocate_local_zone_descriptors(lzd_old,subname)

  ! DEBUG - plot in global box - CHECK WITH REFORMAT ETC IN LRs
  ind=1
  allocate (gpsi(tmb%Lzd%glr%wfd%nvctr_c+7*tmb%Lzd%glr%wfd%nvctr_f),stat=i_stat)
  call memocc(i_stat,gpsi,'gpsi',subname)
  do iorbp=1,tmb%orbs%norbp
     iiorb=iorbp+tmb%orbs%isorb
     ilr = tmb%orbs%inwhichlocreg(iiorb)

     call to_zero(tmb%Lzd%glr%wfd%nvctr_c+7*tmb%Lzd%glr%wfd%nvctr_f,gpsi)
     call Lpsi_to_global2(iproc, tmb%Lzd%Llr(ilr)%wfd%nvctr_c+7*tmb%Lzd%Llr(ilr)%wfd%nvctr_f, &
          tmb%Lzd%glr%wfd%nvctr_c+7*tmb%Lzd%glr%wfd%nvctr_f, &
          1, 1, 1, tmb%Lzd%glr, tmb%Lzd%Llr(ilr), tmb%psi(ind), gpsi)

     write(orbname,*) iiorb
     call plot_wf(trim(dir_output)//trim(adjustl(orbname)),1,at,1.0_dp,tmb%Lzd%glr,&
          tmb%Lzd%hgrids(1),tmb%Lzd%hgrids(2),tmb%Lzd%hgrids(3),rxyz,gpsi)
     !call plot_wf(trim(adjustl(orbname)),1,at,1.0_dp,tmb%Lzd%Llr(ilr),&
     !     tmb%Lzd%hgrids(1),tmb%Lzd%hgrids(2),tmb%Lzd%hgrids(3),rxyz,tmb%psi)
  
     ind = ind + tmb%Lzd%Llr(ilr)%wfd%nvctr_c+7*tmb%Lzd%Llr(ilr)%wfd%nvctr_f
  end do
  i_all=-product(shape(gpsi))*kind(gpsi)
  deallocate(gpsi,stat=i_stat)
  call memocc(i_stat,i_all,'gpsi',subname)
  ! END DEBUG 


  ! Read the coefficient file for each fragment and assemble total coeffs
  if(iformat == WF_FORMAT_PLAIN) then
     open(unitwf,file=filename//'_coeff.bin',status='unknown',form='formatted')
  else if(iformat == WF_FORMAT_BINARY) then
     open(unitwf,file=filename//'_coeff.bin',status='unknown',form='unformatted')
  else
     stop 'Coefficient format not implemented'
  end if

  ! coeffs should eventually go into ref_frag array and then point? or be copied to (probably copied as will deallocate frag)
  unitwf=99
  isforb=0
  do ifrag=1,input_frag%nfrag
     ! find reference fragment this corresponds to
     ifrag_ref=input_frag%frag_index(ifrag)

     full_filename=trim(dir_output)//trim(input_frag%dirname(ifrag_ref))//trim(filename)//'_coeff.bin'

     if(iformat == WF_FORMAT_PLAIN) then
        open(unitwf,file=trim(full_filename),status='unknown',form='formatted')
     else if(iformat == WF_FORMAT_BINARY) then
        open(unitwf,file=trim(full_filename),status='unknown',form='unformatted')
     else
        stop 'Coefficient format not implemented'
     end if

     call read_coeff_minbasis(unitwf,(iformat == WF_FORMAT_PLAIN),iproc,ref_frags(ifrag_ref)%fbasis%forbs%norb,&
          ref_frags(ifrag_ref)%nksorb,ref_frags(ifrag_ref)%coeff,&
          tmb%orbs%eval(isforb+1:isforb+ref_frags(ifrag_ref)%fbasis%forbs%norb))
          !tmb%orbs%eval(isforb+1)

     close(unitwf)

     isforb=isforb+ref_frags(ifrag_ref)%fbasis%forbs%norb
  end do

  ! copy from coeff fragment to global coeffs - take occupied ones first, then unoccupied
  ! (ideally reorder these by eval?)
  isforb=0
  jsforb=0
  call to_zero(tmb%orbs%norb*tmb%orbs%norb, tmb%coeff(1,1))
  do ifrag=1,input_frag%nfrag
     ! find reference fragment this corresponds to
     ifrag_ref=input_frag%frag_index(ifrag)

     do itmb=1,ref_frags(ifrag_ref)%fbasis%forbs%norb
        do jtmb=1,ref_frags(ifrag_ref)%nksorb
           tmb%coeff(isforb+itmb,jsforb+jtmb)=ref_frags(ifrag_ref)%coeff(itmb,jtmb)
        end do
     end do

     isforb=isforb+ref_frags(ifrag_ref)%fbasis%forbs%norb
     jsforb=jsforb+ref_frags(ifrag_ref)%nksorb
  end do

  isforb=0
  do ifrag=1,input_frag%nfrag
     ! find reference fragment this corresponds to
     ifrag_ref=input_frag%frag_index(ifrag)
     do itmb=1,ref_frags(ifrag_ref)%fbasis%forbs%norb
        do jtmb=ref_frags(ifrag_ref)%nksorb+1,ref_frags(ifrag_ref)%fbasis%forbs%norb
           tmb%coeff(isforb+itmb,jsforb+jtmb-ref_frags(ifrag_ref)%nksorb)=ref_frags(ifrag_ref)%coeff(itmb,jtmb)
        end do
     end do

     isforb=isforb+ref_frags(ifrag_ref)%fbasis%forbs%norb
     jsforb=jsforb+ref_frags(ifrag_ref)%fbasis%forbs%norb-ref_frags(ifrag_ref)%nksorb
  end do

open(10)
  do itmb=1,tmb%orbs%norb
     do jtmb=1,tmb%orbs%norb
        if (iproc==0) write(10,*) jtmb,itmb,tmb%coeff(jtmb,itmb)
     end do
  end do
close(10)
call mpi_barrier(bigdft_mpi%mpi_comm,ierr)

  call cpu_time(tr1)
  call system_clock(ncount2,ncount_rate,ncount_max)
  tel=dble(ncount2-ncount1)/dble(ncount_rate)

  if (iproc == 0) then
     call yaml_open_sequence('Reading Waves Time')
     call yaml_sequence(advance='no')
     call yaml_open_map(flow=.true.)
     call yaml_map('Process',iproc)
     call yaml_map('Timing',(/ real(tr1-tr0,kind=8),tel /),fmt='(1pe10.3)')
     call yaml_close_map()
     call yaml_close_sequence()
  end if
  !write(*,'(a,i4,2(1x,1pe10.3))') '- READING WAVES TIME',iproc,tr1-tr0,tel

END SUBROUTINE readmywaves_linear_new


! initializes onwhichatom, inwhichlocreg, locrad and locregcenter from file
subroutine initialize_linear_from_file(iproc,nproc,input_frag,astruct,rxyz,orbs,Lzd,iformat,&
     dir_output,filename,ref_frags,orblist)
  use module_base
  use module_types
  use module_defs
  use yaml_output
  use module_fragments
  use module_interfaces, except_this_one => initialize_linear_from_file
  implicit none
  integer, intent(in) :: iproc, nproc, iformat
  type(fragmentInputParameters),intent(in) :: input_frag
  type(atomic_structure), intent(in) :: astruct
  real(gp), dimension(3,astruct%nat), intent(in) :: rxyz
  type(orbitals_data), intent(inout) :: orbs  !< orbs related to the basis functions, inwhichlocreg and onwhichatom generated in this routine
  type(local_zone_descriptors), intent(inout) :: Lzd !< must already contain Glr and hgrids
  type(system_fragment), dimension(input_frag%nfrag_ref), intent(inout) :: ref_frags
  character(len=*), intent(in) :: filename, dir_output
  integer, dimension(orbs%norb), optional :: orblist

  !Local variables
  character(len=*), parameter :: subname='initialize_linear_from_file'
  character(len =256) :: error
  logical :: lstat
  integer :: ilr, ierr, iorb_old, iorb, ispinor, iorb_out, iforb, isforb, isfat, iiorb, iorbp, ifrag, ifrag_ref
  integer, dimension(3) :: n_old, ns_old
  integer :: i_stat, i_all, confPotOrder, iat
  real(gp), dimension(3) :: hgrids_old
  real(kind=8) :: eval, confPotprefac
  real(gp), dimension(orbs%norb):: locrad
  real(gp), dimension(3) :: locregCenter
  real(kind=8), dimension(:), allocatable :: lrad
  real(gp), dimension(:,:), allocatable :: cxyz
  character(len=256) :: full_filename

  ! to be fixed
  if (present(orblist)) then
print*,'present(orblist)',present(orblist)
     stop 'orblist no longer functional in initialize_linear_from_file due to addition of fragment calculation'
  end if

  ! NOTES:
  ! The orbs%norb family must be all constructed before this routine
  ! This can be done from the input.lin since the number of basis functions should be fixed.
  ! Fragment structure must also be fully initialized, even if this is not a fragment calculation

  call to_zero(orbs%norb,locrad(1))
  call to_zero(orbs%norb,orbs%onwhichatom(1))

  if (iformat == WF_FORMAT_ETSF) then
     stop 'Linear scaling with ETSF writing not implemented yet'
  else if (iformat == WF_FORMAT_BINARY .or. iformat == WF_FORMAT_PLAIN) then
     ! loop over fragments in current system (will still be treated as one fragment in calculation, just separate for restart)
     isforb=0
     isfat=0
     do ifrag=1,input_frag%nfrag
        ! find reference fragment this corresponds to
        ifrag_ref=input_frag%frag_index(ifrag)
        ! loop over orbitals of this fragment
        loop_iforb: do iforb=1,ref_frags(ifrag_ref)%fbasis%forbs%norb
           loop_iorb: do iorbp=1,orbs%norbp
              iiorb=iorbp+orbs%isorb
              ! check if this ref frag orbital corresponds to the orbital we want
              if (iiorb/=iforb+isforb) cycle
              do ispinor=1,orbs%nspinor

                 ! if this is a fragment calculation frag%dirname will contain fragment directory, otherwise it will be empty
                 ! bit of a hack to use orbs here not forbs, but different structures so this is necessary - to clean somehow
                 full_filename=trim(dir_output)//trim(input_frag%dirname(ifrag_ref))//trim(filename)

                 call open_filename_of_iorb(99,(iformat == WF_FORMAT_BINARY),full_filename, &
                      & orbs,iorbp,ispinor,iorb_out,iforb)
                      !& ref_frags(ifrag_ref)%fbasis%forbs,iforb,ispinor,iorb_out)
  
                 call io_read_descr_linear(99,(iformat == WF_FORMAT_PLAIN), iorb_old, eval, n_old(1), n_old(2), n_old(3), &
                      ns_old(1), ns_old(2), ns_old(3), hgrids_old, lstat, error, orbs%onwhichatom(iiorb), &
                      locrad(iiorb), locregCenter, confPotOrder, confPotprefac)

                 orbs%onwhichatom(iiorb) = orbs%onwhichatom(iiorb) + isfat

                 if (.not. lstat) then
                    call yaml_warning(trim(error))
                    stop
                 end if
                 if (iorb_old /= iorb_out) then
                    call yaml_warning('Initialize_linear_from_file')
                    stop
                 end if
                 close(99)
              
              end do
           end do loop_iorb
        end do loop_iforb
        isforb=isforb+ref_frags(ifrag_ref)%fbasis%forbs%norb
        isfat=isfat+ref_frags(ifrag_ref)%astruct_frg%nat        
     end do
  else
     call yaml_warning('Unknown wavefunction file format from filename.')
     stop
  end if

  Lzd%nlr = orbs%norb

  ! Communication of the quantities
  if (nproc > 1)  call mpiallred(orbs%onwhichatom(1),orbs%norb,MPI_SUM,bigdft_mpi%mpi_comm,ierr)
  if (nproc > 1)  call mpiallred(locrad(1),orbs%norb,MPI_SUM,bigdft_mpi%mpi_comm,ierr)

  allocate(cxyz(3,Lzd%nlr),stat=i_stat)
  call memocc(i_stat,cxyz,'cxyz',subname)
  allocate(lrad(Lzd%nlr), stat=i_stat)
  call memocc(i_stat, lrad, 'lrad', subname)

  ! Put the llr in posinp order
  ilr=0
  do iat=1,astruct%nat
     do iorb=1,orbs%norb
        if(iat == orbs%onwhichatom(iorb)) then
           ilr = ilr + 1
           cxyz(1,ilr) = rxyz(1,iat)
           cxyz(2,ilr) = rxyz(2,iat)
           cxyz(3,ilr) = rxyz(3,iat)
           lrad(ilr) = locrad(iorb)
           orbs%inwhichlocreg(iorb) = ilr
        end if
     end do
  end do
  
  ! Allocate the array of localisation regions
  allocate(lzd%Llr(lzd%nlr),stat=i_stat)
  do ilr=1,lzd%nlr
     lzd%Llr(ilr)=locreg_null()
  end do
  do ilr=1,lzd%nlr
      lzd%llr(ilr)%locrad=lrad(ilr)
      lzd%llr(ilr)%locregCenter=cxyz(:,ilr)
  end do

  i_all = -product(shape(cxyz))*kind(cxyz)
  deallocate(cxyz,stat=i_stat)
  call memocc(i_stat,i_all,'cxyz',subname)
  i_all = -product(shape(lrad))*kind(lrad)
  deallocate(lrad,stat=i_stat)
  call memocc(i_stat,i_all,'lrad',subname)

END SUBROUTINE initialize_linear_from_file


!> Copy old support functions from phi to phi_old
subroutine copy_old_supportfunctions(orbs,lzd,phi,lzd_old,phi_old)
  use module_base
  use module_types
  implicit none
  type(orbitals_data), intent(in) :: orbs
  type(local_zone_descriptors), intent(in) :: lzd
  type(local_zone_descriptors), intent(inout) :: lzd_old
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
      if (abs(tt-1.d0) > 1.d-3) then
         write(*,*)'wrong phi_old',iiorb,tt
         !stop 
      end if
  end do

  !!!deallocation
  !!i_all=-product(shape(phi))*kind(phi)
  !!deallocate(phi,stat=i_stat)
  !!call memocc(i_stat,i_all,'phi',subname)

END SUBROUTINE copy_old_supportfunctions


subroutine copy_old_coefficients(norb_tmb, coeff, coeff_old)
  use module_base
  implicit none

  ! Calling arguments
  integer,intent(in):: norb_tmb
  real(8),dimension(:,:),pointer:: coeff, coeff_old

  ! Local variables
  character(len=*),parameter:: subname='copy_old_coefficients'
  integer:: istat
!  integer:: iall

  allocate(coeff_old(norb_tmb,norb_tmb),stat=istat)
  call memocc(istat,coeff_old,'coeff_old',subname)

  call vcopy(norb_tmb*norb_tmb, coeff(1,1), 1, coeff_old(1,1), 1)

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
!NB add_derivatives must be false if we are using phi_array_old instead of psi_old and don't have the keys
subroutine reformat_supportfunctions(iproc,at,rxyz_old,rxyz,add_derivatives,tmb,ndim_old,lzd_old,&
       frag_trans,psi_old,phi_array_old)
  use module_base
  use module_types
  use module_fragments
  use module_interfaces, except_this_one=>reformat_supportfunctions
  implicit none
  integer, intent(in) :: iproc,ndim_old
  type(atoms_data), intent(in) :: at
  real(gp), dimension(3,at%astruct%nat), intent(in) :: rxyz,rxyz_old
  type(DFT_wavefunction), intent(inout) :: tmb
  type(local_zone_descriptors), intent(in) :: lzd_old
  type(fragment_transformation), dimension(tmb%orbs%norbp), intent(in) :: frag_trans
  real(wp), dimension(:), pointer :: psi_old
  type(phi_array), dimension(tmb%orbs%norbp), optional, intent(in) :: phi_array_old
  logical, intent(in) :: add_derivatives
  !Local variables
  character(len=*), parameter :: subname='reformatmywaves'
  logical :: reformat
  integer :: iorb,j,i_stat,i_all,jstart,jstart_old,iiorb,ilr,iiat
  integer:: idir,jstart_old_der,ncount,ilr_old
  integer, dimension(3) :: ns_old,ns,n_old,n
  real(gp), dimension(3) :: centre_old_box,centre_new_box,da
  real(gp) :: tt,tol
  real(wp), dimension(:,:,:,:,:,:), pointer :: phigold
  real(wp), dimension(:), allocatable :: phi_old_der
  integer, dimension(0:7) :: reformat_reason
!  real(gp) :: dnrm2
!  integer :: iat

  reformat_reason=0
  tol=1.d-3

  ! Get the derivatives of the support functions
  if (add_derivatives) then
     allocate(phi_old_der(3*ndim_old),stat=i_stat)
     call memocc(i_stat,phi_old_der,'phi_old_der',subname)
     if (.not. associated(psi_old)) stop 'psi_old not associated in reformat_supportfunctions'
     call get_derivative_supportfunctions(ndim_old, lzd_old%hgrids(1), lzd_old, tmb%orbs, psi_old, phi_old_der)
     jstart_old_der=1
  end if

  jstart_old=1
  jstart=1
  do iorb=1,tmb%orbs%norbp
      iiorb=tmb%orbs%isorb+iorb
      ilr=tmb%orbs%inwhichlocreg(iiorb)
      iiat=tmb%orbs%onwhichatom(iiorb)

      ilr_old=ilr

      n_old(1)=lzd_old%Llr(ilr_old)%d%n1
      n_old(2)=lzd_old%Llr(ilr_old)%d%n2
      n_old(3)=lzd_old%Llr(ilr_old)%d%n3
      n(1)=tmb%lzd%Llr(ilr)%d%n1
      n(2)=tmb%lzd%Llr(ilr)%d%n2
      n(3)=tmb%lzd%Llr(ilr)%d%n3
      ns_old(1)=lzd_old%Llr(ilr_old)%ns1
      ns_old(2)=lzd_old%Llr(ilr_old)%ns2
      ns_old(3)=lzd_old%Llr(ilr_old)%ns3
      ns(1)=tmb%lzd%Llr(ilr)%ns1
      ns(2)=tmb%lzd%Llr(ilr)%ns2
      ns(3)=tmb%lzd%Llr(ilr)%ns3

      !theta=frag_trans(iorb)%theta!0.0d0*(4.0_gp*atan(1.d0)/180.0_gp)
      !newz=frag_trans(iorb)%rot_axis!(/1.0_gp,0.0_gp,0.0_gp/)
      !centre_old(:)=frag_trans(iorb)%rot_center(:)!rxyz_old(:,iiat)
      !shift(:)=frag_trans(iorb)%dr(:)!rxyz(:,iiat)
 
      call reformat_check(reformat,reformat_reason,tol,at,lzd_old%hgrids,tmb%lzd%hgrids,&
           lzd_old%llr(ilr_old)%wfd%nvctr_c,lzd_old%llr(ilr_old)%wfd%nvctr_f,&
           tmb%lzd%llr(ilr)%wfd%nvctr_c,tmb%lzd%llr(ilr)%wfd%nvctr_f,&
           n_old,n,ns_old,ns,frag_trans(iorb),centre_old_box,centre_new_box,da)  
   
      ! just copy psi from old to new as reformat not necessary
      if (.not. reformat) then 

          ! copy from phi_array_old, can use new keys as they should be identical to old keys
          if (present(phi_array_old)) then 
             call compress_plain(n(1),n(2),0,n(1),0,n(2),0,n(3), &
                  tmb%lzd%llr(ilr)%wfd%nseg_c,tmb%lzd%llr(ilr)%wfd%nvctr_c,tmb%lzd%llr(ilr)%wfd%keygloc(1,1), &
                  tmb%lzd%llr(ilr)%wfd%keyvloc(1),tmb%lzd%llr(ilr)%wfd%nseg_f,tmb%lzd%llr(ilr)%wfd%nvctr_f,&
                  tmb%lzd%llr(ilr)%wfd%keygloc(1,tmb%lzd%llr(ilr)%wfd%nseg_c+min(1,tmb%lzd%llr(ilr)%wfd%nseg_f)),&
                  tmb%lzd%llr(ilr)%wfd%keyvloc(tmb%lzd%llr(ilr)%wfd%nseg_c+min(1,tmb%lzd%llr(ilr)%wfd%nseg_f)),   &
                  phi_array_old(iorb)%psig,tmb%psi(jstart),&
                  tmb%psi(jstart+tmb%lzd%llr(ilr)%wfd%nvctr_c+min(1,tmb%lzd%llr(ilr)%wfd%nvctr_f)-1))
             jstart=jstart+tmb%lzd%llr(ilr)%wfd%nvctr_c+7*tmb%lzd%llr(ilr)%wfd%nvctr_f

          ! directly copy psi_old to psi, first check psi_old is actually allocated
          else 
             if (.not. associated(psi_old)) stop 'psi_old not associated in reformat_supportfunctions'
             do j=1,lzd_old%llr(ilr_old)%wfd%nvctr_c
                tmb%psi(jstart)=psi_old(jstart_old)
                jstart=jstart+1
                jstart_old=jstart_old+1
             end do
             do j=1,7*lzd_old%llr(ilr_old)%wfd%nvctr_f-6,7
                tmb%psi(jstart+0)=psi_old(jstart_old+0)
                tmb%psi(jstart+1)=psi_old(jstart_old+1)
                tmb%psi(jstart+2)=psi_old(jstart_old+2)
                tmb%psi(jstart+3)=psi_old(jstart_old+3)
                tmb%psi(jstart+4)=psi_old(jstart_old+4)
                tmb%psi(jstart+5)=psi_old(jstart_old+5)
                tmb%psi(jstart+6)=psi_old(jstart_old+6)
                jstart=jstart+7
                jstart_old=jstart_old+7
            end do
         
         end if
      else
          ! Add the derivatives to the basis functions
          if (add_derivatives) then
             do idir=1,3
                 tt=rxyz(idir,iiat)-rxyz_old(idir,iiat)
                 ncount = lzd_old%llr(ilr_old)%wfd%nvctr_c+7*lzd_old%llr(ilr_old)%wfd%nvctr_f
                 call daxpy(ncount, tt, phi_old_der(jstart_old_der), 1, psi_old(jstart_old), 1)
                 jstart_old_der = jstart_old_der + ncount
             end do
          end if

          ! uncompress or point towards correct phigold as necessary
          if (present(phi_array_old)) then
             phigold=>phi_array_old(iorb)%psig
          else
             allocate(phigold(0:n_old(1),2,0:n_old(2),2,0:n_old(3),2+ndebug),stat=i_stat)
             call memocc(i_stat,phigold,'phigold',subname)
             call psi_to_psig(n_old,lzd_old%llr(ilr_old)%wfd%nvctr_c,lzd_old%llr(ilr_old)%wfd%nvctr_f,&
                  lzd_old%llr(ilr_old)%wfd%nseg_c,lzd_old%llr(ilr_old)%wfd%nseg_f,&
                  lzd_old%llr(ilr_old)%wfd%keyvloc,lzd_old%llr(ilr_old)%wfd%keygloc,&
                  jstart_old,psi_old(jstart_old),phigold)
          end if
   
          !write(100+iproc,*) 'norm phigold ',dnrm2(8*(n1_old+1)*(n2_old+1)*(n3_old+1),phigold,1)
          !write(*,*) 'iproc,norm phigold ',iproc,dnrm2(8*(n1_old+1)*(n2_old+1)*(n3_old+1),phigold,1)

          call reformat_one_supportfunction(tmb%lzd%llr(ilr)%wfd,tmb%lzd%llr(ilr)%geocode,lzd_old%hgrids,&
               n_old,phigold,tmb%lzd%hgrids,n,centre_old_box,centre_new_box,da,&
               frag_trans(iorb),tmb%psi(jstart))

          jstart=jstart+tmb%lzd%llr(ilr)%wfd%nvctr_c+7*tmb%lzd%llr(ilr)%wfd%nvctr_f

          if (present(phi_array_old)) then   
             nullify(phigold)
          else
             i_all=-product(shape(phigold))*kind(phigold)
             deallocate(phigold,stat=i_stat)
             call memocc(i_stat,i_all,'phigold',subname)
          end if

      end if

  end do

  if (add_derivatives) then
     i_all=-product(shape(phi_old_der))*kind(phi_old_der)
     deallocate(phi_old_der,stat=i_stat)
     call memocc(i_stat,i_all,'phi_old_der',subname)
  end if

  call print_reformat_summary(iproc,reformat_reason)

END SUBROUTINE reformat_supportfunctions


!checks whether reformatting is needed based on various criteria and returns final shift and centres needed for reformat
subroutine reformat_check(reformat_needed,reformat_reason,tol,at,hgrids_old,hgrids,nvctr_c_old,nvctr_f_old,&
       nvctr_c,nvctr_f,n_old,n,ns_old,ns,frag_trans,centre_old_box,centre_new_box,da)  
  use module_base
  use module_types
  use module_fragments
  implicit none

  logical, intent(out) :: reformat_needed ! logical telling whether reformat is needed
  integer, dimension(0:7), intent(inout) :: reformat_reason ! array giving reasons for reformatting
  real(gp), intent(in) :: tol ! tolerance for rotations and shifts
  type(atoms_data), intent(in) :: at
  real(gp), dimension(3), intent(in) :: hgrids, hgrids_old
  integer, intent(in) :: nvctr_c, nvctr_f, nvctr_c_old, nvctr_f_old
  integer, dimension(3), intent(in) :: n, n_old, ns, ns_old
  real(gp), dimension(3), intent(out) :: centre_old_box, centre_new_box ! centres of rotation wrt box
  real(gp), dimension(3), intent(out) :: da ! shift to be used in reformat
  type(fragment_transformation), intent(in) :: frag_trans ! includes centres of rotation in global coordinates, shift and angle

  ! local variables 
  real(gp) :: displ, mindist
  integer, dimension(3) :: nb
  logical, dimension(3) :: per

  !conditions for periodicity in the three directions
  per(1)=(at%astruct%geocode /= 'F')
  per(2)=(at%astruct%geocode == 'P')
  per(3)=(at%astruct%geocode /= 'F')

  !buffers related to periodicity
  !WARNING: the boundary conditions are not assumed to change between new and old
  call ext_buffers_coarse(per(1),nb(1))
  call ext_buffers_coarse(per(2),nb(2))
  call ext_buffers_coarse(per(3),nb(3))

  ! centre of rotation with respect to start of box
  centre_old_box(1)=mindist(per(1),at%astruct%cell_dim(1),frag_trans%rot_center(1),hgrids_old(1)*(ns_old(1)-0.5_dp*nb(1)))
  centre_old_box(2)=mindist(per(2),at%astruct%cell_dim(2),frag_trans%rot_center(2),hgrids_old(2)*(ns_old(2)-0.5_dp*nb(2)))
  centre_old_box(3)=mindist(per(3),at%astruct%cell_dim(3),frag_trans%rot_center(3),hgrids_old(3)*(ns_old(3)-0.5_dp*nb(3)))

  centre_new_box(1)=mindist(per(1),at%astruct%cell_dim(1),frag_trans%rot_center_new(1),hgrids(1)*(ns(1)-0.5_dp*nb(1)))
  centre_new_box(2)=mindist(per(2),at%astruct%cell_dim(2),frag_trans%rot_center_new(2),hgrids(2)*(ns(2)-0.5_dp*nb(2)))
  centre_new_box(3)=mindist(per(3),at%astruct%cell_dim(3),frag_trans%rot_center_new(3),hgrids(3)*(ns(3)-0.5_dp*nb(3)))

  !Calculate the shift of the atom to be used in reformat
  da(1)=mindist(per(1),at%astruct%cell_dim(1),centre_new_box(1),centre_old_box(1))
  da(2)=mindist(per(2),at%astruct%cell_dim(2),centre_new_box(2),centre_old_box(2))
  da(3)=mindist(per(3),at%astruct%cell_dim(3),centre_new_box(3),centre_old_box(3))


  !print*,'reformat check',frag_trans%rot_center(2),ns_old(2),centre_old_box(2),&
  !     frag_trans%rot_center_new(2),ns(2),centre_new_box(2),da(2)

  !!write(*,'(a,3(3(f6.3,x),3x))') 'final',centre_old_box,centre_new_box,da

  displ=sqrt(da(1)**2+da(2)**2+da(3)**2)

  !reformatting criterion
  if (hgrids(1) == hgrids_old(1) .and. hgrids(2) == hgrids_old(2) .and. hgrids(3) == hgrids_old(3) &
        .and. nvctr_c  == nvctr_c_old .and. nvctr_f  == nvctr_f_old &
        .and. n_old(1)==n(1)  .and. n_old(2)==n(2) .and. n_old(3)==n(3) &
        .and. abs(frag_trans%theta) <= tol .and. abs(displ) <= tol .and. size(frag_trans%discrete_operations)==0) then
      reformat_reason(0) = reformat_reason(0) + 1
      reformat_needed=.false.
  else
      reformat_needed=.true.
      if (hgrids(1) /= hgrids_old(1) .or. hgrids(2) /= hgrids_old(2) .or. hgrids(3) /= hgrids_old(3)) then 
         reformat_reason(1) = reformat_reason(1) + 1
      end if
      if (nvctr_c  /= nvctr_c_old) then
         reformat_reason(2) = reformat_reason(2) + 1
      end if
      if (nvctr_f  /= nvctr_f_old) then
         reformat_reason(3) = reformat_reason(3) + 1
      end if
      if (n_old(1) /= n(1)  .or. n_old(2) /= n(2) .or. n_old(3) /= n(3) )  then  
         reformat_reason(4) = reformat_reason(4) + 1
      end if
      if (abs(displ) > tol)  then  
         reformat_reason(5) = reformat_reason(5) + 1
      end if
      if (abs(frag_trans%theta) > tol)  then  
         reformat_reason(6) = reformat_reason(6) + 1
      end if
      if (size(frag_trans%discrete_operations) > 0)  then  
         reformat_reason(7) = reformat_reason(7) + 1
      end if
  end if

end subroutine reformat_check

subroutine print_reformat_summary(iproc,reformat_reason)
  use module_base
  use module_types
  implicit none

  integer, intent(in) :: iproc
  integer, dimension(0:7), intent(inout) :: reformat_reason ! array giving reasons for reformatting

  integer :: ierr

  call mpiallred(reformat_reason(0), 7, mpi_sum, bigdft_mpi%mpi_comm, ierr)
  if (iproc==0) then
        write(*,'(1x,a)') 'Overview of the reformatting (several categories may apply):'
        write(*,'(3x,a,i0)') '- No reformatting required: ', reformat_reason(0)
        write(*,'(3x,a,i0)') '- Grid spacing has changed: ', reformat_reason(1)
        write(*,'(3x,a,i0)') '- Number of coarse grid points has changed: ', reformat_reason(2)
        write(*,'(3x,a,i0)') '- Number of fine grid points has changed: ', reformat_reason(3)
        write(*,'(3x,a,i0)') '- Box size has changed: ', reformat_reason(4)
        write(*,'(3x,a,i0)') '- Molecule was shifted: ', reformat_reason(5)
        write(*,'(3x,a,i0)') '- Molecule was rotated: ', reformat_reason(6)
        write(*,'(3x,a,i0)') '- Discrete operations: ', reformat_reason(7)
  end if

end subroutine print_reformat_summary

subroutine psi_to_psig(n,nvctr_c,nvctr_f,nseg_c,nseg_f,keyvloc,keygloc,jstart,psi,psig)
  use module_base
  implicit none

  integer, dimension(3), intent(in) :: n
  integer, intent(in) :: nseg_c, nseg_f, nvctr_c, nvctr_f
  integer, dimension(nseg_c+nseg_f), intent(in) :: keyvloc
  integer, dimension(2,nseg_c+nseg_f), intent(in) :: keygloc
  integer, intent(inout) :: jstart
  real(wp), dimension(jstart:jstart+nvctr_c+7*nvctr_f), intent(in) :: psi
  real(wp), dimension(0:n(1),2,0:n(2),2,0:n(3),2), intent(out) :: psig

  ! local variables
  integer :: iseg, jj, j0, j1, i, ii, i0, i1, i2, i3

  call razero(8*(n(1)+1)*(n(2)+1)*(n(3)+1),psig(0,1,0,1,0,1))

  ! coarse part
  do iseg=1,nseg_c
     jj=keyvloc(iseg)
     j0=keygloc(1,iseg)
     j1=keygloc(2,iseg)
     ii=j0-1
     i3=ii/((n(1)+1)*(n(2)+1))
     ii=ii-i3*(n(1)+1)*(n(2)+1)
     i2=ii/(n(1)+1)
     i0=ii-i2*(n(1)+1)
     i1=i0+j1-j0
     do i=i0,i1
        psig(i,1,i2,1,i3,1) = psi(jstart)
        jstart=jstart+1
     end do
  end do
   
  ! fine part
  do iseg=1,nseg_f
     jj=keyvloc(nseg_c + iseg)
     j0=keygloc(1,nseg_c + iseg)
     j1=keygloc(2,nseg_c + iseg)
     ii=j0-1
     i3=ii/((n(1)+1)*(n(2)+1))
     ii=ii-i3*(n(1)+1)*(n(2)+1)
     i2=ii/(n(1)+1)
     i0=ii-i2*(n(1)+1)
     i1=i0+j1-j0
     do i=i0,i1
        psig(i,2,i2,1,i3,1)=psi(jstart+0)
        psig(i,1,i2,2,i3,1)=psi(jstart+1)
        psig(i,2,i2,2,i3,1)=psi(jstart+2)
        psig(i,1,i2,1,i3,2)=psi(jstart+3)
        psig(i,2,i2,1,i3,2)=psi(jstart+4)
        psig(i,1,i2,2,i3,2)=psi(jstart+5)
        psig(i,2,i2,2,i3,2)=psi(jstart+6)
        jstart=jstart+7
     end do
  end do

end subroutine psi_to_psig



