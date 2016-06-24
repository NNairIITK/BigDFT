!> @file
!!  Routines to do restart of the calculation
!! @author
!!    Copyright (C) 2007-2013 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS


!> Copy old wavefunctions from psi to psi_old
subroutine copy_old_wavefunctions(nproc,orbs,psi,&
     wfd_old,psi_old)
  use module_base
  use module_types
  use yaml_output
  implicit none
  integer, intent(in) :: nproc
  type(orbitals_data), intent(in) :: orbs
  type(wavefunctions_descriptors), intent(in) :: wfd_old
  real(wp), dimension(:), pointer :: psi,psi_old
  !Local variables
  character(len=*), parameter :: subname='copy_old_wavefunctions'
  !real(kind=8), parameter :: eps_mach=1.d-12
  integer :: j,ind1,iorb,oidx,sidx !n(c) nvctrp_old
  real(kind=8) :: tt
  call f_routine(id=subname)

  !add the number of distributed point for the compressed wavefunction
  !tt=dble(wfd_old%nvctr_c+7*wfd_old%nvctr_f)/dble(nproc)
  !n(c) nvctrp_old=int((1.d0-eps_mach*tt) + tt)

!  psi_old=&
!       f_malloc_ptr((wfd_old%nvctr_c+7*wfd_old%nvc_f)*orbs%norbp*orbs%nspinor,!&
!       id='psi_old')
  psi_old = f_malloc_ptr((wfd_old%nvctr_c+7*wfd_old%nvctr_f)*orbs%norbp*orbs%nspinor,id='psi_old')

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
  call f_free_ptr(psi)

  call f_release_routine()

END SUBROUTINE copy_old_wavefunctions


!> Reformat wavefunctions if the mesh have changed (in a restart)
subroutine reformatmywaves(iproc,orbs,at,&
     hx_old,hy_old,hz_old,n1_old,n2_old,n3_old,rxyz_old,wfd_old,psi_old,&
     hx,hy,hz,n1,n2,n3,rxyz,wfd,psi)
  use module_base
  use module_types
  use yaml_output
  use box
  use bounds, only: ext_buffers_coarse
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
  integer :: iat,iorb,j,jj,j0,j1,ii,i0,i1,i2,i3,i,iseg,nb1,nb2,nb3,nvctrcj,n1p1,np,i0jj
  real(gp) :: tx,ty,tz,displ
  type(cell) :: mesh
  real(wp), dimension(:,:,:), allocatable :: psifscf
  real(wp), dimension(:,:,:,:,:,:), allocatable :: psigold


  mesh=cell_new(at%astruct%geocode,[n1,n2,n3],[hx,hy,hz])
  !conditions for periodicity in the three directions
  perx=(at%astruct%geocode /= 'F')
  pery=(at%astruct%geocode == 'P')
  perz=(at%astruct%geocode /= 'F')

  !buffers related to periodicity
  !WARNING: the boundary conditions are not assumed to change between new and old
  call ext_buffers_coarse(perx,nb1)
  call ext_buffers_coarse(pery,nb2)
  call ext_buffers_coarse(perz,nb3)


  psifscf = f_malloc((/ -nb1.to.2*n1+1+nb1, -nb2.to.2*n2+1+nb2, -nb3.to.2*n3+1+nb3 /),id='psifscf')

  tx=0.0_gp
  ty=0.0_gp
  tz=0.0_gp
displ=0.0_gp
  do iat=1,at%astruct%nat
    displ=displ+minimum_distance(mesh,rxyz(:,iat),rxyz_old(:,iat))**2
    !  tx=tx+mindist(perx,at%astruct%cell_dim(1),rxyz(1,iat),rxyz_old(1,iat))**2
    !  ty=ty+mindist(pery,at%astruct%cell_dim(2),rxyz(2,iat),rxyz_old(2,iat))**2
    !  tz=tz+mindist(perz,at%astruct%cell_dim(3),rxyz(3,iat),rxyz_old(3,iat))**2
  enddo
  displ=sqrt(displ)!tx+ty+tz)
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
        call yaml_mapping_open('Reformatting for')
        !write(*,'(1x,a)') 'The wavefunctions need reformatting because:'
        if (hx /= hx_old .or. hy /= hy_old .or. hz /= hz_old) then
           call yaml_mapping_open('hgrid modified',flow=.true.)
              call yaml_map('hgrid_old', (/ hx_old,hy_old,hz_old /),fmt='(1pe20.12)')
              call yaml_map('hgrid', (/ hx,hy,hz /), fmt='(1pe20.12)')
           call yaml_mapping_close()
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
           call yaml_map('Molecule was shifted, norm' , displ , fmt='(1pe19.12)')
           !write(*,"(4x,a,3(1pe19.12))") 'molecule was shifted  ' , tx,ty,tz
        endif
        !write(*,"(1x,a)",advance='NO') 'Reformatting...'
        call yaml_mapping_close()
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
           nvctrcj=wfd%nvctr_c+j
           psi(nvctrcj+0,iorb)=psi_old(nvctrcj+0,iorb)
           psi(nvctrcj+1,iorb)=psi_old(nvctrcj+1,iorb)
           psi(nvctrcj+2,iorb)=psi_old(nvctrcj+2,iorb)
           psi(nvctrcj+3,iorb)=psi_old(nvctrcj+3,iorb)
           psi(nvctrcj+4,iorb)=psi_old(nvctrcj+4,iorb)
           psi(nvctrcj+5,iorb)=psi_old(nvctrcj+5,iorb)
           psi(nvctrcj+6,iorb)=psi_old(nvctrcj+6,iorb)
        enddo

     else

        psigold = f_malloc0((/ 0.to.n1_old, 1.to.2, 0.to.n2_old, 1.to.2, 0.to.n3_old, 1.to.2 /),id='psigold')

        !call f_zero(8*(n1_old+1)*(n2_old+1)*(n3_old+1),psigold)

        n1p1=n1_old+1
        np=n1p1*(n2_old+1)
        ! coarse part
        do iseg=1,wfd_old%nseg_c
           jj=wfd_old%keyvloc(iseg)
           j0=wfd_old%keygloc(1,iseg)
           j1=wfd_old%keygloc(2,iseg)
           ii=j0-1
           i3=ii/np
           ii=ii-i3*np
           i2=ii/n1p1
           i0=ii-i2*n1p1
           i1=i0+j1-j0
           i0jj=jj-i0
           do i=i0,i1
              psigold(i,1,i2,1,i3,1) = psi_old(i+i0jj,iorb)
           enddo
        enddo

        ! fine part
        do iseg=1,wfd_old%nseg_f
           jj=wfd_old%keyvloc(wfd_old%nseg_c + iseg)
           j0=wfd_old%keygloc(1,wfd_old%nseg_c + iseg)
           j1=wfd_old%keygloc(2,wfd_old%nseg_c + iseg)
           ii=j0-1
           i3=ii/np
           ii=ii-i3*np
           i2=ii/n1p1
           i0=ii-i2*n1p1
           i1=i0+j1-j0
           i0jj=jj-i0-1
           do i=i0,i1
              psigold(i,2,i2,1,i3,1)=psi_old(wfd_old%nvctr_c+1+7*(i+i0jj), iorb)
              psigold(i,1,i2,2,i3,1)=psi_old(wfd_old%nvctr_c+2+7*(i+i0jj), iorb)
              psigold(i,2,i2,2,i3,1)=psi_old(wfd_old%nvctr_c+3+7*(i+i0jj), iorb)
              psigold(i,1,i2,1,i3,2)=psi_old(wfd_old%nvctr_c+4+7*(i+i0jj), iorb)
              psigold(i,2,i2,1,i3,2)=psi_old(wfd_old%nvctr_c+5+7*(i+i0jj), iorb)
              psigold(i,1,i2,2,i3,2)=psi_old(wfd_old%nvctr_c+6+7*(i+i0jj), iorb)
              psigold(i,2,i2,2,i3,2)=psi_old(wfd_old%nvctr_c+7+7*(i+i0jj), iorb)
           enddo
        enddo

!write(100+iproc,*) 'norm psigold ',dnrm2(8*(n1_old+1)*(n2_old+1)*(n3_old+1),psigold,1)

        call reformatonewave(displ,wfd,at,hx_old,hy_old,hz_old, & !n(m)
             n1_old,n2_old,n3_old,rxyz_old,psigold,hx,hy,hz,&
             n1,n2,n3,rxyz,psifscf,psi(1,iorb))

        call f_free(psigold)
     end if
  end do

  call f_free(psifscf)

  !if (iproc==0) write(*,"(1x,a)")'done.'

END SUBROUTINE reformatmywaves




!> Reads wavefunction from file and transforms it properly if hgrid or size of simulation cell
!!  have changed
subroutine readmywaves(iproc,filename,iformat,orbs,n1,n2,n3,hx,hy,hz,at,rxyz_old,rxyz,  &
     wfd,psi,orblist)
  use module_base
  use module_types
  use yaml_output
  use module_interfaces, only: open_filename_of_iorb
  use public_enums
  use bounds, only: ext_buffers_coarse
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
  integer :: ncount1,ncount_rate,ncount_max,iorb,ncount2,nb1,nb2,nb3,iorb_out,ispinor,unitwf
  real(kind=4) :: tr0,tr1
  real(kind=8) :: tel
  real(wp), dimension(:,:,:), allocatable :: psifscf
  !integer, dimension(orbs%norb) :: orblist2

  call cpu_time(tr0)
  call system_clock(ncount1,ncount_rate,ncount_max)

  unitwf=99

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

     psifscf = f_malloc((/ -nb1.to.2*n1+1+nb1, -nb2.to.2*n2+1+nb2, -nb3.to.2*n3+1+nb3 /),id='psifscf')

     do iorb=1,orbs%norbp!*orbs%nspinor

        do ispinor=1,orbs%nspinor
           if(present(orblist)) then
              call open_filename_of_iorb(unitwf,(iformat == WF_FORMAT_BINARY),filename, &
                   & orbs,iorb,ispinor,iorb_out, orblist(iorb+orbs%isorb))
           else
              call open_filename_of_iorb(unitwf,(iformat == WF_FORMAT_BINARY),filename, &
                   & orbs,iorb,ispinor,iorb_out)
           end if
           call readonewave(unitwf, (iformat == WF_FORMAT_PLAIN),iorb_out,iproc,n1,n2,n3, &
                & hx,hy,hz,at,wfd,rxyz_old,rxyz,&
                psi(1,ispinor,iorb),orbs%eval(orbs%isorb+iorb),psifscf)
           call f_close(unitwf)
        end do
!!$        do i_all=1,wfd%nvctr_c+7*wfd%nvctr_f
!!$            write(700+iorb,*) i_all, psi(i_all,1,iorb)
!!$        end do

     end do

     call f_free(psifscf)

  else
     call yaml_warning('Unknown wavefunction file format from filename.')
     stop
  end if

  call cpu_time(tr1)
  call system_clock(ncount2,ncount_rate,ncount_max)
  tel=dble(ncount2-ncount1)/dble(ncount_rate)


  if (iproc == 0) then
     call yaml_sequence_open('Reading Waves Time')
     call yaml_sequence(advance='no')
     call yaml_mapping_open(flow=.true.)
     call yaml_map('Process',iproc)
     call yaml_map('Timing',(/ real(tr1-tr0,kind=8),tel /),fmt='(1pe10.3)')
     call yaml_mapping_close()
     call yaml_sequence_close()
  end if
  !write(*,'(a,i4,2(1x,1pe10.3))') '- READING WAVES TIME',iproc,tr1-tr0,tel
END SUBROUTINE readmywaves


!> Verify the presence of a given file
subroutine verify_file_presence(filerad,orbs,iformat,nproc,nforb)
  use module_base
  use module_types
  use public_enums
  use module_interfaces, only: filename_of_iorb
  implicit none
  integer, intent(in) :: nproc
  character(len=*), intent(in) :: filerad
  type(orbitals_data), intent(in) :: orbs
  integer, intent(out) :: iformat
  integer, optional, intent(in) :: nforb
  !local variables
  character(len=500) :: filename
  logical :: onefile,allfiles
  integer :: iorb,ispinor,iorb_out

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
  if (nproc > 1) call mpiallred(allfiles,1,MPI_LAND,comm=bigdft_mpi%mpi_comm)

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
  if (nproc > 1) call mpiallred(allfiles,1,MPI_LAND,comm=bigdft_mpi%mpi_comm)

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
  !Arguments
  character(len=*), intent(in) :: filename
  logical, intent(in) :: lbin
  integer, intent(in) :: iorb,ispinor
  type(orbitals_data), intent(in) :: orbs
  character(len=*), intent(out) :: filename_out
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

!!$  !ternary example, interesting syntax
!!$  realimag = .if. modulo(ispinor,2)==0 .then. 'I' .else. 'R'

  !see if the wavefunction is real or imaginary
  realimag=merge('I','R',modulo(ispinor,2)==0)

!!$  if(modulo(ispinor,2)==0) then
!!$     realimag='I'
!!$  else
!!$     realimag='R'
!!$  end if
  !calculate the spin sector
  spins=orbs%spinsgn(orbs%isorb+iorb)
!!$  spintype=.if. (orbs%nspinor == 4) .then. merge('A','B',ispinor <=2) &
!!$       .else. merge('U','D',spins==1.0_gp)
  if(orbs%nspinor == 4) then
     spintype=merge('A','B',ispinor <=2)
!!$     if (ispinor <=2) then
!!$        spintype='A'
!!$     else
!!$        spintype='B'
!!$     end if
  else
     spintype=merge('U','D',spins==1.0_gp)
!!$     if (spins==1.0_gp) then
!!$        spintype='U'
!!$     else
!!$        spintype='D'
!!$     end if
  end if
  !no spin polarization if nspin=1
  if (orbs%nspin==1) spintype='N'

  !calculate the actual orbital value
  iorb_out=iorb+orbs%isorb-(ikpt-1)*orbs%norb
  !print *,"iorb_out (filename_of_iorb) = ", iorb_out
  !print *,"ikpt (filename_of_iorb) = ", ikpt
  !print *,"orbs%isorb (filename_of_iorb) = ", orbs%isorb

  if(present(iiorb)) iorb_out = iiorb
  !purge the value from the spin sign
  if (spins==-1.0_gp) iorb_out=iorb_out-orbs%norbu

  !value of the orbital
  write(f4,'(a1,i4.4)') "b", iorb_out

  !complete the information in the name of the orbital
  completename='-'//f3//'-'//spintype//realimag
  if (lbin) then
     filename_out = trim(filename)//completename//".bin."//f4
     !print *,'complete name <',trim(filename_out),'> end'
 else
     filename_out = trim(filename)//completename//"."//f4
     !print *,'complete name <',trim(filename_out),'> end'
 end if

  !print *,'filename: ',filename_out
end subroutine filename_of_iorb

subroutine wfn_filename(filename_out,filename,lbin,ikpt,nspinor,nspin,ispinor,spin,iorb)
  use f_ternary
  use yaml_strings
  use module_defs, only: gp
  implicit none
  logical, intent(in) :: lbin
  integer, intent(in) :: ikpt,nspinor,ispinor,nspin,iorb
  real(gp), intent(in) :: spin
  character(len=*), intent(in) :: filename
  character(len=*), intent(out) :: filename_out
  !local variables
  character(len=1) :: spintype

  spintype=.if. (nspinor == 4) .then. merge('A','B',ispinor <=2) .else. merge('U','D',spin==1.0_gp)
  if (nspin==1) spintype='N'  

  call f_strcpy(dest=filename_out,src=&
       filename+'-k'+ikpt**'(i3.3)'+'-'//spintype//&
       merge('I','R',modulo(ispinor,2)==0)+&
       (.if. lbin .then. '.bin.b' .else. '.b')+iorb**'(i4.4)')
end subroutine wfn_filename
  

!> Associate to the absolute value of orbital a filename which depends of the k-point and
!! of the spin sign
subroutine open_filename_of_iorb(unitfile,lbin,filename,orbs,iorb,ispinor,iorb_out,iiorb)
  use module_base
  use module_types
  use module_interfaces, only: filename_of_iorb
  implicit none
  character(len=*), intent(in) :: filename
  logical, intent(in) :: lbin
  integer, intent(in) :: iorb,ispinor
  type(orbitals_data), intent(in) :: orbs
  !>on entry, it suggests the opening unit. On exit, returns the first valid value to which the unit can be associated
  integer, intent(inout) :: unitfile
  integer, intent(out) :: iorb_out
  integer, intent(in), optional :: iiorb
  !local variables
  character(len=500) :: filename_out

  if(present(iiorb)) then
     call filename_of_iorb(lbin,filename,orbs,iorb,ispinor,filename_out,iorb_out,iiorb)
     !restore previous behaviour even though the wannier construction can be compromised
     !call filename_of_iorb(lbin,filename,orbs,iorb,ispinor,filename_out,iorb_out)
  else
     call filename_of_iorb(lbin,filename,orbs,iorb,ispinor,filename_out,iorb_out)
  end if
  call f_open_file(unitfile,file=filename_out,binary=lbin)

end subroutine open_filename_of_iorb


!> Write all my wavefunctions in files by calling writeonewave
subroutine writemywaves(iproc,filename,iformat,orbs,n1,n2,n3,hx,hy,hz,at,rxyz,wfd,psi)
  use module_types
  use module_base
  use yaml_output
  use module_interfaces, only: open_filename_of_iorb
  use public_enums
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
  integer :: ncount1,ncount_rate,ncount_max,iorb,ncount2,iorb_out,ispinor,unitwf
  real(kind=4) :: tr0,tr1
  real(kind=8) :: tel

  unitwf=99

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
           call open_filename_of_iorb(unitwf,(iformat == WF_FORMAT_BINARY),filename, &
                & orbs,iorb,ispinor,iorb_out)
           call writeonewave(unitwf,(iformat == WF_FORMAT_PLAIN),iorb_out,n1,n2,n3,hx,hy,hz, &
                at%astruct%nat,rxyz,wfd%nseg_c,wfd%nvctr_c,wfd%keygloc(1,1),wfd%keyvloc(1),  &
                wfd%nseg_f,wfd%nvctr_f,wfd%keygloc(1,wfd%nseg_c+1),wfd%keyvloc(wfd%nseg_c+1), &
                psi(1,ispinor,iorb),psi(wfd%nvctr_c+1,ispinor,iorb), &
                orbs%eval(iorb+orbs%isorb))
           call f_close(unitwf)
        end do
     enddo

     call cpu_time(tr1)
     call system_clock(ncount2,ncount_rate,ncount_max)
     tel=dble(ncount2-ncount1)/dble(ncount_rate)
     if (iproc == 0) then
        call yaml_sequence_open('Write Waves Time')
        call yaml_sequence(advance='no')
        call yaml_mapping_open(flow=.true.)
        call yaml_map('Process',iproc)
        call yaml_map('Timing',(/ real(tr1-tr0,kind=8),tel /),fmt='(1pe10.3)')
        call yaml_mapping_close()
        call yaml_sequence_close()
     end if
     !write(*,'(a,i4,2(1x,1pe10.3))') '- WRITE WAVES TIME',iproc,tr1-tr0,tel
     !write(*,'(a,1x,i0,a)') '- iproc',iproc,' finished writing waves'
  end if

END SUBROUTINE writemywaves


subroutine read_wave_to_isf(lstat, filename, ln, iorbp, hx, hy, hz, &
     & n1, n2, n3, nspinor, psiscf)
  use module_base
  use module_types
  use module_interfaces, only: readwavetoisf, readwavetoisf_etsf
  use public_enums
  use module_input_keys
  implicit none

  integer, intent(in) :: ln
  character, intent(in) :: filename(ln)
  integer, intent(in) :: iorbp
  integer, intent(out) :: n1, n2, n3, nspinor
  real(gp), intent(out) :: hx, hy, hz
  real(wp), dimension(:,:,:,:), pointer :: psiscf
  logical, intent(out) :: lstat

  character(len = 1024) :: filename_
  integer :: iformat, i

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


!!$subroutine free_wave_to_isf(psiscf)
!!$  use module_base
!!$  implicit none
!!$  real(wp), dimension(:,:,:,:), pointer :: psiscf
!!$
!!$  integer :: i_all, i_stat
!!$
!!$  i_all=-product(shape(psiscf))*kind(psiscf)
!!$  deallocate(psiscf,stat=i_stat)
!!$  call memocc(i_stat,i_all,'psiscf',"free_wave_to_isf_etsf")
!!$END SUBROUTINE free_wave_to_isf


subroutine read_wave_descr(lstat, filename, ln, &
     & norbu, norbd, iorbs, ispins, nkpt, ikpts, nspinor, ispinor)
  use public_enums !module_types
  use module_input_keys
  implicit none
  integer, intent(in) :: ln
  character, intent(in) :: filename(ln)
  integer, intent(out) :: norbu, norbd, nkpt, nspinor
  integer, intent(out) :: iorbs, ispins, ikpts, ispinor
  logical, intent(out) :: lstat

  character(len = 1024) :: filename_
  integer :: iformat, i
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

subroutine tmb_overlap_onsite(iproc, nproc, imethod_overlap, at, tmb, rxyz)

  use module_base
  use module_types
  use locregs, only: copy_locreg_descriptors,allocate_wfd,deallocate_wfd
  use module_interfaces, only: reformat_one_supportfunction
  use module_fragments
  use communications_base, only: comms_linear_null, deallocate_comms_linear, TRANSPOSE_FULL
  use communications_init, only: init_comms_linear
  use communications, only: transpose_localized
  use sparsematrix, only: uncompress_matrix
  use sparsematrix_base, only: matrices, sparse_matrix, &
                               matrices_null, sparse_matrix_null, &
                               deallocate_matrices, deallocate_sparse_matrix, &
                               assignment(=), sparsematrix_malloc_ptr, SPARSE_TASKGROUP
  use sparsematrix_wrappers, only: init_sparse_matrix_wrapper
  use sparsematrix_init, only: init_matrix_taskgroups
  use bigdft_matrices, only: check_local_matrix_extents, init_matrixindex_in_compressed_fortransposed
  use transposed_operations, only: calculate_overlap_transposed, normalize_transposed
  !!use bounds, only: ext_buffers
  !!use locreg_operations
  implicit none

  ! Calling arguments
  integer,intent(in) :: iproc, nproc, imethod_overlap
  type(atoms_data), intent(in) :: at
  type(DFT_wavefunction),intent(inout):: tmb
  real(gp),dimension(3,at%astruct%nat),intent(in) :: rxyz

  ! Local variables
  logical :: reformat
  integer :: iorb,jstart,jstart_tmp
  integer :: iiorb,ilr,iiat,j,iis1,iie1,i1,i
  integer :: ilr_tmp,iiat_tmp,ndim_tmp,ndim,norb_tmp
  integer, dimension(3) :: ns,ns_tmp,n,n_tmp,nglr
  real(gp), dimension(3) :: centre_old_box, centre_new_box, da
  real(wp), dimension(:,:,:,:,:,:), allocatable :: phigold
  real(wp), dimension(:), pointer :: psi_tmp, psit_c_tmp, psit_f_tmp, norm
  integer, dimension(0:7) :: reformat_reason
  type(comms_linear) :: collcom_tmp
  type(local_zone_descriptors) :: lzd_tmp
  real(gp) :: tol
  character(len=*),parameter:: subname='tmb_overlap_onsite'
  type(fragment_transformation) :: frag_trans
  integer :: ierr, ncount, iroot, jproc, ndim_tmp1
  integer,dimension(:),allocatable :: workarray
  type(sparse_matrix) :: smat_tmp
  type(matrices) :: mat_tmp
  integer,dimension(2) :: irow, icol, iirow, iicol
  logical :: wrap_around
  real(gp) :: ddot
  integer :: ind_min, ind_mas, ind_trans_min, ind_trans_max
  !!real(kind=gp), dimension(:,:,:), allocatable :: workarraytmp
  !!real(wp), allocatable, dimension(:,:,:) :: psirold
  !!integer, dimension(3) :: nl, nr
  !!logical, dimension(3) :: per
  !!type(workarr_sumrho) :: w

  ! move all psi into psi_tmp all centred in the same place and calculate overlap matrix
  tol=1.d-3
  reformat_reason=0

  !arbitrarily pick the middle one as assuming it'll be near the centre of structure
  !and therefore have large fine grid
  norb_tmp=tmb%orbs%norb/2
  ilr_tmp=tmb%orbs%inwhichlocreg(norb_tmp)
  iiat_tmp=tmb%orbs%onwhichatom(norb_tmp)

  ! Find out which process handles TMB norb_tmp and get the keys from that process
  do jproc=0,nproc-1
      if (tmb%orbs%isorb_par(jproc)<norb_tmp .and. norb_tmp<=tmb%orbs%isorb_par(jproc)+tmb%orbs%norb_par(jproc,0)) then
          iroot=jproc
          exit
      end if
  end do
  if (iproc/=iroot) then
      ! some processes might already have it allocated
      call deallocate_wfd(tmb%lzd%llr(ilr_tmp)%wfd)
      call allocate_wfd(tmb%lzd%llr(ilr_tmp)%wfd)
  end if
  if (nproc>1) then
      ncount = tmb%lzd%llr(ilr_tmp)%wfd%nseg_c + tmb%lzd%llr(ilr_tmp)%wfd%nseg_f
      workarray = f_malloc(6*ncount,id='workarray')
      if (iproc==iroot) then
          call vcopy(2*ncount, tmb%lzd%llr(ilr_tmp)%wfd%keygloc(1,1), 1, workarray(1), 1)
          call vcopy(2*ncount, tmb%lzd%llr(ilr_tmp)%wfd%keyglob(1,1), 1, workarray(2*ncount+1), 1)
          call vcopy(ncount, tmb%lzd%llr(ilr_tmp)%wfd%keyvloc(1), 1, workarray(4*ncount+1), 1)
          call vcopy(ncount, tmb%lzd%llr(ilr_tmp)%wfd%keyvglob(1), 1, workarray(5*ncount+1), 1)
      end if
      call mpibcast(workarray,root=iroot, comm=bigdft_mpi%mpi_comm)
      if (iproc/=iroot) then
          call vcopy(2*ncount, workarray(1), 1, tmb%lzd%llr(ilr_tmp)%wfd%keygloc(1,1), 1)
          call vcopy(2*ncount, workarray(2*ncount+1), 1, tmb%lzd%llr(ilr_tmp)%wfd%keyglob(1,1), 1)
          call vcopy(ncount, workarray(4*ncount+1), 1, tmb%lzd%llr(ilr_tmp)%wfd%keyvloc(1), 1)
          call vcopy(ncount, workarray(5*ncount+1), 1, tmb%lzd%llr(ilr_tmp)%wfd%keyvglob(1), 1)
      end if
      call f_free(workarray)
  end if

  ! find biggest instead
  !do ilr=1,tmb%lzr%nlr
  !  if (tmb%lzd%llr(ilr)%wfd%nvctr_c
  !end do

  ! Determine size of phi_old and phi
  ndim_tmp=0
  ndim=0
  ndim_tmp1=tmb%lzd%llr(ilr_tmp)%wfd%nvctr_c+7*tmb%lzd%llr(ilr_tmp)%wfd%nvctr_f
  do iorb=1,tmb%orbs%norbp
      iiorb=tmb%orbs%isorb+iorb
      ilr=tmb%orbs%inwhichlocreg(iiorb)
      ndim=ndim+tmb%lzd%llr(ilr)%wfd%nvctr_c+7*tmb%lzd%llr(ilr)%wfd%nvctr_f
      ndim_tmp=ndim_tmp+ndim_tmp1
  end do

  ! should integrate bettwer with existing reformat routines, but restart needs tidying anyway
  psi_tmp = f_malloc_ptr(ndim_tmp,id='psi_tmp')

  jstart=1
  jstart_tmp=1

  n_tmp(1)=tmb%lzd%Llr(ilr_tmp)%d%n1
  n_tmp(2)=tmb%lzd%Llr(ilr_tmp)%d%n2
  n_tmp(3)=tmb%lzd%Llr(ilr_tmp)%d%n3

  ns_tmp(1)=tmb%lzd%Llr(ilr_tmp)%ns1
  ns_tmp(2)=tmb%lzd%Llr(ilr_tmp)%ns2
  ns_tmp(3)=tmb%lzd%Llr(ilr_tmp)%ns3

  nglr(1)=tmb%lzd%glr%d%n1
  nglr(2)=tmb%lzd%glr%d%n2
  nglr(3)=tmb%lzd%glr%d%n3

  do iorb=1,tmb%orbs%norbp
      iiorb=tmb%orbs%isorb+iorb
      ilr=tmb%orbs%inwhichlocreg(iiorb)
      iiat=tmb%orbs%onwhichatom(iiorb)

      n(1)=tmb%lzd%Llr(ilr)%d%n1
      n(2)=tmb%lzd%Llr(ilr)%d%n2
      n(3)=tmb%lzd%Llr(ilr)%d%n3

      ns(1)=tmb%lzd%Llr(ilr)%ns1
      ns(2)=tmb%lzd%Llr(ilr)%ns2
      ns(3)=tmb%lzd%Llr(ilr)%ns3

      !theta=0.d0*(4.0_gp*atan(1.d0)/180.0_gp)
      !newz=(/1.0_gp,0.0_gp,0.0_gp/)
      !centre_old(:)=rxyz(:,iiat)
      !centre_new(:)=rxyz(:,iiat_tmp)
      !shift(:)=centre_new(:)-centre_old(:)

      frag_trans=fragment_transformation_identity()
!!$      frag_trans%theta=0.0d0*(4.0_gp*atan(1.d0)/180.0_gp)
!!$      frag_trans%rot_axis=(/1.0_gp,0.0_gp,0.0_gp/)
      frag_trans%rot_center(:)=rxyz(:,iiat)
      frag_trans%rot_center_new(:)=rxyz(:,iiat_tmp)

      call reformat_check(reformat,reformat_reason,tol,at,tmb%lzd%hgrids,tmb%lzd%hgrids,&
           tmb%lzd%llr(ilr)%wfd%nvctr_c,tmb%lzd%llr(ilr)%wfd%nvctr_f,&
           tmb%lzd%llr(ilr_tmp)%wfd%nvctr_c,tmb%lzd%llr(ilr_tmp)%wfd%nvctr_f,&
           n,n_tmp,ns,ns_tmp,nglr,nglr,at%astruct%geocode,& !tmb%lzd%llr(ilr_tmp)%geocode,&
           frag_trans,centre_old_box,centre_new_box,da,wrap_around)

      if ((.not. reformat) .and. (.not. wrap_around)) then ! copy psi into psi_tmp

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
          phigold = f_malloc((/ 0.to.n(1), 1.to.2, 0.to.n(2), 1.to.2, 0.to.n(3), 1.to.2 /),id='phigold')

          !!!debug
          !!open(3000+iiorb)
          !!do i=1,tmb%lzd%llr(ilr)%wfd%nvctr_c+7*tmb%lzd%llr(ilr)%wfd%nvctr_f
          !!write(3000+iiorb,*) i,tmb%psi(jstart+i)
          !!end do
          !!close(3000+iiorb)

          call psi_to_psig(n,tmb%lzd%llr(ilr)%wfd%nseg_c,tmb%lzd%llr(ilr)%wfd%nvctr_c,&
               tmb%lzd%llr(ilr)%wfd%keygloc,tmb%lzd%llr(ilr)%wfd%keyvloc,&
               tmb%lzd%llr(ilr)%wfd%nseg_f,tmb%lzd%llr(ilr)%wfd%nvctr_f,&
               tmb%lzd%llr(ilr)%wfd%keygloc(1:,tmb%lzd%Llr(ilr)%wfd%nseg_c+1:), &
               tmb%lzd%llr(ilr)%wfd%keyvloc(tmb%lzd%Llr(ilr)%wfd%nseg_c+1:), &
               phigold,tmb%psi(jstart),tmb%psi(jstart+tmb%lzd%llr(ilr)%wfd%nvctr_c))

          if (reformat) then

             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
             !!! reformatting using psir - can only be done if keys are communicated for ilr_tmp
             !!workarraytmp=f_malloc((2*n+31),id='workarraytmp')
             !!psirold=f_malloc0((2*n+31),id='psirold')

             !!!call f_zero((2*n_old(1)+31)*(2*n_old(2)+31)*(2*n_old(3)+31),psirold(1,1,1))
             !!call vcopy((2*n(1)+2)*(2*n(2)+2)*(2*n(3)+2),phigold(0,1,0,1,0,1),1,psirold(1,1,1),1)
             !!call psig_to_psir_free(n(1),n(2),n(3),workarraytmp,psirold)
             !!call f_free(workarraytmp)

             !!per(1)=(at%astruct%geocode /= 'F')
             !!per(2)=(at%astruct%geocode == 'P')
             !!per(3)=(at%astruct%geocode /= 'F')

             !!!buffers related to periodicity
             !!!WARNING: the boundary conditions are not assumed to change between new and old
             !!call ext_buffers(per(1),nl(1),nr(1))
             !!call ext_buffers(per(2),nl(2),nr(2))
             !!call ext_buffers(per(3),nl(3),nr(3))

             !!! centre of rotation with respect to start of box
             !!centre_old_box(1)=frag_trans%rot_center(1)-0.5d0*tmb%lzd%hgrids(1)*(tmb%lzd%llr(ilr)%nsi1-nl(1))
             !!centre_old_box(2)=frag_trans%rot_center(2)-0.5d0*tmb%lzd%hgrids(2)*(tmb%lzd%llr(ilr)%nsi2-nl(2))
             !!centre_old_box(3)=frag_trans%rot_center(3)-0.5d0*tmb%lzd%hgrids(3)*(tmb%lzd%llr(ilr)%nsi3-nl(3))

             !!centre_new_box(1)=frag_trans%rot_center_new(1)-0.5d0*tmb%lzd%hgrids(1)*(tmb%lzd%llr(ilr_tmp)%nsi1-nl(1))
             !!centre_new_box(2)=frag_trans%rot_center_new(2)-0.5d0*tmb%lzd%hgrids(2)*(tmb%lzd%llr(ilr_tmp)%nsi2-nl(2))
             !!centre_new_box(3)=frag_trans%rot_center_new(3)-0.5d0*tmb%lzd%hgrids(3)*(tmb%lzd%llr(ilr_tmp)%nsi3-nl(3))


             !!da=centre_new_box-centre_old_box-(tmb%lzd%hgrids-tmb%lzd%hgrids)*0.5d0

             !!call reformat_one_supportfunction(tmb%lzd%llr(ilr_tmp),tmb%lzd%llr(ilr),tmb%lzd%llr(ilr_tmp)%geocode,&
             !!     & tmb%lzd%hgrids,n,phigold,tmb%lzd%hgrids,n_tmp,centre_old_box,centre_new_box,da,&
             !!     & frag_trans,psi_tmp(jstart_tmp:),psirold)

             !!call f_free(psirold)
             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

             call reformat_one_supportfunction(tmb%lzd%llr(ilr_tmp),tmb%lzd%llr(ilr),at%astruct%geocode,& !,tmb%lzd%llr(ilr_tmp)%geocode,&
                  & tmb%lzd%hgrids,n,phigold,tmb%lzd%hgrids,n_tmp,centre_old_box,centre_new_box,da,&
                  & frag_trans,psi_tmp(jstart_tmp:))

          else

             ! in this case we don't need to reformat, just re-wrap the tmb, so ilr and ilr_tmp should contain same info
             call compress_plain(n_tmp(1),n_tmp(2),0,n_tmp(1),0,n_tmp(2),0,n_tmp(3),  &
                  tmb%lzd%llr(ilr_tmp)%wfd%nseg_c,tmb%lzd%llr(ilr_tmp)%wfd%nvctr_c,&
                  tmb%lzd%llr(ilr_tmp)%wfd%keygloc(1,1),tmb%lzd%llr(ilr_tmp)%wfd%keyvloc(1),   &
                  tmb%lzd%llr(ilr_tmp)%wfd%nseg_f,tmb%lzd%llr(ilr_tmp)%wfd%nvctr_f,&
                  tmb%lzd%llr(ilr_tmp)%wfd%keygloc(1,tmb%lzd%llr(ilr_tmp)%wfd%nseg_c+min(1,tmb%lzd%llr(ilr_tmp)%wfd%nseg_f)),&
                  tmb%lzd%llr(ilr_tmp)%wfd%keyvloc(tmb%lzd%llr(ilr_tmp)%wfd%nseg_c+min(1,tmb%lzd%llr(ilr_tmp)%wfd%nseg_f)), &
                  phigold,psi_tmp(jstart_tmp),psi_tmp(jstart_tmp+tmb%lzd%llr(ilr_tmp)%wfd%nvctr_c))

          end if

          jstart=jstart+tmb%lzd%llr(ilr)%wfd%nvctr_c+7*tmb%lzd%llr(ilr)%wfd%nvctr_f
          jstart_tmp=jstart_tmp+ndim_tmp1

          call f_free(phigold)

      end if

      !!!debug
      !!open(1000+iiorb)
      !!do i=1,ndim_tmp1
      !!   write(1000+iiorb,*) i,psi_tmp(jstart_tmp-ndim_tmp1+i-1),&
      !!        ddot(ndim_tmp1, psi_tmp(jstart_tmp-ndim_tmp1), 1, psi_tmp(jstart_tmp-ndim_tmp1), 1)
      !!end do
      !!close(1000+iiorb)

      !!!debug
      !!write(*,'(a,5(1x,I3),(1x,F8.4),5(2x,3(1x,F6.2)),5(2x,3(1x,I4)),2(1x,L2),4(1x,I5))')'DEBUGnr:',&
      !!     iproc,iiat,iiorb,ilr,ilr_tmp,ddot(ndim_tmp1, psi_tmp(jstart_tmp-ndim_tmp1), 1, psi_tmp(jstart_tmp-ndim_tmp1), 1), &
      !!     da, frag_trans%rot_center, frag_trans%rot_center_new, centre_old_box, centre_new_box, ns, ns_tmp, n, n_tmp, nglr, &
      !!     reformat, wrap_around, tmb%lzd%llr(ilr)%wfd%nvctr_c, tmb%lzd%llr(ilr)%wfd%nvctr_f, &
      !!     tmb%lzd%llr(ilr_tmp)%wfd%nvctr_c, tmb%lzd%llr(ilr_tmp)%wfd%nvctr_f

  end do

  call print_reformat_summary(iproc,nproc,reformat_reason)

  ! now that they are all in one lr, need to calculate overlap matrix
  ! make lzd_tmp contain all identical lrs
  lzd_tmp = local_zone_descriptors_null()
  lzd_tmp%linear=tmb%lzd%linear
  lzd_tmp%nlr=tmb%lzd%nlr
  lzd_tmp%lintyp=tmb%lzd%lintyp
  lzd_tmp%ndimpotisf=tmb%lzd%ndimpotisf
  lzd_tmp%hgrids(:)=tmb%lzd%hgrids(:)

  call nullify_locreg_descriptors(lzd_tmp%glr)
  call copy_locreg_descriptors(tmb%lzd%glr, lzd_tmp%glr)

  iis1=lbound(tmb%lzd%llr,1)
  iie1=ubound(tmb%lzd%llr,1)
  allocate(lzd_tmp%llr(iis1:iie1))

  do i1=iis1,iie1
     call nullify_locreg_descriptors(lzd_tmp%llr(i1))
     call copy_locreg_descriptors(tmb%lzd%llr(ilr_tmp), lzd_tmp%llr(i1))
  end do


  !call nullify_comms_linear(collcom_tmp)
  collcom_tmp=comms_linear_null()
  call init_comms_linear(iproc, nproc, imethod_overlap, ndim_tmp, tmb%orbs, lzd_tmp, &
       tmb%linmat%m%nspin, collcom_tmp)

  smat_tmp = sparse_matrix_null()
  ! Do not initialize the matrix multiplication to save memory. 
  call init_sparse_matrix_wrapper(iproc, nproc, tmb%linmat%s%nspin, tmb%orbs, &
       lzd_tmp, at%astruct, .false., init_matmul=.false., imode=2, smat=smat_tmp)
  call init_matrixindex_in_compressed_fortransposed(iproc, nproc, &
       collcom_tmp, collcom_tmp, collcom_tmp, smat_tmp)
  iirow(1) = smat_tmp%nfvctr
  iirow(2) = 1
  iicol(1) = smat_tmp%nfvctr
  iicol(2) = 1
  call check_local_matrix_extents(iproc, nproc, &
       collcom_tmp, collcom_tmp, tmb%linmat%smmd, smat_tmp, &
       ind_min, ind_mas, ind_trans_min, ind_trans_max, irow, icol)
  iirow(1) = min(irow(1),iirow(1))
  iirow(2) = max(irow(2),iirow(2))
  iicol(1) = min(icol(1),iicol(1))
  iicol(2) = max(icol(2),iicol(2))

  call init_matrix_taskgroups(iproc, nproc, bigdft_mpi%mpi_comm, .false., smat_tmp)

  mat_tmp = matrices_null()
  mat_tmp%matrix_compr = sparsematrix_malloc_ptr(smat_tmp, iaction=SPARSE_TASKGROUP,id='mat_tmp%matrix_compr')

  psit_c_tmp = f_malloc_ptr(sum(collcom_tmp%nrecvcounts_c),id='psit_c_tmp')
  psit_f_tmp = f_malloc_ptr(7*sum(collcom_tmp%nrecvcounts_f),id='psit_f_tmp')

  call transpose_localized(iproc, nproc, ndim_tmp, tmb%orbs, collcom_tmp, &
       TRANSPOSE_FULL, psi_tmp, psit_c_tmp, psit_f_tmp, lzd_tmp)

  ! normalize psi
  !skip the normalize psi step
  !norm = f_malloc_ptr(tmb%orbs%norb,id='norm')
  !call normalize_transposed(iproc, nproc, tmb%orbs, tmb%linmat%s%nspin, collcom_tmp, psit_c_tmp, psit_f_tmp, norm)
  !call f_free_ptr(norm)

  !!call calculate_pulay_overlap(iproc, nproc, tmb%orbs, tmb%orbs, collcom_tmp, collcom_tmp, &
  !!     psit_c_tmp, psit_c_tmp, psit_f_tmp, psit_f_tmp, tmb%linmat%ovrlp_%matrix)
  call calculate_overlap_transposed(iproc, nproc, tmb%orbs, collcom_tmp, &
                 psit_c_tmp, psit_c_tmp, psit_f_tmp, psit_f_tmp, smat_tmp, mat_tmp)
  !call uncompress_matrix(iproc, tmb%linmat%s, mat_tmp%matrix_compr, tmb%linmat%ovrlp_%matrix)
  call uncompress_matrix(iproc, nproc, smat_tmp, mat_tmp%matrix_compr, tmb%linmat%ovrlp_%matrix)

  call deallocate_matrices(mat_tmp)
  call deallocate_sparse_matrix(smat_tmp)

!!!# DEBUG #######
!!call deallocate_local_zone_descriptors(lzd_tmp)
!!call mpi_finalize(i1)
!!stop
!!!# END DEBUG ###

  call deallocate_comms_linear(collcom_tmp)
  call deallocate_local_zone_descriptors(lzd_tmp)

  call f_free_ptr(psit_c_tmp)
  call f_free_ptr(psit_f_tmp)

  call f_free_ptr(psi_tmp)

END SUBROUTINE tmb_overlap_onsite


! we want to compare every tmb to every other tmb, calculating rotations using environment info
subroutine tmb_overlap_onsite_rotate(iproc, nproc, input, at, tmb, rxyz, ref_frags)
  use module_base
  use module_types
  use module_fragments
  use io, only: find_neighbours
  use locregs, only: allocate_wfd, deallocate_wfd
  use module_interfaces, only: reformat_one_supportfunction

  implicit none

  ! Calling arguments
  integer,intent(in) :: iproc, nproc
  type(input_variables),intent(in) :: input
  type(atoms_data), intent(in) :: at
  type(DFT_wavefunction),intent(inout):: tmb
  real(kind=gp),dimension(3,at%astruct%nat),intent(in) :: rxyz
  type(system_fragment), dimension(input%frag%nfrag_ref), intent(in) :: ref_frags

  !local arguments
  integer :: num_neighbours, iat, jat, ntmb_frag_and_env, iiorb, iforb, iiat, iatf, ityp
  integer :: iiat_tmp, ilr_tmp, norb_tmp, jproc, iroot, ncount, j, ndim, ndim_tmp, iorb, ilr
  integer :: jstart, jstart_tmp, ndim_tmp1, istart_tmp, jorb, jlr, jjorb, jjat, iorba, jorba
  integer :: ifr, ifr_ref, ntmb_frag_and_env_dfrag, num_env, i, k
  integer, dimension(3) :: ns, ns_tmp, n, n_tmp, nglr
  integer, dimension(0:7) :: reformat_reason
  integer, dimension(:), allocatable :: workarray, ifrag_ref
  integer, dimension(:,:), allocatable :: map_frag_and_env, frag_map

  real(kind=gp) :: tol, ddot
  real(kind=gp), dimension(3) :: centre_old_box, centre_new_box, da
  real(kind=wp), dimension(:), pointer :: psi_tmp
  real(kind=wp), dimension(:), allocatable :: psi_tmp_i, psi_tmp_j
  real(kind=gp), dimension(:,:), allocatable :: overlap
  real(kind=wp), dimension(:,:,:,:,:,:), allocatable :: phigold

  type(system_fragment), dimension(:), allocatable :: ref_frags_atomic, ref_frags_atomic_dfrag
  type(fragment_transformation), dimension(:,:), pointer :: frag_trans

  logical :: reformat, wrap_around
  !!real(wp) :: dnrm2

  ! check which fragment we're on - if this isn't a fragment calculation then all atoms are in the same fragment, which could be very expensive!

  ifrag_ref=f_malloc(at%astruct%nat,id='ifrag_ref')
  if (input%frag%nfrag>1) then
     iat=0
     do ifr=1,input%frag%nfrag
        ifr_ref=input%frag%frag_index(ifr)
        do iatf=1,ref_frags(ifr_ref)%astruct_frg%nat
           iat=iat+1
           ifrag_ref(iat)=ifr_ref
           !print*,'frag',iat,ifr_ref,input%frag%nfrag
        end do
     end do
  else
     ifrag_ref=1
  end if

  ! allocate and set up ref_frags_atomic, treating each atom as if it is a fragment
  allocate(ref_frags_atomic(at%astruct%nat))
  allocate(ref_frags_atomic_dfrag(at%astruct%nat))

  do iat=1,at%astruct%nat
     call setup_frags_from_astruct(ref_frags_atomic(iat))
     call setup_frags_from_astruct(ref_frags_atomic_dfrag(iat))
  end do


  ! take a default value of 4 if this hasn't been specified in the input
  if (input%lin%frag_num_neighbours==0) then
     num_neighbours=4
  else
     num_neighbours=input%lin%frag_num_neighbours
  end if


  allocate(frag_trans(at%astruct%nat,at%astruct%nat))

  ! pre-fill with identity transformation
  do iat=1,at%astruct%nat
     do jat=1,at%astruct%nat
        frag_trans(iat,jat)=fragment_transformation_identity()
        frag_trans(iat,jat)%rot_center_new=frag_center(1,rxyz(:,jat))
        frag_trans(iat,jat)%Werror=-1.0d0
     end do
     ! diagonal terms have zero error
     frag_trans(iat,jat)%Werror=0.0d0
  end do

  ! get all fragment transformations first, then reformat
  do iat=1,at%astruct%nat
     frag_map=f_malloc0((/ref_frags_atomic(iat)%fbasis%forbs%norb,3/),id='frag_map')

     ! tmbs on the same atom should be consecutive, so we just need to find the first tmb for this atom
     do iiorb=1,tmb%orbs%norb
        iiat=tmb%orbs%onwhichatom(iiorb)
        if (iiat==iat) exit
     end do

     do iforb=1,ref_frags_atomic(iat)%fbasis%forbs%norb
        !tmb frag -> tmb full
        frag_map(iforb,1)=iiorb+iforb-1
        !tmb frag -> atom frag
        frag_map(iforb,2)=1
     end do
     !atom frag -> atom full
     frag_map(1,3)=iat

     ! find environment atoms
     ! don't think we care about keeping track of which atoms they were so we can immediately destroy the mapping array
     map_frag_and_env = f_malloc((/tmb%orbs%norb,3/),id='map_frag_and_env')
     ! version with most extensive matching
     call find_neighbours(num_neighbours,at,rxyz,tmb%orbs,ref_frags_atomic(iat),frag_map,&
          ntmb_frag_and_env,map_frag_and_env,.false.,input%lin%frag_neighbour_cutoff)
     ! with only closest shell
     call find_neighbours(num_neighbours,at,rxyz,tmb%orbs,ref_frags_atomic_dfrag(iat),frag_map,&
          ntmb_frag_and_env_dfrag,map_frag_and_env,.true.,input%lin%frag_neighbour_cutoff)
     call f_free(map_frag_and_env)

     ! we also need nbasis_env
     ref_frags_atomic(iat)%nbasis_env=0
     do iatf=1,ref_frags_atomic(iat)%astruct_env%nat
        ityp=ref_frags_atomic(iat)%astruct_env%iatype(iatf)
        ref_frags_atomic(iat)%nbasis_env=ref_frags_atomic(iat)%nbasis_env+input%lin%norbsPerType(ityp)
     end do

     ! we also need nbasis_env
     ref_frags_atomic_dfrag(iat)%nbasis_env=0
     do iatf=1,ref_frags_atomic_dfrag(iat)%astruct_env%nat
        ityp=ref_frags_atomic_dfrag(iat)%astruct_env%iatype(iatf)
        ref_frags_atomic_dfrag(iat)%nbasis_env=ref_frags_atomic_dfrag(iat)%nbasis_env+input%lin%norbsPerType(ityp)
     end do

     do jat=1,at%astruct%nat
        ! now we treat ref_frags_atomic and environment data therein as if it came from a file
        ! and find transformation between that and current atom
        ! in case where Wahba gives a huge error just do no transformation?****************** (actually maybe want to ignore Si atom i.e. those further away?)

        ! for different atom types don't bother checking for transformation (avoids problems if ntmb/atom doesn't match)
        if (at%astruct%iatype(iat)==at%astruct%iatype(jat)) then
           ! as above, we don't need to keep track of mapping, we just want the transformation
           map_frag_and_env = f_malloc((/ref_frags_atomic(iat)%nbasis_env,3/),id='map_frag_and_env')

           ! if we are an identical fragment then do full matching
           if (ifrag_ref(iat)==ifrag_ref(jat)) then
              call match_environment_atoms(jat-1,at,rxyz,tmb%orbs,ref_frags_atomic(iat),&
                   ref_frags_atomic(iat)%nbasis_env,map_frag_and_env,frag_trans(jat,iat),.false.)
           ! otherwise just look at nearest neighbours, i.e. n closest, not n closest of each type
           else
              call match_environment_atoms(jat-1,at,rxyz,tmb%orbs,ref_frags_atomic_dfrag(iat),&
                   ref_frags_atomic_dfrag(iat)%nbasis_env,map_frag_and_env,frag_trans(jat,iat),.true.)
           end if

           call f_free(map_frag_and_env)

           !if (frag_trans(jat,iat)%Werror > W_tol) call f_increment(itoo_big)
        end if

     end do

     call f_free(frag_map)
  end do

  !debug - check calculated transformations
  if (iproc==0) then
     open(99)
     do iat=1,at%astruct%nat
        do jat=1,at%astruct%nat
           if (at%astruct%iatype(iat)/=at%astruct%iatype(jat)) then
              num_env=-1
           else if (ifrag_ref(iat)==ifrag_ref(jat)) then
              num_env=ref_frags_atomic(iat)%astruct_env%nat-ref_frags_atomic(iat)%astruct_frg%nat
           else
              num_env=ref_frags_atomic_dfrag(iat)%astruct_env%nat-ref_frags_atomic_dfrag(iat)%astruct_frg%nat
           end if
           write(99,'(2(a,1x,I5,1x),F12.6,2x,3(F12.6,1x),6(1x,F18.6),2x,F6.2,3(2x,I5))') &
                trim(at%astruct%atomnames(at%astruct%iatype(iat))),iat,&
                trim(at%astruct%atomnames(at%astruct%iatype(jat))),jat,&
                frag_trans(iat,jat)%theta/(4.0_gp*atan(1.d0)/180.0_gp),frag_trans(iat,jat)%rot_axis,&
                frag_trans(iat,jat)%rot_center,frag_trans(iat,jat)%rot_center_new,&
                frag_trans(iat,jat)%Werror,num_env,ifrag_ref(iat),ifrag_ref(jat)
        end do
     end do
     close(99)
  end if


  ! NOW WE ACTUALLY NEED TO REFORMAT

  ! move all psi into psi_tmp all centred in the same place and calculate overlap matrix
  tol=1.d-3
  reformat_reason=0

  !arbitrarily pick the middle one as assuming it'll be near the centre of structure
  !and therefore have large fine grid
  !might be more efficient to reformat each jorb to iorb's lzd, but that involves communicating lots of keys...
  norb_tmp=tmb%orbs%norb/2
  ilr_tmp=tmb%orbs%inwhichlocreg(norb_tmp)
  iiat_tmp=tmb%orbs%onwhichatom(norb_tmp)

  ! Find out which process handles TMB norb_tmp and get the keys from that process
  do jproc=0,nproc-1
      if (tmb%orbs%isorb_par(jproc)<norb_tmp .and. norb_tmp<=tmb%orbs%isorb_par(jproc)+tmb%orbs%norb_par(jproc,0)) then
          iroot=jproc
          exit
      end if
  end do
  if (iproc/=iroot) then
      ! some processes might already have it allocated
      call deallocate_wfd(tmb%lzd%llr(ilr_tmp)%wfd)
      call allocate_wfd(tmb%lzd%llr(ilr_tmp)%wfd)
  end if
  if (nproc>1) then
      ncount = tmb%lzd%llr(ilr_tmp)%wfd%nseg_c + tmb%lzd%llr(ilr_tmp)%wfd%nseg_f
      workarray = f_malloc(6*ncount,id='workarray')
      if (iproc==iroot) then
          call vcopy(2*ncount, tmb%lzd%llr(ilr_tmp)%wfd%keygloc(1,1), 1, workarray(1), 1)
          call vcopy(2*ncount, tmb%lzd%llr(ilr_tmp)%wfd%keyglob(1,1), 1, workarray(2*ncount+1), 1)
          call vcopy(ncount, tmb%lzd%llr(ilr_tmp)%wfd%keyvloc(1), 1, workarray(4*ncount+1), 1)
          call vcopy(ncount, tmb%lzd%llr(ilr_tmp)%wfd%keyvglob(1), 1, workarray(5*ncount+1), 1)
      end if
      call mpibcast(workarray, root=iroot, comm=bigdft_mpi%mpi_comm)
      if (iproc/=iroot) then
          call vcopy(2*ncount, workarray(1), 1, tmb%lzd%llr(ilr_tmp)%wfd%keygloc(1,1), 1)
          call vcopy(2*ncount, workarray(2*ncount+1), 1, tmb%lzd%llr(ilr_tmp)%wfd%keyglob(1,1), 1)
          call vcopy(ncount, workarray(4*ncount+1), 1, tmb%lzd%llr(ilr_tmp)%wfd%keyvloc(1), 1)
          call vcopy(ncount, workarray(5*ncount+1), 1, tmb%lzd%llr(ilr_tmp)%wfd%keyvglob(1), 1)
      end if
      call f_free(workarray)
  end if


  ndim_tmp1=tmb%lzd%llr(ilr_tmp)%wfd%nvctr_c+7*tmb%lzd%llr(ilr_tmp)%wfd%nvctr_f

  ! Determine size of phi_old and phi
  ndim_tmp=0
  ndim=0
  do iorb=1,tmb%orbs%norbp
      iiorb=tmb%orbs%isorb+iorb
      ilr=tmb%orbs%inwhichlocreg(iiorb)
      ndim=ndim+tmb%lzd%llr(ilr)%wfd%nvctr_c+7*tmb%lzd%llr(ilr)%wfd%nvctr_f
      ndim_tmp=ndim_tmp+ndim_tmp1
  end do

  ! should integrate bettwer with existing reformat routines, but restart needs tidying anyway
  psi_tmp = f_malloc_ptr(ndim_tmp,id='psi_tmp')

  ! since we have a lot of different rotations, store all (reformatted) psi_i but only 1 psi_j at a time
  ! bypass calculate_overlap_transposed and do directly in dense format
  ! think about how to do this better...

  ! get all psi_i - i.e. identity transformation, but reformatting to a single lr
  n_tmp(1)=tmb%lzd%Llr(ilr_tmp)%d%n1
  n_tmp(2)=tmb%lzd%Llr(ilr_tmp)%d%n2
  n_tmp(3)=tmb%lzd%Llr(ilr_tmp)%d%n3
  ns_tmp(1)=tmb%lzd%Llr(ilr_tmp)%ns1
  ns_tmp(2)=tmb%lzd%Llr(ilr_tmp)%ns2
  ns_tmp(3)=tmb%lzd%Llr(ilr_tmp)%ns3
  nglr(1)=tmb%lzd%glr%d%n1
  nglr(2)=tmb%lzd%glr%d%n2
  nglr(3)=tmb%lzd%glr%d%n3

  jstart=1
  jstart_tmp=1
  do iorb=1,tmb%orbs%norbp
      iiorb=tmb%orbs%isorb+iorb
      ilr=tmb%orbs%inwhichlocreg(iiorb)
      iiat=tmb%orbs%onwhichatom(iiorb)

      n(1)=tmb%lzd%Llr(ilr)%d%n1
      n(2)=tmb%lzd%Llr(ilr)%d%n2
      n(3)=tmb%lzd%Llr(ilr)%d%n3
      ns(1)=tmb%lzd%Llr(ilr)%ns1
      ns(2)=tmb%lzd%Llr(ilr)%ns2
      ns(3)=tmb%lzd%Llr(ilr)%ns3

      ! override centres
      frag_trans(iiat,iiat)%rot_center(:)=rxyz(:,iiat)
      frag_trans(iiat,iiat)%rot_center_new(:)=rxyz(:,iiat_tmp)

      call reformat_check(reformat,reformat_reason,tol,at,tmb%lzd%hgrids,tmb%lzd%hgrids,&
           tmb%lzd%llr(ilr)%wfd%nvctr_c,tmb%lzd%llr(ilr)%wfd%nvctr_f,&
           tmb%lzd%llr(ilr_tmp)%wfd%nvctr_c,tmb%lzd%llr(ilr_tmp)%wfd%nvctr_f,&
           n,n_tmp,ns,ns_tmp,nglr,nglr,at%astruct%geocode,& !,tmb%lzd%llr(ilr_tmp)%geocode,&
           frag_trans(iiat,iiat),centre_old_box,centre_new_box,da,wrap_around)  

      if ((.not. reformat) .and. (.not. wrap_around)) then ! copy psi into psi_tmp
          do j=1,tmb%lzd%llr(ilr_tmp)%wfd%nvctr_c
              psi_tmp(jstart_tmp)=tmb%psi(jstart)
              jstart=jstart+1
              jstart_tmp=jstart_tmp+1
          end do
          do j=1,7*tmb%lzd%llr(ilr_tmp)%wfd%nvctr_f-6,7
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
          phigold = f_malloc((/ 0.to.n(1), 1.to.2, 0.to.n(2), 1.to.2, 0.to.n(3), 1.to.2 /),id='phigold')

          call psi_to_psig(n,tmb%lzd%llr(ilr)%wfd%nseg_c,tmb%lzd%llr(ilr)%wfd%nvctr_c,&
               tmb%lzd%llr(ilr)%wfd%keygloc,tmb%lzd%llr(ilr)%wfd%keyvloc,&
               tmb%lzd%llr(ilr)%wfd%nseg_f,tmb%lzd%llr(ilr)%wfd%nvctr_f,&
               tmb%lzd%llr(ilr)%wfd%keygloc(1:,tmb%lzd%Llr(ilr)%wfd%nseg_c+1:), &
               tmb%lzd%llr(ilr)%wfd%keyvloc(tmb%lzd%Llr(ilr)%wfd%nseg_c+1:), &
               phigold,tmb%psi(jstart),tmb%psi(jstart+tmb%lzd%llr(ilr)%wfd%nvctr_c))

          if (reformat) then

             call reformat_one_supportfunction(tmb%lzd%llr(ilr_tmp),tmb%lzd%llr(ilr),at%astruct%geocode,& !,tmb%lzd%llr(ilr_tmp)%geocode,&
                  & tmb%lzd%hgrids,n,phigold,tmb%lzd%hgrids,n_tmp,centre_old_box,centre_new_box,da,&
                  & frag_trans(iiat,iiat),psi_tmp(jstart_tmp:))

          else

             ! in this case we don't need to reformat, just re-wrap the tmb, so ilr and ilr_tmp should contain same info
             call compress_plain(n_tmp(1),n_tmp(2),0,n_tmp(1),0,n_tmp(2),0,n_tmp(3),  &
                  tmb%lzd%llr(ilr_tmp)%wfd%nseg_c,tmb%lzd%llr(ilr_tmp)%wfd%nvctr_c,&
                  tmb%lzd%llr(ilr_tmp)%wfd%keygloc(1,1),tmb%lzd%llr(ilr_tmp)%wfd%keyvloc(1),   &
                  tmb%lzd%llr(ilr_tmp)%wfd%nseg_f,tmb%lzd%llr(ilr_tmp)%wfd%nvctr_f,&
                  tmb%lzd%llr(ilr_tmp)%wfd%keygloc(1,tmb%lzd%llr(ilr_tmp)%wfd%nseg_c+min(1,tmb%lzd%llr(ilr_tmp)%wfd%nseg_f)),&
                  tmb%lzd%llr(ilr_tmp)%wfd%keyvloc(tmb%lzd%llr(ilr_tmp)%wfd%nseg_c+min(1,tmb%lzd%llr(ilr_tmp)%wfd%nseg_f)), &
                  phigold,psi_tmp(jstart_tmp),psi_tmp(jstart_tmp+tmb%lzd%llr(ilr_tmp)%wfd%nvctr_c))
 
          end if

          jstart=jstart+tmb%lzd%llr(ilr)%wfd%nvctr_c+7*tmb%lzd%llr(ilr)%wfd%nvctr_f
          jstart_tmp=jstart_tmp+ndim_tmp1

          call f_free(phigold)

      end if

      !!!debug
      !!write(*,'(a,5(1x,I3),(1x,F8.4),5(2x,3(1x,F6.2)),5(2x,3(1x,I4)),2(1x,L2),4(1x,I5),2(1x,a))')'DEBUGr:',&
      !!     iproc,iiat,iiorb,ilr,ilr_tmp,ddot(ndim_tmp1, psi_tmp(jstart_tmp-ndim_tmp1), 1, psi_tmp(jstart_tmp-ndim_tmp1), 1), &
      !!     da, frag_trans(iiat,iiat)%rot_center, frag_trans(iiat,iiat)%rot_center_new, &
      !!     centre_old_box, centre_new_box, ns, ns_tmp, n, n_tmp, nglr, &
      !!     reformat, wrap_around, tmb%lzd%llr(ilr)%wfd%nvctr_c, tmb%lzd%llr(ilr)%wfd%nvctr_f, &
      !!     tmb%lzd%llr(ilr_tmp)%wfd%nvctr_c, tmb%lzd%llr(ilr_tmp)%wfd%nvctr_f, tmb%lzd%llr(ilr)%geocode, tmb%lzd%glr%geocode

  end do

  !not sure whether to print this - maybe actual reformatting only?
  !call print_reformat_summary(iproc,nproc,reformat_reason)

  psi_tmp_i=f_malloc(ndim_tmp1,id='psi_tmp_i')
  psi_tmp_j=f_malloc(ndim_tmp1,id='psi_tmp_j')

  !maybe make this an argument?
  overlap=f_malloc0((/tmb%orbs%norb,tmb%orbs%norb/),id='overlap')

  istart_tmp=1
  do iiorb=1,tmb%orbs%norb
      !iiorb=tmb%orbs%isorb+iorb
      ilr=tmb%orbs%inwhichlocreg(iiorb)
      iiat=tmb%orbs%onwhichatom(iiorb)

      !these have already been reformatted, so need to communicate each reformatted psi_tmp
      !first check which mpi has this tmb
      do jproc=0,nproc-1
          if (tmb%orbs%isorb_par(jproc)<iiorb .and. iiorb<=tmb%orbs%isorb_par(jproc)+tmb%orbs%norb_par(jproc,0)) then
              iroot=jproc
              exit
          end if
      end do

      if (nproc>1) then
          !workarray = f_malloc(ndim_tmp1,id='workarray')
          if (iproc==iroot) then
              !call vcopy(2*ncount, tmb%lzd%llr(ilr_tmp)%wfd%keygloc(1,1), 1, workarray(1), 1)
              call vcopy(ndim_tmp1, psi_tmp(istart_tmp), 1, psi_tmp_i(1), 1)
          end if
          call mpibcast(psi_tmp_i, root=iroot, comm=bigdft_mpi%mpi_comm)
          if (iproc==iroot) then
              istart_tmp = istart_tmp + ndim_tmp1
          !else
          !    call vcopy(2*ncount, workarray(1), 1, tmb%lzd%llr(ilr_tmp)%wfd%keygloc(1,1), 1)
          !    call vcopy(2*ncount, workarray(2*ncount+1), 1, tmb%lzd%llr(ilr_tmp)%wfd%keyglob(1,1), 1)
          !    call vcopy(ncount, workarray(4*ncount+1), 1, tmb%lzd%llr(ilr_tmp)%wfd%keyvloc(1), 1)
          !    call vcopy(ncount, workarray(5*ncount+1), 1, tmb%lzd%llr(ilr_tmp)%wfd%keyvglob(1), 1)
          end if
          !call f_free(workarray)
      end if

      ! reformat jorb for iorb then calculate overlap
      jstart=1
      do jorb=1,tmb%orbs%norbp
         jjorb=tmb%orbs%isorb+jorb
         jlr=tmb%orbs%inwhichlocreg(jjorb)
         jjat=tmb%orbs%onwhichatom(jjorb)

         ! not sure if this really saves us much, or just makes for poor load balancing
         if (jjat>iiat) then !  .or. (jjat/=7 .and. jjat/=8) .or. (iiat/=7 .and. iiat/=8)   ) then
            jstart=jstart+tmb%lzd%llr(jlr)%wfd%nvctr_c+7*tmb%lzd%llr(jlr)%wfd%nvctr_f
            !overlap(jjorb,iiorb) = 0.0
            cycle
         end if

         n(1)=tmb%lzd%Llr(jlr)%d%n1
         n(2)=tmb%lzd%Llr(jlr)%d%n2
         n(3)=tmb%lzd%Llr(jlr)%d%n3
         ns(1)=tmb%lzd%Llr(jlr)%ns1
         ns(2)=tmb%lzd%Llr(jlr)%ns2
         ns(3)=tmb%lzd%Llr(jlr)%ns3

         ! override centres
         frag_trans(iiat,jjat)%rot_center(:)=rxyz(:,jjat)
         frag_trans(iiat,jjat)%rot_center_new(:)=rxyz(:,iiat_tmp)

         call reformat_check(reformat,reformat_reason,tol,at,tmb%lzd%hgrids,tmb%lzd%hgrids,&
              tmb%lzd%llr(jlr)%wfd%nvctr_c,tmb%lzd%llr(jlr)%wfd%nvctr_f,&
              tmb%lzd%llr(ilr_tmp)%wfd%nvctr_c,tmb%lzd%llr(ilr_tmp)%wfd%nvctr_f,&
              n,n_tmp,ns,ns_tmp,nglr,nglr,at%astruct%geocode,& !,tmb%lzd%llr(ilr_tmp)%geocode,&
              frag_trans(iiat,jjat),centre_old_box,centre_new_box,da,wrap_around)

         if ((.not. reformat) .and. (.not. wrap_around)) then ! copy psi into psi_tmp
             jstart_tmp=1
             do j=1,tmb%lzd%llr(ilr_tmp)%wfd%nvctr_c
                 psi_tmp_j(jstart_tmp)=tmb%psi(jstart)
                 jstart=jstart+1
                 jstart_tmp=jstart_tmp+1
             end do
             do j=1,7*tmb%lzd%llr(ilr_tmp)%wfd%nvctr_f-6,7
                 psi_tmp_j(jstart_tmp+0)=tmb%psi(jstart+0)
                 psi_tmp_j(jstart_tmp+1)=tmb%psi(jstart+1)
                 psi_tmp_j(jstart_tmp+2)=tmb%psi(jstart+2)
                 psi_tmp_j(jstart_tmp+3)=tmb%psi(jstart+3)
                 psi_tmp_j(jstart_tmp+4)=tmb%psi(jstart+4)
                 psi_tmp_j(jstart_tmp+5)=tmb%psi(jstart+5)
                 psi_tmp_j(jstart_tmp+6)=tmb%psi(jstart+6)
                 jstart=jstart+7
                 jstart_tmp=jstart_tmp+7
             end do

         else
             phigold = f_malloc((/ 0.to.n(1), 1.to.2, 0.to.n(2), 1.to.2, 0.to.n(3), 1.to.2 /),id='phigold')

             call psi_to_psig(n,tmb%lzd%llr(jlr)%wfd%nseg_c,tmb%lzd%llr(jlr)%wfd%nvctr_c,&
                  tmb%lzd%llr(jlr)%wfd%keygloc,tmb%lzd%llr(jlr)%wfd%keyvloc,&
                  tmb%lzd%llr(jlr)%wfd%nseg_f,tmb%lzd%llr(jlr)%wfd%nvctr_f,&
                  tmb%lzd%llr(jlr)%wfd%keygloc(1:,tmb%lzd%Llr(jlr)%wfd%nseg_c+1:), &
                  tmb%lzd%llr(jlr)%wfd%keyvloc(tmb%lzd%Llr(jlr)%wfd%nseg_c+1:), &
                  phigold,tmb%psi(jstart),tmb%psi(jstart+tmb%lzd%llr(jlr)%wfd%nvctr_c))

         !!!debug
         !!!if ((iiat==7 .or. iiat==8) .and. (jjat==7 .or. jjat==8)) then
         !!   open(20000+iiorb*100+jjorb)
         !!   do i=1,n(1)
         !!   do j=1,n(2)
         !!   do k=1,n(3)
         !!      write(20000+iiorb*100+jjorb,'(3(I6,1x),1x,2(F12.6,1x))') i,j,k,phigold(i,1,j,1,k,1),&
         !!           dnrm2(8*(n(1)+1)*(n(2)+1)*(n(3)+1),phigold,1)
         !!   end do
         !!   end do
         !!   end do
         !!   close(20000+iiorb*100+jjorb)
         !!!end if

             if (reformat) then

                call reformat_one_supportfunction(tmb%lzd%llr(ilr_tmp),tmb%lzd%llr(jlr),at%astruct%geocode,&  !tmb%lzd%llr(ilr_tmp)%geocode,&
                     & tmb%lzd%hgrids,n,phigold,tmb%lzd%hgrids,n_tmp,centre_old_box,centre_new_box,da,&
                     & frag_trans(iiat,jjat),psi_tmp_j) !,tag=30000+iiorb*100+jjorb)

             else
                ! in this case we don't need to reformat, just re-wrap the tmb, so ilr and ilr_tmp should contain same info
                call compress_plain(n_tmp(1),n_tmp(2),0,n_tmp(1),0,n_tmp(2),0,n_tmp(3),  &
                     tmb%lzd%llr(ilr_tmp)%wfd%nseg_c,tmb%lzd%llr(ilr_tmp)%wfd%nvctr_c,&
                     tmb%lzd%llr(ilr_tmp)%wfd%keygloc(1,1),tmb%lzd%llr(ilr_tmp)%wfd%keyvloc(1),   &
                     tmb%lzd%llr(ilr_tmp)%wfd%nseg_f,tmb%lzd%llr(ilr_tmp)%wfd%nvctr_f,&
                     tmb%lzd%llr(ilr_tmp)%wfd%keygloc(1,tmb%lzd%llr(ilr_tmp)%wfd%nseg_c+min(1,tmb%lzd%llr(ilr_tmp)%wfd%nseg_f)),&
                     tmb%lzd%llr(ilr_tmp)%wfd%keyvloc(tmb%lzd%llr(ilr_tmp)%wfd%nseg_c+min(1,tmb%lzd%llr(ilr_tmp)%wfd%nseg_f)), &
                     phigold,psi_tmp_j(1),psi_tmp_j(tmb%lzd%llr(ilr_tmp)%wfd%nvctr_c+1))
             end if

             jstart=jstart+tmb%lzd%llr(jlr)%wfd%nvctr_c+7*tmb%lzd%llr(jlr)%wfd%nvctr_f   

             call f_free(phigold)
         end if

         overlap(iiorb,jjorb) = ddot(ndim_tmp1, psi_tmp_i(1), 1, psi_tmp_j(1), 1)
         overlap(jjorb,iiorb) = overlap(iiorb,jjorb)

         !!!debug
         !!write(*,'(a,5(1x,I3),1x,2(1x,F8.4),1x,2(1x,L2),1x,2(1x,F7.2))')'DEBUGr2:',&
         !!     iproc,iiat,jjat,iiorb,jjorb,ddot(ndim_tmp1, psi_tmp_j(1), 1, psi_tmp_j(1), 1), &
         !!     ddot(ndim_tmp1, psi_tmp_i(1), 1, psi_tmp_j(1), 1), &
         !!     reformat, wrap_around, frag_trans(iiat,jjat)%Werror, frag_trans(iiat,jjat)%theta/(4.0_gp*atan(1.d0)/180.0_gp)

         !!!debug
         !!write(*,'(a,5(1x,I4),2(2x,F12.6))')'iproc,iat,jat,iorb,jorb,ovrlp',iproc,iiat,jjat,iiorb,jjorb,overlap(jjorb,iiorb),&
         !!      tmb%linmat%ovrlp_%matrix(iiorb,jjorb,1)

         !!!debug
         !!!if ((iiat==7 .or. iiat==8) .and. (jjat==7 .or. jjat==8)) then
         !!   open(10000+iiorb*100+jjorb)
         !!   do i=1,ndim_tmp1
         !!      write(10000+iiorb*100+jjorb,'(I6,1x,2(F12.6,1x),1x,3(F12.6,1x))') i,psi_tmp_i(i),psi_tmp_j(i),&
         !!           ddot(ndim_tmp1, psi_tmp_i(1), 1, psi_tmp_i(1), 1),&
         !!           ddot(ndim_tmp1, psi_tmp_j(1), 1, psi_tmp_j(1), 1),&
         !!           ddot(ndim_tmp1, psi_tmp_i(1), 1, psi_tmp_j(1), 1)
         !!   end do
         !!   close(10000+iiorb*100+jjorb)
         !!!end if

      end do

      if (iproc==0) write(*,'(F6.2,a)') 100.0d0*real(iiorb,dp)/real(tmb%orbs%norb,dp),'%'
  end do

  call mpiallred(overlap, mpi_sum, comm=bigdft_mpi%mpi_comm)

  ! print also various files useful for more direct analysis - eventually tidy this into better formatted outputs
  if (iproc==0) then
     open(98)
     open(97)
     open(96)
     open(95)
     do iat=1,at%astruct%nat
        iorba=0
        do iorb=1,tmb%orbs%norb
           iiat=tmb%orbs%onwhichatom(iorb)
           if (iat/=iiat) cycle
           iorba=iorba+1

           do jat=1,at%astruct%nat
              jorba=0
              do jorb=1,tmb%orbs%norb
                 jjat=tmb%orbs%onwhichatom(jorb)
                 if (jat/=jjat) cycle
                 jorba=jorba+1

                 ! full matrix
                 write(98,'(2(a,1x,I5,1x),4(I5,1x),2x,5(F12.6,1x))') &
                      trim(at%astruct%atomnames(at%astruct%iatype(iat))),iat,&
                      trim(at%astruct%atomnames(at%astruct%iatype(jat))),jat,&
                      iorb,jorb,iorba,jorba,overlap(iorb,jorb),&
                      tmb%linmat%ovrlp_%matrix(iorb,jorb,1),&
                      overlap(iorb,jorb)-tmb%linmat%ovrlp_%matrix(iorb,jorb,1),&
                      frag_trans(iat,jat)%Werror,frag_trans(iat,jat)%theta/(4.0_gp*atan(1.d0)/180.0_gp)

                 ! same 'number' tmb, no rotation
                 if ((iat==jat .or. frag_trans(iat,jat)%theta==0.0) .and. iorba==jorba) then
                    write(97,'(2(a,1x,I5,1x),4(I5,1x),2x,5(F12.6,1x))') &
                         trim(at%astruct%atomnames(at%astruct%iatype(iat))),iat,&
                         trim(at%astruct%atomnames(at%astruct%iatype(jat))),jat,&
                         iorb,jorb,iorba,jorba,overlap(iorb,jorb),&
                         tmb%linmat%ovrlp_%matrix(iorb,jorb,1),&
                         overlap(iorb,jorb)-tmb%linmat%ovrlp_%matrix(iorb,jorb,1),&
                         frag_trans(iat,jat)%Werror,frag_trans(iat,jat)%theta/(4.0_gp*atan(1.d0)/180.0_gp)
                 end if

                 ! same 'number' tmb
                 if (iorba==jorba) then
                    write(96,'(2(a,1x,I5,1x),4(I5,1x),2x,5(F12.6,1x))') &
                         trim(at%astruct%atomnames(at%astruct%iatype(iat))),iat,&
                         trim(at%astruct%atomnames(at%astruct%iatype(jat))),jat,&
                         iorb,jorb,iorba,jorba,overlap(iorb,jorb),&
                         tmb%linmat%ovrlp_%matrix(iorb,jorb,1),&
                         overlap(iorb,jorb)-tmb%linmat%ovrlp_%matrix(iorb,jorb,1),&
                         frag_trans(iat,jat)%Werror,frag_trans(iat,jat)%theta/(4.0_gp*atan(1.d0)/180.0_gp)
                 end if

                 ! first tmb of each atom
                 if (iorba==jorba .and. iorba==1) then
                    write(95,'(2(a,1x,I5,1x),4(I5,1x),2x,5(F12.6,1x))') &
                         trim(at%astruct%atomnames(at%astruct%iatype(iat))),iat,&
                         trim(at%astruct%atomnames(at%astruct%iatype(jat))),jat,&
                         iorb,jorb,iorba,jorba,overlap(iorb,jorb),&
                         tmb%linmat%ovrlp_%matrix(iorb,jorb,1),&
                         overlap(iorb,jorb)-tmb%linmat%ovrlp_%matrix(iorb,jorb,1),&
                         frag_trans(iat,jat)%Werror,frag_trans(iat,jat)%theta/(4.0_gp*atan(1.d0)/180.0_gp)
                 end if

              end do
           end do
        end do
     end do
     close(98)
     close(97)
     close(96)
     close(95)
  end if


  call f_free(overlap)
  call f_free(psi_tmp_i)
  call f_free(psi_tmp_j)
  call f_free_ptr(psi_tmp)

  do iat=1,at%astruct%nat
     ! nullify these as we are pointing to at%astruct versions
     nullify(ref_frags_atomic(iat)%astruct_frg%atomnames)
     nullify(ref_frags_atomic(iat)%astruct_env%atomnames)
     call fragment_free(ref_frags_atomic(iat))

     nullify(ref_frags_atomic_dfrag(iat)%astruct_frg%atomnames)
     nullify(ref_frags_atomic_dfrag(iat)%astruct_env%atomnames)
     call fragment_free(ref_frags_atomic_dfrag(iat))
  end do
  deallocate(ref_frags_atomic_dfrag)
  deallocate(frag_trans)
  call f_free(ifrag_ref)


contains


  subroutine setup_frags_from_astruct(ref_frags)
     implicit none
     type(system_fragment), intent(out) :: ref_frags

     ref_frags=fragment_null()

     ! fill in relevant ref_frags info
     ref_frags%astruct_frg%nat=1
     ref_frags%astruct_env%nat=0

     ! copy some stuff from astruct
     ref_frags%astruct_frg%inputfile_format = at%astruct%inputfile_format
     ref_frags%astruct_frg%units = at%astruct%units
     ref_frags%astruct_frg%geocode = at%astruct%geocode
     ! coordinates can just be zero as we only have 1 atom
     ref_frags%astruct_frg%rxyz = f_malloc0_ptr((/3,1/),id=' ref_frags%astruct_frg%rxyz')

     ! now deal with atom types - easier to keep this coherent with at%astruct rather than having one type
     ref_frags%astruct_frg%ntypes = at%astruct%ntypes
     ref_frags%astruct_frg%atomnames => at%astruct%atomnames
     ref_frags%astruct_frg%iatype = f_malloc_ptr(1,id='ref_frags%astruct_frg%iatype')

     ! polarization etc is irrelevant
     ref_frags%astruct_frg%input_polarization = f_malloc0_ptr(1,&
          id='ref_frags%astruct_frg%input_polarization')
     ref_frags%astruct_frg%ifrztyp = f_malloc0_ptr(1,id='ref_frags%astruct_frg%ifrztyp')
     ref_frags%astruct_frg%iatype(1) = at%astruct%iatype(iat)
     ref_frags%astruct_frg%input_polarization(1) = 100

     call init_minimal_orbitals_data(iproc, nproc, 1, input, ref_frags%astruct_frg, &
          ref_frags%fbasis%forbs, at%astruct)

  end subroutine setup_frags_from_astruct



END SUBROUTINE tmb_overlap_onsite_rotate





!!subroutine tmb_overlap_onsite_rotate(iproc, nproc, at, tmb, rxyz)
!!
!!  use module_base
!!  use module_types
!!  use module_interfaces
!!  use module_fragments
!!  use communications_base, only: comms_linear_null, deallocate_comms_linear
!!  use communications_init, only: init_comms_linear
!!  use communications, only: transpose_localized
!!  implicit none
!!
!!  ! Calling arguments
!!  integer,intent(in) :: iproc, nproc
!!  type(atoms_data), intent(inout) :: at
!!  type(DFT_wavefunction),intent(in):: tmb
!!  real(gp),dimension(3,at%astruct%nat),intent(in) :: rxyz
!!
!!  ! Local variables
!!  logical :: reformat
!!  integer :: iorb,i_stat,i_all,istart,jstart
!!  integer :: iiorb,ilr,iiat,j,iis1,iie1,i1
!!  integer :: jlr,iiat_tmp,ndim_tmp,ndim,norb_tmp,iat,jjat,jjorb
!!  integer, dimension(3) :: ns,nsj,n,nj
!!  real(gp), dimension(3) :: centre_old_box, centre_new_box, da
!!  real(wp), dimension(:,:,:,:,:,:), allocatable :: phigold
!!  real(wp), dimension(:), pointer :: psi_tmp, psit_c_tmp, psit_f_tmp, norm
!!  integer, dimension(0:6) :: reformat_reason
!!  type(comms_linear) :: collcom_tmp
!!  type(local_zone_descriptors) :: lzd_tmp
!!  real(gp) :: tol
!!  character(len=*),parameter:: subname='tmb_overlap_onsite'
!!  type(fragment_transformation) :: frag_trans
!!  real(gp), dimension(:,:), allocatable :: rxyz4_new, rxyz4_ref, rxyz_new, rxyz_ref
!!  real(gp), dimension(:), allocatable :: dist
!!  integer, dimension(:), allocatable :: ipiv
!!
!!  ! move all psi into psi_tmp all centred in the same place and calculate overlap matrix
!!  tol=1.d-3
!!  reformat_reason=0
!!
!!  !arbitrarily pick the middle one as assuming it'll be near the centre of structure
!!  !and therefore have large fine grid
!!  norb_tmp=tmb%orbs%norb/2
!!  jlr=tmb%orbs%inwhichlocreg(norb_tmp)
!!  iiat_tmp=tmb%orbs%onwhichatom(norb_tmp)
!!  jjorb=norb_tmp
!!  jjat=iiat_tmp
!!
!!  ! find biggest instead
!!  !do ilr=1,tmb%lzr%nlr
!!  !  if (tmb%lzd%llr(ilr)%wfd%nvctr_c
!!  !end do
!!
!!  ! Determine size of phi_old and phi
!!  ndim_tmp=0
!!  ndim=0
!!  do iorb=1,tmb%orbs%norbp
!!      iiorb=tmb%orbs%isorb+iorb
!!      ilr=tmb%orbs%inwhichlocreg(iiorb)
!!      ndim=ndim+tmb%lzd%llr(ilr)%wfd%nvctr_c+7*tmb%lzd%llr(ilr)%wfd%nvctr_f
!!      ndim_tmp=ndim_tmp+tmb%lzd%llr(jlr)%wfd%nvctr_c+7*tmb%lzd%llr(jlr)%wfd%nvctr_f
!!  end do
!!
!!  ! should integrate bettwer with existing reformat routines, but restart needs tidying anyway
!!  allocate(psi_tmp(ndim_tmp),stat=i_stat)
!!  call memocc(i_stat,psi_tmp,'psi_tmp',subname)
!!
!!  allocate(rxyz4_ref(3,min(4,at%astruct%nat)), stat=i_stat)
!!  call memocc(i_stat, rxyz4_ref, 'rxyz4_ref', subname)
!!  allocate(rxyz4_new(3,min(4,at%astruct%nat)), stat=i_stat)
!!  call memocc(i_stat, rxyz4_new, 'rxyz4_ref', subname)
!!
!!  allocate(rxyz_ref(3,at%astruct%nat), stat=i_stat)
!!  call memocc(i_stat, rxyz_ref, 'rxyz_ref', subname)
!!  allocate(rxyz_new(3,at%astruct%nat), stat=i_stat)
!!  call memocc(i_stat, rxyz_new, 'rxyz_ref', subname)
!!  allocate(dist(at%astruct%nat), stat=i_stat)
!!  call memocc(i_stat, dist, 'dist', subname)
!!  allocate(ipiv(at%astruct%nat), stat=i_stat)
!!  call memocc(i_stat, ipiv, 'ipiv', subname)
!!
!!  istart=1
!!  jstart=1
!!  do iorb=1,tmb%orbs%norbp
!!      iiorb=tmb%orbs%isorb+iorb
!!      ilr=tmb%orbs%inwhichlocreg(iiorb)
!!      iiat=tmb%orbs%onwhichatom(iiorb)
!!
!!      !do jjorb=1,tmb%orbs%norb
!!      !   !jlr=tmb%orbs%inwhichlocreg(jjorb)
!!      !   jjat=tmb%orbs%onwhichatom(jjorb)
!!
!!         n(1)=tmb%lzd%Llr(ilr)%d%n1
!!         n(2)=tmb%lzd%Llr(ilr)%d%n2
!!         n(3)=tmb%lzd%Llr(ilr)%d%n3
!!         nj(1)=tmb%lzd%Llr(jlr)%d%n1
!!         nj(2)=tmb%lzd%Llr(jlr)%d%n2
!!         nj(3)=tmb%lzd%Llr(jlr)%d%n3
!!         ns(1)=tmb%lzd%Llr(ilr)%ns1
!!         ns(2)=tmb%lzd%Llr(ilr)%ns2
!!         ns(3)=tmb%lzd%Llr(ilr)%ns3
!!         nsj(1)=tmb%lzd%Llr(jlr)%ns1
!!         nsj(2)=tmb%lzd%Llr(jlr)%ns2
!!         nsj(3)=tmb%lzd%Llr(jlr)%ns3
!!
!!         ! find fragment transformation using 3 nearest neighbours
!!         do iat=1,at%astruct%nat
!!            rxyz_new(:,iat)=rxyz(:,iat)
!!            rxyz_ref(:,iat)=rxyz(:,iat)
!!         end do
!!
!!         ! use atom position
!!         frag_trans%rot_center=rxyz(:,iiat)
!!         frag_trans%rot_center_new=rxyz(:,jjat)
!!
!!        ! shift rxyz wrt center of rotation
!!         do iat=1,at%astruct%nat
!!            rxyz_ref(:,iat)=rxyz_ref(:,iat)-frag_trans%rot_center
!!            rxyz_new(:,iat)=rxyz_new(:,iat)-frag_trans%rot_center_new
!!         end do
!!
!!         ! find distances from this atom, sort atoms into neighbour order and take atom and 3 nearest neighbours
!!         do iat=1,at%astruct%nat
!!            dist(iat)=-dsqrt(rxyz_ref(1,iat)**2+rxyz_ref(2,iat)**2+rxyz_ref(3,iat)**2)
!!         end do
!!         call sort_positions(at%astruct%nat,dist,ipiv)
!!         do iat=1,min(4,at%astruct%nat)
!!            rxyz4_ref(:,iat)=rxyz_ref(:,ipiv(iat))
!!         end do
!!
!!         ! find distances from this atom, sort atoms into neighbour order and take atom and 3 nearest neighbours
!!         do iat=1,at%astruct%nat
!!            dist(iat)=-dsqrt(rxyz_new(1,iat)**2+rxyz_new(2,iat)**2+rxyz_new(3,iat)**2)
!!         end do
!!         call sort_positions(at%astruct%nat,dist,ipiv)
!!         do iat=1,min(4,at%astruct%nat)
!!            rxyz4_new(:,iat)=rxyz_new(:,ipiv(iat))
!!         end do
!!
!!         call find_frag_trans(min(4,at%astruct%nat),rxyz4_ref,rxyz4_new,frag_trans)
!!
!!         write(*,'(A,4(I3,1x),3(F12.6,1x),F12.6)') 'iorb,jorb,iat,jat,rot_axis,theta',&
!!              iiorb,jjorb,iiat,jjat,frag_trans%rot_axis,frag_trans%theta/(4.0_gp*atan(1.d0)/180.0_gp)
!!
!!         !frag_trans%theta=0.0d0*(4.0_gp*atan(1.d0)/180.0_gp)
!!         !frag_trans%rot_axis=(/1.0_gp,0.0_gp,0.0_gp/)
!!         !frag_trans%rot_center(:)=rxyz(:,iiat)
!!         !overwrite rot_center_new to account for llr_tmp being in different location
!!         frag_trans%rot_center_new(:)=rxyz(:,iiat_tmp)
!!
!!         call reformat_check(reformat,reformat_reason,tol,at,tmb%lzd%hgrids,tmb%lzd%hgrids,&
!!              tmb%lzd%llr(ilr)%wfd%nvctr_c,tmb%lzd%llr(ilr)%wfd%nvctr_f,&
!!              tmb%lzd%llr(jlr)%wfd%nvctr_c,tmb%lzd%llr(jlr)%wfd%nvctr_f,&
!!              n,nj,ns,nsj,frag_trans,centre_old_box,centre_new_box,da)
!!
!!         if (.not. reformat) then ! copy psi into psi_tmp
!!            do j=1,tmb%lzd%llr(jlr)%wfd%nvctr_c
!!               psi_tmp(jstart)=tmb%psi(istart)
!!               istart=istart+1
!!               jstart=jstart+1
!!            end do
!!            do j=1,7*tmb%lzd%llr(ilr)%wfd%nvctr_f-6,7
!!               psi_tmp(jstart+0)=tmb%psi(istart+0)
!!               psi_tmp(jstart+1)=tmb%psi(istart+1)
!!               psi_tmp(jstart+2)=tmb%psi(istart+2)
!!               psi_tmp(jstart+3)=tmb%psi(istart+3)
!!               psi_tmp(jstart+4)=tmb%psi(istart+4)
!!               psi_tmp(jstart+5)=tmb%psi(istart+5)
!!               psi_tmp(jstart+6)=tmb%psi(istart+6)
!!               istart=istart+7
!!               jstart=jstart+7
!!            end do
!!         else
!!            allocate(phigold(0:n(1),2,0:n(2),2,0:n(3),2+ndebug),stat=i_stat)
!!            call memocc(i_stat,phigold,'phigold',subname)
!!
!!            call psi_to_psig(n,tmb%lzd%llr(ilr)%wfd%nvctr_c,tmb%lzd%llr(ilr)%wfd%nvctr_f,&
!!                 tmb%lzd%llr(ilr)%wfd%nseg_c,tmb%lzd%llr(ilr)%wfd%nseg_f,&
!!                 tmb%lzd%llr(ilr)%wfd%keyvloc,tmb%lzd%llr(ilr)%wfd%keygloc,istart,tmb%psi(jstart),phigold)
!!
!!            call reformat_one_supportfunction(tmb%lzd%llr(jlr),tmb%lzd%llr(ilr),tmb%lzd%llr(jlr)%geocode,&
!!                 tmb%lzd%hgrids,n,phigold,tmb%lzd%hgrids,nj,centre_old_box,centre_new_box,da,&
!!                 frag_trans,psi_tmp(jstart:))
!!
!!            jstart=jstart+tmb%lzd%llr(jlr)%wfd%nvctr_c+7*tmb%lzd%llr(jlr)%wfd%nvctr_f
!!
!!            i_all=-product(shape(phigold))*kind(phigold)
!!            deallocate(phigold,stat=i_stat)
!!            call memocc(i_stat,i_all,'phigold',subname)
!!        end if
!!
!!     !end do
!!  end do
!!
!!
!!  i_all = -product(shape(ipiv))*kind(ipiv)
!!  deallocate(ipiv,stat=i_stat)
!!  call memocc(i_stat,i_all,'ipiv',subname)
!!  i_all = -product(shape(dist))*kind(dist)
!!  deallocate(dist,stat=i_stat)
!!  call memocc(i_stat,i_all,'dist',subname)
!!  i_all = -product(shape(rxyz_ref))*kind(rxyz_ref)
!!  deallocate(rxyz_ref,stat=i_stat)
!!  call memocc(i_stat,i_all,'rxyz_ref',subname)
!!  i_all = -product(shape(rxyz_new))*kind(rxyz_new)
!!  deallocate(rxyz_new,stat=i_stat)
!!  call memocc(i_stat,i_all,'rxyz_new',subname)
!!  i_all = -product(shape(rxyz4_ref))*kind(rxyz4_ref)
!!  deallocate(rxyz4_ref,stat=i_stat)
!!  call memocc(i_stat,i_all,'rxyz4_ref',subname)
!!  i_all = -product(shape(rxyz4_new))*kind(rxyz4_new)
!!  deallocate(rxyz4_new,stat=i_stat)
!!  call memocc(i_stat,i_all,'rxyz4_new',subname)
!!
!!  call print_reformat_summary(iproc,reformat_reason)
!!
!!  ! now that they are all in one lr, need to calculate overlap matrix
!!  ! make lzd_tmp contain all identical lrs
!!  lzd_tmp%linear=tmb%lzd%linear
!!  lzd_tmp%nlr=tmb%lzd%nlr
!!  lzd_tmp%lintyp=tmb%lzd%lintyp
!!  lzd_tmp%ndimpotisf=tmb%lzd%ndimpotisf
!!  lzd_tmp%hgrids(:)=tmb%lzd%hgrids(:)
!!
!!  call nullify_locreg_descriptors(lzd_tmp%glr)
!!  call copy_locreg_descriptors(tmb%lzd%glr, lzd_tmp%glr)
!!
!!  iis1=lbound(tmb%lzd%llr,1)
!!  iie1=ubound(tmb%lzd%llr,1)
!!  allocate(lzd_tmp%llr(iis1:iie1), stat=i_stat)
!!  do i1=iis1,iie1
!!     call nullify_locreg_descriptors(lzd_tmp%llr(i1))
!!     call copy_locreg_descriptors(tmb%lzd%llr(jlr), lzd_tmp%llr(i1))
!!  end do
!!
!!  !call nullify_comms_linear(collcom_tmp)
!!  collcom_tmp=comms_linear_null()
!!  call init_comms_linear(iproc, nproc, ndim_tmp, tmb%orbs, lzd_tmp, collcom_tmp)
!!
!!  allocate(psit_c_tmp(sum(collcom_tmp%nrecvcounts_c)), stat=i_stat)
!!  call memocc(i_stat, psit_c_tmp, 'psit_c_tmp', subname)
!!
!!  allocate(psit_f_tmp(7*sum(collcom_tmp%nrecvcounts_f)), stat=i_stat)
!!  call memocc(i_stat, psit_f_tmp, 'psit_f_tmp', subname)
!!
!!  call transpose_localized(iproc, nproc, ndim_tmp, tmb%orbs, collcom_tmp, &
!!       psi_tmp, psit_c_tmp, psit_f_tmp, lzd_tmp)
!!
!!  ! normalize psi
!!  allocate(norm(tmb%orbs%norb), stat=i_stat)
!!  call memocc(i_stat, norm, 'norm', subname)
!!  call normalize_transposed(iproc, nproc, tmb%orbs, collcom_tmp, psit_c_tmp, psit_f_tmp, norm)
!!  i_all = -product(shape(norm))*kind(norm)
!!  deallocate(norm,stat=i_stat)
!!  call memocc(i_stat,i_all,'norm',subname)
!!
!!  call calculate_pulay_overlap(iproc, nproc, tmb%orbs, tmb%orbs, collcom_tmp, collcom_tmp, &
!!       psit_c_tmp, psit_c_tmp, psit_f_tmp, psit_f_tmp, tmb%linmat%ovrlp%matrix)
!!
!!  call deallocate_comms_linear(collcom_tmp)
!!  call deallocate_local_zone_descriptors(lzd_tmp)
!!
!!  i_all = -product(shape(psit_c_tmp))*kind(psit_c_tmp)
!!  deallocate(psit_c_tmp,stat=i_stat)
!!  call memocc(i_stat,i_all,'psit_c_tmp',subname)
!!
!!  i_all = -product(shape(psit_f_tmp))*kind(psit_f_tmp)
!!  deallocate(psit_f_tmp,stat=i_stat)
!!  call memocc(i_stat,i_all,'psit_f_tmp',subname)
!!
!!  i_all = -product(shape(psi_tmp))*kind(psi_tmp)
!!  deallocate(psi_tmp,stat=i_stat)
!!  call memocc(i_stat,i_all,'psi_tmp',subname)
!!
!!END SUBROUTINE tmb_overlap_onsite_rotate


!!!not correct - not sure if it's worth fixing as only memguess uses it currently
!!subroutine readonewave_linear(unitwf,useFormattedInput,iorb,iproc,n,ns,&
!!     & hgrids,at,llr,rxyz_old,rxyz,locrad,locregCenter,confPotOrder,&
!!     & confPotprefac,psi,eval,onwhichatom,lr,glr,reformat_reason)
!!  use module_base
!!  use module_types
!!  !use internal_io
!!  use module_interfaces, only: reformat_one_supportfunction
!!  use yaml_output
!!  use module_fragments
!!  use io, only: io_read_descr_linear, io_error, read_psi_compress
!!  implicit none
!!  logical, intent(in) :: useFormattedInput
!!  integer, intent(in) :: unitwf,iorb,iproc
!!  integer, dimension(3), intent(in) :: n,ns
!!  !type(wavefunctions_descriptors), intent(in) :: wfd
!!  type(locreg_descriptors), intent(in) :: llr
!!  type(atoms_data), intent(in) :: at
!!  real(gp), dimension(3), intent(in) :: hgrids
!!  real(gp), dimension(3,at%astruct%nat), intent(in) :: rxyz
!!  integer, intent(out) :: confPotOrder
!!  real(gp), intent(out) :: locrad, confPotprefac
!!  real(wp), intent(out) :: eval
!!  real(gp), dimension(3), intent(out) :: locregCenter
!!  real(gp), dimension(3,at%astruct%nat), intent(out) :: rxyz_old
!!  real(wp), dimension(:), pointer :: psi
!!  integer, dimension(*), intent(in) :: onwhichatom
!!  type(locreg_descriptors), intent(in) :: lr, glr
!!  integer, dimension(0:6), intent(out) :: reformat_reason
!!
!!  !local variables
!!  character(len=*), parameter :: subname='readonewave_linear'
!!  character(len = 256) :: error
!!  logical :: lstat,reformat
!!  integer :: iorb_old,nvctr_c_old,nvctr_f_old,iiat
!!  integer :: i1,i2,i3,iel,onwhichatom_tmp
!!  integer, dimension(3) :: ns_old,n_old
!!  real(gp) :: tol
!!  real(gp), dimension(3) :: hgrids_old, centre_old_box, centre_new_box, da
!!  real(gp) :: tt,t1,t2,t3,t4,t5,t6,t7
!!  real(wp), dimension(:,:,:,:,:,:), allocatable :: psigold
!!  type(fragment_transformation) :: frag_trans
!!  ! DEBUG
!!  ! character(len=12) :: orbname
!!  ! real(wp), dimension(:), allocatable :: gpsi
!!
!!
!!  call io_read_descr_linear(unitwf, useFormattedInput, iorb_old, eval, n_old(1), n_old(2), n_old(3), &
!!       ns_old(1), ns_old(2), ns_old(3), hgrids_old, lstat, error, onwhichatom_tmp, locrad, &
!!       locregCenter, confPotOrder, confPotprefac, nvctr_c_old, nvctr_f_old, at%astruct%nat, rxyz_old)
!!
!!  if (.not. lstat) call io_error(trim(error))
!!  if (iorb_old /= iorb) stop 'readonewave_linear'
!!
!!  iiat=onwhichatom(iorb)
!!  tol=1.d-3
!!
!!  !why these hard-coded values?
!!  frag_trans%theta=20.0d0*(4.0_gp*atan(1.d0)/180.0_gp)
!!  frag_trans%rot_axis=(/1.0_gp,0.0_gp,0.0_gp/)
!!  frag_trans%rot_center(:)=(/7.8d0,11.8d0,11.6d0/)
!!  frag_trans%rot_center_new(:)=(/7.8d0,11.2d0,11.8d0/)
!!
!!  call reformat_check(reformat,reformat_reason,tol,at,hgrids,hgrids_old,&
!!       nvctr_c_old,nvctr_f_old,llr%wfd%nvctr_c,llr%wfd%nvctr_f,&
!!       n_old,n,ns_old,ns,frag_trans,centre_old_box,centre_new_box,da)
!!
!!  if (.not. reformat) then
!!     call read_psi_compress(unitwf, useFormattedInput, nvctr_c_old, nvctr_f_old, psi, lstat, error)
!!     if (.not. lstat) call io_error(trim(error))
!!  else
!!     ! add derivative functions at a later date? (needs orbs and lzd)
!!     psigold = f_malloc0((/ 0.to.n_old(1), 1.to.2, 0.to.n_old(2), 1.to.2, 0.to.n_old(3), 1.to.2 /),id='psigold')
!!
!!     !call f_zero(8*(n_old(1)+1)*(n_old(2)+1)*(n_old(3)+1),psigold)
!!     do iel=1,nvctr_c_old
!!        if (useFormattedInput) then
!!           read(unitwf,*) i1,i2,i3,tt
!!        else
!!           read(unitwf) i1,i2,i3,tt
!!        end if
!!        psigold(i1,1,i2,1,i3,1)=tt
!!     enddo
!!     do iel=1,nvctr_f_old
!!        if (useFormattedInput) then
!!           read(unitwf,*) i1,i2,i3,t1,t2,t3,t4,t5,t6,t7
!!        else
!!           read(unitwf) i1,i2,i3,t1,t2,t3,t4,t5,t6,t7
!!        end if
!!        psigold(i1,2,i2,1,i3,1)=t1
!!        psigold(i1,1,i2,2,i3,1)=t2
!!        psigold(i1,2,i2,2,i3,1)=t3
!!        psigold(i1,1,i2,1,i3,2)=t4
!!        psigold(i1,2,i2,1,i3,2)=t5
!!        psigold(i1,1,i2,2,i3,2)=t6
!!        psigold(i1,2,i2,2,i3,2)=t7
!!     enddo
!!
!!     ! NB assuming here geocode is the same in glr and llr
!!     call reformat_one_supportfunction(llr,llr,at%astruct%geocode,hgrids_old,n_old,psigold,hgrids,n, &
!!         centre_old_box,centre_new_box,da,frag_trans,psi)
!!
!!     call f_free(psigold)
!!
!!  endif
!!
!!  !! DEBUG - plot in global box - CHECK WITH REFORMAT ETC IN LRs
!!  !allocate (gpsi(glr%wfd%nvctr_c+7*glr%wfd%nvctr_f),stat=i_stat)
!!  !call memocc(i_stat,gpsi,'gpsi',subname)
!!  !
!!  !call f_zero(glr%wfd%nvctr_c+7*glr%wfd%nvctr_f,gpsi)
!!  !call Lpsi_to_global2(iproc, llr%%wfd%nvctr_c+7*lr%wfd%nvctr_f, glr%wfd%nvctr_c+7*glr%wfd%nvctr_f, &
!!  !     1, 1, 1, glr, lr, psi, gpsi)
!!  !
!!  !write(orbname,*) iorb
!!  !call plot_wf(trim(adjustl(orbname)),1,at,1.0_dp,glr,hgrids(1),hgrids(2),hgrids(3),rxyz,gpsi)
!!  !!call plot_wf(trim(adjustl(orbname)),1,at,1.0_dp,lr,hx,hy,hz,rxyz,psi)
!!  !
!!  !i_all=-product(shape(gpsi))*kind(gpsi)
!!  !deallocate(gpsi,stat=i_stat)
!!  !call memocc(i_stat,i_all,'gpsi',subname)
!!  !! END DEBUG
!!
!!END SUBROUTINE readonewave_linear

!> Reads wavefunction from file and transforms it properly if hgrid or size of simulation cell
!! have changed
subroutine readmywaves_linear_new(iproc,nproc,dir_output,filename,iformat,at,tmb,rxyz,&
       ref_frags,input_frag,frag_calc,kernel_restart,max_nbasis_env,frag_env_mapping,orblist)
  use module_base
  use module_types
  use yaml_output
  use module_fragments
  !use internal_io
  use module_interfaces, only: open_filename_of_iorb, reformat_supportfunctions
  use io, only: read_coeff_minbasis, io_read_descr_linear, read_psig, io_error, read_dense_matrix
  use locreg_operations, only: lpsi_to_global2
  use public_enums
  implicit none
  integer, intent(in) :: iproc, nproc
  integer, intent(in) :: iformat
  type(atoms_data), intent(in) :: at
  type(DFT_wavefunction), intent(inout) :: tmb
  real(gp), dimension(3,at%astruct%nat), intent(in) :: rxyz
  !real(gp), dimension(3,at%astruct%nat), intent(out) :: rxyz_old
  character(len=*), intent(in) :: dir_output, filename
  type(fragmentInputParameters), intent(in) :: input_frag
  type(system_fragment), dimension(input_frag%nfrag_ref), intent(inout) :: ref_frags
  logical, intent(in) :: frag_calc, kernel_restart
  integer, intent(in) :: max_nbasis_env
  integer, dimension(input_frag%nfrag,max_nbasis_env,3), intent(inout) :: frag_env_mapping
  integer, dimension(tmb%orbs%norb), intent(in), optional :: orblist
  !Local variables
  real(gp), parameter :: W_tol=1.e-3_gp !< wahba's tolerance
  integer :: ncount1,ncount_rate,ncount_max,ncount2
  integer :: iorb_out,ispinor,ilr,iorb_old
  integer :: confPotOrder,onwhichatom_tmp,unitwf,itoo_big
  real(gp) :: confPotprefac
!!$ real(gp), dimension(3) :: mol_centre, mol_centre_new
  real(kind=4) :: tr0,tr1
  real(kind=8) :: tel,eval
  character(len=256) :: error, full_filename
  logical :: lstat
  character(len=*), parameter :: subname='readmywaves_linear_new'
  ! to eventually be part of the fragment structure?
  integer :: ndim_old, iiorb, ifrag, ifrag_ref, isfat, iorbp, iforb, isforb, iiat, iat, i, np, c, minperm, iorb, ind
  integer :: iatt, iatf, ityp, num_env
  type(local_zone_descriptors) :: lzd_old
  real(wp), dimension(:), pointer :: psi_old
  type(phi_array), dimension(:), pointer :: phi_array_old
  type(fragment_transformation), dimension(:), pointer :: frag_trans_orb, frag_trans_frag
  integer, dimension(:), allocatable :: ipiv
  real(gp), dimension(:,:), allocatable :: rxyz_new, rxyz4_ref, rxyz4_new, rxyz_ref
  real(gp), dimension(:,:), allocatable :: rxyz_old !<this is read from the disk and not needed
  real(kind=gp), dimension(:), allocatable :: dist
  real(gp) :: max_shift, dtol

  logical :: skip, binary
  integer :: itmb, jtmb, jat
  !!$ integer :: ierr

  ! DEBUG
  ! character(len=12) :: orbname
  real(wp), dimension(:), allocatable :: gpsi

  call cpu_time(tr0)
  call system_clock(ncount1,ncount_rate,ncount_max)

  ! check file format
  if (iformat == WF_FORMAT_ETSF) then
     call f_err_throw('Linear scaling with ETSF writing not implemented yet')
  else if (iformat /= WF_FORMAT_BINARY .and. iformat /= WF_FORMAT_PLAIN) then
     call f_err_throw('Unknown wavefunction file format from filename.')
  end if

  rxyz_old=f_malloc([3,at%astruct%nat],id='rxyz_old')

  ! to be fixed
  if (present(orblist)) then
     call f_err_throw('orblist no longer functional in initialize_linear_from_file due to addition of fragment calculation')
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
  allocate(lzd_old%Llr(lzd_old%nlr))
  do ilr=1,lzd_old%nlr
     call nullify_locreg_descriptors(lzd_old%llr(ilr))
  end do

  ! has size of new orbs, will possibly point towards the same tmb multiple times
  allocate(phi_array_old(tmb%orbs%norbp))
  do iorbp=1,tmb%orbs%norbp
     nullify(phi_array_old(iorbp)%psig)
  end do

  !allocate(frag_trans_orb(tmb%orbs%norbp))

  unitwf=99
  isforb=0
  isfat=0
  call timing(iproc,'tmbrestart','ON')
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
              phi_array_old(iorbp)%psig = f_malloc_ptr((/ 0.to.Lzd_old%Llr(ilr)%d%n1 , 1.to.2 , &
                   0.to.Lzd_old%Llr(ilr)%d%n2 , 1.to.2 , &
                   0.to.Lzd_old%Llr(ilr)%d%n3 , 1.to.2 /),id='phi_array_old(iorbp)%psig')
              call timing(iproc,'tmbrestart','OF')

              !read phig directly
              call timing(iproc,'readtmbfiles','ON')
              call read_psig(unitwf, (iformat == WF_FORMAT_PLAIN), Lzd_old%Llr(ilr)%wfd%nvctr_c, Lzd_old%Llr(ilr)%wfd%nvctr_f, &
                   Lzd_old%Llr(ilr)%d%n1, Lzd_old%Llr(ilr)%d%n2, Lzd_old%Llr(ilr)%d%n3, phi_array_old(iorbp)%psig, lstat, error)
              if (.not. lstat) call io_error(trim(error))
              call timing(iproc,'readtmbfiles','OF')

              call timing(iproc,'tmbrestart','ON')

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
  itoo_big=0
  fragment_if: if (input_frag%nfrag>1) then
     ! Find fragment transformations for each fragment, then put in frag_trans array for each orb
     allocate(frag_trans_frag(input_frag%nfrag))

     isfat=0
     isforb=0
     do ifrag=1,input_frag%nfrag
        ! find reference fragment this corresponds to
        ifrag_ref=input_frag%frag_index(ifrag)

        ! check if we need this fragment transformation on this proc
        ! if this is an environment calculation we need mapping on all mpi, so easier to just calculate all transformations on all procs
        if (ref_frags(ifrag_ref)%astruct_env%nat/=0) then
           skip=.false.
        else
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
        end if

        if (skip) then
           isfat=isfat+ref_frags(ifrag_ref)%astruct_frg%nat
           isforb=isforb+ref_frags(ifrag_ref)%fbasis%forbs%norb
           cycle
        end if

        if (ref_frags(ifrag_ref)%astruct_env%nat==0) then
           rxyz_ref = f_malloc((/ 3, ref_frags(ifrag_ref)%astruct_frg%nat /),id='rxyz_ref')
           rxyz_new = f_malloc((/ 3, ref_frags(ifrag_ref)%astruct_frg%nat /),id='rxyz_new')

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

           call f_free(rxyz_ref)
           call f_free(rxyz_new)

        ! take into account environment coordinates
        else
           call match_environment_atoms(isfat,at,rxyz,tmb%orbs,ref_frags(ifrag_ref),&
                max_nbasis_env,frag_env_mapping(ifrag,:,:),frag_trans_frag(ifrag),.false.)
        end if

        ! in environment case we're calculating all transformations on each MPI, so no need to incrememnt on each
        if (frag_trans_frag(ifrag)%Werror > W_tol .and. ((ref_frags(ifrag_ref)%astruct_env%nat/=0 .and. iproc==0) &
             .or. ref_frags(ifrag_ref)%astruct_env%nat==0))then
           call f_increment(itoo_big)
        end if

        ! useful for identifying which fragments are problematic
        if (iproc==0 .and. frag_trans_frag(ifrag)%Werror>W_tol) then
           write(*,'(A,1x,I3,1x,I3,1x,3(F12.6,1x),2(F12.6,1x),2(I8,1x))') 'ifrag,ifrag_ref,rot_axis,theta,error',&
                ifrag,ifrag_ref,frag_trans_frag(ifrag)%rot_axis,frag_trans_frag(ifrag)%theta/(4.0_gp*atan(1.d0)/180.0_gp),&
                frag_trans_frag(ifrag)%Werror,itoo_big,iproc
           write(*,*) ''
        end if

        isfat=isfat+ref_frags(ifrag_ref)%astruct_frg%nat
        isforb=isforb+ref_frags(ifrag_ref)%fbasis%forbs%norb
     end do

     !if (bigdft_mpi%nproc > 1) then
     !   call mpiallred(frag_env_mapping, mpi_sum, comm=bigdft_mpi%mpi_comm)
     !end if

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

              frag_trans_orb(iorbp)%rot_center=frag_trans_frag(ifrag)%rot_center
              frag_trans_orb(iorbp)%rot_center_new=frag_trans_frag(ifrag)%rot_center_new

              !!!!!!!!!!!!!
              iiat=tmb%orbs%onwhichatom(iiorb)

              ! use atom position
              frag_trans_orb(iorbp)%rot_center=rxyz_old(:,iiat)
              frag_trans_orb(iorbp)%rot_center_new=rxyz(:,iiat)
              !!!!!!!!!!!!!

              frag_trans_orb(iorbp)%rot_axis=(frag_trans_frag(ifrag)%rot_axis)
              frag_trans_orb(iorbp)%theta=frag_trans_frag(ifrag)%theta
              frag_trans_orb(iorbp)%Rmat=frag_trans_frag(ifrag)%Rmat

              !write(*,'(a,x,2(i2,x),4(f5.2,x),6(f7.3,x))'),'trans2',ifrag,iiorb,frag_trans_orb(iorbp)%theta,&
              !     frag_trans_orb(iorbp)%rot_axis, &
              !     frag_trans_orb(iorbp)%rot_center,frag_trans_orb(iorbp)%rot_center_new
           end do
        end do
        isforb=isforb+ref_frags(ifrag_ref)%fbasis%forbs%norb
        isfat=isfat+ref_frags(ifrag_ref)%astruct_frg%nat
     end do

     deallocate(frag_trans_frag)
  else
     ! only 1 'fragment', calculate rotation/shift atom wise, using nearest neighbours
     allocate(frag_trans_orb(tmb%orbs%norbp))

     rxyz4_ref = f_malloc((/ 3, min(4, ref_frags(ifrag_ref)%astruct_frg%nat) /),id='rxyz4_ref')
     rxyz4_new = f_malloc((/ 3, min(4, ref_frags(ifrag_ref)%astruct_frg%nat) /),id='rxyz4_new')

     isforb=0
     isfat=0
     do ifrag=1,input_frag%nfrag
        ! find reference fragment this corresponds to
        ifrag_ref=input_frag%frag_index(ifrag)

        rxyz_ref = f_malloc((/ 3, ref_frags(ifrag_ref)%astruct_frg%nat /),id='rxyz_ref')
        rxyz_new = f_malloc((/ 3, ref_frags(ifrag_ref)%astruct_frg%nat /),id='rxyz_new')
        dist = f_malloc(ref_frags(ifrag_ref)%astruct_frg%nat,id='dist')
        ipiv = f_malloc(ref_frags(ifrag_ref)%astruct_frg%nat,id='ipiv')

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

              call find_frag_trans(min(4,ref_frags(ifrag_ref)%astruct_frg%nat),rxyz4_ref,rxyz4_new,&
                   frag_trans_orb(iorbp))
              if (frag_trans_orb(iorbp)%Werror > W_tol) call f_increment(itoo_big)

!!$              print *,'transformation of the fragment, iforb',iforb
!!$              write(*,'(A,I3,1x,I3,1x,3(F12.6,1x),F12.6)') 'ifrag,iorb,rot_axis,theta',&
!!$                   ifrag,iiorb,frag_trans_orb(iorbp)%rot_axis,frag_trans_orb(iorbp)%theta/(4.0_gp*atan(1.d0)/180.0_gp)

           end do
        end do
        isforb=isforb+ref_frags(ifrag_ref)%fbasis%forbs%norb
        isfat=isfat+ref_frags(ifrag_ref)%astruct_frg%nat

        call f_free(ipiv)
        call f_free(dist)
        call f_free(rxyz_ref)
        call f_free(rxyz_new)

     end do

     call f_free(rxyz4_ref)
     call f_free(rxyz4_new)

  end if fragment_if

  !reduce the number of warnings
  if (nproc >1) call mpiallred(itoo_big,1,op=MPI_SUM,comm=bigdft_mpi%mpi_comm)

  if (itoo_big > 0 .and. iproc==0) call yaml_warning('Found '//itoo_big//' warning of high Wahba cost functions')


  !!!debug - check calculated transformations
  !!if (iproc==0) then
  !!   open(109)
  !!   do ifrag=1,input_frag%nfrag
  !!      ! find reference fragment this corresponds to
  !!      ifrag_ref=input_frag%frag_index(ifrag)
  !!      num_env=ref_frags(ifrag_ref)%astruct_env%nat-ref_frags(ifrag_ref)%astruct_frg%nat
  !!      write(109,'((a,1x,I5,1x),F12.6,2x,3(F12.6,1x),6(1x,F18.6),2x,F6.2,2(2x,I5))') &
  !!           trim(input_frag%label(ifrag_ref)),ifrag,&
  !!           frag_trans_frag(ifrag)%theta/(4.0_gp*atan(1.d0)/180.0_gp),frag_trans_frag(ifrag)%rot_axis,&
  !!           frag_trans_frag(ifrag)%rot_center,frag_trans_frag(ifrag)%rot_center_new,&
  !!           Werror(ifrag),num_env,ifrag_ref
  !!   end do
  !!   close(109)
  !!end if

  call timing(iproc,'tmbrestart','OF')
  call reformat_supportfunctions(iproc,nproc,&
       at,rxyz_old,rxyz,.false.,tmb,ndim_old,lzd_old,frag_trans_orb,&
       psi_old,trim(dir_output),input_frag,ref_frags,max_shift,phi_array_old)
  call timing(iproc,'tmbrestart','ON')

  deallocate(frag_trans_orb)

  do iorbp=1,tmb%orbs%norbp
     !nullify/deallocate here as appropriate, in future may keep
     call f_free_ptr(phi_array_old(iorbp)%psig)
  end do

  deallocate(phi_array_old)
  call f_free(rxyz_old)
  call deallocate_local_zone_descriptors(lzd_old)

  !! DEBUG - plot in global box - CHECK WITH REFORMAT ETC IN LRs
  !ind=1
  !gpsi=f_malloc(tmb%Lzd%glr%wfd%nvctr_c+7*tmb%Lzd%glr%wfd%nvctr_f,id='gpsi')
  !do iorbp=1,tmb%orbs%norbp
  !   iiorb=iorbp+tmb%orbs%isorb
  !   ilr = tmb%orbs%inwhichlocreg(iiorb)
  !   call f_zero(gpsi)
  !   call Lpsi_to_global2(iproc, tmb%Lzd%Llr(ilr)%wfd%nvctr_c+7*tmb%Lzd%Llr(ilr)%wfd%nvctr_f, &
  !        tmb%Lzd%glr%wfd%nvctr_c+7*tmb%Lzd%glr%wfd%nvctr_f, &
  !        1, 1, 1, tmb%Lzd%glr, tmb%Lzd%Llr(ilr), &
  !        tmb%psi(ind:ind+tmb%Lzd%Llr(ilr)%wfd%nvctr_c+7*tmb%Lzd%Llr(ilr)%wfd%nvctr_f), gpsi)
  !   call plot_wf(.false.,trim(dir_output)//trim(adjustl(yaml_toa(iiorb))),1,at,1.0_dp,tmb%Lzd%glr,&
  !        tmb%Lzd%hgrids(1),tmb%Lzd%hgrids(2),tmb%Lzd%hgrids(3),rxyz,gpsi)
  !   ind = ind + tmb%Lzd%Llr(ilr)%wfd%nvctr_c+7*tmb%Lzd%Llr(ilr)%wfd%nvctr_f
  !end do
  !call f_free(gpsi)
  !! END DEBUG

  ! Read the coefficient file for each fragment and assemble total coeffs
  ! coeffs should eventually go into ref_frag array and then point? or be copied to (probably copied as will deallocate frag)
  unitwf=99
  isforb=0
  ! can directly loop over reference fragments here
  !do ifrag=1,input_frag%nfrag
  do ifrag_ref=1,input_frag%nfrag_ref
     ! find reference fragment this corresponds to
     !ifrag_ref=input_frag%frag_index(ifrag)

     ! read coeffs/kernel
     if (kernel_restart) then
        full_filename=trim(dir_output)//trim(input_frag%dirname(ifrag_ref))//'density_kernel.bin'
        !should fragments have some knowledge of spin?
        !assume kernel is in binary if tmbs are...
        binary=(iformat == WF_FORMAT_BINARY)
        call read_dense_matrix(full_filename, binary, tmb%orbs%nspinor, ref_frags(ifrag_ref)%fbasis%forbs%norb, &
             ref_frags(ifrag_ref)%kernel, ref_frags(ifrag_ref)%astruct_frg%nat)


        if (ref_frags(ifrag_ref)%nbasis_env/=0) then
           full_filename=trim(dir_output)//trim(input_frag%dirname(ifrag_ref))//'density_kernel_env.bin'
           !should fragments have some knowledge of spin?
           !assume kernel is in binary if tmbs are...
           binary=(iformat == WF_FORMAT_BINARY)
           call read_dense_matrix(full_filename, binary, tmb%orbs%nspinor, ref_frags(ifrag_ref)%nbasis_env, &
                ref_frags(ifrag_ref)%kernel_env, ref_frags(ifrag_ref)%astruct_env%nat)
        end if

        !!if (iproc==0) then
        !!   open(32)
        !!   do itmb=1,tmb%orbs%norb
        !!      do jtmb=1,tmb%orbs%norb
        !!         write(32,*) itmb,jtmb,tmb%coeff(itmb,jtmb),ref_frags(ifrag_ref)%kernel(itmb,jtmb,1)
        !!      end do
        !!   end do
        !!   write(32,*) ''
        !!   close(32)
        !!end if

     else
        full_filename=trim(dir_output)//trim(input_frag%dirname(ifrag_ref))//trim(filename)//'_coeff.bin'
        call f_open_file(unitwf,file=trim(full_filename),binary=iformat == WF_FORMAT_BINARY)
        call read_coeff_minbasis(unitwf,(iformat == WF_FORMAT_PLAIN),iproc,ref_frags(ifrag_ref)%fbasis%forbs%norb,&
             ref_frags(ifrag_ref)%nelec,ref_frags(ifrag_ref)%fbasis%forbs%norb,ref_frags(ifrag_ref)%coeff,ref_frags(ifrag_ref)%eval)
        call f_close(unitwf)
     end if

     !isforb=isforb+ref_frags(ifrag_ref)%fbasis%forbs%norb
  end do


  call cpu_time(tr1)
  call system_clock(ncount2,ncount_rate,ncount_max)
  tel=dble(ncount2-ncount1)/dble(ncount_rate)

  if (iproc == 0) then
     call yaml_sequence_open('Reading Waves Time')
     call yaml_sequence(advance='no')
     call yaml_mapping_open(flow=.true.)
     call yaml_map('Process',iproc)
     call yaml_map('Timing',(/ real(tr1-tr0,kind=8),tel /),fmt='(1pe10.3)')
     call yaml_mapping_close()
     call yaml_sequence_close()
  end if
  !write(*,'(a,i4,2(1x,1pe10.3))') '- READING WAVES TIME',iproc,tr1-tr0,tel
  call timing(iproc,'tmbrestart','OF')


END SUBROUTINE readmywaves_linear_new


   !> Matches neighbouring atoms from environment file to those in full system
   !! returns atom mapping information and 'best' fragment transformation and corresponding Wahba error
   subroutine match_environment_atoms(isfat,at,rxyz,orbs,ref_frag,max_nbasis_env,frag_env_mapping,frag_trans,ignore_species)
      use module_base
      use module_types
      use module_fragments
      use io, only: dist_and_shift
      implicit none
      integer, intent(in) :: isfat
      type(atoms_data), intent(in) :: at
      real(gp), dimension(3,at%astruct%nat), intent(in) :: rxyz
      type(orbitals_data), intent(in) :: orbs
      type(system_fragment), intent(in) :: ref_frag
      integer, intent(in) :: max_nbasis_env
      integer, dimension(max_nbasis_env,3), intent(out) :: frag_env_mapping
      type(fragment_transformation), intent(out) :: frag_trans
      logical, intent(in) :: ignore_species

      !local variables
      integer :: iat, ityp, ipiv_shift, iatt, iatf, num_env, np, c, itmb, iorb, i, minperm, num_neighbours
      integer, allocatable, dimension(:) :: ipiv, array_tmp, num_neighbours_type
      integer, dimension(:,:), allocatable :: permutations

      real(kind=gp) :: minerror, mintheta, err_tol, rot_tol
      real(kind=gp), allocatable, dimension(:) :: dist
      real(kind=gp), dimension(:,:), allocatable :: rxyz_new_all, rxyz_frg_new, rxyz_new_trial, rxyz_ref, rxyz_new

      logical :: perx, pery, perz, wrong_atom

      character(len=2) :: atom_ref, atom_trial


      !from _env file - includes fragment and environment
      rxyz_ref = f_malloc((/ 3,ref_frag%astruct_env%nat /),id='rxyz_ref')
      !all coordinates in new system, except those in fragment
      rxyz_new_all = f_malloc((/ 3,at%astruct%nat-ref_frag%astruct_frg%nat /),id='rxyz_new_all')
      dist = f_malloc(at%astruct%nat-ref_frag%astruct_frg%nat,id='dist')
      ipiv = f_malloc(at%astruct%nat-ref_frag%astruct_frg%nat,id='ipiv')
      !just the fragment in the new system
      rxyz_frg_new = f_malloc((/ 3,ref_frag%astruct_frg%nat /),id='rxyz_frg_new')

      !just a straightforward copy
      do iat=1,ref_frag%astruct_env%nat
         rxyz_ref(:,iat)=ref_frag%astruct_env%rxyz(:,iat)
      end do

      !also add up how many atoms of each type (don't include fragment itself)
      num_neighbours_type=f_malloc0(at%astruct%ntypes,id='num_neighbours_type')
      do iat=ref_frag%astruct_frg%nat+1,ref_frag%astruct_env%nat
         !be careful here as atom types not necessarily in same order in env file as in main file (or even same number thereof)
         do ityp=1,at%astruct%ntypes
            if (trim(at%astruct%atomnames(ityp))==trim(ref_frag%astruct_env%atomnames(ref_frag%astruct_env%iatype(iat)))) then
               num_neighbours_type(ityp) = num_neighbours_type(ityp)+1
               exit
            end if
         end do
      end do
      num_neighbours=sum(num_neighbours_type)
      if (sum(num_neighbours_type)/=ref_frag%astruct_env%nat-ref_frag%astruct_frg%nat) &
           stop 'Error with num_neighbours_type in fragment environment restart'
      !print*,ifrag,ifrag_ref,sum(num_neighbours_type),num_neighbours_type

      !take all atoms not in this fragment (might be overcomplicating this...)
      do iat=1,isfat
         rxyz_new_all(:,iat)=rxyz(:,iat)
      end do

      do iat=isfat+ref_frag%astruct_frg%nat+1,at%astruct%nat
         rxyz_new_all(:,iat-ref_frag%astruct_frg%nat)=rxyz(:,iat)
      end do

      !just take those in the fragment
      do iat=1,ref_frag%astruct_frg%nat
         rxyz_frg_new(:,iat)=rxyz(:,isfat+iat)
      end do

      !this should just be the fragment centre - fragment xyz comes first in rxyz_env so this should be ok
      frag_trans%rot_center=frag_center(ref_frag%astruct_frg%nat,rxyz_ref(:,1:ref_frag%astruct_frg%nat))
      !the fragment centre in new coordinates
      frag_trans%rot_center_new=frag_center(ref_frag%astruct_frg%nat,rxyz_frg_new)

      ! shift rxyz wrt center of rotation
      do iat=1,ref_frag%astruct_env%nat
         rxyz_ref(:,iat)=rxyz_ref(:,iat)-frag_trans%rot_center
      end do

      ! find distances from this atom BEFORE shifting
      perx=(at%astruct%geocode /= 'F')
      pery=(at%astruct%geocode == 'P')
      perz=(at%astruct%geocode /= 'F')

      !if coordinates wrap around (in periodic), correct before shifting
      !assume that the fragment itself doesn't, just the environment...
      !think about other periodic cases that might need fixing...
      do iat=1,at%astruct%nat-ref_frag%astruct_frg%nat
         dist(iat) = dist_and_shift(perx,at%astruct%cell_dim(1),frag_trans%rot_center_new(1),rxyz_new_all(1,iat))**2
         dist(iat) = dist(iat) + dist_and_shift(pery,at%astruct%cell_dim(2),frag_trans%rot_center_new(2),rxyz_new_all(2,iat))**2
         dist(iat) = dist(iat) + dist_and_shift(perz,at%astruct%cell_dim(3),frag_trans%rot_center_new(3),rxyz_new_all(3,iat))**2
         dist(iat) = -dsqrt(dist(iat))
!!$         write(*,'(A,2(I3,2x),F12.6,3x,2(3(F12.6,1x),2x))') 'ifrag,iat,dist',ifrag,iat,dist(iat),&
!!$              at%astruct%cell_dim(:),rxyz_new_all(:,iat)

         rxyz_new_all(:,iat) = rxyz_new_all(:,iat)-frag_trans%rot_center_new
      end do

      do iat=1,ref_frag%astruct_frg%nat
         rxyz_frg_new(:,iat)=rxyz_frg_new(:,iat)-frag_trans%rot_center_new
      end do

      ! sort atoms into neighbour order
      call sort_positions(at%astruct%nat-ref_frag%astruct_frg%nat,dist,ipiv)

      rxyz_new = f_malloc((/ 3,ref_frag%astruct_env%nat /),id='rxyz_new')

      ! take fragment and closest neighbours (assume that environment atoms were originally the closest)
      ! put mapping in column 2 to avoid overwriting later
      do iat=1,ref_frag%astruct_frg%nat
         rxyz_new(:,iat)=rxyz_frg_new(:,iat)
         frag_env_mapping(iat,2) = isfat+iat
      end do

      iatf=0
      do ityp=1,at%astruct%ntypes
         iatt=0
         if (num_neighbours_type(ityp)==0 .and. (.not. ignore_species)) cycle
         do iat=1,at%astruct%nat-ref_frag%astruct_frg%nat
            !ipiv_shift needed for quantities which reference full rxyz, not rxyz_new_all which has already eliminated frag atoms
            if (ipiv(iat)<=isfat) then
               ipiv_shift=ipiv(iat)
            else
               ipiv_shift=ipiv(iat)+ref_frag%astruct_frg%nat
            end if
            if (at%astruct%iatype(ipiv_shift)/=ityp .and. (.not. ignore_species)) cycle
            iatf=iatf+1
            iatt=iatt+1
            rxyz_new(:,iatf+ref_frag%astruct_frg%nat)=rxyz_new_all(:,ipiv(iat))
            frag_env_mapping(iatf+ref_frag%astruct_frg%nat,2) = ipiv_shift
            if ((ignore_species.and.iatt==num_neighbours) &
                 .or. ((.not. ignore_species) .and. iatt==num_neighbours_type(ityp))) exit
         end do
         !write(*,'(a,7(i3,2x))')'ityp',ityp,at%astruct%ntypes,iatt,iatf,num_neighbours_type(ityp),ifrag,ifrag_ref
         if (ignore_species) exit
      end do

      !print*,'iatf,sum',iatf,sum(num_neighbours_type),num_neighbours,ignore_species
      if (((.not. ignore_species) .and. iatf/=sum(num_neighbours_type)) &
           .or. (ignore_species .and. iatf/=num_neighbours)) stop 'Error num_neighbours/=iatf in match_environment_atoms'

      call f_free(num_neighbours_type)
      call f_free(dist)
      call f_free(ipiv)
      call f_free(rxyz_frg_new)
      call f_free(rxyz_new_all)

      !# don't sort rxyz_ref - just check Wahba permutations for all atoms
      !# assuming small number of neighbours so saves generalizing things and makes it easier for mapping env -> full

!!$      do iat=1,ref_frag%astruct_env%nat
!!$         write(*,'(A,3(I3,2x),2x,2(3(F12.6,1x),2x))') 'ifrag,ifrag_ref,iat,rxyz_ref,rxyz_new',&
!!$              ifrag,ifrag_ref,iat,rxyz_ref_sorted(:,iat),rxyz_new(:,iat)
!!$      end do

      !ADD CHECKING OF ATOM TYPE TO ABOVE SORTING PROCEDURE, for the moment assuming identical atom types
      !if error is above some threshold and we have some degenerate distances
      !then try to find ordering which gives lowest Wahba error
      !also give preference to zero rotation
      !write(*,'(A)') 'Problem matching environment atoms to new coordinates, attempting to find correct order'
      !write(*,'(A)') 'Checking for ordering giving a more accurate transformation/no rotation'

      num_env=ref_frag%astruct_env%nat-ref_frag%astruct_frg%nat

      !assume that we have only a small number of identical distances, or this would become expensive...
      array_tmp=f_malloc(num_env,id='array_tmp')
      do i=1,num_env
         array_tmp(i)=i
      end do

      np=fact(num_env)
      permutations=f_malloc((/num_env,np/),id='permutations')
      c=0
      call reorder(num_env,num_env,c,np,array_tmp,permutations)
      call f_free(array_tmp)

      rxyz_new_trial = f_malloc((/ 3,ref_frag%astruct_env%nat /),id='rxyz_new_trial')

      !the fragment part doesn't change
      do iat=1,ref_frag%astruct_frg%nat
         rxyz_new_trial(:,iat) = rxyz_new(:,iat)
      end do

      minerror=1.d100
      minperm=-1
      mintheta=-1

      !test each permutation
      do i=1,np
         wrong_atom=.false.
         !first check that the atom types are coherent - if not reject this transformation
         do iat=ref_frag%astruct_frg%nat+1,ref_frag%astruct_env%nat
            atom_ref = trim(ref_frag%astruct_env%atomnames(ref_frag%astruct_env%iatype(iat)))
            atom_trial = trim(at%astruct%atomnames(at%astruct%iatype&
                 (frag_env_mapping(permutations(iat-ref_frag%astruct_frg%nat,i)+ref_frag%astruct_frg%nat,2))))
            !write(*,'(a,4(i3,2x),2(a2,2x),3(i3,2x))') 'ifrag,ifrag_ref,i,iat,atom_ref,atom_trial',ifrag,ifrag_ref,i,iat,&
            !     trim(atom_ref),trim(atom_trial),&
            !      frag_env_mapping(iat,2),permutations(iat-ref_frag%astruct_frg%nat,i),&
            !      frag_env_mapping(permutations(iat-ref_frag%astruct_frg%nat,i)+ref_frag%astruct_frg%nat,2)
            if (trim(atom_ref)/=trim(atom_trial) .and. (.not. ignore_species)) then
               wrong_atom=.true.
               exit
            end if
         end do
         if (wrong_atom) cycle

         do iat=ref_frag%astruct_frg%nat+1,ref_frag%astruct_env%nat
           rxyz_new_trial(:,iat) &
                 = rxyz_new(:,ref_frag%astruct_frg%nat &
                   + permutations(iat-ref_frag%astruct_frg%nat,i))
         end do
         call find_frag_trans(ref_frag%astruct_env%nat,rxyz_ref,&
              rxyz_new_trial,frag_trans)
         !if (frag_trans%Werror > W_tol) call f_increment(itoo_big)

         !do iat=1,ref_frag%astruct_env%nat
         !   write(*,'(A,3(I3,2x),2x,2(3(F12.6,1x),2x))') 'ifrag,ifrag_ref,iat,rxyz_new,rxyz_ref',&
         !        ifrag,ifrag_ref,iat,rxyz_new_trial(:,iat),rxyz_ref(:,iat)
         !end do
         !write(*,'(A,I3,2x,3(I3,1x),1x,F12.6)') 'i,perms,error: ',i,permutations(:,i),frag_trans%Werror
         !prioritize no rotation, and if not possible 180 degrees
         !could improve logic/efficiency here, i.e. stop checking once below some threshold
         !if ((frag_trans%Werror < minerror .and. (mintheta/=0 .or. minerror-frag_trans%Werror>1e-6)) &
         !     .or. (frag_trans%Werror-minerror<1e-6.and.frag_trans%theta==0.0d0) then

         err_tol = 1e-3 !1e-6
         rot_tol = 1e-3 !1e-6
         ! less than minerror by more than some tol
         ! or ~same error and zero rotation (wrt tol)
         ! or ~same error, 180 rotation (wrt tol) and not already zero rotation
         if ( (frag_trans%Werror < minerror - err_tol) &
              .or. (abs(frag_trans%Werror - minerror) < err_tol .and. abs(frag_trans%theta - 0.0d0) < rot_tol)) then ! &
!              .or. (abs(frag_trans%Werror - minerror) < err_tol .and. mintheta /= 0.0d0 &
!                   .and. abs(frag_trans%theta - 4.0_gp*atan(1.d0)) < rot_tol) ) then

            mintheta = frag_trans%theta
            minerror = frag_trans%Werror
            minperm = i
         end if
      end do

      ! use this as final transformation
      if (minperm/=-1) then
         !LG: commented it out, maybe it might be useful for debugging
         !write(*,'(A,I3,2x,2(F12.6,2x))') 'Final value of cost function:',&
         !     minperm,minerror,mintheta/(4.0_gp*atan(1.d0)/180.0_gp)
         do iat=ref_frag%astruct_frg%nat+1,ref_frag%astruct_env%nat
            rxyz_new_trial(:,iat) = rxyz_new(:,ref_frag%astruct_frg%nat + permutations(iat-ref_frag%astruct_frg%nat,minperm))
         end do
         call find_frag_trans(ref_frag%astruct_env%nat,rxyz_ref,rxyz_new_trial,frag_trans)
         !if (frag_trans%Werror > W_tol) call f_increment(itoo_big)

         do iat=1,ref_frag%astruct_frg%nat
            frag_env_mapping(iat,3) = frag_env_mapping(iat,2)
         end do
         do iat=ref_frag%astruct_frg%nat+1,ref_frag%astruct_env%nat
            frag_env_mapping(iat,3) = frag_env_mapping(ref_frag%astruct_frg%nat &
                 + permutations(iat-ref_frag%astruct_frg%nat,minperm),2)
         end do

         ! fill in 1st and 2nd columns of env_mapping
         itmb = 0
         do iat=1,ref_frag%astruct_env%nat
            do iorb=1,orbs%norb
               if (orbs%onwhichatom(iorb) == frag_env_mapping(iat,3)) then
                  itmb = itmb+1
                  frag_env_mapping(itmb,1) = iorb
                  frag_env_mapping(itmb,2) = iat
               end if
           end do
         end do
         if (itmb /= ref_frag%nbasis_env) stop 'Error with nbasis_env'

         !do iorb=1,ref_frag%nbasis_env
         !   write(*,'(A,5(1x,I4))') 'mapping: ',ifrag,ifrag_ref,frag_env_mapping(iorb,:)
         !end do

         !debug
         !do iat=1,ref_frag%astruct_env%nat
         !   write(*,'(a,4(i3,2x),a,2(3(f8.2,1x),4x))') 'if,ifr,ia,m,t',ifrag,ifrag_ref,iat,frag_env_mapping(iat,3),&
         !        trim(at%astruct%atomnames(at%astruct%iatype(frag_env_mapping(iat,3)))),rxyz_new_trial(:,iat),rxyz_ref(:,iat)
         !end do

      else
         stop 'Error finding environment transformation'
      end if

      call f_free(rxyz_new_trial)
      call f_free(permutations)
      call f_free(rxyz_ref)
      call f_free(rxyz_new)

   contains

      recursive subroutine reorder(nf,n,c,np,array_in,permutations)
        implicit none

        integer, intent(in) :: nf,n,np
        integer, intent(inout) :: c
        integer, dimension(1:nf), intent(inout) :: array_in
        integer, dimension(1:nf,1:np), intent(inout) :: permutations

        integer :: i, tmp
        integer, dimension(1:nf) :: array_out

        if (n>1) then
           do i=n,1,-1
              array_out=array_in
              tmp=array_in(n)
              array_out(n)=array_in(i)
              array_out(i)=tmp
              !print*,'i',i,n,'in',array_in,'out',array_out
              call reorder(nf,n-1,c,np,array_out,permutations)
           end do
        else
           c=c+1
           !print*,c,array_in
           permutations(:,c)=array_in(:)
           return
        end if

      end subroutine reorder

      function fact(n)
        implicit none

        integer, intent(in) :: n
        integer :: fact

        integer :: i

        fact=1
        do i=1,n
           fact = fact * i
        end do

      end function

   end subroutine match_environment_atoms


!> Initializes onwhichatom, inwhichlocreg, locrad and locregcenter from file
subroutine initialize_linear_from_file(iproc,nproc,input_frag,astruct,rxyz,orbs,Lzd,iformat,&
     dir_output,filename,ref_frags,orblist)
  use module_base
  use module_types
  use module_defs
  use yaml_output
  use module_fragments
  use module_interfaces, only: open_filename_of_iorb
  use locregs, only: locreg_null
  use io, only: io_read_descr_linear
  use public_enums
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
  integer :: ilr, iorb_old, iorb, ispinor, iorb_out, iforb, isforb, isfat, iiorb, iorbp, ifrag, ifrag_ref
  integer, dimension(3) :: n_old, ns_old
  integer :: confPotOrder, iat,unitwf
  real(gp), dimension(3) :: hgrids_old
  real(kind=8) :: eval, confPotprefac
  real(gp), dimension(orbs%norb):: locrad
  real(gp), dimension(3) :: locregCenter
  real(kind=8), dimension(:), allocatable :: lrad
  real(gp), dimension(:,:), allocatable :: cxyz
  character(len=256) :: full_filename


  unitwf=99
  ! to be fixed
  if (present(orblist)) then
     stop 'orblist no longer functional in initialize_linear_from_file due to addition of fragment calculation'
  end if

  ! NOTES:
  ! The orbs%norb family must be all constructed before this routine
  ! This can be done from the input.lin since the number of basis functions should be fixed.
  ! Fragment structure must also be fully initialized, even if this is not a fragment calculation

  call f_zero(locrad)
  call f_zero(orbs%norb,orbs%onwhichatom(1))

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

                 call open_filename_of_iorb(unitwf,(iformat == WF_FORMAT_BINARY),full_filename, &
                      & orbs,iorbp,ispinor,iorb_out,iforb)
                      !& ref_frags(ifrag_ref)%fbasis%forbs,iforb,ispinor,iorb_out)

                 !print *,'before crash',iorbp,trim(full_filename),iiorb,iforb

                 call io_read_descr_linear(unitwf,(iformat == WF_FORMAT_PLAIN), iorb_old, eval, n_old(1), n_old(2), n_old(3), &
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
                 call f_close(unitwf)

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
  if (nproc > 1)  call mpiallred(orbs%onwhichatom(1),orbs%norb,MPI_SUM,comm=bigdft_mpi%mpi_comm)
  if (nproc > 1)  call mpiallred(locrad,MPI_SUM,comm=bigdft_mpi%mpi_comm)

  cxyz = f_malloc((/ 3, Lzd%nlr /),id='cxyz')
  lrad = f_malloc(Lzd%nlr,id='lrad')

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
  allocate(lzd%Llr(lzd%nlr))
  do ilr=1,lzd%nlr
     lzd%Llr(ilr)=locreg_null()
  end do
  do ilr=1,lzd%nlr
      lzd%llr(ilr)%locrad=lrad(ilr)
      lzd%llr(ilr)%locregCenter=cxyz(:,ilr)
  end do

  call f_free(cxyz)
  call f_free(lrad)

END SUBROUTINE initialize_linear_from_file


!> Copy old support functions from phi to phi_old
subroutine copy_old_supportfunctions(iproc,orbs,lzd,phi,lzd_old,phi_old)
  use module_base
  use module_types
  use locregs, only: copy_locreg_descriptors
  use yaml_output
  implicit none
  integer,intent(in) :: iproc
  type(orbitals_data), intent(in) :: orbs
  type(local_zone_descriptors), intent(in) :: lzd
  type(local_zone_descriptors), intent(inout) :: lzd_old
  real(wp), dimension(:), pointer :: phi,phi_old
  !Local variables
  character(len=*), parameter :: subname='copy_old_supportfunctions'
  integer :: iseg,j,ind1,iorb,ii,iiorb,ilr
  real(kind=8) :: tt

  ! First copy global quantities
  call nullify_locreg_descriptors(lzd_old%glr)

  call copy_locreg_descriptors(lzd%glr,lzd_old%glr)

!!$  lzd_old%glr%wfd%nvctr_c = lzd%glr%wfd%nvctr_c
!!$  lzd_old%glr%wfd%nvctr_f = lzd%glr%wfd%nvctr_f
!!$  lzd_old%glr%wfd%nseg_c  = lzd%glr%wfd%nseg_c
!!$  lzd_old%glr%wfd%nseg_f  = lzd%glr%wfd%nseg_f

!!$  !allocations
!!$  call allocate_wfd(lzd_old%glr%wfd)
!!$
!!$  do iseg=1,lzd_old%glr%wfd%nseg_c+lzd_old%glr%wfd%nseg_f
!!$     lzd_old%glr%wfd%keyglob(1,iseg)    = lzd%glr%wfd%keyglob(1,iseg)
!!$     lzd_old%glr%wfd%keyglob(2,iseg)    = lzd%glr%wfd%keyglob(2,iseg)
!!$     lzd_old%glr%wfd%keygloc(1,iseg)    = lzd%glr%wfd%keygloc(1,iseg)
!!$     lzd_old%glr%wfd%keygloc(2,iseg)    = lzd%glr%wfd%keygloc(2,iseg)
!!$     lzd_old%glr%wfd%keyvloc(iseg)      = lzd%glr%wfd%keyvloc(iseg)
!!$     lzd_old%glr%wfd%keyvglob(iseg)     = lzd%glr%wfd%keyvglob(iseg)
!!$  enddo
  !!!deallocation
  !!call deallocate_wfd(lzd%glr%wfd,subname)

  !!lzd_old%glr%d%n1 = lzd%glr%d%n1
  !!lzd_old%glr%d%n2 = lzd%glr%d%n2
  !!lzd_old%glr%d%n3 = lzd%glr%d%n3
!!$  call copy_grid_dimensions(lzd%glr%d, lzd_old%glr%d)


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
      call copy_locreg_descriptors(lzd%llr(ilr), lzd_old%llr(ilr))
      ii = ii + lzd_old%llr(ilr)%wfd%nvctr_c + 7*lzd_old%llr(ilr)%wfd%nvctr_f
  end do

  phi_old = f_malloc_ptr(ii,id='phi_old')

  ! Now copy the suport functions
  if (iproc==0) call yaml_map('Check the normalization of the support functions, tolerance',1.d-3,fmt='(1es12.4)')
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
         !write(*,*)'wrong phi_old',iiorb,tt
         call yaml_warning('support function, value:'//trim(yaml_toa(iiorb,fmt='(i6)'))//trim(yaml_toa(tt,fmt='(1es18.9)')))
         !stop
      end if
  end do
!  if (iproc==0) call yaml_mapping_close()

  !!!deallocation
  !!i_all=-product(shape(phi))*kind(phi)
  !!deallocate(phi,stat=i_stat)
  !!call memocc(i_stat,i_all,'phi',subname)

END SUBROUTINE copy_old_supportfunctions


subroutine copy_old_coefficients(norb_tmb, nfvctr, coeff, coeff_old)
  use module_base
  implicit none

  ! Calling arguments
  integer,intent(in):: norb_tmb, nfvctr
  real(8),dimension(:,:),pointer:: coeff, coeff_old

  ! Local variables
  character(len=*),parameter:: subname='copy_old_coefficients'
!  integer:: istat,iall

  coeff_old = f_malloc_ptr((/ nfvctr, norb_tmb /),id='coeff_old')

  call vcopy(nfvctr*norb_tmb, coeff(1,1), 1, coeff_old(1,1), 1)

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
  !integer :: istat, iall

  inwhichlocreg_old = f_malloc_ptr(norb_tmb,id='inwhichlocreg_old')
  call vcopy(norb_tmb, inwhichlocreg(1), 1, inwhichlocreg_old(1), 1)
  !!iall=-product(shape(inwhichlocreg))*kind(inwhichlocreg)
  !!deallocate(inwhichlocreg,stat=istat)
  !!call memocc(istat,iall,'inwhichlocreg',subname)


  onwhichatom_old = f_malloc_ptr(norb_tmb,id='onwhichatom_old')
  call vcopy(norb_tmb, onwhichatom(1), 1, onwhichatom_old(1), 1)
  !!iall=-product(shape(onwhichatom))*kind(onwhichatom)
  !!deallocate(onwhichatom,stat=istat)
  !!call memocc(istat,iall,'onwhichatom',subname)

END SUBROUTINE copy_old_inwhichlocreg


!> Reformat wavefunctions if the mesh have changed (in a restart)
!! NB add_derivatives must be false if we are using phi_array_old instead of psi_old and don't have the keys
subroutine reformat_supportfunctions(iproc,nproc,at,rxyz_old,rxyz,add_derivatives,tmb,ndim_old,lzd_old,&
       frag_trans,psi_old,input_dir,input_frag,ref_frags,max_shift,phi_array_old)
  use module_base
  use module_types
  use module_fragments
  use module_interfaces, only: reformat_one_supportfunction
  use yaml_output
  use bounds, only: ext_buffers
  implicit none
  integer, intent(in) :: iproc,nproc
  integer, intent(in) :: ndim_old
  type(atoms_data), intent(in) :: at
  real(gp), dimension(3,at%astruct%nat), intent(in) :: rxyz,rxyz_old
  type(DFT_wavefunction), intent(inout) :: tmb
  type(local_zone_descriptors), intent(inout) :: lzd_old
  type(fragment_transformation), dimension(tmb%orbs%norbp), intent(in) :: frag_trans
  real(wp), dimension(:), pointer :: psi_old
  type(phi_array), dimension(tmb%orbs%norbp), optional, intent(in) :: phi_array_old
  logical, intent(in) :: add_derivatives
  character(len=*), intent(in) :: input_dir
  type(fragmentInputParameters), intent(in) :: input_frag
  type(system_fragment), dimension(input_frag%nfrag_ref), intent(in) :: ref_frags
  real(gp),intent(out) :: max_shift
  !Local variables
  character(len=*), parameter :: subname='reformatmywaves'
  logical :: reformat
  integer :: iorb,j,jstart,jstart_old,iiorb,ilr,iiat
  integer:: idir,jstart_old_der,ncount,ilr_old
  !!integer :: i
  integer, dimension(3) :: ns_old,ns,n_old,n,nglr_old,nglr
  real(gp), dimension(3) :: centre_old_box,centre_new_box,da
  real(gp) :: tt,tol
  real(wp), dimension(:,:,:,:,:,:), pointer :: phigold
  real(wp), dimension(:), allocatable :: phi_old_der
  integer, dimension(0:7) :: reformat_reason
  character(len=12) :: orbname!, dummy
  real(wp), allocatable, dimension(:,:,:) :: psirold
  logical :: psirold_ok
  integer, dimension(3) :: nl, nr
  logical, dimension(3) :: per
  character(len=100) :: fragdir
  integer :: ifrag, ifrag_ref, iforb, isforb
  real(kind=gp), dimension(:,:,:), allocatable :: workarraytmp
  logical :: gperx, gpery, gperz, lperx, lpery, lperz, wrap_around
  integer :: gnbl1, gnbr1, gnbl2, gnbr2, gnbl3, gnbr3, lnbl1, lnbr1, lnbl2, lnbr2, lnbl3, lnbr3

  real(gp), external :: dnrm2
!  integer :: iat

  reformat_reason=0
  tol=1.d-3
  max_shift = 0.d0

  ! Get the derivatives of the support functions
  if (add_derivatives) then
     phi_old_der = f_malloc(3*ndim_old,id='phi_old_der')
     if (.not. associated(psi_old)) stop 'psi_old not associated in reformat_supportfunctions'
     call get_derivative_supportfunctions(ndim_old, lzd_old%hgrids(1), lzd_old, tmb%orbs, psi_old, phi_old_der)
     jstart_old_der=1
  end if

  nglr_old(1)=lzd_old%glr%d%n1
  nglr_old(2)=lzd_old%glr%d%n1
  nglr_old(3)=lzd_old%glr%d%n1
  nglr(1)=tmb%lzd%glr%d%n1
  nglr(2)=tmb%lzd%glr%d%n1
  nglr(3)=tmb%lzd%glr%d%n1

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
           n_old,n,ns_old,ns,nglr_old,nglr,at%astruct%geocode,& !lzd_old%llr(ilr)%geocode,&
           frag_trans(iorb),centre_old_box,centre_new_box,da,wrap_around)
      max_shift = max(max_shift,sqrt(da(1)**2+da(2)**2+da(3)**2))

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
             if (.not. wrap_around) then

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

            else

               ! THIS CASE NEEDS TESTING
               call psi_to_psig(n_old,lzd_old%llr(ilr)%wfd%nseg_c,lzd_old%llr(ilr)%wfd%nvctr_c,&
                    lzd_old%llr(ilr)%wfd%keygloc,lzd_old%llr(ilr)%wfd%keyvloc,&
                    lzd_old%llr(ilr)%wfd%nseg_f,lzd_old%llr(ilr)%wfd%nvctr_f,&
                    lzd_old%llr(ilr)%wfd%keygloc(1:,lzd_old%Llr(ilr)%wfd%nseg_c+1:), &
                    lzd_old%llr(ilr)%wfd%keyvloc(lzd_old%Llr(ilr)%wfd%nseg_c+1:), &
                    phigold,psi_old(jstart_old),psi_old(jstart_old+lzd_old%llr(ilr)%wfd%nvctr_c))

               call compress_plain(n(1),n(2),0,n(1),0,n(2),0,n(3), &
                    tmb%lzd%llr(ilr)%wfd%nseg_c,tmb%lzd%llr(ilr)%wfd%nvctr_c,tmb%lzd%llr(ilr)%wfd%keygloc(1,1), &
                    tmb%lzd%llr(ilr)%wfd%keyvloc(1),tmb%lzd%llr(ilr)%wfd%nseg_f,tmb%lzd%llr(ilr)%wfd%nvctr_f,&
                    tmb%lzd%llr(ilr)%wfd%keygloc(1,tmb%lzd%llr(ilr)%wfd%nseg_c+min(1,tmb%lzd%llr(ilr)%wfd%nseg_f)),&
                    tmb%lzd%llr(ilr)%wfd%keyvloc(tmb%lzd%llr(ilr)%wfd%nseg_c+min(1,tmb%lzd%llr(ilr)%wfd%nseg_f)),   &
                    phigold,tmb%psi(jstart),&
                    tmb%psi(jstart+tmb%lzd%llr(ilr)%wfd%nvctr_c+min(1,tmb%lzd%llr(ilr)%wfd%nvctr_f)-1))

               jstart_old=jstart_old+lzd_old%llr(ilr)%wfd%nvctr_c+7*lzd_old%llr(ilr)%wfd%nvctr_f
               jstart=jstart+tmb%lzd%llr(ilr)%wfd%nvctr_c+7*tmb%lzd%llr(ilr)%wfd%nvctr_f

            end if
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

!!$          print*, 'norm of read psi ',dnrm2(lzd_old%llr(ilr_old)%wfd%nvctr_c+7*lzd_old%llr(ilr_old)%wfd%nvctr_f,&
!!$               psi_old(jstart_old),1),&
!!$               lzd_old%llr(ilr_old)%wfd%nvctr_c,lzd_old%llr(ilr_old)%wfd%nvctr_f
!!$          print*, 'norm of read psic ',dnrm2(lzd_old%llr(ilr_old)%wfd%nvctr_c,psi_old(jstart_old),1)
!!$          print*, 'norm of read psif ',dnrm2(lzd_old%llr(ilr_old)%wfd%nvctr_f*7,&
!!$               psi_old(jstart_old-1+lzd_old%llr(ilr_old)%wfd%nvctr_c+min(1,lzd_old%llr(ilr_old)%wfd%nvctr_f)),1)


          ! uncompress or point towards correct phigold as necessary
          if (present(phi_array_old)) then
             phigold=>phi_array_old(iorb)%psig
          else
             phigold = f_malloc_ptr((/ 0.to.n_old(1), 1.to.2, 0.to.n_old(2), 1.to.2, 0.to.n_old(3), 1.to.2 /),id='phigold')

             call psi_to_psig(n_old,lzd_old%llr(ilr)%wfd%nseg_c,lzd_old%llr(ilr)%wfd%nvctr_c,&
                  lzd_old%llr(ilr)%wfd%keygloc,lzd_old%llr(ilr)%wfd%keyvloc,&
                  lzd_old%llr(ilr)%wfd%nseg_f,lzd_old%llr(ilr)%wfd%nvctr_f,&
                  lzd_old%llr(ilr)%wfd%keygloc(1:,lzd_old%Llr(ilr)%wfd%nseg_c+1:), &
                  lzd_old%llr(ilr)%wfd%keyvloc(lzd_old%Llr(ilr)%wfd%nseg_c+1:), &
                  phigold,psi_old(jstart_old),psi_old(jstart_old+lzd_old%llr(ilr)%wfd%nvctr_c))

             jstart_old=jstart_old+lzd_old%llr(ilr)%wfd%nvctr_c+7*lzd_old%llr(ilr)%wfd%nvctr_f

          end if

          !write(100+iproc,*) 'norm phigold ',dnrm2(8*(n1_old+1)*(n2_old+1)*(n3_old+1),phigold,1)
          !write(*,*) 'iproc,norm phigold ',iproc,dnrm2(8*product(n_old+1),phigold,1)

          ! read psir_old directly from files (don't have lzd_old to rebuild it)
          psirold_ok=.true.

          ! allow for fragment calculation
          if (input_frag%nfrag>1) then
             isforb=0
             do ifrag=1,input_frag%nfrag
                ! find reference fragment this corresponds to
                ifrag_ref=input_frag%frag_index(ifrag)
                ! loop over orbitals of this fragment
                do iforb=1,ref_frags(ifrag_ref)%fbasis%forbs%norb
                   if (iiorb==iforb+isforb) exit
                end do
                if (iforb/=ref_frags(ifrag_ref)%fbasis%forbs%norb+1) exit
                isforb=isforb+ref_frags(ifrag_ref)%fbasis%forbs%norb
             end do
             write(orbname,*) iforb
             fragdir=trim(input_frag%dirname(ifrag_ref))
          else
             write(orbname,*) iiorb
             fragdir=trim(input_frag%dirname(1))
          end if

          !! first check if file exists
          !inquire(file=trim(input_dir)//trim(fragdir)//'tmbisf'//trim(adjustl(orbname))//'.dat',exist=psirold_ok)
          !if (.not. psirold_ok) print*,"psirold doesn't exist for reformatting",iiorb,&
          !     trim(input_dir)//'tmbisf'//trim(adjustl(orbname))//'.dat'

          !! read in psirold
          !if (psirold_ok) then
          !   call timing(iproc,'readisffiles','ON')
          !   open(99,file=trim(input_dir)//trim(fragdir)//'tmbisf'//trim(adjustl(orbname))//'.dat',&
          !        form="unformatted",status='unknown')
          !   read(99) dummy
          !   read(99) lzd_old%llr(ilr_old)%d%n1i,lzd_old%llr(ilr_old)%d%n2i,lzd_old%llr(ilr_old)%d%n3i
          !   read(99) lzd_old%llr(ilr_old)%nsi1,lzd_old%llr(ilr_old)%nsi2,lzd_old%llr(ilr_old)%nsi3
          !   psirold=f_malloc((/lzd_old%llr(ilr_old)%d%n1i,lzd_old%llr(ilr_old)%d%n2i,lzd_old%llr(ilr_old)%d%n3i/),id='psirold')
          !   do k=1,lzd_old%llr(ilr_old)%d%n3i
          !      do j=1,lzd_old%llr(ilr_old)%d%n2i
          !         do i=1,lzd_old%llr(ilr_old)%d%n1i
          !            read(99) psirold(i,j,k)
          !         end do
          !      end do
          !   end do
          !   close(99)
          !   call timing(iproc,'readisffiles','OF')
          !end if


          ! Periodicity in the three directions
          gperx=(tmb%lzd%glr%geocode /= 'F')
          gpery=(tmb%lzd%glr%geocode == 'P')
          gperz=(tmb%lzd%glr%geocode /= 'F')

          ! Set the conditions for ext_buffers (conditions for buffer size)
          lperx=(lzd_old%llr(ilr)%geocode /= 'F')
          lpery=(lzd_old%llr(ilr)%geocode == 'P')
          lperz=(lzd_old%llr(ilr)%geocode /= 'F')

          !calculate the size of the buffers of interpolating function grid
          call ext_buffers(gperx,gnbl1,gnbr1)
          call ext_buffers(gpery,gnbl2,gnbr2)
          call ext_buffers(gperz,gnbl3,gnbr3)
          call ext_buffers(lperx,lnbl1,lnbr1)
          call ext_buffers(lpery,lnbl2,lnbr2)
          call ext_buffers(lperz,lnbl3,lnbr3)


          lzd_old%llr(ilr_old)%nsi1=2*lzd_old%llr(ilr_old)%ns1 - (Lnbl1 - Gnbl1)
          lzd_old%llr(ilr_old)%nsi2=2*lzd_old%llr(ilr_old)%ns2 - (Lnbl2 - Gnbl2)
          lzd_old%llr(ilr_old)%nsi3=2*lzd_old%llr(ilr_old)%ns3 - (Lnbl3 - Gnbl3)

          !lzd_old%llr(ilr_old)%d%n1i=2*n_old(1)+31
          !lzd_old%llr(ilr_old)%d%n2i=2*n_old(2)+31
          !lzd_old%llr(ilr_old)%d%n3i=2*n_old(3)+31
          !dimensions of the interpolating scaling functions grid (reduce to +2 for periodic)
          if(lzd_old%llr(ilr)%geocode == 'F') then
             lzd_old%llr(ilr)%d%n1i=2*n_old(1)+31
             lzd_old%llr(ilr)%d%n2i=2*n_old(2)+31
             lzd_old%llr(ilr)%d%n3i=2*n_old(3)+31
          else if(lzd_old%llr(ilr)%geocode == 'S') then
             lzd_old%llr(ilr)%d%n1i=2*n_old(1)+2
             lzd_old%llr(ilr)%d%n2i=2*n_old(2)+31
             lzd_old%llr(ilr)%d%n3i=2*n_old(3)+2
          else
             lzd_old%llr(ilr)%d%n1i=2*n_old(1)+2
             lzd_old%llr(ilr)%d%n2i=2*n_old(2)+2
             lzd_old%llr(ilr)%d%n3i=2*n_old(3)+2
          end if


          psirold_ok=.true.
          workarraytmp=f_malloc((2*n_old+31),id='workarraytmp')
          psirold=f_malloc0((2*n_old+31),id='psirold')

          !call f_zero((2*n_old(1)+31)*(2*n_old(2)+31)*(2*n_old(3)+31),psirold(1,1,1))
          call vcopy((2*n_old(1)+2)*(2*n_old(2)+2)*(2*n_old(3)+2),phigold(0,1,0,1,0,1),1,psirold(1,1,1),1)
          call psig_to_psir_free(n_old(1),n_old(2),n_old(3),workarraytmp,psirold)
          call f_free(workarraytmp)

          !write(*,*) 'iproc,norm psirold ',iproc,dnrm2(product(2*n_old+31),psirold,1),2*n_old+31

          call timing(iproc,'Reformatting ','ON')
          if (psirold_ok) then
             !print*,'using psirold to reformat',iiorb
             ! recalculate centres
             !conditions for periodicity in the three directions
             per(1)=(at%astruct%geocode /= 'F')
             per(2)=(at%astruct%geocode == 'P')
             per(3)=(at%astruct%geocode /= 'F')

             !buffers related to periodicity
             !WARNING: the boundary conditions are not assumed to change between new and old
             call ext_buffers(per(1),nl(1),nr(1))
             call ext_buffers(per(2),nl(2),nr(2))
             call ext_buffers(per(3),nl(3),nr(3))

             ! centre of rotation with respect to start of box
             centre_old_box(1)=frag_trans(iorb)%rot_center(1)-0.5d0*lzd_old%hgrids(1)*(lzd_old%llr(ilr_old)%nsi1-nl(1))
             centre_old_box(2)=frag_trans(iorb)%rot_center(2)-0.5d0*lzd_old%hgrids(2)*(lzd_old%llr(ilr_old)%nsi2-nl(2))
             centre_old_box(3)=frag_trans(iorb)%rot_center(3)-0.5d0*lzd_old%hgrids(3)*(lzd_old%llr(ilr_old)%nsi3-nl(3))

             centre_new_box(1)=frag_trans(iorb)%rot_center_new(1)-0.5d0*tmb%lzd%hgrids(1)*(tmb%lzd%llr(ilr)%nsi1-nl(1))
             centre_new_box(2)=frag_trans(iorb)%rot_center_new(2)-0.5d0*tmb%lzd%hgrids(2)*(tmb%lzd%llr(ilr)%nsi2-nl(2))
             centre_new_box(3)=frag_trans(iorb)%rot_center_new(3)-0.5d0*tmb%lzd%hgrids(3)*(tmb%lzd%llr(ilr)%nsi3-nl(3))

             da=centre_new_box-centre_old_box-(lzd_old%hgrids-tmb%lzd%hgrids)*0.5d0

             call reformat_one_supportfunction(tmb%lzd%llr(ilr),lzd_old%llr(ilr_old),at%astruct%geocode,& !,tmb%lzd%llr(ilr)%geocode,&
                  lzd_old%hgrids,n_old,phigold,tmb%lzd%hgrids,n,centre_old_box,centre_new_box,da,&
                  frag_trans(iorb),tmb%psi(jstart:),psirold)
             call f_free(psirold)
          else ! don't have psirold from file, so reformat using old way
             call reformat_one_supportfunction(tmb%lzd%llr(ilr),lzd_old%llr(ilr_old),at%astruct%geocode,& !,tmb%lzd%llr(ilr)%geocode,&
                  lzd_old%hgrids,n_old,phigold,tmb%lzd%hgrids,n,centre_old_box,centre_new_box,da,&
                  frag_trans(iorb),tmb%psi(jstart:))
          end if
          call timing(iproc,'Reformatting ','OF')
          jstart=jstart+tmb%lzd%llr(ilr)%wfd%nvctr_c+7*tmb%lzd%llr(ilr)%wfd%nvctr_f

          if (present(phi_array_old)) then
             nullify(phigold)
          else
             call f_free_ptr(phigold)
          end if

      end if

   end do

  ! Get the maximal shift among all tasks
  if (nproc>1) then
      call mpiallred(max_shift, 1, mpi_max, comm=bigdft_mpi%mpi_comm)
  end if
  if (iproc==0) call yaml_map('max shift of a locreg center',max_shift,fmt='(es9.2)')

  ! Determine the dumping factor for the confinement. In the limit where the atoms
  ! have not moved, it goes to zero; in the limit where they have moved a lot, it goes to one.
  tt = exp(max_shift*3.465735903d0) - 1.d0 !exponential which is 0 at 0.0 and 1 at 0.2
  tt = min(tt,1.d0) !make sure that the value is not larger than 1.0
  tmb%damping_factor_confinement = tt


  if (add_derivatives) then
     call f_free(phi_old_der)
  end if

  call print_reformat_summary(iproc,nproc,reformat_reason)

END SUBROUTINE reformat_supportfunctions


!> Checks whether reformatting is needed based on various criteria and returns final shift and centres needed for reformat
subroutine reformat_check(reformat_needed,reformat_reason,tol,at,hgrids_old,hgrids,nvctr_c_old,nvctr_f_old,&
       nvctr_c,nvctr_f,n_old,n,ns_old,ns,nglr_old,nglr,geocode,frag_trans,centre_old_box,centre_new_box,da,wrap_around)
  use module_base
  use module_types
  use module_fragments
  use yaml_output
  use box
  use bounds, only: ext_buffers_coarse
  implicit none

  logical, intent(out) :: reformat_needed ! logical telling whether reformat is needed
  integer, dimension(0:7), intent(inout) :: reformat_reason ! array giving reasons for reformatting
  real(gp), intent(in) :: tol ! tolerance for rotations and shifts
  type(atoms_data), intent(in) :: at
  real(gp), dimension(3), intent(in) :: hgrids, hgrids_old
  integer, intent(in) :: nvctr_c, nvctr_f, nvctr_c_old, nvctr_f_old
  integer, dimension(3), intent(in) :: n, n_old, ns, ns_old, nglr_old, nglr
  character(len=1), intent(in) :: geocode !< @copydoc poisson_solver::doc::geocode
  real(gp), dimension(3), intent(out) :: centre_old_box, centre_new_box ! centres of rotation wrt box
  real(gp), dimension(3), intent(out) :: da ! shift to be used in reformat
  type(fragment_transformation), intent(in) :: frag_trans ! includes centres of rotation in global coordinates, shift and angle
  logical, intent(out) :: wrap_around ! periodic case where no reformatting is needed but tmb needs 'unwrapping'


  ! local variables
  real(gp) :: displ !, mindist
  integer, dimension(3) :: nb
  logical, dimension(3) :: per
  integer :: i
  !real(gp), dimension(3) :: centre_new
  type(cell) :: mesh


  mesh=cell_new(at%astruct%geocode,n,hgrids)
  !conditions for periodicity in the three directions
  per(1)=(geocode /= 'F')
  per(2)=(geocode == 'P')
  per(3)=(geocode /= 'F')

  !buffers related to periodicity
  !WARNING: the boundary conditions are not assumed to change between new and old
  call ext_buffers_coarse(per(1),nb(1))
  call ext_buffers_coarse(per(2),nb(2))
  call ext_buffers_coarse(per(3),nb(3))

  !use new (internal) version of mindist, mindist doesn't do the right thing in this case
  ! centre of rotation with respect to start of box
  centre_old_box(1)=mindist_new(per(1),at%astruct%cell_dim(1),frag_trans%rot_center(1),hgrids_old(1)*(ns_old(1)-0.5_dp*nb(1)))
  centre_old_box(2)=mindist_new(per(2),at%astruct%cell_dim(2),frag_trans%rot_center(2),hgrids_old(2)*(ns_old(2)-0.5_dp*nb(2)))
  centre_old_box(3)=mindist_new(per(3),at%astruct%cell_dim(3),frag_trans%rot_center(3),hgrids_old(3)*(ns_old(3)-0.5_dp*nb(3)))
  centre_new_box(1)=mindist_new(per(1),at%astruct%cell_dim(1),frag_trans%rot_center_new(1),hgrids(1)*(ns(1)-0.5_dp*nb(1)))
  centre_new_box(2)=mindist_new(per(2),at%astruct%cell_dim(2),frag_trans%rot_center_new(2),hgrids(2)*(ns(2)-0.5_dp*nb(2)))
  centre_new_box(3)=mindist_new(per(3),at%astruct%cell_dim(3),frag_trans%rot_center_new(3),hgrids(3)*(ns(3)-0.5_dp*nb(3)))
  !centre_old_box=closest_r(mesh,frag_trans%rot_center,center=hgrids_old*(ns_old-0.5_dp*nb))
  !centre_new_box=closest_r(mesh,frag_trans%rot_center_new,center=hgrids*(ns-0.5_dp*nb))

  ! centre_new_box(1)=mindist(per(1),at%astruct%cell_dim(1),frag_trans%rot_center_new(1),hgrids(1)*(ns(1)-0.5_dp*nb(1)))
  ! centre_new_box(2)=mindist(per(2),at%astruct%cell_dim(2),frag_trans%rot_center_new(2),hgrids(2)*(ns(2)-0.5_dp*nb(2)))
  ! centre_new_box(3)=mindist(per(3),at%astruct%cell_dim(3),frag_trans%rot_center_new(3),hgrids(3)*(ns(3)-0.5_dp*nb(3)))

  !print*,'rotated nb',trim(yaml_toa(rotate_vector(frag_trans%rot_axis,frag_trans%theta,hgrids*-0.5_dp*nb),fmt='(f12.8)'))
  !print*,'rotated centre old',trim(yaml_toa(rotate_vector(frag_trans%rot_axis,frag_trans%theta,centre_old_box),fmt='(f12.8)'))
  !print*,'rotated centre new',trim(yaml_toa(rotate_vector(frag_trans%rot_axis,-frag_trans%theta,centre_new_box),fmt='(f12.8)'))
  !Calculate the shift of the atom to be used in reformat
  !da(1)=mindist(per(1),at%astruct%cell_dim(1),centre_new_box(1),centre_old_box(1))
  !da(2)=mindist(per(2),at%astruct%cell_dim(2),centre_new_box(2),centre_old_box(2))
  !da(3)=mindist(per(3),at%astruct%cell_dim(3),centre_new_box(3),centre_old_box(3))
  da=centre_new_box-centre_old_box-(hgrids_old-hgrids)*0.5d0

  !print*,'reformat check',frag_trans%rot_center(2),ns_old(2),centre_old_box(2),&
  !     frag_trans%rot_center_new(2),ns(2),centre_new_box(2),da(2)
  !write(*,'(a,15I4)')'nb,ns_old,ns,n_old,n',nb,ns_old,ns,n_old,n
  !write(*,'(a,3(3(f12.8,x),3x))') 'final centre box',centre_old_box,centre_new_box,da
  !write(*,'(a,3(3(f12.8,x),3x))') 'final centre',frag_trans%rot_center,frag_trans%rot_center_new

  displ=square(mesh,da)!sqrt(da(1)**2+da(2)**2+da(3)**2)

  !reformatting criterion
  if (hgrids(1) == hgrids_old(1) .and. hgrids(2) == hgrids_old(2) .and. hgrids(3) == hgrids_old(3) &
        .and. nvctr_c  == nvctr_c_old .and. nvctr_f  == nvctr_f_old &
        .and. n_old(1)==n(1)  .and. n_old(2)==n(2) .and. n_old(3)==n(3) &
        .and. abs(frag_trans%theta) <= tol .and. abs(displ) <= tol) then
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
  end if

  ! check to make sure we don't need to 'unwrap' (or wrap) tmb in periodic case
  wrap_around=.false.
  !if (.not. reformat_needed) then
     do i=1,3
        if (tmb_wrap(per(i),ns(i),n(i),nglr(i)) .or. tmb_wrap(per(i),ns_old(i),n_old(i),nglr_old(i))) then
           if (ns(i) /= ns_old(i)) then
              wrap_around = .true.
              exit
           end if
        end if
     end do

     if (wrap_around) then
         reformat_reason(7) = reformat_reason(7) + 1
     end if
  !end if

  !write(*,'(a,1(1x,I4),3(2x,F12.6),6(2x,I4),9(2x,F12.6))')'DEBUG:iproc,rc,ns,nb,alat,h,cob',&
  !     bigdft_mpi%iproc, frag_trans%rot_center, ns_old, nb, at%astruct%cell_dim, hgrids_old, centre_old_box, wrap_around

contains

  !> Checks to see if a tmb wraps around the the unit cell
  !! knowing that there could have been a modulo operation
  function tmb_wrap(periodic,ns,n,nglr)
    use module_base
    implicit none
    logical, intent(in) :: periodic
    integer, intent(in) :: ns, n, nglr
    logical :: tmb_wrap

    tmb_wrap = .false.
    if (periodic) then
       !<=?
       if (ns<0 .or. ns+n>nglr) then
          tmb_wrap = .true.
       end if
    end if

  end function tmb_wrap


  !> Calculates the minimum difference between two coordinates
  !! knowing that there could have been a modulo operation
  function mindist_new(periodic,alat,r,r_old)
    use module_base
    implicit none
    logical, intent(in) :: periodic
    real(gp), intent(in) :: r,r_old,alat
    real(gp) :: mindist_new

    mindist_new = r - r_old
    if (periodic) then
       if (mindist_new > 0.5d0*alat) then
          mindist_new = mindist_new - alat
       else if (mindist_new < -0.5d0*alat) then
          mindist_new = mindist_new + alat
       end if
    end if

  end function mindist_new


end subroutine reformat_check


!> Print information about the reformatting due to restart
subroutine print_reformat_summary(iproc,nproc,reformat_reason)
  use module_base
  use module_types
  use yaml_output
  implicit none

  integer, intent(in) :: iproc,nproc
  integer, dimension(0:7), intent(inout) :: reformat_reason ! array giving reasons for reformatting

  if (nproc > 1) call mpiallred(reformat_reason, mpi_sum, comm=bigdft_mpi%mpi_comm)

  if (iproc==0) then
        call yaml_mapping_open('Overview of the reformatting (several categories may apply)')
        call yaml_map('No reformatting required', reformat_reason(0))
        call yaml_map('Grid spacing has changed', reformat_reason(1))
        call yaml_map('Number of coarse grid points has changed', reformat_reason(2))
        call yaml_map('Number of fine grid points has changed', reformat_reason(3))
        call yaml_map('Box size has changed', reformat_reason(4))
        call yaml_map('Molecule was shifted', reformat_reason(5))
        call yaml_map('Molecule was rotated', reformat_reason(6))
        call yaml_map('Wrapping/unwrapping', reformat_reason(7))
        call yaml_mapping_close()
  end if

end subroutine print_reformat_summary


! re-written to be closer to compress_plain and version in writeonewave_linear
subroutine psi_to_psig(n,nseg_c,nvctr_c,keyg_c,keyv_c,nseg_f,nvctr_f,keyg_f,keyv_f,psig,psi_c,psi_f)
  use module_base
  implicit none

  integer, dimension(3), intent(in) :: n
  integer, intent(in) :: nseg_c, nseg_f, nvctr_c, nvctr_f
  integer, dimension(nseg_c), intent(in) :: keyv_c
  integer, dimension(nseg_f), intent(in) :: keyv_f
  integer, dimension(2,nseg_c), intent(in) :: keyg_c
  integer, dimension(2,nseg_f), intent(in) :: keyg_f
  real(wp), dimension(0:n(1),2,0:n(2),2,0:n(3),2), intent(out) :: psig
  real(wp), dimension(nvctr_c), intent(in) :: psi_c
  real(wp), dimension(7,nvctr_f), intent(in) :: psi_f

  ! local variables
  integer :: iseg, jj, j0, j1, ii, i1, i2, i3, i0, i, n1p1, np

  call f_zero(psig)

  n1p1=n(1)+1
  np=n1p1*(n(2)+1)

  !$omp parallel default(private) shared(keyv_c,keyv_f,keyg_c,keyg_f,psig,psi_c,psi_f) &
  !$omp shared(n1p1,np,nseg_c,nseg_f)
  ! coarse part
  !$omp do
  do iseg=1,nseg_c
     jj=keyv_c(iseg)
     j0=keyg_c(1,iseg)
     j1=keyg_c(2,iseg)
     ii=j0-1
     i3=ii/np
     ii=ii-i3*np
     i2=ii/n1p1
     i0=ii-i2*n1p1
     i1=i0+j1-j0
     do i=i0,i1
        psig(i,1,i2,1,i3,1) = psi_c(i-i0+jj)
     enddo
  enddo
  !$omp enddo
  ! fine part
  !$omp do
  do iseg=1,nseg_f
     jj=keyv_f(iseg)
     j0=keyg_f(1,iseg)
     j1=keyg_f(2,iseg)
     ii=j0-1
     i3=ii/np
     ii=ii-i3*np
     i2=ii/n1p1
     i0=ii-i2*n1p1
     i1=i0+j1-j0
     do i=i0,i1
        psig(i,2,i2,1,i3,1) = psi_f(1,i-i0+jj)
        psig(i,1,i2,2,i3,1) = psi_f(2,i-i0+jj)
        psig(i,2,i2,2,i3,1) = psi_f(3,i-i0+jj)
        psig(i,1,i2,1,i3,2) = psi_f(4,i-i0+jj)
        psig(i,2,i2,1,i3,2) = psi_f(5,i-i0+jj)
        psig(i,1,i2,2,i3,2) = psi_f(6,i-i0+jj)
        psig(i,2,i2,2,i3,2) = psi_f(7,i-i0+jj)
     enddo
  enddo
  !$omp enddo
  !$omp end parallel


end subroutine psi_to_psig

!> Associate to the absolute value of orbital a filename which depends of the k-point and
!! of the spin sign
subroutine filename_of_proj(lbin,filename,ikpt,iat,iproj,icplx,filename_out)
  use module_base
  use module_types
  implicit none
  character(len=*), intent(in) :: filename
  logical, intent(in) :: lbin
  integer, intent(in) :: ikpt,iat,iproj,icplx
  character(len=*), intent(out) :: filename_out
  !local variables
  character(len=1) :: realimag
  character(len=3) :: f2
  character(len=4) :: f3
  character(len=5) :: f4
  character(len=13) :: completename

  !calculate k-point
  write(f3,'(a1,i3.3)') "k", ikpt !not more than 999 kpts

  !see if the wavefunction is real or imaginary
  if(icplx==2) then
     realimag='I'
  else
     realimag='R'
  end if

  !value of the atom
  write(f4,'(a1,i4.4)') "a", iat

  !value of the atom
  write(f2,'(i3.3)') iproj

  !complete the information in the name of the orbital
  completename='-'//f3//'-'//f4//'-'//realimag
  if (lbin) then
     filename_out = trim(filename)//completename//".bin."//f2
     !print *,'complete name <',trim(filename_out),'> end'
 else
     filename_out = trim(filename)//completename//"."//f2
     !print *,'complete name <',trim(filename_out),'> end'
 end if

  !print *,'filename: ',filename_out
end subroutine filename_of_proj

!> Write all projectors
subroutine writemyproj(filename,iformat,orbs,hx,hy,hz,at,rxyz,nl)
  use module_types
  use module_base
  use yaml_output
  use gaussians
  use public_enums, only: WF_FORMAT_ETSF, WF_FORMAT_BINARY
  implicit none
  integer, intent(in) :: iformat
  real(gp), intent(in) :: hx,hy,hz
  type(atoms_data), intent(in) :: at
  type(orbitals_data), intent(in) :: orbs
  type(DFT_PSP_projectors), intent(in) :: nl
  real(gp), dimension(3,at%astruct%nat), intent(in) :: rxyz
  character(len=*), intent(in) :: filename
  !Local variables
  type(gaussian_basis_iter) :: iter
  integer :: ncount1,ncount2,ncount_rate,ncount_max
  integer :: iat,ikpt,iproj,iskpt,iekpt,istart,ncplx_k,icplx,l
  integer :: mbseg_c,mbseg_f,mbvctr_c,mbvctr_f
  real(kind=4) :: tr0,tr1
  real(kind=8) :: tel
  character(len=500) :: filename_out
  logical :: lbin

  call yaml_map('Write projectors to file', trim(filename) // '.*')
  !if (iproc == 0) write(*,"(1x,A,A,a)") "Write wavefunctions to file: ", trim(filename),'.*'
  if (iformat == WF_FORMAT_ETSF) then
     !call write_waves_etsf(iproc,filename,orbs,n1,n2,n3,hx,hy,hz,at,rxyz,wfd,psi)
     stop "not implemented proj in ETSF"
  else
     call cpu_time(tr0)
     call system_clock(ncount1,ncount_rate,ncount_max)

     !create projectors for any of the k-point hosted by the processor
     !starting kpoint
     if (orbs%norbp > 0) then
        iskpt=orbs%iokpt(1)
        iekpt=orbs%iokpt(orbs%norbp)
     else
        iskpt=1
        iekpt=1
     end if
     lbin = (iformat == WF_FORMAT_BINARY)

     do ikpt=iskpt,iekpt
        ncplx_k = 2
        if (orbs%kpts(1,ikpt) == 0 .and. orbs%kpts(2,ikpt) == 0 .and. &
             & orbs%kpts(3,ikpt) == 0) ncplx_k = 1
        do iat=1,at%astruct%nat

           call plr_segs_and_vctrs(nl%pspd(iat)%plr,mbseg_c,mbseg_f,mbvctr_c,mbvctr_f)
           ! Start a gaussian iterator.
           call gaussian_iter_start(nl%proj_G, iat, iter)
           iproj = 0
           istart = 1
           do
              if (.not. gaussian_iter_next_shell(nl%proj_G, iter)) exit
              do l = 1, 2*iter%l-1
                 iproj = iproj + 1
                 do icplx = 1, ncplx_k
                    call filename_of_proj(lbin,filename,&
                         & ikpt,iat,iproj,icplx,filename_out)
                    if (lbin) then
                       open(unit=99,file=trim(filename_out),&
                            & status='unknown',form="unformatted")
                    else
                       open(unit=99,file=trim(filename_out),status='unknown')
                    end if
                    call writeonewave(99,.not.lbin,iproj,&
                         & nl%pspd(iat)%plr%d%n1, &
                         & nl%pspd(iat)%plr%d%n2, &
                         & nl%pspd(iat)%plr%d%n3, &
                         & hx,hy,hz, at%astruct%nat,rxyz, &
                         & mbseg_c, mbvctr_c, &
                         & nl%pspd(iat)%plr%wfd%keyglob(1,1), &
                         & nl%pspd(iat)%plr%wfd%keyvglob(1), &
                         & mbseg_f, mbvctr_f, &
                         & nl%pspd(iat)%plr%wfd%keyglob(1,mbseg_c+1), &
                         & nl%pspd(iat)%plr%wfd%keyvglob(mbseg_c+1), &
                         & nl%proj(istart), nl%proj(istart + mbvctr_c), &
                         & UNINITIALIZED(1._wp))

                    close(99)
                    istart = istart + (mbvctr_c+7*mbvctr_f)
                 end do
              end do
           end do
        end do
     enddo

     call cpu_time(tr1)
     call system_clock(ncount2,ncount_rate,ncount_max)
     tel=dble(ncount2-ncount1)/dble(ncount_rate)
     call yaml_sequence_open('Write Proj Time')
     call yaml_sequence(advance='no')
     call yaml_mapping_open(flow=.true.)
     call yaml_map('Timing',(/ real(tr1-tr0,kind=8),tel /),fmt='(1pe10.3)')
     call yaml_mapping_close()
     call yaml_sequence_close()
     !write(*,'(a,i4,2(1x,1pe10.3))') '- WRITE WAVES TIME',iproc,tr1-tr0,tel
     !write(*,'(a,1x,i0,a)') '- iproc',iproc,' finished writing waves'
  end if

END SUBROUTINE writemyproj
