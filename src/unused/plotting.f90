!> @file
!!   These routines are not used anymore.
!! @author
!!    Copyright (C) 2007-2013 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS


!> Routine to plot wavefunctions (old version)
subroutine plot_wf_old(kindplot,orbname,nexpo,at,lr,hx,hy,hz,rxyz,psi,comment)
  use module_base
  use module_types
  implicit none
  character(len=*) :: kindplot
  character(len=10) :: comment
  character(len=11) :: orbname
  integer, intent(in) :: nexpo
  real(gp), intent(in) :: hx,hy,hz
  type(atoms_data), intent(in) :: at
  real(gp), dimension(3,at%nat), intent(in) :: rxyz
  type(locreg_descriptors), intent(in) :: lr
  real(wp), dimension(*) :: psi!wfd%nvctr_c+7*wfd%nvctr_f
  !local variables
  character(len=*), parameter :: subname='plot_wf_old'
  integer :: i_stat,i_all
  integer :: nl1,nl2,nl3,n1i,n2i,n3i,n1,n2,n3
  type(workarr_sumrho) :: w
  real(wp), dimension(:), allocatable :: psir

  n1=lr%d%n1
  n2=lr%d%n2
  n3=lr%d%n3
  n1i=lr%d%n1i
  n2i=lr%d%n2i
  n3i=lr%d%n3i

  call initialize_work_arrays_sumrho(lr,w)
 
  allocate(psir(lr%d%n1i*lr%d%n2i*lr%d%n3i+ndebug),stat=i_stat)
  call memocc(i_stat,psir,'psir',subname)
  !initialisation
  if (lr%geocode == 'F') then
     call to_zero(lr%d%n1i*lr%d%n2i*lr%d%n3i,psir)
  end if
 
  call daub_to_isf(lr,w,psi,psir)

  if (lr%geocode /= 'F') then
     nl1=1
     nl3=1
  else
     nl1=14
     nl3=14
  end if
  !value of the buffer in the y direction
  if (lr%geocode == 'P') then
     nl2=1
  else
     nl2=14
  end if

  if (trim(kindplot)=='POT') then
     !call plot_pot(rx,ry,rz,hx,hy,hz,n1,n2,n3,n1i,n2i,n3i,nl1,nl2,nl3,iounit,psir)
     call plot_pot_full(nexpo,hx,hy,hz,n1,n2,n3,n1i,n2i,n3i,&
          nl1,nl2,nl3,orbname,psir,comment)
  else if (trim(kindplot)=='CUBE') then
     call plot_cube_full(nexpo,at,rxyz,hx,hy,hz,n1,n2,n3,n1i,n2i,n3i,&
          nl1,nl2,nl3,orbname,psir,comment)
  else
     stop 'ERROR: plotting format not recognized'
  end if
  i_all=-product(shape(psir))*kind(psir)
  deallocate(psir,stat=i_stat)
  call memocc(i_stat,i_all,'psir',subname)

  call deallocate_work_arrays_sumrho(w)

END SUBROUTINE plot_wf_old
!!$
!!$
subroutine plot_wf_cube(orbname,at,lr,hx,hy,hz,rxyz,psi,comment)
  use module_base
  use module_types
  implicit none
!HU  character(len=10) :: orbname,comment 
  character(len=10) :: comment
  character(len=11) :: orbname
  type(atoms_data), intent(in) :: at
  real(gp), intent(in) :: hx,hy,hz
  real(gp), dimension(3,at%nat), intent(in) :: rxyz
  type(locreg_descriptors), intent(in) :: lr
  real(wp), dimension(*) :: psi!wfd%nvctr_c+7*wfd%nvctr_f
  !local variables
  character(len=*), parameter :: subname='plot_wf'
  integer :: nw1,nw2,i_stat,i_all,i,n1,n2,n3,nfl1,nfl2,nfl3,nfu1,nfu2,nfu3
  integer :: nxc,nxf,nl1,nl2,nl3,n1i,n2i,n3i
  !n(c) real(gp) :: rx, ry, rz
  real(wp), dimension(0:3) :: scal
  real(wp), dimension(:), allocatable :: psir,w1,w2,x_c_psifscf,x_f_psig

  !n(c) rx=rxyz(1,1)
  !n(c) ry=rxyz(2,1)
  !n(c) rz=rxyz(3,1)

  n1=lr%d%n1
  n2=lr%d%n2
  n3=lr%d%n3
  n1i=lr%d%n1i
  n2i=lr%d%n2i
  n3i=lr%d%n3i
  nfl1=lr%d%nfl1
  nfl2=lr%d%nfl2
  nfl3=lr%d%nfl3
  nfu1=lr%d%nfu1
  nfu2=lr%d%nfu2
  nfu3=lr%d%nfu3

  do i=0,3
     scal(i)=1.d0
  enddo

  select case(lr%geocode)
  case('F')
     !dimension of the work arrays
     ! shrink convention: nw1>nw2
     nw1=max((n3+1)*(2*n1+31)*(2*n2+31),& 
          (n1+1)*(2*n2+31)*(2*n3+31),&
          2*(nfu1-nfl1+1)*(2*(nfu2-nfl2)+31)*(2*(nfu3-nfl3)+31),&
          2*(nfu3-nfl3+1)*(2*(nfu1-nfl1)+31)*(2*(nfu2-nfl2)+31))
     nw2=max(4*(nfu2-nfl2+1)*(nfu3-nfl3+1)*(2*(nfu1-nfl1)+31),&
          4*(nfu1-nfl1+1)*(nfu2-nfl2+1)*(2*(nfu3-nfl3)+31),&
          (n1+1)*(n2+1)*(2*n3+31),&
          (2*n1+31)*(n2+1)*(n3+1))
     nxc=(n1+1)*(n2+1)*(n3+1)!(2*n1+2)*(2*n2+2)*(2*n3+2)
     nxf=7*(nfu1-nfl1+1)*(nfu2-nfl2+1)*(nfu3-nfl3+1)
  case('S')
     !dimension of the work arrays
     nw1=1
     nw2=1
     nxc=(2*n1+2)*(2*n2+31)*(2*n3+2)
     nxf=1
  case('P')
     !dimension of the work arrays, fully periodic case
     nw1=1
     nw2=1
     nxc=(2*n1+2)*(2*n2+2)*(2*n3+2)
     nxf=1
  end select
  !work arrays
  allocate(x_c_psifscf(nxc+ndebug),stat=i_stat)
  call memocc(i_stat,x_c_psifscf,'x_c_psifscf',subname)
  allocate(x_f_psig(nxf+ndebug),stat=i_stat)
  call memocc(i_stat,x_f_psig,'x_f_psig',subname)
  allocate(w1(nw1+ndebug),stat=i_stat)
  call memocc(i_stat,w1,'w1',subname)
  allocate(w2(nw2+ndebug),stat=i_stat)
  call memocc(i_stat,w2,'w2',subname)

  allocate(psir(n1i*n2i*n3i+ndebug),stat=i_stat)
  call memocc(i_stat,psir,'psir',subname)
  !initialisation
  if (lr%geocode == 'F') then
     call to_zero(nxc,x_c_psifscf)
     call to_zero(nxf,x_f_psig)
     call to_zero(n1i*n2i*n3i,psir)
  end if
  select case(lr%geocode)
  case('F')
     call uncompress_forstandard_short(n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3, & 
          lr%wfd%nseg_c,lr%wfd%nvctr_c,lr%wfd%keyg(1,1),lr%wfd%keyv(1),  & 
          lr%wfd%nseg_f,lr%wfd%nvctr_f,&
          lr%wfd%keyg(1,lr%wfd%nseg_c+1),lr%wfd%keyv(lr%wfd%nseg_c+1), &
          scal,psi(1),psi(lr%wfd%nvctr_c+1),x_c_psifscf,x_f_psig)

     call comb_grow_all(n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,w1,w2,&
          x_c_psifscf,x_f_psig,  & 
          psir(1),lr%bounds%kb%ibyz_c,lr%bounds%gb%ibzxx_c,lr%bounds%gb%ibxxyy_c,&
          lr%bounds%gb%ibyz_ff,lr%bounds%gb%ibzxx_f,lr%bounds%gb%ibxxyy_f)
     
  case('P')
     call uncompress_per(n1,n2,n3,lr%wfd%nseg_c,&
          lr%wfd%nvctr_c,lr%wfd%keyg(1,1),lr%wfd%keyv(1),&
          lr%wfd%nseg_f,lr%wfd%nvctr_f,&
          lr%wfd%keyg(1,lr%wfd%nseg_c+1),lr%wfd%keyv(lr%wfd%nseg_c+1),&
          psi(1),psi(lr%wfd%nvctr_c+1),x_c_psifscf,psir(1))

     call convolut_magic_n_per_self(2*n1+1,2*n2+1,2*n3+1,&
          x_c_psifscf,psir(1))
     
  case('S')
     call uncompress_slab(n1,n2,n3,lr%wfd%nseg_c,lr%wfd%nvctr_c,&
          lr%wfd%keyg(1,1),lr%wfd%keyv(1),&
          lr%wfd%nseg_f,lr%wfd%nvctr_f,lr%wfd%keyg(1,lr%wfd%nseg_c+1),&
          lr%wfd%keyv(lr%wfd%nseg_c+1),   &
          psi(1),psi(lr%wfd%nvctr_c+1),x_c_psifscf,psir(1))

     call convolut_magic_n_slab_self(2*n1+1,2*n2+15,2*n3+1,x_c_psifscf,psir(1))
  end select

  i_all=-product(shape(x_c_psifscf))*kind(x_c_psifscf)
  deallocate(x_c_psifscf,stat=i_stat)
  call memocc(i_stat,i_all,'x_c_psifscf',subname)
  i_all=-product(shape(x_f_psig))*kind(x_f_psig)
  deallocate(x_f_psig,stat=i_stat)
  call memocc(i_stat,i_all,'x_f_psig',subname)
  i_all=-product(shape(w1))*kind(w1)
  deallocate(w1,stat=i_stat)
  call memocc(i_stat,i_all,'w1',subname)
  i_all=-product(shape(w2))*kind(w2)
  deallocate(w2,stat=i_stat)
  call memocc(i_stat,i_all,'w2',subname)

  if (lr%geocode /= 'F') then
     nl1=1
     nl3=1
  else
     nl1=14
     nl3=14
  end if
  !value of the buffer in the y direction
  if (lr%geocode == 'P') then
     nl2=1
  else
     nl2=14
  end if


  !call plot_pot(rx,ry,rz,hx,hy,hz,n1,n2,n3,n1i,n2i,n3i,nl1,nl2,nl3,iounit,psir)
!HU  call plot_pot_full(rx,ry,rz,hx,hy,hz,n1,n2,n3,n1i,n2i,n3i,&
!HU       nl1,nl2,nl3,orbname,psir,comment)
  call plot_cube_full(1,at,rxyz,hx,hy,hz,n1,n2,n3,n1i,n2i,n3i,&
       nl1,nl2,nl3,orbname,psir,comment)

  i_all=-product(shape(psir))*kind(psir)
  deallocate(psir,stat=i_stat)
  call memocc(i_stat,i_all,'psir',subname)

END SUBROUTINE plot_wf_cube
!!$
!!$
subroutine plot_pot(rx,ry,rz,hx,hy,hz,n1,n2,n3,n1i,n2i,n3i,nl1,nl2,nl3,iounit,pot)
  use module_base
  implicit none
  integer, intent(in) :: n1,n2,n3,iounit,n1i,n2i,n3i,nl1,nl2,nl3
  real(gp), intent(in) :: rx,ry,rz,hx,hy,hz
  real(dp), dimension(*), intent(in) :: pot
  !local variables
  integer :: i1,i2,i3,ind
  real(gp) :: hxh,hyh,hzh

  hxh=.5_gp*hx
  hyh=.5_gp*hy
  hzh=.5_gp*hz

  open(iounit) 
  open(iounit+1) 
  open(iounit+2) 

  i3=nint(rz/hzh)
  i2=nint(ry/hyh)
  write(*,*) 'plot_p, i2,i3,n2,n3 ',i2,i3,n2,n3
  do i1=0,2*n1
     ind=i1+nl1+(i2+nl2-1)*n1i+(i3+nl3-1)*n1i*n2i
     write(iounit,*) real(i1,gp)*hxh,pot(ind)
  enddo

  i1=nint(rx/hxh)
  i2=nint(ry/hyh)
  write(*,*) 'plot_p, i1,i2 ',i1,i2
  do i3=0,2*n3
     ind=i1+nl1+(i2+nl2-1)*n1i+(i3+nl3-1)*n1i*n2i
     write(iounit+1,*) real(i3,gp)*hzh,pot(ind)
  enddo

  i1=nint(rx/hxh)
  i3=nint(rz/hzh)
  write(*,*) 'plot_p, i1,i3 ',i1,i3
  do i2=0,2*n2
     ind=i1+nl1+(i2+nl2-1)*n1i+(i3+nl3-1)*n1i*n2i
     write(iounit+2,*) real(i2,gp)*hyh,pot(ind)
  enddo

  close(iounit) 
  close(iounit+1) 
  close(iounit+2) 

END SUBROUTINE plot_pot


subroutine plot_pot_full(nexpo,hx,hy,hz,n1,n2,n3,n1i,n2i,n3i,&
     nl1,nl2,nl3,orbname,pot,comment)
  use module_base
  implicit none
  character(len=10), intent(in) :: comment
  character(len=11), intent(in) :: orbname
  integer, intent(in) :: n1,n2,n3,nl1,nl2,nl3,n1i,n2i,n3i,nexpo
  real(gp), intent(in) :: hx,hy,hz
  real(dp), dimension(*), intent(in) :: pot
  !local variables
  integer :: i1,i2,i3,ind
  real(gp) :: hxh,hyh,hzh

!virtual orbitals are identified by their name
!  write(*,'(1x,a,i0)')'printing orbital number= ',iounit
  write(*,'(A)')'printing '//orbname
  hxh=.5_gp*hx
  hyh=.5_gp*hy
  hzh=.5_gp*hz
!  write(orbname,'(i0)')iounit
!  open(unit=22,file='psi'//orbname//'.pot',status='unknown')
  open(unit=22,file=orbname//'.pot',status='unknown')
!  write(22,*)'orbital'//orbname
  write(22,*)orbname
  write(22,*) 2*n1+2,2*n2+2,2*n3+2
  !here there can be shifts in visualising periodic objects (n1 -> n1+1)
  write(22,*) n1*hx,' 0. ',n2*hy
  write(22,*) ' 0. ',' 0. ',n3*hz
  write(22,*)'xyz periodic '//comment
  do i3=0,2*n3+1
     do i2=0,2*n2+1
        do i1=0,2*n1+1
           ind=i1+nl1+(i2+nl2-1)*n1i+(i3+nl3-1)*n1i*n2i
           write(22,*)pot(ind)**nexpo
        end do
     end do
  end do
  close(unit=22) 
!  close(unit=23) 

END SUBROUTINE plot_pot_full


subroutine plot_cube_full(nexpo,at,rxyz,hx,hy,hz,n1,n2,n3,n1i,n2i,n3i,&
     nl1,nl2,nl3,orbname,pot,comment)
  use module_base
  use module_types
  implicit none
!HU  character(len=10), intent(in) :: orbname,comment
  character(len=10), intent(in) :: comment
  character(len=11), intent(in) :: orbname
  type(atoms_data), intent(in) :: at
  real(gp), dimension(3,at%nat), intent(in) :: rxyz
  integer, intent(in) :: n1,n2,n3,nl1,nl2,nl3,n1i,n2i,n3i,nexpo
  real(gp), intent(in) :: hx,hy,hz
  real(dp), dimension(*), intent(in) :: pot
  !local variables
  integer :: i1,i2,i3,ind,iat,j,icount
  real(gp) :: hxh,hyh,hzh
  character(len=3) :: advancestring

!virtual orbitals are identified by their name
!  write(*,'(1x,a,i0)')'printing orbital number= ',iounit
  write(*,'(A)')'printing '//orbname
  hxh=.5_gp*hx
  hyh=.5_gp*hy
  hzh=.5_gp*hz

  open(unit=22,file=trim(orbname)//'.cube',status='unknown')
  write(22,*)orbname
  write(22,*)'CUBE file for orbital wavefunction'
  !number of atoms
!  if (at%geocode=='P') then
     write(22,'(i5,3(f12.6),a)') at%nat,0.0_gp,0.0_gp,0.0_gp,' modified origin'
!  else if (at%geocode=='S') then
!     write(22,'(i5,3(f12.6))') at%nat,0.0_gp,-hyh,0.0_gp
!  else if (at%geocode=='F') then
!     write(22,'(i5,3(f12.6))') at%nat,-hxh,-hyh,-hzh
!  end if
  !grid and grid spacings
  write(22,'(i5,3(f19.12))') 2*n1+2,hxh,0.0_gp,0.0_gp
  write(22,'(i5,3(f19.12))') 2*n2+2,0.0_gp,hyh,0.0_gp
  write(22,'(i5,3(f19.12))') 2*n3+2,0.0_gp,0.0_gp,hzh
  !atomic number and positions
  do iat=1,at%nat
     write(22,'(i5,4(f12.6))') at%nzatom(at%iatype(iat)),at%nelpsp(at%iatype(iat)),(rxyz(j,iat),j=1,3)
  end do

  !the loop is reverted for a cube file
  !charge normalised to the total charge

  do i1=0,2*n1+1
     do i2=0,2*n2+1
        icount=0
        do i3=0,2*n3+1
           icount=icount+1
           if (icount == 6 .or. i3==2*n3+1) then
              advancestring='yes'
              icount=0
           else
              advancestring='no'
           end if

           ind=i1+nl1+(i2+nl2-1)*n1i+(i3+nl3-1)*n1i*n2i
           write(22,'(1x,1pe13.6)',advance=advancestring) pot(ind)**nexpo
        end do
     end do
  end do
  close(22)

END SUBROUTINE plot_cube_full


subroutine plot_psifscf(iunit,hgrid,n1,n2,n3,psifscf)
  use module_base
  implicit none
  integer, intent(in) :: iunit,n1,n2,n3
  real(gp), intent(in) :: hgrid
  real(wp), dimension(-7:2*n1+8,-7:2*n2+8,-7:2*n3+8), intent(in) :: psifscf
  !local variables
  integer :: i1,i2,i3,i
  real(gp) :: hgridh

  hgridh=.5d0*hgrid

  ! along x-axis
  i3=n3
  i2=n2
  do i1=-7,2*n1+8
     write(iunit,'(3(1x,e10.3),1x,e12.5)') &
        real(i1,gp)*hgridh,real(i2,gp)*hgridh,real(i3,gp)*hgridh,psifscf(i1,i2,i3)
  enddo

  ! 111 diagonal
  do i=-7,min(2*n1+8,2*n2+8,2*n3+8)
     i1=i ; i2=i ; i3=i
     write(iunit,'(3(1x,e10.3),1x,e12.5)') &
          real(i1,gp)*hgridh,real(i2,gp)*hgridh,real(i3,gp)*hgridh,psifscf(i1,i2,i3)
  enddo

  ! 1-1-1 diagonal
  do i=-7,min(2*n1+8,2*n2+8,2*n3+8)
     i1=i ; i2=-i ; i3=-i
     write(iunit,'(3(1x,e10.3),1x,e12.5)') &
          real(i1,gp)*hgridh,real(i2,gp)*hgridh,real(i3,gp)*hgridh,psifscf(i1,i2,i3)
  enddo

  ! -11-1 diagonal
  do i=-7,min(2*n1+8,2*n2+8,2*n3+8)
     i1=-i ; i2=i ; i3=-i
     write(iunit,'(3(1x,e10.3),1x,e12.5)') &
          real(i1,gp)*hgridh,real(i2,gp)*hgridh,real(i3,gp)*hgridh,psifscf(i1,i2,i3)
  enddo

  ! -1-11 diagonal
  do i=-7,min(2*n1+8,2*n2+8,2*n3+8)
     i1=-i ; i2=-i ; i3=i
     write(iunit,'(3(1x,e10.3),1x,e12.5)') &
          real(i1,gp)*hgridh,real(i2,gp)*hgridh,real(i3,gp)*hgridh,psifscf(i1,i2,i3)
  enddo

END SUBROUTINE plot_psifscf


subroutine read_potfile(geocode,filename,n1,n2,n3,n1i,n2i,n3i,n3d,i3s,rho)
  use module_base
  implicit none
  character(len=1), intent(in) :: geocode !< @copydoc poisson_solver::doc::geocode
  character(len=*), intent(in) :: filename
  integer, intent(in) :: n1i,n2i,n3i,n3d,n1,n2,n3,i3s
  real(dp), dimension(n1i*n2i*n3d), intent(out) :: rho
  !local variables
  integer :: nl1,nl2,nl3,i1,i2,i3,ind
  real(dp) :: value

  open(unit=22,file=filename,status='unknown')
  read(22,*)!'normalised density'
  read(22,*)!2*n1+2,2*n2+2,2*n3+2
  read(22,*)!alat1,' 0. ',alat2
  read(22,*)!' 0. ',' 0. ',alat3
  read(22,*)!xyz   periodic' !not true in general but needed in the case

  !conditions for periodicity in the three directions
  !value of the buffer in the x and z direction
  if (geocode /= 'F') then
     nl1=1
     nl3=1
  else
     nl1=15
     nl3=15
  end if
  !value of the buffer in the y direction
  if (geocode == 'P') then
     nl2=1
  else
     nl2=15
  end if

  call to_zero(max(n1i*n2i*n3d,1),rho)

  do i3=0,2*n3+1
     do i2=0,2*n2+1
        do i1=0,2*n1+1
           ind=i1+nl1+(i2+nl2-1)*n1i+(i3+nl3-i3s)*n1i*n2i
           read(22,*)value
           if (i3+nl3 >= i3s .and. i3+nl3 <= i3s+n3d-1) then
              rho(ind)=value
           end if
        end do
     end do
  end do
  close(22)

END SUBROUTINE read_potfile


subroutine plot_density_old(geocode,filename,iproc,nproc,n1,n2,n3,n1i,n2i,n3i,n3p,&
     alat1,alat2,alat3,ngatherarr,rho)
  use module_base
  implicit none
  character(len=1), intent(in) :: geocode !< @copydoc poisson_solver::doc::geocode
  character(len=*), intent(in) :: filename
  integer, intent(in) :: iproc,n1i,n2i,n3i,n3p,n1,n2,n3,nproc
  real(gp), intent(in) :: alat1,alat2,alat3
  integer, dimension(0:nproc-1,2), intent(in) :: ngatherarr
  real(dp), dimension(n1i*n2i*n3p), target, intent(in) :: rho
  !local variables
  character(len=*), parameter :: subname='plot_density'
  integer :: nl1,nl2,nl3,i_all,i_stat,i1,i2,i3,ind,ierr,nbxz,nby
  real(dp) :: later_avg
  real(dp), dimension(:), pointer :: pot_ion

  !conditions for periodicity in the three directions
  !value of the buffer in the x and z direction
  if (geocode /= 'F') then
     nl1=1
     nl3=1
     nbxz = 1
  else
     nl1=15
     nl3=15
     nbxz = 0
  end if
  !value of the buffer in the y direction
  if (geocode == 'P') then
     nl2=1
     nby = 1
  else
     nl2=15
     nby = 0
  end if

  if (iproc == 0) then
     open(unit=22,file=filename,status='unknown')
     write(22,*)'normalised density'
     write(22,*) 2*(n1+nbxz),2*(n2+nby),2*(n3+nbxz)
     write(22,*) alat1,' 0. ',alat2
     write(22,*) ' 0. ',' 0. ',alat3
     
     write(22,*)'xyz   periodic' !not true in general but needed in the case
  end if

  if (nproc > 1) then
     !allocate full density in pot_ion array
     allocate(pot_ion(n1i*n2i*n3i+ndebug),stat=i_stat)
     call memocc(i_stat,pot_ion,'pot_ion',subname)

     call MPI_ALLGATHERV(rho,n1i*n2i*n3p,&
          mpidtypd,pot_ion,ngatherarr(0,1),&
          ngatherarr(0,2),mpidtypd,MPI_COMM_WORLD,ierr)
  else
     pot_ion => rho
  end if

  if (iproc == 0) then
     do i3=0,2*(n3+nbxz)-1
        do i2=0,2*(n2+nby)-1
           do i1=0,2*(n1+nbxz)-1
              ind=i1+nl1+(i2+nl2-1)*n1i+(i3+nl3-1)*n1i*n2i
              write(22,*)pot_ion(ind)
           end do
        end do
     end do
     close(22)
     !laterally averaged potential following the isolated direction in surfaces case
     !if (geocode == 'S') then
        open(unit=23,file=trim(filename)//'_avg'//geocode,status='unknown')
        do i2=0,2*n2+1
           later_avg=0.0_dp
           do i3=0,2*n3+1
              do i1=0,2*n1+1
                 ind=i1+nl1+(i2+nl2-1)*n1i+(i3+nl3-1)*n1i*n2i
                 later_avg=later_avg+pot_ion(ind)
              end do
           end do
           later_avg=later_avg/real((2*n1+2)*(2*n3+2),dp) !2D integration/2D Volume
           !problem with periodic/isolated BC
           write(23,*)i2,alat2/real(2*n2+2,dp)*i2,later_avg
        end do
        close(23)
     !end if
  end if

  if (nproc > 1) then
       i_all=-product(shape(pot_ion))*kind(pot_ion)
     deallocate(pot_ion,stat=i_stat)
     call memocc(i_stat,i_all,'pot_ion',subname)
  end if
END SUBROUTINE plot_density_old
