!> @file
!!  Routines to plot wavefunctions
!! @author
!!    Copyright (C) 2010-2015 BigDFT group 
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 
subroutine plot_density_cube_old(filename,iproc,nproc,n1,n2,n3,n1i,n2i,n3i,n3p,nspin,&
     hxh,hyh,hzh,at,rxyz,ngatherarr,rho)
  use module_base
  use module_types
  implicit none
  character(len=*), intent(in) :: filename
  integer, intent(in) :: iproc,n1i,n2i,n3i,n3p,n1,n2,n3,nspin,nproc
  real(gp), intent(in) :: hxh,hyh,hzh
  type(atoms_data), intent(in) :: at
  integer, dimension(0:nproc-1,2), intent(in) :: ngatherarr
  real(gp), dimension(3,at%astruct%nat), intent(in) :: rxyz
  real(dp), dimension(n1i*n2i*n3p,nspin), target, intent(in) :: rho
  !local variables
  character(len=*), parameter :: subname='plot_density_cube_old'
  character(len=3) :: advancestring
  character(len=5) :: suffix
  character(len=15) :: message
  integer :: nl1,nl2,nl3,i1,i2,i3,ind,ierr,icount,j,iat,ia,ib,nbxz,nby
  real(dp) :: a,b
  real(dp), dimension(:,:), pointer :: pot_ion

  !conditions for periodicity in the three directions
  !value of the buffer in the x and z direction
  if (at%astruct%geocode /= 'F') then
     nl1=1
     nl3=1
     nbxz = 1
  else
     nl1=15
     nl3=15
     nbxz = 0
  end if
  !value of the buffer in the y direction
  if (at%astruct%geocode == 'P') then
     nl2=1
     nby = 1
  else
     nl2=15
     nby = 0
  end if

  if (nproc > 1) then
     !allocate full density in pot_ion array
     pot_ion = f_malloc_ptr((/ n1i*n2i*n3i, nspin /),id='pot_ion')

     call MPI_ALLGATHERV(rho(1,1),n1i*n2i*n3p,&
          mpidtypd,pot_ion(1,1),ngatherarr(0,1),&
          ngatherarr(0,2),mpidtypd,bigdft_mpi%mpi_comm,ierr)

     !case for npspin==2
     if (nspin==2) then
        call MPI_ALLGATHERV(rho(1,2),n1i*n2i*n3p,&
             mpidtypd,pot_ion(1,2),ngatherarr(0,1),&
             ngatherarr(0,2),mpidtypd,bigdft_mpi%mpi_comm,ierr)
     end if

  else
     pot_ion => rho
  end if

  if (iproc == 0) then

     if (nspin /=2) then
        suffix=''
        message='total spin'
        a=1.0_dp
        ia=1
        b=0.0_dp
        ib=1
        call cubefile_write
     else
        suffix='-up'
        message='spin up'
        a=1.0_dp
        ia=1
        b=0.0_dp
        ib=2
        call cubefile_write

        suffix='-down'
        message='spin down'
        a=0.0_dp
        ia=1
        b=1.0_dp
        ib=2
        call cubefile_write

        suffix=''
        message='total spin'
        a=1.0_dp
        ia=1
        b=1.0_dp
        ib=2
        call cubefile_write

        suffix='-u-d'
        message='spin difference'
        a=1.0_dp
        ia=1
        b=-1.0_dp
        ib=2
        call cubefile_write
     end if
  end if

  if (nproc > 1) then
     call f_free_ptr(pot_ion)
  end if

contains

  subroutine cubefile_write
        open(unit=22,file=trim(filename)//trim(suffix)//'.cube',status='unknown')
        write(22,*)'CUBE file for charge density'
        write(22,*)'Case for '//trim(message)
        !number of atoms
!        if (geocode=='P') then
           write(22,'(i5,3(f12.6),a)') at%astruct%nat,0.0_gp,0.0_gp,0.0_gp,' modified'
!        else if (geocode=='S') then
!           write(22,'(i5,3(f12.6))') at%astruct%nat,0.0_gp,-hyh,0.0_gp
!        else if (geocode=='F') then
!           write(22,'(i5,3(f12.6))') at%astruct%nat,-hxh,-hyh,-hzh
!        end if
        !grid and grid spacings
        write(22,'(i5,3(f12.6))') 2*(n1+nbxz),hxh,0.0_gp,0.0_gp
        write(22,'(i5,3(f12.6))') 2*(n2+nby) ,0.0_gp,hyh,0.0_gp
        write(22,'(i5,3(f12.6))') 2*(n3+nbxz),0.0_gp,0.0_gp,hzh
        !atomic number and positions
        do iat=1,at%astruct%nat
           write(22,'(i5,4(f12.6))') at%nzatom(at%astruct%iatype(iat)),0.0_gp,(rxyz(j,iat),j=1,3)
        end do

        !the loop is reverted for a cube file
        !charge normalised to the total charge
        do i1=0,2*(n1+nbxz)-1
           do i2=0,2*(n2+nby)-1
              icount=0
              do i3=0,2*(n3+nbxz)-1
                 icount=icount+1
                 if (icount == 6 .or. i3==2*(n3+nbxz) - 1) then
                    advancestring='yes'
                    icount=0
                 else
                    advancestring='no'
                 end if

                 ind=i1+nl1+(i2+nl2-1)*n1i+(i3+nl3-1)*n1i*n2i
                 write(22,'(1x,1pe13.6)',advance=advancestring)&
                      a*pot_ion(ind,ia)+b*pot_ion(ind,ib)
              end do
           end do
        end do
        close(22)
  END SUBROUTINE cubefile_write

END SUBROUTINE plot_density_cube_old


subroutine read_density_cube_old(filename, n1i,n2i,n3i, nspin, hxh,hyh,hzh, nat, rxyz,  rho)
  use module_base
  use module_types
  implicit none
  character(len=*), intent(in) :: filename
  integer, intent(out) ::  n1i,n2i,n3i
  integer, intent(in) :: nspin
  real(gp), intent(out) :: hxh,hyh,hzh
  real(gp), pointer :: rxyz(:,:)
  real(dp), dimension(:), pointer :: rho
  integer, intent(out) ::  nat
 
  !local variables
  character(len=*), parameter :: subname='read_density_cube_old'
  character(len=5) :: suffix
  character(len=15) :: message
  character(len=3) :: advancestring
  integer :: i1,i2,i3,ind,icount,j,iat,ia

  if (nspin /=2) then
     suffix=''
     message='total spin'
     ia=1
     call cubefile_read
  else
     suffix='-up'
     message='spin up'
     ia=1
     call cubefile_read
     
     suffix='-down'
     message='spin down'
     ia=2
     call cubefile_read
     
  end if

contains

  subroutine cubefile_read
       real(dp) dum1,dum2, dum3
        integer idum
        open(unit=22,file=trim(filename)//trim(suffix)//'.cube',status='old')
        read(22,*)! 'CUBE file for charge density'
        read(22,*)! 'Case for '//trim(message)

        read(22,'(i5,3(f12.6),a)')  nat , dum1, dum2, dum3 ! ,0.0_gp,0.0_gp,0.0_gp,' modified'
           

        read(22,'(i5,3(f12.6))') n1i , hxh,   dum1 ,dum2
        read(22,'(i5,3(f12.6))') n2i ,dum1 , hyh  ,  dum2
        read(22,'(i5,3(f12.6))') n3i ,dum1 , dum2 , hzh
        !atomic number and positions

        if( associated(rxyz) ) then
           call f_free_ptr(rxyz)
        end if
        
        rxyz = f_malloc_ptr((/ 3, nat /),id='rxyz')
        

        if( associated(rho ).and. ia==1 ) then
           call f_free_ptr(rho)
        end if
        if(ia==1) then
           rho = f_malloc_ptr(n1i*n2i*n3i,id='rho')
        endif

        do iat=1,nat
           read(22,'(i5,4(f12.6))') idum , dum1 , (rxyz(j,iat),j=1,3)
           ! write(22,'(i5,4(f12.6))') at%nzatom(at%astruct%iatype(iat)),0.0_gp,(rxyz(j,iat),j=1,3)
        end do

        do i1=1,n1i
           do i2=1,n2i
              icount=0
              do i3=1,n3i
                 icount=icount+1
                 if (icount == 6  .or. i3==n3i  ) then
                    advancestring='yes'
                    icount=0
                 else
                    advancestring='no'
                 end if

                 ind=i1+(i2-1)*n1i+(i3-1)*n1i*n2i

                 read(22,'(1x,1pe13.6)',advance=advancestring)  rho(ind) 
              end do
           end do
        end do
        close(22)
  END SUBROUTINE cubefile_read

END SUBROUTINE read_density_cube_old


!> the API of this function has to be modified to make it readable
subroutine plot_wf(units_provided,orbname,nexpo,at,factor,lr,hx,hy,hz,rxyz,psi, &
           unit0_, unitx_, unity_, unitz_)
  use module_base
  use locregs, only: locreg_descriptors
  use module_types, only: atoms_data
  use locreg_operations
  use IObox, only: dump_field
  implicit none
  !Arguments
  logical,intent(in) :: units_provided
  character(len=*) :: orbname
  integer, intent(in) :: nexpo
  real(dp), intent(in) :: factor
  real(gp), intent(in) :: hx,hy,hz
  type(atoms_data), intent(in) :: at
  real(gp), dimension(3,at%astruct%nat), intent(in) :: rxyz
  type(locreg_descriptors), intent(in) :: lr
  real(wp), dimension(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f), intent(in) :: psi
  integer,intent(in),optional :: unit0_, unitx_, unity_, unitz_
  !Local variables
  character(len=*), parameter :: subname='plot_wf'
  integer :: n1i,n2i,n3i,n1,n2,n3,n1s,n2s,n3s
  type(workarr_sumrho) :: w
  real(wp), dimension(:), allocatable :: psir
  integer, dimension(3) :: ndims
  real(gp), dimension(3) :: hgrids
  integer :: iunit0, iunitx, iunity, iunitz
  integer,parameter :: unit0 = 22
  integer,parameter :: unitx = 23
  integer,parameter :: unity = 24
  integer,parameter :: unitz = 25

  ! Check whether all units are present
  if (units_provided) then
      if (.not.present(unit0_)) call f_err_throw('unit0_ is not present',err_name='BIGDFT_RUNTIME_ERROR')
      if (.not.present(unitx_)) call f_err_throw('unitx_ is not present',err_name='BIGDFT_RUNTIME_ERROR')
      if (.not.present(unity_)) call f_err_throw('unity_ is not present',err_name='BIGDFT_RUNTIME_ERROR')
      if (.not.present(unitz_)) call f_err_throw('unitz_ is not present',err_name='BIGDFT_RUNTIME_ERROR')
  end if

  n1=lr%d%n1
  n2=lr%d%n2
  n3=lr%d%n3
  n1s=lr%ns1
  n2s=lr%ns2
  n3s=lr%ns3
  n1i=lr%d%n1i
  n2i=lr%d%n2i
  n3i=lr%d%n3i


  ndims(1)=n1i
  ndims(2)=n2i
  ndims(3)=n3i
  hgrids(1)=0.5_gp*hx
  hgrids(2)=0.5_gp*hy
  hgrids(3)=0.5_gp*hz
  call initialize_work_arrays_sumrho(lr,.true.,w)

  psir = f_malloc(lr%d%n1i*lr%d%n2i*lr%d%n3i,id='psir')
  !initialisation
  if (lr%geocode == 'F') call f_zero(psir)

  call daub_to_isf(lr,w,psi,psir)

  if (units_provided) then
      iunit0 = unit0_
      iunitx = unitx_
      iunity = unity_
      iunitz = unitz_
  else
      iunit0 = unit0
      iunitx = unitx
      iunity = unity
      iunitz = unitz
      !open(unit=iunit0,file=trim(orbname)//'.cube',status='unknown')
      !open(unit=iunitx,file=trim(orbname)//'_avg_x',status='unknown')
      !open(unit=iunity,file=trim(orbname)//'_avg_y',status='unknown')
      !open(unit=iunitz,file=trim(orbname)//'_avg_z',status='unknown')
  end if
!!$
!!$  call write_cube_fields(iunit0,iunitx,iunity,iunitz,' ',&
!!$       at,factor,rxyz,n1i,n2i,n3i,n1s,n2s,n3s,0.5_gp*hx,0.5_gp*hy,0.5_gp*hz,&
!!$       1.0_gp,psir,nexpo,0.0_gp,psir)

  call dump_field(trim(orbname)//'.cube',lr%geocode,ndims,hgrids,1,psir,&
       rxyz,at%astruct%iatype,at%nzatom,at%nelpsp)


  if (units_provided) then
      !close(unit=iunit0)
      !close(unit=iunitx)
      !close(unit=iunity)
      !close(unit=iunitz)
  end if

  call f_free(psir)

  call deallocate_work_arrays_sumrho(w)

END SUBROUTINE plot_wf


