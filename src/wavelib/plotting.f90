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
  call initialize_work_arrays_sumrho(1,[lr],.true.,w)

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

!> Calculate the dipole of a Field given in the rho array.
!! The parallel distribution used is the one of the potential
subroutine calc_dipole(dpbox,nspin,at,rxyz,rho,calculate_quadropole)
  use module_base
  use module_dpbox, only: denspot_distribution
  use module_types
  use yaml_output
  implicit none
  integer, intent(in) :: nspin
  type(denspot_distribution), intent(in) :: dpbox
  type(atoms_data), intent(in) :: at
  real(gp), dimension(3,at%astruct%nat), intent(in) :: rxyz
  real(dp), dimension(dpbox%ndims(1),dpbox%ndims(2),max(dpbox%n3p, 1),nspin), target, intent(in) :: rho
  logical,intent(in) :: calculate_quadropole

  character(len=*), parameter :: subname='calc_dipole'
  integer :: ierr,n3p,nc1,nc2,nc3
  real(gp) :: q,qtot, delta_term,x,y,z,ri,rj
  integer  :: iat,i1,i2,i3, nl1,nl2,nl3, ispin,n1i,n2i,n3i, i, j
  real(gp), dimension(3) :: dipole_el,dipole_cores,tmpdip,charge_center_cores
  real(gp),dimension(3,nspin) :: charge_center_elec
  real(gp), dimension(3,3) :: quadropole_el,quadropole_cores,tmpquadrop
  real(dp), dimension(:,:,:,:), pointer :: ele_rho
!!$  real(dp), dimension(:,:,:,:), pointer :: rho_buf
  
  n1i=dpbox%ndims(1)
  n2i=dpbox%ndims(2)
  n3i=dpbox%ndims(3)
  n3p=dpbox%n3p

  if (dpbox%mpi_env%nproc > 1) then
     !allocate full density in pot_ion array
     ele_rho = f_malloc_ptr((/ n1i, n2i, n3i, nspin /),id='ele_rho')

!Commented out, it is enough to allocate the rho at 1
!!$     ! rho_buf is used instead of rho for avoiding the case n3p=0 in 
!!$     ! some procs which makes MPI_ALLGATHERV failed.
!!$     if (n3p.eq.0) then
!!$       allocate(rho_buf(n1i,n2i,n3p+1,nspin),stat=i_stat)
!!$       call memocc(i_stat,rho_buf,'rho_buf',subname)
!!$       rho_buf = 0.0_dp
!!$     else
!!$       allocate(rho_buf(n1i,n2i,n3p,nspin),stat=i_stat)
!!$       call memocc(i_stat,rho_buf,'rho_buf',subname)
!!$       rho_buf = rho
!!$     endif  


     do ispin=1,nspin
        call MPI_ALLGATHERV(rho(1,1,1,ispin),n1i*n2i*n3p,&
             mpidtypd,ele_rho(1,1,1,ispin),dpbox%ngatherarr(0,1),&
             dpbox%ngatherarr(0,2),mpidtypd,dpbox%mpi_env%mpi_comm,ierr)
     end do

  else
     ele_rho => rho
  end if

  if (at%astruct%geocode /= 'F') then
     nl1=1
     nl3=1
     nc1=n1i
     nc3=n3i
  else
     nl1=15
     nl3=15
     nc1=n1i-31
     nc3=n3i-31
  end if
  !value of the buffer in the y direction
  if (at%astruct%geocode == 'P') then
     nl2=1
     nc2=n2i
  else
     nl2=15
     nc2=n2i-31
  end if

  qtot=0.d0
  dipole_cores(1:3)=0_gp
  do iat=1,at%astruct%nat
     dipole_cores(1:3)=dipole_cores(1:3)+at%nelpsp(at%astruct%iatype(iat)) * rxyz(1:3,iat)
  end do

  dipole_el   (1:3)=0_gp
  do ispin=1,nspin
     do i3=0,nc3 - 1
        do i2=0,nc2 - 1
           do i1=0,nc1 - 1
              !ind=i1+nl1+(i2+nl2-1)*n1i+(i3+nl3-1)*n1i*n2i
              !q= ( ele_rho(ind,ispin) ) * hxh*hyh*hzh 
              q= - ele_rho(i1+nl1,i2+nl2,i3+nl3,ispin) * product(dpbox%hgrids)
              qtot=qtot+q
              dipole_el(1)=dipole_el(1)+ q* at%astruct%cell_dim(1)/real(nc1,dp)*i1 
              dipole_el(2)=dipole_el(2)+ q* at%astruct%cell_dim(2)/real(nc2,dp)*i2
              dipole_el(3)=dipole_el(3)+ q* at%astruct%cell_dim(3)/real(nc3,dp)*i3
           end do
        end do
     end do

  end do


  if (calculate_quadropole) then

      ! charge center
      charge_center_cores(1:3)=0.d0
      qtot=0.d0
      do iat=1,at%astruct%nat
          q=at%nelpsp(at%astruct%iatype(iat))
          charge_center_cores(1:3) = charge_center_cores(1:3) + q*at%astruct%rxyz(1:3,iat)
          qtot=qtot+q
      end do
      charge_center_cores=charge_center_cores/qtot


      quadropole_cores(1:3,1:3)=0._gp
      do iat=1,at%astruct%nat
          do i=1,3
              select case (i)
              case (1)
                  !ri=at%astruct%rxyz(1,iat)-charge_center_cores(1)
                  ri=at%astruct%rxyz(1,iat)
              case (2)
                  !ri=at%astruct%rxyz(2,iat)-charge_center_cores(2)
                  ri=at%astruct%rxyz(2,iat)
              case (3)
                  !ri=at%astruct%rxyz(3,iat)-charge_center_cores(3)
                  ri=at%astruct%rxyz(3,iat)
              case default
                  stop 'wrong value of i'
              end select
              do j=1,3
                  select case (j)
                  case (1)
                      !rj=at%astruct%rxyz(1,iat)-charge_center_cores(1)
                      rj=at%astruct%rxyz(1,iat)
                  case (2)
                      !rj=at%astruct%rxyz(2,iat)-charge_center_cores(2)
                      rj=at%astruct%rxyz(2,iat)
                  case (3)
                      !rj=at%astruct%rxyz(3,iat)-charge_center_cores(3)
                      rj=at%astruct%rxyz(3,iat)
                  case default
                      stop 'wrong value of j'
                  end select
                  if (i==j) then
                      !delta_term = (at%astruct%rxyz(1,iat)-charge_center_cores(1))**2 + &
                      !             (at%astruct%rxyz(2,iat)-charge_center_cores(2))**2 + &
                      !             (at%astruct%rxyz(3,iat)-charge_center_cores(3))**2
                      delta_term = at%astruct%rxyz(1,iat)**2 + &
                                   at%astruct%rxyz(2,iat)**2 + &
                                   at%astruct%rxyz(3,iat)**2
                  else
                      delta_term=0.d0
                  end if
                  q=at%nelpsp(at%astruct%iatype(iat))
                  quadropole_cores(j,i) = quadropole_cores(j,i) + q*(3.d0*rj*ri-delta_term)
              end do
          end do
      end do

      !!ele_rho=0.d0
      !!do iat=1,at%astruct%nat
      !!    i1=nint((at%astruct%rxyz(1,iat)/at%astruct%cell_dim(1))*real(nc1,dp))+nl1
      !!    i2=nint((at%astruct%rxyz(2,iat)/at%astruct%cell_dim(2))*real(nc2,dp))+nl2
      !!    i3=nint((at%astruct%rxyz(3,iat)/at%astruct%cell_dim(3))*real(nc3,dp))+nl3
      !!    if (bigdft_mpi%iproc==0) write(*,*) 'iat,i1,i2,i3',iat,i1,i2,i3
      !!    ele_rho(i1,i2,i3,1)=real(at%nelpsp(at%astruct%iatype(iat)))/product(dpbox%hgrids)
      !!end do


      ! charge center
      charge_center_elec(1:3,1:nspin)=0.d0
      do ispin=1,nspin
          qtot=0.d0
          do i3=0,nc3 - 1
              do i2=0,nc2 - 1
                  do i1=0,nc1 - 1
                      q= - ele_rho(i1+nl1,i2+nl2,i3+nl3,ispin) * product(dpbox%hgrids)
                      x=at%astruct%cell_dim(1)/real(nc1,dp)*i1
                      y=at%astruct%cell_dim(2)/real(nc2,dp)*i2
                      z=at%astruct%cell_dim(3)/real(nc3,dp)*i3
                      charge_center_elec(1,ispin) = charge_center_elec(1,ispin) + q*x
                      charge_center_elec(2,ispin) = charge_center_elec(2,ispin) + q*y
                      charge_center_elec(3,ispin) = charge_center_elec(3,ispin) + q*z
                      qtot=qtot+q
                  end do
              end do
          end do
          charge_center_elec(1:3,ispin)=charge_center_elec(1:3,ispin)/qtot
      end do

      quadropole_el(1:3,1:3)=0._gp
      do ispin=1,nspin
          do i3=0,nc3 - 1
              do i2=0,nc2 - 1
                  do i1=0,nc1 - 1
                      q= - ele_rho(i1+nl1,i2+nl2,i3+nl3,ispin) * product(dpbox%hgrids)
                      x=at%astruct%cell_dim(1)/real(nc1,dp)*i1
                      y=at%astruct%cell_dim(2)/real(nc2,dp)*i2
                      z=at%astruct%cell_dim(3)/real(nc3,dp)*i3
                      do i=1,3
                          select case (i)
                          case (1)
                              !ri=x-charge_center_cores(1)
                              ri=x+(charge_center_cores(1)-charge_center_elec(1,ispin))
                          case (2)
                              !ri=y-charge_center_cores(2)
                              ri=y+(charge_center_cores(2)-charge_center_elec(2,ispin))
                          case (3)
                              !ri=z-charge_center_cores(3)
                              ri=z+(charge_center_cores(3)-charge_center_elec(3,ispin))
                          case default
                              stop 'wrong value of i'
                          end select
                          do j=1,3
                              select case (j)
                              case (1)
                                  !rj=x-charge_center_cores(1)
                                  rj=x+(charge_center_cores(1)-charge_center_elec(1,ispin))
                              case (2)
                                  !rj=y-charge_center_cores(2)
                                  rj=y+(charge_center_cores(2)-charge_center_elec(2,ispin))
                              case (3)
                                  !rj=z-charge_center_cores(3)
                                  rj=z+(charge_center_cores(3)-charge_center_elec(3,ispin))
                              case default
                                  stop 'wrong value of j'
                              end select
                              if (i==j) then
                                  !delta_term = (x-charge_center_cores(1))**2 + &
                                  !             (y-charge_center_cores(2))**2 + &
                                  !             (z-charge_center_cores(3))**2
                                  delta_term = (x+(charge_center_cores(1)-charge_center_elec(1,ispin)))**2 + &
                                               (y+(charge_center_cores(2)-charge_center_elec(2,ispin)))**2 + &
                                               (z+(charge_center_cores(3)-charge_center_elec(3,ispin)))**2
                              else
                                  delta_term=0.d0
                              end if
                              quadropole_el(j,i) = quadropole_el(j,i) + q*(3.d0*rj*ri-delta_term)
                          end do
                      end do
                  end do
              end do
          end do
      end do

      tmpquadrop=quadropole_cores+quadropole_el

  end if

  if(bigdft_mpi%iproc==0) then
     !dipole_el=dipole_el        !/0.393430307_gp  for e.bohr to Debye2or  /0.20822678_gp  for e.A2Debye
     !dipole_cores=dipole_cores  !/0.393430307_gp  for e.bohr to Debye2or  /0.20822678_gp  for e.A2Debye
     !write(*,*) 'dipole_cores', dipole_cores
     !write(*,*) 'dipole_el', dipole_el
     tmpdip=dipole_cores+dipole_el
     call yaml_mapping_open('Electric Dipole Moment (AU)')
       call yaml_map('P vector',tmpdip(1:3),fmt='(1pe13.4)')
       call yaml_map('norm(P)',sqrt(sum(tmpdip**2)),fmt='(1pe14.6)')
     call yaml_mapping_close()
     tmpdip=tmpdip/0.393430307_gp  ! au2debye              
     call yaml_mapping_open('Electric Dipole Moment (Debye)')
       call yaml_map('P vector',tmpdip(1:3),fmt='(1pe13.4)')
       call yaml_map('norm(P)',sqrt(sum(tmpdip**2)),fmt='(1pe14.6)')
     call yaml_mapping_close()


      if (calculate_quadropole) then
          !call yaml_sequence_open('core quadropole')
          !do i=1,3
          !   call yaml_sequence(trim(yaml_toa(quadropole_cores(i,1:3),fmt='(es15.8)')))
          !end do
          !call yaml_sequence_close()

          !call yaml_sequence_open('electronic quadropole')
          !do i=1,3
          !   call yaml_sequence(trim(yaml_toa(quadropole_el(i,1:3),fmt='(es15.8)')))
          !end do
          !call yaml_sequence_close()

          call yaml_sequence_open('Quadrupole Moment (AU)')
          do i=1,3
             call yaml_sequence(trim(yaml_toa(tmpquadrop(i,1:3),fmt='(es15.8)')))
          end do
          call yaml_map('trace',tmpquadrop(1,1)+tmpquadrop(2,2)+tmpquadrop(3,3),fmt='(es12.2)')
          call yaml_sequence_close()
      end if

  endif

  if (dpbox%mpi_env%nproc > 1) then
     call f_free_ptr(ele_rho)
  else
     nullify(ele_rho)
  end if

END SUBROUTINE calc_dipole

subroutine plot_density(iproc,nproc,filename,at,rxyz,kernel,nspin,rho)
  use module_defs, only: gp,dp
  use PStypes, only: coulomb_operator
  use IObox, only: dump_field
  use PSbox, only: PS_gather
  use module_atoms, only: atoms_data
  use dynamic_memory
  implicit none
  integer, intent(in) :: iproc,nproc,nspin
  type(atoms_data), intent(in) :: at
  type(coulomb_operator), intent(in) :: kernel
  character(len=*), intent(in) :: filename
  real(gp), dimension(3,at%astruct%nat), intent(in) :: rxyz
  real(dp), dimension(*) :: rho !< intent(in)
  !local variables
  integer :: ispin
  real(dp), dimension(:,:,:,:), allocatable :: pot_ion

  pot_ion = f_malloc((/kernel%ndims(1),kernel%ndims(2),kernel%ndims(3), nspin /),id='pot_ion')
  if (nproc > 1) then

     !here we might add an extra interface
     call PS_gather(src=rho,dest=pot_ion,kernel=kernel,nsrc=nspin)
  else
     call f_memcpy(n=size(pot_ion),src=rho(1),dest=pot_ion(1,1,1,1))
  end if

  call dump_field(filename,at%astruct%geocode,kernel%ndims,kernel%hgrids,nspin,pot_ion,&
       rxyz,at%astruct%iatype,at%nzatom,at%nelpsp)

  call f_free(pot_ion)

!!$  character(len=*), parameter :: subname='plot_density'
!!$  character(len=5) :: suffix
!!$  character(len=65) :: message
!!$  integer :: ierr,ia,ib,isuffix,fformat,n1i,n2i,n3i
!!$  real(dp) :: a,b
!!$  real(gp) :: hxh,hyh,hzh
!!$  real(dp), dimension(:,:), pointer :: pot_ion
!!$  integer,parameter :: unit0 = 22
!!$  integer,parameter :: unitx = 23
!!$  integer,parameter :: unity = 24
!!$  integer,parameter :: unitz = 25
!!$
!!$  n1i=box%ndims(1)
!!$  n2i=box%ndims(2)
!!$  n3i=box%ndims(3)
!!$
!!$  hxh=box%hgrids(1)
!!$  hyh=box%hgrids(2)
!!$  hzh=box%hgrids(3)
!!$
!!$  if (nproc > 1) then
!!$     !allocate full density in pot_ion array
!!$     pot_ion = f_malloc_ptr((/ box%ndimgrid, nspin /),id='pot_ion')
!!$     
!!$     call mpiallgather(sendbuf=rho(1,1),sendcount=box%ndimpot,&
!!$          recvbuf=pot_ion(1,1),recvcounts=box%ngatherarr(:,1),&
!!$          displs=box%ngatherarr(:,2),comm=box%mpi_env%mpi_comm)
!!$
!!$     !case for npspin==2
!!$     if (nspin==2) then
!!$        call mpiallgather(sendbuf=rho(1,2),sendcount=box%ndimpot,&
!!$             recvbuf=pot_ion(1,2),recvcounts=box%ngatherarr(:,1),&
!!$             displs=box%ngatherarr(:,2),comm=box%mpi_env%mpi_comm)
!!$     end if
!!$
!!$  else
!!$     !pot_ion => rho
!!$     pot_ion = f_malloc_ptr(shape(rho),id='pot_ion')
!!$     call f_memcpy(dest=pot_ion,src=rho)
!!$  end if
!!$
!!$  ! Format = 1 -> cube (default)
!!$  ! Format = 2 -> ETSF
!!$  ! ...
!!$  fformat = 1
!!$  isuffix = index(filename, ".cube", back = .true.)
!!$  if (isuffix > 0) then
!!$     isuffix = isuffix - 1
!!$     fformat = 1
!!$  else
!!$     isuffix = index(filename, ".etsf", back = .true.)
!!$     if (isuffix <= 0) isuffix = index(filename, ".etsf.nc", back = .true.)
!!$     if (isuffix > 0) then
!!$        isuffix = isuffix - 1
!!$        fformat = 2
!!$     else
!!$        isuffix = len(trim(filename))
!!$     end if
!!$  end if
!!$  if (iproc == 0) then
!!$
!!$
!!$     open(unit=unit0,file=trim(filename(:isuffix))//'.cube',status='unknown')
!!$     open(unit=unitx,file=trim(filename(:isuffix))//'_avg_x',status='unknown')
!!$     open(unit=unity,file=trim(filename(:isuffix))//'_avg_y',status='unknown')
!!$     open(unit=unitz,file=trim(filename(:isuffix))//'_avg_z',status='unknown')
!!$
!!$     if (nspin /=2) then
!!$        message='total spin'
!!$        if (fformat == 1) then
!!$           suffix=''
!!$           a=1.0_dp
!!$           ia=1
!!$           b=0.0_dp
!!$           ib=1
!!$           call write_cube_fields(unit0,unitx,unity,unitz,message,&
!!$                at,1.d0,rxyz,n1i,n2i,n3i,0,0,0,hxh,hyh,hzh,&
!!$                a,pot_ion(1,ia),1,b,pot_ion(1,ib))
!!$        else
!!$           call write_etsf_density(filename(:isuffix),message,&
!!$                at,rxyz,n1i,n2i,n3i,hxh,hyh,hzh,&
!!$                pot_ion, 1)
!!$        end if
!!$     else
!!$        if (fformat == 1) then
!!$           suffix=''
!!$           message='total spin'
!!$           a=1.0_dp
!!$           ia=1
!!$           b=0.0_dp
!!$           ib=2
!!$           call write_cube_fields(unit0,unitx,unity,unitz,message,&
!!$                at,1.d0,rxyz,n1i,n2i,n3i,0,0,0,hxh,hyh,hzh,&
!!$                a,pot_ion(1,ia),1,b,pot_ion(1,ib))
!!$
!!$           suffix='-down'
!!$           message='spin down'
!!$           a=0.0_dp
!!$           ia=1
!!$           b=1.0_dp
!!$           ib=2
!!$           call write_cube_fields(unit0,unitx,unity,unitz,message,&
!!$                at,1.d0,rxyz,n1i,n2i,n3i,0,0,0,hxh,hyh,hzh,&
!!$                a,pot_ion(1,ia),1,b,pot_ion(1,ib))
!!$
!!$           suffix='-u-d'
!!$           message='spin difference'
!!$           a=1.0_dp
!!$           ia=1
!!$           b=-2.0_dp
!!$           ib=2
!!$           call write_cube_fields(unit0,unitx,unity,unitz,message,&
!!$                at,1.d0,rxyz,n1i,n2i,n3i,0,0,0,hxh,hyh,hzh,&
!!$                a,pot_ion(1,ia),1,b,pot_ion(1,ib))
!!$
!!$           suffix='-up'
!!$           message='spin up'
!!$           a=1.0_dp
!!$           ia=1
!!$           b=-1.0_dp
!!$           ib=2
!!$           call write_cube_fields(unit0,unitx,unity,unitz,message,&
!!$                at,1.d0,rxyz,n1i,n2i,n3i,0,0,0,hxh,hyh,hzh,&
!!$                a,pot_ion(1,ia),1,b,pot_ion(1,ib))
!!$        else
!!$           message = 'spin up, down, total, difference'
!!$           call write_etsf_density(filename(:isuffix),message,&
!!$                at,rxyz,n1i,n2i,n3i,hxh,hyh,hzh,&
!!$                pot_ion, 2)
!!$        end if
!!$
!!$     end if
!!$
!!$     close(unit=unit0)
!!$     close(unit=unitx)
!!$     close(unit=unity)
!!$     close(unit=unitz)
!!$
!!$  end if
!!$
!!$
!!$  !if (nproc > 1) then
!!$     call f_free_ptr(pot_ion)
!!$  !end if

END SUBROUTINE plot_density
