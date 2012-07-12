!> @file
!!  Routines to create descriptor arrays for density and potential
!! @author
!!    Copyright (C) 2007-2011 BigDFT group (LG)
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 

!> Denspot initialization
subroutine initialize_DFT_local_fields(denspot)
  use module_base
  use module_types
  implicit none
  type(DFT_local_fields), intent(inout) :: denspot

  denspot%rhov_is = EMPTY
  nullify(denspot%rho_C,denspot%V_ext,denspot%Vloc_KS,denspot%rho_psi)
  nullify(denspot%V_XC)
  nullify(denspot%f_XC,denspot%rho_work,denspot%pot_work,denspot%rhov)

  denspot%psoffset=0.0_gp

  if (verbose >1) then
     denspot%PSquiet='NO '
  else
     denspot%PSquiet='YES'
  end if

  call initialize_coulomb_operator(denspot%pkernel)
  call initialize_coulomb_operator(denspot%pkernelseq)
  call initialize_rho_descriptors(denspot%rhod)
  call initialize_denspot_distribution(denspot%dpbox)

  nullify(denspot%mix)
end subroutine initialize_DFT_local_fields

subroutine initialize_coulomb_operator(kernel)
  use module_base
  use module_types
  implicit none
  type(coulomb_operator), intent(out) :: kernel
  
  nullify(kernel%kernel)

  
end subroutine initialize_coulomb_operator


subroutine initialize_denspot_distribution(dpbox)
  use module_base
  use module_types
  implicit none
  type(denspot_distribution), intent(out) :: dpbox
  !local variables
  integer :: i
  
  dpbox%n3d      =uninitialized(dpbox%n3d)      
  dpbox%n3p      =uninitialized(dpbox%n3p)      
  dpbox%n3pi     =uninitialized(dpbox%n3pi)     
  dpbox%i3xcsh   =uninitialized(dpbox%i3xcsh)   
  dpbox%i3s      =uninitialized(dpbox%i3s)      
  dpbox%nrhodim  =uninitialized(dpbox%nrhodim)  
  dpbox%i3rho_add=uninitialized(dpbox%i3rho_add)
  do i=1,3
     dpbox%hgrids(i)=uninitialized(dpbox%hgrids(i))
     dpbox%ndims(i)=uninitialized(dpbox%ndims(i))
  end do

  nullify(dpbox%nscatterarr,dpbox%ngatherarr)
  dpbox%iproc=0
  dpbox%iproc_world=0
  dpbox%nproc=1
  dpbox%mpi_comm=MPI_COMM_WORLD
  
end subroutine initialize_denspot_distribution

subroutine initialize_rho_descriptors(rhod)
  use module_base
  use module_types
  implicit none
  type(rho_descriptors), intent(out) :: rhod

  rhod%geocode='X' !fake value
  rhod%icomm=1 !< lda case
  rhod%nrhotot=uninitialized(rhod%nrhotot)
  rhod%n_csegs=uninitialized(rhod%n_csegs)
  rhod%n_fsegs=uninitialized(rhod%n_fsegs)
  rhod%dp_size=uninitialized(rhod%dp_size)
  rhod%sp_size=uninitialized(rhod%sp_size)
  
  nullify(rhod%spkey,rhod%dpkey,rhod%cseg_b,rhod%fseg_b)

end subroutine initialize_rho_descriptors

subroutine dpbox_set_box(dpbox,Lzd)
  use module_base
  use module_types
  implicit none
  type(local_zone_descriptors), intent(in) :: Lzd
  type(denspot_distribution), intent(inout) :: dpbox
  
  dpbox%hgrids(1)=0.5_gp*Lzd%hgrids(1)
  dpbox%hgrids(2)=0.5_gp*Lzd%hgrids(2)
  dpbox%hgrids(3)=0.5_gp*Lzd%hgrids(3)
  dpbox%ndims(1)=Lzd%Glr%d%n1i
  dpbox%ndims(2)=Lzd%Glr%d%n2i
  dpbox%ndims(3)=Lzd%Glr%d%n3i

end subroutine dpbox_set_box

!>todo: remove n1i and n2i
subroutine denspot_set_history(denspot, iscf, nspin, &
     & n1i, n2i) !to be removed arguments when denspot has dimensions
  use module_base
  use module_types
  use m_ab6_mixing
  implicit none
  type(DFT_local_fields), intent(inout) :: denspot
  integer, intent(in) :: iscf, n1i, n2i, nspin
  
  integer :: potden, npoints, ierr
  character(len=500) :: errmess

  if (iscf < 10) then
     potden = AB6_MIXING_POTENTIAL
     npoints = n1i*n2i*denspot%dpbox%n3p
     if (denspot%dpbox%n3p==0) npoints=1
  else
     potden = AB6_MIXING_DENSITY
     npoints = n1i*n2i*denspot%dpbox%n3d
     if (denspot%dpbox%n3d==0) npoints=1
  end if
  if (iscf > SCF_KIND_DIRECT_MINIMIZATION) then
     allocate(denspot%mix)
     call ab6_mixing_new(denspot%mix, modulo(iscf, 10), potden, &
          AB6_MIXING_REAL_SPACE, npoints, nspin, 0, &
          ierr, errmess, useprec = .false.)
     call ab6_mixing_eval_allocate(denspot%mix)
  else
     nullify(denspot%mix)
  end if
end subroutine denspot_set_history


subroutine denspot_communications(iproc_world,nproc_world,iproc,nproc,mpi_comm,&
     ixc,nspin,geocode,SICapproach,dpbox)
  use module_base
  use module_types
  use module_interfaces, except_this_one => denspot_communications
  implicit none
  integer, intent(in) :: iproc_world,iproc,nproc,mpi_comm,ixc,nspin,nproc_world
  character(len=1), intent(in) :: geocode
  character(len=4), intent(in) :: SICapproach
  type(denspot_distribution), intent(inout) :: dpbox
  !local variables
  character(len = *), parameter :: subname = 'denspot_communications' 
  integer :: i_stat

  dpbox%iproc_world=iproc_world
  dpbox%nproc_world=nproc_world
  dpbox%iproc=iproc
  dpbox%nproc=nproc
  dpbox%mpi_comm=mpi_comm

  ! Create descriptors for density and potentials.
  ! ------------------
  !these arrays should be included in the comms descriptor
  !allocate values of the array for the data scattering in sumrho
  !its values are ignored in the datacode='G' case
  allocate(dpbox%nscatterarr(0:nproc-1,4+ndebug),stat=i_stat)
  call memocc(i_stat,dpbox%nscatterarr,'nscatterarr',subname)
  !allocate array for the communications of the potential
  !also used for the density
  allocate(dpbox%ngatherarr(0:nproc-1,3+ndebug),stat=i_stat)
  call memocc(i_stat,dpbox%ngatherarr,'ngatherarr',subname)

  call dpbox_repartition(iproc,nproc,geocode,'D',ixc,dpbox)

!!$  !create the descriptors for the density and the potential
!!$  !these descriptors should take into account the localisation regions
!!$  call createDensPotDescriptors(iproc,nproc,atoms,grid,hxh,hyh,hzh, &
!!$       rxyz,in%crmult,in%frmult,radii_cf,in%nspin,'D',in%ixc,in%rho_commun, &
!!$       dpbox%n3d,dpbox%n3p,&
!!$       dpbox%n3pi,dpbox%i3xcsh,dpbox%i3s, &
!!$       dpbox%nscatterarr,dpbox%ngatherarr,rhod)

  !Allocate Charge density / Potential in real space
  !here the full_density treatment should be put
  dpbox%nrhodim=nspin
  dpbox%i3rho_add=0
  if (trim(SICapproach)=='NK') then
     dpbox%nrhodim=2*dpbox%nrhodim !to be eliminated with a orbital-dependent potential
     dpbox%i3rho_add=dpbox%ndims(1)*dpbox%ndims(2)*dpbox%i3xcsh+1
  end if

  !fill the full_local_potential dimension
  dpbox%ndimpot=dpbox%ndims(1)*dpbox%ndims(2)*dpbox%n3p
  dpbox%ndimgrid=dpbox%ndims(1)*dpbox%ndims(2)*dpbox%ndims(3)
  dpbox%ndimrhopot=dpbox%ndims(1)*dpbox%ndims(2)*dpbox%n3d*&
       dpbox%nrhodim
end subroutine denspot_communications

subroutine denspot_set_rhov_status(denspot, status, istep, iproc, nproc)
  use module_base
  use module_types
  implicit none
  type(DFT_local_fields), intent(inout) :: denspot
  integer, intent(in) :: status, istep, iproc, nproc

  denspot%rhov_is = status
  
  if (denspot%c_obj /= 0) then
     call denspot_emit_rhov(denspot, istep, iproc, nproc)
  end if
end subroutine denspot_set_rhov_status

subroutine denspot_full_density(denspot, rho_full, iproc, new)
  use module_base
  use module_types
  use m_profiling
  implicit none
  type(DFT_local_fields), intent(in) :: denspot
  integer, intent(in) :: iproc
  integer, intent(out) :: new
  real(gp), dimension(:), pointer :: rho_full

  character(len = *), parameter :: subname = "denspot_full_density"
  integer :: i_stat, nslice, ierr, irhodim, irhoxcsh

  new = 0
  nslice = max(denspot%dpbox%ndimpot, 1)
  if (nslice < denspot%dpbox%ndimgrid) then
     if (iproc == 0) then
        !allocate full density in pot_ion array
        allocate(rho_full(denspot%dpbox%ndimgrid*denspot%dpbox%nrhodim+ndebug),stat=i_stat)
        call memocc(i_stat,rho_full,'rho_full',subname)
        new = 1
        
        ! Ask to gather density to other procs.
        call MPI_BCAST(0, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
     end if

     if (denspot%dpbox%ndimrhopot > 0) then
        irhoxcsh = nslice / denspot%dpbox%n3p * denspot%dpbox%i3xcsh
     else
        irhoxcsh = 0
     end if     
     do irhodim = 1, denspot%dpbox%nrhodim, 1
        if (iproc == 0) then
           call MPI_GATHERV(denspot%rhov(nslice * (irhodim - 1) + irhoxcsh + 1),&
                nslice,mpidtypd,rho_full(denspot%dpbox%ndimgrid * (irhodim - 1) + 1),&
                denspot%dpbox%ngatherarr(0,1),denspot%dpbox%ngatherarr(0,2),&
                mpidtypd,0,MPI_COMM_WORLD,ierr)
        else
           call MPI_GATHERV(denspot%rhov(nslice * (irhodim - 1) + irhoxcsh + 1),&
                nslice,mpidtypd,rho_full(1),&
                denspot%dpbox%ngatherarr(0,1),denspot%dpbox%ngatherarr(0,2),&
                mpidtypd,0,MPI_COMM_WORLD,ierr)
        end if
     end do
  else
     rho_full => denspot%rhov
  end if
END SUBROUTINE denspot_full_density
subroutine denspot_full_v_ext(denspot, pot_full, iproc, new)
  use module_base
  use module_types
  use m_profiling
  implicit none
  type(DFT_local_fields), intent(in) :: denspot
  integer, intent(in) :: iproc
  integer, intent(out) :: new
  real(gp), pointer :: pot_full(:)

  character(len = *), parameter :: subname = "localfields_full_potential"
  integer :: i_stat, ierr

  new = 0
  if (denspot%dpbox%ndimpot < denspot%dpbox%ndimgrid) then
     if (iproc == 0) then
        !allocate full density in pot_ion array
        allocate(pot_full(denspot%dpbox%ndimgrid+ndebug),stat=i_stat)
        call memocc(i_stat,pot_full,'pot_full',subname)
        new = 1
      
        ! Ask to gather density to other procs.
        call MPI_BCAST(1, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
     end if

     call MPI_GATHERV(denspot%v_ext(1,1,1,1),max(denspot%dpbox%ndimpot, 1),&
          mpidtypd,pot_full(1),denspot%dpbox%ngatherarr(0,1),&
          denspot%dpbox%ngatherarr(0,2),mpidtypd,0,MPI_COMM_WORLD,ierr)
  else
     pot_full => denspot%rhov
  end if
END SUBROUTINE denspot_full_v_ext
subroutine denspot_emit_rhov(denspot, iter, iproc, nproc)
  use module_base
  use module_types
  implicit none
  type(DFT_local_fields), intent(in) :: denspot
  integer, intent(in) :: iter, iproc, nproc

  character(len = *), parameter :: subname = "denspot_emit_rhov"
  integer, parameter :: SIGNAL_DONE = -1
  integer, parameter :: SIGNAL_DENSITY = 0
  integer :: message, ierr, i_stat, i_all, new
  real(gp), pointer :: full_dummy(:)
  interface
     subroutine denspot_full_density(denspot, rho_full, iproc, new)
       use module_types
       implicit none
       type(DFT_local_fields), intent(in) :: denspot
       integer, intent(in) :: iproc
       integer, intent(out) :: new
       real(gp), dimension(:), pointer :: rho_full
     END SUBROUTINE denspot_full_density
  end interface

  call timing(iproc,'rhov_signals  ','ON')
  if (iproc == 0) then
     ! Only iproc 0 emit the signal. This call is blocking.
     ! All other procs are blocked by the bcast to wait for
     ! possible transfer to proc 0.
     call localfields_emit_rhov(denspot%c_obj, iter)
     if (nproc > 1) then
        ! After handling the signal, iproc 0 broadcasts to other
        ! proc to continue (jproc == -1).
        message = SIGNAL_DONE
        call MPI_BCAST(message, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
     end if
  else
     do
        call MPI_BCAST(message, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
        if (message == SIGNAL_DONE) then
           exit
        else if (message == SIGNAL_DENSITY) then
           allocate(full_dummy(denspot%dpbox%nrhodim+ndebug),stat=i_stat)
           call memocc(i_stat,full_dummy,'full_dummy',subname)
           ! Gather density to iproc 0
           call denspot_full_density(denspot, full_dummy, iproc, new)
           i_all=-product(shape(full_dummy))*kind(full_dummy)
           deallocate(full_dummy,stat=i_stat)
           call memocc(i_stat,i_all,'full_dummy',subname)
        end if
     end do
  end if
  call timing(iproc,'rhov_signals  ','OF')
END SUBROUTINE denspot_emit_rhov
subroutine denspot_emit_v_ext(denspot, iproc, nproc)
  use module_base
  use module_types
  implicit none
  type(DFT_local_fields), intent(in) :: denspot
  integer, intent(in) :: iproc, nproc

  character(len = *), parameter :: subname = "denspot_emit_v_ext"
  integer, parameter :: SIGNAL_DONE = -1
  integer :: message, ierr, i_stat, i_all, new
  real(gp), pointer :: full_dummy(:)
  interface
     subroutine denspot_full_v_ext(denspot, pot_full, iproc, new)
       use module_types
       implicit none
       type(DFT_local_fields), intent(in) :: denspot
       integer, intent(in) :: iproc
       integer, intent(out) :: new
       real(gp), pointer :: pot_full(:)
     END SUBROUTINE denspot_full_v_ext
  end interface

  call timing(iproc,'rhov_signals  ','ON')
  if (iproc == 0) then
     ! Only iproc 0 emit the signal. This call is blocking.
     ! All other procs are blocked by the bcast to wait for
     ! possible transfer to proc 0.
     call localfields_emit_v_ext(denspot%c_obj)
     if (nproc > 1) then
        ! After handling the signal, iproc 0 broadcasts to other
        ! proc to continue (jproc == -1).
        message = SIGNAL_DONE
        call MPI_BCAST(message, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
     end if
  else
     do
        call MPI_BCAST(message, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
        if (message == SIGNAL_DONE) then
           exit
        else
           allocate(full_dummy(1+ndebug),stat=i_stat)
           call memocc(i_stat,full_dummy,'full_dummy',subname)
           ! Gather density to iproc 0
           call denspot_full_v_ext(denspot, full_dummy, iproc, new)
           i_all=-product(shape(full_dummy))*kind(full_dummy)
           deallocate(full_dummy,stat=i_stat)
           call memocc(i_stat,i_all,'full_dummy',subname)
        end if
     end do
  end if
  call timing(iproc,'rhov_signals  ','OF')
END SUBROUTINE denspot_emit_v_ext

subroutine allocateRhoPot(iproc,Glr,nspin,atoms,rxyz,denspot)
  use module_base
  use module_types
  use module_interfaces, except_this_one => allocateRhoPot
  use m_ab6_mixing
  implicit none
  integer, intent(in) :: iproc,nspin
  type(locreg_descriptors), intent(in) :: Glr
  type(atoms_data), intent(in) :: atoms
  real(gp), dimension(3,atoms%nat), intent(in) :: rxyz
  type(DFT_local_fields), intent(inout) :: denspot

  character(len = *), parameter :: subname = "allocateRhoPot"
  integer :: i_stat

  ! Allocate density and potentials.
  ! --------
  !allocate ionic potential
  if (denspot%dpbox%n3pi > 0) then
     allocate(denspot%V_ext(Glr%d%n1i,Glr%d%n2i,denspot%dpbox%n3pi,1+ndebug),stat=i_stat)
     call memocc(i_stat,denspot%V_ext,'V_ext',subname)
  else
     allocate(denspot%V_ext(1,1,1,1+ndebug),stat=i_stat)
     call memocc(i_stat,denspot%V_ext,'pot_ion',subname)
  end if
  !Allocate XC potential
  if (denspot%dpbox%n3p >0) then
     allocate(denspot%V_XC(Glr%d%n1i,Glr%d%n2i,denspot%dpbox%n3p,nspin+ndebug),stat=i_stat)
     call memocc(i_stat,denspot%V_XC,'V_XC',subname)
  else
     allocate(denspot%V_XC(1,1,1,nspin+ndebug),stat=i_stat)
     call memocc(i_stat,denspot%V_XC,'V_XC',subname)
  end if

  if (denspot%dpbox%n3d >0) then
     allocate(denspot%rhov(Glr%d%n1i*Glr%d%n2i*denspot%dpbox%n3d*&
          denspot%dpbox%nrhodim+ndebug),stat=i_stat)
     call memocc(i_stat,denspot%rhov,'rhov',subname)
  else
     allocate(denspot%rhov(denspot%dpbox%nrhodim+ndebug),stat=i_stat)
     call memocc(i_stat,denspot%rhov,'rhov',subname)
  end if
  !check if non-linear core correction should be applied, and allocate the 
  !pointer if it is the case
  !print *,'i3xcsh',denspot%dpbox%i3s,denspot%dpbox%i3xcsh,denspot%dpbox%n3d
  call calculate_rhocore(iproc,atoms,Glr%d,rxyz,&
       denspot%dpbox%hgrids(1),denspot%dpbox%hgrids(2),denspot%dpbox%hgrids(3),&
       denspot%dpbox%i3s,denspot%dpbox%i3xcsh,&
       denspot%dpbox%n3d,denspot%dpbox%n3p,denspot%rho_C)

!!$  !calculate the XC energy of rhocore
!!$  call xc_init_rho(denspot%dpbox%nrhodim,denspot%rhov,1)
!!$  call XC_potential(atoms%geocode,'D',iproc,nproc,&
!!$       Lzd%Glr%d%n1i,Lzd%Glr%d%n2i,Lzd%Glr%d%n3i,ixc,hxh,hyh,hzh,&
!!$       denspot%rhov,eexcu,vexcu,orbs%nspin,denspot%rho_C,denspot%V_XC,xcstr)


END SUBROUTINE allocateRhoPot

!!$!> Create the descriptors for the density and the potential
!!$subroutine createDensPotDescriptors(iproc,nproc,atoms,gdim,hxh,hyh,hzh,&
!!$     rxyz,crmult,frmult,radii_cf,nspin,datacode,ixc,rho_commun,&
!!$     n3d,n3p,n3pi,i3xcsh,i3s,nscatterarr,ngatherarr,rhodsc)
!> Create the descriptors for the density and the potential
subroutine dpbox_repartition(iproc,nproc,geocode,datacode,ixc,dpbox)

  use module_base
  use module_types
  use Poisson_Solver
  use module_xc
  implicit none
  !Arguments
  integer, intent(in) :: iproc,nproc,ixc
  character(len=1), intent(in) :: geocode,datacode
  type(denspot_distribution), intent(inout) :: dpbox
  !Local variables
  integer :: jproc,n3d,n3p,n3pi,i3xcsh,i3s

  if (datacode == 'D') then
     do jproc=0,nproc-1
        call PS_dim4allocation(geocode,datacode,jproc,nproc,&
             dpbox%ndims(1),dpbox%ndims(2),dpbox%ndims(3),ixc,&
             n3d,n3p,n3pi,i3xcsh,i3s)
        dpbox%nscatterarr(jproc,1)=n3d            !number of planes for the density
        dpbox%nscatterarr(jproc,2)=n3p            !number of planes for the potential
        dpbox%nscatterarr(jproc,3)=i3s+i3xcsh-1   !starting offset for the potential
        dpbox%nscatterarr(jproc,4)=i3xcsh         !GGA XC shift between density and potential
     end do
  end if

  if (iproc< nproc) then
     dpbox%n3d=dpbox%nscatterarr(iproc,1)
     dpbox%n3p=dpbox%nscatterarr(iproc,2)
     dpbox%i3xcsh=dpbox%nscatterarr(iproc,4)
     dpbox%i3s=dpbox%nscatterarr(iproc,3)-dpbox%i3xcsh+1
     dpbox%n3pi=dpbox%n3p
  else
     dpbox%n3d=0
     dpbox%n3p=0
     dpbox%i3xcsh=0
     dpbox%i3s=1
     dpbox%n3pi=dpbox%n3p
  end if

  dpbox%ngatherarr(:,1)=dpbox%ndims(1)*dpbox%ndims(2)*dpbox%nscatterarr(:,2)
  dpbox%ngatherarr(:,2)=dpbox%ndims(1)*dpbox%ndims(2)*dpbox%nscatterarr(:,3)
  !for the density
  dpbox%ngatherarr(:,3)=dpbox%ndims(1)*dpbox%ndims(2)*dpbox%nscatterarr(:,1)

end subroutine dpbox_repartition

!!$  !calculate dimensions of the complete array to be allocated before the reduction procedure
!!$  if (rhodsc%icomm==1) then
!!$     rhodsc%nrhotot=0
!!$     do jproc=0,nproc-1
!!$        rhodsc%nrhotot=rhodsc%nrhotot+nscatterarr(jproc,1)
!!$     end do
!!$  else
!!$     rhodsc%nrhotot=ndims(3)
!!$  end if

!END SUBROUTINE createDensPotDescriptors

subroutine density_descriptors(iproc,nproc,nspin,crmult,frmult,atoms,dpbox,&
     rho_commun,rxyz,radii_cf,rhodsc)
  use module_base
  use module_types
  use module_xc
  use module_interfaces
  implicit none
  integer, intent(in) :: iproc,nproc,nspin
  real(gp), intent(in) :: crmult,frmult
  type(atoms_data), intent(in) :: atoms
  type(denspot_distribution), intent(in) :: dpbox
  character(len=3), intent(in) :: rho_commun
  real(gp), dimension(3,atoms%nat), intent(in) :: rxyz
  real(gp), dimension(atoms%ntypes,3), intent(in) :: radii_cf
  type(rho_descriptors), intent(out) :: rhodsc
  !local variables
  integer :: ierr

  !decide rho communication strategy
  !old way
!!$  if (rho_commun=='MIX' .and. (atoms%geocode.eq.'F') .and. (nproc > 1)) then
!!$     rhodsc%icomm=2
!!$  else
!!$     if (.not.xc_isgga()) then
!!$        rhodsc%icomm=1
!!$     else
!!$        rhodsc%icomm=0
!!$     endif
!!$  end if
  !deafult
  if (.not.xc_isgga()) then
     rhodsc%icomm=1
  else
     rhodsc%icomm=0
  endif
  if ((atoms%geocode.eq.'F') .and. (nproc > 1)) then
     rhodsc%icomm=2
  end if
  !override the  default
  if (rho_commun=='DBL') then
     rhodsc%icomm=0
  else if (rho_commun == 'RSC') then
     rhodsc%icomm=1
  end if

  !in the case of taskgroups the RSC scheme should be overrided
  if (rhodsc%icomm==1 .and. size(dpbox%nscatterarr,1) < nproc) then
     if (atoms%geocode.eq.'F') then
        rhodsc%icomm=2
     else
        rhodsc%icomm=0
     end if
  end if
rhodsc%icomm=2
  !write (*,*) 'hxh,hyh,hzh',hgrids(1),hgrids(2),hgrids(3)
  !create rhopot descriptors
  !allocate rho_descriptors if the density repartition is activated
 
  if (rhodsc%icomm==2) then !rho_commun=='MIX' .and. (atoms%geocode.eq.'F') .and. (nproc > 1)) then! .and. xc_isgga()) then
     call rho_segkey(iproc,atoms,rxyz,crmult,frmult,radii_cf,&
          dpbox%ndims(1),dpbox%ndims(2),dpbox%ndims(3),&
          dpbox%hgrids(1),dpbox%hgrids(2),dpbox%hgrids(3),nspin,rhodsc,.false.)
  else
     !nullify rhodsc pointers
     nullify(rhodsc%spkey)
     nullify(rhodsc%dpkey)
     nullify(rhodsc%cseg_b)
     nullify(rhodsc%fseg_b)
  end if
  
  !calculate dimensions of the complete array to be allocated before the reduction procedure
  if (rhodsc%icomm==1) then
     rhodsc%nrhotot=sum(dpbox%nscatterarr(:,1))
  else
     rhodsc%nrhotot=dpbox%ndims(3)
  end if
 
end subroutine density_descriptors


!> routine which initialised the potential data
subroutine default_confinement_data(confdatarr,norbp)
  use module_base
  use module_types
  implicit none
  integer, intent(in) :: norbp
  type(confpot_data), dimension(norbp), intent(out) :: confdatarr
  !local variables
  integer :: iorb

  !initialize the confdatarr
  do iorb=1,norbp
     confdatarr(iorb)%potorder=0
     !the rest is not useful
     confdatarr(iorb)%prefac     =UNINITIALIZED(confdatarr(iorb)%prefac)     
     confdatarr(iorb)%hh(1)      =UNINITIALIZED(confdatarr(iorb)%hh(1))      
     confdatarr(iorb)%hh(2)      =UNINITIALIZED(confdatarr(iorb)%hh(2))      
     confdatarr(iorb)%hh(3)      =UNINITIALIZED(confdatarr(iorb)%hh(3))      
     confdatarr(iorb)%rxyzConf(1)=UNINITIALIZED(confdatarr(iorb)%rxyzConf(1))
     confdatarr(iorb)%rxyzConf(2)=UNINITIALIZED(confdatarr(iorb)%rxyzConf(2))
     confdatarr(iorb)%rxyzConf(3)=UNINITIALIZED(confdatarr(iorb)%rxyzConf(3))
     confdatarr(iorb)%ioffset(1) =UNINITIALIZED(confdatarr(iorb)%ioffset(1)) 
     confdatarr(iorb)%ioffset(2) =UNINITIALIZED(confdatarr(iorb)%ioffset(2)) 
     confdatarr(iorb)%ioffset(3) =UNINITIALIZED(confdatarr(iorb)%ioffset(3)) 

  end do
end subroutine default_confinement_data

subroutine define_confinement_data(confdatarr,orbs,rxyz,at,hx,hy,hz,&
           confpotorder,potentialprefac,Lzd,confinementCenter)
  use module_base
  use module_types
  implicit none
  real(gp), intent(in) :: hx,hy,hz
  type(atoms_data), intent(in) :: at
  type(orbitals_data), intent(in) :: orbs
  !!type(linearParameters), intent(in) :: lin
  integer,intent(in):: confpotorder
  real(gp),dimension(at%ntypes),intent(in):: potentialprefac
  type(local_zone_descriptors), intent(in) :: Lzd
  real(gp), dimension(3,at%nat), intent(in) :: rxyz
  integer, dimension(orbs%norb), intent(in) :: confinementCenter
  type(confpot_data), dimension(orbs%norbp), intent(out) :: confdatarr
  !local variables
  integer :: iorb,nl1,nl2,nl3,ilr,icenter

  !initialize the confdatarr
  do iorb=1,orbs%norbp
     ilr=orbs%inWhichlocreg(orbs%isorb+iorb)
     icenter=confinementCenter(orbs%isorb+iorb)
     !!confdatarr(iorb)%potorder=lin%confpotorder
     !!confdatarr(iorb)%prefac=lin%potentialprefac(at%iatype(icenter))
     confdatarr(iorb)%potorder=confpotorder
     confdatarr(iorb)%prefac=potentialprefac(at%iatype(icenter))
     confdatarr(iorb)%hh(1)=.5_gp*hx
     confdatarr(iorb)%hh(2)=.5_gp*hy
     confdatarr(iorb)%hh(3)=.5_gp*hz
     confdatarr(iorb)%rxyzConf(1:3)=rxyz(1:3,icenter)!Lzd%Llr(ilr)%locregCenter(1:3)
     call geocode_buffers(Lzd%Llr(ilr)%geocode,nl1,nl2,nl3)
     confdatarr(iorb)%ioffset(1)=lzd%llr(ilr)%nsi1-nl1-1
     confdatarr(iorb)%ioffset(2)=lzd%llr(ilr)%nsi2-nl2-1
     confdatarr(iorb)%ioffset(3)=lzd%llr(ilr)%nsi3-nl3-1
  end do

contains

    subroutine geocode_buffers(geocode,nl1,nl2,nl3)
      implicit none
      integer, intent(in) :: nl1,nl2,nl3
      character(len=1), intent(in) :: geocode
      !local variables
      logical :: perx,pery,perz
      integer :: nr1,nr2,nr3

      !conditions for periodicity in the three directions
      perx=(geocode /= 'F')
      pery=(geocode == 'P')
      perz=(geocode /= 'F')

      call ext_buffers(perx,nl1,nr1)
      call ext_buffers(pery,nl2,nr2)
      call ext_buffers(perz,nl3,nr3)

    end subroutine geocode_buffers
  
end subroutine define_confinement_data


!> Print the distribution schemes
subroutine print_distribution_schemes(unit,nproc,nkpts,norb_par,nvctr_par)
  use module_base
  implicit none
  !Arguments
  integer, intent(in) :: nproc,nkpts,unit
  integer, dimension(0:nproc-1,nkpts), intent(in) :: norb_par,nvctr_par
  !local variables
  integer :: jproc,ikpt,norbp,isorb,ieorb,isko,ieko,nvctrp,ispsi,iepsi,iekc,iskc
  integer :: iko,ikc,nko,nkc

  write(unit,'(1x,a,a)')repeat('-',46),'Direct and transposed data repartition'
  write(unit,'(1x,8(a))')'| proc |',' N. Orbitals | K-pt |  Orbitals  ',&
       '|| N. Components | K-pt |    Components   |'
  do jproc=0,nproc-1
     call start_end_distribution(nproc,nkpts,jproc,norb_par,isko,ieko,norbp)
     call start_end_distribution(nproc,nkpts,jproc,nvctr_par,iskc,iekc,nvctrp)
     iko=isko
     ikc=iskc
     nko=ieko-isko+1
     nkc=iekc-iskc+1
     !print total number of orbitals and components
     write(unit,'(1x,a,i4,a,i8,a,i13,a)')'| ',jproc,' |',norbp,&
          repeat(' ',5)//'|'//repeat('-',6)//'|'//repeat('-',12)//'||',&
          nvctrp,&
          repeat(' ',2)//'|'//repeat('-',6)//'|'//repeat('-',17)//'|'
     !change the values to zero if there is no orbital
     do ikpt=1,min(nko,nkc)
        call start_end_comps(nproc,jproc,norb_par(0,iko),isorb,ieorb)
        call start_end_comps(nproc,jproc,nvctr_par(0,ikc),ispsi,iepsi)
        if (norbp/=0) then
           write(unit,'(a,i4,a,i5,a,i5,a,i4,a,i8,a,i8,a)')&
                ' |'//repeat(' ',6)//'|'//repeat(' ',13)//'|',&
                iko,'  |',isorb,'-',ieorb,&
                ' ||'//repeat(' ',15)//'|',&
                ikc,'  |',ispsi,'-',iepsi,'|'
        else
           write(unit,'(a,i4,a,i5,a,i5,a,i4,a,i8,a,i8,a)')&
                ' |'//repeat(' ',6)//'|'//repeat(' ',13)//'|',&
                0,'  |',0,'-',-1,&
                ' ||'//repeat(' ',15)//'|',&
                ikc,'  |',ispsi,'-',iepsi,'|'
        end if
        iko=iko+1
        ikc=ikc+1
     end do
     if (nko > nkc) then
        do ikpt=nkc+1,nko
           if (norbp/=0) then
              call start_end_comps(nproc,jproc,norb_par(0,iko),isorb,ieorb)
              write(unit,'(a,i4,a,i5,a,i5,2a)') &
                   & ' |'//repeat(' ',6)//'|'//repeat(' ',13)//'|',&
                   & iko,'  |',isorb,'-',ieorb, ' ||'//repeat(' ',15)//'|',&
                   & '      |                 |'
           else
              write(unit,'(a,i4,a,i5,a,i5,2a)') &
                   & ' |'//repeat(' ',6)//'|'//repeat(' ',13)//'|',&
                   & 0,'  |',0,'-',-1, ' ||'//repeat(' ',15)//'|',&
                   & '      |                 |'
           end if
           iko=iko+1
        end do
     else if (nkc > nko) then
        do ikpt=nko+1,nkc
           call start_end_comps(nproc,jproc,nvctr_par(0,ikc),ispsi,iepsi)
           write(unit,'(a,i4,a,i8,a,i8,a)')&
                ' |'//repeat(' ',6)//'|'//repeat(' ',13)//'|'//repeat(' ',4)//'  |'//&
                repeat(' ',12)//'||'//repeat(' ',15)//'|',&
                ikc,'  |',ispsi,'-',iepsi,'|'
           ikc=ikc+1
        end do
     end if
  end do
  
END SUBROUTINE print_distribution_schemes


subroutine start_end_distribution(nproc,nkpts,jproc,ndist,is,ie,norbp)
  implicit none
  integer, intent(in) :: nproc,nkpts,jproc
  integer, dimension(0:nproc-1,nkpts), intent(in) :: ndist
  integer, intent(out) :: is,ie,norbp
  !local variables
  integer :: ikpt
  norbp=0
  do ikpt=1,nkpts
     norbp=norbp+ndist(jproc,ikpt)
  end do
  if (norbp == 0) then
     is=nkpts
     ie=nkpts
  end if
  loop_is: do ikpt=1,nkpts
     if (ndist(jproc,ikpt) /= 0) then
        is=ikpt
        exit loop_is
     end if
  end do loop_is
  loop_ie: do ikpt=nkpts,1,-1
     if (ndist(jproc,ikpt) /= 0) then
        ie=ikpt
        exit loop_ie
     end if
  end do loop_ie
END SUBROUTINE start_end_distribution


subroutine start_end_comps(nproc,jproc,ndist,is,ie)
  implicit none
  integer, intent(in) :: nproc,jproc
  integer, dimension(0:nproc-1), intent(in) :: ndist
  integer, intent(out) :: is,ie
  !local variables
  integer :: kproc

  is=1
  do kproc=0,jproc-1
     is=is+ndist(kproc)
  end do
  ie=is+ndist(jproc)-1
  
END SUBROUTINE start_end_comps
