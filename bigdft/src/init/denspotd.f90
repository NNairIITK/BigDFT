!> @file
!!  Routines to create descriptor arrays for density and potential
!! @author
!!    Copyright (C) 2007-2015 BigDFT group (LG)
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 


!> Denspot initialization
subroutine initialize_DFT_local_fields(denspot, ixc, nspden)
  use module_base
  use module_dpbox, only: dpbox_null
  use module_types
  use module_xc
  use public_enums
  implicit none
  type(DFT_local_fields), intent(inout) :: denspot
  integer, intent(in) :: ixc, nspden

  denspot%rhov_is = EMPTY
  nullify(denspot%rho_C)
  nullify(denspot%rhohat)
  nullify(denspot%V_ext)
  nullify(denspot%rho_ion)
  nullify(denspot%Vloc_KS)
  nullify(denspot%rho_psi)
  nullify(denspot%V_XC)
  nullify(denspot%f_XC)
  nullify(denspot%rho_work)
  nullify(denspot%pot_work)
  nullify(denspot%rhov)

  denspot%psoffset=0.0_gp

  if (verbose >1) then
     denspot%PSquiet='NO '
  else
     denspot%PSquiet='YES'
  end if

  call initialize_coulomb_operator(denspot%pkernel)
  call initialize_coulomb_operator(denspot%pkernelseq)
  call initialize_rho_descriptors(denspot%rhod)
  denspot%dpbox=dpbox_null()

  nullify(denspot%mix)

  if (ixc < 0) then
     call xc_init(denspot%xc, ixc, XC_MIXED, nspden)
  else
     call xc_init(denspot%xc, ixc, XC_ABINIT, nspden)
  end if
end subroutine initialize_DFT_local_fields


subroutine initialize_coulomb_operator(kernel)
  use module_base
  use module_types
  implicit none
  type(coulomb_operator), intent(out) :: kernel
  
  nullify(kernel%kernel)

  
end subroutine initialize_coulomb_operator


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


!> Initialize dpbox from the local zone descriptors
subroutine dpbox_set(dpbox,Lzd,xc,iproc,nproc,mpi_comm,PS_groupsize,SICapproach,geocode,nspin,igpu)
  use module_base
  use module_dpbox, only: denspot_distribution,dpbox_null
  use module_types
  use module_xc
  implicit none
  integer, intent(in) :: iproc,nproc,mpi_comm,PS_groupsize,nspin,igpu
  character(len=1), intent(in) :: geocode !< @copydoc poisson_solver::doc::geocode
  character(len=4), intent(in) :: SICapproach
  type(local_zone_descriptors), intent(in) :: Lzd
  type(xc_info), intent(in) :: xc
  type(denspot_distribution), intent(out) :: dpbox
  !local variables
  integer :: npsolver_groupsize

  dpbox=dpbox_null()

  call dpbox_set_box(dpbox,Lzd)

  !if the taskgroup size is not a divisor of nproc do not create taskgroups
  if (nproc > 1 .and. PS_groupsize > 0 .and. &
       PS_groupsize < nproc .and.&
       mod(nproc,PS_groupsize)==0) then
     npsolver_groupsize=PS_groupsize
  else
     npsolver_groupsize=nproc
  end if
  call mpi_environment_set(dpbox%mpi_env,iproc,nproc,mpi_comm,npsolver_groupsize)

  call denspot_communications(dpbox%mpi_env%iproc,dpbox%mpi_env%nproc,igpu,xc,&
                              nspin,geocode,SICapproach,dpbox)

end subroutine dpbox_set


!> Initialize dpbox (density pot distribution) i.e. the parameters defining the grid
subroutine dpbox_set_box(dpbox,Lzd)
  use module_base
  use module_dpbox, only: denspot_distribution
  use module_types
  implicit none
  type(local_zone_descriptors), intent(in) :: Lzd
  type(denspot_distribution), intent(inout) :: dpbox
 
  !The grid for the potential is twice finer
  dpbox%hgrids(1)=0.5_gp*Lzd%hgrids(1)
  dpbox%hgrids(2)=0.5_gp*Lzd%hgrids(2)
  dpbox%hgrids(3)=0.5_gp*Lzd%hgrids(3)
  !Same dimension
  dpbox%ndims(1)=Lzd%Glr%d%n1i
  dpbox%ndims(2)=Lzd%Glr%d%n2i
  dpbox%ndims(3)=Lzd%Glr%d%n3i
  dpbox%geocode=Lzd%Glr%geocode

end subroutine dpbox_set_box


subroutine denspot_set_history(denspot, iscf, &
     npulayit)
  use module_base
  use module_types
  use module_mixing
  use public_enums, only: SCF_KIND_DIRECT_MINIMIZATION
  implicit none
  type(DFT_local_fields), intent(inout) :: denspot
  integer, intent(in) :: iscf
  integer,intent(in),optional :: npulayit
  
  integer :: potden, npoints, ierr
  character(len=500) :: errmess

  if (iscf < 10) then
     ! Mixing over potential so use dimension of pot (n3p)
     potden = AB7_MIXING_POTENTIAL
     npoints = denspot%dpbox%ndims(1)*denspot%dpbox%ndims(2)*denspot%dpbox%n3p
!!!     npoints = n1i*n2i*denspot%dpbox%n3p
  else
     ! Mixing over density so use dimension of density (n3d)
     potden = AB7_MIXING_DENSITY
     npoints = denspot%dpbox%ndims(1)*denspot%dpbox%ndims(2)*denspot%dpbox%n3d
!!!     npoints = n1i*n2i*denspot%dpbox%n3d
  end if
  if (iscf > SCF_KIND_DIRECT_MINIMIZATION) then
     allocate(denspot%mix)
     if (present(npulayit)) then
         call ab7_mixing_new(denspot%mix, modulo(iscf, 10), potden, &
              AB7_MIXING_REAL_SPACE, npoints, denspot%dpbox%nrhodim, 0, &
              ierr, errmess, npulayit=npulayit, useprec = .false.)
     else
         call ab7_mixing_new(denspot%mix, modulo(iscf, 10), potden, &
              AB7_MIXING_REAL_SPACE, npoints, denspot%dpbox%nrhodim, 0, &
              ierr, errmess, useprec = .false.)
     end if
     call ab7_mixing_eval_allocate(denspot%mix)
  else
     nullify(denspot%mix)
  end if
end subroutine denspot_set_history


subroutine denspot_free_history(denspot)
  use module_types
  use module_mixing
  implicit none
  type(DFT_local_fields), intent(inout) :: denspot
  
  if (associated(denspot%mix)) then
      call ab7_mixing_deallocate(denspot%mix)
      deallocate(denspot%mix)
  end if
end subroutine denspot_free_history


!> Create descriptors for density and potentials (parallel distribution)
subroutine denspot_communications(iproc,nproc,igpu,xc,nspin,geocode,SICapproach,dpbox)
  use module_base
  use module_dpbox, only: denspot_distribution
  use module_types
  use module_xc
  implicit none
  integer, intent(in) :: nspin,iproc,nproc,igpu
  type(xc_info), intent(in) :: xc
  character(len=1), intent(in) :: geocode !< @copydoc poisson_solver::doc::geocode
  character(len=4), intent(in) :: SICapproach
  type(denspot_distribution), intent(inout) :: dpbox

  ! Create descriptors for density and potentials.

  ! these arrays should be included in the comms descriptor
  ! allocate values of the array for the data scattering in sumrho
  ! its values are ignored in the datacode='G' case
  dpbox%nscatterarr = f_malloc_ptr((/ 0.to.nproc-1, 1.to.4 /),id='dpbox%nscatterarr')
  !allocate array for the communications of the potential
  !also used for the density
  dpbox%ngatherarr = f_malloc_ptr((/ 0.to.nproc-1, 1.to.3 /),id='dpbox%ngatherarr')

  call dpbox_repartition(iproc,nproc,igpu,geocode,'D',xc,dpbox)

  ! Allocate Charge density / Potential in real space
  ! here the full_density treatment should be put
  dpbox%nrhodim=nspin
  dpbox%i3rho_add=0
  if (trim(SICapproach)=='NK') then
     dpbox%nrhodim=2*dpbox%nrhodim !to be eliminated with an orbital-dependent potential
     dpbox%i3rho_add=dpbox%ndims(1)*dpbox%ndims(2)*dpbox%i3xcsh+1
  end if

  !fill the full_local_potential dimension
  dpbox%ndimpot=dpbox%ndims(1)*dpbox%ndims(2)*dpbox%n3p
  dpbox%ndimgrid=dpbox%ndims(1)*dpbox%ndims(2)*dpbox%ndims(3)
  dpbox%ndimrhopot=dpbox%ndims(1)*dpbox%ndims(2)*dpbox%n3d*dpbox%nrhodim
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
  use memory_profiling
  implicit none
  type(DFT_local_fields), intent(in) :: denspot
  integer, intent(in) :: iproc
  integer, intent(out) :: new
  real(gp), dimension(:), pointer :: rho_full

  character(len = *), parameter :: subname = "denspot_full_density"
  integer :: nslice, ierr, irhodim, irhoxcsh

  new = 0
  nslice = max(denspot%dpbox%ndimpot, 1)
  if (nslice < denspot%dpbox%ndimgrid) then
     if (iproc == 0) then
        !allocate full density in pot_ion array
        rho_full = f_malloc_ptr(denspot%dpbox%ndimgrid*denspot%dpbox%nrhodim,id='rho_full')
        new = 1
        
        ! Ask to gather density to other procs.
        !LG: wtf is that? call MPI_BCAST(0, 1, MPI_INTEGER, 0, bigdft_mpi%mpi_comm, ierr)
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
                mpidtypd,0,bigdft_mpi%mpi_comm,ierr)
        else
           call MPI_GATHERV(denspot%rhov(nslice * (irhodim - 1) + irhoxcsh + 1),&
                nslice,mpidtypd,rho_full(1),&
                denspot%dpbox%ngatherarr(0,1),denspot%dpbox%ngatherarr(0,2),&
                mpidtypd,0,bigdft_mpi%mpi_comm,ierr)
        end if
     end do
  else
     rho_full => denspot%rhov
  end if
END SUBROUTINE denspot_full_density


subroutine denspot_full_v_ext(denspot, pot_full, iproc, new)
  use module_base
  use module_types
  use memory_profiling
  implicit none
  type(DFT_local_fields), intent(in) :: denspot
  integer, intent(in) :: iproc
  integer, intent(out) :: new
  real(gp), pointer :: pot_full(:)

  character(len = *), parameter :: subname = "localfields_full_potential"
  integer :: ierr

  new = 0
  if (denspot%dpbox%ndimpot < denspot%dpbox%ndimgrid) then
     if (iproc == 0) then
        !allocate full density in pot_ion array
        pot_full = f_malloc_ptr(denspot%dpbox%ndimgrid,id='pot_full')
        new = 1
      
        ! Ask to gather density to other procs.
        !!!call MPI_BCAST(1, 1, MPI_INTEGER, 0, bigdft_mpi%mpi_comm, ierr)
     end if

     call MPI_GATHERV(denspot%v_ext(1,1,1,1),max(denspot%dpbox%ndimpot, 1),&
          mpidtypd,pot_full(1),denspot%dpbox%ngatherarr(0,1),&
          denspot%dpbox%ngatherarr(0,2),mpidtypd,0,bigdft_mpi%mpi_comm,ierr)
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
  integer :: message, ierr, new
  real(gp), pointer :: full_dummy(:)
  interface
     subroutine denspot_full_density(denspot, rho_full, iproc, new)
       use module_defs, only: gp
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
        call mpibcast(message, 1,comm=bigdft_mpi%mpi_comm)
     end if
  else
     do
        call mpibcast(message, 1,comm=bigdft_mpi%mpi_comm)
        if (message == SIGNAL_DONE) then
           exit
        else if (message == SIGNAL_DENSITY) then
           full_dummy = f_malloc_ptr(denspot%dpbox%nrhodim,id='full_dummy')
           ! Gather density to iproc 0
           call denspot_full_density(denspot, full_dummy, iproc, new)
           call f_free_ptr(full_dummy)
        end if
     end do
  end if
  call timing(iproc,'rhov_signals  ','OF')
END SUBROUTINE denspot_emit_rhov


subroutine denspot_emit_v_ext(denspot, iproc, nproc)
  use module_base
  use module_types
  implicit none
  !Arguments
  type(DFT_local_fields), intent(in) :: denspot
  integer, intent(in) :: iproc, nproc
  !Local variables
  character(len = *), parameter :: subname = "denspot_emit_v_ext"
  integer, parameter :: SIGNAL_DONE = -1
  integer :: message, ierr, new
  real(gp), pointer :: full_dummy(:)
  interface
     subroutine denspot_full_v_ext(denspot, pot_full, iproc, new)
       use module_defs, only: gp
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
        call mpibcast(message, 1,comm=bigdft_mpi%mpi_comm)
     end if
  else
     do
        call mpibcast(message, 1,comm=bigdft_mpi%mpi_comm)
        !call MPI_BCAST(message, 1, MPI_INTEGER, 0, bigdft_mpi%mpi_comm, ierr)
        if (message == SIGNAL_DONE) then
           exit
        else
           full_dummy = f_malloc_ptr(1,id='full_dummy')
           ! Gather density to iproc 0
           call denspot_full_v_ext(denspot, full_dummy, iproc, new)
           call f_free_ptr(full_dummy)
        end if
     end do
  end if
  call timing(iproc,'rhov_signals  ','OF')
END SUBROUTINE denspot_emit_v_ext


!> Allocate density and potentials.
subroutine allocateRhoPot(Glr,nspin,atoms,rxyz,denspot)
  use module_base
  use module_types
  use module_interfaces, only: calculate_rhocore
  implicit none
  integer, intent(in) :: nspin
  type(locreg_descriptors), intent(in) :: Glr
  type(atoms_data), intent(in) :: atoms
  real(gp), dimension(3,atoms%astruct%nat), intent(in) :: rxyz
  type(DFT_local_fields), intent(inout) :: denspot

  !allocate ionic potential
  if (denspot%dpbox%n3pi > 0) then
     denspot%V_ext = f_malloc_ptr((/ Glr%d%n1i , Glr%d%n2i , denspot%dpbox%n3pi , 1 /),id='denspot%V_ext')
  else
     denspot%V_ext = f_malloc_ptr((/ 1 , 1 , 1 , 1 /),id='denspot%V_ext')
  end if
  !Allocate XC potential
  if (denspot%dpbox%n3p >0) then
     denspot%V_XC = f_malloc_ptr((/ Glr%d%n1i , Glr%d%n2i , denspot%dpbox%n3p , nspin /),id='denspot%V_XC')
  else
     denspot%V_XC = f_malloc_ptr((/ 1 , 1 , 1 , nspin /),id='denspot%V_XC')
  end if

  !allocate ionic density in the case of a cavity calculation
  if (denspot%pkernel%method /= 'VAC') then
     if (denspot%dpbox%n3pi > 0) then
        denspot%rho_ion = f_malloc_ptr([ Glr%d%n1i , Glr%d%n2i , denspot%dpbox%n3pi , 1 ],id='denspot%rho_ion')
     else
        denspot%rho_ion = f_malloc_ptr([ 1 , 1 , 1 , 1 ],id='denspot%rho_ion')
     end if
  else
     denspot%rho_ion = f_malloc_ptr([ 1 , 1 , 1 , 1 ],id='denspot%rho_ion')
  end if

  if (denspot%dpbox%n3d >0) then
     denspot%rhov = f_malloc_ptr(Glr%d%n1i*Glr%d%n2i*denspot%dpbox%n3d*&
          denspot%dpbox%nrhodim,id='denspot%rhov')
  else
     denspot%rhov = f_malloc0_ptr(denspot%dpbox%nrhodim,id='denspot%rhov')
  end if
  !check if non-linear core correction should be applied, and allocate the 
  !pointer if it is the case
  !print *,'i3xcsh',denspot%dpbox%i3s,denspot%dpbox%i3xcsh,denspot%dpbox%n3d
  call calculate_rhocore(atoms,rxyz,denspot%dpbox,denspot%rho_C)

END SUBROUTINE allocateRhoPot


!!$!> Create the descriptors for the density and the potential
!!$subroutine createDensPotDescriptors(iproc,nproc,atoms,gdim,hxh,hyh,hzh,&
!!$     rxyz,crmult,frmult,radii_cf,nspin,datacode,ixc,rho_commun,&
!!$     n3d,n3p,n3pi,i3xcsh,i3s,nscatterarr,ngatherarr,rhodsc)


!> Do the parallel distribution and the descriptors for the density and the potential
subroutine dpbox_repartition(iproc,nproc,igpu,geocode,datacode,xc,dpbox)

  use module_base
  use module_dpbox, only: denspot_distribution
  use module_types
  use Poisson_Solver
  use module_xc
  implicit none
  !Arguments
  integer, intent(in) :: iproc,nproc,igpu
  type(xc_info), intent(in) :: xc
  character(len=1), intent(in) :: geocode  !< @copydoc poisson_solver::doc::geocode
  character(len=1), intent(in) :: datacode !< @copydoc poisson_solver::doc::datacode
  type(denspot_distribution), intent(inout) :: dpbox
  !Local variables
  integer :: jproc,n3d,n3p,n3pi,i3xcsh,i3s

  if (datacode == 'D') then
     do jproc=0,nproc-1
        call PS_dim4allocation(geocode,datacode,jproc,nproc,&
             dpbox%ndims(1),dpbox%ndims(2),dpbox%ndims(3),xc_isgga(xc),(xc%ixc/=13),&
             igpu,n3d,n3p,n3pi,i3xcsh,i3s)
        dpbox%nscatterarr(jproc,1)=n3d            !number of planes for the density
        dpbox%nscatterarr(jproc,2)=n3p            !number of planes for the potential
        dpbox%nscatterarr(jproc,3)=i3s+i3xcsh-1   !starting offset for the potential
        dpbox%nscatterarr(jproc,4)=i3xcsh         !GGA XC shift between density and potential
     end do
  end if

  if (iproc < nproc) then
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


subroutine density_descriptors(iproc,nproc,xc,nspin,crmult,frmult,atoms,dpbox,&
     rho_commun,rxyz,rhodsc)
  use module_base
  use module_dpbox, only:  denspot_distribution
  use module_types
  use module_xc
  implicit none
  integer, intent(in) :: iproc,nproc,nspin
  type(xc_info), intent(in) :: xc
  real(gp), intent(in) :: crmult,frmult
  type(atoms_data), intent(in) :: atoms
  type(denspot_distribution), intent(in) :: dpbox
  character(len=3), intent(in) :: rho_commun
  real(gp), dimension(3,atoms%astruct%nat), intent(in) :: rxyz
  !real(gp), dimension(atoms%astruct%ntypes,3), intent(in) :: radii_cf
  type(rho_descriptors), intent(out) :: rhodsc
  !local variables

  if (.not.xc_isgga(xc)) then
     rhodsc%icomm=1
  else
     rhodsc%icomm=0
  endif

  !decide rho communication strategy
  !old way
  !override the  default
  if (rho_commun=='DBL') then
     rhodsc%icomm=0
  else if (rho_commun == 'RSC') then
     rhodsc%icomm=1
  else if (rho_commun=='MIX' .and. (atoms%astruct%geocode.eq.'F') .and. (nproc > 1)) then
     rhodsc%icomm=2
  end if
  
!!$  !recent way
!!$  if ((atoms%astruct%geocode.eq.'F') .and. (nproc > 1)) then
!!$     rhodsc%icomm=2
!!$  end if
!!$  !override the  default
!!$  if (rho_commun=='DBL') then
!!$     rhodsc%icomm=0
!!$  else if (rho_commun == 'RSC') then
!!$     rhodsc%icomm=1
!!$  end if

  !in the case of taskgroups the RSC scheme should be overridden
  if (rhodsc%icomm==1 .and. size(dpbox%nscatterarr,1) < nproc) then
     if (atoms%astruct%geocode.eq.'F') then
        rhodsc%icomm=2
     else
        rhodsc%icomm=0
     end if
  end if
  !write (*,*) 'hxh,hyh,hzh',hgrids(1),hgrids(2),hgrids(3)
  !create rhopot descriptors
  !allocate rho_descriptors if the density repartition is activated

  if (rhodsc%icomm==2) then !rho_commun=='MIX' .and. (atoms%astruct%geocode.eq.'F') .and. (nproc > 1)) then! .and. xc_isgga()) then
     call rho_segkey(iproc,atoms,rxyz,crmult,frmult,&
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
  use locreg_operations, only: confpot_data,nullify_confpot_data
  implicit none
  integer, intent(in) :: norbp
  type(confpot_data), dimension(norbp), intent(out) :: confdatarr
  !local variables
  integer :: iorb

  !initialize the confdatarr
  do iorb=1,norbp
     call nullify_confpot_data(confdatarr(iorb))
  end do
end subroutine default_confinement_data




subroutine define_confinement_data(confdatarr,orbs,rxyz,at,hx,hy,hz,&
           confpotorder,potentialprefac,Lzd,confinementCenter)
  use module_base
  use module_types
  use locreg_operations, only: confpot_data
  use bounds, only: geocode_buffers
  implicit none
  real(gp), intent(in) :: hx,hy,hz
  type(atoms_data), intent(in) :: at
  type(orbitals_data), intent(in) :: orbs
  !!type(linearParameters), intent(in) :: lin
  integer,intent(in):: confpotorder
  real(gp),dimension(at%astruct%ntypes),intent(in):: potentialprefac
  type(local_zone_descriptors), intent(in) :: Lzd
  real(gp), dimension(3,at%astruct%nat), intent(in) :: rxyz
  integer, dimension(orbs%norb), intent(in) :: confinementCenter
  type(confpot_data), dimension(orbs%norbp), intent(out) :: confdatarr
  !local variables
  integer :: iorb,nl1,nl2,nl3,ilr,icenter

  !initialize the confdatarr
  do iorb=1,orbs%norbp
     ilr=orbs%inWhichlocreg(orbs%isorb+iorb)
     icenter=confinementCenter(orbs%isorb+iorb)
     !!confdatarr(iorb)%potorder=lin%confpotorder
     !!confdatarr(iorb)%prefac=lin%potentialprefac(at%astruct%iatype(icenter))
     confdatarr(iorb)%potorder=confpotorder
     confdatarr(iorb)%prefac=potentialprefac(at%astruct%iatype(icenter))
     confdatarr(iorb)%hh(1)=.5_gp*hx
     confdatarr(iorb)%hh(2)=.5_gp*hy
     confdatarr(iorb)%hh(3)=.5_gp*hz
     confdatarr(iorb)%rxyzConf(1:3)=rxyz(1:3,icenter)!Lzd%Llr(ilr)%locregCenter(1:3)
     call geocode_buffers(Lzd%Llr(ilr)%geocode, lzd%glr%geocode, nl1, nl2, nl3)
     confdatarr(iorb)%ioffset(1)=lzd%llr(ilr)%nsi1-nl1-1
     confdatarr(iorb)%ioffset(2)=lzd%llr(ilr)%nsi2-nl2-1
     confdatarr(iorb)%ioffset(3)=lzd%llr(ilr)%nsi3-nl3-1
     !confdatarr(iorb)%ioffset(1)=lzd%llr(ilr)%nsi1-1
     !confdatarr(iorb)%ioffset(2)=lzd%llr(ilr)%nsi2-1
     !confdatarr(iorb)%ioffset(3)=lzd%llr(ilr)%nsi3-1
     !confdatarr(iorb)%ioffset(1)=modulo(lzd%llr(ilr)%nsi1-1,lzd%glr%d%n1i)+1-nl1-1
     !confdatarr(iorb)%ioffset(2)=modulo(lzd%llr(ilr)%nsi2-1,lzd%glr%d%n2i)+1-nl2-1
     !confdatarr(iorb)%ioffset(3)=modulo(lzd%llr(ilr)%nsi3-1,lzd%glr%d%n3i)+1-nl3-1
     confdatarr(iorb)%damping   =1.0_gp
  end do

end subroutine define_confinement_data

  


!> Print the distribution schemes
subroutine print_distribution_schemes(nproc,nkpts,norb_par,nvctr_par)
  use module_base
  use yaml_output
  implicit none
  !Arguments
  integer, intent(in) :: nproc,nkpts
  integer, dimension(0:nproc-1,nkpts), intent(in) :: norb_par,nvctr_par
  !local variables
  integer :: jproc,ikpt,norbp,isorb,ieorb,isko,ieko,nvctrp,ispsi,iepsi,iekc,iskc
  integer :: iko,ikc,nko,nkc
  integer :: indentlevel

  call yaml_sequence_open('Direct and transposed data repartition')
     do jproc=0,nproc-1
        call start_end_distribution(nproc,nkpts,jproc,norb_par,isko,ieko,norbp)
        call start_end_distribution(nproc,nkpts,jproc,nvctr_par,iskc,iekc,nvctrp)
        iko=isko
        ikc=iskc
        nko=ieko-isko+1
        nkc=iekc-iskc+1
        !print total number of orbitals and components
        call yaml_mapping_open('Process'//trim(yaml_toa(jproc)))

           call yaml_map('Orbitals and Components', (/ norbp, nvctrp /))
           if (norbp /= 0) then
              call yaml_stream_attributes(indent=indentlevel)
              call yaml_sequence_open('Distribution',flow=.true.)
              call yaml_comment('Orbitals: [From, To], Components: [From, To]')
                 call yaml_newline()
                 do ikpt=1,min(nko,nkc)
                    call start_end_comps(nproc,jproc,norb_par(0,iko),isorb,ieorb)
                    call start_end_comps(nproc,jproc,nvctr_par(0,ikc),ispsi,iepsi)
                    call yaml_newline()
                    call yaml_sequence_open(repeat(' ', max(indentlevel+1,0)) // &
                         & "Kpt"//trim(yaml_toa(iko,fmt='(i4.4)')),flow=.true.)
                       call yaml_map("Orbitals",(/ isorb, ieorb /),fmt='(i5)')
                       call yaml_map("Components",(/ ispsi, iepsi /),fmt='(i8)')
                    call yaml_sequence_close()
                    iko=iko+1
                    ikc=ikc+1
                 end do
                 if (nko > nkc) then
                    do ikpt=nkc+1,nko
                       call start_end_comps(nproc,jproc,norb_par(0,iko),isorb,ieorb)
                       call yaml_sequence_open("Kpt"//trim(yaml_toa(iko,fmt='(i4.4)')),flow=.true.)
                       call yaml_map("Orbitals",(/ isorb, ieorb /),fmt='(i5)')
                       call yaml_sequence_close()
                       call yaml_newline()
                       iko=iko+1
                    end do
                 else if (nkc > nko) then
                    do ikpt=nko+1,nkc
                       call start_end_comps(nproc,jproc,nvctr_par(0,ikc),ispsi,iepsi)
                       call yaml_sequence_open("Kpt"//trim(yaml_toa(iko,fmt='(i4.4)')),flow=.true.)
                       call yaml_map("Components",(/ ispsi, iepsi /),fmt='(i8)')
                       call yaml_sequence_close()
                       call yaml_newline()
                    end do
                 end if
              call yaml_sequence_close()
           end if

        call yaml_mapping_close() ! for Process jproc
     end do
  call yaml_sequence_close()  ! for Data distribution
  
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
