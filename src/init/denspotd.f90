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
  type(DFT_local_fields), intent(out) :: denspot
  !local variables
  integer :: i

  nullify(denspot%rho_C,denspot%V_ext,denspot%Vloc_KS,denspot%rho_psi)
  nullify(denspot%V_XC,denspot%pkernel,denspot%pkernelseq)
  nullify(denspot%f_XC,denspot%rho_full,denspot%pot_full,denspot%rhov)

  denspot%psoffset=0.0_gp

  do i=1,3
     denspot%hgrids(i)=uninitialized(denspot%hgrids(i))
  end do
  if (verbose >1) then
     denspot%PSquiet='NO '
  else
     denspot%PSquiet='YES'
  end if

  call initialize_rho_descriptors(denspot%rhod)
  call initialize_denspot_distribution(denspot%dpcom)

end subroutine initialize_DFT_local_fields

subroutine initialize_denspot_distribution(dpcom)
  use module_base
  use module_types
  implicit none
  type(denspot_distribution), intent(out) :: dpcom
  
  dpcom%n3d      =uninitialized(dpcom%n3d)      
  dpcom%n3p      =uninitialized(dpcom%n3p)      
  dpcom%n3pi     =uninitialized(dpcom%n3pi)     
  dpcom%i3xcsh   =uninitialized(dpcom%i3xcsh)   
  dpcom%i3s      =uninitialized(dpcom%i3s)      
  dpcom%nrhodim  =uninitialized(dpcom%nrhodim)  
  dpcom%i3rho_add=uninitialized(dpcom%i3rho_add)

  nullify(dpcom%nscatterarr,dpcom%ngatherarr)
  
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

subroutine denspot_communications(iproc,nproc,grid,hxh,hyh,hzh,in,atoms,rxyz,radii_cf,dpcom,rhod)
  use module_base
  use module_types
  use module_interfaces, except_this_one => denspot_communications
  implicit none
  integer, intent(in) :: iproc, nproc
  type(grid_dimensions), intent(in) :: grid
  real(gp), intent(in) :: hxh, hyh, hzh
  type(input_variables), intent(in) :: in
  type(atoms_data), intent(in) :: atoms
  real(gp), dimension(3,atoms%nat), intent(in) :: rxyz
  real(gp), dimension(atoms%ntypes,3), intent(in) :: radii_cf
  type(denspot_distribution), intent(inout) :: dpcom
  type(rho_descriptors), intent(out) :: rhod
  !local variables
  character(len = *), parameter :: subname = 'denspot_communications' 
  integer :: i_stat

  ! Create descriptors for density and potentials.
  ! ------------------
  !these arrays should be included in the comms descriptor
  !allocate values of the array for the data scattering in sumrho
  !its values are ignored in the datacode='G' case
  allocate(dpcom%nscatterarr(0:nproc-1,4+ndebug),stat=i_stat)
  call memocc(i_stat,dpcom%nscatterarr,'nscatterarr',subname)
  !allocate array for the communications of the potential
  allocate(dpcom%ngatherarr(0:nproc-1,2+ndebug),stat=i_stat)
  call memocc(i_stat,dpcom%ngatherarr,'ngatherarr',subname)

  !create the descriptors for the density and the potential
  !these descriptors should take into account the localisation regions
  call createDensPotDescriptors(iproc,nproc,atoms,grid,hxh,hyh,hzh, &
       rxyz,in%crmult,in%frmult,radii_cf,in%nspin,'D',in%ixc,in%rho_commun, &
       dpcom%n3d,dpcom%n3p,&
       dpcom%n3pi,dpcom%i3xcsh,dpcom%i3s, &
       dpcom%nscatterarr,dpcom%ngatherarr,rhod)

  !Allocate Charge density / Potential in real space
  !here the full_density treatment should be put
  dpcom%nrhodim=in%nspin
  dpcom%i3rho_add=0
  if (trim(in%SIC%approach)=='NK') then
     dpcom%nrhodim=2*dpcom%nrhodim
     dpcom%i3rho_add=grid%n1i*grid%n2i*dpcom%i3xcsh+1
  end if

  !fill the full_local_potential dimension
  dpcom%ndimpot=grid%n1i*grid%n2i*dpcom%n3p
  dpcom%ndimgrid=grid%n1i*grid%n2i*grid%n3i
  dpcom%ndimrhopot=grid%n1i*grid%n2i*dpcom%n3d*&
       dpcom%nrhodim

end subroutine denspot_communications

subroutine allocateRhoPot(iproc,nproc,Glr,hxh,hyh,hzh,in,atoms,rxyz,radii_cf,denspot)
  use module_base
  use module_types
  use module_interfaces, except_this_one => allocateRhoPot
  implicit none
  integer, intent(in) :: iproc, nproc
  type(locreg_descriptors), intent(in) :: Glr
  real(gp), intent(in) :: hxh, hyh, hzh
  type(input_variables), intent(in) :: in
  type(atoms_data), intent(in) :: atoms
  real(gp), dimension(3,atoms%nat), intent(in) :: rxyz
  real(gp), dimension(atoms%ntypes,3), intent(in) :: radii_cf
  type(DFT_local_fields), intent(inout) :: denspot

  character(len = *), parameter :: subname = "allocateRhoPot"
  integer :: i_stat

  ! Allocate density and potentials.
  ! --------
  !allocate ionic potential
  if (denspot%dpcom%n3pi > 0) then
     allocate(denspot%V_ext(Glr%d%n1i,Glr%d%n2i,denspot%dpcom%n3pi,1+ndebug),stat=i_stat)
     call memocc(i_stat,denspot%V_ext,'V_ext',subname)
  else
     allocate(denspot%V_ext(1,1,1,1+ndebug),stat=i_stat)
     call memocc(i_stat,denspot%V_ext,'pot_ion',subname)
  end if
  !Allocate XC potential
  if (denspot%dpcom%n3p >0) then
     allocate(denspot%V_XC(Glr%d%n1i,Glr%d%n2i,denspot%dpcom%n3p,in%nspin+ndebug),stat=i_stat)
     call memocc(i_stat,denspot%V_XC,'V_XC',subname)
  else
     allocate(denspot%V_XC(1,1,1,in%nspin+ndebug),stat=i_stat)
     call memocc(i_stat,denspot%V_XC,'V_XC',subname)
  end if

  if (denspot%dpcom%n3d >0) then
     allocate(denspot%rhov(Glr%d%n1i*Glr%d%n2i*denspot%dpcom%n3d*&
          denspot%dpcom%nrhodim+ndebug),stat=i_stat)
     call memocc(i_stat,denspot%rhov,'rhov',subname)
  else
     allocate(denspot%rhov(denspot%dpcom%nrhodim+ndebug),stat=i_stat)
     call memocc(i_stat,denspot%rhov,'rhov',subname)
  end if
  !check if non-linear core correction should be applied, and allocate the 
  !pointer if it is the case
  call calculate_rhocore(iproc,atoms,Glr%d,rxyz,hxh,hyh,hzh, &
       denspot%dpcom%i3s,denspot%dpcom%i3xcsh,&
       denspot%dpcom%n3d,denspot%dpcom%n3p,denspot%rho_C)
  
END SUBROUTINE allocateRhoPot

!> Create the descriptors for the density and the potential
subroutine createDensPotDescriptors(iproc,nproc,atoms,gdim,hxh,hyh,hzh,&
     rxyz,crmult,frmult,radii_cf,nspin,datacode,ixc,rho_commun,&
     n3d,n3p,n3pi,i3xcsh,i3s,nscatterarr,ngatherarr,rhodsc)
  use module_base
  use module_types
  use Poisson_Solver
  use module_xc
  implicit none
  !Arguments
  character(len=1), intent(in) :: datacode
  character(len=3), intent(in) :: rho_commun
  integer, intent(in) :: iproc,nproc,ixc,nspin
  real(gp), intent(in) :: crmult,frmult,hxh,hyh,hzh
  type(atoms_data), intent(in) :: atoms
  type(grid_dimensions), intent(in) :: gdim
  real(gp), dimension(atoms%ntypes,3), intent(in) :: radii_cf
  real(gp), dimension(3,atoms%nat), intent(in) :: rxyz
  integer, intent(out) ::  n3d,n3p,n3pi,i3xcsh,i3s
  type(rho_descriptors), intent(out) :: rhodsc
  integer, dimension(0:nproc-1,4), intent(out) :: nscatterarr
  integer, dimension(0:nproc-1,2), intent(out) :: ngatherarr
  !Local variables
  integer :: jproc

  if (datacode == 'D') then
     do jproc=0,iproc-1
        call PS_dim4allocation(atoms%geocode,datacode,jproc,nproc,&
             gdim%n1i,gdim%n2i,gdim%n3i,ixc,&
             n3d,n3p,n3pi,i3xcsh,i3s)
        nscatterarr(jproc,1)=n3d            !number of planes for the density
        nscatterarr(jproc,2)=n3p            !number of planes for the potential
        nscatterarr(jproc,3)=i3s+i3xcsh-1   !starting offset for the potential
        nscatterarr(jproc,4)=i3xcsh         !GGA XC shift between density and potential
     end do
     do jproc=iproc+1,nproc-1
        call PS_dim4allocation(atoms%geocode,datacode,jproc,nproc,&
             gdim%n1i,gdim%n2i,gdim%n3i,ixc,&
             n3d,n3p,n3pi,i3xcsh,i3s)
        nscatterarr(jproc,1)=n3d
        nscatterarr(jproc,2)=n3p
        nscatterarr(jproc,3)=i3s+i3xcsh-1
        nscatterarr(jproc,4)=i3xcsh
     end do
  end if

  call PS_dim4allocation(atoms%geocode,datacode,iproc,nproc,&
       gdim%n1i,gdim%n2i,gdim%n3i,ixc,&
       n3d,n3p,n3pi,i3xcsh,i3s)
  nscatterarr(iproc,1)=n3d
  nscatterarr(iproc,2)=n3p
  nscatterarr(iproc,3)=i3s+i3xcsh-1
  nscatterarr(iproc,4)=i3xcsh

  ngatherarr(:,1)=gdim%n1i*gdim%n2i*nscatterarr(:,2)
  ngatherarr(:,2)=gdim%n1i*gdim%n2i*nscatterarr(:,3)

!write (*,*) 'hxh,hyh,hzh',hxh,hyh,hzh
  !create rhopot descriptors
    !allocate rho_descriptors if the density repartition is activated
  !decide rho communication strategy
  if (rho_commun=='MIX' .and. (atoms%geocode.eq.'F') .and. (nproc > 1) .and. xc_isgga()) then
     call rho_segkey(iproc,atoms,rxyz,crmult,frmult,radii_cf,&
          gdim%n1i,gdim%n2i,gdim%n3i,&
          hxh,hyh,hzh,nspin,rhodsc,.false.)
     rhodsc%icomm=2
  else
     !nullify rhodsc pointers
     nullify(rhodsc%spkey)
     nullify(rhodsc%dpkey)
     nullify(rhodsc%cseg_b)
     nullify(rhodsc%fseg_b)
     if (.not.xc_isgga()) then
        rhodsc%icomm=1
     else
        rhodsc%icomm=0
     endif
  end if

  !calculate dimensions of the complete array to be allocated before the reduction procedure
  if (rhodsc%icomm==1) then
     rhodsc%nrhotot=0
     do jproc=0,nproc-1
        rhodsc%nrhotot=rhodsc%nrhotot+nscatterarr(jproc,1)
     end do
  else
     rhodsc%nrhotot=gdim%n3i
  end if

END SUBROUTINE createDensPotDescriptors

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
  integer :: iorb,nl1,nl2,nl3,icenter,ilr

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
     confdatarr(iorb)%rxyzConf(1:3)=rxyz(1:3,icenter)
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


!> Partition the orbitals between processors to ensure load balancing
!! the criterion will depend on GPU computation
!! and/or on the sizes of the different localisation region.
!!
!! Calculate the number of elements to be sent to each process
!! and the array of displacements.
!! Cubic strategy: 
!!    - the components are equally distributed among the wavefunctions
!!    - each processor has all the orbitals in transposed form
!!    - each wavefunction is equally distributed in its transposed form
!!    - this holds for each k-point, which regroups different processors
subroutine orbitals_communicators(iproc,nproc,lr,orbs,comms,basedist)
  use module_base
  use module_types
  implicit none
  integer, intent(in) :: iproc,nproc
  type(locreg_descriptors), intent(in) :: lr
  type(orbitals_data), intent(inout) :: orbs
  type(communications_arrays), intent(out) :: comms
  integer, dimension(0:nproc-1,orbs%nkpts), intent(in), optional :: basedist
  !local variables
  character(len=*), parameter :: subname='orbitals_communicators'
  logical :: yesorb,yescomp
  integer :: jproc,nvctr_tot,ikpts,iorbp,jorb,norb_tot,ikpt,i_stat,i_all
  integer :: nkptsp,ierr,kproc,jkpts,jkpte,jsorb,lubo,lubc,info,jkpt,nB,nKB,nMB
  integer, dimension(:), allocatable :: mykpts
  logical, dimension(:), allocatable :: GPU_for_comp
  integer, dimension(:,:), allocatable :: nvctr_par,norb_par !<for all the components and orbitals (with k-pts)
  
  !check of allocation of important arrays
  if (.not. associated(orbs%norb_par)) then
     write(*,*)'ERROR: norb_par array not allocated'
     stop
  end if

  allocate(nvctr_par(0:nproc-1,0:orbs%nkpts+ndebug),stat=i_stat)
  call memocc(i_stat,nvctr_par,'nvctr_par',subname)
  allocate(norb_par(0:nproc-1,0:orbs%nkpts+ndebug),stat=i_stat)
  call memocc(i_stat,norb_par,'norb_par',subname)
  allocate(mykpts(orbs%nkpts+ndebug),stat=i_stat)
  call memocc(i_stat,mykpts,'mykpts',subname)

  !initialise the arrays
  do ikpts=0,orbs%nkpts
     do jproc=0,nproc-1
        nvctr_par(jproc,ikpts)=0 
        norb_par(jproc,ikpts)=0 
     end do
  end do

  !calculate the same k-point distribution for the orbitals
  !assign the k-point to the given orbital, counting one orbital after each other
  jorb=1
  ikpts=1
  do jproc=0,nproc-1
     do iorbp=1,orbs%norb_par(jproc,0)
        norb_par(jproc,ikpts)=norb_par(jproc,ikpts)+1
        if (mod(jorb,orbs%norb)==0) then
           ikpts=ikpts+1
        end if
        jorb=jorb+1
     end do
  end do
  !some checks
  if (orbs%norb /= 0) then
     !check the distribution
     do ikpts=1,orbs%nkpts
        !print *,'partition',ikpts,orbs%nkpts,'ikpts',norb_par(:,ikpts)
        norb_tot=0
        do jproc=0,nproc-1
           norb_tot=norb_tot+norb_par(jproc,ikpts)
        end do
        if(norb_tot /= orbs%norb) then
           write(*,*)'ERROR: partition of orbitals incorrect, kpoint:',ikpts
           stop
        end if
     end do
  end if


  !balance the components between processors
  !in the most symmetric way
  !here the components are taken into account for all the k-points

  !create an array which indicate which processor has a GPU associated 
  !from the viewpoint of the BLAS routines (deprecated, not used anymore)
  allocate(GPU_for_comp(0:nproc-1+ndebug),stat=i_stat)
  call memocc(i_stat,GPU_for_comp,'GPU_for_comp',subname)

  if (nproc > 1 .and. .not. GPUshare) then
     call MPI_ALLGATHER(GPUblas,1,MPI_LOGICAL,GPU_for_comp(0),1,MPI_LOGICAL,&
          MPI_COMM_WORLD,ierr)
  else
     GPU_for_comp(0)=GPUblas
  end if

  i_all=-product(shape(GPU_for_comp))*kind(GPU_for_comp)
  deallocate(GPU_for_comp,stat=i_stat)
  call memocc(i_stat,i_all,'GPU_for_comp',subname)

  !old k-point repartition
!!$  !decide the repartition for the components in the same way as the orbitals
!!$  call parallel_repartition_with_kpoints(nproc,orbs%nkpts,(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f),nvctr_par)

!!$  ikpts=1
!!$  ncomp_res=(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f)
!!$  do jproc=0,nproc-1
!!$     loop_comps: do
!!$        if (nvctr_par(jproc,0) >= ncomp_res) then
!!$           nvctr_par(jproc,ikpts)= ncomp_res
!!$           ikpts=ikpts+1
!!$           nvctr_par(jproc,0)=nvctr_par(jproc,0)-ncomp_res
!!$           ncomp_res=(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f)
!!$        else
!!$           nvctr_par(jproc,ikpts)= nvctr_par(jproc,0)
!!$           ncomp_res=ncomp_res-nvctr_par(jproc,0)
!!$           nvctr_par(jproc,0)=0
!!$           exit loop_comps
!!$        end if
!!$        if (nvctr_par(jproc,0) == 0 ) then
!!$           ncomp_res=(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f)
!!$           exit loop_comps
!!$        end if
!!$
!!$     end do loop_comps
!!$  end do

  !new k-point repartition
  if (present(basedist)) then
     do jkpt=1,orbs%nkpts
        do jproc=0,nproc-1
           nvctr_par(jproc,jkpt)=basedist(jproc,jkpt)
        end do
     end do
  else
     !first try the naive repartition
     call kpts_to_procs_via_obj(nproc,orbs%nkpts,(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f),nvctr_par(0,1))
  end if
  !then silently check whether the distribution agree
  info=-1
  call check_kpt_distributions(nproc,orbs%nkpts,orbs%norb,(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f),&
       norb_par(0,1),nvctr_par(0,1),info,lubo,lubc)
  if (info/=0 .and. .not. present(basedist)) then !redo the distribution based on the orbitals scheme
     info=-1
     call components_kpt_distribution(nproc,orbs%nkpts,orbs%norb,(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f),norb_par(0,1),nvctr_par(0,1))
     call check_kpt_distributions(nproc,orbs%nkpts,orbs%norb,(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f),&
          norb_par(0,1),nvctr_par(0,1),info,lubo,lubc)
  end if
  if (info /=0) then
     if (iproc==0) then
        write(*,*)'ERROR for nproc,nkpts,norb,nvctr',nproc,orbs%nkpts,orbs%norb,(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f)
        call print_distribution_schemes(6,nproc,orbs%nkpts,norb_par(0,1),nvctr_par(0,1))
     end if
     call MPI_BARRIER(MPI_COMM_WORLD,ierr)
     stop
  end if

!write(*,'(a,i2,3x,8i7,i10)') 'iproc, nvctr_par(jproc), sum', iproc, (nvctr_par(jproc,1), jproc=0,nproc-1), sum(nvctr_par(:,1))
!write(*,*) 'iproc, (lr%wfd%nvctr_c+7*lr%wfd%nvctr_f)*orbs%norbp', iproc, (lr%wfd%nvctr_c+7*lr%wfd%nvctr_f)*orbs%norbp
  !some checks
  !check the distribution
  do ikpts=1,orbs%nkpts
     !print *,'iproc,cpts:',lr%wfd%nvctr_c+7*lr%wfd%nvctr_f,nvctr_par(:,ikpts)
     nvctr_tot=0
     do jproc=0,nproc-1
        nvctr_tot=nvctr_tot+nvctr_par(jproc,ikpts)
     end do
     if(nvctr_tot /= lr%wfd%nvctr_c+7*lr%wfd%nvctr_f) then
        write(*,*)'ERROR: partition of components incorrect, kpoint:',ikpts
        stop
     end if
  end do

  !this function which associates a given k-point to a processor in the component distribution
  !the association is chosen such that each k-point is associated to only
  !one processor
  !if two processors treat the same k-point the processor which highest rank is chosen
  do ikpts=1,orbs%nkpts
     loop_jproc: do jproc=nproc-1,0,-1
        if (nvctr_par(jproc,ikpts) /= 0) then
           orbs%ikptproc(ikpts)=jproc
           exit loop_jproc
        end if
     end do loop_jproc
  end do
  
  !print*,'check',orbs%ikptproc(:)

!write(*,*) 'orbs%norb_par',orbs%norb_par

  !calculate the number of k-points treated by each processor in both
  ! the component distribution and the orbital distribution.
  !to have a correct distribution, a k-point should be divided between the same processors
  nkptsp=0
  orbs%iskpts=-1
  do ikpts=1,orbs%nkpts
     if (nvctr_par(iproc,ikpts) /= 0 .or. norb_par(iproc,ikpts) /= 0) then
        if (orbs%iskpts == -1) orbs%iskpts=ikpts-1
        nkptsp=nkptsp+1
        mykpts(nkptsp) = ikpts
     end if
  end do
  orbs%nkptsp=nkptsp

!!$  allocate(orbs%ikptsp(orbs%nkptsp+ndebug),stat=i_stat)
!!$  call memocc(i_stat,orbs%ikptsp,'orbs%ikptsp',subname)
!!$  orbs%ikptsp(1:orbs%nkptsp)=mykpts(1:orbs%nkptsp)

  !print the distribution scheme used for this set of orbital
  !in the case of multiple k-points
  if (iproc == 0 .and. verbose > 1 .and. orbs%nkpts > 1) then
     call print_distribution_schemes(6,nproc,orbs%nkpts,norb_par(0,1),nvctr_par(0,1))
  end if

  !print *,iproc,orbs%nkptsp,orbs%norbp,orbs%norb,orbs%nkpts
  !call MPI_BARRIER(MPI_COMM_WORLD,ierr)
  !call MPI_FINALIZE(ierr)
  !stop
  !check that for any processor the orbital k-point repartition is contained into the components
  do jproc=0,nproc-1
     jsorb=0
     do kproc=0,jproc-1
        jsorb=jsorb+orbs%norb_par(kproc,0)
     end do
     jkpts=min(jsorb/orbs%norb+1,orbs%nkpts)
     if (nvctr_par(jproc,jkpts) == 0 .and. orbs%norb_par(jproc,0) /=0 ) then
        if (iproc ==0) write(*,*)'ERROR, jproc: ',jproc,' the orbital k-points distribution starts before the components one'
        !print *,jsorb,jkpts,jproc,orbs%iskpts,nvctr_par(jproc,jkpts)
        stop
     end if
     jkpte=min((jsorb+orbs%norb_par(jproc,0)-1)/orbs%norb+1,orbs%nkpts)
     if (nvctr_par(jproc,jkpte) == 0 .and. orbs%norb_par(jproc,0) /=0) then
        if (iproc ==0) write(*,*)'ERROR, jproc: ',jproc,&
             ' the orbital k-points distribution ends after the components one'
        print *,jsorb,jkpte,jproc,orbs%iskpts,orbs%nkptsp,nvctr_par(jproc,jkpte)
        stop
     end if
  end do

  !before printing the distribution schemes, check that the two distributions contain
  !the same k-points
  yesorb=.false.
  kpt_components: do ikpts=1,orbs%nkptsp
     ikpt=orbs%iskpts+ikpts
     do jorb=1,orbs%norbp
        if (orbs%iokpt(jorb) == ikpt) yesorb=.true.
     end do
     if (.not. yesorb .and. orbs%norbp /= 0) then
        write(*,*)' ERROR: processor ', iproc,' kpt ',ikpt,&
             ' not found in the orbital distribution'
        call MPI_ABORT(MPI_COMM_WORLD, ierr)
     end if
  end do kpt_components

  yescomp=.false.
  kpt_orbitals: do jorb=1,orbs%norbp
     ikpt=orbs%iokpt(jorb)   
     do ikpts=1,orbs%nkptsp
        if (orbs%iskpts+ikpts == ikpt) yescomp=.true.
     end do
     if (.not. yescomp) then
        write(*,*)' ERROR: processor ', iproc,' kpt,',ikpt,&
             'not found in the component distribution'
        call MPI_ABORT(MPI_COMM_WORLD, ierr)
     end if
  end do kpt_orbitals

  !print *,'AAAAiproc',iproc,orbs%iskpts,orbs%iskpts+orbs%nkptsp

  !allocate communication arrays
  allocate(comms%nvctr_par(0:nproc-1,0:orbs%nkpts+ndebug),stat=i_stat)
  call memocc(i_stat,comms%nvctr_par,'nvctr_par',subname)

  allocate(comms%ncntd(0:nproc-1+ndebug),stat=i_stat)
  call memocc(i_stat,comms%ncntd,'ncntd',subname)

  allocate(comms%ncntt(0:nproc-1+ndebug),stat=i_stat)
  call memocc(i_stat,comms%ncntt,'ncntt',subname)
  allocate(comms%ndspld(0:nproc-1+ndebug),stat=i_stat)
  call memocc(i_stat,comms%ndspld,'ndspld',subname)
  allocate(comms%ndsplt(0:nproc-1+ndebug),stat=i_stat)
  call memocc(i_stat,comms%ndsplt,'ndsplt',subname)

  !assign the partition of the k-points to the communication array
  !calculate the number of componenets associated to the k-point
  do jproc=0,nproc-1
     comms%nvctr_par(jproc,0)=0
     do ikpt=1,orbs%nkpts
        comms%nvctr_par(jproc,0)=comms%nvctr_par(jproc,0)+&
             nvctr_par(jproc,ikpt) 
        comms%nvctr_par(jproc,ikpt)=nvctr_par(jproc,ikpt)
     end do
  end do
!!$  do ikpts=1,orbs%nkptsp
!!$     ikpt=orbs%iskpts+ikpts!orbs%ikptsp(ikpts)
!!$     do jproc=0,nproc-1
!!$        comms%nvctr_par(jproc,ikpts)=nvctr_par(jproc,ikpt) 
!!$     end do
!!$  end do

  !with this distribution the orbitals and the components are ordered following k-points
  !there must be no overlap for the components
  !here we will print out the k-points components distribution, in the transposed and in the direct way

  do jproc=0,nproc-1
     comms%ncntd(jproc)=0
     do ikpts=1,orbs%nkpts
        comms%ncntd(jproc)=comms%ncntd(jproc)+&
             nvctr_par(jproc,ikpts)*norb_par(iproc,ikpts)*orbs%nspinor
     end do
  end do
  comms%ndspld(0)=0
  do jproc=1,nproc-1
     comms%ndspld(jproc)=comms%ndspld(jproc-1)+comms%ncntd(jproc-1)
  end do
  !receive buffer
  do jproc=0,nproc-1
     comms%ncntt(jproc)=0
     do ikpts=1,orbs%nkpts
        comms%ncntt(jproc)=comms%ncntt(jproc)+&
             nvctr_par(iproc,ikpts)*norb_par(jproc,ikpts)*orbs%nspinor
     end do
  end do
  comms%ndsplt(0)=0
  do jproc=1,nproc-1
     comms%ndsplt(jproc)=comms%ndsplt(jproc-1)+comms%ncntt(jproc-1)
  end do

  !print *,'iproc,comms',iproc,comms%ncntd,comms%ndspld,comms%ncntt,comms%ndsplt

  i_all=-product(shape(nvctr_par))*kind(nvctr_par)
  deallocate(nvctr_par,stat=i_stat)
  call memocc(i_stat,i_all,'nvctr_par',subname)
  i_all=-product(shape(norb_par))*kind(norb_par)
  deallocate(norb_par,stat=i_stat)
  call memocc(i_stat,i_all,'norb_par',subname)
  i_all=-product(shape(mykpts))*kind(mykpts)
  deallocate(mykpts,stat=i_stat)
  call memocc(i_stat,i_all,'mykpts',subname)

  !calculate the dimension of the wavefunction
  !for the given processor (this is only the cubic strategy)
  orbs%npsidim_orbs=(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f)*orbs%norb_par(iproc,0)*orbs%nspinor
  orbs%npsidim_comp=sum(comms%ncntt(0:nproc-1))
    
!!$  orbs%npsidim=max((lr%wfd%nvctr_c+7*lr%wfd%nvctr_f)*orbs%norb_par(iproc,0)*orbs%nspinor,&
!!$       sum(comms%ncntt(0:nproc-1)))


  nB=max(orbs%npsidim_orbs,orbs%npsidim_comp)*8
  nMB=nB/1024/1024
  nKB=(nB-nMB*1024*1024)/1024
  nB=modulo(nB,1024)

  if (iproc == 0) write(*,'(1x,a,3(i5,a))') &
       'Wavefunctions memory occupation for root MPI process: ',&
       nMB,' MB ',nKB,' KB ',nB,' B'

END SUBROUTINE orbitals_communicators

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
