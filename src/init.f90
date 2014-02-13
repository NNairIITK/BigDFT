!> @file
!!  Routines to initialize the information about localisation regions
!! @author
!!    Copyright (C) 2007-2013 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 


!> Calculates the descriptor arrays and nvctrp
!! Calculates also the bounds arrays needed for convolutions
!! Refers this information to the global localisation region descriptor
subroutine createWavefunctionsDescriptors(iproc,hx,hy,hz,atoms,rxyz,radii_cf,&
     &   crmult,frmult,Glr,output_denspot)
  use module_base
  use module_types
  use yaml_output
  use module_interfaces, except_this_one => createWavefunctionsDescriptors
  implicit none
  !Arguments
  type(atoms_data), intent(in) :: atoms
  integer, intent(in) :: iproc
  real(gp), intent(in) :: hx,hy,hz,crmult,frmult
  real(gp), dimension(3,atoms%astruct%nat), intent(in) :: rxyz
  real(gp), dimension(atoms%astruct%ntypes,3), intent(in) :: radii_cf
  type(locreg_descriptors), intent(inout) :: Glr
  logical, intent(in), optional :: output_denspot
  !local variables
  character(len=*), parameter :: subname='createWavefunctionsDescriptors'
  integer :: i_all,i_stat,i1,i2,i3,iat
  integer :: n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3
  logical :: output_denspot_
  logical, dimension(:,:,:), pointer :: logrid_c,logrid_f

  call f_routine(id=subname)
  call timing(iproc,'CrtDescriptors','ON')
  

  !assign the dimensions to improve (a little) readability
  n1=Glr%d%n1
  n2=Glr%d%n2
  n3=Glr%d%n3
  nfl1=Glr%d%nfl1
  nfl2=Glr%d%nfl2
  nfl3=Glr%d%nfl3
  nfu1=Glr%d%nfu1
  nfu2=Glr%d%nfu2
  nfu3=Glr%d%nfu3

  !assign geocode and the starting points
  Glr%geocode=atoms%astruct%geocode

  ! determine localization region for all orbitals, but do not yet fill the descriptor arrays
  allocate(logrid_c(0:n1,0:n2,0:n3+ndebug),stat=i_stat)
  call memocc(i_stat,logrid_c,'logrid_c',subname)
  allocate(logrid_f(0:n1,0:n2,0:n3+ndebug),stat=i_stat)
  call memocc(i_stat,logrid_f,'logrid_f',subname)

  ! coarse/fine grid quantities
  if (atoms%astruct%ntypes >0) then
     call fill_logrid(atoms%astruct%geocode,n1,n2,n3,0,n1,0,n2,0,n3,0,atoms%astruct%nat,&
          &   atoms%astruct%ntypes,atoms%astruct%iatype,rxyz,radii_cf(1,1),crmult,hx,hy,hz,logrid_c)
     call fill_logrid(atoms%astruct%geocode,n1,n2,n3,0,n1,0,n2,0,n3,0,atoms%astruct%nat,&
          &   atoms%astruct%ntypes,atoms%astruct%iatype,rxyz,radii_cf(1,2),frmult,hx,hy,hz,logrid_f)
  else
     logrid_c=.true.
     logrid_f=.true.
  end if
  call wfd_from_grids(logrid_c,logrid_f,Glr)

  if (atoms%astruct%geocode == 'P' .and. .not. Glr%hybrid_on .and. Glr%wfd%nvctr_c /= (n1+1)*(n2+1)*(n3+1) ) then
     if (iproc ==0) then
        call yaml_warning('The coarse grid does not fill the entire periodic box')
        call yaml_comment('Errors due to translational invariance breaking may occur')
        !write(*,*) ' ERROR: the coarse grid does not fill the entire periodic box'
        !write(*,*) '          errors due to translational invariance breaking may occur'
        !stop
     end if
     if (GPUconv) then
        !        if (iproc ==0)then
        call yaml_warning('The code should be stopped for a GPU calculation')
        call yaml_comment('Since density is not initialised to 10^-20')
        !write(*,*) '          The code should be stopped for a GPU calculation     '
        !write(*,*) '          since density is not initialised to 10^-20               '
        !        end if
        stop
     end if
  end if

  output_denspot_ = .false.
  if (present(output_denspot)) output_denspot_ = output_denspot
  if (output_denspot_) then
     call export_grids("grid.xyz", atoms, rxyz, hx, hy, hz, n1, n2, n3, logrid_c, logrid_f)
  end if

  i_all=-product(shape(logrid_c))*kind(logrid_c)
  deallocate(logrid_c,stat=i_stat)
  call memocc(i_stat,i_all,'logrid_c',subname)
  i_all=-product(shape(logrid_f))*kind(logrid_f)
  deallocate(logrid_f,stat=i_stat)
  call memocc(i_stat,i_all,'logrid_f',subname)

  call timing(iproc,'CrtDescriptors','OF')
  call f_release_routine()
END SUBROUTINE createWavefunctionsDescriptors


subroutine wfd_from_grids(logrid_c, logrid_f, Glr)
   use module_base
   use module_types
   use yaml_output
   implicit none
   !Arguments
   type(locreg_descriptors), intent(inout) :: Glr
   logical, dimension(0:Glr%d%n1,0:Glr%d%n2,0:Glr%d%n3), intent(in) :: logrid_c,logrid_f
   !local variables
   character(len=*), parameter :: subname='wfd_from_grids'
   integer :: i_stat, i_all
   integer :: n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3

   !assign the dimensions to improve (a little) readability
   n1=Glr%d%n1
   n2=Glr%d%n2
   n3=Glr%d%n3
   nfl1=Glr%d%nfl1
   nfl2=Glr%d%nfl2
   nfl3=Glr%d%nfl3
   nfu1=Glr%d%nfu1
   nfu2=Glr%d%nfu2
   nfu3=Glr%d%nfu3

   !allocate kinetic bounds, only for free BC
   if (Glr%geocode == 'F' ) then
      allocate(Glr%bounds%kb%ibyz_c(2,0:n2,0:n3+ndebug),stat=i_stat)
      call memocc(i_stat,Glr%bounds%kb%ibyz_c,'Glr%bounds%kb%ibyz_c',subname)
      allocate(Glr%bounds%kb%ibxz_c(2,0:n1,0:n3+ndebug),stat=i_stat)
      call memocc(i_stat,Glr%bounds%kb%ibxz_c,'Glr%bounds%kb%ibxz_c',subname)
      allocate(Glr%bounds%kb%ibxy_c(2,0:n1,0:n2+ndebug),stat=i_stat)
      call memocc(i_stat,Glr%bounds%kb%ibxy_c,'Glr%bounds%kb%ibxy_c',subname)
      allocate(Glr%bounds%kb%ibyz_f(2,0:n2,0:n3+ndebug),stat=i_stat)
      call memocc(i_stat,Glr%bounds%kb%ibyz_f,'Glr%bounds%kb%ibyz_f',subname)
      allocate(Glr%bounds%kb%ibxz_f(2,0:n1,0:n3+ndebug),stat=i_stat)
      call memocc(i_stat,Glr%bounds%kb%ibxz_f,'Glr%bounds%kb%ibxz_f',subname)
      allocate(Glr%bounds%kb%ibxy_f(2,0:n1,0:n2+ndebug),stat=i_stat)
      call memocc(i_stat,Glr%bounds%kb%ibxy_f,'Glr%bounds%kb%ibxy_f',subname)
   end if

   ! Do the coarse region.
   call num_segkeys(n1,n2,n3,0,n1,0,n2,0,n3,logrid_c,Glr%wfd%nseg_c,Glr%wfd%nvctr_c)
   if (Glr%wfd%nseg_c == 0) then
      ! Check if the number of seg_c (Glr%wfd%nseg_c) > 0
      call yaml_warning('There is no coarse grid points (nseg_c=0)!')
      !write( *,*) ' ERROR: there is no coarse grid points (nseg_c=0)!'
      stop
   end if

   if (Glr%geocode == 'F') then
      call make_bounds(n1,n2,n3,logrid_c,Glr%bounds%kb%ibyz_c,Glr%bounds%kb%ibxz_c,Glr%bounds%kb%ibxy_c)
   end if

   ! Do the fine region.
   call num_segkeys(n1,n2,n3,0,n1,0,n2,0,n3,logrid_f,Glr%wfd%nseg_f,Glr%wfd%nvctr_f)
   if (Glr%geocode == 'F') then
      call make_bounds(n1,n2,n3,logrid_f,Glr%bounds%kb%ibyz_f,Glr%bounds%kb%ibxz_f,Glr%bounds%kb%ibxy_f)
   end if

   ! allocations for arrays holding the wavefunctions and their data descriptors
   call allocate_wfd(Glr%wfd)

   ! now fill the wavefunction descriptor arrays
   ! coarse grid quantities
   call segkeys(n1,n2,n3,0,n1,0,n2,0,n3,logrid_c,Glr%wfd%nseg_c, &
        & Glr%wfd%keyglob(1,1),Glr%wfd%keyvglob(1))
   ! fine grid quantities
   if (Glr%wfd%nseg_f > 0) then
      call segkeys(n1,n2,n3,0,n1,0,n2,0,n3,logrid_f,Glr%wfd%nseg_f, &
           & Glr%wfd%keyglob(1,Glr%wfd%nseg_c+1), Glr%wfd%keyvglob(Glr%wfd%nseg_c+1))
   end if
   !that is the point where the association is given
   !one should consider the possiblity of associating the 
   !arrays with f_associate
!!$   i_all = -product(shape(Glr%wfd%keygloc))*kind(Glr%wfd%keygloc)
!!$   deallocate(Glr%wfd%keygloc,stat=i_stat)
!!$   call memocc(i_stat,i_all,'Glr%wfd%keygloc',subname)
   call f_free_ptr(Glr%wfd%keygloc)
   Glr%wfd%keygloc => Glr%wfd%keyglob

!!$   i_all = -product(shape(Glr%wfd%keyvloc))*kind(Glr%wfd%keyvloc)
!!$   deallocate(Glr%wfd%keyvloc,stat=i_stat)
!!$   call memocc(i_stat,i_all,'Glr%wfd%keyvloc',subname)
   call f_free_ptr(Glr%wfd%keyvloc)
   Glr%wfd%keyvloc => Glr%wfd%keyvglob
 
   ! Copy the information of keyglob to keygloc for Glr (just pointing leads to problem during the deallocation of wfd)
!!$   do i = lbound(Glr%wfd%keyglob,1),ubound(Glr%wfd%keyglob,1)
!!$      do j = lbound(Glr%wfd%keyglob,2),ubound(Glr%wfd%keyglob,2)
!!$         Glr%wfd%keygloc(i,j) = Glr%wfd%keyglob(i,j)
!!$      end do
!!$   end do
   
 !for free BC admits the bounds arrays
   if (Glr%geocode == 'F' ) then
      !allocate grow, shrink and real bounds
      allocate(Glr%bounds%gb%ibzxx_c(2,0:n3,-14:2*n1+16+ndebug),stat=i_stat)
      call memocc(i_stat,Glr%bounds%gb%ibzxx_c,'Glr%bounds%gb%ibzxx_c',subname)
      allocate(Glr%bounds%gb%ibxxyy_c(2,-14:2*n1+16,-14:2*n2+16+ndebug),stat=i_stat)
      call memocc(i_stat,Glr%bounds%gb%ibxxyy_c,'Glr%bounds%gb%ibxxyy_c',subname)
      allocate(Glr%bounds%gb%ibyz_ff(2,nfl2:nfu2,nfl3:nfu3+ndebug),stat=i_stat)
      call memocc(i_stat,Glr%bounds%gb%ibyz_ff,'Glr%bounds%gb%ibyz_ff',subname)
      allocate(Glr%bounds%gb%ibzxx_f(2,nfl3:nfu3,2*nfl1-14:2*nfu1+16+ndebug),stat=i_stat)
      call memocc(i_stat,Glr%bounds%gb%ibzxx_f,'Glr%bounds%gb%ibzxx_f',subname)
      allocate(Glr%bounds%gb%ibxxyy_f(2,2*nfl1-14:2*nfu1+16,2*nfl2-14:2*nfu2+16+ndebug),stat=i_stat)
      call memocc(i_stat,Glr%bounds%gb%ibxxyy_f,'Glr%bounds%gb%ibxxyy_f',subname)

      allocate(Glr%bounds%sb%ibzzx_c(2,-14:2*n3+16,0:n1+ndebug),stat=i_stat)
      call memocc(i_stat,Glr%bounds%sb%ibzzx_c,'Glr%bounds%sb%ibzzx_c',subname)
      allocate(Glr%bounds%sb%ibyyzz_c(2,-14:2*n2+16,-14:2*n3+16+ndebug),stat=i_stat)
      call memocc(i_stat,Glr%bounds%sb%ibyyzz_c,'Glr%bounds%sb%ibyyzz_c',subname)
      allocate(Glr%bounds%sb%ibxy_ff(2,nfl1:nfu1,nfl2:nfu2+ndebug),stat=i_stat)
      call memocc(i_stat,Glr%bounds%sb%ibxy_ff,'Glr%bounds%sb%ibxy_ff',subname)
      allocate(Glr%bounds%sb%ibzzx_f(2,-14+2*nfl3:2*nfu3+16,nfl1:nfu1+ndebug),stat=i_stat)
      call memocc(i_stat,Glr%bounds%sb%ibzzx_f,'Glr%bounds%sb%ibzzx_f',subname)
      allocate(Glr%bounds%sb%ibyyzz_f(2,-14+2*nfl2:2*nfu2+16,-14+2*nfl3:2*nfu3+16+ndebug),stat=i_stat)
      call memocc(i_stat,Glr%bounds%sb%ibyyzz_f,'Glr%bounds%sb%ibyyzz_f',subname)

      allocate(Glr%bounds%ibyyzz_r(2,-14:2*n2+16,-14:2*n3+16+ndebug),stat=i_stat)
      call memocc(i_stat,Glr%bounds%ibyyzz_r,'Glr%bounds%ibyyzz_r',subname)

      call make_all_ib(n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,&
         &   Glr%bounds%kb%ibxy_c,Glr%bounds%sb%ibzzx_c,Glr%bounds%sb%ibyyzz_c,&
         &   Glr%bounds%kb%ibxy_f,Glr%bounds%sb%ibxy_ff,Glr%bounds%sb%ibzzx_f,Glr%bounds%sb%ibyyzz_f,&
         &   Glr%bounds%kb%ibyz_c,Glr%bounds%gb%ibzxx_c,Glr%bounds%gb%ibxxyy_c,&
         &   Glr%bounds%kb%ibyz_f,Glr%bounds%gb%ibyz_ff,Glr%bounds%gb%ibzxx_f,Glr%bounds%gb%ibxxyy_f,&
         &   Glr%bounds%ibyyzz_r)

   end if

   if (Glr%geocode == 'P' .and. Glr%hybrid_on) then
      call make_bounds_per(n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,Glr%bounds,Glr%wfd)
      call make_all_ib_per(n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,&
         &   Glr%bounds%kb%ibxy_f,Glr%bounds%sb%ibxy_ff,Glr%bounds%sb%ibzzx_f,Glr%bounds%sb%ibyyzz_f,&
         &   Glr%bounds%kb%ibyz_f,Glr%bounds%gb%ibyz_ff,Glr%bounds%gb%ibzxx_f,Glr%bounds%gb%ibxxyy_f)
   endif


end subroutine wfd_from_grids

!> Determine localization region for all projectors, but do not yet fill the descriptor arrays
subroutine createProjectorsArrays(lr,rxyz,at,orbs,&
         & radii_cf,cpmult,fpmult,hx,hy,hz,dry_run,nl,proj_G)
  use module_base
  use psp_projectors
  use module_types
  use gaussians, only: gaussian_basis
  implicit none
  real(gp), intent(in) :: cpmult,fpmult,hx,hy,hz
  type(locreg_descriptors),intent(in) :: lr
  type(atoms_data), intent(in) :: at
  type(orbitals_data), intent(in) :: orbs
  real(gp), dimension(3,at%astruct%nat), intent(in) :: rxyz
  real(gp), dimension(at%astruct%ntypes,3), intent(in) :: radii_cf
  type(gaussian_basis),dimension(at%astruct%ntypes),intent(in) :: proj_G
  type(DFT_PSP_projectors), intent(out) :: nl
  logical, intent(in) :: dry_run !< .true. to compute the size only and don't allocate
  !local variables
  character(len=*), parameter :: subname='createProjectorsArrays'
  integer :: n1,n2,n3,nl1,nl2,nl3,nu1,nu2,nu3,mseg,mproj
  integer :: iat,i_stat,i_all,iseg
  logical, dimension(:,:,:), allocatable :: logrid
  call f_routine(id=subname)

  !start from a null structure
  nl=DFT_PSP_projectors_null()

  !allocate the different localization regions of the projectors
  nl%natoms=at%astruct%nat
  !for a structure let the allocator crash when allocating
  allocate(nl%pspd(at%astruct%nat))
  do iat=1,at%astruct%nat
     nl%pspd(iat)=nonlocal_psp_descriptors_null()
  end do

  ! define the region dimensions
  n1 = lr%d%n1
  n2 = lr%d%n2
  n3 = lr%d%n3

  ! determine localization region for all projectors, but do not yet fill the descriptor arrays
  logrid=f_malloc((/0.to.n1,0.to.n2,0.to.n3/),id='logrid')

  call localize_projectors(n1,n2,n3,hx,hy,hz,cpmult,fpmult,&
       rxyz,radii_cf,logrid,at,orbs,nl,proj_G)

  if (dry_run) then
     call f_free(logrid)
     return
  end if

  !here the allocation is possible
  do iat=1,nl%natoms
     !also the fact of allocating pointers with size zero has to be discussed
     !for the moments the bounds are not needed for projectors
     call allocate_wfd(nl%pspd(iat)%plr%wfd)
  end do
  nl%proj=f_malloc0_ptr(nl%nprojel,id='proj')

  ! After having determined the size of the projector descriptor arrays fill them
  do iat=1,at%astruct%nat
     if (nl%pspd(iat)%mproj > 0) then 

        call bounds_to_plr_limits(.false.,1,nl%pspd(iat)%plr,&
             nl1,nl2,nl3,nu1,nu2,nu3)         

!!$        !most likely the call can here be replaced by
!!$        call fill_logrid(at%astruct%geocode,n1,n2,n3,nl1,nu1,nl2,nu2,nl3,nu3,0,1,  &
!!$             1,(/1/),rxyz(1,iat),radii_cf(at%astruct%iatype(iat),3),&
!!$             cpmult,hx,hy,hz,logrid)
!!$        !which make radiicf and rxyz the only external data needed

        call fill_logrid(at%astruct%geocode,n1,n2,n3,nl1,nu1,nl2,nu2,nl3,nu3,0,1,  &
             at%astruct%ntypes,at%astruct%iatype(iat),rxyz(1,iat),radii_cf(1,3),&
             cpmult,hx,hy,hz,logrid)

        call segkeys(n1,n2,n3,nl1,nu1,nl2,nu2,nl3,nu3,logrid,&
             nl%pspd(iat)%plr%wfd%nseg_c,&
             nl%pspd(iat)%plr%wfd%keyglob(1,1),nl%pspd(iat)%plr%wfd%keyvglob(1))

        call transform_keyglob_to_keygloc(lr,nl%pspd(iat)%plr,nl%pspd(iat)%plr%wfd%nseg_c,&
             nl%pspd(iat)%plr%wfd%keyglob(1,1),nl%pspd(iat)%plr%wfd%keygloc(1,1))

        ! fine grid quantities
        call bounds_to_plr_limits(.false.,2,nl%pspd(iat)%plr,&
             nl1,nl2,nl3,nu1,nu2,nu3)         

        call fill_logrid(at%astruct%geocode,n1,n2,n3,nl1,nu1,nl2,nu2,nl3,nu3,0,1,  &
             & at%astruct%ntypes,at%astruct%iatype(iat),rxyz(1,iat),radii_cf(1,2),&
             fpmult,hx,hy,hz,logrid)

        mseg=nl%pspd(iat)%plr%wfd%nseg_f
        iseg=nl%pspd(iat)%plr%wfd%nseg_c+1

        if (mseg > 0) then
           call segkeys(n1,n2,n3,nl1,nu1,nl2,nu2,nl3,nu3,  & 
                logrid,mseg,nl%pspd(iat)%plr%wfd%keyglob(1,iseg),&
                nl%pspd(iat)%plr%wfd%keyvglob(iseg))

           call transform_keyglob_to_keygloc(lr,nl%pspd(iat)%plr,mseg,nl%pspd(iat)%plr%wfd%keyglob(1,iseg),&
                nl%pspd(iat)%plr%wfd%keygloc(1,iseg)) 
        end if
     endif
  enddo

  call f_free(logrid)
  !fill the projectors if the strategy is a distributed calculation
  if (.not. nl%on_the_fly) then
     !calculate the wavelet expansion of projectors
     call fill_projectors(lr,hx,hy,hz,at,orbs,rxyz,nl,0)
  end if

  call f_release_routine()
END SUBROUTINE createProjectorsArrays


        !!$subroutine initRhoPot(iproc, nproc, Glr, hxh, hyh, hzh, atoms, rxyz, crmult, frmult, radii, nspin, ixc, rho_commun, rhodsc, nscatterarr, ngatherarr, pot_ion)
        !!$  use module_base
        !!$  use module_types
        !!$
        !!$  implicit none
        !!$
        !!$  integer, intent(in) :: iproc, nproc
        !!$
        !!$  integer :: i_stat
        !!$
        !!$END SUBROUTINE initRhoPot

        subroutine input_wf_empty(iproc, nproc, psi, hpsi, psit, orbs, &
              & band_structure_filename, input_spin, atoms, d, denspot)
          use module_defs
          use module_types
          use yaml_output
          use module_interfaces, except_this_one => input_wf_empty
          implicit none
          integer, intent(in) :: iproc, nproc
          type(orbitals_data), intent(in) :: orbs
          character(len = *), intent(in) :: band_structure_filename
          integer, intent(in) :: input_spin
          type(atoms_data), intent(in) :: atoms
          type(grid_dimensions), intent(in) :: d
          type(DFT_local_fields), intent(inout) :: denspot
          real(wp), dimension(:), pointer :: psi
          real(kind=8), dimension(:), pointer :: hpsi, psit

          character(len = *), parameter :: subname = "input_wf_empty"
          integer :: i_stat, i_all, nspin, n1i, n2i, n3i, ispin, ierr
          real(gp) :: hxh, hyh, hzh

          !allocate fake psit and hpsi
          allocate(hpsi(max(orbs%npsidim_comp,orbs%npsidim_orbs)+ndebug),stat=i_stat)
          call memocc(i_stat,hpsi,'hpsi',subname)
          if (nproc > 1) then
             allocate(psit(max(orbs%npsidim_comp,orbs%npsidim_orbs)+ndebug),stat=i_stat)
             call memocc(i_stat,psit,'psit',subname)
          else
             psit => psi
          end if
          !fill the rhopot array with the read potential if needed
          if (trim(band_structure_filename) /= '') then
             !only the first processor should read this
             if (iproc == 0) then
        call yaml_map('Reading local potential from file:',trim(band_structure_filename))
        !write(*,'(1x,a)')'Reading local potential from file:'//trim(band_structure_filename)
        call read_density(trim(band_structure_filename),atoms%astruct%geocode,&
                     n1i,n2i,n3i,nspin,hxh,hyh,hzh,denspot%Vloc_KS)
                if (nspin /= input_spin) stop
             else
                allocate(denspot%Vloc_KS(1,1,1,input_spin+ndebug),stat=i_stat)
                call memocc(i_stat,denspot%Vloc_KS,'Vloc_KS',subname)
             end if

             if (nproc > 1) then
                do ispin=1,input_spin
                   call MPI_SCATTERV(denspot%Vloc_KS(1,1,1,ispin),&
                        denspot%dpbox%ngatherarr(0,1),denspot%dpbox%ngatherarr(0,2),&
                        mpidtypw,denspot%rhov((ispin-1)*&
                        d%n1i*d%n2i*denspot%dpbox%n3p+1),&
                        d%n1i*d%n2i*denspot%dpbox%n3p,mpidtypw,0,&
                        bigdft_mpi%mpi_comm,ierr)
                end do
             else
                call vcopy(d%n1i*d%n2i*d%n3i*input_spin,&
                     denspot%Vloc_KS(1,1,1,1),1,denspot%rhov(1),1)
             end if
             !now the meaning is KS potential
             call denspot_set_rhov_status(denspot, KS_POTENTIAL, 0, iproc, nproc)

             i_all=-product(shape(denspot%Vloc_KS))*kind(denspot%Vloc_KS)
             deallocate(denspot%Vloc_KS,stat=i_stat)
             call memocc(i_stat,i_all,'Vloc_KS',subname)

             !add pot_ion potential to the local_potential
             !do ispin=1,in%nspin
             !   !spin up and down together with the XC part
             !   call axpy(Lzd%Glr%d%n1i*Lzd%Glr%d%n2i*n3p,1.0_dp,pot_ion(1),1,&
             !        rhopot((ispin-1)*Lzd%Glr%d%n1i*Lzd%Glr%d%n2i*n3p+1),1)
             !end do
          end if
        END SUBROUTINE input_wf_empty


        !> Random initialisation of the wavefunctions
        !! The initialization of only the scaling function coefficients should be considered
        subroutine input_wf_random(psi, orbs)
          use module_defs
          use module_types
          implicit none

          type(orbitals_data), intent(inout) :: orbs
          real(wp), dimension(:), pointer :: psi

          integer :: icoeff,jorb,iorb,nvctr
          integer :: idum=0
          real(kind=4) :: tt,builtin_rand

          if (max(orbs%npsidim_comp,orbs%npsidim_orbs)>1) &
               call to_zero(max(orbs%npsidim_comp,orbs%npsidim_orbs),psi(1))

          !Fill randomly the wavefunctions coefficients for the orbitals considered
          nvctr=orbs%npsidim_orbs/(orbs%nspinor*orbs%norbp)
          do icoeff=1,nvctr !tt not dependent of iproc
             !Be sure to call always a different random number, per orbital
             do jorb=1,orbs%isorb*orbs%nspinor
                tt=builtin_rand(idum) !call random_number(tt)
             end do
             do iorb=1,orbs%norbp*orbs%nspinor
                tt=builtin_rand(idum) !call random_number(tt)
                psi(icoeff+(iorb-1)*nvctr)=real(tt,wp)
             end do
             do iorb=(orbs%isorb+orbs%norbp)*orbs%nspinor+1,orbs%norb*orbs%nkpts*orbs%nspinor
                tt=builtin_rand(idum) !call random_number(tt)
             end do
          end do

          orbs%eval(1:orbs%norb*orbs%nkpts)=-0.5d0

        END SUBROUTINE input_wf_random


        !> Initialisation of the wavefunctions via import gaussians from CP2K
        subroutine input_wf_cp2k(iproc, nproc, nspin, atoms, rxyz, Lzd, &
             & psi, orbs)
          use module_defs
          use module_types
  use yaml_output
          use module_interfaces, except_this_one => input_wf_cp2k
          implicit none

          integer, intent(in) :: iproc, nproc, nspin
          type(atoms_data), intent(in) :: atoms
  real(gp), dimension(3, atoms%astruct%nat), intent(in) :: rxyz
          type(local_zone_descriptors), intent(in) :: Lzd
          type(orbitals_data), intent(inout) :: orbs
          real(wp), dimension(:), pointer :: psi

          character(len = *), parameter :: subname = "input_wf_cp2k"
          integer :: i_stat, i_all
          type(gaussian_basis) :: gbd
          real(wp), dimension(:,:), pointer :: gaucoeffs

          !import gaussians form CP2K (data in files gaubasis.dat and gaucoeff.dat)
          !and calculate eigenvalues
          if (nspin /= 1) then
             if (iproc==0) then
        call yaml_warning('Gaussian importing is possible only for non-spin polarised calculations')
        call yaml_comment('The reading rules of CP2K files for spin-polarised orbitals are not implemented')
        !write(*,'(1x,a)') 'Gaussian importing is possible only for non-spin polarised calculations'
        !write(*,'(1x,a)') 'The reading rules of CP2K files for spin-polarised orbitals are not implemented'
             end if
             stop
          end if

          call parse_cp2k_files(iproc,'gaubasis.dat','gaucoeff.dat',&
       atoms%astruct%nat,atoms%astruct%ntypes,orbs,atoms%astruct%iatype,rxyz,gbd,gaucoeffs)

          call gaussians_to_wavelets_new(iproc,nproc,Lzd,orbs,gbd,gaucoeffs,psi)

          !deallocate gaussian structure and coefficients
          call deallocate_gwf(gbd,subname)
          i_all=-product(shape(gaucoeffs))*kind(gaucoeffs)
          deallocate(gaucoeffs,stat=i_stat)
          call memocc(i_stat,i_all,'gaucoeffs',subname)
          nullify(gbd%rxyz)

          !call dual_gaussian_coefficients(orbs%norbp,gbd,gaucoeffs)
          orbs%eval(1:orbs%norb*orbs%nkpts)=-0.5d0

        END SUBROUTINE input_wf_cp2k

subroutine input_wf_memory_history(iproc,orbs,atoms,wfn_history,istep_history,oldpsis,rxyz,Lzd,psi)
  use module_base
  use module_types
  use module_interfaces
  use yaml_output
  implicit none
  integer, intent(in) :: iproc,wfn_history
  type(atoms_data), intent(in) :: atoms
  real(gp), dimension(3, atoms%astruct%nat), intent(in) :: rxyz
  type(local_zone_descriptors),intent(in) :: lzd
  type(orbitals_data), intent(in) :: orbs
  type(old_wavefunction), dimension(0:wfn_history+1), intent(inout) :: oldpsis
  integer, intent(inout) :: istep_history
  real(wp), dimension(Lzd%Glr%wfd%nvctr_c+7*Lzd%Glr%wfd%nvctr_f,orbs%nspinor*orbs%norbp), intent(out) :: psi
  !local variables
  character(len=*), parameter :: subname='input_wf_memory_history'
  integer :: i_stat,i_all,istep,jstep,nvctr
  real(wp), dimension(:,:), allocatable :: psi_tmp
  real(gp), dimension(3:9) :: kappa,alpha
  real(gp), dimension(0:9,3:9) :: c

  !set the coefficients
  c=0.0_gp
  kappa(3)=1.69_gp
  kappa(4)=1.75_gp
  kappa(5)=1.82_gp
  kappa(6)=1.84_gp
  kappa(7)=1.86_gp
  kappa(8)=1.88_gp  
  kappa(9)=1.89_gp

  alpha(3)=150.e-3_gp
  alpha(4)=57.e-3_gp
  alpha(5)=18.e-3_gp
  alpha(6)=5.5e-3_gp
  alpha(7)=1.6e-3_gp
  alpha(8)=.44e-3_gp  
  alpha(9)=.12e-3_gp

  c(0:3,3)=alpha(3)*(/-2._gp,3._gp,0._gp,-1._gp /)
  c(0:4,4)=alpha(4)*(/-3._gp,6._gp,-2._gp,-2._gp,1._gp /)
  c(0:5,5)=alpha(5)*(/-6._gp,14._gp,-8._gp,-3._gp,4._gp,-1._gp /)
  c(0:6,6)=alpha(6)*(/-14._gp,36._gp,-27._gp,-2._gp,12._gp,-6._gp,1._gp /)
  c(0:7,7)=alpha(7)*(/-36._gp,99._gp,-88._gp,11._gp,32._gp,-25._gp,8._gp,-1._gp /)
  c(0:8,8)=alpha(8)*(/-99._gp,286._gp,-286._gp,78._gp,78._gp,-90._gp,42._gp,-10._gp,1._gp /)  
  c(0:9,9)=alpha(9)*(/-286._gp,858._gp,-936._gp,364._gp,168._gp,-300._gp,184._gp,-63._gp,12._gp,-1._gp /)  
  !rework the coefficients for the first two elements
  do istep=3,9
     c(0,istep)=c(0,istep)+2._gp-kappa(istep)
     c(1,istep)=c(1,istep)-1._gp
  end do
  !number of componenets
  nvctr=(Lzd%Glr%wfd%nvctr_c+7*Lzd%Glr%wfd%nvctr_f)*orbs%nspinor*orbs%norbp
  !check if history has not yet been filled
  if (istep_history <= wfn_history) then
     !if so, copy the SCF wfn, which is in the last position, in the corresponding history place

     call old_wavefunction_set(oldpsis(istep_history),&
          atoms%astruct%nat,orbs%norbp*orbs%nspinor,&
          oldpsis(wfn_history+1)%Lzd,oldpsis(wfn_history+1)%rxyz,&
          oldpsis(wfn_history+1)%psi)
     !check if it is the first restart
     if (istep_history == 0) then
        do istep=1,wfn_history
                call old_wavefunction_set(oldpsis(istep),&
                atoms%astruct%nat,orbs%norbp*orbs%nspinor,&
                oldpsis(wfn_history+1)%Lzd,oldpsis(wfn_history+1)%rxyz,&
                oldpsis(wfn_history+1)%psi)
        end do
     end if
  end if
if (iproc==0)call yaml_map('Previous SCF wfn copied',.true.)   
  !put to zero the wavefunction
  if (nvctr>0) call to_zero(nvctr,psi(1,1))

  !calculate the reformat with history
  allocate(psi_tmp(Lzd%Glr%wfd%nvctr_c+7*Lzd%Glr%wfd%nvctr_f,orbs%nspinor*orbs%norbp+ndebug),stat=i_stat)
  call memocc(i_stat,psi_tmp,'psi_tmp',subname)

  !first reformat the previous SCF step
  istep=wfn_history+1
  call reformatmywaves(iproc,orbs,atoms,&
       oldpsis(istep)%Lzd%hgrids(1),oldpsis(istep)%Lzd%hgrids(2),oldpsis(istep)%Lzd%hgrids(3),&
       oldpsis(istep)%Lzd%Glr%d%n1,oldpsis(istep)%Lzd%Glr%d%n2,oldpsis(istep)%Lzd%Glr%d%n3,&
       oldpsis(istep)%rxyz,oldpsis(istep)%Lzd%Glr%wfd,&
       oldpsis(istep)%psi,Lzd%hgrids(1),Lzd%hgrids(2),Lzd%hgrids(3),&
       Lzd%Glr%d%n1,Lzd%Glr%d%n2,Lzd%Glr%d%n3,rxyz,Lzd%Glr%wfd,psi_tmp)
  if (nvctr>0) call axpy(nvctr,kappa(wfn_history),psi_tmp(1,1),1,psi(1,1),1)
  call yaml_map('Reformat Previous SCF wfn',.true.)   
  !then the reformatting step based on history
  do jstep=0,wfn_history
     istep=modulo(modulo(istep_history,wfn_history+1)-jstep,wfn_history+1)
     call reformatmywaves(iproc,orbs,atoms,&
          oldpsis(istep)%Lzd%hgrids(1),oldpsis(istep)%Lzd%hgrids(2),oldpsis(istep)%Lzd%hgrids(3),&
          oldpsis(istep)%Lzd%Glr%d%n1,oldpsis(istep)%Lzd%Glr%d%n2,oldpsis(istep)%Lzd%Glr%d%n3,&
          oldpsis(istep)%rxyz,oldpsis(istep)%Lzd%Glr%wfd,&
          oldpsis(istep)%psi,Lzd%hgrids(1),Lzd%hgrids(2),Lzd%hgrids(3),&
          Lzd%Glr%d%n1,Lzd%Glr%d%n2,Lzd%Glr%d%n3,rxyz,Lzd%Glr%wfd,psi_tmp)
     if (nvctr>0) call axpy(nvctr,c(jstep,wfn_history),psi_tmp(1,1),1,psi(1,1),1)
     if (iproc==0)call yaml_map('Reformat Input wfn of Iter.',jstep,advance='no')   
     if (iproc==0)call yaml_comment('Position:'//trim(yaml_toa(istep))//', Step'//trim(yaml_toa(istep_history)))
  end do
  i_all=-product(shape(psi_tmp))*kind(psi_tmp)
  deallocate(psi_tmp,stat=i_stat)
  call memocc(i_stat,i_all,'psi_tmp',subname)

  !increase the iteration step
  istep_history=istep_history+1
  if (istep_history > wfn_history+1) then
     istep=modulo(istep_history,wfn_history+1)
     !and save the input wfn in the history
     call old_wavefunction_set(oldpsis(istep),&
          atoms%astruct%nat,orbs%norbp*orbs%nspinor,&
          Lzd,rxyz,psi)
  end if

end subroutine input_wf_memory_history

subroutine input_wf_memory(iproc, atoms, &
     & rxyz_old, hx_old, hy_old, hz_old, d_old, wfd_old, psi_old, &
     & rxyz, hx, hy, hz, d, wfd, psi, orbs)
  use module_defs
  use module_types
  use module_interfaces, except_this_one => input_wf_memory
  implicit none

  integer, intent(in) :: iproc
  type(atoms_data), intent(in) :: atoms
  real(gp), dimension(3, atoms%astruct%nat), intent(in) :: rxyz, rxyz_old
  real(gp), intent(in) :: hx, hy, hz, hx_old, hy_old, hz_old
  type(grid_dimensions), intent(in) :: d, d_old
  type(wavefunctions_descriptors), intent(in) :: wfd
  type(wavefunctions_descriptors), intent(inout) :: wfd_old
  type(orbitals_data), intent(in) :: orbs
  real(wp), dimension(:), pointer :: psi, psi_old

  character(len = *), parameter :: subname = "input_wf_memory"
  integer :: i_stat, i_all

  !these parts should be reworked for the non-collinear spin case
  call reformatmywaves(iproc,orbs,atoms,hx_old,hy_old,hz_old,&
       d_old%n1,d_old%n2,d_old%n3,rxyz_old,wfd_old,psi_old,hx,hy,hz,&
       & d%n1,d%n2,d%n3,rxyz,wfd,psi)


  i_all=-product(shape(psi_old))*kind(psi_old)
  deallocate(psi_old,stat=i_stat)
  call memocc(i_stat,i_all,'psi_old',subname)
END SUBROUTINE input_wf_memory



subroutine input_memory_linear(iproc, nproc, at, KSwfn, tmb, tmb_old, denspot, input, &
           rxyz_old, rxyz, denspot0, energs, nlpsp, GPU, ref_frags)

  use module_base
  use module_types
  use module_interfaces, except_this_one => input_memory_linear
  use module_fragments
  use yaml_output
  implicit none

  ! Calling arguments
  integer,intent(in) :: iproc, nproc
  type(atoms_data), intent(inout) :: at
  type(DFT_wavefunction),intent(inout):: KSwfn
  type(DFT_wavefunction),intent(inout):: tmb, tmb_old
  type(DFT_local_fields), intent(inout) :: denspot
  type(input_variables),intent(in):: input
  real(gp),dimension(3,at%astruct%nat),intent(in) :: rxyz_old, rxyz
  real(8),dimension(max(denspot%dpbox%ndims(1)*denspot%dpbox%ndims(2)*denspot%dpbox%n3p,1)),intent(out):: denspot0
  type(energy_terms),intent(inout):: energs
  type(DFT_PSP_projectors), intent(inout) :: nlpsp
  type(GPU_pointers), intent(inout) :: GPU
  type(system_fragment), dimension(:), intent(in) :: ref_frags

  ! Local variables
  integer :: ndim_old, ndim, iorb, iiorb, ilr, i_stat, i_all, ilr_old, iiat
         logical:: overlap_calculated
  real(wp), allocatable, dimension(:) :: norm
  type(fragment_transformation), dimension(:), pointer :: frag_trans
  character(len=*),parameter:: subname='input_memory_linear'

  ! Determine size of phi_old and phi
  ndim_old=0
  ndim=0
  do iorb=1,tmb%orbs%norbp
      iiorb=tmb%orbs%isorb+iorb
      ilr=tmb%orbs%inwhichlocreg(iiorb)
      !ilr_old=tmb_old%orbs%inwhichlocreg(iiorb)
      ilr_old=ilr
              !!write(*,*) '###### input_memory_linear: iiorb, ilr', iiorb, ilr
      ndim_old=ndim_old+tmb_old%lzd%llr(ilr_old)%wfd%nvctr_c+7*tmb_old%lzd%llr(ilr_old)%wfd%nvctr_f
      ndim=ndim+tmb%lzd%llr(ilr)%wfd%nvctr_c+7*tmb%lzd%llr(ilr)%wfd%nvctr_f
  end do

  ! Reformat the support functions if we are not using FOE. Otherwise an AO
  ! input guess wil be done below.
  if (input%lin%scf_mode/=LINEAR_FOE) then

     ! define fragment transformation - should eventually be done automatically...
     allocate(frag_trans(tmb%orbs%norbp))

     do iorb=1,tmb%orbs%norbp
         iiat=tmb%orbs%onwhichatom(iorb+tmb%orbs%isorb)
         frag_trans(iorb)%theta=0.0d0*(4.0_gp*atan(1.d0)/180.0_gp)
         frag_trans(iorb)%rot_axis=(/1.0_gp,0.0_gp,0.0_gp/)
         frag_trans(iorb)%rot_center(:)=rxyz_old(:,iiat)
         frag_trans(iorb)%rot_center_new(:)=rxyz(:,iiat)
     end do

     call reformat_supportfunctions(iproc,at,rxyz_old,rxyz,.true.,tmb,ndim_old,tmb_old%lzd,frag_trans,&
          tmb_old%psi,input%dir_output,input%frag,ref_frags)

     deallocate(frag_trans)
  end if
          !!write(*,*) 'after reformat_supportfunctions, iproc',iproc

  ! need the input guess eval for preconditioning as they won't be recalculated
  !if(input%lin%scf_mode==LINEAR_DIRECT_MINIMIZATION) then
  ! needed for Pulay as now using tmb rather Kswfn evals
     do iorb=1,tmb%orbs%norb
        tmb%orbs%eval(iorb) = tmb_old%orbs%eval(iorb)
     end do
  !end if

  call deallocate_local_zone_descriptors(tmb_old%lzd, subname)
  call deallocate_orbitals_data(tmb_old%orbs, subname)

  i_all = -product(shape(tmb_old%psi))*kind(tmb_old%psi)
  deallocate(tmb_old%psi,stat=i_stat)
  call memocc(i_stat,i_all,'tmb_old%psi',subname)

  !!call deallocate_wfd(tmb_old%lzd%glr%wfd,subname)
  !!do ilr=1,tmb_old%lzd%nlr
  !!    call deallocate_wfd(tmb_old%lzd%llr(ilr)%wfd,subname)
          !!end do

          ! Copy the coefficients
  if (input%lin%scf_mode/=LINEAR_FOE) then
      call dcopy(tmb%orbs%norb*tmb%orbs%norb, tmb_old%coeff(1,1), 1, tmb%coeff(1,1), 1)
  end if
          !!write(*,*) 'after dcopy, iproc',iproc

  if (associated(tmb_old%coeff)) then
      i_all=-product(shape(tmb_old%coeff))*kind(tmb_old%coeff)
      deallocate(tmb_old%coeff, stat=i_stat)
      call memocc(i_stat, i_all, 'tmb_old%coeff', subname)
  end if

  ! MOVE LATER 
  !!if (associated(tmb_old%linmat%denskern%matrix_compr)) then
  !!   i_all=-product(shape(tmb_old%linmat%denskern%matrix_compr))*kind(tmb_old%linmat%denskern%matrix_compr)
  !!   deallocate(tmb_old%linmat%denskern%matrix_compr, stat=i_stat)
  !!   call memocc(i_stat, i_all, 'tmb_old%linmat%denskern%matrix_compr', subname)
  !!end if
  if (associated(tmb_old%linmat%denskern_large%matrix_compr)) then
     i_all=-product(shape(tmb_old%linmat%denskern_large%matrix_compr))*kind(tmb_old%linmat%denskern_large%matrix_compr)
     deallocate(tmb_old%linmat%denskern_large%matrix_compr, stat=i_stat)
     call memocc(i_stat, i_all, 'tmb_old%linmat%denskern_large%matrix_compr', subname)
  end if

  ! destroy it all together here - don't have all comms arrays
  !call destroy_DFT_wavefunction(tmb_old)

          !!write(*,*) 'after deallocate, iproc', iproc

   ! normalize tmbs - only really needs doing if we reformatted, but will need to calculate transpose after anyway

   ! Normalize the input guess. If FOE is used, the input guess will be generated below.
   if (input%lin%scf_mode/=LINEAR_FOE) then
       tmb%can_use_transposed=.true.
       overlap_calculated=.false.
       allocate(tmb%psit_c(sum(tmb%collcom%nrecvcounts_c)), stat=i_stat)
       call memocc(i_stat, tmb%psit_c, 'tmb%psit_c', subname)

       allocate(tmb%psit_f(7*sum(tmb%collcom%nrecvcounts_f)), stat=i_stat)
       call memocc(i_stat, tmb%psit_f, 'tmb%psit_f', subname)


       call transpose_localized(iproc, nproc, tmb%npsidim_orbs, tmb%orbs, tmb%collcom, &
            tmb%psi, tmb%psit_c, tmb%psit_f, tmb%lzd)

       ! normalize psi
       allocate(norm(tmb%orbs%norb), stat=i_stat)
       call memocc(i_stat, norm, 'norm', subname)

       call normalize_transposed(iproc, nproc, tmb%orbs, tmb%collcom, tmb%psit_c, tmb%psit_f, norm)

       call untranspose_localized(iproc, nproc, tmb%npsidim_orbs, tmb%orbs, tmb%collcom, &
            tmb%psit_c, tmb%psit_f, tmb%psi, tmb%lzd)

       i_all = -product(shape(norm))*kind(norm)
       deallocate(norm,stat=i_stat)
       call memocc(i_stat,i_all,'norm',subname)
   end if


          ! Update the kernel
  if (input%lin%scf_mode/=LINEAR_FOE) then
      !tmb%can_use_transposed=.false.
      !overlap_calculated = .false.
      !nullify(tmb%psit_c)
      !nullify(tmb%psit_f)
      call reconstruct_kernel(iproc, nproc, tmb%orthpar%methTransformOverlap, tmb%orthpar%blocksize_pdsyev, &
           tmb%orthpar%blocksize_pdgemm, KSwfn%orbs, tmb, overlap_calculated)
      i_all = -product(shape(tmb%psit_c))*kind(tmb%psit_c)
      deallocate(tmb%psit_c,stat=i_stat)
      call memocc(i_stat,i_all,'tmb%psit_c',subname)
      i_all = -product(shape(tmb%psit_f))*kind(tmb%psit_f)
      deallocate(tmb%psit_f,stat=i_stat)
      call memocc(i_stat,i_all,'tmb%psit_f',subname)
  else
     ! By doing an LCAO input guess
     tmb%can_use_transposed=.false.
     tmb%ham_descr%can_use_transposed=.false.
     ! the following subroutine will overwrite phi, therefore store in a temporary array...
     !!allocate(phi_tmp(size(tmb%psi)), stat=i_stat)
     !!call memocc(i_stat, phi_tmp, 'phi_tmp', subname)
     !!call dcopy(size(tmb%psi), tmb%psi, 1, phi_tmp, 1)
     call inputguessConfinement(iproc, nproc, at, input, KSwfn%Lzd%hgrids(1),KSwfn%Lzd%hgrids(2),KSwfn%Lzd%hgrids(3), &
          rxyz,nlpsp,GPU,KSwfn%orbs, kswfn, tmb,denspot,denspot0,energs)
     !!call dcopy(size(tmb%psi), phi_tmp, 1, tmb%psi, 1)
     !!i_all=-product(shape(phi_tmp))*kind(phi_tmp)
     !!deallocate(phi_tmp, stat=i_stat)
     !!call memocc(i_stat, i_all, 'phi_tmp', subname)
     if(tmb%can_use_transposed) then
         i_all=-product(shape(tmb%psit_c))*kind(tmb%psit_c)
         deallocate(tmb%psit_c, stat=i_stat)
         call memocc(i_stat, i_all, 'tmb%psit_c', subname)
         i_all=-product(shape(tmb%psit_f))*kind(tmb%psit_f)
         deallocate(tmb%psit_f, stat=i_stat)
         call memocc(i_stat, i_all, 'tmb%psit_f', subname)
     end if
  end if

  !!if (iproc==0) then
  !!  do i_stat=1,size(tmb%linmat%denskern%matrix_compr)
  !!    write(*,'(a,i8,es20.10)') 'i_stat, tmb%linmat%denskern%matrix_compr(i_stat)', i_stat, tmb%linmat%denskern%matrix_compr(i_stat)
  !!             denspot%rhov)
  !!  end do
  !!end if

          ! Must initialize rhopotold (FOR NOW... use the trivial one)
  if (input%lin%scf_mode/=LINEAR_FOE) then
      call communicate_basis_for_density_collective(iproc, nproc, tmb%lzd, max(tmb%npsidim_orbs,tmb%npsidim_comp), &
           tmb%orbs, tmb%psi, tmb%collcom_sr)
      call sumrho_for_TMBs(iproc, nproc, KSwfn%Lzd%hgrids(1), KSwfn%Lzd%hgrids(2), KSwfn%Lzd%hgrids(3), &
           tmb%collcom_sr, tmb%linmat%denskern_large, KSwfn%Lzd%Glr%d%n1i*KSwfn%Lzd%Glr%d%n2i*denspot%dpbox%n3d, denspot%rhov)
      call dcopy(max(denspot%dpbox%ndims(1)*denspot%dpbox%ndims(2)*denspot%dpbox%n3p,1)*input%nspin, &
           denspot%rhov(1), 1, denspot0(1), 1)
      call updatePotential(input%ixc,input%nspin,denspot,energs%eh,energs%exc,energs%evxc)
      call local_potential_dimensions(tmb%lzd,tmb%orbs,denspot%dpbox%ngatherarr(0,1))
  end if

  ! Orthonormalize the input guess if necessary
  if (input%experimental_mode .and. input%lin%scf_mode/=LINEAR_FOE) then                 
      if (iproc==0) call yaml_map('orthonormalization of input guess','standard')        
      tmb%can_use_transposed = .false.                                         
      call orthonormalizeLocalized(iproc, nproc, -1, tmb%npsidim_orbs, tmb%orbs, tmb%lzd, &
           tmb%linmat%ovrlp, tmb%linmat%inv_ovrlp_large, &
           tmb%collcom, tmb%orthpar, tmb%psi, tmb%psit_c, tmb%psit_f, tmb%can_use_transposed)
  end if  


END SUBROUTINE input_memory_linear

        subroutine input_wf_disk(iproc, nproc, input_wf_format, d, hx, hy, hz, &
             & in, atoms, rxyz, rxyz_old, wfd, orbs, psi)
          use module_defs
          use module_types
          use module_interfaces, except_this_one => input_wf_disk
          implicit none

          integer, intent(in) :: iproc, nproc, input_wf_format
          type(grid_dimensions), intent(in) :: d
          real(gp), intent(in) :: hx, hy, hz
          type(input_variables), intent(in) :: in
          type(atoms_data), intent(in) :: atoms
  real(gp), dimension(3, atoms%astruct%nat), intent(in) :: rxyz
  real(gp), dimension(3, atoms%astruct%nat), intent(out) :: rxyz_old
          type(wavefunctions_descriptors), intent(in) :: wfd
          type(orbitals_data), intent(inout) :: orbs
          real(wp), dimension(:), pointer :: psi

          integer :: ierr

          !restart from previously calculated wavefunctions, on disk
          !since each processor read only few eigenvalues, initialise them to zero for all
          call to_zero(orbs%norb*orbs%nkpts,orbs%eval(1))

          call readmywaves(iproc,trim(in%dir_output) // "wavefunction", input_wf_format, &
               & orbs,d%n1,d%n2,d%n3,hx,hy,hz,atoms,rxyz_old,rxyz,wfd,psi)

          !reduce the value for all the eigenvectors
          if (nproc > 1) call mpiallred(orbs%eval(1),orbs%norb*orbs%nkpts,MPI_SUM,bigdft_mpi%mpi_comm,ierr)

          if (in%iscf > SCF_KIND_DIRECT_MINIMIZATION) then
             !recalculate orbitals occupation numbers
             call evaltoocc(iproc,nproc,.false.,in%Tel,orbs,in%occopt)
             !read potential depending of the mixing scheme
             !considered as optional in the mixing case
             !inquire(file=trim(in%dir_output)//'local_potential.cube',exist=potential_from_disk)
             !if (potential_from_disk)  then
             !   call read_potential_from_disk(iproc,nproc,trim(in%dir_output)//'local_potential.cube',&
     !        atoms%astruct%geocode,ngatherarr,Lzd%Glr%d%n1i,Lzd%Glr%d%n2i,Lzd%Glr%d%n3i,n3p,in%nspin,hxh,hyh,hzh,rhopot)
             !end if
          end if
        END SUBROUTINE input_wf_disk

        !> Input guess wavefunction diagonalization
        subroutine input_wf_diag(iproc,nproc,at,denspot,&
             orbs,nvirt,comms,Lzd,energs,rxyz,&
             nlpsp,ixc,psi,hpsi,psit,G,&
             nspin,symObj,GPU,input,onlywf,proj_G,paw)
           ! Input wavefunctions are found by a diagonalization in a minimal basis set
           ! Each processors write its initial wavefunctions into the wavefunction file
           ! The files are then read by readwave
           ! @todo pass GPU to be a local variable of this routine (initialized and freed here)
           use module_base
           use module_interfaces, except_this_one => input_wf_diag
           use module_types
   use Poisson_Solver, except_dp => dp, except_gp => gp, except_wp => wp
           use yaml_output
           use gaussians
           implicit none
           !Arguments
           integer, intent(in) :: iproc,nproc,ixc
           integer, intent(inout) :: nspin,nvirt
   logical, intent(in) :: onlywf  !if .true. finds only the WaveFunctions and return
           type(atoms_data), intent(in) :: at
           type(DFT_PSP_projectors), intent(inout) :: nlpsp
           type(local_zone_descriptors), intent(inout) :: Lzd
           type(communications_arrays), intent(in) :: comms
           type(energy_terms), intent(inout) :: energs
           type(orbitals_data), intent(inout) :: orbs
           type(DFT_local_fields), intent(inout) :: denspot
           type(GPU_pointers), intent(in) :: GPU
           type(input_variables), intent(in) :: input
           type(symmetry_data), intent(in) :: symObj
   real(gp), dimension(3,at%astruct%nat), intent(in) :: rxyz
           type(gaussian_basis), intent(out) :: G !basis for davidson IG
           real(wp), dimension(:), pointer :: psi,hpsi,psit
   type(gaussian_basis),dimension(at%astruct%ntypes),optional,intent(in)::proj_G
   type(paw_objects),optional,intent(inout)::paw
           !local variables
           character(len=*), parameter :: subname='input_wf_diag'
           logical :: switchGPUconv,switchOCLconv
   integer :: ii,jj
           integer :: i_stat,i_all,nspin_ig,ncplx,irhotot_add,irho_add,ispin,ikpt
           real(gp) :: hxh,hyh,hzh,etol,accurex,eks
           type(orbitals_data) :: orbse
           type(communications_arrays) :: commse
           integer, dimension(:,:), allocatable :: norbsc_arr
           real(wp), dimension(:), allocatable :: passmat
           !real(wp), dimension(:,:,:), allocatable :: mom_vec
           real(gp), dimension(:), allocatable :: locrad
        !   real(wp), dimension(:), pointer :: pot,pot1
           real(wp), dimension(:,:,:), pointer :: psigau
   real(wp),dimension(:),allocatable::psi_
           type(confpot_data), dimension(:), allocatable :: confdatarr
           type(local_zone_descriptors) :: Lzde
           type(GPU_pointers) :: GPUe
        !!$   integer :: idum=0
        !!$   real(kind=4) :: tt,builtin_rand
!!$   real(wp), dimension(:), allocatable :: ovrlp
!!$   real(wp), dimension(:,:), allocatable :: smat,tmp

        !yk
        !  integer :: i!,iorb,jorb,icplx

           allocate(norbsc_arr(at%natsc+1,nspin+ndebug),stat=i_stat)
           call memocc(i_stat,norbsc_arr,'norbsc_arr',subname)
   allocate(locrad(at%astruct%nat+ndebug),stat=i_stat)
           call memocc(i_stat,locrad,'locrad',subname)

           if (iproc == 0) then
              !yaml_output
              !call yaml_newline()
           end if
           !spin for inputguess orbitals
           if (nspin == 4) then
              nspin_ig=1
           else
              nspin_ig=nspin
           end if

           call inputguess_gaussian_orbitals(iproc,nproc,at,rxyz,nvirt,nspin_ig,&
                orbs,orbse,norbsc_arr,locrad,G,psigau,eks,1)

           !allocate communications arrays for inputguess orbitals
           !call allocate_comms(nproc,orbse,commse,subname)
           call orbitals_communicators(iproc,nproc,Lzd%Glr,orbse,commse,basedist=comms%nvctr_par(0:,1:))  

           !use the eval array of orbse structure to save the original values
           allocate(orbse%eval(orbse%norb*orbse%nkpts+ndebug),stat=i_stat)
           call memocc(i_stat,orbse%eval,'orbse%eval',subname)

           hxh=.5_gp*Lzd%hgrids(1)
           hyh=.5_gp*Lzd%hgrids(2)
           hzh=.5_gp*Lzd%hgrids(3)

           !check the communication distribution
          !call check_communications(iproc,nproc,orbse,Lzd%Glr,commse)

           !once the wavefunction coefficients are known perform a set 
           !of nonblocking send-receive operations to calculate overlap matrices

           !!!  !create mpirequests array for controlling the success of the send-receive operation
           !!!  allocate(mpirequests(nproc-1+ndebug),stat=i_stat)
           !!!  call memocc(i_stat,mpirequests,'mpirequests',subname)
           !!!
           !!!  call nonblocking_transposition(iproc,nproc,G%ncoeff,orbse%isorb+orbse%norbp,&
           !!!       orbse%nspinor,psigau,orbse%norb_par,mpirequests)

        ! ###################################################################
        !!experimental part for building the localisation regions
        ! ###################################################################
           call nullify_local_zone_descriptors(Lzde)
           call create_LzdLIG(iproc,nproc,orbs%nspin,input%linear,&
                Lzd%hgrids(1),Lzd%hgrids(2),Lzd%hgrids(3),Lzd%Glr,at,orbse,rxyz,Lzde)

           if(iproc==0 .and. Lzde%linear) call yaml_comment('Entering the Linear IG')
           !write(*,'(1x,A)') 'Entering the Linear IG'

           ! determine the wavefunction dimension
           call wavefunction_dimension(Lzde,orbse)

           !allocate the wavefunction in the transposed way to avoid allocations/deallocations
             allocate(psi(max(orbse%npsidim_orbs,orbse%npsidim_comp)+ndebug),stat=i_stat)
             call memocc(i_stat,psi,'psi',subname)

             !allocate arrays for the GPU if a card is present
             GPUe = GPU
             switchGPUconv=.false.
             switchOCLconv=.false.
             if (GPUconv) then
                call prepare_gpu_for_locham(Lzde%Glr%d%n1,Lzde%Glr%d%n2,Lzde%Glr%d%n3,nspin_ig,&
                     Lzd%hgrids(1),Lzd%hgrids(2),Lzd%hgrids(3),Lzde%Glr%wfd,orbse,GPUe)
        if (iproc == 0) call yaml_comment('GPU data allocated')
             else if (OCLconv) then
        call allocate_data_OCL(Lzde%Glr%d%n1,Lzde%Glr%d%n2,Lzde%Glr%d%n3,at%astruct%geocode,&
                     nspin_ig,Lzde%Glr%wfd,orbse,GPUe)
        if (iproc == 0) call yaml_comment('GPU data allocated')
        !if (iproc == 0) write(*,*) 'GPU data allocated'
             end if

            call timing(iproc,'wavefunction  ','ON')   
           !use only the part of the arrays for building the hamiltonian matrix
             call gaussians_to_wavelets_new(iproc,nproc,Lzde,orbse,G,&
                  psigau(1,1,min(orbse%isorb+1,orbse%norb)),psi)
            call timing(iproc,'wavefunction  ','OF')
             i_all=-product(shape(locrad))*kind(locrad)
             deallocate(locrad,stat=i_stat)
             call memocc(i_stat,i_all,'locrad',subname)

! IF onlywf return
  if(onlywf) then

     !for testing
     !application of the hamiltonian for gaussian based treatment
     !orbse%nspin=nspin
     !call sumrho(iproc,nproc,orbse,Lzd,hxh,hyh,hzh,denspot%dpcom%nscatterarr,&
     !GPU,symObj,denspot%rhod,psi,denspot%rho_psi)
     !call communicate_density(iproc,nproc,orbse%nspin,hxh,hyh,hzh,Lzd,&
     !   denspot%rhod,denspot%dpcom%nscatterarr,denspot%rho_psi,denspot%rhov)
     !orbse%nspin=nspin_ig
     !   
     !!-- if spectra calculation uses a energy dependent potential
     !!    input_wf_diag will write (to be used in abscalc)
     !!    the density to the file electronic_density.cube
     !!  The writing is activated if  5th bit of  in%potshortcut is on.
     !   call plot_density_cube_old('electronic_density',&
     !        iproc,nproc,Lzd%Glr%d%n1,Lzd%Glr%d%n2,Lzd%Glr%d%n3,&
     !        Lzd%Glr%d%n1i,Lzd%Glr%d%n2i,Lzd%Glr%d%n3i,denspot%dpcom%nscatterarr(iproc,2),&
     !        nspin,hxh,hyh,hzh,at,rxyz,denspot%dpcom%ngatherarr,&
     !        denspot%rhov(1+denspot%dpcom%nscatterarr(iproc,4)*Lzd%Glr%d%n1i*Lzd%Glr%d%n2i))
     !---
    !reallocate psi, with good dimensions:
    ii=max(1,max(orbse%npsidim_orbs,orbse%npsidim_comp))+ndebug
    jj=max(1,max(orbs%npsidim_orbs,orbs%npsidim_comp))+ndebug
    if(ii .ne. jj) then
      allocate(psi_(jj),stat=i_stat)
      call memocc(i_stat,psi_,'psi_',subname)
      if(jj<=ii) psi_=psi(1:jj)
      if(jj>ii) then
        psi_(1:ii)=psi(1:ii)
        psi_(ii+1:jj)=1.0d0
      end if
      i_all=-product(shape(psi))*kind(psi)
      deallocate(psi,stat=i_stat)
      call memocc(i_stat,i_all,'psi',subname)
      allocate(psi(jj),stat=i_stat)
      call memocc(i_stat,psi,'psi',subname)
      psi=psi_
      i_all=-product(shape(psi_))*kind(psi_)
      deallocate(psi_,stat=i_stat)
      call memocc(i_stat,i_all,'psi_',subname)
    end if


    !allocate the wavefunction in the transposed way to avoid allocations/deallocations
    allocate(hpsi(max(1,max(orbs%npsidim_orbs,orbs%npsidim_comp))+ndebug),stat=i_stat)
    call memocc(i_stat,hpsi,'hpsi',subname)
     
    if(present(paw)) then
      allocate(paw%spsi(max(1,max(orbs%npsidim_orbs,orbs%npsidim_comp))+ndebug),stat=i_stat)
      call memocc(i_stat,paw%spsi,'spsi',subname)
    end if

    !The following lines are copied from LDiagHam:
    nullify(psit)
    !
    !in the case of minimal basis allocate now the transposed wavefunction
    !otherwise do it only in parallel
    if ( nproc > 1) then
       allocate(psit(max(orbs%npsidim_orbs,orbs%npsidim_comp)+ndebug),stat=i_stat)
       call memocc(i_stat,psit,'psit',subname)
    else
       psit => hpsi
    end if

    !transpose the psi wavefunction
    call transpose_v2(iproc,nproc,orbs,Lzd,comms,psi,work=hpsi,outadd=psit)

    nullify(G%rxyz)

    !Set orbs%eval=-0.5.
    !This will be done in LDiagHam
    !For the moment we skip this, since hpsi is not yet calculated
    !(hpsi is an input argument in LDiagHam)
    orbs%eval(:)=-0.5_wp

    call deallocate_input_wfs()
    return 
  end if

           !check the size of the rhopot array related to NK SIC
        !!$   nrhodim=nspin
        !!$   i3rho_add=0
        !!$   if (input%SIC%approach=='NK') then
        !!$      nrhodim=2*nrhodim
        !!$     i3rho_add=Lzd%Glr%d%n1i*Lzd%Glr%d%n2i*nscatterarr(iproc,4)+1
        !!$   end if

           !application of the hamiltonian for gaussian based treatment
           !if(.false.) then
           !   call sumrho(iproc,nproc,orbse,Lzd%Glr,hxh,hyh,hzh,psi,rhopot,&
           !        nscatterarr,nspin,GPU,symObj,irrzon,phnons,rhodsc)
           !end if

          ! test merging of the cubic and linear code
          !call sumrhoLinear(iproc,nproc,Lzd,orbse,hxh,hyh,hzh,psi,rhopot,nscatterarr,nspin,GPU,symObj, irrzon, phnons, rhodsc)    

           !spin adaptation for the IG in the spinorial case
           orbse%nspin=nspin
           call sumrho(denspot%dpbox,orbse,Lzde,GPUe,symObj,denspot%rhod,psi,denspot%rho_psi)
           call communicate_density(denspot%dpbox,orbse%nspin,denspot%rhod,denspot%rho_psi,denspot%rhov,.false.)
           call denspot_set_rhov_status(denspot, ELECTRONIC_DENSITY, 0, iproc, nproc)

           orbse%nspin=nspin_ig

           !before creating the potential, save the density in the second part 
           !if the case of NK SIC, so that the potential can be created afterwards
           !copy the density contiguously since the GGA is calculated inside the NK routines
           if (input%SIC%approach=='NK') then
              irhotot_add=Lzde%Glr%d%n1i*Lzde%Glr%d%n2i*denspot%dpbox%nscatterarr(iproc,4)+1
              irho_add=Lzde%Glr%d%n1i*Lzde%Glr%d%n2i*denspot%dpbox%nscatterarr(iproc,1)*input%nspin+1
              do ispin=1,input%nspin
                call dcopy(Lzde%Glr%d%n1i*Lzde%Glr%d%n2i*denspot%dpbox%nscatterarr(iproc,2),&
                     denspot%rhov(irhotot_add),1,denspot%rhov(irho_add),1)
                irhotot_add=irhotot_add+Lzde%Glr%d%n1i*Lzde%Glr%d%n2i*denspot%dpbox%nscatterarr(iproc,1)
                irho_add=irho_add+Lzde%Glr%d%n1i*Lzde%Glr%d%n2i*denspot%dpbox%nscatterarr(iproc,2)
              end do
           end if
           call updatePotential(ixc,nspin,denspot,energs%eh,energs%exc,energs%evxc)

!!$   !experimental
!!$   if (nproc == 1) then
!!$
!!$
        !!$     !calculate the overlap matrix as well as the kinetic overlap
        !!$     !in view of complete gaussian calculation
        !!$     allocate(ovrlp(G%ncoeff*G%ncoeff),stat=i_stat)
        !!$     call memocc(i_stat,ovrlp,'ovrlp',subname)
        !!$     allocate(tmp(G%ncoeff,orbse%norb),stat=i_stat)
        !!$     call memocc(i_stat,tmp,'tmp',subname)
        !!$     allocate(smat(orbse%norb,orbse%norb),stat=i_stat)
        !!$     call memocc(i_stat,smat,'smat',subname)
        !!$
        !!$     !overlap calculation of the gaussian matrix
        !!$     call gaussian_overlap(G,G,ovrlp)
        !!$     call dsymm('L','U',G%ncoeff,orbse%norb,1.0_gp,ovrlp(1),G%ncoeff,&
!!$          psigau(1,1,1),G%ncoeff,0.d0,tmp(1,1),G%ncoeff)
        !!$
        !!$     call gemm('T','N',orbse%norb,orbse%norb,G%ncoeff,1.0_gp,&
!!$          psigau(1,1,1),G%ncoeff,tmp(1,1),G%ncoeff,0.0_wp,smat(1,1),orbse%norb)
!!$
        !!$     !print overlap matrices
!!$     print *,'OVERLAP' 
        !!$     do i=1,orbse%norb
!!$        write(*,'(i4,30(1pe10.2))')i,(smat(i,iorb),iorb=1,orbse%norb)
!!$        !write(*,'(i4,30(1pe10.2))')i,(ovrlp(i+(iorb-1)*orbse%norb),&
!!$        !     iorb=1,orbse%norb)
        !!$     end do
!!$     
!!$     !overlap calculation of the kinetic operator
!!$     call kinetic_overlap(G,G,ovrlp)
!!$     call dsymm('L','U',G%ncoeff,orbse%norb,1.0_gp,ovrlp(1),G%ncoeff,&
!!$          psigau(1,1,1),G%ncoeff,0.d0,tmp(1,1),G%ncoeff)
!!$
!!$     call gemm('T','N',orbse%norb,orbse%norb,G%ncoeff,1.0_gp,&
!!$          psigau(1,1,1),G%ncoeff,tmp(1,1),G%ncoeff,0.0_wp,smat(1,1),orbse%norb)
!!$
!!$     !print overlap matrices
!!$     print *,'HAMILTONIAN' 
!!$     tt=0.0_wp
!!$     do i=1,orbse%norb
!!$        write(*,'(i4,30(1pe10.2))')i,(smat(i,iorb),iorb=1,orbse%norb)
!!$        !write(12,'(i5,30(1pe15.8))')i,(smat(i,iorb),iorb=1,orbse%norb)
!!$        tt=tt+smat(i,i)
!!$     end do
!!$     print *,'trace',tt
!!$stop

           !!!
           !!!     !overlap calculation of the kinetic operator
           !!!     call cpu_time(t0)
           !!!     call potential_overlap(G,G,rhopot,Glr%d%n1i,Glr%d%n2i,Glr%d%n3i,hxh,hyh,hzh,&
           !!!          ovrlp)
           !!!     call cpu_time(t1)
           !!!     call dsymm('L','U',G%ncoeff,orbse%norb,1.0_gp,ovrlp(1),G%ncoeff,&
   !!!          psigau(1,1),G%ncoeff,0.d0,tmp(1,1),G%ncoeff)
           !!!
           !!!     call gemm('T','N',orbse%norb,orbse%norb,G%ncoeff,1.0_gp,&
   !!!          psigau(1,1),G%ncoeff,tmp(1,1),G%ncoeff,0.0_wp,smat(1,1),orbse%norb)
           !!!
           !!!     !print overlap matrices
           !!!     tt=0.0_wp
           !!!     do i=1,orbse%norb
           !!!        write(*,'(i5,30(1pe15.8))')i,(smat(i,iorb),iorb=1,orbse%norb)
           !!!        !write(12,'(i5,30(1pe15.8))')i,(smat(i,iorb),iorb=1,orbse%norb)
           !!!        tt=tt+smat(i,i)
           !!!     end do
           !!!     print *,'trace',tt
           !!!     print *, 'time',t1-t0
           !!!
!!$     i_all=-product(shape(ovrlp))*kind(ovrlp)
!!$     deallocate(ovrlp,stat=i_stat)
!!$     call memocc(i_stat,i_all,'ovrlp',subname)
!!$     i_all=-product(shape(tmp))*kind(tmp)
!!$     deallocate(tmp,stat=i_stat)
!!$     call memocc(i_stat,i_all,'tmp',subname)
!!$     i_all=-product(shape(smat))*kind(smat)
!!$     deallocate(smat,stat=i_stat)
!!$     call memocc(i_stat,i_all,'smat',subname)
!!$  end if
           
           
           !allocate the wavefunction in the transposed way to avoid allocations/deallocations
           allocate(hpsi(max(1,max(orbse%npsidim_orbs,orbse%npsidim_comp))+ndebug),stat=i_stat)
           call memocc(i_stat,hpsi,'hpsi',subname)
           
             !call dcopy(orbse%npsidim,psi,1,hpsi,1)
           if (input%exctxpar == 'OP2P') then
              energs%eexctX = UNINITIALIZED(1.0_gp)
           else
              energs%eexctX=0.0_gp
           end if
           
           !change temporarily value of Lzd%npotddim
           allocate(confdatarr(orbse%norbp)) !no stat so tho make it crash
           call local_potential_dimensions(Lzde,orbse,denspot%dpbox%ngatherarr(0,1))
        !   print *,'here',iproc   
           call default_confinement_data(confdatarr,orbse%norbp)

           !spin adaptation for the IG in the spinorial case
           orbse%nspin=nspin
           call full_local_potential(iproc,nproc,orbse,Lzde,Lzde%lintyp,denspot%dpbox,denspot%rhov,denspot%pot_work)
           orbse%nspin=nspin_ig

           !write(*,*) 'size(denspot%pot_work)', size(denspot%pot_work)
           call FullHamiltonianApplication(iproc,nproc,at,orbse,rxyz,&
                Lzde,nlpsp,confdatarr,denspot%dpbox%ngatherarr,denspot%pot_work,psi,hpsi,&
                energs,input%SIC,GPUe,&
                pkernel=denspot%pkernelseq)
!!$   if (orbse%npsidim_orbs > 0) call to_zero(orbse%npsidim_orbs,hpsi(1))
!!$   call  LocalHamiltonianApplication(iproc,nproc,at,orbse,&
!!$        Lzde,confdatarr,denspot%dpbox%ngatherarr,denspot%pot_work,psi,hpsi,&
!!$        energs,input%SIC,GPUe,3,pkernel=denspot%pkernelseq)

           call denspot_set_rhov_status(denspot, KS_POTENTIAL, 0, iproc, nproc)
            !restore the good value
            call local_potential_dimensions(Lzde,orbs,denspot%dpbox%ngatherarr(0,1))

             !deallocate potential
             call free_full_potential(denspot%dpbox%mpi_env%nproc,Lzde%lintyp,denspot%pot_work,subname)

             i_all=-product(shape(orbse%ispot))*kind(orbse%ispot)
             deallocate(orbse%ispot,stat=i_stat)
             call memocc(i_stat,i_all,'orbse%ispot',subname)

             deallocate(confdatarr)
         
           !!!  !calculate the overlap matrix knowing that the original functions are gaussian-based
           !!!  allocate(thetaphi(2,G%nat+ndebug),stat=i_stat)
           !!!  call memocc(i_stat,thetaphi,'thetaphi',subname)
           !!!  thetaphi=0.0_gp
           !!!
           !!!  !calculate the scalar product between the hamiltonian and the gaussian basis
           !!!  allocate(hpsigau(G%ncoeff,orbse%norbp+ndebug),stat=i_stat)
           !!!  call memocc(i_stat,hpsigau,'hpsigau',subname)
           !!!
           !!!
   !!!  call wavelets_to_gaussians(at%astruct%geocode,orbse%norbp,Glr%d%n1,Glr%d%n2,Glr%d%n3,G,&
           !!!       thetaphi,hx,hy,hz,Glr%wfd,hpsi,hpsigau)
           !!!
           !!!  i_all=-product(shape(thetaphi))*kind(thetaphi)
           !!!  deallocate(thetaphi,stat=i_stat)
           !!!  call memocc(i_stat,i_all,'thetaphi',subname)
           
   accurex=abs(eks-energs%ekin)
   !tolerance for comparing the eigenvalues in the case of degeneracies
   etol=accurex/real(orbse%norbu,gp)

   !if (iproc == 0 .and. verbose > 1 .and. at%astruct%geocode=='F') write(*,'(1x,a,2(f19.10))') 'done. ekin_sum,eks:',energs%ekin,eks
   if (iproc == 0 .and. verbose > 1 .and. at%astruct%geocode=='F') call yaml_map('Expected kinetic energy',eks,fmt='(f19.10)')
   if (iproc==0) call yaml_newline()

   call total_energies(energs, 0, iproc)

           if (iproc==0) then
              !yaml output
              !call write_energies(0,0,energs,0.0_gp,0.0_gp,'Input Guess')
              call write_energies(0,0,energs,0.0_gp,0.0_gp,'')
   endif
          
           !!!  call Gaussian_DiagHam(iproc,nproc,at%natsc,nspin,orbs,G,mpirequests,&
           !!!       psigau,hpsigau,orbse,etol,norbsc_arr)
           
           
           !!!  i_all=-product(shape(mpirequests))*kind(mpirequests)
           !!!  deallocate(mpirequests,stat=i_stat)
           !!!  call memocc(i_stat,i_all,'mpirequests',subname)
           
           !!!  i_all=-product(shape(hpsigau))*kind(hpsigau)
           !!!  deallocate(hpsigau,stat=i_stat)
           !!!  call memocc(i_stat,i_all,'hpsigau',subname)
           
             !free GPU if it is the case
             if (GPUconv) then
                call free_gpu(GPUe,orbse%norbp)
        if (iproc == 0) call yaml_comment('GPU data deallocated')
             else if (OCLconv) then
                call free_gpu_OCL(GPUe,orbse,nspin_ig)
        if (iproc == 0) call yaml_comment('GPU data deallocated')
             end if

             !if (iproc == 0 .and. verbose > 1) write(*,'(1x,a)')&
             !     'Input Wavefunctions Orthogonalization:'
          
             !nullify psit (will be created in DiagHam)
             nullify(psit)

             !psivirt can be eliminated here, since it will be allocated before davidson
             !with a gaussian basis
           !!$  call DiagHam(iproc,nproc,at%natsc,nspin_ig,orbs,Glr%wfd,comms,&
           !!$       psi,hpsi,psit,orbse,commse,etol,norbsc_arr,orbsv,psivirt)

            !allocate the passage matrix for transforming the LCAO wavefunctions in the IG wavefucntions
             ncplx=1
             if (orbs%nspinor > 1) ncplx=2
             allocate(passmat(ncplx*orbs%nkptsp*(orbse%norbu*orbs%norbu+orbse%norbd*orbs%norbd)+ndebug),stat=i_stat)
             call memocc(i_stat,passmat,'passmat',subname)
          !!print '(a,10i5)','iproc,passmat',iproc,ncplx*orbs%nkptsp*(orbse%norbu*orbs%norbu+orbse%norbd*orbs%norbd),&
          !!     orbs%nspinor,orbs%nkptsp,orbse%norbu,orbse%norbd,orbs%norbu,orbs%norbd

            if (iproc==0) call yaml_newline()

           !test merging of Linear and cubic
             call LDiagHam(iproc,nproc,at%natsc,nspin_ig,orbs,Lzd,Lzde,comms,&
                 psi,hpsi,psit,input%orthpar,passmat,input%iscf,input%Tel,input%occopt,&
                 orbse,commse,etol,norbsc_arr)

             i_all=-product(shape(passmat))*kind(passmat)
             deallocate(passmat,stat=i_stat)
             call memocc(i_stat,i_all,'passmat',subname)

           if (input%iscf > SCF_KIND_DIRECT_MINIMIZATION .or. input%Tel > 0.0_gp) then
              
              !restore the occupations as they are extracted from DiagHam
              !use correct copying due to k-points
              do ikpt=1,orbs%nkpts
                 call vcopy(orbs%norbu,orbse%occup((ikpt-1)*orbse%norb+1),1,&
                      orbs%occup((ikpt-1)*orbs%norb+1),1)
                 if (orbs%norbd > 0) then
                    call vcopy(orbs%norbd,orbse%occup((ikpt-1)*orbse%norb+orbse%norbu+1),1,&
                         orbs%occup((ikpt-1)*orbs%norb+orbs%norbu+1),1)
                 end if
              end do
              !call dcopy(orbs%norb*orbs%nkpts,orbse%occup(1),1,orbs%occup(1),1) !this is not good with k-points
              !associate the entropic energy contribution
              orbs%eTS=orbse%eTS
              
           end if

        !!$   !yaml output
        !!$   if (iproc ==0) then
        !!$      if(orbse%nspinor==4) then
        !!$         allocate(mom_vec(4,orbse%norb,min(nproc,2)+ndebug),stat=i_stat)
        !!$         call memocc(i_stat,mom_vec,'mom_vec',subname)
        !!$         call to_zero(4*orbse%norb*min(nproc,2),mom_vec(1,1,1))
        !!$      end if
        !!$
        !!$      !experimental part to show the actual occupation numbers which will be put in the inputguess
           !!put the occupation numbers of the normal orbitals
           !call vcopy(orbs%norb*orbs%nkpts,orbs%occup(1),1,orbse%occup(1),1)
           !!put to zero the other values
   !call to_zero(orbse%norb*orbse%nkpts-orbs%norb*orbs%nkpts,&
   !     orbse%occup(min(orbse%norb*orbse%nkpts,orbs%norb*orbs%nkpts+1)))
!!$
!!$      call write_eigenvalues_data(orbse,mom_vec)
!!$      yaml_indent=yaml_indent-2
!!$
!!$      if (orbs%nspinor ==4) then
!!$         i_all=-product(shape(mom_vec))*kind(mom_vec)
!!$         deallocate(mom_vec,stat=i_stat)
!!$         call memocc(i_stat,i_all,'mom_vec',subname)
!!$      end if
!!$   end if


   call deallocate_input_wfs()

contains

subroutine deallocate_input_wfs()

   call deallocate_comms(commse,subname)

   i_all=-product(shape(norbsc_arr))*kind(norbsc_arr)
   deallocate(norbsc_arr,stat=i_stat)
   call memocc(i_stat,i_all,'norbsc_arr',subname)

   if (iproc == 0) then
      !gaussian estimation valid only for Free BC
      if (at%astruct%geocode == 'F') then
         call yaml_newline()
         call yaml_open_map('Accuracy estimation for this run')
         call yaml_map('Energy',accurex,fmt='(1pe9.2)')
         call yaml_map('Convergence Criterion',accurex/real(orbs%norb,kind=8),fmt='(1pe9.2)')
         call yaml_close_map()
            !write(*,'(1x,a,1pe9.2)') 'expected accuracy in energy ',accurex
      !write(*,'(1x,a,1pe9.2)') &
      !&   'expected accuracy in energy per orbital ',accurex/real(orbs%norb,kind=8)
         !write(*,'(1x,a,1pe9.2)') &
         !     'suggested value for gnrm_cv ',accurex/real(orbs%norb,kind=8)
      end if
   endif

   !here we can define the subroutine which generates the coefficients for the virtual orbitals
   call deallocate_gwf(G,subname)
   call deallocate_local_zone_descriptors(Lzde, subname)

   i_all=-product(shape(psigau))*kind(psigau)
   deallocate(psigau,stat=i_stat)
   call memocc(i_stat,i_all,'psigau',subname)

   call deallocate_orbs(orbse,subname)
   i_all=-product(shape(orbse%eval))*kind(orbse%eval)
   deallocate(orbse%eval,stat=i_stat)
   call memocc(i_stat,i_all,'orbse%eval',subname)

end subroutine deallocate_input_wfs

END SUBROUTINE input_wf_diag


subroutine input_wf(iproc,nproc,in,GPU,atoms,rxyz,&
     denspot,denspot0,nlpsp,KSwfn,tmb,energs,inputpsi,input_wf_format,norbv,&
     lzd_old,wfd_old,psi_old,d_old,hx_old,hy_old,hz_old,rxyz_old,tmb_old,ref_frags,cdft,&
     locregcenters)
  use module_defs
  use module_types
  use module_interfaces, except_this_one => input_wf
  use module_fragments
  use constrained_dft
  use dynamic_memory
  use yaml_output
  use gaussians, only:gaussian_basis
  implicit none

  integer, intent(in) :: iproc, nproc, inputpsi, input_wf_format
  type(input_variables), intent(in) :: in
  type(GPU_pointers), intent(inout) :: GPU
  real(gp), intent(in) :: hx_old,hy_old,hz_old
  type(atoms_data), intent(inout) :: atoms
  real(gp), dimension(3, atoms%astruct%nat), target, intent(in) :: rxyz
  type(DFT_local_fields), intent(inout) :: denspot
  type(DFT_wavefunction), intent(inout) :: KSwfn,tmb,tmb_old !<input wavefunctions
  real(gp), dimension(*), intent(out) :: denspot0 !< Initial density / potential, if needed
  type(energy_terms), intent(inout) :: energs !<energies of the system
  !real(wp), dimension(:), pointer :: psi,hpsi,psit
  real(wp), dimension(:), pointer :: psi_old
  integer, intent(out) :: norbv
  type(DFT_PSP_projectors), intent(inout) :: nlpsp
  !type(gaussian_basis), intent(inout) :: gbd
  !real(wp), dimension(:,:), pointer :: gaucoeffs
  type(grid_dimensions), intent(in) :: d_old
  real(gp), dimension(3, atoms%astruct%nat), intent(inout) :: rxyz_old
  type(local_zone_descriptors),intent(inout):: lzd_old
  type(wavefunctions_descriptors), intent(inout) :: wfd_old
  type(system_fragment), dimension(:), pointer :: ref_frags
  type(cdft_data), intent(out) :: cdft
  real(kind=8),dimension(3,atoms%astruct%nat),intent(in),optional :: locregcenters
  !local variables
  character(len = *), parameter :: subname = "input_wf"
  integer :: i_stat, nspin, i_all, iat
  type(gaussian_basis) :: Gvirt
  real(wp), allocatable, dimension(:) :: norm
  !wvl+PAW objects
  integer :: iatyp
  type(gaussian_basis),dimension(atoms%astruct%ntypes)::proj_G
  type(paw_objects)::paw
  logical :: overlap_calculated, perx,pery,perz
  real(gp) :: tx,ty,tz,displ,mindist
  real(gp), dimension(:), pointer :: in_frag_charge
  integer :: infoCoeff, iorb, nstates_max
  real(kind=8) :: pnrm
  type(sparseMatrix) :: ham_small
  !!real(gp), dimension(:,:), allocatable :: ks, ksk
  !!real(gp) :: nonidem

  call f_routine(id='input_wf')

  !nullify paw objects:
  do iatyp=1,atoms%astruct%ntypes
  call nullify_gaussian_basis(proj_G(iatyp))
  end do
  paw%usepaw=0 !Not using PAW
  call nullify_paw_objects(paw)



 !determine the orthogonality parameters
  KSwfn%orthpar = in%orthpar
  if (inputpsi == INPUT_PSI_LINEAR_AO .or. inputpsi == INPUT_PSI_DISK_LINEAR &
      .or. inputpsi == INPUT_PSI_MEMORY_LINEAR) then
     tmb%orthpar%methTransformOverlap = in%lin%methTransformOverlap
     tmb%orthpar%nItOrtho = 1
     tmb%orthpar%blocksize_pdsyev = in%lin%blocksize_pdsyev
     tmb%orthpar%blocksize_pdgemm = in%lin%blocksize_pdgemm
     tmb%orthpar%nproc_pdsyev = in%lin%nproc_pdsyev
  end if

  !SIC parameters
  KSwfn%SIC=in%SIC
  !exact exchange parallelization parameter
  KSwfn%exctxpar=in%exctxpar

  !avoid allocation of the eigenvalues array in case of restart
  if ( inputpsi /= INPUT_PSI_MEMORY_WVL .and. &
     & inputpsi /= INPUT_PSI_MEMORY_GAUSS .and. &
       & inputpsi /= INPUT_PSI_MEMORY_LINEAR .and. &
       & inputpsi /= INPUT_PSI_DISK_LINEAR) then
     allocate(KSwfn%orbs%eval(KSwfn%orbs%norb*KSwfn%orbs%nkpts+ndebug),stat=i_stat)
     call memocc(i_stat,KSwfn%orbs%eval,'eval',subname)
  end if
  ! Still do it for linear restart, to be check...
  if (inputpsi == INPUT_PSI_DISK_LINEAR) then
     if(iproc==0) call yaml_comment('ALLOCATING KSwfn%orbs%eval... is this correct?')
     allocate(KSwfn%orbs%eval(KSwfn%orbs%norb*KSwfn%orbs%nkpts+ndebug),stat=i_stat)
     call memocc(i_stat,KSwfn%orbs%eval,'eval',subname)
  end if

  !all the input formats need to allocate psi except the LCAO input_guess
  ! WARNING: at the moment the linear scaling version allocates psi in the same
  ! way as the LCAO input guess, so it is not necessary to allocate it here.
  ! Maybe to be changed later.
  !if (inputpsi /= 0) then

  if (inputpsi /= INPUT_PSI_LCAO .and. inputpsi /= INPUT_PSI_LINEAR_AO .and. inputpsi /= INPUT_PSI_DISK_LINEAR &
     .and. inputpsi /= INPUT_PSI_MEMORY_LINEAR) then
     allocate(KSwfn%psi(max(KSwfn%orbs%npsidim_comp,KSwfn%orbs%npsidim_orbs)+ndebug),stat=i_stat)
     call memocc(i_stat,KSwfn%psi,'psi',subname)
  end if
  if (inputpsi == INPUT_PSI_LINEAR_AO .or. inputpsi == INPUT_PSI_DISK_LINEAR &
      .or. inputpsi == INPUT_PSI_MEMORY_LINEAR) then
     allocate(tmb%psi(max(tmb%npsidim_comp,tmb%npsidim_orbs)), stat=i_stat)
     call memocc(i_stat, tmb%psi, 'tmb%psi', subname)
     !allocate(tmb%confdatarr(tmb%orbs%norbp))
     !call define_confinement_data(tmb%confdatarr,tmb%orbs,rxyz,atoms,&
     !     KSwfn%Lzd%hgrids(1),KSwfn%Lzd%hgrids(2),KSwfn%Lzd%hgrids(3),4,&
     !     in%lin%potentialprefac_lowaccuracy,tmb%lzd,tmb%orbs%onwhichatom)
  else
     allocate(KSwfn%confdatarr(KSwfn%orbs%norbp))
     call default_confinement_data(KSwfn%confdatarr,KSwfn%orbs%norbp)
     call local_potential_dimensions(KSwfn%Lzd,KSwfn%orbs,denspot%dpbox%ngatherarr(0,1))
  end if

  norbv=abs(in%norbv)

  ! INPUT WAVEFUNCTIONS, added also random input guess
  select case(inputpsi)

  case(INPUT_PSI_EMPTY)
     if (iproc == 0) then
        !write( *,'(1x,a)')&
        !     &   '------------------------------------------------- Empty wavefunctions initialization'
        call yaml_comment('Empty wavefunctions initialization',hfill='-')
        call yaml_open_map("Input Hamiltonian")
     end if

     call input_wf_empty(iproc, nproc,KSwfn%psi, KSwfn%hpsi, KSwfn%psit, KSwfn%orbs, &
          in%band_structure_filename, in%nspin, atoms, KSwfn%Lzd%Glr%d, denspot)

  case(INPUT_PSI_RANDOM)
     if (iproc == 0) then
        !write( *,'(1x,a)')&
        !     &   '------------------------------------------------ Random wavefunctions initialization'
        call yaml_comment('Random wavefunctions Initialization',hfill='-')
        call yaml_open_map("Input Hamiltonian")
     end if

     call input_wf_random(KSwfn%psi, KSwfn%orbs)

  case(INPUT_PSI_CP2K)
     if (iproc == 0) then
        !write(*,'(1x,a)')&
        !     &   '--------------------------------------------------------- Import Gaussians from CP2K'
        call yaml_comment('Import Gaussians from CP2K',hfill='-')
        call yaml_open_map("Input Hamiltonian")
     end if

     call input_wf_cp2k(iproc, nproc, in%nspin, atoms, rxyz, KSwfn%Lzd, &
          KSwfn%psi,KSwfn%orbs)

  case(INPUT_PSI_LCAO)
     if (iproc == 0) then
        !write(*,'(1x,a)')&
        !     &   '------------------------------------------------------- Input Wavefunctions Creation'
        call yaml_comment('Wavefunctions from PSP Atomic Orbitals Initialization',hfill='-')
        call yaml_open_map('Input Hamiltonian')
     end if
     nspin=in%nspin
     !calculate input guess from diagonalisation of LCAO basis (written in wavelets)
     call input_wf_diag(iproc,nproc, atoms,denspot,&
          KSwfn%orbs,norbv,KSwfn%comms,KSwfn%Lzd,energs,rxyz,&
          nlpsp,in%ixc,KSwfn%psi,KSwfn%hpsi,KSwfn%psit,&
          Gvirt,nspin,atoms%astruct%sym,GPU,in,.false.)

  case(INPUT_PSI_MEMORY_WVL)
     !restart from previously calculated wavefunctions, in memory
     if (iproc == 0) then
        !write( *,'(1x,a)')&
        !     &   '-------------------------------------------------------------- Wavefunctions Restart'
        call yaml_comment('Wavefunctions Restart',hfill='-')
        call yaml_open_map("Input Hamiltonian")
     end if
     perx=(atoms%astruct%geocode /= 'F')
     pery=(atoms%astruct%geocode == 'P')
     perz=(atoms%astruct%geocode /= 'F')
 
     tx=0.0_gp
     ty=0.0_gp
     tz=0.0_gp
 
     do iat=1,atoms%astruct%nat
        tx=tx+mindist(perx,atoms%astruct%cell_dim(1),rxyz(1,iat),rxyz_old(1,iat))**2
        ty=ty+mindist(pery,atoms%astruct%cell_dim(2),rxyz(2,iat),rxyz_old(2,iat))**2
        tz=tz+mindist(perz,atoms%astruct%cell_dim(3),rxyz(3,iat),rxyz_old(3,iat))**2
     enddo
     displ=sqrt(tx+ty+tz)
 
     if(displ.eq.0d0 .or. in%inguess_geopt == 0) then
        if (in%wfn_history <= 2) then
           call timing(iproc,'restart_wvl   ','ON')
           call input_wf_memory(iproc, atoms, &
                rxyz_old, hx_old, hy_old, hz_old, d_old, wfd_old, psi_old, &
                rxyz,KSwfn%Lzd%hgrids(1),KSwfn%Lzd%hgrids(2),KSwfn%Lzd%hgrids(3),&
                KSwfn%Lzd%Glr%d,KSwfn%Lzd%Glr%wfd,KSwfn%psi, KSwfn%orbs)
              call timing(iproc,'restart_wvl   ','OF')
        else
           call input_wf_memory_history(iproc,KSwfn%orbs,atoms,in%wfn_history,&
                Kswfn%istep_history,KSwfn%oldpsis,rxyz,Kswfn%Lzd,KSwfn%psi)
        end if
     else if(in%inguess_geopt == 1) then
        call timing(iproc,'restart_rsp   ','ON')
        call input_wf_memory_new(nproc, iproc, atoms, &
             rxyz_old, hx_old, hy_old, hz_old, d_old, wfd_old, psi_old,lzd_old, &
             rxyz,KSwfn%Lzd%hgrids(1),KSwfn%Lzd%hgrids(2),KSwfn%Lzd%hgrids(3),&
             KSwfn%Lzd%Glr%d,KSwfn%Lzd%Glr%wfd,KSwfn%psi, KSwfn%orbs,KSwfn%lzd,displ)
        call timing(iproc,'restart_rsp   ','OF')
    else
        stop 'Wrong value of inguess_geopt in input.perf'
    end if

     if (in%iscf > SCF_KIND_DIRECT_MINIMIZATION) &
           call evaltoocc(iproc,nproc,.false.,in%Tel,KSwfn%orbs,in%occopt)
  case(INPUT_PSI_MEMORY_LINEAR)
     if (iproc == 0) then
        call yaml_comment('Support functions Restart',hfill='-')
        call yaml_open_map("Input Hamiltonian")
     end if
      call input_memory_linear(iproc, nproc, atoms, KSwfn, tmb, tmb_old, denspot, in, &
           rxyz_old, rxyz, denspot0, energs, nlpsp, GPU, ref_frags)
  case(INPUT_PSI_DISK_WVL)
     if (iproc == 0) then
        !write( *,'(1x,a)')&
        !     &   '---------------------------------------------------- Reading Wavefunctions from disk'
        call yaml_comment('Reading Wavefunctions from disk',hfill='-')
        call yaml_open_map("Input Hamiltonian")
     end if
     call input_wf_disk(iproc, nproc, input_wf_format, KSwfn%Lzd%Glr%d,&
          KSwfn%Lzd%hgrids(1),KSwfn%Lzd%hgrids(2),KSwfn%Lzd%hgrids(3),&
          in, atoms, rxyz, rxyz_old, KSwfn%Lzd%Glr%wfd, KSwfn%orbs, KSwfn%psi)

  case(INPUT_PSI_MEMORY_GAUSS)
     !restart from previously calculated gaussian coefficients
     if (iproc == 0) then
        !write( *,'(1x,a)')&
        !     &   '--------------------------------------- Quick Wavefunctions Restart (Gaussian basis)'
        call yaml_comment('Quick Wavefunctions Restart (Gaussian basis)',hfill='-')
        call yaml_open_map("Input Hamiltonian")
     end if
     call restart_from_gaussians(iproc,nproc,KSwfn%orbs,KSwfn%Lzd,&
          KSwfn%Lzd%hgrids(1),KSwfn%Lzd%hgrids(2),KSwfn%Lzd%hgrids(3),&
          KSwfn%psi,KSwfn%gbd,KSwfn%gaucoeffs)

  case(INPUT_PSI_DISK_GAUSS)
     !reading wavefunctions from gaussian file
     if (iproc == 0) then
        !write( *,'(1x,a)')&
        !     &   '------------------------------------------- Reading Wavefunctions from gaussian file'
        call yaml_comment('Reading Wavefunctions from gaussian file',hfill='-')
        call yaml_open_map("Input Hamiltonian")
     end if
     call read_gaussian_information(KSwfn%orbs,KSwfn%gbd,KSwfn%gaucoeffs,&
          trim(in%dir_output)//'wavefunctions.gau')
     !associate the new positions, provided that the atom number is good
     if (KSwfn%gbd%nat == atoms%astruct%nat) then
        KSwfn%gbd%rxyz=>rxyz
     else
        !        if (iproc == 0) then
        call yaml_warning('The atom number does not coincide with the number of gaussian centers')
        !write( *,*)&
        !     &   ' ERROR: the atom number does not coincide with the number of gaussian centers'
        !        end if
        stop
     end if
     call restart_from_gaussians(iproc,nproc,KSwfn%orbs,KSwfn%Lzd,&
          KSwfn%Lzd%hgrids(1),KSwfn%Lzd%hgrids(2),KSwfn%Lzd%hgrids(3),&
          KSwfn%psi,KSwfn%gbd,KSwfn%gaucoeffs)

  case (INPUT_PSI_LINEAR_AO)
     if (iproc == 0) then
        !write(*,'(1x,a)')&
        !     '------------------------------------------------------- Input Wavefunctions Creation'
        call yaml_comment('Input Wavefunctions Creation',hfill='-')
        call yaml_open_map("Input Hamiltonian")
     end if

     ! By doing an LCAO input guess
     tmb%can_use_transposed=.false.
     if (.not.present(locregcenters)) stop 'locregcenters not present!'
     call inputguessConfinement(iproc,nproc,atoms,in,KSwfn%Lzd%hgrids(1),KSwfn%Lzd%hgrids(2),KSwfn%Lzd%hgrids(3), &
          rxyz,nlpsp,GPU,KSwfn%orbs,kswfn,tmb,denspot,denspot0,energs,locregcenters)
     if(tmb%can_use_transposed) then
         i_all=-product(shape(tmb%psit_c))*kind(tmb%psit_c)
         deallocate(tmb%psit_c, stat=i_stat)
         call memocc(i_stat, i_all, 'tmb%psit_c', subname)
         i_all=-product(shape(tmb%psit_f))*kind(tmb%psit_f)
         deallocate(tmb%psit_f, stat=i_stat)
         call memocc(i_stat, i_all, 'tmb%psit_f', subname)
     end if
  case (INPUT_PSI_DISK_LINEAR)
     if (iproc == 0) then
        !write( *,'(1x,a)')&
        !     &   '---------------------------------------------------- Reading Wavefunctions from disk'
        call yaml_comment('Reading Wavefunctions from disk',hfill='-')
        call yaml_open_map("Input Hamiltonian")
     end if

     !if (in%lin%scf_mode==LINEAR_FOE) then
     !    stop 'INPUT_PSI_DISK_LINEAR not allowed with LINEAR_FOE!'
     !end if

     ! By reading the basis functions and coefficients from file
     !call readmywaves_linear(iproc,trim(in%dir_output)//'minBasis',&
     !     & input_wf_format,tmb%npsidim_orbs,tmb%lzd,tmb%orbs, &
     !     & atoms,rxyz_old,rxyz,tmb%psi,tmb%coeff)

     call readmywaves_linear_new(iproc,trim(in%dir_output),'minBasis',input_wf_format,&
          atoms,tmb,rxyz_old,rxyz,ref_frags,in%frag,in%lin%fragment_calculation)

     ! normalize tmbs - only really needs doing if we reformatted, but will need to calculate transpose after anyway
     !nullify(tmb%psit_c)                                                                
     !nullify(tmb%psit_f)  

     tmb%can_use_transposed=.true.
     overlap_calculated=.false.
     allocate(tmb%psit_c(sum(tmb%collcom%nrecvcounts_c)), stat=i_stat)
     call memocc(i_stat, tmb%psit_c, 'tmb%psit_c', subname)
     allocate(tmb%psit_f(7*sum(tmb%collcom%nrecvcounts_f)), stat=i_stat)
     call memocc(i_stat, tmb%psit_f, 'tmb%psit_f', subname)

     call transpose_localized(iproc, nproc, tmb%npsidim_orbs, tmb%orbs, tmb%collcom, &
          tmb%psi, tmb%psit_c, tmb%psit_f, tmb%lzd)

     ! normalize psi
     allocate(norm(tmb%orbs%norb), stat=i_stat)
     call memocc(i_stat, norm, 'norm', subname)

     call normalize_transposed(iproc, nproc, tmb%orbs, tmb%collcom, tmb%psit_c, tmb%psit_f, norm)

     call untranspose_localized(iproc, nproc, tmb%npsidim_orbs, tmb%orbs, tmb%collcom, &
          tmb%psit_c, tmb%psit_f, tmb%psi, tmb%lzd)

     i_all = -product(shape(norm))*kind(norm)
     deallocate(norm,stat=i_stat)
     call memocc(i_stat,i_all,'norm',subname)

     !!allocate(tmb%linmat%denskern%matrix(tmb%orbs%norb,tmb%orbs%norb), stat=i_stat)
     !!call memocc(i_stat, tmb%linmat%denskern%matrix, 'tmb%linmat%denskern%matrix', subname)
     !!call calculate_density_kernel(iproc, nproc, .true., KSwfn%orbs, tmb%orbs, tmb%coeff, tmb%linmat%denskern%matrix)
     !!call compress_matrix_for_allreduce(iproc,tmb%linmat%denskern)
     !!do itmb=1,tmb%orbs%norb
     !!   do jtmb=1,tmb%orbs%norb
     !!      write(20,*) itmb,jtmb,tmb%linmat%denskern%matrix(itmb,jtmb)
     !!   end do
     !!end do
     !!i_all=-product(shape(tmb%linmat%denskern%matrix))*kind(tmb%linmat%denskern%matrix)
     !!deallocate(tmb%linmat%denskern%matrix,stat=i_stat)
     !!call memocc(i_stat,i_all,'tmb%linmat%denskern%matrix',subname)           

     ! CDFT: need to do this here to correct fragment charges in case of constrained transfer integral calculation
     call nullify_cdft_data(cdft)
     nullify(in_frag_charge)
     if (in%lin%constrained_dft) then
        call cdft_data_init(cdft,in%frag,KSwfn%Lzd%Glr%d%n1i*KSwfn%Lzd%Glr%d%n2i*denspot%dpbox%n3d,&
             in%lin%calc_transfer_integrals)
        if (in%lin%calc_transfer_integrals) then
           in_frag_charge=f_malloc_ptr(in%frag%nfrag,id='in_frag_charge')
           call dcopy(in%frag%nfrag,in%frag%charge(1),1,in_frag_charge(1),1)
           ! assume all other fragments neutral, use total system charge to get correct charge for the other fragment
           in_frag_charge(cdft%ifrag_charged(2))=in%ncharge - in_frag_charge(cdft%ifrag_charged(1))
           ! want the difference in number of electrons here, rather than explicitly the charge
           ! actually need this to be more general - perhaps change constraint to be charge rather than number of electrons
           cdft%charge=ref_frags(in%frag%frag_index(cdft%ifrag_charged(1)))%nelec-in_frag_charge(cdft%ifrag_charged(1))&
                -(ref_frags(in%frag%frag_index(cdft%ifrag_charged(2)))%nelec-in_frag_charge(cdft%ifrag_charged(2)))
           !DEBUG
           if (iproc==0) then
              print*,'???????????????????????????????????????????????????????'
              print*,'ifrag_charged1&2,in_frag_charge1&2,ncharge,cdft%charge',cdft%ifrag_charged(1:2),&
              in_frag_charge(cdft%ifrag_charged(1)),in_frag_charge(cdft%ifrag_charged(2)),in%ncharge,cdft%charge
              print*,'??',ref_frags(in%frag%frag_index(cdft%ifrag_charged(1)))%nelec,in_frag_charge(cdft%ifrag_charged(1)),&
                              ref_frags(in%frag%frag_index(cdft%ifrag_charged(2)))%nelec,in_frag_charge(cdft%ifrag_charged(2))
              print*,'???????????????????????????????????????????????????????'
           end if
           !END DEBUG
        else
           in_frag_charge=>in%frag%charge
        end if
     else
        in_frag_charge=>in%frag%charge
     end if

     ! we have to copy the coeffs from the fragment structure to the tmb structure and reconstruct each 'mini' kernel
     ! this is overkill as we are recalculating the kernel anyway - fix at some point
     ! or just put into fragment structure to save recalculating for CDFT
     if (in%lin%fragment_calculation) then
        call fragment_coeffs_to_kernel(iproc,in,in_frag_charge,ref_frags,tmb,KSwfn%orbs,overlap_calculated,&
             nstates_max,in%lin%constrained_dft)
        if (in%lin%calc_transfer_integrals.and.in%lin%constrained_dft) then
           call f_free_ptr(in_frag_charge)
        else
           nullify(in_frag_charge)
        end if
     else
        call dcopy(tmb%orbs%norb**2,ref_frags(1)%coeff(1,1),1,tmb%coeff(1,1),1)
        call dcopy(tmb%orbs%norb,ref_frags(1)%eval(1),1,tmb%orbs%eval(1),1)
        if (associated(ref_frags(1)%coeff)) call f_free_ptr(ref_frags(1)%coeff)
        if (associated(ref_frags(1)%eval)) call f_free_ptr(ref_frags(1)%eval)
     end if

     ! hack occup to make density neutral with full occupations, then unhack after extra diagonalization (using nstates max)
     ! use nstates_max - tmb%orbs%occup set in fragment_coeffs_to_kernel
     if (in%lin%diag_start) then
        ! not worrying about this case as not currently used anyway
        call reconstruct_kernel(iproc, nproc, tmb%orthpar%methTransformOverlap, tmb%orthpar%blocksize_pdsyev, &
             tmb%orthpar%blocksize_pdgemm, tmb%orbs, tmb, overlap_calculated)  
     else
        ! come back to this - reconstruct kernel too expensive with exact version, but Taylor needs to be done ~ 3 times here...
        call reconstruct_kernel(iproc, nproc, tmb%orthpar%methTransformOverlap, tmb%orthpar%blocksize_pdsyev, &
             tmb%orthpar%blocksize_pdgemm, KSwfn%orbs, tmb, overlap_calculated)
     end if
     !!tmb%linmat%ovrlp%matrix=f_malloc_ptr((/tmb%orbs%norb,tmb%orbs%norb/),id='tmb%linmat%ovrlp%matrix')
     !!tmb%linmat%denskern%matrix=f_malloc_ptr((/tmb%orbs%norb,tmb%orbs%norb/),id='tmb%linmat%denskern%matrix')
     !!ks=f_malloc((/tmb%orbs%norb,tmb%orbs%norb/),id='ks')
     !!ksk=f_malloc((/tmb%orbs%norb,tmb%orbs%norb/),id='ksk')
     !!call uncompressMatrix(bigdft_mpi%iproc,tmb%linmat%ovrlp)
     !!call uncompressMatrix(bigdft_mpi%iproc,tmb%linmat%denskern)
     !!call dgemm('n', 't', tmb%orbs%norb, tmb%orbs%norb, tmb%orbs%norb, 1.d0, tmb%linmat%denskern%matrix(1,1), tmb%orbs%norb, &
     !!           tmb%linmat%ovrlp%matrix(1,1), tmb%orbs%norb, 0.d0, ks(1,1), tmb%orbs%norb) 
     !!call dgemm('n', 't', tmb%orbs%norb, tmb%orbs%norb, tmb%orbs%norb, 1.d0, ks(1,1), tmb%orbs%norb, &
     !!           tmb%linmat%denskern%matrix(1,1), tmb%orbs%norb, 0.d0, ksk(1,1), tmb%orbs%norb)

     !!nonidem=0
     !!do itmb=1,tmb%orbs%norb
     !!   do jtmb=1,tmb%orbs%norb
     !!      write(61,*) itmb,jtmb,tmb%linmat%denskern%matrix(itmb,jtmb),ksk(itmb,jtmb),&
     !!           tmb%linmat%denskern%matrix(itmb,jtmb)-ksk(itmb,jtmb),tmb%linmat%ovrlp%matrix(itmb,jtmb)
     !!      nonidem=nonidem+tmb%linmat%denskern%matrix(itmb,jtmb)-ksk(itmb,jtmb)
     !!   end do
     !!end do
     !!print*,'non idempotency',nonidem/tmb%orbs%norb**2

     !!call f_free(ks) 
     !!call f_free(ksk) 
     !!call f_free_ptr(tmb%linmat%ovrlp%matrix)   
     !!call f_free_ptr(tmb%linmat%denskern%matrix)   

     tmb%can_use_transposed=.false. ! - do we really need to deallocate here?
     i_all = -product(shape(tmb%psit_c))*kind(tmb%psit_c)                               
     deallocate(tmb%psit_c,stat=i_stat)                                                 
     call memocc(i_stat,i_all,'tmb%psit_c',subname)                                     
     i_all = -product(shape(tmb%psit_f))*kind(tmb%psit_f)                               
     deallocate(tmb%psit_f,stat=i_stat)                                                 
     call memocc(i_stat,i_all,'tmb%psit_f',subname)
     nullify(tmb%psit_c)
     nullify(tmb%psit_f)

     ! Now need to calculate the charge density and the potential related to this inputguess
     call communicate_basis_for_density_collective(iproc, nproc, tmb%lzd, max(tmb%npsidim_orbs,tmb%npsidim_comp), &
          tmb%orbs, tmb%psi, tmb%collcom_sr)

     call sumrho_for_TMBs(iproc, nproc, KSwfn%Lzd%hgrids(1), KSwfn%Lzd%hgrids(2), KSwfn%Lzd%hgrids(3), &
          tmb%collcom_sr, tmb%linmat%denskern_large, KSwfn%Lzd%Glr%d%n1i*KSwfn%Lzd%Glr%d%n2i*denspot%dpbox%n3d, denspot%rhov)

     ! CDFT: calculate w(r) and w_ab, define some initial guess for V and initialize other cdft_data stuff
     call timing(iproc,'constraineddft','ON')
     if (in%lin%constrained_dft) then
        call cdft_data_allocate(cdft,tmb%linmat%ham)
        if (trim(cdft%method)=='fragment_density') then ! fragment density approach
           if (in%lin%calc_transfer_integrals) stop 'Must use Lowdin for CDFT transfer integral calculations for now'
           if (in%lin%diag_start) stop 'Diag at start probably not working for fragment_density'
           cdft%weight_function=f_malloc_ptr(cdft%ndim_dens,id='cdft%weight_function')
           call calculate_weight_function(in,ref_frags,cdft,&
                KSwfn%Lzd%Glr%d%n1i*KSwfn%Lzd%Glr%d%n2i*denspot%dpbox%n3d,denspot%rhov,tmb,atoms,rxyz,denspot)
           call calculate_weight_matrix_using_density(cdft,tmb,atoms,in,GPU,denspot)
           call f_free_ptr(cdft%weight_function)
        else if (trim(cdft%method)=='lowdin') then ! direct weight matrix approach
           call calculate_weight_matrix_lowdin_wrapper(cdft,tmb,in,ref_frags,.false.,tmb%orthpar%methTransformOverlap)
           ! debug
           !call plot_density(iproc,nproc,'initial_density.cube', &
           !     atoms,rxyz,denspot%dpbox,1,denspot%rhov)
           ! debug
        else 
           stop 'Error invalid method for calculating CDFT weight matrix'
        end if
     end if

     call timing(iproc,'constraineddft','OF')

    !call plot_density(bigdft_mpi%iproc,bigdft_mpi%nproc,'density.cube', &
    !     atoms,rxyz,denspot%dpbox,1,denspot%rhov)

     ! Must initialize rhopotold (FOR NOW... use the trivial one)
     call dcopy(max(denspot%dpbox%ndims(1)*denspot%dpbox%ndims(2)*denspot%dpbox%n3p,1)*in%nspin, &
          denspot%rhov(1), 1, denspot0(1), 1)
     !!call deallocateCommunicationbufferSumrho(tmb%comsr, subname)

     call updatePotential(in%ixc,in%nspin,denspot,energs%eh,energs%exc,energs%evxc)
     call local_potential_dimensions(tmb%lzd,tmb%orbs,denspot%dpbox%ngatherarr(0,1))

    !call plot_density(bigdft_mpi%iproc,bigdft_mpi%nproc,'potential.cube', &
    !     atoms,rxyz,denspot%dpbox,1,denspot%rhov)

    !call plot_density(bigdft_mpi%iproc,bigdft_mpi%nproc,'vext.cube', &
    !     atoms,rxyz,denspot%dpbox,1,denspot%V_ext)

     !! if we want to ignore read in coeffs and diag at start - EXPERIMENTAL
     if (in%lin%diag_start) then
        !if (iproc==0) then
        !print*,'coeffs before extra diag:'
        !do iorb=1,KSwfn%orbs%norb
        !write(*,'(I4,3(F8.2,2x),4(F8.4,2x),2x,4(F8.4,2x))') iorb,KSwfn%orbs%occup(iorb),tmb%orbs%occup(iorb),&
        !tmb%orbs%eval(iorb),tmb%coeff(1:4,iorb),tmb%coeff(5:8,iorb)
        !end do
        !do iorb=KSwfn%orbs%norb+1,tmb%orbs%norb
        !write(*,'(I4,3(F8.2,2x),4(F8.4,2x),2x,4(F8.4,2x))') iorb,0.d0,tmb%orbs%occup(iorb),tmb%orbs%eval(iorb),&
        !tmb%coeff(1:4,iorb),tmb%coeff(5:8,iorb)
        !end do
        !end if
        call get_coeff(iproc,nproc,LINEAR_MIXDENS_SIMPLE,KSwfn%orbs,atoms,rxyz,denspot,GPU,&
             infoCoeff,energs,nlpsp,in%SIC,tmb,pnrm,.false.,.false.,&
             .true.,ham_small,0,0,0,0,in%lin%order_taylor,in%calculate_KS_residue) !in%lin%extra_states) - assume no extra states as haven't set occs for this yet

        !if (iproc==0) then
        !print*,'coeffs after extra diag:'
        !do iorb=1,KSwfn%orbs%norb
        !write(*,'(I4,3(F8.2,2x),4(F8.4,2x),2x,4(F8.4,2x))') iorb,KSwfn%orbs%occup(iorb),tmb%orbs%occup(iorb),tmb%orbs%eval(iorb),&
        !tmb%coeff(1:4,iorb),tmb%coeff(5:8,iorb)
        !end do
        !do iorb=KSwfn%orbs%norb+1,tmb%orbs%norb
        !write(*,'(I4,3(F8.2,2x),4(F8.4,2x),2x,4(F8.4,2x))') iorb,0.d0,tmb%orbs%occup(iorb),tmb%orbs%eval(iorb),&
        !tmb%coeff(1:4,iorb),tmb%coeff(5:8,iorb)
        !end do
        !end if

        !reset occ
        call razero(tmb%orbs%norb,tmb%orbs%occup(1))
        do iorb=1,kswfn%orbs%norb
          tmb%orbs%occup(iorb)=Kswfn%orbs%occup(iorb)
        end do

        ! use the coeffs from the more balanced kernel to get the correctly occupied kernel (maybe reconstruct not needed?)
        call reconstruct_kernel(iproc, nproc, tmb%orthpar%methTransformOverlap, tmb%orthpar%blocksize_pdsyev, &
            tmb%orthpar%blocksize_pdgemm, KSwfn%orbs, tmb, overlap_calculated)     
        !then redo density and potential with correct charge? - for ease doing in linear scaling
     end if

  case default
     !     if (iproc == 0) then
     !write( *,'(1x,a,I0,a)')'ERROR: illegal value of inputPsiId (', in%inputPsiId, ').'
     call yaml_warning('Illegal value of inputPsiId (' // trim(yaml_toa(in%inputPsiId,fmt='(i0)')) // ')')
     call input_psi_help()
     stop
     !     end if

  end select

  !save the previous potential if the rho_work is associated
  if (denspot%rhov_is==KS_POTENTIAL .and. in%iscf==SCF_KIND_GENERALIZED_DIRMIN) then
     if (associated(denspot%rho_work)) then
        call yaml_warning('The reference potential should be empty to correct the hamiltonian!')
        !write(*,*)'ERROR: the reference potential should be empty to correct the hamiltonian!'
        stop
     end if
     call yaml_newline()
     call yaml_comment('Saving the KS potential obtained from IG')
     allocate(denspot%rho_work(denspot%dpbox%ndimpot*denspot%dpbox%nrhodim+ndebug),stat=i_stat)
     call dcopy(denspot%dpbox%ndimpot*denspot%dpbox%nrhodim,&
          denspot%rhov(1),1,denspot%rho_work(1),1)
  end if

  !all the input format need first_orthon except the LCAO input_guess
  ! WARNING: at the momemt the linear scaling version does not need first_orthon.
  ! hpsi and psit have been allocated during the LCAO input guess.
  ! Maybe to be changed later.
  !if (inputpsi /= 0 .and. inputpsi /=-1000) then
  if ( inputpsi /= INPUT_PSI_LCAO .and. inputpsi /= INPUT_PSI_LINEAR_AO .and. &
        inputpsi /= INPUT_PSI_EMPTY .and. inputpsi /= INPUT_PSI_DISK_LINEAR .and. &
        inputpsi /= INPUT_PSI_MEMORY_LINEAR) then
    
     !orthogonalise wavefunctions and allocate hpsi wavefunction (and psit if parallel)
     call first_orthon(iproc,nproc,KSwfn%orbs,KSwfn%Lzd%Glr%wfd,KSwfn%comms,&
          KSwfn%psi,KSwfn%hpsi,KSwfn%psit,in%orthpar)
  end if

  !if (iproc==0 .and. inputpsi /= INPUT_PSI_LINEAR_AO) call yaml_close_map() !input hamiltonian
  if (iproc==0) call yaml_close_map() !input hamiltonian

  if(inputpsi /= INPUT_PSI_LINEAR_AO .and. inputpsi /= INPUT_PSI_DISK_LINEAR .and. &
     inputpsi /= INPUT_PSI_MEMORY_LINEAR) then
     !allocate arrays for the GPU if a card is present
     if (GPUconv) then
        call prepare_gpu_for_locham(KSwfn%Lzd%Glr%d%n1,KSwfn%Lzd%Glr%d%n2,KSwfn%Lzd%Glr%d%n3,&
             in%nspin,&
             KSwfn%Lzd%hgrids(1),KSwfn%Lzd%hgrids(2),KSwfn%Lzd%hgrids(3),&
             KSwfn%Lzd%Glr%wfd,KSwfn%orbs,GPU)
     end if
     !the same with OpenCL, but they cannot exist at same time
     if (OCLconv) then
        call allocate_data_OCL(KSwfn%Lzd%Glr%d%n1,KSwfn%Lzd%Glr%d%n2,KSwfn%Lzd%Glr%d%n3,&
             atoms%astruct%geocode,&
             in%nspin,KSwfn%Lzd%Glr%wfd,KSwfn%orbs,GPU)
        if (iproc == 0) call yaml_comment('GPU data allocated')
        !if (iproc == 0) write(*,*)'GPU data allocated'
     end if
  end if

   ! Emit that new wavefunctions are ready.
   if (inputpsi /= INPUT_PSI_LINEAR_AO .and. inputpsi /= INPUT_PSI_DISK_LINEAR &
        & .and. inputpsi /= INPUT_PSI_MEMORY_LINEAR .and. KSwfn%c_obj /= 0) then
      call kswfn_emit_psi(KSwfn, 0, 0, iproc, nproc)
   end if
   if ((inputpsi == INPUT_PSI_LINEAR_AO .or.&
        inputpsi == INPUT_PSI_DISK_LINEAR .or. &
        inputpsi == INPUT_PSI_MEMORY_LINEAR ).and. tmb%c_obj /= 0) then
      call kswfn_emit_psi(tmb, 0, 0, iproc, nproc)
   end if

   call f_release_routine()

END SUBROUTINE input_wf


!> Check for the input psi (wavefunctions)
!! @param inputpsi            
!!    INPUT_PSI_DISK_WVL      : psi on the disk (wavelets), check if the wavefunctions are all present
!!                              otherwise switch to normal input guess
!!    INPUT_PSI_DISK_LINEAR : psi on memory (linear version)
!!    INPUT_PSI_LCAO          : Use normal input guess (Linear Combination of Atomic Orbitals)
subroutine input_check_psi_id(inputpsi, input_wf_format, dir_output, orbs, lorbs, iproc, nproc, nfrag, frag_dir, ref_frags)
  use module_types
  use yaml_output
  use module_fragments
  use module_interfaces, except_this_one=>input_check_psi_id
  implicit none
  integer, intent(out) :: input_wf_format         !< (out) Format of WF
  integer, intent(inout) :: inputpsi              !< (in) indicate how check input psi, (out) give how to build psi
  integer, intent(in) :: iproc                    !< (in)  id proc
  integer, intent(in) :: nproc                    !< (in)  #proc
  integer, intent(in) :: nfrag                    !< number of fragment directories which need checking
  type(system_fragment), dimension(:), pointer :: ref_frags  !< number of orbitals for each fragment
  character(len=100), dimension(nfrag), intent(in) :: frag_dir !< label for fragment subdirectories (blank if not a fragment calculation)
  character(len = *), intent(in) :: dir_output
  type(orbitals_data), intent(in) :: orbs, lorbs

  logical :: onefile
  integer :: ifrag

  input_wf_format=WF_FORMAT_NONE !default value
  !for the inputPsi == WF_FORMAT_NONE case, check 
  !if the wavefunctions are all present
  !otherwise switch to normal input guess
  if (inputpsi == INPUT_PSI_DISK_WVL) then
     ! Test ETSF file.
     inquire(file=trim(dir_output)//"wavefunction.etsf",exist=onefile)
     if (onefile) then
        input_wf_format = WF_FORMAT_ETSF
     else
        call verify_file_presence(trim(dir_output)//"wavefunction",orbs,input_wf_format,nproc)
     end if
     if (input_wf_format == WF_FORMAT_NONE) then
        if (iproc==0) call yaml_warning('Missing wavefunction files, switch to normal input guess')
        !if (iproc==0) write(*,*)''
        !if (iproc==0) write(*,*)'*********************************************************************'
        !if (iproc==0) write(*,*)'* WARNING: Missing wavefunction files, switch to normal input guess *'
        !if (iproc==0) write(*,*)'*********************************************************************'
        !if (iproc==0) write(*,*)''
        inputpsi=INPUT_PSI_LCAO
     end if
  end if
  ! Test if the files are there for initialization via reading files
  if (inputpsi == INPUT_PSI_DISK_LINEAR) then
     do ifrag=1,nfrag
        ! Test ETSF file.
        inquire(file=trim(dir_output)//"minBasis.etsf",exist=onefile)
        if (onefile) then
           input_wf_format = WF_FORMAT_ETSF
        else
           call verify_file_presence(trim(dir_output)//trim(frag_dir(ifrag))//"minBasis",lorbs,input_wf_format,&
                nproc,ref_frags(ifrag)%fbasis%forbs%norb)
        end if
        if (input_wf_format == WF_FORMAT_NONE) then
           if (iproc==0) call yaml_warning('Missing wavefunction files, switch to normal input guess')
           !if (iproc==0) write(*,*)''
           !if (iproc==0) write(*,*)'*********************************************************************'
           !if (iproc==0) write(*,*)'* WARNING: Missing wavefunction files, switch to normal input guess *'
           !if (iproc==0) write(*,*)'*********************************************************************'
           !if (iproc==0) write(*,*)''
           inputpsi=INPUT_PSI_LINEAR_AO
           ! if one directoy doesn't exist, exit
           exit
        end if
     end do
  end if
END SUBROUTINE input_check_psi_id


subroutine input_wf_memory_new(nproc, iproc, atoms, &
           rxyz_old, hx_old, hy_old, hz_old, d_old, wfd_old, psi_old,lzd_old, &
           rxyz,hx,hy,hz,d,wfd,psi,orbs,lzd,displ)

  use module_defs
  use ao_inguess, only: atomic_info
  use module_types
  use module_interfaces, except_this_one => input_wf_memory_new
  implicit none

  !Global Variables  
  integer, intent(in) :: nproc, iproc
  type(atoms_data), intent(in) :: atoms
  real(gp), dimension(3, atoms%astruct%nat), intent(in) :: rxyz, rxyz_old
  real(gp), intent(in) :: hx, hy, hz, hx_old, hy_old, hz_old,displ
  type(grid_dimensions), intent(in) :: d, d_old
  type(wavefunctions_descriptors), intent(in) :: wfd
  type(wavefunctions_descriptors), intent(inout) :: wfd_old
  type(orbitals_data), intent(in) :: orbs
  type(local_zone_descriptors), intent(inout) :: lzd_old
  type(local_zone_descriptors), intent(in) :: lzd
  real(wp), dimension(:), pointer :: psi, psi_old

  !Local Variables
  character(len = *), parameter :: subname = "input_wf_memory"
  integer :: i_stat, i_all,iorb,nbox,npsir,ist,i,l,k,i1,i2,i3,l1,l2,l3,p1,p2,p3,ii1,ii2,ii3
  type(workarr_sumrho) :: w
  real(wp), dimension(:,:,:), allocatable :: psir,psir_old
  real(wp) :: hhx_old,hhy_old,hhz_old,hhx,hhy,hhz,dgrid1,dgrid2,dgrid3,expfct,x,y,z,s1,s2,s3
  real(wp) :: s1d1,s1d2,s1d3,s2d1,s2d2,s2d3,s3d1,s3d2,s3d3,norm_1,norm_2,norm_3,norm,radius,jacdet
  real(wp), dimension(-1:1) :: coeff,ipv,ipv2


  !To reduce the size, use real(kind=4)
  real(kind=4), dimension(:,:), allocatable :: shift
  real(wp) :: s1_new, s2_new, s3_new,xz,yz,zz,recnormsqr,exp_val, exp_cutoff
  real(wp) :: k1,k2,k3,distance,cutoff

  integer :: istart,irange,iend,rest,ierr, gridx,gridy,gridz,xbox,ybox,zbox,iy,iz

  !Atom description (needed for call to eleconf)
  integer ::nzatom,nvalelec!,nsccode,mxpl,mxchg
  real(wp) :: rcov
  !character(len=2) :: symbol

  if (lzd_old%Glr%geocode .ne. 'F') then
     write(*,*) 'Not implemented for boundary conditions other than free'
     stop
  end if

 ! Daubechies to ISF
  npsir=1
  call initialize_work_arrays_sumrho(Lzd_old%Glr,w)
  nbox = lzd_old%Glr%d%n1i*Lzd_old%Glr%d%n2i*Lzd_old%Glr%d%n3i

  allocate(psir_old(nbox,npsir,orbs%norbp),stat=i_stat)
  call memocc(i_stat,psir_old,'psir_old', subname)
  
  allocate(psir(lzd%Glr%d%n1i*Lzd%Glr%d%n2i*Lzd%Glr%d%n3i,npsir,orbs%norbp),stat=i_stat)
  call memocc(i_stat,psir,'psir',subname)

  allocate(shift(lzd%glr%d%n1i*lzd%glr%d%n2i*lzd%glr%d%n3i,5),stat=i_stat)
  call memocc(i_stat,shift,'shift',subname)  
  
  call to_zero(max(orbs%npsidim_comp,orbs%npsidim_orbs),psi(1)) 
  call to_zero(lzd%Glr%d%n1i*Lzd%Glr%d%n2i*Lzd%Glr%d%n3i*npsir*orbs%norbp,psir(1,1,1)) 
  call to_zero(nbox*npsir*orbs%norbp,psir_old(1,1,1)) 

  call to_zero(lzd%glr%d%n1i*lzd%glr%d%n2i*lzd%glr%d%n3i*5, shift(1,1))

  ist=1
  loop_orbs: do iorb=1,orbs%norbp
     call daub_to_isf(Lzd_old%Glr,w,psi_old(ist),psir_old(1,1,iorb))
     ist=ist+Lzd_old%Glr%wfd%nvctr_c+7*Lzd_old%Glr%wfd%nvctr_f
  end do loop_orbs
  call deallocate_work_arrays_sumrho(w)

  hhx_old = 0.5*hx_old
  hhy_old = 0.5*hy_old
  hhz_old = 0.5*hz_old  

  hhx = 0.5*hx
  hhy = 0.5*hy
  hhz = 0.5*hz  
  
  jacdet = 0.d0
  expfct = 0.d0
  recnormsqr = 0d0

  irange = int(atoms%astruct%nat/nproc)

  rest = atoms%astruct%nat - nproc*irange

  do i = 0,rest-1
     if(iproc .eq. i) irange = irange + 1 
  end do

  if (iproc <= rest-1) then
     istart = iproc*irange + 1
     iend = istart + irange - 1
  else  
     istart = iproc*irange + 1 + rest
     iend = istart + irange - 1
  end if
   
  do k = istart,iend

     !determine sigma of gaussian (sigma is taken as the covalent radius of the atom,rcov)
     nzatom = atoms%nzatom(atoms%astruct%iatype(k))
     nvalelec =  atoms%nelpsp(atoms%astruct%iatype(k))
     call atomic_info(nzatom, nvalelec,rcov=rcov)
     !call eleconf(nzatom, nvalelec,symbol,rcov,rprb,ehomo,neleconf,nsccode,mxpl,mxchg,amu)
     
     radius = 1.0/((rcov)**2)
     cutoff = 3*rcov
     
     !dimensions for box around atom 
     xbox = nint(cutoff/hhx) ; ybox = nint(cutoff/hhy) ; zbox = nint(cutoff/hhz)
  
     gridx = nint(rxyz(1,k)/hhx)+14
     gridy = nint(rxyz(2,k)/hhy)+14
     gridz = nint(rxyz(3,k)/hhz)+14
  
     !$OMP PARALLEL DO DEFAULT(PRIVATE) SHARED(rxyz_old,rxyz,shift,gridx,gridy,gridz,xbox,ybox,zbox, &
     !$OMP& hhz,hhy,hhx,lzd,k) FIRSTPRIVATE(radius,cutoff)
     do i3 = max(1, gridz-zbox), min(lzd%glr%d%n3i,gridz+zbox)
        ii3 = i3-14
        zz = ii3*hhz
        iz = (i3-1)*lzd%glr%d%n1i*lzd%glr%d%n2i
        do i2 = max(1,gridy-ybox), min(lzd%glr%d%n2i,gridy+ybox)
           ii2 = i2 - 14
           yz = ii2*hhy
           iy = (i2-1)*lzd%glr%d%n1i
           do i1 = max(1,gridx-xbox), min(lzd%glr%d%n1i,gridx+xbox)
               
              ii1 = i1 - 14
              xz = ii1*hhx
              
              norm_1 = 0d0 ; norm_2 = 0d0 ; norm_3 = 0d0 
              s1_new = 0d0 ; s2_new = 0d0 ; s3_new = 0d0

              s1d1 = 0d0 ; s1d2 = 0d0 ; s1d3 = 0d0
              s2d1 = 0d0 ; s2d2 = 0d0 ; s2d3 = 0d0
              s3d1 = 0d0 ; s3d2 = 0d0 ; s3d3 = 0d0 
              
              distance = sqrt((xz-rxyz(1,k))**2+(yz-rxyz(2,k))**2+(zz-rxyz(3,k))**2)
              if(distance > cutoff) cycle 
              
              exp_val = 0.5*(distance**2)*radius
              exp_cutoff = 0.5*(cutoff**2)*radius
              
              expfct = ex(exp_val, exp_cutoff) 

              norm = expfct
              recnormsqr = 1/expfct**2

              s1_new = (rxyz(1,k) - rxyz_old(1,k))*expfct
              s2_new = (rxyz(2,k) - rxyz_old(2,k))*expfct
              s3_new = (rxyz(3,k) - rxyz_old(3,k))*expfct

              norm_1 =  expfct*((xz-rxyz(1,k))*radius)
              norm_2 =  expfct*((yz-rxyz(2,k))*radius)
              norm_3 =  expfct*((zz-rxyz(3,k))*radius)

              s1d1 = s1_new*((xz-rxyz(1,k))*radius)
              s1d2 = s1_new*((yz-rxyz(2,k))*radius)
              s1d3 = s1_new*((zz-rxyz(3,k))*radius)
              
              s2d1 = s2_new*((xz-rxyz(1,k))*radius)
              s2d2 = s2_new*((yz-rxyz(2,k))*radius)
              s2d3 = s2_new*((zz-rxyz(3,k))*radius)
     
              s3d1 = s3_new*((xz-rxyz(1,k))*radius)
              s3d2 = s3_new*((yz-rxyz(2,k))*radius)
              s3d3 = s3_new*((zz-rxyz(3,k))*radius)

              s1d1 =       (s1d1*expfct - s1_new*norm_1)*recnormsqr
              s1d2 =       (s1d2*expfct - s1_new*norm_2)*recnormsqr
              s1d3 =       (s1d3*expfct - s1_new*norm_3)*recnormsqr
         
              s2d1 =       (s2d1*expfct - s2_new*norm_1)*recnormsqr
              s2d2 =       (s2d2*expfct - s2_new*norm_2)*recnormsqr
              s2d3 =       (s2d3*expfct - s2_new*norm_3)*recnormsqr
         
              s3d1 =       (s3d1*expfct - s3_new*norm_1)*recnormsqr
              s3d2 =       (s3d2*expfct - s3_new*norm_2)*recnormsqr
              s3d3 =       (s3d3*expfct - s3_new*norm_3)*recnormsqr
        
 
              jacdet = s1d1*s2d2*s3d3 + s1d2*s2d3*s3d1 + s1d3*s2d1*s3d2 - s1d3*s2d2*s3d1-s1d2*s2d1*s3d3 - s1d1*s2d3*s3d2&
                        &  + s1d1 + s2d2 + s3d3 +s1d1*s2d2+s3d3*s1d1+s3d3*s2d2 - s1d2*s2d1 - s3d2*s2d3 - s3d1*s1d3

              shift(i1+iy+iz,1) = real(s1_new,kind=4) +  shift(i1+iy+iz,1)  
              shift(i1+iy+iz,2) = real(s2_new,kind=4) +  shift(i1+iy+iz,2)  
              shift(i1+iy+iz,3) = real(s3_new,kind=4) +  shift(i1+iy+iz,3)  
              shift(i1+iy+iz,4) = real(expfct,kind=4) +  shift(i1+iy+iz,4)  
              shift(i1+iy+iz,5) = real(jacdet,kind=4) +  shift(i1+iy+iz,5)  

           end do     
        end do
      end do
    !$OMP END PARALLEL DO
  end do

 call MPIALLRED(shift(1,1),lzd%glr%d%n1i*lzd%glr%d%n2i*lzd%glr%d%n3i*5, MPI_SUM,bigdft_mpi%mpi_comm,ierr) 

!Interpolation
 do iorb = 1,orbs%norbp
  !$OMP PARALLEL DO DEFAULT(PRIVATE) SHARED(shift,hhx,hhy,hhz,hhx_old,hhy_old,hhz_old,lzd_old,lzd,atoms, &
  !$OMP& psir_old,psir,rxyz_old,rxyz,iorb) 
  do i3 = 1, lzd%glr%d%n3i
      iz = (i3-1)*lzd%glr%d%n2i*lzd%glr%d%n1i
    do i2 = 1, lzd%glr%d%n2i
       iy = (i2-1)*lzd%glr%d%n1i
      do i1 = 1, lzd%glr%d%n1i

         s1 = shift(i1+iy+iz,1) 
         s2 = shift(i1+iy+iz,2) 
         s3 = shift(i1+iy+iz,3) 
         norm = shift(i1+iy+iz,4) 
         jacdet = shift(i1+iy+iz,5)

         jacdet = jacdet + 1.0

         if(norm.eq.0d0) norm = 1d0

         s1 = s1/norm
         s2 = s2/norm
         s3 = s3/norm

         k1 = (i1-14)*hhx_old 
         k2 = (i2-14)*hhy_old 
         k3 = (i3-14)*hhz_old 

         x = k1 - s1
         y = k2 - s2
         z = k3 - s3
     
         p1 = nint(x/hhx_old) + 14
         p2 = nint(y/hhy_old) + 14 
         p3 = nint(z/hhz_old) + 14
        
         dgrid1 = x - (p1-14)*hhx_old
         dgrid2 = y - (p2-14)*hhy_old
         dgrid3 = z - (p3-14)*hhz_old

         if(p1 < 2 .or. p1 > lzd_old%glr%d%n1i-1 .or. p2 < 2 .or. p2 > lzd_old%glr%d%n2i-1 .or. &
            p3 < 2 .or. p3 > lzd_old%glr%d%n3i-1 ) then
            psir(i1+iy+iz,1,iorb) = 0.d0
            cycle
         end if

         do i = -1,1
           do l = -1,1

              l1 = (p1-1) 
              l2 = (p2+i-1)*lzd_old%glr%d%n1i          
              l3 = (p3+l-1)*lzd_old%glr%d%n2i*lzd_old%glr%d%n1i
              
              coeff(-1) = psir_old(l1+l2+l3,1,iorb)

              l1 = p1  
            
              coeff(0) = (psir_old(l1+l2+l3,1,iorb)-coeff(-1))/hhx_old

              l1 = p1+1 
            
              coeff(1) = (psir_old(l1+l2+l3,1,iorb)-coeff(-1)-coeff(0)*2*hhx_old)/(2*hhx_old*hhx_old)
 
              ipv(l) = coeff(-1) + coeff(0)*(hhx_old+dgrid1) + coeff(1)*(hhx_old+dgrid1)*dgrid1

           end do

           coeff(-1) = ipv(-1)
           coeff(0) = (ipv(0) - coeff(-1))/hhz_old
           coeff(1) = (ipv(1) - coeff(-1) - coeff(0)*2*hhz_old)/(2*hhz_old*hhz_old)

           ipv2(i) = coeff(-1) + coeff(0)*(hhz_old+dgrid3) + coeff(1)*(hhz_old+dgrid3)*dgrid3

         end do

         coeff(-1) = ipv2(-1)
         coeff(0) = (ipv2(0) - coeff(-1))/hhy_old
         coeff(1) = (ipv2(1) - coeff(-1) - coeff(0)*2*hhy_old)/(2*hhy_old*hhy_old)
 
         psir(i1+iy+iz,1,iorb) = &
         & (coeff(-1) + coeff(0)*(dgrid2+hhy_old) + coeff(1)*(dgrid2+hhy_old)*dgrid2) /sqrt(abs(jacdet))
        
       end do
      end do
     end do
  !$OMP END PARALLEL DO
  end do

  call initialize_work_arrays_sumrho(Lzd%Glr,w) 
 
  ist=1
  loop_orbs_back: do iorb=1,orbs%norbp
        call isf_to_daub(Lzd%Glr,w,psir(1,1,iorb),psi(ist))
        ist=ist+Lzd%Glr%wfd%nvctr_c+7*Lzd%Glr%wfd%nvctr_f
  end do loop_orbs_back

  call deallocate_work_arrays_sumrho(w)
  
  i_all = -product(shape(psir_old))*kind(psir_old)
  deallocate(psir_old,stat=i_stat)
  call memocc(i_stat,i_all, 'psir_old', subname)

  i_all = -product(shape(psir))*kind(psir)
  deallocate(psir,stat=i_stat)
  call memocc(i_stat,i_all, 'psir', subname)

  i_all = -product(shape(shift))*kind(shift)
  deallocate(shift,stat=i_stat)
  call memocc(i_stat,i_all, 'shift', subname)

  i_all=-product(shape(psi_old))*kind(psi_old)
  deallocate(psi_old,stat=i_stat)
  call memocc(i_stat,i_all,'psi_old',subname)

contains

  real(wp) function ex(x,m)
     implicit none
     real(wp),intent(in) :: x,m

     ex = (1.0 - x/m)**m

  end function ex

END SUBROUTINE input_wf_memory_new

