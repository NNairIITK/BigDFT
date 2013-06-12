!> @file
!!  Routines to initialize the information about localisation regions
!! @author
!!    Copyright (C) 2007-2011 BigDFT group
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
  logical :: my_output_denspot
  logical, dimension(:,:,:), allocatable :: logrid_c,logrid_f

  call timing(iproc,'CrtDescriptors','ON')

  if (iproc == 0) then
     call yaml_open_map('Wavefunctions Descriptors, full simulation domain')
     !write(*,'(1x,a)')&
     !   &   '------------------------------------------------- Wavefunctions Descriptors Creation'
  end if

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

   if (iproc == 0) then
      !write(*,'(2(1x,a,i10))') &
      !     &   'Coarse resolution grid: Number of segments= ',Glr%wfd%nseg_c,'points=',Glr%wfd%nvctr_c
      call yaml_open_map('Coarse resolution grid')!,flow=.true.)
      call yaml_map('No. of segments',Glr%wfd%nseg_c)
      call yaml_map('No. of points',Glr%wfd%nvctr_c)
      call yaml_close_map()
   end if

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

   if (iproc == 0) then
      !write(*,'(2(1x,a,i10))')
      !'  Fine resolution grid: Number of segments= ',Glr%wfd%nseg_f,'points=',Glr%wfd%nvctr_f
      call yaml_open_map('Fine resolution grid')!,flow=.true.)
      call yaml_map('No. of segments',Glr%wfd%nseg_f)
      call yaml_map('No. of points',Glr%wfd%nvctr_f)
      call yaml_close_map()
      call yaml_close_map()
   end if

   ! Create the file grid.xyz to visualize the grid of functions
   my_output_denspot = .false.
   if (present(output_denspot)) my_output_denspot = output_denspot
   if (my_output_denspot) then
      open(unit=22,file='grid.xyz',status='unknown')
      write(22,*) Glr%wfd%nvctr_c+Glr%wfd%nvctr_f+atoms%astruct%nat,' atomic'
      if (atoms%astruct%geocode=='F') then
         write(22,*)'complete simulation grid with low and high resolution points'
      else if (atoms%astruct%geocode =='S') then
         write(22,'(a,2x,3(1x,1pe24.17))')'surface',atoms%astruct%cell_dim(1),atoms%astruct%cell_dim(2),atoms%astruct%cell_dim(3)
      else if (atoms%astruct%geocode =='P') then
         write(22,'(a,2x,3(1x,1pe24.17))')'periodic',atoms%astruct%cell_dim(1),atoms%astruct%cell_dim(2),&
              atoms%astruct%cell_dim(3)
      end if
      do iat=1,atoms%astruct%nat
         write(22,'(a6,2x,3(1x,e12.5),3x)') &
            &   trim(atoms%astruct%atomnames(atoms%astruct%iatype(iat))),rxyz(1,iat),rxyz(2,iat),rxyz(3,iat)
      enddo
      do i3=0,n3  
         do i2=0,n2  
            do i1=0,n1
               if (logrid_c(i1,i2,i3))&
                  &   write(22,'(a4,2x,3(1x,e10.3))') &
                  &   '  g ',real(i1,kind=8)*hx,real(i2,kind=8)*hy,real(i3,kind=8)*hz
            enddo
         enddo
      end do
      do i3=0,n3 
         do i2=0,n2 
            do i1=0,n1
               if (logrid_f(i1,i2,i3))&
                  &   write(22,'(a4,2x,3(1x,e10.3))') &
                  &   '  G ',real(i1,kind=8)*hx,real(i2,kind=8)*hy,real(i3,kind=8)*hz
            enddo
         enddo
      enddo
      close(22)
   endif

   i_all=-product(shape(logrid_c))*kind(logrid_c)
   deallocate(logrid_c,stat=i_stat)
   call memocc(i_stat,i_all,'logrid_c',subname)
   i_all=-product(shape(logrid_f))*kind(logrid_f)
   deallocate(logrid_f,stat=i_stat)
   call memocc(i_stat,i_all,'logrid_f',subname)

   call timing(iproc,'CrtDescriptors','OF')
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
   if (Glr%geocode == 'F') then
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
   call allocate_wfd(Glr%wfd,subname)

   ! now fill the wavefunction descriptor arrays
   ! coarse grid quantities
   call segkeys(n1,n2,n3,0,n1,0,n2,0,n3,logrid_c,Glr%wfd%nseg_c, &
        & Glr%wfd%keyglob(1,1),Glr%wfd%keyvglob(1))
   ! fine grid quantities
   if (Glr%wfd%nseg_f > 0) then
      call segkeys(n1,n2,n3,0,n1,0,n2,0,n3,logrid_f,Glr%wfd%nseg_f, &
           & Glr%wfd%keyglob(1,Glr%wfd%nseg_c+1), Glr%wfd%keyvglob(Glr%wfd%nseg_c+1))
   end if
   i_all = -product(shape(Glr%wfd%keygloc))*kind(Glr%wfd%keygloc)
   deallocate(Glr%wfd%keygloc,stat=i_stat)
   call memocc(i_stat,i_all,'Glr%wfd%keygloc',subname)
   Glr%wfd%keygloc => Glr%wfd%keyglob
   i_all = -product(shape(Glr%wfd%keyvloc))*kind(Glr%wfd%keyvloc)
   deallocate(Glr%wfd%keyvloc,stat=i_stat)
   call memocc(i_stat,i_all,'Glr%wfd%keyvloc',subname)
   Glr%wfd%keyvloc => Glr%wfd%keyvglob
 
   ! Copy the information of keyglob to keygloc for Glr (just pointing leads to problem during the deallocation of wfd)
!!$   do i = lbound(Glr%wfd%keyglob,1),ubound(Glr%wfd%keyglob,1)
!!$      do j = lbound(Glr%wfd%keyglob,2),ubound(Glr%wfd%keyglob,2)
!!$         Glr%wfd%keygloc(i,j) = Glr%wfd%keyglob(i,j)
!!$      end do
!!$   end do

   !for free BC admits the bounds arrays
   if (Glr%geocode == 'F') then
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
subroutine createProjectorsArrays(iproc,lr,rxyz,at,orbs,&
      &   radii_cf,cpmult,fpmult,hx,hy,hz,nlpspd,proj_G,proj)
   use module_base
   use module_types
   use gaussians, only: gaussian_basis
   implicit none
   integer, intent(in) :: iproc
   real(gp), intent(in) :: cpmult,fpmult,hx,hy,hz
   type(locreg_descriptors),intent(in) :: lr
   type(atoms_data), intent(in) :: at
   type(orbitals_data), intent(in) :: orbs
   real(gp), dimension(3,at%astruct%nat), intent(in) :: rxyz
   real(gp), dimension(at%astruct%ntypes,3), intent(in) :: radii_cf
   type(nonlocal_psp_descriptors), intent(out) :: nlpspd
   real(wp), dimension(:), pointer :: proj
   type(gaussian_basis),dimension(at%astruct%ntypes),intent(in) :: proj_G
   !local variables
   character(len=*), parameter :: subname='createProjectorsArrays'
   integer :: n1,n2,n3,nl1,nl2,nl3,nu1,nu2,nu3,mseg,mproj
   integer :: iat,i_stat,i_all,iseg
   logical, dimension(:,:,:), allocatable :: logrid

   !allocate the different localization regions of the projectors
   nlpspd%natoms=at%astruct%nat
   allocate(nlpspd%plr(at%astruct%nat),stat=i_stat)
   do iat=1,at%astruct%nat
      nlpspd%plr(iat)=locreg_null() !< set in types
   end do

  ! define the region dimensions
    n1 = lr%d%n1
    n2 = lr%d%n2
    n3 = lr%d%n3

   ! determine localization region for all projectors, but do not yet fill the descriptor arrays
   allocate(logrid(0:n1,0:n2,0:n3+ndebug),stat=i_stat)
   call memocc(i_stat,logrid,'logrid',subname)

   call localize_projectors(iproc,n1,n2,n3,hx,hy,hz,cpmult,fpmult,&
        rxyz,radii_cf,logrid,at,orbs,nlpspd,proj_G)

   !here the allocation is possible
   do iat=1,nlpspd%natoms
      !for the moments the bounds are not needed for projectors
      call allocate_wfd(nlpspd%plr(iat)%wfd,subname)
   end do

   allocate(proj(nlpspd%nprojel+ndebug),stat=i_stat)
   call memocc(i_stat,proj,'proj',subname)
   if (nlpspd%nprojel >0) call to_zero(nlpspd%nprojel,proj(1))

   ! After having determined the size of the projector descriptor arrays fill them
   do iat=1,at%astruct%nat
      if(at%npspcode(at%astruct%iatype(iat))==7) then  
        call numb_proj_paw_tr(at%astruct%iatype(iat),at%astruct%ntypes,proj_G(at%astruct%iatype(iat)),mproj)
      else
        call numb_proj(at%astruct%iatype(iat),at%astruct%ntypes,at%psppar,at%npspcode,mproj)
      end if
      if (mproj.ne.0) then 

         call bounds_to_plr_limits(.false.,1,nlpspd%plr(iat),&
              nl1,nl2,nl3,nu1,nu2,nu3)         

         call fill_logrid(at%astruct%geocode,n1,n2,n3,nl1,nu1,nl2,nu2,nl3,nu3,0,1,  &
              at%astruct%ntypes,at%astruct%iatype(iat),rxyz(1,iat),radii_cf(1,3),&
              cpmult,hx,hy,hz,logrid)

         call segkeys(n1,n2,n3,nl1,nu1,nl2,nu2,nl3,nu3,logrid,&
              nlpspd%plr(iat)%wfd%nseg_c,&
              nlpspd%plr(iat)%wfd%keyglob(1,1),nlpspd%plr(iat)%wfd%keyvglob(1))

        call transform_keyglob_to_keygloc(lr,nlpspd%plr(iat),nlpspd%plr(iat)%wfd%nseg_c,&
             nlpspd%plr(iat)%wfd%keyglob(1,1),nlpspd%plr(iat)%wfd%keygloc(1,1))

         ! fine grid quantities
         call bounds_to_plr_limits(.false.,2,nlpspd%plr(iat),&
              nl1,nl2,nl3,nu1,nu2,nu3)         

         call fill_logrid(at%astruct%geocode,n1,n2,n3,nl1,nu1,nl2,nu2,nl3,nu3,0,1,  &
            &   at%astruct%ntypes,at%astruct%iatype(iat),rxyz(1,iat),radii_cf(1,2),&
            fpmult,hx,hy,hz,logrid)

         mseg=nlpspd%plr(iat)%wfd%nseg_f
         iseg=nlpspd%plr(iat)%wfd%nseg_c+1

         if (mseg > 0) then
            call segkeys(n1,n2,n3,nl1,nu1,nl2,nu2,nl3,nu3,  & 
                 logrid,mseg,nlpspd%plr(iat)%wfd%keyglob(1,iseg),&
                 nlpspd%plr(iat)%wfd%keyvglob(iseg))

            call transform_keyglob_to_keygloc(lr,nlpspd%plr(iat),mseg,nlpspd%plr(iat)%wfd%keyglob(1,iseg),&
                 nlpspd%plr(iat)%wfd%keygloc(1,iseg)) 
         end if
      endif
   enddo

   i_all=-product(shape(logrid))*kind(logrid)
   deallocate(logrid,stat=i_stat)
   call memocc(i_stat,i_all,'logrid',subname)
   !fill the projectors if the strategy is a distributed calculation
   if (.not. DistProjApply) then
      !calculate the wavelet expansion of projectors
     call fill_projectors(iproc,lr,hx,hy,hz,at,orbs,rxyz,nlpspd,proj,0)
   end if

END SUBROUTINE createProjectorsArrays


!> Fill the preconditioning projectors for a given atom 
subroutine fillPcProjOnTheFly(PPD, Glr, iat, at, hx,hy,hz,startjorb,ecut_pc,   initial_istart_c ) 
   use module_interfaces
   use module_base
   use module_types
   implicit none
   type(pcproj_data_type),  intent(in) ::PPD
   type(locreg_descriptors),  intent(in):: Glr
   integer, intent(in)  ::iat, startjorb
   real(gp), intent(in) ::  ecut_pc, hx,hy,hz
   !! real(gp), pointer :: gaenes(:)
   integer, intent(in) :: initial_istart_c
   type(atoms_data), intent(in) :: at

   ! local variables  
   type(locreg_descriptors) :: Plr
   real(gp) kx, ky, kz
   integer :: jorb, ncplx, istart_c
   real(wp), dimension(PPD%G%ncoeff ) :: Gocc
   character(len=*), parameter :: subname='fillPcProjOnTheFly'

   istart_c=initial_istart_c

   Plr%d%n1 = Glr%d%n1
   Plr%d%n2 = Glr%d%n2
   Plr%d%n3 = Glr%d%n3
   Plr%geocode = at%astruct%geocode


   call plr_segs_and_vctrs(PPD%pc_nlpspd%plr(iat),&
        Plr%wfd%nseg_c,Plr%wfd%nseg_f,Plr%wfd%nvctr_c,Plr%wfd%nvctr_f)
!!$   Plr%wfd%nvctr_c  =PPD%pc_nlpspd%nvctr_p(2*iat-1)-PPD%pc_nlpspd%nvctr_p(2*iat-2)
!!$   Plr%wfd%nvctr_f  =PPD%pc_nlpspd%nvctr_p(2*iat  )-PPD%pc_nlpspd%nvctr_p(2*iat-1)
!!$   Plr%wfd%nseg_c   =PPD%pc_nlpspd%nseg_p(2*iat-1)-PPD%pc_nlpspd%nseg_p(2*iat-2)
!!$   Plr%wfd%nseg_f   =PPD%pc_nlpspd%nseg_p(2*iat  )-PPD%pc_nlpspd%nseg_p(2*iat-1)

   call allocate_wfd(Plr%wfd,subname)

   call vcopy(Plr%wfd%nseg_c+Plr%wfd%nseg_f,&
        PPD%pc_nlpspd%plr(iat)%wfd%keyvglob(1),1,Plr%wfd%keyvglob(1),1)
   call vcopy(2*(Plr%wfd%nseg_c+Plr%wfd%nseg_f),&
        PPD%pc_nlpspd%plr(iat)%wfd%keyglob(1,1),1,Plr%wfd%keyglob(1,1),1)

!!$   Plr%wfd%keyv(:)  = &
!!$        PPD%pc_nlpspd%keyv_p(  PPD%pc_nlpspd%nseg_p(2*iat-2)+1:  PPD%pc_nlpspd%nseg_p(2*iat)   )
!!$   Plr%wfd%keyg(1:2, :)  = &
!!$        PPD%pc_nlpspd%keyg_p( 1:2,  PPD%pc_nlpspd%nseg_p(2*iat-2)+1:  PPD%pc_nlpspd%nseg_p(2*iat)   )

   kx=0.0_gp
   ky=0.0_gp
   kz=0.0_gp

   Gocc=0.0_wp

   jorb=startjorb

   do while( jorb<=PPD%G%ncoeff .and. PPD%iorbtolr(jorb)== iat) 
      if( PPD%gaenes(jorb)<ecut_pc) then

         Gocc(jorb)=1.0_wp
         ncplx=1
         call gaussians_to_wavelets_orb(ncplx,Plr,hx,hy,hz,kx,ky,kz,PPD%G,&
              Gocc(1),PPD%pc_proj(istart_c))
         Gocc(jorb)=0.0_wp

         !! ---------------  use this to plot projectors
         !!$              write(orbname,'(A,i4.4)')'pc_',iproj
         !!$              Plr%bounds = Glr%bounds
         !!$              Plr%d          = Glr%d
         !!$              call plot_wf_cube(orbname,at,Plr,hx,hy,hz,PPD%G%rxyz, PPD%pc_proj(istart_c) ,"1234567890" ) 

         istart_c=istart_c + (   Plr%wfd%nvctr_c    +   7*Plr%wfd%nvctr_f   )


      endif
      jorb=jorb+1

      if(jorb> PPD%G%ncoeff) exit

   end do

   call deallocate_wfd(Plr%wfd,subname)

END SUBROUTINE fillPcProjOnTheFly


!> Fill the preconditioning projectors for a given atom 
subroutine fillPawProjOnTheFly(PAWD, Glr, iat,  hx,hy,hz,kx,ky,kz,startjorb,   initial_istart_c, geocode, at, iatat) 
   use module_interfaces
   use module_base
   use module_types
   implicit none
   type(pawproj_data_type),  intent(in) ::PAWD
   type(locreg_descriptors),  intent(in):: Glr
   integer, intent(in)  ::iat, startjorb
   real(gp), intent(in) ::   hx,hy,hz,kx,ky,kz
   integer, intent(in) :: initial_istart_c
   character(len=1), intent(in) :: geocode
   type(atoms_data) :: at
   integer :: iatat

   ! local variables  
   type(locreg_descriptors) :: Plr
   integer :: jorb, ncplx, istart_c
   real(wp), dimension(PAWD%G%ncoeff ) :: Gocc
   character(len=*), parameter :: subname='fillPawProjOnTheFly'


   !!Just for extracting the covalent radius and rprb
   integer :: nsccode,mxpl,mxchg
   real(gp) ::amu,rprb,ehomo,rcov, cutoff
   character(len=2) :: symbol
   real(kind=8), dimension(6,4) :: neleconf

   istart_c=initial_istart_c

   Plr%d%n1 = Glr%d%n1
   Plr%d%n2 = Glr%d%n2
   Plr%d%n3 = Glr%d%n3
   Plr%geocode = geocode

   call plr_segs_and_vctrs(PAWD%paw_nlpspd%plr(iat),&
        Plr%wfd%nseg_c,Plr%wfd%nseg_f,Plr%wfd%nvctr_c,Plr%wfd%nvctr_f)

!!$   Plr%wfd%nvctr_c  =PAWD%paw_nlpspd%nvctr_p(2*iat-1)-PAWD%paw_nlpspd%nvctr_p(2*iat-2)
!!$   Plr%wfd%nvctr_f  =PAWD%paw_nlpspd%nvctr_p(2*iat  )-PAWD%paw_nlpspd%nvctr_p(2*iat-1)
!!$   Plr%wfd%nseg_c   =PAWD%paw_nlpspd%nseg_p(2*iat-1)-PAWD%paw_nlpspd%nseg_p(2*iat-2)
!!$   Plr%wfd%nseg_f   =PAWD%paw_nlpspd%nseg_p(2*iat  )-PAWD%paw_nlpspd%nseg_p(2*iat-1)

   call allocate_wfd(Plr%wfd,subname)

   call vcopy(Plr%wfd%nseg_c+Plr%wfd%nseg_f,&
        PAWD%paw_nlpspd%plr(iat)%wfd%keyvglob(1),1,Plr%wfd%keyvglob(1),1)
   call vcopy(2*(Plr%wfd%nseg_c+Plr%wfd%nseg_f),&
        PAWD%paw_nlpspd%plr(iat)%wfd%keyglob(1,1),1,Plr%wfd%keyglob(1,1),1)

!!$   Plr%wfd%keyv(:)  = PAWD%paw_nlpspd%keyv_p(  PAWD%paw_nlpspd%nseg_p(2*iat-2)+1:  PAWD%paw_nlpspd%nseg_p(2*iat)   )
!!$   Plr%wfd%keyg(1:2, :)  = PAWD%paw_nlpspd%keyg_p( 1:2,  PAWD%paw_nlpspd%nseg_p(2*iat-2)+1:  PAWD%paw_nlpspd%nseg_p(2*iat)   )

   if (kx**2 + ky**2 + kz**2 == 0.0_gp) then
      ncplx=1
   else
      ncplx=2
   end if

   Gocc=0.0_wp

   jorb=startjorb

   !!Just for extracting the covalent radius 
   call eleconf(at%nzatom( at%astruct%iatype(iatat)), at%nelpsp(at%astruct%iatype(iatat)) ,  &
      &   symbol, rcov, rprb, ehomo,neleconf, nsccode, mxpl, mxchg, amu)

   cutoff=rcov*1.5_gp

   do while( jorb<=PAWD%G%ncoeff         .and. PAWD%iorbtolr(jorb)== iat)      
      Gocc(jorb)=1.0_wp

      call gaussians_c_to_wavelets_orb(ncplx,Plr,hx,hy,hz,kx,ky,kz,PAWD%G,&
         &   Gocc(1),  PAWD%paw_proj(istart_c), cutoff  )

      Gocc(jorb)=0.0_wp
      !!$     !! ---------------  use this to plot projectors
      !!$              write(orbname,'(A,i4.4)')'paw2_',jorb
      !!$              Plr%bounds = Glr%bounds
      !!$              Plr%d          = Glr%d
      !!$              call plot_wf_cube(orbname,PAWD%at,Plr,hx,hy,hz,PAWD%G%rxyz, PAWD%paw_proj(istart_c) ,"1234567890" ) 

      istart_c=istart_c + (   Plr%wfd%nvctr_c    +   7*Plr%wfd%nvctr_f   ) * ncplx


      jorb=jorb+1

      if(jorb> PAWD%G%ncoeff) exit

   end do

   call deallocate_wfd(Plr%wfd,subname)

END SUBROUTINE fillPawProjOnTheFly


!>   Determine localization region for all preconditioning projectors, but do not yet fill the descriptor arrays
subroutine createPcProjectorsArrays(iproc,n1,n2,n3,rxyz,at,orbs,&
      &   radii_cf,cpmult,fpmult,hx,hy,hz, ecut_pc, &
      &   PPD, Glr)
   use module_interfaces, except_this_one => createPcProjectorsArrays
   use module_base
   use module_types
   implicit none
   integer, intent(in) :: iproc,n1,n2,n3
   real(gp), intent(in) :: cpmult,fpmult,hx,hy,hz
   type(atoms_data), intent(in) :: at
   type(orbitals_data), intent(in) :: orbs

   real(gp), dimension(3,at%astruct%nat), intent(in) :: rxyz
   real(gp), dimension(at%astruct%ntypes,3), intent(in) :: radii_cf
   real(gp), intent(in):: ecut_pc

   type(pcproj_data_type) ::PPD

   type(locreg_descriptors),  intent(in):: Glr


   !local variables
   character(len=*), parameter :: subname='createPcProjectorsArrays'
   integer :: nl1,nl2,nl3,nu1,nu2,nu3,mseg,mproj, mvctr
   integer :: iat,i_stat,i_all,iseg, istart_c
   logical, dimension(:,:,:), allocatable :: logrid


   integer :: ng
   logical :: enlargerprb
   real(wp), dimension(:), pointer :: Gocc

   integer, pointer :: iorbto_l(:)
   integer, pointer :: iorbto_m(:)
   integer, pointer :: iorbto_ishell(:)
   integer, pointer :: iorbto_iexpobeg(:)

   integer :: nspin
   integer ::  jorb
   integer :: iproj, startjorb
   real(gp) :: Pcpmult
   integer :: mprojtot, nvctr_c, nvctr_f
   integer :: nprojel_tmp

   Pcpmult=1.5*cpmult

   ng=21
   enlargerprb = .false.
   nspin=1


   nullify(PPD%G%rxyz)
   call gaussian_pswf_basis(ng,enlargerprb,iproc,nspin,at,rxyz,PPD%G,Gocc, PPD%gaenes, &
      &   PPD%iorbtolr,iorbto_l, iorbto_m,  iorbto_ishell,iorbto_iexpobeg  )  


   ! allocated  : gaenes, Gocc , PPD%iorbtolr,iorbto_l, iorbto_m,  iorbto_ishell,iorbto_iexpobeg


   !!$ ========================================================================================



   PPD%pc_nlpspd%natoms=at%astruct%nat
   allocate(PPD%pc_nlpspd%plr(at%astruct%nat),stat=i_stat)

!!$   allocate(PPD%pc_nlpspd%nseg_p(0:2*at%astruct%nat+ndebug),stat=i_stat)
!!$   call memocc(i_stat,PPD%pc_nlpspd%nseg_p,'pc_nlpspd%nseg_p',subname)
!!$   allocate(PPD%pc_nlpspd%nvctr_p(0:2*at%astruct%nat+ndebug),stat=i_stat)
!!$   call memocc(i_stat,PPD%pc_nlpspd%nvctr_p,'pc_nlpspd%nvctr_p',subname)
!!$   allocate(PPD%pc_nlpspd%nboxp_c(2,3,at%astruct%nat+ndebug),stat=i_stat)
!!$   call memocc(i_stat,PPD%pc_nlpspd%nboxp_c,'pc_nlpspd%nboxp_c',subname)
!!$   allocate(PPD%pc_nlpspd%nboxp_f(2,3,at%astruct%nat+ndebug),stat=i_stat)
!!$   call memocc(i_stat,PPD%pc_nlpspd%nboxp_f,'pc_nlpspd%nboxp_f',subname)

   allocate(logrid(0:n1,0:n2,0:n3+ndebug),stat=i_stat)
   call memocc(i_stat,logrid,'logrid',subname)


   call localize_projectors(iproc,n1,n2,n3,hx,hy,hz,Pcpmult,fpmult,rxyz,radii_cf,&
      &   logrid,at,orbs,PPD%pc_nlpspd)

   ! the above routine counts atomic projector and the number of their element for psp
   ! We must therefore correct , later, nlpspd%nprojel  and nlpspd%nproj
   !-------------------

   ! allocations for arrays holding the projectors and their data descriptors
   !here the allocation is possible
   do iat=1,PPD%pc_nlpspd%natoms
      !for the moments the bounds are not needed for projectors
      call allocate_wfd(PPD%pc_nlpspd%plr(iat)%wfd,subname)
   end do

!!$   allocate(PPD%pc_nlpspd%keyg_p(2,PPD%pc_nlpspd%nseg_p(2*at%astruct%nat)+ndebug),stat=i_stat)
!!$   call memocc(i_stat,PPD%pc_nlpspd%keyg_p,'pc_nlpspd%keyg_p',subname)
!!$
!!$
!!$   allocate(PPD%pc_nlpspd%keyv_p(PPD%pc_nlpspd%nseg_p(2*at%astruct%nat)+ndebug),stat=i_stat)
!!$   call memocc(i_stat,PPD%pc_nlpspd%keyv_p,'pc_nlpspd%keyv_p',subname)



   !!$  -- this one delayed, waiting for the correct pc_nlpspd%nprojel, pc_nlpspd%nproj
   !!$  --
   !!$  allocate(pc_proj(pc_nlpspd%nprojel+ndebug),stat=i_stat)
   !!$  call memocc(i_stat,pc_proj,'pc_proj',subname)
   PPD%ecut_pc=ecut_pc

   PPD%pc_nlpspd%nprojel=0
   PPD%pc_nlpspd%nproj  =0

   !!$ =========================================================================================  

   mprojtot=0
   jorb=1  
   ! After having determined the size of the projector descriptor arrays fill them
   do iat=1,at%astruct%nat

      mproj=0

      do while( jorb<=PPD%G%ncoeff         .and. PPD%iorbtolr(jorb)== iat)

         if( PPD%gaenes(jorb)<ecut_pc) then
            mproj=mproj+1
         endif
         if(jorb==PPD%G%ncoeff) exit
         jorb=jorb+1
      end do

      mprojtot=mprojtot+mproj

      PPD%pc_nlpspd%nproj=PPD%pc_nlpspd%nproj+mproj


      if (mproj.ne.0) then 

         nprojel_tmp=0

         call bounds_to_plr_limits(.false.,1,PPD%pc_nlpspd%plr(iat),nl1,nl2,nl3,nu1,nu2,nu3)
!!$         ! coarse grid quantities
!!$         nl1=PPD%pc_nlpspd%nboxp_c(1,1,iat) 
!!$         nl2=PPD%pc_nlpspd%nboxp_c(1,2,iat) 
!!$         nl3=PPD%pc_nlpspd%nboxp_c(1,3,iat) 
!!$
!!$         nu1=PPD%pc_nlpspd%nboxp_c(2,1,iat)
!!$         nu2=PPD%pc_nlpspd%nboxp_c(2,2,iat)
!!$         nu3=PPD%pc_nlpspd%nboxp_c(2,3,iat)

         call fill_logrid(at%astruct%geocode,n1,n2,n3,nl1,nu1,nl2,nu2,nl3,nu3,0,1,  &
            &   at%astruct%ntypes,at%astruct%iatype(iat),rxyz(1,iat),radii_cf(1,3),Pcpmult,hx,hy,hz,logrid)

!!$         iseg=PPD%pc_nlpspd%nseg_p(2*iat-2)+1
!!$         mseg=PPD%pc_nlpspd%nseg_p(2*iat-1)-PPD%pc_nlpspd%nseg_p(2*iat-2)
         mseg=PPD%pc_nlpspd%plr(iat)%wfd%nseg_c

         call segkeys(n1,n2,n3,nl1,nu1,nl2,nu2,nl3,nu3,  & 
!!$         logrid,mseg,PPD%pc_nlpspd%keyg_p(1,iseg),PPD%pc_nlpspd%keyv_p(iseg))
         logrid,mseg,PPD%pc_nlpspd%plr(iat)%wfd%keyglob(1,1),PPD%pc_nlpspd%plr(iat)%wfd%keyvglob(1))

!!$         mvctr =PPD%pc_nlpspd%nvctr_p(2*iat-1)-PPD%pc_nlpspd%nvctr_p(2*iat-2)
         mvctr =PPD%pc_nlpspd%plr(iat)%wfd%nvctr_c

         nprojel_tmp =nprojel_tmp +mproj*mvctr

         call bounds_to_plr_limits(.false.,2,PPD%pc_nlpspd%plr(iat),nl1,nl2,nl3,nu1,nu2,nu3)
!!$         ! fine grid quantities
!!$         nl1=PPD%pc_nlpspd%nboxp_f(1,1,iat)
!!$         nl2=PPD%pc_nlpspd%nboxp_f(1,2,iat)
!!$         nl3=PPD%pc_nlpspd%nboxp_f(1,3,iat)
!!$
!!$         nu1=PPD%pc_nlpspd%nboxp_f(2,1,iat)
!!$         nu2=PPD%pc_nlpspd%nboxp_f(2,2,iat)
!!$         nu3=PPD%pc_nlpspd%nboxp_f(2,3,iat)
         call fill_logrid(at%astruct%geocode,n1,n2,n3,nl1,nu1,nl2,nu2,nl3,nu3,0,1,  &
            &   at%astruct%ntypes,at%astruct%iatype(iat),rxyz(1,iat),radii_cf(1,2),fpmult,hx,hy,hz,logrid)
!!$         iseg=PPD%pc_nlpspd%nseg_p(2*iat-1)+1
!!$         mseg=PPD%pc_nlpspd%nseg_p(2*iat)-PPD%pc_nlpspd%nseg_p(2*iat-1)
         iseg=PPD%pc_nlpspd%plr(iat)%wfd%nseg_c+1
         mseg=PPD%pc_nlpspd%plr(iat)%wfd%nseg_f

         if (mseg > 0) then
            call segkeys(n1,n2,n3,nl1,nu1,nl2,nu2,nl3,nu3,  & 
                 logrid,mseg,&
!!$                 PPD%pc_nlpspd%keyg_p(1,iseg),PPD%pc_nlpspd%keyv_p(iseg))
                 PPD%pc_nlpspd%plr(iat)%wfd%keyglob(1,iseg),&
                 PPD%pc_nlpspd%plr(iat)%wfd%keyvglob(iseg))

            mvctr =PPD%pc_nlpspd%plr(iat)%wfd%nvctr_f!PPD%pc_nlpspd%nvctr_p(2*iat)-PPD%pc_nlpspd%nvctr_p(2*iat-1)

            nprojel_tmp=nprojel_tmp+mproj*mvctr*7

         end if

         if( PPD%DistProjApply)  then
            PPD%pc_nlpspd%nprojel=max(PPD%pc_nlpspd%nprojel,nprojel_tmp   )
         else
            PPD%pc_nlpspd%nprojel= PPD%pc_nlpspd%nprojel+nprojel_tmp 
         endif


      endif

   enddo


   allocate(PPD%pc_proj(PPD%pc_nlpspd%nprojel+ndebug),stat=i_stat)
   call memocc(i_stat,PPD%pc_proj,'pc_proj',subname)

   allocate(PPD%ilr_to_mproj(at%astruct%nat  +ndebug ) , stat=i_stat)
   call memocc(i_stat,PPD%ilr_to_mproj,'ilr_to_mproj',subname)

   allocate(PPD%iproj_to_ene(mprojtot +ndebug ) , stat=i_stat)
   call memocc(i_stat ,PPD%iproj_to_ene,'iproj_to_ene',subname)

   allocate(PPD%iproj_to_factor(mprojtot +ndebug ) , stat=i_stat)
   call memocc(i_stat ,PPD%iproj_to_factor,'iproj_to_factor',subname)

   allocate(PPD%iproj_to_l(mprojtot +ndebug ) , stat=i_stat)
   call memocc(i_stat ,PPD%iproj_to_l,'iproj_to_l',subname)

   PPD%mprojtot=mprojtot


   startjorb=1
   jorb=1
   istart_c=1
   Gocc(:)=0.0_wp


   iproj=0
   do iat=1,at%astruct%nat

      mproj=0
      do while( jorb<=PPD%G%ncoeff         .and. PPD%iorbtolr(jorb)== iat)
         if( PPD%gaenes(jorb)<ecut_pc) then
            mproj=mproj+1
         endif
         if(jorb==PPD%G%ncoeff) exit
         jorb=jorb+1
      end do

      PPD%ilr_to_mproj(iat)=mproj
      if( mproj>0) then
         nvctr_c  =PPD%pc_nlpspd%plr(iat)%wfd%nvctr_c!PPD%pc_nlpspd%nvctr_p(2*iat-1)-PPD%pc_nlpspd%nvctr_p(2*iat-2)
         nvctr_f  =PPD%pc_nlpspd%plr(iat)%wfd%nvctr_f!PPD%pc_nlpspd%nvctr_p(2*iat  )-PPD%pc_nlpspd%nvctr_p(2*iat-1)

         jorb=startjorb
         do while( jorb<=PPD%G%ncoeff         .and. PPD%iorbtolr(jorb)== iat) 
            if( PPD%gaenes(jorb)<ecut_pc) then
               iproj=iproj+1
               PPD%iproj_to_ene(iproj) = PPD%gaenes(jorb)
               PPD%iproj_to_l(iproj)   = iorbto_l(jorb)

               istart_c=istart_c + (   nvctr_c    +   7*nvctr_f   )
            endif
            jorb=jorb+1

            if(jorb> PPD%G%ncoeff) exit
         end do


         if( .not. PPD%DistProjApply) then
            istart_c= istart_c-mproj*(nvctr_c+7*nvctr_f)

            call fillPcProjOnTheFly(PPD, Glr, iat, at, hx,hy,hz, startjorb,ecut_pc ,  istart_c ) 
            istart_c= istart_c+mproj*(nvctr_c+7*nvctr_f)

            !!$
            !!$           ncplx=1
            !!$           rdum=0.0_gp
            !!$
            !!$           mbvctr_c=PPD%pc_nlpspd%nvctr_p(2*iat-1)-PPD%pc_nlpspd%nvctr_p(2*iat-2)
            !!$           mbvctr_f=PPD%pc_nlpspd%nvctr_p(2*iat  )-PPD%pc_nlpspd%nvctr_p(2*iat-1)
            !!$           
            !!$           mbseg_c=PPD%pc_nlpspd%nseg_p(2*iat-1)-PPD%pc_nlpspd%nseg_p(2*iat-2)
            !!$           mbseg_f=PPD%pc_nlpspd%nseg_p(2*iat  )-PPD%pc_nlpspd%nseg_p(2*iat-1)
            !!$           jseg_c=PPD%pc_nlpspd%nseg_p(2*iat-2)+1
            !!$              
            !!$           do idum=1, 9
            !!$              call wpdot_wrap(ncplx,  &
            !!$                   mbvctr_c,mbvctr_f,mbseg_c,mbseg_f,PPD%pc_nlpspd%keyv_p(jseg_c),&
            !!$                   PPD%pc_nlpspd%keyg_p(1,jseg_c),PPD%pc_proj(istart_c-idum*(nvctr_c+7*nvctr_f)),& 
            !!$                   mbvctr_c,mbvctr_f,mbseg_c,mbseg_f,PPD%pc_nlpspd%keyv_p(jseg_c),&
            !!$                   PPD%pc_nlpspd%keyg_p(1,jseg_c),&
            !!$                   PPD%pc_proj(istart_c-idum*(nvctr_c+7*nvctr_f)),&
            !!$                   rdum)
            !!$           end do
         endif

      end if

      !! aggiunger condizione su istartc_c per vedere se e nelprj

      startjorb=jorb

   enddo

   if( .not. PPD%DistProjApply) then
      call deallocate_gwf(PPD%G,subname)
   endif

   i_all=-product(shape(logrid))*kind(logrid)
   deallocate(logrid,stat=i_stat)
   call memocc(i_stat,i_all,'logrid',subname)

   i_all=-product(shape(Gocc))*kind(Gocc)
   deallocate(Gocc,stat=i_stat)
   call memocc(i_stat,i_all,'Gocc',subname)

   !!$  i_all=-product(shape(iorbtolr))*kind(iorbtolr)
   !!$  deallocate(iorbtolr,stat=i_stat)
   !!$  call memocc(i_stat,i_all,'iorbtolr',subname)

   i_all=-product(shape(iorbto_l))*kind(iorbto_l)
   deallocate(iorbto_l,stat=i_stat)
   call memocc(i_stat,i_all,'iorbto_l',subname)

   i_all=-product(shape(iorbto_m))*kind(iorbto_m)
   deallocate(iorbto_m,stat=i_stat)
   call memocc(i_stat,i_all,'iorbto_m',subname)

   i_all=-product(shape(iorbto_ishell))*kind(iorbto_ishell)
   deallocate(iorbto_ishell,stat=i_stat)
   call memocc(i_stat,i_all,'iorbto_ishell',subname)


   i_all=-product(shape(iorbto_iexpobeg))*kind(iorbto_iexpobeg)
   deallocate(iorbto_iexpobeg,stat=i_stat)
   call memocc(i_stat,i_all,'iorbto_iexpobeg',subname)

END SUBROUTINE createPcProjectorsArrays


!> Determine localization region for all preconditioning projectors, but do not yet fill the descriptor arrays
subroutine createPawProjectorsArrays(iproc,n1,n2,n3,rxyz,at,orbs,&
      &   radii_cf,cpmult,fpmult,hx,hy,hz, &
      &   PAWD, Glr)
   use module_interfaces
   use module_base
   use module_types
   implicit none
   integer, intent(in) :: iproc,n1,n2,n3
   real(gp), intent(in) :: cpmult,fpmult,hx,hy,hz
   type(atoms_data), intent(in) :: at
   type(orbitals_data), intent(in) :: orbs

   real(gp), dimension(3,at%astruct%nat), intent(in) :: rxyz
   real(gp), dimension(at%astruct%ntypes,3), intent(in) :: radii_cf

   type(PAWproj_data_type) ::PAWD

   type(locreg_descriptors),  intent(in):: Glr

   !local variables
   character(len=*), parameter :: subname='createPawProjectorsArrays'

   integer :: nl1,nl2,nl3,nu1,nu2,nu3,mseg,mproj, mvctr
   integer :: iat,i_stat,i_all,iseg, istart_c
   logical, dimension(:,:,:), allocatable :: logrid

   real(wp), dimension(:), pointer :: Gocc

   integer, pointer :: iorbto_l(:)
   integer, pointer :: iorbto_paw_nchannels(:)
   integer, pointer :: iorbto_m(:)
   integer, pointer :: iorbto_ishell(:)
   integer, pointer :: iorbto_iexpobeg(:)

   integer :: ncplx
   real(gp) :: kx,ky,kz
   integer ::  jorb
   integer :: iproj, startjorb
   real(gp) :: Pcpmult
   integer :: nvctr_c, nvctr_f
   integer :: iatat

   integer :: ikpt,iskpt,iekpt

   Pcpmult=1.0*cpmult


   nullify(PAWD%G%rxyz)

   call gaussian_pswf_basis_for_paw(at,rxyz,PAWD%G, &
      &   PAWD%iorbtolr,iorbto_l, iorbto_m,  iorbto_ishell,iorbto_iexpobeg  ,&
      &   iorbto_paw_nchannels, PAWD%iprojto_imatrixbeg )  


   allocate(Gocc(PAWD%G%ncoeff+ndebug),stat=i_stat)
   call memocc(i_stat,Gocc,'Gocc',subname)
   call razero(PAWD%G%ncoeff,Gocc)

   ! allocated  : gaenes, Gocc , PAWD%iorbtolr,iorbto_l, iorbto_m,  iorbto_ishell,iorbto_iexpobeg, iorbto_paw_nchannels

   !!$ ========================================================================================
   !---------

   PAWD%paw_nlpspd%natoms=PAWD%G%nat
   allocate(PAWD%paw_nlpspd%plr(PAWD%paw_nlpspd%natoms))


!!$   allocate(PAWD%paw_nlpspd%nseg_p(0:2*PAWD%G%nat+ndebug),stat=i_stat)
!!$   call memocc(i_stat,PAWD%paw_nlpspd%nseg_p,'pc_nlpspd%nseg_p',subname)
!!$   allocate(PAWD%paw_nlpspd%nvctr_p(0:2*PAWD%G%nat+ndebug),stat=i_stat)
!!$   call memocc(i_stat,PAWD%paw_nlpspd%nvctr_p,'pc_nlpspd%nvctr_p',subname)
!!$   allocate(PAWD%paw_nlpspd%nboxp_c(2,3,PAWD%G%nat+ndebug),stat=i_stat)
!!$   call memocc(i_stat,PAWD%paw_nlpspd%nboxp_c,'pc_nlpspd%nboxp_c',subname)
!!$   allocate(PAWD%paw_nlpspd%nboxp_f(2,3,PAWD%G%nat+ndebug),stat=i_stat)
!!$   call memocc(i_stat,PAWD%paw_nlpspd%nboxp_f,'pc_nlpspd%nboxp_f',subname)

   allocate(logrid(0:n1,0:n2,0:n3+ndebug),stat=i_stat)
   call memocc(i_stat,logrid,'logrid',subname)

   call localize_projectors_paw(iproc,n1,n2,n3,hx,hy,hz,Pcpmult,1*fpmult,rxyz,radii_cf,&
      &   logrid,at,orbs,PAWD)

   ! the above routine counts atomic projector and the number of their element for psp
   ! We must therefore correct , later, nlpspd%nprojel  and nlpspd%nproj
   !-------------------

   ! allocations for arrays holding the projectors and their data descriptors
   do iat=1,PAWD%paw_nlpspd%natoms
      !for the moments the bounds are not needed for projectors
      call allocate_wfd(PAWD%paw_nlpspd%plr(iat)%wfd,subname)
   end do

!!$   allocate(PAWD%paw_nlpspd%keyg_p(2,PAWD%paw_nlpspd%nseg_p(2*PAWD%G%nat)+ndebug),stat=i_stat)
!!$   call memocc(i_stat,PAWD%paw_nlpspd%keyg_p,'pc_nlpspd%keyg_p',subname)
!!$
!!$   allocate(PAWD%paw_nlpspd%keyv_p(PAWD%paw_nlpspd%nseg_p(2*PAWD%G%nat)+ndebug),stat=i_stat)
!!$   call memocc(i_stat,PAWD%paw_nlpspd%keyv_p,'pc_nlpspd%keyv_p',subname)

   !!$  -- this one delayed, waiting for the correct pc_nlpspd%nprojel, pc_nlpspd%nproj
   !!$  --
   !!$  allocate(pc_proj(pc_nlpspd%nprojel+ndebug),stat=i_stat)
   !!$  call memocc(i_stat,pc_proj,'pc_proj',subname)
   allocate(PAWD%paw_proj(PAWD%paw_nlpspd%nprojel+ndebug),stat=i_stat)
   call memocc(i_stat,PAWD%paw_proj,'paw_proj',subname)

   allocate(PAWD%ilr_to_mproj(PAWD%G%nat  +ndebug ) , stat=i_stat)
   call memocc(i_stat,PAWD%ilr_to_mproj,'ilr_to_mproj',subname)

   allocate(PAWD%iproj_to_l(PAWD%paw_nlpspd%nproj +ndebug ) , stat=i_stat)
   call memocc(i_stat ,PAWD%iproj_to_l,'iproj_to_l',subname)

   allocate(PAWD%iproj_to_paw_nchannels( PAWD%paw_nlpspd%nproj+ndebug ) , stat=i_stat)
   call memocc(i_stat ,PAWD%iproj_to_paw_nchannels,'iproj_to_paw_nchannels',subname)

   !!$ =========================================================================================  

   jorb=1  
   ! After having determined the size of the projector descriptor arrays fill them
   iat=0
   do iatat=1, at%astruct%nat
      if (  at%paw_NofL(at%astruct%iatype(iatat)).gt.0  ) then
         iat=iat+1
         mproj=0
         do while( jorb<=PAWD%G%ncoeff         .and. PAWD%iorbtolr(jorb)== iat)
            mproj=mproj+1
            if(jorb==PAWD%G%ncoeff) exit
            jorb=jorb+1
         end do

         PAWD%paw_nlpspd%nproj=PAWD%paw_nlpspd%nproj+mproj
         if (mproj.ne.0) then 

            call bounds_to_plr_limits(.false.,1,PAWD%paw_nlpspd%plr(iat),nl1,nl2,nl3,nu1,nu2,nu3)
!!$            ! coarse grid quantities
!!$            nl1=PAWD%paw_nlpspd%nboxp_c(1,1,iat) 
!!$            nl2=PAWD%paw_nlpspd%nboxp_c(1,2,iat) 
!!$            nl3=PAWD%paw_nlpspd%nboxp_c(1,3,iat) 
!!$
!!$            nu1=PAWD%paw_nlpspd%nboxp_c(2,1,iat)
!!$            nu2=PAWD%paw_nlpspd%nboxp_c(2,2,iat)
!!$            nu3=PAWD%paw_nlpspd%nboxp_c(2,3,iat)

            call fill_logrid(at%astruct%geocode,n1,n2,n3,nl1,nu1,nl2,nu2,nl3,nu3,0,1,  &
               &   at%astruct%ntypes,at%astruct%iatype(iatat),rxyz(1,iatat),radii_cf(1,3),Pcpmult,hx,hy,hz,logrid)

!!$            iseg=PAWD%paw_nlpspd%nseg_p(2*iat-2)+1
!!$            mseg=PAWD%paw_nlpspd%nseg_p(2*iat-1)-PAWD%paw_nlpspd%nseg_p(2*iat-2)
            mseg=PAWD%paw_nlpspd%plr(iat)%wfd%nseg_c

            call segkeys(n1,n2,n3,nl1,nu1,nl2,nu2,nl3,nu3,  & 
                 logrid,mseg,&
!!$                 PAWD%paw_nlpspd%keyg_p(1,iseg),PAWD%paw_nlpspd%keyv_p(iseg))
                 PAWD%paw_nlpspd%plr(iat)%wfd%keyglob(1,1),&
                 PAWD%paw_nlpspd%plr(iat)%wfd%keyvglob(1))

            mvctr =PAWD%paw_nlpspd%plr(iat)%wfd%nvctr_c!PAWD%paw_nlpspd%nvctr_p(2*iat-1)-PAWD%paw_nlpspd%nvctr_p(2*iat-2)

            call bounds_to_plr_limits(.false.,2,PAWD%paw_nlpspd%plr(iat),&
                 nl1,nl2,nl3,nu1,nu2,nu3)
!!$            ! fine grid quantities
!!$            nl1=PAWD%paw_nlpspd%nboxp_f(1,1,iat)
!!$            nl2=PAWD%paw_nlpspd%nboxp_f(1,2,iat)
!!$            nl3=PAWD%paw_nlpspd%nboxp_f(1,3,iat)
!!$
!!$            nu1=PAWD%paw_nlpspd%nboxp_f(2,1,iat)
!!$            nu2=PAWD%paw_nlpspd%nboxp_f(2,2,iat)
!!$            nu3=PAWD%paw_nlpspd%nboxp_f(2,3,iat)
            call fill_logrid(at%astruct%geocode,n1,n2,n3,nl1,nu1,nl2,nu2,nl3,nu3,0,1,  &
               &   at%astruct%ntypes,at%astruct%iatype(iatat),rxyz(1,iatat),radii_cf(1,2),1*fpmult,hx,hy,hz,logrid)
            
!!$            iseg=PAWD%paw_nlpspd%nseg_p(2*iat-1)+1
!!$            mseg=PAWD%paw_nlpspd%nseg_p(2*iat)-PAWD%paw_nlpspd%nseg_p(2*iat-1)
            iseg=PAWD%paw_nlpspd%plr(iat)%wfd%nseg_c+1
            mseg=PAWD%paw_nlpspd%plr(iat)%wfd%nseg_f

            if (mseg > 0) then
               call segkeys(n1,n2,n3,nl1,nu1,nl2,nu2,nl3,nu3,  & 
                    logrid,mseg,&
!!$               PAWD%paw_nlpspd%keyg_p(1,iseg),PAWD%paw_nlpspd%keyv_p(iseg))
                    PAWD%paw_nlpspd%plr(iat)%wfd%keyglob(1,iseg),&
                    PAWD%paw_nlpspd%plr(iat)%wfd%keyvglob(iseg))
               
               mvctr =PAWD%paw_nlpspd%plr(iat)%wfd%nvctr_f!PAWD%paw_nlpspd%nvctr_p(2*iat)-PAWD%paw_nlpspd%nvctr_p(2*iat-1)
            end if


         endif
      endif
   enddo

   if (orbs%norbp > 0) then
      iskpt=orbs%iokpt(1)
      iekpt=orbs%iokpt(orbs%norbp)
   else
      iskpt=1
      iekpt=1
   end if

   istart_c=1
   do ikpt=iskpt,iekpt     

      !features of the k-point ikpt
      kx=orbs%kpts(1,ikpt)
      ky=orbs%kpts(2,ikpt)
      kz=orbs%kpts(3,ikpt)
      !!  write( *, '(A,i4,1x,A,3(1x,d20.10))') " IKPT , " , ikpt, " K " , orbs%kpts(:,ikpt)
      !evaluate the complexity of the k-point
      if (kx**2 + ky**2 + kz**2 == 0.0_gp) then
         ncplx=1
      else
         ncplx=2
      end if

      startjorb=1
      jorb=1
      Gocc(:)=0.0_wp
      iproj=0

      iat=0
      do iatat=1, at%astruct%nat
         if (  at%paw_NofL(at%astruct%iatype(iatat)).gt.0  ) then
            iat=iat+1
            mproj=0
            do while( jorb<=PAWD%G%ncoeff         .and. PAWD%iorbtolr(jorb)== iat)
               mproj=mproj+1
               if(jorb==PAWD%G%ncoeff) exit
               jorb=jorb+1
            end do

            PAWD%ilr_to_mproj(iat)=mproj
            if( mproj>0) then
               nvctr_c  =PAWD%paw_nlpspd%plr(iat)%wfd%nvctr_c!PAWD%paw_nlpspd%nvctr_p(2*iat-1)-PAWD%paw_nlpspd%nvctr_p(2*iat-2)
               nvctr_f  =PAWD%paw_nlpspd%plr(iat)%wfd%nvctr_f!PAWD%paw_nlpspd%nvctr_p(2*iat  )-PAWD%paw_nlpspd%nvctr_p(2*iat-1)

               jorb=startjorb
               do while( jorb<=PAWD%G%ncoeff  .and. PAWD%iorbtolr(jorb)== iat) 
                  iproj=iproj+1
                  PAWD%iproj_to_l(iproj)   = iorbto_l(jorb)
                  PAWD%iproj_to_paw_nchannels(iproj)   = iorbto_paw_nchannels(jorb)
                  istart_c=istart_c + (   nvctr_c    +   7*nvctr_f   )*ncplx
                  jorb=jorb+1
                  if(jorb> PAWD%G%ncoeff) exit
               end do
               if( .not. PAWD%DistProjApply) then
                  istart_c= istart_c-mproj*(nvctr_c+7*nvctr_f)*ncplx
                  call fillPawProjOnTheFly(PAWD, Glr, iat,  hx,hy,hz, kx,ky,kz, startjorb,&
                     &   istart_c, at%astruct%geocode , at, iatat) 
                  istart_c= istart_c+mproj*(nvctr_c+7*nvctr_f)*ncplx
               endif
            end if
            startjorb=jorb
         end if
      enddo
   enddo
   if (istart_c-1 /= PAWD%paw_nlpspd%nprojel) stop 'incorrect once-and-for-all psp generation'


   if( .not. PAWD%DistProjApply) then
      call deallocate_gwf_c(PAWD%G,subname)
   endif

   i_all=-product(shape(logrid))*kind(logrid)
   deallocate(logrid,stat=i_stat)
   call memocc(i_stat,i_all,'logrid',subname)

   i_all=-product(shape(Gocc))*kind(Gocc)
   deallocate(Gocc,stat=i_stat)
   call memocc(i_stat,i_all,'Gocc',subname)

   !!$  i_all=-product(shape(iorbtolr))*kind(iorbtolr)
   !!$  deallocate(iorbtolr,stat=i_stat)
   !!$  call memocc(i_stat,i_all,'iorbtolr',subname)

   i_all=-product(shape(iorbto_l))*kind(iorbto_l)
   deallocate(iorbto_l,stat=i_stat)
   call memocc(i_stat,i_all,'iorbto_l',subname)

   i_all=-product(shape(iorbto_paw_nchannels))*kind(iorbto_paw_nchannels)
   deallocate(iorbto_paw_nchannels,stat=i_stat)
   call memocc(i_stat,i_all,'iorbto_paw_nchannels',subname)





   i_all=-product(shape(iorbto_m))*kind(iorbto_m)
   deallocate(iorbto_m,stat=i_stat)
   call memocc(i_stat,i_all,'iorbto_m',subname)

   i_all=-product(shape(iorbto_ishell))*kind(iorbto_ishell)
   deallocate(iorbto_ishell,stat=i_stat)
   call memocc(i_stat,i_all,'iorbto_ishell',subname)


   i_all=-product(shape(iorbto_iexpobeg))*kind(iorbto_iexpobeg)
   deallocate(iorbto_iexpobeg,stat=i_stat)
   call memocc(i_stat,i_all,'iorbto_iexpobeg',subname)


END SUBROUTINE createPawProjectorsArrays

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

  call deallocate_wfd(wfd_old,subname)

  i_all=-product(shape(psi_old))*kind(psi_old)
  deallocate(psi_old,stat=i_stat)
  call memocc(i_stat,i_all,'psi_old',subname)
END SUBROUTINE input_wf_memory



subroutine input_memory_linear(iproc, nproc, at, KSwfn, tmb, tmb_old, denspot, input, &
           rxyz_old, rxyz, denspot0, energs, nlpspd, proj, GPU)

  use module_base
  use module_types
  use module_interfaces, except_this_one => input_memory_linear
  use module_fragments
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
  type(nonlocal_psp_descriptors), intent(in) :: nlpspd
  real(kind=8), dimension(:), pointer :: proj
  type(GPU_pointers), intent(inout) :: GPU

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
         allocate(frag_trans(iorb)%discrete_operations(0),stat=i_stat)
         call memocc(i_stat,frag_trans(iorb)%discrete_operations,'frag_trans(iorb)%discrete_operations',subname)
     end do

     call reformat_supportfunctions(iproc,at,rxyz_old,rxyz,.true.,tmb,ndim_old,tmb_old%lzd,frag_trans,tmb_old%psi)

     do iorb=1,tmb%orbs%norbp
        i_all = -product(shape(frag_trans(iorb)%discrete_operations)*kind(frag_trans(iorb)%discrete_operations))
        deallocate(frag_trans(iorb)%discrete_operations,stat=i_stat)
        call memocc(i_stat,i_all,'frag_trans(iorb)%discrete_operations',subname)
     end do
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
  if (associated(tmb_old%linmat%denskern%matrix_compr)) then
     i_all=-product(shape(tmb_old%linmat%denskern%matrix_compr))*kind(tmb_old%linmat%denskern%matrix_compr)
     deallocate(tmb_old%linmat%denskern%matrix_compr, stat=i_stat)
     call memocc(i_stat, i_all, 'tmb_old%linmat%denskern%matrix_compr', subname)
  end if

  ! destroy it all together here - don't have all comms arrays
  !call destroy_DFT_wavefunction(tmb_old)

  !!write(*,*) 'after deallocate, iproc', iproc

   ! normalize tmbs - only really needs doing if we reformatted, but will need to calculate transpose after anyway

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


  ! Update the kernel
  if (input%lin%scf_mode/=LINEAR_FOE) then
      !tmb%can_use_transposed=.false.
      !overlap_calculated = .false.
      !nullify(tmb%psit_c)
      !nullify(tmb%psit_f)
      call reconstruct_kernel(iproc, nproc, 0, tmb%orthpar%blocksize_pdsyev, tmb%orthpar%blocksize_pdgemm, &
           KSwfn%orbs, tmb, overlap_calculated)
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
          rxyz,nlpspd,proj,GPU,KSwfn%orbs,tmb,denspot,denspot0,energs)
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
  !!  end do
  !!end if

  ! Must initialize rhopotold (FOR NOW... use the trivial one)
  if (input%lin%scf_mode/=LINEAR_FOE) then
      call communicate_basis_for_density_collective(iproc, nproc, tmb%lzd, max(tmb%npsidim_orbs,tmb%npsidim_comp), &
           tmb%orbs, tmb%psi, tmb%collcom_sr)
      call sumrho_for_TMBs(iproc, nproc, KSwfn%Lzd%hgrids(1), KSwfn%Lzd%hgrids(2), KSwfn%Lzd%hgrids(3), &
           tmb%collcom_sr, tmb%linmat%denskern, KSwfn%Lzd%Glr%d%n1i*KSwfn%Lzd%Glr%d%n2i*denspot%dpbox%n3d, denspot%rhov)
      call dcopy(max(denspot%dpbox%ndims(1)*denspot%dpbox%ndims(2)*denspot%dpbox%n3p,1)*input%nspin, &
           denspot%rhov(1), 1, denspot0(1), 1)
      call updatePotential(input%ixc,input%nspin,denspot,energs%eh,energs%exc,energs%evxc)
      call local_potential_dimensions(tmb%lzd,tmb%orbs,denspot%dpbox%ngatherarr(0,1))
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
     nlpspd,proj,ixc,psi,hpsi,psit,G,&
     nspin,symObj,GPU,input,onlywf,proj_G,paw)
   ! Input wavefunctions are found by a diagonalization in a minimal basis set
   ! Each processors write its initial wavefunctions into the wavefunction file
   ! The files are then read by readwave
   ! @todo pass GPU to be a local variable of this routine (initialized and freed here)
   use module_base
   use module_interfaces, except_this_one => input_wf_diag
   use module_types
   use Poisson_Solver, except_dp => dp, except_gp => gp, except_wp => wp
   use libxc_functionals
   use yaml_output
   use gaussians
   implicit none
   !Arguments
   integer, intent(in) :: iproc,nproc,ixc
   integer, intent(inout) :: nspin,nvirt
   logical, intent(in) :: onlywf  !if .true. finds only the WaveFunctions and return
   type(atoms_data), intent(in) :: at
   type(nonlocal_psp_descriptors), intent(in) :: nlpspd
   type(local_zone_descriptors), intent(in) :: Lzd
   type(communications_arrays), intent(in) :: comms
   type(energy_terms), intent(inout) :: energs
   type(orbitals_data), intent(inout) :: orbs
   type(DFT_local_fields), intent(inout) :: denspot
   type(GPU_pointers), intent(in) :: GPU
   type(input_variables), intent(in) :: input
   type(symmetry_data), intent(in) :: symObj
   real(gp), dimension(3,at%astruct%nat), intent(in) :: rxyz
   real(wp), dimension(nlpspd%nprojel), intent(inout) :: proj
   type(gaussian_basis), intent(out) :: G !basis for davidson IG
   type(gaussian_basis),dimension(at%astruct%ntypes),optional,intent(in)::proj_G
   type(paw_objects),optional,intent(inout)::paw
   real(wp), dimension(:), pointer :: psi,hpsi,psit
   !local variables
   character(len=*), parameter :: subname='input_wf_diag'
   logical :: switchGPUconv,switchOCLconv
   integer :: ii,jj
   integer :: i_stat,i_all,nspin_ig,ncplx,irhotot_add,irho_add,ispin
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
        orbs,orbse,norbsc_arr,locrad,G,psigau,eks)

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
        proj,Lzde,nlpspd,confdatarr,denspot%dpbox%ngatherarr,denspot%pot_work,psi,hpsi,&
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
    if (.false.) then
       call DiagHam(iproc,nproc,at%natsc,nspin_ig,orbs,Lzde%Glr%wfd,comms,&
          psi,hpsi,psit,input%orthpar,passmat,orbse,commse,etol,norbsc_arr)
    end if

    if (iproc==0) call yaml_newline()

   !test merging of Linear and cubic
     call LDiagHam(iproc,nproc,at%natsc,nspin_ig,orbs,Lzd,Lzde,comms,&
         psi,hpsi,psit,input%orthpar,passmat,input%iscf,input%Tel,input%occopt,&
         orbse,commse,etol,norbsc_arr)

     i_all=-product(shape(passmat))*kind(passmat)
     deallocate(passmat,stat=i_stat)
     call memocc(i_stat,i_all,'passmat',subname)

   if (input%iscf > SCF_KIND_DIRECT_MINIMIZATION .or. input%Tel > 0.0_gp) then

!commented out, this part has already been done in LDiagHam     
!!$      !clean the array of the IG eigenvalues
!!$      call to_zero(orbse%norb*orbse%nkpts,orbse%eval(1))
!!$      !put the actual values on it
!!$      call dcopy(orbs%norb*orbs%nkpts,orbs%eval(1),1,orbse%eval(1),1)
!!$
!!$      !add a small displacement in the eigenvalues
!!$      do iorb=1,orbs%norb*orbs%nkpts
!!$         tt=builtin_rand(idum)
!!$         orbs%eval(iorb)=orbs%eval(iorb)*(1.0_gp+max(input%Tel,1.0e-3_gp)*real(tt,gp))
!!$      end do
!!$
!!$      !correct the occupation numbers wrt fermi level
!!$      call evaltoocc(iproc,nproc,.false.,input%Tel,orbs,input%occopt)

      !restore the occupations 
      call dcopy(orbs%norb*orbs%nkpts,orbse%occup(1),1,orbs%occup(1),1)
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
!!$      call write_eigenvalues_data(nproc,orbse,mom_vec)
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
     denspot,denspot0,nlpspd,proj,KSwfn,tmb,energs,inputpsi,input_wf_format,norbv,&
     wfd_old,psi_old,d_old,hx_old,hy_old,hz_old,rxyz_old,tmb_old,ref_frags,cdft)
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
  type(nonlocal_psp_descriptors), intent(in) :: nlpspd
  real(kind=8), dimension(:), pointer :: proj
  !type(gaussian_basis), intent(inout) :: gbd
  !real(wp), dimension(:,:), pointer :: gaucoeffs
  type(grid_dimensions), intent(in) :: d_old
  real(gp), dimension(3, atoms%astruct%nat), intent(inout) :: rxyz_old
  type(wavefunctions_descriptors), intent(inout) :: wfd_old
  type(system_fragment), dimension(:), pointer :: ref_frags
  type(cdft_data), intent(out) :: cdft
  !local variables
  character(len = *), parameter :: subname = "input_wf"
  integer :: i_stat, nspin, i_all, ifrag, iorb, itmb, jtmb, ierr
  type(gaussian_basis) :: Gvirt
  real(wp), allocatable, dimension(:) :: norm
  !wvl+PAW objects
  integer :: iatyp
  type(gaussian_basis),dimension(atoms%astruct%ntypes)::proj_G
  type(paw_objects)::paw
  logical :: overlap_calculated
  !!real(gp), dimension(:,:), allocatable :: ks, ksk
  !!real(gp) :: nonidem

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
        call yaml_comment('Random wavefunctions initialization',hfill='-')
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
        call yaml_comment('Wavefunctions from PSP Atomic Orbitals initialization',hfill='-')
        call yaml_open_map('Input Hamiltonian')
     end if
     nspin=in%nspin
     !calculate input guess from diagonalisation of LCAO basis (written in wavelets)
     call input_wf_diag(iproc,nproc, atoms,denspot,&
          KSwfn%orbs,norbv,KSwfn%comms,KSwfn%Lzd,energs,rxyz,&
          nlpspd,proj,in%ixc,KSwfn%psi,KSwfn%hpsi,KSwfn%psit,&
          Gvirt,nspin,atoms%astruct%sym,GPU,in,.false.)

  case(INPUT_PSI_MEMORY_WVL)
     !restart from previously calculated wavefunctions, in memory
     if (iproc == 0) then
        !write( *,'(1x,a)')&
        !     &   '-------------------------------------------------------------- Wavefunctions Restart'
        call yaml_comment('Wavefunctions Restart',hfill='-')
        call yaml_open_map("Input Hamiltonian")
     end if
     if (in%wfn_history <= 2) then
        call input_wf_memory(iproc, atoms, &
             rxyz_old, hx_old, hy_old, hz_old, d_old, wfd_old, psi_old, &
             rxyz,KSwfn%Lzd%hgrids(1),KSwfn%Lzd%hgrids(2),KSwfn%Lzd%hgrids(3),&
             KSwfn%Lzd%Glr%d,KSwfn%Lzd%Glr%wfd,KSwfn%psi, KSwfn%orbs)
     else
        call input_wf_memory_history(iproc,KSwfn%orbs,atoms,in%wfn_history,&
             Kswfn%istep_history,KSwfn%oldpsis,rxyz,Kswfn%Lzd,KSwfn%psi)
     end if

     if (in%iscf > SCF_KIND_DIRECT_MINIMIZATION) &
          call evaltoocc(iproc,nproc,.false.,in%Tel,KSwfn%orbs,in%occopt)
  case(INPUT_PSI_MEMORY_LINEAR)
     if (iproc == 0) then
        call yaml_comment('Support functions Restart',hfill='-')
        call yaml_open_map("Input Hamiltonian")
     end if
      call input_memory_linear(iproc, nproc, atoms, KSwfn, tmb, tmb_old, denspot, in, &
           rxyz_old, rxyz, denspot0, energs, nlpspd, proj, GPU)
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
     call inputguessConfinement(iproc,nproc,atoms,in,KSwfn%Lzd%hgrids(1),KSwfn%Lzd%hgrids(2),KSwfn%Lzd%hgrids(3), &
          rxyz,nlpspd,proj,GPU,KSwfn%orbs,tmb,denspot,denspot0,energs)
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

     ! we have to copy the coeffs from the fragment structure to the tmb structure and reconstruct each 'mini' kernel
     ! this is overkill as we are recalculating the kernel anyway - fix at some point
     ! or just put into fragment structure to save recalculating for CDFT
     if (in%lin%fragment_calculation) then
        call fragment_coeffs_to_kernel(in%frag,ref_frags,tmb,KSwfn%orbs,overlap_calculated)
     end if

     call reconstruct_kernel(iproc, nproc, 0, tmb%orthpar%blocksize_pdsyev, tmb%orthpar%blocksize_pdgemm, &
          KSwfn%orbs, tmb, overlap_calculated)     

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

     ! Now need to calculate the charge density and the potential related to this inputguess
     call communicate_basis_for_density_collective(iproc, nproc, tmb%lzd, max(tmb%npsidim_orbs,tmb%npsidim_comp), &
          tmb%orbs, tmb%psi, tmb%collcom_sr)

     call sumrho_for_TMBs(iproc, nproc, KSwfn%Lzd%hgrids(1), KSwfn%Lzd%hgrids(2), KSwfn%Lzd%hgrids(3), &
          tmb%collcom_sr, tmb%linmat%denskern, KSwfn%Lzd%Glr%d%n1i*KSwfn%Lzd%Glr%d%n2i*denspot%dpbox%n3d, denspot%rhov)

     ! CDFT: calculate w(r), define some initial guess for V and initialize other cdft_data stuff
     ! CDFT: could move to linearscaling but use rhopotold rather than rhov, but keep here in case we want to deallocate fragment stuff

     call timing(iproc,'constraineddft','ON')
     call nullify_cdft_data(cdft)
     if (in%lin%constrained_dft) then
        call cdft_data_init(cdft,in%frag,KSwfn%Lzd%Glr%d%n1i*KSwfn%Lzd%Glr%d%n2i*denspot%dpbox%n3d)
        call cdft_data_allocate(cdft,tmb%linmat%ham)
        if (trim(cdft%method)=='fragment_density') then
           call calculate_weight_function(in,ref_frags,cdft,&
                KSwfn%Lzd%Glr%d%n1i*KSwfn%Lzd%Glr%d%n2i*denspot%dpbox%n3d,denspot%rhov,tmb,atoms,rxyz,denspot)
        end if
     end if
     call timing(iproc,'constraineddft','OF')

     ! Must initialize rhopotold (FOR NOW... use the trivial one)
     call dcopy(max(denspot%dpbox%ndims(1)*denspot%dpbox%ndims(2)*denspot%dpbox%n3p,1)*in%nspin, &
          denspot%rhov(1), 1, denspot0(1), 1)
     !!call deallocateCommunicationbufferSumrho(tmb%comsr, subname)
     call updatePotential(in%ixc,in%nspin,denspot,energs%eh,energs%exc,energs%evxc)
     call local_potential_dimensions(tmb%lzd,tmb%orbs,denspot%dpbox%ngatherarr(0,1))

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
           call yaml_warning('Missing wavefunction files, switch to normal input guess')
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
