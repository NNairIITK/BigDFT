!!****f* BigDFT/createWavefunctionsDescriptors
!! FUNCTION
!!   Calculates the descriptor arrays and nvctrp
!!   Calculates also the bounds arrays needed for convolutions
!!   Refers this information to the global localisation region descriptor
!!
!! COPYRIGHT
!!    Copyright (C) 2007-2010 bigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 
!!
!! SOURCE
!!
subroutine createWavefunctionsDescriptors(iproc,hx,hy,hz,atoms,rxyz,radii_cf,&
     crmult,frmult,Glr)
  use module_base
  use module_types
  implicit none
  !Arguments
  type(atoms_data), intent(in) :: atoms
  integer, intent(in) :: iproc
  real(gp), intent(in) :: hx,hy,hz,crmult,frmult
  real(gp), dimension(3,atoms%nat), intent(in) :: rxyz
  real(gp), dimension(atoms%ntypes,3), intent(in) :: radii_cf
  type(locreg_descriptors), intent(inout) :: Glr
  !local variables
  character(len=*), parameter :: subname='createWavefunctionsDescriptors'
  integer :: i_all,i_stat
  integer :: n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3
  logical, dimension(:,:,:), allocatable :: logrid_c,logrid_f

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
  if (atoms%geocode == 'F') then
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

  if (iproc == 0) then
     write(*,'(1x,a)')&
          '------------------------------------------------- Wavefunctions Descriptors Creation'
  end if

  ! determine localization region for all orbitals, but do not yet fill the descriptor arrays
  allocate(logrid_c(0:n1,0:n2,0:n3+ndebug),stat=i_stat)
  call memocc(i_stat,logrid_c,'logrid_c',subname)
  allocate(logrid_f(0:n1,0:n2,0:n3+ndebug),stat=i_stat)
  call memocc(i_stat,logrid_f,'logrid_f',subname)

  ! coarse grid quantities
  call fill_logrid(atoms%geocode,n1,n2,n3,0,n1,0,n2,0,n3,0,atoms%nat,&
       atoms%ntypes,atoms%iatype,rxyz,radii_cf(1,1),crmult,hx,hy,hz,logrid_c)
  call num_segkeys(n1,n2,n3,0,n1,0,n2,0,n3,logrid_c,Glr%wfd%nseg_c,Glr%wfd%nvctr_c)
  if (iproc == 0) write(*,'(2(1x,a,i10))') &
       'Coarse resolution grid: Number of segments= ',Glr%wfd%nseg_c,'points=',Glr%wfd%nvctr_c

  if (atoms%geocode == 'F') then
     call make_bounds(n1,n2,n3,logrid_c,Glr%bounds%kb%ibyz_c,Glr%bounds%kb%ibxz_c,Glr%bounds%kb%ibxy_c)
  end if

  if (atoms%geocode == 'P' .and. .not. Glr%hybrid_on .and. Glr%wfd%nvctr_c /= (n1+1)*(n2+1)*(n3+1) ) then
     if (iproc ==0)then
        write(*,*)&
          ' ERROR: the coarse grid does not fill the entire periodic box'
        write(*,*)&
          '          errors due to translational invariance breaking may occur'
        !stop
     end if
     if (GPUconv) then
!        if (iproc ==0)then
           write(*,*)&
                '          The code should be stopped for a GPU calculation     '
           write(*,*)&
                '          since density is not initialised to 10^-20               '
!        end if
        stop
     end if
  end if

  call fill_logrid(atoms%geocode,n1,n2,n3,0,n1,0,n2,0,n3,0,atoms%nat,&
       atoms%ntypes,atoms%iatype,rxyz,radii_cf(1,2),frmult,hx,hy,hz,logrid_f)
  call num_segkeys(n1,n2,n3,0,n1,0,n2,0,n3,logrid_f,Glr%wfd%nseg_f,Glr%wfd%nvctr_f)
  if (iproc == 0) write(*,'(2(1x,a,i10))') & 
       '  Fine resolution grid: Number of segments= ',Glr%wfd%nseg_f,'points=',Glr%wfd%nvctr_f
  if (atoms%geocode == 'F') then
     call make_bounds(n1,n2,n3,logrid_f,Glr%bounds%kb%ibyz_f,Glr%bounds%kb%ibxz_f,Glr%bounds%kb%ibxy_f)
  end if

  ! allocations for arrays holding the wavefunctions and their data descriptors
  call allocate_wfd(Glr%wfd,subname)

  ! now fill the wavefunction descriptor arrays
  ! coarse grid quantities
  call segkeys(n1,n2,n3,0,n1,0,n2,0,n3,logrid_c,Glr%wfd%nseg_c,Glr%wfd%keyg(1,1),Glr%wfd%keyv(1))

  ! fine grid quantities
  if (Glr%wfd%nseg_f > 0) then
     call segkeys(n1,n2,n3,0,n1,0,n2,0,n3,logrid_f,Glr%wfd%nseg_f,Glr%wfd%keyg(1,Glr%wfd%nseg_c+1), &
          & Glr%wfd%keyv(Glr%wfd%nseg_c+1))
  end if

  i_all=-product(shape(logrid_c))*kind(logrid_c)
  deallocate(logrid_c,stat=i_stat)
  call memocc(i_stat,i_all,'logrid_c',subname)
  i_all=-product(shape(logrid_f))*kind(logrid_f)
  deallocate(logrid_f,stat=i_stat)
  call memocc(i_stat,i_all,'logrid_f',subname)

  !for free BC admits the bounds arrays
  if (atoms%geocode == 'F') then

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
          Glr%bounds%kb%ibxy_c,Glr%bounds%sb%ibzzx_c,Glr%bounds%sb%ibyyzz_c,&
          Glr%bounds%kb%ibxy_f,Glr%bounds%sb%ibxy_ff,Glr%bounds%sb%ibzzx_f,Glr%bounds%sb%ibyyzz_f,&
          Glr%bounds%kb%ibyz_c,Glr%bounds%gb%ibzxx_c,Glr%bounds%gb%ibxxyy_c,&
          Glr%bounds%kb%ibyz_f,Glr%bounds%gb%ibyz_ff,Glr%bounds%gb%ibzxx_f,Glr%bounds%gb%ibxxyy_f,&
          Glr%bounds%ibyyzz_r)

  end if

  if ( atoms%geocode == 'P' .and. Glr%hybrid_on) then
     call make_bounds_per(n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,Glr%bounds,Glr%wfd)
     call make_all_ib_per(n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,&
          Glr%bounds%kb%ibxy_f,Glr%bounds%sb%ibxy_ff,Glr%bounds%sb%ibzzx_f,Glr%bounds%sb%ibyyzz_f,&
          Glr%bounds%kb%ibyz_f,Glr%bounds%gb%ibyz_ff,Glr%bounds%gb%ibzxx_f,Glr%bounds%gb%ibxxyy_f)
  endif

  !assign geocode and the starting points
  Glr%geocode=atoms%geocode

END SUBROUTINE createWavefunctionsDescriptors
!!***


!!****f* BigDFT/createProjectorsArrays
!! FUNCTION
!!   Determine localization region for all projectors, but do not yet fill the descriptor arrays
!! SOURCE
!!
subroutine createProjectorsArrays(iproc,n1,n2,n3,rxyz,at,orbs,&
     radii_cf,cpmult,fpmult,hx,hy,hz,nlpspd,proj)
  use module_base
  use module_types
  implicit none
  integer, intent(in) :: iproc,n1,n2,n3
  real(gp), intent(in) :: cpmult,fpmult,hx,hy,hz
  type(atoms_data), intent(in) :: at
  type(orbitals_data), intent(in) :: orbs
  real(gp), dimension(3,at%nat), intent(in) :: rxyz
  real(gp), dimension(at%ntypes,3), intent(in) :: radii_cf
  type(nonlocal_psp_descriptors), intent(out) :: nlpspd
  real(wp), dimension(:), pointer :: proj
  !local variables
  character(len=*), parameter :: subname='createProjectorsArrays'
  integer :: nl1,nl2,nl3,nu1,nu2,nu3,mseg,mproj
  integer :: iat,i_stat,i_all,iseg
  logical, dimension(:,:,:), allocatable :: logrid
  
  allocate(nlpspd%nseg_p(0:2*at%nat+ndebug),stat=i_stat)
  call memocc(i_stat,nlpspd%nseg_p,'nlpspd%nseg_p',subname)
  allocate(nlpspd%nvctr_p(0:2*at%nat+ndebug),stat=i_stat)
  call memocc(i_stat,nlpspd%nvctr_p,'nlpspd%nvctr_p',subname)
  allocate(nlpspd%nboxp_c(2,3,at%nat+ndebug),stat=i_stat)
  call memocc(i_stat,nlpspd%nboxp_c,'nlpspd%nboxp_c',subname)
  allocate(nlpspd%nboxp_f(2,3,at%nat+ndebug),stat=i_stat)
  call memocc(i_stat,nlpspd%nboxp_f,'nlpspd%nboxp_f',subname)

  ! determine localization region for all projectors, but do not yet fill the descriptor arrays
  allocate(logrid(0:n1,0:n2,0:n3+ndebug),stat=i_stat)
  call memocc(i_stat,logrid,'logrid',subname)

  call localize_projectors(iproc,n1,n2,n3,hx,hy,hz,cpmult,fpmult,rxyz,radii_cf,&
       logrid,at,orbs,nlpspd)

  ! allocations for arrays holding the projectors and their data descriptors
  allocate(nlpspd%keyg_p(2,nlpspd%nseg_p(2*at%nat)+ndebug),stat=i_stat)
  call memocc(i_stat,nlpspd%keyg_p,'nlpspd%keyg_p',subname)
  allocate(nlpspd%keyv_p(nlpspd%nseg_p(2*at%nat)+ndebug),stat=i_stat)
  call memocc(i_stat,nlpspd%keyv_p,'nlpspd%keyv_p',subname)
  allocate(proj(nlpspd%nprojel+ndebug),stat=i_stat)
  call memocc(i_stat,proj,'proj',subname)

  ! After having determined the size of the projector descriptor arrays fill them
  do iat=1,at%nat
     call numb_proj(at%iatype(iat),at%ntypes,at%psppar,at%npspcode,mproj)
     if (mproj.ne.0) then 

        ! coarse grid quantities
        nl1=nlpspd%nboxp_c(1,1,iat) 
        nl2=nlpspd%nboxp_c(1,2,iat) 
        nl3=nlpspd%nboxp_c(1,3,iat) 

        nu1=nlpspd%nboxp_c(2,1,iat)
        nu2=nlpspd%nboxp_c(2,2,iat)
        nu3=nlpspd%nboxp_c(2,3,iat)
        call fill_logrid(at%geocode,n1,n2,n3,nl1,nu1,nl2,nu2,nl3,nu3,0,1,  &
             at%ntypes,at%iatype(iat),rxyz(1,iat),radii_cf(1,3),cpmult,hx,hy,hz,logrid)

        iseg=nlpspd%nseg_p(2*iat-2)+1
        mseg=nlpspd%nseg_p(2*iat-1)-nlpspd%nseg_p(2*iat-2)

        call segkeys(n1,n2,n3,nl1,nu1,nl2,nu2,nl3,nu3,  & 
             logrid,mseg,nlpspd%keyg_p(1,iseg),nlpspd%keyv_p(iseg))

        ! fine grid quantities
        nl1=nlpspd%nboxp_f(1,1,iat)
        nl2=nlpspd%nboxp_f(1,2,iat)
        nl3=nlpspd%nboxp_f(1,3,iat)

        nu1=nlpspd%nboxp_f(2,1,iat)
        nu2=nlpspd%nboxp_f(2,2,iat)
        nu3=nlpspd%nboxp_f(2,3,iat)
        call fill_logrid(at%geocode,n1,n2,n3,nl1,nu1,nl2,nu2,nl3,nu3,0,1,  &
             at%ntypes,at%iatype(iat),rxyz(1,iat),radii_cf(1,2),fpmult,hx,hy,hz,logrid)
        iseg=nlpspd%nseg_p(2*iat-1)+1
        mseg=nlpspd%nseg_p(2*iat)-nlpspd%nseg_p(2*iat-1)
        if (mseg > 0) then
           call segkeys(n1,n2,n3,nl1,nu1,nl2,nu2,nl3,nu3,  & 
                logrid,mseg,nlpspd%keyg_p(1,iseg),nlpspd%keyv_p(iseg))
        end if
     endif
  enddo

  i_all=-product(shape(logrid))*kind(logrid)
  deallocate(logrid,stat=i_stat)
  call memocc(i_stat,i_all,'logrid',subname)

  !fill the projectors if the strategy is a distributed calculation
  if (.not. DistProjApply) then
     !calculate the wavelet expansion of projectors
     call fill_projectors(iproc,n1,n2,n3,hx,hy,hz,at,orbs,rxyz,nlpspd,proj,0)
  end if

END SUBROUTINE createProjectorsArrays
!!***

!!****f* BigDFT/input_wf_diag
!! FUNCTION
!!   input guess wavefunction diagonalization
!! SOURCE
!!
subroutine input_wf_diag(iproc,nproc,at,&
     orbs,nvirt,comms,Glr,hx,hy,hz,rxyz,rhopot,rhocore,pot_ion,&
     nlpspd,proj,pkernel,pkernelseq,ixc,psi,hpsi,psit,G,&
     nscatterarr,ngatherarr,nspin,potshortcut,symObj,irrzon,phnons,GPU,input)
  ! Input wavefunctions are found by a diagonalization in a minimal basis set
  ! Each processors write its initial wavefunctions into the wavefunction file
  ! The files are then read by readwave
  use module_base
  use module_interfaces, except_this_one => input_wf_diag
  use module_types
  use Poisson_Solver
  implicit none
  !Arguments
  integer, intent(in) :: iproc,nproc,ixc,symObj
  integer, intent(inout) :: nspin,nvirt
  real(gp), intent(in) :: hx,hy,hz
  type(atoms_data), intent(in) :: at
  type(orbitals_data), intent(inout) :: orbs
  type(nonlocal_psp_descriptors), intent(in) :: nlpspd
  type(locreg_descriptors), intent(in) :: Glr
  type(communications_arrays), intent(in) :: comms
  type(GPU_pointers), intent(inout) :: GPU
  type(input_variables):: input
  integer, dimension(0:nproc-1,4), intent(in) :: nscatterarr !n3d,n3p,i3s+i3xcsh-1,i3xcsh
  integer, dimension(0:nproc-1,2), intent(in) :: ngatherarr 
  real(gp), dimension(3,at%nat), intent(in) :: rxyz
  real(wp), dimension(nlpspd%nprojel), intent(in) :: proj
  real(dp), dimension(*), intent(inout) :: rhopot,pot_ion
  type(gaussian_basis), intent(out) :: G !basis for davidson IG
  real(wp), dimension(:), pointer :: psi,hpsi,psit,rhocore
  real(dp), dimension(:), pointer :: pkernel,pkernelseq
  integer, intent(in) ::potshortcut
  integer, dimension(*), intent(in) :: irrzon
  real(dp), dimension(*), intent(in) :: phnons
  !local variables
  character(len=*), parameter :: subname='input_wf_diag'
  logical :: switchGPUconv,switchOCLconv
  integer :: i_stat,i_all,iat,nspin_ig,iorb,idum=0
  real(kind=4) :: tt,builtin_rand
  real(gp) :: hxh,hyh,hzh,eks,eexcu,vexcu,epot_sum,ekin_sum,ehart,eexctX,eproj_sum,etol,accurex
  type(orbitals_data) :: orbse
  type(communications_arrays) :: commse
  integer, dimension(:,:), allocatable :: norbsc_arr
  real(wp), dimension(:), allocatable :: potxc
  real(gp), dimension(:), allocatable :: locrad
  type(locreg_descriptors), dimension(:), allocatable :: Llr
  real(wp), dimension(:), pointer :: pot
  real(wp), dimension(:,:,:), pointer :: psigau
type(orbitals_data):: orbsLIN
type(communications_arrays):: commsLIN
real(8),dimension(:),allocatable:: phi, hphi
real(8),dimension(:,:),allocatable:: HamSmall
real(8),dimension(:),allocatable:: eval
integer:: istat
real(8),dimension(:),pointer:: phiWorkPointer

  allocate(norbsc_arr(at%natsc+1,nspin+ndebug),stat=i_stat)
  call memocc(i_stat,norbsc_arr,'norbsc_arr',subname)
  allocate(locrad(at%nat+ndebug),stat=i_stat)
  call memocc(i_stat,locrad,'locrad',subname)

  if (iproc == 0) then
     write(*,'(1x,a)')&
          '------------------------------------------------------- Input Wavefunctions Creation'
  end if

  !spin for inputguess orbitals
  if (nspin == 4) then
     nspin_ig=1
  else
     nspin_ig=nspin
  end if

  call inputguess_gaussian_orbitals(iproc,nproc,at,rxyz,Glr,nvirt,nspin_ig,&
       orbs,orbse,norbsc_arr,locrad,G,psigau,eks)

  !allocate communications arrays for inputguess orbitals
  !call allocate_comms(nproc,orbse,commse,subname)
  call orbitals_communicators(iproc,nproc,Glr,orbse,commse)  

  hxh=.5_gp*hx
  hyh=.5_gp*hy
  hzh=.5_gp*hz

  !check the communication distribution
  !call check_communications(iproc,nproc,orbse,Glr,commse)

  !once the wavefunction coefficients are known perform a set 
  !of nonblocking send-receive operations to calculate overlap matrices

!!!  !create mpirequests array for controlling the success of the send-receive operation
!!!  allocate(mpirequests(nproc-1+ndebug),stat=i_stat)
!!!  call memocc(i_stat,mpirequests,'mpirequests',subname)
!!!
!!!  call nonblocking_transposition(iproc,nproc,G%ncoeff,orbse%isorb+orbse%norbp,&
!!!       orbse%nspinor,psigau,orbse%norb_par,mpirequests)

  !experimental part for building the localisation regions
  if (at%geocode == 'F') then
     !allocate the array of localisation regions
     allocate(Llr(at%nat+ndebug),stat=i_stat)
     !call memocc(i_stat,Llr,'Llr',subname)

     !print *,'locrad',locrad

     call determine_locreg(at%nat,rxyz,locrad,hx,hy,hz,Glr,Llr)

     do iat=1,at%nat
        call deallocate_lr(Llr(iat),subname)
!!$        call deallocate_wfd(Llr(iat)%wfd,subname)
!!$        if (Llr(iat)%geocode=='F') then
!!$           call deallocate_bounds(Llr(iat)%bounds,subname)
!!$        end if
     end do

     !i_all=-product(shape(Llr))*kind(Llr)
     deallocate(Llr,stat=i_stat) !these allocation are special
     !call memocc(i_stat,i_all,'Llr',subname)
  end if

  !allocate the wavefunction in the transposed way to avoid allocations/deallocations
  allocate(psi(orbse%npsidim+ndebug),stat=i_stat)
  call memocc(i_stat,psi,'psi',subname)

  !allocate arrays for the GPU if a card is present
  switchGPUconv=.false.
  switchOCLconv=.false.
  if (GPUconv .and. potshortcut ==0 ) then
     call prepare_gpu_for_locham(Glr%d%n1,Glr%d%n2,Glr%d%n3,nspin_ig,&
          hx,hy,hz,Glr%wfd,orbse,GPU)
  else if (OCLconv .and. potshortcut ==0) then
     call allocate_data_OCL(Glr%d%n1,Glr%d%n2,Glr%d%n3,at%geocode,&
          nspin_ig,hx,hy,hz,Glr%wfd,orbse,GPU)
     if (iproc == 0) write(*,*)&
          'GPU data allocated'
  else if (GPUconv .and. potshortcut >0 ) then
     switchGPUconv=.true.
     GPUconv=.false.
  else if (OCLconv .and. potshortcut >0 ) then
     switchOCLconv=.true.
     OCLconv=.false.
  end if


  !use only the part of the arrays for building the hamiltonian matrix
  call gaussians_to_wavelets_new(iproc,nproc,Glr,orbse,hx,hy,hz,G,&
       psigau(1,1,min(orbse%isorb+1,orbse%norb)),psi)


  i_all=-product(shape(locrad))*kind(locrad)
  deallocate(locrad,stat=i_stat)
  call memocc(i_stat,i_all,'locrad',subname)

  !application of the hamiltonian for gaussian based treatment
  call sumrho(iproc,nproc,orbse,Glr,ixc,hxh,hyh,hzh,psi,rhopot,&
       & Glr%d%n1i*Glr%d%n2i*nscatterarr(iproc,1),nscatterarr,nspin,GPU, &
       & symObj, irrzon, phnons)
     
  !-- if spectra calculation uses a energy dependent potential
  !    input_wf_diag will write (to be used in abscalc)
  !    the density to the file electronic_density.cube
  !  The writing is activated if  5th bit of  in%potshortcut is on.
  if( iand( potshortcut,16)==0 .and. potshortcut /= 0) then
     call plot_density_cube_old(at%geocode,'electronic_density',&
          iproc,nproc,Glr%d%n1,Glr%d%n2,Glr%d%n3,Glr%d%n1i,Glr%d%n2i,Glr%d%n3i,nscatterarr(iproc,2),  & 
          nspin,hxh,hyh,hzh,at,rxyz,ngatherarr,rhopot(1+nscatterarr(iproc,4)*Glr%d%n1i*Glr%d%n2i))
  endif
  !---
  
  if(orbs%nspinor==4) then
     !this wrapper can be inserted inside the poisson solver 
     call PSolverNC(at%geocode,'D',iproc,nproc,Glr%d%n1i,Glr%d%n2i,Glr%d%n3i,&
          nscatterarr(iproc,1),& !this is n3d
          ixc,hxh,hyh,hzh,&
          rhopot,pkernel,pot_ion,ehart,eexcu,vexcu,0.d0,.true.,4)
  else
     !Allocate XC potential
     if (nscatterarr(iproc,2) >0) then
        allocate(potxc(Glr%d%n1i*Glr%d%n2i*nscatterarr(iproc,2)*nspin+ndebug),stat=i_stat)
        call memocc(i_stat,potxc,'potxc',subname)
     else
        allocate(potxc(1+ndebug),stat=i_stat)
        call memocc(i_stat,potxc,'potxc',subname)
     end if

     call XC_potential(at%geocode,'D',iproc,nproc,&
          Glr%d%n1i,Glr%d%n2i,Glr%d%n3i,ixc,hxh,hyh,hzh,&
          rhopot,eexcu,vexcu,nspin,rhocore,potxc)


     if( iand(potshortcut,4)==0) then
        call H_potential(at%geocode,'D',iproc,nproc,&
             Glr%d%n1i,Glr%d%n2i,Glr%d%n3i,hxh,hyh,hzh,&
             rhopot,pkernel,pot_ion,ehart,0.0_dp,.true.)
     endif


     !sum the two potentials in rhopot array
     !fill the other part, for spin, polarised
     if (nspin == 2) then
        call dcopy(Glr%d%n1i*Glr%d%n2i*nscatterarr(iproc,2),rhopot(1),1,&
             rhopot(Glr%d%n1i*Glr%d%n2i*nscatterarr(iproc,2)+1),1)
     end if
     !spin up and down together with the XC part
     call axpy(Glr%d%n1i*Glr%d%n2i*nscatterarr(iproc,2)*nspin,1.0_dp,potxc(1),1,&
          rhopot(1),1)


     i_all=-product(shape(potxc))*kind(potxc)
     deallocate(potxc,stat=i_stat)
     call memocc(i_stat,i_all,'potxc',subname)

  end if

!!!  if (nproc == 1) then
!!!     !calculate the overlap matrix as well as the kinetic overlap
!!!     !in view of complete gaussian calculation
!!!     allocate(ovrlp(G%ncoeff*G%ncoeff),stat=i_stat)
!!!     call memocc(i_stat,ovrlp,'ovrlp',subname)
!!!     allocate(tmp(G%ncoeff,orbse%norb),stat=i_stat)
!!!     call memocc(i_stat,tmp,'tmp',subname)
!!!     allocate(smat(orbse%norb,orbse%norb),stat=i_stat)
!!!     call memocc(i_stat,smat,'smat',subname)
!!!
!!!     !overlap calculation of the gaussian matrix
!!!     call gaussian_overlap(G,G,ovrlp)
!!!     call dsymm('L','U',G%ncoeff,orbse%norb,1.0_gp,ovrlp(1),G%ncoeff,&
!!!          gaucoeff(1,1),G%ncoeff,0.d0,tmp(1,1),G%ncoeff)
!!!
!!!     call gemm('T','N',orbse%norb,orbse%norb,G%ncoeff,1.0_gp,&
!!!          gaucoeff(1,1),G%ncoeff,tmp(1,1),G%ncoeff,0.0_wp,smat(1,1),orbse%norb)
!!!
!!!     !print overlap matrices
!!!     do i=1,orbse%norb
!!!        write(*,'(i5,30(1pe15.8))')i,(smat(i,iorb),iorb=1,orbse%norb)
!!!     end do
!!!
!!!     !overlap calculation of the kinetic operator
!!!     call kinetic_overlap(G,G,ovrlp)
!!!     call dsymm('L','U',G%ncoeff,orbse%norb,1.0_gp,ovrlp(1),G%ncoeff,&
!!!          gaucoeff(1,1),G%ncoeff,0.d0,tmp(1,1),G%ncoeff)
!!!
!!!     call gemm('T','N',orbse%norb,orbse%norb,G%ncoeff,1.0_gp,&
!!!          gaucoeff(1,1),G%ncoeff,tmp(1,1),G%ncoeff,0.0_wp,smat(1,1),orbse%norb)
!!!
!!!     !print overlap matrices
!!!     tt=0.0_wp
!!!     do i=1,orbse%norb
!!!        write(*,'(i5,30(1pe15.8))')i,(smat(i,iorb),iorb=1,orbse%norb)
!!!        !write(12,'(i5,30(1pe15.8))')i,(smat(i,iorb),iorb=1,orbse%norb)
!!!        tt=tt+smat(i,i)
!!!     end do
!!!     print *,'trace',tt
!!!
!!!     !overlap calculation of the kinetic operator
!!!     call cpu_time(t0)
!!!     call potential_overlap(G,G,rhopot,Glr%d%n1i,Glr%d%n2i,Glr%d%n3i,hxh,hyh,hzh,&
!!!          ovrlp)
!!!     call cpu_time(t1)
!!!     call dsymm('L','U',G%ncoeff,orbse%norb,1.0_gp,ovrlp(1),G%ncoeff,&
!!!          gaucoeff(1,1),G%ncoeff,0.d0,tmp(1,1),G%ncoeff)
!!!
!!!     call gemm('T','N',orbse%norb,orbse%norb,G%ncoeff,1.0_gp,&
!!!          gaucoeff(1,1),G%ncoeff,tmp(1,1),G%ncoeff,0.0_wp,smat(1,1),orbse%norb)
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
!!!     i_all=-product(shape(ovrlp))*kind(ovrlp)
!!!     deallocate(ovrlp,stat=i_stat)
!!!     call memocc(i_stat,i_all,'ovrlp',subname)
!!!     i_all=-product(shape(tmp))*kind(tmp)
!!!     deallocate(tmp,stat=i_stat)
!!!     call memocc(i_stat,i_all,'tmp',subname)
!!!     i_all=-product(shape(smat))*kind(smat)
!!!     deallocate(smat,stat=i_stat)
!!!     call memocc(i_stat,i_all,'smat',subname)
!!!  end if

  if(potshortcut>0) then
!!$    if (GPUconv) then
!!$       call free_gpu(GPU,orbs%norbp)
!!$    end if
     if (switchGPUconv) then
        GPUconv=.true.
     end if
     if (switchOCLconv) then
        OCLconv=.true.
     end if

     call deallocate_orbs(orbse,subname)
     
     !deallocate the gaussian basis descriptors
     call deallocate_gwf(G,subname)
    
     i_all=-product(shape(psigau))*kind(psigau)
     deallocate(psigau,stat=i_stat)
     call memocc(i_stat,i_all,'psigau',subname)
     call deallocate_comms(commse,subname)
     i_all=-product(shape(norbsc_arr))*kind(norbsc_arr)
     deallocate(norbsc_arr,stat=i_stat)
     call memocc(i_stat,i_all,'norbsc_arr',subname)
    return 
  end if

  !allocate the wavefunction in the transposed way to avoid allocations/deallocations
  allocate(hpsi(orbse%npsidim+ndebug),stat=i_stat)
  call memocc(i_stat,hpsi,'hpsi',subname)

  !call dcopy(orbse%npsidim,psi,1,hpsi,1)
  if (input%exctxpar == 'OP2P') eexctX = -99.0_gp

  call full_local_potential(iproc,nproc,Glr%d%n1i*Glr%d%n2i*nscatterarr(iproc,2),Glr%d%n1i*Glr%d%n2i*Glr%d%n3i,nspin,&
       orbse%norb,orbse%norbp,ngatherarr,rhopot,pot)



! THIS IS NOW DONE IN CLUSTER
!call initializeParameters(iproc, nproc, Glr, orbs, orbsLIN, commsLIN, at, phi)
!call improveOrbitals(iproc, nproc, nspin, Glr, orbs, orbsLIN, commsLIN, at, rxyz, nscatterarr, ngatherarr, nlpspd, proj, rhopot, GPU, input, pkernelseq, phi)


!!!write(*,*) 'calling getLocalizedBasis'
!!!call initializeParameters(iproc, nproc, Glr, orbs, orbsLIN, commsLIN, at)
!!!allocate(phi(orbsLIN%npsidim), stat=istat)
!!!call initRandomSeed(iproc, 1)
!!!call random_number(phi)
!!!allocate(hphi(orbsLIN%npsidim), stat=istat)
!!!! Initialize phi at random
!!!call getLocalizedBasis(iproc, nproc, at, orbs, Glr, input, orbsLIN, commsLIN, rxyz, nspin, nlpspd, proj, nscatterarr, ngatherarr, rhopot, GPU, pkernelseq, phi, hphi)
!!!allocate(phiWorkPointer(size(phi)), stat=istat)
!!!call transpose_v(iproc, nproc, orbsLIN, Glr%wfd, commsLIN, phi, work=phiWorkPointer)
!!!call transpose_v(iproc, nproc, orbsLIN, Glr%wfd, commsLIN, hphi, work=phiWorkPointer)
!!!allocate(HamSmall(orbsLIN%norb,orbsLIN%norb), stat=istat)
!!!call transformHam(iproc, nproc, orbsLIN, commsLIN, phi, hphi, HamSmall)
!!!allocate(eval(orbsLIN%norb), stat=istat)
!!!call diagonalizeHamiltonian(iproc, nproc, orbsLIN, HamSmall, eval)
!!!
!!!! Store the new wave function in hphi as a temporary array
!!!call buildWavefunction(iproc, nproc, orbsLIN, commsLIN, phi, hphi, HamSmall)
!!!call dcopy(orbsLIN%npsidim, hphi(1), 1, phi(1), 1)
!!!call untranspose_v(iproc, nproc, orbsLIN, Glr%wfd, commsLIN, phi, work=phiWorkPointer)
!!!! Not necessary to untranspose hphi (no longer needed)
!!!!call untranspose_v(iproc, nproc, orbsLIN, Glr%wfd, commsLIN, hphi, work=phiWorkPointer)
!!!deallocate(phiWorkPointer, stat=istat)

  call HamiltonianApplication(iproc,nproc,at,orbse,hx,hy,hz,rxyz,&
       nlpspd,proj,Glr,ngatherarr,pot,&
       psi,hpsi,ekin_sum,epot_sum,eexctX,eproj_sum,nspin,GPU,pkernel=pkernelseq)

  !deallocate potential
  call free_full_potential(nproc,pot,subname)

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
!!!  call wavelets_to_gaussians(at%geocode,orbse%norbp,Glr%d%n1,Glr%d%n2,Glr%d%n3,G,&
!!!       thetaphi,hx,hy,hz,Glr%wfd,hpsi,hpsigau)
!!!
!!!  i_all=-product(shape(thetaphi))*kind(thetaphi)
!!!  deallocate(thetaphi,stat=i_stat)
!!!  call memocc(i_stat,i_all,'thetaphi',subname)

  accurex=abs(eks-ekin_sum)
  !tolerance for comparing the eigenvalues in the case of degeneracies
  etol=accurex/real(orbse%norbu,gp)
  if (iproc == 0 .and. verbose > 1) write(*,'(1x,a,2(f19.10))') 'done. ekin_sum,eks:',ekin_sum,eks
  if (iproc == 0) then
     write(*,'(1x,a,3(1x,1pe18.11))') 'ekin_sum,epot_sum,eproj_sum',  & 
          ekin_sum,epot_sum,eproj_sum
     write(*,'(1x,a,3(1x,1pe18.11))') '   ehart,   eexcu,    vexcu',ehart,eexcu,vexcu
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
     call free_gpu(GPU,orbse%norbp)
  else if (OCLconv) then
     call free_gpu_OCL(GPU,orbse,nspin_ig)
  end if

  if (iproc == 0 .and. verbose > 1) write(*,'(1x,a)',advance='no')&
       'Input Wavefunctions Orthogonalization:'

  !psivirt can be eliminated here, since it will be allocated before davidson
  !with a gaussian basis
!!$  call DiagHam(iproc,nproc,at%natsc,nspin_ig,orbs,Glr%wfd,comms,&
!!$       psi,hpsi,psit,orbse,commse,etol,norbsc_arr,orbsv,psivirt)

  call DiagHam(iproc,nproc,at%natsc,nspin_ig,orbs,Glr%wfd,comms,&
       psi,hpsi,psit,input,orbse,commse,etol,norbsc_arr)

  if (input%itrpmax > 1) then
     !use the eval array of orbse structure to save the original values
     allocate(orbse%eval(orbs%norb*orbs%nkpts+ndebug),stat=i_stat)
     call memocc(i_stat,orbse%eval,'orbse%eval',subname)
     
     call dcopy(orbs%norb*orbs%nkpts,orbs%eval(1),1,orbse%eval(1),1)

     !add a small displacement in the eigenvalues
     do iorb=1,orbs%norb*orbs%nkpts
        if (iorb <= orbs%norb) then
           if (orbs%efermi == -.1_gp .and. orbs%occup(iorb) < 1.0_gp) then
              orbs%efermi=orbs%eval(iorb)
           end if
        end if
        tt=builtin_rand(idum)
        orbs%eval(iorb)=orbs%eval(iorb)*(1.0_gp+input%Tel*real(tt,gp))
        !use the first k-point to guess fermi energy input
     end do

     !correct the occupation numbers wrt fermi level
     call evaltoocc(iproc,.false.,input%Tel,orbs)

     !restore the occupation numbers
     call dcopy(orbs%norb*orbs%nkpts,orbse%eval(1),1,orbs%eval(1),1)

     i_all=-product(shape(orbse%eval))*kind(orbse%eval)
     deallocate(orbse%eval,stat=i_stat)
     call memocc(i_stat,i_all,'orbse%eval',subname)
  end if

  call deallocate_comms(commse,subname)

  i_all=-product(shape(norbsc_arr))*kind(norbsc_arr)
  deallocate(norbsc_arr,stat=i_stat)
  call memocc(i_stat,i_all,'norbsc_arr',subname)

  if (iproc == 0) then
     if (verbose > 1) write(*,'(1x,a)')'done.'
     !gaussian estimation valid only for Free BC
     if (at%geocode == 'F') then
        write(*,'(1x,a,1pe9.2)') 'expected accuracy in energy ',accurex
        write(*,'(1x,a,1pe9.2)') &
          'expected accuracy in energy per orbital ',accurex/real(orbs%norb,kind=8)
        !write(*,'(1x,a,1pe9.2)') &
        !     'suggested value for gnrm_cv ',accurex/real(orbs%norb,kind=8)
     end if
  endif

  !here we can define the subroutine which generates the coefficients for the virtual orbitals
  call deallocate_gwf(G,subname)

  i_all=-product(shape(psigau))*kind(psigau)
  deallocate(psigau,stat=i_stat)
  call memocc(i_stat,i_all,'psigau',subname)

  call deallocate_orbs(orbse,subname)



END SUBROUTINE input_wf_diag
!!***

!contains





!!subroutine getLocalizedBasis(iproc, nproc, at, orbs, Glr, input, orbsLIN, commsLIN, rxyz, nspin, nlpspd, &
!!    proj, nscatterarr, ngatherarr, rhopot, GPU, pkernelseq, phi, hphi, trH, rxyzParabola, &
!!    idsxMin, idsxMax, infoBasisFunctions)
!!!
!!! Purpose:
!!! ========
!!!   Calculates the localized basis functions phi. These basis functions are eigenfunctions of the ordinary Hamiltonian
!!!   with an additional parabolic potential centered at the atoms. The eigenfunctions are determined by minimizing the trace.
!!!
!!! Calling arguments:
!!!   Input arguments
!!!   Output arguments
!!!    phi   the localized basis functions
!!!
!!use module_base
!!use module_types
!!use module_interfaces, except_this_one => getLocalizedBasis
!!  use Poisson_Solver
!!!use allocModule
!!implicit none
!!
!!! Calling arguments
!!integer:: iproc, nproc, idsxMin, idsxMax, infoBasisFunctions
!!type(atoms_data), intent(in) :: at
!!type(orbitals_data):: orbs
!!type(locreg_descriptors), intent(in) :: Glr
!!type(input_variables):: input
!!type(orbitals_data):: orbsLIN
!!type(communications_arrays):: commsLIN
!!real(8),dimension(3,at%nat):: rxyz, rxyzParabola
!!integer:: nspin
!!type(nonlocal_psp_descriptors), intent(in) :: nlpspd
!!real(wp), dimension(nlpspd%nprojel), intent(in) :: proj
!!integer, dimension(0:nproc-1,4), intent(in) :: nscatterarr !n3d,n3p,i3s+i3xcsh-1,i3xcsh
!!integer, dimension(0:nproc-1,2), intent(in) :: ngatherarr 
!!real(dp), dimension(*), intent(inout) :: rhopot
!!type(GPU_pointers), intent(inout) :: GPU
!!real(dp), dimension(:), pointer :: pkernelseq
!!real(8),dimension(orbsLIN%npsidim):: phi, hphi
!!real(8):: trH
!!
!!! Local variables
!!real(8) ::epot_sum,ekin_sum,eexctX,eproj_sum
!!integer:: iorb, jorb, j, icountSDSatur, iat, icountSwitch, idsx, icountDIISFailureTot, icountDIISFailureCons, itBest
!!integer:: istat, istart, jstart, ierr, i, i0, i1, i2, i3, j0, j1, ii, jj, ix0, iy0, iz0, iseg, jproc, it
!!integer:: itMax, nbasisPerAtForDebug, icount, ishift, ncong, ix, iy, iz, iiAt, info, lwork, iall, nvctrp
!!real(8),dimension(:),allocatable:: hphiold, alpha
!!real(8),dimension(:),allocatable:: phir, eval, work
!!real(8),dimension(:,:),allocatable:: HamSmall
!!real(8):: hxh, hyh, hzh, dis, tt, rr, trace, ddot, dnrm2, fnrm, fnrmMax, meanAlpha, gnrm, gnrm_zero
!!real(8):: kx, ky, kz, tt1, tt2, tt3, parabShift, gnrmMax, trHm1, trHm2, d2trH
!!type(workarr_sumrho) :: w
!!type(workarr_locham) :: w2
!!!integer,dimension(:),allocatable:: onWhichAtom
!!real(8),dimension(:),pointer:: phiWorkPointer
!!real(8),dimension(:),allocatable:: phiWorkPointer2
!!logical:: debug, precond, quiet, allowDIIS
!!logical,dimension(:),allocatable:: move
!!character(len=*),parameter:: subname='getLocalizedBasis'
!!real(8),dimension(:),allocatable:: diag, fnrmOldArr
!!real(8),dimension(:,:),allocatable:: fnrmArr, fnrmOvrlpArr
!!type(diis_objects):: diisLIN
!!type(diis_objects),dimension(:),allocatable:: diisArr
!!
!!
!!
!!call allocateLocalArrays()
!!
!!
!!
!!icountSDSatur=0
!!icountSwitch=0
!!icountDIISFailureTot=0
!!icountDIISFailureCons=0
!!
!!
!!! No DIIS in the beginning
!!call initializeDIISParameters(idsxMax)
!!allowDIIS=.true.
!!if(allowDIIS) then
!!else
!!    diisLIN%idsx=0
!!    call deallocate_diis_objects(diisLIN,subname)
!!end if
!!
!!
!!allocate(diisArr(orbsLIN%norb), stat=istat)
!!do iorb=1,orbsLIN%norb
!!    ! parameters for DIIS
!!    diisArr(iorb)%switchSD=.false.
!!    diisArr(iorb)%idiistol=0
!!    diisArr(iorb)%mids=1
!!    diisArr(iorb)%ids=0
!!    diisArr(iorb)%idsx=idsxMax
!!    diisArr(iorb)%energy_min=1.d10
!!    diisArr(iorb)%energy_old=1.d10
!!    diisArr(iorb)%energy=1.d10
!!    diisArr(iorb)%alpha=2.d0
!!    ! allocate only matrix, no history array (therefore '0')
!!    call allocate_diis_objects(diisArr(iorb)%idsx, 0, 1, diisArr(iorb), subname) ! 1 for k-points
!!end do
!!
!!
!!
!!
!!if(iproc==0) write(*,'(a)') '============================ basis functions creation... ============================'
!!itMax=10000
!!!alpha=1.d-3
!!alpha=1.d-2
!!precond=.true.
!!trHm1=0.d0
!!trHm2=0.d0
!!if(iproc==0) write(*,'(a,i0)') 'using DIIS with history length ', diisLIN%idsx
!!if(iproc==0) write(*,'(a,es12.5)') 'convergence criterion is', orbsLIN%convCrit
!!if(iproc==0) write(*,'(a,i0)') 'minimal number of iterations: ', orbsLIN%nItMin
!!if(iproc==0) write(*,'(a,i0)') 'maximal number of iterations: ', orbsLIN%nItMax
!!iterLoop: do it=1,itMax
!!    trace=0.d0
!!    fnrmMax=0.d0
!!    fnrm=0.d0
!!
!!    if (iproc==0) then
!!        write( *,'(1x,a,i0)') repeat('-',77 - int(log(real(it))/log(10.))) // ' iter=', it
!!    endif
!!
!!
!!    ! Orthonormalize the orbitals.
!!    if(iproc==0) then
!!        write(*,'(x,a)', advance='no') 'Orthonormalization... '
!!    end if
!!    call orthogonalize(iproc, nproc, orbsLIN, commsLIN, Glr%wfd, phi, input)
!!
!!    ! Untranspose phi
!!    call untranspose_v(iproc, nproc, orbsLIN, Glr%wfd, commsLIN, phi, work=phiWorkPointer)
!!
!!
!!    ! Calculate the unconstrained gradient.
!!    if(iproc==0) then
!!        write(*,'(a)', advance='no') 'Hamiltonian application... '
!!    end if
!!    call HamiltonianApplicationParabola(iproc,nproc,at,orbsLIN,input%hx,input%hy,input%hz,rxyz,&
!!         nlpspd,proj,Glr,ngatherarr,Glr%d%n1i*Glr%d%n2i*nscatterarr(iproc,2),&
!!         rhopot(1),&
!!         phi(1),hphi(1),ekin_sum,epot_sum,eexctX,eproj_sum,nspin,GPU, orbsLIN%onWhichAtom, rxyzParabola, pkernel=pkernelseq)
!!
!!
!!    ! Apply the orthoconstraint to the gradient. This subroutine also calculates the trace trH.
!!    if(iproc==0) then
!!        write(*,'(a)', advance='no') 'orthoconstraint... '
!!    end if
!!    call transpose_v(iproc, nproc, orbsLIN, Glr%wfd, commsLIN, hphi, work=phiWorkPointer)
!!    call transpose_v(iproc, nproc, orbsLIN, Glr%wfd, commsLIN, phi, work=phiWorkPointer)
!!    call orthoconstraintNotSymmetric(iproc, nproc, orbsLIN, commsLIN, Glr%wfd, phi, hphi, trH, diag)
!!
!!
!!    ! Calculate the norm of the gradient (fnrmArr) and determine the angle between the current gradient and that
!!    ! of the previous iteration (fnrmOvrlpArr).
!!    nvctrp=commsLIN%nvctr_par(iproc,1) ! 1 for k-point
!!    istart=1
!!    do iorb=1,orbsLIN%norb
!!        if(it>1) fnrmOvrlpArr(iorb,2)=ddot(nvctrp*orbs%nspinor, hphi(istart), 1, hphiold(istart), 1)
!!        fnrmArr(iorb,2)=ddot(nvctrp*orbs%nspinor, hphi(istart), 1, hphi(istart), 1)
!!        istart=istart+nvctrp*orbs%nspinor
!!    end do
!!    call mpi_allreduce(fnrmArr(1,2), fnrmArr(1,1), orbsLIN%norb, mpi_double_precision, mpi_sum, mpi_comm_world, ierr)
!!    call mpi_allreduce(fnrmOvrlpArr(1,2), fnrmOvrlpArr(1,1), orbsLIN%norb, mpi_double_precision, mpi_sum, mpi_comm_world, ierr)
!!
!!    ! Keep the gradient for the next iteration.
!!    if(it>1) then
!!        call dcopy(orbsLIN%norb, fnrmArr(1,1), 1, fnrmOldArr(1), 1)
!!    end if
!!
!!    ! Determine the gradient norm and its maximal component. In addition, adapt the
!!    !  step size for the steepest descent minimization (depending on the angle 
!!    ! between the current gradient and the one from the previous iteration).
!!    ! This is of course only necessary if we are using steepest descent and not DIIS.
!!    do iorb=1,orbsLIN%norb
!!        fnrm=fnrm+fnrmArr(iorb,1)
!!        if(fnrmArr(iorb,1)>fnrmMax) fnrmMax=fnrmArr(iorb,1)
!!        if(it>1 .and. diisLIN%idsx==0 .and. .not.diisLIN%switchSD) then
!!        ! Adapt step size for the steepest descent minimization.
!!            tt=fnrmOvrlpArr(iorb,1)/sqrt(fnrmArr(iorb,1)*fnrmOldArr(iorb))
!!            if(tt>.7d0) then
!!                alpha(iorb)=alpha(iorb)*1.05d0
!!            else
!!                alpha(iorb)=alpha(iorb)*.5d0
!!            end if
!!        end if
!!    end do
!!    fnrm=sqrt(fnrm)
!!    fnrmMax=sqrt(fnrmMax)
!!    ! Copy the gradient (will be used in the next iteration to adapt the step size).
!!    call dcopy(orbsLIN%norb*nvctrp*orbs%nspinor, hphi(1), 1, hphiold(1), 1)
!!
!!    ! Untranspose hphi.
!!    call untranspose_v(iproc, nproc, orbsLIN, Glr%wfd, commsLIN, hphi, work=phiWorkPointer)
!!    !deallocate(phiWorkPointer, stat=istat)
!!
!!
!!    ! Precondition the gradient
!!    if(iproc==0) then
!!        write(*,'(a)') 'preconditioning. '
!!    end if
!!    ncong=10
!!    gnrm=1.d3 ; gnrm_zero=1.d3
!!    if(precond) call preconditionallLIN(iproc, nproc, orbsLIN, Glr, input%hx, input%hy, input%hz, &
!!        ncong, hphi, gnrm, gnrm_zero, gnrmMax,  at%nat, rxyz, orbsLIN%onWhichAtom, at, it)
!!    !if(iproc==0) then
!!    !    write(*,'(a)') 'done. '
!!    !end if
!!
!!    tt=gnrm
!!    call mpi_allreduce(tt, gnrm, 1, mpi_double_precision, mpi_sum, mpi_comm_world, ierr)
!!    gnrm=sqrt(gnrm)
!!    tt=gnrmMax
!!    call mpi_allreduce(tt, gnrmMax, 1, mpi_double_precision, mpi_max, mpi_comm_world, ierr)
!!    gnrmMax=sqrt(gnrmMax)
!!
!!
!!    tt=sum(alpha)
!!    meanAlpha=tt/dble(orbsLIN%norb)
!!
!!    ! Keep the values of the two previous iterations
!!    trHm2=trHm1
!!    trHm1=trH 
!!    d2trH=trHm2-2.d0*trHm1+trH
!!
!!
!!
!!    if(iproc==0) write(*,'(x,a,i6,2es15.7,f14.7)') 'iter, fnrm, fnrmMax, trace', it, fnrm, fnrmMax, trH
!!    if(iproc==0) write(1000,'(i6,2es15.7,f15.7,es12.4)') it, fnrm, fnrmMax, trH, meanAlpha
!!    !if(fnrmMax<1.d0) allowDIIS=.true.
!!    !if(fnrmMax<1.d-2) then
!!    !if(fnrmMax<1.d-2 .and. it>=15) then
!!    if((fnrmMax<orbsLIN%convCrit .and. it>=orbsLIN%nItMin) .or. it>=orbsLIN%nItMax) then
!!        if(it>=orbsLIN%nItMax) then
!!            if(iproc==0) write(*,'(a,i0,a)') 'WARNING: not converged within ', it, &
!!                ' iterations! Exiting loop due to limitations of iterations.'
!!            if(iproc==0) write(*,'(a,2es15.7,f12.7)') 'Final values for fnrm, fnrmMax, trace: ', fnrm, fnrmMax, trH
!!            infoBasisFunctions=1
!!        else
!!            if(iproc==0) then
!!                write(*,'(a,i0,a,2es15.7,f12.7)') 'converged in ', it, ' iterations.'
!!                write (*,'(a,2es15.7,f12.7)') 'Final values for fnrm, fnrmMax, trace: ', fnrm, fnrmMax, trH
!!            end if
!!            infoBasisFunctions=0
!!        end if
!!        if(iproc==0) write(*,'(a)') '============================== basis functions created =============================='
!!        !allocate(phiWorkPointer(size(phi)), stat=istat)
!!        call untranspose_v(iproc, nproc, orbsLIN, Glr%wfd, commsLIN, phi, work=phiWorkPointer)
!!        !deallocate(phiWorkPointer, stat=istat)
!!        call plotOrbitals(iproc, orbsLIN, Glr, phi, at%nat, rxyz, orbsLIN%onWhichAtom, .5d0*input%hx, &
!!            .5d0*input%hy, .5d0*input%hz, 1)
!!        exit iterLoop
!!    end if
!!    if(fnrmMax<1.d2 .and. .not.precond) then
!!        if(iproc==0) write(*,'(a)') 'starting preconditioner...'
!!        alpha=10.d0*alpha
!!        precond=.true.
!!    end if
!!
!!
!!    call DIISorSD()
!!    if(iproc==0) then
!!        if(diisLIN%idsx>0) then
!!            write(*,'(x,3(a,i0))') 'DIIS informations: history length=',diisLIN%idsx, ', consecutive failures=', &
!!                icountDIISFailureCons, ', total failures=', icountDIISFailureTot
!!        else
!!            write(*,'(x,a,es9.3,a,i0)') 'steepest descent informations: mean alpha=', meanAlpha, &
!!            ', consecutive successes=', icountSDSatur
!!        end if
!!    end if
!!    if(.not. diisLIN%switchSD) call improve()
!!
!!
!!
!! 
!!
!!end do iterLoop
!!
!!
!!
!!
!!
!!
!!
!!
!!contains
!!
!!    subroutine initializeDIISParameters(idsxHere)
!!    ! Purpose:
!!    ! ========
!!    !   Initializes all parameters needed for the DIIS procedure.
!!    !
!!    ! Calling arguments
!!    !   idsx    DIIS history length
!!    !
!!    implicit none
!!    
!!    ! Calling arguments
!!    integer:: idsxHere
!!
!!    diisLIN%switchSD=.false.
!!    diisLIN%idiistol=0
!!    diisLIN%mids=1
!!    diisLIN%ids=0
!!    diisLIN%idsx=idsxHere
!!    diisLIN%energy_min=1.d10
!!    diisLIN%energy_old=1.d10
!!    diisLIN%energy=1.d10
!!    diisLIN%alpha=2.d0
!!    call allocate_diis_objects(diisLIN%idsx, orbsLIN%npsidim, 1, diisLIN, subname) ! 1 for k-points
!!    end subroutine initializeDIISParameters
!!
!!
!!    subroutine DIISorSD()
!!    if(diisLIN%switchSD) diisLIN%switchSD=.false.
!!    ! Determine wheter to use DIIS or SD
!!    if(trH<=diisLIN%energy_min) then
!!        ! everything ok
!!        diisLIN%energy_min=trH
!!        diisLIN%switchSD=.false.
!!        itBest=it
!!        icountSDSatur=icountSDSatur+1
!!        icountDIISFailureCons=0
!!        if(icountSDSatur>=10 .and. diisLIN%idsx==0 .and. allowDIIS) then
!!            ! switch back to DIIS 
!!            icountSwitch=icountSwitch+1
!!            !diisLIN%idsx=idsx
!!            idsx=max(idsxMin,idsxMax-icountSwitch)
!!            !diisLIN%idsx=idsxMax
!!            if(iproc==0) write(*,'(a,i0)') 'switch to DIIS with new history length ', idsx
!!            call initializeDIISParameters(idsx)
!!            !diisLIN%ids=0
!!            !diisLIN%mids=1
!!            !call allocate_diis_objects(diisLIN%idsx, orbsLIN%npsidim, 1, diisLIN, subname) ! 1 for k-points
!!            icountDIISFailureTot=0
!!            icountDIISFailureCons=0
!!        end if
!!    else
!!        ! the trace is growing.
!!        ! Count how many times this occurs and switch to SD after 3 times.
!!        icountDIISFailureCons=icountDIISFailureCons+1
!!        icountDIISFailureTot=icountDIISFailureTot+1
!!        icountSDSatur=0
!!        if((icountDIISFailureCons>=2 .or. icountDIISFailureTot>=3) .and. diisLIN%idsx>0) then
!!            alpha=1.d0
!!            if(iproc==0) then
!!                if(icountDIISFailureCons>=2) write(*,'(x,a,i0,a,es10.3)') 'DIIS failed ', &
!!                    icountDIISFailureCons, ' times consecutievly. Switch to SD with stepsize', alpha(1)
!!                if(icountDIISFailureTot>=3) write(*,'(x,a,i0,a,es10.3)') 'DIIS failed ', &
!!                    icountDIISFailureTot, ' times in total. Switch to SD with stepsize', alpha(1)
!!            end if
!!            ! Try to get back the orbitals of the best iteration. This is possible if
!!            ! these orbitals are still present in the DIIS history.
!!            if(it-itBest<diisLIN%idsx) then
!!               if(iproc==0) then
!!                   if(iproc==0) write(*,'(x,a,i0,a)')  'Recover the orbitals from iteration ', &
!!                       itBest, ' which are the best so far.'
!!               end if
!!               ii=modulo(diisLIN%mids-(it-itBest),diisLIN%mids)
!!               !write(*,'(a,2i5)') 'diisLIN%mids, ii', diisLIN%mids, ii
!!               nvctrp=commsLIN%nvctr_par(iproc,1) ! 1 for k-point
!!               call dcopy(orbsLIN%norb*nvctrp, diisLIN%psidst(ii*nvctrp*orbsLIN%norb+1), 1, phi(1), 1)
!!            end if
!!            call deallocate_diis_objects(diisLIN, subname)
!!            diisLIN%idsx=0
!!            diisLIN%switchSD=.true.
!!        end if
!!    end if
!!    end subroutine DIISorSD
!!
!!
!!    subroutine improve()
!!    ! For DIIS 
!!    if (diisLIN%idsx > 0) then
!!       diisLIN%mids=mod(diisLIN%ids,diisLIN%idsx)+1
!!       diisLIN%ids=diisLIN%ids+1
!!       do iorb=1,orbsLIN%norb
!!           diisArr(iorb)%mids=mod(diisArr(iorb)%ids,diisArr(iorb)%idsx)+1
!!           diisArr(iorb)%ids=diisArr(iorb)%ids+1
!!       end do
!!    end if
!!
!!    ! Follow the gradient using steepest descent.
!!    ! The same, but transposed
!!    !allocate(phiWorkPointer(size(phi)), stat=istat)
!!    call transpose_v(iproc, nproc, orbsLIN, Glr%wfd, commsLIN, hphi, work=phiWorkPointer)
!!    
!!    ! steepest descent
!!    if(diisLIN%idsx==0) then
!!        istart=1
!!        nvctrp=commsLIN%nvctr_par(iproc,1) ! 1 for k-point
!!        do iorb=1,orbsLIN%norb
!!            call daxpy(nvctrp*orbs%nspinor, -alpha(iorb), hphi(istart), 1, phi(istart), 1)
!!            istart=istart+nvctrp*orbs%nspinor
!!        end do
!!    else
!!        ! DIIS
!!        quiet=.true. ! less output
!!        call psimix(iproc, nproc, orbsLIN, commsLIN, diisLIN, hphi, phi, quiet)
!!        !call psimixVariable(iproc, nproc, orbsLIN, commsLIN, diisLIN, diisArr, hphi, phi, quiet)
!!    end if
!!    !deallocate(phiWorkPointer, stat=istat)
!!    end subroutine improve
!!
!!
!!
!!    subroutine allocateLocalArrays()
!!    !
!!    ! Purpose:
!!    ! ========
!!    !   This subroutine allocates all local arrays.
!!    !
!!
!!    allocate(hphiold(orbsLIN%npsidim), stat=istat)
!!    call memocc(istat, orbsLIN%npsidim, 'hphiold', subname)
!!
!!    allocate(alpha(orbsLIN%norb), stat=istat)
!!    call memocc(istat, orbsLIN%norbp, 'alpha', subname)
!!
!!    allocate(fnrmArr(orbsLIN%norb,2), stat=istat)
!!    call memocc(istat, orbsLIN%norb*2, 'fnrmArr', subname)
!!
!!    allocate(fnrmOldArr(orbsLIN%norb), stat=istat)
!!    call memocc(istat, orbsLIN%norb, 'fnrmOldArr', subname)
!!
!!    allocate(fnrmOvrlpArr(orbsLIN%norb,2), stat=istat)
!!    call memocc(istat, orbsLIN%norb*2, 'fnrmOvrlpArr', subname)
!!
!!    allocate(phiWorkPointer(size(phi)), stat=istat)
!!    call memocc(istat, size(phi), 'phiWorkPointer', subname)
!!    
!!    allocate(diag(orbsLIN%norb), stat=istat)
!!    
!!    ! Allocate the work arrays which will be used for the preconditioning.
!!    call initialize_work_arrays_locham(Glr,orbs%nspinor,w2)
!!
!!    end subroutine allocateLocalArrays
!!
!!
!!    subroutine deallocateLocalArrays()
!!    !
!!    ! Purpose:
!!    ! ========
!!    !   This subroutine deallocate all local arrays.
!!    !
!!
!!    call deallocate_work_arrays_locham(Glr,w2)
!!    
!!    iall=-product(shape(hphiold))*kind(hphiold)
!!    deallocate(hphiold, stat=istat)
!!    call memocc(istat, iall, 'hphiold', subname)
!!    
!!    iall=-product(shape(alpha))*kind(alpha)
!!    deallocate(alpha, stat=istat)
!!    call memocc(istat, iall, 'alpha', subname)
!!
!!    iall=-product(shape(phiWorkPointer))*kind(phiWorkPointer)
!!    deallocate(phiWorkPointer, stat=istat)
!!    call memocc(istat, iall, 'phiWorkPointer', subname)
!!    
!!    ! if diisLIN%idsx==0, these arrays have already been deallocated
!!    if(diisLIN%idsx>0 .and. idsxMax>0) call deallocate_diis_objects(diisLIN,subname)
!!
!!    end subroutine deallocateLocalArrays
!!
!!
!!end subroutine getLocalizedBasis





!END SUBROUTINE input_wf_diag






!subroutine plotOrbitals(iproc, orbs, Glr, phi, nat, rxyz, onWhichAtom, hxh, hyh, hzh, it)
subroutine plotOrbitals(iproc, orbs, Glr, phi, nat, rxyz, onWhichAtom, hxh, hyh, hzh, it)
!
! Plots the orbitals
!
use module_base
use module_types
implicit none

! Calling arguments
integer:: iproc
type(orbitals_data), intent(inout) :: orbs
type(locreg_descriptors), intent(in) :: Glr
real(8),dimension((Glr%wfd%nvctr_c+7*Glr%wfd%nvctr_f)*orbs%nspinor*orbs%norbp):: phi
integer:: nat
real(8),dimension(3,nat):: rxyz
integer,dimension(orbs%norbp):: onWhichAtom
real(8):: hxh, hyh, hzh
integer:: it

integer:: ix, iy, iz, ix0, iy0, iz0, iiAt, jj, iorb, i1, i2, i3, istart, ii, istat
real(8),dimension(:),allocatable:: phir
type(workarr_sumrho) :: w

allocate(phir(Glr%d%n1i*Glr%d%n2i*Glr%d%n3i), stat=istat)

call initialize_work_arrays_sumrho(Glr,w)

istart=0
   
!write(*,*) 'write, orbs%nbasisp', orbs%norbp
    orbLoop: do iorb=1,orbs%norbp
        call daub_to_isf(Glr,w,phi(istart+1),phir(1))
        iiAt=onWhichAtom(iorb)
        ix0=nint(rxyz(1,iiAt)/hxh)
        iy0=nint(rxyz(2,iiAt)/hyh)
        iz0=nint(rxyz(3,iiAt)/hzh)

        jj=0
        open(unit=(iproc+1)*1000000+it*1000+iorb*10+7)
        open(unit=(iproc+1)*1000000+it*1000+iorb*10+8)
        open(unit=(iproc+1)*1000000+it*1000+iorb*10+9)
        do i3=1,Glr%d%n3i
            do i2=1,Glr%d%n2i
                do i1=1,Glr%d%n1i
                   jj=jj+1
                   ! z component of point jj
                   iz=jj/(Glr%d%n2i*Glr%d%n1i)
                   ! Subtract the 'lower' xy layers
                   ii=jj-iz*(Glr%d%n2i*Glr%d%n1i)
                   ! y component of point jj
                   iy=ii/Glr%d%n1i
                   ! Subtract the 'lower' y rows
                   ii=ii-iy*Glr%d%n1i
                   ! x component
                   ix=ii

                   if(iy==ix0 .and. iz==iz0) write((iproc+1)*1000000+it*1000+iorb*10+7,*) ix, phir(jj)
                   ! Write along y-axis
                   if(ix==ix0 .and. iz==iz0) write((iproc+1)*1000000+it*1000+iorb*10+8,*) iy, phir(jj)
                   ! Write along z-axis
                   if(ix==ix0 .and. iy==iy0) write((iproc+1)*1000000+it*1000+iorb*10+9,*) iz, phir(jj)


                end do
            end do
        end do
        close(unit=(iproc+1)*1000000+it*1000+iorb*10+7)
        close(unit=(iproc+1)*1000000+it*1000+iorb*10+8)
        close(unit=(iproc+1)*1000000+it*1000+iorb*10+9)

        istart=istart+(Glr%wfd%nvctr_c+7*Glr%wfd%nvctr_f)*orbs%nspinor

    end do orbLoop

call deallocate_work_arrays_sumrho(w)
deallocate(phir, stat=istat)


end subroutine plotOrbitals







!!subroutine transformHam(iproc, nproc, orbs, comms, phi, hphi, HamSmall)
!!!
!!! Purpose:
!!! =======
!!!   Builds the Hamiltonian in the basis of the localized basis functions phi. To do so, it gets all basis
!!!   functions |phi_i> and H|phi_i> and then calculates H_{ij}=<phi_i|H|phi_j>. The basis functions phi are
!!!   provided in the transposed form.
!!!
!!! Calling arguments:
!!! ==================
!!!   Input arguments
!!!     HamLarge   Hamiltonian in large basis
!!!     phi        small basis set
!!!   Output arguments
!!!     HamSmall   Hamiltonian in small basis
!!!
!!use module_base
!!use module_types
!!implicit none
!!
!!! Calling arguments
!!integer:: iproc, nproc
!!type(orbitals_data), intent(inout) :: orbs
!!type(communications_arrays), intent(in) :: comms
!!real(8),dimension(sum(comms%nvctr_par(iproc,1:orbs%nkptsp))*orbs%nspinor,orbs%norb), intent(in) :: phi, hphi
!!real(8),dimension(orbs%norb,orbs%norb),intent(out):: HamSmall
!!
!!! Local variables
!!integer:: iorb, jorb, istat, ierr, nvctrp, iall
!!real(8),dimension(:,:),allocatable:: HamTemp
!!character(len=*),parameter:: subname='transformHam'
!!
!!
!!! Allocate a temporary array if there are several MPI processes
!!if(nproc>1) then
!!    allocate(HamTemp(orbs%norb,orbs%norb), stat=istat)
!!    call memocc(istat, orbs%norb*orbs%norb, 'HamTemp', subname)
!!end if
!!
!!! nvctrp is the amount of each phi hold by the current process
!!nvctrp=sum(comms%nvctr_par(iproc,1:orbs%nkptsp))*orbs%nspinor
!!
!!! Build the Hamiltonian. In the parallel case, each process writes its Hamiltonian in HamTemp
!!! and a mpi_allreduce sums up the contribution from all processes.
!!if(nproc==1) then
!!    call dgemm('t', 'n', orbs%norb, orbs%norb, nvctrp, 1.d0, phi(1,1), nvctrp, &
!!               hphi(1,1), nvctrp, 0.d0, HamSmall(1,1), orbs%norb)
!!else
!!    call dgemm('t', 'n', orbs%norb, orbs%norb, nvctrp, 1.d0, phi(1,1), nvctrp, &
!!               hphi(1,1), nvctrp, 0.d0, HamTemp(1,1), orbs%norb)
!!end if
!!if(nproc>1) then
!!    call mpi_allreduce(HamTemp(1,1), HamSmall(1,1), orbs%norb**2, mpi_double_precision, mpi_sum, mpi_comm_world, ierr)
!!end if
!!
!!if(nproc>1) then
!!   iall=-product(shape(HamTemp))*kind(HamTemp)
!!   deallocate(HamTemp,stat=istat)
!!   call memocc(istat, iall, 'HamTemp', subname)
!!end if
!!
!!end subroutine transformHam
!!
!!
!!
!!
!!
!!subroutine initializeParameters(iproc, nproc, Glr, orbs, orbsLIN, commsLIN, at, phi, input, rxyz, occupForInguess)
!!!
!!! Purpose:
!!! ========
!!!   This subroutine initializes all parameters needed for the linear scaling version.
!!!
!!! Calling arguments:
!!! ==================
!!!   Input arguments
!!!   Output arguments
!!!
!!use module_base
!!use module_types
!!use module_interfaces, except_this_one => initializeParameters
!!implicit none
!!
!!! Calling arguments
!!integer:: iproc, nproc
!!type(locreg_descriptors), intent(in) :: Glr
!!type(orbitals_data), intent(inout) :: orbs
!!type(orbitals_data), intent(out) :: orbsLIN
!!type(communications_arrays), intent(in) :: commsLIN
!!type(atoms_data), intent(in) :: at
!!real(8),dimension(:),allocatable:: phi
!!type(input_variables), intent(in) :: input
!!real(8),dimension(3,at%nat):: rxyz
!!real(8),dimension(32,at%nat):: occupForInguess
!!
!!! Local variables
!!integer:: ii, jproc, jj, istat, iorb, i, jorb, ierr, ii2, j, istart, jstart, iat, ityp, iall, norb_tot, iiOrb
!!real(8):: tt, ddot
!!integer,dimension(:),allocatable:: norb_par
!!real(8),dimension(:),pointer:: phiWorkPointer
!!integer,dimension(:),allocatable:: nbasisPerAt, nbasisArr, orbsPerAt
!!character(len=*),parameter:: subname='initializeParameters'
!!character(len=20):: atomname
!!
!!! Copy everything
!!orbsLin=orbs
!!
!!! Decide wheter quadratic or quartic potential
!!orbsLIN%power=4  ! quartic
!!!orbsLIN%power=2  ! quadratic
!!
!!if(iproc==0) then
!!    if(orbsLIN%power==2) then
!!        write(*,'(a)') 'The basis functions are created using a parabolic potential.'
!!    else if(orbsLIN%power==4) then
!!        write(*,'(a)') 'The basis functions are created using a quartic potential.'
!!    else
!!        write(*,'(a,i0,a)') 'ERROR: can only use parabolic or quartic potentials, but found power ', orbsLIN%power, '!' 
!!        stop
!!    end if
!!end if
!!
!!! Number of basis functions
!!allocate(nbasisPerAt(at%nat), stat=istat)
!!call memocc(istat, at%nat, 'nbasisPerAt', subname)
!!
!!! Read in the number of basis functions per atom type and save this information
!!! in the array nbasisPerAt.
!!allocate(orbsPerAt(at%ntypes), stat=istat)
!!allocate(orbsLIN%parabPrefacArr(at%ntypes), stat=istat)
!!open(unit=999, file='orbitalsValues')
!!read(999,*) orbsLIN%nItMin, orbsLIN%nItMax
!!read(999,*) orbsLIN%convCritInit, orbsLIN%convCritFinal
!!read(999,*) orbsLIN%DIISHistMin, orbsLIN%DIISHistMax
!!if(orbsLIN%DIISHistMin>orbsLIN%DIISHistMax) then
!!    if(iproc==0) write(*,'(a,i0,a,i0,a)') 'ERROR: DIISHistMin must not be larger than &
!!    & DIISHistMax, but you chose ', orbsLIN%DIISHistMin, ' and ', orbsLIN%DIISHistMax, '!'
!!    stop
!!end if
!!do iat=1,at%ntypes
!!    read(999,*) atomname, orbsPerAt(iat), orbsLIN%parabPrefacArr(iat)
!!    if(mod(orbsPerAt(iat),7)/=0) then
!!        !write(*,'(a,a,a,i0)') 'ERROR: the number of orbitals per atom must be a multiple of 7, but for ',trim(atomname),' we found ', orbsPerAt(iat)
!!        !stop
!!    end if
!!    if(iproc==0) write(*,'(a,a,a,i0,a,es9.3)') 'parameters for atom ', trim(atomname), &
!!        ': number of basis functions=', orbsPerAt(iat), ', prefactor=', orbsLIN%parabPrefacArr(iat)
!!    !do
!!    !    ishell=ii+occupForInguess(iat, ishell)
!!    !end do
!!    !occupForInguess
!!end do
!!close(unit=999) 
!!!open(unit=999, file='parabolaValues')
!!!do iat=1,at%ntypes
!!!    read(999,*) atomname, orbsLIN%parabPrefacArr(iat)
!!!    if(iproc==0) write(*,'(a,a,a,es12.4)') 'parabola value for ', trim(atomname),' is:', orbsLIN%parabPrefacArr(iat)
!!!end do
!!!close(unit=999)
!!
!!
!!! Count how many basis functions we have.
!!do iat=1,at%nat
!!    ityp=at%iatype(iat)
!!    nbasisPerAt(iat)=orbsPerAt(at%iatype(iat))
!!end do
!!orbs%nbasis=sum(nbasisPerAt) ! Choose such that they can be evenly distributed
!!deallocate(orbsPerAt, stat=istat)
!!
!!
!!! Determine how many orbitals shall be treated by each process.
!!! tt gives the ideal split. Unfortunately this will likely not be an integer.
!!tt=dble(orbs%nbasis)/dble(nproc)
!!orbs%nbasisp=floor(tt)
!!! Each process will handle floor(tt) orbitals.
!!! Since some orbitals are remaining if tt was not an integer, calculate the number of orbitals
!!! that some processors have to handle on top. This number is given by ii.
!!tt=tt-dble(floor(tt))
!!ii=int(tt*nproc)
!!do jproc=0,ii-1
!!    if(iproc==jproc) orbs%nbasisp=orbs%nbasisp+1
!!end do
!!if(ii==0) then
!!    if(iproc==0) write(*,'(a,2(i0,a))') 'Processes from 0 to ',nproc-1,' treat ',orbs%nbasisp,' orbitals.'
!!else
!!    if(iproc==0) write(*,'(a,5(i0,a))') 'Processes from 0 to ',ii-1,' treat ',orbs%nbasisp, &
!!        ' orbitals, processes from ',ii,' to ',nproc-1,' treat ',orbs%nbasisp-1,' orbitals.'
!!end if
!!
!!! The number of orbitals handled by each preocessor is stored in orbsLIN%norb_par.
!!! The local array norb_par is used to collect the data and is then deallocated
!!allocate(norb_par(0:nproc-1), stat=istat)
!!call memocc(istat, nproc, 'norb_par', subname)
!!norb_par=0
!!norb_par(iproc)=orbs%nbasisp
!!call mpi_allreduce(norb_par(0), orbsLIN%norb_par, nproc, mpi_integer, mpi_sum, mpi_comm_world, ierr)
!!
!!iall=-product(shape(norb_par))*kind(norb_par)
!!deallocate(norb_par, stat=istat)
!!call memocc(istat, iall, 'norb_par', subname)
!!
!!if(.not.allocated(orbsLIN%onWhichAtom)) allocate(orbsLIN%onWhichAtom(orbs%nbasisp), &
!!    stat=istat) ; if(istat/=0) stop 'ERROR in allocating onWhichAtom'
!!if(.not.allocated(orbsLIN%orbitalNumber)) allocate(orbsLIN%orbitalNumber(orbs%nbasisp), &
!!    stat=istat) ; if(istat/=0) stop 'ERROR in allocating onWhichAtom'
!!if(.not.allocated(orbsLIN%parabolaShift)) allocate(orbsLIN%parabolaShift(3,orbs%nbasisp), &
!!    stat=istat) ; if(istat/=0) stop 'ERROR in allocating onWhichAtom'
!!
!!
!!! Distribute the centers of the parabolic potential among the MPI processes.
!!! There are four counters:
!!!   jproc: indicates which MPI process is handling the basis function which is being treated
!!!   ii: counts the atom numbers
!!!   jj: counts the orbitals handled by a given process
!!!   iiOrb: counts the number of orbitals for a given atoms thas has already been assigned
!!!   ii2: counts the number of basis functions for a given atom which are on the 'previous' MPI process
!!!        and then stores this number in orbsLIN%nbasisOnPreviousMPI.
!!!        (example: MPI 0 handles 4 orbitals centered on atom i and MPI 1 handles the remaining 8, then
!!!         orbsLIN%nbasisOnPreviousMPI will have the value 4 for MPI 1)
!!jproc=0
!!ii=1
!!jj=0
!!ii2=0
!!iiOrb=0
!!orbsLIN%onWhichAtom=-1
!!orbsLIN%orbitalNumber=-1
!!!write(*,'(a,i2,3x,20i4)') 'iproc, nbasisPerAt', iproc, nbasisPerAt
!!
!!! THERE IS A PROBLEM WHEN ONLY 1 ORBITAL PER ATOM
!!do iorb=1,orbs%nbasis
!!call mpi_barrier(mpi_comm_world, ierr)
!!
!!    ! Switch to the next MPI process if the numbers of orbitals for a given
!!    ! MPI process is reached.
!!    if(jj==orbsLIN%norb_par(jproc)) then
!!        jproc=jproc+1
!!  !if(iproc==0) write(*,'(a,4i5)') 'iorb, jj, ii, ii2', iorb, jj, ii, ii2
!!        if(iproc==jproc .and. iiorb/=nbasisPerAt(max(ii,1))) orbsLIN%nbasisOnPreviousMPI=ii2
!!        jj=0
!!        ii2=0
!!    end if
!!    
!!    ! Switch to the next atom if the number of basis functions for this atom is reached.
!!    !if(mod(iorb,nbasisPerAt(max(ii,1)))==1 .or. nbasisPerAt(max(ii,1))==1) then
!!    !if((ii==0 .or. (ii2>1 .and. mod(ii2,nbasisPerAt(max(ii,1)))==1)) .or. nbasisPerAt(max(ii,1))==1) then
!!
!!    !if(iiOrb==nbasisPerAt(max(ii,1)) .or. nbasisPerAt(max(ii,1))==1) then
!!    if(iiOrb==nbasisPerAt(max(ii,1))) then
!!        ii=ii+1
!!        ii2=0
!!        iiOrb=0
!!    end if
!!    ii2=ii2+1
!!    jj=jj+1
!!    iiOrb=iiOrb+1
!!    !if(iproc==0) write(*,'(a,4i5)') 'iorb, ii, jj, jproc', iorb, ii, jj, jproc
!!    if(iproc==jproc) orbsLIN%onWhichAtom(jj)=ii
!!    if(iproc==jproc) orbsLIN%orbitalNumber(jj)=iiOrb
!!    !write(*,'(a,i2,i3,3x,20i4)') 'iproc, iorb, orbsLIN%onWhichAtom, orbsLIN%orbitalNumber', iproc, iorb, orbsLIN%onWhichAtom, orbsLIN%orbitalNumber
!!
!!end do    
!!orbsLIN%orbitalNumber=orbsLIN%orbitalNumber!+orbsLIN%nbasisOnPreviousMPI
!!do iorb=1,orbs%nbasisp
!!    !ii=mod(orbsLIN%orbitalNumber(iorb)-1,7)+1
!!    !tt=dble(ceiling(dble(orbsLIN%orbitalNumber(iorb))/7.d0))
!!    !tt=1.d0
!!    !select case(ii)
!!    !case(1)
!!    !    orbsLIN%parabolaShift(1,iorb)=0.d0
!!    !    orbsLIN%parabolaShift(2,iorb)=0.d0
!!    !    orbsLIN%parabolaShift(3,iorb)=0.d0
!!    !case(2)
!!    !    orbsLIN%parabolaShift(1,iorb)=-1.d-1*tt
!!    !    orbsLIN%parabolaShift(2,iorb)=0.d0
!!    !    orbsLIN%parabolaShift(3,iorb)=0.d0
!!    !case(3)
!!    !    orbsLIN%parabolaShift(1,iorb)=1.d-1*tt
!!    !    orbsLIN%parabolaShift(2,iorb)=0.d0
!!    !    orbsLIN%parabolaShift(3,iorb)=0.d0
!!    !case(4)
!!    !    orbsLIN%parabolaShift(1,iorb)=0.d0
!!    !    orbsLIN%parabolaShift(2,iorb)=-1.d-1*tt
!!    !    orbsLIN%parabolaShift(3,iorb)=0.d0
!!    !case(5)
!!    !    orbsLIN%parabolaShift(1,iorb)=0.d0
!!    !    orbsLIN%parabolaShift(2,iorb)=1.d-1*tt
!!    !    orbsLIN%parabolaShift(3,iorb)=0.d0
!!    !case(6)
!!    !    orbsLIN%parabolaShift(1,iorb)=0.d0
!!    !    orbsLIN%parabolaShift(2,iorb)=0.d0
!!    !    orbsLIN%parabolaShift(3,iorb)=-1.d-1*tt
!!    !case(7)
!!    !    orbsLIN%parabolaShift(1,iorb)=0.d0
!!    !    orbsLIN%parabolaShift(2,iorb)=0.d0
!!    !    orbsLIN%parabolaShift(3,iorb)=1.d-1*tt
!!    !end select
!!    orbsLIN%parabolaShift(1,iorb)=0.d0
!!    orbsLIN%parabolaShift(2,iorb)=0.d0
!!    orbsLIN%parabolaShift(3,iorb)=0.d0
!!!write(*,'(a,2i5,3es9.2)') 'iproc, iorb, orbsLIN%parabolaShift(:,iorb)', iproc, iorb, orbsLIN%parabolaShift(:,iorb)
!!end do
!!
!!!write(*,'(a,i2,3x,200i4)') 'iproc, orbsLIN%onWhichAtom', iproc, orbsLIN%onWhichAtom
!!!write(*,'(a,i2,3x,200i4)') 'iproc, orbsLIN%orbitalNumber', iproc, orbsLIN%orbitalNumber
!!!write(*,'(a,i2,3x,200i4)') 'iproc, orbsLIN%nbasisOnPreviousMPI', iproc, orbsLIN%nbasisOnPreviousMPI
!!iall=-product(shape(nbasisPerAt))*kind(nbasisPerAt)
!!deallocate(nbasisPerAt, stat=istat)
!!call memocc(istat, iall, 'nbasisPerAt', subname)
!!
!!
!!! Now adapt a few parameters in orbsLin%norb
!!
!!! The number of orbitals in orbsLIN is the number of basis functions
!!orbsLin%norb=orbs%nbasis
!!! At the moment only working for closed shell!
!!orbsLin%norbu=orbs%nbasis
!!orbsLin%norbp=orbs%nbasisp
!!if(orbsLIN%norbd>0) stop 'ERROR: not yet implemented for nspin=2!'
!!
!!! READ IN THE PREFACTOR FOR THE PARABOLIC POTENTIAL
!!!open(unit=999, file='parabPrefac')
!!!read(999,*) orbsLIN%parabPrefac
!!!read(999,*) orbsLIN%parabMaxVal
!!!close(unit=999)
!!
!!! Read in the prefactors for the parabola centered at the atoms. Each atom type
!!! has it own value.
!!!allocate(orbsLIN%parabPrefacArr(at%ntypes), stat=istat)
!!!open(unit=999, file='parabolaValues')
!!!do iat=1,at%ntypes
!!!    read(999,*) atomname, orbsLIN%parabPrefacArr(iat)
!!!    if(iproc==0) write(*,'(a,a,a,es12.4)') 'parabola value for ', trim(atomname),' is:', orbsLIN%parabPrefacArr(iat)
!!!end do
!!!close(unit=999)
!!
!!! orbsLIN%isorb is the 'first' orbital for a given MPI process.
!!norb_tot=0
!!do jproc=0,iproc-1
!!   norb_tot=norb_tot+orbsLIN%norb_par(jproc)
!!end do
!!!reference orbital for process
!!orbsLIN%isorb=norb_tot
!!
!!
!!! I don't know what this means, but it seems that it is necessary.
!!nullify(orbsLIN%spinsgn)
!!allocate(orbsLIN%spinsgn(orbs%nbasis), stat=istat)
!!orbsLIN%spinsgn=1.d0  ! WHY LIKE THIS?
!!
!!! I don't know what this means, but it seems that it is necessary.
!!nullify(orbsLIN%iokpt)
!!allocate(orbsLIN%iokpt(orbs%nbasisp), stat=istat)
!!orbsLIN%iokpt=1 ! WHY LIKE THIS?
!!
!!! I don't know what this means, but it seems that it is necessary.
!!nullify(orbsLIN%occup)
!!allocate(orbsLIN%occup(orbs%nbasis), stat=istat)
!!orbsLIN%occup=2.d0  ! WHY LIKE THIS?
!!
!!! The eigenvalues from the input guess. They are used for the preconditioning.
!!nullify(orbsLIN%eval)
!!allocate(orbsLIN%eval(orbs%nbasis), stat=istat)
!!orbsLIN%eval=-.5d0
!!
!!
!!! Assign the parameters needed for the communication to commsLIN.
!!call orbitals_communicators(iproc,nproc,Glr,orbsLIN,commsLIN)
!!
!!
!!! Allocate phi and initialize it at random
!!allocate(phi(orbsLIN%npsidim), stat=istat)
!!call initRandomSeed(iproc, 1)
!!call random_number(phi)
!!
!!!write(*,*) 'calling createInputGuess'
!!!call createInputGuess(iproc, orbsLIN, Glr, input, at, rxyz, phi)
!!
!!! Orthonormalize phi.
!!allocate(phiWorkPointer(size(phi)), stat=istat)
!!call transpose_v(iproc, nproc, orbsLIN, Glr%wfd, commsLIN, phi, work=phiWorkPointer)
!!call orthogonalize(iproc, nproc, orbsLIN, commsLIN, Glr%wfd, phi, input)
!!call untranspose_v(iproc, nproc, orbsLIN, Glr%wfd, commsLIN, phi, work=phiWorkPointer)
!!deallocate(phiWorkPointer, stat=istat)
!!
!!
!!end subroutine initializeParameters
!!
!!
!!
!!
!!
!!
!!
!!subroutine diagonalizeHamiltonian(iproc, nproc, orbsLIN, HamSmall, eval)
!!!
!!! Purpose:
!!! ========
!!!   Diagonalizes the Hamiltonian HamSmall and makes sure that all MPI processes give
!!!   the same result. This is done by requiring that the first entry of each vector
!!!   is positive.
!!!
!!! Calling arguments:
!!! ==================
!!!   Input arguments
!!!     iproc     process ID
!!!     nproc     number of MPI processes
!!!     orbsLIN   type containing many parameters
!!!   Input / Putput arguments
!!!     HamSmall  on input: the Hamiltonian
!!!               on exit: the eigenvectors
!!!   Output arguments
!!!     eval      the associated eigenvalues 
!!!
!!use module_base
!!use module_types
!!implicit none
!!
!!! Calling arguments
!!integer:: iproc, nproc
!!type(orbitals_data), intent(inout) :: orbsLIN
!!real(8),dimension(orbsLIN%norb, orbsLIN%norb):: HamSmall
!!real(8),dimension(orbsLIN%norb):: eval
!!
!!! Local variables
!!integer:: lwork, info, istat, i, iorb, jorb
!!real(8),dimension(:),allocatable:: work
!!
!!
!!! Diagonalize the Hamiltonian 
!!lwork=-1 
!!allocate(work(1), stat=istat) ; if(istat/=0) stop 'ERROR in allocating work' 
!!call dsyev('v', 'l', orbsLIN%norb, HamSmall(1,1), orbsLIN%norb, eval(1), work(1), lwork, info) 
!!lwork=work(1) 
!!deallocate(work, stat=istat) ; if(istat/=0) stop 'ERROR in deallocating work' 
!!allocate(work(lwork), stat=istat) ; if(istat/=0) stop 'ERROR in allocating work' 
!!call dsyev('v', 'l', orbsLIN%norb, HamSmall(1,1), orbsLIN%norb, eval(1), work(1), lwork, info) 
!!
!!! Make sure that the eigenvectors are the same for all MPI processes. To do so, require that 
!!! the first entry of each vector is positive.
!!do iorb=1,orbsLIN%norb
!!    if(HamSmall(1,iorb)<0.d0) then
!!        do jorb=1,orbsLIN%norb
!!            HamSmall(jorb,iorb)=-HamSmall(jorb,iorb)
!!        end do
!!    end if
!!end do
!!
!!
!!!! Write the eigenvalues.
!!!if(iproc==0) write(*,'(/,a)') 'The eigenvalues:'
!!!do i=1,orbsLIN%norb
!!!    !if(i==p%norb) then
!!!    !    message=' <-- HOMO'
!!!    !else if(i==p%norb+1) then
!!!    !    message=' <-- LUMO'
!!!    !else
!!!    !    message=''
!!!    !end if
!!!    !write(*,'(a,i0,a,es10.3,a)') 'eval(',i,') = ', eval(i), message
!!!    if(iproc==0) write(*,'(a,i0,a,es10.3)') 'eval(',i,') = ', eval(i)
!!!end do
!!
!!
!!end subroutine diagonalizeHamiltonian
!!
!!
!!
!!
!!
!!subroutine buildWavefunction(iproc, nproc, orbs, orbsLIN, comms, commsLIN, phiOld, phiNew, HamSmall)
!!
!!use module_base
!!use module_types
!!implicit none
!!
!!! Calling arguments
!!integer:: iproc, nproc
!!type(orbitals_data), intent(in) :: orbs
!!type(orbitals_data), intent(in) :: orbsLIN
!!type(communications_arrays), intent(in) :: comms
!!type(communications_arrays), intent(in) :: commsLIN
!!real(8),dimension(sum(commsLIN%nvctr_par(iproc,1:orbsLIN%nkptsp))*orbsLIN%nspinor,orbsLIN%norb) :: phiOld
!!!real(8),dimension(sum(commsLIN%nvctr_par(iproc,1:orbsLIN%nkptsp))*orbsLIN%nspinor,orbs%norb) :: phiNew
!!real(8),dimension(sum(comms%nvctr_par(iproc,1:orbs%nkptsp))*orbs%nspinor,orbs%norb) :: phiNew
!!real(8),dimension(orbsLIN%norb,orbsLIN%norb):: HamSmall
!!
!!! Local variables
!!integer:: nvctrp
!!
!!
!!! Is this necessary??
!!phiNew=0.d0
!!
!!!nvctrp=sum(commsLIN%nvctr_par(iproc,1:orbsLIN%nkptsp))*orbsLIN%nspinor
!!nvctrp=sum(comms%nvctr_par(iproc,1:orbs%nkptsp))*orbs%nspinor
!!call dgemm('n', 'n', nvctrp, orbs%norb, orbsLIN%norb, 1.d0, phiOld(1,1), nvctrp, HamSmall(1,1), &
!!           orbsLIN%norb, 0.d0, phiNew(1,1), nvctrp)
!!
!!
!!
!!end subroutine buildWavefunction





!!!subroutine improveOrbitals(iproc, nproc, nspin, Glr, orbs, orbsLIN, comms, commsLIN, at, rxyz, rxyzParab, &
!!!    nscatterarr, ngatherarr, nlpspd, proj, sizeRhopot, rhopot, GPU, input, pkernelseq, phi, psi, psit, &
!!!    iter, infoBasisFunctions, n3p, pulayAt, pulayDir, shift, ebs_mod)
!!!!
!!!! Purpose:
!!!! ========
!!!!   Improves the eigenvectors according to the updated electronic density.
!!!!
!!!! Calling arguments:
!!!! ==================
!!!!
!!!use module_base
!!!use module_types
!!!use module_interfaces, except_this_one => improveOrbitals
!!!implicit none
!!!
!!!! Calling arguments
!!!integer:: iproc, nproc, nspin, sizeRhopot, infoBasisFunctions
!!!type(locreg_descriptors), intent(in) :: Glr
!!!type(orbitals_data), intent(inout) :: orbs, orbsLIN
!!!type(communications_arrays), intent(in) :: comms
!!!type(communications_arrays), intent(in) :: commsLIN
!!!type(atoms_data), intent(in) :: at
!!!real(8),dimension(3,at%nat):: rxyz, rxyzParab
!!!integer, dimension(0:nproc-1,4), intent(in) :: nscatterarr !n3d,n3p,i3s+i3xcsh-1,i3xcsh
!!!integer, dimension(0:nproc-1,2), intent(in) :: ngatherarr 
!!!type(nonlocal_psp_descriptors), intent(in) :: nlpspd
!!!real(wp), dimension(nlpspd%nprojel), intent(in) :: proj
!!!!real(dp), dimension(*), intent(inout) :: rhopot
!!!real(dp), dimension(sizeRhopot), intent(inout) :: rhopot
!!!type(GPU_pointers), intent(inout) :: GPU
!!!type(input_variables):: input
!!!real(dp), dimension(:), pointer :: pkernelseq
!!!real(8),dimension(orbsLIN%npsidim):: phi
!!!real(8),dimension(orbs%npsidim):: psi, psit
!!!integer:: iter, n3p
!!!integer,optional:: pulayAt, pulayDir
!!!real(8),optional:: shift, ebs_mod
!!!
!!!! Local variables
!!!integer:: istat, i, j, istart, jstart, ierr
!!!real(8),dimension(:),allocatable:: hphi, eval
!!!real(8),dimension(:,:),allocatable:: HamSmall
!!!real(8),dimension(:,:,:),allocatable:: matrixElements
!!!real(8),dimension(:),pointer:: phiWorkPointer
!!!real(8):: epot_sum,ekin_sum,eexctX,eproj_sum, ddot, trace, dnrm2, tt
!!!character(len=1)::num,num2
!!!character(len=20):: filename
!!!logical:: modifiedEnergy
!!!real(8),dimension(:),allocatable:: phiOld, psiOld
!!!integer:: iorb, jorb, korb, nvctrp, idsx, idsxStart, idsxMin, idsxMax
!!!real(wp), dimension(:), pointer :: potential
!!!character(len=*),parameter:: subname='improveOrbitals'
!!!
!!!allocate(hphi(orbsLIN%npsidim), stat=istat)
!!!
!!!
!!!!if(iter<=0) then
!!!!    write(num,'(i1)') iproc
!!!!    filename='phiIproc'//num
!!!!    open(unit=iproc+1, file=trim(filename))
!!!!    do i=1,size(phi)
!!!!        !write(iproc+1,*) phi(i)
!!!!        read(iproc+1,*) phi(i)
!!!!    end do
!!!!    close(unit=iproc+1)
!!!!else 
!!!    !call random_number(phi)
!!!if(iter>=0) then
!!!    !read(3000000+iproc,*) phi
!!!    !read(1000000+iproc,*) phi
!!!    !read(5000000+iproc,*) phi
!!!    call getLocalizedBasis(iproc, nproc, at, orbs, Glr, input, orbsLIN, commsLIN, rxyz, nspin, nlpspd, proj, &
!!!        nscatterarr, ngatherarr, rhopot, GPU, pkernelseq, phi, hphi, trace, rxyzParab, &
!!!        orbsLIN%DIISHistMin, orbsLIN%DIISHistMax, infoBasisFunctions)
!!!    !write(4000000+iproc,*) phi
!!!    !write(5000000+iproc,*) phi
!!!else
!!!    ! idsx is the maximal DIIS history length, idsxStart is the history length with
!!!    ! which DIIS start (0 means SD). If everything works fine, the program switches
!!!    ! automatically to DIIS with histoy length idsx.
!!!    !if(.not.present(pulayAt) .or. .not.present(pulayDir) .or. .not.present(shift)) stop 'pulayAt and/or pulayDir and/or shift not present!'
!!!    allocate(phiOld(size(phi)))
!!!    phiOld=phi
!!!
!!!    !!&! Take the 'old' HamSmall
!!!    !!&call HamiltonianApplication(iproc,nproc,at,orbsLIN,input%hx,input%hy,input%hz,rxyz,&
!!!    !!&     nlpspd,proj,Glr,ngatherarr,Glr%d%n1i*Glr%d%n2i*nscatterarr(iproc,2),&
!!!    !!&     rhopot(1),&
!!!    !!&     phi(1),hphi(1),ekin_sum,epot_sum,eexctX,eproj_sum,nspin,GPU,pkernel=pkernelseq)
!!!    !!&allocate(phiWorkPointer(max(size(phi),size(psi))), stat=istat)
!!!    !!&call transpose_v(iproc, nproc, orbsLIN, Glr%wfd, commsLIN, phi, work=phiWorkPointer)
!!!    !!&call transpose_v(iproc, nproc, orbsLIN, Glr%wfd, commsLIN, hphi, work=phiWorkPointer)
!!!    !!&deallocate(phiWorkPointer, stat=istat)
!!!    !!&allocate(HamSmall(orbsLIN%norb,orbsLIN%norb), stat=istat)
!!!    !!&call transformHam(iproc, nproc, orbsLIN, commsLIN, phi, hphi, HamSmall)
!!!    !!&allocate(eval(orbsLIN%norb), stat=istat)
!!!    !!&call diagonalizeHamiltonian(iproc, nproc, orbsLIN, HamSmall, eval)
!!!    !!&allocate(phiWorkPointer(max(size(phi),size(psi))), stat=istat)
!!!    !!&call untranspose_v(iproc, nproc, orbsLIN, Glr%wfd, commsLIN, phi, work=phiWorkPointer)
!!!    !!&call untranspose_v(iproc, nproc, orbsLIN, Glr%wfd, commsLIN, hphi, work=phiWorkPointer)
!!!    !!&deallocate(phiWorkPointer, stat=istat)
!!!
!!!    
!!!    !call pulayNew(iproc, nproc, at, orbs, Glr, input, orbsLIN, commsLIN, rxyz, nspin, nlpspd, proj, nscatterarr, ngatherarr, rhopot, GPU, pkernelseq, phi, hphi, rxyzParab, pulayAt, pulayDir, shift)
!!!    call getLocalizedBasis(iproc, nproc, at, orbs, Glr, input, orbsLIN, commsLIN, rxyz, &
!!!        nspin, nlpspd, proj, nscatterarr, ngatherarr, rhopot, GPU, pkernelseq, phi, &
!!!        hphi, trace, rxyzParab, orbsLIN%DIISHistMin, orbsLIN%DIISHistMax, infoBasisFunctions)
!!!    !istart=1
!!!    !do iorb=1,orbsLIN%norbp
!!!    !    write(*,'(a,i4,i3,i5,2es22.12)') 'pulayAt, iproc, iorb, <phi|phiOld>', pulayAt, iproc, iorb, ddot(Glr%wfd%nvctr_c+7*Glr%wfd%nvctr_f, phi(istart), 1, phiOld(istart), 1)
!!!    !    istart=istart+Glr%wfd%nvctr_c+7*Glr%wfd%nvctr_f
!!!    !end do
!!!
!!!    !!&allocate(phiWorkPointer(max(size(phi),size(psi))), stat=istat)
!!!    !!&call transpose_v(iproc, nproc, orbsLIN, Glr%wfd, commsLIN, phi, work=phiWorkPointer)
!!!    !!&call transpose_v(iproc, nproc, orbsLIN, Glr%wfd, commsLIN, hphi, work=phiWorkPointer)
!!!    !!&deallocate(phiWorkPointer, stat=istat)
!!!end if
!!!    !call getLocalizedBasis2(iproc, nproc, at, orbs, Glr, input, orbsLIN, commsLIN, rxyz, nspin, nlpspd, proj, nscatterarr, ngatherarr, rhopot, GPU, pkernelseq, phi, hphi, trace, rxyzParab)
!!!!    if(iter<=9) then
!!!!        write(num,'(i1)') iproc
!!!!        write(num2,'(i1)') iter
!!!!        filename='phiIproc'//num//'_'//num2
!!!!        open(unit=iproc+1, file=trim(filename))
!!!!        do i=1,size(phi)
!!!!            write(iproc+1,*) phi(i)
!!!!        end do
!!!!        close(unit=iproc+1)
!!!!    end if
!!!!end if
!!!
!!!!! Estimate the pulay forces by finite differences
!!!!do iat=1,at%nat
!!!!    do icomp=1,3
!!!!        ! shift the parabolic center
!!!!    end do
!!!!end do
!!!
!!!
!!!!call HamiltonianApplicationParabola(iproc,nproc,at,orbsLIN,input%hx,input%hy,input%hz,rxyz,&
!!!!     nlpspd,proj,Glr,ngatherarr,Glr%d%n1i*Glr%d%n2i*nscatterarr(iproc,2),&
!!!!     rhopot(1),&
!!!!     phi(1),hphi(1),ekin_sum,epot_sum,eexctX,eproj_sum,nspin,GPU, orbsLIN%onWhichAtom, pkernel=pkernelseq)
!!!!if(iter>=0) then
!!!
!!!
!!!modifiedEnergy=.true.
!!!if(modifiedEnergy) then
!!!    ! Calculate <phi_i|H|phi_j>
!!!    allocate(matrixElements(orbsLIN%norb, orbsLIN%norb,2))
!!!    call HamiltonianApplicationParabola(iproc,nproc,at,orbsLIN,input%hx,input%hy,input%hz,rxyz,&
!!!         nlpspd,proj,Glr,ngatherarr,Glr%d%n1i*Glr%d%n2i*nscatterarr(iproc,2),&
!!!         rhopot(1),&
!!!         phi(1),hphi(1),ekin_sum,epot_sum,eexctX,eproj_sum,nspin,GPU, orbsLIN%onWhichAtom, rxyzParab, pkernel=pkernelseq)
!!!    allocate(phiWorkPointer(max(size(phi),size(psi))), stat=istat)
!!!    call transpose_v(iproc, nproc, orbsLIN, Glr%wfd, commsLIN, phi, work=phiWorkPointer)
!!!    call transpose_v(iproc, nproc, orbsLIN, Glr%wfd, commsLIN, hphi, work=phiWorkPointer)
!!!    nvctrp=sum(commsLIN%nvctr_par(iproc,1:orbsLIN%nkptsp))*orbs%nspinor
!!!    jstart=1
!!!    do jorb=1,orbsLIN%norb
!!!        istart=1
!!!        do iorb=1,orbsLIN%norb
!!!            matrixElements(iorb,jorb,2)=ddot(nvctrp, phi(istart), 1, hphi(jstart), 1)
!!!            istart=istart+nvctrp
!!!        end do
!!!        jstart=jstart+nvctrp
!!!    end do
!!!    call mpi_allreduce(matrixElements(1,1,2), matrixElements(1,1,1), orbsLIN%norb**2, &
!!!        mpi_double_precision, mpi_sum, mpi_comm_world, ierr)
!!!    !if(iproc==0) then
!!!    !    write(*,*) 'matrix Elements'
!!!    !    do iorb=1,orbsLIN%norb
!!!    !        write(*,'(80es9.2)') (matrixElements(iorb,jorb,1), jorb=1,orbsLIN%norb)
!!!    !    end do
!!!    !end if
!!!    call untranspose_v(iproc, nproc, orbsLIN, Glr%wfd, commsLIN, phi, work=phiWorkPointer)
!!!    call untranspose_v(iproc, nproc, orbsLIN, Glr%wfd, commsLIN, hphi, work=phiWorkPointer)
!!!    deallocate(phiWorkPointer)
!!!end if
!!!
!!!if(.true.) then
!!!    if(iproc==0) write(*,'(a)', advance='no') 'calculation of physical orbitals... '
!!!
!!!    !allocate the potential in the full box
!!!    call full_local_potential(iproc,nproc,Glr%d%n1i*Glr%d%n2i*n3p,Glr%d%n1i*Glr%d%n2i*Glr%d%n3i,input%nspin,&
!!!         orbs%norb,orbs%norbp,ngatherarr,rhopot,potential)
!!!
!!!    !call HamiltonianApplication(iproc,nproc,atoms,orbs,hx,hy,hz,rxyz,&
!!!    !     nlpspd,proj,Glr,ngatherarr,potential,psi,hpsi,ekin_sum,epot_sum,eexctX,eproj_sum,&
!!!    !     in%nspin,GPU,pkernel=pkernelseq)
!!!
!!!    call HamiltonianApplication(iproc,nproc,at,orbsLIN,input%hx,input%hy,input%hz,rxyz,&
!!!         nlpspd,proj,Glr,ngatherarr,potential,&
!!!         phi(1),hphi(1),ekin_sum,epot_sum,eexctX,eproj_sum,nspin,GPU,pkernel=pkernelseq)
!!!
!!!    !call HamiltonianApplication(iproc,nproc,at,orbsLIN,input%hx,input%hy,input%hz,rxyz,&
!!!    !     nlpspd,proj,Glr,ngatherarr,Glr%d%n1i*Glr%d%n2i*nscatterarr(iproc,2),&
!!!    !     rhopot(1),&
!!!    !     phi(1),hphi(1),ekin_sum,epot_sum,eexctX,eproj_sum,nspin,GPU,pkernel=pkernelseq)
!!!
!!!    !deallocate potential
!!!    call free_full_potential(nproc,potential,subname)
!!!
!!!
!!!
!!!    allocate(phiWorkPointer(max(size(phi),size(psi))), stat=istat)
!!!    call transpose_v(iproc, nproc, orbsLIN, Glr%wfd, commsLIN, phi, work=phiWorkPointer)
!!!    call transpose_v(iproc, nproc, orbsLIN, Glr%wfd, commsLIN, hphi, work=phiWorkPointer)
!!!    deallocate(phiWorkPointer, stat=istat)
!!!    allocate(HamSmall(orbsLIN%norb,orbsLIN%norb), stat=istat)
!!!    call transformHam(iproc, nproc, orbsLIN, commsLIN, phi, hphi, HamSmall)
!!!    !if(iter>=0) then
!!!    !    do iorb=1,orbsLIN%norb
!!!    !        write(110+iproc,'(100f7.3)') (HamSmall(iorb,jorb), jorb=1,orbsLIN%norb)
!!!    !    end do
!!!    !else
!!!    !    do iorb=1,orbsLIN%norb
!!!    !        write(120+iproc,'(100f7.3)') (HamSmall(iorb,jorb), jorb=1,orbsLIN%norb)
!!!    !    end do
!!!    !end if
!!!    if(iproc==0) write(*,'(a)', advance='no') 'Linear Algebra... '
!!!    allocate(eval(orbsLIN%norb), stat=istat)
!!!    call diagonalizeHamiltonian(iproc, nproc, orbsLIN, HamSmall, eval)
!!!    !if(iproc==0) then
!!!    !    write(*,*) 'diagonalized HamSmall'
!!!    !    do iorb=1,orbsLIN%norb
!!!    !        write(*,'(80es9.2)') (HamSmall(iorb,jorb), jorb=1,orbsLIN%norb)
!!!    !    end do
!!!    !end if
!!!end if
!!!
!!!! Calculate the modified band structure energy
!!!if(modifiedEnergy) then
!!!    tt=0.d0
!!!    do iorb=1,orbs%norb
!!!        do jorb=1,orbsLIN%norb
!!!            do korb=1,orbsLIN%norb
!!!                tt=tt+HamSmall(korb,iorb)*HamSmall(jorb,iorb)*matrixElements(korb,jorb,1)
!!! !if(iproc==0) write(*,'(a,3i5,3es12.4,es16.5)') 'iorb, jorb, korb, HamSmall(korb,iorb), HamSmall(jorb,iorb), matrixElements(korb,jorb,1), tt', iorb, jorb, korb, HamSmall(korb,iorb), HamSmall(jorb,iorb), matrixElements(korb,jorb,1), tt
!!!            end do
!!!        end do
!!!    end do
!!!    if(present(ebs_mod)) then
!!!        if(nspin==1) ebs_mod=2.d0*tt ! 2 for closed shell
!!!    end if
!!!end if
!!!
!!!!write(*,*) 'after diagonalizeHamiltonian, iproc',iproc
!!!
!!!!if(iter>=0) then
!!!!    do iorb=1,orbsLIN%norb
!!!!        !write(80+iproc,'(100f7.3)') (HamSmall(iorb,jorb), jorb=1,orbsLIN%norb)
!!!!    end do
!!!!else
!!!!    do iorb=1,orbsLIN%norb
!!!!        !write(90+iproc,'(100f7.3)') (HamSmall(iorb,jorb), jorb=1,orbsLIN%norb)
!!!!    end do
!!!!end if
!!!nvctrp=sum(commsLIN%nvctr_par(iproc,1:orbsLIN%nkptsp))*orbs%nspinor
!!!!if(iter>=0) then
!!!!    do istat=1,orbsLIN%norb*nvctrp
!!!!        if(mod(istat-1,nvctrp)==0) write(130+iproc,'(i0,a)') nint(dble(istat)/dble(nvctrp))+1,' ---------------------------'
!!!!        !write(130+iproc,*) phi(istat)
!!!!    end do
!!!!else
!!!!    do istat=1,orbsLIN%norb*nvctrp
!!!!        if(mod(istat-1,nvctrp)==0) write(140+iproc,'(i0,a)') nint(dble(istat)/dble(nvctrp))+1,' ---------------------------'
!!!!        !write(140+iproc,*) phi(istat)
!!!!    end do
!!!!end if
!!!
!!!
!!!! Store the new wave function in psi
!!!allocate(psiOld(size(psi)))
!!!psiOld=psi
!!!allocate(phiWorkPointer(max(size(phi),size(psi))), stat=istat)
!!!call untranspose_v(iproc, nproc, orbs, Glr%wfd, comms, psiOld, work=phiWorkPointer)
!!!deallocate(phiWorkPointer)
!!!
!!!! THIS IS A TEST
!!!psiOld=0.d0
!!!nvctrp=sum(comms%nvctr_par(iproc,1:orbs%nkptsp))*orbs%nspinor
!!!jstart=1
!!!do jorb=1,orbs%norb
!!!    istart=1
!!!    do iorb=1,orbsLIN%norb
!!!        call daxpy(nvctrp, HamSmall(iorb,jorb), phi(istart), 1, psiOld(jstart), 1)
!!!    istart=istart+nvctrp
!!!    end do
!!!    jstart=jstart+nvctrp
!!!end do
!!!!if(iter>=0) then
!!!!    do istat=1,orbs%norb*nvctrp
!!!!        if(mod(istat-1,nvctrp)==0) write(150+iproc,'(i0,a)') nint(dble(istat)/dble(nvctrp))+1,' ---------------------------'
!!!!        write(150+iproc,*) psiOld(istat)
!!!!    end do
!!!!else
!!!!    do istat=1,orbs%norb*nvctrp
!!!!        if(mod(istat-1,nvctrp)==0) write(160+iproc,'(i0,a)') nint(dble(istat)/dble(nvctrp))+1,' ---------------------------'
!!!!        write(160+iproc,*) psiOld(istat)
!!!!    end do
!!!!end if
!!!
!!!call buildWavefunction(iproc, nproc, orbs, orbsLIN, comms, commsLIN, phi, psi, HamSmall)
!!!nvctrp=sum(comms%nvctr_par(iproc,1:orbs%nkptsp))*orbs%nspinor
!!!!if(iter>=0) then
!!!!    do istat=1,orbs%norb*nvctrp
!!!!        if(mod(istat-1,nvctrp)==0) write(40+iproc,'(i0,a)') nint(dble(istat)/dble(nvctrp))+1,' ---------------------------'
!!!!        write(40+iproc,*) psi(istat)
!!!!    end do
!!!!else
!!!!    do istat=1,orbs%norb*nvctrp
!!!!        if(mod(istat-1,nvctrp)==0) write(50+iproc,'(i0,a)') nint(dble(istat)/dble(nvctrp))+1,' ---------------------------'
!!!!        write(50+iproc,*) psi(istat)
!!!!    end do
!!!! !call mpi_barrier(mpi_comm_world, iorb)
!!!! !stop
!!!!end if
!!!!istart=1
!!!!do iorb=1,orbs%norbp
!!!!    !write(*,'(a,i4,i3,i5,2es12.5)') 'pulayAt, iproc, iorb, <psi|psiOld>', pulayAt, iproc, iorb, ddot(Glr%wfd%nvctr_c+7*Glr%wfd%nvctr_f, psi(istart), 1, psiOld(istart), 1)
!!!!    !write(*,'(a,i4,i3,i5,2es12.5)') 'pulayAt, iproc, iorb, <psi|psi>', pulayAt, iproc, iorb, ddot(Glr%wfd%nvctr_c+7*Glr%wfd%nvctr_f, psi(istart), 1, psi(istart), 1)
!!!!    write(*,'(a,i4,i3,i5,2es12.5)') 'pulayAt, iproc, iorb, <psiOld|psiOld>', pulayAt, iproc, iorb, ddot(Glr%wfd%nvctr_c+7*Glr%wfd%nvctr_f, psiOld(istart), 1, psiOld(istart), 1)
!!!!    istart=istart+Glr%wfd%nvctr_c+7*Glr%wfd%nvctr_f
!!!!end do
!!!call dcopy(orbs%npsidim, psi, 1, psit, 1)
!!!if(iproc==0) write(*,'(a)') 'done.'
!!!
!!!!!! Get the Pulay forces
!!!!write(*,*) 'iproc, iter', iproc, iter
!!!!write(*,*) 'before calling pulay, iproc', iproc
!!!!call mpi_barrier(mpi_comm_world, i)
!!!!if(iter>=0) call pulay(iproc, nproc, Glr, orbs, orbsLIN, comms, commsLIN, input, at, rxyz, phi, hphi, psi)
!!!!!write(*,*) 'ATTENTION --  FIRST TEST!'
!!!!!call HamiltonianApplicationParabola(iproc,nproc,at,orbsLIN,input%hx,input%hy,input%hz,rxyz,&
!!!!!     nlpspd,proj,Glr,ngatherarr,Glr%d%n1i*Glr%d%n2i*nscatterarr(iproc,2),&
!!!!!     rhopot(1),&
!!!!!     phi(1),hphi(1),ekin_sum,epot_sum,eexctX,eproj_sum,nspin,GPU, orbsLIN%onWhichAtom, rxyz, pkernel=pkernelseq)
!!!!call pulay(iproc, nproc, Glr, orbs, orbsLIN, comms, commsLIN, input, at, rxyz, phi, hphi, psi, nscatterarr, ngatherarr, nlpspd, proj, sizeRhopot, rhopot, GPU, pkernelseq)
!!!!call mpi_barrier(mpi_comm_world, i)
!!!
!!!! Orthonormalization necessary?
!!!!call initRandomSeed(iproc, 1)
!!!!call random_number(psi)
!!!!call transpose_v(iproc, nproc, orbs, Glr%wfd, comms, psi, work=phiWorkPointer)
!!!!call orthogonalize(iproc, nproc, orbs, comms, Glr%wfd, psi, input)
!!!
!!!allocate(phiWorkPointer(max(size(phi),size(psi))), stat=istat)
!!!call untranspose_v(iproc, nproc, orbsLIN, Glr%wfd, commsLIN, phi, work=phiWorkPointer)
!!!call untranspose_v(iproc, nproc, orbs, Glr%wfd, comms, psi, work=phiWorkPointer)
!!!! Not necessary to untranspose hphi (no longer needed)
!!!!call untranspose_v(iproc, nproc, orbsLIN, Glr%wfd, commsLIN, hphi, work=phiWorkPointer)
!!!
!!!!if(iter>=0) then
!!!!    do istat=1,size(psi)
!!!!        write(60+iproc,*) psi(istat)
!!!!    end do
!!!!    write(*,*) iproc, 'writes', size(psi)
!!!!else
!!!!    do istat=1,size(psi)
!!!!        write(70+iproc,*) psi(istat)
!!!!    end do
!!!!end if
!!!
!!!!allocate(psiOld(size(psi)))
!!!!psiOld=psi
!!!!istart=1
!!!!write(*,'(a,2i6,2i12)') 'iproc, istart, istart+Glr%wfd%nvctr_c+7*Glr%wfd%nvctr_f, size(psiOld)', iproc, istart, istart+Glr%wfd%nvctr_c+7*Glr%wfd%nvctr_f, size(psiOld)
!!!!do iorb=1,orbs%norbp
!!!!    !write(*,'(a,i4,i3,i5,2es12.5)') 'pulayAt, iproc, iorb, <psi|psiOld>', pulayAt, iproc, iorb, ddot(Glr%wfd%nvctr_c+7*Glr%wfd%nvctr_f, psi(istart), 1, psiOld(istart), 1)
!!!!    write(*,'(a,i4,i3,i5,2es12.5)') 'pulayAt, iproc, iorb, <psi|psi>', pulayAt, iproc, iorb, ddot(Glr%wfd%nvctr_c+7*Glr%wfd%nvctr_f, psi(istart), 1, psi(istart), 1)
!!!!    !write(*,'(a,i4,i3,i5,2es12.5)') 'pulayAt, iproc, iorb, <psiOld|psiOld>', pulayAt, iproc, iorb, ddot(Glr%wfd%nvctr_c+7*Glr%wfd%nvctr_f, psiOld(istart), 1, psiOld(istart), 1)
!!!!    istart=istart+Glr%wfd%nvctr_c+7*Glr%wfd%nvctr_f
!!!!    write(*,'(a,2i5,2i12)') 'iproc, iorb, istart, size(psiOld)', iproc, iorb, istart, size(psiOld)
!!!!end do
!!!!if(iter<0) then
!!!!    write(*,*) iproc, 'reads', size(psi)
!!!!    ! waste time
!!!!    do istat=1,10000000
!!!!        tt=exp(dble(istat))
!!!!    end do
!!!!    do istat=1,size(psi)
!!!!        read(60+iproc,*) psiOld(istat)
!!!!    end do
!!!!    do iorb=1,orbs%norbp
!!!!        write(*,*) 'iproc, ddot2', ddot(Glr%wfd%nvctr_c+7*Glr%wfd%nvctr_f, psi(1), 1, psiOld(1), 1)/(dnrm2(Glr%wfd%nvctr_c+7*Glr%wfd%nvctr_f, psi(1), 1)*dnrm2(Glr%wfd%nvctr_c+7*Glr%wfd%nvctr_f, psiOld(1), 1))
!!!!    end do
!!!!end if
!!!
!!!
!!!!!! ONLY A TEST
!!!!!!call initRandomSeed(iproc, 1)
!!!!!!call random_number(phi)
!!!!!call transpose_v(iproc, nproc, orbsLIN, Glr%wfd, commsLIN, phi, work=phiWorkPointer)
!!!!!call orthogonalize(iproc, nproc, orbsLIN, commsLIN, Glr%wfd, phi, input)
!!!!!call untranspose_v(iproc, nproc, orbsLIN, Glr%wfd, commsLIN, phi, work=phiWorkPointer)
!!!deallocate(phiWorkPointer, stat=istat)
!!!
!!!
!!!!!istart=1
!!!!!do i=1,orbs%norbp
!!!!!    jstart=1
!!!!!    do j=1,i
!!!!!        write(*,*) 'i, j, ddot', i, j, ddot(Glr%wfd%nvctr_c+7*Glr%wfd%nvctr_f, psi(istart), 1, psi(jstart), 1)
!!!!!        !write(*,*) 'i, j, ddot', i, j, ddot(Glr%wfd%nvctr_c+7*Glr%wfd%nvctr_f, phi(istart), 1, phi(jstart), 1)
!!!!!        jstart=jstart+Glr%wfd%nvctr_c+7*Glr%wfd%nvctr_f
!!!!!    end do
!!!!!    istart=istart+Glr%wfd%nvctr_c+7*Glr%wfd%nvctr_f
!!!!!end do
!!!
!!!
!!!end subroutine improveOrbitals




subroutine createInputGuess(iproc, orbsLIN, Glr, input, at, rxyz, phi)
!
! Purpose:
! ========
!   Crates an input guess for the the phi orbitals based on the solutions
!   of the quantum mechanical harmonic oscillator.
! ATTENTION: Parabola is not defined with factor 1/2, so this input guess is maybe not good
!
! Input arguments:
! ================
!
use module_base
use module_types
implicit none

! Calling arguments
integer:: iproc
type(orbitals_data), intent(inout) :: orbsLIN
type(locreg_descriptors), intent(in) :: Glr
type(input_variables), intent(in) :: input
type(atoms_data), intent(in) :: at
real(8),dimension(3,at%nat):: rxyz
real(8),dimension(orbsLIN%npsidim):: phi

! Local variables
integer:: iorb, ideg, ndeg, n1, n2, n3, ninguess, istat, istart, ii, itot, iiAt
integer,dimension(:,:),allocatable:: quantumNumbers
real(8),dimension(:),allocatable:: phir
real(8):: kx, ky, kz, tt, dnrm2
type(workarr_locham) :: w2




allocate(phir(Glr%d%n1i*Glr%d%n2i*Glr%d%n3i),stat=istat)
call initialize_work_arrays_locham(Glr,orbsLIN%nspinor,w2)

kx=0.d0
ky=0.d0
kz=0.d0




! nexc means that we hav to calculate the nexc-excited state. Since there is degeneracy in the 3-dimensional case, 
! we have to find out which state belongs to nexc.

! We have to find orbsLIN%norbp input wavefunctions.
ninguess=0
istart=0
itot=0
iorb=0
!outLoop: do iorb=1,orbsLIN%norbp
outLoop: do
    iorb=iorb+1
    ! ndeg gives the degeneracy of this level
    ndeg=int(dble((iorb+1)*(iorb+2))/2.d0)
    allocate(quantumNumbers(3,ndeg), stat=istat)
    call getQuantumNumbers(iorb, ndeg, quantumNumbers)
    do ideg=1,ndeg
        ! It is possible that some orbitals for a given atoms are generated by another (the 'previous') MPI process.
        ! Therefore skip these levels until we reach a level that has not been used by another MPI process.
        itot=itot+1
        if(itot<=orbsLIN%nbasisOnPreviousMPI) cycle

        ninguess=ninguess+1
        ii=orbsLIN%onWhichAtom(ninguess)
        iiAt=at%iatype(ii)
        call getEigenfunction(orbsLIN, Glr, input, quantumNumbers(1,ideg), rxyz(1,ii), orbsLIN%parabPrefacArr(iiAt), phir)
!do istat=1,Glr%d%n1i*Glr%d%n2i*Glr%d%n3i
!    if(iproc==0) write(900+ninguess,*) istat, phir(istat)
!end do
        tt=input%hx*.5d0*dnrm2(Glr%d%n1i*Glr%d%n2i*Glr%d%n3i, phir(1), 1)
 write(*,*) 'iproc, iorb, tt', iproc, iorb, tt
        !call dscal(Glr%d%n1i*Glr%d%n2i*Glr%d%n3i, 1/tt, phir(1), 1)
        call isf_to_daub(input%hx, input%hy, input%hz, kx, ky, kz, orbsLIN%nspinor, Glr, w2, phir(1), phi(istart+1), tt)
!do istat=1,(Glr%wfd%nvctr_c+7*Glr%wfd%nvctr_f)*orbsLIN%nspinor
!    if(iproc==0) write(1100+ninguess,*) istat, phi(istart+istat)
!end do
        istart=istart+(Glr%wfd%nvctr_c+7*Glr%wfd%nvctr_f)*orbsLIN%nspinor
        if(ninguess==orbsLIN%norbp) exit outLoop
    end do
    deallocate(quantumNumbers, stat=istat)
end do outLoop





end subroutine createInputGuess




subroutine getQuantumNumbers(maxNum, ndeg, quantumNumbers)
!
! Purpose:
! ========
!   Calculates the possible quantum numbers n1, n2, n3 for the 3-dimensional harmonic 
!   oscillator for a given energy level. This energy level is characterized by its degeneracy.
!   The level of degeneracy is (maxNum+1)*(maxNum+2)/2. The possible quantum numbers n1, n2, n3
!   have to fulfill n1+n2+n3=MaxNum.
!
! Calling arguments:
! ==================
!   Input arguments: 
!     maxNum           each quantum number n1, n2, n3 cannot be larger than maxNum
!     ndeg             the level of degeneracy for this maxNum
!   Output arguments:
!     quantumNumbers   array containing the allowed quantum numbers n1, n2, n3
!
implicit none

! Calling arguments
integer,intent(in):: maxNum, ndeg
integer,dimension(3,ndeg),intent(out):: quantumNumbers

! Local variables
integer:: i1, i2, i3, ii


ii=0
do i3=0,maxNum
    do i2=0,maxNum
        do i1=0,maxNum
            if(i1+i2+i3==maxNum) then
                ii=ii+1
                quantumNumbers(1,ii)=i1
                quantumNumbers(2,ii)=i2
                quantumNumbers(3,ii)=i3
            end if
        end do
    end do
end do

! Check whether the correct numbers of quantum numbers was created.
if(ii/=ndeg) then
    write(*,'(a)') 'ERROR: the degeneracy is ', ndeg, 'but ', ii, ' levels were created!'
end if

end subroutine getQuantumNumbers




subroutine getEigenfunction(orbsLIN, Glr, input, quantumNumbers, rxyz, parabPrefac, phir)
!
! Purpose:
! ========
!   Creates the eigenfunction for the 3-dimensional harmonic oscillator for given quantum
!   numbers n1, n2, n3. Since we have m=1 and hbar=1, this eigenfunction is
!   psi_n(x) = (omega/pi)^(1/4)*1/sqrt(2^n*n!)*H_n(sqrt(omega)*x)*exp(-.5*omega*x^2)
!   Here omega=sqrt(k) (the parabolic potential is defined as 1/2*k*x^2) and H_n is the
!   Hermite polynomial of N-th order.
! Calling arguments:
! ==================
!   Input arguments:
!   Output arguments:
!   
use module_base
use module_types
implicit none

! Calling arguments
type(orbitals_data), intent(inout) :: orbsLIN
type(locreg_descriptors), intent(in) :: Glr
type(input_variables), intent(in) :: input
integer,dimension(3):: quantumNumbers
real(8),dimension(3):: rxyz
real(8):: parabPrefac
real(8),dimension(Glr%d%n1i*Glr%d%n2i*Glr%d%n3i),intent(out):: phir

! Local variables
integer:: ix, iy, iz, factorial, ii, N, ix0, iy0, iz0
real(8):: fac, chi, Hermite, phirx, phiry, phirz, tt, omega
real(8),parameter:: pi=3.141592653589793d0


!write(*,'(a,3i4)') 'quantumNumbers:', quantumNumbers(1), quantumNumbers(2), quantumNumbers(3)

omega=sqrt(2.d0*parabPrefac*5.d0)
!omega=sqrt(2.d0*parabPrefac*50.d0)

ix0=nint(rxyz(1)/(input%hx*.5d0))+15
iy0=nint(rxyz(2)/(input%hy*.5d0))+15
iz0=nint(rxyz(3)/(input%hz*.5d0))+15


!write(*,'(a,3es14.6)') 'rxyz(1), rxyz(2), rxyz(3)', rxyz(1), rxyz(2), rxyz(3)
ii=1
do iz=1,Glr%d%n3i
    do iy=1,Glr%d%n2i
        do ix=1,Glr%d%n1i

            !chi=orbsLIN%parabPrefac**.25d0*(ix-rxyz(1))*input%hx*.5d0
            !chi=sqrt(omega)*(ix*input%hx*.5d0-rxyz(1))
            chi=sqrt(omega)*(ix-ix0)*input%hx*.5d0
            N=quantumNumbers(1)
!if(N>1) write(999,*) 'N, factorial(N)', N, factorial(N)
            fac=(omega/pi)**(.25d0)/sqrt(dble(factorial(N))*2.d0**N)
            phirx=fac*exp(-.5d0*chi**2)*Hermite(N, chi)

            !chi=orbsLIN%parabPrefac**.25d0*(iy-rxyz(2))*input%hy*.5d0
            !chi=sqrt(omega)*(iy*input%hy*.5d0-rxyz(2))
            chi=sqrt(omega)*(iy-iy0)*input%hx*.5d0
            N=quantumNumbers(2)
            fac=(omega/pi)**(.25d0)/sqrt(dble(factorial(N))*2.d0**N)
            phiry=fac*exp(-.5d0*chi**2)*Hermite(N, chi)

            !chi=orbsLIN%parabPrefac**.25d0*(iz-rxyz(3))*input%hz*.5d0
            !chi=sqrt(omega)*(iz*input%hz*.5d0-rxyz(3))
            chi=sqrt(omega)*(iz-iz0)*input%hx*.5d0
            N=quantumNumbers(3)
            fac=(omega/pi)**(.25d0)/sqrt(dble(factorial(N))*2.d0**N)
            phirz=fac*exp(-.5d0*chi**2)*Hermite(N, chi)

            !!! Add some random noise
            !!call random_number(tt)
            !!tt=tt*.1d0
            !!tt=tt+.95d0
            tt=1.d0
            phir(ii)=phirx*phiry*phirz*tt
!if(phir(ii)>5.d-1) write(*,*) 'phir(ii)',phir(ii)
            ii=ii+1
        end do
    end do
end do



end subroutine getEigenfunction




function factorial(N)
!
! Purpose:
! ========
!   Calculates the factorial of N.
! Calling arguments:
! ==================
!   Input arguments:
!     N           value for which the factorial shall be calculated
!   Output arguments:
!     factorial   the factorial of N
!
implicit none
integer:: N
integer:: factorial

integer:: i

factorial=1
do i=1,N
    factorial=factorial*i
end do

end function factorial




function Hermite(N, x)
!
! Purpose:
! ========
!   Calculate the value of the Hermite polynomial of N-th order at the
!   point x. A recursion formula is used, namely
!   H_n(x) = 2*x*H_{n-1}(x)-2*(n-1)*H_{n-2}(x)  with H_0(x) = 1, H_1(x) = 2*x
!
! Calling arguments:
! ==================
!   Input arguments:
!     N         order of the polynomial
!     x         argument of the polynomial
!   Output arguments:
!     Hermite   the value of the polynomial at position x
!
implicit none

! Calling arguments
integer:: N
real(8):: x
real(8):: Hermite

! Local variables
integer:: i, istat
real(8),dimension(:),allocatable:: HerArr

if(N==0) then
    Hermite=1.d0
else if(N==1) then
    Hermite=2.d0*x
else
    allocate(HerArr(0:N), stat=istat)
    HerArr(0)=1.d0
    HerArr(1)=2.d0*x
    do i=2,N
        HerArr(i)=2.d0*x*HerArr(i-1)-2.d0*dble(i-1)*HerArr(i-2)
    end do
    Hermite=HerArr(N)
    deallocate(HerArr, stat=istat)
end if

end function



subroutine pulayNew(iproc, nproc, at, orbs, lr, input, orbsLIN, commsLIN, rxyz, nspin, &
    nlpspd, proj, nscatterarr, ngatherarr, rhopot, GPU, pkernelseq, phi, hphi, &
    rxyzParabola, pulayAt, pulayDir, shift)
!
! Purpose:
! ========
!   Calculates the localized basis functions phi. These basis functions are eigenfunctions of the ordinary Hamiltonian
!   with an additional parabolic potential centered at the atoms. The eigenfunctions are determined by minimizing the trace.
!
! Calling arguments:
!   Input arguments
!   Output arguments
!    phi   the localized basis functions
!
use module_base
use module_types
use module_interfaces, except_this_one => pulayNew
  use Poisson_Solver
!use allocModule
implicit none

! Calling arguments
integer:: iproc, nproc
type(atoms_data), intent(in) :: at
type(orbitals_data):: orbs
type(locreg_descriptors), intent(in) :: lr
type(input_variables):: input
type(orbitals_data):: orbsLIN
type(communications_arrays):: commsLIN
real(8),dimension(3,at%nat):: rxyz, rxyzParabola
integer:: nspin
type(nonlocal_psp_descriptors), intent(in) :: nlpspd
real(wp), dimension(nlpspd%nprojel), intent(in) :: proj
integer, dimension(0:nproc-1,4), intent(in) :: nscatterarr !n3d,n3p,i3s+i3xcsh-1,i3xcsh
integer, dimension(0:nproc-1,2), intent(in) :: ngatherarr 
real(dp), dimension(*), intent(inout) :: rhopot
type(GPU_pointers), intent(inout) :: GPU
real(dp), dimension(:), pointer :: pkernelseq
real(8),dimension(orbsLIN%npsidim):: phi, hphi
integer:: pulayAt, pulayDir
real(8):: shift
! Local variables
integer:: iorb, jorb, korb, lorb, istat, istart, nvctrp, ierr, iat, ii, jj, kk, ik, il, ll, jproc, kproc, lproc, id, norbOnAt, info
integer:: jstart, lwork, ncplx, it, iiAt, idir
real(8),dimension(:),allocatable:: b, a, dvec, eval, work, phiw, hphiw, phi1, phiw2
real(8),dimension(:,:),allocatable:: eps, dTemp, d
real(8),dimension(:,:),allocatable:: emat, U, HU, HtildeSmall
real(8),dimension(:,:,:),allocatable:: Htilde
real(8):: tt, hphiNrm, dnrm2, fracTot, parabPrefac, ddot, ekin_sum, epot_sum, eexctx, eproj_sum
integer,dimension(:,:),allocatable:: orbsOnAtom
integer,dimension(:,:,:),allocatable:: tempArr
integer,dimension(:),allocatable:: lorbArr, displs, ipiv
real(8),dimension(:),pointer:: phiWorkPointer
character(len=1):: direction
!real(8),dimension(:,:),allocatable:: rxyzParabola


! Calculate the unconstrained gradient.
call HamiltonianApplicationParabola(iproc,nproc,at,orbsLIN,input%hx,input%hy,input%hz,rxyz,&
     nlpspd,proj,lr,ngatherarr,lr%d%n1i*lr%d%n2i*nscatterarr(iproc,2),&
     rhopot(1),&
     phi(1),hphi(1),ekin_sum,epot_sum,eexctX,eproj_sum,nspin,GPU, orbsLIN%onWhichAtom, rxyzParabola, pkernel=pkernelseq)

allocate(orbsOnAtom(orbsLIN%norb,2), stat=istat)
allocate(lorbArr(0:nproc-1), stat=istat)
allocate(displs(0:nproc-1), stat=istat)
atomsLoop: do iat=pulayAt,pulayAt
    orbsOnAtom=0
    nvctrp=sum(commsLIN%nvctr_par(iproc,1:orbsLIN%nkptsp))*orbsLIN%nspinor


    ! First find out which atoms are centered on atom iat
    korb=0
    lorb=0
    do jproc=0,nproc-1
        do jorb=1,orbsLIN%norb_par(jproc)
            korb=korb+1
            if(jproc==iproc) then
                if(orbsLIN%onWhichAtom(jorb)==iat) then
                    lorb=lorb+1
                    orbsOnAtom(lorb,2)=korb
                end if
            end if
        end do
    end do
    call mpi_gather(lorb, 1, mpi_integer, lorbArr(0), 1, mpi_integer, 0, mpi_comm_world, ierr)
    displs(0)=0
    do jproc=1,nproc-1
        displs(jproc)=displs(jproc-1)+lorbArr(jproc-1)
    end do
    call mpi_gatherv(orbsOnAtom(1,2), lorb, mpi_integer, orbsOnAtom(1,1), lorbArr(0), displs(0), &
        mpi_integer, 0, mpi_comm_world, ierr)
    do istat=1,sum(lorbArr)
        !if(iproc==0) write(*,*)'istat, orbsOnAtom(istat,1)', istat, orbsOnAtom(istat,1)
    end do

    ! Send orbsOnAtom to all processes
    if(iproc==0) norbOnAt=sum(lorbArr)
    call mpi_bcast(norbOnAt, 1, mpi_integer, 0, mpi_comm_world, ierr)
    call mpi_barrier(mpi_comm_world, ierr)
    call mpi_bcast(orbsOnAtom(1,1), norbOnat, mpi_integer, 0, mpi_comm_world, ierr)
    call mpi_barrier(mpi_comm_world, ierr)
 
    ! Calculate <phi|H|phi> for all orbitals. This requires to tranpose them first.
    allocate(Htilde(orbsLIN%norb,orbsLIN%norb,2), stat=istat)
    Htilde=0.d0
    allocate(phiWorkPointer(size(phi)), stat=istat)
    call transpose_v(iproc, nproc, orbsLIN, lr%wfd, commsLIN, phi, work=phiWorkPointer)
    call transpose_v(iproc, nproc, orbsLIN, lr%wfd, commsLIN, hphi, work=phiWorkPointer)
    deallocate(phiWorkPointer, stat=istat)
    nvctrp=commsLIN%nvctr_par(iproc,1)
    istart=1
    do iorb=1,orbsLIN%norb
        jstart=1
        do jorb=1,orbsLIN%norb
            Htilde(iorb,jorb,2)=ddot(nvctrp, phi(istart), 1, hphi(jstart), 1)
            jstart=jstart+nvctrp
        end do
        istart=istart+nvctrp
    end do
    call mpi_allreduce (Htilde(1,1,2), Htilde(1,1,1), orbsLIN%norb**2, mpi_double_precision, mpi_sum, mpi_comm_world, ierr)
    if(iproc==0) write(*,*) '<phi|H|phi>'
    do iorb=1,orbsLIN%norb
        if(iproc==0) write(*,'(100f8.4)') (Htilde(iorb,jorb,1), jorb=1,orbsLIN%norb)
    end do

    ! Cut out the part of the matrix containing all orbitals centered on the current atom
    allocate(HtildeSmall(norbOnAt,norbOnAt), stat=istat)
    if(iproc==0) then
        do iorb=1,norbOnAt
            do jorb=1,norbOnAt
                HtildeSmall(jorb,iorb)=Htilde(orbsOnAtom(jorb,1),orbsOnAtom(iorb,1),1)
            end do
        end do 
        !do iorb=1,norbOnAt
        !    write(*,'(100f8.4)') (HtildeSmall(iorb,jorb), jorb=1,norbOnAt)
        !end do
        
        ! Diagonalize Htilde on root only
        allocate(eval(norbOnAt), stat=istat)
        lwork=1000*norbOnAt
        allocate(work(lwork), stat=istat)
        call dsyev('v', 'l', norbOnAt, HtildeSmall, norbOnAt, eval, work, lwork, info)
        if(info/=0) then
            write(*,'(a,i0)') 'ERROR in dsyev, info= ',info
            stop
        end if
        deallocate(work, stat=istat)
        !write(*,*) '--------------'
        !do iorb=1,norbOnAt
        !    write(*,'(100f8.4)') (HtildeSmall(iorb,jorb), jorb=1,norbOnAt)
        !end do
    end if

    ! Broadcast HtildeSmall
    call mpi_bcast(HtildeSmall(1,1), norbOnAt**2, mpi_double_precision, 0, mpi_comm_world, ierr)
    call mpi_barrier(mpi_comm_world, ierr)
    !write(*,*) '--------------', iproc, norbOnAt
    !do iorb=1,norbOnAt
    !    write(*,'(100f8.4)') (HtildeSmall(iorb,jorb), jorb=1,norbOnAt)
    !end do

    ! Build new linear combination
    allocate(phiw(orbsLIN%npsidim), stat=istat)
    allocate(hphiw(orbsLIN%npsidim), stat=istat)
    phiw=0.d0 ; hphiw=0.d0
    nvctrp=commsLIN%nvctr_par(iproc,1)
    do iorb=1,norbOnAt
        istart=(orbsOnAtom(iorb,1)-1)*nvctrp+1
        do jorb=1,norbOnAt
            jstart=(orbsOnAtom(jorb,1)-1)*nvctrp+1
            call daxpy(nvctrp, HtildeSmall(jorb,iorb), phi(jstart), 1, phiw(istart), 1)
        end do
    end do

    !!! Undo this transformations only a test)
    !!allocate(phiw2(size(phiw)), stat=istat)
    !!phiw2=phiw
    !!phiw=0.d0
    !!do iorb=1,norbOnAt
    !!    istart=(orbsOnAtom(iorb,1)-1)*nvctrp+1
    !!    do jorb=1,norbOnAt
    !!        jstart=(orbsOnAtom(jorb,1)-1)*nvctrp+1
    !!        call daxpy(nvctrp, HtildeSmall(iorb,jorb), phiw2(jstart), 1, phiw(istart), 1)
    !!    end do
    !!end do

    allocate(phiWorkPointer(size(phi)), stat=istat)
    call untranspose_v(iproc, nproc, orbsLIN, lr%wfd, commsLIN, phi, work=phiWorkPointer)
    call untranspose_v(iproc, nproc, orbsLIN, lr%wfd, commsLIN, phiw, work=phiWorkPointer)
    call untranspose_v(iproc, nproc, orbsLIN, lr%wfd, commsLIN, hphi, work=phiWorkPointer)
    deallocate(phiWorkPointer, stat=istat)

    ! Check whether the subspace diagonalization was successful
    ! This should at the moment always be done to get the diagonal elements
    if(.true.) then
        nspin=1
        !call HamiltonianApplicationParabola(iproc,nproc,at,orbsLIN,input%hx,input%hy,input%hz,rxyz,&
        !     nlpspd,proj,lr,ngatherarr,lr%d%n1i*lr%d%n2i*nscatterarr(iproc,2),&
        !     rhopot(1),&
        !     phiw(1),hphiw(1),ekin_sum,epot_sum,eexctX,eproj_sum,nspin,GPU, orbsLIN%onWhichAtom, rxyzParabola, pkernel=pkernelseq)
        call HamiltonianApplicationParabola(iproc,nproc,at,orbsLIN,input%hx,input%hy,input%hz,rxyz,&
            nlpspd,proj,lr,ngatherarr,lr%d%n1i*lr%d%n2i*nscatterarr(iproc,2),&
            rhopot(1),&
            phiw(1),hphiw(1),ekin_sum,epot_sum,eexctX,eproj_sum,nspin,GPU, orbsLIN%onWhichAtom, rxyzParabola, pkernel=pkernelseq)

        Htilde=0.d0
        allocate(phiWorkPointer(size(phi)), stat=istat)
        call transpose_v(iproc, nproc, orbsLIN, lr%wfd, commsLIN, phiw, work=phiWorkPointer)
        call transpose_v(iproc, nproc, orbsLIN, lr%wfd, commsLIN, hphiw, work=phiWorkPointer)
        deallocate(phiWorkPointer, stat=istat)
        nvctrp=commsLIN%nvctr_par(iproc,1)
        istart=1
        do iorb=1,orbsLIN%norb
            jstart=1
            do jorb=1,orbsLIN%norb
                Htilde(iorb,jorb,2)=ddot(nvctrp, phiw(istart), 1, hphiw(jstart), 1)
                jstart=jstart+nvctrp
            end do
            istart=istart+nvctrp
        end do
        call mpi_allreduce (Htilde(1,1,2), Htilde(1,1,1), orbsLIN%norb**2 ,mpi_double_precision, mpi_sum, mpi_comm_world, ierr)
        !if(iproc==0) write(*,*) 'after subspace diag: <phi|H|phi>'
        !do iorb=1,orbsLIN%norb
        !    if(iproc==0) write(*,'(100f8.4)') (Htilde(iorb,jorb,1), jorb=1,orbsLIN%norb)
        !end do
    end if       

    directionLoop: do idir=pulayDir,pulayDir
           if(idir==1) then
               direction='x'
           else if(idir==2) then
               direction='y'
           else if(idir==3) then
               direction='z'
           end if
           ! Calculate the matrix elements <phiw|V|phiw>
           ncplx=1
           it=1
           ! First apply the perturbation to all orbitals belonging to iproc
           allocate(phiWorkPointer(size(phi)), stat=istat)
           call untranspose_v(iproc, nproc, orbsLIN, lr%wfd, commsLIN, phiw, work=phiWorkPointer)
           hphiw=phiw
           istart=1
           do iorb=1,orbsLIN%norbp
               iiAt=orbsLIN%onWhichAtom(iorb)
               parabPrefac=orbsLIN%parabPrefacArr(at%iatype(orbsLIN%onWhichAtom(iorb)))
               !write(*,'(a,2i6,es12.4)') 'before: iproc, iorb, dnrm2(hphiw)', iproc, iorb, dnrm2(nvctrp,hphiw(istart),1)
               call subTemp(lr,ncplx,&
                    input%hx,input%hy,input%hz,hphiw(istart), rxyzParabola(3,iat), orbsLIN, parabPrefac, it, direction, shift)
               !write(*,'(a,2i6,es12.4)') 'after: iproc, iorb, dnrm2(hphiw)', iproc, iorb, dnrm2(nvctrp,hphiw(istart),1)
               istart=istart+(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f)*orbsLIN%nspinor
           end do
    
           call transpose_v(iproc, nproc, orbsLIN, lr%wfd, commsLIN, phiw, work=phiWorkPointer)
           call transpose_v(iproc, nproc, orbsLIN, lr%wfd, commsLIN, hphiw, work=phiWorkPointer)
    
    
           ! nvctrp is the amount of each phi hold by the current process
           nvctrp=sum(commsLIN%nvctr_par(iproc,1:orbsLIN%nkptsp))*orbsLIN%nspinor
           
           allocate(dTemp(orbsLIN%norb,orbsLIN%norb), stat=istat)
           allocate(d(orbsLIN%norb,orbsLIN%norb), stat=istat)
           dTemp=0.d0
           istart=1
           do iorb=1,orbsLIN%norb
               jstart=1
               do jorb=1,orbsLIN%norb
                   dTemp(iorb,jorb)=ddot(nvctrp, phiw(istart), 1, hphiw(jstart), 1)
                   !if(iproc==0) write(*,'(a,2i6,2es12.4)') 'iorb, jorb, dnrm2(phiw), dnrm2(hphiw)', iorb, jorb, dnrm2(nvctrp,phiw(istart),1), dnrm2(nvctrp,hphiw(jstart),1)
                   jstart=jstart+nvctrp
               end do
               istart=istart+nvctrp
           end do
           d=0.d0
           call mpi_allreduce(dTemp(1,1), d(1,1), orbsLIN%norb**2, mpi_double_precision, mpi_sum, mpi_comm_world, ierr)
           ! Add the diagonal element shift**2
           do iorb=1,orbsLIN%norb
               d(iorb,iorb)=d(iorb,iorb)+parabPrefac*dble(shift)**2
           end do
           
           if(iproc==0) write(*,'(a,es12.4)') '<phi|V|phi> with shift=',shift
           do iorb=1,orbsLIN%norb
               if(iproc==0) write(*,'(100es15.8)') (d(iorb,jorb), jorb=1,orbsLIN%norb)
           end do
    
           ! Make a linear combination according to perturbation theory
           !!do iorb=1,orbsLIN%norbp
           !!    write(*,'(a,i5,i3,i5,es12.5)') 'iat, iproc, iorb, dnrm2(phiw)', iat, iproc, iorb, dnrm2(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f, phiw(istart), 1)
           !!    istart=istart+lr%wfd%nvctr_c+7*lr%wfd%nvctr_f
           !!end do
           allocate(phi1(orbsLIN%npsidim), stat=istat)
           phi1=0.d0
           do iorb=1,norbOnAt
               !istart=(iorb-1)*nvctrp+1
               istart=(orbsOnAtom(iorb,1)-1)*nvctrp+1
               do jorb=1,norbOnAt
                   jstart=(orbsOnAtom(jorb,1)-1)*nvctrp+1
                   if(iorb==jorb) cycle
                   call daxpy(nvctrp, d(orbsOnATom(iorb,1),orbsOnAtom(jorb,1))/ &
                   (Htilde(orbsOnAtom(iorb,1),orbsOnAtom(iorb,1),1)-Htilde(orbsOnAtom(jorb,1),orbsOnAtom(jorb,1),1)), &
                   phiw(jstart), 1, phi1(istart), 1)
               end do
           end do
           call untranspose_v(iproc, nproc, orbsLIN, lr%wfd, commsLIN, phi1, work=phiWorkPointer)
           !istart=1
           !do iorb=1,orbsLIN%norbp
           !    write(*,'(a,i3,i5,2es12.5)') 'iproc, iorb, dnrm2(phi1), dnrm2(phi)', iproc, iorb, dnrm2(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f, phi1(istart), 1), dnrm2(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f, phi(istart), 1)
           !    istart=istart+lr%wfd%nvctr_c+7*lr%wfd%nvctr_f
           !end do
           call transpose_v(iproc, nproc, orbsLIN, lr%wfd, commsLIN, phi1, work=phiWorkPointer)

           ! Add the correction phi1 to the old orbitals phiw
           istart=1
           do iorb=1,norbOnAt
               jstart=(orbsOnAtom(iorb,1)-1)*nvctrp+1
               call daxpy(nvctrp, 1.d0, phi1(jstart), 1, phiw(jstart), 1)
               !call daxpy(nvctrp, 1.d0, phi1(istart), 1, phiw(jstart), 1)
               istart=istart+nvctrp
           end do

           call transpose_v(iproc, nproc, orbsLIN, lr%wfd, commsLIN, phi, work=phiWorkPointer)
           ! Transform back to 'old' orbitals
           allocate(phiw2(size(phiw)), stat=istat)
           phiw2=phiw
           phiw=0.d0
           do iorb=1,norbOnAt
               istart=(orbsOnAtom(iorb,1)-1)*nvctrp+1
               do jorb=1,norbOnAt
                   jstart=(orbsOnAtom(jorb,1)-1)*nvctrp+1
                   call daxpy(nvctrp, HtildeSmall(iorb,jorb), phiw2(jstart), 1, phiw(istart), 1)
               end do
           end do

           ! Now mix the orbitals
           phiw2=phi
           istart=1
           do iorb=1,norbOnAt
               jstart=(orbsOnAtom(iorb,1)-1)*nvctrp+1
               !call dcopy(nvctrp, phiw(istart), 1, phiw2(jstart), 1)
               call dcopy(nvctrp, phiw(jstart), 1, phiw2(jstart), 1)
               istart=istart+nvctrp
           end do

           ! This is not necessary, only for debugging
           call untranspose_v(iproc, nproc, orbsLIN, lr%wfd, commsLIN, phiw2, work=phiWorkPointer)
           !istart=1
           !do iorb=1,orbsLIN%norbp
           !    write(*,'(a,i4,i3,i5,2es12.5)') 'iat, iproc, iorb, dnrm2(phiw2)', iat, iproc, iorb, dnrm2(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f, phiw2(istart), 1)
           !    istart=istart+lr%wfd%nvctr_c+7*lr%wfd%nvctr_f
           !end do
           call transpose_v(iproc, nproc, orbsLIN, lr%wfd, commsLIN, phiw2, work=phiWorkPointer)


           ! Now orthonormalize
           call orthogonalize(iproc, nproc, orbsLIN, commsLIN, lr%wfd, phiw2, input)
    
           call untranspose_v(iproc, nproc, orbsLIN, lr%wfd, commsLIN, phi, work=phiWorkPointer)
           call untranspose_v(iproc, nproc, orbsLIN, lr%wfd, commsLIN, phiw, work=phiWorkPointer)
           call untranspose_v(iproc, nproc, orbsLIN, lr%wfd, commsLIN, hphiw, work=phiWorkPointer)
           call untranspose_v(iproc, nproc, orbsLIN, lr%wfd, commsLIN, phiw2, work=phiWorkPointer)
           deallocate(phiWorkPointer, stat=istat)

           ! Now phiw2 contains the perturbed basis functions.
           ! Copy them back
           !istart=1
           !do iorb=1,orbsLIN%norbp
           !    write(*,'(a,i4,i3,i5,2es24.15)') 'iat, iproc, iorb, <phi|phiw2>', iat, iproc, iorb, ddot(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f, phi(istart), 1, phiw2(istart), 1)
           !    istart=istart+lr%wfd%nvctr_c+7*lr%wfd%nvctr_f
           !end do
           phi=phiw2
    
        !call mpi_barrier(mpi_comm_world, ierr)
        !stop
    end do directionLoop
end do atomsLoop

end subroutine pulayNew


subroutine pulay(iproc, nproc, lr, orbs, orbsLIN, comms, commsLIN, input, at, rxyz, phi, hphi, &
    psi, nscatterarr, ngatherarr, nlpspd, proj, sizeRhopot, rhopot, GPU, pkernelseq)
!
! Purpose:
! ========
!   Calculates the Pulay forces
!
! Calling arguments:
! ==================
!
use module_base
use module_types
use module_interfaces, except_this_one => pulay
implicit none

! Calling arguments
integer:: iproc, nproc, sizeRhopot
type(locreg_descriptors), intent(in) :: lr
type(orbitals_data), intent(inout) :: orbs, orbsLIN
type(communications_arrays), intent(in) :: comms
type(communications_arrays), intent(in) :: commsLIN
type(atoms_data), intent(in) :: at
real(8),dimension(3,at%nat):: rxyz
integer, dimension(0:nproc-1,4), intent(in) :: nscatterarr !n3d,n3p,i3s+i3xcsh-1,i3xcsh
integer, dimension(0:nproc-1,2), intent(in) :: ngatherarr 
type(nonlocal_psp_descriptors), intent(in) :: nlpspd
real(wp), dimension(nlpspd%nprojel), intent(in) :: proj
!real(dp), dimension(*), intent(inout) :: rhopot
real(dp), dimension(sizeRhopot), intent(inout) :: rhopot
type(GPU_pointers), intent(inout) :: GPU
type(input_variables):: input
real(dp), dimension(:), pointer :: pkernelseq
real(8),dimension(orbsLIN%npsidim):: phi, hphi
real(8),dimension(orbs%npsidim):: psi
!integer:: iter
!!
! Local variables
integer:: iorb, jorb, korb, lorb, istat, istart, nvctrp, ierr, iat, ii, jj, kk, ik, il, ll, jproc, kproc, lproc, id, norbOnAt, info
integer:: jstart, lwork, nspin
real(8),dimension(:),allocatable:: b, d, a, dvec, eval, work, phiw, hphiw
real(8),dimension(:,:),allocatable:: eps
real(8),dimension(:,:),allocatable:: emat, U, HU, HtildeSmall
real(8),dimension(:,:,:),allocatable:: Htilde
real(8):: tt, hphiNrm, dnrm2, fracTot, parabPrefac, ddot, ekin_sum, epot_sum, eexctx, eproj_sum
integer,dimension(:,:),allocatable:: orbsOnAtom
integer,dimension(:,:,:),allocatable:: tempArr
integer,dimension(:),allocatable:: lorbArr, displs, ipiv
real(8),dimension(:),pointer:: phiWorkPointer
real(8),dimension(:,:),allocatable:: rxyzParabola
!!
 write(*,*) 'ATTENTION -- TEST!'
 write(*,'(a,2i9)') 'size(rhopot), lr%d%n1i*lr%d%n2i*nscatterarr(iproc,2)', size(rhopot), lr%d%n1i*lr%d%n2i*nscatterarr(iproc,2)
       allocate(rxyzParabola(3,at%nat), stat=istat)
       rxyzParabola=rxyz
       call HamiltonianApplicationParabola(iproc,nproc,at,orbsLIN,input%hx,input%hy,input%hz,rxyz,&
            nlpspd,proj,lr,ngatherarr,lr%d%n1i*lr%d%n2i*nscatterarr(iproc,2),&
            rhopot(1),&
            phi(1),hphi(1),ekin_sum,epot_sum,eexctX,eproj_sum,nspin,GPU, orbsLIN%onWhichAtom, rxyzParabola, pkernel=pkernelseq)

fracTot=0.d0
allocate(eps(orbsLIN%norb,orbsLIN%norb), stat=istat)
allocate(d(orbsLIN%norb**2), stat=istat)
allocate(orbsOnAtom(orbsLIN%norb,2), stat=istat)
allocate(lorbArr(0:nproc-1), stat=istat)
allocate(displs(0:nproc-1), stat=istat)
orbsOnAtom=0
nvctrp=sum(commsLIN%nvctr_par(iproc,1:orbsLIN%nkptsp))*orbsLIN%nspinor
!!psiLoop: do korb=1,orbs%norb
    istart=1
    ! First calculate the epsilons
    !call getEpsilon(iproc, nproc, orbsLIN, commsLIN, eps, phi(1), hphi(1))
    !do iorb=1,orbslIN%norb
    !    !if(iproc==0) write(*,'(100f9.5)') (eps((jorb-1)*orbsLIN%norb+iorb), jorb=1,orbsLIN%norb)
    !    if(iproc==0) write(*,'(100f9.5)') (eps(iorb,jorb), jorb=1,orbsLIN%norb)
    !end do
    !phi1Loop: do iorb=1,orbsLIN%norb
    atomsLoop: do iat=1,at%nat

        call mpi_barrier(mpi_comm_world, ierr)  

        !iiAt=orbsLIN%onWhichAtom(iorb)
        !parabPrefac=orbsLIN%parabPrefacArr(at%iatype(orbsLIN%onWhichAtom(iorb)))
        parabPrefac=orbsLIN%parabPrefacArr(at%iatype(iat))
        !call getD(iproc,nproc,orbsLIN,lr,input%hx,input%hy,input%hz,phi, hphi, at%nat, rxyz(1,iat), parabPrefac, at, commsLIN, d)
        !if(iproc==0) write(*,*) '===================================================='
        !if(iproc==0) write(*,*) '===================================================='
        !do iorb=1,orbsLIN%norb
        !    if(iproc==0) write(*,'(100f11.5)') (d((jorb-1)*orbsLIN%norb+iorb), jorb=1,orbsLIN%norb)
        !end do

        ! First find out which atoms are centered on atom iat
        korb=0
        lorb=0
        do jproc=0,nproc-1
            do jorb=1,orbsLIN%norb_par(jproc)
                korb=korb+1
                if(jproc==iproc) then
                    if(orbsLIN%onWhichAtom(jorb)==iat) then
                        lorb=lorb+1
                        orbsOnAtom(lorb,2)=korb
                    end if
                end if
            end do
        end do
        call mpi_gather(lorb, 1, mpi_integer, lorbArr(0), 1, mpi_integer, 0, mpi_comm_world, ierr)
        !if(iproc==0) then
            displs(0)=0
            do jproc=1,nproc-1
                displs(jproc)=displs(jproc-1)+lorbArr(jproc-1)
            end do
        !end if
        if(iproc==0) write(*,'(a,100i4)') 'lorbArr', lorbArr
        if(iproc==0) write(*,'(a,100i4)') 'displs', displs
        call mpi_barrier(mpi_comm_world, ierr)
        write(*,*) 'before mpi_gatherv, iproc', iproc
        call mpi_gatherv(orbsOnAtom(1,2), lorb, mpi_integer, orbsOnAtom(1,1), lorbArr(0), displs(0), &
            mpi_integer, 0, mpi_comm_world, ierr)
        write(*,*) 'after mpi_gatherv, iproc', iproc
        do istat=1,sum(lorbArr)
            if(iproc==0) write(*,*)'istat, orbsOnAtom(istat,1)', istat, orbsOnAtom(istat,1)
        end do
        call mpi_barrier(mpi_comm_world, ierr)
 write(*,*) 'here, iproc', iproc
        call mpi_barrier(mpi_comm_world, ierr)
  write(*,*) 'after barrier, iproc', iproc

        ! Send orbsOnAtom to all processes
        if(iproc==0) norbOnAt=sum(lorbArr)
        call mpi_bcast(norbOnAt, 1, mpi_integer, 0, mpi_comm_world, ierr)
        call mpi_barrier(mpi_comm_world, ierr)
          write(*,*) 'norbOnAt, iproc', norbOnAt, iproc
        call mpi_bcast(orbsOnAtom(1,1), norbOnat, mpi_integer, 0, mpi_comm_world, ierr)
        call mpi_barrier(mpi_comm_world, ierr)
   write(*,*) 'after bcast, iproc', iproc
 
        !if(iproc==0) then
            !norbOnAt=sum(lorbArr)
            !allocate(U(orbsLIN%npsidim,norbOnAt), stat=istat)
            !allocate(HU(orbsLIN%npsidim,norbOnAt), stat=istat)
            allocate(Htilde(orbsLIN%norb,orbsLIN%norb,2), stat=istat)
            Htilde=0.d0
            allocate(phiWorkPointer(size(phi)), stat=istat)
            call transpose_v(iproc, nproc, orbsLIN, lr%wfd, commsLIN, phi, work=phiWorkPointer)
            call transpose_v(iproc, nproc, orbsLIN, lr%wfd, commsLIN, hphi, work=phiWorkPointer)
            deallocate(phiWorkPointer, stat=istat)
            nvctrp=commsLIN%nvctr_par(iproc,1)
            istart=1
            do iorb=1,orbsLIN%norb
                jstart=1
                do jorb=1,orbsLIN%norb
 !write(*,*) iproc, iorb, jorb
                    Htilde(iorb,jorb,2)=ddot(nvctrp, phi(istart), 1, hphi(jstart), 1)
                    jstart=jstart+nvctrp
                end do
                istart=istart+nvctrp
            end do
            call mpi_allreduce (Htilde(1,1,2), Htilde(1,1,1), orbsLIN%norb**2, mpi_double_precision, mpi_sum, mpi_comm_world, ierr)
            do iorb=1,orbsLIN%norb
                if(iproc==0) write(*,'(100f8.4)') (Htilde(iorb,jorb,1), jorb=1,orbsLIN%norb)
            end do
        !end if

        allocate(HtildeSmall(norbOnAt,norbOnAt), stat=istat)
        if(iproc==0) then
            do iorb=1,norbOnAt
                do jorb=1,norbOnAt
                    HtildeSmall(jorb,iorb)=Htilde(orbsOnAtom(jorb,1),orbsOnAtom(iorb,1),1)
                end do
            end do 
            do iorb=1,norbOnAt
                write(*,'(100f8.4)') (HtildeSmall(iorb,jorb), jorb=1,norbOnAt)
            end do
            
            ! Diagonalize Htilde
             allocate(eval(norbOnAt), stat=istat)
             lwork=1000*norbOnAt
             allocate(work(lwork), stat=istat)
            call dsyev('v', 'l', norbOnAt, HtildeSmall, norbOnAt, eval, work, lwork, info)
            deallocate(work, stat=istat)
            write(*,*) '--------------'
            do iorb=1,norbOnAt
                write(*,'(100f8.4)') (HtildeSmall(iorb,jorb), jorb=1,norbOnAt)
            end do
        end if
        ! Broadcast HtildeSmall
        call mpi_bcast(HtildeSmall(1,1), norbOnAt**2, mpi_double_precision, 0, mpi_comm_world, ierr)
        call mpi_barrier(mpi_comm_world, ierr)
            write(*,*) '--------------', iproc, norbOnAt
            do iorb=1,norbOnAt
                write(*,'(100f8.4)') (HtildeSmall(iorb,jorb), jorb=1,norbOnAt)
            end do

       ! Build new linear combination
       allocate(phiw(orbsLIN%npsidim), stat=istat)
       allocate(hphiw(orbsLIN%npsidim), stat=istat)
       nvctrp=commsLIN%nvctr_par(iproc,1)
       jstart=1
       do iorb=1,norbOnAt
           istart=(orbsOnAtom(iorb,1)-1)*nvctrp+1
           call dcopy(nvctrp, phi(istart), 1, phiw(jstart), 1)
           jstart=jstart+nvctrp
       end do

       allocate(phiWorkPointer(size(phi)), stat=istat)
       call untranspose_v(iproc, nproc, orbsLIN, lr%wfd, commsLIN, phi, work=phiWorkPointer)
       call untranspose_v(iproc, nproc, orbsLIN, lr%wfd, commsLIN, phiw, work=phiWorkPointer)
       call untranspose_v(iproc, nproc, orbsLIN, lr%wfd, commsLIN, hphi, work=phiWorkPointer)
       deallocate(phiWorkPointer, stat=istat)

       nspin=1
       !call HamiltonianApplicationParabola(iproc,nproc,at,orbsLIN,input%hx,input%hy,input%hz,rxyz,&
       !     nlpspd,proj,lr,ngatherarr,lr%d%n1i*lr%d%n2i*nscatterarr(iproc,2),&
       !     rhopot(1),&
       !     phiw(1),hphiw(1),ekin_sum,epot_sum,eexctX,eproj_sum,nspin,GPU, orbsLIN%onWhichAtom, rxyzParabola, pkernel=pkernelseq)
       call HamiltonianApplicationParabola(iproc,nproc,at,orbsLIN,input%hx,input%hy,input%hz,rxyz,&
            nlpspd,proj,lr,ngatherarr,lr%d%n1i*lr%d%n2i*nscatterarr(iproc,2),&
            rhopot(1),&
            phi(1),hphi(1),ekin_sum,epot_sum,eexctX,eproj_sum,nspin,GPU, orbsLIN%onWhichAtom, rxyzParabola, pkernel=pkernelseq)

       Htilde=0.d0
       allocate(phiWorkPointer(size(phi)), stat=istat)
       call transpose_v(iproc, nproc, orbsLIN, lr%wfd, commsLIN, phi, work=phiWorkPointer)
       call transpose_v(iproc, nproc, orbsLIN, lr%wfd, commsLIN, hphi, work=phiWorkPointer)
       deallocate(phiWorkPointer, stat=istat)
       nvctrp=commsLIN%nvctr_par(iproc,1)
       istart=1
       do iorb=1,orbsLIN%norb
           jstart=1
           do jorb=1,orbsLIN%norb
               Htilde(iorb,jorb,2)=ddot(nvctrp, phi(istart), 1, hphi(jstart), 1)
               jstart=jstart+nvctrp
           end do
           istart=istart+nvctrp
       end do
       call mpi_allreduce (Htilde(1,1,2), Htilde(1,1,1), orbsLIN%norb**2 ,mpi_double_precision, mpi_sum, mpi_comm_world, ierr)
       do iorb=1,orbsLIN%norb
           if(iproc==0) write(*,'(100f8.4)') (Htilde(iorb,jorb,1), jorb=1,orbsLIN%norb)
       end do

    call mpi_barrier(mpi_comm_world, ierr)
    stop


        if(iproc==0) then
            if(allocated(emat)) deallocate(emat)
            norbOnAt=sum(lorbArr)
            allocate(emat(norbOnAt**2,norbOnAt**2))
            emat=0.d0
            !emat=1.d-100
        end if
        
        if(iproc==0) allocate(tempArr(2,norbOnAt**2,norbOnAt**2))
        if(iproc==0) tempArr=0
        ! Build the matrix emat
        if(iproc==0) then
            do jorb=1,norbOnAt**2
                do iorb=1,norbOnAt**2
                    ! the 'blocks'
                    !write(*,'(a,3i7)') 'iorb, sum(lorbArr), mod(iorb-1,sum(lorbArr))+1',iorb, sum(lorbArr), mod(iorb-1,sum(lorbArr))+1
                    !if(abs(iorb-jorb)<norbOnAt) emat(iorb,jorb)=eps(orbsOnAtom(mod(iorb-1,norbOnAt)+1,1),orbsOnAtom(mod(jorb-1,norbOnAt)+1,1))
                    !if(abs(iorb-jorb)<norbOnAt) tempArr(1,iorb,jorb)=orbsOnAtom(mod(iorb-1,norbOnAt)+1,1) ; tempArr(2,iorb,jorb)=orbsOnAtom(mod(jorb-1,norbOnAt)+1,1)
                    if(ceiling(dble(iorb)/dble(norbOnAt))==ceiling(dble(jorb)/dble(norbOnAt))) then
                        emat(iorb,jorb)=eps(orbsOnAtom(mod(iorb-1,norbOnAt)+1,1), &
                        orbsOnAtom(mod(jorb-1,norbOnAt)+1,1))
                        !tempArr(1,iorb,jorb)=orbsOnAtom(mod(iorb-1,norbOnAt)+1,1) ; tempArr(2,iorb,jorb)=orbsOnAtom(mod(jorb-1,norbOnAt)+1,1)
                    end if
                    ! the 'pseudodiagonal' elements
                    if(ceiling(dble(jorb)/dble(norbOnAt))==mod(iorb-1,norbOnAt)+1) then
                        emat(iorb,jorb)=emat(iorb,jorb)+eps(orbsOnAtom(ceiling(dble(iorb)/dble(norbOnAt)),1),&
                        orbsOnAtom(mod(jorb-1,norbOnAt)+1,1))
                        tempArr(1,iorb,jorb)=orbsOnAtom(ceiling(dble(iorb)/dble(norbOnAt)),1) ; &
                        tempArr(2,iorb,jorb)=orbsOnAtom(mod(jorb-1,norbOnAt)+1,1)
                    end if
                end do
            end do
        end if
        !!! Delete the elements belonging to a_ii
        !!if(iproc==0) then
        !!    do jorb=1,norbOnAt**2
        !!        do iorb=1,norbOnAt**2
        !!            if(mod(iorb-1,norbOnAt)+1==ceiling(dble(iorb)/dble(norbOnAt)) .and. iorb==jorb) then
        !!                emat(iorb,jorb)=0.d0
        !!            end if
        !!        end do
        !!    end do
        !!end if
        if(iproc==0) then
            do iorb=1,norbOnAt**2
                write(*,'(100f9.5)') (emat(iorb,jorb), jorb=1,norbOnAt**2)
            end do
        end if
        if(iproc==0) write(*,*) '--------------------------------'
        if(iproc==0) then
            do iorb=1,norbOnAt**2
                write(*,'(100(2i3,2x))') (tempArr(1:2,iorb,jorb), jorb=1,norbOnAt**2)
            end do
        end if

        ! Build the vector dvec
        if(allocated(dvec)) deallocate(dvec)
        if(iproc==0) then
            allocate(dvec(norbOnAt**2))
            dvec=0.d0
            ii=1
            do iorb=1,norbOnAt**2
                if(iorb==orbsOnAtom(ii,1)) then
                    jj=1
                    do jorb=1,norbOnAt**2
                        if(jorb==orbsOnAtom(jj,1)) then
                            tempArr(1,(ii-1)*norbOnAt+jj,1)=iorb ; tempArr(2,(ii-1)*norbOnAt+jj,1)=jorb
                            dvec((ii-1)*norbOnAt+jj)=d((iorb-1)*norbOnAt+jorb)
                            jj=jj+1
                        end if
                    end do
                    ii=ii+1
                end if
            end do
  
            do iorb=1,norbOnAt**2
                write(*,'(2i4)') tempArr(1,iorb,1), tempArr(2,iorb,1)
            end do
        end if


        ! Solve emat*a=dvec
        if(iproc==0) then
            allocate(a(norbOnAt**2))
            allocate(ipiv(norbOnAt**2))
            call dgesv(norbOnAt**2, 1, emat(1,1), norbOnAt**2, ipiv(1), dvec(1), norbOnAt**2, info)
            if(info/=0) then
                write(*,'(a,i0)') 'ERROR in dgesv, errorcode=', info
                do iorb=1,norbOnAt**2
                    write(*,'(100f9.5)') (emat(iorb,jorb), jorb=1,norbOnAt**2)
                end do
                stop
            end if
            do iorb=1,norbOnAt**2
                write(*,'(a,i4,es12.5)') 'iorb, dvec(iorb)', iorb, dvec(iorb)
            end do
        end if
        

        !! Build the matrix epsilon and the vector d only for those orbitals centerd in atom iat.
        !jj=0
        !ik=0
        !do jproc=0,nproc-1
        !    do jorb=1,orbsLIN%norb_par(jproc)
        !        jj=jj+1
        !        if(jproc==iproc) then
        !            if(orbsLIN%onWhichAtom(jorb)==iat) then
        !                kk=0
        !                il=0
        !                do kproc=0,nproc-1
        !                    do korb=1,orbsLIN%norb_par(kproc)
        !                        kk=kk+1
        !                        if(kproc==iproc) then
        !                            if(orbsLIN%onWhichAtom(korb)==iat) then
        !                                ik=ik+1
        !                                dvec(ik)=-2.d0*d(jj+kk)
        !                                ll=0
        !                                do lproc=0,nproc-1
        !                                    do lorb=1,orbsLIN%norb_par(lproc)
        !                                        ll=ll+1
        !                                        if(lproc==iproc) then
        !                                            if(orbsLIN%onWhichAtom(lorb)==iat) then
        !                                                il=il+1
        !                                                emat(ik,il)=eps((ik-1)*orbsLIN%norbp+il)
        !                                            end if
        !                                        end if
        !                                    end do
        !                                end do
        !                            end if
        !                        end if
        !                    end do
        !                end do
        !            end if
        !        end if
        !    end do
        !end do
        !do iorb=1,14
        !    write(*,'(20f9.4)') (emat(iorb,jorb), jorb=1,14)
        !end do
        !if(iproc==0) write(*,*) '****************************************************'
        !if(iproc==0) write(*,*) '****************************************************'
        !do iorb=1,14
        !    write(*,'(20f9.4)') dvec(iorb)
        !end do

        !do iorb=1,orbsLIN%norb**2
        !    if(iproc==0) write(*,'(f9.5)') d(iorb)
        !end do
        !tt=dnrm2(nvctrp, hphi(istart), 1)
        !call mpi_allreduce(tt, hphiNrm, 1, mpi_double_precision, mpi_sum, mpi_comm_world, ierr)
        !tt=0.d0
        !do jorb=1,orbsLIN%norb
        !    write(*,'(a,3i4,2es14.6)') 'iproc, iorb, jorb, eps(jorb), hphiNrm', iproc, iorb, jorb, eps(jorb), hphiNrm
        !    tt=tt+eps(jorb)**2
        !end do
        !write(*,'(a,2i4,2es14.6,f8.3)') 'iproc, iorb, epsSum, hphiNrm, sqrt(tt)/hphiNrm', iproc, iorb, sqrt(tt), hphiNrm, sqrt(tt)/hphiNrm
        !fracTot=fracTot+sqrt(tt)/hphiNrm
        ! Then calculate the b coefficients
        !!call getB(iproc, nproc, orbsLIN, commsLIN, eps, phi, hphi, iorb)
!!        phi2Loop: do jorb=1,orbsLIN%norb
!!
!!            
!!            
!!            ! Then calculate the d coefficients
!!            call getD()
!!  
!!            ! Calculate the force component
!!            call getForce()
!!        
!!        end do phi2Loop
          istart=istart+nvctrp
    !end do phi1Loop
    end do atomsLoop
    fracTot=fracTot/orbsLIN%norb
    write(*,'(a,i4,f8.3)') 'iproc, fracTot', iproc, fracTot
   tt=fracTot
   call mpi_allreduce(tt, fracTot, 1, mpi_double_precision, mpi_sum, mpi_comm_world, ierr)
   if(iproc==0) write(*,'(a,f8.3)') 'mean totFrac', fracTot/dble(nproc)
!!end do psiLoop
!!
!!
end subroutine pulay
!!
!!
!!
!!
subroutine getEpsilon(iproc, nproc, orbsLIN, commsLIN, eps, phi, hphi)

use module_base
use module_types
implicit none

! Calling arguments
integer:: iproc, nproc
type(orbitals_data), intent(inout):: orbsLIN
type(communications_arrays), intent(in):: commsLIN
real(8),dimension(orbsLIN%norb,orbsLIN%norb):: eps
real(8),dimension(orbsLIN%npsidim):: phi
real(8),dimension(sum(commsLIN%nvctr_par(iproc,1:orbsLIN%nkptsp))*orbsLIN%nspinor*orbsLIN%norb):: hphi

! Local variables
integer:: nvctrp, istat, ierr, istart, jstart, iorb, jorb
real(8),dimension(:,:),allocatable:: epsTemp 
real(8):: ddot

! nvctrp is the amount of each phi hold by the current process
nvctrp=sum(commsLIN%nvctr_par(iproc,1:orbsLIN%nkptsp))*orbsLIN%nspinor

allocate(epsTemp(orbsLIN%norb,orbsLIN%norb), stat=istat)
epsTemp=0.d0
istart=1
do iorb=1,orbsLIN%norb
    jstart=1
    do jorb=1,orbsLIN%norb
        epsTemp(iorb,jorb)=ddot(nvctrp, phi(istart), 1, hphi(jstart), 1)
        epsTemp(iorb,jorb)=epsTemp(iorb,jorb)+ddot(nvctrp, phi(jstart), 1, hphi(istart), 1)
        jstart=jstart+nvctrp
    end do
    istart=istart+nvctrp
end do
eps=0.d0
call mpi_allreduce(epsTemp(1,1), eps(1,1), orbsLIN%norb**2, mpi_double_precision, mpi_sum, mpi_comm_world, ierr)

end subroutine getEpsilon






!!!!subroutine getD(iproc, nproc, orbsLIN, commsLIN, d, phi)
!!!!
!!!!use module_base
!!!!use module_types
!!!!implicit none
!!!!
!!!!! Calling arguments
!!!!integer:: iproc, nproc
!!!!type(orbitals_data), intent(inout):: orbsLIN
!!!!type(communications_arrays), intent(in):: commsLIN
!!!!real(8),dimension(orbsLIN%norb*orbsLIN%norb):: d
!!!!real(8),dimension(orbsLIN%npsidim):: phi
!!!!
!!!!! Local variables
!!!!integer:: nvctrp, istat, ierr, istart, jstart, iorb, jorb
!!!!real(8),dimension(:),allocatable:: dTemp 
!!!!real(8):: ddot
!!!!
!!!!! nvctrp is the amount of each phi hold by the current process
!!!!nvctrp=sum(commsLIN%nvctr_par(iproc,1:orbsLIN%nkptsp))*orbsLIN%nspinor
!!!!
!!!!allocate(dTemp(orbsLIN%norb*orbsLIN%norb), stat=istat)
!!!!istart=1
!!!!do iorb=1,orbsLIN%norb
!!!!    jstart=1
!!!!    do jorb=1,orbsLIN%norb
!!!!        !call calculateMatrixElement()
!!!!        jstart=jstart+nvctrp
!!!!    end do
!!!!    istart=istart+nvctrp
!!!!end do
!!!!call mpi_allreduce(dTemp(1), d(1), orbsLIN%norb, mpi_double_precision, mpi_sum, mpi_comm_world, ierr)
!!!!
!!!!end subroutine getD



subroutine getB(iproc, nproc, orbsLIN, commsLIN, e, phi, hphi, iorb)

use module_base
use module_types
implicit none

! Calling arguments
integer:: iproc, nproc, iorb
type(orbitals_data), intent(inout):: orbsLIN
type(communications_arrays), intent(in):: commsLIN
real(8),dimension(orbsLIN%norb):: e
real(8),dimension(orbsLIN%npsidim):: phi, hphi

! Local variables
integer:: jorb

do jorb=1,orbsLIN%norb
    if(jorb==iorb) cycle
    e(jorb)=e(iorb)/(e(iorb)-e(jorb))
end do

end subroutine getB





!!subroutine getE
!!implicit none
!!
!!
!!
!!call HamiltonianApplicationParabola(iproc,nproc,at,orbsLIN,input%hx,input%hy,input%hz,rxyz,&
!!     nlpspd,proj,Glr,ngatherarr,Glr%d%n1i*Glr%d%n2i*nscatterarr(iproc,2),&
!!     rhopot(1),&
!!     phi(1),hphi(1),ekin_sum,epot_sum,eexctX,eproj_sum,nspin,GPU, orbsLIN%onWhichAtom, rxyzParabola, pkernel=pkernelseq)
!!
!!
!!
!!end subroutine getE





subroutine inputOrbitals(iproc,nproc,at,&
     orbs,nvirt,comms,Glr,hx,hy,hz,rxyz,rhopot,rhocore,pot_ion,&
     nlpspd,proj,pkernel,pkernelseq,ixc,psi,hpsi,psit,G,&
     nscatterarr,ngatherarr,nspin,potshortcut,symObj,irrzon,phnons,GPU,input)
  ! Input wavefunctions are found by a diagonalization in a minimal basis set
  ! Each processors write its initial wavefunctions into the wavefunction file
  ! The files are then read by readwave
  use module_base
  use module_interfaces, except_this_one => inputOrbitals
  use module_types
  use Poisson_Solver
  implicit none
  !Arguments
  integer, intent(in) :: iproc,nproc,ixc,symObj
  integer, intent(inout) :: nspin,nvirt
  real(gp), intent(in) :: hx,hy,hz
  type(atoms_data), intent(in) :: at
  type(orbitals_data), intent(inout) :: orbs
  type(nonlocal_psp_descriptors), intent(in) :: nlpspd
  type(locreg_descriptors), intent(in) :: Glr
  type(communications_arrays), intent(in) :: comms
  type(GPU_pointers), intent(inout) :: GPU
  type(input_variables):: input
  integer, dimension(0:nproc-1,4), intent(in) :: nscatterarr !n3d,n3p,i3s+i3xcsh-1,i3xcsh
  integer, dimension(0:nproc-1,2), intent(in) :: ngatherarr 
  real(gp), dimension(3,at%nat), intent(in) :: rxyz
  real(wp), dimension(nlpspd%nprojel), intent(in) :: proj
  real(dp), dimension(*), intent(inout) :: rhopot,pot_ion
  type(gaussian_basis), intent(out) :: G !basis for davidson IG
  real(wp), dimension(:), pointer :: hpsi,psit,rhocore
  real(8),dimension(orbs%npsidim):: psi
  real(dp), dimension(:), pointer :: pkernel,pkernelseq
  integer, intent(in) ::potshortcut
  integer, dimension(*), intent(in) :: irrzon
  real(dp), dimension(*), intent(in) :: phnons
  !local variables
  character(len=*), parameter :: subname='input_wf_diag'
  logical :: switchGPUconv,switchOCLconv
  integer :: i_stat,i_all,iat,nspin_ig,iorb,idum=0
  real(kind=4) :: tt,builtin_rand
  real(gp) :: hxh,hyh,hzh,eks,eexcu,vexcu,epot_sum,ekin_sum,ehart,eexctX,eproj_sum,etol,accurex
  type(orbitals_data) :: orbse
  type(communications_arrays) :: commse
  integer, dimension(:,:), allocatable :: norbsc_arr
  real(wp), dimension(:), allocatable :: potxc
  real(gp), dimension(:), allocatable :: locrad
  type(locreg_descriptors), dimension(:), allocatable :: Llr
  real(wp), dimension(:,:,:), pointer :: psigau
type(orbitals_data):: orbsLIN
type(communications_arrays):: commsLIN
real(8),dimension(:),allocatable:: phi, hphi
real(8),dimension(:,:),allocatable:: HamSmall
real(8),dimension(:),allocatable:: eval
integer:: istat
real(8),dimension(:),pointer:: phiWorkPointer

  allocate(norbsc_arr(at%natsc+1,nspin+ndebug),stat=i_stat)
  call memocc(i_stat,norbsc_arr,'norbsc_arr',subname)
  allocate(locrad(at%nat+ndebug),stat=i_stat)
  call memocc(i_stat,locrad,'locrad',subname)

  write(*,*) 'in inputOrbitals'

  if (iproc == 0) then
     write(*,'(1x,a)')&
          '------------------------------------------------------- Input Wavefunctions Creation'
  end if

  !spin for inputguess orbitals
  if (nspin == 4) then
     nspin_ig=1
  else
     nspin_ig=nspin
  end if

  call inputguess_gaussian_orbitals(iproc,nproc,at,rxyz,Glr,nvirt,nspin_ig,&
       orbs,orbse,norbsc_arr,locrad,G,psigau,eks)

  !allocate communications arrays for inputguess orbitals
  !call allocate_comms(nproc,orbse,commse,subname)
  call orbitals_communicators(iproc,nproc,Glr,orbse,commse)  

  hxh=.5_gp*hx
  hyh=.5_gp*hy
  hzh=.5_gp*hz

  !check the communication distribution
  !call check_communications(iproc,nproc,orbse,Glr,commse)

  !once the wavefunction coefficients are known perform a set 
  !of nonblocking send-receive operations to calculate overlap matrices

!!!  !create mpirequests array for controlling the success of the send-receive operation
!!!  allocate(mpirequests(nproc-1+ndebug),stat=i_stat)
!!!  call memocc(i_stat,mpirequests,'mpirequests',subname)
!!!
!!!  call nonblocking_transposition(iproc,nproc,G%ncoeff,orbse%isorb+orbse%norbp,&
!!!       orbse%nspinor,psigau,orbse%norb_par,mpirequests)

  !experimental part for building the localisation regions
  if (at%geocode == 'F') then
     !allocate the array of localisation regions
     allocate(Llr(at%nat+ndebug),stat=i_stat)
     !call memocc(i_stat,Llr,'Llr',subname)

     !print *,'locrad',locrad

     call determine_locreg(at%nat,rxyz,locrad,hx,hy,hz,Glr,Llr)

     do iat=1,at%nat
        call deallocate_lr(Llr(iat),subname)
!!$        call deallocate_wfd(Llr(iat)%wfd,subname)
!!$        if (Llr(iat)%geocode=='F') then
!!$           call deallocate_bounds(Llr(iat)%bounds,subname)
!!$        end if
     end do

     !i_all=-product(shape(Llr))*kind(Llr)
     deallocate(Llr,stat=i_stat) !these allocation are special
     !call memocc(i_stat,i_all,'Llr',subname)
  end if

!!!  !allocate the wavefunction in the transposed way to avoid allocations/deallocations
!!!  allocate(psi(orbse%npsidim+ndebug),stat=i_stat)
!!!  call memocc(i_stat,psi,'psi',subname)

  !allocate arrays for the GPU if a card is present
  switchGPUconv=.false.
  switchOCLconv=.false.
  if (GPUconv .and. potshortcut ==0 ) then
     call prepare_gpu_for_locham(Glr%d%n1,Glr%d%n2,Glr%d%n3,nspin,&
          hx,hy,hz,Glr%wfd,orbse,GPU)
  else if (OCLconv .and. potshortcut ==0) then
     call allocate_data_OCL(Glr%d%n1,Glr%d%n2,Glr%d%n3,at%geocode,&
          nspin,hx,hy,hz,Glr%wfd,orbse,GPU)
     if (iproc == 0) write(*,*)&
          'GPU data allocated'
  else if (GPUconv .and. potshortcut >0 ) then
     switchGPUconv=.true.
     GPUconv=.false.
  else if (OCLconv .and. potshortcut >0 ) then
     switchOCLconv=.true.
     OCLconv=.false.
  end if

write(50+iproc,*) psigau(:,:,orbse%isorb+1:orbse%isorb+orbse%norbp)

  !use only the part of the arrays for building the hamiltonian matrix
  call gaussians_to_wavelets_new(iproc,nproc,Glr,orbse,hx,hy,hz,G,&
       psigau(1,1,min(orbse%isorb+1,orbse%norb)),psi)


  i_all=-product(shape(locrad))*kind(locrad)
  deallocate(locrad,stat=i_stat)
  call memocc(i_stat,i_all,'locrad',subname)

!!!!  !application of the hamiltonian for gaussian based treatment
!!!!  call sumrho(iproc,nproc,orbse,Glr,ixc,hxh,hyh,hzh,psi,rhopot,&
!!!!       & Glr%d%n1i*Glr%d%n2i*nscatterarr(iproc,1),nscatterarr,nspin,GPU, &
!!!!       & symObj, irrzon, phnons)
!!!!     
!!!!  !-- if spectra calculation uses a energy dependent potential
!!!!  !    input_wf_diag will write (to be used it in abscalc)
!!!!  !    the density to the file electronic_density.cube
!!!!  !  The writing is activated if  5th bit of  in%potshortcut is on.
!!!!  if( iand( potshortcut,16)==0 .and. potshortcut /= 0) then
!!!!     call plot_density_cube_old(at%geocode,'electronic_density',&
!!!!          iproc,nproc,Glr%d%n1,Glr%d%n2,Glr%d%n3,Glr%d%n1i,Glr%d%n2i,Glr%d%n3i,nscatterarr(iproc,2),  & 
!!!!          nspin,hxh,hyh,hzh,at,rxyz,ngatherarr,rhopot(1+nscatterarr(iproc,4)*Glr%d%n1i*Glr%d%n2i))
!!!!  endif
!!!!  !---
!!!!  
!!!!  if(orbs%nspinor==4) then
!!!!     !this wrapper can be inserted inside the poisson solver 
!!!!     call PSolverNC(at%geocode,'D',iproc,nproc,Glr%d%n1i,Glr%d%n2i,Glr%d%n3i,&
!!!!          nscatterarr(iproc,1),& !this is n3d
!!!!          ixc,hxh,hyh,hzh,&
!!!!          rhopot,pkernel,pot_ion,ehart,eexcu,vexcu,0.d0,.true.,4)
!!!!  else
!!!!     !Allocate XC potential
!!!!     if (nscatterarr(iproc,2) >0) then
!!!!        allocate(potxc(Glr%d%n1i*Glr%d%n2i*nscatterarr(iproc,2)*nspin+ndebug),stat=i_stat)
!!!!        call memocc(i_stat,potxc,'potxc',subname)
!!!!     else
!!!!        allocate(potxc(1+ndebug),stat=i_stat)
!!!!        call memocc(i_stat,potxc,'potxc',subname)
!!!!     end if
!!!!
!!!!     call XC_potential(at%geocode,'D',iproc,nproc,&
!!!!          Glr%d%n1i,Glr%d%n2i,Glr%d%n3i,ixc,hxh,hyh,hzh,&
!!!!          rhopot,eexcu,vexcu,nspin,rhocore,potxc)
!!!!
!!!!
!!!!     if( iand(potshortcut,4)==0) then
!!!!        call H_potential(at%geocode,'D',iproc,nproc,&
!!!!             Glr%d%n1i,Glr%d%n2i,Glr%d%n3i,hxh,hyh,hzh,&
!!!!             rhopot,pkernel,pot_ion,ehart,0.0_dp,.true.)
!!!!     endif
!!!!
!!!!
!!!!     !sum the two potentials in rhopot array
!!!!     !fill the other part, for spin, polarised
!!!!     if (nspin == 2) then
!!!!        call dcopy(Glr%d%n1i*Glr%d%n2i*nscatterarr(iproc,2),rhopot(1),1,&
!!!!             rhopot(Glr%d%n1i*Glr%d%n2i*nscatterarr(iproc,2)+1),1)
!!!!     end if
!!!!     !spin up and down together with the XC part
!!!!     call axpy(Glr%d%n1i*Glr%d%n2i*nscatterarr(iproc,2)*nspin,1.0_dp,potxc(1),1,&
!!!!          rhopot(1),1)
!!!!
!!!!
!!!!     i_all=-product(shape(potxc))*kind(potxc)
!!!!     deallocate(potxc,stat=i_stat)
!!!!     call memocc(i_stat,i_all,'potxc',subname)
!!!!
!!!!  end if
!!!!
!!!!!!!  if (nproc == 1) then
!!!!!!!     !calculate the overlap matrix as well as the kinetic overlap
!!!!!!!     !in view of complete gaussian calculation
!!!!!!!     allocate(ovrlp(G%ncoeff*G%ncoeff),stat=i_stat)
!!!!!!!     call memocc(i_stat,ovrlp,'ovrlp',subname)
!!!!!!!     allocate(tmp(G%ncoeff,orbse%norb),stat=i_stat)
!!!!!!!     call memocc(i_stat,tmp,'tmp',subname)
!!!!!!!     allocate(smat(orbse%norb,orbse%norb),stat=i_stat)
!!!!!!!     call memocc(i_stat,smat,'smat',subname)
!!!!!!!
!!!!!!!     !overlap calculation of the gaussian matrix
!!!!!!!     call gaussian_overlap(G,G,ovrlp)
!!!!!!!     call dsymm('L','U',G%ncoeff,orbse%norb,1.0_gp,ovrlp(1),G%ncoeff,&
!!!!!!!          gaucoeff(1,1),G%ncoeff,0.d0,tmp(1,1),G%ncoeff)
!!!!!!!
!!!!!!!     call gemm('T','N',orbse%norb,orbse%norb,G%ncoeff,1.0_gp,&
!!!!!!!          gaucoeff(1,1),G%ncoeff,tmp(1,1),G%ncoeff,0.0_wp,smat(1,1),orbse%norb)
!!!!!!!
!!!!!!!     !print overlap matrices
!!!!!!!     do i=1,orbse%norb
!!!!!!!        write(*,'(i5,30(1pe15.8))')i,(smat(i,iorb),iorb=1,orbse%norb)
!!!!!!!     end do
!!!!!!!
!!!!!!!     !overlap calculation of the kinetic operator
!!!!!!!     call kinetic_overlap(G,G,ovrlp)
!!!!!!!     call dsymm('L','U',G%ncoeff,orbse%norb,1.0_gp,ovrlp(1),G%ncoeff,&
!!!!!!!          gaucoeff(1,1),G%ncoeff,0.d0,tmp(1,1),G%ncoeff)
!!!!!!!
!!!!!!!     call gemm('T','N',orbse%norb,orbse%norb,G%ncoeff,1.0_gp,&
!!!!!!!          gaucoeff(1,1),G%ncoeff,tmp(1,1),G%ncoeff,0.0_wp,smat(1,1),orbse%norb)
!!!!!!!
!!!!!!!     !print overlap matrices
!!!!!!!     tt=0.0_wp
!!!!!!!     do i=1,orbse%norb
!!!!!!!        write(*,'(i5,30(1pe15.8))')i,(smat(i,iorb),iorb=1,orbse%norb)
!!!!!!!        !write(12,'(i5,30(1pe15.8))')i,(smat(i,iorb),iorb=1,orbse%norb)
!!!!!!!        tt=tt+smat(i,i)
!!!!!!!     end do
!!!!!!!     print *,'trace',tt
!!!!!!!
!!!!!!!     !overlap calculation of the kinetic operator
!!!!!!!     call cpu_time(t0)
!!!!!!!     call potential_overlap(G,G,rhopot,Glr%d%n1i,Glr%d%n2i,Glr%d%n3i,hxh,hyh,hzh,&
!!!!!!!          ovrlp)
!!!!!!!     call cpu_time(t1)
!!!!!!!     call dsymm('L','U',G%ncoeff,orbse%norb,1.0_gp,ovrlp(1),G%ncoeff,&
!!!!!!!          gaucoeff(1,1),G%ncoeff,0.d0,tmp(1,1),G%ncoeff)
!!!!!!!
!!!!!!!     call gemm('T','N',orbse%norb,orbse%norb,G%ncoeff,1.0_gp,&
!!!!!!!          gaucoeff(1,1),G%ncoeff,tmp(1,1),G%ncoeff,0.0_wp,smat(1,1),orbse%norb)
!!!!!!!
!!!!!!!     !print overlap matrices
!!!!!!!     tt=0.0_wp
!!!!!!!     do i=1,orbse%norb
!!!!!!!        write(*,'(i5,30(1pe15.8))')i,(smat(i,iorb),iorb=1,orbse%norb)
!!!!!!!        !write(12,'(i5,30(1pe15.8))')i,(smat(i,iorb),iorb=1,orbse%norb)
!!!!!!!        tt=tt+smat(i,i)
!!!!!!!     end do
!!!!!!!     print *,'trace',tt
!!!!!!!     print *, 'time',t1-t0
!!!!!!!
!!!!!!!     i_all=-product(shape(ovrlp))*kind(ovrlp)
!!!!!!!     deallocate(ovrlp,stat=i_stat)
!!!!!!!     call memocc(i_stat,i_all,'ovrlp',subname)
!!!!!!!     i_all=-product(shape(tmp))*kind(tmp)
!!!!!!!     deallocate(tmp,stat=i_stat)
!!!!!!!     call memocc(i_stat,i_all,'tmp',subname)
!!!!!!!     i_all=-product(shape(smat))*kind(smat)
!!!!!!!     deallocate(smat,stat=i_stat)
!!!!!!!     call memocc(i_stat,i_all,'smat',subname)
!!!!!!!  end if
!!!!
!!!!  if(potshortcut>0) then
!!!!!!$    if (GPUconv) then
!!!!!!$       call free_gpu(GPU,orbs%norbp)
!!!!!!$    end if
!!!!     if (switchGPUconv) then
!!!!        GPUconv=.true.
!!!!     end if
!!!!     if (switchOCLconv) then
!!!!        OCLconv=.true.
!!!!     end if
!!!!
!!!!     call deallocate_orbs(orbse,subname)
!!!!     
!!!!     !deallocate the gaussian basis descriptors
!!!!     call deallocate_gwf(G,subname)
!!!!    
!!!!     i_all=-product(shape(psigau))*kind(psigau)
!!!!     deallocate(psigau,stat=i_stat)
!!!!     call memocc(i_stat,i_all,'psigau',subname)
!!!!     call deallocate_comms(commse,subname)
!!!!     i_all=-product(shape(norbsc_arr))*kind(norbsc_arr)
!!!!     deallocate(norbsc_arr,stat=i_stat)
!!!!     call memocc(i_stat,i_all,'norbsc_arr',subname)
!!!!    return 
!!!!  end if
!!!!
!!!!  !allocate the wavefunction in the transposed way to avoid allocations/deallocations
!!!!  allocate(hpsi(orbse%npsidim+ndebug),stat=i_stat)
!!!!  call memocc(i_stat,hpsi,'hpsi',subname)
!!!!
!!!!  !call dcopy(orbse%npsidim,psi,1,hpsi,1)
!!!!  if (input%exctxpar == 'OP2P') eexctX = -99.0_gp
!!!!
!!!!
!!!!
!!!!! THIS IS NOW DONE IN CLUSTER
!!!!!call initializeParameters(iproc, nproc, Glr, orbs, orbsLIN, commsLIN, at, phi)
!!!!!call improveOrbitals(iproc, nproc, nspin, Glr, orbs, orbsLIN, commsLIN, at, rxyz, nscatterarr, ngatherarr, nlpspd, proj, rhopot, GPU, input, pkernelseq, phi)
!!!!
!!!!
!!!!!!!write(*,*) 'calling getLocalizedBasis'
!!!!!!!call initializeParameters(iproc, nproc, Glr, orbs, orbsLIN, commsLIN, at)
!!!!!!!allocate(phi(orbsLIN%npsidim), stat=istat)
!!!!!!!call initRandomSeed(iproc, 1)
!!!!!!!call random_number(phi)
!!!!!!!allocate(hphi(orbsLIN%npsidim), stat=istat)
!!!!!!!! Initialize phi at random
!!!!!!!call getLocalizedBasis(iproc, nproc, at, orbs, Glr, input, orbsLIN, commsLIN, rxyz, nspin, nlpspd, proj, nscatterarr, ngatherarr, rhopot, GPU, pkernelseq, phi, hphi)
!!!!!!!allocate(phiWorkPointer(size(phi)), stat=istat)
!!!!!!!call transpose_v(iproc, nproc, orbsLIN, Glr%wfd, commsLIN, phi, work=phiWorkPointer)
!!!!!!!call transpose_v(iproc, nproc, orbsLIN, Glr%wfd, commsLIN, hphi, work=phiWorkPointer)
!!!!!!!allocate(HamSmall(orbsLIN%norb,orbsLIN%norb), stat=istat)
!!!!!!!call transformHam(iproc, nproc, orbsLIN, commsLIN, phi, hphi, HamSmall)
!!!!!!!allocate(eval(orbsLIN%norb), stat=istat)
!!!!!!!call diagonalizeHamiltonian(iproc, nproc, orbsLIN, HamSmall, eval)
!!!!!!!
!!!!!!!! Store the new wave function in hphi as a temporary array
!!!!!!!call buildWavefunction(iproc, nproc, orbsLIN, commsLIN, phi, hphi, HamSmall)
!!!!!!!call dcopy(orbsLIN%npsidim, hphi(1), 1, phi(1), 1)
!!!!!!!call untranspose_v(iproc, nproc, orbsLIN, Glr%wfd, commsLIN, phi, work=phiWorkPointer)
!!!!!!!! Not necessary to untranspose hphi (no longer needed)
!!!!!!!!call untranspose_v(iproc, nproc, orbsLIN, Glr%wfd, commsLIN, hphi, work=phiWorkPointer)
!!!!!!!deallocate(phiWorkPointer, stat=istat)
!!!!
!!!!  call HamiltonianApplication(iproc,nproc,at,orbse,hx,hy,hz,rxyz,&
!!!!       nlpspd,proj,Glr,ngatherarr,Glr%d%n1i*Glr%d%n2i*nscatterarr(iproc,2),&
!!!!       rhopot,&
!!!!       psi,hpsi,ekin_sum,epot_sum,eexctX,eproj_sum,nspin,GPU,pkernel=pkernelseq)
!!!!
!!!!!!!  !calculate the overlap matrix knowing that the original functions are gaussian-based
!!!!!!!  allocate(thetaphi(2,G%nat+ndebug),stat=i_stat)
!!!!!!!  call memocc(i_stat,thetaphi,'thetaphi',subname)
!!!!!!!  thetaphi=0.0_gp
!!!!!!!
!!!!!!!  !calculate the scalar product between the hamiltonian and the gaussian basis
!!!!!!!  allocate(hpsigau(G%ncoeff,orbse%norbp+ndebug),stat=i_stat)
!!!!!!!  call memocc(i_stat,hpsigau,'hpsigau',subname)
!!!!!!!
!!!!!!!
!!!!!!!  call wavelets_to_gaussians(at%geocode,orbse%norbp,Glr%d%n1,Glr%d%n2,Glr%d%n3,G,&
!!!!!!!       thetaphi,hx,hy,hz,Glr%wfd,hpsi,hpsigau)
!!!!!!!
!!!!!!!  i_all=-product(shape(thetaphi))*kind(thetaphi)
!!!!!!!  deallocate(thetaphi,stat=i_stat)
!!!!!!!  call memocc(i_stat,i_all,'thetaphi',subname)
!!!!
!!!!  accurex=abs(eks-ekin_sum)
!!!!  !tolerance for comparing the eigenvalues in the case of degeneracies
!!!!  etol=accurex/real(orbse%norbu,gp)
!!!!  if (iproc == 0 .and. verbose > 1) write(*,'(1x,a,2(f19.10))') 'done. ekin_sum,eks:',ekin_sum,eks
!!!!  if (iproc == 0) then
!!!!     write(*,'(1x,a,3(1x,1pe18.11))') 'ekin_sum,epot_sum,eproj_sum',  & 
!!!!          ekin_sum,epot_sum,eproj_sum
!!!!     write(*,'(1x,a,3(1x,1pe18.11))') '   ehart,   eexcu,    vexcu',ehart,eexcu,vexcu
!!!!  endif
!!!!
!!!!!!!  call Gaussian_DiagHam(iproc,nproc,at%natsc,nspin,orbs,G,mpirequests,&
!!!!!!!       psigau,hpsigau,orbse,etol,norbsc_arr)
!!!!
!!!!
!!!!!!!  i_all=-product(shape(mpirequests))*kind(mpirequests)
!!!!!!!  deallocate(mpirequests,stat=i_stat)
!!!!!!!  call memocc(i_stat,i_all,'mpirequests',subname)
!!!!
!!!!!!!  i_all=-product(shape(hpsigau))*kind(hpsigau)
!!!!!!!  deallocate(hpsigau,stat=i_stat)
!!!!!!!  call memocc(i_stat,i_all,'hpsigau',subname)
!!!!
!!!!  !free GPU if it is the case
!!!!  if (GPUconv) then
!!!!     call free_gpu(GPU,orbse%norbp)
!!!!  else if (OCLconv) then
!!!!     call free_gpu_OCL(GPU,orbse%norbp)
!!!!  end if
!!!!
!!!!  if (iproc == 0 .and. verbose > 1) write(*,'(1x,a)',advance='no')&
!!!!       'Input Wavefunctions Orthogonalization:'
!!!!
!!!!  !psivirt can be eliminated here, since it will be allocated before davidson
!!!!  !with a gaussian basis
!!!!!!$  call DiagHam(iproc,nproc,at%natsc,nspin_ig,orbs,Glr%wfd,comms,&
!!!!!!$       psi,hpsi,psit,orbse,commse,etol,norbsc_arr,orbsv,psivirt)
!!!!
!!!!  call DiagHam(iproc,nproc,at%natsc,nspin_ig,orbs,Glr%wfd,comms,&
!!!!       psi,hpsi,psit,input,orbse,commse,etol,norbsc_arr)
!!!!
!!!!  if (input%itrpmax > 1) then
!!!!     !use the eval array of orbse structure to save the original values
!!!!     allocate(orbse%eval(orbs%norb*orbs%nkpts+ndebug),stat=i_stat)
!!!!     call memocc(i_stat,orbse%eval,'orbse%eval',subname)
!!!!     
!!!!     call dcopy(orbs%norb*orbs%nkpts,orbs%eval(1),1,orbse%eval(1),1)
!!!!
!!!!     !add a small displacement in the eigenvalues
!!!!     do iorb=1,orbs%norb*orbs%nkpts
!!!!        tt=builtin_rand(idum)
!!!!        orbs%eval(iorb)=orbs%eval(iorb)*(1.0_gp+0.05_gp*real(tt,gp))
!!!!     end do
!!!!
!!!!     !correct the occupation numbers wrt fermi level
!!!!     call Fermilevel(.false.,1.e-2_gp,orbs)
!!!!
!!!!     !restore the occupation numbers
!!!!     call dcopy(orbs%norb*orbs%nkpts,orbse%eval(1),1,orbs%eval(1),1)
!!!!
!!!!     i_all=-product(shape(orbse%eval))*kind(orbse%eval)
!!!!     deallocate(orbse%eval,stat=i_stat)
!!!!     call memocc(i_stat,i_all,'orbse%eval',subname)
!!!!  end if
!!!!
!!!!  call deallocate_comms(commse,subname)
!!!!
!!!!  i_all=-product(shape(norbsc_arr))*kind(norbsc_arr)
!!!!  deallocate(norbsc_arr,stat=i_stat)
!!!!  call memocc(i_stat,i_all,'norbsc_arr',subname)
!!!!
!!!!  if (iproc == 0) then
!!!!     if (verbose > 1) write(*,'(1x,a)')'done.'
!!!!     !gaussian estimation valid only for Free BC
!!!!     if (at%geocode == 'F') then
!!!!        write(*,'(1x,a,1pe9.2)') 'expected accuracy in energy ',accurex
!!!!        write(*,'(1x,a,1pe9.2)') &
!!!!          'expected accuracy in energy per orbital ',accurex/real(orbs%norb,kind=8)
!!!!        !write(*,'(1x,a,1pe9.2)') &
!!!!        !     'suggested value for gnrm_cv ',accurex/real(orbs%norb,kind=8)
!!!!     end if
!!!!  endif
!!!!
!!!!  !here we can define the subroutine which generates the coefficients for the virtual orbitals
!!!!  call deallocate_gwf(G,subname)
!!!!
!!!!  i_all=-product(shape(psigau))*kind(psigau)
!!!!  deallocate(psigau,stat=i_stat)
!!!!  call memocc(i_stat,i_all,'psigau',subname)
!!!!
!!!!  call deallocate_orbs(orbse,subname)



END SUBROUTINE inputOrbitals



!! THIS IS NOT WORKING !!
subroutine getD(iproc,nproc,orbsLIN,lr,hx,hy,hz,phi, hphi, nat, rxyzParab, parabPrefac, at, commsLIN, d)
  use module_base
  use module_types
  use module_interfaces
  implicit none
  integer, intent(in) :: iproc,nproc
  real(gp), intent(in) :: hx,hy,hz
  type(locreg_descriptors), intent(in) :: lr
  type(orbitals_data), intent(in) :: orbsLIN
  !real(wp), dimension(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f,orbs%nspinor,orbs%norbp), intent(inout) :: hpsi
  real(wp), dimension(orbsLIN%npsidim), intent(inout) :: phi, hphi
integer,intent(in):: nat
real(8),dimension(3),intent(in):: rxyzParab
real(8):: parabPrefac
type(atoms_data), intent(in) :: at
type(communications_arrays), intent(in):: commsLIN
real(8),dimension(orbsLIN%norb,orbsLIN%norb):: d
  !local variables
  integer :: iorb,inds,ncplx,ikpt,ierr
  real(wp) :: cprecr,scpr,eval_zero,evalmax 
  real(gp) :: kx,ky,kz
real(8),dimension(:),pointer:: phiWorkPointer
integer:: iiAt, istat, nvctrp, istart, jstart, jorb, it
real(8):: ddot
real(8),dimension(:,:),allocatable:: dTemp

  ncplx=1
  it=1
  ! First apply the perturbation to all orbitals belonging to iproc
  allocate(phiWorkPointer(size(phi)), stat=istat)
  call untranspose_v(iproc, nproc, orbsLIN, lr%wfd, commsLIN, phi, work=phiWorkPointer)
  hphi=phi
  istart=1
  do iorb=1,orbsLIN%norbp
      !iiAt=onWhichAtom(iorb)
      !parabPrefac=orbsLIN%parabPrefacArr(at%iatype(onWhichAtom(iorb)))


 !! THIS IS NOT WORKING !!
      !!call subTemp(lr,ncplx,&
      !!     hx,hy,hz,hphi(istart), rxyzParab, orbsLIN, parabPrefac, it)


      istart=istart+(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f)*orbsLIN%nspinor
  end do

  call transpose_v(iproc, nproc, orbsLIN, lr%wfd, commsLIN, phi, work=phiWorkPointer)
  call transpose_v(iproc, nproc, orbsLIN, lr%wfd, commsLIN, hphi, work=phiWorkPointer)


  ! nvctrp is the amount of each phi hold by the current process
  nvctrp=sum(commsLIN%nvctr_par(iproc,1:orbsLIN%nkptsp))*orbsLIN%nspinor
  
  allocate(dTemp(orbsLIN%norb,orbsLIN%norb), stat=istat)
  dTemp=0.d0
  istart=1
  do iorb=1,orbsLIN%norb
      jstart=1
      do jorb=1,orbsLIN%norb
          dTemp(iorb,jorb)=ddot(nvctrp, phi(istart), 1, hphi(jstart), 1)
          jstart=jstart+nvctrp
      end do
      istart=istart+nvctrp
  end do
  d=0.d0
  call mpi_allreduce(dTemp(1,1), d(1,1), orbsLIN%norb**2, mpi_double_precision, mpi_sum, mpi_comm_world, ierr)

  call untranspose_v(iproc, nproc, orbsLIN, lr%wfd, commsLIN, phi, work=phiWorkPointer)
  call untranspose_v(iproc, nproc, orbsLIN, lr%wfd, commsLIN, hphi, work=phiWorkPointer)
  deallocate(phiWorkPointer, stat=istat)


END SUBROUTINE getD
!!***




subroutine subTemp(lr,ncplx,&
     hx,hy,hz,x,  rxyzParab, orbs, parabPrefac, it, direction, shift)
  use module_base
  use module_types
  ! Solves (KE+cprecr*I)*xx=yy by conjugate gradient method
  ! x is the right hand side on input and the solution on output
  implicit none
  integer, intent(in) :: ncplx
  real(gp), intent(in) :: hx,hy,hz
  type(locreg_descriptors), intent(in) :: lr
  real(wp), dimension((lr%wfd%nvctr_c+7*lr%wfd%nvctr_f)*ncplx), intent(inout) :: x
real(8),dimension(3),intent(in):: rxyzParab
type(orbitals_data), intent(in) :: orbs
real(8):: parabPrefac
integer:: it
character(len=1):: direction
real(8):: shift
  ! local variables
  character(len=*), parameter :: subname='precondition_residue'
  real(gp), dimension(0:7) :: scal
  real(wp) :: rmr_old,rmr_new,alpha,beta
  integer :: i_stat,i_all,icong
  type(workarr_precond) :: w
  real(wp), dimension(:), allocatable :: b,r,d
real(8):: dnrm2

  !arrays for the CG procedure
  allocate(b(ncplx*(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f)+ndebug),stat=i_stat)
  call memocc(i_stat,b,'b',subname)
  allocate(r(ncplx*(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f)+ndebug),stat=i_stat)
  call memocc(i_stat,r,'r',subname)
  allocate(d(ncplx*(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f)+ndebug),stat=i_stat)
  call memocc(i_stat,d,'d',subname)

  call allocate_work_arrays(lr%geocode,lr%hybrid_on,ncplx,lr%d,w)

  !call precondition_preconditioner(lr,ncplx,hx,hy,hz,scal,cprecr,w,x,b)

  scal=1.d0
!write(*,*) 'dnrm2(x) 1' , dnrm2(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f, x, 1)
  call callApplyLinearOperator(ncplx,lr,hx,hy,hz,x,d,w,scal, rxyzParab, orbs, parabPrefac, it, direction, shift)
  x=d
!write(*,*) 'dnrm2(x)', dnrm2(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f, x, 1)

  i_all=-product(shape(b))*kind(b)
  deallocate(b,stat=i_stat)
  call memocc(i_stat,i_all,'b',subname)
  i_all=-product(shape(r))*kind(r)
  deallocate(r,stat=i_stat)
  call memocc(i_stat,i_all,'r',subname)
  i_all=-product(shape(d))*kind(d)
  deallocate(d,stat=i_stat)
  call memocc(i_stat,i_all,'d',subname)

  call deallocate_work_arrays(lr%geocode,lr%hybrid_on,ncplx,w)

END SUBROUTINE subTemp






subroutine callApplyLinearOperator(ncplx,lr,hx,hy,hz,&
     x,y,w,scal, rxyzParab, orbs, parabPrefac, it, direction, shift)! y:=Ax
  use module_base
  use module_types
  implicit none
  integer, intent(in) :: ncplx
  real(gp), intent(in) :: hx,hy,hz
  type(locreg_descriptors), intent(in) :: lr
  real(gp), dimension(0:7), intent(in) :: scal
  real(wp), dimension(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f,ncplx), intent(in) ::  x
  type(workarr_precond), intent(inout) :: w
  real(wp), dimension(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f,ncplx), intent(out) ::  y
real(8),dimension(3),intent(in):: rxyzParab
type(orbitals_data), intent(in) :: orbs
real(8):: parabPrefac
integer:: it
character(len=1):: direction
real(8):: shift
  !local variables
  integer :: idx,nf

     do idx=1,ncplx
        call applyLinearOperator(lr%d%n1,lr%d%n2,lr%d%n3,&
             lr%d%nfl1,lr%d%nfu1,lr%d%nfl2,lr%d%nfu2,lr%d%nfl3,lr%d%nfu3, &
             lr%wfd%nseg_c,lr%wfd%nvctr_c,lr%wfd%keyg,lr%wfd%keyv,&
             lr%wfd%nseg_f,lr%wfd%nvctr_f,&
             lr%wfd%keyg(1,lr%wfd%nseg_c+min(1,lr%wfd%nseg_f)),&
             lr%wfd%keyv(lr%wfd%nseg_c+min(1,lr%wfd%nseg_f)), &
             scal,hx,&
             lr%bounds%kb%ibyz_c,lr%bounds%kb%ibxz_c,lr%bounds%kb%ibxy_c,&
             lr%bounds%kb%ibyz_f,lr%bounds%kb%ibxz_f,lr%bounds%kb%ibxy_f,&
             x(1,idx),x(lr%wfd%nvctr_c+min(1,lr%wfd%nvctr_f),idx),&
             y(1,idx),y(lr%wfd%nvctr_c+min(1,lr%wfd%nvctr_f),idx),&
             w%xpsig_c,w%xpsig_f,w%ypsig_c,w%ypsig_f,&
             w%x_f1,w%x_f2,w%x_f3, rxyzParab, orbs, lr, parabPrefac, it, direction, shift)
     end do
END SUBROUTINE callApplyLinearOperator

! ypsi = (1/2) \Nabla^2 xpsi + cprecr xpsi
subroutine applyLinearOperator(n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3, & 
     nseg_c,nvctr_c,keyg_c,keyv_c,nseg_f,nvctr_f,keyg_f,keyv_f, &
     scal,hgrid,ibyz_c,ibxz_c,ibxy_c,ibyz_f,ibxz_f,ibxy_f,&
     xpsi_c,xpsi_f,ypsi_c,ypsi_f,&
     xpsig_c,xpsig_f,ypsig_c,ypsig_f,x_f1,x_f2,x_f3, rxyzParab, orbs, lr, parabPrefac, it, direction, shift)
  use module_base
  use module_types
  implicit none
  integer, intent(in) :: n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3
  integer, intent(in) :: nseg_c,nvctr_c,nseg_f,nvctr_f
  real(gp), intent(in) :: hgrid
  integer, dimension(nseg_c), intent(in) :: keyv_c
  integer, dimension(nseg_f), intent(in) :: keyv_f
  integer, dimension(2,nseg_c), intent(in) :: keyg_c
  integer, dimension(2,nseg_f), intent(in) :: keyg_f
  integer, dimension(2,0:n2,0:n3), intent(in) :: ibyz_c,ibyz_f
  integer, dimension(2,0:n1,0:n3), intent(in) :: ibxz_c,ibxz_f
  integer, dimension(2,0:n1,0:n2), intent(in) :: ibxy_c,ibxy_f
  real(wp), dimension(0:3), intent(in) :: scal
  real(wp), dimension(nvctr_c), intent(in) :: xpsi_c
  real(wp), dimension(7,nvctr_f), intent(in) :: xpsi_f
  real(wp), dimension(nvctr_c), intent(out) :: ypsi_c
  real(wp), dimension(7,nvctr_f), intent(out) :: ypsi_f
  real(wp), dimension(0:n1,0:n2,0:n3), intent(inout) :: xpsig_c,ypsig_c
  real(wp), dimension(7,nfl1:nfu1,nfl2:nfu2,nfl3:nfu3), intent(inout) :: xpsig_f,ypsig_f
  real(wp), dimension(nfl1:nfu1,nfl2:nfu2,nfl3:nfu3), intent(inout) :: x_f1
  real(wp), dimension(nfl2:nfu2,nfl1:nfu1,nfl3:nfu3), intent(inout) :: x_f2
  real(wp), dimension(nfl3:nfu3,nfl1:nfu1,nfl2:nfu2), intent(inout) :: x_f3
real(8),dimension(3),intent(in):: rxyzParab
type(orbitals_data), intent(in) :: orbs
type(locreg_descriptors), intent(in) :: lr
real(8):: parabPrefac
integer:: it
character(len=1):: direction
real(8):: shift

! Local variables
integer:: ix, iy, iz, ix0, iy0, iz0, i1, i2, i3, istat, ii, jj
real(8):: hxh, hyh, hzh, tt, dis, kx, ky, kz
real(8),dimension(:,:,:),allocatable:: xpsig_cTemp
real(8),dimension(:,:,:,:),allocatable:: xpsig_fTemp
real(8),dimension(:),allocatable:: phi, phir
type(workarr_locham) :: w2

real(8):: dnrm2


!write(*,*) 'dnrm2(xpsi_c)', dnrm2(nvctr_c, xpsi_c, 1)
!write(*,*) 'dnrm2(xpsi_f)', dnrm2(7*nvctr_f, xpsi_f, 1)

xpsig_c=0.d0
xpsig_f=0.d0

  call uncompress_forstandard(n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,  &
       nseg_c,nvctr_c,keyg_c,keyv_c,  & 
       nseg_f,nvctr_f,keyg_f,keyv_f,  & 
       scal,xpsi_c,xpsi_f,xpsig_c,xpsig_f,x_f1,x_f2,x_f3)

!write(*,*) 'dnrm2(xpsig_c)', dnrm2((n1+1)*(n2+1)*(n3+1), xpsig_c, 1)

  call ConvolLinear(n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3, & 
       hgrid,ibyz_c,ibxz_c,ibxy_c,ibyz_f,ibxz_f,ibxy_f,xpsig_c,&
       xpsig_f,ypsig_c,ypsig_f,x_f1,x_f2,x_f3, rxyzParab(1), parabPrefac, it, direction, shift)
!write(*,*) 'dnrm2(ypsig_c)', dnrm2((n1+1)*(n2+1)*(n3+1), ypsig_c, 1)
!write(*,*) 'parabPrefac',parabPrefac

  call compress_forstandard(n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,  &
       nseg_c,nvctr_c,keyg_c,keyv_c,  & 
       nseg_f,nvctr_f,keyg_f,keyv_f,  & 
       scal,ypsig_c,ypsig_f,ypsi_c,ypsi_f)

!write(*,*) 'dnrm2(ypsi_c)', dnrm2(nvctr_c, ypsi_c, 1)


END SUBROUTINE applyLinearOperator




subroutine estimatePerturbedOrbitals(iproc, nproc, at, orbs, lr, input, orbsLIN, commsLIN, rxyz,&
     nspin, nlpspd, proj, nscatterarr, ngatherarr, rhopot, GPU, pkernelseq, phi, rxyzParabola, perturbation)
!
! Purpose:
! ========
!   Calculates the localized basis functions phi. These basis functions are eigenfunctions of the ordinary Hamiltonian
!   with an additional parabolic potential centered at the atoms. The eigenfunctions are determined by minimizing the trace.
!
! Calling arguments:
!   Input arguments
!   Output arguments
!    phi   the localized basis functions
!
use module_base
use module_types
use module_interfaces, except_this_one => estimatePerturbedOrbitals
implicit none

! Calling arguments
integer:: iproc, nproc
type(atoms_data), intent(in) :: at
type(orbitals_data):: orbs
type(locreg_descriptors), intent(in) :: lr
type(input_variables):: input
type(orbitals_data):: orbsLIN
type(communications_arrays):: commsLIN
real(8),dimension(3,at%nat):: rxyz, rxyzParabola
integer:: nspin
type(nonlocal_psp_descriptors), intent(in) :: nlpspd
real(wp), dimension(nlpspd%nprojel), intent(in) :: proj
integer, dimension(0:nproc-1,4), intent(in) :: nscatterarr !n3d,n3p,i3s+i3xcsh-1,i3xcsh
integer, dimension(0:nproc-1,2), intent(in) :: ngatherarr 
real(dp), dimension(*), intent(inout) :: rhopot
type(GPU_pointers), intent(inout) :: GPU
real(dp), dimension(:), pointer :: pkernelseq
real(8),dimension(orbsLIN%npsidim):: phi
real(8),dimension(3,at%nat):: perturbation
! Local variables
integer:: iorb, jorb, korb, lorb, istat, istart, nvctrp, ierr, iat, ii, jj, kk, ik, il, ll, jproc, kproc, lproc, id, norbOnAt, info
integer:: jstart, lwork, ncplx, it, iiAt, idir
real(8),dimension(:),allocatable:: b, a, dvec, eval, work, phiw, hphiw, phi1, phiw2, phiPerturbed
real(8),dimension(:,:),allocatable:: eps, dTemp, d
real(8),dimension(:,:),allocatable:: emat, U, HU, HtildeSmall
real(8),dimension(:,:,:),allocatable:: Htilde
real(8):: tt, hphiNrm, dnrm2, fracTot, parabPrefac, ddot, ekin_sum, epot_sum, eexctx, eproj_sum
integer,dimension(:,:),allocatable:: orbsOnAtom
integer,dimension(:,:,:),allocatable:: tempArr
integer,dimension(:),allocatable:: lorbArr, displs, ipiv
real(8),dimension(:),pointer:: phiWorkPointer
character(len=1):: direction
real(8):: shift
real(8),dimension(:),allocatable:: hphi
!real(8),dimension(:,:),allocatable:: rxyzParabola


! Calculate the unconstrained gradient.
allocate(hphi(orbsLIN%npsidim), stat=istat)
call HamiltonianApplicationParabola(iproc,nproc,at,orbsLIN,input%hx,input%hy,input%hz,rxyz,&
     nlpspd,proj,lr,ngatherarr,lr%d%n1i*lr%d%n2i*nscatterarr(iproc,2),&
     rhopot(1),&
     phi(1),hphi(1),ekin_sum,epot_sum,eexctX,eproj_sum,nspin,GPU, orbsLIN%onWhichAtom, rxyzParabola, pkernel=pkernelseq)

allocate(orbsOnAtom(orbsLIN%norb,2), stat=istat)
allocate(lorbArr(0:nproc-1), stat=istat)
allocate(displs(0:nproc-1), stat=istat)
allocate(phiPerturbed(size(phi)), stat=istat)

! First copy phi to phiPerturbed
allocate(phiWorkPointer(size(phi)), stat=istat)
call transpose_v(iproc, nproc, orbsLIN, lr%wfd, commsLIN, phi, work=phiWorkPointer)
phiPerturbed=phi
call untranspose_v(iproc, nproc, orbsLIN, lr%wfd, commsLIN, phi, work=phiWorkPointer)
deallocate(phiWorkPointer, stat=istat)

atomsLoop: do iat=1,at%nat
    orbsOnAtom=0
    nvctrp=sum(commsLIN%nvctr_par(iproc,1:orbsLIN%nkptsp))*orbsLIN%nspinor


    ! First find out which atoms are centered on atom iat
    korb=0
    lorb=0
    do jproc=0,nproc-1
        do jorb=1,orbsLIN%norb_par(jproc)
            korb=korb+1
            if(jproc==iproc) then
                if(orbsLIN%onWhichAtom(jorb)==iat) then
                    lorb=lorb+1
                    orbsOnAtom(lorb,2)=korb
                end if
            end if
        end do
    end do
    call mpi_gather(lorb, 1, mpi_integer, lorbArr(0), 1, mpi_integer, 0, mpi_comm_world, ierr)
    displs(0)=0
    do jproc=1,nproc-1
        displs(jproc)=displs(jproc-1)+lorbArr(jproc-1)
    end do
    call mpi_gatherv(orbsOnAtom(1,2), lorb, mpi_integer, orbsOnAtom(1,1), lorbArr(0), displs(0), &
        mpi_integer, 0, mpi_comm_world, ierr)
    do istat=1,sum(lorbArr)
        if(iproc==0) write(*,*)'istat, orbsOnAtom(istat,1)', istat, orbsOnAtom(istat,1)
    end do

    ! Send orbsOnAtom to all processes
    if(iproc==0) norbOnAt=sum(lorbArr)
    call mpi_bcast(norbOnAt, 1, mpi_integer, 0, mpi_comm_world, ierr)
    call mpi_barrier(mpi_comm_world, ierr)
    call mpi_bcast(orbsOnAtom(1,1), norbOnat, mpi_integer, 0, mpi_comm_world, ierr)
    call mpi_barrier(mpi_comm_world, ierr)
 
    ! Calculate <phi|H|phi> for all orbitals. This requires to tranpose them first.
    allocate(Htilde(orbsLIN%norb,orbsLIN%norb,2), stat=istat)
    Htilde=0.d0
    allocate(phiWorkPointer(size(phi)), stat=istat)
    call transpose_v(iproc, nproc, orbsLIN, lr%wfd, commsLIN, phi, work=phiWorkPointer)
    call transpose_v(iproc, nproc, orbsLIN, lr%wfd, commsLIN, hphi, work=phiWorkPointer)
    deallocate(phiWorkPointer, stat=istat)
    nvctrp=commsLIN%nvctr_par(iproc,1)
    istart=1
    do iorb=1,orbsLIN%norb
        jstart=1
        do jorb=1,orbsLIN%norb
            Htilde(iorb,jorb,2)=ddot(nvctrp, phi(istart), 1, hphi(jstart), 1)
            jstart=jstart+nvctrp
        end do
        istart=istart+nvctrp
    end do
    call mpi_allreduce (Htilde(1,1,2), Htilde(1,1,1), orbsLIN%norb**2, mpi_double_precision, mpi_sum, mpi_comm_world, ierr)
    !if(iproc==0) write(*,*) '<phi|H|phi>'
    !do iorb=1,orbsLIN%norb
    !    if(iproc==0) write(*,'(100f8.4)') (Htilde(iorb,jorb,1), jorb=1,orbsLIN%norb)
    !end do

    ! Cut out the part of the matrix containing all orbitals centered on the current atom
    allocate(HtildeSmall(norbOnAt,norbOnAt), stat=istat)
    if(iproc==0) then
        do iorb=1,norbOnAt
            do jorb=1,norbOnAt
                HtildeSmall(jorb,iorb)=Htilde(orbsOnAtom(jorb,1),orbsOnAtom(iorb,1),1)
            end do
        end do 
        !do iorb=1,norbOnAt
        !    write(*,'(100f8.4)') (HtildeSmall(iorb,jorb), jorb=1,norbOnAt)
        !end do
        
        ! Diagonalize Htilde on root only
        allocate(eval(norbOnAt), stat=istat)
        lwork=1000*norbOnAt
        allocate(work(lwork), stat=istat)
        call dsyev('v', 'l', norbOnAt, HtildeSmall, norbOnAt, eval, work, lwork, info)
        if(info/=0) then
            write(*,'(a,i0)') 'ERROR in dsyev, info= ',info
            stop
        end if
        deallocate(work, stat=istat)
        !write(*,*) '--------------'
        !do iorb=1,norbOnAt
        !    write(*,'(100f8.4)') (HtildeSmall(iorb,jorb), jorb=1,norbOnAt)
        !end do
    end if

    ! Broadcast HtildeSmall
    call mpi_bcast(HtildeSmall(1,1), norbOnAt**2, mpi_double_precision, 0, mpi_comm_world, ierr)
    call mpi_barrier(mpi_comm_world, ierr)
    !write(*,*) '--------------', iproc, norbOnAt
    !do iorb=1,norbOnAt
    !    write(*,'(100f8.4)') (HtildeSmall(iorb,jorb), jorb=1,norbOnAt)
    !end do

    ! Build new linear combination
    allocate(phiw(orbsLIN%npsidim), stat=istat)
    allocate(hphiw(orbsLIN%npsidim), stat=istat)
    phiw=0.d0 ; hphiw=0.d0
    nvctrp=commsLIN%nvctr_par(iproc,1)
    do iorb=1,norbOnAt
        istart=(orbsOnAtom(iorb,1)-1)*nvctrp+1
        do jorb=1,norbOnAt
            jstart=(orbsOnAtom(jorb,1)-1)*nvctrp+1
            call daxpy(nvctrp, HtildeSmall(jorb,iorb), phi(jstart), 1, phiw(istart), 1)
        end do
    end do

    !!! Undo this transformations only a test)
    !!allocate(phiw2(size(phiw)), stat=istat)
    !!phiw2=phiw
    !!phiw=0.d0
    !!do iorb=1,norbOnAt
    !!    istart=(orbsOnAtom(iorb,1)-1)*nvctrp+1
    !!    do jorb=1,norbOnAt
    !!        jstart=(orbsOnAtom(jorb,1)-1)*nvctrp+1
    !!        call daxpy(nvctrp, HtildeSmall(iorb,jorb), phiw2(jstart), 1, phiw(istart), 1)
    !!    end do
    !!end do

    allocate(phiWorkPointer(size(phi)), stat=istat)
    call untranspose_v(iproc, nproc, orbsLIN, lr%wfd, commsLIN, phi, work=phiWorkPointer)
    call untranspose_v(iproc, nproc, orbsLIN, lr%wfd, commsLIN, phiw, work=phiWorkPointer)
    call untranspose_v(iproc, nproc, orbsLIN, lr%wfd, commsLIN, hphi, work=phiWorkPointer)
    deallocate(phiWorkPointer, stat=istat)

    ! Check whether the subspace diagonalization was successful
    ! This should at the moment always be done to get the diagonal elements
    if(.true.) then
        nspin=1
        !call HamiltonianApplicationParabola(iproc,nproc,at,orbsLIN,input%hx,input%hy,input%hz,rxyz,&
        !     nlpspd,proj,lr,ngatherarr,lr%d%n1i*lr%d%n2i*nscatterarr(iproc,2),&
        !     rhopot(1),&
        !     phiw(1),hphiw(1),ekin_sum,epot_sum,eexctX,eproj_sum,nspin,GPU, orbsLIN%onWhichAtom, rxyzParabola, pkernel=pkernelseq)
        call HamiltonianApplicationParabola(iproc,nproc,at,orbsLIN,input%hx,input%hy,input%hz,rxyz,&
            nlpspd,proj,lr,ngatherarr,lr%d%n1i*lr%d%n2i*nscatterarr(iproc,2),&
            rhopot(1),&
            phiw(1),hphiw(1),ekin_sum,epot_sum,eexctX,eproj_sum,nspin,GPU, orbsLIN%onWhichAtom, rxyzParabola, pkernel=pkernelseq)

        Htilde=0.d0
        allocate(phiWorkPointer(size(phi)), stat=istat)
        call transpose_v(iproc, nproc, orbsLIN, lr%wfd, commsLIN, phiw, work=phiWorkPointer)
        call transpose_v(iproc, nproc, orbsLIN, lr%wfd, commsLIN, hphiw, work=phiWorkPointer)
        deallocate(phiWorkPointer, stat=istat)
        nvctrp=commsLIN%nvctr_par(iproc,1)
        istart=1
        do iorb=1,orbsLIN%norb
            jstart=1
            do jorb=1,orbsLIN%norb
                Htilde(iorb,jorb,2)=ddot(nvctrp, phiw(istart), 1, hphiw(jstart), 1)
                jstart=jstart+nvctrp
            end do
            istart=istart+nvctrp
        end do
        call mpi_allreduce (Htilde(1,1,2), Htilde(1,1,1), orbsLIN%norb**2 ,mpi_double_precision, mpi_sum, mpi_comm_world, ierr)
        !if(iproc==0) write(*,*) 'after subspace diag: <phi|H|phi>', iat
        !do iorb=1,orbsLIN%norb
        !    if(iproc==0) write(*,'(100f8.4)') (Htilde(iorb,jorb,1), jorb=1,orbsLIN%norb)
        !end do
    end if       

     if(.not.allocated(phiw2)) allocate(phiw2(size(phiw)), stat=istat)
     !phiw2=0.d0

    ! After this operation phiw2 contains the transformed phi in the transposed way.
    phiw2=phiw  ! here it is transposed?
    directionLoop: do idir=1,3
           if(idir==1) then
               direction='x'
           else if(idir==2) then
               direction='y'
           else if(idir==3) then
               direction='z'
           end if
           ! Calculate the matrix elements <phiw|V|phiw>
           ncplx=1
           it=1
           ! First apply the perturbation to all orbitals belonging to iproc
           allocate(phiWorkPointer(size(phi)), stat=istat)
           call untranspose_v(iproc, nproc, orbsLIN, lr%wfd, commsLIN, phiw, work=phiWorkPointer)
           hphiw=phiw
           istart=1
           shift=perturbation(idir,iat)
!write(*,*) 'ATTENTION: CHANGING SHIFT!!'
!SHIFT=1.d2
           do iorb=1,orbsLIN%norbp
               iiAt=orbsLIN%onWhichAtom(iorb)
               parabPrefac=orbsLIN%parabPrefacArr(at%iatype(orbsLIN%onWhichAtom(iorb)))
               !write(*,'(a,2i6,es12.4)') 'before: iproc, iorb, dnrm2(hphiw)', iproc, iorb, dnrm2(nvctrp,hphiw(istart),1)
               call subTemp(lr,ncplx,&
                    input%hx,input%hy,input%hz,hphiw(istart), rxyzParabola(3,iat), orbsLIN, parabPrefac, it, direction, shift)
               !write(*,'(a,2i6,es12.4)') 'after: iproc, iorb, dnrm2(hphiw)', iproc, iorb, dnrm2(nvctrp,hphiw(istart),1)
               istart=istart+(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f)*orbsLIN%nspinor
           end do
    
           call transpose_v(iproc, nproc, orbsLIN, lr%wfd, commsLIN, phiw, work=phiWorkPointer)
           call transpose_v(iproc, nproc, orbsLIN, lr%wfd, commsLIN, hphiw, work=phiWorkPointer)
    
    
           ! nvctrp is the amount of each phi hold by the current process
           nvctrp=sum(commsLIN%nvctr_par(iproc,1:orbsLIN%nkptsp))*orbsLIN%nspinor
           
           allocate(dTemp(orbsLIN%norb,orbsLIN%norb), stat=istat)
           allocate(d(orbsLIN%norb,orbsLIN%norb), stat=istat)
           dTemp=0.d0
           istart=1
           do iorb=1,orbsLIN%norb
               jstart=1
               do jorb=1,orbsLIN%norb
                   dTemp(iorb,jorb)=ddot(nvctrp, phiw(istart), 1, hphiw(jstart), 1)
                   !if(iproc==0) write(*,'(a,2i6,2es12.4)') 'iorb, jorb, dnrm2(phiw), dnrm2(hphiw)', iorb, jorb, dnrm2(nvctrp,phiw(istart),1), dnrm2(nvctrp,hphiw(jstart),1)
                   jstart=jstart+nvctrp
               end do
               istart=istart+nvctrp
           end do
           d=0.d0
           call mpi_allreduce(dTemp(1,1), d(1,1), orbsLIN%norb**2, mpi_double_precision, mpi_sum, mpi_comm_world, ierr)
           ! Add the diagonal element shift**2
           do iorb=1,orbsLIN%norb
               d(iorb,iorb)=d(iorb,iorb)+parabPrefac*dble(shift)**2
           end do
           
           !if(iproc==0) write(*,'(a,es12.4)') '<phi|V|phi> with shift=',shift
           !do iorb=1,orbsLIN%norb
           !    if(iproc==0) write(*,'(100es15.8)') (d(iorb,jorb), jorb=1,orbsLIN%norb)
           !end do
    
           ! Make a linear combination according to perturbation theory
           !!do iorb=1,orbsLIN%norbp
           !!    write(*,'(a,i5,i3,i5,es12.5)') 'iat, iproc, iorb, dnrm2(phiw)', iat, iproc, iorb, dnrm2(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f, phiw(istart), 1)
           !!    istart=istart+lr%wfd%nvctr_c+7*lr%wfd%nvctr_f
           !!end do
           allocate(phi1(orbsLIN%npsidim), stat=istat)
           phi1=0.d0
           do iorb=1,norbOnAt
               !istart=(iorb-1)*nvctrp+1
               istart=(orbsOnAtom(iorb,1)-1)*nvctrp+1
               do jorb=1,norbOnAt
                   jstart=(orbsOnAtom(jorb,1)-1)*nvctrp+1
                   if(iorb==jorb) cycle
                   call daxpy(nvctrp, d(orbsOnATom(iorb,1),orbsOnAtom(jorb,1))/&
                   (Htilde(orbsOnAtom(iorb,1),orbsOnAtom(iorb,1),1)-Htilde(orbsOnAtom(jorb,1),orbsOnAtom(jorb,1),1)), &
                   phiw(jstart), 1, phi1(istart), 1)
               end do
           end do
           !call untranspose_v(iproc, nproc, orbsLIN, lr%wfd, commsLIN, phi1, work=phiWorkPointer)
           !istart=1
           !do iorb=1,orbsLIN%norbp
           !    write(*,'(a,i3,i5,2es12.5)') 'iproc, iorb, dnrm2(phi1), dnrm2(phi)', iproc, iorb, dnrm2(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f, phi1(istart), 1), dnrm2(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f, phi(istart), 1)
           !    istart=istart+lr%wfd%nvctr_c+7*lr%wfd%nvctr_f
           !end do
           !call transpose_v(iproc, nproc, orbsLIN, lr%wfd, commsLIN, phi1, work=phiWorkPointer)

           ! Add the correction phi1 to the old orbitals phiw2
           istart=1
           do iorb=1,norbOnAt
               jstart=(orbsOnAtom(iorb,1)-1)*nvctrp+1
               call daxpy(nvctrp, 1.d0, phi1(jstart), 1, phiw2(jstart), 1)
               !call daxpy(nvctrp, 1.d0, phi1(jstart), 1, phiw(jstart), 1)
               !call daxpy(nvctrp, 1.d0, phi1(istart), 1, phiw(jstart), 1)
               istart=istart+nvctrp
        
           end do
        ! Try to stop the loop here
        end do directionLoop

        ! Now phiw2 contains the perturbed orbitals.

           !call transpose_v(iproc, nproc, orbsLIN, lr%wfd, commsLIN, phi, work=phiWorkPointer)
           ! Transform back to 'old' orbitals
           !if(.not.allocated(phiw2)) allocate(phiw2(size(phiw)), stat=istat)
           !phiw2=phiw
           phiw=0.d0
           do iorb=1,norbOnAt
               istart=(orbsOnAtom(iorb,1)-1)*nvctrp+1
               do jorb=1,norbOnAt
                   jstart=(orbsOnAtom(jorb,1)-1)*nvctrp+1
                   call daxpy(nvctrp, HtildeSmall(iorb,jorb), phiw2(jstart), 1, phiw(istart), 1)
               end do
           end do

           !! Now mix the orbitals
           !phiw2=phi
           !istart=1
           !do iorb=1,norbOnAt
           !    jstart=(orbsOnAtom(iorb,1)-1)*nvctrp+1
           !    !call dcopy(nvctrp, phiw(istart), 1, phiw2(jstart), 1)
           !    call dcopy(nvctrp, phiw(jstart), 1, phiw2(jstart), 1)
           !    istart=istart+nvctrp
           !end do
           ! Mix the orbitals in phiPerturbed. phiPerturbed was initialized to be =phi in the very beginning,
           ! So we can just replace the perturbed orbitals without doing anything else
           istart=1
           do iorb=1,norbOnAt
               jstart=(orbsOnAtom(iorb,1)-1)*nvctrp+1
               !call dcopy(nvctrp, phiw(istart), 1, phiw2(jstart), 1)
               call dcopy(nvctrp, phiw(jstart), 1, phiPerturbed(jstart), 1)
               istart=istart+nvctrp
           end do

           ! This is not necessary, only for debugging
           !call untranspose_v(iproc, nproc, orbsLIN, lr%wfd, commsLIN, phiw2, work=phiWorkPointer)
           !istart=1
           !do iorb=1,orbsLIN%norbp
           !    write(*,'(a,i4,i3,i5,2es12.5)') 'iat, iproc, iorb, dnrm2(phiw2)', iat, iproc, iorb, dnrm2(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f, phiw2(istart), 1)
           !    istart=istart+lr%wfd%nvctr_c+7*lr%wfd%nvctr_f
           !end do
           !call transpose_v(iproc, nproc, orbsLIN, lr%wfd, commsLIN, phiw2, work=phiWorkPointer)


           ! Now orthonormalize
           call orthogonalize(iproc, nproc, orbsLIN, commsLIN, lr%wfd, phiw2, input)
    
           !call untranspose_v(iproc, nproc, orbsLIN, lr%wfd, commsLIN, phi, work=phiWorkPointer)
           call untranspose_v(iproc, nproc, orbsLIN, lr%wfd, commsLIN, phiw, work=phiWorkPointer)
           call untranspose_v(iproc, nproc, orbsLIN, lr%wfd, commsLIN, hphiw, work=phiWorkPointer)
           call untranspose_v(iproc, nproc, orbsLIN, lr%wfd, commsLIN, phiw2, work=phiWorkPointer)
           deallocate(phiWorkPointer, stat=istat)

           ! Now phiw2 contains the perturbed basis functions.
           ! Copy them back
           !istart=1
           !do iorb=1,orbsLIN%norbp
           !    write(*,'(a,i4,i3,i5,2es24.15)') 'iat, iproc, iorb, <phi|phiw2>', iat, iproc, iorb, ddot(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f, phi(istart), 1, phiw2(istart), 1)
           !    istart=istart+lr%wfd%nvctr_c+7*lr%wfd%nvctr_f
           !end do

           ! DO NOT COPY HERE
           !phi=phiw2
    
        !call mpi_barrier(mpi_comm_world, ierr)
        !stop
    !end do directionLoop
end do atomsLoop

! Replace phi with the perturbed phiPerturbed
call orthogonalize(iproc, nproc, orbsLIN, commsLIN, lr%wfd, phiPerturbed, input)
allocate(phiWorkPointer(size(phi)), stat=istat)
call untranspose_v(iproc, nproc, orbsLIN, lr%wfd, commsLIN, phiPerturbed, work=phiWorkPointer)
deallocate(phiWorkPointer, stat=istat)
write(*,*) 'check', ddot(size(phi), phi(1), 1, phiPerturbed(1), 1)/&
    (dnrm2(size(phi), phi(1), 1)*dnrm2(size(phiPerturbed), phiPerturbed(1), 1))
phi=phiPerturbed

end subroutine estimatePerturbedOrbitals
