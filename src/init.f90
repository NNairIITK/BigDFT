!calculates the descriptor arrays
!calculates also the bounds arrays needed for convolutions
!refers this information to the global localisation region descriptor
subroutine createWavefunctionsDescriptors(iproc,hx,hy,hz,atoms,rxyz,radii_cf,&
     crmult,frmult,Glr,orbs)
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
  type(orbitals_data), intent(inout) :: orbs
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
        if (iproc ==0)then
           write(*,*)&
                '          The code should be stopped for a GPU calculation     '
           write(*,*)&
                '          since density is not initialised to 10^-20               '
        end if
        stop
     end if
  end if

  ! fine grid quantities
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
  call segkeys(n1,n2,n3,0,n1,0,n2,0,n3,logrid_f,Glr%wfd%nseg_f,Glr%wfd%keyg(1,Glr%wfd%nseg_c+1), &
       & Glr%wfd%keyv(Glr%wfd%nseg_c+1))

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

end subroutine createWavefunctionsDescriptors

!pass to implicit none while inserting types on this routine
subroutine createProjectorsArrays(iproc,n1,n2,n3,rxyz,at,&
     radii_cf,cpmult,fpmult,hx,hy,hz,nlpspd,proj)
  use module_base
  use module_types
  implicit none
  type(atoms_data), intent(in) :: at
  integer, intent(in) :: iproc,n1,n2,n3
  real(gp), intent(in) :: cpmult,fpmult,hx,hy,hz
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
       logrid,at,nlpspd)

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
        call segkeys(n1,n2,n3,nl1,nu1,nl2,nu2,nl3,nu3,  & 
             logrid,mseg,nlpspd%keyg_p(1,iseg),nlpspd%keyv_p(iseg))

     endif
  enddo

  i_all=-product(shape(logrid))*kind(logrid)
  deallocate(logrid,stat=i_stat)
  call memocc(i_stat,i_all,'logrid',subname)

  !fill the projectors if the strategy is a distributed calculation
  if (.not. DistProjApply) then
     !calculate the wavelet expansion of projectors
     call fill_projectors(iproc,n1,n2,n3,hx,hy,hz,cpmult,fpmult,at,rxyz,radii_cf,&
          nlpspd,proj,0)
  end if

END SUBROUTINE createProjectorsArrays

subroutine import_gaussians(iproc,nproc,cpmult,fpmult,radii_cf,at,orbs,comms,&
     Glr,hx,hy,hz,rxyz,rhopot,pot_ion,nlpspd,proj,& 
     pkernel,ixc,psi,psit,hpsi,nscatterarr,ngatherarr,nspin)
  use module_base
  use module_interfaces, except_this_one => import_gaussians
  use module_types
  use Poisson_Solver
  implicit none
  integer, intent(in) :: iproc,nproc,ixc,nspin
  real(gp), intent(in) :: hx,hy,hz,cpmult,fpmult
  type(atoms_data), intent(in) :: at
  type(nonlocal_psp_descriptors), intent(in) :: nlpspd
  type(locreg_descriptors), intent(in) :: Glr
  type(communications_arrays), intent(in) :: comms
  integer, dimension(0:nproc-1,4), intent(in) :: nscatterarr !n3d,n3p,i3s+i3xcsh-1,i3xcsh
  integer, dimension(0:nproc-1,2), intent(in) :: ngatherarr 
  real(gp), dimension(at%ntypes,3), intent(in) :: radii_cf  
  real(gp), dimension(3,at%nat), intent(in) :: rxyz
  real(wp), dimension(nlpspd%nprojel), intent(in) :: proj
  real(dp), dimension(*), intent(in) :: pkernel
  type(orbitals_data), intent(inout) :: orbs
  real(dp), dimension(*), intent(inout) :: rhopot
  real(wp), dimension(*), intent(inout) :: pot_ion
  real(wp), dimension(:), pointer :: psi,psit,hpsi
  !local variables
  character(len=*), parameter :: subname='import_gaussians'
  integer :: i_stat,i_all
  real(gp) :: hxh,hyh,hzh,eexcu,vexcu,epot_sum,ekin_sum,ehart,eproj_sum
  type(gaussian_basis) :: CP2K
  type(GPU_pointers) :: GPU !added for interface compatibility, not working here
  real(wp), dimension(:,:), pointer :: wfn_cp2k

  if (iproc.eq.0) then
     write(*,'(1x,a)')&
          '--------------------------------------------------------- Import Gaussians from CP2K'
  end if

  if (nspin /= 1) then
     if (iproc==0) then
        write(*,'(1x,a)')&
             'Gaussian importing is possible only for non-spin polarised calculations'
        write(*,'(1x,a)')&
             'The writing rules of CP2K files for spin-polarised orbitals are not implemented'
     end if
     stop
  end if

  hxh=.5_gp*hx
  hyh=.5_gp*hy
  hzh=.5_gp*hz

  call parse_cp2k_files(iproc,'gaubasis.dat','gaucoeff.dat',&
       at%nat,at%ntypes,orbs,at%iatype,rxyz,CP2K,wfn_cp2k)

  !allocate the wavefunction in the transposed way to avoid allocations/deallocations
  allocate(psi(orbs%npsidim+ndebug),stat=i_stat)
  call memocc(i_stat,psi,'psi',subname)


!!$  !put to zero the value of psi function to control what happens for a complex case
!!$  call razero(orbs%npsidim,psi)

  !calculate the overlap matrix for debugging and testing
!!$  allocate(ovrlp(CP2K%ncoeff*CP2K%ncoeff),stat=i_stat)
!!$  call memocc(i_stat,ovrlp,'ovrlp',subname)
!!$  allocate(tmp(CP2K%ncoeff,norb),stat=i_stat)
!!$  call memocc(i_stat,tmp,'tmp',subname)
!!$  allocate(smat(norb,norb),stat=i_stat)
!!$  call memocc(i_stat,smat,'smat',subname)
!!$  !overlap calculation of the gaussian matrix, to be done in view of quick restart
!!$  call gaussian_overlap(CP2K,CP2K,ovrlp)
!!$  call dsymm('L','U',CP2K%ncoeff,norb,1.0_gp,ovrlp(1),CP2K%ncoeff,wfn_cp2k(1,1),CP2K%ncoeff,&
!!$       0.d0,tmp(1,1),CP2K%ncoeff)
!!$  i_all=-product(shape(ovrlp))*kind(ovrlp)
!!$  deallocate(ovrlp,stat=i_stat)
!!$  call memocc(i_stat,i_all,'ovrlp',subname)
!!$  call gemm('T','N',norb,norb,CP2K%ncoeff,1.0_gp,wfn_cp2k(1,1),CP2K%ncoeff,tmp(1,1),CP2K%ncoeff,&
!!$       0.0_wp,smat(1,1),norb)
!!$  !print overlap matrices
!!$  do i=1,norb
!!$     write(*,'(i5,30(1pe19.12))')i,(smat(i,iorb),iorb=1,norb)
!!$  end do
!!$  i_all=-product(shape(tmp))*kind(tmp)
!!$  deallocate(tmp,stat=i_stat)
!!$  call memocc(i_stat,i_all,'tmp',subname)
!!$  i_all=-product(shape(smat))*kind(smat)
!!$  deallocate(smat,stat=i_stat)
!!$  call memocc(i_stat,i_all,'smat',subname)

  call gaussians_to_wavelets(iproc,nproc,at%geocode,orbs,Glr%d,&
       hx,hy,hz,Glr%wfd,CP2K,wfn_cp2k,psi)

  !deallocate CP2K variables
  call deallocate_gwf(CP2K,subname)
  !nullify gaussian centers
  nullify(CP2K%rxyz)
  i_all=-product(shape(wfn_cp2k))*kind(wfn_cp2k)
  deallocate(wfn_cp2k,stat=i_stat)
  call memocc(i_stat,i_all,'wfn_cp2k',subname)

  call sumrho(iproc,nproc,orbs,Glr,ixc,hxh,hyh,hzh,psi,rhopot,&
       Glr%d%n1i*Glr%d%n2i*nscatterarr(iproc,1),nscatterarr,1,GPU)

  call PSolver(at%geocode,'D',iproc,nproc,Glr%d%n1i,Glr%d%n2i,Glr%d%n3i,ixc,hxh,hyh,hzh,&
       rhopot,pkernel,pot_ion,ehart,eexcu,vexcu,0.0_dp,.true.,1)

  !allocate the wavefunction in the transposed way to avoid allocations/deallocations
  allocate(hpsi(orbs%npsidim+ndebug),stat=i_stat)
  call memocc(i_stat,hpsi,'hpsi',subname)

  call HamiltonianApplication(iproc,nproc,at,orbs,hx,hy,hz,rxyz,cpmult,fpmult,radii_cf,&
       nlpspd,proj,Glr,ngatherarr,&
       Glr%d%n1i*Glr%d%n2i*nscatterarr(iproc,2),&
       rhopot(1+Glr%d%n1i*Glr%d%n2i*nscatterarr(iproc,4)),&
       psi,hpsi,ekin_sum,epot_sum,eproj_sum,1,GPU)

  if (iproc == 0 .and. verbose > 1) write(*,'(1x,a,(f19.10))') 'done. ekin_sum',ekin_sum

  if (iproc == 0) then
     write(*,'(1x,a,3(1x,1pe18.11))') 'ekin_sum,epot_sum,eproj_sum',  & 
          ekin_sum,epot_sum,eproj_sum
     write(*,'(1x,a,3(1x,1pe18.11))') '   ehart,   eexcu,    vexcu',ehart,eexcu,vexcu
  endif

  !after having applied the hamiltonian to all the atomic orbitals
  !we split the semicore orbitals from the valence ones
  !this is possible since the semicore orbitals are the first in the 
  !order, so the linear algebra on the transposed wavefunctions 
  !may be splitted

  if (iproc == 0 .and. verbose > 1) write(*,'(1x,a)',advance='no')&
       'Imported Wavefunctions Orthogonalization:'

  call DiagHam(iproc,nproc,at%natsc,nspin,orbs,Glr%wfd,comms,psi,hpsi,psit)

  if (iproc == 0 .and. verbose > 1) write(*,'(1x,a)')'done.'

END SUBROUTINE import_gaussians

subroutine input_wf_diag(iproc,nproc,cpmult,fpmult,radii_cf,at,&
     orbs,orbsv,nvirt,comms,Glr,hx,hy,hz,rxyz,rhopot,pot_ion,&
     nlpspd,proj,pkernel,ixc,psi,hpsi,psit,psivirt,&
     nscatterarr,ngatherarr,nspin)
  ! Input wavefunctions are found by a diagonalization in a minimal basis set
  ! Each processors write its initial wavefunctions into the wavefunction file
  ! The files are then read by readwave
  use module_base
  use module_interfaces, except_this_one => input_wf_diag
  use module_types
  use Poisson_Solver
  implicit none
  integer, intent(in) :: iproc,nproc,ixc
  integer, intent(inout) :: nspin,nvirt
  real(gp), intent(in) :: hx,hy,hz,cpmult,fpmult
  type(atoms_data), intent(in) :: at
  type(orbitals_data), intent(inout) :: orbs
  type(nonlocal_psp_descriptors), intent(in) :: nlpspd
  type(locreg_descriptors), intent(in) :: Glr
  type(communications_arrays), intent(in) :: comms
  integer, dimension(0:nproc-1,4), intent(in) :: nscatterarr !n3d,n3p,i3s+i3xcsh-1,i3xcsh
  integer, dimension(0:nproc-1,2), intent(in) :: ngatherarr 
  real(gp), dimension(at%ntypes,3), intent(in) :: radii_cf  
  real(gp), dimension(3,at%nat), intent(in) :: rxyz
  real(wp), dimension(nlpspd%nprojel), intent(in) :: proj
  real(dp), dimension(*), intent(in) :: pkernel
  real(dp), dimension(*), intent(inout) :: rhopot,pot_ion
  type(orbitals_data), intent(out) :: orbsv
  real(wp), dimension(:), pointer :: psi,hpsi,psit,psivirt
  !local variables
  character(len=*), parameter :: subname='input_wf_diag'
  integer, parameter :: ngx=31
  integer :: i_stat,i_all,iat,nspin_ig
  real(gp) :: hxh,hyh,hzh,eks,eexcu,vexcu,epot_sum,ekin_sum,ehart,eproj_sum,etol,accurex
  type(gaussian_basis) :: G
  type(orbitals_data) :: orbse
  type(communications_arrays) :: commse
  type(GPU_pointers) :: GPU
  integer, dimension(:,:), allocatable :: norbsc_arr
  real(gp), dimension(:), allocatable :: locrad
  type(locreg_descriptors), dimension(:), allocatable :: Llr
  real(wp), dimension(:,:,:), pointer :: psigau

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
       orbs,orbse,orbsv,norbsc_arr,locrad,G,psigau,eks)

  !allocate communications arrays for inputguess orbitals
  call allocate_comms(nproc,commse,subname)
  call orbitals_communicators(iproc,nproc,Glr,orbse,commse)  

  i_all=-product(shape(orbse%norb_par))*kind(orbse%norb_par)
  deallocate(orbse%norb_par,stat=i_stat)
  call memocc(i_stat,i_all,'orbse%norb_par',subname)

  hxh=.5_gp*hx
  hyh=.5_gp*hy
  hzh=.5_gp*hz

  !check the communication distribution
  call check_communications(iproc,nproc,orbse,Glr,commse)

  !once the wavefunction coefficients are known perform a set 
  !of nonblocking send-receive operations to calculate overlap matrices

!!$  !create mpirequests array for controlling the success of the send-receive operation
!!$  allocate(mpirequests(nproc-1+ndebug),stat=i_stat)
!!$  call memocc(i_stat,mpirequests,'mpirequests',subname)
!!$
!!$  call nonblocking_transposition(iproc,nproc,G%ncoeff,orbse%isorb+orbse%norbp,&
!!$       orbse%nspinor,psigau,orbse%norb_par,mpirequests)


  !experimental part for building the localisation regions
  if (at%geocode == 'F') then
     !allocate the array of localisation regions
     allocate(Llr(at%nat+ndebug),stat=i_stat)
     !call memocc(i_stat,Llr,'Llr',subname)

     !print *,'locrad',locrad

     call determine_locreg(at%nat,rxyz,locrad,hx,hy,hz,Glr,Llr)

     do iat=1,at%nat
        call deallocate_wfd(Llr(iat)%wfd,subname)
        if (Llr(iat)%geocode=='F') then
           call deallocate_bounds(Llr(iat)%bounds,subname)
        end if
     end do

     !i_all=-product(shape(Llr))*kind(Llr)
     deallocate(Llr,stat=i_stat) !these allocation are special
     !call memocc(i_stat,i_all,'Llr',subname)
  end if

  !allocate the wavefunction in the transposed way to avoid allocations/deallocations
  allocate(psi(orbse%npsidim+ndebug),stat=i_stat)
  call memocc(i_stat,psi,'psi',subname)

  !allocate arrays for the GPU if a card is present
  if (GPUconv) then
       call prepare_gpu_for_locham(Glr%d%n1,Glr%d%n2,Glr%d%n3,nspin,&
            hx,hy,hz,Glr%wfd,orbse,GPU)
  end if

  !use only the part of the arrays for building the hamiltonian matrix
  call gaussians_to_wavelets(iproc,nproc,at%geocode,orbse,Glr%d,&
       hx,hy,hz,Glr%wfd,G,psigau(1,1,min(orbse%isorb+1,orbse%norb)),psi)

  i_all=-product(shape(locrad))*kind(locrad)
  deallocate(locrad,stat=i_stat)
  call memocc(i_stat,i_all,'locrad',subname)

  !application of the hamiltonian for gaussian based treatment
  call sumrho(iproc,nproc,orbse,Glr,ixc,hxh,hyh,hzh,psi,rhopot,&
       Glr%d%n1i*Glr%d%n2i*nscatterarr(iproc,1),nscatterarr,nspin,GPU)
  
  if(orbs%nspinor==4) then
     !this wrapper can be inserted inside the poisson solver 
     call PSolverNC(at%geocode,'D',iproc,nproc,Glr%d%n1i,Glr%d%n2i,Glr%d%n3i,&
          nscatterarr(iproc,1),& !this is n3d
          ixc,hxh,hyh,hzh,&
          rhopot,pkernel,pot_ion,ehart,eexcu,vexcu,0.d0,.true.,4)
  else
     call PSolver(at%geocode,'D',iproc,nproc,Glr%d%n1i,Glr%d%n2i,Glr%d%n3i,&
          ixc,hxh,hyh,hzh,&
          rhopot,pkernel,pot_ion,ehart,eexcu,vexcu,0.d0,.true.,nspin)
  end if


!!$  if (nproc == 1) then
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
!!$          gaucoeff(1,1),G%ncoeff,0.d0,tmp(1,1),G%ncoeff)
!!$
!!$     call gemm('T','N',orbse%norb,orbse%norb,G%ncoeff,1.0_gp,&
!!$          gaucoeff(1,1),G%ncoeff,tmp(1,1),G%ncoeff,0.0_wp,smat(1,1),orbse%norb)
!!$
!!$     !print overlap matrices
!!$     do i=1,orbse%norb
!!$        write(*,'(i5,30(1pe15.8))')i,(smat(i,iorb),iorb=1,orbse%norb)
!!$     end do
!!$
!!$     !overlap calculation of the kinetic operator
!!$     call kinetic_overlap(G,G,ovrlp)
!!$     call dsymm('L','U',G%ncoeff,orbse%norb,1.0_gp,ovrlp(1),G%ncoeff,&
!!$          gaucoeff(1,1),G%ncoeff,0.d0,tmp(1,1),G%ncoeff)
!!$
!!$     call gemm('T','N',orbse%norb,orbse%norb,G%ncoeff,1.0_gp,&
!!$          gaucoeff(1,1),G%ncoeff,tmp(1,1),G%ncoeff,0.0_wp,smat(1,1),orbse%norb)
!!$
!!$     !print overlap matrices
!!$     tt=0.0_wp
!!$     do i=1,orbse%norb
!!$        write(*,'(i5,30(1pe15.8))')i,(smat(i,iorb),iorb=1,orbse%norb)
!!$        !write(12,'(i5,30(1pe15.8))')i,(smat(i,iorb),iorb=1,orbse%norb)
!!$        tt=tt+smat(i,i)
!!$     end do
!!$     print *,'trace',tt
!!$
!!$     !overlap calculation of the kinetic operator
!!$     call cpu_time(t0)
!!$     call potential_overlap(G,G,rhopot,Glr%d%n1i,Glr%d%n2i,Glr%d%n3i,hxh,hyh,hzh,&
!!$          ovrlp)
!!$     call cpu_time(t1)
!!$     call dsymm('L','U',G%ncoeff,orbse%norb,1.0_gp,ovrlp(1),G%ncoeff,&
!!$          gaucoeff(1,1),G%ncoeff,0.d0,tmp(1,1),G%ncoeff)
!!$
!!$     call gemm('T','N',orbse%norb,orbse%norb,G%ncoeff,1.0_gp,&
!!$          gaucoeff(1,1),G%ncoeff,tmp(1,1),G%ncoeff,0.0_wp,smat(1,1),orbse%norb)
!!$
!!$     !print overlap matrices
!!$     tt=0.0_wp
!!$     do i=1,orbse%norb
!!$        write(*,'(i5,30(1pe15.8))')i,(smat(i,iorb),iorb=1,orbse%norb)
!!$        !write(12,'(i5,30(1pe15.8))')i,(smat(i,iorb),iorb=1,orbse%norb)
!!$        tt=tt+smat(i,i)
!!$     end do
!!$     print *,'trace',tt
!!$     print *, 'time',t1-t0
!!$
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
  allocate(hpsi(orbse%npsidim+ndebug),stat=i_stat)
  call memocc(i_stat,hpsi,'hpsi',subname)
  
  call HamiltonianApplication(iproc,nproc,at,orbse,hx,hy,hz,rxyz,cpmult,fpmult,radii_cf,&
       nlpspd,proj,Glr,ngatherarr,Glr%d%n1i*Glr%d%n2i*nscatterarr(iproc,2),&
       rhopot(1+Glr%d%n1i*Glr%d%n2i*nscatterarr(iproc,4)),&
       psi,hpsi,ekin_sum,epot_sum,eproj_sum,nspin,GPU)

!!$  !calculate the overlap matrix knowing that the original functions are gaussian-based
!!$  allocate(thetaphi(2,G%nat+ndebug),stat=i_stat)
!!$  call memocc(i_stat,thetaphi,'thetaphi',subname)
!!$  thetaphi=0.0_gp
!!$
!!$  !calculate the scalar product between the hamiltonian and the gaussian basis
!!$  allocate(hpsigau(G%ncoeff,orbse%norbp+ndebug),stat=i_stat)
!!$  call memocc(i_stat,hpsigau,'hpsigau',subname)
!!$
!!$
!!$  call wavelets_to_gaussians(at%geocode,orbse%norbp,Glr%d%n1,Glr%d%n2,Glr%d%n3,G,&
!!$       thetaphi,hx,hy,hz,Glr%wfd,hpsi,hpsigau)
!!$
!!$  i_all=-product(shape(thetaphi))*kind(thetaphi)
!!$  deallocate(thetaphi,stat=i_stat)
!!$  call memocc(i_stat,i_all,'thetaphi',subname)


  accurex=abs(eks-ekin_sum)
  !tolerance for comparing the eigenvalues in the case of degeneracies
  etol=accurex/real(orbse%norbu,gp)
  if (iproc == 0 .and. verbose > 1) write(*,'(1x,a,2(f19.10))') 'done. ekin_sum,eks:',ekin_sum,eks
  if (iproc == 0) then
     write(*,'(1x,a,3(1x,1pe18.11))') 'ekin_sum,epot_sum,eproj_sum',  & 
          ekin_sum,epot_sum,eproj_sum
     write(*,'(1x,a,3(1x,1pe18.11))') '   ehart,   eexcu,    vexcu',ehart,eexcu,vexcu
  endif

!!$  call Gaussian_DiagHam(iproc,nproc,at%natsc,nspin,orbs,G,mpirequests,&
!!$       psigau,hpsigau,orbse,etol,norbsc_arr)

  !deallocate the gaussian basis descriptors
  call deallocate_gwf(G,subname)


!!$  i_all=-product(shape(mpirequests))*kind(mpirequests)
!!$  deallocate(mpirequests,stat=i_stat)
!!$  call memocc(i_stat,i_all,'mpirequests',subname)

  i_all=-product(shape(psigau))*kind(psigau)
  deallocate(psigau,stat=i_stat)
  call memocc(i_stat,i_all,'psigau',subname)
!!$  i_all=-product(shape(hpsigau))*kind(hpsigau)
!!$  deallocate(hpsigau,stat=i_stat)
!!$  call memocc(i_stat,i_all,'hpsigau',subname)



  !after having applied the hamiltonian to all the atomic orbitals
  !we split the semicore orbitals from the valence ones
  !this is possible since the semicore orbitals are the first in the 
  !order, so the linear algebra on the transposed wavefunctions 
  !may be splitted


  !free GPU if it is the case
  if (GPUconv) then
     call free_gpu(GPU,orbse%norbp)
  end if

  if (iproc == 0 .and. verbose > 1) write(*,'(1x,a)',advance='no')&
       'Input Wavefunctions Orthogonalization:'

  call DiagHam(iproc,nproc,at%natsc,nspin_ig,orbs,Glr%wfd,comms,&
       psi,hpsi,psit,orbse,commse,etol,norbsc_arr,orbsv,psivirt)
 
  call deallocate_comms(commse,subname)

  i_all=-product(shape(norbsc_arr))*kind(norbsc_arr)
  deallocate(norbsc_arr,stat=i_stat)
  call memocc(i_stat,i_all,'norbsc_arr',subname)

  if (iproc == 0 .and. verbose > 1) then
     write(*,'(1x,a)')'done.'
     !gaussian estimation valid only for Free BC
     if (at%geocode == 'F') then
        write(*,'(1x,a,1pe9.2)') &
          'expected accuracy in kinetic energy due to grid size',accurex
        write(*,'(1x,a,1pe9.2)') &
             'suggested value for gnrm_cv ',accurex/real(orbs%norb,kind=8)
     end if
  endif


  if (nvirt == 0) then
     i_all=-product(shape(orbsv%occup))*kind(orbsv%occup)
     deallocate(orbsv%occup,stat=i_stat)
     call memocc(i_stat,i_all,'orbsv%occup',subname)
     i_all=-product(shape(orbsv%spinsgn))*kind(orbsv%spinsgn)
     deallocate(orbsv%spinsgn,stat=i_stat)
     call memocc(i_stat,i_all,'orbsv%spinsgn',subname)
     i_all=-product(shape(orbsv%kpts))*kind(orbsv%kpts)
     deallocate(orbsv%kpts,stat=i_stat)
     call memocc(i_stat,i_all,'orbsv%kpts',subname)
     i_all=-product(shape(orbsv%kwgts))*kind(orbsv%kwgts)
     deallocate(orbsv%kwgts,stat=i_stat)
     call memocc(i_stat,i_all,'orbsv%kwgts',subname)
     i_all=-product(shape(orbsv%iokpt))*kind(orbsv%iokpt)
     deallocate(orbsv%iokpt,stat=i_stat)
     call memocc(i_stat,i_all,'orbsv%iokpt',subname)
  end if

  i_all=-product(shape(orbse%occup))*kind(orbse%occup)
  deallocate(orbse%occup,stat=i_stat)
  call memocc(i_stat,i_all,'orbse%occup',subname)
  i_all=-product(shape(orbse%spinsgn))*kind(orbse%spinsgn)
  deallocate(orbse%spinsgn,stat=i_stat)
  call memocc(i_stat,i_all,'orbse%spinsgn',subname)
  i_all=-product(shape(orbse%kpts))*kind(orbse%kpts)
  deallocate(orbse%kpts,stat=i_stat)
  call memocc(i_stat,i_all,'orbse%kpts',subname)
  i_all=-product(shape(orbse%kwgts))*kind(orbse%kwgts)
  deallocate(orbse%kwgts,stat=i_stat)
  call memocc(i_stat,i_all,'orbse%kwgts',subname)
  i_all=-product(shape(orbse%iokpt))*kind(orbse%iokpt)
  deallocate(orbse%iokpt,stat=i_stat)
  call memocc(i_stat,i_all,'orbse%iokpt',subname)

     
end subroutine input_wf_diag

 
