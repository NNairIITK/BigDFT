subroutine createWavefunctionsDescriptors(iproc,nproc,n1,n2,n3,output_grid,&
     hx,hy,hz,atoms,rxyz,radii_cf,crmult,frmult,&
     wfd,nvctrp,norb,norbp,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,bounds,nspinor,hybrid_on)
  !calculates the descriptor arrays keyg and keyv as well as nseg_c,nseg_f,nvctr_c,nvctr_f,nvctrp
  !calculates also the bounds arrays needed for convolutions
  use module_base
  use module_types
  implicit none
  !Arguments
  type(atoms_data), intent(in) :: atoms
  logical, intent(in) :: hybrid_on
  integer, intent(in) :: iproc,nproc,n1,n2,n3,norb,norbp,nspinor,output_grid
  integer, intent(in) :: nfl1,nfu1,nfl2,nfu2,nfl3,nfu3
  real(gp), intent(in) :: hx,hy,hz,crmult,frmult
  real(gp), dimension(3,atoms%nat), intent(in) :: rxyz
  real(gp), dimension(atoms%ntypes,3), intent(in) :: radii_cf
  integer, intent(out) :: nvctrp
  type(wavefunctions_descriptors) , intent(out) :: wfd
  type(convolutions_bounds), intent(out) :: bounds
  !local variables
  character(len=*), parameter :: subname='createWavefunctionsDescriptors'
  real(kind=8), parameter :: eps_mach=1.d-12
  integer :: iat,i1,i2,i3,norbme,norbyou,jpst,jproc,i_all,i_stat
  real(kind=8) :: tt
  logical, dimension(:,:,:), allocatable :: logrid_c,logrid_f

  !allocate kinetic bounds, only for free BC
  if (atoms%geocode == 'F') then
     allocate(bounds%kb%ibyz_c(2,0:n2,0:n3+ndebug),stat=i_stat)
     call memocc(i_stat,bounds%kb%ibyz_c,'bounds%kb%ibyz_c',subname)
     allocate(bounds%kb%ibxz_c(2,0:n1,0:n3+ndebug),stat=i_stat)
     call memocc(i_stat,bounds%kb%ibxz_c,'bounds%kb%ibxz_c',subname)
     allocate(bounds%kb%ibxy_c(2,0:n1,0:n2+ndebug),stat=i_stat)
     call memocc(i_stat,bounds%kb%ibxy_c,'bounds%kb%ibxy_c',subname)
     allocate(bounds%kb%ibyz_f(2,0:n2,0:n3+ndebug),stat=i_stat)
     call memocc(i_stat,bounds%kb%ibyz_f,'bounds%kb%ibyz_f',subname)
     allocate(bounds%kb%ibxz_f(2,0:n1,0:n3+ndebug),stat=i_stat)
     call memocc(i_stat,bounds%kb%ibxz_f,'bounds%kb%ibxz_f',subname)
     allocate(bounds%kb%ibxy_f(2,0:n1,0:n2+ndebug),stat=i_stat)
     call memocc(i_stat,bounds%kb%ibxy_f,'bounds%kb%ibxy_f',subname)
  end if

  if (iproc.eq.0) then
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
  call num_segkeys(n1,n2,n3,0,n1,0,n2,0,n3,logrid_c,wfd%nseg_c,wfd%nvctr_c)
  if (iproc.eq.0) write(*,'(2(1x,a,i10))') &
       'Coarse resolution grid: Number of segments= ',wfd%nseg_c,'points=',wfd%nvctr_c

  if (atoms%geocode == 'F') then
     call make_bounds(n1,n2,n3,logrid_c,bounds%kb%ibyz_c,bounds%kb%ibxz_c,bounds%kb%ibxy_c)
  end if

  if (atoms%geocode == 'P' .and. wfd%nvctr_c /= (n1+1)*(n2+1)*(n3+1) ) then
     if (iproc ==0)then
        write(*,*)&
          ' WARNING: the coarse grid does not fill the entire periodic box'
        write(*,*)&
          '          errors due to translational invariance breaking may occur'

     end if
  end if

  ! fine grid quantities
  call fill_logrid(atoms%geocode,n1,n2,n3,0,n1,0,n2,0,n3,0,atoms%nat,&
       atoms%ntypes,atoms%iatype,rxyz,radii_cf(1,2),frmult,hx,hy,hz,logrid_f)
  call num_segkeys(n1,n2,n3,0,n1,0,n2,0,n3,logrid_f,wfd%nseg_f,wfd%nvctr_f)
  if (iproc.eq.0) write(*,'(2(1x,a,i10))') &
       '  Fine resolution grid: Number of segments= ',wfd%nseg_f,'points=',wfd%nvctr_f
  if (atoms%geocode == 'F') then
     call make_bounds(n1,n2,n3,logrid_f,bounds%kb%ibyz_f,bounds%kb%ibxz_f,bounds%kb%ibxy_f)
  end if

  ! Create the file grid.xyz to visualize the grid of functions
  if (iproc ==0 .and. output_grid==1) then
     open(unit=22,file='grid.xyz',status='unknown')
     write(22,*) wfd%nvctr_c+wfd%nvctr_f+atoms%nat,' atomic'
     write(22,*)'complete simulation grid with low and high resolution points'
     do iat=1,atoms%nat
        write(22,'(a6,2x,3(1x,e12.5),3x)') &
             trim(atoms%atomnames(atoms%iatype(iat))),rxyz(1,iat),rxyz(2,iat),rxyz(3,iat)
     enddo
     do i3=0,n3  
        do i2=0,n2  
           do i1=0,n1
              if (logrid_c(i1,i2,i3))&
                   write(22,'(a4,2x,3(1x,e10.3))') &
                   '  g ',real(i1,kind=8)*hx,real(i2,kind=8)*hy,real(i3,kind=8)*hz
           enddo
        enddo
     end do
     do i3=0,n3 
        do i2=0,n2 
           do i1=0,n1
              if (logrid_f(i1,i2,i3))&
                   write(22,'(a4,2x,3(1x,e10.3))') &
                   '  G ',real(i1,kind=8)*hx,real(i2,kind=8)*hy,real(i3,kind=8)*hz
           enddo
        enddo
     enddo
     close(22)
  endif

  ! allocations for arrays holding the wavefunctions and their data descriptors
  call allocate_wfd(wfd,subname)

  ! now fill the wavefunction descriptor arrays
  ! coarse grid quantities
  call segkeys(n1,n2,n3,0,n1,0,n2,0,n3,logrid_c,wfd%nseg_c,wfd%keyg(1,1),wfd%keyv(1))

  ! fine grid quantities
  call segkeys(n1,n2,n3,0,n1,0,n2,0,n3,logrid_f,wfd%nseg_f,wfd%keyg(1,wfd%nseg_c+1), &
       & wfd%keyv(wfd%nseg_c+1))

  i_all=-product(shape(logrid_c))*kind(logrid_c)
  deallocate(logrid_c,stat=i_stat)
  call memocc(i_stat,i_all,'logrid_c',subname)
  i_all=-product(shape(logrid_f))*kind(logrid_f)
  deallocate(logrid_f,stat=i_stat)
  call memocc(i_stat,i_all,'logrid_f',subname)

  !distribution of wavefunction arrays between processors
  norbme=max(min((iproc+1)*norbp,norb)-iproc*norbp,0)
  !write(*,'(a,i0,a,i0,a)') '- iproc ',iproc,' treats ',norbme,' orbitals '
  if (iproc == 0 .and. nproc>1) then
     jpst=0
     do jproc=0,nproc-2
        norbme=max(min((jproc+1)*norbp,norb)-jproc*norbp,0)
        norbyou=max(min((jproc+2)*norbp,norb)-(jproc+1)*norbp,0)
        if (norbme /= norbyou) then
           !this is a screen output that must be modified
           write(*,'(3(a,i0),a)')&
                ' Processes from ',jpst,' to ',jproc,' treat ',norbme,' orbitals '
           jpst=jproc+1
        end if
     end do
     write(*,'(3(a,i0),a)')&
          ' Processes from ',jpst,' to ',nproc-1,' treat ',norbyou,' orbitals '
  end if

  tt=dble(wfd%nvctr_c+7*wfd%nvctr_f)/dble(nproc)
  nvctrp=int((1.d0-eps_mach*tt) + tt)

  if (iproc.eq.0) write(*,'(1x,a,i0)') &
       'Wavefunction memory occupation per orbital (Bytes): ',&
       nvctrp*nproc*8*nspinor !

  !for free BC admits the bounds arrays
  if (atoms%geocode == 'F') then

     !allocate grow, shrink and real bounds
     allocate(bounds%gb%ibzxx_c(2,0:n3,-14:2*n1+16+ndebug),stat=i_stat)
     call memocc(i_stat,bounds%gb%ibzxx_c,'bounds%gb%ibzxx_c',subname)
     allocate(bounds%gb%ibxxyy_c(2,-14:2*n1+16,-14:2*n2+16+ndebug),stat=i_stat)
     call memocc(i_stat,bounds%gb%ibxxyy_c,'bounds%gb%ibxxyy_c',subname)
     allocate(bounds%gb%ibyz_ff(2,nfl2:nfu2,nfl3:nfu3+ndebug),stat=i_stat)
     call memocc(i_stat,bounds%gb%ibyz_ff,'bounds%gb%ibyz_ff',subname)
     allocate(bounds%gb%ibzxx_f(2,nfl3:nfu3,2*nfl1-14:2*nfu1+16+ndebug),stat=i_stat)
     call memocc(i_stat,bounds%gb%ibzxx_f,'bounds%gb%ibzxx_f',subname)
     allocate(bounds%gb%ibxxyy_f(2,2*nfl1-14:2*nfu1+16,2*nfl2-14:2*nfu2+16+ndebug),stat=i_stat)
     call memocc(i_stat,bounds%gb%ibxxyy_f,'bounds%gb%ibxxyy_f',subname)

     allocate(bounds%sb%ibzzx_c(2,-14:2*n3+16,0:n1+ndebug),stat=i_stat)
     call memocc(i_stat,bounds%sb%ibzzx_c,'bounds%sb%ibzzx_c',subname)
     allocate(bounds%sb%ibyyzz_c(2,-14:2*n2+16,-14:2*n3+16+ndebug),stat=i_stat)
     call memocc(i_stat,bounds%sb%ibyyzz_c,'bounds%sb%ibyyzz_c',subname)
     allocate(bounds%sb%ibxy_ff(2,nfl1:nfu1,nfl2:nfu2+ndebug),stat=i_stat)
     call memocc(i_stat,bounds%sb%ibxy_ff,'bounds%sb%ibxy_ff',subname)
     allocate(bounds%sb%ibzzx_f(2,-14+2*nfl3:2*nfu3+16,nfl1:nfu1+ndebug),stat=i_stat)
     call memocc(i_stat,bounds%sb%ibzzx_f,'bounds%sb%ibzzx_f',subname)
     allocate(bounds%sb%ibyyzz_f(2,-14+2*nfl2:2*nfu2+16,-14+2*nfl3:2*nfu3+16+ndebug),stat=i_stat)
     call memocc(i_stat,bounds%sb%ibyyzz_f,'bounds%sb%ibyyzz_f',subname)

     allocate(bounds%ibyyzz_r(2,-14:2*n2+16,-14:2*n3+16+ndebug),stat=i_stat)
     call memocc(i_stat,bounds%ibyyzz_r,'bounds%ibyyzz_r',subname)

     call make_all_ib(n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,&
          bounds%kb%ibxy_c,bounds%sb%ibzzx_c,bounds%sb%ibyyzz_c,&
          bounds%kb%ibxy_f,bounds%sb%ibxy_ff,bounds%sb%ibzzx_f,bounds%sb%ibyyzz_f,&
          bounds%kb%ibyz_c,bounds%gb%ibzxx_c,bounds%gb%ibxxyy_c,&
          bounds%kb%ibyz_f,bounds%gb%ibyz_ff,bounds%gb%ibzxx_f,bounds%gb%ibxxyy_f,&
          bounds%ibyyzz_r)

  end if

!*************Added by Alexey************************************************************
  if ( atoms%geocode == 'P' .and. hybrid_on) then
	  call make_bounds_per(n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,bounds,wfd)
	  call make_all_ib_per(n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,&
     bounds%kb%ibxy_f,bounds%sb%ibxy_ff,bounds%sb%ibzzx_f,bounds%sb%ibyyzz_f,&
     bounds%kb%ibyz_f,bounds%gb%ibyz_ff,bounds%gb%ibzxx_f,bounds%gb%ibxxyy_f)
  endif	
!****************************************************************************************  
END SUBROUTINE createWavefunctionsDescriptors

subroutine make_bounds_per(n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,bounds,wfd)
use module_base
use module_types
implicit none
type(wavefunctions_descriptors), intent(in) :: wfd
type(convolutions_bounds),intent(out):: bounds

integer,intent(in)::n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3

logical,allocatable,dimension(:,:,:)::logrid
character(len=*), parameter :: subname='make_bounds'
integer::i_stat,i_all,nseg_c,i2,i3

allocate(bounds%kb%ibyz_f(2,0:n2,0:n3+ndebug),stat=i_stat)
call memocc(i_stat,bounds%kb%ibyz_f,'bounds%kb%ibyz_f',subname)
allocate(bounds%kb%ibxz_f(2,0:n1,0:n3+ndebug),stat=i_stat)
call memocc(i_stat,bounds%kb%ibxz_f,'bounds%kb%ibxz_f',subname)
allocate(bounds%kb%ibxy_f(2,0:n1,0:n2+ndebug),stat=i_stat)
call memocc(i_stat,bounds%kb%ibxy_f,'bounds%kb%ibxy_f',subname)

allocate(bounds%gb%ibyz_ff(2,nfl2:nfu2,nfl3:nfu3+ndebug),stat=i_stat)
call memocc(i_stat,bounds%gb%ibyz_ff,'bounds%gb%ibyz_ff',subname)
allocate(bounds%gb%ibzxx_f(2,nfl3:nfu3,0:2*n1+1+ndebug),stat=i_stat)
call memocc(i_stat,bounds%gb%ibzxx_f,'bounds%gb%ibzxx_f',subname)
allocate(bounds%gb%ibxxyy_f(2,0:2*n1+1,0:2*n2+1+ndebug),stat=i_stat)
call memocc(i_stat,bounds%gb%ibxxyy_f,'bounds%gb%ibxxyy_f',subname)

allocate(bounds%sb%ibxy_ff(2,nfl1:nfu1,nfl2:nfu2+ndebug),stat=i_stat)
call memocc(i_stat,bounds%sb%ibxy_ff,'bounds%sb%ibxy_ff',subname)
allocate(bounds%sb%ibzzx_f(2,0:2*n3+1,nfl1:nfu1+ndebug),stat=i_stat)
call memocc(i_stat,bounds%sb%ibzzx_f,'bounds%sb%ibzzx_f',subname)
allocate(bounds%sb%ibyyzz_f(2,0:2*n2+1,0:2*n3+1+ndebug),stat=i_stat)
call memocc(i_stat,bounds%sb%ibyyzz_f,'bounds%sb%ibyyzz_f',subname)

allocate(logrid(0:n1,0:n2,0:n3+ndebug),stat=i_stat)
call memocc(i_stat,logrid,'logrid',subname)

nseg_c=wfd%nseg_c
call make_logrid_f(n1,n2,n3, & 
     wfd%nseg_c,wfd%nvctr_c,wfd%keyg(1,1),wfd%keyv(1),  & 
     wfd%nseg_f,wfd%nvctr_f,wfd%keyg(1,nseg_c+1),wfd%keyv(nseg_c+1),  & 
     logrid)
	 
call make_bounds(n1,n2,n3,logrid,bounds%kb%ibyz_f,bounds%kb%ibxz_f,bounds%kb%ibxy_f)

i_all=-product(shape(logrid))*kind(logrid)
deallocate(logrid,stat=i_stat)
call memocc(i_stat,i_all,'logrid',subname)

end subroutine make_bounds_per


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
  integer :: nl1,nl2,nl3,nu1,nu2,nu3,mseg,mvctr,mproj,istart
  integer :: iat,i_stat,i_all,ityp,iseg,natyp
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
     nvctrp,Glr,hx,hy,hz,rxyz,rhopot,pot_ion,nlpspd,proj,& 
     pkernel,ixc,psi,psit,hpsi,nscatterarr,ngatherarr,nspin)
  use module_base
  use module_interfaces, except_this_one => import_gaussians
  use module_types
  use Poisson_Solver
  implicit none
  integer, intent(in) :: iproc,nproc,ixc,nvctrp,nspin
  real(gp), intent(in) :: hx,hy,hz,cpmult,fpmult
  type(atoms_data), intent(in) :: at
  type(orbitals_data), intent(in) :: orbs
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
  integer :: i,iorb,i_stat,i_all,ierr,info,jproc,n_lp,jorb,j
  real(kind=4) :: t1,t0
  real(gp) :: hxh,hyh,hzh,eks,eexcu,vexcu,epot_sum,ekin_sum,ehart,eproj_sum,accurex,maxdiff
  type(gaussian_basis) :: CP2K
  integer, dimension(:), allocatable :: iwork
  real(gp), dimension(:), allocatable :: ones,ovrlp,work
  real(gp), dimension(:,:), allocatable :: tmp,smat
  real(wp), dimension(:,:), pointer :: wfn_cp2k

  if (iproc.eq.0) then
     write(*,'(1x,a)')&
          '--------------------------------------------------------- Import Gaussians from CP2K'
  end if

  if (nspin /= 1) then
     if (iproc==0) write(*,'(1x,a)')&
          'Gaussian importing is possible only for non-spin polarised calculations'
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

  call gaussians_to_wavelets(at%geocode,iproc,nproc,orbs,Glr%d%n1,Glr%d%n2,Glr%d%n3,&
     Glr%d%nfl1,Glr%d%nfu1,Glr%d%nfl2,Glr%d%nfu2,Glr%d%nfl3,Glr%d%nfu3,&
     hx,hy,hz,Glr%wfd,CP2K,wfn_cp2k,psi)

  !deallocate CP2K variables
  call deallocate_gwf(CP2K,subname)
  !nullify gaussian centers
  nullify(CP2K%rxyz)
  i_all=-product(shape(wfn_cp2k))*kind(wfn_cp2k)
  deallocate(wfn_cp2k,stat=i_stat)
  call memocc(i_stat,i_all,'wfn_cp2k',subname)

  call sumrho(iproc,nproc,orbs,Glr,ixc,hxh,hyh,hzh,psi,rhopot,&
       Glr%d%n1i*Glr%d%n2i*nscatterarr(iproc,1),nscatterarr,1)

  call PSolver(at%geocode,'D',iproc,nproc,Glr%d%n1i,Glr%d%n2i,Glr%d%n3i,ixc,hxh,hyh,hzh,&
       rhopot,pkernel,pot_ion,ehart,eexcu,vexcu,0.0_dp,.true.,1)

  !allocate the wavefunction in the transposed way to avoid allocations/deallocations
  allocate(hpsi(npsidim+ndebug),stat=i_stat)
  call memocc(i_stat,hpsi,'hpsi',subname)

  call HamiltonianApplication(iproc,nproc,at,orbs,hx,hy,hz,rxyz,cpmult,fpmult,radii_cf,&
       nlpspd,proj,Glr,ngatherarr,&
       Glr%d%n1i*Glr%d%n2i*nscatterarr(iproc,2),&
       rhopot(1+Glr%d%n1i*Glr%d%n2i*nscatterarr(iproc,4)),&
       psi,hpsi,ekin_sum,epot_sum,eproj_sum,1)

  accurex=abs(eks-ekin_sum)
  if (iproc.eq.0) write(*,'(1x,a,2(f19.10))') 'done. ekin_sum,eks:',ekin_sum,eks

  if (iproc.eq.0) then
     write(*,'(1x,a,3(1x,1pe18.11))') 'ekin_sum,epot_sum,eproj_sum',  & 
          ekin_sum,epot_sum,eproj_sum
     write(*,'(1x,a,3(1x,1pe18.11))') '   ehart,   eexcu,    vexcu',ehart,eexcu,vexcu
  endif

  !after having applied the hamiltonian to all the atomic orbitals
  !we split the semicore orbitals from the valence ones
  !this is possible since the semicore orbitals are the first in the 
  !order, so the linear algebra on the transposed wavefunctions 
  !may be splitted

  if (iproc.eq.0) write(*,'(1x,a)',advance='no')&
       'Imported Wavefunctions Orthogonalization:'

  call DiagHam(iproc,nproc,at%natsc,nspin,orbs,nvctrp,Glr%wfd,comms,psi,hpsi,psit)

  if (iproc.eq.0) write(*,'(1x,a)')'done.'

END SUBROUTINE import_gaussians

subroutine input_wf_diag(iproc,nproc,cpmult,fpmult,radii_cf,at,&
     nvirte,nvirtep,nvirt,nvctrp,Glr,hx,hy,hz,rxyz,rhopot,pot_ion,&
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
  !!!!!!!!!!SONO ARRIVATO QUI CON LA CONVERSIONE
  integer, intent(in) :: iproc,nproc,norb,norbp,ixc,nvctrp
  integer, intent(inout) :: nspin,nvirte,nvirtep,nvirt
  real(gp), intent(in) :: hx,hy,hz,cpmult,fpmult
  type(atoms_data), intent(in) :: at
  type(nonlocal_psp_descriptors), intent(in) :: nlpspd
  type(locreg_descriptors), intent(in) :: Glr
  integer, dimension(0:nproc-1,4), intent(in) :: nscatterarr !n3d,n3p,i3s+i3xcsh-1,i3xcsh
  integer, dimension(0:nproc-1,2), intent(in) :: ngatherarr 
  real(gp), dimension(norb), intent(in) :: spinsgn
  real(gp), dimension(at%ntypes,3), intent(in) :: radii_cf  
  real(gp), dimension(3,at%nat), intent(in) :: rxyz
  real(wp), dimension(nlpspd%nprojel), intent(in) :: proj
  real(dp), dimension(*), intent(in) :: pkernel
  real(dp), dimension(*), intent(inout) :: rhopot,pot_ion
  real(wp), dimension(norb), intent(out) :: eval
  real(wp), dimension(:), pointer :: psi,hpsi,psit,psivirt
  !local variables
  character(len=*), parameter :: subname='input_wf_diag'
  real(kind=8), parameter :: eps_mach=1.d-12
  integer, parameter :: ngx=31
  integer :: i,iorb,iorbsc,imatrsc,iorbst,imatrst,i_stat,i_all,ierr,info
  integer :: norbe,norbep,norbi,norbj,norbeme,ndim_hamovr,n_lp,norbsc,jproc,jpst,norbeyou
  integer :: ispin,norbu,norbd,iorbst2,ist,n2hamovr,nsthamovr,nspinor,iat
  real(gp) :: hxh,hyh,hzh,tt,eks,eexcu,vexcu,epot_sum,ekin_sum,ehart,eproj_sum,etol,accurex
  type(gaussian_basis) :: G
  logical, dimension(:,:,:), allocatable :: scorb
  integer, dimension(:), allocatable :: ng,iorbtolr
  integer, dimension(:,:), allocatable :: nl,norbsc_arr
  real(wp), dimension(:,:), allocatable :: gaucoeff
  real(gp), dimension(:), allocatable :: occupe,spinsgne,locrad
  real(gp), dimension(:,:), allocatable :: xp,occupat
  real(gp), dimension(:,:,:), allocatable :: psiat
  type(locreg_descriptors), dimension(:), allocatable :: Llr

  !Calculate no. up and down orbitals for spin-polarized starting guess
  norbu=0
  norbd=0
  do iorb=1,norb
     if(spinsgn(iorb)>0.0_gp) norbu=norbu+1
     if(spinsgn(iorb)<0.0_gp) norbd=norbd+1
  end do
  if(nspin==4) then
     nspinor=4
     nspin=2
  else
     nspinor=1
  end if
!!$  !calculate dimension of the interpolating scaling function grid
!!$  select case(at%geocode)
!!$     case('F')
!!$        n1i=2*n1+31
!!$        n2i=2*n2+31
!!$        n3i=2*n3+31
!!$     case('S')
!!$        n1i=2*n1+2
!!$        n2i=2*n2+31
!!$        n3i=2*n3+2
!!$     case('P')
!!$        n1i=2*n1+2
!!$        n2i=2*n2+2
!!$        n3i=2*n3+2
!!$  end select

  allocate(xp(ngx,at%ntypes+ndebug),stat=i_stat)
  call memocc(i_stat,xp,'xp',subname)
  allocate(psiat(ngx,5,at%ntypes+ndebug),stat=i_stat)
  call memocc(i_stat,psiat,'psiat',subname)
  allocate(occupat(5,at%ntypes+ndebug),stat=i_stat)
  call memocc(i_stat,occupat,'occupat',subname)
  allocate(ng(at%ntypes+ndebug),stat=i_stat)
  call memocc(i_stat,ng,'ng',subname)
  allocate(nl(4,at%ntypes+ndebug),stat=i_stat)
  call memocc(i_stat,nl,'nl',subname)
  allocate(scorb(4,2,at%natsc+ndebug),stat=i_stat)
  call memocc(i_stat,scorb,'scorb',subname)
  allocate(norbsc_arr(at%natsc+1,nspin+ndebug),stat=i_stat)
  call memocc(i_stat,norbsc_arr,'norbsc_arr',subname)
  allocate(locrad(at%nat+ndebug),stat=i_stat)
  call memocc(i_stat,locrad,'locrad',subname)

  if (iproc.eq.0) then
     write(*,'(1x,a)')&
          '------------------------------------------------------- Input Wavefunctions Creation'
  end if

  !Generate the input guess via the inguess_generator
  !here we should allocate the gaussian basis descriptors 
  !the prescriptions can be found in the creation of psp basis

  call readAtomicOrbitals(iproc,ngx,xp,psiat,occupat,ng,nl,at,norbe,norbsc,nspin,&
       scorb,norbsc_arr,locrad)

  !Check for max number of virtual orbitals
  nvirte=norbe-max(norbu,norbd)!the unoccupied orbitals available as a LCAO
  if(nvirt==nvirte .and. nvirt/=0 .and. iproc==0) then
     write(*,'(1x,a)')&
          "WARNING: A smaller number of virtual orbitals may be needed for better convergence."
     write(*,'(1x,a,i0)')'         Put nvirte= ',nvirte
  end if
  if(nvirte<nvirt)then
     nvirt=nvirte 
     if(iproc==0)write(*,'(1x,a,i3)')&
          "WARNING: Number of virtual orbitals is too large. New value: ",nvirt
  end if
  
  !no Davidson calculation if nvirt=0
  if (nvirt==0) nvirte=0

  !  allocate wavefunctions and their occupation numbers
  allocate(occupe(nspin*norbe+ndebug),stat=i_stat)
  call memocc(i_stat,occupe,'occupe',subname)
  !the number of orbitals to be considered is doubled in the case of a spin-polarised calculation
  tt=dble(nspin*norbe)/dble(nproc)
  norbep=int((1.d0-eps_mach*tt) + tt)


  !this is the distribution procedure for cubic code
  if (iproc == 0 .and. nproc>1) then
     jpst=0
     do jproc=0,nproc-2
        norbeme=max(min((jproc+1)*norbep,nspin*norbe)-jproc*norbep,0)
        norbeyou=max(min((jproc+2)*norbep,nspin*norbe)-(jproc+1)*norbep,0)
        if (norbeme /= norbeyou) then
           !this is a screen output that must be modified
           write(*,'(3(a,i0),a)')&
                ' Processes from ',jpst,' to ',jproc,' treat ',norbeme,' inguess orbitals '
           jpst=jproc+1
        end if
     end do
     write(*,'(3(a,i0),a)')&
          ' Processes from ',jpst,' to ',nproc-1,' treat ',norbeyou,' inguess orbitals '
  end if
  

  hxh=.5_gp*hx
  hyh=.5_gp*hy
  hzh=.5_gp*hz

  !allocate the wavefunction in the transposed way to avoid allocations/deallocations
!  allocate(psi(nvctrp,norbep*nproc*nspinor),stat=i_stat)
  allocate(psi(nvctrp*norbep*nproc+ndebug),stat=i_stat)
  call memocc(i_stat,psi,'psi',subname)
  
  ! Create input guess orbitals
  !call createAtomicOrbitals(iproc,nproc,at,rxyz,norbe,norbep,norbsc,occupe,occupat,&
  !     ngx,xp,psiat,ng,nl,wfd,n1,n2,n3,hx,hy,hz,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,nspin,psi,eks,scorb)
  allocate(gaucoeff(norbe,norbep+ndebug),stat=i_stat)
  call memocc(i_stat,gaucoeff,'gaucoeff',subname)
  allocate(iorbtolr(norbep+ndebug),stat=i_stat)
  call memocc(i_stat,iorbtolr,'iorbtolr',subname)

  
  call AtomicOrbitals(iproc,nproc,at,rxyz,norbe,norbep,norbsc,occupe,occupat,&
     ngx,xp,psiat,ng,nl,nspin,eks,scorb,G,gaucoeff,iorbtolr)!,&
     !wfd,n1,n2,n3,hx,hy,hz,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,psi)
  !createAtomicOrbitals should generate the gaussian basis set and the data needed to generate
  !the wavefunctions inside a given localisation region.

  !create the localisation region which are associated to the gaussian extensions and plot them
!!$  Glr%geocode=at%geocode
!!$  Glr%ns1=0
!!$  Glr%ns2=0
!!$  Glr%ns3=0
!!$  Glr%d%n1=n1
!!$  Glr%d%n2=n2
!!$  Glr%d%n3=n3
!!$  Glr%d%nfl1=nfl1
!!$  Glr%d%nfl2=nfl2
!!$  Glr%d%nfl3=nfl3
!!$  Glr%d%nfu1=nfu1
!!$  Glr%d%nfu2=nfu2
!!$  Glr%d%nfu3=nfu3
!!$  Glr%wfd=wfd !to be tested

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

  call gaussians_to_wavelets(at%geocode,iproc,nproc,norbe*nspin,norbep,&
     Glr%d%n1,Glr%d%n2,Glr%d%n3,&
     Glr%d%nfl1,Glr%d%nfu1,Glr%d%nfl2,Glr%d%nfu2,Glr%d%nfl3,Glr%d%nfu3,&
     hx,hy,hz,Glr%wfd,G,gaucoeff,psi)

  
!!$  !!plot the initial LCAO wavefunctions
!!$  !do i=2*iproc+1,2*iproc+2
!!$  !   iounit=15+3*(i-1)
!!$  !   print *,'iounit',iounit,'-',iounit+2
!!$  !   call plot_wf(iounit,n1,n2,n3,hgrid,nseg_c,nvctr_c,keyg,keyv,nseg_f,nvctr_f,  & 
!!$  !        rxyz(1,1),rxyz(2,1),rxyz(3,1),psi(:,i-2*iproc:i-2*iproc), &
!!$  !        ibyz_c,ibzxx_c,ibxxyy_c,ibyz_ff,ibzxx_f,ibxxyy_f,ibyyzz_r,&
!!$  !        nfl1,nfu1,nfl2,nfu2,nfl3,nfu3)
!!$  !end do

  !deallocate the gaussian basis descriptors
  call deallocate_gwf(G,subname)


  i_all=-product(shape(gaucoeff))*kind(gaucoeff)
  deallocate(gaucoeff,stat=i_stat)
  call memocc(i_stat,i_all,'gaucoeff',subname)
  i_all=-product(shape(locrad))*kind(locrad)
  deallocate(locrad,stat=i_stat)
  call memocc(i_stat,i_all,'locrad',subname)
  i_all=-product(shape(iorbtolr))*kind(iorbtolr)
  deallocate(iorbtolr,stat=i_stat)
  call memocc(i_stat,i_all,'iorbtolr',subname)


  
  i_all=-product(shape(scorb))*kind(scorb)
  deallocate(scorb,stat=i_stat)
  call memocc(i_stat,i_all,'scorb',subname)
  i_all=-product(shape(xp))*kind(xp)
  deallocate(xp,stat=i_stat)
  call memocc(i_stat,i_all,'xp',subname)
  i_all=-product(shape(psiat))*kind(psiat)
  deallocate(psiat,stat=i_stat)
  call memocc(i_stat,i_all,'psiat',subname)
  i_all=-product(shape(occupat))*kind(occupat)
  deallocate(occupat,stat=i_stat)
  call memocc(i_stat,i_all,'occupat',subname)
  i_all=-product(shape(ng))*kind(ng)
  deallocate(ng,stat=i_stat)
  call memocc(i_stat,i_all,'ng',subname)
  i_all=-product(shape(nl))*kind(nl)
  deallocate(nl,stat=i_stat)
  call memocc(i_stat,i_all,'nl',subname)


  ! resulting charge density and potential
  allocate(spinsgne(nspin*norbe+ndebug),stat=i_stat)
  call memocc(i_stat,spinsgne,'spinsgne',subname)
  ist=1
  do ispin=1,nspin
     spinsgne(ist:ist+norbe-1)=real(1-2*(ispin-1),gp)
     ist=norbe+1
  end do

!!$  !call the gaussian basis structure associated to the input guess
!!$  call gaussian_pswf_basis(iproc,at,rxyz,G)
!!$  !create the density starting from input guess gaussians
    
  call sumrho(iproc,nproc,nspin*norbe,norbep,Glr,ixc,hxh,hyh,hzh,occupe,  & 
       psi,rhopot,Glr%d%n1i*Glr%d%n2i*nscatterarr(iproc,1),nscatterarr,&
       nspin,1,spinsgne,hybrid_on)

  call PSolver(at%geocode,'D',iproc,nproc,Glr%d%n1i,Glr%d%n2i,Glr%d%n3i,ixc,hxh,hyh,hzh,&
       rhopot,pkernel,pot_ion,ehart,eexcu,vexcu,0.0_dp,.true.,nspin)

  !allocate the wavefunction in the transposed way to avoid allocations/deallocations
  allocate(hpsi(nvctrp*norbep*nproc+ndebug),stat=i_stat)
  call memocc(i_stat,hpsi,'hpsi',subname)
  
  call HamiltonianApplication(iproc,nproc,at,hx,hy,hz,rxyz,cpmult,fpmult,radii_cf,&
       nspin*norbe,norbep,occupe,nlpspd,proj,Glr,&
       ngatherarr,Glr%d%n1i*Glr%d%n2i*nscatterarr(iproc,2),&
       rhopot(1+Glr%d%n1i*Glr%d%n2i*nscatterarr(iproc,4)),&
       psi,hpsi,ekin_sum,epot_sum,eproj_sum,nspin,1,spinsgne,hybrid_on)

  i_all=-product(shape(spinsgne))*kind(spinsgne)
  deallocate(spinsgne,stat=i_stat)
  call memocc(i_stat,i_all,'spinsgne',subname)

  i_all=-product(shape(occupe))*kind(occupe)
  deallocate(occupe,stat=i_stat)
  call memocc(i_stat,i_all,'occupe',subname)

  accurex=abs(eks-ekin_sum)
  !tolerance for comparing the eigenvalues in the case of degeneracies
  etol=accurex/real(norbe,gp)
  if (iproc.eq.0) write(*,'(1x,a,2(f19.10))') 'done. ekin_sum,eks:',ekin_sum,eks
  if (iproc.eq.0) then
     write(*,'(1x,a,3(1x,1pe18.11))') 'ekin_sum,epot_sum,eproj_sum',  & 
          ekin_sum,epot_sum,eproj_sum
     write(*,'(1x,a,3(1x,1pe18.11))') '   ehart,   eexcu,    vexcu',ehart,eexcu,vexcu
  endif


  !after having applied the hamiltonian to all the atomic orbitals
  !we split the semicore orbitals from the valence ones
  !this is possible since the semicore orbitals are the first in the 
  !order, so the linear algebra on the transposed wavefunctions 
  !may be splitted


  if (iproc.eq.0) write(*,'(1x,a)',advance='no')&
       'Input Wavefunctions Orthogonalization:'

  call DiagHam(iproc,nproc,at%natsc,nspin,nspinor,norbu,norbd,norb,norbp,nvctrp,Glr%wfd,&
       psi,hpsi,psit,eval,norbe,norbep,etol,norbsc_arr,nvirte,nvirtep,psivirt)
 

  i_all=-product(shape(norbsc_arr))*kind(norbsc_arr)
  deallocate(norbsc_arr,stat=i_stat)
  call memocc(i_stat,i_all,'norbsc_arr',subname)

  if(nspinor==4) nspin=4

!  if(nspin==4) then
!     call psitospi(iproc,nproc,norbe,norbep,norbsc,nat,&
!          wfd%nvctr_c,wfd%nvctr_f,at%iatype,at%ntypes,&
!          at%iasctype,at%natsc,at%natpol,nspin,spinsgne,psi)
!  end if

  if (iproc.eq.0) then
     write(*,'(1x,a)')'done.'
     !gaussian estimation valid only for Free BC
     if (at%geocode == 'F') then
        write(*,'(1x,a,1pe9.2)') &
          'expected accuracy in kinetic energy due to grid size',accurex
        write(*,'(1x,a,1pe9.2)') &
             'suggested value for gnrm_cv ',accurex/real(norb,kind=8)
     end if
  endif
     
end subroutine input_wf_diag

!!****f* BigDFT/DiagHam
!! DESCRIPTION
!!    Diagonalise the hamiltonian in a basis set of norbe orbitals and select the first
!!    norb eigenvectors. Works also with the spin-polarisation case and perform also the 
!!    treatment of semicore atoms. 
!!    In the absence of norbe parameters, it simply diagonalize the hamiltonian in the given
!!    orbital basis set.
!! COPYRIGHT
!!    Copyright (C) 2008 CEA Grenoble
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!! INPUT VARIABLES
!!    iproc  process id
!!    nproc  number of mpi processes
!!    natsc  number of semicore atoms for the orthogonalisation treatment
!!           used as a dimension for the array of semicore atoms
!!    nspin  spin polarised id; 1 => non spin-polarised; 2 => spin-polarised (collinear)
!!    norbu  number of up orbitals in the spin-polarised case; for non spin-pol equal to norb
!!    norbd  number of down orbitals in the spin-polarised case; for non spin-pol equal to 0
!!    norb   total number of orbitals of the resulting eigenfunctions
!!    norbp  number of orbitals in parallel. For nproc=1 norbp=norb
!!    nvirte  number of virtual orbitals to be saved as input guess for the Davidson method
!!    nvctrp number of points of the wavefunctions for each orbital in the transposed sense
!!    wfd    data structure of the wavefunction descriptors
!!    norbe  (optional) number of orbitals of the initial set of wavefunction, to be reduced
!!    etol   tolerance for which a degeneracy should be printed. Set to zero if absent
!! INPUT-OUTPUT VARIABLES
!!    psi    wavefunctions. 
!!           If norbe is absent: on input, set of norb wavefunctions, 
!!                               on output eigenfunctions
!!           If norbe is present: on input, set of norbe wavefunctions, 
!!                                on output the first norb eigenfunctions
!!    hpsi   hamiltonian on the wavefunctions
!!           If norbe is absent: on input, set of norb arrays, 
!!                               destroyed on output
!!           If norbe is present: on input, set of norbe wavefunctions, 
!!                                destroyed on output
!! OUTPUT VARIABLES
!!    psit   wavefunctions in the transposed form.
!!           On input: nullified
!!           on Output: transposed wavefunction but only if nproc>1, nullified otherwise
!!    psivirt wavefunctions for input guess of the Davidson method in the transposed form.
!!           On input: nullified
!!           if nvirte >0: on Output transposed wavefunction (if nproc>1), direct otherwise
!!           if nvirte=0: nullified
!!    eval   array of the first norb eigenvalues       
!! AUTHOR
!!    Luigi Genovese
!! CREATION DATE
!!    February 2008
!! SOURCE
!! 
subroutine DiagHam(iproc,nproc,natsc,nspin,orbs,nvctrp,wfd,comms,&
     psi,hpsi,psit,& !mandatory
     norbe,norbep,norbsc_arr,nvirte,nvirtep,psivirt) !optional
  use module_base
  use module_types
  use module_interfaces, except_this_one => DiagHam
  implicit none
  integer, intent(in) :: iproc,nproc,natsc,nspin,nvctrp
  type(wavefunctions_descriptors), intent(in) :: wfd
  type(communications_arrays), intent(in) :: comms
  type(orbitals_data), intent(inout) :: orbs
  real(wp), dimension(:), pointer :: psi,hpsi,psit
  !optional arguments
  integer, optional, intent(in) :: norbe,norbep,nvirte
  integer, optional, intent(out) :: nvirtep
  real(gp), optional, intent(in) :: etol
  integer, optional, dimension(natsc+1,nspin), intent(in) :: norbsc_arr
  real(wp), dimension(:), pointer, optional :: psivirt
   !real(kind=8), optional, dimension(:,:), pointer :: ppsi
  !local variables
  character(len=*), parameter :: subname='DiagHam'
  real(kind=8), parameter :: eps_mach=1.d-12
  logical :: semicore,minimal
  integer :: i,ndim_hamovr,i_all,i_stat,n2hamovr,nsthamovr,ierr,norbi_max,j
  integer :: norbtot,norbtotp,natsceff,norbsc,ndh1,ispin,nvctr
  real(gp) :: tolerance
  real(kind=8) :: tt
  integer, dimension(:,:), allocatable :: norbgrp
  real(wp), dimension(:,:), allocatable :: hamovr
  real(wp), dimension(:), pointer :: psiw

  !performs some check of the arguments
  if (present(etol)) then
     tolerance=etol
  else
     tolerance=0.0_gp
  end if

  if (present(norbe) .neqv. present(norbep)) then
     if (iproc ==0) write(*,'(1x,a)')&
          'ERROR (DiagHam): the variables norbe and norbep must be present at the same time'
     stop
  else
     minimal=present(norbe)
  end if

  semicore=present(norbsc_arr)

  !define the grouping of the orbitals: for the semicore case, follow the semocore atoms,
  !otherwise use the number of orbitals, separated in the spin-polarised case
  !fro the spin polarised case it is supposed that the semicore orbitals are disposed equally
  if (semicore) then
     norbi_max=max(maxval(norbsc_arr),nvirte)

     !calculate the dimension of the overlap matrix
     !take the maximum as the two spin dimensions
     ndim_hamovr=0
     do ispin=1,nspin
        ndh1=0
        norbsc=0
        do i=1,natsc+1
           ndh1=ndh1+norbsc_arr(i,ispin)**2
        end do
        ndim_hamovr=max(ndim_hamovr,ndh1)
     end do
     if (natsc > 0) then
        if (nspin == 2) then
           if (sum(norbsc_arr(1:natsc,1)) /= sum(norbsc_arr(1:natsc,2))) then
              write(*,'(1x,a)')&
                'ERROR (DiagHam): The number of semicore orbitals must be the same for both spins'
              stop
           end if
        end if
        norbsc=sum(norbsc_arr(1:natsc,1))
     else
        norbsc=0
     end if

     natsceff=natsc
     allocate(norbgrp(natsceff+1,nspin+ndebug),stat=i_stat)
     call memocc(i_stat,norbgrp,'norbgrp',subname)

     !assign the grouping of the orbitals
     do j=1,nspin
        do i=1,natsceff+1
           norbgrp(i,j)=norbsc_arr(i,j)
        end do
     end do
  else
     norbi_max=max(norbu,norbd) !this works also for non spin-polarised since there norbu=norb
     ndim_hamovr=norbi_max**2

     natsceff=0
     allocate(norbgrp(1,nspin+ndebug),stat=i_stat)
     call memocc(i_stat,norbgrp,'norbgrp',subname)

     norbsc=0
     norbgrp(1,1)=norbu
     if (nspin == 2) norbgrp(1,2)=norbd

  end if

  !assign total orbital number for calculating the overlap matrix and diagonalise the system
  if(minimal) then
     norbtot=nspin*norbe !beware that norbe is equal both for spin up and down
     norbtotp=norbep !this is coherent with nspin*norbe
  else
     norbtot=norb
     norbtotp=norbp
  end if
  if (nproc > 1) then
     allocate(psiw(nvctrp*norbtotp*nproc+ndebug),stat=i_stat)
     call memocc(i_stat,psiw,'psiw',subname)
  else
     psiw => null()
  end if

  !transpose all the wavefunctions for having a piece of all the orbitals 
  !for each processor
  call transpose_v(iproc,nproc,norbtotp,1,wfd,nvctrp,comms,psi,work=psiw)
  call transpose_v(iproc,nproc,norbtotp,1,wfd,nvctrp,comms,hpsi,work=psiw)

 if (nproc > 1) then
     i_all=-product(shape(psiw))*kind(psiw)
     deallocate(psiw,stat=i_stat)
     call memocc(i_stat,i_all,'psiw',subname)

     n2hamovr=4
     nsthamovr=3
  else
     !allocation values
     n2hamovr=2
     nsthamovr=1
  end if

  allocate(hamovr(nspin*ndim_hamovr,n2hamovr+ndebug),stat=i_stat)
  call memocc(i_stat,hamovr,'hamovr',subname)

  if (iproc.eq.0) write(*,'(1x,a)',advance='no')&
       'Overlap Matrix...'

  call overlap_matrices(norbtot,nvctrp,natsceff,nspin,ndim_hamovr,norbgrp,&
       hamovr(1,nsthamovr),psi,hpsi)

  if (minimal) then
     !deallocate hpsi in the case of a minimal basis
     i_all=-product(shape(hpsi))*kind(hpsi)
     deallocate(hpsi,stat=i_stat)
     call memocc(i_stat,i_all,'hpsi',subname)
  end if

  if (nproc > 1) then
     !reduce the overlap matrix between all the processors
     call MPI_ALLREDUCE(hamovr(1,3),hamovr(1,1),2*nspin*ndim_hamovr,&
          mpidtypw,MPI_SUM,MPI_COMM_WORLD,ierr)
  end if

  call solve_eigensystem(iproc,orbs%norb,orbs%norbu,orbs%norbd,norbi_max,&
       ndim_hamovr,natsceff,nspin,tolerance,norbgrp,hamovr,orbs%eval)

  !in the case of minimal basis allocate now the transposed wavefunction
  !otherwise do it only in parallel
  if (minimal .or. nproc > 1) then
        allocate(psit(orbs%npsidim+ndebug),stat=i_stat)
        call memocc(i_stat,psit,'psit',subname)
  else
     psit => hpsi
  end if
  
  !allocate the pointer for virtual orbitals
  if(present(nvirte) .and. present(psivirt) .and. nvirte > 0) then
     allocate(psivirt(orbsv%npsidim+ndebug),stat=i_stat)
     call memocc(i_stat,psivirt,'psivirt',subname)
  end if

  if (iproc.eq.0) write(*,'(1x,a)',advance='no')'Building orthogonal Wavefunctions...'
  nvctr=wfd%nvctr_c+7*wfd%nvctr_f

  if (.not. present(nvirte)) then
     call build_eigenvectors(orbs%norbu,orbs%norbd,orbs%norb,norbtot,nvctrp,nvctr,&
          natsceff,nspin,orbs%nspinor,ndim_hamovr,norbgrp,hamovr,psi,psit)
  else
     call build_eigenvectors(orbs%norbu,orbs%norbd,orbs%norb,norbtot,nvctrp,nvctr,&
          natsceff,nspin,orbs%nspinor,ndim_hamovr,norbgrp,hamovr,psi,psit,nvirte,psivirt)
  end if
  
  !if(nproc==1.and.nspinor==4) call psitransspi(nvctrp,norbu+norbd,psit,.false.)
     
  i_all=-product(shape(hamovr))*kind(hamovr)
  deallocate(hamovr,stat=i_stat)
  call memocc(i_stat,i_all,'hamovr',subname)
  i_all=-product(shape(norbgrp))*kind(norbgrp)
  deallocate(norbgrp,stat=i_stat)
  call memocc(i_stat,i_all,'norbgrp',subname)

  if (minimal) then
     !deallocate the old psi
     i_all=-product(shape(psi))*kind(psi)
     deallocate(psi,stat=i_stat)
     call memocc(i_stat,i_all,'psi',subname)
  else if (nproc == 1) then
     !reverse objects for the normal diagonalisation in serial
     !at this stage hpsi is the eigenvectors and psi is the old wavefunction
     !this will restore the correct identification
     nullify(hpsi)
     hpsi => psi
!     if(nspinor==4) call psitransspi(nvctrp,norb,psit,.false.) 
    nullify(psi)
     psi => psit
  end if

  !orthogonalise the orbitals in the case of semi-core atoms
  if (norbsc > 0) then
!!$     if(nspin==1) then
!!$        call orthon_p(iproc,nproc,norb,nvctrp,wfd%nvctr_c+7*wfd%nvctr_f,psit,nspinor) 
!!$     else
     call orthon_p(iproc,nproc,orbs%norbu,nvctrp,wfd%nvctr_c+7*wfd%nvctr_f,psit,&
          orbs%nspinor) 
     if(orbs%norbd > 0) then
        call orthon_p(iproc,nproc,orbs%norbd,nvctrp,wfd%nvctr_c+7*wfd%nvctr_f,&
             psit(1+nvctrp*norbu),orbs%nspinor) 
!!$        end if
     end if
  end if


  if (minimal) then
     allocate(hpsi(orbs%npsidim+ndebug),stat=i_stat)
     call memocc(i_stat,hpsi,'hpsi',subname)
!     hpsi=0.0d0
     if (nproc > 1) then
        !allocate the direct wavefunction
        allocate(psi(orbs%npsidim+ndebug),stat=i_stat)
        call memocc(i_stat,psi,'psi',subname)
     else
        psi => psit
     end if
  end if

  !this untranspose also the wavefunctions 
  call untranspose_v(iproc,nproc,orbs%norbp,orbs%nspinor,wfd,nvctrp,comms,&
       psit,work=hpsi,outadd=psi(1))

  if (nproc == 1) then
     nullify(psit)
  end if

end subroutine DiagHam
!!***

subroutine overlap_matrices(norbe,nvctrp,natsc,nspin,ndim_hamovr,norbsc_arr,hamovr,psi,hpsi)
  use module_base
  implicit none
  integer, intent(in) :: norbe,nvctrp,natsc,ndim_hamovr,nspin
  integer, dimension(natsc+1,nspin), intent(in) :: norbsc_arr
  real(wp), dimension(nspin*ndim_hamovr,2), intent(out) :: hamovr
  real(wp), dimension(nvctrp,norbe), intent(in) :: psi,hpsi
  !local variables
  integer :: iorbst,imatrst,norbi,i,ispin,j
  real(kind=4) :: t0,t1

  !calculate the overlap matrix for each group of the semicore atoms
  !       hamovr(jorb,iorb,3)=+psit(k,jorb)*hpsit(k,iorb)
  !       hamovr(jorb,iorb,4)=+psit(k,jorb)* psit(k,iorb)
  iorbst=1
  imatrst=1
  do ispin=1,nspin !this construct assumes that the semicore is identical for both the spins
     do i=1,natsc+1
        norbi=norbsc_arr(i,ispin)
        call gemm('T','N',norbi,norbi,nvctrp,1.0_wp,psi(1,iorbst),nvctrp,hpsi(1,iorbst),nvctrp,&
             0.0_wp,hamovr(imatrst,1),norbi)
        call gemm('T','N',norbi,norbi,nvctrp,1.0_wp,psi(1,iorbst),nvctrp,psi(1,iorbst),nvctrp,&
             0.0_wp,hamovr(imatrst,2),norbi)
        iorbst=iorbst+norbi
        imatrst=imatrst+norbi**2
     end do
  end do

end subroutine overlap_matrices

subroutine solve_eigensystem(iproc,norb,norbu,norbd,norbi_max,ndim_hamovr,natsc,nspin,etol,&
     norbsc_arr,hamovr,eval)
  use module_base
  implicit none
  integer, intent(in) :: iproc,norb,norbi_max,ndim_hamovr,natsc,nspin,norbu,norbd
  integer, dimension(natsc+1,nspin), intent(in) :: norbsc_arr
  real(gp), intent(in) :: etol
  real(wp), dimension(nspin*ndim_hamovr,2), intent(inout) :: hamovr
  real(wp), dimension(norb), intent(out) :: eval
  !local variables
  character(len=*), parameter :: subname='solve_eigensystem'
  character(len=64) :: message
  integer :: iorbst,imatrst,norbi,n_lp,info,i_all,i_stat,iorb,i,ndegen,nwrtmsg,jorb,istart,norbj
  integer :: jjorb,jiorb
  real(wp), dimension(2) :: preval
  real(wp), dimension(:), allocatable :: work_lp,evale

  !find the eigenfunctions for each group
  n_lp=max(10,4*norbi_max)
  allocate(work_lp(n_lp+ndebug),stat=i_stat)
  call memocc(i_stat,work_lp,'work_lp',subname)
  allocate(evale(nspin*norbi_max+ndebug),stat=i_stat)
  call memocc(i_stat,evale,'evale',subname)

  if (iproc.eq.0) write(*,'(1x,a)')'Linear Algebra...'

  nwrtmsg=0
  ndegen=0

  preval=0.0_wp
  iorbst=1
  imatrst=1
  do i=1,natsc+1
     norbi=norbsc_arr(i,1)

!!$     if (iproc == 0) then
!!$        !write the matrices on a file
!!$        open(31+2*(i-1))
!!$        do jjorb=1,norbi
!!$           write(31+2*(i-1),'(2000(1pe10.2))')&
!!$                (hamovr(imatrst-1+jiorb+(jjorb-1)*norbi,1),jiorb=1,norbi)
!!$        end do
!!$        close(31+2*(i-1))
!!$        open(32+2*(i-1))
!!$        do jjorb=1,norbi
!!$           write(32+2*(i-1),'(2000(1pe10.2))')&
!!$                (hamovr(imatrst-1+jiorb+(jjorb-1)*norbi,2),jiorb=1,norbi)
!!$        end do
!!$        close(32+2*(i-1))
!!$
!!$     end if

     !write(11,*)hamovr(:,1:2)

     call sygv(1,'V','U',norbi,hamovr(imatrst,1),norbi,hamovr(imatrst,2),&
          norbi,evale(1),work_lp(1),n_lp,info)
     if (info.ne.0) write(*,*) 'SYGV ERROR',info,i,natsc+1

     !do the diagonalisation separately in case of spin polarization     
     if (nspin==2) then
        norbj=norbsc_arr(i,2)
        call sygv(1,'V','U',norbj,hamovr(imatrst+ndim_hamovr,1),&
             norbj,hamovr(imatrst+ndim_hamovr,2),norbj,evale(norbi+1),work_lp(1),n_lp,info)
        if (info.ne.0) write(*,*) 'SYGV ERROR',info,i,natsc+1
     end if

!!$     if (iproc == 0) then
!!$        !write the matrices on a file
!!$        open(12)
!!$        do jjorb=1,norbi
!!$           do jiorb=1,norbi
!!$              write(12,'(1x,2(i0,1x),2(1pe24.17,1x))')jjorb,jiorb,&
!!$                   hamovr(jjorb+norbi*(jiorb-1),1),hamovr(jjorb+norbi*(jiorb-1),2)
!!$           end do
!!$        end do
!!$        close(12)
!!$        !open(33+2*(i-1))
!!$        !write(33+2*(i-1),'(2000(1pe10.2))')&
!!$        !        (hamovr(imatrst-1+jiorb+(jjorb-1)*norbi,1),jiorb=1,norbi)
!!$        !end do
!!$        !close(33+2*(i-1))
!!$        !open(34+2*(i-1))
!!$        !do jjorb=1,norbi
!!$        !   write(34+2*(i-1),'(2000(1pe10.2))')&
!!$        !        (hamovr(imatrst-1+jiorb+(jjorb-1)*norbi,2),jiorb=1,norbi)
!!$        !end do
!!$        !close(34+2*(i-1))
!!$
!!$     end if

     !writing rules, control if the last eigenvector is degenerate
     !do this for each spin
     !for each spin it is supposed that only the last group is not completely passed
     !and also that the components of each of the group but the last are the same for up and 
     !down polarisation. Do not work properly in the other cases
     do iorb=1,norbi
        if (nspin==1) then
           if (nwrtmsg==1) then
              if (abs(evale(iorb)-preval(1)) <= etol) then
                 !degeneracy found
                 message='  <- found degeneracy'
                 ndegen=ndegen+1
              else
                 nwrtmsg=0
              end if
           end if
           if (iorb+iorbst-1 == norb) then
              nwrtmsg=1
              message=' <- Last eigenvalue for input wavefunctions'
              preval(1)=evale(iorb)
           end if
           if (iproc.eq.0) then
              if (nwrtmsg == 1) then
                 write(*,'(1x,a,i0,a,1x,1pe21.14,a)') &
                      'evale(',iorb+iorbst-1,')=',evale(iorb),trim(message)
              else
                 write(*,'(1x,a,i0,a,1x,1pe21.14)') &
                      'evale(',iorb+iorbst-1,')=',evale(iorb)
              end if
           end if
        else
           if (nwrtmsg==1) then
              if (abs(evale(iorb)-preval(1)) <= etol .and. &
                   abs(evale(iorb+norbi)-preval(2)) <= etol) then
                 !degeneracy found
                 message='  <-deg->  '
                 !ndegen=ndegen+1 removed, only for non magnetized cases
              else if (abs(evale(iorb)-preval(1)) <= etol) then
                 !degeneracy found
                 message='  <-deg    '
              else if (abs(evale(iorb+norbi)-preval(2)) <= etol) then
                 !degeneracy found
                 message='    deg->  '
              else
                 nwrtmsg=0
              end if
           end if
           if (iorb+iorbst-1 == norbu .and. iorb+iorbst-1 == norbd) then
              nwrtmsg=1
              message='  <-Last-> ' 
              preval(1)=evale(iorb)
              preval(2)=evale(iorb+norbi)
           else if (iorb+iorbst-1 == norbu) then
              nwrtmsg=1
              message='  <-Last   '
              preval(1)=evale(iorb)
           else if (iorb+iorbst-1 == norbd) then
              nwrtmsg=1
              message='    Last-> '
              preval(2)=evale(iorb+norbi)
           end if
           if (iproc == 0) then
              if (nwrtmsg==1) then
                 write(*,'(1x,a,i4,a,1x,1pe21.14,a12,a,i4,a,1x,1pe21.14)') &
                      'evale(',iorb+iorbst-1,',u)=',evale(iorb),message,&
                      'evale(',iorb+iorbst-1,',d)=',evale(iorb+norbi)
              else
                 write(*,'(1x,a,i4,a,1x,1pe21.14,12x,a,i4,a,1x,1pe21.14)') &
                      'evale(',iorb+iorbst-1,',u)=',evale(iorb),&
                      'evale(',iorb+iorbst-1,',d)=',evale(iorb+norbi)
              end if
           end if
        end if
     end do
     if (nspin==1) then
        do iorb=iorbst,min(norbi+iorbst-1,norb)
           eval(iorb)=evale(iorb-iorbst+1)
        end do
     else
        do iorb=iorbst,min(norbi+iorbst-1,norbu)
           eval(iorb)=evale(iorb-iorbst+1)
        end do
        do iorb=iorbst,min(norbi+iorbst-1,norbd)
           eval(iorb+norbu)=evale(iorb-iorbst+1+norbi)
        end do
     end if
     iorbst=iorbst+norbi
     imatrst=imatrst+norbi**2
  end do

  i_all=-product(shape(work_lp))*kind(work_lp)
  deallocate(work_lp,stat=i_stat)
  call memocc(i_stat,i_all,'work_lp',subname)
  i_all=-product(shape(evale))*kind(evale)
  deallocate(evale,stat=i_stat)
  call memocc(i_stat,i_all,'evale',subname)

end subroutine solve_eigensystem

subroutine build_eigenvectors(norbu,norbd,norb,norbe,nvctrp,natsc,nspin,nspinor,&
     ndim_hamovr,norbsc_arr,hamovr,psi,ppsit,nvirte,psivirt)
  use module_base
  implicit none
  !Arguments
  integer, intent(in) :: norbu,norbd,norb,norbe,nvctrp,natsc,nspin,nspinor,ndim_hamovr
  integer, dimension(natsc+1,nspin), intent(in) :: norbsc_arr
  real(wp), dimension(nspin*ndim_hamovr), intent(in) :: hamovr
  real(wp), dimension(nvctrp,norbe), intent(in) :: psi
  real(wp), dimension(nvctrp*nspinor,norb), intent(out) :: ppsit
  integer, intent(in), optional :: nvirte
  real(wp), dimension(:), pointer, optional :: psivirt
  !Local variables
  character(len=*), parameter :: subname='build_eigenvectors'
  integer, parameter :: iunit=1978
  integer :: ispin,iorbst,iorbst2,imatrst,norbsc,norbi,norbj,i,iorb,i_stat,i_all
  logical :: exists
  real(gp) :: mx,my,mz,mnorm,fac,ma,mb,mc,md
  real(wp), dimension(:,:), allocatable :: tpsi

  !perform the vector-matrix multiplication for building the input wavefunctions
  ! ppsit(k,iorb)=+psit(k,jorb)*hamovr(jorb,iorb,1)
  !!     iorbst=1
  !!     imatrst=1
  !!     do i=1,natsc
  !!        norbi=norbsc_arr(i)
  !!        call DGEMM('N','N',nvctrp,norbi,norbi,1.d0,psi(1,iorbst),nvctrp,&
  !!             hamovr(imatrst,1),norbi,0.d0,ppsit(1,iorbst),nvctrp)
  !!        iorbst=iorbst+norbi
  !!        imatrst=imatrst+norbi**2
  !!     end do
  !!     norbi=norbsc_arr(natsc+1)
  !!     norbj=norb-norbsc
  !!     call DGEMM('N','N',nvctrp,norbj,norbi,1.d0,psi(1,iorbst),nvctrp,&
  !!          hamovr(imatrst,1),norbi,0.d0,ppsit(1,iorbst),nvctrp)

  !ppsi(k,iorb)=+psi(k,jorb)*hamovr(jorb,iorb,1)

  !allocate the pointer for virtual orbitals

  if(nspinor==1) then
     iorbst=1
     iorbst2=1
     imatrst=1
     do ispin=1,nspin
        norbsc=0
        do i=1,natsc
           norbi=norbsc_arr(i,ispin)
           norbsc=norbsc+norbi
           call gemm('N','N',nvctrp,norbi,norbi,1.0_wp,psi(1,iorbst),nvctrp,&
                hamovr(imatrst),norbi,0.0_wp,ppsit(1,iorbst2),nvctrp)
           iorbst=iorbst+norbi
           iorbst2=iorbst2+norbi
           imatrst=imatrst+norbi**2
        end do
        norbi=norbsc_arr(natsc+1,ispin)
        if(ispin==1) norbj=norbu-norbsc
        if(ispin==2) norbj=norbd-norbsc
        !        write(*,'(1x,a,5i4)') "DIMS:",norbi,norbj,iorbst,imatrst
        !        norbj=norb-norbsc
        if(norbj>0) then
           call gemm('N','N',nvctrp,norbj,norbi,1.0_wp,psi(1,iorbst),nvctrp,&
                hamovr(imatrst),norbi,0.0_wp,ppsit(1,iorbst2),nvctrp)
        end if

        !now store the input wavefunctions for the Davidson treatment
        !we take the rest of the orbitals which are not assigned
        !from the group of non-semicore orbitals
        !the results are orthogonal with each other by construction
        !in the case of semicore atomes the orthogonality is not guaranteed
        if (present(nvirte) .and. nvirte >0) then
           call gemm('N','N',nvctrp,nvirte,norbi,1.0_wp,psi(1,iorbst),nvctrp,&
                hamovr(imatrst+norbi*norbj),norbi,0.0_wp,psivirt(1),nvctrp)
        end if
        iorbst=norbi+norbsc+1 !this is equal to norbe+1
        iorbst2=norbu+1
        imatrst=ndim_hamovr+1
     end do

  else
     allocate(tpsi(nvctrp,norbe+ndebug),stat=i_stat)
     call memocc(i_stat,tpsi,'tpsi',subname)
     iorbst=1
     iorbst2=1
     imatrst=1
     do ispin=1,nspin
        norbsc=0
        do i=1,natsc
           norbi=norbsc_arr(i,ispin)
           norbsc=norbsc+norbi
           call gemm('N','N',nvctrp,norbi,norbi,1.0_wp,psi(1,iorbst),nvctrp,&
                hamovr(imatrst),norbi,0.0_wp,tpsi(1,iorbst2),nvctrp)
           iorbst=iorbst+norbi
           iorbst2=iorbst2+norbi
           imatrst=imatrst+norbi**2
        end do
        norbi=norbsc_arr(natsc+1,ispin)
        if(ispin==1) norbj=norbu-norbsc
        if(ispin==2) norbj=norbd-norbsc
        !        write(*,'(1x,a,5i4)') "DIMS:",norbi,norbj,iorbst,imatrst
        !        norbj=norb-norbsc
        if(norbj>0) then
           call gemm('N','N',nvctrp,norbj,norbi,1.0_wp,psi(1,iorbst),nvctrp,&
                hamovr(imatrst),norbi,0.0_wp,tpsi(1,iorbst2),nvctrp)
        end if
        iorbst=norbi+norbsc+1 !this is equal to norbe+1
        iorbst2=norbu+1
        imatrst=ndim_hamovr+1
     end do
     !here we should put razero
     ppsit=0.0_wp
     inquire(file='moments',exist=exists)
     if (.not.exists) then
        stop 'The file "moments does not exist!'
     endif
     open(unit=iunit,file='moments',form='formatted',action='read',status='old')
     fac=0.5_gp
     do iorb=1,norbu+norbd
        read(unit=iunit,fmt=*,iostat=i_stat) mx,my,mz
        if (i_stat /= 0) then
           write(unit=*,fmt='(a,i0,a,i0,a)') &
                'The file "moments" has the line ',iorb,&
                ' which have not 3 numbers for the orbital ',iorb,'.'
           stop 'The file "moments" is not correct!'
        end if
        mnorm=sqrt(mx**2+my**2+mz**2)
        mx=mx/mnorm
        my=my/mnorm
        mz=mz/mnorm

        ma=0.0_gp
        mb=0.0_gp
        mc=0.0_gp
        md=0.0_gp

        if(mz > 0.0_gp .and. iorb<=norbu) then 
           ma=ma+mz
        else
           mc=mc+abs(mz)
        end if
        if(mx > 0.0_gp .and. iorb<=norbu) then 
           ma=ma+fac*mx
           mb=mb+fac*mx
           mc=mc+fac*mx
           md=md+fac*mx
        else
           ma=ma-fac*abs(mx)
           mb=mb-fac*abs(mx)
           mc=mc+fac*abs(mx)
           md=md+fac*abs(mx)
        end if
        if(my > 0.0_gp .and. iorb<=norbu) then 
           ma=ma+fac*my
           mb=mb-fac*my
           mc=mc+fac*my
           md=md+fac*my
        else
           ma=ma-fac*abs(my)
           mb=mb+fac*abs(my)
           mc=mc+fac*abs(my)
           md=md+fac*abs(my)
        end if
        if(mx==0.0_gp .and. my==0.0_gp .and. mz==0.0_gp) then
           ma=1.0_gp/sqrt(2.0_gp)
           mb=0.0_gp
           mc=1.0_gp/sqrt(2.0_gp)
           md=0.0_gp
        end if
        do i=1,nvctrp
           ppsit(2*i-1,iorb)=real(ma,wp)*tpsi(i,iorb)
           ppsit(2*i,iorb)=real(mb,wp)*tpsi(i,iorb)
           ppsit(2*i+2*nvctrp-1,iorb)=real(mc,wp)*tpsi(i,iorb)
           ppsit(2*i+2*nvctrp,iorb)=real(md,wp)*tpsi(i,iorb)
        end do
     end do
     close(unit=iunit)
    
     i_all=-product(shape(tpsi))*kind(tpsi)
     deallocate(tpsi,stat=i_stat)
     call memocc(i_stat,i_all,'tpsi',subname)
  end if

end subroutine build_eigenvectors

!  call psitospi(iproc,nproc,norbe,norbep,norbsc,nat,&
!       wfd%nvctr_c,wfd%nvctr_f,at%iatype,at%ntypes,&
!       at%iasctype,at%natsc,at%natpol,nspin,spinsgne,otoa,psi)
! Reads magnetic moments from file ('moments') and transforms the
! atomic orbitals to spinors 
! warning: Does currently not work for mx<0
!
subroutine psitospi(iproc,nproc,norbe,norbep,norbsc,&
     & nvctr_c,nvctr_f,nat,iatype,ntypes, &
     iasctype,natsc,natpol,nspin,spinsgne,otoa,psi)
  use module_base
  implicit none
  integer, intent(in) :: norbe,norbep,iproc,nproc,nat
  integer, intent(in) :: nvctr_c,nvctr_f
  integer, intent(in) :: ntypes
  integer, intent(in) :: norbsc,natsc,nspin
  integer, dimension(ntypes), intent(in) :: iasctype
  integer, dimension(norbep), intent(in) :: otoa
  integer, dimension(nat), intent(in) :: iatype,natpol
  integer, dimension(norbe*nspin), intent(in) :: spinsgne
  real(kind=8), dimension(nvctr_c+7*nvctr_f,4*norbep), intent(out) :: psi
  !local variables
  character(len=*), parameter :: subname='psitospi'
  logical :: myorbital,polarised
  integer :: iatsc,i_all,i_stat,ispin,ipolres,ipolorb,nvctr
  integer :: iorb,jorb,iat,ity,i,ictot,inl,l,m,nctot,nterm
  real(kind=8) :: facu,facd
  real(kind=8) :: mx,my,mz,mnorm,fac
  real(kind=8), dimension(:,:), allocatable :: mom,psi_o
  logical, dimension(4) :: semicore
  integer, dimension(2) :: iorbsc,iorbv

  !initialise the orbital counters
  iorbsc(1)=0
  iorbv(1)=norbsc
  !used in case of spin-polarisation, ignored otherwise
  iorbsc(2)=norbe
  iorbv(2)=norbsc+norbe


  if (iproc ==0) then
     write(*,'(1x,a)',advance='no')'Transforming AIO to spinors...'
  end if
  
  nvctr=nvctr_c+7*nvctr_f

  allocate(mom(3,nat+ndebug),stat=i_stat)
  call memocc(i_stat,mom,'mom',subname)

  open(unit=1978,file='moments')
  do i=1,nat
     read(1978,*) mx,my,mz
     mnorm=sqrt(mx**2+my**2+mz**2)
     mom(1,iat)=mx/mnorm
     mom(2,iat)=my/mnorm
     mom(3,iat)=mz/mnorm
  end do
  close(1978)
  fac=0.5d0
  do iorb=norbep*nproc,1,-1
     jorb=iorb-iproc*norbep
!     print *,'Kolla', shape(psi),4*iorb,shape(spinsgne),iorb
     if (myorbital(iorb,nspin*norbe,iproc,nproc)) then
        mx=mom(1,otoa(iorb))
        my=mom(2,otoa(iorb))
        mz=mom(3,otoa(iorb))
        if(spinsgne(jorb)>0.0d0) then
           do i=1,nvctr
              psi(i,iorb*4-3) = (mz+fac*(my+mx))*psi(i,iorb)
              psi(i,iorb*4-2) = fac*(my-mx)*psi(i,iorb)
              psi(i,iorb*4-1) = (fac*(mx-my))*psi(i,iorb)
              psi(i,iorb*4)   = fac*(my-mx)*psi(i,iorb)
           end do
        else
           do i=1,nvctr
              psi(i,iorb*4-3) = (fac*(mx+my))*psi(i,iorb)
              psi(i,iorb*4-2) = -fac*(my+mx)*psi(i,iorb)
              psi(i,iorb*4-1) = -(mz+fac*(my+mx))*psi(i,iorb)
              psi(i,iorb*4)   = -fac*(my-mx)*psi(i,iorb)
           end do
        end if
     end if
!     print *,'OtoA',(otoa(iorb),iorb=1,norbe)

  end do
     i_all=-product(shape(mom))*kind(mom)
     deallocate(mom,stat=i_stat)
     call memocc(i_stat,i_all,'mom',subname)

  if (iproc ==0) then
     write(*,'(1x,a)')'done.'
  end if

END SUBROUTINE psitospi
 
