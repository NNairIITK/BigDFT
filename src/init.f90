subroutine createWavefunctionsDescriptors(iproc,nproc,idsx,n1,n2,n3,output_grid,&
     hgrid,nat,ntypes,iatype,atomnames,alat1,alat2,alat3,rxyz,radii_cf,crmult,frmult,&
     wfd,nvctrp,norb,norbp,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,bounds)
  !calculates the descriptor arrays keyg and keyv as well as nseg_c,nseg_f,nvctr_c,nvctr_f,nvctrp
  !calculates also the bounds arrays needed for convolutions

  use module_types

  implicit none
  !Arguments
  integer, intent(in) :: iproc,nproc,idsx,n1,n2,n3,nat,ntypes,norb,norbp
  integer, intent(out) :: nvctrp
  logical, intent(in) :: output_grid
  integer, intent(in) :: iatype(nat)
  real(kind=8), intent(in) :: hgrid,crmult,frmult,alat1,alat2,alat3
  real(kind=8) :: rxyz(3, nat), radii_cf(ntypes, 2)
  character(len=20), intent(in) :: atomnames(100)
  integer,intent(in):: nfl1,nfu1,nfl2,nfu2,nfl3,nfu3
  type(wavefunctions_descriptors) , intent(out) :: wfd
  !boundary arrays
  type(convolutions_bounds), intent(out) :: bounds

  !Local variables
  real(kind=8), parameter :: eps_mach=1.d-12
  integer :: iat,i1,i2,i3,norbme,norbyou,jpst,jproc,i_all,i_stat
  real(kind=8) :: tt
  logical, allocatable :: logrid_c(:,:,:), logrid_f(:,:,:)

  !allocate kinetic bounds
  allocate(bounds%kb%ibyz_c(2,0:n2,0:n3),stat=i_stat)
  call memocc(i_stat,product(shape(bounds%kb%ibyz_c))*kind(bounds%kb%ibyz_c),'ibyz_c','crtwvfnctsdescriptors')
  allocate(bounds%kb%ibxz_c(2,0:n1,0:n3),stat=i_stat)
  call memocc(i_stat,product(shape(bounds%kb%ibxz_c))*kind(bounds%kb%ibxz_c),'ibxz_c','crtwvfnctsdescriptors')
  allocate(bounds%kb%ibxy_c(2,0:n1,0:n2),stat=i_stat)
  call memocc(i_stat,product(shape(bounds%kb%ibxy_c))*kind(bounds%kb%ibxy_c),'ibxy_c','crtwvfnctsdescriptors')
  allocate(bounds%kb%ibyz_f(2,0:n2,0:n3),stat=i_stat)
  call memocc(i_stat,product(shape(bounds%kb%ibyz_f))*kind(bounds%kb%ibyz_f),'ibyz_f','crtwvfnctsdescriptors')
  allocate(bounds%kb%ibxz_f(2,0:n1,0:n3),stat=i_stat)
  call memocc(i_stat,product(shape(bounds%kb%ibxz_f))*kind(bounds%kb%ibxz_f),'ibxz_f','crtwvfnctsdescriptors')
  allocate(bounds%kb%ibxy_f(2,0:n1,0:n2),stat=i_stat)
  call memocc(i_stat,product(shape(bounds%kb%ibxy_f))*kind(bounds%kb%ibxy_f),'ibxy_f','crtwvfnctsdescriptors')

  if (iproc.eq.0) then
     write(*,'(1x,a)')&
          '------------------------------------------------- Wavefunctions Descriptors Creation'
  end if

  ! Create the file grid.ascii to visualize the grid of functions
  if (iproc.eq.0 .and. output_grid) then
     open(unit=22,file='grid.ascii',status='unknown')
     write(22,*) nat
     write(22,*) alat1,' 0. ',alat2
     write(22,*) ' 0. ',' 0. ',alat3
     do iat=1,nat
        write(22,'(3(1x,e12.5),3x,a20)') &
             rxyz(1,iat),rxyz(2,iat),rxyz(3,iat),atomnames(iatype(iat))
     enddo
  endif

  ! determine localization region for all orbitals, but do not yet fill the descriptor arrays
  allocate(logrid_c(0:n1,0:n2,0:n3),stat=i_stat)
  call memocc(i_stat,product(shape(logrid_c))*kind(logrid_c),'logrid_c','crtwvfnctsdescriptors')
  allocate(logrid_f(0:n1,0:n2,0:n3),stat=i_stat)
  call memocc(i_stat,product(shape(logrid_f))*kind(logrid_f),'logrid_f','crtwvfnctsdescriptors')

  ! coarse grid quantities
  call fill_logrid(n1,n2,n3,0,n1,0,n2,0,n3,0,nat,ntypes,iatype,rxyz, & 
       radii_cf(1,1),crmult,hgrid,logrid_c)
  if (iproc.eq.0 .and. output_grid) then
     do i3=0,n3  
        do i2=0,n2  
           do i1=0,n1
              if (logrid_c(i1,i2,i3))&
                   write(22,'(3(1x,e10.3),1x,a4)') &
                   real(i1,kind=8)*hgrid,real(2,kind=8)*hgrid,real(i3,kind=8)*hgrid,'  g '
           enddo
        enddo
     end do
  endif
  call num_segkeys(n1,n2,n3,0,n1,0,n2,0,n3,logrid_c,wfd%nseg_c,wfd%nvctr_c)
  if (iproc.eq.0) write(*,'(2(1x,a,i10))') &
       'Coarse resolution grid: Number of segments= ',wfd%nseg_c,'points=',wfd%nvctr_c
  call make_bounds(n1,n2,n3,logrid_c,bounds%kb%ibyz_c,bounds%kb%ibxz_c,bounds%kb%ibxy_c)

  ! fine grid quantities
  call fill_logrid(n1,n2,n3,0,n1,0,n2,0,n3,0,nat,ntypes,iatype,rxyz, & 
       radii_cf(1,2),frmult,hgrid,logrid_f)
  if (iproc.eq.0 .and. output_grid) then
     do i3=0,n3 
        do i2=0,n2 
           do i1=0,n1
              if (logrid_f(i1,i2,i3))&
                   write(22,'(3(1x,e10.3),1x,a4)') &
                   real(i1,kind=8)*hgrid,real(i2,kind=8)*hgrid,real(i3,kind=8)*hgrid,'  G '
           enddo
        enddo
     enddo
  endif
  call num_segkeys(n1,n2,n3,0,n1,0,n2,0,n3,logrid_f,wfd%nseg_f,wfd%nvctr_f)
  if (iproc.eq.0) write(*,'(2(1x,a,i10))') &
       '  Fine resolution grid: Number of segments= ',wfd%nseg_f,'points=',wfd%nvctr_f
  call make_bounds(n1,n2,n3,logrid_f,bounds%kb%ibyz_f,bounds%kb%ibxz_f,bounds%kb%ibxy_f)

  if (iproc.eq.0 .and. output_grid) close(22)

  ! allocations for arrays holding the wavefunctions and their data descriptors
  call allocate_wfd(wfd,'crtwvfnctsdescriptors')

  ! now fill the wavefunction descriptor arrays
  ! coarse grid quantities
  call segkeys(n1,n2,n3,0,n1,0,n2,0,n3,logrid_c,wfd%nseg_c,wfd%keyg(1,1),wfd%keyv(1))

  ! fine grid quantities
  call segkeys(n1,n2,n3,0,n1,0,n2,0,n3,logrid_f,wfd%nseg_f,wfd%keyg(1,wfd%nseg_c+1), &
       & wfd%keyv(wfd%nseg_c+1))

  i_all=-product(shape(logrid_c))*kind(logrid_c)
  deallocate(logrid_c,stat=i_stat)
  call memocc(i_stat,i_all,'logrid_c','crtwvfnctsdescriptors')
  i_all=-product(shape(logrid_f))*kind(logrid_f)
  deallocate(logrid_f,stat=i_stat)
  call memocc(i_stat,i_all,'logrid_f','crtwvfnctsdescriptors')

  ! allocate wavefunction arrays
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
       nvctrp*nproc*8

  !allocate grow, shrink and real bounds
  allocate(bounds%gb%ibzxx_c(2,0:n3,-14:2*n1+16),stat=i_stat)
  call memocc(i_stat,product(shape(bounds%gb%ibzxx_c))*kind(bounds%gb%ibzxx_c),'ibzxx_c','crtwvfnctsdescriptors')
  allocate(bounds%gb%ibxxyy_c(2,-14:2*n1+16,-14:2*n2+16),stat=i_stat)
  call memocc(i_stat,product(shape(bounds%gb%ibxxyy_c))*kind(bounds%gb%ibxxyy_c),'ibxxyy_c','crtwvfnctsdescriptors')
  allocate(bounds%gb%ibyz_ff(2,nfl2:nfu2,nfl3:nfu3),stat=i_stat)
  call memocc(i_stat,product(shape(bounds%gb%ibyz_ff))*kind(bounds%gb%ibyz_ff),'ibyz_ff','crtwvfnctsdescriptors')
  allocate(bounds%gb%ibzxx_f(2,nfl3:nfu3,2*nfl1-14:2*nfu1+16),stat=i_stat)
  call memocc(i_stat,product(shape(bounds%gb%ibzxx_f))*kind(bounds%gb%ibzxx_f),'ibzxx_f','crtwvfnctsdescriptors')
  allocate(bounds%gb%ibxxyy_f(2,2*nfl1-14:2*nfu1+16,2*nfl2-14:2*nfu2+16),stat=i_stat)
  call memocc(i_stat,product(shape(bounds%gb%ibxxyy_f))*kind(bounds%gb%ibxxyy_f),'ibxxyy_f','crtwvfnctsdescriptors')

  allocate(bounds%sb%ibzzx_c(2,-14:2*n3+16,0:n1),stat=i_stat)
  call memocc(i_stat,product(shape(bounds%sb%ibzzx_c))*kind(bounds%sb%ibzzx_c),'ibzzx_c','crtwvfnctsdescriptors')
  allocate(bounds%sb%ibyyzz_c(2,-14:2*n2+16,-14:2*n3+16),stat=i_stat)
  call memocc(i_stat,product(shape(bounds%sb%ibyyzz_c))*kind(bounds%sb%ibyyzz_c),'ibyyzz_c','crtwvfnctsdescriptors')
  allocate(bounds%sb%ibxy_ff(2,nfl1:nfu1,nfl2:nfu2),stat=i_stat)
  call memocc(i_stat,product(shape(bounds%sb%ibxy_ff))*kind(bounds%sb%ibxy_ff),'ibxy_ff','crtwvfnctsdescriptors')
  allocate(bounds%sb%ibzzx_f(2,-14+2*nfl3:2*nfu3+16,nfl1:nfu1),stat=i_stat)
  call memocc(i_stat,product(shape(bounds%sb%ibzzx_f))*kind(bounds%sb%ibzzx_f),'ibzzx_f','crtwvfnctsdescriptors')
  allocate(bounds%sb%ibyyzz_f(2,-14+2*nfl2:2*nfu2+16,-14+2*nfl3:2*nfu3+16),stat=i_stat)
  call memocc(i_stat,product(shape(bounds%sb%ibyyzz_f))*kind(bounds%sb%ibyyzz_f),'ibyyzz_f','crtwvfnctsdescriptors')

  allocate(bounds%ibyyzz_r(2,-14:2*n2+16,-14:2*n3+16),stat=i_stat)
  call memocc(i_stat,product(shape(bounds%ibyyzz_r))*kind(bounds%ibyyzz_r),'ibyyzz_r','crtwvfnctsdescriptors')


  call make_all_ib(n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,&
       bounds%kb%ibxy_c,bounds%sb%ibzzx_c,bounds%sb%ibyyzz_c,&
       bounds%kb%ibxy_f,bounds%sb%ibxy_ff,bounds%sb%ibzzx_f,bounds%sb%ibyyzz_f,&
       bounds%kb%ibyz_c,bounds%gb%ibzxx_c,bounds%gb%ibxxyy_c,&
       bounds%kb%ibyz_f,bounds%gb%ibyz_ff,bounds%gb%ibzxx_f,bounds%gb%ibxxyy_f,&
       bounds%ibyyzz_r)

  !***********************************************************************************************
END SUBROUTINE createWavefunctionsDescriptors

!pass to implicit none while inserting types on this routine
subroutine createProjectorsArrays(iproc,n1,n2,n3,rxyz,nat,ntypes,iatype,atomnames,&
     & psppar,npspcode,radii_cf,cpmult,fpmult,hgrid,nlpspd,proj)

  use module_types

  implicit none
  type(nonlocal_psp_descriptors), intent(out) :: nlpspd
  character(len=20), dimension(100),intent(in) :: atomnames
  integer, intent(in) :: iproc,n1,n2,n3,nat,ntypes
  real(kind=8), intent(in) :: cpmult,fpmult,hgrid
  integer, dimension(nat), intent(in) :: iatype
  integer, dimension(ntypes), intent(in) :: npspcode
  real(kind=8), dimension(3,nat), intent(in) :: rxyz
  real(kind=8), dimension(ntypes,2), intent(in) :: radii_cf
  real(kind=8), dimension(0:4,0:6,ntypes), intent(in) :: psppar
  real(kind=8), dimension(:), pointer :: proj
  !local variables
  integer :: nl1,nl2,nl3,nu1,nu2,nu3,mseg,mvctr,mproj,istart,istart_c,istart_f,mvctr_c,mvctr_f
  integer :: nl1_c,nl1_f,nl2_c,nl2_f,nl3_c,nl3_f,nu1_c,nu1_f,nu2_c,nu2_f,nu3_c,nu3_f
  integer :: iat,i_stat,i_all,nterm_max,i,l,m,iproj,ityp,nterm,iseg
  real(kind=8) :: fpi,factor,scpr,gau_a,rx,ry,rz
  logical, dimension(:,:,:), allocatable :: logrid
  real(kind=8), dimension(:), allocatable :: fac_arr
  integer, dimension(:), allocatable :: lx,ly,lz


  if (iproc.eq.0) then
     write(*,'(1x,a)')&
          '------------------------------------------------------------ PSP Projectors Creation'
     write(*,'(1x,a4,4x,a4,1x,a)')&
          'Atom','Name','Number of projectors'
  end if

  allocate(nlpspd%nseg_p(0:2*nat),stat=i_stat)
  call memocc(i_stat,product(shape(nlpspd%nseg_p))*kind(nlpspd%nseg_p),&
       'nseg_p','createprojectorsarrays')
  allocate(nlpspd%nvctr_p(0:2*nat),stat=i_stat)
  call memocc(i_stat,product(shape(nlpspd%nvctr_p))*kind(nlpspd%nvctr_p),&
       'nvctr_p','createprojectorsarrays')
  allocate(nlpspd%nboxp_c(2,3,nat),stat=i_stat)
  call memocc(i_stat,product(shape(nlpspd%nboxp_c))*kind(nlpspd%nboxp_c),&
       'nboxp_c','createprojectorsarrays')
  allocate(nlpspd%nboxp_f(2,3,nat),stat=i_stat)
  call memocc(i_stat,product(shape(nlpspd%nboxp_f))*kind(nlpspd%nboxp_f),&
       'nboxp_f','createprojectorsarrays')

  ! determine localization region for all projectors, but do not yet fill the descriptor arrays
  allocate(logrid(0:n1,0:n2,0:n3),stat=i_stat)
  call memocc(i_stat,product(shape(logrid))*kind(logrid),'logrid','createprojectorsarrays')

  nlpspd%nseg_p(0)=0 
  nlpspd%nvctr_p(0)=0 

  istart=1
  nlpspd%nproj=0
  do iat=1,nat

     call numb_proj(iatype(iat),ntypes,psppar,npspcode,mproj)
     if (mproj.ne.0) then 

        if (iproc.eq.0) write(*,'(1x,i4,2x,a6,1x,i20)')&
             iat,trim(atomnames(iatype(iat))),mproj


        !if (iproc.eq.0) write(*,'(1x,a,2(1x,i0))')&
        !     'projector descriptors for atom with mproj ',iat,mproj
        nlpspd%nproj=nlpspd%nproj+mproj

        ! coarse grid quantities
        call  pregion_size(rxyz(1,iat),radii_cf(1,2),cpmult,iatype(iat),ntypes, &
             hgrid,n1,n2,n3,nl1,nu1,nl2,nu2,nl3,nu3)
        !if (iproc.eq.0) write(*,'(a,6(i4))') 'coarse grid',nl1,nu1,nl2,nu2,nl3,nu3
        nlpspd%nboxp_c(1,1,iat)=nl1 
        nlpspd%nboxp_c(1,2,iat)=nl2 
        nlpspd%nboxp_c(1,3,iat)=nl3 

        nlpspd%nboxp_c(2,1,iat)=nu1
        nlpspd%nboxp_c(2,2,iat)=nu2
        nlpspd%nboxp_c(2,3,iat)=nu3

        call fill_logrid(n1,n2,n3,nl1,nu1,nl2,nu2,nl3,nu3,0,1,  &
             ntypes,iatype(iat),rxyz(1,iat),radii_cf(1,2),cpmult,hgrid,logrid)
        call num_segkeys(n1,n2,n3,nl1,nu1,nl2,nu2,nl3,nu3,logrid,mseg,mvctr)
        !if (iproc.eq.0) write(*,'(1x,a,2(1x,i0))') 'mseg,mvctr,coarse projectors ',mseg,mvctr

        nlpspd%nseg_p(2*iat-1)=nlpspd%nseg_p(2*iat-2) + mseg
        nlpspd%nvctr_p(2*iat-1)=nlpspd%nvctr_p(2*iat-2) + mvctr
        istart=istart+mvctr*mproj

        ! fine grid quantities
        call  pregion_size(rxyz(1,iat),radii_cf(1,2),fpmult,iatype(iat),ntypes, &
             hgrid,n1,n2,n3,nl1,nu1,nl2,nu2,nl3,nu3)
        !if (iproc.eq.0) write(*,'(a,6(i4))') 'fine   grid',nl1,nu1,nl2,nu2,nl3,nu3
        nlpspd%nboxp_f(1,1,iat)=nl1
        nlpspd%nboxp_f(1,2,iat)=nl2
        nlpspd%nboxp_f(1,3,iat)=nl3

        nlpspd%nboxp_f(2,1,iat)=nu1
        nlpspd%nboxp_f(2,2,iat)=nu2
        nlpspd%nboxp_f(2,3,iat)=nu3

        call fill_logrid(n1,n2,n3,nl1,nu1,nl2,nu2,nl3,nu3,0,1,  &
             ntypes,iatype(iat),rxyz(1,iat),radii_cf(1,2),fpmult,hgrid,logrid)
        call num_segkeys(n1,n2,n3,nl1,nu1,nl2,nu2,nl3,nu3,logrid,mseg,mvctr)
        !if (iproc.eq.0) write(*,'(1x,a,2(1x,i0))') 'mseg,mvctr, fine  projectors ',mseg,mvctr
        nlpspd%nseg_p(2*iat)=nlpspd%nseg_p(2*iat-1) + mseg
        nlpspd%nvctr_p(2*iat)=nlpspd%nvctr_p(2*iat-1) + mvctr

        istart=istart+7*mvctr*mproj

     else  !(atom has no nonlocal PSP, e.g. H)
        nlpspd%nseg_p(2*iat-1)=nlpspd%nseg_p(2*iat-2) 
        nlpspd%nvctr_p(2*iat-1)=nlpspd%nvctr_p(2*iat-2) 
        nlpspd%nseg_p(2*iat)=nlpspd%nseg_p(2*iat-1) 
        nlpspd%nvctr_p(2*iat)=nlpspd%nvctr_p(2*iat-1) 
     endif
  enddo

  if (iproc.eq.0) then
     write(*,'(28x,a)') '------'
     write(*,'(1x,a,i5)') 'Total number of projectors =',nlpspd%nproj
  end if

  ! allocations for arrays holding the projectors and their data descriptors
  allocate(nlpspd%keyg_p(2,nlpspd%nseg_p(2*nat)),stat=i_stat)
  call memocc(i_stat,product(shape(nlpspd%keyg_p))*kind(nlpspd%keyg_p),&
       'keyg_p','createprojectorsarrays')
  allocate(nlpspd%keyv_p(nlpspd%nseg_p(2*nat)),stat=i_stat)
  call memocc(i_stat,product(shape(nlpspd%keyv_p))*kind(nlpspd%keyv_p),&
       'keyv_p','createprojectorsarrays')
  nlpspd%nprojel=istart-1
  allocate(proj(nlpspd%nprojel),stat=i_stat)
  call memocc(i_stat,product(shape(proj))*kind(proj),'proj','createprojectorsarrays')


  ! After having determined the size of the projector descriptor arrays fill them
  istart_c=1
  do iat=1,nat
     call numb_proj(iatype(iat),ntypes,psppar,npspcode,mproj)
     if (mproj.ne.0) then 

        ! coarse grid quantities
        nl1=nlpspd%nboxp_c(1,1,iat) 
        nl2=nlpspd%nboxp_c(1,2,iat) 
        nl3=nlpspd%nboxp_c(1,3,iat) 

        nu1=nlpspd%nboxp_c(2,1,iat)
        nu2=nlpspd%nboxp_c(2,2,iat)
        nu3=nlpspd%nboxp_c(2,3,iat)
        call fill_logrid(n1,n2,n3,nl1,nu1,nl2,nu2,nl3,nu3,0,1,  &
             ntypes,iatype(iat),rxyz(1,iat),radii_cf(1,2),cpmult,hgrid,logrid)

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
        call fill_logrid(n1,n2,n3,nl1,nu1,nl2,nu2,nl3,nu3,0,1,  &
             ntypes,iatype(iat),rxyz(1,iat),radii_cf(1,2),fpmult,hgrid,logrid)
        iseg=nlpspd%nseg_p(2*iat-1)+1
        mseg=nlpspd%nseg_p(2*iat)-nlpspd%nseg_p(2*iat-1)
        call segkeys(n1,n2,n3,nl1,nu1,nl2,nu2,nl3,nu3,  & 
             logrid,mseg,nlpspd%keyg_p(1,iseg),nlpspd%keyv_p(iseg))

     endif
  enddo

  if (iproc.eq.0) write(*,'(1x,a)',advance='no') &
       'Calculating wavelets expansion of projectors...'
  !allocate these vectors up to the maximum size we can get
  nterm_max=10 !if GTH nterm_max=3
  allocate(fac_arr(nterm_max),stat=i_stat)
  call memocc(i_stat,product(shape(fac_arr))*kind(fac_arr),'fac_arr','createprojectorsarrays')
  allocate(lx(nterm_max),stat=i_stat)
  call memocc(i_stat,product(shape(lx))*kind(lx),'lx','createprojectorsarrays')
  allocate(ly(nterm_max),stat=i_stat)
  call memocc(i_stat,product(shape(ly))*kind(ly),'ly','createprojectorsarrays')
  allocate(lz(nterm_max),stat=i_stat)
  call memocc(i_stat,product(shape(lz))*kind(lz),'lz','createprojectorsarrays')

  iproj=0
  fpi=(4.d0*atan(1.d0))**(-.75d0)
  do iat=1,nat
     rx=rxyz(1,iat) 
     ry=rxyz(2,iat) 
     rz=rxyz(3,iat)
     ityp=iatype(iat)

     !decide the loop bounds
     do l=1,4 !generic case, also for HGHs (for GTH it will stop at l=2)
        do i=1,3 !generic case, also for HGHs (for GTH it will stop at i=2)
           if (psppar(l,i,ityp).ne.0.d0) then
              gau_a=psppar(l,0,ityp)
              factor=sqrt(2.d0)*fpi/(sqrt(gau_a)**(2*(l-1)+4*i-1))
              do m=1,2*l-1
                 mvctr_c=nlpspd%nvctr_p(2*iat-1)-nlpspd%nvctr_p(2*iat-2)
                 mvctr_f=nlpspd%nvctr_p(2*iat  )-nlpspd%nvctr_p(2*iat-1)
                 istart_f=istart_c+mvctr_c
                 nl1_c=nlpspd%nboxp_c(1,1,iat)
                 nl2_c=nlpspd%nboxp_c(1,2,iat)
                 nl3_c=nlpspd%nboxp_c(1,3,iat)
                 nl1_f=nlpspd%nboxp_f(1,1,iat)
                 nl2_f=nlpspd%nboxp_f(1,2,iat)
                 nl3_f=nlpspd%nboxp_f(1,3,iat)

                 nu1_c=nlpspd%nboxp_c(2,1,iat)
                 nu2_c=nlpspd%nboxp_c(2,2,iat)
                 nu3_c=nlpspd%nboxp_c(2,3,iat)
                 nu1_f=nlpspd%nboxp_f(2,1,iat)
                 nu2_f=nlpspd%nboxp_f(2,2,iat)
                 nu3_f=nlpspd%nboxp_f(2,3,iat)

                 call calc_coeff_proj(l,i,m,nterm_max,nterm,lx,ly,lz,fac_arr)

                 fac_arr(1:nterm)=factor*fac_arr(1:nterm)

                 call crtproj(iproc,nterm,n1,n2,n3,nl1_c,nu1_c,nl2_c,nu2_c,nl3_c,nu3_c, &
                      & nl1_f,nu1_f,nl2_f,nu2_f,nl3_f,nu3_f,radii_cf(iatype(iat),2), & 
                      & cpmult,fpmult,hgrid,gau_a,fac_arr,rx,ry,rz,lx,ly,lz, & 
                      & mvctr_c,mvctr_f,proj(istart_c), &
                      & proj(istart_f))

                 iproj=iproj+1
                 ! testing
                 call wnrm(mvctr_c,mvctr_f,proj(istart_c), &
                      & proj(istart_f),scpr)
                 if (abs(1.d0-scpr).gt.1.d-2) then
                    print *,'norm projector for atom ',trim(atomnames(iatype(iat))),&
                         'iproc,l,i,rl,scpr=',iproc,l,i,gau_a,scpr
                    stop 'norm projector'
                 end if

                 ! testing end
                 istart_c=istart_f+7*mvctr_f
                 if (istart_c.gt.istart) stop 'istart_c > istart'

                 !do iterm=1,nterm
                 !   if (iproc.eq.0) write(*,'(1x,a,i0,1x,a,1pe10.3,3(1x,i0))') &
                 !        'projector: iat,atomname,gau_a,lx,ly,lz ', & 
                 !        iat,trim(atomnames(iatype(iat))),gau_a,lx(iterm),ly(iterm),lz(iterm)
                 !enddo


              enddo
           endif
        enddo
     enddo
  enddo
  if (iproj.ne.nlpspd%nproj) stop 'incorrect number of projectors created'
  ! projector part finished
  if (iproc ==0) write(*,'(1x,a)')'done.'

  i_all=-product(shape(logrid))*kind(logrid)
  deallocate(logrid,stat=i_stat)
  call memocc(i_stat,i_all,'logrid','createprojectorsarrays')
  i_all=-product(shape(fac_arr))*kind(fac_arr)
  deallocate(fac_arr,stat=i_stat)
  call memocc(i_stat,i_all,'fac_arr','createprojectorsarrays')
  i_all=-product(shape(lx))*kind(lx)
  deallocate(lx,stat=i_stat)
  call memocc(i_stat,i_all,'lx','createprojectorsarrays')
  i_all=-product(shape(ly))*kind(ly)
  deallocate(ly,stat=i_stat)
  call memocc(i_stat,i_all,'ly','createprojectorsarrays')
  i_all=-product(shape(lz))*kind(lz)
  deallocate(lz,stat=i_stat)
  call memocc(i_stat,i_all,'lz','createprojectorsarrays')

END SUBROUTINE createProjectorsArrays

subroutine import_gaussians(parallel,iproc,nproc,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3, & 
     nat,norb,norbp,occup,n1,n2,n3,nvctrp,hgrid,rxyz, & 
     rhopot,pot_ion,wfd,bounds,nlpspd,proj,  &
     atomnames,ntypes,iatype,pkernel,psppar,npspcode,ixc,&
     psi,psit,hpsi,eval,accurex,datacode,nscatterarr,ngatherarr,nspin,spinar)

  use module_types
  use Poisson_Solver

  implicit none
  include 'mpif.h'
  type(wavefunctions_descriptors), intent(in) :: wfd
  type(convolutions_bounds), intent(in) :: bounds
  type(nonlocal_psp_descriptors), intent(in) :: nlpspd
  logical, intent(in) :: parallel
  character(len=20), dimension(100), intent(in) :: atomnames
  character(len=1), intent(in) :: datacode
  integer, intent(in) :: iproc,nproc,nat,ntypes,norb,norbp,n1,n2,n3,ixc
  integer, intent(in) :: nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,nvctrp
  integer, intent(in) :: nspin
  real(kind=8), dimension(norb), intent(in) :: spinar
  real(kind=8), intent(in) :: hgrid
  real(kind=8), intent(out) :: accurex
  integer, dimension(nat), intent(in) :: iatype
  integer, dimension(ntypes), intent(in) :: npspcode
  integer, dimension(0:nproc-1,4), intent(in) :: nscatterarr !n3d,n3p,i3s+i3xcsh-1,i3xcsh
  integer, dimension(0:nproc-1,2), intent(in) :: ngatherarr 
  real(kind=8), dimension(norb), intent(in) :: occup
  real(kind=8), dimension(3,nat), intent(in) :: rxyz
  real(kind=8), dimension(0:4,0:6,ntypes), intent(in) :: psppar
  real(kind=8), dimension(nlpspd%nprojel), intent(in) :: proj
  real(kind=8), dimension(*), intent(in) :: pkernel
  real(kind=8), dimension(*), intent(inout) :: rhopot,pot_ion
  real(kind=8), dimension(norb), intent(out) :: eval
  real(kind=8), dimension(:,:), pointer :: psi,psit,hpsi

  !local variables
  integer :: i,iorb,i_stat,i_all,ierr,info,jproc,n_lp,jorb
  real(kind=8) :: hgridh,tt,eks,eexcu,vexcu,epot_sum,ekin_sum,ehart,eproj_sum
  real(kind=8), dimension(:), allocatable :: work_lp,pot,ones
  real(kind=8), dimension(:,:), allocatable :: hamovr

  if (iproc.eq.0) then
     write(*,'(1x,a)')&
          '--------------------------------------------------------- Import Gaussians from CP2K'
  end if

  hgridh=.5d0*hgrid

  if (parallel) then
     !allocate the wavefunction in the transposed way to avoid allocations/deallocations
     allocate(psi(nvctrp,norbp*nproc),stat=i_stat)
     call memocc(i_stat,product(shape(psi))*kind(psi),'psi','import_gaussians')
  else
     allocate(psi(nvctrp,norb),stat=i_stat)
     call memocc(i_stat,product(shape(psi))*kind(psi),'psi','import_gaussians')
  end if

  !read the values for the gaussian code and insert them on psi 
  call gautowav(iproc,nproc,nat,ntypes,norb,norbp,n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,&
       wfd%nvctr_c,wfd%nvctr_f,wfd%nseg_c,wfd%nseg_f,wfd%keyg,wfd%keyv,&
       iatype,occup,rxyz,hgrid,psi,eks)

!!$  !!plot the initial gaussian wavefunctions
!!$  !do i=2*iproc+1,2*iproc+2
!!$  !   iounit=15+3*(i-1)
!!$  !   print *,'iounit',iounit,'-',iounit+2
!!$  !   call plot_wf(iounit,n1,n2,n3,hgrid,nseg_c,nvctr_c,keyg,keyv,nseg_f,nvctr_f,  & 
!!$  !        rxyz(1,1),rxyz(2,1),rxyz(3,1),psi(:,i-2*iproc:i-2*iproc), & 
!!$          ibyz_c,ibzxx_c,ibxxyy_c,ibyz_ff,ibzxx_f,ibxxyy_f,ibyyzz_r,&
!!$          nfl1,nfu1,nfl2,nfu2,nfl3,nfu3)
!!$  !end do


  ! resulting charge density and potential
  allocate(ones(norb),stat=i_stat)
  call memocc(i_stat,product(shape(ones))*kind(ones),'ones','import_gaussians')
  ones(:)=1.0d0

  call sumrho(parallel,iproc,nproc,norb,norbp,n1,n2,n3,hgrid,occup,  & 
       wfd,psi,rhopot,(2*n1+31)*(2*n2+31)*nscatterarr(iproc,1),nscatterarr,1,ones,&
       nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,bounds)

  call PSolver('F',datacode,iproc,nproc,2*n1+31,2*n2+31,2*n3+31,ixc,hgridh,hgridh,hgridh,&
       rhopot,pkernel,pot_ion,ehart,eexcu,vexcu,0.d0,.true.,1)

  if (parallel) then
     !allocate the wavefunction in the transposed way to avoid allocations/deallocations
     allocate(hpsi(nvctrp,norbp*nproc),stat=i_stat)
     call memocc(i_stat,product(shape(hpsi))*kind(hpsi),'hpsi','import_gaussians')
  else
     allocate(hpsi(nvctrp,norbp),stat=i_stat)
     call memocc(i_stat,product(shape(hpsi))*kind(hpsi),'hpsi','import_gaussians')
  end if


  call HamiltonianApplication(parallel,datacode,iproc,nproc,nat,ntypes,iatype,hgrid,&
       psppar,npspcode,norb,norbp,occup,n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,&
       wfd,bounds,nlpspd,proj,&
       ngatherarr,nscatterarr(iproc,2),rhopot(1+(2*n1+31)*(2*n2+31)*nscatterarr(iproc,4)),&
       psi,hpsi,ekin_sum,epot_sum,eproj_sum,1,ones)

  i_all=-product(shape(ones))*kind(ones)
  deallocate(ones,stat=i_stat)
  call memocc(i_stat,i_all,'ones','import_gaussians')


  accurex=abs(eks-ekin_sum)
  if (iproc.eq.0) write(*,'(1x,a,2(f19.10))') 'done. ekin_sum,eks:',ekin_sum,eks


  !after having applied the hamiltonian to all the atomic orbitals
  !we split the semicore orbitals from the valence ones
  !this is possible since the semicore orbitals are the first in the 
  !order, so the linear algebra on the transposed wavefunctions 
  !may be splitted

  if (iproc.eq.0) write(*,'(1x,a)',advance='no')&
       'Imported Wavefunctions Orthogonalization:'

  if (parallel) then

     !transpose all the wavefunctions for having a piece of all the orbitals 
     !for each processor
     !allocate the wavefunction in the transposed way to avoid allocations/deallocations
     allocate(psit(nvctrp,norbp*nproc),stat=i_stat)
     call memocc(i_stat,product(shape(psit))*kind(psit),'psit','import_gaussians')

     call timing(iproc,'Un-TransSwitch','ON')
     call switch_waves(iproc,nproc,norb,norbp,wfd%nvctr_c,wfd%nvctr_f,nvctrp,psi,psit)
     call timing(iproc,'Un-TransSwitch','OF')
     call timing(iproc,'Un-TransComm  ','ON')
     call MPI_ALLTOALL(psit,nvctrp*norbp,MPI_DOUBLE_PRECISION,  &
          psi,nvctrp*norbp,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
     call timing(iproc,'Un-TransComm  ','OF')

     call timing(iproc,'Un-TransSwitch','ON')
     call switch_waves(iproc,nproc,norb,norbp,wfd%nvctr_c,wfd%nvctr_f,nvctrp,hpsi,psit)
     call timing(iproc,'Un-TransSwitch','OF')
     call timing(iproc,'Un-TransComm  ','ON')
     call MPI_ALLTOALL(psit,nvctrp*norbp,MPI_DOUBLE_PRECISION,  &
          hpsi,nvctrp*norbp,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
     call timing(iproc,'Un-TransComm  ','OF')

     !end of transposition

     allocate(hamovr(norb**2,4),stat=i_stat)
     call memocc(i_stat,product(shape(hamovr))*kind(hamovr),'hamovr','import_gaussians')

     !calculate the overlap matrix for each group of the semicore atoms
     !       hamovr(jorb,iorb,3)=+psit(k,jorb)*hpsit(k,iorb)
     !       hamovr(jorb,iorb,4)=+psit(k,jorb)* psit(k,iorb)

     if (iproc.eq.0) write(*,'(1x,a)',advance='no')&
          'Overlap Matrix...'

     call DGEMM('T','N',norb,norb,nvctrp,1.d0,psi,nvctrp,hpsi,nvctrp,&
          0.d0,hamovr(1,3),norb)

     call DGEMM('T','N',norb,norb,nvctrp,1.d0,psi,nvctrp,psi,nvctrp,&
          0.d0,hamovr(1,4),norb)

     !reduce the overlap matrix between all the processors
     call MPI_ALLREDUCE(hamovr(1,3),hamovr(1,1),2*norb**2,&
          MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)

     !print the overlap matrix in the wavelet case
     print *,norb
     open(33)
     do iorb=1,norb
        write(33,'(2000(1pe10.2))')&
             (hamovr(jorb+(iorb-1)*norb,2),jorb=1,norb)
     end do
     close(33)

     !stop

     !found the eigenfunctions for each group
     n_lp=max(10,4*norb)
     allocate(work_lp(n_lp),stat=i_stat)
     call memocc(i_stat,product(shape(work_lp))*kind(work_lp),'work_lp','import_gaussians')

     if (iproc.eq.0) write(*,'(1x,a)')'Linear Algebra...'

     call DSYGV(1,'V','U',norb,hamovr(1,1),norb,hamovr(1,2),norb,eval,work_lp,n_lp,info)

     if (info.ne.0) write(*,*) 'DSYGV ERROR',info
!!$        !!write the matrices on a file
!!$        !open(33+2*(i-1))
!!$        !do jjorb=1,norbi
!!$        !   write(33+2*(i-1),'(2000(1pe10.2))')&
!!$        !        (hamovr(imatrst-1+jiorb+(jjorb-1)*norbi,1),jiorb=1,norbi)
!!$        !end do
!!$        !close(33+2*(i-1))
!!$        !open(34+2*(i-1))
!!$        !do jjorb=1,norbi
!!$        !   write(34+2*(i-1),'(2000(1pe10.2))')&
!!$        !        (hamovr(imatrst-1+jjorb+(jiorb-1)*norbi,1),jiorb=1,norbi)
!!$        !end do
!!$        !close(34+2*(i-1))

     if (iproc.eq.0) then
        do iorb=1,norb
           write(*,'(1x,a,i0,a,1x,1pe21.14)') 'eval(',iorb,')=',eval(iorb)
        enddo
     endif

     i_all=-product(shape(work_lp))*kind(work_lp)
     deallocate(work_lp,stat=i_stat)
     call memocc(i_stat,i_all,'work_lp','import_gaussians')

     if (iproc.eq.0) write(*,'(1x,a)',advance='no')'Building orthogonal Imported Wavefunctions...'

     !perform the vector-matrix multiplication for building the input wavefunctions
     ! ppsit(k,iorb)=+psit(k,jorb)*hamovr(jorb,iorb,1)

     call DGEMM('N','N',nvctrp,norb,norb,1.d0,psi,nvctrp,&
          hamovr(1,1),norb,0.d0,psit,nvctrp)

     i_all=-product(shape(hamovr))*kind(hamovr)
     deallocate(hamovr,stat=i_stat)
     call memocc(i_stat,i_all,'hamovr','import_gaussians')

     !retranspose the wavefunctions
     call timing(iproc,'Un-TransComm  ','ON')
     call MPI_ALLTOALL(psit,nvctrp*norbp,MPI_DOUBLE_PRECISION,  &
          hpsi,nvctrp*norbp,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
     call timing(iproc,'Un-TransComm  ','OF')
     call timing(iproc,'Un-TransSwitch','ON')
     call unswitch_waves(iproc,nproc,norb,norbp,wfd%nvctr_c,wfd%nvctr_f,nvctrp,hpsi,psi)
     call timing(iproc,'Un-TransSwitch','OF')

  else !serial case

     write(*,'(1x,a)',advance='no')'Overlap Matrix...'

     allocate(hamovr(norb**2,4),stat=i_stat)
     call memocc(i_stat,product(shape(hamovr))*kind(hamovr),'hamovr','import_gaussians')
     !hamovr(jorb,iorb,3)=+psi(k,jorb)*hpsi(k,iorb)
     call DGEMM('T','N',norb,norb,nvctrp,1.d0,psi,nvctrp,hpsi,nvctrp,&
          0.d0,hamovr(1,1),norb)
     call DGEMM('T','N',norb,norb,nvctrp,1.d0,psi,nvctrp,psi,nvctrp,&
          0.d0,hamovr(1,2),norb)

     n_lp=max(10,4*norb)
     allocate(work_lp(n_lp),stat=i_stat)
     call memocc(i_stat,product(shape(work_lp))*kind(work_lp),'work_lp','import_gaussians')

     write(*,'(1x,a)')'Linear Algebra...'
     call DSYGV(1,'V','U',norb,hamovr(1,1),norb,hamovr(1,2),norb,eval,work_lp,n_lp,info)

     if (info.ne.0) write(*,*) 'DSYGV ERROR',info
     if (iproc.eq.0) then
        do iorb=1,norb
           write(*,'(1x,a,i0,a,1x,1pe21.14)') 'evale(',iorb,')=',eval(iorb)
        enddo
     endif

     i_all=-product(shape(work_lp))*kind(work_lp)
     deallocate(work_lp,stat=i_stat)
     call memocc(i_stat,i_all,'work_lp','import_gaussians')

     write(*,'(1x,a)',advance='no')'Building orthogonal Imported Wavefunctions...'

     !copy the values into hpsi
     do iorb=1,norb
        do i=1,wfd%nvctr_c+7*wfd%nvctr_f
           hpsi(i,iorb)=psi(i,iorb)
        end do
     end do
     !ppsi(k,iorb)=+psi(k,jorb)*hamovr(jorb,iorb,1)
     call DGEMM('N','N',nvctrp,norb,norb,1.d0,hpsi,nvctrp,hamovr(1,1),norb,0.d0,psi,nvctrp)

     i_all=-product(shape(hamovr))*kind(hamovr)
     deallocate(hamovr,stat=i_stat)
     call memocc(i_stat,i_all,'hamovr','import_gaussians')

  endif

  if (iproc.eq.0) write(*,'(1x,a)')'done.'

END SUBROUTINE import_gaussians

subroutine input_wf_diag(parallel,iproc,nproc,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,&
     nat,natsc,norb,norbp,n1,n2,n3,nvctrp,hgrid,rxyz, & 
     rhopot,pot_ion,wfd,bounds,nlpspd,proj,  &
     atomnames,ntypes,iatype,iasctype,pkernel,nzatom,nelpsp,psppar,npspcode,ixc,&
     ppsi,ppsit,eval,accurex,datacode,nscatterarr,ngatherarr,nspin,spinar)
  ! Input wavefunctions are found by a diagonalization in a minimal basis set
  ! Each processors write its initial wavefunctions into the wavefunction file
  ! The files are then read by readwave

  use module_types
  use Poisson_Solver

  implicit none
  include 'mpif.h'
  type(wavefunctions_descriptors), intent(in) :: wfd
  type(nonlocal_psp_descriptors), intent(in) :: nlpspd
  type(convolutions_bounds), intent(in) :: bounds
  logical, intent(in) :: parallel
  character(len=20), dimension(100), intent(in) :: atomnames
  character(len=1), intent(in) :: datacode
  integer, intent(in) :: iproc,nproc,nat,natsc,ntypes,norb,norbp,n1,n2,n3,ixc
  integer, intent(in) :: nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,nvctrp
  integer, intent(in) :: nspin
  real(kind=8), dimension(norb), intent(in) :: spinar
  real(kind=8), intent(in) :: hgrid
  real(kind=8), intent(out) :: accurex
  integer, dimension(nat), intent(in) :: iatype
  integer, dimension(ntypes), intent(in) :: iasctype,npspcode,nzatom,nelpsp
  integer, dimension(0:nproc-1,4), intent(in) :: nscatterarr !n3d,n3p,i3s+i3xcsh-1,i3xcsh
  integer, dimension(0:nproc-1,2), intent(in) :: ngatherarr 
  real(kind=8), dimension(3,nat), intent(in) :: rxyz
  real(kind=8), dimension(0:4,0:6,ntypes), intent(in) :: psppar
  real(kind=8), dimension(nlpspd%nprojel), intent(in) :: proj
  real(kind=8), dimension(*), intent(in) :: pkernel
  real(kind=8), dimension(*), intent(inout) :: rhopot,pot_ion
  real(kind=8), dimension(norb), intent(out) :: eval
  real(kind=8), dimension(:,:), pointer :: ppsi,ppsit

  !local variables
  real(kind=8), parameter :: eps_mach=1.d-12
  integer, parameter :: ngx=31
  integer :: i,iorb,iorbsc,imatrsc,iorbst,imatrst,i_stat,i_all,ierr,info,jproc,jpst,norbeyou
  integer :: norbe,norbep,norbi,norbj,norbeme,ndim_hamovr,n_lp,norbi_max,norbsc
  integer :: ispin,norbu,norbd,iorbst2
  real(kind=8) :: hgridh,tt,eks,eexcu,vexcu,epot_sum,ekin_sum,ehart,eproj_sum
  logical, dimension(:,:), allocatable :: scorb
  integer, dimension(:), allocatable :: norbsc_arr,ng
  integer, dimension(:,:), allocatable :: nl
  real(kind=8), dimension(:), allocatable :: work_lp,pot,evale,occupe,ones
  real(kind=8), dimension(:,:), allocatable :: xp,occupat,hamovr,psi,hpsi
  real(kind=8), dimension(:,:,:), allocatable :: psiw,psiat

  !Calculate no. up and down orbitals for spin-polarized starting guess
  norbu=0
  norbd=0
  do iorb=1,norb
     if(spinar(iorb)>0.0d0) norbu=norbu+1
     if(spinar(iorb)<0.0d0) norbd=norbd+1
  end do

  allocate(xp(ngx,ntypes),stat=i_stat)
  call memocc(i_stat,product(shape(xp))*kind(xp),'xp','input_wf_diag')
  allocate(psiat(ngx,5,ntypes),stat=i_stat)
  call memocc(i_stat,product(shape(psiat))*kind(psiat),'psiat','input_wf_diag')
  allocate(occupat(5,ntypes),stat=i_stat)
  call memocc(i_stat,product(shape(occupat))*kind(occupat),'occupat','input_wf_diag')
  allocate(ng(ntypes),stat=i_stat)
  call memocc(i_stat,product(shape(ng))*kind(ng),'ng','input_wf_diag')
  allocate(nl(4,ntypes),stat=i_stat)
  call memocc(i_stat,product(shape(nl))*kind(nl),'nl','input_wf_diag')
  allocate(scorb(4,natsc),stat=i_stat)
  call memocc(i_stat,product(shape(scorb))*kind(scorb),'scorb','input_wf_diag')

  allocate(norbsc_arr(natsc+1),stat=i_stat)
  call memocc(i_stat,product(shape(norbsc_arr))*kind(norbsc_arr),'norbsc_arr','input_wf_diag')

  if (iproc.eq.0) then
     write(*,'(1x,a)')&
          '------------------------------------------------------- Input Wavefunctions Creation'
  end if

  ! Read the inguess.dat file or generate the input guess via the inguess_generator
  call readAtomicOrbitals(iproc,ngx,xp,psiat,occupat,ng,nl,nzatom,nelpsp,psppar,&
       & npspcode,norbe,norbsc,atomnames,ntypes,iatype,iasctype,nat,natsc,scorb,&
       & norbsc_arr)

  !  allocate wavefunctions and their occupation numbers
  allocate(occupe(norbe),stat=i_stat)
  call memocc(i_stat,product(shape(occupe))*kind(occupe),'occupe','input_wf_diag')
  tt=dble(norbe)/dble(nproc)
  norbep=int((1.d0-eps_mach*tt) + tt)

  if (iproc == 0 .and. nproc>1) then
     jpst=0
     do jproc=0,nproc-2
        norbeme=max(min((jproc+1)*norbep,norbe)-jproc*norbep,0)
        norbeyou=max(min((jproc+2)*norbep,norbe)-(jproc+1)*norbep,0)
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

  hgridh=.5d0*hgrid

  if (parallel) then
     !allocate the wavefunction in the transposed way to avoid allocations/deallocations
     allocate(psi(nvctrp,norbep*nproc),stat=i_stat)
     call memocc(i_stat,product(shape(psi))*kind(psi),'psi','input_wf_diag')
  else
     allocate(psi(nvctrp,norbe),stat=i_stat)
     call memocc(i_stat,product(shape(psi))*kind(psi),'psi','input_wf_diag')
  end if

  ! Create input guess orbitals
  call createAtomicOrbitals(iproc,nproc,atomnames,&
       nat,rxyz,norbe,norbep,norbsc,occupe,occupat,ngx,xp,psiat,ng,nl,&
       wfd%nvctr_c,wfd%nvctr_f,n1,n2,n3,hgrid,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,&
       wfd%nseg_c,wfd%nseg_f,wfd%keyg,wfd%keyv,iatype,ntypes,iasctype,natsc,psi,eks,scorb)

!!$  !!plot the initial LCAO wavefunctions
!!$  !do i=2*iproc+1,2*iproc+2
!!$  !   iounit=15+3*(i-1)
!!$  !   print *,'iounit',iounit,'-',iounit+2
!!$  !   call plot_wf(iounit,n1,n2,n3,hgrid,nseg_c,nvctr_c,keyg,keyv,nseg_f,nvctr_f,  & 
!!$  !        rxyz(1,1),rxyz(2,1),rxyz(3,1),psi(:,i-2*iproc:i-2*iproc), &
!!$          ibyz_c,ibzxx_c,ibxxyy_c,ibyz_ff,ibzxx_f,ibxxyy_f,ibyyzz_r,&
!!$          nfl1,nfu1,nfl2,nfu2,nfl3,nfu3)
!!$  !end do


  i_all=-product(shape(scorb))*kind(scorb)
  deallocate(scorb,stat=i_stat)
  call memocc(i_stat,i_all,'scorb','input_wf_diag')
  i_all=-product(shape(xp))*kind(xp)
  deallocate(xp,stat=i_stat)
  call memocc(i_stat,i_all,'xp','input_wf_diag')
  i_all=-product(shape(psiat))*kind(psiat)
  deallocate(psiat,stat=i_stat)
  call memocc(i_stat,i_all,'psiat','input_wf_diag')
  i_all=-product(shape(occupat))*kind(occupat)
  deallocate(occupat,stat=i_stat)
  call memocc(i_stat,i_all,'occupat','input_wf_diag')
  i_all=-product(shape(ng))*kind(ng)
  deallocate(ng,stat=i_stat)
  call memocc(i_stat,i_all,'ng','input_wf_diag')
  i_all=-product(shape(nl))*kind(nl)
  deallocate(nl,stat=i_stat)
  call memocc(i_stat,i_all,'nl','input_wf_diag')


  ! resulting charge density and potential
  allocate(ones(norbe),stat=i_stat)
  call memocc(i_stat,product(shape(ones))*kind(ones),'ones','input_wf_diag')
  ones(:)=1.0d0

  call sumrho(parallel,iproc,nproc,norbe,norbep,n1,n2,n3,hgrid,occupe,  & 
       wfd,psi,rhopot,(2*n1+31)*(2*n2+31)*nscatterarr(iproc,1),nscatterarr,1,ones, &
       nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,bounds)

  call PSolver('F',datacode,iproc,nproc,2*n1+31,2*n2+31,2*n3+31,ixc,hgridh,hgridh,hgridh,&
       rhopot,pkernel,pot_ion,ehart,eexcu,vexcu,0.d0,.true.,1)

  if (parallel) then
     !allocate the wavefunction in the transposed way to avoid allocations/deallocations
     allocate(hpsi(nvctrp,norbep*nproc),stat=i_stat)
     call memocc(i_stat,product(shape(hpsi))*kind(hpsi),'hpsi','input_wf_diag')
  else
     allocate(hpsi(nvctrp,norbe),stat=i_stat)
     call memocc(i_stat,product(shape(hpsi))*kind(hpsi),'hpsi','input_wf_diag')
  end if

  call HamiltonianApplication(parallel,datacode,iproc,nproc,nat,ntypes,iatype,hgrid,&
       psppar,npspcode,norbe,norbep,occupe,n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,&
       wfd,bounds,nlpspd,proj,&
       ngatherarr,nscatterarr(iproc,2),rhopot(1+(2*n1+31)*(2*n2+31)*nscatterarr(iproc,4)),&
       psi,hpsi,ekin_sum,epot_sum,eproj_sum,1,ones)

  i_all=-product(shape(ones))*kind(ones)
  deallocate(ones,stat=i_stat)
  call memocc(i_stat,i_all,'ones','input_wf_diag')

  i_all=-product(shape(occupe))*kind(occupe)
  deallocate(occupe,stat=i_stat)
  call memocc(i_stat,i_all,'occupe','input_wf_diag')

  accurex=abs(eks-ekin_sum)
  if (iproc.eq.0) write(*,'(1x,a,2(f19.10))') 'done. ekin_sum,eks:',ekin_sum,eks

  !after having applied the hamiltonian to all the atomic orbitals
  !we split the semicore orbitals from the valence ones
  !this is possible since the semicore orbitals are the first in the 
  !order, so the linear algebra on the transposed wavefunctions 
  !may be splitted

  if (iproc.eq.0) write(*,'(1x,a)',advance='no')&
       'Input Wavefunctions Orthogonalization:'

  norbi_max=maxval(norbsc_arr)

  !calculate the dimension of the overlap matrix
  ndim_hamovr=0
  do i=1,natsc+1
     ndim_hamovr=ndim_hamovr+norbsc_arr(i)**2
  end do

  if (parallel) then

     !transpose all the wavefunctions for having a piece of all the orbitals 
     !for each processor
     !here the timing is related to the input guess part
     allocate(psiw(nvctrp,norbep,nproc),stat=i_stat)
     call memocc(i_stat,product(shape(psiw))*kind(psiw),'psiw','input_wf_diag')

     call switch_waves(iproc,nproc,norbe,norbep,wfd%nvctr_c,wfd%nvctr_f,nvctrp,psi,psiw)
     call MPI_ALLTOALL(psiw,nvctrp*norbep,MPI_DOUBLE_PRECISION,  &
          psi,nvctrp*norbep,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)

     call switch_waves(iproc,nproc,norbe,norbep,wfd%nvctr_c,wfd%nvctr_f,nvctrp,hpsi,psiw)
     call MPI_ALLTOALL(psiw,nvctrp*norbep,MPI_DOUBLE_PRECISION,  &
          hpsi,nvctrp*norbep,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)

     i_all=-product(shape(psiw))*kind(psiw)
     deallocate(psiw,stat=i_stat)
     call memocc(i_stat,i_all,'psiw','input_wf_diag')
     !end of transposition

     allocate(hamovr(ndim_hamovr,4),stat=i_stat)
     call memocc(i_stat,product(shape(hamovr))*kind(hamovr),'hamovr','input_wf_diag')

     !calculate the overlap matrix for each group of the semicore atoms
     !       hamovr(jorb,iorb,3)=+psit(k,jorb)*hpsit(k,iorb)
     !       hamovr(jorb,iorb,4)=+psit(k,jorb)* psit(k,iorb)
     iorbst=1
     imatrst=1
     !print *,'norbi',norbi,natsc,norbsc_arr(natsc+1)

     if (iproc.eq.0) write(*,'(1x,a)',advance='no')&
          'Overlap Matrix...'

     do i=1,natsc
        norbi=norbsc_arr(i)
        call DGEMM('T','N',norbi,norbi,nvctrp,1.d0,psi(1,iorbst),nvctrp,hpsi(1,iorbst),nvctrp,&
             0.d0,hamovr(imatrst,3),norbi)
        call DGEMM('T','N',norbi,norbi,nvctrp,1.d0,psi(1,iorbst),nvctrp,psi(1,iorbst),nvctrp,&
             0.d0,hamovr(imatrst,4),norbi)
        iorbst=iorbst+norbi
        imatrst=imatrst+norbi**2
     end do
     norbi=norbsc_arr(natsc+1)
     call DGEMM('T','N',norbi,norbi,nvctrp,1.d0,psi(1,iorbst),nvctrp,hpsi(1,iorbst),nvctrp,&
          0.d0,hamovr(imatrst,3),norbi)

     i_all=-product(shape(hpsi))*kind(hpsi)
     deallocate(hpsi,stat=i_stat)
     call memocc(i_stat,i_all,'hpsi','input_wf_diag')

     call DGEMM('T','N',norbi,norbi,nvctrp,1.d0,psi(1,iorbst),nvctrp,psi(1,iorbst),nvctrp,&
          0.d0,hamovr(imatrst,4),norbi)

     !reduce the overlap matrix between all the processors
     call MPI_ALLREDUCE(hamovr(1,3),hamovr(1,1),2*ndim_hamovr,&
          MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)

     !found the eigenfunctions for each group
     n_lp=max(10,4*norbi_max)
     allocate(work_lp(n_lp),stat=i_stat)
     call memocc(i_stat,product(shape(work_lp))*kind(work_lp),'work_lp','input_wf_diag')
     allocate(evale(norbi_max),stat=i_stat)
     call memocc(i_stat,product(shape(evale))*kind(evale),'evale','input_wf_diag')

     if (iproc.eq.0) write(*,'(1x,a)')'Linear Algebra...'

     iorbst=1
     imatrst=1
     do i=1,natsc+1
        norbi=norbsc_arr(i)
        call DSYGV(1,'V','U',norbi,hamovr(imatrst,1),norbi,hamovr(imatrst,2),&
             norbi,evale,work_lp,n_lp,info)

        if (info.ne.0) write(*,*) 'DSYGV ERROR',info,i,natsc+1
!!$        !write the matrices on a file
!!$        !open(33+2*(i-1))
!!$        !do jjorb=1,norbi
!!$        !   write(33+2*(i-1),'(2000(1pe10.2))')&
!!$        !        (hamovr(imatrst-1+jiorb+(jjorb-1)*norbi,1),jiorb=1,norbi)
!!$        !end do
!!$        !close(33+2*(i-1))
!!$        !open(34+2*(i-1))
!!$        !do jjorb=1,norbi
!!$        !   write(34+2*(i-1),'(2000(1pe10.2))')&
!!$        !        (hamovr(imatrst-1+jjorb+(jiorb-1)*norbi,1),jiorb=1,norbi)
!!$        !end do
!!$        !close(34+2*(i-1))

        if (iproc.eq.0) then
           do iorb=1,norbi
              if (nspin==1) then
                 if (iorb+iorbst-1 == norb) then
                    write(*,'(1x,a,i0,a,1x,1pe21.14,a)') &
                         'evale(',iorb+iorbst-1,')=',evale(iorb),&
                         ' <- Last eigenvalue for input wavefncts'
                 else
                    write(*,'(1x,a,i0,a,1x,1pe21.14)') &
                         'evale(',iorb+iorbst-1,')=',evale(iorb)
                 end if
              else
                 if (norbu==norbd) then
                    if (iorb+iorbst-1 == norbu) then
                       write(*,'(1x,a,i0,a,1x,1pe21.14,a)') &
                            'evale(',iorb+iorbst-1,')=',evale(iorb),&
                            ' <- Last eigenvalue for input wavefncts of both spins'
                    else
                       write(*,'(1x,a,i0,a,1x,1pe21.14)') &
                            'evale(',iorb+iorbst-1,')=',evale(iorb)
                    end if
                 else
                    if (iorb+iorbst-1 == norbu) then
                       write(*,'(1x,a,i0,a,1x,1pe21.14,a)') &
                            'evale(',iorb+iorbst-1,')=',evale(iorb),&
                            ' <- Last eigenvalue for input wavefncts of spin up'
                    else if (iorb+iorbst-1 == norbd) then
                       write(*,'(1x,a,i0,a,1x,1pe21.14,a)') &
                            'evale(',iorb+iorbst-1,')=',evale(iorb),&
                            ' <- Last eigenvalue for input wavefncts of spin down'
                    else
                       write(*,'(1x,a,i0,a,1x,1pe21.14)') &
                            'evale(',iorb+iorbst-1,')=',evale(iorb)
                    end if
                 end if
              end if
           enddo
        endif
        do iorb=iorbst,min(norbi+iorbst-1,norb)
           eval(iorb)=evale(iorb-iorbst+1)
        enddo
        iorbst=iorbst+norbi
        imatrst=imatrst+norbi**2
     end do

     !     write(*,*) "NORBO",norb,norbe,norbu,norbd,norbu+norbd,norbp
     ! Copy eigenvalues from NM to spin-polarized channels
     if(nspin>1) then
        !        do iorb=1,norbu
        !!           write(*,*) 'jorb:',iorb,norbi_max,norbe,norb,product(shape(eval)),product(shape(evale))
        !           evale(iorb)=eval(iorb)
        !        end do
        do iorb=1,norbd
           !           write(*,*) 'korb:',iorb,iorb+norbu,norbi_max,norbe,norb,product(shape(evale))
           eval(iorb+norbu)=eval(iorb)
        end do
     end if

     i_all=-product(shape(work_lp))*kind(work_lp)
     deallocate(work_lp,stat=i_stat)
     call memocc(i_stat,i_all,'work_lp','input_wf_diag')
     i_all=-product(shape(evale))*kind(evale)
     deallocate(evale,stat=i_stat)
     call memocc(i_stat,i_all,'evale','input_wf_diag')

     !allocate the transposed wavefunction
     allocate(ppsit(nvctrp,norbp*nproc),stat=i_stat)
     call memocc(i_stat,product(shape(ppsit))*kind(ppsit),'ppsit','input_wf_diag')

     if (iproc.eq.0) write(*,'(1x,a)',advance='no')'Building orthogonal Input Wavefunctions...'

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
     do ispin=1,nspin
        iorbst=1
        iorbst2=0
        if(ispin==2) iorbst2=norbu
        imatrst=1
        norbsc=0
        do i=1,natsc
           norbi=norbsc_arr(i)
           norbsc=norbsc+norbi
           call DGEMM('N','N',nvctrp,norbi,norbi,1.d0,psi(1,iorbst),nvctrp,&
                hamovr(imatrst,1),norbi,0.d0,ppsit(1,iorbst2+iorbst),nvctrp)
           iorbst=iorbst+norbi
           imatrst=imatrst+norbi**2
        end do
        norbi=norbsc_arr(natsc+1)
        if(ispin==1) norbj=norbu-norbsc
        if(ispin==2) norbj=norbd-norbsc
        !        write(*,'(1x,a,5i4)') "DIMS:",norbi,norbj,iorbst,imatrst
        !        norbj=norb-norbsc
        if(norbj>0) then
           call DGEMM('N','N',nvctrp,norbj,norbi,1.d0,psi(1,iorbst),nvctrp,&
                hamovr(imatrst,1),norbi,0.d0,ppsit(1,iorbst2+iorbst),nvctrp)
        end if
     end do


     i_all=-product(shape(psi))*kind(psi)
     deallocate(psi,stat=i_stat)
     call memocc(i_stat,i_all,'psi','input_wf_diag')
     i_all=-product(shape(hamovr))*kind(hamovr)
     deallocate(hamovr,stat=i_stat)
     call memocc(i_stat,i_all,'hamovr','input_wf_diag')
     i_all=-product(shape(norbsc_arr))*kind(norbsc_arr)
     deallocate(norbsc_arr,stat=i_stat)
     call memocc(i_stat,i_all,'norbsc_arr','input_wf_diag')

     !orthogonalise the orbitals in the case of semi-core atoms
     if (norbsc > 0) then
        if(nspin==1) then
           call orthon_p(iproc,nproc,norb,nvctrp,wfd%nvctr_c+7*wfd%nvctr_f,ppsit) 
        else
           call orthon_p(iproc,nproc,norbu,nvctrp,wfd%nvctr_c+7*wfd%nvctr_f,ppsit) 
           if(norbd>0) then
              call orthon_p(iproc,nproc,norbd,nvctrp,wfd%nvctr_c+7*wfd%nvctr_f,ppsit(1,norbu+1)) 
           end if
        end if
     end if

     !retranspose the wavefunctions
     allocate(psiw(nvctrp,norbp,nproc),stat=i_stat)
     call memocc(i_stat,product(shape(psiw))*kind(psiw),'psiw','input_wf_diag')

     call timing(iproc,'Un-TransComm  ','ON')
     call MPI_ALLTOALL(ppsit,nvctrp*norbp,MPI_DOUBLE_PRECISION,  &
          psiw,nvctrp*norbp,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
     call timing(iproc,'Un-TransComm  ','OF')

     !allocate the direct wavefunction
     allocate(ppsi(nvctrp,norbp*nproc),stat=i_stat)
     call memocc(i_stat,product(shape(ppsi))*kind(ppsi),'ppsi','input_wf_diag')

     call timing(iproc,'Un-TransSwitch','ON')
     call unswitch_waves(iproc,nproc,norb,norbp,wfd%nvctr_c,wfd%nvctr_f,nvctrp,psiw,ppsi)
     call timing(iproc,'Un-TransSwitch','OF')

     i_all=-product(shape(psiw))*kind(psiw)
     deallocate(psiw,stat=i_stat)
     call memocc(i_stat,i_all,'psiw','input_wf_diag')

  else !serial case

     write(*,'(1x,a)',advance='no')'Overlap Matrix...'

     allocate(hamovr(ndim_hamovr,2),stat=i_stat)
     call memocc(i_stat,product(shape(hamovr))*kind(hamovr),'hamovr','input_wf_diag')
     !hamovr(jorb,iorb,3)=+psi(k,jorb)*hpsi(k,iorb)
     iorbst=1
     imatrst=1
     do i=1,natsc+1
        norbi=norbsc_arr(i)
        call DGEMM('T','N',norbi,norbi,nvctrp,1.d0,psi(1,iorbst),nvctrp,hpsi(1,iorbst),nvctrp,&
             0.d0,hamovr(imatrst,1),norbi)
        call DGEMM('T','N',norbi,norbi,nvctrp,1.d0,psi(1,iorbst),nvctrp,psi(1,iorbst),nvctrp,&
             0.d0,hamovr(imatrst,2),norbi)
        iorbst=iorbst+norbi
        imatrst=imatrst+norbi**2
     end do

     i_all=-product(shape(hpsi))*kind(hpsi)
     deallocate(hpsi,stat=i_stat)
     call memocc(i_stat,i_all,'hpsi','input_wf_diag')


     n_lp=max(10,4*norbi_max)
     allocate(work_lp(n_lp),stat=i_stat)
     call memocc(i_stat,product(shape(work_lp))*kind(work_lp),'work_lp','input_wf_diag')
     allocate(evale(max(norbi_max,norb)),stat=i_stat)
     call memocc(i_stat,product(shape(evale))*kind(evale),'evale','input_wf_diag')

     write(*,'(1x,a)')'Linear Algebra...'
     iorbst=1
     imatrst=1
     do i=1,natsc+1
        norbi=norbsc_arr(i)
        call DSYGV(1,'V','U',norbi,hamovr(imatrst,1),norbi,hamovr(imatrst,2),&
             norbi,evale,work_lp,n_lp,info)

        if (info.ne.0) write(*,*) 'DSYGV ERROR',info,i,natsc+1
        if (iproc.eq.0) then
           do iorb=1,norbi
              if (nspin==1) then
                 if (iorb+iorbst-1 == norb) then
                    write(*,'(1x,a,i0,a,1x,1pe21.14,a)') &
                         'evale(',iorb+iorbst-1,')=',evale(iorb),&
                         ' <- Last eigenvalue for input wavefncts'
                 else
                    write(*,'(1x,a,i0,a,1x,1pe21.14)') &
                         'evale(',iorb+iorbst-1,')=',evale(iorb)
                 end if
              else
                 if (norbu==norbd) then
                    if (iorb+iorbst-1 == norbu) then
                       write(*,'(1x,a,i0,a,1x,1pe21.14,a)') &
                            'evale(',iorb+iorbst-1,')=',evale(iorb),&
                            ' <- Last eigenvalue for input wavefncts of both spins'
                    else
                       write(*,'(1x,a,i0,a,1x,1pe21.14)') &
                            'evale(',iorb+iorbst-1,')=',evale(iorb)
                    end if
                 else
                    if (iorb+iorbst-1 == norbu) then
                       write(*,'(1x,a,i0,a,1x,1pe21.14,a)') &
                            'evale(',iorb+iorbst-1,')=',evale(iorb),&
                            ' <- Last eigenvalue for input wavefncts of spin up'
                    else if (iorb+iorbst-1 == norbd) then
                       write(*,'(1x,a,i0,a,1x,1pe21.14,a)') &
                            'evale(',iorb+iorbst-1,')=',evale(iorb),&
                            ' <- Last eigenvalue for input wavefncts of spin down'
                    else
                       write(*,'(1x,a,i0,a,1x,1pe21.14)') &
                            'evale(',iorb+iorbst-1,')=',evale(iorb)
                    end if
                 end if
              end if
           enddo
        endif
        do iorb=iorbst,min(norbi+iorbst-1,norb)
           !           write(*,*) 'iorb:',iorb,iorb-iorbst+1,norbi_max,norbe,norb,product(shape(evale))
           eval(iorb)=evale(iorb-iorbst+1)
        enddo
        iorbst=iorbst+norbi
        imatrst=imatrst+norbi**2
     end do
     !     write(*,*) "NORBO",norb,norbe,norbu,norbd,norbu+norbd,norbp
     ! Copy eigenvalues from NM to spin-polarized channels
     if(nspin>1) then
        !        do iorb=1,norbu
        !!           write(*,*) 'jorb:',iorb,norbi_max,norbe,norb,product(shape(eval)),product(shape(evale))
        !           evale(iorb)=eval(iorb)
        !        end do
        do iorb=1,norbd
           !           write(*,*) 'korb:',iorb,iorb+norbu,norbi_max,norbe,norb,product(shape(evale))
           eval(iorb+norbu)=eval(iorb)
        end do
     end if
     !
     i_all=-product(shape(work_lp))*kind(work_lp)
     deallocate(work_lp,stat=i_stat)
     call memocc(i_stat,i_all,'work_lp','input_wf_diag')
     i_all=-product(shape(evale))*kind(evale)
     deallocate(evale,stat=i_stat)
     call memocc(i_stat,i_all,'evale','input_wf_diag')

     write(*,'(1x,a)',advance='no')'Building orthogonal Input Wavefunctions...'

     !allocate the wavefunction
     allocate(ppsi(nvctrp,norb),stat=i_stat)
     call memocc(i_stat,product(shape(ppsi))*kind(ppsi),'ppsi','input_wf_diag')

     !ppsi(k,iorb)=+psi(k,jorb)*hamovr(jorb,iorb,1)
     do ispin=1,nspin
        iorbst=1
        iorbst2=0
        if(ispin==2) iorbst2=norbu
        imatrst=1
        norbsc=0
        do i=1,natsc
           norbi=norbsc_arr(i)
           norbsc=norbsc+norbi
           call DGEMM('N','N',nvctrp,norbi,norbi,1.d0,psi(1,iorbst),nvctrp,&
                hamovr(imatrst,1),norbi,0.d0,ppsi(1,iorbst2+iorbst),nvctrp)
           iorbst=iorbst+norbi
           imatrst=imatrst+norbi**2
        end do
        norbi=norbsc_arr(natsc+1)
        if(ispin==1) norbj=norbu-norbsc
        if(ispin==2) norbj=norbd-norbsc
        !        write(*,'(1x,a,5i4)') "DIMS:",norbi,norbj,iorbst,imatrst
        !        norbj=norb-norbsc
        if(norbj>0) then
           call DGEMM('N','N',nvctrp,norbj,norbi,1.d0,psi(1,iorbst),nvctrp,&
                hamovr(imatrst,1),norbi,0.d0,ppsi(1,iorbst2+iorbst),nvctrp)
        end if
     end do

     i_all=-product(shape(psi))*kind(psi)
     deallocate(psi,stat=i_stat)
     call memocc(i_stat,i_all,'psi','input_wf_diag')
     i_all=-product(shape(hamovr))*kind(hamovr)
     deallocate(hamovr,stat=i_stat)
     call memocc(i_stat,i_all,'hamovr','input_wf_diag')
     i_all=-product(shape(norbsc_arr))*kind(norbsc_arr)
     deallocate(norbsc_arr,stat=i_stat)
     call memocc(i_stat,i_all,'norbsc_arr','input_wf_diag')


     !orthogonalise the orbitals in the case of semi-core atoms
     if (norbsc > 0) then
        if(nspin==1) then
           call orthon(norb,nvctrp,ppsi)
        else
           call orthon(norbu,nvctrp,ppsi)
           if(norbd>0) then
              call orthon(norbd,nvctrp,ppsi(1,norbu+1))
           end if
        end if
     end if

  endif

  if (iproc.eq.0) write(*,'(1x,a)')'done.'

END SUBROUTINE input_wf_diag
