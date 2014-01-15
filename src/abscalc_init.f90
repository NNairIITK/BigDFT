
!> Fill the preconditioning projectors for a given atom 
subroutine fillPcProjOnTheFly(PPD, Glr, iat, at, hx,hy,hz,startjorb,ecut_pc,   initial_istart_c ) 
  use module_interfaces
  use module_base
  use module_types
  use module_abscalc
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


  call plr_segs_and_vctrs(PPD%pc_nl%pspd(iat)%plr,&
       Plr%wfd%nseg_c,Plr%wfd%nseg_f,Plr%wfd%nvctr_c,Plr%wfd%nvctr_f)
!!$  Plr%wfd%nvctr_c  =PPD%pc_nlpspd%nvctr_p(2*iat-1)-PPD%pc_nlpspd%nvctr_p(2*iat-2)
!!$  Plr%wfd%nvctr_f  =PPD%pc_nlpspd%nvctr_p(2*iat  )-PPD%pc_nlpspd%nvctr_p(2*iat-1)
!!$  Plr%wfd%nseg_c   =PPD%pc_nlpspd%nseg_p(2*iat-1)-PPD%pc_nlpspd%nseg_p(2*iat-2)
!!$  Plr%wfd%nseg_f   =PPD%pc_nlpspd%nseg_p(2*iat  )-PPD%pc_nlpspd%nseg_p(2*iat-1)

  call allocate_wfd(Plr%wfd,subname)

  call vcopy(Plr%wfd%nseg_c+Plr%wfd%nseg_f,&
       PPD%pc_nl%pspd(iat)%plr%wfd%keyvglob(1),1,Plr%wfd%keyvglob(1),1)
  call vcopy(2*(Plr%wfd%nseg_c+Plr%wfd%nseg_f),&
       PPD%pc_nl%pspd(iat)%plr%wfd%keyglob(1,1),1,Plr%wfd%keyglob(1,1),1)

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
  use module_abscalc
  implicit none
  type(pawproj_data_type),  intent(in) ::PAWD
  type(locreg_descriptors),  intent(in):: Glr
  integer, intent(in)  ::iat, startjorb
  real(gp), intent(in) ::   hx,hy,hz,kx,ky,kz
  integer, intent(in) :: initial_istart_c
  character(len=1), intent(in) :: geocode !< @copydoc poisson_solver::doc::geocode
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

  call plr_segs_and_vctrs(PAWD%paw_nl%pspd(iat)%plr,&
       Plr%wfd%nseg_c,Plr%wfd%nseg_f,Plr%wfd%nvctr_c,Plr%wfd%nvctr_f)

  call allocate_wfd(Plr%wfd,subname)

  call vcopy(Plr%wfd%nseg_c+Plr%wfd%nseg_f,&
       PAWD%paw_nl%pspd(iat)%plr%wfd%keyvglob(1),1,Plr%wfd%keyvglob(1),1)
  call vcopy(2*(Plr%wfd%nseg_c+Plr%wfd%nseg_f),&
       PAWD%paw_nl%pspd(iat)%plr%wfd%keyglob(1,1),1,Plr%wfd%keyglob(1,1),1)

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


!rm?!>   Determine localization region for all preconditioning projectors, but do not yet fill the descriptor arrays
!rm?subroutine createPcProjectorsArrays(iproc,n1,n2,n3,rxyz,at,orbs,&
!rm?     &   radii_cf,cpmult,fpmult,hx,hy,hz, ecut_pc, &
!rm?     &   PPD, Glr)
!rm?  use module_interfaces, except_this_one => createPcProjectorsArrays
!rm?  use module_base
!rm?  use module_types
!rm?  implicit none
!rm?  integer, intent(in) :: iproc,n1,n2,n3
!rm?  real(gp), intent(in) :: cpmult,fpmult,hx,hy,hz
!rm?  type(atoms_data), intent(in) :: at
!rm?  type(orbitals_data), intent(in) :: orbs
!rm?
!rm?  real(gp), dimension(3,at%astruct%nat), intent(in) :: rxyz
!rm?  real(gp), dimension(at%astruct%ntypes,3), intent(in) :: radii_cf
!rm?  real(gp), intent(in):: ecut_pc
!rm?
!rm?  type(pcproj_data_type) ::PPD
!rm?
!rm?  type(locreg_descriptors),  intent(in):: Glr
!rm?
!rm?
!rm?  !local variables
!rm?  character(len=*), parameter :: subname='createPcProjectorsArrays'
!rm?  integer :: nl1,nl2,nl3,nu1,nu2,nu3,mseg,mproj, mvctr
!rm?  integer :: iat,i_stat,i_all,iseg, istart_c
!rm?  logical, dimension(:,:,:), allocatable :: logrid
!rm?
!rm?
!rm?  integer :: ng
!rm?  logical :: enlargerprb
!rm?  real(wp), dimension(:), pointer :: Gocc
!rm?
!rm?  integer, pointer :: iorbto_l(:)
!rm?  integer, pointer :: iorbto_m(:)
!rm?  integer, pointer :: iorbto_ishell(:)
!rm?  integer, pointer :: iorbto_iexpobeg(:)
!rm?
!rm?  integer :: nspin
!rm?  integer ::  jorb
!rm?  integer :: iproj, startjorb
!rm?  real(gp) :: Pcpmult
!rm?  integer :: mprojtot, nvctr_c, nvctr_f
!rm?  integer :: nprojel_tmp
!rm?
!rm?  Pcpmult=1.5*cpmult
!rm?
!rm?  ng=21
!rm?  enlargerprb = .false.
!rm?  nspin=1
!rm?
!rm?
!rm?  nullify(PPD%G%rxyz)
!rm?  call gaussian_pswf_basis(ng,enlargerprb,iproc,nspin,at,rxyz,PPD%G,Gocc, PPD%gaenes, &
!rm?       &   PPD%iorbtolr,iorbto_l, iorbto_m,  iorbto_ishell,iorbto_iexpobeg  )  
!rm?
!rm?
!rm?  ! allocated  : gaenes, Gocc , PPD%iorbtolr,iorbto_l, iorbto_m,  iorbto_ishell,iorbto_iexpobeg
!rm?
!rm?
!rm?!!$ ========================================================================================
!rm?
!rm?
!rm?
!rm?  PPD%pc_nlpspd%natoms=at%astruct%nat
!rm?  allocate(PPD%pc_nlpspd%plr(at%astruct%nat),stat=i_stat)
!rm?
!rm?!!$   allocate(PPD%pc_nlpspd%nseg_p(0:2*at%astruct%nat+ndebug),stat=i_stat)
!rm?!!$   call memocc(i_stat,PPD%pc_nlpspd%nseg_p,'pc_nlpspd%nseg_p',subname)
!rm?!!$   allocate(PPD%pc_nlpspd%nvctr_p(0:2*at%astruct%nat+ndebug),stat=i_stat)
!rm?!!$   call memocc(i_stat,PPD%pc_nlpspd%nvctr_p,'pc_nlpspd%nvctr_p',subname)
!rm?!!$   allocate(PPD%pc_nlpspd%nboxp_c(2,3,at%astruct%nat+ndebug),stat=i_stat)
!rm?!!$   call memocc(i_stat,PPD%pc_nlpspd%nboxp_c,'pc_nlpspd%nboxp_c',subname)
!rm?!!$   allocate(PPD%pc_nlpspd%nboxp_f(2,3,at%astruct%nat+ndebug),stat=i_stat)
!rm?!!$   call memocc(i_stat,PPD%pc_nlpspd%nboxp_f,'pc_nlpspd%nboxp_f',subname)
!rm?
!rm?  allocate(logrid(0:n1,0:n2,0:n3+ndebug),stat=i_stat)
!rm?  call memocc(i_stat,logrid,'logrid',subname)
!rm?
!rm?
!rm?  call localize_projectors(iproc,n1,n2,n3,hx,hy,hz,Pcpmult,fpmult,rxyz,radii_cf,&
!rm?       &   logrid,at,orbs,PPD%pc_nlpspd)
!rm?
!rm?  ! the above routine counts atomic projector and the number of their element for psp
!rm?  ! We must therefore correct , later, nlpspd%nprojel  and nlpspd%nproj
!rm?  !-------------------
!rm?
!rm?  ! allocations for arrays holding the projectors and their data descriptors
!rm?  !here the allocation is possible
!rm?  do iat=1,PPD%pc_nlpspd%natoms
!rm?     !for the moments the bounds are not needed for projectors
!rm?     call allocate_wfd(PPD%pc_nl%pspd(iat)%plr%wfd,subname)
!rm?  end do
!rm?
!rm?!!$   allocate(PPD%pc_nlpspd%keyg_p(2,PPD%pc_nlpspd%nseg_p(2*at%astruct%nat)+ndebug),stat=i_stat)
!rm?!!$   call memocc(i_stat,PPD%pc_nlpspd%keyg_p,'pc_nlpspd%keyg_p',subname)
!rm?!!$
!rm?!!$
!rm?!!$   allocate(PPD%pc_nlpspd%keyv_p(PPD%pc_nlpspd%nseg_p(2*at%astruct%nat)+ndebug),stat=i_stat)
!rm?!!$   call memocc(i_stat,PPD%pc_nlpspd%keyv_p,'pc_nlpspd%keyv_p',subname)
!rm?
!rm?
!rm?
!rm?!!$  -- this one delayed, waiting for the correct pc_nlpspd%nprojel, pc_nlpspd%nproj
!rm?!!$  --
!rm?!!$  allocate(pc_proj(pc_nlpspd%nprojel+ndebug),stat=i_stat)
!rm?!!$  call memocc(i_stat,pc_proj,'pc_proj',subname)
!rm?  PPD%ecut_pc=ecut_pc
!rm?
!rm?  PPD%pc_nlpspd%nprojel=0
!rm?  PPD%pc_nlpspd%nproj  =0
!rm?
!rm?!!$ =========================================================================================  
!rm?
!rm?  mprojtot=0
!rm?  jorb=1  
!rm?  ! After having determined the size of the projector descriptor arrays fill them
!rm?  do iat=1,at%astruct%nat
!rm?
!rm?     mproj=0
!rm?
!rm?     do while( jorb<=PPD%G%ncoeff         .and. PPD%iorbtolr(jorb)== iat)
!rm?
!rm?        if( PPD%gaenes(jorb)<ecut_pc) then
!rm?           mproj=mproj+1
!rm?        endif
!rm?        if(jorb==PPD%G%ncoeff) exit
!rm?        jorb=jorb+1
!rm?     end do
!rm?
!rm?     mprojtot=mprojtot+mproj
!rm?
!rm?     PPD%pc_nlpspd%nproj=PPD%pc_nlpspd%nproj+mproj
!rm?
!rm?
!rm?     if (mproj.ne.0) then 
!rm?
!rm?        nprojel_tmp=0
!rm?
!rm?        call bounds_to_plr_limits(.false.,1,PPD%pc_nl%pspd(iat)%plr,nl1,nl2,nl3,nu1,nu2,nu3)
!rm?!!$         ! coarse grid quantities
!rm?!!$         nl1=PPD%pc_nlpspd%nboxp_c(1,1,iat) 
!rm?!!$         nl2=PPD%pc_nlpspd%nboxp_c(1,2,iat) 
!rm?!!$         nl3=PPD%pc_nlpspd%nboxp_c(1,3,iat) 
!rm?!!$
!rm?!!$         nu1=PPD%pc_nlpspd%nboxp_c(2,1,iat)
!rm?!!$         nu2=PPD%pc_nlpspd%nboxp_c(2,2,iat)
!rm?!!$         nu3=PPD%pc_nlpspd%nboxp_c(2,3,iat)
!rm?
!rm?        call fill_logrid(at%astruct%geocode,n1,n2,n3,nl1,nu1,nl2,nu2,nl3,nu3,0,1,  &
!rm?             &   at%astruct%ntypes,at%astruct%iatype(iat),rxyz(1,iat),radii_cf(1,3),Pcpmult,hx,hy,hz,logrid)
!rm?
!rm?!!$         iseg=PPD%pc_nlpspd%nseg_p(2*iat-2)+1
!rm?!!$         mseg=PPD%pc_nlpspd%nseg_p(2*iat-1)-PPD%pc_nlpspd%nseg_p(2*iat-2)
!rm?        mseg=PPD%pc_nl%pspd(iat)%plr%wfd%nseg_c
!rm?
!rm?        call segkeys(n1,n2,n3,nl1,nu1,nl2,nu2,nl3,nu3,  & 
!rm?                                !!$         logrid,mseg,PPD%pc_nlpspd%keyg_p(1,iseg),PPD%pc_nlpspd%keyv_p(iseg))
!rm?             logrid,mseg,PPD%pc_nl%pspd(iat)%plr%wfd%keyglob(1,1),PPD%pc_nl%pspd(iat)%plr%wfd%keyvglob(1))
!rm?
!rm?!!$         mvctr =PPD%pc_nlpspd%nvctr_p(2*iat-1)-PPD%pc_nlpspd%nvctr_p(2*iat-2)
!rm?        mvctr =PPD%pc_nl%pspd(iat)%plr%wfd%nvctr_c
!rm?
!rm?        nprojel_tmp =nprojel_tmp +mproj*mvctr
!rm?
!rm?        call bounds_to_plr_limits(.false.,2,PPD%pc_nl%pspd(iat)%plr,nl1,nl2,nl3,nu1,nu2,nu3)
!rm?!!$         ! fine grid quantities
!rm?!!$         nl1=PPD%pc_nlpspd%nboxp_f(1,1,iat)
!rm?!!$         nl2=PPD%pc_nlpspd%nboxp_f(1,2,iat)
!rm?!!$         nl3=PPD%pc_nlpspd%nboxp_f(1,3,iat)
!rm?!!$
!rm?!!$         nu1=PPD%pc_nlpspd%nboxp_f(2,1,iat)
!rm?!!$         nu2=PPD%pc_nlpspd%nboxp_f(2,2,iat)
!rm?!!$         nu3=PPD%pc_nlpspd%nboxp_f(2,3,iat)
!rm?        call fill_logrid(at%astruct%geocode,n1,n2,n3,nl1,nu1,nl2,nu2,nl3,nu3,0,1,  &
!rm?             &   at%astruct%ntypes,at%astruct%iatype(iat),rxyz(1,iat),radii_cf(1,2),fpmult,hx,hy,hz,logrid)
!rm?!!$         iseg=PPD%pc_nlpspd%nseg_p(2*iat-1)+1
!rm?!!$         mseg=PPD%pc_nlpspd%nseg_p(2*iat)-PPD%pc_nlpspd%nseg_p(2*iat-1)
!rm?        iseg=PPD%pc_nl%pspd(iat)%plr%wfd%nseg_c+1
!rm?        mseg=PPD%pc_nl%pspd(iat)%plr%wfd%nseg_f
!rm?
!rm?        if (mseg > 0) then
!rm?           call segkeys(n1,n2,n3,nl1,nu1,nl2,nu2,nl3,nu3,  & 
!rm?                logrid,mseg,&
!rm?                                !!$                 PPD%pc_nlpspd%keyg_p(1,iseg),PPD%pc_nlpspd%keyv_p(iseg))
!rm?                PPD%pc_nl%pspd(iat)%plr%wfd%keyglob(1,iseg),&
!rm?                PPD%pc_nl%pspd(iat)%plr%wfd%keyvglob(iseg))
!rm?
!rm?           mvctr =PPD%pc_nl%pspd(iat)%plr%wfd%nvctr_f!PPD%pc_nlpspd%nvctr_p(2*iat)-PPD%pc_nlpspd%nvctr_p(2*iat-1)
!rm?
!rm?           nprojel_tmp=nprojel_tmp+mproj*mvctr*7
!rm?
!rm?        end if
!rm?
!rm?        if( PPD%DistProjApply)  then
!rm?           PPD%pc_nlpspd%nprojel=max(PPD%pc_nlpspd%nprojel,nprojel_tmp   )
!rm?        else
!rm?           PPD%pc_nlpspd%nprojel= PPD%pc_nlpspd%nprojel+nprojel_tmp 
!rm?        endif
!rm?
!rm?
!rm?     endif
!rm?
!rm?  enddo
!rm?
!rm?
!rm?  allocate(PPD%pc_proj(PPD%pc_nlpspd%nprojel+ndebug),stat=i_stat)
!rm?  call memocc(i_stat,PPD%pc_proj,'pc_proj',subname)
!rm?
!rm?  allocate(PPD%ilr_to_mproj(at%astruct%nat  +ndebug ) , stat=i_stat)
!rm?  call memocc(i_stat,PPD%ilr_to_mproj,'ilr_to_mproj',subname)
!rm?
!rm?  allocate(PPD%iproj_to_ene(mprojtot +ndebug ) , stat=i_stat)
!rm?  call memocc(i_stat ,PPD%iproj_to_ene,'iproj_to_ene',subname)
!rm?
!rm?  allocate(PPD%iproj_to_factor(mprojtot +ndebug ) , stat=i_stat)
!rm?  call memocc(i_stat ,PPD%iproj_to_factor,'iproj_to_factor',subname)
!rm?
!rm?  allocate(PPD%iproj_to_l(mprojtot +ndebug ) , stat=i_stat)
!rm?  call memocc(i_stat ,PPD%iproj_to_l,'iproj_to_l',subname)
!rm?
!rm?  PPD%mprojtot=mprojtot
!rm?
!rm?
!rm?  startjorb=1
!rm?  jorb=1
!rm?  istart_c=1
!rm?  Gocc(:)=0.0_wp
!rm?
!rm?
!rm?  iproj=0
!rm?  do iat=1,at%astruct%nat
!rm?
!rm?     mproj=0
!rm?     do while( jorb<=PPD%G%ncoeff         .and. PPD%iorbtolr(jorb)== iat)
!rm?        if( PPD%gaenes(jorb)<ecut_pc) then
!rm?           mproj=mproj+1
!rm?        endif
!rm?        if(jorb==PPD%G%ncoeff) exit
!rm?        jorb=jorb+1
!rm?     end do
!rm?
!rm?     PPD%ilr_to_mproj(iat)=mproj
!rm?     if( mproj>0) then
!rm?        nvctr_c  =PPD%pc_nl%pspd(iat)%plr%wfd%nvctr_c!PPD%pc_nlpspd%nvctr_p(2*iat-1)-PPD%pc_nlpspd%nvctr_p(2*iat-2)
!rm?        nvctr_f  =PPD%pc_nl%pspd(iat)%plr%wfd%nvctr_f!PPD%pc_nlpspd%nvctr_p(2*iat  )-PPD%pc_nlpspd%nvctr_p(2*iat-1)
!rm?
!rm?        jorb=startjorb
!rm?        do while( jorb<=PPD%G%ncoeff         .and. PPD%iorbtolr(jorb)== iat) 
!rm?           if( PPD%gaenes(jorb)<ecut_pc) then
!rm?              iproj=iproj+1
!rm?              PPD%iproj_to_ene(iproj) = PPD%gaenes(jorb)
!rm?              PPD%iproj_to_l(iproj)   = iorbto_l(jorb)
!rm?
!rm?              istart_c=istart_c + (   nvctr_c    +   7*nvctr_f   )
!rm?           endif
!rm?           jorb=jorb+1
!rm?
!rm?           if(jorb> PPD%G%ncoeff) exit
!rm?        end do
!rm?
!rm?
!rm?        if( .not. PPD%DistProjApply) then
!rm?           istart_c= istart_c-mproj*(nvctr_c+7*nvctr_f)
!rm?
!rm?           call fillPcProjOnTheFly(PPD, Glr, iat, at, hx,hy,hz, startjorb,ecut_pc ,  istart_c ) 
!rm?           istart_c= istart_c+mproj*(nvctr_c+7*nvctr_f)
!rm?
!rm?!!$
!rm?!!$           ncplx=1
!rm?!!$           rdum=0.0_gp
!rm?!!$
!rm?!!$           mbvctr_c=PPD%pc_nlpspd%nvctr_p(2*iat-1)-PPD%pc_nlpspd%nvctr_p(2*iat-2)
!rm?!!$           mbvctr_f=PPD%pc_nlpspd%nvctr_p(2*iat  )-PPD%pc_nlpspd%nvctr_p(2*iat-1)
!rm?!!$           
!rm?!!$           mbseg_c=PPD%pc_nlpspd%nseg_p(2*iat-1)-PPD%pc_nlpspd%nseg_p(2*iat-2)
!rm?!!$           mbseg_f=PPD%pc_nlpspd%nseg_p(2*iat  )-PPD%pc_nlpspd%nseg_p(2*iat-1)
!rm?!!$           jseg_c=PPD%pc_nlpspd%nseg_p(2*iat-2)+1
!rm?!!$              
!rm?!!$           do idum=1, 9
!rm?!!$              call wpdot_wrap(ncplx,  &
!rm?!!$                   mbvctr_c,mbvctr_f,mbseg_c,mbseg_f,PPD%pc_nlpspd%keyv_p(jseg_c),&
!rm?!!$                   PPD%pc_nlpspd%keyg_p(1,jseg_c),PPD%pc_proj(istart_c-idum*(nvctr_c+7*nvctr_f)),& 
!rm?!!$                   mbvctr_c,mbvctr_f,mbseg_c,mbseg_f,PPD%pc_nlpspd%keyv_p(jseg_c),&
!rm?!!$                   PPD%pc_nlpspd%keyg_p(1,jseg_c),&
!rm?!!$                   PPD%pc_proj(istart_c-idum*(nvctr_c+7*nvctr_f)),&
!rm?!!$                   rdum)
!rm?!!$           end do
!rm?        endif
!rm?
!rm?     end if
!rm?
!rm?     !! aggiunger condizione su istartc_c per vedere se e nelprj
!rm?
!rm?     startjorb=jorb
!rm?
!rm?  enddo
!rm?
!rm?  if( .not. PPD%DistProjApply) then
!rm?     call deallocate_gwf(PPD%G,subname)
!rm?  endif
!rm?
!rm?  i_all=-product(shape(logrid))*kind(logrid)
!rm?  deallocate(logrid,stat=i_stat)
!rm?  call memocc(i_stat,i_all,'logrid',subname)
!rm?
!rm?  i_all=-product(shape(Gocc))*kind(Gocc)
!rm?  deallocate(Gocc,stat=i_stat)
!rm?  call memocc(i_stat,i_all,'Gocc',subname)
!rm?
!rm?!!$  i_all=-product(shape(iorbtolr))*kind(iorbtolr)
!rm?!!$  deallocate(iorbtolr,stat=i_stat)
!rm?!!$  call memocc(i_stat,i_all,'iorbtolr',subname)
!rm?
!rm?  i_all=-product(shape(iorbto_l))*kind(iorbto_l)
!rm?  deallocate(iorbto_l,stat=i_stat)
!rm?  call memocc(i_stat,i_all,'iorbto_l',subname)
!rm?
!rm?  i_all=-product(shape(iorbto_m))*kind(iorbto_m)
!rm?  deallocate(iorbto_m,stat=i_stat)
!rm?  call memocc(i_stat,i_all,'iorbto_m',subname)
!rm?
!rm?  i_all=-product(shape(iorbto_ishell))*kind(iorbto_ishell)
!rm?  deallocate(iorbto_ishell,stat=i_stat)
!rm?  call memocc(i_stat,i_all,'iorbto_ishell',subname)
!rm?
!rm?
!rm?  i_all=-product(shape(iorbto_iexpobeg))*kind(iorbto_iexpobeg)
!rm?  deallocate(iorbto_iexpobeg,stat=i_stat)
!rm?  call memocc(i_stat,i_all,'iorbto_iexpobeg',subname)
!rm?
!rm?END SUBROUTINE createPcProjectorsArrays
!rm?

!> Determine localization region for all preconditioning projectors, but do not yet fill the descriptor arrays
subroutine createPawProjectorsArrays(iproc,n1,n2,n3,rxyz,at,orbs,&
     &   radii_cf,cpmult,fpmult,hx,hy,hz, &
     &   PAWD, Glr)
  use module_interfaces
  use module_base
  use module_types
  use module_abscalc
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
  !start from a null structure
  PAWD%paw_nl=DFT_PSP_projectors_null()

  PAWD%paw_nl%natoms=PAWD%G%nat
  allocate(PAWD%paw_nl%pspd(PAWD%paw_nlpspd%natoms))
  do iat=1,PAWD%paw_nl%natoms
     PAWD%paw_nl%pspd(iat)=nonlocal_psp_descriptors_null()
  end do

  allocate(logrid(0:n1,0:n2,0:n3+ndebug),stat=i_stat)
  call memocc(i_stat,logrid,'logrid',subname)

  call localize_projectors_paw(iproc,n1,n2,n3,hx,hy,hz,Pcpmult,1*fpmult,rxyz,radii_cf,&
       &   logrid,at,orbs,PAWD)

  ! the above routine counts atomic projector and the number of their element for psp
  ! We must therefore correct , later, nlpspd%nprojel  and nlpspd%nproj
  !-------------------

  ! allocations for arrays holding the projectors and their data descriptors
  do iat=1,PAWD%paw_nl%natoms
     !for the moments the bounds are not needed for projectors
     call allocate_wfd(PAWD%paw_nl%pspd(iat)%plr%wfd)
  end do
  PAWD%paw_nl%proj=f_malloc0_ptr(PAWD%paw_nl%nprojel,id='paw_proj')

!!$  allocate(PAWD%paw_proj(PAWD%paw_nlpspd%nprojel+ndebug),stat=i_stat)
!!$  call memocc(i_stat,PAWD%paw_proj,'paw_proj',subname)

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

           call bounds_to_plr_limits(.false.,1,PAWD%paw_nl%pspd(iat)%plr,nl1,nl2,nl3,nu1,nu2,nu3)

           call fill_logrid(at%astruct%geocode,n1,n2,n3,nl1,nu1,nl2,nu2,nl3,nu3,0,1,  &
                &   at%astruct%ntypes,at%astruct%iatype(iatat),rxyz(1,iatat),radii_cf(1,3),Pcpmult,hx,hy,hz,logrid)
2)
           mseg=PAWD%paw_nl%pspd(iat)%plr%wfd%nseg_c

           call segkeys(n1,n2,n3,nl1,nu1,nl2,nu2,nl3,nu3,  & 
                logrid,mseg,&
                PAWD%paw_nl%pspd(iat)%plr%wfd%keyglob(1,1),&
                PAWD%paw_nl%pspd(iat)%plr%wfd%keyvglob(1))

           mvctr =PAWD%paw_nl%pspd(iat)%plr%wfd%nvctr_c!PAWD%paw_nlpspd%nvctr_p(2*iat-1)-PAWD%paw_nlpspd%nvctr_p(2*iat-2)

           call bounds_to_plr_limits(.false.,2,PAWD%paw_nl%pspd(iat)%plr,&
                nl1,nl2,nl3,nu1,nu2,nu3)
           call fill_logrid(at%astruct%geocode,n1,n2,n3,nl1,nu1,nl2,nu2,nl3,nu3,0,1,  &
                &   at%astruct%ntypes,at%astruct%iatype(iatat),rxyz(1,iatat),radii_cf(1,2),1*fpmult,hx,hy,hz,logrid)

           iseg=PAWD%paw_nl%pspd(iat)%plr%wfd%nseg_c+1
           mseg=PAWD%paw_nl%pspd(iat)%plr%wfd%nseg_f

           if (mseg > 0) then
              call segkeys(n1,n2,n3,nl1,nu1,nl2,nu2,nl3,nu3,  & 
                   logrid,mseg,&
                   PAWD%paw_nl%pspd(iat)%plr%wfd%keyglob(1,iseg),&
                   PAWD%paw_nl%pspd(iat)%plr%wfd%keyvglob(iseg))

              mvctr =PAWD%paw_nl%pspd(iat)%plr%wfd%nvctr_f!PAWD%paw_nlpspd%nvctr_p(2*iat)-PAWD%paw_nlpspd%nvctr_p(2*iat-1)
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
              nvctr_c  =PAWD%paw_nl%pspd(iat)%plr%wfd%nvctr_c!PAWD%paw_nlpspd%nvctr_p(2*iat-1)-PAWD%paw_nlpspd%nvctr_p(2*iat-2)
              nvctr_f  =PAWD%paw_nl%pspd(iat)%plr%wfd%nvctr_f!PAWD%paw_nlpspd%nvctr_p(2*iat  )-PAWD%paw_nlpspd%nvctr_p(2*iat-1)

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

subroutine localize_projectors_paw(iproc,n1,n2,n3,hx,hy,hz,cpmult,fpmult,rxyz,radii_cf,&
     logrid,at,orbs,PAWD)
  use module_base
  use module_types
  use module_abscalc
  implicit none
  integer, intent(in) :: iproc,n1,n2,n3
  real(gp), intent(in) :: cpmult,fpmult,hx,hy,hz
  type(atoms_data), intent(in) :: at
  type(orbitals_data), intent(in) :: orbs

  real(gp), dimension(3,at%astruct%nat), intent(in) :: rxyz
  real(gp), dimension(at%astruct%ntypes,3), intent(in) :: radii_cf
  logical, dimension(0:n1,0:n2,0:n3), intent(inout) :: logrid
  type(PAWproj_data_type) ::PAWD

  !Local variables
  integer :: istart,ityp,natyp,iat,mproj,nl1,nu1,nl2,nu2,nl3,nu3,mvctr,mseg,nprojelat,i,l
  integer :: ikpt,nkptsproj,ikptp,i_stat
  real(gp) :: maxfullvol,totfullvol,totzerovol,zerovol,fullvol,maxrad,maxzerovol,rad
  integer :: natpaw

  if (iproc.eq.0) then
     write(*,'(1x,a)')&
          '------------------------------------------------------------ PSP Projectors Creation'
     write(*,'(1x,a4,4x,a4,2(1x,a))')&
          'Type','Name','Number of atoms','Number of paw projectors per atom'
  end if
  
!!$  PAWD%paw_nlpspd%nseg_p(0)=0 
!!$  PAWD%paw_nlpspd%nvctr_p(0)=0 

  istart=1
  PAWD%paw_nlpspd%nproj=0
  PAWD%paw_nlpspd%nprojel=0

  if (iproc ==0) then
     !print the number of projectors to be created
     do ityp=1,at%astruct%ntypes
        natyp=0
        mproj=0
        if(at%paw_NofL(ityp).gt.0) then
           do iat=1,at%astruct%nat
              if (at%astruct%iatype(iat) == ityp) then
                 if(natyp.eq.0) then
                    call numb_proj_paw(ityp,mproj)                    
                 endif
                 natyp=natyp+1
              endif
           end do
           write(*,'(1x,i4,2x,a6,1x,i15,i21)')&
                ityp,trim(at%astruct%atomnames(ityp)),natyp,mproj
        end if
     end do
  end if

  !count number of PAW projectors
  natpaw=0
  do iat=1,at%astruct%nat
     if(  at%paw_NofL(at%astruct%iatype(iat)).gt.0) then
        call numb_proj_paw(at%astruct%iatype(iat),mproj)
        if (mproj /= 0) then 
           natpaw=natpaw+1
        end if
     end if
  end do
  PAWD%paw_nlpspd%natoms=natpaw
  allocate(PAWD%paw_nlpspd%plr(PAWD%paw_nlpspd%natoms),stat=i_stat)

  natpaw=0
  do iat=1,at%astruct%nat

     if(  at%paw_NofL(at%astruct%iatype(iat)).gt.0) then

        call numb_proj_paw(at%astruct%iatype(iat),mproj)

        if (mproj /= 0) then 
           natpaw=natpaw+1
           PAWD%paw_nlpspd%nproj=PAWD%paw_nlpspd%nproj+mproj



           ! coarse grid quantities
           call pregion_size(at%astruct%geocode,rxyz(1,iat),radii_cf(at%astruct%iatype(iat),3),cpmult, &
                hx,hy,hz,n1,n2,n3,nl1,nu1,nl2,nu2,nl3,nu3)

           PAWD%paw_nlpspd%plr(natpaw)%ns1=nl1     
           PAWD%paw_nlpspd%plr(natpaw)%ns2=nl2                                   
           PAWD%paw_nlpspd%plr(natpaw)%ns3=nl3                                   
                                                
           PAWD%paw_nlpspd%plr(natpaw)%d%n1=nu1-nl1
           PAWD%paw_nlpspd%plr(natpaw)%d%n2=nu2-nl2
           PAWD%paw_nlpspd%plr(natpaw)%d%n3=nu3-nl3

           call fill_logrid(at%astruct%geocode,n1,n2,n3,nl1,nu1,nl2,nu2,nl3,nu3,0,1,  &
                at%astruct%ntypes,at%astruct%iatype(iat),rxyz(1,iat),radii_cf(1,3),cpmult,hx,hy,hz,logrid)
           call num_segkeys(n1,n2,n3,nl1,nu1,nl2,nu2,nl3,nu3,logrid,mseg,mvctr)

           PAWD%paw_nlpspd%plr(natpaw)%wfd%nseg_c=mseg
           PAWD%paw_nlpspd%plr(natpaw)%wfd%nvctr_c=mvctr

           istart=istart+mvctr*mproj

           nprojelat=mvctr*mproj

           !print *,'iat,mvctr',iat,mvctr,mseg,mproj

           ! fine grid quantities

           call pregion_size(at%astruct%geocode,rxyz(1,iat),radii_cf(at%astruct%iatype(iat),2),&
                fpmult,&
                hx,hy,hz,n1,n2,n3,nl1,nu1,nl2,nu2,nl3,nu3)

           PAWD%paw_nlpspd%plr(natpaw)%d%nfl1=nl1-&
                PAWD%paw_nlpspd%plr(natpaw)%ns1
           PAWD%paw_nlpspd%plr(natpaw)%d%nfl2=nl2-&
                PAWD%paw_nlpspd%plr(natpaw)%ns2
           PAWD%paw_nlpspd%plr(natpaw)%d%nfl3=nl3-&
                PAWD%paw_nlpspd%plr(natpaw)%ns3
           PAWD%paw_nlpspd%plr(natpaw)%d%nfu1=nu1-&
                PAWD%paw_nlpspd%plr(natpaw)%ns1
           PAWD%paw_nlpspd%plr(natpaw)%d%nfu2=nu2-&
                PAWD%paw_nlpspd%plr(natpaw)%ns2
           PAWD%paw_nlpspd%plr(natpaw)%d%nfu3=nu3-&
                PAWD%paw_nlpspd%plr(natpaw)%ns3

           call fill_logrid(at%astruct%geocode,n1,n2,n3,nl1,nu1,nl2,nu2,nl3,nu3,0,1,  &
                at%astruct%ntypes,at%astruct%iatype(iat),rxyz(1,iat),radii_cf(1,2),fpmult,hx,hy,hz,logrid)
           call num_segkeys(n1,n2,n3,nl1,nu1,nl2,nu2,nl3,nu3,logrid,mseg,mvctr)
           !if (iproc.eq.0) write(*,'(1x,a,2(1x,i0))') 'mseg,mvctr, fine  projectors ',mseg,mvctr
           PAWD%paw_nlpspd%plr(natpaw)%wfd%nseg_f=mseg
           PAWD%paw_nlpspd%plr(natpaw)%wfd%nvctr_f=mvctr

           istart=istart+7*mvctr*mproj
           nprojelat=nprojelat+7*mvctr*mproj

           PAWD%paw_nlpspd%nprojel=max(PAWD%paw_nlpspd%nprojel,nprojelat)

           !print *,'iat,nprojelat',iat,nprojelat,mvctr,mseg

        else  !(atom has no nonlocal PSP, e.g. H)

           PAWD%paw_nlpspd%plr(natpaw)%wfd%nseg_c=0
           PAWD%paw_nlpspd%plr(natpaw)%wfd%nvctr_c=0
           PAWD%paw_nlpspd%plr(natpaw)%wfd%nseg_f=0
           PAWD%paw_nlpspd%plr(natpaw)%wfd%nvctr_f=0

           !! the following is necessary to the creati of preconditioning projectors

           ! coarse grid quantities
           call pregion_size(at%astruct%geocode,rxyz(1,iat),radii_cf(at%astruct%iatype(iat),3),cpmult, &
                hx,hy,hz,n1,n2,n3,nl1,nu1,nl2,nu2,nl3,nu3)
           
           PAWD%paw_nlpspd%plr(natpaw)%ns1=nl1     
           PAWD%paw_nlpspd%plr(natpaw)%ns2=nl2                                   
           PAWD%paw_nlpspd%plr(natpaw)%ns3=nl3                                   

           PAWD%paw_nlpspd%plr(natpaw)%d%n1=nu1-nl1
           PAWD%paw_nlpspd%plr(natpaw)%d%n2=nu2-nl2
           PAWD%paw_nlpspd%plr(natpaw)%d%n3=nu3-nl3

           ! fine grid quantities
           call pregion_size(at%astruct%geocode,rxyz(1,iat),radii_cf(at%astruct%iatype(iat),2),fpmult,&
                hx,hy,hz,n1,n2,n3,nl1,nu1,nl2,nu2,nl3,nu3)

           PAWD%paw_nlpspd%plr(natpaw)%d%nfl1=nl1-&
                PAWD%paw_nlpspd%plr(natpaw)%ns1
           PAWD%paw_nlpspd%plr(natpaw)%d%nfl2=nl2-&
                PAWD%paw_nlpspd%plr(natpaw)%ns2
           PAWD%paw_nlpspd%plr(natpaw)%d%nfl3=nl3-&
                PAWD%paw_nlpspd%plr(natpaw)%ns3
           PAWD%paw_nlpspd%plr(natpaw)%d%nfu1=nu1-&
                PAWD%paw_nlpspd%plr(natpaw)%ns1
           PAWD%paw_nlpspd%plr(natpaw)%d%nfu2=nu2-&
                PAWD%paw_nlpspd%plr(natpaw)%ns2
           PAWD%paw_nlpspd%plr(natpaw)%d%nfu3=nu3-&
                PAWD%paw_nlpspd%plr(natpaw)%ns3

        endif
     endif
  enddo
  
  
  !   if (memorylimit /= 0.e0 .and. .not. DistProjApply .and. &
  !        real(istart-1,kind=4) > memorylimit*134217728.0e0) then
  !      if (iproc == 0) then
  !         write(*,'(44x,a)') '------ On-the-fly paw projectors application'
  !      end if
  !      DistProjApply =.true.
  !   end if
  
  !calculate the fraction of the projector array used for allocate zero values
  !control the hardest and the softest gaussian
  totzerovol=0.0_gp
  maxfullvol=0.0_gp
  totfullvol=0.0_gp
  do iat=1,at%astruct%nat
     if(  at%paw_NofL(at%astruct%iatype(iat)).gt.0) then
        ityp=at%astruct%iatype(iat)
        maxrad=min(maxval(at%psppar(1:4,0,ityp)),cpmult/15.0_gp*radii_cf(ityp,3))
        zerovol=0.0_gp
        fullvol=0.0_gp
        do l=1,4
           do i=1,3
              if (at%psppar(l,i,ityp) /= 0.0_gp) then
                 rad=min(at%psppar(l,0,ityp),cpmult/15.0_gp*radii_cf(ityp,3))
                 zerovol=zerovol+(maxrad**3-rad**3)
                 fullvol=fullvol+maxrad**3
              end if
           end do
        end do
        if (fullvol >= maxfullvol .and. fullvol > 0.0_gp) then
           maxzerovol=zerovol/fullvol
           maxfullvol=fullvol
        end if
        totzerovol=totzerovol+zerovol
        totfullvol=totfullvol+fullvol
     endif
  end do

  !assign the total quantity per atom
  zerovol=0.d0
  if (totfullvol /= 0.0_gp) then
     if (PAWD%DistProjApply) then
        zerovol=maxzerovol
     else
        zerovol=totzerovol/totfullvol
     end if
  end if

  !here is the point in which the projector strategy should be decided
  !DistProjApply shoud never change after this point

  !number of elements of the projectors
  if (.not. PAWD%DistProjApply) PAWD%paw_nlpspd%nprojel=istart-1

  nkptsproj=1
  if ((.not.PAWD%DistProjApply) .and. orbs%norbp > 0) then
     nkptsproj = 0
     !the new solution did not work when there is no orbital on the processor
     do ikptp=1,orbs%nkptsp! orbs%iokpt(1), orbs%iokpt(orbs%norbp)
        ikpt=orbs%iskpts+ikptp
!!$         print *, " k points ", orbs%kpts

        if (orbs%kpts(1,ikpt)**2+orbs%kpts(2,ikpt)**2+orbs%kpts(3,ikpt)**2 >0 .and. &
             &  orbs%nspinor > 1) then
           nkptsproj = nkptsproj + 2
        else
           nkptsproj = nkptsproj + 1
        end if
     end do
  else if (PAWD%DistProjApply) then
     !the new solution did not work when there is no orbital on the processor
     do ikptp=1,orbs%nkptsp! orbs%iokpt(1), orbs%iokpt(orbs%norbp)
        ikpt=orbs%iskpts+ikptp
        if (orbs%kpts(1,ikpt)**2+orbs%kpts(2,ikpt)**2+orbs%kpts(3,ikpt)**2 >0 .and. &
             &  orbs%nspinor > 1) then
           nkptsproj = max(nkptsproj, 2)
        end if
     end do
  end if
  !   print *, " nkptsproj EST    ", nkptsproj
  !   print *, " PAWD%paw_nlpspd%nprojel EST  ", PAWD%paw_nlpspd%nprojel

  PAWD%paw_nlpspd%nprojel=nkptsproj*PAWD%paw_nlpspd%nprojel
  if (iproc == 0) then
     if (PAWD%DistProjApply) then
        write(*,'(44x,a)') '------  PAWD: On-the-fly projectors application'
     else
        write(*,'(44x,a)') '------'
     end if
     write(*,'(1x,a,i21)') 'Total number of projectors =',PAWD%paw_nlpspd%nproj
     write(*,'(1x,a,i21)') 'Total number of components =',PAWD%paw_nlpspd%nprojel
     write(*,'(1x,a,i21)') 'Percent of zero components =',nint(100.0_gp*zerovol)
  end if
contains
  
subroutine numb_proj_paw(ityp,mproj)
  integer , intent(in):: ityp
  integer, intent(out):: mproj
  

  integer :: il,jtyp

  mproj=0
  il=0
  do jtyp=1,ityp-1
     il=il+at%paw_NofL(jtyp)
  enddo
  do i =1, at%paw_NofL(ityp)
     il=il+1
     if( at%paw_l(il).ge.0) then
        mproj=mproj+at%paw_nofchannels(il)*(2*at%paw_l(il) +1)
     else
        mproj=mproj+at%paw_nofchannels(il)*(-2*at%paw_l(il) -1)        
     endif
  enddo
end subroutine numb_proj_paw

END subroutine localize_projectors_paw
