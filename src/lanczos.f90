!> @file
!!  Lanczos diagonalisation used by XANES calculation
!! @author
!!    Copyright (C) 2009-2011 BigDFT group (AM, LG)
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 

!> Lanczos diagonalization
subroutine xabs_lanczos(iproc,nproc,at,hx,hy,hz,rxyz,&
     radii_cf,nlpspd,proj,Lzd,dpcom,potential,&
     energs,nspin,GPU,in_iat_absorber,&
     in , PAWD , orbs )! add to interface
   use module_base
   use module_types
   use lanczos_interface
   use lanczos_base
   use module_interfaces ,except_this_one => xabs_lanczos
   implicit none
   integer, intent(in) :: iproc,nproc,nspin
   real(gp), intent(in) :: hx,hy,hz
   type(atoms_data), intent(in), target :: at
   type(nonlocal_psp_descriptors), intent(in), target :: nlpspd
   type(local_zone_descriptors), intent(inout), target :: Lzd
   type(denspot_distribution), intent(in), target :: dpcom
   real(gp), dimension(3,at%astruct%nat), intent(in), target :: rxyz
   real(gp), dimension(at%astruct%ntypes,3), intent(in), target ::  radii_cf
   real(wp), dimension(nlpspd%nprojel), intent(in), target :: proj
   real(wp), dimension(max(dpcom%ndimpot,1),nspin), target :: potential
   type(energy_terms), intent(inout) :: energs
   type(GPU_pointers), intent(inout) , target :: GPU
   integer, intent(in) :: in_iat_absorber

   type(input_variables),intent(in), target :: in
   type(pawproj_data_type), target ::PAWD
   type(orbitals_data), intent(inout), target :: orbs

   !local variables
   character(len=*), parameter :: subname='lanczos'
   integer :: i_stat,i_all
   type(lanczos_args) :: ha
   integer :: i

   real(wp), pointer :: Gabs_coeffs(:)
   real(wp), dimension(:), pointer :: pot

   if(iproc==0) write(*,*) " IN ROUTINE LANCZOS "

   if (GPUconv) then
      call prepare_gpu_for_locham(Lzd%Glr%d%n1,Lzd%Glr%d%n2,Lzd%Glr%d%n3,in%nspin,&
         &   hx,hy,hz,Lzd%Glr%wfd,orbs,GPU)
   end if
   GPU%full_locham=.true.
   if (OCLconv) then
      call allocate_data_OCL(Lzd%Glr%d%n1,Lzd%Glr%d%n2,Lzd%Glr%d%n3,at%astruct%geocode,&
         &   in%nspin,Lzd%Glr%wfd,orbs,GPU)
      if (iproc == 0) write(*,*) 'GPU data allocated'
   end if


   allocate(orbs%eval(orbs%norb+ndebug),stat=i_stat)
   call memocc(i_stat,orbs%eval,'orbs%eval',subname)
   orbs%occup(1:orbs%norb)=1.0_gp
   orbs%spinsgn(1:orbs%norb)=1.0_gp
   orbs%eval(1:orbs%norb)=1.0_gp
   !call allocate_comms(nproc,ha%comms,subname)
   call orbitals_communicators(iproc,nproc,Lzd%Glr,orbs,ha%comms)  

   call local_potential_dimensions(Lzd,orbs,dpcom%ngatherarr(0,1))

   allocate(Gabs_coeffs(2*in%L_absorber+1+ndebug),stat=i_stat)
   call memocc(i_stat,Gabs_coeffs,'Gabs_coeffs',subname)


   if(   at%paw_NofL( at%astruct%iatype(   in_iat_absorber ) ) .gt. 0   ) then     
      Gabs_coeffs(:)=in%Gabs_coeffs(:)
   else     
      print * ," You are asking for a spectra for atom " , in_iat_absorber
      print *, " but at%paw_NofL( at%astruct%iatype(   in_iat_absorber ) )=0 " 
      print *, " this mean that the pseudopotential file is not pawpatched. "
      print *, " You'll have to generated the patch with pseudo"
      STOP     
   endif

   call full_local_potential(iproc,nproc,orbs,Lzd,0,dpcom,potential,pot)

!!$   call full_local_potential(iproc,nproc,ndimpot,Lzd%Glr%d%n1i*Lzd%Glr%d%n2i*Lzd%Glr%d%n3i,&
!!$        in%nspin,Lzd%Glr%d%n1i*Lzd%Glr%d%n2i*Lzd%Glr%d%n3i*in%nspin,0,&
!!$        orbs,Lzd,0,ngatherarr,potential,pot)

   ha%in_iat_absorber=in_iat_absorber
   ha%Labsorber  = in%L_absorber
   ha%iproc=iproc
   ha%nproc=nproc
   ha%at=>at !!
   ha%hx=hx 
   ha%hy=hy
   ha%hz=hz
   ha%rxyz=>rxyz

   ha%radii_cf=>radii_cf
   ha%nlpspd=>nlpspd !!
   ha%proj=>proj !!
   ha%Lzd=>Lzd !!!
   ha%ngatherarr=>dpcom%ngatherarr
   ha%ndimpot=dpcom%ndimpot
   ha%potential=>pot
   ha%energs = energs
   ha%nspin=nspin
   ha%GPU=>GPU !!
   ha%Gabs_coeffs=>Gabs_coeffs
   ha%PAWD=> PAWD
   ha%SIC=>in%SIC
   ha%orbs=>orbs

   call EP_inizializza(ha) 

   if(.true.) then
      LB_nsteps =in%nsteps
      call LB_allocate_for_lanczos( )
      call EP_allocate_for_eigenprob(LB_nsteps)
      call EP_make_dummy_vectors(10)


      call LB_passeggia(0,LB_nsteps,      get_EP_dim, EP_initialize_start , EP_normalizza,&
         &   EP_Moltiplica, EP_GramSchmidt ,EP_set_all_random, EP_copy,   EP_mat_mult, &
         &   EP_scalare,EP_add_from_vect_with_fact     )


      if(ha%iproc==0) then
         open(unit=22,file="alphabeta")
         write(22,*) LB_nsteps, EP_norma2_initialized_state
         print *, " alpha and beta from lanczos "
         WRITE(*,'(I5,1000(1ES23.16))')    LB_nsteps, EP_norma2_initialized_state
         do i=0, LB_nsteps-1
            write(22,*) LB_alpha(i), LB_beta(i)
            WRITE(*,'(2ES23.16)')  LB_alpha(i), LB_beta(i)
         enddo

         close(unit=22)
      endif

      call LB_de_allocate_for_lanczos( )

   endif

   call deallocate_comms(ha%comms,subname)

   call EP_free(ha%iproc)

   if (GPUconv) then
      call free_gpu(GPU,orbs%norbp)
   else if (OCLconv) then
      call free_gpu_OCL(GPU,orbs,in%nspin)
   end if

   i_all=-product(shape(orbs%eval))*kind(orbs%eval)
   deallocate(orbs%eval,stat=i_stat)
   call memocc(i_stat,i_all,'orbs%eval',subname)

   i_all=-product(shape(Gabs_coeffs))*kind(Gabs_coeffs)
   deallocate(Gabs_coeffs,stat=i_stat)
   call memocc(i_stat,i_all,'Gabs_coeffs',subname)

   call free_full_potential(dpcom%mpi_env%nproc,0,pot,subname)

END SUBROUTINE xabs_lanczos


!> Chebychev polynomials to calculate the density of states
subroutine xabs_chebychev(iproc,nproc,at,hx,hy,hz,rxyz,&
     radii_cf,nlpspd,proj,Lzd,dpcom,potential,&
     energs,nspin,GPU,in_iat_absorber,in, PAWD , orbs  )

   use module_base
   use module_types
   use lanczos_interface
   use lanczos_base
   ! per togliere il bug 
   use module_interfaces, except_this_one => xabs_chebychev

   implicit none
   integer  :: iproc,nproc,nspin
   real(gp)  :: hx,hy,hz
   type(atoms_data), target :: at
   type(nonlocal_psp_descriptors), target :: nlpspd
   type(local_zone_descriptors), target :: Lzd
   type(denspot_distribution), intent(in), target :: dpcom
   real(gp), dimension(3,at%astruct%nat), target :: rxyz
   real(gp), dimension(at%astruct%ntypes,3), intent(in), target ::  radii_cf
   real(wp), dimension(nlpspd%nprojel), target :: proj
   real(wp), dimension(max(dpcom%ndimpot,1),nspin), target :: potential
   type(energy_terms), intent(inout) :: energs
   type(GPU_pointers), intent(inout) , target :: GPU
   integer, intent(in) :: in_iat_absorber 

   type(input_variables),intent(in), target :: in
   type(pawproj_data_type), target ::PAWD
   type(orbitals_data), intent(inout), target :: orbs

   !Local variables
   character(len=*), parameter :: subname='chebychev'
   integer :: i_stat,i_all
   type(lanczos_args) :: ha
   integer :: i

   real(wp), dimension(:), pointer  :: pot

   real(gp) :: eval_min, eval_max, fact_cheb, cheb_shift
   real(gp) :: Pi
   logical:: dopaw
   real(gp) :: GetBottom

   if (iproc==0) print *, " IN ROUTINE  chebychev  "

   Pi=acos(-1.0_gp)

   if (GPUconv) then
      call prepare_gpu_for_locham(Lzd%Glr%d%n1,Lzd%Glr%d%n2,Lzd%Glr%d%n3,in%nspin,&
         &   hx,hy,hz,Lzd%Glr%wfd,orbs,GPU)
   end if

   GPU%full_locham=.true.

   if (OCLconv) then
      call allocate_data_OCL(Lzd%Glr%d%n1,Lzd%Glr%d%n2,Lzd%Glr%d%n3,at%astruct%geocode,&
         &   in%nspin,Lzd%Glr%wfd,orbs,GPU)
      if (iproc == 0) write(*,*)&
         &   'GPU data allocated'
   end if

   allocate(orbs%eval(orbs%norb *orbs%nkpts  +ndebug),stat=i_stat)
   call memocc(i_stat,orbs%eval,'orbs%eval',subname)
   orbs%occup(1:orbs%norb*orbs%nkpts )=1.0_gp
   orbs%spinsgn(1:orbs%norb*orbs%nkpts )=1.0_gp
   orbs%eval(1:orbs%norb*orbs%nkpts )=1.0_gp

   call orbitals_communicators(iproc,nproc,Lzd%Glr,orbs,ha%comms)  

   call local_potential_dimensions(Lzd,orbs,dpcom%ngatherarr(0,1))

   if(   at%paw_NofL( at%astruct%iatype(   in_iat_absorber ) ) .gt. 0   ) then     
   else
      print * ," You are asking for a spactra for atom " , in_iat_absorber
      print *, " but at%paw_NofL( at%astruct%iatype(   in_iat_absorber ) )=0 " 
      print *, " this mean that the pseudopotential file is not pawpatched. "
      print *, " You'll have to generated the patch with pseudo"
      STOP     
   endif

   call full_local_potential(iproc,nproc,orbs,Lzd,0,dpcom,potential,pot)
!!$   call full_local_potential(iproc,nproc,ndimpot,Lzd%Glr%d%n1i*Lzd%Glr%d%n2i*Lzd%Glr%d%n3i,&
!!$        in%nspin,Lzd%Glr%d%n1i*Lzd%Glr%d%n2i*Lzd%Glr%d%n3i*in%nspin,0,&
!!$        orbs,Lzd,0,ngatherarr,potential,pot)

   !associate hamapp_arg pointers
   ha%in_iat_absorber=in_iat_absorber
   ha%Labsorber  = in%L_absorber
   ha%iproc=iproc
   ha%nproc=nproc
   ha%at=>at !!
   ha%hx=hx 
   ha%hy=hy
   ha%hz=hz
   ha%rxyz=>rxyz

   ha%radii_cf=>radii_cf
   ha%nlpspd=>nlpspd !!
   ha%proj=>proj !!
   ha%Lzd=>Lzd !!!
   ha%ngatherarr=>dpcom%ngatherarr
   ha%ndimpot=dpcom%ndimpot
   ha%potential=>pot
   ha%energs = energs
   ha%nspin=nspin
   ha%GPU=>GPU !!
   ha%Gabs_coeffs=>in%Gabs_coeffs
   ha%PAWD=> PAWD
   ha%SIC=>in%SIC
   ha%orbs=>orbs

   call EP_inizializza(ha)  

   !!$  if(.false.) then
   !!$
   !!$     ! trova il valore massimo 
   !!$     shift =-0.0
   !!$     tol   =1.0D-8
   !!$     accontentati_di=1
   !!$     
   !!$     cercacount = LB_cerca( 10, shift, tol, set_EP_shift, EP_allocate_for_eigenprob, EP_make_dummy_vectors, &
!!$          get_EP_dim, EP_initialize_start , EP_normalizza, EP_Moltiplica, EP_GramSchmidt, &
!!$          EP_set_all_random, EP_copy , EP_mat_mult,  EP_scalare,EP_add_from_vect_with_fact,accontentati_di)
   !!$     
   !!$     if(iproc==0) then
   !!$        print *, " maximal eigenvalues " 
   !!$        print *, LB_eval
   !!$     endif
   !!$     eval_max = LB_eval(0)
   !!$     
   !!$     ! trova il valore minimo 
   !!$     shift =-10000
   !!$     
   !!$     
   !!$     cercacount = LB_cerca( 10, shift, tol, set_EP_shift, EP_allocate_for_eigenprob, EP_make_dummy_vectors, &
         !!$          get_EP_dim, EP_initialize_start , EP_normalizza, EP_Moltiplica, EP_GramSchmidt, &
         !!$          EP_set_all_random, EP_copy , EP_mat_mult,  EP_scalare,EP_add_from_vect_with_fact,accontentati_di)
   !!$     
   !!$     if(iproc==0) then
   !!$        print *, " minima eigenvalues" 
   !!$        print *, LB_eval
   !!$     endif
   !!$     eval_min = LB_eval(0)+10000
   !!$
   !!$  else

   eval_min = GetBottom(  at, nspin)-1.0 - in%abscalc_bottomshift
   eval_max = 4.0*Pi*Pi*(1.0/hx/hx + 1.0/hy/hy + 1.0/hz/hz  )/2.0*1.1 +2

   !   endif

   cheb_shift=0.5*(eval_min+ eval_max) 
   fact_cheb = (2-0.0001)/(eval_max-eval_min)

   LB_nsteps = in%nsteps
   LB_norbp  = orbs%norbp
   LB_nproc=nproc
   LB_iproc=iproc

   call LB_allocate_for_chebychev( )
   call EP_allocate_for_eigenprob(6) 
   call EP_make_dummy_vectors(4)


   !! call set_EP_shift(-cheb_shift) 
   call set_EP_shift(0.0_gp) 

   if( sum( ha%at%paw_NofL ).gt.0 ) then
      dopaw=.true.
   else
      dopaw=.false.
   endif
   if(ha%iproc==0) then
      print *, "weigths ", orbs%kwgts
   endif

   call LB_passeggia_Chebychev (LB_nsteps, cheb_shift,  fact_cheb,     get_EP_dim, EP_initialize_start , EP_normalizza,&
      &   EP_Moltiplica, EP_GramSchmidt ,EP_set_all_random, EP_copy,   EP_mat_mult, &
      &   EP_scalare_multik,EP_add_from_vect_with_fact  , EP_multbyfact, EP_ApplySinv, EP_ApplyS, dopaw, &
      &   in%abscalc_S_do_cg,  in%abscalc_Sinv_do_cg, in%xabs_res_prefix, orbs%nkpts , orbs%norb_par, &
      &   orbs%kwgts(   orbs%iskpts+1    :    orbs%iskpts + orbs%norb_par(ha%iproc,0))   )

   if(ha%iproc==0) then
      print *, "coefficients from Chebychev "
      WRITE(*,'(I5,2ES23.16)')   2*LB_nsteps, cheb_shift,  fact_cheb
      print *,"... " 
      do i=0, 2*LB_nsteps-1
         if(i>2*LB_nsteps-1 -10) then
            WRITE(*,'(1ES23.16)')   LB_alpha_cheb(1, i)
         endif
      enddo
   endif

   call free_full_potential(dpcom%mpi_env%nproc,0,pot,subname)
   nullify(ha%potential)


   !deallocate communication and orbitals descriptors
   call deallocate_comms(ha%comms,subname)

   call EP_free(ha%iproc)
   call  LB_de_allocate_for_cheb( )

   !!$ this free is already executed by bigdft
   !!$
   if (GPUconv) then
      call free_gpu(GPU,orbs%norbp)
   else if (OCLconv) then
      call free_gpu_OCL(GPU,orbs,in%nspin)
   end if

   i_all=-product(shape(orbs%eval))*kind(orbs%eval)
   deallocate(orbs%eval,stat=i_stat)
   call memocc(i_stat,i_all,'orbs%eval',subname)


   !!$  i_all=-product(shape(Gabsorber%nshell))*kind(Gabsorber%nshell)
   !!$  deallocate(Gabsorber%nshell,stat=i_stat)
   !!$  call memocc(i_stat,i_all,'Gabsorber%nshell',subname)
   !!$
   !!$  i_all=-product(shape(Gabsorber%nam))*kind(Gabsorber%nam)
   !!$  deallocate(Gabsorber%nam,stat=i_stat)
   !!$  call memocc(i_stat,i_all,'Gabsorber%nam',subname)
   !!$
   !!$  i_all=-product(shape(Gabsorber%ndoc))*kind(Gabsorber%ndoc)
   !!$  deallocate(Gabsorber%ndoc,stat=i_stat)
   !!$  call memocc(i_stat,i_all,'Gabsorber%ndoc',subname)
   !!$
   !!$  i_all=-product(shape(Gabsorber%xp))*kind(Gabsorber%xp)
   !!$  deallocate(Gabsorber%xp,stat=i_stat)
   !!$  call memocc(i_stat,i_all,'Gabsorber%xp',subname)
   !!$
   !!$  i_all=-product(shape(Gabsorber%psiat))*kind(Gabsorber%psiat)
   !!$  deallocate(Gabsorber%psiat,stat=i_stat)
   !!$  call memocc(i_stat,i_all,'Gabsorber%psiat',subname)
   !!$
   !!$  if( associated(Gabsorber%rxyz)) then
   !!$     i_all=-product(shape(Gabsorber%rxyz))*kind(Gabsorber%rxyz)
   !!$     deallocate(Gabsorber%rxyz,stat=i_stat)
   !!$     call memocc(i_stat,i_all,'Gabsorber%rxyz',subname)
   !!$  endif

END SUBROUTINE xabs_chebychev


!> Finds the spectra solving  (H-omega)x=b
subroutine xabs_cg(iproc,nproc,at,hx,hy,hz,rxyz,&
      &   radii_cf,nlpspd,proj,Lzd,dpcom,potential,&
      &   energs,nspin,GPU,in_iat_absorber,&
      &   in , rhoXanes, PAWD , PPD, orbs )
   use module_base
   use module_types
   use lanczos_interface
   use lanczos_base
   ! per togliere il bug 
   use module_interfaces ,except_this_one => xabs_cg

   implicit none

   integer  :: iproc,nproc,nspin
   real(gp)  :: hx,hy,hz
   type(atoms_data), target :: at
   type(nonlocal_psp_descriptors), target :: nlpspd
   type(local_zone_descriptors), target :: Lzd
   type(pcproj_data_type), target ::PPD
   type(denspot_distribution), intent(in), target :: dpcom
   real(gp), dimension(3,at%astruct%nat), target :: rxyz
   real(gp), dimension(at%astruct%ntypes,3), intent(in), target ::  radii_cf
   real(wp), dimension(nlpspd%nprojel), target :: proj
   real(wp), dimension(max(dpcom%ndimpot,1),nspin), target :: potential
   real(wp), dimension(max(dpcom%ndimpot,1),nspin), target :: rhoXanes
   type(energy_terms), intent(inout) :: energs
   type(GPU_pointers), intent(inout) , target :: GPU
   integer, intent(in) :: in_iat_absorber
   type(pawproj_data_type), target ::PAWD
   type(input_variables),intent(in), target :: in
   type(orbitals_data), intent(inout), target :: orbs

   !local variables
   character(len=*), parameter :: subname='xabs_cg'
   integer :: i_stat,i_all
   type(lanczos_args) :: ha
   integer :: i,j
   real(gp) Ene,gamma,  res 

   real(wp),   pointer  :: Gabs_coeffs(:)
   real(wp), dimension(:), pointer  :: pot

   logical:: useold
   real(gp) , pointer ::potentialclone(:,:)

   if( iand( in%potshortcut,16)>0) then
      allocate(potentialclone(max(dpcom%ndimpot,1),nspin+ndebug),stat=i_stat)
      call memocc(i_stat,potentialclone,'potentialclone',subname)
      potentialclone=potential
   endif

   if(iproc==0) print *, " IN ROUTINE xabs_cg "

   if (GPUconv) then
      call prepare_gpu_for_locham(Lzd%Glr%d%n1,Lzd%Glr%d%n2,Lzd%Glr%d%n3,in%nspin,&
         &   hx,hy,hz,Lzd%Glr%wfd,orbs,GPU)
   end if
   GPU%full_locham=.true.
   if (OCLconv) then
      call allocate_data_OCL(Lzd%Glr%d%n1,Lzd%Glr%d%n2,Lzd%Glr%d%n3,at%astruct%geocode,&
         &   in%nspin,Lzd%Glr%wfd,orbs,GPU)
      if (iproc == 0) write(*,*)&
         &   'GPU data allocated'
   end if

   allocate(orbs%eval(orbs%norb+ndebug),stat=i_stat)
   call memocc(i_stat,orbs%eval,'orbs%eval',subname)

   orbs%occup(1:orbs%norb)=1.0_gp
   orbs%spinsgn(1:orbs%norb)=1.0_gp
   orbs%eval(1:orbs%norb)=1.0_gp
   !call allocate_comms(nproc,ha%comms,subname)
   call orbitals_communicators(iproc,nproc,Lzd%Glr,orbs,ha%comms)  

   call local_potential_dimensions(Lzd,orbs,dpcom%ngatherarr(0,1))

   allocate(Gabs_coeffs(2*in%L_absorber+1+ndebug),stat=i_stat)
   call memocc(i_stat,Gabs_coeffs,'Gabs_coeffs',subname)

   if(   at%paw_NofL( at%astruct%iatype(   in_iat_absorber ) ) .gt. 0   ) then     
      Gabs_coeffs(:)=in%Gabs_coeffs(:)
   else
      print * ," You are asking for a spactra for atom " , in_iat_absorber
      print *, " but at%paw_NofL( at%astruct%iatype(   in_iat_absorber ) )=0 " 
      print *, " this mean that the pseudopotential file is not pawpatched. "
      print *, " You'll have to generated the patch with pseudo"
      STOP     
   endif
   
   call full_local_potential(iproc,nproc,orbs,Lzd,0,dpcom,potential,pot)
!!$   call full_local_potential(iproc,nproc,ndimpot,Lzd%Glr%d%n1i*Lzd%Glr%d%n2i*Lzd%Glr%d%n3i,&
!!$        in%nspin,Lzd%Glr%d%n1i*Lzd%Glr%d%n2i*Lzd%Glr%d%n3i*in%nspin,0,&
!!$        orbs,Lzd,0,ngatherarr,potential,pot)

   ha%in_iat_absorber=in_iat_absorber
   ha%Labsorber  = in%L_absorber
   ha%iproc=iproc
   ha%nproc=nproc
   ha%at=>at !!
   ha%hx=hx 
   ha%hy=hy
   ha%hz=hz
   ha%rxyz=>rxyz

   ha%radii_cf=>radii_cf
   ha%nlpspd=>nlpspd !!
   ha%proj=>proj !!
   ha%Lzd=>Lzd !!!
   ha%ngatherarr=>dpcom%ngatherarr
   ha%ndimpot=dpcom%ndimpot
   ha%potential=>pot
   ha%energs = energs
   ha%nspin=nspin
   ha%GPU=>GPU !!
   ha%Gabs_coeffs=>Gabs_coeffs
   ha%PAWD=> PAWD 
   ha%PPD=> PPD
   ha%SIC=>in%SIC
   ha%orbs=>orbs


   call EP_inizializza(ha) 


   if(.true.) then
      LB_nsteps =in%nsteps
      call LB_allocate_for_lanczos( )
      call EP_allocate_for_eigenprob(10)
      call EP_make_dummy_vectors(10)


      do i=0,0

         ene = 0.22_gp + i*0.03_gp

         if( iand( in%potshortcut,16)>0) then
            potential=potentialclone
            do j=1, dpcom%ndimpot

               if( mod(j-1,100)==0) then
                  print *, " dirac_hara punto",j
               endif
               call dirac_hara (rhoXanes(j,1), ene , potential(j,1))
            enddo
         endif

         gamma = 0.03_gp

         if(i==0) then
            useold=.false.
         else
            useold=.true.
         endif

         res =  LB_cg(    get_EP_dim, EP_initialize_start , EP_normalizza,&
            &   EP_Moltiplica4spectra,  EP_copy,  &
            &   EP_scalare,EP_add_from_vect_with_fact   , EP_multbyfact  ,EP_precondition, Ene, gamma, 1.0D-2, useold )

         print *, ene, res
         open(unit=22,file="cgspectra.dat", position="append")
         write(22,*) ene, res
         close(unit=22)

      enddo

      call LB_de_allocate_for_lanczos( )

   endif

   call deallocate_comms(ha%comms,subname)

   call EP_free(ha%iproc)

   if (GPUconv) then
      call free_gpu(GPU,orbs%norbp)
   else if (OCLconv) then
      call free_gpu_OCL(GPU,orbs,in%nspin)
   end if

   i_all=-product(shape(orbs%eval))*kind(orbs%eval)
   deallocate(orbs%eval,stat=i_stat)
   call memocc(i_stat,i_all,'orbs%eval',subname)

   i_all=-product(shape(Gabs_coeffs))*kind(Gabs_coeffs)
   deallocate(Gabs_coeffs,stat=i_stat)
   call memocc(i_stat,i_all,'Gabs_coeffs',subname)

   if( iand( in%potshortcut,16)>0) then
      i_all=-product(shape(potentialclone))*kind(potentialclone)
      deallocate(potentialclone,stat=i_stat)
      call memocc(i_stat,i_all,'potentialclone',subname)
   endif

   call free_full_potential(dpcom%mpi_env%nproc,0,pot,subname)
   nullify(ha%potential)

END SUBROUTINE xabs_cg


subroutine dirac_hara (rho, E , V)
   use module_base

   implicit none

   !Arguments
   real(gp), intent(in   ) :: rho, E
   real(gp), intent(inout) :: V
   !Local variables
   real(gp), parameter :: f= 1.919158d0  !( 4/(9*pi)**(1/3)  )
   real(gp), parameter :: pi=3.141592653589793_gp

   real(gp) :: Vcorr, rs, xk, EV,x
   integer :: i

   if(rho>1.0e-4) then
      rs = (3.0_gp / (4.0_gp*pi*rho)) ** (1.0_gp/3.0_gp)
   else
      rs=1000.0_gp
   endif

   Vcorr=V

   EV=E-Vcorr
   if(EV<=0) then
      return
   endif

   do i=1,10

      EV=E-Vcorr

      if(EV<=0) then
         return
      endif

      xk =sqrt(2*EV  )

      x = xk *rs/ f


      if ((x-1)  < 1.0D-6) return
      Vcorr =V - (f/pi/rs) * ( log(abs((1+x) / (1-x))) * (1-x**2) / (2*x))

   end do
   V=Vcorr
   return
END SUBROUTINE dirac_hara


function GetBottom(atoms,nspin)

   use module_base
   use module_types
   use module_interfaces

   implicit none
   !Arguments
   real(gp) :: GetBottom
   type(atoms_data), intent(in) :: atoms
   integer, intent(in) :: nspin
   !Local variables
   character(len=*), parameter :: subname='GetBottom'
   integer, parameter :: noccmax=2,lmax=4,nmax=6,nelecmax=32, ng=21

   integer :: ity,  i_all
   real(gp) , pointer :: expo(:),  occup(:,:)
   real(gp)   psi(ng,5)

   integer :: i_stat
   real(gp) :: gaenes_aux(5)
   integer, dimension(lmax) :: nl
   integer nspinor, iat, noncoll

   ! if (in_iat_absorber.ne.0) then

   allocate(expo(ng +ndebug  ), stat=i_stat)
   call memocc(i_stat,expo,'expo',subname)

   allocate(occup ( noccmax  ,lmax+1+ndebug ), stat=i_stat)
   call memocc(i_stat,occup,'occup',subname)

   GetBottom=1.0e4_gp

   !for the moment, only collinear
   nspinor=1
   !if non-collinear it is like nspin=1 but with the double of orbitals
   if (nspinor == 4) then
      noncoll=2
   else
      noncoll=1
   end if

   do ity=1, atoms%astruct%ntypes
      do iat=1, atoms%astruct%nat
         if (ity.eq.atoms%astruct%iatype(iat)) exit
      end do
      call count_atomic_shells(lmax,noccmax,nelecmax,nspin,nspinor,atoms%aocc(1,iat),occup,nl)

      call iguess_generator_modified(atoms%nzatom(ity),atoms%nelpsp(ity),&
         &   real(atoms%nelpsp(ity),gp),atoms%psppar(0,0,ity),&
         &   atoms%npspcode(ity),  &
         &   atoms%nlcc_ngv,atoms%nlcc_ngc,atoms%nlccpar,&
         &   ng-1,nl,5,noccmax,lmax,occup,expo,&
         &   psi,.false., gaenes_aux  )

      if( minval(gaenes_aux ) < GetBottom) GetBottom=minval(gaenes_aux )
   enddo

   i_all=-product(shape(occup))*kind(occup)
   deallocate(occup,stat=i_stat)
   call memocc(i_stat,i_all,'occup',subname)

   i_all=-product(shape(expo))*kind(expo)
   deallocate(expo,stat=i_stat)
   call memocc(i_stat,i_all,'expo',subname)

END FUNCTION GetBottom
