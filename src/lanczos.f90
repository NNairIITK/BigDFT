
subroutine lanczos(iproc,nproc,at,hx,hy,hz,rxyz,&
     cpmult,fpmult,radii_cf,nlpspd,proj,lr,ngatherarr,ndimpot,potential,&
      ekin_sum,epot_sum,eproj_sum,nspin,GPU,in_iat_absorber, doorthoocc, Occ_norb, Occ_psit, Occ_eval,&
  in  )! aggiunger a interface
  use module_base
  use module_types
  use lanczos_interface
  use lanczos_base
  ! per togliere il bug 
  use module_interfaces ,except_this_one => lanczos


  implicit none
  integer  :: iproc,nproc,ndimpot,nspin
  real(gp)  :: hx,hy,hz,cpmult,fpmult
  type(atoms_data), target :: at
  type(nonlocal_psp_descriptors), target :: nlpspd
  type(locreg_descriptors), target :: lr
  integer, dimension(0:nproc-1,2), target :: ngatherarr 
  real(gp), dimension(3,at%nat), target :: rxyz
  real(gp), dimension(at%ntypes,3), target :: radii_cf  
  real(wp), dimension(nlpspd%nprojel), target :: proj
  real(wp), dimension(max(ndimpot,1),nspin), target :: potential

  real(gp) :: ekin_sum,epot_sum,eproj_sum
  type(GPU_pointers), intent(inout) , target :: GPU
  integer, intent(in) :: in_iat_absorber
  
  logical, intent(in) :: doorthoocc
  integer, intent(in)  :: Occ_norb
  real(wp), dimension(:), pointer :: Occ_psit
  real(wp), dimension(:), pointer :: Occ_eval

  type(input_variables),intent(in) :: in


  !local variables
  character(len=*), parameter :: subname='lanczos'
  integer :: i_stat,i_all
  type(lanczos_args) :: ha
  real(8) shift, tol
  real(8) dum
  integer i, cercacount, ierr

  type(gaussian_basis), target ::  Gabsorber
  real(wp),   pointer  :: Gabs_coeffs(:)
  real(wp),  pointer :: dum_coeffs(:,:)
  character(len=80) :: filename
  logical :: projeexists


  !print *, " IN ROUTINE LANCZOS "

  !create the orbitals descriptors, for virtual and inputguess orbitals
  call orbitals_descriptors(iproc,nproc,1,1,0,1,in%nkpt,in%kpt,in%wkpt,ha%orbs)

  if (GPUconv) then
     call prepare_gpu_for_locham(lr%d%n1,lr%d%n2,lr%d%n3,in%nspin,&
          hx,hy,hz,lr%wfd,ha%orbs,GPU)
  end if
  GPU%full_locham=.true.


  allocate(ha%orbs%eval(ha%orbs%norb+ndebug),stat=i_stat)
  call memocc(i_stat,ha%orbs%eval,'ha%orbs%eval',subname)
  ha%orbs%occup(1:ha%orbs%norb)=1.0_gp
  ha%orbs%spinsgn(1:ha%orbs%norb)=1.0_gp
  ha%orbs%eval(1:ha%orbs%norb)=1.0_gp
  !call allocate_comms(nproc,ha%comms,subname)
  call orbitals_communicators(iproc,nproc,lr,ha%orbs,ha%comms)  
  allocate(Gabs_coeffs(2*in%L_absorber+1+ndebug),stat=i_stat)
  call memocc(i_stat,Gabs_coeffs,'Gabs_coeffs',subname)
 
!!$  
  !call allocate_comms(nproc,ha%orbs,ha%comms,subname)
     ! write(filename,'(A,A,A)') "gproje_", at%atomnames(at%iatype(  in_iat_absorber )) , "_1s_dipole"
     write(filename,'(A,A,A,I1)') "gproje_", at%atomnames(at%iatype(  in_iat_absorber )) , "_1s_",  in%L_absorber

     inquire(FILE=filename,EXIST=projeexists)
     if( projeexists .and. .not. in%abscalc_eqdiff   ) then


        print *, "leggo " ,  in%abscalc_eqdiff  
        nullify( dum_coeffs  ) 

        call read_gaussian_information(0 ,1 ,ha%orbs,Gabsorber,dum_coeffs , filename, .true. )

        Gabsorber%rxyz(:,1)=rxyz(:, in_iat_absorber )

        Gabs_coeffs(:)=in%Gabs_coeffs(:)

     else
          
           
           call GetExcitedOrbitalAsG(in_iat_absorber ,Gabsorber,&
                at,rxyz,nproc,iproc,1,    in%L_absorber ,  in%abscalc_eqdiff )
           print * , " uscito get "
           
           allocate(dum_coeffs(2*in%L_absorber+1,1+ndebug),stat=i_stat)
           call memocc(i_stat,dum_coeffs,'dum_coeffs',subname)
           
           dum_coeffs(:,1) = in%Gabs_coeffs(:)
           
           Gabs_coeffs(:)=in%Gabs_coeffs(:)
           
        if( iproc.eq.0) then
           print *, " chiamo write " 
           call write_gaussian_information( 0 ,1 ,ha%orbs,Gabsorber,dum_coeffs ,filename)
        endif
     endif
     
  

  !associate hamapp_arg pointers
  ha%iproc=iproc
  ha%nproc=nproc
  ha%at=>at !!
  ha%hx=hx 
  ha%hy=hy
  ha%hz=hz
  ha%rxyz=>rxyz
  ha%cpmult=cpmult
  ha%fpmult=fpmult


  ha%radii_cf=>radii_cf
  ha%nlpspd=>nlpspd !!
  ha%proj=>proj !!
  ha%lr=>lr !!!
  ha%ngatherarr=>ngatherarr
  ha%ndimpot=ndimpot
  ha%potential=>potential
  ha%ekin_sum=ekin_sum
  ha%epot_sum=epot_sum
  ha%eproj_sum=eproj_sum
  ha%nspin=nspin
  ha%GPU=>GPU !!
  ha%Gabsorber=>Gabsorber 
  ha%Gabs_coeffs=>Gabs_coeffs
  
  !initialise the arguments for HamiltonianApplication


  call EP_inizializza(ha) 
  

  call EP_memorizza_stato(Gabsorber) ! se uno stato e' memorizzato EP_initialize_start usa quello, se no random
     


  if(.true.) then
     LB_nsteps =2000
     call LB_allocate_for_lanczos( )
     call EP_allocate_for_eigenprob(LB_nsteps)
     call EP_make_dummy_vectors(10)

     
     
     if(doorthoocc) then
        print *, "  in lanczos associated(Occ_psit) ", associated(Occ_psit) 
        call EP_store_occupied_orbitals(iproc, nproc, Occ_norb, Occ_psit ) 
     endif
     
     call LB_passeggia(0,LB_nsteps,      get_EP_dim, EP_initialize_start , EP_normalizza,&
          EP_Moltiplica, EP_GramSchmidt ,EP_set_all_random, EP_copy,   EP_mat_mult, &
          EP_scalare,EP_add_from_vect_with_fact     )


     

     open(unit=22,file="alphabeta")
     write(22,*) LB_nsteps, EP_norma2_initialized_state
     do i=0, LB_nsteps-1
        write(22,*) LB_alpha(i), LB_beta(i)
     enddo
     if(doorthoocc) then
        write(22,*) Occ_norb
        do i=1, Occ_norb
           write(22,*) Occ_eval(i), EP_occprojections
        enddo
     endif
     close(unit=22)
  endif


  if(0.eq.1) then
     print *, " chiamo cerca  " , get_EP_dim()
     
     shift =-300.0
     tol   =1.0D-14


     LB_iproc=iproc
     LB_nproc=nproc
     
     cercacount = LB_cerca( 10, shift, tol, set_EP_shift, EP_allocate_for_eigenprob, EP_make_dummy_vectors, &
          get_EP_dim, EP_initialize_start , EP_normalizza, EP_Moltiplica, EP_GramSchmidt, &
          EP_set_all_random, EP_copy , EP_mat_mult,  EP_scalare,EP_add_from_vect_with_fact)
     
     
     stop


     ! ---------------------------------------------g

     call EP_allocate_for_eigenprob(10)
     
     print *, " chiamo EP_make_dummy_vectors "
     call  EP_make_dummy_vectors(10)    

     if ( iproc.eq.0)    print *, " chiamo initialize start " 
     call EP_initialize_start()
     call EP_normalizza(0)
     ! call EP_initialize_start_0( Gabsorber)
     
!!$     if ( iproc.eq.0)    print *, " copio psi "
!!$     call EP_copia_per_prova(psi)
     call  EP_Moltiplica(-1,0) 
     dum= EP_scalare(0,-1)
     ! if (iproc.eq.0) then
     print *, " il prodotto scalare psi H psi est " , dum 
     print * , " uscito da EP_initialize_start "
     ! endif
     stop
     
     !------------------------------------------
     
    
     print *, " ciao " 
     print *, LB_eval(:9)
     call LB_de_allocate_for_lanczos( )
  endif



  if(0.eq.1) then
     print *, " chiamo allocate for eigen prob"
     call EP_allocate_for_eigenprob(10)
     
     print *, " chiamo EP_make_dummy_vectors "
     call  EP_make_dummy_vectors(10)    

     if ( iproc.eq.0)    print *, " chiamo initialize start " 
     call EP_initialize_start()
     call EP_normalizza(0)
     ! call EP_initialize_start_0( Gabsorber)
    
!!$     if ( iproc.eq.0)    print *, " copio psi "
!!$     call EP_copia_per_prova(psi)
     call  EP_Moltiplica(-1,0) 
     dum= EP_scalare(0,-1)
     ! if (iproc.eq.0) then
        print *, " il prodotto scalare psi H psi est " , dum 
        print * , " uscito da EP_initialize_start "
     ! endif
  endif

  !deallocate communication and orbitals descriptors
  call deallocate_comms(ha%comms,subname)


  call EP_free()

  if (GPUconv) then
     call free_gpu(GPU,ha%orbs%norbp)
  end if


  call deallocate_orbs(ha%orbs,subname)

  i_all=-product(shape(ha%orbs%eval))*kind(ha%orbs%eval)
  deallocate(ha%orbs%eval,stat=i_stat)
  call memocc(i_stat,i_all,'ha%orbs%spinsgn',subname)

!!$  i_all=-product(shape(Gabs_coeffs))*kind(Gabs_coeffs)
!!$  deallocate(Gabs_coeffs,stat=i_stat)
!!$  call memocc(i_stat,i_all,'Gabs_coeffs',subname)



end subroutine lanczos

subroutine chebychev(iproc,nproc,at,hx,hy,hz,rxyz,&
     cpmult,fpmult,radii_cf,nlpspd,proj,lr,ngatherarr,ndimpot,potential,&
     ekin_sum,epot_sum,eproj_sum,nspin,GPU,in_iat_absorber,in  )! aggiunger a interface
  use module_base
  use module_types
  use lanczos_interface
  use lanczos_base
  ! per togliere il bug 
  use module_interfaces, except_this_one => chebychev
  

  implicit none
  integer  :: iproc,nproc,ndimpot,nspin
  real(gp)  :: hx,hy,hz,cpmult,fpmult
  type(atoms_data), target :: at
  type(nonlocal_psp_descriptors), target :: nlpspd
  type(locreg_descriptors), target :: lr
  integer, dimension(0:nproc-1,2), target :: ngatherarr 
  real(gp), dimension(3,at%nat), target :: rxyz
  real(gp), dimension(at%ntypes,3), target :: radii_cf  
  real(wp), dimension(nlpspd%nprojel), target :: proj
  real(wp), dimension(max(ndimpot,1),nspin), target :: potential

  real(gp) :: ekin_sum,epot_sum,eproj_sum
  type(GPU_pointers), intent(inout) , target :: GPU
  integer, intent(in) :: in_iat_absorber
  

  type(input_variables),intent(in) :: in


  !local variables
  character(len=*), parameter :: subname='chebychev'
  integer :: i_stat,i_all
  type(lanczos_args) :: ha
  real(8) shift, tol
  real(8) dum
  integer i, cercacount, ierr

  type(gaussian_basis), target ::  Gabsorber
  real(wp),   pointer  :: Gabs_coeffs(:)
  real(wp),  pointer :: dum_coeffs(:,:)
  character(len=80) :: filename
  logical :: projeexists
  real(gp) eval_min, eval_max, fact_cheb, cheb_shift
  integer accontentati_di

  print *, " IN ROUTINE  chebychev  "

  !create the orbitals descriptors, for virtual and inputguess orbitals
  call orbitals_descriptors(iproc,nproc,1,1,0,1,in%nkpt,in%kpt,in%wkpt,ha%orbs)

  if (GPUconv) then
     call prepare_gpu_for_locham(lr%d%n1,lr%d%n2,lr%d%n3,in%nspin,&
          hx,hy,hz,lr%wfd,ha%orbs,GPU)
  end if
  GPU%full_locham=.true.


  allocate(ha%orbs%eval(ha%orbs%norb+ndebug),stat=i_stat)
  call memocc(i_stat,ha%orbs%eval,'ha%orbs%eval',subname)
  ha%orbs%occup(1:ha%orbs%norb)=1.0_gp
  ha%orbs%spinsgn(1:ha%orbs%norb)=1.0_gp
  ha%orbs%eval(1:ha%orbs%norb)=1.0_gp
  !call allocate_comms(nproc,ha%comms,subname)
  call orbitals_communicators(iproc,nproc,lr,ha%orbs,ha%comms)  
  allocate(Gabs_coeffs(2*in%L_absorber+1+ndebug),stat=i_stat)
  call memocc(i_stat,Gabs_coeffs,'Gabs_coeffs',subname)
 
!!$  
  !call allocate_comms(nproc,ha%orbs,ha%comms,subname)
  write(filename,'(A,A,A,I1)') "gproje_", at%atomnames(at%iatype(  in_iat_absorber )) , "_1s_",  in%L_absorber
     
  inquire(FILE=filename,EXIST=projeexists)
  if( projeexists .and. .not. in%abscalc_eqdiff   ) then


     print *, "leggo " ,  in%abscalc_eqdiff  
     nullify( dum_coeffs  ) 
     
     call read_gaussian_information(0 ,1 ,ha%orbs,Gabsorber,dum_coeffs , filename, .true. )
     
     Gabsorber%rxyz(:,1)=rxyz(:, in_iat_absorber )
     
     Gabs_coeffs(:)=in%Gabs_coeffs(:)
     
  else
     call GetExcitedOrbitalAsG(in_iat_absorber ,Gabsorber,&
          at,rxyz,nproc,iproc,1,    in%L_absorber ,  in%abscalc_eqdiff )
     print * , " uscito get "
           
           allocate(dum_coeffs(2*in%L_absorber+1,1+ndebug),stat=i_stat)
           call memocc(i_stat,dum_coeffs,'dum_coeffs',subname)
           
           dum_coeffs(:,1) = in%Gabs_coeffs(:)
           
           Gabs_coeffs(:)=in%Gabs_coeffs(:)
           
        if( iproc.eq.0) then
           print *, " chiamo write " 
           call write_gaussian_information( 0 ,1 ,ha%orbs,Gabsorber,dum_coeffs ,filename)
        endif
     endif
     


  

  !associate hamapp_arg pointers
  ha%iproc=iproc
  ha%nproc=nproc
  ha%at=>at !!
  ha%hx=hx 
  ha%hy=hy
  ha%hz=hz
  ha%rxyz=>rxyz
  ha%cpmult=cpmult
  ha%fpmult=fpmult


  ha%radii_cf=>radii_cf
  ha%nlpspd=>nlpspd !!
  ha%proj=>proj !!
  ha%lr=>lr !!!
  ha%ngatherarr=>ngatherarr
  ha%ndimpot=ndimpot
  ha%potential=>potential
  ha%ekin_sum=ekin_sum
  ha%epot_sum=epot_sum
  ha%eproj_sum=eproj_sum
  ha%nspin=nspin
  ha%GPU=>GPU !!
  ha%Gabsorber=>Gabsorber 
  ha%Gabs_coeffs=>Gabs_coeffs
  
  !initialise the arguments for HamiltonianApplication


  call EP_inizializza(ha) 
  
  call  EP_memorizza_stato(Gabsorber) 


  if(.true.) then

     ! trova il valore massimo 
     shift =-0.0
     tol   =1.0D-8
     accontentati_di=1
     
     cercacount = LB_cerca( 10, shift, tol, set_EP_shift, EP_allocate_for_eigenprob, EP_make_dummy_vectors, &
          get_EP_dim, EP_initialize_start , EP_normalizza, EP_Moltiplica, EP_GramSchmidt, &
          EP_set_all_random, EP_copy , EP_mat_mult,  EP_scalare,EP_add_from_vect_with_fact,accontentati_di)
     
     if(iproc==0) then
        print *, " valori propri massimi " 
        print *, LB_eval
     endif
     eval_max = LB_eval(0)
     
     
     ! trova il valore minimo 
     shift =-10000
     
     
     cercacount = LB_cerca( 10, shift, tol, set_EP_shift, EP_allocate_for_eigenprob, EP_make_dummy_vectors, &
          get_EP_dim, EP_initialize_start , EP_normalizza, EP_Moltiplica, EP_GramSchmidt, &
          EP_set_all_random, EP_copy , EP_mat_mult,  EP_scalare,EP_add_from_vect_with_fact,accontentati_di)
     
     if(iproc==0) then
        print *, " valori propri minimi " 
        print *, LB_eval
     endif
     eval_min = LB_eval(0)+10000

  else
     eval_min=-3.0
     eval_max=380.0
  endif

  cheb_shift=0.5*(eval_min+ eval_max) 
  fact_cheb = (2-0.0001)/(eval_max-eval_min)
     
        
     
  if(.true.) then
     
     call EP_memorizza_stato(Gabsorber) ! se uno stato e' memorizzato EP_initialize_start usa quello, se no random
     LB_nsteps =2000
     
     call LB_allocate_for_chebychev( )
     call EP_allocate_for_eigenprob(3) ! invece di nsteps, giusto qualche vettore per fare i calcoli
     call EP_make_dummy_vectors(2)
  

     call set_EP_shift(-cheb_shift) 
     
     call LB_passeggia_Chebychev (LB_nsteps, cheb_shift,  fact_cheb,     get_EP_dim, EP_initialize_start , EP_normalizza,&
          EP_Moltiplica, EP_GramSchmidt ,EP_set_all_random, EP_copy,   EP_mat_mult, &
          EP_scalare,EP_add_from_vect_with_fact  , EP_multbyfact  )

     open(unit=22,file="alphabeta")
     write(22,*) LB_nsteps, EP_norma2_initialized_state
     do i=0, LB_nsteps-1
        write(22,*) LB_alpha(i)
     enddo
     close(unit=22)
  endif

  !deallocate communication and orbitals descriptors
  call deallocate_comms(ha%comms,subname)

  call EP_free()

  if (GPUconv) then
     call free_gpu(GPU,ha%orbs%norbp)
  end if


  call deallocate_orbs(ha%orbs,subname)

  i_all=-product(shape(ha%orbs%eval))*kind(ha%orbs%eval)
  deallocate(ha%orbs%eval,stat=i_stat)
  call memocc(i_stat,i_all,'ha%orbs%spinsgn',subname)

!!$  i_all=-product(shape(Gabs_coeffs))*kind(Gabs_coeffs)
!!$  deallocate(Gabs_coeffs,stat=i_stat)
!!$  call memocc(i_stat,i_all,'Gabs_coeffs',subname)



end subroutine chebychev

