!> @file
!!  Lanczos diagonalisation used by XANES calculation
!! @author
!!    Copyright (C) 2009-2011 BigDFT group (AM, LG)
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 

!>   Lanczos diagonalization
subroutine xabs_lanczos(iproc,nproc,at,hx,hy,hz,rxyz,&
     radii_cf,nlpspd,proj,lr,ngatherarr,ndimpot,potential,&
     ekin_sum,epot_sum,eproj_sum,nspin,GPU,in_iat_absorber,&
     in  )! aggiunger a interface
  use module_base
  use module_types
  use lanczos_interface
  use lanczos_base
  ! per togliere il bug 
  use module_interfaces ,except_this_one => xabs_lanczos

  implicit none

  integer, intent(in) :: iproc,nproc,ndimpot,nspin
  real(gp), intent(in) :: hx,hy,hz
  type(atoms_data), target :: at
  type(nonlocal_psp_descriptors), target :: nlpspd
  type(locreg_descriptors), target :: lr
  integer, dimension(0:nproc-1,2), target :: ngatherarr 
  real(gp), dimension(3,at%nat), target :: rxyz
  real(gp), dimension(at%ntypes,3), intent(in), target ::  radii_cf
  real(wp), dimension(nlpspd%nprojel), target :: proj
  real(wp), dimension(max(ndimpot,1),nspin), target :: potential

  real(gp) :: ekin_sum,epot_sum,eproj_sum
  type(GPU_pointers), intent(inout) , target :: GPU
  integer, intent(in) :: in_iat_absorber
  

  type(input_variables),intent(in) :: in

  !local variables
  character(len=*), parameter :: subname='lanczos'
  integer :: i_stat,i_all
  type(lanczos_args) :: ha
  integer :: i

  type(gaussian_basis), target ::  Gabsorber
  real(wp),   pointer  :: Gabs_coeffs(:)
  real(wp), dimension(:), pointer  :: pot
  real(wp),  pointer, dimension(:,:)  :: dum_coeffs
  character(len=800) :: filename
  logical :: projeexists

  if(iproc==0) print *, " IN ROUTINE LANCZOS "

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
 
  write(filename,'(A,A,A,I1)') "gproje_", trim(at%atomnames(at%iatype(  in_iat_absorber ))) , "_1s_",  in%L_absorber
  
  inquire(FILE=trim(filename),EXIST=projeexists)
  
  if( projeexists .and. .not. in%abscalc_eqdiff   ) then
     
     if(iproc==0) then
        print *, "reading  precalculated  projection on pseudofunctions"
        print *, "for 1s of atom number ", in_iat_absorber, " ( atomname = ", at%atomnames(at%iatype(  in_iat_absorber ))," )"
        print *, "After application of a 2*L-pole  with L= ", in%L_absorber
        print *," from file " , filename
     endif

     nullify( dum_coeffs  ) 
     nullify(Gabsorber%nshell)
     nullify(Gabsorber%nam)
     nullify(Gabsorber%ndoc)
     nullify(Gabsorber%xp)
     nullify(Gabsorber%psiat)
     
     call read_gaussian_information (ha%orbs,Gabsorber,dum_coeffs , filename, .true. )    
     Gabsorber%rxyz(:,1)=rxyz(:, in_iat_absorber )
     
     i_all=-product(shape(dum_coeffs))*kind(dum_coeffs)
     deallocate(dum_coeffs,stat=i_stat)
     call memocc(i_stat,i_all,'coeffs',subname)
     
     Gabs_coeffs(:)=in%Gabs_coeffs(:)
     
  else
     
        if(iproc==0) then
           print *, "calculating  projection on pseudofunctions"
           print *, "for 1s of atom number ", in_iat_absorber, " ( atomname = ", at%atomnames(at%iatype(  in_iat_absorber ))," )"
           print *, "After application of a 2*L-pole  with L= ", in%L_absorber
        endif
           call GetExcitedOrbitalAsG(in_iat_absorber ,Gabsorber,&
                at,rxyz,nproc,iproc,   in%L_absorber ,  in%abscalc_eqdiff )
           
           allocate(dum_coeffs(2*in%L_absorber+1,1+ndebug),stat=i_stat)
           call memocc(i_stat,dum_coeffs,'dum_coeffs',subname)
           
           dum_coeffs(:,1) = in%Gabs_coeffs(:)
           
           Gabs_coeffs(:)=in%Gabs_coeffs(:)
           
        if( iproc.eq.0) then
           print *,"writing them on file " , filename
           call write_gaussian_information( 0 ,1 ,ha%orbs,Gabsorber,dum_coeffs ,filename)
        endif
     
        i_all=-product(shape(dum_coeffs))*kind(dum_coeffs)
        deallocate(dum_coeffs,stat=i_stat)
        call memocc(i_stat,i_all,'coeffs',subname)
        
     endif

  !allocate the potential in the full box
  call full_local_potential(iproc,nproc,ndimpot,lr%d%n1i*lr%d%n2i*lr%d%n3i,in%nspin,&
       ha%orbs%norb,ha%orbs%norbp,ngatherarr,potential,pot)


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
  ha%lr=>lr !!!
  ha%ngatherarr=>ngatherarr
  ha%ndimpot=ndimpot
  ha%potential=>pot
  ha%ekin_sum=ekin_sum
  ha%epot_sum=epot_sum
  ha%eproj_sum=eproj_sum
  ha%nspin=nspin
  ha%GPU=>GPU !!
  ha%Gabsorber=>Gabsorber 
  ha%Gabs_coeffs=>Gabs_coeffs

  call EP_inizializza(ha) 

  call EP_memorizza_stato(Gabsorber) 
     
  if(.true.) then
     LB_nsteps =in%nsteps
     call LB_allocate_for_lanczos( )
     call EP_allocate_for_eigenprob(LB_nsteps)
     call EP_make_dummy_vectors(10)
     
     
     call LB_passeggia(0,LB_nsteps,      get_EP_dim, EP_initialize_start , EP_normalizza,&
          EP_Moltiplica, EP_GramSchmidt ,EP_set_all_random, EP_copy,   EP_mat_mult, &
          EP_scalare,EP_add_from_vect_with_fact     )


     if(ha%iproc==0) then
        open(unit=22,file="alphabeta")
        write(22,*) LB_nsteps, EP_norma2_initialized_state
        print *, " alpha and beta from lanczos "
        WRITE(*,'(I5,1ES23.16)')    LB_nsteps, EP_norma2_initialized_state
        do i=0, LB_nsteps-1
           write(22,*) LB_alpha(i), LB_beta(i)
           WRITE(*,'(2ES23.16)')  LB_alpha(i), LB_beta(i)
        enddo

        close(unit=22)
     endif

     call LB_de_allocate_for_lanczos( )


  endif


  
  call deallocate_comms(ha%comms,subname)

  call EP_free()

  if (GPUconv) then
     call free_gpu(GPU,ha%orbs%norbp)
  end if

  call deallocate_orbs(ha%orbs,subname)

  i_all=-product(shape(Gabsorber%rxyz))*kind(Gabsorber%rxyz)
  deallocate(Gabsorber%rxyz,stat=i_stat)
  call memocc(i_stat,i_all,'Gabsorber%rxyz',subname)

  call deallocate_gwf(Gabsorber, subname)

  i_all=-product(shape(ha%orbs%eval))*kind(ha%orbs%eval)
  deallocate(ha%orbs%eval,stat=i_stat)
  call memocc(i_stat,i_all,'ha%orbs%eval',subname)

  i_all=-product(shape(Gabs_coeffs))*kind(Gabs_coeffs)
  deallocate(Gabs_coeffs,stat=i_stat)
  call memocc(i_stat,i_all,'Gabs_coeffs',subname)


  call deallocate_abscalc_input(in, subname)


END SUBROUTINE xabs_lanczos


!>   Chebychev polynomials to calculate the density of states
subroutine xabs_chebychev(iproc,nproc,at,hx,hy,hz,rxyz,&
     radii_cf,nlpspd,proj,lr,ngatherarr,ndimpot,potential,&
     ekin_sum,epot_sum,eproj_sum,nspin,GPU,in_iat_absorber,in  )! aggiunger a interface

  use module_base
  use module_types
  use lanczos_interface
  use lanczos_base
  ! per togliere il bug 
  use module_interfaces, except_this_one => xabs_chebychev
  
  implicit none
  integer  :: iproc,nproc,ndimpot,nspin
  real(gp)  :: hx,hy,hz
  type(atoms_data), target :: at
  type(nonlocal_psp_descriptors), target :: nlpspd
  type(locreg_descriptors), target :: lr
  integer, dimension(0:nproc-1,2), target :: ngatherarr 
  real(gp), dimension(3,at%nat), target :: rxyz
  real(gp), dimension(at%ntypes,3), intent(in), target ::  radii_cf
  real(wp), dimension(nlpspd%nprojel), target :: proj
  real(wp), dimension(max(ndimpot,1),nspin), target :: potential

  real(gp) :: ekin_sum,epot_sum,eproj_sum
  type(GPU_pointers), intent(inout) , target :: GPU
  integer, intent(in) :: in_iat_absorber
  
  type(input_variables),intent(in) :: in

  !Local variables
  character(len=*), parameter :: subname='chebychev'
  integer :: i_stat,i_all
  type(lanczos_args) :: ha
  real(8) :: shift, tol
  integer :: i, cercacount

  type(gaussian_basis), target ::  Gabsorber
  real(wp), pointer :: Gabs_coeffs(:)
  real(wp), dimension(:), pointer :: pot
  real(wp), pointer, dimension (:,:) :: dum_coeffs
  character(len=80) :: filename
  logical :: projeexists
  real(gp) :: eval_min, eval_max, fact_cheb, cheb_shift
  integer :: accontentati_di
  real(gp) :: Pi

  if (iproc==0) print *, " IN ROUTINE  chebychev  "

  Pi=acos(-1.0_gp)

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

  call orbitals_communicators(iproc,nproc,lr,ha%orbs,ha%comms)  

  allocate(Gabs_coeffs(2*in%L_absorber+1+ndebug),stat=i_stat)
  call memocc(i_stat,Gabs_coeffs,'Gabs_coeffs',subname)
 

  write(filename,'(A,A,A,I1)') "gproje_", trim(at%atomnames(at%iatype(  in_iat_absorber ))) , "_1s_",  in%L_absorber

  inquire(FILE=trim(filename),EXIST=projeexists)

  if( projeexists .and. .not. in%abscalc_eqdiff   ) then
     
     if(iproc==0) then
        print *, "reading  precalculated  projection on pseudofunctions"
        print *, "for 1s of atom number ", in_iat_absorber, " ( atomname = ", at%atomnames(at%iatype(  in_iat_absorber ))," )"
        print *, "After application of a 2*L-pole  with L= ", in%L_absorber
        print *," from file " , filename
     endif
     
     nullify( dum_coeffs  ) 

     nullify(Gabsorber%nshell)
     nullify(Gabsorber%nam)
     nullify(Gabsorber%ndoc)
     nullify(Gabsorber%xp)
     nullify(Gabsorber%psiat)
   
     call read_gaussian_information(ha%orbs,Gabsorber,dum_coeffs , filename, .true. )
     Gabsorber%rxyz(:,1)=rxyz(:, in_iat_absorber )
     
     i_all=-product(shape(dum_coeffs))*kind(dum_coeffs)
     deallocate(dum_coeffs,stat=i_stat)
     call memocc(i_stat,i_all,'coeffs',subname)

     Gabs_coeffs(:)=in%Gabs_coeffs(:)
  else
     if(iproc==0) then
        print *, "calculating  projection on pseudofunctions"
        print *, "for 1s of atom number ", in_iat_absorber, " ( atomname = ", at%atomnames(at%iatype(  in_iat_absorber ))," )"
        print *, "After application of a 2*L-pole  with L= ", in%L_absorber
     endif
     
     call GetExcitedOrbitalAsG(in_iat_absorber ,Gabsorber,&
          at,rxyz,nproc,iproc,   in%L_absorber ,  in%abscalc_eqdiff )
     
     allocate(dum_coeffs(2*in%L_absorber+1,1+ndebug),stat=i_stat)
     call memocc(i_stat,dum_coeffs,'dum_coeffs',subname)
     
     dum_coeffs(:,1) = in%Gabs_coeffs(:)
     
     Gabs_coeffs(:)=in%Gabs_coeffs(:)
     
     if( iproc.eq.0) then
        print *,"writing them on file " , filename
        call write_gaussian_information( 0 ,1 ,ha%orbs,Gabsorber,dum_coeffs ,trim(filename))
     endif
   
     
     i_all=-product(shape(dum_coeffs))*kind(dum_coeffs)
     deallocate(dum_coeffs,stat=i_stat)
     call memocc(i_stat,i_all,'coeffs',subname)
     
  endif

  !allocate the potential in the full box
  call full_local_potential(iproc,nproc,ndimpot,lr%d%n1i*lr%d%n2i*lr%d%n3i,in%nspin,&
       ha%orbs%norb,ha%orbs%norbp,ngatherarr,potential,pot)
  
  print *, "OK "
  !associate hamapp_arg pointers
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
  ha%lr=>lr !!!
  ha%ngatherarr=>ngatherarr
  ha%ndimpot=ndimpot
  ha%potential=>pot
  ha%ekin_sum=ekin_sum
  ha%epot_sum=epot_sum
  ha%eproj_sum=eproj_sum
  ha%nspin=nspin
  ha%GPU=>GPU !!
  ha%Gabsorber=>Gabsorber 
  ha%Gabs_coeffs=>Gabs_coeffs
 
  call EP_inizializza(ha)  
  print *, "OK 1"
 
  call  EP_memorizza_stato(Gabsorber) 

  print *, "OK 2"
 
  if(.false.) then

     ! trova il valore massimo 
     shift =-0.0
     tol   =1.0D-8
     accontentati_di=1
     
     cercacount = LB_cerca( 10, shift, tol, set_EP_shift, EP_allocate_for_eigenprob, EP_make_dummy_vectors, &
          get_EP_dim, EP_initialize_start , EP_normalizza, EP_Moltiplica, EP_GramSchmidt, &
          EP_set_all_random, EP_copy , EP_mat_mult,  EP_scalare,EP_add_from_vect_with_fact,accontentati_di)
     
     if(iproc==0) then
        print *, " maximal eigenvalues " 
        print *, LB_eval
     endif
     eval_max = LB_eval(0)
     
     ! trova il valore minimo 
     shift =-10000
     
     
     cercacount = LB_cerca( 10, shift, tol, set_EP_shift, EP_allocate_for_eigenprob, EP_make_dummy_vectors, &
          get_EP_dim, EP_initialize_start , EP_normalizza, EP_Moltiplica, EP_GramSchmidt, &
          EP_set_all_random, EP_copy , EP_mat_mult,  EP_scalare,EP_add_from_vect_with_fact,accontentati_di)
     
     if(iproc==0) then
        print *, " minima eigenvalues" 
        print *, LB_eval
     endif
     eval_min = LB_eval(0)+10000

  else
     eval_min = GetBottom(  at, iproc)-1.0
     eval_max = 4.0*Pi*Pi*(1.0/hx/hx + 1.0/hy/hy + 1.0/hz/hz  )/2.0*1.01+1
  endif



  cheb_shift=0.5*(eval_min+ eval_max) 
  fact_cheb = (2-0.0001)/(eval_max-eval_min)
     
  print *, "OK 4"


  if(.true.) then
     
     call EP_memorizza_stato(Gabsorber) ! se uno stato e' memorizzato EP_initialize_start usa quello, se no random
     LB_nsteps = in%nsteps
     
     print *, "OK 45"
 
     call LB_allocate_for_chebychev( )
     call EP_allocate_for_eigenprob(3) ! invece di nsteps, giusto qualche vettore per fare i calcoli
     call EP_make_dummy_vectors(2)
  

     call set_EP_shift(-cheb_shift) 
     print *, "OK 5"
 
     call LB_passeggia_Chebychev (LB_nsteps, cheb_shift,  fact_cheb,     get_EP_dim, EP_initialize_start , EP_normalizza,&
          EP_Moltiplica, EP_GramSchmidt ,EP_set_all_random, EP_copy,   EP_mat_mult, &
          EP_scalare,EP_add_from_vect_with_fact  , EP_multbyfact  )

     if(ha%iproc==0) then
        print *, "coefficients from Chebychev "
        WRITE(*,'(I5,2ES23.16)')   2*LB_nsteps, cheb_shift,  fact_cheb
        print *,"... " 
        do i=0, 2*LB_nsteps-1
           if(i>2*LB_nsteps-1 -10) then
              WRITE(*,'(1ES23.16)')   LB_alpha(i)
           endif
        enddo
     endif
  endif

  call free_full_potential(nproc,pot,subname)

  nullify(ha%potential)

  !deallocate communication and orbitals descriptors
  call deallocate_comms(ha%comms,subname)

  call EP_free()
  call  LB_de_allocate_for_lanczos( )

!!$ this free is already executed by bigdft
!!$
  if (GPUconv) then
     call free_gpu(GPU,ha%orbs%norbp)
  end if

  i_all=-product(shape(ha%orbs%eval))*kind(ha%orbs%eval)
  deallocate(ha%orbs%eval,stat=i_stat)
  call memocc(i_stat,i_all,'ha%orbs%eval',subname)
  
  i_all=-product(shape(Gabs_coeffs))*kind(Gabs_coeffs)
  deallocate(Gabs_coeffs,stat=i_stat)
  call memocc(i_stat,i_all,'Gabs_coeffs',subname)

  call deallocate_orbs(ha%orbs,subname)

  i_all=-product(shape(Gabsorber%rxyz))*kind(Gabsorber%rxyz)
  deallocate(Gabsorber%rxyz,stat=i_stat)
  call memocc(i_stat,i_all,'Gabsorber%rxyz',subname)

  call deallocate_gwf(Gabsorber, subname)

  call deallocate_abscalc_input(in, subname)


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


!>   finds the spectra solving  (H-omega)x=b
subroutine xabs_cg(iproc,nproc,at,hx,hy,hz,rxyz,&
     radii_cf,nlpspd,proj,lr,ngatherarr,ndimpot,potential,&
     ekin_sum,epot_sum,eproj_sum,nspin,GPU,in_iat_absorber,&
     in , rhoXanes )
  use module_base
  use module_types
  use lanczos_interface
  use lanczos_base
  ! per togliere il bug 
  use module_interfaces ,except_this_one => xabs_lanczos

  implicit none

  integer  :: iproc,nproc,ndimpot,nspin
  real(gp)  :: hx,hy,hz
  type(atoms_data), target :: at
  type(nonlocal_psp_descriptors), target :: nlpspd
  type(locreg_descriptors), target :: lr
  integer, dimension(0:nproc-1,2), target :: ngatherarr 
  real(gp), dimension(3,at%nat), target :: rxyz
  real(gp), dimension(at%ntypes,3), intent(in), target ::  radii_cf
  real(wp), dimension(nlpspd%nprojel), target :: proj
  real(wp), dimension(max(ndimpot,1),nspin), target :: potential
  real(wp), dimension(max(ndimpot,1),nspin), target :: rhoXanes



  real(gp) :: ekin_sum,epot_sum,eproj_sum
  type(GPU_pointers), intent(inout) , target :: GPU
  integer, intent(in) :: in_iat_absorber

  type(input_variables),intent(in) :: in

  !local variables
  character(len=*), parameter :: subname='xabs_cg'
  integer :: i_stat,i_all
  type(lanczos_args) :: ha
  integer :: i,j
  real(gp) Ene,gamma,  res


  type(gaussian_basis), target ::  Gabsorber
  real(wp),   pointer  :: Gabs_coeffs(:)
  real(wp), dimension(:), pointer  :: pot
  real(wp),  pointer, dimension(:,:)  :: dum_coeffs
  character(len=800) :: filename
  logical :: projeexists
  logical:: useold
  real(gp) , pointer ::potentialclone(:,:)

  if( iand( in%potshortcut,16)>0) then
     allocate(potentialclone(max(ndimpot,1),nspin+ndebug),stat=i_stat)
     call memocc(i_stat,potentialclone,'potentialclone',subname)
     potentialclone=potential
  endif

  if(iproc==0) print *, " IN ROUTINE xabs_cg "

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
 
  write(filename,'(A,A,A,I1)') "gproje_", trim(at%atomnames(at%iatype(  in_iat_absorber ))) , "_1s_",  in%L_absorber
  
  inquire(FILE=trim(filename),EXIST=projeexists)
  
  if( projeexists .and. .not. in%abscalc_eqdiff   ) then
     
     if(iproc==0) then
        print *, "reading  precalculated  projection on pseudofunctions"
        print *, "for 1s of atom number ", in_iat_absorber, " ( atomname = ", at%atomnames(at%iatype(  in_iat_absorber ))," )"
        print *, "After application of a 2*L-pole  with L= ", in%L_absorber
        print *," from file " , filename
     endif

     nullify( dum_coeffs  ) 
     nullify(Gabsorber%nshell)
     nullify(Gabsorber%nam)
     nullify(Gabsorber%ndoc)
     nullify(Gabsorber%xp)
     nullify(Gabsorber%psiat)
     
     call read_gaussian_information (ha%orbs,Gabsorber,dum_coeffs , filename, .true. )    
     Gabsorber%rxyz(:,1)=rxyz(:, in_iat_absorber )
     
     i_all=-product(shape(dum_coeffs))*kind(dum_coeffs)
     deallocate(dum_coeffs,stat=i_stat)
     call memocc(i_stat,i_all,'coeffs',subname)
     
     Gabs_coeffs(:)=in%Gabs_coeffs(:)
     
  else
     
        if(iproc==0) then
           print *, "calculating  projection on pseudofunctions"
           print *, "for 1s of atom number ", in_iat_absorber, " ( atomname = ", at%atomnames(at%iatype(  in_iat_absorber ))," )"
           print *, "After application of a 2*L-pole  with L= ", in%L_absorber
        endif
           call GetExcitedOrbitalAsG(in_iat_absorber ,Gabsorber,&
                at,rxyz,nproc,iproc,   in%L_absorber ,  in%abscalc_eqdiff )
           
           allocate(dum_coeffs(2*in%L_absorber+1,1+ndebug),stat=i_stat)
           call memocc(i_stat,dum_coeffs,'dum_coeffs',subname)
           
           dum_coeffs(:,1) = in%Gabs_coeffs(:)
           
           Gabs_coeffs(:)=in%Gabs_coeffs(:)
           
        if( iproc.eq.0) then
           print *,"writing them on file " , filename
           call write_gaussian_information( 0 ,1 ,ha%orbs,Gabsorber,dum_coeffs ,filename)
        endif
     
        i_all=-product(shape(dum_coeffs))*kind(dum_coeffs)
        deallocate(dum_coeffs,stat=i_stat)
        call memocc(i_stat,i_all,'coeffs',subname)
        
     endif

     !allocate the potential in the full box
     call full_local_potential(iproc,nproc,ndimpot,lr%d%n1i*lr%d%n2i*lr%d%n3i,in%nspin,&
          ha%orbs%norb,ha%orbs%norbp,ngatherarr,potential,pot)

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
  ha%lr=>lr !!!
  ha%ngatherarr=>ngatherarr
  ha%ndimpot=ndimpot
  ha%potential=>pot
  ha%ekin_sum=ekin_sum
  ha%epot_sum=epot_sum
  ha%eproj_sum=eproj_sum
  ha%nspin=nspin
  ha%GPU=>GPU !!
  ha%Gabsorber=>Gabsorber 
  ha%Gabs_coeffs=>Gabs_coeffs

  call EP_inizializza(ha) 

  call EP_memorizza_stato(Gabsorber) 
     
  if(.true.) then
     LB_nsteps =in%nsteps
     call LB_allocate_for_lanczos( )
     call EP_allocate_for_eigenprob(10)
     call EP_make_dummy_vectors(10)
     
         
     do i=0,70

        ene = 0.22_gp + i*0.03_gp

        if( iand( in%potshortcut,16)>0) then
           potential=potentialclone
           do j=1, ndimpot

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

        print *, "chiamo CG "
        res =  LB_cg(    get_EP_dim, EP_initialize_start , EP_normalizza,&
             EP_Moltiplica4spectra,  EP_copy,  &
             EP_scalare,EP_add_from_vect_with_fact   , EP_multbyfact  ,EP_precondition, Ene, gamma, 1.0D-4, useold )
        
        print *, ene, res
        open(unit=22,file="cgspectra", position='append')
        write(22,*) ene, res
        close(unit=22)
        
     enddo


     call LB_de_allocate_for_lanczos( )


  endif
  
  call deallocate_comms(ha%comms,subname)

  call EP_free()

  if (GPUconv) then
     call free_gpu(GPU,ha%orbs%norbp)
  end if

  call deallocate_orbs(ha%orbs,subname)

  i_all=-product(shape(Gabsorber%rxyz))*kind(Gabsorber%rxyz)
  deallocate(Gabsorber%rxyz,stat=i_stat)
  call memocc(i_stat,i_all,'Gabsorber%rxyz',subname)

  call deallocate_gwf(Gabsorber, subname)

  i_all=-product(shape(ha%orbs%eval))*kind(ha%orbs%eval)
  deallocate(ha%orbs%eval,stat=i_stat)
  call memocc(i_stat,i_all,'ha%orbs%eval',subname)

  i_all=-product(shape(Gabs_coeffs))*kind(Gabs_coeffs)
  deallocate(Gabs_coeffs,stat=i_stat)
  call memocc(i_stat,i_all,'Gabs_coeffs',subname)


  if( iand( in%potshortcut,16)>0) then
     i_all=-product(shape(potentialclone))*kind(potentialclone)
     deallocate(potentialclone,stat=i_stat)
     call memocc(i_stat,i_all,'potentialclone',subname)
  endif

  call free_full_potential(nproc,pot,subname)

  call deallocate_abscalc_input(in, subname)


END SUBROUTINE xabs_cg


subroutine dirac_hara (rho, E , V)
  use module_base
    
  implicit none

  real(gp), intent(in   ) :: rho, E
  real(gp), intent(inout) :: V
  real(gp), parameter :: f= 1.919158d0  !( 4/(9*pi)**(1/3)  )
  real(gp), parameter :: pi=3.141592653589793_gp
  
  real(gp) Vcorr, rs, xk, EV,x
  integer i

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
