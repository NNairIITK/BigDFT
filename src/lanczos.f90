
subroutine lanczos(iproc,nproc,at,hx,hy,hz,rxyz,Gabsorber,Gabs_coeffs,&
     cpmult,fpmult,radii_cf,nlpspd,proj,lr,ngatherarr,ndimpot,potential,&
     ekin_sum,epot_sum,eproj_sum,nspin,GPU)
  use module_base
  use module_types
  use lanczos_interface
  use lanczos_base
  implicit none
  integer  :: iproc,nproc,ndimpot,nspin
  real(gp)  :: hx,hy,hz,cpmult,fpmult
  type(atoms_data), target :: at
  type(nonlocal_psp_descriptors), target :: nlpspd
  type(locreg_descriptors), target :: lr
  type(gaussian_basis), target :: Gabsorber
  integer, dimension(0:nproc-1,2), target :: ngatherarr 
  real(gp), dimension(3,at%nat), target :: rxyz
  real(gp), dimension(at%ntypes,3), target :: radii_cf  
  real(wp), dimension(nlpspd%nprojel), target :: proj
  real(wp), dimension(max(ndimpot,1),nspin), target :: potential
  real(wp), dimension(3), target :: Gabs_coeffs
  real(gp) :: ekin_sum,epot_sum,eproj_sum
  type(GPU_pointers), intent(inout) , target :: GPU
  !local variables
  character(len=*), parameter :: subname='lanczos'
  integer :: i_stat,i_all
  type(lanczos_args) :: ha



 
  !create the orbitals descriptors, for virtual and inputguess orbitals



  allocate(ha%orbs%norb_par(0:nproc-1+ndebug),stat=i_stat)
  call memocc(i_stat,ha%orbs%norb_par,'ha%orbs%norb_par',subname)
  !davidson treatment for spin-pol case should be reworked

  call orbitals_descriptors(iproc,nproc,1,1,0,1,ha%orbs)

  !allocate communications arrays for virtual orbitals
  !warning: here the aim is just to calculate npsidim, should be fixed
  !call allocate_comms(nproc,ha%orbs,ha%comms,subname)
  call orbitals_communicators(iproc,nproc,lr,ha%orbs,ha%comms)  


  i_all=-product(shape(ha%orbs%norb_par))*kind(ha%orbs%norb_par)
  deallocate(ha%orbs%norb_par,stat=i_stat)
  call memocc(i_stat,i_all,'ha%orbs%norb_par',subname)
    

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

  print *, " chiamo inizializza " 
  call inizializza(ha)

  print *, " chiamo allocate for eigen prob"
  call EP_allocate_for_eigenprob(1)

  print *, " chiamo initialize start " 
  call EP_initialize_start(Gabsorber)

  print * , " uscito da EP_initialize_start "

  !deallocate communication and orbitals descriptors
  call deallocate_comms(ha%comms,subname)
  i_all=-product(shape(ha%orbs%occup))*kind(ha%orbs%occup)
  deallocate(ha%orbs%occup,stat=i_stat)
  call memocc(i_stat,i_all,'ha%orbs%occup',subname)
  i_all=-product(shape(ha%orbs%spinsgn))*kind(ha%orbs%spinsgn)
  deallocate(ha%orbs%spinsgn,stat=i_stat)
  call memocc(i_stat,i_all,'ha%orbs%spinsgn',subname)
  i_all=-product(shape(ha%orbs%kpts))*kind(ha%orbs%kpts)
  deallocate(ha%orbs%kpts,stat=i_stat)
  call memocc(i_stat,i_all,'ha%orbs%kpts',subname)
  i_all=-product(shape(ha%orbs%kwgts))*kind(ha%orbs%kwgts)
  deallocate(ha%orbs%kwgts,stat=i_stat)
  call memocc(i_stat,i_all,'ha%orbs%kwgts',subname)
  i_all=-product(shape(ha%orbs%iokpt))*kind(ha%orbs%iokpt)
  deallocate(ha%orbs%iokpt,stat=i_stat)
  call memocc(i_stat,i_all,'ha%orbs%iokpt',subname)

  
end subroutine lanczos
