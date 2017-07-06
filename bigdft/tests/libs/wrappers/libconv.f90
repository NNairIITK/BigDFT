!> test program for the convolution library
program libconv
  use futile
  use box
  use f_functions
  use locregs
  use locreg_operations
  use f_trees
  use BigDFT_API, only: bigdft_init_errors,bigdft_init_timing_categories
  use numerics
  implicit none
  integer :: i
  real(f_double) :: ekin,crmult,frmult,maxdiff
  type(cell) :: mesh
  type(box_iterator) :: bit
  type(locreg_descriptors) :: lr
  type(workarr_sumrho) :: w
  type(workarr_locham) :: wrk_lh
  logical, dimension(3) :: pers
  real(f_double), dimension(3) :: kpoint,oxyz,angrad,hgrids
  real(f_double), dimension(6) :: k_strten
  type(f_function), dimension(3) :: funcs
  type(f_tree) :: dict_posinp
  real(f_double), dimension(:), allocatable :: psir,psi,tpsi,psir_work,tpsir
  
  call f_lib_initialize()
 
  call bigdft_init_errors()
  call bigdft_init_timing_categories()

  hgrids=0.333_f_double
  dict_posinp=f_tree_load('{positions: [{ H: [0.0, 0.0, 0.0]}], cell: [10,15,11]}')
  crmult=50.0_f_double
  frmult=80.0_f_double
  angrad=onehalf*pi
  oxyz=5.0_f_double
  kpoint=0.0_f_double

  call define_lr(lr,dict_posinp,crmult,frmult,hgrids)

  call f_tree_free(dict_posinp)
  
  mesh=cell_new(lr%geocode,[lr%d%n1i,lr%d%n2i,lr%d%n3i],0.5_f_double*hgrids,&
       alpha_bc=angrad(1),beta_ac=angrad(2),gamma_ab=angrad(3))

  pers=cell_periodic_dims(mesh)
  do i=1,3
     if (pers(i)) then
        funcs(i)=f_function_new(f_exp_cosine,&
             length=mesh%ndims(i)*mesh%hgrids(i),frequency=2.0_f_double)
     else
        funcs(i)=f_function_new(f_shrink_gaussian,&
             length=mesh%ndims(1)*mesh%hgrids(1))
!!$        funcs(i)=f_function_new(f_gaussian,&
!!$             exponent=1.0_f_double/(10.0_f_double*mesh%hgrids(i)**2))
     end if
     !funcs(i)=f_function_new(f_constant,prefactor=11.0_f_double)
  end do

  call yaml_map('Number of grid points in the high-resolution domain',mesh%ndims)

  call yaml_map('Periodicity of the x,y,z dimensions',pers)

  psir=f_malloc(mesh%ndim,id='psir')
  tpsir=f_malloc(mesh%ndim,id='tpsir')
  psir_work=f_malloc0(mesh%ndim,id='psir_work')
  psi=f_malloc0(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f,id='psi')
  tpsi=f_malloc0(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f,id='tpsi')

  bit=box_iter(mesh,centered=.true.)
  !take the reference functions
  call separable_3d_function(bit,funcs,1.0_f_double,psir)
  call separable_3d_laplacian(bit,funcs,-2.0_f_double,tpsir)

  !initalize workspaces
  call initialize_work_arrays_sumrho(lr,.true.,w)
  call initialize_work_arrays_locham(lr,1,.true.,wrk_lh)

  !from real space to wavelet
  call isf_to_daub(lr,w,psir,psi)

  !from wavelets to real space and fine scaling functions space
  call daub_to_isf_locham(1,lr,wrk_lh,psi,psir_work)

  !multiply by the potential (zero in this case)
  call f_zero(psir_work)

  !calculate results of the laplacian
  call isf_to_daub_kinetic_clone(mesh%hgrids(1),mesh%hgrids(2),mesh%hgrids(3),&
       kpoint(1),kpoint(2),kpoint(3),1,lr,wrk_lh,&
       psir_work,tpsi,ekin,k_strten)

  !from wavelets to real space
  call daub_to_isf(lr,w,tpsi,psir_work)

  !free work arrays
  call deallocate_work_arrays_sumrho(w)
  call deallocate_work_arrays_locham(wrk_lh)
  call deallocate_locreg_descriptors(lr)

  !compare the results (tpsir with psir_work) 
  call f_diff(size(tpsir),tpsir,psir_work,maxdiff,ind=i)

  call yaml_map('Final difference',maxdiff)
  call yaml_map('At point',i)
!  call yaml_map('Corresponding to value',ind_to_iarr(mesh%ndims,i))
  call yaml_map('Calculated value',psir_work(i))
  call yaml_map('Reference value',tpsir(i))

  call f_free(psir,tpsir,psir_work)
  call f_free(psi,tpsi)

!!$  !list of routines which have to be duplicated and modified
!!$  call convolut_kinetic_slab_T_k(2*lr%d%n1+1,2*lr%d%n2+15,2*lr%d%n3+1,&
!!$                    hgridh,w%x_c(1,idx),w%y_c(1,idx),ekino,kx,ky,kz)
!!$
!!$ ! compute the kinetic part and add  it to psi_out
!!$ ! the kinetic energy is calculated at the same time
!!$ call convolut_kinetic_slab_T(2*lr%d%n1+1,2*lr%d%n2+15,2*lr%d%n3+1,&
!!$      hgridh,w%x_c(1,idx),w%y_c(1,idx),ekino)
!!$
!!$ call convolut_kinetic_hyb_T(lr%d%n1,lr%d%n2,lr%d%n3, &
!!$      lr%d%nfl1,lr%d%nfu1,lr%d%nfl2,lr%d%nfu2,lr%d%nfl3,lr%d%nfu3,  &
!!$      hgridh,w%x_c(1,idx),w%x_f(1,idx),w%y_c(1,idx),w%y_f(1,idx),kstrteno,&
!!$      w%x_f1(1,idx),w%x_f2(1,idx),w%x_f3(1,idx),lr%bounds%kb%ibyz_f,&
!!$      lr%bounds%kb%ibxz_f,lr%bounds%kb%ibxy_f)
!!$
!!$ call convolut_kinetic_per_T_k(2*lr%d%n1+1,2*lr%d%n2+1,2*lr%d%n3+1,&
!!$      hgridh,w%x_c(1,idx),w%y_c(1,idx),kstrteno,kx,ky,kz)
!!$ 
!!$ call convolut_kinetic_per_t(2*lr%d%n1+1,2*lr%d%n2+1,2*lr%d%n3+1,&
!!$      hgridh,w%x_c(1,idx),w%y_c(1,idx),kstrteno)
  
  call f_lib_finalize()

contains

  function ind_to_iarr(ndims,ind) result(iarr)
    implicit none
    integer, dimension(:), intent(in) :: ndims
    integer, intent(in) :: ind
    integer, dimension(size(ndims)) :: iarr
    !local variables
    integer :: j,jnd,n,m

    n=size(ndims)   
    jnd=ind
    m=1
    do j=1,n-1
       iarr(j)=modulo(jnd-1,ndims(j+1))+1
       jnd=jnd-(iarr(j)-1)*m
       if (j<n-1) jnd=jnd/ndims(j+1)
       m=m*ndims(j)
    end do
    iarr(n)=modulo(jnd-1,ndims(n))+1
  end function ind_to_iarr


end program libconv


subroutine define_lr(lr,tree_posinp,crmult,frmult,hgrids)
  use module_atoms
  use locregs
  use locregs_init
  use futile
  use f_trees
  use pseudopotentials
  implicit none
  type(f_tree), intent(in) :: tree_posinp
  type(locreg_descriptors), intent(out) :: lr
  real(f_double), intent(in) :: crmult,frmult
  real(f_double), dimension(3), intent(in) :: hgrids
  !local variables
  type(atoms_data) :: atoms
  type(dictionary), pointer :: dict_targets,dict_dm,types,var
  real(f_double), dimension(:,:), allocatable :: rxyz

  call dict_init(dict_dm)
  call dict_init(dict_targets)
  call nullify_atoms_data(atoms)

  !> fill the atomic structure datatype
  call astruct_set(atoms%astruct,tree_posinp%d,0.0_f_double,.true.,1.e-8_f_double,&
       [0.0_f_double,0.0_f_double,0.0_f_double],1,.true.)

  ! Generate the dict of types for later use.
  call astruct_dict_get_types(tree_posinp%d, types)

  nullify(var)
  do while(iterating(var,on=types))
     call psp_dict_fill_all(dict_dm, trim(dict_key(var)), 1, 8.0_f_double,crmult,frmult)
  end do
  call dict_free(types)

  call atoms_fill(atoms,dict_dm,frmult,1,.false.,16,0,0.0_f_double)

  rxyz=f_malloc(src=atoms%astruct%rxyz,id='rxyz')

  call lr_set(lr,0,.false.,.true.,crmult,frmult,hgrids,rxyz,atoms,.true.,.false.)

  call f_free(rxyz)

  call deallocate_atoms_data(atoms)
  call dict_free(dict_dm,dict_targets)

end subroutine define_lr


subroutine isf_to_daub_kinetic_clone(hx,hy,hz,kx,ky,kz,nspinor,lr,w,psir,hpsi,ekin,k_strten)
  use locregs
  use locreg_operations
  use module_defs
  implicit none
  integer, intent(in) :: nspinor
  real(gp), intent(in) :: hx,hy,hz,kx,ky,kz
  type(locreg_descriptors), intent(in) :: lr
  type(workarr_locham), intent(inout) :: w
  real(wp), dimension(lr%d%n1i*lr%d%n2i*lr%d%n3i,nspinor), intent(inout) :: psir
  real(gp), intent(out) :: ekin
  real(wp), dimension(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f,nspinor), intent(inout) :: hpsi
  real(wp), dimension(6) :: k_strten
  !Local variables
  logical :: usekpts
  integer :: idx,i,i_f,iseg_f,ipsif,isegf
  real(gp) :: ekino
  real(wp), dimension(0:3) :: scal
  real(gp), dimension(3) :: hgridh
  real(wp), dimension(6) :: kstrten,kstrteno


  !control whether the k points are to be used
  !real k-point different from Gamma still not implemented
  usekpts = kx**2+ky**2+kz**2 > 0.0_gp .or. nspinor == 2

  hgridh(1)=hx*.5_gp
  hgridh(2)=hy*.5_gp
  hgridh(3)=hz*.5_gp

  do i=0,3
     scal(i)=1.0_wp
  enddo

  !starting point for the fine degrees, to avoid boundary problems
  i_f=min(1,lr%wfd%nvctr_f)
  iseg_f=min(1,lr%wfd%nseg_f)
  ipsif=lr%wfd%nvctr_c+i_f
  isegf=lr%wfd%nseg_c+iseg_f

  !call MPI_COMM_RANK(bigdft_mpi%mpi_comm,iproc,ierr)
  ekin=0.0_gp

  kstrten=0.0_wp
  select case(lr%geocode)
  case('F')

     !here kpoints cannot be used (for the moment, to be activated for the 
     !localisation region scheme
     if (usekpts) stop 'K points not allowed for Free BC locham'

     do idx=1,nspinor

        call comb_shrink(lr%d%n1,lr%d%n2,lr%d%n3,&
             lr%d%nfl1,lr%d%nfu1,lr%d%nfl2,lr%d%nfu2,lr%d%nfl3,lr%d%nfu3,&
             w%w1,w%w2,psir(1,idx),&
             lr%bounds%kb%ibxy_c,lr%bounds%sb%ibzzx_c,lr%bounds%sb%ibyyzz_c,&
             lr%bounds%sb%ibxy_ff,lr%bounds%sb%ibzzx_f,lr%bounds%sb%ibyyzz_f,&
             w%y_c(1,idx),w%y_f(1,idx))

        call ConvolkineticT(lr%d%n1,lr%d%n2,lr%d%n3,&
             lr%d%nfl1,lr%d%nfu1,lr%d%nfl2,lr%d%nfu2,lr%d%nfl3,lr%d%nfu3,  &
             hx,hy,hz,&      !here the grid spacings are supposed to be equal.  SM: not any more
             lr%bounds%kb%ibyz_c,lr%bounds%kb%ibxz_c,lr%bounds%kb%ibxy_c,&
             lr%bounds%kb%ibyz_f,lr%bounds%kb%ibxz_f,lr%bounds%kb%ibxy_f, &
             w%x_c(1,idx),w%x_f(1,idx),&
             w%y_c(1,idx),w%y_f(1,idx),ekino, &
             w%x_f1(1,idx),w%x_f2(1,idx),w%x_f3(1,idx),111)
        ekin=ekin+ekino

        !new compression routine in standard form
        call compress_and_accumulate_standard(lr%d,lr%wfd,&
             lr%wfd%keyvloc(1),lr%wfd%keyvloc(isegf),&
             lr%wfd%keygloc(1,1),lr%wfd%keygloc(1,isegf),&
             w%y_c(1,idx),w%y_f(1,idx),&
             hpsi(1,idx),hpsi(ipsif,idx))
!!$        call compress_forstandard(lr%d%n1,lr%d%n2,lr%d%n3,&
!!$             lr%d%nfl1,lr%d%nfu1,lr%d%nfl2,lr%d%nfu2,lr%d%nfl3,lr%d%nfu3,  &
!!$             lr%wfd%nseg_c,lr%wfd%nvctr_c,&
!!$             lr%wfd%keygloc(1,1),lr%wfd%keyv(1),&
!!$             lr%wfd%nseg_f,lr%wfd%nvctr_f,&
!!$             lr%wfd%keygloc(1,lr%wfd%nseg_c+iseg_f),lr%wfd%keyv(lr%wfd%nseg_c+iseg_f),   &
!!$             scal,w%y_c(1,idx),w%y_f(1,idx),hpsi(1,idx),hpsi(lr%wfd%nvctr_c+i_f,idx))

     end do

  case('S')

     if (usekpts) then
        !first calculate the proper arrays then transpose them before passing to the
        !proper routine
        do idx=1,nspinor
           call convolut_magic_t_slab_self(2*lr%d%n1+1,2*lr%d%n2+15,2*lr%d%n3+1,&
                psir(1,idx),w%y_c(1,idx))
        end do

        !Transposition of the work arrays (use psir as workspace)
        call transpose_for_kpoints(nspinor,2*lr%d%n1+2,2*lr%d%n2+31,2*lr%d%n3+2,&
             w%x_c,psir,.true.)
        call transpose_for_kpoints(nspinor,2*lr%d%n1+2,2*lr%d%n2+31,2*lr%d%n3+2,&
             w%y_c,psir,.true.)

        ! compute the kinetic part and add  it to psi_out
        ! the kinetic energy is calculated at the same time
        ! do this thing for both components of the spinors
        do idx=1,nspinor,2
           call convolut_kinetic_slab_T_k(2*lr%d%n1+1,2*lr%d%n2+15,2*lr%d%n3+1,&
                hgridh,w%x_c(1,idx),w%y_c(1,idx),ekino,kx,ky,kz)
           ekin=ekin+ekino        
        end do

        !re-Transposition of the work arrays (use psir as workspace)
        call transpose_for_kpoints(nspinor,2*lr%d%n1+2,2*lr%d%n2+31,2*lr%d%n3+2,&
             w%y_c,psir,.false.)

        do idx=1,nspinor
           !new compression routine in mixed form
           call analyse_slab_self(lr%d%n1,lr%d%n2,lr%d%n3,&
                w%y_c(1,idx),psir(1,idx))
           call compress_and_accumulate_mixed(lr%d,lr%wfd,&
                lr%wfd%keyvloc(1),lr%wfd%keyvloc(isegf),&
                lr%wfd%keygloc(1,1),lr%wfd%keygloc(1,isegf),&
                psir(1,idx),hpsi(1,idx),hpsi(ipsif,idx))

!!$           call compress_slab(lr%d%n1,lr%d%n2,lr%d%n3,&
!!$                lr%wfd%nseg_c,lr%wfd%nvctr_c,&
!!$                lr%wfd%keygloc(1,1),lr%wfd%keyv(1),   & 
!!$                lr%wfd%nseg_f,lr%wfd%nvctr_f,&
!!$                lr%wfd%keygloc(1,lr%wfd%nseg_c+iseg_f),lr%wfd%keyv(lr%wfd%nseg_c+iseg_f),   & 
!!$                w%y_c(1,idx),hpsi(1,idx),hpsi(lr%wfd%nvctr_c+i_f,idx),psir(1,idx))
        end do

     else
        do idx=1,nspinor
           call convolut_magic_t_slab_self(2*lr%d%n1+1,2*lr%d%n2+15,2*lr%d%n3+1,&
                psir(1,idx),w%y_c(1,idx))

           ! compute the kinetic part and add  it to psi_out
           ! the kinetic energy is calculated at the same time
           call convolut_kinetic_slab_T(2*lr%d%n1+1,2*lr%d%n2+15,2*lr%d%n3+1,&
                hgridh,w%x_c(1,idx),w%y_c(1,idx),ekino)
           ekin=ekin+ekino

           !new compression routine in mixed form
           call analyse_slab_self(lr%d%n1,lr%d%n2,lr%d%n3,&
                w%y_c(1,idx),psir(1,idx))
           call compress_and_accumulate_mixed(lr%d,lr%wfd,&
                lr%wfd%keyvloc(1),lr%wfd%keyvloc(isegf),&
                lr%wfd%keygloc(1,1),lr%wfd%keygloc(1,isegf),&
                psir(1,idx),hpsi(1,idx),hpsi(ipsif,idx))

!!$           call compress_slab(lr%d%n1,lr%d%n2,lr%d%n3,&
!!$                lr%wfd%nseg_c,lr%wfd%nvctr_c,&
!!$                lr%wfd%keygloc(1,1),lr%wfd%keyv(1),   & 
!!$                lr%wfd%nseg_f,lr%wfd%nvctr_f,&
!!$                lr%wfd%keygloc(1,lr%wfd%nseg_c+iseg_f),lr%wfd%keyv(lr%wfd%nseg_c+iseg_f),   & 
!!$                w%y_c(1,idx),hpsi(1,idx),hpsi(lr%wfd%nvctr_c+i_f,idx),psir(1,idx))
        end do
     end if

  case('P')

     if (lr%hybrid_on) then

        !here kpoints cannot be used, such BC are used in general to mimic the Free BC
        if (usekpts) stop 'K points not allowed for hybrid BC locham'

        !here the grid spacing is not halved
        hgridh(1)=hx
        hgridh(2)=hy
        hgridh(3)=hz
        do idx=1,nspinor
           call comb_shrink_hyb(lr%d%n1,lr%d%n2,lr%d%n3,&
                lr%d%nfl1,lr%d%nfu1,lr%d%nfl2,lr%d%nfu2,lr%d%nfl3,lr%d%nfu3,&
                w%w2,w%w1,psir(1,idx),w%y_c(1,idx),w%y_f(1,idx),lr%bounds%sb)

           call convolut_kinetic_hyb_T(lr%d%n1,lr%d%n2,lr%d%n3, &
                lr%d%nfl1,lr%d%nfu1,lr%d%nfl2,lr%d%nfu2,lr%d%nfl3,lr%d%nfu3,  &
                hgridh,w%x_c(1,idx),w%x_f(1,idx),w%y_c(1,idx),w%y_f(1,idx),kstrteno,&
                w%x_f1(1,idx),w%x_f2(1,idx),w%x_f3(1,idx),lr%bounds%kb%ibyz_f,&
                lr%bounds%kb%ibxz_f,lr%bounds%kb%ibxy_f)
           kstrten=kstrten+kstrteno
           !ekin=ekin+ekino

           call compress_and_accumulate_standard(lr%d,lr%wfd,&
                lr%wfd%keyvloc(1),lr%wfd%keyvloc(isegf),&
                lr%wfd%keygloc(1,1),lr%wfd%keygloc(1,isegf),&
                w%y_c(1,idx),w%y_f(1,idx),hpsi(1,idx),hpsi(ipsif,idx))

!!$           call compress_per_f(lr%d%n1,lr%d%n2,lr%d%n3,&
!!$                lr%wfd%nseg_c,lr%wfd%nvctr_c,&
!!$                lr%wfd%keygloc(1,1),lr%wfd%keyv(1),& 
!!$                lr%wfd%nseg_f,lr%wfd%nvctr_f,&
!!$                lr%wfd%keygloc(1,lr%wfd%nseg_c+iseg_f),lr%wfd%keyv(lr%wfd%nseg_c+iseg_f), & 
!!$                w%y_c(1,idx),w%y_f(1,idx),hpsi(1,idx),hpsi(lr%wfd%nvctr_c+i_f,idx),&
!!$                lr%d%nfl1,lr%d%nfl2,lr%d%nfl3,lr%d%nfu1,lr%d%nfu2,lr%d%nfu3)
        end do
     else
        if (usekpts) then
           !first calculate the proper arrays then transpose them before passing to the
           !proper routine
           do idx=1,nspinor
              call convolut_magic_t_per_self(2*lr%d%n1+1,2*lr%d%n2+1,2*lr%d%n3+1,&
                   psir(1,idx),w%y_c(1,idx))
           end do

           !Transposition of the work arrays (use psir as workspace)
           call transpose_for_kpoints(nspinor,2*lr%d%n1+2,2*lr%d%n2+2,2*lr%d%n3+2,&
                w%x_c,psir,.true.)
           call transpose_for_kpoints(nspinor,2*lr%d%n1+2,2*lr%d%n2+2,2*lr%d%n3+2,&
                w%y_c,psir,.true.)


           ! compute the kinetic part and add  it to psi_out
           ! the kinetic energy is calculated at the same time
           do idx=1,nspinor,2
              !print *,'AAA',2*lr%d%n1+1,2*lr%d%n2+1,2*lr%d%n3+1,hgridh

              call convolut_kinetic_per_T_k(2*lr%d%n1+1,2*lr%d%n2+1,2*lr%d%n3+1,&
                   hgridh,w%x_c(1,idx),w%y_c(1,idx),kstrteno,kx,ky,kz)
              kstrten=kstrten+kstrteno
              !ekin=ekin+ekino
           end do

           !Transposition of the work arrays (use psir as workspace)
           call transpose_for_kpoints(nspinor,2*lr%d%n1+2,2*lr%d%n2+2,2*lr%d%n3+2,&
                w%y_c,psir,.false.)

           do idx=1,nspinor

              call analyse_per_self(lr%d%n1,lr%d%n2,lr%d%n3,&
                   w%y_c(1,idx),psir(1,idx))
              call compress_and_accumulate_mixed(lr%d,lr%wfd,&
                   lr%wfd%keyvloc(1),lr%wfd%keyvloc(isegf),&
                   lr%wfd%keygloc(1,1),lr%wfd%keygloc(1,isegf),&
                   psir(1,idx),hpsi(1,idx),hpsi(ipsif,idx))

!!$              call compress_per(lr%d%n1,lr%d%n2,lr%d%n3,&
!!$                   lr%wfd%nseg_c,lr%wfd%nvctr_c,&
!!$                   lr%wfd%keygloc(1,1),lr%wfd%keyv(1),& 
!!$                   lr%wfd%nseg_f,lr%wfd%nvctr_f,&
!!$                   lr%wfd%keygloc(1,lr%wfd%nseg_c+iseg_f),lr%wfd%keyv(lr%wfd%nseg_c+iseg_f),&
!!$                   w%y_c(1,idx),hpsi(1,idx),hpsi(lr%wfd%nvctr_c+i_f,idx),psir(1,idx))
           end do
        else
           !first calculate the proper arrays then transpose them before passing to the
           !proper routine
           do idx=1,nspinor
!!$              call convolut_magic_t_per_self(2*lr%d%n1+1,2*lr%d%n2+1,2*lr%d%n3+1,&
!!$                   psir(1,idx),w%y_c(1,idx))

              call d_s0s0_1d_sym8_imd(3,0,[2*lr%d%n1+2,2*lr%d%n2+2,2*lr%d%n3+2],0,&
                   [2*lr%d%n1+2,2*lr%d%n2+2,2*lr%d%n3+2],&
                   [2*lr%d%n1+2,2*lr%d%n2+2,2*lr%d%n3+2],1,&
                   psir(1,idx),w%y_c(1,idx),&
                   1.0_wp, 1.0_wp,0.0_wp)
              call d_s0s0_1d_sym8_imd(3,1,[2*lr%d%n1+2,2*lr%d%n2+2,2*lr%d%n3+2],0,&
                   [2*lr%d%n1+2,2*lr%d%n2+2,2*lr%d%n3+2],&
                   [2*lr%d%n1+2,2*lr%d%n2+2,2*lr%d%n3+2],1,&
                   w%y_c(1,idx),psir(1,idx),&
                   1.0_wp, 1.0_wp,0.0_wp)
              call d_s0s0_1d_sym8_imd(3,2,[2*lr%d%n1+2,2*lr%d%n2+2,2*lr%d%n3+2],0,&
                   [2*lr%d%n1+2,2*lr%d%n2+2,2*lr%d%n3+2],&
                   [2*lr%d%n1+2,2*lr%d%n2+2,2*lr%d%n3+2],1,&
                   psir(1,idx),w%y_c(1,idx),&
                   1.0_wp, 1.0_wp,0.0_wp)

              

!!$              call d_s0s0_1d_sym8_imd(3,0,&
!!$                   2*lr%d%n1+2,2*lr%d%n2+2,2*lr%d%n3+2,&
!!$                   [0,0,0],&
!!$                   2*lr%d%n1+2,2*lr%d%n2+2,2*lr%d%n3+2,1,&
!!$                   psir(1,idx),w%y_c(1,idx),&
!!$                   1.0_wp,1.0_wp,0.0_wp)

              ! compute the kinetic part and add  it to psi_out
              ! the kinetic energy is calculated at the same time
              call convolut_kinetic_per_t(2*lr%d%n1+1,2*lr%d%n2+1,2*lr%d%n3+1,&
                   hgridh,w%x_c(1,idx),w%y_c(1,idx),kstrteno)
              kstrten=kstrten+kstrteno

              call switch_s0_for_libconv(lr%d%n1+1,lr%d%n2+1,lr%d%n3+1,w%y_c,psir)

              call d_s0s1_1d_sym8(3,0,[lr%d%n1+1,lr%d%n2+1,lr%d%n3+1],0,&
                   [lr%d%n1+1,lr%d%n2+1,lr%d%n3+1],&
                   [lr%d%n1+1,lr%d%n2+1,lr%d%n3+1],1,&
                   w%y_c(1,idx),psir(1,idx),&
                   1.0_wp,0.0_wp)
              call d_s0s1_1d_sym8(3,1,[lr%d%n1+1,lr%d%n2+1,lr%d%n3+1],0,&
                   [lr%d%n1+1,lr%d%n2+1,lr%d%n3+1],&
                   [lr%d%n1+1,lr%d%n2+1,lr%d%n3+1],1,&
                   psir(1,idx),w%y_c(1,idx),&
                   1.0_wp,0.0_wp)
              call d_s0s1_1d_sym8(3,2,[lr%d%n1+1,lr%d%n2+1,lr%d%n3+1],0,&
                   [lr%d%n1+1,lr%d%n2+1,lr%d%n3+1],&
                   [lr%d%n1+1,lr%d%n2+1,lr%d%n3+1],1,&
                   w%y_c(1,idx),psir(1,idx),&
                   1.0_wp,0.0_wp)

!!$              call analyse_per_self(lr%d%n1,lr%d%n2,lr%d%n3,&
!!$                   w%y_c(1,idx),psir(1,idx))
              call compress_and_accumulate_mixed(lr%d,lr%wfd,&
                   lr%wfd%keyvloc(1),lr%wfd%keyvloc(isegf),&
                   lr%wfd%keygloc(1,1),lr%wfd%keygloc(1,isegf),&
                   psir(1,idx),hpsi(1,idx),hpsi(ipsif,idx))

!!$              call compress_per(lr%d%n1,lr%d%n2,lr%d%n3,&
!!$                   lr%wfd%nseg_c,lr%wfd%nvctr_c,&
!!$                   lr%wfd%keygloc(1,1),lr%wfd%keyv(1),& 
!!$                   lr%wfd%nseg_f,lr%wfd%nvctr_f,&
!!$                   lr%wfd%keygloc(1,lr%wfd%nseg_c+iseg_f),lr%wfd%keyv(lr%wfd%nseg_c+iseg_f),& 
!!$                   w%y_c(1,idx),hpsi(1,idx),hpsi(lr%wfd%nvctr_c+i_f,idx),psir(1,idx))
           end do
        end if

     end if
     ekin=ekin+kstrten(1)+kstrten(2)+kstrten(3)
     k_strten=kstrten 

  end select

END SUBROUTINE isf_to_daub_kinetic_clone

subroutine switch_s0_for_libconv(n1,n2,n3,x,y)
  use f_precisions, only: f_double
  use dynamic_memory, only: f_memcpy
  implicit none
  integer, intent(in) :: n1,n2,n3
  real(f_double), dimension(2*n1,2*n2,2*n3), intent(inout) :: x
  real(f_double), dimension(2*n1,2*n2,2*n3), intent(inout) :: y
  !local variables
  integer :: i1,i2,i3

  do i3=1,2*n3
     do i2=1,2*n2
        do i1=1,2*n1
           y(i1,i2,i3)=x(modulo(i1,2*n1)+1,modulo(i2,2*n2)+1,modulo(i3,2*n3)+1)
        end do
     end do
  end do
  call f_memcpy(src=y,dest=x)

end subroutine switch_s0_for_libconv
