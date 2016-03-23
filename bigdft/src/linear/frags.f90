!needs cleaning once we stabilize which options are useful
subroutine fragment_coeffs_to_kernel(iproc,input,input_frag_charge,ref_frags,tmb,ksorbs,overlap_calculated,&
  nstates_max,cdft,restart_mode,rmax)
  use yaml_output
  use module_base
  use module_types
  use module_fragments
  use communications_base, only: TRANSPOSE_FULL
  use communications, only: transpose_localized
  use sparsematrix_base, only: sparsematrix_malloc_ptr, DENSE_FULL, assignment(=)
  use sparsematrix, only: uncompress_matrix, gather_matrix_from_taskgroups_inplace, &
                          uncompress_matrix2
  use transposed_operations, only: calculate_overlap_transposed
  use module_interfaces, only: write_eigenvalues_data
  use public_enums
  implicit none
  type(DFT_wavefunction), intent(inout) :: tmb
  type(input_variables), intent(in) :: input
  type(system_fragment), dimension(input%frag%nfrag_ref), intent(inout) :: ref_frags
  type(orbitals_data), intent(inout) :: ksorbs
  logical, intent(inout) :: overlap_calculated
  real(kind=gp), dimension(input%frag%nfrag), intent(in) :: input_frag_charge
  integer, intent(in) :: iproc
  integer, intent(out) :: nstates_max ! number of states in total if we consider all partially occupied fragment states to be fully occupied
  logical, intent(in) :: cdft
  integer, intent(in) :: restart_mode
  real(kind=gp), intent(in) :: rmax

  integer :: iorb, isforb, jsforb, ifrag, ifrag_ref, itmb, jtmb, num_extra_per_frag, linstate, jf, pm, ortho_size, s, nelecfrag
  integer, allocatable, dimension(:) :: ipiv
  real(gp), dimension(:,:), allocatable :: coeff_final
  real(gp) :: nelecorbs, nelecfrag_tot, jstate_max, homo_diff, lag_mult, fac, cdft_charge
  real(gp), dimension(:), allocatable :: eval_tmp, eval_tmp2
  character(len=*), parameter :: subname='fragment_coeffs_to_kernel'

  integer :: rand_size
  integer, allocatable, dimension(:) :: rand_seed
  real(kind=dp) :: rtime, random_noise
  character(len=10) :: sys_time
  logical :: random, lincombm, lincombp !

  real(wp), dimension(:,:,:), pointer :: mom_vec_fake
  logical :: completely_random, use_tmbs_as_coeffs


  call timing(iproc,'kernel_init','ON')
  call f_routine(id='fragment_coeffs_to_kernel')

  ! need to do this properly/rearrange routines
  if (cdft) then
     ! otherwise doesn't make sense
     if (input%frag%nfrag_ref==2) homo_diff=(ref_frags(1)%eval(ceiling(ref_frags(1)%nelec/2.0_gp))&
          -ref_frags(2)%eval(ceiling(ref_frags(2)%nelec/2.0_gp)))/2.0d0
     !if (cdft%charge<0) lag_mult=-0.5, otherwise +0.5
     !lag_mult=-0.05d0

     !kind of unfortunate that we have to recreate this information here -
     !should maybe move the cdft initialization?
     cdft_charge=input%frag%charge(1)!(cdft%ifrag_charged(1))
     if (cdft_charge<0) then
        lag_mult=-abs(input%lin%cdft_lag_mult_init)
     else
        lag_mult=abs(input%lin%cdft_lag_mult_init)
     end if
  else
     homo_diff=0.0d0
     lag_mult=0.0d0
  end if

  ! adding random noise to starting to help with local minima problem
  completely_random=.false. ! completely random start for coeffs
  use_tmbs_as_coeffs=.false.
  if (restart_mode==LIN_RESTART_RANDOM) then
     completely_random=.true.
  else if (restart_mode==LIN_RESTART_DIAG_KERNEL) then
     use_tmbs_as_coeffs=.true.
     stop 'Error in restart: cannot construct coefficients using diagonal kernel'
  else if (restart_mode==LIN_RESTART_KERNEL) then
     stop 'Error in restart: annot construct coefficients using kernel'
  end if

  random=(rmax>0.0d0) ! add a bit of noise
  random_noise=0.0d0
  rtime=0.0d0
  if (random .or. completely_random) then
     call random_seed(size=rand_size)
     allocate(rand_seed(1:rand_size))
     call date_and_time(time=sys_time)
     ! coeffs need to be the same across processors
     if (iproc==0) read(sys_time,*) rtime
     if (bigdft_mpi%nproc>1) call mpiallred(rtime, 1, mpi_sum, comm=bigdft_mpi%mpi_comm)
     rand_seed=int(rtime*1000.0_dp)
     call random_seed(put=rand_seed)
     deallocate(rand_seed) 
  end if

  nelecfrag_tot=0
  do ifrag=1,input%frag%nfrag
     ifrag_ref=input%frag%frag_index(ifrag)
     nelecfrag_tot=nelecfrag_tot+ref_frags(ifrag_ref)%nelec-input_frag_charge(ifrag)
  end do

  if (completely_random) then
     if (bigdft_mpi%iproc==0) print*,'Starting coeffs are replaced with a random guess'
  else if (random) then
     if (bigdft_mpi%iproc==0) print*,'Random noise added to starting coeffs'
  end if

  ! in theory we could add/remove states depending on their energies, but for now we force the user to specify
  ! need to include occupations as we actually want to compare number of electrons here?
  nelecorbs=0
  do iorb=1,ksorbs%norb
     nelecorbs=nelecorbs+ksorbs%occup(iorb)
  end do

  lincombm=.false.
  lincombp=.false.
  !if (nint(nelecorbs)/=nelecfrag_tot) then
  if (abs(nelecorbs-nelecfrag_tot)>0.01d0) then
     !EXPERIMENTAL
     !use an average (lin. combination) of the LUMOs from each fragment
     !take the easiest case of 1 extra or missing electron, to be generalized/improved
     lincombm=.false.
     lincombp=.false.
     if (nelecorbs-nelecfrag_tot==1.0d0) then
        lincombm=.true. 
     else if (nelecorbs-nelecfrag_tot==-1.0d0) then
        lincombp=.true.
     else
        print*,'User should specify which fragments charges are added to/removed from in charged fragment calculation',&
             nelecfrag_tot,nelecorbs,ksorbs%norb
        stop
     end if
     if (iproc==0) print*,'Warning, experimental guess for unconstrained charged fragment calculation, proceed with caution',&
             nelecfrag_tot,nelecorbs,ksorbs%norb
  end if

  if (mod(input%norbsempty,input%frag%nfrag)/=0) then
     if (bigdft_mpi%iproc==0) call yaml_warning('Number of extra bands does not divide evenly among fragments')
     !if (bigdft_mpi%iproc==0) print*,'Warning, number of extra bands does not divide evenly among fragments'
     num_extra_per_frag=(input%norbsempty-mod(input%norbsempty,input%frag%nfrag))/input%frag%nfrag
  else
     num_extra_per_frag=input%norbsempty/input%frag%nfrag
  end if

  ! Calculate the overlap matrix between the TMBs.
  if(.not. overlap_calculated) then
     call timing(iproc,'kernel_init','OF')
     if(.not.tmb%can_use_transposed) then
         if(associated(tmb%psit_c)) then
             call f_free_ptr(tmb%psit_c)
         end if
         if(associated(tmb%psit_f)) then
             call f_free_ptr(tmb%psit_f)
         end if
         tmb%psit_c = f_malloc_ptr(sum(tmb%collcom%nrecvcounts_c),id='tmb%psit_c')
         tmb%psit_f = f_malloc_ptr(7*sum(tmb%collcom%nrecvcounts_f),id='tmb%psit_f')
         call transpose_localized(bigdft_mpi%iproc, bigdft_mpi%nproc, tmb%npsidim_orbs, tmb%orbs, tmb%collcom, &
              TRANSPOSE_FULL, tmb%psi, tmb%psit_c, tmb%psit_f, tmb%lzd)
         tmb%can_use_transposed=.true.
     end if

     call calculate_overlap_transposed(bigdft_mpi%iproc, bigdft_mpi%nproc, tmb%orbs, tmb%collcom, &
          tmb%psit_c, tmb%psit_c, tmb%psit_f, tmb%psit_f, tmb%linmat%s, tmb%linmat%ovrlp_)
     !!call gather_matrix_from_taskgroups_inplace(bigdft_mpi%iproc, bigdft_mpi%nproc, tmb%linmat%s, tmb%linmat%ovrlp_)
     ! This can then be deleted if the transition to the new type has been completed.
     !tmb%linmat%ovrlp%matrix_compr=tmb%linmat%ovrlp_%matrix_compr

     overlap_calculated=.true.
     call timing(iproc,'kernel_init','ON')
  end if

  ! set coeffs to tmbs with occupancies so each is equally weighted
  ! can only easily be used to get kernel (via FOE) or using diag
  ! occupancies should be reset afterwards
  if (.not. use_tmbs_as_coeffs .or. completely_random) then
     call fill_occupied_coeffs()

     ! still experimental and only for specific cases
     if (lincombp .or. lincombm) then
        call fill_inbetween_coeffs()
     end if

     call fill_unoccupied_coeffs()

     if (.not. completely_random) then
        call reorder_and_print_coeffs()
     end if
     call f_free(eval_tmp)
  else
     call set_coeffs_to_tmbs()

     if (random) call add_random_noise()
  end if


  call f_release_routine()
  call timing(iproc,'kernel_init','OF')

contains

  ! used for diagonal case only - in other case we orthonormalize in sections so easier to keep it internal
  subroutine add_random_noise()
    implicit none

    do itmb=1,tmb%orbs%norb
       do jtmb=1,tmb%orbs%norb
          call random_number(random_noise)
          random_noise = ((random_noise-0.5d0)*2.0d0)*rmax
          tmb%coeff(itmb,jtmb) = tmb%coeff(itmb,jtmb) + random_noise
       end do
    end do

  end subroutine add_random_noise


  !still assuming neutral/correct charge distribution given
  !might delete this option eventually, as now doing via kernel - unless can think of a way to make this work in direct min case?
  subroutine set_coeffs_to_tmbs()
    implicit none

    !real(kind=dp) :: sumo, sumof
    call f_zero(tmb%orbs%norb*tmb%orbs%norb, tmb%coeff(1,1))

    coeff_final=f_malloc0((/tmb%orbs%norb,tmb%orbs%norb/),id='coeff_final')

    jsforb=0
    !sumo=0.0d0
    do ifrag=1,input%frag%nfrag
       !sumof=0.0d0
       ! find reference fragment this corresponds to
       ifrag_ref=input%frag%frag_index(ifrag)
       call f_zero(tmb%orbs%norb*tmb%orbs%norb, tmb%coeff(1,1))
       nelecfrag=int(ref_frags(ifrag_ref)%nelec-input_frag_charge(ifrag))
       do jtmb=1,ref_frags(ifrag_ref)%fbasis%forbs%norb
          tmb%coeff(jsforb+jtmb,jtmb)=1.0d0
          tmb%orbs%occup(jsforb+jtmb)=real(nelecfrag,dp)/real(ref_frags(ifrag_ref)%fbasis%forbs%norb,dp) !ref_frags(ifrag_ref)%coeff(jtmb,jtmb) !
          tmb%orbs%eval(jsforb+jtmb)=-0.5d0
          !sumo=sumo+ref_frags(ifrag_ref)%coeff(jtmb,jtmb)
          !sumof=sumof+ref_frags(ifrag_ref)%coeff(jtmb,jtmb)
          !write(*,'(a,6(2x,I4),4(2x,F6.2))') 'iproc,if,ifr,ntmbf,nef,it,n,nw',iproc,ifrag,ifrag_ref,&
          !     ref_frags(ifrag_ref)%fbasis%forbs%norb,nelecfrag,jtmb,&
          !     real(nelecfrag,dp)/real(ref_frags(ifrag_ref)%fbasis%forbs%norb,dp),&
          !     ref_frags(ifrag_ref)%coeff(jtmb,jtmb),sumof,sumo
       end do
       !if (iproc==0) print*,ifrag,ifrag_ref,ref_frags(ifrag_ref)%fbasis%forbs%norb,nelecfrag,&
       !      real(nelecfrag,dp)/real(ref_frags(ifrag_ref)%fbasis%forbs%norb,dp),tmb%coeff(jtmb,jtmb)

       !ortho_size=ref_frags(ifrag_ref)%fbasis%forbs%norb
       !tmb%linmat%ovrlp_%matrix = sparsematrix_malloc_ptr(tmb%linmat%s, &
       !                           iaction=DENSE_FULL, id='tmb%linmat%ovrlp_%matrix')
       !call timing(iproc,'kernel_init','OF')
       !call uncompress_matrix2(iproc, bigdft_mpi%nproc, tmb%linmat%s, &
       !     tmb%linmat%ovrlp_%matrix_compr, tmb%linmat%ovrlp_%matrix)
       !call reorthonormalize_coeff(bigdft_mpi%iproc, bigdft_mpi%nproc, ortho_size, &
       !     tmb%orthpar%blocksize_pdsyev, tmb%orthpar%blocksize_pdgemm, input%lin%order_taylor, &
       !     tmb%orbs, tmb%linmat%s, tmb%linmat%ks, tmb%linmat%ovrlp_, tmb%coeff, ksorbs)
       !call timing(iproc,'kernel_init','ON')
       !call f_free_ptr(tmb%linmat%ovrlp_%matrix)

       ! update coeff_final matrix following coeff reorthonormalization
       do jtmb=1,ref_frags(ifrag_ref)%fbasis%forbs%norb
          ! renormalize occupations to compensate for charge leaking
          !tmb%orbs%occup(jsforb+jtmb)=tmb%orbs%occup(jsforb+jtmb)*nelecfrag/sumof
          do itmb=1,tmb%orbs%norb
             coeff_final(itmb,jsforb+jtmb)=tmb%coeff(itmb,jtmb)
          end do
       end do

       jsforb=jsforb+ref_frags(ifrag_ref)%fbasis%forbs%norb
    end do

    call vcopy(tmb%orbs%norb*tmb%orbs%norb,coeff_final(1,1),1,tmb%coeff(1,1),1)
  
    call f_free(coeff_final)

    nstates_max=tmb%orbs%norb

    if (nstates_max/=ksorbs%norb) then
       if (bigdft_mpi%iproc==0) print*,'Warning, number of states with non-zero occupation in fragments (',nstates_max,&
            ') differs from number of KS states (',ksorbs%norb,') - might have convergence problems'
    end if

  end subroutine set_coeffs_to_tmbs

  subroutine fill_occupied_coeffs()
    implicit none
  
    ! copy from coeff fragment to global coeffs - occupied states only
    isforb=0
    jsforb=0
    nstates_max=0
    eval_tmp=f_malloc0(tmb%orbs%norb,id='eval_tmp')
    coeff_final=f_malloc0((/tmb%orbs%norb,tmb%orbs%norb/),id='coeff_final')

    !initialize to zero
    if (lincombm .or. lincombp) then
       call f_zero(tmb%coeff)
       call f_zero(tmb%orbs%eval)
    end if
  
    do ifrag=1,input%frag%nfrag
       ! find reference fragment this corresponds to
       ifrag_ref=input%frag%frag_index(ifrag)
       call f_zero(tmb%orbs%norb*tmb%orbs%norb, tmb%coeff(1,1))
  
       jstate_max=(ref_frags(ifrag_ref)%nelec-input_frag_charge(ifrag))/2.0_gp+num_extra_per_frag
       if (lincombp) jstate_max=jstate_max-1
       do jtmb=1,ceiling(jstate_max)
  
          ! want random mixing across fragments in both cases
          if (random.or.completely_random) then
             do itmb=1,isforb
                call random_number(random_noise)
                random_noise=((random_noise-0.5d0)*2.0d0)*rmax
                tmb%coeff(itmb,jtmb)=random_noise
             end do
          end if
  
          do itmb=1,ref_frags(ifrag_ref)%fbasis%forbs%norb
             if (random.or.completely_random) then
                call random_number(random_noise)
                random_noise=((random_noise-0.5d0)*2.0d0)*rmax
             end if
              if (.not. completely_random) then
                tmb%coeff(isforb+itmb,jtmb)=ref_frags(ifrag_ref)%coeff(itmb,jtmb)+random_noise
             else
                tmb%coeff(isforb+itmb,jtmb)=random_noise
             end if
          end do
          tmb%orbs%eval(jsforb+jtmb)=ref_frags(ifrag_ref)%eval(jtmb)-((-1)**(ifrag))*lag_mult-homo_diff
  
          ! want random mixing across fragments in both cases
          if (random .or. completely_random) then
             do itmb=isforb+ref_frags(ifrag_ref)%fbasis%forbs%norb+1,tmb%orbs%norb
                call random_number(random_noise)
                random_noise=((random_noise-0.5d0)*2.0d0)*rmax
                tmb%coeff(itmb,jtmb)=random_noise
             end do
          end if
  
          if (ceiling(jstate_max)/=jstate_max.and.jtmb==ceiling(jstate_max)) then
             tmb%orbs%occup(jtmb+jsforb)=(jstate_max*2.0d0)-2*(ceiling(jstate_max)-1)
             ! want partly occupied states to be last of unoccupied
             eval_tmp(jsforb+jtmb)=tmb%orbs%eval(jsforb+jtmb)+10.0d0
          else
             tmb%orbs%occup(jtmb+jsforb)=2.0d0
             eval_tmp(jsforb+jtmb)=tmb%orbs%eval(jsforb+jtmb)
          end if
       end do
       nstates_max=nstates_max+ceiling(jstate_max)
  
       ! should correct the occupation for kernel here, but as we replace the smaller kernel with the correct bigger kernel
       ! don't worry about this for now
  
       ! reorthonormalize the coeffs for each fragment - don't need unoccupied states here
       ! not sure why we're not including the extra states here but for now just go with it
       ortho_size=ceiling((ref_frags(ifrag_ref)%nelec-input_frag_charge(ifrag))/2.0_gp)
       if (lincombp) ortho_size=ortho_size-1
       tmb%linmat%ovrlp_%matrix = sparsematrix_malloc_ptr(tmb%linmat%s, &
                                  iaction=DENSE_FULL, id='tmb%linmat%ovrlp_%matrix')
       call timing(iproc,'kernel_init','OF')
       call uncompress_matrix2(iproc, bigdft_mpi%nproc, tmb%linmat%s, &
            tmb%linmat%ovrlp_%matrix_compr, tmb%linmat%ovrlp_%matrix)
       call reorthonormalize_coeff(bigdft_mpi%iproc, bigdft_mpi%nproc, ortho_size, &
            tmb%orthpar%blocksize_pdsyev, tmb%orthpar%blocksize_pdgemm, input%lin%order_taylor, &
            tmb%orbs, tmb%linmat%s, tmb%linmat%ks, tmb%linmat%ovrlp_, tmb%coeff, ksorbs)
       call timing(iproc,'kernel_init','ON')
       call f_free_ptr(tmb%linmat%ovrlp_%matrix)
  
       !! debug
       !!tmb%linmat%denskern%matrix=f_malloc_ptr((/tmb%orbs%norb,tmb%orbs%norb/),id='tmb%linmat%denskern%matrix')
       !!call uncompress_matrix(bigdft_mpi%iproc,tmb%linmat%denskern)
       !!do itmb=1,tmb%orbs%norb
       !!   do jtmb=1,tmb%orbs%norb
       !!      write(30+ifrag,*) itmb,jtmb,tmb%linmat%denskern%matrix(itmb,jtmb),tmb%coeff(itmb,jtmb)
       !!   end do
       !!end do
       !!call f_free_ptr(tmb%linmat%denskern%matrix)    
       !! end debug
  
       ! update coeff_final matrix following coeff reorthonormalization
       do jtmb=1,ceiling(jstate_max)
          do itmb=1,tmb%orbs%norb
             coeff_final(itmb,jsforb+jtmb)=tmb%coeff(itmb,jtmb)
          end do
       end do
  
       isforb=isforb+ref_frags(ifrag_ref)%fbasis%forbs%norb
       jsforb=jsforb+ceiling(jstate_max)
    end do
  
    call vcopy(tmb%orbs%norb*tmb%orbs%norb,coeff_final(1,1),1,tmb%coeff(1,1),1)
  
    call f_free(coeff_final)
  
  end subroutine fill_occupied_coeffs

  subroutine fill_inbetween_coeffs()
    implicit none

    ! add states with evenly distributed charge (i.e. last few occupied/first few unoccupied)
    ! assuming each fragment has at least a LUMO and HOMO-1 available (i.e. not fully occupied)
    isforb=0
    do ifrag=1,input%frag%nfrag
       ! find reference fragment this corresponds to
       ifrag_ref=input%frag%frag_index(ifrag)
       jstate_max=(ref_frags(ifrag_ref)%nelec-input_frag_charge(ifrag))/2.0_gp+num_extra_per_frag
       if (lincombp) jstate_max=jstate_max-1

       jtmb=ceiling(jstate_max)+1
       linstate=jsforb+jtmb-ceiling(jstate_max)

       do jf=0,input%frag%nfrag-1
          fac=1.0d0/dsqrt(real(input%frag%nfrag,gp))
          if (jf==ifrag-1.and.ifrag/=1) then
             pm=-1
          else
             pm=1
          end if
          !version where occupancies are hacked instead
          !if (jf==ifrag-1) then
          !   pm=1
          !   fac=1.0d0
          !else
          !   pm=0
          !   fac=0.0d0
          !end if

          !!to probably be deleted
          !!if (jf==0) then
          !!   pm=1
          !!else
          !!   if ((jf==0.and.lincombm).or.(jf==input%frag%nfrag-1.and.lincombp)) then
          !!      pm=1
          !!      fac=1.0d0/dsqrt(real(input%frag%nfrag,gp))
          !!   else if (jf==ifrag-1) then
          !!      pm=1
          !!      fac=dsqrt(real(input%frag%nfrag,gp)-1.0d0)/dsqrt(real(input%frag%nfrag,gp))
          !!   else if ((ifrag==1.and.lincombm).or.(ifrag==input%frag%nfrag.and.lincombp)) then
          !!      !pm=-1
          !!      pm=1
          !!      fac=1.0d0/dsqrt(real(input%frag%nfrag,gp))
          !!   else
          !!      pm=0
          !!      fac=0.0d0
          !!   end if
          !!end if

          if (random .or. completely_random) then ! want random mixing across fragments in both cases
             do itmb=1,isforb
                call random_number(random_noise)
                random_noise=((random_noise-0.5d0)*2.0d0)*rmax
                tmb%coeff(itmb,jsforb+jtmb-ceiling(jstate_max)+jf)=tmb%coeff(itmb,jsforb+jtmb-ceiling(jstate_max)+jf)&
                     +random_noise*fac
             end do
          end if

          do itmb=1,ref_frags(ifrag_ref)%fbasis%forbs%norb
             if (random.or.completely_random) then
                call random_number(random_noise)
                random_noise=((random_noise-0.5d0)*2.0d0)*rmax
             end if
             if (.not. completely_random) then
                tmb%coeff(isforb+itmb,jsforb+jtmb-ceiling(jstate_max)+jf)&
                     =tmb%coeff(isforb+itmb,jsforb+jtmb-ceiling(jstate_max)+jf)&
                     +pm*(ref_frags(ifrag_ref)%coeff(itmb,jtmb)+random_noise)*fac
             else
                tmb%coeff(isforb+itmb,jsforb+jtmb-ceiling(jstate_max)+jf)&
                     =tmb%coeff(isforb+itmb,jsforb+jtmb-ceiling(jstate_max)+jf)&
                     +random_noise*fac
             end if
          end do
          tmb%orbs%eval(jsforb+jtmb-ceiling(jstate_max)+jf)=tmb%orbs%eval(jsforb+jtmb-ceiling(jstate_max)+jf)+&
               (ref_frags(ifrag_ref)%eval(jtmb))/real(input%frag%nfrag,gp)!&
               !-((-1)**(ifrag))*lag_mult-homo_diff
          eval_tmp(jsforb+jtmb-ceiling(jstate_max)+jf)=eval_tmp(jsforb+jtmb-ceiling(jstate_max)+jf)&
               +(tmb%orbs%eval(jsforb+jtmb-ceiling(jstate_max)+jf)+5.0d0)/real(input%frag%nfrag,gp)
          !print*,'eval',jsforb+jtmb-ceiling(jstate_max)+jf,eval_tmp(jsforb+jtmb-ceiling(jstate_max)+jf),&
          !     tmb%orbs%eval(jsforb+jtmb-ceiling(jstate_max)+jf)
       end do

       if (random .or. completely_random) then ! want random mixing across fragments in both cases
          do jf=0,input%frag%nfrag-1
             do itmb=isforb+ref_frags(ifrag_ref)%fbasis%forbs%norb+1,tmb%orbs%norb
                call random_number(random_noise)
                random_noise=((random_noise-0.5d0)*2.0d0)*rmax
                tmb%coeff(itmb,jsforb+jtmb-ceiling(jstate_max)+jf)=tmb%coeff(itmb,jsforb+jtmb-ceiling(jstate_max)+jf)&
                     +random_noise*fac
             end do
          end do
       end if

       !GENERALIZE HERE FOR OTHER CHARGE STATES
       !eval_tmp(jsforb+jtmb-ceiling(jstate_max)+jf)=eval_tmp(jsforb+jtmb-ceiling(jstate_max)+jf)+5.0d0/real(input%frag%nfrag,gp)
       tmb%orbs%occup(jsforb+jtmb-ceiling(jstate_max)+ifrag-1)=0.0d0
       !version where occupancies are fixed
       !if (lincombp) tmb%orbs%occup(jsforb+jtmb-ceiling(jstate_max)+ifrag)=real(2*input%frag%nfrag-1,gp)/real(input%frag%nfrag,gp)

       !if (bigdft_mpi%iproc==0) print*,'ifrag,jtmb,occ,iorb',ifrag,jtmb,0.0,jsforb+jtmb-ceiling(jstate_max)
       !end do

       isforb=isforb+ref_frags(ifrag_ref)%fbasis%forbs%norb
       !jsforb=jsforb+ref_frags(ifrag_ref)%fbasis%forbs%norb-ceiling(jstate_max)
    end do

    !first/last should be occupied with occupancy 1 - needs generalizing
    if (lincombm) then
       tmb%orbs%occup(jsforb+jtmb-ceiling(jstate_max))=1.0d0
    else if (lincombp) then
       tmb%orbs%occup(jsforb+jtmb-ceiling(jstate_max)+input%frag%nfrag-1)=1.0d0
    end if

    jsforb=jsforb+input%frag%nfrag

  end subroutine fill_inbetween_coeffs

  subroutine fill_unoccupied_coeffs()
    implicit none

    isforb=0
    do ifrag=1,input%frag%nfrag
       ! find reference fragment this corresponds to
       ifrag_ref=input%frag%frag_index(ifrag)
       jstate_max=(ref_frags(ifrag_ref)%nelec-input_frag_charge(ifrag))/2.0_gp+num_extra_per_frag
       if (lincombm) jstate_max=jstate_max+1
       do jtmb=ceiling(jstate_max)+1,ref_frags(ifrag_ref)%fbasis%forbs%norb
          if (random .or. completely_random) then ! want random mixing across fragments in both cases
             do itmb=1,isforb
                call random_number(random_noise)
                random_noise=((random_noise-0.5d0)*2.0d0)*rmax
                tmb%coeff(itmb,jsforb+jtmb-ceiling(jstate_max))=random_noise
             end do
          end if

          do itmb=1,ref_frags(ifrag_ref)%fbasis%forbs%norb
             if (random.or.completely_random) then
                call random_number(random_noise)
                random_noise=((random_noise-0.5d0)*2.0d0)*rmax
             end if
             if (.not. completely_random) then
                tmb%coeff(isforb+itmb,jsforb+jtmb-ceiling(jstate_max))=ref_frags(ifrag_ref)%coeff(itmb,jtmb)+random_noise
             else
                tmb%coeff(isforb+itmb,jsforb+jtmb-ceiling(jstate_max))=random_noise
             end if
          end do
          tmb%orbs%eval(jsforb+jtmb-ceiling(jstate_max))=ref_frags(ifrag_ref)%eval(jtmb)-((-1)**(ifrag))*lag_mult-homo_diff
          eval_tmp(jsforb+jtmb-ceiling(jstate_max))=tmb%orbs%eval(jsforb+jtmb-ceiling(jstate_max))+20.0d0

          if (random .or. completely_random) then ! want random mixing across fragments in both cases
             do itmb=isforb+ref_frags(ifrag_ref)%fbasis%forbs%norb+1,tmb%orbs%norb
                call random_number(random_noise)
                random_noise=((random_noise-0.5d0)*2.0d0)*rmax
                tmb%coeff(itmb,jsforb+jtmb-ceiling(jstate_max))=random_noise
             end do
          end if

          tmb%orbs%occup(jsforb+jtmb-ceiling(jstate_max))=0.0d0
          !if (bigdft_mpi%iproc==0) print*,'ifrag,jtmb,occ,iorb',ifrag,jtmb,0.0,jsforb+jtmb-ceiling(jstate_max)
       end do

       isforb=isforb+ref_frags(ifrag_ref)%fbasis%forbs%norb
       jsforb=jsforb+ref_frags(ifrag_ref)%fbasis%forbs%norb-ceiling(jstate_max)
    end do

  end subroutine fill_unoccupied_coeffs

  subroutine reorder_and_print_coeffs()
    use module_interfaces, only: write_eigenvalues_data
    implicit none

    !!print*,'nstates_max:',nstates_max,ksorbs%norb,tmb%orbs%norb
    if (bigdft_mpi%iproc==0) then
        call yaml_map('nstates_max',nstates_max)
        call yaml_map('ksorbs%norb',ksorbs%norb)
        call yaml_map('tmb%orbs%norb',tmb%orbs%norb)
    end if

    ! reorder unoccupied states so that extra states functions correctly
    ! still needed just in case number of empty bands doesn't divide by number of fragments
    ipiv=f_malloc(tmb%orbs%norb,id='ipiv')
    call order_coeffs_by_energy(tmb%orbs%norb-nstates_max,tmb%orbs%norb,tmb%coeff(1,nstates_max+1),&
         eval_tmp(nstates_max+1),ipiv(1))!,tmb%orbs%eval(nstates_max+1))
    eval_tmp2=f_malloc(tmb%orbs%norb-nstates_max,id='eval_tmp2')
    call vcopy(tmb%orbs%norb-nstates_max,tmb%orbs%eval(nstates_max+1),1,eval_tmp2(1),1)
    do itmb=nstates_max+1,tmb%orbs%norb
       tmb%orbs%eval(itmb)=eval_tmp2(ipiv(itmb-nstates_max))
    end do
    call vcopy(tmb%orbs%norb-nstates_max,tmb%orbs%occup(nstates_max+1),1,eval_tmp2(1),1)
    do itmb=nstates_max+1,tmb%orbs%norb
       tmb%orbs%occup(itmb)=eval_tmp2(ipiv(itmb-nstates_max))
    end do
    call f_free(eval_tmp2)
    ! reorder ksorbs%norb states by energy - no longer taking charge as input
    call order_coeffs_by_energy(ksorbs%norb,tmb%orbs%norb,tmb%coeff(1,1),eval_tmp(1),ipiv(1))!,tmb%orbs%eval(1))
            !eval_tmp2=f_malloc(tmb%orbs%norb,id='eval_tmp2')
            !call vcopy(tmb%orbs%norb,tmb%orbs%occup(1),1,eval_tmp2(1),1)
            !do itmb=1,ksorbs%norb
            !   tmb%orbs%occup(itmb)=eval_tmp2(ipiv(itmb))
            !end do
            !call vcopy(tmb%orbs%norb,tmb%orbs%eval(1),1,eval_tmp2(1),1)
            !call vcopy(tmb%orbs%norb,eval_tmp(1),1,tmb%orbs%eval(1),1)
            !nullify(mom_vec_fake)
            !if (bigdft_mpi%iproc==0) then 
            !   call write_eigenvalues_data(0.1d0,tmb%orbs,mom_vec_fake)
            !end if
            !call vcopy(tmb%orbs%norb,eval_tmp2(1),1,tmb%orbs%eval(1),1)
            !call f_free(eval_tmp2)
    call vcopy(ksorbs%norb,tmb%orbs%eval(1),1,eval_tmp(1),1)
    do itmb=1,ksorbs%norb
       tmb%orbs%eval(itmb)=eval_tmp(ipiv(itmb))
    end do
    call f_free(ipiv)
    ! debug
    !do itmb=1,tmb%orbs%norb
    !   do jtmb=1,tmb%orbs%norb
    !      write(50,*) itmb,jtmb,tmb%coeff(itmb,jtmb)
    !   end do
    !end do
    ! end debug
    ! print starting eigenvalues
    if(bigdft_mpi%iproc==0) then
       call yaml_sequence_open('TMB eigenvalues',flow=.true.)
       call yaml_newline()
       do iorb=1,tmb%orbs%norb
           call yaml_mapping_open(flow=.true.)
           call yaml_map('index',iorb)
           call yaml_map('value',tmb%orbs%eval(iorb),fmt='(es20.12)')
           call yaml_mapping_close()
           if(iorb==ksorbs%norb) then
               call yaml_comment('  <-- last occupied orbital')
           else if(iorb==ksorbs%norb+1) then
               call yaml_comment('  <-- first virtual orbital')
           else
           end if
           call yaml_newline()
       end do
       call yaml_sequence_close()
    end if

    if (nstates_max/=ksorbs%norb) then
       if (bigdft_mpi%iproc==0) print*,'Warning, number of states with non-zero occupation in fragments (',nstates_max,&
            ') differs from number of KS states (',ksorbs%norb,') - might have convergence problems'
    end if

    ! shouldn't we actually check if nstates_max is (less than?) equal to ks_e?
    ! don't think it needs to be an output still...

    !!!!!!!!!!!!!!!
    ! need the eigenvalues to be in ksorbs%eval
    call vcopy(ksorbs%norb,tmb%orbs%eval(1),1,ksorbs%eval(1),1)
    call evaltoocc(bigdft_mpi%iproc,bigdft_mpi%nproc,.false.,input%tel,ksorbs,input%occopt)

    nullify(mom_vec_fake)
    if (bigdft_mpi%iproc ==0) then 
       call write_eigenvalues_data(0.1d0,ksorbs,mom_vec_fake)
    end if
    !!!!!!!!!!!!!!!
  end subroutine reorder_and_print_coeffs

end subroutine fragment_coeffs_to_kernel





!think about cdft and charged systems...
!also still need to activate completely random case, but need to think about purification first as will definitely be necessary
subroutine fragment_kernels_to_kernel(iproc,input,input_frag_charge,ref_frags,tmb,ksorbs,&
  overlap_calculated,cdft,diagonal_kernel,max_nbasis_env,frag_env_mapping,rmax)
  use yaml_output
  use module_base
  use module_types
  use module_fragments
  use communications_base, only: TRANSPOSE_FULL
  use communications, only: transpose_localized
  use sparsematrix_base, only: sparsematrix_malloc0_ptr, DENSE_FULL, assignment(=)
  use sparsematrix, only: uncompress_matrix, gather_matrix_from_taskgroups_inplace, &
                          uncompress_matrix2, compress_matrix
  use transposed_operations, only: calculate_overlap_transposed
  implicit none
  type(DFT_wavefunction), intent(inout) :: tmb
  type(input_variables), intent(in) :: input
  type(system_fragment), dimension(input%frag%nfrag_ref), intent(inout) :: ref_frags
  type(orbitals_data), intent(inout) :: ksorbs
  logical, intent(inout) :: overlap_calculated
  real(kind=gp), dimension(input%frag%nfrag), intent(in) :: input_frag_charge
  integer, intent(in) :: iproc
  logical, intent(in) :: cdft
  logical, intent(in) :: diagonal_kernel
  integer, intent(in) :: max_nbasis_env
  integer, dimension(input%frag%nfrag,max_nbasis_env,3), intent(in) :: frag_env_mapping
  real(kind=8), intent(in) :: rmax

  integer :: iorb, isforb, jsforb, ifrag, ifrag_ref, itmb, jtmb, num_extra_per_frag, linstate, jf, pm, ortho_size, s
  integer, allocatable, dimension(:) :: ipiv
  real(gp), dimension(:,:), allocatable :: coeff_final
  real(gp) :: nelecorbs, nelecfrag_tot, jstate_max, homo_diff, lag_mult, fac, nelecfrag
  real(gp), dimension(:), allocatable :: eval_tmp, eval_tmp2
  character(len=*), parameter :: subname='fragment_coeffs_to_kernel'

  integer :: rand_size
  integer, allocatable, dimension(:) :: rand_seed
  real(kind=dp) :: rtime, random_noise
  character(len=10) :: sys_time
  logical :: random, completely_random, env_exists, neutral

  real(wp), dimension(:,:,:), pointer :: mom_vec_fake

  call timing(iproc,'kernel_init','ON')
  call f_routine(id='fragment_coeffs_to_kernel')

  !! need to do this properly/rearrange routines
  !if (cdft) then
  !   ! otherwise doesn't make sense
  !   if (input%frag%nfrag_ref==2) homo_diff=(ref_frags(1)%eval(ceiling(ref_frags(1)%nelec/2.0_gp))&
  !        -ref_frags(2)%eval(ceiling(ref_frags(2)%nelec/2.0_gp)))/2.0d0
  !   !if (cdft%charge<0) lag_mult=-0.5, otherwise +0.5
  !   lag_mult=-0.05d0
  !else
  !   homo_diff=0.0d0
  !   lag_mult=0.0d0
  !end if

  ! adding random noise to kernel to help with local minima problem
  random=(rmax>0.0d0) ! add a bit of noise
  completely_random=.false. ! completely random start for coeffs

  random_noise=0.0d0
  rtime=0.0d0
  if (random .or. completely_random) then
     call random_seed(size=rand_size)
     allocate(rand_seed(1:rand_size))
     call date_and_time(time=sys_time)
     ! coeffs need to be the same across processors
     if (iproc==0) read(sys_time,*) rtime
     if (bigdft_mpi%nproc>1) call mpiallred(rtime, 1, mpi_sum, comm=bigdft_mpi%mpi_comm)
     rand_seed=int(rtime*1000.0_dp)
     call random_seed(put=rand_seed)
     deallocate(rand_seed) 
  end if

  !can we add charge to kernel for cdft?
  !nelecfrag_tot=0
  !do ifrag=1,input%frag%nfrag
  !   ifrag_ref=input%frag%frag_index(ifrag)
  !   nelecfrag_tot=nelecfrag_tot+ref_frags(ifrag_ref)%nelec-input_frag_charge(ifrag)
  !end do

  if (completely_random) then
     if (bigdft_mpi%iproc==0) print*,'Starting kernel is replaced with a random guess'
  else if (random) then
     if (bigdft_mpi%iproc==0) print*,'Random noise added to starting kernel'
  end if

  !!in theory we could add/remove states depending on their energies, but for now we force the user to specify
  !!need to include occupations as we actually want to compare number of electrons here?
  !nelecorbs=0
  !do iorb=1,ksorbs%norb
  !   nelecorbs=nelecorbs+ksorbs%occup(iorb)
  !end do

  !for the moment just carry on even if charge is incorrect
  !if (nint(nelecorbs)/=nelecfrag_tot) then
  !   print*,'User should specify which fragments charges are added to/removed from in charged fragment calculation',&
  !        nelecfrag_tot,nelecorbs,ksorbs%norb
  !end if

  !likewise, ignore extra states
  !if (mod(input%norbsempty,input%frag%nfrag)/=0) then
  !   if (bigdft_mpi%iproc==0) call yaml_warning('Number of extra bands does not divide evenly among fragments')
  !   num_extra_per_frag=(input%norbsempty-mod(input%norbsempty,input%frag%nfrag))/input%frag%nfrag
  !else
  !   num_extra_per_frag=input%norbsempty/input%frag%nfrag
  !end if

  ! Calculate the overlap matrix between the TMBs - why? - I guess for purification/reorthonormalization...
  !if(.not. overlap_calculated) then
  !   call timing(iproc,'kernel_init','OF')
  !   if(.not.tmb%can_use_transposed) then
  !       if(associated(tmb%psit_c)) then
  !           call f_free_ptr(tmb%psit_c)
  !       end if
  !       if(associated(tmb%psit_f)) then
  !           call f_free_ptr(tmb%psit_f)
  !       end if
  !       tmb%psit_c = f_malloc_ptr(sum(tmb%collcom%nrecvcounts_c),id='tmb%psit_c')
  !       tmb%psit_f = f_malloc_ptr(7*sum(tmb%collcom%nrecvcounts_f),id='tmb%psit_f')
  !       call transpose_localized(bigdft_mpi%iproc, bigdft_mpi%nproc, tmb%npsidim_orbs, tmb%orbs, tmb%collcom, &
  !            TRANSPOSE_FULL, tmb%psi, tmb%psit_c, tmb%psit_f, tmb%lzd)
  !       tmb%can_use_transposed=.true.
  !   end if
  !
  !   call calculate_overlap_transposed(bigdft_mpi%iproc, bigdft_mpi%nproc, tmb%orbs, tmb%collcom, &
  !        tmb%psit_c, tmb%psit_c, tmb%psit_f, tmb%psit_f, tmb%linmat%s, tmb%linmat%ovrlp_)
  !   !!call gather_matrix_from_taskgroups_inplace(bigdft_mpi%iproc, bigdft_mpi%nproc, tmb%linmat%s, tmb%linmat%ovrlp_)
  !   ! This can then be deleted if the transition to the new type has been completed.
  !   !tmb%linmat%ovrlp%matrix_compr=tmb%linmat%ovrlp_%matrix_compr
  !
  !   overlap_calculated=.true.
  !   call timing(iproc,'kernel_init','ON')
  !end if

  ! set coeffs to tmbs with occupancies so each is equally weighted
  ! can only easily be used to get kernel (via FOE) or using diag
  ! occupancies should be reset afterwards

  ! for now working in dense format
  tmb%linmat%kernel_%matrix = sparsematrix_malloc0_ptr(tmb%linmat%l, DENSE_FULL, id='tmb%linmat%kernel__%matrix')

  ! check if environment exists
  env_exists=.false.
  do ifrag_ref=1,input%frag%nfrag_ref
     if (ref_frags(ifrag_ref)%astruct_env%nat/=0) then
        env_exists=.true.
        exit
     end if
  end do

  if ((.not. diagonal_kernel) .or. completely_random) then
     if (env_exists) then
        call fill_kernel_from_frag_env()
     else
        call fill_kernel_from_fragments()
     end if

     !it seems to work better to just start with the neutral guess, at least for constraining charge differences
     !should check if adding purification improves the situation
     !distribute extra charge across tmbs if necessary
     !neutral=.true.
     !do ifrag=1,input%frag%nfrag
     !   if (input_frag_charge(ifrag)/=0) then
     !      neutral=.false.
     !      exit
     !   end if
     !end do

     !if (.not. neutral) call add_charge_to_diagonal()

  else if (completely_random) then
      call fill_random_kernel()
  else if (diagonal_kernel) then
     call fill_diagonal_kernel()
  end if

  if (random .and. (.not. completely_random)) then
     call add_noise_to_kernel()
  end if


  !if (iproc==0) then
  !    open(27)
  !    do itmb=1,tmb%orbs%norb
  !      do jtmb=1,tmb%orbs%norb
  !         write(27,*) itmb,jtmb,tmb%linmat%kernel_%matrix(itmb,jtmb,1)
  !      end do
  !    end do
  !   write(27,*) ''
  !   close(27)
  !end if 


  call compress_matrix(iproc,tmb%linmat%l,inmat=tmb%linmat%kernel_%matrix,outmat=tmb%linmat%kernel_%matrix_compr)  
  call f_free_ptr(tmb%linmat%kernel_%matrix) 

  call f_release_routine()
  call timing(iproc,'kernel_init','OF')

contains

  !add excess/deficit of electrons to diagonal elements (for fragment kernel approach) - assume kernels being read in were from neutral systems!
  subroutine add_charge_to_diagonal()
    implicit none

    isforb=0
    do ifrag=1,input%frag%nfrag
       ifrag_ref=input%frag%frag_index(ifrag)

       !nelecfrag=ref_frags(ifrag_ref)%nelec-input_frag_charge(ifrag)
       do itmb=1,ref_frags(ifrag_ref)%fbasis%forbs%norb
          tmb%linmat%kernel_%matrix(isforb+itmb,isforb+itmb,1) = tmb%linmat%kernel_%matrix(isforb+itmb,isforb+itmb,1) &
               - input_frag_charge(ifrag)/real(ref_frags(ifrag_ref)%fbasis%forbs%norb,dp)
       end do

       isforb=isforb+ref_frags(ifrag_ref)%fbasis%forbs%norb
    end do

  end subroutine add_charge_to_diagonal


  !assuming neutral/correct charge distribution given
  subroutine fill_diagonal_kernel()
    implicit none

    isforb=0
    do ifrag=1,input%frag%nfrag
       ifrag_ref=input%frag%frag_index(ifrag)

       nelecfrag=ref_frags(ifrag_ref)%nelec-input_frag_charge(ifrag)
       do itmb=1,ref_frags(ifrag_ref)%fbasis%forbs%norb
          tmb%linmat%kernel_%matrix(isforb+itmb,isforb+itmb,1)=nelecfrag/real(ref_frags(ifrag_ref)%fbasis%forbs%norb,dp)
       end do

       isforb=isforb+ref_frags(ifrag_ref)%fbasis%forbs%norb
    end do

  end subroutine fill_diagonal_kernel

  subroutine fill_kernel_from_frag_env()
    implicit none

    integer :: itmb_full, jtmb_full, ntmb
    logical, parameter :: nn_only=.false.

    ! add loop over spin

    do ifrag=1,input%frag%nfrag
       ifrag_ref=input%frag%frag_index(ifrag)

       ! use frag_env_mapping to figure out where to put elements
       ! assuming complete symmetry between frag+env i.e. overwriting elements rather than averaging

       ! include frag-environment only
       if (nn_only) then
          ntmb=ref_frags(ifrag_ref)%fbasis%forbs%norb
       ! also include env-env terms
       else
          ntmb=ref_frags(ifrag_ref)%nbasis_env
       end if

       do itmb=1,ntmb
          itmb_full = frag_env_mapping(ifrag,itmb,1)

          do jtmb=itmb,ref_frags(ifrag_ref)%nbasis_env
             jtmb_full = frag_env_mapping(ifrag,jtmb,1)

             tmb%linmat%kernel_%matrix(itmb_full,jtmb_full,1) = ref_frags(ifrag_ref)%kernel_env(itmb,jtmb,1)
             tmb%linmat%kernel_%matrix(jtmb_full,itmb_full,1) = ref_frags(ifrag_ref)%kernel_env(itmb,jtmb,1)
             !write(*,'(A,3(2(1x,I4),2x),F12.6)') 'Kij',ifrag,ifrag_ref,itmb,jtmb,itmb_full,jtmb_full,ref_frags(ifrag_ref)%kernel_env(itmb,jtmb,1)
          end do
       end do
    end do
  
    ! should we purify kernel?
  
  end subroutine fill_kernel_from_frag_env

  subroutine fill_kernel_from_fragments()
    implicit none

    ! add loop over spin

    isforb=0
    do ifrag=1,input%frag%nfrag
       ifrag_ref=input%frag%frag_index(ifrag)
       do itmb=1,ref_frags(ifrag_ref)%fbasis%forbs%norb
          call vcopy(ref_frags(ifrag_ref)%fbasis%forbs%norb,ref_frags(ifrag_ref)%kernel(1,itmb,1),1,&
               tmb%linmat%kernel_%matrix(1+isforb,itmb+isforb,1),1)
       end do
  
       isforb=isforb+ref_frags(ifrag_ref)%fbasis%forbs%norb
    end do
  
    ! should we purify kernel?
  
  end subroutine fill_kernel_from_fragments

  subroutine fill_random_kernel()
    implicit none

    do itmb=1,tmb%linmat%l%nfvctr
       do jtmb=1,tmb%linmat%l%nfvctr
          call random_number(random_noise)
          random_noise=((random_noise-0.5d0)*2.0d0)*rmax
          tmb%linmat%kernel_%matrix(itmb,jtmb,1)=random_noise
       end do
    end do

    ! almost certainly need purification here - add in all cases?

  end subroutine fill_random_kernel

  subroutine add_noise_to_kernel()
    implicit none

    do itmb=1,tmb%linmat%l%nfvctr
       do jtmb=1,tmb%linmat%l%nfvctr
          call random_number(random_noise)
          random_noise=((random_noise-0.5d0)*2.0d0)*rmax
          tmb%linmat%kernel_%matrix(itmb,jtmb,1)=tmb%linmat%kernel_%matrix(itmb,jtmb,1)+random_noise
       end do
    end do

  end subroutine add_noise_to_kernel


end subroutine fragment_kernels_to_kernel



