subroutine fragment_coeffs_to_kernel(iproc,input,input_frag_charge,ref_frags,tmb,ksorbs,overlap_calculated,&
  nstates_max,cdft)
  use yaml_output
  use module_base
  use module_types
  use module_interfaces, except_this_one => fragment_coeffs_to_kernel
  use module_fragments
  use communications_base, only: TRANSPOSE_FULL
  use communications, only: transpose_localized
  use sparsematrix_base, only: sparsematrix_malloc_ptr, DENSE_FULL, assignment(=)
  use sparsematrix, only: uncompress_matrix, gather_matrix_from_taskgroups_inplace, &
                          uncompress_matrix2
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

  integer :: iorb, isforb, jsforb, ifrag, ifrag_ref, itmb, jtmb, num_extra_per_frag
  integer, allocatable, dimension(:) :: ipiv
  real(gp), dimension(:,:), allocatable :: coeff_final
  !real(gp), dimension(:,:), allocatable :: ks, ksk
  !*real(gp), dimension(:), allocatable :: kernel_final
  real(gp) :: nelecorbs, nelecfrag_tot, jstate_max, homo_diff, lag_mult
  real(gp), dimension(:), allocatable :: eval_tmp, eval_tmp2
  character(len=*), parameter :: subname='fragment_coeffs_to_kernel'

  integer :: rand_size
  integer, allocatable, dimension(:) :: rand_seed
  real(kind=dp) :: rtime, random_noise, rmax
  character(len=10) :: sys_time
  logical :: random, completely_random

  real(wp), dimension(:,:,:), pointer :: mom_vec_fake

  call timing(iproc,'kernel_init','ON')
  call f_routine(id='fragment_coeffs_to_kernel')

  ! need to do this properly/rearrange routines
  if (cdft) then
     ! otherwise doesn't make sense
     if (input%frag%nfrag_ref==2) homo_diff=(ref_frags(1)%eval(ceiling(ref_frags(1)%nelec/2.0_gp))&
          -ref_frags(2)%eval(ceiling(ref_frags(2)%nelec/2.0_gp)))/2.0d0
     !if (cdft%charge<0) lag_mult=-0.5, otherwise +0.5
     lag_mult=-0.05d0
  else
     homo_diff=0.0d0
     lag_mult=0.0d0
  end if

  ! adding random noise to starting to help with local minima problem
  random=.false. ! add a bit of noise
  completely_random=.false. ! completely random start for coeffs

  rmax=0.2d0
  random_noise=0.0d0
  rtime=0.0d0
  if (random .or. completely_random) then
     call random_seed(size=rand_size)
     allocate(rand_seed(1:rand_size))
     call date_and_time(time=sys_time)
     ! coeffs need to be the same across processors
     if (iproc==0) read(sys_time,*) rtime
     if (bigdft_mpi%nproc>1) call mpiallred(rtime, 1, mpi_sum, bigdft_mpi%mpi_comm)
     rand_seed=int(rtime*1000.0_dp)
     call random_seed(put=rand_seed)
     deallocate(rand_seed) 
  end if
  nstates_max=0
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

  if (nint(nelecorbs)/=nelecfrag_tot) then
     print*,'User must specify which fragments charges are added to/removed from in charged fragment calculation',&
          nelecfrag_tot,nelecorbs,ksorbs%norb
     stop
  end if

  if (mod(input%norbsempty,input%frag%nfrag)/=0) then
     if (bigdft_mpi%iproc==0) call yaml_warning('Number of extra bands does not divide evenly among fragments')
     !if (bigdft_mpi%iproc==0) print*,'Warning, number of extra bands does not divide evenly among fragments'
     num_extra_per_frag=(input%norbsempty-mod(input%norbsempty,input%frag%nfrag))/input%frag%nfrag
  else
     num_extra_per_frag=input%norbsempty/input%frag%nfrag
  end if


  eval_tmp=f_malloc(tmb%orbs%norb,id='eval_tmp')
  coeff_final=f_malloc((/tmb%orbs%norb,tmb%orbs%norb/),id='coeff_final')
  !*kernel_final=f_malloc(tmb%linmat%denskern%nvctr,id='kernel_final')
  !ref_frags(ifrag_ref)%kernel=f_malloc((/tmb%orbs%norb,tmb%orbs%norb/),id='ref_frags(ifrag_ref)%kernel')

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
     !call timing(iproc,'renormCoefComp','OF')

     call calculate_overlap_transposed(bigdft_mpi%iproc, bigdft_mpi%nproc, tmb%orbs, tmb%collcom, &
          tmb%psit_c, tmb%psit_c, tmb%psit_f, tmb%psit_f, tmb%linmat%s, tmb%linmat%ovrlp_)
     !!call gather_matrix_from_taskgroups_inplace(bigdft_mpi%iproc, bigdft_mpi%nproc, tmb%linmat%s, tmb%linmat%ovrlp_)
     ! This can then be deleted if the transition to the new type has been completed.
     !tmb%linmat%ovrlp%matrix_compr=tmb%linmat%ovrlp_%matrix_compr


     !call timing(iproc,'renormCoefComp','ON')
     overlap_calculated=.true.
     call timing(iproc,'kernel_init','ON')
  end if

  ! copy from coeff fragment to global coeffs - occupied states only
  isforb=0
  jsforb=0
  call f_zero(coeff_final)
  !*call f_zero(tmb%linmat%denskern%nvctr,kernel_final(1))
  !!tmb%linmat%ovrlp%matrix=f_malloc_ptr((/tmb%orbs%norb,tmb%orbs%norb/),id='tmb%ovrlp%matrix')
  !!call uncompress_matrix(iproc,tmb%linmat%ovrlp)
  do ifrag=1,input%frag%nfrag
     ! find reference fragment this corresponds to
     ifrag_ref=input%frag%frag_index(ifrag)
     call f_zero(tmb%orbs%norb*tmb%orbs%norb, tmb%coeff(1,1))

     jstate_max=(ref_frags(ifrag_ref)%nelec-input_frag_charge(ifrag))/2.0_gp+num_extra_per_frag
     !jstate_max=ref_frags(ifrag_ref)%nelec/2.0_gp+num_extra_per_frag
     do jtmb=1,ceiling(jstate_max)

        if (random .or. completely_random) then ! want random mixing across fragments in both cases
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
           tmb%orbs%eval(jsforb+jtmb)=ref_frags(ifrag_ref)%eval(jtmb)-((-1)**(ifrag))*lag_mult-homo_diff
        end do

        if (random .or. completely_random) then ! want random mixing across fragments in both cases
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
        !if (bigdft_mpi%iproc==0) print*,'ifrag,jtmb,occ,iorb',ifrag,jtmb,tmb%orbs%occup(jtmb+jsforb),jtmb+jsforb
     end do
     nstates_max=nstates_max+ceiling(jstate_max)

     ! debug
     !do itmb=1,tmb%orbs%norb
     !   do jtmb=1,tmb%orbs%norb
     !      write(40+ifrag,*) itmb,jtmb,tmb%coeff(itmb,jtmb)
     !  end do
     !end do
     ! end debug

     !call f_zero(tmb%linmat%denskern%nvctr,tmb%linmat%denskern%matrix_compr(1))

     ! should correct the occupation for kernel here, but as we replace the smaller kernel with the correct bigger kernel
     ! don't worry about this for now

     ! reorthonormalize the coeffs for each fragment - don't need unoccupied states here
     tmb%linmat%ovrlp_%matrix = sparsematrix_malloc_ptr(tmb%linmat%s, &
                                iaction=DENSE_FULL, id='tmb%linmat%ovrlp_%matrix')
     call timing(iproc,'kernel_init','OF')
     call uncompress_matrix2(iproc, bigdft_mpi%nproc, tmb%linmat%s, &
          tmb%linmat%ovrlp_%matrix_compr, tmb%linmat%ovrlp_%matrix)
     call reorthonormalize_coeff(bigdft_mpi%iproc, bigdft_mpi%nproc, &
          ceiling((ref_frags(ifrag_ref)%nelec-input_frag_charge(ifrag))/2.0_gp), &
          tmb%orthpar%blocksize_pdsyev, tmb%orthpar%blocksize_pdgemm, input%lin%order_taylor, &
          tmb%orbs, tmb%linmat%s, tmb%linmat%ks, tmb%linmat%ovrlp_, tmb%coeff, ksorbs)
     call timing(iproc,'kernel_init','ON')
     call f_free_ptr(tmb%linmat%ovrlp_%matrix)

     !! debug
     !!output final kernel
     !! 20 - if just calculate, 21 if reconstruct total, 22 if reconstruct then sum
     !tmb%linmat%denskern%matrix=f_malloc_ptr((/tmb%orbs%norb,tmb%orbs%norb/),id='tmb%linmat%denskern%matrix')
     !call uncompress_matrix(bigdft_mpi%iproc,tmb%linmat%denskern)
     !!do itmb=1,tmb%orbs%norb
     !!   do jtmb=1,tmb%orbs%norb
     !!      write(30+ifrag,*) itmb,jtmb,tmb%linmat%denskern%matrix(itmb,jtmb),tmb%coeff(itmb,jtmb)
     !!   end do
     !!end do
     !do itmb=1,ref_frags(ifrag_ref)%fbasis%forbs%norb
     !   do jtmb=1,ref_frags(ifrag_ref)%fbasis%forbs%norb
     !      ref_frags(ifrag_ref)%kernel(itmb,jtmb)=tmb%linmat%denskern%matrix(itmb,jtmb).....
     !   end do
     !end do
     !call f_free_ptr(tmb%linmat%denskern%matrix)    
     !! end debug

     ! assemble complete kernel from separate fragment kernels
     !call daxpy(tmb%linmat%denskern%nvctr,1.0d0,tmb%linmat%denskern%matrix_compr(1),1,kernel_final(1),1)

     ! update coeff_final matrix following coeff reorthonormalization
     do jtmb=1,ceiling(jstate_max)
        do itmb=1,tmb%orbs%norb
           coeff_final(itmb,jsforb+jtmb)=tmb%coeff(itmb,jtmb)
        end do
     end do

     isforb=isforb+ref_frags(ifrag_ref)%fbasis%forbs%norb
     jsforb=jsforb+ceiling(jstate_max)
  end do
  !!call f_free_ptr(tmb%linmat%ovrlp%matrix)

  !*call vcopy(tmb%linmat%denskern%nvctr,kernel_final(1),1,tmb%linmat%denskern%matrix_compr(1),1)
  call vcopy(tmb%orbs%norb*tmb%orbs%norb,coeff_final(1,1),1,tmb%coeff(1,1),1)

  !*call f_free(kernel_final)
  call f_free(coeff_final)

  ! debug
  !output final kernel
  ! 20 - if just calculate, 21 if reconstruct total, 22 if reconstruct then sum
  !tmb%linmat%denskern%matrix=f_malloc_ptr((/tmb%orbs%norb,tmb%orbs%norb/),id='tmb%linmat%denskern%matrix')
  !call uncompress_matrix(bigdft_mpi%iproc,tmb%linmat%denskern)
  !do itmb=1,tmb%orbs%norb
  !   do jtmb=1,tmb%orbs%norb
  !      write(22,*) itmb,jtmb,tmb%linmat%denskern%matrix(itmb,jtmb),tmb%coeff(itmb,jtmb)
  !   end do
  !end do

  ! check final kernel is idempotent
  !tmb%linmat%ovrlp%matrix=f_malloc_ptr((/tmb%orbs%norb,tmb%orbs%norb/),id='tmb%linmat%ovrlp%matrix')
  !ks=f_malloc((/tmb%orbs%norb,tmb%orbs%norb/),id='ks')
  !ksk=f_malloc((/tmb%orbs%norb,tmb%orbs%norb/),id='ksk')
  !call uncompress_matrix(bigdft_mpi%iproc,tmb%linmat%ovrlp)
  !call dgemm('n', 't', tmb%orbs%norb, tmb%orbs%norb, tmb%orbs%norb, 1.d0, tmb%linmat%denskern%matrix(1,1), tmb%orbs%norb, &
  !           tmb%linmat%ovrlp%matrix(1,1), tmb%orbs%norb, 0.d0, ks(1,1), tmb%orbs%norb) 
  !call dgemm('n', 't', tmb%orbs%norb, tmb%orbs%norb, tmb%orbs%norb, 1.d0, ks(1,1), tmb%orbs%norb, &
  !           tmb%linmat%denskern%matrix(1,1), tmb%orbs%norb, 0.d0, ksk(1,1), tmb%orbs%norb)

  !nonidem=0
  !do itmb=1,tmb%orbs%norb
  !   do jtmb=1,tmb%orbs%norb
  !      write(60,*) itmb,jtmb,tmb%linmat%denskern%matrix(itmb,jtmb),ksk(itmb,jtmb),&
  !           tmb%linmat%denskern%matrix(itmb,jtmb)-ksk(itmb,jtmb)
  !      nonidem=nonidem+tmb%linmat%denskern%matrix(itmb,jtmb)-ksk(itmb,jtmb)
  !   end do
  !end do
  !print*,'non idempotency',nonidem/tmb%orbs%norb**2

  !call f_free(ks) 
  !call f_free(ksk) 
  !call f_free_ptr(tmb%linmat%ovrlp%matrix)   
  !call f_free_ptr(tmb%linmat%denskern%matrix)    
  ! end debug

  ! add unoccupied states to complete coeffs
  isforb=0
  do ifrag=1,input%frag%nfrag
     ! find reference fragment this corresponds to
     ifrag_ref=input%frag%frag_index(ifrag)
     jstate_max=(ref_frags(ifrag_ref)%nelec-input_frag_charge(ifrag))/2.0_gp+num_extra_per_frag
     !jstate_max=ref_frags(ifrag_ref)%nelec/2.0_gp+num_extra_per_frag
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
           tmb%orbs%eval(jsforb+jtmb-ceiling(jstate_max))=ref_frags(ifrag_ref)%eval(jtmb)-((-1)**(ifrag))*lag_mult-homo_diff
           eval_tmp(jsforb+jtmb-ceiling(jstate_max))=tmb%orbs%eval(jsforb+jtmb-ceiling(jstate_max))+20.0d0
        end do

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
  if (.not. completely_random) then
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
     call f_free(eval_tmp)
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
        !!write(*,'(1x,a)') '-------------------------------------------------'
        !!write(*,'(1x,a)') 'some selected eigenvalues:'
        !!do iorb=1,tmb%orbs%norb!max(ksorbs%norb-8,1),min(ksorbs%norb+8,tmb%orbs%norb)
        !!    if(iorb==ksorbs%norb) then
        !!        write(*,'(3x,a,i0,a,es20.12,a)') 'eval(',iorb,')= ',tmb%orbs%eval(iorb),'  <-- last occupied orbital'
        !!    else if(iorb==ksorbs%norb+1) then
        !!        write(*,'(3x,a,i0,a,es20.12,a)') 'eval(',iorb,')= ',tmb%orbs%eval(iorb),'  <-- first virtual orbital'
        !!    else
        !!        write(*,'(3x,a,i0,a,es20.12)') 'eval(',iorb,')= ',tmb%orbs%eval(iorb)
        !!    end if
        !!end do
        !!write(*,'(1x,a)') '-------------------------------------------------'
        !!write(*,'(1x,a,2es24.16)') 'lowest, highest ev:',tmb%orbs%eval(1),tmb%orbs%eval(tmb%orbs%norb)

        call yaml_sequence_open('TMB eigenvalues',flow=.true.)
        call yaml_newline()
        do iorb=1,tmb%orbs%norb
            call yaml_mapping_open(flow=.true.)
            call yaml_map('index',iorb)
            call yaml_map('value',tmb%orbs%eval(iorb),fmt='(es20.12)')
            call yaml_mapping_close()
            if(iorb==ksorbs%norb) then
                !!write(*,'(3x,a,i0,a,es20.12,a)') 'eval(',iorb,')= ',tmb%orbs%eval(iorb),'  <-- last occupied orbital'
                call yaml_comment('  <-- last occupied orbital')
            else if(iorb==ksorbs%norb+1) then
                !!write(*,'(3x,a,i0,a,es20.12,a)') 'eval(',iorb,')= ',tmb%orbs%eval(iorb),'  <-- first virtual orbital'
                call yaml_comment('  <-- first virtual orbital')
            else
                !!write(*,'(3x,a,i0,a,es20.12)') 'eval(',iorb,')= ',tmb%orbs%eval(iorb)
            end if
            call yaml_newline()
        end do
        call yaml_sequence_close()
        !!write(*,'(1x,a)') '-------------------------------------------------'
        !!write(*,'(1x,a,2es24.16)') 'lowest, highest ev:',tmb%orbs%eval(1),tmb%orbs%eval(tmb%orbs%norb)

     end if

     if (nstates_max/=ksorbs%norb) then
        if (bigdft_mpi%iproc==0) print*,'Warning, number of states with non-zero occupation in fragments (',nstates_max,&
             ') differs from number of KS states (',ksorbs%norb,') - might have convergence problems'
     end if

     !!!!!!!!!!!!!!!
     ! need the eigenvalues to be in ksorbs%eval
     call vcopy(ksorbs%norb,tmb%orbs%eval(1),1,ksorbs%eval(1),1)
     call evaltoocc(bigdft_mpi%iproc,bigdft_mpi%nproc,.false.,input%tel,ksorbs,input%occopt)

     nullify(mom_vec_fake)
     if (bigdft_mpi%iproc ==0) then 
        call write_eigenvalues_data(0.1d0,ksorbs,mom_vec_fake)
     end if
     !!!!!!!!!!!!!!!
  end if ! completely random
  call f_release_routine()
  call timing(iproc,'kernel_init','OF')

end subroutine fragment_coeffs_to_kernel
