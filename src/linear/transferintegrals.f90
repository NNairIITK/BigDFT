! calculation of cSc and cHc using original coeffs (HOMO and LUMO only) and new Hamiltonian and overlap matrices
subroutine calc_transfer_integrals_old(iproc,nproc,input_frag,ref_frags,orbs,ham,ovrlp)
  use module_base
  use module_types
  use yaml_output
  use module_fragments
  use internal_io
  use module_interfaces
  use sparsematrix_base, only: sparse_matrix
  implicit none

  integer, intent(in) :: iproc, nproc
  type(fragmentInputParameters), intent(in) :: input_frag
  type(system_fragment), dimension(input_frag%nfrag_ref), intent(in) :: ref_frags
  type(orbitals_data), intent(in) :: orbs
  type(sparse_matrix), intent(inout) :: ham, ovrlp
  !Local variables
  character(len=*), parameter :: subname='calc_transfer_integrals'
  integer :: i_stat, i_all, ifrag, jfrag, ntmb_tot, ind, itmb, ifrag_ref, ierr, ih, jh
  !integer :: jfrag_ref, jtmb
  integer, allocatable, dimension(:) :: homo
  real(gp), allocatable, dimension(:,:) :: homo_coeffs

  allocate(homo(input_frag%nfrag), stat=i_stat)
  call memocc(i_stat, homo, 'homo', subname)

  do ifrag=1,input_frag%nfrag
     ifrag_ref=input_frag%frag_index(ifrag)
     homo(ifrag)=ceiling(ref_frags(ifrag_ref)%nelec/2.0_gp)
  end do

  ntmb_tot=ham%nfvctr!=orbs%norb
  allocate(homo_coeffs(ntmb_tot,input_frag%nfrag), stat=i_stat)
  call memocc(i_stat, homo_coeffs, 'homo_coeffs', subname)

  if (input_frag%nfrag/=2) stop 'Error, only 2 fragments may currently be considered for transfer integral calculation'
  ! activate site energies only in case of more fragments

  if (iproc==0) write(*,*) 'HOMO and LUMO are defined as those of the neutral fragment'

  ! combine individual homo coeffs into a big ntmb_tot x input_frag%nfrag array

  ind=ref_frags(input_frag%frag_index(1))%fbasis%forbs%norb
  do ih=-1,2
     if (homo(input_frag%frag_index(1))+ih>ref_frags(input_frag%frag_index(1))%fbasis%forbs%norb) cycle

     call to_zero(ntmb_tot, homo_coeffs(1,1))

     do itmb=1,ref_frags(input_frag%frag_index(1))%fbasis%forbs%norb
        homo_coeffs(itmb,1)=ref_frags(input_frag%frag_index(1))%coeff(itmb,homo(input_frag%frag_index(1))+ih)
     end do

     do jh=-1,2
        if (homo(input_frag%frag_index(2))+jh>ref_frags(input_frag%frag_index(2))%fbasis%forbs%norb) cycle     

        call to_zero(ntmb_tot, homo_coeffs(1,2))

        do itmb=1,ref_frags(input_frag%frag_index(2))%fbasis%forbs%norb
           homo_coeffs(ind+itmb,2)=ref_frags(input_frag%frag_index(2))%coeff(itmb,homo(input_frag%frag_index(2))+jh)
        end do

        if (iproc==0) then
           if (ih<0) then
              write(*,'(a,I2)',advance='NO') 'Fragment 1 HOMO-',abs(ih)
           else if (ih==0) then
              write(*,'(a)',advance='NO') 'Fragment 1 HOMO'
           else if (ih==1) then
              write(*,'(a)',advance='NO') 'Fragment 1 LUMO'
           else
              write(*,'(a,I2)',advance='NO') 'Fragment 1 LUMO+',ih-1
           end if
        end if

        if (iproc==0) then
           if (jh<0) then
              write(*,'(a,I2,a)') ', fragment 2 HOMO-',abs(jh),'.  '
           else if (jh==0) then
              write(*,'(a)') ', fragment 2 HOMO.  '
           else if (jh==1) then
              write(*,'(a)') ', fragment 2 LUMO.  '
           else
              write(*,'(a,I2)') ', fragment 2 LUMO+',jh-1,'.  '
           end if
        end if

        call calc_transfer_integral_old(iproc,nproc,input_frag,orbs,ham,ovrlp,homo_coeffs)

     end do
  end do

  i_all = -product(shape(homo_coeffs))*kind(homo_coeffs)
  deallocate(homo_coeffs,stat=i_stat)
  call memocc(i_stat,i_all,'homo_coeffs',subname)

  i_all = -product(shape(homo))*kind(homo)
  deallocate(homo,stat=i_stat)
  call memocc(i_stat,i_all,'homo',subname)

end subroutine calc_transfer_integrals_old





! calculation of cSc and cHc using original coeffs (HOMO and LUMO only) and new Hamiltonian and overlap matrices
subroutine calc_transfer_integral_old(iproc,nproc,input_frag,orbs,ham,ovrlp,homo_coeffs)
  use module_base
  use module_types
  use yaml_output
  use module_fragments
  use internal_io
  use module_interfaces
  use sparsematrix_base, only: sparse_matrix
  use sparsematrix, only: uncompress_matrix
  implicit none

  integer, intent(in) :: iproc, nproc
  type(fragmentInputParameters), intent(in) :: input_frag
  type(orbitals_data), intent(in) :: orbs
  type(sparse_matrix), intent(inout) :: ham, ovrlp
  real(kind=gp), dimension(ovrlp%nfvctr,input_frag%nfrag), intent(in) :: homo_coeffs
  !Local variables
  character(len=*), parameter :: subname='calc_transfer_integral'
  integer :: i_stat, i_all, ifrag, jfrag, ntmb_tot, ind, itmb, ierr, i, j
  !integer :: jfrag_ref, jtmb
  real(gp), allocatable, dimension(:,:) :: homo_ham, homo_ovrlp, coeff_tmp
  real(gp) :: orthog_energy


  ! make the coeff copies more efficient?

  allocate(coeff_tmp(orbs%norbp,input_frag%nfrag), stat=i_stat)
  call memocc(i_stat, coeff_tmp, 'coeff_tmp', subname)

  allocate(homo_ham(input_frag%nfrag,input_frag%nfrag), stat=i_stat)
  call memocc(i_stat, homo_ham, 'homo_ham', subname)
  allocate(homo_ovrlp(input_frag%nfrag,input_frag%nfrag), stat=i_stat)
  call memocc(i_stat, homo_ovrlp, 'homo_ovrlp', subname)

  !!allocate(ham%matrix(ham%nfvctr,ham%nfvctr), stat=i_stat)
  !!call memocc(i_stat, ham%matrix, 'ham%matrix', subname)
  ham%matrix=f_malloc_ptr((/ham%nfvctr,ham%nfvctr/),id='ham%matrix')
  
  call uncompress_matrix(iproc,ham)

  !DGEMM(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
  !rows op(a) and c, cols op(b) and c, cols op(a) and rows op(b)
  call to_zero(input_frag%nfrag**2, homo_ham(1,1))
  if (orbs%norbp>0) then
     call dgemm('n', 'n', orbs%norbp, input_frag%nfrag, orbs%norb, 1.d0, &
          ham%matrix(orbs%isorb+1,1),orbs%norb, &
          homo_coeffs(1,1), orbs%norb, 0.d0, &
          coeff_tmp, orbs%norbp)
     call dgemm('t', 'n', input_frag%nfrag, input_frag%nfrag, orbs%norbp, 1.d0, homo_coeffs(orbs%isorb+1,1), &
          orbs%norb, coeff_tmp, orbs%norbp, 0.d0, homo_ham, input_frag%nfrag)
  end if


  if (nproc>1) then
      call mpiallred(homo_ham(1,1), input_frag%nfrag**2, mpi_sum, bigdft_mpi%mpi_comm)
  end if

  !!i_all=-product(shape(ham%matrix))*kind(ham%matrix)
  !!deallocate(ham%matrix, stat=i_stat)
  !!call memocc(i_stat, i_all, 'ham%matrix', subname)
  call f_free_ptr(ham%matrix)

  !!allocate(ovrlp%matrix(ovrlp%nfvctr,ovrlp%nfvctr), stat=i_stat)
  !!call memocc(i_stat, ovrlp%matrix, 'ovrlp%matrix', subname)
  ovrlp%matrix=f_malloc_ptr((/ovrlp%nfvctr,ovrlp%nfvctr/),id='ovrlp%matrix')
  call uncompress_matrix(iproc,ovrlp)

  call to_zero(input_frag%nfrag**2, homo_ovrlp(1,1))
  if (orbs%norbp>0) then
     call dgemm('n', 'n', orbs%norbp, input_frag%nfrag, orbs%norb, 1.d0, ovrlp%matrix(orbs%isorb+1,1), &
          orbs%norb, homo_coeffs(1,1), orbs%norb, 0.d0, coeff_tmp, orbs%norbp)
     call dgemm('t', 'n', input_frag%nfrag, input_frag%nfrag, orbs%norbp, 1.d0, homo_coeffs(orbs%isorb+1,1), &
          orbs%norb, coeff_tmp, orbs%norbp, 0.d0, homo_ovrlp, input_frag%nfrag)
  end if

  if (nproc>1) then
      call mpiallred(homo_ovrlp(1,1), input_frag%nfrag**2, mpi_sum, bigdft_mpi%mpi_comm)
  end if

  !!i_all=-product(shape(ovrlp%matrix))*kind(ovrlp%matrix)
  !!deallocate(ovrlp%matrix, stat=i_stat)
  !!call memocc(i_stat, i_all, 'ovrlp%matrix', subname)
  call f_free_ptr(ovrlp%matrix)

  i_all = -product(shape(coeff_tmp))*kind(coeff_tmp)
  deallocate(coeff_tmp,stat=i_stat)
  call memocc(i_stat,i_all,'coeff_tmp',subname)

  ! output results
  !if (iproc==0) write(*,'(a)') '-----------------------------------------------------------------------------------------'
  if (input_frag%nfrag/=2) then
     !!if (iproc==0) write(*,*) 'Transfer integrals and site energies:'
     !!if (iproc==0) write(*,*) 'frag i, frag j, energy, overlap'
     if (iproc==0) then
         call yaml_open_map('Transfer integrals and site energies')
     end if
     do jfrag=1,input_frag%nfrag
        do ifrag=1,input_frag%nfrag
           !!if (iproc==0) write(*,'(2(I5,1x),1x,2(F16.12,1x))') jfrag, ifrag, homo_ham(jfrag,ifrag), homo_ovrlp(jfrag,ifrag)
           if (iproc==0) then
               call yaml_map('frag j',jfrag)
               call yaml_map('frag i',ifrag)
               call yaml_map('energy',homo_ham(jfrag,ifrag),fmt='(f16.12)')
               call yaml_map('overlap',homo_ovrlp(jfrag,ifrag),fmt='(f16.12)')
           end if
        end do
     end do
     call yaml_close_map()
  else ! include orthogonalized results as well
     !if (iproc==0) write(*,*) 'Site energies:'
     !if (iproc==0) write(*,*) 'frag i, energy, overlap, orthog energy'
     !i=1
     !j=2
     !orthog_energy= (0.5_gp/(1.0_gp-homo_ovrlp(i,j)**2)) &
     !             * ( (homo_ham(i,i)+homo_ham(j,j)) - 2.0_gp*homo_ham(i,j)*homo_ovrlp(i,j) &
     !             + (homo_ham(i,i)-homo_ham(j,j))*dsqrt(1.0_gp-homo_ovrlp(i,j)**2) )
     !if (iproc==0) write(*,'((I5,1x),1x,3(F16.12,1x))') 1, homo_ham(1,1), homo_ovrlp(1,1), orthog_energy
     !orthog_energy= (0.5_gp/(1.0_gp-homo_ovrlp(i,j)**2)) &
     !             * ( (homo_ham(i,i)+homo_ham(j,j)) - 2.0_gp*homo_ham(i,j)*homo_ovrlp(i,j) &
     !             - (homo_ham(i,i)-homo_ham(j,j))*dsqrt(1.0_gp-homo_ovrlp(i,j)**2) )
     !if (iproc==0) write(*,'((I5,1x),1x,3(F16.12,1x))') 2, homo_ham(2,2), homo_ovrlp(2,2), orthog_energy

     !!if (iproc==0) write(*,*) 'Transfer integrals:'
     !!if (iproc==0) write(*,*) 'frag i, frag j, energy, overlap, orthog energy'
     i=1
     j=2
     orthog_energy=(homo_ham(i,j)-0.5_gp*(homo_ham(i,i)+homo_ham(j,j))*homo_ovrlp(i,j))/(1.0_gp-homo_ovrlp(i,j)**2)
     !!if (iproc==0) write(*,'(2(I5,1x),1x,3(F16.12,1x))') 1, 2, homo_ham(1,2), homo_ovrlp(1,2),orthog_energy
     i=2
     j=1
     orthog_energy=(homo_ham(i,j)-0.5_gp*(homo_ham(i,i)+homo_ham(j,j))*homo_ovrlp(i,j))/(1.0_gp-homo_ovrlp(i,j)**2)
     !!if (iproc==0) write(*,'(2(I5,1x),1x,3(F16.12,1x))') 1, 2, homo_ham(2,1), homo_ovrlp(2,1),orthog_energy
     if (iproc==0) then
         call yamL_open_map('Transfer integrals')
         call yaml_map('frag i',1)
         call yaml_map('frag j',2)
         call yaml_map('energy',homo_ham(1,2))
         call yaml_map('overlap',homo_ovrlp(1,2))
         call yaml_map('orthog energy',orthog_energy)
         call yaml_map('frag i',1)
         call yaml_map('frag j',2)
         call yaml_map('energy',homo_ham(2,1))
         call yaml_map('overlap',homo_ovrlp(2,1))
         call yaml_map('orthog energy',orthog_energy)
         call yaml_close_map()
     end if

  end if
  !!if (iproc==0) write(*,'(a)') '-----------------------------------------------------------------------------------------'

  i_all = -product(shape(homo_ham))*kind(homo_ham)
  deallocate(homo_ham,stat=i_stat)
  call memocc(i_stat,i_all,'homo_ham',subname)
  i_all = -product(shape(homo_ovrlp))*kind(homo_ovrlp)
  deallocate(homo_ovrlp,stat=i_stat)
  call memocc(i_stat,i_all,'homo_ovrlp',subname)

end subroutine calc_transfer_integral_old


! calculation of cSc and cHc using original coeffs and new Hamiltonian and overlap matrices
! parallelization to be improved
! also have already uncompressed and recompressed ovrlp, so could change this
subroutine calc_transfer_integral(iproc,nproc,nstates,orbs,ham,ovrlp,homo_coeffs1,homo_coeffs2,homo_ham,homo_ovrlp)
  use module_base
  use module_types
  use yaml_output
  use module_fragments
  use internal_io
  use module_interfaces
  use sparsematrix_base, only: sparse_matrix
  implicit none

  integer, intent(in) :: iproc, nproc, nstates
  type(orbitals_data), intent(in) :: orbs
  type(sparse_matrix), intent(inout) :: ham, ovrlp
  real(kind=gp), dimension(ovrlp%nfvctr,nstates), intent(in) :: homo_coeffs1, homo_coeffs2
  real(kind=gp), dimension(nstates), intent(inout) :: homo_ham, homo_ovrlp

  !Local variables
  integer :: i_stat, i_all, ifrag, jfrag, ntmb_tot, ind, itmb, ierr, i, j, istate
  real(gp), allocatable, dimension(:,:) :: coeff_tmp
  real(gp) :: orthog_energy

  coeff_tmp=f_malloc((/orbs%norbp,nstates/), id='coeff_tmp')

  !DGEMM(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
  !rows op(a) and c, cols op(b) and c, cols op(a) and rows op(b)
  !ham%matrix=f_malloc_ptr((/ham%nfvctr,ham%nfvctr/), id='ham%matrix')
  !call uncompress_matrix(iproc,ham)
  if (orbs%norbp>0) then
     do istate=1,nstates
        call dgemm('n', 'n', orbs%norbp, 1, orbs%norb, 1.d0, &
             ham%matrix(orbs%isorb+1,1),orbs%norb, &
             homo_coeffs1(1,istate), orbs%norb, 0.d0, &
             coeff_tmp(1,istate), orbs%norbp)
        call dgemm('t', 'n', 1, 1, orbs%norbp, 1.d0, homo_coeffs2(orbs%isorb+1,istate), &
             orbs%norb, coeff_tmp(1,istate), orbs%norbp, 0.d0, homo_ham(istate), 1)
     end do
  else
     call to_zero(nstates,homo_ham(1))
  end if

  if (nproc>1) then
      call mpiallred(homo_ham(1), nstates, mpi_sum, bigdft_mpi%mpi_comm)
  end if

  !call f_free_ptr(ham%matrix)

  !ovrlp%matrix=f_malloc_ptr((/ovrlp%nfvctr,ovrlp%nfvctr/), id='ovrlp%matrix')
  !call uncompress_matrix(iproc,ovrlp)

  if (orbs%norbp>0) then
     do istate=1,nstates
        call dgemm('n', 'n', orbs%norbp, 1, orbs%norb, 1.d0, ovrlp%matrix(orbs%isorb+1,1), &
             orbs%norb, homo_coeffs1(1,istate), orbs%norb, 0.d0, coeff_tmp(1,istate), orbs%norbp)
        call dgemm('t', 'n', 1, 1, orbs%norbp, 1.d0, homo_coeffs2(orbs%isorb+1,istate), &
             orbs%norb, coeff_tmp(1,istate), orbs%norbp, 0.d0, homo_ovrlp(istate), 1)
     end do
  else
     call to_zero(nstates,homo_ovrlp(1))
  end if

  if (nproc>1) then
      call mpiallred(homo_ovrlp(1), nstates, mpi_sum, bigdft_mpi%mpi_comm)
  end if

  !call f_free_ptr(ovrlp%matrix)
  call f_free(coeff_tmp)

end subroutine calc_transfer_integral


! calculation of cSc and cHc using original coeffs and new Hamiltonian and overlap matrices
! parallelization to be improved
! only calculates transfer integrals if we have two fragments
! occs are for neutral reference fragments...
subroutine calc_site_energies_transfer_integrals(iproc,nproc,meth_overlap,input_frag,ref_frags,orbs,ham,ovrlp)
  use module_base
  use module_types
  use yaml_output
  use module_fragments
  use internal_io
  use module_interfaces
  use sparsematrix_base, only: sparse_matrix
  use sparsematrix, only: uncompress_matrix
  implicit none

  integer, intent(in) :: iproc, nproc, meth_overlap
  type(fragmentInputParameters), intent(in) :: input_frag
  type(orbitals_data), intent(in) :: orbs
  type(sparse_matrix), intent(inout) :: ham, ovrlp
  type(system_fragment), dimension(input_frag%nfrag_ref), intent(in) :: ref_frags
  !Local variables
  integer :: i_stat, i_all, ifrag, jfrag, ntmb_tot, ind, itmb, ierr, i, j, nstates, istate, ih, ifrag_ref
  integer :: ifrag_ref1, ifrag_ref2, homo1, homo2, jh, above_lumo, iind, jind, norb_tmp
  !integer :: jfrag_ref, jtmb
  real(gp), allocatable, dimension(:,:) :: coeffs_tmp, homo_coeffs, coeffs_orthog
  real(gp), allocatable, dimension(:) :: frag_sum, homo_ham, homo_ovrlp
  real(gp), allocatable, dimension(:) :: frag_sum_orthog, homo_ham_orthog, homo_ovrlp_orthog, eval_sum
  real(gp) :: frag_sum_tot, frag_sum_tot_orthog, eval_sum_tot, orthog_energy
  real(gp), dimension(1) :: trans_int_energy, trans_int_energy_orthog, trans_int_ovrlp, trans_int_ovrlp_orthog
  character(len=8) :: str
  logical, parameter :: separate_site_energies=.false.

  call timing(iproc,'transfer_int','ON')
  call f_routine(id='calc_site_energies_transfer_integrals')

  nstates=0
  above_lumo=0
  do ifrag=1,input_frag%nfrag
     ifrag_ref= input_frag%frag_index(ifrag)
     nstates=nstates+min(ceiling((ref_frags(ifrag_ref)%nelec+1)/2.0_gp)+above_lumo,ref_frags(ifrag_ref)%fbasis%forbs%norb)
  end do

  homo_ham=f_malloc(nstates,id='homo_ham')
  homo_ovrlp=f_malloc(nstates,id='homo_ovrlp')
  homo_coeffs=f_malloc0((/ovrlp%nfvctr,nstates/), id='homo_coeffs')
  !coeffs_tmp=f_malloc((/ovrlp%nfvctr,ovrlp%nfvctr/), id='coeffs_tmp')
  ovrlp%matrix=f_malloc_ptr((/ovrlp%nfvctr,ovrlp%nfvctr/), id='ovrlp%matrix')
  call uncompress_matrix(iproc,ovrlp)

  istate=1
  ind=1
  do ifrag=1,input_frag%nfrag
     ifrag_ref=input_frag%frag_index(ifrag)

     norb_tmp=min(ceiling((ref_frags(ifrag_ref)%nelec+1)/2.0_gp)+above_lumo,ref_frags(ifrag_ref)%fbasis%forbs%norb)

     !call to_zero(ovrlp%nfvctr**2,coeffs_tmp(1,1),1)
     do ih=1,norb_tmp
        call vcopy(ref_frags(ifrag_ref)%fbasis%forbs%norb,ref_frags(ifrag_ref)%coeff(1,ih),1,homo_coeffs(ind,istate+ih-1),1)
        !call vcopy(ref_frags(ifrag_ref)%fbasis%forbs%norb,ref_frags(ifrag_ref)%coeff(1,ih),1,coeffs_tmp(ind,ih),1)
     end do

     !call reorthonormalize_coeff(iproc, nproc, norb_tmp, -8, -8, meth_overlap, orbs, ovrlp, coeffs_tmp(1,1))
     !call vcopy(orbs%norb*norb_tmp,coeffs_tmp(1,1),1,homo_coeffs(1,istate),1)

     istate=istate+norb_tmp
     ind=ind+ref_frags(ifrag_ref)%fbasis%forbs%norb
  end do
  !call f_free(coeffs_tmp)

  ham%matrix=f_malloc_ptr((/ham%nfvctr,ham%nfvctr/), id='ham%matrix')
  call uncompress_matrix(iproc,ham)
  if (separate_site_energies .or. input_frag%nfrag==1) call calc_transfer_integral(iproc,nproc,nstates,&
       orbs,ham,ovrlp,homo_coeffs,homo_coeffs,homo_ham,homo_ovrlp)

  ! orthogonalize
  coeffs_tmp=f_malloc0((/orbs%norb,orbs%norb/), id='coeffs_tmp')
  call vcopy(orbs%norb*nstates,homo_coeffs(1,1),1,coeffs_tmp(1,1),1)
  call reorthonormalize_coeff(iproc, nproc, nstates, -8, -8, meth_overlap, orbs, ovrlp, coeffs_tmp(1,1), orbs)
  coeffs_orthog=f_malloc((/orbs%norb,nstates/), id='coeffs_orthog')
  call vcopy(orbs%norb*nstates,coeffs_tmp(1,1),1,coeffs_orthog(1,1),1)
  call f_free(coeffs_tmp)

  ! only calculate site energies separately if specified or if not calculating them below
  if (separate_site_energies .or. input_frag%nfrag==1) then
     homo_ham_orthog=f_malloc(nstates, id='homo_ham_orthog')
     homo_ovrlp_orthog=f_malloc(nstates, id='homo_ovrlp_orthog')

     call calc_transfer_integral(iproc,nproc,nstates,orbs,ham,ovrlp,coeffs_orthog,coeffs_orthog,&
          homo_ham_orthog,homo_ovrlp_orthog)

     frag_sum=f_malloc0(nstates, id='frag_sum')
     frag_sum_orthog=f_malloc0(nstates, id='frag_sum_orthog')
     eval_sum=f_malloc0(nstates, id='eval_sum')

     !if (iproc==0) write(*,'(a)') '-------------------------------------------------------------------------------------------------'
     !if (iproc==0) write(*,*) 'Site energies:'
     if (iproc==0) call yaml_open_sequence('Site energies',flow=.true.)

     !!if (iproc==0) write(*,*) 'state, energy, orthog energy, frag eval, overlap, orthog overlap, occ'
     if (iproc==0) call yaml_comment('state, energy, orthog energy, frag eval, overlap, orthog overlap, occ')

     istate=1
     frag_sum_tot=0
     frag_sum_tot_orthog=0
     eval_sum_tot=0
     do ifrag=1,input_frag%nfrag
        ifrag_ref=input_frag%frag_index(ifrag)
        !if (iproc==0) write(*,'(a,i3)') trim(input_frag%label(ifrag_ref)),ifrag
        if (iproc==0) call yaml_open_map(flow=.true.)
        if (iproc==0) call yaml_map('label',trim(input_frag%label(ifrag_ref)))
        do ih=1,min(ceiling((ref_frags(ifrag_ref)%nelec+1)/2.0_gp)+above_lumo,ref_frags(ifrag_ref)%fbasis%forbs%norb)
           !!if (iproc==0) call yaml_open_map(flow=.true.)
           if (iproc==0) call yaml_newline()
           if (ih<ceiling(ref_frags(ifrag_ref)%nelec/2.0_gp)) then
              write(str,'(I2)') abs(ih-ceiling(ref_frags(ifrag_ref)%nelec/2.0_gp))
              !if (iproc==0) write(*,'(a8)',advance='NO') ' HOMO-'//trim(adjustl(str))
              if (iproc==0) call yaml_map('state','HOMO-'//trim(adjustl(str)))
           else if (ih==ceiling(ref_frags(ifrag_ref)%nelec/2.0_gp)) then
              !if (iproc==0) write(*,'(a8)',advance='NO') ' HOMO'
              if (iproc==0) call yaml_map('state','HOMO')
           else if (ih==ceiling(ref_frags(ifrag_ref)%nelec/2.0_gp)+1) then
              !if (iproc==0) write(*,'(a8)',advance='NO') ' LUMO'
              if (iproc==0) call yaml_map('state','LUMO')
           else
              write(str,'(I2)') ih-1-ceiling(ref_frags(ifrag_ref)%nelec/2.0_gp)
              !if (iproc==0) write(*,'(a8)',advance='NO') ' LUMO+'//trim(adjustl(str))
              if (iproc==0) call yaml_map('state','LUMO+'//trim(adjustl(str)))
           end if
           !if (iproc==0) write(*,'(1x,5(F20.12,1x))',advance='NO') homo_ham(istate), homo_ham_orthog(istate), &
           !     ref_frags(ifrag_ref)%eval(ih), homo_ovrlp(istate), homo_ovrlp_orthog(istate)
           if (iproc==0) then
               call yaml_map('energy',homo_ham(istate),fmt='(es16.8)')
               call yaml_map('orthog energy',homo_ham_orthog(istate),fmt='(es16.8)')
               call yaml_map('frag eval',ref_frags(ifrag_ref)%eval(ih),fmt='(es16.8)')
               call yaml_map('overlap',homo_ovrlp(istate),fmt='(es14.6)')
               !call yaml_map('orthog overlap',homo_ovrlp_orthog(istate))
           end if     
           if (ih<ceiling(ref_frags(ifrag_ref)%nelec/2.0_gp)) then
              frag_sum(ifrag)=frag_sum(ifrag)+homo_ham(istate)
              frag_sum_orthog(ifrag)=frag_sum_orthog(ifrag)+homo_ham_orthog(istate)
              eval_sum(ifrag)=eval_sum(ifrag)+ref_frags(ifrag_ref)%eval(ih)
              !if (iproc==0) write(*,'(1x,F4.2)') 2.0_gp
              if (iproc==0) call yaml_map('occ',2.0_gp,fmt='(f5.2)')
           else if (ih==ceiling(ref_frags(ifrag_ref)%nelec/2.0_gp)) then
              if (mod(real(ref_frags(ifrag_ref)%nelec,gp),2.0_gp)/=0.0_gp) then
                 frag_sum(ifrag)=frag_sum(ifrag)+0.5_gp*homo_ham(istate)
                 frag_sum_orthog(ifrag)=frag_sum_orthog(ifrag)+0.5_gp*homo_ham_orthog(istate)
                 eval_sum(ifrag)=eval_sum(ifrag)+0.5_gp*ref_frags(ifrag_ref)%eval(ih)
                 !if (iproc==0) write(*,'(1x,F4.2)') 1.0_gp
                 if (iproc==0) call yaml_map('occ',1.0_gp,fmt='(f5.2)')
              else
                 frag_sum(ifrag)=frag_sum(ifrag)+homo_ham(istate)
                 frag_sum_orthog(ifrag)=frag_sum_orthog(ifrag)+homo_ham_orthog(istate)
                 eval_sum(ifrag)=eval_sum(ifrag)+ref_frags(ifrag_ref)%eval(ih)
                 !if (iproc==0) write(*,'(1x,F4.2)') 2.0_gp
                 if (iproc==0) call yaml_map('occ',2.0_gp,fmt='(f5.2)')
              end if
           else
              !if (iproc==0) write(*,'(1x,F4.2)') 0.0_gp
              if (iproc==0) call yaml_map('occ',0.0_gp,fmt='(f5.2)')
           end if
           istate=istate+1
           !!if (iproc==0) call yaml_close_map()
        end do
        !if (iproc==0) write(*,'(9x,3(F20.12,1x))') 2.0_gp*frag_sum(ifrag),&
        !     2.0_gp*frag_sum_orthog(ifrag),2.0_gp*eval_sum(ifrag)
          !if (iproc==0) write(*,'(a)') '------------------------------------------------------------------------'//&
          !     '-------------------------'
        if(iproc==0) then
            call yaml_newline
            call yaml_map('2*frag sum',2.0_gp*frag_sum(ifrag))
            call yaml_map('2*frag sum orthog',2.0_gp*frag_sum_orthog(ifrag))
            call yaml_map('2*eval sum',2.0_gp*eval_sum(ifrag))
            call yaml_close_map()
            call yaml_newline()
        end if
        frag_sum_tot=frag_sum_tot+frag_sum(ifrag)
        frag_sum_tot_orthog=frag_sum_tot_orthog+frag_sum_orthog(ifrag)
        eval_sum_tot=eval_sum_tot+eval_sum(ifrag)
     end do
     if (iproc==0) call yaml_close_sequence()

     if (iproc==0) then
         call yaml_map('2.0_gp*frag_sum_tot',2.0_gp*frag_sum_tot)
         call yaml_map('2.0_gp*frag_sum_tot_orthog',2.0_gp*frag_sum_tot_orthog)
         call yaml_map('2.0_gp*eval_sum_tot',2.0_gp*eval_sum_tot)
     end if

     !if (iproc==0) write(*,'(9x,3(F20.12,1x))') 2.0_gp*frag_sum_tot, 2.0_gp*frag_sum_tot_orthog,2.0_gp*eval_sum_tot
     !if (iproc==0) write(*,'(a)') '-------------------------------------------------------------------------------------------------'

     call f_free(eval_sum)
     call f_free(frag_sum)
     call f_free(frag_sum_orthog)
     call f_free(homo_ham_orthog)
     call f_free(homo_ovrlp_orthog)
  else
     orthog_energy=0.0d0
  end if

  if (input_frag%nfrag>=2) then
     !if (iproc==0) write(*,*) 'Transfer integrals (HOMO and LUMO are defined as those of the neutral fragment):'
     if (iproc==0) call yaml_open_sequence('Transfer integrals &
         &(HOMO and LUMO are defined as those of the neutral fragment)',flow=.true.)
     if (iproc==0) call yaml_newline()
     !if (iproc==0) write(*,*) 'state1, state2, energy, orthog energy, orthog energy2, overlap, orthog overlap, occ1, occ2'
     if (iproc==0) call yaml_comment('state1, state2, energy, orthog energy, &
         &orthog energy2, overlap, orthog overlap, occ1, occ2')
     iind=0
     do ifrag=1,input_frag%nfrag
        ifrag_ref1=input_frag%frag_index(ifrag)
        homo1=ceiling((ref_frags(ifrag_ref1)%nelec)/2.0_gp)

        jind=0
        do jfrag=1,ifrag
           ifrag_ref2=input_frag%frag_index(jfrag)
           homo2=ceiling((ref_frags(ifrag_ref2)%nelec)/2.0_gp)

           do jh=-above_lumo,1+above_lumo
              if (homo2+jh>ref_frags(ifrag_ref2)%fbasis%forbs%norb) cycle  
              if (homo2+jh<1) cycle  
              do ih=-above_lumo,1+above_lumo
                 if (homo1+ih>ref_frags(ifrag_ref1)%fbasis%forbs%norb) cycle
                 if (homo1+ih<1) cycle  

                 i=homo1+ih+iind
                 j=homo2+jh+jind
                      
                 if (iproc==0) then
                    call yaml_open_map()!flow=.true.)
                    if (iproc==0) call yaml_newline()
                    if (ih<0) then
                       write(str,'(I2)') abs(ih)
                       !write(*,'(a,I3,a8)',advance='NO') trim(input_frag%label(ifrag_ref1)),ifrag,' HOMO-'//trim(adjustl(str))
                       call yaml_map('label',trim(input_frag%label(ifrag_ref1)))
                       call yaml_map('i',ifrag)
                       call yaml_map('s1','HOMO-'//trim(adjustl(str)))
                    else if (ih==0) then
                       !write(*,'(a,I3,a8)',advance='NO') trim(input_frag%label(ifrag_ref1)),ifrag,' HOMO  '
                       call yaml_map('label',trim(input_frag%label(ifrag_ref1)))
                       call yaml_map('i',ifrag)
                       call yaml_map('s1','HOMO')
                    else if (ih==1) then
                       !write(*,'(a,I3,a8)',advance='NO') trim(input_frag%label(ifrag_ref1)),ifrag,' LUMO  '
                       call yaml_map('label',trim(input_frag%label(ifrag_ref1)))
                       call yaml_map('i',ifrag)
                       call yaml_map('s1','LUMO')
                    else
                       write(str,'(I2)') ih-1
                       !write(*,'(a,I3,a8)',advance='NO') trim(input_frag%label(ifrag_ref1)),ifrag,' LUMO+'//trim(adjustl(str))
                       call yaml_map('label',trim(input_frag%label(ifrag_ref1)))
                       call yaml_map('i',ifrag)
                       call yaml_map('s1','LUMO+'//trim(adjustl(str)))
                    end if
                 end if

                 if (iproc==0) then
                    if (jh<0) then
                       write(str,'(I2)') abs(jh)
                       !write(*,'(3x,a,I3,a8)',advance='NO') trim(input_frag%label(ifrag_ref2)),jfrag,&
                       !     ' HOMO-'//trim(adjustl(str))
                       call yaml_map('label',trim(input_frag%label(ifrag_ref2)))
                       call yaml_map('j',jfrag)
                       call yaml_map('s1','HOMO-'//trim(adjustl(str)))
                    else if (jh==0) then
                       !write(*,'(3x,a,I3,a8)',advance='NO') trim(input_frag%label(ifrag_ref2)),jfrag,' HOMO  '
                       call yaml_map('label',trim(input_frag%label(ifrag_ref2)))
                       call yaml_map('j',jfrag)
                       call yaml_map('s1','HOMO')
                    else if (jh==1) then
                       !write(*,'(3x,a,I3,a8)',advance='NO') trim(input_frag%label(ifrag_ref2)),jfrag,' LUMO  '
                       call yaml_map('label',trim(input_frag%label(ifrag_ref2)))
                       call yaml_map('j',jfrag)
                       call yaml_map('s1','LUMO')
                    else
                       !write(str,'(I2)') jh-1
                       !write(*,'(3x,a,I3,a8)',advance='NO') trim(input_frag%label(ifrag_ref2)),jfrag,&
                       !     ' LUMO+'//trim(adjustl(str))
                       call yaml_map('label',trim(input_frag%label(ifrag_ref2)))
                       call yaml_map('j',jfrag)
                       call yaml_map('s1','LUMO+'//trim(adjustl(str)))
                    end if
                 end if

                 call calc_transfer_integral(iproc,nproc,1,orbs,ham,ovrlp,homo_coeffs(1,i),homo_coeffs(1,j),&
                      trans_int_energy(1),trans_int_ovrlp(1))
                 call calc_transfer_integral(iproc,nproc,1,orbs,ham,ovrlp,coeffs_orthog(1,i),coeffs_orthog(1,j),&
                      trans_int_energy_orthog(1),trans_int_ovrlp_orthog(1))

                 !orthog_energy=(trans_int_energy(1)-0.5_gp*(homo_ham(i)+homo_ham(j))*trans_int_ovrlp(1))&
                 !     /(1.0_gp-trans_int_ovrlp(1)**2)
      
                 !if (iproc==0) write(*,'(2x,5(F16.12,1x))',advance='NO') trans_int_energy(1), &
                 !     trans_int_energy_orthog(1), orthog_energy, trans_int_ovrlp(1), trans_int_ovrlp_orthog(1)
                 if (iproc==0) then
                     call yaml_map('trans_int_energy',trans_int_energy(1),fmt='(f16.12)')
                     call yaml_map('trans_int_energy_orthog',trans_int_energy_orthog(1),fmt='(f16.12)')
                     !call yamL_map('orthog_energy',orthog_energy,fmt='(f16.12)')
                     call yaml_map('trans_int_ovrlp',trans_int_ovrlp(1),fmt='(f16.12)')
                     call yaml_map('trans_int_ovrlp_orthog',trans_int_ovrlp_orthog(1),fmt='(f16.12)')
                 end if

                 if (homo1+ih<ceiling(ref_frags(ifrag_ref1)%nelec/2.0_gp)) then
                    !if (iproc==0) write(*,'(1x,F4.2)',advance='NO') 2.0_gp
                    if (iproc==0) call yaml_map('occ1',2.0_gp,fmt='(f4.2)')
                 else if (homo1+ih==ceiling(ref_frags(ifrag_ref1)%nelec/2.0_gp)) then
                    if (mod(real(ref_frags(ifrag_ref1)%nelec,gp),2.0_gp)/=0.0_gp) then
                       !if (iproc==0) write(*,'(1x,F4.2)',advance='NO') 1.0_gp
                       if (iproc==0) call yaml_map('occ1',1.0_gp,fmt='(f4.2)')
                    else
                       !if (iproc==0) write(*,'(1x,F4.2)',advance='NO') 2.0_gp
                       if (iproc==0) call yaml_map('occ1',2.0_gp,fmt='(f4.2)')
                    end if
                 else
                    !if (iproc==0) write(*,'(1x,F4.2)',advance='NO') 0.0_gp
                    if (iproc==0) call yaml_map('occ1',0.0_gp,fmt='(f4.2)')
                 end if

                 if (homo2+jh<ceiling(ref_frags(ifrag_ref2)%nelec/2.0_gp)) then
                    !if (iproc==0) write(*,'(1x,F4.2)') 2.0_gp
                    if (iproc==0) call yaml_map('occ2',2.0_gp,fmt='(f4.2)')
                 else if (homo2+jh==ceiling(ref_frags(ifrag_ref2)%nelec/2.0_gp)) then
                    if (mod(real(ref_frags(ifrag_ref2)%nelec,gp),2.0_gp)/=0.0_gp) then
                       !if (iproc==0) write(*,'(1x,F4.2)') 1.0_gp
                       if (iproc==0) call yaml_map('occ2',1.0_gp,fmt='(f4.2)')
                    else
                       !if (iproc==0) write(*,'(1x,F4.2)') 2.0_gp
                       if (iproc==0) call yaml_map('occ2',2.0_gp,fmt='(f4.2)')
                    end if
                 else
                    !if (iproc==0) write(*,'(1x,F4.2)') 0.0_gp
                    if (iproc==0) call yaml_map('occ2',0.0_gp,fmt='(f4.2)')
                 end if

                 if (iproc==0) call yaml_close_map()
                 if (iproc==0) call yaml_newline()

              end do
           end do
           !if (iproc==0) write(*,'(a)') '------------------------------------------------------------------------'//&
           !    '-------------------------'
           jind=jind+min(ceiling((ref_frags(ifrag_ref2)%nelec+1)/2.0_gp)+above_lumo,ref_frags(ifrag_ref2)%fbasis%forbs%norb)
        end do
        iind=iind+min(ceiling((ref_frags(ifrag_ref1)%nelec+1)/2.0_gp)+above_lumo,ref_frags(ifrag_ref1)%fbasis%forbs%norb)
     end do

     if (iproc==0) call yaml_close_sequence()
  end if

  call f_free_ptr(ham%matrix)
  call f_free_ptr(ovrlp%matrix)

  call f_free(homo_ham)
  call f_free(homo_ovrlp)
  call f_free(homo_coeffs)
  call f_free(coeffs_orthog)

  call f_release_routine()
  call timing(iproc,'transfer_int','OF')

end subroutine calc_site_energies_transfer_integrals
