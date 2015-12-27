module matrix_operations
    use module_base
    implicit none

    private

    !> Public routines
    public :: overlapPowerGeneral
    public :: overlap_minus_one_half_serial
    public :: overlap_plus_minus_one_half_exact
    public :: deviation_from_unity_parallel, deviation_from_unity_parallel_new
    public :: overlap_power_minus_one_half_parallel
    public :: check_taylor_order
    public :: calculate_S_minus_one_half_onsite


    contains


      !> S^-1 exact only works for symmetric matrices
      !! BOTH sparse matrices must be present together and inv_ovrlp should be nullified pointer, NOT inv_ovrlp_smat%matrix
      !! when sparse matrices present, check is performed to see whether %matrix is allocated so that its allocated status remains unchanged
      !! contents of %matrix not guaranteed to be correct though
      !! power: -2 -> S^-1/2, 2 -> S^1/2, 1 -> S^-1
      subroutine overlapPowerGeneral(iproc, nproc, iorder, ncalc, power, blocksize, imode, &
                 ovrlp_smat, inv_ovrlp_smat, ovrlp_mat, inv_ovrlp_mat, check_accur, &
                 max_error, mean_error, nspinx)
           !!foe_nseg, foe_kernel_nsegline, foe_istsegline, foe_keyg)
        use module_base
        use sparsematrix_base, only: sparse_matrix, matrices, &
                                sparsematrix_malloc_ptr, sparsematrix_malloc, sparsematrix_malloc0, sparsematrix_malloc0_ptr, &
                                assignment(=), &
                                SPARSE_FULL, DENSE_PARALLEL, SPARSE_MATMUL_LARGE, &
                                DENSE_MATMUL, DENSE_FULL, SPARSEMM_SEQ, SPARSE_TASKGROUP
        use sparsematrix, only: compress_matrix, uncompress_matrix, &
                                transform_sparse_matrix, transform_sparse_matrix_local, &
                                compress_matrix_distributed_wrapper, &
                                uncompress_matrix_distributed2, &
                                sequential_acces_matrix_fast2, sequential_acces_matrix_fast, &
                                gather_matrix_from_taskgroups, gather_matrix_from_taskgroups_inplace, &
                                uncompress_matrix2, transform_sparsity_pattern, &
                                sparsemm_new, matrix_matrix_mult_wrapper
        use parallel_linalg, only: dpotrf_parallel, dpotri_parallel
        use ice, only: inverse_chebyshev_expansion
        use yaml_output
        implicit none
        
        ! Calling arguments
        integer,intent(in) :: iproc, nproc, iorder, blocksize, ncalc
        integer,dimension(ncalc),intent(in) :: power
        integer,intent(in) :: imode
        type(sparse_matrix),intent(in) :: ovrlp_smat, inv_ovrlp_smat
        type(matrices),intent(in) :: ovrlp_mat
        type(matrices),dimension(ncalc),intent(inout) :: inv_ovrlp_mat
        logical,intent(in) :: check_accur
        real(kind=8),intent(out),optional :: max_error, mean_error
        integer,intent(in),optional :: nspinx !< overwrite the default spin value
        
        ! Local variables
        integer :: iorb, jorb, info, iiorb, isorb, norbp, ii, ii_inv, iii, ierr, i, its, maxits
        integer :: matrixindex_in_compressed, nmaxvalk, icalc
        real(kind=8), dimension(:,:), pointer :: inv_ovrlpp, ovrlppowerp
        real(kind=8), dimension(:,:), pointer :: inv_ovrlp_half_tmp
        real(kind=8), dimension(:), pointer :: ovrlpminonep_new
        real(kind=8), dimension(:,:,:), pointer :: ovrlpminone, ovrlp_local, inv_ovrlp_local, ovrlppoweroldp, ovrlpminonep
        real(kind=8) :: factor
        logical :: ovrlp_allocated, inv_ovrlp_allocated
      
        ! new for sparse taylor
        integer :: nout, nseq, ispin, ishift, ishift2, isshift, ilshift, ilshift2, nspin, iline, icolumn, ist, j, iorder_taylor
        integer,dimension(:,:,:),allocatable :: istindexarr
        real(kind=8),dimension(:),pointer :: ovrlpminone_sparse
        real(kind=8),dimension(:),allocatable :: ovrlp_compr_seq, ovrlpminone_sparse_seq, ovrlp_large_compr, tmparr, resmat
        real(kind=8),dimension(:),allocatable :: invovrlp_compr_seq, ovrlpminoneoldp_new
        real(kind=8),dimension(:,:),allocatable :: ovrlpminoneoldp, invovrlpp, ovrlp_largep
        real(kind=8),dimension(:,:,:),allocatable :: invovrlpp_arr
        real(kind=8),dimension(:,:),allocatable :: invovrlpp_arr_new
        real(kind=8),dimension(:,:),allocatable :: Amat12p, Amat21p, Amat21
        real(kind=8),dimension(:,:),pointer :: Amat12, Amat11p, Amat22p
        real(kind=8),dimension(:),pointer :: Amat12_compr
        real(kind=8),dimension(:),allocatable :: Amat21_compr, Amat12_seq, Amat21_seq, tmpmat
        integer,parameter :: SPARSE=1
        integer,parameter :: DENSE=2
        real(kind=8) :: ex, max_error_p, mean_error_p, tt1, tt2
        real(kind=8),dimension(:),allocatable :: factor_arr
        real(kind=8),dimension(:),allocatable :: ovrlp_largep_new, invovrlpp_new
        real(kind=8),dimension(:),allocatable :: Amat12p_new, Amat21p_new
        real(kind=8),dimension(:),pointer :: Amat11p_new, Amat22p_new
        real(kind=8),dimension(:),allocatable :: rpower
      
      
        !!write(*,*) 'iorder',iorder

      
        call f_routine(id='overlapPowerGeneral')
        call timing(iproc,'lovrlp^-1     ','ON')
      
        ! several calculations at one are at the moment only possible for sparse exact, Taylor or ICE
        if (ncalc>1) then
            if (imode/=SPARSE .or. iorder<0) stop 'non-compliant arguments for ncalc>0'
        end if
      
        ! Indicate which power is calculated
        do icalc=1,ncalc
            select case (power(icalc))
            case (-2)
                inv_ovrlp_mat(icalc)%power = -0.5d0
            case (2)
                inv_ovrlp_mat(icalc)%power = 0.5d0
            case (1)
                inv_ovrlp_mat(icalc)%power = -1.d0
            case default
                stop 'wrong value of power(icalc)'
            end select
        end do
      
        ! Overwrite the default spin value is present. Usefull to manipulate spinless
        ! matrices even in a polarized calculation
        if (present(nspinx)) then
            nspin=nspinx
        else
            nspin=ovrlp_smat%nspin
        end if


        if (iproc==0) then
            call yaml_newline()
            call yaml_mapping_open('calculate S^x')
            if (imode==SPARSE) then
                call yaml_map('mode','sparse')
            else if (imode==DENSE) then
                call yaml_map('mode','dense')
            end if
            !call yaml_map('power(1)',power(1))
            call yaml_mapping_open('powers')
            do icalc=1,ncalc
                select case (power(icalc))
                case (-2)
                    call yaml_map('x','-1/2')
                case (2)
                    call yaml_map('x','1/2')
                case (1)
                    call yaml_map('x','-1')
                case default
                    stop 'wrong power(icalc)'
                end select
            end do
            call yaml_mapping_close()
            call yaml_map('order',iorder)
        end if
      
      
        ! Perform a check of the arguments
      
        if (imode/=SPARSE .and. imode/=DENSE) stop 'wrong imode'
      
        if (imode==DENSE) then
            if (.not.associated(ovrlp_mat%matrix)) stop 'ovrlp_mat%matrix not associated'
            if (.not.associated(inv_ovrlp_mat(1)%matrix)) stop 'inv_ovrlp_mat(1)%matrix not associated'
        end if
        
        if (check_accur) then
            if (.not.present(max_error)) stop 'max_error not present'
            if (.not.present(mean_error)) stop 'mean_error not present'
        end if
      
        if (power(1)/=-2 .and. power(1)/=1 .and. power(1)/=2) stop 'wrong value of power(1)'
      
        if (nproc/=1 .and. nproc/=bigdft_mpi%nproc) stop 'wrong value of nproc'
      
        ! Decide whether this routine is called in parallel or in serial.
        ! If parallel, take the default values from orbs, otherwise adjust them.
        !if (nproc>1) then
            norbp=ovrlp_smat%nfvctrp
            isorb=ovrlp_smat%isfvctr
        !else
        !    norbp=norb
        !    isorb=0
        !end if
      
        sparse_dense: if (imode==DENSE) then
            if (iorder==0) then
                call vcopy(ovrlp_smat%nfvctr*ovrlp_smat%nfvctr*nspin,ovrlp_mat%matrix(1,1,1),1,inv_ovrlp_mat(1)%matrix(1,1,1),1)
                if (power(1)==1) then
                   if (blocksize<0) then
                       do ispin=1,nspin
                           call overlap_minus_one_exact_serial(ovrlp_smat%nfvctr,inv_ovrlp_mat(1)%matrix(1:,1:,ispin))
                       end do
                   else
                      stop 'check if working - upper half may not be filled'
                      call dpotrf_parallel(iproc, nproc, blocksize, bigdft_mpi%mpi_comm, 'l', &
                           ovrlp_smat%nfvctr, inv_ovrlp_mat(1)%matrix(1:,1:,1), ovrlp_smat%nfvctr)
                      call dpotri_parallel(iproc, nproc, blocksize, bigdft_mpi%mpi_comm, 'l', &
                           ovrlp_smat%nfvctr, inv_ovrlp_mat(1)%matrix(1:,1:,1), ovrlp_smat%nfvctr)
                   end if
                else if (power(1)==2) then
                    do ispin=1,nspin
                        call overlap_plus_minus_one_half_exact(bigdft_mpi%nproc,ovrlp_smat%nfvctr, &
                             blocksize,.true.,inv_ovrlp_mat(1)%matrix(1:,1:,ispin),inv_ovrlp_smat)
                    end do
                else if (power(1)==-2) then
                    do ispin=1,nspin 
                        call overlap_plus_minus_one_half_exact(bigdft_mpi%nproc,ovrlp_smat%nfvctr, &
                             blocksize,.false.,inv_ovrlp_mat(1)%matrix(1:,1:,ispin),inv_ovrlp_smat)
                    end do
                end if
            else if (iorder<0) then
                ! sign approach as used in CP2K
                ! use 4 submatrices
                !if (nproc>1) then
                    !Amat12p = sparsematrix_malloc(inv_ovrlp_smat, iaction=DENSE_PARALLEL, id='Amat12p')
                    !Amat21p = sparsematrix_malloc(inv_ovrlp_smat, iaction=DENSE_PARALLEL, id='Amat21p')
                    !Amat11p = sparsematrix_malloc_ptr(inv_ovrlp_smat, iaction=DENSE_PARALLEL, id='Amat11p')
                !else
                    Amat12p = sparsematrix_malloc(ovrlp_smat,iaction=DENSE_PARALLEL,id='Amat12p')
                    Amat21p = sparsematrix_malloc(ovrlp_smat,iaction=DENSE_PARALLEL,id='Amat21p')
                    Amat11p = sparsematrix_malloc_ptr(ovrlp_smat,iaction=DENSE_PARALLEL,id='Amat11p')
                !end if
                ! save some memory but keep code clear - Amat22 and Amat11 should be identical as only combining S and I
                Amat22p=>Amat11p
                !if (nproc>1) then
                !    Amat21=sparsematrix_malloc0(inv_ovrlp_smat, iaction=DENSE_FULL, id='Amat21')
                !else
                    Amat21=sparsematrix_malloc0(ovrlp_smat,iaction=DENSE_FULL,id='Amat21')
                !end if
      
                do ispin=1,nspin
      
                    Amat12=>inv_ovrlp_mat(1)%matrix(:,:,ispin)
      
                    call vcopy(ovrlp_smat%nfvctr*ovrlp_smat%nfvctr,ovrlp_mat%matrix(1,1,ispin),1,Amat12(1,1),1)
                    do iorb=1,ovrlp_smat%nfvctr
                        Amat21(iorb,iorb)=1.0d0
                    end do
      
                    ! calculate Xn+1=0.5*Xn*(3I-Xn**2)
                    do its=1,abs(iorder)
                        if (norbp>0) call dgemm('n', 'n', ovrlp_smat%nfvctr, norbp, ovrlp_smat%nfvctr, -0.5d0, Amat12(1,1), &
                             ovrlp_smat%nfvctr, Amat21(1,isorb+1), ovrlp_smat%nfvctr, 0.0d0, Amat11p(1,1), ovrlp_smat%nfvctr)
                        !call dgemm('n', 'n', ovrlp_smat%nfvctr, norbp, ovrlp_smat%nfvctr, -0.5d0, Amat21(1,1), &
                        !     ovrlp_smat%nfvctr, Amat12(1,isorb+1), ovrlp_smat%nfvctr, 0.0d0, Amat22p(1,1), ovrlp_smat%nfvctr)
                        do iorb=1,norbp
                            Amat11p(iorb+isorb,iorb)=Amat11p(iorb+isorb,iorb)+1.5d0
                        !    Amat22p(iorb+isorb,iorb)=Amat22p(iorb+isorb,iorb)+1.5d0
                        end do
                        if (norbp>0) call dgemm('n', 'n', ovrlp_smat%nfvctr, norbp, ovrlp_smat%nfvctr, 1.0d0, Amat12(1,1), &
                             ovrlp_smat%nfvctr, Amat22p(1,1), ovrlp_smat%nfvctr, 0.0d0, Amat12p(1,1), ovrlp_smat%nfvctr)
                        if (norbp>0) call dgemm('n', 'n', ovrlp_smat%nfvctr, norbp, ovrlp_smat%nfvctr, 1.0d0, Amat21(1,1), &
                             ovrlp_smat%nfvctr, Amat11p(1,1), ovrlp_smat%nfvctr, 0.0d0, Amat21p(1,1), ovrlp_smat%nfvctr)
                        if(nproc > 1) then
                            call timing(iproc,'lovrlp^-1     ','OF')
                            call timing(iproc,'lovrlp_comm   ','ON')
                            call mpi_allgatherv(Amat12p, ovrlp_smat%nfvctr*norbp, mpi_double_precision, Amat12, &
                                 ovrlp_smat%nfvctr*ovrlp_smat%nfvctr_par(:), ovrlp_smat%nfvctr*ovrlp_smat%isfvctr_par, &
                                 mpi_double_precision, bigdft_mpi%mpi_comm, ierr)
                            call mpi_allgatherv(Amat21p, ovrlp_smat%nfvctr*norbp, mpi_double_precision, Amat21, &
                                 ovrlp_smat%nfvctr*ovrlp_smat%nfvctr_par(:), ovrlp_smat%nfvctr*ovrlp_smat%isfvctr_par, &
                                 mpi_double_precision, bigdft_mpi%mpi_comm, ierr)
                            call timing(iproc,'lovrlp_comm   ','OF')
                            call timing(iproc,'lovrlp^-1     ','ON')
                        else
                            call vcopy(ovrlp_smat%nfvctr**2,Amat12p(1,1),1,Amat12(1,1),1)
                            call vcopy(ovrlp_smat%nfvctr**2,Amat21p(1,1),1,Amat21(1,1),1)
                        end if
                    end do
      
      
                    if (power(1)==1) then
                        if (norbp>0) call dgemm('n', 'n', ovrlp_smat%nfvctr, norbp, ovrlp_smat%nfvctr, 1.0d0, Amat21(1,1), &
                             ovrlp_smat%nfvctr, Amat21p(1,1), ovrlp_smat%nfvctr, 0.0d0, Amat12p(1,1), ovrlp_smat%nfvctr)
                        if (nproc>1) then
                            call timing(iproc,'lovrlp^-1     ','OF')
                            call timing(iproc,'lovrlp_comm   ','ON')
                            call mpi_allgatherv(Amat12p, ovrlp_smat%nfvctr*norbp, &
                                 mpi_double_precision, inv_ovrlp_mat(1)%matrix(1,1,ispin), &
                                 ovrlp_smat%nfvctr*ovrlp_smat%nfvctr_par(:), ovrlp_smat%nfvctr*ovrlp_smat%isfvctr_par, &
                                 mpi_double_precision, bigdft_mpi%mpi_comm, ierr)
                            call timing(iproc,'lovrlp_comm   ','OF')
                            call timing(iproc,'lovrlp^-1     ','ON')
                        else
                            call vcopy(ovrlp_smat%nfvctr**2, Amat12p(1,1), 1, inv_ovrlp_mat(1)%matrix(1,1,ispin), 1)
                        end if
                    !else if (power(1)==2) then
                    !   call vcopy(ovrlp_smat%nfvctr**2,Amat12(1,1),1,inv_ovrlp_mat(1)%matrix(1,1),1)
                    else if (power(1)==-2) then
                        call vcopy(ovrlp_smat%nfvctr**2,Amat21(1,1),1,inv_ovrlp_mat(1)%matrix(1,1,ispin),1)
                    end if
      
                end do
      
                nullify(Amat22p)
                call f_free_ptr(Amat11p)
                call f_free(Amat12p)
                call f_free(Amat21p)
                nullify(Amat12)
                call f_free(Amat21)
      
            else
                ! we're missing cholesky method here... to avoid going to ridiculously high order Taylor expansion instead subtract 1000
                iorder_taylor=iorder
                if (iorder>1000) iorder_taylor=iorder-1000
                if (iorder_taylor>1) then
                    if (nproc>1) then
                        ovrlpminone => inv_ovrlp_mat(1)%matrix(:,:,:)
                        !ovrlpminonep = sparsematrix_malloc_ptr(ovrlp_smat,iaction=DENSE_PARALLEL,id='ovrlpminonep')
                        ovrlpminonep = f_malloc_ptr((/ovrlp_smat%nfvctr,max(ovrlp_smat%nfvctrp,1),nspin/),id='ovrlpminonep')
                    else
                        ovrlpminone = sparsematrix_malloc_ptr(ovrlp_smat,iaction=DENSE_FULL,id='ovrlpminone')
                        ovrlpminonep => ovrlpminone
                    end if
      
                    do ispin=1,nspin
                        if (ovrlp_smat%nfvctrp>0) call matrix_minus_identity_dense(ovrlp_smat%nfvctr,&
                                          ovrlp_smat%isfvctr,ovrlp_smat%nfvctrp, &
                                          ovrlp_mat%matrix(1:,ovrlp_smat%isfvctr+1:,ispin),ovrlpminonep(1:,1:,ispin))
      
      
                        !!if (iproc==0) write(*,*) 'isorb, ovrlp_mat%matrix(1,isorb+1,ispin)',isorb, ovrlp_mat%matrix(1,isorb+1,ispin)
                        !!do iorb=1,norbp
                        !!    do jorb=1,ovrlp_smat%nfvctr
                        !!        write(2800+10*iproc+ispin,'(a,3i8,3es14.6)') 'ispin, iorb, jorb, vals', &
                        !!             ispin, iorb, jorb, ovrlpminonep(jorb,iorb,ispin), ovrlp_mat%matrix(jorb,isorb+iorb,ispin)
                        !!    end do
                        !!end do
      
      
                        if(nproc > 1) then
                            call timing(iproc,'lovrlp^-1     ','OF')
                            call timing(iproc,'lovrlp_comm   ','ON')
                            call mpi_allgatherv(ovrlpminonep(1,1,ispin), ovrlp_smat%nfvctr*norbp, &
                                 mpi_double_precision, ovrlpminone(1,1,ispin), &
                                 ovrlp_smat%nfvctr*ovrlp_smat%nfvctr_par(:), ovrlp_smat%nfvctr*ovrlp_smat%isfvctr_par, &
                                 mpi_double_precision, bigdft_mpi%mpi_comm, ierr)
                            call timing(iproc,'lovrlp_comm   ','OF')
                            call timing(iproc,'lovrlp^-1     ','ON')
                        end if
      
                    end do
      
                    ovrlppoweroldp = sparsematrix_malloc_ptr(ovrlp_smat,iaction=DENSE_PARALLEL,id='ovrlppoweroldp')
                    if (norbp>0) call vcopy(ovrlp_smat%nfvctr*norbp*nspin,ovrlpminonep(1,1,1),1,ovrlppoweroldp(1,1,1),1)
      
                    if (nproc>1) then
                        call f_free_ptr(ovrlpminonep)
                    else
                        nullify(ovrlpminonep)
                    end if
      
                end if
      
                ovrlppowerp = sparsematrix_malloc_ptr(ovrlp_smat,iaction=DENSE_PARALLEL,id='ovrlppowerp')

      
                do ispin=1,nspin
      
                    if (power(1)==1) then
                        factor=-1.0d0
                    else if (power(1)==2) then
                        factor=0.5d0
                    else if (power(1)==-2) then
                        factor=-0.5d0
                    end if
      
                    if (nproc>1) then
                        inv_ovrlpp = sparsematrix_malloc_ptr(ovrlp_smat,iaction=DENSE_PARALLEL,id='inv_ovrlpp')
                    else
                        inv_ovrlpp => inv_ovrlp_mat(1)%matrix(:,:,ispin)
                    end if
      
                    if (norbp>0) call first_order_taylor_dense(ovrlp_smat%nfvctr,isorb,norbp,power(1),&
                        ovrlp_mat%matrix(1:,isorb+1:,ispin),inv_ovrlpp)
                    !!do iorb=1,norbp
                    !!    do jorb=1,ovrlp_smat%nfvctr
                    !!        write(2900+10*iproc+ispin,'(a,3i8,3es14.6)') 'ispin, iorb, jorb, vals', &
                    !!             ispin, iorb, jorb, inv_ovrlpp(jorb,iorb), ovrlp_mat%matrix(jorb,isorb+iorb,ispin)
                    !!    end do
                    !!end do
      
                    do i=2,iorder_taylor
                        if (norbp>0) call dgemm('n', 'n', ovrlp_smat%nfvctr, norbp, ovrlp_smat%nfvctr, &
                                          1.d0, ovrlpminone(1,1,ispin), &
                                          ovrlp_smat%nfvctr, ovrlppoweroldp(1,1,ispin), ovrlp_smat%nfvctr, &
                                          0.d0, ovrlppowerp(1,1), ovrlp_smat%nfvctr)
                       factor=newfactor(power(1),i,factor)
                        !!do iorb=1,norbp
                        !!    do jorb=1,ovrlp_smat%nfvctr
                        !!        write(3000+10*iproc+ispin,'(a,3i8,3es14.6)') 'ispin, iorb, jorb, vals', &
                        !!             ispin, iorb, jorb, ovrlppowerp(jorb,iorb), inv_ovrlpp(jorb,iorb), ovrlp_mat%matrix(jorb,isorb+iorb,ispin)
                        !!    end do
                        !!end do
                        call daxpy(ovrlp_smat%nfvctr*norbp,factor,ovrlppowerp,1,inv_ovrlpp,1)
                        !!do iorb=1,norbp
                        !!    do jorb=1,ovrlp_smat%nfvctr
                        !!        write(3100+10*iproc+ispin,'(a,4i8,3es14.6)') 'ispin, i, iorb, jorb, vals', &
                        !!             ispin, i, iorb, jorb, factor, ovrlppowerp(jorb,iorb), inv_ovrlpp(jorb,iorb)
                        !!    end do
                        !!end do
                        !!if (iproc==0) write(*,'(a,2i8,es16.9)') 'ispin, i, sum(inv_ovrlpp)', ispin, i, sum(inv_ovrlpp)
                        if (i/=iorder_taylor.and.norbp>0) then
                            call vcopy(ovrlp_smat%nfvctr*norbp,ovrlppowerp(1,1),1,ovrlppoweroldp(1,1,ispin),1)
                        end if
                    end do
      
      
                    !!write(*,'(a,2i8,es15.6)') 'iproc, ispin, sum(inv_ovrlpp)', iproc, ispin, sum(inv_ovrlpp)
                    if(nproc > 1) then
                        call timing(iproc,'lovrlp^-1     ','OF')
                        call timing(iproc,'lovrlp_comm   ','ON')
                        call mpi_allgatherv(inv_ovrlpp, ovrlp_smat%nfvctr*norbp, mpi_double_precision, &
                             inv_ovrlp_mat(1)%matrix(1,1,ispin), &
                             ovrlp_smat%nfvctr*ovrlp_smat%nfvctr_par(:), ovrlp_smat%nfvctr*ovrlp_smat%isfvctr_par, &
                             mpi_double_precision, bigdft_mpi%mpi_comm, ierr)
                        call timing(iproc,'lovrlp_comm   ','OF')
                        call timing(iproc,'lovrlp^-1     ','ON')
                        call f_free_ptr(inv_ovrlpp)
                    end if
                    !!if (iproc==0) write(*,'(a,2i8,es15.6)') 'iproc, ispin, sum(inv_ovrlp_mat(1)%matrix(:,:,ispin))', iproc, ispin, sum(inv_ovrlp_mat(1)%matrix(:,:,ispin))
                end do
      
      
                call f_free_ptr(ovrlppowerp)
      
                if (iorder_taylor>1) then
                    if(nproc > 1) then
                        nullify(ovrlpminone)
                    else
                        call f_free_ptr(ovrlpminone)
                    end if
                end if
      
                if (iorder_taylor>1) then
                    call f_free_ptr(ovrlppoweroldp)
                else
                    nullify(inv_ovrlpp)
                end if
            end if

            if (check_accur) then
                do ispin=1,nspin
                    call check_accur_overlap_minus_one(iproc,nproc,ovrlp_smat%nfvctr,&
                         ovrlp_smat%nfvctrp,ovrlp_smat%isfvctr,power(1),&
                         ovrlp_mat%matrix(:,:,ispin),inv_ovrlp_mat(1)%matrix(:,:,ispin),ovrlp_smat,max_error,mean_error)
                    if (iproc==0) then
                        call yaml_newline()
                        if (nspin==1) then
                            call yaml_map('max / mean error',(/max_error,mean_error/),fmt='(es8.2)')
                        else
                            if (ispin==1) then
                                call yaml_map('spin up, max / mean error',(/max_error,mean_error/),fmt='(es8.2)')
                            else if (ispin==2) then
                                call yaml_map('spin down, max / mean error',(/max_error,mean_error/),fmt='(es8.2)')
                            end if
                        end if
                    end if
                end do
            end if
        else if (imode==SPARSE) then
            if (iorder==0) then
                ovrlp_local = sparsematrix_malloc_ptr(inv_ovrlp_smat, iaction=DENSE_FULL, id='ovrlp_local')
                inv_ovrlp_local = sparsematrix_malloc_ptr(inv_ovrlp_smat, iaction=DENSE_FULL, id='inv_ovrlp_local')
                do icalc=1,ncalc
                    call timing(iproc,'lovrlp^-1     ','OF')
                    !!tmpmat = sparsematrix_malloc(ovrlp_smat,iaction=SPARSE_FULL,id='tmpmat')
                    !!call gather_matrix_from_taskgroups(iproc, nproc, ovrlp_smat, ovrlp_mat%matrix_compr, tmpmat)
                    call uncompress_matrix2(iproc, nproc, ovrlp_smat, ovrlp_mat%matrix_compr, ovrlp_local)
                    !!call f_free(tmpmat)
                    call timing(iproc,'lovrlp^-1     ','ON')
                    do ispin=1,nspin
                        !!write(*,*) 'sum(ovrlp_local(:,:,ispin))',sum(ovrlp_local(:,:,ispin))
                        call vcopy(ovrlp_smat%nfvctr*ovrlp_smat%nfvctr,ovrlp_local(1,1,ispin),1,inv_ovrlp_local(1,1,ispin),1)
                        if (power(icalc)==1) then
                           if (blocksize<0) then
                              call overlap_minus_one_exact_serial(ovrlp_smat%nfvctr,inv_ovrlp_local(1:,1:,ispin))
                              !if (iproc==0) then
                              !    do i=1,ovrlp_smat%nfvctr
                              !        do j=1,ovrlp_smat%nfvctr
                              !            write(400,'(2(i6,1x),es19.12)') i, j, inv_ovrlp_local(i,j,ispin)
                              !        end do
                              !    end do
                              !end if
                           else
                              !stop 'check if working - upper half may not be filled'
                              call dpotrf_parallel(iproc, nproc, blocksize, bigdft_mpi%mpi_comm, 'l', &
                                   ovrlp_smat%nfvctr, inv_ovrlp_local(1:,1:,ispin), ovrlp_smat%nfvctr)
                              call dpotri_parallel(iproc, nproc, blocksize, bigdft_mpi%mpi_comm, 'l', &
                                   ovrlp_smat%nfvctr, inv_ovrlp_local(1:,1:,ispin), ovrlp_smat%nfvctr)
                              !fill upper half...
                              do i=1,ovrlp_smat%nfvctr
                                  do j=i,ovrlp_smat%nfvctr
                                      inv_ovrlp_local(i,j,ispin)=inv_ovrlp_local(j,i,ispin)
                                  end do
                              end do
                           end if
                        else if (power(icalc)==2) then
                            call overlap_plus_minus_one_half_exact(bigdft_mpi%nproc,ovrlp_smat%nfvctr, &
                                 blocksize,.true.,inv_ovrlp_local(1:,1:,ispin),inv_ovrlp_smat)
                        else if (power(icalc)==-2) then
                            call overlap_plus_minus_one_half_exact(bigdft_mpi%nproc,ovrlp_smat%nfvctr, &
                                 blocksize,.false.,inv_ovrlp_local(1:,1:,ispin),inv_ovrlp_smat)
                        end if
                        call timing(iproc,'lovrlp^-1     ','OF')
                        call timing(iproc,'lovrlp^-1     ','ON')
                        !!write(*,*) 'sum(inv_ovrlp_local(:,:,ispin))',sum(inv_ovrlp_local(:,:,ispin))
                    end do
                    !call compress_matrix(iproc, inv_ovrlp_smat, inmat=inv_ovrlp_local, outmat=inv_ovrlp_mat(icalc)%matrix_compr)
                    do ispin=1,nspin
                        ishift=(ispin-1)*inv_ovrlp_smat%nvctr
                        ishift2=(ispin-1)*inv_ovrlp_smat%nvctrp_tg
                        call compress_matrix_distributed_wrapper(iproc, nproc, inv_ovrlp_smat, layout=DENSE_PARALLEL, &
                             matrixp=&
                               inv_ovrlp_local(1:,inv_ovrlp_smat%isfvctr+1:inv_ovrlp_smat%isfvctr+inv_ovrlp_smat%nfvctrp,ispin), &
                             matrix_compr=inv_ovrlp_mat(icalc)%matrix_compr(ishift2+1:))
                    end do
                end do
                call f_free_ptr(ovrlp_local)
                call f_free_ptr(inv_ovrlp_local)
                ! #############################################################################
            else if (iorder<0) then ! could improve timing for checking, but for now just making sure it works
                ! use 4 submatrices
                !!if (nproc>0) then
                !!    Amat12p = sparsematrix_malloc(inv_ovrlp_smat, iaction=DENSE_MATMUL, id='Amat12p')
                !!    Amat21p = sparsematrix_malloc0(inv_ovrlp_smat, iaction=DENSE_MATMUL, id='Amat21p')
                !!    Amat11p = sparsematrix_malloc0_ptr(inv_ovrlp_smat, iaction=DENSE_MATMUL, id='Amat11p')
                !!else
                !!    Amat12p = sparsematrix_malloc(inv_ovrlp_smat, iaction=DENSE_FULL, id='Amat12p')
                !!    Amat21p = sparsematrix_malloc0(inv_ovrlp_smat, iaction=DENSE_FULL, id='Amat21p')
                !!    Amat11p = sparsematrix_malloc0_ptr(inv_ovrlp_smat, iaction=DENSE_FULL, id='Amat11p')
                !!end if
                Amat12p_new = f_malloc(inv_ovrlp_smat%smmm%nvctrp,id='Amat12p_new')
                Amat21p_new = f_malloc0(inv_ovrlp_smat%smmm%nvctrp,id='Amat21p_new')
                Amat11p_new = f_malloc0_ptr(inv_ovrlp_smat%smmm%nvctrp,id='Amat11p_new')
                ! save some memory but keep code clear - Amat22 and Amat11 should be identical as only combining S and I
                !!Amat22p=>Amat11p
                Amat22p_new=>Amat11p_new
                Amat12_compr=>inv_ovrlp_mat(1)%matrix_compr(1:)
      
                call transform_sparse_matrix_local(ovrlp_smat, inv_ovrlp_smat, 'small_to_large', &
                     smatrix_compr_in=ovrlp_mat%matrix_compr, lmatrix_compr_out=Amat12_compr)
                Amat12_seq = sparsematrix_malloc(inv_ovrlp_smat, iaction=SPARSEMM_SEQ, id='Amat12_seq')
                Amat21_seq = sparsematrix_malloc(inv_ovrlp_smat, iaction=SPARSEMM_SEQ, id='Amat21_seq')
                Amat21_compr = sparsematrix_malloc(inv_ovrlp_smat, iaction=SPARSE_FULL, id='Amat21_compr')
      
                do ispin=1,nspin
      
                    ishift=(ispin-1)*inv_ovrlp_smat%nvctr
                    ishift2=(ispin-1)*inv_ovrlp_smat%nvctrp_tg
      
                    call sequential_acces_matrix_fast2(inv_ovrlp_smat, &
                         Amat12_compr(ishift2+1:), Amat12_seq)
                    call timing(iproc,'lovrlp^-1     ','OF')
                    !!call uncompress_matrix_distributed2(iproc, inv_ovrlp_smat, DENSE_MATMUL, &
                    !!     Amat12_compr(ishift2+1:), Amat12p)
      
                    !SM: This might not work with taskgroups
                    if (inv_ovrlp_smat%ntaskgroup/=1) stop 'overlapPowerGeneral: inv_ovrlp_smat%ntaskgroup/=1'
                    call transform_sparsity_pattern(inv_ovrlp_smat%nfvctr, &
                         inv_ovrlp_smat%smmm%nvctrp_mm, inv_ovrlp_smat%smmm%isvctr_mm, &
                         inv_ovrlp_smat%nseg, inv_ovrlp_smat%keyv, inv_ovrlp_smat%keyg, &
                         inv_ovrlp_smat%smmm%line_and_column_mm, &
                         inv_ovrlp_smat%smmm%nvctrp, inv_ovrlp_smat%smmm%isvctr, &
                         inv_ovrlp_smat%smmm%nseg, inv_ovrlp_smat%smmm%keyv, inv_ovrlp_smat%smmm%keyg, &
                         inv_ovrlp_smat%smmm%istsegline, 'small_to_large', &
                         Amat12_compr(ishift2+inv_ovrlp_smat%smmm%isvctr_mm+1:), &
                         Amat12p_new)
                    call timing(iproc,'lovrlp^-1     ','ON')
      
                    !!do iorb=1,inv_ovrlp_smat%smmm%nfvctrp
                    !!    Amat21p(iorb+inv_ovrlp_smat%smmm%isfvctr,iorb)=1.0d0
                    !!end do
                    do iorb=1,inv_ovrlp_smat%smmm%nvctrp
                        ii = inv_ovrlp_smat%smmm%isvctr + iorb
                        !!call get_line_and_column(ii, inv_ovrlp_smat%smmm%nseg, inv_ovrlp_smat%smmm%keyv, &
                        !!     inv_ovrlp_smat%smmm%keyg, iline, icolumn)
                        iline = inv_ovrlp_smat%smmm%line_and_column(1,iorb)
                        icolumn = inv_ovrlp_smat%smmm%line_and_column(2,iorb)
                        if (iline==icolumn) then
                            !write(*,*) 'iorb, ii, iline, icolumn', iorb, ii, iline, icolumn
                            Amat21p_new(iorb)=1.0d0
                        end if
                    end do
                    !write(*,*) 'init: sum(Amat21p_new)', sum(Amat21p_new)
                    call timing(iproc,'lovrlp^-1     ','OF')
                    !call compress_matrix_distributed(iproc, nproc, inv_ovrlp_smat, DENSE_MATMUL, &
                    !     Amat21p, Amat21_compr(inv_ovrlp_smat%isvctrp_tg+1:))
                    call compress_matrix_distributed_wrapper(iproc, nproc, inv_ovrlp_smat, SPARSE_MATMUL_LARGE, &
                         Amat21p_new, Amat21_compr(inv_ovrlp_smat%isvctrp_tg+1:))
                    !!write(*,*) 'after compr, sum(Amat21_compr)', sum(Amat21_compr)
                    call timing(iproc,'lovrlp^-1     ','ON')
                    call sequential_acces_matrix_fast(inv_ovrlp_smat, Amat21_compr, Amat21_seq)
      
                    ! calculate Xn+1=0.5*Xn*(3I-Xn**2)
                    do its=1,abs(iorder)
                        call timing(iproc,'lovrlp^-1     ','OF')
                        !call sparsemm(inv_ovrlp_smat, Amat12_seq, Amat21p, Amat11p)
                        !!write(*,*) 'sum(Amat21p_new)', sum(Amat21p_new)
                        call sparsemm_new(inv_ovrlp_smat, Amat12_seq, Amat21p_new, Amat11p_new)
                        call timing(iproc,'lovrlp^-1     ','ON')
                        !!write(*,*) 'after matmul: sum(Amat11p_new)', sum(Amat11p_new)
      
                        if (inv_ovrlp_smat%smmm%nvctrp>0) then
                            !call vscal(inv_ovrlp_smat%nfvctr*inv_ovrlp_smat%smmm%nfvctrp,-0.5d0,Amat11p(1,1),1)
                            call vscal(inv_ovrlp_smat%smmm%nvctrp,-0.5d0,Amat11p_new(1),1)
                        end if
                        !!write(*,*) 'after scal: sum(Amat11p_new)', sum(Amat11p_new)
                        !call vscal(ovrlp_smat%nfvctr*norbp,-0.5d0,Amat22p(1,1),1)
                        !!do iorb=1,inv_ovrlp_smat%smmm%nfvctrp
                        !!    Amat11p(iorb+inv_ovrlp_smat%smmm%isfvctr,iorb)=Amat11p(iorb+inv_ovrlp_smat%smmm%isfvctr,iorb)+1.5d0
                        !!!    Amat22p(iorb+isorb,iorb)=Amat22p(iorb+isorb,iorb)+1.5d0
                        !!end do
                        do iorb=1,inv_ovrlp_smat%smmm%nvctrp
                            ii = inv_ovrlp_smat%smmm%isvctr + iorb
                            !!call get_line_and_column(ii, inv_ovrlp_smat%smmm%nseg, &
                            !!     inv_ovrlp_smat%smmm%keyv, inv_ovrlp_smat%smmm%keyg, iline, icolumn)
                            iline = inv_ovrlp_smat%smmm%line_and_column(1,iorb)
                            icolumn = inv_ovrlp_smat%smmm%line_and_column(2,iorb)
                            if (iline==icolumn) then
                                Amat11p_new(iorb)=Amat11p_new(iorb)+1.5d0
                            end if
                        !    Amat22p(iorb+isorb,iorb)=Amat22p(iorb+isorb,iorb)+1.5d0
                        end do
      
                        call timing(iproc,'lovrlp^-1     ','OF')
                        !!call sparsemm(inv_ovrlp_smat, Amat12_seq, Amat22p, Amat12p)
                        !!call sparsemm(inv_ovrlp_smat, Amat21_seq, Amat11p, Amat21p)
                        !!write(*,*) 'sum(Amat11p_new)',sum(Amat11p_new)
                        call sparsemm_new(inv_ovrlp_smat, Amat12_seq, Amat22p_new, Amat12p_new)
                        call sparsemm_new(inv_ovrlp_smat, Amat21_seq, Amat11p_new, Amat21p_new)
                        call timing(iproc,'lovrlp^-1     ','ON')
      
                        if (its/=abs(iorder).or.power(1)/=2) then
                            call timing(iproc,'lovrlp^-1     ','OF')
                            !!call compress_matrix_distributed(iproc, nproc, inv_ovrlp_smat, DENSE_MATMUL, &
                            !!     Amat21p, Amat21_compr(inv_ovrlp_smat%isvctrp_tg+1:))
                            !!write(*,*) 'sum(Amat21p_new)',sum(Amat21p_new)
                            call compress_matrix_distributed_wrapper(iproc, nproc, inv_ovrlp_smat, SPARSE_MATMUL_LARGE, &
                                 Amat21p_new, Amat21_compr(inv_ovrlp_smat%isvctrp_tg+1:))
                            call timing(iproc,'lovrlp^-1     ','ON')
                        end if
                        if (its/=abs(iorder).or.power(1)==1) then
                            call sequential_acces_matrix_fast(inv_ovrlp_smat, Amat21_compr, Amat21_seq)
                        end if
                        if (its/=abs(iorder).or.power(1)==2) then
                            call timing(iproc,'lovrlp^-1     ','OF')
                            !!call compress_matrix_distributed(iproc, nproc, inv_ovrlp_smat, DENSE_MATMUL, &
                            !!     Amat12p, Amat12_compr)
                            call compress_matrix_distributed_wrapper(iproc, nproc, inv_ovrlp_smat, SPARSE_MATMUL_LARGE, &
                                 Amat12p_new, Amat12_compr)
                            call timing(iproc,'lovrlp^-1     ','ON')
                        end if
                        if (its/=abs(iorder)) then
                            call sequential_acces_matrix_fast2(inv_ovrlp_smat, Amat12_compr, Amat12_seq)
                        end if
                    end do
      
      
                    if (power(1)==1) then
                        call timing(iproc,'lovrlp^-1     ','OF')
                        !call sparsemm(inv_ovrlp_smat, Amat21_seq, Amat21p, Amat12p)
                        call sparsemm_new(inv_ovrlp_smat, Amat21_seq, Amat21p_new, Amat12p_new)
                        !!call compress_matrix_distributed(iproc, nproc, inv_ovrlp_smat, DENSE_MATMUL, Amat12p, &
                        !!     inv_ovrlp_mat(1)%matrix_compr(ishift2+1:))
                        call compress_matrix_distributed_wrapper(iproc, nproc, inv_ovrlp_smat, SPARSE_MATMUL_LARGE, Amat12p_new, &
                             inv_ovrlp_mat(1)%matrix_compr(ishift2+1:))
                        call timing(iproc,'lovrlp^-1     ','ON')
                    !else if (power(1)==2) then
                    !    call vcopy(inv_ovrlp_smat%nvctr,Amat12_compr(1),1,inv_ovrlp_smat%matrix_compr(1),1)
                    else if (power(1)==-2) then
                        call vcopy(inv_ovrlp_smat%nvctrp_tg,Amat21_compr(inv_ovrlp_smat%isvctrp_tg+1),1,&
                             inv_ovrlp_mat(1)%matrix_compr(ishift2+1),1)
                        !!write(*,*) 'sum(inv_ovrlp_mat(1)%matrix_compr)',sum(inv_ovrlp_mat(1)%matrix_compr)
                    end if
      
                end do
      
                call f_free(Amat12_seq)
                nullify(Amat22p)
                !!call f_free_ptr(Amat11p)
      
                nullify(Amat12_compr)
                call f_free(Amat21_compr)
                !!call f_free(Amat12p)
                !!call f_free(Amat21p)
                call f_free(Amat21_seq)
                call f_free(Amat12p_new)
                call f_free(Amat21p_new)
                call f_free_ptr(Amat11p_new)
      
      
            else
                if (iorder<1000) then
                    ovrlp_large_compr = sparsematrix_malloc(inv_ovrlp_smat, iaction=SPARSE_FULL, id='ovrlp_large_compr')
                    ! This is a bit quick and dirty
                    tmparr = sparsematrix_malloc(ovrlp_smat,iaction=SPARSE_FULL,id='tmparr')
                    call gather_matrix_from_taskgroups(iproc, nproc, ovrlp_smat, ovrlp_mat%matrix_compr, tmparr)
                    call transform_sparse_matrix(ovrlp_smat, inv_ovrlp_smat, &
                         tmparr, ovrlp_large_compr, 'small_to_large')
                    !!write(500+bigdft_mpi%iproc,'(a,2es16.8)') 'tmparr, large', tmparr(1), ovrlp_large_compr(1)
                    call f_free(tmparr)
      
                    !ovrlpminonep = sparsematrix_malloc_ptr(inv_ovrlp_smat, iaction=DENSE_MATMUL, id='ovrlpminonep')
                    ovrlpminonep_new = f_malloc_ptr(inv_ovrlp_smat%smmm%nvctrp,id='ovrlpminonep_new')
                    !invovrlpp_arr = f_malloc((/inv_ovrlp_smat%nfvctr,inv_ovrlp_smat%smmm%nfvctrp,ncalc/),id='invovrlpp_arr')
                    invovrlpp_arr_new = f_malloc((/inv_ovrlp_smat%smmm%nvctrp,ncalc/),id='invovrlpp_arr_new')
      
                    if (iorder>1) then
                        ovrlpminone_sparse_seq = sparsematrix_malloc(inv_ovrlp_smat, iaction=SPARSEMM_SEQ, &
                             id='ovrlpminone_sparse_seq')
                        ovrlpminone_sparse = sparsematrix_malloc_ptr(inv_ovrlp_smat, iaction=SPARSE_FULL, &
                             id='ovrlpminone_sparse')
                        !ovrlpminoneoldp = sparsematrix_malloc(inv_ovrlp_smat, iaction=DENSE_MATMUL, id='ovrlpminoneoldp')
                        ovrlpminoneoldp_new = f_malloc(inv_ovrlp_smat%smmm%nvctrp, id='ovrlpminoneoldp_new')
                    end if
      
                    factor_arr = f_malloc(ncalc,id='factor_arr')
      
                    do ispin=1,nspin
      
                        isshift=(ispin-1)*ovrlp_smat%nvctr
                        ilshift=(ispin-1)*inv_ovrlp_smat%nvctr
                        ilshift2=(ispin-1)*inv_ovrlp_smat%nvctrp_tg
      
                        if (iorder>1) then
                            !!ovrlpminone_sparse_seq = sparsematrix_malloc(inv_ovrlp_smat, iaction=SPARSEMM_SEQ, &
                            !!     id='ovrlpminone_sparse_seq')
                            !!ovrlpminone_sparse = sparsematrix_malloc_ptr(inv_ovrlp_smat, iaction=SPARSE_FULL, &
                            !!     id='ovrlpminone_sparse')
                            !!ovrlpminoneoldp = sparsematrix_malloc(inv_ovrlp_smat, iaction=DENSE_PARALLEL, id='ovrlpminoneoldp')
      
                            call matrix_minus_identity_sparse(ovrlp_smat%nfvctr, inv_ovrlp_smat, &
                                 ovrlp_large_compr(isshift+1), ovrlpminone_sparse)
                            !!write(500+bigdft_mpi%iproc,'(a,2es16.8)') 'large, sparse', ovrlp_large_compr(1), ovrlpminone_sparse(1)
                            call sequential_acces_matrix_fast(inv_ovrlp_smat, ovrlpminone_sparse, ovrlpminone_sparse_seq)
                            call timing(iproc,'lovrlp^-1     ','OF')
                            !call uncompress_matrix_distributed(iproc, inv_ovrlp_smat, DENSE_MATMUL, ovrlpminone_sparse, ovrlpminoneoldp)
                            call transform_sparsity_pattern(inv_ovrlp_smat%nfvctr, &
                                 inv_ovrlp_smat%smmm%nvctrp_mm, inv_ovrlp_smat%smmm%isvctr_mm, &
                                 inv_ovrlp_smat%nseg, inv_ovrlp_smat%keyv, inv_ovrlp_smat%keyg, &
                                 inv_ovrlp_smat%smmm%line_and_column_mm, &
                                 inv_ovrlp_smat%smmm%nvctrp, inv_ovrlp_smat%smmm%isvctr, &
                                 inv_ovrlp_smat%smmm%nseg, inv_ovrlp_smat%smmm%keyv, inv_ovrlp_smat%smmm%keyg, &
                                 inv_ovrlp_smat%smmm%istsegline, 'small_to_large', &
                                 ovrlpminone_sparse(ilshift2+inv_ovrlp_smat%smmm%isvctr_mm+1:), &
                                 ovrlpminoneoldp_new)
                            call timing(iproc,'lovrlp^-1     ','ON')
      
                            !!call f_free_ptr(ovrlpminone_sparse)
      
                            do icalc=1,ncalc
                                select case (power(icalc))
                                case(1)
                                    factor_arr(icalc)=-1.0d0
                                case(2)
                                    factor_arr(icalc)=0.5d0
                                case(-2)
                                    factor_arr(icalc)=-0.5d0
                                end select
                            end do
                        end if
      
      
                        if (inv_ovrlp_smat%smmm%nfvctrp>0) then
                            call timing(iproc,'lovrlp^-1     ','OF')
                            !!call uncompress_matrix_distributed(iproc, inv_ovrlp_smat, DENSE_MATMUL, &
                            !!     ovrlp_large_compr, ovrlpminonep(:,:,1))
                            !!write(500+bigdft_mpi%iproc,'(a,2es16.8)') 'BEF large, new', &
                            !!    ovrlp_large_compr(ilshift2+inv_ovrlp_smat%smmm%isvctr_mm+1), ovrlpminonep_new(1)
                            call transform_sparsity_pattern(inv_ovrlp_smat%nfvctr, &
                                 inv_ovrlp_smat%smmm%nvctrp_mm, inv_ovrlp_smat%smmm%isvctr_mm, &
                                 inv_ovrlp_smat%nseg, inv_ovrlp_smat%keyv, inv_ovrlp_smat%keyg, &
                                 inv_ovrlp_smat%smmm%line_and_column_mm, &
                                 inv_ovrlp_smat%smmm%nvctrp, inv_ovrlp_smat%smmm%isvctr, &
                                 inv_ovrlp_smat%smmm%nseg, inv_ovrlp_smat%smmm%keyv, inv_ovrlp_smat%smmm%keyg, &
                                 inv_ovrlp_smat%smmm%istsegline, 'small_to_large', &
                                 ovrlp_large_compr(ilshift2+inv_ovrlp_smat%smmm%isvctr_mm+1), &
                                 ovrlpminonep_new)
                            !!write(500+bigdft_mpi%iproc,'(a,2es16.8)') 'large, new', &
                            !!    ovrlp_large_compr(ilshift2+inv_ovrlp_smat%smmm%isvctr_mm+1), ovrlpminonep_new(1)
                            call timing(iproc,'lovrlp^-1     ','ON')
                            !!if (.not.check_accur) call f_free(ovrlp_large_compr)
                            do icalc=1,ncalc
                                !!call first_order_taylor_dense(inv_ovrlp_smat%nfvctr,inv_ovrlp_smat%smmm%isfvctr, &
                                !!     inv_ovrlp_smat%smmm%nfvctrp,power(icalc),ovrlpminonep,invovrlpp_arr(1,1,icalc))
                                call first_order_taylor_sparse_new(power(icalc), inv_ovrlp_smat, &
                                     ovrlpminonep_new, invovrlpp_arr_new(1,icalc))
                                !!write(500+bigdft_mpi%iproc,'(a,2es16.8)') 'FIRST ovrlpminonep_new(1), invovrlpp_arr_new(1,icalc)', ovrlpminonep_new(1), invovrlpp_arr_new(1,icalc)
                            end do
                        end if
                        call f_free(ovrlp_large_compr)
      
                        do i=2,iorder
                            call timing(iproc,'lovrlp^-1     ','OF')
                            !!call sparsemm(inv_ovrlp_smat, ovrlpminone_sparse_seq, ovrlpminoneoldp, ovrlpminonep)
                            !!write(500+bigdft_mpi%iproc,*) 'ovrlpminone_sparse_seq(1)', ovrlpminone_sparse_seq(1)
                            call sparsemm_new(inv_ovrlp_smat, ovrlpminone_sparse_seq, ovrlpminoneoldp_new, ovrlpminonep_new)
                            !!write(500+bigdft_mpi%iproc,'(a,2es16.8)') 'ovrlpminoneoldp_new(1), ovrlpminonep_new(1)', ovrlpminoneoldp_new(1), ovrlpminonep_new(1)
                            call timing(iproc,'lovrlp^-1     ','ON')
                            do icalc=1,ncalc
                                factor_arr(icalc)=newfactor(power(icalc),i,factor_arr(icalc))
                                if (inv_ovrlp_smat%smmm%nvctrp>0) then
                                    !!call daxpy(inv_ovrlp_smat%nfvctr*inv_ovrlp_smat%smmm%nfvctrp,factor_arr(icalc), &
                                    !!     ovrlpminonep,1,invovrlpp_arr(1,1,icalc),1)
                                    call daxpy(inv_ovrlp_smat%smmm%nvctrp,factor_arr(icalc), &
                                         ovrlpminonep_new,1,invovrlpp_arr_new(1,icalc),1)
                                    !!write(500+bigdft_mpi%iproc,'(a,2es16.8)') 'ovrlpminonep_new(1), invovrlpp_arr_new(1,icalc)', ovrlpminonep_new(1), invovrlpp_arr_new(1,icalc)
                                end if
                            end do
                            if (i/=iorder.and.inv_ovrlp_smat%smmm%nvctrp>0) then
                                !!call vcopy(inv_ovrlp_smat%nfvctr*inv_ovrlp_smat%smmm%nfvctrp,&
                                !!ovrlpminonep(1,1,1),1,ovrlpminoneoldp(1,1),1)
                                call vcopy(inv_ovrlp_smat%smmm%nvctrp,&
                                ovrlpminonep_new(1),1,ovrlpminoneoldp_new(1),1)
                            end if
                        end do
                        !!call to_zero(inv_ovrlp_smat%nvctr, inv_ovrlp_smat%matrix_compr(1))
                        call timing(iproc,'lovrlp^-1     ','OF')
                        do icalc=1,ncalc
                            !!call compress_matrix_distributed(iproc, nproc, inv_ovrlp_smat, &
                            !!     DENSE_MATMUL, invovrlpp_arr(1:,1:,icalc), &
                            !!     inv_ovrlp_mat(icalc)%matrix_compr(ilshift2+1:))
                            call compress_matrix_distributed_wrapper(iproc, nproc, inv_ovrlp_smat, &
                                 SPARSE_MATMUL_LARGE, invovrlpp_arr_new(1:,icalc), &
                                 inv_ovrlp_mat(icalc)%matrix_compr(ilshift2+1:))
                            !!write(500+bigdft_mpi%iproc,'(a,2es16.8)') 'inv_mat, inv_compr', invovrlpp_arr_new(1,icalc), inv_ovrlp_mat(icalc)%matrix_compr(ilshift2+1)
                        end do
                        call timing(iproc,'lovrlp^-1     ','ON')
      
                    end do
      
                    call f_free(factor_arr)
      
                    if (iorder>1) then
                        call f_free(ovrlpminone_sparse_seq)
                        !!call f_free(ovrlpminoneoldp)
                        call f_free(ovrlpminoneoldp_new)
                        call f_free_ptr(ovrlpminone_sparse)
                    end if
      
                    !!call f_free(invovrlpp_arr)
                    call f_free(invovrlpp_arr_new)
                    !!call f_free_ptr(ovrlpminonep)
                    call f_free_ptr(ovrlpminonep_new)
      
                else
      
                    ! @ NEW: ICE ##########################
                    !!select case (power(1))
                    !!case (-2)
                    !!    ex=-0.5d0
                    !!case (2)
                    !!    ex=0.5d0
                    !!case (1)
                    !!    ex=-1.d0
                    !!case default
                    !!    stop 'wrong power(1)'
                    !!end select
                    !!tmpmat = sparsematrix_malloc(ovrlp_smat,iaction=SPARSE_FULL,id='tmpmat')
                    !!call vcopy(ovrlp_smat%nvctr, ovrlp_mat%matrix_compr(1), 1, tmpmat(1), 1)
                    !!call gather_matrix_from_taskgroups_inplace(iproc, nproc, ovrlp_smat, ovrlp_mat)
                    rpower = f_malloc(ncalc,id='rpower')
                    do icalc=1,ncalc
                        select case (power(icalc))
                        case (-2)
                            rpower(icalc) = -0.5d0
                        case (2)
                            rpower(icalc) = 0.5d0
                        case (1)
                            rpower(icalc) = -1.d0
                        case default
                            stop 'wrong value of power(icalc)'
                        end select
                    end do
                    call inverse_chebyshev_expansion(iproc, nproc, iorder-1000, &
                         ovrlp_smat, inv_ovrlp_smat, ncalc, rpower, ovrlp_mat, inv_ovrlp_mat)
                    call f_free(rpower)
                    !!call vcopy(ovrlp_smat%nvctr, tmpmat(1), 1, ovrlp_mat%matrix_compr(1), 1)
                    !!call f_free(tmpmat)
                    ! #####################################
                end if
            end if
      
            if (check_accur) then
                ! HERE STARTS LINEAR CHECK ##########################
                !!invovrlpp = sparsematrix_malloc(inv_ovrlp_smat, iaction=DENSE_MATMUL, id='invovrlpp')
                invovrlpp_new = f_malloc(inv_ovrlp_smat%smmm%nvctrp, id='invovrlpp_new')
                !!if (iorder<1 .or. iorder>=1000) then
                    ovrlp_large_compr = sparsematrix_malloc(inv_ovrlp_smat, iaction=SPARSE_TASKGROUP, id='ovrlp_large_compr')
                    call transform_sparse_matrix_local(ovrlp_smat, inv_ovrlp_smat, 'small_to_large', &
                         smatrix_compr_in=ovrlp_mat%matrix_compr, lmatrix_compr_out=ovrlp_large_compr)
                !!end if
                invovrlp_compr_seq = sparsematrix_malloc(inv_ovrlp_smat, iaction=SPARSEMM_SEQ, id='ovrlp_large_compr_seq')
                !!ovrlp_largep = sparsematrix_malloc(inv_ovrlp_smat, iaction=DENSE_MATMUL, id='ovrlp_largep')
                ovrlp_largep_new = f_malloc(inv_ovrlp_smat%smmm%nvctrp,id='ovrlp_largep')

                !!if (iproc==0) write(*,*) 'TEST ##############################################'
                !!do icalc=1,ncalc
                !!    resmat = sparsematrix_malloc(inv_ovrlp_smat, iaction=SPARSE_TASKGROUP, id='resmat')
                !!    call matrix_matrix_mult_wrapper(iproc, nproc, inv_ovrlp_smat, ovrlp_large_compr, inv_ovrlp_mat(icalc)%matrix_compr, resmat)
                !!    call deviation_from_unity_parallel_new(iproc, nproc, inv_ovrlp_smat%nfvctr, inv_ovrlp_smat%nfvctrp, inv_ovrlp_smat%isfvctr, &
                !!         resmat, inv_ovrlp_smat, tt1, tt2)
                !!    if (iproc==0) then
                !!        call yaml_map('NEW max mean error',(/tt1,tt2/))
                !!        do i=1,size(resmat,1)
                !!            write(500,*) i, resmat(i)
                !!        end do
                !!        do i=1,size(ovrlp_large_compr)
                !!            write(600,*) i, ovrlp_large_compr(i)
                !!        end do
                !!        do i=1,size(inv_ovrlp_mat(icalc)%matrix_compr)
                !!            write(700,*) i, inv_ovrlp_mat(icalc)%matrix_compr(i)
                !!        end do
                !!    end if
                !!    call f_free(resmat)
                !!end do
                !!if (iproc==0) write(*,*) 'END TEST ##########################################'
      
                if (iproc==0) then
                    call yaml_newline()
                    call yaml_sequence_open('error estimation')
                end if
                do icalc=1,ncalc
                    do ispin=1,nspin
                        isshift=(ispin-1)*ovrlp_smat%nvctr
                        ilshift=(ispin-1)*inv_ovrlp_smat%nvctrp_tg
                        ilshift2=(ispin-1)*inv_ovrlp_smat%nvctrp_tg
                        call timing(iproc,'lovrlp^-1     ','OF')
                        !!call uncompress_matrix_distributed2(iproc, inv_ovrlp_smat, &
                        !!     DENSE_MATMUL, ovrlp_large_compr(ilshift+1), ovrlp_largep)
                        if (inv_ovrlp_smat%smmm%nvctrp_mm>0) then
                            call transform_sparsity_pattern(inv_ovrlp_smat%nfvctr, &
                                 inv_ovrlp_smat%smmm%nvctrp_mm, inv_ovrlp_smat%smmm%isvctr_mm, &
                                 inv_ovrlp_smat%nseg, inv_ovrlp_smat%keyv, inv_ovrlp_smat%keyg, &
                                 inv_ovrlp_smat%smmm%line_and_column_mm, &
                                 inv_ovrlp_smat%smmm%nvctrp, inv_ovrlp_smat%smmm%isvctr, &
                                 inv_ovrlp_smat%smmm%nseg, inv_ovrlp_smat%smmm%keyv, inv_ovrlp_smat%smmm%keyg, &
                                 inv_ovrlp_smat%smmm%istsegline, 'small_to_large', &
                                 ovrlp_large_compr(ilshift+inv_ovrlp_smat%smmm%isvctr_mm-inv_ovrlp_smat%isvctrp_tg+1), &
                                 ovrlp_largep_new)
                        end if
      
                        call timing(iproc,'lovrlp^-1     ','ON')
                        call sequential_acces_matrix_fast2(inv_ovrlp_smat, &
                             inv_ovrlp_mat(icalc)%matrix_compr(ilshift2+1:), invovrlp_compr_seq)
                        !!write(500+bigdft_mpi%iproc,'(a,2es16.8)') 'inv, inv_seq', inv_ovrlp_mat(icalc)%matrix_compr(ilshift2+1), invovrlp_compr_seq(1)
                        !!write(*,*) 'sum(inv_ovrlp_mat(1)%matrix_compr(ilshift+1:ilshift+inv_ovrlp_smat%nvctr)', sum(inv_ovrlp_mat(1)%matrix_compr(ilshift+1:ilshift+inv_ovrlp_smat%nvctr))
      
                        if (power(icalc)==1) then
                            call check_accur_overlap_minus_one_sparse_new(iproc, nproc, inv_ovrlp_smat, ovrlp_smat%nfvctr, &
                                 inv_ovrlp_smat%smmm%nfvctrp, inv_ovrlp_smat%smmm%isfvctr, &
                                 inv_ovrlp_smat%smmm%nseq, inv_ovrlp_smat%smmm%nout, &
                                 invovrlp_compr_seq, ovrlp_largep_new, power(icalc), &
                                 max_error, mean_error)
                        else if (power(icalc)==2) then
                            call timing(iproc,'lovrlp^-1     ','OF')
                            !!call uncompress_matrix_distributed2(iproc, inv_ovrlp_smat, DENSE_MATMUL, &
                            !!     inv_ovrlp_mat(icalc)%matrix_compr(ilshift2+1:), invovrlpp)
                            if (inv_ovrlp_smat%smmm%nvctrp_mm>0) then
                                ist = ilshift2+inv_ovrlp_smat%smmm%isvctr_mm-inv_ovrlp_smat%isvctrp_tg+1
                                call transform_sparsity_pattern(inv_ovrlp_smat%nfvctr, &
                                     inv_ovrlp_smat%smmm%nvctrp_mm, inv_ovrlp_smat%smmm%isvctr_mm, &
                                     inv_ovrlp_smat%nseg, inv_ovrlp_smat%keyv, inv_ovrlp_smat%keyg, &
                                     inv_ovrlp_smat%smmm%line_and_column_mm, &
                                     inv_ovrlp_smat%smmm%nvctrp, inv_ovrlp_smat%smmm%isvctr, &
                                     inv_ovrlp_smat%smmm%nseg, inv_ovrlp_smat%smmm%keyv, inv_ovrlp_smat%smmm%keyg, &
                                     inv_ovrlp_smat%smmm%istsegline, 'small_to_large', &
                                     inv_ovrlp_mat(icalc)%matrix_compr(ist:), &
                                     invovrlpp_new)
                            end if
                            call timing(iproc,'lovrlp^-1     ','ON')
                            call check_accur_overlap_minus_one_sparse_new(iproc, nproc, inv_ovrlp_smat, ovrlp_smat%nfvctr, &
                                 inv_ovrlp_smat%smmm%nfvctrp, inv_ovrlp_smat%smmm%isfvctr, &
                                 inv_ovrlp_smat%smmm%nseq, inv_ovrlp_smat%smmm%nout, &
                                 invovrlp_compr_seq, invovrlpp_new, power(icalc), &
                                 max_error, mean_error, cmatp=ovrlp_largep_new)
                        else if (power(icalc)==-2) then
                            ovrlp_compr_seq = sparsematrix_malloc(inv_ovrlp_smat, iaction=SPARSEMM_SEQ, id='ovrlp_compr_seq') 
                            call sequential_acces_matrix_fast2(inv_ovrlp_smat, ovrlp_large_compr(ilshift+1), ovrlp_compr_seq)
                            call timing(iproc,'lovrlp^-1     ','OF')
                            !!call uncompress_matrix_distributed2(iproc, inv_ovrlp_smat, DENSE_MATMUL, &
                            !!     inv_ovrlp_mat(icalc)%matrix_compr(ilshift2+1:), invovrlpp)
                            if (inv_ovrlp_smat%smmm%nvctrp_mm>0) then
                                ist = ilshift2+inv_ovrlp_smat%smmm%isvctr_mm-inv_ovrlp_smat%isvctrp_tg+1
                                call transform_sparsity_pattern(inv_ovrlp_smat%nfvctr, &
                                     inv_ovrlp_smat%smmm%nvctrp_mm, inv_ovrlp_smat%smmm%isvctr_mm, &
                                     inv_ovrlp_smat%nseg, inv_ovrlp_smat%keyv, inv_ovrlp_smat%keyg, &
                                     inv_ovrlp_smat%smmm%line_and_column_mm, &
                                     inv_ovrlp_smat%smmm%nvctrp, inv_ovrlp_smat%smmm%isvctr, &
                                     inv_ovrlp_smat%smmm%nseg, inv_ovrlp_smat%smmm%keyv, inv_ovrlp_smat%smmm%keyg, &
                                     inv_ovrlp_smat%smmm%istsegline, 'small_to_large', &
                                     inv_ovrlp_mat(icalc)%matrix_compr(ist:), &
                                     invovrlpp_new)
                            end if
                            call timing(iproc,'lovrlp^-1     ','ON')
                            call check_accur_overlap_minus_one_sparse_new(iproc, nproc, inv_ovrlp_smat, ovrlp_smat%nfvctr, &
                                 inv_ovrlp_smat%smmm%nfvctrp, inv_ovrlp_smat%smmm%isfvctr, &
                                 inv_ovrlp_smat%smmm%nseq, inv_ovrlp_smat%smmm%nout, &
                                 invovrlp_compr_seq, invovrlpp_new, power(icalc), &
                                 max_error, mean_error, &
                                 ovrlp_compr_seq)
                            call f_free(ovrlp_compr_seq)
                        else
                            stop 'wrong power(icalc)'
                        end if
                        if (iproc==0) then
                            call yaml_newline()
                            if (nspin==1) then
                                call yaml_map('max / mean error',(/max_error,mean_error/),fmt='(es8.2)')
                            else
                                if (ispin==1) then
                                    call yaml_map('spin up, max / mean error',(/max_error,mean_error/),fmt='(es8.2)')
                                else if (ispin==2) then
                                    call yaml_map('spin down, max / mean error',(/max_error,mean_error/),fmt='(es8.2)')
                                end if
                            end if
                        end if
                    end do
                end do
                call yaml_sequence_close()
                call f_free(invovrlp_compr_seq)
                !!call f_free(ovrlp_largep)
                call f_free(ovrlp_largep_new)
                !!call f_free(invovrlpp)
                call f_free(invovrlpp_new)
                call f_free(ovrlp_large_compr)
                !HERE ENDS LINEAR CHECK #############################
            end if
        end if sparse_dense
      
        if (iproc==0) then
            call yaml_mapping_close()
            call yaml_newline()
        end if
      
        call timing(iproc,'lovrlp^-1     ','OF')
        call f_release_routine()

      
      end subroutine overlapPowerGeneral


      subroutine check_accur_overlap_minus_one_sparse_new(iproc, nproc, smat, norb, norbp, isorb, nseq, nout, &
                 amat_seq, bmatp, power, &
                 max_error, mean_error, dmat_seq, cmatp)
        use module_base
        use sparsematrix_base, only: sparse_matrix
        use sparsematrix, only: sparsemm_new
        implicit none
        integer,intent(in) :: iproc, nproc, norb, norbp, isorb, nseq, nout, power
        type(sparse_matrix) :: smat
        real(kind=8),dimension(nseq),intent(in) :: amat_seq
        real(kind=8),dimension(smat%smmm%nvctrp),intent(in) :: bmatp
        real(kind=8),intent(out) :: max_error, mean_error
        real(kind=8),dimension(nseq),intent(in),optional :: dmat_seq
        real(kind=8),dimension(smat%smmm%nvctrp),intent(in),optional :: cmatp
      
        real(kind=8), allocatable, dimension(:,:) :: tmp, tmp2
        real(kind=8), allocatable, dimension(:) :: tmpp, tmp2p
        integer :: ierr, i,j
      
        call f_routine(id='check_accur_overlap_minus_one_sparse_new')
      
        tmpp=f_malloc0((smat%smmm%nvctrp),id='tmpp')
        if (power==1) then
           !!call dgemm('n', 'n', norb, norbp, norb, 1.d0, inv_ovrlp(1,1), &
           !!     norb, ovrlp(1,isorb+1), norb, 0.d0, tmpp(1,1), norb)
           call sparsemm_new(smat, amat_seq, bmatp, tmpp)
           call deviation_from_unity_parallel_new(iproc, nproc, tmpp, smat, max_error, mean_error)
        else if (power==2) then
            if (.not.present(cmatp)) stop 'cmatp not present'
           !!call dgemm('n', 'n', norb, norbp, norb, 1.d0, inv_ovrlp(1,1), &
           !!     norb, inv_ovrlp(1,isorb+1), norb, 0.d0, tmpp(1,1), norb)
           !write(*,*) 'iproc, sum(amat_seq), sum(bmatp)', iproc, sum(amat_seq), sum(bmatp)
           call sparsemm_new(smat, amat_seq, bmatp, tmpp)
           call max_matrix_diff_parallel_new(iproc, norb, norbp, isorb, tmpp, cmatp, smat, max_error, mean_error)
           !max_error=0.5d0*max_error
           !mean_error=0.5d0*mean_error
        else if (power==-2) then
           if (.not.present(dmat_seq)) stop 'dmat_seq not present'
           !!call dgemm('n', 'n', norb, norbp, norb, 1.d0, inv_ovrlp(1,1), &
           !!     norb, inv_ovrlp(1,isorb+1), norb, 0.d0, tmpp(1,1), norb)
           !!write(500+bigdft_mpi%iproc,*) 'amat_seq(1)', amat_seq(1)
           call sparsemm_new(smat, amat_seq, bmatp, tmpp)
           tmp2p=f_malloc0(smat%smmm%nvctrp,id='tmp2p')
           !!call dgemm('n', 'n', norb, norbp, norb, 1.d0, ovrlp(1,1), &
           !!     norb, tmpp(1,1), norb, 0.d0, tmp2p(1,1), norb)
           call sparsemm_new(smat, dmat_seq, tmpp, tmp2p)
           call deviation_from_unity_parallel_new(iproc, nproc, tmp2p, smat, max_error, mean_error)
           !max_error=0.5d0*max_error
           !mean_error=0.5d0*mean_error
           call f_free(tmp2p)
        else
           stop 'Error in check_accur_overlap_minus_one_sparse_new'
        end if
        call f_free(tmpp)
      
        call f_release_routine()
      
      end subroutine check_accur_overlap_minus_one_sparse_new


      subroutine overlap_minus_one_half_serial(iproc, nproc, iorder, power, blocksize, &
                 norb, ovrlp_matrix, inv_ovrlp_matrix, check_accur, &
                 smat, max_error, mean_error)
        use module_base
        use sparsematrix_base, only: sparse_matrix
        use yaml_output
        implicit none
        
        ! Calling arguments
        integer,intent(in) :: iproc, nproc, iorder, blocksize, power, norb
        real(kind=8),dimension(norb,norb),intent(in) :: ovrlp_matrix
        real(kind=8),dimension(:,:),pointer,intent(inout) :: inv_ovrlp_matrix
        type(sparse_matrix),intent(in) :: smat
        logical,intent(in) :: check_accur
        real(kind=8),intent(out),optional :: max_error, mean_error
        
        ! Local variables
        integer :: iorb, jorb, info, iiorb, isorb, norbp, ii, ii_inv, iii, ierr, i, its, maxits
        integer :: matrixindex_in_compressed, nmaxvalk
        real(kind=8), dimension(:,:), pointer :: ovrlpminonep, ovrlpminone, inv_ovrlpp, ovrlppowerp, ovrlppoweroldp
        real(kind=8), dimension(:,:), pointer :: inv_ovrlp_half_tmp, ovrlp_local, inv_ovrlp_local
        real(kind=8) :: factor
        logical :: ovrlp_allocated, inv_ovrlp_allocated
      
        ! new for sparse taylor
        integer :: nout, nseq
        integer,dimension(:,:,:),allocatable :: istindexarr
        real(kind=8),dimension(:),pointer :: ovrlpminone_sparse
        real(kind=8),dimension(:),allocatable :: ovrlp_compr_seq, ovrlpminone_sparse_seq, ovrlp_large_compr
        real(kind=8),dimension(:),allocatable :: invovrlp_compr_seq
        real(kind=8),dimension(:,:),allocatable :: ovrlpminoneoldp, invovrlpp, ovrlp_largep
        real(kind=8),dimension(:,:),allocatable :: Amat12p, Amat21p, Amat21
        real(kind=8),dimension(:,:),pointer :: Amat12, Amat11p, Amat22p
        real(kind=8),dimension(:),pointer :: Amat12_compr
        real(kind=8),dimension(:),allocatable :: Amat21_compr, Amat12_seq, Amat21_seq
        integer,parameter :: SPARSE=1
        integer,parameter :: DENSE=2
      
      
        !!write(*,*) 'iorder',iorder
      
      
        call f_routine(id='overlap_minus_one_half_serial')
        call timing(iproc,'lovrlp^-1     ','ON')
      
        if (nproc>1) then
            stop 'this routine only works in serial'
        end if
      
        
        if (check_accur) then
            if (.not.present(max_error)) stop 'max_error not present'
            if (.not.present(mean_error)) stop 'mean_error not present'
        end if
      
        if (power/=-2 .and. power/=1 .and. power/=2) stop 'wrong value of power'
      
      
            if (iorder==0) then
                call vcopy(norb*norb,ovrlp_matrix(1,1),1,inv_ovrlp_matrix(1,1),1)
                if (power==1) then
                   call overlap_minus_one_exact_serial(norb,inv_ovrlp_matrix)
                else if (power==2) then
                   call overlap_plus_minus_one_half_exact(1,norb,blocksize,.true.,inv_ovrlp_matrix,smat)
                else if (power==-2) then
                   call overlap_plus_minus_one_half_exact(1,norb,blocksize,.false.,inv_ovrlp_matrix,smat)
                end if
            else if (iorder<0) then
                Amat12p = f_malloc((/norb,norb/), id='Amat12p')
                Amat21p = f_malloc((/norb,norb/), id='Amat21p')
                Amat11p = f_malloc_ptr((/norb,norb/), id='Amat11p')
                ! save some memory but keep code clear - Amat22 and Amat11 should be identical as only combining S and I
                Amat22p=>Amat11p
                Amat12=>inv_ovrlp_matrix
                Amat21=f_malloc0((/norb,norb/), id='Amat21')
      
                call vcopy(norb*norb,ovrlp_matrix(1,1),1,Amat12(1,1),1)
                do iorb=1,norb
                    Amat21(iorb,iorb)=1.0d0
                end do
      
                ! calculate Xn+1=0.5*Xn*(3I-Xn**2)
                do its=1,abs(iorder)
                    call dgemm('n', 'n', norb, norb, norb, -0.5d0, Amat12(1,1), &
                         norb, Amat21(1,1), norb, 0.0d0, Amat11p(1,1), norb)
                    do iorb=1,norb
                        Amat11p(iorb,iorb)=Amat11p(iorb,iorb)+1.5d0
                    end do
                    call dgemm('n', 'n', norb, norb, norb, 1.0d0, Amat12(1,1), &
                         norb, Amat22p(1,1), norb, 0.0d0, Amat12p(1,1), norb)
                    call dgemm('n', 'n', norb, norb, norb, 1.0d0, Amat21(1,1), &
                         norb, Amat11p(1,1), norb, 0.0d0, Amat21p(1,1), norb)
                    call vcopy(norb**2,Amat12p(1,1),1,Amat12(1,1),1)
                    call vcopy(norb**2,Amat21p(1,1),1,Amat21(1,1),1)
                end do
      
                nullify(Amat22p)
                call f_free_ptr(Amat11p)
      
                if (power==1) then
                    call dgemm('n', 'n', norb, norb, norb, 1.0d0, Amat21(1,1), &
                         norb, Amat21p(1,1), norb, 0.0d0, Amat12p(1,1), norb)
                    call vcopy(norb**2, Amat12p(1,1), 1, inv_ovrlp_matrix(1,1), 1)
                else if (power==-2) then
                    call vcopy(norb**2,Amat21(1,1),1,inv_ovrlp_matrix(1,1),1)
                end if
      
                call f_free(Amat12p)
                call f_free(Amat21p)
                nullify(Amat12)
                call f_free(Amat21)
      
            else
                if (iorder>1) then
                    ovrlpminone = f_malloc_ptr((/norb,norb/), id='ovrlpminone')
                    ovrlpminonep => ovrlpminone
      
                    call matrix_minus_identity_dense(norb,0,norb,ovrlp_matrix(1,1),ovrlpminonep)
      
                    ovrlppoweroldp = f_malloc_ptr((/norb,norb/), id='ovrlppoweroldp')
      
                    call vcopy(norb*norb,ovrlpminonep(1,1),1,ovrlppoweroldp(1,1),1)
      
                    nullify(ovrlpminonep)
      
                    ovrlppowerp = f_malloc_ptr((/norb,norb/), id='ovrlppowerp')
      
                    if (power==1) then
                        factor=-1.0d0
                    else if (power==2) then
                        factor=0.5d0
                    else if (power==-2) then
                        factor=-0.5d0
                    end if
                end if
      
                if (nproc>1) then
                    inv_ovrlpp = f_malloc_ptr((/norb,norb/), id='inv_ovrlpp')
                else
                    inv_ovrlpp => inv_ovrlp_matrix
                end if
      
                call first_order_taylor_dense(norb,0,norb,power,ovrlp_matrix(1,1),inv_ovrlpp)
      
                do i=2,iorder
                    call dgemm('n', 'n', norb, norb, norb, 1.d0, ovrlpminone(1,1), &
                         norb, ovrlppoweroldp(1,1), norb, 0.d0, ovrlppowerp(1,1), norb)
                    factor=newfactor(power,i,factor)
                    call daxpy(norb*norb,factor,ovrlppowerp,1,inv_ovrlpp,1)
                    if (i/=iorder) call vcopy(norb*norb,ovrlppowerp(1,1),1,ovrlppoweroldp(1,1),1)
                end do
      
                if (iorder>1) then
                    if(nproc > 1) then
                        nullify(ovrlpminone)
                    else
                        call f_free_ptr(ovrlpminone)
                    end if
      
                    call f_free_ptr(ovrlppowerp)
                    call f_free_ptr(ovrlppoweroldp)
      
                end if
      
                nullify(inv_ovrlpp)
            end if
      
            if (check_accur) then
                call check_accur_overlap_minus_one(iproc,nproc,norb,norb,0,power,ovrlp_matrix,inv_ovrlp_matrix,&
                     smat, max_error,mean_error)
            end if
      
      
        call f_release_routine()
        call timing(iproc,'lovrlp^-1     ','OF')
      
      
      end subroutine overlap_minus_one_half_serial


      subroutine check_accur_overlap_minus_one(iproc,nproc,norb,norbp,isorb,power,ovrlp,inv_ovrlp,&
                 smat,max_error,mean_error)
        use module_base
        use sparsematrix_base, only: sparse_matrix
        implicit none
        integer,intent(in) :: iproc, nproc, norb, norbp, isorb, power
        real(kind=8),dimension(norb,norb),intent(in) :: ovrlp, inv_ovrlp
        type(sparse_matrix),intent(in) :: smat
        real(kind=8),intent(out) :: max_error, mean_error
      
        real(kind=8), allocatable, dimension(:,:) :: tmp, tmp2
        real(kind=8), allocatable, dimension(:,:) :: tmpp, tmp2p
        integer :: ierr, i,j
      
        call f_routine(id='check_accur_overlap_minus_one')
      
        tmpp=f_malloc((/norb,norbp/),id='tmpp')
        if (power==1) then
           if (norbp>0) then
              call dgemm('n', 'n', norb, norbp, norb, 1.d0, inv_ovrlp(1,1), &
                   norb, ovrlp(1,isorb+1), norb, 0.d0, tmpp(1,1), norb)
           end if
           call deviation_from_unity_parallel(iproc, nproc, norb, norbp, isorb, tmpp, smat, max_error, mean_error)
        else if (power==2) then
           if (norbp>0) then
               call dgemm('n', 'n', norb, norbp, norb, 1.d0, inv_ovrlp(1,1), &
                   norb, inv_ovrlp(1,isorb+1), norb, 0.d0, tmpp(1,1), norb)
               call max_matrix_diff_parallel(iproc, norb, norbp, isorb, tmpp, ovrlp(1,isorb+1), smat, max_error, mean_error)
               max_error=0.5d0*max_error
               mean_error=0.5d0*mean_error
           else
               !still have to call routine so that allreduce is called by all mpi, but avoid explicit out of bounds reference for ovrlp
               call max_matrix_diff_parallel(iproc, norb, norbp, isorb, tmpp, ovrlp, smat, max_error, mean_error)
           end if
        else if (power==-2) then
           if (norbp>0) then
              call dgemm('n', 'n', norb, norbp, norb, 1.d0, inv_ovrlp(1,1), &
                   norb, inv_ovrlp(1,isorb+1), norb, 0.d0, tmpp(1,1), norb)
           end if
           tmp2p=f_malloc((/norb,norbp/),id='tmp2p')
           if (norbp>0) then
              call dgemm('n', 'n', norb, norbp, norb, 1.d0, ovrlp(1,1), &
                   norb, tmpp(1,1), norb, 0.d0, tmp2p(1,1), norb)
           end if
           call deviation_from_unity_parallel(iproc, nproc, norb, norbp, isorb, tmp2p, smat, max_error, mean_error)
           max_error=0.5d0*max_error
           mean_error=0.5d0*mean_error
           call f_free(tmp2p)
        else
           stop 'Error in check_accur_overlap_minus_one'
        end if
        call f_free(tmpp)
      
        call f_release_routine()
      
      end subroutine check_accur_overlap_minus_one


      subroutine deviation_from_unity_parallel_new(iproc, nproc, ovrlp, smat, max_deviation, mean_deviation)
        use module_base
        use sparsematrix_base, only: sparse_matrix
        use sparsematrix_init, only: matrixindex_in_compressed
        implicit none
      
        ! Calling arguments
        integer,intent(in):: iproc, nproc
        type(sparse_matrix),intent(in) :: smat
        real(8),dimension(smat%smmm%nvctrp),intent(in):: ovrlp
        real(8),intent(out):: max_deviation, mean_deviation
      
        ! Local variables
        integer:: iorb, iiorb, jorb, ierr, ind, iline, icolumn, i, ii
        real(8):: error, num
        real(kind=8),dimension(2) :: reducearr
      
        call f_routine(id='deviation_from_unity_parallel')
      
        call timing(iproc,'dev_from_unity','ON') 
        max_deviation=0.d0
        mean_deviation=0.d0
        num=0.d0
        !!do iorb=1,norbp
        !!   iiorb=iorb+isorb
        !!   !$omp parallel default(private) shared(norb, iiorb, ovrlp, iorb, max_deviation, mean_deviation, num, smat)
        !!   !$omp do reduction(max:max_deviation) reduction(+:mean_deviation,num)
        !!   do jorb=1,norb
        !!      ind=matrixindex_in_compressed(smat,jorb,iiorb)
        !!      if (ind>0) then
        !!          ! This entry is within the sparsity pattern, i.e. it matters for the error.
        !!          if(iiorb==jorb) then
        !!             error=(ovrlp(jorb,iorb)-1.d0)**2
        !!          else
        !!             error=ovrlp(jorb,iorb)**2
        !!          end if
        !!          max_deviation=max(error,max_deviation)
        !!          mean_deviation=mean_deviation+error
        !!          num=num+1.d0
        !!      end if
        !!   end do
        !!   !$omp end do
        !!   !$omp end parallel
        !!end do
      
      
        do i=1,smat%smmm%nvctrp
            ii = smat%smmm%isvctr + i
            !!call get_line_and_column(ii, smat%smmm%nseg, smat%smmm%keyv, smat%smmm%keyg, iline, icolumn)
            iline = smat%smmm%line_and_column(1,i)
            icolumn = smat%smmm%line_and_column(2,i)
            ind=matrixindex_in_compressed(smat,icolumn,iline)
            if (ind>0) then
                ! This entry is within the sparsity pattern, i.e. it matters for the error.
                if(icolumn==iline) then
                   error=(ovrlp(i)-1.d0)**2
                else
                   error=ovrlp(i)**2
                end if
                !write(*,'(a,3i8,2es16.7)') 'i, iline, icolumn, ovrlp(i), error', i, iline, icolumn, ovrlp(i), error
                max_deviation=max(error,max_deviation)
                mean_deviation=mean_deviation+error
                num=num+1.d0
            end if
        end do
      
      
      
        if (nproc>1) then
            reducearr(1)=mean_deviation
            reducearr(2)=num
            call mpiallred(max_deviation, 1, mpi_max, comm=bigdft_mpi%mpi_comm)
            call mpiallred(reducearr(1), 2, mpi_sum, comm=bigdft_mpi%mpi_comm)
            mean_deviation=reducearr(1)
            num=reducearr(2)
        end if
      
        mean_deviation=mean_deviation/num
        mean_deviation=sqrt(mean_deviation)
        max_deviation=sqrt(max_deviation)
      
        call timing(iproc,'dev_from_unity','OF') 
      
        call f_release_routine()
      
      end subroutine deviation_from_unity_parallel_new


      subroutine first_order_taylor_dense(norb,isorb,norbp,power,ovrlpp,inv_ovrlpp)
      use module_base
        implicit none
        integer,intent(in) :: norb, isorb, norbp, power
        real(kind=8),dimension(norb,norbp),intent(in) :: ovrlpp
        real(kind=8),dimension(norb,norbp),intent(out) :: inv_ovrlpp
      
        integer :: iorb, jorb, iiorb
      
        if (power==1) then
           !$omp parallel do default(private) shared(inv_ovrlpp,ovrlpp,norb,isorb,norbp)
           do iorb=1,norbp
              iiorb=isorb+iorb
              do jorb=1,norb
                 if(iiorb==jorb) then
                    inv_ovrlpp(jorb,iorb) = 2.d0 - ovrlpp(jorb,iorb)
                 else
                    inv_ovrlpp(jorb,iorb) = -ovrlpp(jorb,iorb)
                 end if
              end do
           end do
           !$omp end parallel do
        else if (power==2) then
           !$omp parallel do default(private) shared(inv_ovrlpp,ovrlpp,norb,isorb,norbp)
           do iorb=1,norbp
              iiorb=isorb+iorb
              do jorb=1,norb
                 if(iiorb==jorb) then
                    inv_ovrlpp(jorb,iorb) = 0.5d0 + 0.5d0*ovrlpp(jorb,iorb)
                 else
                    inv_ovrlpp(jorb,iorb) = 0.5d0*ovrlpp(jorb,iorb)
                 end if
              end do
           end do
           !$omp end parallel do
        else if (power==-2) then
           !$omp parallel do default(private) shared(inv_ovrlpp,ovrlpp,norb,isorb,norbp)
           do iorb=1,norbp
              iiorb=isorb+iorb
              do jorb=1,norb
                 if(iiorb==jorb) then
                    inv_ovrlpp(jorb,iorb) = 1.5d0 - 0.5d0*ovrlpp(jorb,iorb)
                 else
                    inv_ovrlpp(jorb,iorb) = -0.5d0*ovrlpp(jorb,iorb)
                 end if
              end do
           end do
           !$omp end parallel do
        else
           stop 'Error in first_order_taylor_dense'
        end if
      
      end subroutine first_order_taylor_dense


      subroutine matrix_minus_identity_sparse(norb, smat, ovrlp_compr, ovrlpminone_compr)
        use module_base
        use sparsematrix_base, only: sparse_matrix
        implicit none

        ! Calling arguments
        integer,intent(in) :: norb
        type(sparse_matrix),intent(in) :: smat
        real(kind=8),dimension(smat%nvctr),intent(in) :: ovrlp_compr
        real(kind=8),dimension(smat%nvctr),intent(out) :: ovrlpminone_compr

        ! Local variables
        integer :: ii, iseg, jorb, iiorb, jjorb

        !$omp parallel do default(private) &
        !$omp shared(smat, norb, ovrlpminone_compr, ovrlp_compr)
        do iseg=1,smat%nseg
            ii=smat%keyv(iseg)-1
             ! A segment is always on one line, therefore no double loop
            do jorb=smat%keyg(1,1,iseg),smat%keyg(2,1,iseg)
                iiorb = smat%keyg(1,2,iseg)
                jjorb = jorb
                ii=ii+1
                if (iiorb==jjorb) then
                    ovrlpminone_compr(ii)=ovrlp_compr(ii)-1.d0
                else
                    ovrlpminone_compr(ii)=ovrlp_compr(ii)
                end if
            end do
        end do
        !$omp end parallel do

      end subroutine matrix_minus_identity_sparse


      pure function newfactor(power,order,factor)
        implicit none
        integer, intent(in) :: power, order
        real(kind=8), intent(in) :: factor
        real(kind=8) :: newfactor
      
        if (power==1) then
           newfactor=-factor
        else if (power==2) then
           newfactor=factor*(1.5d0-real(order,kind=8))/real(order,kind=8)
        else if (power==-2) then
           newfactor=factor*(0.5d0-real(order,kind=8))/real(order,kind=8)
        end if
      
      end function newfactor


      subroutine matrix_minus_identity_dense(norb,isorb,norbp,matinp,matoutp)
        implicit none
        integer,intent(in) :: norb, isorb, norbp
        real(kind=8),dimension(norb,norbp),intent(in) :: matinp
        real(kind=8),dimension(norb,norbp),intent(out) :: matoutp
      
        integer :: iorb, jorb, iiorb
      
        !$omp parallel do default(private) shared(matinp,matoutp,norb,isorb,norbp)
        do iorb=1,norbp
           iiorb=isorb+iorb
           do jorb=1,norb
              if(iiorb==jorb) then
                 matoutp(jorb,iorb) = matinp(jorb,iorb) - 1.0d0
              else
                 matoutp(jorb,iorb) = matinp(jorb,iorb)
              end if
           end do
        end do
        !$omp end parallel do
      
      end subroutine matrix_minus_identity_dense


      subroutine first_order_taylor_sparse_new(power, smat, ovrlpp, inv_ovrlpp)
      use module_base
      use sparsematrix_base, only: sparse_matrix
        implicit none
        !!integer,intent(in) :: norb, isorb, norbp, power
        !!real(kind=8),dimension(norb,norbp),intent(in) :: ovrlpp
        !!real(kind=8),dimension(norb,norbp),intent(out) :: inv_ovrlpp
        integer,intent(in) :: power
        type(sparse_matrix),intent(in) :: smat
        real(kind=8),dimension(smat%smmm%nvctrp),intent(in) :: ovrlpp
        real(kind=8),dimension(smat%smmm%nvctrp),intent(out) :: inv_ovrlpp
      
        integer :: i, ii, iline, icolumn
      
        if (power==1) then
           !!!$omp parallel do default(private) shared(inv_ovrlpp,ovrlpp,norb,isorb,norbp)
           !!do iorb=1,norbp
           !!   iiorb=isorb+iorb
           !!   do jorb=1,norb
           !!      if(iiorb==jorb) then
           !!         inv_ovrlpp(jorb,iorb) = 2.d0 - ovrlpp(jorb,iorb)
           !!      else
           !!         inv_ovrlpp(jorb,iorb) = -ovrlpp(jorb,iorb)
           !!      end if
           !!   end do
           !!end do
           !!!$omp end parallel do
      
           do i=1,smat%smmm%nvctrp
               ii = smat%smmm%isvctr + i
               !!call get_line_and_column(ii, smat%smmm%nseg, smat%smmm%keyv, smat%smmm%keyg, iline, icolumn)
               iline = smat%smmm%line_and_column(1,i)
               icolumn = smat%smmm%line_and_column(2,i)
               if (iline==icolumn) then
                   inv_ovrlpp(i) = 2.d0 - ovrlpp(i)
               else
                   inv_ovrlpp(i) = -ovrlpp(i)
               end if
           end do
        else if (power==2) then
           !!!$omp parallel do default(private) shared(inv_ovrlpp,ovrlpp,norb,isorb,norbp)
           !!do iorb=1,norbp
           !!   iiorb=isorb+iorb
           !!   do jorb=1,norb
           !!      if(iiorb==jorb) then
           !!         inv_ovrlpp(jorb,iorb) = 0.5d0 + 0.5d0*ovrlpp(jorb,iorb)
           !!      else
           !!         inv_ovrlpp(jorb,iorb) = 0.5d0*ovrlpp(jorb,iorb)
           !!      end if
           !!   end do
           !!end do
           !!!$omp end parallel do
      
           do i=1,smat%smmm%nvctrp
               ii = smat%smmm%isvctr + i
               !!call get_line_and_column(ii, smat%smmm%nseg, smat%smmm%keyv, smat%smmm%keyg, iline, icolumn)
               iline = smat%smmm%line_and_column(1,i)
               icolumn = smat%smmm%line_and_column(2,i)
               if (iline==icolumn) then
                   inv_ovrlpp(i) = 0.5d0 + 0.5d0*ovrlpp(i)
               else
                   inv_ovrlpp(i) = 0.5d0*ovrlpp(i)
               end if
           end do
        else if (power==-2) then
           !!!$omp parallel do default(private) shared(inv_ovrlpp,ovrlpp,norb,isorb,norbp)
           !!do iorb=1,norbp
           !!   iiorb=isorb+iorb
           !!   do jorb=1,norb
           !!      if(iiorb==jorb) then
           !!         inv_ovrlpp(jorb,iorb) = 1.5d0 - 0.5d0*ovrlpp(jorb,iorb)
           !!      else
           !!         inv_ovrlpp(jorb,iorb) = -0.5d0*ovrlpp(jorb,iorb)
           !!      end if
           !!   end do
           !!end do
           !!!$omp end parallel do
      
           do i=1,smat%smmm%nvctrp
               ii = smat%smmm%isvctr + i
               !!call get_line_and_column(ii, smat%smmm%nseg, smat%smmm%keyv, smat%smmm%keyg, iline, icolumn)
               iline = smat%smmm%line_and_column(1,i)
               icolumn = smat%smmm%line_and_column(2,i)
               if (iline==icolumn) then
                   inv_ovrlpp(i) = 1.5d0 - 0.5d0*ovrlpp(i)
               else
                   inv_ovrlpp(i) = -0.5d0*ovrlpp(i)
               end if
           end do
        else
           stop 'Error in first_order_taylor_dense'
        end if
      
      end subroutine first_order_taylor_sparse_new


      subroutine overlap_minus_one_exact_serial(norb,inv_ovrlp)
        use module_base
        implicit none
        integer,intent(in) :: norb
        real(kind=8),dimension(norb,norb),intent(inout) :: inv_ovrlp
      
        integer :: info, iorb, jorb
      
        call f_routine(id='overlap_minus_one_exact_serial')
      
        call dpotrf('u', norb, inv_ovrlp(1,1), norb, info)
        if(info/=0) then
           write(*,'(1x,a,i0)') 'ERROR in dpotrf, info=',info
           stop
        end if
        call dpotri('u', norb, inv_ovrlp(1,1), norb, info)
        if(info/=0) then
           write(*,'(1x,a,i0)') 'ERROR in dpotri, info=',info
           stop
        end if
      
        ! fill lower half
        !$omp parallel do default(private) shared(inv_ovrlp,norb)
        do iorb=1,norb
           do jorb=1,iorb-1
              inv_ovrlp(iorb,jorb)=inv_ovrlp(jorb,iorb)
           end do
        end do
        !$omp end parallel do 
      
        call f_release_routine()
      
      end subroutine overlap_minus_one_exact_serial


      subroutine overlap_plus_minus_one_half_exact(nproc,norb,blocksize,plusminus,inv_ovrlp_half,smat)
        use module_base
        use sparsematrix_base, only: sparse_matrix
        use parallel_linalg, only: dgemm_parallel, dsyev_parallel
        implicit none
        integer,intent(in) :: nproc,norb,blocksize
        real(kind=8),dimension(norb,norb) :: inv_ovrlp_half
        logical, intent(in) :: plusminus
        type(sparse_matrix),intent(in) :: smat
      
        integer :: info, iorb, jorb, ierr, iiorb, isorb, norbp, lwork, jjorb,ninetynine
        real(kind=8),dimension(:),allocatable :: eval, work
        real(kind=8),dimension(:,:),allocatable :: tempArr, orig_ovrlp
        real(kind=8),dimension(:,:),pointer :: inv_ovrlp_halfp
        real(kind=8),dimension(:,:), allocatable :: vr,vl ! for non-symmetric LAPACK
        real(kind=8),dimension(:),allocatable:: eval1 ! for non-symmetric LAPACK
        real(dp) :: temp, max_error, mean_error
        real(dp), allocatable, dimension(:) :: temp_vec
        logical, parameter :: symmetric=.true.
        logical, parameter :: check_lapack=.true.
        integer :: korb
      
      
        call f_routine(id='overlap_plus_minus_one_half_exact')
      
        !!!if (nproc>1) then
        !!    if (.not.present(smat)) then 
        !!        call f_err_throw('overlap_plus_minus_one_half_exact: for nproc>1, smat must be present!', &
        !!             err_name='BIGDFT_RUNTIME_ERROR')
        !!    end if
        !!!end if
                 
      
        eval=f_malloc(norb,id='eval')
        if(blocksize>0) then
           call dsyev_parallel(bigdft_mpi%iproc, bigdft_mpi%nproc, min(blocksize,norb), bigdft_mpi%mpi_comm, 'v', 'l', norb, &
                inv_ovrlp_half(1,1), norb, eval(1), info)
           if(info/=0) then
              write(*,'(a,i0)') 'ERROR in dsyev_parallel, info=', info
           end if
        else
           if (symmetric) then
              if (check_lapack) then
                 orig_ovrlp=f_malloc((/norb,norb/),id='orig_ovrlp')
                 call vcopy(norb*norb,inv_ovrlp_half(1,1),1,orig_ovrlp(1,1),1)
                 !the original overlap has to be symmetrized
                 do iorb=1,norb
                    do jorb=1,iorb-1
                       orig_ovrlp(jorb,iorb)=orig_ovrlp(iorb,jorb)
                    end do
                 end do
              end if
              work=f_malloc(1000,id='work')
              call dsyev('v', 'l', norb, inv_ovrlp_half(1,1), norb, eval, work, -1, info)
              lwork = int(work(1))
              call f_free(work)
              work=f_malloc(lwork,id='work')
                 !!do iorb=1,norb
                 !!   do jorb=1,norb
                 !!      write(2000+bigdft_mpi%iproc,'(a,3i8,es16.7)') 'iproc, iorb, jorb, val', bigdft_mpi%iproc, iorb, jorb, inv_ovrlp_half(jorb,iorb)
                 !!   end do
                 !!end do
              !!do jorb=1,norb
              !!    do korb=1,norb
              !!        write(910,'(a,2i8,es14.5)') 'jorb, korb, inv_ovrlp_half(korb,jorb)', jorb, korb, inv_ovrlp_half(korb,jorb)
              !!    end do
              !!end do
              call dsyev('v', 'l', norb, inv_ovrlp_half(1,1), norb, eval, work, lwork, info)
              !!do jorb=1,norb
              !!    do korb=1,norb
              !!        write(920,'(a,2i8,es14.5)') 'jorb, korb, inv_ovrlp_half(korb,jorb)', jorb, korb, inv_ovrlp_half(korb,jorb)
              !!    end do
              !!end do
              if (check_lapack) then
                 tempArr=f_malloc((/norb,norb/), id='tempArr')
                 do iorb=1,norb
                    do jorb=1,norb
                       tempArr(jorb,iorb)=inv_ovrlp_half(jorb,iorb)*eval(iorb)
                    end do
                 end do
                 inv_ovrlp_halfp=f_malloc_ptr((/norb,norb/), id='inv_ovrlp_halfp')
                 call dgemm('n', 't', norb, norb, norb, 1.d0, inv_ovrlp_half, &
                      norb, tempArr, norb, 0.d0, inv_ovrlp_halfp, norb)
                 call f_free(tempArr)
                 call max_matrix_diff(bigdft_mpi%iproc, norb, inv_ovrlp_halfp, orig_ovrlp, smat, max_error, mean_error)
                 if (abs(max_error)>1.0d-8) then
                    if (bigdft_mpi%iproc==0) then
                       ninetynine=f_get_free_unit(99)
                       open(ninetynine,file='dsyev_input.txt')
                       do iorb=1,norb
                          do jorb=1,norb
                             write(ninetynine,*) iorb,jorb,orig_ovrlp(iorb,jorb)
                          end do
                       end do
                       call f_close(ninetynine)
                       open(ninetynine,file='dsyev_output.txt')
                       do iorb=1,norb
                          do jorb=1,norb
                             write(ninetynine,*) iorb,jorb,inv_ovrlp_halfp(iorb,jorb)
                          end do
                       end do
                       call f_close(ninetynine)
                    end if
                    call f_err_throw('LAPACK error for dsyev in overlap_plus_minus_one_half_exact (maxerr='//&
                         trim(yaml_toa(max_error))//'), erroneous matrices dumped in files "dsyev_(in/out)put.txt"',&
                         err_name='BIGDFT_RUNTIME_ERROR')
                 end if
                 call f_free_ptr(inv_ovrlp_halfp)
              end if
           else
              if (check_lapack) then
                 orig_ovrlp=f_malloc((/norb,norb/),id='orig_ovrlp')
                 call vcopy(norb*norb,inv_ovrlp_half(1,1),1,orig_ovrlp(1,1),1)
              end if
              work=f_malloc(1000,id='work')
              eval1=f_malloc(norb,id='eval1')
              vl=f_malloc((/norb,norb/),id='vl')
              vr=f_malloc((/norb,norb/),id='vr')
              call dgeev( 'v','v', norb, inv_ovrlp_half(1,1), norb, eval, eval1, VL, norb, VR,&
                   norb, WORK, -1, info )
              lwork = nint(work(1))
              call f_free(work)
              work=f_malloc(lwork,id='work')
              call DGEEV( 'v','v', norb, inv_ovrlp_half(1,1), norb, eval, eval1, VL, norb, VR,&
                   norb, WORK, LWORK, info )
              call vcopy(norb*norb,vl(1,1),1,inv_ovrlp_half(1,1),1)
              call f_free(eval1)
              call f_free(vr)
              call f_free(vl)
              temp_vec=f_malloc(norb,id='temp_vec')
              do iorb=1,norb
                 do jorb=iorb+1,norb
                    if (eval(jorb) < eval(iorb)) then
                       temp = eval(iorb)
                       temp_vec = inv_ovrlp_half(:,iorb)
                       eval(iorb) = eval(jorb)
                       eval(jorb) = temp
                       inv_ovrlp_half(:,iorb) = inv_ovrlp_half(:,jorb)
                       inv_ovrlp_half(:,jorb) = temp_vec
                    end if
                 end do
              end do
              call f_free(temp_vec)
              if (check_lapack) then
                 tempArr=f_malloc((/norb,norb/), id='tempArr')
                 do iorb=1,norb
                    do jorb=1,norb
                       tempArr(jorb,iorb)=inv_ovrlp_half(jorb,iorb)*eval(iorb)
                    end do
                 end do
                 inv_ovrlp_halfp=f_malloc_ptr((/norb,norb/), id='inv_ovrlp_halfp')
                 call dgemm('n', 't', norb, norb, norb, 1.d0, inv_ovrlp_half, &
                      norb, tempArr, norb, 0.d0, inv_ovrlp_halfp, norb)
                 call f_free(tempArr)
                 call max_matrix_diff(bigdft_mpi%iproc, norb, inv_ovrlp_halfp, orig_ovrlp, smat, max_error, mean_error)
                 if (bigdft_mpi%iproc==0.and.abs(max_error)>1.0d-8) then
                    print*,'LAPACK error for dgeev in overlap_plus_minus_one_half_exact',max_error
                    open(99,file='dgeev_input.txt')
                    do iorb=1,norb
                       do jorb=1,norb
                         write(99,*) iorb,jorb,orig_ovrlp(iorb,jorb)
                       end do
                    end do
                    close(99)
                    open(99,file='dgeev_output.txt')
                    do iorb=1,norb
                       do jorb=1,norb
                         write(99,*) iorb,jorb,inv_ovrlp_halfp(iorb,jorb)
                       end do
                    end do
                    close(99)
                    call mpi_finalize(bigdft_mpi%mpi_comm)
                    stop
                 end if
                 call f_free_ptr(inv_ovrlp_halfp)
              end if
           end if
      
           if(info/=0) then
              write(*,'(a,i0,2x,i0)') 'ERROR in dsyev (overlap_plus_minus_one_half_exact), info, norb=', info, norb
              open(99,file='dsyev_overlap.txt')
              do iorb=1,norb
                 do jorb=1,norb
                    write(99,*) iorb,jorb,orig_ovrlp(iorb,jorb)
                 end do
              end do
              close(99)
              stop
           end if
           call f_free(orig_ovrlp)
           call f_free(work)
        end if
      
        ! Calculate S^{-1/2}. 
        ! First calculate ovrlp*diag(1/sqrt(eval)) (ovrlp is the diagonalized overlap
        ! matrix and diag(1/sqrt(evall)) the diagonal matrix consisting of the inverse square roots of the eigenvalues...
        !write(*,*) 'iproc, present(orbs), norb, orbs%norb', bigdft_mpi%iproc, present(orbs), norb, orbs%norb
        !if (present(orbs).and.bigdft_mpi%nproc>1.and.blocksize<0) then
           !if (norb/=orbs%norb) stop 'Error with orbs%norb in overlap_plus_minus_one_half_exact'
           if (nproc>1) then
               norbp=smat%nfvctrp
               isorb=smat%isfvctr
           else
               norbp=norb
               isorb=0
           end if
        !else
        !   norbp=norb
        !   isorb=0
        !end if
        tempArr=f_malloc((/norbp,norb/), id='tempArr')
        !$omp parallel do default(private) shared(tempArr,inv_ovrlp_half,eval,plusminus,norb,norbp,isorb)
        do iorb=1,norb
           do jorb=1,norbp
              jjorb=jorb+isorb
              if (plusminus) then
                 tempArr(jorb,iorb)=inv_ovrlp_half(jjorb,iorb)*sqrt(abs(eval(iorb)))
              else
                 tempArr(jorb,iorb)=inv_ovrlp_half(jjorb,iorb)*1.d0/sqrt(abs(eval(iorb)))
              end if
           end do
        end do
        !$omp end parallel do
      
 
        call f_free(eval)
        inv_ovrlp_halfp=f_malloc_ptr((/norb,norbp/), id='inv_ovrlp_halfp')
        ! ...and now apply the diagonalized overlap matrix to the matrix constructed above.
        ! This will give S^{+/-1/2}.
        !!if(blocksize<0) then
           if (norbp>0) call dgemm('n', 't', norb, norbp, norb, 1.d0, inv_ovrlp_half, &
                norb, tempArr, norbp, 0.d0, inv_ovrlp_halfp, norb)
           !if (present(orbs).and.bigdft_mpi%nproc>1) then
           if (nproc>1) then
              call mpi_allgatherv(inv_ovrlp_halfp, norb*norbp, mpi_double_precision, inv_ovrlp_half, &
                         norb*smat%nfvctr_par(:), norb*smat%isfvctr_par, mpi_double_precision, bigdft_mpi%mpi_comm, ierr)
           else
              call vcopy(norb*norbp,inv_ovrlp_halfp(1,1),1,inv_ovrlp_half(1,1),1)
           end if
        !!else
        !!   call dgemm_parallel(bigdft_mpi%iproc, bigdft_mpi%nproc, blocksize, bigdft_mpi%mpi_comm, 'n', 't', norb, norb, norb, &
        !!        1.d0, inv_ovrlp_half, norb, tempArr, norb, 0.d0, inv_ovrlp_halfp, norb)
        !!   call vcopy(norb*norbp,inv_ovrlp_halfp(1,1),1,inv_ovrlp_half(1,1),1)
        !!end if
  
        call f_free_ptr(inv_ovrlp_halfp)
        call f_free(tempArr)
      
        call f_release_routine()
      
      end subroutine overlap_plus_minus_one_half_exact


      subroutine deviation_from_unity_parallel(iproc, nproc, norb, norbp, isorb, ovrlp, smat, max_deviation, mean_deviation)
        use module_base
        use sparsematrix_base, only: sparse_matrix
        use sparsematrix_init, only: matrixindex_in_compressed
        implicit none
      
        ! Calling arguments
        integer,intent(in):: iproc, nproc, norb, norbp, isorb
        real(8),dimension(norb,norbp),intent(in):: ovrlp
        type(sparse_matrix),intent(in) :: smat
        real(8),intent(out):: max_deviation, mean_deviation
      
        ! Local variables
        integer:: iorb, iiorb, jorb, ierr, ind
        real(8):: error, num
        real(kind=8),dimension(2) :: reducearr
      
        call f_routine(id='deviation_from_unity_parallel')
      
        call timing(iproc,'dev_from_unity','ON') 
        max_deviation=0.d0
        mean_deviation=0.d0
        num=0.d0
        do iorb=1,norbp
           iiorb=iorb+isorb
           !$omp parallel default(private) shared(norb, iiorb, ovrlp, iorb, max_deviation, mean_deviation, num, smat)
           !$omp do reduction(max:max_deviation) reduction(+:mean_deviation,num)
           do jorb=1,norb
              ind=matrixindex_in_compressed(smat,jorb,iiorb)
              if (ind>0) then
                  ! This entry is within the sparsity pattern, i.e. it matters for the error.
                  if(iiorb==jorb) then
                     error=(ovrlp(jorb,iorb)-1.d0)**2
                  else
                     error=ovrlp(jorb,iorb)**2
                  end if
                  max_deviation=max(error,max_deviation)
                  mean_deviation=mean_deviation+error
                  num=num+1.d0
              end if
           end do
           !$omp end do
           !$omp end parallel
        end do
        if (nproc>1) then
            reducearr(1)=mean_deviation
            reducearr(2)=num
            call mpiallred(max_deviation, 1, mpi_max, comm=bigdft_mpi%mpi_comm)
            call mpiallred(reducearr(1), 2, mpi_sum, comm=bigdft_mpi%mpi_comm)
            mean_deviation=reducearr(1)
            num=reducearr(2)
        end if
      
        mean_deviation=mean_deviation/num
        mean_deviation=sqrt(mean_deviation)
        max_deviation=sqrt(max_deviation)
      
        call timing(iproc,'dev_from_unity','OF') 
      
        call f_release_routine()
      
      end subroutine deviation_from_unity_parallel


      subroutine max_matrix_diff(iproc, norb, mat1, mat2, smat, max_deviation, mean_deviation)
        use module_base
        use sparsematrix_base, only: sparse_matrix
        use sparsematrix_init, only: matrixindex_in_compressed
        implicit none
      
        ! Calling arguments
        integer,intent(in):: iproc, norb
        real(8),dimension(norb,norb),intent(in):: mat1, mat2
        type(sparse_matrix),intent(in) :: smat
        real(8),intent(out):: max_deviation, mean_deviation
      
        ! Local variables
        integer:: iorb, jorb, ind
        real(8):: error, num
      
        call timing(iproc,'dev_from_unity','ON') 
        max_deviation=0.d0
        mean_deviation=0.d0
        num=0.d0
        do iorb=1,norb
           do jorb=1,norb
              ind=matrixindex_in_compressed(smat,jorb,iorb)
              if (ind>0) then
                  error=(mat1(jorb,iorb)-mat2(jorb,iorb))**2
                  max_deviation=max(error,max_deviation)
                  mean_deviation=mean_deviation+error
                  num=num+1.d0
              end if
           end do
        end do
        mean_deviation=mean_deviation/num
        mean_deviation=sqrt(mean_deviation)
        max_deviation=sqrt(max_deviation)
        call timing(iproc,'dev_from_unity','OF') 
      
      end subroutine max_matrix_diff


      subroutine max_matrix_diff_parallel(iproc, norb, norbp, isorb, mat1, mat2, &
                 smat, max_deviation, mean_deviation)
        use module_base
        use sparsematrix_base, only: sparse_matrix
        use sparsematrix_init, only: matrixindex_in_compressed
        implicit none
      
        ! Calling arguments
        integer,intent(in):: iproc, norb, norbp, isorb
        real(8),dimension(norb,norbp),intent(in):: mat1, mat2
        type(sparse_matrix),intent(in) :: smat
        real(8),intent(out):: max_deviation, mean_deviation
      
        ! Local variables
        integer:: iorb, iiorb, jorb, ierr, ind
        real(8):: error, num
        real(kind=8),dimension(2) :: reducearr
      
        call timing(iproc,'dev_from_unity','ON') 
        max_deviation=0.d0
        mean_deviation=0.d0
        num=0.d0
        do iorb=1,norbp
           iiorb=iorb+isorb
           !$omp parallel default(private) shared(iorb, iiorb, norb, mat1, mat2, max_deviation, mean_deviation, num, smat)
           !$omp do reduction(max:max_deviation) reduction(+:mean_deviation,num)
           do jorb=1,norb
              ind=matrixindex_in_compressed(smat,jorb,iiorb)
              if (ind>0) then
                  ! This entry is within the sparsity pattern, i.e. it matters for the error.
                  error=(mat1(jorb,iorb)-mat2(jorb,iorb))**2
                  max_deviation=max(error,max_deviation)
                  mean_deviation=mean_deviation+error
                  num=num+1.d0
              end if
           end do
           !$omp end do
           !$omp end parallel
        end do
      
        if (bigdft_mpi%nproc>1) then
            reducearr(1)=mean_deviation
            reducearr(2)=num
            call mpiallred(max_deviation, 1, mpi_max, comm=bigdft_mpi%mpi_comm)
            call mpiallred(reducearr, mpi_sum, comm=bigdft_mpi%mpi_comm)
            mean_deviation=reducearr(1)
            num=reducearr(2)
        end if
      
        mean_deviation=mean_deviation/num
        mean_deviation=sqrt(mean_deviation)
        max_deviation=sqrt(max_deviation)
      
        call timing(iproc,'dev_from_unity','OF') 
      
      end subroutine max_matrix_diff_parallel


      subroutine max_matrix_diff_parallel_new(iproc, norb, norbp, isorb, mat1, mat2, &
                 smat, max_deviation, mean_deviation)
        use module_base
        use sparsematrix_base, only: sparse_matrix
        use sparsematrix_init, only: matrixindex_in_compressed
        implicit none
      
        ! Calling arguments
        integer,intent(in):: iproc, norb, norbp, isorb
        type(sparse_matrix),intent(in) :: smat
        real(8),dimension(smat%smmm%nvctrp),intent(in):: mat1, mat2
        real(8),intent(out):: max_deviation, mean_deviation
      
        ! Local variables
        integer:: iorb, iiorb, jorb, ierr, ind, iline, icolumn, i, ii
        real(8):: error, num
        real(kind=8),dimension(2) :: reducearr
      
        call timing(iproc,'dev_from_unity','ON') 
        max_deviation=0.d0
        mean_deviation=0.d0
        num=0.d0
        !!do iorb=1,norbp
        !!   iiorb=iorb+isorb
        !!   !$omp parallel default(private) shared(iorb, iiorb, norb, mat1, mat2, max_deviation, mean_deviation, num, smat)
        !!   !$omp do reduction(max:max_deviation) reduction(+:mean_deviation,num)
        !!   do jorb=1,norb
        !!      ind=matrixindex_in_compressed(smat,jorb,iiorb)
        !!      if (ind>0) then
        !!          ! This entry is within the sparsity pattern, i.e. it matters for the error.
        !!          error=(mat1(jorb,iorb)-mat2(jorb,iorb))**2
        !!          max_deviation=max(error,max_deviation)
        !!          mean_deviation=mean_deviation+error
        !!          num=num+1.d0
        !!      end if
        !!   end do
        !!   !$omp end do
        !!   !$omp end parallel
        !!end do
      
        do i=1,smat%smmm%nvctrp
            ii = smat%smmm%isvctr + i
            !!call get_line_and_column(ii, smat%smmm%nseg, smat%smmm%keyv, smat%smmm%keyg, iline, icolumn)
            iline = smat%smmm%line_and_column(1,i)
            icolumn = smat%smmm%line_and_column(2,i)
            ind=matrixindex_in_compressed(smat,icolumn,iline)
            if (ind>0) then
                ! This entry is within the sparsity pattern, i.e. it matters for the error.
                error=(mat1(i)-mat2(i))**2
                max_deviation=max(error,max_deviation)
                mean_deviation=mean_deviation+error
                num=num+1.d0
            end if
        end do
      
        if (bigdft_mpi%nproc>1) then
            reducearr(1)=mean_deviation
            reducearr(2)=num
            call mpiallred(max_deviation, 1, mpi_max, comm=bigdft_mpi%mpi_comm)
            call mpiallred(reducearr, mpi_sum, comm=bigdft_mpi%mpi_comm)
            mean_deviation=reducearr(1)
            num=reducearr(2)
        end if
      
        mean_deviation=mean_deviation/num
        mean_deviation=sqrt(mean_deviation)
        max_deviation=sqrt(max_deviation)
      
        call timing(iproc,'dev_from_unity','OF') 
      
      end subroutine max_matrix_diff_parallel_new


      subroutine overlap_power_minus_one_half_parallel(iproc, nproc, meth_overlap, ovrlp, ovrlp_mat, &
                 inv_ovrlp_half, inv_ovrlp_half_)
        use module_base
        use sparsematrix_base, only: sparse_matrix, matrices, matrices_null, &
                                     allocate_matrices, deallocate_matrices
        use sparsematrix_init, only: matrixindex_in_compressed
        use sparsematrix, only: synchronize_matrix_taskgroups
        implicit none
      
        ! Calling arguments
        integer,intent(in) :: iproc, nproc, meth_overlap
        type(sparse_matrix),intent(in) :: ovrlp
        type(matrices),intent(inout) :: ovrlp_mat
        type(sparse_matrix),intent(in) :: inv_ovrlp_half
        type(matrices),intent(inout) :: inv_ovrlp_half_
      
        ! Local variables
        integer(kind=8) :: ii, iend
        integer :: i, iorb, n, istat, iall, jorb, korb, jjorb, kkorb!, ilr
        integer :: iiorb, ierr, iseg, ind, ishift_ovrlp, ishift_inv_ovrlp, ispin
        real(kind=8) :: error
        real(kind=8),dimension(:,:),pointer :: ovrlp_tmp, ovrlp_tmp_inv_half
        logical,dimension(:),allocatable :: in_neighborhood
        character(len=*),parameter :: subname='overlap_power_minus_one_half_parallel'
        !type(matrices) :: inv_ovrlp_half_
        !!integer :: itaskgroups, iitaskgroup, imin, imax
      
        !!imin=ovrlp%nvctr
        !!imax=0
        !!do itaskgroups=1,ovrlp%ntaskgroupp
        !!    iitaskgroup = ovrlp%taskgroupid(itaskgroups)
        !!    imin = min(imin,ovrlp%taskgroup_startend(1,1,iitaskgroup))
        !!    imax = max(imax,ovrlp%taskgroup_startend(2,1,iitaskgroup))
        !!end do
      
      
        call f_routine('overlap_power_minus_one_half_parallel')
        call timing(iproc,'lovrlp^-1/2par','ON')
      
        in_neighborhood = f_malloc(ovrlp%nfvctr,id='in_neighborhood')
      
        !inv_ovrlp_half_ = matrices_null()
        !call allocate_matrices(inv_ovrlp_half, allocate_full=.false., matname='inv_ovrlp_half_', mat=inv_ovrlp_half_)
        call f_zero(inv_ovrlp_half%nvctrp_tg*inv_ovrlp_half%nspin, inv_ovrlp_half_%matrix_compr(1))
      
        !DEBUG
        !if (iproc==0) then
        !   jjorb=0
        !   do jorb=1,orbs%norb
        !      do korb=1,orbs%norb
        !         ind = ovrlp%matrixindex_in_compressed(korb, jorb)
        !         if (ind>0) then
        !            print*, korb,jorb,ovrlp%matrix_compr(ind)
        !         else
        !            print*, korb,jorb,0.0d0
        !         end if
        !         !write(1200+iproc,'(2i8,es20.10)') kkorb, jjorb, ovrlp_tmp(kkorb,jjorb)
        !      end do
        !   end do
        !end if
        !call mpi_barrier(bigdft_mpi%mpi_comm,istat)
      
        !!ii=0
        !!do ispin=1,ovrlp%nspin
        !!    do i=1,ovrlp%nvctr
        !!        ii=ii+1
        !!        write(950+iproc,*) 'ii, i, val', ii, i, ovrlp_mat%matrix_compr(ii)
        !!    end do
        !!end do
      
      
        spin_loop: do ispin=1,ovrlp%nspin
      
            ishift_ovrlp=(ispin-1)*ovrlp%nvctr
            ishift_inv_ovrlp=(ispin-1)*inv_ovrlp_half%nvctr
      
            do iorb=1,ovrlp%nfvctrp
               iiorb=ovrlp%isfvctr+iorb
               !ilr=orbs%inwhichlocreg(iiorb)
               ! We are at the start of a new atom
               ! Count all orbitals that are in the neighborhood
      
               iseg=ovrlp%istsegline(iiorb)
               iend=int(iiorb,kind=8)*int(ovrlp%nfvctr,kind=8)
               n=0
               in_neighborhood(:)=.false.
               do 
                  do i=ovrlp%keyg(1,1,iseg),ovrlp%keyg(2,1,iseg)
                     in_neighborhood(i)=.true.
                     !if (iproc==0) write(*,*) 'iiorb, iseg, i, n', iiorb, iseg, i, n
                     n=n+1
                  end do
                  iseg=iseg+1
                  if (iseg>ovrlp%nseg) exit
                  ii = int((ovrlp%keyg(1,2,iseg)-1),kind=8)*int(ovrlp%nfvctr,kind=8) + &
                       int(ovrlp%keyg(1,1,iseg),kind=8)
                  if (ii>iend) exit
               end do
               !if (iproc==0) write(*,*) 'iiorb, n', iiorb, n
      
               ovrlp_tmp = f_malloc0_ptr((/n,n/),id='ovrlp_tmp')
      
               jjorb=0
               do jorb=1,ovrlp%nfvctr
                  if (.not.in_neighborhood(jorb)) cycle
                  jjorb=jjorb+1
                  kkorb=0
                  do korb=1,ovrlp%nfvctr
                     if (.not.in_neighborhood(korb)) cycle
                     kkorb=kkorb+1
                     ind = matrixindex_in_compressed(ovrlp,korb,jorb)
                     if (ind>0) then
                        !!if (ind<imin) then
                        !!    write(*,*) 'ind,imin',ind,imin
                        !!    stop 'ind<imin'
                        !!end if
                        !!if (ind>imax) then
                        !!    write(*,*) 'ind,imax',ind,imax
                        !!    stop 'ind>imax'
                        !!end if
                        ind=ind+ishift_ovrlp
                        ovrlp_tmp(kkorb,jjorb)=ovrlp_mat%matrix_compr(ind-ovrlp%isvctrp_tg)
                     else
                        ovrlp_tmp(kkorb,jjorb)=0.d0
                     end if
                     !!write(1200+iproc,'(2i8,es20.10)') kkorb, jjorb, ovrlp_tmp(kkorb,jjorb)
                  end do
               end do
                    
               ovrlp_tmp_inv_half = f_malloc_ptr((/n,n/),id='ovrlp_tmp_inv_half')
               call vcopy(n*n, ovrlp_tmp(1,1), 1, ovrlp_tmp_inv_half(1,1), 1)
               !!do jorb=1,n
               !!    do korb=1,n
               !!        write(900,'(a,2i8,es14.5)') 'jorb, korb, ovrlp_tmp(korb,jorb)', jorb, korb, ovrlp_tmp(korb,jorb)
               !!    end do
               !!end do
      
               !if (iiorb==orbs%norb) then
               !print*,''
               !print*,'ovrlp_tmp',n,iiorb
               !do jorb=1,n
               !print*,jorb,ovrlp_tmp(:,jorb)
               !end do
               !end if
      
      
               ! Calculate S^-1/2 for the small overlap matrix
               !!call overlapPowerGeneral(iproc, nproc, meth_overlap, -2, -8, n, orbs, imode=2, check_accur=.true.,&
               !!     ovrlp=ovrlp_tmp, inv_ovrlp=ovrlp_tmp_inv_half, error=error)
               !!call overlapPowerGeneral(iproc, 1, meth_overlap, -2, -8, n, orbs, imode=2, &
               !!     ovrlp_smat=ovrlp, inv_ovrlp_smat=inv_ovrlp_half, &
               !!     ovrlp_mat=ovrlp_mat, inv_ovrlp_mat=inv_ovrlp_half_, check_accur=.true., &
               !!     ovrlp=ovrlp_tmp, inv_ovrlp=ovrlp_tmp_inv_half, error=error)
               call overlap_plus_minus_one_half_exact(1, n, -8, .false., ovrlp_tmp_inv_half,inv_ovrlp_half)
      
      
               !if (iiorb==orbs%norb) then
               !print*,''
               !print*,'inv_ovrlp_tmp',n,iiorb,error
               !do jorb=1,n
               !print*,jorb,ovrlp_tmp_inv_half(:,jorb)
               !end do
               !end if
      
               jjorb=0
               do jorb=1,ovrlp%nfvctr
                  if (.not.in_neighborhood(jorb)) cycle
                  jjorb=jjorb+1
                  kkorb=0
                  if (jorb==iiorb) then
                     do korb=1,ovrlp%nfvctr
                        if (.not.in_neighborhood(korb)) cycle
                        kkorb=kkorb+1
                        ind = matrixindex_in_compressed(inv_ovrlp_half,korb,jorb)
                        if (ind>0) then
                           ind=ind+ishift_inv_ovrlp
                           inv_ovrlp_half_%matrix_compr(ind-inv_ovrlp_half%isvctrp_tg)=ovrlp_tmp_inv_half(kkorb,jjorb)
                           !if (iproc==0) write(*,'(a,6i8,es16.8)') 'ind, inv_ovrlp_half%isvctrp_tg, jorb, korb, jjorb, kkorb, val', &
                           !                  ind, inv_ovrlp_half%isvctrp_tg, jorb, korb, jjorb, kkorb, ovrlp_tmp_inv_half(kkorb,jjorb)
                           !if (iiorb==orbs%norb) print*,'problem here?!',iiorb,kkorb,jjorb,korb,jorb,ind,ovrlp_tmp_inv_half(kkorb,jjorb)
                        end if
                        !write(1300+iproc,'(2i8,es20.10)') kkorb, jjorb, ovrlp_tmp(kkorb,jjorb)
                     end do
                     exit !no need to keep looping
                  end if
               end do
      
      
               call f_free_ptr(ovrlp_tmp_inv_half)
               call f_free_ptr(ovrlp_tmp)
      
            end do
      
            !!if (nproc>1)then
            !!    call mpiallred(inv_ovrlp_half_%matrix_compr(1), inv_ovrlp_half%nvctr*inv_ovrlp_half%nspin, mpi_sum, bigdft_mpi%mpi_comm)
            !!end if
            call synchronize_matrix_taskgroups(iproc, nproc, inv_ovrlp_half, inv_ovrlp_half_)
      
        end do spin_loop
      
        call f_free(in_neighborhood)
      
        !if (iproc==0) then
        !   jjorb=0
        !   do jorb=1,orbs%norb
        !      do korb=1,orbs%norb
        !         ind = inv_ovrlp_half%matrixindex_in_compressed(korb, jorb)
        !         if (ind>0) then
        !            print*, korb,jorb,inv_ovrlp_half%matrix_compr(ind)
        !         else
        !            print*, korb,jorb,0.0d0
        !         end if
        !         !write(1200+iproc,'(2i8,es20.10)') kkorb, jjorb, ovrlp_tmp(kkorb,jjorb)
        !      end do
        !   end do
        !end if
        !call mpi_finalize(ind)
        !stop
      
        !call deallocate_matrices(inv_ovrlp_half_)
      
      
        call f_release_routine
        call timing(iproc,'lovrlp^-1/2par','OF')
      
      end subroutine overlap_power_minus_one_half_parallel


    subroutine check_taylor_order(error, max_error, order_taylor)
      use module_base
      use yaml_output
      implicit none
    
      ! Calling arguments
      real(kind=8),intent(in) :: error, max_error
      integer,intent(inout) :: order_taylor
    
      ! Local variables
      character(len=12) :: act
      integer,parameter :: max_order_positive=50
      integer,parameter :: max_order_negative=-20
      logical :: is_ice
    
      if (order_taylor>=1000) then
          order_taylor=order_taylor-1000
          is_ice=.true.
      else
          is_ice=.false.
      end if
    
      if (order_taylor/=0) then
          ! only do this if approximations (Taylor or "negative thing") are actually used
          if (error<=1.d-1*max_error) then
              !! error is very small, so decrease the order of the polynomial
              !if (order_taylor>20) then
              !    ! always keep a minimum of 20
              !    act=' (decreased)'
              !    if (order_taylor>0) then
              !        order_taylor = floor(0.9d0*real(order_taylor,kind=8))
              !    else
              !        order_taylor = ceiling(0.9d0*real(order_taylor,kind=8))
              !    end if
              !end if
          else if (error>max_error) then
              ! error is too big, increase the order of the Taylor series by 10%
              act=' (increased)'
              if (order_taylor>0) then
                  order_taylor = ceiling(max(1.1d0*real(order_taylor,kind=8),real(order_taylor+5,kind=8)))
              else
                  order_taylor = floor(min(1.1d0*real(order_taylor,kind=8),real(order_taylor-5,kind=8)))
              end if
          else
              ! error is small enough, so do nothing
              act=' (unchanged)'
          end if
          !if (bigdft_mpi%iproc==0) call yaml_map('new Taylor order',trim(yaml_toa(order_taylor,fmt='(i0)'))//act)
      end if
    
      if (order_taylor>0) then
          if (order_taylor>max_order_positive) then
              order_taylor=max_order_positive
              if (bigdft_mpi%iproc==0) call yaml_warning('Taylor order reached maximum')
          end if
      else
          if (order_taylor<max_order_negative) then
              order_taylor=max_order_negative
              if (bigdft_mpi%iproc==0) call yaml_warning('Taylor order reached maximum')
          end if
      end if
    
      if (is_ice) then
          order_taylor=order_taylor+1000
      end if
    
    end subroutine check_taylor_order


    !< Calculate "local versions" of S^{-1/2}, i.e. take only the small subblocks of all support functions
    !! on one atom and calculate S^{-1/2} for this subblock. The remaining parts of teh matrix are empty.
    subroutine calculate_S_minus_one_half_onsite(iproc, nproc, norb, onwhichatom, smats, smatl, ovrlp_, inv_ovrlp_)
      use module_base
      use sparsematrix_base, only: sparse_matrix, matrices
      use sparsematrix_init, only: matrixindex_in_compressed
      implicit none

      ! Calling arguments
      integer,intent(in) :: iproc, nproc, norb
      integer,dimension(norb),intent(in) :: onwhichatom
      type(sparse_matrix),intent(in) :: smats, smatl
      type(matrices),intent(inout) :: ovrlp_, inv_ovrlp_

      ! Local variables
      integer :: nat, natp, isat, ii, iorb, iiat, n, jorb, jjat, ind, korb
      real(kind=8),dimension(:,:),allocatable :: matrix

      call f_routine(id='calculate_S_minus_one_half_onsite')

      ! Parallelize over the atoms
      nat = maxval(onwhichatom)
      if (iproc<nat) then
          natp = nat/nproc
          ii = nat-nproc*natp
          if (iproc<ii) then
              natp = natp + 1
              isat = iproc*natp
          else
              isat = ii*(natp+1) + (iproc-ii)*natp
          end if
      else
          natp = 0
          isat = nat
      end if


      call f_zero(inv_ovrlp_%matrix_compr)

      iorb = 1
      do
          iiat = onwhichatom(iorb)
          n = 0
          !write(*,*) 'iproc, iorb, iiat', iproc, iorb, iiat
          if (iiat>=isat+1 .and. iiat<=isat+natp) then
              do jorb=iorb,norb
                  jjat = onwhichatom(jorb)
                  if (jjat/=iiat) exit
                  n = n + 1
              end do
              matrix = f_malloc((/n,n/),id='matrix')
              do jorb=iorb,iorb+n-1
                  do korb=iorb,iorb+n-1
                      ind = matrixindex_in_compressed(smats, korb, jorb)
                      matrix(korb-iorb+1,jorb-iorb+1) = ovrlp_%matrix_compr(ind)
                  end do
              end do
              call  overlap_plus_minus_one_half_exact(1,n,-1,.false.,matrix,smats)
              do jorb=iorb,iorb+n-1
                  do korb=iorb,iorb+n-1
                      ind = matrixindex_in_compressed(smatl, korb, jorb)
                      inv_ovrlp_%matrix_compr(ind) = matrix(korb-iorb+1,jorb-iorb+1)
                  end do
              end do
              call f_free(matrix)
          end if
          iorb = iorb + max(n,1)
          if (iorb>norb) exit
      end do

      call mpiallred(inv_ovrlp_%matrix_compr, mpi_sum, comm=bigdft_mpi%mpi_comm)

      call f_release_routine()

    end subroutine calculate_S_minus_one_half_onsite


end module matrix_operations
