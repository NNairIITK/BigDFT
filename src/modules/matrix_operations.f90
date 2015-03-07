module matrix_operations
    use module_base
    implicit none

    private

    !> Public routines
    public :: overlapPowerGeneral


    contains


      !> S^-1 exact only works for symmetric matrices
      !! BOTH sparse matrices must be present together and inv_ovrlp should be nullified pointer, NOT inv_ovrlp_smat%matrix
      !! when sparse matrices present, check is performed to see whether %matrix is allocated so that its allocated status remains unchanged
      !! contents of %matrix not guaranteed to be correct though - inv_ovrlp_smat%can_use_dense set accordingly
      !! power: -2 -> S^-1/2, 2 -> S^1/2, 1 -> S^-1
      subroutine overlapPowerGeneral(iproc, nproc, iorder, ncalc, power, blocksize, imode, &
                 ovrlp_smat, inv_ovrlp_smat, ovrlp_mat, inv_ovrlp_mat, check_accur, &
                 max_error, mean_error, nspinx)
           !!foe_nseg, foe_kernel_nsegline, foe_istsegline, foe_keyg)
        use module_base
        use module_types
        use module_interfaces
        use sparsematrix_base, only: sparse_matrix, &
                                sparsematrix_malloc_ptr, sparsematrix_malloc, sparsematrix_malloc0, sparsematrix_malloc0_ptr, &
                                assignment(=), &
                                SPARSE_FULL, DENSE_PARALLEL, SPARSE_MATMUL_LARGE, &
                                DENSE_MATMUL, DENSE_FULL, SPARSEMM_SEQ, SPARSE_TASKGROUP
        use sparsematrix_init, only: get_line_and_column
        use sparsematrix, only: compress_matrix, uncompress_matrix, &
                                transform_sparse_matrix, transform_sparse_matrix_local, &
                                compress_matrix_distributed_wrapper, &
                                uncompress_matrix_distributed2, &
                                sequential_acces_matrix_fast2, sequential_acces_matrix_fast, &
                                gather_matrix_from_taskgroups, gather_matrix_from_taskgroups_inplace, &
                                uncompress_matrix2, transform_sparsity_pattern, &
                                sparsemm_new
        use yaml_output
        implicit none
        
        ! Calling arguments
        integer,intent(in) :: iproc, nproc, iorder, blocksize, ncalc
        integer,dimension(ncalc),intent(in) :: power
        integer,intent(in) :: imode
        type(sparse_matrix),intent(inout) :: ovrlp_smat, inv_ovrlp_smat
        type(matrices),intent(inout) :: ovrlp_mat
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
        real(kind=8) :: factor, newfactor
        logical :: ovrlp_allocated, inv_ovrlp_allocated
      
        ! new for sparse taylor
        integer :: nout, nseq, ispin, ishift, ishift2, isshift, ilshift, ilshift2, nspin, iline, icolumn, ist
        integer,dimension(:,:,:),allocatable :: istindexarr
        real(kind=8),dimension(:),pointer :: ovrlpminone_sparse
        real(kind=8),dimension(:),allocatable :: ovrlp_compr_seq, ovrlpminone_sparse_seq, ovrlp_large_compr, tmparr
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
        real(kind=8) :: ex, max_error_p, mean_error_p
        real(kind=8),dimension(:),allocatable :: factor_arr
        real(kind=8),dimension(:),allocatable :: ovrlp_largep_new, invovrlpp_new
        real(kind=8),dimension(:),allocatable :: Amat12p_new, Amat21p_new
        real(kind=8),dimension(:),pointer :: Amat11p_new, Amat22p_new
      
      
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
                           call overlap_minus_one_exact_serial(ovrlp_smat%nfvctr,inv_ovrlp_mat(1)%matrix(1,1,ispin))
                       end do
                   else
                      stop 'check if working - upper half may not be filled'
                      call dpotrf_parallel(iproc, nproc, blocksize, bigdft_mpi%mpi_comm, 'l', &
                           ovrlp_smat%nfvctr, inv_ovrlp_mat(1)%matrix(1,1,1), ovrlp_smat%nfvctr)
                      call dpotri_parallel(iproc, nproc, blocksize, bigdft_mpi%mpi_comm, 'l', &
                           ovrlp_smat%nfvctr, inv_ovrlp_mat(1)%matrix(1,1,1), ovrlp_smat%nfvctr)
                   end if
                else if (power(1)==2) then
                    do ispin=1,nspin
                        call overlap_plus_minus_one_half_exact(bigdft_mpi%nproc,ovrlp_smat%nfvctr, &
                             blocksize,.true.,inv_ovrlp_mat(1)%matrix(1,1,ispin),inv_ovrlp_smat)
                    end do
                else if (power(1)==-2) then
                    do ispin=1,nspin 
                        call overlap_plus_minus_one_half_exact(bigdft_mpi%nproc,ovrlp_smat%nfvctr, &
                             blocksize,.false.,inv_ovrlp_mat(1)%matrix(1,1,ispin),inv_ovrlp_smat)
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
                if (iorder>1) then
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
                                          ovrlp_mat%matrix(1,ovrlp_smat%isfvctr+1,ispin),ovrlpminonep(1,1,ispin))
      
      
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
                        ovrlp_mat%matrix(1,isorb+1,ispin),inv_ovrlpp)
                    !!do iorb=1,norbp
                    !!    do jorb=1,ovrlp_smat%nfvctr
                    !!        write(2900+10*iproc+ispin,'(a,3i8,3es14.6)') 'ispin, iorb, jorb, vals', &
                    !!             ispin, iorb, jorb, inv_ovrlpp(jorb,iorb), ovrlp_mat%matrix(jorb,isorb+iorb,ispin)
                    !!    end do
                    !!end do
      
                    do i=2,iorder
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
                        if (i/=iorder.and.norbp>0) then
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
      
                if (iorder>1) then
                    if(nproc > 1) then
                        nullify(ovrlpminone)
                    else
                        call f_free_ptr(ovrlpminone)
                    end if
                end if
      
                if (iorder>1) then
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
                              call overlap_minus_one_exact_serial(ovrlp_smat%nfvctr,inv_ovrlp_local(1,1,ispin))
                           else
                              stop 'check if working - upper half may not be filled'
                              call dpotrf_parallel(iproc, nproc, blocksize, bigdft_mpi%mpi_comm, 'l', &
                                   ovrlp_smat%nfvctr, inv_ovrlp_local(1,1,ispin), ovrlp_smat%nfvctr)
                              call dpotri_parallel(iproc, nproc, blocksize, bigdft_mpi%mpi_comm, 'l', &
                                   ovrlp_smat%nfvctr, inv_ovrlp_local(1,1,ispin), ovrlp_smat%nfvctr)
                           end if
                        else if (power(icalc)==2) then
                            call overlap_plus_minus_one_half_exact(bigdft_mpi%nproc,ovrlp_smat%nfvctr, &
                                 blocksize,.true.,inv_ovrlp_local(1,1,ispin),inv_ovrlp_smat)
                        else if (power(icalc)==-2) then
                            call overlap_plus_minus_one_half_exact(bigdft_mpi%nproc,ovrlp_smat%nfvctr, &
                                 blocksize,.false.,inv_ovrlp_local(1,1,ispin),inv_ovrlp_smat)
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
      
                call transform_sparse_matrix_local(ovrlp_smat, inv_ovrlp_smat, &
                     ovrlp_mat%matrix_compr, Amat12_compr, 'small_to_large')
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
                    call ice(iproc, nproc, iorder-1000, ovrlp_smat, inv_ovrlp_smat, ncalc, power, ovrlp_mat, inv_ovrlp_mat)
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
                    call transform_sparse_matrix_local(ovrlp_smat, inv_ovrlp_smat, &
                         ovrlp_mat%matrix_compr, ovrlp_large_compr, 'small_to_large')
                !!end if
                invovrlp_compr_seq = sparsematrix_malloc(inv_ovrlp_smat, iaction=SPARSEMM_SEQ, id='ovrlp_large_compr_seq')
                !!ovrlp_largep = sparsematrix_malloc(inv_ovrlp_smat, iaction=DENSE_MATMUL, id='ovrlp_largep')
                ovrlp_largep_new = f_malloc(inv_ovrlp_smat%smmm%nvctrp,id='ovrlp_largep')
      
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

end module matrix_operations
