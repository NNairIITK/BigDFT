module selinv

    private

    public :: selinv_wrapper

    contains

      !	 Copyright (c) 2012 The Regents of the University of California,
      !	 through Lawrence Berkeley National Laboratory.  
      !
      !  Author: Lin Lin
      !	 
      !  This file is part of PEXSI. All rights reserved.
      !
      !	 Redistribution and use in source and binary forms, with or without
      !	 modification, are permitted provided that the following conditions are met:
      !
      !	 (1) Redistributions of source code must retain the above copyright notice, this
      !	 list of conditions and the following disclaimer.
      !	 (2) Redistributions in binary form must reproduce the above copyright notice,
      !	 this list of conditions and the following disclaimer in the documentation
      !	 and/or other materials provided with the distribution.
      !	 (3) Neither the name of the University of California, Lawrence Berkeley
      !	 National Laboratory, U.S. Dept. of Energy nor the names of its contributors may
      !	 be used to endorse or promote products derived from this software without
      !	 specific prior written permission.
      !
      !	 THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
      !	 ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
      !	 WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
      !	 DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
      !	 ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
      !	 (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
      !	 LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
      !	 ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
      !	 (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
      !	 SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
      !
      !	 You are under no obligation whatsoever to provide any bug fixes, patches, or
      !	 upgrades to the features, functionality or performance of the source code
      !	 ("Enhancements") to anyone; however, if you choose to make your Enhancements
      !	 available either publicly, or directly to Lawrence Berkeley National
      !	 Laboratory, without imposing a separate written license agreement for such
      !	 Enhancements, then you hereby grant the following license: a non-exclusive,
      !	 royalty-free perpetual license to install, use, modify, prepare derivative
      !	 works, incorporate into other computer software, distribute, and sublicense
      !	 such enhancements or derivative works thereof, in binary and source code form.
      !
      !> @file f_driver_ksdft.f90
      !> @brief FORTRAN version of the driver for solving KSDFT.
      !> @date 2014-04-02
      subroutine selinv_driver(iproc, nproc, comm, nfvctr, nvctr, row_ind, col_ptr, mat_s, inv_ovrlp)
      !use module_base
      use pexsi_base, only: f_ppexsi_options
      use pexsi_interfaces
      use yaml_output
      use iso_c_binding
      !use sparsematrix_base
      use wrapper_mpi
      !use time_profiling
      use dynamic_memory
      use futile
      implicit none
      !include 'mpif.h'
      
      ! Calling arguments
      integer,intent(in) :: iproc, nproc, comm
      integer,intent(in) :: nfvctr !< number of row/columns
      integer,intent(in) :: nvctr ! number of non-zero entries
      integer,dimension(nvctr),intent(in) :: row_ind !< row_ind from the CCS format
      integer,dimension(nfvctr),intent(in) :: col_ptr !< col_ptr from the CCS format
      real(kind=8),dimension(nvctr),intent(in) :: mat_s
      real(kind=8),dimension(nvctr),intent(out) :: inv_ovrlp
      
      integer(c_int) :: nrows, nnz, nnzLocal, numColLocal
      integer(c_int), allocatable, dimension(:) ::  colptrLocal, rowindLocal
      real(c_double), allocatable, dimension(:) ::  HnzvalLocal, SnzvalLocal, invSnzvalLocal
      
      integer(c_int):: nprow, npcol, npSymbFact, outputFileIndex, ic
      integer :: ierr, ii
      double precision:: timeSta, timeEnd
      integer(c_int):: info
      integer(c_intptr_t) :: plan
      type(f_ppexsi_options) :: options
      real(kind=8),dimension(:),pointer :: mat_h_local, mat_s_local
      integer,dimension(:),pointer :: col_ptr_local, row_ind_local
      integer :: nfvctr_local_min, nfvctr_local_max
      real(kind=8) :: nfvctr_local_avg
      integer :: nvctr_local_min, nvctr_local_max
      real(kind=8) :: nvctr_local_avg
      
      integer:: i, j , nfvctr_local, nvctr_local, isvctr_local, maxproc
      integer:: numColLocalFirst, firstCol
      integer:: irow, jcol
      
      integer:: readComm
      integer:: isProcRead
      
      call f_routine(id='pexsi_driver')

      if (iproc==0) call yaml_comment('SelInv calculation of inverse',hfill='~')
      

      ! Determine the number of processes used.
      ! Ideal number of processes per pole:
      ii = nproc! / npoles
      ! Number of processor rows / columns
      ii = max(1,floor(sqrt(real(ii,kind=8))))
      nprow = int(ii,kind=c_int)
      npcol = int(ii,kind=c_int)
      ! "Highest" MPI task working
      !maxproc = npoles*nprow*npcol
      maxproc = nprow*npcol

      if (iproc==0) then
          !write(*,*) 'nproc, npoles, nprow, npcol, maxproc', nproc, npoles, nprow, npcol, maxproc
          call yaml_mapping_open('PEXSI parallelization')
          call yaml_map('Total number of cores',nproc)
          call yaml_map('Number of cores used',maxproc)
          call yaml_map('Processor grid dimensions',(/nprow,npcol/))
          call yaml_mapping_close()
      end if
      
      
          
          ! Split the processors to read matrix
          if( iproc < maxproc ) then
            isProcRead = 1
          else
            isProcRead = 0
          endif
          
          call mpi_comm_split( comm, isProcRead, iproc, readComm, ierr )
          
      if (iproc<maxproc) then
          

            call distribute_matrix(iproc, maxproc, nfvctr, nvctr, col_ptr, row_ind, mat_s, &
                 nfvctr_local, nvctr_local, isvctr_local, col_ptr_local, row_ind_local, mat_s_local)
            !Because they are C integers..?
            nrows = int(nfvctr,kind=c_int)
            nnz = int(nvctr,kind=c_int)
            nnzLocal = int(nvctr_local,kind=c_int)
            numColLocal = int(nfvctr_local,kind=c_int)

            nfvctr_local_min = nfvctr_local
            nfvctr_local_max = nfvctr_local
            call mpiallred(nfvctr_local_min,count=1,op=mpi_min,comm=readComm)
            call mpiallred(nfvctr_local_max,count=1,op=mpi_max,comm=readComm)
            nfvctr_local_avg = real(nfvctr,kind=8)/real(maxproc,kind=8)
            nvctr_local_min = nvctr_local
            nvctr_local_max = nvctr_local
            call mpiallred(nvctr_local_min,count=1,op=mpi_min,comm=readComm)
            call mpiallred(nvctr_local_max,count=1,op=mpi_max,comm=readComm)
            nvctr_local_avg = real(nvctr,kind=8)/real(maxproc,kind=8)
          
            if(iproc== 0) then
              call yaml_mapping_open('Dimensions of the matrices')
              call yaml_mapping_open('Global')
              call yaml_map('number of columns',nfvctr)
              call yaml_map('number of non-zero elements)',nvctr)
              call yaml_mapping_close()
              call yaml_mapping_open('Distributed')
              call yaml_map('number of columns (min/max/avg)', &
                   (/real(nfvctr_local_min,kind=8),real(nfvctr_local_max,kind=8),nfvctr_local_avg/),fmt='(f7.1)')
              call yaml_map('number of non-zero elements (min/max/avg)', &
                   (/real(nvctr_local_min,kind=8),real(nvctr_local_max,kind=8),nvctr_local_avg/),fmt='(f7.1)')
              call yaml_mapping_close()
              call yaml_mapping_close()
            endif
          
            ! Allocate memory
            colptrLocal = f_malloc(numColLocal+1,id='colptrLocal')
            rowindLocal = f_malloc(nnzLocal,id='rowindLocal')
            HnzvalLocal = f_malloc(nnzLocal,id='HnzvalLocal')
            SnzvalLocal = f_malloc(nnzLocal,id='SnzvalLocal')
            invSnzvalLocal = f_malloc(nnzLocal,id='invSnzvalLocal')
          

            do i=1,nfvctr_local+1
                ic = int(i,kind=c_int)
                colptrLocal(ic) = int(col_ptr_local(i),kind=c_int)
            end do
            do i=1,nvctr_local
                ic = int(i,kind=c_int)
                rowIndLocal(ic) = int(row_ind_local(i),kind=c_int)
                SnzvalLocal(ic) = real(mat_s_local(i),kind=c_double)
            end do

            !write(*,*) 'OLD: HnzvalLocal',HnzvalLocal
            !write(*,*) 'OLD: iproc, SnzvalLocal',iproc,SnzvalLocal
          
          !endif
          
          ! Step 1. Initialize
          
          
          if( mod( iproc, nprow * npcol ) .eq. 0 ) then
            outputFileIndex = iproc / (nprow * npcol);
          else
            outputFileIndex = -1;
          endif
          
          
          plan = f_ppexsi_plan_initialize(&
            readComm,&
            nprow,&
            npcol,&
            outputFileIndex,&
            info )
          
          if( info .ne. 0 ) then
              call mpi_finalize( ierr )
              call exit(info)
          endif
          
          call f_ppexsi_set_default_options(&
            options )
          
          
          call f_ppexsi_load_real_symmetric_hs_matrix(&
                plan,&       
                options,&
                nrows,&
                nnz,&
                nnzLocal,&
                numColLocal,&
                colptrLocal,& 
                rowindLocal,&
                HnzvalLocal,&
                int(0,kind=c_int),&
                SnzvalLocal,&
                info ) 
          
          if( info .ne. 0 ) then
            call mpi_finalize( ierr )
            call exit(info)
          endif
          
          
            
          call f_ppexsi_symbolic_factorize_real_symmetric_matrix(&                                              
            plan,&                                                                                              
            options,&
            info)                                                                                               
          
          if( info .ne. 0 ) then                                                                                
              call mpi_finalize( ierr )                                                                         
              call exit(info)
          endif
            
            
            
          call f_ppexsi_selinv_real_symmetric_matrix(&                                                          
            plan,& 
            options,&
            SnzvalLocal,&
            invSnzvalLocal,&                                                                                    
            info)                                                                                               
                                                                                                                
          if( info .ne. 0 ) then
              call mpi_finalize( ierr )
              call exit(info)                                                                                   
          endif                                                                                                 
          
            
          
          
          ! Step 3. Clean up */
          
          call f_ppexsi_plan_finalize( plan, info )
          
          
          ! Gather the local copies of the kernel
          call f_zero(inv_ovrlp)
          do i=1,nvctr_local
              !write(*,*) 'iproc, i, isvctr_local+i-1', iproc, i, isvctr_local+i-1
              inv_ovrlp(isvctr_local+i-1) = invSnzvalLocal(i)
          end do
          call mpiallred(inv_ovrlp,mpi_sum,comm=readComm)

          
          if( isProcRead == 1 ) then
            call f_free(colptrLocal)
            call f_free(rowindLocal)
            call f_free(HnzvalLocal)
            call f_free(SnzvalLocal)
            call f_free(invSnzvalLocal)
          endif

        call f_free_ptr(col_ptr_local)
        call f_free_ptr(row_ind_local)
        call f_free_ptr(mat_h_local)
        call f_free_ptr(mat_s_local)

      end if

          call mpi_comm_free( readComm, ierr )

      call mpibcast(inv_ovrlp, root=0, comm=comm) 

      if (iproc==0) call yaml_comment('SelInv calculation of inverse finished',hfill='~')
      
      call f_release_routine()
      
      end subroutine selinv_driver



      subroutine distribute_matrix(iproc, nproc, nfvctr, nvctr, col_ptr, row_ind, mat_S, &
                 nfvctr_local, nvctr_local, isvctr_local, col_ptr_local, row_ind_local, mat_s_local)
        !use module_base
        use dynamic_memory
        use wrapper_mpi
        implicit none

        ! Calling arguments
        integer,intent(in) :: iproc !< ID of MPI task
        integer,intent(in) :: nproc !< number of MPI tasks
        integer,intent(in) :: nfvctr !< number of rows/columns (global)
        integer,intent(in) :: nvctr !< number of nonzero elements (global)
        integer,dimension(nfvctr),intent(in) :: col_ptr !< col_ptr from the CCS format
        integer,dimension(nvctr),intent(in) :: row_ind !< row_ind from the CCS format
        real(kind=8),dimension(nvctr),intent(in) :: mat_s !< values of S
        integer,intent(out) :: nfvctr_local !< number of columns (local)
        integer,intent(out) :: nvctr_local !< number of nonzero elements (local)
        integer,intent(out) :: isvctr_local !< global index of the first local non-zero element for each MPI task
        integer,dimension(:),pointer,intent(out) :: col_ptr_local !< local part of col_ptr
        integer,dimension(:),pointer,intent(out) :: row_ind_local !< local part of row_ind
        real(kind=8),dimension(:),pointer,intent(out) :: mat_s_local !< local entries of S

        ! Local variables
        integer :: ncol, icol, ii, i
        integer,dimension(:),allocatable :: col_ptr_tmp

        call f_routine(id='distribute_matrix')
        
        ! Distribute the columns, the last process takes what remains
        ncol = nfvctr/nproc
        if (iproc<nproc-1) then
            nfvctr_local = ncol
        else
            nfvctr_local = nfvctr - ncol*(nproc-1)
        end if

        ! Temporary copy of col_ptr, with an additional entry at the end
        col_ptr_tmp = f_malloc(nfvctr+1,id='col_ptr_tmp')
        do icol=1,nfvctr
            col_ptr_tmp(icol) = col_ptr(icol)
        end do
        col_ptr_tmp(nfvctr+1) = nvctr + 1

        ! Distribute the number non-zero entries
        ii = ncol*iproc+1
        nvctr_local = col_ptr_tmp(ii+nfvctr_local) - col_ptr_tmp(ii)

        col_ptr_local = f_malloc_ptr(nfvctr_local+1,id='col_ptr_local')
        row_ind_local = f_malloc_ptr(nvctr_local,id='row_ind_local')
        mat_s_local = f_malloc_ptr(nvctr_local,id='mat_s_local')

        ! Local copy of col_ptr, with respect to the first local element
        do i=1,nfvctr_local+1
            col_ptr_local(i) = col_ptr_tmp(iproc*ncol+i) - col_ptr_tmp(iproc*ncol+1) + 1
        end do

        ! Local copy of row_ind
        ii = col_ptr_tmp(iproc*ncol+1)-1
        do i=1,nvctr_local
            row_ind_local(i) = row_ind(ii+i)
            mat_s_local(i) = mat_s(ii+i)
        end do

        ! Global index of the first local element
        isvctr_local = col_ptr_tmp(ncol*iproc+1)


        call f_free(col_ptr_tmp)
        
        call f_release_routine()

      end subroutine distribute_matrix



   !!   subroutine distribute_matrix

   !!     ! Calling arguments
   !!     integer, intent(in) :: iproc, nproc, nfvctr, nvctr, nfvctr_local, nvctr_local
   !!     integer,dimension(nfvctr),intent(in) :: col_ptr
   !!     integer,dimension(nvctr),intent(in) :: row_ind
   !!     real(kind=8),dimesion(nvctr),intent(in) ::  mat_h, mat_s
   !!     integer,dimension(nfvctr_local+1),intent(out) :: col_ptr_local
   !!     integer,dimension(nvctr_local),intent(out) :: row_ind_local
   !!     real(kind=8),dimension(nvctr_local),intent(out) :: mat_h_local, mat_s_local

   !!     ! Local variables


   !!     if (iproc<nproc-1) then
   !!         do i=1,nfvctr_local+1
   !!             col_ptr_local(i) = col_ptr(iproc*nfvctr_local+i)
   !!         end do
   !!     else
   !!         do i=1,nfvctr_local
   !!             col_ptr_local(i) = col_ptr(iproc*nfvctr_local+i)
   !!         end do
   !!         col_ptr_local(nfvctr_local+1) = nvctr
   !!     end if

   !!     do i=1,nvctr_local
   !!         row_ind_local(i) = row_ind(iproc*nvctr_local)
   !!     end do

   !!   end subroutine distribute_matrix


    subroutine selinv_wrapper(iproc, nproc, comm, smats, smatl, ovrlp, inv_ovrlp)
      use futile
      use sparsematrix_base
      use sparsematrix, only: transform_sparse_matrix
      use sparsematrix_init, only: sparsebigdft_to_ccs
      implicit none
    
      ! Calling arguments
      integer,intent(in) :: iproc, nproc, comm
      type(sparse_matrix),intent(in) :: smats, smatl
      type(matrices),intent(in) :: ovrlp
      type(matrices),intent(out) :: inv_ovrlp
    
      ! Local variables
      integer,dimension(:),allocatable :: row_ind, col_ptr
      real(mp),dimension(:),allocatable :: ovrlp_large
    
      call f_routine(id='selinv_wrapper')
    
      !call write_pexsi_matrices(iproc, nproc, smatm, smats, ham%matrix_compr, ovrlp%matrix_compr)
      row_ind = f_malloc(smatl%nvctr,id='row_ind')
      col_ptr = f_malloc(smatl%nfvctr,id='col_ptr')
      call sparsebigdft_to_ccs(smatl%nfvctr, smatl%nvctr, smatl%nseg, smatl%keyg, row_ind, col_ptr)
      ! At the moment not working for nspin>1
      ovrlp_large = sparsematrix_malloc(smatl, iaction=SPARSE_FULL, id='ovrlp_large')
      if (smats%ntaskgroup/=1 .or. smatl%ntaskgroup/=1) then
          call f_err_throw('PEXSI is not yet tested with matrix taskgroups', err_name='SPARSEMATRIX_RUNTIME_ERROR')
      end if
      call transform_sparse_matrix(iproc, smats, smatl, SPARSE_FULL, 'small_to_large', &
           smat_in=ovrlp%matrix_compr, lmat_out=ovrlp_large)
      call selinv_driver(iproc, nproc, comm, smatl%nfvctr, smatl%nvctr, row_ind, col_ptr, &
           ovrlp_large, inv_ovrlp%matrix_compr)
    
      call f_free(ovrlp_large)
      call f_free(row_ind)
      call f_free(col_ptr)
    
      call f_release_routine()
    
    end subroutine selinv_wrapper


end module selinv
