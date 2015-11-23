module pexsi

    private

    public :: pexsi_driver

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
      subroutine pexsi_driver(Hfile, Sfile, charge, npoles, mumin, mumax, mu, temperature, tol_charge, &
                 nsize_kernel, kernel, energy)
      use module_base
      use f_ppexsi_interface
      use iso_c_binding
      implicit none
      !include 'mpif.h'
      
      ! Calling arguments
      character(len=*),intent(in) :: Hfile, Sfile
      real(kind=8),intent(in) :: charge
      integer,intent(in) :: npoles, nsize_kernel
      real(kind=8),intent(in) :: mumin, mumax, mu, temperature, tol_charge
      real(kind=8),dimension(nsize_kernel),intent(out) :: kernel
      real(kind=8),intent(out) :: energy
      
      integer(c_int) :: nrows, nnz, nnzLocal, numColLocal
      integer(c_int), allocatable, dimension(:) ::  colptrLocal, rowindLocal
      real(c_double), allocatable, dimension(:) ::  &
        HnzvalLocal, SnzvalLocal
      real(c_double), allocatable, dimension(:) ::  &
        DMnzvalLocal, EDMnzvalLocal, FDMnzvalLocal
      real(c_double) :: numElectronExact, muPEXSI, numElectronPEXSI, &
        muMinInertia, muMaxInertia
      integer(c_int) :: numTotalInertiaIter, numTotalPEXSIIter
      real(c_double) :: totalEnergyH, totalEnergyS, totalFreeEnergy
      
      integer(c_int):: nprow, npcol, npSymbFact, outputFileIndex
      integer :: iproc, nproc, ierr
      double precision:: timeSta, timeEnd
      integer(c_int):: info
      integer(c_intptr_t) :: plan
      type(f_ppexsi_options) :: options
      
      integer:: i, j 
      integer:: numColLocalFirst, firstCol
      integer:: irow, jcol
      
      integer:: readComm
      integer:: isProcRead
      
      call f_routine(id='f_driver_ksdft')
      
      !call mpi_init( ierr )
      call mpi_comm_rank( MPI_COMM_WORLD, iproc, ierr )
      call mpi_comm_size( MPI_COMM_WORLD, nproc, ierr )
      
      if (nproc/=1) then
          call f_err_throw('pexsi not yet ready in parallel')
      end if
      
      !Hfile            = "lap2dr.matrix"
      !Sfile            = "lap2dr.matrix"
      !Hfile            = "hamiltonian_sparse_PEXSI.bin"
      !Sfile            = "overlap_sparse_PEXSI.bin"
      
      ! Only use 4 processors in this example
      nprow = 1
      npcol = 1
      
      ! Split the processors to read matrix
      if( iproc < nprow * npcol ) then
      	isProcRead = 1
      else
      	isProcRead = 0
      endif
      
      call mpi_comm_split( MPI_COMM_WORLD, isProcRead, iproc, readComm, ierr )
      
      
      if( isProcRead == 1 ) then
        call f_read_distsparsematrix_formatted_head( &
          trim(Hfile)//char(0),&
          nrows,&
          nnz,&
          nnzLocal,&
          numColLocal,&
          readComm )
      
        if( iproc .eq. 0 ) then
          write(*,*) "Matrix size (local data on proc 0):" 
          write(*,*) "size = ", nrows
          write(*,*) "nnz  = ", nnz
          write(*,*) "nnzLocal = ", nnzLocal
          write(*,*) "numColLocal = ", numColLocal
        endif
      
        ! Allocate memory
        allocate( colptrLocal( numColLocal + 1 ) )
        allocate( rowindLocal( nnzLocal ) )
        allocate( HnzvalLocal( nnzLocal ) )
        allocate( SnzvalLocal( nnzLocal ) )
        allocate( DMnzvalLocal( nnzLocal ) )
        allocate( EDMnzvalLocal( nnzLocal ) )
        allocate( FDMnzvalLocal( nnzLocal ) )
      
        call f_read_distsparsematrix_formatted (&
          trim(Hfile)//char(0),&
          nrows,&
          nnz,&
          nnzLocal,&
          numColLocal,&
          colptrLocal,&
          rowindLocal,&
          HnzvalLocal,&
          readComm )
      
        call f_read_distsparsematrix_formatted (&
          trim(Sfile)//char(0),&
          nrows,&
          nnz,&
          nnzLocal,&
          numColLocal,&
          colptrLocal,&
          rowindLocal,&
          SnzvalLocal,&
          readComm )
      
      endif
      
      ! Step 1. Initialize PEXSI 
      
      ! Set the outputFileIndex to be the pole index.
      ! The first processor for each pole outputs information
      
      if( mod( iproc, nprow * npcol ) .eq. 0 ) then
        outputFileIndex = iproc / (nprow * npcol);
      else
        outputFileIndex = -1;
      endif
      
      
      plan = f_ppexsi_plan_initialize(&
        MPI_COMM_WORLD,&
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
      
      numElectronExact = charge
      options%muMin0   = mumin
      options%muMax0   = mumax
      options%mu0      = mu
      options%deltaE   = 20d0 
      options%numPole  = npoles
      options%temperature = temperature
      options%muPEXSISafeGuard = 0.2d0
      options%numElectronPEXSITolerance = tol_charge
      
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
            0,&
            SnzvalLocal,&
            info ) 
      
      if( info .ne. 0 ) then
      	call mpi_finalize( ierr )
      	call exit(info)
      endif
      
      
      if( iproc == 0 ) then
        write(*,*)  "Finish setting up the matrix."
      endif
      
      ! Step 2. PEXSI Solve
      
      call f_ppexsi_dft_driver(&
        plan,&
        options,&
        numElectronExact,&
        muPEXSI,&
        numElectronPEXSI,&
        muMinInertia,&
        muMaxInertia,&
        numTotalInertiaIter,&
        numTotalPEXSIIter,&
        info)
      
      if( info .ne. 0 ) then
      	call mpi_finalize( ierr )
      	call exit(info)
      endif
      
      
      if( iproc == 0 ) then
        write(*,*)  "Finish DFT driver."
      endif
      
      if( isProcRead == 1 ) then
        call f_ppexsi_retrieve_real_symmetric_dft_matrix(&
          plan,&
          DMnzvalLocal,&
          EDMnzvalLocal,&
          FDMnzvalLocal,&
          totalEnergyH,&
          totalEnergyS,&
          totalFreeEnergy,&
          info)
      
        if( iproc == 0 ) then
          write(*,*) "Output from the main program."
          write(*,*) "Total energy (H*DM)         = ", totalEnergyH
          write(*,*) "Total energy (S*EDM)        = ", totalEnergyS
          write(*,*) "Total free energy           = ", totalFreeEnergy
        endif
      endif
      
      energy = totalEnergyH
      
      
      ! Step 3. Clean up */
      
      call f_ppexsi_plan_finalize( plan, info )
      
      call mpi_comm_free( readComm, ierr )
      !call mpi_finalize( ierr )
      
      kernel(:) = DMnzvalLocal
      
      if( isProcRead == 1 ) then
        deallocate( colptrLocal )
        deallocate( rowindLocal )
        deallocate( HnzvalLocal )
        deallocate( DMnzvalLocal )
        deallocate( EDMnzvalLocal )
        deallocate( FDMnzvalLocal )
      endif
      
      call f_routine(id='f_driver_ksdft')
      
      end subroutine pexsi_driver

end module pexsi
