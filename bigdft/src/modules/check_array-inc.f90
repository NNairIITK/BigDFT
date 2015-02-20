
      maxdiff=0.0_gp

      if (nproc == 1 .or. ndims == 0) return

      !check that the positions are identical for all the processes
      rxyz_glob=f_malloc((/ndims,nproc/),id='rxyz_glob')

      !gather the results for all the processors
      call MPI_GATHER(array,ndims,mpidtypg,&
           rxyz_glob,ndims,mpidtypg,0,mpi_comm,ierr)
      do jproc=2,nproc
         do i=1,ndims
            maxdiff=max(maxdiff,&
                 abs(rxyz_glob(i,jproc)-rxyz_glob(i,1)))
         end do
      end do

      call f_free(rxyz_glob)
