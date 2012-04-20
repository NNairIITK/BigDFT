subroutine glr_get_psi_size(glr, psisize)
  use module_types
  implicit none
  type(locreg_descriptors), intent(in) :: glr
  integer, intent(out) :: psisize

  psisize = glr%wfd%nvctr_c + 7 * glr%wfd%nvctr_f
END SUBROUTINE glr_get_psi_size

subroutine kswfn_free_scf_data(KSwfn, freePsit)
  use module_base
  use module_types
  use m_profiling
  implicit none
  type(DFT_wavefunction), intent(inout) :: KSwfn
  logical, intent(in) :: freePsit
  
  character(len = *), parameter :: subname = "kswfn_free_scf_data"
  integer :: i_all, i_stat

  ! Clean KSwfn parts only needed in the SCF loop.
  call deallocate_diis_objects(KSwfn%diis,subname)
  i_all=-product(shape(KSwfn%hpsi))*kind(KSwfn%hpsi)
  deallocate(KSwfn%hpsi,stat=i_stat)
  call memocc(i_stat,i_all,'hpsi',subname)
  if (freePsit) then
     i_all=-product(shape(KSwfn%psit))*kind(KSwfn%psit)
     deallocate(KSwfn%psit,stat=i_stat)
     call memocc(i_stat,i_all,'psit',subname)
  else
     nullify(KSwfn%psit)
  end if
end subroutine kswfn_free_scf_data

subroutine kswfn_emit_psi(KSwfn, iter, iproc, nproc)
  use module_base
  use module_types
  implicit none
  type(DFT_wavefunction), intent(in) :: KSwfn
  integer, intent(in) :: iter, iproc, nproc

  integer, parameter :: SIGNAL_DONE = -1
  integer, parameter :: SIGNAL_WAIT = -2
  integer :: message, ierr, data(2), orbSize
  integer :: status(MPI_STATUS_SIZE)

  call timing(iproc,'wf_signals    ','ON')
  if (iproc == 0) then
     ! Only iproc 0 emit the signal. This call is blocking.
     ! All other procs are blocked by the bcast to wait for
     ! possible transfer to proc 0.
     call wf_emit_psi(KSwfn%c_obj, iter)
     if (nproc > 1) then
        ! After handling the signal, iproc 0 broadcasts to other
        ! proc to continue (jproc == -1).
        message = SIGNAL_DONE
        call MPI_BCAST(message, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
     end if
  else
     message = SIGNAL_WAIT
     do
        if (message == SIGNAL_DONE) then
           exit
        end if
        call MPI_BCAST(message, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
        
        if (message > 0 .and. iproc == message) then
           ! Will have to send to iproc 0 some of psi.
           call MPI_RECV(data, 2, MPI_INTEGER, 0, 123, MPI_COMM_WORLD, status, ierr)
           call glr_get_psi_size(KSwfn%Lzd%Glr, orbSize)
           call MPI_SEND(KSwfn%psi(1 + data(1) * orbSize), data(2), MPI_DOUBLE_PRECISION, &
                & 0, 123, MPI_COMM_WORLD, ierr)
        end if
     end do
  end if
  call timing(iproc,'wf_signals    ','OF')
END SUBROUTINE kswfn_emit_psi

subroutine kswfn_mpi_copy(psic, jproc, iorbp, psiSize)
  use module_base
  use module_types
  implicit none
  integer, intent(in) :: psiSize, jproc, iorbp
  real(wp), intent(inout) :: psic(psiSize)

  integer :: ierr
  integer :: status(MPI_STATUS_SIZE)

  if (jproc == 0) return

  call MPI_BCAST(jproc, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

  call MPI_SEND((/ iorbp, psiSize /), 2, MPI_INTEGER, jproc, 123, MPI_COMM_WORLD, ierr)
  call MPI_RECV(psic, psiSize, MPI_DOUBLE_PRECISION, jproc, 123, MPI_COMM_WORLD, status, ierr)
END SUBROUTINE kswfn_mpi_copy
