!> @file
!!  Routines related to the initialization of the Kohn-Sham wavefunctions
!! @author
!!    Copyright (C) 2013-2013 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS



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
  use memory_profiling
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


subroutine kswfn_emit_psi(Wfn, iter, psi_or_hpsi, iproc, nproc)
  use module_base
  use module_types
  implicit none
  type(DFT_wavefunction), intent(in) :: Wfn
  integer, intent(in) :: iter, iproc, nproc, psi_or_hpsi

  integer, parameter :: SIGNAL_DONE = -1
  integer, parameter :: SIGNAL_WAIT = -2
  integer :: message, ierr, data(2)
  integer :: status(MPI_STATUS_SIZE)

  call timing(iproc,'wf_signals    ','ON')
  if (iproc == 0) then
     ! Only iproc 0 emit the signal. This call is blocking.
     ! All other procs are blocked by the bcast to wait for
     ! possible transfer to proc 0.
     if (psi_or_hpsi == 0) then
        call wf_emit_psi(Wfn%c_obj, iter)
     else
        call wf_emit_hpsi(Wfn%c_obj, iter)
     end if
     if (nproc > 1) then
        ! After handling the signal, iproc 0 broadcasts to other
        ! proc to continue (jproc == -1).
        message = SIGNAL_DONE
        call MPI_BCAST(message, 1, MPI_INTEGER, 0, bigdft_mpi%mpi_comm, ierr)
     end if
  else
     message = SIGNAL_WAIT
     do
        if (message == SIGNAL_DONE) then
           exit
        end if
        call MPI_BCAST(message, 1, MPI_INTEGER, 0, bigdft_mpi%mpi_comm, ierr)
        
        if (message > 0 .and. iproc == message) then
           ! Will have to send to iproc 0 some of psi.
           call MPI_RECV(data, 2, MPI_INTEGER, 0, 123, bigdft_mpi%mpi_comm, status, ierr)
           if (psi_or_hpsi == 0) then
              call MPI_SEND(Wfn%psi(1 + data(1)), data(2), MPI_DOUBLE_PRECISION, &
                   & 0, 123, bigdft_mpi%mpi_comm, ierr)
           else
              call MPI_SEND(Wfn%hpsi(1 + data(1)), data(2), MPI_DOUBLE_PRECISION, &
                   & 0, 123, bigdft_mpi%mpi_comm, ierr)
           end if
        end if
     end do
  end if
  call timing(iproc,'wf_signals    ','OF')
END SUBROUTINE kswfn_emit_psi


subroutine kswfn_mpi_copy(psic, jproc, psiStart, psiSize)
  use module_base
  use module_types
  implicit none
  integer, intent(in) :: psiSize, jproc, psiStart
  real(wp), intent(inout) :: psic(psiSize)

  integer :: ierr
  integer :: status(MPI_STATUS_SIZE)

  if (jproc == 0) return

  call MPI_BCAST(jproc, 1, MPI_INTEGER, 0, bigdft_mpi%mpi_comm, ierr)

  call MPI_SEND((/ psiStart, psiSize /), 2, MPI_INTEGER, jproc, 123, bigdft_mpi%mpi_comm, ierr)
  call MPI_RECV(psic, psiSize, MPI_DOUBLE_PRECISION, jproc, 123, bigdft_mpi%mpi_comm, status, ierr)
END SUBROUTINE kswfn_mpi_copy


subroutine kswfn_init_comm(wfn, in, atoms, dpbox, iproc, nproc)
  use module_types
  use module_interfaces, except_this_one => kswfn_init_comm
  use communications, only: collective_comms_null, init_collective_comms, init_collective_comms_sumrho
  implicit none
  integer, intent(in) :: iproc, nproc
  type(DFT_wavefunction), intent(inout) :: wfn
  type(input_variables), intent(in) :: in
  type(atoms_data),intent(in) :: atoms
  type(denspot_distribution), intent(in) :: dpbox

  ! Nullify all pointers
  nullify(wfn%psi)
  nullify(wfn%hpsi)
  nullify(wfn%psit)
  nullify(wfn%psit_c)
  nullify(wfn%psit_f)
  nullify(wfn%spsi)
  nullify(wfn%gaucoeffs)

  call initialize_communication_potential(iproc, nproc, dpbox%nscatterarr, &
       & wfn%orbs, wfn%lzd, wfn%comgp)

  !call nullify_collective_comms(wfn%collcom)
  !call nullify_collective_comms(wfn%collcom_sr)
  wfn%collcom=collective_comms_null()
  wfn%collcom_sr=collective_comms_null()

  call init_collective_comms(iproc, nproc, wfn%npsidim_orbs, wfn%orbs, wfn%lzd, wfn%collcom)
  call init_collective_comms_sumrho(iproc, nproc, wfn%lzd, wfn%orbs, dpbox%nscatterarr, wfn%collcom_sr)

END SUBROUTINE kswfn_init_comm


subroutine kswfn_emit_lzd(Wfn, iproc, nproc)
  use module_base
  use module_types
  implicit none
  type(DFT_wavefunction), intent(in) :: Wfn
  integer, intent(in) :: iproc, nproc

  call timing(iproc,'wf_signals    ','ON')
  if (iproc == 0) then
     call wf_emit_lzd(Wfn%c_obj)
  end if
  call timing(iproc,'wf_signals    ','OF')
END SUBROUTINE kswfn_emit_lzd
