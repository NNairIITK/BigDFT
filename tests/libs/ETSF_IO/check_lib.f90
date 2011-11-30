program schtroumpf_lib
  use BigDFT_API

  implicit none

  character(len = *), parameter :: filename = "data/wavefunction.etsf"
  integer, parameter :: iorbp = 23
  real(wp), dimension(:,:,:,:), pointer :: psiscf
  real(gp) :: hx, hy, hz
  integer :: n1, n2, n3, nspinor
  
  real(wp) :: nrm
  integer :: ierr, i, j, k  

  call MPI_INIT(ierr)

  write(*,"(3A)") " --- Test read_wave_to_isf_etsf() from ", filename ," ---"
  call read_wave_to_isf_etsf(filename, len(filename), iorbp, hx, hy, hz, n1, n2, n3, nspinor, psiscf)

  nrm = real(0, wp)
  do k = 1, n3, 1
     do j = 1, n2, 1
        do i = 1, n1, 1
           nrm = nrm + psiscf(i, j, k, 1) * psiscf(i, j, k, 1)
           if (nspinor == 2) then
              nrm = nrm + psiscf(i, j, k, 2) * psiscf(i, j, k, 2)
           end if
        end do
     end do
  end do
  write(*,"(A,3F10.6)") " hgrid values for iscf representation:    ", hx, hy, hz
  write(*,"(A,3I10)")   " number of points in iscf representation: ", n1, n2, n3
  write(*,"(A,I2,A,22x,F12.8)")  " norm of orbital ", iorbp, ":", nrm

  call free_wave_to_isf_etsf(psiscf)

  call memocc(0,0,'count','stop')

  call MPI_FINALIZE(ierr)
end program schtroumpf_lib
