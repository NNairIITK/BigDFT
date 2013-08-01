program schtroumpf_lib
  use BigDFT_API

  implicit none

!!$  character(len = *), parameter :: filename = "data/wavefunction.etsf"
  character(len = 1024) :: filename
  integer, parameter :: iorbp = 23
  real(wp), dimension(:,:,:,:), pointer :: psiscf
  real(gp) :: hx, hy, hz
  integer :: n1, n2, n3, nspinor, norbu, norbd, nkpt, iorb, ispin, ikpt, ispinor
  logical :: lstat
  
  real(wp) :: nrm
  integer :: ierr, i, j, k  

  call MPI_INIT(ierr)
  call f_malloc_set_status(memory_limit=0.e0)

  call get_command_argument(1, value = filename)

  write(*,"(3A)") " --- Test read_wave_to_isf_etsf() from ", trim(filename) ," ---"
  call read_wave_descr(lstat, trim(filename), len(trim(filename)), &
       & norbu, norbd, iorb, ispin, nkpt, ikpt, nspinor, ispinor)
  write(*, "(A,2x,4I7)") " ETSF wavefunction file (nou, nod, nk, sp):", &
       norbu, norbd, nkpt, nspinor
  if (.not. lstat) stop

  call read_wave_to_isf(lstat, trim(filename), len(trim(filename)), iorbp, hx, hy, hz, &
       & n1, n2, n3, nspinor, psiscf)
  if (.not. lstat) stop

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

  call free_wave_to_isf(psiscf)

  call memocc(0,0,'count','stop')

  call MPI_FINALIZE(ierr)
end program schtroumpf_lib
