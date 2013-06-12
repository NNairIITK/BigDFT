!> @file
!!  Routines to estimate the use of memory
!! @author
!!    Copyright (C) Luigi Genovese, CEA Grenoble, France, 2007-2011
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 


!>   Estimation of the used memory
subroutine MemoryEstimator(nproc,idsx,lr,nat,norb,nspinor,nkpt,nprojel,nspin,itrpmax,iscf,peakmem)

  use module_base
  use module_types
  use Poisson_Solver
  use yaml_output

  implicit none

  !Arguments
  integer, intent(in) :: nproc,idsx,nat,norb,nspin,nprojel
  integer, intent(in) :: nkpt,nspinor,itrpmax,iscf
  type(locreg_descriptors), intent(in) :: lr
  real(kind=8), intent(out) :: peakmem
  !Local variables
  character(len=*), parameter :: subname='MemoryEstimator'
  real(kind=8), parameter :: eps_mach=1.d-12
  integer :: norbp,nvctrp,n1,n2,n3
  integer :: n01,n02,n03,m1,m2,m3,md1,md2,md3,nd1,nd2,nd3
  integer(kind=8) :: mworkham, mworkrho
  real(kind=8) :: omemwf,omemker,omemden,omempot,omemproj,nden,npotden,npotham,narr
  real(kind=8) :: tt,tmemker,tmemden,tmemps,tmemha
!!$ real(kind=8) :: timinamount

  n1=lr%d%n1
  n2=lr%d%n2
  n3=lr%d%n3

  !here we must add the estimation for the projectors

  tt=dble(norb * nkpt)/dble(nproc)
  norbp=int((1.d0-eps_mach*tt) + tt)
  tt=dble(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f)/dble(nproc)
  nvctrp=int((1.d0-eps_mach*tt) + tt)

  ! Multiply the size of one wavefunction per its spinor value.
!!$  if(nspin==4) then !quadruple size for complex spinors
  norbp=norbp ! do not multiply also the number of bands.
  nvctrp=nvctrp*nspinor
!!$  end if

  !wavefunction memory per orbitals
  omemwf=real(nvctrp*nproc*8,kind=8)
  
  if (lr%geocode == 'P') then
     call P_FFT_dimensions(2*n1+2,2*n2+2,2*n3+2,m1,m2,m3,n01,n02,n03,md1,md2,md3,nd1,nd2,nd3,nproc,.false.)
     n01=2*n1+2
     n02=2*n2+2
     n03=2*n3+2
  else if (lr%geocode == 'S') then
     call S_FFT_dimensions(2*n1+2,2*n2+31,2*n3+2,m1,m2,m3,n01,n02,n03,md1,md2,md3,nd1,nd2,nd3,nproc,0,.false.)
     n01=2*n1+2
     n02=2*n2+31
     n03=2*n3+2
  else if (lr%geocode == 'F') then
     call F_FFT_dimensions(2*n1+31,2*n2+31,2*n3+31,m1,m2,m3,n01,n02,n03,md1,md2,md3,nd1,nd2,nd3,nproc,0,.false.)
     n01=2*n1+31
     n02=2*n2+31
     n03=2*n3+31
  end if
  tt = 8.d0*real(n1*n2*n3,kind=8)/real(n01*n02*n03,kind=8)

  !density memory
  omemden=real(md3*md2/nproc,kind=8)*8.d0*real(md1*nspin,kind=8)
  !kernel memory
  omemker=real(nd2*nd3/nproc,kind=8)*8.d0*real(nd1,kind=8)
  !memory of full grid arrays
  omempot=real(n02*n03,kind=8)*8.d0*real(n01*nspin,kind=8)
  !memory of nonlocal pseudopotential arrays
  omemproj=real(nprojel,kind=8)*8.d0

  ! Work arrays.
  call memspace_work_arrays_sumrho(lr, mworkrho)
  call memspace_work_arrays_locham(lr, mworkham) !n(m)
  ! pot_ion, rhopot, potxc
  nden=3.d0
  ! In Hamiltonian application: pot + psir + work arrays
  npotham=1.d0+nspinor+real(mworkham * 8 * nspinor, kind=8) / omempot
  ! In sumrho: Rho_p + psir + work arrays
  npotden=1.d0+1.d0+real(mworkrho * 8, kind=8) / omempot
  ! Mixing arrays.
  if (itrpmax /= 1) then
     if (mod(iscf, 10) == 1) narr = 5
     if (mod(iscf, 10) == 2) narr = 2
     if (mod(iscf, 10) == 3) narr = 3
     if (mod(iscf, 10) == 4) narr = 5
     if (mod(iscf, 10) == 5) narr = 10
     if (mod(iscf, 10) == 6) narr = 10
     if (mod(iscf, 10) == 7) narr = 1 + 2 * 7
     nden = nden + narr
  end if

  call yaml_comment('Estimation of Memory Consumption',hfill='-')
  call yaml_open_map('Memory requirements for principal quantities (MiB.KiB)')
    call yaml_map('Subspace Matrix',trim(MibdotKib(real(norb,kind=8)**2)),advance='no')
      call yaml_comment('(Number of Orbitals:'//trim(yaml_toa(norb))//')',tabbing=50)
    call yaml_map('Single orbital',trim(MibdotKib(omemwf)),advance='no')
      call yaml_comment('(Number of Components:'//trim(yaml_toa(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f))//')',tabbing=50)
      if(nproc > 1 ) omemwf=24.d0*real(norbp*nvctrp*nproc,kind=8)  !takes into account psit
      if(nproc == 1 ) omemwf=16.d0*real(norbp*nvctrp*nproc,kind=8)
    call yaml_map('All (distributed) orbitals',trim(MibdotKib(omemwf)),advance='no')
      call yaml_comment('(Number of Orbitals per MPI task:'//trim(yaml_toa(norbp))//')',tabbing=50)
      if(nproc > 1 ) omemwf=8.d0*real(2*idsx+3,kind=8)*real(norbp*nvctrp*nproc,kind=8)
      if(nproc == 1 ) omemwf=8.d0*real(2*idsx+2,kind=8)*real(norbp*nvctrp*nproc,kind=8)
    call yaml_map('Wavefunction storage size',trim(MibdotKib(omemwf)),advance='no')
      call yaml_comment('(DIIS/SD workspaces included)',tabbing=50)
    call yaml_map('Nonlocal Pseudopotential Arrays',trim(MibdotKib(omemproj)))
    call yaml_map('Full Uncompressed (ISF) grid',trim(MibdotKib(omempot)))
    call yaml_map('Workspaces storage size',trim(MibdotKib(real(max(mworkrho,mworkham),kind=8))))
  call yaml_close_map()

!!$  write(*,'(1x,a)')&
!!$       '------------------------------------------------------------------ Memory Estimation'
!!$  write(*,'(1x,a,i5,a,i6,a,3(i5))')&
!!$       'Number of atoms=',nat,' Number of orbitals=',norb,' Sim. Box Dimensions= ', n1,n2,n3
!!$  write(*,'(1x,a,i0,a)')&
!!$       'Estimation performed for ',nproc,' processors.'
!!$  write(*,'(1x,a)')&
!!$       'Memory occupation for principal arrays:'
!!$  write(*,'(1x,a,2(i6,a3))')&
!!$       '              Poisson Solver Kernel (K):',mega(omemker),'MB',kappa(omemker),'KB'  
!!$  write(*,'(1x,a,2(i6,a3))')&
!!$       '             Poisson Solver Density (D):',mega(omemden),'MB',kappa(omemden),'KB'  
!!$  write(*,'(1x,a,2(i6,a3))')&
!!$       '    Single Wavefunction for one orbital:',mega(omemwf),'MB',kappa(omemwf),'KB'  
!!$  if(nproc > 1 ) omemwf=24.d0*real(norbp*nvctrp*nproc,kind=8)  !takes into account psit
!!$  if(nproc == 1 ) omemwf=16.d0*real(norbp*nvctrp*nproc,kind=8)
!!$  write(*,'(1x,a,2(i6,a3))')&
!!$       '   All Wavefunctions for each processor:',mega(omemwf),'MB',kappa(omemwf),'KB'  
!!$  if(nproc > 1 ) omemwf=8.d0*real(2*idsx+3,kind=8)*real(norbp*nvctrp*nproc,kind=8)
!!$  if(nproc == 1 ) omemwf=8.d0*real(2*idsx+2,kind=8)*real(norbp*nvctrp*nproc,kind=8)
!!$  write(*,'(1x,a,2(i6,a3))')&
!!$       '      Wavefunctions + DIIS per proc (W):',mega(omemwf),'MB',kappa(omemwf),'KB'  
!!$  write(*,'(1x,a,2(i6,a3))')&
!!$       '    Nonlocal Pseudopotential Arrays (P):',mega(omemproj),'MB',kappa(omemproj),'KB'  
!!$  write(*,'(1x,a,2(i6,a3))')&
!!$       '   Arrays of full uncompressed grid (U):',mega(omempot),'MB',kappa(omempot),'KB' 
!!$
!!$  write(*,'(1x,a)')&
!!$       'Estimation of Memory requirements for principal code sections:'
!!$  write(*,'(1x,a)')&
!!$       ' Kernel calculation | Density Construction | Poisson Solver | Hamiltonian application'
  if (nproc > 1) then 
!!$     write(*,'(1x,a,I0,a,I2,a,I0,a,I2,a)')&
!!$       '      ~19*K         |   W+~',nint(npotden),'*U+~',nint(nden),&
!!$       & '*D+K+P   |   ~12*D+K+W+P  |   W+~',nint(npotham),'*U+~',nint(nden),'*D+K+P '
     tmemker=19.d0*omemker
     tmemden=omemwf+nden*omemden+npotden*omempot+omemker+omemproj
     tmemps=12.d0*omemden+omemwf+omemker+omemproj
     tmemha=nden*omemden+npotham*omempot+omemwf+omemker+omemproj
  else
!!$     write(*,'(1x,a,I0,a,I2,a,I0,a,I2,a)')&
!!$       '      ~11*K         |   W+~',nint(npotden - 1.d0),'*U+~',nint(nden),&
!!$       & '*D+K+P   |   ~8*D+K+W+P   |   W+~',nint(npotham - 1.d0),'*U+~',nint(nden),'*D+K+P '
     tmemker=11.d0*omemker
     tmemden=omemwf+nden*omemden+(npotden-1.d0)*omempot+omemker+omemproj
     tmemps=8.d0*omemden+omemwf+omemker+omemproj
     tmemha=nden*omemden+(npotham-1.d0)*omempot+omemwf+omemker+omemproj
  end if
!!$  write(*,'(1x,4(1x,i8,a))')&
!!$       mega(tmemker),'MB         | ',mega(tmemden),'MB          |',&
!!$       mega(tmemps),'MB     |     ',mega(tmemha),'MB'
  !estimation of the memory peak
  peakmem=max(tmemker,tmemden,tmemps,tmemha)

  call yaml_open_map('Accumulated memory requirements during principal run stages (MiB.KiB)')
     call yaml_map('Kernel calculation',trim(MibdotKib(tmemker)))
     call yaml_map('Density Construction',trim(MibdotKib(tmemden)))
     call yaml_map('Poisson Solver',trim(MibdotKib(tmemps)))
     call yaml_map('Hamiltonian application',trim(MibdotKib(tmemha)))
!           call yaml_comment('Wfn, Work, Den, Ker ',tabbing=50)
  call yaml_close_map()
  call yaml_map('Estimated Memory Peak (MB)',yaml_toa(mega(peakmem)))


!!$  write(*,'(1x,a,i0,a)')&
!!$       'The overall memory requirement needed for this calculation is thus: ',&
!!$       mega(peakmem),' MB'
!!$  tminamount=real(3*(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f)*8,kind=8)+3.d0*real(n01*n02,kind=8)+&
!!$       (3.d0+2.d0*tt)*omempot+omemproj
!!$  write(*,'(1x,a)')&
!!$       'By reducing the DIIS history and/or increasing the number of processors the amount of'
!!$  write(*,'(1x,a,i0,a)')&
!!$       ' memory can be reduced but for this system it will never be less than ',&
!!$       mega(tminamount),' MB'

contains

  function mega(omemory)
    implicit none
    real(kind=8), intent(in) :: omemory
    integer(kind=8) :: mega
    mega=int(omemory/1048576.d0,kind=8)
  end function mega

  function kappa(omemory)
    implicit none
    real(kind=8), intent(in) :: omemory
    integer :: kappa
    kappa=ceiling((omemory-aint(omemory/1048576.d0)*1048576.d0)/1024.d0)
  end function kappa

  function MiBdotKiB(omemory)
    implicit none
    real(kind=8), intent(in) :: omemory
    character(len=50) MiBdotKiB

    MiBdotKiB=repeat(' ',len(MiBdotKiB))

    MiBdotKiB=trim(adjustl(yaml_toa(int(mega(omemory)))))//'.'//&
         trim(adjustl(yaml_toa(int(kappa(omemory)))))
    
  end function MiBdotKiB

END SUBROUTINE MemoryEstimator
