!!subroutine check_cutoff(iproc, nproc, orbs, lzd, hx, hy, hz, locrad, confdatarr, psi)
!!use module_base
!!use module_types
!!use module_interfaces, except_this_one => apply_position_operators
!!implicit none
!!
!!! Calling arguments
!!integer,intent(in):: iproc, nproc
!!type(orbitals_data),intent(in):: orbs
!!type(local_zone_descriptors),intent(in):: lzd
!!real(8),intent(in):: hx, hy, hz, locrad
!!type(confpot_data),dimension(orbs%norbp),intent(in):: confdatarr
!!real(8),dimension(max(orbs%npsidim_orbs,orbs%npsidim_comp)),intent(in):: psi
!!
!!! Local variables
!!integer:: oidx, iorb, ilr, npot, icenter, i_stat, i_all, ist_c, ist_f, ist, iiorb, iall, ierr
!!real(8):: hxh, hyh, hzh, ddot, tt, t1, t2, time, weight_in, weight_out
!!real(8),dimension(:,:),allocatable:: psir, psirx, psiry, psirz
!!type(workarr_sumrho):: work_sr
!!real(8),dimension(0:3),parameter:: scal=1.d0
!!real(8),dimension(:,:,:),allocatable:: ypsitemp_c
!!real(8),dimension(:,:,:,:),allocatable:: ypsitemp_f
!!character(len=*),parameter:: subname='apply_orbitaldependent_potential'
!!integer, dimension(3) :: ishift !temporary variable in view of wavefunction creation
!!
!!  ishift=(/0,0,0/)
!!
!!  oidx = 0
!!  do iorb=1,orbs%norbp
!!     ilr = orbs%inwhichlocreg(iorb+orbs%isorb)
!!
!!  
!!     !initialise the work arrays
!!     call initialize_work_arrays_sumrho(lzd%llr(ilr), work_sr)
!!
!!     ! Wavefunction in real space
!!     allocate(psir(Lzd%Llr(ilr)%d%n1i*Lzd%Llr(ilr)%d%n2i*Lzd%Llr(ilr)%d%n3i,orbs%nspinor+ndebug),stat=i_stat)
!!     call memocc(i_stat,psir,'psir',subname)
!!     call razero(Lzd%Llr(ilr)%d%n1i*Lzd%Llr(ilr)%d%n2i*Lzd%Llr(ilr)%d%n3i*orbs%nspinor,psir)
!!
!!     call daub_to_isf(lzd%llr(ilr), work_sr, psi(1+oidx), psir)
!!     !apply the potential to the psir wavefunction and calculate potential energy
!!     !icenter=confinementCenter(iorb)
!!     !components of the potential
!!     npot=orbs%nspinor
!!     if (orbs%nspinor == 2) npot=1
!!
!!     write(*,'(a,2i8,3es16.7)') 'iproc, iorb, confdatarr(iorb)%rxyzConf', iproc, iorb, confdatarr(iorb)%rxyzConf
!!
!!     call get_cutoff_weight(lzd%llr(ilr)%d%n1i, lzd%llr(ilr)%d%n2i, lzd%llr(ilr)%d%n3i, &
!!                             lzd%llr(ilr)%d%n1i, lzd%llr(ilr)%d%n2i, lzd%llr(ilr)%d%n3i, &
!!                             ishift, lzd%llr(ilr)%d%n2, lzd%llr(ilr)%d%n3, orbs%nspinor, &
!!                             psir, locrad, weight_in, weight_out, &
!!                             confdatarr(iorb), lzd%llr(ilr)%bounds%ibyyzz_r) !optional
!!
!!     write(*,'(a,2i8,3es16.6)') 'iproc, iorb, weight_in, weight_out, ratio', &
!!         iproc, iorb, weight_in, weight_out, weight_in/(weight_in+weight_out)
!!
!!     i_all=-product(shape(psir))*kind(psir)
!!     deallocate(psir,stat=i_stat)
!!     call memocc(i_stat,i_all,'psir',subname)
!!
!!     call deallocate_work_arrays_sumrho(work_sr)
!!
!!     oidx = oidx + (Lzd%Llr(ilr)%wfd%nvctr_c+7*Lzd%Llr(ilr)%wfd%nvctr_f)*orbs%nspinor
!!
!!  enddo
!!
!!
!!end subroutine check_cutoff
!!
!!
!!
!!subroutine get_cutoff_weight(n1i,n2i,n3i,n1ip,n2ip,n3ip,ishift,n2,n3,nspinor,psir,&
!!     cutoff, weight_in, weight_out, &
!!     confdata,ibyyzz_r) !optional
!!  use module_base
!!  use module_types
!!  implicit none
!!  integer, intent(in) :: n1i,n2i,n3i,n1ip,n2ip,n3ip,n2,n3,nspinor
!!  integer, dimension(3), intent(in) :: ishift !<offset of potential box in wfn box coords.
!!  real(wp), dimension(n1i,n2i,n3i,nspinor), intent(in) :: psir !< real-space wfn in lr
!!  real(8),intent(in):: cutoff
!!  real(8),intent(out):: weight_in, weight_out
!!  type(confpot_data), intent(in), optional :: confdata !< data for the confining potential
!!  integer, dimension(2,-14:2*n2+16,-14:2*n3+16), intent(in), optional :: ibyyzz_r !< bounds in lr
!!  !local variables
!!  integer :: i1,i2,i3,ispinor,i1s,i1e,i2s,i2e,i3s,i3e,i1st,i1et
!!  real(wp) :: tt11,tt22,tt33,tt44,tt13,tt14,tt23,tt24,tt31,tt32,tt41,tt42,tt
!!  real(wp) :: psir1,psir2,psir3,psir4,pot1,pot2,pot3,pot4
!!  real(wp):: ttx, tty, ttz, cutoff2
!!
!!  ! Square of the cutoff radius
!!  cutoff2=cutoff**2
!!
!!  ! Initialize return values
!!  weight_in=0.d0
!!  weight_out=0.d0
!!
!!
!!  !loop on wavefunction
!!  !calculate the limits in all the directions
!!  !regions in which both the potential and wavefunctions are defined
!!  i3s=max(1,ishift(3)+1)
!!  i3e=min(n3i,n3ip+ishift(3))
!!  i2s=max(1,ishift(2)+1)
!!  i2e=min(n2i,n2ip+ishift(2))
!!  i1s=max(1,ishift(1)+1)
!!  i1e=min(n1i,n1ip+ishift(1))
!!
!!
!!  !$omp parallel default(none)&
!!  !$omp shared(psir,n1i,n2i,n3i,n1ip,n2ip,n3ip,n2,n3,ibyyzz_r,nspinor)&
!!  !$omp shared(i1s,i1e,i2s,i2e,i3s,i3e,ishift,confdata,weight_in,weight_out)&
!!  !$omp private(ispinor,i1,i2,i3,i1st,i1et)&
!!  !$omp private(tt11,tt22,tt33,tt44,tt13,tt14,tt23,tt24,tt31,tt32,tt41,tt42,tt)&
!!  !$omp private(psir1,psir2,psir3,psir4,pot1,pot2,pot3,pot4,ttx,tty,ttz,cutoff2)
!!
!!!!$  !$omp parallel default(private)&
!!!!$  !$omp shared(pot,psir,n1i,n2i,n3i,n1ip,n2ip,n3ip,n2,n3,ibyyzz_r,nspinor)&
!!!!$  !$omp shared(i1s,i1e,i2s,i2e,i3s,i3e,ishift)
!!  !case without bounds
!!
!!
!!  !!!put to zero the external part of psir if the potential is more little than the wavefunction
!!  !!!first part of the array
!!  !!do ispinor=1,nspinor
!!  !!   !$omp do 
!!  !!   do i3=1,i3s-1
!!  !!      do i2=1,n2i
!!  !!         do i1=1,n1i
!!  !!           psirx(i1,i2,i3,ispinor)=0.0_wp 
!!  !!           psiry(i1,i2,i3,ispinor)=0.0_wp 
!!  !!           psirz(i1,i2,i3,ispinor)=0.0_wp 
!!  !!         end do
!!  !!      end do
!!  !!   end do
!!  !!   !$omp end do
!!  !!end do
!!
!!  !!!central part of the array
!!  !!do ispinor=1,nspinor
!!  !!   !$omp do 
!!  !!   do i3=i3s,i3e
!!
!!  !!      !first part
!!  !!      do i2=1,i2s-1
!!  !!         do i1=1,n1i
!!  !!            psirx(i1,i2,i3,ispinor)=0.0_wp 
!!  !!            psiry(i1,i2,i3,ispinor)=0.0_wp 
!!  !!            psirz(i1,i2,i3,ispinor)=0.0_wp 
!!  !!         end do
!!  !!      end do
!!  !!      !central part
!!  !!      do i2=i2s,i2e
!!  !!         do i1=1,i1s-1
!!  !!            psirx(i1,i2,i3,ispinor)=0.0_wp 
!!  !!            psiry(i1,i2,i3,ispinor)=0.0_wp 
!!  !!            psirz(i1,i2,i3,ispinor)=0.0_wp 
!!  !!         end do
!!  !!         do i1=i1e+1,n1i
!!  !!            psirx(i1,i2,i3,ispinor)=0.0_wp 
!!  !!            psiry(i1,i2,i3,ispinor)=0.0_wp 
!!  !!            psirz(i1,i2,i3,ispinor)=0.0_wp 
!!  !!         end do
!!  !!      end do
!!  !!      !last part
!!  !!      do i2=i2e+1,n2i
!!  !!         do i1=1,n1i
!!  !!            psirx(i1,i2,i3,ispinor)=0.0_wp 
!!  !!            psiry(i1,i2,i3,ispinor)=0.0_wp 
!!  !!            psirz(i1,i2,i3,ispinor)=0.0_wp 
!!  !!         end do
!!  !!      end do
!!
!!  !!   end do
!!  !!   !$omp end do
!!  !!end do
!!
!!
!!  !!!last part of the array
!!  !!do ispinor=1,nspinor
!!  !!   !$omp do 
!!  !!   do i3=i3e+1,n3i
!!  !!      do i2=1,n2i
!!  !!         do i1=1,n1i
!!  !!            psirx(i1,i2,i3,ispinor)=0.0_wp 
!!  !!            psiry(i1,i2,i3,ispinor)=0.0_wp 
!!  !!            psirz(i1,i2,i3,ispinor)=0.0_wp 
!!  !!         end do
!!  !!      end do
!!  !!   end do
!!  !!   !$omp end do
!!  !!end do
!!
!!
!!  !important part of the array
!!  if (nspinor==4) then
!!      stop 'not yet implemented for nspinor==4!'
!!     !!!$omp do
!!     !!do i3=i3s,i3e
!!     !!   do i2=i2s,i2e
!!     !!      !thanks to the optional argument the conditional is done at compile time
!!     !!      if (present(ibyyzz_r)) then
!!     !!         i1st=max(i1s,ibyyzz_r(1,i2-15,i3-15)+1) !in bounds coordinates
!!     !!         i1et=min(i1e,ibyyzz_r(2,i2-15,i3-15)+1) !in bounds coordinates
!!     !!      else
!!     !!         i1st=i1s
!!     !!         i1et=i1e
!!     !!      end if
!!     !!      !no need of setting up to zero values outside wavefunction bounds
!!     !!      do i1=i1st,i1et
!!     !!         !wavefunctions
!!     !!         psir1=psir(i1,i2,i3,1)
!!     !!         psir2=psir(i1,i2,i3,2)
!!     !!         psir3=psir(i1,i2,i3,3)
!!     !!         psir4=psir(i1,i2,i3,4)
!!     !!         !potentials + confining term
!!     !!         pot1=pot(i1-ishift(1),i2-ishift(2),i3-ishift(3),1)+cp(i1,i2,i3)
!!     !!         pot2=pot(i1-ishift(1),i2-ishift(2),i3-ishift(3),2)+cp(i1,i2,i3)
!!     !!         pot3=pot(i1-ishift(1),i2-ishift(2),i3-ishift(3),3)+cp(i1,i2,i3)
!!     !!         pot4=pot(i1-ishift(1),i2-ishift(2),i3-ishift(3),4)+cp(i1,i2,i3)
!!
!!     !!         !diagonal terms
!!     !!         tt11=pot1*psir1 !p1
!!     !!         tt22=pot1*psir2 !p2
!!     !!         tt33=pot4*psir3 !p3
!!     !!         tt44=pot4*psir4 !p4
!!     !!         !Rab*Rb
!!     !!         tt13=pot2*psir3 !p1
!!     !!         !Iab*Ib
!!     !!         tt14=pot3*psir4 !p1
!!     !!         !Rab*Ib
!!     !!         tt23=pot2*psir4 !p2
!!     !!         !Iab*Rb
!!     !!         tt24=pot3*psir3 !p2
!!     !!         !Rab*Ra
!!     !!         tt31=pot2*psir1 !p3
!!     !!         !Iab*Ia
!!     !!         tt32=pot3*psir2 !p3
!!     !!         !Rab*Ia
!!     !!         tt41=pot2*psir2 !p4
!!     !!         !Iab*Ra
!!     !!         tt42=pot3*psir1 !p4
!!
!!     !!         !value of the potential energy
!!     !!         epot_p=epot_p+tt11*psir1+tt22*psir2+tt33*psir3+tt44*psir4+&
!!     !!              2.0_gp*tt31*psir3-2.0_gp*tt42*psir4+2.0_gp*tt41*psir4+2.0_gp*tt32*psir3
!!
!!     !!         !wavefunction update
!!     !!         !p1=h1p1+h2p3-h3p4
!!     !!         !p2=h1p2+h2p4+h3p3
!!     !!         !p3=h2p1+h3p2+h4p3
!!     !!         !p4=h2p2-h3p1+h4p4
!!     !!         psir(i1,i2,i3,1)=tt11+tt13-tt14
!!     !!         psir(i1,i2,i3,2)=tt22+tt23+tt24
!!     !!         psir(i1,i2,i3,3)=tt33+tt31+tt32
!!     !!         psir(i1,i2,i3,4)=tt44+tt41-tt42
!!     !!      end do
!!     !!   end do
!!     !!end do
!!     !!!$omp end do
!!
!!  else !case with nspinor /=4
!!     do ispinor=1,nspinor
!!        !$omp do
!!        do i3=i3s,i3e
!!           do i2=i2s,i2e
!!              !thanks to the optional argument the conditional is done at compile time
!!              if (present(ibyyzz_r)) then
!!                 i1st=max(i1s,ibyyzz_r(1,i2-15,i3-15)+1) !in bounds coordinates
!!                 i1et=min(i1e,ibyyzz_r(2,i2-15,i3-15)+1) !in bounds coordinates
!!              else
!!                 i1st=i1s
!!                 i1et=i1e
!!              end if
!!              !no need of setting up to zero values outside wavefunction bounds
!!              do i1=i1st,i1et
!!                 psir1=psir(i1,i2,i3,ispinor)
!!                 !the local potential is always real (npot=1) + confining term
!!                 !!pot1=pot(i1-ishift(1),i2-ishift(2),i3-ishift(3),1)+cp(i1,i2,i3)
!!                 ttx=(confdata%hh(1)*real(i1+confdata%ioffset(1),wp)-confdata%rxyzConf(1))**2
!!                 tty=(confdata%hh(2)*real(i2+confdata%ioffset(2),wp)-confdata%rxyzConf(2))**2
!!                 ttz=(confdata%hh(3)*real(i3+confdata%ioffset(3),wp)-confdata%rxyzConf(3))**2
!!                 tt=ttx+tty+ttz
!!                 !write(1000,*) tt, cutoff2, psir1**2
!!                 if(tt>cutoff2) then
!!                     weight_out=weight_out+psir1**2
!!                 else
!!                     weight_in=weight_in+psir1**2
!!                 end if
!!              end do
!!           end do
!!        end do
!!        !$omp end do
!!     end do
!!  end if
!!  
!!  
!!  !$omp end parallel
!!
!!
!!END SUBROUTINE get_cutoff_weight



!!!subroutine diagonalizeHamiltonian(iproc, nproc, orbs, HamSmall, eval)
!!!!
!!!! Purpose:
!!!! ========
!!!!   Diagonalizes the Hamiltonian HamSmall and makes sure that all MPI processes give
!!!!   the same result. This is done by requiring that the first entry of each vector
!!!!   is positive.
!!!!
!!!! Calling arguments:
!!!! ==================
!!!!   Input arguments:
!!!!   ----------------
!!!!     iproc     process ID
!!!!     nproc     number of MPI processes
!!!!     orbs      type describing the physical orbitals psi
!!!!   Input / Putput arguments
!!!!     HamSmall  on input: the Hamiltonian
!!!!               on exit: the eigenvectors
!!!!   Output arguments
!!!!     eval      the associated eigenvalues 
!!!!
!!!use module_base
!!!use module_types
!!!implicit none
!!!
!!!! Calling arguments
!!!integer:: iproc, nproc
!!!type(orbitals_data), intent(inout) :: orbs
!!!real(8),dimension(orbs%norb, orbs%norb):: HamSmall
!!!real(8),dimension(orbs%norb):: eval
!!!
!!!! Local variables
!!!integer:: lwork, info, istat, iall, i, iorb, jorb
!!!real(8),dimension(:),allocatable:: work
!!!character(len=*),parameter:: subname='diagonalizeHamiltonian'
!!!
!!!  ! Get the optimal work array size
!!!  lwork=-1 
!!!  allocate(work(1), stat=istat)
!!!  call memocc(istat, work, 'work', subname)
!!!  call dsyev('v', 'l', orbs%norb, HamSmall(1,1), orbs%norb, eval(1), work(1), lwork, info) 
!!!  lwork=work(1) 
!!!
!!!  ! Deallocate the work array ane reallocate it with the optimal size
!!!  iall=-product(shape(work))*kind(work)
!!!  deallocate(work, stat=istat) ; if(istat/=0) stop 'ERROR in deallocating work' 
!!!  call memocc(istat, iall, 'work', subname)
!!!  allocate(work(lwork), stat=istat) ; if(istat/=0) stop 'ERROR in allocating work' 
!!!  call memocc(istat, work, 'work', subname)
!!!
!!!  ! Diagonalize the Hamiltonian
!!!  call dsyev('v', 'l', orbs%norb, HamSmall(1,1), orbs%norb, eval(1), work(1), lwork, info) 
!!!
!!!  ! Deallocate the work array.
!!!  iall=-product(shape(work))*kind(work)
!!!  deallocate(work, stat=istat) ; if(istat/=0) stop 'ERROR in deallocating work' 
!!!  call memocc(istat, iall, 'work', subname)
!!!  
!!!  ! Make sure that the eigenvectors are the same for all MPI processes. To do so, require that 
!!!  ! the first entry of each vector is positive.
!!!  do iorb=1,orbs%norb
!!!      if(HamSmall(1,iorb)<0.d0) then
!!!          do jorb=1,orbs%norb
!!!              HamSmall(jorb,iorb)=-HamSmall(jorb,iorb)
!!!          end do
!!!      end if
!!!  end do
!!!
!!!
!!!end subroutine diagonalizeHamiltonian



!!!subroutine get_potential_matrices(iproc, nproc, at, orbs, lzd, op, comon, mad, rxyz, &
!!!           confdatarr, hx, psi, potmat)
!!!use module_base
!!!use module_types
!!!use module_interfaces, eccept_this_one => get_potential_matrices
!!!implicit none
!!!
!!!! Calling arguments
!!!integer,intent(in):: iproc, nproc
!!!type(atoms_data),intent(in):: at
!!!type(orbitals_data),intent(in):: orbs
!!!type(local_zone_descriptors),intent(in):: lzd
!!!type(overlapParameters),intent(inout):: op
!!!type(p2pComms),intent(inout):: comon
!!!type(matrixDescriptors),intent(in):: mad
!!!real(8),dimension(3,at%nat),intent(in):: rxyz
!!!real(8),intent(in):: hx
!!!type(confpot_data),dimension(orbs%norbp),intent(in):: confdatarr
!!!real(8),dimension(max(orbs%npsidim_orbs,orbs%npsidim_comp)),intent(inout):: psi
!!!real(8),dimension(orbs%norb,orbs%norb,at%nat),intent(out):: potmat
!!!
!!!! Local variables
!!!integer:: iorb, ilr, ilrold, istat, iall
!!!real(8),dimension(:,:),allocatable:: ttmat
!!!real(8):: tt1, tt2, tt3, tt4, tt5
!!!real(8),dimension(:),allocatable:: vpsi
!!!character(len=*),parameter:: subname='get_potential_matrices'
!!!
!!!allocate(vpsi(max(orbs%npsidim_orbs,orbs%npsidim_comp)), stat=istat)
!!!call memocc(istat, vpsi, 'vpsi', subname)
!!!
!!!
!!!ilrold=-1
!!!do iorb=1,orbs%norb
!!!    ilr=orbs%inwhichlocreg(iorb)
!!!    if(ilr==ilrold) cycle
!!!    call apply_orbitaldependent_potential(iproc, nproc, at, orbs, lzd, rxyz, &
!!!         confdatarr, hx, psi, ilr, vpsi)
!!!
!!!    !call extractOrbital3(iproc, nproc, orbs, orbs%npsidim, orbs%inWhichLocreg, lzd, op, vpsi, comon%nsendBuf, comon%sendBuf)
!!!    !call postCommsOverlapNew(iproc, nproc, orbs, op, lzd, vpsi, comon, tt1, tt2)
!!!    !allocate(ttmat(lin%orbs%norb,lin%orbs%norb))
!!!    !call collectnew(iproc, nproc, comon, lin%mad,lin%op, lin%orbs, input, lin%lzd, comon%nsendbuf, &
!!!    !     comon%sendbuf, comon%nrecvbuf, comon%recvbuf, ttmat, tt3, tt4, tt5)
!!!    !deallocate(ttmat)
!!!    call getMatrixElements2(iproc, nproc, lzd, orbs, op, comon, psi, vpsi, mad, potmat(1,1,ilr))
!!!    ilrold=ilr
!!!    
!!!end do
!!!
!!!iall=-product(shape(vpsi))*kind(vpsi)
!!!deallocate(vpsi, stat=istat)
!!!call memocc(istat, iall, 'vpsi', subname)
!!!
!!!
!!!
!!!end subroutine get_potential_matrices


!!!subroutine buildWavefunction(iproc, nproc, orbs, orbsLIN, comms, commsLIN, phi, psi, HamSmall)
!!!!
!!!! Purpose:
!!!! =======
!!!!   Builds the physical orbitals psi as a linear combination of the basis functions phi. The coefficients
!!!!   for this linear combination are obtained by diagonalizing the Hamiltonian matrix HamSmall.
!!!!
!!!! Calling arguments:
!!!! ==================
!!!!   Input arguments:
!!!!   ----------------
!!!!     iproc      process ID
!!!!     nproc      total number of processes
!!!!     orbs       type describing the physical orbitals psi
!!!!     orbsLIN    type describing the basis functions phi
!!!!     comms      type containing the communication parameters for the physical orbitals psi
!!!!     commsLIN   type containing the communication parameters for the basis functions phi
!!!!     phi        the basis functions 
!!!!     HamSmall   the  Hamiltonian matrix
!!!!   Output arguments:
!!!!   -----------------
!!!!     psi        the physical orbitals 
!!!!
!!!
!!!use module_base
!!!use module_types
!!!implicit none
!!!
!!!! Calling arguments
!!!integer:: iproc, nproc
!!!type(orbitals_data), intent(in) :: orbs
!!!type(orbitals_data), intent(in) :: orbsLIN
!!!type(comms_cubic), intent(in) :: comms
!!!type(comms_cubic), intent(in) :: commsLIN
!!!real(8),dimension(sum(commsLIN%nvctr_par(iproc,1:orbsLIN%nkptsp))*orbsLIN%nspinor,orbsLIN%norb) :: phi
!!!real(8),dimension(sum(comms%nvctr_par(iproc,1:orbs%nkptsp))*orbs%nspinor,orbs%norb) :: psi
!!!real(8),dimension(orbsLIN%norb,orbsLIN%norb):: HamSmall
!!!
!!!! Local variables
!!!integer:: nvctrp
!!!
!!!
!!!  nvctrp=sum(comms%nvctr_par(iproc,1:orbs%nkptsp))*orbs%nspinor
!!!  call dgemm('n', 'n', nvctrp, orbs%norb, orbsLIN%norb, 1.d0, phi(1,1), nvctrp, HamSmall(1,1), &
!!!             orbsLIN%norb, 0.d0, psi(1,1), nvctrp)
!!!  
!!!
!!!end subroutine buildWavefunction
!!
!!
!!

!!
!!!subroutine buildWavefunctionModified(iproc, nproc, orbs, orbsLIN, comms, commsLIN, phi, psi, coeff)
!!!
!!!!
!!!! Purpose:
!!!! =======
!!!!   Builds the physical orbitals psi as a linear combination of the basis functions phi. The coefficients
!!!!   for this linear combination are obtained by diagonalizing the Hamiltonian matrix HamSmall.
!!!!
!!!! Calling arguments:
!!!! ==================
!!!!   Input arguments:
!!!!   ----------------
!!!!     iproc      process ID
!!!!     nproc      total number of processes
!!!!     orbs       type describing the physical orbitals psi
!!!!     orbsLIN    type describing the basis functions phi
!!!!     comms      type containing the communication parameters for the physical orbitals psi
!!!!     commsLIN   type containing the communication parameters for the basis functions phi
!!!!     phi        the basis functions 
!!!!     coeff      the coefficients for the linear combination
!!!!   Output arguments:
!!!!   -----------------
!!!!     psi        the physical orbitals 
!!!!
!!!
!!!use module_base
!!!use module_types
!!!implicit none
!!!
!!!! Calling arguments
!!!integer:: iproc, nproc
!!!type(orbitals_data), intent(in) :: orbs
!!!type(orbitals_data), intent(in) :: orbsLIN
!!!type(comms_cubic), intent(in) :: comms
!!!type(comms_cubic), intent(in) :: commsLIN
!!!real(8),dimension(sum(commsLIN%nvctr_par(iproc,1:orbsLIN%nkptsp))*orbsLIN%nspinor,orbsLIN%norb) :: phi
!!!real(8),dimension(sum(comms%nvctr_par(iproc,1:orbs%nkptsp))*orbs%nspinor,orbs%norb) :: psi
!!!real(8),dimension(orbsLIN%norb,orbs%norb):: coeff
!!!
!!!! Local variables
!!!integer:: nvctrp
!!!
!!!
!!!  nvctrp=sum(comms%nvctr_par(iproc,1:orbs%nkptsp))*orbs%nspinor
!!!  call dgemm('n', 'n', nvctrp, orbs%norb, orbsLIN%norb, 1.d0, phi(1,1), nvctrp, coeff(1,1), &
!!!             orbsLIN%norb, 0.d0, psi(1,1), nvctrp)
!!!  
!!!
!!!end subroutine buildWavefunctionModified


!!subroutine cut_at_boundaries(lzd, orbs, psi)
!!  use module_base
!!  use module_types
!!  use module_interfaces
!!  implicit none
!!
!!  ! Calling arguments
!!  type(local_zone_descriptors),intent(in):: lzd
!!  type(orbitals_data),intent(in):: orbs
!!  real(8),dimension(orbs%npsidim_orbs),intent(inout):: psi
!!
!!  ! Local variables
!!  integer:: istc, istf, iorb, iiorb, ilr, i1, i2, i3, istat, iall, ii1, ii2, ii3
!!  real(8):: r0, r1, r2, r3, rr, tt
!!  real(8),dimension(:,:,:,:,:,:),allocatable:: psig
!!  character(len=*),parameter:: subname='flatten_at_boundaries'
!!  
!!
!!  istc=1
!!  istf=1
!!  do iorb=1,orbs%norbp
!!      iiorb=orbs%isorb+iorb
!!      ilr=orbs%inwhichlocreg(iiorb)
!!
!!      allocate(psig(0:lzd%llr(ilr)%d%n1,2,0:lzd%llr(ilr)%d%n2,2,0:lzd%llr(ilr)%d%n3,2), stat=istat)
!!      call memocc(istat, psig, 'psig', subname)
!!      call to_zero(8*(lzd%llr(ilr)%d%n1+1)*(lzd%llr(ilr)%d%n2+1)*(lzd%llr(ilr)%d%n3+1), psig(0,1,0,1,0,1))
!!
!!      istf = istf + lzd%llr(ilr)%wfd%nvctr_c
!!      call uncompress(lzd%llr(ilr)%d%n1, lzd%llr(ilr)%d%n2, lzd%llr(ilr)%d%n3, &
!!           lzd%llr(ilr)%wfd%nseg_c, lzd%llr(ilr)%wfd%nvctr_c, lzd%llr(ilr)%wfd%keygloc, lzd%llr(ilr)%wfd%keyvloc,  &
!!           lzd%llr(ilr)%wfd%nseg_f, lzd%llr(ilr)%wfd%nvctr_f, &
!!           lzd%llr(ilr)%wfd%keygloc(1,lzd%llr(ilr)%wfd%nseg_c+min(1,lzd%llr(ilr)%wfd%nseg_f)), &
!!           lzd%llr(ilr)%wfd%keyvloc(lzd%llr(ilr)%wfd%nseg_c+min(1,lzd%llr(ilr)%wfd%nseg_f)), &
!!           psi(istc), psi(istf), psig)
!!
!!      r0 = (lzd%llr(ilr)%locrad-8.d0*lzd%hgrids(1))**2
!!      !r0 = (lzd%llr(ilr)%locrad-2.d0)**2
!!      do i3=0,lzd%llr(ilr)%d%n3
!!          ii3=lzd%llr(ilr)%ns3+i3
!!          r3 = (dble(ii3)*lzd%hgrids(3)-lzd%llr(ilr)%locregCenter(3))**2
!!          do i2=0,lzd%llr(ilr)%d%n2
!!              ii2=lzd%llr(ilr)%ns2+i2
!!              r2 = (dble(ii2)*lzd%hgrids(2)-lzd%llr(ilr)%locregCenter(2))**2
!!              do i1=0,lzd%llr(ilr)%d%n1
!!                  ii1=lzd%llr(ilr)%ns1+i1
!!                  r1 = (dble(ii1)*lzd%hgrids(1)-lzd%llr(ilr)%locregCenter(1))**2
!!                  rr=r1+r2+r3
!!                  !write(999,'(5es14.3)') r1, r2, r3, rr, r0
!!                  if(rr>=r0) then
!!                      psig(i1,1,i2,1,i3,1)=0.d0
!!                      psig(i1,2,i2,1,i3,1)=0.d0
!!                      psig(i1,1,i2,2,i3,1)=0.d0
!!                      psig(i1,2,i2,2,i3,1)=0.d0
!!                      psig(i1,1,i2,1,i3,2)=0.d0
!!                      psig(i1,2,i2,1,i3,2)=0.d0
!!                      psig(i1,1,i2,2,i3,2)=0.d0
!!                      psig(i1,2,i2,2,i3,2)=0.d0
!!                  else
!!                      !write(*,*) 'not zero'
!!                  end if
!!              end do
!!          end do
!!      end do
!!
!!      call compress(lzd%llr(ilr)%d%n1, lzd%llr(ilr)%d%n2, &
!!           0, lzd%llr(ilr)%d%n1, 0, lzd%llr(ilr)%d%n2, 0, lzd%llr(ilr)%d%n3, &
!!           lzd%llr(ilr)%wfd%nseg_c, lzd%llr(ilr)%wfd%nvctr_c, lzd%llr(ilr)%wfd%keygloc, lzd%llr(ilr)%wfd%keyvloc,  &
!!           lzd%llr(ilr)%wfd%nseg_f, lzd%llr(ilr)%wfd%nvctr_f, &
!!           lzd%llr(ilr)%wfd%keygloc(1,lzd%llr(ilr)%wfd%nseg_c+min(1,lzd%llr(ilr)%wfd%nseg_f)), &
!!           lzd%llr(ilr)%wfd%keyvloc(lzd%llr(ilr)%wfd%nseg_c+min(1,lzd%llr(ilr)%wfd%nseg_f)),  &
!!           psig, psi(istc), psi(istf))
!!
!!      istf = istf + 7*lzd%llr(ilr)%wfd%nvctr_f
!!      istc = istc + lzd%llr(ilr)%wfd%nvctr_c + 7*lzd%llr(ilr)%wfd%nvctr_f
!!
!!      iall=-product(shape(psig))*kind(psig)
!!      deallocate(psig,stat=istat)
!!      call memocc(istat,iall,'psig',subname)
!!
!!
!!  end do
!!
!!end subroutine cut_at_boundaries






!!subroutine cut_at_boundaries2(lr, orbs, hx, hy, hz, psi)
!!  use module_base
!!  use module_types
!!  use module_interfaces
!!  implicit none
!!
!!  ! Calling arguments
!!  type(locreg_descriptors),intent(in):: lr
!!  type(orbitals_data),intent(in):: orbs
!!  real(8),intent(in):: hx, hy, hz
!!  real(8),dimension(orbs%npsidim_orbs),intent(inout):: psi
!!
!!  ! Local variables
!!  integer:: istc, istf, iorb, i1, i2, i3, istat, iall, ii1, ii2, ii3
!!  real(8):: r0, r1, r2, r3, rr, tt
!!  real(8),dimension(:,:,:,:,:,:),allocatable:: psig
!!  character(len=*),parameter:: subname='flatten_at_boundaries'
!!  
!!
!!  istc=1
!!  istf=1
!!
!!  allocate(psig(0:lr%d%n1,2,0:lr%d%n2,2,0:lr%d%n3,2), stat=istat)
!!  call memocc(istat, psig, 'psig', subname)
!!  call to_zero(8*(lr%d%n1+1)*(lr%d%n2+1)*(lr%d%n3+1), psig(0,1,0,1,0,1))
!!
!!  istf = istf + lr%wfd%nvctr_c
!!  call uncompress(lr%d%n1, lr%d%n2, lr%d%n3, &
!!       lr%wfd%nseg_c, lr%wfd%nvctr_c, lr%wfd%keygloc, lr%wfd%keyvloc,  &
!!       lr%wfd%nseg_f, lr%wfd%nvctr_f, &
!!       lr%wfd%keygloc(1,lr%wfd%nseg_c+min(1,lr%wfd%nseg_f)), &
!!       lr%wfd%keyvloc(lr%wfd%nseg_c+min(1,lr%wfd%nseg_f)), &
!!       psi(istc), psi(istf), psig)
!!
!!  r0 = (lr%locrad-8.d0*hx)**2
!!  !r0 = (lr%locrad-2.d0)**2
!!  do i3=0,lr%d%n3
!!      ii3=lr%ns3+i3
!!      r3 = (dble(ii3)*hz-lr%locregCenter(3))**2
!!      do i2=0,lr%d%n2
!!          ii2=lr%ns2+i2
!!          r2 = (dble(ii2)*hy-lr%locregCenter(2))**2
!!          do i1=0,lr%d%n1
!!              ii1=lr%ns1+i1
!!              r1 = (dble(ii1)*hx-lr%locregCenter(1))**2
!!              rr=r1+r2+r3
!!              !write(999,'(5es14.3)') r1, r2, r3, rr, r0
!!              if(rr>=r0) then
!!                  psig(i1,1,i2,1,i3,1)=0.d0
!!                  psig(i1,2,i2,1,i3,1)=0.d0
!!                  psig(i1,1,i2,2,i3,1)=0.d0
!!                  psig(i1,2,i2,2,i3,1)=0.d0
!!                  psig(i1,1,i2,1,i3,2)=0.d0
!!                  psig(i1,2,i2,1,i3,2)=0.d0
!!                  psig(i1,1,i2,2,i3,2)=0.d0
!!                  psig(i1,2,i2,2,i3,2)=0.d0
!!              else
!!                  !write(*,*) 'not zero'
!!              end if
!!          end do
!!      end do
!!  end do
!!
!!  call compress(lr%d%n1, lr%d%n2, &
!!       0, lr%d%n1, 0, lr%d%n2, 0, lr%d%n3, &
!!       lr%wfd%nseg_c, lr%wfd%nvctr_c, lr%wfd%keygloc, lr%wfd%keyvloc,  &
!!       lr%wfd%nseg_f, lr%wfd%nvctr_f, &
!!       lr%wfd%keygloc(1,lr%wfd%nseg_c+min(1,lr%wfd%nseg_f)), &
!!       lr%wfd%keyvloc(lr%wfd%nseg_c+min(1,lr%wfd%nseg_f)),  &
!!       psig, psi(istc), psi(istf))
!!
!!  istf = istf + 7*lr%wfd%nvctr_f
!!  istc = istc + lr%wfd%nvctr_c + 7*lr%wfd%nvctr_f
!!
!!  iall=-product(shape(psig))*kind(psig)
!!  deallocate(psig,stat=istat)
!!  call memocc(istat,iall,'psig',subname)
!!
!!
!!
!!end subroutine cut_at_boundaries2










!!subroutine flatten_at_boundaries2(lr, orbs, hx, hy, hz, psi)
!!  use module_base
!!  use module_types
!!  use module_interfaces
!!  implicit none
!!
!!  ! Calling arguments
!!  type(locreg_descriptors),intent(in):: lr
!!  type(orbitals_data),intent(in):: orbs
!!  real(8),intent(in):: hx, hy, hz
!!  real(8),dimension(orbs%npsidim_orbs),intent(inout):: psi
!!
!!  ! Local variables
!!  integer:: istc, istf, iorb, i1, i2, i3, istat, iall, ii1, ii2, ii3
!!  real(8):: r0, r1, r2, r3, rr, tt
!!  real(8),dimension(:,:,:,:,:,:),allocatable:: psig
!!  character(len=*),parameter:: subname='flatten_at_boundaries'
!!  
!!
!!  istc=1
!!  istf=1
!!
!!  allocate(psig(0:lr%d%n1,2,0:lr%d%n2,2,0:lr%d%n3,2), stat=istat)
!!  call memocc(istat, psig, 'psig', subname)
!!  call to_zero(8*(lr%d%n1+1)*(lr%d%n2+1)*(lr%d%n3+1), psig(0,1,0,1,0,1))
!!
!!  istf = istf + lr%wfd%nvctr_c
!!  call uncompress(lr%d%n1, lr%d%n2, lr%d%n3, &
!!       lr%wfd%nseg_c, lr%wfd%nvctr_c, lr%wfd%keygloc, lr%wfd%keyvloc,  &
!!       lr%wfd%nseg_f, lr%wfd%nvctr_f, &
!!       lr%wfd%keygloc(1,lr%wfd%nseg_c+min(1,lr%wfd%nseg_f)), &
!!       lr%wfd%keyvloc(lr%wfd%nseg_c+min(1,lr%wfd%nseg_f)), &
!!       psi(istc), psi(istf), psig)
!!
!!  r0=lr%locrad**2/5.d0
!!  !r0 = (lr%locrad-2.d0)**2
!!  do i3=0,lr%d%n3
!!      ii3=lr%ns3+i3
!!      r3 = (dble(ii3)*hz-lr%locregCenter(3))**2
!!      do i2=0,lr%d%n2
!!          ii2=lr%ns2+i2
!!          r2 = (dble(ii2)*hy-lr%locregCenter(2))**2
!!          do i1=0,lr%d%n1
!!              ii1=lr%ns1+i1
!!              r1 = (dble(ii1)*hx-lr%locregCenter(1))**2
!!              rr=r1+r2+r3
!!              tt=exp(-(rr-lr%locrad**2/2.d0)/r0)
!!              if(tt<1.d0) then
!!                  psig(i1,1,i2,1,i3,1)=tt*psig(i1,1,i2,1,i3,1)
!!                  psig(i1,2,i2,1,i3,1)=tt*psig(i1,2,i2,1,i3,1)
!!                  psig(i1,1,i2,2,i3,1)=tt*psig(i1,1,i2,2,i3,1)
!!                  psig(i1,2,i2,2,i3,1)=tt*psig(i1,2,i2,2,i3,1)
!!                  psig(i1,1,i2,1,i3,2)=tt*psig(i1,1,i2,1,i3,2)
!!                  psig(i1,2,i2,1,i3,2)=tt*psig(i1,2,i2,1,i3,2)
!!                  psig(i1,1,i2,2,i3,2)=tt*psig(i1,1,i2,2,i3,2)
!!                  psig(i1,2,i2,2,i3,2)=tt*psig(i1,2,i2,2,i3,2)
!!              end if
!!          end do
!!      end do
!!  end do
!!
!!  call compress(lr%d%n1, lr%d%n2, &
!!       0, lr%d%n1, 0, lr%d%n2, 0, lr%d%n3, &
!!       lr%wfd%nseg_c, lr%wfd%nvctr_c, lr%wfd%keygloc, lr%wfd%keyvloc,  &
!!       lr%wfd%nseg_f, lr%wfd%nvctr_f, &
!!       lr%wfd%keygloc(1,lr%wfd%nseg_c+min(1,lr%wfd%nseg_f)), &
!!       lr%wfd%keyvloc(lr%wfd%nseg_c+min(1,lr%wfd%nseg_f)),  &
!!       psig, psi(istc), psi(istf))
!!
!!  istf = istf + 7*lr%wfd%nvctr_f
!!  istc = istc + lr%wfd%nvctr_c + 7*lr%wfd%nvctr_f
!!
!!  iall=-product(shape(psig))*kind(psig)
!!  deallocate(psig,stat=istat)
!!  call memocc(istat,iall,'psig',subname)
!!
!!
!!
!!end subroutine flatten_at_boundaries2





!!subroutine get_both_gradients(iproc, nproc, lzd, orbs, psi, gnrm_in, gnrm_out)
!!  use module_base
!!  use module_types
!!  use module_interfaces
!!  implicit none
!!
!!  ! Calling arguments
!!  integer,intent(in):: iproc, nproc
!!  type(local_zone_descriptors),intent(in):: lzd
!!  type(orbitals_data),intent(in):: orbs
!!  real(8),dimension(orbs%npsidim_orbs),intent(in):: psi
!!  real(8),intent(out):: gnrm_out, gnrm_in
!!
!!  ! Local variables
!!  integer:: istc, istf, iorb, iiorb, ilr, i1, i2, i3, istat, iall, ierr, ii1, ii2, ii3, ipts_out, ipts_in
!!  real(8):: r0, r1, r2, r3, rr, tt, pts_in, pts_out
!!  real(8),dimension(:,:,:,:,:,:),allocatable:: psig
!!  character(len=*),parameter:: subname='flatten_at_boundaries'
!!  
!!
!!  istc=1
!!  istf=1
!!  gnrm_out=0.d0
!!  gnrm_in=0.d0
!!  pts_in=0.d0
!!  pts_out=0.d0
!!  do iorb=1,orbs%norbp
!!      iiorb=orbs%isorb+iorb
!!      ilr=orbs%inwhichlocreg(iiorb)
!!
!!      allocate(psig(0:lzd%llr(ilr)%d%n1,2,0:lzd%llr(ilr)%d%n2,2,0:lzd%llr(ilr)%d%n3,2), stat=istat)
!!      call memocc(istat, psig, 'psig', subname)
!!      call to_zero(8*(lzd%llr(ilr)%d%n1+1)*(lzd%llr(ilr)%d%n2+1)*(lzd%llr(ilr)%d%n3+1), psig(0,1,0,1,0,1))
!!
!!      istf = istf + lzd%llr(ilr)%wfd%nvctr_c
!!      call uncompress(lzd%llr(ilr)%d%n1, lzd%llr(ilr)%d%n2, lzd%llr(ilr)%d%n3, &
!!           lzd%llr(ilr)%wfd%nseg_c, lzd%llr(ilr)%wfd%nvctr_c, lzd%llr(ilr)%wfd%keygloc, lzd%llr(ilr)%wfd%keyvloc,  &
!!           lzd%llr(ilr)%wfd%nseg_f, lzd%llr(ilr)%wfd%nvctr_f, &
!!           lzd%llr(ilr)%wfd%keygloc(1,lzd%llr(ilr)%wfd%nseg_c+min(1,lzd%llr(ilr)%wfd%nseg_f)), &
!!           lzd%llr(ilr)%wfd%keyvloc(lzd%llr(ilr)%wfd%nseg_c+min(1,lzd%llr(ilr)%wfd%nseg_f)), &
!!           psi(istc), psi(istf), psig)
!!
!!      ipts_in=0
!!      ipts_out=0
!!      !r0=(lzd%llr(ilr)%locrad-3.d0)**2
!!      r0=(lzd%llr(ilr)%locrad-8.d0*lzd%hgrids(1))**2
!!      do i3=0,lzd%llr(ilr)%d%n3
!!          ii3=lzd%llr(ilr)%ns3+i3
!!          r3 = (dble(ii3)*lzd%hgrids(3)-lzd%llr(ilr)%locregCenter(3))**2
!!          do i2=0,lzd%llr(ilr)%d%n2
!!              ii2=lzd%llr(ilr)%ns2+i2
!!              r2 = (dble(ii2)*lzd%hgrids(2)-lzd%llr(ilr)%locregCenter(2))**2
!!              do i1=0,lzd%llr(ilr)%d%n1
!!                  ii1=lzd%llr(ilr)%ns1+i1
!!                  r1 = (dble(ii1)*lzd%hgrids(1)-lzd%llr(ilr)%locregCenter(1))**2
!!                  rr=r1+r2+r3
!!                  if(rr>r0) then
!!                      gnrm_out=gnrm_out+psig(i1,1,i2,1,i3,1)**2
!!                      gnrm_out=gnrm_out+psig(i1,2,i2,1,i3,1)**2
!!                      gnrm_out=gnrm_out+psig(i1,1,i2,2,i3,1)**2
!!                      gnrm_out=gnrm_out+psig(i1,2,i2,2,i3,1)**2
!!                      gnrm_out=gnrm_out+psig(i1,1,i2,1,i3,2)**2
!!                      gnrm_out=gnrm_out+psig(i1,2,i2,1,i3,2)**2
!!                      gnrm_out=gnrm_out+psig(i1,1,i2,2,i3,2)**2
!!                      gnrm_out=gnrm_out+psig(i1,2,i2,2,i3,2)**2
!!                      if(psig(i1,1,i2,1,i3,1)/=0.d0) then
!!                          ! point carries scaling function
!!                          pts_out=pts_out+1.d0
!!                          ipts_out=ipts_out+1
!!                      end if
!!                      if(psig(i1,2,i2,1,i3,1)/=0.d0) then
!!                          ! point carries wavelets
!!                          pts_out=pts_out+7.d0
!!                          ipts_out=ipts_out+7
!!                      end if
!!                  else
!!                      gnrm_in=gnrm_in+psig(i1,1,i2,1,i3,1)**2
!!                      gnrm_in=gnrm_in+psig(i1,2,i2,1,i3,1)**2
!!                      gnrm_in=gnrm_in+psig(i1,1,i2,2,i3,1)**2
!!                      gnrm_in=gnrm_in+psig(i1,2,i2,2,i3,1)**2
!!                      gnrm_in=gnrm_in+psig(i1,1,i2,1,i3,2)**2
!!                      gnrm_in=gnrm_in+psig(i1,2,i2,1,i3,2)**2
!!                      gnrm_in=gnrm_in+psig(i1,1,i2,2,i3,2)**2
!!                      gnrm_in=gnrm_in+psig(i1,2,i2,2,i3,2)**2
!!                      if(psig(i1,1,i2,1,i3,1)/=0.d0) then
!!                          ! point carries scaling function
!!                          pts_in=pts_in+1.d0
!!                          ipts_in=ipts_in+1
!!                      end if
!!                      if(psig(i1,2,i2,1,i3,1)/=0.d0) then
!!                          ! point carries wavelets
!!                          pts_in=pts_in+7.d0
!!                          ipts_in=ipts_in+7
!!                      end if
!!                  end if
!!              end do
!!          end do
!!      end do
!!
!!      if (ipts_in+ipts_out /= lzd%llr(ilr)%wfd%nvctr_c+7*lzd%llr(ilr)%wfd%nvctr_f) then
!!          write(*,'(2(a,i0))') 'ERROR: ',ipts_in+ipts_out,&
!!                      ' = ipts_in+ipts_out /= lzd%llr(ilr)%wfd%nvctr_c+7*lzd%llr(ilr)%wfd%nvctr_f = ',&
!!                      lzd%llr(ilr)%wfd%nvctr_c+7*lzd%llr(ilr)%wfd%nvctr_f
!!          stop
!!      end if
!!
!!      !!call compress(lzd%llr(ilr)%d%n1, lzd%llr(ilr)%d%n2, &
!!      !!     0, lzd%llr(ilr)%d%n1, 0, lzd%llr(ilr)%d%n2, 0, lzd%llr(ilr)%d%n3, &
!!      !!     lzd%llr(ilr)%wfd%nseg_c, lzd%llr(ilr)%wfd%nvctr_c, lzd%llr(ilr)%wfd%keygloc, lzd%llr(ilr)%wfd%keyvloc,  &
!!      !!     lzd%llr(ilr)%wfd%nseg_f, lzd%llr(ilr)%wfd%nvctr_f, &
!!      !!     lzd%llr(ilr)%wfd%keygloc(1,lzd%llr(ilr)%wfd%nseg_c+min(1,lzd%llr(ilr)%wfd%nseg_f)), &
!!      !!     lzd%llr(ilr)%wfd%keyvloc(lzd%llr(ilr)%wfd%nseg_c+min(1,lzd%llr(ilr)%wfd%nseg_f)),  &
!!      !!     psig, psi(istc), psi(istf))
!!
!!      istf = istf + 7*lzd%llr(ilr)%wfd%nvctr_f
!!      istc = istc + lzd%llr(ilr)%wfd%nvctr_c + 7*lzd%llr(ilr)%wfd%nvctr_f
!!
!!      iall=-product(shape(psig))*kind(psig)
!!      deallocate(psig,stat=istat)
!!      call memocc(istat,iall,'psig',subname)
!!
!!
!!  end do
!!
!!  !if(pts_in>0) gnrm_in=gnrm_in/pts_in
!!  !if(pts_out>0) gnrm_out=gnrm_out/pts_out
!!  call mpiallred(gnrm_out, 1, mpi_sum, mpi_comm_world, ierr)
!!  call mpiallred(gnrm_in, 1, mpi_sum, mpi_comm_world, ierr)
!!  call mpiallred(pts_in, 1, mpi_sum, mpi_comm_world, ierr)
!!  call mpiallred(pts_out, 1, mpi_sum, mpi_comm_world, ierr)
!!  gnrm_out=sqrt(gnrm_out/dble(orbs%norb))
!!  gnrm_in=sqrt(gnrm_in/dble(orbs%norb))
!!  if(iproc==0) write(*,'(a,5es14.4)') 'pts_in, pts_out, gnrm_in, gnrm_out, gnrm_in/gnrm_out', &
!!      pts_in, pts_out, gnrm_in, gnrm_out, gnrm_in/gnrm_out
!!
!!end subroutine get_both_gradients





!!subroutine create_penalty_basis_function(iproc, nproc, orbs, lzd, hx, hy, hz, confdatarr, psi)
!!use module_base
!!use module_types
!!use module_interfaces, except_this_one => apply_position_operators
!!implicit none
!!
!!! Calling arguments
!!integer,intent(in):: iproc, nproc
!!type(orbitals_data),intent(in):: orbs
!!type(local_zone_descriptors),intent(in):: lzd
!!real(8),intent(in):: hx, hy, hz
!!type(confpot_data),dimension(orbs%norbp),intent(in):: confdatarr
!!real(8),dimension(max(orbs%npsidim_orbs,orbs%npsidim_comp)),intent(inout):: psi
!!
!!! Local variables
!!integer:: oidx, iorb, ilr, npot, icenter, i_stat, i_all, ist_c, ist_f, ist, iiorb, iall, ierr, owa, owanext
!!real(8):: hxh, hyh, hzh, ddot, tt, t1, t2, time
!!real(8),dimension(:,:),allocatable:: psir
!!type(workarr_sumrho):: work_sr
!!real(8),dimension(0:3),parameter:: scal=1.d0
!!real(8),dimension(:,:,:),allocatable:: ypsitemp_c
!!real(8),dimension(:,:,:,:),allocatable:: ypsitemp_f
!!character(len=*),parameter:: subname='apply_position_operators'
!!integer, dimension(3) :: ishift !temporary variable in view of wavefunction creation
!!!!interface
!!!!subroutine position_operator(iproc, n1, n2, n3, nl1, nl2, nl3, nbuf, nspinor, psir, &
!!!!     hxh, hyh, hzh, dir, &
!!!!     ibyyzz_r) !optional
!!!!use module_base
!!!!implicit none
!!!!integer, intent(in) :: iproc, n1,n2,n3,nl1,nl2,nl3,nbuf,nspinor
!!!!real(wp), dimension(-14*nl1:2*n1+1+15*nl1,-14*nl2:2*n2+1+15*nl2,-14*nl3:2*n3+1+15*nl3,nspinor), intent(inout) :: psir
!!!!real(8),intent(in):: hxh, hyh, hzh
!!!!character(len=1),intent(in):: dir
!!!!integer, dimension(2,-14:2*n2+16,-14:2*n3+16), intent(in), optional :: ibyyzz_r
!!!!end subroutine
!!!!end interface
!!
!!  ishift=(/0,0,0/)
!!
!!  oidx = 0
!!  do iorb=1,orbs%norbp
!!     iiorb=orbs%isorb+iorb
!!     ilr = orbs%inwhichlocreg(iiorb)
!!     owa = orbs%onwhichatom(iiorb)
!!     if(iiorb<orbs%norb) then
!!         owanext=orbs%onwhichatom(iiorb+1)
!!     else
!!         owanext=lzd%nlr+1
!!     end if
!!
!!     if(owa/=owanext) then
!!         ! last basis function
!!
!!  
!!         !initialise the work arrays
!!         call initialize_work_arrays_sumrho(lzd%llr(ilr), work_sr)
!!
!!         ! Wavefunction in real space
!!         allocate(psir(Lzd%Llr(ilr)%d%n1i*Lzd%Llr(ilr)%d%n2i*Lzd%Llr(ilr)%d%n3i,orbs%nspinor+ndebug),stat=i_stat)
!!         call memocc(i_stat,psir,'psir',subname)
!!         call razero(Lzd%Llr(ilr)%d%n1i*Lzd%Llr(ilr)%d%n2i*Lzd%Llr(ilr)%d%n3i*orbs%nspinor,psir)
!!
!!         call daub_to_isf(lzd%llr(ilr), work_sr, psi(1+oidx), psir)
!!
!!         !!do i_stat=1,Lzd%Llr(ilr)%d%n1i*Lzd%Llr(ilr)%d%n2i*Lzd%Llr(ilr)%d%n3i
!!         !!    write(1000+iproc,'(i9,es18.7,i9)') i_stat, psir(i_stat,1), Lzd%Llr(ilr)%d%n1i*Lzd%Llr(ilr)%d%n2i*Lzd%Llr(ilr)%d%n3i
!!         !!end do
!!         !apply the potential to the psir wavefunction and calculate potential energy
!!         hxh=.5d0*hx
!!         hyh=.5d0*hy
!!         hzh=.5d0*hz
!!         !icenter=confinementCenter(iorb)
!!         !components of the potential
!!         npot=orbs%nspinor
!!         if (orbs%nspinor == 2) npot=1
!!
!!         !!call apply_confinement(iproc, lzd%llr(ilr)%d%n1,lzd%llr(ilr)%d%n2,lzd%llr(ilr)%d%n3,1,1,1,0,orbs%nspinor, psir, &
!!         !!     rxyz(1,icenter), hxh, hyh, hzh, lin%potentialprefac(at%iatype(icenter)), lin%confpotorder, &
!!         !!     lzd%llr(ilr)%nsi1, lzd%llr(ilr)%nsi2, lzd%llr(ilr)%nsi3,  &
!!         !!     lzd%llr(ilr)%bounds%ibyyzz_r) !optional
!!         write(*,*) 'confdatarr(iorb)%prefac',confdatarr(iorb)%prefac
!!        call penalty_basis_function(lzd%llr(ilr)%d%n1i,lzd%llr(ilr)%d%n2i,lzd%llr(ilr)%d%n3i,&
!!             lzd%llr(ilr)%d%n1i,lzd%llr(ilr)%d%n2i,lzd%llr(ilr)%d%n3i,&
!!             ishift,lzd%llr(ilr)%d%n2,lzd%llr(ilr)%d%n3,&
!!             orbs%nspinor,psir(1,1),&
!!             confdata=confdatarr(iorb),ibyyzz_r=lzd%llr(ilr)%bounds%ibyyzz_r)
!!
!!         call isf_to_daub(lzd%llr(ilr), work_sr, psir, psi(1+oidx))
!!
!!
!!
!!         i_all=-product(shape(psir))*kind(psir)
!!         deallocate(psir,stat=i_stat)
!!         call memocc(i_stat,i_all,'psir',subname)
!!
!!
!!         call deallocate_work_arrays_sumrho(work_sr)
!!
!!     end if
!!
!!     oidx = oidx + (Lzd%Llr(ilr)%wfd%nvctr_c+7*Lzd%Llr(ilr)%wfd%nvctr_f)*orbs%nspinor
!!
!!  enddo
!!
!!
!!end subroutine create_penalty_basis_function



!!subroutine penalty_basis_function(n1i,n2i,n3i,n1ip,n2ip,n3ip,ishift,n2,n3,nspinor,psir,&
!!     confdata,ibyyzz_r) !optional
!!  use module_base
!!  use module_types
!!  implicit none
!!  integer, intent(in) :: n1i,n2i,n3i,n1ip,n2ip,n3ip,n2,n3,nspinor
!!  integer, dimension(3), intent(in) :: ishift !<offset of potential box in wfn box coords.
!!  real(wp), dimension(n1i,n2i,n3i,nspinor), intent(inout) :: psir !< real-space wfn in lr
!!  type(confpot_data), intent(in), optional :: confdata !< data for the confining potential
!!  integer, dimension(2,-14:2*n2+16,-14:2*n3+16), intent(in), optional :: ibyyzz_r !< bounds in lr
!!  !local variables
!!  integer :: i1,i2,i3,ispinor,i1s,i1e,i2s,i2e,i3s,i3e,i1st,i1et
!!  real(wp) :: tt11,tt22,tt33,tt44,tt13,tt14,tt23,tt24,tt31,tt32,tt41,tt42,tt
!!  real(wp) :: psir1,psir2,psir3,psir4,pot1,pot2,pot3,pot4
!!
!!
!!  !write(*,*) 'present(confdata)', present(confdata)
!!  !write(*,*) 'confdata%prefac, confdata%potorder', confdata%prefac, confdata%potorder
!!  !write(*,*) 'n1ip*n2ip*n3ip', n1ip*n2ip*n3ip
!!
!!  !loop on wavefunction
!!  !calculate the limits in all the directions
!!  !regions in which both the potential and wavefunctions are defined
!!  i3s=max(1,ishift(3)+1)
!!  i3e=min(n3i,n3ip+ishift(3))
!!  i2s=max(1,ishift(2)+1)
!!  i2e=min(n2i,n2ip+ishift(2))
!!  i1s=max(1,ishift(1)+1)
!!  i1e=min(n1i,n1ip+ishift(1))
!!
!!
!!  !$omp parallel default(none)&
!!  !$omp shared(pot,psir,n1i,n2i,n3i,n1ip,n2ip,n3ip,n2,n3,epot,ibyyzz_r,nspinor)&
!!  !!$omp shared(psir,n1i,n2i,n3i,n1ip,n2ip,n3ip,n2,n3,ibyyzz_r,nspinor)& ! for jaguar
!!  !$omp shared(i1s,i1e,i2s,i2e,i3s,i3e,ishift)&
!!  !$omp private(ispinor,i1,i2,i3,epot_p,i1st,i1et)&
!!  !!$omp private(ispinor,i1,i2,i3,i1st,i1et)& for jaguar
!!  !$omp private(tt11,tt22,tt33,tt44,tt13,tt14,tt23,tt24,tt31,tt32,tt41,tt42,tt)&
!!  !$omp private(psir1,psir2,psir3,psir4,pot1,pot2,pot3,pot4)
!!
!!!!$  !$omp parallel default(private)&
!!!!$  !$omp shared(pot,psir,n1i,n2i,n3i,n1ip,n2ip,n3ip,n2,n3,epot,ibyyzz_r,nspinor)&
!!!!$  !$omp shared(i1s,i1e,i2s,i2e,i3s,i3e,ishift)
!!  !case without bounds
!!
!!
!!  !put to zero the external part of psir if the potential is more little than the wavefunction
!!  !first part of the array
!!  do ispinor=1,nspinor
!!     !$omp do 
!!     do i3=1,i3s-1
!!        do i2=1,n2i
!!           do i1=1,n1i
!!             psir(i1,i2,i3,ispinor)=0.0_wp 
!!           end do
!!        end do
!!     end do
!!     !$omp end do
!!  end do
!!
!!  !central part of the array
!!  do ispinor=1,nspinor
!!     !$omp do 
!!     do i3=i3s,i3e
!!
!!        !first part
!!        do i2=1,i2s-1
!!           do i1=1,n1i
!!              psir(i1,i2,i3,ispinor)=0.0_wp 
!!           end do
!!        end do
!!        !central part
!!        do i2=i2s,i2e
!!           do i1=1,i1s-1
!!              psir(i1,i2,i3,ispinor)=0.0_wp 
!!           end do
!!           do i1=i1e+1,n1i
!!              psir(i1,i2,i3,ispinor)=0.0_wp 
!!           end do
!!        end do
!!        !last part
!!        do i2=i2e+1,n2i
!!           do i1=1,n1i
!!              psir(i1,i2,i3,ispinor)=0.0_wp 
!!           end do
!!        end do
!!
!!     end do
!!     !$omp end do
!!  end do
!!
!!
!!  !last part of the array
!!  do ispinor=1,nspinor
!!     !$omp do 
!!     do i3=i3e+1,n3i
!!        do i2=1,n2i
!!           do i1=1,n1i
!!              psir(i1,i2,i3,ispinor)=0.0_wp 
!!           end do
!!        end do
!!     end do
!!     !$omp end do
!!  end do
!!
!!
!!  !important part of the array
!!  if (nspinor==4) then
!!     stop 'not yet implemented for nspinor==4!'
!!     !!!!$omp do
!!     !!!do i3=i3s,i3e
!!     !!!   do i2=i2s,i2e
!!     !!!      !thanks to the optional argument the conditional is done at compile time
!!     !!!      if (present(ibyyzz_r)) then
!!     !!!         i1st=max(i1s,ibyyzz_r(1,i2-15,i3-15)+1) !in bounds coordinates
!!     !!!         i1et=min(i1e,ibyyzz_r(2,i2-15,i3-15)+1) !in bounds coordinates
!!     !!!      else
!!     !!!         i1st=i1s
!!     !!!         i1et=i1e
!!     !!!      end if
!!     !!!      !no need of setting up to zero values outside wavefunction bounds
!!     !!!      do i1=i1st,i1et
!!     !!!         !wavefunctions
!!     !!!         psir1=psir(i1,i2,i3,1)
!!     !!!         psir2=psir(i1,i2,i3,2)
!!     !!!         psir3=psir(i1,i2,i3,3)
!!     !!!         psir4=psir(i1,i2,i3,4)
!!     !!!         !potentials + confining term
!!     !!!         pot1=pot(i1-ishift(1),i2-ishift(2),i3-ishift(3),1)+cp(i1,i2,i3)
!!     !!!         pot2=pot(i1-ishift(1),i2-ishift(2),i3-ishift(3),2)+cp(i1,i2,i3)
!!     !!!         pot3=pot(i1-ishift(1),i2-ishift(2),i3-ishift(3),3)+cp(i1,i2,i3)
!!     !!!         pot4=pot(i1-ishift(1),i2-ishift(2),i3-ishift(3),4)+cp(i1,i2,i3)
!!
!!     !!!         !diagonal terms
!!     !!!         tt11=pot1*psir1 !p1
!!     !!!         tt22=pot1*psir2 !p2
!!     !!!         tt33=pot4*psir3 !p3
!!     !!!         tt44=pot4*psir4 !p4
!!     !!!         !Rab*Rb
!!     !!!         tt13=pot2*psir3 !p1
!!     !!!         !Iab*Ib
!!     !!!         tt14=pot3*psir4 !p1
!!     !!!         !Rab*Ib
!!     !!!         tt23=pot2*psir4 !p2
!!     !!!         !Iab*Rb
!!     !!!         tt24=pot3*psir3 !p2
!!     !!!         !Rab*Ra
!!     !!!         tt31=pot2*psir1 !p3
!!     !!!         !Iab*Ia
!!     !!!         tt32=pot3*psir2 !p3
!!     !!!         !Rab*Ia
!!     !!!         tt41=pot2*psir2 !p4
!!     !!!         !Iab*Ra
!!     !!!         tt42=pot3*psir1 !p4
!!
!!     !!!         !value of the potential energy
!!     !!!         epot_p=epot_p+tt11*psir1+tt22*psir2+tt33*psir3+tt44*psir4+&
!!     !!!              2.0_gp*tt31*psir3-2.0_gp*tt42*psir4+2.0_gp*tt41*psir4+2.0_gp*tt32*psir3
!!
!!     !!!         !wavefunction update
!!     !!!         !p1=h1p1+h2p3-h3p4
!!     !!!         !p2=h1p2+h2p4+h3p3
!!     !!!         !p3=h2p1+h3p2+h4p3
!!     !!!         !p4=h2p2-h3p1+h4p4
!!     !!!         psir(i1,i2,i3,1)=tt11+tt13-tt14
!!     !!!         psir(i1,i2,i3,2)=tt22+tt23+tt24
!!     !!!         psir(i1,i2,i3,3)=tt33+tt31+tt32
!!     !!!         psir(i1,i2,i3,4)=tt44+tt41-tt42
!!     !!!      end do
!!     !!!   end do
!!     !!!end do
!!     !!!!$omp end do
!!
!!  else !case with nspinor /=4
!!     do ispinor=1,nspinor
!!        !$omp do
!!        do i3=i3s,i3e
!!           do i2=i2s,i2e
!!              !thanks to the optional argument the conditional is done at compile time
!!              if (present(ibyyzz_r)) then
!!                 i1st=max(i1s,ibyyzz_r(1,i2-15,i3-15)+1) !in bounds coordinates
!!                 i1et=min(i1e,ibyyzz_r(2,i2-15,i3-15)+1) !in bounds coordinates
!!              else
!!                 i1st=i1s
!!                 i1et=i1e
!!              end if
!!              !no need of setting up to zero values outside wavefunction bounds
!!              !write(*,'(a,6i9)') 'i1st, i1et, i2s, i2e, i3s, i3e', i1st, i1et, i2s, i2e, i3s, i3e
!!              do i1=i1st,i1et
!!                 pot1=cp(i1,i2,i3)
!!                 psir(i1,i2,i3,ispinor)=pot1
!!              end do
!!           end do
!!        end do
!!        !$omp end do
!!     end do
!!  end if
!!  
!!  
!!  !$omp end parallel
!!
!!contains
!!  
!!  !inline the definition of the confining potential
!!  real(wp) function cp(i1,i2,i3)
!!    implicit none
!!    integer, intent(in) :: i1,i2,i3
!!    !local variables
!!    real(wp) :: r2
!!    !to be sure that the conditional is executed at compile time
!!    if (present(confdata)) then
!!       r2=(confdata%hh(1)*real(i1+confdata%ioffset(1),wp)-confdata%rxyzConf(1))**2 +&
!!            (confdata%hh(2)*real(i2+confdata%ioffset(2),wp)-confdata%rxyzConf(2))**2 +&
!!            (confdata%hh(3)*real(i3+confdata%ioffset(3),wp)-confdata%rxyzConf(3))**2 
!!       !if(r2>=81.d0) write(*,'(6i8,3es11.2,es13.4)') i1, i2, i3, confdata%ioffset(1), confdata%ioffset(2), confdata%ioffset(3), confdata%rxyzConf(1), confdata%rxyzConf(2), confdata%rxyzConf(3), r2 
!!
!!       cp=confdata%prefac*r2**(confdata%potorder/2)
!!    else
!!       cp=0.0_wp
!!    end if
!!
!!  end function cp
!!
!!END SUBROUTINE penalty_basis_function


!!subroutine flatten_at_boundaries(lzd, orbs, psi)
!!  use module_base
!!  use module_types
!!  use module_interfaces
!!  implicit none
!!
!!  ! Calling arguments
!!  type(local_zone_descriptors),intent(in):: lzd
!!  type(orbitals_data),intent(in):: orbs
!!  real(8),dimension(orbs%npsidim_orbs),intent(inout):: psi
!!
!!  ! Local variables
!!  integer:: istc, istf, iorb, iiorb, ilr, i1, i2, i3, istat, iall, ii1, ii2, ii3
!!  real(8):: r0, r1, r2, r3, rr, tt
!!  real(8),dimension(:,:,:,:,:,:),allocatable:: psig
!!  character(len=*),parameter:: subname='flatten_at_boundaries'
!!  
!!
!!  istc=1
!!  istf=1
!!  do iorb=1,orbs%norbp
!!      iiorb=orbs%isorb+iorb
!!      ilr=orbs%inwhichlocreg(iiorb)
!!
!!      allocate(psig(0:lzd%llr(ilr)%d%n1,2,0:lzd%llr(ilr)%d%n2,2,0:lzd%llr(ilr)%d%n3,2), stat=istat)
!!      call memocc(istat, psig, 'psig', subname)
!!      call to_zero(8*(lzd%llr(ilr)%d%n1+1)*(lzd%llr(ilr)%d%n2+1)*(lzd%llr(ilr)%d%n3+1), psig(0,1,0,1,0,1))
!!
!!      istf = istf + lzd%llr(ilr)%wfd%nvctr_c
!!      call uncompress(lzd%llr(ilr)%d%n1, lzd%llr(ilr)%d%n2, lzd%llr(ilr)%d%n3, &
!!           lzd%llr(ilr)%wfd%nseg_c, lzd%llr(ilr)%wfd%nvctr_c, lzd%llr(ilr)%wfd%keygloc, lzd%llr(ilr)%wfd%keyvloc,  &
!!           lzd%llr(ilr)%wfd%nseg_f, lzd%llr(ilr)%wfd%nvctr_f, &
!!           lzd%llr(ilr)%wfd%keygloc(1,lzd%llr(ilr)%wfd%nseg_c+min(1,lzd%llr(ilr)%wfd%nseg_f)), &
!!           lzd%llr(ilr)%wfd%keyvloc(lzd%llr(ilr)%wfd%nseg_c+min(1,lzd%llr(ilr)%wfd%nseg_f)), &
!!           psi(istc), psi(istf), psig)
!!
!!      r0=lzd%llr(ilr)%locrad**2/5.d0
!!      do i3=0,lzd%llr(ilr)%d%n3
!!          ii3=lzd%llr(ilr)%ns3+i3
!!          r3 = (dble(ii3)*lzd%hgrids(3)-lzd%llr(ilr)%locregCenter(3))**2
!!          do i2=0,lzd%llr(ilr)%d%n2
!!              ii2=lzd%llr(ilr)%ns2+i2
!!              r2 = (dble(ii2)*lzd%hgrids(2)-lzd%llr(ilr)%locregCenter(2))**2
!!              do i1=0,lzd%llr(ilr)%d%n1
!!                  ii1=lzd%llr(ilr)%ns1+i1
!!                  r1 = (dble(ii1)*lzd%hgrids(1)-lzd%llr(ilr)%locregCenter(1))**2
!!                  rr=r1+r2+r3
!!                  tt=exp(-(rr-lzd%llr(ilr)%locrad**2/2.d0)/r0)
!!                  tt=min(tt,1.d0)
!!                  psig(i1,1,i2,1,i3,1)=tt*psig(i1,1,i2,1,i3,1)
!!                  psig(i1,2,i2,1,i3,1)=tt*psig(i1,2,i2,1,i3,1)
!!                  psig(i1,1,i2,2,i3,1)=tt*psig(i1,1,i2,2,i3,1)
!!                  psig(i1,2,i2,2,i3,1)=tt*psig(i1,2,i2,2,i3,1)
!!                  psig(i1,1,i2,1,i3,2)=tt*psig(i1,1,i2,1,i3,2)
!!                  psig(i1,2,i2,1,i3,2)=tt*psig(i1,2,i2,1,i3,2)
!!                  psig(i1,1,i2,2,i3,2)=tt*psig(i1,1,i2,2,i3,2)
!!                  psig(i1,2,i2,2,i3,2)=tt*psig(i1,2,i2,2,i3,2)
!!              end do
!!          end do
!!      end do
!!
!!      call compress(lzd%llr(ilr)%d%n1, lzd%llr(ilr)%d%n2, &
!!           0, lzd%llr(ilr)%d%n1, 0, lzd%llr(ilr)%d%n2, 0, lzd%llr(ilr)%d%n3, &
!!           lzd%llr(ilr)%wfd%nseg_c, lzd%llr(ilr)%wfd%nvctr_c, lzd%llr(ilr)%wfd%keygloc, lzd%llr(ilr)%wfd%keyvloc,  &
!!           lzd%llr(ilr)%wfd%nseg_f, lzd%llr(ilr)%wfd%nvctr_f, &
!!           lzd%llr(ilr)%wfd%keygloc(1,lzd%llr(ilr)%wfd%nseg_c+min(1,lzd%llr(ilr)%wfd%nseg_f)), &
!!           lzd%llr(ilr)%wfd%keyvloc(lzd%llr(ilr)%wfd%nseg_c+min(1,lzd%llr(ilr)%wfd%nseg_f)),  &
!!           psig, psi(istc), psi(istf))
!!
!!      istf = istf + 7*lzd%llr(ilr)%wfd%nvctr_f
!!      istc = istc + lzd%llr(ilr)%wfd%nvctr_c + 7*lzd%llr(ilr)%wfd%nvctr_f
!!
!!      iall=-product(shape(psig))*kind(psig)
!!      deallocate(psig,stat=istat)
!!      call memocc(istat,iall,'psig',subname)
!!
!!
!!  end do
!!
!!end subroutine flatten_at_boundaries



!!subroutine get_weighted_gradient(iproc, nproc, lzd, orbs, psi)
!!  use module_base
!!  use module_types
!!  use module_interfaces
!!  implicit none
!!
!!  ! Calling arguments
!!  integer,intent(in):: iproc, nproc
!!  type(local_zone_descriptors),intent(in):: lzd
!!  type(orbitals_data),intent(in):: orbs
!!  real(8),dimension(orbs%npsidim_orbs),intent(in):: psi
!!
!!  ! Local variables
!!  integer:: istc, istf, iorb, iiorb, ilr, i1, i2, i3, istat, iall, ierr, ii1, ii2, ii3
!!  real(8):: r0, r1, r2, r3, rr, tt, gnrm
!!  real(8),dimension(:,:,:,:,:,:),allocatable:: psig
!!  character(len=*),parameter:: subname='flatten_at_boundaries'
!!  
!!
!!  istc=1
!!  istf=1
!!  gnrm=0.d0
!!  do iorb=1,orbs%norbp
!!      iiorb=orbs%isorb+iorb
!!      ilr=orbs%inwhichlocreg(iiorb)
!!
!!      allocate(psig(0:lzd%llr(ilr)%d%n1,2,0:lzd%llr(ilr)%d%n2,2,0:lzd%llr(ilr)%d%n3,2), stat=istat)
!!      call memocc(istat, psig, 'psig', subname)
!!      call to_zero(8*(lzd%llr(ilr)%d%n1+1)*(lzd%llr(ilr)%d%n2+1)*(lzd%llr(ilr)%d%n3+1), psig(0,1,0,1,0,1))
!!
!!      istf = istf + lzd%llr(ilr)%wfd%nvctr_c
!!      call uncompress(lzd%llr(ilr)%d%n1, lzd%llr(ilr)%d%n2, lzd%llr(ilr)%d%n3, &
!!           lzd%llr(ilr)%wfd%nseg_c, lzd%llr(ilr)%wfd%nvctr_c, lzd%llr(ilr)%wfd%keygloc, lzd%llr(ilr)%wfd%keyvloc,  &
!!           lzd%llr(ilr)%wfd%nseg_f, lzd%llr(ilr)%wfd%nvctr_f, &
!!           lzd%llr(ilr)%wfd%keygloc(1,lzd%llr(ilr)%wfd%nseg_c+min(1,lzd%llr(ilr)%wfd%nseg_f)), &
!!           lzd%llr(ilr)%wfd%keyvloc(lzd%llr(ilr)%wfd%nseg_c+min(1,lzd%llr(ilr)%wfd%nseg_f)), &
!!           psi(istc), psi(istf), psig)
!!
!!      r0=lzd%llr(ilr)%locrad**2/3.d0
!!      do i3=0,lzd%llr(ilr)%d%n3
!!          ii3=lzd%llr(ilr)%ns3+i3
!!          r3 = (dble(ii3)*lzd%hgrids(3)-lzd%llr(ilr)%locregCenter(3))**2
!!          do i2=0,lzd%llr(ilr)%d%n2
!!              ii2=lzd%llr(ilr)%ns2+i2
!!              r2 = (dble(ii2)*lzd%hgrids(2)-lzd%llr(ilr)%locregCenter(2))**2
!!              do i1=0,lzd%llr(ilr)%d%n1
!!                  ii1=lzd%llr(ilr)%ns1+i1
!!                  r1 = (dble(ii1)*lzd%hgrids(1)-lzd%llr(ilr)%locregCenter(1))**2
!!                  rr=r1+r2+r3
!!                  !tt=exp(-(rr-lzd%llr(ilr)%locrad**2/2.d0)/r0)
!!                  tt=exp(-rr/r0)
!!                  !tt=min(tt,1.d0)
!!                  gnrm=gnrm+tt*psig(i1,1,i2,1,i3,1)**2
!!                  gnrm=gnrm+tt*psig(i1,2,i2,1,i3,1)**2
!!                  gnrm=gnrm+tt*psig(i1,1,i2,2,i3,1)**2
!!                  gnrm=gnrm+tt*psig(i1,2,i2,2,i3,1)**2
!!                  gnrm=gnrm+tt*psig(i1,1,i2,1,i3,2)**2
!!                  gnrm=gnrm+tt*psig(i1,2,i2,1,i3,2)**2
!!                  gnrm=gnrm+tt*psig(i1,1,i2,2,i3,2)**2
!!                  gnrm=gnrm+tt*psig(i1,2,i2,2,i3,2)**2
!!              end do
!!          end do
!!      end do
!!
!!      !!call compress(lzd%llr(ilr)%d%n1, lzd%llr(ilr)%d%n2, &
!!      !!     0, lzd%llr(ilr)%d%n1, 0, lzd%llr(ilr)%d%n2, 0, lzd%llr(ilr)%d%n3, &
!!      !!     lzd%llr(ilr)%wfd%nseg_c, lzd%llr(ilr)%wfd%nvctr_c, lzd%llr(ilr)%wfd%keygloc, lzd%llr(ilr)%wfd%keyvloc,  &
!!      !!     lzd%llr(ilr)%wfd%nseg_f, lzd%llr(ilr)%wfd%nvctr_f, &
!!      !!     lzd%llr(ilr)%wfd%keygloc(1,lzd%llr(ilr)%wfd%nseg_c+min(1,lzd%llr(ilr)%wfd%nseg_f)), &
!!      !!     lzd%llr(ilr)%wfd%keyvloc(lzd%llr(ilr)%wfd%nseg_c+min(1,lzd%llr(ilr)%wfd%nseg_f)),  &
!!      !!     psig, psi(istc), psi(istf))
!!
!!      istf = istf + 7*lzd%llr(ilr)%wfd%nvctr_f
!!      istc = istc + lzd%llr(ilr)%wfd%nvctr_c + 7*lzd%llr(ilr)%wfd%nvctr_f
!!
!!      iall=-product(shape(psig))*kind(psig)
!!      deallocate(psig,stat=istat)
!!      call memocc(istat,iall,'psig',subname)
!!
!!
!!  end do
!!
!!  call mpiallred(gnrm, 1, mpi_sum, mpi_comm_world, ierr)
!!  gnrm=gnrm/dble(orbs%norb)
!!  if(iproc==0) write(*,*) 'weighted grnm',gnrm
!!
!!end subroutine get_weighted_gradient


!!subroutine update_kernel(norb, Umat, kernel)
!!  use module_base
!!  use module_types
!!  implicit none
!!  
!!  ! Calling arguments
!!  integer,intent(in):: norb
!!  real(8),dimension(norb,norb),intent(in):: Umat
!!  real(8),dimension(norb,norb),intent(inout):: kernel
!!  
!!  ! Local variables
!!  integer:: iorb, jorb, korb, lorb, istat, iall
!!  real(8):: tt
!!  real(8),dimension(:,:),allocatable:: kernelold
!!  character(len=*),parameter:: subname='update_kernel'
!!  
!!  allocate(kernelold(norb,norb), stat=istat)
!!  call memocc(istat, kernelold, 'kernelold', subname)
!!  
!!  call dcopy(norb**2, kernel(1,1), 1, kernelold(1,1), 1)
!!  do iorb=1,norb
!!      do jorb=1,norb
!!          tt=0.d0
!!          do korb=1,norb
!!              do lorb=1,norb
!!                  tt=tt+kernelold(korb,lorb)*Umat(korb,iorb)*Umat(lorb,jorb)
!!                  !tt=tt+kernelold(korb,lorb)*Umat(iorb,korb)*Umat(jorb,lorb)
!!              end do
!!          end do
!!          kernel(jorb,iorb)=tt
!!      end do
!!  end do
!!  
!!  iall=-product(shape(kernelold))*kind(kernelold)
!!  deallocate(kernelold, stat=istat)
!!  call memocc(istat, iall, 'kernelold', subname)
!!
!!end subroutine update_kernel


!!subroutine check_locregCenters(iproc, lzd, locregCenter, hx, hy, hz)
!!  use module_base
!!  use module_types
!!  implicit none
!!  
!!  ! Calling arguments
!!  integer,intent(in):: iproc
!!  type(local_zone_descriptors),intent(in):: lzd
!!  real(8),dimension(3,lzd%nlr),intent(in):: locregCenter
!!  real(8),intent(in):: hx, hy, hz
!!  
!!  ! Local variables
!!  integer:: ilr, ierr
!!  
!!  do ilr=1,lzd%nlr
!!      if( floor(locregCenter(1,ilr)/hx) < 0 .or. ceiling(locregCenter(1,ilr)/hx) > lzd%glr%d%n1 ) then
!!          if(iproc==0) then
!!              write(*,'(1x,a,i0,a,i0,1x,i0,a,i0,1x,i0)') 'ERROR: new center for locreg ',ilr,&
!!                  ' is outside of box in x direction! Box limits=',0,lzd%glr%d%n1,&
!!                  ', center=',floor(locregCenter(1,ilr)/hx),ceiling(locregCenter(1,ilr)/hx)
!!          end if
!!          call mpi_barrier(mpi_comm_world, ierr)
!!          stop
!!      end if
!!      if( floor(locregCenter(2,ilr)/hy) < 0 .or. ceiling(locregCenter(2,ilr)/hy) > lzd%glr%d%n2 ) then
!!          if(iproc==0) then
!!              write(*,'(1x,a,i0,a,i0,1x,i0,a,i0,1x,i0)') 'ERROR: new center for locreg ',ilr,&
!!                  'is outside of box in y direction! Box limits=',0,lzd%glr%d%n2,&
!!                  ', center=',floor(locregCenter(2,ilr)/hy),ceiling(locregCenter(2,ilr)/hy)
!!          end if
!!          call mpi_barrier(mpi_comm_world, ierr)
!!          stop
!!      end if
!!      if( floor(locregCenter(3,ilr)/hz) < 0 .or. ceiling(locregCenter(3,ilr)/hz) > lzd%glr%d%n3 ) then
!!          if(iproc==0) then
!!              write(*,'(1x,a,i0,a,i0,1x,i0,a,i0,1x,i0)') 'ERROR: new center for locreg ',ilr,&
!!                  'is outside of box in z direction! Box limits=',0,lzd%glr%d%n3,&
!!                  ', center=',floor(locregCenter(3,ilr)/hz),ceiling(locregCenter(3,ilr)/hz)
!!          end if
!!          call mpi_barrier(mpi_comm_world, ierr)
!!          stop
!!      end if
!!  end do
!!
!!end subroutine check_locregCenters           




!!subroutine apply_position_operators(iproc, nproc, orbs, lzd, hx, hy, hz, confdatarr, psi, order, xpsi, ypsi, zpsi)
!!use module_base
!!use module_types
!!use module_interfaces, except_this_one => apply_position_operators
!!implicit none
!!
!!! Calling arguments
!!integer,intent(in):: iproc, nproc, order
!!type(orbitals_data),intent(in):: orbs
!!type(local_zone_descriptors),intent(in):: lzd
!!real(8),intent(in):: hx, hy, hz
!!type(confpot_data),dimension(orbs%norbp),intent(in):: confdatarr
!!real(8),dimension(max(orbs%npsidim_orbs,orbs%npsidim_comp)),intent(in):: psi
!!real(8),dimension(max(orbs%npsidim_orbs,orbs%npsidim_comp)),intent(out):: xpsi, ypsi, zpsi
!!
!!! Local variables
!!integer:: oidx, iorb, ilr, npot, icenter, i_stat, i_all, ist_c, ist_f, ist, iiorb, iall, ierr
!!real(8):: hxh, hyh, hzh, ddot, tt, t1, t2, time
!!real(8),dimension(:,:),allocatable:: psir, psirx, psiry, psirz
!!type(workarr_sumrho):: work_sr
!!real(8),dimension(0:3),parameter:: scal=1.d0
!!real(8),dimension(:,:,:),allocatable:: ypsitemp_c
!!real(8),dimension(:,:,:,:),allocatable:: ypsitemp_f
!!character(len=*),parameter:: subname='apply_position_operators'
!!integer, dimension(3) :: ishift !temporary variable in view of wavefunction creation
!!!!interface
!!!!subroutine position_operator(iproc, n1, n2, n3, nl1, nl2, nl3, nbuf, nspinor, psir, &
!!!!     hxh, hyh, hzh, dir, &
!!!!     ibyyzz_r) !optional
!!!!use module_base
!!!!implicit none
!!!!integer, intent(in) :: iproc, n1,n2,n3,nl1,nl2,nl3,nbuf,nspinor
!!!!real(wp), dimension(-14*nl1:2*n1+1+15*nl1,-14*nl2:2*n2+1+15*nl2,-14*nl3:2*n3+1+15*nl3,nspinor), intent(inout) :: psir
!!!!real(8),intent(in):: hxh, hyh, hzh
!!!!character(len=1),intent(in):: dir
!!!!integer, dimension(2,-14:2*n2+16,-14:2*n3+16), intent(in), optional :: ibyyzz_r
!!!!end subroutine
!!!!end interface
!!
!!  ishift=(/0,0,0/)
!!
!!  call to_zero(max(orbs%npsidim_orbs,orbs%npsidim_comp), xpsi(1))
!!  call to_zero(max(orbs%npsidim_orbs,orbs%npsidim_comp), ypsi(1))
!!  call to_zero(max(orbs%npsidim_orbs,orbs%npsidim_comp), zpsi(1))
!!  oidx = 0
!!  do iorb=1,orbs%norbp
!!     ilr = orbs%inwhichlocreg(iorb+orbs%isorb)
!!
!!     iiorb=orbs%isorb+iorb
!!     !!write(*,'(a,4i8,4x,3i6)') 'iproc, iorb, iiorb, ilr, confdatarr(iorb)%ioffset(:)', &
!!     !!    iproc, iorb, iiorb, ilr, confdatarr(iorb)%ioffset(:)
!!     !!write(*,'(a,3i8,6i6)') 'iproc, iiorb, ilr, is1, ie1, is2, ie2, is3, ie3', &
!!     !!    1, lzd%llr(ilr)%d%n1i, 1, lzd%llr(ilr)%d%n2i, 1, lzd%llr(ilr)%d%n3i
!!  
!!     !initialise the work arrays
!!     call initialize_work_arrays_sumrho(lzd%llr(ilr), work_sr)
!!
!!     ! Wavefunction in real space
!!     allocate(psir(Lzd%Llr(ilr)%d%n1i*Lzd%Llr(ilr)%d%n2i*Lzd%Llr(ilr)%d%n3i,orbs%nspinor+ndebug),stat=i_stat)
!!     call memocc(i_stat,psir,'psir',subname)
!!     call razero(Lzd%Llr(ilr)%d%n1i*Lzd%Llr(ilr)%d%n2i*Lzd%Llr(ilr)%d%n3i*orbs%nspinor,psir)
!!
!!     allocate(psirx(Lzd%Llr(ilr)%d%n1i*Lzd%Llr(ilr)%d%n2i*Lzd%Llr(ilr)%d%n3i,orbs%nspinor+ndebug),stat=i_stat)
!!     call memocc(i_stat,psirx,'psirx',subname)
!!     call razero(Lzd%Llr(ilr)%d%n1i*Lzd%Llr(ilr)%d%n2i*Lzd%Llr(ilr)%d%n3i*orbs%nspinor,psirx)
!!
!!     allocate(psiry(Lzd%Llr(ilr)%d%n1i*Lzd%Llr(ilr)%d%n2i*Lzd%Llr(ilr)%d%n3i,orbs%nspinor+ndebug),stat=i_stat)
!!     call memocc(i_stat,psiry,'psiry',subname)
!!     call razero(Lzd%Llr(ilr)%d%n1i*Lzd%Llr(ilr)%d%n2i*Lzd%Llr(ilr)%d%n3i*orbs%nspinor,psiry)
!!
!!     allocate(psirz(Lzd%Llr(ilr)%d%n1i*Lzd%Llr(ilr)%d%n2i*Lzd%Llr(ilr)%d%n3i,orbs%nspinor+ndebug),stat=i_stat)
!!     call memocc(i_stat,psirz,'psirz',subname)
!!     call razero(Lzd%Llr(ilr)%d%n1i*Lzd%Llr(ilr)%d%n2i*Lzd%Llr(ilr)%d%n3i*orbs%nspinor,psirz)
!!
!!     !transform the wavefunction in Daubechies basis to the wavefunction in ISF basis
!!     !the psir wavefunction is given in the spinorial form
!!
!!     !psi(1+oidx+lzd%llr(ilr)%wfd%nvctr_c:1+oidx+lzd%llr(ilr)%wfd%nvctr_c+7*lzd%llr(ilr)%wfd%nvctr_f-1)=0.d0
!!
!!     call daub_to_isf(lzd%llr(ilr), work_sr, psi(1+oidx), psir)
!!
!!     !!do i_stat=1,Lzd%Llr(ilr)%d%n1i*Lzd%Llr(ilr)%d%n2i*Lzd%Llr(ilr)%d%n3i
!!     !!    write(1000+iproc,'(i9,es18.7,i9)') i_stat, psir(i_stat,1), Lzd%Llr(ilr)%d%n1i*Lzd%Llr(ilr)%d%n2i*Lzd%Llr(ilr)%d%n3i
!!     !!end do
!!     !apply the potential to the psir wavefunction and calculate potential energy
!!     hxh=.5d0*hx
!!     hyh=.5d0*hy
!!     hzh=.5d0*hz
!!     !icenter=confinementCenter(iorb)
!!     !components of the potential
!!     npot=orbs%nspinor
!!     if (orbs%nspinor == 2) npot=1
!!
!!     !!call apply_confinement(iproc, lzd%llr(ilr)%d%n1,lzd%llr(ilr)%d%n2,lzd%llr(ilr)%d%n3,1,1,1,0,orbs%nspinor, psir, &
!!     !!     rxyz(1,icenter), hxh, hyh, hzh, lin%potentialprefac(at%iatype(icenter)), lin%confpotorder, &
!!     !!     lzd%llr(ilr)%nsi1, lzd%llr(ilr)%nsi2, lzd%llr(ilr)%nsi3,  &
!!     !!     lzd%llr(ilr)%bounds%ibyyzz_r) !optional
!!     if(lzd%llr(ilr)%geocode == 'F')then
!!        call position_operators(lzd%Glr%d%n1i,lzd%Glr%d%n2i,lzd%Glr%d%n3i, &
!!                             lzd%llr(ilr)%d%n1i, lzd%llr(ilr)%d%n2i, lzd%llr(ilr)%d%n3i, &
!!                             lzd%llr(ilr)%d%n1i, lzd%llr(ilr)%d%n2i, lzd%llr(ilr)%d%n3i, &
!!                             ishift, lzd%llr(ilr)%d%n2, lzd%llr(ilr)%d%n3, orbs%nspinor, &
!!                             psir, order, psirx, psiry, psirz, &
!!                             confdatarr(iorb), lzd%llr(ilr)%bounds%ibyyzz_r) !optional
!!     else
!!        call position_operators(lzd%Glr%d%n1i,lzd%Glr%d%n2i,lzd%Glr%d%n3i, &
!!                             lzd%llr(ilr)%d%n1i, lzd%llr(ilr)%d%n2i, lzd%llr(ilr)%d%n3i, &
!!                             lzd%llr(ilr)%d%n1i, lzd%llr(ilr)%d%n2i, lzd%llr(ilr)%d%n3i, &
!!                             ishift, lzd%llr(ilr)%d%n2, lzd%llr(ilr)%d%n3, orbs%nspinor, &
!!                             psir, order, psirx, psiry, psirz, &
!!                             confdatarr(iorb)) !optional
!!     end if
!!
!!     call isf_to_daub(lzd%llr(ilr), work_sr, psirx, xpsi(1+oidx))
!!     call isf_to_daub(lzd%llr(ilr), work_sr, psiry, ypsi(1+oidx))
!!     call isf_to_daub(lzd%llr(ilr), work_sr, psirz, zpsi(1+oidx))
!!
!!     !!call dcopy(Lzd%Llr(ilr)%d%n1i*Lzd%Llr(ilr)%d%n2i*Lzd%Llr(ilr)%d%n3i*orbs%nspinor, psir(1,1), 1, vpsir(1,1), 1)
!!     !!call position_operator(iproc, lzd%llr(ilr)%d%n1,lzd%llr(ilr)%d%n2,lzd%llr(ilr)%d%n3,1,1,1,0,orbs%nspinor, vpsir, &
!!     !!     hxh, hyh, hzh, confdatarr(iorb)%ioffset, 'x', &
!!     !!     lzd%llr(ilr)%bounds%ibyyzz_r) !optional
!!     !!call isf_to_daub(lzd%llr(ilr), work_sr, vpsir, xpsi(1+oidx))
!!
!!     !!call dcopy(Lzd%Llr(ilr)%d%n1i*Lzd%Llr(ilr)%d%n2i*Lzd%Llr(ilr)%d%n3i*orbs%nspinor, psir(1,1), 1, vpsir(1,1), 1)
!!     !!call position_operator(iproc, lzd%llr(ilr)%d%n1,lzd%llr(ilr)%d%n2,lzd%llr(ilr)%d%n3,1,1,1,0,orbs%nspinor, vpsir, &
!!     !!     hxh, hyh, hzh, confdatarr(iorb)%ioffset, 'y', &
!!     !!     lzd%llr(ilr)%bounds%ibyyzz_r) !optional
!!     !!call isf_to_daub(lzd%llr(ilr), work_sr, vpsir, ypsi(1+oidx))
!!
!!     !!call dcopy(Lzd%Llr(ilr)%d%n1i*Lzd%Llr(ilr)%d%n2i*Lzd%Llr(ilr)%d%n3i*orbs%nspinor, psir(1,1), 1, vpsir(1,1), 1)
!!     !!call position_operator(iproc, lzd%llr(ilr)%d%n1,lzd%llr(ilr)%d%n2,lzd%llr(ilr)%d%n3,1,1,1,0,orbs%nspinor, vpsir, &
!!     !!     hxh, hyh, hzh, confdatarr(iorb)%ioffset, 'z', &
!!     !!     lzd%llr(ilr)%bounds%ibyyzz_r) !optional
!!     !!call isf_to_daub(lzd%llr(ilr), work_sr, vpsir, zpsi(1+oidx))
!!
!!     !!iall=(Lzd%Llr(ilr)%wfd%nvctr_c+7*Lzd%Llr(ilr)%wfd%nvctr_f)*orbs%nspinor
!!     !!write(*,'(a,i5,es16.4)') 'iorb, ddot x', iorb, ddot(iall, xpsi(1+oidx), 1, psi(1+oidx), 1)
!!     !!write(*,'(a,i5,es16.4)') 'iorb, ddot y', iorb, ddot(iall, ypsi(1+oidx), 1, psi(1+oidx), 1)
!!     !!write(*,'(a,i5,es16.4)') 'iorb, ddot z', iorb, ddot(iall, zpsi(1+oidx), 1, psi(1+oidx), 1)
!!
!!
!!     i_all=-product(shape(psir))*kind(psir)
!!     deallocate(psir,stat=i_stat)
!!     call memocc(i_stat,i_all,'psir',subname)
!!
!!     i_all=-product(shape(psirx))*kind(psirx)
!!     deallocate(psirx,stat=i_stat)
!!     call memocc(i_stat,i_all,'psirx',subname)
!!
!!     i_all=-product(shape(psiry))*kind(psiry)
!!     deallocate(psiry,stat=i_stat)
!!     call memocc(i_stat,i_all,'psiry',subname)
!!
!!     i_all=-product(shape(psirz))*kind(psirz)
!!     deallocate(psirz,stat=i_stat)
!!     call memocc(i_stat,i_all,'psirz',subname)
!!
!!     call deallocate_work_arrays_sumrho(work_sr)
!!
!!     oidx = oidx + (Lzd%Llr(ilr)%wfd%nvctr_c+7*Lzd%Llr(ilr)%wfd%nvctr_f)*orbs%nspinor
!!
!!  enddo
!!
!!
!!end subroutine apply_position_operators





!!subroutine position_operators(Gn1i,Gn2i,Gn3i,n1i,n2i,n3i,n1ip,n2ip,n3ip,ishift,n2,n3,nspinor,psir,order,&
!!     psirx, psiry, psirz, &
!!     confdata,ibyyzz_r) !optional
!!  use module_base
!!  use module_types
!!  implicit none
!!  integer, intent(in) :: Gn1i,Gn2i,Gn3i,n1i,n2i,n3i,n1ip,n2ip,n3ip,n2,n3,nspinor,order
!!  integer, dimension(3), intent(in) :: ishift !<offset of potential box in wfn box coords.
!!  real(wp), dimension(n1i,n2i,n3i,nspinor), intent(in) :: psir !< real-space wfn in lr
!!  real(wp), dimension(n1i,n2i,n3i,nspinor), intent(out) :: psirx, psiry, psirz !< x,y,z operator applied to real-space wfn in lr
!!  type(confpot_data), intent(in), optional :: confdata !< data for the confining potential
!!  integer, dimension(2,-14:2*n2+16,-14:2*n3+16), intent(in), optional :: ibyyzz_r !< bounds in lr
!!  !local variables
!!  integer :: ii1,ii2,ii3,i1,i2,i3,ispinor,i1s,i1e,i2s,i2e,i3s,i3e,i1st,i1et
!!  real(wp) :: tt11,tt22,tt33,tt44,tt13,tt14,tt23,tt24,tt31,tt32,tt41,tt42,tt
!!  real(wp) :: psir1,psir2,psir3,psir4,pot1,pot2,pot3,pot4
!!  real(wp):: ttx, tty, ttz, potx, poty, potz
!!
!!  write(*,*) 'in position_operators'
!!
!!
!!  !loop on wavefunction
!!  !calculate the limits in all the directions
!!  !regions in which both the potential and wavefunctions are defined
!!  i3s=max(1,ishift(3)+1)
!!  i3e=min(n3i,n3ip+ishift(3))
!!  i2s=max(1,ishift(2)+1)
!!  i2e=min(n2i,n2ip+ishift(2))
!!  i1s=max(1,ishift(1)+1)
!!  i1e=min(n1i,n1ip+ishift(1))
!!
!!
!!  !$omp parallel default(none)&
!!  !$omp shared(psir,psirx,psiry,psirz,n1i,n2i,n3i,n1ip,n2ip,n3ip,n2,n3,ibyyzz_r,nspinor)&
!!  !$omp shared(i1s,i1e,i2s,i2e,i3s,i3e,ishift,confdata,order,Gn1i,Gn2i,Gn3i)&
!!  !$omp private(ispinor,i1,i2,i3,i1st,i1et)&
!!  !$omp private(tt11,tt22,tt33,tt44,tt13,tt14,tt23,tt24,tt31,tt32,tt41,tt42,tt)&
!!  !$omp private(psir1,psir2,psir3,psir4,pot1,pot2,pot3,pot4,ttx,tty,ttz,potx,poty,potz)
!!
!!!!$  !$omp parallel default(private)&
!!!!$  !$omp shared(pot,psir,n1i,n2i,n3i,n1ip,n2ip,n3ip,n2,n3,ibyyzz_r,nspinor)&
!!!!$  !$omp shared(i1s,i1e,i2s,i2e,i3s,i3e,ishift)
!!  !case without bounds
!!
!!
!!  !put to zero the external part of psir if the potential is more little than the wavefunction
!!  !first part of the array
!!  do ispinor=1,nspinor
!!     !$omp do 
!!     do i3=1,i3s-1
!!        do i2=1,n2i
!!           do i1=1,n1i
!!             psirx(i1,i2,i3,ispinor)=0.0_wp 
!!             psiry(i1,i2,i3,ispinor)=0.0_wp 
!!             psirz(i1,i2,i3,ispinor)=0.0_wp 
!!           end do
!!        end do
!!     end do
!!     !$omp end do
!!  end do
!!
!!  !central part of the array
!!  do ispinor=1,nspinor
!!     !$omp do 
!!     do i3=i3s,i3e
!!
!!        !first part
!!        do i2=1,i2s-1
!!           do i1=1,n1i
!!              psirx(i1,i2,i3,ispinor)=0.0_wp 
!!              psiry(i1,i2,i3,ispinor)=0.0_wp 
!!              psirz(i1,i2,i3,ispinor)=0.0_wp 
!!           end do
!!        end do
!!        !central part
!!        do i2=i2s,i2e
!!           do i1=1,i1s-1
!!              psirx(i1,i2,i3,ispinor)=0.0_wp 
!!              psiry(i1,i2,i3,ispinor)=0.0_wp 
!!              psirz(i1,i2,i3,ispinor)=0.0_wp 
!!           end do
!!           do i1=i1e+1,n1i
!!              psirx(i1,i2,i3,ispinor)=0.0_wp 
!!              psiry(i1,i2,i3,ispinor)=0.0_wp 
!!              psirz(i1,i2,i3,ispinor)=0.0_wp 
!!           end do
!!        end do
!!        !last part
!!        do i2=i2e+1,n2i
!!           do i1=1,n1i
!!              psirx(i1,i2,i3,ispinor)=0.0_wp 
!!              psiry(i1,i2,i3,ispinor)=0.0_wp 
!!              psirz(i1,i2,i3,ispinor)=0.0_wp 
!!           end do
!!        end do
!!
!!     end do
!!     !$omp end do
!!  end do
!!
!!
!!  !last part of the array
!!  do ispinor=1,nspinor
!!     !$omp do 
!!     do i3=i3e+1,n3i
!!        do i2=1,n2i
!!           do i1=1,n1i
!!              psirx(i1,i2,i3,ispinor)=0.0_wp 
!!              psiry(i1,i2,i3,ispinor)=0.0_wp 
!!              psirz(i1,i2,i3,ispinor)=0.0_wp 
!!           end do
!!        end do
!!     end do
!!     !$omp end do
!!  end do
!!
!!
!!  !important part of the array
!!  if (nspinor==4) then
!!      stop 'not yet implemented for nspinor==4!'
!!     !!!$omp do
!!     !!do i3=i3s,i3e
!!     !!   do i2=i2s,i2e
!!     !!      !thanks to the optional argument the conditional is done at compile time
!!     !!      if (present(ibyyzz_r)) then
!!     !!         i1st=max(i1s,ibyyzz_r(1,i2-15,i3-15)+1) !in bounds coordinates
!!     !!         i1et=min(i1e,ibyyzz_r(2,i2-15,i3-15)+1) !in bounds coordinates
!!     !!      else
!!     !!         i1st=i1s
!!     !!         i1et=i1e
!!     !!      end if
!!     !!      !no need of setting up to zero values outside wavefunction bounds
!!     !!      do i1=i1st,i1et
!!     !!         !wavefunctions
!!     !!         psir1=psir(i1,i2,i3,1)
!!     !!         psir2=psir(i1,i2,i3,2)
!!     !!         psir3=psir(i1,i2,i3,3)
!!     !!         psir4=psir(i1,i2,i3,4)
!!     !!         !potentials + confining term
!!     !!         pot1=pot(i1-ishift(1),i2-ishift(2),i3-ishift(3),1)+cp(i1,i2,i3)
!!     !!         pot2=pot(i1-ishift(1),i2-ishift(2),i3-ishift(3),2)+cp(i1,i2,i3)
!!     !!         pot3=pot(i1-ishift(1),i2-ishift(2),i3-ishift(3),3)+cp(i1,i2,i3)
!!     !!         pot4=pot(i1-ishift(1),i2-ishift(2),i3-ishift(3),4)+cp(i1,i2,i3)
!!
!!     !!         !diagonal terms
!!     !!         tt11=pot1*psir1 !p1
!!     !!         tt22=pot1*psir2 !p2
!!     !!         tt33=pot4*psir3 !p3
!!     !!         tt44=pot4*psir4 !p4
!!     !!         !Rab*Rb
!!     !!         tt13=pot2*psir3 !p1
!!     !!         !Iab*Ib
!!     !!         tt14=pot3*psir4 !p1
!!     !!         !Rab*Ib
!!     !!         tt23=pot2*psir4 !p2
!!     !!         !Iab*Rb
!!     !!         tt24=pot3*psir3 !p2
!!     !!         !Rab*Ra
!!     !!         tt31=pot2*psir1 !p3
!!     !!         !Iab*Ia
!!     !!         tt32=pot3*psir2 !p3
!!     !!         !Rab*Ia
!!     !!         tt41=pot2*psir2 !p4
!!     !!         !Iab*Ra
!!     !!         tt42=pot3*psir1 !p4
!!
!!     !!         !value of the potential energy
!!     !!         epot_p=epot_p+tt11*psir1+tt22*psir2+tt33*psir3+tt44*psir4+&
!!     !!              2.0_gp*tt31*psir3-2.0_gp*tt42*psir4+2.0_gp*tt41*psir4+2.0_gp*tt32*psir3
!!
!!     !!         !wavefunction update
!!     !!         !p1=h1p1+h2p3-h3p4
!!     !!         !p2=h1p2+h2p4+h3p3
!!     !!         !p3=h2p1+h3p2+h4p3
!!     !!         !p4=h2p2-h3p1+h4p4
!!     !!         psir(i1,i2,i3,1)=tt11+tt13-tt14
!!     !!         psir(i1,i2,i3,2)=tt22+tt23+tt24
!!     !!         psir(i1,i2,i3,3)=tt33+tt31+tt32
!!     !!         psir(i1,i2,i3,4)=tt44+tt41-tt42
!!     !!      end do
!!     !!   end do
!!     !!end do
!!     !!!$omp end do
!!
!!  else !case with nspinor /=4
!!     write(*,'(a,3i9)') 'confdata%ioffset(1), confdata%ioffset(2), confdata%ioffset(3)', &
!!                confdata%ioffset(1), confdata%ioffset(2), confdata%ioffset(3)
!!     do ispinor=1,nspinor
!!        !$omp do
!!        do ii3=i3s,i3e
!!           i3=mod(ii3+confdata%ioffset(3)-1,Gn3i)+1
!!           do ii2=i2s,i2e
!!              i2=mod(ii2+confdata%ioffset(2)-1,Gn2i)+1
!!              !thanks to the optional argument the conditional is done at compile time
!!              if (present(ibyyzz_r)) then
!!                 i1st=max(i1s,ibyyzz_r(1,ii2-15,ii3-15)+1) !in bounds coordinates
!!                 i1et=min(i1e,ibyyzz_r(2,ii2-15,ii3-15)+1) !in bounds coordinates
!!              else
!!                 i1st=i1s
!!                 i1et=i1e
!!              end if
!!              !no need of setting up to zero values outside wavefunction bounds
!!              do ii1=i1st,i1et
!!                 i1=mod(ii1+confdata%ioffset(1)-1,Gn1i)+1
!!                 psir1=psir(ii1,ii2,ii3,ispinor)
!!                 !the local potential is always real (npot=1) + confining term
!!                 !!pot1=pot(i1-ishift(1),i2-ishift(2),i3-ishift(3),1)+cp(i1,i2,i3)
!!                 potx=(confdata%hh(1)*real(i1,wp))**order
!!                 poty=(confdata%hh(2)*real(i2,wp))**order
!!                 potz=(confdata%hh(3)*real(i3,wp))**order
!!
!!                 ttx=potx*psir1
!!                 tty=poty*psir1
!!                 ttz=potz*psir1
!!
!!                 psirx(ii1,ii2,ii3,ispinor)=ttx
!!                 psiry(ii1,ii2,ii3,ispinor)=tty
!!                 psirz(ii1,ii2,ii3,ispinor)=ttz
!!              end do
!!           end do
!!        end do
!!        !$omp end do
!!     end do
!!  end if
!!  
!!  
!!  !$omp end parallel
!!
!!
!!END SUBROUTINE position_operators







!!subroutine commutator(norb, A, B, res)
!!implicit none
!!
!!! Calling arguments
!!integer,intent(in):: norb
!!real(8),dimension(norb,norb),intent(in):: A, B
!!real(8),dimension(norb,norb),intent(out):: res
!!
!!! Local variables
!!real(8),dimension(norb,norb):: AB, BA
!!
!!call dgemm('n', 'n', norb, norb, norb, 1.d0, A, norb, B, norb, 0.d0, AB, norb)
!!call dgemm('n', 'n', norb, norb, norb, 1.d0, B, norb, A, norb, 0.d0, BA, norb)
!!res=AB-BA
!!
!!end subroutine commutator


!!subroutine apply_r_operators(iproc, nproc, orbs, lzd, hx, hy, hz, confdatarr, psi, order, vpsi)
!!use module_base
!!use module_types
!!use module_interfaces, except_this_one => apply_r_operators
!!implicit none
!!
!!! Calling arguments
!!integer,intent(in):: iproc, nproc, order
!!type(orbitals_data),intent(in):: orbs
!!type(local_zone_descriptors),intent(in):: lzd
!!real(8),intent(in):: hx, hy, hz
!!type(confpot_data),dimension(orbs%norbp),intent(in):: confdatarr
!!real(8),dimension(max(orbs%npsidim_orbs,orbs%npsidim_comp)),intent(in):: psi
!!real(8),dimension(max(orbs%npsidim_orbs,orbs%npsidim_comp)),intent(out):: vpsi
!!
!!! Local variables
!!integer:: oidx, iorb, ilr, npot, icenter, i_stat, i_all, ist_c, ist_f, ist, iiorb, iall, ierr
!!real(8):: hxh, hyh, hzh, ddot, tt, t1, t2, time
!!real(8),dimension(:,:),allocatable:: psir, psirx, psiry, psirz
!!type(workarr_sumrho):: work_sr
!!real(8),dimension(0:3),parameter:: scal=1.d0
!!real(8),dimension(:,:,:),allocatable:: ypsitemp_c
!!real(8),dimension(:,:,:,:),allocatable:: ypsitemp_f
!!character(len=*),parameter:: subname='apply_orbitaldependent_potential'
!!integer, dimension(3) :: ishift !temporary variable in view of wavefunction creation
!!
!!  ishift=(/0,0,0/)
!!
!!  !xpsi=0.d0
!!  !ypsi=0.d0
!!  call to_zero(max(orbs%npsidim_orbs,orbs%npsidim_comp), vpsi(1))
!!  oidx = 0
!!  do iorb=1,orbs%norbp
!!     ilr = orbs%inwhichlocreg(iorb+orbs%isorb)
!!
!!     iiorb=orbs%isorb+iorb
!!     !!write(*,'(a,4i8,4x,3i6)') 'iproc, iorb, iiorb, ilr, confdatarr(iorb)%ioffset(:)', &
!!     !!    iproc, iorb, iiorb, ilr, confdatarr(iorb)%ioffset(:)
!!     !!write(*,'(a,3i8,6i6)') 'iproc, iiorb, ilr, is1, ie1, is2, ie2, is3, ie3', &
!!     !!    1, lzd%llr(ilr)%d%n1i, 1, lzd%llr(ilr)%d%n2i, 1, lzd%llr(ilr)%d%n3i
!!  
!!     !initialise the work arrays
!!     call initialize_work_arrays_sumrho(lzd%llr(ilr), work_sr)
!!
!!     ! Wavefunction in real space
!!     allocate(psir(Lzd%Llr(ilr)%d%n1i*Lzd%Llr(ilr)%d%n2i*Lzd%Llr(ilr)%d%n3i,orbs%nspinor+ndebug),stat=i_stat)
!!     call memocc(i_stat,psir,'psir',subname)
!!     call razero(Lzd%Llr(ilr)%d%n1i*Lzd%Llr(ilr)%d%n2i*Lzd%Llr(ilr)%d%n3i*orbs%nspinor,psir)
!!
!!     !!allocate(psirx(Lzd%Llr(ilr)%d%n1i*Lzd%Llr(ilr)%d%n2i*Lzd%Llr(ilr)%d%n3i,orbs%nspinor+ndebug),stat=i_stat)
!!     !!call memocc(i_stat,psirx,'psirx',subname)
!!     !!call razero(Lzd%Llr(ilr)%d%n1i*Lzd%Llr(ilr)%d%n2i*Lzd%Llr(ilr)%d%n3i*orbs%nspinor,psirx)
!!
!!     !!allocate(psiry(Lzd%Llr(ilr)%d%n1i*Lzd%Llr(ilr)%d%n2i*Lzd%Llr(ilr)%d%n3i,orbs%nspinor+ndebug),stat=i_stat)
!!     !!call memocc(i_stat,psiry,'psiry',subname)
!!     !!call razero(Lzd%Llr(ilr)%d%n1i*Lzd%Llr(ilr)%d%n2i*Lzd%Llr(ilr)%d%n3i*orbs%nspinor,psiry)
!!
!!     !!allocate(psirz(Lzd%Llr(ilr)%d%n1i*Lzd%Llr(ilr)%d%n2i*Lzd%Llr(ilr)%d%n3i,orbs%nspinor+ndebug),stat=i_stat)
!!     !!call memocc(i_stat,psirz,'psirz',subname)
!!     !!call razero(Lzd%Llr(ilr)%d%n1i*Lzd%Llr(ilr)%d%n2i*Lzd%Llr(ilr)%d%n3i*orbs%nspinor,psirz)
!!
!!     !transform the wavefunction in Daubechies basis to the wavefunction in ISF basis
!!     !the psir wavefunction is given in the spinorial form
!!
!!     !psi(1+oidx+lzd%llr(ilr)%wfd%nvctr_c:1+oidx+lzd%llr(ilr)%wfd%nvctr_c+7*lzd%llr(ilr)%wfd%nvctr_f-1)=0.d0
!!
!!     call daub_to_isf(lzd%llr(ilr), work_sr, psi(1+oidx), psir)
!!     !!!write(*,*) 'WARNING DEBUG in r_operator'
!!     !!!psir=1.d0/sqrt(dble(Lzd%Llr(ilr)%d%n1i*Lzd%Llr(ilr)%d%n2i*Lzd%Llr(ilr)%d%n3i))
!!     !apply the potential to the psir wavefunction and calculate potential energy
!!     hxh=.5d0*hx
!!     hyh=.5d0*hy
!!     hzh=.5d0*hz
!!     !icenter=confinementCenter(iorb)
!!     !components of the potential
!!     npot=orbs%nspinor
!!     if (orbs%nspinor == 2) npot=1
!!
!!     !!call apply_confinement(iproc, lzd%llr(ilr)%d%n1,lzd%llr(ilr)%d%n2,lzd%llr(ilr)%d%n3,1,1,1,0,orbs%nspinor, psir, &
!!     !!     rxyz(1,icenter), hxh, hyh, hzh, lin%potentialprefac(at%iatype(icenter)), lin%confpotorder, &
!!     !!     lzd%llr(ilr)%nsi1, lzd%llr(ilr)%nsi2, lzd%llr(ilr)%nsi3,  &
!!     !!     lzd%llr(ilr)%bounds%ibyyzz_r) !optional
!!     if(lzd%llr(ilr)%geocode == 'F') then
!!        call r_operator(lzd%Glr%d%n1i, lzd%Glr%d%n2i, lzd%Glr%d%n3i, &
!!                        lzd%llr(ilr)%d%n1i, lzd%llr(ilr)%d%n2i, lzd%llr(ilr)%d%n3i, &
!!                        lzd%llr(ilr)%d%n1i, lzd%llr(ilr)%d%n2i, lzd%llr(ilr)%d%n3i, &
!!                        ishift, lzd%llr(ilr)%d%n2, lzd%llr(ilr)%d%n3, orbs%nspinor, &
!!                        psir, order, &
!!                        confdatarr(iorb), lzd%llr(ilr)%bounds%ibyyzz_r) !optional
!!     else
!!        call r_operator(lzd%Glr%d%n1i, lzd%Glr%d%n2i, lzd%Glr%d%n3i, &
!!                        lzd%llr(ilr)%d%n1i, lzd%llr(ilr)%d%n2i, lzd%llr(ilr)%d%n3i, &
!!                        lzd%llr(ilr)%d%n1i, lzd%llr(ilr)%d%n2i, lzd%llr(ilr)%d%n3i, &
!!                        ishift, lzd%llr(ilr)%d%n2, lzd%llr(ilr)%d%n3, orbs%nspinor, &
!!                        psir, order, &
!!                        confdatarr(iorb)) !optional
!!     end if
!!
!!     call isf_to_daub(lzd%llr(ilr), work_sr, psir, vpsi(1+oidx))
!!     !!call isf_to_daub(lzd%llr(ilr), work_sr, psiry, ypsi(1+oidx))
!!     !!call isf_to_daub(lzd%llr(ilr), work_sr, psirz, zpsi(1+oidx))
!!
!!     !!call dcopy(Lzd%Llr(ilr)%d%n1i*Lzd%Llr(ilr)%d%n2i*Lzd%Llr(ilr)%d%n3i*orbs%nspinor, psir(1,1), 1, vpsir(1,1), 1)
!!     !!call position_operator(iproc, lzd%llr(ilr)%d%n1,lzd%llr(ilr)%d%n2,lzd%llr(ilr)%d%n3,1,1,1,0,orbs%nspinor, vpsir, &
!!     !!     hxh, hyh, hzh, confdatarr(iorb)%ioffset, 'x', &
!!     !!     lzd%llr(ilr)%bounds%ibyyzz_r) !optional
!!     !!call isf_to_daub(lzd%llr(ilr), work_sr, vpsir, xpsi(1+oidx))
!!
!!     !!call dcopy(Lzd%Llr(ilr)%d%n1i*Lzd%Llr(ilr)%d%n2i*Lzd%Llr(ilr)%d%n3i*orbs%nspinor, psir(1,1), 1, vpsir(1,1), 1)
!!     !!call position_operator(iproc, lzd%llr(ilr)%d%n1,lzd%llr(ilr)%d%n2,lzd%llr(ilr)%d%n3,1,1,1,0,orbs%nspinor, vpsir, &
!!     !!     hxh, hyh, hzh, confdatarr(iorb)%ioffset, 'y', &
!!     !!     lzd%llr(ilr)%bounds%ibyyzz_r) !optional
!!     !!call isf_to_daub(lzd%llr(ilr), work_sr, vpsir, ypsi(1+oidx))
!!
!!     !!call dcopy(Lzd%Llr(ilr)%d%n1i*Lzd%Llr(ilr)%d%n2i*Lzd%Llr(ilr)%d%n3i*orbs%nspinor, psir(1,1), 1, vpsir(1,1), 1)
!!     !!call position_operator(iproc, lzd%llr(ilr)%d%n1,lzd%llr(ilr)%d%n2,lzd%llr(ilr)%d%n3,1,1,1,0,orbs%nspinor, vpsir, &
!!     !!     hxh, hyh, hzh, confdatarr(iorb)%ioffset, 'z', &
!!     !!     lzd%llr(ilr)%bounds%ibyyzz_r) !optional
!!     !!call isf_to_daub(lzd%llr(ilr), work_sr, vpsir, zpsi(1+oidx))
!!
!!     !!iall=(Lzd%Llr(ilr)%wfd%nvctr_c+7*Lzd%Llr(ilr)%wfd%nvctr_f)*orbs%nspinor
!!     !!write(*,'(a,i5,es16.4)') 'iorb, ddot x', iorb, ddot(iall, xpsi(1+oidx), 1, psi(1+oidx), 1)
!!     !!write(*,'(a,i5,es16.4)') 'iorb, ddot y', iorb, ddot(iall, ypsi(1+oidx), 1, psi(1+oidx), 1)
!!     !!write(*,'(a,i5,es16.4)') 'iorb, ddot z', iorb, ddot(iall, zpsi(1+oidx), 1, psi(1+oidx), 1)
!!
!!
!!     i_all=-product(shape(psir))*kind(psir)
!!     deallocate(psir,stat=i_stat)
!!     call memocc(i_stat,i_all,'psir',subname)
!!
!!     !!i_all=-product(shape(psirx))*kind(psirx)
!!     !!deallocate(psirx,stat=i_stat)
!!     !!call memocc(i_stat,i_all,'psirx',subname)
!!
!!     !!i_all=-product(shape(psiry))*kind(psiry)
!!     !!deallocate(psiry,stat=i_stat)
!!     !!call memocc(i_stat,i_all,'psiry',subname)
!!
!!     !!i_all=-product(shape(psirz))*kind(psirz)
!!     !!deallocate(psirz,stat=i_stat)
!!     !!call memocc(i_stat,i_all,'psirz',subname)
!!
!!     call deallocate_work_arrays_sumrho(work_sr)
!!
!!     oidx = oidx + (Lzd%Llr(ilr)%wfd%nvctr_c+7*Lzd%Llr(ilr)%wfd%nvctr_f)*orbs%nspinor
!!
!!  enddo
!!
!!
!!end subroutine apply_r_operators








!!subroutine r_operator(Gn1i,Gn2i,Gn3i,n1i,n2i,n3i,n1ip,n2ip,n3ip,ishift,n2,n3,nspinor,psir,order,&
!!     confdata,ibyyzz_r) !optional
!!  use module_base
!!  use module_types
!!  implicit none
!!  integer, intent(in) :: Gn1i,Gn2i,Gn3i,n1i,n2i,n3i,n1ip,n2ip,n3ip,n2,n3,nspinor,order
!!  integer, dimension(3), intent(in) :: ishift !<offset of potential box in wfn box coords.
!!  real(wp), dimension(n1i,n2i,n3i,nspinor), intent(inout) :: psir !< real-space wfn in lr
!!  type(confpot_data), intent(in), optional :: confdata !< data for the confining potential
!!  integer, dimension(2,-14:2*n2+16,-14:2*n3+16), intent(in), optional :: ibyyzz_r !< bounds in lr
!!  !local variables
!!  integer :: i1,i2,i3,ii1,ii2,ii3,ispinor,i1s,i1e,i2s,i2e,i3s,i3e,i1st,i1et
!!  real(wp) :: tt11,tt22,tt33,tt44,tt13,tt14,tt23,tt24,tt31,tt32,tt41,tt42,tt
!!  real(wp) :: psir1,psir2,psir3,psir4,pot1,pot2,pot3,pot4
!!  real(wp):: ttx, tty, ttz, potx, poty, potz
!!
!!  if(order/=1 .and. order/=2) stop 'wrong order'
!!
!!  !loop on wavefunction
!!  !calculate the limits in all the directions
!!  !regions in which both the potential and wavefunctions are defined
!!  i3s=max(1,ishift(3)+1)
!!  i3e=min(n3i,n3ip+ishift(3))
!!  i2s=max(1,ishift(2)+1)
!!  i2e=min(n2i,n2ip+ishift(2))
!!  i1s=max(1,ishift(1)+1)
!!  i1e=min(n1i,n1ip+ishift(1))
!!
!!
!!  !$omp parallel default(none)&
!!  !$omp shared(psir,n1i,n2i,n3i,n1ip,n2ip,n3ip,n2,n3,ibyyzz_r,nspinor)&
!!  !$omp shared(i1s,i1e,i2s,i2e,i3s,i3e,ishift,confdata,order,Gn1i,Gn2i,Gn3i)&
!!  !$omp private(ispinor,i1,i2,i3,i1st,i1et)&
!!  !$omp private(tt11,tt22,tt33,tt44,tt13,tt14,tt23,tt24,tt31,tt32,tt41,tt42,tt)&
!!  !$omp private(psir1,psir2,psir3,psir4,pot1,pot2,pot3,pot4,ttx,tty,ttz,potx,poty,potz)
!!
!!!!$  !$omp parallel default(private)&
!!!!$  !$omp shared(pot,psir,n1i,n2i,n3i,n1ip,n2ip,n3ip,n2,n3,ibyyzz_r,nspinor)&
!!!!$  !$omp shared(i1s,i1e,i2s,i2e,i3s,i3e,ishift)
!!  !case without bounds
!!
!!
!!  !put to zero the external part of psir if the potential is more little than the wavefunction
!!  !first part of the array
!!  do ispinor=1,nspinor
!!     !$omp do 
!!     do i3=1,i3s-1
!!        do i2=1,n2i
!!           do i1=1,n1i
!!             psir(i1,i2,i3,ispinor)=0.0_wp 
!!             !!psiry(i1,i2,i3,ispinor)=0.0_wp 
!!             !!psirz(i1,i2,i3,ispinor)=0.0_wp 
!!           end do
!!        end do
!!     end do
!!     !$omp end do
!!  end do
!!
!!  !central part of the array
!!  do ispinor=1,nspinor
!!     !$omp do 
!!     do i3=i3s,i3e
!!
!!        !first part
!!        do i2=1,i2s-1
!!           do i1=1,n1i
!!              psir(i1,i2,i3,ispinor)=0.0_wp 
!!              !!psiry(i1,i2,i3,ispinor)=0.0_wp 
!!              !!psirz(i1,i2,i3,ispinor)=0.0_wp 
!!           end do
!!        end do
!!        !central part
!!        do i2=i2s,i2e
!!           do i1=1,i1s-1
!!              psir(i1,i2,i3,ispinor)=0.0_wp 
!!              !!psiry(i1,i2,i3,ispinor)=0.0_wp 
!!              !!psirz(i1,i2,i3,ispinor)=0.0_wp 
!!           end do
!!           do i1=i1e+1,n1i
!!              psir(i1,i2,i3,ispinor)=0.0_wp 
!!              !!psiry(i1,i2,i3,ispinor)=0.0_wp 
!!              !!psirz(i1,i2,i3,ispinor)=0.0_wp 
!!           end do
!!        end do
!!        !last part
!!        do i2=i2e+1,n2i
!!           do i1=1,n1i
!!              psir(i1,i2,i3,ispinor)=0.0_wp 
!!              !!psiry(i1,i2,i3,ispinor)=0.0_wp 
!!              !!psirz(i1,i2,i3,ispinor)=0.0_wp 
!!           end do
!!        end do
!!
!!     end do
!!     !$omp end do
!!  end do
!!
!!
!!  !last part of the array
!!  do ispinor=1,nspinor
!!     !$omp do 
!!     do i3=i3e+1,n3i
!!        do i2=1,n2i
!!           do i1=1,n1i
!!              psir(i1,i2,i3,ispinor)=0.0_wp 
!!              !!psiry(i1,i2,i3,ispinor)=0.0_wp 
!!              !!psirz(i1,i2,i3,ispinor)=0.0_wp 
!!           end do
!!        end do
!!     end do
!!     !$omp end do
!!  end do
!!
!!
!!  !important part of the array
!!  if (nspinor==4) then
!!      stop 'not yet implemented for nspinor==4!'
!!     !!!$omp do
!!     !!do i3=i3s,i3e
!!     !!   do i2=i2s,i2e
!!     !!      !thanks to the optional argument the conditional is done at compile time
!!     !!      if (present(ibyyzz_r)) then
!!     !!         i1st=max(i1s,ibyyzz_r(1,i2-15,i3-15)+1) !in bounds coordinates
!!     !!         i1et=min(i1e,ibyyzz_r(2,i2-15,i3-15)+1) !in bounds coordinates
!!     !!      else
!!     !!         i1st=i1s
!!     !!         i1et=i1e
!!     !!      end if
!!     !!      !no need of setting up to zero values outside wavefunction bounds
!!     !!      do i1=i1st,i1et
!!     !!         !wavefunctions
!!     !!         psir1=psir(i1,i2,i3,1)
!!     !!         psir2=psir(i1,i2,i3,2)
!!     !!         psir3=psir(i1,i2,i3,3)
!!     !!         psir4=psir(i1,i2,i3,4)
!!     !!         !potentials + confining term
!!     !!         pot1=pot(i1-ishift(1),i2-ishift(2),i3-ishift(3),1)+cp(i1,i2,i3)
!!     !!         pot2=pot(i1-ishift(1),i2-ishift(2),i3-ishift(3),2)+cp(i1,i2,i3)
!!     !!         pot3=pot(i1-ishift(1),i2-ishift(2),i3-ishift(3),3)+cp(i1,i2,i3)
!!     !!         pot4=pot(i1-ishift(1),i2-ishift(2),i3-ishift(3),4)+cp(i1,i2,i3)
!!
!!     !!         !diagonal terms
!!     !!         tt11=pot1*psir1 !p1
!!     !!         tt22=pot1*psir2 !p2
!!     !!         tt33=pot4*psir3 !p3
!!     !!         tt44=pot4*psir4 !p4
!!     !!         !Rab*Rb
!!     !!         tt13=pot2*psir3 !p1
!!     !!         !Iab*Ib
!!     !!         tt14=pot3*psir4 !p1
!!     !!         !Rab*Ib
!!     !!         tt23=pot2*psir4 !p2
!!     !!         !Iab*Rb
!!     !!         tt24=pot3*psir3 !p2
!!     !!         !Rab*Ra
!!     !!         tt31=pot2*psir1 !p3
!!     !!         !Iab*Ia
!!     !!         tt32=pot3*psir2 !p3
!!     !!         !Rab*Ia
!!     !!         tt41=pot2*psir2 !p4
!!     !!         !Iab*Ra
!!     !!         tt42=pot3*psir1 !p4
!!
!!     !!         !value of the potential energy
!!     !!         epot_p=epot_p+tt11*psir1+tt22*psir2+tt33*psir3+tt44*psir4+&
!!     !!              2.0_gp*tt31*psir3-2.0_gp*tt42*psir4+2.0_gp*tt41*psir4+2.0_gp*tt32*psir3
!!
!!     !!         !wavefunction update
!!     !!         !p1=h1p1+h2p3-h3p4
!!     !!         !p2=h1p2+h2p4+h3p3
!!     !!         !p3=h2p1+h3p2+h4p3
!!     !!         !p4=h2p2-h3p1+h4p4
!!     !!         psir(i1,i2,i3,1)=tt11+tt13-tt14
!!     !!         psir(i1,i2,i3,2)=tt22+tt23+tt24
!!     !!         psir(i1,i2,i3,3)=tt33+tt31+tt32
!!     !!         psir(i1,i2,i3,4)=tt44+tt41-tt42
!!     !!      end do
!!     !!   end do
!!     !!end do
!!     !!!$omp end do
!!
!!  else !case with nspinor /=4
!!     do ispinor=1,nspinor
!!        !$omp do
!!        do ii3=i3s,i3e
!!           i3=mod(ii3+confdata%ioffset(3)-1,Gn3i)+1
!!           do ii2=i2s,i2e
!!              i2=mod(ii2+confdata%ioffset(2)-1,Gn2i)+1
!!              !thanks to the optional argument the conditional is done at compile time
!!              if (present(ibyyzz_r)) then
!!                 i1st=max(i1s,ibyyzz_r(1,ii2-15,ii3-15)+1) !in bounds coordinates
!!                 i1et=min(i1e,ibyyzz_r(2,ii2-15,ii3-15)+1) !in bounds coordinates
!!              else
!!                 i1st=i1s
!!                 i1et=i1e
!!              end if
!!              !no need of setting up to zero values outside wavefunction bounds
!!              do ii1=i1st,i1et
!!                 i1=mod(ii1+confdata%ioffset(1)-1,Gn1i)+1
!!                 psir1=psir(ii1,ii2,ii3,ispinor)
!!                 !the local potential is always real (npot=1) + confining term
!!                 !!pot1=pot(i1-ishift(1),i2-ishift(2),i3-ishift(3),1)+cp(i1,i2,i3)
!!                 ttx=(confdata%hh(1)*real(i1,wp))**2
!!                 tty=(confdata%hh(2)*real(i2,wp))**2
!!                 ttz=(confdata%hh(3)*real(i3,wp))**2
!!
!!                 tt = ttx+tty+ttz
!!
!!                 if(order==1) then
!!                     tt=sqrt(tt)
!!                 end if
!!
!!                 psir(ii1,ii2,ii3,ispinor)=tt*psir1
!!              end do
!!           end do
!!        end do
!!        !$omp end do
!!     end do
!!  end if
!!  
!!  
!!  !$omp end parallel
!!
!!
!!END SUBROUTINE r_operator




!!subroutine update_confdatarr(lzd, orbs, locregCenter, confdatarr)
!!  use module_base
!!  use module_types
!!  implicit none
!!  
!!  ! Calling arguments
!!  type(local_zone_descriptors),intent(in):: lzd
!!  type(orbitals_data),intent(in):: orbs
!!  real(8),dimension(3,lzd%nlr),intent(in):: locregCenter
!!  type(confpot_data),dimension(orbs%norbp),intent(inout):: confdatarr
!!  
!!  ! Local variables
!!  integer:: iorb, iiorb, ilr, icenter, nl1, nl2, nl3
!!  
!!  ! Update confdatarr...
!!  do iorb=1,orbs%norbp
!!     iiorb=orbs%isorb+iorb
!!     ilr=orbs%inWhichlocreg(iiorb)
!!     icenter=orbs%inwhichlocreg(iiorb)
!!     !confdatarr(iorb)%potorder=lin%confpotorder
!!     !confdatarr(iorb)%prefac=lin%potentialprefac(at%iatype(icenter))
!!     !confdatarr(iorb)%hh(1)=.5_gp*hx
!!     !confdatarr(iorb)%hh(2)=.5_gp*hy
!!     !confdatarr(iorb)%hh(3)=.5_gp*hz
!!     confdatarr(iorb)%rxyzConf(1:3)=locregCenter(1:3,icenter)
!!     call my_geocode_buffers(lzd%Llr(ilr)%geocode,nl1,nl2,nl3)
!!     confdatarr(iorb)%ioffset(1)=lzd%llr(ilr)%nsi1-nl1-1
!!     confdatarr(iorb)%ioffset(2)=lzd%llr(ilr)%nsi2-nl2-1
!!     confdatarr(iorb)%ioffset(3)=lzd%llr(ilr)%nsi3-nl3-1
!!  end do
!!
!!end subroutine update_confdatarr


!!function dfactorial(n)
!!implicit none
!!
!!! Calling arguments
!!integer,intent(in):: n
!!real(8):: dfactorial
!!
!!! Local variables
!!integer:: i
!!
!!  dfactorial=1.d0
!!  do i=1,n
!!      dfactorial=dfactorial*dble(i)
!!  end do
!!
!!end function dfactorial



!!subroutine apply_orbitaldependent_potential(iproc, nproc, at, orbs, lzd, rxyz, confdatarr, hx, &
!!           psi, centralLocreg, vpsi)
!!use module_base
!!use module_types
!!use module_interfaces, except_this_one => apply_orbitaldependent_potential
!!implicit none
!!
!!! Calling arguments
!!integer,intent(in):: iproc, nproc, centralLocreg
!!type(atoms_data),intent(in):: at
!!type(orbitals_data),intent(in):: orbs
!!type(local_zone_descriptors),intent(in):: lzd
!!real(8),dimension(3,at%nat),intent(in):: rxyz
!!type(confpot_data),dimension(orbs%norbp),intent(in):: confdatarr
!!real(8),intent(in):: hx
!!!real(8),dimension(lzd%lpsidimtot),intent(in):: psi  !!!! ATENTION, intent should be in !
!!!real(8),dimension(lzd%lpsidimtot),intent(inout):: psi
!!real(8),dimension(max(orbs%npsidim_orbs,orbs%npsidim_comp)),intent(inout):: psi
!!real(8),dimension(max(orbs%npsidim_orbs,orbs%npsidim_comp)),intent(out):: vpsi
!!
!!! Local variables
!!integer:: oidx, iorb, ilr, npot, icenter, i_stat, i_all, ist_c, ist_f, ist, iiorb, iall, ierr
!!real(8):: hxh, hyh, hzh, ddot, tt, t1, t2, time
!!real(8),dimension(:,:),allocatable:: psir, vpsir
!!type(workarr_precond) :: work
!!type(workarrays_quartic_convolutions):: work_conv
!!real(8),dimension(0:3),parameter:: scal=1.d0
!!real(8),dimension(:,:,:),allocatable:: ypsitemp_c
!!real(8),dimension(:,:,:,:),allocatable:: ypsitemp_f
!!character(len=*),parameter:: subname='apply_orbitaldependent_potential'
!!
!!
!!
!!  call to_zero(max(orbs%npsidim_orbs,orbs%npsidim_comp), vpsi(1))
!!  ist_c=1
!!  ist_f=1
!!  do iorb=1,orbs%norbp
!!      iiorb=iorb+orbs%isorb
!!      ilr = orbs%inwhichlocreg(iiorb)
!!      if(centralLocreg<0) then
!!          !icenter=lin%orbs%inWhichLocregp(iorb)
!!          icenter=orbs%inWhichLocreg(iiorb)
!!      else
!!          icenter=centralLocreg
!!      end if
!!      !components of the potential
!!      npot=orbs%nspinor
!!      if (orbs%nspinor == 2) npot=1
!!      ist_f=ist_f+lzd%llr(ilr)%wfd%nvctr_c
!!      call allocate_workarrays_quartic_convolutions(lzd%llr(ilr), subname, work_conv)
!!
!!      call uncompress_for_quartic_convolutions(lzd%llr(ilr)%d%n1, lzd%llr(ilr)%d%n2, lzd%llr(ilr)%d%n3, &
!!           lzd%llr(ilr)%d%nfl1, lzd%llr(ilr)%d%nfu1, &
!!           lzd%llr(ilr)%d%nfl2, lzd%llr(ilr)%d%nfu2, &
!!           lzd%llr(ilr)%d%nfl3, lzd%llr(ilr)%d%nfu3, &
!!           lzd%llr(ilr)%wfd%nseg_c, lzd%llr(ilr)%wfd%nvctr_c, &
!!           lzd%llr(ilr)%wfd%keygloc, lzd%llr(ilr)%wfd%keyvloc,  & 
!!           lzd%llr(ilr)%wfd%nseg_f, lzd%llr(ilr)%wfd%nvctr_f, &
!!           lzd%llr(ilr)%wfd%keygloc(1,lzd%llr(ilr)%wfd%nseg_c+min(1,lzd%llr(ilr)%wfd%nseg_f)), &
!!           lzd%llr(ilr)%wfd%keyvloc(lzd%llr(ilr)%wfd%nseg_c+min(1,lzd%llr(ilr)%wfd%nseg_f)),  & 
!!           scal, psi(ist_c), psi(ist_f), &
!!           work_conv)
!!
!!      if(confdatarr(iorb)%potorder==4) then
!!          call ConvolQuartic4(iproc,nproc,lzd%llr(ilr)%d%n1, lzd%llr(ilr)%d%n2, lzd%llr(ilr)%d%n3, &
!!               lzd%llr(ilr)%d%nfl1, lzd%llr(ilr)%d%nfu1, &
!!               lzd%llr(ilr)%d%nfl2, lzd%llr(ilr)%d%nfu2, &
!!               lzd%llr(ilr)%d%nfl3, lzd%llr(ilr)%d%nfu3, & 
!!               hx, lzd%llr(ilr)%ns1, lzd%llr(ilr)%ns2, lzd%llr(ilr)%ns3, &
!!               lzd%llr(ilr)%bounds%kb%ibyz_c, lzd%llr(ilr)%bounds%kb%ibxz_c, lzd%llr(ilr)%bounds%kb%ibxy_c, &
!!               lzd%llr(ilr)%bounds%kb%ibyz_f, lzd%llr(ilr)%bounds%kb%ibxz_f, lzd%llr(ilr)%bounds%kb%ibxy_f, &
!!               rxyz(1,ilr), confdatarr(iorb)%prefac, .false., 0.d0, &
!!               work_conv%xx_c, work_conv%xx_f1, work_conv%xx_f, &
!!               work_conv%xy_c, work_conv%xy_f2, work_conv%xy_f, &
!!               work_conv%xz_c, work_conv%xz_f4, work_conv%xz_f, &
!!               work_conv%y_c, work_conv%y_f)
!!      else if(confdatarr(iorb)%potorder==6) then
!!          call ConvolSextic(lzd%llr(ilr)%d%n1, lzd%llr(ilr)%d%n2, lzd%llr(ilr)%d%n3, &
!!               lzd%llr(ilr)%d%nfl1, lzd%llr(ilr)%d%nfu1, &
!!               lzd%llr(ilr)%d%nfl2, lzd%llr(ilr)%d%nfu2, &
!!               lzd%llr(ilr)%d%nfl3, lzd%llr(ilr)%d%nfu3, & 
!!               hx, lzd%llr(ilr)%ns1, lzd%llr(ilr)%ns2, lzd%llr(ilr)%ns3, &
!!               lzd%llr(ilr)%bounds%kb%ibyz_c, lzd%llr(ilr)%bounds%kb%ibxz_c, lzd%llr(ilr)%bounds%kb%ibxy_c, &
!!               lzd%llr(ilr)%bounds%kb%ibyz_f, lzd%llr(ilr)%bounds%kb%ibxz_f, lzd%llr(ilr)%bounds%kb%ibxy_f, &
!!               rxyz(1,ilr), confdatarr(iorb)%prefac, .false., 0.d0, &
!!               work_conv%xx_c, work_conv%xx_f1, work_conv%xx_f, &
!!               work_conv%xy_c, work_conv%xy_f2, work_conv%xy_f, &
!!               work_conv%xz_c, work_conv%xz_f4, work_conv%xz_f, &
!!               work_conv%y_c, work_conv%y_f)
!!      else
!!          stop 'wronf conf pot'
!!      end if
!!
!!      call compress_forstandard(lzd%llr(ilr)%d%n1, lzd%llr(ilr)%d%n2, lzd%llr(ilr)%d%n3, &
!!           lzd%llr(ilr)%d%nfl1, lzd%llr(ilr)%d%nfu1, &
!!           lzd%llr(ilr)%d%nfl2, lzd%llr(ilr)%d%nfu2, &
!!           lzd%llr(Ilr)%d%nfl3, lzd%llr(ilr)%d%nfu3, &
!!           lzd%llr(ilr)%wfd%nseg_c, lzd%llr(ilr)%wfd%nvctr_c, &
!!           lzd%llr(ilr)%wfd%keygloc, lzd%llr(ilr)%wfd%keyvloc,  & 
!!           lzd%llr(ilr)%wfd%nseg_f, lzd%llr(ilr)%wfd%nvctr_f, &
!!           lzd%llr(ilr)%wfd%keygloc(1,lzd%llr(ilr)%wfd%nseg_c+min(1,lzd%llr(ilr)%wfd%nseg_f)), &
!!           lzd%llr(ilr)%wfd%keyvloc(lzd%llr(ilr)%wfd%nseg_c+min(1,lzd%llr(ilr)%wfd%nseg_f)),  & 
!!           scal, work_conv%y_c, work_conv%y_f, vpsi(ist_c), vpsi(ist_f))
!!
!!      call deallocate_workarrays_quartic_convolutions(lzd%llr(ilr), subname, work_conv)
!!
!!      ist_f = ist_f + 7*lzd%llr(ilr)%wfd%nvctr_f
!!      ist_c = ist_c + lzd%llr(ilr)%wfd%nvctr_c + 7*lzd%llr(ilr)%wfd%nvctr_f
!!
!!  end do
!!
!!
!!
!!end subroutine apply_orbitaldependent_potential



!!!plot gradient
!!allocate(phiplot(tmb%lzd%glr%wfd%nvctr_c+7*tmb%lzd%glr%wfd%nvctr_f))
!!ist=1
!!do iorb=1,tmbopt%orbs%norbp
!!    iiorb=tmbopt%orbs%isorb+iorb
!!    ilr=tmbopt%orbs%inwhichlocreg(iiorb)
!!    sdim=tmbopt%lzd%llr(ilr)%wfd%nvctr_c+7*tmbopt%lzd%llr(ilr)%wfd%nvctr_f
!!    ldim=tmbopt%lzd%glr%wfd%nvctr_c+7*tmbopt%lzd%glr%wfd%nvctr_f
!!    call to_zero(tmb%lzd%glr%wfd%nvctr_c+7*tmb%lzd%glr%wfd%nvctr_f, phiplot(1))
!!    call Lpsi_to_global2(iproc, nproc, sdim, ldim, tmbopt%orbs%norb, tmbopt%orbs%nspinor, 1, tmbopt%lzd%glr, &
!!         tmbopt%lzd%llr(ilr), lhphiopt(ist), phiplot(1))
!!    !!do istat=1,sdim
!!    !!    write(300,*) lhphiopt(ist+istat-1)
!!    !!end do
!!    !!do istat=1,ldim
!!    !!    write(400,*) phiplot(istat)
!!    !!end do
!!    !!call small_to_large_locreg(iproc, nproc, tmbopt%lzd, tmblarge2%lzd, tmbopt%orbs, tmblarge2%orbs, &
!!    !!     tmbopt%psi, tmblarge2%psi)
!!    write(num,'(i3.3)') iiorb
!!    call plot_wf('gradient'//num,2,at,1.d0,tmbopt%lzd%glr,tmb%lzd%hgrids(1),tmb%lzd%hgrids(2),tmb%lzd%hgrids(3),rxyz,phiplot)
!!    ncount=tmbopt%lzd%llr(ilr)%wfd%nvctr_c+7*tmbopt%lzd%llr(ilr)%wfd%nvctr_f
!!    ist = ist + ncount
!!end do
!!deallocate(phiplot)


!!!subroutine diagonalizeHamiltonian2(iproc, nproc, orbs, nsubmax, HamSmall, ovrlp, eval)
!!!!
!!!! Purpose:
!!!! ========
!!!!   Diagonalizes the Hamiltonian HamSmall and makes sure that all MPI processes give
!!!!   the same result. This is done by requiring that the first entry of each vector
!!!!   is positive.
!!!!
!!!! Calling arguments:
!!!! ==================
!!!!   Input arguments:
!!!!   ----------------
!!!!     iproc     process ID
!!!!     nproc     number of MPI processes
!!!!     orbs      type describing the physical orbitals psi
!!!!   Input / Putput arguments
!!!!     HamSmall  on input: the Hamiltonian
!!!!               on exit: the eigenvectors
!!!!   Output arguments
!!!!     eval      the associated eigenvalues 
!!!!
!!!use module_base
!!!use module_types
!!!implicit none
!!!
!!!! Calling arguments
!!!integer :: iproc, nproc, nsubmax
!!!type(orbitals_data), intent(inout) :: orbs
!!!real(kind=8),dimension(orbs%norb, orbs%norb),intent(inout) :: HamSmall
!!!real(kind=8),dimension(orbs%norb, orbs%norb),intent(in) :: ovrlp
!!!real(kind=8),dimension(orbs%norb),intent(out) :: eval
!!!
!!!! Local variables
!!!integer :: lwork, info, istat, iall, iorb, jorb
!!!real(kind=8),dimension(:),allocatable :: work
!!!character(len=*),parameter :: subname='diagonalizeHamiltonian'
!!!
!!!  ! temp change
!!!  real(8),dimension(:),allocatable:: eval1,beta
!!!  real(8),dimension(:,:), allocatable :: vr,vl,ovrlp_copy
!!!  !real(8),dimension(:,:), allocatable :: inv_ovrlp,ks
!!!  integer :: ierr
!!!  real(8) :: temp, tt, ddot
!!!
!!!  call timing(iproc,'diagonal_seq  ','ON')
!!!
!!!  !! OLD VERSION #####################################################################################################
!!!  ! Get the optimal work array size
!!!  lwork=-1 
!!!  allocate(work(1), stat=istat)
!!!  call memocc(istat, work, 'work', subname)
!!!  call dsygv(1, 'v', 'l', orbs%norb, HamSmall(1,1), orbs%norb, ovrlp(1,1), orbs%norb, eval(1), work(1), lwork, info) 
!!!  lwork=int(work(1))
!!!
!!!  !!! find inverse overlap and premultiply Hamiltonian
!!!!!  allocate(inv_ovrlp(1:orbs%norb,1:orbs%norb))
!!!  !!call dcopy(orbs%norb**2,ovrlp(1,1),1,inv_ovrlp(1,1),1)
!!!  !!! Exact inversion
!!!  !!call dpotrf('l', orbs%norb, inv_ovrlp(1,1), orbs%norb, info)
!!!  !!if(info/=0) then
!!!  !!   write(*,'(1x,a,i0)') 'ERROR in dpotrf, info=',info
!!!  !!   stop
!!!  !!end if
!!!  !!call dpotri('l', orbs%norb, inv_ovrlp(1,1), orbs%norb, info)
!!!  !!if(info/=0) then
!!!  !!   write(*,'(1x,a,i0)') 'ERROR in dpotri, info=',info
!!!  !!   stop
!!!  !!end if
!!!
!!!  !!! fill the upper triangle
!!!  !!do iorb=1,orbs%norb
!!!  !!   do jorb=1,iorb-1
!!!  !!      inv_ovrlp(jorb,iorb)=inv_ovrlp(iorb,jorb)
!!!  !!   end do
!!!  !!end do
!!!
!!!  !!allocate(ks(1:orbs%norb,1:orbs%norb))
!!!  !!call dgemm('n','n', orbs%norb,orbs%norb,orbs%norb,1.d0,inv_ovrlp(1,1),orbs%norb,&
!!!  !!     HamSmall(1,1),orbs%norb,0.d0,ks(1,1),orbs%norb)
!!!  !!call dcopy(orbs%norb**2,ks(1,1),1,HamSmall(1,1),1)
!!!  !!deallocate(ks)
!!!  !!deallocate(inv_ovrlp)
!!!  !!!!!!!!!!!
!!!  !!allocate(ham_copy(1:orbs%norb,1:orbs%norb))
!!!
!!!
!!!  !allocate(ovrlp_copy(1:orbs%norb,1:orbs%norb), stat=istat)
!!!  !call memocc(istat, ovrlp_copy, 'ovrlp_copy', subname)
!!!  !allocate(vl(1:orbs%norb,1:orbs%norb), stat=istat)
!!!  !call memocc(istat, vl, 'vl', subname)
!!!  !allocate(vr(1:orbs%norb,1:orbs%norb), stat=istat)
!!!  !call memocc(istat, vr, 'vr', subname)
!!!  !allocate(eval1(1:orbs%norb), stat=istat)
!!!  !call memocc(istat, eval1, 'eval1', subname)
!!!  !allocate(beta(1:orbs%norb), stat=istat)
!!!  !call memocc(istat, beta, 'beta', subname)
!!!
!!!  !call dcopy(orbs%norb**2, ovrlp(1,1), 1, ovrlp_copy(1,1), 1)
!!!
!!!  !!$call dggev('v', 'v',orbs%norb,&
!!!  !!$      HamSmall(1,1), orbs%norb, ovrlp(1,1), orbs%norb, eval, eval1, beta, &
!!!  !!$      vl,orbs%norb,vr,orbs%norb,work, lwork, ierr)
!!!  !!call DGEEV( 'v','v', orbs%norb, HamSmall(1,1), orbs%norb, eval, eval1, VL, orbs%norb, VR,&
!!!  !!     orbs%norb, WORK, LWORK, ierr )
!!!  !!$lwork=work(1) 
!!!
!!!  ! Deallocate the work array and reallocate it with the optimal size
!!!  iall=-product(shape(work))*kind(work)
!!!  deallocate(work, stat=istat) ; if(istat/=0) stop 'ERROR in deallocating work' 
!!!  call memocc(istat, iall, 'work', subname)
!!!  allocate(work(lwork), stat=istat) ; if(istat/=0) stop 'ERROR in allocating work' 
!!!  call memocc(istat, work, 'work', subname)
!!!
!!!  ! Diagonalize the Hamiltonian
!!!!!  call dcopy(orbs%norb**2, HamSmall(1,1), 1, vl(1,1), 1)
!!!!!  call dcopy(orbs%norb**2, ovrlp(1,1), 1, vr(1,1), 1)
!!!!!  call dcopy(orbs%norb**2, HamSmall(1,1), 1, inv_ovrlp(1,1), 1)
!!!!!$  call dcopy(orbs%norb**2, ovrlp(1,1), 1, ks(1,1), 1)
!!!  call dsygv(1, 'v', 'l', orbs%norb, HamSmall(1,1), orbs%norb, ovrlp(1,1), orbs%norb, eval(1), work(1), lwork, info) 
!!!  iall=-product(shape(work))*kind(work)
!!!  deallocate(work, stat=istat) ; if(istat/=0) stop 'ERROR in deallocating work' 
!!!  call memocc(istat, iall, 'work', subname)
!!!!!  do iorb=1,orbs%norb
!!!!!    do jorb=1,orbs%norb
!!!!!      write(200+iproc,*) iorb,jorb,hamsmall(jorb,iorb)
!!!!!    end do
!!!!!    write(250+iproc,*) iorb,eval(iorb)
!!!!!  end do
!!!!!  call dcopy(orbs%norb**2, vl(1,1), 1, HamSmall(1,1), 1)
!!!!!  call dcopy(orbs%norb**2, vr(1,1), 1, ovrlp(1,1), 1)
!!!
!!!
!!!  !do iorb=1,orbs%norb
!!!  !  eval(iorb) = eval(iorb) / beta(iorb)
!!!  !end do
!!!
!!!!!$$$  lwork=-1
!!!!!$$$  call dggev('v', 'v',orbs%norb,&
!!!!!$$$        HamSmall(1,1), orbs%norb, ovrlp(1,1), orbs%norb, eval, eval1, beta, &
!!!!!$$$        vl,orbs%norb,vr,orbs%norb,work, lwork, ierr)
!!!!!$$$  lwork=work(1) 
!!!!!$$$  iall=-product(shape(work))*kind(work)
!!!!!$$$  deallocate(work, stat=istat) ; if(istat/=0) stop 'ERROR in deallocating work' 
!!!!!$$$  call memocc(istat, iall, 'work', subname)
!!!!!$$$  allocate(work(lwork), stat=istat) ; if(istat/=0) stop 'ERROR in allocating work' 
!!!!!$$$  call memocc(istat, work, 'work', subname)
!!!!!$$$  call dggev('v', 'v',orbs%norb,&
!!!!!$$$        HamSmall(1,1), orbs%norb, ovrlp(1,1), orbs%norb, eval, eval1, beta, &
!!!!!$$$        vl,orbs%norb,vr,orbs%norb,work, lwork, ierr)
!!!!!$$$
!!!!!$$$        hamsmall=vl
!!!!!$$$  do iorb=1,orbs%norb
!!!!!$$$     do jorb=iorb,orbs%norb
!!!!!$$$        if (eval(jorb)/beta(jorb) < eval(iorb)/beta(iorb)) then
!!!!!$$$           temp = eval(iorb)
!!!!!$$$           temp_vec = HamSmall(:,iorb)
!!!!!$$$           eval(iorb) = eval(jorb)
!!!!!$$$           eval(jorb) = temp
!!!!!$$$           HamSmall(:,iorb) = HamSmall(:,jorb)
!!!!!$$$           HamSmall(:,jorb) = temp_vec
!!!!!$$$           temp=beta(iorb)
!!!!!$$$           beta(iorb)=beta(jorb)
!!!!!$$$           beta(jorb)=temp
!!!!!$$$        end if
!!!!!$$$     end do
!!!!!$$$  end do
!!!!!$$$
!!!!!$$$
!!!!!$$$
!!!!!$$$  call dcopy(orbs%norb**2, ks(1,1), 1, ovrlp(1,1), 1)
!!!!!$$$  do iorb=1,orbs%norb
!!!!!$$$      call dgemv('n', orbs%norb, orbs%norb, 1.d0, ovrlp(1,1), &
!!!!!$$$           orbs%norb, hamsmall(1,iorb), 1, 0.d0, vl(1,iorb), 1)
!!!!!$$$      tt=ddot(orbs%norb, hamsmall(1,iorb),  1, vl(1,iorb), 1)
!!!!!$$$      call dscal(orbs%norb, 1/sqrt(tt), hamsmall(1,iorb), 1)
!!!!!$$$  end do
!!!!!$$$
!!!!!$$$
!!!!!$$$
!!!!!$$$
!!!!!$$$!!  do iorb=1,orbs%norb
!!!!!$$$!!    do jorb=1,orbs%norb
!!!!!$$$!!      write(300+iproc,*) iorb,jorb,hamsmall(jorb,iorb)
!!!!!$$$!!    end do
!!!!!$$$!!    write(350+iproc,*) iorb,eval(iorb)/beta(iorb)
!!!!!$$$!!  end do
!!!!!$$$!!
!!!!!$$$!!  lwork=-1
!!!!!$$$!!  call dcopy(orbs%norb**2, inv_ovrlp(1,1), 1, HamSmall(1,1), 1)
!!!!!$$$!!  call dcopy(orbs%norb**2, ks(1,1), 1, ovrlp(1,1), 1)
!!!!!$$$!!  call DGEEV( 'v','v', orbs%norb, HamSmall(1,1), orbs%norb, eval, eval1, VL, orbs%norb, VR,&
!!!!!$$$!!       orbs%norb, WORK, LWORK, ierr )
!!!!!$$$!!  ! Deallocate the work array and reallocate it with the optimal size
!!!!!$$$!!  lwork=work(1) 
!!!!!$$$!!  iall=-product(shape(work))*kind(work)
!!!!!$$$!!  deallocate(work, stat=istat) ; if(istat/=0) stop 'ERROR in deallocating work' 
!!!!!$$$!!  call memocc(istat, iall, 'work', subname)
!!!!!$$$!!  allocate(work(lwork), stat=istat) ; if(istat/=0) stop 'ERROR in allocating work' 
!!!!!$$$!!  call memocc(istat, work, 'work', subname)
!!!!!$$$!!
!!!!!$$$!!  call DGEEV( 'v','v', orbs%norb, HamSmall(1,1), orbs%norb, eval, eval1, VL, orbs%norb, VR,&
!!!!!$$$!!       orbs%norb, WORK, LWORK, ierr )
!!!!!$$$!!
!!!!!$$$!!  HamSmall=vl
!!!!!$$$!!  do iorb=1,orbs%norb
!!!!!$$$!!     do jorb=iorb,orbs%norb
!!!!!$$$!!        if (eval(jorb) < eval(iorb)) then
!!!!!$$$!!           temp = eval(iorb)
!!!!!$$$!!           temp_vec = HamSmall(:,iorb)
!!!!!$$$!!           eval(iorb) = eval(jorb)
!!!!!$$$!!           eval(jorb) = temp
!!!!!$$$!!           HamSmall(:,iorb) = HamSmall(:,jorb)
!!!!!$$$!!           HamSmall(:,jorb) = temp_vec
!!!!!$$$!!        end if
!!!!!$$$!!     end do
!!!!!$$$!!  end do
!!!!!$$$!!
!!!!!$$$!!  do iorb=1,orbs%norb
!!!!!$$$!!    do jorb=1,orbs%norb
!!!!!$$$!!      write(400+iproc,*) iorb,jorb,hamsmall(jorb,iorb)
!!!!!$$$!!    end do
!!!!!$$$!!    write(450+iproc,*) iorb,eval(iorb)
!!!!!$$$!!  end do
!!!!!$$$!!
!!!!!$$$!!!  do iorb=1,orbs%norb
!!!!!$$$!!!    write(36,*) vl(:,iorb)
!!!!!$$$!!!    write(37,*) vr(:,iorb)
!!!!!$$$!!!  end do
!!!!!$$$!!!  write(36,*) ''
!!!!!$$$!!!  write(37,*) ''
!!!!!$$$!!!  write(38,*) 'eval',eval
!!!!!$$$!!!  write(38,*) 'eval1',eval1
!!!!!$$$!!!  !write(38,*) 'beta',beta
!!!!!$$$
!!!!!$$$  ! Deallocate the work array.
!!!!!$$$  iall=-product(shape(work))*kind(work)
!!!!!$$$  deallocate(work, stat=istat) ; if(istat/=0) stop 'ERROR in deallocating work' 
!!!!!$$$  call memocc(istat, iall, 'work', subname)
!!!!!$$$  
!!!!!$$$  ! Make sure that the eigenvectors are the same for all MPI processes. To do so, require that 
!!!!!$$$  ! the first entry of each vector is positive.
!!!!!$$$  do iorb=1,orbs%norb
!!!!!$$$      if(HamSmall(1,iorb)<0.d0) then
!!!!!$$$          do jorb=1,orbs%norb
!!!!!$$$              HamSmall(jorb,iorb)=-HamSmall(jorb,iorb)
!!!!!$$$          end do
!!!!!$$$      end if
!!!!!$$$  end do
!!!!!$$$  !! #################################################################################################################
!!!!!$$$
!!!!!$$$  deallocate(vl)
!!!!!$$$  deallocate(vr)
!!!!!$$$  deallocate(eval1)
!!!!!$$$  deallocate(beta)
!!!!!$$$
!!!!!$$$if (.false.) then
!!!!!$$$  allocate(vl(1:orbs%norb,1:orbs%norb))
!!!!!$$$  allocate(vr(1:orbs%norb,1:orbs%norb))
!!!!!$$$  allocate(eval1(1:orbs%norb))
!!!!!$$$  allocate(eval2(1:orbs%norb))
!!!!!$$$
!!!!!$$$
!!!!!$$$  allocate(work(lwork), stat=istat) ; if(istat/=0) stop 'ERROR in allocating work' 
!!!!!$$$  ! check eigenvalues of overlap matrix
!!!!!$$$      call DGEEV( 'v','v', orbs%norb, ovrlp(1,1), orbs%norb, eval2, eval1, VL, orbs%norb, VR,&
!!!!!$$$                  orbs%norb, WORK, LWORK, ierr )
!!!!!$$$!  write(40,*) 'eval',eval2
!!!!!$$$!  write(40,*) 'eval1',eval1
!!!!!$$$!  write(40,*) 'sum',sum(eval2)
!!!!!$$$
!!!!!$$$!  do iorb=1,orbs%norb
!!!!!$$$!    write(44,*) vl(:,iorb)
!!!!!$$$!    write(45,*) vr(:,iorb)
!!!!!$$$!  end do
!!!!!$$$!  write(44,*) ''
!!!!!$$$!  write(45,*) ''
!!!!!$$$
!!!!!$$$  write(41,*) 'sum olap eigs',sum(eval2)
!!!!!$$$
!!!!!$$$  deallocate(work)
!!!!!$$$  deallocate(vl)
!!!!!$$$  deallocate(vr)
!!!!!$$$  deallocate(eval1)
!!!!!$$$  deallocate(eval2)
!!!!!$$$end if
!!!
!!!  !!!! NEW VERSION #####################################################################################################
!!!  !!! Determine the maximal number of non-zero subdiagonals
!!!  !!!!nsubmax=0
!!!  !!!!do iorb=1,orbs%norb
!!!  !!!!    nsub=0
!!!  !!!!    do jorb=orbs%norb,iorb+1,-1
!!!  !!!!        if(Hamsmall(jorb,iorb)/=0.d0) then
!!!  !!!!            nsub=jorb-iorb
!!!  !!!!            exit
!!!  !!!!        end if
!!!  !!!!    end do
!!!  !!!!    if(iproc==0) write(*,*) 'iorb,nsub',iorb,nsub
!!!  !!!!    nsubmax=max(nsub,nsubmax)
!!!  !!!!end do
!!!  !!!!if(iproc==0) write(*,*) 'nsubmax',nsubmax
!!!  !!!!if(iproc==0) then
!!!  !!!!      do iorb=1,orbs%norb
!!!  !!!!           write(*,'(14es10.3)') (hamsmall(iorb,jorb), jorb=1,orbs%norb)
!!!  !!!!      end do
!!!  !!!!end if
!!!
!!!  !!! Copy to banded format
!!!  !!allocate(ham_band(nsubmax+1,orbs%norb), stat=istat)
!!!  !!call memocc(istat, ham_band, 'ham_band', subname)
!!!  !!allocate(ovrlp_band(nsubmax+1,orbs%norb), stat=istat)
!!!  !!call memocc(istat, ovrlp_band, 'ovrlp_band', subname)
!!!  !!do iorb=1,orbs%norb
!!!  !!    do jorb=iorb,min(iorb+nsubmax,orbs%norb)
!!!  !!        ham_band(1+jorb-iorb,iorb)=HamSmall(jorb,iorb)
!!!  !!        ovrlp_band(1+jorb-iorb,iorb)=ovrlp(jorb,iorb)
!!!  !!    end do
!!!  !!end do
!!!  !!!!if(iproc==0) then
!!!  !!!!      write(*,*) '+++++++++++++++++++++++++++++'
!!!  !!!!      do iorb=1,nsubmax+1
!!!  !!!!           write(*,'(14es10.3)') (ham_band(iorb,jorb), jorb=1,orbs%norb)
!!!  !!!!      end do
!!!  !!!!end if
!!!
!!!
!!!  !!!!! Get the optimal work array size
!!!  !!!!lwork=-1 
!!!  !!!!allocate(work(1), stat=istat)
!!!  !!!!call memocc(istat, work, 'work', subname)
!!!  !!!!call dsygv(1, 'v', 'l', orbs%norb, HamSmall(1,1), orbs%norb, ovrlp(1,1), orbs%norb, eval(1), work(1), lwork, info) 
!!!  !!!!lwork=work(1) 
!!!
!!!  !!!!! Deallocate the work array ane reallocate it with the optimal size
!!!  !!!!iall=-product(shape(work))*kind(work)
!!!  !!!!deallocate(work, stat=istat) ; if(istat/=0) stop 'ERROR in deallocating work' 
!!!  !!!!call memocc(istat, iall, 'work', subname)
!!!  !!allocate(work(3*orbs%norb), stat=istat) ; if(istat/=0) stop 'ERROR in allocating work' 
!!!  !!call memocc(istat, work, 'work', subname)
!!!
!!!  !!! Diagonalize the Hamiltonian
!!!  !!!call dsygv(1, 'v', 'l', orbs%norb, HamSmall(1,1), orbs%norb, ovrlp(1,1), orbs%norb, eval(1), work(1), lwork, info) 
!!!  !!call dsbgv('v', 'l', orbs%norb, nsubmax, nsubmax, ham_band(1,1), nsubmax+1, ovrlp_band(1,1), nsubmax+1, &
!!!  !!     eval(1), HamSmall(1,1), orbs%norb, work, info)
!!!
!!!  !!! Deallocate the work array.
!!!  !!iall=-product(shape(work))*kind(work)
!!!  !!deallocate(work, stat=istat) ; if(istat/=0) stop 'ERROR in deallocating work' 
!!!  !!call memocc(istat, iall, 'work', subname)
!!!  !!
!!!  !!! Make sure that the eigenvectors are the same for all MPI processes. To do so, require that 
!!!  !!! the first entry of each vector is positive.
!!!  !!do iorb=1,orbs%norb
!!!  !!    if(HamSmall(1,iorb)<0.d0) then
!!!  !!        do jorb=1,orbs%norb
!!!  !!            HamSmall(jorb,iorb)=-HamSmall(jorb,iorb)
!!!  !!        end do
!!!  !!    end if
!!!  !!end do
!!!
!!!
!!!  !!iall=-product(shape(ham_band))*kind(ham_band)
!!!  !!deallocate(ham_band, stat=istat) ; if(istat/=0) stop 'ERROR in deallocating ham_band' 
!!!  !!call memocc(istat, iall, 'ham_band', subname)
!!!
!!!  !!iall=-product(shape(ovrlp_band))*kind(ovrlp_band)
!!!  !!deallocate(ovrlp_band, stat=istat) ; if(istat/=0) stop 'ERROR in deallocating ovrlp_band' 
!!!  !!call memocc(istat, iall, 'ovrlp_band', subname)
!!!
!!!  call timing(iproc,'diagonal_seq  ','OF')
!!!
!!!end subroutine diagonalizeHamiltonian2



!!subroutine plot_gradient(iproc, nproc, num, lzd, orbs, psi)
!!  use module_base
!!  use module_types
!!  use module_interfaces
!!  implicit none
!!
!!  ! Calling arguments
!!  integer,intent(in):: iproc, nproc, num
!!  type(local_zone_descriptors),intent(in):: lzd
!!  type(orbitals_data),intent(in):: orbs
!!  real(8),dimension(orbs%npsidim_orbs),intent(in):: psi
!!
!!  ! Local variables
!!  integer:: istc, istf, iorb, iiorb, ilr, i1, i2, i3, istat, iall, ierr, ii2, ii3
!!  real(8):: r0, r1, r2, r3, rr, tt, gnrm
!!  real(8),dimension(:,:,:,:,:,:),allocatable:: psig
!!  character(len=*),parameter:: subname='flatten_at_boundaries'
!!  
!!
!!  istc=1
!!  istf=1
!!  do iorb=1,orbs%norbp
!!      iiorb=orbs%isorb+iorb
!!      ilr=orbs%inwhichlocreg(iiorb)
!!
!!      allocate(psig(0:lzd%llr(ilr)%d%n1,2,0:lzd%llr(ilr)%d%n2,2,0:lzd%llr(ilr)%d%n3,2), stat=istat)
!!      call memocc(istat, psig, 'psig', subname)
!!      call to_zero(8*(lzd%llr(ilr)%d%n1+1)*(lzd%llr(ilr)%d%n2+1)*(lzd%llr(ilr)%d%n3+1), psig(0,1,0,1,0,1))
!!
!!      istf = istf + lzd%llr(ilr)%wfd%nvctr_c
!!      call uncompress(lzd%llr(ilr)%d%n1, lzd%llr(ilr)%d%n2, lzd%llr(ilr)%d%n3, &
!!           lzd%llr(ilr)%wfd%nseg_c, lzd%llr(ilr)%wfd%nvctr_c, lzd%llr(ilr)%wfd%keygloc, lzd%llr(ilr)%wfd%keyvloc,  &
!!           lzd%llr(ilr)%wfd%nseg_f, lzd%llr(ilr)%wfd%nvctr_f, &
!!           lzd%llr(ilr)%wfd%keygloc(1,lzd%llr(ilr)%wfd%nseg_c+min(1,lzd%llr(ilr)%wfd%nseg_f)), &
!!           lzd%llr(ilr)%wfd%keyvloc(lzd%llr(ilr)%wfd%nseg_c+min(1,lzd%llr(ilr)%wfd%nseg_f)), &
!!           psi(istc), psi(istf), psig)
!!
!!      ii2=nint(lzd%llr(ilr)%d%n2/2.d0)
!!      ii3=nint(lzd%llr(ilr)%d%n3/2.d0)
!!      do i1=0,lzd%llr(ilr)%d%n1
!!          write(num+iiorb,*) i1, psig(i1,1,ii2,1,ii3,1)
!!      end do
!!
!!      !!call compress(lzd%llr(ilr)%d%n1, lzd%llr(ilr)%d%n2, &
!!      !!     0, lzd%llr(ilr)%d%n1, 0, lzd%llr(ilr)%d%n2, 0, lzd%llr(ilr)%d%n3, &
!!      !!     lzd%llr(ilr)%wfd%nseg_c, lzd%llr(ilr)%wfd%nvctr_c, lzd%llr(ilr)%wfd%keygloc, lzd%llr(ilr)%wfd%keyvloc,  &
!!      !!     lzd%llr(ilr)%wfd%nseg_f, lzd%llr(ilr)%wfd%nvctr_f, &
!!      !!     lzd%llr(ilr)%wfd%keygloc(1,lzd%llr(ilr)%wfd%nseg_c+min(1,lzd%llr(ilr)%wfd%nseg_f)), &
!!      !!     lzd%llr(ilr)%wfd%keyvloc(lzd%llr(ilr)%wfd%nseg_c+min(1,lzd%llr(ilr)%wfd%nseg_f)),  &
!!      !!     psig, psi(istc), psi(istf))
!!
!!      istf = istf + 7*lzd%llr(ilr)%wfd%nvctr_f
!!      istc = istc + lzd%llr(ilr)%wfd%nvctr_c + 7*lzd%llr(ilr)%wfd%nvctr_f
!!
!!      iall=-product(shape(psig))*kind(psig)
!!      deallocate(psig,stat=istat)
!!      call memocc(istat,iall,'psig',subname)
!!
!!
!!  end do
!!
!!end subroutine plot_gradient

               
              
!!subroutine project_gradient(iproc, nproc, tmb, lphi, lhphi)
!!  use module_base
!!  use module_types
!!  implicit none
!!
!!  integer,intent(in):: iproc, nproc
!!  ! Calling arguments
!!  type(DFT_wavefunction),intent(in):: tmb
!!  real(8),dimension(tmb%orbs%npsidim_orbs),intent(in):: lphi, lhphi
!!
!!  ! Local variables
!!  real(8),dimension(:),allocatable:: psit_c, psit_f, hpsit_c, hpsit_f, hpsittmp_c, hpsittmp_f, lhphitmp
!!  real(8),dimension(:,:),allocatable:: lagmat
!!  real(8):: fnrm, ddot
!!  integer:: istat, iall, istart, iorb, iiorb, ilr, ncount, ierr, jorb
!!  character(len=*),parameter:: subname='project_gradient'
!!
!!  ! Calculate matrix <gradient|TMBs>
!!  allocate(psit_c(sum(tmb%collcom%nrecvcounts_c)), stat=istat)
!!  call memocc(istat, psit_c, 'psit_c', subname)
!!  allocate(psit_f(7*sum(tmb%collcom%nrecvcounts_f)), stat=istat)
!!  call memocc(istat, psit_f, 'psit_f', subname)
!!  call transpose_localized(iproc, nproc, tmb%orbs, tmb%collcom, lphi, psit_c, psit_f, tmb%lzd)
!!  allocate(hpsit_c(sum(tmb%collcom%nrecvcounts_c)), stat=istat)
!!  call memocc(istat, hpsit_c, 'hpsit_c', subname)
!!  allocate(hpsit_f(7*sum(tmb%collcom%nrecvcounts_f)), stat=istat)
!!  call memocc(istat, hpsit_f, 'hpsit_f', subname)
!!  call transpose_localized(iproc, nproc, tmb%orbs, tmb%collcom, lhphi, hpsit_c, hpsit_f, tmb%lzd)
!!
!!  allocate(lagmat(tmb%orbs%norb,tmb%orbs%norb), stat=istat)
!!  call memocc(istat, lagmat, 'lagmat', subname)
!!  call calculate_overlap_transposed(iproc, nproc, tmb%orbs, tmb%mad, tmb%collcom, psit_c, hpsit_c, psit_f, hpsit_f, lagmat)
!!  if(iproc==0) then
!!      do iorb=1,tmb%orbs%norb
!!          do jorb=1,tmb%orbs%norb
!!              write(300,*) iorb, jorb, lagmat(jorb,iorb)
!!          end do
!!      end do
!!  end if
!!
!!  allocate(hpsittmp_c(sum(tmb%collcom%nrecvcounts_c)), stat=istat)
!!  call memocc(istat, hpsittmp_c, 'hpsittmp_c', subname)
!!  allocate(hpsittmp_f(7*sum(tmb%collcom%nrecvcounts_f)), stat=istat)
!!  call memocc(istat, hpsittmp_f, 'hpsittmp_f', subname)
!!  call build_linear_combination_transposed(tmb%orbs%norb, lagmat, tmb%collcom, &
!!       hpsit_c, hpsit_f, .true., hpsittmp_c, hpsittmp_f, iproc)
!!
!!  allocate(lhphitmp(tmb%orbs%npsidim_orbs), stat=istat)
!!  call memocc(istat, lhphitmp, 'lhphitmp', subname)
!!  call untranspose_localized(iproc, nproc, tmb%orbs, tmb%collcom, hpsittmp_c, hpsittmp_f, lhphitmp, tmb%lzd)
!!
!!
!!  fnrm=0.d0
!!  istart=1
!!  do iorb=1,tmb%orbs%norbp
!!      iiorb=tmb%orbs%isorb+iorb
!!      ilr=tmb%orbs%inwhichlocreg(iiorb)
!!      ncount=tmb%lzd%llr(ilr)%wfd%nvctr_c+7*tmb%lzd%llr(ilr)%wfd%nvctr_f
!!      fnrm=fnrm+ddot(ncount, lhphitmp(istart), 1, lhphitmp(istart), 1)
!!      istart=istart+ncount
!!  end do
!!  call mpiallred(fnrm, 1, mpi_sum, mpi_comm_world, ierr)
!!  fnrm=sqrt(fnrm/dble(tmb%orbs%norb))
!!  if(iproc==0) write(*,*) 'projected gradient:',fnrm
!!
!!
!!
!!  iall=-product(shape(psit_c))*kind(psit_c)
!!  deallocate(psit_c,stat=istat)
!!  call memocc(istat,iall,'psit_c',subname)
!!
!!  iall=-product(shape(psit_f))*kind(psit_f)
!!  deallocate(psit_f,stat=istat)
!!  call memocc(istat,iall,'psit_f',subname)
!!
!!  iall=-product(shape(hpsit_c))*kind(hpsit_c)
!!  deallocate(hpsit_c,stat=istat)
!!  call memocc(istat,iall,'hpsit_c',subname)
!!
!!  iall=-product(shape(hpsit_f))*kind(hpsit_f)
!!  deallocate(hpsit_f,stat=istat)
!!  call memocc(istat,iall,'hpsit_f',subname)
!!
!!  iall=-product(shape(hpsittmp_c))*kind(hpsittmp_c)
!!  deallocate(hpsittmp_c,stat=istat)
!!  call memocc(istat,iall,'hpsittmp_c',subname)
!!
!!  iall=-product(shape(hpsittmp_f))*kind(hpsittmp_f)
!!  deallocate(hpsittmp_f,stat=istat)
!!  call memocc(istat,iall,'hpsittmp_f',subname)
!!
!!  iall=-product(shape(lhphitmp))*kind(lhphitmp)
!!  deallocate(lhphitmp,stat=istat)
!!  call memocc(istat,iall,'lhphitmp',subname)
!!
!!  iall=-product(shape(lagmat))*kind(lagmat)
!!  deallocate(lagmat,stat=istat)
!!  call memocc(istat,iall,'lagmat',subname)
!!
!!endsubroutine project_gradient


subroutine build_new_linear_combinations(iproc, lzd, orbs, op, nrecvbuf, recvbuf, omat, reset, lphi)
use module_base
use module_types
implicit none

!Calling arguments
integer,intent(in) :: iproc
type(local_zone_descriptors),intent(in) :: lzd
type(orbitals_data),intent(in) :: orbs
type(overlapParameters),intent(in) :: op
integer,intent(in) :: nrecvbuf
real(kind=8),dimension(nrecvbuf),intent(in) :: recvbuf
real(kind=8),dimension(orbs%norb,orbs%norb),intent(in) :: omat
logical,intent(in) :: reset
real(kind=8),dimension(max(orbs%npsidim_orbs,orbs%npsidim_comp)),intent(out) :: lphi

! Local variables
integer :: jst, ilrold, iorb, iiorb, ilr, ncount, jorb, jjorb, ldim, indout, gdim, iorbref
integer :: istart, iend, iseg, start, kseg, kold, kstart, kend 
real(kind=8) :: tt

   call timing(iproc,'build_lincomb ','ON')

      ! Build new lphi
      if(reset) then
          !!lphi=0.d0
          call to_zero(max(orbs%npsidim_orbs,orbs%npsidim_comp), lphi(1))
      end if

      indout=1
      ilrold=-1
      do iorb=1,orbs%norbp
          iiorb=orbs%isorb+iorb
          ilr=orbs%inwhichlocreg(iiorb)
          !if(ilr>ilrold) then
          if(ilr/=ilrold) then
              iorbref=iorb
          end if
          gdim=lzd%llr(ilr)%wfd%nvctr_c+7*lzd%llr(ilr)%wfd%nvctr_f
          do jorb=1,op%noverlaps(iiorb)
              jjorb=op%overlaps(jorb,iiorb)
              !jjorb=op%overlaps(jorb,ilr)
              jst=op%indexInRecvBuf(iorbref,jjorb)
              ldim=op%wfd_overlap(jorb,iorbref)%nvctr_c+7*op%wfd_overlap(jorb,iorbref)%nvctr_f
              tt=omat(jjorb,iiorb)
              !!tt=tt*lzd%cutoffweight(jjorb,iiorb)
              !Coarse part
              kold=1
              do iseg=1,op%wfd_overlap(jorb,iorbref)%nseg_c
                  istart=op%wfd_overlap(jorb,iorbref)%keyglob(1,iseg)
                  iend=op%wfd_overlap(jorb,iorbref)%keyglob(2,iseg)
                  ncount=iend-istart+1
                  inner_loop: do kseg=kold,lzd%llr(ilr)%wfd%nseg_c
                     kstart = lzd%llr(ilr)%wfd%keyglob(1,kseg)
                     kend   = lzd%llr(ilr)%wfd%keyglob(2,kseg)
                     if(kstart <= iend .and. kend >= istart) then 
                        kold = kseg + 1
                        start = lzd%llr(ilr)%wfd%keyvglob(kseg) + max(0,istart-kstart)
                        call daxpy(ncount, tt, recvBuf(jst), 1, lphi(indout+start-1), 1)
                        jst=jst+ncount
                        exit inner_loop
                     end if
                  end do inner_loop
              end do
              ! Fine part
              kold = 1
              jst=op%indexInRecvBuf(iorbref,jjorb)              
              do iseg=1,op%wfd_overlap(jorb,iorbref)%nseg_f
                 istart=op%wfd_overlap(jorb,iorbref)%keyglob(1,iseg+op%wfd_overlap(jorb,iorbref)%nseg_c)
                 iend=op%wfd_overlap(jorb,iorbref)%keyglob(2,iseg+op%wfd_overlap(jorb,iorbref)%nseg_c)
                 start = op%wfd_overlap(jorb,iorbref)%keyvglob(iseg+op%wfd_overlap(jorb,iorbref)%nseg_c)
                 ncount=7*(iend-istart+1)
                 inner_loop2: do kseg=kold,lzd%llr(ilr)%wfd%nseg_f
                    kstart = lzd%llr(ilr)%wfd%keyglob(1,kseg+lzd%llr(ilr)%wfd%nseg_c)
                    kend   = lzd%llr(ilr)%wfd%keyglob(2,kseg+lzd%llr(ilr)%wfd%nseg_c)
                    if(kstart <= iend .and. kend >= istart) then 
                       kold = kseg + 1
                       start = lzd%llr(ilr)%wfd%nvctr_c+(lzd%llr(ilr)%wfd%keyvglob(kseg+lzd%llr(ilr)%wfd%nseg_c) +&
                              max(0,istart-kstart)-1)*7
                       call daxpy(ncount, tt, recvBuf(jst+op%wfd_overlap(jorb,iorbref)%nvctr_c), 1,&
                               lphi(indout+start), 1)
                       jst=jst+ncount
                       exit inner_loop2
                    end if
                  end do inner_loop2
              end do
          end do
          indout=indout+gdim
          ilrold=ilr

      end do

   call timing(iproc,'build_lincomb ','OF')
          

end subroutine build_new_linear_combinations


!!subroutine communicate_basis_for_density(iproc, nproc, lzd, llborbs, lphi, comsr)
!!  use module_base
!!  use module_types
!!  use module_interfaces, except_this_one => communicate_basis_for_density
!!  implicit none
!!  
!!  ! Calling arguments
!!  integer,intent(in) :: iproc, nproc
!!  type(local_zone_descriptors),intent(in) :: lzd
!!  type(orbitals_data),intent(in) :: llborbs
!!  real(kind=8),dimension(llborbs%npsidim_orbs),intent(in) :: lphi
!!  type(p2pComms),intent(inout) :: comsr
!!  
!!  ! Local variables
!!  integer :: ist, istr, iorb, iiorb, ilr
!!  type(workarr_sumrho) :: w
!!
!!  call timing(iproc,'commbasis4dens','ON') !lr408t
!!
!!  ! Allocate the communication buffers for the calculation of the charge density.
!!  !call allocateCommunicationbufferSumrho(iproc, comsr, subname)
!!  ! Transform all orbitals to real space.
!!  ist=1
!!  istr=1
!!  do iorb=1,llborbs%norbp
!!      iiorb=llborbs%isorb+iorb
!!      ilr=llborbs%inWhichLocreg(iiorb)
!!      call initialize_work_arrays_sumrho(lzd%Llr(ilr), w)
!!      call daub_to_isf(lzd%Llr(ilr), w, lphi(ist), comsr%sendBuf(istr))
!!      call deallocate_work_arrays_sumrho(w)
!!      ist = ist + lzd%Llr(ilr)%wfd%nvctr_c + 7*lzd%Llr(ilr)%wfd%nvctr_f
!!      istr = istr + lzd%Llr(ilr)%d%n1i*lzd%Llr(ilr)%d%n2i*lzd%Llr(ilr)%d%n3i
!!  end do
!!  if(istr/=comsr%nsendBuf+1) then
!!      write(*,'(a,i0,a)') 'ERROR on process ',iproc,' : istr/=comsr%nsendBuf+1'
!!      stop
!!  end if
!!  
!!  ! Post the MPI messages for the communication of sumrho. Since we use non blocking point
!!  ! to point communication, the program will continue immediately. The messages will be gathered
!!  ! in the subroutine sumrhoForLocalizedBasis2.
!!  !!call postCommunicationSumrho2(iproc, nproc, comsr, comsr%sendBuf, comsr%recvBuf)
!!  call post_p2p_communication(iproc, nproc, comsr%nsendbuf, comsr%sendbuf, comsr%nrecvbuf, comsr%recvbuf, comsr)
!!
!!  call timing(iproc,'commbasis4dens','OF') !lr408t
!!
!!end subroutine communicate_basis_for_density




subroutine check_idempotency(iproc, nproc, tmb, diff)
  use module_base
  use module_types
  use sparsematrix_base, only: sparsematrix_malloc_ptr, DENSE_FULL, assignment(=)
  use sparsematrix, only: uncompress_matrix
  implicit none

  ! Calling variables
  integer,intent(in) :: iproc, nproc
  type(DFT_wavefunction),intent(inout) :: tmb
  real(kind=8),intent(out) :: diff

  ! Local variables
  integer :: iorb, iiorb, jsegstart, jsegend, jseg, jorb, jjorb, ierr
  real(kind=8),dimension(:,:),allocatable :: ks, ksk, ksksk


  tmb%linmat%ovrlp%matrix=f_malloc_ptr((/tmb%orbs%norb,tmb%orbs%norb/),id='tmb%linmat%ovrlp%matrix')
  tmb%linmat%kernel_%matrix = sparsematrix_malloc_ptr(tmb%linmat%l, iaction=DENSE_FULL, id='tmb%linmat%kernel_%matrix')
  call uncompress_matrix(iproc,tmb%linmat%ovrlp)
  call uncompress_matrix(iproc,tmb%linmat%denskern_large)


  call dscal(tmb%orbs%norb**2, 0.5d0, tmb%linmat%kernel_%matrix, 1)

  ks=f_malloc((/tmb%orbs%norb,tmb%orbs%norb/))
  ksk=f_malloc((/tmb%orbs%norb,tmb%orbs%norbp/))
  ksksk=f_malloc((/tmb%orbs%norb,tmb%orbs%norb/))

  call to_zero(tmb%orbs%norb**2, ks(1,1))
  if (tmb%orbs%norbp>0) then
      call dgemm('n', 'n', tmb%orbs%norb, tmb%orbs%norbp, tmb%orbs%norb, &
                 1.d0, tmb%linmat%kernel_%matrix(1,1), tmb%orbs%norb, &
                 tmb%linmat%ovrlp%matrix(1,tmb%orbs%isorb+1), tmb%orbs%norb, &
                 0.d0, ks(1,tmb%orbs%isorb+1), tmb%orbs%norb) 
  end if
  call mpiallred(ks(1,1), tmb%orbs%norb**2, mpi_sum, bigdft_mpi%mpi_comm, ierr)
  if (tmb%orbs%norbp>0) then
      call dgemm('n', 'n', tmb%orbs%norb, tmb%orbs%norbp, tmb%orbs%norb, &
                 1.d0, ks(1,1), tmb%orbs%norb, &
                 tmb%linmat%kernel_%matrix(1,tmb%orbs%isorb+1), tmb%orbs%norb, &
                 0.d0, ksk(1,1), tmb%orbs%norb)
  end if
  !!if (tmb%orbs%norbp>0) then
  !!    call dgemm('n', 'n', tmb%orbs%norb, tmb%orbs%norbp, tmb%orbs%norb, 1.d0, ks(1,1), tmb%orbs%norb, &
  !!               ksk(1,1), tmb%orbs%norb, 0.d0, ksksk(1,1), tmb%orbs%norb)
  !!end if


  diff=0.d0
  do iorb=tmb%orbs%isorb+1,tmb%orbs%isorb+tmb%orbs%norbp
      iiorb=iorb-tmb%orbs%isorb
      jsegstart=tmb%linmat%denskern_large%istsegline(iorb)
      if (iorb<tmb%orbs%norb) then
          jsegend=tmb%linmat%denskern_large%istsegline(iorb+1)-1
      else
          jsegend=tmb%linmat%denskern_large%nseg
      end if
      do jseg=jsegstart,jsegend
          do jorb=tmb%linmat%denskern_large%keyg(1,jseg),tmb%linmat%denskern_large%keyg(2,jseg)
              jjorb=jorb-(iorb-1)*tmb%orbs%norb
              diff = diff + (ksk(jjorb,iiorb)-tmb%linmat%kernel_%matrix(jjorb,iorb))**2
          end do
      end do
  end do
  call mpiallred(diff, 1, mpi_sum, bigdft_mpi%mpi_comm, ierr)
  diff=sqrt(diff)


  call f_free(ks)
  call f_free(ksk)
  call f_free(ksksk)
  call f_free_ptr(tmb%linmat%ovrlp%matrix)
  call f_free_ptr(tmb%linmat%kernel_%matrix)

end subroutine check_idempotency



subroutine my_geocode_buffers(geocode,nl1,nl2,nl3)
  implicit none
  integer, intent(out) :: nl1,nl2,nl3
  character(len=1), intent(in) :: geocode !< @copydoc poisson_solver::doc::geocode
  !local variables
  logical :: perx,pery,perz
  integer :: nr1,nr2,nr3

  !conditions for periodicity in the three directions
  perx=(geocode /= 'F')
  pery=(geocode == 'P')
  perz=(geocode /= 'F')

  call ext_buffers(perx,nl1,nr1)
  call ext_buffers(pery,nl2,nr2)
  call ext_buffers(perz,nl3,nr3)

end subroutine my_geocode_buffers

subroutine getLocalizedBasis(iproc,nproc,at,orbs,rxyz,denspot,GPU,trH,trH_old,&
    fnrm,infoBasisFunctions,nlpsp,scf_mode,ldiis,SIC,tmb,energs_base,&
    nit_precond,target_function,&
    correction_orthoconstraint,nit_basis,&
    ratio_deltas,ortho_on,extra_states,itout,conv_crit,experimental_mode,early_stop,&
    gnrm_dynamic, min_gnrm_for_dynamic, can_use_ham, order_taylor, kappa_conv, method_updatekernel,&
    purification_quickreturn, correction_co_contra)
  !
  ! Purpose:
  ! ========
  !   Calculates the localized basis functions phi. These basis functions are obtained by adding a
  !   quartic potential centered on the atoms to the ordinary Hamiltonian. The eigenfunctions are then
  !   determined by minimizing the trace until the gradient norm is below the convergence criterion.
  use module_base
  use module_types
  use yaml_output
  use module_interfaces, except_this_one => getLocalizedBasis, except_this_one_A => writeonewave
  use communications, only: transpose_localized, start_onesided_communication
  use sparsematrix_base, only: assignment(=), sparsematrix_malloc, SPARSE_FULL
  !  use Poisson_Solver
  !use allocModule
  implicit none

  ! Calling arguments
  integer,intent(in) :: iproc, nproc, order_taylor
  integer,intent(out) :: infoBasisFunctions
  type(atoms_data), intent(in) :: at
  type(orbitals_data) :: orbs
  real(kind=8),dimension(3,at%astruct%nat) :: rxyz
  type(DFT_local_fields), intent(inout) :: denspot
  type(GPU_pointers), intent(inout) :: GPU
  real(kind=8),intent(out) :: trH, fnrm
  real(kind=8),intent(inout) :: trH_old
  type(DFT_PSP_projectors),intent(inout) :: nlpsp
  integer,intent(in) :: scf_mode
  type(localizedDIISParameters),intent(inout) :: ldiis
  type(DFT_wavefunction),target,intent(inout) :: tmb
  type(SIC_data) :: SIC !<parameters for the SIC methods
  type(energy_terms),intent(in) :: energs_base
  integer, intent(in) :: nit_precond, target_function, correction_orthoconstraint, nit_basis
  real(kind=8),intent(out) :: ratio_deltas
  logical, intent(inout) :: ortho_on
  integer, intent(in) :: extra_states
  integer,intent(in) :: itout
  real(kind=8),intent(in) :: conv_crit, early_stop, gnrm_dynamic, min_gnrm_for_dynamic, kappa_conv
  logical,intent(in) :: experimental_mode, purification_quickreturn
  logical,intent(out) :: can_use_ham
  integer,intent(in) :: method_updatekernel
  logical,intent(in) :: correction_co_contra
 
  ! Local variables
  real(kind=8) :: fnrmMax, meanAlpha, ediff_best, alpha_max, delta_energy, delta_energy_prev, ediff
  integer :: iorb, it, it_tot, ncount, jorb, ncharge
  real(kind=8),dimension(:),allocatable :: alpha,fnrmOldArr,alphaDIIS, hpsit_c_tmp, hpsit_f_tmp, hpsi_noconf, psidiff
  real(kind=8),dimension(:),allocatable :: delta_energy_arr
  real(kind=8),dimension(:),allocatable :: hpsi_noprecond, occup_tmp, kernel_compr_tmp, philarge
  logical :: energy_increased, overlap_calculated
  real(kind=8),dimension(:),pointer :: lhphiold, lphiold, hpsit_c, hpsit_f, hpsi_small
  type(energy_terms) :: energs
  real(kind=8), dimension(2):: reducearr
  real(gp) :: econf, dynamic_convcrit, kappa_mean
  integer :: i, ist, iiorb, ilr, ii, kappa_satur
  real(kind=8) :: energy_first, hxh, hyh, hzh, trH_ref, charge
  real(kind=8),dimension(:),allocatable :: kernel_best
  integer ::  correction_orthoconstraint_local, npsidim_large, ists, istl, sdim, ldim, nspin, nit_exit
  logical :: energy_diff, energy_increased_previous, complete_reset, even
  real(kind=8),dimension(3),save :: kappa_history
  integer,save :: nkappa_history
  logical,save :: has_already_converged
  logical,dimension(7) :: exit_loop

  call f_routine(id='getLocalizedBasis')

  delta_energy_arr=f_malloc(nit_basis+6,id='delta_energy_arr')
  kernel_best=f_malloc(tmb%linmat%l%nvctr,id='kernel_best')
  energy_diff=.false.


  ! Allocate all local arrays.
  call allocateLocalArrays()


  call timing(iproc,'getlocbasinit','ON')
  tmb%can_use_transposed=.false.
  !!if(iproc==0) write(*,'(1x,a)') '======================== Creation of the basis functions... ========================'

  alpha=ldiis%alphaSD
  alphaDIIS=ldiis%alphaDIIS
  ldiis%resetDIIS=.false.
  ldiis%immediateSwitchToSD=.false.
 
  call timing(iproc,'getlocbasinit','OF')

  overlap_calculated=.false.
  it=0
  it_tot=0
  !ortho=.true.
  call local_potential_dimensions(iproc,tmb%ham_descr%lzd,tmb%orbs,denspot%xc,denspot%dpbox%ngatherarr(0,1))
  call start_onesided_communication(iproc, nproc, max(denspot%dpbox%ndimpot,1), denspot%rhov, &
       tmb%ham_descr%comgp%nrecvbuf, tmb%ham_descr%comgp%recvbuf, tmb%ham_descr%comgp, tmb%ham_descr%lzd)

  delta_energy_prev=1.d100

  energy_increased_previous=.false.
  ratio_deltas=1.d0
  ediff_best=1.d0
  ediff=1.d0
  delta_energy_prev=1.d0
  delta_energy_arr=1.d0
  trH_ref=trH_old
  dynamic_convcrit=1.d-100
  kappa_satur=0


  ! Count whether there is an even or an odd number of electrons
  charge=0.d0
  do iorb=1,orbs%norb
      charge=charge+orbs%occup(iorb)
  end do
  ncharge=nint(charge)
  even=(mod(ncharge,2)==0)

  ! Purify the initial kernel (only when necessary and if there is an even number of electrons)
  if (target_function/=TARGET_FUNCTION_IS_TRACE .and. even .and. scf_mode==LINEAR_FOE) then
      if (iproc==0) then
          call yaml_sequence(advance='no')
          call yaml_mapping_open(flow=.true.)
          call yaml_map('Initial kernel purification',.true.)
      end if
      overlap_calculated=.true.
      call purify_kernel(iproc, nproc, tmb, overlap_calculated, 1, 30, order_taylor, purification_quickreturn)
      if (iproc==0) call yaml_mapping_close()
  end if

  if (itout==0) then
      nkappa_history=0
      kappa_history=0.d0
      has_already_converged=.false.
  end if

  iterLoop: do
      it=it+1
      it=max(it,1) !since it could become negative (2 is subtracted if the loop cycles)
      it_tot=it_tot+1

      fnrmMax=0.d0
      fnrm=0.d0
  
      if (iproc==0) then
          call yaml_sequence(advance='no')
          call yaml_mapping_open(flow=.true.)
          call yaml_comment('iter:'//yaml_toa(it,fmt='(i6)'),hfill='-')
          if (target_function==TARGET_FUNCTION_IS_TRACE) then
              call yaml_map('target function','TRACE')
          else if (target_function==TARGET_FUNCTION_IS_ENERGY) then
              call yaml_map('target function','ENERGY')
          else if (target_function==TARGET_FUNCTION_IS_HYBRID) then
              call yaml_map('target function','HYBRID')
          end if
      end if


      ! Calculate the unconstrained gradient by applying the Hamiltonian.
      !!if (iproc==0) write(*,*) 'tmb%psi(1)',tmb%psi(1)
      if (tmb%ham_descr%npsidim_orbs > 0)  call to_zero(tmb%ham_descr%npsidim_orbs,tmb%hpsi(1))
      call small_to_large_locreg(iproc, tmb%npsidim_orbs, tmb%ham_descr%npsidim_orbs, tmb%lzd, tmb%ham_descr%lzd, &
           tmb%orbs, tmb%psi, tmb%ham_descr%psi)

      call NonLocalHamiltonianApplication(iproc,at,tmb%ham_descr%npsidim_orbs,tmb%orbs,rxyz,&
           tmb%ham_descr%lzd,nlpsp,tmb%ham_descr%psi,tmb%hpsi,energs%eproj)
      ! only kinetic because waiting for communications
      call LocalHamiltonianApplication(iproc,nproc,at,tmb%ham_descr%npsidim_orbs,tmb%orbs,&
           tmb%ham_descr%lzd,tmb%confdatarr,denspot%dpbox%ngatherarr,denspot%pot_work,&
           & tmb%ham_descr%psi,tmb%hpsi,energs,SIC,GPU,3,denspot%xc,&
           & pkernel=denspot%pkernelseq,dpbox=denspot%dpbox,&
           & potential=denspot%rhov,comgp=tmb%ham_descr%comgp)
      call full_local_potential(iproc,nproc,tmb%orbs,tmb%ham_descr%lzd,2,denspot%dpbox,&
           & denspot%xc,denspot%rhov,denspot%pot_work,tmb%ham_descr%comgp)
      ! only potential
      if (target_function==TARGET_FUNCTION_IS_HYBRID) then
          call vcopy(tmb%ham_descr%npsidim_orbs, tmb%hpsi(1), 1, hpsi_noconf(1), 1)
          call LocalHamiltonianApplication(iproc,nproc,at,tmb%ham_descr%npsidim_orbs,tmb%orbs,&
               tmb%ham_descr%lzd,tmb%confdatarr,denspot%dpbox%ngatherarr,denspot%pot_work,&
               & tmb%ham_descr%psi,tmb%hpsi,energs,SIC,GPU,2,denspot%xc,&
               & pkernel=denspot%pkernelseq,dpbox=denspot%dpbox,&
               & potential=denspot%rhov,comgp=tmb%ham_descr%comgp,&
               hpsi_noconf=hpsi_noconf,econf=econf)

          if (nproc>1) then
              call mpiallred(econf, 1, mpi_sum, bigdft_mpi%mpi_comm)
          end if

      else
          call LocalHamiltonianApplication(iproc,nproc,at,tmb%ham_descr%npsidim_orbs,tmb%orbs,&
               tmb%ham_descr%lzd,tmb%confdatarr,denspot%dpbox%ngatherarr,&
               & denspot%pot_work,tmb%ham_descr%psi,tmb%hpsi,energs,SIC,GPU,2,denspot%xc,&
               & pkernel=denspot%pkernelseq,dpbox=denspot%dpbox,&
               & potential=denspot%rhov,comgp=tmb%ham_descr%comgp)
      end if


      !!if (target_function==TARGET_FUNCTION_IS_HYBRID .and. iproc==0) then
      !!    write(*,*) 'econf, econf/tmb%orbs%norb',econf, econf/tmb%orbs%norb
      !!end if

      call timing(iproc,'glsynchham2','ON')
      call SynchronizeHamiltonianApplication(nproc,tmb%ham_descr%npsidim_orbs,tmb%orbs,tmb%ham_descr%lzd,GPU,denspot%xc,tmb%hpsi,&
           energs%ekin,energs%epot,energs%eproj,energs%evsic,energs%eexctX)
      call timing(iproc,'glsynchham2','OF')

      if (iproc==0) then
          call yaml_map('Hamiltonian Applied',.true.)
      end if

      ! Use this subroutine to write the energies, with some fake number
      ! to prevent it from writing too much
      if (iproc==0) then
          call write_energies(0,0,energs,0.d0,0.d0,'',.true.)
      end if

      !if (iproc==0) write(*,'(a,5es16.6)') 'ekin, eh, epot, eproj, eex', &
      !              energs%ekin, energs%eh, energs%epot, energs%eproj, energs%exc

      if (iproc==0) then
          call yaml_map('Orthoconstraint',.true.)
      end if


      if (target_function==TARGET_FUNCTION_IS_HYBRID) then
          call transpose_localized(iproc, nproc, tmb%ham_descr%npsidim_orbs, tmb%orbs, tmb%ham_descr%collcom, &
               hpsi_noconf, hpsit_c, hpsit_f, tmb%ham_descr%lzd)
          if (method_updatekernel==UPDATE_BY_FOE) then
              !@NEW
              !!if(associated(tmb%ham_descr%psit_c)) then
              !!    call f_free_ptr(tmb%ham_descr%psit_c)
              !!    associated_psitlarge_c=.true.
              !!else
              !!    associated_psitlarge_c=.false.
              !!end if
              !!if(associated(tmb%ham_descr%psit_f)) then
              !!    call f_free_ptr(tmb%ham_descr%psit_f)
              !!    associated_psitlarge_f=.true.
              !!else
              !!    associated_psitlarge_f=.false.
              !!end if

              !!tmb%ham_descr%psit_c = f_malloc_ptr(tmb%ham_descr%collcom%ndimind_c,id='tmb%ham_descr%psit_c')
              !!tmb%ham_descr%psit_f = f_malloc_ptr(7*tmb%ham_descr%collcom%ndimind_f,id='tmb%ham_descr%psit_f')
              call transpose_localized(iproc, nproc, tmb%ham_descr%npsidim_orbs, tmb%orbs, tmb%ham_descr%collcom, &
                   tmb%ham_descr%psi, tmb%ham_descr%psit_c, tmb%ham_descr%psit_f, tmb%ham_descr%lzd)
              call calculate_overlap_transposed(iproc, nproc, tmb%orbs, tmb%ham_descr%collcom, &
                   tmb%ham_descr%psit_c, hpsit_c, tmb%ham_descr%psit_f, hpsit_f, tmb%linmat%m, tmb%linmat%ham_)

              !!if (.not.associated(tmb%psit_c)) then
              !!    tmb%psit_c = f_malloc_ptr(tmb%collcom%ndimind_c,id='tmb%psit_c')
              !!end if
              !!if (.not.associated(tmb%psit_f)) then
              !!    tmb%psit_f = f_malloc_ptr(7*tmb%collcom%ndimind_f,id='tmb%psit_f')
              !!end if
              call transpose_localized(iproc, nproc, tmb%npsidim_orbs, tmb%orbs, tmb%collcom, &
                   tmb%psi, tmb%psit_c, tmb%psit_f, tmb%lzd)
              call calculate_overlap_transposed(iproc, nproc, tmb%orbs, tmb%collcom, &
                   tmb%psit_c, tmb%psit_c, tmb%psit_f, tmb%psit_f, tmb%linmat%s, tmb%linmat%ovrlp_)
              if (iproc==0) call yaml_newline()
              if (iproc==0) call yaml_sequence_open('kernel update by FOE')
              call foe(iproc, nproc, 0.d0, &
                   energs%ebs, -1, -10, order_taylor, purification_quickreturn, 0, &
                   FOE_FAST, tmb, tmb%foe_obj)
              if (iproc==0) call yaml_sequence_close()
              !if (.not.associated_psitlarge_c) then
              !    call f_free_ptr(tmb%ham_descr%psit_c)
              !end if
              !if (.not.associated_psitlarge_f) then
              !    call f_free_ptr(tmb%ham_descr%psit_f)
              !end if
              !if (associated_psitlarge_c .and. associated_psitlarge_f) then
                  tmb%ham_descr%can_use_transposed=.true.
              !end if
              !@ENDNEW
          end if
      else
          call transpose_localized(iproc, nproc, tmb%ham_descr%npsidim_orbs, tmb%orbs, tmb%ham_descr%collcom, &
               tmb%hpsi, hpsit_c, hpsit_f, tmb%ham_descr%lzd)
      end if

      ncount=sum(tmb%ham_descr%collcom%nrecvcounts_c)
      if(ncount>0) call vcopy(ncount, hpsit_c(1), 1, hpsit_c_tmp(1), 1)
      ncount=7*sum(tmb%ham_descr%collcom%nrecvcounts_f)
      if(ncount>0) call vcopy(ncount, hpsit_f(1), 1, hpsit_f_tmp(1), 1)

      ! optimize the tmbs for a few extra states
      if (target_function==TARGET_FUNCTION_IS_ENERGY.and.extra_states>0) then
          kernel_compr_tmp = sparsematrix_malloc(tmb%linmat%l, iaction=SPARSE_FULL, id='kernel_compr_tmp')
          call vcopy(tmb%linmat%l%nvctr, tmb%linmat%kernel_%matrix_compr(1), 1, kernel_compr_tmp(1), 1)
          !allocate(occup_tmp(tmb%orbs%norb), stat=istat)
          !call memocc(istat, occup_tmp, 'occup_tmp', subname)
          !call vcopy(tmb%orbs%norb, tmb%orbs%occup(1), 1, occup_tmp(1), 1)
          !call to_zero(tmb%orbs%norb,tmb%orbs%occup(1))
          !call vcopy(orbs%norb, orbs%occup(1), 1, tmb%orbs%occup(1), 1)
          !! occupy the next few states - don't need to preserve the charge as only using for support function optimization
          !do iorb=1,tmb%orbs%norb
          !   if (tmb%orbs%occup(iorb)==1.0_gp) then
          !      tmb%orbs%occup(iorb)=2.0_gp
          !   else if (tmb%orbs%occup(iorb)==0.0_gp) then
          !      do jorb=iorb,min(iorb+extra_states-1,tmb%orbs%norb)
          !         tmb%orbs%occup(jorb)=2.0_gp
          !      end do
          !      exit
          !   end if
          !end do
          call calculate_density_kernel(iproc, nproc, .true., tmb%orbs, tmb%orbs, tmb%coeff, &
               tmb%linmat%l, tmb%linmat%kernel_)
          !tmb%linmat%denskern_large%matrix_compr = tmb%linmat%kernel_%matrix_compr
          !call transform_sparse_matrix(tmb%linmat%denskern, tmb%linmat%denskern_large, 'large_to_small')
      end if

      correction_orthoconstraint_local=correction_orthoconstraint
      !if (target_function==TARGET_FUNCTION_IS_HYBRID) then
      !    correction_orthoconstraint_local=2
      !end if
      !if(.not.ortho_on) then
      !    correction_orthoconstraint_local=2
      !end if
      !write(*,*) 'correction_orthoconstraint, correction_orthoconstraint_local',correction_orthoconstraint, correction_orthoconstraint_local


      !!! PLOT ###########################################################################
      !!hxh=0.5d0*tmb%lzd%hgrids(1)      
      !!hyh=0.5d0*tmb%lzd%hgrids(2)      
      !!hzh=0.5d0*tmb%lzd%hgrids(3)      
      !!npsidim_large=tmb%lzd%glr%wfd%nvctr_c+7*tmb%lzd%glr%wfd%nvctr_f
      !!allocate(philarge((tmb%lzd%glr%wfd%nvctr_c+7*tmb%lzd%glr%wfd%nvctr_f)*tmb%orbs%norbp))
      !!philarge=0.d0
      !!ists=1
      !!istl=1
      !!do iorb=1,tmb%orbs%norbp
      !!    ilr = tmb%orbs%inWhichLocreg(tmb%orbs%isorb+iorb)
      !!    sdim=tmb%lzd%llr(ilr)%wfd%nvctr_c+7*tmb%lzd%llr(ilr)%wfd%nvctr_f
      !!    ldim=tmb%lzd%glr%wfd%nvctr_c+7*tmb%lzd%glr%wfd%nvctr_f
      !!    nspin=1 !this must be modified later
      !!    call Lpsi_to_global2(iproc, sdim, ldim, tmb%orbs%norb, tmb%orbs%nspinor, nspin, tmb%lzd%glr, &
      !!         tmb%lzd%llr(ilr), tmb%psi(ists), philarge(istl))
      !!    ists=ists+tmb%lzd%llr(ilr)%wfd%nvctr_c+7*tmb%lzd%llr(ilr)%wfd%nvctr_f
      !!    istl=istl+tmb%lzd%glr%wfd%nvctr_c+7*tmb%lzd%glr%wfd%nvctr_f
      !!end do
      !!call plotOrbitals(iproc, tmb, philarge, at%astruct%nat, rxyz, hxh, hyh, hzh, 100*itout+it, 'orbs')
      !!deallocate(philarge)
      !!! END PLOT #######################################################################


      !if (iproc==0) write(*,*) 'tmb%linmat%denskern%matrix_compr(1)',tmb%linmat%denskern%matrix_compr(1)
      call calculate_energy_and_gradient_linear(iproc, nproc, it, ldiis, fnrmOldArr, alpha, trH, trH_old, fnrm, fnrmMax, &
           meanAlpha, alpha_max, energy_increased, tmb, lhphiold, overlap_calculated, energs_base, &
           hpsit_c, hpsit_f, nit_precond, target_function, correction_orthoconstraint_local, hpsi_small, &
           experimental_mode, correction_co_contra, hpsi_noprecond, order_taylor, method_updatekernel)


      !!! PLOT ###########################################################################
      !!hxh=0.5d0*tmb%lzd%hgrids(1)      
      !!hyh=0.5d0*tmb%lzd%hgrids(2)      
      !!hzh=0.5d0*tmb%lzd%hgrids(3)      
      !!npsidim_large=tmb%lzd%glr%wfd%nvctr_c+7*tmb%lzd%glr%wfd%nvctr_f
      !!allocate(philarge((tmb%lzd%glr%wfd%nvctr_c+7*tmb%lzd%glr%wfd%nvctr_f)*tmb%orbs%norbp))
      !!philarge=0.d0
      !!ists=1
      !!istl=1
      !!do iorb=1,tmb%orbs%norbp
      !!    ilr = tmb%orbs%inWhichLocreg(tmb%orbs%isorb+iorb)
      !!    sdim=tmb%lzd%llr(ilr)%wfd%nvctr_c+7*tmb%lzd%llr(ilr)%wfd%nvctr_f
      !!    ldim=tmb%lzd%glr%wfd%nvctr_c+7*tmb%lzd%glr%wfd%nvctr_f
      !!    nspin=1 !this must be modified later
      !!    !call Lpsi_to_global2(iproc, sdim, ldim, tmb%orbs%norb, tmb%orbs%nspinor, nspin, tmb%lzd%glr, &
      !!    !      tmb%lzd%llr(ilr), hpsi_small(ists), philarge(istl))
      !!    call Lpsi_to_global2(iproc, sdim, ldim, tmb%orbs%norb, tmb%orbs%nspinor, nspin, tmb%lzd%glr, &
      !!          tmb%lzd%llr(ilr), tmb%psi(ists), philarge(istl))
      !!    ists=ists+tmb%lzd%llr(ilr)%wfd%nvctr_c+7*tmb%lzd%llr(ilr)%wfd%nvctr_f
      !!    istl=istl+tmb%lzd%glr%wfd%nvctr_c+7*tmb%lzd%glr%wfd%nvctr_f
      !!end do
      !!!call plotOrbitals(iproc, tmb, philarge, at%astruct%nat, rxyz, hxh, hyh, hzh, 100*itout+it, 'grad')
      !!call plotOrbitals(iproc, tmb, philarge, at%astruct%nat, rxyz, hxh, hyh, hzh, 100*itout+it, 'tmbs')
      !!deallocate(philarge)
      !!! END PLOT #######################################################################



      if (experimental_mode) then
          if (it_tot==1) then
              energy_first=trH
          end if
          !!if (iproc==0) write(*,'(a,3es16.7)') 'trH, energy_first, (trH-energy_first)/energy_first', &
          !!                                      trH, energy_first, (trH-energy_first)/energy_first
          if (iproc==0) call yaml_map('rel D',(trH-energy_first)/energy_first,fmt='(es9.2)')
          if ((trH-energy_first)/energy_first>early_stop .and. itout>0) then
              energy_diff=.true.
              !!if (iproc==0) write(*,'(a,3es16.7)') 'new stopping crit: trH, energy_first, (trH-energy_first)/energy_first', &
              !!                                      trH, energy_first, (trH-energy_first)/energy_first
          end if
      end if

      if (target_function==TARGET_FUNCTION_IS_ENERGY.and.extra_states>0) then
          call vcopy(tmb%linmat%l%nvctr, kernel_compr_tmp(1), 1, tmb%linmat%kernel_%matrix_compr(1), 1)
          call f_free(kernel_compr_tmp)
      end if

      ediff=trH-trH_old
      ediff_best=trH-trH_ref
      !!if (iproc==0) write(*,*) 'trH, trH_ref', trH, trH_ref

      if (it>1 .and. (target_function==TARGET_FUNCTION_IS_HYBRID .or. experimental_mode)) then
          if (.not.energy_increased .and. .not.energy_increased_previous) then
              if (.not.ldiis%switchSD) then
                  ratio_deltas=ediff_best/delta_energy_prev
              else
                  ratio_deltas=ediff_best/delta_energy_arr(ldiis%itBest)
              end if
          else
              ! use a default value
              if (iproc==0) then
                  call yaml_warning('use a fake value for kappa')
                  call yaml_newline()
              end if
              ratio_deltas=0.5d0
          end if
          if (ldiis%switchSD) then
              !!ratio_deltas=0.5d0
              !!if (iproc==0) write(*,*) 'WARNING: TEMPORARY FIX for ratio_deltas!'
          end if
          if (iproc==0) call yaml_map('kappa',ratio_deltas,fmt='(es10.3)')
          if (target_function==TARGET_FUNCTION_IS_HYBRID) then
              !if (ratio_deltas>0.d0) then
              if (ratio_deltas>1.d-12) then
                  if (iproc==0) call yaml_map('kappa to history',.true.)
                  nkappa_history=nkappa_history+1
                  ii=mod(nkappa_history-1,3)+1
                  kappa_history(ii)=ratio_deltas
              end if
              !!if (nkappa_history>=3) then
              !!    kappa_mean=sum(kappa_history)/3.d0
              !!    if (iproc==0) call yaml_map('mean kappa',kappa_mean,fmt='(es10.3)')
              !!    dynamic_convcrit=conv_crit/kappa_mean
              !!    if (iproc==0) call yaml_map('dynamic conv crit',dynamic_convcrit,fmt='(es9.2)')
              !!end if
          end if
      end if
      if (target_function==TARGET_FUNCTION_IS_HYBRID) then
          if (nkappa_history>=3) then
              kappa_mean=sum(kappa_history)/3.d0
              if (iproc==0) call yaml_map('mean kappa',kappa_mean,fmt='(es10.3)')
              !dynamic_convcrit=conv_crit/kappa_mean
              dynamic_convcrit=gnrm_dynamic/kappa_mean
              if (iproc==0) call yaml_map('dynamic conv crit',dynamic_convcrit,fmt='(es9.2)')
          end if
      end if

      if (energy_increased) then
          energy_increased_previous=.true.
      else
          energy_increased_previous=.false.
      end if



      !!delta_energy_prev=delta_energy

      if (energy_increased) then
          !if (iproc==0) write(*,*) 'WARNING: ENERGY INCREASED'
          !if (iproc==0) call yaml_warning('The target function increased, D='&
          !              //trim(adjustl(yaml_toa(trH-ldiis%trmin,fmt='(es10.3)'))))
          if (iproc==0) then
              call yaml_newline()
              call yaml_map('iter',it,fmt='(i5)')
              call yaml_map('fnrm',fnrm,fmt='(es9.2)')
              call yaml_map('Omega',trH,fmt='(es22.15)')
              call yaml_map('D',ediff,fmt='(es9.2)')
              call yaml_map('D best',ediff_best,fmt='(es9.2)')
          end if
          tmb%ham_descr%can_use_transposed=.false.
          call vcopy(tmb%npsidim_orbs, lphiold(1), 1, tmb%psi(1), 1)
          can_use_ham=.false.
          call vcopy(tmb%linmat%l%nvctr, kernel_best(1), 1, tmb%linmat%kernel_%matrix_compr(1), 1)
          trH_old=0.d0
          it=it-2 !go back one iteration (minus 2 since the counter was increased)
          !if(associated(tmb%ham_descr%psit_c)) then
          !    call f_free_ptr(tmb%ham_descr%psit_c)
          !end if
          !if(associated(tmb%ham_descr%psit_f)) then
          !    call f_free_ptr(tmb%ham_descr%psit_f)
          !end if
          !!if(iproc==0) write(*,*) 'it_tot',it_tot
          overlap_calculated=.false.
          ! print info here anyway for debugging
          if (it_tot<2*nit_basis) then ! just in case the step size is the problem
              call yaml_mapping_close()
              call bigdft_utils_flush(unit=6)
             cycle
          else if(it_tot<3*nit_basis) then ! stop orthonormalizing the tmbs
             if (iproc==0) call yaml_newline()
             if (iproc==0) call yaml_warning('Energy increasing, switching off orthonormalization of tmbs')
             ortho_on=.false.
             alpha=alpha*5.0d0/3.0d0 ! increase alpha to make up for decrease from previous iteration
          end if
      else
          can_use_ham=.true.
      end if 


      ! information on the progress of the optimization
      if (iproc==0) then
          call yaml_newline()
          call yaml_map('iter',it,fmt='(i5)')
          call yaml_map('fnrm',fnrm,fmt='(es9.2)')
          call yaml_map('Omega',trH,fmt='(es22.15)')
          call yaml_map('D',ediff,fmt='(es9.2)')
          call yaml_map('D best',ediff_best,fmt='(es9.2)')
      end if

      ! Add some extra iterations if DIIS failed (max 6 failures are allowed before switching to SD)
      nit_exit=min(nit_basis+ldiis%icountDIISFailureTot,nit_basis+6)

      ! Determine whether the loop should be exited
      exit_loop(1) = (it>=nit_exit)
      exit_loop(2) = (it_tot>=3*nit_basis)
      exit_loop(3) = energy_diff
      exit_loop(4) = (fnrm<conv_crit .and. experimental_mode)
      exit_loop(5) = (experimental_mode .and. fnrm<dynamic_convcrit .and. fnrm<min_gnrm_for_dynamic &
                     .and. (it>1 .or. has_already_converged)) ! first overall convergence not allowed in a first iteration
      exit_loop(6) = (itout==0 .and. it>1 .and. ratio_deltas<kappa_conv .and.  ratio_deltas>0.d0)
      if (ratio_deltas>0.d0 .and. ratio_deltas<1.d-1) then
          kappa_satur=kappa_satur+1
      else
          kappa_satur=0
      end if
      exit_loop(7) = (.false. .and. itout>0 .and. kappa_satur>=2)

      if(any(exit_loop)) then
          if(exit_loop(1)) then
              infoBasisFunctions=-1
              if(iproc==0) call yaml_map('exit criterion','net number of iterations')
          end if
          if (exit_loop(2)) then
              infoBasisFunctions=-2
              if (iproc==0) call yaml_map('exit criterion','total number of iterations')
          end if
          if (exit_loop(3)) then
              infoBasisFunctions=it
              if (iproc==0) call yaml_map('exit criterion','energy difference')
          end if
          if (exit_loop(4)) then
              if (iproc==0) call yaml_map('exit criterion','gradient')
              infoBasisFunctions=it
          end if
          if (exit_loop(5)) then
              if (iproc==0) call yaml_map('exit criterion','dynamic gradient')
              infoBasisFunctions=it
              has_already_converged=.true.
          end if
          if (exit_loop(6)) then
              infoBasisFunctions=it
              if (iproc==0) call yaml_map('exit criterion','extended input guess')
          end if
          if (exit_loop(7)) then
              infoBasisFunctions=it
              if (iproc==0) call yaml_map('exit criterion','kappa')
          end if
          if (can_use_ham) then
              ! Calculate the Hamiltonian matrix, since we have all quantities ready. This matrix can then be used in the first
              ! iteration of get_coeff.
              call calculate_overlap_transposed(iproc, nproc, tmb%orbs, tmb%ham_descr%collcom, &
                   tmb%ham_descr%psit_c, hpsit_c_tmp, tmb%ham_descr%psit_f, hpsit_f_tmp, tmb%linmat%m, tmb%linmat%ham_)
              ! This can then be deleted if the transition to the new type has been completed.
              !tmb%linmat%ham%matrix_compr=tmb%linmat%ham_%matrix_compr
          end if

          if (iproc==0) then
              !yaml output
              call yaml_mapping_close() !iteration
              call bigdft_utils_flush(unit=6)
          end if

          exit iterLoop
      end if
      trH_old=trH

      if (ldiis%isx>0) then
          ldiis%mis=mod(ldiis%is,ldiis%isx)+1 !to store the energy at the correct location in the history
      end if
      call hpsitopsi_linear(iproc, nproc, it, ldiis, tmb, &
           lphiold, alpha, trH, meanAlpha, alpha_max, alphaDIIS, hpsi_small, ortho_on, psidiff, &
           experimental_mode, order_taylor, trH_ref, kernel_best, complete_reset)
      !if (iproc==0) write(*,*) 'kernel_best(1)',kernel_best(1)
      !if (iproc==0) write(*,*) 'tmb%linmat%denskern%matrix_compr(1)',tmb%linmat%denskern%matrix_compr(1)


      overlap_calculated=.false.
      ! It is now not possible to use the transposed quantities, since they have changed.
      if(tmb%ham_descr%can_use_transposed) then
          !call f_free_ptr(tmb%ham_descr%psit_c)
          !call f_free_ptr(tmb%ham_descr%psit_f)
          tmb%ham_descr%can_use_transposed=.false.
      end if

      ! Estimate the energy change, that is to be expected in the next optimization
      ! step, given by the product of the force and the "displacement" .
      if (target_function==TARGET_FUNCTION_IS_HYBRID .or. experimental_mode) then
          call estimate_energy_change(tmb%npsidim_orbs, tmb%orbs, tmb%lzd, psidiff, hpsi_noprecond, delta_energy)
          ! This is a hack...
          if (energy_increased) then
              delta_energy=1.d100
              !ratio_deltas=1.d100
          end if
          !if (iproc==0) write(*,*) 'delta_energy', delta_energy
          delta_energy_prev=delta_energy
          delta_energy_arr(max(it,1))=delta_energy !max since the counter was decreased if there are problems, might lead to wrong results otherwise
      end if


      ! Only need to reconstruct the kernel if it is actually used.
      if ((target_function/=TARGET_FUNCTION_IS_TRACE .or. scf_mode==LINEAR_DIRECT_MINIMIZATION) &
           .and. .not.complete_reset ) then
          if(scf_mode/=LINEAR_FOE) then
              call reconstruct_kernel(iproc, nproc, order_taylor, tmb%orthpar%blocksize_pdsyev, &
                   tmb%orthpar%blocksize_pdgemm, orbs, tmb, overlap_calculated)
              if (iproc==0) call yaml_map('reconstruct kernel',.true.)
          else if (experimental_mode .and. .not.complete_reset) then
              if (method_updatekernel==UPDATE_BY_PURIFICATION) then
                  if (iproc==0) then
                      call yaml_map('purify kernel',.true.)
                      call yaml_newline()
                  end if
                  call purify_kernel(iproc, nproc, tmb, overlap_calculated, 1, 30, order_taylor, purification_quickreturn)
                  !tmb%linmat%denskern_large%matrix_compr = tmb%linmat%kernel_%matrix_compr
              else if (method_updatekernel==UPDATE_BY_FOE) then
                  if (iproc==0) then
                      call yaml_map('purify kernel',.false.)
                  end if
              end if
          end if
      end if

      if (iproc==0) then
          !yaml output
          call yaml_mapping_close() !iteration
          call bigdft_utils_flush(unit=6)
      end if


  end do iterLoop

  ! Write the final results
  if (iproc==0) then
      call yaml_sequence(label='final_supfun'//trim(adjustl(yaml_toa(itout,fmt='(i3.3)'))),advance='no')
      call yaml_mapping_open(flow=.true.)
      call yaml_comment('iter:'//yaml_toa(it,fmt='(i6)'),hfill='-')
      if (target_function==TARGET_FUNCTION_IS_TRACE) then
          call yaml_map('target function','TRACE')
      else if (target_function==TARGET_FUNCTION_IS_ENERGY) then
          call yaml_map('target function','ENERGY')
      else if (target_function==TARGET_FUNCTION_IS_HYBRID) then
          call yaml_map('target function','HYBRID')
      end if
      call write_energies(0,0,energs,0.d0,0.d0,'',.true.)
      call yaml_newline()
      call yaml_map('iter',it,fmt='(i5)')
      call yaml_map('fnrm',fnrm,fmt='(es9.2)')
      call yaml_map('Omega',trH,fmt='(es22.15)')
      call yaml_map('D',ediff,fmt='(es9.2)')
      call yaml_map('D best',ediff_best,fmt='(es9.2)')
      call yaml_mapping_close() !iteration
      call bigdft_utils_flush(unit=6)
  end if


  if (iproc==0) then
      call yaml_comment('Support functions created')
  end if


  ! Deallocate potential
  call f_free_ptr(denspot%pot_work)


  ! Keep the values for the next iteration
  reducearr(1)=0.d0
  reducearr(2)=0.d0
  do iorb=1,tmb%orbs%norbp
      reducearr(1)=reducearr(1)+alpha(iorb)
      reducearr(2)=reducearr(2)+alphaDIIS(iorb)
  end do

  if (nproc > 1) then
      call mpiallred(reducearr(1), 2, mpi_sum, bigdft_mpi%mpi_comm)
  end if

  reducearr(1)=reducearr(1)/dble(tmb%orbs%norb)
  reducearr(2)=reducearr(2)/dble(tmb%orbs%norb)

  ldiis%alphaSD=reducearr(1)
  ldiis%alphaDIIS=reducearr(2)


  ! Deallocate all local arrays.
  call deallocateLocalArrays()
  call f_free(delta_energy_arr)
  call f_free(kernel_best)

  call f_release_routine()

contains


    subroutine allocateLocalArrays()
    !
    ! Purpose:
    ! ========
    !   This subroutine allocates all local arrays.
    !
      alpha = f_malloc(tmb%orbs%norbp,id='alpha')
      alphaDIIS = f_malloc(tmb%orbs%norbp,id='alphaDIIS')
      fnrmOldArr = f_malloc(tmb%orbs%norb,id='fnrmOldArr')
      hpsi_small = f_malloc_ptr(max(tmb%npsidim_orbs, tmb%npsidim_comp),id='hpsi_small')
      lhphiold = f_malloc_ptr(max(tmb%npsidim_orbs, tmb%npsidim_comp),id='lhphiold')
      lphiold = f_malloc_ptr(size(tmb%psi),id='lphiold')
      hpsit_c = f_malloc_ptr(sum(tmb%ham_descr%collcom%nrecvcounts_c),id='hpsit_c')
      hpsit_f = f_malloc_ptr(7*sum(tmb%ham_descr%collcom%nrecvcounts_f),id='hpsit_f')
      hpsit_c_tmp = f_malloc(sum(tmb%ham_descr%collcom%nrecvcounts_c),id='hpsit_c_tmp')
      hpsit_f_tmp = f_malloc(7*sum(tmb%ham_descr%collcom%nrecvcounts_f),id='hpsit_f_tmp')
      hpsi_noconf = f_malloc(tmb%ham_descr%npsidim_orbs,id='hpsi_noconf')
      psidiff = f_malloc(tmb%npsidim_orbs,id='psidiff')
      hpsi_noprecond = f_malloc(tmb%npsidim_orbs,id='hpsi_noprecond')

    end subroutine allocateLocalArrays


    subroutine deallocateLocalArrays()
    !
    ! Purpose:
    ! ========
    !   This subroutine deallocates all local arrays.
    !
    call f_free(alpha)
    call f_free(alphaDIIS)
    call f_free(fnrmOldArr)
    call f_free_ptr(hpsi_small)
    call f_free_ptr(lhphiold)
    call f_free_ptr(lphiold)
    call f_free_ptr(hpsit_c)
    call f_free_ptr(hpsit_f)
    call f_free(hpsit_c_tmp)
    call f_free(hpsit_f_tmp)
    call f_free(hpsi_noconf)
    call f_free(psidiff)
    call f_free(hpsi_noprecond)

    end subroutine deallocateLocalArrays


end subroutine getLocalizedBasis


subroutine purify_kernel(iproc, nproc, tmb, overlap_calculated, it_shift, it_opt, order_taylor, &
           max_inversion_error, purification_quickreturn, ispin)
  use module_base
  use module_types
  use yaml_output
  use communications_base, only: TRANSPOSE_FULL
  use communications, only: transpose_localized
  use sparsematrix_base, only: sparsematrix_malloc_ptr, DENSE_FULL, assignment(=), matrices, &
                               matrices_null, allocate_matrices, deallocate_matrices
  use sparsematrix, only: uncompress_matrix, gather_matrix_from_taskgroups_inplace, &
                          uncompress_matrix2, compress_matrix2, trace_sparse
  use foe_base, only: foe_data_get_real
  use transposed_operations, only: calculate_overlap_transposed
  use matrix_operations, only: overlapPowerGeneral, check_taylor_order
  implicit none

  ! Calling arguments
  integer,intent(in) :: iproc, nproc
  integer,intent(inout) :: order_taylor
  real(kind=8),intent(in) :: max_inversion_error
  type(DFT_wavefunction),intent(inout):: tmb
  logical,intent(inout):: overlap_calculated
  integer,intent(in) :: it_shift, it_opt
  logical,intent(in) :: purification_quickreturn
  integer,intent(in) :: ispin

  ! Local variables
  integer :: it, iorb, jorb, jsegstart, jsegend, jseg, jjorb, iiorb !info, lwork, 
  integer :: ishift, isshift, ilshift
  real(kind=8) :: alpha, shift
  real(kind=8),dimension(:,:),allocatable :: ks, ksk, ksksk, kernel_prime
  !real(kind=8),dimension(:),allocatable :: eval, work
  character(len=*),parameter :: subname='purify_kernel'
  real(kind=8) :: diff, tr_KS, chargediff, max_error, mean_error
  !logical :: overlap_associated, inv_ovrlp_associated
  real(kind=8),dimension(2) :: bisec_bounds
  logical,dimension(2) :: bisec_bounds_ok
  !real(kind=8),dimension(:,:),pointer :: ovrlp_onehalf, ovrlp_minusonehalf
  type(matrices),dimension(1) :: ovrlppowers_
  type(matrices),dimension(1) :: ovrlp_onehalf_



  if (purification_quickreturn) then
      if (iproc==0) call yaml_warning('quick return in purification')
      if (iproc==0) call yaml_newline()
      return
  end if

  call f_routine(id='purify_kernel')

  isshift=(ispin-1)*tmb%linmat%s%nvctrp_tg
  ilshift=(ispin-1)*tmb%linmat%l%nvctrp_tg

  ovrlp_onehalf_(1) = matrices_null()
  call allocate_matrices(tmb%linmat%l, allocate_full=.true., matname='ovrlp_onehalf_', mat=ovrlp_onehalf_(1))
  ovrlppowers_(1) = matrices_null()
  call allocate_matrices(tmb%linmat%l, allocate_full=.true., matname='ovrlppowers_', mat=ovrlppowers_(1))


  ! Calculate the overlap matrix between the TMBs.
  if(.not. overlap_calculated) then
     if(.not.tmb%can_use_transposed) then
         !!if(associated(tmb%psit_c)) then
         !!    call f_free_ptr(tmb%psit_c)
         !!end if
         !!if(associated(tmb%psit_f)) then
         !!    call f_free_ptr(tmb%psit_f)
         !!end if
         !!tmb%psit_c = f_malloc_ptr(sum(tmb%collcom%nrecvcounts_c),id='tmb%psit_c')
         !!tmb%psit_f = f_malloc_ptr(7*sum(tmb%collcom%nrecvcounts_f),id='tmb%psit_f')
         call transpose_localized(iproc, nproc, tmb%npsidim_orbs, tmb%orbs, tmb%collcom, &
              TRANSPOSE_FULL, tmb%psi, tmb%psit_c, tmb%psit_f, tmb%lzd)
         tmb%can_use_transposed=.true.
     end if
     !call timing(iproc,'renormCoefComp','OF')

     call calculate_overlap_transposed(iproc, nproc, tmb%orbs, tmb%collcom, &
          tmb%psit_c, tmb%psit_c, tmb%psit_f, tmb%psit_f, tmb%linmat%s, tmb%linmat%ovrlp_)
     !!call gather_matrix_from_taskgroups_inplace(iproc, nproc, tmb%linmat%s, tmb%linmat%ovrlp_)

     !call timing(iproc,'renormCoefComp','ON')
     overlap_calculated=.true.
  end if


  tmb%linmat%ovrlp_%matrix = sparsematrix_malloc_ptr(tmb%linmat%s, iaction=DENSE_FULL, id='tmb%linmat%ovrlp_%matrix')
  tmb%linmat%kernel_%matrix = sparsematrix_malloc_ptr(tmb%linmat%l,iaction=DENSE_FULL,id='tmb%linmat%kernel_%matrix')
  call uncompress_matrix2(iproc, nproc, tmb%linmat%s, &
       tmb%linmat%ovrlp_%matrix_compr, tmb%linmat%ovrlp_%matrix)
  call uncompress_matrix2(iproc, nproc, tmb%linmat%l, &
       tmb%linmat%kernel_%matrix_compr, tmb%linmat%kernel_%matrix)

  ks=f_malloc((/tmb%linmat%l%nfvctr,tmb%linmat%l%nfvctr/),id='ks')
  ksk=f_malloc((/tmb%linmat%l%nfvctr,tmb%linmat%l%nfvctrp/),id='ksk')
  ksksk=f_malloc((/tmb%linmat%l%nfvctr,tmb%linmat%l%nfvctr/),id='ksksk')
  kernel_prime=f_malloc([tmb%linmat%l%nfvctr,tmb%linmat%l%nfvctr],id='kernel_prime')

  !ovrlp_onehalf=f_malloc_ptr((/tmb%linmat%l%nfvctr,tmb%linmat%l%nfvctr/),id='ovrlp_onehalf')
  !ovrlp_minusonehalf=f_malloc_ptr((/tmb%linmat%l%nfvctr,tmb%linmat%l%nfvctr/),id='ovrlp_minusonehalf')




  call timing(iproc,'purify_kernel ','ON') 

  if (tmb%linmat%l%nspin==1) then
      call dscal(tmb%linmat%l%nfvctr**2, 0.5d0, tmb%linmat%kernel_%matrix, 1)
  end if



  !!tmb%linmat%ovrlp_%matrix_compr = tmb%linmat%ovrlp%matrix_compr
  tr_KS=trace_sparse(iproc, nproc, tmb%linmat%s, tmb%linmat%l, &
        tmb%linmat%ovrlp_%matrix_compr(isshift+1:), &
        tmb%linmat%kernel_%matrix_compr(ilshift+1:), ispin)
  if (iproc==0) then
      call yaml_map('tr(KS) before purification',tr_KS)
      call yaml_newline
  end if


  alpha=1.d-4
  chargediff=0.d0
  
  !!if (.not.associated(tmb%linmat%inv_ovrlp_large%matrix_compr)) then
  !!    inv_ovrlp_associated=.false.
  !!    !!allocate(tmb%linmat%inv_ovrlp_large%matrix_compr(tmb%linmat%inv_ovrlp_large%nvctr),stat=istat)
  !!    !!call memocc(istat,tmb%linmat%inv_ovrlp_large%matrix_compr,'tmb%linmat%inv_ovrlp_large%matrix_compr',subname)
  !!    tmb%linmat%inv_ovrlp_large%matrix_compr=f_malloc_ptr(tmb%linmat%inv_ovrlp_large%nvctr,&
  !!        id='tmb%linmat%inv_ovrlp_large%matrix_compr')
  !!else
  !!    inv_ovrlp_associated=.true.
  !!end if

  if (it_shift>1) then
      call calculate_overlap_onehalf()
      call f_zero(kernel_prime)
      if (tmb%linmat%l%nfvctrp>0) then
          !SM: need to fix the spin here
          call dgemm('n', 'n', tmb%linmat%l%nfvctr, tmb%linmat%l%nfvctrp, tmb%linmat%l%nfvctr, &
                     1.d0, tmb%linmat%kernel_%matrix, tmb%linmat%l%nfvctr, &
                     ovrlp_onehalf_(1)%matrix(1,tmb%linmat%l%isfvctr+1,1), tmb%linmat%l%nfvctr, &
                     0.d0, ksksk, tmb%linmat%l%nfvctr) 
          call dgemm('n', 'n', tmb%linmat%l%nfvctr, tmb%linmat%l%nfvctrp, tmb%linmat%l%nfvctr, &
                     1.d0, ovrlp_onehalf_(1)%matrix, tmb%linmat%l%nfvctr, &
                     ksksk, tmb%linmat%l%nfvctr, &
                     0.d0, kernel_prime(1,tmb%linmat%l%isfvctr+1), tmb%linmat%l%nfvctr) 
      end if

      if (nproc > 1) then
          call mpiallred(kernel_prime, mpi_sum, comm=bigdft_mpi%mpi_comm)
      end if
  end if


  shift=0.d0
  bisec_bounds=0.d0
  bisec_bounds_ok=.false.


  shift_loop: do ishift=1,it_shift

  if (iproc==0) call yaml_newline()
  if (iproc==0) call yaml_map('shift of eigenvalues',shift,fmt='(es10.3)')

  if (iproc==0) call yaml_sequence_open('purification process')

      ! shift the eigenvalues of the density kernel, using ks as temporary variable
      if (shift/=0.d0) then
          if (ishift==1) stop 'eigenvalue shift not allowed for first iteration'
          do iorb=1,tmb%linmat%l%nfvctr
              do jorb=1,tmb%linmat%l%nfvctr
                  if (jorb==iorb) then
                      ks(jorb,iorb)=kernel_prime(jorb,iorb)+shift
                  else
                      ks(jorb,iorb)=kernel_prime(jorb,iorb)
                  end if
              end do
          end do
          !SM: need to fix the spin here
          call f_zero(tmb%linmat%l%nfvctr**2, tmb%linmat%kernel_%matrix(1,1,1))
          if (tmb%linmat%l%nfvctrp>0) then
              call dgemm('n', 'n', tmb%linmat%l%nfvctr, tmb%linmat%l%nfvctrp, tmb%linmat%l%nfvctr, &
                         1.d0, ks, tmb%linmat%l%nfvctr, &
                         ovrlppowers_(1)%matrix(1,tmb%linmat%l%isfvctr+1,1), tmb%linmat%l%nfvctr, &
                         0.d0, ksksk, tmb%linmat%l%nfvctr) 
              call dgemm('n', 'n', tmb%linmat%l%nfvctr, tmb%linmat%l%nfvctrp, tmb%linmat%l%nfvctr, &
                         1.d0, ovrlppowers_(1)%matrix, tmb%linmat%l%nfvctr, &
                         ksksk, tmb%linmat%l%nfvctr, &
                         0.d0, tmb%linmat%kernel_%matrix(1,tmb%linmat%l%isfvctr+1,1), tmb%linmat%l%nfvctr) 
          end if
    

          if (nproc > 1) then
             !SM: need to fix the spin here
             call mpiallred(tmb%linmat%kernel_%matrix(1,1,1), tmb%linmat%l%nfvctr**2, mpi_sum, comm=bigdft_mpi%mpi_comm)
          end if
      end if


      do it=1,it_opt

          call f_zero(tmb%linmat%l%nfvctr**2, ks(1,1))
          if (tmb%linmat%l%nfvctrp>0) then
              call dgemm('n', 'n', tmb%linmat%l%nfvctr, tmb%linmat%l%nfvctrp, tmb%linmat%l%nfvctr, &
                         1.d0, tmb%linmat%kernel_%matrix(1,1,1), tmb%linmat%l%nfvctr, &
                         tmb%linmat%ovrlp_%matrix(1,tmb%linmat%l%isfvctr+1,1), tmb%linmat%l%nfvctr, &
                         0.d0, ks(1,tmb%linmat%l%isfvctr+1), tmb%linmat%l%nfvctr) 
          end if

          if (nproc > 1) then
              call mpiallred(ks, mpi_sum, comm=bigdft_mpi%mpi_comm)
          end if

          if (tmb%linmat%l%nfvctrp>0) then
              call dgemm('n', 'n', tmb%linmat%l%nfvctr, tmb%linmat%l%nfvctrp, tmb%linmat%l%nfvctr, &
                         1.d0, ks(1,1), tmb%linmat%l%nfvctr, &
                         tmb%linmat%kernel_%matrix(1,tmb%linmat%l%isfvctr+1,1), tmb%linmat%l%nfvctr, &
                         0.d0, ksk(1,1), tmb%linmat%l%nfvctr)
          end if
          if (tmb%linmat%l%nfvctrp>0) then
              call dgemm('n', 'n', tmb%linmat%l%nfvctr, tmb%linmat%l%nfvctrp, tmb%linmat%l%nfvctr, &
                         1.d0, ks(1,1), tmb%linmat%l%nfvctr, &
                         ksk(1,1), tmb%linmat%l%nfvctr, 0.d0, ksksk(1,1), tmb%linmat%l%nfvctr)
          end if


          diff=0.d0
          do iorb=tmb%linmat%l%isfvctr+1,tmb%linmat%l%isfvctr+tmb%linmat%l%nfvctrp
              iiorb=iorb-tmb%linmat%l%isfvctr
              jsegstart=tmb%linmat%l%istsegline(iorb)
              if (iorb<tmb%linmat%l%nfvctr) then
                  jsegend=tmb%linmat%l%istsegline(iorb+1)-1
              else
                  jsegend=tmb%linmat%l%nseg
              end if
              do jseg=jsegstart,jsegend
                  ! A segment is always on one line, therefore no double loop
                  do jorb=tmb%linmat%l%keyg(1,1,jseg),tmb%linmat%l%keyg(2,1,jseg)
                      jjorb=jorb
                      diff = diff + (ksk(jjorb,iiorb)-tmb%linmat%kernel_%matrix(jjorb,iorb,1))**2
                  end do
              end do
          end do

          call compress_matrix2(iproc,tmb%linmat%l, &
               inmat=tmb%linmat%kernel_%matrix, outmat=tmb%linmat%kernel_%matrix_compr)
          !!tmb%linmat%ovrlp_%matrix_compr = tmb%linmat%ovrlp%matrix_compr
          tr_KS=trace_sparse(iproc, nproc, tmb%linmat%s, tmb%linmat%l, &
                tmb%linmat%ovrlp_%matrix_compr(isshift+1:), &
                tmb%linmat%kernel_%matrix_compr(ilshift+1:), ispin)
          if (tmb%linmat%l%nspin==2) then
              chargediff=tr_KS-foe_data_get_real(tmb%foe_obj,"charge",ispin)
          else if (tmb%linmat%l%nspin==1) then
              chargediff=2.d0*tr_KS-foe_data_get_real(tmb%foe_obj,"charge",ispin)
          end if

          if (nproc > 1) then
              call mpiallred(diff, 1, mpi_sum, comm=bigdft_mpi%mpi_comm)
          end if

          diff=sqrt(diff)
          if (iproc==0) then
              call yaml_newline()
              call yaml_sequence(advance='no')
              call yaml_mapping_open(flow=.true.)
              call yaml_map('iter',it)
              call yaml_map('diff from idempotency',diff,fmt='(es9.3)')
              call yaml_map('charge diff',chargediff,fmt='(es10.3)')
              !call yaml_map('alpha',alpha,fmt='(es8.2)')
              call yaml_mapping_close()
          end if

          call f_zero(tmb%linmat%l%nfvctr**2, tmb%linmat%kernel_%matrix(1,1,1))
          do iorb=1,tmb%linmat%l%nfvctrp
              iiorb=iorb+tmb%linmat%l%isfvctr
              do jorb=1,tmb%linmat%l%nfvctr
                  tmb%linmat%kernel_%matrix(jorb,iiorb,1) = 3.d0*ksk(jorb,iorb) - 2.d0*ksksk(jorb,iorb)
              end do
          end do

          if (nproc > 1) then
              call mpiallred(tmb%linmat%kernel_%matrix(1,1,1), &
                   tmb%linmat%l%nfvctr**2, mpi_sum, comm=bigdft_mpi%mpi_comm)
          end if

          if (diff<1.d-10) exit

      end do

      call compress_matrix2(iproc,tmb%linmat%l, &
           inmat=tmb%linmat%kernel_%matrix, outmat=tmb%linmat%kernel_%matrix_compr)
      !!tmb%linmat%ovrlp_%matrix_compr = tmb%linmat%ovrlp%matrix_compr
      tr_KS=trace_sparse(iproc, nproc, tmb%linmat%s, tmb%linmat%l, &
            tmb%linmat%ovrlp_%matrix_compr(isshift+1:), &
            tmb%linmat%kernel_%matrix_compr(ilshift+1:), ispin)
      if (tmb%linmat%l%nspin==2) then
          chargediff=tr_KS-foe_data_get_real(tmb%foe_obj,"charge",ispin)
      else if (tmb%linmat%l%nspin==1) then
          chargediff=2.d0*tr_KS-foe_data_get_real(tmb%foe_obj,"charge",ispin)
      end if

      if (iproc==0) call yaml_sequence_close

      if (abs(chargediff)<1.d-6) exit shift_loop

      if (chargediff>0) then
          ! make this the new upper bound for the bisection
          bisec_bounds(2)=shift
          ! choose new shift, based on whether the lower bound is known or not
          if (bisec_bounds_ok(1)) then
              shift=0.5d0*(bisec_bounds(1)+bisec_bounds(2))
          else
              shift=bisec_bounds(2)-0.01d0
          end if
          bisec_bounds_ok(2)=.true.
      end if
      if (chargediff<0) then
          ! make this the new lower bound for the bisection
          bisec_bounds(1)=shift
          ! choose new shift, based on whether the upper bound is known or not
          if (bisec_bounds_ok(2)) then
              shift=0.5d0*(bisec_bounds(1)+bisec_bounds(2))
          else
              shift=bisec_bounds(1)+0.01d0
          end if
          bisec_bounds_ok(1)=.true.
      end if



  end do shift_loop


  !if (iproc==0) call yaml_sequence_close

  if (tmb%linmat%l%nspin==1) then
      call dscal(tmb%linmat%l%nfvctr**2, 2.0d0, tmb%linmat%kernel_%matrix, 1)
  end if

  call timing(iproc,'purify_kernel ','OF') 

  call f_free(ks)
  call f_free(ksk)
  call f_free(ksksk)
  call f_free(kernel_prime)

  !call f_free_ptr(ovrlp_onehalf)
  !call f_free_ptr(ovrlp_minusonehalf)

  !if (.not.inv_ovrlp_associated) then
  !    call f_free_ptr(tmb%linmat%inv_ovrlp_large%matrix_compr)
  !end if


  call compress_matrix2(iproc, tmb%linmat%l, inmat=tmb%linmat%kernel_%matrix, outmat=tmb%linmat%kernel_%matrix_compr)

  !!tmb%linmat%ovrlp_%matrix_compr = tmb%linmat%ovrlp%matrix_compr
  tr_KS=trace_sparse(iproc, nproc, tmb%linmat%s, tmb%linmat%l, &
        tmb%linmat%ovrlp_%matrix_compr(isshift+1:), &
        tmb%linmat%kernel_%matrix_compr(ilshift+1:), ispin)
  if (iproc==0) then
      call yaml_newline()
      call yaml_map('tr(KS) after purification',tr_KS)
  end if


  call f_free_ptr(tmb%linmat%ovrlp_%matrix)
  call f_free_ptr(tmb%linmat%kernel_%matrix)

  call deallocate_matrices(ovrlp_onehalf_(1))
  call deallocate_matrices(ovrlppowers_(1))

  call f_release_routine()



      contains

        subroutine calculate_overlap_onehalf()
          ! Taylor approximation of S^1/2 and S^-1/2 up to higher order

          call overlapPowerGeneral(iproc, nproc, order_taylor, 1, (/2/), -1, &
               imode=2, ovrlp_smat=tmb%linmat%s, inv_ovrlp_smat=tmb%linmat%l, &
               ovrlp_mat=tmb%linmat%ovrlp_, inv_ovrlp_mat=ovrlp_onehalf_, check_accur=.true., &
               max_error=max_error, mean_error=mean_error)
          call check_taylor_order(mean_error, max_inversion_error, order_taylor)
          call overlapPowerGeneral(iproc, nproc, order_taylor, 1, (/-2/), -1, &
               imode=2, ovrlp_smat=tmb%linmat%s, inv_ovrlp_smat=tmb%linmat%l, &
               ovrlp_mat=tmb%linmat%ovrlp_, inv_ovrlp_mat=ovrlppowers_, check_accur=.true., &
               max_error=max_error, mean_error=mean_error)
          call check_taylor_order(mean_error, max_inversion_error, order_taylor)
          !if (iproc==0) then
          !    call yaml_map('max error of S^-1/2',max_error,fmt='(es9.2)')
          !    call yaml_map('mean error of S^-1/2',mean_error,fmt='(es9.2)')
          !end if
      end subroutine calculate_overlap_onehalf

end subroutine purify_kernel
