!> @file
!!   Routines to calculate the exact exchange potential
!! @author
!!    Copyright (C) 2010-2013 BigDFT group 
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 


!> Calculate the exact exchange potential
subroutine exact_exchange_potential(iproc,nproc,geocode,nspin,lr,orbs,n3parr,n3p,&
     hxh,hyh,hzh,pkernel,psi,psir,eexctX)

  use module_base
  use module_types
  use Poisson_Solver, except_dp => dp, except_gp => gp, except_wp => wp
  use module_xc
  use yaml_output

  implicit none
  character(len=1), intent(in) :: geocode !< @copydoc poisson_solver::coulomb_operator::geocode             !< Determine Boundary conditions
  integer, intent(in) :: iproc,nproc                  !< MPI information
  integer, intent(in) :: n3p,nspin                    !< spin and ...
  real(gp), intent(in) :: hxh,hyh,hzh                 !< hgrid
  type(locreg_descriptors), intent(in) :: lr          !< Local region descriptors
  type(orbitals_data), intent(in) :: orbs             !< Orbitals
  integer, dimension(0:nproc-1), intent(in) :: n3parr
  real(wp), dimension(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f,orbs%nspinor,orbs%norbp), intent(in) :: psi
  type(coulomb_operator), intent(in) :: pkernel
  real(gp), intent(out) :: eexctX
  real(wp), dimension(max(lr%d%n1i*lr%d%n2i*lr%d%n3i*orbs%norbp,n3parr(0)*orbs%norb)), intent(out) :: psir
  !local variables
  character(len=*), parameter :: subname='exact_exchange_potential'
  integer :: i_all,i_stat,ierr,ispinor,ispsiw,ispin,norb,ncall
  integer :: i1,i2,i3p,iorb,iorbs,jorb,jorbs,ispsir,ind3,ind2,ind1i,ind1j,jproc,igran,ngran
  real(gp) :: ehart,hfac,exctXfac,sign,sfac,hfaci,hfacj
  type(workarr_sumrho) :: w
  integer, dimension(:,:), allocatable :: ncommarr
  real(wp), dimension(:), allocatable :: psiw
  real(wp), dimension(:,:,:,:), allocatable :: rp_ij

  !call timing(iproc,'Exchangecorr  ','ON')

  exctXfac = xc_exctXfac()

  eexctX=0.0_gp

  call initialize_work_arrays_sumrho(lr,w)
  
  !save the value of the previous offset of the kernel
  !kerneloff=pkernel(1)
  !put to szero the offset to subtract the energy of the momentum
  !pkernel(1)=0.0_dp

  !the granularity of the calculation is set by ngran
  !for the moment it is irrelevant but if the poisson solver is modified
  !we may increase this value
  ngran=1

  !partial densities with a given granularity
  allocate(rp_ij(lr%d%n1i,lr%d%n2i,n3p,ngran+ndebug),stat=i_stat)
  call memocc(i_stat,rp_ij,'rp_ij',subname)
  allocate(psiw(max(max(lr%d%n1i*lr%d%n2i*lr%d%n3i*orbs%norbp,n3parr(0)*orbs%norb),1)+ndebug),stat=i_stat)
  call memocc(i_stat,psiw,'psiw',subname)

  if (geocode == 'F') then
     call razero(lr%d%n1i*lr%d%n2i*lr%d%n3i*orbs%norbp,psiw)
  end if

  !uncompress the wavefunction in the real grid
  !and switch the values of the function
  ispinor=1
  ispsiw=1
  do iorb=1,orbs%norbp
     call daub_to_isf(lr,w,psi(1,ispinor,iorb),psiw(ispsiw))
     ispsir=1+(iorb-1)*n3parr(0)
     do jproc=0,nproc-1
        !write(*,'(a,1x,8(i10))'),'iproc,jproc',iproc,jproc,iorb,orbs%norbp,ispsir,ispsiw,&
        !     lr%d%n1i*lr%d%n2i*max(lr%d%n3i*orbs%norbp,n3p*orbs%norb),n3parr(jproc)
        call dcopy(n3parr(jproc),psiw(ispsiw),1,psir(ispsir),1)
        ispsiw=ispsiw+n3parr(jproc)
        if (jproc /= nproc-1) then
           do jorb=iorb,orbs%norbp
              ispsir=ispsir+n3parr(jproc)
           end do
           do jorb=1,iorb-1
              ispsir=ispsir+n3parr(jproc+1)
           end do
        end if
     end do
  end do
  call deallocate_work_arrays_sumrho(w)

  !communicate them between processors
  if (nproc > 1) then
     !arrays for the communication between processors
     !valid only for one k-point for the moment
     !and only real functions (nspinor=1)
     !this distribution is in principle valid also for k-points

     allocate(ncommarr(0:nproc-1,4+ndebug),stat=i_stat)
     call memocc(i_stat,ncommarr,'ncommarr',subname)

     !count array for orbitals => components
     do jproc=0,nproc-1
        ncommarr(jproc,1)=n3parr(jproc)*orbs%norb_par(iproc,0)
     end do
     !displacement array for orbitals => components
     ncommarr(0,2)=0
     do jproc=1,nproc-1
        ncommarr(jproc,2)=ncommarr(jproc-1,2)+ncommarr(jproc-1,1)
     end do
     !count array for components => orbitals
     do jproc=0,nproc-1
        ncommarr(jproc,3)=n3parr(iproc)*orbs%norb_par(jproc,0)
     end do
     !displacement array for components => orbitals
     ncommarr(0,4)=0
     do jproc=1,nproc-1
        ncommarr(jproc,4)=ncommarr(jproc-1,4)+ncommarr(jproc-1,3)
     end do

     call MPI_ALLTOALLV(psir,ncommarr(0,1),ncommarr(0,2),mpidtypw, &
          psiw,ncommarr(0,3),ncommarr(0,4),mpidtypw,bigdft_mpi%mpi_comm,ierr)

  else
     call dcopy(lr%d%n1i*lr%d%n2i*n3p*orbs%norb,psir,1,psiw,1)
  end if

  !this is the array of the actions of the X potential on psi
  call razero(lr%d%n1i*lr%d%n2i*n3p*orbs%norb,psir)

  !build the partial densities for the poisson solver, calculate the partial potential
  !and accumulate the result
  !do it for different spins
  !for non spin-polarised systems there is a factor of two
  !non-collinear spin not yet implemented
  if (nspin==2) then
     sfac=1.0_gp
  else 
     sfac=0.5_gp
  end if

  !number of orbitals, all quantum numbers
  norb=orbs%norb!*orbs%nkpts, for the future
  iorb=1
  jorb=1
  ncall=0
  do ispin=1,nspin
     if (ispin==1) then
        iorb=1
        jorb=1
        norb=orbs%norbu
        sign=1.0_gp
     else
        iorb=orbs%norbu+1
        jorb=orbs%norbu+1
        norb=orbs%norb
        sign=-1.0_gp
     end if
     orbital_loop: do
        iorbs=iorb
        jorbs=jorb
        hfac=1/(hxh*hyh*hzh)
        do igran=1,ngran
           if (iorb > norb) exit orbital_loop
           if (orbs%spinsgn(iorb) == sign .and. orbs%spinsgn(jorb) == sign) then
              !calculate partial density (real functions), no spin-polarisation
              do i3p=1,n3p
                 ind3=(i3p-1)*lr%d%n1i*lr%d%n2i
                 do i2=1,lr%d%n2i
                    ind2=(i2-1)*lr%d%n1i+ind3
                    do i1=1,lr%d%n1i
                       ind1i=i1+ind2+(iorb-1)*lr%d%n1i*lr%d%n2i*n3p
                       ind1j=i1+ind2+(jorb-1)*lr%d%n1i*lr%d%n2i*n3p
                       rp_ij(i1,i2,i3p,igran)=hfac*psiw(ind1i)*psiw(ind1j)
                    end do
                 end do
              end do
           end if
           jorb=jorb+1
           if (jorb > norb) then
              iorb=iorb+1
              jorb=iorb
           end if
        end do
        jorb=jorbs
        iorb=iorbs
        do igran=1,ngran
           if (iorb > norb) exit orbital_loop
           if (orbs%spinsgn(iorb) == sign .and. orbs%spinsgn(jorb) == sign) then
              !this factor is only valid with one k-point
              !can be easily generalised to the k-point case
              hfac=sfac*orbs%occup(iorb)*orbs%occup(jorb)

              !print *,'test',iproc,iorb,jorb,sum(rp_ij(:,:,:,igran))
              !partial exchange term for each partial density
              ncall=ncall+1
              if (iproc == 0 .and. verbose > 1) then
                 !write(*,*)'Exact exchange calculation: spin, orbitals:',ispin,iorb,jorb
                 call yaml_comment('Exact exchange calculation: ' // trim(yaml_toa( nint(real(ncall,gp)/ &
                 &    real(orbs%norbu*(orbs%norbu+1)/2+orbs%norbd*(orbs%norbd+1)/2,gp)*100.0_gp),fmt='(i3)')) // '%')
                 !write(*,'(1x,a,i3,a2)')'Exact exchange calculation:',&
                 !     nint(real(ncall,gp)/real(orbs%norbu*(orbs%norbu+1)/2+orbs%norbd*(orbs%norbd+1)/2,gp)*100.0_gp),'% '
              end if
              
              call H_potential('D',pkernel,rp_ij(1,1,1,igran),rp_ij,ehart,0.0_dp,.false.,&
                   quiet='YES')

!!$              call PSolver(geocode,'D',iproc,nproc,lr%d%n1i,lr%d%n2i,lr%d%n3i,&
!!$                   0,hxh,hyh,hzh,rp_ij(1,1,1,igran),pkernel,rp_ij,ehart,zero,zero,&
!!$                   0.d0,.false.,1,quiet='YES')
              if (iorb==jorb) then
                 eexctX=eexctX+hfac*real(ehart,gp)
              else
                 eexctX=eexctX+2.0_gp*hfac*real(ehart,gp)
              end if
              !print *,'PSOLVER,ehart,iproc',iproc,ehart,hfac
           end if
           jorb=jorb+1
           if (jorb > norb) then
              iorb=iorb+1
              jorb=iorb
           end if
        end do
        jorb=jorbs
        iorb=iorbs
        do igran=1,ngran
           if (iorb > norb) exit orbital_loop
           if (orbs%spinsgn(iorb) == sign .and. orbs%spinsgn(jorb) == sign) then
              !this factor is only valid with one k-point
              !we have to correct with the kwgts if we want more than one k-point
              hfaci=-sfac*orbs%occup(jorb)
              hfacj=-sfac*orbs%occup(iorb)

              if (iorb /= jorb) then
                 !accumulate the results for each of the wavefunctions concerned
                 do i3p=1,n3p
                    ind3=(i3p-1)*lr%d%n1i*lr%d%n2i
                    do i2=1,lr%d%n2i
                       ind2=(i2-1)*lr%d%n1i+ind3
                       do i1=1,lr%d%n1i
                          ind1i=i1+ind2+(iorb-1)*lr%d%n1i*lr%d%n2i*n3p
                          ind1j=i1+ind2+(jorb-1)*lr%d%n1i*lr%d%n2i*n3p
                          psir(ind1i)=psir(ind1i)+hfaci*rp_ij(i1,i2,i3p,igran)*psiw(ind1j)
                          psir(ind1j)=psir(ind1j)+hfacj*rp_ij(i1,i2,i3p,igran)*psiw(ind1i)
                       end do
                    end do
                 end do
              else
                 !accumulate the results for each of the wavefunctions concerned
                 do i3p=1,n3p
                    ind3=(i3p-1)*lr%d%n1i*lr%d%n2i
                    do i2=1,lr%d%n2i
                       ind2=(i2-1)*lr%d%n1i+ind3
                       do i1=1,lr%d%n1i
                          ind1i=i1+ind2+(iorb-1)*lr%d%n1i*lr%d%n2i*n3p
                          ind1j=i1+ind2+(jorb-1)*lr%d%n1i*lr%d%n2i*n3p
                          psir(ind1i)=psir(ind1i)+hfacj*rp_ij(i1,i2,i3p,igran)*psiw(ind1j)
                       end do
                    end do
                 end do
              end if
           end if
           jorb=jorb+1
           if (jorb > norb) then
              iorb=iorb+1
              jorb=iorb
           end if
        end do
     end do orbital_loop
  end do

  !the exact exchange energy is half the Hartree energy (which already has another half)
  eexctX=-exctXfac*eexctX

  if (iproc == 0) call yaml_map('Exact Exchange Energy',eexctX,fmt='(1pe18.1)')
  !if (iproc == 0) write(*,'(1x,a,1x,1pe18.11)')'Exact Exchange Energy:',eexctX

  !assign the potential for each function
  if (nproc > 1) then
     !call dcopy(lr%d%n1i*lr%d%n2i*n3p*orbs%norb,psir,1,psirt,1)
     !recommunicate the values in the psir array
     call MPI_ALLTOALLV(psir,ncommarr(0,3),ncommarr(0,4),mpidtypw, &
          psiw,ncommarr(0,1),ncommarr(0,2),mpidtypw,bigdft_mpi%mpi_comm,ierr)
     !redress the potential
     ispsiw=1
     do iorb=1,orbs%norbp
        ispsir=1+(iorb-1)*n3parr(0)
        do jproc=0,nproc-1
           call dcopy(n3parr(jproc),psiw(ispsir),1,psir(ispsiw),1)
           ispsiw=ispsiw+n3parr(jproc)
           if (jproc /= nproc-1) then
              do jorb=iorb,orbs%norbp
                 ispsir=ispsir+n3parr(jproc)
              end do
              do jorb=1,iorb-1
                 ispsir=ispsir+n3parr(jproc+1)
              end do
           end if
        end do
     end do
  end if

  i_all=-product(shape(rp_ij))*kind(rp_ij)
  deallocate(rp_ij,stat=i_stat)
  call memocc(i_stat,i_all,'rp_ij',subname)
  
  i_all=-product(shape(psiw))*kind(psiw)
  deallocate(psiw,stat=i_stat)
  call memocc(i_stat,i_all,'psiw',subname)


  if (nproc > 1) then
     i_all=-product(shape(ncommarr))*kind(ncommarr)
     deallocate(ncommarr,stat=i_stat)
     call memocc(i_stat,i_all,'ncommarr',subname)
  end if

  !call timing(iproc,'Exchangecorr  ','OF')

END SUBROUTINE exact_exchange_potential


subroutine prepare_psirocc(iproc,nproc,lr,orbsocc,n3p,n3parr,psiocc,psirocc)
  use module_base
  use module_types
  implicit none
  integer, intent(in) :: iproc,nproc,n3p
  type(locreg_descriptors), intent(in) :: lr
  type(orbitals_data), intent(in) :: orbsocc
  integer, dimension(0:nproc-1), intent(in) :: n3parr
  real(wp), dimension(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f,orbsocc%nspinor,orbsocc%norbp), intent(in) :: psiocc
  real(wp), dimension(max(lr%d%n1i*lr%d%n2i*lr%d%n3i*orbsocc%norbp,n3parr(0)*orbsocc%norb)), intent(out) :: psirocc
  !local variables
  character(len=*), parameter :: subname='prepare_psirocc'
  integer :: i_all,i_stat,ierr,ispinor,ispsiw
  integer :: iorb,jorb,ispsir,jproc
  type(workarr_sumrho) :: w
  integer, dimension(:,:), allocatable :: ncommocc
  real(wp), dimension(:), allocatable :: psiwocc

  call initialize_work_arrays_sumrho(lr,w)

  allocate(psiwocc(max(max(lr%d%n1i*lr%d%n2i*lr%d%n3i*orbsocc%norbp,n3parr(0)*orbsocc%norb),1)+ndebug),stat=i_stat)
  call memocc(i_stat,psiwocc,'psiwocc',subname)

  call razero(max(max(lr%d%n1i*lr%d%n2i*lr%d%n3i*orbsocc%norbp,n3parr(0)*orbsocc%norb),1),psiwocc)

  if (lr%geocode == 'F') then
     call razero(lr%d%n1i*lr%d%n2i*lr%d%n3i*orbsocc%norbp,psirocc)
  end if

  !uncompress the wavefunction in the real grid
  !and switch the values of the function

  !occupied orbitals
  ispinor=1
  ispsiw=1
  do iorb=1,orbsocc%norbp
     call daub_to_isf(lr,w,psiocc(1,ispinor,iorb),psirocc(ispsiw))
     ispsir=1+(iorb-1)*n3parr(0)
     do jproc=0,nproc-1
        !write(*,'(a,1x,8(i10))'),'iproc,jproc',iproc,jproc,iorb,orbs%norbp,ispsir,ispsiw,&
        !     lr%d%n1i*lr%d%n2i*max(lr%d%n3i*orbs%norbp,n3p*orbs%norb),n3parr(jproc)
        call dcopy(n3parr(jproc),psirocc(ispsiw),1,psiwocc(ispsir),1)
        ispsiw=ispsiw+n3parr(jproc)
        if (jproc /= nproc-1) then
           do jorb=iorb,orbsocc%norbp
              ispsir=ispsir+n3parr(jproc)
           end do
           do jorb=1,iorb-1
              ispsir=ispsir+n3parr(jproc+1)
           end do
        end if
     end do
  end do
  call deallocate_work_arrays_sumrho(w)

  !communicate them between processors
  !occupied orbitals
  if (nproc > 1) then
     !arrays for the communication between processors
     !valid only for one k-point for the moment
     !and only real functions (nspinor=1)
     !this distribution is in principle valid also for k-points
     allocate(ncommocc(0:nproc-1,4+ndebug),stat=i_stat)
     call memocc(i_stat,ncommocc,'ncommocc',subname)

     !count occay for orbitals => components
     do jproc=0,nproc-1
        ncommocc(jproc,1)=n3parr(jproc)*orbsocc%norb_par(iproc,0)
     end do
     !displacement array for orbitals => components
     ncommocc(0,2)=0
     do jproc=1,nproc-1
        ncommocc(jproc,2)=ncommocc(jproc-1,2)+ncommocc(jproc-1,1)
     end do
     !count occay for components => orbitals
     do jproc=0,nproc-1
        ncommocc(jproc,3)=n3parr(iproc)*orbsocc%norb_par(jproc,0)
     end do
     !displacement array for components => orbitals
     ncommocc(0,4)=0
     do jproc=1,nproc-1
        ncommocc(jproc,4)=ncommocc(jproc-1,4)+ncommocc(jproc-1,3)
     end do

     call MPI_ALLTOALLV(psiwocc,ncommocc(0,1),ncommocc(0,2),mpidtypw, &
          psirocc,ncommocc(0,3),ncommocc(0,4),mpidtypw,bigdft_mpi%mpi_comm,ierr)
  else
     call dcopy(lr%d%n1i*lr%d%n2i*n3p*orbsocc%norb,psiwocc,1,psirocc,1)
  end if
  i_all=-product(shape(psiwocc))*kind(psiwocc)
  deallocate(psiwocc,stat=i_stat)
  call memocc(i_stat,i_all,'psiwocc',subname)

  if (nproc > 1) then
     i_all=-product(shape(ncommocc))*kind(ncommocc)
     deallocate(ncommocc,stat=i_stat)
     call memocc(i_stat,i_all,'ncommocc',subname)
  end if

END SUBROUTINE prepare_psirocc


!> Calculate the exact exchange potential only on virtual orbitals
!! by knowing the occupied orbitals and their distribution
!! both sets of orbitals are to be 
subroutine exact_exchange_potential_virt(iproc,nproc,geocode,nspin,lr,orbsocc,orbsvirt,n3parr,n3p,&
     hxh,hyh,hzh,pkernel,psirocc,psivirt,psirvirt)
  use module_base
  use module_types
  use Poisson_Solver, except_dp => dp, except_gp => gp, except_wp => wp
  use yaml_output
  implicit none
  character(len=1), intent(in) :: geocode !< @copydoc poisson_solver::coulomb_operator::geocode
  integer, intent(in) :: iproc,nproc,n3p,nspin
  real(gp), intent(in) :: hxh,hyh,hzh
  type(locreg_descriptors), intent(in) :: lr
  type(orbitals_data), intent(in) :: orbsocc,orbsvirt
  integer, dimension(0:nproc-1), intent(in) :: n3parr
  real(wp), dimension(max(lr%d%n1i*lr%d%n2i*lr%d%n3i*orbsocc%norbp,n3parr(0)*orbsocc%norb)), intent(in) :: psirocc
  real(wp), dimension(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f,orbsvirt%nspinor,orbsvirt%norbp), intent(in) :: psivirt
  type(coulomb_operator), intent(in) :: pkernel
  real(wp), dimension(max(lr%d%n1i*lr%d%n2i*lr%d%n3i*orbsvirt%norbp,n3parr(0)*orbsvirt%norb)), intent(out) :: psirvirt
  !local variables
  character(len=*), parameter :: subname='exact_exchange_potential_virt'
  integer :: i_all,i_stat,ierr,ispinor,ispsiw,ispin,norbocc,norbvirt
  integer :: i1,i2,i3p,iorb,iorbs,jorb,jorbs,ispsir,ind3,ind2,ind1i,ind1j,jproc,igran,ngran
  real(gp) :: ehart,hfac,sign,sfac,hfaci
  type(workarr_sumrho) :: w
  integer, dimension(:,:), allocatable :: ncommvirt
  real(wp), dimension(:), allocatable :: psiwvirt
  real(wp), dimension(:,:,:,:), allocatable :: rp_ij

  !call timing(iproc,'Exchangecorr  ','ON')

  call initialize_work_arrays_sumrho(lr,w)
  
  !the granularity of the calculation is set by ngran
  !for the moment it is irrelevant but if the poisson solver is modified
  !we may increase this value
  ngran=1

  !partial densities with a given granularity
  allocate(rp_ij(lr%d%n1i,lr%d%n2i,n3p,ngran+ndebug),stat=i_stat)
  call memocc(i_stat,rp_ij,'rp_ij',subname)
  allocate(psiwvirt(max(max(lr%d%n1i*lr%d%n2i*lr%d%n3i*orbsvirt%norbp,n3parr(0)*orbsvirt%norb),1)+ndebug),stat=i_stat)
  call memocc(i_stat,psiwvirt,'psiwvirt',subname)

  if (geocode == 'F') then
     call razero(lr%d%n1i*lr%d%n2i*lr%d%n3i*orbsvirt%norbp,psiwvirt)
  end if

  !uncompress the wavefunction in the real grid
  !and switch the values of the function
  !unoccupied orbitals
  ispinor=1
  ispsiw=1
  do iorb=1,orbsvirt%norbp
     call daub_to_isf(lr,w,psivirt(1,ispinor,iorb),psiwvirt(ispsiw))
     ispsir=1+(iorb-1)*n3parr(0)
     do jproc=0,nproc-1
        !write(*,'(a,1x,8(i10))'),'iproc,jproc',iproc,jproc,iorb,orbs%norbp,ispsir,ispsiw,&
        !     lr%d%n1i*lr%d%n2i*max(lr%d%n3i*orbs%norbp,n3p*orbs%norb),n3parr(jproc)
        call dcopy(n3parr(jproc),psiwvirt(ispsiw),1,psirvirt(ispsir),1)
        ispsiw=ispsiw+n3parr(jproc)
        if (jproc /= nproc-1) then
           do jorb=iorb,orbsvirt%norbp
              ispsir=ispsir+n3parr(jproc)
           end do
           do jorb=1,iorb-1
              ispsir=ispsir+n3parr(jproc+1)
           end do
        end if
     end do
  end do

  call deallocate_work_arrays_sumrho(w)

  !communicate them between processors
  !unoccupied orbitals
  if (nproc > 1) then
     !arrays for the communication between processors
     !valid only for one k-point for the moment
     !and only real functions (nspinor=1)
     !this distribution is in principle valid also for k-points

     allocate(ncommvirt(0:nproc-1,4+ndebug),stat=i_stat)
     call memocc(i_stat,ncommvirt,'ncommvirt',subname)

     !count array for orbitals => components
     do jproc=0,nproc-1
        ncommvirt(jproc,1)=n3parr(jproc)*orbsvirt%norb_par(iproc,0)
     end do
     !displacement array for orbitals => components
     ncommvirt(0,2)=0
     do jproc=1,nproc-1
        ncommvirt(jproc,2)=ncommvirt(jproc-1,2)+ncommvirt(jproc-1,1)
     end do
     !count array for components => orbitals
     do jproc=0,nproc-1
        ncommvirt(jproc,3)=n3parr(iproc)*orbsvirt%norb_par(jproc,0)
     end do
     !displacement array for components => orbitals
     ncommvirt(0,4)=0
     do jproc=1,nproc-1
        ncommvirt(jproc,4)=ncommvirt(jproc-1,4)+ncommvirt(jproc-1,3)
     end do

     call MPI_ALLTOALLV(psirvirt,ncommvirt(0,1),ncommvirt(0,2),mpidtypw, &
          psiwvirt,ncommvirt(0,3),ncommvirt(0,4),mpidtypw,bigdft_mpi%mpi_comm,ierr)

  else
     call dcopy(lr%d%n1i*lr%d%n2i*n3p*orbsvirt%norb,psirvirt,1,psiwvirt,1)
  end if

  !this is the array of the actions of the X potential on psi
  call razero(lr%d%n1i*lr%d%n2i*n3p*orbsvirt%norb,psirvirt)

  !build the partial densities for the poisson solver, calculate the partial potential
  !and accumulate the result
  !do it for different spins
  !for non spin-polarised systems there is a factor of two
  !non-collinear spin not yet implemented
  if (nspin==2) then
     sfac=1.0_gp
  else 
     sfac=0.5_gp
  end if

  !number of orbitals, all quantum numbers
  norbocc=orbsocc%norb!*orbs%nkpts, for the future
  norbvirt=orbsvirt%norb
  iorb=1 !index on virtual orbitals
  jorb=1 !index on occupied orbitals
  
  do ispin=1,nspin
     if (ispin==1) then
        iorb=1
        jorb=1
        norbocc=orbsocc%norbu
        norbvirt=orbsvirt%norbu
        sign=1.0_gp
     else
        iorb=orbsvirt%norbu+1
        jorb=orbsocc%norbu+1
        norbocc=orbsocc%norb
        norbvirt=orbsvirt%norb
        sign=-1.0_gp
     end if
     orbital_loop: do
        iorbs=iorb
        jorbs=jorb
        hfac=1/(hxh*hyh*hzh)
        do igran=1,ngran
           if (iorb > norbvirt) exit orbital_loop
           if (orbsvirt%spinsgn(iorb) == sign .and. orbsocc%spinsgn(jorb) == sign) then
              !calculate partial density (real functions), no spin-polarisation
              do i3p=1,n3p
                 ind3=(i3p-1)*lr%d%n1i*lr%d%n2i
                 do i2=1,lr%d%n2i
                    ind2=(i2-1)*lr%d%n1i+ind3
                    do i1=1,lr%d%n1i
                       ind1i=i1+ind2+(iorb-1)*lr%d%n1i*lr%d%n2i*n3p
                       ind1j=i1+ind2+(jorb-1)*lr%d%n1i*lr%d%n2i*n3p
                       rp_ij(i1,i2,i3p,igran)=hfac*psiwvirt(ind1i)*psirocc(ind1j)
                    end do
                 end do
              end do
           end if
           jorb=jorb+1
           if (jorb > norbocc) then
              iorb=iorb+1
              jorb=1
           end if
        end do
        jorb=jorbs
        iorb=iorbs
        do igran=1,ngran
           if (iorb > norbvirt) exit orbital_loop
           if (orbsvirt%spinsgn(iorb) == sign .and. orbsocc%spinsgn(jorb) == sign) then

              !partial exchange term for each partial density
              if (iproc == 0 .and. verbose > 1) then
                 call yaml_map('Exact exchange calculation (spin, orbitals)', (/ispin,iorb,jorb /))
                 !write(*,*)'Exact exchange calculation: spin, orbitals:',ispin,iorb,jorb
              end if

              call H_potential('D',pkernel,rp_ij(1,1,1,igran),rp_ij,ehart,0.0_dp,.false.,&
                   quiet='YES')

!!$              call PSolver(geocode,'D',iproc,nproc,lr%d%n1i,lr%d%n2i,lr%d%n3i,&
!!$                   0,hxh,hyh,hzh,rp_ij(1,1,1,igran),pkernel,rp_ij,ehart,zero,zero,&
!!$                   0.d0,.false.,1,quiet='YES')

           end if
           jorb=jorb+1
           if (jorb > norbocc) then
              iorb=iorb+1
              jorb=1
           end if
        end do
        jorb=jorbs
        iorb=iorbs
        do igran=1,ngran
           if (iorb > norbvirt) exit orbital_loop
           if (orbsvirt%spinsgn(iorb) == sign .and. orbsocc%spinsgn(jorb) == sign) then
              !this factor is only valid with one k-point
              !we have to correct with the kwgts if we want more than one k-point
              hfaci=-sfac*orbsocc%occup(jorb)

              !accumulate the results for each of the wavefunctions concerned
              do i3p=1,n3p
                 ind3=(i3p-1)*lr%d%n1i*lr%d%n2i
                 do i2=1,lr%d%n2i
                    ind2=(i2-1)*lr%d%n1i+ind3
                    do i1=1,lr%d%n1i
                       ind1i=i1+ind2+(iorb-1)*lr%d%n1i*lr%d%n2i*n3p
                       ind1j=i1+ind2+(jorb-1)*lr%d%n1i*lr%d%n2i*n3p
                       psirvirt(ind1i)=psirvirt(ind1i)+hfaci*rp_ij(i1,i2,i3p,igran)*psirocc(ind1j)
                    end do
                 end do
              end do
           end if
           jorb=jorb+1
           if (jorb > norbocc) then
              iorb=iorb+1
              jorb=1
           end if
        end do
     end do orbital_loop
  end do

  !assign the potential for each function
  if (nproc > 1) then
     !call dcopy(lr%d%n1i*lr%d%n2i*n3p*orbs%norb,psir,1,psirt,1)
     !recommunicate the values in the psir array
     call MPI_ALLTOALLV(psirvirt,ncommvirt(0,3),ncommvirt(0,4),mpidtypw, &
          psiwvirt,ncommvirt(0,1),ncommvirt(0,2),mpidtypw,bigdft_mpi%mpi_comm,ierr)
     !redress the potential
     ispsiw=1
     do iorb=1,orbsvirt%norbp
        ispsir=1+(iorb-1)*n3parr(0)
        do jproc=0,nproc-1
           call dcopy(n3parr(jproc),psiwvirt(ispsir),1,psirvirt(ispsiw),1)
           ispsiw=ispsiw+n3parr(jproc)
           if (jproc /= nproc-1) then
              do jorb=iorb,orbsvirt%norbp
                 ispsir=ispsir+n3parr(jproc)
              end do
              do jorb=1,iorb-1
                 ispsir=ispsir+n3parr(jproc+1)
              end do
           end if
        end do
     end do
  end if

  i_all=-product(shape(rp_ij))*kind(rp_ij)
  deallocate(rp_ij,stat=i_stat)
  call memocc(i_stat,i_all,'rp_ij',subname)

  i_all=-product(shape(psiwvirt))*kind(psiwvirt)
  deallocate(psiwvirt,stat=i_stat)
  call memocc(i_stat,i_all,'psiwvirt',subname)


  if (nproc > 1) then
     i_all=-product(shape(ncommvirt))*kind(ncommvirt)
     deallocate(ncommvirt,stat=i_stat)
     call memocc(i_stat,i_all,'ncommvirt',subname)
  end if

END SUBROUTINE exact_exchange_potential_virt


!> Calculate the exact exchange potential on occupied orbitals
!! within the symmetric round-robin scheme
!! the psi is already given in the real-space form
subroutine exact_exchange_potential_round(iproc,nproc,nspin,lr,orbs,&
     hxh,hyh,hzh,pkernel,psi,dpsir,eexctX)
  use module_base
  use module_types
  use Poisson_Solver, except_dp => dp, except_gp => gp, except_wp => wp
  use module_xc
  use yaml_output
  implicit none
  integer, intent(in) :: iproc,nproc,nspin
  real(gp), intent(in) :: hxh,hyh,hzh
  type(locreg_descriptors), intent(in) :: lr
  type(orbitals_data), intent(in) :: orbs
  real(wp), dimension(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f,orbs%nspinor,orbs%norbp), intent(in) :: psi
  type(coulomb_operator), intent(in) :: pkernel
  real(gp), intent(out) :: eexctX
  real(wp), dimension(lr%d%n1i*lr%d%n2i*lr%d%n3i,orbs%norbp), intent(out) :: dpsir
  !local variables
  character(len=*), parameter :: subname='exact_exchange_potential_round'
  logical :: doit
  integer :: i_all,i_stat,ierr,ncommsstep,ncommsstep2,isnow,irnow,isnow2,irnow2,jsorb,kproc,norbp
  integer :: i,iorb,jorb,jproc,igroup,ngroup,ngroupp,nend,isorb,iorbs,jorbs
  integer :: icount,nprocgr,iprocgrs,iprocgrr,itestproc,norbi,norbj,ncalltot,icountmax,iprocref,ncalls
  real(gp) :: ehart,hfac,exctXfac,sfac,hfaci,hfacj,hfac2
  integer, dimension(4) :: mpireq,mpireq2
  integer, dimension(MPI_STATUS_SIZE,4) :: mpistat,mpistat2
  type(workarr_sumrho) :: w
  integer, dimension(:), allocatable :: igrpr,ndatac
  integer, dimension(:,:), allocatable :: nvctr_par
  integer, dimension(:,:,:), allocatable :: jprocsr,iprocpm1,ndatas,iorbgr
  real(wp), dimension(:), allocatable :: rp_ij
  real(wp), dimension(:,:), allocatable :: psir
  real(wp), dimension(:,:,:,:), allocatable :: psiw,dpsiw

  !call timing(iproc,'Exchangecorr  ','ON')

  exctXfac = xc_exctXfac()

  eexctX=0.0_gp

  !build the partial densities for the poisson solver, calculate the partial potential
  !and accumulate the result
  !do it for different spins
  !for non spin-polarised systems there is a factor of two
  !non-collinear spin not yet implemented
  if (nspin==2) then
     sfac=1.0_gp
     ngroup=2
  else 
     sfac=0.5_gp
     ngroup=1
  end if

  hfac=1.0_gp/(hxh*hyh*hzh)

  !here we can start with the round-robin scheme
  !since the orbitals are all occupied we have to use the symmetric scheme
  !we have first to define the number of groups, which correspond to the repartition 
  !of spin up and spin down orbitals
  allocate(nvctr_par(0:nproc-1,ngroup+ndebug),stat=i_stat)
  call memocc(i_stat,nvctr_par,'nvctr_par',subname)

  allocate(iorbgr(2,0:nproc-1,ngroup+ndebug),stat=i_stat)
  call memocc(i_stat,iorbgr,'iorbgr',subname)

  !test array for data sending
  allocate(ndatas(2,0:nproc-1,ngroup+ndebug),stat=i_stat)
  call memocc(i_stat,ndatas,'ndatas',subname)

  
  if (ngroup==2) then
     isorb=0
     do jproc=0,nproc-1
        iorbgr(1,jproc,1)=isorb
        iorbgr(2,jproc,1)=1
        norbp=max(min(isorb+orbs%norb_par(jproc,0),orbs%norbu)-isorb,0)
        if (norbp == 0) then
           iorbgr(1,jproc,1)=0
           iorbgr(1,jproc,2)=0
        end if
        nvctr_par(jproc,1)=norbp*lr%d%n1i*lr%d%n2i*lr%d%n3i
        iorbgr(1,jproc,2)=isorb
        iorbgr(2,jproc,2)=norbp+1
        norbp=max(isorb+orbs%norb_par(jproc,0)-max(orbs%norbu,isorb),0)
        if (norbp == 0) then
           iorbgr(1,jproc,2)=0
           iorbgr(2,jproc,2)=1
        end if
        nvctr_par(jproc,2)=norbp*lr%d%n1i*lr%d%n2i*lr%d%n3i
        isorb=isorb+orbs%norb_par(jproc,0)
     end do
!!$     if (iproc ==0) then
!!$        print '(a,10(1x,i8))','iproc,nvctr_parA',iproc,nvctr_par(:,1)/(lr%d%n1i*lr%d%n2i*lr%d%n3i)
!!$        print '(a,10(1x,i8))','iproc,nvctr_parB',iproc,nvctr_par(:,2)/(lr%d%n1i*lr%d%n2i*lr%d%n3i)
!!$        print '(a,10(1x,i8))','iproc,iorbgr',iproc,iorbgr
!!$     end if
  else
     isorb=0
     do jproc=0,nproc-1
        iorbgr(1,jproc,1)=isorb
        iorbgr(2,jproc,1)=1
        nvctr_par(jproc,1)=orbs%norb_par(jproc,0)*lr%d%n1i*lr%d%n2i*lr%d%n3i
        isorb=isorb+orbs%norb_par(jproc,0)
     end do
  end if


  !here we can allocate the working arrays giving the maximum
  !between the components for each group
  ngroupp=0
  do igroup=1,ngroup
     if (nvctr_par(iproc,igroup) > 0) then
        ngroupp=ngroupp+1
     end if
  end do

  !determine the array of the groups which are of interest for this processor
  allocate(igrpr(ngroupp+ndebug),stat=i_stat)
  call memocc(i_stat,igrpr,'igrpr',subname)
  allocate(iprocpm1(2,0:nproc-1,ngroupp+ndebug),stat=i_stat)
  call memocc(i_stat,iprocpm1,'iprocpm1',subname)

  !test array for data calculation
  allocate(ndatac(ngroupp+ndebug),stat=i_stat)
  call memocc(i_stat,ndatac,'ndatac',subname)


  !determine for each processor the groups which has to be used
  icount=0
  do igroup=1,ngroup
     if (nvctr_par(iproc,igroup) > 0) then
        icount=icount+1
        igrpr(icount)=igroup
     end if
  end do

  !calculate the processor which lies after and before the present in the list
  iprocpm1=-1
  do igroup=1,ngroupp
     iprocgrs=-1
     iprocgrr=-1
     !define the number of data to calculate in total
     ndatac(igroup)=0
     do kproc=0,nproc-1
        ndatac(igroup)=ndatac(igroup)-nvctr_par(kproc,igrpr(igroup))
        if (nvctr_par(modulo(iproc+kproc,nproc),igrpr(igroup)) > 0) then
           iprocgrs=iprocgrs+1
           iprocpm1(1,iprocgrs,igroup)=modulo(iproc+kproc,nproc)
        end if
        if (nvctr_par(modulo(iproc-kproc,nproc),igrpr(igroup)) > 0) then
           iprocgrr=iprocgrr+1
           iprocpm1(2,iprocgrr,igroup)=modulo(iproc-kproc,nproc)
        end if
     end do
  end do

  !find the processor whih has the maximum number of groups
  icountmax=0
  do kproc=0,nproc-1
     icount=0
     do igroup=1,ngroup
        if (nvctr_par(kproc,igroup) > 0) then
           icount=icount+1
        end if
     end do
     if (icount > icountmax) then
        iprocref=kproc
        icountmax=icount
     end if
  end do

  !calculate the list of send-receive operations which have to be performed per group
  !allocate it at the maximum size needed
  allocate(jprocsr(4,0:nproc/2+1,ngroupp+ndebug),stat=i_stat)
  call memocc(i_stat,jprocsr,'jprocsr',subname)
  !initalise array to minus one
  jprocsr=-1

  ncalltot=0
  do igroup=1,ngroupp
  ncalltot=ncalltot+&
       (nvctr_par(iproc,igrpr(igroup))/(lr%d%n1i*lr%d%n2i*lr%d%n3i))*&
       (nvctr_par(iproc,igrpr(igroup))/(lr%d%n1i*lr%d%n2i*lr%d%n3i)+1)/2
       !calculate the number of processors per group
     nprocgr=0
     do kproc=0,nproc-1
        if (nvctr_par(kproc,igrpr(igroup)) > 0) nprocgr=nprocgr+1
     end do
     
     !do not send anything if there is only one member in the group
     if (nprocgr > 1) then
        do kproc=0,(nprocgr-1)/2-1
           !define the arrays for send-receive of data
           jprocsr(1,kproc,igroup)=iprocpm1(2,kproc,igroup)
           jprocsr(2,kproc,igroup)=iprocpm1(2,kproc+1,igroup)
           if (iproc == iprocref) then
              ncalltot=ncalltot+&
                   (nvctr_par(jprocsr(2,kproc,igroup),igrpr(igroup))/(lr%d%n1i*lr%d%n2i*lr%d%n3i))*&
                   (nvctr_par(iproc,igrpr(igroup))/(lr%d%n1i*lr%d%n2i*lr%d%n3i))
             end if
           if (kproc > 0) then
              jprocsr(3,kproc,igroup)=iprocpm1(2,kproc,igroup)
              jprocsr(4,kproc,igroup)=iprocpm1(1,kproc,igroup)
           end if
        end do
        kproc=(nprocgr-1)/2
        !the last step behaves differently if the group number is odd or even
        if (modulo(nprocgr,2) == 0) then
           jprocsr(1,kproc,igroup)=iprocpm1(2,kproc,igroup)
           jprocsr(2,kproc,igroup)=iprocpm1(2,kproc+1,igroup)
           if (iproc == iprocref) then
              ncalltot=ncalltot+&
                   (nvctr_par(jprocsr(2,kproc,igroup),igrpr(igroup))/(lr%d%n1i*lr%d%n2i*lr%d%n3i))*&
                   (nvctr_par(iproc,igrpr(igroup))/(lr%d%n1i*lr%d%n2i*lr%d%n3i))
             end if
           if (kproc > 0) then
              jprocsr(3,kproc,igroup)=iprocpm1(2,kproc,igroup)
              jprocsr(4,kproc,igroup)=iprocpm1(1,kproc,igroup)
           end if
        else
           jprocsr(3,kproc,igroup)=iprocpm1(2,kproc,igroup)
           jprocsr(4,kproc,igroup)=iprocpm1(1,kproc,igroup)
        end if
     end if
  end do

  itestproc=-1 !-1=no debug verbosity
  if (itestproc > -1) then
     !simulation of communication
     isnow=1
     isnow2=1
     nend=(nproc-1)/2+1
     ncommsstep2=0
     ndatas=0

     do jproc=0,nend
        irnow=3-isnow
        ncommsstep=0
        !sending receiving data
        do igroup=1,ngroupp
           if (jprocsr(1,jproc,igroup) /= -1) then
              ncommsstep=ncommsstep+1
              !send the fixed array to the processor which comes in the list
              if (iprocpm1(1,1,igroup) == itestproc) then
                 print *,'step',jproc+1,': sending',nvctr_par(jprocsr(1,jproc,igroup),igrpr(igroup)),&
                      'elements from',iproc,'to',iprocpm1(1,1,igroup)
              end if
              ndatas(1,iprocpm1(1,1,igroup),igrpr(igroup))=ndatas(1,iprocpm1(1,1,igroup),igrpr(igroup))+&
                   nvctr_par(jprocsr(1,jproc,igroup),igrpr(igroup))
           end if
           if (jprocsr(2,jproc,igroup) /= -1) then
              ncommsstep=ncommsstep+1
              if (iproc == itestproc) then
                 print *,'step',jproc+1,': receiving',nvctr_par(jprocsr(2,jproc,igroup),igrpr(igroup)),&
                      'elements from',iprocpm1(2,1,igroup),'to',iproc
              end if
              ndatas(1,iproc,igrpr(igroup))=ndatas(1,iproc,igrpr(igroup))-&
                   nvctr_par(jprocsr(2,jproc,igroup),igrpr(igroup))
           end if
        end do

        !calculation for orbitals to be performed
        do igroup=1,ngroupp
           if (jproc==0) then
              ndatac(igroup)=ndatac(igroup)+nvctr_par(iproc,igrpr(igroup))
           else
              if (jprocsr(2,jproc-1,igroup) /=-1) then
                 ndatac(igroup)=ndatac(igroup)+&
                      nvctr_par(jprocsr(2,jproc-1,igroup),igrpr(igroup))  
                 if (iproc == itestproc) then
                    print '(5(1x,a,i8))','step',jproc+1,'group:',igrpr(igroup),&
                         ':processing',nvctr_par(jprocsr(2,jproc-1,igroup),igrpr(igroup)),&
                         'elements in',iproc,'from',jprocsr(2,jproc-1,igroup)
                 end if
              end if
           end if
        end do

        !copy the results which have been received
        if (ncommsstep2 > 0) then
           do igroup=1,ngroupp
              if (jprocsr(4,jproc-1,igroup) /= -1) then
                 if (iproc == itestproc) then
                    print '(5(1x,a,i8))','step',jproc+1,'group:',igrpr(igroup),&
                         ':copying',nvctr_par(jprocsr(4,jproc-1,igroup),igrpr(igroup)),&
                         'elements from',jprocsr(4,jproc-1,igroup),'in',iproc
                 end if
                 ndatac(igroup)=ndatac(igroup)+&
                      nvctr_par(jprocsr(4,jproc-1,igroup),igrpr(igroup)) 
              end if
           end do
        end if

        !send-receive of the results
        ncommsstep2=0
        do igroup=1,ngroupp
           if (jprocsr(3,jproc,igroup) /= -1) then
              ncommsstep2=ncommsstep2+1
              if (jprocsr(3,jproc,igroup) == itestproc) then
                 print '(5(1x,a,i8))','step',jproc+1,'group:',igrpr(igroup),&
                      ': sending',nvctr_par(jprocsr(3,jproc,igroup),igrpr(igroup)),&
                      'elements from',iproc,'to',jprocsr(3,jproc,igroup)
              end if
              ndatas(2,jprocsr(3,jproc,igroup),igrpr(igroup))=ndatas(2,jprocsr(3,jproc,igroup),igrpr(igroup))+&
                   nvctr_par(jprocsr(3,jproc,igroup),igrpr(igroup))
           end if
           if (jprocsr(4,jproc,igroup) /= -1) then
              ncommsstep2=ncommsstep2+1
              if (iproc == itestproc) then
                 print '(5(1x,a,i8))','step',jproc+1,'group:',igrpr(igroup),&
                      ': receiving',nvctr_par(iproc,igrpr(igroup)),&
                      'elements from',jprocsr(4,jproc,igroup),'to',iproc
              end if
              ndatas(2,iproc,igrpr(igroup))=ndatas(2,iproc,igrpr(igroup))-&
                   nvctr_par(iproc,igrpr(igroup))
           end if
        end do

     end do

     if (nproc > 1) call mpiallred(ndatas(1,0,1),2*nproc*ngroup,MPI_SUM,bigdft_mpi%mpi_comm,ierr)
     !if(iproc ==0)print *,'iproc,datas',iproc,ndatas

     do igroup=1,ngroupp
        if (ndatac(igroup) /=0) then
           write(*,*)'ERROR: OP2P communication simulation failed: processor',iproc,&
                ' has calculated',ndatac(igroup),' data more than needed'
           stop
        end if
        if (ndatas(1,iproc,igrpr(igroup)) /=0 .or. ndatas(2,iproc,igrpr(igroup)) /=0) then
           write(*,*)'ERROR: OP2P communication simulation failed: processor',iproc,&
                ' has not a zero balance of send-receive calls',ndatas(1:2,iproc,igrpr(igroup))
           stop
        end if
     end do
  end if
  !stop
  !open(100+iproc)  
  
  call initialize_work_arrays_sumrho(lr,w)
  allocate(psir(lr%d%n1i*lr%d%n2i*lr%d%n3i,orbs%norbp+ndebug),stat=i_stat)
  call memocc(i_stat,psir,'psir',subname)
  
  call razero(lr%d%n1i*lr%d%n2i*lr%d%n3i*orbs%norbp,psir)
  
  !uncompress the wavefunction in the real grid
  do iorb=1,orbs%norbp
     !here ispinor is equal to one
     call daub_to_isf(lr,w,psi(1,1,iorb),psir(1,iorb))
  end do
  
  call deallocate_work_arrays_sumrho(w)
  
  allocate(psiw(lr%d%n1i*lr%d%n2i*lr%d%n3i,maxval(orbs%norb_par(:,0)),2,ngroupp+ndebug),stat=i_stat)
  call memocc(i_stat,psiw,'psiw',subname)
  allocate(dpsiw(lr%d%n1i*lr%d%n2i*lr%d%n3i,maxval(orbs%norb_par(:,0)),3,ngroupp+ndebug),stat=i_stat)
  call memocc(i_stat,dpsiw,'dpsiw',subname)
  !partial densities and potentials
  allocate(rp_ij(lr%d%n1i*lr%d%n2i*lr%d%n3i+ndebug),stat=i_stat)
  call memocc(i_stat,rp_ij,'rp_ij',subname)
  
  !this is the array of the actions of the X potential on psi
  call razero(lr%d%n1i*lr%d%n2i*lr%d%n3i*maxval(orbs%norb_par,1)*2*ngroupp,psiw)
  call razero(lr%d%n1i*lr%d%n2i*lr%d%n3i*maxval(orbs%norb_par,1)*3*ngroupp,dpsiw)
  
  call razero(lr%d%n1i*lr%d%n2i*lr%d%n3i*orbs%norbp,dpsir)

  ncalls=0
  !real communication
  isnow=1
  isnow2=1
  nend=(nproc-1)/2+1
  ncommsstep2=0
  do jproc=0,nend
     irnow=3-isnow
     ncommsstep=0
     !sending receiving data
     do igroup=1,ngroupp
        if (jprocsr(1,jproc,igroup) /= -1) then
           ncommsstep=ncommsstep+1
           if (iprocpm1(1,1,igroup) == itestproc) then
              print *,'step',jproc+1,': sending',nvctr_par(jprocsr(1,jproc,igroup),igrpr(igroup)),&
                   'elements from',iproc,'to',iprocpm1(1,1,igroup)
           end if
           
           if (jproc == 0) then
              call MPI_ISEND(psir(1,iorbgr(2,iproc,igrpr(igroup))),nvctr_par(jprocsr(1,jproc,igroup),igrpr(igroup)),&
                   mpidtypw,iprocpm1(1,1,igroup),&
                   iproc+2*nproc*jproc,bigdft_mpi%mpi_comm,mpireq(ncommsstep),ierr)
           else
              call MPI_ISEND(psiw(1,1,isnow,igroup),nvctr_par(jprocsr(1,jproc,igroup),igrpr(igroup)),&
                   mpidtypw,iprocpm1(1,1,igroup),&
                   iproc+2*nproc*jproc,bigdft_mpi%mpi_comm,mpireq(ncommsstep),ierr)
           end if
        end if
        if (jprocsr(2,jproc,igroup) /= -1) then
           ncommsstep=ncommsstep+1
           if (iproc == itestproc) then
              print *,'step',jproc+1,': receiving',nvctr_par(jprocsr(2,jproc,igroup),igrpr(igroup)),&
                   'elements from',iprocpm1(2,1,igroup),'to',iproc
           end if
           
           call MPI_IRECV(psiw(1,1,irnow,igroup),nvctr_par(jprocsr(2,jproc,igroup),igrpr(igroup)),&
                mpidtypw,iprocpm1(2,1,igroup),&
                iprocpm1(2,1,igroup)+2*nproc*jproc,bigdft_mpi%mpi_comm,mpireq(ncommsstep),ierr)
        end if
     end do
     
     do igroup=1,ngroupp
        if (jproc /= 0 .and. jprocsr(3,jproc,igroup) /= -1) then
           !put to zero the sending element
           call razero(lr%d%n1i*lr%d%n2i*lr%d%n3i*maxval(orbs%norb_par,1),dpsiw(1,1,3,igroup))
        end if
     end do
     
     !calculation for orbitals to be performed
     do igroup=1,ngroupp
        if (jproc == 0) then
           doit=.true.
        else if (jprocsr(2,jproc-1,igroup) /=-1) then
           doit=.true.
        else
           doit=.false.
        end if
        if (doit) then
           !calculation of the partial densities and potentials
           !starting point of the loop
           !here there is the calculation routine
           !number of orbitals to be treated locally
           norbi=nvctr_par(iproc,igrpr(igroup))/(lr%d%n1i*lr%d%n2i*lr%d%n3i)
           if (jproc == 0) then
              norbj=norbi
           else
              norbj=nvctr_par(jprocsr(2,jproc-1,igroup),igrpr(igroup))/(lr%d%n1i*lr%d%n2i*lr%d%n3i)
           end if
           !calculating the starting orbitals locally
           iorbs=iorbgr(2,iproc,igrpr(igroup))
           if (jproc == 0) then
              jorbs=iorbs
           else
              jorbs=iorbgr(2,jprocsr(2,jproc-1,igroup),igrpr(igroup))
           end if
           !calculate the starting orbital globally
           isorb=iorbgr(1,iproc,igrpr(igroup))
           if (jproc==0) then
              jsorb=isorb
           else       
              jsorb=iorbgr(1,jprocsr(2,jproc-1,igroup),igrpr(igroup))
              !if (igrpr(igroup) == 2) jsorb=orbs%norbu+jsorb
           end if
           
           !loop over all the orbitals
           !for the first step do only the upper triangular part
           do iorb=iorbs,iorbs+norbi-1
              hfacj=-sfac*orbs%occup(iorb+isorb)
              do jorb=jorbs,jorbs+norbj-1
                 !first cross-check whether the spin indices are the same
                 if (orbs%spinsgn(isorb+iorb) /= orbs%spinsgn(jsorb+jorb)) then
                    write(*,*)'ERROR in partitioning the orbitals',&
                         iorb+isorb,jorb+jsorb,igroup,jsorb,iproc
                    stop
                 end if
                 hfaci=-sfac*orbs%occup(jorb+jsorb)
                 !do it only for upper triangular results
                 if (jproc /= 0 .or. jorb+jsorb >= iorb+isorb) then
                    if (jproc == 0 ) then
                       do i=1,lr%d%n1i*lr%d%n2i*lr%d%n3i
                          rp_ij(i)=hfac*psir(i,iorb)*psir(i,jorb)
                       end do
                    else
                       do i=1,lr%d%n1i*lr%d%n2i*lr%d%n3i
                          rp_ij(i)=hfac*psir(i,iorb)*psiw(i,jorb-jorbs+1,isnow,igroup)
                       end do
                    end if
                    ncalls=ncalls+1                    
                    !Poisson solver in sequential
                    if (iproc == iprocref .and. verbose > 1) then
                          call yaml_comment('Exact exchange calculation: ' // trim(yaml_toa( &
                               nint(real(ncalls,gp)/real(ncalltot,gp)*100.0_gp),fmt='(i3)')) //'%')
                          !write(*,'(1x,a,i3,a2)')'Exact exchange calculation: ',&
                          !     nint(real(ncalls,gp)/real(ncalltot,gp)*100.0_gp),' %'
                       !write(*,'(1x,a,2(1x,i5))')'Exact exchange calculation: ',ncalls,ncalltot
                       !write(*,*)'Exact exchange calculation: spin, orbitals:',igrpr(igroup),iorb,jorb
                    end if

                    call H_potential('D',pkernel,rp_ij,rp_ij,ehart,0.0_dp,.false.,&
                         quiet='YES')
                    
                    !this factor is only valid with one k-point
                    !can be easily generalised to the k-point case
                    hfac2=sfac*orbs%occup(iorb+isorb)*orbs%occup(jorb+jsorb)
                    
                    !exact exchange energy
                    if (iorb+isorb == jorb+jsorb) then
                       eexctX=eexctX+hfac2*real(ehart,gp)
                    else
                       !if the result has to be sent away
                       if (jprocsr(3,jproc,igroup) /= -1 .or. jproc==0) then
                          eexctX=eexctX+2.0_gp*hfac2*real(ehart,gp)
                       else !otherwise other processors are already calculating it
                          eexctX=eexctX+hfac2*real(ehart,gp)
                       end if
                    end if
                    !accumulate the results for each of the wavefunctions concerned
                    if (jproc == 0) then
                       do i=1,lr%d%n1i*lr%d%n2i*lr%d%n3i
                          dpsir(i,iorb)=dpsir(i,iorb)+&
                               hfaci*rp_ij(i)*psir(i,jorb)
                       end do
                       if (jorb+jsorb /= iorb+isorb) then
                          do i=1,lr%d%n1i*lr%d%n2i*lr%d%n3i
                             dpsir(i,jorb)=dpsir(i,jorb)+&
                                  hfacj*rp_ij(i)*psir(i,iorb)
                          end do
                          !write(100+iproc,*)jorb+jsorb,iorb+isorb,igrpr(igroup) 
                       end if
                    else
                       do i=1,lr%d%n1i*lr%d%n2i*lr%d%n3i
                          dpsir(i,iorb)=dpsir(i,iorb)+&
                               hfaci*rp_ij(i)*psiw(i,jorb-jorbs+1,isnow,igroup)
                       end do
                    end if
                    !write(100+iproc,*)iorb+isorb,jorb+jsorb,igrpr(igroup)
                 end if
                 
                 !fill the set of the vector to be sent to the other processes
                 !in the first step the results are self-contained
                 if (jproc /= 0 .and. jprocsr(3,jproc,igroup) /= -1) then
                    !write(100+iproc,*)jorb+jsorb,iorb+isorb,igrpr(igroup) 
                    do i=1,lr%d%n1i*lr%d%n2i*lr%d%n3i
                       dpsiw(i,jorb-jorbs+1,3,igroup)=dpsiw(i,jorb-jorbs+1,3,igroup)+&
                            hfacj*rp_ij(i)*psir(i,iorb)
                    end do
                 end if
              end do
           end do
        end if
     end do
         
     if (ncommsstep2 > 0) then
        !verify that the messages have been passed
        call MPI_WAITALL(ncommsstep2,mpireq2,mpistat2,ierr)
        if (ierr /=0)  print *,'step2,ierr',jproc+1,iproc,ierr,mpistat2 !,MPI_STATUSES_IGNORE
        
        !copy the results which have been received (the messages sending are after)
        do igroup=1,ngroupp
           if (jprocsr(4,jproc-1,igroup) /= -1) then
              if (iproc == itestproc) then
                 print '(5(1x,a,i8))','step',jproc+1,'group:',igrpr(igroup),&
                      ':copying',nvctr_par(jprocsr(4,jproc-1,igroup),igrpr(igroup)),&
                      'processed elements from',jprocsr(4,jproc-1,igroup),'in',iproc
              end if
              
              call axpy(nvctr_par(iproc,igrpr(igroup)),1.0_wp,dpsiw(1,1,irnow2,igroup),1,&
                   dpsir(1,iorbgr(2,iproc,igrpr(igroup))),1)
           end if
        end do
     end if
     
     ncommsstep2=0
     !meanwhile, we can receive the result from the processor which has the psi 
     irnow2=3-isnow2
     do igroup=1,ngroupp
        if (jprocsr(3,jproc,igroup) /= -1) then
           ncommsstep2=ncommsstep2+1
           if (jprocsr(3,jproc,igroup) == itestproc) then
              print '(5(1x,a,i8))','step',jproc+1,'group:',igrpr(igroup),&
                   ': sending',nvctr_par(jprocsr(3,jproc,igroup),igrpr(igroup)),&
                   'elements from',iproc,'to',jprocsr(3,jproc,igroup)
           end if
           call dcopy(nvctr_par(jprocsr(3,jproc,igroup),igrpr(igroup)),&
                dpsiw(1,1,3,igroup),1,dpsiw(1,1,isnow2,igroup),1)
           
           call MPI_ISEND(dpsiw(1,1,isnow2,igroup),&
                nvctr_par(jprocsr(3,jproc,igroup),igrpr(igroup)),mpidtypw,&
                jprocsr(3,jproc,igroup),&
                iproc+nproc+2*nproc*jproc,bigdft_mpi%mpi_comm,mpireq2(ncommsstep2),ierr)
        end if
        if (jprocsr(4,jproc,igroup) /= -1) then
           ncommsstep2=ncommsstep2+1
           if (iproc == itestproc) then
              print '(5(1x,a,i8))','step',jproc+1,'group:',igrpr(igroup),&
                   ': receiving',nvctr_par(iproc,igrpr(igroup)),&
                   'elements from',jprocsr(4,jproc,igroup),'to',iproc
           end if
           call MPI_IRECV(dpsiw(1,1,irnow2,igroup),&
                nvctr_par(iproc,igrpr(igroup)),mpidtypw,jprocsr(4,jproc,igroup),&
                jprocsr(4,jproc,igroup)+nproc+2*nproc*jproc,bigdft_mpi%mpi_comm,mpireq2(ncommsstep2),ierr)
           
        end if
     end do
     if (jproc>1) isnow2=3-isnow2
     
     if (ncommsstep /=0) then
        !verify that the messages have been passed
        !print *,'waiting,iproc',iproc
        call MPI_WAITALL(ncommsstep,mpireq,mpistat,ierr)
        if (ierr /=0) print *,'step,ierr',jproc+1,iproc,ierr,mpistat !,MPI_STATUSES_IGNORE
        !print *,'done,iproc',iproc
     end if
     isnow=3-isnow
     ncommsstep=0
  end do
  
  !call MPI_BARRIER(bigdft_mpi%mpi_comm,ierr)
  if (nproc>1) call mpiallred(eexctX,1,MPI_SUM,bigdft_mpi%mpi_comm,ierr)
  
  !the exact exchange energy is half the Hartree energy (which already has another half)
  eexctX=-exctXfac*eexctX
  
  if (iproc == 0) call yaml_map('Exact Exchange Energy',eexctX,fmt='(1pe18.11)')
  !if (iproc == 0) write(*,'(1x,a,1x,1pe18.11)')'Exact Exchange Energy:',eexctX

  !close(100+iproc)
  i_all=-product(shape(nvctr_par))*kind(nvctr_par)
  deallocate(nvctr_par,stat=i_stat)
  call memocc(i_stat,i_all,'nvctr_par',subname)
  
  i_all=-product(shape(iorbgr))*kind(iorbgr)
  deallocate(iorbgr,stat=i_stat)
  call memocc(i_stat,i_all,'iorbgr',subname)


  i_all=-product(shape(ndatas))*kind(ndatas)
  deallocate(ndatas,stat=i_stat)
  call memocc(i_stat,i_all,'ndatas',subname)

  i_all=-product(shape(ndatac))*kind(ndatac)
  deallocate(ndatac,stat=i_stat)
  call memocc(i_stat,i_all,'ndatac',subname)

  i_all=-product(shape(rp_ij))*kind(rp_ij)
  deallocate(rp_ij,stat=i_stat)
  call memocc(i_stat,i_all,'rp_ij',subname)
  
  i_all=-product(shape(psiw))*kind(psiw)
  deallocate(psiw,stat=i_stat)
  call memocc(i_stat,i_all,'psiw',subname)

  i_all=-product(shape(dpsiw))*kind(dpsiw)
  deallocate(dpsiw,stat=i_stat)
  call memocc(i_stat,i_all,'dpsiw',subname)


  i_all=-product(shape(psir))*kind(psir)
  deallocate(psir,stat=i_stat)
  call memocc(i_stat,i_all,'psir',subname)

  i_all=-product(shape(igrpr))*kind(igrpr)
  deallocate(igrpr,stat=i_stat)
  call memocc(i_stat,i_all,'igrpr',subname)

  i_all=-product(shape(iprocpm1))*kind(iprocpm1)
  deallocate(iprocpm1,stat=i_stat)
  call memocc(i_stat,i_all,'iprocpm1',subname)

  i_all=-product(shape(jprocsr))*kind(jprocsr)
  deallocate(jprocsr,stat=i_stat)
  call memocc(i_stat,i_all,'jprocsr',subname)

  !call timing(iproc,'Exchangecorr  ','OF')

END SUBROUTINE exact_exchange_potential_round


!!$!> Calculate the exact exchange potential on occupied orbitals
!!$!! within the symmetric round-robin scheme
!!$!! the psi is already given in the real-space form
!!$subroutine exact_exchange_potential_round_new(iproc,nproc,geocode,nspin,lr,orbs,&
!!$     hxh,hyh,hzh,pkernel,psi,dpsir,eexctX)
!!$  use module_base
!!$  use module_types
!!$  use Poisson_Solver
!!$  use module_xc
!!$  implicit none
!!$  character(len=1), intent(in) :: geocode !< @copydoc poisson_solver::coulomb_operator::geocode
!!$  integer, intent(in) :: iproc,nproc,nspin
!!$  real(gp), intent(in) :: hxh,hyh,hzh
!!$  type(locreg_descriptors), intent(in) :: lr
!!$  type(orbitals_data), intent(in) :: orbs
!!$  real(wp), dimension(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f,orbs%nspinor,orbs%norbp), intent(in) :: psi
!!$  real(dp), dimension(*), intent(in) :: pkernel
!!$  real(gp), intent(out) :: eexctX
!!$  real(wp), dimension(lr%d%n1i*lr%d%n2i*lr%d%n3i,orbs%norbp), intent(out) :: dpsir
!!$  !local variables
!!$  character(len=*), parameter :: subname='exact_exchange_potential_round_new'
!!$  logical :: doit
!!$  integer :: i_all,i_stat,ierr,ispin,ncommsstep,ncommsstep2,isnow,irnow,isnow2,irnow2,jsorb,kproc,norbp,jgroup
!!$  integer :: i,iorb,jorb,jproc,igroup,ngroup,ngroupp,jprocsend,jprocrecv,jprocrecv2,nend,isorb,iorbs,iorbe,jorbs,jorbe
!!$  integer :: icount,nprocgr,iprocgrs,iprocgrr,itestproc,norbi,norbj,iproclast,ncalltot,icountmax,iprocref,ncalls
!!$  real(gp) :: ehart,hfac,exctXfac,sfac,hfaci,hfacj,hfac2
!!$  integer, dimension(4) :: mpireq,mpireq2
!!$  integer, dimension(MPI_STATUS_SIZE,4) :: mpistat,mpistat2
!!$  type(workarr_sumrho) :: w
!!$  integer, dimension(:), allocatable :: igrpr,ndatac
!!$  integer, dimension(:,:), allocatable :: nvctr_par
!!$  integer, dimension(:,:,:), allocatable :: jprocsr,iprocpm1,ndatas,iorbgr
!!$  real(wp), dimension(:), allocatable :: rp_ij
!!$  real(wp), dimension(:,:), allocatable :: psir
!!$  real(wp), dimension(:,:,:,:), allocatable :: psiw,dpsiw
!!$
!!$  !call timing(iproc,'Exchangecorr  ','ON')
!!$
!!$  exctXfac = xc_exctXfac()
!!$
!!$  eexctX=0.0_gp
!!$
!!$  !build the partial densities for the poisson solver, calculate the partial potential
!!$  !and accumulate the result
!!$  !do it for different spins
!!$  !for non spin-polarised systems there is a factor of two
!!$  !non-collinear spin not yet implemented
!!$  if (nspin==2) then
!!$     sfac=1.0_gp
!!$     ngroup=2
!!$  else 
!!$     sfac=0.5_gp
!!$     ngroup=1
!!$  end if
!!$
!!$  hfac=1.0_gp/(hxh*hyh*hzh)
!!$
!!$  !here we can start with the round-robin scheme
!!$  !since the orbitlas are all occupied we have to use the symmetric scheme
!!$  !we have first to define the number of groups, which correspond to the repartition 
!!$  !of spin up and spin down orbitals
!!$  allocate(nvctr_par(0:nproc-1,ngroup+ndebug),stat=i_stat)
!!$  call memocc(i_stat,nvctr_par,'nvctr_par',subname)
!!$
!!$  allocate(iorbgr(2,0:nproc-1,ngroup+ndebug),stat=i_stat)
!!$  call memocc(i_stat,iorbgr,'iorbgr',subname)
!!$
!!$
!!$  
!!$  if (ngroup==2) then
!!$     isorb=0
!!$     do jproc=0,nproc-1
!!$        iorbgr(1,jproc,1)=isorb
!!$        iorbgr(2,jproc,1)=1
!!$        norbp=max(min(isorb+orbs%norb_par(jproc),orbs%norbu)-isorb,0)
!!$        if (norbp == 0) then
!!$           iorbgr(1,jproc,1)=0
!!$           iorbgr(1,jproc,2)=0
!!$        end if
!!$        nvctr_par(jproc,1)=norbp*lr%d%n1i*lr%d%n2i*lr%d%n3i
!!$        iorbgr(1,jproc,2)=isorb
!!$        iorbgr(2,jproc,2)=norbp+1
!!$        norbp=max(isorb+orbs%norb_par(jproc)-max(orbs%norbu,isorb),0)
!!$        if (norbp == 0) then
!!$           iorbgr(1,jproc,2)=0
!!$           iorbgr(2,jproc,2)=1
!!$        end if
!!$        nvctr_par(jproc,2)=norbp*lr%d%n1i*lr%d%n2i*lr%d%n3i
!!$        isorb=isorb+orbs%norb_par(jproc)
!!$     end do
!!$     if (iproc ==0) then
!!$        print '(a,10(1x,i8))','iproc,nvctr_parA',iproc,nvctr_par(:,1)/(lr%d%n1i*lr%d%n2i*lr%d%n3i)
!!$        print '(a,10(1x,i8))','iproc,nvctr_parB',iproc,nvctr_par(:,2)/(lr%d%n1i*lr%d%n2i*lr%d%n3i)
!!$        print '(a,10(1x,i8))','iproc,iorbgr',iproc,iorbgr
!!$     end if
!!$  else
!!$     isorb=0
!!$     do jproc=0,nproc-1
!!$        iorbgr(1,jproc,1)=isorb
!!$        iorbgr(2,jproc,1)=1
!!$        nvctr_par(jproc,1)=orbs%norb_par(jproc)*lr%d%n1i*lr%d%n2i*lr%d%n3i
!!$        isorb=isorb+orbs%norb_par(jproc)
!!$     end do
!!$  end if
!!$
!!$
!!$subroutine op2p_communication_descriptors(iproc,nproc,ngroup,nvctr_par,op2p)
!!$  use module_base
!!$  use module_types
!!$  implicit none
!!$  integer, intent(in) :: iproc,nproc,ngroup
!!$  integer, dimension(0:nproc-1,ngroup), intent(in) :: nvctr_par
!!$  type(op2p_descriptors), intent(out) :: op2p
!!$  !local variables
!!$  integer :: igroup
!!$
!!$  !allocate and copy nvctr_par
!!$  allocate(op2p%nvctr_par(0:nproc-1,ngroup+ndebug),stat=i_stat)
!!$  call memocc(i_stat,op2p%nvctr_par,'nvctr_par',subname)
!!$  
!!$  do igroup=1,ngroup
!!$     do kproc=0,nproc-1
!!$        op2p%nvctr_par(kproc,igroup)=nvctr_par(kproc,igroup)
!!$     end do
!!$  end do
!!$
!!$  !here we can allocate the working arrays giving the maximum
!!$  !between the components for each group
!!$  op2p%ngroupp=0
!!$  do igroup=1,ngroup
!!$     if (nvctr_par(iproc,igroup) > 0) then
!!$        op2p%ngroupp=op2p%ngroupp+1
!!$     end if
!!$  end do
!!$
!!$  !determine the array of the groups which are of interest for this processor
!!$  allocate(op2p%igrpr(op2p%ngroupp+ndebug),stat=i_stat)
!!$  call memocc(i_stat,igrpr,'igrpr',subname)
!!$  allocate(op2p%iprocpm1(2,0:nproc-1,op2p%ngroupp+ndebug),stat=i_stat)
!!$  call memocc(i_stat,iprocpm1,'iprocpm1',subname)
!!$
!!$
!!$  !determine for each processor the groups which has to be used
!!$  icount=0
!!$  do igroup=1,ngroup
!!$     if (nvctr_par(iproc,igroup) > 0) then
!!$        icount=icount+1
!!$        op2p%igrpr(icount)=igroup
!!$     end if
!!$  end do
!!$
!!$  !calculate the processor which lies after and before the present in the list
!!$  op2p%iprocpm1=-1
!!$  do igroup=1,op2p%ngroupp
!!$     iprocgrs=-1
!!$     iprocgrr=-1
!!$     do kproc=0,nproc-1
!!$        if (nvctr_par(modulo(iproc+kproc,nproc),op2p%igrpr(igroup)) > 0) then
!!$           iprocgrs=iprocgrs+1
!!$           op2p%iprocpm1(1,iprocgrs,igroup)=modulo(iproc+kproc,nproc)
!!$        end if
!!$        if (op2p%nvctr_par(modulo(iproc-kproc,nproc),op2p%igrpr(igroup)) > 0) then
!!$           iprocgrr=iprocgrr+1
!!$           op2p%iprocpm1(2,iprocgrr,igroup)=modulo(iproc-kproc,nproc)
!!$        end if
!!$     end do
!!$  end do
!!$
!!$  !find the processor whih has the maximum number of groups
!!$  icountmax=0
!!$  do kproc=0,nproc-1
!!$     icount=0
!!$     do igroup=1,ngroup
!!$        if (nvctr_par(kproc,igroup) > 0) then
!!$           icount=icount+1
!!$        end if
!!$     end do
!!$     if (icount > icountmax) then
!!$        op2p%iprocref=kproc
!!$        icountmax=icount
!!$     end if
!!$  end do
!!$
!!$  !calculate the list of send-receive operations which have to be performed per group
!!$  !allocate it at the maximum size needed
!!$  allocate(op2p%jprocsr(4,0:nproc/2+1,op2p%ngroupp+ndebug),stat=i_stat)
!!$  call memocc(i_stat,jprocsr,'jprocsr',subname)
!!$  !initalise array to minus one
!!$  jprocsr=-1
!!$
!!$  do igroup=1,op2p%ngroupp
!!$     !calculate the number of processors per group
!!$     nprocgr=0
!!$     do kproc=0,nproc-1
!!$        if (nvctr_par(kproc,op2p%igrpr(igroup)) > 0) nprocgr=nprocgr+1
!!$     end do
!!$     
!!$     !do not sent anything if there is only one member in the group
!!$     if (nprocgr > 1) then
!!$        do kproc=0,(nprocgr-1)/2-1
!!$           !define the arrays for send-receive of data
!!$           op2p%jprocsr(1,kproc,igroup)=op2p%iprocpm1(2,kproc,igroup)
!!$           op2p%jprocsr(2,kproc,igroup)=op2p%iprocpm1(2,kproc+1,igroup)
!!$           if (iproc == iprocref) then
!!$              ncalltot=ncalltot+&
!!$                   (nvctr_par(jprocsr(2,kproc,igroup),igrpr(igroup))/(lr%d%n1i*lr%d%n2i*lr%d%n3i))*&
!!$                   (nvctr_par(iproc,igrpr(igroup))/(lr%d%n1i*lr%d%n2i*lr%d%n3i))
!!$             end if
!!$           if (kproc > 0) then
!!$              op2p%jprocsr(3,kproc,igroup)=op2p%iprocpm1(2,kproc,igroup)
!!$              op2p%jprocsr(4,kproc,igroup)=op2p%iprocpm1(1,kproc,igroup)
!!$           end if
!!$        end do
!!$        kproc=(nprocgr-1)/2
!!$        !the last step behaves differently if the group number is odd or even
!!$        if (modulo(nprocgr,2) == 0) then
!!$           jprocsr(1,kproc,igroup)=iprocpm1(2,kproc,igroup)
!!$           jprocsr(2,kproc,igroup)=iprocpm1(2,kproc+1,igroup)
!!$           if (iproc == iprocref) then
!!$              ncalltot=ncalltot+&
!!$                   (nvctr_par(jprocsr(2,kproc,igroup),igrpr(igroup))/(lr%d%n1i*lr%d%n2i*lr%d%n3i))*&
!!$                   (nvctr_par(iproc,igrpr(igroup))/(lr%d%n1i*lr%d%n2i*lr%d%n3i))
!!$             end if
!!$           if (kproc > 0) then
!!$              jprocsr(3,kproc,igroup)=iprocpm1(2,kproc,igroup)
!!$              jprocsr(4,kproc,igroup)=iprocpm1(1,kproc,igroup)
!!$           end if
!!$        else
!!$           jprocsr(3,kproc,igroup)=iprocpm1(2,kproc,igroup)
!!$           jprocsr(4,kproc,igroup)=iprocpm1(1,kproc,igroup)
!!$        end if
!!$     end if
!!$  end do
!!$
!!$  !stop
!!$  !open(100+iproc)  
!!$  
!!$  call initialize_work_arrays_sumrho(lr,w)
!!$  allocate(psir(lr%d%n1i*lr%d%n2i*lr%d%n3i,orbs%norbp+ndebug),stat=i_stat)
!!$  call memocc(i_stat,psir,'psir',subname)
!!$  
!!$  call razero(lr%d%n1i*lr%d%n2i*lr%d%n3i*orbs%norbp,psir)
!!$  
!!$  !uncompress the wavefunction in the real grid
!!$  do iorb=1,orbs%norbp
!!$     !here ispinor is equal to one
!!$     call daub_to_isf(lr,w,psi(1,1,iorb),psir(1,iorb))
!!$  end do
!!$  
!!$  call deallocate_work_arrays_sumrho(w)
!!$  
!!$  allocate(psiw(lr%d%n1i*lr%d%n2i*lr%d%n3i,maxval(orbs%norb_par),2,ngroupp+ndebug),stat=i_stat)
!!$  call memocc(i_stat,psiw,'psiw',subname)
!!$  allocate(dpsiw(lr%d%n1i*lr%d%n2i*lr%d%n3i,maxval(orbs%norb_par),3,ngroupp+ndebug),stat=i_stat)
!!$  call memocc(i_stat,dpsiw,'dpsiw',subname)
!!$  !partial densities and potentials
!!$  allocate(rp_ij(lr%d%n1i*lr%d%n2i*lr%d%n3i+ndebug),stat=i_stat)
!!$  call memocc(i_stat,rp_ij,'rp_ij',subname)
!!$  
!!$  !this is the array of the actions of the X potential on psi
!!$  call razero(lr%d%n1i*lr%d%n2i*lr%d%n3i*maxval(orbs%norb_par)*2*ngroupp,psiw)
!!$  call razero(lr%d%n1i*lr%d%n2i*lr%d%n3i*maxval(orbs%norb_par)*3*ngroupp,dpsiw)
!!$  
!!$  call razero(lr%d%n1i*lr%d%n2i*lr%d%n3i*orbs%norbp,dpsir)
!!$
!!$  ncalls=0
!!$  !real communication
!!$  isnow=1
!!$  isnow2=1
!!$  nend=(nproc-1)/2+1
!!$  ncommsstep2=0
!!$  do jproc=0,nend
!!$     irnow=3-isnow
!!$     ncommsstep=0
!!$     !sending receiving data
!!$     do igroup=1,ngroupp
!!$        if (jprocsr(1,jproc,igroup) /= -1) then
!!$           ncommsstep=ncommsstep+1
!!$           if (iprocpm1(1,1,igroup) == itestproc) then
!!$              print *,'step',jproc+1,': sending',nvctr_par(jprocsr(1,jproc,igroup),igrpr(igroup)),&
!!$                   'elements from',iproc,'to',iprocpm1(1,1,igroup)
!!$           end if
!!$           
!!$           if (jproc == 0) then
!!$              call MPI_ISEND(psir(1,iorbgr(2,iproc,igrpr(igroup))),nvctr_par(jprocsr(1,jproc,igroup),igrpr(igroup)),&
!!$                   mpidtypw,iprocpm1(1,1,igroup),&
!!$                   iproc+2*nproc*jproc,bigdft_mpi%mpi_comm,mpireq(ncommsstep),ierr)
!!$           else
!!$              call MPI_ISEND(psiw(1,1,isnow,igroup),nvctr_par(jprocsr(1,jproc,igroup),igrpr(igroup)),&
!!$                   mpidtypw,iprocpm1(1,1,igroup),&
!!$                   iproc+2*nproc*jproc,bigdft_mpi%mpi_comm,mpireq(ncommsstep),ierr)
!!$           end if
!!$        end if
!!$        if (jprocsr(2,jproc,igroup) /= -1) then
!!$           ncommsstep=ncommsstep+1
!!$           if (iproc == itestproc) then
!!$              print *,'step',jproc+1,': receiving',nvctr_par(jprocsr(2,jproc,igroup),igrpr(igroup)),&
!!$                   'elements from',iprocpm1(2,1,igroup),'to',iproc
!!$           end if
!!$           
!!$           call MPI_IRECV(psiw(1,1,irnow,igroup),nvctr_par(jprocsr(2,jproc,igroup),igrpr(igroup)),&
!!$                mpidtypw,iprocpm1(2,1,igroup),&
!!$                iprocpm1(2,1,igroup)+2*nproc*jproc,bigdft_mpi%mpi_comm,mpireq(ncommsstep),ierr)
!!$        end if
!!$     end do
!!$     
!!$     do igroup=1,ngroupp
!!$        if (jproc /= 0 .and. jprocsr(3,jproc,igroup) /= -1) then
!!$           !put to zero the sending element
!!$           call razero(lr%d%n1i*lr%d%n2i*lr%d%n3i*maxval(orbs%norb_par),dpsiw(1,1,3,igroup))
!!$        end if
!!$     end do
!!$     
!!$     !calculation for orbitals to be performed
!!$     do igroup=1,ngroupp
!!$        if (jproc == 0) then
!!$           doit=.true.
!!$        else if (jprocsr(2,jproc-1,igroup) /=-1) then
!!$           doit=.true.
!!$        else
!!$           doit=.false.
!!$        end if
!!$        if (doit) then
!!$           !calculation of the partial densities and potentials
!!$           !starting point of the loop
!!$           !here there is the calculation routine
!!$           !number of orbitals to be treated locally
!!$           norbi=nvctr_par(iproc,igrpr(igroup))/(lr%d%n1i*lr%d%n2i*lr%d%n3i)
!!$           if (jproc == 0) then
!!$              norbj=norbi
!!$           else
!!$              norbj=nvctr_par(jprocsr(2,jproc-1,igroup),igrpr(igroup))/(lr%d%n1i*lr%d%n2i*lr%d%n3i)
!!$           end if
!!$           !calculating the starting orbitals locally
!!$           iorbs=iorbgr(2,iproc,igrpr(igroup))
!!$           if (jproc == 0) then
!!$              jorbs=iorbs
!!$           else
!!$              jorbs=iorbgr(2,jprocsr(2,jproc-1,igroup),igrpr(igroup))
!!$           end if
!!$           !calculate the starting orbital globally
!!$           isorb=iorbgr(1,iproc,igrpr(igroup))
!!$           if (jproc==0) then
!!$              jsorb=isorb
!!$           else       
!!$              jsorb=iorbgr(1,jprocsr(2,jproc-1,igroup),igrpr(igroup))
!!$              !if (igrpr(igroup) == 2) jsorb=orbs%norbu+jsorb
!!$           end if
!!$           
!!$           !loop over all the orbitals
!!$           !for the first step do only the upper triangular part
!!$           do iorb=iorbs,iorbs+norbi-1
!!$              hfacj=-sfac*orbs%occup(iorb+isorb)
!!$              do jorb=jorbs,jorbs+norbj-1
!!$                 !first cross-check whether the spin indices are the same
!!$                 if (orbs%spinsgn(isorb+iorb) /= orbs%spinsgn(jsorb+jorb)) then
!!$                    write(*,*)'ERROR in partitioning the orbitals',&
!!$                         iorb+isorb,jorb+jsorb,igroup,jsorb,iproc
!!$                    stop
!!$                 end if
!!$                 hfaci=-sfac*orbs%occup(jorb+jsorb)
!!$                 !do it only for upper triangular results
!!$                 if (jproc /= 0 .or. jorb+jsorb >= iorb+isorb) then
!!$                    if (jproc == 0 ) then
!!$                       do i=1,lr%d%n1i*lr%d%n2i*lr%d%n3i
!!$                          rp_ij(i)=hfac*psir(i,iorb)*psir(i,jorb)
!!$                       end do
!!$                    else
!!$                       do i=1,lr%d%n1i*lr%d%n2i*lr%d%n3i
!!$                          rp_ij(i)=hfac*psir(i,iorb)*psiw(i,jorb-jorbs+1,isnow,igroup)
!!$                       end do
!!$                    end if
!!$                    ncalls=ncalls+1                    
!!$                    !Poisson solver in sequential
!!$                    if (iproc == iprocref .and. verbose > 1) then
!!$                          write(*,'(1x,a,i3,a2)')'Exact exchange calculation: ',&
!!$                               nint(real(ncalls,gp)/real(ncalltot,gp)*100.0_gp),' %'
!!$                       !write(*,'(1x,a,2(1x,i5))')'Exact exchange calculation: ',ncalls,ncalltot
!!$                       !write(*,*)'Exact exchange calculation: spin, orbitals:',igrpr(igroup),iorb,jorb
!!$                    end if
!!$
!!$                    call H_potential(geocode,'D',0,1,&
!!$                         lr%d%n1i,lr%d%n2i,lr%d%n3i,hxh,hyh,hzh,&
!!$                         rp_ij,pkernel,rp_ij,ehart,0.0_dp,.false.,&
!!$                         quiet='YES')
!!$                    
!!$                    !this factor is only valid with one k-point
!!$                    !can be easily generalised to the k-point case
!!$                    hfac2=sfac*orbs%occup(iorb+isorb)*orbs%occup(jorb+jsorb)
!!$                    
!!$                    !exact exchange energy
!!$                    if (iorb+isorb == jorb+jsorb) then
!!$                       eexctX=eexctX+hfac2*real(ehart,gp)
!!$                    else
!!$                       !if the result has to be sent away
!!$                       if (jprocsr(3,jproc,igroup) /= -1 .or. jproc==0) then
!!$                          eexctX=eexctX+2.0_gp*hfac2*real(ehart,gp)
!!$                       else !otherwise other processors are already calculating it
!!$                          eexctX=eexctX+hfac2*real(ehart,gp)
!!$                       end if
!!$                    end if
!!$                    !accumulate the results for each of the wavefunctions concerned
!!$                    if (jproc == 0) then
!!$                       do i=1,lr%d%n1i*lr%d%n2i*lr%d%n3i
!!$                          dpsir(i,iorb)=dpsir(i,iorb)+&
!!$                               hfaci*rp_ij(i)*psir(i,jorb)
!!$                       end do
!!$                       if (jorb+jsorb /= iorb+isorb) then
!!$                          do i=1,lr%d%n1i*lr%d%n2i*lr%d%n3i
!!$                             dpsir(i,jorb)=dpsir(i,jorb)+&
!!$                                  hfacj*rp_ij(i)*psir(i,iorb)
!!$                          end do
!!$                          !write(100+iproc,*)jorb+jsorb,iorb+isorb,igrpr(igroup) 
!!$                       end if
!!$                    else
!!$                       do i=1,lr%d%n1i*lr%d%n2i*lr%d%n3i
!!$                          dpsir(i,iorb)=dpsir(i,iorb)+&
!!$                               hfaci*rp_ij(i)*psiw(i,jorb-jorbs+1,isnow,igroup)
!!$                       end do
!!$                    end if
!!$                    !write(100+iproc,*)iorb+isorb,jorb+jsorb,igrpr(igroup)
!!$                 end if
!!$                 
!!$                 !fill the set of the vector to be sent to the other processes
!!$                 !in the first step the results are self-contained
!!$                 if (jproc /= 0 .and. jprocsr(3,jproc,igroup) /= -1) then
!!$                    !write(100+iproc,*)jorb+jsorb,iorb+isorb,igrpr(igroup) 
!!$                    do i=1,lr%d%n1i*lr%d%n2i*lr%d%n3i
!!$                       dpsiw(i,jorb-jorbs+1,3,igroup)=dpsiw(i,jorb-jorbs+1,3,igroup)+&
!!$                            hfacj*rp_ij(i)*psir(i,iorb)
!!$                    end do
!!$                 end if
!!$              end do
!!$           end do
!!$        end if
!!$     end do
!!$         
!!$     if (ncommsstep2 > 0) then
!!$        !verify that the messages have been passed
!!$        call MPI_WAITALL(ncommsstep2,mpireq2,mpistat2,ierr)
!!$        if (ierr /=0)  print *,'step2,ierr',jproc+1,iproc,ierr,mpistat,MPI_STATUSES_IGNORE
!!$        
!!$        !copy the results which have been received (the messages sending are after)
!!$        do igroup=1,ngroupp
!!$           if (jprocsr(4,jproc-1,igroup) /= -1) then
!!$              if (iproc == itestproc) then
!!$                 print '(5(1x,a,i8))','step',jproc+1,'group:',igrpr(igroup),&
!!$                      ':copying',nvctr_par(jprocsr(4,jproc-1,igroup),igrpr(igroup)),&
!!$                      'processed elements from',jprocsr(4,jproc-1,igroup),'in',iproc
!!$              end if
!!$              
!!$              call axpy(nvctr_par(iproc,igrpr(igroup)),1.0_wp,dpsiw(1,1,irnow2,igroup),1,&
!!$                   dpsir(1,iorbgr(2,iproc,igrpr(igroup))),1)
!!$           end if
!!$        end do
!!$     end if
!!$     
!!$     ncommsstep2=0
!!$     !meanwhile, we can receive the result from the processor which has the psi 
!!$     irnow2=3-isnow2
!!$     do igroup=1,ngroupp
!!$        if (jprocsr(3,jproc,igroup) /= -1) then
!!$           ncommsstep2=ncommsstep2+1
!!$           if (jprocsr(3,jproc,igroup) == itestproc) then
!!$              print '(5(1x,a,i8))','step',jproc+1,'group:',igrpr(igroup),&
!!$                   ': sending',nvctr_par(jprocsr(3,jproc,igroup),igrpr(igroup)),&
!!$                   'elements from',iproc,'to',jprocsr(3,jproc,igroup)
!!$           end if
!!$           call dcopy(nvctr_par(jprocsr(3,jproc,igroup),igrpr(igroup)),&
!!$                dpsiw(1,1,3,igroup),1,dpsiw(1,1,isnow2,igroup),1)
!!$           
!!$           call MPI_ISEND(dpsiw(1,1,isnow2,igroup),&
!!$                nvctr_par(jprocsr(3,jproc,igroup),igrpr(igroup)),mpidtypw,&
!!$                jprocsr(3,jproc,igroup),&
!!$                iproc+nproc+2*nproc*jproc,bigdft_mpi%mpi_comm,mpireq2(ncommsstep2),ierr)
!!$        end if
!!$        if (jprocsr(4,jproc,igroup) /= -1) then
!!$           ncommsstep2=ncommsstep2+1
!!$           if (iproc == itestproc) then
!!$              print '(5(1x,a,i8))','step',jproc+1,'group:',igrpr(igroup),&
!!$                   ': receiving',nvctr_par(iproc,igrpr(igroup)),&
!!$                   'elements from',jprocsr(4,jproc,igroup),'to',iproc
!!$           end if
!!$           call MPI_IRECV(dpsiw(1,1,irnow2,igroup),&
!!$                nvctr_par(iproc,igrpr(igroup)),mpidtypw,jprocsr(4,jproc,igroup),&
!!$                jprocsr(4,jproc,igroup)+nproc+2*nproc*jproc,bigdft_mpi%mpi_comm,mpireq2(ncommsstep2),ierr)
!!$           
!!$        end if
!!$     end do
!!$     if (jproc>1) isnow2=3-isnow2
!!$     
!!$     if (ncommsstep /=0) then
!!$        !verify that the messages have been passed
!!$        !print *,'waiting,iproc',iproc
!!$        call MPI_WAITALL(ncommsstep,mpireq,mpistat,ierr)
!!$        if (ierr /=0) print *,'step,ierr',jproc+1,iproc,ierr,mpistat,MPI_STATUSES_IGNORE
!!$        !print *,'done,iproc',iproc
!!$     end if
!!$     isnow=3-isnow
!!$     ncommsstep=0
!!$  end do
!!$  
!!$  !call MPI_BARRIER(bigdft_mpi%mpi_comm,ierr)
!!$  call mpiallred(eexctX,1,MPI_SUM,bigdft_mpi%mpi_comm,ierr)
!!$  
!!$  !the exact exchange energy is half the Hartree energy (which already has another half)
!!$  eexctX=-exctXfac*eexctX
!!$  
!!$  if (iproc == 0) write(*,'(a,1x,1pe18.11)')'Exact Exchange Energy:',eexctX
!!$  !close(100+iproc)
!!$  i_all=-product(shape(nvctr_par))*kind(nvctr_par)
!!$  deallocate(nvctr_par,stat=i_stat)
!!$  call memocc(i_stat,i_all,'nvctr_par',subname)
!!$  
!!$  i_all=-product(shape(iorbgr))*kind(iorbgr)
!!$  deallocate(iorbgr,stat=i_stat)
!!$  call memocc(i_stat,i_all,'iorbgr',subname)
!!$
!!$
!!$
!!$  i_all=-product(shape(rp_ij))*kind(rp_ij)
!!$  deallocate(rp_ij,stat=i_stat)
!!$  call memocc(i_stat,i_all,'rp_ij',subname)
!!$  
!!$  i_all=-product(shape(psiw))*kind(psiw)
!!$  deallocate(psiw,stat=i_stat)
!!$  call memocc(i_stat,i_all,'psiw',subname)
!!$
!!$  i_all=-product(shape(dpsiw))*kind(dpsiw)
!!$  deallocate(dpsiw,stat=i_stat)
!!$  call memocc(i_stat,i_all,'dpsiw',subname)
!!$
!!$
!!$  i_all=-product(shape(psir))*kind(psir)
!!$  deallocate(psir,stat=i_stat)
!!$  call memocc(i_stat,i_all,'psir',subname)
!!$
!!$  i_all=-product(shape(igrpr))*kind(igrpr)
!!$  deallocate(igrpr,stat=i_stat)
!!$  call memocc(i_stat,i_all,'igrpr',subname)
!!$
!!$  i_all=-product(shape(iprocpm1))*kind(iprocpm1)
!!$  deallocate(iprocpm1,stat=i_stat)
!!$  call memocc(i_stat,i_all,'iprocpm1',subname)
!!$
!!$  i_all=-product(shape(jprocsr))*kind(jprocsr)
!!$  deallocate(jprocsr,stat=i_stat)
!!$  call memocc(i_stat,i_all,'jprocsr',subname)
!!$
!!$  !call timing(iproc,'Exchangecorr  ','OF')
!!$
!!$END SUBROUTINE exact_exchange_potential_round_new
!!$
!!$subroutine OP2P_comm_simulation(iproc,nproc,op2p)
!!$  use module_base
!!$  use module_types
!!$  implicit none
!!$  integer, intent(in) :: iproc,nproc
!!$  type(op2p_descriptors), intent(in) :: op2p
!!$  !local variables
!!$  integer :: itestproc,nend,jproc,kproc
!!$  integer, dimension(:), allocatable :: ndatac
!!$  integer, dimension(:,:,:), allocatable :: ndatas
!!$
!!$  !test array for data sending
!!$  allocate(ndatas(2,0:nproc-1,op2p%ngroup+ndebug),stat=i_stat)
!!$  call memocc(i_stat,ndatas,'ndatas',subname)
!!$
!!$  !test array for data calculation
!!$  allocate(ndatac(op2p%ngroupp+ndebug),stat=i_stat)
!!$  call memocc(i_stat,ndatac,'ndatac',subname)
!!$
!!$  do igroup=1,op2p%ngroupp
!!$     !define the number of data to calculate in total
!!$     ndatac(igroup)=0
!!$     do kproc=0,nproc-1
!!$        ndatac(igroup)=ndatac(igroup)-op2p%nvctr_par(kproc,op2p%igrpr(igroup))
!!$     end do
!!$  end do
!!$
!!$  itestproc=-1 !-1=no debug verbosity
!!$  if (itestproc > -1) then
!!$     !simulation of communication
!!$     nend=(nproc-1)/2+1
!!$     ndatas=0
!!$     do jproc=0,nend
!!$        !sending receiving data
!!$        do igroup=1,op2p%ngroupp
!!$           if (op2p%jprocsr(1,jproc,igroup) /= -1) then
!!$              !send the fixed array to the processor which comes in the list
!!$              if (op2p%iprocpm1(1,1,igroup) == itestproc) then
!!$                 print *,'step',jproc+1,': sending',op2p%nvctr_par(op2p%jprocsr(1,jproc,igroup),op2p%igrpr(igroup)),&
!!$                      'elements from',iproc,'to',op2p%iprocpm1(1,1,igroup)
!!$              end if
!!$              ndatas(1,op2p%iprocpm1(1,1,igroup),op2p%igrpr(igroup))=ndatas(1,op2p%iprocpm1(1,1,igroup),op2p%igrpr(igroup))+&
!!$                   op2p%nvctr_par(op2p%jprocsr(1,jproc,igroup),op2p%igrpr(igroup))
!!$           end if
!!$           if (op2p%jprocsr(2,jproc,igroup) /= -1) then
!!$              if (iproc == itestproc) then
!!$                 print *,'step',jproc+1,': receiving',op2p%nvctr_par(op2p%jprocsr(2,jproc,igroup),op2p%igrpr(igroup)),&
!!$                      'elements from',op2p%iprocpm1(2,1,igroup),'to',iproc
!!$              end if
!!$              ndatas(1,iproc,op2p%igrpr(igroup))=ndatas(1,iproc,op2p%igrpr(igroup))-&
!!$                   op2p%nvctr_par(op2p%jprocsr(2,jproc,igroup),op2p%igrpr(igroup))
!!$           end if
!!$        end do
!!$
!!$        !calculation for orbitals to be performed
!!$        do igroup=1,op2p%ngroupp
!!$           if (jproc==0) then
!!$              ndatac(igroup)=ndatac(igroup)+op2p%nvctr_par(iproc,op2p%igrpr(igroup))
!!$           else
!!$              if (op2p%jprocsr(2,jproc-1,igroup) /=-1) then
!!$                 ndatac(igroup)=ndatac(igroup)+&
!!$                      op2p%nvctr_par(jprocsr(2,jproc-1,igroup),op2p%igrpr(igroup))  
!!$                 if (iproc == itestproc) then
!!$                    print '(5(1x,a,i8))','step',jproc+1,'group:',op2p%igrpr(igroup),&
!!$                         ':processing',op2p%nvctr_par(op2p%jprocsr(2,jproc-1,igroup),op2p%igrpr(igroup)),&
!!$                         'elements in',iproc,'from',op2p%jprocsr(2,jproc-1,igroup)
!!$                 end if
!!$              end if
!!$           end if
!!$        end do
!!$
!!$        !copy the results which have been received
!!$        if (ncommsstep2 > 0) then
!!$           do igroup=1,op2p%ngroupp
!!$              if (op2p%jprocsr(4,jproc-1,igroup) /= -1) then
!!$                 if (iproc == itestproc) then
!!$                    print '(5(1x,a,i8))','step',jproc+1,'group:',op2p%igrpr(igroup),&
!!$                         ':copying',op2p%nvctr_par(op2p%jprocsr(4,jproc-1,igroup),op2p%igrpr(igroup)),&
!!$                         'elements from',op2p%jprocsr(4,jproc-1,igroup),'in',iproc
!!$                 end if
!!$                 ndatac(igroup)=ndatac(igroup)+&
!!$                      op2p%nvctr_par(op2p%jprocsr(4,jproc-1,igroup),op2p%igrpr(igroup)) 
!!$              end if
!!$           end do
!!$        end if
!!$
!!$        !send-receive of the results
!!$        ncommsstep2=0
!!$        do igroup=1,op2p%ngroupp
!!$           if (op2p%jprocsr(3,jproc,igroup) /= -1) then
!!$              ncommsstep2=ncommsstep2+1
!!$              if (op2p%jprocsr(3,jproc,igroup) == itestproc) then
!!$                 print '(5(1x,a,i8))','step',jproc+1,'group:',op2p%igrpr(igroup),&
!!$                      ': sending',op2p%nvctr_par(op2p%jprocsr(3,jproc,igroup),op2p%igrpr(igroup)),&
!!$                      'elements from',iproc,'to',op2p%jprocsr(3,jproc,igroup)
!!$              end if
!!$              ndatas(2,op2p%jprocsr(3,jproc,igroup),op2p%igrpr(igroup))=ndatas(2,op2p%jprocsr(3,jproc,igroup),op2p%igrpr(igroup))+&
!!$                   op2p%nvctr_par(op2p%jprocsr(3,jproc,igroup),op2p%igrpr(igroup))
!!$           end if
!!$           if (op2p%jprocsr(4,jproc,igroup) /= -1) then
!!$              ncommsstep2=ncommsstep2+1
!!$              if (iproc == itestproc) then
!!$                 print '(5(1x,a,i8))','step',jproc+1,'group:',op2p%igrpr(igroup),&
!!$                      ': receiving',op2p%nvctr_par(iproc,op2p%igrpr(igroup)),&
!!$                      'elements from',op2p%jprocsr(4,jproc,igroup),'to',iproc
!!$              end if
!!$              ndatas(2,iproc,op2p%igrpr(igroup))=ndatas(2,iproc,op2p%igrpr(igroup))-&
!!$                   op2p%nvctr_par(iproc,op2p%igrpr(igroup))
!!$           end if
!!$        end do
!!$
!!$     end do
!!$
!!$     call mpiallred(ndatas(1,0,1),2*nproc*op2p%ngroup,MPI_SUM,bigdft_mpi%mpi_comm,ierr)
!!$     !if(iproc ==0)print *,'iproc,datas',iproc,ndatas
!!$
!!$     do igroup=1,op2p%ngroupp
!!$        if (ndatac(igroup) /=0) then
!!$           write(*,*)'ERROR: OP2P communication simulation failed: processor',iproc,&
!!$                ' has calculated',ndatac(igroup),' data more than needed'
!!$           stop
!!$        end if
!!$        if (ndatas(1,iproc,op2p%igrpr(igroup)) /=0 .or. ndatas(2,iproc,op2p%igrpr(igroup)) /=0) then
!!$           write(*,*)'ERROR: OP2P communication simulation failed: processor',iproc,&
!!$                ' has not a zero balance of send-receive calls',ndatas(1:2,iproc,op2p%igrpr(igroup))
!!$           stop
!!$        end if
!!$     end do
!!$  end if
!!$
!!$  i_all=-product(shape(ndatas))*kind(ndatas)
!!$  deallocate(ndatas,stat=i_stat)
!!$  call memocc(i_stat,i_all,'ndatas',subname)
!!$
!!$  i_all=-product(shape(ndatac))*kind(ndatac)
!!$  deallocate(ndatac,stat=i_stat)
!!$  call memocc(i_stat,i_all,'ndatac',subname)
!!$
!!$
!!$END SUBROUTINE OP2P_comm_simulation
