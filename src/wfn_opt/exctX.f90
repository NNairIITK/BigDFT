!> @file
!!   Routines to calculate the exact exchange potential
!! @author
!!    Copyright (C) 2010-2013 BigDFT group 
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 


!> Calculate the exact exchange potential
subroutine exact_exchange_potential(iproc,nproc,geocode,xc,nspin,lr,orbs,n3parr,n3p,&
     hxh,hyh,hzh,pkernel,psi,psir,eexctX)

  use module_base
  use module_types
  use Poisson_Solver, except_dp => dp, except_gp => gp
  use module_xc
  use yaml_output
  use locreg_operations
  implicit none
  character(len=1), intent(in) :: geocode !< @copydoc poisson_solver::doc::geocode             !< Determine Boundary conditions
  integer, intent(in) :: iproc,nproc                  !< MPI information
  integer, intent(in) :: n3p,nspin                    !< spin and ...
  real(gp), intent(in) :: hxh,hyh,hzh                 !< hgrid
  type(xc_info), intent(in) :: xc
  type(locreg_descriptors), intent(in) :: lr          !< Local region descriptors
  type(orbitals_data), intent(in) :: orbs             !< Orbitals
  integer, dimension(0:nproc-1), intent(in) :: n3parr
  real(wp), dimension(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f,orbs%nspinor,orbs%norbp), intent(in) :: psi
  type(coulomb_operator), intent(inout) :: pkernel
  real(gp), intent(out) :: eexctX
  real(wp), dimension(max(lr%d%n1i*lr%d%n2i*lr%d%n3i*orbs%norbp,n3parr(0)*orbs%norb)), intent(out) :: psir
  !local variables
  character(len=*), parameter :: subname='exact_exchange_potential'
  integer :: ierr,ispinor,ispsiw,ispin,norb,ncall
  integer :: i1,i2,i3p,iorb,iorbs,jorb,jorbs,ispsir,ind3,ind2,ind1i,ind1j,jproc,igran,ngran
  real(gp) :: ehart,hfac,exctXfac,sign,sfac,hfaci,hfacj
  type(workarr_sumrho) :: w
  integer, dimension(:,:), allocatable :: ncommarr
  real(wp), dimension(:), allocatable :: psiw
  real(wp), dimension(:,:,:,:), allocatable :: rp_ij

  !call timing(iproc,'Exchangecorr  ','ON')

  exctXfac = xc_exctXfac(xc)

  eexctX=0.0_gp

  call initialize_work_arrays_sumrho(1,[lr],.true.,w)
  
  !save the value of the previous offset of the kernel
  !kerneloff=pkernel(1)
  !put to szero the offset to subtract the energy of the momentum
  !pkernel(1)=0.0_dp

  !the granularity of the calculation is set by ngran
  !for the moment it is irrelevant but if the poisson solver is modified
  !we may increase this value
  ngran=1

  !partial densities with a given granularity
  rp_ij = f_malloc((/ lr%d%n1i, lr%d%n2i, n3p, ngran /),id='rp_ij')
  psiw = f_malloc(max(max(lr%d%n1i*lr%d%n2i*lr%d%n3i*orbs%norbp, n3parr(0)*orbs%norb), 1),id='psiw')

  if (geocode == 'F') then
     !call to_zero(lr%d%n1i*lr%d%n2i*lr%d%n3i*orbs%norbp,psiw)
     call f_zero(psiw)
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
        call vcopy(n3parr(jproc),psiw(ispsiw),1,psir(ispsir),1)
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

     ncommarr = f_malloc((/ 0.to.nproc-1, 1.to.4 /),id='ncommarr')

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
     call vcopy(lr%d%n1i*lr%d%n2i*n3p*orbs%norb,psir(1),1,psiw(1),1)
  end if

  !this is the array of the actions of the X potential on psi
  call f_zero(lr%d%n1i*lr%d%n2i*n3p*orbs%norb,psir(1))

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
     !call vcopy(lr%d%n1i*lr%d%n2i*n3p*orbs%norb,psir,1,psirt,1)
     !recommunicate the values in the psir array
     call MPI_ALLTOALLV(psir,ncommarr(0,3),ncommarr(0,4),mpidtypw, &
          psiw,ncommarr(0,1),ncommarr(0,2),mpidtypw,bigdft_mpi%mpi_comm,ierr)
     !redress the potential
     ispsiw=1
     do iorb=1,orbs%norbp
        ispsir=1+(iorb-1)*n3parr(0)
        do jproc=0,nproc-1
           call vcopy(n3parr(jproc),psiw(ispsir),1,psir(ispsiw),1)
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

  call f_free(rp_ij)
  call f_free(psiw)


  if (nproc > 1) then
     call f_free(ncommarr)
  end if

  !call timing(iproc,'Exchangecorr  ','OF')

END SUBROUTINE exact_exchange_potential


subroutine prepare_psirocc(iproc,nproc,lr,orbsocc,n3p,n3parr,psiocc,psirocc)
  use module_base
  use module_types
  use locreg_operations
  implicit none
  integer, intent(in) :: iproc,nproc,n3p
  type(locreg_descriptors), intent(in) :: lr
  type(orbitals_data), intent(in) :: orbsocc
  integer, dimension(0:nproc-1), intent(in) :: n3parr
  real(wp), dimension(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f,orbsocc%nspinor,orbsocc%norbp), intent(in) :: psiocc
  real(wp), dimension(max(lr%d%n1i*lr%d%n2i*lr%d%n3i*orbsocc%norbp,n3parr(0)*orbsocc%norb)), intent(out) :: psirocc
  !local variables
  character(len=*), parameter :: subname='prepare_psirocc'
  integer :: ierr,ispinor,ispsiw
  integer :: iorb,jorb,ispsir,jproc
  type(workarr_sumrho) :: w
  integer, dimension(:,:), allocatable :: ncommocc
  real(wp), dimension(:), allocatable :: psiwocc

  call initialize_work_arrays_sumrho(1,[lr],.true.,w)

  psiwocc = f_malloc0(max(max(lr%d%n1i*lr%d%n2i*lr%d%n3i*orbsocc%norbp, n3parr(0)*orbsocc%norb), 1),id='psiwocc')

  !call to_zero(max(max(lr%d%n1i*lr%d%n2i*lr%d%n3i*orbsocc%norbp,n3parr(0)*orbsocc%norb),1),psiwocc)

  if (lr%geocode == 'F') then
     !call f_zero(lr%d%n1i*lr%d%n2i*lr%d%n3i*orbsocc%norbp,psirocc)
     call f_zero(psirocc)
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
        call vcopy(n3parr(jproc),psirocc(ispsiw),1,psiwocc(ispsir),1)
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
     ncommocc = f_malloc((/ 0.to.nproc-1, 1.to.4 /),id='ncommocc')

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
     call vcopy(lr%d%n1i*lr%d%n2i*n3p*orbsocc%norb,psiwocc(1),1,psirocc(1),1)
  end if
  call f_free(psiwocc)

  if (nproc > 1) then
     call f_free(ncommocc)
  end if

END SUBROUTINE prepare_psirocc


!> Calculate the exact exchange potential only on virtual orbitals
!! by knowing the occupied orbitals and their distribution
!! both sets of orbitals are to be 
subroutine exact_exchange_potential_virt(iproc,nproc,geocode,nspin,lr,orbsocc,orbsvirt,n3parr,n3p,&
     hxh,hyh,hzh,pkernel,psirocc,psivirt,psirvirt)
  use module_base
  use module_types
  use Poisson_Solver, except_dp => dp, except_gp => gp
  use yaml_output
  use locreg_operations
  implicit none
  character(len=1), intent(in) :: geocode !< @copydoc poisson_solver::doc::geocode
  integer, intent(in) :: iproc,nproc,n3p,nspin
  real(gp), intent(in) :: hxh,hyh,hzh
  type(locreg_descriptors), intent(in) :: lr
  type(orbitals_data), intent(in) :: orbsocc,orbsvirt
  integer, dimension(0:nproc-1), intent(in) :: n3parr
  real(wp), dimension(max(lr%d%n1i*lr%d%n2i*lr%d%n3i*orbsocc%norbp,n3parr(0)*orbsocc%norb)), intent(in) :: psirocc
  real(wp), dimension(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f,orbsvirt%nspinor,orbsvirt%norbp), intent(in) :: psivirt
  type(coulomb_operator), intent(inout) :: pkernel
  real(wp), dimension(max(lr%d%n1i*lr%d%n2i*lr%d%n3i*orbsvirt%norbp,n3parr(0)*orbsvirt%norb)), intent(out) :: psirvirt
  !local variables
  character(len=*), parameter :: subname='exact_exchange_potential_virt'
  integer :: ierr,ispinor,ispsiw,ispin,norbocc,norbvirt
  integer :: i1,i2,i3p,iorb,iorbs,jorb,jorbs,ispsir,ind3,ind2,ind1i,ind1j,jproc,igran,ngran
  real(gp) :: ehart,hfac,sign,sfac,hfaci
  type(workarr_sumrho) :: w
  integer, dimension(:,:), allocatable :: ncommvirt
  real(wp), dimension(:), allocatable :: psiwvirt
  real(wp), dimension(:,:,:,:), allocatable :: rp_ij

  !call timing(iproc,'Exchangecorr  ','ON')

  call initialize_work_arrays_sumrho(1,[lr],.true.,w)
  
  !the granularity of the calculation is set by ngran
  !for the moment it is irrelevant but if the poisson solver is modified
  !we may increase this value
  ngran=1

  !partial densities with a given granularity
  rp_ij = f_malloc((/ lr%d%n1i, lr%d%n2i, n3p, ngran /),id='rp_ij')
  psiwvirt = f_malloc(max(max(lr%d%n1i*lr%d%n2i*lr%d%n3i*orbsvirt%norbp, n3parr(0)*orbsvirt%norb), 1),id='psiwvirt')

  if (geocode == 'F') then
     call f_zero(psiwvirt)
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
        call vcopy(n3parr(jproc),psiwvirt(ispsiw),1,psirvirt(ispsir),1)
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

     ncommvirt = f_malloc((/ 0.to.nproc-1, 1.to.4 /),id='ncommvirt')

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
     call vcopy(lr%d%n1i*lr%d%n2i*n3p*orbsvirt%norb,psirvirt(1),1,psiwvirt(1),1)
  end if

  !this is the array of the actions of the X potential on psi
  call f_zero(psirvirt)

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
     !call vcopy(lr%d%n1i*lr%d%n2i*n3p*orbs%norb,psir,1,psirt,1)
     !recommunicate the values in the psir array
     call MPI_ALLTOALLV(psirvirt,ncommvirt(0,3),ncommvirt(0,4),mpidtypw, &
          psiwvirt,ncommvirt(0,1),ncommvirt(0,2),mpidtypw,bigdft_mpi%mpi_comm,ierr)
     !redress the potential
     ispsiw=1
     do iorb=1,orbsvirt%norbp
        ispsir=1+(iorb-1)*n3parr(0)
        do jproc=0,nproc-1
           call vcopy(n3parr(jproc),psiwvirt(ispsir),1,psirvirt(ispsiw),1)
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

  call f_free(rp_ij)
  call f_free(psiwvirt)


  if (nproc > 1) then
     call f_free(ncommvirt)
  end if

END SUBROUTINE exact_exchange_potential_virt


!> Calculate the exact exchange potential on occupied orbitals
!! within the symmetric round-robin scheme
!! the psi is already given in the real-space form
subroutine exact_exchange_potential_round(iproc,nproc,xc,nspin,lr,orbs,&
     hxh,hyh,hzh,pkernel,psi,dpsir,eexctX)
  use module_base
  use module_types
  use Poisson_Solver, except_dp => dp, except_gp => gp
  use module_xc
  use yaml_output
  use locreg_operations
  implicit none
  integer, intent(in) :: iproc,nproc,nspin
  real(gp), intent(in) :: hxh,hyh,hzh
  type(xc_info), intent(in) :: xc
  type(locreg_descriptors), intent(in) :: lr
  type(orbitals_data), intent(in) :: orbs
  real(wp), dimension(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f,orbs%nspinor,orbs%norbp), intent(in) :: psi
  type(coulomb_operator), intent(inout) :: pkernel
  real(gp), intent(out) :: eexctX
  real(wp), dimension(lr%d%n1i*lr%d%n2i*lr%d%n3i,orbs%norbp), intent(out) :: dpsir
  !local variables
  character(len=*), parameter :: subname='exact_exchange_potential_round'
  logical :: doit, use_mpi_get,new_mpi_pattern
  logical :: get_post,get_start,acc_post,acc_start
  integer :: ierr,ncommsstep,ncommsstep2,isnow,irnow,isnow2,irnow2,jsorb,kproc,norbp
  integer :: i,iorb,jorb,jproc,igroup,ngroup,ngroupp,nend,isorb,iorbs,jorbs,ii,iproc_totake,iproc_toput
  integer :: icount,nprocgr,iprocgrs,iprocgrr,itestproc,norbi,norbj,ncalltot,icountmax,iprocref,ncalls
  real(gp) :: ehart,hfac,exctXfac,sfac,hfaci,hfacj,hfac2
  integer, dimension(4) :: mpireq,mpireq2
  integer, dimension(MPI_STATUS_SIZE,4) :: mpistat,mpistat2
  type(workarr_sumrho) :: w
  integer, dimension(:), allocatable :: igrpr,ndatac
  integer, dimension(:,:), allocatable :: nvctr_par,igrprarr
  integer, dimension(:,:,:), allocatable :: jprocsr,iprocpm1,ndatas,iorbgr
  real(wp), dimension(:), allocatable :: rp_ij
  real(wp), dimension(:,:), allocatable :: psir
  real(wp), dimension(:,:,:,:), allocatable :: psiw,dpsiw
  integer :: win, win2, win3, size_of_double,win4,base_group
  integer, dimension(2) :: grp_acc_start,grp_acc_post
  integer :: grp_put,grp_get

  !call timing(iproc,'Exchangecorr  ','ON')
  use_mpi_get = .false.
  new_mpi_pattern=.true.

  exctXfac = xc_exctXfac(xc)

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
  nvctr_par = f_malloc((/ 0.to.nproc-1, 1.to.ngroup /),id='nvctr_par')

  iorbgr = f_malloc((/ 1.to.2, 0.to.nproc-1, 1.to.ngroup /),id='iorbgr')

  !test array for data sending
  ndatas = f_malloc((/ 1.to.2, 0.to.nproc-1, 1.to.ngroup /),id='ndatas')


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
  igrpr = f_malloc(ngroupp,id='igrpr')
  igrprarr= f_malloc0((/ 1.to.ngroup, 0.to.nproc-1 /),id='igrprarr')
  iprocpm1 = f_malloc((/ 1.to.2, 0.to.nproc-1, 1.to.ngroupp /),id='iprocpm1')

  !test array for data calculation
  ndatac = f_malloc(ngroupp,id='ndatac')


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
           igrprarr(igroup,kproc)=icount
        end if
     end do
     if (icount > icountmax) then
        iprocref=kproc
        icountmax=icount
     end if
  end do

  !calculate the list of send-receive operations which have to be performed per group
  !allocate it at the maximum size needed
  jprocsr = f_malloc((/ 1.to.4, 0.to.nproc/2+1, 1.to.ngroupp /),id='jprocsr')
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
           if(new_mpi_pattern) then
               jprocsr(1,kproc,igroup)= iprocpm1(1,kproc+1,igroup)
               jprocsr(2,kproc,igroup)= iprocpm1(2,kproc+1,igroup)
           else
               jprocsr(1,kproc,igroup)=iprocpm1(2,kproc,igroup)
               jprocsr(2,kproc,igroup)=iprocpm1(2,kproc+1,igroup)
           end if 
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
           if(new_mpi_pattern) then
               jprocsr(1,kproc,igroup)= iprocpm1(1,kproc+1,igroup)
               jprocsr(2,kproc,igroup)= iprocpm1(2,kproc+1,igroup)
           else
               jprocsr(1,kproc,igroup)=iprocpm1(2,kproc,igroup)
               jprocsr(2,kproc,igroup)=iprocpm1(2,kproc+1,igroup)
           end if
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

  itestproc=-1! -1=no debug verbosity
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

     if (nproc > 1) call mpiallred(ndatas,MPI_SUM,comm=bigdft_mpi%mpi_comm)
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

  call initialize_work_arrays_sumrho(1,[lr],.true.,w)
  psir = f_malloc0((/ lr%d%n1i*lr%d%n2i*lr%d%n3i, orbs%norbp /),id='psir')

  !call to_zero(lr%d%n1i*lr%d%n2i*lr%d%n3i*orbs%norbp,psir(1,1))

  !uncompress the wavefunction in the real grid
  do iorb=1,orbs%norbp
     !here ispinor is equal to one
     call daub_to_isf(lr,w,psi(1,1,iorb),psir(1,iorb))
  end do

  call deallocate_work_arrays_sumrho(w)

  psiw = f_malloc0((/ 1.to.lr%d%n1i*lr%d%n2i*lr%d%n3i, 1.to.maxval(orbs%norb_par(:,0)), 1.to.2, 1.to.ngroupp /),id='psiw')


  if (nproc>1) then
     call mpi_type_size(mpi_double_precision, size_of_double, ierr)
  else
     size_of_double = 8
  end if



  dpsiw = f_malloc0((/ 1.to.lr%d%n1i*lr%d%n2i*lr%d%n3i, 1.to.maxval(orbs%norb_par(:,0)), 1.to.3, 1.to.ngroupp /),id='dpsiw')
  !partial densities and potentials
  rp_ij = f_malloc(lr%d%n1i*lr%d%n2i*lr%d%n3i,id='rp_ij')

  if(use_mpi_get) then 
     if (.not. new_mpi_pattern) then
        win = mpiwindow(lr%d%n1i*lr%d%n2i*lr%d%n3i*maxval(orbs%norb_par(:,0))*2*ngroupp,&
             psiw(1,1,1,1),bigdft_mpi%mpi_comm)
        win2 = mpiwindow(lr%d%n1i*lr%d%n2i*lr%d%n3i*maxval(orbs%norb_par(:,0))*3*ngroupp,dpsiw(1,1,1,1), bigdft_mpi%mpi_comm)
     end if
     ! To use MPI_get, we need a window covering psir, which will be used only once
     win3 = mpiwindow(lr%d%n1i*lr%d%n2i*lr%d%n3i*orbs%norbp,psir(1,1), bigdft_mpi%mpi_comm)
     if (new_mpi_pattern) then
        win=MPI_WIN_NULL
        win2=MPI_WIN_NULL
        win4 = mpiwindow(lr%d%n1i*lr%d%n2i*lr%d%n3i*orbs%norbp,dpsir(1,1), bigdft_mpi%mpi_comm)
        base_group=mpigroup(bigdft_mpi%mpi_comm)
     end if
  end if

  !this is the array of the actions of the X potential on psi
  !ii=lr%d%n1i*lr%d%n2i*lr%d%n3i*maxval(orbs%norb_par(:,0))*2*ngroupp
  !call to_zero(lr%d%n1i*lr%d%n2i*lr%d%n3i*maxval(orbs%norb_par,1)*2*ngroupp,psiw(1,1,1,1))
  !call to_zero(ii,psiw(1,1,1,1))
  !ii=lr%d%n1i*lr%d%n2i*lr%d%n3i*maxval(orbs%norb_par(:,0))*3*ngroupp
  !call to_zero(lr%d%n1i*lr%d%n2i*lr%d%n3i*maxval(orbs%norb_par,1)*3*ngroupp,dpsiw(1,1,1,1))
  !call to_zero(ii,dpsiw(1,1,1,1))

  call f_zero(dpsir)
  ncalls=0
  !real communication
  isnow=1
  isnow2=1
  nend=(nproc-1)/2+1
  ncommsstep2=0

  get_post=.false.
  get_start=.false.
  acc_post=.false.
  acc_start=.false.
  grp_get=mpigroup_null()
  grp_put=mpigroup_null()
  grp_acc_start=mpigroup_null()
  grp_acc_post=mpigroup_null()
  do jproc=0,nend
     irnow=3-isnow
     ncommsstep=0
     !sending receiving data
     do igroup=1,ngroupp

        if (jprocsr(1,jproc,igroup) /= -1) then
           ncommsstep=ncommsstep+1
           if (iprocpm1(1,1,igroup) == itestproc) then
              print *,'step',jproc+1,': sending',nvctr_par(iproc,igrpr(igroup)),&
                   'elements from',iproc,'to',jprocsr(1,jproc,igroup)
           end if
           if( .not. use_mpi_get) then    
             if (.not. new_mpi_pattern) then
               if (jproc == 0) then
                 call MPI_ISEND(psir(1,iorbgr(2,iproc,igrpr(igroup))),nvctr_par(jprocsr(1,jproc,igroup),igrpr(igroup)),&
                      mpidtypw,iprocpm1(1,1,igroup),&
                      iproc+2*nproc*jproc,bigdft_mpi%mpi_comm,mpireq(ncommsstep),ierr)
               else
                 call MPI_ISEND(psiw(1,1,isnow,igroup),nvctr_par(jprocsr(1,jproc,igroup),igrpr(igroup)),&
                      mpidtypw,iprocpm1(1,1,igroup),&
                      iproc+2*nproc*jproc,bigdft_mpi%mpi_comm,mpireq(ncommsstep),ierr)
               end if
             else
             call MPI_ISEND(psir(1,iorbgr(2,iproc,igrpr(igroup))),nvctr_par(iproc,igrpr(igroup)),&
                      mpidtypw,jprocsr(1,jproc,igroup),&
                      iproc,bigdft_mpi%mpi_comm,mpireq(ncommsstep),ierr)
             end if
           else if (new_mpi_pattern) then
              !create exposure epoch for the window win3
              !create the passive group, that should match the creation of an active group on a remote proc
              if (ngroupp==2 .and. jprocsr(1,jproc,ngroupp) /= -1) then
                 if (igroup==1) grp_put=p2p_group(base_group,p1=iproc,&
                                        p2=jprocsr(1,jproc,igroup),p3=jprocsr(1,jproc,2))

              else

                 grp_put=p2p_group(base_group,p1=iproc,p2=jprocsr(1,jproc,igroup))
              end if
              if (igroup==1) then
                 call mpiwinpost(grp_put,win3,MPI_MODE_NOCHECK+MPI_MODE_NOSTORE+MPI_MODE_NOPUT)
                 get_post=.true.
              end if
           end if
        end if

        if (jprocsr(2,jproc,igroup) == -1 .and. igroup==1 .and. new_mpi_pattern .and. use_mpi_get) then
           if (ngroupp==2 .and. jprocsr(2,jproc,ngroupp) /= -1) then
              grp_get=p2p_group(base_group,p1=iproc,p2=jprocsr(2,jproc,2))
              call mpiwinstart(grp_get,win3,MPI_MODE_NOCHECK)
              get_start=.true.
           end if
        end if
        
        if (jprocsr(2,jproc,igroup) /= -1) then
           ncommsstep=ncommsstep+1
           if (iproc == itestproc) then
              print *,'step',jproc+1,': receiving',nvctr_par(jprocsr(2,jproc,igroup),igrpr(igroup)),&
                   'elements from',jprocsr(2,jproc,igroup),'to',iproc
           end if

           if(use_mpi_get) then 
              if (new_mpi_pattern) then
                 iproc_totake=jprocsr(2,jproc,igroup) !this should always be the same
                 !create the active group, that should match the creation of a passive group on a remote proc
                 if (ngroupp==2 .and. jprocsr(2,jproc,ngroupp) /= -1) then
                    if (igroup==1) grp_get=p2p_group(base_group,p1=iproc,p2=jprocsr(2,jproc,igroup),p3=jprocsr(2,jproc,2))
                 else
                    if (.not. get_start) grp_get=p2p_group(base_group,p1=iproc,p2=jprocsr(2,jproc,igroup))
                 end if
                 if (.not. get_start) then
                    call mpiwinstart(grp_get,win3,MPI_MODE_NOCHECK)
                    get_start=.true.
                 end if
                 call mpiget(origin=psiw(1,1,irnow,igroup), &
                      count=nvctr_par(iproc_totake,igrpr(igroup)),&
                      target_rank=iproc_totake,&
                      target_disp=int((iorbgr(2,iproc_totake,igrprarr(igrpr(igroup), iproc_totake))-1)&
                      *(lr%d%n1i*lr%d%n2i*lr%d%n3i), kind=mpi_address_kind),window=win3)
              else
                 if (jproc == 0) then
                    call mpiget(origin=psiw(1,1,irnow,igroup), &
                         count=nvctr_par(jprocsr(2,jproc,igroup),igrpr(igroup)),&
                         target_rank=iprocpm1(2,1,igroup),&
                         target_disp=int((iorbgr(2,iprocpm1(2,1,igroup),igrprarr(igrpr(igroup), iprocpm1(2,1,igroup)))-1)&
                         *(lr%d%n1i*lr%d%n2i*lr%d%n3i), kind=mpi_address_kind),window=win3)
                 else             
                    call mpiget(origin=psiw(1,1,irnow,igroup), &
                         count=nvctr_par(jprocsr(2,jproc,igroup),igrpr(igroup)),&
                         target_rank=iprocpm1(2,1,igroup),&
                         target_disp=int((igrprarr(igrpr(igroup), iprocpm1(2,1,igroup))-1)*&
                         (lr%d%n1i*lr%d%n2i*lr%d%n3i*maxval(orbs%norb_par(:,0)*2))&
                         + (isnow-1)*(lr%d%n1i*lr%d%n2i*lr%d%n3i*maxval(orbs%norb_par(:,0))), kind=mpi_address_kind),&
                         window=win)
                 end if
              end if
           else
             if (.not. new_mpi_pattern) then
              call MPI_IRECV(psiw(1,1,irnow,igroup),nvctr_par(jprocsr(2,jproc,igroup),igrpr(igroup)),&
                   mpidtypw,iprocpm1(2,1,igroup),&
                   iprocpm1(2,1,igroup)+2*nproc*jproc,bigdft_mpi%mpi_comm,mpireq(ncommsstep),ierr)
             else 
              call MPI_IRECV(psiw(1,1,irnow,igroup),nvctr_par(jprocsr(2,jproc,igroup),igrpr(igroup)),&
                   mpidtypw,jprocsr(2,jproc,igroup),&
                   jprocsr(2,jproc,igroup),bigdft_mpi%mpi_comm,mpireq(ncommsstep),ierr)

             end if
           end if
        end if
     end do

     do igroup=1,ngroupp
        if (jproc /= 0 .and. jprocsr(3,jproc,igroup) /= -1) then
           !put to zero the sending element
           ii=lr%d%n1i*lr%d%n2i*lr%d%n3i*maxval(orbs%norb_par(:,0))
           !call to_zero(lr%d%n1i*lr%d%n2i*lr%d%n3i*maxval(orbs%norb_par,1),dpsiw(1,1,3,igroup))
           call f_zero(ii,dpsiw(1,1,3,igroup))
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
                       !$omp parallel do default(shared) private(i)
                       do i=1,lr%d%n1i*lr%d%n2i*lr%d%n3i
                          rp_ij(i)=hfac*psir(i,iorb)*psir(i,jorb)
                       end do
                       !$omp end parallel do
                    else
                       !$omp parallel do default(shared) private(i)
                       do i=1,lr%d%n1i*lr%d%n2i*lr%d%n3i
                          rp_ij(i)=hfac*psir(i,iorb)*psiw(i,jorb-jorbs+1,isnow,igroup)
                       end do
                       !$omp end parallel do
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
                       !$omp parallel do default(shared) private(i)
                       do i=1,lr%d%n1i*lr%d%n2i*lr%d%n3i
                          dpsir(i,iorb)=dpsir(i,iorb)+&
                               hfaci*rp_ij(i)*psir(i,jorb)
                       end do
                       !$omp end parallel do
                       if (jorb+jsorb /= iorb+isorb) then
                          !$omp parallel do default(shared) private(i)
                          do i=1,lr%d%n1i*lr%d%n2i*lr%d%n3i
                             dpsir(i,jorb)=dpsir(i,jorb)+&
                                  hfacj*rp_ij(i)*psir(i,iorb)
                          end do
                          !$omp end parallel do
                          !write(100+iproc,*)jorb+jsorb,iorb+isorb,igrpr(igroup) 
                       end if
                    else
                       !this part is summed on the win4 in the new version,
                       !to be controlled if it conflicts with the mpi_accumulate
                       !$omp parallel do default(shared) private(i)
                       do i=1,lr%d%n1i*lr%d%n2i*lr%d%n3i
                          dpsir(i,iorb)=dpsir(i,iorb)+&
                               hfaci*rp_ij(i)*psiw(i,jorb-jorbs+1,isnow,igroup)
                       end do
                       !$omp end parallel do
                    end if
                    !write(100+iproc,*)iorb+isorb,jorb+jsorb,igrpr(igroup)
                 end if

                 !fill the set of the vector to be sent to the other processes
                 !in the first step the results are self-contained
                 if (jproc /= 0 .and. jprocsr(3,jproc,igroup) /= -1) then
                    !write(100+iproc,*)jorb+jsorb,iorb+isorb,igrpr(igroup)
                    !$omp parallel do default(shared) private(i)
                    do i=1,lr%d%n1i*lr%d%n2i*lr%d%n3i
                       dpsiw(i,jorb-jorbs+1,3,igroup)=dpsiw(i,jorb-jorbs+1,3,igroup)+&
                            hfacj*rp_ij(i)*psir(i,iorb)
                    end do
                    !$omp end parallel do
                 end if
              end do
           end do
        end if
     end do
     if(use_mpi_get) then
        if (new_mpi_pattern) then
           !here at the first passage the group has not yet been created
           !if (jproc/=0) then 
           !here assert=MPI_MODE_NOPRECEDE for jproc==0
           !call mpi_fence(win4,rma_grp=grp_active(isnow))
           !then free the group, as it will be recreated by the mpi_accumulate
           if (acc_start) then
              call mpiwincomplete(win4)
              call mpigroup_free(grp_acc_start(3-isnow2)) !as the irnow2 is toggled after
              acc_start=.false.
           end if
           if (acc_post) then
              call mpiwinwait(win4)
              call mpigroup_free(grp_acc_post(isnow2))
              acc_post=.false.
           end if
        else
           call mpi_fence(win2)
        end if
     end if
     if (ncommsstep2 > 0) then
        !verify that the messages have been passed
        if(.not. use_mpi_get) call MPI_WAITALL(ncommsstep2,mpireq2,mpistat2,ierr)
        !copy the results which have been received (the messages sending are after)
        !this part is already done by the mpi_accumulate
        if (.not. (use_mpi_get .and. new_mpi_pattern)) then
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
     end if

     ncommsstep2=0
     !meanwhile, we can receive the result from the processor which has the psi 

     irnow2=3-isnow2
     do igroup=1,ngroupp
        if (jprocsr(3,jproc,igroup) /= -1) then
           ncommsstep2=ncommsstep2+1
           if (jprocsr(3,jproc,igroup) == itestproc) then
              print *,'step',jproc+1,'group:',igrpr(igroup),&
                   ': accum',nvctr_par(jprocsr(3,jproc,igroup),igrpr(igroup)),&
                   'elements from',iproc,'to',jprocsr(3,jproc,igroup)
           end if
           call vcopy(nvctr_par(jprocsr(3,jproc,igroup),igrpr(igroup)),&
                dpsiw(1,1,3,igroup),1,dpsiw(1,1,isnow2,igroup),1)
           if(.not. use_mpi_get) then
              call MPI_ISEND(dpsiw(1,1,isnow2,igroup),&
                   nvctr_par(jprocsr(3,jproc,igroup),igrpr(igroup)),mpidtypw,&
                   jprocsr(3,jproc,igroup),&
                   iproc+nproc+2*nproc*jproc,bigdft_mpi%mpi_comm,mpireq2(ncommsstep2),ierr)
           else if (new_mpi_pattern) then
              !version with accumulate
              !print *,'XXXXXXXXXXXXXXXXhere',jprocsr(3,jproc,igroup),jproc,igroup,
              if (.not. acc_start) then
                 call mpiwinstart(grp_acc_start(isnow2),win4,MPI_MODE_NOCHECK)
                 acc_start=.true.
              end if
              iproc_toput=jprocsr(3,jproc,igroup)
              call mpiaccumulate(origin=dpsiw(1,1,isnow2,igroup),&
                   count=nvctr_par(iproc_toput,igrpr(igroup)),& !this one has to be changed for the version with put
                   target_rank=iproc_toput,&
                   target_disp=int((iorbgr(2,iproc_toput,igrpr(igroup))-1)*lr%d%n1i*lr%d%n2i*lr%d%n3i, kind=mpi_address_kind),&
                   op=MPI_SUM,window=win4)
           end if
        end if
     end do

     do igroup=1,ngroupp
        if (jprocsr(4,jproc,igroup) /= -1) then
           ncommsstep2=ncommsstep2+1
           if(use_mpi_get .and. .not. new_mpi_pattern) then
              call mpiget(origin=dpsiw(1,1,irnow2,igroup),&
                   count=nvctr_par(iproc,igrpr(igroup)),&
                   target_rank=jprocsr(4,jproc,igroup),&
                   target_disp=int((igrprarr(igrpr(igroup), jprocsr(4,jproc,igroup))-1)*&
                   (lr%d%n1i*lr%d%n2i*lr%d%n3i*maxval(orbs%norb_par(:,0)*3))&
                   + (isnow2-1)*(lr%d%n1i*lr%d%n2i*lr%d%n3i*maxval(orbs%norb_par(:,0))), kind=mpi_address_kind) ,&
                   window=win2)
           else if (new_mpi_pattern .and. use_mpi_get) then
              if (.not. acc_post) then
                 call mpiwinpost(grp_acc_post(irnow2),win4,MPI_MODE_NOCHECK+MPI_MODE_NOSTORE)
                 acc_post=.true.
              end if
           else
              call MPI_IRECV(dpsiw(1,1,irnow2,igroup),&
                   nvctr_par(iproc,igrpr(igroup)),mpidtypw,jprocsr(4,jproc,igroup),&
                   jprocsr(4,jproc,igroup)+nproc+2*nproc*jproc,bigdft_mpi%mpi_comm,mpireq2(ncommsstep2),ierr)
           end if
        end if
     end do

     if(use_mpi_get) then
        if (new_mpi_pattern) then
           !here assert=MPI_MODE_NOPUT
           !call mpi_fence(win3)
           if (get_post) then
              call mpiwincomplete(win3)
              get_post=.false.
              !the get group will now become the put group
              if (jproc > 1) then
                 grp_acc_post(isnow2)=grp_put
              else
                 grp_acc_post(irnow2)=grp_put
              end if
              grp_put=mpigroup_null()
           end if
           if (get_start) then 
              call mpiwinwait(win3)
              get_start=.false.
              !the put group will now become the get group
              if (jproc > 1) then
                 grp_acc_start(irnow2)=grp_get
              else
                 grp_acc_start(isnow2)=grp_get
              end if
              grp_get=mpigroup_null()
           end if

        else
           if (jproc == 0) then
              call mpi_fence(win3)
           else
              call mpi_fence(win)
           endif
        end if
     else
        if (ncommsstep /=0) then
           !verify that the messages have been passed
           !print *,'waiting,iproc',iproc
           call MPI_WAITALL(ncommsstep,mpireq,mpistat,ierr)
           if (ierr /=0) print *,'step,ierr',jproc+1,iproc,ierr,mpistat !,MPI_STATUSES_IGNORE
           !print *,'done,iproc',iproc
        end if
     end if
     if (jproc>1) isnow2=3-isnow2
     isnow=3-isnow
     ncommsstep=0
  end do

  !call MPI_BARRIER(bigdft_mpi%mpi_comm,ierr)
  if (nproc>1) call mpiallred(eexctX,1,MPI_SUM,comm=bigdft_mpi%mpi_comm)

  !the exact exchange energy is half the Hartree energy (which already has another half)
  eexctX=-exctXfac*eexctX

  if (iproc == 0) call yaml_map('Exact Exchange Energy',eexctX,fmt='(1pe18.11)')
  !if (iproc == 0) write(*,'(1x,a,1x,1pe18.11)')'Exact Exchange Energy:',eexctX
  if(use_mpi_get) then
     if (.not. new_mpi_pattern) then
        call mpi_win_free(win, ierr)
        call mpi_win_free(win2, ierr)
     end if
     call mpi_win_free(win3, ierr)
     if (new_mpi_pattern) then
        call mpigroup_free(grp_get)
        call mpigroup_free(grp_put)
        call mpigroup_free(grp_acc_start(1))
        call mpigroup_free(grp_acc_start(2))
        call mpigroup_free(grp_acc_post(1))
        call mpigroup_free(grp_acc_post(2))
        call mpi_win_free(win4, ierr)
        call mpigroup_free(base_group)
     end if
  end if
  !close(100+iproc)
  call f_free(nvctr_par)
  call f_free(iorbgr)
  call f_free(ndatas)
  call f_free(ndatac)
  call f_free(rp_ij)
  call f_free(psiw)
  call f_free(dpsiw)
  call f_free(psir)
  call f_free(igrpr)
  call f_free(igrprarr)
  call f_free(iprocpm1)
  call f_free(jprocsr)
  !call timing(iproc,'Exchangecorr  ','OF')

END SUBROUTINE exact_exchange_potential_round

!> Calculate the exact exchange potential on occupied orbitals
!! within the symmetric round-robin scheme
!! the psi is already given in the real-space form
subroutine exact_exchange_potential_round_clean(iproc,nproc,xc,nspin,ndim,orbs,&
     pkernel,psir,dpsir,eexctX)
  use module_base
  use module_types
  use Poisson_Solver, except_dp => dp, except_gp => gp
  use module_xc
  use yaml_output
  use locreg_operations
  use overlap_point_to_point
  implicit none
  integer, intent(in) :: iproc,nproc,nspin,ndim
  !real(gp), intent(in) :: hxh,hyh,hzh
  type(xc_info), intent(in) :: xc
  type(orbitals_data), intent(in) :: orbs
!  real(wp), dimension(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f,orbs%nspinor,orbs%norbp), intent(in) :: psi
  real(wp), dimension(ndim,orbs%norbp), intent(in) :: psir
    
  type(coulomb_operator), intent(inout) :: pkernel
  real(gp), intent(out) :: eexctX
  real(wp), dimension(ndim,orbs%norbp), intent(out) :: dpsir
  !local variables
  character(len=*), parameter :: subname='exact_exchange_potential_round'
  integer, parameter :: LOCAL_=2,GLOBAL_=1
  integer, parameter :: SEND_DATA=1,RECV_DATA=2,SEND_RES=3,RECV_RES=4
  integer, parameter :: DATA_=1,RES_=2
  logical :: doit,symmetric
  integer :: ierr,ndata_comms,nres_comms,isnow,irnow,isnow2,irnow2,jsorb,norbp,source,dest,nstep_max,nsteps
  integer :: i,igroup,ngroup,ngroupp,nend,isorb,iorbs,jorbs,ii,count,istep,jproc,iorb,jorbs_tmp,norbpu,norbpd
  integer :: icount,nprocgr,iprocgrs,iprocgrr,itestproc,norbi,norbj,ncalltot,icountmax,iprocref,ncalls
  real(gp) :: ehart,exctXfac,sfac
  integer, dimension(4) :: mpireq,mpireq2
  type(workarr_sumrho) :: w
  integer, dimension(:), allocatable :: igrpr,ndatac
  integer, dimension(:,:), allocatable :: nvctr_par
  integer, dimension(:,:,:), allocatable :: jprocsr,iprocpm1,ndatas,iorbgr
  real(wp), dimension(:), allocatable :: rp_ij
!  real(wp), dimension(:,:), allocatable :: psir
  real(wp), dimension(:,:,:,:), allocatable :: psiw,dpsiw
  type(local_data) :: phi_i,phi_j


  !call timing(iproc,'Exchangecorr  ','ON')
  !decide the strategy for the communication
  symmetric=.true.

  exctXfac = xc_exctXfac(xc)

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

  !construct the OP2P scheme and test it
  !use temporaryly tyhe nvrct_par array
  nvctr_par = f_malloc0((/ 0.to.nproc-1, 1.to.ngroup /),id='nvctr_par')
  isorb=0
  do jproc=0,nproc-1
     norbp=orbs%norb_par(jproc,0)
     !transition region
     if (isorb+norbp > orbs%norbu .and. isorb < orbs%norbu) then
        nvctr_par(jproc,1)=orbs%norbu-isorb
        if (ngroup==2) nvctr_par(jproc,2)=isorb+norbp-orbs%norbu
     else if (isorb >= orbs%norbu .and. ngroup==2) then
        nvctr_par(jproc,2)=norbp
     else
        nvctr_par(jproc,1)=norbp
     end if
     isorb=isorb+norbp
  end do
!!$  call yaml_map('Orbital repartition'+yaml_toa(iproc),[orbs%norbu,orbs%norbd,sum(orbs%norb_par(:,0))])
  !if (iproc==0)call yaml_map('Orbital repartition'+yaml_toa(iproc),nvctr_par)
  if (any(sum(nvctr_par,dim=1) /= [orbs%norbu,orbs%norbd])) &
       call f_err_throw('Error in orbital repartition'+yaml_toa(iproc)+';'+yaml_toa(sum(nvctr_par,dim=1)),&
       err_name='BIGDFT_RUNTIME_ERROR')
  call OP2P_unitary_test(bigdft_mpi%mpi_comm,iproc,nproc,ngroup,ndim,nvctr_par,.true.)


  !here we can start with the round-robin scheme
  !since the orbitals are all occupied we have to use the symmetric scheme
  !we have first to define the number of groups, which correspond to the repartition 
  !of spin up and spin down orbitals
!  nvctr_par = f_malloc((/ 0.to.nproc-1, 1.to.ngroup /),id='nvctr_par')

  iorbgr = f_malloc((/ 1.to.2, 0.to.nproc-1, 1.to.ngroup /),id='iorbgr')

  if (ngroup==2) then
     isorb=0
     do jproc=0,nproc-1
        iorbgr(GLOBAL_,jproc,1)=isorb
        iorbgr(LOCAL_,jproc,1)=1
        norbp=max(min(isorb+orbs%norb_par(jproc,0),orbs%norbu)-isorb,0)
        if (norbp == 0) then
           iorbgr(GLOBAL_,jproc,1)=0
           !iorbgr(GLOBAL_,jproc,2)=0
        end if
        nvctr_par(jproc,1)=norbp*ndim
        iorbgr(GLOBAL_,jproc,2)=isorb
        iorbgr(LOCAL_,jproc,2)=norbp+1
        norbp=max(isorb+orbs%norb_par(jproc,0)-max(orbs%norbu,isorb),0)
        if (norbp == 0) then
           iorbgr(GLOBAL_,jproc,2)=0
           iorbgr(LOCAL_,jproc,2)=1
        end if
        nvctr_par(jproc,2)=norbp*ndim
        isorb=isorb+orbs%norb_par(jproc,0)
     end do
  else
     isorb=0
     do jproc=0,nproc-1
        iorbgr(GLOBAL_,jproc,1)=isorb
        iorbgr(LOCAL_,jproc,1)=1
        nvctr_par(jproc,1)=orbs%norb_par(jproc,0)*ndim
        isorb=isorb+orbs%norb_par(jproc,0)
     end do
  end if


  !> initialization

  !here we can allocate the working arrays giving the maximum
  !between the components for each group
  ngroupp=0
  do igroup=1,ngroup
     if (nvctr_par(iproc,igroup) > 0) then
        ngroupp=ngroupp+1
     end if
  end do

  !determine the array of the groups which are of interest for this processor
  igrpr = f_malloc(ngroupp,id='igrpr')
  iprocpm1 = f_malloc((/ 1.to.2, 0.to.nproc-1, 1.to.ngroupp /),id='iprocpm1')

  !test array for data calculation
  ndatac = f_malloc0(ngroupp,id='ndatac')

  !determine for each processor the groups which has to be used
  icount=0
  do igroup=1,ngroup
     if (nvctr_par(iproc,igroup) > 0) then
        icount=icount+1
        igrpr(icount)=igroup
     end if
  end do


  !find the processor whih has the maximum number of groups
  icountmax=0
  do jproc=0,nproc-1
     icount=0
     do igroup=1,ngroup
        if (nvctr_par(jproc,igroup) > 0) then
           icount=icount+1
        end if
     end do
     if (icount > icountmax) then
        iprocref=jproc
        icountmax=icount
     end if
  end do

  if (symmetric) then
     nstep_max=(nproc-1)/2+1!nproc/2+1
  else
     nstep_max=nproc-1
  end if

  !calculate the processor which lies after and before the present in the list
  iprocpm1=mpirank_null()
  do igroup=1,ngroupp
     iprocgrs=-1
     iprocgrr=-1
     !define the number of data to calculate in total
     do jproc=0,nproc-1
        ndatac(igroup)=ndatac(igroup)-nvctr_par(jproc,igrpr(igroup))
        if (nvctr_par(modulo(iproc+jproc,nproc),igrpr(igroup)) > 0 .and. .true.) then
           iprocgrs=iprocgrs+1
           iprocpm1(1,iprocgrs,igroup)=modulo(iproc+jproc,nproc)
        end if
        if (nvctr_par(modulo(iproc-jproc,nproc),igrpr(igroup)) > 0 .and. .true.) then
           iprocgrr=iprocgrr+1
           iprocpm1(2,iprocgrr,igroup)=modulo(iproc-jproc,nproc)
        end if
     end do
  end do

  !calculate the list of send-receive operations which have to be performed per group
  !allocate it at the maximum size needed
  jprocsr = f_malloc([1.to.4, 1.to.ngroupp, 0.to.nstep_max],id='jprocsr')
  !initalise array to rank_null
  jprocsr=mpirank_null()

  ncalltot=0
  do igroup=1,ngroupp
     ncalltot=ncalltot+&
          (nvctr_par(iproc,igrpr(igroup))/(ndim))*&
          (nvctr_par(iproc,igrpr(igroup))/(ndim)+1)/2
     !calculate the number of processors per group
     nprocgr=0
     do jproc=0,nproc-1
        if (nvctr_par(jproc,igrpr(igroup)) > 0) nprocgr=nprocgr+1
     end do
     !do not send anything if there is only one member
     if (nprocgr > 1) then
        if (symmetric) then
           nsteps=(nprocgr-1)/2
        else
           nsteps=nprocgr-1
        end if
        do istep=0,nsteps-1
           !define the arrays for send-receive of data
           jprocsr(SEND_DATA,igroup,istep)= iprocpm1(1,istep+1,igroup)
           jprocsr(RECV_DATA,igroup,istep)= iprocpm1(2,istep+1,igroup)
           if (iproc == iprocref .and. jprocsr(RECV_DATA,igroup,istep) /= mpirank_null()) then
              ncalltot=ncalltot+&
                   (nvctr_par(jprocsr(RECV_DATA,igroup,istep),igrpr(igroup))/(ndim))*&
                   (nvctr_par(iproc,igrpr(igroup))/(ndim))
           end if
           if (istep > 0 .and. symmetric) then
              jprocsr(SEND_RES,igroup,istep)=iprocpm1(2,istep,igroup)
              jprocsr(RECV_RES,igroup,istep)=iprocpm1(1,istep,igroup)
           end if
        end do
        !last case
        istep=nsteps!(nprocgr-1)/2
        !the last step behaves differently if the number of members is odd or even
        if (modulo(nprocgr,2) == 0 .or. .not. symmetric) then
           jprocsr(SEND_DATA,igroup,istep)= iprocpm1(1,istep+1,igroup)
           jprocsr(RECV_DATA,igroup,istep)= iprocpm1(2,istep+1,igroup)
           if (iproc == iprocref .and. jprocsr(RECV_DATA,igroup,istep) /= mpirank_null()) then
              ncalltot=ncalltot+&
                   (nvctr_par(jprocsr(RECV_DATA,igroup,istep),igrpr(igroup))/(ndim))*&
                   (nvctr_par(iproc,igrpr(igroup))/(ndim))
           end if
           if (istep > 0 .and. symmetric) then
              jprocsr(SEND_RES,igroup,istep)=iprocpm1(2,istep,igroup)
              jprocsr(RECV_RES,igroup,istep)=iprocpm1(1,istep,igroup)
           end if
        else
           jprocsr(SEND_RES,igroup,istep)=iprocpm1(2,istep,igroup)
           jprocsr(RECV_RES,igroup,istep)=iprocpm1(1,istep,igroup)
        end if
     end if
  end do


  !partial densities and potentials
  rp_ij = f_malloc(ndim,id='rp_ij')


  call f_zero(dpsir)

  if (ngroupp>0) then
     phi_i=local_data_init(orbs%norbp,ndim)
     call set_local_data(phi_i,iorbgr(GLOBAL_,iproc,igrpr(1)),&
          psir,dpsir)
  end if

  ncalls=0
  !real communication
  itestproc=mpirank_null()-1! =no debug verbosity
    
  psiw = f_malloc0([ ndim, maxval(orbs%norb_par(:,0)), ngroupp,2],id='psiw')
  dpsiw = f_malloc0([ndim, maxval(orbs%norb_par(:,0)), ngroupp,3],id='dpsiw')

  !test array for data sending
  ndatas = f_malloc0((/ 1.to.2, 0.to.nproc-1, 1.to.ngroup /),id='ndatas')

  isnow=1
  isnow2=1
  nend=(nproc-1)/2+1
  nres_comms=0


  do istep=0,nstep_max!nend
     irnow=3-isnow
     ndata_comms=0
     !sending receiving data
     do igroup=1,ngroupp
        dest=jprocsr(SEND_DATA,igroup,istep)
        if (dest /= mpirank_null()) then
           count=nvctr_par(iproc,igrpr(igroup))
           ndata_comms=ndata_comms+1
           !send the fixed array to the processor which comes in the list
           ndatas(DATA_,dest,igrpr(igroup))=ndatas(DATA_,dest,igrpr(igroup))+nvctr_par(dest,igrpr(igroup))
           call mpisend(psir(1,iorbgr(LOCAL_,iproc,igrpr(igroup))),count,&
                dest=dest,tag=iproc,comm=bigdft_mpi%mpi_comm,request=mpireq(ndata_comms),&
                verbose= dest==itestproc)
        end if

        source=jprocsr(RECV_DATA,igroup,istep)
        if (source /= mpirank_null()) then
           count=nvctr_par(source,igrpr(igroup))
           ndata_comms=ndata_comms+1
           ndatas(DATA_,iproc,igrpr(igroup))=ndatas(DATA_,iproc,igrpr(igroup))-count
           call mpirecv(psiw(1,1,igroup,irnow),count,&
                source=source,tag=source,comm=bigdft_mpi%mpi_comm,request=mpireq(ndata_comms),verbose= source == itestproc)
        end if
     end do

     do igroup=1,ngroupp
        if (istep /= 0 .and. jprocsr(SEND_RES,igroup,istep) /= mpirank_null()) then
           !put to zero the sending element
           ii=ndim*maxval(orbs%norb_par(:,0))
           call f_zero(ii,dpsiw(1,1,igroup,3))
        end if
     end do

     !calculation for orbitals to be performed
     loop_nocomm: do igroup=1,ngroupp
        if (istep == 0) then
           doit=.true.
           source=iproc
        else if (jprocsr(RECV_DATA,igroup,istep-1) /= mpirank_null()) then
           doit=.true.
           source=jprocsr(RECV_DATA,igroup,istep-1)
           if (iproc == itestproc) then
              print '(5(1x,a,i8))','step',istep+1,'group:',igrpr(igroup),&
                   ':processing',nvctr_par(source,igrpr(igroup)),&
                   'elements in',iproc,'from',source
           end if
        else
           doit=.false.
        end if
        if (doit) then
           ndatac(igroup)=ndatac(igroup)+nvctr_par(source,igrpr(igroup))
           !calculation of the partial densities and potentials
           !starting point of the loop
           !here there is the calculation routine
           !number of orbitals to be treated locally
           norbi=nvctr_par(iproc,igrpr(igroup))/(ndim)
           norbj=nvctr_par(source,igrpr(igroup))/(ndim)

           !calculating the starting orbitals locally
           iorbs=iorbgr(LOCAL_,iproc,igrpr(igroup))
           jorbs=iorbgr(LOCAL_,source,igrpr(igroup))

           !calculate the starting orbital globally
           isorb=iorbgr(GLOBAL_,iproc,igrpr(igroup))
           jsorb=iorbgr(GLOBAL_,source,igrpr(igroup))

!!$           phi_i=local_data_init(norbi,ndim)
!!$           call set_local_data(phi_i,isorb+iorbs-1,&
!!$                psir(1,iorbs),dpsir(1,iorbs))

           if (istep/=0) then
              phi_j=local_data_init(norbj,ndim)
              call set_local_data(phi_j,jsorb+jorbs-1,&
                   psiw(1,1,igroup,isnow),dpsiw(1,1,igroup,3))
              jorbs_tmp=1
           else
              phi_j=phi_i
              jorbs_tmp=iorbs
           end if

           call internal_calculation_exctx(istep,sfac,pkernel,orbs%norb,orbs%occup,orbs%spinsgn,&
                jprocsr(SEND_RES,igroup,istep) /= mpirank_null(),norbi,norbj,iorbs,jorbs_tmp,&
                phi_i,phi_j,eexctX,rp_ij)
           if (iproc == iprocref .and. verbose > 1 .and. igroup==1) then
              if (istep == 0) then
                 ncalls=ncalls+((norbi-iorbs+1)*(norbj-jorbs+1+1))/2
              else
                 ncalls=ncalls+(norbi-iorbs+1)*(norbj-jorbs+1)
              end if
              call yaml_comment('Exact exchange calculation: '+nint(real(ncalls,gp)/real(ncalltot,gp)*100.0_gp)**'(i3)'+'%')
           end if

!!$           call free_local_data(phi_i)
           if (istep/=0) call free_local_data(phi_j)

        end if
     end do loop_nocomm


     if (nres_comms > 0) then
        !verify that the messages have been passed
        call mpiwaitall(nres_comms,mpireq2)
        !copy the results which have been received (the messages sending are after)
        !this part is already done by the mpi_accumulate
        do igroup=1,ngroupp
           source=jprocsr(RECV_RES,igroup,istep-1)
           if (source /= mpirank_null()) then
              if (iproc == itestproc) then
                 print '(5(1x,a,i8))','step',istep+1,'group:',igrpr(igroup),&
                      ':copying',nvctr_par(source,igrpr(igroup)),&
                      'processed elements from',source,'in',iproc
              end if
              ndatac(igroup)=ndatac(igroup)+nvctr_par(source,igrpr(igroup))
              !WARNING: should here source == iproc?
              call axpy(nvctr_par(iproc,igrpr(igroup)),1.0_wp,dpsiw(1,1,igroup,irnow2),1,&
                   dpsir(1,iorbgr(LOCAL_,iproc,igrpr(igroup))),1)
           end if
        end do
     end if
     nres_comms=0
     !meanwhile, we can receive the result from the processor which has the psi 

     irnow2=3-isnow2
     do igroup=1,ngroupp
        dest=jprocsr(SEND_RES,igroup,istep)
        if (dest /= mpirank_null()) then
           nres_comms=nres_comms+1
           call f_memcpy(n=nvctr_par(dest,igrpr(igroup)),src=dpsiw(1,1,igroup,3),&
                dest=dpsiw(1,1,igroup,isnow2))
           call mpisend(dpsiw(1,1,igroup,isnow2),&
                nvctr_par(dest,igrpr(igroup)),&
                dest=dest,tag=iproc+nproc+2*nproc*istep,comm=bigdft_mpi%mpi_comm,&
                request=mpireq2(nres_comms))
           ndatas(RES_,dest,igrpr(igroup))=ndatas(RES_,dest,igrpr(igroup))+nvctr_par(dest,igrpr(igroup))
        end if
     end do
     do igroup=1,ngroupp
        source=jprocsr(RECV_RES,igroup,istep)
        if (source /= mpirank_null()) then
           nres_comms=nres_comms+1
           call mpirecv(dpsiw(1,1,igroup,irnow2),&
                nvctr_par(iproc,igrpr(igroup)),source=source,&
                tag=source+nproc+2*nproc*istep,comm=bigdft_mpi%mpi_comm,request=mpireq2(nres_comms))
           ndatas(RES_,iproc,igrpr(igroup))=ndatas(RES_,iproc,igrpr(igroup))-nvctr_par(iproc,igrpr(igroup))
        end if
     end do

     !verify that the messages have been passed
     call mpiwaitall(ndata_comms,mpireq)

     if (istep>1) isnow2=3-isnow2
     isnow=3-isnow
     ndata_comms=0
  end do

  if (ngroupp>0) call free_local_data(phi_i)

  !unitary tests of the calculation
  if (itestproc > -1) then
     if (nproc > 1) call mpiallred(ndatas,MPI_SUM,comm=bigdft_mpi%mpi_comm)
     !if(iproc ==0)print *,'iproc,datas',iproc,ndatas
     do igroup=1,ngroupp
        if (ndatac(igroup) /=0) then
           write(*,*)'ERROR: OP2P communication simulation failed: processor',iproc,&
                ' has calculated',ndatac(igroup),' data more than needed'
           stop
        end if
        if (ndatas(DATA_,iproc,igrpr(igroup)) /=0 .or. ndatas(RES_,iproc,igrpr(igroup)) /=0) then
           write(*,*)'ERROR: OP2P communication simulation failed: processor',iproc,&
                ' has not a zero balance of send-receive calls',ndatas(:,iproc,igrpr(igroup))
           stop
        end if
     end do
  end if

  !call MPI_BARRIER(bigdft_mpi%mpi_comm,ierr)
  if (nproc>1) call mpiallred(eexctX,1,MPI_SUM,comm=bigdft_mpi%mpi_comm)

  !the exact exchange energy is half the Hartree energy (which already has another half)
  eexctX=-exctXfac*eexctX

  if (iproc == 0) call yaml_map('Exact Exchange Energy',eexctX,fmt='(1pe18.11)')
  !if (iproc == 0) write(*,'(1x,a,1x,1pe18.11)')'Exact Exchange Energy:',eexctX
  !close(100+iproc)
  call f_free(nvctr_par)
  call f_free(iorbgr)
  call f_free(ndatas)
  call f_free(ndatac)
  call f_free(rp_ij)
  call f_free(psiw)
  call f_free(dpsiw)
!  call f_free(psir)
  call f_free(igrpr)
  call f_free(iprocpm1)
  call f_free(jprocsr)
  !call timing(iproc,'Exchangecorr  ','OF')

END SUBROUTINE exact_exchange_potential_round_clean


!!$subroutine rma_group(base_group,ranks,grp)
!!$  implicit none
!!$  integer, intent(in) :: base_group
!!$  integer, dimension(3), intent(in) :: ranks
!!$  integer, intent(out) :: grp
!!$  !local variables
!!$  integer :: nlist,ilist,center,irank
!!$  integer, dimension(1) :: minl,maxl
!!$  integer, dimension(3) :: list,ipiv
!!$
!!$  !order the list
!!$  list=ranks
!!$  minl=minloc(list)
!!$  list(minl(1))=-2 ! exclude min
!!$  maxl=maxloc(list)
!!$  center=6-minl(1)-maxl(1) !take the other
!!$  nlist=0
!!$  ipiv=[minl(1),center,maxl(1)]
!!$  do ilist=1,3
!!$     irank=ranks(ipiv(ilist))
!!$     if (irank >= 0) then
!!$        nlist=nlist+1
!!$        list(nlist)=irank
!!$     end if
!!$  end do
!!$
!!$  if (nlist > 0) then
!!$     grp=mpigroupincl(base_group,nlist,list)
!!$  else
!!$     grp=mpigroup_null()
!!$  end if
!!$
!!$end subroutine rma_group

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
!!$  character(len=1), intent(in) :: geocode !< @copydoc poisson_solver::doc::geocode
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
!!$  integer :: ierr,ispin,ncommsstep,ncommsstep2,isnow,irnow,isnow2,irnow2,jsorb,kproc,norbp,jgroup
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
!!$  call initialize_work_arrays_sumrho(1,lr,.true.,w)
!!$  allocate(psir(lr%d%n1i*lr%d%n2i*lr%d%n3i,orbs%norbp+ndebug),stat=i_stat)
!!$  call memocc(i_stat,psir,'psir',subname)
!!$  
!!$  call to_zero(lr%d%n1i*lr%d%n2i*lr%d%n3i*orbs%norbp,psir)
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
!!$  call to_zero(lr%d%n1i*lr%d%n2i*lr%d%n3i*maxval(orbs%norb_par)*2*ngroupp,psiw)
!!$  call to_zero(lr%d%n1i*lr%d%n2i*lr%d%n3i*maxval(orbs%norb_par)*3*ngroupp,dpsiw)
!!$  
!!$  call to_zero(lr%d%n1i*lr%d%n2i*lr%d%n3i*orbs%norbp,dpsir)
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
!!$           call to_zero(lr%d%n1i*lr%d%n2i*lr%d%n3i*maxval(orbs%norb_par),dpsiw(1,1,3,igroup))
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
!!$           call vcopy(nvctr_par(jprocsr(3,jproc,igroup),igrpr(igroup)),&
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
