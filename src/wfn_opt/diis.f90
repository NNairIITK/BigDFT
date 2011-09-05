!> @file
!!   Routines for density mixing and wavefunction update
!! @author
!!    Copyright (C) 2007-2011 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 

!> Extract the energy (the quantity which has to be minimised by the wavefunction)
!! and calculate the corresponding gradient.
!! The energy can be the actual Kohn-Sham energy or the trace of the hamiltonian, 
!! depending of the functional we want to calculate. The gradient wrt the wavefucntion
!! Is then put in hpsi accordingly to the functional
subroutine calculate_energy_and_gradient_new(iter,iproc,nproc,orbs,comms,GPU,lr,orthpar,hx,hy,hz,ncong,iscf,&
     energs,psi,psit,hpsi,gnrm,gnrm_zero,energy)
  use module_base
  use module_types
  use module_interfaces!, except_this_one => calculate_energy_and_gradient_new
  implicit none
  integer, intent(in) :: iproc,nproc,ncong,iscf,iter
  real(gp), intent(in) :: hx,hy,hz
  type(orbitals_data), intent(inout) :: orbs
  type(communications_arrays), intent(in) :: comms
  type(locreg_descriptors), intent(in) :: lr
  type(GPU_pointers), intent(in) :: GPU
  type(orthon_data), intent(in) :: orthpar
  type(energy_terms), intent(out) :: energs
  real(gp), intent(out) :: gnrm,gnrm_zero,energy
  real(wp), dimension(:), pointer :: psi,psit,hpsi
  !local variables
  character(len=*), parameter :: subname='calculate_energy_and_gradient' 
  logical :: lcs
  integer :: ierr,ikpt,iorb,i_all,i_stat,k,ncplx,jorb
  real(gp) :: energybs,trH,rzeroorbs,tt,energyKS
  real(wp), dimension(:), allocatable :: passmat
  real(wp), dimension(:,:,:), allocatable :: mom_vec
  real(wp), dimension(:), pointer :: psiw !< fake pointer for the calculation of the hamiltonian matrix

  !band structure energy calculated with occupation numbers
  energs%ebs=energs%ekin+energs%epot+energs%eproj !the potential energy contains also exctX
  !this is the Kohn-Sham energy
  energs%eKS=energs%ebs-energs%eh+energs%exc-energs%vxc-energs%eexctX+energs%eion+energs%edisp

  !calculate orbital poloarisation directions
  if(orbs%nspinor==4) then
     allocate(mom_vec(4,orbs%norb,min(nproc,2)+ndebug),stat=i_stat)
     call memocc(i_stat,mom_vec,'mom_vec',subname)

     call calc_moments(iproc,nproc,orbs%norb,orbs%norb_par,&
          lr%wfd%nvctr_c+7*lr%wfd%nvctr_f,orbs%nspinor,psi,mom_vec)
  end if


  if (iproc==0 .and. verbose > 1) then
     write(*,'(1x,a)',advance='no')&
          'done,  orthoconstraint...'
  end if

!  !transpose the hpsi wavefunction
!  call transpose_v(iproc,nproc,orbs,lr%wfd,comms,hpsi,work=psi)
!
  !transpose the psit wavefunction for the non-collinear case
  if (nproc == 1) then
     !associate psit pointer for orthoconstraint and transpose it (for the non-collinear case)
     psit => psi
     call transpose_v(iproc,nproc,orbs,lr%wfd,comms,psit)
  end if

  !calculate overlap matrix of the psi-hamiltonian to find the passage matrix for rho
  !wavefunctions are in orbital distribution at the beginning and at the end
  !the eigenvalues in orbs are overwritten with the diagonal matrix
  !psit is the new transposed wavefunction 
  !hpsi is destroyed

  !allocate the passage matrix for transforming the LCAO wavefunctions in the IG wavefucntions
  ncplx=1
  if (orbs%nspinor > 1) ncplx=2
  allocate(passmat(ncplx*orbs%nkptsp*(orbs%norbu*orbs%norbu+orbs%norbd*orbs%norbd)+ndebug),stat=i_stat)
  call memocc(i_stat,passmat,'passmat',subname)

  allocate(psiw(orbs%npsidim+ndebug),stat=i_stat)
  call memocc(i_stat,psiw,'psiw',subname)

  !put the hpsi wavefunction in the work array
  call dcopy(orbs%npsidim,hpsi,1,psiw,1)

  call DiagHam(iproc,nproc,0,orbs%nspin,orbs,lr%wfd,comms,&
     psi,psiw,psit,orthpar,passmat)

  !print out the passage matrix (valid for one k-point only and ncplx=1)
  do iorb=1,orbs%norbu
     write(*,'(100(1pe14.7))')(passmat(iorb+orbs%norbu*(jorb-1)),jorb=1,orbs%norbu)
  end do

  !DEBUG: test 
  !this passage matrix should be applied to the 
  !the second passage should give a passage matrix of identity
  !the matrix

!!$
!!$
!!$  norbi_max=max(orbs%norbu,orbs%norbd) 
!!$  ndim_hamovr=norbi_max**2
!!$  !for complex matrices the dimension is doubled
!!$  if (nspinor /=1) then
!!$     ndim_hamovr=2*ndim_hamovr
!!$  end if
!!$  allocate(norbgrp(1,nspin+ndebug),stat=i_stat)
!!$  call memocc(i_stat,norbgrp,'norbgrp',subname)
!!$  norbgrp(1,1)=orbs%norbu
!!$  norbgrp(1,2)=orbs%norbd
!!$
!!$  allocate(hamovr(nspin*ndim_hamovr,2,orbsu%nkpts+ndebug),stat=i_stat)
!!$  call memocc(i_stat,hamovr,'hamovr',subname)
!!$
!!$  !initialise hamovr
!!$  call razero(nspin*ndim_hamovr*2*orbsu%nkpts,hamovr)
!!$  ispsi=1
!!$  do ikptp=1,orbsu%nkptsp
!!$     ikpt=orbsu%iskpts+ikptp!orbsu%ikptsp(ikptp)
!!$     
!!$     nvctrp=commu%nvctr_par(iproc,ikptp)
!!$     if (nvctrp == 0) cycle
!!$     
!!$     !print *,'iproc,nvctrp,nspin,norb,ispsi,ndimovrlp',iproc,nvctrp,nspin,norb,ispsi,ndimovrlp(ispin,ikpt-1)
!!$     call overlap_matrices(norbtot,nvctrp,natsceff,nspin,nspinor,&
!!$          & ndim_hamovr,norbgrp,hamovr(1,1,ikpt),psi(ispsi),hpsi(ispsi))
!!$     
!!$     ispsi=ispsi+nvctrp*norbtot*orbsu%nspinor
!!$  end do
!!$  if (nproc > 1) then
!!$     !reduce the overlap matrix between all the processors
!!$     call mpiallred(hamovr(1,1,1),2*nspin*ndim_hamovr*orbsu%nkpts,&
!!$          MPI_SUM,MPI_COMM_WORLD,ierr)
!!$  end if
!!$
!!$  !diagonalize hamovr
  


!!$  i_all=-product(shape(hamovr))*kind(hamovr)
!!$  deallocate(hamovr,stat=i_stat)
!!$  call memocc(i_stat,i_all,'hamovr',subname)
!!$  i_all=-product(shape(norbgrp))*kind(norbgrp)
!!$  deallocate(norbgrp,stat=i_stat)
!!$  call memocc(i_stat,i_all,'norbgrp',subname)


  !after having applied the hamiltonian to all the atomic orbitals
  !we split the semicore orbitals from the valence ones
  !this is possible since the semicore orbitals are the first in the 
  !order, so the linear algebra on the transposed wavefunctions 
  !may be splitted

!!$  ispsi=1
!!$  do ikptp=1,orbsu%nkptsp
!!$     ikpt=orbsu%iskpts+ikptp!orbsu%ikptsp(ikptp)
!!$     
!!$     nvctrp=commu%nvctr_par(iproc,ikptp)
!!$     if (nvctrp == 0) cycle
!!$     
!!$     !print *,'iproc,nvctrp,nspin,norb,ispsi,ndimovrlp',iproc,nvctrp,nspin,norb,ispsi,ndimovrlp(ispin,ikpt-1)
!!$     call overlap_matrices(norbtot,nvctrp,natsceff,nspin,nspinor,&
!!$          & ndim_hamovr,norbgrp,hamovr(1,1,ikpt),psi(ispsi),hpsi(ispsi))
!!$     
!!$     ispsi=ispsi+nvctrp*norbtot*orbsu%nspinor
!!$  end do


  ! Apply  orthogonality constraints to all orbitals belonging to iproc
  !takes also into account parallel k-points distribution
  !here the orthogonality with respect to other occupied functions should be 
  !passed as an optional argument
  call orthoconstraint(iproc,nproc,orbs,comms,lr%wfd,psit,hpsi,trH)

  !retranspose the hpsi wavefunction
  call untranspose_v(iproc,nproc,orbs,lr%wfd,comms,hpsi,work=psi)

  !after having calcutated the trace of the hamiltonian, the functional have to be defined
  !new value without the trace, to be added in hpsitopsi
  if (iscf >1) then
     energy=energs%trH
  else
     energy=energs%trH-energs%eh+energs%exc-energs%vxc-energs%eexctX+energs%eion+energs%edisp
  end if

  !check that the trace of the hamiltonian is compatible with the 
  !band structure energy 
  !this can be done only if the occupation numbers are all equal
  tt=(energs%ebs-energs%trH)/energs%trH
  if (((abs(tt) > 1.d-10 .and. .not. GPUconv) .or.&
       (abs(tt) > 1.d-8 .and. GPUconv)) .and. iproc==0) then 
     !write this warning only if the system is closed shell
     call check_closed_shell(orbs,lcs)
     if (lcs) then
        write( *,'(1x,a,1pe9.2,2(1pe22.14))') &
             'ERROR: inconsistency between gradient and energy',tt,energs%ebs,energs%trH
     end if
  endif

  call timing(iproc,'Precondition  ','ON')
  if (iproc==0 .and. verbose > 1) then
     write(*,'(1x,a)',advance='no')&
          'done,  preconditioning...'
  end if

  !Preconditions all orbitals belonging to iproc
  !and calculate the partial norm of the residue
  !switch between CPU and GPU treatment
  if (GPUconv) then
     call preconditionall_GPU(iproc,nproc,orbs,lr,hx,hy,hz,ncong,&
          hpsi,gnrm,gnrm_zero,GPU)
  else if (OCLconv) then
     call preconditionall_OCL(iproc,nproc,orbs,lr,hx,hy,hz,ncong,&
          hpsi,gnrm,gnrm_zero,GPU)
  else
     call preconditionall(iproc,nproc,orbs,lr,hx,hy,hz,ncong,hpsi,gnrm,gnrm_zero)
  end if

  !sum over all the partial residues
  if (nproc > 1) then
     call mpiallred(gnrm,1,MPI_SUM,MPI_COMM_WORLD,ierr)
     call mpiallred(gnrm_zero,1,MPI_SUM,MPI_COMM_WORLD,ierr)
  endif

  !count the number of orbitals which have zero occupation number
  !weight this with the corresponding k point weight
  rzeroorbs=0.0_gp
  do ikpt=1,orbs%nkpts
     do iorb=1,orbs%norb
        if (orbs%occup(iorb+(ikpt-1)*orbs%norb) == 0.0_gp) then
           rzeroorbs=rzeroorbs+orbs%kwgts(ikpt)
        end if
     end do
  end do
  !commented out, the kwgts sum already to one
  !if (orbs%nkpts > 1) nzeroorbs=nint(real(nzeroorbs,gp)/real(orbs%nkpts,gp))

  gnrm=sqrt(gnrm/(real(orbs%norb,gp)-rzeroorbs))

  if (rzeroorbs /= 0.0_gp) then
     gnrm_zero=sqrt(gnrm_zero/rzeroorbs)
  else
     gnrm_zero=0.0_gp
  end if

  if (iproc==0 .and. verbose > 1) then
     write(*,'(1x,a)')&
          'done.'
  end if
  call timing(iproc,'Precondition  ','OF')

  i_all=-product(shape(passmat))*kind(passmat)
  deallocate(passmat,stat=i_stat)
  call memocc(i_stat,i_all,'passmat',subname)

  i_all=-product(shape(psiw))*kind(psiw)
  deallocate(psiw,stat=i_stat)
  call memocc(i_stat,i_all,'psiw',subname)



  if (orbs%nspinor == 4) then
     !only the root process has the correct array
     if(iproc==0 .and. verbose > 0) then
        write(*,'(1x,a)')&
             'Magnetic polarization per orbital'
        write(*,'(1x,a)')&
             '  iorb    m_x       m_y       m_z'
        do iorb=1,orbs%norb
           write(*,'(1x,i5,3f10.5)') &
                iorb,(mom_vec(k,iorb,1)/mom_vec(1,iorb,1),k=2,4)
        end do
     end if
     i_all=-product(shape(mom_vec))*kind(mom_vec)
     deallocate(mom_vec,stat=i_stat)
     call memocc(i_stat,i_all,'mom_vec',subname)
  end if

  !write the energy information
  if (iproc == 0) then
     if (verbose > 0 .and. iscf<1) then
        write( *,'(1x,a,3(1x,1pe18.11))') 'ekin_sum,epot_sum,eproj_sum',  & 
             energs%ekin,energs%epot,energs%eproj
        write( *,'(1x,a,3(1x,1pe18.11))') '   ehart,   eexcu,    vexcu',energs%eh,energs%exc,energs%vxc
     end if
     if (iscf > 1) then
        if (gnrm_zero == 0.0_gp) then
           write( *,'(1x,a,i6,2x,1pe24.17,1x,1pe9.2)') 'iter, tr(H),gnrm',iter,energs%trH,gnrm
        else
           write( *,'(1x,a,i6,2x,1pe24.17,2(1x,1pe9.2))') 'iter, tr(H),gnrm,gnrm_zero',iter,energs%trH,gnrm,gnrm_zero
        end if
     else
        if (gnrm_zero == 0.0_gp) then
           write( *,'(1x,a,i6,2x,1pe24.17,1x,1pe9.2)') 'iter,total energy,gnrm',iter,energs%eKS,gnrm
        else
           write( *,'(1x,a,i6,2x,1pe24.17,2(1x,1pe9.2))') 'iter,total energy,gnrm,gnrm_zero',iter,energs%eKS,gnrm,gnrm_zero
        end if
     end if
  endif

end subroutine calculate_energy_and_gradient_new


!> Mix the electronic density or the potential using DIIS
subroutine mix_rhopot(iproc,nproc,npoints,alphamix,mix,rhopot,istep,&
     & n1,n2,n3,ucvol,rpnrm,nscatterarr)
  use module_base
  use defs_basis, only: AB6_NO_ERROR
  use m_ab6_mixing
  implicit none
  integer, intent(in) :: npoints, istep, n1, n2, n3, nproc, iproc
  real(gp), intent(in) :: alphamix, ucvol
  integer, intent(in) :: nscatterarr(0:nproc-1,4)
  type(ab6_mixing_object), intent(inout) :: mix
  real(dp), dimension(npoints), intent(inout) :: rhopot
  real(gp), intent(out) :: rpnrm
  !local variables
  integer :: ierr,ie,ii,i_stat,i_all
  character(len = *), parameter :: subname = "mix_rhopot"
  character(len = 500) :: errmess
  integer, allocatable :: user_data(:)

  ! Calculate the residue and put it in rhopot
  if (istep > 1) then
     ! rhopot = vin
     call axpy(npoints, -1.d0, mix%f_fftgr(1,1, mix%i_vrespc(1)), 1, &
          & rhopot(1), 1)
     call dscal(npoints, 1.d0 - alphamix, rhopot(1), 1)
     ! rhopot = alpha(vin - v(out-1))
  else
     mix%f_fftgr(:,:, mix%i_vrespc(1)) = 0.d0
  end if
  ! rhopot = v(out-1) and fftgr = alpha(vin - v(out-1))
  call dswap(npoints, mix%f_fftgr(1,1, mix%i_vrespc(1)), 1, &
       & rhopot(1), 1)

  ! Store the scattering of rho in user_data
  allocate(user_data(2 * nproc), stat = i_stat)
  call memocc(i_stat,user_data,'user_data',subname)
  do ii = 1, nproc, 1
     user_data(1 + (ii - 1 ) * 2:ii * 2) = &
          & n1 * n2 * (/ nscatterarr(iproc, 2), nscatterarr(iproc, 4) /)
  end do

  ! Do the mixing
  call ab6_mixing_eval(mix, rhopot, istep, n1 * n2 * n3, ucvol, &
       & MPI_COMM_WORLD, (nproc > 1), ierr, errmess, resnrm = rpnrm, &
       & fnrm = fnrm_denpot, fdot = fdot_denpot, user_data = user_data)
  if (ierr /= AB6_NO_ERROR) then
     if (iproc == 0) write(0,*) errmess
     call MPI_ABORT(MPI_COMM_WORLD, ierr, ie)
  end if
  rpnrm = sqrt(rpnrm) / real(n1 * n2 * n3, gp)
  rpnrm = rpnrm / (1.d0 - alphamix)

  i_all=-product(shape(user_data))*kind(user_data)
  deallocate(user_data,stat=i_stat)
  call memocc(i_stat,i_all,'user_data',subname)

  ! Copy new in vrespc
  call dcopy(npoints, rhopot(1), 1, mix%f_fftgr(1,1, mix%i_vrespc(1)), 1)

END SUBROUTINE mix_rhopot


subroutine psimix(iproc,nproc,orbs,comms,diis,hpsit,psit, quiet)
  use module_base
  use module_types
  use module_interfaces, except_this_one => psimix
  implicit none
  integer, intent(in) :: iproc,nproc
  type(orbitals_data), intent(in) :: orbs
  type(communications_arrays), intent(in) :: comms
  type(diis_objects), intent(inout) :: diis
  real(wp), dimension(sum(comms%ncntt(0:nproc-1))), intent(inout) :: psit,hpsit
  logical, optional:: quiet ! to avoid that the DIIS weights are written
  !real(wp), dimension(:), pointer :: psit,hpsit
  !local variables
  integer :: ikptp,nvctrp,ispsi,ispsidst,jj

  if (diis%idsx > 0) then
     !do not transpose the hpsi wavefunction into the diis array
     !for compatibility with the k-points distribution
     ispsi=1
     ispsidst=1
     do ikptp=1,orbs%nkptsp
        nvctrp=comms%nvctr_par(iproc,ikptp)
        if (nvctrp == 0) cycle
     !here we can choose to store the DIIS arrays with single precision
     !psidst=psit
        call dcopy(nvctrp*orbs%norb*orbs%nspinor,&
             psit(ispsi),1,&
             diis%psidst(ispsidst+nvctrp*orbs%nspinor*orbs%norb*(diis%mids-1)),1)

     !hpsidst=hpsi
     !   call dcopy(nvctrp*orbs%norb*orbs%nspinor,&
     !        hpsit(ispsi),1,&
     !        hpsidst(ispsidst+nvctrp*orbs%nspinor*orbs%norb*(mids-1)),1)

     do jj=0,nvctrp*orbs%norb*orbs%nspinor-1
        diis%hpsidst(ispsidst+nvctrp*orbs%nspinor*orbs%norb*(diis%mids-1)+jj)&
             =real(hpsit(ispsi+jj),tp) !diis precision conversion
     end do
        ispsi=ispsi+nvctrp*orbs%norb*orbs%nspinor
        ispsidst=ispsidst+nvctrp*orbs%norb*orbs%nspinor*diis%idsx
     end do

     !here we should separate between up and down spin orbitals, maybe
     if(present(quiet)) then
         call diisstp(iproc,nproc,orbs,comms,diis,psit,quiet)
     else
         call diisstp(iproc,nproc,orbs,comms,diis,psit)
     end if

  else
     ! update all wavefunctions with the preconditioned gradient
     if (diis%energy > diis%energy_old) then
        diis%alpha=max(5.d-2,.5_wp*diis%alpha)
        if (diis%alpha == 5.d-2) write(*,*) ' WARNING: Convergence problem or limit'
     else
        diis%alpha=min(1.05_wp*diis%alpha,diis%alpha_max)
     endif
     if (iproc == 0 .and. verbose > 0) write(*,'(1x,a,1pe11.3)') 'alpha=',diis%alpha

!!     do iorb=1,orbs%norb*orbs%nspinor
!!        call axpy(comms%nvctr_par(iproc),&
!!             -alpha,hpsi(1+comms%nvctr_par(iproc)*(iorb-1)),1,&
!!             psit(1+comms%nvctr_par(iproc)*(iorb-1)),1)
!!     enddo

     call axpy(sum(comms%ncntt(0:nproc-1)),-diis%alpha,hpsit(1),1,psit(1),1)

  endif

END SUBROUTINE psimix


subroutine diis_or_sd(iproc,idsx,nkptsp,diis)
  use module_base
  use module_types
  implicit none
  integer, intent(in) :: iproc,idsx,nkptsp
  type(diis_objects), intent(inout) :: diis

  !here we should add the DIIS/SD switch
  !add switch between DIIS and SD
  if (diis%energy == diis%energy_min .and. .not. diis%switchSD) diis%idiistol=0
  if (diis%energy > diis%energy_min .and. idsx >0 .and. .not. diis%switchSD) then
     diis%idiistol=diis%idiistol+1
  end if
  if (diis%idiistol > idsx .and. .not. diis%switchSD) then
     !the energy has not decreasing for too much steps, switching to SD for next steps
     if (iproc ==0) write(*,'(1x,a,1pe9.2,a)')&
          'WARNING: The energy value is growing (delta=',diis%energy-diis%energy_min,') switch to SD'
     diis%switchSD=.true.
     !il is not necessary to deallocate these arrays
     !i_all=-product(shape(psidst))*kind(psidst)
     !deallocate(psidst,stat=i_stat)
     !call memocc(i_stat,i_all,'psidst',subname)
     !i_all=-product(shape(hpsidst_sp))*kind(hpsidst_sp)
     !deallocate(hpsidst_sp,stat=i_stat)
     !call memocc(i_stat,i_all,'hpsidst_sp',subname)
     !i_all=-product(shape(ads))*kind(ads)
     !deallocate(ads,stat=i_stat)
     !call memocc(i_stat,i_all,'ads',subname)
     diis%idsx=0
     diis%idiistol=0
  end if

  if ((diis%energy == diis%energy_min) .and. diis%switchSD) then
     diis%idiistol=diis%idiistol+1
  end if
!  if (diis%idiistol > idsx .and. diis%switchSD) then
   if (diis%idiistol > 10000*idsx .and. diis%switchSD) then
     !restore the original DIIS
     if (iproc ==0) write(*,'(1x,a,1pe9.2)')&
          'WARNING: The energy value is now decreasing again, coming back to DIIS'
     diis%switchSD=.false.
     diis%idsx=idsx
     diis%ids=0
     diis%idiistol=0

     !no need to reallocate
     !allocate(psidst(sum(comms%ncntt(0:nproc-1))*idsx+ndebug),stat=i_stat)
     !call memocc(i_stat,psidst,'psidst',subname)
     !allocate(hpsidst_sp(sum(comms%ncntt(0:nproc-1))*idsx+ndebug),stat=i_stat)
     !call memocc(i_stat,hpsidst_sp,'hpsidst_sp',subname)
     !allocate(ads(idsx+1,idsx+1,orbs%nkptsp*3+ndebug),stat=i_stat)
     !call memocc(i_stat,ads,'ads',subname)

     call razero(nkptsp*3*(idsx+1)**2,diis%ads)
  end if

END SUBROUTINE diis_or_sd


!> calculates the DIIS extrapolated solution psit in the ids-th DIIS step 
!! using  the previous iteration points psidst and the associated error 
!! vectors (preconditioned gradients) hpsidst
subroutine diisstp(iproc,nproc,orbs,comms,diis,psit,quiet)
  use module_base
  use module_types
  implicit none
! Arguments
  integer, intent(in) :: nproc,iproc
  type(orbitals_data), intent(in) :: orbs
  type(communications_arrays), intent(in) :: comms
  type(diis_objects), intent(inout) :: diis
  real(wp), dimension(sum(comms%ncntt(0:nproc-1))), intent(out) :: psit
  logical, optional:: quiet ! to avoid that the DIIS weights are written
! Local variables
  character(len=*), parameter :: subname='diisstp'
  integer :: i,j,ist,jst,mi,iorb,info,jj,mj,k,i_all,i_stat,ierr
  integer :: ikptp,ikpt,ispsi,ispsidst,nvctrp
  integer, dimension(:), allocatable :: ipiv
  real(dp), dimension(:,:), allocatable :: rds
  logical:: quietHere
  ! This decides whether the DIIS weights are written to the screen or not.
  ! If quietHere is false, then they are written, otherwise not.
  if(present(quiet) .and. quiet) then
      quietHere=.true.
  else
      quietHere=.false.
  end if

  allocate(ipiv(diis%idsx+1+ndebug),stat=i_stat)
  call memocc(i_stat,ipiv,'ipiv',subname)
  allocate(rds(diis%idsx+1,orbs%nkpts+ndebug),stat=i_stat)
  call memocc(i_stat,rds,'rds',subname)

  call razero((diis%idsx+1)*orbs%nkpts,rds)

  ispsidst=1
  do ikptp=1,orbs%nkptsp
     ikpt=orbs%iskpts+ikptp!orbs%ikptsp(ikptp)
     nvctrp=comms%nvctr_par(iproc,ikptp)
     if (nvctrp == 0) cycle

     ! set up DIIS matrix (upper triangle)
     if (diis%ids > diis%idsx) then
        ! shift left up matrix
        do i=1,diis%idsx-1
           do j=1,i
              diis%ads(j,i,ikptp,1)=diis%ads(j+1,i+1,ikptp,1)
           end do
        end do
     end if

     ! calculate new line, use rds as work array for summation
     ist=max(1,diis%ids-diis%idsx+1)
     do i=ist,diis%ids
        mi=mod(i-1,diis%idsx)+1
!!     do iorb=1,norb*nspinor
!!        tt=dot(nvctrp,hpsidst(1,iorb,mids),1,hpsidst(1,iorb,mi),1)
!!        rds(i-ist+1)=rds(i-ist+1)+tt
!!     end do
        !to be corrected for complex wavefunctions
        rds(i-ist+1,ikpt)=dot(nvctrp*orbs%norb*orbs%nspinor,&
             diis%hpsidst(ispsidst+(diis%mids-1)*nvctrp*orbs%norb*orbs%nspinor),1,&
             diis%hpsidst(ispsidst+(mi-1)*nvctrp*orbs%norb*orbs%nspinor),1)
        !this has to be inserted in module_base
        !call ds_dot(nvctrp*orbs%norb*orbs%nspinor,&
        !     hpsidst_sp,ispsidst+(mids-1)*nvctrp*orbs%norb*orbs%nspinor,1,&
        !     hpsidst_sp,ispsidst+(mi  -1)*nvctrp*orbs%norb*orbs%nspinor,1,&
        !     rds(i-ist+1,ikpt))
     end do
     ispsidst=ispsidst+nvctrp*orbs%norb*orbs%nspinor*diis%idsx
  end do

  if (nproc > 1) then
     call mpiallred(rds(1,1),(diis%idsx+1)*orbs%nkpts,MPI_SUM,MPI_COMM_WORLD,ierr)
!     call MPI_ALLREDUCE(MPI_IN_PLACE,rds,(diis%idsx+1)*orbs%nkpts,  & 
!                 mpidtypw,MPI_SUM,MPI_COMM_WORLD,ierr)
!
!!     call MPI_ALLREDUCE(rds,ads(1,min(diis%idsx,ids),1),min(ids,diis%idsx),  & 
!!                 MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
  endif
  
  ispsi=1
  ispsidst=1
  do ikptp=1,orbs%nkptsp
     ikpt=orbs%iskpts+ikptp!orbs%ikptsp(ikptp)
     nvctrp=comms%nvctr_par(iproc,ikptp)
     if (nvctrp == 0) cycle

     do i=1,min(diis%ids,diis%idsx)
        diis%ads(i,min(diis%idsx,diis%ids),ikptp,1)=rds(i,ikpt)
     end do

     ! copy to work array, right hand side, boundary elements
     do j=1,min(diis%idsx,diis%ids)
        diis%ads(j,min(diis%idsx,diis%ids)+1,ikptp,2)=1.0_wp
        rds(j,ikpt)=0.d0
        do i=j,min(diis%idsx,diis%ids)
           diis%ads(j,i,ikptp,2)=diis%ads(j,i,ikptp,1)
        end do
     end do
     diis%ads(min(diis%idsx,diis%ids)+1,min(diis%idsx,diis%ids)+1,ikptp,2)=0.0_dp
     rds(min(diis%idsx,diis%ids)+1,ikpt)=1.0_dp
     
     !if(iproc==0)  write(6,*) 'DIIS matrix'
     !do i=1,min(diis%idsx,ids)+1
     !  if(iproc==0)  write(6,'(i3,12(1x,e9.2))') iproc,(ads(i,j,2),j=1,min(diis%idsx,ids)+1),rds(i)
     !enddo
     if (diis%ids > 1) then
        ! solve linear system:(LAPACK)
        call DSYSV('U',min(diis%idsx,diis%ids)+1,1,diis%ads(1,1,ikptp,2),diis%idsx+1,  & 
             ipiv,rds(1,ikpt),diis%idsx+1,diis%ads(1,1,ikptp,3),(diis%idsx+1)**2,info)
        
        if (info /= 0) then
           print*, 'diisstp: DSYSV',info
        end if
     else
        rds(1,ikpt)=1.0_dp
     endif

! new guess
     do iorb=1,orbs%norb
        call razero(nvctrp*orbs%nspinor,psit(ispsi+(iorb-1)*nvctrp*orbs%nspinor))
        
        jst=max(1,diis%ids-diis%idsx+1)
        jj=0
        do j=jst,diis%ids
           jj=jj+1
           mj=mod(j-1,diis%idsx)+1
           do k=1,nvctrp*orbs%nspinor
              psit(ispsi+(iorb-1)*nvctrp*orbs%nspinor+k-1)=&
                   psit(ispsi+(iorb-1)*nvctrp*orbs%nspinor+k-1)+&
                   rds(jj,ikpt)*(&
                   diis%psidst(ispsidst+k-1+(iorb-1)*nvctrp*orbs%nspinor+&
                   (mj-1)*orbs%norb*orbs%nspinor*nvctrp)&
                   -real(diis%hpsidst(ispsidst+k-1+(iorb-1)*nvctrp*orbs%nspinor+&
                   (mj-1)*orbs%norb*orbs%nspinor*nvctrp),wp))
           end do
        end do
     end do
     ispsi=ispsi+nvctrp*orbs%norb*orbs%nspinor
     ispsidst=ispsidst+nvctrp*orbs%norb*orbs%nspinor*diis%idsx
  end do
  ! Output to screen, depending on policy.
  if (verbose >= 10) then
     call broadcast_kpt_objects(nproc, orbs%nkpts, diis%idsx+1, rds, orbs%ikptproc)
  end if
  if (iproc == 0 .and. verbose > 0 .and. .not.quietHere) then 
     if (verbose < 10) then
        !we restrict the printing to the first k point only.
        write(*,'(1x,a,2x,12(1x,1pe9.2))')'DIIS weights',(rds(j,1),j=1,min(diis%idsx,diis%ids)+1)
     else
        do ikpt = 1, orbs%nkpts
           write(*,'(1x,a,I3.3,a,2x,12(1x,1pe9.2))')'DIIS weights (kpt #', ikpt, &
                & ')', (rds(j,0),j=1,min(diis%idsx,diis%ids)+1)
        end do
     end if
  endif
  
  i_all=-product(shape(ipiv))*kind(ipiv)
  deallocate(ipiv,stat=i_stat)
  call memocc(i_stat,i_all,'ipiv',subname)
  i_all=-product(shape(rds))*kind(rds)
  deallocate(rds,stat=i_stat)
  call memocc(i_stat,i_all,'rds',subname)

END SUBROUTINE diisstp


!> compute a dot product of two single precision vectors 
!! returning a double precision result
subroutine ds_dot(ndim,x,x0,dx,y,y0,dy,dot_out)
  use module_base
  use module_types
  implicit none
  integer, intent(in) :: ndim,x0,dx,y0,dy
  real(tp), dimension(ndim), intent(in) :: x,y
  real(wp), intent(out) :: dot_out
  integer :: jj,ix,iy
  
  dot_out=0.0d0
  ix=x0
  iy=y0
  do jj=1,ndim
    dot_out=dot_out + real(x(ix),wp)*real(y(iy),wp)
    ix=ix+dx
    iy=iy+dy
  end do
END SUBROUTINE ds_dot

function s2d_dot(ndim,x,dx,y,dy)
  implicit none
  integer, intent(in) :: ndim,dx,dy
  real(kind=4), dimension(ndim), intent(in) :: x,y
  real(kind=8) :: s2d_dot
  !local variables
  integer :: jj,ix,iy
  real(kind=8) :: dot_out

  !quick return if possible
  if (ndim <= 0) then
     s2d_dot=0.0d0
     return
  end if

  dot_out=0.0d0
  ix=1
  iy=1
  do jj=1,ndim
     dot_out=dot_out + real(x(ix),kind=8)*real(y(iy),kind=8)
     ix=ix+dx
     iy=iy+dy
  end do

  s2d_dot=dot_out

end function s2d_dot








!!****f* BigDFT/psimix
!! FUNCTION
!! COPYRIGHT
!!    Copyright (C) 2007-2010 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 
!! SOURCE
!!
subroutine psimixVariable(iproc,nproc,orbs,comms,diis,diisArr, hpsit,psit, quiet)
  use module_base
  use module_types
  use module_interfaces, excpet_this_one => psimixVariable
  implicit none
  integer, intent(in) :: iproc,nproc
  type(orbitals_data), intent(in) :: orbs
  type(communications_arrays), intent(in) :: comms
  type(diis_objects), intent(inout) :: diis
  type(diis_objects),dimension(orbs%norb),intent(in out):: diisArr
  real(wp), dimension(sum(comms%ncntt(0:nproc-1))), intent(inout) :: psit,hpsit
  logical, optional:: quiet ! to avoid that the DIIS weights are written
  !real(wp), dimension(:), pointer :: psit,hpsit
  !local variables
  integer :: ikptp,nvctrp,ispsi,ispsidst,jj,ierr


  if (diis%idsx > 0) then
     !do not transpose the hpsi wavefunction into the diis array
     !for compatibility with the k-points distribution
     ispsi=1
     ispsidst=1
     do ikptp=1,orbs%nkptsp
        nvctrp=comms%nvctr_par(iproc,ikptp)
        if (nvctrp == 0) cycle
        
     !here we can choose to store the DIIS arrays with single precision
     !psidst=psit
        call dcopy(nvctrp*orbs%norb*orbs%nspinor,&
             psit(ispsi),1,&
             diis%psidst(ispsidst+nvctrp*orbs%nspinor*orbs%norb*(diis%mids-1)),1)
     !hpsidst=hpsi
     !   call dcopy(nvctrp*orbs%norb*orbs%nspinor,&
     !        hpsit(ispsi),1,&
     !        hpsidst(ispsidst+nvctrp*orbs%nspinor*orbs%norb*(mids-1)),1)

     do jj=0,nvctrp*orbs%norb*orbs%nspinor-1
        diis%hpsidst(ispsidst+nvctrp*orbs%nspinor*orbs%norb*(diis%mids-1)+jj)&
             =real(hpsit(ispsi+jj),tp) !diis precision conversion
     end do
        ispsi=ispsi+nvctrp*orbs%norb*orbs%nspinor
        ispsidst=ispsidst+nvctrp*orbs%norb*orbs%nspinor*diis%idsx
     end do

     !here we should separate between up and down spin orbitals, maybe
     if(present(quiet)) then
         call diisstpVariable(iproc,nproc,orbs,comms,diis,diisArr,psit,quiet)
     else
         call diisstpVariable(iproc,nproc,orbs,comms,diis,diisArr,psit)
     end if

  else
     ! update all wavefunctions with the preconditioned gradient
     if (diis%energy > diis%energy_old) then
        diis%alpha=max(.125_wp,.5_wp*diis%alpha)
        if (diis%alpha == .125_wp) write(*,*) ' WARNING: Convergence problem or limit'
     else
        diis%alpha=min(1.05_wp*diis%alpha,1._wp)
     endif
     if (iproc == 0 .and. verbose > 0) write(*,'(1x,a,1pe11.3)') 'alpha=',diis%alpha

!!     do iorb=1,orbs%norb*orbs%nspinor
!!        call axpy(comms%nvctr_par(iproc),&
!!             -alpha,hpsi(1+comms%nvctr_par(iproc)*(iorb-1)),1,&
!!             psit(1+comms%nvctr_par(iproc)*(iorb-1)),1)
!!     enddo

     call axpy(sum(comms%ncntt(0:nproc-1)),-diis%alpha,hpsit(1),1,psit(1),1)

  endif

END SUBROUTINE psimixVariable
!!***



! diis subroutine:
! calculates the DIIS extrapolated solution psit in the ids-th DIIS step 
! using  the previous iteration points psidst and the associated error 
! vectors (preconditioned gradients) hpsidst
subroutine diisstpVariable(iproc,nproc,orbs,comms,diis,diisArr,psit,quiet)
  use module_base
  use module_types
  implicit none
! Arguments
  integer, intent(in) :: nproc,iproc
  type(orbitals_data), intent(in) :: orbs
  type(communications_arrays), intent(in) :: comms
  type(diis_objects), intent(inout) :: diis
  type(diis_objects),dimension(orbs%norb),intent(in out):: diisArr
  real(wp), dimension(sum(comms%ncntt(0:nproc-1))), intent(out) :: psit
  logical, optional:: quiet ! to avoid that the DIIS weights are written
! Local variables
  character(len=*), parameter :: subname='diisstp'
  integer :: i,j,ist,jst,mi,iorb,info,jj,mj,k,i_all,i_stat,ierr
  integer :: ikptp,ikpt,ispsi,ispsidst,nvctrp
  integer, dimension(:), allocatable :: ipiv
  real(dp), dimension(:,:), allocatable :: rds
  logical:: quietHere

  ! This decides whether the DIIS weights are written to the screen or not.
  ! If quietHere is false, then they are written, otherwise not.
  if(present(quiet) .and. quiet) then
      quietHere=.true.
  else
      quietHere=.false.
  end if

  allocate(ipiv(diisArr(1)%idsx+1+ndebug),stat=i_stat)
  call memocc(i_stat,ipiv,'ipiv',subname)
  allocate(rds(diisArr(1)%idsx+1,orbs%nkpts+ndebug),stat=i_stat)
  call memocc(i_stat,rds,'rds',subname)


  orbsLoop: do iorb=1,orbs%norb
      call razero((diisArr(iorb)%idsx+1)*orbs%nkpts,rds)

      ispsidst=1
      do ikptp=1,orbs%nkptsp
         ikpt=orbs%iskpts+ikptp!orbs%ikptsp(ikptp)
         nvctrp=comms%nvctr_par(iproc,ikptp)
         if (nvctrp == 0) cycle

         ! set up DIIS matrix (upper triangle)
         if (diisArr(iorb)%ids > diisArr(iorb)%idsx) then
            ! shift left up matrix
            do i=1,diisArr(iorb)%idsx-1
               do j=1,i
                  diisArr(iorb)%ads(j,i,ikptp,1)=diisArr(iorb)%ads(j+1,i+1,ikptp,1)
               end do
            end do
         end if

         ! calculate new line, use rds as work array for summation
         ist=max(1,diisArr(iorb)%ids-diisArr(iorb)%idsx+1)
         do i=ist,diisArr(iorb)%ids
            mi=mod(i-1,diisArr(iorb)%idsx)+1
!!         do iorb=1,norb*nspinor
!!            tt=dot(nvctrp,hpsidst(1,iorb,mids),1,hpsidst(1,iorb,mi),1)
!!            rds(i-ist+1)=rds(i-ist+1)+tt
!!         end do
            !to be corrected for complex wavefunctions
            rds(i-ist+1,ikpt)=dot(nvctrp*orbs%nspinor,&
                 diis%hpsidst(ispsidst+(iorb-1)*nvctrp*orbs%nspinor+(diisArr(iorb)%mids-1)*nvctrp*orbs%norb*orbs%nspinor),1,&
                 diis%hpsidst(ispsidst+(iorb-1)*nvctrp*orbs%nspinor+(mi-1)*nvctrp*orbs%norb*orbs%nspinor),1)
            !this has to be inserted in module_base
            !call ds_dot(nvctrp*orbs%norb*orbs%nspinor,&
            !     hpsidst_sp,ispsidst+(mids-1)*nvctrp*orbs%norb*orbs%nspinor,1,&
            !     hpsidst_sp,ispsidst+(mi  -1)*nvctrp*orbs%norb*orbs%nspinor,1,&
            !     rds(i-ist+1,ikpt))
         end do
         ispsidst=ispsidst+nvctrp*orbs%norb*orbs%nspinor*diis%idsx
      end do

      if (nproc > 1) then
         call mpiallred(rds(1,1),(diisArr(iorb)%idsx+1)*orbs%nkpts,MPI_SUM,MPI_COMM_WORLD,ierr)
!         call MPI_ALLREDUCE(MPI_IN_PLACE,rds,(diis%idsx+1)*orbs%nkpts,  & 
!                     mpidtypw,MPI_SUM,MPI_COMM_WORLD,ierr)
!
!!         call MPI_ALLREDUCE(rds,ads(1,min(diis%idsx,ids),1),min(ids,diis%idsx),  & 
!!                     MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
      endif
      
      ispsi=1
      ispsidst=1
      do ikptp=1,orbs%nkptsp
         ikpt=orbs%iskpts+ikptp!orbs%ikptsp(ikptp)
         nvctrp=comms%nvctr_par(iproc,ikptp)
         if (nvctrp == 0) cycle

         do i=1,min(diisArr(iorb)%ids,diisArr(iorb)%idsx)
            diisArr(iorb)%ads(i,min(diisArr(iorb)%idsx,diisArr(iorb)%ids),ikptp,1)=rds(i,ikpt)
         end do

         ! copy to work array, right hand side, boundary elements
         do j=1,min(diisArr(iorb)%idsx,diisArr(iorb)%ids)
            diisArr(iorb)%ads(j,min(diisArr(iorb)%idsx,diisArr(iorb)%ids)+1,ikptp,2)=1.0_wp
            rds(j,ikpt)=0.d0
            do i=j,min(diisArr(iorb)%idsx,diisArr(iorb)%ids)
               diisArr(iorb)%ads(j,i,ikptp,2)=diisArr(iorb)%ads(j,i,ikptp,1)
            end do
         end do
         diisArr(iorb)%ads(min(diisArr(iorb)%idsx,diisArr(iorb)%ids)+1,min(diisArr(iorb)%idsx,diisArr(iorb)%ids)+1,ikptp,2)=0.0_dp
         rds(min(diisArr(iorb)%idsx,diisArr(iorb)%ids)+1,ikpt)=1.0_dp
         
         !if(iproc==0)  write(6,*) 'DIIS matrix'
         !do i=1,min(diis%idsx,ids)+1
         !  if(iproc==0)  write(6,'(i3,12(1x,e9.2))') iproc,(ads(i,j,2),j=1,min(diis%idsx,ids)+1),rds(i)
         !enddo
         if (diisArr(iorb)%ids > 1) then
            ! solve linear system:(LAPACK)
            call DSYSV('U',min(diisArr(iorb)%idsx,diisArr(iorb)%ids)+1,1,diisArr(iorb)%ads(1,1,ikptp,2),diisArr(iorb)%idsx+1,  & 
                 ipiv,rds(1,ikpt),diisArr(iorb)%idsx+1,diisArr(iorb)%ads(1,1,ikptp,3),(diisArr(iorb)%idsx+1)**2,info)
            
            if (info /= 0) then
               print*, 'diisstp: DSYSV',info
            end if
         else
            rds(1,ikpt)=1.0_dp
         endif
     !end do orbsLoop

! new guess
     !do iorb=1,orbs%norb
         call razero(nvctrp*orbs%nspinor,psit(ispsi+(iorb-1)*nvctrp*orbs%nspinor))
         
         jst=max(1,diisArr(iorb)%ids-diisArr(iorb)%idsx+1)
         jj=0
         do j=jst,diisArr(iorb)%ids
            jj=jj+1
            mj=mod(j-1,diisArr(iorb)%idsx)+1
            do k=1,nvctrp*orbs%nspinor
               psit(ispsi+(iorb-1)*nvctrp*orbs%nspinor+k-1)=&
                    psit(ispsi+(iorb-1)*nvctrp*orbs%nspinor+k-1)+&
                    rds(jj,ikpt)*(&
                    diis%psidst(ispsidst+k-1+(iorb-1)*nvctrp*orbs%nspinor+&
                    (mj-1)*orbs%norb*orbs%nspinor*nvctrp)&
                    -real(diis%hpsidst(ispsidst+k-1+(iorb-1)*nvctrp*orbs%nspinor+&
                    (mj-1)*orbs%norb*orbs%nspinor*nvctrp),wp))
            end do
         end do
      !end do
      ispsi=ispsi+nvctrp*orbs%norb*orbs%nspinor
      ispsidst=ispsidst+nvctrp*orbs%norb*orbs%nspinor*diis%idsx
    end do
 end do orbsLoop
  ! Output to screen, depending on policy.
  if (verbose >= 10) then
     call broadcast_kpt_objects(nproc, orbs%nkpts, diis%idsx+1, rds, orbs%ikptproc)
  end if
  if (iproc == 0 .and. verbose > 0 .and. .not.quietHere) then 
     if (verbose < 10) then
        !we restrict the printing to the first k point only.
        write(*,'(1x,a,2x,12(1x,1pe9.2))')'DIIS weights',(rds(j,1),j=1,min(diis%idsx,diis%ids)+1)
     else
        do ikpt = 1, orbs%nkpts
           write(*,'(1x,a,I3.3,a,2x,12(1x,1pe9.2))')'DIIS weights (kpt #', ikpt, &
                & ')', (rds(j,0),j=1,min(diis%idsx,diis%ids)+1)
        end do
     end if
  endif
  
  i_all=-product(shape(ipiv))*kind(ipiv)
  deallocate(ipiv,stat=i_stat)
  call memocc(i_stat,i_all,'ipiv',subname)
  i_all=-product(shape(rds))*kind(rds)
  deallocate(rds,stat=i_stat)
  call memocc(i_stat,i_all,'rds',subname)

END SUBROUTINE diisstpVariable
