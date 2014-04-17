!> @file
!!   Routines to orthogonalize the wavefunctions
!! @author
!!    Copyright (C) 2007-2013 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 


!> Orthogonality routine, for all the orbitals
!! Uses wavefunctions in their transposed form
subroutine orthogonalize(iproc,nproc,orbs,comms,psi,orthpar,paw)
  use module_base
  use module_types
  use module_interfaces, except_this_one_A => orthogonalize
  use communications_base, only: comms_cubic
  implicit none
  integer, intent(in) :: iproc,nproc
  type(orbitals_data), intent(inout) :: orbs
  type(comms_cubic), intent(in) :: comms
  !>Choose which orthogonalization method shall be used:
  !!    orthpar%methOrtho==0: Cholesky orthonormalization (i.e. a pseudo Gram-Schmidt)
  !!    orthpar%methOrtho==1: hybrid Gram-Schmidt/Cholesky orthonormalization
  !!    orthpar%methOrtho==2: Loewdin orthonormalization
  type(orthon_data), intent(in) :: orthpar
  real(wp), dimension(comms%nvctr_par(iproc,0)*orbs%nspinor*orbs%norb), intent(inout) :: psi
  type(paw_objects),optional,intent(inout) :: paw
  !local variables
  character(len=*), parameter :: subname='orthogonalize'
  integer :: i_stat,i_all
  !integer :: i,idx
  integer :: ispin,nspin,nspinor,usepaw=0
  integer, dimension(:,:), allocatable :: ndim_ovrlp
  real(wp), dimension(:), allocatable :: ovrlp
  integer,dimension(:),allocatable:: norbArr
  character(len=20):: category
  !real(kind=8) :: test,test_max

  ! Determine whether we have close shell (nspin=1) or spin polarized (nspin=2)
  if (orbs%norbd>0) then
     nspin=2
  else 
     nspin=1
  end if

  !Determine whether we are in a paw calculation:
  if(present(paw)) usepaw=paw%usepaw

  ! ndim_ovrlp describes the shape of the overlap matrix.
  allocate(ndim_ovrlp(nspin,0:orbs%nkpts+ndebug),stat=i_stat)
  call memocc(i_stat,ndim_ovrlp,'ndim_ovrlp',subname)
  
  ! Allocate norbArr which contains the number of up and down orbitals.
  allocate(norbArr(nspin), stat=i_stat)
  call memocc(i_stat,norbArr,'norbArr',subname)
  do ispin=1,nspin
     if(ispin==1) norbArr(ispin)=orbs%norbu
     if(ispin==2) norbArr(ispin)=orbs%norbd
  end do

  ! Choose which orthogonalization method shall be used:
  ! methOrtho==0: Cholesky orthonormalization (i.e. a pseudo Gram-Schmidt)
  ! methOrtho==1: hybrid Gram-Schmidt/Cholesky orthonormalization
  ! methOrtho==2: Loewdin orthonormalization
  if(orthpar%methOrtho==0) then
     category='Chol'
     call timing(iproc, trim(category)//'_comput', 'ON')

     call dimension_ovrlp(nspin,orbs,ndim_ovrlp)

     ! Allocate the overlap matrix
     allocate(ovrlp(ndim_ovrlp(nspin,orbs%nkpts)+ndebug),stat=i_stat)
     call memocc(i_stat,ovrlp,'ovrlp',subname)

     ! Make a loop over npsin; calculate the overlap matrix (for up/down, resp.) and orthogonalize (again for up/down, resp.).
     do ispin=1,nspin
        if(usepaw==1) then
           call getOverlap_paw(iproc,nproc,nspin,norbArr(ispin),orbs,comms,&
                psi(1),paw%spsi(1),ndim_ovrlp,ovrlp,norbArr,1,ispin,category)

           call cholesky(iproc,nspin,norbArr(ispin),psi(1),orbs,comms,&
                ndim_ovrlp,ovrlp(1),norbArr,1,ispin,paw)
        else
           call getOverlap(iproc,nproc,nspin,norbArr(ispin),orbs,comms,&
                psi(1),ndim_ovrlp,ovrlp,norbArr,1,ispin,category)

           call cholesky(iproc,nspin,norbArr(ispin),psi(1),orbs,comms,&
                ndim_ovrlp,ovrlp(1),norbArr,1,ispin)
        end if
        !write(*,*)'orthogonality l80 erase me:'
        !write(*,*)'ovrlp',ovrlp
      end do

!             test = 0.d0 ; test_max = 0.d0 ; idx = 0.d0
!
!             do i = 1,orbs%norb**2
!
!                idx = i/orbs%norb
!
!
!                if((idx+1)-(i-idx*orbs%norb).eq.0) then
!                test = abs(ovrlp(i) - 1.d0)
!                else
!                test = abs(ovrlp(i) - 0.d0)
!                end if
!
!                !if(i.eq.1) test = abs(ovrlp(i)-1.d0)
!                if(i.eq.orbs%norb**2) test = abs(ovrlp(i)-1.d0)
!            
!
!                if(test.gt.test_max) test_max = test
!
!            end do
!
 !           write(*,*) 'Ovrlp-Difference', test_max


             ! Deallocate the arrays.
             i_all=-product(shape(ovrlp))*kind(ovrlp)
             deallocate(ovrlp,stat=i_stat)
             call memocc(i_stat,i_all,'ovrlp',subname)


          else if(orthpar%methOrtho==1) then
               category='GS/Chol'
               call timing(iproc, trim(category)//'_comput', 'ON')
               
               ! Make a hybrid Gram-Schmidt/Cholesky orthonormalization.
       if(usepaw==1) then
         call gsChol(iproc,nproc,psi(1),orthpar,nspinor,orbs,nspin,ndim_ovrlp,norbArr,comms,paw)
       else
         call gsChol(iproc,nproc,psi(1),orthpar,nspinor,orbs,nspin,ndim_ovrlp,norbArr,comms)
       end if

          else if(orthpar%methOrtho==2) then
             category='Loewdin'
             call timing(iproc,trim(category)//'_comput','ON')

             call dimension_ovrlp(nspin,orbs,ndim_ovrlp)
             
             ! Allocate the overlap matrix
             allocate(ovrlp(ndim_ovrlp(nspin,orbs%nkpts)+ndebug),stat=i_stat)
             call memocc(i_stat,ovrlp,'ovrlp',subname)
                  
             ! Make a loop over npsin; calculate the overlap matrix (for up/down,resp.) and orthogonalize (again for up/down,resp.).
             do ispin=1,nspin
        if(usepaw==1) then
           call getOverlap_paw(iproc,nproc,nspin,norbArr(ispin),orbs,comms,&
                psi(1),paw%spsi(1),ndim_ovrlp,ovrlp,norbArr,1,ispin,category)
           call loewdin(iproc,norbArr(ispin),orbs%nspinor,1,ispin,orbs,comms,&
                nspin,psi,ovrlp,ndim_ovrlp,norbArr,paw)
        else
           call getOverlap(iproc,nproc,nspin,norbArr(ispin),orbs,comms,psi(1),ndim_ovrlp,ovrlp,norbArr,1,ispin,category)
           call loewdin(iproc,norbArr(ispin),orbs%nspinor,1,ispin,orbs,comms,&
                nspin,psi,ovrlp,ndim_ovrlp,norbArr)
        end if
        !write(*,*)'orthogonality l117 erase me:'
        !write(*,*)'ovrlp',ovrlp
             end do
             
             ! Deallocate the arrays.
             i_all=-product(shape(ovrlp))*kind(ovrlp)
             deallocate(ovrlp,stat=i_stat)
             call memocc(i_stat,i_all,'ovrlp',subname)
                  
          else
             if(iproc==0) write(*,'(a,i0)') 'ERROR: invalid choice for methOrtho:',orthpar%methOrtho
             if(iproc==0) write(*,'(a)') "Change it in 'input.perf' to 0, 1 or 2!"
             stop
          end if

          ! Deallocate the remaining arrays.
          i_all=-product(shape(norbArr))*kind(norbArr)
          deallocate(norbArr, stat=i_stat)
          call memocc(i_stat,i_all,'norbArr',subname)

          i_all=-product(shape(ndim_ovrlp))*kind(ndim_ovrlp)
          deallocate(ndim_ovrlp, stat=i_stat)
          call memocc(i_stat,i_all,'ndim_ovrlp',subname)


          call timing(iproc,trim(category)//'_comput','OF')
          
        END SUBROUTINE orthogonalize


        subroutine check_closed_shell(orbs,lcs)
          use module_base
          use module_types
          implicit none
          type(orbitals_data), intent(in) :: orbs
          logical, intent(out) :: lcs
          !local variables
          integer :: iorb
          lcs=.true.
          do iorb=orbs%norb*orbs%nkpts,1,-1
             if ( orbs%occup(iorb) /= real(3-orbs%nspin,gp)) then
                lcs=.false.
                exit
             end if
          end do
        END SUBROUTINE check_closed_shell


!> Orthogonality constraint routine, for all the orbitals
!! Uses wavefunctions in their transposed form
subroutine orthoconstraint(iproc,nproc,orbs,comms,symm,psi,hpsi,scprsum,spsi) !n(c) wfd (arg:5)
  use module_base
  use module_types
  use module_interfaces, except_this_one => orthoconstraint
  use communications_base, only: comms_cubic
  implicit none
          logical, intent(in) :: symm !< symmetrize the lagrange multiplier after calculation
          integer, intent(in) :: iproc,nproc
          type(orbitals_data), intent(in) :: orbs
          type(comms_cubic), intent(in) :: comms
          !n(c) type(wavefunctions_descriptors), intent(in) :: wfd
          real(wp), dimension(orbs%npsidim_comp), intent(in) :: psi
          real(wp), dimension(orbs%npsidim_comp), optional, intent(in) :: spsi
          real(wp), dimension(orbs%npsidim_comp), intent(inout) :: hpsi
          real(dp), intent(out) :: scprsum
          !local variables
          character(len=*), parameter :: subname='orthoconstraint'
          integer :: i_stat,i_all,ierr,iorb,ialag,jorb !n(c) ise
          integer :: ispin,nspin,ikpt,norb,norbs,ncomp,nvctrp,ispsi,ikptp,nspinor
          real(dp) :: occ !n(c) tt
          real(gp), dimension(2) :: aij,aji
          integer, dimension(:,:), allocatable :: ndim_ovrlp
          real(wp), dimension(:), allocatable :: alag,paw_ovrlp

          !separate the orthogonalisation procedure for up and down orbitals 
          !and for different k-points
          call timing(iproc,'LagrM_comput  ','ON')

          !number of components of the overlap matrix for parallel case
          !calculate the dimension of the overlap matrix for each k-point
          if (orbs%norbd > 0) then
             nspin=2
          else
             nspin=1
          end if

          !number of components for the overlap matrix in wp-kind real numbers

          allocate(ndim_ovrlp(nspin,0:orbs%nkpts+ndebug),stat=i_stat)
          call memocc(i_stat,ndim_ovrlp,'ndim_ovrlp',subname)

          call dimension_ovrlp(nspin,orbs,ndim_ovrlp)

          allocate(alag(ndim_ovrlp(nspin,orbs%nkpts)+ndebug),stat=i_stat)
          call memocc(i_stat,alag,'alag',subname)
  
  !Allocate ovrlp for PAW: 
  if(present(spsi)) then
    norb=max(orbs%norbu,orbs%norbd,1)
    allocate(paw_ovrlp(norb+ndebug),stat=i_stat)
    call memocc(i_stat,paw_ovrlp,'paw_ovrlp',subname)
  end if

          !put to zero all the k-points which are not needed
          call to_zero(ndim_ovrlp(nspin,orbs%nkpts),alag)

          !do it for each of the k-points and separate also between up and down orbitals in the non-collinear case
          ispsi=1
          do ikptp=1,orbs%nkptsp
             ikpt=orbs%iskpts+ikptp!orbs%ikptsp(ikptp)

             do ispin=1,nspin

                call orbitals_and_components(iproc,ikpt,ispin,orbs,comms,&
                     nvctrp,norb,norbs,ncomp,nspinor)
                if (nvctrp == 0) cycle
        ialag=ndim_ovrlp(ispin,ikpt-1)+1

                if(nspinor==1) then
           if (symm) then
              !call gemmsy('T','N',norb,norb,nvctrp,1.0_wp,psi(ispsi),&
              call gemm('T','N',norb,norb,nvctrp,1.0_wp,psi(ispsi),&
                   max(1,nvctrp),hpsi(ispsi),max(1,nvctrp),0.0_wp,&
                   alag(ndim_ovrlp(ispin,ikpt-1)+1),norb)
           else
              call gemm('T','N',norb,norb,nvctrp,1.0_wp,psi(ispsi),&
                   ! TEMPORARYcall gemmsy('T','N',norb,norb,nvctrp,1.0_wp,psi(ispsi),&
                   max(1,nvctrp),hpsi(ispsi),max(1,nvctrp),0.0_wp,&
                   alag(ndim_ovrlp(ispin,ikpt-1)+1),norb)
           end if
        else
        !this part should be recheck in the case of nspinor == 2
        call c_gemm('C','N',norb,norb,ncomp*nvctrp,(1.0_wp,0.0_wp),psi(ispsi),&
             max(1,ncomp*nvctrp), &
             hpsi(ispsi),max(1,ncomp*nvctrp),(0.0_wp,0.0_wp),&
             alag(ndim_ovrlp(ispin,ikpt-1)+1),norb)
        end if
        !In PAW we should do <psi|H|psi>/<psi|S|psi>
        !However when the overlap is too large (usually when we have a bad initial guess)
        !Dividing by <psi|S|psi> is not a good idea since |gnrm> might get too large.
        !Hence the following part is not done:
        if (present(spsi) .and. 1==2) then
          if(nspinor==1) then
             !dgemmsy desactivated for the moment due to SIC
             !call gemmsy('T','N',norb,norb,nvctrp,1.0_wp,psi(ispsi),&
             call gemm('T','N',norb,norb,nvctrp,1.0_wp,psi(ispsi),&
                  max(1,nvctrp),psi(ispsi),max(1,nvctrp),0.0_wp,&
                  paw_ovrlp(1),norb)
                  !write(*,*)'orthoconstraint l260, erase me:'
                  !write(*,*)'<psi|psi>',paw_ovrlp
             call gemm('T','N',norb,norb,nvctrp,1.0_wp,psi(ispsi),&
                  max(1,nvctrp),spsi(ispsi),max(1,nvctrp),1.0_wp,&
                  paw_ovrlp(1),norb)
          else
             !this part should be recheck in the case of nspinor == 2
             call c_gemm('C','N',norb,norb,ncomp*nvctrp,(1.0_wp,0.0_wp),psi(ispsi),&
                  max(1,ncomp*nvctrp), &
                  psi(ispsi),max(1,ncomp*nvctrp),(0.0_wp,0.0_wp),&
                  paw_ovrlp(1),norb)
             call c_gemm('C','N',norb,norb,ncomp*nvctrp,(1.0_wp,0.0_wp),psi(ispsi),&
                  max(1,ncomp*nvctrp), &
                  spsi(ispsi),max(1,ncomp*nvctrp),(1.0_wp,0.0_wp),&
                  paw_ovrlp(1),norb)
          end if
          if(nproc>1) call mpiallred(paw_ovrlp(1),1,MPI_SUM,MPI_COMM_WORLD,ierr)
          alag(ialag:ialag+norb)=alag(ialag:ialag+norb)/paw_ovrlp(1:norb)
          !write(*,*)'orthoconstraint l268, erase me:'
          !write(*,*)'<psi|S|psi>',paw_ovrlp
        end if
        ispsi=ispsi+nvctrp*norb*nspinor
     end do
  end do

if (nproc > 1) then
  call timing(iproc,'LagrM_comput  ','OF')
  call timing(iproc,'LagrM_commun  ','ON')
  call mpiallred(alag(1),ndim_ovrlp(nspin,orbs%nkpts),MPI_SUM,bigdft_mpi%mpi_comm,ierr)
  call timing(iproc,'LagrM_commun  ','OF')
  call timing(iproc,'LagrM_comput  ','ON')
end if

!now each processors knows all the overlap matrices for each k-point
!even if it does not handle it.
!this is somehow redundant but it is one way of reducing the number of communications
!without defining group of processors

!calculate the sum of the diagonal of the overlap matrix, for each k-point
scprsum=0.0_dp
!for each k-point calculate the gradient
ispsi=1
do ikptp=1,orbs%nkptsp
  ikpt=orbs%iskpts+ikptp!orbs%ikptsp(ikptp)

  do ispin=1,nspin
     !n(c) if (ispin==1) ise=0
     call orbitals_and_components(iproc,ikpt,ispin,orbs,comms,&
          nvctrp,norb,norbs,ncomp,nspinor)
     if (nvctrp == 0) cycle

        !!$        !correct the orthogonality constraint if there are some orbitals which have zero occupation number
     if (symm) then
        do iorb=1,norb
           do jorb=iorb+1,norb
        !!$              if (orbs%occup((ikpt-1)*orbs%norb+iorb+ise) /= 0.0_gp .and. &
        !!$                   orbs%occup((ikpt-1)*orbs%norb+jorb+ise) == 0.0_gp) then
              aij(1)=alag(ndim_ovrlp(ispin,ikpt-1)+iorb+(jorb-1)*norbs)
              aji(1)=alag(ndim_ovrlp(ispin,ikpt-1)+jorb+(iorb-1)*norbs)
              aij(2)=0.0_gp
              aji(2)=0.0_gp
              alag(ndim_ovrlp(ispin,ikpt-1)+iorb+(jorb-1)*norbs) = 0.5_gp*(aij(1)+aji(1))!0.0_wp
              alag(ndim_ovrlp(ispin,ikpt-1)+jorb+(iorb-1)*norbs) = 0.5_gp*(aij(1)+aji(1))!0.0_wp
              if (norbs == 2*norb) then !imaginary part, if present
                 aij(2)=alag(ndim_ovrlp(ispin,ikpt-1)+iorb+(jorb-1)*norbs+1)
                 aji(2)=alag(ndim_ovrlp(ispin,ikpt-1)+jorb+(iorb-1)*norbs+1)
                 alag(ndim_ovrlp(ispin,ikpt-1)+iorb+(jorb-1)*norbs+1)=0.5_gp*(aij(2)-aji(2))
                 alag(ndim_ovrlp(ispin,ikpt-1)+jorb+(iorb-1)*norbs+1)=0.5_gp*(aji(2)-aij(2))
                 end if
        !!$                 !if (iproc ==0) print *,'i,j',iorb,jorb,alag(ndim_ovrlp(ispin,ikpt-1)+iorb+(jorb-1)*norbs)
        !!$              end if
              end do
           end do
        end if
        !!$        if (iproc ==0) print *,'matrix'
        !!$        do iorb=1,norb
        !!$           if (iproc ==0) print '(a,i3,100(1pe14.4))','i,j',iorb,&
        !!$                (alag(ndim_ovrlp(ispin,ikpt-1)+iorb+(jorb-1)*norbs),jorb=1,norb)
        !!$        end do

                !calculate the scprsum if the k-point is associated to this processor
                !the scprsum always coincide with the trace of the hamiltonian
                if (orbs%ikptproc(ikpt) == iproc) then
                   occ=real(orbs%kwgts(ikpt),dp)*real(3-orbs%nspin,gp)
                   if (nspinor == 4) occ=real(orbs%kwgts(ikpt),dp)
                   if(nspinor == 1) then
                      do iorb=1,norb
                         !write(*,'(a,2i5,2es14.7)') 'iproc, iorb, occ, lagMatVal', iproc, iorb, occ, real(alag(ndim_ovrlp(ispin,ikpt-1)+iorb+(iorb-1)*norbs),dp)
                         scprsum=scprsum+&
                              occ*real(alag(ndim_ovrlp(ispin,ikpt-1)+iorb+(iorb-1)*norbs),dp)
                      enddo
                   else if (nspinor == 4 .or. nspinor == 2) then
                      !not sure about the imaginary part of the diagonal (should be zero if H is hermitian)
                      do iorb=1,norb
                         scprsum=scprsum+&
                              occ*real(alag(ndim_ovrlp(ispin,ikpt-1)+2*iorb-1+(iorb-1)*norbs),dp)
                         scprsum=scprsum+&
                              occ*real(alag(ndim_ovrlp(ispin,ikpt-1)+2*iorb+(iorb-1)*norbs),dp)
                      enddo
                   end if
                end if
                !n(c) ise=norb

                if(nspinor==1 .and. nvctrp /= 0) then
                   call gemm('N','N',nvctrp,norb,norb,-1.0_wp,psi(ispsi),max(1,nvctrp),&
                        alag(ndim_ovrlp(ispin,ikpt-1)+1),norb,1.0_wp,&
                        hpsi(ispsi),max(1,nvctrp))
                else if (nvctrp /= 0) then
                   call c_gemm('N','N',ncomp*nvctrp,norb,norb,(-1.0_wp,0.0_wp),psi(ispsi),max(1,ncomp*nvctrp),&
                        alag(ndim_ovrlp(ispin,ikpt-1)+1),norb,(1.0_wp,0.0_wp),hpsi(ispsi),max(1,ncomp*nvctrp))
                end if

        !Only for PAW:
        if (present(spsi)) then
          if(nspinor==1 .and. nvctrp /= 0) then
             call gemm('N','N',nvctrp,norb,norb,-1.0_wp,spsi(ispsi),max(1,nvctrp),&
                  alag(ndim_ovrlp(ispin,ikpt-1)+1),norb,1.0_wp,&
                  hpsi(ispsi),max(1,nvctrp))
          else if (nvctrp /= 0) then
             call c_gemm('N','N',ncomp*nvctrp,norb,norb,(-1.0_wp,0.0_wp),spsi(ispsi),max(1,ncomp*nvctrp),&
                  alag(ndim_ovrlp(ispin,ikpt-1)+1),norb,(1.0_wp,0.0_wp),hpsi(ispsi),max(1,ncomp*nvctrp))
          end if
        end if

                ispsi=ispsi+nvctrp*norb*nspinor
             end do
          end do

          if (nproc > 1) then
             !n(c) tt=scprsum
             call mpiallred(scprsum,1,MPI_SUM,bigdft_mpi%mpi_comm,ierr)
             !call MPI_ALLREDUCE(tt,scprsum,1,mpidtypd,MPI_SUM,bigdft_mpi%mpi_comm,ierr)
          end if

          i_all=-product(shape(alag))*kind(alag)
          deallocate(alag,stat=i_stat)
          call memocc(i_stat,i_all,'alag',subname)

  if(present(spsi)) then
    i_all=-product(shape(paw_ovrlp))*kind(paw_ovrlp)
    deallocate(paw_ovrlp,stat=i_stat)
    call memocc(i_stat,i_all,'paw_ovrlp',subname)
  end if

  i_all=-product(shape(ndim_ovrlp))*kind(ndim_ovrlp)
  deallocate(ndim_ovrlp,stat=i_stat)
  call memocc(i_stat,i_all,'ndim_ovrlp',subname)

  call timing(iproc,'LagrM_comput  ','OF')

END SUBROUTINE orthoconstraint


!> Found the linear combination of the wavefunctions which diagonalises
!! the overlap matrix
subroutine subspace_diagonalisation(iproc,nproc,orbs,comms,psi,hpsi,evsum)
  use module_base
  use module_types
  use yaml_output
  use communications_base, only: comms_cubic
  implicit none
  integer, intent(in) :: iproc,nproc
  type(orbitals_data), intent(inout) :: orbs !eval is updated
  type(comms_cubic), intent(in) :: comms
  real(wp), dimension(comms%nvctr_par(iproc,0)*orbs%nspinor*orbs%norb), intent(in) :: hpsi
  real(wp), dimension(comms%nvctr_par(iproc,0)*orbs%nspinor*orbs%norb), intent(inout) :: psi
  real(wp), intent(out) :: evsum
  !local variables
  character(len=*), parameter :: subname='subspace_diagonalisation'
  integer :: i_stat,i_all,ierr,info,iorb,n_lp,n_rp,npsiw,isorb,ise,jorb,ncplx
  integer :: ispin,nspin,ikpt,norb,norbs,ncomp,nvctrp,ispsi,ikptp,nspinor
  !integer :: istart
  real(wp) :: occ,asymm
  real(gp), dimension(2) :: aij,aji
  integer, dimension(:,:), allocatable :: ndim_ovrlp
  real(wp), dimension(:), allocatable :: work_lp,work_rp,psiw
  real(wp), dimension(:), allocatable :: hamks

  !separate the diagonalisation procedure for up and down orbitals 
  !and for different k-points

!!$  !number of components of the overlap matrix for parallel case
!!$  istart=2
!!$  if (nproc == 1) istart=1

  !calculate the dimension of the overlap matrix for each k-point
  if (orbs%norbd > 0) then
     nspin=2
  else
     nspin=1
  end if

  !number of components for the overlap matrix in wp-kind real numbers
  allocate(ndim_ovrlp(nspin,0:orbs%nkpts+ndebug),stat=i_stat)
  call memocc(i_stat,ndim_ovrlp,'ndim_ovrlp',subname)

  call dimension_ovrlp(nspin,orbs,ndim_ovrlp)

  allocate(hamks(ndim_ovrlp(nspin,orbs%nkpts)+ndebug),stat=i_stat)
  call memocc(i_stat,hamks,'hamks',subname)

  !put to zero all the k-points which are not needed
  call to_zero(ndim_ovrlp(nspin,orbs%nkpts),hamks)

  !dimension of the work arrays
  n_lp=0
  n_rp=0
  npsiw=0

  !do it for each of the k-points and separate also between up and down orbitals in the non-collinear case
  ispsi=1
  do ikptp=1,orbs%nkptsp
     ikpt=orbs%iskpts+ikptp!orbs%ikptsp(ikptp)

     do ispin=1,nspin

        call orbitals_and_components(iproc,ikpt,ispin,orbs,comms,&
             nvctrp,norb,norbs,ncomp,nspinor)
        if (nvctrp == 0) cycle

        if(nspinor==1) then
           call gemm('T','N',norb,norb,nvctrp,1.0_wp,psi(ispsi),max(1,nvctrp),hpsi(ispsi),&
                max(1,nvctrp),0.0_wp,&
                hamks(ndim_ovrlp(ispin,ikpt-1)+1),norb)
        else
           !this part should be recheck in the case of nspinor == 2
           call c_gemm('C','N',norb,norb,ncomp*nvctrp,(1.0_wp,0.0_wp),psi(ispsi),&
                max(1,ncomp*nvctrp), &
                hpsi(ispsi),max(1,ncomp*nvctrp),(0.0_wp,0.0_wp),&
                hamks(ndim_ovrlp(ispin,ikpt-1)+1),norb)
        end if
        ispsi=ispsi+nvctrp*norb*nspinor

        !dimensions of the work arrays
        n_lp=max(4*norbs,1000,n_lp)
        n_rp=max(3*norb+1,n_rp)
        npsiw=max(nvctrp*orbs%norb*nspinor,npsiw)

     end do
  end do

  if (nproc > 1) then
     call mpiallred(hamks(1),ndim_ovrlp(nspin,orbs%nkpts),MPI_SUM,bigdft_mpi%mpi_comm,ierr)
  end if

  !now each processors knows all the overlap matrices for each k-point
  !even if it does not handle it.
  !this is somehow redundant but it is one way of reducing the number of communications
  !without defining group of processors

  allocate(work_lp(n_lp*2+ndebug),stat=i_stat)
  call memocc(i_stat,work_lp,'work_lp',subname)
  allocate(work_rp(n_rp+ndebug),stat=i_stat)
  call memocc(i_stat,work_rp,'work_rp',subname)

  allocate(psiw(npsiw+ndebug),stat=i_stat)
  call memocc(i_stat,psiw,'psiw',subname)

  !!if(iproc==0) then
  !!    ierr=0
  !!    do i_all=1,norb
  !!        do i_stat=1,norb
  !!            ierr=ierr+1
  !!            write(13000+iproc,*) i_all, i_stat, hamks(ierr,1)
  !!        end do
  !!    end do
  !!    write(13000+iproc,*) '=============================='
  !!end if

  !for each k-point now reorthogonalise wavefunctions
  !assume the hamiltonian is a hermitian matrix in the subspace.
  !Evaluate the non-symmetricity of the hamiltonian in the subspace
  if (iproc==0) then
     asymm=0.0_dp
     do ikpt=1,orbs%nkpts
        do ispin=1,nspin
           call orbitals_and_components(iproc,ikpt,ispin,orbs,comms,&
                nvctrp,norb,norbs,ncomp,nspinor)
           if (norbs == 2*norb) then
              ncplx=2
           else
              ncplx=1
           end if
           do iorb=1,norb
              do jorb=iorb+1,norb
                 aij(1)=hamks(ndim_ovrlp(ispin,ikpt-1)+1+(iorb-1)*ncplx+(jorb-1)*norbs)
                 aji(1)=hamks(ndim_ovrlp(ispin,ikpt-1)+1+(jorb-1)*ncplx+(iorb-1)*norbs)
                 aij(2)=0.0_gp
                 aji(2)=0.0_gp
                 if (norbs == 2*norb) then !imaginary part, if present
                    aij(2)=hamks(ndim_ovrlp(ispin,ikpt-1)+2+(iorb-1)*ncplx+(jorb-1)*norbs)
                    aji(2)=hamks(ndim_ovrlp(ispin,ikpt-1)+2+(jorb-1)*ncplx+(iorb-1)*norbs)
                 end if
                 asymm=max(asymm,(aij(1)-aji(1))**2+(aij(2)+aji(2))**2)
              end do
           end do
        end do
     end do
     call yaml_map('Non-Hermiticity of Hamiltonian in the Subspace',asymm,fmt='(1pe9.2)')
     if (asymm > 1.d-10) then
        call yaml_warning('KS Hamiltonian is not Hermitian in the subspace, diff:'//&
             trim(yaml_toa(asymm,fmt='(1pe9.2)')))
        if (verbose >= 3) then
           call yaml_open_sequence('KS Hamiltonian Matrix(ces)',advance='no')
           call yaml_comment('Rank of the matrix: '//adjustl(trim(yaml_toa(norb,fmt='(i6)'))))
           do ikpt=1,orbs%nkpts
              do ispin=1,nspin
                 if (orbs%nkpts > 1) then
                    call yaml_comment('Kpt No. '//adjustl(trim(yaml_toa(ikpt,fmt='(i6)')))//&
                         ', Spin No.'//adjustl(trim(yaml_toa(ispin,fmt='(i6)'))))
                 else
                    call yaml_comment('Spin No.'//adjustl(trim(yaml_toa(ispin,fmt='(i6)'))))
                 end if
                 call yaml_sequence(advance='no')
                 call orbitals_and_components(iproc,ikpt,ispin,orbs,comms,&
                      nvctrp,norb,norbs,ncomp,nspinor)
                 !call yaml_stream_attributes(indent=indentlevel)
                 call yaml_open_sequence(flow=.true.)
                 do iorb=1,norb
                    call yaml_sequence()
                    call yaml_open_sequence()
                    do jorb=1,norb
                       if (norbs == 2*norb) then
                          ncplx=2
                          call yaml_sequence(trim(yaml_toa(&
                               cmplx(&
                               hamks(ndim_ovrlp(ispin,ikpt-1)+1+(jorb-1)*ncplx+(iorb-1)*norbs),&
                               hamks(ndim_ovrlp(ispin,ikpt-1)+2+(jorb-1)*ncplx+(iorb-1)*norbs),kind=8),&
                               fmt='(1pe9.2)')))
                       else         
                          call yaml_sequence(trim(yaml_toa(hamks(ndim_ovrlp(ispin,ikpt-1)+jorb+(iorb-1)*norbs),&
                               fmt='(1pe9.2)')))
                       end if
                    end do
                    call yaml_close_sequence()
                    if (iorb < norb) call yaml_newline()
                 end do
                 call yaml_close_sequence()
              end do
           end do
           call yaml_close_sequence()
        end if
     end if
  end if
  ispsi=1
  evsum=0.0_wp
  do ikptp=1,orbs%nkptsp
     ikpt=orbs%iskpts+ikptp!orbs%ikptsp(ikptp)
     isorb=1
     do ispin=1,nspin

        call orbitals_and_components(iproc,ikpt,ispin,orbs,comms,&
             nvctrp,norb,norbs,ncomp,nspinor)
        if (nvctrp == 0) cycle

        if(nspinor==1) then

                   call syev('V','U',norb,hamks(ndim_ovrlp(ispin,ikpt-1)+1),norb,&
                orbs%eval(isorb+(ikpt-1)*orbs%norb),work_lp(1),n_lp,info)
           if (info /= 0) write(*,*) 'SYEV ERROR',info

        else

           call  heev('V','U',norb,hamks(ndim_ovrlp(ispin,ikpt-1)+1),norb,&
                orbs%eval(isorb+(ikpt-1)*orbs%norb),work_lp(1),n_lp,work_rp(1),info)
           if (info /= 0) write(*,*) 'HEEV ERROR',info

        end if

        !calculate the evsum if the k-point is associated to this processor
        if (orbs%ikptproc(ikpt) == iproc) then
           if (ispin==1) ise=0
           do iorb=1,norb
              occ=real(orbs%kwgts(ikpt)*orbs%occup((ikpt-1)*orbs%norb+iorb+ise),dp)
              evsum=evsum+orbs%eval(isorb+iorb-1+(ikpt-1)*orbs%norb)*occ
           enddo
           ise=norb
        end if

        ispsi=ispsi+nvctrp*norb*nspinor
        isorb=isorb+norb
     end do
  end do

        !!do iorb=1,norb
        !!   occ=real(orbs%kwgts(ikpt)*orbs%occup((ikpt-1)*orbs%norb+iorb),wp)
        !!   evsum=evsum+orbs%eval(isorb+iorb-1+(ikpt-1)*orbs%norb)*occ
        !!   !if (iproc.eq.0) write(*,'(1x,a,i0,a,1x,1pe21.14)') 'eval(',iorb,')=',eval(iorb)
        !!enddo

        !!        ! Transform to KS orbitals
        !!        ! dgemm can be used instead of daxpy
        !!        if(nspinor==1) then
        !!           do iorb=1,norb
        !!              call to_zero(nvctrp,psitt(1,iorb))
        !!              do jorb=1,norb
        !!                 alpha=hamks(jorb,iorb,1)
        !!                 call axpy(nvctrp,alpha,psit(1,jorb),1,psitt(1,iorb),1)
        !!              enddo
        !!           enddo
        !!        else
        !!           do iorb=1,norb
        !!              call to_zero(nvctrp*nspinor,psitt(1,iorb))
        !!              do jorb=1,norb
        !!                 call c_axpy(ncomp*nvctrp,hamks(2*jorb-1,iorb,1),psit(1,jorb),1,psitt(1,iorb),1)
        !!              enddo
        !!           enddo
        !!        end if

        !the matrix which is applied here is the passage matrix which can be exported outside of the routine

  if (nproc > 1) then
     do ikpt = 1, orbs%nkpts, 1
        do ispin=1,nspin
           call orbitals_and_components(iproc,ikpt,ispin,orbs,comms,&
                nvctrp,norb,norbs,ncomp,nspinor)
           if (iproc /= orbs%ikptproc(ikpt)) &
                & hamks(ndim_ovrlp(ispin,ikpt-1)+1:ndim_ovrlp(ispin,ikpt-1)+norbs*norb) = 0._wp
        end do
     end do
     call mpiallred(hamks(1),ndim_ovrlp(nspin,orbs%nkpts),MPI_SUM,bigdft_mpi%mpi_comm,ierr)
  end if

  ispsi=1
  do ikptp=1,orbs%nkptsp
     ikpt=orbs%iskpts+ikptp!orbs%ikptsp(ikptp)
     isorb=1
     do ispin=1,nspin

        call orbitals_and_components(iproc,ikpt,ispin,orbs,comms,&
             nvctrp,norb,norbs,ncomp,nspinor)
        if (nvctrp == 0) cycle

        !sample of dgemm
        if (nspinor == 1) then
           call gemm('N','N',nvctrp,norb,norb,1.0_wp,psi(ispsi),max(1,nvctrp),&
                hamks(ndim_ovrlp(ispin,ikpt-1)+1),norb,0.0_wp,psiw(1),max(1,nvctrp))
        else
           call c_gemm('N','N',ncomp*nvctrp,norb,norb,(1.0_wp,0.0_wp),&
                psi(ispsi),max(1,ncomp*nvctrp),hamks(ndim_ovrlp(ispin,ikpt-1)+1),norb,&
                (0.0_wp,0.0_wp),psiw(1),max(1,ncomp*nvctrp))
        end if

        call vcopy(nvctrp*norb*nspinor,psiw(1),1,psi(ispsi),1)

        !here we should add the same transformation for hpsi if required

        ispsi=ispsi+nvctrp*norb*nspinor
        isorb=isorb+norb
     end do
  end do

  if (nproc > 1) then
     !evsumtmp=evsum
     call mpiallred(evsum,1,MPI_SUM,bigdft_mpi%mpi_comm,ierr)
  end if

  i_all=-product(shape(psiw))*kind(psiw)
  deallocate(psiw,stat=i_stat)
  call memocc(i_stat,i_all,'psiw',subname)

  i_all=-product(shape(work_lp))*kind(work_lp)
  deallocate(work_lp,stat=i_stat)
  call memocc(i_stat,i_all,'work_lp',subname)
  i_all=-product(shape(work_rp))*kind(work_rp)
  deallocate(work_rp,stat=i_stat)
  call memocc(i_stat,i_all,'work_rp',subname)

  i_all=-product(shape(hamks))*kind(hamks)
  deallocate(hamks,stat=i_stat)
  call memocc(i_stat,i_all,'hamks',subname)

  i_all=-product(shape(ndim_ovrlp))*kind(ndim_ovrlp)
  deallocate(ndim_ovrlp,stat=i_stat)
  call memocc(i_stat,i_all,'ndim_ovrlp',subname)

END SUBROUTINE subspace_diagonalisation


!> Makes sure all psivirt/gradients are othogonal to the occupied states psi.
!! This routine is almost the same as orthoconstraint_p. Difference:
!! hpsi(:,norb) -->  psivirt(:,nvirte) , therefore rectangular alag.
!! @warning
!!   Orthogonality to spin polarized channels is achieved in two calls,
subroutine orthon_virt_occup(iproc,nproc,orbs,orbsv,comms,commsv,psi_occ,psi_virt,msg)
  use module_base
  use module_types
          use communications_base, only: comms_cubic
  implicit none
  logical, intent(in) :: msg
  integer, intent(in) :: iproc,nproc
  type(orbitals_data), intent(in) :: orbs,orbsv
  type(comms_cubic), intent(in) :: comms,commsv
  real(wp), dimension(comms%nvctr_par(iproc,0)*orbs%nspinor*orbs%norb), intent(in) :: psi_occ
  real(wp), dimension(commsv%nvctr_par(iproc,0)*orbsv%nspinor*orbsv%norb), intent(inout) :: psi_virt
  !local variables
  character(len=*), parameter :: subname='orthon_virt_occup'
  integer :: i_stat,i_all,ierr,ispsiv,iorb,jorb,isorb
  integer :: ispin,nspin,ikpt,norb,norbs,ncomp,nvctrp,ispsi,ikptp,nspinor
  integer :: norbv,norbsv,ncompv,nvctrpv,nspinorv
  real(wp) :: scprsum,tt
  integer, dimension(:,:), allocatable :: ndim_ovrlp
  real(wp), dimension(:), allocatable :: alag

  !separate the orthogonalisation procedure for up and down orbitals 
  !and for different k-points
  call timing(iproc,'LagrM_comput  ','ON')

  !calculate the dimension of the overlap matrix for each k-point
  if (orbs%norbd > 0) then
     nspin=2
  else
     nspin=1
  end if

  !number of components for the overlap matrix in wp-kind real numbers

  allocate(ndim_ovrlp(nspin,0:orbs%nkpts+ndebug),stat=i_stat)
  call memocc(i_stat,ndim_ovrlp,'ndim_ovrlp',subname)

  call dimension_ovrlp_virt(nspin,orbs,orbsv,ndim_ovrlp)

  allocate(alag(ndim_ovrlp(nspin,orbs%nkpts)+ndebug),stat=i_stat)
  call memocc(i_stat,alag,'alag',subname)

  !put to zero all the k-points which are not needed
  call to_zero(ndim_ovrlp(nspin,orbs%nkpts),alag)

  !differentiate between real and complex wavefunctions
  !Lower triangle of overlap matrix using BLAS
  !     ovrlp(iorb,jorb)=psit(k,iorb)*psit(k,jorb) ; lower triangle


          !do it for each of the k-points and separate also between up and down orbitals in the non-collinear case
          ispsi=1
          ispsiv=1
          do ikptp=1,orbs%nkptsp
             ikpt=orbs%iskpts+ikptp!orbs%ikptsp(ikptp)

             do ispin=1,nspin

                call orbitals_and_components(iproc,ikpt,ispin,orbs,comms,&
                     nvctrp,norb,norbs,ncomp,nspinor)
                call orbitals_and_components(iproc,ikpt,ispin,orbsv,commsv,&
                     nvctrpv,norbv,norbsv,ncompv,nspinorv)
                !there checks ensure that the component distribution scheme of virtual and occupied states is the same
                if (nvctrpv /= nvctrp) stop 'nvctrp'
                if (ncompv /= ncomp) stop 'ncomp'
                if (nspinorv /= nspinor) stop 'nspinor'

                if (nvctrp == 0) cycle
                
                norbv=orbsv%norbu
                if (ispin==2) norbv=orbsv%norbd

                !print *,'nvctrp',iproc,commsv%nvctr_par(iproc,ikptp),nvctrp,ikpt,orbs%nkpts,orbsv%nkpts,norbs,norbv,orbsv%nkptsp,orbs%nkptsp

                !print *,'iproc,nvctrp,nspin,norb,ispsi,ndim_ovrlp',iproc,nvctrp,nspin,norb,ispsi,ndim_ovrlp(ispin,ikpt-1)

                if(nspinor==1) then
                   call gemm('T','N',norb,norbv,nvctrp,1.0_wp,psi_occ(ispsi),max(1,nvctrp),&
                        psi_virt(ispsiv),max(1,nvctrp),0.0_wp,&
                        alag(ndim_ovrlp(ispin,ikpt-1)+1),norb)
                else
                   !this part should be recheck in the case of nspinor == 2
                   call c_gemm('C','N',norb,norbv,ncomp*nvctrp,(1.0_wp,0.0_wp),psi_occ(ispsi),&
                        max(1,ncomp*nvctrp), &
                        psi_virt(ispsiv),max(1,ncomp*nvctrp),(0.0_wp,0.0_wp),&
                        alag(ndim_ovrlp(ispin,ikpt-1)+1),norb)
                end if

                ispsi=ispsi+nvctrp*norb*nspinor
                ispsiv=ispsiv+nvctrp*norbv*nspinor
             end do
          end do

          if (nproc > 1) then
             call timing(iproc,'LagrM_comput  ','OF')
             call timing(iproc,'LagrM_commun  ','ON')
             call mpiallred(alag(1),ndim_ovrlp(nspin,orbs%nkpts),MPI_SUM,bigdft_mpi%mpi_comm,ierr)
             !call MPI_ALLREDUCE (alag(1,2),alag(1,1),ndim_ovrlp(nspin,orbs%nkpts),&
             !     mpidtypw,MPI_SUM,bigdft_mpi%mpi_comm,ierr)
             call timing(iproc,'LagrM_commun  ','OF')
             call timing(iproc,'LagrM_comput  ','ON')
          end if

          !now each processors knows all the overlap matrices for each k-point
          !even if it does not handle it.
          !this is somehow redundant but it is one way of reducing the number of communications
          !without defining group of processors

          !for each k-point now reorthogonalise wavefunctions
          ispsi=1
          ispsiv=1
          isorb=1
          do ikptp=1,orbs%nkptsp
             ikpt=orbs%iskpts+ikptp!orbs%ikptsp(ikptp)

             do ispin=1,nspin

                call orbitals_and_components(iproc,ikpt,ispin,orbs,comms,&
                     nvctrp,norb,norbs,ncomp,nspinor)
                if (nvctrp == 0) cycle

                norbv=orbsv%norbu
                if (ispin==2) norbv=orbsv%norbd


                if (msg .and. .false.) then
                   write(*,'(1x,a)')'scalar products are'
                   write(*,'(1x,a)')'iocc  ivirt       value'!               zero if<1d-12'

                   scprsum=0.0_wp

                   do iorb=1,norb
                      do jorb=1,norbv
                         tt=alag(isorb+iorb-1+(jorb-1)*norbs)
                         write(*,'(1x,2i3,1pe21.14)')iorb,jorb,tt
                         scprsum=scprsum+tt**2
                         !if(abs(tt)<1d-12)alag(iorb,jorb,1)=0d0
                         !if(msg)write(*,'(2(i3),7x,2(1pe21.14))')iorb,jorb,tt,alag(iorb,jorb,1)
                      end do
                   enddo
                   scprsum=sqrt(scprsum/real(norb,wp)/real(norbv,wp))
                   write(*,'(1x,a,1pe21.14)')'sqrt sum squares is',scprsum
                   write(*,'(1x)')
                end if

                if(nspinor==1 .and. nvctrp /= 0) then
                   call gemm('N','N',nvctrp,norbv,norb,-1.0_wp,psi_occ(ispsi),max(1,nvctrp),&
                        alag(ndim_ovrlp(ispin,ikpt-1)+1),norb,1.0_wp,&
                        psi_virt(ispsiv),max(1,nvctrp))
                else if (nvctrp /= 0) then
                   call c_gemm('N','N',ncomp*nvctrp,norbv,norb,(-1.0_wp,0.0_wp),psi_occ(ispsi),&
                        max(1,ncomp*nvctrp),&
                        alag(ndim_ovrlp(ispin,ikpt-1)+1),norb,(1.0_wp,0.0_wp),psi_virt(ispsi),&
                        max(1,ncomp*nvctrp))
                end if
                ispsi=ispsi+nvctrp*norb*nspinor
                ispsiv=ispsiv+nvctrp*norbv*nspinor
                isorb=isorb+norbs*norbv
             end do
          end do


          i_all=-product(shape(alag))*kind(alag)
          deallocate(alag,stat=i_stat)
          call memocc(i_stat,i_all,'alag',subname)

          i_all=-product(shape(ndim_ovrlp))*kind(ndim_ovrlp)
          deallocate(ndim_ovrlp,stat=i_stat)
          call memocc(i_stat,i_all,'ndim_ovrlp',subname)

          call timing(iproc,'LagrM_comput  ','OF')

        END SUBROUTINE orthon_virt_occup


        subroutine complex_components(nspinor,norb,norbs,ncomp)
          implicit none
          integer, intent(in) :: nspinor,norb
          integer, intent(out) :: norbs,ncomp

          if(nspinor == 1) then
             norbs=norb
             ncomp=1 !useless
          else if (nspinor == 2) then
             norbs=2*norb
             ncomp=1
          else if (nspinor == 4) then
             norbs=2*norb
             ncomp=2
          end if
          
        END SUBROUTINE complex_components


        subroutine orbitals_and_components(iproc,ikpt,ispin,orbs,comms,nvctrp,norb,norbs,ncomp,nspinor)
          use module_base
          use module_types
          use communications_base, only: comms_cubic
          implicit none
          integer, intent(in) :: iproc,ikpt,ispin
          type(orbitals_data), intent(in) :: orbs
          type(comms_cubic), intent(in) :: comms
          integer, intent(out) :: nvctrp,norb,norbs,ncomp,nspinor
          
          nvctrp=comms%nvctr_par(iproc,ikpt)
          norb=orbs%norbu
          nspinor=orbs%nspinor
          if (ispin==2) norb=orbs%norbd

          call complex_components(nspinor,norb,norbs,ncomp)

        END SUBROUTINE orbitals_and_components


        subroutine dimension_ovrlp(nspin,orbs,ndim_ovrlp)
          use module_base
          use module_types
          implicit none
          integer, intent(in) :: nspin
          type(orbitals_data), intent(in) :: orbs
          integer, dimension(nspin,0:orbs%nkpts), intent(out) :: ndim_ovrlp
          !local variables
          integer :: norb,norbs,ncomp,ikpt

          ndim_ovrlp(1,0)=0
          if (nspin == 2) then
             norb=orbs%norbu

             !this is first k-point
             call complex_components(orbs%nspinor,norb,norbs,ncomp)

             ndim_ovrlp(2,0)=norbs*norb
          end if

          do ikpt=1,orbs%nkpts
             !this part should be enhanced for real k-points
             norb=orbs%norbu
             if (nspin == 2) norb = orbs%norbd
             !this is ikpt k-point
             call complex_components(orbs%nspinor,norb,norbs,ncomp)

             ndim_ovrlp(1,ikpt)=ndim_ovrlp(nspin,ikpt-1)+norbs*norb
             if (orbs%norbd > 0) then
                norb=orbs%norbu
                !this is ikpt+1
                call complex_components(orbs%nspinor,norb,norbs,ncomp)
                if (ikpt == orbs%nkpts) then
                   ndim_ovrlp(2,ikpt)=ndim_ovrlp(1,ikpt)
                else
                   ndim_ovrlp(2,ikpt)=ndim_ovrlp(1,ikpt)+norbs*norb
                end if
             end if
          end do

END SUBROUTINE dimension_ovrlp


subroutine dimension_ovrlp_virt(nspin,orbs,orbsv,ndim_ovrlp)
  use module_base
  use module_types
  implicit none
  integer, intent(in) :: nspin
  type(orbitals_data), intent(in) :: orbs,orbsv
  integer, dimension(nspin,0:orbs%nkpts), intent(out) :: ndim_ovrlp
  !local variables
  integer :: norb,norbs,ncomp,ikpt,norbv

  ndim_ovrlp(1,0)=0
  if (nspin == 2) then
     norb=orbs%norbu
     norbv=orbsv%norbu
     call complex_components(orbs%nspinor,norb,norbs,ncomp)

     ndim_ovrlp(2,0)=norbs*norbv
  end if

  do ikpt=1,orbs%nkpts
     !this part should be enhanced for real k-points
     norb=orbs%norbu
     norbv=orbsv%norbu 
     if (nspin == 2) then
        norb=orbs%norbd
        norbv=orbsv%norbd
     end if


     call complex_components(orbs%nspinor,norb,norbs,ncomp)

     ndim_ovrlp(1,ikpt)=ndim_ovrlp(nspin,ikpt-1)+norbs*norbv
     if (orbs%norbd > 0) then

        norb=orbs%norbu
        norbv=orbsv%norbu
        
        call complex_components(orbs%nspinor,norb,norbs,ncomp)

        if (ikpt == orbs%nkpts) then
           ndim_ovrlp(2,ikpt)=ndim_ovrlp(1,ikpt)
        else
           ndim_ovrlp(2,ikpt)=ndim_ovrlp(1,ikpt)+norbs*norbv
        end if
     end if
  end do

END SUBROUTINE dimension_ovrlp_virt


!> Effect of orthogonality constraints on gradient 
subroutine orthoconstraint_p(iproc,nproc,norb,occup,nvctrp,psit,hpsit,scprsum,nspinor)
  use module_base
  use module_types
  implicit none
  integer, intent(in) :: iproc,nproc,norb,nvctrp,nspinor
  real(gp), dimension(norb), intent(in) :: occup
  real(wp), dimension(nspinor*nvctrp,norb), intent(in) :: psit 
  real(dp), intent(out) :: scprsum
  real(wp), dimension(nspinor*nvctrp,norb), intent(out) :: hpsit
  !local variables
  character(len=*), parameter :: subname='orthoconstraint_p'
  integer :: i_stat,i_all,istart,iorb,ierr,norbs,ncomp
  real(dp) :: occ
  real(wp), dimension(:,:,:), allocatable :: alag

  call timing(iproc,'LagrM_comput  ','ON')

  istart=2
  if (nproc == 1) istart=1

  if(nspinor == 1) then
     norbs=norb
  else if(nspinor == 2) then
     norbs=2*norb
     ncomp=1
  else if (nspinor == 4) then
     norbs=2*norb
     ncomp=2
  end if

  allocate(alag(norbs,norb,istart+ndebug),stat=i_stat)
  call memocc(i_stat,alag,'alag',subname)

  !initialise if nvctrp=0
  if (nvctrp == 0) then
     call to_zero(norbs*norb*istart,alag(1,1,1))
  end if

  !     alag(jorb,iorb,istart)=+psit(k,jorb)*hpsit(k,iorb)
  if(nspinor==1) then
     call GEMM('T','N',norb,norb,nvctrp,1.0_wp,psit(1,1),max(1,nvctrp),hpsit(1,1),max(1,nvctrp),0.0_wp,&
          alag(1,1,istart),norb)
  else
     !this part should be recheck in the case of nspinor == 2
     call C_GEMM('C','N',norb,norb,ncomp*nvctrp,(1.0_wp,0.0_wp),psit(1,1),max(1,ncomp*nvctrp), &
          hpsit(1,1),max(1,ncomp*nvctrp),(0.0_wp,0.0_wp),alag(1,1,istart),norb)
  end if

  if (nproc > 1) then
     call timing(iproc,'LagrM_comput  ','OF')
     call timing(iproc,'LagrM_commun  ','ON')
     call MPI_ALLREDUCE(alag(1,1,2),alag(1,1,1),norbs*norb,&
          mpidtypd,MPI_SUM,bigdft_mpi%mpi_comm,ierr)
     call timing(iproc,'LagrM_commun  ','OF')
     call timing(iproc,'LagrM_comput  ','ON')
  end if
!          if (iproc.eq.0) then
!          write(*,*) 'ALAG',iproc,norb,norbs
!          do iorb=1,norb
!          write(*,'(10(1x,1pe10.3))') (alag(jorb,iorb,1),jorb=1,norbs)
!          enddo
!          endif



  scprsum=0.0_dp
  if(nspinor == 1) then
     do iorb=1,norb
        occ=real(occup(iorb),dp)
        scprsum=scprsum+occ*real(alag(iorb,iorb,1),dp)
     enddo
  else if (nspinor == 4 .or. nspinor == 2) then
     !not sure about the imaginary part of the diagonal
    do iorb=1,norb
       occ=real(occup(iorb),dp)
       scprsum=scprsum+occ*real(alag(2*iorb-1,iorb,1),dp)
       scprsum=scprsum+occ*real(alag(2*iorb,iorb,1),dp)
     enddo
  end if


!  if(iproc==0) print *,'ortho_p',scprsum

  ! hpsit(k,iorb)=-psit(k,jorb)*alag(jorb,iorb,1)
  if(nspinor==1) then
     call GEMM('N','N',nvctrp,norb,norb,-1.0_wp,psit(1,1),max(1,nvctrp),alag(1,1,1),norb,1.0_wp,&
          hpsit(1,1),max(1,nvctrp))
  else
     call C_GEMM('N','N',ncomp*nvctrp,norb,norb,(-1.0_wp,0.0_wp),psit(1,1),max(1,ncomp*nvctrp),&
          alag(1,1,1),norb,(1.0_wp,0.0_wp),hpsit(1,1),max(1,ncomp*nvctrp))
  end if


  i_all=-product(shape(alag))*kind(alag)
  deallocate(alag,stat=i_stat)
  call memocc(i_stat,i_all,'alag',subname)



  call timing(iproc,'LagrM_comput  ','OF')

END SUBROUTINE orthoconstraint_p


!> Gram-Schmidt orthogonalisation
subroutine orthon_p(iproc,nproc,norb,nvctrp,psit,nspinor)
  use module_base
  use module_types
  implicit none
  integer, intent(in) :: iproc,nproc,norb,nvctrp,nspinor
  real(wp), dimension(nspinor*nvctrp,norb), intent(inout) :: psit
  !local variables
  character(len=*), parameter :: subname='orthon_p'
  integer :: info,i_all,i_stat,nvctr_eff,ierr,istart,norbs,ncomp
  real(wp) :: tt,ttLOC
  real(wp), dimension(:,:,:), allocatable :: ovrlp
  integer :: volta

  call timing(iproc,'GramS_comput  ','ON')

  do volta=1,2

  if (norb == 1) then 

     !for the inhomogeneous distribution this should  be changed
     nvctr_eff=nvctrp!min(nvctr_tot-iproc*nvctrp,nvctrp)

     if (nvctr_eff > 0) then
     !parallel treatment of a run with only one orbital
        if(nspinor==1) then
           tt=nrm2(nvctr_eff,psit(1,1),1)
        else
           !print *,'for one orbital the norm of the spinor must be calculated'
           !stop
           tt=nrm2(nvctr_eff*nspinor,psit(1,1),1) !NOT CORRECT
        end if
        ttLOC=tt**2
     else
        ttLOC=0.0_wp
     end if
     
     if (nproc > 1) then
        call MPI_ALLREDUCE(ttLOC,tt,1,mpidtypd,MPI_SUM,bigdft_mpi%mpi_comm,ierr)
     else
        tt=ttLOC
     end if

     tt=1.0_wp/sqrt(tt)
     if(nspinor==1) then 
        !correct normalisation
        call vscal(nvctr_eff,tt,psit(1,1),1)
     else
        !not correct, to be adjusted
        call vscal(nvctr_eff*nspinor,tt,psit(1,1),1)
     end if

  else

     istart=2
     if (nproc == 1) istart=1

     if(nspinor==1) then
        norbs=norb
     else if (nspinor ==2) then
        norbs=2*norb
        ncomp=1
     else if (nspinor ==4) then
        norbs=2*norb
        ncomp=2
     end if

     allocate(ovrlp(norbs,norb,istart+ndebug),stat=i_stat)
     call memocc(i_stat,ovrlp,'ovrlp',subname)

     call to_zero(norbs*norb*istart,ovrlp(1,1,1))

     ! Upper triangle of overlap matrix using BLAS
     !     ovrlp(iorb,jorb)=psit(k,iorb)*psit(k,jorb) ; upper triangle
     if(nspinor==1) then
        call syrk('L','T',norb,nvctrp,1.0_wp,psit(1,1),max(1,nvctrp),0.0_wp,ovrlp(1,1,istart),norb)
     else
!!        ovrlp=0.0d0
!!        do iorb=1,norb
!!           do jorb=1,norb
!!              ttr=ddot(nvctrp*nspinor,psit(1,iorb),1,psit(1,jorb),1)
!!              tti=ddot(nvctrp,psit(1,iorb),1,psit(nvctrp+1,jorb),1)
!!              tti=tti-ddot(nvctrp,psit(nvctrp+1,iorb),1,psit(1,jorb),1)
!!              tti=tti-ddot(nvctrp,psit(3*nvctrp+1,iorb),1,psit(2*nvctrp+1,jorb),1)
!!              tti=tti+ddot(nvctrp,psit(2*nvctrp+1,iorb),1,psit(3*nvctrp+1,jorb),1)
!!              ovrlp(2*iorb-1,jorb,1)=ttr
!!              ovrlp(2*iorb,jorb,1)=tti*0.0d0
!!              print *,iorb,norb,ttr,tti
!!           end do
!!        end do
!!        stop
        call herk('L','C',norb,ncomp*nvctrp,1.0_wp,psit(1,1),max(1,ncomp*nvctrp),&
             0.0_wp,ovrlp(1,1,istart),norb)
     end if

     if (nproc > 1) then
        call timing(iproc,'GramS_comput  ','OF')
        call timing(iproc,'GramS_commun  ','ON')
        call MPI_ALLREDUCE (ovrlp(1,1,2),ovrlp(1,1,1),norbs*norb,&
             MPI_DOUBLE_PRECISION,MPI_SUM,bigdft_mpi%mpi_comm,ierr)
        call timing(iproc,'GramS_commun  ','OF')
        call timing(iproc,'GramS_comput  ','ON')
     end if
     
!!     if (iproc==0) then
!!        write(*,*) 'parallel ovrlp'
!!        do i=1,norbs,min(2,nspinor)
!!           write(*,'(10(1x,1pe10.3))') (ovrlp(i,j,1),j=1,norb)
!!        enddo
!!     end if

     !to be excluded if nvctrp==0
     if(nspinor==1) then
        
        ! Cholesky factorization
        call potrf( 'L',norb,ovrlp(1,1,1),norb,info)
        if (info /= 0) then
           write(*,*) 'info Cholesky factorization',info
        end if
        
        ! calculate L^{-1}
        call trtri( 'L','N',norb,ovrlp(1,1,1),norb,info)
        if (info.ne.0) write(6,*) 'info L^-1',info
        
        ! new vectors   
        call trmm ('R','L','T','N',nvctrp,norb,1.0_wp,ovrlp(1,1,1),norb,psit(1,1),max(1,nvctrp))

     else

       ! Cholesky factorization
!!        do i=1,norb
!!           if(iproc==0) then
!!              write(*,*) 'parallel ovrlp',i
!!              write(*,'(10f10.3)') (ovrlp(j,i,1), j=1,norbs)
!!           end if
!!        end do
        call c_potrf( 'L',norb,ovrlp(1,1,1),norb,info )
        if (info /= 0) then
           write(*,*) 'info Cholesky factorization',info
        end if
        
        ! calculate L^{-1}
!!         do i=1,norb
!!           if(iproc==0) then
!!              write(*,*) 'parallel ovrlp2',i
!!              write(*,'(10f10.3)') (ovrlp(j,i,1), j=1,norbs)
!!           end if
!!        end do
       call c_trtri( 'L','N',norb,ovrlp(1,1,1),norb,info)
        if (info.ne.0) write(6,*) 'info L^-1',info
        
!!        do i=1,norb
!!           if(iproc==0) then
!!              write(*,'(10f10.3)') (ovrlp(j,i,1), j=1,norbs)
!!           end if
!!        end do
       ! new vectors   !!check if third argument should be transpose or conjugate
        call c_trmm ('R','L','C','N',ncomp*nvctrp,norb,(1.0_wp,0.0_wp),&
             ovrlp(1,1,1),norb,psit(1,1),max(1,ncomp*nvctrp))

        !if(nproc==1) call psitransspi(nvctrp,norb,psit,.true.)


     end if

     i_all=-product(shape(ovrlp))*kind(ovrlp)
     deallocate(ovrlp,stat=i_stat)
     call memocc(i_stat,i_all,'ovrlp',subname)

  end if

  enddo

  call timing(iproc,'GramS_comput  ','OF')

END SUBROUTINE orthon_p


!> Loewdin orthogonalisation
!! @todo
!!  The loewe routines must be uniformised serial/parallel and nspinor should be added
subroutine loewe_p(iproc,norb,ndim,nvctrp,nvctr_tot,psit)
  use module_base
  use module_types
  implicit none
  !Argumments
  integer, intent(in) :: iproc,norb,ndim,nvctrp,nvctr_tot
  real(dp), dimension(nvctrp,ndim), intent(inout) :: psit
  !Local variables
  character(len=*), parameter :: subname='loewe_p'
  real(kind=8), allocatable :: ovrlp(:,:,:),evall(:),psitt(:,:)
  real(kind=8) :: dnrm2
  real(kind=8) :: tt,ttLOC
  integer :: nvctr_eff,i_all,i_stat,ierr,info,jorb,lorb

  if (norb == 1) then

     nvctr_eff=min(nvctr_tot-iproc*nvctrp,nvctrp)

     if (nvctr_eff > 0) then
     !parallel treatment of a run with only one orbital
     tt=dnrm2(nvctr_eff,psit,1)     
     ttLOC=tt**2

     else
        ttLOC =0.d0
     end if
     
     call MPI_ALLREDUCE(ttLOC,tt,1,MPI_DOUBLE_PRECISION,MPI_SUM,bigdft_mpi%mpi_comm,ierr)

     tt=1.d0/sqrt(tt)
     call vscal(nvctr_eff,tt,psit(1,1),1)

     !stop 'more than one orbital needed for a parallel run'

  else

     allocate(ovrlp(norb,norb,3+ndebug),stat=i_stat)
     call memocc(i_stat,ovrlp,'ovrlp',subname)
     allocate(evall(norb+ndebug),stat=i_stat)
     call memocc(i_stat,evall,'evall',subname)

     ! Upper triangle of overlap matrix using BLAS
     !     ovrlp(iorb,jorb)=psit(k,iorb)*psit(k,jorb) ; upper triangle
     call DSYRK('U','T',norb,nvctrp,1.d0,psit,nvctrp,0.d0,ovrlp(1,1,2),norb)

     ! Full overlap matrix using  BLAS
     !     ovrlap(jorb,iorb,2)=+psit(k,jorb)*psit(k,iorb)
     !      call DGEMM('T','N',norb,norb,nvctrp,1.d0,psit,&
     !                     nvctrp,psit,nvctrp,0.d0,ovrlp(1,1,2),norb)

     call MPI_ALLREDUCE(ovrlp(1,1,2),ovrlp(1,1,1),norb**2,&
          MPI_DOUBLE_PRECISION,MPI_SUM,bigdft_mpi%mpi_comm,ierr)

     !       write(*,*) 'OVERLAP',iproc
     !       do i=1,norb
     !       write(*,'(10(x,e17.10))') (ovrlp(i,j,1),j=1,norb)
     !       enddo

     ! LAPACK
     call DSYEV('V','U',norb,ovrlp(1,1,1),norb,evall,ovrlp(1,1,3),norb**2,info)
     if (info.ne.0) write(6,*) 'info loewe', info
     !        if (iproc.eq.0) then 
     !          write(6,*) 'overlap eigenvalues'
     !77        format(8(1x,e10.3))
     !          if (norb.le.16) then
     !          write(6,77) evall
     !          else
     !          write(6,77) (evall(i),i=1,4), (evall(i),i=norb-3,norb)
     !          endif
     !        endif

     ! calculate S^{-1/2} ovrlp(*,*,3)
     do lorb=1,norb
        do jorb=1,norb
           ovrlp(jorb,lorb,2)=ovrlp(jorb,lorb,1)*sqrt(1.d0/evall(lorb))
        end do
     end do
     !        do 3985,j=1,norb
     !        do 3985,i=1,norb
     !        ovrlp(i,j,3)=0.d0
     !        do 3985,l=1,norb
     !3985    ovrlp(i,j,3)=ovrlp(i,j,3)+ovrlp(i,l,1)*ovrlp(j,l,2)
     ! BLAS:
     call DGEMM('N','T',norb,norb,norb,1.d0,ovrlp(1,1,1),norb,&
          ovrlp(1,1,2),norb,0.d0,ovrlp(1,1,3),norb)

     allocate(psitt(nvctrp,ndim+ndebug),stat=i_stat)
     call memocc(i_stat,psitt,'psitt',subname)
     ! new eigenvectors
     !   psitt(i,iorb)=psit(i,jorb)*ovrlp(jorb,iorb,3)
     call DGEMM('N','N',nvctrp,norb,norb,1.d0,psit,nvctrp,ovrlp(1,1,3),norb,0.d0,psitt,nvctrp)
     call vcopy(nvctrp*ndim,psitt(1,1),1,psit(1,1),1)
     i_all=-product(shape(psitt))*kind(psitt)
     deallocate(psitt,stat=i_stat)
     call memocc(i_stat,i_all,'psitt',subname)

     i_all=-product(shape(ovrlp))*kind(ovrlp)
     deallocate(ovrlp,stat=i_stat)
     call memocc(i_stat,i_all,'ovrlp',subname)
     i_all=-product(shape(evall))*kind(evall)
     deallocate(evall,stat=i_stat)
     call memocc(i_stat,i_all,'evall',subname)

  end if

END SUBROUTINE loewe_p


!> Loewdin orthogonalisation
subroutine loewe(norb,nvctrp,psi)
  use module_base
  implicit real(kind=8) (a-h,o-z)
  dimension psi(nvctrp,norb)
  character(len=*), parameter :: subname='loewe'
  real(kind=8), allocatable :: ovrlp(:,:,:),evall(:),tpsi(:,:)

  if (norb.eq.1) then
     tt=0.d0
     do i=1,nvctrp
        tt=tt+psi(i,1)**2
     enddo
     tt=1.d0/sqrt(tt)
     do i=1,nvctrp
        psi(i,1)=psi(i,1)*tt
     enddo

  else

     allocate(ovrlp(norb,norb,3+ndebug),stat=i_stat)
     call memocc(i_stat,ovrlp,'ovrlp',subname)
     allocate(evall(norb+ndebug),stat=i_stat)
     call memocc(i_stat,evall,'evall',subname)

     ! Overlap matrix using BLAS
     !     ovrlp(iorb,jorb)=psi(k,iorb)*psi(k,jorb) ; upper triangle
     call DSYRK('U','T',norb,nvctrp,1.d0,psi,nvctrp,0.d0,ovrlp(1,1,1),norb)

     !       write(*,*) 'OVERLAP'
     !       do i=1,norb
     !       write(*,'(10(1x,1pe17.10))') (ovrlp(i,j,1),j=1,norb)
     !       enddo


     ! LAPACK
     call DSYEV('V','U',norb,ovrlp(1,1,1),norb,evall,ovrlp(1,1,3),norb**2,info)
     if (info.ne.0) write(6,*) 'info loewe', info
     !          write(6,*) 'overlap eigenvalues'
     !77        format(8(1x,e10.3))
     !          if (norb.le.16) then
     !          write(6,77) evall
     !          else
     !          write(6,77) (evall(i),i=1,4), (evall(i),i=norb-3,norb)
     !          endif

     ! calculate S^{-1/2} ovrlp(*,*,3)
     do lorb=1,norb
        do jorb=1,norb
           ovrlp(jorb,lorb,2)=ovrlp(jorb,lorb,1)*sqrt(1.d0/evall(lorb))
        end do
     end do
     !        do 3985,j=1,norb
     !        do 3985,i=1,norb
     !        ovrlp(i,j,3)=0.d0
     !        do 3985,l=1,norb
     !3985    ovrlp(i,j,3)=ovrlp(i,j,3)+ovrlp(i,l,1)*ovrlp(j,l,2)
     ! BLAS:
     call DGEMM('N','T',norb,norb,norb,1.d0,ovrlp(1,1,1),norb,ovrlp(1,1,2),norb,0.d0,ovrlp(1,1,3),norb)

     ! new eigenvectors
     allocate(tpsi(nvctrp,norb+ndebug),stat=i_stat)
     call memocc(i_stat,tpsi,'tpsi',subname)
     !   tpsi(i,iorb)=psi(i,jorb)*ovrlp(jorb,iorb,3)
     call DGEMM('N','N',nvctrp,norb,norb,1.d0,psi(1,1),nvctrp,ovrlp(1,1,3),norb,0.d0,tpsi,nvctrp)
     call vcopy(nvctrp*norb,tpsi(1,1),1,psi(1,1),1)
     i_all=-product(shape(tpsi))*kind(tpsi)
     deallocate(tpsi,stat=i_stat)
     call memocc(i_stat,i_all,'tpsi',subname)

     i_all=-product(shape(ovrlp))*kind(ovrlp)
     deallocate(ovrlp,stat=i_stat)
     call memocc(i_stat,i_all,'ovrlp',subname)
     i_all=-product(shape(evall))*kind(evall)
     deallocate(evall,stat=i_stat)
     call memocc(i_stat,i_all,'evall',subname)

  endif

END SUBROUTINE loewe

subroutine checkortho_paw(iproc,norb,nvctrp,psit,spsi)
  use module_base
  implicit real(kind=8) (a-h,o-z)
  integer :: iproc,norb,nvctrp
  dimension psit(nvctrp,norb)
  dimension spsi(nvctrp,norb)
  character(len=*), parameter :: subname='checkortho_paw'
  real(kind=8), allocatable :: ovrlp(:,:,:)

  allocate(ovrlp(norb,norb,2+ndebug),stat=i_stat)
  call memocc(i_stat,ovrlp,'ovrlp',subname)
  ovrlp=0.d0

  do iorb=1,norb
     do jorb=1,norb
        !<psi|S|psi>
        ovrlp(iorb,jorb,2)=ddot(nvctrp,psit(1,iorb),1,spsi(1,jorb),1)
        !add <psi|psi>
        ovrlp(iorb,jorb,2)=ovrlp(iorb,jorb,2)+&
         & ddot(nvctrp,psit(1,iorb),1,psit(1,jorb),1)
     end do
  end do

  call MPI_ALLREDUCE(ovrlp(1,1,2),ovrlp(1,1,1),norb**2,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)

  toler=1.d-10
  dev=0.d0
  do iorb=1,norb
     do jorb=1,norb
        scpr=ovrlp(iorb,jorb,1)
        if (iorb.eq.jorb) then
           dev=dev+(scpr-1.d0)**2
        else
           dev=dev+scpr**2
        endif
        if (iproc == 0) then
           if (iorb.eq.jorb .and. abs(scpr-1.d0).gt.toler) write(*,'(1x,a,2(1x,i0),1x,1pe12.6)')&
                'ERROR ORTHO',iorb,jorb,scpr
           if (iorb.ne.jorb .and. abs(scpr).gt.toler)      write(*,'(1x,a,2(1x,i0),1x,1pe12.6)')&
                'ERROR ORTHO',iorb,jorb,scpr
        end if
     end do
  end do

  if (dev.gt.toler) write(*,'(1x,a,i0,1pe13.5)') 'Deviation from orthogonality ',iproc,dev

  i_all=-product(shape(ovrlp))*kind(ovrlp)
  deallocate(ovrlp,stat=i_stat)
  call memocc(i_stat,i_all,'ovrlp',subname)

END SUBROUTINE checkortho_paw

subroutine checkortho_p(iproc,norb,nvctrp,psit)
  use module_base
  use module_types
  implicit real(kind=8) (a-h,o-z)
  dimension psit(nvctrp,norb)
  character(len=*), parameter :: subname='checkortho_p'
  real(kind=8), allocatable :: ovrlp(:,:,:)

  allocate(ovrlp(norb,norb,2+ndebug),stat=i_stat)
  call memocc(i_stat,ovrlp,'ovrlp',subname)

  do iorb=1,norb
     do jorb=1,norb
        ovrlp(iorb,jorb,2)=ddot(nvctrp,psit(1,iorb),1,psit(1,jorb),1)
     end do
  end do

  call MPI_ALLREDUCE(ovrlp(1,1,2),ovrlp(1,1,1),norb**2,MPI_DOUBLE_PRECISION,MPI_SUM,bigdft_mpi%mpi_comm,ierr)

  toler=1.d-10
  dev=0.d0
  do iorb=1,norb
     do jorb=1,norb
        scpr=ovrlp(iorb,jorb,1)
        if (iorb.eq.jorb) then
           dev=dev+(scpr-1.d0)**2
        else
           dev=dev+scpr**2
        endif
        if (iproc == 0) then
           if (iorb.eq.jorb .and. abs(scpr-1.d0).gt.toler) write(*,'(1x,a,2(1x,i0),1x,1pe12.6)')&
                'ERROR ORTHO',iorb,jorb,scpr
           if (iorb.ne.jorb .and. abs(scpr).gt.toler)      write(*,'(1x,a,2(1x,i0),1x,1pe12.6)')&
                'ERROR ORTHO',iorb,jorb,scpr
        end if
     end do
  end do

  if (dev.gt.toler) write(*,'(1x,a,i0,1pe13.5)') 'Deviation from orthogonality ',iproc,dev

  i_all=-product(shape(ovrlp))*kind(ovrlp)
  deallocate(ovrlp,stat=i_stat)
  call memocc(i_stat,i_all,'ovrlp',subname)

END SUBROUTINE checkortho_p


subroutine checkortho(norb,nvctrp,psi)
  use module_base
  implicit real(kind=8) (a-h,o-z)
  dimension psi(nvctrp,norb)
  character(len=*), parameter :: subname='checkortho'
  real(kind=8), allocatable :: ovrlp(:,:,:)

  allocate(ovrlp(norb,norb,1+ndebug),stat=i_stat)
  call memocc(i_stat,ovrlp,'ovrlp',subname)

  do iorb=1,norb
     do jorb=1,norb
        ovrlp(iorb,jorb,1)=ddot(nvctrp,psi(1,iorb),1,psi(1,jorb),1)
     enddo
  enddo

  toler=1.d-10
  dev=0.d0
  do iorb=1,norb
     do jorb=1,norb
        scpr=ovrlp(iorb,jorb,1)
        if (iorb.eq.jorb) then
           dev=dev+(scpr-1.d0)**2
        else
           dev=dev+scpr**2
        endif
        if (iorb.eq.jorb .and. abs(scpr-1.d0).gt.toler) write(*,'(1x,a,2(1x,i0),1x,1pe12.6)')&
             'ERROR ORTHO',iorb,jorb,scpr
        if (iorb.ne.jorb .and. abs(scpr).gt.toler)      write(*,'(1x,a,2(1x,i0),1x,1pe12.6)')&
             'ERROR ORTHO',iorb,jorb,scpr
     enddo
  enddo

  if (dev.gt.1.d-10) write(*,'(1x,a,i0,1pe13.5)') 'Deviation from orthogonality ',0,dev

  i_all=-product(shape(ovrlp))*kind(ovrlp)
  deallocate(ovrlp,stat=i_stat)
  call memocc(i_stat,i_all,'ovrlp',subname)


END SUBROUTINE checkortho


!> At the start each processor has all the Psi's but only its part of the HPsi's
!! at the end each processor has only its part of the Psi's
subroutine KStrans_p(nproc,norb,nvctrp,occup,  & 
     hpsit,psit,evsum,eval,nspinor)
  use module_base
  use module_types
  implicit none
  integer, intent(in) :: nproc,norb,nvctrp,nspinor
  real(wp), intent(out) :: evsum
  real(gp), dimension(norb), intent(in) :: occup
  real(wp), dimension(nvctrp*nspinor,norb), intent(in) :: hpsit
  real(wp), dimension(norb), intent(out) :: eval
  real(wp), dimension(nvctrp*nspinor,norb), intent(out) :: psit
  !local variables
  character(len=*), parameter :: subname='KStrans_p'
  integer :: i_all,i_stat,ierr,iorb,jorb,n_lp,istart,info,norbs,ncomp
  real(wp) :: alpha
  ! arrays for KS orbitals
  real(wp), dimension(:), allocatable :: work_lp,work_rp
  real(wp), dimension(:,:), allocatable :: psitt
  real(wp), dimension(:,:,:), allocatable :: hamks

  if(nspinor==4) then
     norbs=2*norb
     ncomp=2
  else if (nspinor==2) then
     norbs=2*norb
     ncomp=1
  else
     norbs=norb
  end if
  ! set up Hamiltonian matrix
  allocate(hamks(norbs,norb,2+ndebug),stat=i_stat)
  call memocc(i_stat,hamks,'hamks',subname)

  do jorb=1,norb
     do iorb=1,norbs
        hamks(iorb,jorb,2)=0.0_wp
     enddo
  enddo
  if (nproc > 1) then
     istart=2
  else
     istart=1
  end if

  if(nspinor==1) then
!     do iorb=1,norb
!        do jorb=1,norb
!           scpr=ddot(nvctrp,psit(1,jorb),1,hpsit(1,iorb),1)
!           hamks(iorb,jorb,istart)=scpr
!        enddo
!     enddo
     call gemm('T','N',norb,norb,nvctrp,1.0_wp,psit(1,1),max(1,nvctrp),hpsit(1,1),max(1,nvctrp),0.0_wp,&
          hamks(1,1,istart),norb)
  else
     call c_gemm('C','N',norb,norb,ncomp*nvctrp,(1.0_wp,0.0_wp),psit(1,1),max(1,ncomp*nvctrp), &
          hpsit(1,1),max(1,ncomp*nvctrp),(0.0_wp,0.0_wp),hamks(1,1,istart),norb)
  end if

  if (nproc > 1) then
     call MPI_ALLREDUCE(hamks(1,1,2),hamks(1,1,1),norbs*norb,mpidtypw,&
          MPI_SUM,bigdft_mpi%mpi_comm,ierr)
  end if
!  do iorb=1,norb
!     if(iproc==0) write(*,'(30f10.5)')(hamks(jorb,iorb,2),jorb=1,norbs)
!  end do
  !        write(*,*) 'KS Hamiltonian',iproc
  !        do iorb=1,norb
  !        write(*,'(10(1x,e10.3))') (hamks(iorb,jorb,1),jorb=1,norb)
  !        enddo

  n_lp=max(4*norbs,1000)
  allocate(work_lp(n_lp*2+ndebug),stat=i_stat)
  call memocc(i_stat,work_lp,'work_lp',subname)
  if(nspinor==1) then
     call  syev('V','U',norb,hamks(1,1,1),norb,eval(1),work_lp(1),n_lp,info)
  else
     allocate(work_rp(3*norb+1+ndebug),stat=i_stat)
     call memocc(i_stat,work_rp,'work_rp',subname)
     
     call  heev('V','U',norb,hamks(1,1,1),norb,eval(1),work_lp(1),n_lp,work_rp(1),info)
     
     i_all=-product(shape(work_rp))*kind(work_rp)
     deallocate(work_rp,stat=i_stat)
     call memocc(i_stat,i_all,'work_rp',subname)
  end if

  evsum=0.0_wp
  do iorb=1,norb
     evsum=evsum+eval(iorb)*real(occup(iorb),wp)
     !if (iproc.eq.0) write(*,'(1x,a,i0,a,1x,1pe21.14)') 'eval(',iorb,')=',eval(iorb)
  enddo
  i_all=-product(shape(work_lp))*kind(work_lp)
  deallocate(work_lp,stat=i_stat)
  call memocc(i_stat,i_all,'work_lp',subname)
  if (info.ne.0) write(*,*) 'DSYEV ERROR',info

  allocate(psitt(nvctrp*nspinor,norb+ndebug),stat=i_stat)
  call memocc(i_stat,psitt,'psitt',subname)
  ! Transform to KS orbitals
  ! dgemm can be used instead of daxpy
  if(nspinor==1) then
     do iorb=1,norb
        call to_zero(nvctrp,psitt(1,iorb))
        do jorb=1,norb
           alpha=hamks(jorb,iorb,1)
           call axpy(nvctrp,alpha,psit(1,jorb),1,psitt(1,iorb),1)
        enddo
     enddo
  else
     do iorb=1,norb
        call to_zero(nvctrp*nspinor,psitt(1,iorb))
        do jorb=1,norb
           call c_axpy(ncomp*nvctrp,hamks(2*jorb-1,iorb,1),psit(1,jorb),1,psitt(1,iorb),1)
        enddo
     enddo
  end if
  i_all=-product(shape(hamks))*kind(hamks)
  deallocate(hamks,stat=i_stat)
  call memocc(i_stat,i_all,'hamks',subname)

  call vcopy(nvctrp*norb*nspinor,psitt(1,1),1,psit(1,1),1)
  i_all=-product(shape(psitt))*kind(psitt)
  deallocate(psitt,stat=i_stat)
  call memocc(i_stat,i_all,'psitt',subname)

END SUBROUTINE KStrans_p


!> This subroutine orthonormalizes the orbitals psi in a parallel way. To do so, it first transposes the orbitals to all
!! processors using mpi_alltoallv. The orthonomalization is then done in this data layout using a combination of blockwise Gram-Schmidt
!! and Cholesky orthonomalization. At the end the vectors are again untransposed.
!!
!! Input arguments:
!!  @param  iproc      process ID
!!  @param  nproc      total number of processes
!!  @param  norb       total number of vectors that have to be orthonomalized, shared over all processes
!!  @param  orthpar    data type containing many parameters
!!  @param  nspinor    size of spinor
!!  @param  nspin      spin components
!!  @param  ndilmovrlp dimension of overlap
!!  @param  norbArr
!!  @param  comms      Communication arrays
!! Input/Output arguments:
!!  @param  psi
!!      - on input: the vectors to be orthonormalized
!!      - on output: the orthonomalized vectors
subroutine gsChol(iproc, nproc, psi, orthpar, nspinor, orbs, nspin,ndim_ovrlp,norbArr,comms,paw)
  use module_base
  use module_types
  use module_interfaces, except_this_one_A => gsChol
  use communications_base, only: comms_cubic
  implicit none

  ! Calling arguments
  !integer, intent(in) :: ikpt
  integer, intent(in) :: iproc, nproc, nspin
  integer, intent(inout) :: nspinor
  type(orthon_data), intent(in):: orthpar
  type(orbitals_data):: orbs
  type(comms_cubic), intent(in) :: comms
  integer, dimension(nspin), intent(in) :: norbArr
  integer, dimension(nspin,0:orbs%nkpts), intent(inout) :: ndim_ovrlp
  real(wp),dimension(comms%nvctr_par(iproc,0)*orbs%nspinor*orbs%norb),intent(inout):: psi
  type(paw_objects),optional,intent(inout)::paw
  
  ! Local variables
  integer:: iblock, jblock, ist, jst, iter, iter2, gcd, blocksize, blocksizeSmall, i_stat, i_all,usepaw=0
  integer:: getBlocksize, ispin
  real(wp),dimension(:), allocatable :: ovrlp
  character(len=*), parameter:: subname='gsChol',category='GS/Chol'
  
  if(present(paw))usepaw=paw%usepaw  

  ! Make a loop over spin up/down.
  do ispin=1,nspin
     ! Get the blocksize.
     blocksize=getBlocksize(orthpar, norbArr(ispin))
     
     ! There are two orthonormalization subroutines: gramschmidt orthogonalizes a given bunch of vectors to another bunch
     ! of already orthonormal vectors, and the subroutine cholesky orthonormalizes the given bunch.
     ! First determine how many bunches can be created for the given blocksize.
     iter=floor(real(norbArr(ispin))/real(blocksize))
     
     ! Get the dimensions of the overlap matrix for handling blocksize orbitals.
     call dimension_ovrlpFixedNorb(nspin,orbs,ndim_ovrlp,blocksize)
     allocate(ovrlp(ndim_ovrlp(nspin,orbs%nkpts)+ndebug),stat=i_stat)
     call memocc(i_stat,ovrlp,'ovrlp',subname)
     
     ! Make a loop over all blocks.
     do iblock=1,iter
        ! ist is the starting orbitals of the current bunch of vectors.
        ist=(iblock-1)*blocksize+1
        ! Now orthogonalize this bunch to all previous ones.
        do jblock=1,iblock-1
           ! jst is the starting vector of the bunch to which the current bunch has to be orthogonalized.
           jst=blocksize*(jblock-1)+1
           if(usepaw==1) then
              call getOverlapDifferentPsi_paw(iproc, nproc, nspin, blocksize,orbs, &
                   comms, psi(1),paw%spsi(1), ndim_ovrlp, ovrlp, norbArr, ist, jst, ispin, category)
              call gramschmidt(iproc, blocksize, psi(1), ndim_ovrlp, ovrlp, &
                   orbs, nspin, nspinor, comms, norbArr, ist, jst, ispin,paw)
           else
              call getOverlapDifferentPsi(iproc, nproc, nspin, blocksize,orbs, &
                   comms, psi(1), ndim_ovrlp, ovrlp, norbArr, ist, jst, ispin, category)
              call gramschmidt(iproc, blocksize, psi(1), ndim_ovrlp, ovrlp, &
                   orbs, nspin, nspinor, comms, norbArr, ist, jst, ispin)
           end if
        end do
    
        ! Orthonormalize the current bunch of vectors.
        if(usepaw==1) then
           call getOverlap_paw(iproc, nproc, nspin, blocksize, orbs, comms, psi(1), &
                paw%spsi(1),ndim_ovrlp, ovrlp, norbArr, ist, ispin, category)
           call cholesky(iproc, nspin,blocksize, psi(1), orbs, &
                comms, ndim_ovrlp, ovrlp(1), norbArr, ist, ispin,paw)
        else
           call getOverlap(iproc, nproc, nspin, blocksize, orbs, comms, psi(1), &
                ndim_ovrlp, ovrlp, norbArr, ist, ispin, category)
           call cholesky(iproc, nspin, blocksize, psi(1), orbs, &
                comms, ndim_ovrlp, ovrlp(1), norbArr, ist, ispin)
        end if
    
    end do

    i_all=-product(shape(ovrlp))*kind(ovrlp)
    deallocate(ovrlp, stat=i_stat)
    call memocc(i_stat,i_all,'ovrlp',subname)
    

    ! Orthonormalize the remaining vectors, if there are any.
    remainingIf: if(blocksize*iter/=norbArr(ispin)) then
        ! ist is the starting vector of the bunch that still havs to be orthonormalized.
        ist=blocksize*iter+1

        ! We have to find a new block size that matches both the remaining vectors and the already orthonomalized ones. This is done by determining
        ! the greatest common divisor of these two numbers.
        blocksizeSmall=gcd(blocksize*iter,norbArr(ispin)-ist+1)

        ! Get the dimensions of the overlap matrix for handling blocksize orbitals.
        call dimension_ovrlpFixedNorb(nspin,orbs,ndim_ovrlp,blocksizeSmall)
        allocate(ovrlp(ndim_ovrlp(nspin,orbs%nkpts)+ndebug),stat=i_stat)
        call memocc(i_stat,ovrlp,'ovrlp',subname)

        ! Determine how many blocks can be created with this new block size.
        iter2=(norbArr(ispin)-ist+1)/blocksizeSmall
        ! Now make a loop over all these blocks.
        do iblock=1,iter2
            ! ist is the starting vector of the current bunch.
            ist=iter*blocksize+blocksizeSmall*(iblock-1)+1
            ! Now orthogonalize this bunch to all previous ones.
            do jblock=1,(blocksize*iter)/blocksizeSmall+iblock-1
                ! jst is the starting vector of the bunch to which the current bunch has to be orthogonalized.
                jst=blocksizeSmall*(jblock-1)+1
                if(usepaw==1) then
                   call getOverlapDifferentPsi_paw(iproc, nproc, nspin, blocksizeSmall, &
                        orbs, comms, psi(1), paw%spsi(1),ndim_ovrlp, ovrlp, norbArr, ist, jst, ispin, category)
                   call gramschmidt(iproc, blocksizeSmall, psi(1), ndim_ovrlp, &
                        ovrlp, orbs, nspin, nspinor, comms, norbArr, ist, jst, ispin,paw)
                else
                   call getOverlapDifferentPsi(iproc, nproc, nspin, blocksizeSmall, &
                        orbs, comms, psi(1), ndim_ovrlp, ovrlp, norbArr, ist, jst, ispin, category)
                   call gramschmidt(iproc, blocksizeSmall, psi(1), ndim_ovrlp, &
                        ovrlp, orbs, nspin, nspinor, comms, norbArr, ist, jst, ispin)
                end if

            end do
            ! Orthonormalize the current bunch of vectors.
            if(usepaw==1) then
               call getOverlap_paw(iproc, nproc, nspin, blocksizeSmall, orbs, comms,&
                    psi(1), paw%spsi(1),ndim_ovrlp, ovrlp, norbArr, ist, ispin, category)
               call cholesky(iproc, nspin, blocksizeSmall, psi(1), &
                    orbs, comms, ndim_ovrlp, ovrlp(1), norbArr, ist, ispin,paw)
            else
               call getOverlap(iproc, nproc, nspin, blocksizeSmall, orbs, comms,&
                    psi(1), ndim_ovrlp, ovrlp, norbArr, ist, ispin, category)
               call cholesky(iproc, nspin, blocksizeSmall, psi(1), &
                    orbs, comms, ndim_ovrlp, ovrlp(1), norbArr, ist, ispin)
            end if
        end do
        i_all=-product(shape(ovrlp))*kind(ovrlp)
        deallocate(ovrlp, stat=i_stat)
        call memocc(i_stat,i_all,'ovrlp',subname)
    end if remainingIf
    
end do

END SUBROUTINE gsChol


!>  This subroutine orthogonalizes a given bunch of vectors in psit to another bunch of equal size. These other vectors
!!  are assumed to be orthonomal themselves. The starting indices of the two bunches are given by block1 and block2.
!!  The orthonormalization is done in parallel, assuming that each process holds a small portion of each vector.
!!
!!  Input arguments:
!!   @param  iproc       process ID
!!   @param  norbIn     number of orbitals to be orthonormalized
!!   @param  ndim_ovrlp      describes the shape of the overlap matrix
!!   @param  orbs       type that contains many parameters concerning the orbitals
!!   @param  nspin      closed shell -> nspin=1 ; spin polarised -> nspin=2
!!   @param  nspinor    real wavefunction -> nspinor=1, complex wavefunction -> nspinor>1
!!   @param  comms      type containing parameters for communicating the wavefunstion between processors
!!   @param  norbTot    total number of orbitals (if nspin=2:
!!                 norbTot(1)=total number of up orbitals
!!                 norbTot(2)=total number of down orbitals)
!!   @param  block1     gives the starting orbital of the orbitals to be orthogonalized
!!   @param  block2     gives the starting orbital of the orbitals to which they shall be orthogonalized
!!   @param  ispinIn    indicates whether the up or down orbitals shall be handled
!!  Input/Output arguments:
!!   @param  psit       the vectors that shall be orthonormalized
!!   @param  ovrlp      the overlap matrix which will be destroyed during this subroutine
subroutine gramschmidt(iproc, norbIn, psit, ndim_ovrlp, ovrlp, orbs, nspin,&
     nspinor, comms, norbTot, block1, block2, ispinIn, paw)
use module_base
use module_types
use communications_base, only: comms_cubic
implicit none

! Calling arguments
integer,intent(in):: iproc, norbIn, nspin, block1, block2, ispinIn
integer,intent(out) :: nspinor
type(orbitals_data):: orbs
type(comms_cubic), intent(in) :: comms
type(paw_objects),optional,intent(inout)::paw
real(wp),dimension(comms%nvctr_par(iproc,0)*orbs%nspinor*orbs%norb),intent(inout):: psit
integer,dimension(nspin,0:orbs%nkpts):: ndim_ovrlp
real(wp),dimension(ndim_ovrlp(nspin,orbs%nkpts)):: ovrlp
integer,dimension(nspin):: norbTot

! Local arguments
integer:: nvctrp, i_stat, i_all, ncomp, ikptp, ikpt, ispin, norb, norbs, istThis, istOther,usepaw=0
!real(kind=8),allocatable::raux(:)
real(kind=8),dimension(:),allocatable:: A1D
character(len=*),parameter:: subname='gramschmidt'

if(present(paw))usepaw=paw%usepaw

! Initialize the starting indices. istThis is the starting index of the orbitals that shall be orthogonalized,
! istOther is the starting index of the orbitals to which they shall be orthogonalized.
istThis=1
istOther=1

! Make a loop over the number of k-points handled by the process.
do ikptp=1,orbs%nkptsp
    ! ikpt is the number of the k-point.
    ikpt=orbs%iskpts+ikptp
    ! Now make a loop over spin up and down.
    do ispin=1,nspin
        ! This subroutine gives essentially back nvctrp, i.e. the length of the vectors for.
        ! In addition it sets the value of nspinor to orbs%nspinor.
        call orbitals_and_components(iproc,ikpt,ispin,orbs,comms,&
            nvctrp,norb,norbs,ncomp,nspinor)
        ! The subroutine also overwrite the variable norb with the total number of orbitals.
        ! However we want to keep the value of norbIn (since we possibly treat only a part of the orbitals).
        norb=norbIn

        ! Allocate the matrix A which will hold some partial results.
        allocate(A1D(nvctrp*norb*nspinor), stat=i_stat)
        call memocc(i_stat, A1D, 'A1D', subname)

        ! Count up the starting indices.
        istThis=istThis+nvctrp*(block1-1)*nspinor
        istOther=istOther+nvctrp*(block2-1)*nspinor

        if(ispin==ispinIn) then
            ! Calculate matrix product psit*ovrlp=A. This will give the components that will be projected out of psit.
            ! We actually calculate -psit*ovrlp=-A, since this is better for further processing with daxpy.
            if(nspinor==1) then
                call dgemm('n', 'n', nvctrp, norb, norb, -1.d0, psit(istOther), nvctrp,&
                     ovrlp(ndim_ovrlp(ispin,ikpt-1)+1), norb, 0.d0, A1D(1), nvctrp)
            else
                call zgemm('n', 'n', nvctrp, norb, norb, (-1.d0,0.d0), psit(istOther), nvctrp, &
                     ovrlp(ndim_ovrlp(ispin,ikpt-1)+1), norb, (0.d0,0.d0), A1D(1), nvctrp)
            end if
            ! Now project out: psit=psit-A.
            ! Since we calculated -A, we have to put psit=psit+A and can use daxpy to perform psit=A+psit
            if(nspinor==1) then
                call daxpy(nvctrp*norb*nspinor,1.d0,A1D(1),1,psit(istThis),1)
            else
                call daxpy(nvctrp*norb*nspinor,1.d0,A1D(1),1,psit(istThis),1)
            end if

           if(paw%usepaw==1) then
           !Do the same for SPSI and cprj:
           !Pending: This is not yet coded
              stop

              ! We actually calculate -psit*ovrlp=-A, since this is better for further processing with daxpy.
              if(nspinor==1) then
                  call dgemm('n', 'n', nvctrp, norb, norb, -1.d0, paw%spsi(istOther), nvctrp,&
                       ovrlp(ndim_ovrlp(ispin,ikpt-1)+1), norb, 0.d0, A1D(1), nvctrp)
              else
                  call zgemm('n', 'n', nvctrp, norb, norb, (-1.d0,0.d0), paw%spsi(istOther), nvctrp, &
                       ovrlp(ndim_ovrlp(ispin,ikpt-1)+1), norb, (0.d0,0.d0), A1D(1), nvctrp)
              end if
              call daxpy(nvctrp*norb*nspinor,1.d0,A1D(1),1,paw%spsi(istThis),1)

              !update cprj
              !allocate(raux(2*paw%lmnmax*norb*nspinor))
              !call memocc(i_stat,raux,'raux',subname)
              !!
              !do iat=1,paw%natom
              !   raux=0.d0
              !   !copy cprj%cp object to a simple array 'raux'
              !   call cprj_to_array(paw%cprj(iat,:),raux,norb,nspinor,istThis-1,1)
              !   !
              !   !PENDING: 
              !   stop
              !   !if(nspinor==1) then
              !   !    call dgemm('n', 'n', 2*paw%lmnmax, norb, norb, -1.d0, raux(istOther), nvctrp,&
              !   !         ovrlp(ndim_ovrlp(ispin,ikpt-1)+1), norb, 0.d0, A1D(1), nvctrp)
              !   !else
              !   !    call zgemm('n', 'n', nvctrp, norb, norb, (-1.d0,0.d0), paw%spsi(istOther), nvctrp, &
              !   !         ovrlp(ndim_ovrlp(ispin,ikpt-1)+1), norb, (0.d0,0.d0), A1D(1), nvctrp)
              !   !end if
              !   !call daxpy(2*paw%lmnmax*norb*nspinor,1.d0,A1D(1),1,raux,1)
              !end do
              !!
              !i_all=-product(shape(raux))*kind(raux)
              !deallocate(raux,stat=i_stat)
              !call memocc(i_stat,i_all,'raux',subname)
           end if

        end if


        ! Increase the starting indices. This will bring the starting index to the start of the the next spin case (up/down) and k-point.
        istThis=istThis+nvctrp*(norbTot(ispin)-block1+1)*nspinor
        istOther=istOther+nvctrp*(norbTot(ispin)-block2+1)*nspinor

        i_all=-product(shape(A1D))*kind(A1D)
        deallocate(A1D)
        call memocc(i_stat,i_all,'A1D',subname)
    end do
end do

END SUBROUTINE gramschmidt


!>  This subroutine orthonormalizes a given bunch of vectors psi.
!!  It first calculates the Cholesky factorization S=L*L^T of the overlap matrix S.  This matrix L is then
!!  inverted to get L^{-1} and the orthonormal vectors are finally given by psi=psi*L^{-1}.
!!
!!  Input arguments:
!!   @param  iproc      process ID
!!   @param  norbIn     number of orbitals to be orthonormalized
!!   @param  nspinor    real wavefunction -> nspinor=1, complex wavefunction -> nspinor>1
!!   @param  nspin      closed shell -> nspin=1 ; spin polarised -> nspin=2
!!   @param  orbs       type that contains many parameters concerning the orbitals
!!   @param  comms      type containing parameters for communicating the wavefunction between processors
!!   @param  ndim_ovrlp describes the shape of the overlap matrix
!!   @param  norbTot    total number of orbitals (if nspin=2:
!!                      - norbTot(1)=total number of up orbitals
!!                      - norbTot(2)=total number of down orbitals)
!!   @param  block1     gives the starting orbital of the orbitals to be orthonormalized
!!   @param  ispinIn    indicates whether the up or down orbitals shall be handled
!!  Input/Output arguments:
!!   @param  psi        the vectors that shall be orthonormalized
!!   @param  Lc      the overlap matrix which will be destroyed during this subroutine
subroutine cholesky(iproc,nspin, norbIn, psi, orbs, comms, ndim_ovrlp, ovrlp, norbTot, block1, ispinIn,paw)

use module_base
use module_types
use communications_base, only: comms_cubic
implicit none

! Calling arguments
!integer:: iproc,nvctrp,norbIn, nspinor, nspin, norbTot, block1, ispinIn
integer:: iproc,nvctrp,norbIn, block1, ispinIn,nspin
type(orbitals_data):: orbs
type(comms_cubic):: comms
real(kind=8),dimension(orbs%npsidim_comp),intent(inout):: psi
integer,dimension(nspin,0:orbs%nkpts):: ndim_ovrlp
real(kind=8),dimension(ndim_ovrlp(nspin,orbs%nkpts),1):: ovrlp
integer,dimension(orbs%nspin):: norbTot
type(paw_objects),optional,intent(inout)::paw

! Local variables
integer:: ist, info, ispin, ikptp, ikpt, ncomp, norbs, norb,nspinor
integer:: i_all,i_stat,iat
integer:: usepaw=0
real(kind=8),dimension(:,:),allocatable::raux
character(len=*),parameter:: subname='cholesky'

if(present(paw))usepaw=paw%usepaw
 
! Set the starting index to 1.
ist=1
! Make a loop over the number of k-points handled by the process.
do ikptp=1,orbs%nkptsp
    ! ikpt is the number of the k-point.
    ikpt=orbs%iskpts+ikptp
    ! Now make a loop over spin up and down.
    do ispin=1,nspin
        ! This subroutine gives essentially back nvctrp, i.e. the length of the vectors for.
        ! In addition it sets the value of nspinor to orbs%nspinor.
        call orbitals_and_components(iproc,ikpt,ispin,orbs,comms,&
            nvctrp,norb,norbs,ncomp,nspinor)
        ! The subroutine also overwrite the variable norb with the total number of orbitals.
        ! However we want to keep the value of norbIn (since we possibly treat only a part of the orbitals).

        norb=norbIn
        ! Count up the starting index
        ist=ist+nvctrp*(block1-1)*nspinor
        
 
        ! The following part is only executed if ispin==ispinIn. Otherwise only the starting index ist
        ! is increased.
        if(ispin==ispinIn) then
            
            ! Make a Cholesky factorization of L.
            if(nspinor==1) then
                call dpotrf('l', norb, ovrlp(ndim_ovrlp(ispin,ikpt-1)+1,1), norb, info)
            else
                call zpotrf('l', norb, ovrlp(ndim_ovrlp(ispin,ikpt-1)+1,1), norb, info)
            end if
            !print *,'info',info
            ! Invert the Cholesky matrix: L^{-1}.
            if(nspinor==1) then
                call dtrtri('l', 'n', norb, ovrlp(ndim_ovrlp(ispin,ikpt-1)+1,1), norb, info)
            else
                call ztrtri('l', 'n', norb, ovrlp(ndim_ovrlp(ispin,ikpt-1)+1,1), norb, info)
            end if
            !print *,'info',info
            ! Calculate the matrix product psi*L^{-1}=psi. This will give the orthonormal orbitals.
            ! For PAW: update spsi, and cprj (below)
            if(nspinor==1) then
                call dtrmm('r', 'l', 't', 'n', nvctrp, norb, 1.d0, &
                     ovrlp(ndim_ovrlp(ispin,ikpt-1)+1,1), norb, psi(ist), nvctrp)
                if(usepaw==1) then
                   call dtrmm('r', 'l', 't', 'n', nvctrp, norb, 1.d0, &
                        ovrlp(ndim_ovrlp(ispin,ikpt-1)+1,1), norb, paw%spsi(ist), nvctrp)
                end if
            else
                call ztrmm('r', 'l', 'c', 'n', ncomp*nvctrp, norb, (1.d0,0.d0),&
                     ovrlp(ndim_ovrlp(ispin,ikpt-1)+1,1), norb, psi(ist), ncomp*nvctrp)
                if(usepaw==1) then
                   call ztrmm('r', 'l', 'c', 'n', ncomp*nvctrp, norb, (1.d0,0.d0),&
                        ovrlp(ndim_ovrlp(ispin,ikpt-1)+1,1), norb, paw%spsi(ist), ncomp*nvctrp)
                end if
            end if

            if(usepaw==1) then
              !Pending: check that this works in parallel, and with nspinor=2
              !update cprj
              allocate(raux(2*paw%lmnmax,norb*nspinor),stat=i_stat)
              call memocc(i_stat,raux,'raux',subname)
              do iat=1,paw%natom
                raux=0.d0
                !copy cprj%cp objet to a simple array 'raux'
                call cprj_to_array(paw%cprj(iat,:),raux,norb,nspinor,ndim_ovrlp(ispin,ikpt-1),1)
                ! Calculate the matrix product cprj*L^{-1}=cprj.
                if(nspinor==1) then
                   call dtrmm('r', 'l', 't', 'n', 2*paw%lmnmax, norb, 1.d0, &
                        ovrlp(ndim_ovrlp(ispin,ikpt-1)+1,1), norb, raux, 2*paw%lmnmax)
                else
                   call ztrmm('r', 'l', 'c', 'n', ncomp*2*paw%lmnmax, norb, (1.d0,0.d0),&
                        ovrlp(ndim_ovrlp(ispin,ikpt-1)+1,1), norb, raux, ncomp*2*paw%lmnmax)
                end if
                !
                !copy back raux to cprj%cp
                call cprj_to_array(paw%cprj(iat,:),raux,norb,nspinor,ndim_ovrlp(ispin,ikpt-1),2)
              end do
              i_all=-product(shape(raux))*kind(raux)
              deallocate(raux,stat=i_stat)
              call memocc(i_stat,i_all,'raux',subname)
 
            end if !usepaw
        end if !InSpin


 
        ! Increase the starting index.
        ist=ist+nvctrp*(norbTot(ispin)-block1+1)*nspinor

    end do
end do         

END SUBROUTINE cholesky


!> Orthonormalizes the vectors provided in psit by a loewdin orthonormalization.
!!  Input arguments:
!!   @param  iproc      process ID
!!   @param  norbIn     number of orbitals to be orthonormalized
!!   @param  nspinor    real wavefunction -> nspinor=1, complex wavefunction -> nspinor>1
!!   @param  block1     gives the starting orbital of the orbitals to be orthonormalized
!!   @param  ispinIn    indicates whether the up or down orbitals shall be handled
!!   @param  orbs       type that contains many parameters concerning the orbitals
!!   @param  comms      type containing parameters for communicating the wavefunstion between processors
!!   @param  nspin      closed shell -> nspin=1 ; spin polarised -> nspin=2
!!   @param  ndim_ovrlp describes the shape of the overlap matrix
!!   @param  norbTot    total number of orbitals (if nspin=2:
!!                      - norbTot(1)=total number of up orbitals
!!                      - norbTot(2)=total number of down orbitals)
!!  Input/output Arguments
!!   @param  psit       the orbitals to be orthonormalized
!!   @param  ovrlp      the overlap matrix which will be destroyed during this subroutine
subroutine loewdin(iproc, norbIn, nspinor, block1, ispinIn, orbs, comms, nspin, psit, ovrlp, ndim_ovrlp, norbTot,paw)

use module_base
use module_types
use communications_base, only: comms_cubic
implicit none

! Calling arguments
integer,intent(in):: iproc,norbIn, nspin, block1, ispinIn
type(paw_objects),optional,intent(inout)::paw
integer, intent(inout) :: nspinor
type(orbitals_data),intent(in):: orbs
type(comms_cubic),intent(in):: comms
real(kind=8),dimension(comms%nvctr_par(iproc,0)*orbs%nspinor*orbs%norb),intent(in out):: psit
integer,dimension(nspin,0:orbs%nkpts):: ndim_ovrlp
real(kind=8),dimension(ndim_ovrlp(nspin,orbs%nkpts)):: ovrlp
integer,dimension(nspin):: norbTot

! Local variables
integer:: jorb, lorb, i_stat, i_all, info, nvctrp, ispin, ist, ikptp, ikpt, ncomp, norbs, norb, lwork,usepaw=0
integer:: ii,iat,jj,shift,ispinor,iorb,ilmn
real(kind=8),allocatable::raux(:,:,:,:)
real(kind=8),dimension(:),allocatable:: evall, psitt
real(kind=8),dimension(:,:),allocatable:: tempArr
character(len=*), parameter :: subname='loewdin'

if(present(paw))usepaw=paw%usepaw

! Allocate the work arrays.
lwork=nspinor*norbIn**2+10
allocate(tempArr(norbIn**2*nspinor,2), stat=i_stat)
call memocc(i_stat,tempArr,'tempArr',subname)

allocate(evall(norbIn), stat=i_stat)
call memocc(i_stat,evall,'evall',subname)

ist=1
! Make a loop over the number of k-points handled by the process.
do ikptp=1,orbs%nkptsp
    ! ikpt is the number of the k-point.
    ikpt=orbs%iskpts+ikptp
    ! Now make a loop over spin up and down.
    do ispin=1,nspin
        ! This subroutine gives essentially back nvctrp, i.e. the length of the vectors for.
        ! In addition it sets the value of nspinor to orbs%nspinor.
        call orbitals_and_components(iproc,ikpt,ispin,orbs,comms,&
            nvctrp,norb,norbs,ncomp,nspinor)
        ! The subroutine also overwrite the variable norb with the total number of orbitals.
        ! However we want to keep the value of norbIn (since we possibly treat only a part of the orbitals).
        norb=norbIn
        ! Count up the starting index
        ist=ist+nvctrp*(block1-1)*nspinor
 
        ! The following part is only executed if ispin==ispinIn. Otherwise only the starting index ist
        ! is increased.
        if(ispin==ispinIn) then

            ! Diagonalize the overlap matrix.
            if(nspinor==1) then
                call dsyev('v', 'l', norb, ovrlp(ndim_ovrlp(ispin,ikpt-1)+1), norb,&
                     evall, tempArr(1,1), lwork, info)
  
           
            else
                call zheev('v', 'l', norb,ovrlp(ndim_ovrlp(ispin,ikpt-1)+1), norb,&
                     evall, tempArr(1,1), lwork, tempArr(1,2), info)
            end if
            if (info/=0) then
                write(*,'(a,i0)') 'ERROR in dsyev (Loewdin); info=',info
                stop
            end if

            ! Calculate S^{-1/2}. 
            ! First calulate ovrlp*diag(evall) (ovrlp is the diagonalized overlap
            ! matrix and diag(evall) the diagonal matrix consisting of the eigenvalues...
            do lorb=1,norb
               do jorb=1,norb*nspinor
                  tempArr((lorb-1)*norb*nspinor+jorb,1)=&
                       ovrlp(ndim_ovrlp(ispin,ikpt-1)+(lorb-1)*norb*nspinor+jorb)*sqrt(1.d0/evall(lorb))
               end do
            end do

            ! ...and now apply the diagonalized overlap matrix to the matrix constructed above.
            ! This will give S^{-1/2}.
            if(nspinor==1) then
                call dgemm('n', 't', norb, norb, norb, 1.d0, ovrlp(ndim_ovrlp(ispin,ikpt-1)+1), norb,&
                     tempArr(1,1), norb, 0.d0, tempArr(1,2), norb)
            else
                call zgemm('n', 't', norb, norb, norb, (1.d0,0.d0), ovrlp(ndim_ovrlp(ispin,ikpt-1)+1), norb,&
                     tempArr(1,1), norb, (0.d0,0.d0), tempArr(1,2), norb)
            end if

            ! Now calculate the orthonormal orbitals by applying S^{-1/2} to the orbitals.
            ! This requires the use of a temporary variable psitt.
            allocate(psitt(nvctrp*norb*nspinor),stat=i_stat)
            call memocc(i_stat,psitt,'psitt',subname)
            if(nspinor==1) then
                call dgemm('n', 'n', nvctrp, norb, norb, 1.d0, psit(ist), &
                     nvctrp, tempArr(1,2), norb, 0.d0, psitt, nvctrp)
            else
                call zgemm('n', 'n', nvctrp, norb, norb, (1.d0,0.d0), &
                     psit(ist), nvctrp, tempArr(1,2), norb, (0.d0,0.d0), psitt, nvctrp)
            end if

            ! Now copy the orbitals from the temporary variable to psit.
            call vcopy(nvctrp*norb*nspinor, psitt(1), 1, psit(ist), 1)

            ! For PAW: upgrade also spsi and cprj
            if(usepaw==1) then

               if(nspinor==1) then
                   call dgemm('n', 'n', nvctrp, norb, norb, 1.d0, paw%spsi(ist), &
                        nvctrp, tempArr(1,2), norb, 0.d0, psitt, nvctrp)
               else
                   call zgemm('n', 'n', nvctrp, norb, norb, (1.d0,0.d0), &
                        paw%spsi(ist), nvctrp, tempArr(1,2), norb, (0.d0,0.d0), psitt, nvctrp)
               end if

               ! Now copy the orbitals from the temporary variable to psit.
               call vcopy(nvctrp*norb*nspinor, psitt(1), 1, paw%spsi(ist), 1)
               
               !Now upgrade cprj:
               !Pending: check that this works for more than 1 orbital, and in parallel
               !update cprj
               !icprj=icprj+(block1-1)*nspinor
               allocate(raux(2,paw%lmnmax,paw%natom,norb*nspinor))
               call memocc(i_stat,raux,'raux',subname)
               raux=0.d0
               ii=0
               do iorb=1,norb
                 jj=0
                 do jorb=1,norb
                   shift=(iorb-1)*norb*nspinor+jorb
                   do ispinor=1,nspinor
                     ii=ii+1
                     jj=jj+1
                     do iat=1,paw%natom
                       do ilmn=1,paw%cprj(iat,jj)%nlmn
                       raux(:,ilmn,iat,iorb)=raux(:,ilmn,iat,iorb)&
                         +tempArr(ii,2)*paw%cprj(iat,jj)%cp(:,ilmn)
                       end do
                     end do
                   end do
                 end do
               end do
               jj=0
               do iorb=1,norb
                 do ispinor=1,nspinor
                   jj=jj+1
                   do iat=1,paw%natom
                     do ilmn=1,paw%cprj(iat,jj)%nlmn
                       paw%cprj(iat,jj)%cp(:,ilmn)=raux(:,ilmn,iat,jj)
                     end do
                   end do
                 end do
               end do
               !
               i_all=-product(shape(raux))*kind(raux)
               deallocate(raux,stat=i_stat)
               call memocc(i_stat,i_all,'raux',subname)
            end if !usepaw

            ! Deallocate the temporary variable psitt.
            i_all=-product(shape(psitt))*kind(psitt)
            deallocate(psitt,stat=i_stat)
            call memocc(i_stat,i_all,'psitt',subname)

        end if
        ! Increase the starting index.
        ist=ist+nvctrp*(norbTot(ispin)-block1+1)*nspinor

    end do

end do         


! Deallocate the remaining arrays.
i_all=-product(shape(tempArr))*kind(tempArr)
deallocate(tempArr,stat=i_stat)
call memocc(i_stat,i_all,'tempArr',subname)

i_all=-product(shape(evall))*kind(evall)
deallocate(evall,stat=i_stat)
call memocc(i_stat,i_all,'evall',subname)

END SUBROUTINE loewdin

!>  This subroutine calculates the overlap matrix for a given bunch of orbitals. It also takes into 
!!  account k-points and spin.
!!
!!  Input arguments:
!!   @param  iproc      process ID
!!   @param  nproc      total number of processes
!!   @param  nspin      closed shell -> nspin=1 ; spin polarised -> nspin=2
!!   @param  norbIn     number of orbitals to be orthonormalized
!!   @param  orbs       type that contains many parameters concerning the orbitals
!!   @param  comms      type containing parameters for communicating the wavefunstion between processors
!!   @param  ndim_ovrlp      describes the shape of the overlap matrix
!!   @param  norbTot    total number of orbitals (if nspin=2:
!!               - norbTot(1)=total number of up orbitals
!!               - norbTot(2)=total number of down orbitals)
!!   @param  block1     gives the starting orbital of the orbitals to be orthonormalized
!!   @param  ispinIn    indicates whether the up or down orbitals shall be handled
!!   @param  category   gives the category for the timing
!!  Output arguments:
!   @param  ovrlp      the overlap matrix of the orbitals given in psip


subroutine getOverlap(iproc,nproc,nspin,norbIn,orbs,comms,&
     psi,ndim_ovrlp,ovrlp,norbTot,block1,ispinIn,category)

  use module_base
  use module_types
  use communications_base, only: comms_cubic
  implicit none

  ! Calling arguments
  character(len=*), intent(in) :: category
  integer,intent(in):: iproc,nproc,nspin,norbIn,block1,ispinIn
  type(orbitals_data),intent(in):: orbs
  type(comms_cubic),intent(in) :: comms
  real(wp),dimension(comms%nvctr_par(iproc,0)*orbs%nspinor*orbs%norb),intent(in) :: psi
  integer,dimension(nspin,0:orbs%nkpts),intent(in):: ndim_ovrlp
  real(wp),dimension(ndim_ovrlp(nspin,orbs%nkpts)),intent(out):: ovrlp
  integer,dimension(nspin),intent(in):: norbTot

  ! Local variables
  integer:: ispsi,ikptp,ikpt,ispin,nspinor,ncomp,norbs,ierr,nvctrp,norb



  ! Set the whole overlap matrix to zero. This is necessary since each process treats only a part
  ! of the matrix.
  call to_zero(ndim_ovrlp(nspin,orbs%nkpts),ovrlp)


  ispsi=1
  ! First make a loop over the k points handled by this process.
  do ikptp=1,orbs%nkptsp
     ! ikpt is the index of the k point.
     ikpt=orbs%iskpts+ikptp

     ! Now make also a loop over spin up/down.
     do ispin=1,nspin

        ! This subroutine gives essentially back nvctrp, i.e. the length of the vectors for which the overlap
        ! matrix shall be calculated. In addition it sets the value of nspinor to orbs%nspinor.
        call orbitals_and_components(iproc,ikpt,ispin,orbs,comms,&
             nvctrp,norb,norbs,ncomp,nspinor)
        ! The subroutine also overwrite the variable norb with the total number of orbitals.
        ! However we want to keep the value of norbIn (since we treat only a part of the orbitals).
        norb=norbIn
        ! Put the starting index to the right place. The current block of vector starts at the block1-th vector.
        ispsi=ispsi+nvctrp*(block1-1)*nspinor
        if(ispin==ispinIn) then
           if (nvctrp == 0) cycle

           ! Now calclulate one part of the overlap matrix. The starting index of this part is given by ndim_ovrlp(ispin,ikpt-1)+1.
           if(nspinor==1) then
              call syrk('L','T',norb,nvctrp,1.0_wp,psi(ispsi),max(1,nvctrp),&
                   0.0_wp,ovrlp(ndim_ovrlp(ispin,ikpt-1)+1),norb)
           else
              call herk('L','C',norb,ncomp*nvctrp,1.0_wp,psi(ispsi),&
                   max(1,ncomp*nvctrp),0.0_wp,ovrlp(ndim_ovrlp(ispin,ikpt-1)+1),norb)
           end if
        end if
        ! Move the starting indices to the end of the actual k point. This is necessary since nvctrp is
        ! different for the next k point and we cannot jump directly to the starting indices of our block for 
        ! the next k point.
        ispsi=ispsi+nvctrp*(norbTot(ispin)-block1+1)*nspinor
     end do
  end do

  !call MPI_BARRIER(bigdft_mpi%mpi_comm,ierr)
  !print *,'here',iproc

  if (nproc > 1) then
     !call timing(iproc,'GramS_comput  ','OF')
     !call timing(iproc,'GramS_commun  ','ON')
     call timing(iproc, trim(category)//'_comput', 'OF')
     call timing(iproc, trim(category)//'_commun', 'ON')
     call mpiallred(ovrlp(1),ndim_ovrlp(nspin,orbs%nkpts),MPI_SUM,bigdft_mpi%mpi_comm,ierr)
     call timing(iproc, trim(category)//'_commun', 'OF')
     call timing(iproc, trim(category)//'_comput', 'ON')
     !call timing(iproc,'GramS_commun  ','OF')
     !call timing(iproc,'GramS_comput  ','ON')
  end if

  ! Now each processors knows all the overlap matrices for each k-point
  ! even if it does not handle it.
  ! This is somehow redundant but it is one way of reducing the number of communications
  ! without defining group of processors.

END SUBROUTINE getOverlap

!>  This subroutine calculates the overlap matrix for a given bunch of orbitals. It also takes into 
!!  account k-points and spin.
!!
!!  Input arguments:
!!   @param  iproc      process ID
!!   @param  nproc      total number of processes
!!   @param  nspin      closed shell -> nspin=1 ; spin polarised -> nspin=2
!!   @param  norbIn     number of orbitals to be orthonormalized
!!   @param  orbs       type that contains many parameters concerning the orbitals
!!   @param  comms      type containing parameters for communicating the wavefunstion between processors
!!   @param  ndim_ovrlp      describes the shape of the overlap matrix
!!   @param  norbTot    total number of orbitals (if nspin=2:
!!               - norbTot(1)=total number of up orbitals
!!               - norbTot(2)=total number of down orbitals)
!!   @param  block1     gives the starting orbital of the orbitals to be orthonormalized
!!   @param  ispinIn    indicates whether the up or down orbitals shall be handled
!!   @param  category   gives the category for the timing
!!  Output arguments:
!!   @param  ovrlp      the overlap matrix of the orbitals given in psi
subroutine getOverlap_paw(iproc,nproc,nspin,norbIn,orbs,comms,&
     psi,spsi,ndim_ovrlp,ovrlp,norbTot,block1,ispinIn,category)

  use module_base
  use module_types
  use communications_base, only: comms_cubic
  implicit none

  ! Calling arguments
  character(len=*), intent(in) :: category
  integer,intent(in):: iproc,nproc,nspin,norbIn,block1,ispinIn
  type(orbitals_data),intent(in):: orbs
  type(comms_cubic),intent(in) :: comms
  real(wp),dimension(comms%nvctr_par(iproc,0)*orbs%nspinor*orbs%norb),intent(in) :: psi
  real(wp),dimension(comms%nvctr_par(iproc,0)*orbs%nspinor*orbs%norb),intent(in) :: spsi
  integer,dimension(nspin,0:orbs%nkpts),intent(in):: ndim_ovrlp
  real(wp),dimension(ndim_ovrlp(nspin,orbs%nkpts)),intent(out):: ovrlp
  integer,dimension(nspin),intent(in):: norbTot

  ! Local variables
  integer:: ispsi,ikptp,ikpt,ispin,nspinor,ncomp,norbs,ierr,nvctrp,norb
  real(wp),dimension(ndim_ovrlp(nspin,orbs%nkpts)):: ovrlp_pw


  ! Set the whole overlap matrix to zero. This is necessary since each process treats only a part
  ! of the matrix.
  call to_zero(ndim_ovrlp(nspin,orbs%nkpts),ovrlp)
  call to_zero(ndim_ovrlp(nspin,orbs%nkpts),ovrlp_pw)


  ispsi=1
  ! First make a loop over the k points handled by this process.
  do ikptp=1,orbs%nkptsp
     ! ikpt is the index of the k point.
     ikpt=orbs%iskpts+ikptp

     ! Now make also a loop over spin up/down.
     do ispin=1,nspin

        ! This subroutine gives essentially back nvctrp, i.e. the length of the vectors for which the overlap
        ! matrix shall be calculated. In addition it sets the value of nspinor to orbs%nspinor.
        call orbitals_and_components(iproc,ikpt,ispin,orbs,comms,&
             nvctrp,norb,norbs,ncomp,nspinor)
        ! The subroutine also overwrite the variable norb with the total number of orbitals.
        ! However we want to keep the value of norbIn (since we treat only a part of the orbitals).
        norb=norbIn
        ! Put the starting index to the right place. The current block of vector starts at the block1-th vector.
        ispsi=ispsi+nvctrp*(block1-1)*nspinor
        if(ispin==ispinIn) then
           if (nvctrp == 0) cycle

           ! Now calclulate one part of the overlap matrix. The starting index of this part is given by ndim_ovrlp(ispin,ikpt-1)+1.
           !Notice that two overlaps are computed:
           ! overlap_pw= <psi|psi>
           ! overlap   = <psi|S|psi>
           if(nspinor==1) then
              call syrk('L','T',norb,nvctrp,1.0_wp,psi(ispsi),max(1,nvctrp),&
                   0.0_wp,ovrlp_pw(ndim_ovrlp(ispin,ikpt-1)+1),norb)
              !for nspinor==1, ncomp==1
              call gemm('t','n',norb,norb,ncomp*nvctrp,1.0_wp,psi(ispsi),&
                   ncomp*nvctrp,spsi(ispsi),ncomp*nvctrp,0.d0,ovrlp(ndim_ovrlp(ispin,ikpt-1)+1),norb)
           else
              call herk('L','C',norb,ncomp*nvctrp,1.0_wp,psi(ispsi),&
                   max(1,ncomp*nvctrp),0.0_wp,ovrlp_pw(ndim_ovrlp(ispin,ikpt-1)+1),norb)
              !
              call c_gemm('c','n',norb,norb,ncomp*nvctrp,(1.0_wp,0.0_wp),psi(ispsi),&
                   ncomp*nvctrp,spsi(ispsi),ncomp*nvctrp,(0.d0,0.d0),ovrlp(ndim_ovrlp(ispin,ikpt-1)+1),norb)
           end if
        end if
        ! Move the starting indices to the end of the actual k point. This is necessary since nvctrp is
        ! different for the next k point and we cannot jump directly to the starting indices of our block for 
        ! the next k point.
        ispsi=ispsi+nvctrp*(norbTot(ispin)-block1+1)*nspinor
     end do
  end do

  !Sum the two overlaps:
  !overlap matrix in paw: O=1+S
  !<psi|O|psi> =  <psi|psi> + <psi|S|psi>
  ovrlp=ovrlp_pw + ovrlp

  !call MPI_BARRIER(MPI_COMM_WORLD,ierr)
  !print *,'here',iproc

  if (nproc > 1) then
     !call timing(iproc,'GramS_comput  ','OF')
     !call timing(iproc,'GramS_commun  ','ON')
     call timing(iproc, trim(category)//'_comput', 'OF')
     call timing(iproc, trim(category)//'_commun', 'ON')
     call mpiallred(ovrlp(1),ndim_ovrlp(nspin,orbs%nkpts),MPI_SUM,MPI_COMM_WORLD,ierr)
     !call MPI_ALLREDUCE (ovrlp(1,2),ovrlp(1,1),ndim_ovrlp(nspin,orbs%nkpts),mpidtypw,MPI_SUM,MPI_COMM_WORLD,ierr)
     call timing(iproc, trim(category)//'_commun', 'OF')
     call timing(iproc, trim(category)//'_comput', 'ON')
     !call timing(iproc,'GramS_commun  ','OF')
     !call timing(iproc,'GramS_comput  ','ON')
  end if

  ! Now each processors knows all the overlap matrices for each k-point
  ! even if it does not handle it.
  ! This is somehow redundant but it is one way of reducing the number of communications
  ! without defining group of processors.

END SUBROUTINE getOverlap_paw

!>  This subroutine calculates the overlap matrix for a given bunch of orbitals. It also takes into 
!!  account k-points and spin.
!!
!!  Input arguments:
!!   @param  iproc      process ID
!!   @param  nproc      total number of processes
!!   @param  nspin      closed shell -> nspin=1 ; spin polarised -> nspin=2
!!   @param  norbIn     number of orbitals to be orthonormalized
!!   @param  istart     second dimension of the overlpa matrix
!!   @param  orbs       type that contains many parameters concerning the orbitals
!!   @param  comms      type containing parameters for communicating the wavefunstion between processors
!!   @param  psit    the orbitals 
!!   @param  ndim_ovrlp  describes the shape of the overlap matrix
!!   @param  norbTot    total number of orbitals (if nspin=2:
!!               - norbTot(1)=total number of up orbitals
!!               - norbTot(2)=total number of down orbitals)
!!   @param  block1     gives the starting orbital of the orbitals to be orthogonalized
!!   @param  block2     gives the starting orbital of the orbitals to which the orbitals shall orthogonalized
!!   @param  ispinIn    indicates whether the up or down orbitals shall be handled
!!   @param  category   gives the category for the timing
!!  Output arguments:
!!   @param  ovrlp      the overlap matrix of the orbitals given in psi
subroutine getOverlapDifferentPsi(iproc, nproc, nspin, norbIn, orbs, comms,&
     psit, ndim_ovrlp, ovrlp, norbTot, block1, block2, ispinIn, category)

  use module_base
  use module_types
  use communications_base, only: comms_cubic
  implicit none

  ! Calling arguments
  !integer,intent(in):: iproc, nproc, nspin, norbIn,  istart, norbTot, block1, block2
  character(len=*), intent(in) :: category
  integer,intent(in):: iproc, nproc, nspin, norbIn, block1, block2, ispinIn
  type(orbitals_data),intent(in):: orbs
  type(comms_cubic),intent(in) :: comms
  real(kind=8),dimension(comms%nvctr_par(iproc,0)*orbs%nspinor*orbs%norb),intent(in) :: psit
  integer,dimension(nspin,0:orbs%nkpts),intent(in):: ndim_ovrlp
  real(kind=8),dimension(ndim_ovrlp(nspin,orbs%nkpts)):: ovrlp
  integer,dimension(nspin):: norbTot
  ! Local variables
  integer:: ikptp, ikpt, ispin, nspinor, ncomp, norbs, ierr, nvctrp, norb, ispsi1, ispsi2
  
  ! Set the whole overlap matrix to zero. This is necessary since each process treats only a part
  ! of the matrix.
  call to_zero(ndim_ovrlp(nspin,orbs%nkpts),ovrlp)

  ispsi1=1
  ispsi2=1
  ! First make a loop over the k points handled by this process.
  do ikptp=1,orbs%nkptsp
     ! ikpt is the index of the k point.
     ikpt=orbs%iskpts+ikptp
     
     ! Now make also a loop over spin up/down.
     do ispin=1,nspin
        
        ! This subroutine gives essentially back nvctrp, i.e. the length of the vectors for which the overlap
        ! matrix shall be calculated. In addition it sets the value of nspinor to orbs%nspinor.
        call orbitals_and_components(iproc,ikpt,ispin,orbs,comms,&
             nvctrp,norb,norbs,ncomp,nspinor)
        ! The subroutine also overwrite the variable norb with the total number of orbitals.
        ! However we want to keep the value of norbIn (since we treat only a part of the orbitals).
        norb=norbIn
        
        ! Put the starting index to the right place. The current block of vector starts at the block1-th and
        ! block2-th vector, respectively. 
        ispsi1=ispsi1+nvctrp*(block1-1)*nspinor
        ispsi2=ispsi2+nvctrp*(block2-1)*nspinor
        if(ispin==ispinIn) then
            if (nvctrp == 0) cycle
       
            ! Now calclulate one part of the overlap matrix. The starting index of this part is given by ndim_ovrlp(ispin,ikpt-1)+1.
            if(nspinor==1) then
               call gemm('t','n',norb,norb,ncomp*nvctrp,1.0_wp,psit(ispsi2),&
                    ncomp*nvctrp,psit(ispsi1),ncomp*nvctrp,0.d0,ovrlp(ndim_ovrlp(ispin,ikpt-1)+1),norb)
            else
               call c_gemm('c','n',norb,norb,ncomp*nvctrp,(1.0_wp,0.0_wp),psit(ispsi2),&
                    ncomp*nvctrp,psit(ispsi1),ncomp*nvctrp,(0.d0,0.d0),ovrlp(ndim_ovrlp(ispin,ikpt-1)+1),norb)
            end if

        end if
        ! Move the starting indices to the end of the actual k point. This is necessary since nvctrp is
        ! different for the next k point and we cannot jump directly to the starting indices of our block for 
        ! the next k point.
        ispsi1=ispsi1+nvctrp*(norbTot(ispin)-block1+1)*nspinor
        ispsi2=ispsi2+nvctrp*(norbTot(ispin)-block2+1)*nspinor

     end do
  end do

  ! Sum up the overlap matrices from all processes.
  if (nproc > 1) then
     !call timing(iproc,'GramS_comput  ','OF')
     !call timing(iproc,'GramS_commun  ','ON')
     call timing(iproc,trim(category)//'_comput','OF')
     call timing(iproc,trim(category)//'_commun','ON')
     call mpiallred(ovrlp(1),ndim_ovrlp(nspin,orbs%nkpts),MPI_SUM,bigdft_mpi%mpi_comm,ierr)
     !call mpi_allreduce(ovrlp(1,2),ovrlp(1,1),ndim_ovrlp(nspin,orbs%nkpts),mpi_double_precision,mpi_sum,bigdft_mpi%mpi_comm,ierr)
     call timing(iproc,trim(category)//'_commun','OF')
     call timing(iproc,trim(category)//'_comput','ON')
     !call timing(iproc,'GramS_commun  ','OF')
     !call timing(iproc,'GramS_comput  ','ON')
  end if
  
  ! Now each processors knows all the overlap matrices for each k-point even if it does not handle it.
  
END SUBROUTINE getOverlapDifferentPsi

!>  This subroutine calculates the overlap matrix for a given bunch of orbitals. It also takes into 
!!  account k-points and spin.
!!
!!  Input arguments:
!!   @param  iproc      process ID
!!   @param  nproc      total number of processes
!!   @param  nspin      closed shell -> nspin=1 ; spin polarised -> nspin=2
!!   @param  norbIn     number of orbitals to be orthonormalized
!!   @param  istart     second dimension of the overlpa matrix
!!   @param  orbs       type that contains many parameters concerning the orbitals
!!   @param  comms      type containing parameters for communicating the wavefunstion between processors
!!   @param  psit    the orbitals 
!!   @param  ndim_ovrlp  describes the shape of the overlap matrix
!!   @param  norbTot    total number of orbitals (if nspin=2:
!!               - norbTot(1)=total number of up orbitals
!!               - norbTot(2)=total number of down orbitals)
!!   @param  block1     gives the starting orbital of the orbitals to be orthogonalized
!!   @param  block2     gives the starting orbital of the orbitals to which the orbitals shall orthogonalized
!!   @param  ispinIn    indicates whether the up or down orbitals shall be handled
!!   @param  category   gives the category for the timing
!!  Output arguments:
!!   @param  ovrlp      the overlap matrix of the orbitals given in psi
subroutine getOverlapDifferentPsi_paw(iproc, nproc, nspin, norbIn, orbs, comms,&
     psit, spsit,ndim_ovrlp, ovrlp, norbTot, block1, block2, ispinIn, category)

  use module_base
  use module_types
  use communications_base, only: comms_cubic
  implicit none

  ! Calling arguments
  !integer,intent(in):: iproc, nproc, nspin, norbIn,  istart, norbTot, block1, block2
  character(len=*), intent(in) :: category
  integer,intent(in):: iproc, nproc, nspin, norbIn, block1, block2, ispinIn
  type(orbitals_data),intent(in):: orbs
  type(comms_cubic),intent(in) :: comms
  real(kind=8),dimension(comms%nvctr_par(iproc,0)*orbs%nspinor*orbs%norb),intent(in) :: psit,spsit
  integer,dimension(nspin,0:orbs%nkpts),intent(in):: ndim_ovrlp
  real(kind=8),dimension(ndim_ovrlp(nspin,orbs%nkpts)):: ovrlp
  integer,dimension(nspin):: norbTot
  ! Local variables
  integer:: ikptp, ikpt, ispin, nspinor, ncomp, norbs, ierr, nvctrp, norb, ispsi1, ispsi2
  real(kind=8),dimension(ndim_ovrlp(nspin,orbs%nkpts)):: ovrlp_pw
  
  ! Set the whole overlap matrix to zero. This is necessary since each process treats only a part
  ! of the matrix.
  call to_zero(ndim_ovrlp(nspin,orbs%nkpts),ovrlp)
  call to_zero(ndim_ovrlp(nspin,orbs%nkpts),ovrlp_pw)

  ispsi1=1
  ispsi2=1
  ! First make a loop over the k points handled by this process.
  do ikptp=1,orbs%nkptsp
     ! ikpt is the index of the k point.
     ikpt=orbs%iskpts+ikptp
     
     ! Now make also a loop over spin up/down.
     do ispin=1,nspin
        
        ! This subroutine gives essentially back nvctrp, i.e. the length of the vectors for which the overlap
        ! matrix shall be calculated. In addition it sets the value of nspinor to orbs%nspinor.
        call orbitals_and_components(iproc,ikpt,ispin,orbs,comms,&
             nvctrp,norb,norbs,ncomp,nspinor)
        ! The subroutine also overwrite the variable norb with the total number of orbitals.
        ! However we want to keep the value of norbIn (since we treat only a part of the orbitals).
        norb=norbIn
        
        ! Put the starting index to the right place. The current block of vector starts at the block1-th and
        ! block2-th vector, respectively. 
        ispsi1=ispsi1+nvctrp*(block1-1)*nspinor
        ispsi2=ispsi2+nvctrp*(block2-1)*nspinor
        if(ispin==ispinIn) then
            if (nvctrp == 0) cycle
       
            ! Now calclulate one part of the overlap matrix. The starting index of this part is given by ndim_ovrlp(ispin,ikpt-1)+1.
            !
            !Notice that two overlaps are computed:
            ! overlap_pw= <psi|psi>
            ! overlap   = <psi|S|psi>
            if(nspinor==1) then
               call gemm('t','n',norb,norb,ncomp*nvctrp,1.0_wp,psit(ispsi2),&
                    ncomp*nvctrp,spsit(ispsi1),ncomp*nvctrp,0.d0,ovrlp(ndim_ovrlp(ispin,ikpt-1)+1),norb)
               call gemm('t','n',norb,norb,ncomp*nvctrp,1.0_wp,psit(ispsi2),&
                    ncomp*nvctrp,psit(ispsi1),ncomp*nvctrp,0.d0,ovrlp_pw(ndim_ovrlp(ispin,ikpt-1)+1),norb)
            else
               call c_gemm('c','n',norb,norb,ncomp*nvctrp,(1.0_wp,0.0_wp),psit(ispsi2),&
                    ncomp*nvctrp,spsit(ispsi1),ncomp*nvctrp,(0.d0,0.d0),ovrlp(ndim_ovrlp(ispin,ikpt-1)+1),norb)
               call c_gemm('c','n',norb,norb,ncomp*nvctrp,(1.0_wp,0.0_wp),psit(ispsi2),&
                    ncomp*nvctrp,psit(ispsi1),ncomp*nvctrp,(0.d0,0.d0),ovrlp_pw(ndim_ovrlp(ispin,ikpt-1)+1),norb)
            end if

        end if
        ! Move the starting indices to the end of the actual k point. This is necessary since nvctrp is
        ! different for the next k point and we cannot jump directly to the starting indices of our block for 
        ! the next k point.
        ispsi1=ispsi1+nvctrp*(norbTot(ispin)-block1+1)*nspinor
        ispsi2=ispsi2+nvctrp*(norbTot(ispin)-block2+1)*nspinor

     end do
  end do

  !Sum the two overlaps:
  !overlap matrix in paw: O=1+S
  !<psi|O|psi> =  <psi|psi> + <psi|S|psi>
  ovrlp=ovrlp_pw + ovrlp


  ! Sum up the overlap matrices from all processes.
  if (nproc > 1) then
     !call timing(iproc,'GramS_comput  ','OF')
     !call timing(iproc,'GramS_commun  ','ON')
     call timing(iproc,trim(category)//'_comput','OF')
     call timing(iproc,trim(category)//'_commun','ON')
     call mpiallred(ovrlp(1),ndim_ovrlp(nspin,orbs%nkpts),MPI_SUM,MPI_COMM_WORLD,ierr)
     !call mpi_allreduce(ovrlp(1,2),ovrlp(1,1),ndim_ovrlp(nspin,orbs%nkpts),mpi_double_precision,mpi_sum,mpi_comm_world,ierr)
     call timing(iproc,trim(category)//'_commun','OF')
     call timing(iproc,trim(category)//'_comput','ON')
     !call timing(iproc,'GramS_commun  ','OF')
     !call timing(iproc,'GramS_comput  ','ON')
  end if
  
  ! Now each processors knows all the overlap matrices for each k-point even if it does not handle it.
  
END SUBROUTINE getOverlapDifferentPsi_paw




subroutine dimension_ovrlpFixedNorb(nspin,orbs,ndim_ovrlp,norb)
  use module_base
  use module_types
  implicit none
  integer, intent(in) :: nspin,norb
  type(orbitals_data), intent(in) :: orbs
  integer, dimension(nspin,0:orbs%nkpts), intent(inout) :: ndim_ovrlp
  !local variables
  integer :: norbs,ncomp,ikpt

  ndim_ovrlp(1,0)=0
  if (nspin == 2) then
     !norb=orbs%norbu

     !this is first k-point
     call complex_components(orbs%nspinor,norb,norbs,ncomp)

     ndim_ovrlp(2,0)=norbs*norb
  end if

  do ikpt=1,orbs%nkpts
     !this part should be enhanced for real k-points
     !norb=orbs%norbu
     !if (nspin == 2) norb = orbs%norbd
     !this is ikpt k-point
     call complex_components(orbs%nspinor,norb,norbs,ncomp)

     ndim_ovrlp(1,ikpt)=ndim_ovrlp(nspin,ikpt-1)+norbs*norb
     if (orbs%norbd > 0) then
        !norb=orbs%norbu
        !this is ikpt+1
        call complex_components(orbs%nspinor,norb,norbs,ncomp)
        if (ikpt == orbs%nkpts) then
           ndim_ovrlp(2,ikpt)=ndim_ovrlp(1,ikpt)
        else
           ndim_ovrlp(2,ikpt)=ndim_ovrlp(1,ikpt)+norbs*norb
        end if
     end if
  end do

END SUBROUTINE dimension_ovrlpFixedNorb







!!****f* BigDFT/orthoconstraint
!! FUNCTION
!!   Orthogonality routine, for all the orbitals
!!   Uses wavefunctions in their transposed form
!! SOURCE
!!
!!subroutine orthoconstraintNotSymmetric(iproc,nproc,ndim,orbs,comms,wfd,psi,hpsi,scprsum,lagMatDiag)
!!  use module_base
!!  use module_types
!!  implicit none
!!  integer, intent(in) :: iproc,nproc,ndim
!!  type(orbitals_data), intent(in) :: orbs
!!  type(comms_cubic), intent(in) :: comms
!!  type(wavefunctions_descriptors), intent(in) :: wfd
!!  !real(wp), dimension(sum(comms%nvctr_par(iproc,1:orbs%nkptsp))*orbs%nspinor*orbs%norb), intent(in) :: psi
!!  !real(wp), dimension(sum(comms%nvctr_par(iproc,1:orbs%nkptsp))*orbs%nspinor*orbs%norb), intent(out) :: hpsi
!!  real(wp), dimension(ndim*orbs%nspinor*orbs%norb), intent(in) :: psi
!!  real(wp), dimension(ndim*orbs%nspinor*orbs%norb), intent(out) :: hpsi
!!  real(dp), intent(out) :: scprsum
!!  real(dp),dimension(orbs%norb),intent(out):: lagMatDiag
!!  !local variables
!!  character(len=*), parameter :: subname='orthoconstraintNotSymmetric'
!!  integer :: i_stat,i_all,ierr,iorb,ise,jorb
!!  integer :: ispin,nspin,ikpt,norb,norbs,ncomp,nvctrp,ispsi,ikptp,nspinor
!!  real(dp) :: occ,tt
!!  integer, dimension(:,:), allocatable :: ndim_ovrlp
!!  real(wp), dimension(:), allocatable :: alag
!!
!!
!!integer:: istart, jstart
!!
!!
!!  !separate the orthogonalisation procedure for up and down orbitals 
!!  !and for different k-points
!!  call timing(iproc,'LagrM_comput  ','ON')
!!
!!  !number of components of the overlap matrix for parallel case
!!  !calculate the dimension of the overlap matrix for each k-point
!!  if (orbs%norbd > 0) then
!!     nspin=2
!!  else
!!     nspin=1
!!  end if
!!
!!  !number of components for the overlap matrix in wp-kind real numbers
!!
!!  allocate(ndim_ovrlp(nspin,0:orbs%nkpts+ndebug),stat=i_stat)
!!  call memocc(i_stat,ndim_ovrlp,'ndim_ovrlp',subname)
!!
!!  call dimension_ovrlp(nspin,orbs,ndim_ovrlp)
!!
!!  allocate(alag(ndim_ovrlp(nspin,orbs%nkpts)+ndebug),stat=i_stat)
!!  call memocc(i_stat,alag,'alag',subname)
!!
!!  !put to zero all the k-points which are not needed
!!  call to_zero(ndim_ovrlp(nspin,orbs%nkpts),alag)
!!
!!  !do it for each of the k-points and separate also between up and down orbitals in the non-collinear case
!!  ispsi=1
!!  do ikptp=1,orbs%nkptsp
!!     ikpt=orbs%iskpts+ikptp!orbs%ikptsp(ikptp)
!!
!!     do ispin=1,nspin
!!
!!        call orbitals_and_components(iproc,ikptp,ispin,orbs,comms,&
!!             nvctrp,norb,norbs,ncomp,nspinor)
!!        if (nvctrp == 0) cycle
!!
!!        if(nspinor==1) then
!!           call gemm('T','N',norb,norb,nvctrp,1.0_wp,psi(ispsi),&
!!                max(1,nvctrp),hpsi(ispsi),max(1,nvctrp),0.0_wp,&
!!                alag(ndim_ovrlp(ispin,ikpt-1)+1),norb)
!!        else
!!           !this part should be recheck in the case of nspinor == 2
!!           call c_gemm('C','N',norb,norb,ncomp*nvctrp,(1.0_wp,0.0_wp),psi(ispsi),&
!!                max(1,ncomp*nvctrp), &
!!                hpsi(ispsi),max(1,ncomp*nvctrp),(0.0_wp,0.0_wp),&
!!                alag(ndim_ovrlp(ispin,ikpt-1)+1),norb)
!!        end if
!!        ispsi=ispsi+nvctrp*norb*nspinor
!!     end do
!!  end do
!!
!!  if (nproc > 1) then
!!     call timing(iproc,'LagrM_comput  ','OF')
!!     call timing(iproc,'LagrM_commun  ','ON')
!!     call mpiallred(alag(1),ndim_ovrlp(nspin,orbs%nkpts),MPI_SUM,bigdft_mpi%mpi_comm,ierr)
!!     call timing(iproc,'LagrM_commun  ','OF')
!!     call timing(iproc,'LagrM_comput  ','ON')
!!  end if
!!
!!  ! Copy the diagonal of the matrix
!!  do iorb=1,orbs%norb
!!      lagMatDiag(iorb)=alag((iorb-1)*orbs%norb+iorb)
!!  end do
!!! Lagrange multiplier matrix
!!!!if(iproc==0) write(*,*) 'Lagrange multiplier matrix'
!!!!do iorb=1,norb
!!!!    !if(iproc==0) write(*,'(80f8.4)') (alag(norb*jorb+iorb), jorb=0,norb-1)
!!!!    do jorb=1,norb
!!!!        write(1100+iproc,*) iorb, jorb, alag(norb*(iorb-1)+jorb)
!!!!    end do
!!!!end do
!!
!!
!!  !now each processors knows all the overlap matrices for each k-point
!!  !even if it does not handle it.
!!  !this is somehow redundant but it is one way of reducing the number of communications
!!  !without defining group of processors
!!
!!  !calculate the sum of the diagonal of the overlap matrix, for each k-point
!!  scprsum=0.0_dp
!!  !for each k-point calculate the gradient
!!  ispsi=1
!!  do ikptp=1,orbs%nkptsp
!!     ikpt=orbs%iskpts+ikptp!orbs%ikptsp(ikptp)
!!
!!     do ispin=1,nspin
!!        if (ispin==1) ise=0
!!        call orbitals_and_components(iproc,ikptp,ispin,orbs,comms,&
!!             nvctrp,norb,norbs,ncomp,nspinor)
!!        if (nvctrp == 0) cycle
!!
!!!!$        !correct the orthogonality constraint if there are some orbitals which have zero occupation number
!!!!$        do iorb=1,norb
!!!!$           do jorb=iorb+1,norb
!!!!$              if (orbs%occup((ikpt-1)*orbs%norb+iorb+ise) /= 0.0_gp .and. &
!!!!$                   orbs%occup((ikpt-1)*orbs%norb+jorb+ise) == 0.0_gp) then
!!!!$                 alag(ndim_ovrlp(ispin,ikpt-1)+iorb+(jorb-1)*norbs) = 0.0_wp
!!!!$                 alag(ndim_ovrlp(ispin,ikpt-1)+jorb+(iorb-1)*norbs) = 0.0_wp
!!!!$                 !if (iproc ==0) print *,'i,j',iorb,jorb,alag(ndim_ovrlp(ispin,ikpt-1)+iorb+(jorb-1)*norbs)
!!!!$              end if
!!!!$           end do
!!!!$        end do
!!
!!        !calculate the scprsum if the k-point is associated to this processor
!!        !the scprsum always coincide with the trace of the hamiltonian
!!        if (orbs%ikptproc(ikpt) == iproc) then
!!           occ=real(orbs%kwgts(ikpt),dp)
!!           if(nspinor == 1) then
!!              do iorb=1,norb
!!                 scprsum=scprsum+occ*real(alag(ndim_ovrlp(ispin,ikpt-1)+iorb+(iorb-1)*norbs),dp)
!!              enddo
!!           else if (nspinor == 4 .or. nspinor == 2) then
!!              !not sure about the imaginary part of the diagonal (should be zero if H is hermitian)
!!              do iorb=1,norb
!!                 scprsum=scprsum+&
!!                      occ*real(alag(ndim_ovrlp(ispin,ikpt-1)+2*iorb-1+(iorb-1)*norbs),dp)
!!                 scprsum=scprsum+&
!!                      occ*real(alag(ndim_ovrlp(ispin,ikpt-1)+2*iorb+(iorb-1)*norbs),dp)
!!              enddo
!!           end if
!!        end if
!!        ise=norb
!!
!!        if(nspinor==1 .and. nvctrp /= 0) then
!!           !call gemm('N','N',nvctrp,norb,norb,-1.0_wp,psi(ispsi),max(1,nvctrp),&
!!           !     alag(ndim_ovrlp(ispin,ikpt-1)+1),norb,1.0_wp,&
!!           !     hpsi(ispsi),max(1,nvctrp))
!!           call gemm('N','N',nvctrp,norb,norb,-.5_wp,psi(ispsi),max(1,nvctrp),&
!!                alag(ndim_ovrlp(ispin,ikpt-1)+1),norb,1.0_wp,&
!!                hpsi(ispsi),max(1,nvctrp))
!!           call gemm('N','T',nvctrp,norb,norb,-.5_wp,psi(ispsi),max(1,nvctrp),&
!!                alag(ndim_ovrlp(ispin,ikpt-1)+1),norb,1.0_wp,&
!!                hpsi(ispsi),max(1,nvctrp))
!!        else if (nvctrp /= 0) then
!!           stop 'not implemented for nspinor/=1!'
!!           call c_gemm('N','N',ncomp*nvctrp,norb,norb,(-1.0_wp,0.0_wp),psi(ispsi),max(1,ncomp*nvctrp),&
!!                alag(ndim_ovrlp(ispin,ikpt-1)+1),norb,(1.0_wp,0.0_wp),hpsi(ispsi),max(1,ncomp*nvctrp))
!!        end if
!!        ispsi=ispsi+nvctrp*norb*nspinor
!!     end do
!!  end do
!!
!!  if (nproc > 1) then
!!     tt=scprsum
!!     call mpiallred(scprsum,1,MPI_SUM,bigdft_mpi%mpi_comm,ierr)
!!     !call MPI_ALLREDUCE(tt,scprsum,1,mpidtypd,MPI_SUM,bigdft_mpi%mpi_comm,ierr)
!!  end if
!!
!!  i_all=-product(shape(alag))*kind(alag)
!!  deallocate(alag,stat=i_stat)
!!  call memocc(i_stat,i_all,'alag',subname)
!!
!!  i_all=-product(shape(ndim_ovrlp))*kind(ndim_ovrlp)
!!  deallocate(ndim_ovrlp,stat=i_stat)
!!  call memocc(i_stat,i_all,'ndim_ovrlp',subname)
!!
!!  call timing(iproc,'LagrM_comput  ','OF')
!!
!!END SUBROUTINE orthoconstraintNotSymmetric
!!***







!!!!!!! FUNCTION
!!!!!!!    Orthogonality routine, for all the orbitals
!!!!!!!    Uses wavefunctions in their transposed form
!!!!!!!
!!!!!!! COPYRIGHT
!!!!!!!    Copyright (C) 2007-2010 BigDFT group
!!!!!!!    This file is distributed under the terms of the
!!!!!!!    GNU General Public License, see ~/COPYING file
!!!!!!!    or http://www.gnu.org/copyleft/gpl.txt .
!!!!!!!    For the list of contributors, see ~/AUTHORS 
!!!!!!!
!!!!!!! SOURCE
!!!!!!!
!!!!!!subroutine orthogonalizeLIN(iproc,nproc,orbs,comms,wfd,psi,input)
!!!!!subroutine orthogonalizeLIN(iproc, lproc, uproc, norbPerGroup, orbs, comms, psi, input, newComm)
!!!!!  use module_base
!!!!!  use module_types
!!!!!  implicit none
!!!!!  integer, intent(in) :: iproc,lproc, uproc, norbPerGroup, newComm
!!!!!  type(orbitals_data), intent(in) :: orbs
!!!!!  type(comms_cubic), intent(in) :: comms
!!!!!  type(input_variables), intent(in) :: input
!!!!!  !real(wp), dimension(sum(comms%nvctr_par(iproc,1:orbs%nkptsp))*orbs%nspinor*orbs%norb), intent(inout) :: psi
!!!!!  real(wp), dimension(orbs%npsidim), intent(inout) :: psi
!!!!!  !local variables
!!!!!  character(len=*), parameter :: subname='orthogonalize'
!!!!!  integer :: i_stat,i_all,ierr,info
!!!!!  integer :: ispin,nspin,ikpt,norb,norbs,ncomp,nvctrp,ispsi,ikptp,nspinor
!!!!!  integer, dimension(:,:), allocatable :: ndim_ovrlp
!!!!!  real(wp), dimension(:), allocatable :: ovrlp
!!!!!  integer,dimension(:),allocatable:: norbArr
!!!!!  character(len=20):: category
!!!!!integer:: nproc
!!!!!
!!!!!  nproc=uproc-lproc+1
!!!!!
!!!!!  ! Determine wheter we have close shell (nspin=1) or spin polarized (nspin=2)
!!!!!  if (orbs%norbd>0) then 
!!!!!     nspin=2 
!!!!!  else 
!!!!!     nspin=1 
!!!!!  end if
!!!!!
!!!!!  ! ndim_ovrlp describes the shape of the overlap matrix.
!!!!!  allocate(ndim_ovrlp(nspin,0:orbs%nkpts+ndebug),stat=i_stat)
!!!!!  call memocc(i_stat,ndim_ovrlp,'ndim_ovrlp',subname)
!!!!!  
!!!!!  ! Allocate norbArr which contains the number of up and down orbitals.
!!!!!  allocate(norbArr(nspin), stat=i_stat)
!!!!!  call memocc(i_stat,norbArr,'norbArr',subname)
!!!!!  do ispin=1,nspin
!!!!!     !if(ispin==1) norbArr(ispin)=orbs%norbu
!!!!!     !if(ispin==2) norbArr(ispin)=orbs%norbd
!!!!!     if(ispin==1) norbArr(ispin)=norbPerGroup
!!!!!     if(ispin==2) stop 'ERROR: not yet implemented for ispin==2!!'
!!!!!  end do
!!!!!
!!!!!  ! Choose which orthogonalization method shall be used:
!!!!!  ! methOrtho==0: Cholesky orthonormalization (i.e. a pseudo Gram-Schmidt)
!!!!!  ! methOrtho==1: hybrid Gram-Schmidt/Cholesky orthonormalization
!!!!!  ! methOrtho==2: Loewdin orthonormalization
!!!!!  if(input%methOrtho==0) then
!!!!!     category='Chol'
!!!!!     call timing(iproc, trim(category)//'_comput', 'ON')
!!!!!
!!!!!     !call dimension_ovrlp(nspin,orbs,ndim_ovrlp)
!!!!!     call dimension_ovrlpFixedNorb(nspin,orbs,ndim_ovrlp,norbPerGroup)
!!!!!
!!!!!     ! Allocate the overlap matrix
!!!!!     allocate(ovrlp(ndim_ovrlp(nspin,orbs%nkpts)+ndebug),stat=i_stat)
!!!!!     call memocc(i_stat,ovrlp,'ovrlp',subname)
!!!!!
!!!!!     !print *,'there',iproc
!!!!!
!!!!!     ! Make a loop over npsin; calculate the overlap matrix (for up/down, resp.) and orthogonalize (again for up/down, resp.).
!!!!!     do ispin=1,nspin
!!!!!        !call getOverlap(iproc,nproc,nspin,norbArr(ispin),orbs,comms,&
!!!!!        !     psi(1),ndim_ovrlp,ovrlp,norbArr,1,ispin,category)
!!!!!        call getOverlapLIN(iproc,nproc,nspin,norbArr(ispin),orbs,comms,&
!!!!!             psi(1),ndim_ovrlp,ovrlp,norbArr,1,ispin,category, newComm)
!!!!!        call cholesky(iproc,nproc,norbArr(ispin),psi(1),nspinor,nspin,orbs,comms,&
!!!!!             ndim_ovrlp,ovrlp(1),norbArr,1,ispin)
!!!!!!call cholesky(iproc,nproc,norbArr(ispin),psi(1),nspinor,nspin,orbs,comms,&
!!!!!!     ndim_ovrlp,ovrlp(1),norbArr,1,ispin)
!!!!!!do i_stat=1,size(ovrlp)
!!!!!!    write(2000+iproc,*) i_stat, ovrlp(i_stat)
!!!!!!end do
!!!!!     end do
!!!!!
!!!!!     ! Deallocate the arrays.
!!!!!     i_all=-product(shape(ovrlp))*kind(ovrlp)
!!!!!     deallocate(ovrlp,stat=i_stat)
!!!!!     call memocc(i_stat,i_all,'ovrlp',subname)
!!!!!
!!!!!  else if(input%methOrtho==1) then
!!!!!       category='GS/Chol'
!!!!!       call timing(iproc, trim(category)//'_comput', 'ON')
!!!!!       
!!!!!       ! Make a hybrid Gram-Schmidt/Cholesky orthonormalization.
!!!!!       call gsChol(iproc,nproc,psi(1),input,nspinor,orbs,nspin,ndim_ovrlp,norbArr,comms)
!!!!!  else if(input%methOrtho==2) then
!!!!!     category='Loewdin'
!!!!!     call timing(iproc,trim(category)//'_comput','ON')
!!!!!
!!!!!     call dimension_ovrlp(nspin,orbs,ndim_ovrlp)
!!!!!     
!!!!!     ! Allocate the overlap matrix
!!!!!     allocate(ovrlp(ndim_ovrlp(nspin,orbs%nkpts)+ndebug),stat=i_stat)
!!!!!     call memocc(i_stat,ovrlp,'ovrlp',subname)
!!!!!          
!!!!!     ! Make a loop over npsin; calculate the overlap matrix (for up/down,resp.) and orthogonalize (again for up/down,resp.).
!!!!!     do ispin=1,nspin
!!!!!        call getOverlap(iproc,nproc,nspin,norbArr(ispin),orbs,comms,psi(1),ndim_ovrlp,ovrlp,norbArr,1,ispin,category)
!!!!!        call loewdin(iproc,nproc,norbArr(ispin),orbs%nspinor,1,ispin,orbs,comms,nspin,psi,ovrlp,ndim_ovrlp,norbArr)
!!!!!     end do
!!!!!     
!!!!!     ! Deallocate the arrays.
!!!!!     i_all=-product(shape(ovrlp))*kind(ovrlp)
!!!!!     deallocate(ovrlp,stat=i_stat)
!!!!!     call memocc(i_stat,i_all,'ovrlp',subname)
!!!!!          
!!!!!  else
!!!!!     if(iproc==0) write(*,'(a)') 'ERROR: invalid choice for methOrtho.'
!!!!!     if(iproc==0) write(*,'(a)') "Change it in 'input.perf' to 0, 1 or 2!"
!!!!!     stop
!!!!!  end if
!!!!!
!!!!!  ! Deallocate the remaining arrays.
!!!!!  i_all=-product(shape(norbArr))*kind(norbArr)
!!!!!  deallocate(norbArr, stat=i_stat)
!!!!!  call memocc(i_stat,i_all,'norbArr',subname)
!!!!!
!!!!!  i_all=-product(shape(ndim_ovrlp))*kind(ndim_ovrlp)
!!!!!  deallocate(ndim_ovrlp, stat=i_stat)
!!!!!  call memocc(i_stat,i_all,'ndim_ovrlp',subname)
!!!!!
!!!!!
!!!!!  call timing(iproc,trim(category)//'_comput','OF')
!!!!!  
!!!!!END SUBROUTINE orthogonalizeLIN
!!!!!
!!!!!
!!!!!subroutine getOverlapLIN(iproc,nproc,nspin,norbIn,orbs,comms,&
!!!!!     psi,ndim_ovrlp,ovrlp,norbTot,block1,ispinIn,category, newComm)
!!!!!  !
!!!!!  ! Purpose:
!!!!!  ! =======
!!!!!  !  This subroutine calculates the overlap matrix for a given bunch of orbitals. It also takes into 
!!!!!  !  account k-points and spin.
!!!!!  !
!!!!!  ! Calling arguments:
!!!!!  ! =================
!!!!!  !  Input arguments:
!!!!!  !    iproc      process ID
!!!!!  !    nproc      total number of processes
!!!!!  !    nspin      closed shell -> nspin=1 ; spin polarised -> nspin=2
!!!!!  !    norbIn     number of orbitals to be orthonormalized
!!!!!  !    orbs       type that contains many parameters concerning the orbitals
!!!!!  !    comms      type containing parameters for communicating the wavefunstion between processors
!!!!!  !    ndim_ovrlp      describes the shape of the overlap matrix
!!!!!  !    norbTot    total number of orbitals (if nspin=2:
!!!!!  !                 norbTot(1)=total number of up orbitals
!!!!!  !                 norbTot(2)=total number of down orbitals)
!!!!!  !    block1     gives the starting orbital of the orbitals to be orthonormalized
!!!!!  !    ispinIn    indicates whether the up or down orbitals shall be handled
!!!!!  !    catgeory   gives the category for the timing
!!!!!  !  Output arguments:
!!!!!  !    ovrlp      the overlap matrix of the orbitals given in psi
!!!!!  !
!!!!!  use module_base
!!!!!  use module_types
!!!!!  implicit none
!!!!!
!!!!!  ! Calling arguments
!!!!!  character(len=*), intent(in) :: category
!!!!!  integer,intent(in):: iproc,nproc,nspin,norbIn,block1,ispinIn, newComm
!!!!!  type(orbitals_data),intent(in):: orbs
!!!!!  type(comms_cubic),intent(in) :: comms
!!!!!  real(wp),dimension(sum(comms%nvctr_par(iproc,1:orbs%nkptsp))*orbs%nspinor*orbs%norb),intent(in) :: psi
!!!!!  integer,dimension(nspin,0:orbs%nkpts),intent(in):: ndim_ovrlp
!!!!!  real(wp),dimension(ndim_ovrlp(nspin,orbs%nkpts)),intent(out):: ovrlp
!!!!!  integer,dimension(nspin),intent(in):: norbTot
!!!!!
!!!!!  ! Local variables
!!!!!  integer:: ispsi,ikptp,ikpt,ispin,nspinor,ncomp,norbs,ierr,nvctrp,norb
!!!!!
!!!!!
!!!!!
!!!!!  ! Set the whole overlap matrix to zero. This is necessary since each process treats only a part
!!!!!  ! of the matrix.
!!!!!  call to_zero(ndim_ovrlp(nspin,orbs%nkpts),ovrlp)
!!!!!
!!!!!
!!!!!  ispsi=1
!!!!!  ! First make a loop over the k points handled by this process.
!!!!!  do ikptp=1,orbs%nkptsp
!!!!!     ! ikpt is the index of the k point.
!!!!!     ikpt=orbs%iskpts+ikptp
!!!!!
!!!!!     ! Now make also a loop over spin up/down.
!!!!!     do ispin=1,nspin
!!!!!        do jproc=lproc,uproc
!!!!!           do iorb=1,orbs%norb_par(jproc)
!!!!!
!!!!!              ! This subroutine gives essentially back nvctrp, i.e. the length of the vectors for which the overlap
!!!!!              ! matrix shall be calculated. In addition it sets the value of nspinor to orbs%nspinor.
!!!!!              call orbitals_and_components(iproc,ikptp,ispin,orbs,comms,&
!!!!!                   nvctrp,norb,norbs,ncomp,nspinor, jproc, iorb)
!!!!!              ! The subroutine also overwrite the variable norb with the total number of orbitals.
!!!!!              ! However we want to keep the value of norbIn (since we treat only a part of the orbitals).
!!!!!              norb=norbIn
!!!!!              ! Put the starting index to the right place. The current block of vector starts at the block1-th vector.
!!!!!              ispsi=ispsi+nvctrp*(block1-1)*nspinor
!!!!!              if(ispin==ispinIn) then
!!!!!                 if (nvctrp == 0) cycle
!!!!!
!!!!!                 ! Now calclulate one part of the overlap matrix. The starting index of this part is given by ndim_ovrlp(ispin,ikpt-1)+1.
!!!!!                 if(nspinor==1) then
!!!!!                    call syrk('L','T',norb,nvctrp,1.0_wp,psi(ispsi),max(1,nvctrp),&
!!!!!                         0.0_wp,ovrlp(ndim_ovrlp(ispin,ikpt-1)+1),norb)
!!!!!                 else
!!!!!                    call herk('L','C',norb,ncomp*nvctrp,1.0_wp,psi(ispsi),&
!!!!!                         max(1,ncomp*nvctrp),0.0_wp,ovrlp(ndim_ovrlp(ispin,ikpt-1)+1),norb)
!!!!!                 end if
!!!!!              end if
!!!!!              ! Move the starting indices to the end of the actual k point. This is necessary since nvctrp is
!!!!!              ! different for the next k point and we cannot jump directly to the starting indices of our block for 
!!!!!              ! the next k point.
!!!!!              ispsi=ispsi+nvctrp*(norbTot(ispin)-block1+1)*nspinor
!!!!!
!!!!!           end do
!!!!!        end do
!!!!!     end do
!!!!!  end do
!!!!!
!!!!!  !call MPI_BARRIER(bigdft_mpi%mpi_comm,ierr)
!!!!!  !print *,'here',iproc
!!!!!
!!!!!  if (nproc > 1) then
!!!!!     !call timing(iproc,'GramS_comput  ','OF')
!!!!!     !call timing(iproc,'GramS_commun  ','ON')
!!!!!     call timing(iproc, trim(category)//'_comput', 'OF')
!!!!!     call timing(iproc, trim(category)//'_commun', 'ON')
!!!!!     !call mpiallred(ovrlp(1),ndim_ovrlp(nspin,orbs%nkpts),MPI_SUM,bigdft_mpi%mpi_comm,ierr)
!!!!!     call mpiallred(ovrlp(1),ndim_ovrlp(nspin,orbs%nkpts),MPI_SUM,newComm,ierr)
!!!!!     !call MPI_ALLREDUCE (ovrlp(1,2),ovrlp(1,1),ndim_ovrlp(nspin,orbs%nkpts),mpidtypw,MPI_SUM,bigdft_mpi%mpi_comm,ierr)
!!!!!     call timing(iproc, trim(category)//'_commun', 'OF')
!!!!!     call timing(iproc, trim(category)//'_comput', 'ON')
!!!!!     !call timing(iproc,'GramS_commun  ','OF')
!!!!!     !call timing(iproc,'GramS_comput  ','ON')
!!!!!  end if
!!!!!
!!!!!  ! Now each processors knows all the overlap matrices for each k-point
!!!!!  ! even if it does not handle it.
!!!!!  ! This is somehow redundant but it is one way of reducing the number of communications
!!!!!  ! without defining group of processors.
!!!!!
!!!!!end subroutine getOverlapLIN
!!!!!
!!!!!
!!!!!subroutine orbitals_and_componentsLIN(iproc,ikptp,ispin,orbs,comms,nvctrp,norb,norbs,ncomp,nspinor, jproc, iorb)
!!!!!  use module_base
!!!!!  use module_types
!!!!!  implicit none
!!!!!  integer, intent(in) :: iproc,ikptp,ispin, jproc, iorb
!!!!!  type(orbitals_data), intent(in) :: orbs
!!!!!  type(comms_cubic), intent(in) :: comms
!!!!!  integer, intent(out) :: nvctrp,norb,norbs,ncomp,nspinor
!!!!!
!!!!!  !nvctrp=comms%nvctr_par(iproc,ikptp)
!!!!!  nvctrp=comms%nvctr_parLIN(iorb,jproc,iproc,ikptp)
!!!!!  norb=orbs%norbu
!!!!!  nspinor=orbs%nspinor
!!!!!  if (ispin==2) norb=orbs%norbd
!!!!!
!!!!!  call complex_components(nspinor,norb,norbs,ncomp)

!!!!!END SUBROUTINE orbitals_and_componentsLIN
