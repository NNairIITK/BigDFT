subroutine orthonormalizeLocalized(iproc, nproc, lin, input, lphi, ovrlp)
use module_base
use module_types
implicit none

! Calling arguments
integer,intent(in):: iproc, nproc
type(linearParameters),intent(inout):: lin
type(input_variables),intent(in):: input
real(8),dimension(lin%lorbs%npsidim),intent(inout):: lphi
real(8),dimension(lin%lzd%orbs%norb,lin%lzd%orbs%norb),intent(out):: ovrlp

! Local variables
integer:: it, istat, iall, iorb, jorb, ierr
real(8),dimension(:),allocatable:: lphiovrlp
character(len=*),parameter:: subname='orthonormalize'
logical:: converged
real(8):: maxError, t1, t2, timeCommun, timeComput, timeCalcOvrlp, t3, t4, timeExpand, timeLoewdin, timeTransform, timeExtract

  allocate(lphiovrlp(lin%op%ndim_lphiovrlp), stat=istat)
  call memocc(istat, lphiovrlp, 'lphiovrlp',subname)


  timeComput=0.d0
  timeCommun=0.d0
  timeCalcOvrlp=0.d0
  timeExpand=0.d0
  timeLoewdin=0.d0
  timeTransform=0.d0
  timeExtract=0.d0
  converged=.false.
  do it=1,lin%nItOrtho
      !if(iproc==0) write(*,'(a,i0)') 'at it=',it
      call cpu_time(t1)
      !call extractOrbital(iproc, nproc, lin%orbs, lin%lorbs%npsidim, lin%onWhichAtomAll, lin%lzd, lin%op, lphi, lin%comon)
      call extractOrbital2(iproc, nproc, lin%orbs, lin%lorbs%npsidim, lin%onWhichAtomAll, lin%lzd, lin%op, lphi, lin%comon)
      call cpu_time(t2)
      timeExtract=timeExtract+t2-t1
      timeComput=timeComput+t2-t1
      call cpu_time(t1)
      call postCommsOverlap(iproc, nproc, lin%comon)
      !call gatherOrbitals(iproc, nproc, lin%comon)
      call gatherOrbitals2(iproc, nproc, lin%comon)
      call cpu_time(t2)
      timeCommun=timeCommun+t2-t1
      call cpu_time(t1)
      call calculateOverlapMatrix2(iproc, nproc, lin%lzd%orbs, lin%op, lin%comon, lin%onWhichAtomAll, ovrlp)
      call cpu_time(t2)
      timeCalcOvrlp=timeCalcOvrlp+t2-t1
      !!do iorb=1,lin%orbs%norb
      !!    do jorb=1,lin%orbs%norb
      !!        write(500+iproc,*) iorb, jorb, ovrlp(jorb,iorb)
      !!    end do
      !!end do
      call checkUnity(iproc, lin%lzd%orbs%norb, ovrlp, maxError)
      if(iproc==0) write(*,'(3x,a,es12.4)') 'maximal deviation from unity:', maxError
      if(maxError<lin%convCritOrtho) then
          converged=.true.
          call cpu_time(t2)
          timeComput=timeComput+t2-t1
          exit
      else if(it==lin%nItOrtho) then
          call cpu_time(t2)
          timeComput=timeComput+t2-t1
          exit
      end if
      call cpu_time(t3)
      !call transformOverlapMatrix(iproc, nproc, lin%lzd%orbs, ovrlp)
      call transformOverlapMatrix(iproc, nproc, lin%lzd%orbs%norb, ovrlp)
      call cpu_time(t4)
      timeTransform=timeTransform+t4-t3
      !call transformOverlapMatrix2(iproc, nproc, lin%lzd%orbs, ovrlp)
      call cpu_time(t3)
      !call expandOrbital(iproc, nproc, lin%lzd%orbs, input, lin%onWhichAtomAll, lin%lzd, lin%op, lin%comon, lphiovrlp)
      call expandOrbital2(iproc, nproc, lin%lzd%orbs, input, lin%onWhichAtomAll, lin%lzd, lin%op, lin%comon, lphiovrlp)
      call cpu_time(t4)
      timeExpand=timeExpand+t4-t3
      call cpu_time(t3)
      call globalLoewdin(iproc, nproc, lin%lzd%orbs, lin%lorbs, lin%onWhichAtomAll, lin%lzd, lin%op, ovrlp, lphiovrlp, lphi)
      call cpu_time(t4)
      timeLoewdin=timeLoewdin+t4-t3
      call cpu_time(t2)
      timeComput=timeComput+t2-t1
  end do

  if(converged) then
      if(iproc==0) write(*,'(3x,a,i0,a)') 'done in ', it, ' iterations.'
  else 
      if(iproc==0) write(*,'(3x,a,i0,a)') 'WARNING: orthonormalization not converged within ', lin%nItOrtho, ' iterations.'
  end if

  call mpiallred(timeComput, 1, mpi_sum, mpi_comm_world, ierr)
  call mpiallred(timeCommun, 1, mpi_sum, mpi_comm_world, ierr)
  call mpiallred(timeCalcOvrlp, 1, mpi_sum, mpi_comm_world, ierr)
  call mpiallred(timeExpand, 1, mpi_sum, mpi_comm_world, ierr)
  call mpiallred(timeLoewdin, 1, mpi_sum, mpi_comm_world, ierr)
  call mpiallred(timeTransform, 1, mpi_sum, mpi_comm_world, ierr)
  call mpiallred(timeExtract, 1, mpi_sum, mpi_comm_world, ierr)
  timeComput=timeComput/dble(nproc)
  timeCommun=timeCommun/dble(nproc)
  timeCalcOvrlp=timeCalcOvrlp/dble(nproc)
  timeExpand=timeExpand/dble(nproc)
  timeLoewdin=timeLoewdin/dble(nproc)
  timeTransform=timeTransform/dble(nproc)
  timeExtract=timeExtract/dble(nproc)
  if(iproc==0) write(*,'(3x,a,es9.3,a,f4.1,a)') 'time for computation:', timeComput, '=',&
                100.d0*timeComput/(timeComput+timeCommun), '%'
  if(iproc==0) write(*,'(3x,a,es9.3,a,f4.1,a)') 'time for communication:', timeCommun, '=',&
               100.d0*timeCommun/(timeComput+timeCommun), '%'
  if(iproc==0) write(*,'(3x,a,es9.3,a,f4.1,a)') 'time for calculating overlap:', timeCalcOvrlp, '=',&
               100.d0*timeCalcOvrlp/(timeComput+timeCommun), '%'
  if(iproc==0) write(*,'(3x,a,es9.3,a,f4.1,a)') 'time for expansion:', timeExpand, '=',&
               100.d0*timeExpand/(timeComput+timeCommun), '%'
  if(iproc==0) write(*,'(3x,a,es9.3,a,f4.1,a)') 'time for Loewdin:', timeLoewdin, '=',&
               100.d0*timeLoewdin/(timeComput+timeCommun), '%'
  if(iproc==0) write(*,'(3x,a,es9.3,a,f4.1,a)') 'time for transform:', timeTransform, '=',&
               100.d0*timeTransform/(timeComput+timeCommun), '%'
  if(iproc==0) write(*,'(3x,a,es9.3,a,f4.1,a)') 'time for extract:', timeExtract, '=',&
               100.d0*timeExtract/(timeComput+timeCommun), '%'

  iall=-product(shape(lphiovrlp))*kind(lphiovrlp)
  deallocate(lphiovrlp, stat=istat)
  call memocc(istat, iall, 'lphiovrlp', subname)


end subroutine orthonormalizeLocalized


subroutine orthoconstraintLocalized(iproc, nproc, lin, input, lphi, lhphi, trH)
use module_base
use module_types
implicit none

! Calling arguments
integer,intent(in):: iproc, nproc
type(linearParameters),intent(inout):: lin
type(input_variables),intent(in):: input
real(8),dimension(lin%lorbs%npsidim),intent(in):: lphi
real(8),dimension(lin%lorbs%npsidim),intent(inout):: lhphi
real(8),intent(out):: trH

! Local variables
integer:: it, istat, iall, iorb
real(8),dimension(:),allocatable:: lphiovrlp
real(8),dimension(:,:),allocatable:: lagmat
character(len=*),parameter:: subname='orthoconstraintLocalized'



  allocate(lagmat(lin%orbs%norb,lin%orbs%norb), stat=istat)
  call memocc(istat, lagmat, 'lagmat',subname)
  allocate(lphiovrlp(lin%op%ndim_lphiovrlp), stat=istat)
  call memocc(istat, lphiovrlp, 'lphiovrlp',subname)


  ! Put lphi in the sendbuffer, i.e. lphi will be sent to other processes' receive buffer.
  !call extractOrbital(iproc, nproc, lin%orbs, lin%lorbs%npsidim, lin%onWhichAtomAll, lin%lzd, lin%op, lphi, lin%comon)
  call extractOrbital2(iproc, nproc, lin%orbs, lin%lorbs%npsidim, lin%onWhichAtomAll, lin%lzd, lin%op, lphi, lin%comon)
  call postCommsOverlap(iproc, nproc, lin%comon)
  !call gatherOrbitals(iproc, nproc, lin%comon)
  call gatherOrbitals2(iproc, nproc, lin%comon)
  ! Put lhphi to the sendbuffer, so we can the calculate <lphi|lhphi>
  !call extractOrbital(iproc, nproc, lin%orbs, lin%lorbs%npsidim, lin%onWhichAtomAll, lin%lzd, lin%op, lhphi, lin%comon)
  call extractOrbital2(iproc, nproc, lin%orbs, lin%lorbs%npsidim, lin%onWhichAtomAll, lin%lzd, lin%op, lhphi, lin%comon)
  call calculateOverlapMatrix2(iproc, nproc, lin%lzd%orbs, lin%op, lin%comon, lin%onWhichAtomAll, lagmat)
  trH=0.d0
  do iorb=1,lin%orbs%norb
      trH=trH+lagmat(iorb,iorb)
  end do
  ! Expand the receive buffer, i.e. lphi
  !call expandOrbital(iproc, nproc, lin%lzd%orbs, input, lin%onWhichAtomAll, lin%lzd, lin%op, lin%comon, lphiovrlp)
  call expandOrbital2(iproc, nproc, lin%lzd%orbs, input, lin%onWhichAtomAll, lin%lzd, lin%op, lin%comon, lphiovrlp)
  call applyOrthoconstraint(iproc, nproc, lin%lzd%orbs, lin%lorbs, lin%onWhichAtomAll, lin%lzd, lin%op, lagmat, lphiovrlp, lhphi)

  iall=-product(shape(lagmat))*kind(lagmat)
  deallocate(lagmat, stat=istat)
  call memocc(istat, iall, 'lagmat', subname)
  iall=-product(shape(lphiovrlp))*kind(lphiovrlp)
  deallocate(lphiovrlp, stat=istat)
  call memocc(istat, iall, 'lphiovrlp', subname)


end subroutine orthoconstraintLocalized


subroutine orthoconstraintNonorthogonal(iproc, nproc, lin, input, ovrlp, lphi, lhphi, trH)
use module_base
use module_types
implicit none

! Calling arguments
integer,intent(in):: iproc, nproc
type(linearParameters),intent(inout):: lin
type(input_variables),intent(in):: input
real(8),dimension(lin%lzd%orbs%norb,lin%lzd%orbs%norb),intent(in):: ovrlp
real(8),dimension(lin%lorbs%npsidim),intent(in):: lphi
real(8),dimension(lin%lorbs%npsidim),intent(inout):: lhphi
real(8),intent(out):: trH

! Local variables
integer:: it, istat, iall, iorb
real(8),dimension(:),allocatable:: lphiovrlp, lhphiovrlp
real(8),dimension(:,:),allocatable:: lagmat
character(len=*),parameter:: subname='orthoconstraintLocalized'



  allocate(lagmat(lin%orbs%norb,lin%orbs%norb), stat=istat)
  call memocc(istat, lagmat, 'lagmat',subname)
  allocate(lphiovrlp(lin%op%ndim_lphiovrlp), stat=istat)
  call memocc(istat, lphiovrlp, 'lphiovrlp',subname)
  allocate(lhphiovrlp(lin%op%ndim_lphiovrlp), stat=istat)
  call memocc(istat, lhphiovrlp, 'lhphiovrlp',subname)


  ! Put lphi in the sendbuffer, i.e. lphi will be sent to other processes' receive buffer.
  call extractOrbital2(iproc, nproc, lin%orbs, lin%lorbs%npsidim, lin%onWhichAtomAll, lin%lzd, lin%op, lphi, lin%comon)
  call postCommsOverlap(iproc, nproc, lin%comon)
  call gatherOrbitals2(iproc, nproc, lin%comon)
  ! Put lhphi to the sendbuffer, so we can the calculate <lphi|lhphi>
  call extractOrbital2(iproc, nproc, lin%orbs, lin%lorbs%npsidim, lin%onWhichAtomAll, lin%lzd, lin%op, lhphi, lin%comon)
  call calculateOverlapMatrix2(iproc, nproc, lin%lzd%orbs, lin%op, lin%comon, lin%onWhichAtomAll, lagmat)
  trH=0.d0
  do iorb=1,lin%orbs%norb
      trH=trH+lagmat(iorb,iorb)
  end do
  ! Expand the receive buffer, i.e. lphi
  !call expandOrbital(iproc, nproc, lin%lzd%orbs, input, lin%onWhichAtomAll, lin%lzd, lin%op, lin%comon, lphiovrlp)
  call expandOrbital2(iproc, nproc, lin%lzd%orbs, input, lin%onWhichAtomAll, lin%lzd, lin%op, lin%comon, lphiovrlp)

  ! Now we also have to send lhphi
  call extractOrbital2(iproc, nproc, lin%orbs, lin%lorbs%npsidim, lin%onWhichAtomAll, lin%lzd, lin%op, lhphi, lin%comon)
  call postCommsOverlap(iproc, nproc, lin%comon)
  call gatherOrbitals2(iproc, nproc, lin%comon)
  ! Expand the receive buffer, i.e. lhphi
  call expandOrbital2(iproc, nproc, lin%lzd%orbs, input, lin%onWhichAtomAll, lin%lzd, lin%op, lin%comon, lhphiovrlp)
  !call applyOrthoconstraintNonorthogonal(iproc, nproc, lin%lzd%orbs, lin%lorbs, lin%onWhichAtomAll, lin%lzd, lin%op, lagmat, ovrlp, lphiovrlp, lhphiovrlp, lhphi)
  call applyOrthoconstraintNonorthogonal2(iproc, nproc, lin%lzd%orbs, lin%lorbs, lin%onWhichAtomAll, lin%lzd, lin%op,&
       lagmat, ovrlp, lphiovrlp, lhphiovrlp, lhphi)

  iall=-product(shape(lagmat))*kind(lagmat)
  deallocate(lagmat, stat=istat)
  call memocc(istat, iall, 'lagmat', subname)
  iall=-product(shape(lphiovrlp))*kind(lphiovrlp)
  deallocate(lphiovrlp, stat=istat)
  call memocc(istat, iall, 'lphiovrlp', subname)
  iall=-product(shape(lhphiovrlp))*kind(lhphiovrlp)
  deallocate(lhphiovrlp, stat=istat)
  call memocc(istat, iall, 'lhphiovrlp', subname)


end subroutine orthoconstraintNonorthogonal




subroutine orthoconstraintLocalized2(iproc,nproc,orbs,comms,wfd,psi,hpsi,scprsum,lagMatDiag)
  use module_base
  use module_types
  implicit none
  integer, intent(in) :: iproc,nproc
  type(orbitals_data), intent(in) :: orbs
  type(communications_arrays), intent(in) :: comms
  type(wavefunctions_descriptors), intent(in) :: wfd
  real(wp), dimension(sum(comms%nvctr_par(iproc,1:orbs%nkptsp))*orbs%nspinor*orbs%norb), intent(in) :: psi
  real(wp), dimension(sum(comms%nvctr_par(iproc,1:orbs%nkptsp))*orbs%nspinor*orbs%norb), intent(out) :: hpsi
  real(dp), intent(out) :: scprsum
  real(dp),dimension(orbs%norb),intent(out):: lagMatDiag
  !local variables
  character(len=*), parameter :: subname='orthoconstraintNotSymmetric'
  integer :: i_stat,i_all,ierr,iorb,ise,jorb
  integer :: ispin,nspin,ikpt,norb,norbs,ncomp,nvctrp,ispsi,ikptp,nspinor
  real(dp) :: occ,tt
  integer, dimension(:,:), allocatable :: ndimovrlp
  real(wp), dimension(:), allocatable :: alag


integer:: istart, jstart


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

  allocate(ndimovrlp(nspin,0:orbs%nkpts+ndebug),stat=i_stat)
  call memocc(i_stat,ndimovrlp,'ndimovrlp',subname)

  call dimension_ovrlp(nspin,orbs,ndimovrlp)

  allocate(alag(ndimovrlp(nspin,orbs%nkpts)+ndebug),stat=i_stat)
  call memocc(i_stat,alag,'alag',subname)

  !put to zero all the k-points which are not needed
  call razero(ndimovrlp(nspin,orbs%nkpts),alag)

  !do it for each of the k-points and separate also between up and down orbitals in the non-collinear case
  ispsi=1
  do ikptp=1,orbs%nkptsp
     ikpt=orbs%iskpts+ikptp!orbs%ikptsp(ikptp)

     do ispin=1,nspin

        call orbitals_and_components(iproc,ikptp,ispin,orbs,comms,&
             nvctrp,norb,norbs,ncomp,nspinor)
        if (nvctrp == 0) cycle

        if(nspinor==1) then
           call gemm('T','N',norb,norb,nvctrp,1.0_wp,psi(ispsi),&
                max(1,nvctrp),hpsi(ispsi),max(1,nvctrp),0.0_wp,&
                alag(ndimovrlp(ispin,ikpt-1)+1),norb)
        else
           !this part should be recheck in the case of nspinor == 2
           call c_gemm('C','N',norb,norb,ncomp*nvctrp,(1.0_wp,0.0_wp),psi(ispsi),&
                max(1,ncomp*nvctrp), &
                hpsi(ispsi),max(1,ncomp*nvctrp),(0.0_wp,0.0_wp),&
                alag(ndimovrlp(ispin,ikpt-1)+1),norb)
        end if
        ispsi=ispsi+nvctrp*norb*nspinor
     end do
  end do

  if (nproc > 1) then
     call timing(iproc,'LagrM_comput  ','OF')
     call timing(iproc,'LagrM_commun  ','ON')
     call mpiallred(alag(1),ndimovrlp(nspin,orbs%nkpts),MPI_SUM,MPI_COMM_WORLD,ierr)
     call timing(iproc,'LagrM_commun  ','OF')
     call timing(iproc,'LagrM_comput  ','ON')
  end if

  ! Copy the diagonal of the matrix
  do iorb=1,orbs%norb
      lagMatDiag(iorb)=alag((iorb-1)*orbs%norb+iorb)
  end do
! Lagrange multiplier matrix
!!if(iproc==0) write(*,*) 'Lagrange multiplier matrix'
!!do iorb=1,norb
!!    !if(iproc==0) write(*,'(80f8.4)') (alag(norb*jorb+iorb), jorb=0,norb-1)
!!    do jorb=1,norb
!!        write(1100+iproc,*) iorb, jorb, alag(norb*(iorb-1)+jorb)
!!    end do
!!end do


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
        if (ispin==1) ise=0
        call orbitals_and_components(iproc,ikptp,ispin,orbs,comms,&
             nvctrp,norb,norbs,ncomp,nspinor)
        if (nvctrp == 0) cycle

!!$        !correct the orthogonality constraint if there are some orbitals which have zero occupation number
!!$        do iorb=1,norb
!!$           do jorb=iorb+1,norb
!!$              if (orbs%occup((ikpt-1)*orbs%norb+iorb+ise) /= 0.0_gp .and. &
!!$                   orbs%occup((ikpt-1)*orbs%norb+jorb+ise) == 0.0_gp) then
!!$                 alag(ndimovrlp(ispin,ikpt-1)+iorb+(jorb-1)*norbs) = 0.0_wp
!!$                 alag(ndimovrlp(ispin,ikpt-1)+jorb+(iorb-1)*norbs) = 0.0_wp
!!$                 !if (iproc ==0) print *,'i,j',iorb,jorb,alag(ndimovrlp(ispin,ikpt-1)+iorb+(jorb-1)*norbs)
!!$              end if
!!$           end do
!!$        end do

        !calculate the scprsum if the k-point is associated to this processor
        !the scprsum always coincide with the trace of the hamiltonian
        if (orbs%ikptproc(ikpt) == iproc) then
           occ=real(orbs%kwgts(ikpt),dp)
           if(nspinor == 1) then
              do iorb=1,norb
                 scprsum=scprsum+occ*real(alag(ndimovrlp(ispin,ikpt-1)+iorb+(iorb-1)*norbs),dp)
              enddo
           else if (nspinor == 4 .or. nspinor == 2) then
              !not sure about the imaginary part of the diagonal (should be zero if H is hermitian)
              do iorb=1,norb
                 scprsum=scprsum+&
                      occ*real(alag(ndimovrlp(ispin,ikpt-1)+2*iorb-1+(iorb-1)*norbs),dp)
                 scprsum=scprsum+&
                      occ*real(alag(ndimovrlp(ispin,ikpt-1)+2*iorb+(iorb-1)*norbs),dp)
              enddo
           end if
        end if
        ise=norb

        if(nspinor==1 .and. nvctrp /= 0) then
           !call gemm('N','N',nvctrp,norb,norb,-1.0_wp,psi(ispsi),max(1,nvctrp),&
           !     alag(ndimovrlp(ispin,ikpt-1)+1),norb,1.0_wp,&
           !     hpsi(ispsi),max(1,nvctrp))
           call gemm('N','N',nvctrp,norb,norb,-.5_wp,psi(ispsi),max(1,nvctrp),&
                alag(ndimovrlp(ispin,ikpt-1)+1),norb,1.0_wp,&
                hpsi(ispsi),max(1,nvctrp))
           call gemm('N','T',nvctrp,norb,norb,-.5_wp,psi(ispsi),max(1,nvctrp),&
                alag(ndimovrlp(ispin,ikpt-1)+1),norb,1.0_wp,&
                hpsi(ispsi),max(1,nvctrp))
        else if (nvctrp /= 0) then
           stop 'not implemented for nspinor/=1!'
           call c_gemm('N','N',ncomp*nvctrp,norb,norb,(-1.0_wp,0.0_wp),psi(ispsi),max(1,ncomp*nvctrp),&
                alag(ndimovrlp(ispin,ikpt-1)+1),norb,(1.0_wp,0.0_wp),hpsi(ispsi),max(1,ncomp*nvctrp))
        end if
        ispsi=ispsi+nvctrp*norb*nspinor
     end do
  end do

  if (nproc > 1) then
     tt=scprsum
     call mpiallred(scprsum,1,MPI_SUM,MPI_COMM_WORLD,ierr)
     !call MPI_ALLREDUCE(tt,scprsum,1,mpidtypd,MPI_SUM,MPI_COMM_WORLD,ierr)
  end if

  i_all=-product(shape(alag))*kind(alag)
  deallocate(alag,stat=i_stat)
  call memocc(i_stat,i_all,'alag',subname)

  i_all=-product(shape(ndimovrlp))*kind(ndimovrlp)
  deallocate(ndimovrlp,stat=i_stat)
  call memocc(i_stat,i_all,'ndimovrlp',subname)

  call timing(iproc,'LagrM_comput  ','OF')

END SUBROUTINE orthoconstraintLocalized2


subroutine getOverlapMatrix(iproc, nproc, lin, input, lphi, ovrlp)
use module_base
use module_types
implicit none

! Calling arguments
integer,intent(in):: iproc, nproc
type(linearParameters),intent(inout):: lin
type(input_variables),intent(in):: input
real(8),dimension(lin%lorbs%npsidim),intent(inout):: lphi
real(8),dimension(lin%lzd%orbs%norb,lin%lzd%orbs%norb),intent(out):: ovrlp

! Local variables
integer:: it, istat, iall, iorb, jorb, ierr
real(8),dimension(:),allocatable:: lphiovrlp
character(len=*),parameter:: subname='orthonormalize'
logical:: converged

  allocate(lphiovrlp(lin%op%ndim_lphiovrlp), stat=istat)
  call memocc(istat, lphiovrlp, 'lphiovrlp',subname)

  call extractOrbital2(iproc, nproc, lin%orbs, lin%lorbs%npsidim, lin%onWhichAtomAll, lin%lzd, lin%op, lphi, lin%comon)
  call postCommsOverlap(iproc, nproc, lin%comon)
  call gatherOrbitals2(iproc, nproc, lin%comon)
  call calculateOverlapMatrix2(iproc, nproc, lin%lzd%orbs, lin%op, lin%comon, lin%onWhichAtomAll, ovrlp)

  iall=-product(shape(lphiovrlp))*kind(lphiovrlp)
  deallocate(lphiovrlp, stat=istat)
  call memocc(istat, iall, 'lphiovrlp', subname)


end subroutine getOverlapMatrix


subroutine getOverlapMatrix2(iproc, nproc, lin, input, lphi, ovrlp)
use module_base
use module_types
implicit none

! Calling arguments
integer,intent(in):: iproc, nproc
type(linearParameters),intent(inout):: lin
type(input_variables),intent(in):: input
real(8),dimension(lin%lb%lzd%orbs%npsidim),intent(inout):: lphi
real(8),dimension(lin%lb%lzd%orbs%norb,lin%lb%lzd%orbs%norb),intent(out):: ovrlp

! Local variables
integer:: it, istat, iall, iorb, jorb, ierr
real(8),dimension(:),allocatable:: lphiovrlp
character(len=*),parameter:: subname='orthonormalize'
logical:: converged

  !!allocate(lphiovrlp(lin%op_lb%ndim_lphiovrlp), stat=istat)
  !!call memocc(istat, lphiovrlp, 'lphiovrlp',subname)

  call extractOrbital2(iproc,nproc,lin%lb%orbs,lin%lb%lorbs%npsidim,lin%lb%onWhichAtomAll,lin%lb%lzd,lin%op_lb,lphi,lin%comon_lb)
  call postCommsOverlap(iproc, nproc, lin%comon_lb)
  call gatherOrbitals2(iproc, nproc, lin%comon_lb)
  call calculateOverlapMatrix2(iproc, nproc, lin%lb%lzd%orbs, lin%op_lb, lin%comon_lb, lin%lb%onWhichAtomAll, ovrlp)

  !!iall=-product(shape(lphiovrlp))*kind(lphiovrlp)
  !!deallocate(lphiovrlp, stat=istat)
  !!call memocc(istat, iall, 'lphiovrlp', subname)


end subroutine getOverlapMatrix2



subroutine initCommsOrtho(iproc, nproc, lzd, onWhichAtomAll, input, op, comon, tag)
use module_base
use module_types
use module_interfaces, exceptThisOne => initCommsOrtho
implicit none

! Calling arguments
integer,intent(in):: iproc, nproc
type(linear_zone_descriptors),intent(in):: lzd
integer,dimension(lzd%orbs%norb),intent(in):: onWhichAtomAll
type(input_variables),intent(in):: input
type(overlapParameters),intent(out):: op
type(p2pCommsOrthonormality),intent(out):: comon
integer,intent(inout):: tag

! Local variables
integer:: iorb, jorb, iiorb, jproc, ioverlaporb, ioverlapMPI, ilr, jlr
integer:: ilrold, is1, ie1, is2, ie2, is3, ie3, js1, je1, js2, je2, js3
integer::  je3, istat
logical:: ovrlpx, ovrlpy, ovrlpz
character(len=*),parameter:: subname='initCommsOrtho'

allocate(comon%noverlaps(0:nproc-1), stat=istat)
call memocc(istat, comon%noverlaps, 'comon%noverlaps',subname)
allocate(op%noverlaps(lzd%orbs%norb), stat=istat)
call memocc(istat, op%noverlaps, 'op%noverlaps',subname)

! Count how many overlaping regions each orbital / process has.
call countOverlaps(iproc, nproc, lzd%orbs, lzd, onWhichAtomAll, op, comon)

allocate(op%overlaps(maxval(op%noverlaps),lzd%orbs%norb), stat=istat)
call memocc(istat, op%overlaps, 'op%overlaps', subname)
allocate(comon%overlaps(maxval(comon%noverlaps),0:nproc-1), stat=istat)
call memocc(istat, comon%overlaps, 'comon%overlaps', subname)
allocate(op%indexInRecvBuf(lzd%orbs%norbp,lzd%orbs%norb), stat=istat)
call memocc(istat, op%indexInRecvBuf, 'op%indexInRecvBuf', subname)
allocate(op%indexInSendBuf(lzd%orbs%norbp,lzd%orbs%norb), stat=istat)
call memocc(istat, op%indexInSendBuf, 'op%indexInSendBuf', subname)

! Determine the overlapping orbitals.
call determineOverlaps(iproc, nproc, lzd%orbs, lzd, onWhichAtomAll, op, comon)

allocate(op%olr(maxval(op%noverlaps(lzd%orbs%isorb+1:lzd%orbs%isorb+lzd%orbs%norbp)),lzd%orbs%norbp), stat=istat)

! Set the orbital descriptors for the overlap regions.
call determineOverlapDescriptors(iproc, nproc, lzd%orbs, lzd, lzd%Glr, onWhichAtomAll, op)


allocate(comon%comarr(8,maxval(comon%noverlaps),0:nproc-1), stat=istat)
call memocc(istat, comon%comarr, 'comon%comarr', subname)
allocate(comon%communComplete(maxval(comon%noverlaps),0:nproc-1), stat=istat)
call memocc(istat, comon%communComplete, 'comun%communComplete', subname)
call setCommsOrtho(iproc, nproc, lzd%orbs, onWhichAtomAll, lzd, op, comon, tag)

allocate(comon%recvBuf(comon%nrecvBuf), stat=istat)
call memocc(istat, comon%recvBuf, 'comon%recvBuf', subname)
allocate(comon%sendBuf(comon%nsendBuf), stat=istat)
call memocc(istat, comon%sendBuf, 'comon%sendBuf', subname)
comon%recvBuf=555.55d0
!write(*,'(a,2i11)') 'iproc, comon%nsendBuf', iproc, comon%nsendBuf

!call postCommsOverlap(iproc, nproc, comon)

! Initialize the index arrays for the transformations from overlap region
! to ordinary localization region.
allocate(op%indexExpand(comon%nrecvBuf), stat=istat)
call memocc(istat, op%indexExpand, 'op%indexExpand',subname)
call indicesForExpansion(iproc, nproc, lzd%orbs, input, onWhichAtomAll, lzd, op, comon)

! Initialize the index arrays for the transformations from the ordinary localization region
! to the overlap region.
allocate(op%indexExtract(comon%nsendBuf), stat=istat)
call memocc(istat, op%indexExtract, 'op%indexExtract',subname)
call indicesForExtraction(iproc, nproc, lzd%orbs, lzd%orbs%npsidim, onWhichAtomAll, lzd, op, comon)


!!do istat=1,comon%nrecvBuf
!!    write(400+iproc,*) istat, op%indexExpand(istat)
!!end do
!!do istat=1,comon%nsendBuf
!!    write(450+iproc,*) istat, op%indexExtract(istat)
!!end do


!!call mpi_barrier(mpi_comm_world, istat)
!!stop






end subroutine initCommsOrtho







!> Returns the starting and ending indices (on the coarse grid) of a given localization region.
subroutine getIndices(lr, is1, ie1, is2, ie2, is3, ie3)
use module_base
use module_types
implicit none

! Calling arguments
type(locreg_descriptors),intent(in):: lr
integer,intent(out):: is1, ie1, is2, ie2, is3, ie3

  is1=lr%ns1+1
  ie1=lr%ns1+lr%d%n1
  is2=lr%ns2+1
  ie2=lr%ns2+lr%d%n2
  is3=lr%ns3+1
  ie3=lr%ns3+lr%d%n3

end subroutine getIndices




! Count for each orbital and each process the number of overlapping orbitals.
subroutine countOverlaps(iproc, nproc, orbs, lzd, onWhichAtom, op, comon)
use module_base
use module_types
implicit none

! Calling arguments
integer,intent(in):: iproc, nproc
type(orbitals_data),intent(in):: orbs
type(linear_zone_descriptors),intent(in):: lzd
integer,dimension(orbs%norb),intent(in):: onWhichAtom
type(overlapParameters),intent(out):: op
type(p2pCommsOrthonormality),intent(out):: comon

! Local variables
integer:: jproc, iorb, jorb, ioverlapMPI, ioverlaporb, ilr, jlr, ilrold, is1, ie1, is2, ie2, is3, ie3
integer:: js1, je1, js2, je2, js3, je3, iiorb
logical:: ovrlpx, ovrlpy, ovrlpz

iiorb=0
do jproc=0,nproc-1
    ioverlapMPI=0 ! counts the overlaps for the given MPI process.
    ilrold=-1
    do iorb=1,orbs%norb_par(jproc)
        ioverlaporb=0 ! counts the overlaps for the given orbital.
        iiorb=iiorb+1 ! counts the total orbitals
        ilr=onWhichAtom(iiorb)
        call getIndices(lzd%llr(ilr), is1, ie1, is2, ie2, is3, ie3)
        do jorb=1,orbs%norb
            jlr=onWhichAtom(jorb)
            call getIndices(lzd%llr(jlr), js1, je1, js2, je2, js3, je3)
            ovrlpx = ( is1<=je1 .and. ie1>=js1 )
            ovrlpy = ( is2<=je2 .and. ie2>=js2 )
            ovrlpz = ( is3<=je3 .and. ie3>=js3 )
            !if(iproc==0) write(*,'(a,6i5,5x,6i5,5x,3l)') 'is1, ie1, is2, ie2, is3, ie3   js1, je1, js2, je2, js3, je3  ovrlpx, ovrlpy, ovrlpz', &
            !  is1, ie1, is2, ie2, is3, ie3, js1, je1, js2, je2, js3, je3, ovrlpx, ovrlpy, ovrlpz
            if(ovrlpx .and. ovrlpy .and. ovrlpz) then
                ioverlaporb=ioverlaporb+1
                if(ilr/=ilrold) then
                    ! if ilr==ilrold, we are in the same localization region, so the MPI prosess
                    ! would get the same orbitals again. Therefore the counter is not increased
                    ! in that case.
                    ioverlapMPI=ioverlapMPI+1
                end if
            end if
        end do 
        op%noverlaps(iiorb)=ioverlaporb
        !!if(iproc==0) write(*,'(a,2i8)') 'iiorb, op%noverlaps(iiorb)', iiorb, op%noverlaps(iiorb)
        ilrold=ilr
    end do
    comon%noverlaps(jproc)=ioverlapMpi
    !if(iproc==0) write(*,'(a,2i8)') 'jproc, comon%noverlaps(jproc)', jproc, comon%noverlaps(jproc)
end do

end subroutine countOverlaps




subroutine determineOverlaps(iproc, nproc, orbs, lzd, onWhichAtom, op, comon)
use module_base
use module_types
implicit none

! Calling arguments
integer,intent(in):: iproc, nproc
type(orbitals_data),intent(in):: orbs
type(linear_zone_descriptors),intent(in):: lzd
integer,dimension(orbs%norb),intent(in):: onWhichAtom
type(overlapParameters),intent(out):: op
type(p2pCommsOrthonormality),intent(out):: comon

! Local variables
integer:: jproc, iorb, jorb, ioverlapMPI, ioverlaporb, ilr, jlr, ilrold, is1, ie1, is2, ie2, is3, ie3
integer:: js1, je1, js2, je2, js3, je3, iiorb
logical:: ovrlpx, ovrlpy, ovrlpz

  ! Initialize to some value which will never be used.
  op%overlaps=-1
  comon%overlaps=-1

  iiorb=0
  do jproc=0,nproc-1
      ioverlapMPI=0 ! counts the overlaps for the given MPI process.
      ilrold=-1
      do iorb=1,orbs%norb_par(jproc)
          ioverlaporb=0 ! counts the overlaps for the given orbital.
          iiorb=iiorb+1 ! counts the total orbitals
          ilr=onWhichAtom(iiorb)
          call getIndices(lzd%llr(ilr), is1, ie1, is2, ie2, is3, ie3)
          do jorb=1,orbs%norb
              jlr=onWhichAtom(jorb)
              call getIndices(lzd%llr(jlr), js1, je1, js2, je2, js3, je3)
              ovrlpx = ( is1<=je1 .and. ie1>=js1 )
              ovrlpy = ( is2<=je2 .and. ie2>=js2 )
              ovrlpz = ( is3<=je3 .and. ie3>=js3 )
              !if(iproc==0) write(*,'(a,6i5,5x,6i5,5x,3l)') 'is1, ie1, is2, ie2, is3, ie3   js1, je1, js2, je2, js3, je3  ovrlpx, ovrlpy, ovrlpz', &
              !  is1, ie1, is2, ie2, is3, ie3, js1, je1, js2, je2, js3, je3, ovrlpx, ovrlpy, ovrlpz
              if(ovrlpx .and. ovrlpy .and. ovrlpz) then
                  ioverlaporb=ioverlaporb+1
                  op%overlaps(ioverlaporb,iiorb)=jorb
                  if(ilr/=ilrold) then
                      ! if ilr==ilrold, we are in th same localization region, so the MPI prosess
                      ! would get the same orbitals again. Therefore the counter is not increased
                      ! in that case.
                      ioverlapMPI=ioverlapMPI+1
                      comon%overlaps(ioverlapMPI,jproc)=jorb
                  end if
              end if
          end do 
          !if(iproc==0) write(*,'(a,i3,5x,100i5)') 'iiorb, op%overlaps', iiorb, op%overlaps(:,iiorb) 
          ilrold=ilr
      end do
      !if(iproc==0) write(*,'(a,i3,5x,100i5)') 'jproc, comon%overlaps', jproc, comon%overlaps(:,jproc) 
  end do

end subroutine determineOverlaps



subroutine determineOverlapDescriptors(iproc, nproc, orbs, lzd, Glr, onWhichAtom, op)
use module_base
use module_types
implicit none

! Calling arguments
integer,intent(in):: iproc, nproc
type(orbitals_data),intent(in):: orbs
type(linear_zone_descriptors),intent(in):: lzd
type(locreg_descriptors),intent(in):: Glr
integer,dimension(orbs%norb),intent(in):: onWhichAtom
type(overlapParameters),intent(inout):: op

! Local variables
integer:: iorb, jorb, jjorb, ilr, jlr, iiorb


do iorb=1,orbs%norbp
    iiorb=lzd%orbs%isorb_par(iproc)+iorb
    ilr=onWhichAtom(iiorb)
!if(iproc==0) write(*,'(a,2i10)') 'iorb, op%noverlaps(iorb)', iorb, op%noverlaps(iorb)
    do jorb=1,op%noverlaps(iiorb)
        jjorb=op%overlaps(jorb,iiorb)
        jlr=onWhichAtom(jjorb)
!write(*,*) 'calling get_overlap_region_periodic'
        !call get_overlap_region_periodic(ilr, jlr, Glr, 1, lzd%llr, lzd%nlr, op%olr(jorb,iorb))
        call get_overlap_region_periodic2(ilr, jlr, Glr, 1, lzd%llr, lzd%nlr, op%olr(jorb,iorb))
        !write(*,'(a,13i8)') 'iproc, iorb, jorb, iiorb, jjorb, ilr, jlr, nvctr_c, nvctr_f, ncount, n1, n2, n3', iproc, iorb, jorb, iiorb, jjorb, ilr, jlr, &
        !    op%olr(jorb,iorb)%wfd%nvctr_c, op%olr(jorb,iorb)%wfd%nvctr_f, op%olr(jorb,iorb)%wfd%nvctr_c+7*op%olr(jorb,iorb)%wfd%nvctr_f, op%olr(jorb,iorb)%d%n1, op%olr(jorb,iorb)%d%n2, op%olr(jorb,iorb)%d%n3
    end do
end do

end subroutine determineOverlapDescriptors




subroutine setCommsOrtho(iproc, nproc, orbs, onWhichAtom, lzd, op, comon, tag)
use module_base
use module_types
implicit none

! Calling arguments
integer,intent(in):: iproc, nproc
type(orbitals_data),intent(in):: orbs
integer,dimension(orbs%norb),intent(in):: onWhichAtom
type(linear_zone_descriptors),intent(in):: lzd
type(overlapParameters),intent(inout):: op
type(p2pCommsOrthonormality),intent(out):: comon
integer,intent(inout):: tag

! Local variables
integer:: jproc, iorb, jorb, iiorb, jjorb, mpisource, mpidest, istsource, istdest, ncount, istat, iall, ijorb
integer:: ilr, ilrold, jprocold, ildim, ierr
integer,dimension(:),allocatable:: istsourceArr, istdestArr
character(len=*),parameter:: subname='setCommsOrtho'
logical,dimension(:),allocatable:: receivedOrbital

allocate(istsourceArr(0:nproc-1), stat=istat)
call memocc(istat, istsourceArr, 'istsourceArr',subname)
allocate(istdestArr(0:nproc-1), stat=istat)
call memocc(istat, istdestArr, 'istdestArr',subname)
allocate(receivedOrbital(orbs%norb), stat=istat)
call memocc(istat, receivedOrbital, 'receivedOrbital', subname)

istsourceArr=1
istdestArr=1

comon%nsendBuf=0
comon%nrecvBuf=0


op%indexInRecvBuf=0

iiorb=0
jprocold=-1
do jproc=0,nproc-1
    receivedOrbital=.false.
    ijorb=0
    ilrold=-1
    do iorb=1,orbs%norb_par(jproc)
       iiorb=iiorb+1 
       ilr=onWhichAtom(iiorb)
       ! Check whether process jproc has already received orbital jjorb.
       if(ilr==ilrold) cycle
       ildim=lzd%llr(ilr)%wfd%nvctr_c+7*lzd%llr(ilr)%wfd%nvctr_f
       do jorb=1,op%noverlaps(iiorb)
           jjorb=op%overlaps(jorb,iiorb)
           !write(*,'(a,7i8)') 'iproc, iiorb, jjorb, ilr, ilrold, jproc, jprocold', iproc, iiorb, jjorb, ilr, ilrold, jproc, jprocold
           ijorb=ijorb+1
           mpisource=orbs%onWhichMPI(jjorb)
           mpidest=jproc
           istsource=istsourceArr(mpisource)
           istdest=istdestArr(mpidest)
           if(iproc==jproc) then
               ncount=op%olr(jorb,iorb)%wfd%nvctr_c+7*op%olr(jorb,iorb)%wfd%nvctr_f
               !write(*,'(a,4i9)') 'iproc, iorb, jorb, ncount', iproc, iorb, jorb, ncount
           end if
           call mpi_bcast(ncount, 1, mpi_integer, jproc, mpi_comm_world, ierr)
           tag=tag+1
           receivedOrbital(jjorb)=.true.
           call setCommsParameters(mpisource, mpidest, istsource, istdest, ncount, tag, comon%comarr(1,ijorb,jproc))
           !if(iproc==0) write(*,'(6(a,i0))') 'process ',mpisource,' sends ',ncount,' elements from position ',istsource,' to position ',istdest,' on process ',mpidest,', tag=',tag
           if(iproc==mpisource) then
               !write(*,'(4(a,i0))') 'adding ',ncount,' elements for orbital ',iiorb,' to nsendBuf, iproc=',iproc,', jproc=',jproc
               comon%nsendBuf=comon%nsendBuf+ncount
           end if
           if(iproc==mpidest) then
               !write(*,'(3(a,i0))') 'process ',iproc,' will get orbital ',jjorb,' at position ',istdest
               op%indexInRecvBuf(iorb,jjorb)=istdest
               comon%nrecvBuf=comon%nrecvBuf+ncount
               op%ndim_lphiovrlp=op%ndim_lphiovrlp+ildim
           end if
           istsourceArr(mpisource) = istsourceArr(mpisource) + ncount
           istdestArr(mpidest) = istdestArr(mpidest) + ncount
       end do
       ilrold=ilr
    end do
    jprocold=jproc
end do
   
iall = -product(shape(istsourceArr))*kind(istsourceArr)
deallocate(istsourceArr, stat=istat)
call memocc(istat, iall, 'istsourceArr',subname)
iall = -product(shape(istdestArr))*kind(istdestArr)
deallocate(istdestArr, stat=istat)
call memocc(istat, iall, 'istdestArr',subname)
iall = -product(shape(receivedOrbital))*kind(receivedOrbital)
deallocate(receivedOrbital, stat=istat)
call memocc(istat, iall, 'receivedOrbital',subname)


end subroutine setCommsOrtho



subroutine setCommsParameters(mpisource, mpidest, istsource, istdest, ncount, tag, comarr)
use module_base
use module_types
implicit none

! Calling arguments
integer,intent(in):: mpisource, mpidest, istsource, istdest, ncount, tag
integer,dimension(8),intent(out):: comarr


! From which MPI process shall the orbital be sent.
comarr(1)=mpisource

! Starting index on the sending process.
comarr(2)=istsource

! Amount of datat to be sent
comarr(3)=ncount

! To which MPI process shall the orbital be sent.
comarr(4)=mpidest

! Starting index on the receiving process.
comarr(5)=istdest

! Tag for the communication
comarr(6)=tag

! comarr(7): this entry is used as request for the mpi_isend.

! comarr(8): this entry is used as request for the mpi_irecv.



end subroutine setCommsParameters




subroutine postCommsOverlap(iproc, nproc, comon)
use module_base
use module_types
implicit none

! Calling arguments
integer,intent(in):: iproc, nproc
type(p2pCommsOrthonormality),intent(inout):: comon

! Local variables
integer:: jproc, iorb, mpisource, istsource, ncount, mpidest, istdest, tag, nsends, nreceives, ierr


nsends=0
nreceives=0
comon%communComplete=.false.
do jproc=0,nproc-1
    do iorb=1,comon%noverlaps(jproc)
        mpisource=comon%comarr(1,iorb,jproc)
        istsource=comon%comarr(2,iorb,jproc)
        ncount=comon%comarr(3,iorb,jproc)
        mpidest=comon%comarr(4,iorb,jproc)
        istdest=comon%comarr(5,iorb,jproc)
        tag=comon%comarr(6,iorb,jproc)
        if(mpisource/=mpidest) then
            ! The orbitals are on different processes, so we need a point to point communication.
            if(iproc==mpisource) then
                !write(*,'(6(a,i0))') 'process ', mpisource, ' sends ', ncount, ' elements from position ', istsource, ' to position ', istdest, ' on process ', mpidest, ', tag=',tag
                call mpi_isend(comon%sendBuf(istsource), ncount, mpi_double_precision, mpidest, tag,&
                     mpi_comm_world, comon%comarr(7,iorb,jproc), ierr)
                !call mpi_isend(sendBuf(istsource), ncount, mpi_double_precision, mpidest, tag, mpi_comm_world, lin%comsr%comarr(8,iorb,jproc), ierr)
                comon%comarr(8,iorb,jproc)=mpi_request_null !is this correct?
                nsends=nsends+1
            else if(iproc==mpidest) then
                !write(*,'(6(a,i0))') 'process ', mpidest, ' receives ', ncount, ' elements at position ', istdest, ' from position ', istsource, ' on process ', mpisource, ', tag=',tag
                call mpi_irecv(comon%recvBuf(istdest), ncount, mpi_double_precision, mpisource, tag,&
                     mpi_comm_world, comon%comarr(8,iorb,jproc), ierr)
                comon%comarr(7,iorb,jproc)=mpi_request_null !is this correct?
                nreceives=nreceives+1
            else
                comon%comarr(7,iorb,jproc)=mpi_request_null
                comon%comarr(8,iorb,jproc)=mpi_request_null
            end if
        else
            ! The orbitals are on the same process, so simply copy them.
            if(iproc==mpisource) then
                call dcopy(ncount, comon%sendBuf(istsource), 1, comon%recvBuf(istdest), 1)
                !write(*,'(6(a,i0),a,2es15.7)') 'process ', iproc, ' copies ', ncount, ' elements from position ', istsource, ' to position ', istdest, ' on process ', iproc, ', tag=',tag,'. comon%sendBuf(istsource), comon%recvBuf(istdest)', comon%sendBuf(istsource), comon%recvBuf(istdest)
                comon%comarr(7,iorb,jproc)=mpi_request_null
                comon%comarr(8,iorb,jproc)=mpi_request_null
                nsends=nsends+1
                nreceives=nreceives+1
                comon%communComplete(iorb,iproc)=.true.
            else
                comon%comarr(7,iorb,jproc)=mpi_request_null
                comon%comarr(8,iorb,jproc)=mpi_request_null
            end if

        end if
    end do
end do


end subroutine postCommsOverlap



subroutine extractOrbital(iproc, nproc, orbs, sizePhi, onWhichAtom, lzd, op, phi, comon)
use module_base
use module_types
implicit none

! Calling arguments
integer,intent(in):: iproc, nproc, sizePhi
type(orbitals_data),intent(in):: orbs
integer,dimension(orbs%norb),intent(in):: onWhichAtom
type(linear_zone_descriptors),intent(in):: lzd
type(overlapParameters),intent(inout):: op
real(8),dimension(sizePhi),intent(in):: phi
type(p2pCommsOrthonormality),intent(out):: comon

! Local variables
integer:: iorb, jorb, korb, ind, indovrlp, ilr, klr, ilrold, jjorb, jjlr, jjproc, iiproc, iiprocold, gdim, ldim, kkorb, lorb
integer:: i

indovrlp=1
op%indexInSendBuf=0

ilrold=-1
do iorb=1,orbs%norb
    ilr=onWhichAtom(iorb)
    iiproc=orbs%onWhichMPI(iorb)
    if(ilr==ilrold .and. iiproc==iiprocold) cycle ! otherwise we would extract the same again
    do jorb=1,op%noverlaps(iorb)
        jjorb=op%overlaps(jorb,iorb)
        jjlr=onWhichAtom(jjorb)
        jjproc=orbs%onWhichMPI(jjorb)
        if(iproc==jjproc) then
            ! Get the correct descriptors
            korb=jjorb-orbs%isorb
            !write(*,'(a,5i8)') 'iorb, jorb, jjorb, jjproc, korb', iorb, jorb, jjorb, jjproc, korb
            do i=1,op%noverlaps(jjorb)
                !write(*,'(a,5i8)') 'iproc, iorb, korb, i, op%overlaps(i,korb)', iproc, iorb, korb, i, op%overlaps(i,korb)
                if(op%overlaps(i,jjorb)==iorb) then
                    lorb=i
                    exit
                end if
            end do
            !write(*,'(a,5i9)') 'iproc, iorb, jorb, korb, lorb', iproc, iorb, jorb, korb, lorb
            gdim=lzd%llr(jjlr)%wfd%nvctr_c+7*lzd%llr(jjlr)%wfd%nvctr_f
            ldim=op%olr(lorb,korb)%wfd%nvctr_c+7*op%olr(lorb,korb)%wfd%nvctr_f
            ind=1
            do kkorb=orbs%isorb+1,jjorb-1
                klr=onWhichAtom(kkorb)
                ind = ind + lzd%llr(klr)%wfd%nvctr_c + 7*lzd%llr(klr)%wfd%nvctr_f
            end do
            !write(*,'(5(a,i0))') 'process ',iproc,' adds ',op%olr(lorb,korb)%wfd%nvctr_c+7*op%olr(lorb,korb)%wfd%nvctr_f,' elements at position ',indovrlp,' from orbital ',jjorb,' for orbital ', iorb
            call psi_to_locreg2(iproc, nproc, ldim, gdim, op%olr(lorb,korb), lzd%llr(jjlr), phi(ind), comon%sendBuf(indovrlp))
            op%indexInSendBuf(jjorb-orbs%isorb,iorb)=indovrlp
            indovrlp=indovrlp+op%olr(lorb,korb)%wfd%nvctr_c+7*op%olr(lorb,korb)%wfd%nvctr_f
        end if
    end do
    ilrold=ilr
    iiprocold=iiproc
end do

if(indovrlp/=comon%nsendBuf+1) then
    write(*,'(x,a,i0,a,3x,i0,2x,i0)') 'ERROR on process ', iproc, ': indovrlp/=comon%nsendBuf+1', indovrlp, comon%nsendBuf+1
    stop
end if



end subroutine extractOrbital




subroutine extractOrbital2(iproc, nproc, orbs, sizePhi, onWhichAtom, lzd, op, phi, comon)
use module_base
use module_types
implicit none

! Calling arguments
integer,intent(in):: iproc, nproc, sizePhi
type(orbitals_data),intent(in):: orbs
integer,dimension(orbs%norb),intent(in):: onWhichAtom
type(linear_zone_descriptors),intent(in):: lzd
type(overlapParameters),intent(inout):: op
real(8),dimension(sizePhi),intent(in):: phi
type(p2pCommsOrthonormality),intent(out):: comon

! Local variables
integer:: iorb, jorb, korb, ind, indovrlp, ilr, klr, ilrold, jjorb, jjlr, jjproc, iiproc, iiprocold, gdim, ldim, kkorb, lorb
integer:: i, indSource

indovrlp=1
op%indexInSendBuf=0

ilrold=-1
iiprocold=-1
do iorb=1,orbs%norb
    ilr=onWhichAtom(iorb)
    iiproc=orbs%onWhichMPI(iorb)
    if(ilr==ilrold .and. iiproc==iiprocold) cycle ! otherwise we would extract the same again
    do jorb=1,op%noverlaps(iorb)
        jjorb=op%overlaps(jorb,iorb)
        jjlr=onWhichAtom(jjorb)
        jjproc=orbs%onWhichMPI(jjorb)
        if(iproc==jjproc) then
            ! Get the correct descriptors
            korb=jjorb-orbs%isorb
            !write(*,'(a,5i8)') 'iorb, jorb, jjorb, jjproc, korb', iorb, jorb, jjorb, jjproc, korb
            do i=1,op%noverlaps(jjorb)
                !write(*,'(a,5i8)') 'iproc, iorb, korb, i, op%overlaps(i,korb)', iproc, iorb, korb, i, op%overlaps(i,korb)
                if(op%overlaps(i,jjorb)==iorb) then
                    lorb=i
                    exit
                end if
            end do
            !write(*,'(a,5i9)') 'iproc, iorb, jorb, korb, lorb', iproc, iorb, jorb, korb, lorb
            gdim=lzd%llr(jjlr)%wfd%nvctr_c+7*lzd%llr(jjlr)%wfd%nvctr_f
            ldim=op%olr(lorb,korb)%wfd%nvctr_c+7*op%olr(lorb,korb)%wfd%nvctr_f
            ind=1
            do kkorb=orbs%isorb+1,jjorb-1
                klr=onWhichAtom(kkorb)
                ind = ind + lzd%llr(klr)%wfd%nvctr_c + 7*lzd%llr(klr)%wfd%nvctr_f
            end do
            !write(*,'(5(a,i0))') 'process ',iproc,' adds ',op%olr(lorb,korb)%wfd%nvctr_c+7*op%olr(lorb,korb)%wfd%nvctr_f,' elements at position ',indovrlp,' from orbital ',jjorb,' for orbital ', iorb
            !call psi_to_locreg2(iproc, nproc, ldim, gdim, op%olr(lorb,korb), lzd%llr(jjlr), phi(ind), comon%sendBuf(indovrlp))
            do i=0,ldim-1
                indSource=ind+op%indexExtract(indovrlp+i)-1
                comon%sendBuf(indovrlp+i)=phi(indSource)
            end do
            op%indexInSendBuf(jjorb-orbs%isorb,iorb)=indovrlp
            indovrlp=indovrlp+op%olr(lorb,korb)%wfd%nvctr_c+7*op%olr(lorb,korb)%wfd%nvctr_f
        end if
    end do
    ilrold=ilr
    iiprocold=iiproc
end do

if(indovrlp/=comon%nsendBuf+1) then
    write(*,'(x,a,i0,a,3x,i0,2x,i0)') 'ERROR on process ', iproc, ': indovrlp/=comon%nsendBuf+1', indovrlp, comon%nsendBuf+1
    stop
end if

end subroutine extractOrbital2





subroutine gatherOrbitals(iproc, nproc, comon)
use module_base
use module_types
implicit none

! Calling arguments
integer,intent(in):: iproc, nproc
type(p2pCommsOrthonormality),intent(inout):: comon

! Local variables
integer:: kproc, mpisource, mpidest, nfast, nslow, nsameproc, ierr
integer,dimension(mpi_status_size):: stat
logical:: sendComplete, receiveComplete


! Check whether the communications have completed.
nfast=0
nsameproc=0
testLoop: do
    do kproc=1,comon%noverlaps(iproc)
        if(comon%communComplete(kproc,iproc)) cycle
        call mpi_test(comon%comarr(7,kproc,iproc), sendComplete, stat, ierr)
        call mpi_test(comon%comarr(8,kproc,iproc), receiveComplete, stat, ierr)
        if(sendComplete .and. receiveComplete) comon%communComplete(kproc,iproc)=.true.
        if(comon%communComplete(kproc,iproc)) then
            !write(*,'(2(a,i0))') 'fast communication; process ', iproc, ' has received orbital ', kproc
            mpisource=comon%comarr(1,kproc,iproc)
            mpidest=comon%comarr(4,kproc,iproc)
            if(mpisource/=mpidest) then
                nfast=nfast+1
            else
                nsameproc=nsameproc+1
            end if
        end if
    end do
    ! If we made it until here, either all all the communication is
    ! complete or we better wait for each single orbital.
    exit testLoop
end do testLoop


! Wait for the communications that have not completed yet
nslow=0
do kproc=1,comon%noverlaps(iproc)
    if(comon%communComplete(kproc,iproc)) cycle
    !write(*,'(2(a,i0))') 'process ', iproc, ' is waiting for orbital ', korb
    nslow=nslow+1
    call mpi_wait(comon%comarr(7,kproc,iproc), stat, ierr)   !COMMENTED BY PB
    call mpi_wait(comon%comarr(8,kproc,iproc), stat, ierr)   !COMMENTED BY PB
    comon%communComplete(kproc,iproc)=.true.
end do

!call mpiallred(nreceives, 1, mpi_sum, mpi_comm_world, ierr)
call mpiallred(nfast, 1, mpi_sum, mpi_comm_world, ierr)
call mpiallred(nslow, 1, mpi_sum, mpi_comm_world, ierr)
call mpiallred(nsameproc, 1, mpi_sum, mpi_comm_world, ierr)
!if(iproc==0) write(*,'(x,2(a,i0),a)') 'statistics: - ', nfast+nslow, ' point to point communications, of which ', &
!                       nfast, ' could be overlapped with computation.'
!if(iproc==0) write(*,'(x,a,i0,a)') '            - ', nsameproc, ' copies on the same processor.'



end subroutine gatherOrbitals





subroutine gatherOrbitals2(iproc, nproc, comon)
use module_base
use module_types
implicit none

! Calling arguments
integer,intent(in):: iproc, nproc
type(p2pCommsOrthonormality),intent(inout):: comon

! Local variables
integer:: jorb, mpisource, mpidest, nfast, nslow, nsameproc, ierr, jproc
integer,dimension(mpi_status_size):: stat
logical:: sendComplete, receiveComplete


! Check whether the communications have completed.
nfast=0
nsameproc=0
testLoop: do
    do jproc=0,nproc-1
        do jorb=1,comon%noverlaps(jproc)
            if(comon%communComplete(jorb,jproc)) cycle
            call mpi_test(comon%comarr(7,jorb,jproc), sendComplete, stat, ierr)
            call mpi_test(comon%comarr(8,jorb,jproc), receiveComplete, stat, ierr)
            if(sendComplete .and. receiveComplete) comon%communComplete(jorb,jproc)=.true.
            if(comon%communComplete(jorb,jproc)) then
                !write(*,'(2(a,i0))') 'fast communication; process ', iproc, ' has received orbital ', jorb
                mpisource=comon%comarr(1,jorb,jproc)
                mpidest=comon%comarr(4,jorb,jproc)
                if(mpisource/=mpidest) then
                    nfast=nfast+1
                else
                    nsameproc=nsameproc+1
                end if
            end if
        end do
    end do
    ! If we made it until here, either all all the communication is
    ! complete or we better wait for each single orbital.
    exit testLoop
end do testLoop


! Wait for the communications that have not completed yet
nslow=0
do jproc=0,nproc-1
    do jorb=1,comon%noverlaps(jproc)
        if(comon%communComplete(jorb,jproc)) cycle
        !write(*,'(2(a,i0))') 'process ', iproc, ' is waiting for orbital ', korb
        nslow=nslow+1
        call mpi_wait(comon%comarr(7,jorb,jproc), stat, ierr)   !COMMENTED BY PB
        call mpi_wait(comon%comarr(8,jorb,jproc), stat, ierr)   !COMMENTED BY PB
        comon%communComplete(jorb,jproc)=.true.
    end do
end do

!call mpiallred(nreceives, 1, mpi_sum, mpi_comm_world, ierr)
call mpiallred(nfast, 1, mpi_sum, mpi_comm_world, ierr)
call mpiallred(nslow, 1, mpi_sum, mpi_comm_world, ierr)
call mpiallred(nsameproc, 1, mpi_sum, mpi_comm_world, ierr)
!if(iproc==0) write(*,'(x,2(a,i0),a)') 'statistics: - ', nfast+nslow, ' point to point communications, of which ', &
!                       nfast, ' could be overlapped with computation.'
!if(iproc==0) write(*,'(x,a,i0,a)') '            - ', nsameproc, ' copies on the same processor.'





end subroutine gatherOrbitals2





subroutine calculateOverlapMatrix(iproc, nproc, orbs, op, comon, onWhichAtom, lovrlp)
use module_base
use module_types
implicit none

! Calling arguments
integer,intent(in):: iproc, nproc
type(orbitals_data),intent(in):: orbs
type(overlapParameters),intent(in):: op
type(p2pCommsOrthonormality),intent(inout):: comon
integer,dimension(orbs%norb),intent(in):: onWhichAtom
real(8),dimension(maxval(op%noverlaps),orbs%norbp),intent(out):: lovrlp

! Local variables
integer:: iorb, jorb, iiorb, jjorb, jjproc, ii, ist, jst, jjlr, istat, ncount
real(8):: ddot

lovrlp=0.d0


do iorb=1,orbs%norbp
    iiorb=orbs%isorb+iorb
    do jorb=1,op%noverlaps(iiorb)
        jjorb=op%overlaps(jorb,iiorb)
        ! Starting index of orbital iorb, already transformed to overlap region with jjorb.
        ! We have to find the first orbital on the same MPI and in the same locreg as jjorb.
        jjlr=onWhichAtom(jjorb)
        jjproc=orbs%onWhichMPI(jjorb)
        ii=orbs%isorb_par(jjproc)+1
        do
            if(onWhichAtom(ii)==jjlr) exit
            ii=ii+1
        end do
        ist=op%indexInSendBuf(iorb,ii)
        ! Starting index of orbital jjorb
        jst=op%indexInRecvBuf(iorb,jjorb)
        !write(*,'(5(a,i0))') 'process ',iproc,' calculates overlap of ',iiorb,' and ',jjorb,'. ist=',ist,' jst=',jst 
        !ncount=op%olr(jorb,iiorb)%wfd%nvctr_c+7*op%olr(jorb,iiorb)%wfd%nvctr_f
        ncount=op%olr(jorb,iorb)%wfd%nvctr_c+7*op%olr(jorb,iorb)%wfd%nvctr_f
        !write(*,'(a,4i8)') 'iproc, iiorb, jjorb, ncount', iproc, iiorb, jjorb, ncount
        !ovrlp(iiorb,jjorb)=ddot(ncount, comon%sendBuf(ist), 1, comon%recvBuf(jst), 1)
        lovrlp(jorb,iorb)=ddot(ncount, comon%sendBuf(ist), 1, comon%recvBuf(jst), 1)
    end do
end do

!!do iorb=1,orbs%norbp
!!    iiorb=orbs%isorb+iorb
!!    do jorb=1,op%noverlaps(iiorb)
!!        write(500+iproc,*) iorb, jorb, lovrlp(jorb,iorb)
!!    end do
!!end do



end subroutine calculateOverlapMatrix




subroutine calculateOverlapMatrix2(iproc, nproc, orbs, op, comon, onWhichAtom, ovrlp)
use module_base
use module_types
implicit none

! Calling arguments
integer,intent(in):: iproc, nproc
type(orbitals_data),intent(in):: orbs
type(overlapParameters),intent(in):: op
type(p2pCommsOrthonormality),intent(inout):: comon
integer,dimension(orbs%norb),intent(in):: onWhichAtom
real(8),dimension(orbs%norb,orbs%norb),intent(out):: ovrlp

! Local variables
integer:: iorb, jorb, iiorb, jjorb, jjproc, ii, ist, jst, jjlr, istat, ncount, ierr, jj, korb, kkorb, iilr, iiproc
real(8):: ddot
real(8),dimension(orbs%norb**2):: buffer


ovrlp=0.d0

do iorb=1,orbs%norbp
    iiorb=orbs%isorb+iorb
    do jorb=1,op%noverlaps(iiorb)
        jjorb=op%overlaps(jorb,iiorb)
        ! Starting index of orbital iorb, already transformed to overlap region with jjorb.
        ! We have to find the first orbital on the same MPI and in the same locreg as jjorb.
        jjlr=onWhichAtom(jjorb)
        jjproc=orbs%onWhichMPI(jjorb)
        jj=orbs%isorb_par(jjproc)+1
        do
            if(onWhichAtom(jj)==jjlr) exit
            jj=jj+1
        end do
        ist=op%indexInSendBuf(iorb,jj)
        ! Starting index of orbital jjorb.
        ! We have to find the first orbital on the same MPI and in the same locreg as iiorb.
        iilr=onWhichAtom(iiorb)
        iiproc=orbs%onWhichMPI(iiorb)
        do korb=1,orbs%norbp
            kkorb=orbs%isorb_par(iiproc)+korb
            if(onWhichAtom(kkorb)==iilr) then
                ii=korb
                exit
            end if
        end do
        jst=op%indexInRecvBuf(ii,jjorb)
        !write(*,'(5(a,i0))') 'process ',iproc,' calculates overlap of ',iiorb,' and ',jjorb,'. ist=',ist,' jst=',jst 
        !ncount=op%olr(jorb,iiorb)%wfd%nvctr_c+7*op%olr(jorb,iiorb)%wfd%nvctr_f
        ncount=op%olr(jorb,iorb)%wfd%nvctr_c+7*op%olr(jorb,iorb)%wfd%nvctr_f
        !write(*,'(a,4i8)') 'iproc, iiorb, jjorb, ncount', iproc, iiorb, jjorb, ncount
        !ovrlp(iiorb,jjorb)=ddot(ncount, comon%sendBuf(ist), 1, comon%recvBuf(jst), 1)
        ovrlp(iiorb,jjorb)=ddot(ncount, comon%sendBuf(ist), 1, comon%recvBuf(jst), 1)
    end do
end do


call mpiallred(ovrlp(1,1), orbs%norb**2, mpi_sum, mpi_comm_world, ierr)

!!call dcopy(orbs%norb**2, ovrlp(1,1), 1, buffer(1), 1)
!!call mpi_allreduce(buffer(1), ovrlp(1,1), orbs%norb**2, mpi_double_precision, &
!!     mpi_sum, mpi_comm_world, ierr)



end subroutine calculateOverlapMatrix2






subroutine transformOverlapMatrix(iproc, nproc, norb, ovrlp)
use module_base
use module_types
implicit none

! Calling arguments
integer,intent(in):: iproc, nproc, norb
real(8),dimension(norb,norb),intent(inout):: ovrlp

! Local variables
integer:: lwork, istat, iall, iorb, jorb, info
real(8),dimension(:),allocatable:: work, eval
real(8),dimension(:,:,:),allocatable:: tempArr
character(len=*),parameter:: subname='transformOverlapMatrix'


allocate(eval(norb), stat=istat)
call memocc(istat, eval, 'eval', subname)
allocate(tempArr(norb,norb,2), stat=istat)
call memocc(istat, tempArr, 'tempArr', subname)

!!write(450+iproc,*) '---------------------------------------'
!!do iorb=1,norb
!!  do jorb=1,norb
!!      write(450+iproc,*) iorb, jorb, ovrlp(jorb,iorb)
!!      if(abs(ovrlp(iorb,jorb)-ovrlp(jorb,iorb))>1.d-10 ) then
!!          if(iproc==0) write(*,'(a,3i7,3es20.12)') 'ERROR: not symmetric, iproc, iorb, jorb, ovrlp(iorb,jorb), &
!!               &ovrlp(jorb,iorb), abs(ovrlp(iorb,jorb)-ovrlp(jorb,iorb))', iproc, iorb, jorb, ovrlp(iorb,jorb), &
!!               ovrlp(jorb,iorb), abs(ovrlp(iorb,jorb)-ovrlp(jorb,iorb))
!!          !stop
!!      end if
!!  end do
!!end do
!!do iorb=1,norb
!!    do jorb=1,norb
!!        if(iproc==0) write(1401,'(2i6,es26.17)') iorb, jorb, ovrlp(iorb,jorb)
!!    end do
!!end do


lwork=1000*norb
allocate(work(lwork), stat=istat)
call memocc(istat, work, 'work', subname)
call dsyev('v', 'l', norb, ovrlp(1,1), norb, eval, work, lwork, info)
if(info/=0) then
    write(*,'(a,i0)') 'ERROR in dsyev, info=', info
    stop
end if
iall=-product(shape(work))*kind(work)
deallocate(work, stat=istat)
call memocc(istat, iall, 'work', subname)
!!do iorb=1,norb
!!    do jorb=1,norb
!!        if(iproc==0) write(1402,'(2i6,es26.17)') iorb, jorb, ovrlp(iorb,jorb)
!!    end do
!!    if(iproc==0) write(1450,*) iorb, eval(iorb)
!!end do

! Calculate S^{-1/2}. 
! First calulate ovrlp*diag(1/sqrt(evall)) (ovrlp is the diagonalized overlap
! matrix and diag(1/sqrt(evall)) the diagonal matrix consisting of the inverse square roots of the eigenvalues...
do iorb=1,norb
    do jorb=1,norb
        tempArr(jorb,iorb,1)=ovrlp(jorb,iorb)*1.d0/sqrt(eval(iorb))
    end do
end do
!!do iorb=1,norb
!!    do jorb=1,norb
!!        if(iproc==0) write(1403,'(2i6,es26.17)') iorb, jorb, temparr(iorb,jorb,1)
!!    end do
!!end do

! ...and now apply the diagonalized overlap matrix to the matrix constructed above.
! This will give S^{-1/2}.
call dgemm('n', 't', norb, norb, norb, 1.d0, ovrlp(1,1), &
     norb, tempArr(1,1,1), norb, 0.d0, &
     tempArr(1,1,2), norb)
call dcopy(norb**2, tempArr(1,1,2), 1, ovrlp(1,1), 1)
!!do iorb=1,norb
!!    do jorb=1,norb
!!        if(iproc==0) write(1405,'(2i6,es26.17)') iorb, jorb, ovrlp(iorb,jorb)
!!    end do
!!end do


iall=-product(shape(eval))*kind(eval)
deallocate(eval, stat=istat)
call memocc(istat, iall, 'eval', subname)
iall=-product(shape(tempArr))*kind(tempArr)
deallocate(tempArr, stat=istat)
call memocc(istat, iall, 'tempArr', subname)


end subroutine transformOverlapMatrix



subroutine expandOrbital(iproc, nproc, orbs, input, onWhichAtom, lzd, op, comon, lphiovrlp)
use module_base
use module_types
implicit none

! Calling arguments
integer,intent(in):: iproc, nproc
type(orbitals_data),intent(in):: orbs
type(input_variables),intent(in):: input
integer,dimension(orbs%norb),intent(in):: onWhichAtom
type(linear_zone_descriptors),intent(in):: lzd
type(overlapParameters),intent(in):: op
type(p2pCommsOrthonormality),intent(in):: comon
real(8),dimension(op%ndim_lphiovrlp),intent(out):: lphiovrlp

! Local variables
integer:: ind, iorb, iiorb, ilr, gdim, ldim, jorb, jjorb, jst, ilrold


lphiovrlp=0.d0

ind=1
ilrold=-1
do iorb=1,orbs%norbp
    iiorb=orbs%isorb+iorb
    ilr=onWhichAtom(iiorb)
    if(ilr==ilrold) cycle
    gdim=lzd%llr(ilr)%wfd%nvctr_c+7*lzd%llr(ilr)%wfd%nvctr_f
    do jorb=1,op%noverlaps(iiorb)
        jjorb=op%overlaps(jorb,iiorb)
        ! Starting index of orbital jjorb
        jst=op%indexInRecvBuf(iorb,jjorb)
        !ldim=op%olr(jorb,iiorb)%wfd%nvctr_c+7*op%olr(jorb,iiorb)%wfd%nvctr_f
        ldim=op%olr(jorb,iorb)%wfd%nvctr_c+7*op%olr(jorb,iorb)%wfd%nvctr_f
        !call Lpsi_to_global2(iproc, nproc, ldim, gdim, orbs%norbp, orbs%nspinor, input%nspin, lzd%llr(ilr), op%olr(jorb,iiorb), comon%recvBuf(jst), lphiovrlp(ind))
        call Lpsi_to_global2(iproc, nproc, ldim, gdim, orbs%norbp, orbs%nspinor, input%nspin, lzd%llr(ilr),&
             op%olr(jorb,iorb), comon%recvBuf(jst), lphiovrlp(ind))
        ind=ind+gdim
    end do
    ilrold=ilr
end do

end subroutine expandOrbital




subroutine expandOrbital2(iproc, nproc, orbs, input, onWhichAtom, lzd, op, comon, lphiovrlp)
use module_base
use module_types
implicit none

! Calling arguments
integer,intent(in):: iproc, nproc
type(orbitals_data),intent(in):: orbs
type(input_variables),intent(in):: input
integer,dimension(orbs%norb),intent(in):: onWhichAtom
type(linear_zone_descriptors),intent(in):: lzd
type(overlapParameters),intent(in):: op
type(p2pCommsOrthonormality),intent(in):: comon
real(8),dimension(op%ndim_lphiovrlp),intent(out):: lphiovrlp

! Local variables
integer:: ind, iorb, iiorb, ilr, gdim, ldim, jorb, jjorb, jst, ilrold, i, indDest


lphiovrlp=0.d0

ind=1
ilrold=-1
do iorb=1,orbs%norbp
    iiorb=orbs%isorb+iorb
    ilr=onWhichAtom(iiorb)
    if(ilr==ilrold) cycle
    gdim=lzd%llr(ilr)%wfd%nvctr_c+7*lzd%llr(ilr)%wfd%nvctr_f
    do jorb=1,op%noverlaps(iiorb)
        jjorb=op%overlaps(jorb,iiorb)
        ! Starting index of orbital jjorb
        jst=op%indexInRecvBuf(iorb,jjorb)
        !ldim=op%olr(jorb,iiorb)%wfd%nvctr_c+7*op%olr(jorb,iiorb)%wfd%nvctr_f
        ldim=op%olr(jorb,iorb)%wfd%nvctr_c+7*op%olr(jorb,iorb)%wfd%nvctr_f
        !call Lpsi_to_global2(iproc, nproc, ldim, gdim, orbs%norbp, orbs%nspinor, input%nspin, lzd%llr(ilr), op%olr(jorb,iiorb), comon%recvBuf(jst), lphiovrlp(ind))
        !call Lpsi_to_global2(iproc, nproc, ldim, gdim, orbs%norbp, orbs%nspinor, input%nspin, lzd%llr(ilr), op%olr(jorb,iorb), comon%recvBuf(jst), lphiovrlp(ind))
        do i=0,ldim-1
            indDest=ind+op%indexExpand(jst+i)-1
            lphiovrlp(indDest)=comon%recvBuf(jst+i)
        end do
        ind=ind+gdim
    end do
    ilrold=ilr
end do

end subroutine expandOrbital2





subroutine checkUnity(iproc, norb, ovrlp, maxError)
  use module_base
  use module_types
  implicit none
  
  ! Calling arguments
  integer,intent(in):: iproc, norb
  real(8),dimension(norb,norb),intent(in):: ovrlp
  real(8),intent(out):: maxError
  
  ! Local variables
  integer:: iorb, jorb
  real(8):: error
  
  maxError=0.d0
  do iorb=1,norb
      do jorb=1,norb
          if(iorb==jorb) then
              error=abs(ovrlp(jorb,iorb)-1.d0)
          else
              error=abs(ovrlp(jorb,iorb))
          end if
          if(error>maxError) then
              maxError=error
          end if
      end do
  end do

end subroutine checkUnity



subroutine localGramschmidt(iproc, nproc, orbs, lorbs, onWhichAtom, lzd, op, comon, lovrlp, lphiovrlp, lphi)
use module_base
use module_types
implicit none

! Calling arguments
integer,intent(in):: iproc, nproc
type(orbitals_data),intent(in):: orbs, lorbs
integer,dimension(orbs%norb),intent(in):: onWhichAtom
type(linear_zone_descriptors),intent(in):: lzd
type(overlapParameters),intent(in):: op
type(p2pCommsOrthonormality),intent(in):: comon
real(8),dimension(maxval(op%noverlaps),orbs%norbp),intent(in):: lovrlp
real(8),dimension(op%ndim_lphiovrlp),intent(in):: lphiovrlp
real(8),dimension(lorbs%npsidim),intent(inout):: lphi

! Local variables
integer:: ist, jst, ilr, ilrold, iorb, iiorb, jorb, ncount, jjorb
real(8):: dnrm2, ddot, tt



ist=1
jst=1
ilrold=-1
do iorb=1,orbs%norbp
    iiorb=orbs%isorb+iorb
    ilr=onWhichAtom(iiorb)
    if(ilr==ilrold) then
        ! Set back the index of lphiovrlp, since we again need the same orbitals.
        jst=jst-op%noverlaps(iiorb)*ncount
    end if
    ncount=lzd%llr(ilr)%wfd%nvctr_c+7*lzd%llr(ilr)%wfd%nvctr_f
    do jorb=1,op%noverlaps(iiorb)
        jjorb=op%overlaps(jorb,iiorb)
        if(iiorb==jjorb) then
            jst=jst+ncount
            cycle
        end if
        !write(*,'(a,3i7,3x,es15.7)') 'before: iproc, iorb, jorb, ddot', iproc, iorb, jorb, ddot(ncount, lphiovrlp(jst), 1, lphi(ist), 1) 
        tt=ddot(ncount, lphiovrlp(jst), 1, lphiovrlp(jst), 1)
        call daxpy(ncount, -lovrlp(jorb,iorb)/tt, lphiovrlp(jst), 1, lphi(ist), 1)
        !write(*,'(a,3i7,3x,es15.7)') 'after: iproc, iorb, jorb, ddot', iproc, iorb, jorb, ddot(ncount, lphiovrlp(jst), 1, lphi(ist), 1) 
        jst=jst+ncount
    end do
    ! Normalize
    tt=dnrm2(ncount, lphi(ist), 1)
    call dscal(ncount, 1/tt, lphi(ist), 1)
    ist=ist+ncount
    ilrold=ilr
end do


end subroutine localGramschmidt




subroutine globalLoewdin(iproc, nproc, orbs, lorbs, onWhichAtom, lzd, op, ovrlp, lphiovrlp, lphi)
use module_base
use module_types
implicit none

! Calling arguments
integer,intent(in):: iproc, nproc
type(orbitals_data),intent(in):: orbs, lorbs
integer,dimension(orbs%norb),intent(in):: onWhichAtom
type(linear_zone_descriptors),intent(in):: lzd
type(overlapParameters),intent(in):: op
real(8),dimension(orbs%norb,orbs%norb),intent(in):: ovrlp
real(8),dimension(op%ndim_lphiovrlp),intent(in):: lphiovrlp
real(8),dimension(lorbs%npsidim),intent(out):: lphi

! Local variables
integer:: iorb, jorb, iiorb, ilr, ist, jst, ilrold, jjorb, ncount
real(8):: tt, dnrm2

lphi=0.d0

ist=1
jst=1
ilrold=-1
do iorb=1,orbs%norbp
    iiorb=orbs%isorb+iorb
    ilr=onWhichAtom(iiorb)
    if(ilr==ilrold) then
        ! Set back the index of lphiovrlp, since we again need the same orbitals.
        jst=jst-op%noverlaps(iiorb)*ncount
    end if
    ncount=lzd%llr(ilr)%wfd%nvctr_c+7*lzd%llr(ilr)%wfd%nvctr_f
    do jorb=1,op%noverlaps(iiorb)
        jjorb=op%overlaps(jorb,iiorb)
        call daxpy(ncount, ovrlp(jjorb,iiorb), lphiovrlp(jst), 1, lphi(ist), 1)
        jst=jst+ncount
    end do

    ! Normalize
    tt=dnrm2(ncount, lphi(ist), 1)
    call dscal(ncount, 1/tt, lphi(ist), 1)

    ist=ist+ncount
    ilrold=ilr

end do


end subroutine globalLoewdin




subroutine applyOrthoconstraint(iproc, nproc, orbs, lorbs, onWhichAtom, lzd, op, lagmat, lphiovrlp, lhphi)
use module_base
use module_types
implicit none

! Calling arguments
integer,intent(in):: iproc, nproc
type(orbitals_data),intent(in):: orbs, lorbs
integer,dimension(orbs%norb),intent(in):: onWhichAtom
type(linear_zone_descriptors),intent(in):: lzd
type(overlapParameters),intent(in):: op
real(8),dimension(orbs%norb,orbs%norb),intent(in):: lagmat
real(8),dimension(op%ndim_lphiovrlp),intent(in):: lphiovrlp
real(8),dimension(lorbs%npsidim),intent(out):: lhphi

! Local variables
integer:: iorb, jorb, iiorb, ilr, ist, jst, ilrold, jjorb, ncount



ist=1
jst=1
ilrold=-1
do iorb=1,orbs%norbp
    iiorb=orbs%isorb+iorb
    ilr=onWhichAtom(iiorb)
    if(ilr==ilrold) then
        ! Set back the index of lphiovrlp, since we again need the same orbitals.
        jst=jst-op%noverlaps(iiorb)*ncount
    end if
    ncount=lzd%llr(ilr)%wfd%nvctr_c+7*lzd%llr(ilr)%wfd%nvctr_f
    do jorb=1,op%noverlaps(iiorb)
        jjorb=op%overlaps(jorb,iiorb)
        call daxpy(ncount, -.5d0*lagmat(jjorb,iiorb), lphiovrlp(jst), 1, lhphi(ist), 1)
        call daxpy(ncount, -.5d0*lagmat(iiorb,jjorb), lphiovrlp(jst), 1, lhphi(ist), 1)
        jst=jst+ncount
    end do
    ist=ist+ncount
    ilrold=ilr
end do



end subroutine applyOrthoconstraint



subroutine applyOrthoconstraintNonorthogonal(iproc,nproc,orbs,lorbs,onWhichAtom,lzd,op,lagmat,ovrlp,lphiovrlp,lhphiovrlp,lhphi)
use module_base
use module_types
implicit none

! Calling arguments
integer,intent(in):: iproc, nproc
type(orbitals_data),intent(in):: orbs, lorbs
integer,dimension(orbs%norb),intent(in):: onWhichAtom
type(linear_zone_descriptors),intent(in):: lzd
type(overlapParameters),intent(in):: op
real(8),dimension(orbs%norb,orbs%norb),intent(in):: lagmat, ovrlp
real(8),dimension(op%ndim_lphiovrlp),intent(in):: lphiovrlp, lhphiovrlp
real(8),dimension(lorbs%npsidim),intent(out):: lhphi

! Local variables
integer:: iorb, jorb, iiorb, ilr, ist, jst, ilrold, jjorb, ncount



ist=1
jst=1
ilrold=-1
do iorb=1,orbs%norbp
    iiorb=orbs%isorb+iorb
    ilr=onWhichAtom(iiorb)
    if(ilr==ilrold) then
        ! Set back the index of lphiovrlp, since we again need the same orbitals.
        jst=jst-op%noverlaps(iiorb)*ncount
    end if
    ncount=lzd%llr(ilr)%wfd%nvctr_c+7*lzd%llr(ilr)%wfd%nvctr_f
    call dscal(ncount, 1.5d0, lhphi(ist), 1)
    do jorb=1,op%noverlaps(iiorb)
        jjorb=op%overlaps(jorb,iiorb)
        call daxpy(ncount, -.5d0*ovrlp(jjorb,iiorb), lhphiovrlp(jst), 1, lhphi(ist), 1)
        call daxpy(ncount, -.5d0*lagmat(jjorb,iiorb), lphiovrlp(jst), 1, lhphi(ist), 1)
        call daxpy(ncount, -.5d0*lagmat(iiorb,jjorb), lphiovrlp(jst), 1, lhphi(ist), 1)
        jst=jst+ncount
    end do
    ist=ist+ncount
    ilrold=ilr
end do



end subroutine applyOrthoconstraintNonorthogonal




subroutine applyOrthoconstraintNonorthogonal2(iproc,nproc,orbs,lorbs,onWhichAtom,lzd,op,lagmat,ovrlp,lphiovrlp,lhphiovrlp,lhphi)
use module_base
use module_types
implicit none

! Calling arguments
integer,intent(in):: iproc, nproc
type(orbitals_data),intent(in):: orbs, lorbs
integer,dimension(orbs%norb),intent(in):: onWhichAtom
type(linear_zone_descriptors),intent(in):: lzd
type(overlapParameters),intent(in):: op
real(8),dimension(orbs%norb,orbs%norb),intent(in):: lagmat, ovrlp
real(8),dimension(op%ndim_lphiovrlp),intent(in):: lphiovrlp, lhphiovrlp
real(8),dimension(lorbs%npsidim),intent(out):: lhphi

! Local variables
integer:: iorb, jorb, iiorb, ilr, ist, jst, ilrold, jjorb, ncount, info, i, istat, iall
real(8):: tt
real(8),dimension(:,:),allocatable:: ovrlp2
real(8),dimension(:,:),allocatable:: ovrlp_minus_one_lagmat, ovrlp_minus_one_lagmat_trans
character(len=*),parameter:: subname='applyOrthoconstraintNonorthogonal2'

allocate(ovrlp_minus_one_lagmat(orbs%norb,orbs%norb), stat=istat)
call memocc(istat, ovrlp_minus_one_lagmat, 'ovrlp_minus_one_lagmat', subname)
allocate(ovrlp_minus_one_lagmat_trans(orbs%norb,orbs%norb), stat=istat)
call memocc(istat, ovrlp_minus_one_lagmat_trans, 'ovrlp_minus_one_lagmat_trans', subname)
allocate(ovrlp2(orbs%norb,orbs%norb), stat=istat)
call memocc(istat, orbs%norb, 'orbs%norb', subname)

call dcopy(orbs%norb**2, ovrlp(1,1), 1, ovrlp2(1,1), 1)
! Invert the overlap matrix
call dpotrf('l', orbs%norb, ovrlp2(1,1), orbs%norb, info)
if(info/=0) then
    write(*,'(x,a,i0)') 'ERROR in dpotrf, info=',info
    stop
end if
call dpotri('l', orbs%norb, ovrlp2(1,1), orbs%norb, info)
if(info/=0) then
    write(*,'(x,a,i0)') 'ERROR in dpotri, info=',info
    stop
end if


! Multiply the Lagrange multiplier matrix with S^-1/2.
! First fill the upper triangle.
do iorb=1,orbs%norb
    do jorb=1,iorb-1
        ovrlp2(jorb,iorb)=ovrlp2(iorb,jorb)
    end do
end do
call dgemm('n', 'n', orbs%norb, orbs%norb, orbs%norb, 1.d0, ovrlp2(1,1), orbs%norb, lagmat(1,1), orbs%norb, &
     0.d0, ovrlp_minus_one_lagmat(1,1), orbs%norb)
call dgemm('n', 't', orbs%norb, orbs%norb, orbs%norb, 1.d0, ovrlp2(1,1), orbs%norb, lagmat(1,1), orbs%norb, &
     0.d0, ovrlp_minus_one_lagmat_trans(1,1), orbs%norb)


ist=1
jst=1
ilrold=-1
do iorb=1,orbs%norbp
    iiorb=orbs%isorb+iorb
    ilr=onWhichAtom(iiorb)
    if(ilr==ilrold) then
        ! Set back the index of lphiovrlp, since we again need the same orbitals.
        jst=jst-op%noverlaps(iiorb)*ncount
    end if
    ncount=lzd%llr(ilr)%wfd%nvctr_c+7*lzd%llr(ilr)%wfd%nvctr_f
    !call dscal(ncount, 1.5d0, lhphi(ist), 1)
    do jorb=1,op%noverlaps(iiorb)
        jjorb=op%overlaps(jorb,iiorb)
        !call daxpy(ncount, -.5d0*ovrlp(jjorb,iiorb), lhphiovrlp(jst), 1, lhphi(ist), 1)
        call daxpy(ncount, -.5d0*ovrlp_minus_one_lagmat(jjorb,iiorb), lphiovrlp(jst), 1, lhphi(ist), 1)
        call daxpy(ncount, -.5d0*ovrlp_minus_one_lagmat_trans(jjorb,iiorb), lphiovrlp(jst), 1, lhphi(ist), 1)
        jst=jst+ncount
    end do
    ist=ist+ncount
    ilrold=ilr
end do


iall=-product(shape(ovrlp_minus_one_lagmat))*kind(ovrlp_minus_one_lagmat)
deallocate(ovrlp_minus_one_lagmat, stat=istat)
call memocc(istat, iall, 'ovrlp_minus_one_lagmat', subname)
iall=-product(shape(ovrlp_minus_one_lagmat_trans))*kind(ovrlp_minus_one_lagmat_trans)
deallocate(ovrlp_minus_one_lagmat_trans, stat=istat)
call memocc(istat, iall, 'ovrlp_minus_one_lagmat_trans', subname)
iall=-product(shape(ovrlp2))*kind(ovrlp2)
deallocate(ovrlp2, stat=istat)
call memocc(istat, iall, 'ovrlp2', subname)


end subroutine applyOrthoconstraintNonorthogonal2






subroutine transformOverlapMatrix2(iproc, nproc, orbs, ovrlp)
use module_base
use module_types
implicit none

! Calling arguments
integer,intent(in):: iproc, nproc
type(orbitals_data),intent(in):: orbs
real(8),dimension(orbs%norb,orbs%norb),intent(inout):: ovrlp

! Local variables
integer:: lwork, istat, iall, iorb, jorb, info
real(8),dimension(:),allocatable:: work, eval
real(8),dimension(:,:,:),allocatable:: tempArr
character(len=*),parameter:: subname='transformOverlapMatrix'


allocate(eval(orbs%norb), stat=istat)
call memocc(istat, eval, 'eval', subname)
allocate(tempArr(orbs%norb,orbs%norb,2), stat=istat)
call memocc(istat, tempArr, 'tempArr', subname)

do iorb=1,orbs%norb
    do jorb=1,orbs%norb
        if(iorb==jorb) then
            ovrlp(iorb,jorb)=1.5d0-.5d0*ovrlp(iorb,jorb)
        else
            ovrlp(iorb,jorb)=-.5d0*ovrlp(iorb,jorb)
        end if
    end do
end do



endsubroutine transformOverlapMatrix2





subroutine indicesForExpansion(iproc, nproc, orbs, input, onWhichAtom, lzd, op, comon)
use module_base
use module_types
implicit none

! Calling arguments
integer,intent(in):: iproc, nproc
type(orbitals_data),intent(in):: orbs
type(input_variables),intent(in):: input
integer,dimension(orbs%norb),intent(in):: onWhichAtom
type(linear_zone_descriptors),intent(in):: lzd
type(overlapParameters),intent(in):: op
type(p2pCommsOrthonormality),intent(in):: comon

! Local variables
integer:: ind, iorb, iiorb, ilr, gdim, ldim, jorb, jjorb, jst, ilrold



ind=1
ilrold=-1
do iorb=1,orbs%norbp
    iiorb=orbs%isorb+iorb
    ilr=onWhichAtom(iiorb)
    if(ilr==ilrold) cycle
    gdim=lzd%llr(ilr)%wfd%nvctr_c+7*lzd%llr(ilr)%wfd%nvctr_f
    do jorb=1,op%noverlaps(iiorb)
        jjorb=op%overlaps(jorb,iiorb)
        ! Starting index of orbital jjorb
        jst=op%indexInRecvBuf(iorb,jjorb)
        !ldim=op%olr(jorb,iiorb)%wfd%nvctr_c+7*op%olr(jorb,iiorb)%wfd%nvctr_f
        ldim=op%olr(jorb,iorb)%wfd%nvctr_c+7*op%olr(jorb,iorb)%wfd%nvctr_f
        !call Lpsi_to_global2(iproc, nproc, ldim, gdim, orbs%norbp, orbs%nspinor, input%nspin, lzd%llr(ilr), op%olr(jorb,iiorb), comon%recvBuf(jst), lphiovrlp(ind))
        !call Lpsi_to_global2(iproc, nproc, ldim, gdim, orbs%norbp, orbs%nspinor, input%nspin, lzd%llr(ilr), op%olr(jorb,iorb), comon%recvBuf(jst), lphiovrlp(ind))
        call index_of_Lpsi_to_global2(iproc, nproc, ldim, gdim, orbs%norbp, orbs%nspinor, input%nspin, lzd%llr(ilr),&
             op%olr(jorb,iorb), op%indexExpand(jst))
        ind=ind+gdim
    end do
    ilrold=ilr
end do

end subroutine indicesForExpansion





subroutine indicesForExtraction(iproc, nproc, orbs, sizePhi, onWhichAtom, lzd, op, comon)
use module_base
use module_types
implicit none

! Calling arguments
integer,intent(in):: iproc, nproc, sizePhi
type(orbitals_data),intent(in):: orbs
integer,dimension(orbs%norb),intent(in):: onWhichAtom
type(linear_zone_descriptors),intent(in):: lzd
type(overlapParameters),intent(inout):: op
type(p2pCommsOrthonormality),intent(out):: comon

! Local variables
integer:: iorb, jorb, korb, ind, indovrlp, ilr, klr, ilrold, jjorb, jjlr, jjproc, iiproc, iiprocold, gdim, ldim, kkorb, lorb
integer:: i

indovrlp=1
op%indexInSendBuf=0

ilrold=-1
do iorb=1,orbs%norb
    ilr=onWhichAtom(iorb)
    iiproc=orbs%onWhichMPI(iorb)
    if(ilr==ilrold .and. iiproc==iiprocold) cycle ! otherwise we would extract the same again
    do jorb=1,op%noverlaps(iorb)
        jjorb=op%overlaps(jorb,iorb)
        jjlr=onWhichAtom(jjorb)
        jjproc=orbs%onWhichMPI(jjorb)
        if(iproc==jjproc) then
            ! Get the correct descriptors
            korb=jjorb-orbs%isorb
            !write(*,'(a,5i8)') 'iorb, jorb, jjorb, jjproc, korb', iorb, jorb, jjorb, jjproc, korb
            do i=1,op%noverlaps(jjorb)
                !write(*,'(a,5i8)') 'iproc, iorb, korb, i, op%overlaps(i,korb)', iproc, iorb, korb, i, op%overlaps(i,korb)
                if(op%overlaps(i,jjorb)==iorb) then
                    lorb=i
                    exit
                end if
            end do
            !write(*,'(a,5i9)') 'iproc, iorb, jorb, korb, lorb', iproc, iorb, jorb, korb, lorb
            gdim=lzd%llr(jjlr)%wfd%nvctr_c+7*lzd%llr(jjlr)%wfd%nvctr_f
            ldim=op%olr(lorb,korb)%wfd%nvctr_c+7*op%olr(lorb,korb)%wfd%nvctr_f
            ind=1
            do kkorb=orbs%isorb+1,jjorb-1
                klr=onWhichAtom(kkorb)
                ind = ind + lzd%llr(klr)%wfd%nvctr_c + 7*lzd%llr(klr)%wfd%nvctr_f
            end do
            !write(*,'(5(a,i0))') 'process ',iproc,' adds ',op%olr(lorb,korb)%wfd%nvctr_c+7*op%olr(lorb,korb)%wfd%nvctr_f,' elements at position ',indovrlp,' from orbital ',jjorb,' for orbital ', iorb
            !call psi_to_locreg2(iproc, nproc, ldim, gdim, op%olr(lorb,korb), lzd%llr(jjlr), phi(ind), comon%sendBuf(indovrlp))
            call index_of_psi_to_locreg2(iproc, nproc, ldim, gdim, op%olr(lorb,korb), lzd%llr(jjlr), op%indexExtract(indovrlp))
            op%indexInSendBuf(jjorb-orbs%isorb,iorb)=indovrlp
            indovrlp=indovrlp+op%olr(lorb,korb)%wfd%nvctr_c+7*op%olr(lorb,korb)%wfd%nvctr_f
        end if
    end do
    ilrold=ilr
    iiprocold=iiproc
end do

if(indovrlp/=comon%nsendBuf+1) then
    write(*,'(x,a,i0,a,3x,i0,2x,i0)') 'ERROR on process ', iproc, ': indovrlp/=comon%nsendBuf+1', indovrlp, comon%nsendBuf+1
    stop
end if



end subroutine indicesForExtraction






subroutine getMatrixElements2(iproc, nproc, lin, lphi, lhphi, matrixElements)
use module_base
use module_types
use module_interfaces, exceptThisOne => getMatrixElements2
implicit none

! Calling arguments
integer,intent(in):: iproc, nproc
type(linearParameters),intent(inout):: lin
real(8),dimension(lin%lb%lzd%orbs%npsidim),intent(in):: lphi, lhphi
real(8),dimension(lin%lb%lzd%orbs%norb,lin%lb%lzd%orbs%norb),intent(out):: matrixElements

! Local variables
integer:: it, istat, iall, iorb
real(8),dimension(:),allocatable:: lphiovrlp
character(len=*),parameter:: subname='getMatrixElements'



  allocate(lphiovrlp(lin%op_lb%ndim_lphiovrlp), stat=istat)
  call memocc(istat, lphiovrlp, 'lphiovrlp',subname)


  ! Put lphi in the sendbuffer, i.e. lphi will be sent to other processes' receive buffer.
  !call extractOrbital(iproc, nproc, lin%orbs, lin%lorbs%npsidim, lin%onWhichAtomAll, lin%lzd, lin%op, lphi, lin%comon)
  call extractOrbital2(iproc,nproc,lin%lb%lorbs,lin%lb%lorbs%npsidim,lin%lb%onWhichAtomAll,lin%lb%lzd,lin%op_lb,lphi,lin%comon_lb)
  call postCommsOverlap(iproc, nproc, lin%comon_lb)
  !call gatherOrbitals(iproc, nproc, lin%comon)
  call gatherOrbitals2(iproc, nproc, lin%comon_lb)
  ! Put lhphi to the sendbuffer, so we can the calculate <lphi|lhphi>
  !call extractOrbital(iproc, nproc, lin%orbs, lin%lorbs%npsidim, lin%onWhichAtomAll, lin%lzd, lin%op, lhphi, lin%comon)
  call extractOrbital2(iproc,nproc,lin%lb%lorbs,lin%lb%lorbs%npsidim,lin%lb%onWhichAtomAll,lin%lb%lzd,lin%op_lb,lhphi,lin%comon_lb)
  call calculateOverlapMatrix2(iproc, nproc, lin%lb%lzd%orbs, lin%op_lb, lin%comon_lb, lin%lb%onWhichAtomAll, matrixElements)

  !!!allocate(lphiovrlp(lin%op_lb%ndim_lphiovrlp), stat=istat)
  !!!call memocc(istat, lphiovrlp, 'lphiovrlp',subname)

  !!!call extractOrbital2(iproc, nproc, lin%lb%orbs, lin%lb%lorbs%npsidim, lin%lb%onWhichAtomAll, lin%lb%lzd, lin%op_lb, lphi, lin%comon_lb)
  !!!call postCommsOverlap(iproc, nproc, lin%comon_lb)
  !!!call gatherOrbitals2(iproc, nproc, lin%comon_lb)
  !!!call calculateOverlapMatrix2(iproc, nproc, lin%lb%lzd%orbs, lin%op_lb, lin%comon_lb, lin%lb%onWhichAtomAll, ovrlp)

  !!!iall=-product(shape(lphiovrlp))*kind(lphiovrlp)
  !!!deallocate(lphiovrlp, stat=istat)
  !!!call memocc(istat, iall, 'lphiovrlp', subname)

end subroutine getMatrixElements2
