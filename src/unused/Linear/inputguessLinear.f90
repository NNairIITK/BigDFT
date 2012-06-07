!!subroutine orthonormalizeCoefficients(orbs, orbsig, coeff)
!!use module_base
!!use module_types
!!implicit none
!!
!!! Calling arguments
!!type(orbitals_data),intent(in):: orbs, orbsig
!!real(8),dimension(orbsig%norb,orbs%norb),intent(inout):: coeff
!!
!!! Local variables
!!integer:: iorb, jorb, istat, iall, lwork, info
!!real(8),dimension(:),allocatable:: work, eval
!!real(8),dimension(:,:),allocatable:: ovrlp, coeffTemp
!!real(8),dimension(:,:,:),allocatable:: tempArr
!!character(len=*),parameter:: subname='orthonormalizeCoefficients'
!!real(8):: ddot
!!
!!        allocate(ovrlp(orbs%norb,orbs%norb), stat=istat)
!!        call memocc(istat, ovrlp, 'ovrlp', subname)
!!        allocate(eval(orbs%norb), stat=istat)
!!        call memocc(istat, eval, 'eval', subname)
!!        allocate(tempArr(orbs%norb,orbs%norb,2), stat=istat)
!!        call memocc(istat, tempArr, 'tempArr', subname)
!!        allocate(coeffTemp(orbsig%norb,orbs%norb), stat=istat)
!!        call memocc(istat, coeffTemp, 'coeffTemp', subname)
!!
!!
!!        !!! Orthonormalize the coefficient vectors (Gram-Schmidt).
!!        !!do iorb=1,orbs%norb
!!        !!    do jorb=1,iorb-1
!!        !!        tt=ddot(orbsig%norb, coeff(1,iorb), 1, coeff(1,jorb), 1)
!!        !!        call daxpy(orbsig%norb, -tt, coeff(1,jorb), 1, coeff(1,iorb), 1)
!!        !!    end do
!!        !!    tt=dnrm2(orbsig%norb, coeff(1,iorb), 1)
!!        !!    call dscal(orbsig%norb, 1/tt, coeff(1,iorb), 1)
!!        !!end do
!!
!!        !!! Orthonormalize the coefficient vectors (Loewdin).
!!        !!do iorb=1,orbs%norb
!!        !!    do jorb=1,orbs%norb
!!        !!        ovrlp(iorb,jorb)=ddot(orbsig%norb, coeff(1,iorb), 1, coeff(1,jorb), 1)
!!        !!    end do
!!        !!end do
!!
!!        allocate(work(1), stat=istat)
!!        call memocc(istat, work, 'work', subname)
!!        call dsyev('v', 'l', orbs%norb, ovrlp(1,1), orbs%norb, eval, work, -1, info)
!!        lwork=work(1)
!!        iall=-product(shape(work))*kind(work)
!!        deallocate(work, stat=istat)
!!        call memocc(istat, iall, 'work', subname)
!!        allocate(work(lwork), stat=istat)
!!        call memocc(istat, work, 'work', subname)
!!        call dsyev('v', 'l', orbs%norb, ovrlp(1,1), orbs%norb, eval, work, lwork, info)
!!        iall=-product(shape(work))*kind(work)
!!        deallocate(work, stat=istat)
!!        call memocc(istat, iall, 'work', subname)
!!
!!        ! Calculate S^{-1/2}. 
!!        ! First calulate ovrlp*diag(1/sqrt(evall)) (ovrlp is the diagonalized overlap
!!        ! matrix and diag(1/sqrt(evall)) the diagonal matrix consisting of the inverse square roots of the eigenvalues...
!!        do iorb=1,orbs%norb
!!            do jorb=1,orbs%norb
!!                tempArr(jorb,iorb,1)=ovrlp(jorb,iorb)*1.d0/sqrt(eval(iorb))
!!            end do
!!        end do
!!
!!        ! ...and now apply the diagonalized overlap matrix to the matrix constructed above.
!!        ! This will give S^{-1/2}.
!!        call dgemm('n', 't', orbs%norb, orbs%norb, orbs%norb, 1.d0, ovrlp(1,1), &
!!        orbs%norb, tempArr(1,1,1), orbs%norb, 0.d0, &
!!        tempArr(1,1,2), orbs%norb)
!!
!!        ! Now calculate the orthonormal orbitals by applying S^{-1/2} to the orbitals.
!!        ! This requires the use of a temporary variable phidTemp.
!!        call dgemm('n', 'n', orbsig%norb, orbs%norb, orbs%norb, 1.d0, coeff(1,1), &
!!             orbsig%norb, tempArr(1,1,2), orbs%norb, 0.d0, &
!!             coeffTemp(1,1), orbsig%norb)
!!        
!!        ! Now copy the orbitals from the temporary variable to phid.
!!        call dcopy(orbs%norb*orbsig%norb, coeffTemp(1,1), 1, coeff(1,1), 1)
!!
!!        iall=-product(shape(ovrlp))*kind(ovrlp)
!!        deallocate(ovrlp, stat=istat)
!!        call memocc(istat, iall, 'ovrlp', subname)
!!
!!        iall=-product(shape(eval))*kind(eval)
!!        deallocate(eval, stat=istat)
!!        call memocc(istat, iall, 'eval', subname)
!!
!!        iall=-product(shape(tempArr))*kind(tempArr)
!!        deallocate(tempArr, stat=istat)
!!        call memocc(istat, iall, 'tempArr', subname)
!!
!!        iall=-product(shape(coeffTemp))*kind(coeffTemp)
!!        deallocate(coeffTemp, stat=istat)
!!        call memocc(istat, iall, 'coeffTemp', subname)
!!
!!end subroutine orthonormalizeCoefficients

!!!subroutine postCommsVectorOrthonormalization(iproc, nproc, newComm, comom)
!!!use module_base
!!!use module_types
!!!implicit none
!!!
!!!! Calling arguments
!!!integer,intent(in):: iproc, nproc, newComm
!!!type(p2pCommsOrthonormalityMatrix),intent(inout):: comom
!!!
!!!! Local variables
!!!integer:: nsends, nreceives, jproc, iorb, mpisource, istsource, ncount, mpidest, istdest, tag, ierr
!!!
!!!nsends=0
!!!nreceives=0
!!!comom%communComplete=.false.
!!!do jproc=0,nproc-1
!!!  do iorb=1,comom%noverlapProc(jproc)
!!!     mpisource=comom%comarr(1,iorb,jproc)
!!!     istsource=comom%comarr(2,iorb,jproc)
!!!     ncount=comom%comarr(3,iorb,jproc)
!!!     mpidest=comom%comarr(4,iorb,jproc)
!!!     istdest=comom%comarr(5,iorb,jproc)
!!!     tag=comom%comarr(6,iorb,jproc)
!!!     !if(iproc==0) write(*,'(a,4i9)') 'jproc, iorb, mpisource, mpidest', jproc, iorb, mpisource, mpidest
!!!     if(mpisource/=mpidest) then
!!!        ! The orbitals are on different processes, so we need a point to point communication.
!!!        if(iproc==mpisource) then
!!!           !write(*,'(6(a,i0))') 'process ', mpisource, ' sends ', ncount, ' elements from position ', istsource, ' to position ', istdest, ' on process ', mpidest, ', tag=',tag
!!!           call mpi_isend(comom%sendBuf(istsource), ncount, mpi_double_precision, mpidest, tag, newComm,&
!!!                comom%comarr(7,iorb,jproc), ierr)
!!!           !call mpi_isend(sendBuf(istsource), ncount, mpi_double_precision, mpidest, tag, mpi_comm_world, lin%comsr%comarr(8,iorb,jproc), ierr)
!!!           comom%comarr(8,iorb,jproc)=mpi_request_null !is this correct?
!!!           nsends=nsends+1
!!!        else if(iproc==mpidest) then
!!!           !write(*,'(6(a,i0))') 'Aprocess ', mpidest, ' receives ', ncount,&
!!!           !     ' elements at position ', istdest, ' from position ', istsource, ' on process ', mpisource, ', tag=',tag
!!!           call mpi_irecv(comom%recvBuf(istdest), ncount, mpi_double_precision, mpisource, tag, newComm,&
!!!                comom%comarr(8,iorb,jproc), ierr)
!!!           comom%comarr(7,iorb,jproc)=mpi_request_null !is this correct?
!!!           nreceives=nreceives+1
!!!        else
!!!           comom%comarr(7,iorb,jproc)=mpi_request_null
!!!           comom%comarr(8,iorb,jproc)=mpi_request_null
!!!        end if
!!!     else
!!!        ! The orbitals are on the same process, so simply copy them.
!!!        if(iproc==mpisource) then
!!!           call dcopy(ncount, comom%sendBuf(istsource), 1, comom%recvBuf(istdest), 1)
!!!           !write(*,'(6(a,i0))') 'process ', iproc, ' copies ', ncount, ' elements from position ', istsource, ' to position ', istdest, ' on process ', iproc, ', tag=',tag
!!!           comom%comarr(7,iorb,jproc)=mpi_request_null
!!!           comom%comarr(8,iorb,jproc)=mpi_request_null
!!!           nsends=nsends+1
!!!           nreceives=nreceives+1
!!!           comom%communComplete(iorb,iproc)=.true.
!!!        else
!!!           comom%comarr(7,iorb,jproc)=mpi_request_null
!!!           comom%comarr(8,iorb,jproc)=mpi_request_null
!!!        end if
!!!
!!!     end if
!!!  end do
!!!end do
!!!
!!!
!!!
!!!end subroutine postCommsVectorOrthonormalization

!!!!subroutine gatherVectors(iproc, nproc, newComm, comom)
!!!!use module_base
!!!!use module_types
!!!!implicit none
!!!!
!!!!! Calling arguments
!!!!integer,intent(in):: iproc, nproc, newComm
!!!!type(p2pCommsOrthonormalityMatrix),intent(inout):: comom
!!!!
!!!!! Local variables
!!!!integer:: jorb, mpisource, mpidest, nfast, nslow, nsameproc, ierr, jproc
!!!!integer,dimension(mpi_status_size):: stat
!!!!logical:: sendComplete, receiveComplete
!!!!
!!!!!!! Check whether the communications have completed.
!!!!!!nfast=0
!!!!!!nsameproc=0
!!!!!!testLoop: do
!!!!!!    do jproc=0,nproc-1
!!!!!!        do jorb=1,comom%noverlapProc(jproc)
!!!!!!            if(comom%communComplete(jorb,jproc)) cycle
!!!!!!            call mpi_test(comom%comarr(7,jorb,jproc), sendComplete, stat, ierr)
!!!!!!            call mpi_test(comom%comarr(8,jorb,jproc), receiveComplete, stat, ierr)
!!!!!!            if(sendComplete .and. receiveComplete) comom%communComplete(jorb,jproc)=.true.
!!!!!!            if(comom%communComplete(jorb,jproc)) then
!!!!!!                !write(*,'(2(a,i0))') 'fast communication; process ', iproc, ' has received orbital ', jorb
!!!!!!                mpisource=comom%comarr(1,jorb,jproc)
!!!!!!                mpidest=comom%comarr(4,jorb,jproc)
!!!!!!                if(mpisource/=mpidest) then
!!!!!!                    nfast=nfast+1
!!!!!!                else
!!!!!!                    nsameproc=nsameproc+1
!!!!!!                end if
!!!!!!            end if
!!!!!!        end do
!!!!!!    end do
!!!!!!    ! If we made it until here, either all all the communication is
!!!!!!    ! complete or we better wait for each single orbital.
!!!!!!    exit testLoop
!!!!!!end do testLoop
!!!!
!!!!
!!!!! Wait for the communications that have not completed yet
!!!!nslow=0
!!!!do jproc=0,nproc-1
!!!!  do jorb=1,comom%noverlapProc(jproc)
!!!!     if(comom%communComplete(jorb,jproc)) cycle
!!!!     !write(*,'(2(a,i0))') 'process ', iproc, ' is waiting for orbital ', korb
!!!!     nslow=nslow+1
!!!!     call mpi_wait(comom%comarr(7,jorb,jproc), stat, ierr)   !COMMENTED BY PB
!!!!     call mpi_wait(comom%comarr(8,jorb,jproc), stat, ierr)   !COMMENTED BY PB
!!!!     comom%communComplete(jorb,jproc)=.true.
!!!!  end do
!!!!end do
!!!!
!!!!!call mpiallred(nreceives, 1, mpi_sum, mpi_comm_world, ierr)
!!!!!call mpiallred(nfast, 1, mpi_sum, newComm, ierr)
!!!!!call mpiallred(nslow, 1, mpi_sum, newComm, ierr)
!!!!!call mpiallred(nsameproc, 1, mpi_sum, newComm, ierr)
!!!!!if(iproc==0) write(*,'(1x,2(a,i0),a)') 'statistics: - ', nfast+nslow, ' point to point communications, of which ', &
!!!!!                       nfast, ' could be overlapped with computation.'
!!!!!if(iproc==0) write(*,'(1x,a,i0,a)') '            - ', nsameproc, ' copies on the same processor.'
!!!!
!!!!
!!!!end subroutine gatherVectors

!!!!!!subroutine applyOrthoconstraintVectors(iproc, nproc, methTransformOverlap, correctionOrthoconstraint, blocksize_pdgemm, &
!!!!!!           comm, norb, norbmax, norbp, isorb, nlr, noverlaps, onWhichAtom, ovrlp, &
!!!!!!           lagmat, comom, mlr, mad, orbs, grad, ovrlp_minus_one_lagmat, ovrlp_minus_one_lagmat_trans)
!!!!!!use module_base
!!!!!!use module_types
!!!!!!use module_interfaces, exceptThisOne => applyOrthoconstraintVectors
!!!!!!implicit none
!!!!!!
!!!!!!! Calling arguments
!!!!!!integer,intent(in):: iproc, nproc, methTransformOverlap, correctionOrthoconstraint, blocksize_pdgemm
!!!!!!integer,intent(in):: comm, norb, norbmax, norbp, isorb, nlr, noverlaps
!!!!!!integer,dimension(norb),intent(in):: onWhichAtom
!!!!!!real(8),dimension(norb,norb),intent(in):: ovrlp
!!!!!!real(8),dimension(norb,norb),intent(inout):: lagmat
!!!!!!type(p2pCommsOrthonormalityMatrix),intent(in):: comom
!!!!!!type(matrixLocalizationRegion),dimension(nlr),intent(in):: mlr
!!!!!!type(matrixDescriptors),intent(in):: mad
!!!!!!type(orbitals_data),intent(in):: orbs
!!!!!!real(8),dimension(norbmax,norbp),intent(inout):: grad
!!!!!!real(8),dimension(norb,norb),intent(out):: ovrlp_minus_one_lagmat, ovrlp_minus_one_lagmat_trans
!!!!!!
!!!!!!! Local variables
!!!!!!integer:: info, iorb, ilrold, iiorb, jjorb, ilr, ncount, jorb, ijorb, istat, iall
!!!!!!real(8),dimension(:,:),allocatable:: ovrlp2
!!!!!!real(8):: tt
!!!!!!character(len=*),parameter:: subname='applyOrthoconstraintVectors'
!!!!!!
!!!!!!
!!!!!!allocate(ovrlp2(norb,norb), stat=istat)
!!!!!!call memocc(istat, ovrlp2, 'ovrlp2', subname)
!!!!!!
!!!!!!call dcopy(norb**2, ovrlp(1,1), 1, ovrlp2(1,1), 1)
!!!!!!
!!!!!!correctionIf: if(correctionOrthoconstraint==0) then
!!!!!!    ! Invert the overlap matrix
!!!!!!    call overlapPowerMinusOne(iproc, nproc, methTransformOverlap, norb, mad, orbs, ovrlp2)
!!!!!!    
!!!!!!    ! Multiply the Lagrange multiplier matrix with S^-1/2.
!!!!!!    ! First fill the upper triangle.
!!!!!!    do iorb=1,norb
!!!!!!        do jorb=1,iorb-1
!!!!!!            ovrlp2(jorb,iorb)=ovrlp2(iorb,jorb)
!!!!!!        end do
!!!!!!    end do
!!!!!!    if(blocksize_pdgemm<0) then
!!!!!!        !!call dgemm('n', 'n', norb, norb, norb, 1.d0, ovrlp2(1,1), norb, lagmat(1,1), norb, &
!!!!!!        !!     0.d0, ovrlp_minus_one_lagmat(1,1), norb)
!!!!!!        !!call dgemm('n', 't', norb, norb, norb, 1.d0, ovrlp2(1,1), norb, lagmat(1,1), norb, &
!!!!!!        !!     0.d0, ovrlp_minus_one_lagmat_trans(1,1), norb)
!!!!!!        !!call dsymm('l', 'l', norb, norb, 1.d0, ovrlp2(1,1), norb, lagmat(1,1), norb, &
!!!!!!        !!     0.d0, ovrlp_minus_one_lagmat(1,1), norb)
!!!!!!        ovrlp_minus_one_lagmat=0.d0
!!!!!!        !!call dgemm_compressed2(iproc, nproc, norb, mad%nsegline, mad%nseglinemax, mad%keygline, mad%nsegmatmul, mad%keygmatmul, &
!!!!!!        !!     ovrlp2, lagmat, ovrlp_minus_one_lagmat)
!!!!!!        call dgemm_compressed_parallel(iproc, nproc, norb, mad%nsegline, mad%nseglinemax,&
!!!!!!             mad%keygline, mad%nsegmatmul, mad%keygmatmul, &
!!!!!!             orbs%norb_par, orbs%isorb_par, orbs%norbp, ovrlp2, lagmat, ovrlp_minus_one_lagmat)
!!!!!!        !ovrlp_minus_one_lagmat=lagmat
!!!!!!        ! Transpose lagmat
!!!!!!        do iorb=1,norb
!!!!!!            do jorb=iorb+1,norb
!!!!!!                tt=lagmat(jorb,iorb)
!!!!!!                lagmat(jorb,iorb)=lagmat(iorb,jorb)
!!!!!!                lagmat(iorb,jorb)=tt
!!!!!!            end do
!!!!!!        end do
!!!!!!        !!call dsymm('l', 'l', norb, norb, 1.d0, ovrlp2(1,1), norb, lagmat(1,1), norb, &
!!!!!!        !!     0.d0, ovrlp_minus_one_lagmat_trans(1,1), norb)
!!!!!!        !ovrlp_minus_one_lagmat_trans=lagmat
!!!!!!        ovrlp_minus_one_lagmat_trans=0.d0
!!!!!!        !!call dgemm_compressed2(iproc, nproc, norb, mad%nsegline, mad%nseglinemax, mad%keygline, mad%nsegmatmul, mad%keygmatmul, &
!!!!!!        !!     ovrlp2, lagmat, ovrlp_minus_one_lagmat_trans)
!!!!!!        call dgemm_compressed_parallel(iproc, nproc, norb, mad%nsegline, mad%nseglinemax,&
!!!!!!             mad%keygline, mad%nsegmatmul, mad%keygmatmul, &
!!!!!!             orbs%norb_par, orbs%isorb_par, orbs%norbp, ovrlp2, lagmat, ovrlp_minus_one_lagmat_trans)
!!!!!!    
!!!!!!    else
!!!!!!        call dsymm_parallel(iproc, nproc, blocksize_pdgemm, comm, 'l', 'l', norb, norb, 1.d0, ovrlp2(1,1), &
!!!!!!             norb, lagmat(1,1), norb, &
!!!!!!             0.d0, ovrlp_minus_one_lagmat(1,1), norb)
!!!!!!        ! Transpose lagmat
!!!!!!        do iorb=1,norb
!!!!!!            do jorb=iorb+1,norb
!!!!!!                tt=lagmat(jorb,iorb)
!!!!!!                lagmat(jorb,iorb)=lagmat(iorb,jorb)
!!!!!!                lagmat(iorb,jorb)=tt
!!!!!!            end do
!!!!!!        end do
!!!!!!        call dsymm_parallel(iproc, nproc, blocksize_pdgemm, comm, 'l', 'l', norb, norb, 1.d0, ovrlp2(1,1), norb, &
!!!!!!            lagmat(1,1), norb, 0.d0, ovrlp_minus_one_lagmat_trans(1,1), norb)
!!!!!!    end if
!!!!!!else if(correctionOrthoconstraint==1) then correctionIf
!!!!!!    do iorb=1,norb
!!!!!!        do jorb=1,norb
!!!!!!            ovrlp_minus_one_lagmat(jorb,iorb)=lagmat(jorb,iorb)
!!!!!!            ovrlp_minus_one_lagmat_trans(jorb,iorb)=lagmat(iorb,jorb)
!!!!!!        end do
!!!!!!    end do
!!!!!!end if correctionIf
!!!!!!
!!!!!!
!!!!!!!!ilrold=-1
!!!!!!!!ijorb=0
!!!!!!!!do iorb=1,norbp
!!!!!!!!    iiorb=isorb+iorb
!!!!!!!!    ilr=onWhichAtom(iiorb)
!!!!!!!!    if(ilr==ilrold) then
!!!!!!!!        ! Set back the index of lphiovrlp, since we again need the same orbitals.
!!!!!!!!        !ijorb=ijorb-comom%noverlap(ilr)
!!!!!!!!        ijorb=ijorb-comom%noverlap(iiorb)
!!!!!!!!    end if
!!!!!!!!    ncount=mlr(ilr)%norbinlr
!!!!!!!!    !do jorb=1,comom%noverlap(ilr)
!!!!!!!!    do jorb=1,comom%noverlap(iiorb)
!!!!!!!!        ijorb=ijorb+1
!!!!!!!!        !jjorb=comom%overlaps(jorb,ilr)
!!!!!!!!        jjorb=comom%overlaps(jorb,iiorb)
!!!!!!!!        call daxpy(ncount, -.5d0*ovrlp_minus_one_lagmat(jjorb,iiorb), vecOvrlp(1,ijorb), 1, grad(1,iorb), 1)
!!!!!!!!        call daxpy(ncount, -.5d0*ovrlp_minus_one_lagmat_trans(jjorb,iiorb), vecOvrlp(1,ijorb), 1, grad(1,iorb), 1)
!!!!!!!!    end do
!!!!!!!!    ilrold=ilr
!!!!!!!!end do
!!!!!!
!!!!!!
!!!!!!!!iall=-product(shape(ovrlp_minus_one_lagmat))*kind(ovrlp_minus_one_lagmat)
!!!!!!!!deallocate(ovrlp_minus_one_lagmat, stat=istat)
!!!!!!!!call memocc(istat, iall, 'ovrlp_minus_one_lagmat', subname)
!!!!!!!!iall=-product(shape(ovrlp_minus_one_lagmat_trans))*kind(ovrlp_minus_one_lagmat_trans)
!!!!!!!!deallocate(ovrlp_minus_one_lagmat_trans, stat=istat)
!!!!!!!!call memocc(istat, iall, 'ovrlp_minus_one_lagmat_trans', subname)
!!!!!!iall=-product(shape(ovrlp2))*kind(ovrlp2)
!!!!!!deallocate(ovrlp2, stat=istat)
!!!!!!call memocc(istat, iall, 'ovrlp2', subname)
!!!!!!
!!!!!!
!!!!!!end subroutine applyOrthoconstraintVectors




!!!subroutine buildLinearCombinations(iproc, nproc, lzdig, lzd, orbsig, orbs, input, coeff, lchi, locregShape, &
!!!           tag, comonig, opig, madig, lphi)
!!!use module_base
!!!use module_types
!!!use module_interfaces, exceptThisOne => buildLinearCombinations
!!!implicit none
!!!
!!!! Calling arguments
!!!integer,intent(in):: iproc, nproc
!!!type(local_zone_descriptors),intent(in):: lzdig, lzd
!!!type(orbitals_data),intent(in):: orbsig, orbs
!!!type(input_variables),intent(in):: input
!!!real(8),dimension(orbsig%norb,orbs%norb),intent(in):: coeff
!!!real(8),dimension(orbsig%npsidim_orbs),intent(in):: lchi
!!!character(len=1),intent(in):: locregShape
!!!integer,intent(inout):: tag
!!!type(p2pComms):: comonig
!!!type(overlapParameters):: opig
!!!type(matrixDescriptors):: madig
!!!real(8),dimension(orbs%npsidim_orbs),intent(out):: lphi
!!!
!!!! Local variables
!!!integer:: istat, iall, ist, jst, ilr, ilrold, iorb, iiorb, ncount, jorb, jjorb
!!!!type(overlapParameters):: op
!!!!type(p2pCommsOrthonormality):: comon
!!!real(8),dimension(:),allocatable:: lchiovrlp
!!!character(len=*),parameter:: subname='buildLinearCombinations'
!!!!type(matrixDescriptors):: mad !just for calling collectnew, not really needed
!!!real(8),dimension(:,:),allocatable:: ttmat
!!!real(8):: tt1, tt2, tt3
!!!
!!!!tag=10000
!!!!call initCommsOrtho(iproc, nproc, lzdig, orbsig, orbsig%inWhichLocreg, input, locregShape, op, comon, tag)
!!!allocate(lchiovrlp(opig%ndim_lphiovrlp), stat=istat)
!!!call memocc(istat, lchiovrlp, 'lchiovrlp',subname)
!!!
!!!call allocateCommuncationBuffersOrtho(comonig, subname)
!!!!call extractOrbital2(iproc,nproc,orbsig,orbsig%npsidim,orbsig%inWhichLocreg,lzdig,op,lchi,comon)
!!!call extractOrbital3(iproc,nproc,orbsig,orbsig,orbsig%npsidim_orbs,orbsig%inWhichLocreg,&
!!!     lzdig,lzdig,opig,opig,lchi,comonig%nsendBuf,comonig%sendBuf)
!!!!call postCommsOverlap(iproc, nproc, comon)
!!!call postCommsOverlapNew(iproc, nproc, orbsig, opig, lzdig, lchi, comonig, tt1, tt2)
!!!!call gatherOrbitals2(iproc, nproc, comon)
!!!!!allocate(ttmat(orbsig%norb,orbsig%norb))
!!!call collectnew(iproc, nproc, comonig, madig, opig, orbsig, lzdig, comonig%nsendbuf, &
!!!     comonig%sendbuf, comonig%nrecvbuf, comonig%recvbuf, tt1, tt2, tt3)
!!!!!deallocate(ttmat)
!!!call expandOrbital2(iproc, nproc, orbsig, input, orbsig%inWhichLocreg, lzdig, opig, comonig, lchiovrlp)
!!!call deallocateCommuncationBuffersOrtho(comonig, subname)
!!!
!!!
!!!
!!!lphi=0.d0
!!!
!!!ist=1
!!!jst=1
!!!ilrold=-1
!!!do iorb=1,orbs%norbp
!!!    iiorb=orbs%isorb+iorb
!!!    ilr=orbs%inWhichLocreg(iiorb)
!!!    if(ilr==ilrold) then
!!!        ! Set back the index of lphiovrlp, since we again need the same orbitals.
!!!        jst=jst-opig%noverlaps(iiorb-1)*ncount
!!!    end if
!!!    !write(*,'(a,6i13)') 'iproc, iorb, iiorb, op%noverlaps(iiorb), ilr, jst', iproc, iorb, iiorb, op%noverlaps(iiorb), ilr, jst
!!!    ncount=lzd%llr(ilr)%wfd%nvctr_c+7*lzd%llr(ilr)%wfd%nvctr_f
!!!    do jorb=1,opig%noverlaps(iiorb)
!!!        jjorb=opig%overlaps(jorb,iiorb)
!!!        !call daxpy(ncount, ovrlp(jjorb,iiorb), lphiovrlp(jst), 1, lphi(ist), 1)
!!!        call daxpy(ncount, coeff(jjorb,iiorb), lchiovrlp(jst), 1, lphi(ist), 1)
!!!        jst=jst+ncount
!!!    end do
!!!
!!!    ist=ist+ncount
!!!    ilrold=ilr
!!!
!!!end do
!!!
!!!!!if(ist/=orbs%npsidim+1) then
!!!!!    write(*,'(a,i0,a,2(2x,i0))') 'ERROR on process ',iproc,': ist/=orbs%npsidim+1',ist,orbs%npsidim+1
!!!!!    stop
!!!!!end if
!!!if(ist>orbs%npsidim_orbs+1) then
!!!    write(*,'(a,i0,a,2(2x,i0))') 'ERROR on process ',iproc,': ist/=orbs%npsidim_orbs+1',ist,orbs%npsidim_orbs+1
!!!    stop
!!!end if
!!!
!!!
!!!
!!!!call deallocate_overlapParameters(op, subname)
!!!!call deallocate_p2pCommsOrthonormality(comon, subname)
!!!
!!!
!!!iall=-product(shape(lchiovrlp))*kind(lchiovrlp)
!!!deallocate(lchiovrlp, stat=istat)
!!!call memocc(istat, iall, 'lchiovrlp', subname)
!!!
!!!
!!!
!!!end subroutine buildLinearCombinations

!!!!subroutine buildLinearCombinationsVariable(iproc, nproc, lzdig, lzd, orbsig, orbs, input, coeff, lchi, tag, lphi)
!!!!use module_base
!!!!use module_types
!!!!use module_interfaces, exceptThisOne => buildLinearCombinationsVariable
!!!!implicit none
!!!!
!!!!! Calling arguments
!!!!integer,intent(in):: iproc, nproc
!!!!type(local_zone_descriptors),intent(in):: lzdig, lzd
!!!!type(orbitals_data),intent(in):: orbsig, orbs
!!!!type(input_variables),intent(in):: input
!!!!real(8),dimension(orbsig%norb,orbs%norb),intent(in):: coeff
!!!!real(8),dimension(orbsig%npsidim_orbs),intent(in):: lchi
!!!!integer,intent(inout):: tag
!!!!real(8),dimension(orbs%npsidim_orbs),intent(out):: lphi
!!!!
!!!!! Local variables
!!!!integer:: istat, iall, ist, jst, ilr, ilrold, iorb, iiorb, ncount, jorb, jjorb, ii
!!!!type(overlapParameters):: op
!!!!type(p2pComms):: comon
!!!!real(8),dimension(:),allocatable:: lchiovrlp
!!!!character(len=*),parameter:: subname='buildLinearCombinationsVariable'
!!!!
!!!!!tag=10000
!!!!call initCommsOrthoVariable(iproc, nproc, lzdig, orbs, orbsig, orbsig%inWhichLocreg, input, op, comon, tag)
!!!!allocate(lchiovrlp(op%ndim_lphiovrlp), stat=istat)
!!!!call memocc(istat, lchiovrlp, 'lchiovrlp',subname)
!!!!
!!!!call allocateCommuncationBuffersOrtho(comon, subname)
!!!!call extractOrbital2Variable(iproc, nproc, orbs, orbsig, orbsig%npsidim_orbs, lzdig, op, lchi, comon)
!!!!call postCommsOverlap(iproc, nproc, comon)
!!!!call gatherOrbitals2(iproc, nproc, comon)
!!!!!call mpi_barrier(mpi_comm_world, ist)
!!!!!stop
!!!!call expandOrbital2(iproc, nproc, orbs, input, orbs%inWhichLocreg, lzdig, op, comon, lchiovrlp)
!!!!
!!!!call deallocateCommuncationBuffersOrtho(comon, subname)
!!!!
!!!!
!!!!
!!!!lphi=0.d0
!!!!
!!!!ist=1
!!!!jst=1
!!!!ilrold=-1
!!!!do iorb=1,orbs%norbp
!!!!    iiorb=orbs%isorb+iorb
!!!!    ilr=orbs%inWhichLocreg(iiorb)
!!!!    if(ilr==ilrold) then
!!!!        ! Set back the index of lphiovrlp, since we again need the same orbitals.
!!!!        jst=jst-op%noverlaps(iiorb-1)*ncount
!!!!    end if
!!!!    ncount=lzd%llr(ilr)%wfd%nvctr_c+7*lzd%llr(ilr)%wfd%nvctr_f
!!!!    do jorb=1,op%noverlaps(iiorb)
!!!!        jjorb=op%overlaps(jorb,iiorb)
!!!!        !call daxpy(ncount, ovrlp(jjorb,iiorb), lphiovrlp(jst), 1, lphi(ist), 1)
!!!!        call daxpy(ncount, coeff(jjorb,iiorb), lchiovrlp(jst), 1, lphi(ist), 1)
!!!!        jst=jst+ncount
!!!!    end do
!!!!
!!!!    ist=ist+ncount
!!!!    ilrold=ilr
!!!!
!!!!end do
!!!!
!!!!if(ist/=orbs%npsidim_orbs+1) then
!!!!    write(*,'(a,i0,a,2(2x,i0))') 'ERROR on process ',iproc,&
!!!!         ': ist/=orbsig%npsidim+1',ist,orbsig%npsidim_orbs+1
!!!!    stop
!!!!end if
!!!!
!!!!
!!!!
!!!!call deallocate_overlapParameters(op, subname)
!!!!call deallocate_p2pComms(comon, subname)
!!!!
!!!!
!!!!iall=-product(shape(lchiovrlp))*kind(lchiovrlp)
!!!!deallocate(lchiovrlp, stat=istat)
!!!!call memocc(istat, iall, 'lchiovrlp', subname)
!!!!
!!!!
!!!!
!!!!end subroutine buildLinearCombinationsVariable


!!!subroutine initializeInguessParameters(iproc, nproc, orbs, orbsig, newComm, ip)
!!!use module_base
!!!use module_types
!!!implicit none
!!!
!!!! Calling arguments
!!!integer,intent(in):: iproc, nproc
!!!type(orbitals_data),intent(in):: orbs, orbsig
!!!integer,intent(in):: newComm
!!!type(inguessParameters),intent(inout):: ip
!!!
!!!! Local variables
!!!integer:: ii, kk, jproc, istat, ierr, norbTarget, iorb, iiorb
!!!real(8):: tt
!!!character(len=*),parameter:: subname='initializeInguessParameters'
!!!
!!!
!!!  !!ip%norb=orbs%norb
!!!  !!ip%norbtot=orbsig%norb
!!!
!!!  ! In order to symplify the transposing/untransposing, the orbitals are padded with zeros such that 
!!!  ! they can be distributed evenly over all processes when being transposed. The new length of the 
!!!  ! orbitals after this padding is then given by ip%norbtotPad.
!!!  !!ip%norbtotPad=orbsig%norb
!!!  !!do
!!!  !!    if(mod(ip%norbtotPad, nproc)==0) exit
!!!  !!    ip%norbtotPad=ip%norbtotPad+1
!!!  !!end do
!!!
!!!
!!!
!!!  !!!! Calculate the number of elements that each process has when the vectors are transposed.
!!!  !!!! nvctrp is the total number, nvctrp_nz is the nonzero numbers.
!!!  !!!allocate(ip%nvctrp_nz(0:nproc-1), stat=istat)
!!!  !!!call memocc(istat, ip%nvctrp_nz, 'ip%nvctrp_nz', subname)
!!!  !!!tt=ip%norbtot/dble(nproc)
!!!  !!!ii=floor(tt)
!!!  !!!! ii is now the number of elements that every process has. Distribute the remaining ones.
!!!  !!!ip%nvctrp_nz=ii
!!!  !!!kk=ip%norbtot-nproc*ii
!!!  !!!ip%nvctrp_nz(0:kk-1)=ii+1
!!!  !!!! Check wheter this distribution is correct
!!!  !!!ii=0
!!!  !!!do jproc=0,nproc-1
!!!  !!!   ii=ii+ip%nvctrp_nz(jproc)
!!!  !!!end do
!!!  !!!if(ii/=ip%norbtot) then
!!!  !!!   if(iproc==0) write(*,'(3x,a)') 'ERROR: wrong partition of ip%norbtot!'
!!!  !!!   call mpi_barrier(newComm, ierr)
!!!  !!!   stop
!!!  !!!end if
!!!
!!!  ! With the padded zeros, the elements can be distributed evenly.
!!!  !!!ip%nvctrp=ip%norbtotPad/nproc
!!!
!!!  ! Define the values for the mpi_alltoallv.
!!!  ! sendcounts: number of elements that a given  process sends to another process.
!!!  ! senddispls: offset of the starting index on a given process for the send operation to another process.
!!!  !!!allocate(ip%sendcounts(0:nproc-1), stat=istat)
!!!  !!!call memocc(istat, ip%sendcounts, 'ip%sendcounts', subname)
!!!  !!!allocate(ip%senddispls(0:nproc-1), stat=istat)
!!!  !!!call memocc(istat, ip%senddispls, 'ip%senddispls', subname)
!!!  !!!ii=0
!!!  !!!do jproc=0,nproc-1
!!!  !!!    ip%sendcounts(jproc)=ip%nvctrp*orbs%norb_par(iproc,0)
!!!  !!!    ip%senddispls(jproc)=ii
!!!  !!!    ii=ii+ip%sendcounts(jproc)
!!!  !!!end do
!!!  !!!! recvcounts: number of elements that a given process receives from another process.
!!!  !!!! recvdispls: offset of the starting index on a given process for the receive operation from another process.
!!!  !!!allocate(ip%recvcounts(0:nproc-1), stat=istat)
!!!  !!!call memocc(istat, ip%recvcounts, 'ip%recvcounts', subname)
!!!  !!!allocate(ip%recvdispls(0:nproc-1), stat=istat)
!!!  !!!call memocc(istat, ip%recvdispls, 'ip%recvdispls', subname)
!!!  !!!ii=0
!!!  !!!do jproc=0,nproc-1
!!!  !!!    ip%recvcounts(jproc)=ip%nvctrp*orbs%norb_par(jproc,0)
!!!  !!!    ip%recvdispls(jproc)=ii
!!!  !!!    ii=ii+ip%recvcounts(jproc)
!!!  !!!end do
!!!
!!!  !!!! Determine the size of the work array needed for the transposition.
!!!  !!!ip%sizeWork=max(ip%norbtotPad*orbs%norb_par(iproc,0),sum(ip%recvcounts(:)))
!!!
!!!end subroutine initializeInguessParameters



