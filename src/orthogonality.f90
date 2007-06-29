subroutine orthoconstraint_p(iproc,nproc,norb,norbp,occup,nvctrp,psit,hpsit,scprsum)
  !Effect of orthogonality constraints on gradient 
  implicit real*8 (a-h,o-z)
  logical, parameter :: parallel=.true.
  dimension psit(nvctrp,norbp*nproc),hpsit(nvctrp,norbp*nproc),occup(norb)
  allocatable :: alag(:,:,:)
  include 'mpif.h'

  call timing(iproc,'LagrM_comput  ','ON')

  allocate(alag(norb,norb,2),stat=i_stat)
  call memocc(i_stat,product(shape(alag))*kind(alag),'alag','orthoconstraint_p')
  !     alag(jorb,iorb,2)=+psit(k,jorb)*hpsit(k,iorb)
  call DGEMM('T','N',norb,norb,nvctrp,1.d0,psit,nvctrp,hpsit,nvctrp,0.d0,alag(1,1,2),norb)

  call timing(iproc,'LagrM_comput  ','OF')
  call timing(iproc,'LagrM_commun  ','ON')
  call MPI_ALLREDUCE(alag(1,1,2),alag(1,1,1),norb**2,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
  call timing(iproc,'LagrM_commun  ','OF')
  call timing(iproc,'LagrM_comput  ','ON')
  !        if (iproc.eq.0) then
  !        write(*,*) 'ALAG',iproc
  !        do iorb=1,norb
  !        write(*,'(10(1x,1pe10.3))') (alag(iorb,jorb,1),jorb=1,norb)
  !        enddo
  !        endif

  scprsum=0.d0
  do iorb=1,norb
     scprsum=scprsum+occup(iorb)*alag(iorb,iorb,1)
  enddo

  ! hpsit(k,iorb)=-psit(k,jorb)*alag(jorb,iorb,1)
  call DGEMM('N','N',nvctrp,norb,norb,-1.d0,psit,nvctrp,alag,norb,1.d0,hpsit,nvctrp)
  i_all=-product(shape(alag))*kind(alag)
  deallocate(alag,stat=i_stat)
  call memocc(i_stat,i_all,'alag','orthoconstraint_p')

  call timing(iproc,'LagrM_comput  ','OF')

END SUBROUTINE orthoconstraint_p

subroutine orthoconstraint(norb,norbp,occup,nvctrp,psi,hpsi,scprsum)
  !Effect of orthogonality constraints on gradient 
  implicit real*8 (a-h,o-z)
  logical, parameter :: parallel=.false.
  dimension psi(nvctrp,norbp),hpsi(nvctrp,norbp),occup(norb)
  allocatable :: alag(:,:,:)

  call timing(iproc,'LagrM_comput  ','ON')

  allocate(alag(norb,norb,2),stat=i_stat)
  call memocc(i_stat,product(shape(alag))*kind(alag),'alag','orthoconstraint')

  !     alag(jorb,iorb,2)=+psi(k,jorb)*hpsi(k,iorb)
  call DGEMM('T','N',norb,norb,nvctrp,1.d0,psi,nvctrp,hpsi,nvctrp,0.d0,alag(1,1,1),norb)

  scprsum=0.d0
  do iorb=1,norb
     scprsum=scprsum+occup(iorb)*alag(iorb,iorb,1)
  enddo

  ! hpsit(k,iorb)=-psit(k,jorb)*alag(jorb,iorb,1)
  call DGEMM('N','N',nvctrp,norb,norb,-1.d0,psi,nvctrp,alag,norb,1.d0,hpsi,nvctrp)

  i_all=-product(shape(alag))*kind(alag)
  deallocate(alag,stat=i_stat)
  call memocc(i_stat,i_all,'alag','orthoconstraint')

  call timing(iproc,'LagrM_comput  ','OF')

END SUBROUTINE orthoconstraint
