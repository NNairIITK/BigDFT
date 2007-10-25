subroutine orthoconstraint_p(iproc,nproc,norb,occup,nvctrp,psit,hpsit,scprsum)
  !Effect of orthogonality constraints on gradient 
  implicit real(kind=8) (a-h,o-z)
  logical, parameter :: parallel=.true.
  dimension psit(nvctrp,norb),hpsit(nvctrp,norb),occup(norb)
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


subroutine orthoconstraint(norb,occup,nvctrp,psi,hpsi,scprsum)
  !Effect of orthogonality constraints on gradient 
  implicit real(kind=8) (a-h,o-z)
  logical, parameter :: parallel=.false.
  dimension psi(nvctrp,norb),hpsi(nvctrp,norb),occup(norb)
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


subroutine orthon_p(iproc,nproc,norb,nvctrp,psit)
  ! Gram-Schmidt orthogonalisation
  implicit real(kind=8) (a-h,o-z)
  logical, parameter :: parallel=.true.
  dimension psit(nvctrp,norb)
  real(kind=8), allocatable :: ovrlp(:,:,:)
  include 'mpif.h'

  call timing(iproc,'GramS_comput  ','ON')

  if (norb.eq.1) then

!!$         if (iproc .eq. 0) then
!!$            tt=dnrm2(nvctrp,psi,1)
!!$            tt=1.d0/tt
!!$            call dscal(nvctrp,tt,psi,1)
!!$         end if

     stop 'more than one orbital needed for a parallel run'

  else

     allocate(ovrlp(norb,norb,2),stat=i_stat)
     call memocc(i_stat,product(shape(ovrlp))*kind(ovrlp),'ovrlp','orthon_p')

     ! Upper triangle of overlap matrix using BLAS
     !     ovrlp(iorb,jorb)=psit(k,iorb)*psit(k,jorb) ; upper triangle
     call DSYRK('L','T',norb,nvctrp,1.d0,psit,nvctrp,0.d0,ovrlp(1,1,2),norb)

     call timing(iproc,'GramS_comput  ','OF')
     call timing(iproc,'GramS_commun  ','ON')
     call MPI_ALLREDUCE (ovrlp(1,1,2),ovrlp(1,1,1),norb**2,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
     call timing(iproc,'GramS_commun  ','OF')
     call timing(iproc,'GramS_comput  ','ON')

     !  write(*,*) 'parallel ovrlp'
     !  do i=1,norb
     !  write(*,'(10(1x,e10.3))') (ovrlp(i,j,1),j=1,norb)
     !  enddo


     ! Cholesky factorization
     call dpotrf( 'L', norb, ovrlp, norb, info )
     if (info.ne.0) write(6,*) 'info Cholesky factorization', info

     ! calculate L^{-1}
     call DTRTRI( 'L', 'N', norb, ovrlp, norb, info )
     if (info.ne.0) write(6,*) 'info L^-1', info

     ! new vectors   
     call DTRMM ('R', 'L', 'T', 'N', nvctrp, norb, 1.d0, ovrlp, norb, psit, nvctrp)

     i_all=-product(shape(ovrlp))*kind(ovrlp)
     deallocate(ovrlp,stat=i_stat)
     call memocc(i_stat,i_all,'ovrlp','orthon_p')

  end if

  call timing(iproc,'GramS_comput  ','OF')

END SUBROUTINE orthon_p


subroutine orthon(norb,nvctrp,psi)
  ! Gram-Schmidt orthogonalisation
  implicit real(kind=8) (a-h,o-z)
  logical, parameter :: parallel=.false.
  dimension psi(nvctrp,norb)
  real(kind=8), allocatable :: ovrlp(:,:)

  call timing(iproc,'GramS_comput  ','ON')

  if (norb.eq.1) then
     tt=dnrm2(nvctrp,psi,1)
     tt=1.d0/tt
     call dscal(nvctrp,tt,psi,1)

  else

     allocate(ovrlp(norb,norb),stat=i_stat)
     call memocc(i_stat,product(shape(ovrlp))*kind(ovrlp),'ovrlp','orthon')

     ! Overlap matrix using BLAS
     !     ovrlp(iorb,jorb)=psi(k,iorb)*psi(k,jorb) ; upper triangle
     call DSYRK('L','T',norb,nvctrp,1.d0,psi,nvctrp,0.d0,ovrlp,norb)

     !  write(*,*) 'ovrlp'
     !  do i=1,norb
     !  write(*,'(10(1x,e10.3))') (ovrlp(i,j),j=1,norb)
     !  enddo

     ! Cholesky factorization
     call dpotrf( 'L', norb, ovrlp, norb, info )
     if (info.ne.0) write(6,*) 'info Cholesky factorization', info

     ! calculate L^{-1}
     call DTRTRI( 'L', 'N', norb, ovrlp, norb, info )
     if (info.ne.0) write(6,*) 'info L^-1', info

     ! new vectors   
     call DTRMM ('R', 'L', 'T', 'N', nvctrp, norb, 1.d0, ovrlp, norb, psi, nvctrp)

     i_all=-product(shape(ovrlp))*kind(ovrlp)
     deallocate(ovrlp,stat=i_stat)
     call memocc(i_stat,i_all,'ovrlp','orthon')

  endif

  call timing(iproc,'GramS_comput  ','OF')

end subroutine orthon


subroutine loewe_p(iproc,nproc,norb,ndim,nvctrp,psit)
  ! loewdin orthogonalisation
  implicit real(kind=8) (a-h,o-z)
  logical, parameter :: parallel=.true.
  dimension psit(nvctrp,ndim)
  real(kind=8), allocatable :: ovrlp(:,:,:),evall(:),psitt(:,:)
  include 'mpif.h'

  if (norb.eq.1) stop 'more than one orbital needed for a parallel run'

  allocate(ovrlp(norb,norb,3),stat=i_stat)
  call memocc(i_stat,product(shape(ovrlp))*kind(ovrlp),'ovrlp','loewe_p')
  allocate(evall(norb),stat=i_stat)
  call memocc(i_stat,product(shape(evall))*kind(evall),'evall','loewe_p')

  ! Upper triangle of overlap matrix using BLAS
  !     ovrlp(iorb,jorb)=psit(k,iorb)*psit(k,jorb) ; upper triangle
  call DSYRK('U','T',norb,nvctrp,1.d0,psit,nvctrp,0.d0,ovrlp(1,1,2),norb)

  ! Full overlap matrix using  BLAS
  !     ovrlap(jorb,iorb,2)=+psit(k,jorb)*psit(k,iorb)
  !      call DGEMM('T','N',norb,norb,nvctrp,1.d0,psit,nvctrp,psit,nvctrp,0.d0,ovrlp(1,1,2),norb)

  call MPI_ALLREDUCE (ovrlp(1,1,2),ovrlp(1,1,1),norb**2,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)

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
  call DGEMM('N','T',norb,norb,norb,1.d0,ovrlp(1,1,1),norb,ovrlp(1,1,2),norb,0.d0,ovrlp(1,1,3),norb)

  allocate(psitt(nvctrp,ndim),stat=i_stat)
  call memocc(i_stat,product(shape(psitt))*kind(psitt),'psitt','loewe_p')
  ! new eigenvectors
  !   psitt(i,iorb)=psit(i,jorb)*ovrlp(jorb,iorb,3)
  call DGEMM('N','N',nvctrp,norb,norb,1.d0,psit,nvctrp,ovrlp(1,1,3),norb,0.d0,psitt,nvctrp)
  call DCOPY(nvctrp*ndim,psitt,1,psit,1)
  i_all=-product(shape(psitt))*kind(psitt)
  deallocate(psitt,stat=i_stat)
  call memocc(i_stat,i_all,'psitt','loewe_p')

  i_all=-product(shape(ovrlp))*kind(ovrlp)
  deallocate(ovrlp,stat=i_stat)
  call memocc(i_stat,i_all,'ovrlp','loewe_p')
  i_all=-product(shape(evall))*kind(evall)
  deallocate(evall,stat=i_stat)
  call memocc(i_stat,i_all,'evall','loewe_p')

END SUBROUTINE loewe_p


subroutine loewe(norb,nvctrp,psi)
  ! loewdin orthogonalisation
  implicit real(kind=8) (a-h,o-z)
  dimension psi(nvctrp,norb)
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

     allocate(ovrlp(norb,norb,3),stat=i_stat)
     call memocc(i_stat,product(shape(ovrlp))*kind(ovrlp),'ovrlp','loewe')
     allocate(evall(norb),stat=i_stat)
     call memocc(i_stat,product(shape(evall))*kind(evall),'evall','loewe')

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
     allocate(tpsi(nvctrp,norb),stat=i_stat)
     call memocc(i_stat,product(shape(tpsi))*kind(tpsi),'tpsi','loewe')
     !   tpsi(i,iorb)=psi(i,jorb)*ovrlp(jorb,iorb,3)
     call DGEMM('N','N',nvctrp,norb,norb,1.d0,psi,nvctrp,ovrlp(1,1,3),norb,0.d0,tpsi,nvctrp)
     call DCOPY(nvctrp*norb,tpsi,1,psi,1)
     i_all=-product(shape(tpsi))*kind(tpsi)
     deallocate(tpsi,stat=i_stat)
     call memocc(i_stat,i_all,'tpsi','loewe')

     i_all=-product(shape(ovrlp))*kind(ovrlp)
     deallocate(ovrlp,stat=i_stat)
     call memocc(i_stat,i_all,'ovrlp','loewe')
     i_all=-product(shape(evall))*kind(evall)
     deallocate(evall,stat=i_stat)
     call memocc(i_stat,i_all,'evall','loewe')

  endif

END SUBROUTINE loewe


subroutine checkortho_p(iproc,nproc,norb,nvctrp,psit)
  implicit real(kind=8) (a-h,o-z)
  dimension psit(nvctrp,norb)
  real(kind=8), allocatable :: ovrlp(:,:,:)
  include 'mpif.h'

  allocate(ovrlp(norb,norb,2),stat=i_stat)
  call memocc(i_stat,product(shape(ovrlp))*kind(ovrlp),'ovrlp','checkortho_p')

  do iorb=1,norb
     do jorb=1,norb
        ovrlp(iorb,jorb,2)=ddot(nvctrp,psit(1,iorb),1,psit(1,jorb),1)
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
  call memocc(i_stat,i_all,'ovrlp','checkortho_p')

END SUBROUTINE checkortho_p


subroutine checkortho(norb,nvctrp,psi)
  implicit real(kind=8) (a-h,o-z)
  dimension psi(nvctrp,norb)
  real(kind=8), allocatable :: ovrlp(:,:,:)

  allocate(ovrlp(norb,norb,1),stat=i_stat)
  call memocc(i_stat,product(shape(ovrlp))*kind(ovrlp),'ovrlp','checkortho')

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
  call memocc(i_stat,i_all,'ovrlp','checkortho')


END SUBROUTINE checkortho


subroutine KStrans_p(iproc,nproc,norb,ndim,nvctrp,occup,  & 
     hpsit,psit,evsum,eval)
  ! at the start each processor has all the Psi's but only its part of the HPsi's
  ! at the end each processor has only its part of the Psi's
  implicit real(kind=8) (a-h,o-z)
  dimension occup(norb),eval(norb)
  dimension psit(nvctrp,ndim),hpsit(nvctrp,ndim)
  ! arrays for KS orbitals
  allocatable :: hamks(:,:,:),work_lp(:),psitt(:,:)
  include 'mpif.h'

  ! set up Hamiltonian matrix
  allocate(hamks(norb,norb,2),stat=i_stat)
  call memocc(i_stat,product(shape(hamks))*kind(hamks),'hamks','kstrans_p')
  do jorb=1,norb
     do iorb=1,norb
        hamks(iorb,jorb,2)=0.d0
     enddo
  enddo
  do iorb=1,norb
     do jorb=1,norb
        scpr=ddot(nvctrp,psit(1,jorb),1,hpsit(1,iorb),1)
        hamks(iorb,jorb,2)=scpr
     enddo
  enddo

  call MPI_ALLREDUCE(hamks(1,1,2),hamks(1,1,1),norb**2,MPI_DOUBLE_PRECISION,&
       MPI_SUM,MPI_COMM_WORLD,ierr)

  !        write(*,*) 'KS Hamiltonian',iproc
  !        do iorb=1,norb
  !        write(*,'(10(1x,e10.3))') (hamks(iorb,jorb,1),jorb=1,norb)
  !        enddo

  n_lp=max(4*norb,1000)
  allocate(work_lp(n_lp),stat=i_stat)
  call memocc(i_stat,product(shape(work_lp))*kind(work_lp),'work_lp','kstrans_p')
  call  DSYEV('V','U',norb,hamks,norb,eval, work_lp, n_lp, info )
  evsum=0.d0
  do iorb=1,norb
     evsum=evsum+eval(iorb)*occup(iorb)
     !if (iproc.eq.0) write(*,'(1x,a,i0,a,1x,1pe21.14)') 'eval(',iorb,')=',eval(iorb)
  enddo
  i_all=-product(shape(work_lp))*kind(work_lp)
  deallocate(work_lp,stat=i_stat)
  call memocc(i_stat,i_all,'work_lp','kstrans_p')
  if (info.ne.0) write(*,*) 'DSYEV ERROR',info

  allocate(psitt(nvctrp,ndim),stat=i_stat)
  call memocc(i_stat,product(shape(psitt))*kind(psitt),'psitt','kstrans_p')
  ! Transform to KS orbitals
  do iorb=1,norb
     call razero(nvctrp,psitt(1,iorb))
     do jorb=1,norb
        alpha=hamks(jorb,iorb,1)
        call daxpy(nvctrp,alpha,psit(1,jorb),1,psitt(1,iorb),1)
     enddo
  enddo
  i_all=-product(shape(hamks))*kind(hamks)
  deallocate(hamks,stat=i_stat)
  call memocc(i_stat,i_all,'hamks','kstrans_p')

  call DCOPY(nvctrp*ndim,psitt,1,psit,1)
  i_all=-product(shape(psitt))*kind(psitt)
  deallocate(psitt,stat=i_stat)
  call memocc(i_stat,i_all,'psitt','kstrans_p')

END SUBROUTINE KStrans_p


subroutine KStrans(norb,nvctrp,occup,hpsi,psi,evsum,eval)
  ! at the start each processor has all the Psi's but only its part of the HPsi's
  ! at the end each processor has only its part of the Psi's
  implicit real(kind=8) (a-h,o-z)
  dimension occup(norb),eval(norb)
  dimension psi(nvctrp,norb),hpsi(nvctrp,norb)
  ! arrays for KS orbitals
  allocatable :: hamks(:,:,:),work_lp(:),psitt(:,:)

  ! set up Hamiltonian matrix
  allocate(hamks(norb,norb,2),stat=i_stat)
  call memocc(i_stat,product(shape(hamks))*kind(hamks),'hamks','kstrans')
  do jorb=1,norb
     do iorb=1,norb
        hamks(iorb,jorb,2)=0.d0
     enddo
  enddo
  do iorb=1,norb
     do jorb=1,norb
        scpr=ddot(nvctrp,psi(1,jorb),1,hpsi(1,iorb),1)
        hamks(iorb,jorb,1)=scpr
     enddo
  enddo

  !        write(*,*) 'KS Hamiltonian',0
  !        do iorb=1,norb
  !        write(*,'(10(1x,e10.3))') (hamks(iorb,jorb,1),jorb=1,norb)
  !        enddo

  n_lp=max(4*norb,1000)
  allocate(work_lp(n_lp),stat=i_stat)
  call memocc(i_stat,product(shape(work_lp))*kind(work_lp),'work_lp','kstrans')
  call  DSYEV('V','U',norb,hamks,norb,eval, work_lp, n_lp, info )
  evsum=0.d0
  do iorb=1,norb
     evsum=evsum+eval(iorb)*occup(iorb)
     !write(*,'(1x,a,i0,a,1x,1pe21.14)') 'eval(',iorb,')=',eval(iorb)
  enddo
  i_all=-product(shape(work_lp))*kind(work_lp)
  deallocate(work_lp,stat=i_stat)
  call memocc(i_stat,i_all,'work_lp','kstrans')
  if (info.ne.0) write(*,*) 'DSYEV ERROR',info

  allocate(psitt(nvctrp,norb),stat=i_stat)
  call memocc(i_stat,product(shape(psitt))*kind(psitt),'psitt','kstrans')
  ! Transform to KS orbitals
  do iorb=1,norb
     call razero(nvctrp,psitt(1,iorb))
     do jorb=1,norb
        alpha=hamks(jorb,iorb,1)
        call daxpy(nvctrp,alpha,psi(1,jorb),1,psitt(1,iorb),1)
     enddo
  enddo
  i_all=-product(shape(hamks))*kind(hamks)
  deallocate(hamks,stat=i_stat)
  call memocc(i_stat,i_all,'hamks','kstrans')

  call DCOPY(nvctrp*norb,psitt,1,psi,1)
  i_all=-product(shape(psitt))*kind(psitt)
  deallocate(psitt,stat=i_stat)
  call memocc(i_stat,i_all,'psitt','kstrans')

END SUBROUTINE KStrans
