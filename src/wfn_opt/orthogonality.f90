!orthogonality routine, for all the orbitals
!uses wavefunctions in their transposed form
subroutine orthogonalize(iproc,nproc,orbs,comms,wfd,psi)
  use module_base
  use module_types
  implicit none
  integer, intent(in) :: iproc,nproc
  type(orbitals_data), intent(in) :: orbs
  type(communications_arrays), intent(in) :: comms
  type(wavefunctions_descriptors), intent(in) :: wfd
  real(wp), dimension(sum(comms%nvctr_par(iproc,1:orbs%nkptsp))*orbs%nspinor*orbs%norb), intent(inout) :: psi
  !local variables
  character(len=*), parameter :: subname='orthogonalize'
  integer :: i_stat,i_all,i_err,info
  integer :: istart,ispin,nspin,ikpt,norb,norbs,ncomp,nvctrp,ispsi
  integer, dimension(:), allocatable :: ndimovrlp
  real(wp), dimension(:,:), allocatable :: ovrlp

  !separate the orthogonalisation procedure for up and down orbitals 
  !and for different k-points
  call timing(iproc,'GramS_comput  ','ON')

  !number of components of the overlap matrix for parallel case
  istart=2
  if (nproc == 1) istart=1

  !calculate the dimension of the overlap matrix for each k-point
  if (orbs%norbd > 0) then
     nspin=2
  else
     nspin=1
  end if

  !number of components for the overlap matrix in wp-kind real numbers

  allocate(ndimovrlp(nspin,0:orbs%nkpts+ndebug),stat=i_stat)
  call memocc(i_stat,ndimovrlp,'ndimovrlp',subname)

  ndimovrlp(1:2,0)=0

  do ikpt=1,orbs%nkpts
     !this part should be enhanced for real k-points
     norb=orbs%norbu
     if(orbs%nspinor==1) then
        norbs=norb
        ncomp=1 !useless
     else if (orbs%nspinor ==2) then
        norbs=2*norb
        ncomp=1
     else if (orbs%nspinor ==4) then
        norbs=2*norb
        ncomp=2
     end if
     ndimovrlp(1,ikpt)=ndimovrlp(nspin,ikpt-1)+norbs*norb
     if (orbs%norbd > 0) then
        norb=orbs%norbd
        if(orbs%nspinor==1) then
           norbs=norb
           ncomp=1 !useless
        else if (orbs%nspinor ==2) then
           norbs=2*norb
           ncomp=1
        else if (orbs%nspinor ==4) then
           norbs=2*norb
           ncomp=2
        end if
        ndimovrlp(2,ikpt)=ndimovrlp(1,ikpt)+norbs*norb
     end if
  end do

  allocate(ovrlp(ndimovrlp(nspin,orbs%nkpts),istart+ndebug),stat=i_stat)
  call memocc(i_stat,ovrlp,'ovrlp',subname)

  !put to zero all the k-points which are not needed
  call razero(ndimovrlp(nspin,orbs%nkpts)*istart,ovrlp)

  !differentiate between real and complex wavefunctions
  !Lower triangle of overlap matrix using BLAS
  !     ovrlp(iorb,jorb)=psit(k,iorb)*psit(k,jorb) ; lower triangle


  !do it for each of the k-points and separate also between up and down orbitals in the non-collinear case
  ispsi=1
  do ikptp=1,orbs%nkptsp
     ikpt=orbs%iskpts+ikptp

     do ispin=1,nspin
        nvctrp=comms%nvctr_par(iproc,ikptp)
        norb=orbs%norbu
        if (ispin==2) norb=orbs%norbd

        if(orbs%nspinor==1) then
           norbs=norb
           ncomp=1 !useless
        else if (orbs%nspinor ==2) then
           norbs=2*norb
           ncomp=1
        else if (orbs%nspinor ==4) then
           norbs=2*norb
           ncomp=2
        end if

        if(orbs%nspinor==1) then
           call syrk('L','T',norb,nvctrp,1.0_wp,psit(ispsi),max(1,nvctrp),&
                0.0_wp,ovrlp(ndimovrlp(ispin,ikpt-1)+1,istart),norb)
        else
           call herk('L','C',norb,ncomp*nvctrp,1.0_wp,psit(ispsi),max(1,ncomp*nvctrp),&
                0.0_wp,ovrlp(ndimovrlp(ispin,ikpt-1)+1,istart),norb)
        end if
        ispsi=ispsi+nvctrp*norb*orbs%nspinor
     end do
  end do

  if (nproc > 1) then
     call timing(iproc,'GramS_comput  ','OF')
     call timing(iproc,'GramS_commun  ','ON')
     call MPI_ALLREDUCE (ovrlp(1,2),ovrlp(1,1),ndimovrlp(nspin,orbs%nkpts),&
          mpidtypw,MPI_SUM,MPI_COMM_WORLD,ierr)
     call timing(iproc,'GramS_commun  ','OF')
     call timing(iproc,'GramS_comput  ','ON')
  end if

  !now each processors knows all the overlap matrices for each k-point
  !even if it does not handle it.
  !this is somehow redundant but it is the only way of reducing the number of communications
  !without defining group of processors

  !for each k-point now reorthogonalise wavefunctions
  ispsi=1
  do ikptp=1,orbs%nkptsp
     ikpt=orbs%iskpts+ikptp

     do ispin=1,nspin
        nvctrp=comms%nvctr_par(iproc,ikptp)
        norb=orbs%norbu
        if (ispin==2) norb=orbs%norbd

        if(orbs%nspinor==1) then
           norbs=norb
           ncomp=1 !useless
        else if (orbs%nspinor ==2) then
           norbs=2*norb
           ncomp=1
        else if (orbs%nspinor ==4) then
           norbs=2*norb
           ncomp=2
        end if

        !to be excluded if nvctrp==0
        if(orbs%nspinor==1 .and. nvctrp /= 0) then

           ! Cholesky factorization
           call potrf( 'L',norb,ovrlp(ndimovrlp(ispin,ikpt-1)+1,1),norb,info)
           if (info /= 0) then
              write(*,*) 'info Cholesky factorization',info
           end if

           ! calculate L^{-1}
           call trtri( 'L','N',norb,ovrlp(ndimovrlp(ispin,ikpt-1)+1,1),norb,info)
           if (info.ne.0) write(6,*) 'info L^-1',info

           ! new vectors   
           call trmm('R','L','T','N',nvctrp,norb,1.0_wp,ovrlp(ndimovrlp(ispin,ikpt-1)+1,1),&
                norb,psit(ispsi),max(1,nvctrp))

        else

           ! Cholesky factorization
           call c_potrf( 'L',norb,ovrlp(ndimovrlp(ispin,ikpt-1)+1,1),norb,info )
           if (info /= 0) then
              write(*,*) 'info Cholesky factorization',info
           end if

           ! calculate L^{-1}
           call c_trtri( 'L','N',norb,ovrlp(ndimovrlp(ispin,ikpt-1)+1,1),norb,info)
           if (info /= 0) write(6,*) 'info L^-1',info
           ! new vectors   !!check if third argument should be transpose or conjugate
           call c_trmm('R','L','C','N',ncomp*nvctrp,norb,(1.0_wp,0.0_wp),&
                ovrlp(ndimovrlp(ispin,ikpt-1)+1,1),norb,psit(ispsi),max(1,ncomp*nvctrp))

        end if
        ispsi=ispsi+nvctrp*norb*orbs%nspinor
     end do
  end do


  i_all=-product(shape(ovrlp))*kind(ovrlp)
  deallocate(ovrlp,stat=i_stat)
  call memocc(i_stat,i_all,'ovrlp',subname)

  i_all=-product(shape(ndimovrlp))*kind(ndimovrlp)
  deallocate(ndimovrlp,stat=i_stat)
  call memocc(i_stat,i_all,'ndimovrlp',subname)


  call timing(iproc,'GramS_comput  ','OF')

!!$
!!$call orthon_p(iproc,nproc,orbs%norbu,comms%nvctr_par(iproc),wfd%nvctr_c+7*wfd%nvctr_f,&
!!$    psi,orbs%nspinor) 
!!$if(orbs%norbd > 0) then
!!$  call orthon_p(iproc,nproc,orbs%norbd,comms%nvctr_par(iproc),wfd%nvctr_c+7*wfd%nvctr_f,&
!!$       psi(1+comms%nvctr_par(iproc)*orbs%norbu),orbs%nspinor) 
!!$end if

end subroutine orthogonalize


subroutine orthoconstraint_p(iproc,nproc,norb,occup,nvctrp,psit,hpsit,scprsum,nspinor)
  !Effect of orthogonality constraints on gradient 
  use module_base
  implicit none
  integer, intent(in) :: iproc,nproc,norb,nvctrp,nspinor
  real(gp), dimension(norb), intent(in) :: occup
  real(wp), dimension(nspinor*nvctrp,norb), intent(in) :: psit 
  real(dp), intent(out) :: scprsum
  real(wp), dimension(nspinor*nvctrp,norb), intent(out) :: hpsit
  !local variables
  character(len=*), parameter :: subname='orthoconstraint_p'
  integer :: i_stat,i_all,istart,iorb,jorb,ierr,norbs,i,j,ncomp
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
     call razero(norbs*norb*istart,alag)
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
          mpidtypd,MPI_SUM,MPI_COMM_WORLD,ierr)
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

end subroutine orthoconstraint_p

subroutine orthon_p(iproc,nproc,norb,nvctrp,nvctr_tot,psit,nspinor)
  ! Gram-Schmidt orthogonalisation
  use module_base
  implicit none
  integer, intent(in) :: iproc,nproc,norb,nvctrp,nvctr_tot,nspinor
  real(wp), dimension(nspinor*nvctrp,norb), intent(inout) :: psit
  !local variables
  character(len=*), parameter :: subname='orthon_p'
  integer :: info,i_all,i_stat,nvctr_eff,ierr,istart,i,j,norbs,iorb,jorb,ncomp
  real(wp) :: tt,ttLOC,ttr,tti
  real(wp), dimension(:,:,:), allocatable :: ovrlp

  call timing(iproc,'GramS_comput  ','ON')

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
        call MPI_ALLREDUCE(ttLOC,tt,1,mpidtypd,MPI_SUM,MPI_COMM_WORLD,ierr)
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

     call razero(norbs*norb*istart,ovrlp)

     ! Upper triangle of overlap matrix using BLAS
     !     ovrlp(iorb,jorb)=psit(k,iorb)*psit(k,jorb) ; upper triangle
     if(nspinor==1) then
        call syrk('L','T',norb,nvctrp,1.0_wp,psit(1,1),max(1,nvctrp),0.0_wp,ovrlp(1,1,istart),norb)
     else
!!$        ovrlp=0.0d0
!!$        do iorb=1,norb
!!$           do jorb=1,norb
!!$              ttr=ddot(nvctrp*nspinor,psit(1,iorb),1,psit(1,jorb),1)
!!$              tti=ddot(nvctrp,psit(1,iorb),1,psit(nvctrp+1,jorb),1)
!!$              tti=tti-ddot(nvctrp,psit(nvctrp+1,iorb),1,psit(1,jorb),1)
!!$              tti=tti-ddot(nvctrp,psit(3*nvctrp+1,iorb),1,psit(2*nvctrp+1,jorb),1)
!!$              tti=tti+ddot(nvctrp,psit(2*nvctrp+1,iorb),1,psit(3*nvctrp+1,jorb),1)
!!$              ovrlp(2*iorb-1,jorb,1)=ttr
!!$              ovrlp(2*iorb,jorb,1)=tti*0.0d0
!!$              print *,iorb,norb,ttr,tti
!!$           end do
!!$        end do
!!$        stop
        call herk('L','C',norb,ncomp*nvctrp,1.0_wp,psit(1,1),max(1,ncomp*nvctrp),&
             0.0_wp,ovrlp(1,1,istart),norb)
     end if

     if (nproc > 1) then
        call timing(iproc,'GramS_comput  ','OF')
        call timing(iproc,'GramS_commun  ','ON')
        call MPI_ALLREDUCE (ovrlp(1,1,2),ovrlp(1,1,1),norbs*norb,&
             MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
        call timing(iproc,'GramS_commun  ','OF')
        call timing(iproc,'GramS_comput  ','ON')
     end if
     
!!$     if (iproc==0) then
!!$        write(*,*) 'parallel ovrlp'
!!$        do i=1,norbs,min(2,nspinor)
!!$           write(*,'(10(1x,1pe10.3))') (ovrlp(i,j,1),j=1,norb)
!!$        enddo
!!$     end if

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
!!$        do i=1,norb
!!$           if(iproc==0) then
!!$              write(*,*) 'parallel ovrlp',i
!!$              write(*,'(10f10.3)') (ovrlp(j,i,1), j=1,norbs)
!!$           end if
!!$        end do
        call c_potrf( 'L',norb,ovrlp(1,1,1),norb,info )
        if (info /= 0) then
           write(*,*) 'info Cholesky factorization',info
        end if
        
        ! calculate L^{-1}
!!$         do i=1,norb
!!$           if(iproc==0) then
!!$              write(*,*) 'parallel ovrlp2',i
!!$              write(*,'(10f10.3)') (ovrlp(j,i,1), j=1,norbs)
!!$           end if
!!$        end do
       call c_trtri( 'L','N',norb,ovrlp(1,1,1),norb,info)
        if (info.ne.0) write(6,*) 'info L^-1',info
        
!!$        do i=1,norb
!!$           if(iproc==0) then
!!$              write(*,'(10f10.3)') (ovrlp(j,i,1), j=1,norbs)
!!$           end if
!!$        end do
       ! new vectors   !!check if third argument should be transpose or conjugate
        call c_trmm ('R','L','C','N',ncomp*nvctrp,norb,(1.0_wp,0.0_wp),&
             ovrlp(1,1,1),norb,psit(1,1),max(1,ncomp*nvctrp))

        !if(nproc==1) call psitransspi(nvctrp,norb,psit,.true.)


     end if

     i_all=-product(shape(ovrlp))*kind(ovrlp)
     deallocate(ovrlp,stat=i_stat)
     call memocc(i_stat,i_all,'ovrlp',subname)

  end if

  call timing(iproc,'GramS_comput  ','OF')

END SUBROUTINE orthon_p

!the loewe routines must be uniformised serial/parallel and nspinor should be added
subroutine loewe_p(iproc,nproc,norb,ndim,nvctrp,nvctr_tot,psit)
  ! loewdin orthogonalisation
  use module_base
  implicit real(kind=8) (a-h,o-z)
  logical, parameter :: parallel=.true.
  dimension psit(nvctrp,ndim)
  character(len=*), parameter :: subname='loewe_p'
  real(kind=8), allocatable :: ovrlp(:,:,:),evall(:),psitt(:,:)

  if (norb.eq.1) then

     nvctr_eff=min(nvctr_tot-iproc*nvctrp,nvctrp)

     if (nvctr_eff > 0) then
     !parallel treatment of a run with only one orbital
     tt=dnrm2(nvctr_eff,psit,1)     
     ttLOC=tt**2

     else
        ttLOC =0.d0
     end if
     
     call MPI_ALLREDUCE(ttLOC,tt,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)

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
          MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)

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
     call DCOPY(nvctrp*ndim,psitt,1,psit,1)
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


subroutine loewe(norb,nvctrp,psi)
  ! loewdin orthogonalisation
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
     call DCOPY(nvctrp*norb,tpsi,1,psi,1)
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


subroutine checkortho_p(iproc,nproc,norb,nvctrp,psit)
  use module_base
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


subroutine KStrans_p(iproc,nproc,norb,nvctrp,occup,  & 
     hpsit,psit,evsum,eval,nspinor)
  ! at the start each processor has all the Psi's but only its part of the HPsi's
  ! at the end each processor has only its part of the Psi's
  !implicit real(kind=8) (a-h,o-z)
  use module_base
  implicit none
  integer, intent(in) :: iproc,nproc,norb,nvctrp,nspinor
  real(wp), intent(out) :: evsum
  real(gp), dimension(norb), intent(in) :: occup
  real(wp), dimension(nvctrp*nspinor,norb), intent(in) :: hpsit
  real(wp), dimension(norb), intent(out) :: eval
  real(wp), dimension(nvctrp*nspinor,norb), intent(out) :: psit
  !local variables
  character(len=*), parameter :: subname='KStrans_p'
  integer :: i_all,i_stat,ierr,iorb,jorb,n_lp,istart,info,norbs,ncomp
  real(dp) :: scpr,alpha
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
          MPI_SUM,MPI_COMM_WORLD,ierr)
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
  if(nspinor==1) then
     do iorb=1,norb
        call razero(nvctrp,psitt(1,iorb))
        do jorb=1,norb
           alpha=hamks(jorb,iorb,1)
           call axpy(nvctrp,alpha,psit(1,jorb),1,psitt(1,iorb),1)
        enddo
     enddo
  else
     do iorb=1,norb
        call razero(nvctrp*nspinor,psitt(1,iorb))
        do jorb=1,norb
           call c_axpy(ncomp*nvctrp,hamks(2*jorb-1,iorb,1),psit(1,jorb),1,psitt(1,iorb),1)
        enddo
     enddo
  end if
  i_all=-product(shape(hamks))*kind(hamks)
  deallocate(hamks,stat=i_stat)
  call memocc(i_stat,i_all,'hamks',subname)

  call DCOPY(nvctrp*norb*nspinor,psitt,1,psit,1)
  i_all=-product(shape(psitt))*kind(psitt)
  deallocate(psitt,stat=i_stat)
  call memocc(i_stat,i_all,'psitt',subname)

END SUBROUTINE KStrans_p
