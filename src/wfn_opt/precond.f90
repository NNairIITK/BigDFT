
subroutine preconditionall(iproc,nproc,norb,norbp,n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,hgrid,&
     ncong,nseg_c,nseg_f,nvctr_c,nvctr_f,keyg,keyv,eval,&
     ibyz_c,ibxz_c,ibxy_c,ibyz_f,ibxz_f,ibxy_f,hpsi)
  ! Calls the preconditioner for each orbital treated by the processor
  implicit real(kind=8) (a-h,o-z)
  dimension ibyz_c(2,0:n2,0:n3),ibxz_c(2,0:n1,0:n3),ibxy_c(2,0:n1,0:n2)
  dimension ibyz_f(2,0:n2,0:n3),ibxz_f(2,0:n1,0:n3),ibxy_f(2,0:n1,0:n2)
  dimension keyg(2,nseg_c+nseg_f),keyv(nseg_c+nseg_f)
  dimension hpsi(nvctr_c+7*nvctr_f,norbp),eval(norb)
  
  do iorb=iproc*norbp+1,min((iproc+1)*norbp,norb)
     
     cprecr=-eval(iorb)
     !       write(*,*) 'cprecr',iorb,cprecr!
     call precong(iorb,n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3, &
          nseg_c,nvctr_c,nseg_f,nvctr_f,keyg,keyv, &
          ncong,cprecr,hgrid,ibyz_c,ibxz_c,ibxy_c,ibyz_f,ibxz_f,ibxy_f,hpsi(1,iorb-iproc*norbp))
     
  enddo

end subroutine preconditionall

subroutine precong(iorb,n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3, &
     nseg_c,nvctr_c,nseg_f,nvctr_f,keyg,keyv, &
     ncong,cprecr,hgrid,ibyz_c,ibxz_c,ibxy_c,ibyz_f,ibxz_f,ibxy_f,hpsi)
  ! Solves (KE+cprecr*I)*xx=yy by conjugate gradient method
  ! hpsi is the right hand side on input and the solution on output
  implicit real(kind=8) (a-h,o-z)
  dimension ibyz_c(2,0:n2,0:n3),ibxz_c(2,0:n1,0:n3),ibxy_c(2,0:n1,0:n2)
  dimension ibyz_f(2,0:n2,0:n3),ibxz_f(2,0:n1,0:n3),ibxy_f(2,0:n1,0:n2)
  dimension hpsi(nvctr_c+7*nvctr_f),scal(0:3),residues(ncong)
  dimension keyg(2,nseg_c+nseg_f),keyv(nseg_c+nseg_f)
  allocatable rpsi(:),ppsi(:),wpsi(:),spsi(:)
  logical, parameter :: newmethod=.true.

!***********************************************************************************************
  allocatable :: xpsig_c(:,:,:), ypsig_c(:,:,:)
  allocatable :: xpsig_f(:,:,:,:), ypsig_f(:,:,:,:)
  real(kind=8), allocatable, dimension(:,:,:) :: x_f1,x_f2,x_f3 ! input
!***********************************************************************************************
  !       WAVELET AND SCALING FUNCTION SECOND DERIVATIVE FILTERS
  PARAMETER(B2=24.8758460293923314D0,A2=3.55369228991319019D0)
  LOGICAL,PARAMETER::INGUESS_ON=.TRUE.

  !      The input guess consists of diagonal preconditioning of the original gradient.
  !       In contrast to older version, not only the wavelet part and the scfunction
  !       part are multiplied by different factors, but the scfunction part is 
  !       subjected to wavelet analysis with periodic boundaries. Then the wavelets
  !       on different scales are multiplied by different factors and backward wavelet !       transformed to scaling functions.
  !
  !       The new input guess is turned on if the parameter INGUESS_ON
  !       has value .TRUE.
  !       
  tt=sum(hpsi)
  if(abs(tt)/=0.0d0) then
!     print *,"PREcond",abs(tt)
 

  allocate(rpsi(nvctr_c+7*nvctr_f),stat=i_stat)
  call memocc(i_stat,product(shape(rpsi))*kind(rpsi),'rpsi','precong')
  allocate(ppsi(nvctr_c+7*nvctr_f),stat=i_stat)
  call memocc(i_stat,product(shape(ppsi))*kind(ppsi),'ppsi','precong')
  allocate(wpsi(nvctr_c+7*nvctr_f),stat=i_stat)
  call memocc(i_stat,product(shape(wpsi))*kind(wpsi),'wpsi','precong')


  !! check initial residue of original equation
  allocate(spsi(nvctr_c+7*nvctr_f),stat=i_stat)
  call memocc(i_stat,product(shape(spsi))*kind(spsi),'spsi','precong')

!***********************************************************************************************
  do i=1,nvctr_c+7*nvctr_f
     spsi(i)=hpsi(i)
  enddo


  FAC_H=1.d0/hgrid**2
  H0=    1.5D0*A2*FAC_H;    H1=(A2+B2*.5D0)*FAC_H
  H2=(A2*.5D0+B2)*FAC_H;    H3=    1.5D0*B2*FAC_H

  scal(0)=sqrt(1.D0/(H0+cprecr)) ;  scal(1)=sqrt(1.D0/(H1+cprecr)) 
  scal(2)=sqrt(1.D0/(H2+cprecr)) ;  scal(3)=sqrt(1.D0/(H3+cprecr))


  IF (INGUESS_ON) THEN
     !          The right hand side is temporarily stored in the rpsi array        
     rpsi=hpsi           
     !          and preconditioned with D^{-1/2} as usual:
     call  wscalv(nvctr_c,nvctr_f,scal,rpsi,rpsi(nvctr_c+1))

     !          hpsi is now diagonally preconditioned with alexey's old preconditioner;
     !          inside the diagonal preconditioner a factor of D^{1/2} was added
     !          to make the overall factor D^{-1/2} again
     CALL prec_diag(n1,n2,n3,hgrid,nseg_c,nvctr_c,nvctr_f,&
          keyg,keyv,hpsi,hpsi(nvctr_c+1),cprecr,scal,A2,B2)
  ELSE
     !          assume as input guess x=y
     !          hpsi is preconditioned with D^{-1/2} as usual
     call  wscalv(nvctr_c,nvctr_f,scal,hpsi,hpsi(nvctr_c+1))
  ENDIF

  if (newmethod) then
     allocate(xpsig_c(0:n1,0:n2,0:n3),stat=i_stat)
     call memocc(i_stat,product(shape(xpsig_c))*kind(xpsig_c),'xpsig_c','calc_grad_reza')
     allocate(xpsig_f(7,nfl1:nfu1,nfl2:nfu2,nfl3:nfu3),stat=i_stat)
     call memocc(i_stat,product(shape(xpsig_f))*kind(xpsig_f),'xpsig_f','calc_grad_reza')
     allocate(ypsig_c(0:n1,0:n2,0:n3),stat=i_stat)
     call memocc(i_stat,product(shape(ypsig_c))*kind(ypsig_c),'ypsig_c','calc_grad_reza')
     allocate(ypsig_f(7,nfl1:nfu1,nfl2:nfu2,nfl3:nfu3),stat=i_stat)
     call memocc(i_stat,product(shape(ypsig_f))*kind(ypsig_f),'ypsig_f','calc_grad_reza')

     allocate(x_f1(nfl1:nfu1,nfl2:nfu2,nfl3:nfu3),stat=i_stat)
     call memocc(i_stat,product(shape(x_f1))*kind(x_f1),'x_f1','applylocpotkinall')
     allocate(x_f2(nfl2:nfu2,nfl1:nfu1,nfl3:nfu3),stat=i_stat)
     call memocc(i_stat,product(shape(x_f2))*kind(x_f2),'x_f2','applylocpotkinall')
     allocate(x_f3(nfl3:nfu3,nfl1:nfu1,nfl2:nfu2),stat=i_stat)
     call memocc(i_stat,product(shape(x_f3))*kind(x_f3),'x_f3','applylocpotkinall')


     call razero((nfu1-nfl1+1)*(nfu2-nfl2+1)*(nfu3-nfl3+1),x_f1)
     call razero((nfu1-nfl1+1)*(nfu2-nfl2+1)*(nfu3-nfl3+1),x_f2)
     call razero((nfu1-nfl1+1)*(nfu2-nfl2+1)*(nfu3-nfl3+1),x_f3)

     call razero((n1+1)*(n2+1)*(n3+1),xpsig_c)
     call razero(7*(nfu1-nfl1+1)*(nfu2-nfl2+1)*(nfu3-nfl3+1),xpsig_f)

     call razero((n1+1)*(n2+1)*(n3+1),ypsig_c)
     call razero(7*(nfu1-nfl1+1)*(nfu2-nfl2+1)*(nfu3-nfl3+1),ypsig_f)

     call CALC_GRAD_REZA(n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3, &
          nseg_c,nvctr_c,keyg,keyv,nseg_f,nvctr_f,keyg(1,nseg_c+1),keyv(nseg_c+1), &
          scal,cprecr,hgrid,ibyz_c,ibxz_c,ibxy_c,ibyz_f,ibxz_f,ibxy_f,hpsi,&
          hpsi(nvctr_c+1),wpsi,wpsi(nvctr_c+1),&
          xpsig_c,xpsig_f,ypsig_c,ypsig_f,&
          x_f1,x_f2,x_f3)
  else
      call CALC_GRAD_REZA_prev(n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3, &
           nseg_c,nvctr_c,keyg,keyv,nseg_f,nvctr_f,keyg(1,nseg_c+1),keyv(nseg_c+1), &
           scal,cprecr,hgrid,ibyz_c,ibxz_c,ibxy_c,ibyz_f,ibxz_f,ibxy_f,&
           hpsi,hpsi(nvctr_c+1),wpsi,wpsi(nvctr_c+1))

  end if

  IF (INGUESS_ON) THEN 
     do i=1,nvctr_c+7*nvctr_f
        tt=wpsi(i)-rpsi(i)  ! rpsi instead of hpsi: alexey
        rpsi(i)=tt
        ppsi(i)=tt
     enddo

  ELSE
     do i=1,nvctr_c+7*nvctr_f
        tt=wpsi(i)-hpsi(i)  ! normal
        rpsi(i)=tt
        ppsi(i)=tt
     enddo
  ENDIF


  do icong=2,ncong
     if (newmethod) then
        call CALC_GRAD_REZA(n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3, &
             nseg_c,nvctr_c,keyg,keyv,nseg_f,nvctr_f,keyg(1,nseg_c+1),keyv(nseg_c+1), &
             scal,cprecr,hgrid,ibyz_c,ibxz_c,ibxy_c,ibyz_f,ibxz_f,&
             ibxy_f,ppsi,ppsi(nvctr_c+1),wpsi,wpsi(nvctr_c+1),&
             xpsig_c,xpsig_f,ypsig_c,ypsig_f,&
             x_f1,x_f2,x_f3)
     else
        call CALC_GRAD_REZA_prev(n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3, &
          nseg_c,nvctr_c,keyg,keyv,nseg_f,nvctr_f,keyg(1,nseg_c+1),keyv(nseg_c+1), &
          scal,cprecr,hgrid,ibyz_c,ibxz_c,ibxy_c,ibyz_f,ibxz_f,ibxy_f,ppsi,ppsi(nvctr_c+1),&
          wpsi,wpsi(nvctr_c+1))
     end if
     
     alpha1=0.d0 ; alpha2=0.d0
     do i=1,nvctr_c+7*nvctr_f
        alpha1=alpha1+rpsi(i)*rpsi(i)
        alpha2=alpha2+rpsi(i)*wpsi(i)
     enddo
     residues(icong)=alpha1
     alpha=alpha1/alpha2
     
     do i=1,nvctr_c+7*nvctr_f
        hpsi(i)=hpsi(i)-alpha*ppsi(i)
        rpsi(i)=rpsi(i)-alpha*wpsi(i)
     end do

!    if (alpha1.lt.tol) goto 1010
     if (icong.ge.ncong) goto 1010
     
     beta1=0.d0 ; beta2=0.d0
     do i=1,nvctr_c+7*nvctr_f
        beta1=beta1+rpsi(i)*wpsi(i)
        beta2=beta2+ppsi(i)*wpsi(i)
     enddo
     beta=beta1/beta2
     
     do i=1,nvctr_c+7*nvctr_f
        ppsi(i)=rpsi(i)-beta*ppsi(i)
     end do

  end do
1010 continue

  !  D^{-1/2} times solution
  call wscalv(nvctr_c,nvctr_f,scal,hpsi,hpsi(nvctr_c+1))
  
  !write(*,'(i4,(100(1x,e8.2)))') iorb,(residues(icong),icong=2,ncong)
  
  ! check final residue of original equation
  do i=0,3
     scal(i)=1.d0
  enddo
  if (newmethod) then
     call  CALC_GRAD_REZA(n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3, &
          nseg_c,nvctr_c,keyg,keyv,nseg_f,nvctr_f,keyg(1,nseg_c+1),keyv(nseg_c+1), &
          scal,cprecr,hgrid,ibyz_c,ibxz_c,ibxy_c,ibyz_f,ibxz_f,&
          ibxy_f,hpsi,hpsi(nvctr_c+1),wpsi,wpsi(nvctr_c+1),&
          xpsig_c,xpsig_f,ypsig_c,ypsig_f,&
          x_f1,x_f2,x_f3)
  else
      call CALC_GRAD_REZA_prev(n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3, &
       nseg_c,nvctr_c,keyg,keyv,nseg_f,nvctr_f,keyg(1,nseg_c+1),keyv(nseg_c+1), &
       scal,cprecr,hgrid,ibyz_c,ibxz_c,ibxy_c,ibyz_f,ibxz_f,ibxy_f,&
       hpsi,hpsi(nvctr_c+1),wpsi,wpsi(nvctr_c+1))
  end if

     
  tt=0.d0
  do i=1,nvctr_c+7*nvctr_f
     tt=tt+(wpsi(i)-spsi(i))**2
  enddo
  !write(*,'(1x,a,1x,i0,1x,1pe13.6)') 'Precond, final residue',iorb,sqrt(tt)
  i_all=-product(shape(spsi))*kind(spsi)
  deallocate(spsi,stat=i_stat)
  call memocc(i_stat,i_all,'spsi','precong')
  ! checkend

  i_all=-product(shape(rpsi))*kind(rpsi)
  deallocate(rpsi,stat=i_stat)
  call memocc(i_stat,i_all,'rpsi','precong')
  i_all=-product(shape(ppsi))*kind(ppsi)
  deallocate(ppsi,stat=i_stat)
  call memocc(i_stat,i_all,'ppsi','precong')
  i_all=-product(shape(wpsi))*kind(wpsi)
  deallocate(wpsi,stat=i_stat)
  call memocc(i_stat,i_all,'wpsi','precong')

  if (newmethod) then
     i_all=-product(shape(xpsig_c))*kind(xpsig_c)
     deallocate(xpsig_c,stat=i_stat)
     call memocc(i_stat,i_all,'xpsig_c','calc_grad_reza')

     i_all=-product(shape(ypsig_c))*kind(ypsig_c)
     deallocate(ypsig_c,stat=i_stat)
     call memocc(i_stat,i_all,'ypsig_c','calc_grad_reza')

     i_all=-product(shape(xpsig_f))*kind(xpsig_f)
     deallocate(xpsig_f,stat=i_stat)
     call memocc(i_stat,i_all,'xpsig_f','calc_grad_reza')

     i_all=-product(shape(ypsig_f))*kind(ypsig_f)
     deallocate(ypsig_f,stat=i_stat)
     call memocc(i_stat,i_all,'ypsig_f','calc_grad_reza')

     i_all=-product(shape(x_f1))*kind(x_f1)
     deallocate(x_f1,stat=i_stat)
     call memocc(i_stat,i_all,'x_f1','applylocpotkinall')

     i_all=-product(shape(x_f2))*kind(x_f2)
     deallocate(x_f2,stat=i_stat)
     call memocc(i_stat,i_all,'x_f2','applylocpotkinall')

     i_all=-product(shape(x_f3))*kind(x_f3)
     deallocate(x_f3,stat=i_stat)
     call memocc(i_stat,i_all,'x_f3','applylocpotkinall')
  end if

!  print *,"AFTcond",sum(hpsi)
  else
!     print *,'no preconditioning on empty orbital'
  end if

end subroutine precong

SUBROUTINE CALC_GRAD_REZA(n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3, & 
     nseg_c,nvctr_c,keyg_c,keyv_c,nseg_f,nvctr_f,keyg_f,keyv_f, &
     scal,cprecr,hgrid,ibyz_c,ibxz_c,ibxy_c,ibyz_f,ibxz_f,ibxy_f,xpsi_c,xpsi_f,ypsi_c,ypsi_f,&
 xpsig_c,xpsig_f,ypsig_c,ypsig_f,&
 x_f1,x_f2,x_f3)
  ! ypsi = (1/2) \Nabla^2 xpsi
  implicit real(kind=8) (a-h,o-z)        
  dimension keyg_c(2,nseg_c),keyv_c(nseg_c),keyg_f(2,nseg_f),keyv_f(nseg_f)
  dimension xpsi_c(nvctr_c),xpsi_f(7,nvctr_f),scal(0:3)
  dimension ypsi_c(nvctr_c),ypsi_f(7,nvctr_f)
  dimension ibyz_c(2,0:n2,0:n3),ibxz_c(2,0:n1,0:n3),ibxy_c(2,0:n1,0:n2)
  dimension ibyz_f(2,0:n2,0:n3),ibxz_f(2,0:n1,0:n3),ibxy_f(2,0:n1,0:n2)
!***********************************************************************************************
   dimension  xpsig_c(0:n1,0:n2,0:n3),xpsig_f(7,nfl1:nfu1,nfl2:nfu2,nfl3:nfu3)
   dimension  ypsig_c(0:n1,0:n2,0:n3),ypsig_f(7,nfl1:nfu1,nfl2:nfu2,nfl3:nfu3)
   dimension  x_f1(nfl1:nfu1,nfl2:nfu2,nfl3:nfu3),x_f2(nfl2:nfu2,nfl1:nfu1,nfl3:nfu3)
   dimension  x_f3(nfl3:nfu3,nfl1:nfu1,nfl2:nfu2)
!***********************************************************************************************

  call uncompress_forstandard(n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,  &
       nseg_c,nvctr_c,keyg_c,keyv_c,  & 
       nseg_f,nvctr_f,keyg_f,keyv_f,  & 
       scal,xpsi_c,xpsi_f,xpsig_c,xpsig_f,x_f1,x_f2,x_f3)

!!$  ypsig_c=xpsig_c
!!$  ypsig_f=xpsig_f
  call Convolkinetic(n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3, & 
       cprecr,hgrid,ibyz_c,ibxz_c,ibxy_c,ibyz_f,ibxz_f,ibxy_f,xpsig_c,&
       xpsig_f,ypsig_c,ypsig_f,x_f1,x_f2,x_f3)

  call compress_forstandard(n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,  &
       nseg_c,nvctr_c,keyg_c,keyv_c,  & 
       nseg_f,nvctr_f,keyg_f,keyv_f,  & 
       scal,ypsig_c,ypsig_f,ypsi_c,ypsi_f)

END SUBROUTINE CALC_GRAD_REZA


subroutine prec_diag(n1,n2,n3,hgrid,nseg_c,nvctr_c,nvctr_f,&
     keyg_c,keyv_c,hpsi_c,hpsi_f,C,scal,A2,B2)
  ! 
  !
  implicit real(kind=8) (a-h,o-z)
  dimension keyg_c(2,nseg_c),keyv_c(nseg_c),hpsi_c(nvctr_c),hpsi_f(7,nvctr_f)
  real(kind=8), allocatable, dimension(:,:,:) :: hpsip
  real(kind=8)::scal(0:3) 
  real(kind=8),parameter::atomic_length=2.d0,FAC_LEN=2.D0

  !      Number of sweeps in wavelet transformation
  !      THE BIGGEST SCALING FUNCTION STEP: atomic_length*FAC_LEN
  !      (NOT JUST ATOMIC_LENGTH, BECAUSE SO IT IS BETTER IN PRACTICE) 
  NUM_TRANS=NINT(log(atomic_length*FAC_LEN/hgrid)/log(2.d0))
  N2_NT=2**NUM_TRANS
  !write(*,'(1x,a)') 'NUMBER OF WAVELET TRANSFORMS (sweeps)',NUM_TRANS

  ! Find right leading dimensions for array


  !       ND1+1 IS THE MULTIPLE OF N2_N
  !       WHICH IS CLOSEST TO N1+1 FROM ABOVE. 
  ND1=CEILING( REAL(N1+1,KIND=8)/REAL(N2_NT,KIND=8) ) *N2_NT-1
  !       THE SAME FOR ND2,ND3.
  ND2=CEILING( REAL(N2+1,KIND=8)/REAL(N2_NT,KIND=8) ) *N2_NT-1
  ND3=CEILING( REAL(N3+1,KIND=8)/REAL(N2_NT,KIND=8) ) *N2_NT-1

  !write(*,'(3(1x,a,i0))')'ND1=',ND1,'ND2=',ND2,'ND3=',ND3

  allocate(hpsip(0:nd1,0:nd2,0:nd3),stat=i_stat)
  call memocc(i_stat,product(shape(hpsip))*kind(hpsip),'hpsip','prec_diag')

  HPSIP=0.D0

  ! coarse part
  do iseg=1,nseg_c
     jj=keyv_c(iseg)
     j0=keyg_c(1,iseg)
     j1=keyg_c(2,iseg)
     ii=j0-1
     i3=ii/((n1+1)*(n2+1))
     ii=ii-i3*(n1+1)*(n2+1)
     i2=ii/(n1+1)
     i0=ii-i2*(n1+1)
     i1=i0+j1-j0
     do i=i0,i1
        hpsip(i,i2,i3)=hpsi_c(i-i0+jj)
     enddo
  enddo

  FAC_H=1.D0/((HGRID*REAL(N2_NT,KIND=8))**2)

  H0=    1.5D0*A2*FAC_H;    H1=(A2+B2*.5D0)*FAC_H
  H2=(A2*.5D0+B2)*FAC_H;    H3=    1.5D0*B2*FAC_H

  !       FORWARD TRANSFORM THE COARSE SCALING FUNCTIONS NUM_TRANS TIMES
  CALL ANA_REPEATED_PER(ND1,ND2,ND3,HPSIP,NUM_TRANS,NN1,NN2,NN3) 

  NNN1=NN1; NNN2=NN2; NNN3=NN3 

  !       DIAGONALLY PRECONDITION THE RESULTING COARSE WAVELETS
  CALL PRECOND_PROPER(ND1,ND2,ND3,HPSIP,NUM_TRANS,NNN1,NNN2,NNN3,H0,H1,H2,H3,C)

  HPSIP=HPSIP/SCAL(0) ! apply (wscal)^(-1)

  !       BACKWARD TRANSFORM THE COARSE SCALING FUNCTIONS NUM_TRANS TIMES
  CALL SYN_REPEATED_PER(ND1,ND2,ND3,HPSIP,NUM_TRANS,NN1,NN2,NN3)

  !       DIAGONALLY PRECONDITION THE FINE WAVELETS
  DO I=1,NVCTR_F
     HPSI_F(1,I)=HPSI_F(1,I)*scal(1)
     HPSI_F(2,I)=HPSI_F(2,I)*scal(1)
     HPSI_F(4,I)=HPSI_F(4,I)*scal(1)

     HPSI_F(3,I)=HPSI_F(3,I)*scal(2)
     HPSI_F(5,I)=HPSI_F(5,I)*scal(2)
     HPSI_F(6,I)=HPSI_F(6,I)*scal(2)

     HPSI_F(7,I)=HPSI_F(7,I)*scal(3)
  ENDDO

  ! coarse part
  do iseg=1,nseg_c
     jj=keyv_c(iseg)
     j0=keyg_c(1,iseg)
     j1=keyg_c(2,iseg)
     ii=j0-1
     i3=ii/((n1+1)*(n2+1))
     ii=ii-i3*(n1+1)*(n2+1)
     i2=ii/(n1+1)
     i0=ii-i2*(n1+1)
     i1=i0+j1-j0
     do i=i0,i1
        hpsi_c(i-i0+jj)=hpsip(i,i2,i3)
     enddo
  enddo

  i_all=-product(shape(hpsip))*kind(hpsip)
  deallocate(hpsip,stat=i_stat)
  call memocc(i_stat,i_all,'hpsip','prec_diag')

end subroutine prec_diag

SUBROUTINE PRECOND_PROPER(nd1,nd2,nd3,x,NUM_TRANS,N1,N2,N3,H0,H1,H2,H3,EPS)
  implicit real(kind=8) (a-h,o-z)
  dimension  x(0:nd1,0:nd2,0:nd3)


  DO I_TRANS=1,NUM_TRANS
     N1P=2*(N1+1)-1
     N2P=2*(N2+1)-1
     N3P=2*(N3+1)-1

     IF (N1P.GT.ND1) STOP 'N1 BEYOND BORDERS'
     IF (N2P.GT.ND2) STOP 'N2 BEYOND BORDERS'
     IF (N3P.GT.ND3) STOP 'N3 BEYOND BORDERS'

     N1PP=N1+1
     N2PP=N2+1
     N3PP=N3+1

     F1=1.D0/(H1+EPS); F2=1.D0/(H2+EPS);  F3=1.D0/(H3+EPS)       


     IF (I_TRANS.EQ.1) THEN 

        F0=1.D0/(H0+EPS)

        DO I3=0,N3
           I3P=I3+N3PP
           DO I2=0,N2
              I2P=I2+N2PP
              DO I1=0,N1
                 I1P=I1+N1PP

                 X(I1,I2,I3)=X(I1,I2,I3)*F0

                 X(I1P,I2,I3)=X(I1P,I2,I3)*F1
                 X(I1,I2P,I3)=X(I1,I2P,I3)*F1
                 X(I1,I2,I3P)=X(I1,I2,I3P)*F1

                 X(I1P,I2P,I3)=X(I1P,I2P,I3)*F2
                 X(I1,I2P,I3P)=X(I1,I2P,I3P)*F2
                 X(I1P,I2,I3P)=X(I1P,I2,I3P)*F2

                 X(I1P,I2P,I3P)=X(I1P,I2P,I3P)*F3

              ENDDO
           ENDDO
        ENDDO

     ELSE

        DO I3=0,N3
           I3P=I3+N3PP
           DO I2=0,N2
              I2P=I2+N2PP
              DO I1=0,N1
                 I1P=I1+N1PP

                 X(I1P,I2,I3)=X(I1P,I2,I3)*F1
                 X(I1,I2P,I3)=X(I1,I2P,I3)*F1
                 X(I1,I2,I3P)=X(I1,I2,I3P)*F1

                 X(I1P,I2P,I3)=X(I1P,I2P,I3)*F2
                 X(I1,I2P,I3P)=X(I1,I2P,I3P)*F2
                 X(I1P,I2,I3P)=X(I1P,I2,I3P)*F2

                 X(I1P,I2P,I3P)=X(I1P,I2P,I3P)*F3

              ENDDO
           ENDDO
        ENDDO

     ENDIF

     N1=N1P
     N2=N2P
     N3=N3P

     H1=H1*4.D0
     H2=H2*4.D0
     H3=H3*4.D0

  ENDDO

END SUBROUTINE PRECOND_PROPER

