
subroutine preconditionall(iproc,nproc,norb,norbp,n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,hgrid,&
     ncong,nseg_c,nseg_f,nvctr_c,nvctr_f,keyg,keyv,eval,&
     ibyz_c,ibxz_c,ibxy_c,ibyz_f,ibxz_f,ibxy_f,hpsi)
  ! Calls the preconditioner for each orbital treated by the processor
  implicit real*8 (a-h,o-z)
  dimension ibyz_c(2,0:n2,0:n3),ibxz_c(2,0:n1,0:n3),ibxy_c(2,0:n1,0:n2)
  dimension ibyz_f(2,0:n2,0:n3),ibxz_f(2,0:n1,0:n3),ibxy_f(2,0:n1,0:n2)
  dimension keyg(2,nseg_c+nseg_f),keyv(nseg_c+nseg_f)
  dimension hpsi(nvctr_c+7*nvctr_f,norbp),eval(norb)
  
  call timing(iproc,'Precondition  ','ON')
  
  do iorb=iproc*norbp+1,min((iproc+1)*norbp,norb)
     
     cprecr=-eval(iorb)
     !       write(*,*) 'cprecr',iorb,cprecr!
     call precong(iorb,n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3, &
          nseg_c,nvctr_c,nseg_f,nvctr_f,keyg,keyv, &
          ncong,cprecr,hgrid,ibyz_c,ibxz_c,ibxy_c,ibyz_f,ibxz_f,ibxy_f,hpsi(1,iorb-iproc*norbp))
     
  enddo
  
  call timing(iproc,'Precondition  ','OF')
  
  return
end subroutine preconditionall

subroutine precong(iorb,n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3, &
     nseg_c,nvctr_c,nseg_f,nvctr_f,keyg,keyv, &
     ncong,cprecr,hgrid,ibyz_c,ibxz_c,ibxy_c,ibyz_f,ibxz_f,ibxy_f,hpsi)
  ! Solves (KE+cprecr*I)*xx=yy by conjugate gradient method
  ! hpsi is the right hand side on input and the solution on output
  implicit real*8 (a-h,o-z)
  dimension ibyz_c(2,0:n2,0:n3),ibxz_c(2,0:n1,0:n3),ibxy_c(2,0:n1,0:n2)
  dimension ibyz_f(2,0:n2,0:n3),ibxz_f(2,0:n1,0:n3),ibxy_f(2,0:n1,0:n2)
  dimension hpsi(nvctr_c+7*nvctr_f),scal(0:3),residues(ncong)
  dimension keyg(2,nseg_c+nseg_f),keyv(nseg_c+nseg_f)
  allocatable rpsi(:),ppsi(:),wpsi(:),spsi(:)

!***********************************************************************************************
  allocatable xpsig_c(:,:,:),ypsig_c(:,:,:)
  allocatable xpsig_f(:,:,:,:),ypsig_f(:,:,:,:)
  real*8,allocatable,dimension(:,:,:)::x_f1,x_f2,x_f3 ! input
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
!***********************************************************************************************
  do i=1,nvctr_c+7*nvctr_f
     spsi(i)=hpsi(i)
  enddo
  !        do i=0,3
  !        scal(i)=1.d0
  !        enddo
  !        call  CALC_GRAD_REZA(n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3, &
  !              nseg_c,nvctr_c,keyg,keyv,nseg_f,nvctr_f,keyg(1,nseg_c+1),keyv(nseg_c+1), &
  !              scal,cprecr,hgrid,ibyz_c,ibxz_c,ibxy_c,ibyz_f,ibxz_f,ibxy_f,hpsi,hpsi(nvctr_c+1),wpsi,wpsi(nvctr_c+1))
  !        tt=0.d0
  !	do i=1,nvctr_c+7*nvctr_f
  !	tt=tt+(wpsi(i)-spsi(i))**2
  !        enddo
  !        write(*,*) 'initial residue',iorb,sqrt(tt)
  !! checkend



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

  call  CALC_GRAD_REZA(n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3, &
       nseg_c,nvctr_c,keyg,keyv,nseg_f,nvctr_f,keyg(1,nseg_c+1),keyv(nseg_c+1), &
       scal,cprecr,hgrid,ibyz_c,ibxz_c,ibxy_c,ibyz_f,ibxz_f,ibxy_f,hpsi,&
	   hpsi(nvctr_c+1),wpsi,wpsi(nvctr_c+1),&
		 xpsig_c,xpsig_f,ypsig_c,ypsig_f,&
		 x_f1,x_f2,x_f3)

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
     call CALC_GRAD_REZA(n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3, &
          nseg_c,nvctr_c,keyg,keyv,nseg_f,nvctr_f,keyg(1,nseg_c+1),keyv(nseg_c+1), &
          scal,cprecr,hgrid,ibyz_c,ibxz_c,ibxy_c,ibyz_f,ibxz_f,&
		  ibxy_f,ppsi,ppsi(nvctr_c+1),wpsi,wpsi(nvctr_c+1),&
		 xpsig_c,xpsig_f,ypsig_c,ypsig_f,&
		 x_f1,x_f2,x_f3)
     
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

!	if (alpha1.lt.tol) goto 1010
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
  call  CALC_GRAD_REZA(n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3, &
       nseg_c,nvctr_c,keyg,keyv,nseg_f,nvctr_f,keyg(1,nseg_c+1),keyv(nseg_c+1), &
       scal,cprecr,hgrid,ibyz_c,ibxz_c,ibxy_c,ibyz_f,ibxz_f,&
	   ibxy_f,hpsi,hpsi(nvctr_c+1),wpsi,wpsi(nvctr_c+1),&
	 xpsig_c,xpsig_f,ypsig_c,ypsig_f,&
	 x_f1,x_f2,x_f3)

	 
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

!***********************************************************************************************
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
!***********************************************************************************************
  return
end subroutine precong

SUBROUTINE CALC_GRAD_REZA(n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3, & 
     nseg_c,nvctr_c,keyg_c,keyv_c,nseg_f,nvctr_f,keyg_f,keyv_f, &
     scal,cprecr,hgrid,ibyz_c,ibxz_c,ibxy_c,ibyz_f,ibxz_f,ibxy_f,xpsi_c,xpsi_f,ypsi_c,ypsi_f,&
 xpsig_c,xpsig_f,ypsig_c,ypsig_f,&
 x_f1,x_f2,x_f3)
  ! ypsi = (1/2) \Nabla^2 xpsi
  implicit real*8 (a-h,o-z)        
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

  call Convolkinetic(n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3, & 
       cprecr,hgrid,ibyz_c,ibxz_c,ibxy_c,ibyz_f,ibxz_f,ibxy_f,xpsig_c,&
	   xpsig_f,ypsig_c,ypsig_f,x_f1,x_f2,x_f3)

  call compress_forstandard(n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,  &
       nseg_c,nvctr_c,keyg_c,keyv_c,  & 
       nseg_f,nvctr_f,keyg_f,keyv_f,  & 
       scal,ypsig_c,ypsig_f,ypsi_c,ypsi_f)

END SUBROUTINE CALC_GRAD_REZA
