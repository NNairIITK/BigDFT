! Calls the preconditioner for each orbital treated by the processor
subroutine preconditionall(geocode,iproc,nproc,norb,norbp,n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,&
     hx,hy,hz,ncong,nspinor,wfd,eval,kb,hpsi,gnrm)
  use module_base
  use module_types
  implicit none
  type(wavefunctions_descriptors), intent(in) :: wfd
  type(kinetic_bounds), intent(in) :: kb
  character(len=1), intent(in) :: geocode
  integer, intent(in) :: iproc,nproc,norb,norbp,n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3
  integer, intent(in) :: nspinor,ncong
  real(gp), intent(in) :: hx,hy,hz
  real(wp), dimension(norb), intent(in) :: eval
  real(dp), intent(out) :: gnrm
  real(wp), dimension(wfd%nvctr_c+7*wfd%nvctr_f,norbp*nspinor), intent(inout) :: hpsi
  !local variables
  integer :: iorb,inds,indo
  real(wp) :: cprecr
  real(dp) :: scpr
  real(kind=8), external :: dnrm2

  ! Preconditions all orbitals belonging to iproc
  !and calculate the norm of the residue

  ! norm of gradient
  gnrm=0.0_dp
  do iorb=iproc*norbp+1,min((iproc+1)*norbp,norb)
     indo=(iorb-1)*nspinor+1-iproc*norbp*nspinor
     !loop over the spinorial components
     do inds=indo,indo+nspinor-1

        scpr=dnrm2(wfd%nvctr_c+7*wfd%nvctr_f,hpsi(1,inds),1)
        gnrm=gnrm+scpr**2

        select case(geocode)
        case('F')
           !in this case the grid spacings are uniform
           cprecr=-eval(iorb)
           if(scpr /=0.0_dp) then
              call precong(iorb,n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3, &
                   wfd%nseg_c,wfd%nvctr_c,wfd%nseg_f,wfd%nvctr_f,wfd%keyg,wfd%keyv, &
                   ncong,cprecr,hx,kb%ibyz_c,kb%ibxz_c,kb%ibxy_c,&
                   kb%ibyz_f,kb%ibxz_f,kb%ibxy_f,hpsi(1,inds))
           end if
        case('P')
           cprecr=0.5_wp
           !           cprecr=abs(eval(iorb))
           !		   if (cprecr.lt..1_wp) cprecr=.5_wp
           if (ncong.eq.0) then
              call prec_fft(n1,n2,n3, &
                   wfd%nseg_c,wfd%nvctr_c,wfd%nseg_f,wfd%nvctr_f,wfd%keyg,wfd%keyv, &
                   cprecr,hx,hy,hz,hpsi(1,inds))
           else
              call precong_per(n1,n2,n3, &
                   wfd%nseg_c,wfd%nvctr_c,wfd%nseg_f,wfd%nvctr_f,wfd%keyg,wfd%keyv, &
                   ncong,cprecr,hx,hy,hz,hpsi(1,inds))
           endif
        case('S')
           cprecr=0.5_wp
			if (ncong.eq.0) then
				call prec_fft_slab(n1,n2,n3, &
				wfd%nseg_c,wfd%nvctr_c,wfd%nseg_f,wfd%nvctr_f,wfd%keyg,wfd%keyv, &
                cprecr,hx,hy,hz,hpsi(1,inds))
			else
        		call precong_slab(n1,n2,n3, &
				wfd%nseg_c,wfd%nvctr_c,wfd%nseg_f,wfd%nvctr_f,wfd%keyg,wfd%keyv, &
				ncong,cprecr,hx,hy,hz,hpsi(1,inds))
			endif
        end select
     end do
  enddo

end subroutine preconditionall

subroutine precong(iorb,n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3, &
     nseg_c,nvctr_c,nseg_f,nvctr_f,keyg,keyv, &
     ncong,cprecr,hgrid,ibyz_c,ibxz_c,ibxy_c,ibyz_f,ibxz_f,ibxy_f,hpsi)
  ! Solves (KE+cprecr*I)*xx=yy by conjugate gradient method
  ! hpsi is the right hand side on input and the solution on output
  use module_base
  implicit none
  !implicit real(kind=8) (a-h,o-z)
  integer, intent(in) :: iorb,n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3
  integer, intent(in) :: nseg_c,nvctr_c,nseg_f,nvctr_f,ncong
  real(gp), intent(in) :: hgrid
  real(dp), intent(in) :: cprecr
  integer, dimension(nseg_c+nseg_f), intent(in) :: keyv
  integer, dimension(2,nseg_c+nseg_f), intent(in) :: keyg
  integer, dimension(2,0:n2,0:n3), intent(in) :: ibyz_c,ibyz_f
  integer, dimension(2,0:n1,0:n3), intent(in) :: ibxz_c,ibxz_f
  integer, dimension(2,0:n1,0:n2), intent(in) :: ibxy_c,ibxy_f
  real(wp), dimension(nvctr_c+7*nvctr_f), intent(inout) :: hpsi
  !local variables
  character(len=*), parameter :: subname='precong'
  logical, parameter :: inguess_on=.true.
  !       wavelet and scaling function second derivative filters
  real(wp), parameter :: b2=24.8758460293923314_wp, a2=3.55369228991319019_wp
  integer :: i,icong,i_stat,i_all
  real(wp) :: fac_h,h0,h1,h2,h3,tt,alpha1,alpha2,alpha,beta1,beta2,beta
  real(wp), dimension(0:3) :: scal
  real(wp), dimension(:), allocatable :: rpsi,ppsi,wpsi,spsi
  real(wp), dimension(:,:,:,:), allocatable :: xpsig_f,ypsig_f
  real(wp), dimension(:,:,:), allocatable :: xpsig_c,ypsig_c,x_f1,x_f2,x_f3


  ! The input guess consists of diagonal preconditioning of the original gradient.
  ! In contrast to older version, not only the wavelet part and the scfunction
  ! part are multiplied by different factors, but the scfunction part is 
  ! subjected to wavelet analysis with periodic boundaries. Then the wavelets
  ! on different scales are multiplied by different factors and backward wavelet 
  ! transformed to scaling functions.
  !
  ! The new input guess is turned on if the parameter INGUESS_ON
  ! has value .TRUE.
  ! 
  
  allocate(rpsi(nvctr_c+7*nvctr_f+ndebug),stat=i_stat)
  call memocc(i_stat,rpsi,'rpsi',subname)
  allocate(ppsi(nvctr_c+7*nvctr_f+ndebug),stat=i_stat)
  call memocc(i_stat,ppsi,'ppsi',subname)
  allocate(wpsi(nvctr_c+7*nvctr_f+ndebug),stat=i_stat)
  call memocc(i_stat,wpsi,'wpsi',subname)

!!$  !array of initial wavefunction
!!$  allocate(spsi(nvctr_c+7*nvctr_f),stat=i_stat)
!!$  call memocc(i_stat,spsi,'spsi',subname)
!!$  do i=1,nvctr_c+7*nvctr_f
!!$     spsi(i)=hpsi(i)
!!$  enddo

  fac_h=1.0_wp/real(hgrid,wp)**2
  h0=    1.5_wp*a2*fac_h
  h1=(a2+b2*.5_wp)*fac_h
  h2=(a2*.5_wp+b2)*fac_h
  h3=    1.5_wp*b2*fac_h

  scal(0)=sqrt(1.0_wp/(h0+cprecr)) 
  scal(1)=sqrt(1.0_wp/(h1+cprecr)) 
  scal(2)=sqrt(1.0_wp/(h2+cprecr)) 
  scal(3)=sqrt(1.0_wp/(h3+cprecr))

  if (inguess_on) then
     !          the right hand side is temporarily stored in the rpsi array        
     rpsi=hpsi           
     !          and preconditioned with d^{-1/2} as usual:
     call  wscalv(nvctr_c,nvctr_f,scal,rpsi,rpsi(nvctr_c+1))

     !          hpsi is now diagonally preconditioned with alexey's old preconditioner;
     !          inside the diagonal preconditioner a factor of d^{1/2} was added
     !          to make the overall factor d^{-1/2} again
     call prec_diag(n1,n2,n3,hgrid,nseg_c,nvctr_c,nvctr_f,&
          keyg,keyv,hpsi,hpsi(nvctr_c+1),cprecr,scal,a2,b2)
  else
     !          assume as input guess x=y
     !          hpsi is preconditioned with d^{-1/2} as usual
     call  wscalv(nvctr_c,nvctr_f,scal,hpsi,hpsi(nvctr_c+1))
  endif

  !allocate work arrays
  allocate(xpsig_c(0:n1,0:n2,0:n3+ndebug),stat=i_stat)
  call memocc(i_stat,xpsig_c,'xpsig_c',subname)
  allocate(xpsig_f(7,nfl1:nfu1,nfl2:nfu2,nfl3:nfu3+ndebug),stat=i_stat)
  call memocc(i_stat,xpsig_f,'xpsig_f',subname)
  allocate(ypsig_c(0:n1,0:n2,0:n3+ndebug),stat=i_stat)
  call memocc(i_stat,ypsig_c,'ypsig_c',subname)
  allocate(ypsig_f(7,nfl1:nfu1,nfl2:nfu2,nfl3:nfu3+ndebug),stat=i_stat)
  call memocc(i_stat,ypsig_f,'ypsig_f',subname)

  allocate(x_f1(nfl1:nfu1,nfl2:nfu2,nfl3:nfu3+ndebug),stat=i_stat)
  call memocc(i_stat,x_f1,'x_f1',subname)
  allocate(x_f2(nfl2:nfu2,nfl1:nfu1,nfl3:nfu3+ndebug),stat=i_stat)
  call memocc(i_stat,x_f2,'x_f2',subname)
  allocate(x_f3(nfl3:nfu3,nfl1:nfu1,nfl2:nfu2+ndebug),stat=i_stat)
  call memocc(i_stat,x_f3,'x_f3',subname)
  
  !initalize to zero the work arrays, probably not needed
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

  loop_precond: do icong=2,ncong

     call CALC_GRAD_REZA(n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3, &
          nseg_c,nvctr_c,keyg,keyv,nseg_f,nvctr_f,keyg(1,nseg_c+1),keyv(nseg_c+1), &
          scal,cprecr,hgrid,ibyz_c,ibxz_c,ibxy_c,ibyz_f,ibxz_f,&
          ibxy_f,ppsi,ppsi(nvctr_c+1),wpsi,wpsi(nvctr_c+1),&
          xpsig_c,xpsig_f,ypsig_c,ypsig_f,&
          x_f1,x_f2,x_f3)

     alpha1=0.0_wp 
     alpha2=0.0_wp
     do i=1,nvctr_c+7*nvctr_f
        alpha1=alpha1+rpsi(i)*rpsi(i)
        alpha2=alpha2+rpsi(i)*wpsi(i)
     enddo

     !residues(icong)=alpha1
     alpha=alpha1/alpha2        

     !write(10+iorb,'(1x,i0,3(1x,1pe24.17))')icong,alpha1,alpha2,alpha

     do i=1,nvctr_c+7*nvctr_f
        hpsi(i)=hpsi(i)-alpha*ppsi(i)
        rpsi(i)=rpsi(i)-alpha*wpsi(i)
     end do

     if (icong >= ncong) exit loop_precond

     beta1=0.0_wp 
     beta2=0.0_wp

     do i=1,nvctr_c+7*nvctr_f
        beta1=beta1+rpsi(i)*wpsi(i)
        beta2=beta2+ppsi(i)*wpsi(i)
     enddo

     beta=beta1/beta2        

     do i=1,nvctr_c+7*nvctr_f
        ppsi(i)=rpsi(i)-beta*ppsi(i)
     end do

  end do loop_precond

  !  D^{-1/2} times solution
  call wscalv(nvctr_c,nvctr_f,scal,hpsi,hpsi(nvctr_c+1))

  !write(*,'(i4,(100(1x,e8.2)))') iorb,(residues(icong),icong=2,ncong)

!!$  ! check final residue of original equation
!!$  do i=0,3
!!$     scal(i)=1.d0
!!$  enddo
!!$
!!$  call CALC_GRAD_REZA(n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3, &
!!$       nseg_c,nvctr_c,keyg,keyv,nseg_f,nvctr_f,keyg(1,nseg_c+1),keyv(nseg_c+1), &
!!$       scal,cprecr,hgrid,ibyz_c,ibxz_c,ibxy_c,ibyz_f,ibxz_f,&
!!$       ibxy_f,hpsi,hpsi(nvctr_c+1),wpsi,wpsi(nvctr_c+1),&
!!$       xpsig_c,xpsig_f,ypsig_c,ypsig_f,&
!!$       x_f1,x_f2,x_f3)
!!$     
!!$  tt=0.d0
!!$  do i=1,nvctr_c+7*nvctr_f
!!$     tt=tt+(wpsi(i)-spsi(i))**2
!!$  enddo
!!$  !write(*,'(1x,a,1x,i0,1x,1pe13.6)') 'Precond, final residue',iorb,sqrt(tt)
!!$  i_all=-product(shape(spsi))*kind(spsi)
!!$  deallocate(spsi,stat=i_stat)
!!$  call memocc(i_stat,i_all,'spsi',subname)
!!$  ! checkend

  i_all=-product(shape(rpsi))*kind(rpsi)
  deallocate(rpsi,stat=i_stat)
  call memocc(i_stat,i_all,'rpsi',subname)
  i_all=-product(shape(ppsi))*kind(ppsi)
  deallocate(ppsi,stat=i_stat)
  call memocc(i_stat,i_all,'ppsi',subname)
  i_all=-product(shape(wpsi))*kind(wpsi)
  deallocate(wpsi,stat=i_stat)
  call memocc(i_stat,i_all,'wpsi',subname)


  i_all=-product(shape(xpsig_c))*kind(xpsig_c)
  deallocate(xpsig_c,stat=i_stat)
  call memocc(i_stat,i_all,'xpsig_c',subname)

  i_all=-product(shape(ypsig_c))*kind(ypsig_c)
  deallocate(ypsig_c,stat=i_stat)
  call memocc(i_stat,i_all,'ypsig_c',subname)

  i_all=-product(shape(xpsig_f))*kind(xpsig_f)
  deallocate(xpsig_f,stat=i_stat)
  call memocc(i_stat,i_all,'xpsig_f',subname)

  i_all=-product(shape(ypsig_f))*kind(ypsig_f)
  deallocate(ypsig_f,stat=i_stat)
  call memocc(i_stat,i_all,'ypsig_f',subname)

  i_all=-product(shape(x_f1))*kind(x_f1)
  deallocate(x_f1,stat=i_stat)
  call memocc(i_stat,i_all,'x_f1',subname)

  i_all=-product(shape(x_f2))*kind(x_f2)
  deallocate(x_f2,stat=i_stat)
  call memocc(i_stat,i_all,'x_f2',subname)

  i_all=-product(shape(x_f3))*kind(x_f3)
  deallocate(x_f3,stat=i_stat)
  call memocc(i_stat,i_all,'x_f3',subname)
     
end subroutine precong

! ypsi = (1/2) \Nabla^2 xpsi + cprecr xpsi
subroutine calc_grad_reza(n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3, & 
     nseg_c,nvctr_c,keyg_c,keyv_c,nseg_f,nvctr_f,keyg_f,keyv_f, &
     scal,cprecr,hgrid,ibyz_c,ibxz_c,ibxy_c,ibyz_f,ibxz_f,ibxy_f,&
     xpsi_c,xpsi_f,ypsi_c,ypsi_f,&
     xpsig_c,xpsig_f,ypsig_c,ypsig_f,x_f1,x_f2,x_f3)
  use module_base
  implicit none
  integer, intent(in) :: n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3
  integer, intent(in) :: nseg_c,nvctr_c,nseg_f,nvctr_f
  real(wp), intent(in) :: cprecr
  real(gp), intent(in) :: hgrid
  integer, dimension(nseg_c), intent(in) :: keyv_c
  integer, dimension(nseg_f), intent(in) :: keyv_f
  integer, dimension(2,nseg_c), intent(in) :: keyg_c
  integer, dimension(2,nseg_f), intent(in) :: keyg_f
  integer, dimension(2,0:n2,0:n3), intent(in) :: ibyz_c,ibyz_f
  integer, dimension(2,0:n1,0:n3), intent(in) :: ibxz_c,ibxz_f
  integer, dimension(2,0:n1,0:n2), intent(in) :: ibxy_c,ibxy_f
  real(wp), dimension(0:3), intent(in) :: scal
  real(wp), dimension(nvctr_c), intent(in) :: xpsi_c
  real(wp), dimension(7,nvctr_f), intent(in) :: xpsi_f
  real(wp), dimension(nvctr_c), intent(out) :: ypsi_c
  real(wp), dimension(7,nvctr_f), intent(out) :: ypsi_f
  real(wp), dimension(0:n1,0:n2,0:n3), intent(inout) :: xpsig_c,ypsig_c
  real(wp), dimension(7,nfl1:nfu1,nfl2:nfu2,nfl3:nfu3), intent(inout) :: xpsig_f,ypsig_f
  real(wp), dimension(nfl1:nfu1,nfl2:nfu2,nfl3:nfu3), intent(inout) :: x_f1
  real(wp), dimension(nfl2:nfu2,nfl1:nfu1,nfl3:nfu3), intent(inout) :: x_f2
  real(wp), dimension(nfl3:nfu3,nfl1:nfu1,nfl2:nfu2), intent(inout) :: x_f3

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

end subroutine calc_grad_reza


subroutine prec_diag(n1,n2,n3,hgrid,nseg_c,nvctr_c,nvctr_f,&
     keyg_c,keyv_c,hpsi_c,hpsi_f,c,scal,a2,b2)
  use module_base
  implicit none
  integer, intent(in) :: n1,n2,n3,nseg_c,nvctr_c,nvctr_f
  real(wp), intent(in) :: c,a2,b2
  real(gp), intent(in) :: hgrid
  integer, dimension(nseg_c), intent(in) :: keyv_c
  integer, dimension(2,nseg_c), intent(in) :: keyg_c
  real(wp), dimension(0:3), intent(in) :: scal
  real(wp), dimension(nvctr_c), intent(inout) :: hpsi_c
  real(wp), dimension(7,nvctr_f), intent(inout) :: hpsi_f
  !local variables
  character(len=*), parameter :: subname='prec_diag'
  real(gp), parameter ::atomic_length=2.0_gp,fac_len=2.0_gp
  integer :: num_trans,n2_nt,nd1,nd2,nd3,iseg,jj,j0,ii,i3,i2,i
  integer :: nn1,nn2,nn3,nnn1,nnn2,nnn3,i0,i_all,i_stat,i1,j1
  real(wp) :: h0,h1,h2,h3,fac_h
  real(wp), dimension(:,:,:), allocatable :: hpsip


  !      number of sweeps in wavelet transformation
  !      the biggest scaling function step: atomic_length*fac_len
  !      (not just atomic_length, because so it is better in practice) 
  num_trans=nint(log(atomic_length*fac_len/hgrid)/log(2.0_gp))
  n2_nt=2**num_trans
  !write(*,'(1x,a)') 'number of wavelet transforms (sweeps)',num_trans

  ! find right leading dimensions for array


  !       nd1+1 is the multiple of n2_n
  !       which is closest to n1+1 from above. 
  nd1=ceiling( real(n1+1,kind=8)/real(n2_nt,kind=8)) *n2_nt-1
  !       the same for nd2,nd3.
  nd2=ceiling( real(n2+1,kind=8)/real(n2_nt,kind=8)) *n2_nt-1
  nd3=ceiling( real(n3+1,kind=8)/real(n2_nt,kind=8)) *n2_nt-1

  !write(*,'(3(1x,a,i0))')'nd1=',nd1,'nd2=',nd2,'nd3=',nd3

  allocate(hpsip(0:nd1,0:nd2,0:nd3+ndebug),stat=i_stat)
  call memocc(i_stat,hpsip,'hpsip',subname)

  hpsip=0.0_wp

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

  fac_h=real(1.0_gp/((hgrid*real(n2_nt,gp))**2),wp)

  h0=    1.5_wp*a2*fac_h
  h1=(a2+b2*.5d0)*fac_h
  h2=(a2*.5_wp+b2)*fac_h
  h3=    1.5_wp*b2*fac_h

  !       forward transform the coarse scaling functions num_trans times
  call ana_repeated_per(nd1,nd2,nd3,hpsip,num_trans,nn1,nn2,nn3) 

  nnn1=nn1
  nnn2=nn2
  nnn3=nn3 

  !       diagonally precondition the resulting coarse wavelets
  call precond_proper(nd1,nd2,nd3,hpsip,num_trans,nnn1,nnn2,nnn3,h0,h1,h2,h3,c)

  hpsip=hpsip/scal(0) ! apply (wscal)^(-1)

  !       backward transform the coarse scaling functions num_trans times
  call syn_repeated_per(nd1,nd2,nd3,hpsip,num_trans,nn1,nn2,nn3)

  !       diagonally precondition the fine wavelets
  do i=1,nvctr_f
     hpsi_f(1,i)=hpsi_f(1,i)*scal(1)
     hpsi_f(2,i)=hpsi_f(2,i)*scal(1)
     hpsi_f(4,i)=hpsi_f(4,i)*scal(1)

     hpsi_f(3,i)=hpsi_f(3,i)*scal(2)
     hpsi_f(5,i)=hpsi_f(5,i)*scal(2)
     hpsi_f(6,i)=hpsi_f(6,i)*scal(2)

     hpsi_f(7,i)=hpsi_f(7,i)*scal(3)
  enddo

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
  call memocc(i_stat,i_all,'hpsip',subname)

end subroutine prec_diag

subroutine precond_proper(nd1,nd2,nd3,x,num_trans,n1,n2,n3,h0,h1,h2,h3,eps)
  use module_base
  implicit none
  integer, intent(in) :: nd1,nd2,nd3,num_trans
  integer, intent(inout) :: n1,n2,n3
  real(wp), intent(in) :: eps,h0
  real(wp), intent(inout) :: h1,h2,h3
  real(wp), dimension(0:nd1,0:nd2,0:nd3), intent(inout) :: x
  !local variables
  integer :: i_trans,n1p,n2p,n3p,n1pp,n2pp,n3pp,i1,i2,i3,i1p,i2p,i3p
  real(wp) :: f0,f1,f2,f3


  do i_trans=1,num_trans
     n1p=2*(n1+1)-1
     n2p=2*(n2+1)-1
     n3p=2*(n3+1)-1

     if (n1p.gt.nd1) stop 'n1 beyond borders'
     if (n2p.gt.nd2) stop 'n2 beyond borders'
     if (n3p.gt.nd3) stop 'n3 beyond borders'

     n1pp=n1+1
     n2pp=n2+1
     n3pp=n3+1

     f1=1.0_wp/(h1+eps)
     f2=1.0_wp/(h2+eps)
     f3=1.0_wp/(h3+eps)       

     if (i_trans == 1) then 

        f0=1.d0/(h0+eps)

        do i3=0,n3
           i3p=i3+n3pp
           do i2=0,n2
              i2p=i2+n2pp
              do i1=0,n1
                 i1p=i1+n1pp

                 x(i1,i2,i3)=x(i1,i2,i3)*f0

                 x(i1p,i2,i3)=x(i1p,i2,i3)*f1
                 x(i1,i2p,i3)=x(i1,i2p,i3)*f1
                 x(i1,i2,i3p)=x(i1,i2,i3p)*f1

                 x(i1p,i2p,i3)=x(i1p,i2p,i3)*f2
                 x(i1,i2p,i3p)=x(i1,i2p,i3p)*f2
                 x(i1p,i2,i3p)=x(i1p,i2,i3p)*f2

                 x(i1p,i2p,i3p)=x(i1p,i2p,i3p)*f3

              enddo
           enddo
        enddo

     else

        do i3=0,n3
           i3p=i3+n3pp
           do i2=0,n2
              i2p=i2+n2pp
              do i1=0,n1
                 i1p=i1+n1pp

                 x(i1p,i2,i3)=x(i1p,i2,i3)*f1
                 x(i1,i2p,i3)=x(i1,i2p,i3)*f1
                 x(i1,i2,i3p)=x(i1,i2,i3p)*f1

                 x(i1p,i2p,i3)=x(i1p,i2p,i3)*f2
                 x(i1,i2p,i3p)=x(i1,i2p,i3p)*f2
                 x(i1p,i2,i3p)=x(i1p,i2,i3p)*f2

                 x(i1p,i2p,i3p)=x(i1p,i2p,i3p)*f3

              enddo
           enddo
        enddo

     endif

     n1=n1p
     n2=n2p
     n3=n3p

     h1=h1*4.0_wp
     h2=h2*4.0_wp
     h3=h3*4.0_wp

  enddo

end subroutine precond_proper

