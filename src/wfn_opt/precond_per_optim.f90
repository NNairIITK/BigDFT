subroutine precong_per(n1,n2,n3,nseg_c,nvctr_c,nseg_f,nvctr_f,keyg,keyv, &
     ncong,cprecr,hx,hy,hz,x)
  use module_base
  ! Solves (KE+cprecr*I)*xx=yy by conjugate gradient method
  ! x is the right hand side on input and the solution on output
  implicit none
  integer, intent(in) :: n1,n2,n3,ncong
  integer, intent(in) :: nseg_c,nvctr_c,nseg_f,nvctr_f
  real(gp), intent(in) :: hx,hy,hz,cprecr
  integer, dimension(2,nseg_c+nseg_f), intent(in) :: keyg
  integer, dimension(nseg_c+nseg_f), intent(in) :: keyv
  real(wp), intent(inout) ::  x(nvctr_c+7*nvctr_f)
  ! local variables
  integer, parameter :: lowfil=-14,lupfil=14
  real(gp)::scal(0:8)
  real(wp)::rmr,rmr_new,alpha,beta
  integer i,i_stat,i_all
  real(wp),allocatable::b(:),r(:),d(:)
  real(wp),allocatable::psifscf(:),ww(:)
  integer::nd1,nd2,nd3
  integer::n1f,n3f,n1b,n3b,nd1f,nd3f,nd1b,nd3b	
  real(gp), allocatable,dimension(:,:) :: af,bf,cf,ef
  integer,allocatable,dimension(:)::modul1,modul2,modul3

  ! Array sizes for the real-to-complex FFT: note that n1(there)=n1(here)+1
  ! and the same for n2,n3.
  call dimensions_fft(n1,n2,n3,nd1,nd2,nd3,n1f,n3f,n1b,n3b,nd1f,nd3f,nd1b,nd3b)

  call allocate_all

  call prepare_sdc(n1,n2,n3,modul1,modul2,modul3,af,bf,cf,ef,hx,hy,hz)
  !	initializes the wavelet scaling coefficients	
  call wscal_init_per(scal,hx,hy,hz,cprecr)
  !b=x
  call dcopy(nvctr_c+7*nvctr_f,x,1,b,1) 

  !	compute the input guess x via a Fourier transform in a cubic box.
  !	Arrays psifscf and ww serve as work arrays for the Fourier
  call prec_fft_fast(n1,n2,n3,nseg_c,nvctr_c,nseg_f,nvctr_f,keyg,keyv, &
       cprecr,hx,hy,hz,x,&
       psifscf(1),psifscf(n1+2),psifscf(n1+n2+3),ww(1),ww(nd1b*nd2*nd3*4+1),&
       ww(nd1b*nd2*nd3*4+nd1*nd2*nd3f*4+1),&
       nd1,nd2,nd3,n1f,n1b,n3f,n3b,nd1f,nd1b,nd3f,nd3b)

  call apply_hp_sd(n1,n2,n3,nseg_c,nvctr_c,nseg_f,nvctr_f,keyg,keyv, &
       cprecr,hx,hy,hz,x,d,psifscf,ww,modul1,modul2,modul3,af,bf,cf,ef) ! d:=Ax
  r=b-d


  call wscal_per(nvctr_c,nvctr_f,scal,r(1),r(nvctr_c+1),d(1),d(nvctr_c+1))
  !rmr=dot_product(r,d)
  rmr=dot(nvctr_c+7*nvctr_f,r(1),1,d(1),1)
  
  do i=1,ncong 
     !		write(*,*)i,sqrt(rmr)

     call apply_hp_sd(n1,n2,n3,nseg_c,nvctr_c,nseg_f,nvctr_f,keyg,keyv, &
          cprecr,hx,hy,hz,d,b,psifscf,ww,modul1,modul2,modul3,af,bf,cf,ef) ! b:=Ad

     !alpha=rmr/dot_product(d,b)
     alpha=rmr/dot(nvctr_c+7*nvctr_f,d(1),1,b(1),1)
     x=x+alpha*d
     r=r-alpha*b

     call wscal_per(nvctr_c,nvctr_f,scal,r(1),r(nvctr_c+1),b(1),b(nvctr_c+1))
     !rmr_new=dot_product(r,b)
     rmr_new=dot(nvctr_c+7*nvctr_f,r(1),1,b(1),1)

     beta=rmr_new/rmr
     d=b+beta*d
     rmr=rmr_new
  enddo

  call deallocate_all

contains
  subroutine allocate_all
    allocate(modul1(lowfil:n1+lupfil+ndebug),stat=i_stat)
    call memocc(i_stat,modul1,'modul1','precong_per')
    allocate(modul2(lowfil:n2+lupfil+ndebug),stat=i_stat)
    call memocc(i_stat,modul2,'modul2','precong_per')
    allocate(modul3(lowfil:n3+lupfil+ndebug),stat=i_stat)
    call memocc(i_stat,modul3,'modul3','precong_per')
    allocate(af(lowfil:lupfil,3+ndebug),stat=i_stat)
    call memocc(i_stat,af,'af','precong_per')
    allocate(bf(lowfil:lupfil,3+ndebug),stat=i_stat)
    call memocc(i_stat,bf,'bf','precong_per')
    allocate(cf(lowfil:lupfil,3+ndebug),stat=i_stat)
    call memocc(i_stat,cf,'cf','precong_per')
    allocate(ef(lowfil:lupfil,3+ndebug),stat=i_stat)
    call memocc(i_stat,ef,'ef','precong_per')

    allocate(b(nvctr_c+7*nvctr_f+ndebug),stat=i_stat)
    call memocc(i_stat,b,'b','precong_per')
    allocate(r(nvctr_c+7*nvctr_f+ndebug),stat=i_stat)
    call memocc(i_stat,r,'r','precong_per')
    allocate(d(nvctr_c+7*nvctr_f+ndebug),stat=i_stat)
    call memocc(i_stat,d,'','precong_per')
    allocate( psifscf((2*n1+2)*(2*n2+2)*(2*n3+2)+ndebug),stat=i_stat )
    call memocc(i_stat,psifscf,'psifscf','precong_per')
    allocate( ww((2*n1+2)*(2*n2+2)*(2*n3+2)+ndebug) ,stat=i_stat)
    call memocc(i_stat,ww,'ww','precong_per')
  end subroutine allocate_all

  subroutine deallocate_all

    i_all=-product(shape(modul1))*kind(modul1)
    deallocate(modul1,stat=i_stat)
    call memocc(i_stat,i_all,'modul1','last_orthon')

    i_all=-product(shape(modul2))*kind(modul2)
    deallocate(modul2,stat=i_stat)
    call memocc(i_stat,i_all,'modul2','last_orthon')

    i_all=-product(shape(modul3))*kind(modul3)
    deallocate(modul3,stat=i_stat)
    call memocc(i_stat,i_all,'modul3','last_orthon')

    i_all=-product(shape(af))*kind(af)
    deallocate(af,stat=i_stat)
    call memocc(i_stat,i_all,'af','last_orthon')

    i_all=-product(shape(bf))*kind(bf)
    deallocate(bf,stat=i_stat)
    call memocc(i_stat,i_all,'bf','last_orthon')

    i_all=-product(shape(cf))*kind(cf)
    deallocate(cf,stat=i_stat)
    call memocc(i_stat,i_all,'cf','last_orthon')

    i_all=-product(shape(ef))*kind(ef)
    deallocate(ef,stat=i_stat)
    call memocc(i_stat,i_all,'ef','last_orthon')


    i_all=-product(shape(psifscf))*kind(psifscf)
    deallocate(psifscf,stat=i_stat)
    call memocc(i_stat,i_all,'psifscf','last_orthon')

    i_all=-product(shape(ww))*kind(ww)
    deallocate(ww,stat=i_stat)
    call memocc(i_stat,i_all,'ww','last_orthon')

    i_all=-product(shape(b))*kind(b)
    deallocate(b,stat=i_stat)
    call memocc(i_stat,i_all,'b','last_orthon')

    i_all=-product(shape(r))*kind(r)
    deallocate(r,stat=i_stat)
    call memocc(i_stat,i_all,'r','last_orthon')

    i_all=-product(shape(d))*kind(d)
    deallocate(d,stat=i_stat)
    call memocc(i_stat,i_all,'d','last_orthon')
  end subroutine deallocate_all
end subroutine precong_per


subroutine prec_fft_fast(n1,n2,n3,nseg_c,nvctr_c,nseg_f,nvctr_f,keyg,keyv, &
     cprecr,hx,hy,hz,hpsi,&
     kern_k1,kern_k2,kern_k3,z1,z3,x_c,&
     nd1,nd2,nd3,n1f,n1b,n3f,n3b,nd1f,nd1b,nd3f,nd3b)
  ! Solves (KE+cprecr*I)*xx=yy by FFT in a cubic box 
  ! x_c is the right hand side on input and the solution on output
  ! This version uses work arrays kern_k1-kern_k3 and z allocated elsewhere
  use module_base
  implicit none 
  integer, intent(in) :: n1,n2,n3
  integer,intent(in)::nd1,nd2,nd3
  integer,intent(in)::n1f,n3f,n1b,n3b,nd1f,nd3f,nd1b,nd3b	
  integer, intent(in) :: nseg_c,nvctr_c,nseg_f,nvctr_f
  real(gp), intent(in) :: hx,hy,hz,cprecr
  integer, dimension(2,nseg_c+nseg_f), intent(in) :: keyg
  integer, dimension(nseg_c+nseg_f), intent(in) :: keyv
  real(wp), intent(inout) ::  hpsi(nvctr_c+7*nvctr_f) 

  !work arrays
  real(gp):: kern_k1(0:n1),kern_k2(0:n2),kern_k3(0:n3)
  real(wp),dimension(0:n1,0:n2,0:n3):: x_c! in and out of Fourier preconditioning
  real(wp)::z1(2,nd1b,nd2,nd3,2)! work array
  real(wp)::z3(2,nd1,nd2,nd3f,2)! work array

  call wscal_f(nvctr_f,hpsi(nvctr_c+1),hx,hy,hz,cprecr)

  call make_kernel(n1,hx,kern_k1)
  call make_kernel(n2,hy,kern_k2)
  call make_kernel(n3,hz,kern_k3)

  call uncompress_c(hpsi,x_c,keyg(1,1),keyv(1),nseg_c,nvctr_c,n1,n2,n3)

  call  hit_with_kernel(x_c,z1,z3,kern_k1,kern_k2,kern_k3,n1+1,n2+1,n3+1,nd1,nd2,nd3,&
       n1f,n1b,n3f,n3b,nd1f,nd1b,nd3f,nd3b,cprecr)

  call   compress_c(hpsi,x_c,keyg(1,1),keyv(1),nseg_c,nvctr_c,n1,n2,n3)

end subroutine prec_fft_fast

subroutine prec_fft(n1,n2,n3, &
     nseg_c,nvctr_c,nseg_f,nvctr_f,keyg,keyv, &
     cprecr,hx,hy,hz,hpsi)
  ! Solves (KE+cprecr*I)*xx=yy by conjugate gradient method
  ! hpsi is the right hand side on input and the solution on output
  use module_base
  implicit none 
  integer, intent(in) :: n1,n2,n3
  integer, intent(in) :: nseg_c,nvctr_c,nseg_f,nvctr_f
  real(gp), intent(in) :: hx,hy,hz,cprecr
  integer, dimension(2,nseg_c+nseg_f), intent(in) :: keyg
  integer, dimension(nseg_c+nseg_f), intent(in) :: keyv
  real(wp), intent(inout) ::  hpsi(nvctr_c+7*nvctr_f)
  !local variables
  integer nd1,nd2,nd3,i_stat,i_all
  real(gp), dimension(:), allocatable :: kern_k1,kern_k2,kern_k3
  real(wp), dimension(:,:,:), allocatable :: x_c! in and out of Fourier preconditioning
  real(wp), dimension(:,:,:,:,:), allocatable::z1,z3 ! work array for FFT
  integer::n1f,n3f,n1b,n3b,nd1f,nd3f,nd1b,nd3b	

  ! Array sizes for the real-to-complex FFT: note that n1(there)=n1(here)+1
  ! and the same for n2,n3.
  call dimensions_fft(n1,n2,n3,nd1,nd2,nd3,n1f,n3f,n1b,n3b,nd1f,nd3f,nd1b,nd3b)

  call allocate_all

  ! diagonally precondition the wavelet part  
  call wscal_f(nvctr_f,hpsi(nvctr_c+1),hx,hy,hz,cprecr)

  call make_kernel(n1,hx,kern_k1)
  call make_kernel(n2,hy,kern_k2)
  call make_kernel(n3,hz,kern_k3)

  call uncompress_c(hpsi,x_c,keyg(1,1),keyv(1),nseg_c,nvctr_c,n1,n2,n3)

  !	solve the helmholtz equation for the scfunction part  
  call  hit_with_kernel(x_c,z1,z3,kern_k1,kern_k2,kern_k3,n1+1,n2+1,n3+1,nd1,nd2,nd3,&
       n1f,n1b,n3f,n3b,nd1f,nd1b,nd3f,nd3b,cprecr)

  call   compress_c(hpsi,x_c,keyg(1,1),keyv(1),nseg_c,nvctr_c,n1,n2,n3)

  call deallocate_all

contains
  subroutine allocate_all
    allocate(kern_k1(0:n1+ndebug),stat=i_stat)
    call memocc(i_stat,kern_k1,'kern_k1','prec_fft')
    allocate(kern_k2(0:n2+ndebug),stat=i_stat)
    call memocc(i_stat,kern_k2,'kern_k2','prec_fft')
    allocate(kern_k3(0:n3+ndebug),stat=i_stat)
    call memocc(i_stat,kern_k3,'kern_k3','prec_fft')
    allocate(z1(2,nd1b,nd2,nd3,2+ndebug),stat=i_stat) ! work array for fft
    call memocc(i_stat,z1,'z1','prec_fft')
    allocate(z3(2,nd1,nd2,nd3f,2+ndebug),stat=i_stat) ! work array for fft
    call memocc(i_stat,z3,'z3','prec_fft')
    allocate(x_c(0:n1,0:n2,0:n3+ndebug),stat=i_stat)
    call memocc(i_stat,x_c,'x_c','prec_fft')
  end subroutine allocate_all
  subroutine deallocate_all
    i_all=-product(shape(z1))*kind(z1)
    deallocate(z1,stat=i_stat)
    call memocc(i_stat,i_all,'z1','last_orthon')

    i_all=-product(shape(z3))*kind(z3)
    deallocate(z3,stat=i_stat)
    call memocc(i_stat,i_all,'z3','last_orthon')

    i_all=-product(shape(kern_k1))*kind(kern_k1)
    deallocate(kern_k1,stat=i_stat)
    call memocc(i_stat,i_all,'kern_k1','last_orthon')

    i_all=-product(shape(kern_k2))*kind(kern_k2)
    deallocate(kern_k2,stat=i_stat)
    call memocc(i_stat,i_all,'kern_k2','last_orthon')

    i_all=-product(shape(kern_k3))*kind(kern_k3)
    deallocate(kern_k3,stat=i_stat)
    call memocc(i_stat,i_all,'kern_k3','last_orthon')

    i_all=-product(shape(x_c))*kind(x_c)
    deallocate(x_c,stat=i_stat)
    call memocc(i_stat,i_all,'x_c','last_orthon')

  end subroutine deallocate_all
end subroutine prec_fft


subroutine apply_hp(n1,n2,n3, &
     nseg_c,nvctr_c,nseg_f,nvctr_f,keyg,keyv, &
     cprecr,hx,hy,hz,x,y,psifscf,ww)
  !	Applies the operator (KE+cprecr*I)*x=y
  !	array x is input, array y is output
  use module_base
  implicit none
  integer, intent(in) :: n1,n2,n3
  integer, intent(in) :: nseg_c,nvctr_c,nseg_f,nvctr_f
  real(gp), intent(in) :: hx,hy,hz,cprecr
  integer, dimension(2,nseg_c+nseg_f), intent(in) :: keyg
  integer, dimension(nseg_c+nseg_f), intent(in) :: keyv
  real(wp), intent(in) ::  x(nvctr_c+7*nvctr_f)  
  real(wp), intent(out) ::  y(nvctr_c+7*nvctr_f)

  real(gp) hgridh(3)	
  real(wp),dimension((2*n1+2)*(2*n2+2)*(2*n3+2))::ww,psifscf
  call uncompress_per(n1,n2,n3,nseg_c,nvctr_c,keyg(1,1),keyv(1),   &
       nseg_f,nvctr_f,keyg(1,nseg_c+1),keyv(nseg_c+1),   &
       x(1),x(nvctr_c+1),psifscf,ww)

  hgridh(1)=hx*.5_gp
  hgridh(2)=hy*.5_gp
  hgridh(3)=hz*.5_gp
  call convolut_kinetic_per_c(2*n1+1,2*n2+1,2*n3+1,hgridh,psifscf,ww,cprecr)

  call compress_per(n1,n2,n3,nseg_c,nvctr_c,keyg(1,1),keyv(1),   & 
       nseg_f,nvctr_f,keyg(1,nseg_c+1),keyv(nseg_c+1),   & 
       ww,y(1),y(nvctr_c+1),psifscf)
end subroutine apply_hp


subroutine wscal_f(mvctr_f,psi_f,hx,hy,hz,c)
  ! multiplies a wavefunction psi_c,psi_f (in vector form) with a scaling vector (scal)
  use module_base
  implicit none
  integer,intent(in)::mvctr_f
  real(gp),intent(in)::c,hx,hy,hz

  real(wp)::psi_f(7,mvctr_f)
  real(gp)::scal(7),hh(3)
  !WAVELET AND SCALING FUNCTION SECOND DERIVATIVE FILTERS, diagonal elements
  real(gp),PARAMETER::B2=24.8758460293923314_gp,A2=3.55369228991319019_gp

  integer i

  hh(1)=.5_gp/hx**2
  hh(2)=.5_gp/hy**2
  hh(3)=.5_gp/hz**2

  scal(1)=1._gp/(b2*hh(1)+a2*hh(2)+a2*hh(3)+c)       !  2 1 1
  scal(2)=1._gp/(a2*hh(1)+b2*hh(2)+a2*hh(3)+c)       !  1 2 1
  scal(3)=1._gp/(b2*hh(1)+b2*hh(2)+a2*hh(3)+c)       !  2 2 1
  scal(4)=1._gp/(a2*hh(1)+a2*hh(2)+b2*hh(3)+c)       !  1 1 2
  scal(5)=1._gp/(b2*hh(1)+a2*hh(2)+b2*hh(3)+c)       !  2 1 2
  scal(6)=1._gp/(a2*hh(1)+b2*hh(2)+b2*hh(3)+c)       !  1 2 2
  scal(7)=1._gp/(b2*hh(1)+b2*hh(2)+b2*hh(3)+c)       !  2 2 2

  do i=1,mvctr_f
     psi_f(1,i)=psi_f(1,i)*scal(1)       !  2 1 1
     psi_f(2,i)=psi_f(2,i)*scal(2)       !  1 2 1
     psi_f(3,i)=psi_f(3,i)*scal(3)       !  2 2 1
     psi_f(4,i)=psi_f(4,i)*scal(4)       !  1 1 2
     psi_f(5,i)=psi_f(5,i)*scal(5)       !  2 1 2
     psi_f(6,i)=psi_f(6,i)*scal(6)       !  1 2 2
     psi_f(7,i)=psi_f(7,i)*scal(7)       !  2 2 2
  enddo

end subroutine wscal_f


subroutine wscal_per(mvctr_c,mvctr_f,scal,psi_c_in,psi_f_in,psi_c_out,psi_f_out)
  ! multiplies a wavefunction psi_c,psi_f (in vector form) with a scaling vector (scal)
  use module_base
  implicit none
  integer,intent(in)::mvctr_c,mvctr_f
  real(gp),intent(in)::scal(0:8)
  real(wp),intent(in)::psi_c_in(mvctr_c),psi_f_in(7,mvctr_f)
  real(wp),intent(out)::psi_c_out(mvctr_c),psi_f_out(7,mvctr_f)

  integer i

  do i=1,mvctr_c
     psi_c_out(i)=psi_c_in(i)*scal(0)           !  1 1 1
  enddo

  do i=1,mvctr_f
     psi_f_out(1,i)=psi_f_in(1,i)*scal(1)       !  2 1 1
     psi_f_out(2,i)=psi_f_in(2,i)*scal(2)       !  1 2 1
     psi_f_out(3,i)=psi_f_in(3,i)*scal(3)       !  2 2 1
     psi_f_out(4,i)=psi_f_in(4,i)*scal(4)       !  1 1 2
     psi_f_out(5,i)=psi_f_in(5,i)*scal(5)       !  2 1 2
     psi_f_out(6,i)=psi_f_in(6,i)*scal(6)       !  1 2 2
     psi_f_out(7,i)=psi_f_in(7,i)*scal(7)       !  2 2 2
  enddo

end subroutine wscal_per


subroutine wscal_init_per(scal,hx,hy,hz,c)
  !	initialization for the array scal in the subroutine wscal_per 	
  use module_base
  implicit none
  real(wp),intent(in)::c,hx,hy,hz
  real(wp),PARAMETER::B2=24.8758460293923314D0,A2=3.55369228991319019D0
  real(wp),intent(out)::scal(0:8)
  real(gp)::hh(3)

  hh(1)=.5_wp/hx**2
  hh(2)=.5_wp/hy**2
  hh(3)=.5_wp/hz**2

  scal(0)=1._wp/(a2*hh(1)+a2*hh(2)+a2*hh(3)+c)       !  1 1 1
  scal(1)=1._wp/(b2*hh(1)+a2*hh(2)+a2*hh(3)+c)       !  2 1 1
  scal(2)=1._wp/(a2*hh(1)+b2*hh(2)+a2*hh(3)+c)       !  1 2 1
  scal(3)=1._wp/(b2*hh(1)+b2*hh(2)+a2*hh(3)+c)       !  2 2 1
  scal(4)=1._wp/(a2*hh(1)+a2*hh(2)+b2*hh(3)+c)       !  1 1 2
  scal(5)=1._wp/(b2*hh(1)+a2*hh(2)+b2*hh(3)+c)       !  2 1 2
  scal(6)=1._wp/(a2*hh(1)+b2*hh(2)+b2*hh(3)+c)       !  1 2 2
  scal(7)=1._wp/(b2*hh(1)+b2*hh(2)+b2*hh(3)+c)       !  2 2 2

end subroutine wscal_init_per

