subroutine prec_fft_slab(n1,n2,n3, &
     nseg_c,nvctr_c,nseg_f,nvctr_f,keyg,keyv, &
     cprecr,hx,hy,hz,hpsi)
  ! Solves (KE+cprecr*I)*xx=yy by conjugate gradient method
  ! hpsi is the right hand side on input and the solution on output
  use module_base
  implicit none 
  integer, intent(in) :: n1,n2,n3
  integer, intent(in) :: nseg_c,nvctr_c,nseg_f,nvctr_f
  real(kind=8), intent(in) :: hx,hy,hz,cprecr
  integer, dimension(2,nseg_c+nseg_f), intent(in) :: keyg
  integer, dimension(nseg_c+nseg_f), intent(in) :: keyv
  real(kind=8), intent(inout) ::  hpsi(nvctr_c+7*nvctr_f)
  !local variables
  integer nd1,nd2,nd3,i_stat,i_all
  real(kind=8)::hgrid(3)
  real(kind=8), dimension(:), allocatable :: kern_k1,kern_k3
  real(kind=8), dimension(:,:,:), allocatable :: x_c! in and out of Fourier preconditioning
  real(kind=8), dimension(:,:,:,:,:), allocatable::z(:,:,:,:,:) ! work array for FFT
  nd1=n1+2;        nd2=n2+2;    nd3=n3+2

  call allocate_all
  
! diagonally precondition the wavelet part  
  call wscal_f(nvctr_f,hpsi(nvctr_c+1),hx,hy,hz,cprecr)

  call make_kernel_slab(n1,hx,kern_k1)
  call make_kernel_slab(n3,hz,kern_k3)

  call uncompress_c(hpsi,x_c,keyg(1,1),keyv(1),nseg_c,nvctr_c,n1,n2,n3)

  hgrid(1)=hx
  hgrid(2)=hy
  hgrid(3)=hz
  
!    solve the helmholtz equation for the scfunction part  
!  call  hit_with_kernel(x_c,z,kern_k1,kern_k2,kern_k3,n1,n2,n3,nd1,nd2,nd3,cprecr)
  call hit_with_kernel_slab(x_c,z,kern_k1,kern_k3,n1,n2,n3,cprecr,hgrid)    

  call compress_c(hpsi,x_c,keyg(1,1),keyv(1),nseg_c,nvctr_c,n1,n2,n3)

  call deallocate_all

  contains

      subroutine allocate_all
          implicit none
          allocate(kern_k1(0:n1),stat=i_stat)
          call memocc(i_stat,product(shape(kern_k1))*kind(kern_k1),'kern_k1','prec_fft')
          allocate(kern_k3(0:n3),stat=i_stat)
          call memocc(i_stat,product(shape(kern_k3))*kind(kern_k3),'kern_k3','prec_fft')
          allocate(z(2,nd1,nd2,nd3,2),stat=i_stat) ! work array for fft
          call memocc(i_stat,product(shape(z))*kind(z),'z','prec_fft')
          allocate(x_c(0:n1,0:n2,0:n3),stat=i_stat)
          call memocc(i_stat,product(shape(x_c))*kind(x_c),'x_c','prec_fft')
      end subroutine allocate_all

      subroutine deallocate_all
          implicit none
          integer :: i_all,i_stat
          i_all=-product(shape(z))*kind(z)
          deallocate(z,stat=i_stat)
          call memocc(i_stat,i_all,'z','last_orthon')
          
          i_all=-product(shape(kern_k1))*kind(kern_k1)
          deallocate(kern_k1,stat=i_stat)
          call memocc(i_stat,i_all,'kern_k1','last_orthon')
          
          i_all=-product(shape(kern_k3))*kind(kern_k3)
          deallocate(kern_k3,stat=i_stat)
          call memocc(i_stat,i_all,'kern_k3','last_orthon')
          
          i_all=-product(shape(x_c))*kind(x_c)
          deallocate(x_c,stat=i_stat)
          call memocc(i_stat,i_all,'x_c','last_orthon')
        
      end subroutine deallocate_all
end subroutine prec_fft_slab

