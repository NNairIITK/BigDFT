! Compresses a psig wavefunction into psi_c,psi_f form
subroutine compress(n1,n2,n3,nl1,nu1,nl2,nu2,nl3,nu3, & 
     mseg_c,mvctr_c,keyg_c,keyv_c,  & 
     mseg_f,mvctr_f,keyg_f,keyv_f,  & 
     psig,psi_c,psi_f)
  use module_base
  implicit none
  integer, intent(in) :: n1,n2,n3,nl1,nu1,nl2,nu2,nl3,nu3
  integer, intent(in) :: mseg_c,mvctr_c,mseg_f,mvctr_f
  integer, dimension(mseg_c), intent(in) :: keyv_c
  integer, dimension(mseg_f), intent(in) :: keyv_f
  integer, dimension(2,mseg_c), intent(in) :: keyg_c
  integer, dimension(2,mseg_f), intent(in) :: keyg_f
  real(wp), dimension(nl1:nu1,2,nl2:nu2,2,nl3:nu3,2), intent(in) :: psig
  real(wp), dimension(mvctr_c), intent(out) :: psi_c
  real(wp), dimension(7,mvctr_f), intent(out) :: psi_f
  !local variables
  integer :: iseg,jj,j0,j1,ii,i1,i2,i3,i0,i

  ! coarse part
  do iseg=1,mseg_c
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
        psi_c(i-i0+jj)=psig(i,1,i2,1,i3,1)
     enddo
  enddo

  ! fine part
  do iseg=1,mseg_f
     jj=keyv_f(iseg)
     j0=keyg_f(1,iseg)
     j1=keyg_f(2,iseg)
     ii=j0-1
     i3=ii/((n1+1)*(n2+1))
     ii=ii-i3*(n1+1)*(n2+1)
     i2=ii/(n1+1)
     i0=ii-i2*(n1+1)
     i1=i0+j1-j0
     do i=i0,i1
        psi_f(1,i-i0+jj)=psig(i,2,i2,1,i3,1)
        psi_f(2,i-i0+jj)=psig(i,1,i2,2,i3,1)
        psi_f(3,i-i0+jj)=psig(i,2,i2,2,i3,1)
        psi_f(4,i-i0+jj)=psig(i,1,i2,1,i3,2)
        psi_f(5,i-i0+jj)=psig(i,2,i2,1,i3,2)
        psi_f(6,i-i0+jj)=psig(i,1,i2,2,i3,2)
        psi_f(7,i-i0+jj)=psig(i,2,i2,2,i3,2)
     enddo
  enddo

END SUBROUTINE compress

! Expands the compressed wavefunction in vector form (psi_c,psi_f) into the psig format
subroutine uncompress_forstandard_short(n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,  & 
     mseg_c,mvctr_c,keyg_c,keyv_c,  & 
     mseg_f,mvctr_f,keyg_f,keyv_f,  & 
     scal,psi_c,psi_f,psig_c,psig_f)
  use module_base
  implicit none
  integer, intent(in) :: n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,mseg_c,mvctr_c,mseg_f,mvctr_f
  integer, dimension(mseg_c), intent(in) :: keyv_c
  integer, dimension(mseg_f), intent(in) :: keyv_f
  integer, dimension(2,mseg_c), intent(in) :: keyg_c
  integer, dimension(2,mseg_f), intent(in) :: keyg_f
  real(wp), dimension(0:3), intent(in) :: scal
  real(wp), dimension(mvctr_c), intent(in) :: psi_c
  real(wp), dimension(7,mvctr_f), intent(in) :: psi_f
  real(wp), dimension(0:n1,0:n2,0:n3), intent(out) :: psig_c
  real(wp), dimension(7,nfl1:nfu1,nfl2:nfu2,nfl3:nfu3), intent(out) :: psig_f
  !local variables
  integer :: iseg,jj,j0,j1,ii,i1,i2,i3,i0,i

  ! coarse part
  do iseg=1,mseg_c
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
        psig_c(i,i2,i3)=psi_c(i-i0+jj)*scal(0)
     enddo
  enddo

  ! fine part
  do iseg=1,mseg_f
     jj=keyv_f(iseg)
     j0=keyg_f(1,iseg)
     j1=keyg_f(2,iseg)
     ii=j0-1
     i3=ii/((n1+1)*(n2+1))
     ii=ii-i3*(n1+1)*(n2+1)
     i2=ii/(n1+1)
     i0=ii-i2*(n1+1)
     i1=i0+j1-j0
     do i=i0,i1
        psig_f(1,i,i2,i3)=psi_f(1,i-i0+jj)*scal(1)
        
        psig_f(2,i,i2,i3)=psi_f(2,i-i0+jj)*scal(1)
        
        psig_f(3,i,i2,i3)=psi_f(3,i-i0+jj)*scal(2)
        psig_f(4,i,i2,i3)=psi_f(4,i-i0+jj)*scal(1)
        
        psig_f(5,i,i2,i3)=psi_f(5,i-i0+jj)*scal(2)
        psig_f(6,i,i2,i3)=psi_f(6,i-i0+jj)*scal(2)
        psig_f(7,i,i2,i3)=psi_f(7,i-i0+jj)*scal(3)
     enddo
  enddo

end subroutine uncompress_forstandard_short

! Expands the compressed wavefunction in vector form (psi_c,psi_f) into the psig format
subroutine uncompress_forstandard(n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,  & 
     mseg_c,mvctr_c,keyg_c,keyv_c,  & 
     mseg_f,mvctr_f,keyg_f,keyv_f,  & 
     scal,psi_c,psi_f,psig_c,psig_f,&
     x_f1,x_f2,x_f3)
  use module_base
  implicit none
  integer, intent(in) :: n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,mseg_c,mvctr_c,mseg_f,mvctr_f
  integer, dimension(mseg_c), intent(in) :: keyv_c
  integer, dimension(mseg_f), intent(in) :: keyv_f
  integer, dimension(2,mseg_c), intent(in) :: keyg_c
  integer, dimension(2,mseg_f), intent(in) :: keyg_f
  real(wp), dimension(0:3), intent(in) :: scal
  real(wp), dimension(mvctr_c), intent(in) :: psi_c
  real(wp), dimension(7,mvctr_f), intent(in) :: psi_f
  real(wp), dimension(0:n1,0:n2,0:n3), intent(out) :: psig_c
  real(wp), dimension(7,nfl1:nfu1,nfl2:nfu2,nfl3:nfu3), intent(out) :: psig_f
  real(wp), dimension(nfl1:nfu1,nfl2:nfu2,nfl3:nfu3), intent(out) :: x_f1
  real(wp), dimension(nfl2:nfu2,nfl1:nfu1,nfl3:nfu3), intent(out) :: x_f2
  real(wp), dimension(nfl3:nfu3,nfl1:nfu1,nfl2:nfu2), intent(out) :: x_f3
  !local variables
  integer :: iseg,jj,j0,j1,ii,i1,i2,i3,i0,i
 

  ! coarse part
  do iseg=1,mseg_c
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
        psig_c(i,i2,i3)=psi_c(i-i0+jj)*scal(0)
     enddo
  enddo

  ! fine part
  do iseg=1,mseg_f
     jj=keyv_f(iseg)
     j0=keyg_f(1,iseg)
     j1=keyg_f(2,iseg)
     ii=j0-1
     i3=ii/((n1+1)*(n2+1))
     ii=ii-i3*(n1+1)*(n2+1)
     i2=ii/(n1+1)
     i0=ii-i2*(n1+1)
     i1=i0+j1-j0
     do i=i0,i1
        psig_f(1,i,i2,i3)=psi_f(1,i-i0+jj)*scal(1)
        x_f1(i,i2,i3)=psig_f(1,i,i2,i3)
        
        psig_f(2,i,i2,i3)=psi_f(2,i-i0+jj)*scal(1)
        x_f2(i2,i,i3)=psig_f(2,i,i2,i3)
        
        psig_f(3,i,i2,i3)=psi_f(3,i-i0+jj)*scal(2)
        psig_f(4,i,i2,i3)=psi_f(4,i-i0+jj)*scal(1)
        x_f3(i3,i,i2)=psig_f(4,i,i2,i3)
        
        psig_f(5,i,i2,i3)=psi_f(5,i-i0+jj)*scal(2)
        psig_f(6,i,i2,i3)=psi_f(6,i-i0+jj)*scal(2)
        psig_f(7,i,i2,i3)=psi_f(7,i-i0+jj)*scal(3)
     enddo
  enddo


end subroutine uncompress_forstandard

! Compresses a psig wavefunction into psi_c,psi_f form
subroutine compress_forstandard(n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,  & 
     mseg_c,mvctr_c,keyg_c,keyv_c,  & 
     mseg_f,mvctr_f,keyg_f,keyv_f,  & 
     scal,psig_c,psig_f,psi_c,psi_f)
  use module_base
  implicit none
  integer, intent(in) :: n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,mseg_c,mvctr_c,mseg_f,mvctr_f
  integer, dimension(mseg_c), intent(in) :: keyv_c
  integer, dimension(mseg_f), intent(in) :: keyv_f
  integer, dimension(2,mseg_c), intent(in) :: keyg_c
  integer, dimension(2,mseg_f), intent(in) :: keyg_f
  real(wp), dimension(0:3), intent(in) :: scal
  real(wp), dimension(0:n1,0:n2,0:n3), intent(in) :: psig_c
  real(wp), dimension(7,nfl1:nfu1,nfl2:nfu2,nfl3:nfu3), intent(in) :: psig_f
  real(wp), dimension(mvctr_c), intent(out) :: psi_c
  real(wp), dimension(7,mvctr_f), intent(out) :: psi_f
  !local variables
  integer :: iseg,jj,j0,j1,ii,i1,i2,i3,i0,i

  ! coarse part
  do iseg=1,mseg_c
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
        psi_c(i-i0+jj)=psig_c(i,i2,i3)*scal(0)
     enddo
  enddo

  ! fine part
  do iseg=1,mseg_f
     jj=keyv_f(iseg)
     j0=keyg_f(1,iseg)
     j1=keyg_f(2,iseg)
     ii=j0-1
     i3=ii/((n1+1)*(n2+1))
     ii=ii-i3*(n1+1)*(n2+1)
     i2=ii/(n1+1)
     i0=ii-i2*(n1+1)
     i1=i0+j1-j0
     do i=i0,i1
        psi_f(1,i-i0+jj)=psig_f(1,i,i2,i3)*scal(1)
        psi_f(2,i-i0+jj)=psig_f(2,i,i2,i3)*scal(1)
        psi_f(3,i-i0+jj)=psig_f(3,i,i2,i3)*scal(2)
        psi_f(4,i-i0+jj)=psig_f(4,i,i2,i3)*scal(1)
        psi_f(5,i-i0+jj)=psig_f(5,i,i2,i3)*scal(2)
        psi_f(6,i-i0+jj)=psig_f(6,i,i2,i3)*scal(2)
        psi_f(7,i-i0+jj)=psig_f(7,i,i2,i3)*scal(3)
     enddo
  enddo

end subroutine compress_forstandard

! Compresses a wavefunction that is given in terms of fine scaling functions (psifscf) into 
! the retained coarse scaling functions and wavelet coefficients (psi_c,psi_f)
subroutine compress_per(n1,n2,n3,nseg_c,nvctr_c,keyg_c,keyv_c,  & 
     nseg_f,nvctr_f,keyg_f,keyv_f,  & 
     psifscf,psi_c,psi_f,psig)
  use module_base
  implicit none
  integer, intent(in) :: n1,n2,n3,nseg_c,nvctr_c,nseg_f,nvctr_f
  integer, dimension(nseg_c), intent(in) :: keyv_c
  integer, dimension(nseg_f), intent(in) :: keyv_f
  integer, dimension(2,nseg_c), intent(in) :: keyg_c
  integer, dimension(2,nseg_f), intent(in) :: keyg_f
  real(wp), dimension((2*n1+2)*(2*n2+2)*(2*n3+2)), intent(in) :: psifscf
  real(wp), dimension(0:n1,2,0:n2,2,0:n3,2), intent(inout) :: psig
  real(wp), dimension(nvctr_c), intent(out) :: psi_c
  real(wp), dimension(7,nvctr_f), intent(out) :: psi_f
  !local variables
  integer :: iseg,jj,j0,j1,ii,i1,i2,i3,i0,i

  ! decompose wavelets into coarse scaling functions and wavelets
  call analyse_per_self(n1,n2,n3,psifscf,psig)

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
        psi_c(i-i0+jj)=psig(i,1,i2,1,i3,1)
     enddo
  enddo

  ! fine part
  do iseg=1,nseg_f
     jj=keyv_f(iseg)
     j0=keyg_f(1,iseg)
     j1=keyg_f(2,iseg)
     ii=j0-1
     i3=ii/((n1+1)*(n2+1))
     ii=ii-i3*(n1+1)*(n2+1)
     i2=ii/(n1+1)
     i0=ii-i2*(n1+1)
     i1=i0+j1-j0
     do i=i0,i1
        psi_f(1,i-i0+jj)=psig(i,2,i2,1,i3,1)
        psi_f(2,i-i0+jj)=psig(i,1,i2,2,i3,1)
        psi_f(3,i-i0+jj)=psig(i,2,i2,2,i3,1)
        psi_f(4,i-i0+jj)=psig(i,1,i2,1,i3,2)
        psi_f(5,i-i0+jj)=psig(i,2,i2,1,i3,2)
        psi_f(6,i-i0+jj)=psig(i,1,i2,2,i3,2)
        psi_f(7,i-i0+jj)=psig(i,2,i2,2,i3,2)
     enddo
  enddo

end subroutine compress_per

subroutine uncompress_per(n1,n2,n3,nseg_c,nvctr_c,keyg_c,keyv_c,  & 
     nseg_f,nvctr_f,keyg_f,keyv_f,  & 
     psi_c,psi_f,psifscf,psig)
  ! Expands the compressed wavefunction in vector form (psi_c,psi_f) 
  ! into fine scaling functions (psifscf)
  use module_base
  implicit none
  integer, intent(in) :: n1,n2,n3,nseg_c,nvctr_c,nseg_f,nvctr_f
  integer, dimension(nseg_c), intent(in) :: keyv_c
  integer, dimension(nseg_f), intent(in) :: keyv_f
  integer, dimension(2,nseg_c), intent(in) :: keyg_c
  integer, dimension(2,nseg_f), intent(in) :: keyg_f
  real(wp), dimension(nvctr_c), intent(in) :: psi_c
  real(wp), dimension(7,nvctr_f), intent(in) :: psi_f
  real(wp), dimension(0:n1,2,0:n2,2,0:n3,2), intent(inout) :: psig
  real(wp), dimension((2*n1+2)*(2*n2+2)*(2*n3+2)), intent(out) :: psifscf
  !local variables
  integer :: iseg,jj,j0,j1,ii,i1,i2,i3,i0,i

  call razero(8*(n1+1)*(n2+1)*(n3+1),psig)

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
        psig(i,1,i2,1,i3,1)=psi_c(i-i0+jj)
     enddo
  enddo

  ! fine part
  do iseg=1,nseg_f
     jj=keyv_f(iseg)
     j0=keyg_f(1,iseg)
     j1=keyg_f(2,iseg)
     ii=j0-1
     i3=ii/((n1+1)*(n2+1))
     ii=ii-i3*(n1+1)*(n2+1)
     i2=ii/(n1+1)
     i0=ii-i2*(n1+1)
     i1=i0+j1-j0
     do i=i0,i1
        psig(i,2,i2,1,i3,1)=psi_f(1,i-i0+jj)
        psig(i,1,i2,2,i3,1)=psi_f(2,i-i0+jj)
        psig(i,2,i2,2,i3,1)=psi_f(3,i-i0+jj)
        psig(i,1,i2,1,i3,2)=psi_f(4,i-i0+jj)
        psig(i,2,i2,1,i3,2)=psi_f(5,i-i0+jj)
        psig(i,1,i2,2,i3,2)=psi_f(6,i-i0+jj)
        psig(i,2,i2,2,i3,2)=psi_f(7,i-i0+jj)
     enddo
  enddo

  ! calculate fine scaling functions
  call synthese_per_self(n1,n2,n3,psig,psifscf)

end subroutine uncompress_per

! Compresses a wavefunction that is given in terms of fine scaling functions (psifscf) into 
! the retained coarse scaling functions and wavelet coefficients (psi_c,psi_f)
subroutine compress_slab(n1,n2,n3,nseg_c,nvctr_c,keyg_c,keyv_c,  & 
     nseg_f,nvctr_f,keyg_f,keyv_f,  & 
     psifscf,psi_c,psi_f,psig)
  use module_base
  implicit none
  integer, intent(in) :: n1,n2,n3,nseg_c,nvctr_c,nseg_f,nvctr_f
  integer, dimension(nseg_c), intent(in) :: keyv_c
  integer, dimension(nseg_f), intent(in) :: keyv_f
  integer, dimension(2,nseg_c), intent(in) :: keyg_c
  integer, dimension(2,nseg_f), intent(in) :: keyg_f
  real(wp), dimension((2*n1+2)*(2*n2+16)*(2*n3+2)), intent(in) :: psifscf
  real(wp), dimension(0:n1,2,0:n2,2,0:n3,2), intent(inout) :: psig
  real(wp), dimension(nvctr_c), intent(out) :: psi_c
  real(wp), dimension(7,nvctr_f), intent(out) :: psi_f
  !local variables
  integer :: iseg,jj,j0,j1,ii,i1,i2,i3,i0,i
  real(wp),allocatable:: ww(:)

  ! decompose wavelets into coarse scaling functions and wavelets

  call analyse_slab_self(n1,n2,n3,psifscf,psig)

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
        psi_c(i-i0+jj)=psig(i,1,i2,1,i3,1)
     enddo
  enddo

  ! fine part
  do iseg=1,nseg_f
     jj=keyv_f(iseg)
     j0=keyg_f(1,iseg)
     j1=keyg_f(2,iseg)
     ii=j0-1
     i3=ii/((n1+1)*(n2+1))
     ii=ii-i3*(n1+1)*(n2+1)
     i2=ii/(n1+1)
     i0=ii-i2*(n1+1)
     i1=i0+j1-j0
     do i=i0,i1
        psi_f(1,i-i0+jj)=psig(i,2,i2,1,i3,1)
        psi_f(2,i-i0+jj)=psig(i,1,i2,2,i3,1)
        psi_f(3,i-i0+jj)=psig(i,2,i2,2,i3,1)
        psi_f(4,i-i0+jj)=psig(i,1,i2,1,i3,2)
        psi_f(5,i-i0+jj)=psig(i,2,i2,1,i3,2)
        psi_f(6,i-i0+jj)=psig(i,1,i2,2,i3,2)
        psi_f(7,i-i0+jj)=psig(i,2,i2,2,i3,2)
     enddo
  enddo

end subroutine compress_slab

subroutine uncompress_slab(n1,n2,n3,nseg_c,nvctr_c,keyg_c,keyv_c,  & 
     nseg_f,nvctr_f,keyg_f,keyv_f,  & 
     psi_c,psi_f,psifscf,psig)
  ! Expands the compressed wavefunction in vector form (psi_c,psi_f) 
  ! into fine scaling functions (psifscf)
  use module_base
  implicit none
  integer, intent(in) :: n1,n2,n3,nseg_c,nvctr_c,nseg_f,nvctr_f
  integer, dimension(nseg_c), intent(in) :: keyv_c
  integer, dimension(nseg_f), intent(in) :: keyv_f
  integer, dimension(2,nseg_c), intent(in) :: keyg_c
  integer, dimension(2,nseg_f), intent(in) :: keyg_f
  real(wp), dimension(nvctr_c), intent(in) :: psi_c
  real(wp), dimension(7,nvctr_f), intent(in) :: psi_f
  real(wp), dimension(0:n1,2,0:n2,2,0:n3,2), intent(inout) :: psig
  real(wp), dimension((2*n1+2)*(2*n2+16)*(2*n3+2)), intent(out) :: psifscf
  !local variables
  integer :: iseg,jj,j0,j1,ii,i1,i2,i3,i0,i
  real(wp),allocatable:: ww(:)

  call razero(8*(n1+1)*(n2+1)*(n3+1),psig)

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
        psig(i,1,i2,1,i3,1)=psi_c(i-i0+jj)
     enddo
  enddo

  ! fine part
  do iseg=1,nseg_f
     jj=keyv_f(iseg)
     j0=keyg_f(1,iseg)
     j1=keyg_f(2,iseg)
     ii=j0-1
     i3=ii/((n1+1)*(n2+1))
     ii=ii-i3*(n1+1)*(n2+1)
     i2=ii/(n1+1)
     i0=ii-i2*(n1+1)
     i1=i0+j1-j0
     do i=i0,i1
        psig(i,2,i2,1,i3,1)=psi_f(1,i-i0+jj)
        psig(i,1,i2,2,i3,1)=psi_f(2,i-i0+jj)
        psig(i,2,i2,2,i3,1)=psi_f(3,i-i0+jj)
        psig(i,1,i2,1,i3,2)=psi_f(4,i-i0+jj)
        psig(i,2,i2,1,i3,2)=psi_f(5,i-i0+jj)
        psig(i,1,i2,2,i3,2)=psi_f(6,i-i0+jj)
        psig(i,2,i2,2,i3,2)=psi_f(7,i-i0+jj)
     enddo
  enddo

  ! calculate fine scaling functions
  
  call synthese_slab_self(n1,n2,n3,psig,psifscf)

end subroutine uncompress_slab
