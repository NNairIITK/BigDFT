!> @file
!!  Routines to apply the local part of the hamiltonian
!! @author
!!    Copyright (C) 2009-2011 BigDFT group 
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 
!!
!>  Applies the local potential and kinetic energy operator to one wavefunction 
!! Input: pot,psi
!! Output: hpsi,epot,ekin
subroutine applylocpotkinone_per(n1,n2,n3, & 
     hx,hy,hz,nseg_c,nseg_f,nvctr_c,nvctr_f,keyg,keyv,  & 
     psir,psi_in,psi_out,psi,pot,hpsi,epot,ekin,npot,nspinor)
  use module_base
  use module_interfaces
  implicit none
  integer, intent(in) :: n1,n2,n3,nseg_c,nseg_f,nvctr_c,nvctr_f,npot,nspinor
  real(gp), intent(in) :: hx,hy,hz
  integer, dimension(nseg_c+nseg_f), intent(in) :: keyv
  integer, dimension(2,nseg_c+nseg_f), intent(in) :: keyg
  real(wp), dimension((2*n1+2)*(2*n2+2)*(2*n3+2),npot), intent(in) :: pot
  real(wp), dimension(nvctr_c+7*nvctr_f,nspinor), intent(in) :: psi
  real(wp), dimension((2*n1+2)*(2*n2+2)*(2*n3+2),nspinor), intent(inout) :: psir,psi_in,psi_out
  real(gp), intent(out) :: epot,ekin
  real(wp), dimension(nvctr_c+7*nvctr_f,nspinor), intent(out) :: hpsi
  !local variables
  character(len=*), parameter :: subname='applylocpotkinone_per'
  integer :: i,i_stat,i_all,idx
  real(wp) :: tt
  real(gp) :: v,p,epot_p
  real(gp), dimension(3) :: hgridh

  ! Initialisation of potential energy  
!!  epot=0.0_gp
  ekin=0.0_gp

  hgridh(1)=hx*.5_gp
  hgridh(2)=hy*.5_gp
  hgridh(3)=hz*.5_gp


  ! Wavefunction expressed everywhere in fine scaling functions (for potential and kinetic energy)
  do idx=1,nspinor
     call uncompress_per(n1,n2,n3,nseg_c,nvctr_c,keyg(1,1),keyv(1),   &
          nseg_f,nvctr_f,keyg(1,nseg_c+1),keyv(nseg_c+1),   &
          psi(1,idx),psi(nvctr_c+1,idx),psi_in(1,idx),psir(1,idx))

     ! psir serves as a work array	   
     call convolut_magic_n_per(2*n1+1,2*n2+1,2*n3+1,psi_in(1,idx),psir(1,idx),psi_out(1,idx)) 
  end do

  call apply_potential(n1,n2,n3,0,0,0,0,nspinor,npot,psir,pot,epot)

!!  !$omp parallel default(private)&
!!  !$omp shared(pot,psir,n1,n2,n3,epot)
!!
!!  epot_p=0._gp
!!  !$omp do
!!  do i=1,(2*n1+2)*(2*n2+2)*(2*n3+2)
!!     v=real(pot(i),gp)
!!     p=real(psir(i),gp)
!!     tt=pot(i)*psir(i)
!!     epot_p=epot_p+p*v*p
!!     psir(i)=tt
!!  enddo
!!  !$omp end do
!!
!!  !$omp critical
!!  epot=epot+epot_p
!!  !$omp end critical
!!
!!  !$omp end parallel

  do idx=1,nspinor
     call convolut_magic_t_per_self(2*n1+1,2*n2+1,2*n3+1,psir(1,idx),psi_out(1,idx))

     ! compute the kinetic part and add  it to psi_out
     ! the kinetic energy is calculated at the same time
     !here we should insert the treatment for k-points
     call convolut_kinetic_per_t(2*n1+1,2*n2+1,2*n3+1,hgridh,psi_in(1,idx),psi_out(1,idx),ekin)

     call compress_per(n1,n2,n3,nseg_c,nvctr_c,keyg(1,1),keyv(1),   & 
          nseg_f,nvctr_f,keyg(1,nseg_c+1),keyv(nseg_c+1),   & 
          psi_out(1,idx),hpsi(1,idx),hpsi(nvctr_c+1,idx),psir(1,idx))
  end do


END SUBROUTINE applylocpotkinone_per


!>  Applies the local potential and kinetic energy operator to one wavefunction 
!! Input: pot,psi
!! Output: hpsi,epot,ekin
subroutine applylocpotkinone_hyb(n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3, & 
     hx,hy,hz,nseg_c,nseg_f,nvctr_c,nvctr_f,keyg,keyv,  & 
     psir,psi,pot,hpsi,epot,ekin,bounds,nspinor,npot)
  use module_base
  use module_types
  use module_interfaces
  implicit none
  type(convolutions_bounds),intent(in):: bounds
  integer, intent(in) :: n1,n2,n3,nseg_c,nseg_f,nvctr_c,nvctr_f,nspinor,npot
  integer,intent(in):: nfl1,nfl2,nfl3,nfu1,nfu2,nfu3
  real(gp), intent(in) :: hx,hy,hz
  integer, dimension(nseg_c+nseg_f), intent(in) :: keyv
  integer, dimension(2,nseg_c+nseg_f), intent(in) :: keyg
  real(wp), dimension((2*n1+2)*(2*n2+2)*(2*n3+2),npot), intent(in) :: pot
  real(wp), dimension(nvctr_c+7*nvctr_f), intent(in) :: psi
  real(wp), dimension((2*n1+2)*(2*n2+2)*(2*n3+2)), intent(inout) :: psir
  real(gp), intent(out) :: epot,ekin
  real(wp), dimension(nvctr_c+7*nvctr_f), intent(out) :: hpsi
  !local variables
  integer :: i
  integer nf
  real(wp) :: tt
  real(gp) :: v,p
  real(gp), dimension(3) :: hgridh,hgrid
  real(wp),allocatable::x_f(:,:,:,:),x_c(:,:,:)
  real(wp),allocatable,dimension(:)::x_f1,x_f2,x_f3
  real(wp),allocatable,dimension(:,:,:)::y_c
  real(wp),allocatable,dimension(:,:,:,:)::y_f
  real(wp),allocatable,dimension(:)::w,ww
 

  integer i_stat,i_all
  integer nw,nww,i2,i3
  

  !these allocation should be rised up by a level
  !for the moment do not allow hybrid BC for nspinor/=1/=npot
  if (nspinor /= 1 .or. npot /=1) stop 'Complex functions not allowed for hybrid BC'

  ! Wavefunction expressed everywhere in fine scaling functions (for potential and kinetic energy)
  nf=7*(nfu1-nfl1+1)*(nfu2-nfl2+1)*(nfu3-nfl3+1)

  nw=max(4*(nfu2-nfl2+1)*(nfu3-nfl3+1)*(2*n1+2),(2*n1+2)*(n2+2)*(n3+2))
  nw=max(nw,2*(n3+1)*(n1+1)*(n2+1))	   ! for the comb_shrink_hyb_c
  nw=max(nw,4*(2*n3+2)*(nfu1-nfl1+1)*(nfu2-nfl2+1)) ! for the _f
  allocate(w(nw+ndebug),stat=i_stat)
  call memocc(i_stat,w,'w','applylocpotkinone_hyb')
  
  nww=max(2*(nfu3-nfl3+1)*(2*n1+2)*(2*n2+2),(n3+1)*(2*n1+2)*(2*n2+2))
  nww=max(nww,4*(n2+1)*(n3+1)*(n1+1))	! for the comb_shrink_hyb_c   
  nww=max(nww,2*(2*n2+2)*(2*n3+2)*(nfu1-nfl1+1)) ! for the _f
  allocate(ww(nww+ndebug),stat=i_stat)
  call memocc(i_stat,ww,'ww','applylocpotkinone_hyb')

   allocate(x_f(7,nfl1:nfu1,nfl2:nfu2,nfl3:nfu3+ndebug),stat=i_stat)
   call memocc(i_stat,x_f,'x_f','applylocpotkinone_hyb')
   allocate(x_c(0:n1,0:n2,0:n3+ndebug),stat=i_stat)
   call memocc(i_stat,x_c,'x_c ','applylocpotkinone_hyb')
   allocate(x_f1(nf+ndebug),stat=i_stat)
   call memocc(i_stat,x_f1,'x_f1','applylocpotkinone_hyb')
   allocate(x_f2(nf+ndebug),stat=i_stat)
   call memocc(i_stat,x_f2,'x_f2','applylocpotkinone_hyb')
   allocate(x_f3(nf+ndebug),stat=i_stat)
   call memocc(i_stat,x_f3,'x_f3','applylocpotkinone_hyb')
 
   allocate(y_f(7,nfl1:nfu1,nfl2:nfu2,nfl3:nfu3+ndebug),stat=i_stat)
   call memocc(i_stat,y_f,'y_f','applylocpotkinone_hyb')
   allocate(y_c(0:n1,0:n2,0:n3+ndebug),stat=i_stat)
   call memocc(i_stat,y_c,'y_c','applylocpotkinone_hyb')

  call uncompress_per_f(n1,n2,n3,nseg_c,nvctr_c,keyg(1,1),keyv(1),   &
       nseg_f,nvctr_f,keyg(1,nseg_c+1),keyv(nseg_c+1),   &
       psi(1),psi(nvctr_c+1),x_c,x_f,x_f1,x_f2,x_f3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3)

! x_c: input, psir1: output
! psir: work array
  call comb_grow_all_hybrid(n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,nw,nww,&
       w,ww,x_c,x_f,psir,bounds%gb)
  
  call apply_potential(n1,n2,n3,0,0,0,0,nspinor,npot,psir,pot,epot)

!!  epot=0.0_gp
!!  do i=1,(2*n1+2)*(2*n2+2)*(2*n3+2)
!!     v=real(pot(i),gp)
!!     p=real(psir(i),gp)
!!     tt=pot(i)*psir(i)
!!     epot=epot+p*v*p
!!     psir(i)=tt
!!  enddo

! y_c has the scfunction output of the kinetic energy operator  
!psir  : input, y_c: output, psi_in:work
 call comb_shrink_hyb(n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,ww,w,psir,y_c,y_f,bounds%sb)

  hgrid(1)=hx
  hgrid(2)=hy
  hgrid(3)=hz

! compute the kinetic part and add  it to psi_out
! the kinetic energy is calculated at the same time
  
  call convolut_kinetic_hyb_T(n1,n2,n3, &
     nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,  &
     hgrid,x_c,x_f,y_c,y_f,ekin,x_f1,x_f2,x_f3,bounds%kb%ibyz_f,&
	bounds%kb%ibxz_f,bounds%kb%ibxy_f)

  call compress_per_f(n1,n2,n3,nseg_c,nvctr_c,keyg(1,1),keyv(1),   & 
       nseg_f,nvctr_f,keyg(1,nseg_c+1),keyv(nseg_c+1), & 
     y_c,y_f,hpsi(1),hpsi(nvctr_c+1),nfl1,nfl2,nfl3,nfu1,nfu2,nfu3)

     i_all=-product(shape(y_c))*kind(y_c)
     deallocate(y_c,stat=i_stat)
     call memocc(i_stat,i_all,'y_c','applylocpotkinone_hyb')
     i_all=-product(shape(x_c))*kind(x_c)
     deallocate(x_c,stat=i_stat)
     call memocc(i_stat,i_all,'x_c','applylocpotkinone_hyb')

     i_all=-product(shape(x_f1))*kind(x_f1)
     deallocate(x_f1,stat=i_stat)
     call memocc(i_stat,i_all,'x_f1','applylocpotkinone_hyb')
     i_all=-product(shape(x_f2))*kind(x_f2)
     deallocate(x_f2,stat=i_stat)
     call memocc(i_stat,i_all,'x_f2','applylocpotkinone_hyb')
     i_all=-product(shape(x_f3))*kind(x_f3)
     deallocate(x_f3,stat=i_stat)
     call memocc(i_stat,i_all,'x_f3','applylocpotkinone_hyb')
     i_all=-product(shape(y_f))*kind(y_f)
     deallocate(y_f,stat=i_stat)
     call memocc(i_stat,i_all,'y_f','applylocpotkinone_hyb')
     i_all=-product(shape(x_f))*kind(x_f)
     deallocate(x_f,stat=i_stat)
     call memocc(i_stat,i_all,'x_f','applylocpotkinone_hyb')
     i_all=-product(shape(w))*kind(w)
     deallocate(w,stat=i_stat)
     call memocc(i_stat,i_all,'w','applylocpotkinone_hyb')
     i_all=-product(shape(ww))*kind(ww)
     deallocate(ww,stat=i_stat)
     call memocc(i_stat,i_all,'ww','applylocpotkinone_hyb')
 
END SUBROUTINE applylocpotkinone_hyb


!>  Applies the local potential and kinetic energy operator to one wavefunction 
!! Input: pot,psi
!! Output: hpsi,epot,ekin
subroutine applylocpotkinone_slab(n1,n2,n3, & 
     hx,hy,hz,nseg_c,nseg_f,nvctr_c,nvctr_f,keyg,keyv,  & 
     psir,psi_in,psi_out,psi,pot,hpsi,epot,ekin,nspinor,npot)
  use module_base
  use module_interfaces
  implicit none
  integer, intent(in) :: n1,n2,n3,nseg_c,nseg_f,nvctr_c,nvctr_f,nspinor,npot
  real(gp), intent(in) :: hx,hy,hz
  integer, dimension(nseg_c+nseg_f), intent(in) :: keyv
  integer, dimension(2,nseg_c+nseg_f), intent(in) :: keyg
  real(wp), dimension((2*n1+2)*(2*n2+31)*(2*n3+2),npot), intent(in) :: pot
  real(wp), dimension(nvctr_c+7*nvctr_f,nspinor), intent(in) :: psi
  real(wp), dimension((2*n1+2)*(2*n2+16)*(2*n3+2),nspinor), intent(inout) :: psi_in
  real(wp), dimension((2*n1+2)*(2*n2+31)*(2*n3+2),nspinor), intent(inout) :: psir,psi_out
  real(gp), intent(out) :: epot,ekin
  real(wp), dimension(nvctr_c+7*nvctr_f,nspinor), intent(out) :: hpsi
  !local variables
  integer :: i,idx
  real(wp) :: tt
  real(gp) :: v,p
  real(gp), dimension(3) :: hgridh

! Wavefunction expressed everywhere in fine scaling functions (for potential and kinetic energy)
!	psir serves as a work array	   

  hgridh(1)=hx*.5_gp
  hgridh(2)=hy*.5_gp
  hgridh(3)=hz*.5_gp

  do idx=1,nspinor
     call uncompress_slab(n1,n2,n3,nseg_c,nvctr_c,keyg(1,1),keyv(1),   &
          nseg_f,nvctr_f,keyg(1,nseg_c+1),keyv(nseg_c+1),   &
          psi(1,idx),psi(nvctr_c+1,idx),psi_in(1,idx),psir(1,idx))
     
     !	psi_out serves as a work array	   
     call convolut_magic_n_slab(2*n1+1,2*n2+15,2*n3+1,psi_in(1,idx),psir(1,idx),psi_out(1,idx)) 
  end do
  
  call apply_potential(n1,n2,n3,0,1,0,0,nspinor,npot,psir,pot,epot)

!!  epot=0.0_gp
!!  do i=1,(2*n1+2)*(2*n2+31)*(2*n3+2)
!!     v=real(pot(i),gp)
!!     p=real(psir(i),gp)
!!     tt=pot(i)*psir(i)
!!     epot=epot+p*v*p
!!     psir(i)=tt
!!  enddo

  do idx=1,nspinor
     call convolut_magic_t_slab_self(2*n1+1,2*n2+15,2*n3+1,psir(1,idx),psi_out(1,idx))

     ! compute the kinetic part and add  it to psi_out
     ! the kinetic energy is calculated at the same time
     call convolut_kinetic_slab_T(2*n1+1,2*n2+15,2*n3+1,hgridh,psi_in(1,idx),psi_out(1,idx),ekin)
  
     call compress_slab(n1,n2,n3,nseg_c,nvctr_c,keyg(1,1),keyv(1),   & 
          nseg_f,nvctr_f,keyg(1,nseg_c+1),keyv(nseg_c+1),   & 
          psi_out(1,idx),hpsi(1,idx),hpsi(nvctr_c+1,idx),psir(1,idx))
  end do

END SUBROUTINE applylocpotkinone_slab


!>  Applies the local potential and kinetic energy operator to one wavefunction 
!! Input: pot,psi
!! Output: hpsi,epot,ekin
subroutine applylocpotkinone_prev(n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,nbuf, & 
     hgrid,nseg_c,nseg_f,nvctr_c,nvctr_f,keyg,keyv,  & 
     ibyz_c,ibxz_c,ibxy_c,ibyz_f,ibxz_f,ibxy_f, & 
     y_c,y_f,psir,  &
     psi,pot,hpsi,epot,ekin,x_c,x_fc,x_f,w1,w2,&
     ibzzx_c,ibyyzz_c,ibxy_ff,ibzzx_f,ibyyzz_f,&
     ibzxx_c,ibxxyy_c,ibyz_ff,ibzxx_f,ibxxyy_f,nw1,nw2,ibyyzz_r)
  use module_base
  implicit real(kind=8) (a-h,o-z)
  dimension ibyz_c(2,0:n2,0:n3),ibxz_c(2,0:n1,0:n3),ibxy_c(2,0:n1,0:n2)
  dimension ibyz_f(2,0:n2,0:n3),ibxz_f(2,0:n1,0:n3),ibxy_f(2,0:n1,0:n2)
  dimension pot((2*n1+31)*(2*n2+31)*(2*n3+31))
  dimension keyg(2,nseg_c+nseg_f),keyv(nseg_c+nseg_f)
  dimension psi(nvctr_c+7*nvctr_f),scal(0:3)
  dimension hpsi(nvctr_c+7*nvctr_f)
  dimension y_c(0:n1,0:n2,0:n3)
  dimension y_f(7,0:n1,0:n2,0:n3)
  dimension psir((2*n1+31)*(2*n2+31)*(2*n3+31))
  !********************Alexey***************************************************************
  ! for shrink:
  integer ibzzx_c(2,-14:2*n3+16,0:n1) 
  integer ibyyzz_c(2,-14:2*n2+16,-14:2*n3+16)

  integer ibxy_ff(2,nfl1:nfu1,nfl2:nfu2)
  integer ibzzx_f(2,-14+2*nfl3:2*nfu3+16,nfl1:nfu1) 
  integer ibyyzz_f(2,-14+2*nfl2:2*nfu2+16,-14+2*nfl3:2*nfu3+16)

  ! for grow:
  integer ibzxx_c(2,0:n3,-14:2*n1+16) ! extended boundary arrays
  integer ibxxyy_c(2,-14:2*n1+16,-14:2*n2+16)

  integer ibyz_ff(2,nfl2:nfu2,nfl3:nfu3)
  integer ibzxx_f(2,nfl3:nfu3,2*nfl1-14:2*nfu1+16)
  integer ibxxyy_f(2,2*nfl1-14:2*nfu1+16,2*nfl2-14:2*nfu2+16)

  ! for real space:
  integer,intent(in):: ibyyzz_r(2,-14:2*n2+16,-14:2*n3+16)
  !*****************************************************************************************
  real(kind=8) x_c(0:n1,0:n2,0:n3),  x_fc(0:n1,0:n2,0:n3,3), x_f(7,nfl1:nfu1,nfl2:nfu2,nfl3:nfu3)! input
  real(kind=8) w1(nw1),w2(nw2) ! work
  !***********************************************************************************************
  character(len=*), parameter :: subname='applylocpotkinone_prev'
  do i=0,3
     scal(i)=1.d0
  enddo

  call uncompress_forstandard_prev(n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,  & 
       nseg_c,nvctr_c,keyg(1,1),keyv(1),  & 
       nseg_f,nvctr_f,keyg(1,nseg_c+1),keyv(nseg_c+1),   &
       scal,psi(1),psi(nvctr_c+1),x_c,x_fc,x_f)

  call comb_grow_all_prev(n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,w1,w2,x_c,x_f,  & 
       psir,ibyz_c,ibzxx_c,ibxxyy_c,ibyz_ff,ibzxx_f,ibxxyy_f)

  if (nbuf.eq.0) then
     call realspace_prev(ibyyzz_r,pot,psir,epot,n1,n2,n3)
  else
     call realspace_nbuf_prev(ibyyzz_r,pot,psir,epot,n1,n2,n3,nbuf)
  endif

  call comb_shrink(n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,w1,w2,psir,&
       ibxy_c,ibzzx_c,ibyyzz_c,ibxy_ff,ibzzx_f,ibyyzz_f,y_c,y_f)!,ibyz_c,ibyz_f)

  call ConvolkineticT_prev(n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,  &
       hgrid,ibyz_c,ibxz_c,ibxy_c,ibyz_f,ibxz_f,ibxy_f,x_c,x_fc,x_f,y_c,y_f,ekin)

  call compress_forstandard(n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,  &
       nseg_c,nvctr_c,keyg(1,1),       keyv(1),   &
       nseg_f,nvctr_f,keyg(1,nseg_c+1),keyv(nseg_c+1),   &
       scal,y_c,y_f,hpsi(1),hpsi(nvctr_c+1))

END SUBROUTINE applylocpotkinone_prev


subroutine realspace_prev(ibyyzz_r,pot,psir,epot,n1,n2,n3)
  implicit none
  integer,intent(in)::n1,n2,n3
  integer,intent(in)::ibyyzz_r(2,-14:2*n2+16,-14:2*n3+16)

  real(kind=8),intent(in)::pot(-14:2*n1+16,-14:2*n2+16,-14:2*n3+16)
  real(kind=8),intent(inout)::psir(-14:2*n1+16,-14:2*n2+16,-14:2*n3+16)

  real(kind=8),intent(out)::epot
  real(kind=8) tt
  integer i1,i2,i3

  epot=0.d0
  do i3=-14,2*n3+16
     do i2=-14,2*n2+16
        do i1=max(ibyyzz_r(1,i2,i3),-14),min(ibyyzz_r(2,i2,i3),2*n1+16)
           tt=pot(i1,i2,i3)*psir(i1,i2,i3)
           epot=epot+tt*psir(i1,i2,i3)
           psir(i1,i2,i3)=tt
        enddo
     enddo
  enddo

END SUBROUTINE realspace_prev


subroutine realspace_nbuf_prev(ibyyzz_r,pot,psir,epot,nb1,nb2,nb3,nbuf)
  implicit none
  integer,intent(in)::nb1,nb2,nb3,nbuf
  integer,intent(in)::ibyyzz_r(2,-14:2*nb2+16,-14:2*nb3+16)
  real(kind=8),intent(in)::pot(-14:2*nb1+16-4*nbuf,-14:2*nb2+16-4*nbuf,-14:2*nb3+16-4*nbuf)
  real(kind=8),intent(inout)::psir(-14:2*nb1+16,-14:2*nb2+16,-14:2*nb3+16)

  real(kind=8),intent(out)::epot
  real(kind=8) tt,dnrm2
  integer i1,i2,i3

  epot=0.d0
  do i3=-14,2*nb3+16
     if (i3.ge.-14+2*nbuf .and. i3.le.2*nb3+16-2*nbuf) then
        do i2=-14,2*nb2+16
           if (i2.ge.-14+2*nbuf .and. i2.le.2*nb2+16-2*nbuf) then
              do i1=-14+2*nbuf,ibyyzz_r(1,i2,i3)-1
                 psir(i1,i2,i3)=0.d0
              enddo
              do i1=max(ibyyzz_r(1,i2,i3),-14+2*nbuf),min(ibyyzz_r(2,i2,i3),2*nb1+16-2*nbuf)
                 tt=pot(i1-2*nbuf,i2-2*nbuf,i3-2*nbuf)*psir(i1,i2,i3)
                 epot=epot+tt*psir(i1,i2,i3)
                 psir(i1,i2,i3)=tt
              enddo
              do i1=ibyyzz_r(2,i2,i3)+1,2*nb1+16-2*nbuf
                 psir(i1,i2,i3)=0.d0
              enddo
           else
              do i1=-14,2*nb1+16
                 psir(i1,i2,i3)=0.d0
              enddo
           endif
        enddo
     else
        do i2=-14,2*nb2+16
           do i1=-14,2*nb1+16
              psir(i1,i2,i3)=0.d0
           enddo
        enddo
     endif
  enddo

END SUBROUTINE realspace_nbuf_prev


!> ypsi = @f$(1/2) \nabla^2 xpsi@f$
SUBROUTINE CALC_GRAD_REZA_prev(n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3, & 
     nseg_c,nvctr_c,keyg_c,keyv_c,nseg_f,nvctr_f,keyg_f,keyv_f, &
     scal,cprecr,hgrid,ibyz_c,ibxz_c,ibxy_c,ibyz_f,ibxz_f,ibxy_f,xpsi_c,xpsi_f,ypsi_c,ypsi_f)
  use module_base
  implicit none
  integer :: keyg_c(2,nseg_c),keyv_c(nseg_c),keyg_f(2,nseg_f),keyv_f(nseg_f)
  real(kind=8) :: xpsi_c(nvctr_c),xpsi_f(7,nvctr_f),scal(0:3)
  real(kind=8) :: ypsi_c(nvctr_c),ypsi_f(7,nvctr_f)
  integer :: ibyz_c(2,0:n2,0:n3),ibxz_c(2,0:n1,0:n3),ibxy_c(2,0:n1,0:n2)
  integer :: ibyz_f(2,0:n2,0:n3),ibxz_f(2,0:n1,0:n3),ibxy_f(2,0:n1,0:n2)
  character(len=*), parameter :: subname='CALC_GRAD_REZA_prev'
  real(kind=8),allocatable :: xpsig_c(:,:,:),ypsig_c(:,:,:)
  real(kind=8),allocatable :: xpsig_f(:,:,:,:),ypsig_f(:,:,:,:)
  real(kind=8),allocatable :: xpsig_fc(:,:,:,:)
  real(kind=8) :: cprecr,hgrid
  integer :: i_all,i_stat
  integer :: n1,n2,n3,nfl1,nfl2,nfl3,nfu1,nfu2,nfu3,nseg_c,nseg_f,nvctr_c,nvctr_f

  allocate(xpsig_c(0:n1,0:n2,0:n3+ndebug),stat=i_stat)
  call memocc(i_stat,xpsig_c,'xpsig_c',subname)
  allocate(xpsig_fc(0:n1,0:n2,0:n3,3+ndebug),stat=i_stat)
  call memocc(i_stat,xpsig_fc,'xpsig_fc',subname)
  allocate(xpsig_f(7,nfl1:nfu1,nfl2:nfu2,nfl3:nfu3+ndebug),stat=i_stat)
  call memocc(i_stat,xpsig_f,'xpsig_f',subname)
  allocate(ypsig_c(0:n1,0:n2,0:n3+ndebug),stat=i_stat)
  call memocc(i_stat,ypsig_c,'ypsig_c',subname)
  allocate(ypsig_f(7,nfl1:nfu1,nfl2:nfu2,nfl3:nfu3+ndebug),stat=i_stat)
  call memocc(i_stat,ypsig_f,'ypsig_f',subname)

  call uncompress_forstandard_prev(n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,  &
       nseg_c,nvctr_c,keyg_c,keyv_c,  & 
       nseg_f,nvctr_f,keyg_f,keyv_f,  & 
       scal,xpsi_c,xpsi_f,xpsig_c,xpsig_fc,xpsig_f)

  call Convolkinetic_prev(n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3, & 
       cprecr,hgrid,ibyz_c,ibxz_c,ibxy_c,ibyz_f,ibxz_f,ibxy_f,xpsig_c,xpsig_fc,xpsig_f,ypsig_c,ypsig_f)

  call compress_forstandard(n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,  &
       nseg_c,nvctr_c,keyg_c,keyv_c,  & 
       nseg_f,nvctr_f,keyg_f,keyv_f,  & 
       scal,ypsig_c,ypsig_f,ypsi_c,ypsi_f)

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
  i_all=-product(shape(xpsig_fc))*kind(xpsig_fc)
  deallocate(xpsig_fc,stat=i_stat)
  call memocc(i_stat,i_all,'xpsig_fc',subname)

END SUBROUTINE CALC_GRAD_REZA_prev
