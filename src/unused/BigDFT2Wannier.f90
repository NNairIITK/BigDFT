!!>
!subroutine angularpart(l, mr, np, nx, ny, nz, ix, iy, iz, &
!                    xx, yy, zz, n_proj, ylm)
!
!   ! This routine returns the angular part of the spherical harmonic identified by indices (l,mr)
!   ! Calcutations are made in Cartesian coordinates
!
!   implicit none
!
!   ! I/O variables
!   integer, intent(in) :: l(n_proj), mr(n_proj)
!   integer, intent(in) :: np, nx, ny, nz, ix, iy, iz, n_proj
!   real(kind=8), intent(in) :: xx, yy, zz
!   real(kind=8), dimension(nx,ny,nz), intent(out) :: ylm
!
!   ! local variables
!   real(kind=8), parameter :: pi=3.141592653589793238462643383279d0
!   real(kind=8), parameter :: eps8  = 1.0e-8
!   real(kind=8), external :: s, pz, px, py, dz2, dxz, dyz, dx2my2, dxy, &
!                          fz3, fxz2, fyz2, fzx2my2, fxyz, fxx2m3y2, fy3x2my2
!   real(kind=8) :: rr
!   real(kind=8) :: bs2, bs3, bs6, bs12
!
!   bs2 = 1.d0/sqrt(2.d0)
!   bs3 = 1.d0/sqrt(3.d0)
!   bs6 = 1.d0/sqrt(6.d0)
!   bs12 = 1.d0/sqrt(12.d0)
!
!
!   if (l(np) > 3 .OR. l(np) < -5 ) then 
!      write(*,*) 'error, l out of range '
!   else
!      if (l(np)>=0) then
!         if (mr(np) < 1 .OR. mr(np) > 2*l(np)+1) then
!	    write(*,*) 'error, mr out of range'
!	 end if
!      else
!         if (mr(np) < 1 .OR. mr(np) > abs(l(np))+1 ) then 
!	    write(*,*) 'error, mr out of range'
!         end if
!      end if
!   end if
!
!   rr = sqrt( xx*xx + yy*yy + zz*zz )
!    
!   if (l(np)==0) then   ! s orbital
!      ylm(ix,iy,iz) = s(xx,yy,zz,rr)  
!   end if
!
!   if (l(np)==1) then   ! p orbitals
!      if (mr(np)==1) ylm(ix,iy,iz) = pz(xx,yy,zz,rr) 
!      if (mr(np)==2) ylm(ix,iy,iz) = px(xx,yy,zz,rr)
!      if (mr(np)==3) ylm(ix,iy,iz) = py(xx,yy,zz,rr)
!   end if
!
!   if (l(np)==2) then   ! d orbitals
!      if (mr(np)==1) ylm(ix,iy,iz) = dz2(xx,yy,zz,rr)
!      if (mr(np)==2) ylm(ix,iy,iz) = dxz(xx,yy,zz,rr)
!      if (mr(np)==3) ylm(ix,iy,iz) = dyz(xx,yy,zz,rr)
!      if (mr(np)==4) ylm(ix,iy,iz) = dx2my2(xx,yy,zz,rr)
!      if (mr(np)==5) ylm(ix,iy,iz) = dxy(xx,yy,zz,rr)
!   endif
!
!   if (l(np)==3) then   ! f orbitals
!      if (mr(np)==1) ylm(ix,iy,iz) = fz3(xx,yy,zz,rr)
!      if (mr(np)==2) ylm(ix,iy,iz) = fxz2(xx,yy,zz,rr)
!      if (mr(np)==3) ylm(ix,iy,iz) = fyz2(xx,yy,zz,rr)
!      if (mr(np)==4) ylm(ix,iy,iz) = fzx2my2(xx,yy,zz,rr)
!      if (mr(np)==5) ylm(ix,iy,iz) = fxyz(xx,yy,zz,rr)
!      if (mr(np)==6) ylm(ix,iy,iz) = fxx2m3y2(xx,yy,zz,rr)
!      if (mr(np)==7) ylm(ix,iy,iz) = fy3x2my2(xx,yy,zz,rr)
!   endif
!
!   if (l(np)==-1) then  !  sp hybrids
!      if (mr(np)==1) ylm(ix,iy,iz) = bs2 * ( s(xx,yy,zz,rr) + px(xx,yy,zz,rr) ) 
!      if (mr(np)==2) ylm(ix,iy,iz) = bs2 * ( s(xx,yy,zz,rr) - px(xx,yy,zz,rr) ) 
!   end if
!
!   if (l(np)==-2) then  !  sp2 hybrids 
!      if (mr(np)==1) ylm(ix,iy,iz) = bs3*s(xx,yy,zz,rr)-bs6*px(xx,yy,zz,rr)+bs2*py(xx,yy,zz,rr)
!      if (mr(np)==2) ylm(ix,iy,iz) = bs3*s(xx,yy,zz,rr)-bs6*px(xx,yy,zz,rr)-bs2*py(xx,yy,zz,rr)
!      if (mr(np)==3) ylm(ix,iy,iz) = bs3*s(xx,yy,zz,rr) +2.d0*bs6*px(xx,yy,zz,rr) 
!   end if
!
!   if (l(np)==-3) then  !  sp3 hybrids
!      if (mr(np)==1) ylm(ix,iy,iz) = 0.5d0*(s(xx,yy,zz,rr)+px(xx,yy,zz,rr)+py(xx,yy,zz,rr)+pz(xx,yy,zz,rr))
!      if (mr(np)==2) ylm(ix,iy,iz) = 0.5d0*(s(xx,yy,zz,rr)+px(xx,yy,zz,rr)-py(xx,yy,zz,rr)-pz(xx,yy,zz,rr))
!      if (mr(np)==3) ylm(ix,iy,iz) = 0.5d0*(s(xx,yy,zz,rr)-px(xx,yy,zz,rr)+py(xx,yy,zz,rr)-pz(xx,yy,zz,rr))
!      if (mr(np)==4) ylm(ix,iy,iz) = 0.5d0*(s(xx,yy,zz,rr)-px(xx,yy,zz,rr)-py(xx,yy,zz,rr)+pz(xx,yy,zz,rr))
!   end if
!
!   if (l(np)==-4) then  !  sp3d hybrids
!      if (mr(np)==1) ylm(ix,iy,iz) = bs3*s(xx,yy,zz,rr)-bs6*px(xx,yy,zz,rr)+bs2*py(xx,yy,zz,rr)
!      if (mr(np)==2) ylm(ix,iy,iz) = bs3*s(xx,yy,zz,rr)-bs6*px(xx,yy,zz,rr)-bs2*py(xx,yy,zz,rr)
!      if (mr(np)==3) ylm(ix,iy,iz) = bs3*s(xx,yy,zz,rr) +2.d0*bs6*px(xx,yy,zz,rr) 
!      if (mr(np)==4) ylm(ix,iy,iz) = bs2*pz(xx,yy,zz,rr)+bs2*dz2(xx,yy,zz,rr)
!      if (mr(np)==5) ylm(ix,iy,iz) =-bs2*pz(xx,yy,zz,rr)+bs2*dz2(xx,yy,zz,rr)
!   end if
!
!   if (l(np)==-5) then  ! sp3d2 hybrids
!      if (mr(np)==1) ylm(ix,iy,iz) = bs6*s(xx,yy,zz,rr)-bs2*px(xx,yy,zz,rr)-bs12*dz2(xx,yy,zz,rr)+.5d0*dx2my2(xx,yy,zz,rr)
!      if (mr(np)==2) ylm(ix,iy,iz) = bs6*s(xx,yy,zz,rr)+bs2*px(xx,yy,zz,rr)-bs12*dz2(xx,yy,zz,rr)+.5d0*dx2my2(xx,yy,zz,rr)
!      if (mr(np)==3) ylm(ix,iy,iz) = bs6*s(xx,yy,zz,rr)-bs2*py(xx,yy,zz,rr)-bs12*dz2(xx,yy,zz,rr)-.5d0*dx2my2(xx,yy,zz,rr)
!      if (mr(np)==4) ylm(ix,iy,iz) = bs6*s(xx,yy,zz,rr)+bs2*py(xx,yy,zz,rr)-bs12*dz2(xx,yy,zz,rr)-.5d0*dx2my2(xx,yy,zz,rr)
!      if (mr(np)==5) ylm(ix,iy,iz) = bs6*s(xx,yy,zz,rr)-bs2*pz(xx,yy,zz,rr)+bs3*dz2(xx,yy,zz,rr)
!      if (mr(np)==6) ylm(ix,iy,iz) = bs6*s(xx,yy,zz,rr)+bs2*pz(xx,yy,zz,rr)+bs3*dz2(xx,yy,zz,rr)
!   end if
!
!END SUBROUTINE angularpart
!
!
!! The following functions are used to calculate angular parts of the spherical harmonics
!!======== l = 0 =====================================================================
!function s(xx,yy,zz,rr)
!   implicit none
!   real(kind=8), parameter :: pi=3.141592653589793238462643383279d0
!   real(kind=8) :: s, xx, yy, zz, rr
!   s = 1.d0/ sqrt(4*pi)
!END FUNCTION s
!
!
!!======== l = 1 =====================================================================
!function pz(xx,yy,zz,rr)
!   implicit none
!   real(kind=8), parameter :: pi=3.141592653589793238462643383279d0
!   real(kind=8) ::pz, xx, yy, zz, rr
!   pz =  sqrt(3.d0/(4*pi)) * (zz/rr)
!END FUNCTION pz
!
!function px(xx,yy,zz,rr)
!   implicit none
!   real(kind=8), parameter :: pi=3.141592653589793238462643383279d0
!   real(kind=8) ::px, xx, yy, zz, rr
!   px =  sqrt(3.d0/(4*pi)) * (xx/rr)
!END FUNCTION px
!
!function py(xx,yy,zz,rr)
!   implicit none
!   real(kind=8), parameter :: pi=3.141592653589793238462643383279d0
!   real(kind=8) ::py, xx, yy, zz, rr
!   py =  sqrt(3.d0/(4*pi)) * (yy/rr)
!END FUNCTION py
!
!
!!======== l = 2 =====================================================================
!function dz2(xx,yy,zz,rr)
!   implicit none
!   real(kind=8), parameter :: pi=3.141592653589793238462643383279d0
!   real(kind=8) ::dz2, xx, yy, zz, rr
!   dz2 =  sqrt(1.25d0/(4*pi)) * (-xx*xx-yy*yy+2.d0*zz*zz)/(rr*rr)
!END FUNCTION dz2
!
!function dxz(xx,yy,zz,rr)
!   implicit none
!   real(kind=8), parameter :: pi=3.141592653589793238462643383279d0
!   real(kind=8) ::dxz, xx, yy, zz, rr
!   dxz =  sqrt(15.d0/(4*pi)) * (xx*zz)/(rr*rr)
!END FUNCTION dxz
!
!function dyz(xx,yy,zz,rr)
!   implicit none
!   real(kind=8), parameter :: pi=3.141592653589793238462643383279d0
!   real(kind=8) ::dyz, xx, yy, zz, rr
!   dyz =  sqrt(15.d0/(4*pi)) * (yy*zz)/(rr*rr)
!END FUNCTION dyz
!
!function dx2my2(xx,yy,zz,rr)
!   implicit none
!   real(kind=8), parameter :: pi=3.141592653589793238462643383279d0
!   real(kind=8) ::dx2my2, xx, yy, zz, rr
!   dx2my2 =  sqrt(3.75d0/(4*pi)) * (xx*xx-yy*yy)/(rr*rr)
!END FUNCTION dx2my2
!
!function dxy(xx,yy,zz,rr)
!   implicit none
!   real(kind=8), parameter :: pi=3.141592653589793238462643383279d0
!   real(kind=8) ::dxy, xx, yy, zz, rr
!   dxy =  sqrt(3.75d0/(4*pi)) * (xx*yy)/(rr*rr)
!END FUNCTION dxy
!
!
!!======== l = 3 =====================================================================
!function fz3(xx,yy,zz,rr)
!   implicit none
!   real(kind=8), parameter :: pi=3.141592653589793238462643383279d0
!   real(kind=8) ::fz3, xx, yy, zz, rr
!   fz3 =  0.25d0*sqrt(7.d0/pi) * (zz*(2.d0*zz*zz-3.d0*xx*xx-3.d0*yy*yy)) / (rr*rr*rr)
!END FUNCTION fz3
!
!function fxz2(xx,yy,zz,rr)
!   implicit none
!   real(kind=8), parameter :: pi=3.141592653589793238462643383279d0
!   real(kind=8) ::fxz2, xx, yy, zz, rr
!   fxz2 =  0.25d0*sqrt(10.5d0/pi) * (xx*(4.d0*zz*zz-xx*xx-yy*yy)) / (rr*rr*rr)
!END FUNCTION fxz2
!
!function fyz2(xx,yy,zz,rr)
!   implicit none
!   real(kind=8), parameter :: pi=3.141592653589793238462643383279d0
!   real(kind=8) ::fyz2, xx, yy, zz, rr
!   fyz2 =  0.25d0*sqrt(10.5d0/pi) * (yy*(4.d0*zz*zz-xx*xx-yy*yy)) / (rr*rr*rr)
!END FUNCTION fyz2
!
!function fzx2my2(xx,yy,zz,rr)
!   implicit none
!   real(kind=8), parameter :: pi=3.141592653589793238462643383279d0
!   real(kind=8) ::fzx2my2, xx, yy, zz, rr
!   fzx2my2 =  0.25d0*sqrt(105d0/pi) * (zz*(xx*xx-yy*yy)) / (rr*rr*rr)
!END FUNCTION fzx2my2
!
!function fxyz(xx,yy,zz,rr)
!   implicit none
!   real(kind=8), parameter :: pi=3.141592653589793238462643383279d0
!   real(kind=8) ::fxyz, xx, yy, zz, rr
!   fxyz =  0.25d0*sqrt(105d0/pi) * (xx*yy*zz) / (rr*rr*rr)
!END FUNCTION fxyz
!
!function fxx2m3y2(xx,yy,zz,rr)
!   implicit none
!   real(kind=8), parameter :: pi=3.141592653589793238462643383279d0
!   real(kind=8) ::fxx2m3y2, xx, yy, zz, rr
!   fxx2m3y2 =  0.25d0*sqrt(17.5d0/pi) * (xx*(xx*xx-3.d0*yy*yy)) / (rr*rr*rr)
!END FUNCTION fxx2m3y2
!
!function fy3x2my2(xx,yy,zz,rr)
!   implicit none
!   real(kind=8), parameter :: pi=3.141592653589793238462643383279d0
!   real(kind=8) ::fy3x2my2, xx, yy, zz, rr
!   fy3x2my2 =  0.25d0*sqrt(17.5d0/pi) * (yy*(3.d0*xx*xx-yy*yy)) / (rr*rr*rr)
!END FUNCTION fy3x2my2


!!subroutine make_precheck(iproc,nproc,input,Glr,orbsv,commsv,orbsp,commsp,atoms,w,rxyz,n_proj,ctr_proj,&
!!           x_proj,y_proj,z_proj,l,mr,rvalue,zona,amnk_bands_sorted,sph_daub)
!!   use BigDFT_API
!!   use Poisson_Solver
!!   implicit none
!!   integer, intent(in) :: iproc, nproc, n_proj
!!   type(input_variables),intent(in) :: input
!!   type(locreg_descriptors), intent(in) :: Glr
!!   type(orbitals_data), intent(inout) :: orbsv,orbsp
!!   type(communications_arrays), target :: commsv,commsp
!!   type(atoms_data), intent(in) :: atoms
!!   type(workarr_sumrho), intent(in) :: w
!!   real(gp), dimension(3,atoms%nat), intent(in) :: rxyz
!!   real(kind=8), dimension (n_proj,3), intent(in) :: ctr_proj, x_proj, y_proj, z_proj
!!   integer, dimension (n_proj),intent(in) :: l, mr, rvalue
!!   real, dimension (n_proj), intent(in) :: zona
!!   integer, dimension (:), pointer :: amnk_bands_sorted
!!   real(wp),dimension(:),pointer :: sph_daub
!!   !local variables
!!   character(len=*), parameter :: subname='make_precheck'
!!   integer :: i_stat, i_all, npsidim, i, j, k, np, npp, pshft
!!   integer :: ind, nb, ierr, npsidim2,nvctrp
!!   real(kind=8) :: b1, b2, b3, r0x, r0y, r0z, zz, yy, xx
!!   real(wp), allocatable :: psi_etsfv(:,:),sph_har_etsf(:)!,sph_daub(:)
!!   real(wp), allocatable :: psi_etsf2(:)
!!   real(wp), pointer :: pwork(:)
!!   character(len=60) :: filename
!!   real(kind=8), allocatable :: ylm(:,:,:), func_r(:,:,:)
!!   real(kind=8), allocatable :: amnk(:,:), amnk_guess(:)
!!   integer :: n_virt, n_virt_tot
!!   real(kind=8), allocatable, dimension(:) :: amnk_guess_sorted
!!   real(gp), dimension(3,atoms%nat) :: rxyz_old
!!
!!   call timing(iproc,'CrtProjectors ','ON')
!!
!!   ! Read wavefunction from file and transforms it properly if hgrid or size of simulation cell have changed
!!   allocate(psi_etsfv(Glr%wfd%nvctr_c+7*Glr%wfd%nvctr_f,max(orbsv%norbp*orbsv%nspinor,1)),stat=i_stat)
!!   call memocc(i_stat,psi_etsfv,'psi_etsfv',subname)
!!   if(associated(orbsv%eval)) nullify(orbsv%eval)
!!   allocate(orbsv%eval(orbsv%norb*orbsv%nkpts), stat=i_stat)
!!   call memocc(i_stat,orbsv%eval,'orbsv%eval',subname)
!!
!!   filename= trim(input%dir_output) // 'virtuals'// trim(wfformat_read)
!!   call readmywaves(iproc,filename,orbsv,Glr%d%n1,Glr%d%n2,Glr%d%n3,input%hx,input%hy,input%hz,atoms,rxyz_old,rxyz,  & 
!!      Glr%wfd,psi_etsfv)
!!   i_all = -product(shape(orbsv%eval))*kind(orbsv%eval)
!!   deallocate(orbsv%eval,stat=i_stat)
!!   nullify(orbsv%eval)
!!   call memocc(i_stat,i_all,'orbsv%eval',subname)
!!
!!   ! Tranposition of the distribution of the BigDFT wavefunctions : orbitals -> components.
!!   npsidim=max((Glr%wfd%nvctr_c+7*Glr%wfd%nvctr_f)*orbsv%norbp*orbsv%nspinor,sum(commsv%ncntt(0:nproc-1)))
!!   npsidim2=max((Glr%wfd%nvctr_c+7*Glr%wfd%nvctr_f)*orbsp%norbp,sum(commsp%ncntt(0:nproc-1)))
!!   allocate(psi_etsf2(npsidim),stat=i_stat) !!doing this because psi_etsfv does not incorporate enough space for transpose
!!   call memocc(i_stat,psi_etsf2,'psi_etsf2',subname)
!!
!!   call razero(npsidim,psi_etsf2)
!!   if(nproc > 1) then
!!     allocate(pwork(npsidim),stat=i_stat)
!!     call memocc(i_stat,pwork,'pwork',subname)
!!     call transpose_v(iproc,nproc,orbsv,Glr%wfd,commsv,psi_etsfv(1,1),work=pwork,outadd=psi_etsf2(1))
!!     i_all = -product(shape(pwork))*kind(pwork)
!!     deallocate(pwork,stat=i_stat)
!!     call memocc(i_stat,i_all,'pwork',subname)
!!   else
!!      ! just copy the wavefunctions 
!!      k=0
!!      do j=1,orbsv%norbp*orbsv%nspinor
!!      do i=1,Glr%wfd%nvctr_c+7*Glr%wfd%nvctr_f
!!         k=k+1
!!         psi_etsf2(k) = psi_etsfv(i,j)
!!      end do
!!      end do
!!   end if
!!   i_all=-product(shape(psi_etsfv))*kind(psi_etsfv)
!!   deallocate(psi_etsfv,stat=i_stat)
!!   call memocc(i_stat,i_all,'psi_etsfv',subname)
!!
!!   ! - b1, b2 and b3 are the norm of the lattice parameters.
!!   b1=atoms%alat1
!!   b2=atoms%alat2
!!   b3=atoms%alat3
!!   ! - Allocations
!!   allocate(amnk(orbsv%norb,orbsp%norb),stat=i_stat)
!!   call memocc(i_stat,amnk,'amnk',subname)
!!   allocate(amnk_guess(orbsv%norb),stat=i_stat)
!!   call memocc(i_stat,amnk_guess,'amnk_guess',subname)
!!   allocate(sph_daub(npsidim2), stat=i_stat)
!!   call memocc(i_stat,sph_daub,'sph_daub',subname)
!!
!!   ! Begining of the algorithm to compute the scalar product in order to find the best unoccupied orbitals to use to compute the actual Amnk matrix :
!!   if (iproc==0) then
!!      write(*,*) '!==================================!'
!!      write(*,*) '! Calculating amnk=<virt|sph_har>  !'
!!      write(*,*) '!       in pre-check mode :        !'
!!      write(*,*) '!==================================!'
!!      write(*,'(A12,4x,A15)') 'Virtual band', 'amnk_guess(nb)='
!!   end if
!!
!!   ! Calculation of the spherical harmonics in parallel.
!!   ! It is done in the real space and then converted in the Daubechies representation.
!!   pshft = 0
!!   do npp=1, orbsp%norbp
!!      np = npp + orbsp%isorb
!!      ! Convolution buffer : n1i=2*n1+31 -> explains the '13*input%hx*0.5' term
!!      r0x=ctr_proj(np,1)*b1+13*input%hx*0.5
!!      r0y=ctr_proj(np,2)*b2+13*input%hy*0.5
!!      r0z=ctr_proj(np,3)*b3+13*input%hz*0.5
!!      do k=1,Glr%d%n3i
!!         zz=(k-1)*input%hz*0.5-r0z
!!         do j=1,Glr%d%n2i
!!            yy=(j-1)*input%hy*0.5-r0y
!!            do i=1,Glr%d%n1i
!!               ind=(k-1)*Glr%d%n2i*Glr%d%n1i+(j-1)*Glr%d%n1i+i
!!               xx=(i-1)*input%hx*0.5-r0x
!!               call angularpart(l, mr, np, Glr%d%n1i, Glr%d%n2i, Glr%d%n3i, i, j, k, &
!!                     xx, yy, zz, x_proj, y_proj, z_proj, n_proj, ylm)
!!               call radialpart(rvalue, zona, np, Glr%d%n1i, Glr%d%n2i, Glr%d%n3i, i, j, k, &
!!                     xx, yy, zz, n_proj, func_r)
!!               ! The 'sqrt(input%hx*0.5*input%hy*0.5*input%hz*0.5)' term is here to normalize spherical harmonics
!!               sph_har_etsf(ind)=func_r(i,j,k)*ylm(i,j,k)*sqrt(input%hx*0.5*input%hy*0.5*input%hz*0.5)
!!            end do
!!         end do
!!      end do
!!      call isf_to_daub(Glr,w,sph_har_etsf(1),sph_daub(1+pshft))
!!      pshft=pshft + max(Glr%wfd%nvctr_c+7*Glr%wfd%nvctr_f,commsp%ncntt(iproc)/orbsp%norbp)
!!   end do
!!
!!   call timing(iproc,'CrtProjectors ','OF')
!!
!!   ! Tranposition of the distribution of the spherical harmonics: orbitals -> components.
!!   allocate(pwork(npsidim2),stat=i_stat)
!!   call memocc(i_stat,pwork,'pwork',subname)
!!   call transpose_v(iproc,nproc,orbsp,Glr%wfd,commsp,sph_daub,work=pwork)
!!   i_all = -product(shape(pwork))*kind(pwork)
!!   deallocate(pwork,stat=i_stat)
!!   call memocc(i_stat,i_all,'pwork',subname)
!!   call timing(iproc,'ApplyProj     ','ON')
!!
!!   ! Scalar product of amnk=<sph_daub|psi> in parallel.
!!   call razero(orbsp%norb*orbsv%norb,amnk)
!!   nvctrp=commsv%nvctr_par(iproc,1)
!!   call gemm('T','N',orbsv%norb,orbsp%norb,nvctrp,1.0_wp,psi_etsf2(1),max(1,nvctrp),&
!!        sph_daub(1),max(1,nvctrp),0.0_wp,amnk(1,1),orbsv%norb)
!!      
!!   ! Construction of the whole Amnk_guess matrix.
!!   call mpiallred(amnk(1,1),orbsv%norb*orbsp%norb,MPI_SUM,MPI_COMM_WORLD,ierr)
!!
!!   ! For each unoccupied orbitals, check how they project on spherical harmonics.
!!   ! The greater amnk_guess(nb) is, the more they project on spherical harmonics.
!!   do nb=1,orbsv%norb
!!      amnk_guess(nb)=0.0
!!      do np=1,orbsp%norb
!!         amnk_guess(nb)=amnk_guess(nb)+(amnk(nb,np))**2
!!      end do
!!      if (iproc==0) write(*,'(I4,11x,F12.6)') nb, sqrt(amnk_guess(nb))
!!   end do
!!
!!   ! Choice of the unoccupied orbitals to calculate the Amnk matrix
!!   if (iproc==0) then
!!      write(*,*) 
!!      write(*,'(1a)') 'These are the virtual bands to use to construct the actual Amn and Mmn matrices :'
!!      write(*,'(A4,4x,A17)') 'Band', 'sqrt(amnk_guess)='
!!   end if
!!   allocate(amnk_guess_sorted(n_virt),stat=i_stat)
!!   call memocc(i_stat,amnk_guess_sorted,'amnk_guess_sorted',subname)
!!   do nb=1,n_virt
!!      amnk_guess_sorted(nb)=maxval(amnk_guess,1)
!!      amnk_bands_sorted(nb)=maxloc(amnk_guess,1)
!!      amnk_guess(amnk_bands_sorted(nb))=0.d0
!!   if (iproc==0) write(*,'(I4,3x,F12.6)') amnk_bands_sorted(nb), sqrt(amnk_guess_sorted(nb))
!!   end do
!!
!!   ! End of the pre-check mode
!!   i_all = -product(shape(psi_etsf2))*kind(psi_etsf2)
!!   deallocate(psi_etsf2,stat=i_stat)
!!   call memocc(i_stat,i_all,'psi_etsf2',subname)
!!   i_all = -product(shape(amnk))*kind(amnk)
!!   deallocate(amnk,stat=i_stat)
!!   call memocc(i_stat,i_all,'amnk',subname)
!!   i_all = -product(shape(amnk_guess_sorted))*kind(amnk_guess_sorted)
!!   deallocate(amnk_guess_sorted,stat=i_stat)
!!   call memocc(i_stat,i_all,'amnk_guess_sorted',subname)
!!   i_all = -product(shape(amnk_guess))*kind(amnk_guess)
!!   deallocate(amnk_guess,stat=i_stat)
!!   call memocc(i_stat,i_all,'amnk_guess',subname)
!!
!!   if (iproc==0) then
!!      write(*,*) '!==================================!'
!!      write(*,*) '! Calculating amnk=<virt|sph_har>  !'
!!      write(*,*) '!     in pre-check mode done       !'
!!      write(*,*) '!==================================!'
!!      write(*,*)
!!      write(*,*)
!!   end if
!!
!!   ! Rewrite the input.inter file to add the chosen unoccupied states.
!!   if (iproc==0) call write_inter(n_virt, amnk_bands_sorted)
!!
!!   call timing(iproc,'ApplyProj     ','OF')
!!
!!END SUBROUTINE make_precheck
