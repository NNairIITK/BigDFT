
subroutine get_derivative_supportfunctions(ndim, hgrid, lzd, lorbs, phi, phid)
  use module_base
  use module_types
  use module_interfaces, except_this_one => get_derivative_supportfunctions
  implicit none
  
  ! Calling arguments
  integer,intent(in):: ndim
  real(kind=8),intent(in) :: hgrid
  type(local_zone_descriptors),intent(in) :: lzd
  type(orbitals_data),intent(in) :: lorbs
  real(kind=8),dimension(ndim),intent(in) :: phi !< Basis functions
  real(kind=8),dimension(3*ndim),intent(inout) :: phid  !< Derivative basis functions
  
  ! Local variables
  integer :: ist1, iorb, ist, ilr, iiorb

  ist=1
  ist1=1
  do iorb=1,lorbs%norbp
     iiorb=lorbs%isorb+iorb
     ilr=lorbs%inWhichLocreg(iiorb)

     call get_one_derivative_supportfunction(ilr,hgrid,lzd,phi(ist),phid(ist1))

     ist = ist + lzd%llr(ilr)%wfd%nvctr_c + 7*lzd%llr(ilr)%wfd%nvctr_f
     ist1 = ist1 + 3*(lzd%llr(ilr)%wfd%nvctr_c + 7*lzd%llr(ilr)%wfd%nvctr_f)
  end do

end subroutine get_derivative_supportfunctions


subroutine get_one_derivative_supportfunction(ilr,hgrid,lzd,phi,phid)
   use module_base
   use module_types
   use module_interfaces
   implicit none
   
   ! Calling arguments
   integer, intent(in) :: ilr
   real(kind=8),intent(in) :: hgrid
   type(local_zone_descriptors),intent(in) :: lzd
   real(kind=8),dimension(lzd%llr(ilr)%wfd%nvctr_c+7*lzd%llr(ilr)%wfd%nvctr_f),intent(in) :: phi !< Basis functions
   real(kind=8),dimension(3*(lzd%llr(ilr)%wfd%nvctr_c+7*lzd%llr(ilr)%wfd%nvctr_f)),intent(inout) :: phid  !< Derivative basis functions

   ! Local variables
   integer :: nf, istat, iall
   integer :: isty_c, istz_c, isty_f, istz_f
   real(kind=8),dimension(0:3),parameter :: scal=1.d0
   real(kind=8),dimension(:),allocatable :: w_f1, w_f2, w_f3
   real(kind=8),dimension(:,:,:),allocatable :: w_c, phix_c, phiy_c, phiz_c
   real(kind=8),dimension(:,:,:,:),allocatable :: w_f, phix_f, phiy_f, phiz_f
   character(len=*),parameter :: subname='get_one_derivative_supportfunction'

   call allocateWorkarrays()

   ! Uncompress the wavefunction.
   call uncompress_forstandard(lzd%llr(ilr)%d%n1, lzd%llr(ilr)%d%n2, lzd%llr(ilr)%d%n3, &
        lzd%llr(ilr)%d%nfl1, lzd%llr(ilr)%d%nfu1, & 
        lzd%llr(ilr)%d%nfl2, lzd%llr(ilr)%d%nfu2, lzd%llr(ilr)%d%nfl3, lzd%llr(ilr)%d%nfu3,  &
        lzd%llr(ilr)%wfd%nseg_c, lzd%llr(ilr)%wfd%nvctr_c, lzd%llr(ilr)%wfd%keygloc, lzd%llr(ilr)%wfd%keyvloc,  &
        lzd%llr(ilr)%wfd%nseg_f, lzd%llr(ilr)%wfd%nvctr_f, &
        lzd%llr(ilr)%wfd%keygloc(1,lzd%llr(ilr)%wfd%nseg_c+min(1,lzd%llr(ilr)%wfd%nseg_f)), &
        lzd%llr(ilr)%wfd%keyvloc(lzd%llr(ilr)%wfd%nseg_c+min(1,lzd%llr(ilr)%wfd%nseg_f)),  &
        scal, phi(1), phi(1+lzd%llr(ilr)%wfd%nvctr_c), w_c, w_f, w_f1, w_f2, w_f3)


   call createDerivativeBasis(lzd%llr(ilr)%d%n1, lzd%llr(ilr)%d%n2, lzd%llr(ilr)%d%n3, &
        lzd%llr(ilr)%d%nfl1, lzd%llr(ilr)%d%nfu1, lzd%llr(ilr)%d%nfl2, lzd%llr(ilr)%d%nfu2, &
        lzd%llr(ilr)%d%nfl3, lzd%llr(ilr)%d%nfu3,  &
        hgrid, lzd%llr(ilr)%bounds%kb%ibyz_c, lzd%llr(ilr)%bounds%kb%ibxz_c, lzd%llr(ilr)%bounds%kb%ibxy_c, &
        lzd%llr(ilr)%bounds%kb%ibyz_f, lzd%llr(ilr)%bounds%kb%ibxz_f, lzd%llr(ilr)%bounds%kb%ibxy_f, &
        w_c, w_f, w_f1, w_f2, w_f3, phix_c, phix_f, phiy_c, phiy_f, phiz_c, phiz_f)

   ! Compress the x wavefunction.
   call compress_forstandard(lzd%llr(ilr)%d%n1, lzd%llr(ilr)%d%n2, lzd%llr(ilr)%d%n3, &
        lzd%llr(ilr)%d%nfl1, lzd%llr(ilr)%d%nfu1, &
        lzd%llr(ilr)%d%nfl2, lzd%llr(ilr)%d%nfu2, lzd%llr(ilr)%d%nfl3, lzd%llr(ilr)%d%nfu3, &
        lzd%llr(ilr)%wfd%nseg_c, lzd%llr(ilr)%wfd%nvctr_c, lzd%llr(ilr)%wfd%keygloc, lzd%llr(ilr)%wfd%keyvloc, &
        lzd%llr(ilr)%wfd%nseg_f, lzd%llr(ilr)%wfd%nvctr_f, &
        lzd%llr(ilr)%wfd%keygloc(1,lzd%llr(ilr)%wfd%nseg_c+min(1,lzd%llr(ilr)%wfd%nseg_f)), &
        lzd%llr(ilr)%wfd%keyvloc(lzd%llr(ilr)%wfd%nseg_c+min(1,lzd%llr(ilr)%wfd%nseg_f)),  &
        scal, phix_c, phix_f, phid(1), phid(1+lzd%llr(ilr)%wfd%nvctr_c))

   ! Compress the y wavefunction.
   isty_c = 1 + lzd%llr(ilr)%wfd%nvctr_c + 7*lzd%llr(ilr)%wfd%nvctr_f
   isty_f = 1 + 2*lzd%llr(ilr)%wfd%nvctr_c + 7*lzd%llr(ilr)%wfd%nvctr_f
   call compress_forstandard(lzd%llr(ilr)%d%n1, lzd%llr(ilr)%d%n2, lzd%llr(ilr)%d%n3, &
        lzd%llr(ilr)%d%nfl1, lzd%llr(ilr)%d%nfu1, &
        lzd%llr(ilr)%d%nfl2, lzd%llr(ilr)%d%nfu2, lzd%llr(ilr)%d%nfl3, lzd%llr(ilr)%d%nfu3, &
        lzd%llr(ilr)%wfd%nseg_c, lzd%llr(ilr)%wfd%nvctr_c, lzd%llr(ilr)%wfd%keygloc, lzd%llr(ilr)%wfd%keyvloc, &
        lzd%llr(ilr)%wfd%nseg_f, lzd%llr(ilr)%wfd%nvctr_f, &
        lzd%llr(ilr)%wfd%keygloc(1,lzd%llr(ilr)%wfd%nseg_c+min(1,lzd%llr(ilr)%wfd%nseg_f)), &
        lzd%llr(ilr)%wfd%keyvloc(lzd%llr(ilr)%wfd%nseg_c+min(1,lzd%llr(ilr)%wfd%nseg_f)),  &
        scal, phiy_c, phiy_f, phid(isty_c), phid(isty_f))

   ! Compress the z wavefunction.
   istz_c = 1 + 2*(lzd%llr(ilr)%wfd%nvctr_c + 7*lzd%llr(ilr)%wfd%nvctr_f)
   istz_f = 1 + 3*lzd%llr(ilr)%wfd%nvctr_c + 2*7*lzd%llr(ilr)%wfd%nvctr_f
   call compress_forstandard(lzd%llr(ilr)%d%n1, lzd%llr(ilr)%d%n2, lzd%llr(ilr)%d%n3, &
        lzd%llr(ilr)%d%nfl1, lzd%llr(ilr)%d%nfu1, &
        lzd%llr(ilr)%d%nfl2, lzd%llr(ilr)%d%nfu2, lzd%llr(ilr)%d%nfl3, lzd%llr(ilr)%d%nfu3, &
        lzd%llr(ilr)%wfd%nseg_c, lzd%llr(ilr)%wfd%nvctr_c, lzd%llr(ilr)%wfd%keygloc, lzd%llr(ilr)%wfd%keyvloc, &
        lzd%llr(ilr)%wfd%nseg_f, lzd%llr(ilr)%wfd%nvctr_f, &
        lzd%llr(ilr)%wfd%keygloc(1,lzd%llr(ilr)%wfd%nseg_c+min(1,lzd%llr(ilr)%wfd%nseg_f)), &
        lzd%llr(ilr)%wfd%keyvloc(lzd%llr(ilr)%wfd%nseg_c+min(1,lzd%llr(ilr)%wfd%nseg_f)),  &
        scal, phiz_c, phiz_f, phid(istz_c), phid(istz_f))

   call deallocateWorkarrays()

contains

  subroutine allocateWorkarrays()
  
    ! THIS IS COPIED FROM allocate_work_arrays. Works only for free boundary.
    nf=(lzd%llr(ilr)%d%nfu1-lzd%llr(ilr)%d%nfl1+1)*(lzd%llr(ilr)%d%nfu2-lzd%llr(ilr)%d%nfl2+1)* &
       (lzd%llr(ilr)%d%nfu3-lzd%llr(ilr)%d%nfl3+1)

    ! Allocate work arrays
    allocate(w_c(0:lzd%llr(ilr)%d%n1,0:lzd%llr(ilr)%d%n2,0:lzd%llr(ilr)%d%n3+ndebug), stat=istat)
    call memocc(istat, w_c, 'w_c', subname)
    !!w_c=0.d0
    call to_zero((lzd%llr(ilr)%d%n1+1)*(lzd%llr(ilr)%d%n2+1)*(lzd%llr(ilr)%d%n3+1), w_c(0,0,0))

    allocate(w_f(7,lzd%llr(ilr)%d%nfl1:lzd%llr(ilr)%d%nfu1,lzd%llr(ilr)%d%nfl2:lzd%llr(ilr)%d%nfu2, &
                 lzd%llr(ilr)%d%nfl3:lzd%llr(ilr)%d%nfu3+ndebug), stat=istat)
    call memocc(istat, w_f, 'w_f', subname)
    !!w_f=0.d0
    call to_zero(7*nf, w_f(1,lzd%llr(ilr)%d%nfl1,lzd%llr(ilr)%d%nfl2,lzd%llr(ilr)%d%nfl3))

  
    allocate(w_f1(nf+ndebug), stat=istat)
    call memocc(istat, w_f1, 'w_f1', subname)
    !!w_f1=0.d0
    call to_zero(nf, w_f1(1))
    
    allocate(w_f2(nf+ndebug), stat=istat)
    call memocc(istat, w_f2, 'w_f2', subname)
    !!w_f2=0.d0
    call to_zero(nf, w_f2(1))

    allocate(w_f3(nf+ndebug), stat=istat)
    call memocc(istat, w_f3, 'w_f3', subname)
    !!w_f3=0.d0
    call to_zero(nf, w_f3(1))
  
  
    allocate(phix_f(7,lzd%llr(ilr)%d%nfl1:lzd%llr(ilr)%d%nfu1,lzd%llr(ilr)%d%nfl2:lzd%llr(ilr)%d%nfu2, &
                    lzd%llr(ilr)%d%nfl3:lzd%llr(ilr)%d%nfu3), stat=istat)
    call memocc(istat, phix_f, 'phix_f', subname)
    !!phix_f=0.d0
    call to_zero(7*nf, phix_f(1,lzd%llr(ilr)%d%nfl1,lzd%llr(ilr)%d%nfl2,lzd%llr(ilr)%d%nfl3))

    allocate(phix_c(0:lzd%llr(ilr)%d%n1,0:lzd%llr(ilr)%d%n2,0:lzd%llr(ilr)%d%n3), stat=istat)
    call memocc(istat, phix_c, 'phix_c', subname)
    !!phix_c=0.d0
    call to_zero((lzd%llr(ilr)%d%n1+1)*(lzd%llr(ilr)%d%n2+1)*(lzd%llr(ilr)%d%n3+1), phix_c(0,0,0))

    allocate(phiy_f(7,lzd%llr(ilr)%d%nfl1:lzd%llr(ilr)%d%nfu1,lzd%llr(ilr)%d%nfl2:lzd%llr(ilr)%d%nfu2, &
                    lzd%llr(ilr)%d%nfl3:lzd%llr(ilr)%d%nfu3), stat=istat)
    call memocc(istat, phiy_f, 'phiy_f', subname)
    !!phiy_f=0.d0
    call to_zero(7*nf, phiy_f(1,lzd%llr(ilr)%d%nfl1,lzd%llr(ilr)%d%nfl2,lzd%llr(ilr)%d%nfl3))

    allocate(phiy_c(0:lzd%llr(ilr)%d%n1,0:lzd%llr(ilr)%d%n2,0:lzd%llr(ilr)%d%n3), stat=istat)
    call memocc(istat, phiy_c, 'phiy_c', subname)
    !!phiy_c=0.d0
    call to_zero((lzd%llr(ilr)%d%n1+1)*(lzd%llr(ilr)%d%n2+1)*(lzd%llr(ilr)%d%n3+1), phiy_c(0,0,0))

    allocate(phiz_f(7,lzd%llr(ilr)%d%nfl1:lzd%llr(ilr)%d%nfu1,lzd%llr(ilr)%d%nfl2:lzd%llr(ilr)%d%nfu2, &
                    lzd%llr(ilr)%d%nfl3:lzd%llr(ilr)%d%nfu3), stat=istat)
    call memocc(istat, phiz_f, 'phiz_f', subname)
    !!phiz_f=0.d0
    call to_zero(7*nf, phiz_f(1,lzd%llr(ilr)%d%nfl1,lzd%llr(ilr)%d%nfl2,lzd%llr(ilr)%d%nfl3))

    allocate(phiz_c(0:lzd%llr(ilr)%d%n1,0:lzd%llr(ilr)%d%n2,0:lzd%llr(ilr)%d%n3), stat=istat)
    call memocc(istat, phiz_c, 'phiz_c', subname)
    !!phiz_c=0.d0
    call to_zero((lzd%llr(ilr)%d%n1+1)*(lzd%llr(ilr)%d%n2+1)*(lzd%llr(ilr)%d%n3+1), phiz_c(0,0,0))
  
  end subroutine allocateWorkarrays


  subroutine deallocateWorkarrays

    iall=-product(shape(w_c))*kind(w_c)
    deallocate(w_c, stat=istat)
    call memocc(istat, iall, 'w_c', subname)

    iall=-product(shape(w_f))*kind(w_f)
    deallocate(w_f, stat=istat)
    call memocc(istat, iall, 'w_f', subname)

    iall=-product(shape(w_f1))*kind(w_f1)
    deallocate(w_f1, stat=istat)
    call memocc(istat, iall, 'w_f1', subname)

    iall=-product(shape(w_f2))*kind(w_f2)
    deallocate(w_f2, stat=istat)
    call memocc(istat, iall, 'w_f2', subname)

    iall=-product(shape(w_f3))*kind(w_f3)
    deallocate(w_f3, stat=istat)
    call memocc(istat, iall, 'w_f3', subname)

    iall=-product(shape(phix_f))*kind(phix_f)
    deallocate(phix_f, stat=istat)
    call memocc(istat, iall, 'phix_f', subname)

    iall=-product(shape(phix_c))*kind(phix_c)
    deallocate(phix_c, stat=istat)
    call memocc(istat, iall, 'phix_c', subname)

    iall=-product(shape(phiy_f))*kind(phiy_f)
    deallocate(phiy_f, stat=istat)
    call memocc(istat, iall, 'phiy_f', subname)

    iall=-product(shape(phiy_c))*kind(phiy_c)
    deallocate(phiy_c, stat=istat)
    call memocc(istat, iall, 'phiy_c', subname)

    iall=-product(shape(phiz_f))*kind(phiz_f)
    deallocate(phiz_f, stat=istat)
    call memocc(istat, iall, 'phiz_f', subname)

    iall=-product(shape(phiz_c))*kind(phiz_c)
    deallocate(phiz_c, stat=istat)
    call memocc(istat, iall, 'phiz_c', subname)

  end subroutine deallocateWorkarrays

end subroutine get_one_derivative_supportfunction

!Experimenting the calculation of outward flux to determine the energy error associated with locrad
 subroutine correction_locrad(iproc, nproc, tmb, orbs, coeff)
 use module_base
 use module_types
 use module_interfaces
 implicit none

 integer, intent(in) :: iproc, nproc
 type(DFT_wavefunction), intent(in) :: tmb
 type(orbitals_data), intent(in) :: orbs
 real(kind=8), dimension(tmb%orbs%norb,orbs%norb), intent(in) :: coeff
 !Local variables
 integer :: ist, iorb, iiorb, ilr, ndim, ndimr, istat, iall
 integer :: i1 ,i2, i3, ipt, jjorb, kkorb, ierr
 real(kind=8) :: x, y, z, factor
 real(kind=8), allocatable, dimension(:) :: phider, psirX, psirY, psirZ, dphi, phidr, dE 
 real(kind=8), allocatable, dimension(:) :: psit_c, psit_f, hpsit_c, hpsit_f, phidr_c, phidr_f
 real(kind=8), allocatable, dimension(:) :: matrix_compr, overlap_compr
 real(kind=8), allocatable, dimension(:,:) :: matrix, overlap 
 type(workarr_sumrho) :: w
 character(len=*),parameter:: subname='correction_locrad'


  allocate(phidr(max(tmb%orbs%npsidim_orbs, tmb%orbs%npsidim_comp)), stat=istat)
  call memocc(istat,phidr,'phidr',subname)
  call to_zero(max(tmb%orbs%npsidim_orbs, tmb%orbs%npsidim_comp),phidr(1))
  
  ! First construct the radial derivatives
  ist = 1
  do iorb = 1, tmb%orbs%norbp
     iiorb = iorb + tmb%orbs%isorb
     ilr = tmb%orbs%inwhichlocreg(iiorb)
     ndim = tmb%lzd%llr(ilr)%wfd%nvctr_c + 7*tmb%lzd%llr(ilr)%wfd%nvctr_f
     ndimr= tmb%lzd%llr(ilr)%d%n1i*tmb%lzd%llr(ilr)%d%n2i*tmb%lzd%llr(ilr)%d%n3i

     allocate(phider(3*ndim),stat=istat)
     call memocc(istat,phider,'phider',subname)
     call to_zero(3*ndim, phider(1))

     !call get_derivative_supportfunctions(ndim, tmb%lzd%hgrids(1), tmb%lzd, tmb%orbs, tmb%psi(ist), phider)
     call get_one_derivative_supportfunction(ilr,tmb%lzd%hgrids(1),tmb%lzd,tmb%psi(ist),phider)    

     ! transform the derivatives to real space 
     call initialize_work_arrays_sumrho(tmb%lzd%llr(ilr),w)
     allocate(psirX(ndimr),stat=istat)
     call memocc(istat,psirX,'psirX',subname)
     allocate(psirY(ndimr),stat=istat)
     call memocc(istat,psirY,'psirY',subname)
     allocate(psirZ(ndimr),stat=istat)
     call memocc(istat,psirZ,'psirZ',subname)

     call daub_to_isf(tmb%lzd%llr(ilr),w,phider(1),psirX)
     call daub_to_isf(tmb%lzd%llr(ilr),w,phider(1+ndim),psirY)
     call daub_to_isf(tmb%lzd%llr(ilr),w,phider(1+2*ndim),psirZ)
     
     !Construct radial derivative
     ! Note the addition of 1.d-4 is to avoid division by zero (introduces a 1.d-8 error in the DE/dL)
     ! -16 because of buffers (redo with correct function for periodic)
     allocate(dphi(ndimr),stat=istat)
     call memocc(istat,dphi,'dphi',subname)
     call to_zero(ndimr, dphi(1))
     do i3= 1, tmb%lzd%llr(ilr)%d%n3i
        z = (tmb%lzd%llr(ilr)%nsi3 + i3-15)*0.5*tmb%lzd%hgrids(3) - tmb%lzd%llr(ilr)%locregCenter(3) 
        do i2= 1, tmb%lzd%llr(ilr)%d%n2i
           y= (tmb%lzd%llr(ilr)%nsi2 + i2-15)*0.5*tmb%lzd%hgrids(2) - tmb%lzd%llr(ilr)%locregCenter(2) 
           do i1= 1, tmb%lzd%llr(ilr)%d%n1i
              x = (tmb%lzd%llr(ilr)%nsi1 + i1-15)*0.5*tmb%lzd%hgrids(1) - tmb%lzd%llr(ilr)%locregCenter(1) 
              ipt = (i3-1)*tmb%lzd%llr(ilr)%d%n2i*tmb%lzd%llr(ilr)%d%n1i + (i2-1)*tmb%lzd%llr(ilr)%d%n1i + i1
              factor = sqrt(x**2+y**2+z**2)
              if (factor /=0.0_wp) &
               dphi(ipt) = dphi(ipt) + psirX(ipt)*x/factor + psirY(ipt)*y/factor + psirZ(ipt)*z/factor
           end do
        end do
     end do

     ! transform back to daub
     call isf_to_daub(tmb%lzd%llr(ilr),w,dphi,phidr(ist))

     ist = ist + ndim
     call deallocate_work_arrays_sumrho(w)        
     iall = -product(shape(psirX))*kind(psirX)
     deallocate(psirX,stat=istat)
     call memocc(istat,iall,'psirX',subname)
     iall = -product(shape(psirY))*kind(psirY)
     deallocate(psirY,stat=istat)
     call memocc(istat,iall,'psirY',subname)
     iall = -product(shape(psirZ))*kind(psirZ)
     deallocate(psirZ,stat=istat)
     call memocc(istat,iall,'psirZ',subname)
     iall = -product(shape(phider))*kind(phider)
     deallocate(phider,stat=istat)
     call memocc(istat,iall,'phider',subname)
     iall = -product(shape(dphi))*kind(dphi)
     deallocate(dphi,stat=istat)
     call memocc(istat,iall,'dphi',subname)
    
  end do

  ! Suppose we have the h|psi> in tmb%hpsi
  allocate(phidr_c(tmb%collcom%ndimind_c),stat=istat)
  call memocc(istat,phidr_c,'phidr_c',subname)
  allocate(phidr_f(7*tmb%collcom%ndimind_f),stat=istat)
  call memocc(istat,phidr_f,'phidr_f',subname)
  call transpose_localized(iproc,nproc,tmb%orbs,tmb%collcom, &
       phidr, phidr_c, phidr_f, tmb%lzd)
  iall = -product(shape(phidr))*kind(phidr)
  deallocate(phidr,stat=istat)
  call memocc(istat,iall,'phidr',subname)
  
  
  allocate(psit_c(tmb%collcom%ndimind_c),stat=istat)
  call memocc(istat,psit_c,'psit_c',subname)
  allocate(psit_f(7*tmb%collcom%ndimind_f),stat=istat)
  call memocc(istat,psit_f,'psit_f',subname)
  allocate(hpsit_c(tmb%collcom%ndimind_c),stat=istat)
  call memocc(istat,hpsit_c,'hpsit_c',subname)
  allocate(hpsit_f(7*tmb%collcom%ndimind_f),stat=istat)
  call memocc(istat,hpsit_f,'hpsit_f',subname)

  call transpose_localized(iproc, nproc, tmb%orbs,  tmb%collcom, &
       tmb%psi, psit_c, psit_f, tmb%lzd)
  call transpose_localized(iproc, nproc, tmb%orbs,  tmb%collcom, &
       tmb%hpsi, hpsit_c, hpsit_f, tmb%lzd)


  allocate(matrix(tmb%orbs%norb,tmb%orbs%norb),stat=istat)
  call memocc(istat,matrix,'matrix',subname)  
  call to_zero(tmb%orbs%norb**2,matrix(1,1))
  allocate(overlap(tmb%orbs%norb,tmb%orbs%norb),stat=istat)
  call memocc(istat,overlap,'overlap',subname)  
  call to_zero(tmb%orbs%norb**2,overlap(1,1))

  allocate(matrix_compr(tmb%mad%nvctr),stat=istat)
  call memocc(istat,matrix_compr,'matrix_compr',subname)
  allocate(overlap_compr(tmb%mad%nvctr),stat=istat)
  call memocc(istat,overlap_compr,'overlap_compr',subname)

  call calculate_overlap_transposed(iproc, nproc, tmb%orbs, tmb%mad, &
       tmb%collcom, hpsit_c, phidr_c, hpsit_f, phidr_f, matrix_compr)
  call calculate_overlap_transposed(iproc, nproc, tmb%orbs, tmb%mad, &
       tmb%collcom, psit_c, phidr_c, psit_f, phidr_f, overlap_compr)

  call uncompressMatrix(tmb%orbs%norb, tmb%mad, matrix_compr, matrix)
  call uncompressMatrix(tmb%orbs%norb, tmb%mad, overlap_compr, overlap)

  iall = -product(shape(matrix_compr))*kind(matrix_compr)
  deallocate(matrix_compr,stat=istat)
  call memocc(istat,iall,'matrix_compr',subname)

  iall = -product(shape(overlap_compr))*kind(overlap_compr)
  deallocate(overlap_compr,stat=istat)
  call memocc(istat,iall,'overlap_compr',subname)

  allocate(dE(tmb%orbs%norb),stat=istat)
  call memocc(istat,dE,'dE',subname)  
  call to_zero(tmb%orbs%norb, dE(1))
  do iorb=1,orbs%norbp
      iiorb=orbs%isorb+iorb
      do jjorb=1,tmb%orbs%norb
          do kkorb=1,tmb%orbs%norb
             dE(jjorb) = dE(jjorb) - &
              2*orbs%occup(iiorb)*coeff(jjorb,iiorb)*coeff(kkorb,iiorb)* &
              (matrix(jjorb,kkorb) - orbs%eval(iiorb)*overlap(jjorb,kkorb))
          end do
      end do
  end do
  call mpiallred(dE(1), tmb%orbs%norb, mpi_sum, bigdft_mpi%mpi_comm, ierr)
  !!if(iproc==0) then
  !!     do iorb=1,tmb%orbs%norb
  !!         print *,'Basis function ',iorb,'on atom',tmb%orbs%onwhichatom(iorb)
  !!         print *,' has an energy error of ', dE(iorb)
  !!     end do
  !!     print *,'Total variation is of :',sum(dE)
  !!end if

  iall = -product(shape(dE))*kind(dE)
  deallocate(dE,stat=istat)
  call memocc(istat,iall,'dE',subname)
  iall = -product(shape(overlap))*kind(overlap)
  deallocate(overlap,stat=istat)
  call memocc(istat,iall,'overlap',subname)
  iall = -product(shape(matrix))*kind(matrix)
  deallocate(matrix,stat=istat)
  call memocc(istat,iall,'matrix',subname)
  iall = -product(shape(phidr_c))*kind(phidr_c)
  deallocate(phidr_c,stat=istat)
  call memocc(istat,iall,'phidr_c',subname)
  iall = -product(shape(phidr_f))*kind(phidr_f)
  deallocate(phidr_f,stat=istat)
  call memocc(istat,iall,'phidr_f',subname)
  iall = -product(shape(psit_c))*kind(psit_c)
  deallocate(psit_c,stat=istat)
  call memocc(istat,iall,'psit_c',subname)
  iall = -product(shape(psit_f))*kind(psit_f)
  deallocate(psit_f,stat=istat)
  call memocc(istat,iall,'psit_f',subname)
  iall = -product(shape(hpsit_c))*kind(hpsit_c)
  deallocate(hpsit_c,stat=istat)
  call memocc(istat,iall,'hpsit_c',subname)
  iall = -product(shape(hpsit_f))*kind(hpsit_f)
  deallocate(hpsit_f,stat=istat)
  call memocc(istat,iall,'hpsit_f',subname)

  !First going to need the divergence of the derivative basis, to do this we need to derive it again.
  !allocate(psidiv(tmb%orbs_shamop%npsidim_orbs), stat=istat)
  !call memocc(istat, psidiv, 'psidiv', subname)
  !allocate(outflux(tmb%orbs_shamop%norb),stat=istat)
  !call memocc(istat, outflux, 'outflux', subname)

  !call get_divergence(tmb%orbs_shamop%npsidim_orbs,tmb%lzd_shamop%hgrids(1), tmb%lzd_shamop, tmb%orbs_shamop, tmbder%psi, psidiv)

  ! Now integrate the divergence only in the outer region (corresponds to a shell of 32 isf points because of large region).
  !!ldir = 1
  !!outflux = 0.0_dp
  !!do iorb = 1, tmb%orbs_shamop%norbp
  !!   iiorb = iorb + tmb%orbs_shamop%isorb
  !!   ilr = tmb%orbs_shamop%inwhichlocreg(iiorb)
  !!   call initialize_work_arrays_sumrho(tmb%lzd_shamop%llr(ilr),w)
  !!   allocate(psir(tmb%lzd_shamop%llr(ilr)%d%n1i*tmb%lzd_shamop%llr(ilr)%d%n2i*tmb%lzd_shamop%llr(ilr)%d%n3i),stat=istat)
  !!   call memocc(istat,psir,'psir',subname)
  !!   call daub_to_isf(tmb%lzd_shamop%llr(ilr),w,psidiv(ldir),psir)
  !!   !call daub_to_isf(tmb%lzd_shamop%llr(ilr),w,tmb%psi_shamop(ldir),psir)
  !!   do i1 = 1, tmb%lzd_shamop%llr(ilr)%d%n1i
  !!      if(i1 > 32 .and. tmb%lzd_shamop%llr(ilr)%d%n1i-i1 > 32) cycle
  !!      do i2 = 1 , tmb%lzd_shamop%llr(ilr)%d%n2i
  !!         if(i2 > 32 .and. tmb%lzd_shamop%llr(ilr)%d%n2i-i2 > 32) cycle
  !!         do i3 = 1, tmb%lzd_shamop%llr(ilr)%d%n3i 
  !!            if(i3 > 32 .and. tmb%lzd_shamop%llr(ilr)%d%n3i-i3 > 32) cycle
  !!            ipt = (i3-1)*tmb%lzd_shamop%llr(ilr)%d%n2i*tmb%lzd_shamop%llr(ilr)%d%n1i + (i2-1)*tmb%lzd_shamop%llr(ilr)%d%n1i + i1
  !!            outflux(iiorb) = outflux(iiorb) + psir(ipt)!*psir(ipt)
  !!         end do
  !!      end do 
  !!   end do
  !!   call deallocate_work_arrays_sumrho(w)
  !!   iall = -product(shape(psir))*kind(psir)
  !!   deallocate(psir,stat=istat)
  !!   call memocc(istat,iall,'psir',subname)
  !!   ldir = ldir + tmb%lzd_shamop%llr(ilr)%wfd%nvctr_c + 7*tmb%lzd_shamop%llr(ilr)%wfd%nvctr_f
  !!end do

  !!call mpiallred(outflux(1),tmb%orbs_shamop%norb,MPI_SUM,bigdft_mpi%mpi_comm,ierr)

  !!if(iproc == 0) then
  !!   factor = 0.5*tmb%lzd_shamop%hgrids(1)*0.5*tmb%lzd_shamop%hgrids(2)*0.5*tmb%lzd_shamop%hgrids(3)
  !!   do iiorb = 1, tmb%orbs_shamop%norb
  !!      print *,'Basis function ',iiorb,'on atom',tmb%orbs_shamop%onwhichatom(iiorb)
  !!      print *,' has an outward flux of ', outflux(iiorb)!*factor
  !!   end do
  !!end if

  !!iall = -product(shape(psidiv))*kind(psidiv)
  !!deallocate(psidiv, stat=istat)
  !!call memocc(istat, iall, 'psidiv', subname)

end subroutine correction_locrad


subroutine get_derivative(idir, ndim, hgrid, orbs, lzd, phi, phider)
  use module_base
  use module_types
  use module_interfaces!, except_this_one => get_divergence
  implicit none
  
  ! Calling arguments
  integer,intent(in):: ndim, idir
  real(kind=8),intent(in) :: hgrid
  type(orbitals_data),intent(in) :: orbs
  type(local_zone_descriptors),intent(in) :: lzd
  real(kind=8),dimension(ndim),intent(in) :: phi !< Basis functions
  real(kind=8),dimension(ndim),intent(inout) :: phider  !< Derivative of basis functions
  
  ! Local variables
  integer :: nf, istat, iall, iorb, iiorb, ilr, istrt
  real(kind=8),dimension(0:3),parameter :: scal=1.d0
  real(kind=8),dimension(:),allocatable :: w_f1, w_f2, w_f3
  real(kind=8),dimension(:,:,:),allocatable :: w_c, phider_c
  real(kind=8),dimension(:,:,:,:),allocatable :: w_f, phider_f
  character(len=*),parameter :: subname='get_derivative'

   call to_zero(ndim,phider(1))

   istrt = 1
   do iorb=1, orbs%norbp
     iiorb = iorb+orbs%isorb 
     ilr = orbs%inwhichlocreg(iiorb)

     call allocateWorkarrays()
 
     ! Uncompress the wavefunction.
     call uncompress_forstandard(lzd%llr(ilr)%d%n1, lzd%llr(ilr)%d%n2, lzd%llr(ilr)%d%n3, &
          lzd%llr(ilr)%d%nfl1, lzd%llr(ilr)%d%nfu1, & 
          lzd%llr(ilr)%d%nfl2, lzd%llr(ilr)%d%nfu2, lzd%llr(ilr)%d%nfl3, lzd%llr(ilr)%d%nfu3,  &
          lzd%llr(ilr)%wfd%nseg_c, lzd%llr(ilr)%wfd%nvctr_c, lzd%llr(ilr)%wfd%keygloc, lzd%llr(ilr)%wfd%keyvloc,  &
          lzd%llr(ilr)%wfd%nseg_f, lzd%llr(ilr)%wfd%nvctr_f, &
          lzd%llr(ilr)%wfd%keygloc(1,lzd%llr(ilr)%wfd%nseg_c+min(1,lzd%llr(ilr)%wfd%nseg_f)), &
          lzd%llr(ilr)%wfd%keyvloc(lzd%llr(ilr)%wfd%nseg_c+min(1,lzd%llr(ilr)%wfd%nseg_f)),  &
          scal, phi(istrt), phi(istrt+lzd%llr(ilr)%wfd%nvctr_c), w_c, w_f, w_f1, w_f2, w_f3)
  
     if(idir==1) then
         call createDerivativeX(lzd%llr(ilr)%d%n1, lzd%llr(ilr)%d%n2, lzd%llr(ilr)%d%n3, &
              lzd%llr(ilr)%d%nfl1, lzd%llr(ilr)%d%nfu1, lzd%llr(ilr)%d%nfl2, lzd%llr(ilr)%d%nfu2, &
              lzd%llr(ilr)%d%nfl3, lzd%llr(ilr)%d%nfu3,  &
              hgrid, lzd%llr(ilr)%bounds%kb%ibyz_c, lzd%llr(ilr)%bounds%kb%ibyz_f, &
              w_c, w_f, w_f1, phider_c, phider_f)
         ! Compress the x wavefunction.
         call compress_forstandard(lzd%llr(ilr)%d%n1, lzd%llr(ilr)%d%n2, lzd%llr(ilr)%d%n3, &
              lzd%llr(ilr)%d%nfl1, lzd%llr(ilr)%d%nfu1, &
              lzd%llr(ilr)%d%nfl2, lzd%llr(ilr)%d%nfu2, lzd%llr(ilr)%d%nfl3, lzd%llr(ilr)%d%nfu3, &
              lzd%llr(ilr)%wfd%nseg_c, lzd%llr(ilr)%wfd%nvctr_c, lzd%llr(ilr)%wfd%keygloc, lzd%llr(ilr)%wfd%keyvloc, &
              lzd%llr(ilr)%wfd%nseg_f, lzd%llr(ilr)%wfd%nvctr_f, &
              lzd%llr(ilr)%wfd%keygloc(1,lzd%llr(ilr)%wfd%nseg_c+min(1,lzd%llr(ilr)%wfd%nseg_f)), &
              lzd%llr(ilr)%wfd%keyvloc(lzd%llr(ilr)%wfd%nseg_c+min(1,lzd%llr(ilr)%wfd%nseg_f)),  &
              scal, phider_c, phider_f, phider(istrt), phider(istrt+lzd%llr(ilr)%wfd%nvctr_c))
     else if (idir==2) then
        call createDerivativeY(lzd%llr(ilr)%d%n1, lzd%llr(ilr)%d%n2, lzd%llr(ilr)%d%n3, &
              lzd%llr(ilr)%d%nfl1, lzd%llr(ilr)%d%nfu1, lzd%llr(ilr)%d%nfl2, lzd%llr(ilr)%d%nfu2, &
              lzd%llr(ilr)%d%nfl3, lzd%llr(ilr)%d%nfu3,  &
              hgrid, lzd%llr(ilr)%bounds%kb%ibxz_c, lzd%llr(ilr)%bounds%kb%ibxz_f, &
              w_c, w_f, w_f2, phider_c, phider_f)
        ! Compress the y wavefunction.
        call compress_forstandard(lzd%llr(ilr)%d%n1, lzd%llr(ilr)%d%n2, lzd%llr(ilr)%d%n3, &
             lzd%llr(ilr)%d%nfl1, lzd%llr(ilr)%d%nfu1, &
             lzd%llr(ilr)%d%nfl2, lzd%llr(ilr)%d%nfu2, lzd%llr(ilr)%d%nfl3, lzd%llr(ilr)%d%nfu3, &
             lzd%llr(ilr)%wfd%nseg_c, lzd%llr(ilr)%wfd%nvctr_c, lzd%llr(ilr)%wfd%keygloc, lzd%llr(ilr)%wfd%keyvloc, &
             lzd%llr(ilr)%wfd%nseg_f, lzd%llr(ilr)%wfd%nvctr_f, &
             lzd%llr(ilr)%wfd%keygloc(1,lzd%llr(ilr)%wfd%nseg_c+min(1,lzd%llr(ilr)%wfd%nseg_f)), &
             lzd%llr(ilr)%wfd%keyvloc(lzd%llr(ilr)%wfd%nseg_c+min(1,lzd%llr(ilr)%wfd%nseg_f)),  &
             scal, phider_c, phider_f, phider(istrt), phider(istrt+lzd%llr(ilr)%wfd%nvctr_c))
     else
        call createDerivativeZ(lzd%llr(ilr)%d%n1, lzd%llr(ilr)%d%n2, lzd%llr(ilr)%d%n3, &
              lzd%llr(ilr)%d%nfl1, lzd%llr(ilr)%d%nfu1, lzd%llr(ilr)%d%nfl2, lzd%llr(ilr)%d%nfu2, &
              lzd%llr(ilr)%d%nfl3, lzd%llr(ilr)%d%nfu3,  &
              hgrid, lzd%llr(ilr)%bounds%kb%ibxy_c, lzd%llr(ilr)%bounds%kb%ibxy_f, &
              w_c, w_f, w_f3, phider_c, phider_f)
        ! Compress the z wavefunction.
        call compress_forstandard(lzd%llr(ilr)%d%n1, lzd%llr(ilr)%d%n2, lzd%llr(ilr)%d%n3, &
             lzd%llr(ilr)%d%nfl1, lzd%llr(ilr)%d%nfu1, &
             lzd%llr(ilr)%d%nfl2, lzd%llr(ilr)%d%nfu2, lzd%llr(ilr)%d%nfl3, lzd%llr(ilr)%d%nfu3, &
             lzd%llr(ilr)%wfd%nseg_c, lzd%llr(ilr)%wfd%nvctr_c, lzd%llr(ilr)%wfd%keygloc, lzd%llr(ilr)%wfd%keyvloc, &
             lzd%llr(ilr)%wfd%nseg_f, lzd%llr(ilr)%wfd%nvctr_f, &
             lzd%llr(ilr)%wfd%keygloc(1,lzd%llr(ilr)%wfd%nseg_c+min(1,lzd%llr(ilr)%wfd%nseg_f)), &
             lzd%llr(ilr)%wfd%keyvloc(lzd%llr(ilr)%wfd%nseg_c+min(1,lzd%llr(ilr)%wfd%nseg_f)),  &
             scal, phider_c, phider_f, phider(istrt),phider(istrt+lzd%llr(ilr)%wfd%nvctr_c ))
     end if

     istrt = istrt + lzd%llr(ilr)%wfd%nvctr_c + 7*lzd%llr(ilr)%wfd%nvctr_f
  
     call deallocateWorkarrays()                                

   end do

contains
  subroutine allocateWorkarrays()

    ! THIS IS COPIED FROM allocate_work_arrays. Works only for free boundary.
    nf=(lzd%llr(ilr)%d%nfu1-lzd%llr(ilr)%d%nfl1+1)*(lzd%llr(ilr)%d%nfu2-lzd%llr(ilr)%d%nfl2+1)* &
       (lzd%llr(ilr)%d%nfu3-lzd%llr(ilr)%d%nfl3+1)

    ! Allocate work arrays
    allocate(w_c(0:lzd%llr(ilr)%d%n1,0:lzd%llr(ilr)%d%n2,0:lzd%llr(ilr)%d%n3+ndebug), stat=istat)
    call memocc(istat, w_c, 'w_c', subname)
    !!w_c=0.d0
    call to_zero((lzd%llr(ilr)%d%n1+1)*(lzd%llr(ilr)%d%n2+1)*(lzd%llr(ilr)%d%n3+1), w_c(0,0,0))

    allocate(w_f(7,lzd%llr(ilr)%d%nfl1:lzd%llr(ilr)%d%nfu1,lzd%llr(ilr)%d%nfl2:lzd%llr(ilr)%d%nfu2, &
                 lzd%llr(ilr)%d%nfl3:lzd%llr(ilr)%d%nfu3+ndebug), stat=istat)
    call memocc(istat, w_f, 'w_f', subname)
    !!w_f=0.d0
    call to_zero(7*nf, w_f(1,lzd%llr(ilr)%d%nfl1,lzd%llr(ilr)%d%nfl2,lzd%llr(ilr)%d%nfl3))


    allocate(w_f1(nf+ndebug), stat=istat)
    call memocc(istat, w_f1, 'w_f1', subname)
    !!w_f1=0.d0
    call to_zero(nf, w_f1(1))

    allocate(w_f2(nf+ndebug), stat=istat)
    call memocc(istat, w_f2, 'w_f2', subname)
    !!w_f2=0.d0
    call to_zero(nf, w_f2(1))

    allocate(w_f3(nf+ndebug), stat=istat)
    call memocc(istat, w_f3, 'w_f3', subname)
    !!w_f3=0.d0
    call to_zero(nf, w_f3(1))


    allocate(phider_f(7,lzd%llr(ilr)%d%nfl1:lzd%llr(ilr)%d%nfu1,lzd%llr(ilr)%d%nfl2:lzd%llr(ilr)%d%nfu2, &
                    lzd%llr(ilr)%d%nfl3:lzd%llr(ilr)%d%nfu3), stat=istat)
    call memocc(istat, phider_f, 'phider_f', subname)
    !!phix_f=0.d0
    call to_zero(7*nf, phider_f(1,lzd%llr(ilr)%d%nfl1,lzd%llr(ilr)%d%nfl2,lzd%llr(ilr)%d%nfl3))

    allocate(phider_c(0:lzd%llr(ilr)%d%n1,0:lzd%llr(ilr)%d%n2,0:lzd%llr(ilr)%d%n3), stat=istat)
    call memocc(istat, phider_c, 'phider_c', subname)
    !!phix_c=0.d0
    call to_zero((lzd%llr(ilr)%d%n1+1)*(lzd%llr(ilr)%d%n2+1)*(lzd%llr(ilr)%d%n3+1), phider_c(0,0,0))

  end subroutine allocateWorkarrays


  subroutine deallocateWorkarrays

    iall=-product(shape(w_c))*kind(w_c)
    deallocate(w_c, stat=istat)
    call memocc(istat, iall, 'w_c', subname)

    iall=-product(shape(w_f))*kind(w_f)
    deallocate(w_f, stat=istat)
    call memocc(istat, iall, 'w_f', subname)

    iall=-product(shape(w_f1))*kind(w_f1)
    deallocate(w_f1, stat=istat)
    call memocc(istat, iall, 'w_f1', subname)

    iall=-product(shape(w_f2))*kind(w_f2)
    deallocate(w_f2, stat=istat)
    call memocc(istat, iall, 'w_f2', subname)

    iall=-product(shape(w_f3))*kind(w_f3)
    deallocate(w_f3, stat=istat)
    call memocc(istat, iall, 'w_f3', subname)

    iall=-product(shape(phider_f))*kind(phider_f)
    deallocate(phider_f, stat=istat)
    call memocc(istat, iall, 'phider_f', subname)

    iall=-product(shape(phider_c))*kind(phider_c)
    deallocate(phider_c, stat=istat)
    call memocc(istat, iall, 'phider_c', subname)

  end subroutine deallocateWorkarrays
end subroutine get_derivative

