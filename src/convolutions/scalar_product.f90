subroutine NL_HGH_application(npspcode,psppar,ncplx_p,n_p,wfd_p,proj,&
     ncplx_w,n_w,wfd_w,psi,hpsi,eproj)
  use module_base
  use module_types, only: wavefunctions_descriptors
  implicit none
  integer, intent(in) :: npspcode
  integer, intent(in) :: ncplx_p,n_p,ncplx_w,n_w
  type(wavefunctions_descriptors), intent(in) :: wfd_p !< descriptors of projectors
  type(wavefunctions_descriptors), intent(in) :: wfd_w !< descriptors of wavefunction
  real(gp), dimension(0:4,0:6), intent(in) :: psppar !< parameters of HGH PSP
  real(wp), dimension(wfd_p%nvctr_c+7*wfd_p%nvctr_f,ncplx_p,n_p), intent(in) :: proj !< components of the projectors, real and imaginary parts
  real(wp), dimension(wfd_w%nvctr_c+7*wfd_w%nvctr_f,ncplx_w,n_w), intent(in) :: psi !< components of wavefunctions, real and imaginary parts

  real(wp), dimension(wfd_w%nvctr_c+7*wfd_w%nvctr_f,ncplx_w,n_w), intent(inout) :: hpsi !< components of wavefunctions, real and imaginary parts
  real(gp), intent(out) :: eproj

  !local variables
  character(len=*), parameter :: subname='NL_HGH_application'
  integer :: nmseg_c,nmseg_f
  integer, dimension(:), allocatable :: nbsegs_cf,keyag_lin_cf
  integer, dimension(:,:), allocatable :: psi_mask
  real(wp), dimension(:,:), allocatable :: psi_pack
  real(wp), dimension(:,:,:), allocatable :: pdpsi,hpdpsi
  real(wp), dimension(:,:,:,:), allocatable :: scpr

  call f_routine(id=subname)

  eproj=0.0_gp

  !allocate temporary workspace
  keyag_lin_cf=f_malloc(wfd_w%nseg_c+wfd_w%nseg_f,id='keyag_lin_cf')


  nbsegs_cf=f_malloc0(wfd_p%nseg_c+wfd_p%nseg_f,id='nbsegs_cf')  
  !find the dimension for masking array, coarse case
  call vcopy(wfd_w%nseg_c+wfd_w%nseg_f,wfd_w%keyglob(1,1),2,keyag_lin_cf(1),1)

  call count_wblas_segs(wfd_w%nseg_c,wfd_p%nseg_c,keyag_lin_cf(1),&
       wfd_w%keyglob(1,1),wfd_p%keyglob(1,1),nbsegs_cf(1))
!  print *,'no of points',sum(nbsegs_cf),wfd_w%nseg_c,wfd_p%nseg_c
  call integrate_nseg(wfd_p%nseg_c,nbsegs_cf(1),nmseg_c)
!  print *,'no of points',nmseg_c

  if (wfd_w%nseg_f >0 .and. wfd_p%nseg_f > 0 ) then
     call count_wblas_segs(wfd_w%nseg_f,wfd_p%nseg_f,keyag_lin_cf(wfd_w%nseg_c+1),&
          wfd_w%keyglob(1,wfd_w%nseg_c+1),wfd_p%keyglob(1,wfd_p%nseg_c+1),&
          nbsegs_cf(wfd_p%nseg_c+1))
     call integrate_nseg(wfd_p%nseg_f,nbsegs_cf(wfd_p%nseg_c+1),nmseg_f)
  else
     nmseg_f=0
  end if
  !the masking array can be allocated

  !print *,'no of points',nmseg_c,nmseg_f
  psi_mask=f_malloc0((/3,nmseg_c+nmseg_f/),id='psi_mask')
  !and filled
  call fill_wblas_segs(wfd_w%nseg_c,wfd_p%nseg_c,nmseg_c,&
       nbsegs_cf(1),keyag_lin_cf(1),wfd_w%keyglob(1,1),wfd_p%keyglob(1,1),&
       wfd_w%keyvglob(1),wfd_p%keyvglob(1),psi_mask(1,1))
  if (nmseg_f > 0) then
     call fill_wblas_segs(wfd_w%nseg_f,wfd_p%nseg_f,nmseg_f,&
          nbsegs_cf(wfd_p%nseg_c+1),keyag_lin_cf(wfd_w%nseg_c+1),&
          wfd_w%keyglob(1,wfd_w%nseg_c+1),wfd_p%keyglob(1,wfd_p%nseg_c+1),&
          wfd_w%keyvglob(wfd_w%nseg_c+1),wfd_p%keyvglob(wfd_p%nseg_c+1),&
          psi_mask(1,nmseg_c+1))
  end if
  call f_free(keyag_lin_cf,nbsegs_cf)

  !create other workspaces  
  psi_pack=f_malloc0((/wfd_p%nvctr_c+7*wfd_p%nvctr_f,n_w*ncplx_w/),id='psi_pack')
  scpr=f_malloc((/ncplx_w,n_w,ncplx_p,n_p/),id='scpr')
  pdpsi=f_malloc((/max(ncplx_w,ncplx_p),n_w,n_p/),id='pdpsi')
  hpdpsi=f_malloc((/max(ncplx_w,ncplx_p),n_w,n_p/),id='hpdpsi')

  call proj_dot_psi(n_p*ncplx_p,wfd_p,proj,n_w*ncplx_w,wfd_w,psi,&
       nmseg_c,nmseg_f,psi_mask,psi_pack,scpr)

  !first create the coefficients for the application of the matrix
  !pdpsi = < p_i | psi >
  call full_coefficients('C',ncplx_p,n_p,'N',ncplx_w,n_w,scpr,'N',pdpsi)

  call apply_hij_coeff(npspcode,psppar,max(ncplx_w,ncplx_p)*n_w,n_p,pdpsi,hpdpsi)

  !then create the coefficients for the evaluation of the projector energy
  !pdpsi= < psi | p_i> = conj(< p_i | psi >)
  call full_coefficients('N',ncplx_p,n_p,'C',ncplx_w,n_w,scpr,'C',pdpsi)

  !the energy can be calculated here
  eproj=dot(max(ncplx_p,ncplx_w)*n_w*n_p,pdpsi(1,1,1),1,hpdpsi(1,1,1),1)

  !then the coefficients have to be transformed for the projectors
  call reverse_coefficients(ncplx_p,n_p,ncplx_w,n_w,hpdpsi,scpr)

  call scpr_proj_p_hpsi(n_p*ncplx_p,wfd_p,proj,n_w*ncplx_w,wfd_w,nmseg_c,nmseg_f,psi_mask,psi_pack,scpr,hpsi)

  !free workspaces
  call f_free(psi_mask)
  call f_free(psi_pack)
  call f_free(scpr)     
  call f_free(pdpsi)
  call f_free(hpdpsi)

  call f_release_routine()

contains

  !> routine for applying the coefficients needed HGH-type PSP to the scalar product
  !! among wavefunctions and projectors. The coefficients are real therefore 
  !! there is no need to separate scpr in its real and imaginary part before
  pure subroutine apply_hij_coeff(npspcode,psppar,n_w,n_p,scpr,hscpr)
    use module_base, only: gp
    implicit none
    integer, intent(in) :: n_p,n_w,npspcode
    real(gp), dimension(0:4,0:6), intent(in) :: psppar
    real(gp), dimension(n_w,n_p), intent(in) :: scpr
    real(gp), dimension(n_w,n_p), intent(out) :: hscpr
    !local variables
    integer :: i,j,l,m,iproj,iw
    real(gp), dimension(3,3,4) :: hij
    real(gp), dimension(7,3,4) :: cproj,dproj 

    !fill the hij matrix
    call hgh_hij_matrix(npspcode,psppar,hij)

    reversed_loop: do iw=1,n_w
       dproj=0.0_gp

       iproj=1
       !fill the complete coefficients
       do l=1,4 !diagonal in l
          do i=1,3
             if (psppar(l,i) /= 0.0_gp) then
                do m=1,2*l-1
                   cproj(m,i,l)=scpr(iw,iproj)
                   iproj=iproj+1
                end do
             else
                do m=1,2*l-1
                   cproj(m,i,l)=0.0_gp
                end do
             end if
          end do
       end do

       !applies the hij matrix
       do l=1,4 !diagonal in l
          do i=1,3
             do j=1,3
                do m=1,2*l-1 !diagonal in m
                   dproj(m,i,l)=dproj(m,i,l)+&
                        hij(i,j,l)*cproj(m,j,l)
                end do
             end do
          end do
       end do

       !copy back the coefficient
       iproj=1
       !fill the complete coefficients
       do l=1,4 !diagonal in l
          do i=1,3
             if (psppar(l,i) /= 0.0_gp) then
                do m=1,2*l-1
                   hscpr(iw,iproj)=dproj(m,i,l)
                   iproj=iproj+1
                end do
             end if
          end do
       end do
    end do reversed_loop

  end subroutine apply_hij_coeff

  pure subroutine hgh_hij_matrix(npspcode,psppar,hij)
    use module_base
    implicit none
    !Arguments
    integer, intent(in) :: npspcode
    real(gp), dimension(0:4,0:6), intent(in) :: psppar
    real(gp), dimension(3,3,4), intent(out) :: hij
    !Local variables
    integer :: l,i,j
    real(gp), dimension(2,2,3) :: offdiagarr

    !enter the coefficients for the off-diagonal terms (HGH case, npspcode=3)
    offdiagarr(1,1,1)=-0.5_gp*sqrt(3._gp/5._gp)
    offdiagarr(2,1,1)=-0.5_gp*sqrt(100._gp/63._gp)
    offdiagarr(1,2,1)=0.5_gp*sqrt(5._gp/21._gp)
    offdiagarr(2,2,1)=0.0_gp !never used
    offdiagarr(1,1,2)=-0.5_gp*sqrt(5._gp/7._gp)  
    offdiagarr(2,1,2)=-7._gp/3._gp*sqrt(1._gp/11._gp)
    offdiagarr(1,2,2)=1._gp/6._gp*sqrt(35._gp/11._gp)
    offdiagarr(2,2,2)=0.0_gp !never used
    offdiagarr(1,1,3)=-0.5_gp*sqrt(7._gp/9._gp)
    offdiagarr(2,1,3)=-9._gp*sqrt(1._gp/143._gp)
    offdiagarr(1,2,3)=0.5_gp*sqrt(63._gp/143._gp)
    offdiagarr(2,2,3)=0.0_gp !never used

    !  call to_zero(3*3*4,hij(1,1,1))
    hij=0.0_gp

    do l=1,4
       !term for all npspcodes
       loop_diag: do i=1,3
          hij(i,i,l)=psppar(l,i) !diagonal term
          if ((npspcode == 3 .and. l/=4 .and. i/=3) .or. &
               ((npspcode == 10 .or. npspcode == 12) .and. i/=3)) then !HGH(-K) case, offdiagonal terms
             loop_offdiag: do j=i+1,3
                if (psppar(l,j) == 0.0_gp) exit loop_offdiag
                !offdiagonal HGH term
                if (npspcode == 3) then !traditional HGH convention
                   hij(i,j,l)=offdiagarr(i,j-i,l)*psppar(l,j)
                else !HGH-K convention
                   hij(i,j,l)=psppar(l,i+j+1)
                end if
                hij(j,i,l)=hij(i,j,l) !symmetrization
             end do loop_offdiag
          end if
       end do loop_diag
    end do

  end subroutine hgh_hij_matrix

  pure subroutine reverse_coefficients(ncplx_p,n_p,ncplx_w,n_w,pdpsi,scpr)
    implicit none
    integer, intent(in) :: ncplx_p,ncplx_w,n_p,n_w
    real(wp), dimension(max(ncplx_w,ncplx_p),n_w,n_p), intent(in) :: pdpsi
    real(wp), dimension(ncplx_w,n_w,ncplx_p,n_p), intent(out) :: scpr
    !local variables
    logical :: cplx_p,cplx_w,cplx_pw
    integer :: iw,ip,icplx

    cplx_p=ncplx_p==2
    cplx_w=ncplx_w==2
    cplx_pw=cplx_p .and. cplx_w

    if (cplx_pw) then
       do ip=1,n_p
          do iw=1,n_w
             scpr(1,iw,1,ip)=pdpsi(1,iw,ip)
             scpr(2,iw,1,ip)=pdpsi(2,iw,ip)
             scpr(1,iw,2,ip)=-pdpsi(2,iw,ip)
             scpr(2,iw,2,ip)=pdpsi(1,iw,ip)
          end do
       end do
       !copy the values, only one of the two might be 2
    else if (cplx_p) then
       do ip=1,n_p
          do icplx=1,ncplx_p
             do iw=1,n_w
                scpr(1,iw,icplx,ip)=pdpsi(icplx,iw,ip)
             end do
          end do
       end do
    else if (cplx_w) then
       do ip=1,n_p
          do iw=1,n_w
             do icplx=1,ncplx_w
                scpr(icplx,iw,1,ip)=pdpsi(icplx,iw,ip)
             end do
          end do
       end do
    else !real case
       do ip=1,n_p
          do iw=1,n_w
             scpr(1,iw,1,ip)=pdpsi(1,iw,ip)
          end do
       end do

    end if
  end subroutine reverse_coefficients

  !> Identify the coefficients
  pure subroutine full_coefficients(trans_p,ncplx_p,n_p,trans_w,ncplx_w,n_w,scpr,trans,pdpsi)
    implicit none
    integer, intent(in) :: ncplx_p,ncplx_w,n_p,n_w
    character(len=1), intent(in) :: trans_p,trans_w,trans
    real(wp), dimension(ncplx_w,n_w,ncplx_p,n_p), intent(in) :: scpr
    real(wp), dimension(max(ncplx_w,ncplx_p),n_w,n_p), intent(out) :: pdpsi
    !local variables
    logical :: cplx_p,cplx_w,cplx_pw
    integer :: iw,ip,ieps_p,ieps_w,ieps
    real(wp) :: prfr,prfi,pifr,pifi

    cplx_p=ncplx_p==2
    cplx_w=ncplx_w==2
    cplx_pw=cplx_p .and. cplx_w

    ieps_p=1
    if (trans_p=='C' .and. cplx_p) ieps_p=-1
    ieps_w=1
    if (trans_w=='C' .and. cplx_w) ieps_w=-1
    ieps=1
    if (trans=='C' .and. (cplx_p .or. cplx_w)) ieps=-1


    !the coefficients have to be transformed to the complex version
    if ((.not. cplx_p) .and. (.not.  cplx_w)) then
       !real case, simply copy the values
       do ip=1,n_p
          do iw=1,n_w
             pdpsi(1,iw,ip)=scpr(1,iw,1,ip)
          end do
       end do
    else
       !complex case, build real and imaginary part when applicable
       prfi=0.0_wp
       pifr=0.0_wp
       pifi=0.0_wp
       do ip=1,n_p
          do iw=1,n_w
             prfr=scpr(1,iw,1,ip)
             if (cplx_p) pifr=scpr(1,iw,2,ip)
             if (cplx_w) prfi=scpr(2,iw,1,ip)
             if (cplx_pw) pifi=scpr(2,iw,2,ip)   
             !real part
             pdpsi(1,iw,ip)=prfr-ieps_p*ieps_w*pifi
             !imaginary part
             pdpsi(2,iw,ip)=ieps*ieps_w*prfi+ieps*ieps_p*pifr
          end do
       end do
    end if

  end subroutine full_coefficients

end subroutine NL_HGH_application

!> Applies one real projector operator in the form |p> hp <p| onto a set of wavefunctions described by the same descriptors
!! accumulate the result on the array hpsi and calculate the energy in the form \sum_w <psi_w|p> hp <p|psi_w>
subroutine apply_oneproj_operator(wfd_p,proj,hp,n_w,wfd_w,psi,hpsi,scpr)
  use module_base
  use module_types, only: wavefunctions_descriptors
  implicit none
  integer, intent(in) :: n_w !< complex components of the wavefunction
  real(wp), intent(in) :: hp !<coefficient of the projector operator
  type(wavefunctions_descriptors), intent(in) :: wfd_p !< descriptors of projectors
  type(wavefunctions_descriptors), intent(in) :: wfd_w !< descriptors of wavefunction
  !  real(gp), dimension(ncplx_o,ncomp_p,ncomp_p,ncomp_w), intent(in) :: hij !< matrix of operator in nonlocal projectors basis
  real(wp), dimension(wfd_p%nvctr_c+7*wfd_p%nvctr_f), intent(in) :: proj !< components of the projector
  real(wp), dimension(wfd_w%nvctr_c+7*wfd_w%nvctr_f,n_w), intent(in) :: psi !< components of wavefunction
  real(wp), dimension(wfd_w%nvctr_c+7*wfd_w%nvctr_f,n_w), intent(inout) :: hpsi !<application of NL operator on psi
  real(wp), dimension(n_w), intent(out) :: scpr !<array of <p|psi_w>, to be used to evaluate energy terms
  !local variables
  character(len=*), parameter :: subname='apply_oneproj'
  integer :: is_w,is_sw,is_p,is_sp,iw
  integer, dimension(:,:), allocatable :: psi_mask
  integer :: proj_count, i_proj
  !routines which are optimized in separate files
  external :: wpdot_keys,wpdot_mask,waxpy_mask

  call f_routine(id=subname)

  !calculate starting points of the fine regions
  !they have to be calculated considering that there could be no fine grid points
  !therefore the array values should not go out of bounds even though their value is actually not used
  is_w=wfd_w%nvctr_c+min(wfd_w%nvctr_f,1)
  is_sw=wfd_w%nseg_c+min(wfd_w%nseg_f,1)

  is_p=wfd_p%nvctr_c+min(wfd_p%nvctr_f,1)
  is_sp=wfd_p%nseg_c+min(wfd_p%nseg_f,1)

  !mask array to avoid multiple calls to bitonic search routines
  psi_mask=f_malloc0((/3,wfd_w%nseg_c+wfd_w%nseg_f/),id='psi_mask')
  call wpdot_keys(wfd_w%nvctr_c,wfd_w%nvctr_f,wfd_w%nseg_c,wfd_w%nseg_f,&
       wfd_w%keyvglob(1),wfd_w%keyvglob(is_sw),wfd_w%keyglob(1,1),wfd_w%keyglob(1,is_sw),&
       psi(1,1),psi(is_w,1),&
       wfd_p%nvctr_c,wfd_p%nvctr_f,wfd_p%nseg_c,wfd_p%nseg_f,&
       wfd_p%keyvglob(1),wfd_p%keyvglob(is_sp),wfd_p%keyglob(1,1),wfd_p%keyglob(1,is_sp),&
       proj(1),proj(is_p),&
       scpr(1))
  !use now mask arrays to calculate the rest of the scalar product
  do iw=2,n_w
  call wpdot_keys(wfd_w%nvctr_c,wfd_w%nvctr_f,wfd_w%nseg_c,wfd_w%nseg_f,&
       wfd_w%keyvglob(1),wfd_w%keyvglob(is_sw),&
       wfd_w%keyglob(1,1),wfd_w%keyglob(1,is_sw),&
       psi(1,iw),psi(is_w,iw),&
       wfd_p%nvctr_c,wfd_p%nvctr_f,wfd_p%nseg_c,wfd_p%nseg_f,&
       wfd_p%keyvglob(1),wfd_p%keyvglob(is_sp),&
       wfd_p%keyglob(1,1),wfd_p%keyglob(1,is_sp),&
       proj(1),proj(is_p),&
       scpr(iw))

!!$     call wpdot_mask(wfd_w%nvctr_c,wfd_w%nvctr_f,wfd_w%nseg_c,wfd_w%nseg_f,&
!!$          psi_mask(1,1),psi_mask(1,is_sw),psi(1,iw),psi(is_w,iw),&
!!$          wfd_p%nvctr_c,wfd_p%nvctr_f,proj(1),proj(is_p),scpr(iw))
  end do

  !then reduce the projector in the wavefunction
  do iw=1,n_w
     call waxpy(hp*scpr(iw),wfd_p%nvctr_c,wfd_p%nvctr_f,&
          wfd_p%nseg_c,wfd_p%nseg_f,&
          wfd_p%keyvglob(1),wfd_p%keyvglob(is_sp),&
          wfd_p%keyglob(1,1),wfd_p%keyglob(1,is_sp),proj(1),proj(is_p),& 
          wfd_w%nvctr_c,wfd_w%nvctr_f,wfd_w%nseg_c,wfd_w%nseg_f,&
          wfd_w%keyvglob(1),wfd_w%keyvglob(is_sw),&
          wfd_w%keyglob(1,1),wfd_w%keyglob(1,is_sw),&
          hpsi(1,iw),hpsi(is_w,iw))
!!$     call waxpy_mask(wfd_w%nvctr_c,wfd_w%nvctr_f,wfd_w%nseg_c,wfd_w%nseg_f,&
!!$          psi_mask(1,1),psi_mask(1,is_sw),hpsi(1,iw),hpsi(is_w,iw),&
!!$          wfd_p%nvctr_c,wfd_p%nvctr_f,proj(1),proj(is_p),&
!!$          hp*scpr(iw))
  end do

  call f_free(psi_mask)

  call f_release_routine()

end subroutine apply_oneproj_operator


!> Performs the scalar product of a projector with a wavefunction each one writeen in Daubechies basis
!! with its own descriptors.
!! A masking array is then calculated to avoid the calculation of bitonic search for the scalar product
!! If the number of projectors is bigger than 1 the wavefunction is also packed in the number of components
!! of the projector to ease its successive application
subroutine proj_dot_psi(n_p,wfd_p,proj,n_w,wfd_w,psi,nmseg_c,nmseg_f,psi_mask,psi_pack,scpr)
  use module_base, only: wp,gemm
  use module_types, only: wavefunctions_descriptors
  implicit none
  integer, intent(in) :: n_p !< number of projectors (real and imaginary part included)
  integer, intent(in) :: n_w !< number of wavefunctions (real and imaginary part included)
  integer, intent(in) :: nmseg_c,nmseg_f !< segments of the masking array
  type(wavefunctions_descriptors), intent(in) :: wfd_p !< descriptors of projectors
  type(wavefunctions_descriptors), intent(in) :: wfd_w !< descriptors of wavefunction
  real(wp), dimension(wfd_p%nvctr_c+7*wfd_p%nvctr_f,n_p), intent(in) :: proj !< components of the projectors
  real(wp), dimension(wfd_w%nvctr_c+7*wfd_w%nvctr_f,n_w), intent(in) :: psi !< components of wavefunction
  integer, dimension(3,nmseg_c+nmseg_f), intent(in) :: psi_mask !<lookup array in the wfn segments
  !indicating the points where data have to be taken for dot product
  ! always produced. Has to be initialized to zero first
  real(wp), dimension(wfd_p%nvctr_c+7*wfd_p%nvctr_f,n_w), intent(inout) :: psi_pack !< packed array of psi in projector form
  !needed only when n_p is bigger than one 
  real(wp), dimension(n_w,n_p), intent(out) :: scpr !< array of the scalar product of all the components
  !local variables
  logical, parameter :: mask=.true.,pack=.true.
  integer :: is_w,is_sw,is_p,is_sp,iw,ip,is_sm
  !intensive routines
  external :: wpdot_keys_pack,wpdot_mask_pack

  !calculate starting points of the fine regions
  !they have to be calculated considering that there could be no fine grid points
  !therefore the array values should not go out of bounds even though their value is actually not used
  is_w=wfd_w%nvctr_c+min(wfd_w%nvctr_f,1)
  is_sw=wfd_w%nseg_c+min(wfd_w%nseg_f,1)

  is_p=wfd_p%nvctr_c+min(wfd_p%nvctr_f,1)
  is_sp=wfd_p%nseg_c+min(wfd_p%nseg_f,1)

  is_sm=nmseg_c+min(nmseg_f,1)

  if (pack) then
     if (.not. mask) then
        do iw=1,n_w
           call wpdot_keys_pack(wfd_w%nvctr_c,wfd_w%nvctr_f,wfd_w%nseg_c,wfd_w%nseg_f,&
                wfd_w%keyvglob(1),wfd_w%keyvglob(is_sw),wfd_w%keyglob(1,1),wfd_w%keyglob(1,is_sw),&
                psi(1,iw),psi(is_w,iw),&
                wfd_p%nvctr_c,wfd_p%nvctr_f,wfd_p%nseg_c,wfd_p%nseg_f,&
                wfd_p%keyvglob(1),wfd_p%keyvglob(is_sp),wfd_p%keyglob(1,1),wfd_p%keyglob(1,is_sp),&
                proj(1,1),proj(is_p,1),&
                psi_pack(1,iw),psi_pack(is_p,iw),scpr(iw,1))
        end do
     else 
        do iw=1,n_w
           call wpdot_mask_pack(wfd_w%nvctr_c,wfd_w%nvctr_f,nmseg_c,nmseg_f,&
                psi_mask(1,1),psi_mask(1,is_sm),psi(1,iw),psi(is_w,iw),&
                wfd_p%nvctr_c,wfd_p%nvctr_f,proj(1,1),proj(is_p,1),&
                psi_pack(1,iw),psi_pack(is_p,iw),scpr(iw,1))
        end do
     end if

     !now that the packed array is constructed linear algebra routine can be used to calculate
     !use multithreaded dgemm or customized ones in the case of no OMP parallelized algebra
     !scpr(iw,ip) = < psi_iw| p_ip >
     if (n_p > 1) then
        call gemm('T','N',n_w,n_p-1,wfd_p%nvctr_c+7*wfd_p%nvctr_f,1.0_wp,psi_pack(1,1),&
             wfd_p%nvctr_c+7*wfd_p%nvctr_f,proj(1,2),wfd_p%nvctr_c+7*wfd_p%nvctr_f,0.0_wp,&
             scpr(1,2),n_w)
     end if

  else
     do ip=1,n_p
        do iw=1,n_w
           call wpdot_keys(wfd_w%nvctr_c,wfd_w%nvctr_f,wfd_w%nseg_c,wfd_w%nseg_f,&
                wfd_w%keyvglob(1),wfd_w%keyvglob(is_sw),wfd_w%keyglob(1,1),wfd_w%keyglob(1,is_sw),&
                psi(1,iw),psi(is_w,iw),&
                wfd_p%nvctr_c,wfd_p%nvctr_f,wfd_p%nseg_c,wfd_p%nseg_f,&
                wfd_p%keyvglob(1),wfd_p%keyvglob(is_sp),wfd_p%keyglob(1,1),wfd_p%keyglob(1,is_sp),&
                proj(1,ip),proj(is_p,ip),&
                scpr(iw,ip))
        end do
     end do
  end if
end subroutine proj_dot_psi

!> Performs the update of a set of wavefunctions with a projector each one written in Daubechies basis
!! with its own descriptors.
!! A masking array is used calculated to avoid the calculation of bitonic search for the scalar product
!! If the number of projectors is bigger than 1 the wavefunction is also given by packing in the number of components
!! of the projector to ease its successive application
subroutine scpr_proj_p_hpsi(n_p,wfd_p,proj,n_w,wfd_w,&
     nmseg_c,nmseg_f,psi_mask,hpsi_pack,scpr,hpsi)
  use module_base, only: wp,gemm,to_zero
  use module_types, only: wavefunctions_descriptors
  implicit none
  integer, intent(in) :: n_p !< number of projectors (real and imaginary part included)
  integer, intent(in) :: n_w !< number of wavefunctions (real and imaginary part included)
  integer, intent(in) :: nmseg_c,nmseg_f !< segments of the masking array
  type(wavefunctions_descriptors), intent(in) :: wfd_p !< descriptors of projectors
  type(wavefunctions_descriptors), intent(in) :: wfd_w !< descriptors of wavefunction
  real(wp), dimension(n_w,n_p), intent(in) :: scpr !< array of the scalar product of all the components
  real(wp), dimension(wfd_p%nvctr_c+7*wfd_p%nvctr_f,n_p), intent(in) :: proj !< components of the projectors
  integer, dimension(3,nmseg_c+nmseg_f), intent(in) :: psi_mask !<lookup array in the wfn segments
  !indicating the points where data have to be taken for dot product
  ! always produced. Has to be initialized to zero first
  real(wp), dimension(wfd_p%nvctr_c+7*wfd_p%nvctr_f,n_w), intent(inout) :: hpsi_pack !< work array of hpsi in projector form
  !needed only when n_p is bigger than one 

  real(wp), dimension(wfd_w%nvctr_c+7*wfd_w%nvctr_f,n_w), intent(inout) :: hpsi !< wavefunction result
  !local variables
  logical, parameter :: mask=.false.,pack=.true.
  external :: waxpy_mask_unpack
  integer :: is_w,is_sw,is_p,is_sp,iw,is_sm

  is_w=wfd_w%nvctr_c+min(wfd_w%nvctr_f,1)
  is_sw=wfd_w%nseg_c+min(wfd_w%nseg_f,1)

  is_p=wfd_p%nvctr_c+min(wfd_p%nvctr_f,1)
  is_sp=wfd_p%nseg_c+min(wfd_p%nseg_f,1)

  is_sm=nmseg_c+min(nmseg_f,1)

  if (pack) then
     !once the coefficients are determined fill the components of the wavefunction with the last projector
     !linear algebra up to the second last projector
     !|psi_iw>=O_iw,jp| p_jp>
     
     if (n_p > 1) then
        call gemm('N','T',wfd_p%nvctr_c+7*wfd_p%nvctr_f,n_w,n_p-1,&
             1.0_wp,proj(1,1),wfd_p%nvctr_c+7*wfd_p%nvctr_f,&
             scpr(1,1),n_w,0.0_wp,&
             hpsi_pack(1,1),wfd_p%nvctr_c+7*wfd_p%nvctr_f)
     else
        call to_zero(n_w*(wfd_p%nvctr_c+7*wfd_p%nvctr_f),hpsi_pack(1,1))
     end if

     !then last projector
     if (mask) then
        do iw=1,n_w
           call waxpy_mask_unpack(wfd_w%nvctr_c,wfd_w%nvctr_f,nmseg_c,nmseg_f,&
                psi_mask(1,1),psi_mask(1,is_sm),hpsi_pack(1,iw),hpsi_pack(is_p,iw),&
                hpsi(1,iw),hpsi(is_w,iw),&
                wfd_p%nvctr_c,wfd_p%nvctr_f,proj(1,n_p),proj(is_p,n_p),&
                scpr(iw,n_p))
        end do
     else
        do iw=1,n_w
           call waxpy_keys_unpack(wfd_w%nvctr_c,wfd_w%nvctr_f,wfd_w%nseg_c,wfd_w%nseg_f,&
                wfd_w%keyvglob(1),wfd_w%keyvglob(is_sw),wfd_w%keyglob(1,1),wfd_w%keyglob(1,is_sw),&
                hpsi(1,iw),hpsi(is_w,iw),&
                wfd_p%nvctr_c,wfd_p%nvctr_f,wfd_p%nseg_c,wfd_p%nseg_f,&
                wfd_p%keyvglob(1),wfd_p%keyvglob(is_sp),&
                wfd_p%keyglob(1,1),wfd_p%keyglob(1,is_sp),&
                proj(1,n_p),proj(is_p,n_p),&
                hpsi_pack(1,iw),hpsi_pack(is_p,iw),scpr(iw,n_p))
        end do
     end if
  else
     
  end if

end subroutine scpr_proj_p_hpsi

!> Calculates the dot product between a wavefunctions apsi and a projector bpsi (both in compressed form)
!! The array mask is then constructed so that successive application of the projector on the same object can
!! be done without bitonic search
subroutine wpdot_keys(  &
     mavctr_c,mavctr_f,maseg_c,maseg_f,keyav_c,keyav_f,keyag_c,keyag_f,apsi_c,apsi_f,&
     mbvctr_c,mbvctr_f,mbseg_c,mbseg_f,keybv_c,keybv_f,keybg_c,keybg_f,bpsi_c,bpsi_f,&
     scpr)
  use module_base, only: wp
  implicit none

  real(wp), dimension(mavctr_c), intent(in) :: apsi_c
  real(wp), dimension(7,mavctr_f), intent(in) :: apsi_f
  real(wp), intent(out) :: scpr
  integer, intent(in) :: mavctr_c,mavctr_f,maseg_c,maseg_f,mbvctr_c,mbvctr_f,mbseg_c,mbseg_f
  integer, dimension(maseg_c), intent(in) :: keyav_c
  integer, dimension(maseg_f), intent(in) :: keyav_f
  integer, dimension(mbseg_c), intent(in) :: keybv_c
  integer, dimension(mbseg_f), intent(in) :: keybv_f
  integer, dimension(2,maseg_c), intent(in) :: keyag_c
  integer, dimension(2,maseg_f), intent(in) :: keyag_f
  integer, dimension(2,mbseg_c), intent(in) :: keybg_c
  integer, dimension(2,mbseg_f), intent(in) :: keybg_f
  real(wp), dimension(mbvctr_c), intent(in) :: bpsi_c
  real(wp), dimension(7,mbvctr_f), intent(in) :: bpsi_f
  !local variables
  integer :: ibseg,jaj,jb1,jb0,jbj,iaoff,iboff,length,ja0,ja1,i,j
  real(wp) :: scpr1,scpr0,tt
  integer :: iaseg0,ibsegs,ibsege
  !these arrays have to be allocatable
  integer, dimension(maseg_c) :: keyag_c_lin !>linear version of second indices of keyag_c
  integer, dimension(maseg_f) :: keyag_f_lin !>linear version of second indices of keyag_f
  !Variables for OpenMP
  !$ integer :: ithread,nthread,nchunk
  !$ integer :: omp_get_thread_num,omp_get_num_threads

  keyag_c_lin = keyag_c(1,:) !speed up access in hunt subroutine by consecutive arrangement in memory
  keyag_f_lin = keyag_f(1,:) !speed up access in hunt subroutine by consecutive arrangement in memory

  scpr=0.0_wp

  !$omp parallel default (none) &
  !$omp shared (maseg_c,keyav_c,keyag_c,keyag_c_lin,keybg_c,mbseg_c,keybv_c,mbseg_f,maseg_f)&
  !$omp shared (apsi_c,bpsi_c,bpsi_f,keybv_f,keybg_f,keyag_f,keyag_f_lin,keyav_f)&
  !$omp shared (apsi_f,scpr) &
!!$  !$omp parallel default(shared) &
  !$omp private(jaj,iaoff,length,ja1,ja0,jb1,jb0,iboff,scpr0,scpr1) &
  !$omp private(jbj,ibseg,iaseg0,i,j,tt)!!!,ithread,nthread,ibsegs,ibsege,nchunk)

  scpr0=0.0_wp
  scpr1=0.0_wp

!!!!start of general region

  !alternative way of parallelizing the loop, to be tested to explore performances
  !LG  ibsegs=1
  !LG  ibsege=mbseg_c
  !LG  !$ ithread=omp_get_thread_num()
  !LG  !$ nthread=omp_get_num_threads() 

  iaseg0=1 

  !coarse part. Loop on the projectors segments
  !LG  !separate in chunks the loop among the threads
  !LG  !$ nchunck=max(mbseg_c/nthread,1)
  !LG  !$ ibsegs=min(ithread*nchunck,mbseg_c+1)
  !LG  !$ ibsege=min((ithread+1)*nchunck,mbseg_c)

  !$omp do schedule(static)
  do ibseg=1,mbseg_c
     jbj=keybv_c(ibseg)
     !     jb0=keybg_c(1,ibseg) !starting point of projector segment
     jb0=max(keybg_c(1,ibseg),keyag_c_lin(1))
     jb1=keybg_c(2,ibseg) !ending point of projector segment
     iboff = max(jb0-keybg_c(1,ibseg),0)

     !find the starting point of the wavefunction segment
     !warning: hunt is assuming that the variable is always found
     !if it is not, iaseg0 is put to maseg + 1 so that the loop is disabled
     call hunt_inline(keyag_c_lin,maseg_c,jb0,iaseg0)
     if (iaseg0==0) then  !segment not belonging to the wavefunctions, go further
        iaseg0=1
        cycle     
     end if
     !now pass through all the wavefunction segments until the end of the segment is 
     !still contained in projector segment
     nonconvex_loop_c: do while(iaseg0 <= maseg_c)
        !length = jb1-jb0
        !iaoff = jb0-keyag_c_lin(iaseg0)!jb0-ja0

        ja0=keyag_c_lin(iaseg0)
        ja1=min(jb1,keyag_c(2,iaseg0)) 
        length = ja1-jb0
        iaoff = max(jb0-ja0,0) !no offset if we are already inside

        jaj=keyav_c(iaseg0)

        do i=0,length
           tt=apsi_c(jaj+iaoff+i)
           scpr0=scpr0+tt*bpsi_c(jbj+i+iboff)
        enddo

        !call op_c_inline(length,jaj+iaoff,jbj+iboff,scpr0)

        if ((ja1<=jb1 .and. length>=0) .or. iaseg0==maseg_c) exit nonconvex_loop_c !segment is finished
        iaseg0=iaseg0+1
        jb0=max(keybg_c(1,ibseg),keyag_c_lin(iaseg0))
        if (keyag_c_lin(iaseg0)>jb1) exit nonconvex_loop_c !segment is not covered
        jbj=jbj+max(jb0-keybg_c(1,ibseg),0)
     end do nonconvex_loop_c
     !disable loop if the end is reached
     if (iaseg0 == maseg_c .and. keybg_c(1,ibseg)> keyag_c_lin(maseg_c)) iaseg0=iaseg0+1
  enddo
  !stop
  !$omp end do nowait

  ! fine part
  !LG  ibsegs=1
  !LG  ibsege=mbseg_f

  iaseg0=1

  !LG  !separate in chunks the loop among the threads
  !LG  !$ nchunck=max(mbseg_f/nthread,1)
  !LG  !$ ibsegs=min(ithread*nchunck,mbseg_f+1)
  !LG  !$ ibsege=min((ithread+1)*nchunck,mbseg_f)

  !$omp do schedule(static)
  do ibseg=1,mbseg_f
     jbj=keybv_f(ibseg)
     !jb0=keybg_f(1,ibseg)
     jb0=max(keybg_f(1,ibseg),keyag_f_lin(1))
     jb1=keybg_f(2,ibseg)
     iboff = max(jb0-keybg_f(1,ibseg),0)
     call hunt_inline(keyag_f_lin,maseg_f,jb0,iaseg0)
     if (iaseg0==0) then  !segment not belonging to the wavefunctions, go further
        iaseg0=1
        cycle     
     end if
     nonconvex_loop_f: do while(iaseg0 <= maseg_f)
!!$     length = jb1-jb0
!!$     iaoff = jb0-keyag_f_lin(iaseg0)

        ja0=keyag_f_lin(iaseg0) !still doubts about copying in automatic array
        ja1=min(jb1,keyag_f(2,iaseg0)) 
        length = ja1-jb0
        iaoff = max(jb0-ja0,0) !no offset if we are already inside
        jaj=keyav_f(iaseg0)

        do i=0,length
           do j=1,7
              tt=apsi_f(j,jaj+iaoff+i)
              scpr1=scpr1+tt*bpsi_f(j,jbj+i+iboff)
           end do
        enddo
        !call op_f_inline(length,jaj+iaoff,jbj+iboff,scpr1)

        if ((ja1<=jb1 .and. length>=0) .or. iaseg0==maseg_f) exit nonconvex_loop_f !segment is finished  
        iaseg0=iaseg0+1
        jb0=max(keybg_f(1,ibseg),keyag_f_lin(iaseg0))
        if (keyag_f_lin(iaseg0)>jb1) exit nonconvex_loop_f !segment is not covered 
        jbj=jbj+max(jb0-keybg_f(1,ibseg),0)
     end do nonconvex_loop_f
     !disable loop if the end is reached
     if (iaseg0 == maseg_f .and. keybg_f(1,ibseg)> keyag_f_lin(maseg_f)) iaseg0=iaseg0+1
  enddo
  !$omp end do !implicit barrier 

!!!!end of general region

  scpr0=scpr0+scpr1

  !$omp critical 
  scpr=scpr+scpr0
  !$omp end critical

  !$omp end parallel

!  include 'wpdot-inc.f90'

contains

!!$  pure subroutine op_c_inline(length,jaj,jbj,scpr0)
!!$    implicit none
!!$    integer, intent(in) :: length,jaj,jbj
!!$    real(wp), intent(inout) :: scpr0
!!$    !local variables
!!$    integer :: i
!!$    real(wp) :: tt
!!$    do i=0,length
!!$       tt=apsi_c(jaj+i)!iaoff+i)
!!$       scpr0=scpr0+tt*bpsi_c(jbj+i)!+iboff)
!!$    enddo
!!$  end subroutine op_c_inline
!!$
!!$  pure subroutine op_f_inline(length,jaj,jbj,scpr1)
!!$    implicit none
!!$    integer, intent(in) :: length,jaj,jbj
!!$    real(wp), intent(inout) :: scpr1
!!$    !local variables
!!$    integer :: i,j
!!$    real(wp) :: tt
!!$    do i=0,length
!!$       do j=1,7
!!$          tt=apsi_f(j,jaj+i)!iaoff+i)
!!$          scpr1=scpr1+tt*bpsi_f(j,jbj+i)!+iboff)
!!$       end do
!!$    enddo
!!$  end subroutine op_f_inline

  include 'scalar_product-inc.f90'
  
END SUBROUTINE wpdot_keys

!> Calculates the dot product between a wavefunctions apsi and a projector bpsi (both in compressed form)
!! The array mask is then constructed so that successive application of the projector on the same object can
!! be done without bitonic search
!! The array of the wavefunction is also compressed
!! so that successive projector application can be performed also with linear algebra routines
subroutine wpdot_keys_pack(  &
     mavctr_c,mavctr_f,maseg_c,maseg_f,keyav_c,keyav_f,keyag_c,keyag_f,apsi_c,apsi_f,  &
     mbvctr_c,mbvctr_f,mbseg_c,mbseg_f,keybv_c,keybv_f,keybg_c,keybg_f,bpsi_c,bpsi_f,&
     apack_c,apack_f,scpr)
  use module_base, only: wp
  implicit none
  real(wp), dimension(mavctr_c), intent(in) :: apsi_c
  real(wp), dimension(7,mavctr_f), intent(in) :: apsi_f
  real(wp), intent(out) :: scpr
  real(wp), dimension(mbvctr_c), intent(inout) :: apack_c
  real(wp), dimension(7,mbvctr_f), intent(inout) :: apack_f

  integer, intent(in) :: mavctr_c,mavctr_f,maseg_c,maseg_f,mbvctr_c,mbvctr_f,mbseg_c,mbseg_f
  integer, dimension(maseg_c), intent(in) :: keyav_c
  integer, dimension(maseg_f), intent(in) :: keyav_f
  integer, dimension(mbseg_c), intent(in) :: keybv_c
  integer, dimension(mbseg_f), intent(in) :: keybv_f
  integer, dimension(2,maseg_c), intent(in) :: keyag_c
  integer, dimension(2,maseg_f), intent(in) :: keyag_f
  integer, dimension(2,mbseg_c), intent(in) :: keybg_c
  integer, dimension(2,mbseg_f), intent(in) :: keybg_f
  real(wp), dimension(mbvctr_c), intent(in) :: bpsi_c
  real(wp), dimension(7,mbvctr_f), intent(in) :: bpsi_f
  !local variables
  integer :: ibseg,jaj,jb1,jb0,jbj,iaoff,iboff,length,ja0,ja1,i,j
  real(wp) :: scpr1,scpr0,tt
  integer :: iaseg0,ibsegs,ibsege
  !these arrays have to be allocatable
  integer, dimension(maseg_c) :: keyag_c_lin !>linear version of second indices of keyag_c
  integer, dimension(maseg_f) :: keyag_f_lin !>linear version of second indices of keyag_f
  !Variables for OpenMP
  !$ integer :: ithread,nthread,nchunk
  !$ integer :: omp_get_thread_num,omp_get_num_threads

  keyag_c_lin = keyag_c(1,:) !speed up access in hunt subroutine by consecutive arrangement in memory
  keyag_f_lin = keyag_f(1,:) !speed up access in hunt subroutine by consecutive arrangement in memory

  scpr=0.0_wp

  !$omp parallel default (none) &
  !$omp shared (maseg_c,keyav_c,keyag_c,keyag_c_lin,keybg_c,mbseg_c,keybv_c,mbseg_f,maseg_f)&
  !$omp shared (apsi_c,bpsi_c,bpsi_f,keybv_f,keybg_f,keyag_f,keyag_f_lin,keyav_f)&
  !$omp shared (apsi_f,scpr,apack_c,apack_f) &
!!$  !$omp parallel default(shared) &
  !$omp private(jaj,iaoff,length,ja1,ja0,jb1,jb0,iboff,scpr0,scpr1) &
  !$omp private(jbj,ibseg,iaseg0,i,j,tt)!!!,ithread,nthread,ibsegs,ibsege,nchunk)

  scpr0=0.0_wp
  scpr1=0.0_wp

!!!!start of general region

  !alternative way of parallelizing the loop, to be tested to explore performances
  !LG  ibsegs=1
  !LG  ibsege=mbseg_c
  !LG  !$ ithread=omp_get_thread_num()
  !LG  !$ nthread=omp_get_num_threads() 

  iaseg0=1 

  !coarse part. Loop on the projectors segments
  !LG  !separate in chunks the loop among the threads
  !LG  !$ nchunck=max(mbseg_c/nthread,1)
  !LG  !$ ibsegs=min(ithread*nchunck,mbseg_c+1)
  !LG  !$ ibsege=min((ithread+1)*nchunck,mbseg_c)

  !$omp do schedule(static)
  do ibseg=1,mbseg_c
     jbj=keybv_c(ibseg)
     !     jb0=keybg_c(1,ibseg) !starting point of projector segment
     jb0=max(keybg_c(1,ibseg),keyag_c_lin(1))
     jb1=keybg_c(2,ibseg) !ending point of projector segment
     iboff = max(jb0-keybg_c(1,ibseg),0)

     !find the starting point of the wavefunction segment
     !warning: hunt is assuming that the variable is always found
     !if it is not, iaseg0 is put to maseg + 1 so that the loop is disabled
     call hunt_inline(keyag_c_lin,maseg_c,jb0,iaseg0)
     if (iaseg0==0) then  !segment not belonging to the wavefunctions, go further
        iaseg0=1
        cycle     
     end if
     !now pass through all the wavefunction segments until the end of the segment is 
     !still contained in projector segment
     nonconvex_loop_c: do while(iaseg0 <= maseg_c)
        !length = jb1-jb0
        !iaoff = jb0-keyag_c_lin(iaseg0)!jb0-ja0

        ja0=keyag_c_lin(iaseg0)
        ja1=min(jb1,keyag_c(2,iaseg0)) 
        length = ja1-jb0
        iaoff = max(jb0-ja0,0) !no offset if we are already inside

        jaj=keyav_c(iaseg0)

        do i=0,length
           tt=apsi_c(jaj+iaoff+i)
           scpr0=scpr0+tt*bpsi_c(jbj+i+iboff)
           apack_c(jbj+i+iboff)=tt
           !apack_c(jbj+i)=tt
        enddo
        !call op_c_inline(length,jaj+iaoff,jbj+iboff,scpr0)

        if ((ja1<=jb1 .and. length>=0) .or. iaseg0==maseg_c) exit nonconvex_loop_c !segment is finished
        iaseg0=iaseg0+1
        jb0=max(keybg_c(1,ibseg),keyag_c_lin(iaseg0))
        if (keyag_c_lin(iaseg0)>jb1) exit nonconvex_loop_c !segment is not covered
        jbj=jbj+max(jb0-keybg_c(1,ibseg),0)
     end do nonconvex_loop_c
     !disable loop if the end is reached
     if (iaseg0 == maseg_c .and. keybg_c(1,ibseg)> keyag_c_lin(maseg_c)) iaseg0=iaseg0+1
  enddo
  !stop
  !$omp end do nowait

  ! fine part
  !LG  ibsegs=1
  !LG  ibsege=mbseg_f

  iaseg0=1

  !LG  !separate in chunks the loop among the threads
  !LG  !$ nchunck=max(mbseg_f/nthread,1)
  !LG  !$ ibsegs=min(ithread*nchunck,mbseg_f+1)
  !LG  !$ ibsege=min((ithread+1)*nchunck,mbseg_f)

  !$omp do schedule(static)
  do ibseg=1,mbseg_f
     jbj=keybv_f(ibseg)
     !jb0=keybg_f(1,ibseg)
     jb0=max(keybg_f(1,ibseg),keyag_f_lin(1))
     jb1=keybg_f(2,ibseg)
     iboff = max(jb0-keybg_f(1,ibseg),0)
     call hunt_inline(keyag_f_lin,maseg_f,jb0,iaseg0)
     if (iaseg0==0) then  !segment not belonging to the wavefunctions, go further
        iaseg0=1
        cycle     
     end if
     nonconvex_loop_f: do while(iaseg0 <= maseg_f)
!!$     length = jb1-jb0
!!$     iaoff = jb0-keyag_f_lin(iaseg0)

        ja0=keyag_f_lin(iaseg0) !still doubts about copying in automatic array
        ja1=min(jb1,keyag_f(2,iaseg0)) 
        length = ja1-jb0
        iaoff = max(jb0-ja0,0) !no offset if we are already inside
        jaj=keyav_f(iaseg0)

        do i=0,length
           do j=1,7
              tt=apsi_f(j,jaj+iaoff+i)
              scpr1=scpr1+tt*bpsi_f(j,jbj+i+iboff)
              apack_f(j,jbj+i+iboff)=tt
              !apack_f(j,jbj+i)=tt
           end do
        enddo
        !call op_f_inline(length,jaj+iaoff,jbj+iboff,scpr1)

        if ((ja1<=jb1 .and. length>=0) .or. iaseg0==maseg_f) exit nonconvex_loop_f !segment is finished  
        iaseg0=iaseg0+1
        jb0=max(keybg_f(1,ibseg),keyag_f_lin(iaseg0))
        if (keyag_f_lin(iaseg0)>jb1) exit nonconvex_loop_f !segment is not covered 
        jbj=jbj+max(jb0-keybg_f(1,ibseg),0)
     end do nonconvex_loop_f
     !disable loop if the end is reached
     if (iaseg0 == maseg_f .and. keybg_f(1,ibseg)> keyag_f_lin(maseg_f)) iaseg0=iaseg0+1
  enddo
  !$omp end do !implicit barrier 

!!!!end of general region

  scpr0=scpr0+scpr1

  !$omp critical 
  scpr=scpr+scpr0
  !$omp end critical

  !$omp end parallel

!  include 'wpdot-inc.f90'

contains

!!$  subroutine op_c_inline(length,jaj,jbj,scpr0)
!!$    implicit none
!!$    integer, intent(in) :: length,jaj,jbj
!!$    real(wp), intent(inout) :: scpr0
!!$    !local variables
!!$    integer :: i
!!$    real(wp) :: tt
!!$
!!$    do i=0,length
!!$       tt=apsi_c(jaj+i)!iaoff+i)
!!$       scpr0=scpr0+tt*bpsi_c(jbj+i)!+iboff)
!!$       !apack_c(jbj+i+iboff)=tt
!!$       apack_c(jbj+i)=tt
!!$    enddo
!!$  end subroutine op_c_inline
!!$
!!$  subroutine op_f_inline(length,jaj,jbj,scpr1)
!!$    implicit none
!!$    integer, intent(in) :: length,jaj,jbj
!!$    real(wp), intent(inout) :: scpr1
!!$    !local variables
!!$    integer :: i,j
!!$    real(wp) :: tt
!!$
!!$    do i=0,length
!!$       do j=1,7
!!$          tt=apsi_f(j,jaj+i)!iaoff+i)
!!$          scpr1=scpr1+tt*bpsi_f(j,jbj+i)!+iboff)
!!$          !apack_f(j,jbj+i+iboff)=tt
!!$          apack_f(j,jbj+i)=tt
!!$       end do
!!$    enddo
!!$  end subroutine op_f_inline

  include 'scalar_product-inc.f90'
  
END SUBROUTINE wpdot_keys_pack

!> Calculates the dot product between a wavefunctions apsi and a projector bpsi (both in compressed form)
!! The array mask is then constructed so that successive application of the projector on the same object can
!! be done without bitonic search
!! The array of the wavefunction is also compressed
!! so that successive projector application can be performed also with linear algebra routines
subroutine waxpy_keys_unpack(  &
     mavctr_c,mavctr_f,maseg_c,maseg_f,keyav_c,keyav_f,keyag_c,keyag_f,apsi_c,apsi_f,  &
     mbvctr_c,mbvctr_f,mbseg_c,mbseg_f,keybv_c,keybv_f,keybg_c,keybg_f,bpsi_c,bpsi_f,&
     apack_c,apack_f,scpr)
  use module_base, only: wp
  implicit none
  real(wp), dimension(mavctr_c), intent(inout) :: apsi_c
  real(wp), dimension(7,mavctr_f), intent(inout) :: apsi_f
  real(wp), intent(in) :: scpr
  real(wp), dimension(mbvctr_c), intent(in) :: apack_c
  real(wp), dimension(7,mbvctr_f), intent(in) :: apack_f
  integer, intent(in) :: mavctr_c,mavctr_f,maseg_c,maseg_f,mbvctr_c,mbvctr_f,mbseg_c,mbseg_f
  integer, dimension(maseg_c), intent(in) :: keyav_c
  integer, dimension(maseg_f), intent(in) :: keyav_f
  integer, dimension(mbseg_c), intent(in) :: keybv_c
  integer, dimension(mbseg_f), intent(in) :: keybv_f
  integer, dimension(2,maseg_c), intent(in) :: keyag_c
  integer, dimension(2,maseg_f), intent(in) :: keyag_f
  integer, dimension(2,mbseg_c), intent(in) :: keybg_c
  integer, dimension(2,mbseg_f), intent(in) :: keybg_f
  real(wp), dimension(mbvctr_c), intent(in) :: bpsi_c
  real(wp), dimension(7,mbvctr_f), intent(in) :: bpsi_f
  !local variables
  integer :: ibseg,jaj,jb1,jb0,jbj,iaoff,iboff,length,ja0,ja1,i,j
  integer :: iaseg0,ibsegs,ibsege
  real(wp) :: tt
  !these arrays have to be allocatable
  integer, dimension(maseg_c) :: keyag_c_lin !>linear version of second indices of keyag_c
  integer, dimension(maseg_f) :: keyag_f_lin !>linear version of second indices of keyag_f
  !Variables for OpenMP
  !$ integer :: ithread,nthread,nchunk
  !$ integer :: omp_get_thread_num,omp_get_num_threads

  keyag_c_lin = keyag_c(1,:) !speed up access in hunt subroutine by consecutive arrangement in memory
  keyag_f_lin = keyag_f(1,:) !speed up access in hunt subroutine by consecutive arrangement in memory

  !$omp parallel default (none) &
  !$omp shared (maseg_c,keyav_c,keyag_c,keyag_c_lin,keybg_c,mbseg_c,keybv_c,mbseg_f,maseg_f)&
  !$omp shared (apsi_c,bpsi_c,bpsi_f,keybv_f,keybg_f,keyag_f,keyag_f_lin,keyav_f)&
  !$omp shared (apsi_f,scpr,apack_c,apack_f) &
!!$  !$omp parallel default(shared) &
  !$omp private(jaj,iaoff,length,ja1,ja0,jb1,jb0,iboff,i,j,tt) &
  !$omp private(jbj,ibseg,iaseg0)!!!,ithread,nthread,ibsegs,ibsege,nchunk)

  !alternative way of parallelizing the loop, to be tested to explore performances
  !LG  ibsegs=1
  !LG  ibsege=mbseg_c
  !LG  !$ ithread=omp_get_thread_num()
  !LG  !$ nthread=omp_get_num_threads() 

  iaseg0=1 

  !coarse part. Loop on the projectors segments
  !LG  !separate in chunks the loop among the threads
  !LG  !$ nchunck=max(mbseg_c/nthread,1)
  !LG  !$ ibsegs=min(ithread*nchunck,mbseg_c+1)
  !LG  !$ ibsege=min((ithread+1)*nchunck,mbseg_c)

  !$omp do schedule(static)
  do ibseg=1,mbseg_c
     jbj=keybv_c(ibseg)
     !     jb0=keybg_c(1,ibseg) !starting point of projector segment
     jb0=max(keybg_c(1,ibseg),keyag_c_lin(1))
     jb1=keybg_c(2,ibseg) !ending point of projector segment
     iboff = max(jb0-keybg_c(1,ibseg),0)

     !find the starting point of the wavefunction segment
     !warning: hunt is assuming that the variable is always found
     !if it is not, iaseg0 is put to maseg + 1 so that the loop is disabled
     call hunt_inline(keyag_c_lin,maseg_c,jb0,iaseg0)
     if (iaseg0==0) then  !segment not belonging to the wavefunctions, go further
        iaseg0=1
        cycle     
     end if
     !now pass through all the wavefunction segments until the end of the segment is 
     !still contained in projector segment
     nonconvex_loop_c: do while(iaseg0 <= maseg_c)
        !length = jb1-jb0
        !iaoff = jb0-keyag_c_lin(iaseg0)!jb0-ja0

        ja0=keyag_c_lin(iaseg0)
        ja1=min(jb1,keyag_c(2,iaseg0)) 
        length = ja1-jb0
        iaoff = max(jb0-ja0,0) !no offset if we are already inside

        jaj=keyav_c(iaseg0)

        do i=0,length
           !tt=bpsi_c(jbj+i)
           tt=apack_c(jbj+i+iboff)+scpr*bpsi_c(jbj+i+iboff)
           apsi_c(jaj+i+iaoff)=apsi_c(jaj+i+iaoff)+tt
        enddo
!        call op_c_inline(length,jaj+iaoff,jbj+iboff,scpr)

        if ((ja1<=jb1 .and. length>=0) .or. iaseg0==maseg_c) exit nonconvex_loop_c !segment is finished
        iaseg0=iaseg0+1
        jb0=max(keybg_c(1,ibseg),keyag_c_lin(iaseg0))
        if (keyag_c_lin(iaseg0)>jb1) exit nonconvex_loop_c !segment is not covered
        jbj=jbj+max(jb0-keybg_c(1,ibseg),0)
     end do nonconvex_loop_c
     !disable loop if the end is reached
     if (iaseg0 == maseg_c .and. keybg_c(1,ibseg)> keyag_c_lin(maseg_c)) iaseg0=iaseg0+1
  enddo
  !stop
  !$omp end do nowait

  ! fine part
  !LG  ibsegs=1
  !LG  ibsege=mbseg_f

  iaseg0=1

  !LG  !separate in chunks the loop among the threads
  !LG  !$ nchunck=max(mbseg_f/nthread,1)
  !LG  !$ ibsegs=min(ithread*nchunck,mbseg_f+1)
  !LG  !$ ibsege=min((ithread+1)*nchunck,mbseg_f)

  !$omp do schedule(static)
  do ibseg=1,mbseg_f
     jbj=keybv_f(ibseg)
     !jb0=keybg_f(1,ibseg)
     jb0=max(keybg_f(1,ibseg),keyag_f_lin(1))
     jb1=keybg_f(2,ibseg)
     iboff = max(jb0-keybg_f(1,ibseg),0)
     call hunt_inline(keyag_f_lin,maseg_f,jb0,iaseg0)
     if (iaseg0==0) then  !segment not belonging to the wavefunctions, go further
        iaseg0=1
        cycle     
     end if
     nonconvex_loop_f: do while(iaseg0 <= maseg_f)
!!$     length = jb1-jb0
!!$     iaoff = jb0-keyag_f_lin(iaseg0)

        ja0=keyag_f_lin(iaseg0) !still doubts about copying in automatic array
        ja1=min(jb1,keyag_f(2,iaseg0)) 
        length = ja1-jb0
        iaoff = max(jb0-ja0,0) !no offset if we are already inside
        jaj=keyav_f(iaseg0)

        do i=0,length
           do j=1,7
              tt=apack_f(j,jbj+i+iboff)+scpr*bpsi_f(j,jbj+i+iboff)
              apsi_f(j,jaj+i+iaoff)=apsi_f(j,jaj+i+iaoff)+tt
           end do
        enddo
        !call op_f_inline(length,jaj+iaoff,jbj+iboff,scpr)

        if ((ja1<=jb1 .and. length>=0) .or. iaseg0==maseg_f) exit nonconvex_loop_f !segment is finished  
        iaseg0=iaseg0+1
        jb0=max(keybg_f(1,ibseg),keyag_f_lin(iaseg0))
        if (keyag_f_lin(iaseg0)>jb1) exit nonconvex_loop_f !segment is not covered 
        jbj=jbj+max(jb0-keybg_f(1,ibseg),0)
     end do nonconvex_loop_f
     !disable loop if the end is reached
     if (iaseg0 == maseg_f .and. keybg_f(1,ibseg)> keyag_f_lin(maseg_f)) iaseg0=iaseg0+1
  enddo
  !$omp end do !implicit barrier 
  !$omp end parallel


contains

!!$  subroutine op_c_inline(length,jaj,jbj,scpr)
!!$    implicit none
!!$    integer, intent(in) :: length,jaj,jbj
!!$    real(wp), intent(in) :: scpr
!!$    !local variables
!!$    integer :: i
!!$    real(wp) :: tt
!!$
!!$     do i=0,length
!!$        !tt=bpsi_c(jbj+i)
!!$        tt=apack_c(jbj+i)+scpr*bpsi_c(jbj+i)
!!$        apsi_c(jaj+i)=apsi_c(jaj+i)+tt
!!$     enddo
!!$
!!$  end subroutine op_c_inline
!!$
!!$  subroutine op_f_inline(length,jaj,jbj,scpr)
!!$    implicit none
!!$    integer, intent(in) :: length,jaj,jbj
!!$    real(wp), intent(in) :: scpr
!!$    !local variables
!!$    integer :: i,j
!!$    real(wp) :: tt
!!$
!!$     do i=0,length
!!$        do j=1,7
!!$           tt=apack_f(j,jbj+i)+scpr*bpsi_f(j,jbj+i)
!!$           apsi_f(j,jaj+i)=apsi_f(j,jaj+i)+tt
!!$        end do
!!$     enddo
!!$  end subroutine op_f_inline

  include 'scalar_product-inc.f90'
  
END SUBROUTINE waxpy_keys_unpack


!> Calculates the dot product between a wavefunctions apsi and a projector bpsi (both in compressed form)
!! The array mask is used so that application of the projector can
!! be done without bitonic search
!! The array of the wavefunction is also compressed
!! so that successive projector application can be performed also with linear algebra routines
subroutine wpdot_mask_pack(  &
     mavctr_c,mavctr_f,mseg_c,mseg_f,amask_c,amask_f,apsi_c,apsi_f,  &
     mbvctr_c,mbvctr_f,bpsi_c,bpsi_f,&
     apack_c,apack_f,scpr)
  use module_base, only: wp
  implicit none
  integer, intent(in) :: mavctr_c,mavctr_f,mseg_c,mseg_f,mbvctr_c,mbvctr_f
  integer, dimension(3,mseg_c), intent(in) :: amask_c
  integer, dimension(3,mseg_f), intent(in) :: amask_f
  real(wp), dimension(mavctr_c), intent(in) :: apsi_c
  real(wp), dimension(mbvctr_c), intent(in) :: bpsi_c
  real(wp), dimension(7,mavctr_f), intent(in) :: apsi_f
  real(wp), dimension(7,mbvctr_f), intent(in) :: bpsi_f
  real(wp), dimension(mbvctr_c), intent(inout) :: apack_c
  real(wp), dimension(7,mbvctr_f), intent(inout) :: apack_f
  real(wp), intent(out) :: scpr
  !local variables
  integer :: iaseg,jaj,jbj,length,i,j
  real(wp) :: scpr1,scpr0,tt

  scpr=0.0_wp

  !$omp parallel default(shared) &
  !$omp private(i,jaj,length,tt,j,scpr0,scpr1) &
  !$omp private(jbj,iaseg)

  scpr0=0.0_wp
  scpr1=0.0_wp

  !$omp do !schedule(static)
  do iaseg=1,mseg_c
     !with the masking array there is no need to perform the bitonic search anymore
     length=amask_c(1,iaseg) !number of elements to be copied
     jaj   =amask_c(2,iaseg) !starting point in original array
     jbj   =amask_c(3,iaseg) !starting point in packed array
     do i=0,length-1 !reduced by one
        tt=apsi_c(jaj+i)
        scpr0=scpr0+tt*bpsi_c(jbj+i)
        apack_c(jbj+i)=tt
     enddo
  end do
  !$omp end do nowait

  !$omp do !schedule(static)
  do iaseg=1,mseg_f
     !with the masking array there is no need to perform the bitonic search anymore
     length=amask_f(1,iaseg) !number of elements to be copied
     jaj   =amask_f(2,iaseg) !starting point in original array
     jbj   =amask_f(3,iaseg) !starting point in packed array
     do i=0,length-1 !reduced by one
        do j=1,7
           tt=apsi_f(j,jaj+i)
           scpr1=scpr1+tt*bpsi_f(j,jbj+i)
           apack_f(j,jbj+i)=tt
        end do
     enddo
  end do
  !$omp end do !implicit barrier 

  scpr0=scpr0+scpr1

  !$omp critical 
  scpr=scpr+scpr0
  !$omp end critical

  !$omp end parallel
  
END SUBROUTINE wpdot_mask_pack

!> Calculates the dot product between a wavefunctions apsi and a projector bpsi (both in compressed form)
!! The array mask is used so that application of the projector can
!! be done without bitonic search
subroutine wpdot_mask(  &
     mavctr_c,mavctr_f,mseg_c,mseg_f,amask_c,amask_f,apsi_c,apsi_f,  &
     mbvctr_c,mbvctr_f,bpsi_c,bpsi_f,&
     scpr)
  use module_base, only: wp
  implicit none
  integer, intent(in) :: mavctr_c,mavctr_f,mseg_c,mseg_f,mbvctr_c,mbvctr_f
  integer, dimension(3,mseg_c), intent(in) :: amask_c
  integer, dimension(3,mseg_f), intent(in) :: amask_f
  real(wp), dimension(mavctr_c), intent(in) :: apsi_c
  real(wp), dimension(mbvctr_c), intent(in) :: bpsi_c
  real(wp), dimension(7,mavctr_f), intent(in) :: apsi_f
  real(wp), dimension(7,mbvctr_f), intent(in) :: bpsi_f
  real(wp), intent(out) :: scpr
  !local variables
  integer :: iaseg,jaj,jbj,length,i,j
  real(wp) :: scpr1,scpr0,tt

  scpr=0.0_wp

  !$omp parallel default(shared) &
  !$omp private(i,jaj,length,tt,j,scpr0,scpr1) &
  !$omp private(jbj,iaseg)

  scpr0=0.0_wp
  scpr1=0.0_wp

  !$omp do !schedule(static)
  do iaseg=1,mseg_c
     !with the masking array there is no need to perform the bitonic search anymore
     length=amask_c(1,iaseg) !number of elements to be copied
     jaj   =amask_c(2,iaseg) !starting point in original array
     jbj   =amask_c(3,iaseg) !starting point in packed array
     do i=0,length-1 !reduced by one
        tt=apsi_c(jaj+i)
        scpr0=scpr0+tt*bpsi_c(jbj+i)
     enddo
  end do
  !$omp end do nowait

  !$omp do !schedule(static)
  do iaseg=1,mseg_f
     !with the masking array there is no need to perform the bitonic search anymore
     length=amask_f(1,iaseg) !number of elements to be copied
     jaj   =amask_f(2,iaseg) !starting point in original array
     jbj   =amask_f(3,iaseg) !starting point in packed array
     do i=0,length-1 !reduced by one
        do j=1,7
           tt=apsi_f(j,jaj+i)
           scpr1=scpr1+tt*bpsi_f(j,jbj+i)
        end do
     enddo
  end do
  !$omp end do !implicit barrier 

  scpr0=scpr0+scpr1

  !$omp critical 
  scpr=scpr+scpr0
  !$omp end critical

  !$omp end parallel
  
END SUBROUTINE wpdot_mask


!> Rank 1 update of wavefunction a with wavefunction b: apsi=apsi+scpr*bpsi
!! The update is only done in the localization region of apsi
!! The array mask is used so that application of the projector can
!! be done without bitonic search
subroutine waxpy_mask(  &
     mavctr_c,mavctr_f,mseg_c,mseg_f,amask_c,amask_f,apsi_c,apsi_f,  &
     mbvctr_c,mbvctr_f,bpsi_c,bpsi_f,scpr)
  use module_base, only: wp
  implicit none
  integer, intent(in) :: mavctr_c,mavctr_f,mseg_c,mseg_f,mbvctr_c,mbvctr_f
  real(wp), intent(in) :: scpr
  integer, dimension(3,mseg_c), intent(in) :: amask_c
  integer, dimension(3,mseg_f), intent(in) :: amask_f
  real(wp), dimension(mbvctr_c), intent(in) :: bpsi_c
  real(wp), dimension(7,mbvctr_f), intent(in) :: bpsi_f
  real(wp), dimension(mavctr_c), intent(inout) :: apsi_c
  real(wp), dimension(7,mavctr_f), intent(inout) :: apsi_f

  !local variables
  integer :: iaseg,jaj,jbj,length,i,j
  real(wp) :: tt

  !quick return if possible
  if (scpr==0.0_wp) return

  !$omp parallel default(shared) &
  !$omp private(i,jaj,length,tt,j) &
  !$omp private(jbj,iaseg)

  !$omp do !schedule(static)
  do iaseg=1,mseg_c
     !with the masking array there is no need to perform the bitonic search anymore
     length=amask_c(1,iaseg) !number of elements to be copied
     jaj   =amask_c(2,iaseg) !starting point in original array
     jbj   =amask_c(3,iaseg) !starting point in packed array
     do i=0,length-1 !reduced by one
        tt=bpsi_c(jbj+i)
        apsi_c(jaj+i)=apsi_c(jaj+i)+scpr*tt
     enddo
  end do
  !$omp end do nowait

  !$omp do !schedule(static)
  do iaseg=1,mseg_f
     !with the masking array there is no need to perform the bitonic search anymore
     length=amask_f(1,iaseg) !number of elements to be copied
     jaj   =amask_f(2,iaseg) !starting point in original array
     jbj   =amask_f(3,iaseg) !starting point in packed array
     do i=0,length-1 !reduced by one
        do j=1,7
           tt=bpsi_f(j,jbj+i)
           apsi_f(j,jaj+i)=apsi_f(j,jaj+i)+scpr*tt
        end do
     enddo
  end do
  !$omp end do !implicit barrier 

  !$omp end parallel
  
END SUBROUTINE waxpy_mask

!> Rank 1 update of wavefunction a with wavefunction b: apsi=apsi+scpr*bpsi
!! The update is only done in the localization region of apsi
!! The array mask is used so that application of the projector can
!! be done without bitonic search
subroutine waxpy_mask_unpack(  &
     mavctr_c,mavctr_f,mseg_c,mseg_f,amask_c,amask_f,apack_c,apack_f,&
     apsi_c,apsi_f, &
     mbvctr_c,mbvctr_f,bpsi_c,bpsi_f,scpr)
  use module_base, only: wp
  implicit none
  integer, intent(in) :: mavctr_c,mavctr_f,mseg_c,mseg_f,mbvctr_c,mbvctr_f
  real(wp), intent(in) :: scpr
  integer, dimension(3,mseg_c), intent(in) :: amask_c
  integer, dimension(3,mseg_f), intent(in) :: amask_f
  real(wp), dimension(mbvctr_c), intent(in) :: bpsi_c
  real(wp), dimension(7,mbvctr_f), intent(in) :: bpsi_f
  real(wp), dimension(mbvctr_c), intent(in) :: apack_c
  real(wp), dimension(7,mbvctr_f), intent(in) :: apack_f
  real(wp), dimension(mavctr_c), intent(inout) :: apsi_c
  real(wp), dimension(7,mavctr_f), intent(inout) :: apsi_f
  !local variables
  integer :: iaseg,jaj,jbj,length,i,j
  real(wp) :: tt

  !$omp parallel default(shared) &
  !$omp private(i,jaj,length,tt,j) &
  !$omp private(jbj,iaseg)

  !$omp do !schedule(static)
  do iaseg=1,mseg_c
     !with the masking array there is no need to perform the bitonic search anymore
     length=amask_c(1,iaseg) !number of elements to be copied
     jaj   =amask_c(2,iaseg) !starting point in original array
     jbj   =amask_c(3,iaseg) !starting point in packed array
     do i=0,length-1 !reduced by one
        !tt=bpsi_c(jbj+i)
        tt=apack_c(jbj+i)+scpr*bpsi_c(jbj+i)
        apsi_c(jaj+i)=apsi_c(jaj+i)+tt
     enddo
  end do
  !$omp end do nowait

  !$omp do !schedule(static)
  do iaseg=1,mseg_f
     !with the masking array there is no need to perform the bitonic search anymore
     length=amask_f(1,iaseg) !number of elements to be copied
     jaj   =amask_f(2,iaseg) !starting point in original array
     jbj   =amask_f(3,iaseg) !starting point in packed array
     do i=0,length-1 !reduced by one
        do j=1,7
           !tt=bpsi_f(j,jbj+i)
           tt=apack_f(j,jbj+i)+scpr*bpsi_f(j,jbj+i)
           apsi_f(j,jaj+i)=apsi_f(j,jaj+i)+tt
        end do
     enddo
  end do
  !$omp end do !implicit barrier 

  !$omp end parallel
  
END SUBROUTINE waxpy_mask_unpack

!> find the number of chunks which are needed to perform blas operations among two compressed wavefunctions
subroutine count_wblas_segs(maseg,mbseg,keyag_lin,keyag,keybg,nbsegs)
  implicit none
  integer, intent(in) :: maseg,mbseg
  integer, dimension(maseg), intent(in) :: keyag_lin !>linear version of second indices of keyag
  integer, dimension(2,maseg), intent(in) :: keyag !>values of the keys ordered for compression a
  integer, dimension(2,mbseg), intent(in) :: keybg !>values of the keys ordered for compression b
  integer, dimension(mbseg), intent(inout) :: nbsegs !>number of common segments for each segment of b
  !local variables
  integer :: ibseg,jb1,jb0,length,ja0,ja1,imask
  integer :: iaseg0,ibsegs,ibsege

  !$omp parallel default (none) &
  !$omp shared (maseg,keyag,keyag_lin,keybg,mbseg,nbsegs)&
!!$  !$omp parallel default(shared) &
  !$omp private(length,ja1,ja0,jb1,jb0,imask) &
  !$omp private(ibseg,iaseg0)!!!,ithread,nthread,ibsegs,ibsege,nchunk)

  !alternative way of parallelizing the loop, to be tested to explore performances
  !LG  ibsegs=1
  !LG  ibsege=mbseg
  !LG  !$ ithread=omp_get_thread_num()
  !LG  !$ nthread=omp_get_num_threads() 

  iaseg0=1 

  !coarse part. Loop on the projectors segments
  !LG  !separate in chunks the loop among the threads
  !LG  !$ nchunck=max(mbseg/nthread,1)
  !LG  !$ ibsegs=min(ithread*nchunck,mbseg+1)
  !LG  !$ ibsege=min((ithread+1)*nchunck,mbseg)

  !$omp do schedule(static)
  do ibseg=1,mbseg
     jb0=max(keybg(1,ibseg),keyag_lin(1)) !starting point of projector segment
     jb1=keybg(2,ibseg) !ending point of projector segment

     !first passage
     imask=0

     !find the starting point of the wavefunction segment
     !warning: hunt is assuming that the variable is always found
     !if it is not, iaseg0 is put to maseg + 1 so that the loop is disabled
     call hunt_inline(keyag_lin,maseg,jb0,iaseg0)
     if (iaseg0==0) then  !segment not belonging to the wavefunctions, go further
        iaseg0=1
        cycle     
     end if
     !now pass through all the wavefunction segments until the end of the segment is 
     !still contained in projector segment
     nonconvex_loop: do while(iaseg0 <= maseg)
        ja0=keyag_lin(iaseg0)
        ja1=min(jb1,keyag(2,iaseg0)) 
        length = ja1-jb0

        !count the active segments for this ibseg
        if (length+1 > 0) then
           imask=imask+1
        end if
        
        if ((ja1<=jb1 .and. length>=0) .or. iaseg0==maseg) exit nonconvex_loop !segment is finished
        iaseg0=iaseg0+1
        jb0=max(keybg(1,ibseg),keyag_lin(iaseg0))
        if (keyag_lin(iaseg0)>jb1) exit nonconvex_loop !segment is not covered
     end do nonconvex_loop

     !first passage, fill the segments array
     nbsegs(ibseg)=imask

     !disable loop if the end is reached
     if (iaseg0 == maseg .and. keybg(1,ibseg)> keyag_lin(maseg)) iaseg0=iaseg0+1
  enddo
  !$omp end do 
  !$omp end parallel

  contains

    include 'scalar_product-inc.f90'

end subroutine count_wblas_segs

!> count the total number of segments and define the integral array of displacements
pure subroutine integrate_nseg(mseg,msegs,nseg_tot)
  implicit none
  integer, intent(in) :: mseg
  integer, dimension(mseg), intent(inout) :: msegs
  integer, intent(out) :: nseg_tot
  !local variables
  integer :: iseg,jseg

  nseg_tot=0
  do iseg=1,mseg
     jseg=msegs(iseg)
     msegs(iseg)=nseg_tot
     nseg_tot=nseg_tot+jseg
  end do
end subroutine integrate_nseg

!> find the number of chunks which are needed to perform blas operations among two compressed wavefunctions
subroutine fill_wblas_segs(maseg,mbseg,mask_segs,isegs_offset,keyag_lin,keyag,keybg,keyav,keybv,amask)
  implicit none
  integer, intent(in) :: maseg,mbseg,mask_segs
  integer, dimension(maseg), intent(in) :: keyav !>position of the segments in compressed a storage
  integer, dimension(mbseg), intent(in) :: keybv !>position of the segments in compressed b storage
  integer, dimension(maseg), intent(in) :: keyag_lin !>linear version of second indices of keyag
  integer, dimension(mbseg), intent(in) :: isegs_offset !>displacement in common segments for each segment of b. 
                                                            !! Integral function of nbsegs vector
  integer, dimension(2,maseg), intent(in) :: keyag !>values of the keys ordered for compression a
  integer, dimension(2,mbseg), intent(in) :: keybg !>values of the keys ordered for compression b
  integer, dimension(3,mask_segs), intent(out) :: amask !>masking array and positions in compressed storage
  !local variables
  integer :: ibseg,jaj,jb1,jb0,jbj,iaoff,iboff,length,ja0,ja1,imask
  integer :: iaseg0,ibsegs,ibsege

  !$omp parallel default (none) &
  !$omp shared (maseg,keyav,keyag,keyag_lin,keybg,mbseg,keybv,isegs_offset,amask)&
!!$  !$omp parallel default(shared) &
  !$omp private(jaj,iaoff,length,ja1,ja0,jb1,jb0,iboff,imask) &
  !$omp private(jbj,ibseg,iaseg0)!!!,ithread,nthread,ibsegs,ibsege,nchunk)

  !alternative way of parallelizing the loop, to be tested to explore performances
  !LG  ibsegs=1
  !LG  ibsege=mbseg
  !LG  !$ ithread=omp_get_thread_num()
  !LG  !$ nthread=omp_get_num_threads() 

  iaseg0=1 

  !coarse part. Loop on the projectors segments
  !LG  !separate in chunks the loop among the threads
  !LG  !$ nchunck=max(mbseg/nthread,1)
  !LG  !$ ibsegs=min(ithread*nchunck,mbseg+1)
  !LG  !$ ibsege=min((ithread+1)*nchunck,mbseg)

  !$omp do schedule(static)
  do ibseg=1,mbseg
     jbj=keybv(ibseg)
     jb0=max(keybg(1,ibseg),keyag_lin(1)) !starting point of projector segment
     jb1=keybg(2,ibseg) !ending point of projector segment
     iboff = max(jb0-keybg(1,ibseg),0)

     !second passage, retrieve starting point
     imask=isegs_offset(ibseg)

     !find the starting point of the wavefunction segment
     !warning: hunt is assuming that the variable is always found
     !if it is not, iaseg0 is put to maseg + 1 so that the loop is disabled
     call hunt_inline(keyag_lin,maseg,jb0,iaseg0)
     if (iaseg0==0) then  !segment not belonging to the wavefunctions, go further
        iaseg0=1
        cycle     
     end if
     !now pass through all the wavefunction segments until the end of the segment is 
     !still contained in projector segment
     nonconvex_loop: do while(iaseg0 <= maseg)
        !length = jb1-jb0
        !iaoff = jb0-keyag_lin(iaseg0)!jb0-ja0

        ja0=keyag_lin(iaseg0)
        ja1=min(jb1,keyag(2,iaseg0)) 
        length = ja1-jb0
        iaoff = max(jb0-ja0,0) !no offset if we are already inside

        jaj=keyav(iaseg0)
        
        !second passage: fill the masking elements as they have to be used
        if (length+1 > 0) then
           imask=imask+1
           !with the masking array there should be no need to perform the bitonic search anymore
           amask(1,imask)=length+1  !number of elements to be copied
           amask(2,imask)=jaj+iaoff !starting point in original array
           amask(3,imask)=jbj+iboff !starting point in packed array
        end if

        if ((ja1<=jb1 .and. length>=0) .or. iaseg0==maseg) exit nonconvex_loop !segment is finished
        iaseg0=iaseg0+1
        jb0=max(keybg(1,ibseg),keyag_lin(iaseg0))
        if (keyag_lin(iaseg0)>jb1) exit nonconvex_loop !segment is not covered
        jbj=jbj+max(jb0-keybg(1,ibseg),0)
     end do nonconvex_loop
     !disable loop if the end is reached
     if (iaseg0 == maseg .and. keybg(1,ibseg)> keyag_lin(maseg)) iaseg0=iaseg0+1
  enddo
  !$omp end do 
  !$omp end parallel
 contains

    include 'scalar_product-inc.f90'

  end subroutine fill_wblas_segs
