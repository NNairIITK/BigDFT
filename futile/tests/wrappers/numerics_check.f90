!> @file
!!  Test of some functionalities of the numeric groups
!! @author
!!    Copyright (C) 2016-2017 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS
program numeric_check
  use futile
  use f_harmonics
  character(len=*), parameter :: input1=&
       "  {name: ndim, shortname: n, default: 30,"//&
       "  help_string: Size of the array for multipoles,"//&
       "  help_dict: {Allowed values: integer}}"
  character(len=*), parameter :: input2=&
       "  {name: boldify, shortname: b, default: None,"//&
       "  help_string: Boldify the string as a test,"//&
       "  help_dict: {Allowed values: string scalar}}"
  character(len=*), parameter :: input3=&
       "  {name: blinkify, shortname: l, default: None,"//&
       "  help_string: Make the string blinking,"//&
       "  help_dict: {Allowed values: string scalar}}"

  character(len=*), parameter :: inputs=&
       '-'//input1//f_cr//&
       '-'//input2//f_cr//&
       '-'//input3


  integer :: n,i
  type(f_multipoles) :: mp
  type(dictionary), pointer :: options
  real(f_double), dimension(3) :: rxyz
  real(f_double), dimension(:), allocatable :: density
  
  call f_lib_initialize()
  call yaml_new_document()
  call yaml_argparse(options,inputs)
  n=options//'ndim'

  density=f_malloc(n,id='density')
  call f_random_number(density)

  rxyz=1.0_f_double

  !create random multipoles
  call f_multipoles_create(mp,2)

  do i=1,n
     call f_multipoles_accumulate(mp%Q,mp%lmax,rxyz,density(i))
  end do

  !here we may print the results of the multipole calculations
  call yaml_mapping_open('Multipoles of the array')
  call yaml_map('q0',sum(density))
  call yaml_mapping_close()

  call yaml_mapping_open('Calculated multipoles')
  call yaml_map('q0',mp%Q(0)%ptr)
  call yaml_map('q1',mp%Q(1)%ptr)
  call yaml_map('q2',mp%Q(2)%ptr)
  call yaml_mapping_close()

  call f_multipoles_free(mp)

  call f_free(density)
  call dict_free(options)

  !test of the multipole preserving routine
  !initialize the work arrays needed to integrate with isf
  !names of the routines to be redefined
  call initialize_real_space_conversion(isf_m=mp_isf_order)

  call finalize_real_space_conversion()

  call f_lib_finalize()

end program numeric_check

!> Creates the charge density of a Gaussian function, to be used for the local part
!! of the pseudopotentials (gives the error function term when later processed by the Poisson solver).
subroutine gaussian_density(perx, pery, perz, n1i, n2i, n3i, nbl1, nbl2, nbl3, i3s, n3pi, hxh, hyh, hzh, rx, ry, rz, &
     rloc, zion, multipole_preserving, use_iterator, mp_isf, &
     dpbox, nmpx, nmpy, nmpz, mpx, mpy, mpz, nrho, pot_ion, rholeaked)
  use module_base
  use module_dpbox, only: denspot_distribution, dpbox_iterator, DPB_POT_ION, dpbox_iter, dpbox_iter_next
  use gaussians, only: mp_exp
  implicit none
  ! Calling arguments
  logical,intent(in) :: perx, pery, perz
  integer,intent(in) :: n1i, n2i, n3i, nrho, i3s, n3pi
  real(kind=8),intent(in) :: rloc, rx, ry, rz, hxh, hyh, hzh
  integer,intent(in) :: nbl1, nbl2, nbl3
  integer,intent(in) :: zion !< ionic charge (integer!)
  logical,intent(in) :: multipole_preserving, use_iterator
  integer,intent(in) :: mp_isf !< interpolating scaling function order for the multipole preserving
  integer,intent(in) :: nmpx, nmpy, nmpz !< sizes of the temporary arrays; if too small the code stops
  real(kind=8),dimension(0:nmpx),intent(inout) :: mpx !< temporary array for the exponetials in x direction
  real(kind=8),dimension(0:nmpy),intent(inout) :: mpy !< temporary array for the exponetials in y direction
  real(kind=8),dimension(0:nmpz),intent(inout) :: mpz !< temporary array for the exponetials in z direction
  type(denspot_distribution),intent(in) :: dpbox
  real(dp),dimension(nrho),intent(inout) :: pot_ion
  real(kind=8),intent(inout) :: rholeaked

  ! Local variables
  real(kind=8) :: rlocinv2sq, charge, cutoff, xp, yp, zp
  integer,dimension(2,3) :: nbox
  integer :: i1, i2, i3, isx, iex, isy, iey, isz, iez, j1, j2, j3, ind
  type(dpbox_iterator) :: boxit
  real(gp),parameter :: mp_tiny = 1.e-30_gp
  logical :: gox, goy, goz

  call f_routine(id='gaussian_density')


  !rloc=at%psppar(0,0,atit%ityp)
  rlocinv2sq=0.5_gp/rloc**2
  charge=real(zion,gp)/(2.0_gp*pi*sqrt(2.0_gp*pi)*rloc**3)

  !write(*,*) 'rloc, charge', rloc, charge

  !cutoff of the range
  cutoff=10.0_gp*rloc
  if (multipole_preserving) then
     !We want to have a good accuracy of the last point rloc*10
     !cutoff=cutoff+max(hxh,hyh,hzh)*real(16,kind=gp)
     cutoff=cutoff+max(hxh,hyh,hzh)*real(mp_isf,kind=gp)
  end if
  isx=floor((rx-cutoff)/hxh)
  isy=floor((ry-cutoff)/hyh)
  isz=floor((rz-cutoff)/hzh)

  iex=ceiling((rx+cutoff)/hxh)
  iey=ceiling((ry+cutoff)/hyh)
  iez=ceiling((rz+cutoff)/hzh)

  ! Check whether the temporary arrays are large enough
  if (iex-isx>nmpx) then
     call f_err_throw('Temporary array in x direction too small',err_name='BIGDFT_RUNTIME_ERROR')
  end if
  if (iey-isy>nmpy) then
     call f_err_throw('Temporary array in y direction too small',err_name='BIGDFT_RUNTIME_ERROR')
  end if
  if (iez-isz>nmpz) then
     call f_err_throw('Temporary array in z direction too small',err_name='BIGDFT_RUNTIME_ERROR')
  end if

  do i1=isx,iex
     mpx(i1-isx) = mp_exp(hxh,rx,rlocinv2sq,i1,0,multipole_preserving)
  end do
  do i2=isy,iey
     mpy(i2-isy) = mp_exp(hyh,ry,rlocinv2sq,i2,0,multipole_preserving)
  end do
  do i3=isz,iez
     mpz(i3-isz) = mp_exp(hzh,rz,rlocinv2sq,i3,0,multipole_preserving)
  end do

  if (use_iterator) then
     nbox(1,1)=isx
     nbox(1,2)=isy
     nbox(1,3)=isz
     nbox(2,1)=iex
     nbox(2,2)=iey
     nbox(2,3)=iez
     boxit = dpbox_iter(dpbox,DPB_POT_ION,nbox)
     do while(dpbox_iter_next(boxit))
        xp = mpx(boxit%ibox(1)-nbox(1,1)) * mpy(boxit%ibox(2)-nbox(1,2)) * mpz(boxit%ibox(3)-nbox(1,3))
        pot_ion(boxit%ind) = pot_ion(boxit%ind) - xp*charge
     end do

  else

     do i3=isz,iez
        zp = mpz(i3-isz)
        if (abs(zp) < mp_tiny) cycle
        !call ind_positions(perz,i3,grid%n3,j3,goz) 
        call ind_positions_new(perz,i3,n3i,j3,goz) 
        j3=j3+nbl3+1
        do i2=isy,iey
           yp = zp*mpy(i2-isy)
           if (abs(yp) < mp_tiny) cycle
           !call ind_positions(pery,i2,grid%n2,j2,goy)
           call ind_positions_new(pery,i2,n2i,j2,goy)
           do i1=isx,iex
              xp = yp*mpx(i1-isx)
              if (abs(xp) < mp_tiny) cycle
              !call ind_positions(perx,i1,grid%n1,j1,gox)
              call ind_positions_new(perx,i1,n1i,j1,gox)
              if (j3 >= i3s .and. j3 <= i3s+n3pi-1  .and. goy  .and. gox ) then
                 ind=j1+1+nbl1+(j2+nbl2)*n1i+(j3-i3s)*n1i*n2i
                 pot_ion(ind)=pot_ion(ind)-xp*charge
              else if (.not. goz ) then
                 rholeaked=rholeaked+xp*charge
              endif
           enddo
        enddo
     enddo


  end if

  call f_release_routine()

end subroutine gaussian_density


subroutine 3d_gaussian
     !Determine the maximal bounds for mpx, mpy, mpz (1D-integral)
     if (at%astruct%nat >0) then
        cutoff=10.0_gp*maxval(at%psppar(0,0,:))
     else
        cutoff=0.0
     end if
     if (at%multipole_preserving) then
        !We want to have a good accuracy of the last point rloc*10
        cutoff=cutoff+max(hxh,hyh,hzh)*real(at%mp_isf,kind=gp)
     end if
     !Separable function: do 1-D integrals before and store it.
     mpx = f_malloc( (/ 0 .to. (ceiling(cutoff/hxh) - floor(-cutoff/hxh)) + 1 /),id='mpx')
     mpy = f_malloc( (/ 0 .to. (ceiling(cutoff/hyh) - floor(-cutoff/hyh)) + 1 /),id='mpy')
     mpz = f_malloc( (/ 0 .to. (ceiling(cutoff/hzh) - floor(-cutoff/hzh)) + 1 /),id='mpz')
     atit = atoms_iter(at%astruct)
     do while(atoms_iter_next(atit))

        !!-     do iat=1,at%astruct%nat
        !!-        ityp=at%astruct%iatype(iat)
        !!-       rx=rxyz(1,iat) 
        !!-       ry=rxyz(2,iat)
        !!-       rz=rxyz(3,iat)
        rx=rxyz(1,atit%iat) 
        ry=rxyz(2,atit%iat)
        rz=rxyz(3,atit%iat)

        !!-        rloc=at%psppar(0,0,ityp)
        !!-        charge=real(at%nelpsp(ityp),gp)/(2.0_gp*pi*sqrt(2.0_gp*pi)*rloc**3)
        rloc=at%psppar(0,0,atit%ityp)
        rlocinv2sq=0.5_gp/rloc**2
        charge=real(at%nelpsp(atit%ityp),gp)/(2.0_gp*pi*sqrt(2.0_gp*pi)*rloc**3)

        !cutoff of the range
        cutoff=10.0_gp*rloc
        if (at%multipole_preserving) then
           !We want to have a good accuracy of the last point rloc*10
           cutoff=cutoff+max(hxh,hyh,hzh)*real(at%mp_isf,kind=gp)
        end if

        if (use_iterator) then
           nbox(1,1)=floor((rx-cutoff)/hxh)
           nbox(1,2)=floor((ry-cutoff)/hyh)
           nbox(1,3)=floor((rz-cutoff)/hzh)
           nbox(2,1)=ceiling((rx+cutoff)/hxh)
           nbox(2,2)=ceiling((ry+cutoff)/hyh)
           nbox(2,3)=ceiling((rz+cutoff)/hzh)

           !Separable function: do 1-D integrals before and store it.
           !mpx = f_malloc( (/ nbox(1,1).to.nbox(2,1) /),id='mpx')
           !mpy = f_malloc( (/ nbox(1,2).to.nbox(2,2) /),id='mpy')
           !mpz = f_malloc( (/ nbox(1,3).to.nbox(2,3) /),id='mpz')
           do i1=nbox(1,1),nbox(2,1)
              mpx(i1-nbox(1,1)) = mp_exp(hxh,rx,rlocinv2sq,i1,0,at%multipole_preserving)
           end do
           do i2=nbox(1,2),nbox(2,2)
              mpy(i2-nbox(1,2)) = mp_exp(hyh,ry,rlocinv2sq,i2,0,at%multipole_preserving)
           end do
           do i3=nbox(1,3),nbox(2,3)
              mpz(i3-nbox(1,3)) = mp_exp(hzh,rz,rlocinv2sq,i3,0,at%multipole_preserving)
           end do

        else
           isx=floor((rx-cutoff)/hxh)
           isy=floor((ry-cutoff)/hyh)
           isz=floor((rz-cutoff)/hzh)

           iex=ceiling((rx+cutoff)/hxh)
           iey=ceiling((ry+cutoff)/hyh)
           iez=ceiling((rz+cutoff)/hzh)

           !Separable function: do 1-D integrals before and store it.
           !call mp_calculate(rx,ry,rz,hxh,hyh,hzh,cutoff,rlocinv2sq,at%multipole_preserving,mpx,mpy,mpz)
           !!mpx = f_malloc( (/ isx.to.iex /),id='mpx')
           !!mpy = f_malloc( (/ isy.to.iey /),id='mpy')
           !!mpz = f_malloc( (/ isz.to.iez /),id='mpz')
           do i1=isx,iex
              mpx(i1-isx) = mp_exp(hxh,rx,rlocinv2sq,i1,0,at%multipole_preserving)
           end do
           do i2=isy,iey
              mpy(i2-isy) = mp_exp(hyh,ry,rlocinv2sq,i2,0,at%multipole_preserving)
           end do
           do i3=isz,iez
              mpz(i3-isz) = mp_exp(hzh,rz,rlocinv2sq,i3,0,at%multipole_preserving)
           end do

        end if

        if ( .not. any(at%npspcode == PSPCODE_PAW) ) then

           if (use_iterator) then
              boxit = dpbox_iter(dpbox,DPB_POT_ION,nbox)
              do while(dpbox_iter_next(boxit))
                 xp = mpx(boxit%ibox(1)-nbox(1,1)) * mpy(boxit%ibox(2)-nbox(1,2)) * mpz(boxit%ibox(3)-nbox(1,3))
                 pot_ion(boxit%ind) = pot_ion(boxit%ind) - xp*charge
                 !write(*,'(4(i0,1x),2(1pe20.10))') boxit%ibox(1),boxit%ibox(2),boxit%ibox(3),boxit%ind,xp,pot_ion(boxit%ind)
              end do
           else


              !Calculate Ionic Density using HGH parameters.
              !Eq. 1.104, T. Deutsch and L. Genovese, JDN. 12, 2011
              do i3=isz,iez
                 zp = mpz(i3-isz)
                 !call ind_positions(perz,i3,n3,j3,goz)
                 call ind_positions_new(perz,i3,n3i,j3,goz) 
                 j3=j3+nbl3+1
                 if ( goz .and. (j3<i3s.or.j3>i3s+n3pi-1) ) cycle
                 indj3=(j3-i3s)*n1i*n2i
                 do i2=isy,iey
                    yp = zp*mpy(i2-isy)
                    !call ind_positions(pery,i2,n2,j2,goy)
                    call ind_positions_new(pery,i2,n2i,j2,goy)
                    if (goz.and.(.not.goy)) cycle
                    indj23=1+nbl1+(j2+nbl2)*n1i+indj3
                    do i1=isx,iex
                       xp = yp*mpx(i1-isx)
                       !call ind_positions(perx,i1,n1,j1,gox)
                       call ind_positions_new(perx,i1,n1i,j1,gox)
                       if (j3 >= i3s .and. j3 <= i3s+n3pi-1 .and. goy .and. gox) then
                          ind=j1+indj23
                          pot_ion(ind)=pot_ion(ind)-xp*charge
                          !write(*,'(4(i0,1x),2(1pe20.10))') i1,i2,i3,ind,xp,pot_ion(ind)
                       else if (.not. goz ) then
                          rholeaked=rholeaked+xp*charge
                       endif
                    enddo
                 enddo
              enddo
           end if

        else
           !r2paw=at%pawtab(ityp)%rpaw**2
           if (use_iterator) then
              !r2paw=at%pawtab(atit%ityp)%rpaw**2
              do while(dpbox_iter_next(boxit))
                 xp = mpx(boxit%ibox(1)-nbox(1,1)) * mpy(boxit%ibox(2)-nbox(1,2)) * mpz(boxit%ibox(3)-nbox(1,3))
                 r2 = (boxit%x-rx)**2 + (boxit%y-ry)**2 + (boxit%z-rz)**2                 
                 if (.not. pawErfCorrection) then
                    rr1(1) = sqrt(r2)
                    !This converges very slowly
                    call paw_splint(at%pawtab(atit%ityp)%wvl%rholoc%msz, &
                         & at%pawtab(atit%ityp)%wvl%rholoc%rad, &
                         & at%pawtab(atit%ityp)%wvl%rholoc%d(:,1), &
                         & at%pawtab(atit%ityp)%wvl%rholoc%d(:,2), &
                         & 1,rr1,raux1,ierr)
                 else
                    !Take the HGH form for rho_L (long range)
                    raux1(1)=-xp*charge
                 end if
                 !raux=-4.d0**(3.0d0/2.0d0)*exp(-4.d0*pi*r2)
                 !Rholeaked is not calculated!!
                 pot_ion(boxit%ind) = pot_ion(boxit%ind) + raux1(1)
              enddo
           else
              !Calculate Ionic Density using splines, PAW case
              !r2paw=at%pawtab(atit%ityp)%rpaw**2
              do i3=isz,iez
                 zp = mpz(i3-isz)
                 if (abs(zp) < mp_tiny) cycle
                 !call ind_positions(perz,i3,n3,j3,goz)
                 call ind_positions_new(perz,i3,n3i,j3,goz) 
                 j3=j3+nbl3+1
                 indj3=(j3-i3s)*n1i*n2i
                 z=real(i3,gp)*hzh-rz
                 zsq=z**2
                 do i2=isy,iey
                    yp = zp*mpy(i2-isy)
                    if (abs(yp) < mp_tiny) cycle
                    !call ind_positions(pery,i2,n2,j2,goy)
                    call ind_positions_new(pery,i2,n2i,j2,goy)
                    indj23=1+nbl1+(j2+nbl2)*n1i+indj3
                    y=real(i2,gp)*hyh-ry
                    yzsq=y**2+zsq
                    do i1=isx,iex
                       xp = yp*mpx(i1-isx)
                       if (abs(xp) < mp_tiny) cycle
                       !call ind_positions(perx,i1,n1,j1,gox)
                       call ind_positions_new(perx,i1,n1i,j1,gox)
                       x=real(i1,gp)*hxh-rx
                       r2=x**2+yzsq
                       !if(r2>r2paw) cycle
                       if(.not. pawErfCorrection) then
                          !This converges very slowly                
                          rr1(1)=sqrt(r2)
                          call paw_splint(at%pawtab(atit%ityp)%wvl%rholoc%msz, &
                               & at%pawtab(atit%ityp)%wvl%rholoc%rad, &
                               & at%pawtab(atit%ityp)%wvl%rholoc%d(:,1), &
                               & at%pawtab(atit%ityp)%wvl%rholoc%d(:,2), &
                               & 1,rr1,raux1,ierr)
                       else
                          !Take the HGH form for rho_L (long range)
                          raux1(1)=-xp*charge
                       end if
                       !raux=-4.d0**(3.0d0/2.0d0)*exp(-4.d0*pi*r2)

                       if (j3 >= i3s .and. j3 <= i3s+n3pi-1  .and. goy  .and. gox ) then
                          ind=j1+indj23
                          pot_ion(ind)=pot_ion(ind)+raux1(1)
                       else if (.not. goz) then
                          rholeaked=rholeaked-raux1(1)
                       endif
                    enddo
                 enddo
              enddo
           end if

        end if

        !De-allocate for multipole preserving
        !call f_free(mpx,mpy,mpz)

     end do
   end subroutine 3d_gaussian
