!> @file
!! module implementing the fingerprints
!!
!! @author NAN
!! @section LICENCE
!!    Copyright (C) 2014 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS

module module_fingerprints

    private

    public :: init_fingerprint
    public :: finalize_fingerprint
    public :: fingerprint
    public :: fpdistance
    public :: equal


    contains

!=====================================================================
subroutine finalize_fingerprint(fp)
    use SPREDbase
    implicit none
    !parameters
    real(gp), allocatable, intent(inout) :: fp(:)

    call f_free(fp) 
end subroutine finalize_fingerprint
!=====================================================================
subroutine init_fingerprint(inputs,nat,geocode,nid,fp)
    use SPREDbase
    use SPREDtypes
   !! use f_enums, f_str => str
    implicit none
    !parameters
    type(SPRED_inputs), intent(in) :: inputs
    integer, intent(in) :: nat
    character(len=*), intent(in) :: geocode
    integer, intent(out) :: nid   
    real(gp), allocatable, intent(out) :: fp(:) !the fingerprint array
    !internal

    if(allocated(fp))then
        call f_err_throw('Array fp is already allocated.')
    endif

    select case(trim(f_char(inputs%fp_method)))
      case('OMF_FP_METHOD')
        if(geocode/='F')then
          call f_err_throw('geocode /= F but a fingerprint (OMF) for free BC is used')
!!          call yaml_warning('geocode /= F but a fingerprint (OMF) for free BC is used')
        endif
        nid=inputs%fp_angmom*nat
        fp = f_malloc((/ 1.to.nid/),id='fp')
      case('OMP_FP_METHOD')
        if(geocode/='P' .and. geocode/='S')then
          call f_err_throw('geocode /= P or S but a fingerprint (OMP) for periodic BC is used')
        endif
        nid=inputs%fp_angmom*inputs%fp_natx_sphere*nat
        fp = f_malloc((/ 1.to.nid/),id='fp')
      case('OMPOLD_FP_METHOD')
        if(geocode/='P')then
          call f_err_throw('geocode /= P but a fingerprint (OMPOLD) for periodic BC is used')
        endif
        nid=inputs%fp_angmom*nat
        fp = f_malloc((/ 1.to.nid/),id='fp')
      case('OMSOLD_FP_METHOD')
        if(geocode/='S')then
          call f_err_throw('geocode /= S but a fingerprint (OMSOLD) for slab BC is used')
        endif
        nid=inputs%fp_angmom*nat
        fp = f_malloc((/ 1.to.nid/),id='fp')
      case DEFAULT
        call f_err_throw('Following FP method is unknown: '//trim(f_char(inputs%fp_method)))
    end select
end subroutine init_fingerprint
!=================================================================
subroutine fingerprint(inputs,nidIn,nat,alat,rcov,rxyz,fp)
    use SPREDbase
    use SPREDtypes
    implicit none
    !parameters
    type(SPRED_inputs), intent(in) :: inputs
    integer, intent(in) :: nidIn
    integer, intent(in) :: nat
    real(gp), intent(in)         :: alat(3)
    real(gp), intent(in) :: rcov(nat)
    real(gp), intent(in) :: rxyz(3,nat)
    real(gp), intent(out) :: fp(nidIn)
    !internal
    integer :: nid
    real(gp) :: alat_dmy(3,3)
 


    select case(trim(f_char(inputs%fp_method)))
      case('OMF_FP_METHOD')
        nid=inputs%fp_angmom*nat
        if(nidIn/=nid) then
          call f_err_throw('Array fp has wrong size')
        endif
        call fingerprint_freebc(nat,nid,alat,'F',rcov,rxyz,fp)
      case('OMP_FP_METHOD')
        nid=inputs%fp_angmom*inputs%fp_natx_sphere*nat
        if(nidIn/=nid) then
          call f_err_throw('Array fp has wrong size')
        endif

        !alat_dmy can be removed as soon as bigdft is able to handle
        !non-orthorombic cells (more precisely, as soon as the lattice vector
        !data structure in bigdft is of dimension (3,3))
        alat_dmy(1,1)=alat(1)
        alat_dmy(2,1)=0.0_gp
        alat_dmy(3,1)=0.0_gp
        alat_dmy(1,2)=0.0_gp
        alat_dmy(2,2)=alat(2)
        alat_dmy(3,2)=0.0_gp
        alat_dmy(1,3)=0.0_gp
        alat_dmy(2,3)=0.0_gp
        alat_dmy(3,3)=alat(3)

        call fingerprint_periodic(nat, inputs%fp_natx_sphere, inputs%fp_angmom, alat_dmy, rxyz, rcov, fp)
      case('OMPOLD_FP_METHOD')
        nid=inputs%fp_angmom*nat
        if(nidIn/=nid) then
          call f_err_throw('Array fp has wrong size')
        endif
        call fingerprint_freebc(nat,nid,alat,'P',rcov,rxyz,fp)
      case('OMSOLD_FP_METHOD')
        nid=inputs%fp_angmom*nat
        if(nidIn/=nid) then
          call f_err_throw('Array fp has wrong size')
        endif
        call fingerprint_freebc(nat,nid,alat,'S',rcov,rxyz,fp)
      case DEFAULT
        call f_err_throw('Following FP method is unknown: '//trim(f_char(inputs%fp_method)))
    end select
end subroutine fingerprint
!=====================================================================
subroutine fingerprint_periodic(nat, natx_sphere, lseg, alat, rxyz, rcov, fpall)
  implicit real*8 (a-h,o-z)
  parameter(nwork=100)
  dimension workalat(nwork) 
  dimension rxyz_sphere(3, natx_sphere),rcov_sphere(natx_sphere)
  dimension fpall(lseg*natx_sphere,nat),fp(lseg*natx_sphere),amplitude(natx_sphere)
  dimension rxyz(3,nat),rcov(nat)
  dimension alat(3, 3),alatalat(3,3),eigalat(3)
  allocatable   :: om(:,:) , work(:)
  integer :: nid
  

! parameters for cutoff function
  width_cutoff=6.d0
  nex_cutoff=3
  radius_cutoff=sqrt(2.d0*nex_cutoff)*width_cutoff
  radius_cutoff2=radius_cutoff**2
  factor_cutoff=1.d0/(2.d0*nex_cutoff*width_cutoff**2)
!  write(*,*) 'width_cutoff,radius_cutoff',width_cutoff,radius_cutoff

  do i=1,3
  do j=1,3
  alatalat(i,j)=alat(1,i)*alat(1,j)+alat(2,i)*alat(2,j)+alat(3,i)*alat(3,j)
  enddo
  enddo
  call dsyev('N', 'L', 3, alatalat, 3, eigalat, workalat, nwork, info)
!  write(*,*) 'alat EVals',eigalat
!  write(*,*) 'ixyzmax',int(sqrt(1.d0/eigalat(1))*radius_cutoff)
  ixyzmax= int(sqrt(1.d0/eigalat(1))*radius_cutoff) + 1

! loop over all center atoms
  natsmax=0
  natsmin=1000000
  do iat = 1, nat
        
     nat_sphere=0
     do jat = 1, nat
         do ix = -ixyzmax,ixyzmax
           do iy = -ixyzmax,ixyzmax
             do iz = -ixyzmax,ixyzmax
                xj = rxyz(1, jat) + ix*alat(1,1)+iy*alat(1,2)+iz*alat(1,3)
                yj = rxyz(2, jat) + ix*alat(2,1)+iy*alat(2,2)+iz*alat(2,3)
                zj = rxyz(3, jat) + ix*alat(3,1)+iy*alat(3,2)+iz*alat(3,3)
                dist2 = (xj-rxyz(1, iat))**2+(yj-rxyz(2, iat))**2+(zj-rxyz(3, iat))**2

                if (dist2.le.radius_cutoff2) then
!                    write(*,*) 'dist ',jat,ix,iy,iz,sqrt(dist2)
!                    write(*,*) xj,yj,zj
!                    write(*,*) rxyz(1, jat),rxyz(2, jat),rxyz(3, jat)
                    nat_sphere=nat_sphere+1
                    if (nat_sphere.gt.natx_sphere) stop 'enlarge natx_sphere'
                    amplitude(nat_sphere)=(1.d0 - dist2*factor_cutoff)**nex_cutoff
                    rxyz_sphere(1,nat_sphere)=xj
                    rxyz_sphere(2,nat_sphere)=yj
                    rxyz_sphere(3,nat_sphere)=zj
                    rcov_sphere(nat_sphere)=rcov(jat)
                endif
             enddo
           enddo
         enddo
     enddo
     natsmin=min(natsmin,nat_sphere)
     natsmax=max(natsmax,nat_sphere)

! set up big overlap matrix
          nid=nat_sphere*lseg
          allocate(om(nid,nid))
!          call  create_om(nat_sphere,rxyz_sphere,rcov_sphere,nid,om)
          if (lseg.eq.1) then
          call  create_om_1(nat_sphere,rxyz_sphere,rcov_sphere,om)
          call mltampl_1(nat_sphere,amplitude,om)
          else if (lseg.eq.4) then
          call  OLDcreate_om_4(nat_sphere,rxyz_sphere,rcov_sphere,om)

!          call Wblock(nat_sphere,lseg,om)
!     stop
          call mltampl_4(nat_sphere,amplitude,om)
          else
              stop 'wrong lseg'
          endif

          tt=0.d0
          do i=1,nid
          do j=i,nid
          tt=max(tt,(om(i,j)-om(j,i))**2)
          if ((om(i,j)-om(j,i))**2.gt.1.d-6) write(*,*) i,j,om(i,j),om(j,i)
          enddo
          enddo
          if (tt.gt.1.d-6) write(*,*) 'max dev symmetry',tt

          lwork=max(1,3*nid-1)
          allocate(work(lwork))
          call DSYEV('N','L',nid,om,nid,fp,work,-1,info)
          if (info.ne.0) stop 'info query'
          lwork=nint(work(1))
          deallocate(work)
        
          allocate(work(lwork))
          call DSYEV('N','L',nid,om,nid,fp,work,lwork,info)
          if (info.ne.0) stop 'info OM diagonalisation'
          deallocate(work,om)


          do i=1,nid
          fpall(i,iat)=fp(nid+1-i)
          enddo
          do i=nid+1,lseg*natx_sphere
          fpall(i,iat)=0.d0
          enddo

!          write(*,*) 'fpall,iat=',iat,nid,nat_sphere,lseg
!          write(*,'(10(2x,e10.3))') (fpall(i,iat),i=1,lseg*natx_sphere)
     
          if (fp(1).lt.-1.d-12) then 
              write(*,*) fp(1)
              stop 'error negative EV'
          endif

  end do
           write(*,*) 'min,max number of atoms in sphere ',natsmin,natsmax

end subroutine fingerprint_periodic
         subroutine create_om_1(nat,rxyz,rcov,om)
         implicit real*8 (a-h,o-z)
         dimension rxyz(3,nat),rcov(nat),om(nat,nat)

           ! Gaussian overlap
           !  <sj|si>
           do iat=1,nat
              xi=rxyz(1,iat)
              yi=rxyz(2,iat)
              zi=rxyz(3,iat)

              do jat=1,nat

                 d2=(rxyz(1,jat) -xi)**2 +(rxyz(2,jat)-yi)**2+(rxyz(3,jat)-zi)**2
                 r=.5d0/(rcov(iat)**2 + rcov(jat)**2)
                 om(jat,iat)= sqrt(4.d0*r*(rcov(iat)*rcov(jat)))**3 * exp(-d2*r)

                 !xji=rxyz(1,jat) - xi
                 !yji=rxyz(2,jat) - yi 
                 !zji=rxyz(3,jat) - zi
                 !d2=xji*xji + yji*yji + zji*zji
                 !scov=rcov(iat)**2 + rcov(jat)**2
                 !factor=(sqrt(2.d0*rcov(iat)*rcov(jat)/scov))**3
                 !facp=1.d0/sqrt(factor*scov)
                 !arg=d2*.5d0/scov
                 !fexp=factor*exp(-arg)
           !  <sj|si>
                 !om(jat,iat)= fexp
              enddo
           enddo
  end subroutine create_om_1
         subroutine create_om_4(nat,rxyz,rcov,om)
         implicit real*8 (a-h,o-z)
         dimension rxyz(3,nat),rcov(nat),om(4,nat,4,nat)


           ! Gaussian overlap
           do iat=1,nat
              xi=rxyz(1,iat)
              yi=rxyz(2,iat)
              zi=rxyz(3,iat)

              do jat=1,nat
                 xji=rxyz(1,jat) - xi
                 yji=rxyz(2,jat) - yi
                 zji=rxyz(3,jat) - zi
                 d2=xji*xji + yji*yji + zji*zji
                 scov=rcov(iat)**2 + rcov(jat)**2
                 factor=(sqrt(2.d0*rcov(iat)*rcov(jat)/scov))**3
                 facp=sqrt(scov)/(rcov(iat)*rcov(jat))
                 arg=d2*.5d0/scov
                 fexp=factor*exp(-arg)
           !  <sj|si>
                 om(1,jat,1,iat)= fexp
           !  <pj|si>
                    om(2,jat,1,iat)= xji*(facp*fexp*rcov(jat)**2)/scov
                    om(3,jat,1,iat)= yji*(facp*fexp*rcov(jat)**2)/scov
                    om(4,jat,1,iat)= zji*(facp*fexp*rcov(jat)**2)/scov
                    om(1,jat,2,iat)=-xji*(facp*fexp*rcov(iat)**2)/scov
                    om(1,jat,3,iat)=-yji*(facp*fexp*rcov(iat)**2)/scov
                    om(1,jat,4,iat)=-zji*(facp*fexp*rcov(iat)**2)/scov
            ! <pj|pi> 
                    om(2,jat,2,iat)= fexp*(1.d0 - xji*xji/scov)
                    om(3,jat,2,iat)= (fexp/scov) * yji*xji
                    om(4,jat,2,iat)= (fexp/scov) * zji*xji
                    om(2,jat,3,iat)= (fexp/scov) * xji*yji
                    om(3,jat,3,iat)= fexp*(1.d0 -  yji*yji/scov)
                    om(4,jat,3,iat)= (fexp/scov) * zji*yji
                    om(2,jat,4,iat)= (fexp/scov) * xji*zji
                    om(3,jat,4,iat)= (fexp/scov) * yji*zji
                    om(4,jat,4,iat)= fexp*(1.d0 - zji*zji/scov)
              enddo
           enddo

  end subroutine create_om_4
         subroutine OLDcreate_om_4(nat,rxyz,rcov,om)
         implicit real*8 (a-h,o-z)
         dimension rxyz(3,nat),rcov(nat),om(4,nat,4,nat)


           ! Gaussian overlap
           do iat=1,nat
              xi=rxyz(1,iat)
              yi=rxyz(2,iat)
              zi=rxyz(3,iat)

              do jat=1,nat
                 xji=rxyz(1,jat) - xi
                 yji=rxyz(2,jat) - yi
                 zji=rxyz(3,jat) - zi
                 d2=xji*xji + yji*yji + zji*zji
                 r=.5d0/(rcov(iat)**2 + rcov(jat)**2)
           !  <sj|si>
                 om(1,jat,1,iat)= sqrt(4.d0*r*(rcov(iat)*rcov(jat)))**3 * exp(-d2*r)
           !  <pj|si>
                    sji= sqrt(4.d0*r*(rcov(jat)*rcov(iat)))**3 * exp(-d2*r)
                    tt= sqrt(8.d0) *rcov(jat)*r * sji
                    om(2,jat,1,iat)= tt * xji
                    om(3,jat,1,iat)= tt * yji
                    om(4,jat,1,iat)= tt * zji
                    tt= sqrt(8.d0) *rcov(iat)*r * sji
                    om(1,jat,2,iat)=-tt * xji
                    om(1,jat,3,iat)=-tt * yji
                    om(1,jat,4,iat)=-tt * zji
            ! <pj|pi> 
                    tt = -8.d0*rcov(iat)*rcov(jat) * r*r * sji
                    om(2,jat,2,iat)= tt *(xji* xji - .5d0/r)
                    om(3,jat,2,iat)= tt *(yji* xji         )
                    om(4,jat,2,iat)= tt *(zji* xji         )
                    om(2,jat,3,iat)= tt *(xji* yji         )
                    om(3,jat,3,iat)= tt *(yji* yji - .5d0/r)
                    om(4,jat,3,iat)= tt *(zji* yji         )
                    om(2,jat,4,iat)= tt *(xji* zji         )
                    om(3,jat,4,iat)= tt *(yji* zji         )
                    om(4,jat,4,iat)= tt *(zji* zji - .5d0/r)
              enddo
           enddo

  end subroutine OLDcreate_om_4

    subroutine mltampl_4(nat,amplitude,om)
    implicit real*8 (a-h,o-z)
    dimension amplitude(nat),om(4,nat,4,nat)
          do i=1,nat
          do j=1,nat
          do il=1,4
          do jl=1,4
          om(il,i,jl,j)=om(il,i,jl,j)*amplitude(i)*amplitude(j)
          enddo
          enddo
          enddo
          enddo
      end subroutine mltampl_4

    subroutine mltampl_1(nat,amplitude,om)
    implicit real*8 (a-h,o-z)
    dimension amplitude(nat),om(nat,nat)
          do i=1,nat
!          write(*,'(10(1x,e9.2))') (om(i,j),j=1,nat)
          do j=1,nat
!          write(*,'(i4,i4,3(e14.7))') i,j,om(i,j),amplitude(i),amplitude(j)
          om(i,j)=om(i,j)*amplitude(i)*amplitude(j)
          enddo
          enddo
      end subroutine mltampl_1


!=====================================================================
subroutine fingerprint_freebc(nat,nid,alat,geocode,rcov,rxyzIn,fp)
    !calculates an overlap matrix for atom centered GTO of the form:
    !s-type: 1/norm_s  exp(-(1/2)*(r/rcov)**2)
    !px type: 1/norm_p exp(-(1/2)*(r/rcov)**2) x/r  and analageously
    !         for py and pz
    use SPREDbase
    implicit none
    !parameters
    integer, intent(in)          :: nat
    integer, intent(in)          :: nid
    real(gp), intent(in)         :: alat(3)
    character(len=1), intent(in) :: geocode
    real(gp),intent(in)          :: rcov(nat)
    real(gp),intent(in)          :: rxyzIn(3,nat)
    real(gp),intent(out)         :: fp(nid)
    !internal
    integer  :: info
    integer  :: lwork
    integer  :: i1,i2,i3, n1, n2, n3
    integer  :: igto,jgto, iat, jat
    real(gp) :: tau(3)
    real(gp) :: cutoff, d2, r
    real(gp) :: sji, xi,yi,zi, xji, yji, zji, tt
    real(gp), parameter :: sqrt8=sqrt(8.0_gp)
    real(gp), allocatable, dimension(:,:) :: om
    real(gp), allocatable, dimension(:) :: workf
    real(gp), allocatable, dimension(:,:) :: rxyz

    ! WARNING! check convergence to ensure that the following
    !cutoff is large enough
    !! exp(-0.5*cutoff^2/rcov^2) = 1E-16  ==>
    !! cutoff^2 = 2*16*log(10)*rcov^2 ==> cutoff ~=8.5 rcov
    cutoff=9*maxval(rcov)

    rxyz = f_malloc((/3,nat/),id='rxyz')
    !with these settings the fingerprints have about 9 correct
    !decimal places
    if (geocode == 'F') then       ! free boundary conditions
        n1=0 ; n2=0 ; n3=0
        do iat=1,nat
            rxyz(1,iat)=rxyzIn(1,iat)
            rxyz(2,iat)=rxyzIn(2,iat)
            rxyz(3,iat)=rxyzIn(3,iat)
        enddo
    else if (geocode == 'S') then  ! surface boundary conditions,
                                   !non-periodic direction i s
        n1=nint(cutoff/alat(1))
        n2=0
        n3=nint(cutoff/alat(3))
        do iat=1,nat
            rxyz(1,iat)=modulo(rxyzIn(1,iat),alat(1))
            rxyz(2,iat)=rxyzIn(2,iat)
            rxyz(3,iat)=modulo(rxyzIn(3,iat),alat(3))
        enddo
    else if (geocode == 'P') then  ! periodic boundary conditions
        n1=nint(cutoff/alat(1))
        n2=nint(cutoff/alat(2))
        n3=nint(cutoff/alat(3))
        do iat=1,nat
            rxyz(1,iat)=modulo(rxyzIn(1,iat),alat(1))
            rxyz(2,iat)=modulo(rxyzIn(2,iat),alat(2))
            rxyz(3,iat)=modulo(rxyzIn(3,iat),alat(3))
        enddo
     else
        stop 'unrecognized BC in fingerprint'
     endif
     if (n1+n2+n3.gt.30) write(*,*) 'Warning n1,n2,n3 too big ',&
                                     n1,n2,n3

    if(nid .ne. nat .and. nid .ne. 4*nat)&
    stop ' nid should be either nat or  4*nat '

    om = f_malloc((/nid,nid/),id='om')
    om(:,:)=0.0_gp

    do i1=-n1,n1
        do i2=-n2,n2
            do i3=-n3,n3

               tau(1)=alat(1)*i1
               tau(2)=alat(2)*i2
               tau(3)=alat(3)*i3

               ! Gaussian overlap
               !  <sj|si>
               do iat=1,nat
                   xi=rxyz(1,iat) + tau(1)
                   yi=rxyz(2,iat) + tau(2)
                   zi=rxyz(3,iat) + tau(3)

                   do jat=iat,nat
                         d2=(rxyz(1,jat) -xi)**2 +(rxyz(2,jat)-yi)**2&
                            +(rxyz(3,jat)-zi)**2
                         r=.5_gp/(rcov(iat)**2 + rcov(jat)**2)
                         om(jat,iat)=om(jat,iat) + sqrt(4.0_gp*r*&
                        (rcov(iat)*rcov(jat)))**3 * exp(-d2*r)
                    enddo
                enddo
            enddo !i3
        enddo !i2
    enddo !i1

    !!  so far only s-s have been calculated
    if(nid == 4*nat) then  ! both s and p (nid = 4nat)

    do i1=-n1,n1
        do i2=-n2,n2
            do i3=-n3,n3
               tau(1)=alat(1)*i1
               tau(2)=alat(2)*i2
               tau(3)=alat(3)*i3
               !  <s|p>
                do iat=1,nat
                    xi=rxyz(1,iat) + tau(1)
                    yi=rxyz(2,iat) + tau(2)
                    zi=rxyz(3,iat) + tau(3)
                    do jat=1,nat   ! NOTE: do not use  jat=iat,nat
                                   !because all elements are on the
                                   !same side of the diagonal
                        xji=rxyz(1,jat) - xi
                        yji=rxyz(2,jat) - yi
                        zji=rxyz(3,jat) - zi

                        d2=xji*xji + yji*yji + zji*zji
                        r=.5_gp/(rcov(jat)**2 + rcov(iat)**2)

                        sji= sqrt(4.0_gp*r*(rcov(jat)*rcov(iat)))**3&
                            * exp(-d2*r)

                        !  <pj|si>
                        tt= sqrt8 *rcov(jat)*r * sji

                        om(1+nat + (jat-1)*3 ,iat )=  om(1+nat +&
                                       (jat-1)*3 ,iat ) + tt * xji
                        om(2+nat + (jat-1)*3 ,iat )=  om(2+nat +&
                                       (jat-1)*3 ,iat ) + tt * yji
                        om(3+nat + (jat-1)*3 ,iat )=  om(3+nat +&
                                       (jat-1)*3 ,iat ) + tt * zji

                       !! !  <sj|pi> no need, because they are on
                       !! !  the other side of the diagonal of the
                       !! !  symmetric matrix
                       !!  tt=-sqrt8 *rcov(iat)*r * sji

                       !!om(jat, 1+nat + (iat-1)*3 )= &
                       !!   om(jat, 1+nat + (iat-1)*3 ) + tt * xji
                       !!om(jat, 2+nat + (iat-1)*3 )= &
                       !!   om(jat, 2+nat + (iat-1)*3 ) + tt * yji
                       !!om(jat, 3+nat + (iat-1)*3 )= &
                       !!   om(jat, 3+nat + (iat-1)*3 ) + tt * zji
                    enddo
                enddo
                ! <pj|pi>
                do iat=1,nat
                    xi=rxyz(1,iat) + tau(1)
                    yi=rxyz(2,iat) + tau(2)
                    zi=rxyz(3,iat) + tau(3)
                    do jat=iat,nat
                        xji=rxyz(1,jat) - xi
                        yji=rxyz(2,jat) - yi
                        zji=rxyz(3,jat) - zi

                        d2=xji*xji + yji*yji + zji*zji
                        r=.5_gp/(rcov(jat)**2 + rcov(iat)**2)

                        sji= sqrt(4.0_gp*r*(rcov(jat)*rcov(iat)))**3&
                             * exp(-d2*r)

                        igto=nat+1 +(iat-1)*3
                        jgto=nat+1 +(jat-1)*3

                        tt = -8.0_gp*rcov(iat)*rcov(jat) * r*r * sji

                        om(jgto   , igto  )=  om(jgto   , igto  ) +&
                                           tt *(xji* xji - .5_gp/r)
                        om(jgto   , igto+1)=  om(jgto   , igto+1) +&
                                           tt *(yji* xji         )
                        om(jgto   , igto+2)=  om(jgto   , igto+2) +&
                                           tt *(zji* xji         )
                        om(jgto+1 , igto  )=  om(jgto+1 , igto  ) +&
                                           tt *(xji* yji         )
                        om(jgto+1 , igto+1)=  om(jgto+1 , igto+1) +&
                                           tt *(yji* yji - .5_gp/r)
                        om(jgto+1 , igto+2)=  om(jgto+1 , igto+2) +&
                                           tt *(zji* yji         )
                        om(jgto+2 , igto  )=  om(jgto+2 , igto  ) +&
                                           tt *(xji* zji         )
                        om(jgto+2 , igto+1)=  om(jgto+2 , igto+1) +&
                                           tt *(yji* zji         )
                        om(jgto+2 , igto+2)=  om(jgto+2 , igto+2) +&
                                           tt *(zji* zji - .5_gp/r)
                    enddo
                enddo

            enddo  ! i3
        enddo  ! i2
    enddo  ! i1
    endif  ! both s and p

    lwork=max(1,3*nid-1)
    workf =  f_malloc((/lwork/),id='workf')
    call DSYEV('N','L',nid,om,nid,fp,workf,-1,info)
    if (info.ne.0) stop 'info query'
    lwork=nint(workf(1))
    call f_free(workf)

    workf =  f_malloc((/lwork/),id='workf')
    call DSYEV('N','L',nid,om,nid,fp,workf,lwork,info)
    if (info.ne.0) stop 'info'

    call f_free(om)
    call f_free(rxyz)
    call f_free(workf)
end subroutine fingerprint_freebc
!=====================================================================
subroutine fpdistance(inputs,nid,nat,fp1,fp2,d)
    use SPREDbase
    use SPREDtypes
    implicit none
    !parameters
    type(SPRED_inputs), intent(in) :: inputs
    integer, intent(in) :: nid
    integer, intent(in) :: nat
    real(gp), intent(in) :: fp1(nid), fp2(nid)
    real(gp), intent(out) :: d 
    !internal
    select case(trim(f_char(inputs%fp_method)))
      case('OMF_FP_METHOD','OMPOLD_FP_METHOD','OMSOLD_FP_METHOD')
        call fpdistance_omf(inputs,nid,nat,fp1,fp2,d)
      case('OMP_FP_METHOD')
        call fpdistance_omp(inputs,nid,nat,fp1,fp2,d)
      case DEFAULT
        call f_err_throw('Following FP method is unknown: '//trim(f_char(inputs%fp_method)))
    end select
end subroutine fpdistance
!=====================================================================
subroutine fpdistance_omp(inputs,nid,nat,fp1,fp2,d)
    use SPREDbase
    use SPREDtypes
    implicit none
    !parameters
    type(SPRED_inputs), intent(in) :: inputs
    integer, intent(in) :: nid
    integer, intent(in) :: nat
    real(gp), intent(in) :: fp1(inputs%fp_angmom*inputs%fp_natx_sphere,nat), fp2(inputs%fp_angmom*inputs%fp_natx_sphere,nat)
    real(gp), intent(out) :: d 
    !internal
    integer :: iat, jat, l
    real(gp) :: tt
    real(gp) :: cost(nat,nat)
    integer :: iassign(nat)


    do iat=1,nat
        do jat=1,nat
            tt=0.d0
            do l=1,inputs%fp_angmom*inputs%fp_natx_sphere
                tt=tt+(fp1(l,iat)-fp2(l,jat))**2
            enddo
            tt=sqrt(tt) !really taking the sqrt here?
            cost(iat,jat)=tt
        enddo
    enddo
    call apc(nat, cost, iassign, d)
end subroutine fpdistance_omp
!=====================================================================
subroutine fpdistance_omf(inputs,nid,nat,fp1,fp2,d)
    use SPREDbase
    use SPREDtypes
    implicit none
    !parameters
    type(SPRED_inputs), intent(in) :: inputs
    integer, intent(in) :: nid
    integer, intent(in) :: nat
    real(gp), intent(in) :: fp1(nid), fp2(nid)
    real(gp), intent(out) :: d
    !internal
    integer :: i

    d=0.0_gp
    do i=1,nid
        d = d + (fp1(i)-fp2(i))**2
    enddo
    d=sqrt(d/(inputs%fp_angmom*nat))
end subroutine fpdistance_omf
!=====================================================================
logical function equal(inputs,nid,nat,iproc,prefix,txt,en_delta,fp_delta,epot1,epot2,fp1,fp2)
    use SPREDbase
    use SPREDtypes
    implicit none
    !parameter
    type(SPRED_inputs), intent(in) :: inputs
    integer, intent(in) :: nid
    integer, intent(in) :: nat
    integer, intent(in) :: iproc
    character(len=*), intent(in) :: prefix
    real(gp), intent(in) :: epot1,epot2
    real(gp), intent(in) :: en_delta, fp_delta
    real(gp), intent(in) :: fp1(nid), fp2(nid)
    character(len=2), intent(in) :: txt
    !internal
    real(gp) :: d

    equal=.false.
    call fpdistance(inputs,nid,nat,fp1,fp2,d)
    if(iproc==0)write(*,'(a,1x,a,1x,es14.7,1x,es14.7)')trim(adjustl(prefix))//'ediff, fpdist ',txt,abs(epot1-epot2),d
    if (abs(epot1-epot2).lt.en_delta) then
        call fpdistance(inputs,nid,nat,fp1,fp2,d)
        if (d.lt.fp_delta) then ! identical
            equal=.true.
        endif
    endif
end function
SUBROUTINE APC(N,A,F,Z)
implicit none
! Modified by Ali Sadeghi to get real*8 matrix A(N,N) and converted to F90 
!
! SOLUTION OF THE LINEAR MIN-SUM ASSIGNMENT PROBLEM.
! HUNGARIAN METHOD. COMPLEXITY O(N**3).
!
! MEANING OF THE INPUT PARAMETERS:
! N      = NUMBER OF ROWS AND COLUMNS OF THE COST MATRIX.
! A(I,J) = COST OF THE ASSIGNMENT OF ROW  I  TO COLUMN  J .
! ON RETURN, THE INPUT PARAMETERS ARE UNCHANGED.
!
! MEANING OF THE OUTPUT PARAMETERS:
! F(I) = COLUMN ASSIGNED TO ROW  I .
! Z    = COST OF THE OPTIMAL ASSIGNMENT =
!      = A(1,F(1)) + A(2,F(2)) + ... + A(N,F(N)) .
!
!
! THE CODE IS BASED ON THE HUNGARIAN METHOD AS DESCRIBED BY
! LAWLER (COMBINATORIAL OPTIMIZATION : NETWORKS AND
! MATROIDS, HOLT, RINEHART AND WINSTON, NEW YORK, 1976).
! THE ALGORITHMIC PASCAL-LIKE DESCRIPTION OF THE CODE IS
! GIVEN IN G.CARPANETO, S.MARTELLO AND P.TOTH, ALGORITHMS AND
! CODES FOR THE ASSIGNMENT PROBLEM, ANNALS OF OPERATIONS
! RESEARCH 7, 1988.
!
! SUBROUTINE APC DETERMINES THE INITIAL DUAL AND PARTIAL
! PRIMAL SOLUTIONS AND THEN SEARCHES FOR AUGMENTING PATHS
! UNTIL ALL ROWS AND COLUMNS ARE ASSIGNED.
!
! MEANING OF THE MAIN INTERNAL VARIABLES:
! FB(J) = ROW ASSIGNED TO COLUMN  J .
! M     = NUMBER OF INITIAL ASSIGNMENTS.
! U(I)  = DUAL VARIABLE ASSOCIATED WITH ROW  I .
! V(J)  = DUAL VARIABLE ASSOCIATED WITH COLUMN  J .
!
! APC NEEDS THE FOLLOWING SUBROUTINES: INCR
!                                      INIT
!                                      PATH
!
! THIS WORK WAS SUPPORTED BY  C.N.R. , ITALY.
      INTEGER n
      REAL(8)  A(n,n),Z,U(n),V(n)
      integer F(N),FB(n), RC(n)
      INTEGER M,I,J
! SEARCH FOR THE INITIAL DUAL AND PARTIAL PRIMAL SOLUTIONS.
      CALL INIT(N,A,F,M,U,V,FB,RC)
      IF ( M .NE. N ) then 
! SOLUTION OF THE REDUCED PROBLEM.
      DO  I=1,N
        IF ( F(I) == 0 ) THEN 
! DETERMINATION OF AN AUGMENTING PATH STARTING FROM ROW  I .
        CALL PATH(N,A,I,F,J,U,V,FB,RC)
! ASSIGNMENT OF ROW  I  AND COLUMN  J .
        CALL INCR(n,F,J,FB,RC)
        ENDIF
      ENDDO    
      ENDIF

! COMPUTATION OF THE SOLUTION COST  Z .
      Z = sum(u(1:N)) + sum(V(1:N))
      END SUBROUTINE


      SUBROUTINE INCR(n,F,J,FB,RC)
!
! ASSIGNMENT OF COLUMN  J .
!
      INTEGER n,I,J,JJ,  F(n),FB(n),RC(n)
   10 I = RC(J)
      FB(J) = I
      JJ = F(I)
      F(I) = J
      J = JJ
      IF ( J > 0 ) GO TO 10
      RETURN
      END SUBROUTINE


      SUBROUTINE INIT(N,A,F,M,U,V,FB,P)
!
! SEARCH FOR THE INITIAL DUAL AND PARTIAL PRIMAL SOLUTIONS.
!
! P(I) = FIRST UNSCANNED COLUMN OF ROW  I .
       IMPLICIT NONE
!
      INTEGER n,m, F(n),FB(n),P(n)
      real(8) A(n,n) , U(n),V(n)
      REAL(8), parameter :: INF = 1.d9
      real(8) min, IA
      integer i,j, k,R, JMIN, KK
! PHASE 1 .
      M = 0
      F(1:N)=0
      FB(1:N)=0
! SCANNING OF THE COLUMNS ( INITIALIZATION OF  V(J) ).
      DO 40 J=1,N
        MIN = INF
        DO 30 I=1,N
          IA = A(I,J)
          IF ( IA .GT. MIN ) GO TO 30
          IF ( IA .LT. MIN ) GO TO 20
          IF ( F(I) .NE. 0 ) GO TO 30
   20     MIN = IA
          R = I
   30   CONTINUE
        V(J) = MIN
        IF ( F(R) .NE. 0 ) GO TO 40
! ASSIGNMENT OF COLUMN  J  TO ROW  R .
        M = M + 1
        FB(J) = R
        F(R) = J
        U(R) = 0.d0
        P(R) = J + 1
   40 CONTINUE
! PHASE 2 .
! SCANNING OF THE UNASSIGNED ROWS ( UPDATING OF  U(I) ).
      DO 110 I=1,N
        IF ( F(I) .NE. 0 ) GO TO 110
        MIN = INF
        DO 60 K=1,N
          IA = A(I,K) - V(K)
          IF ( IA .GT. MIN )  GO TO 60
          IF ( IA .LT. MIN )  GO TO 50
          IF ( FB(K) .NE. 0 ) GO TO 60
          IF ( FB(J) .EQ. 0 ) GO TO 60
   50     MIN = IA
          J = K
   60   CONTINUE
        U(I) = MIN
        JMIN = J
        IF ( FB(J) .EQ. 0 ) GO TO 100
        DO 80 J=JMIN,N
          IF ( A(I,J) - V(J) .GT. MIN ) GO TO 80
          R = FB(J)
          KK = P(R)
          IF ( KK .GT. N ) GO TO 80
          DO 70 K=KK,N
            IF ( FB(K) .GT. 0 ) GO TO 70
            IF ( A(R,K) - U(R) - V(K) .EQ. 0.d0 ) GO TO 90
   70     CONTINUE
          P(R) = N + 1
   80   CONTINUE
        GO TO 110
! REASSIGNMENT OF ROW  R  AND COLUMN  K .
   90   F(R) = K
        FB(K) = R
        P(R) = K + 1
! ASSIGNMENT OF COLUMN  J  TO ROW  I .
  100   M = M + 1
        F(I) = J
        FB(J)= I
        P(I) = J + 1
  110 CONTINUE
      RETURN
      END SUBROUTINE
      SUBROUTINE PATH(N,A,II,F,JJ,U,V,FB,RC)
!
! DETERMINATION OF AN AUGMENTING PATH STARTING FROM
! UNASSIGNED ROW  II  AND TERMINATING AT UNASSIGNED COLUMN
! JJ , WITH UPDATING OF DUAL VARIABLES  U(I)  AND  V(J) .
!
! MEANING OF THE MAIN INTERNAL VARIABLES:
! LR(L) = L-TH LABELLED ROW ( L=1,NLR ).
! PI(J) = MIN ( A(I,J) - U(I) - V(J) , SUCH THAT ROW  I  IS
!         LABELLED AND NOT EQUAL TO  FB(J) ).
! RC(J) = ROW PRECEDING COLUMN  J  IN THE CURRENT
!         ALTERNATING PATH.
! UC(L) = L-TH UNLABELLED COLUMN ( L=1,NUC ).
!
      implicit none
      INTEGER N 
      real(8)  A(n,n),U(n),V(N),PI(n), IA, MIN
      INTEGER F(N),LR(n),UC(n)
      INTEGER FB(n),RC(n)
      REAL(8), parameter :: INF = 1.d9
      integer  i,j,k,L, ii,jj, NUC , NLR, R
! INITIALIZATION.
      LR(1) = II
      DO 10 K=1,N
        PI(K) = A(II,K) - U(II) - V(K)
        RC(K) = II
        UC(K) = K
   10 CONTINUE
      NUC = N
      NLR = 1
      GO TO 40
! SCANNING OF THE LABELLED ROWS.
   20 R = LR(NLR)
      DO 30 L=1,NUC
        J = UC(L)
        IA = A(R,J) - U(R) - V(J)
        IF ( IA .GE. PI(J) ) GO TO 30
        PI(J) = IA
        RC(J) = R
   30 CONTINUE
! SEARCH FOR A ZERO ELEMENT IN AN UNLABELLED COLUMN.
   40 DO 50 L=1,NUC
        J = UC(L)
        IF ( PI(J) .EQ. 0.d0 ) GO TO 100
   50 CONTINUE
! UPDATING OF THE DUAL VARIABLES  U(I)  AND  V(J) .
      MIN = INF
      DO 60 L=1,NUC
        J = UC(L)
        IF ( MIN .GT. PI(J) ) MIN = PI(J)
   60 CONTINUE
      DO 70 L=1,NLR
        R = LR(L)
        U(R) = U(R) + MIN
   70 CONTINUE
      DO 90 J=1,N
        IF ( PI(J) .EQ. 0.d0 ) GO TO 80
        PI(J) = PI(J) - MIN
        GO TO 90
   80   V(J) = V(J) - MIN
   90 CONTINUE
      GO TO 40
  100 IF ( FB(J) .EQ. 0 ) GO TO 110
! LABELLING OF ROW  FB(J)  AND REMOVAL OF THE LABEL  OF COLUMN  J .
      NLR = NLR + 1
      LR(NLR) = FB(J)
      UC(L) = UC(NUC)
      NUC = NUC - 1
      GO TO 20
! DETERMINATION OF THE UNASSIGNED COLUMN  J .
  110 JJ = J
      RETURN
      END SUBROUTINE
 

end module
