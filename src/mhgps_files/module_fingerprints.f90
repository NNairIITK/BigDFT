!> @file
!! module implementing the fingerprints
!!     
!! @author Bastian Schaefer
!! @section LICENCE
!!    Copyright (C) 2014 UNIBAS
!!    This file is not freely distributed.
!!    A licence is necessary from UNIBAS

module module_fingerprints
    implicit none

    private

    public :: fingerprint
    public :: fpdistance
    public :: equal

    contains

!=====================================================================
subroutine fingerprint(nat,nid,alat,geocode,rcov,rxyz,fp)
    !calculates an overlap matrix for atom centered GTO of the form:
    !s-type: 1/norm_s  exp(-(1/2)*(r/rcov)**2)
    !px type: 1/norm_p exp(-(1/2)*(r/rcov)**2) x/r  and analageously
    !         for py and pz
    use module_base
    implicit none
    !parameters
    integer, intent(in)          :: nat
    integer, intent(in)          :: nid
    real(gp), intent(in)         :: alat(3)
    character(len=1), intent(in) :: geocode
    real(gp),intent(in)          :: rcov(nat)
    real(gp),intent(in)          :: rxyz(3,nat)
    real(gp),intent(out)         :: fp(nid)
    !internal
    integer  :: info
    integer  :: i1,i2,i3, n1, n2, n3  
    integer  :: igto,jgto, iat, jat
    real(gp) :: tau(3)
    real(gp) :: cutoff, d2, r
    real(gp) :: sji, xi,yi,zi, xji, yji, zji, tt 
    real(gp), parameter :: sqrt8=sqrt(8.0_gp)
    real(gp), allocatable, dimension(:,:) :: om,workf

    ! WARNING! check convergence to ensure that the folloing
    !cutoff is large enough
    !! exp(-0.5*cutoff^2/rcov^2) = 1E-16  ==>
    !! cutoff^2 = 2*16*log(10)*rcov^2 ==> cutoff ~=8.5 rcov 
    cutoff=9*maxval(rcov)

    !with these settings the fingerprints have about 9 correct
    !decimal places
    if (geocode == 'F') then       ! free boundary conditions
        n1=0 ; n2=0 ; n3=0
    else if (geocode == 'S') then  ! surface boundary conditions,
                                   !non-periodic direction i s
        n1=nint(cutoff/alat(1))
        n2=0
        n3=nint(cutoff/alat(3))
    else if (geocode == 'P') then  ! periodic boundary conditions
        n1=nint(cutoff/alat(1))
        n2=nint(cutoff/alat(2))
        n3=nint(cutoff/alat(3))
     else
        stop 'unrecognized BC in fingerprint'
     endif
     if (n1+n2+n3.gt.30) write(*,*) 'Warning n1,n2,n3 too big ',&
                                     n1,n2,n3

    if(nid .ne. nat .and. nid .ne. 4*nat)&
    stop ' nid should be either nat or  4*nat '

    om = f_malloc((/nid,nid/),id='om')
    workf =  f_malloc((/nid,nid/),id='workf')
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

    call DSYEV('N','L',nid,om,nid,fp,workf,nid**2,info)
    if (info.ne.0) stop 'info'

    call f_free(om)
    call f_free(workf)
end subroutine fingerprint
!=====================================================================
subroutine fpdistance(nid,fp1,fp2,d)
    use module_base
    implicit none
    !parameters
    integer, intent(in) :: nid
    real(gp), intent(in) :: fp1(nid), fp2(nid)
    real(gp), intent(out) :: d
    !internal
    integer :: i

    d=0.0_gp
    do i=1,nid
        d = d + (fp1(i)-fp2(i))**2
    enddo
    d=sqrt(d/nid)
end subroutine fpdistance
!=====================================================================
logical function equal(nid,en_delta,fp_delta,epot1,epot2,fp1,fp2)
    use module_base
    use module_global_variables, only: iproc
    implicit none
    !parameter
    integer, intent(in) :: nid
    real(gp), intent(in) :: epot1,epot2
    real(gp), intent(in) :: en_delta, fp_delta
    real(gp), intent(in) :: fp1(nid), fp2(nid)
    !internal
    real(gp) :: d=1.d100

    equal=.false.
    call fpdistance(nid,fp1,fp2,d)
    if(iproc==0)write(*,*)'(MHGPS) ediff, fpdist',abs(epot1-epot2),d
    if (abs(epot1-epot2).lt.en_delta) then
        call fpdistance(nid,fp1,fp2,d)
        if (d.lt.fp_delta) then ! identical
            equal=.true.
        endif
    endif
end function
!=====================================================================
end module
