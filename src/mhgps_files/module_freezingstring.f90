!> @file
!! module implementing the freezing string technique
!!     
!! @author Bastian Schaefer
!! @section LICENCE
!!    Copyright (C) 2014 UNIBAS
!!    This file is not freely distributed.
!!    A licence is necessary from UNIBAS
module module_freezingstring
    implicit none
    
    contains
!TODO: create function "get_input_guess" which returns
!an inpute guess for a ts (consisting of corrds. and minmode dir.)
!this method simpliy call grow_string
!and then uses cubic splines to find the tangent at the highest energy node
!!=====================================================================
!subroutine grow_string(nat,alat,step,gammainv,perpnrmtol,trust,nstepsmax,nstringmax,nstring,string)
!    use module_base
!    implicit none
!    !parameters
!    integer, intent(in) :: nat,nstepsmax
!    integer, intent(inout) :: nstringmax,nstring
!    real(gp) , intent(in) :: step,gammainv,perpnrmtol,trust,alat(3)
!    real(gp), allocatable, intent(inout) :: string(:,:,:)
!    !constants
!    character(len=10), parameter :: method = 'linlst'
!    !internal
!    integer :: i,j,k,istart
!    integer, parameter :: resize=10
!    real(gp), allocatable :: stringTMP(:,:,:)
!    real(gp) :: tangent(3*nat)
!    real(gp) :: tangentright(3*nat)
!    real(gp) :: dist,dir(3*nat),stepint,perpnrmtol_squared,trust_squared
!    
!    real(gp) :: dist1,dist2
!    integer :: finished=1
!    integer :: nresizes
!         !if finished >0: not finished
!         !if finished =0: finished
!         !if finished <0: one more node in the middle
!
!    if(.not. allocated(string))then
!        write(*,*)'STOP, string not allocated'
!        stop
!    endif
!    perpnrmtol_squared=perpnrmtol**2
!    trust_squared=trust**2
!
!!write(200,*)'#######################'
!
!    nstring=1
!    dir = string(:,2,1)-string(:,1,1)
!    dist = sqrt(dot_product(dir,dir))
!  stepint=0.1_gp*dist
!!  stepint=step
!!write(333,*)dist
!    if(dist<=2._gp*stepint)then!interpolation done
!        return
!    endif
!    nresizes=0
!    do while (.not. finished .or. nresizes>=100)!maximum 100 resizes of string array
!        istart=nstring
!        do i=istart,nstringmax-1
!!            dir = string(:,2,i)-string(:,1,i)
!!            dist = sqrt(dot_product(dir,dir))
!!            stepint=0.05_gp*dist
!            call interpol(method,nat,string(1,1,nstring),&
!                 string(1,2,nstring),stepint,string(1,1,nstring+1),&
!                 string(1,2,nstring+1),tangent,tangentright,finished)
!            dir = string(:,1,nstring+1)-string(:,2,nstring)
!            dist1 = sqrt(dot_product(dir,dir))
!            dir = string(:,1,nstring+1)-string(:,1,nstring)
!            dist2 = sqrt(dot_product(dir,dir))
!            if(dist1<dist2)then
!            dir = string(:,2,nstring)-string(:,1,nstring)
!            dist = sqrt(dot_product(dir,dir))
!            write(*,*)dist,dist1,dist2,j,i,stepint,nstring
!                stop'dist1<dist2'
!            endif
!            call optim_cg(nat,alat,stepint,gammainv,perpnrmtol_squared,trust_squared,nstepsmax,tangent,string(1,1,i+1),string(1,2,i+1))
!            dir = string(:,2,i+1)-string(:,1,i+1)
!            dist = sqrt(dot_product(dir,dir))
!            nstring=nstring+1
!!write(*,*)'nstring',nstring,nstringmax-1
!!write(*,*)'dist',dist
!            if(dist<=2._gp*stepint)then!interpolation done
!                return
!            endif
!        enddo
!        nresizes=nresizes+1
!        if(allocated(stringTmp))then
!            deallocate(stringTmp)
!        endif
!        allocate(stringTmp(3*nat,2,nstringmax))
!        stringTmp=string
!        deallocate(string)
!        nstringmax=nstringmax+resize
!        allocate(string(3*nat,2,nstringmax))
!!write(*,*)'resize',j
!        do k=1,(nstringmax-resize)
!            string(:,:,k)=stringTmp(:,:,k)
!        enddo
!        deallocate(stringTmp)
!    enddo
!
!end subroutine
!!=====================================================================
!subroutine optim_cg(nat,alat,step,gammainv,perpnrmtol_squared,trust_squared,nstepsmax,tangent,rxyz1,rxyz2)
!    use module_base
!    implicit none
!    !parameters
!    integer, intent(in) :: nat,nstepsmax
!    real(gp), intent(in) :: tangent(3*nat),step,gammainv,perpnrmtol_squared,trust_squared,alat(3)
!    real(gp), intent(inout) :: rxyz1(3*nat),rxyz2(3*nat)
!    !internal
!    real(gp) :: fxyz1(3*nat),fxyz2(3*nat)
!    real(gp) :: perp1(3*nat),perp2(3*nat)
!    real(gp) :: epot1, epot2
!    real(gp) :: d0,dmax,dist,dir(3*nat)
!    real(gp) :: dispPrev1(3*nat),disp1(3*nat)
!    real(gp) :: dispPrev2(3*nat),disp2(3*nat)
!    real(gp) :: alpha1,alpha2
!    real(gp) :: fnrmPrev1, fnrm1
!    real(gp) :: fnrmPrev2, fnrm2
!    real(gp) :: perpnrmPrev1_squared, perpnrm1_squared
!    real(gp) :: perpnrmPrev2_squared, perpnrm2_squared
!    real(gp) :: dispnrm_squared
!    integer :: istep
!
!    real(gp) :: tangentint(3*nat),tmp
!
!    tangentint=tangent
!
!
!    dir=rxyz2-rxyz1
!    d0=sqrt(dot_product(dir,dir))
!    dmax=d0+0.5_gp*step
!    call energyandforces(nat, alat, rxyz1,fxyz1,epot1,'cnt_enf_freezing_string')
!    call energyandforces(nat, alat, rxyz2,fxyz2,epot2,'cnt_enf_freezing_string')
!
!    fnrmPrev1=sqrt(dot_product(fxyz1,fxyz1))
!    fnrmPrev2=sqrt(dot_product(fxyz2,fxyz2))
!
!    call perpend(nat,tangentint,fxyz1,perp1)
!    call perpend(nat,tangentint,fxyz2,perp2)
!
!    perpnrmPrev1_squared = dot_product(perp1,perp1)
!    perpnrmPrev2_squared = dot_product(perp2,perp2)
!    perpnrm1_squared=perpnrmPrev1_squared
!    perpnrm2_squared=perpnrmPrev2_squared
!
!    !first steps: steepest descent
!    dispPrev1=gammainv*perp1
!    dispnrm_squared=maxval(dispPrev1**2)
!    if(dispnrm_squared > trust_squared)then
!        dispPrev1=dispPrev1*sqrt(trust_squared/dispnrm_squared)
!    endif
!    rxyz1=rxyz1+dispPrev1
!    dispPrev2=gammainv*perp2
!    dispnrm_squared=maxval(dispPrev2**2)
!    if(dispnrm_squared > trust_squared)then
!        dispPrev2=dispPrev2*sqrt(trust_squared/dispnrm_squared)
!    endif
!    rxyz2=rxyz2+dispPrev2
!
!    call energyandforces(nat, alat, rxyz1,fxyz1,epot1,'cnt_enf_freezing_string')
!    call energyandforces(nat, alat, rxyz2,fxyz2,epot2,'cnt_enf_freezing_string')
!    fnrm1=sqrt(dot_product(fxyz1,fxyz1))
!    fnrm2=sqrt(dot_product(fxyz2,fxyz2))
!    
!    dir=rxyz2-rxyz1
!    dist=sqrt(dot_product(dir,dir))
!    if(dist>dmax)then
!        return
!    endif
!
!!write(200,*)'-------------------'
!    !other steps: cg
!    do istep=2,nstepsmax
!
!    dir=rxyz2-rxyz1
!    dist=sqrt(dot_product(dir,dir))
!    tangentint=dir/dist
!        !move left node
!        call perpend(nat,tangentint,fxyz1,perp1)
!        perpnrm1_squared = dot_product(perp1,perp1)
!!        if(perpnrm1_squared>=perpnrmtol_squared)then
!!!            if(fnrm1>fnrmPrev1)then
!            if(perpnrm1_squared>perpnrmPrev1_squared)then
!                alpha1=1._gp
!            else
!                alpha1 = perpnrm1_squared / perpnrmPrev1_squared
!            endif
!            disp1=gammainv*perp1+ alpha1 * dispPrev1
!            dispnrm_squared=maxval(disp1**2)
!!write(*,*)'trust1',sqrt(dispnrm_squared),sqrt(trust_squared)
!            if(dispnrm_squared > trust_squared)then
!!write(200,*)'Limiting step length!'
!!                disp1 = sqrt(trust_squared) * disp1/sqrt(dot_product(disp1,disp1))
!                 disp1=disp1*sqrt(trust_squared/dispnrm_squared)
!            endif
!            rxyz1=rxyz1+disp1
!            dispPrev1=disp1
!            perpnrmPrev1_squared=perpnrm1_squared
!            fnrmPrev1=fnrm1
!            call energyandforces(nat, alat, rxyz1,fxyz1,epot1,'cnt_enf_freezing_string')
!            fnrm1=sqrt(dot_product(fxyz1,fxyz1))
!!        endif
!
!        !move right node
!            call energyandforces(nat, alat, rxyz2,fxyz2,epot2,'cnt_enf_freezing_string')
!            fnrm2=sqrt(dot_product(fxyz2,fxyz2))
!            call perpend(nat,tangentint,fxyz2,perp2)
!            perpnrm2_squared = dot_product(perp2,perp2)
!!        if(perpnrm2_squared>=perpnrmtol_squared)then
!!            if(fnrm2>fnrmPrev2)then
!            if(perpnrm2_squared>perpnrmPrev2_squared)then
!                alpha2=1._gp
!            else
!                alpha2 = perpnrm2_squared / perpnrmPrev2_squared
!            endif
!            disp2=gammainv*perp2+ alpha2 * dispPrev2
!            dispnrm_squared=maxval(disp2**2)
!!write(*,*)'trust2',sqrt(dispnrm_squared),sqrt(trust_squared)
!            if(dispnrm_squared > trust_squared)then
!!write(200,*)'Limiting step length!'
!!                disp2 = sqrt(trust_squared) * disp2/sqrt(dot_product(disp2,disp2))
!                 disp2=disp2*sqrt(trust_squared/dispnrm_squared)
!            endif
!!write(*,*)'trust2',sqrt(dispnrm_squared),sqrt(trust_squared)
!!write(*,*)'disp2',sqrt(dot_product(disp2,disp2))
!            rxyz2=rxyz2+disp2
!            dispPrev2=disp2
!            perpnrmPrev2_squared=perpnrm2_squared
!            fnrmPrev2=fnrm2
!!        endif
!    
!!write(*,*)alpha1,alpha2
!!write(200,*)'perp',sqrt(perpnrm1_squared), sqrt(perpnrm2_squared)
!!write(*,*)tangent
!!write(*,*)rxyz2
!!write(*,*)perpnrm2_squared
!        dir=rxyz2-rxyz1
!        dist=sqrt(dot_product(dir,dir))
!        if(dist>dmax.or. (perpnrm1_squared<perpnrmtol_squared &
!        !if((perpnrm1_squared<perpnrmtol_squared &
!                    &.and. perpnrm2_squared<perpnrmtol_squared))then
!            if(dist>dmax)then
!                 write(200,*)'exit due to dmax'   
!            endif
!            return
!        endif
!
!
!    enddo
!end subroutine
!=====================================================================
subroutine perpend(nat,tangent,fxyz,perp)
    use module_base
    !returns a vector perp that contains
    !the perpendicular components of fyxz to
    !the tangent vector
    !that is: all components of fxyz in direction
    !of tangent are substracted from xyz
    implicit none
    !parameters
    integer, intent(in) :: nat
    real(gp), intent(in) :: tangent(3*nat),fxyz(3*nat)
    real(gp), intent(out) :: perp(3*nat)

    perp = fxyz - dot_product(tangent,fxyz)*tangent
    
end subroutine
!=====================================================================
subroutine lin_interpol(nat,left, right, step,interleft,interright,&
                       tangent, finished)
    use module_base
    implicit none
    !parameters
    integer, intent(in)  :: nat
    real(gp), intent(in)  :: left(3*nat)
    real(gp), intent(in)  :: right(3*nat)
    real(gp), intent(in)  :: step
    real(gp), intent(out) :: interleft(3*nat)
    real(gp), intent(out) :: interright(3*nat)
    real(gp), intent(out) :: tangent(3*nat)
    integer, intent(out) :: finished
    !internal

    !tangent points from left to right:    
    tangent = right-left
    tangent = tangent / sqrt(dot_product(tangent,tangent))
    interleft = left+step*tangent
    interright = right - step*tangent
end subroutine
!=====================================================================
subroutine lst_interpol(nat,left,right,step,interleft,interright,&
                        tangentleft,tangentright,finished)
    !on return:
    !if finished > 0: interleft and intergiht contain the new nodes
    !if finished < 0: left and right are two close, only one new node
    !                 is returned in interleft.
    !                 interright is meaningless in this case.
    !if finished = 0: left and right are closer than 'step'
    !                 => freezing string search finsihed
    !                 nothing is returned, interleft and interright
    !                 are meaningless 
    use module_base
    use module_interpol
    implicit none
    !parameters
    integer, intent(in)  :: nat
    real(gp), intent(in)  :: left(3,nat)
    real(gp), intent(in)  :: right(3,nat)
    real(gp), intent(inout)  :: step
    real(gp), intent(out) :: interleft(3,nat)
    real(gp), intent(out) :: interright(3,nat)
    real(gp), intent(out) :: tangentleft(3,nat)
    real(gp), intent(out) :: tangentright(3,nat)
    integer, intent(out) :: finished
    !constants
    integer, parameter :: nimages=200
    real(gp), parameter :: stepfrct=0.2_gp
    !internal
    integer :: i,j,tnat
    integer :: iinterleft
    integer :: iinterright
    real(gp) :: lstpath(3,nat,nimages)
    real(gp) :: arc(nimages), arcl,arclh, pthl
    real(gp) :: delta,deltaold
    real(gp) :: diff(3,nat)
    real(gp) :: nimo
    real(gp) :: yp1=huge(1._gp), ypn=huge(1._gp)!natural splines
    real(gp) :: y2vec(3,nat,nimages)
    real(gp) :: tau, ydmy
    real(gp) :: lambda
    !functions
    real(gp) :: dnrm2
    tnat=3*nat

    !create high density lst path
    nimo=1._gp/real(nimages-1,gp)
    do i=1,nimages
        lambda  = real(i-1,gp)*nimo
        call lstpthpnt(nat,left,right,lambda,lstpath(:,:,i))
    enddo
    
    !measure arc length 
    arc(1)=0._gp
    do i=2,nimages
        diff = lstpath(:,:,i) - lstpath(:,:,i-1)
        arc(i)  = arc(i-1) + dnrm2(tnat,diff,1)
    enddo
    arcl=arc(nimages)
    if(step<0._gp)step=stepfrct*arcl

    if(arcl < step)then
        finished=0
        return
    endif

    !compute the spline parameters (y2vec)
    !parametrize curve as a function of the
    !integrated arc length
    do i=1,nimages
        do j=1,nat
            !potentially performance issues since lstpath
            !is not transversed in column-major order in
            !spline_wrapper:
            call spline_wrapper(arc,lstpath(1,j,1),nimages,&
                               yp1,ypn,y2vec(1,j,1))
            call spline_wrapper(arc,lstpath(2,j,1),nimages,&
                               yp1,ypn,y2vec(2,j,1))
            call spline_wrapper(arc,lstpath(3,j,1),nimages,&
                               yp1,ypn,y2vec(3,j,1))
        enddo
    enddo

    if(arcl < 2._gp*step)then!only one more point
        !we have to return the point in the 'middle'    
        tau = 0.5_gp*arcl
        do j=1,nat
            !potentially performance issues since lstpath
            !is not transversed in column-major order in
            !splint_wrapper
            call splint_wrapper(arc,lstpath(1,j,1),y2vec(1,j,1),&
                 nimages,tau,interleft(1,j),tangentleft(1,j))
            call splint_wrapper(arc,lstpath(2,j,1),y2vec(2,j,1),&
                 nimages,tau,interleft(2,j),tangentleft(2,j))
            call splint_wrapper(arc,lstpath(3,j,1),y2vec(3,j,1),&
                 nimages,tau,interleft(3,j),tangentleft(3,j))
        enddo
        finished=-1
    else! standard case
        !we have to return the two points interleft
        !and interright whose distances to left and right
        !are roughly 'step'

        !first left...
        tau = step
        do j=1,nat
            !potentially performance issues since lstpath
            !is not transversed in column-major order in
            !splint_wrapper
            call splint_wrapper(arc,lstpath(1,j,1),y2vec(1,j,1),&
                 nimages,tau,interleft(1,j),tangentleft(1,j))
            call splint_wrapper(arc,lstpath(2,j,1),y2vec(2,j,1),&
                 nimages,tau,interleft(2,j),tangentleft(2,j))
            call splint_wrapper(arc,lstpath(3,j,1),y2vec(3,j,1),&
                 nimages,tau,interleft(3,j),tangentleft(3,j))
        enddo

        !...then right
        tau = arcl-step
        do j=1,nat
            !potentially performance issues since lstpath
            !is not transversed in column-major order in
            !splint_wrapper
            call splint_wrapper(arc,lstpath(1,j,1),y2vec(1,j,1),&
                 nimages,tau,interright(1,j),tangentright(1,j))
            call splint_wrapper(arc,lstpath(2,j,1),y2vec(2,j,1),&
                 nimages,tau,interright(2,j),tangentright(2,j))
            call splint_wrapper(arc,lstpath(3,j,1),y2vec(3,j,1),&
                 nimages,tau,interright(3,j),tangentright(3,j))
        enddo
        finished=2
    endif
end subroutine
!=====================================================================
subroutine interpol(method,nat,left,right,step,interleft,interright,&
                    tangentleft,tangentright,finished)
    use module_base
    implicit none
    !parameters
    character(len=*), intent(in) :: method
    integer, intent(in)  :: nat
    real(gp), intent(in)  :: left(3*nat)
    real(gp), intent(in)  :: right(3*nat)
    real(gp), intent(inout)  :: step
    real(gp), intent(out) :: interleft(3*nat)
    real(gp), intent(out) :: interright(3*nat)
    real(gp), intent(out) :: tangentleft(3*nat)
    real(gp), intent(out) :: tangentright(3*nat)
    integer, intent(out) :: finished

    if(trim(adjustl(method))=='lincat')then
        call lin_interpol(nat,left, right, step,interleft,interright,&
                       tangentleft,finished)
        tangentright=tangentleft
    else if(trim(adjustl(method))=='linlst')then
        call lst_interpol(nat,left, right, step,interleft,interright,&
                       tangentleft,tangentright,finished)
    endif
end subroutine
!=====================================================================
subroutine spline_wrapper(xvec,yvec,ndim,yp1,ypn,y2vec)
    !routine for initializing the spline vectors
    !xvec[1..ndim] and yvec[1..ndim] contain the tabulated function.
    !yi= f(xi), x1 < x2 < ... < xN .
    !yp1, ypn: values of first derivative of the interpolating
    !function at points 1 and ndim
    !y2vec: second derivatives of the interpolating function at the
    !tabulated points
    use module_base
    implicit none
    !parameters
    integer, intent(in)  :: ndim
    real(gp), intent(in) :: xvec(ndim), yvec(ndim)
    real(gp), intent(in) :: yp1, ypn
    real(gp), intent(out) :: y2vec(ndim)
    !internal
    real(gp) :: xt(ndim), yt(ndim)
    real(gp) :: ytp1, ytpn
    if(xvec(1).eq.xvec(ndim)) then
        y2vec=0.0_gp
    elseif(xvec(1).gt.xvec(2)) then
        xt=-xvec
        yt=yvec
        ytp1=-yp1
        ytpn=-ypn
        call spline(xt,yt,ndim,ytp1,ytpn,y2vec)
    else
        call spline(xvec,yvec,ndim,yp1,ypn,y2vec)
    endif
end subroutine
!=====================================================================
subroutine splint_wrapper(xvec,yvec,y2vec,ndim,tau,yval,dy)
    !xvec[1..ndim] and yvec[1..ndim] contain the tabulated function.
    !yi= f(xi), x1 < x2 < ... < xN .
    !y2vec: second derivatives of the interpolating function at the
    !tabulated points
    !tau: spline's parameter
    !yval: cubic spline interpolation value at tay
    !dy: derivative of spline at tau (with respect to 
    !    the parametrization
    use module_base
    implicit none
    !parameters
    integer, intent(in)  :: ndim
    real(gp), intent(in) :: xvec(ndim), yvec(ndim)
    real(gp), intent(in) :: tau
    real(gp), intent(out) :: yval, dy
    real(gp), intent(in) :: y2vec(ndim)
    !internal
    real(gp) :: xt(ndim), yt(ndim), taut
    if(xvec(1).eq.xvec(ndim)) then
        yval=yvec(1)
        dy=0.0_gp
    elseif(xvec(1).gt.xvec(2)) then
        xt=-xvec
        yt=yvec
        taut=-tau
        call splint(xt,yt,y2vec,ndim,taut,yval,dy)
        dy=-dy
    else
        call splint(xvec,yvec,y2vec,ndim,tau,yval,dy)
    endif
end subroutine
!=====================================================================
subroutine spline(xvec,yvec,ndim,yp1,ypn,y2vec)
    !translated to f90 from numerical recipes
    use module_base
    implicit none
    !parameter
    integer, intent(in) :: ndim
    real(gp), intent(in) :: xvec(ndim), yvec(ndim)
    real(gp), intent(in) :: yp1, ypn
    real(gp), intent(out) :: y2vec(ndim)
    !internal
    integer  :: i,k
    real(gp) :: p,qn,sig,un,work(ndim)
    if (yp1 > .99e30_gp) then
        y2vec(1)=0.0_gp
        work(1)=0.0_gp
    else
        y2vec(1)=-0.5_gp
        work(1)=(3./(xvec(2)-xvec(1)))*((yvec(2)-yvec(1))/&
                (xvec(2)-xvec(1))-yp1)
    endif
    do i=2,ndim-1
        sig=(xvec(i)-xvec(i-1))/(xvec(i+1)-xvec(i-1))
        p=sig*y2vec(i-1)+2.0_gp
        y2vec(i)=(sig-1.0_gp)/p
        work(i)=(6.0_gp*((yvec(i+1)-yvec(i))/(xvec(i+1)-xvec(i))-&
                (yvec(i)-yvec(i-1))/(xvec(i)-xvec(i-1)))/&
                (xvec(i+1)-xvec(i-1))-sig*work(i-1))/p  
    enddo
    if(ypn>.99e30_gp) then
        qn=0.0_gp
        un=0.0_gp
    else
        qn=0.5_gp
        un=(3.0_gp/(xvec(ndim)-xvec(ndim-1)))*&
           (ypn-(yvec(ndim)-yvec(ndim-1))/(xvec(ndim)-xvec(ndim-1)))
    endif
    y2vec(ndim)=(un-qn*work(ndim-1))/(qn*y2vec(ndim-1)+1.0_gp)
    do k=ndim-1,1,-1
        y2vec(k)=y2vec(k)*y2vec(k+1)+work(k)
    enddo
end subroutine
!=====================================================================
subroutine splint(xvec,yvec,y2vec,ndim,tau,yval,dy)
    !translated to f90 from numerical recipes
    use module_base
    implicit none
    !parameters
    integer, intent(in) :: ndim
    real(gp), intent(in) :: xvec(ndim), yvec(ndim)
    real(gp), intent(in) :: y2vec(ndim)
    real(gp), intent(in)  :: tau
    real(gp), intent(out) :: yval, dy
    !internal
    integer :: k,khi,klo
    real(gp):: a,b,h,hy
    klo=1
    khi=ndim
    do while(khi-klo>1)
        k=(khi+klo)/2
        if(xvec(k)>tau)then
            khi=k
        else
            klo=k
        endif
    enddo
    h=xvec(khi)-xvec(klo)
    if(almostequal(xvec(khi),xvec(klo),4))&
            stop 'bad xvec input in splint'
    a=(xvec(khi)-tau)/h
    b=(tau-xvec(klo))/h
    yval=a*yvec(klo)+b*yvec(khi)+((a**3-a)*y2vec(klo)+&
         (b**3-b)*y2vec(khi))*(h**2)/6.0_gp  
    
    !compute the derivative at point x with respect to x
    hy=yvec(khi)-yvec(klo)
    dy=hy/h+(-(3.0_gp*a**2-1.0_gp)*y2vec(klo)+&
       (3.0_gp*b**2-1.0_gp)*y2vec(khi))/6.0_gp*h
end subroutine
!=====================================================================
logical function almostequal( x, y, ulp )
    use module_base
    real(gp), intent(in) :: x
    real(gp), intent(in) :: y
    integer, intent(in) :: ulp
    almostequal = abs(x-y)<( real(ulp,gp)*&
                  spacing(max(abs(x),abs(y)))) 
end function 
!=====================================================================
end module
