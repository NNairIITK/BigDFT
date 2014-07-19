!! @file
!! @author Bastian Schaefer
!! @section LICENCE
!!    Copyright (C) 2014 UNIBAS
!!    This file is not freely distributed.
!!    A licence is necessary from UNIBAS
module module_lst
implicit none

!to be removed when moved to bigdft:
integer, parameter :: gp=kind(1.0d0)

contains
subroutine lst_penalty(nat,rxyzR,rxyzP,rxyz,lambda,val,force)
!computes the linear synchronous penalty function
!and gradient
!
!see:
!
!Halgren, T. A., & Lipscomb, W. N. (1977). The synchronous-transit
!method for determining reaction pathways and locating molecular
!transition states. Chemical Physics Letters, 49(2), 225–232.
!doi:10.1016/0009-2614(77)80574-5
!
!and
!
!Behn, A., Zimmerman, P. M., Bell, A. T., & Head-Gordon, M. (2011).
!Incorporating Linear Synchronous Transit Interpolation into
!the Growing String Method: Algorithm and Applications.
!Journal of Chemical Theory and Computation, 7(12), 4019–4025.
!doi:10.1021/ct200654u
!
!and
!
!BS: For computation of gradient see also personal notes in intel
!notebook (paper version) from June 18th, 2014
!
    implicit none
    !parameters
    integer, intent(in)  :: nat
    real(gp), intent(in) :: rxyzR(3,nat) !reactant
    real(gp), intent(in) :: rxyzP(3,nat) !product
    real(gp), intent(in) :: rxyz(3,nat)  !the positon at which
                                         !the penalty fct. is to
                                         !be evaluated
    real(gp), intent(in)  :: lambda !interpolation parameter
                                    !(0<=lambda<=1)
    real(gp), intent(out) :: val  !the value of the penalty function
    real(gp), intent(out) :: force(3,nat) !the negative gradient of
                                          !the penal. fct. at
                                          !rxyz
    !internal
    integer :: b !outer loop
    integer :: a !inner loop
    real(gp) :: rabi !interpolated interatomic distances
    real(gp) :: rabim4 !rim4=ri**(-4)
    real(gp) :: rabC !computed interatomic distances
    real(gp) :: rabR !interatomic dist. of reactant
    real(gp) :: rabP !interatomic dist. of product
    real(gp) :: oml
    real(gp) :: rxRb, ryRb, rzRb
    real(gp) :: rxPb, ryPb, rzPb
    real(gp) :: rxb, ryb, rzb
    real(gp) :: rxd, ryd, rzd
    real(gp) :: rabiMrabC
    real(gp) :: tt, ttx, tty, ttz
    real(gp) :: waxi, wayi, wazi 

    oml=1.0_gp-lambda
    val=0.0_gp
    !first sum
    do b = 1, nat-1
        rxRb = rxyzR(1,b)
        ryRb = rxyzR(2,b)
        rzRb = rxyzR(3,b)
        rxPb = rxyzP(1,b)
        ryPb = rxyzP(2,b)
        rzPb = rxyzP(3,b)
        rxb = rxyz(1,b)
        ryb = rxyz(2,b)
        rzb = rxyz(3,b)
        do a = b+1, nat
            !compute interatomic distances of reactant
            rabR = (rxyzR(1,a)-rxRb)**2+&
                   (rxyzR(2,a)-ryRb)**2+&
                   (rxyzR(3,a)-rzRb)**2
            rabR = sqrt(rabR)
            !compute interatomic distances of product
            rabP = (rxyzP(1,a)-rxPb)**2+&
                   (rxyzP(2,a)-ryPb)**2+&
                   (rxyzP(3,a)-rzPb)**2
            rabP = sqrt(rabP)
            !compute interpolated interatomic distances
            rabi = oml*rabR+lambda*rabP
            rabim4 = 1.0_gp / (rabi**4)
            !compute interatomic distances at rxyz
            rxd = rxyz(1,a)-rxb
            ryd = rxyz(2,a)-ryb
            rzd = rxyz(3,a)-rzb
            rabC = rxd**2+ryd**2+rzd**2
            rabC = sqrt(rabC)
            !compute function value
            rabiMrabC = rabi - rabC
            val = val + rabiMrabC**2 * rabim4 
            !compute negative gradient
            rabiMrabC = -2.0_gp*rabiMrabC * rabim4 / rabC
            ttx = rabiMrabC * rxd
            tty = rabiMrabC * ryd
            ttz = rabiMrabC * rzd
            force(1,b) = force(1,b) + ttx
            force(2,b) = force(2,b) + tty
            force(3,b) = force(3,b) + ttz
            force(1,a) = force(1,a) - ttx
            force(2,a) = force(2,a) - tty
            force(3,a) = force(3,a) - ttz
        enddo
    enddo

    !second sum
    do a = 1, nat
        waxi = oml*rxyzR(1,a) + lambda*rxyzP(1,a)
        wayi = oml*rxyzR(2,a) + lambda*rxyzP(2,a)
        wazi = oml*rxyzR(3,a) + lambda*rxyzP(3,a)
        ttx = waxi-rxyz(1,a)
        tty = wayi-rxyz(2,a)
        ttz = wazi-rxyz(3,a)
        tt   =  ttx**2 + tty**2 + ttz**2
        val = val + 1.e-6_gp * tt
        force(1,a) = force(1,a) + 2.e-6_gp * ttx
        force(2,a) = force(2,a) + 2.e-6_gp * tty
        force(3,a) = force(3,a) + 2.e-6_gp * ttz
    enddo
end subroutine

!subroutine fire_original(nat,rxyz,fxyz,epot,fnrm,count_fr,displ)
!    implicit none
!    !parameters
!    integer, intent(in) :: nat
!    integer :: maxit
!    real(gp)  :: fmax_tol,fnrm_tol,dt_max
!    real(gp), intent(inout) :: rxyz(3,nat),fxyz(3,nat)
!    real(gp)  :: epot,fmax,fnrm
!    logical :: success
!    real(gp), intent(inout) :: count_fr
!    real(gp), intent(out) :: displ
!    !internal
!    real(gp) :: vxyz(3,nat),ff(3,nat),power
!!    real(gp) :: ekin
!    real(gp) :: at1,at2,at3,vxyz_norm,fxyz_norm,alpha
!    integer :: iat,iter,check,cut
!    integer, parameter :: n_min=5
!    real(gp), parameter :: f_inc=1.1_gp, f_dec=0.5_gp, alpha_start=0.1_gp,f_alpha=0.99_gp
!    real(gp) :: dt
!    real(gp) :: ddot,dnrm2
!    character(len=5), allocatable :: xat(:)
!    real(gp) :: rxyzIn(3,nat)
!
!    real(gp) :: epot_prev, anoise
!    real(gp) :: rxyz_prev(3,nat)
!    real(gp) :: counter_reset
!    real(gp) :: disp(3,nat),dmy1,dmy2
!
!    displ=0.0_gp
!    dt_max = 4.5d-2
!    maxit=50000
!
!    rxyzIn=rxyz
!
!    success=.false.
!
!!    dt=0.1_gp*dt_max
!    dt=dt_max
!    alpha=alpha_start
!    vxyz=0.0_gp
!    check=0
!    cut=1
!    call energyandforces(nat,rxyz,ff,epot)
!        count_fr=count_fr+1.0_gp
!!write(*,'(a)')'FIRE         iter     epot            fmax           fnrm
!!dt             alpha'
!    do iter=1,maxit
!        epot_prev=epot
!        rxyz_prev=rxyz
!        1000 continue
!        call daxpy(3*nat,dt,vxyz(1,1),1,rxyz(1,1),1)
!        call daxpy(3*nat,0.5_gp*dt*dt,ff(1,1),1,rxyz(1,1),1)
!!        ekin=ddot(3*nat,vxyz(1,1),1,vxyz(1,1),1)
!!        ekin=0.5_gp*ekin
!        call energyandforces(nat,rxyz,fxyz,epot)
!        count_fr=count_fr+1.0_gp
!        disp=rxyz-rxyz_prev
!        call fmaxfnrm(nat,disp,dmy1,dmy2)
!        displ=displ+dmy2
!        do iat=1,nat
!           at1=fxyz(1,iat)
!           at2=fxyz(2,iat)
!           at3=fxyz(3,iat)
!           !C Evolution of the velocities of the system
!           vxyz(1,iat)=vxyz(1,iat) + (.5_gp*dt) * (at1 + ff(1,iat))
!           vxyz(2,iat)=vxyz(2,iat) + (.5_gp*dt) * (at2 + ff(2,iat))
!           vxyz(3,iat)=vxyz(3,iat) + (.5_gp*dt) * (at3 + ff(3,iat))
!           !C Memorization of old forces
!           ff(1,iat) = at1
!           ff(2,iat) = at2
!           ff(3,iat) = at3
!        end do
!!write(*,*)rxyz(1,1),vxyz(1,1)
!        call fmaxfnrm(nat,fxyz,fmax,fnrm)
!!        call convcheck(fmax,fmax_tol,check)
!!        if(check > 5)then
!         if (fnrm.lt.1.e-3_gp)   then
!!            write(*,*)'FIRE converged, force calls: ',iter+1,epot
!            write(100,'(a,x,i0,5(1x,es14.7))')'FIRE converged # e evals, epot, fmax, fnrm, dt, alpha: ',int(count_fr),epot,fmax,fnrm,dt,alpha
!!            write(*,'(a,x,i0,x,i0,5(1x,es14.7))')'FIRE converged iter, epot,
!!            fmax, fnrm, dt, alpha:
!!            ',iter,iter+int(count_fr)+1+int(counter_reset),epot,fmax,fnrm,dt,alpha
!            success=.true.
!            return
!        endif
!        power = ddot(3*nat,fxyz,1,vxyz,1)
!        vxyz_norm = dnrm2(3*nat,vxyz,1)
!        fxyz_norm = dnrm2(3*nat,fxyz,1)
!        vxyz = (1.0_gp-alpha)*vxyz + alpha * fxyz * vxyz_norm / fxyz_norm
!        if(power<=0)then
!!write(*,*)'freeze'
!            vxyz=0.0_gp
!            cut=iter
!            dt=dt*f_dec
!            alpha=alpha_start
!        else if(power > 0 .and. iter-cut >n_min)then
!            dt= min(dt*f_inc,dt_max)
!            alpha = alpha*f_alpha
!        endif
!    enddo
!    if(fmax > fmax_tol)then
!            write(100,'(a,x,i0,5(1x,es14.7))')'FIRE ERROR not converged iter, epot, fmax, fnrm, dt, alpha: ',iter,epot,fmax,fnrm,dt,alpha
!    endif
!end subroutine


end module
!!program penalty_unittest
!!use module_lst
!!implicit none
!!      include 'sizes.i'
!!      include 'atoms.i'
!!      include 'syntrn.i'
!!integer :: istat,nat,iat
!!real(gp), allocatable :: rxyz1(:,:),rxyz2(:,:),g1(:,:),g2(:,:),rxyz(:,:),rand(:,:)
!!character(len=100) :: filename,line, filename2
!!real(gp) :: lambda,val,rr
!!real(gp) :: transit
!!
!!
!!call get_command_argument(1,filename)
!!call get_command_argument(2,filename2)
!!call get_command_argument(3,line)
!!read(line,*)lambda
!!
!!nat=0
!!open(unit=33,file=trim(adjustl(filename)))
!!read(33,*)
!!read(33,*)
!!read(33,*)
!!read(33,*)
!!do
!!    read(33,'(a)',iostat=istat)line
!!    if(istat/=0)exit
!!    nat=nat+1
!!enddo
!!write(*,*)'Found nat: ',nat
!!rewind(33)
!!read(33,*)
!!read(33,*)
!!read(33,*)
!!read(33,*)
!!allocate(rxyz1(3,nat),rxyz2(3,nat),g1(3,nat),g2(3,nat),rxyz(3,nat),rand(3,nat))
!!do iat=1,nat
!!    read(33,'(a)',iostat=istat)line
!!    if(istat/=0)exit
!!    read(line,*)rxyz1(1,iat),rxyz1(2,iat),rxyz1(3,iat)
!!enddo
!!close(33)
!!open(unit=33,file=trim(adjustl(filename2)))
!!read(33,*)
!!read(33,*)
!!read(33,*)
!!read(33,*)
!!do iat=1,nat
!!    read(33,'(a)',iostat=istat)line
!!    if(istat/=0)exit
!!    read(line,*)rxyz2(1,iat),rxyz2(2,iat),rxyz2(3,iat)
!!enddo
!!close(33)
!!
!!
!!n=nat
!!t=lambda
!!pm=0.0_gp
!!
!!call random_seed
!!call random_number(rand)
!!call random_number(rr)
!!
!!do iat=1,nat
!!xmin1(3*iat-2) = rxyz1(1,iat)
!!xmin1(3*iat-1) = rxyz1(2,iat)
!!xmin1(3*iat) = rxyz1(3,iat)
!!xmin2(3*iat-2) = rxyz2(1,iat)
!!xmin2(3*iat-1) = rxyz2(2,iat)
!!xmin2(3*iat) = rxyz2(3,iat)
!!enddo
!!
!!!rr=lambda
!!rxyz=(1._gp-rr)*rxyz1+rr*rxyz2 + rand
!!call lst_penalty(nat,rxyz1,rxyz2,rxyz,lambda,val,g1)
!!g1=-g1
!!!write(*,*)'mine:   ',val
!!val=transit(rxyz,g2)
!!!write(*,*)'tinker: ',val
!!do iat=1,nat
!!write(*,*)g1(1,iat)-g2(1,iat)
!!write(*,*)g1(2,iat)-g2(2,iat)
!!write(*,*)g1(3,iat)-g2(3,iat)
!!enddo
!!
!!end program
!!
!!      function transit (xx,g)
!!      implicit none
!!      include 'sizes.i'
!!      include 'atoms.i'
!!      include 'syntrn.i'
!!      integer i,j,nvar
!!      integer ix,iy,iz
!!      integer jx,jy,jz
!!      real*8 transit,value
!!      real*8 xci,yci,zci
!!      real*8 xcd,ycd,zcd
!!      real*8 x1i,y1i,z1i
!!      real*8 x1d,y1d,z1d
!!      real*8 x2i,y2i,z2i
!!      real*8 x2d,y2d,z2d
!!      real*8 xmi,ymi,zmi
!!      real*8 xmd,ymd,zmd
!!      real*8 gamma,term
!!      real*8 termx,termy,termz
!!      real*8 cutoff,cutoff2
!!      real*8 r1,r2,rc,rm
!!      real*8 ri,ri4,rd
!!      real*8 wi,wc,wd
!!      real*8 tq,pq
!!      real*8 xx(*)
!!      real*8 g(*)
!!      character*9 mode
!!!
!!!
!!!     zero out the synchronous transit function and gradient
!!!
!!      value = 0.0d0
!!      nvar = 3 * n
!!      do i = 1, nvar
!!         g(i) = 0.0d0
!!      end do
!!      tq = 1.0d0 - t
!!!
!!!     set the cutoff distance for interatomic distances
!!!
!!      cutoff = 1000.0d0
!!      cutoff2 = cutoff**2
!!!
!!!     set the type of synchronous transit path to be used
!!!
!!!      if (pm .eq. 0.0d0) then
!!         mode = 'LINEAR'
!!!      else
!!!         mode = 'QUADRATIC'
!!!         pq = 1.0d0 - pm
!!!      end if
!!!
!!!     portion based on interpolated interatomic distances
!!!
!!      do i = 1, n-1
!!         iz = 3 * i
!!         iy = iz - 1
!!         ix = iz - 2
!!         xci = xx(ix)
!!         yci = xx(iy)
!!         zci = xx(iz)
!!         x1i = xmin1(ix)
!!         y1i = xmin1(iy)
!!         z1i = xmin1(iz)
!!         x2i = xmin2(ix)
!!         y2i = xmin2(iy)
!!         z2i = xmin2(iz)
!!!         if (mode .eq. 'QUADRATIC') then
!!!            xmi = xm(ix)
!!!            ymi = xm(iy)
!!!            zmi = xm(iz)
!!!         end if
!!         do j = i+1, n
!!            jz = 3 * j
!!            jy = jz - 1
!!            jx = jz - 2
!!            xcd = xci - xx(jx)
!!            ycd = yci - xx(jy)
!!            zcd = zci - xx(jz)
!!            x1d = x1i - xmin1(jx)
!!            y1d = y1i - xmin1(jy)
!!            z1d = z1i - xmin1(jz)
!!            x2d = x2i - xmin2(jx)
!!            y2d = y2i - xmin2(jy)
!!            z2d = z2i - xmin2(jz)
!!            rc = xcd**2 + ycd**2 + zcd**2
!!            r1 = x1d**2 + y1d**2 + z1d**2
!!            r2 = x2d**2 + y2d**2 + z2d**2
!!!            if (min(rc,r1,r2) .lt. cutoff2) then
!!               rc = sqrt(rc)
!!               r1 = sqrt(r1)
!!               r2 = sqrt(r2)
!!               ri = tq*r1 + t*r2
!!!               if (mode .eq. 'QUADRATIC') then
!!!                  xmd = xmi - xm(jx)
!!!                  ymd = ymi - xm(jy)
!!!                  zmd = zmi - xm(jz)
!!!                  rm = sqrt(xmd**2+ymd**2+zmd**2)
!!!                  gamma = (rm-pq*r1-pm*r2) / (pm*pq)
!!!                  ri = ri + gamma*t*tq
!!!               end if
!!               ri4 = ri**4
!!               rd = rc - ri
!!               value = value + rd**2/ri4
!!               term = 2.0d0 * rd/(ri4*rc)
!!               termx = term * xcd
!!               termy = term * ycd
!!               termz = term * zcd
!!               g(ix) = g(ix) + termx
!!               g(iy) = g(iy) + termy
!!               g(iz) = g(iz) + termz
!!               g(jx) = g(jx) - termx
!!               g(jy) = g(jy) - termy
!!               g(jz) = g(jz) - termz
!!!            end if
!!         end do
!!      end do
!!!
!!!     portion used to supress rigid rotations and translations
!!!
!!      do i = 1, nvar
!!         wc = xx(i)
!!         wi = tq*xmin1(i) + t*xmin2(i)
!!         wd = wc - wi
!!         value = value + 0.000001d0*wd**2
!!         g(i) = g(i) + 0.000002d0*wd
!!      end do
!!      transit = value
!!      return
!!      end
