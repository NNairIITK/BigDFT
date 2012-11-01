!> @file
!!  Routines to do BFGS geometry optimisation
!! @author
!!    Copyright (C) 2007-2011 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS

!> BFGS driver routine
subroutine bfgsdriver(nproc,iproc,rxyz,fxyz,epot,at,rst,in,ncount_bigdft)
    !n(c) use module_base
    use module_types
    use module_interfaces
    use minpar
    implicit none
    integer, intent(in) :: nproc,iproc
    integer, intent(inout) :: ncount_bigdft
    type(atoms_data), intent(inout) :: at
    type(input_variables), intent(inout) :: in
    type(restart_objects), intent(inout) :: rst
    real(gp), intent(inout) :: epot
    real(gp), dimension(3*at%nat), intent(inout) :: rxyz
    real(gp), dimension(3*at%nat), intent(inout) :: fxyz
    real(gp) :: fluct=0.0_gp,fnrm,fmax,fnoise
    integer :: infocode,i,ixyz,iat,istat,icall,icheck
    character(len=4) :: fn4
    character(len=40) :: comment
    logical :: move_this_coordinate
    integer ::  nr
    integer ::  nwork
    real(gp), dimension(6) :: strten
    real(gp),allocatable:: x(:),f(:),work(:)
    !character(len=4) :: fn4
    !character(len=40) :: comment
    !real(gp), dimension(3*at%nat) :: rxyz0,rxyzwrite
    !character(len=*), parameter :: subname='bfgs'

    in%inputPsiId=1
    icheck=0
    !if(iproc==0) write(*,*) 'EPOT=',epot
    !return

    nr=0
    do i=1,3*at%nat
        iat=(i-1)/3+1
        ixyz=mod(i-1,3)+1
        if(move_this_coordinate(at%ifrztyp(iat),ixyz)) nr=nr+1
    enddo
    parmin%iflag=0
    nwork=nr*nr+3*nr+3*nr*nr+3*nr
    allocate(work(nwork),stat=istat)
    if(istat/=0) stop 'ERROR: failure allocating work.'
    allocate(x(nr),stat=istat)
    if(istat/=0) stop 'ERROR: failure allocating x.'
    allocate(f(nr),stat=istat)
    if(istat/=0) stop 'ERROR: failure allocating f.'
    icall=0
    do 
        !call nebforce(n,np,x,f,fnrmtot,pnow,nproc,iproc,atoms,rst,ll_inputs,ncount_bigdft)
        !do ip=1,np-1
        !    call atomic_copymoving_forward(atoms,n,f(1,ip),nr,fa(1,ip))
        !enddo
        !if(icall/=0) then
            call call_bigdft(nproc,iproc,at,rxyz,in,epot,fxyz,strten,fnoise,rst,infocode)
            ncount_bigdft=ncount_bigdft+1
        !endif
        call atomic_copymoving_forward(at,3*at%nat,fxyz,nr,f)
        call atomic_copymoving_forward(at,3*at%nat,rxyz,nr,x)
        call fnrmandforcemax(fxyz,fnrm,fmax,at%nat)
        if(fmax<3.d-1) call updatefluctsum(fnoise,fluct) !n(m)
        call convcheck(fmax,fluct*in%frac_fluct,in%forcemax,icheck) !n(m)
        if(iproc==0) write(*,*) 'ICHECK ',icheck
        if(icheck>5) parmin%converged=.true.
        !call calmaxforcecomponentanchors(atoms,np,f(1,1),fnrm,fspmax)
        !call checkconvergence(parmin,fspmax)
        !if(ncount_bigdft>in%ncount_cluster_x-1)
        !if(iproc==0) write(*,*) 'nr=',nr,f(1)
        if (iproc == 0) then
           write(fn4,'(i4.4)') ncount_bigdft
           write(comment,'(a,1pe10.3)')'BFGS:fnrm= ',sqrt(fnrm)
           call  write_atomic_file(trim(in%dir_output)//'posout_'//fn4,epot,rxyz,at,trim(comment),forces=fxyz)
        endif
        call bfgs_reza(iproc,in%dir_output,nr,x,epot,f,nwork,work,in%betax,sqrt(fnrm),fmax, &
            ncount_bigdft,fluct*in%frac_fluct,fluct,at)
        !x(1:nr)=x(1:nr)+1.d-2*f(1:nr)
        call atomic_copymoving_backward(at,nr,x,3*at%nat,rxyz)
        if(parmin%converged) then
           if(iproc==0) write(16,'(a,i0,a)') "   BFGS converged in ",icall," iterations"
           if(iproc==0) then
              write(fn4,'(i4.4)') ncount_bigdft
              write(comment,'(a,1pe10.3)')'BFGS:fnrm= ',sqrt(fnrm)
              call  write_atomic_file(trim(in%dir_output)//'posout_'//fn4,epot,rxyz,at,trim(comment),forces=fxyz)
           endif
        endif
        !if(ncount_bigdft>in%ncount_cluster_x-1)
        !do ip=1,np-1
        !    call atomic_copymoving_backward(atoms,nr,xa(1,ip),n,x(1,ip))
        !enddo
        if(parmin%converged) exit
        if(parmin%iflag<=0) exit
        icall=icall+1
        if(icall>in%ncount_cluster_x) exit
    enddo
    deallocate(work,stat=istat);if(istat/=0) stop 'ERROR: failure deallocating work.'
    deallocate(x,stat=istat);if(istat/=0) stop 'ERROR: failure deallocating x.'
    deallocate(f,stat=istat);if(istat/=0) stop 'ERROR: failure deallocating f.'
END SUBROUTINE bfgsdriver


subroutine inithess(iproc,nr,nat,rat,atoms,hess)

    use module_types
    implicit none
    integer :: iproc,nr,nat,iat,jat,nsb,nrsqtwo,i,j,k,info
    real(kind=8) :: rat(3,nat),hess(nr,nr),r0types(4,4),fctypes(4,4),soft,hard
    type(atoms_data), intent(inout) :: atoms
    integer, allocatable::ita(:),isb(:,:)
    real(8), allocatable::r0bonds(:),fcbonds(:),evec(:,:),eval(:),wa(:)
    real(kind=8) :: dx,dy,dz,r,tt
    real(8) :: eval_i,evec_i


    nrsqtwo=2*nr**2
    if(nr/=3*atoms%nat) then
        stop 'ERROR: This subroutine works only for systems without fixed atoms.'
    endif
    allocate(ita(nat),isb(10*nat,2),r0bonds(10*nat),fcbonds(10*nat))
    allocate(evec(nr,nr),eval(nr),wa(nrsqtwo))
    do iat=1,nat
        if(trim(atoms%atomnames(atoms%iatype(iat)))=='H') then
            ita(iat)=1
        elseif(trim(atoms%atomnames(atoms%iatype(iat)))=='C') then
            ita(iat)=2
        elseif(trim(atoms%atomnames(atoms%iatype(iat)))=='N') then
            ita(iat)=3
        elseif(trim(atoms%atomnames(atoms%iatype(iat)))=='O') then
            ita(iat)=4
        else
            if(iproc==0) then
                write(*,'(a)') 'ERROR: This PBFGS is only implemented for systems which '
                write(*,'(a)') '       contain only organic elements, namely H,C,N,O.'
                write(*,'(a)') '       so use BFGS instead.'
            endif
            stop
        endif
    enddo
    call init_parameters(r0types,fctypes)
    !r0types(1:4,1:4)=2.d0 ; fctypes(1:4,1:4)=5.d2
    nsb=0
    do iat=1,nat
        do jat=iat+1,nat
            dx=rat(1,jat)-rat(1,iat)
            dy=rat(2,jat)-rat(2,iat)
            dz=rat(3,jat)-rat(3,iat)
            r=sqrt(dx**2+dy**2+dz**2)
            !if(iat==21 .and. jat==27 .and. iproc==0) then
            !    write(*,*) 'REZA ',r,1.35d0*r0types(ita(iat),ita(jat))
            !endif
            if(r<1.35d0*r0types(ita(iat),ita(jat))) then
                nsb=nsb+1
                if(nsb>10*nat) stop 'ERROR: too many stretching bonds, is everything OK?'
                isb(nsb,1)=iat
                isb(nsb,2)=jat
                r0bonds(nsb)=r0types(ita(iat),ita(jat)) !CAUTION: equil. bond length from amber
                !r0bonds(nsb)=r !CAUTION: current bond length assumed as equil. 
                fcbonds(nsb)=fctypes(ita(iat),ita(jat))
            endif
        enddo
    enddo
    if(iproc==0) write(*,*) 'NSB ',nsb
    !if(iproc==0) then
    !    do i=1,nsb
    !        write(*,'(a,i5,2f20.10,2i4,2(x,a))') 'PAR ', &
    !            i,r0bonds(i),fcbonds(i),isb(i,1),isb(i,2), &
    !            trim(atoms%atomnames(atoms%iatype(isb(i,1)))),trim(atoms%atomnames(atoms%iatype(isb(i,2))))
    !    enddo
    !endif
    call pseudohess(nat,rat,nsb,isb(1,1),isb(1,2),fcbonds,r0bonds,hess)
    evec(1:nr,1:nr)=hess(1:nr,1:nr)
    !if(iproc==0) write(*,*) 'HESS ',hess(:,:)

    call DSYEV('V','L',nr,evec,nr,eval,wa,nrsqtwo,info)


    if(info/=0) stop 'ERROR: DSYEV in inithess failed.'
    if(iproc==0) then
        do i=1,nr
            write(*,'(i5,es20.10)') i,eval(i)
        enddo
    endif
    !stop
    hard=eval(nr)
    soft=eval(nr-nsb+1)
    do k=1,nr
        if(eval(k)<soft) then
            eval(k)=soft
        endif
        eval(k)=1.d0/sqrt(eval(k)**2+soft**2)
    enddo
    do i=1,nr
    do j=i,nr
        tt=0.d0
        do k=1,nr
            !ep=1.d0/max(1.d-5,eval(k))
            !ep=sqrt(ep**2+(20.d0/eval(nr))**2)
            !if(eval(k
            !ep=sqrt(eval(k)**2+constant**2)
            tt=tt+eval(k)*evec(i,k)*evec(j,k)
        enddo
        hess(i,j)=tt
    enddo
    enddo
    do i=1,nr
    do j=1,i-1
        hess(i,j)=hess(j,i)
    enddo
    enddo
    deallocate(ita,isb,r0bonds,fcbonds,evec,eval,wa)
end subroutine inithess


subroutine init_parameters(r0,fc)
    implicit none
    integer :: i,j
    real(kind=8) :: r0(4,4),fc(4,4)
    !((0.0104 / 0.239) / 27.2114) * (0.529177^2) = 0.000447802457
    r0(1,1)=0.80d0/0.529d0
    r0(2,1)=1.09d0/0.529d0 ; r0(2,2)=1.51d0/0.529d0
    r0(3,1)=1.01d0/0.529d0 ; r0(3,2)=1.39d0/0.529d0 ; r0(3,3)=1.10d0/0.529d0
    r0(4,1)=0.96d0/0.529d0 ; r0(4,2)=1.26d0/0.529d0 ; r0(4,3)=1.10d0/0.529d0 ; r0(4,4)=1.10/0.529d0
    do i=1,4
        do j=i+1,4
            r0(i,j)=r0(j,i)
        enddo
    enddo
    fc(1,1)=1.00d3*4.48d-4
    fc(2,1)=3.40d2*4.48d-4 ; fc(2,2)=3.31d2*4.48d-4
    fc(3,1)=4.34d2*4.48d-4 ; fc(3,2)=4.13d2*4.48d-4 ; fc(3,3)=4.56d3*4.48d-4
    fc(4,1)=5.53d2*4.48d-4 ; fc(4,2)=5.43d2*4.48d-4 ; fc(4,3)=4.56d3*4.48d-4 ; fc(4,4)=4.56d3*4.48d-4
    do i=1,4
        do j=i+1,4
            fc(i,j)=fc(j,i)
        enddo
    enddo
end subroutine init_parameters


subroutine pseudohess(nat,rat,nbond,indbond1,indbond2,sprcons,xl0,hess)
    implicit none
    integer :: nat,nbond,indbond1(nbond),indbond2(nbond)
    real(kind=8) :: rat(3,nat),sprcons(nbond),xl0(nbond),hess(3*nat,3*nat)
    integer :: iat,jat,i,j,ibond
    real(kind=8) :: dx,dy,dz,r2,r,rinv! n(c) r3inv !,rinv2,rinv4,rinv8,rinv10,rinv14,rinv16
    real(kind=8) :: dxsq,dysq,dzsq,dxdy,dxdz,dydz,tt1,tt2,tt3
    real(kind=8) :: h11,h22,h33,h12,h13,h23
    do j=1,3*nat
        do i=1,3*nat
            hess(i,j)=0.d0
        enddo
    enddo
    do ibond=1,nbond
        iat=indbond1(ibond)
        jat=indbond2(ibond)
        dx=rat(1,iat)-rat(1,jat)
        dy=rat(2,iat)-rat(2,jat)
        dz=rat(3,iat)-rat(3,jat)
        r2=dx**2+dy**2+dz**2
        r=sqrt(r2) ; rinv=1.d0/r !n(c) ; r3inv=rinv**3
        !rinv2=1.d0/r2
        !rinv4=rinv2*rinv2
        !rinv8=rinv4*rinv4
        !rinv10=rinv8*rinv2
        !rinv14=rinv10*rinv4
        !rinv16=rinv8*rinv8
        dxsq=dx*dx ; dysq=dy*dy ; dzsq=dz*dz
        dxdy=dx*dy ; dxdz=dx*dz ; dydz=dy*dz
        !tt1=672.d0*rinv16
        !tt2=48.d0*rinv14
        !tt3=192.d0*rinv10
        !tt4=24.d0*rinv8
        tt1=sprcons(ibond)
        tt2=xl0(ibond)*rinv
        tt3=tt2*rinv**2
        !calculating the six distinct elements of 6 by 6 block
        !h11=dxsq*tt1-tt2-dxsq*tt3+tt4
        !h22=dysq*tt1-tt2-dysq*tt3+tt4
        !h33=dzsq*tt1-tt2-dzsq*tt3+tt4
        !h12=dxdy*tt1-dxdy*tt3
        !h13=dxdz*tt1-dxdz*tt3
        !h23=dydz*tt1-dydz*tt3

        !k_b*(1-l0/l+l0*(x_i-x_j)^2/l^3)
        h11=tt1*(1.d0-tt2+dxsq*tt3)
        h22=tt1*(1.d0-tt2+dysq*tt3)
        h33=tt1*(1.d0-tt2+dzsq*tt3)
        h12=tt1*dxdy*tt3
        h13=tt1*dxdz*tt3
        h23=tt1*dydz*tt3
        i=3*(iat-1)+1 ; j=3*(jat-1)+1
        !filling upper-left traingle (summing-up is necessary)
        hess(i+0,i+0)=hess(i+0,i+0)+h11
        hess(i+0,i+1)=hess(i+0,i+1)+h12
        hess(i+1,i+1)=hess(i+1,i+1)+h22
        hess(i+0,i+2)=hess(i+0,i+2)+h13
        hess(i+1,i+2)=hess(i+1,i+2)+h23
        hess(i+2,i+2)=hess(i+2,i+2)+h33
        !filling lower-right traingle (summing-up is necessary)
        hess(j+0,j+0)=hess(j+0,j+0)+h11
        hess(j+0,j+1)=hess(j+0,j+1)+h12
        hess(j+1,j+1)=hess(j+1,j+1)+h22
        hess(j+0,j+2)=hess(j+0,j+2)+h13
        hess(j+1,j+2)=hess(j+1,j+2)+h23
        hess(j+2,j+2)=hess(j+2,j+2)+h33
        !filling 3 by 3 block
        !summing-up is not needed but it may be necessary for PBC
        hess(i+0,j+0)=-h11 ; hess(i+1,j+0)=-h12 ; hess(i+2,j+0)=-h13
        hess(i+0,j+1)=-h12 ; hess(i+1,j+1)=-h22 ; hess(i+2,j+1)=-h23
        hess(i+0,j+2)=-h13 ; hess(i+1,j+2)=-h23 ; hess(i+2,j+2)=-h33
        !write(*,'(i3,5es20.10)') ibond,hess(i+0,i+0),tt1,tt2,tt3,xl0(ibond)
    enddo
    !filling the lower triangle of 3Nx3N Hessian matrix
    do i=1,3*nat-1
        do j =i+1,3*nat
            hess(j,i)=hess(i,j)
        enddo
    enddo
    
end subroutine pseudohess


subroutine bfgs_reza(iproc,dir_output,nr,x,epot,f,nwork,work,alphax,fnrm,fmax,ncount_bigdft,flt1,flt2,atoms)
    use minpar, only:parmin
    use module_base
    use module_types
    implicit none
    integer :: iproc,nr,nwork,mf,my,ms,nrsqtwo,iw1,iw2,iw3,iw4,info,i,j,l,mx
    integer :: ncount_bigdft
    character(len=*), intent(in) :: dir_output
    real(kind=8) :: x(nr),f(nr),epot,work(nwork),alphax,flt1,flt2
    type(atoms_data), intent(inout) :: atoms
    !real(8), allocatable::eval(:),umat(:)
    !type(parameterminimization)::parmin
    real(kind=8) :: DDOT,tt1,tt2,de,fnrm,fmax,beta
    real(kind=8) :: tt3,tt4,tt5,tt6
    real(8), save::epotold,alpha,alphamax,zeta
    logical, save::reset
    integer, save::isatur
    if(nwork/=nr*nr+3*nr+3*nr*nr+3*nr) then
        stop 'ERROR: size of work array is insufficient.'
    endif
    nrsqtwo=nr*nr*2
    mf=nr*nr+1       !for force of previous iteration in wiki notation
    my=mf+nr         !for y_k in wiki notation
    ms=my+nr         !for s_k in wiki notation
    iw1=ms+nr        !work array to keep the hessian untouched
    iw2=iw1+nr*nr    !for work array of DSYTRF
    iw3=iw2+nrsqtwo  !for p_k in wiki notation
    mx =iw3+nr       !for position of previous iteration
    iw4=mx+nr        !for eigenvalues of inverse og hessian
    if(parmin%iflag==0) then
        parmin%iflag=1
        parmin%converged=.false.   !! STEFAN Stefan stefan
        parmin%iter=0
        epotold=epot
        alpha=8.d-1
        reset=.false.
        alphamax=0.9d0
        zeta=1.d0
        isatur=0
        if(iproc==0) then
        open(unit=1390,file=trim(dir_output)//'bfgs_eigenvalues.dat',status='replace')
        close(1390)
        endif
    else
        parmin%iter=parmin%iter+1
    endif
    if(fnrm<min(6.d-2,max(1.d-2,2.d-3*sqrt(real(nr,8))))) then
        if(isatur<99) isatur=isatur+1
    else
        isatur=0
    endif
    de=epot-epotold
    !fnrm=calnorm(nr,f);fmax=calmaxforcecomponent(nr,f)
    if(iproc==0) then
    !write(*,'(a10,i5,es23.15,es11.3,2es12.5,2es12.4,i3)') &
    !    'GEOPT_BFGS',parmin%iter,epot,de,fnrm,fmax,zeta,alpha,isatur
    !       '(I5,1x,I5,2x,a10,2x,1pe21.14,2x,e9.2,1(1pe11.3),3(1pe10.2),2x,a,I3,2x,a,1pe8.2E1)'
    write(*,'(i5,1x,i5,2x,a10,2x,1es21.14,2x,es9.2,es11.3,3es10.2,2x,a7,i3)') &
        ncount_bigdft,parmin%iter,'GEOPT_BFGS',epot,de,fmax,fnrm,flt1,flt2,'isatur=',isatur
    write(16,'(i5,1x,i5,2x,a10,2x,1es21.14,2x,es9.2,es11.3,3es10.2,2x,a7,i3)') &
        ncount_bigdft,parmin%iter,'GEOPT_BFGS',epot,de,fmax,fnrm,flt1,flt2,'isatur=',isatur
    endif
    close(16)
    open(unit=16,file=trim(dir_output)//'geopt.mon',status='unknown',position='APPEND')
    !if(parmin%iter==602) then
    !    do i=1,nr/3
    !        write(31,*) x(i*3-2),x(i*3-1),x(i*3-0)
    !    enddo
    !    stop
    !endif
    !if(fmax<parmin%fmaxtol) then
    if(parmin%converged) then
        !parmin%converged=.true.
        parmin%iflag=0
        if(iproc==0) then
        write(*,'(a,i4,es23.15,2es12.5)') &
            'BFGS FINISHED: itfire,epot,fnrm,fmax ',parmin%iter,epot,fnrm,fmax
        endif
        return
    endif

    !if(de>0.d0 .and. zeta>1.d-1) then
    if(de>5.d-2) then
        epot=epotold
        x(1:nr)=work(mx:mx-1+nr)
        f(1:nr)=work(mf:mf-1+nr)
        reset=.true.
        !alpha=max(alpha*0.5d0/1.1d0,1.d-2)
        zeta=max(zeta*2.d-1,1.d-3)
        isatur=0
    else
        !zeta=1.d0
        !if(zeta>1.d-1) zeta=min(zeta*1.1d0,1.d0)
        zeta=min(zeta*1.1d0,1.d0)
        !isatur=isatur+1
    endif
    if(parmin%iter==0 .or. reset) then
        reset=.false.
        !if(isatur>=10) then
        !    reset=.false.
        !    !alpha=5.d-1
        !endif

        if(trim(parmin%approach)=='PBFGS') then
            call inithess(iproc,nr,atoms%nat,x,atoms,work(1))
        else
            work(1:nr*nr)=0.d0
            do i=1,nr
                work(i+(i-1)*nr)=zeta*alphax
            enddo
        endif
        work(iw3:iw3-1+nr)=zeta*alphax*f(1:nr)
    else
        work(ms:ms-1+nr)=x(1:nr)-work(mx:mx-1+nr)
        work(my:my-1+nr)=work(mf:mf-1+nr)-f(1:nr)
        tt1=DDOT(nr,work(my),1,work(ms),1)
        do i=1,nr
            tt2=0.d0
            do j=1,nr
                tt2=tt2+work(i+(j-1)*nr)*work(my-1+j)
            enddo
            work(iw2-1+i)=tt2
        enddo
        tt2=DDOT(nr,work(my),1,work(iw2),1)
        !write(21,*) parmin%iter,tt1,tt2
        !tt1=max(tt1,1.d-2)
        do i=1,nr
            do j=i,nr
                l=i+(j-1)*nr
                work(l)=work(l)+(tt1+tt2)*work(ms-1+i)*work(ms-1+j)/tt1**2- &
                    (work(iw2-1+i)*work(ms-1+j)+work(iw2-1+j)*work(ms-1+i))/tt1
                work(j+(i-1)*nr)=work(l)
            enddo
        enddo
        !do i=1,nr
        !    tt2=0.d0
        !    do j=1,nr
        !        tt2=tt2+work(j+(i-1)*nr)*f(j)
        !    enddo
        !    work(iw3-1+i)=tt2
        !enddo
        !write(31,*) zeta
        work(iw1:iw1-1+nr*nr)=work(1:nr*nr)
        call DSYEV('V','L',nr,work(iw1),nr,work(iw4),work(iw2),nrsqtwo,info)
        if(info/=0) stop 'ERROR: DSYEV in bfgs_reza failed.'
        tt1=work(iw4+0)    ; tt2=work(iw4+1)    ; tt3=work(iw4+2)
        tt4=work(iw4+nr-3) ; tt5=work(iw4+nr-2) ; tt6=work(iw4+nr-1)
        if(iproc==0) then
        open(unit=1390,file=trim(dir_output)//'bfgs_eigenvalues.dat',status='old',position='append')
        write(1390,'(i5,6es15.5)') parmin%iter,tt1,tt2,tt3,tt4,tt5,tt6
        close(1390)
        endif
        work(iw3:iw3-1+nr)=0.d0
        if(isatur<3) then
            beta=1.d-1/alphax
        elseif(isatur<6) then
            beta=1.d-2/alphax
        elseif(isatur<10) then
            beta=1.d-3/alphax
        else
            beta=1.d-3/alphax
        endif
        !do j=1,nr
        !    if(work(iw4-1+j)>0.d0) then
        !        tt3=work(iw4-1+j)
        !        exit
        !    enddo
        !enddo
        tt3=alphax*0.5d0
        do j=1,nr
            tt1=DDOT(nr,work(iw1+nr*(j-1)),1,f,1)
            if(work(iw4-1+j)<tt3) then
                tt4=tt3
            else
                tt4=work(iw4-1+j)
            endif
            tt2=1.d0/sqrt(1.d0/tt4**2+beta**2)
            do i=1,nr
                work(iw3-1+i)=work(iw3-1+i)+tt1*work(iw1-1+i+nr*(j-1))*tt2
            enddo
        enddo
    endif
    epotold=epot
    work(mf:mf-1+nr)=f(1:nr)
    work(mx:mx-1+nr)=x(1:nr)
    alpha=min(alphamax,alpha*1.1d0)
    x(1:nr)=x(1:nr)+alpha*work(iw3:iw3-1+nr)
end subroutine bfgs_reza


!> Driver for the LBFGS routine found on the Nocedal Homepage
!! The subroutines have only been modified slightly, so a VIMDIFF will show all modifications!
!! This is helpfull when we are looking for the source of problems during BFGS runs
subroutine lbfgsdriver(nproc,iproc,rxyz,fxyz,etot,at,rst,in,ncount_bigdft,fail) 
  use module_base
  use module_types
  use module_interfaces
!  use par_driver
  use minpar
  implicit none
!  type(driverparameters)::par
  integer, intent(in) :: nproc,iproc
  integer, intent(inout) :: ncount_bigdft
  type(atoms_data), intent(inout) :: at
  type(input_variables), intent(inout) :: in
  type(restart_objects), intent(inout) :: rst
  real(gp), intent(inout) :: etot
  real(gp), dimension(3*at%nat), intent(inout) :: rxyz
  logical, intent(out) :: fail
  real(gp), dimension(3*at%nat), intent(out) :: fxyz

  !n(c) real(gp), dimension(3*at%nat):: txyz, sxyz
  real(gp) :: fluct,fnrm, fnoise
  real(gp) :: fmax
!  logical :: check
  integer :: check
  integer :: infocode,i,ixyz,iat,nitsd
  real(gp) :: fnormmax_sw,etotprev
  character(len=4) :: fn4
  character(len=40) :: comment
  logical :: move_this_coordinate

  integer ::  n,nr,ndim
  integer ::  NWORK
  real(gp),allocatable:: X(:),G(:),DIAG(:),W(:)
  real(gp):: F,TEPS!,XTOL,GTOL,,STPMIN,STPMAX
  real(gp), dimension(6) :: strten
  real(gp), dimension(3*at%nat) :: rxyz0,rxyzwrite
  integer ::  IPRINT(2),IFLAG,ICALL,M
  character(len=*), parameter :: subname='bfgs'
  integer :: i_stat,i_all

  check=0

!  call init_driver(par)     !Initialize the parameters
  parmin%finstep=0
  parmin%alpha=1.d0
  fail=.false.
  fnrm=1.d10
  nitsd=10!500                 !Maximum number of SD steps before entering BFGS
  fnormmax_sw=in%forcemax!1.e-2_gp      !SD till the max force comp is less than this value
  
  !Dummy variables
  !n(c) txyz=0._gp
  !n(c) sxyz=0._gp
  
  if (iproc==0)    write(*,*) 'Maximum number of SD steps used in the beginning: ',nitsd

  call steepdes(nproc,iproc,at,rxyz,etot,fxyz,rst,ncount_bigdft,fnrm,fnoise,in,&
       fnormmax_sw,nitsd,fluct)
  etotprev=etot
  rxyz0=rxyz     !Save initial positions, since the unconstrained degrees of freedom will be updated upon them
  rxyzwrite=rxyz
  call fnrmandforcemax(fxyz,fnrm,fmax,at%nat)
  !call fnrmandforcemax(fxyz,fnrm,fmax,at)
  !check if the convergence is reached after SD
  call convcheck(fmax,fluct*in%frac_fluct,in%forcemax,check) !n(m)

  if (check.gt.5) then
     if (iproc.eq.0) write(*,*) 'Converged before entering BFGS'
     return
  endif


  !Make a list of all degrees of freedom that should be passed to bfgs
  n=3*at%nat
  nr=0
  do i=1,3*at%nat
     iat=(i-1)/3+1
     ixyz=mod(i-1,3)+1
     if(move_this_coordinate(at%ifrztyp(iat),ixyz)) nr=nr+1
  enddo
  if(iproc==0) write(*,*) 'DOF: n,nr ',n,nr
     NDIM=nr
     NWORK=NDIM*(2*parmin%MSAVE +1)+2*parmin%MSAVE
      
     allocate(X(NDIM),stat=i_stat)
     call memocc(i_stat,X,'X',subname)
     allocate(G(NDIM),stat=i_stat)
     call memocc(i_stat,G,'G',subname)
     allocate(DIAG(NDIM),stat=i_stat)
     call memocc(i_stat,DIAG,'DIAG',subname)
     allocate(W(NWORK),stat=i_stat)
     call memocc(i_stat,W,'W',subname)

     call atomic_copymoving_forward(at,n,rxyz,nr,X)

     N=nr
     M=parmin%MSAVE
     IPRINT(1)= 1
     IPRINT(2)= 0
     F=etot
!     We do not wish to provide the diagonal matrices Hk0, and 
!     therefore set DIAGCO to FALSE.

     TEPS=0.0_gp
     ICALL=0
     IFLAG=0

 20   CONTINUE
        if (parmin%IWRITE) then
           if (iproc == 0) then
              write(fn4,'(i4.4)') ncount_bigdft
              write(comment,'(a,1pe10.3)')'BFGS:fnrm= ',sqrt(fnrm)
              call  write_atomic_file(trim(in%dir_output)//'posout_'//fn4,etot,rxyz,at,trim(comment),forces=fxyz)
           endif
           parmin%IWRITE=.false.
        endif
        rxyzwrite=rxyz

        if (fmax < 3.d-1) call updatefluctsum(fnoise,fluct) !n(m)
   if (iproc==0.and.ICALL.ne.0.and.parmin%verbosity > 0) & 
              &write(16,'(I5,1x,I5,2x,a11,1x,1pe21.14,2x,e9.2,1(1pe11.3),3(1pe10.2),2x,a,I3,2x,a,1pe8.2E1)')&
              &ncount_bigdft,ICALL,"GEOPT_LBFGS",etot,etot-etotprev,fmax,sqrt(fnrm),fluct*in%frac_fluct,fluct&
              &,"BFGS-it=",parmin%finstep,"alpha=",parmin%alpha
   if (iproc==0.and.ICALL.ne.0.and.parmin%verbosity > 0) & 
              & write(* ,'(I5,1x,I5,2x,a11,1x,1pe21.14,2x,e9.2,1(1pe11.3),3(1pe10.2),2x,a,I3,2x,a,1pe8.2E1)')&
              &ncount_bigdft,ICALL,"GEOPT_LBFGS",etot,etot-etotprev,fmax,sqrt(fnrm),fluct*in%frac_fluct,fluct&
              &,"BFGS-it=",parmin%finstep,"alpha=",parmin%alpha
              etotprev=etot
              if (iproc==0.and.ICALL.ne.0.and.parmin%verbosity > 0) write(*,'(1x,a,1pe14.5,2(1x,a,1pe14.5))')&
                           'FORCES norm(Ha/Bohr): maxval=',fmax,'fnrm2=',fnrm,'fluct=', fluct
              call convcheck(fmax,fluct*in%frac_fluct, in%forcemax,check) !n(m)
              if (ncount_bigdft >= in%ncount_cluster_x) goto 50
              close(16)
              open(unit=16,file=trim(in%dir_output)//'geopt.mon',status='unknown',position='APPEND')

      if(check.gt.5) then
         if(iproc==0)  write(16,'(a,i0,a)') "   BFGS converged in ",ICALL," iterations"
         if (iproc == 0) then
            write(fn4,'(i4.4)') ncount_bigdft
            write(comment,'(a,1pe10.3)')'BFGS:fnrm= ',sqrt(fnrm)
            call  write_atomic_file(trim(in%dir_output)//'posout_'//fn4,etot,rxyz,at,trim(comment),forces=fxyz)
         endif
         goto 100
      endif

      
      rxyz=rxyz0
      call atomic_copymoving_backward(at,nr,X,n,rxyz)
!      txyz=rxyz
!      alpha=0._gp
!      call atomic_axpy(at,txyz,alpha,sxyz,rxyz)
      in%inputPsiId=1
!      if(ICALL.ne.0) call call_bigdft(nproc,iproc,at,rxyz,in,F,fxyz,rst,infocode)
      if(ICALL.ne.0) call call_bigdft(nproc,iproc,at,rxyz,in,F,fxyz,strten,fnoise,rst,infocode)
      if(ICALL.ne.0) ncount_bigdft=ncount_bigdft+1
      call atomic_copymoving_forward(at,n,fxyz,nr,G)
      etot=F
      G=-G
      call fnrmandforcemax(fxyz,fnrm,fmax,at%nat)
!      call fnrmandforcemax(fxyz,fnrm,fmax,at)

      CALL LBFGS(IPROC,IN,PARMIN,N,M,X,F,G,DIAG,IPRINT,TEPS,W,IFLAG)
      IF(IFLAG.LE.0) GO TO 50
      ICALL=ICALL + 1
!     We allow at most the given number of evaluations of F and G
      if(ncount_bigdft>in%ncount_cluster_x-1)  then
        goto 100
      endif
      close(16)
      open(unit=16,file=trim(in%dir_output)//'geopt.mon',status='unknown',position='append')
      GO TO 20
  50  CONTINUE
        if (iproc==0) write(*,*) "# Error in BFGS, switching to SD and CG"
        if (iproc==0) write(16,*) "Error in BFGS, switching to SD and CG"
        rxyz(:)=rxyzwrite(:)
        fail=.true.
 100  CONTINUE
        
      i_all=-product(shape(X))*kind(X)
      deallocate(X,stat=i_stat)
      call memocc(i_stat,i_all,'X',subname)
      i_all=-product(shape(G))*kind(G)
      deallocate(G,stat=i_stat)
      call memocc(i_stat,i_all,'G',subname)
      i_all=-product(shape(DIAG))*kind(DIAG)
      deallocate(DIAG,stat=i_stat)
      call memocc(i_stat,i_all,'DIAG',subname)
      i_all=-product(shape(W))*kind(W)
      deallocate(W,stat=i_stat)
      call memocc(i_stat,i_all,'W',subname)

END SUBROUTINE lbfgsdriver

subroutine atomic_copymoving_forward(atoms,n,x,nr,xa)
    use module_types
    implicit none
    type(atoms_data), intent(inout) :: atoms
    integer :: n,nr,i,iat,ixyz,ir
    real(kind=8) :: x(n),xa(nr)
    logical :: move_this_coordinate
    ir=0
    do i=1,3*atoms%nat
        iat=(i-1)/3+1
        ixyz=mod(i-1,3)+1
        if(move_this_coordinate(atoms%ifrztyp(iat),ixyz)) then
            ir=ir+1
            xa(ir)=x(i)
        endif
    enddo
    if(ir/=nr) stop 'ERROR: inconsistent number of relaxing DOF'
END SUBROUTINE atomic_copymoving_forward


subroutine atomic_copymoving_backward(atoms,nr,xa,n,x)
    use module_types
    implicit none
    type(atoms_data), intent(inout) :: atoms
    integer :: n,nr,i,iat,ixyz,ir
    real(kind=8) :: x(n),xa(nr)
    logical :: move_this_coordinate
    ir=0
    do i=1,3*atoms%nat
        iat=(i-1)/3+1
        ixyz=mod(i-1,3)+1
        if(move_this_coordinate(atoms%ifrztyp(iat),ixyz)) then
            ir=ir+1
            x(i)=xa(ir)
        endif
    enddo
    if(ir/=nr) stop 'ERROR: inconsistent number of relaxing DOF'
END SUBROUTINE atomic_copymoving_backward
