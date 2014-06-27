module module_saddle

contains

subroutine findbonds(nat,rcov,pos,nbond,iconnect)
!has to be called before findsad (if operating in biomolecule mode)
use module_base
    implicit none
    !parameters
    integer, intent(in) :: nat
    real(gp), intent(in) :: rcov(nat)
    real(gp), intent(in) :: pos(3,nat)
    integer, intent(out) :: nbond
    integer, intent(out) :: iconnect(2,1000)
    !internal
    integer :: iat,jat
    real(gp) :: dist2
    nbond=0
    do iat=1,nat
        do jat=1,iat-1
            dist2=(pos(1,iat)-pos(1,jat))**2+(pos(2,iat)-pos(2,jat))**2+(pos(3,iat)-pos(3,jat))**2
            if (dist2.le.(1.2d0*(rcov(iat)+rcov(jat)))**2) then
                nbond=nbond+1
                if (nbond.gt.1000) stop 'nbond>1000'
                iconnect(1,nbond)=iat
                iconnect(2,nbond)=jat
            endif
        enddo
    enddo
    write(*,*) 'FOUND ',nbond,' BOND'
end subroutine

subroutine findsad(imode,nat,alat,rcov,alpha0_trans,alpha0_rot,curvgraddiff,nit_trans,&
           nit_rot,nhistx_trans,nhistx_rot,tolc,tolf,tightenfac,rmsdispl0,&
           trustr,wpos,etot,fout,minmode,fnrmtol,count,count_sd,displ,ec,&
           converged,atomnames,nbond,iconnect,alpha_stretch0,recompIfCurvPos,maxcurvrise,cutoffratio)
use module_base
use module_global_variables, only: runObj
!imode=1 for clusters
!imode=2 for biomolecules
        implicit none
        integer, intent(in) :: imode
        integer, intent(in) :: recompIfCurvPos
        real(gp), intent(in) :: alat(3),rcov(nat)
        integer :: nat,nbond
        integer,optional :: iconnect(2,nbond)
        real(gp), dimension(3,nat) :: wpos,fout
        real(gp) :: count
        character(20):: atomnames(nat)
        real(gp), allocatable, dimension(:,:,:) :: rxyz,fxyz,ff,rr,rrr,fff,fstretch,fxyzraw,rxyzraw
        real(gp), allocatable, dimension(:,:) :: aa,dd,dds,delta
        real(gp), allocatable, dimension(:) ::eval,work,res,scpr,wold
        real(gp) :: t1,t2,t3
        real(gp) :: t1raw,t2raw,t3raw
        logical check,steep
        character(2) atomname
        real(gp) :: displ
        real(gp) :: maxd,scl
        real(gp) :: alpha_stretch, fnrm,etotp,etotold,ec,alpha,tt,dt,count_sd,detot,s,fnrmtol,cosangle,etot,st
        integer :: lwork, iat , l, itswitch, nhist, ndim, it, ihist, i
real(gp),optional :: alpha_stretch0
real(gp) :: alpha0_trans,alpha0_rot,curvgraddiff,tolc,tolf,tightenfac,rmsdispl0,trustr
integer:: nit_trans,nit_rot,nhistx_trans,nhistx_rot,recompute
real(gp) :: minmode(3,nat),tol,displold,rmsdispl,curv,displ2,overlap
logical :: converged,optCurvConv
logical :: flag,tooFar,fixfragmented,subspaceSucc
real(gp), allocatable, dimension(:,:) :: gradrot,minmodeold
real(gp) :: tnatdmy(3,nat)
integer :: nputback,lastit
real(gp), intent(in) :: maxcurvrise
real(gp), intent(in) :: cutoffratio

      
integer :: j,jdim,idim,info


real(gp) :: minmode0(3,nat)
character(len=100) :: filename
    !functions
    real(gp) :: ddot,dnrm2

minmode0=minmode

    rmsdispl=rmsdispl0

    flag=.true.
    fixfragmented=.false.
    nputback=0
    converged=.false.
    subspaceSucc=.true.
    lastit=-(10*nhistx_trans)

    check=.true.

    displ=0.d0
    displ2=0.d0
    displold=0.d0
    curv=1000.d0
    !allocate arrays
    lwork=1000+10*nat**2

    alpha_stretch=alpha_stretch0


! allocate arrays
    lwork=1000+10*nat**2
    allocate(rxyz(3,nat,0:nhistx_trans),fxyz(3,nat,0:nhistx_trans),aa(nhistx_trans,nhistx_trans),fxyzraw(3,nat,0:nhistx_trans),rxyzraw(3,nat,0:nhistx_trans))
    allocate(eval(nhistx_trans),res(nhistx_trans),work(lwork))
    allocate(ff(3,nat,0:nhistx_trans),rr(3,nat,0:nhistx_trans),dd(3,nat),dds(3,nat),delta(3,nat))
    allocate(fff(3,nat,0:nhistx_trans),rrr(3,nat,0:nhistx_trans),scpr(nhistx_trans),wold(nbond),fstretch(3,nat,0:nhistx_trans))
    allocate(gradrot(3,nat),minmodeold(3,nat))
    wold=0.d0
    fstretch=0.d0

    do iat=1,nat
        do l=1,3
            rxyz(l,iat,0)=wpos(l,iat)
        enddo
    enddo

    call fixfrag_posvel(nat,rcov,rxyz(1,1,0),tnatdmy,1,fixfragmented)

    count=count+1.d0
!        call energyandforces(nat,rxyz(1,1,0),fxyz(1,1,0),etot)
!        energyandforces(nat,alat,rxyz(1,1,0),fxyz(1,1,0),etot,'cnt_enf_geopt')
    call minenergyandforces(imode,nat,alat,rxyz(1,1,0),rxyzraw(1,1,0),&
    fxyz(1,1,0),fstretch(1,1,0),fxyzraw(1,1,0),etot,iconnect,nbond,atomnames,&
    wold,alpha_stretch0,alpha_stretch)
    ec=ec+1.d0
    if(imode==2)rxyz(:,:,0)=rxyz(:,:,0)+alpha_stretch*fstretch(:,:,0)
    t1=0.d0 ; t2=0.d0 ; t3=0.d0
    t1raw=0.d0 ; t2raw=0.d0 ; t3raw=0.d0
    do iat=1,nat
        t1raw=t1raw+fxyzraw(1,iat,0)**2 ; t2raw=t2raw+fxyzraw(2,iat,0)**2;t3raw=t3raw+fxyzraw(3,iat,0)**2
    enddo
    fnrm=sqrt(t1raw+t2raw+t3raw)
    etotold=etot
    etotp=etot

    itswitch=2
    ndim=0
    nhist=0
    alpha=alpha0_trans
    do it=1,nit_trans
        if (check) write(100,*) 'it:',it,etot,fnrm,itswitch
        nhist=nhist+1

        if ((.not. subspaceSucc) .or. fnrm.gt.100.d0 .or. it.le.itswitch) then
          ndim=0
          steep=.true.
!          if (it.gt.itswitch) itswitch=it+nhistx_trans
          if (check) write(100,*) "STEEP"
      else
          steep=.false.
!          alpha=alpha0_trans
      endif

! make space in the history list
        if (nhist.gt.nhistx_trans) then
            nhist=nhistx_trans
            do ihist=0,nhist-1
                do iat=1,nat
                    do l=1,3
                        rxyz(l,iat,ihist)=rxyz(l,iat,ihist+1)
                        fxyz(l,iat,ihist)=fxyz(l,iat,ihist+1)
                        rxyzraw(l,iat,ihist)=rxyzraw(l,iat,ihist+1)
                        fxyzraw(l,iat,ihist)=fxyzraw(l,iat,ihist+1)
                        fstretch(l,iat,ihist)=fstretch(l,iat,ihist+1)
                     enddo
                enddo
            enddo
        endif

       if(fnrm > 10.d0*fnrmtol)then
            flag=.true.
       endif

!       !START FINDING LOWEST MODE
!       if(fnrm<=tightenfac*fnrmtol .and. flag .and. curv<0.d0)then
!!       if(fnrm<=tightenfac*fnrmtol .and. flag)then
!           if(check) write(100,*)'tighten'
!           write(*,*)'tighten'
!           tol=tolf
!           recompute=it
!           flag=.false.
!       else if(it==1)then
!           tol=tolc
!       else
!           tol=tolc
!       endif
!       !if(abs(displ-displold)>0.029d0*sqrt(dble(3*nat))&
!       tooFar = abs(displ-displold)>rmsdispl*sqrt(dble(3*nat))
!       if(tooFar&
!!       &.or. it==1&
!       &.or. it==1 .or. (curv>=0.d0 .and. mod(it,recompIfCurvPos)==0)&
!       &.or.recompute==it)then
!           call opt_curv(imode,nat,alat,alpha0_rot,curvgraddiff,nit_rot,nhistx_rot,rxyzraw(1,1,nhist-1),fxyzraw(1,1,nhist-1),&
!                        &minmode(1,1),curv,gradrot(1,1),&
!                        &tol,count,count_sd,displ2,ec,check,optCurvConv,iconnect,nbond,atomnames,2.d-4,maxcurvrise,cutoffratio)
!           minmode = minmode / dnrm2(3*nat,minmode(1,1),1)
!           if(.not.optCurvConv)then
!               write(*,*) 'WARNING: opt_curv failed'
!               converged=.false.
!stop 'opt_curv failed'
!               return
!           endif
!           overlap=ddot(3*nat,minmodeold(1,1),1,minmode(1,1),1)
!           write(*,*)'overlap',overlap
!           minmodeold=minmode
!           displold=displ
!           recompute=huge(1)
!           if(tooFar)then
!               flag = .true.
!           endif
!       endif
!       !END FINDING LOWEST MODE

       600 continue
       call modify_gradient(check,nat,ndim,rrr(1,1,1),eval(1),res(1),fxyz(1,1,nhist-1),alpha,dd(1,1))

       !invert gradient in minmode direction
!       dd=dd-2.d0*ddot(3*nat,dd(1,1),1,minmode(1,1),1)*minmode

       tt=0.d0
       dt=0.d0
       maxd=-huge(1.d0)
       do iat=1,nat
           dt=dd(1,iat)**2+dd(2,iat)**2+dd(3,iat)**2
           tt=tt+dt
           maxd=max(maxd,dt)
       enddo
       tt=sqrt(tt)
       maxd=sqrt(maxd)

       !trust radius approach
!       if(maxd>trustr .and. alpha>1.d-1*alpha0_trans)then
       if(maxd>trustr)then
           write(*,*)'step too large',maxd,it
           scl=0.5d0*trustr/maxd
           dd=dd*scl
           tt=tt*scl
        endif
        !do the move
        rxyz(:,:,nhist)=rxyz(:,:,nhist-1)-dd(:,:)
       call fixfrag_posvel(nat,rcov,rxyz(1,1,nhist),tnatdmy,1,fixfragmented)
!       displ=displ+tt

       if (steep) then 
           count_sd=count_sd+1.d0
       else
           count=count+1.d0
       endif

       runObj%inputs%inputPsiId=1 
       call minenergyandforces(imode,nat,alat,rxyz(1,1,nhist),rxyzraw(1,1,nhist),&
            fxyz(1,1,nhist),fstretch(1,1,nhist),fxyzraw(1,1,nhist),etotp&
            ,iconnect,nbond,atomnames,wold,alpha_stretch0,alpha_stretch)
       ec=ec+1.d0
       detot=etotp-etotold
       s=0.d0 ; st=0.d0
       t1=0.d0 ; t2=0.d0 ; t3=0.d0
       t1raw=0.d0 ; t2raw=0.d0 ; t3raw=0.d0
       do iat=1,nat
           t1raw=t1raw+fxyzraw(1,iat,nhist)**2 ; t2raw=t2raw+fxyzraw(2,iat,nhist)**2 ;t3raw=t3raw+fxyzraw(3,iat,nhist)**2
           t1=t1+fxyz(1,iat,nhist)**2 ; t2=t2+fxyz(2,iat,nhist)**2 ;t3=t3+fxyz(3,iat,nhist)**2
           s=s+dd(1,iat)**2+dd(2,iat)**2+dd(3,iat)**2
           st=st+fxyz(1,iat,nhist)*dd(1,iat)+fxyz(2,iat,nhist)*dd(2,iat)+fxyz(3,iat,nhist)*dd(3,iat)
       enddo
       fnrm=sqrt(t1raw+t2raw+t3raw)

       write(*,'(a,i6,2(1x,e21.14),1x,2(1x,es10.3),xi6,2(xes10.3))')&
       &'GEOPT it,etot,etotold,Detot,fnrm,ndim,alpha,alpha_stretch',&
       &it-1,etotp,etotold,detot,fnrm,ndim,alpha,alpha_stretch

       etot=etotp
       etotold=etot
!       if (fnrm.le.fnrmtol .and. curv<0.d0) goto 1000
       if (fnrm.le.fnrmtol) goto 1000

       !now do step in hard directions
       if(imode==2)then
!           fstretch(:,:,nhist)=fstretch(:,:,nhist)-2.d0*ddot(3*nat,fstretch(1,1,nhist),1,minmode(1,1),1)*minmode
           dds=alpha_stretch*(fstretch(:,:,nhist)-2.d0*ddot(3*nat,fstretch(1,1,nhist),1,minmode(1,1),1)*minmode)
           dt=0.d0
           maxd=-huge(1.d0)
           do iat=1,nat
!               dt=fstretch(1,iat,nhist)**2+fstretch(2,iat,nhist)**2+fstretch(3,iat,nhist)**2
               dt=dds(1,iat)**2+dds(2,iat)**2+dds(3,iat)**2
               maxd=max(maxd,dt)
           enddo
           maxd=sqrt(maxd)

           !trust radius approach
           if(maxd>trustr)then
               write(*,'(a,es10.3,xi0,xes10.3)')'hard direction step too large:maxd,it,alpha_stretch',maxd,it,alpha_stretch
               scl=0.5d0*trustr/maxd
               dds=dds*scl
!control
     dt=0.d0
     maxd=-huge(1.d0)
     do iat=1,nat
         dt=dds(1,iat)**2+dds(2,iat)**2+dds(3,iat)**2
         maxd=max(maxd,dt)
     enddo
     maxd=sqrt(maxd)
write(*,*)'control',maxd
           endif

!           rxyz(:,:,nhist)=rxyz(:,:,nhist)+alpha_stretch*fstretch(:,:,nhist)
           rxyz(:,:,nhist)=rxyz(:,:,nhist)+dds
           call fixfrag_posvel(nat,rcov,rxyz(1,1,nhist),tnatdmy,1,fixfragmented)
       endif
!write(437,*)sqrt((rxyz(1,1,nhist)-rxyz(1,2,nhist))**2+(rxyz(2,1,nhist)-rxyz(2,2,nhist))**2+(rxyz(3,1,nhist)-rxyz(3,2,nhist))**2)


       cosangle=-st/sqrt((t1+t2+t3)*s)
       if (cosangle.gt..20d0) then
           alpha=alpha*1.10d0
       else
           alpha=max(alpha*.85d0,alpha0_trans)
       endif
       if (check) write(100,*) 'cosangle ',cosangle,alpha

       call getSubSpaceEvecEval(check,nat,nhist,nhistx_trans,ndim,cutoffratio,lwork,work,rxyz,&
                              &fxyz,aa,rr,ff,rrr,fff,eval,res,subspaceSucc)

       delta=rxyz(:,:,nhist)-rxyz(:,:,nhist-1)
       !write(673,*)dnrm2(3*nat,delta(1,1),1),tt,tt-dnrm2(3*nat,delta(1,1),1)
       displ=displ+dnrm2(3*nat,delta(1,1),1)
  enddo
    write(557,*) nat
    write(557,*) "  "
    do iat=1,nat
        write(557,'(a,3(1x,e24.17))') " LJ ",wpos(1,iat),wpos(2,iat),wpos(3,iat)
    enddo
    write(557,*) "  "
    do iat=1,nat
        write(557,'(3(1x,e24.17))') minmode0(1,iat),minmode0(2,iat),minmode0(3,iat)
    enddo

    write(*,*) "No convergence in findsad"
stop 'no convergence in findsad'
    return

1000 continue
     converged=.true.
     etot=etotp
     write(*,*) "convergence at",it

     do iat=1,nat
     do l=1,3
     wpos(l,iat)= rxyz(l,iat,nhist)
     fout(l,iat)= fxyzraw(l,iat,nhist)
     enddo
     enddo

             write(556,*)
             write(556,*)10.d0,0.d0,10.d0
             write(556,*)0.d0,0.d0,10.d0
             do iat=1,nat
               write(556,*) wpos(1,iat),wpos(2,iat),wpos(3,iat),atomnames(iat)
             enddo
!             write(556,*)nat
!             write(556,*)
!             do iat=1,nat
!               write(556,*)atomnames(iat), wpos(1,iat),wpos(2,iat),wpos(3,iat)
!             enddo

  end subroutine

subroutine modify_gradient(check,nat,ndim,rrr,eval,res,fxyz,alpha,dd)
use module_base
    implicit none
    !parameters
    logical, intent(in) :: check
    integer, intent(in) :: nat
    integer, intent(in) :: ndim
    real(gp), intent(out) :: dd(3,nat)
    real(gp), intent(in) :: fxyz(3,nat)
    real(gp), intent(in) :: rrr(3,nat,ndim)
    real(gp), intent(in) :: eval(ndim)
    real(gp), intent(in) :: res(ndim)
    real(gp), intent(in) :: alpha
    !internal
    integer :: iat,jat,i,l
    real(gp) :: scpr(ndim)
    real(gp) :: tt

    ! decompose gradient
    
    do iat=1,nat
        do l=1,3
            dd(l,iat)=-fxyz(l,iat)
        enddo
    enddo
    do i=1,ndim
        scpr(i)=0.d0
        do iat=1,nat
            do l=1,3
                scpr(i)=scpr(i)-fxyz(l,iat)*rrr(l,iat,i)
            enddo
        enddo
        do iat=1,nat
            do l=1,3
                dd(l,iat)=dd(l,iat)-scpr(i)*rrr(l,iat,i)
            enddo
        enddo
    enddo
    
    if (check) write(100,*) 'alpha=',alpha
    !simple sd in space orthogonal to relevant subspace
    do iat=1,nat
        do l=1,3
            dd(l,iat)=dd(l,iat)*alpha
        enddo
    enddo
    
    do i=1,ndim
    !quasi newton in relevant subspace
        tt=scpr(i)/sqrt(eval(i)**2+res(i)**2)
        if (check) write(100,'(a,i3,3(1x,es10.3))')'i,tt,eval,res',i,&
                   1.d0/sqrt(eval(i)**2+res(i)**2),eval(i),res(i)
        do iat=1,nat
            do l=1,3
                dd(l,iat)=dd(l,iat)+tt*rrr(l,iat,i)
            enddo
        enddo
    enddo

end subroutine


subroutine getSubSpaceEvecEval(check,nat,nhist,nhistx,ndim,cutoffratio,lwork,work,rxyz,&
                              &fxyz,aa,rr,ff,rrr,fff,eval,res,success)
use module_base
    !hard-coded parameters:
    !threshold for linear dependency:
    !if (eval(idim)/eval(nhist).gt.1.d-4) then
    implicit none
    !parameters
    logical :: check
    integer, intent(in) :: nat,nhist,nhistx,lwork
    integer, intent(out) :: ndim
    real(gp), intent(in) :: rxyz(3,nat,0:nhistx),fxyz(3,nat,0:nhistx)
    real(gp), intent(out) :: aa(nhistx,nhistx),eval(nhistx),work(lwork)
    real(gp), intent(out) :: rr(3,nat,0:nhistx), ff(3,nat,0:nhistx)
    real(gp), intent(out) :: rrr(3,nat,0:nhistx), fff(3,nat,0:nhistx)
    real(gp), intent(out) :: res(nhistx)
    real(gp), intent(in) :: cutoffratio
    logical, intent(out) :: success
    !internal
    integer :: i,j,l,iat,info,idim,jdim
    real(gp) :: tt
    real(gp) :: rnorm(nhistx)
!write(444,*)'in getSub'
!write(444,*)rxyz(:,:,0)
!write(444,*)
!write(444,*)rxyz(:,:,1)

    success=.false.

!write(444,*)'nhist',nhist
    ! calculate norms
    do i=1,nhist
        rnorm(i)=0.d0
         do iat=1,nat
             do l=1,3
!                rnorm(i)=rnorm(i) + drxyz(l,iat,i)**2
                rnorm(i)=rnorm(i) + (rxyz(l,iat,i)-rxyz(l,iat,i-1))**2
!write(444,*)(rxyz(l,iat,i)-rxyz(l,iat,i-1))
             enddo
         enddo
         rnorm(i)=1.d0/sqrt(rnorm(i))
!rnorm(i) = (0.5d0*(1.d0+(((dble(i)-0.2d0*dble(nhist))/(dble(nhist)*0.25d0))/sqrt(1.d0+((dble(i)-0.2d0*dble(nhist))/(dble(nhist)*0.25d0))**2))))/sqrt(rnorm(i))
!rnorm(i)=1.d0
    enddo

    do i=1,nhist
        do j=1,nhist
            aa(i,j)=0.d0
            do iat=1,nat
                do l=1,3
                !aa(i,j)=aa(i,j) + drxyz(l,iat,i)*drxyz(l,iat,j)
                aa(i,j)=aa(i,j) + (rxyz(l,iat,i)-rxyz(l,iat,i-1))&
                &*(rxyz(l,iat,j)-rxyz(l,iat,j-1))
                enddo
            enddo
            aa(i,j)=aa(i,j)*rnorm(i)*rnorm(j)
         enddo
    enddo

    call dsyev('V',"L",nhist,aa,nhistx,eval,work,lwork,info)
    if (info.ne.0) then
        if (check) write(100,*) ' Over ', info
        return
!        stop 'info'
    endif
    do i=1,nhist
        if (check) write(100,*) "Overl ",i,eval(i)
    enddo
    
    do idim=1,nhist
        do iat=1,nat
            do l=1,3
                rr(l,iat,idim)=0.d0
                ff(l,iat,idim)=0.d0
            enddo
        enddo
    enddo

    ndim=0
    do idim=1,nhist
        if (eval(idim)/eval(nhist).gt.cutoffratio) then    ! HERE
            ndim=ndim+1
    
            do jdim=1,nhist
                do iat=1,nat
                    do l=1,3
!                        rr(l,iat,ndim)=rr(l,iat,ndim)+&
!                                      &aa(jdim,idim)*rnorm(jdim)*drxyz(l,iat,jdim)
!                        ff(l,iat,ndim)=ff(l,iat,ndim)-&
!                                      &aa(jdim,idim)*rnorm(jdim)*dfxyz(l,iat,jdim)
                         rr(l,iat,ndim)=rr(l,iat,ndim)+&
                                       &aa(jdim,idim)*rnorm(jdim)*(rxyz(l,iat,jdim)-rxyz(l,iat,jdim-1))
                         ff(l,iat,ndim)=ff(l,iat,ndim)-&
                                       &aa(jdim,idim)*rnorm(jdim)*(fxyz(l,iat,jdim)-fxyz(l,iat,jdim-1))

                    enddo
                enddo
            enddo
    
            do iat=1,nat
                do l=1,3
                    rr(l,iat,ndim)=rr(l,iat,ndim)/sqrt(abs(eval(idim)))
                    ff(l,iat,ndim)=ff(l,iat,ndim)/sqrt(abs(eval(idim)))
                enddo
            enddo
        endif
    enddo
    if (check) write(100,'(a,i3)') "ndim= ",ndim


    ! Hessian matrix in significant orthogonal subspace
    do i=1,ndim
        do j=1,ndim
            aa(i,j)=0.d0
            do iat=1,nat
                do l=1,3
                    aa(i,j)=aa(i,j) + .5d0*(rr(l,iat,i)*ff(l,iat,j)+&
                                      &rr(l,iat,j)*ff(l,iat,i))
                enddo
            enddo
        enddo
    enddo

    call dsyev('V',"L",ndim,aa,nhistx,eval,work,lwork,info)
    if (info.ne.0) then
        if(check) write(100,*) ' A ', info
        return
!        stop 'info'
    endif

    ! calculate eigenvectors
    do i=1,ndim
        do iat=1,nat
            do l=1,3
                rrr(l,iat,i)=0.d0
                fff(l,iat,i)=0.d0
            enddo
        enddo
    enddo

    do i=1,ndim
        tt=0.d0
        do j=1,ndim
            do iat=1,nat
                do l=1,3
                    rrr(l,iat,i)=rrr(l,iat,i) + aa(j,i)*rr(l,iat,j)
                    fff(l,iat,i)=fff(l,iat,i) + aa(j,i)*ff(l,iat,j)
                enddo
            enddo
        enddo
        do iat=1,nat
            do l=1,3
                tt=tt+(fff(l,iat,i)-eval(i)*rrr(l,iat,i))**2
            enddo
        enddo
        res(i)=sqrt(tt)
        if (check) write(100,'(a,i3,e14.7,1x,e12.5,2(1x,e9.2))') 'EVAL,RES ',i,eval(i),res(i)
    enddo
    success=.true.
end subroutine
             

subroutine elim_moment_fs(nat,vxyz)
use module_base
    implicit real(gp) (a-h,o-z)
    dimension vxyz(3,nat)
    sx=0.d0 ; sy=0.d0 ; sz=0.d0
    do iat=1,nat
        sx=sx+vxyz(1,iat)
        sy=sy+vxyz(2,iat)
        sz=sz+vxyz(3,iat)
    enddo
    sx=sx/nat ; sy=sy/nat ; sz=sz/nat
    do iat=1,nat
        vxyz(1,iat)=vxyz(1,iat)-sx
        vxyz(2,iat)=vxyz(2,iat)-sy
        vxyz(3,iat)=vxyz(3,iat)-sz
    enddo
end subroutine

subroutine elim_torque_fs(nat,rxyz,vxyz)
use module_base
    implicit real(gp) (a-h,o-z)
    dimension rxyz(3,nat),vxyz(3,nat),t(3)

    ! center of mass
    cmx=0.d0 ; cmy=0.d0 ; cmz=0.d0
    do iat=1,nat
        cmx=cmx+rxyz(1,iat)
        cmy=cmy+rxyz(2,iat)
        cmz=cmz+rxyz(3,iat)
    enddo
    cmx=cmx/nat ; cmy=cmy/nat ; cmz=cmz/nat
    
    do it=1,100

        ! torque and radii in planes
        t(1)=0.d0 ; t(2)=0.d0 ; t(3)=0.d0
        sx=0.d0 ; sy=0.d0 ; sz=0.d0
        do iat=1,nat
            t(1)=t(1)+(rxyz(2,iat)-cmy)*vxyz(3,iat)-&
                &(rxyz(3,iat)-cmz)*vxyz(2,iat)
            t(2)=t(2)+(rxyz(3,iat)-cmz)*vxyz(1,iat)-&
                &(rxyz(1,iat)-cmx)*vxyz(3,iat)
            t(3)=t(3)+(rxyz(1,iat)-cmx)*vxyz(2,iat)-&
                &(rxyz(2,iat)-cmy)*vxyz(1,iat)
            sx=sx+(rxyz(1,iat)-cmx)**2
            sy=sy+(rxyz(2,iat)-cmy)**2
            sz=sz+(rxyz(3,iat)-cmz)**2
        enddo

        if (t(1)**2+t(2)**2+t(3)**2.lt.1.d-22) return

        ii=0
        tmax=0.d0
        do i=1,3
            if (t(i)**2.gt.tmax**2) then
                ii=i
                tmax=t(i)
            endif
        enddo

        ! modify velocities
        if (ii.eq.1) then
            cx=t(1)/(sz+sy)
            do iat=1,nat
                vxyz(2,iat)=vxyz(2,iat)+cx*(rxyz(3,iat)-cmz)
                vxyz(3,iat)=vxyz(3,iat)-cx*(rxyz(2,iat)-cmy)
            enddo
        else if(ii.eq.2) then
            cy=t(2)/(sz+sx)
            do iat=1,nat
                vxyz(1,iat)=vxyz(1,iat)-cy*(rxyz(3,iat)-cmz)
                vxyz(3,iat)=vxyz(3,iat)+cy*(rxyz(1,iat)-cmx)
            enddo
        else if(ii.eq.3) then
            cz=t(3)/(sy+sx)
            do iat=1,nat
                vxyz(1,iat)=vxyz(1,iat)+cz*(rxyz(2,iat)-cmy)
                vxyz(2,iat)=vxyz(2,iat)-cz*(rxyz(1,iat)-cmx)
            enddo
        else
            stop 'wrong ii'
        endif

    enddo
    write(100,'(a,3(e11.3))') 'WARNING REMAINING TORQUE',t
end subroutine

subroutine opt_curv(imode,nat,alat,alpha0,curvgraddiff,nit,nhistx,rxyz_fix,fxyz_fix,dxyzin,curv,fout,fnrmtol&
                   &,count,count_sd,displ,ec,check,converged,iconnect,nbond,atomnames,alpha_stretch0,maxcurvrise,cutoffratio)!,mode)
use module_base
    implicit none
integer, intent(in) :: imode
integer, intent(in) :: nbond
integer, intent(in) :: iconnect(2,nbond)
character(len=20) :: atomnames(nat)
real(gp),intent(in) :: alpha_stretch0
real(gp) :: alpha_stretch
real(gp), intent(in) :: alat(3)
real(gp), intent(in) :: maxcurvrise,cutoffratio
    integer, intent(in) :: nit,nhistx
    real(gp), intent(in) :: alpha0,curvgraddiff
    real(gp), dimension(3,nat) :: dxyzin,fout,rxyz_fix,fxyz_fix!,mode
    real(gp), allocatable, dimension(:,:,:) ::rxyz,fxyz,ff,rr,rrr,fff,fxyzraw,rxyzraw,fstretch
    real(gp), allocatable, dimension(:,:) :: aa,dd
    real(gp), allocatable, dimension(:) :: eval,work,res,scpr
    logical, intent(out) :: converged
    logical check,steep
    integer :: lwork,nhist,i,iat,l,itswitch,ndim
    integer :: ihist,it,nat
    real(gp) :: displ,ec,curv,count
    real(gp) :: count_sd,fnrmtol,curvold,t1,t2,t3,fnrm,curvp,t1raw,t2raw,t3raw
    real(gp) :: alpha,dcurv,s,st,tt,cosangle
    logical :: subspaceSucc
    real(gp),allocatable :: wold(:)
    !functions
    real(gp) :: ddot

    alpha_stretch=alpha_stretch0

    converged =.false.
    subspaceSucc=.true.
    displ=0.d0
    ! allocate arrays
    lwork=1000+10*nat**2
    allocate(rxyz(3,nat,0:nhistx),fxyz(3,nat,0:nhistx),fxyzraw(3,nat,0:nhistx),rxyzraw(3,nat,0:nhistx),aa(nhistx,nhistx))
    allocate(eval(nhistx),res(nhistx),work(lwork))
    allocate(ff(3,nat,0:nhistx),rr(3,nat,0:nhistx),dd(3,nat))
    allocate(fff(3,nat,0:nhistx),rrr(3,nat,0:nhistx),scpr(nhistx),fstretch(3,nat,0:nhistx),wold(nbond))
    wold=0.d0

    do iat=1,nat
       do l=1,3
          rxyz(l,iat,0)=dxyzin(l,iat)
       enddo
    enddo

    count=count+1.d0
!    call minenergyandforces(nat,rxyz(1,1,0),rxyzraw(1,1,0),fxyz(1,1,0),fstretch(1,1,0),fxyzraw(1,1,0),etot,iconnect,nbond,atomnames,wold,alpha_stretch)
!    rxyz(:,:,0)=rxyz(:,:,0)+alpha_stretch*fstretch(:,:,0)
!    ec=ec+1.d0
     call mincurvgrad(imode,nat,alat,curvgraddiff,rxyz_fix(1,1),fxyz_fix(1,1),rxyz(1,1,0),&
          rxyzraw(1,1,0),fxyz(1,1,0),fstretch(1,1,0),fxyzraw(1,1,0),curv,1,ec,&
          iconnect,nbond,atomnames,wold,alpha_stretch0,alpha_stretch)
    if(imode==2)rxyz(:,:,0)=rxyz(:,:,0)+alpha_stretch*fstretch(:,:,0)
    t1=0.d0 ; t2=0.d0 ; t3=0.d0
    t1raw=0.d0 ; t2raw=0.d0 ; t3raw=0.d0
    do iat=1,nat
       t1raw=t1raw+fxyzraw(1,iat,0)**2 ; t2raw=t2raw+fxyzraw(2,iat,0)**2;t3raw=t3raw+fxyzraw(3,iat,0)**2
    enddo
    fnrm=sqrt(t1raw+t2raw+t3raw)
    curvold=curv
    curvp=curv
 
    itswitch=2
    ndim=0
    nhist=0
    alpha=alpha0
    do it=1,nit
        if (check) write(100,*) 'it:',it,curv,fnrm,itswitch
        nhist=nhist+1
 
!        if (fnrm.gt.50.d0 .or. it.le.itswitch ) then
        if ((.not. subspaceSucc) .or. fnrm.gt.2000.d0 .or. it.le.itswitch ) then
            ndim=0
            steep=.true.
!            if (it.gt.itswitch) itswitch=it+nhistx
            if (check) write(100,*) "STEEP"
        else
            steep=.false.
        endif

        ! make space in the history list
        if (nhist.gt.nhistx) then
         nhist=nhistx
         do ihist=0,nhist-1
            do iat=1,nat
               do l=1,3
                  rxyz(l,iat,ihist)=rxyz(l,iat,ihist+1)
                  fxyz(l,iat,ihist)=fxyz(l,iat,ihist+1)
                  rxyzraw(l,iat,ihist)=rxyzraw(l,iat,ihist+1)
                  fxyzraw(l,iat,ihist)=fxyzraw(l,iat,ihist+1)
                  fstretch(l,iat,ihist)=fstretch(l,iat,ihist+1)
                enddo
             enddo
          enddo
      endif


        500 continue
        call modify_gradient(check,nat,ndim,rrr(1,1,1),eval(1),res(1),fxyz(1,1,nhist-1),alpha,dd(1,1))


        tt=0.d0
        do iat=1,nat
            do l=1,3
                tt=tt+dd(l,iat)**2
            enddo
        enddo
        displ=displ+sqrt(tt)
!write(444,*)'dd',nhist

        do iat=1,nat
            rxyz(1,iat,nhist)=rxyz(1,iat,nhist-1)-dd(1,iat)
            rxyz(2,iat,nhist)=rxyz(2,iat,nhist-1)-dd(2,iat)
            rxyz(3,iat,nhist)=rxyz(3,iat,nhist-1)-dd(3,iat)
        enddo
!write(444,*)1,rxyz(1,1,nhist)

        if (steep) then 
            count_sd=count_sd+1.d0
        else
            count=count+1.d0
        endif
     call mincurvgrad(imode,nat,alat,curvgraddiff,rxyz_fix(1,1),fxyz_fix(1,1),rxyz(1,1,nhist),&
          rxyzraw(1,1,nhist),fxyz(1,1,nhist),fstretch(1,1,nhist),fxyzraw(1,1,nhist),&
          curvp,1,ec,iconnect,nbond,atomnames,wold,alpha_stretch0,alpha_stretch)
        dcurv=curvp-curvold


        s=0.d0 ; st=0.d0
        t1=0.d0 ; t2=0.d0 ; t3=0.d0
        t1raw=0.d0 ; t2raw=0.d0 ; t3raw=0.d0
        do iat=1,nat
            t1raw=t1raw+fxyzraw(1,iat,nhist)**2 ; t2raw=t2raw+fxyzraw(2,iat,nhist)**2 ;t3raw=t3raw+fxyzraw(3,iat,nhist)**2
            t1=t1+fxyz(1,iat,nhist)**2 ; t2=t2+fxyz(2,iat,nhist)**2 ;t3=t3+fxyz(3,iat,nhist)**2
            s=s+dd(1,iat)**2+dd(2,iat)**2+dd(3,iat)**2
            st=st+fxyz(1,iat,nhist)*dd(1,iat)+fxyz(2,iat,nhist)*dd(2,iat)+fxyz(3,iat,nhist)*dd(3,iat)
         enddo
         fnrm=sqrt(t1raw+t2raw+t3raw)

         write(*,'(a,i6,2(1x,es21.14),1x,4(1x,es10.3),xi6)')&
        &'CURV  it,curv,curvold,Dcurv,fnrm,alpha,alpha_stretch,ndim',&
        &it-1,curvp,curvold,dcurv,fnrm,alpha,alpha_stretch,ndim

        if (dcurv.gt.maxcurvrise .and. alpha>1.d-1*alpha0) then 
            if (check) write(100,'(a,i4,1x,e9.2)') "WARN: it,dcurv", it,dcurv
            !alpha=min(.5d0*alpha,alpha0)
            alpha=.5d0*alpha
            if (check) write(100,'(a,1x,e9.2)') 'alpha reset ',alpha
            ndim=0
!            if(.not.steep)then
                do iat=1,nat
                    rxyz(1,iat,0)=rxyzraw(1,iat,nhist-1)
                    rxyz(2,iat,0)=rxyzraw(2,iat,nhist-1)
                    rxyz(3,iat,0)=rxyzraw(3,iat,nhist-1)
                    rxyzraw(1,iat,0)=rxyzraw(1,iat,nhist-1)
                    rxyzraw(2,iat,0)=rxyzraw(2,iat,nhist-1)
                    rxyzraw(3,iat,0)=rxyzraw(3,iat,nhist-1)

                    fxyz(1,iat,0)=fxyzraw(1,iat,nhist-1)
                    fxyz(2,iat,0)=fxyzraw(2,iat,nhist-1)
                    fxyz(3,iat,0)=fxyzraw(3,iat,nhist-1)
                    fxyzraw(1,iat,0)=fxyzraw(1,iat,nhist-1)
                    fxyzraw(2,iat,0)=fxyzraw(2,iat,nhist-1)
                    fxyzraw(3,iat,0)=fxyzraw(3,iat,nhist-1)
                enddo
                nhist=1
!            endif
            goto  500
        endif
            curv=curvp
            curvold=curv
        if (fnrm.le.fnrmtol) goto 1000
        if(imode==2)then
            rxyz(:,:,nhist)=rxyz(:,:,nhist)+alpha_stretch*fstretch(:,:,nhist)
        endif

        cosangle=-st/sqrt((t1+t2+t3)*s)
write(*,*)'cosangle',cosangle
        if (cosangle.gt..20d0) then
            alpha=alpha*1.10d0
        else
            alpha=max(alpha*.85d0,0.5d0*alpha0)
        endif
        if (check) write(100,*) 'cosangle ',cosangle,alpha


       call getSubSpaceEvecEval(check,nat,nhist,nhistx,ndim,cutoffratio,lwork,work,rxyz,&
                              &fxyz,aa,rr,ff,rrr,fff,eval,res,subspaceSucc)
    enddo

    write(555,*) nat
    write(555,*) "  "
    do iat=1,nat
        write(555,'(a,3(1x,e24.17))') trim(adjustl(atomnames(iat))),dxyzin(1,iat),dxyzin(2,iat),dxyzin(3,iat)
    enddo
    write(100,*) it,curv,fnrm
    write(*,*) "No convergence in optcurv"
stop 'no convergence in optcurv'
    return
    1000 continue
    converged=.true.
    do iat=1,nat
        do l=1,3
            dxyzin(l,iat)= rxyz(l,iat,nhist)
            fout(l,iat)= fxyzraw(l,iat,nhist)
        enddo
    enddo
    curv=curvp
end subroutine
!subroutine opt_curv(nat,alpha0,curvgraddiff,nit,nhistx,rxyz_fix,fxyz_fix,dxyzin,curv,fout,fnrmtol&
!                   &,count,count_sd,displ,ec,check,converged)!,mode)
!    implicit none
!    integer, intent(in) :: nit,nhistx
!    real(gp), intent(in) :: alpha0,curvgraddiff
!    real(gp), dimension(3,nat) :: dxyzin,fout,rxyz_fix,fxyz_fix!,mode
!    real(gp), allocatable, dimension(:,:,:) :: rxyz,fxyz,ff,rr,rrr,fff
!    real(gp), allocatable, dimension(:,:) :: aa,dd
!    real(gp), allocatable, dimension(:) :: eval,work,res,scpr
!    logical, intent(out) :: converged
!    logical check,steep
!    integer :: lwork,nhist,i,iat,l,itswitch,ndim
!    integer :: ihist,it,nat
!    real(gp) :: displ,ec,curv,count
!    real(gp) :: count_sd,fnrmtol,curvold,t1,t2,t3,fnrm,curvp
!    real(gp) :: alpha,ts,dcurv,s,st,tt,cosangle
!    logical :: subspaceSucc
!    !functions
!    real(gp) :: ddot
!    
!   
!    converged =.false. 
!    subspaceSucc=.true.
!    displ=0.d0
!    ! allocate arrays
!    lwork=1000+10*nat**2
!    allocate(rxyz(3,nat,0:nhistx),fxyz(3,nat,0:nhistx),aa(nhistx,nhistx))
!    allocate(eval(nhistx),res(nhistx),work(lwork))
!    allocate(ff(3,nat,0:nhistx),rr(3,nat,0:nhistx),dd(3,nat))
!    allocate(fff(3,nat,0:nhistx),rrr(3,nat,0:nhistx),scpr(nhistx))
!    
!    do iat=1,nat
!        do l=1,3
!            rxyz(l,iat,0)=dxyzin(l,iat)
!        enddo
!    enddo
!
!    count=count+1.d0
!    call curvgrad(nat,curvgraddiff,rxyz_fix(1,1),fxyz_fix(1,1),rxyz(1,1,0)&
!                 &,curv,fxyz(1,1,0),1,ec)
!    t1=0.d0 ; t2=0.d0 ; t3=0.d0
!    do iat=1,nat
!        t1=t1+fxyz(1,iat,0)**2 ; t2=t2+fxyz(2,iat,0)**2 ; t3=t3+fxyz(3,iat,0)**2
!    enddo
!    fnrm=sqrt(t1+t2+t3)
!    curvold=curv
!    curvp=curv
!
!    itswitch=0
!    ndim=0
!    nhist=0
!    alpha=alpha0
!    do it=1,nit
!        if (check) write(100,*) 'it:',it,curv,fnrm,itswitch
!        nhist=nhist+1
!    
!        !if (fnrm.gt.100.d0 .or. it.le.itswitch ) then
!        if ((.not. subspaceSucc) .or. fnrm.gt.1000.d0 .or. it.le.itswitch ) then
!            ndim=0
!            steep=.true.
!            if (it.gt.itswitch) itswitch=it+nhistx
!            if (check) write(100,*) "STEEP"
!        else
!            steep=.false.
!        endif
!    
!      !make space in the history list
!      if (nhist.gt.nhistx) then
!          nhist=nhistx
!          do ihist=0,nhist-1
!             do iat=1,nat
!                do l=1,3
!                  rxyz(l,iat,ihist)=rxyz(l,iat,ihist+1)
!                  fxyz(l,iat,ihist)=fxyz(l,iat,ihist+1)
!                enddo
!             enddo
!          enddo
!      endif
!
!        ! decompose gradient
!        500 continue
!        do iat=1,nat
!            do l=1,3
!                dd(l,iat)=-fxyz(l,iat,nhist-1)
!            enddo
!        enddo
!        do i=1,ndim
!            scpr(i)=0.d0
!            do iat=1,nat
!                do l=1,3
!                    scpr(i)=scpr(i)-fxyz(l,iat,nhist-1)*rrr(l,iat,i)
!                enddo
!            enddo
!            do iat=1,nat
!                do l=1,3
!                    dd(l,iat)=dd(l,iat)-scpr(i)*rrr(l,iat,i)
!                enddo
!            enddo
!        enddo
!
!        ts=0.d0
!        do iat=1,nat
!            do l=1,3
!                ts=ts+dd(l,iat)**2
!            enddo
!        enddo
!
!        if (check) write(100,*) 'alpha=',alpha
!        do iat=1,nat
!            do l=1,3
!                dd(l,iat)=dd(l,iat)*alpha
!            enddo
!        enddo
!
!        do i=1,ndim
!            tt=scpr(i)/sqrt(eval(i)**2+res(i)**2)
!            if (check) write(100,'(a,i3,3(1x,e10.3))') 'i,tt,eval,res ',&
!                       &i,1.d0/sqrt(eval(i)**2+res(i)**2),eval(i),res(i)
!            do iat=1,nat
!                do l=1,3
!                    dd(l,iat)=dd(l,iat)+tt*rrr(l,iat,i)
!                enddo
!            enddo
!        enddo
!
!        tt=0.d0
!        do iat=1,nat
!            do l=1,3
!                tt=tt+dd(l,iat)**2
!            enddo
!        enddo
!        displ=displ+sqrt(tt)
!
!        do iat=1,nat
!           rxyz(1,iat,nhist)=rxyz(1,iat,nhist-1)-dd(1,iat)
!           rxyz(2,iat,nhist)=rxyz(2,iat,nhist-1)-dd(2,iat)
!           rxyz(3,iat,nhist)=rxyz(3,iat,nhist-1)-dd(3,iat)
!        enddo
!    
!        if (steep) then
!            count_sd=count_sd+1.d0
!        else
!            count=count+1.d0
!        endif
!        call curvgrad(nat,curvgraddiff,rxyz_fix(1,1),fxyz_fix(1,1),rxyz(1,1,nhist),&
!                     &curvp,fxyz(1,1,nhist),1,ec)
!        dcurv=curvp-curvold
!
!         write(*,'(a,i6,2(1x,e21.14),1x,5(1x,e10.3),xi6)')&
!        &'CURV  it,curv,curvold,Dcurv,fnrm,fnrmp/fnrm,dnrm/fnrm,ndim',&
!        &it-1,curvp,curvold,dcurv,fnrm,sqrt(ts)/fnrm,sqrt(tt)/fnrm,alpha,ndim
!
!        s=0.d0 ; st=0.d0
!        t1=0.d0 ; t2=0.d0 ; t3=0.d0
!        do iat=1,nat
!            t1=t1+fxyz(1,iat,nhist)**2 ; t2=t2+fxyz(2,iat,nhist)**2 ; t3=t3+&
!                                                              &fxyz(3,iat,nhist)**2
!            s=s+dd(1,iat)**2+dd(2,iat)**2+dd(3,iat)**2
!            st=st+fxyz(1,iat,nhist)*dd(1,iat)+fxyz(2,iat,nhist)*dd(2,iat)&
!               &+fxyz(3,iat,nhist)*dd(3,iat)
!        enddo
!        fnrm=sqrt(t1+t2+t3)
!
!
!        if (dcurv.gt.1.d-9 .and. alpha>1.d-1*alpha0) then
!            if (check) write(100,'(a,i4,1x,e9.2)') "WARN: it,dcurv", it,dcurv
!            !alpha=min(.5d0*alpha,alpha0)
!            alpha=.5d0*alpha
!            if (check) write(100,'(a,1x,e9.2)') 'alpha reset ',alpha
!            ndim=0
!            if(.not.steep)then
!                do iat=1,nat
!                    rxyz(1,iat,0)=rxyz(1,iat,nhist-1)
!                    rxyz(2,iat,0)=rxyz(2,iat,nhist-1)
!                    rxyz(3,iat,0)=rxyz(3,iat,nhist-1)
!
!                    fxyz(1,iat,0)=fxyz(1,iat,nhist-1)
!                    fxyz(2,iat,0)=fxyz(2,iat,nhist-1)
!                    fxyz(3,iat,0)=fxyz(3,iat,nhist-1)
!                enddo
!                nhist=1
!            endif
!            goto  500
!        else
!            curv=curvp
!            curvold=curv
!        endif
!        if (fnrm.le.fnrmtol) goto 1000
!!if(abs(ddot(3*nat,rxyz1(1,1),1,mode(1,1),1))>.99d0)then
!!goto 1000
!!endif
!
!
!        cosangle=-st/sqrt((t1+t2+t3)*s)
!        if (cosangle.gt..20d0) then
!            alpha=alpha*1.10d0
!        else
!            alpha=max(alpha*.85d0,alpha0)
!        endif
!        if (check) write(100,*) 'cosangle ',cosangle,alpha
!
!
!       call getSubSpaceEvecEval(check,nat,nhist,nhistx,ndim,lwork,work,&
!            &rxyz,fxyz,aa,rr,ff,rrr,fff,eval,res,subspaceSucc)
!    enddo
!
!    write(555,*) nat
!    write(555,*) "  "
!    do iat=1,nat
!        write(555,'(a,3(1x,e24.17))') " LJ ",dxyzin(1,iat),dxyzin(2,iat),dxyzin(3,iat)
!    enddo
!    write(100,*) it,curv,fnrm
!    write(100,*) "No convergence in optcurv"
!    return
!    1000 continue
!    converged = .true.
!    do iat=1,nat
!        do l=1,3
!            dxyzin(l,iat)= rxyz(l,iat,nhist)
!            fout(l,iat)= fxyz(l,iat,nhist)
!        enddo
!    enddo
!    curv=curvp
!end subroutine


subroutine curvgrad(nat,alat,diff,rxyz1,fxyz1,vec,curv,rotforce,imethod,ec)
    !computes the (curvature along vec) = vec^t H vec / (vec^t*vec)
    !vec mus be normalized
use module_base
use module_energyandforces
    implicit none
    !parameters
    integer, intent(in) :: nat,imethod
    real(gp), intent(in) :: rxyz1(3,nat),fxyz1(3,nat),diff
    real(gp), intent(inout) :: vec(3,nat)
    real(gp), intent(out) :: curv
    real(gp), intent(out) :: rotforce(3,nat)
    real(gp), intent(inout) :: ec
    real(gp), intent(in) :: alat(3)
    !internal
    real(gp) :: diffinv, etot2
    real(gp),allocatable :: rxyz2(:,:), fxyz2(:,:)
    real(gp),allocatable :: drxyz(:,:), dfxyz(:,:)
    !functions
    real(gp), external :: dnrm2,ddot
real(gp) :: fluct

    
!    diff=1.d-3 !lennard jones

    diffinv=1.0d0/(diff)

    allocate(rxyz2(3,nat),fxyz2(3,nat))
    allocate(drxyz(3,nat),dfxyz(3,nat))
    call elim_moment_fs(nat,vec(1,1))
!    call elim_torque_fs(nat,rxyz1(1,1),vec(1,1))
    call elim_torque_reza(nat,rxyz1(1,1),vec(1,1))

    vec = vec / dnrm2(3*nat,vec(1,1),1)
    rxyz2 = rxyz1 + diff * vec

    call energyandforces(nat,alat,rxyz2(1,1),fxyz2(1,1),etot2)
!    call energyandforces(nat, alat, rxyz2(1,1),fxyz2(1,1),etot2,'cnt_enf_forcebar_decomp')
    ec=ec+1.d0

    if(imethod==1)then
        dfxyz = fxyz2(:,:)-fxyz1(:,:)
        drxyz = rxyz2(:,:)-rxyz1(:,:)
        drxyz = drxyz * diffinv
        curv  = - diffinv * ddot(3*nat,dfxyz(1,1),1,drxyz(1,1),1)

        rotforce = 2.d0*dfxyz*diffinv + 2.d0 * curv * drxyz
        call elim_moment_fs(nat,rotforce(1,1))
!        call elim_torque_fs(nat,rxyz1(1,1),rotforce(1,1))
        call elim_torque_reza(nat,rxyz1(1,1),rotforce(1,1))

        !remove numerical noise:
!write(*,*)'rotnoise',ddot(3*nat,rotforce(1,1),1,vec(1,1),1)
!        rotforce = rotforce - ddot(3*nat,rotforce(1,1),1,vec(1,1),1)*vec
    else
        stop 'unknown method for curvature computation'
    endif
    
end subroutine


subroutine precondition(nat,check,pos,force)
use module_base
    implicit real(gp) (a-h,o-z)
    dimension pos(3,nat),force(3,nat),dis(3,nat),d(3,nat),spring(nat,nat)
    logical check

    fnrmout=0.d0
    do iat=1,nat
        fnrmout=fnrmout+force(1,iat)**2+force(2,iat)**2+force(3,iat)**2
    enddo
    fnrmout=sqrt(fnrmout)

    nprec=1000
    !beta=.5d-4
    beta=.3d-4
    do iat=1,nat
        dis(1,iat)=0.d0
        dis(2,iat)=0.d0
        dis(3,iat)=0.d0
    enddo

    ! calculate spring constants
    bondlength=1.0d0
    cspring0=1000.d0
    do iat=1,nat
        do jat=1,iat-1
            dd2=((pos(1,iat)-pos(1,jat))**2+(pos(2,iat)-pos(2,jat))**2&
               &+(pos(3,iat)-pos(3,jat))**2)/bondlength**2
            spring(iat,jat)=cspring0*exp(.5d0-.5d0*dd2)
        enddo
    enddo

    do iprec=1,nprec
        do iat=1,nat
            d(1,iat)=beta*force(1,iat)
            d(2,iat)=beta*force(2,iat)
            d(3,iat)=beta*force(3,iat)
            dis(1,iat)=dis(1,iat)+d(1,iat)
            dis(2,iat)=dis(2,iat)+d(2,iat)
            dis(3,iat)=dis(3,iat)+d(3,iat)
        enddo

        do iat=1,nat
            do jat=1,iat-1
                stretchx=d(1,iat)-d(1,jat)
                stretchy=d(2,iat)-d(2,jat)
                stretchz=d(3,iat)-d(3,jat)
                cspring=spring(iat,jat)
                force(1,iat)=force(1,iat)-cspring*stretchx
                force(2,iat)=force(2,iat)-cspring*stretchy
                force(3,iat)=force(3,iat)-cspring*stretchz
                force(1,jat)=force(1,jat)+cspring*stretchx
                force(2,jat)=force(2,jat)+cspring*stretchy
                force(3,jat)=force(3,jat)+cspring*stretchz
            enddo
        enddo
        ! write(24,'(i6,100(1x,e9.2))') iprec,force
        fnrm=0.d0
        do iat=1,nat
            fnrm=fnrm+force(1,iat)**2+force(2,iat)**2+force(3,iat)**2
        enddo
        fnrm=sqrt(fnrm)
        !write(*,*) 'iprec ',iprec,fnrm

        if (fnrm.lt.1.d-2*fnrmout) then
            if(check) write(100,*) 'iprec=',iprec,fnrm,fnrmout
            goto 100
        endif

    enddo
    if(check) write(100,*) 'FAIL iprec=',iprec,fnrm,fnrmout

    100 continue
    do iat=1,nat
        force(1,iat)=dis(1,iat)
        force(2,iat)=dis(2,iat)
        force(3,iat)=dis(3,iat)
    enddo

    return
end subroutine

!subroutine fix_fragmentation_fs(nat,rcov,rxyz,nputback,fixfragmented)
!    implicit real(gp) (a-h,o-z)
!    dimension rxyz(3,nat)
!    logical :: fixfragmented
!    real(gp), intent(in) :: rcov(nat)
!    ! automatic arrays
!    logical belong(nat)
!    fixfragmented=.true.
!
!    nloop=1
!    100  continue
!
!    iat=1
!    belong(iat)=.true.
!    ncluster=1
!    do iat=2,nat
!        belong(iat)=.false.
!    enddo
!
!    !   ic=0
!    do
!        nadd=0
!        do iat=1,nat
!            xi=rxyz(1,iat) ; yi=rxyz(2,iat) ; zi=rxyz(3,iat)
!            if (belong(iat)) then
!                do jat=1,nat
!                    bondlength=rcov(iat)+rcov(jat)
!                    xj=rxyz(1,jat) ; yj=rxyz(2,jat) ; zj=rxyz(3,jat)
!                    if ( (xi-xj)**2+(yi-yj)**2+(zi-zj)**2 .le. (bondlength*1.25d0)**2) then
!
!                        if (belong(jat).eqv. .false.) nadd=nadd+1
!                        belong(jat)=.true.
!                    endif
!                enddo
!            endif
!        enddo
!        ncluster=ncluster+nadd
!        !     ic=ic+1 ; write(100,*) 'nadd,ncluster',ic,nadd,ncluster
!        if (nadd.eq.0) exit
!    enddo
!
!    if (ncluster.eq.nat) then
!        !   write(100,*) 'No fragmentation has occured',nloop
!        fixfragmented=.false.
!        return
!    else
!        nputback=nputback+1
!
!        write(*,*) 'fragmentation occured',nloop,ncluster
!        !write(444,*) nat,iat
!        !write(444,*) ' 50. 0. 50. '
!        !write(444,*) ' 0. 0. 50. '
!        !do kat=1,nat
!        !write(444,*) rxyz(1,kat),rxyz(2,kat),rxyz(3,kat),'  LJ '
!        !enddo
!        
!        ! make sure the part that flew away is smaller tham the cluster
!        if (ncluster.le.nat/2) then
!            !     write(100,*) 'FLIP'
!            do iat=1,nat
!                if (belong(iat)) then
!                    belong(iat)=.false.
!                else
!                    belong(iat)=.true.
!                endif
!            enddo
!         endif
!
!        !find minimum distance between clusters
!        ii=-99999
!                jj=-99999
!        ddmin=1.d100
!        do iat=1,nat
!            if (belong(iat).eqv. .false.) then
!                xi=rxyz(1,iat) ; yi=rxyz(2,iat) ; zi=rxyz(3,iat)
!                do jat=1,nat
!                    if (belong(jat)) then
!                        xj=rxyz(1,jat) ; yj=rxyz(2,jat) ; zj=rxyz(3,jat)
!                        dd= (xi-xj)**2+(yi-yj)**2+(zi-zj)**2
!                        if (dd.lt.ddmin) then
!                            jj=jat
!                            ii=iat
!!write(*,*)'index',ii,jj
!                            ddmin=dd
!                        endif
!                    endif
!                enddo
!            endif
!        enddo
!!write(*,*)'index',ii,jj
!
!        ! pull back the fragment of atoms that flew away
!        d1=rxyz(1,ii)-rxyz(1,jj)
!        d2=rxyz(2,ii)-rxyz(2,jj)
!        d3=rxyz(3,ii)-rxyz(3,jj)
!        bondlength=rcov(ii)+rcov(jj)
!        tt=bondlength/sqrt(d1**2+d2**2+d3**2)
!        do iat=1,nat
!            if (belong(iat).eqv. .false.) then
!                rxyz(1,iat)=rxyz(1,iat)+d1*(tt-1.0d0)
!                rxyz(2,iat)=rxyz(2,iat)+d2*(tt-1.0d0)
!                rxyz(3,iat)=rxyz(3,iat)+d3*(tt-1.0d0)
!            endif
!        enddo
!    
!        !write(444,*) nat, ' fixed configuration'
!        !write(444,*) ' 50. 0. 50. '
!        !write(444,*) ' 0. 0. 50. '
!        !do iat=1,nat
!        !write(444,*) rxyz(1,iat),rxyz(2,iat),rxyz(3,iat),'  LJ '
!        !enddo
!        nloop=nloop+1
!        goto 100
!    endif
!end subroutine



subroutine fixfrag_posvel(nat,rcov,pos,vel,option,occured)
!UNTERSCHIED ZUR MH ROUTINE:
!pos(:,iat)=pos(:,iat)+vec(:)*((mindist-1.0d0*bondlength)/mindist)
!ANSTATT
!pos(:,iat)=pos(:,iat)+vec(:)*((mindist-1.5d0*bondlength)/mindist)
!MACHT GROSSEN PERFORMANCE UNTERSCHIED AUS
!

!This subroutine can perform two tasks.
!ATTENTION: it will only work on free BC!!!
!
!option=1
!The atoms in pos are analyzed and, if there is a fragmentation occuring, the
!main fragment will be identified and all neighboring fragments will be moved towards the nearest
!atom of the main fragment. The array pos will then be updated and returned.
!
!option=2
!The fragments are identified and the center of mass of all fragments are computed separately.
!The center of mass of all cluster is also computed.
!Then, the velocities are modified in such a way that the projection of the velocities 
!along the vector pointing towards the center of mass of all fragments are inverted 
!!use module_base
!!use module_types
!!use m_ab6_symmetry
use module_base
implicit none
integer, intent(in) :: nat
!type(atoms_data), intent(in) :: at
real(gp),dimension(3,nat), INTENT(INOUT) :: pos
real(gp),dimension(3,nat), INTENT(INOUT) :: vel
real(gp),dimension(nat), INTENT(IN) :: rcov
integer, INTENT(IN):: option
integer :: nfrag, nfragold
logical :: occured,niter
real(gp)::  dist, mindist, angle, vec(3), cmass(3), velcm(3), bondlength, bfactor,rnrmi,scpr
real(gp):: ekin,vcm1,vcm2,vcm3,ekin0,scale
real(gp), allocatable:: cm_frags(:,:), vel_frags(:,:)
integer::iat, jat, natfragx(1), imin(2),ifrag
integer, allocatable:: fragcount(:)
integer, allocatable:: nat_frags(:)
integer, dimension(nat):: fragarr
logical, allocatable:: invert(:)

!The bondlength (in atomic units) is read from file input.bondcut
! OPTION 1: System is considered to be fragmented if the minimal distance between two atoms in the fragment is more than 2.0*bondlength
!           The two fragment are then brought together such that the minimal distance equals 1.5*bondlength
! OPTION  : System is considered to be fragmented if the minimal distance between two atoms in the fragment is more than 2.0*bondlength
!           the velocities are then inverted
!open(unit=43,file="input.bondcut")
!read(43,*) bondlength
!close(43)

if (option == 1) then 
   bfactor=1.5d0
else if (option == 2) then 
   bfactor=2.d0
else
   stop 'wrong option'
endif



fragarr(:)=0                     !Array, which atom belongs to which fragment
nfrag=0                       !Number of fragments

!Check for correct input
if (option.ne.1 .and. option.ne.2) stop "Wrong option in fixfrag_refvels"

!Calculate number of fragments and fragmentlist of the atoms
loop_nfrag: do
   nfragold=nfrag
   do iat=1,nat                !Check the first atom that isn't part of a cluster yet
      if(fragarr(iat)==0) then
         nfrag=nfrag+1
         fragarr(iat)=nfrag
         exit 
      endif
   enddo
   if (nfragold==nfrag) exit loop_nfrag
7000 niter=.false.
   do iat=1,nat                !Check if all the other atoms are part of the current cluster
      do jat=1,nat
         bondlength=rcov(iat)+rcov(jat)
         if(nfrag==fragarr(iat) .AND. jat.ne.iat .AND. fragarr(jat)==0) then
            dist=(pos(1,iat)-pos(1,jat))**2+(pos(2,iat)-pos(2,jat))**2+(pos(3,iat)-pos(3,jat))**2
            if(dist<(bfactor*bondlength)**2) then
               fragarr(jat)=nfrag
               niter=.true.
            endif
         endif
      enddo
   enddo
   if(niter) then
      goto 7000
   endif
end do loop_nfrag


!   if(iproc==0) write(*,*) '(MH) nfrag=',nfrag
occured=.false.
if(nfrag.ne.1) then          !"if there is fragmentation..."
   occured=.true.
write(671,*)'fragmentation occured'
!   if(iproc==0) then
!      call yaml_open_map('(MH) FIX')
!      call yaml_map('(MH) Number of Fragments counted with option', (/nfrag,option/))
!   endif
   if (option==1) then !OPTION=1, FIX FRAGMENTATION
      !   if(nfrag.ne.1) then          !"if there is fragmentation..."

      !Find out which fragment is the main cluster
      allocate(fragcount(nfrag))
      fragcount=0
      do ifrag=1,nfrag
         do iat=1,nat
            if(fragarr(iat)==ifrag) then
               fragcount(ifrag)=fragcount(ifrag)+1
            endif
         enddo
      enddo
      natfragx=maxloc(fragcount(:))
      !if(iproc==0) call yaml_map('(MH) The main Fragment index is', natfragx(1))

      !Find the minimum distance between the clusters
      do ifrag=1,nfrag
         mindist=1.d100
         if(ifrag.ne.natfragx(1)) then
            do iat=1,nat
               if(fragarr(iat)==ifrag) then
                  do jat=1,nat
                     if(fragarr(jat)==natfragx(1)) then
                        dist=(pos(1,iat)-pos(1,jat))**2+(pos(2,iat)-pos(2,jat))**2+(pos(3,iat)-pos(3,jat))**2
                        if(dist<mindist**2) then
                           mindist=sqrt(dist)
                           imin(1)=jat  !Atom with minimal distance in main fragment
                           imin(2)=iat   !Atom with minimal distance in fragment ifrag
                        endif
                     endif
                  enddo
               endif
            enddo

!            if (iproc == 0) then
!               write(444,*) nat, 'atomic '
!               write(444,*) 'A fragmented configuration ',imin(1),imin(2)
!               do iat=1,nat
!                  write(444,'(a5,3(e15.7),l1)') ' Mg  ',pos(1,iat),pos(2,iat),pos(3,iat)
!               enddo
!            endif


            vec(:)=pos(:,imin(1))-pos(:,imin(2))
            bondlength=rcov(imin(1))+rcov(imin(2))
            do iat=1,nat        !Move fragments back towards the main fragment 
               if(fragarr(iat)==ifrag) then
                  pos(:,iat)=pos(:,iat)+vec(:)*((mindist-1.0d0*bondlength)/mindist)
                  fragarr(iat)=natfragx(1)
               endif
            enddo


         endif
      enddo
      deallocate(fragcount)
!      if(iproc==0) then
!         call yaml_comment('(MH) FIX: Fragmentation fixed! Keep on hopping...')
!         call yaml_close_map()
!      end if
!      if (iproc == 0) then
!         write(444,*) nat, 'atomic '
!         write(444,*) ' fixed configuration '
!         do iat=1,nat
!            write(444,'(a5,3(e15.7),l1)') ' Mg  ',pos(1,iat),pos(2,iat),pos(3,iat)
!         enddo
!      endif

      !   endif
   elseif(option==2) then !OPTION=2, INVERT VELOCITIES
      !   if(nfrag.ne.1) then          !"if there is fragmentation..."
!      if(iproc==0) call yaml_map('(MH) FIX: Preparing to invert velocities, option:',option)
      !Compute center of mass of all fragments and the collectiove velocity of each fragment
      allocate(cm_frags(3,nfrag),vel_frags(3,nfrag),nat_frags(nfrag))
      allocate(invert(nfrag))
      cm_frags(:,:)=0.d0
      vel_frags(:,:)=0.d0
      nat_frags(:)=0         !number of atoms per fragment
      cmass(:)=0.d0
      velcm(:)=0.d0
      do iat=1,nat
         ifrag=fragarr(iat)
         nat_frags(ifrag)=nat_frags(ifrag)+1
         cm_frags(:,ifrag)=cm_frags(:,ifrag)+pos(:,iat)
         vel_frags(:,ifrag)=vel_frags(:,ifrag)+vel(:,iat)
      enddo

      do ifrag=1,nfrag
         cm_frags(:,ifrag)=cm_frags(:,ifrag)/real(nat_frags(ifrag),8)
         vel_frags(:,ifrag)=vel_frags(:,ifrag)/real(nat_frags(ifrag),8)
         cmass(:)=cmass(:)+cm_frags(:,ifrag)*nat_frags(ifrag)/real(nat,8)
         velcm(:)=velcm(:)+vel_frags(:,ifrag)*nat_frags(ifrag)/real(nat,8)
      enddo
!      if (iproc==0) call yaml_map('(MH) CM VELOCITY',sqrt(velcm(1)**2+velcm(2)**2+velcm(3)**2))
!      if (velcm(1)**2+velcm(2)**2+velcm(3)**2.gt.1.d-24) then
!         if (iproc==0) call yaml_comment('(MH) NONZERO CM VELOCITY')
!      endif


      ! now cm_frags contains the unit vector pointing from the center of mass of the entire system to the center of mass of the fragment
      do ifrag=1,nfrag
         cm_frags(:,ifrag)=cm_frags(:,ifrag)-cmass(:)
         rnrmi=1.d0/sqrt(cm_frags(1,ifrag)**2+cm_frags(2,ifrag)**2+cm_frags(3,ifrag)**2)
         cm_frags(1,ifrag)=cm_frags(1,ifrag)*rnrmi
         cm_frags(2,ifrag)=cm_frags(2,ifrag)*rnrmi
         cm_frags(3,ifrag)=cm_frags(3,ifrag)*rnrmi
         angle=cm_frags(1,ifrag)*vel_frags(1,ifrag)+cm_frags(2,ifrag)*vel_frags(2,ifrag)+cm_frags(3,ifrag)*vel_frags(3,ifrag)
         rnrmi=1.d0/sqrt(vel_frags(1,ifrag)**2+vel_frags(2,ifrag)**2+vel_frags(3,ifrag)**2)
         angle=angle*rnrmi
         if (angle.gt.0.d0) then
            invert(ifrag)=.true.
         else
            invert(ifrag)=.false.
         endif
!         if (iproc==0) then
!           write(*,*) '(MH) ifrag, angle ',ifrag, angle,invert(ifrag)
!           call yaml_open_map('(MH) Frag. Info',flow=.true.)
!            call yaml_map('ifrag',ifrag)
!            call yaml_map('angle',angle)
!            call yaml_map('ifrag inverted',invert(ifrag))
!           call yaml_close_map(advance='yes')
!         endif
      enddo
      !Decompose each atomic velocity into an component parallel and perpendicular to the cm_frags  vector and inter the 
      !paralle part if it point away from the CM

      !Check kinetic energy before inversion
      ekin0=0.d0
      vcm1=0.d0
      vcm2=0.d0
      vcm3=0.d0
      do iat=1,nat
         ekin0=ekin0+vel(1,iat)**2+vel(2,iat)**2+vel(3,iat)**2
         vcm1=vcm1+vel(1,iat)
         vcm2=vcm2+vel(2,iat)
         vcm3=vcm3+vel(3,iat)
      enddo
!      if (iproc==0) then
!          write(*,'(a,e14.7,3(e10.3))') '(MH) EKIN CM before invert',ekin0,vcm1,vcm2,vcm3
!!          call yaml_open_map(,flow=.true.)
!          call yaml_map('(MH) EKIN CM before invert',(/ekin0,vcm1,vcm2,vcm3/),fmt='(e10.3)')
!!          call yaml_close_map(advance='yes')
!      endif
!      if (iproc==0) call torque(nat,pos,vel)
      !Checkend kinetic energy before inversion

      do iat=1,nat
         ! inversions  by fragment group
         ifrag=fragarr(iat)
         if (invert(ifrag)) then
            scpr=cm_frags(1,ifrag)*vel(1,iat)+cm_frags(2,ifrag)*vel(2,iat)+cm_frags(3,ifrag)*vel(3,iat)
            vel(:,iat)=vel(:,iat)-scpr*cm_frags(:,ifrag)*2.d0
         endif
      enddo

      call elim_moment_fs(nat,vel)
!      call elim_torque_fs(nat,pos,vel)
      call elim_torque_reza(nat,pos,vel)

      ! scale velocities to regain initial ekin0
      ekin=0.d0
      do iat=1,nat
         ekin=ekin+vel(1,iat)**2+vel(2,iat)**2+vel(3,iat)**2
      enddo
      scale=sqrt(ekin0/ekin)
      do iat=1,nat
         vel(1,iat)=vel(1,iat)*scale
         vel(2,iat)=vel(2,iat)*scale
         vel(3,iat)=vel(3,iat)*scale
      enddo

      !Check kinetic energy after inversion
      ekin=0.d0
      vcm1=0.d0
      vcm2=0.d0
      vcm3=0.d0
      do iat=1,nat
         ekin=ekin+vel(1,iat)**2+vel(2,iat)**2+vel(3,iat)**2
         vcm1=vcm1+vel(1,iat)
         vcm2=vcm2+vel(2,iat)
         vcm3=vcm3+vel(3,iat)
      enddo
!      if (iproc==0) then
!          write(*,'(a,e14.7,3(e10.3))') '(MH) EKIN CM after  invert',ekin,vcm1,vcm2,vcm3
!          !call yaml_open_map('(MH) EKIN CM after invert',flow=.true.)
!          call yaml_map('(MH) EKIN CM after invert',(/ekin0,vcm1,vcm2,vcm3/),fmt='(e10.3)')
!          !call yaml_close_map(advance='yes')
!      endif
!      if (iproc==0) call torque(nat,pos,vel)
      !Checkend kinetic energy after inversion

      !Check angle  after inversion
      vel_frags(:,:)=0.d0
      do iat=1,nat
         ifrag=fragarr(iat)
         vel_frags(:,ifrag)=vel_frags(:,ifrag)+vel(:,iat)
      enddo
      do ifrag=1,nfrag
         angle=cm_frags(1,ifrag)*vel_frags(1,ifrag)+cm_frags(2,ifrag)*vel_frags(2,ifrag)+cm_frags(3,ifrag)*vel_frags(3,ifrag)
         rnrmi=1.d0/sqrt(vel_frags(1,ifrag)**2+vel_frags(2,ifrag)**2+vel_frags(3,ifrag)**2)
         angle=angle*rnrmi
!         if (iproc==0) then
!           call yaml_open_map('(MH) Frag',flow=.true.)
!            call yaml_map('ifrag',ifrag)
!            call yaml_map('angle a invert',angle)
!           call yaml_close_map(advance='yes')
!         endif
        
      enddo
      !Checkend kinetic energy after inversion


      deallocate(cm_frags,vel_frags,nat_frags)
      deallocate(invert)
      !   endif
      !else
      !   stop "Wrong option within ff-rv"
!      if(iproc==0) write(*,*) "(MH) FIX: Velocity component towards the center of mass inverted! Keep on hopping..."
!      if(iproc==0) call yaml_scalar('(MH) FIX: Velocity component towards the center of mass inverted! Keep on hopping...')
   endif
endif
end subroutine fixfrag_posvel

!*****************************************************************************************
!BELOW HERE ONLY REZAS PRECOND. ROUTINES
subroutine preconditioner_reza(iproc,nr,nat,rat,atomnames,fat)
use module_base
   implicit none
   integer :: iproc,nr,nat,nsb,nrsqtwo,i,j,k,info
   real(gp) :: rat(3,nat),fat(nr),soft,hard
   character(20):: atomnames(nat)
   real(gp), allocatable::evec(:,:),eval(:),wa(:),hess(:,:)
   real(gp) :: dx,dy,dz,r,tt, DDOT

   allocate(hess(nr,nr))
   call make_hessian(iproc,nr,nat,rat,atomnames,hess,nsb)

   nrsqtwo=2*nr**2
   allocate(evec(nr,nr),eval(nr),wa(nrsqtwo))
   evec(1:nr,1:nr)=hess(1:nr,1:nr)
   !if(iproc==0) write(*,*) 'HESS ',hess(:,:)

   call DSYEV('V','L',nr,evec,nr,eval,wa,nrsqtwo,info)
   if(info/=0) stop 'ERROR: DSYEV in inithess failed.'

   if (iproc==0) then
       do i=1,nr
           write(*,'(i5,es20.10)') i,eval(i)
       enddo
   endif

   hard=eval(nr)
   soft=eval(nr-nsb+1)
   write(*,*) 'SOFT in preconditioner_reza',soft,nr,nr-nsb+1
   do k=1,nr
       wa(k)=DDOT(nr,fat,1,evec(1,k),1)
       !write(*,'(i5,2es20.10)') k,eval(k),wa(k)
   enddo
   do i=1,nr
      tt=0.d0
      do j=1,nr
           tt=tt+(wa(j)/sqrt(eval(j)**2+ 66.d0**2))*evec(i,j)
           !tt=tt+(wa(j)/sqrt(eval(j)**2+soft**2))*evec(i,j)
      enddo
      fat(i)=tt
   enddo
   deallocate(evec,eval,wa,hess)
end subroutine preconditioner_reza
subroutine make_hessian(iproc,nr,nat,rat,atomnames,hess,nsb)
use module_base
   implicit none
   integer :: iproc,nr,nat,iat,jat,nsb,i,j,k,info
   real(gp) :: rat(3,nat),hess(nr,nr),r0types(4,4),fctypes(4,4)
   character(20):: atomnames(nat)
   integer, allocatable::ita(:),isb(:,:)
   real(gp), allocatable::r0bonds(:),fcbonds(:)
   real(gp) :: dx,dy,dz,r,tt

   if(nr/=3*nat) then
       write(*,'(a)') 'This subroutine works only for systems without fixed atoms.'
       stop
   endif
   allocate(ita(nat),isb(10*nat,2),r0bonds(10*nat),fcbonds(10*nat))
   do iat=1,nat
   !write(*,*) trim(atomnames(iat))
   !stop
       if(trim(atomnames(iat))=='H') then
           ita(iat)=1
       elseif(trim(atomnames(iat))=='C') then
           ita(iat)=2
       elseif(trim(atomnames(iat))=='N') then
           ita(iat)=3
       elseif(trim(atomnames(iat))=='O') then
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
               r0bonds(nsb)=r0types(ita(iat),ita(jat)) !CAUTION: equil. bond le >  from amber
               !r0bonds(nsb)=r !CAUTION: current bond le >  assumed as equil. 
               fcbonds(nsb)=fctypes(ita(iat),ita(jat))
           endif
       enddo
   enddo
   if(iproc==0) write(*,*) 'NSB ',nsb
   !if(iproc==0) then
   !    do i=1,nsb
   !        write(*,'(a,i5,2f20.10,2i4,2(x,a))') 'PAR ', &
   !            i,r0bonds(i),fcbonds(i),isb(i,1),isb(i,2), &
   !            trim(atoms%astruct%atomnames(atoms%astruct%iatype(isb(i,1)))),trim(atoms%astruct%atomnames(atoms%astruct%iatype(isb(i,2))))
   !    enddo
   !endif
   call pseudohess(nat,rat,nsb,isb(1,1),isb(1,2),fcbonds,r0bonds,hess)
   deallocate(ita,isb,r0bonds,fcbonds)
end subroutine make_hessian
!*****************************************************************************************
subroutine init_parameters(r0,fc)
use module_base
    implicit none
    integer :: i,j
    real(gp) :: r0(4,4),fc(4,4)
    !real(gp):: units_r=1.d0/0.529d0
    real(gp):: units_r=1.d0
    !real(gp):: units_k=4.48d-4
    real(gp):: units_k=1.d0
    !((0.0104 / 0.239) / 27.2114) * (0.529177^2) = 0.000447802457
    r0(1,1)=0.80d0*units_r
    r0(2,1)=1.09d0*units_r ; r0(2,2)=1.51d0*units_r
    r0(3,1)=1.01d0*units_r ; r0(3,2)=1.39d0*units_r ; r0(3,3)=1.10d0*units_r
    r0(4,1)=0.96d0*units_r ; r0(4,2)=1.26d0*units_r ; r0(4,3)=1.10d0*units_r ; r0(4,4)=1.10*units_r
    do i=1,4
        do j=i+1,4
            r0(i,j)=r0(j,i)
        enddo
    enddo
    fc(1,1)=1.00d3*units_k
    fc(2,1)=3.40d2*units_k ; fc(2,2)=3.31d2*units_k
    fc(3,1)=4.34d2*units_k ; fc(3,2)=4.13d2*units_k ; fc(3,3)=4.56d3*units_k
    fc(4,1)=5.53d2*units_k ; fc(4,2)=5.43d2*units_k ; fc(4,3)=4.56d3*units_k ; fc(4,4)=4.56d3*units_k
    do i=1,4
        do j=i+1,4
            fc(i,j)=fc(j,i)
        enddo
    enddo
end subroutine init_parameters
!*****************************************************************************************
subroutine pseudohess(nat,rat,nbond,indbond1,indbond2,sprcons,xl0,hess)
use module_base
    implicit none
    integer :: nat,nbond,indbond1(nbond),indbond2(nbond)
    real(gp) :: rat(3,nat),sprcons(nbond),xl0(nbond),hess(3*nat,3*nat)
    integer :: iat,jat,i,j,ibond
    real(gp) :: dx,dy,dz,r2,r,rinv! n(c) r3inv
!,rinv2,rinv4,rinv8,rinv10,rinv14,rinv16
    real(gp) :: dxsq,dysq,dzsq,dxdy,dxdz,dydz,tt1,tt2,tt3
    real(gp) :: h11,h22,h33,h12,h13,h23
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
!*****************************************************************************************

subroutine minenergyandforces(imode,nat,alat,rat,rxyzraw,fat,fstretch,&
           fxyzraw,epot,iconnect,nbond_,atomnames,wold,alpha_stretch0,alpha_stretch)
use module_base
use module_energyandforces
    implicit real(gp) (a-h,o-z)
    dimension :: rat(3,nat),fat(3,nat),iconnect(2,nbond_),fstretch(3,nat),fxyzraw(3,nat),rxyzraw(3,nat)
    dimension :: ss(nbond_,nbond_),w(nbond_),vv(3,nat,nbond_),wold(nbond_)
    real(gp), intent(in) :: alpha_stretch0
    real(gp), intent(in) :: alat(3)
    character(100) filename
    character(20):: atomnames(nat)

    rxyzraw=rat
    call energyandforces(nat,alat,rat,fat,epot)
    fxyzraw=fat
    fnrmold=dnrm2(3*nat,fat,1)
    fstretch=0.d0

    if(imode==2)then
        call projectbond(nat,nbond_,rat,fat,fstretch,iconnect,atomnames,wold,alpha_stretch0,alpha_stretch)
    endif
 
end subroutine minenergyandforces
subroutine mincurvgrad(imode,nat,alat,diff,rxyz1,fxyz1,vec,vecraw,&
           rotforce,rotfstretch,rotforceraw,curv,imethod,ec,&
           iconnect,nbond_,atomnames,wold,alpha_stretch0,alpha_stretch)
use module_base
    implicit real(gp) (a-h,o-z)
    integer, intent(in) :: imode
    real(gp) :: alat(3)
    real(gp) :: rxyz1(3,nat),fxyz1(3,nat)
    real(gp) :: rxyz2(3,nat),fxyz2(3,nat)
    real(gp) :: vec(3,nat), vecraw(3,nat)
    real(gp) :: rotforce(3,nat),rotforceraw(3,nat)
    real(gp) :: rotfstretch(3,nat)
    real(gp) :: wold(nbond_)
    real(gp) :: curv,ec
    real(gp), intent(in) :: alpha_stretch0
    integer :: imethod,nbond_,iconnect(2,nbond_)
    character(20):: atomnames(nat)

     vecraw=vec
!    call energyandforces(nat,rat,fat,epot)
     call curvgrad(nat,alat,diff,rxyz1,fxyz1,vec,curv,rotforce,imethod,ec)
     rxyz2 =rxyz1+diff*vec
     rotforceraw=rotforce
     rotfstretch=0.d0
     fnrmold=dnrm2(3*nat,rotforce,1)

     if(imode==2)then
         call projectbond(nat,nbond_,rxyz2,rotforce,rotfstretch,iconnect,atomnames,wold,alpha_stretch0,alpha_stretch)
     endif

end subroutine mincurvgrad

subroutine projectbond(nat,nbond,rat,fat,fstretch,iconnect,atomnames,wold,alpha_stretch0,alpha_stretch)
use module_base
    implicit none
    integer, intent(in) :: nat
    integer, intent(in) :: nbond
    real(gp), intent(in) :: rat(3,nat)
    real(gp), intent(inout) :: fat(3,nat)
    real(gp), intent(inout) :: fstretch(3,nat)
    integer, intent(in) :: iconnect(2,nbond)
    character(20), intent(in) :: atomnames(nat)
    real(gp), intent(inout) :: wold(nbond)
    real(gp), intent(in) :: alpha_stretch0
    real(gp), intent(inout) :: alpha_stretch
    !internal
    integer :: iat,jat,ibond,jbond,l,nsame,info
    real(gp) :: ss(nbond,nbond),w(nbond),vv(3,nat,nbond)
    real(gp) :: per
    !functions
    real(gp) :: ddot

    fstretch=0.d0

! set up positional overlap matrix
     vv=0.d0
     do ibond=1,nbond
     iat=iconnect(1,ibond)
     jat=iconnect(2,ibond)
     do l=1,3
     vv(l,iat,ibond)=rat(l,jat)-rat(l,iat)
     vv(l,jat,ibond)=rat(l,iat)-rat(l,jat)
     enddo
     enddo

     do ibond=1,nbond
     do jbond=1,nbond
     ss(ibond,jbond)=ddot(3*nat,vv(1,1,ibond),1,vv(1,1,jbond),1)
     enddo
     w(ibond)=ddot(3*nat,vv(1,1,ibond),1,fat,1)
     enddo

     nsame=0
     do ibond=1,nbond
      if ( wold(ibond)*w(ibond).gt.0.d0) nsame=nsame+1
         wold(ibond)=w(ibond)
     enddo
     per=float(nsame)/nbond
     if (per.gt. .66d0) then
         alpha_stretch=alpha_stretch*1.10d0
!         alpha_stretch=min(100.d0*alpha_stretch0,alpha_stretch*1.10d0)
!         alpha_stretch=min(2.5d0*alpha_stretch0,alpha_stretch*1.10d0)
!         alpha_stretch=min(1.d-3,alpha_stretch*1.10d0)
     else
         alpha_stretch=max(1.d-2*alpha_stretch0,alpha_stretch/1.10d0)
!         alpha_stretch=max(1.d-10,alpha_stretch/1.10d0)
     endif
     write(555,'(f5.1,a,1x,es10.3)') 100*per,' percent of bond directions did not switch sign',alpha_stretch

     call DPOSV('L', nbond, 1, ss, nbond, w, nbond, info )
     if (info.ne.0) then
        stop 'info DPOSV in minenergyforces'
     endif

! calculate projected force
     do iat=1,nat
     do l=1,3
     fstretch(l,iat)=0.d0
     enddo
     enddo

     do ibond=1,nbond
     do iat=1,nat
     do l=1,3
     fstretch(l,iat)=fstretch(l,iat)+w(ibond)*vv(l,iat,ibond)
     enddo
     enddo
     enddo
!     fnrmst=dnrm2(3*nat,fstretch,1)
     do iat=1,nat
     do l=1,3
     fat(l,iat)=fat(l,iat)-fstretch(l,iat)
     enddo
     enddo

!     fnrm=dnrm2(3*nat,fat,1)
!     write(555,'(a,i4,3(e11.3),2x,e17.10)') 'fnrmold,fnrm,fnrmst',istretch,fnrmold,fnrm,fnrmst,epot


end subroutine
subroutine elim_torque_reza(nat,rat0,fat)
use module_base
!  use module_base
  implicit none
  integer, intent(in) :: nat
  real(gp), dimension(3*nat), intent(in) :: rat0
  real(gp), dimension(3*nat), intent(inout) :: fat
  !local variables
  character(len=*), parameter :: subname='elim_torque_reza'
  integer :: i,iat,i_all,i_stat
  real(gp) :: vrotnrm,cmx,cmy,cmz,alpha,totmass
  !this is an automatic array but it should be allocatable
  real(gp), dimension(3) :: evaleria
  real(gp), dimension(3,3) :: teneria
  real(gp), dimension(3*nat) :: rat
  real(gp), dimension(3*nat,3) :: vrot
  real(gp), dimension(:), allocatable :: amass
real(gp) :: dnrm2

  allocate(amass(nat),stat=i_stat)
!  call memocc(i_stat,amass,'amass',subname)

  rat=rat0
  amass(1:nat)=1.0d0
  !project out rotations
  totmass=0.0d0
  cmx=0.0d0
  cmy=0.0d0
  cmz=0.0d0
  do i=1,3*nat-2,3
     iat=(i+2)/3
     cmx=cmx+amass(iat)*rat(i+0)
     cmy=cmy+amass(iat)*rat(i+1)
     cmz=cmz+amass(iat)*rat(i+2)
     totmass=totmass+amass(iat)
  enddo
  cmx=cmx/totmass
  cmy=cmy/totmass
  cmz=cmz/totmass
  do i=1,3*nat-2,3
     rat(i+0)=rat(i+0)-cmx
     rat(i+1)=rat(i+1)-cmy
     rat(i+2)=rat(i+2)-cmz
  enddo

  call moment_of_inertia(nat,rat,teneria,evaleria)
  do iat=1,nat
     i=iat*3-2
     call cross(teneria(1,1),rat(i),vrot(i,1))
     call cross(teneria(1,2),rat(i),vrot(i,2))
     call cross(teneria(1,3),rat(i),vrot(i,3))
  enddo
  call normalizevector(3*nat,vrot(1,1))
  call normalizevector(3*nat,vrot(1,2))
  call normalizevector(3*nat,vrot(1,3))

  do i=1,3*nat-2,3
     rat(i+0)=rat(i+0)+cmx
     rat(i+1)=rat(i+1)+cmy
     rat(i+2)=rat(i+2)+cmz
  enddo

  vrotnrm=dnrm2(3*nat,vrot(1,1),1)
  if (vrotnrm /= 0.0d0) vrot(1:3*nat,1)=vrot(1:3*nat,1)/vrotnrm
  vrotnrm=dnrm2(3*nat,vrot(1,2),1)
  if (vrotnrm /= 0.0d0) vrot(1:3*nat,2)=vrot(1:3*nat,2)/vrotnrm
  vrotnrm=dnrm2(3*nat,vrot(1,3),1)
  if (vrotnrm /= 0.0d0) vrot(1:3*nat,3)=vrot(1:3*nat,3)/vrotnrm

  do i=1,3
     alpha=0.0d0
     if(abs(evaleria(i)).gt.1.d-10) then
        alpha=dot_product(vrot(:,i),fat(:))
        fat(:)=fat(:)-alpha*vrot(:,i)
     endif
  enddo

  i_all=-product(shape(amass))*kind(amass)
  deallocate(amass,stat=i_stat)
!  call memocc(i_stat,i_all,'amass',subname)

END SUBROUTINE elim_torque_reza
subroutine cross(a,b,c)
use module_base
!  use module_base
  implicit none
  real(gp), dimension(3), intent(in) :: a,b
  real(gp), dimension(3), intent(out) :: c

  c(1)=a(2)*b(3)-b(2)*a(3)
  c(2)=a(3)*b(1)-b(3)*a(1)
  c(3)=a(1)*b(2)-b(1)*a(2)
END SUBROUTINE cross


subroutine moment_of_inertia(nat,rat,teneria,evaleria)
use module_base
!  use module_base
  implicit none
  integer, intent(in) :: nat
  real(gp), dimension(3,nat), intent(in) :: rat
  real(gp), dimension(3), intent(out) :: evaleria
  real(gp), dimension(3,3), intent(out) :: teneria
  !local variables
  character(len=*), parameter :: subname='moment_of_inertia'
  integer, parameter::lwork=100
  integer :: iat,info,i_all,i_stat
  real(gp) :: tt
  real(gp), dimension(lwork) :: work
  real(gp), dimension(:), allocatable :: amass

  allocate(amass(nat),stat=i_stat)
!  call memocc(i_stat,amass,'amass',subname)

  !positions relative to center of geometry
  amass(1:nat)=1.0d0
  !calculate inertia tensor
  teneria(1:3,1:3)=0.0d0
  do iat=1,nat
     tt=amass(iat)
     teneria(1,1)=teneria(1,1)+tt*(rat(2,iat)*rat(2,iat)+rat(3,iat)*rat(3,iat))
     teneria(2,2)=teneria(2,2)+tt*(rat(1,iat)*rat(1,iat)+rat(3,iat)*rat(3,iat))
     teneria(3,3)=teneria(3,3)+tt*(rat(1,iat)*rat(1,iat)+rat(2,iat)*rat(2,iat))
     teneria(1,2)=teneria(1,2)-tt*(rat(1,iat)*rat(2,iat))
     teneria(1,3)=teneria(1,3)-tt*(rat(1,iat)*rat(3,iat))
     teneria(2,3)=teneria(2,3)-tt*(rat(2,iat)*rat(3,iat))
     teneria(2,1)=teneria(1,2)
     teneria(3,1)=teneria(1,3)
     teneria(3,2)=teneria(2,3)
  enddo
  !diagonalize inertia tensor
  call DSYEV('V','L',3,teneria,3,evaleria,work,lwork,info)
  i_all=-product(shape(amass))*kind(amass)
  deallocate(amass,stat=i_stat)
!  call memocc(i_stat,i_all,'amass',subname)

END SUBROUTINE moment_of_inertia


subroutine normalizevector(n,v)
use module_base
!  use module_base
  implicit none
  integer, intent(in) :: n
  real(gp), dimension(n), intent(inout) :: v
  !local variables
  integer :: i
  real(gp) :: vnrm

  vnrm=0.0d0
  do i=1,n
     vnrm=vnrm+v(i)**2
  enddo
  vnrm=sqrt(vnrm)
  if (vnrm /= 0.0d0) v(1:n)=v(1:n)/vnrm

END SUBROUTINE normalizevector
end module
