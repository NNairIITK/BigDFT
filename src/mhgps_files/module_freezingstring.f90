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
!=====================================================================
subroutine get_ts_guess(nat,alat,rxyz1,rxyz2,tsguess,minmodeguess)
    use module_base
    use yaml_output
    use module_global_variables,&
            only: iproc, mhgps_verbosity,&
                  gammainv => ts_guess_gammainv,&
                  perpnrmtol => ts_guess_perpnrmtol,&
                  trust => ts_guess_trust,&
                  nstepsmax => ts_guess_nstepsmax
    use module_energyandforces
    implicit none
    !parameters
    integer, intent(in) :: nat
    real(gp), intent(in) :: alat(3)
    real(gp), intent(in) :: rxyz1(3,nat),rxyz2(3,nat)
    real(gp), intent(out) :: tsguess(3,nat),minmodeguess(3,nat)
    !internal
    integer :: finished
    integer :: i,j,iat
    integer :: npath,istring,istringmax,isidemax,ipathmax
    integer :: nnodes,nstring,nstringmax=20
    integer :: ihigh,ncorr
    real(gp), allocatable :: string(:,:,:)
    real(gp), allocatable :: path(:,:,:)
    real(gp), allocatable :: forces(:,:,:)
    real(gp), allocatable :: energies(:)
    real(gp), allocatable :: arc(:)
    real(gp), allocatable :: y2vec(:,:,:)
    real(gp), allocatable :: tangent(:,:,:)
    real(gp) :: fnoise
    real(gp) :: emax
    real(gp) :: tau,rdmy
    real(gp) :: yp1=huge(1._gp), ypn=huge(1._gp)!natural splines
    !functions
    real(gp) :: dnrm2

    string = f_malloc((/ 1.to.3*nat, 1.to.2, 1.to.nstringmax/),&
                     id='string')
    call vcopy(3*nat, rxyz1(1,1), 1, string(1,1,1), 1)
    call vcopy(3*nat, rxyz2(1,1), 1, string(1,2,1), 1)
    
    call grow_string(nat,alat,gammainv,perpnrmtol,trust,nstepsmax,&
                    nstringmax,nstring,string,finished)

    nnodes=2*nstring
    ncorr =0
    if(finished==1)then
    !if(string(1,2,nstring)>=1.e30_gp)then
        nnodes=nnodes-1
        ncorr = 1
    endif
    path = f_malloc((/1.to.nnodes,1.to.3,1.to.nat/),id='path')
    y2vec = f_malloc((/1.to.nnodes,1.to.3,1.to.nat/),id='y2vec')
    forces = f_malloc((/1.to.3,1.to.nat,1.to.nnodes/),id='forces')
    tangent = f_malloc((/1.to.3,1.to.nat,1.to.nnodes/),id='tangent')
    energies = f_malloc((/1.to.nnodes/),id='energies')
    arc = f_malloc((/1.to.nnodes/),id='arc')

    !Copy path in correct order for 
    !spline interpolation later on.
    !Also compute energies at the nodes
    !That's ok, because their energies
    !have NOT yet been computed in optim_cg
    !(which is called by grow_string)
    emax=-huge(1._gp)
    npath=0
    do istring=1,nstring
        npath=npath+1
        do iat=1,nat
            path(npath,1,iat)=string(3*iat-2,1,istring)
            path(npath,2,iat)=string(3*iat-1,1,istring)
            path(npath,3,iat)=string(3*iat  ,1,istring)
        enddo
        !parametrize spline such that the i-th node
        !is at parameter value i:
        arc(npath)=real(npath,gp)
        !due to column major order,
        !pass string() to energy and forces, not path():
        call energyandforces(nat,alat,string(1,1,istring),&
             forces(1,1,npath),fnoise,energies(npath))
        if(energies(npath)>emax)then
            emax       = energies(npath)
            istringmax = istring
            isidemax   = 1
            ipathmax   = npath
        endif
    enddo
    do istring=nstring-ncorr,1,-1
        npath=npath+1
        do iat=1,nat
            path(npath,1,iat)=string(3*iat-2,2,istring)
            path(npath,2,iat)=string(3*iat-1,2,istring)
            path(npath,3,iat)=string(3*iat  ,2,istring)
        enddo
        !parametrize spline such that the i-th node
        !is at parameter value i:
        arc(npath)=real(npath,gp)
        !due to column major order,
        !pass string() to energy and forces, not path():
        call energyandforces(nat,alat,string(1,2,istring),&
             forces(1,1,npath),fnoise,energies(npath))
        if(energies(npath)>emax)then
            emax       = energies(npath)
            istringmax = istring
            isidemax   = 2
            ipathmax   = npath
        endif
    enddo
!DEBUGGING:
if(npath/=nnodes)stop'npat/=nnodes'
write(*,*)'ipathmax',ipathmax
    !tsguess is highest energy node:
    do iat=1,nat
        tsguess(1,iat) = path(ipathmax,1,iat)
        tsguess(2,iat) = path(ipathmax,2,iat)
        tsguess(3,iat) = path(ipathmax,3,iat)
    enddo
   
    !now we have to generate a guess for the minmode: 

    !generate spline parameters for splines used for tangents
    do i=1,nat
        call spline_wrapper(arc,path(1,1,i),npath,&
                           yp1,ypn,y2vec(1,1,i))
        call spline_wrapper(arc,path(1,2,i),npath,&
                           yp1,ypn,y2vec(1,2,i))
        call spline_wrapper(arc,path(1,3,i),npath,&
                           yp1,ypn,y2vec(1,3,i))
    enddo
    !generate tangent
    do j=1,npath
        tau=real(j,gp)
        do i=1,nat
            call splint_wrapper(arc,path(1,1,i),y2vec(1,1,i),&
                 npath,tau,rdmy,tangent(1,i,j))
            call splint_wrapper(arc,path(1,2,i),y2vec(1,2,i),&
                 npath,tau,rdmy,tangent(2,i,j))
            call splint_wrapper(arc,path(1,3,i),y2vec(1,3,i),&
                 npath,tau,rdmy,tangent(3,i,j))
        enddo
        rdmy = dnrm2(3*nat,tangent(1,1,j),1)
        tangent(:,:,j) = tangent(:,:,j) / rdmy
    enddo
    minmodeguess=tangent(:,:,ipathmax)

    if(iproc==0 .and. mhgps_verbosity>=5)then
        call write_path(nat,npath,path,energies,tangent)
    endif     

    call f_free(string) 
    call f_free(path) 
    call f_free(y2vec) 
    call f_free(forces) 
    call f_free(tangent) 
    call f_free(energies) 
    call f_free(arc) 
end subroutine
!=====================================================================
subroutine write_path(nat,npath,path,energies,tangent)
    use module_base
    use yaml_output
    use module_interfaces
    use module_global_variables, only: isadc, atoms,ixyz_int,&
                                       currDir
    implicit none
    !parameters
    integer, intent(in) :: nat, npath
    real(gp), intent(in) :: path(npath,3,nat)
    real(gp), intent(in) :: energies(npath)
    real(gp), intent(in) :: tangent(3,nat,npath)
    !internal
    integer :: ipath,iat
    character(len=4) :: fn4
    character(len=150) :: comment
    real(gp) :: pathint(3,nat,npath)

    do ipath=1,npath
        do iat=1,nat
            pathint(1,iat,ipath)=path(ipath,1,iat)
            pathint(2,iat,ipath)=path(ipath,2,iat)
            pathint(3,iat,ipath)=path(ipath,3,iat)
        enddo
    enddo

    do ipath=1,npath
        write(fn4,'(i4.4)') ipath
        write(comment,'(a)')&
           'ATTENTION! Forces below are no forces&
           but tangents to the guessed reaction path'
        call write_atomic_file(currDir//'/sad'//isadc//'_igpath_'//&
           fn4,energies(ipath),pathint(1,1,ipath),ixyz_int,&
           atoms,trim(comment),forces=tangent(1,1,ipath))
    enddo
end subroutine
!=====================================================================
subroutine grow_string(nat,alat,gammainv,perpnrmtol,trust,&
                       nstepsmax,nstringmax,nstring,string,finished)
    use module_base
    use yaml_output
    use module_global_variables, only: iproc
    implicit none
    !parameters
    integer, intent(in) :: nat
    integer, intent(in) :: nstepsmax
    integer, intent(inout) :: nstringmax
    integer, intent(inout) :: nstring
    integer, intent(out) :: finished
                          !if finished ==2: not finished
                          !if finished ==0: finished
                          !if finished ==1: one more node
                          !                 in the middle
    real(gp) , intent(in) :: gammainv
    real(gp) , intent(in) :: perpnrmtol
    real(gp) , intent(in) :: trust
    real(gp) , intent(in) :: alat(3)
    real(gp), allocatable, intent(inout) :: string(:,:,:)
    !constants
    character(len=10), parameter :: method = 'linlst'
    !internal
    integer :: i,j,k,istart
    integer, parameter :: resize=10
    real(gp), allocatable :: stringTMP(:,:,:)
    real(gp) :: tangentleft(3*nat)
    real(gp) :: tangentright(3*nat)
    real(gp) :: step
    real(gp) :: perpnrmtol_squared
    real(gp) :: trust_squared
    integer :: nresizes
    
    finished=2

    if(iproc==0)call yaml_comment('(MHGPS) entering &
                grow_string')

    if((.not. allocated(string)))then
        if(iproc==0)call yaml_warning('(MHGPS) STOP, string &
                    in grow_string not allocated')
        stop
    endif
    perpnrmtol_squared=perpnrmtol**2
    trust_squared=trust**2

    nstring=1
    step=-1._gp
    nresizes=0
    do!maximum 100 resizes of string array
        istart=nstring
        !actual string growing is done 
        !in the following loop
        do i=istart,nstringmax-1
            if(finished==0 .or. finished==1)return!interpolation done
            call interpol(method,nat,string(1,1,nstring),&
                 string(1,2,nstring),step,string(1,1,nstring+1),&
                 string(1,2,nstring+1),tangentleft,tangentright,&
                 finished)
write(*,*)'basian'
!            if(finished==1)string(1,2,nstring+1)=huge(1.0_gp)
            if(finished/=0)then
if(i/=nstring)stop'DEBUGGING i/=nstring'
               if(perpnrmtol>0)& 
                call optim_cg(nat,alat,finished,step,gammainv,&
                    perpnrmtol_squared,trust_squared,nstepsmax,&
                    tangentleft,tangentright,string(1,1,i+1),&
                    string(1,2,i+1))
                nstring=nstring+1
            endif
        enddo
        if(finished==0.or.finished==1)return
        !What follows is just resizing of string array, 
        !if needed.
        nresizes=nresizes+1
        if(nresizes>100)then
            if(iproc==0)call yaml_warning('(MHGPS) STOP, too&
                        many resizes in grow_string')
            stop
        endif
        if(allocated(stringTmp))then
            deallocate(stringTmp)
        endif
        allocate(stringTmp(3*nat,2,nstringmax))
        stringTmp=string
        deallocate(string)
        nstringmax=nstringmax+resize
        allocate(string(3*nat,2,nstringmax))
        do k=1,(nstringmax-resize)
            string(:,:,k)=stringTmp(:,:,k)
        enddo
        deallocate(stringTmp)
    enddo

end subroutine
!=====================================================================
subroutine optim_cg(nat,alat,finished,step,gammainv,&
           perpnrmtol_squared,trust_squared,nstepsmax,tangent1,&
           tangent2,rxyz1,rxyz2)
    use module_base
    use module_energyandforces
    use module_global_variables, only: iproc
    implicit none
    !parameters
    integer, intent(in)  :: nat
    integer, intent(in)  :: finished
    integer, intent(in)  :: nstepsmax
    real(gp), intent(in) :: tangent1(3*nat)
    real(gp), intent(in) :: tangent2(3*nat)
    real(gp), intent(in) :: step
    real(gp), intent(in) :: gammainv
    real(gp), intent(in) :: perpnrmtol_squared
    real(gp), intent(in) :: trust_squared
    real(gp), intent(in) :: alat(3)
    real(gp), intent(inout) :: rxyz2(3*nat)
    real(gp), intent(inout) :: rxyz1(3*nat)
    !internal
    real(gp) :: epot1
    real(gp) :: epot2
    real(gp) :: d0 !inital distance between the new nodes
    real(gp) :: fxyz1(3*nat),fxyz2(3*nat)
    real(gp) :: perp1(3*nat),perp2(3*nat)
    real(gp) :: dmax,dist,dir(3*nat)
    real(gp) :: dispPrev1(3*nat),disp1(3*nat)
    real(gp) :: dispPrev2(3*nat),disp2(3*nat)
    real(gp) :: alpha1,alpha2
    real(gp) :: perpnrmPrev1_squared, perpnrm1_squared
    real(gp) :: perpnrmPrev2_squared, perpnrm2_squared
    real(gp) :: dispnrm_squared
    integer :: istep
    real(gp) :: fnoise
    !functionals
    real(gp) :: dnrm2, ddot

    dir=rxyz2-rxyz1
    d0=dnrm2(3*nat,dir(1),1)
    dmax=d0+0.5_gp*step

    !first steps: steepest descent
    !left
    call energyandforces(nat,alat,rxyz1,fxyz1,fnoise,epot1)
    call perpend(nat,tangent1,fxyz1,perp1)
    perpnrmPrev1_squared = ddot(3*nat,perp1(1),1,perp1(1),1)
    perpnrm1_squared=perpnrmPrev1_squared
    dispPrev1=gammainv*perp1
    dispnrm_squared=maxval(dispPrev1**2)
    if(dispnrm_squared > trust_squared)then
        dispPrev1=dispPrev1*sqrt(trust_squared/dispnrm_squared)
    endif
    rxyz1=rxyz1+dispPrev1
    !right
    if(finished==2)then
        call energyandforces(nat,alat,rxyz2,fxyz2,fnoise,epot2)
        call perpend(nat,tangent2,fxyz2,perp2)
        perpnrmPrev2_squared = ddot(3*nat,perp2(1),1,perp2(1),1)
        perpnrm2_squared=perpnrmPrev2_squared
        write(*,'(a,i3.3,4(1x,es10.3))')&
        '   (MHGPS) ',1,sqrt(perpnrm1_squared),epot1,&
        sqrt(perpnrm2_squared),epot2
        dispPrev2=gammainv*perp2
        dispnrm_squared=maxval(dispPrev2**2)
        if(dispnrm_squared > trust_squared)then
            dispPrev2=dispPrev2*sqrt(trust_squared/dispnrm_squared)
        endif
        rxyz2=rxyz2+dispPrev2
    
        dir=rxyz2-rxyz1
        dist=dnrm2(3*nat,dir(1),1)
        if(dist>dmax)then
            return
        endif
    endif

    !other steps: cg
    do istep=2,nstepsmax

        !move left node
        call energyandforces(nat,alat,rxyz1,fxyz1,fnoise,epot1)
        call perpend(nat,tangent1,fxyz1,perp1)
        perpnrm1_squared = ddot(3*nat,perp1(1),1,perp1(1),1)
        if(perpnrm1_squared>perpnrmPrev1_squared)then
            alpha1=1._gp
        else
            alpha1 = perpnrm1_squared / perpnrmPrev1_squared
        endif
        disp1=gammainv*perp1+ alpha1 * dispPrev1
        dispnrm_squared=maxval(disp1**2)
        if(dispnrm_squared > trust_squared)then
             disp1=disp1*sqrt(trust_squared/dispnrm_squared)
        endif
        rxyz1=rxyz1+disp1
        dispPrev1=disp1
        perpnrmPrev1_squared=perpnrm1_squared
        
        if(finished==2)then 
            !move right node
            call energyandforces(nat,alat,rxyz2,fxyz2,fnoise,epot2)
            call perpend(nat,tangent2,fxyz2,perp2)
            perpnrm2_squared = ddot(3*nat,perp2(1),1,perp2(1),1)
            write(*,'(a,i3.3,4(1x,es10.3))')&
            '   (MHGPS) ',istep,sqrt(perpnrm1_squared),epot1,&
            sqrt(perpnrm2_squared),epot2
            if(perpnrm2_squared>perpnrmPrev2_squared)then
                alpha2=1._gp
            else
                alpha2 = perpnrm2_squared / perpnrmPrev2_squared
            endif
            disp2=gammainv*perp2+ alpha2 * dispPrev2
            dispnrm_squared=maxval(disp2**2)
            if(dispnrm_squared > trust_squared)then
                 disp2=disp2*sqrt(trust_squared/dispnrm_squared)
            endif
            rxyz2=rxyz2+disp2
            dispPrev2=disp2
            perpnrmPrev2_squared=perpnrm2_squared
        
            dir=rxyz2-rxyz1
            dist=dnrm2(3*nat,dir(1),1)
            !perpnrm1_squared is from last iteration
            !but we accept this, since reevaluation
            !of energies and forces is too expensive
            !just for the purpose of the following
            !comparison of perpnrm1_squared and perpnrmtol_squared
            if(dist>dmax.or. (perpnrm1_squared<perpnrmtol_squared&
            !if((perpnrm1_squared<perpnrmtol_squared &
               &.and. perpnrm2_squared<perpnrmtol_squared))then
            !!!if(dist>dmax)then
            !!!!            if(dist>dmax.and.iproc==0)then
            !!!                write(200,*)'exit due to dmax'   
            !!!            endif
                !we do not compute and do not return energies and
                !forces of the latest rxyz2 and rxyz2. If needed,
                !the must be computed outside.
            !!!            return
            endif
        else if (perpnrm1_squared<perpnrmtol_squared) then
            return
        endif
    enddo
end subroutine
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
    integer, intent(in)    :: nat
    real(gp), intent(in)   :: left(3*nat)
    real(gp), intent(in)   :: right(3*nat)
    real(gp), intent(inout):: step
    real(gp), intent(out)  :: interleft(3*nat)
    real(gp), intent(out)  :: interright(3*nat)
    real(gp), intent(out)  :: tangent(3*nat)
    integer, intent(out)   :: finished
    !constants
    real(gp), parameter :: stepfrct=0.1_gp! freezing string step size
    !internal
    real(gp) :: arcl
    !functions
    real(gp) :: dnrm2
stop 'lin_interpol not debugged yet'
    !tangent points from left to right:    
    tangent = right-left
    arcl = dnrm2(3*nat,tangent(1),1)
    tangent = tangent / arcl
    
    if(step<0._gp)step=stepfrct * arcl

    if(arcl < step)then
        finished=0
    else if(arcl < 2._gp*step)then!only one more point
        interleft  = left + 0.5_gp * arcl * tangent
        finished = 1
    else
        interleft  = left + step * tangent
        interright = right - step * tangent
        finished = 2
    endif
end subroutine
!=====================================================================
subroutine lst_interpol(nat,left,right,step,interleft,interright,&
                        tangentleft,tangentright,finished)
    !Given two distinct structures, lst_interpol interpolates
    !inwards (that is in a direction connecting both strucutres)
    !using the linear synchronous transit (LST) technique.
    !A high density path made from 'nimages' nodes using LST
    !is generated. Then this path is parametrized as a function
    !of its integrated path length using natural cubic splines.
    !In order to avoid uncontinous changes of the tangent direction,
    !a second spline parametrization is done by using
    !nimagesC<<nimages equally spaced nodes from the first spline
    !parameterization. Tangent directions are taken from this 
    !second spline.
    !
    !on return:
    !if finished =2: interleft and intergiht contain the new nodes
    !if finished =1 0: left and right are too close, only one new node
    !                 is returned in interleft.
    !                 interright is meaningless in this case.
    !if finished = 0: left and right are closer than 'step'
    !                 => freezing string search finsihed
    !                 nothing is returned, interleft and interright
    !                 are meaningless
    use module_base
    use module_interpol
    use module_global_variables, only:stepfrct=>lst_interpol_stepfrct
    implicit none
    !parameters
    integer, intent(in)      :: nat
    real(gp), intent(in)     :: left(3,nat)
    real(gp), intent(in)     :: right(3,nat)
    real(gp), intent(inout)  :: step
    real(gp), intent(out)    :: interleft(3,nat)
    real(gp), intent(out)    :: interright(3,nat)
    real(gp), intent(out)    :: tangentleft(3,nat)
    real(gp), intent(out)    :: tangentright(3,nat)
    integer, intent(out)     :: finished
    !constants
    integer, parameter  :: nimages=200 
    integer, parameter  :: nimagesC=5 !setting nimagesC=nimages
                                      !should give similar
                                      !implementation to the freezing
                                      !string publication
!    real(gp), parameter :: stepfrct=0.1_gp! freezing string step size
    !real(gp), parameter :: stepfrct=0.05555555556_gp! freezing string
                                                     ! step size
    !internal
    integer  :: i
    integer  :: j
    integer  :: tnat
    integer  :: iat
    integer  :: nimagestang
    real(gp) :: lstpath(3,nat,nimages)
    real(gp) :: lstpathRM(nimages,3,nat)
    real(gp) :: lstpathCRM(nimagesC,3,nat)
    real(gp) :: arc(nimages)
    real(gp) :: arcl
    real(gp) :: arcC(nimagesC)
    real(gp) :: arclC
    real(gp) :: diff(3,nat)
    real(gp) :: nimo
    real(gp) :: yp1=huge(1._gp), ypn=huge(1._gp)!natural splines
    real(gp) :: y2vec(nimages,3,nat)
    real(gp) :: y2vecC(nimagesC,3,nat)
    real(gp) :: tau
    real(gp) :: rdmy
    real(gp) :: lambda
    !functions
    real(gp) :: dnrm2

!<-DEBUG START------------------------------------------------------>
!character(len=5) :: fc5
!character(len=200) :: filename,line
!integer :: istat
!integer, save :: ic
!real(gp) :: dmy
!character(len=5):: xat(22)
!open(unit=33,file='input001/pos001.ascii')
!read(33,*)
!read(33,*)
!read(33,*)
!read(33,*)
!do iat=1,22
!    read(33,'(a)',iostat=istat)line
!    if(istat/=0)exit
!    read(line,*)dmy,dmy,dmy,xat(iat)
!enddo
!close(33)
!<-DEBUG END-------------------------------------------------------->

    tnat=3*nat

    !create high density lst path
    nimo=1._gp/real(nimages-1,gp)
    do i=1,nimages
        lambda  = real(i-1,gp)*nimo
        call lstpthpnt(nat,left,right,lambda,lstpath(1,1,i))
    enddo

    !measure arc length 
    arc(1)=0._gp
    do i=2,nimages
        diff = lstpath(:,:,i) - lstpath(:,:,i-1)
        arc(i)  = arc(i-1) + dnrm2(tnat,diff(1,1),1)
    enddo
    arcl=arc(nimages)

    if(step<0._gp)step=stepfrct*arcl
write(*,*)'arcl ',arcl,step
    if(arcl < step)then
write(*,*)'lstinterpol finished 0'
        finished=0
        return
    endif
write(*,*)'lstinterpol 1'

    !rewrite lstpath to row major ordering
    !(for faster access in spline routines)
    do iat=1,nat
        do i=1,nimages
            lstpathRM(i,1,iat)=lstpath(1,iat,i)
            lstpathRM(i,2,iat)=lstpath(2,iat,i)
            lstpathRM(i,3,iat)=lstpath(3,iat,i)
        enddo
    enddo

    !compute the spline parameters (y2vec)
    !parametrize curve as a function of the
    !integrated arc length
    do i=1,nat
        call spline_wrapper(arc,lstpathRM(1,1,i),nimages,&
                           yp1,ypn,y2vec(1,1,i))
        call spline_wrapper(arc,lstpathRM(1,2,i),nimages,&
                           yp1,ypn,y2vec(1,2,i))
        call spline_wrapper(arc,lstpathRM(1,3,i),nimages,&
                           yp1,ypn,y2vec(1,3,i))
    enddo

    !generate nodes at which tangents are computed
    nimagestang=min(nimagesC,nimages)
    nimo=1._gp/real(nimagestang-1,gp)
    do j=1,nimagestang
        tau  = arcl*real(j-1,gp)*nimo
        arcC(j)=tau
        do i=1,nat
            call splint_wrapper(arc,lstpathRM(1,1,i),y2vec(1,1,i),&
                 nimages,tau,lstpathCRM(j,1,i),rdmy)
            call splint_wrapper(arc,lstpathRM(1,2,i),y2vec(1,2,i),&
                 nimages,tau,lstpathCRM(j,2,i),rdmy)
            call splint_wrapper(arc,lstpathRM(1,3,i),y2vec(1,3,i),&
                 nimages,tau,lstpathCRM(j,3,i),rdmy)
        enddo
    enddo
   
    !generate spline parameters for splines used for tangents
    do i=1,nat
        call spline_wrapper(arcC,lstpathCRM(1,1,i),nimagestang,&
                           yp1,ypn,y2vecC(1,1,i))
        call spline_wrapper(arcC,lstpathCRM(1,2,i),nimagestang,&
                           yp1,ypn,y2vecC(1,2,i))
        call spline_wrapper(arcC,lstpathCRM(1,3,i),nimagestang,&
                           yp1,ypn,y2vecC(1,3,i))
    enddo
write(*,*)'lstinterpol 2'

!<-DEBUG START------------------------------------------------------>
!!check interpolated path
!do j=1,200
!tau  = arcl*real(j-1,gp)/real(200-1,gp)
!!tau  = arc(j)
!        do i=1,nat
!            call splint_wrapper(arc,lstpathRM(1,1,i),y2vec(1,1,i),&
!                 nimages,tau,interleft(1,i),rdmy)
!            call splint_wrapper(arc,lstpathRM(1,2,i),y2vec(1,2,i),&
!                 nimages,tau,interleft(2,i),rdmy)
!            call splint_wrapper(arc,lstpathRM(1,3,i),y2vec(1,3,i),&
!                 nimages,tau,interleft(3,i),rdmy)
!        enddo
!        do i=1,nat
!            call splint_wrapper(arcC,lstpathCRM(1,1,i),&
!                 y2vecC(1,1,i),nimagestang,tau,rdmy,tangentleft(1,i))
!write(*,*)rdmy-interleft(1,i)
!!write(*,*)y2vecC(:,1,i)-y2vec(:,1,i)
!!write(*,*)lstpathCRM(:,1,i)-lstpathRM(:,1,i)
!            call splint_wrapper(arcC,lstpathCRM(1,2,i),&
!                 y2vecC(1,2,i),nimagestang,tau,rdmy,tangentleft(2,i))
!write(*,*)rdmy-interleft(2,i)
!!write(*,*)y2vecC(:,2,i)-y2vec(:,2,i)
!!write(*,*)lstpathCRM(:,2,i)-lstpathRM(:,2,i)
!            call splint_wrapper(arcC,lstpathCRM(1,3,i),&
!                 y2vecC(1,3,i),nimagestang,tau,rdmy,tangentleft(3,i))
!write(*,*)rdmy-interleft(3,i)
!!write(*,*)y2vecC(:,3,i)-y2vec(:,3,i)
!!write(*,*)lstpathCRM(:,3,i)-lstpathRM(:,3,i)
!        enddo
!!if(mod(j,100)==0)then
!write(fc5,'(i5.5)')j
!write(filename,*)'pospline_'//fc5//'.ascii'
!open(99,file=trim(adjustl((filename))))
!write(99,'(a)')'# BigDFT file'
!write(99,*)10.0 ,0, 10.0 
!write(99,*)0, 0, 10.0 
!do iat=1,nat
!write(99,'(3(1xes24.17),1x,a)')interleft(1,iat)*0.529d0,interleft(2,iat)&
!                           *0.529d0,interleft(3,iat)*0.529d0,xat(iat)
!enddo
!write(99,'(a)')"#metaData: forces (Ha/Bohr) =[ \"
!do iat=1,nat-1
!write(99,'(a,3(1x,es24.17";"),1x,a)')'#',tangentleft(1,iat)*0.529d0,&
!            tangentleft(2,iat)*0.529d0,tangentleft(3,iat)*0.529d0,' \'
!enddo
!iat=nat
!write(99,'(a,3(1x,es24.17";"),1x,a)')'#',tangentleft(1,iat)*0.529d0,&
!            tangentleft(2,iat)*0.529d0,tangentleft(3,iat)*0.529d0,' ]'
!close(99)
!!endif
!
!
!enddo
!stop
!<-DEBUG END-------------------------------------------------------->


    if(arcl < 2._gp*step)then!only one more point
write(*,*)'lstinterpol 3'
        !we have to return the point in the 'middle'    
        tau = 0.5_gp*arcl
        !generate coordinates
        do i=1,nat
            call splint_wrapper(arc,lstpathRM(1,1,i),y2vec(1,1,i),&
                 nimages,tau,interleft(1,i),rdmy)
            call splint_wrapper(arc,lstpathRM(1,2,i),y2vec(1,2,i),&
                 nimages,tau,interleft(2,i),rdmy)
            call splint_wrapper(arc,lstpathRM(1,3,i),y2vec(1,3,i),&
                 nimages,tau,interleft(3,i),rdmy)
        enddo
        !generate tangent
        do i=1,nat
            call splint_wrapper(arcC,lstpathCRM(1,1,i),y2vecC(1,1,i),&
                 nimagestang,tau,rdmy,tangentleft(1,i))
            call splint_wrapper(arcC,lstpathCRM(1,2,i),y2vecC(1,2,i),&
                 nimagestang,tau,rdmy,tangentleft(2,i))
            call splint_wrapper(arcC,lstpathCRM(1,3,i),y2vecC(1,3,i),&
                 nimagestang,tau,rdmy,tangentleft(3,i))
        enddo
        rdmy = dnrm2(tnat,tangentleft(1,1),1)
        tangentleft = tangentleft / rdmy
        !return code: only one more node inserted
        finished=1
    else! standard case
write(*,*)'lstinterpol 4'
        !we have to return the two points interleft
        !and interright whose distances to left and right
        !are roughly 'step'

        !first left...
        tau = step
        do i=1,nat
            !potentially performance issues since lstpath
            !is not transversed in column-major order in
            !splint_wrapper
            call splint_wrapper(arc,lstpathRM(1,1,i),y2vec(1,1,i),&
                 nimages,tau,interleft(1,i),tangentleft(1,i))
            call splint_wrapper(arc,lstpathRM(1,2,i),y2vec(1,2,i),&
                 nimages,tau,interleft(2,i),tangentleft(2,i))
            call splint_wrapper(arc,lstpathRM(1,3,i),y2vec(1,3,i),&
                 nimages,tau,interleft(3,i),tangentleft(3,i))
        enddo
        !generate coordinates for left node
        do i=1,nat
            call splint_wrapper(arc,lstpathRM(1,1,i),y2vec(1,1,i),&
                 nimages,tau,interleft(1,i),rdmy)
            call splint_wrapper(arc,lstpathRM(1,2,i),y2vec(1,2,i),&
                 nimages,tau,interleft(2,i),rdmy)
            call splint_wrapper(arc,lstpathRM(1,3,i),y2vec(1,3,i),&
                 nimages,tau,interleft(3,i),rdmy)
        enddo
        !generate tangent for left node
        do i=1,nat
            call splint_wrapper(arcC,lstpathCRM(1,1,i),y2vecC(1,1,i),&
                 nimagestang,tau,rdmy,tangentleft(1,i))
            call splint_wrapper(arcC,lstpathCRM(1,2,i),y2vecC(1,2,i),&
                 nimagestang,tau,rdmy,tangentleft(2,i))
            call splint_wrapper(arcC,lstpathCRM(1,3,i),y2vecC(1,3,i),&
                 nimagestang,tau,rdmy,tangentleft(3,i))
        enddo
        rdmy = dnrm2(tnat,tangentleft(1,1),1)
        tangentleft = tangentleft / rdmy

        !...then right
        tau = arcl-step
        !generate coordinates for right node
        do i=1,nat
            call splint_wrapper(arc,lstpathRM(1,1,i),y2vec(1,1,i),&
                 nimages,tau,interright(1,i),rdmy)
            call splint_wrapper(arc,lstpathRM(1,2,i),y2vec(1,2,i),&
                 nimages,tau,interright(2,i),rdmy)
            call splint_wrapper(arc,lstpathRM(1,3,i),y2vec(1,3,i),&
                 nimages,tau,interright(3,i),rdmy)
        enddo
        !generate tangent for right node
        do i=1,nat
            call splint_wrapper(arcC,lstpathCRM(1,1,i),y2vecC(1,1,i),&
                 nimagestang,tau,rdmy,tangentright(1,i))
            call splint_wrapper(arcC,lstpathCRM(1,2,i),y2vecC(1,2,i),&
                 nimagestang,tau,rdmy,tangentright(2,i))
            call splint_wrapper(arcC,lstpathCRM(1,3,i),y2vecC(1,3,i),&
                 nimagestang,tau,rdmy,tangentright(3,i))
        enddo
        rdmy = dnrm2(tnat,tangentright(1,1),1)
        tangentright = tangentright / rdmy

        !return code: two more nodes inserted
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
write(*,*)'interpol'
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
    use module_misc
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
end module
