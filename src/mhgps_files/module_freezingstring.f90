module module_freezingstring
    implicit none
    
    contains
subroutine grow_string(nat,alat,step,gammainv,perpnrmtol,trust,nstepsmax,nstringmax,nstring,string)
    implicit none
    !parameters
    integer, intent(in) :: nat,nstepsmax
    integer, intent(inout) :: nstringmax,nstring
    real(8) , intent(in) :: step,gammainv,perpnrmtol,trust,alat(3)
    real(8), allocatable, intent(inout) :: string(:,:,:)
    !internal
    integer :: i,j,k,istart
    integer, parameter :: resize=10
    real(8), allocatable :: stringTMP(:,:,:)
    real(8) :: tangent(3*nat)
    real(8) :: dist,dir(3*nat),stepint,perpnrmtol_squared,trust_squared
    
    real(8) :: dist1,dist2

    if(.not. allocated(string))then
        write(*,*)'STOP, string not allocated'
    endif
    perpnrmtol_squared=perpnrmtol**2
    trust_squared=trust**2

!write(200,*)'#######################'

    nstring=1
    dir = string(:,2,1)-string(:,1,1)
    dist = sqrt(dot_product(dir,dir))
  stepint=0.1d0*dist
!  stepint=step
!write(333,*)dist
    if(dist<=2.d0*stepint)then!interpolation done
        return
    endif
    do j=1,100!maximum 100 resizes of string array
        istart=nstring
        do i=istart,nstringmax-1
!            dir = string(:,2,i)-string(:,1,i)
!            dist = sqrt(dot_product(dir,dir))
!            stepint=0.05d0*dist
            call lin_interpol(nat,string(1,1,nstring),string(1,2,nstring),stepint,string(1,1,nstring+1),string(1,2,nstring+1),tangent)
            dir = string(:,1,nstring+1)-string(:,2,nstring)
            dist1 = sqrt(dot_product(dir,dir))
            dir = string(:,1,nstring+1)-string(:,1,nstring)
            dist2 = sqrt(dot_product(dir,dir))
            if(dist1<dist2)then
            dir = string(:,2,nstring)-string(:,1,nstring)
            dist = sqrt(dot_product(dir,dir))
            write(*,*)dist,dist1,dist2,j,i,stepint,nstring
                stop'dist1<dist2'
            endif
            call optim_cg(nat,alat,stepint,gammainv,perpnrmtol_squared,trust_squared,nstepsmax,tangent,string(1,1,i+1),string(1,2,i+1))
            dir = string(:,2,i+1)-string(:,1,i+1)
            dist = sqrt(dot_product(dir,dir))
            nstring=nstring+1
!write(*,*)'nstring',nstring,nstringmax-1
!write(*,*)'dist',dist
            if(dist<=2.d0*stepint)then!interpolation done
                return
            endif
        enddo
        if(allocated(stringTmp))then
            deallocate(stringTmp)
        endif
        allocate(stringTmp(3*nat,2,nstringmax))
        stringTmp=string
        deallocate(string)
        nstringmax=nstringmax+resize
        allocate(string(3*nat,2,nstringmax))
!write(*,*)'resize',j
        do k=1,(nstringmax-resize)
            string(:,:,k)=stringTmp(:,:,k)
        enddo
        deallocate(stringTmp)
    enddo

end subroutine
subroutine optim_cg(nat,alat,step,gammainv,perpnrmtol_squared,trust_squared,nstepsmax,tangent,rxyz1,rxyz2)
    implicit none
    !parameters
    integer, intent(in) :: nat,nstepsmax
    real(8), intent(in) :: tangent(3*nat),step,gammainv,perpnrmtol_squared,trust_squared,alat(3)
    real(8), intent(inout) :: rxyz1(3*nat),rxyz2(3*nat)
    !internal
    real(8) :: fxyz1(3*nat),fxyz2(3*nat)
    real(8) :: perp1(3*nat),perp2(3*nat)
    real(8) :: epot1, epot2
    real(8) :: d0,dmax,dist,dir(3*nat)
    real(8) :: dispPrev1(3*nat),disp1(3*nat)
    real(8) :: dispPrev2(3*nat),disp2(3*nat)
    real(8) :: alpha1,alpha2
    real(8) :: fnrmPrev1, fnrm1
    real(8) :: fnrmPrev2, fnrm2
    real(8) :: perpnrmPrev1_squared, perpnrm1_squared
    real(8) :: perpnrmPrev2_squared, perpnrm2_squared
    real(8) :: dispnrm_squared
    integer :: istep

    real(8) :: tangentint(3*nat),tmp

    tangentint=tangent


    dir=rxyz2-rxyz1
    d0=sqrt(dot_product(dir,dir))
    dmax=d0+0.5d0*step
    call energyandforces(nat, alat, rxyz1,fxyz1,epot1,'cnt_enf_freezing_string')
    call energyandforces(nat, alat, rxyz2,fxyz2,epot2,'cnt_enf_freezing_string')

    fnrmPrev1=sqrt(dot_product(fxyz1,fxyz1))
    fnrmPrev2=sqrt(dot_product(fxyz2,fxyz2))

    call perpend(nat,tangentint,fxyz1,perp1)
    call perpend(nat,tangentint,fxyz2,perp2)

    perpnrmPrev1_squared = dot_product(perp1,perp1)
    perpnrmPrev2_squared = dot_product(perp2,perp2)
    perpnrm1_squared=perpnrmPrev1_squared
    perpnrm2_squared=perpnrmPrev2_squared

    !first steps: steepest descent
    dispPrev1=gammainv*perp1
    dispnrm_squared=maxval(dispPrev1**2)
    if(dispnrm_squared > trust_squared)then
        dispPrev1=dispPrev1*sqrt(trust_squared/dispnrm_squared)
    endif
    rxyz1=rxyz1+dispPrev1
    dispPrev2=gammainv*perp2
    dispnrm_squared=maxval(dispPrev2**2)
    if(dispnrm_squared > trust_squared)then
        dispPrev2=dispPrev2*sqrt(trust_squared/dispnrm_squared)
    endif
    rxyz2=rxyz2+dispPrev2

    call energyandforces(nat, alat, rxyz1,fxyz1,epot1,'cnt_enf_freezing_string')
    call energyandforces(nat, alat, rxyz2,fxyz2,epot2,'cnt_enf_freezing_string')
    fnrm1=sqrt(dot_product(fxyz1,fxyz1))
    fnrm2=sqrt(dot_product(fxyz2,fxyz2))
    
    dir=rxyz2-rxyz1
    dist=sqrt(dot_product(dir,dir))
    if(dist>dmax)then
        return
    endif

!write(200,*)'-------------------'
    !other steps: cg
    do istep=2,nstepsmax

    dir=rxyz2-rxyz1
    dist=sqrt(dot_product(dir,dir))
    tangentint=dir/dist
        !move left node
        call perpend(nat,tangentint,fxyz1,perp1)
        perpnrm1_squared = dot_product(perp1,perp1)
!        if(perpnrm1_squared>=perpnrmtol_squared)then
!!            if(fnrm1>fnrmPrev1)then
            if(perpnrm1_squared>perpnrmPrev1_squared)then
                alpha1=1.d0
            else
                alpha1 = perpnrm1_squared / perpnrmPrev1_squared
            endif
            disp1=gammainv*perp1+ alpha1 * dispPrev1
            dispnrm_squared=maxval(disp1**2)
!write(*,*)'trust1',sqrt(dispnrm_squared),sqrt(trust_squared)
            if(dispnrm_squared > trust_squared)then
!write(200,*)'Limiting step length!'
!                disp1 = sqrt(trust_squared) * disp1/sqrt(dot_product(disp1,disp1))
                 disp1=disp1*sqrt(trust_squared/dispnrm_squared)
            endif
            rxyz1=rxyz1+disp1
            dispPrev1=disp1
            perpnrmPrev1_squared=perpnrm1_squared
            fnrmPrev1=fnrm1
            call energyandforces(nat, alat, rxyz1,fxyz1,epot1,'cnt_enf_freezing_string')
            fnrm1=sqrt(dot_product(fxyz1,fxyz1))
!        endif

        !move right node
            call energyandforces(nat, alat, rxyz2,fxyz2,epot2,'cnt_enf_freezing_string')
            fnrm2=sqrt(dot_product(fxyz2,fxyz2))
            call perpend(nat,tangentint,fxyz2,perp2)
            perpnrm2_squared = dot_product(perp2,perp2)
!        if(perpnrm2_squared>=perpnrmtol_squared)then
!            if(fnrm2>fnrmPrev2)then
            if(perpnrm2_squared>perpnrmPrev2_squared)then
                alpha2=1.d0
            else
                alpha2 = perpnrm2_squared / perpnrmPrev2_squared
            endif
            disp2=gammainv*perp2+ alpha2 * dispPrev2
            dispnrm_squared=maxval(disp2**2)
!write(*,*)'trust2',sqrt(dispnrm_squared),sqrt(trust_squared)
            if(dispnrm_squared > trust_squared)then
!write(200,*)'Limiting step length!'
!                disp2 = sqrt(trust_squared) * disp2/sqrt(dot_product(disp2,disp2))
                 disp2=disp2*sqrt(trust_squared/dispnrm_squared)
            endif
!write(*,*)'trust2',sqrt(dispnrm_squared),sqrt(trust_squared)
!write(*,*)'disp2',sqrt(dot_product(disp2,disp2))
            rxyz2=rxyz2+disp2
            dispPrev2=disp2
            perpnrmPrev2_squared=perpnrm2_squared
            fnrmPrev2=fnrm2
!        endif
    
!write(*,*)alpha1,alpha2
!write(200,*)'perp',sqrt(perpnrm1_squared), sqrt(perpnrm2_squared)
!write(*,*)tangent
!write(*,*)rxyz2
!write(*,*)perpnrm2_squared
        dir=rxyz2-rxyz1
        dist=sqrt(dot_product(dir,dir))
        if(dist>dmax.or. (perpnrm1_squared<perpnrmtol_squared &
        !if((perpnrm1_squared<perpnrmtol_squared &
                    &.and. perpnrm2_squared<perpnrmtol_squared))then
            if(dist>dmax)then
                 write(200,*)'exit due to dmax'   
            endif
            return
        endif


    enddo
end subroutine
subroutine perpend(nat,tangent,fxyz,perp)
    !returns a vector perp that contains
    !the perpendicular components of fyxz to
    !the tangent vector
    !that is: all components of fxyz in direction
    !of tangent are substracted from xyz
    implicit none
    !parameters
    integer, intent(in) :: nat
    real(8), intent(in) :: tangent(3*nat),fxyz(3*nat)
    real(8), intent(out) :: perp(3*nat)

    perp = fxyz - dot_product(tangent,fxyz)*tangent
    
end subroutine
subroutine lin_interpol(nat,left, right, step,interleft,interright,tangent)
    implicit none
    !parameters
    integer, intent(in) :: nat
    real(8), intent(in) :: left(3*nat), right(3*nat),step
    real(8), intent(out) :: interleft(3*nat), interright(3*nat), tangent(3*nat)
    !internal
    
    tangent = right-left
    tangent = tangent / sqrt(dot_product(tangent,tangent))
    interleft = left+step*tangent
    interright = right - step*tangent
end subroutine
subroutine lst_interpol()
!HIER WEITER
HIER WEITER
end module

