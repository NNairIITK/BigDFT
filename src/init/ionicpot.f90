subroutine IonicEnergyandForces(iproc,nproc,at,hxh,hyh,hzh,&
     rxyz,eion,fion,psoffset,n1,n2,n3,n1i,n2i,n3i,i3s,n3pi,pot_ion,pkernel)
  use module_base
  use module_types
  use Poisson_Solver
  implicit none
  type(atoms_data), intent(in) :: at
  integer, intent(in) :: iproc,nproc,n1,n2,n3,n1i,n2i,n3i,i3s,n3pi
  real(gp), intent(in) :: hxh,hyh,hzh
  real(gp), dimension(3,at%nat), intent(in) :: rxyz
  real(dp), dimension(*), intent(in) :: pkernel
  real(gp), intent(out) :: eion,psoffset
  real(gp), dimension(3,at%nat), intent(out) :: fion
  real(dp), dimension(*), intent(out) :: pot_ion
  !local variables
  character(len=*), parameter :: subname='IonicEnergyandForces'
  logical :: slowion=.false.
  logical :: perx,pery,perz,gox,goy,goz
  integer :: iat,ii,i_all,i_stat,ityp,jat,jtyp,nbl1,nbr1,nbl2,nbr2,nbl3,nbr3
  integer :: isx,iex,isy,iey,isz,iez,i1,i2,i3,j1,j2,j3,ind,ierr
  real(gp) :: ucvol,rloc,twopitothreehalf,pi,atint,shortlength,charge,eself,rx,ry,rz
  real(gp) :: fxion,fyion,fzion,dist,fxslf,fyslf,fzslf,fxerf,fyerf,fzerf,cutoff,zero
  real(gp) :: hxx,hxy,hxz,hyy,hyz,hzz,chgprod
  real(gp) :: x,y,z,xp,Vel,prefactor,r2,arg,ehart,Mz,cmassy
  real(gp), dimension(3,3) :: gmet,rmet,rprimd,gprimd
  !other arrays for the ewald treatment
  real(gp), dimension(:,:), allocatable :: fewald,xred,gion

  pi=4.d0*datan(1.d0)

  if (at%geocode == 'P') then
     !here we insert the calculation of the ewald forces
     allocate(fewald(3,at%nat+ndebug),stat=i_stat)
     call memocc(i_stat,fewald,'fewald',subname)
     allocate(xred(3,at%nat+ndebug),stat=i_stat)
     call memocc(i_stat,xred,'xred',subname)

     !calculate rprimd
     rprimd(:,:)=0.0_gp

     rprimd(1,1)=at%alat1
     rprimd(2,2)=at%alat2
     rprimd(3,3)=at%alat3

     !calculate the metrics and the volume
     call metric(gmet,gprimd,-1,rmet,rprimd,ucvol)

     !calculate reduced coordinates
     do iat=1,at%nat
        do ii=1,3
           xred(ii,iat)= gprimd(1,ii)*rxyz(1,iat)+gprimd(2,ii)*rxyz(2,iat)+&
                gprimd(3,ii)*rxyz(3,iat)
        end do
     end do

     !calculate ewald energy and forces
     call ewald(eion,gmet,fewald,at%nat,at%ntypes,rmet,at%iatype,ucvol,&
          xred,real(at%nelpsp,kind=8))

     !make forces dimensional
     do iat=1,at%nat
        do ii=1,3
           fion(ii,iat)= - (gprimd(ii,1)*fewald(1,iat)+&
                gprimd(ii,2)*fewald(2,iat)+&
                gprimd(ii,3)*fewald(3,iat))
        end do
        !if (nproc==1 .and. slowion) print *,'iat,fion',iat,(fion(j1,iat),j1=1,3)
     end do

     i_all=-product(shape(xred))*kind(xred)
     deallocate(xred,stat=i_stat)
     call memocc(i_stat,i_all,'xred',subname)
     i_all=-product(shape(fewald))*kind(fewald)
     deallocate(fewald,stat=i_stat)
     call memocc(i_stat,i_all,'fewald',subname)

     !now calculate the integral of the local psp
     !this is the offset to be applied in the Poisson Solver to have a neutralizing background
     psoffset=0.0_gp
     shortlength=0.0_gp
     charge=0.0_gp
     twopitothreehalf=2.0_gp*pi*sqrt(2.0_gp*pi)
     do iat=1,at%nat
        ityp=at%iatype(iat)
        rloc=at%psppar(0,0,ityp)
        atint=at%psppar(0,1,ityp)+3.0_gp*at%psppar(0,2,ityp)+&
             15.0_gp*at%psppar(0,3,ityp)+105.0_gp*at%psppar(0,4,ityp)
        psoffset=psoffset+rloc**3*atint
        shortlength=shortlength+real(at%nelpsp(ityp),gp)*rloc**2
        charge=charge+real(at%nelpsp(ityp),gp)
     end do
     psoffset=twopitothreehalf*psoffset
     shortlength=shortlength*2.d0*pi

     !print *,'psoffset',psoffset,'pspcore',(psoffset+shortlength)*charge/(at%alat1*at%alat2*at%alat3)
     !if (iproc ==0) print *,'eion',eion,charge/ucvol*(psoffset+shortlength)
     !correct ionic energy taking into account the PSP core correction
     eion=eion+charge/ucvol*(psoffset+shortlength)

!!$     !in the surfaces case, correct the energy term following (J.Chem.Phys. 111(7)-3155, 1999)
!!$     if (at%geocode == 'S') then
!!$        !calculate the Mz dipole component (which in our case corresponds to y direction)
!!$        !first calculate the center of mass
!!$        cmassy=0.0_gp
!!$        do iat=1,at%nat
!!$           cmassy=cmassy+rxyz(2,iat)
!!$        end do
!!$        
!!$        Mz=0.0_gp
!!$        do iat=1,at%nat
!!$           ityp=at%iatype(iat)
!!$           Mz=Mz+real(at%nelpsp(ityp),gp)*(rxyz(2,iat)-cmassy)
!!$        end do
!!$        
!!$        !correct energy and forces in the y direction
!!$        eion=eion+0.5_gp/ucvol*Mz**2
!!$        do iat=1,at%nat
!!$           ityp=at%iatype(iat)
!!$           fion(2,iat)=fion(2,iat)-real(at%nelpsp(ityp),gp)/ucvol*Mz
!!$           if (nproc==1 .and. slowion) print *,'iat,fion',iat,(fion(j1,iat),j1=1,3)
!!$        end do
!!$
!!$     end if

  else if (at%geocode == 'F') then

     eion=0.0_gp
     do iat=1,at%nat
        ityp=at%iatype(iat)
        rx=rxyz(1,iat) 
        ry=rxyz(2,iat)
        rz=rxyz(3,iat)
        !inizialization of the forces
        fxion=0.0_gp
        fyion=0.0_gp
        fzion=0.0_gp
        !initialisation of the hessian
        hxx=0.0_gp
        hxy=0.0_gp
        hxz=0.0_gp
        hyy=0.0_gp
        hyz=0.0_gp
        hzz=0.0_gp

        !    ion-ion interaction
        do jat=1,iat-1
           dist=sqrt((rx-rxyz(1,jat))**2+(ry-rxyz(2,jat))**2+(rz-rxyz(3,jat))**2)
           jtyp=at%iatype(jat)
           chgprod=real(at%nelpsp(jtyp),gp)*real(at%nelpsp(ityp),gp)
           eion=eion+chgprod/dist
           !forces
           fxion=fxion+chgprod/(dist**3)*(rx-rxyz(1,jat))
           fyion=fyion+chgprod/(dist**3)*(ry-rxyz(2,jat))
           fzion=fzion+chgprod/(dist**3)*(rz-rxyz(3,jat))
           !hessian matrix
           hxx=hxx+3.0_gp*chgprod/(dist**5)*(rx-rxyz(1,jat))**2-chgprod/(dist**3)
           hxy=hxy+3.0_gp*chgprod/(dist**5)*(rx-rxyz(1,jat))*(ry-rxyz(2,jat))
           hxz=hxz+3.0_gp*chgprod/(dist**5)*(rx-rxyz(1,jat))*(rz-rxyz(3,jat))
           hyy=hyy+3.0_gp*chgprod/(dist**5)*(ry-rxyz(2,jat))**2-chgprod/(dist**3)
           hyz=hyz+3.0_gp*chgprod/(dist**5)*(ry-rxyz(2,jat))*(rz-rxyz(3,jat))
           hzz=hzz+3.0_gp*chgprod/(dist**5)*(rz-rxyz(3,jat))**2-chgprod/(dist**3)
        enddo
        do jat=iat+1,at%nat
           dist=sqrt((rx-rxyz(1,jat))**2+(ry-rxyz(2,jat))**2+(rz-rxyz(3,jat))**2)
           jtyp=at%iatype(jat)
           chgprod=real(at%nelpsp(jtyp),gp)*real(at%nelpsp(ityp),gp)
           !forces
           fxion=fxion+chgprod/(dist**3)*(rx-rxyz(1,jat))
           fyion=fyion+chgprod/(dist**3)*(ry-rxyz(2,jat))
           fzion=fzion+chgprod/(dist**3)*(rz-rxyz(3,jat))
           !hessian matrix
           hxx=hxx+3.0_gp*chgprod/(dist**5)*(rx-rxyz(1,jat))**2-chgprod/(dist**3)
           hxy=hxy+3.0_gp*chgprod/(dist**5)*(rx-rxyz(1,jat))*(ry-rxyz(2,jat))
           hxz=hxz+3.0_gp*chgprod/(dist**5)*(rx-rxyz(1,jat))*(rz-rxyz(3,jat))
           hyy=hyy+3.0_gp*chgprod/(dist**5)*(ry-rxyz(2,jat))**2-chgprod/(dist**3)
           hyz=hyz+3.0_gp*chgprod/(dist**5)*(ry-rxyz(2,jat))*(rz-rxyz(3,jat))
           hzz=hzz+3.0_gp*chgprod/(dist**5)*(rz-rxyz(3,jat))**2-chgprod/(dist**3)
        end do

        fion(1,iat)=fxion
        fion(2,iat)=fyion
        fion(3,iat)=fzion

        !if (nproc==1 .and. slowion) print *,'iat,fion',iat,(fion(j1,iat),j1=1,3)
        !energy which comes from the self-interaction of the spread charge
        eself=eself+real(at%nelpsp(ityp)**2,gp)*0.5_gp*sqrt(1.d0/pi)/at%psppar(0,0,ityp)
     end do

     !if (nproc==1 .and. slowion) print *,'eself',eself
     
  end if

  !for the surfaces BC,
  !activate for the moment only the slow calculation of the ionic energy and forces
  !if (at%geocode == 'S' .or. at%geocode == 'P') slowion=.true.
  if (at%geocode == 'S') slowion=.true.
  !slowion=.true.
  
  if (slowion) then

     !case of slow ionic calculation
     !conditions for periodicity in the three directions
     perx=(at%geocode /= 'F')
     pery=(at%geocode == 'P')
     perz=(at%geocode /= 'F')

     call ext_buffers(perx,nbl1,nbr1)
     call ext_buffers(pery,nbl2,nbr2)
     call ext_buffers(perz,nbl3,nbr3)

     !the ions corresponds to gaussian charges disposed in the same way as the pseudopotentials

     !first calculate the self-energy and the forces
     !(the latter are zero for a symmetric grid distribution)

     !self energy initialisation
     eself=0.0_gp
     do iat=1,at%nat

        fion(1,iat)=0.0_gp
        fion(2,iat)=0.0_gp
        fion(3,iat)=0.0_gp

        ityp=at%iatype(iat)
        rloc=at%psppar(0,0,ityp)
        charge=real(at%nelpsp(ityp),gp)/(2.0_gp*pi*sqrt(2.0_gp*pi)*rloc**3)
        prefactor=real(at%nelpsp(ityp),gp)/(2.0_gp*pi*sqrt(2.0_gp*pi)*rloc**5)
        cutoff=10.0_gp*rloc

        !calculate the self energy of the isolated bc
        eself=eself+real(at%nelpsp(ityp),gp)**2/rloc

!!$        if (n3pi >0 ) then
!!$
!!$           call razero(n1i*n2i*n3pi,pot_ion)
!!$           rx=rxyz(1,iat) 
!!$           ry=rxyz(2,iat)
!!$           rz=rxyz(3,iat)
!!$
!!$           isx=floor((rx-cutoff)/hxh)
!!$           isy=floor((ry-cutoff)/hyh)
!!$           isz=floor((rz-cutoff)/hzh)
!!$
!!$           iex=ceiling((rx+cutoff)/hxh)
!!$           iey=ceiling((ry+cutoff)/hyh)
!!$           iez=ceiling((rz+cutoff)/hzh)
!!$
!!$           !these nested loops will be used also for the actual ionic forces, to be recalculated
!!$           do i3=isz,iez
!!$              z=real(i3,gp)*hzh-rz
!!$              call ind_positions(perz,i3,n3,j3,goz) 
!!$              j3=j3+nbl3+1
!!$              do i2=isy,iey
!!$                 y=real(i2,gp)*hyh-ry
!!$                 call ind_positions(pery,i2,n2,j2,goy)
!!$                 do i1=isx,iex
!!$                    x=real(i1,gp)*hxh-rx
!!$                    call ind_positions(perx,i1,n1,j1,gox)
!!$                    r2=x**2+y**2+z**2
!!$                    arg=r2/rloc**2
!!$                    xp=exp(-.5_gp*arg)
!!$                    if (j3 >= i3s .and. j3 <= i3s+n3pi-1  .and. goy  .and. gox ) then
!!$                       ind=j1+1+nbl1+(j2+nbl2)*n1i+(j3-i3s+1-1)*n1i*n2i
!!$                       pot_ion(ind)=pot_ion(ind)-real(xp*charge,dp)
!!$                    endif
!!$                 end do
!!$              end do
!!$           end do
!!$        end if
!!$
!!$       
!!$        !application of the Poisson solver to calculate the self energy and the potential
!!$        !here the value of the datacode must be kept fixed
!!$        call PSolver(at%geocode,'D',iproc,nproc,n1i,n2i,n3i,0,hxh,hyh,hzh,&
!!$             pot_ion,pkernel,pot_ion,ehart,zero,zero,&
!!$             2.0_gp*pi*rloc**2*real(at%nelpsp(ityp),gp),.false.,1)
!!$        eself=eself+ehart
!!$
!!$        !initialise forces calculation
!!$        fxslf=0.0_gp
!!$        fyslf=0.0_gp
!!$        fzslf=0.0_gp
!!$
!!$        if (n3pi >0 ) then
!!$           do i3=isz,iez
!!$              z=real(i3,gp)*hzh-rz
!!$              call ind_positions(perz,i3,n3,j3,goz) 
!!$              j3=j3+nbl3+1
!!$              do i2=isy,iey
!!$                 y=real(i2,gp)*hyh-ry
!!$                 call ind_positions(pery,i2,n2,j2,goy)
!!$                 do i1=isx,iex
!!$                    x=real(i1,gp)*hxh-rx
!!$                    call ind_positions(perx,i1,n1,j1,gox)
!!$                    r2=x**2+y**2+z**2
!!$                    arg=r2/rloc**2
!!$                    xp=exp(-.5_gp*arg)
!!$                    if (j3 >= i3s .and. j3 <= i3s+n3pi-1  .and. goy  .and. gox ) then
!!$                       ind=j1+1+nbl1+(j2+nbl2)*n1i+(j3-i3s+1-1)*n1i*n2i
!!$                       !error function part
!!$                       Vel=pot_ion(ind)
!!$                       fxslf=fxslf+xp*Vel*x
!!$                       fyslf=fyslf+xp*Vel*y
!!$                       fzslf=fzslf+xp*Vel*z
!!$                    endif
!!$                 end do
!!$              end do
!!$           end do
!!$        end if
!!$
!!$        fion(1,iat)=-hxh*hyh*hzh*prefactor*fxslf
!!$        fion(2,iat)=-hxh*hyh*hzh*prefactor*fyslf
!!$        fion(3,iat)=-hxh*hyh*hzh*prefactor*fzslf
!!$
!!$        !if (nproc==1) print *,'iat,fself',iat,fxslf,fyslf,fzslf

     enddo

     eself=0.5_gp/sqrt(pi)*eself

     !if (nproc==1) 
     !print *,'iproc,eself',iproc,eself
     call razero(n1i*n2i*n3pi,pot_ion)

     if (n3pi >0 ) then
        !then calculate the hartree energy and forces of the charge distributions
        !(and save the values for the ionic potential)

        do iat=1,at%nat
           ityp=at%iatype(iat)
           rx=rxyz(1,iat) 
           ry=rxyz(2,iat)
           rz=rxyz(3,iat)

           rloc=at%psppar(0,0,ityp)
           charge=real(at%nelpsp(ityp),gp)/(2.0_gp*pi*sqrt(2.0_gp*pi)*rloc**3)
           cutoff=10.0_gp*rloc

           isx=floor((rx-cutoff)/hxh)
           isy=floor((ry-cutoff)/hyh)
           isz=floor((rz-cutoff)/hzh)

           iex=ceiling((rx+cutoff)/hxh)
           iey=ceiling((ry+cutoff)/hyh)
           iez=ceiling((rz+cutoff)/hzh)

           !these nested loops will be used also for the actual ionic forces, to be recalculated
           do i3=isz,iez
              z=real(i3,gp)*hzh-rz
              call ind_positions(perz,i3,n3,j3,goz) 
              j3=j3+nbl3+1
              do i2=isy,iey
                 y=real(i2,gp)*hyh-ry
                 call ind_positions(pery,i2,n2,j2,goy)
                 do i1=isx,iex
                    x=real(i1,gp)*hxh-rx
                    call ind_positions(perx,i1,n1,j1,gox)
                    r2=x**2+y**2+z**2
                    arg=r2/rloc**2
                    xp=exp(-.5_gp*arg)
                    if (j3 >= i3s .and. j3 <= i3s+n3pi-1  .and. goy  .and. gox ) then
                       ind=j1+1+nbl1+(j2+nbl2)*n1i+(j3-i3s+1-1)*n1i*n2i
                       pot_ion(ind)=pot_ion(ind)-xp*charge
                    endif
                 enddo
              enddo
           enddo

        enddo

     end if

     !now call the Poisson Solver for the global energy forces
     call PSolver(at%geocode,'D',iproc,nproc,n1i,n2i,n3i,0,hxh,hyh,hzh,&
          pot_ion,pkernel,pot_ion,ehart,zero,zero,-2.0_gp*psoffset,.false.,1)

     eion=ehart-eself

     !print *,'ehart,eself',iproc,ehart,eself

     !if (nproc==1) 
     !print *,'iproc,eion',iproc,eion


     do iat=1,at%nat
        ityp=at%iatype(iat)
        !coordinates of the center
        rx=rxyz(1,iat) 
        ry=rxyz(2,iat) 
        rz=rxyz(3,iat)
        !inizialization of the forces
        fxerf=0.0_gp
        fyerf=0.0_gp
        fzerf=0.0_gp

        !local part
        rloc=at%psppar(0,0,ityp)
        prefactor=real(at%nelpsp(ityp),gp)/(2.0_gp*pi*sqrt(2.0_gp*pi)*rloc**5)
        !maximum extension of the gaussian
        cutoff=10.0_gp*rloc

        isx=floor((rx-cutoff)/hxh)
        isy=floor((ry-cutoff)/hyh)
        isz=floor((rz-cutoff)/hzh)

        iex=ceiling((rx+cutoff)/hxh)
        iey=ceiling((ry+cutoff)/hyh)
        iez=ceiling((rz+cutoff)/hzh)

        !calculate the forces near the atom due to the error function part of the potential
        !calculate forces for all atoms only in the distributed part of the simulation box
        if (n3pi >0 ) then
           do i3=isz,iez
              z=real(i3,gp)*hzh-rz
              call ind_positions(perz,i3,n3,j3,goz) 
              j3=j3+nbl3+1
              do i2=isy,iey
                 y=real(i2,gp)*hyh-ry
                 call ind_positions(pery,i2,n2,j2,goy)
                 do i1=isx,iex
                    x=real(i1,gp)*hxh-rx
                    call ind_positions(perx,i1,n1,j1,gox)
                    r2=x**2+y**2+z**2
                    arg=r2/rloc**2
                    xp=exp(-.5_gp*arg)
                    if (j3 >= i3s .and. j3 <= i3s+n3pi-1  .and. goy  .and. gox ) then
                       ind=j1+1+nbl1+(j2+nbl2)*n1i+(j3-i3s+1-1)*n1i*n2i
                       !error function part
                       Vel=pot_ion(ind)
                       fxerf=fxerf+xp*Vel*x
                       fyerf=fyerf+xp*Vel*y
                       fzerf=fzerf+xp*Vel*z
                    endif
                 end do
              end do
           end do
        end if
        !final result of the forces

        fion(1,iat)=fion(1,iat)+(hxh*hyh*hzh*prefactor)*fxerf
        fion(2,iat)=fion(2,iat)+(hxh*hyh*hzh*prefactor)*fyerf
        fion(3,iat)=fion(3,iat)+(hxh*hyh*hzh*prefactor)*fzerf

        !if (nproc==1) print *,'iat,fion',iat,(fion(j1,iat),j1=1,3)

!!$        write(10+iat,'(1x,f8.3,i5,(1x,3(1x,1pe12.5)))',advance='no') &
!!$             hxh,iat,(fion(j1,iat),j1=1,3)


     end do

     if (nproc > 1) then
        allocate(gion(3,at%nat+ndebug),stat=i_stat)
        call memocc(i_stat,gion,'gion',subname)
        do iat=1,at%nat
           gion(1,iat)=fion(1,iat)
           gion(2,iat)=fion(2,iat)
           gion(3,iat)=fion(3,iat)
        end do

        call MPI_ALLREDUCE(gion,fion,3*at%nat,mpidtypg,MPI_SUM,MPI_COMM_WORLD,ierr)

        i_all=-product(shape(gion))*kind(gion)
        deallocate(gion,stat=i_stat)
        call memocc(i_stat,i_all,'gion',subname)

     end if

     !if (iproc ==0) print *,'eion',eion,psoffset,shortlength

  end if

  if (iproc == 0) write(*,'(1x,a,1pe22.14)') 'ion-ion interaction energy',eion
end subroutine IonicEnergyandForces


subroutine createIonicPotential(geocode,iproc,nproc,at,rxyz,&
     hxh,hyh,hzh,ef,n1,n2,n3,n3pi,i3s,n1i,n2i,n3i,pkernel,pot_ion,eion,psoffset)
  use module_base
  use module_types
  use Poisson_Solver
  implicit none
  character(len=1), intent(in) :: geocode
  integer, intent(in) :: iproc,nproc,n1,n2,n3,n3pi,i3s,n1i,n2i,n3i
  real(gp), intent(in) :: hxh,hyh,hzh,psoffset
  type(atoms_data), intent(in) :: at
  real(gp), dimension(3), intent(in) :: ef
  real(gp), dimension(3,at%nat), intent(in) :: rxyz
  real(dp), dimension(*), intent(in) :: pkernel
  real(wp), dimension(*), intent(inout) :: pot_ion
  real(gp), intent(out) :: eion
  !local variables
  logical :: perx,pery,perz,gox,goy,goz,htoobig=.false.,efwrite
  integer :: iat,jat,i1,i2,i3,j1,j2,j3,isx,isy,isz,iex,iey,iez,ierr,ityp,jtyp,nspin
  integer :: ind,i_all,i_stat,nbl1,nbr1,nbl2,nbr2,nbl3,nbr3,nloc,iloc
  real(kind=8) :: hgridh,pi,rholeaked,dist,rloc,charge,cutoff,x,y,z,r2,arg,xp,tt,rx,ry,rz
  real(kind=8) :: tt_tot,rholeaked_tot,eself,potxyz
  real(gp) :: ehart,eexcu,vexcu,elecfield,ystart,yend
  real(kind=8), dimension(4) :: charges_mpi

  call timing(iproc,'CrtLocPot     ','ON')

  if (iproc.eq.0) then
     write(*,'(1x,a)')&
          '----------------------------------------------------------- Ionic Potential Creation'
  end if

  pi=4.d0*atan(1.d0)
  ! Ionic charge (must be calculated for the PS active processes)
  rholeaked=0.d0
  ! Ionic energy (can be calculated for all the processors)

  !Creates charge density arising from the ionic PSP cores
  call razero(n1i*n2i*n3pi,pot_ion)

  !conditions for periodicity in the three directions
  perx=(geocode /= 'F')
  pery=(geocode == 'P')
  perz=(geocode /= 'F')

  call ext_buffers(perx,nbl1,nbr1)
  call ext_buffers(pery,nbl2,nbr2)
  call ext_buffers(perz,nbl3,nbr3)

  if (n3pi >0 .and. .not. htoobig) then

     do iat=1,at%nat
        ityp=at%iatype(iat)
        rx=rxyz(1,iat) 
        ry=rxyz(2,iat)
        rz=rxyz(3,iat)

        rloc=at%psppar(0,0,ityp)
        charge=real(at%nelpsp(ityp),kind=8)/(2.d0*pi*sqrt(2.d0*pi)*rloc**3)
        cutoff=10.d0*rloc

        isx=floor((rx-cutoff)/hxh)
        isy=floor((ry-cutoff)/hyh)
        isz=floor((rz-cutoff)/hzh)

        iex=ceiling((rx+cutoff)/hxh)
        iey=ceiling((ry+cutoff)/hyh)
        iez=ceiling((rz+cutoff)/hzh)

        !these nested loops will be used also for the actual ionic forces, to be recalculated
        do i3=isz,iez
           z=real(i3,kind=8)*hzh-rz
           call ind_positions(perz,i3,n3,j3,goz) 
           j3=j3+nbl3+1
           do i2=isy,iey
              y=real(i2,kind=8)*hyh-ry
              call ind_positions(pery,i2,n2,j2,goy)
              do i1=isx,iex
                 x=real(i1,kind=8)*hxh-rx
                 call ind_positions(perx,i1,n1,j1,gox)
                 r2=x**2+y**2+z**2
                 arg=r2/rloc**2
                 xp=exp(-.5d0*arg)
                 if (j3 >= i3s .and. j3 <= i3s+n3pi-1  .and. goy  .and. gox ) then
                    ind=j1+1+nbl1+(j2+nbl2)*n1i+(j3-i3s+1-1)*n1i*n2i
                    pot_ion(ind)=pot_ion(ind)-xp*charge
                 else if (.not. goz ) then
                    rholeaked=rholeaked+xp*charge
                 endif
              enddo
           enddo
        enddo

     enddo

  end if

  ! Check
  tt=0.d0
  do j3=1,n3pi
     do i2= -nbl2,2*n2+1+nbr2
        do i1= -nbl1,2*n1+1+nbr1
           ind=i1+1+nbl1+(i2+nbl2)*n1i+(j3-1)*n1i*n2i
           tt=tt+pot_ion(ind)
        enddo
     enddo
  enddo

  tt=tt*hxh*hyh*hzh
  rholeaked=rholeaked*hxh*hyh*hzh

  !print *,'test case input_rho_ion',iproc,i3start,i3end,n3pi,2*n3+16,tt

  if (nproc > 1) then
     charges_mpi(1)=tt
     charges_mpi(2)=rholeaked

     call MPI_ALLREDUCE(charges_mpi(1),charges_mpi(3),2,mpidtypd, &
          MPI_SUM,MPI_COMM_WORLD,ierr)

     tt_tot=charges_mpi(3)
     rholeaked_tot=charges_mpi(4)
  else
     tt_tot=tt
     rholeaked_tot=rholeaked
  end if

  if (iproc.eq.0) write(*,'(1x,a,f26.12,2x,1pe10.3)') &
       'total ionic charge, leaked charge ',tt_tot,rholeaked_tot

  if (.not. htoobig) then
     call timing(iproc,'CrtLocPot     ','OF')
     !here the value of the datacode must be kept fixed
     nspin=1
     call PSolver(geocode,'D',iproc,nproc,n1i,n2i,n3i,0,hxh,hyh,hzh,&
          pot_ion,pkernel,pot_ion,ehart,eexcu,vexcu,-psoffset,.false.,nspin)
     call timing(iproc,'CrtLocPot     ','ON')
  end if

  !print *,'ehart',ehart
  !print *,'true eion',ehart-eself

!!$  !calculate the value of the offset to be put
!!$  tt_tot=0.d0
!!$  do ind=1,n1i*n2i*n3i
!!$     tt_tot=tt_tot+pot_ion(ind)
!!$  end do
!!$  print *,'previous offset',tt_tot*hxh*hyh*hzh

  if (n3pi > 0) then
     do iat=1,at%nat
        ityp=at%iatype(iat)

        rx=rxyz(1,iat)
        ry=rxyz(2,iat)
        rz=rxyz(3,iat)

        ! determine number of local terms
        nloc=0
        do iloc=1,4
           if (at%psppar(0,iloc,ityp) /= 0.d0) nloc=iloc
        enddo
        rloc=at%psppar(0,0,ityp)
        cutoff=10.d0*rloc

        isx=floor((rx-cutoff)/hxh)
        isy=floor((ry-cutoff)/hyh)
        isz=floor((rz-cutoff)/hzh)

        iex=ceiling((rx+cutoff)/hxh)
        iey=ceiling((ry+cutoff)/hyh)
        iez=ceiling((rz+cutoff)/hzh)
        
        if (nloc /= 0) then

           do i3=isz,iez
              z=real(i3,kind=8)*hzh-rz
              call ind_positions(perz,i3,n3,j3,goz) 
              j3=j3+nbl3+1
              if (goz .and. j3 >= i3s .and. j3 <=  i3s+n3pi-1) then
                 do i2=isy,iey
                    y=real(i2,kind=8)*hyh-ry
                    call ind_positions(pery,i2,n2,j2,goy)
                    if (goy) then
                       do i1=isx,iex
                          x=real(i1,kind=8)*hxh-rx
                          call ind_positions(perx,i1,n1,j1,gox)
                          if (gox) then
                             r2=x**2+y**2+z**2
                             arg=r2/rloc**2
                             xp=exp(-.5d0*arg)
                             tt=at%psppar(0,nloc,ityp)
                             do iloc=nloc-1,1,-1
                                tt=arg*tt+at%psppar(0,iloc,ityp)
                             enddo
                             ind=j1+1+nbl1+(j2+nbl2)*n1i+(j3-i3s+1-1)*n1i*n2i
                             pot_ion(ind)=pot_ion(ind)+xp*tt
                          end if
                       enddo
                    end if
                 enddo
              end if
           end do

        end if

     enddo


     if (htoobig) then
        !add to pot_ion an explicit error function to correct in the case of big grid spacing
        !for the moment works only in the isolated BC case
        do i3=1,n3pi
           z=real(i3+i3s-1-nbl3-1,gp)*hzh
           do i2=1,n2i
              y=real(i2-nbl2-1,gp)*hyh
              do i1=1,n1i
                 x=real(i1-nbl1-1,gp)*hxh
                 ind=i1+(i2-1)*n1i+(i3-1)*n1i*n2i
                 call sum_erfcr(at%nat,at%ntypes,x,y,z,at%iatype,at%nelpsp,at%psppar,rxyz,potxyz)
                 pot_ion(ind)=pot_ion(ind)+potxyz
              end do
           end do
        end do
     end if
     

  end if

!!$  !calculate the value of the offset to be put
!!$  tt_tot=0.d0
!!$  do ind=1,n1i*n2i*n3i
!!$     tt_tot=tt_tot+pot_ion(ind)
!!$  end do
!!$  print *,'actual offset',tt_tot*hxh*hyh*hzh

  elecfield=ef(1)
  !correct ystart and yend to the whole simulation box 
  ystart=max(ef(2),0.0_gp)
  yend=min(ef(3),at%alat2)

  !use rhopot to calculate the potential from a constant electric field along y direction
  if (elecfield /= 0.0_gp) then
     !constant electric field allowed only for free BC
     if (geocode == 'P') then
     if (iproc == 0) write(*,'(1x,a)') &
          'The constant electric field is allowed only for Free and Surfaces BC'
     stop
     end if
     if (iproc.eq.0) write(*,'(1x,3(a,1pe10.2),a)') &
          'Constant electric field of',elecfield,&
          ' Ha*Bohr for:',ystart,' < y <',yend,' Bohr'

     !write or not electric field in a separate file
     efwrite=.true.

     if (n3pi > 0) then
        do i3=1,n3pi
           !z=real(i3+i3s-1-nbl3-1,gp)*hzh
           do i2=1,n2i
              y=real(i2-nbl2-1,gp)*hyh
              if (y < ystart) then
                 do i1=1,n1i
                    !x=real(i1-nbl1-1,gp)*hxh
                    ind=i1+(i2-1)*n1i+(i3-1)*n1i*n2i
                    pot_ion(ind)=pot_ion(ind)+elecfield*(0.5_gp*(ystart-yend))
                 end do
              else if (y > yend) then
                 do i1=1,n1i
                    !x=real(i1-nbl1-1,gp)*hxh
                    ind=i1+(i2-1)*n1i+(i3-1)*n1i*n2i
                    pot_ion(ind)=pot_ion(ind)-elecfield*(0.5_gp*(ystart-yend))
                 end do
              else
                 do i1=1,n1i
                    !x=real(i1-nbl1-1,gp)*hxh
                    ind=i1+(i2-1)*n1i+(i3-1)*n1i*n2i
                    pot_ion(ind)=pot_ion(ind)+elecfield*(y-0.5_gp*(ystart+yend))
                 end do
              end if
           end do
        end do

        if (efwrite .and. iproc == 0) then
           open(unit=17,file='elecfield_y',status='unknown')
           do i2=nbl2+1,n2i-nbr2-1
              y=real(i2-nbl2-1,gp)*hyh
              if (y < ystart) then
                 write(17,*)i2,y,elecfield*(0.5_gp*(ystart-yend))
              else if (y > yend) then
                 write(17,*)i2,y,-elecfield*(0.5_gp*(ystart-yend))
              else
                 write(17,*)i2,y,elecfield*(y-0.5_gp*(ystart+yend))
              end if
           end do
           close(17)
        end if

     end if
  end if

  call timing(iproc,'CrtLocPot     ','OF')

end subroutine createIonicPotential
!determine the index in which the potential must be inserted, following the BC
!determine also whether the index is inside or outside the box for free BC
subroutine ind_positions(periodic,i,n,j,go)
  implicit none
  logical, intent(in) :: periodic
  integer, intent(in) :: i,n
  logical, intent(out) :: go
  integer, intent(out) :: j

  if (periodic) then
     go=.true.
     j=modulo(i,2*n+2)
  else
     j=i
     if (i >= -14 .and. i <= 2*n+16) then
        go=.true.
     else
        go=.false.
     end if
  end if

end subroutine ind_positions

subroutine sum_erfcr(nat,ntypes,x,y,z,iatype,nelpsp,psppar,rxyz,potxyz)
  use module_base
  implicit none
  integer, intent(in) :: nat,ntypes
  real(gp) :: x,y,z
  integer, dimension(nat), intent(in) :: iatype
  integer, dimension(ntypes), intent(in) :: nelpsp
  real(gp), dimension(0:4,0:6,ntypes), intent(in) :: psppar
  real(gp), dimension(3,nat), intent(in) :: rxyz
  real(wp), intent(out) :: potxyz
  !local variables
  integer :: iat,ityp
  real(wp) :: pi,charge,erf
  real(gp) :: r,sq2rl,rx,ry,rz,derf_val
  

  pi=4.0_wp*atan(1.0_wp)

  potxyz =0.0_wp

  do iat=1,nat

     ityp=iatype(iat)
     sq2rl=sqrt(2.0_gp)*psppar(0,0,ityp)
     charge=real(nelpsp(ityp),wp)

     rx=rxyz(1,iat)-x 
     ry=rxyz(2,iat)-y
     rz=rxyz(3,iat)-z

     r=sqrt(rx**2+ry**2+rz**2)

     if (r == 0.0_gp) then
        potxyz = potxyz - charge*2.0_wp/(sqrt(pi)*real(sq2rl,wp))
     else
        call derf(derf_val,r/sq2rl)
        potxyz = potxyz - charge*real(derf_val/r,wp)
     end if

  end do

end subroutine sum_erfcr

subroutine ext_buffers(periodic,nl,nr)
  implicit none
  logical, intent(in) :: periodic
  integer, intent(out) :: nl,nr

  if (periodic) then
     nl=0
     nr=0
  else
     nl=14
     nr=15
  end if
end subroutine ext_buffers

subroutine derf(derf_yy,yy)

 use module_base
 implicit none
 real(gp),intent(in) :: yy
 real(gp),intent(out) :: derf_yy
 integer          ::  done,ii,isw
 real(gp), parameter :: &
       ! coefficients for 0.0 <= yy < .477
       &  pp(5)=(/ 113.8641541510502e0_gp, 377.4852376853020e0_gp,  &
       &           3209.377589138469e0_gp, .1857777061846032e0_gp,  &
       &           3.161123743870566e0_gp /)
  real(gp), parameter :: &
       &  qq(4)=(/ 244.0246379344442e0_gp, 1282.616526077372e0_gp,  &
       &           2844.236833439171e0_gp, 23.60129095234412e0_gp/)
  ! coefficients for .477 <= yy <= 4.0
  real(gp), parameter :: &
       &  p1(9)=(/ 8.883149794388376e0_gp, 66.11919063714163e0_gp,  &
       &           298.6351381974001e0_gp, 881.9522212417691e0_gp,  &
       &           1712.047612634071e0_gp, 2051.078377826071e0_gp,  &
       &           1230.339354797997e0_gp, 2.153115354744038e-8_gp, &
       &           .5641884969886701e0_gp /)
  real(gp), parameter :: &
       &  q1(8)=(/ 117.6939508913125e0_gp, 537.1811018620099e0_gp,  &
       &           1621.389574566690e0_gp, 3290.799235733460e0_gp,  &
       &           4362.619090143247e0_gp, 3439.367674143722e0_gp,  &
       &           1230.339354803749e0_gp, 15.74492611070983e0_gp/)
  ! coefficients for 4.0 < y,
  real(gp), parameter :: &
       &  p2(6)=(/ -3.603448999498044e-01_gp, -1.257817261112292e-01_gp,   &
       &           -1.608378514874228e-02_gp, -6.587491615298378e-04_gp,   &
       &           -1.631538713730210e-02_gp, -3.053266349612323e-01_gp/)
  real(gp), parameter :: &
       &  q2(5)=(/ 1.872952849923460e0_gp   , 5.279051029514284e-01_gp,    &
       &           6.051834131244132e-02_gp , 2.335204976268692e-03_gp,    &
       &           2.568520192289822e0_gp /)
  real(gp), parameter :: &
       &  sqrpi=.5641895835477563e0_gp, xbig=13.3e0_gp, xlarge=6.375e0_gp, xmin=1.0e-10_gp
  real(gp) ::  res,xden,xi,xnum,xsq,xx

 xx = yy
 isw = 1
!Here change the sign of xx, and keep track of it thanks to isw
 if (xx<0.0e0_gp) then
  isw = -1
  xx = -xx
 end if

 done=0

!Residual value, if yy < -6.375e0_gp
 res=-1.0e0_gp

!abs(yy) < .477, evaluate approximation for erfc
 if (xx<0.477e0_gp) then
! xmin is a very small number
  if (xx<xmin) then
   res = xx*pp(3)/qq(3)
  else
   xsq = xx*xx
   xnum = pp(4)*xsq+pp(5)
   xden = xsq+qq(4)
   do ii = 1,3
    xnum = xnum*xsq+pp(ii)
    xden = xden*xsq+qq(ii)
   end do
   res = xx*xnum/xden
  end if
  if (isw==-1) res = -res
  done=1
 end if

!.477 < abs(yy) < 4.0 , evaluate approximation for erfc
 if (xx<=4.0e0_gp .and. done==0 ) then
  xsq = xx*xx
  xnum = p1(8)*xx+p1(9)
  xden = xx+q1(8)
  do ii=1,7
   xnum = xnum*xx+p1(ii)
   xden = xden*xx+q1(ii)
  end do
  res = xnum/xden
  res = res* exp(-xsq)
  if (isw.eq.-1) res = res-1.0e0_gp
  done=1
 end if

!y > 13.3e0_gp
 if (isw > 0 .and. xx > xbig .and. done==0 ) then
  res = 1.0e0_gp
  done=1
 end if

!4.0 < yy < 13.3e0_gp  .or. -6.375e0_gp < yy < -4.0
!evaluate minimax approximation for erfc
 if ( ( isw > 0 .or. xx < xlarge ) .and. done==0 ) then
  xsq = xx*xx
  xi = 1.0e0_gp/xsq
  xnum= p2(5)*xi+p2(6)
  xden = xi+q2(5)
  do ii = 1,4
   xnum = xnum*xi+p2(ii)
   xden = xden*xi+q2(ii)
  end do
  res = (sqrpi+xi*xnum/xden)/xx
  res = res* exp(-xsq)
  if (isw.eq.-1) res = res-1.0e0_gp
 end if

!All cases have been investigated
 derf_yy = res

end subroutine derf
