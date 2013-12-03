!> @file
!!  Routines for the ionic energy contribution
!! @author
!!    Copyright (C) 2007-2013 BigDFT group (LG)
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 


!> Calculate the ionic contribution to the energy and the forces
subroutine IonicEnergyandForces(iproc,nproc,dpbox,at,elecfield,&
     & rxyz,eion,fion,dispersion,edisp,fdisp,ewaldstr,n1,n2,n3,&
     & pot_ion,pkernel,psoffset)
  use module_base
  use module_types
  use Poisson_Solver, except_dp => dp, except_gp => gp, except_wp => wp
  use vdwcorrection
  use yaml_output
  implicit none
  type(denspot_distribution), intent(in) :: dpbox
  type(atoms_data), intent(in) :: at
  integer, intent(in) :: iproc,nproc,n1,n2,n3,dispersion
  real(gp), dimension(3), intent(in) :: elecfield
  real(gp), dimension(3,at%astruct%nat), intent(in) :: rxyz
  type(coulomb_operator), intent(in) :: pkernel
  real(gp), intent(out) :: eion,edisp,psoffset
  real(dp), dimension(6),intent(out) :: ewaldstr
  real(gp), dimension(:,:), pointer :: fion,fdisp
  real(dp), dimension(*), intent(out) :: pot_ion
  !local variables
  character(len=*), parameter :: subname='IonicEnergyandForces'
  logical :: slowion=.false.
  logical :: perx,pery,perz,gox,goy,goz
  integer :: n1i,n2i,n3i,i3s,n3pi
  integer :: i,iat,ii,i_all,i_stat,ityp,jat,jtyp,nbl1,nbr1,nbl2,nbr2,nbl3,nbr3
  integer :: isx,iex,isy,iey,isz,iez,i1,i2,i3,j1,j2,j3,ind,ierr
  real(gp) :: ucvol,rloc,twopitothreehalf,pi,atint,shortlength,charge,eself,rx,ry,rz
  real(gp) :: fxion,fyion,fzion,dist,fxerf,fyerf,fzerf,cutoff
  real(gp) :: hxh,hyh,hzh
  real(gp) :: hxx,hxy,hxz,hyy,hyz,hzz,chgprod
  real(gp) :: x,y,z,xp,Vel,prefactor,r2,arg,ehart,de
  !real(gp) :: Mz,cmassy
  real(gp), dimension(3,3) :: gmet,rmet,rprimd,gprimd
  !other arrays for the ewald treatment
  real(gp), dimension(:,:), allocatable :: fewald,xred
  real(gp), dimension(3) :: cc

  allocate(fion(3,at%astruct%nat+ndebug),stat=i_stat)
  call memocc(i_stat,fion,'fion',subname)
  allocate(fdisp(3,at%astruct%nat+ndebug),stat=i_stat)
  call memocc(i_stat,fdisp,'fdisp',subname)

  ! Aliasing
  hxh = dpbox%hgrids(1)
  hyh = dpbox%hgrids(2)
  hzh = dpbox%hgrids(3)
  n1i = dpbox%ndims(1)
  n2i = dpbox%ndims(2)
  n3i = dpbox%ndims(3)
  i3s = dpbox%i3s+dpbox%i3xcsh
  n3pi = dpbox%n3pi

  pi=4.d0*datan(1.d0)
  psoffset=0.0_gp
  ewaldstr=0.0_gp
  if (at%astruct%geocode == 'P') then
     !here we insert the calculation of the ewald forces
     allocate(fewald(3,at%astruct%nat+ndebug),stat=i_stat)
     call memocc(i_stat,fewald,'fewald',subname)
     allocate(xred(3,at%astruct%nat+ndebug),stat=i_stat)
     call memocc(i_stat,xred,'xred',subname)

     !calculate rprimd
     rprimd(:,:)=0.0_gp

     rprimd(1,1)=at%astruct%cell_dim(1)
     rprimd(2,2)=at%astruct%cell_dim(2)
     rprimd(3,3)=at%astruct%cell_dim(3)

     !calculate the metrics and the volume
     call metric(gmet,gprimd,-1,rmet,rprimd,ucvol)

     !calculate reduced coordinates
     do iat=1,at%astruct%nat
        do ii=1,3
           xred(ii,iat)= gprimd(1,ii)*rxyz(1,iat)+gprimd(2,ii)*rxyz(2,iat)+&
                gprimd(3,ii)*rxyz(3,iat)
        end do
     end do

     !calculate ewald energy and forces + stress
     call ewald(eion,gmet,fewald,at%astruct%nat,at%astruct%ntypes,rmet,at%astruct%iatype,ucvol,&
          xred,real(at%nelpsp,kind=8))
     ewaldstr=0.0_dp
     call ewald2(gmet,at%astruct%nat,at%astruct%ntypes,rmet,rprimd,ewaldstr,at%astruct%iatype,&
          ucvol,xred,real(at%nelpsp,kind=8))

! our sequence of strten elements : 11 22 33 12 13 23
! abinit output                   : 11 22 33 23 13 12

     !make forces dimensional
     do iat=1,at%astruct%nat
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
     do iat=1,at%astruct%nat
        ityp=at%astruct%iatype(iat)
        rloc=at%psppar(0,0,ityp)
        atint=at%psppar(0,1,ityp)+3.0_gp*at%psppar(0,2,ityp)+&
             15.0_gp*at%psppar(0,3,ityp)+105.0_gp*at%psppar(0,4,ityp)
        psoffset=psoffset+rloc**3*atint
        shortlength=shortlength+real(at%nelpsp(ityp),gp)*rloc**2
        charge=charge+real(at%nelpsp(ityp),gp)
     end do
     psoffset=twopitothreehalf*psoffset
     shortlength=shortlength*2.d0*pi

     !print *,'psoffset',psoffset,'pspcore',(psoffset+shortlength)*charge/(at%astruct%cell_dim(1)*at%astruct%cell_dim(2)*at%astruct%cell_dim(3))
     !if (iproc ==0) print *,'eion',eion,charge/ucvol*(psoffset+shortlength)
     !correct ionic energy taking into account the PSP core correction
     eion=eion+charge/ucvol*(psoffset+shortlength)

     !symmetrization of ewald stress (probably not needed)
     if (at%astruct%sym%symObj >= 0) call symm_stress((iproc==0),ewaldstr,at%astruct%sym%symObj)
     !PSP core correction of the stress tensor (diag.)
     ewaldstr(1:3)=ewaldstr(1:3)-charge*(psoffset+shortlength)/ucvol/ucvol

!!$     if (iproc == 0) then
!!$        write(*,*) 'STRESS TENSOR: EWALD + PSP-CORE'
!!$        write(*,*) ewaldstr(1:3)
!!$        write(*,*) ewaldstr(6),ewaldstr(5),ewaldstr(4)
!!$     end if

!!!     !in the surfaces case, correct the energy term following (J.Chem.Phys. 111(7)-3155, 1999)
!!!     if (at%astruct%geocode == 'S') then
!!!        !calculate the Mz dipole component (which in our case corresponds to y direction)
!!!        !first calculate the center of mass
!!!        cmassy=0.0_gp
!!!        do iat=1,at%astruct%nat
!!!           cmassy=cmassy+rxyz(2,iat)
!!!        end do
!!!        
!!!        Mz=0.0_gp
!!!        do iat=1,at%astruct%nat
!!!           ityp=at%astruct%iatype(iat)
!!!           Mz=Mz+real(at%nelpsp(ityp),gp)*(rxyz(2,iat)-cmassy)
!!!        end do
!!!        
!!!        !correct energy and forces in the y direction
!!!        eion=eion+0.5_gp/ucvol*Mz**2
!!!        do iat=1,at%astruct%nat
!!!           ityp=at%astruct%iatype(iat)
!!!           fion(2,iat)=fion(2,iat)-real(at%nelpsp(ityp),gp)/ucvol*Mz
!!!           if (nproc==1 .and. slowion) print *,'iat,fion',iat,(fion(j1,iat),j1=1,3)
!!!        end do
!!!
!!!     end if

  else if (at%astruct%geocode == 'F') then

     eion=0.0_gp
     eself=0.0_gp
     do iat=1,at%astruct%nat
        ityp=at%astruct%iatype(iat)
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
           jtyp=at%astruct%iatype(jat)
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
        do jat=iat+1,at%astruct%nat
           dist=sqrt((rx-rxyz(1,jat))**2+(ry-rxyz(2,jat))**2+(rz-rxyz(3,jat))**2)
           jtyp=at%astruct%iatype(jat)
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
  !if (at%astruct%geocode == 'S' .or. at%astruct%geocode == 'P') slowion=.true.
  if (at%astruct%geocode == 'S') slowion=.true.
  !slowion=.true.

  if (slowion) then

     !case of slow ionic calculation
     !conditions for periodicity in the three directions
     perx=(at%astruct%geocode /= 'F')
     pery=(at%astruct%geocode == 'P')
     perz=(at%astruct%geocode /= 'F')

     call ext_buffers(perx,nbl1,nbr1)
     call ext_buffers(pery,nbl2,nbr2)
     call ext_buffers(perz,nbl3,nbr3)

     !the ions corresponds to gaussian charges disposed in the same way as the pseudopotentials

     !first calculate the self-energy and the forces
     !(the latter are zero for a symmetric grid distribution)

     !self energy initialisation
     eself=0.0_gp
     do iat=1,at%astruct%nat

        fion(1,iat)=0.0_gp
        fion(2,iat)=0.0_gp
        fion(3,iat)=0.0_gp

        ityp=at%astruct%iatype(iat)
        rloc=at%psppar(0,0,ityp)
        charge=real(at%nelpsp(ityp),gp)/(2.0_gp*pi*sqrt(2.0_gp*pi)*rloc**3)
        prefactor=real(at%nelpsp(ityp),gp)/(2.0_gp*pi*sqrt(2.0_gp*pi)*rloc**5)
        cutoff=10.0_gp*rloc

        !calculate the self energy of the isolated bc
        eself=eself+real(at%nelpsp(ityp),gp)**2/rloc
     enddo

     eself=0.5_gp/sqrt(pi)*eself

     !if (nproc==1) 
     !print *,'iproc,eself',iproc,eself
     call to_zero(n1i*n2i*n3pi,pot_ion(1))

     if (n3pi >0 ) then
        !then calculate the hartree energy and forces of the charge distributions
        !(and save the values for the ionic potential)

        do iat=1,at%astruct%nat
           ityp=at%astruct%iatype(iat)
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
              !call ind_positions(perz,i3,n3,j3,goz) 
              call ind_positions_new(perz,i3,n3i,j3,goz) 
              j3=j3+nbl3+1
              do i2=isy,iey
                 y=real(i2,gp)*hyh-ry
                 !call ind_positions(pery,i2,n2,j2,goy)
                 call ind_positions_new(pery,i2,n2i,j2,goy)
                 do i1=isx,iex
                    x=real(i1,gp)*hxh-rx
                    !call ind_positions(perx,i1,n1,j1,gox)
                    call ind_positions_new(perx,i1,n1i,j1,gox)
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
     call H_potential('D',pkernel,pot_ion,pot_ion,ehart,-2.0_gp*psoffset,.false.)

     eion=ehart-eself

     !print *,'ehart,eself',iproc,ehart,eself

     !if (nproc==1) 
     !print *,'iproc,eion',iproc,eion

     do iat=1,at%astruct%nat
        ityp=at%astruct%iatype(iat)
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
              !call ind_positions(perz,i3,n3,j3,goz) 
              call ind_positions_new(perz,i3,n3i,j3,goz) 
              j3=j3+nbl3+1
              do i2=isy,iey
                 y=real(i2,gp)*hyh-ry
                 !call ind_positions(pery,i2,n2,j2,goy)
                 call ind_positions_new(pery,i2,n2i,j2,goy)
                 do i1=isx,iex
                    x=real(i1,gp)*hxh-rx
                    !call ind_positions(perx,i1,n1,j1,gox)
                    call ind_positions_new(perx,i1,n1i,j1,gox)
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

!!!        write(10+iat,'(1x,f8.3,i5,(1x,3(1x,1pe12.5)))',advance='no') &
!!!             hxh,iat,(fion(j1,iat),j1=1,3)


     end do

     if (pkernel%mpi_env%nproc > 1) then
        call mpiallred(fion(1,1),3*at%astruct%nat,MPI_SUM,pkernel%mpi_env%mpi_comm,ierr)
     end if

     !if (iproc ==0) print *,'eion',eion,psoffset,shortlength

  end if

  ! Add contribution from constant electric field to the forces
  call center_of_charge(at,rxyz,cc)
  do iat=1,at%astruct%nat
     ityp=at%astruct%iatype(iat)
     charge=real(at%nelpsp(ityp),gp)
     fion(1:3,iat)=fion(1:3,iat)+(charge*elecfield(1:3))
     !ry=rxyz(2,iat) 
     !eion=eion-(charge*elecfield)*ry
     de=0.0_gp
     do i=1,3
        de=de+elecfield(i)*(rxyz(i,iat)-cc(i))
     end do
     !eion=eion-charge*sum(elecfield(1:3)*rxyz(1:3,iat))
     eion=eion-charge*de
  enddo

  if (iproc == 0) then
     if(all(elecfield(1:3)==0._gp)) then 
     !      write(*,'(1x,a,1pe22.14)') 'ion-ion interaction energy',eion
        call yaml_map('Ion-Ion interaction energy',eion,fmt='(1pe22.14)')
     else 
     !      write(*,'(1x,a,1pe22.14)') 'ion-ion and ion-electric field interaction energy',eion
        call yaml_map('Ion-electric field interaction energy',eion,fmt='(1pe22.14)')
     endif
  end if

  ! Add empiric correction for Van der Waals forces and energy.
  call vdwcorrection_calculate_energy(edisp,rxyz,at,dispersion)
  if (iproc == 0 .and. edisp /= 0.0_gp) then
!!$     write(*,'(1x,a, e12.5,1x,a)') &
!!$          'Dispersion Correction Energy: ', dispersion_energy, 'Hartree'
     call yaml_map('Dispersion Correction Energy (Ha)',edisp,fmt='(1pe22.14)')
  end if

  call vdwcorrection_calculate_forces(fdisp,rxyz,at,dispersion)
  call vdwcorrection_freeparams() 
END SUBROUTINE IonicEnergyandForces


subroutine createEffectiveIonicPotential(iproc, nproc, verb, in, atoms, rxyz, shift, &
     & Glr, hxh, hyh, hzh, rhopotd, pkernel, pot_ion, elecfield, psoffset,rholoc)
  use module_base
  use module_types

  implicit none

  integer, intent(in) :: iproc,nproc
  logical, intent(in) :: verb
  real(gp), intent(in) :: hxh,hyh,hzh,psoffset
  type(atoms_data), intent(in) :: atoms
  type(locreg_descriptors), intent(in) :: Glr
  type(input_variables), intent(in) :: in
  type(denspot_distribution), intent(in) :: rhopotd
  real(gp), intent(in) :: elecfield(3)
  real(gp), dimension(3), intent(in) :: shift
  real(gp), dimension(3,atoms%astruct%nat), intent(in) :: rxyz
  type(coulomb_operator), intent(in) :: pkernel
  real(wp), dimension(*), intent(inout) :: pot_ion
  type(rholoc_objects),intent(in)::rholoc  !Object used for PAW+WVL
                                           !contains the local density

  character(len = *), parameter :: subname = "createEffectiveIonicPotential"
  logical :: counterions
  integer :: i_stat, i_all
  real(dp), dimension(:), allocatable :: counter_ions

  ! Compute the main ionic potential.
  call createIonicPotential(atoms%astruct%geocode, iproc, nproc, verb, atoms, rxyz, hxh, hyh, hzh, &
       & elecfield, Glr%d%n1, Glr%d%n2, Glr%d%n3, rhopotd%n3pi, rhopotd%i3s + rhopotd%i3xcsh, &
       & Glr%d%n1i, Glr%d%n2i, Glr%d%n3i, pkernel, pot_ion, psoffset,rholoc)

  !inquire for the counter_ion potential calculation (for the moment only xyz format)
  inquire(file='posinp_ci.xyz',exist=counterions)
  if (counterions) then
     if (rhopotd%n3pi > 0) then
        allocate(counter_ions(Glr%d%n1i*Glr%d%n2i*rhopotd%n3pi+ndebug),stat=i_stat)
        call memocc(i_stat,counter_ions,'counter_ions',subname)
     else
        allocate(counter_ions(1+ndebug),stat=i_stat)
        call memocc(i_stat,counter_ions,'counter_ions',subname)
     end if

     call CounterIonPotential(atoms%astruct%geocode,iproc,nproc,in,shift,&
          &   hxh,hyh,hzh,Glr%d,rhopotd%n3pi,rhopotd%i3s + rhopotd%i3xcsh,pkernel,counter_ions)

     !sum that to the ionic potential
     call axpy(Glr%d%n1i*Glr%d%n2i*rhopotd%n3pi,1.0_dp,counter_ions(1),1,&
          &   pot_ion(1),1)

     i_all=-product(shape(counter_ions))*kind(counter_ions)
     deallocate(counter_ions,stat=i_stat)
     call memocc(i_stat,i_all,'counter_ions',subname)
  end if
END SUBROUTINE createEffectiveIonicPotential


subroutine createIonicPotential(geocode,iproc,nproc,verb,at,rxyz,&
     hxh,hyh,hzh,elecfield,n1,n2,n3,n3pi,i3s,n1i,n2i,n3i,pkernel,pot_ion,psoffset,rholoc)
  use module_base
  use module_types
  use yaml_output
  !use gaussians
!  use module_interfaces, except_this_one => createIonicPotential
  use Poisson_Solver, except_dp => dp, except_gp => gp, except_wp => wp
  implicit none
  character(len=1), intent(in) :: geocode !< @copydoc poisson_solver::doc::geocode
  integer, intent(in) :: iproc,nproc,n1,n2,n3,n3pi,i3s,n1i,n2i,n3i
  logical, intent(in) :: verb
  real(gp), intent(in) :: hxh,hyh,hzh,psoffset
  type(atoms_data), intent(in) :: at
  real(gp), dimension(3), intent(in) :: elecfield
  real(gp), dimension(3,at%astruct%nat), intent(in) :: rxyz
  type(coulomb_operator), intent(in) :: pkernel
  real(wp), dimension(*), intent(inout) :: pot_ion
  type(rholoc_objects),intent(in)::rholoc

  !local variables
  character(len=*), parameter :: subname='createIonicPotential'
  character(len = 3) :: quiet
  logical :: perx,pery,perz,gox,goy,goz,htoobig=.false.,efwrite,check_potion=.false.
  integer :: iat,i1,i2,i3,j1,j2,j3,isx,isy,isz,iex,iey,iez,ierr,ityp !n(c) nspin
  integer :: ind,i_all,i_stat,nbl1,nbr1,nbl2,nbr2,nbl3,nbr3,nloc,iloc
  real(kind=8) :: pi,rholeaked,rloc,charge,cutoff,x,y,z,r2,arg,xp,tt,rx,ry,rz
  real(kind=8) :: tt_tot,rholeaked_tot,potxyz
  real(kind=8) :: raux,raux2,rr,r2paw
  real(wp) :: maxdiff
  real(gp) :: ehart
  real(dp), dimension(2) :: charges_mpi
  real(dp), dimension(:), allocatable :: potion_corr
  !real(dp), dimension(:), allocatable :: den_aux

  call timing(iproc,'CrtLocPot     ','ON')

  !call initialize_real_space_conversion() !initialize the work arrays needed to integrate with isf

  pi=4.d0*atan(1.d0)
  ! Ionic charge (must be calculated for the PS active processes)
  rholeaked=0.d0
  ! Ionic energy (can be calculated for all the processors)

  !Creates charge density arising from the ionic PSP cores
  call to_zero(n1i*n2i*n3pi,pot_ion(1))

  !conditions for periodicity in the three directions
  perx=(geocode /= 'F')
  pery=(geocode == 'P')
  perz=(geocode /= 'F')

  call ext_buffers(perx,nbl1,nbr1)
  call ext_buffers(pery,nbl2,nbr2)
  call ext_buffers(perz,nbl3,nbr3)

  if (n3pi >0 .and. .not. htoobig) then

     do iat=1,at%astruct%nat
        ityp=at%astruct%iatype(iat)
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

!       Calculate Ionic Density
!       using HGH parameters.
!       Eq. 1.104, T. Deutsch and L. Genovese, JDN. 12, 2011
        if( .not. any(at%npspcode == 7) ) then

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
                    !use multipole-preserving function
!!$                    xp=mp_exp(hxh,rx,0.5_gp/(rloc**2),i1,0,.true.)*&
!!$                         mp_exp(hyh,ry,0.5_gp/(rloc**2),i2,0,.true.)*&
!!$                         mp_exp(hzh,rz,0.5_gp/(rloc**2),i3,0,.true.)
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

!       Calculate Ionic Density using splines, 
!       PAW case
        else
             
           r2paw=rholoc%radius(ityp)**2
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
                    !if(r2>r2paw) cycle
                    rr=sqrt(r2)
                    if(1==2) then
                      !This converges very slow                
                      call splint(rholoc%msz(ityp),rholoc%rad(1:rholoc%msz(ityp),ityp),&
&                      rholoc%d(1:rholoc%msz(ityp),1,ityp),rholoc%d(1:rholoc%msz(ityp),2,ityp),1,rr,raux,ierr)
                    else
                      !Take the HGH form for rho_L (long range)
                      arg=r2/rloc**2
                      !use multipole-preserving function
!!$                      xp=mp_exp(hxh,rx,0.5_gp/(rloc**2),i1,0,.true.)*&
!!$                           mp_exp(hyh,ry,0.5_gp/(rloc**2),i2,0,.true.)*&
!!$                           mp_exp(hzh,rz,0.5_gp/(rloc**2),i3,0,.true.)
                      xp=exp(-.5d0*arg)
                      raux=-xp*charge
                    end if
                    !raux=-4.d0**(3.0d0/2.0d0)*exp(-4.d0*pi*r2)

                    if (j3 >= i3s .and. j3 <= i3s+n3pi-1  .and. goy  .and. gox ) then
                       ind=j1+1+nbl1+(j2+nbl2)*n1i+(j3-i3s+1-1)*n1i*n2i
                       pot_ion(ind)=pot_ion(ind)+raux
                    else if (.not. goz ) then
                       rholeaked=rholeaked-raux
                    endif
                 enddo
              enddo
           enddo

        end if



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

  if (pkernel%mpi_env%nproc > 1) then
     charges_mpi(1)=tt
     charges_mpi(2)=rholeaked

     call mpiallred(charges_mpi(1),2,MPI_SUM,pkernel%mpi_env%mpi_comm,ierr)

     tt_tot=charges_mpi(1)
     rholeaked_tot=charges_mpi(2)
  else
     tt_tot=tt
     rholeaked_tot=rholeaked
  end if

  if (verb) then
     call yaml_comment('Ionic Potential Creation',hfill='-')
     call yaml_map('Total ionic charge',tt_tot,fmt='(f26.12)')
     if (rholeaked_tot /= 0.0_gp) call yaml_map('Leaked charge',rholeaked_tot,fmt='(1pe10.3)')
     !write(*,'(1x,a)')&
     !     '----------------------------------------------------------- Ionic Potential Creation'
     !write(*,'(1x,a,f26.12,2x,1pe10.3)') 'total ionic charge, leaked charge ',tt_tot,rholeaked_tot
     quiet = "no "
  else
     quiet = "yes"
  end if

  if (.not. htoobig) then
     call timing(iproc,'CrtLocPot     ','OF')
     !here the value of the datacode must be kept fixed
     !n(c) nspin=1

     !if (nproc > 1) call MPI_BARRIER(bigdft_mpi%mpi_env%mpi_comm,ierr)

     call H_potential('D',pkernel,pot_ion,pot_ion,ehart,-psoffset,.false.,quiet=quiet)

     call timing(iproc,'CrtLocPot     ','ON')
     
     if (check_potion) then
        !if (iproc == 0) write(*,'(1x,a)',advance='no') 'Check the ionic potential...'
          
        allocate(potion_corr(n1i*n2i*n3pi+ndebug),stat=i_stat)
        call memocc(i_stat,potion_corr,'potion_corr',subname)

        call razero(n1i*n2i*n3pi,potion_corr)

        !calculate pot_ion with an explicit error function to correct in the case of big grid spacings
        !for the moment works only in the isolated BC case
        do i3=1,n3pi
           z=real(i3+i3s-1-nbl3-1,gp)*hzh
           do i2=1,n2i
              y=real(i2-nbl2-1,gp)*hyh
              do i1=1,n1i
                 x=real(i1-nbl1-1,gp)*hxh
                 ind=i1+(i2-1)*n1i+(i3-1)*n1i*n2i
                 !if (i1==49 .and. i2==46 .and. i3==44) then
                    call sum_erfcr(at%astruct%nat,at%astruct%ntypes,x,y,z,at%astruct%iatype,at%nelpsp,at%psppar,rxyz,potxyz)
                 !   stop
                 !end if
                 potion_corr(ind)=potion_corr(ind)+potxyz
                 !write(18,'(3(i6),i12,3(1x,1pe24.17))')i1,i2,i3,ind,potion_corr(ind),pot_ion(ind)
              end do
           end do
        end do

        !then calculate the maximum difference in the sup norm
        maxdiff=0.0_wp
        do i3=1,n3pi
           do i2=1,n2i
              do i1=1,n1i
                 ind=i1+(i2-1)*n1i+(i3-1)*n1i*n2i
                 maxdiff=max(maxdiff,abs(potion_corr(ind)-pot_ion(ind)))
                 !write(17,'(3(i6),i12,3(1x,1pe24.17))')i1,i2,i3,ind,potion_corr(ind),pot_ion(ind),maxdiff
              end do
           end do
        end do

        call mpiallred(maxdiff,1,MPI_MAX,pkernel%mpi_env%mpi_comm,ierr)

        if (iproc == 0) call yaml_map('Check the ionic potential',maxdiff,fmt='(1pe24.17)')
        !if (iproc == 0) write(*,'(1x,a,1pe24.17)')'...done. MaxDiff=',maxdiff

        stop

        i_all=-product(shape(potion_corr))*kind(potion_corr)
        deallocate(potion_corr,stat=i_stat)
        call memocc(i_stat,i_all,'potion_corr',subname)

     end if

  end if


!!!  !calculate the value of the offset to be put
!!!  tt_tot=0.d0
!!!  do ind=1,n1i*n2i*n3i
!!!     tt_tot=tt_tot+pot_ion(ind)
!!!  end do
!!!  print *,'previous offset',tt_tot*hxh*hyh*hzh

  if (n3pi > 0) then
!    Only for HGH pseudos
     do iat=1,at%astruct%nat
        ityp=at%astruct%iatype(iat)

        rx=rxyz(1,iat)
        ry=rxyz(2,iat)
        rz=rxyz(3,iat)

        rloc=at%psppar(0,0,ityp)
        cutoff=10.d0*rloc

        isx=floor((rx-cutoff)/hxh)
        isy=floor((ry-cutoff)/hyh)
        isz=floor((rz-cutoff)/hzh)

        iex=ceiling((rx+cutoff)/hxh)
        iey=ceiling((ry+cutoff)/hyh)
        iez=ceiling((rz+cutoff)/hzh)
        
        if( at%npspcode(1) .ne.7) then

!          Add the remaining local terms of Eq. (9)
!          in JCP 129, 014109(2008)

           ! determine number of local terms
           nloc=0
           do iloc=1,4
              if (at%psppar(0,iloc,ityp) /= 0.d0) nloc=iloc
           enddo

           !do not add the local part for the vacancy
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

                                !use multipole-preserving function
!!$                                xp=mp_exp(hxh,rx,0.5_gp/(rloc**2),i1,0,.true.)*&
!!$                                     mp_exp(hyh,ry,0.5_gp/(rloc**2),i2,0,.true.)*&
!!$                                     mp_exp(hzh,rz,0.5_gp/(rloc**2),i3,0,.true.)
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
           end if !nloc
        else !HGH or PAW
           ! For PAW, add V^PAW-V_L^HGH
           charge=real(at%nelpsp(ityp),kind=8)
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
                             rr=sqrt(r2)
                             !1) V_L^HGH
                             if(rr>0.01d0) then
                               arg=rr/(sqrt(2.0)*rloc)
                               call derf_ab(tt,arg)
                               raux2=-charge/rr*tt  
                             else
                               !In this case we deduce the values
                               !from a quadratic interpolation (due to 1/rr factor)
                               call interpol_vloc(rr,rloc,charge,raux2)
                             end if
                             !2) V^PAW from splines
                             call splint(rholoc%msz(ityp),rholoc%rad(1:rholoc%msz(ityp),ityp),&
&                              rholoc%d(1:rholoc%msz(ityp),3,ityp),rholoc%d(1:rholoc%msz(ityp),4,ityp),1,rr,raux,ierr)
                             
                             ind=j1+1+nbl1+(j2+nbl2)*n1i+(j3-i3s+1-1)*n1i*n2i
                             pot_ion(ind)=pot_ion(ind)+raux-raux2
                          end if
                       enddo
                    end if
                 enddo
              end if
           end do
        end if ! at%npspcode(iat) .ne.7
     end do !iat
     !debug exit

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
                 call sum_erfcr(at%astruct%nat,at%astruct%ntypes,x,y,z,at%astruct%iatype,at%nelpsp,at%psppar,rxyz,potxyz)
                 pot_ion(ind)=pot_ion(ind)+potxyz
              end do
           end do
        end do
     end if
     
  end if
 

!!!  !calculate the value of the offset to be put
!!!  tt_tot=0.d0
!!!  do ind=1,n1i*n2i*n3i
!!!     tt_tot=tt_tot+pot_ion(ind)
!!!  end do
!!!  print *,'actual offset',tt_tot*hxh*hyh*hzh

  !use rhopot to calculate the potential from a constant electric field along y direction
  if (.not. all(elecfield(1:3) == 0.0_gp)) then
     !constant electric field allowed only for surface and free BC
     if (geocode == 'P') then
     !if (iproc == 0) 
           write(*,'(1x,a)') &
          'The constant electric field is not allowed for Fully Periodic BC.'
          !'The constant electric field is allowed only for Free and Surfaces BC'
     stop
      !constant electric field allowed for surface BC only normal to the surface
     elseif (geocode == 'S' .and. (elecfield(1) /= 0.0_gp .or. elecfield(3) /= 0.0_gp) ) then
     !if (iproc == 0) 
           write(*,'(1x,a)') &
          'Only normal constant electric field (Ex=Ez=0) is allowed for Surface BC.'
     stop
     end if
     if (verb) call yaml_map('Constant electric field (Ha/Bohr)',elecfield(1:3),fmt='(es10.2)')
     !if (verb) write(*,'(1x,a,"(",es10.2,", ",es10.2,", ",es10.2,") ", a)') &
     !     'Constant electric field ',elecfield(1:3),' Ha/Bohr'
!or         'Parabolic confining potential: rprb=',elecfield,&
!           ';  v_conf(r)= 1/(2*rprb**4) * r**2'

     !write or not electric field in a separate file
     efwrite=.false.!true.

     if (n3pi > 0) then
        do i3=1,n3pi
           z=real(i3+i3s-1-nbl3-1,gp)*hzh
           do i2=1,n2i
              y=real(i2-nbl2-1,gp)*hyh
                 do i1=1,n1i
                    x=real(i1-nbl1-1,gp)*hxh
                    ind=i1+(i2-1)*n1i+(i3-1)*n1i*n2i
                    pot_ion(ind)=pot_ion(ind)+elecfield(1)*x+elecfield(2)*y+elecfield(3)*z
!                    parabola: these two lines replace the above line 
!                              comment out the if case and calculate x, z
!                    r2=(x-rx)**2+(y-ry)**2+(z-rz)**2
!                    pot_ion(ind)=pot_ion(ind)+0.5_gp/(elecfield**4)*r2
                 end do
           end do
        end do

        if (efwrite .and. iproc == 0) then
           open(unit=17,file='elecpotential_x',status='unknown')
           write(17,*) "# x , external electric potential(x,y=0,z=0)"
           do i1=nbl1+1,n1i-nbr1-1
              x=real(i1-nbl1-1,gp)*hxh
              write(17,*)x,-elecfield(1)*x
           end do
           close(17)
           open(unit=17,file='elecpotential_y',status='unknown')
           write(17,*) "# y , external electric potential(x=0,y,z=0)"
           do i2=nbl2+1,n2i-nbr2-1
              y=real(i2-nbl2-1,gp)*hyh
              write(17,*)y,-elecfield(2)*y
           end do
           close(17)
           open(unit=17,file='elecpotential_z',status='unknown')
           write(17,*) "# z , external electric potential(x=0,y=0,z)"
           do i3=1,n3pi
              z=real(i3+i3s-1-nbl3-1,gp)*hzh
              write(17,*)z,-elecfield(3)*z
           end do
           close(17)
        end if

     end if
  end if

  !call finalize_real_space_conversion(subname)

  call timing(iproc,'CrtLocPot     ','OF')

contains

! We use a quadratic interpolation to get vloc(x)
! useful for small values of x
  SUBROUTINE interpol_vloc(xx,rloc,charge,yy)
    implicit none
    real(dp),intent(in)  :: xx,rloc,charge
    real(dp),intent(out) :: yy 
!   local variables
    real(dp)::l0,l1,l2,x0,x1,x2,y0,y1,y2

!   Find 3 points (x0,y0), (x1,y1), (x2,y2).
    x0=0.01d0; x1=0.02d0; x2=0.03d0
    call calcVloc(y0,x0,rloc,charge)
    call calcVloc(y1,x1,rloc,charge)
    call calcVloc(y2,x2,rloc,charge)   

!   Find a polynomial of the form:
!   P(x)=y0L0(x) + y1L1(x) + y2L2(x)
 
!   L0(x) = (x-x1)(x-x2)/((x0-x1)(x0-x2))
    l0=(xx-x1)*(xx-x2)/((x0-x1)*(x0-x2))
!   L1(x) = (x-x0)(x-x2)/((x1-x0)(x1-x2))
    l1=(xx-x0)*(xx-x2)/((x1-x0)*(x1-x2))
!   L2(x) = (x-x0)(x-x1)/((x2-x0)(x2-x1))
    l2=(xx-x0)*(xx-x1)/((x2-x0)*(x2-x1))

    yy=y0*l0+y1*l1+y2*l2

  END SUBROUTINE interpol_vloc

  subroutine calcVloc(yy,xx,rloc,Z)
   implicit none
   INTEGER, PARAMETER   :: DP = KIND(1.0D0)          ! double precision
   real(dp),intent(in)  :: xx,rloc,Z
   real(dp),intent(out) :: yy
   real(dp):: arg,tt
  
   arg=xx/(sqrt(2.0)*rloc)
   call derf_ab(tt,arg)
   yy=-Z/xx*tt
  
  
  end subroutine calcVloc


END SUBROUTINE createIonicPotential


!> Determine the index in which the potential must be inserted, following the BC
!! Determine also whether the index is inside or outside the box for free BC
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

END SUBROUTINE ind_positions

!>   Determine the index in which the potential must be inserted, following the BC
!!   Determine also whether the index is inside or outside the box for free BC
subroutine ind_positions_new(periodic,i,ni,j,go)
  implicit none
  logical, intent(in) :: periodic
  integer, intent(in) :: i,ni
  logical, intent(out) :: go
  integer, intent(out) :: j

  if (periodic) then
     go=.true.
     j=modulo(i,ni)
  else
     j=i
     if (i >= -14 .and. i <= ni-15) then
        go=.true.
     else
        go=.false.
     end if
  end if

END SUBROUTINE ind_positions_new


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
  real(wp) :: pi,charge
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
        call derf_ab(derf_val,r/sq2rl)
        potxyz = potxyz - charge*real(derf_val/r,wp)
     end if

     !write(*,'(a,1x,i0,3(1pe24.17))')'iat,r,derf_val,derf_val/r',iat,r/sq2rl,derf_val,derf_val/r

  end do

END SUBROUTINE sum_erfcr


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
END SUBROUTINE ext_buffers


!> Read and initialize counter-ions potentials (read psp files)
subroutine CounterIonPotential(geocode,iproc,nproc,in,shift,&
     hxh,hyh,hzh,grid,n3pi,i3s,pkernel,pot_ion)
  use module_base
  use module_types
  use module_interfaces, except_this_one => CounterIonPotential
  use Poisson_Solver, except_dp => dp, except_gp => gp, except_wp => wp
  use yaml_output
  !use gaussians
  implicit none
  character(len=1), intent(in) :: geocode !< @copydoc poisson_solver::doc::geocode
  integer, intent(in) :: iproc,nproc,n3pi,i3s
  real(gp), intent(in) :: hxh,hyh,hzh
  real(gp), dimension(3), intent(in) :: shift
  type(input_variables), intent(in) :: in
  type(grid_dimensions), intent(in) :: grid
  type(coulomb_operator), intent(in) :: pkernel
  real(wp), dimension(*), intent(inout) :: pot_ion
  !local variables
  character(len=*), parameter :: subname='CounterIonPotential'
  logical :: htoobig=.false.,check_potion=.false.
  logical :: perx,pery,perz,gox,goy,goz
  integer :: iat,i1,i2,i3,j1,j2,j3,isx,isy,isz,iex,iey,iez,ierr,ityp,nspin
  integer :: ind,i_all,i_stat,nbl1,nbr1,nbl2,nbr2,nbl3,nbr3
  real(kind=8) :: pi,rholeaked,rloc,charge,cutoff,x,y,z,r2,arg,xp,tt,rx,ry,rz
  real(kind=8) :: tt_tot,rholeaked_tot,potxyz
  real(wp) :: maxdiff
  real(gp) :: ehart
  type(atoms_data) :: at
  real(dp), dimension(2) :: charges_mpi
  real(dp), dimension(:), allocatable :: potion_corr
  real(gp), dimension(:,:), allocatable :: radii_cf

  call timing(iproc,'CrtLocPot     ','ON')

  !call initialize_real_space_conversion() !initialize the work arrays needed to integrate with isf

  if (iproc.eq.0) then
     write(*,'(1x,a)')&
          '--------------------------------------------------- Counter Ionic Potential Creation'
  end if

  !read the positions of the counter ions from file
  call read_atomic_file('posinp_ci',iproc,at%astruct)
  call allocate_atoms_nat(at, subname)
  call allocate_atoms_ntypes(at, subname)
  ! Read associated pseudo files.
  call init_atomic_values((iproc == 0), at, in%ixc)
  call read_atomic_variables(at, 'input.occup', in%nspin)

  allocate(radii_cf(at%astruct%ntypes,3+ndebug),stat=i_stat)
  call memocc(i_stat,radii_cf,'radii_cf',subname)

  !read the specifications of the counter ions from pseudopotentials
  call read_radii_variables(at, radii_cf, in%crmult, in%frmult, in%projrad)
  if (iproc == 0) call print_atomic_variables(at, radii_cf, max(in%hx,in%hy,in%hz), in%ixc)

  pi=4.d0*atan(1.d0)
  ! Ionic charge (must be calculated for the PS active processes)
  rholeaked=0.d0
  ! Ionic energy (can be calculated for all the processors)

  !Creates charge density arising from the ionic PSP cores
  call razero(grid%n1i*grid%n2i*n3pi,pot_ion)


  !conditions for periodicity in the three directions
  perx=(geocode /= 'F')
  pery=(geocode == 'P')
  perz=(geocode /= 'F')

  call ext_buffers(perx,nbl1,nbr1)
  call ext_buffers(pery,nbl2,nbr2)
  call ext_buffers(perz,nbl3,nbr3)

  if (n3pi >0 .and. .not. htoobig) then

     do iat=1,at%astruct%nat
        ityp=at%astruct%iatype(iat)
        !shift the positions of the counter_ion wrt the box
        rx=at%astruct%rxyz(1,iat)-shift(1)
        ry=at%astruct%rxyz(2,iat)-shift(2)
        rz=at%astruct%rxyz(3,iat)-shift(3)

        if (iproc == 0) then
           write(*,'(1x,a,i6,3(1pe14.7))')'counter ion No. ',iat,rx,ry,rz
        end if

        rloc=at%psppar(0,0,ityp)
        charge=real(at%nelpsp(ityp),kind=8)/(2.d0*pi*sqrt(2.d0*pi)*rloc**3)
        cutoff=10.d0*rloc

        isx=floor((rx-cutoff)/hxh)
        isy=floor((ry-cutoff)/hyh)
        isz=floor((rz-cutoff)/hzh)

        iex=ceiling((rx+cutoff)/hxh)
        iey=ceiling((ry+cutoff)/hyh)
        iez=ceiling((rz+cutoff)/hzh)

        !print *,'rloc,iat,nelpsp',isx,iex,isy,iey,isz,iez,shift(:),iproc

        do i3=isz,iez
           z=real(i3,kind=8)*hzh-rz
           call ind_positions(perz,i3,grid%n3,j3,goz) 
           j3=j3+nbl3+1
           do i2=isy,iey
              y=real(i2,kind=8)*hyh-ry
              call ind_positions(pery,i2,grid%n2,j2,goy)
              do i1=isx,iex
                 x=real(i1,kind=8)*hxh-rx
                 call ind_positions(perx,i1,grid%n1,j1,gox)
                 r2=x**2+y**2+z**2
                 arg=r2/rloc**2

                 !multipole-preserving approach
!!$                 xp=mp_exp(hxh,rx,0.5_gp/(rloc**2),i1,0,.true.)*&
!!$                      mp_exp(hyh,ry,0.5_gp/(rloc**2),i2,0,.true.)*&
!!$                      mp_exp(hzh,rz,0.5_gp/(rloc**2),i3,0,.true.)
                 xp=exp(-.5d0*arg)
                 if (j3 >= i3s .and. j3 <= i3s+n3pi-1  .and. goy  .and. gox ) then
                    ind=j1+1+nbl1+(j2+nbl2)*grid%n1i+(j3-i3s+1-1)*grid%n1i*grid%n2i
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
     do i2= -nbl2,2*grid%n2+1+nbr2
        do i1= -nbl1,2*grid%n1+1+nbr1
           ind=i1+1+nbl1+(i2+nbl2)*grid%n1i+(j3-1)*grid%n1i*grid%n2i
           tt=tt+pot_ion(ind)
        enddo
     enddo
  enddo

  tt=tt*hxh*hyh*hzh
  rholeaked=rholeaked*hxh*hyh*hzh

  if (pkernel%mpi_env%nproc > 1) then
     charges_mpi(1)=tt
     charges_mpi(2)=rholeaked

     call mpiallred(charges_mpi(1),2,MPI_SUM,pkernel%mpi_env%mpi_comm,ierr)

     tt_tot=charges_mpi(1)
     rholeaked_tot=charges_mpi(2)
  else
     tt_tot=tt
     rholeaked_tot=rholeaked
  end if

  if (iproc == 0) write(*,'(1x,a,f26.12,2x,1pe10.3)') &
       'total ionic charge, leaked charge ',tt_tot,rholeaked_tot

  if (.not. htoobig) then
     call timing(iproc,'CrtLocPot     ','OF')
     !here the value of the datacode must be kept fixed
     nspin=1

     call H_potential('D',pkernel,pot_ion,pot_ion,ehart,0.0_gp,.false.)

     call timing(iproc,'CrtLocPot     ','ON')
     
     if (check_potion) then
        !if (iproc == 0) write(*,'(1x,a)',advance='no') 'Check the ionic potential...'
          
        allocate(potion_corr(grid%n1i*grid%n2i*n3pi+ndebug),stat=i_stat)
        call memocc(i_stat,potion_corr,'potion_corr',subname)

        call razero(grid%n1i*grid%n2i*n3pi,potion_corr)

        !calculate pot_ion with an explicit error function to correct in the case of big grid spacings
        !for the moment works only in the isolated BC case
        do i3=1,n3pi
           z=real(i3+i3s-1-nbl3-1,gp)*hzh
           do i2=1,grid%n2i
              y=real(i2-nbl2-1,gp)*hyh
              do i1=1,grid%n1i
                 x=real(i1-nbl1-1,gp)*hxh
                 ind=i1+(i2-1)*grid%n1i+(i3-1)*grid%n1i*grid%n2i
                 !if (i1==49 .and. i2==46 .and. i3==44) then
                    call sum_erfcr(at%astruct%nat,at%astruct%ntypes,x,y,z,at%astruct%iatype,at%nelpsp,at%psppar,&
                          at%astruct%rxyz,potxyz)
                 !   stop
                 !end if
                 potion_corr(ind)=potion_corr(ind)+potxyz
                 !write(18,'(3(i6),i12,3(1x,1pe24.17))')i1,i2,i3,ind,potion_corr(ind),pot_ion(ind)
              end do
           end do
        end do

        !then calculate the maximum difference in the sup norm
        maxdiff=0.0_wp
        do i3=1,n3pi
           do i2=1,grid%n2i
              do i1=1,grid%n1i
                 ind=i1+(i2-1)*grid%n1i+(i3-1)*grid%n1i*grid%n2i
                 maxdiff=max(maxdiff,abs(potion_corr(ind)-pot_ion(ind)))
                 !write(17,'(3(i6),i12,3(1x,1pe24.17))')i1,i2,i3,ind,potion_corr(ind),pot_ion(ind),maxdiff
              end do
           end do
        end do

        call mpiallred(maxdiff,1,MPI_MAX,pkernel%mpi_env%mpi_comm,ierr)

        if (iproc == 0) call yaml_map('Check the ionic potential',maxdiff,fmt='(1pe24.17)')
        !if (iproc == 0) write(*,'(1x,a,1pe24.17)')'...done. MaxDiff=',maxdiff

        stop

        i_all=-product(shape(potion_corr))*kind(potion_corr)
        deallocate(potion_corr,stat=i_stat)
        call memocc(i_stat,i_all,'potion_corr',subname)

     end if

  end if

  if (n3pi > 0 .and. htoobig) then
     !add to pot_ion an explicit error function to correct in the case of big grid spacing
     !for the moment works only in the isolated BC case
     do i3=1,n3pi
        z=real(i3+i3s-1-nbl3-1,gp)*hzh
        do i2=1,grid%n2i
           y=real(i2-nbl2-1,gp)*hyh
           do i1=1,grid%n1i
              x=real(i1-nbl1-1,gp)*hxh
              ind=i1+(i2-1)*grid%n1i+(i3-1)*grid%n1i*grid%n2i
              call sum_erfcr(at%astruct%nat,at%astruct%ntypes,x,y,z,at%astruct%iatype,at%nelpsp,at%psppar,at%astruct%rxyz,potxyz)
              pot_ion(ind)=pot_ion(ind)+potxyz
           end do
        end do
     end do
  end if

  !deallocations
  call deallocate_atoms(at,subname) 

  i_all=-product(shape(radii_cf))*kind(radii_cf)
  deallocate(radii_cf,stat=i_stat)
  call memocc(i_stat,i_all,'radii_cf',subname)

  i_all=-product(shape(at%astruct%rxyz))*kind(at%astruct%rxyz)
  deallocate(at%astruct%rxyz,stat=i_stat)
  call memocc(i_stat,i_all,'at%astruct%rxyz',subname)


  !call finalize_real_space_conversion(subname)
  call timing(iproc,'CrtLocPot     ','OF')

END SUBROUTINE CounterIonPotential
