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
  use module_base, pi => pi_param
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
  type(coulomb_operator), intent(inout) :: pkernel
  real(gp), intent(out) :: eion,edisp,psoffset
  real(dp), dimension(6),intent(out) :: ewaldstr
  real(gp), dimension(:,:), pointer :: fion,fdisp
  real(dp), dimension(*), intent(out) :: pot_ion
  !local variables
  character(len=*), parameter :: subname='IonicEnergyandForces'
  logical :: slowion=.false.
  logical :: perx,pery,perz,gox,goy,goz
  integer :: n1i,n2i,n3i,i3s,n3pi
  integer :: i,iat,ii,ityp,jat,jtyp,nbl1,nbr1,nbl2,nbr2,nbl3,nbr3
  integer :: isx,iex,isy,iey,isz,iez,i1,i2,i3,j1,j2,j3,ind
  real(gp) :: ucvol,rloc,twopitothreehalf,atint,shortlength,charge,eself,rx,ry,rz
  real(gp) :: fxion,fyion,fzion,dist,fxerf,fyerf,fzerf,cutoff
  real(gp) :: hxh,hyh,hzh
  real(gp) :: hxx,hxy,hxz,hyy,hyz,hzz,chgprod
  real(gp) :: x,y,z,xp,Vel,prefactor,r2,arg,ehart,de
  !real(gp) :: Mz,cmassy
  real(gp), dimension(3,3) :: gmet,rmet,rprimd,gprimd
  !other arrays for the ewald treatment
  real(gp), dimension(:,:), allocatable :: fewald,xred
  real(gp), dimension(3) :: cc

  call timing(iproc,'ionic_energy','ON')
  fion = f_malloc_ptr((/ 3, at%astruct%nat /),id='fion')
  fdisp = f_malloc_ptr((/ 3, at%astruct%nat /),id='fdisp')

  ! Aliasing
  hxh = dpbox%hgrids(1)
  hyh = dpbox%hgrids(2)
  hzh = dpbox%hgrids(3)
  n1i = dpbox%ndims(1)
  n2i = dpbox%ndims(2)
  n3i = dpbox%ndims(3)
  i3s = dpbox%i3s+dpbox%i3xcsh
  n3pi = dpbox%n3pi

  psoffset=0.0_gp
  ewaldstr=0.0_gp
  if (at%astruct%geocode == 'P') then
     !here we insert the calculation of the ewald forces
     fewald = f_malloc((/ 3, at%astruct%nat /),id='fewald')
     xred = f_malloc((/ 3, at%astruct%nat /),id='xred')

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
     call ewald(iproc,nproc,eion,gmet,fewald,at%astruct%nat,at%astruct%ntypes,rmet,at%astruct%iatype,ucvol,&
          xred,real(at%nelpsp,kind=8))
     ewaldstr=0.0_dp
     call ewald2(iproc,nproc,gmet,at%astruct%nat,at%astruct%ntypes,rmet,rprimd,ewaldstr,at%astruct%iatype,&
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

     call f_free(xred)
     call f_free(fewald)

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
     if (at%astruct%sym%symObj >= 0) call symm_stress(ewaldstr,at%astruct%sym%symObj)
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

     !LR: commented hessian as not currently using it

     !$omp parallel default(none) &
     !$omp private(iat,ityp,rx,ry,rz,fxion,fyion,fzion,jtyp,chgprod,dist) &
     !$omp shared(at,rxyz,fion,eself,eion)
     !$omp do reduction(+:eself,eion)
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
        !hxx=0.0_gp
        !hxy=0.0_gp
        !hxz=0.0_gp
        !hyy=0.0_gp
        !hyz=0.0_gp
        !hzz=0.0_gp

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
           !hxx=hxx+3.0_gp*chgprod/(dist**5)*(rx-rxyz(1,jat))**2-chgprod/(dist**3)
           !hxy=hxy+3.0_gp*chgprod/(dist**5)*(rx-rxyz(1,jat))*(ry-rxyz(2,jat))
           !hxz=hxz+3.0_gp*chgprod/(dist**5)*(rx-rxyz(1,jat))*(rz-rxyz(3,jat))
           !hyy=hyy+3.0_gp*chgprod/(dist**5)*(ry-rxyz(2,jat))**2-chgprod/(dist**3)
           !hyz=hyz+3.0_gp*chgprod/(dist**5)*(ry-rxyz(2,jat))*(rz-rxyz(3,jat))
           !hzz=hzz+3.0_gp*chgprod/(dist**5)*(rz-rxyz(3,jat))**2-chgprod/(dist**3)
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
           !hxx=hxx+3.0_gp*chgprod/(dist**5)*(rx-rxyz(1,jat))**2-chgprod/(dist**3)
           !hxy=hxy+3.0_gp*chgprod/(dist**5)*(rx-rxyz(1,jat))*(ry-rxyz(2,jat))
           !hxz=hxz+3.0_gp*chgprod/(dist**5)*(rx-rxyz(1,jat))*(rz-rxyz(3,jat))
           !hyy=hyy+3.0_gp*chgprod/(dist**5)*(ry-rxyz(2,jat))**2-chgprod/(dist**3)
           !hyz=hyz+3.0_gp*chgprod/(dist**5)*(ry-rxyz(2,jat))*(rz-rxyz(3,jat))
           !hzz=hzz+3.0_gp*chgprod/(dist**5)*(rz-rxyz(3,jat))**2-chgprod/(dist**3)
        end do

        fion(1,iat)=fxion
        fion(2,iat)=fyion
        fion(3,iat)=fzion

        !if (nproc==1 .and. slowion) print *,'iat,fion',iat,(fion(j1,iat),j1=1,3)
        !energy which comes from the self-interaction of the spread charge
       eself=eself+real(at%nelpsp(ityp)**2,gp)*0.5_gp*sqrt(1.d0/pi)/at%psppar(0,0,ityp)
     end do
     !$omp end do
     !$omp end parallel

     !if (nproc==1 .and. slowion) print *,'eself',eself

  end if

  !for the surfaces BC,
  !activate for the moment only the slow calculation of the ionic energy and forces
  !the slowion command has also to be activated for the cavity calculation
  !if (at%astruct%geocode == 'S' .or. at%astruct%geocode == 'P') slowion=.true.
  if (at%astruct%geocode == 'S' .or. pkernel%method /= 'VAC') slowion=.true.

   slowion_if: if (slowion) then

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

        !calculate the self energy of the isolated bc
        eself=eself+real(at%nelpsp(ityp),gp)**2/rloc
     enddo

     eself=0.5_gp/sqrt(pi)*eself

     !if (nproc==1) 
     !print *,'iproc,eself',iproc,eself
     call f_zero(n1i*n2i*n3pi,pot_ion(1))

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
                       ind=j1+1+nbl1+(j2+nbl2)*n1i+(j3-i3s)*n1i*n2i
                       pot_ion(ind)=pot_ion(ind)-xp*charge
                    endif
                 enddo
              enddo
           enddo

        enddo

     end if

  end if slowion_if

  !in the case of cavity the ionic energy is only considered as the self energy
  nocavity_if: if (pkernel%method /= 'VAC') then
     eion=-eself
  else if (slowion) then
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
                       ind=j1+1+nbl1+(j2+nbl2)*n1i+(j3-i3s)*n1i*n2i
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
        call mpiallred(fion,MPI_SUM,comm=pkernel%mpi_env%mpi_comm)
     end if

     !if (iproc ==0) print *,'eion',eion,psoffset,shortlength

  end if nocavity_if

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

  call timing(iproc,'ionic_energy','OF')

END SUBROUTINE IonicEnergyandForces

!> calculates the value of the dielectric funnction for a smoothed cavity 
!! given a set of centres and radii.
!! Need the epsilon0 as well as the radius of the cavit and its smoothness
subroutine epsilon_rigid_cavity(geocode,ndims,hgrids,nat,rxyz,radii,epsilon0,delta,eps)
  use f_utils
  implicit none
  character(len=1), intent(in) :: geocode !< @copydoc poisson_solver::doc::geocode
  integer, intent(in) :: nat !< number of centres defining the cavity
  real(kind=8), intent(in) :: epsilon0 !< dielectric constant of th solvent
  real(kind=8), intent(in) :: delta !< smoothness factor of the cavity
  integer, dimension(3), intent(in) :: ndims   !< dimensions of the simulation box
  real(kind=8), dimension(3), intent(in) :: hgrids !< grid spacings
  real(kind=8), dimension(nat), intent(in) :: radii !< radii of each of the atoms
  !> position of all the atoms in the grid coordinates
  real(kind=8), dimension(3,nat), intent(in) :: rxyz
  real(kind=8), dimension(ndims(1),ndims(2),ndims(3)), intent(out) :: eps !< dielectric function
  !local variables
  logical :: perx,pery,perz
  integer :: i1,i2,i3,iat,nbl1,nbl2,nbl3,nbr1,nbr2,nbr3,unt
  real(kind=8) :: r2,x,y2,z2,d2,y,z,eps_min,eps1
  !  real(kind=8), dimension(3) :: deps

  !buffers associated to the geocode
  !conditions for periodicity in the three directions
  perx=(geocode /= 'F')
  pery=(geocode == 'P')
  perz=(geocode /= 'F')

  call ext_buffers(perx,nbl1,nbr1)
  call ext_buffers(pery,nbl2,nbr2)
  call ext_buffers(perz,nbl3,nbr3)


  do i3=1,ndims(3)
     z=hgrids(3)*(i3-1-nbl3)
     z2=z*z
     do i2=1,ndims(2)
        y=hgrids(2)*(i2-1-nbl2)
        y2=y*y
        do i1=1,ndims(1)
           x=hgrids(1)*(i1-1-nbl1)
           r2=x*x+y2+z2
           !choose the closest atom
           eps_min=1.d100
           do iat=1,nat
              d2=(x-rxyz(1,iat))**2+(y-rxyz(2,iat))**2+(z-rxyz(3,iat))**2
              if (d2.eq.0.d0) d2=1.0d-30
              eps1=epsl(sqrt(d2),radii(iat),delta,epsilon0)
              if (eps1< eps_min) then
                 !deps(1)=depsoeps(sqrt(d2),radii(iat),delta,epsilon0)*(x-rxyz(1,iat))/sqrt(d2)
                 !deps(2)=depsoeps(sqrt(d2),radii(iat),delta,epsilon0)*(y-rxyz(2,iat))/sqrt(d2)
                 !deps(3)=depsoeps(sqrt(d2),radii(iat),delta,epsilon0)*(z-rxyz(3,iat))/sqrt(d2)
                 eps_min=eps1
              end if
              if (abs(eps_min-1.d0) < epsilon(1.d0)) exit
           end do
           if (nat==0) then
              eps_min=1.d0
              !deps=0.d0
           end if
           eps(i1,i2,i3)=eps_min
           !dlogeps(1:3,i1,i2,i3)=deps(1:3)
        end do
     end do
  end do

  unt=f_get_free_unit(21)
  call f_open_file(unt,file='epsilon.dat')
  i1=1!n03/2
  do i2=1,ndims(2)
     do i3=1,ndims(3)
        write(unt,'(2(1x,I4),2(1x,e14.7))')i2,i3,eps(i1,i2,i3),eps(ndims(1)/2,i2,i3)
     end do
  end do
  call f_close(unt)

  unt=f_get_free_unit(22)
  call f_open_file(unt,file='epsilon_line.dat')
  do i2=1,ndims(2)
   write(unt,'(1x,I8,1(1x,e22.15))')i2,eps(ndims(1)/2,i2,ndims(3)/2)
  end do
  call f_close(unt)

  contains

    pure function epsl(r,rc,delta,epsilon0)
      implicit none
      real(kind=8), intent(in) :: r,rc,delta,epsilon0
      real(kind=8) :: epsl
      !local variables
      real(kind=8) :: d

      d=(r-rc)/delta
      epsl=0.5d0*(epsilon0-1.d0)*(erf(d)+1.d0)+1.d0
    end function epsl

    pure function depsoeps(r,rc,delta,epsilon0)
      implicit none
      real(kind=8), intent(in) :: r,rc,delta,epsilon0
      real(kind=8) :: depsoeps
      !local variables
      real(kind=8) :: d

      d=(r-rc)/delta
      depsoeps=(epsilon0-1.d0)/delta*exp(-d**2)/epsl(r,rc,delta,epsilon0)
    end function depsoeps
 end subroutine epsilon_rigid_cavity

!> calculates the value of the dielectric function for a smoothed cavity 
!! given a set of centres and radii. Based on error function.
!! Need the epsilon0 as well as the radius of the cavit and its smoothness
subroutine epsilon_rigid_cavity_error_multiatoms_bc(geocode,ndims,hgrids,natreal,rxyzreal,radiireal,epsilon0,delta,&
     eps,dlogeps,oneoeps,oneosqrteps,corr,IntSur,IntVol)
  use f_utils
  use module_defs, only : Bohr_Ang,bigdft_mpi
  use f_enums
  use yaml_output
  use dynamic_memory

  implicit none
  character(len=1), intent(in) :: geocode !< @copydoc poisson_solver::doc::geocode
  integer, intent(in) :: natreal !< number of centres defining the cavity
  real(kind=8), intent(in) :: epsilon0 !< dielectric constant of the solvent
  real(kind=8), intent(in) :: delta !< smoothness factor of the cavity
  integer, dimension(3), intent(in) :: ndims   !< dimensions of the simulation box
  real(kind=8), dimension(3), intent(in) :: hgrids !< grid spacings
  real(kind=8), dimension(natreal), intent(in) :: radiireal !< radii of each of the atoms
  !> position of all the atoms in the grid coordinates
  real(kind=8), dimension(3,natreal), intent(in) :: rxyzreal
  real(kind=8), dimension(ndims(1),ndims(2),ndims(3)), intent(out) :: eps !< dielectric function
  real(kind=8), dimension(3,ndims(1),ndims(2),ndims(3)), intent(out) :: dlogeps !< dlogeps
  real(kind=8), dimension(ndims(1),ndims(2),ndims(3)), intent(out) :: oneoeps !< inverse of epsilon. Needed for PI method.
  !> inverse square root of epsilon. Needed for PCG method.
  real(kind=8), dimension(ndims(1),ndims(2),ndims(3)), intent(out) :: oneosqrteps
  real(kind=8), dimension(ndims(1),ndims(2),ndims(3)), intent(out) :: corr !< correction term of the Generalized Laplacian.
  !> Surface and volume integral needed for non-electrostatic contributions to the energy.
  real(kind=8), intent(out) :: IntSur,IntVol 

  !local variables
  logical :: perx,pery,perz
  integer :: i,i1,i2,i3,iat,jat,ii,nat,j,k,l,px,py,pz,unt
  integer :: nbl1,nbl2,nbl3,nbr1,nbr2,nbr3,imin
  real(kind=8) :: r2,x,y2,z2,d,d2,d12,y,z,eps_min,eps1,pi,de2,dde,d1,oneod,h,coeff,dmin,dmax,oneoeps0,oneosqrteps0
  real(kind=8) :: r,t,fact1,fact2,fact3,dd,dtx,curr,value,valuemin
  real(kind=8), dimension(3) :: deps,ddeps,v,rv,shift,sh

  real(kind=8), parameter :: valuebc=1.d0
  real(kind=8), dimension(6) :: plandist
  integer, dimension(6) :: ba
  real(kind=8), dimension(:), allocatable :: ep,ddep
  real(kind=8), dimension(:,:), allocatable :: dep
  real(kind=8), dimension(3,27*natreal) :: rxyztot
  real(kind=8), dimension(27*natreal) :: radiitot
  real(kind=8), dimension(:), allocatable :: radii
  real(kind=8), dimension(:,:), allocatable :: rxyz


  !buffers associated to the geocode
  !conditions for periodicity in the three directions
  perx=(geocode /= 'F')
  pery=(geocode == 'P')
  perz=(geocode /= 'F')

  call ext_buffers(perx,nbl1,nbr1)
  call ext_buffers(pery,nbl2,nbr2)
  call ext_buffers(perz,nbl3,nbr3)

  IntSur=0.d0
  IntVol=0.d0
  pi = 4.d0*datan(1.d0)
  r=0.d0
  t=0.d0
  oneoeps0=1.d0/epsilon0
  oneosqrteps0=1.d0/dsqrt(epsilon0)

  shift(1)=hgrids(1)*ndims(1)
  shift(2)=hgrids(2)*ndims(2)
  shift(3)=hgrids(3)*ndims(3)

!------------------------------------------------------------------------------------------------------
! Depending of Free, Periodic or Surface bc, image atoms are or not included.

  if (bigdft_mpi%iproc==0) then
   do iat=1,natreal
    call yaml_map('real input atoms',iat)
    call yaml_map('radii',radiireal(iat))
    call yaml_map('rxyz',rxyzreal(:,iat))
   end do
  end if

  px=0
  py=0
  pz=0
  if (perx) px=1
  if (pery) py=1
  if (perz) pz=1

  rxyztot(:,:)=0.d0

  i=0
  do iat=1,natreal
   ba(1:6)=0
   ! checking what are the image atoms to include in the calculation of the
   ! cavity.
   rv(1:3)=rxyzreal(1:3,iat)
   plandist(1)=dabs(rv(1))
   plandist(2)=dabs(shift(1)-rv(1))
   plandist(3)=dabs(rv(2))
   plandist(4)=dabs(shift(2)-rv(2))
   plandist(5)=dabs(rv(3))
   plandist(6)=dabs(shift(3)-rv(3))
   do ii=1,6
    valuemin=1.d0
    d=plandist(ii)
    value=epsl(d,radiireal(iat),delta)
    if (value.lt.valuebc) then ! valuebc is the value to check on the box border to accept or refuse an image atom.
     if (abs(value).lt.valuemin) then
      valuemin=abs(value)
      imin=ii
     end if
     select case(ii)
     case (1)
      ba(1)=1*px
     case (2)
      ba(2)=-1*px
     case (3)
      ba(3)=1*py
     case (4)
      ba(4)=-1*py
     case (5)
      ba(5)=1*pz
     case (6)
      ba(6)=-1*pz
     end select
    end if
   end do

   do j=ba(6),ba(5)
    sh(3)=real(j,kind=8)*shift(3)
    do k=ba(4),ba(3)
     sh(2)=real(k,kind=8)*shift(2)
     do l=ba(2),ba(1)
      sh(1)=real(l,kind=8)*shift(1)
      rv(1:3)=rxyzreal(1:3,iat) + sh(1:3)
      i=i+1
      rxyztot(1:3,i)=rv(1:3)
      radiitot(i)=radiireal(iat)
     end do
    end do
   end do

  end do

   nat=i

   ep=f_malloc(nat,id='ep')
   ddep=f_malloc(nat,id='ddep')
   dep=f_malloc([3,nat],id='dep')
   rxyz=f_malloc([3,nat],id='rxyz')
   radii=f_malloc(nat,id='radii')

   rxyz(1:3,1:nat)=rxyztot(1:3,1:nat)
   radii(1:nat)=radiitot(1:nat)

   if (bigdft_mpi%iproc==0) then
    write(*,*)plandist
    write(*,'(1x,a,1x,e14.7,1x,a,1x,i4)')'Value min =',valuemin,'at bc side',imin
    call yaml_map('nat',nat)
    do iat=1,nat
     call yaml_map('atom',iat)
     call yaml_map('radii',radii(iat))
     call yaml_map('rxyz',rxyz(:,iat))
    end do
   end if

!------------------------------------------------------------------------------------------------------
! Starting the cavity building for rxyztot atoms=real+image atoms (total natcurr) for periodic
! and surface boundary conditions or atoms=real for free bc.

  do i3=1,ndims(3)
   z=hgrids(3)*(i3-1-nbl3)
   v(3)=z
   do i2=1,ndims(2)
    y=hgrids(2)*(i2-1-nbl2)
    v(2)=y
    do i1=1,ndims(1)
     x=hgrids(1)*(i1-1-nbl1)
     v(1)=x

     do iat=1,nat
      d2=(x-rxyz(1,iat))**2+(y-rxyz(2,iat))**2+(z-rxyz(3,iat))**2
      d=dsqrt(d2)

      if (d2.eq.0.d0) then
       d2=1.0d-30
       ep(iat)=epsl(d,radii(iat),delta)
       do i=1,3
        dep(i,iat)=0.d0
       end do
       ddep(iat)=0.d0
      else
       oneod=1.d0/d
       ep(iat)=epsl(d,radii(iat),delta)
       d1=d1eps(d,radii(iat),delta)
       coeff=2.d0*((sqrt(d2)-radii(iat))/(delta**2))
       do i=1,3
        h=(v(i)-rxyz(i,iat))*oneod
        dep(i,iat) =d1*h
       end do
       ddep(iat)=d1*(2.d0*oneod-coeff)
      end if

     end do

     IntVol = IntVol + (1.d0-product(ep))

     eps(i1,i2,i3)=(epsilon0-1.d0)*product(ep)+1.d0
     oneoeps(i1,i2,i3)=1.d0/eps(i1,i2,i3)
     oneosqrteps(i1,i2,i3)=1.d0/dsqrt(eps(i1,i2,i3))

     do i=1,3
      deps(i)=0.d0
      do jat=0,nat-1
       curr=dep(i,jat+1)
       do iat=1,nat-1
        curr=curr*ep(modulo(iat+jat,nat)+1)
       end do
        deps(i) = deps(i) + curr
      end do
      deps(i) = deps(i)*(epsilon0-1.d0)
     end do

     d12=0.d0
     do i=1,3
      dlogeps(i,i1,i2,i3)=deps(i)/eps(i1,i2,i3)
      d12 = d12 + deps(i)**2
     end do

     IntSur = IntSur + dsqrt(d12)

     dd=0.d0
     do jat=1,nat
      curr=ddep(jat)
      do iat=1,nat-1
       curr=curr*ep(modulo(iat+jat-1,nat)+1)
      end do
      dd = dd + curr
     end do

      do i=1,3
       do iat=1,nat-1
        do jat=iat+1,nat
         curr=dep(i,iat)*dep(i,jat)
         do ii=1,nat
          if ((ii.eq.iat).or.(ii.eq.jat)) then
          else
           curr=curr*ep(ii)
          end if
         end do
         curr=curr*2.d0
         dd = dd + curr
        end do
       end do
      end do

     dd=dd*(epsilon0-1.d0)
     corr(i1,i2,i3)=(-0.125d0/pi)*(0.5d0*d12/eps(i1,i2,i3)-dd)

    end do
   end do
  end do

  IntSur=IntSur*hgrids(1)*hgrids(2)*hgrids(3)/(epsilon0-1.d0)
  IntVol=IntVol*hgrids(1)*hgrids(2)*hgrids(3)

  unt=f_get_free_unit(21)
  call f_open_file(unt,file='epsilon.dat')
  i1=1!n03/2
  do i2=1,ndims(2)
     do i3=1,ndims(3)
        write(unt,'(2(1x,I4),2(1x,e14.7))')i2,i3,eps(i1,i2,i3),eps(ndims(1)/2,i2,i3)
     end do
  end do
  call f_close(unt)

  unt=f_get_free_unit(22)
  call f_open_file(unt,file='epsilon_line.dat')
  do i2=1,ndims(2)
   write(unt,'(1x,I8,1(1x,e22.15))')i2,eps(ndims(1)/2,i2,ndims(3)/2)
  end do
  call f_close(unt)

  call f_free(ep)
  call f_free(ddep)
  call f_free(dep)
  call f_free(rxyz)
  call f_free(radii)

  contains
    pure function epsl(r,rc,delta)
      implicit none
      real(kind=8), intent(in) :: r,rc,delta
      real(kind=8) :: epsl
      !local variables
      real(kind=8) :: d

      d=(r-rc)/delta
      epsl=0.5d0*(erf(d)+1.d0)
    end function epsl

    pure function d1eps(r,rc,delta)
      use numerics, only: safe_exp
      implicit none
      real(kind=8), intent(in) :: r,rc,delta
      real(kind=8) :: d1eps
      !local variables
      real(kind=8) :: d

      d=(r-rc)/delta
      d1eps=(1.d0/(delta*sqrt(pi)))*max(safe_exp(-d**2),1.0d-24)
    end function d1eps

end subroutine epsilon_rigid_cavity_error_multiatoms_bc

!> calculates the value of the dielectric function for a smoothed cavity 
!! given a set of centres and radii. Based on error function.
!! Need the epsilon0 as well as the radius of the cavit and its smoothness
subroutine epsilon_rigid_cavity_error_multiatoms(geocode,ndims,hgrids,nat,rxyz,radii,epsilon0,delta,&
     eps,dlogeps,oneoeps,oneosqrteps,corr,IntSur,IntVol)
  use f_utils
  use module_defs, only : Bohr_Ang,bigdft_mpi
  use f_enums
  use yaml_output

  implicit none
  character(len=1), intent(in) :: geocode !< @copydoc poisson_solver::doc::geocode
  integer, intent(in) :: nat !< number of centres defining the cavity
  real(kind=8), intent(in) :: epsilon0 !< dielectric constant of the solvent
  real(kind=8), intent(in) :: delta !< smoothness factor of the cavity
  integer, dimension(3), intent(in) :: ndims   !< dimensions of the simulation box
  real(kind=8), dimension(3), intent(in) :: hgrids !< grid spacings
  real(kind=8), dimension(nat), intent(in) :: radii !< radii of each of the atoms
  !> position of all the atoms in the grid coordinates
  real(kind=8), dimension(3,nat), intent(in) :: rxyz
  real(kind=8), dimension(ndims(1),ndims(2),ndims(3)), intent(out) :: eps !< dielectric function
  real(kind=8), dimension(3,ndims(1),ndims(2),ndims(3)), intent(out) :: dlogeps !< dlogeps
  real(kind=8), dimension(ndims(1),ndims(2),ndims(3)), intent(out) :: oneoeps !< inverse of epsilon. Needed for PI method.
  !> inverse square root of epsilon. Needed for PCG method.
  real(kind=8), dimension(ndims(1),ndims(2),ndims(3)), intent(out) :: oneosqrteps
  real(kind=8), dimension(ndims(1),ndims(2),ndims(3)), intent(out) :: corr !< correction term of the Generalized Laplacian.
  !> Surface and volume integral needed for non-electrostatic contributions to the energy.
  real(kind=8), intent(out) :: IntSur,IntVol 

  !local variables
  logical :: perx,pery,perz
  integer :: i,i1,i2,i3,iat,jat,ii,nbl1,nbl2,nbl3,nbr1,nbr2,nbr3,unt
  real(kind=8) :: r2,x,y2,z2,d,d2,d12,y,z,eps_min,eps1,pi,de2,dde,d1,oneod,h,coeff,dmin,dmax,oneoeps0,oneosqrteps0
  real(kind=8) :: r,t,fact1,fact2,fact3,dd,dtx,curr
  real(kind=8), dimension(3) :: deps,ddeps,v
  real(kind=8), dimension(nat) :: ep,ddep
  real(kind=8), dimension(3,nat) :: dep

  !buffers associated to the geocode
  !conditions for periodicity in the three directions
  perx=(geocode /= 'F')
  pery=(geocode == 'P')
  perz=(geocode /= 'F')

  call ext_buffers(perx,nbl1,nbr1)
  call ext_buffers(pery,nbl2,nbr2)
  call ext_buffers(perz,nbl3,nbr3)

  IntSur=0.d0
  IntVol=0.d0
  pi = 4.d0*datan(1.d0)
  r=0.d0
  t=0.d0
  oneoeps0=1.d0/epsilon0
  oneosqrteps0=1.d0/dsqrt(epsilon0)

  if (bigdft_mpi%iproc==0) then
   call yaml_map('nbl1',nbl1)
   call yaml_map('nbl2',nbl2)
   call yaml_map('nbl3',nbl3)
   do iat=1,nat
    call yaml_map('atom',iat)
    call yaml_map('rxyz',rxyz(:,iat))
   end do
  end if

  do i3=1,ndims(3)
   z=hgrids(3)*(i3-1-nbl3)
   v(3)=z
   do i2=1,ndims(2)
    y=hgrids(2)*(i2-1-nbl2)
    v(2)=y
    do i1=1,ndims(1)
     x=hgrids(1)*(i1-1-nbl1)
     v(1)=x

     do iat=1,nat
      d2=(x-rxyz(1,iat))**2+(y-rxyz(2,iat))**2+(z-rxyz(3,iat))**2
      d=dsqrt(d2)

      if (d2.eq.0.d0) then
       d2=1.0d-30
       ep(iat)=epsl(d,radii(iat),delta)
       do i=1,3
        dep(i,iat)=0.d0
       end do
       ddep(iat)=0.d0
      else
       oneod=1.d0/d
       ep(iat)=epsl(d,radii(iat),delta)
       d1=d1eps(d,radii(iat),delta)
       coeff=2.d0*((sqrt(d2)-radii(iat))/(delta**2))
       do i=1,3
        h=(v(i)-rxyz(i,iat))*oneod
        dep(i,iat) =d1*h
       end do
       ddep(iat)=d1*(2.d0*oneod-coeff)
      end if

     end do

     IntVol = IntVol + (1.d0-product(ep))

     eps(i1,i2,i3)=(epsilon0-1.d0)*product(ep)+1.d0
     oneoeps(i1,i2,i3)=1.d0/eps(i1,i2,i3)
     oneosqrteps(i1,i2,i3)=1.d0/dsqrt(eps(i1,i2,i3))

     do i=1,3
      deps(i)=0.d0
      do jat=0,nat-1
       curr=dep(i,jat+1)
       do iat=1,nat-1
        curr=curr*ep(modulo(iat+jat,nat)+1)
       end do
        deps(i) = deps(i) + curr
      end do
      deps(i) = deps(i)*(epsilon0-1.d0)
     end do

     d12=0.d0
     do i=1,3
      dlogeps(i,i1,i2,i3)=deps(i)/eps(i1,i2,i3)
      d12 = d12 + deps(i)**2
     end do

     IntSur = IntSur + dsqrt(d12)

     dd=0.d0
     do jat=1,nat
      curr=ddep(jat)
      do iat=1,nat-1
       curr=curr*ep(modulo(iat+jat-1,nat)+1)
      end do
      dd = dd + curr
     end do

      do i=1,3
       do iat=1,nat-1
        do jat=iat+1,nat
         curr=dep(i,iat)*dep(i,jat)
         do ii=1,nat
          if ((ii.eq.iat).or.(ii.eq.jat)) then
          else
           curr=curr*ep(ii)
          end if
         end do
         curr=curr*2.d0
         dd = dd + curr
        end do
       end do
      end do

     dd=dd*(epsilon0-1.d0)
     corr(i1,i2,i3)=(-0.125d0/pi)*(0.5d0*d12/eps(i1,i2,i3)-dd)

    end do
   end do
  end do

  IntSur=IntSur*hgrids(1)*hgrids(2)*hgrids(3)/(epsilon0-1.d0)
  IntVol=IntVol*hgrids(1)*hgrids(2)*hgrids(3)

  unt=f_get_free_unit(21)
  call f_open_file(unt,file='epsilon.dat')
  i1=1!n03/2
  do i2=1,ndims(2)
     do i3=1,ndims(3)
        write(unt,'(2(1x,I4),2(1x,e14.7))')i2,i3,eps(i1,i2,i3),eps(ndims(1)/2,i2,i3)
     end do
  end do
  call f_close(unt)

  unt=f_get_free_unit(22)
  call f_open_file(unt,file='epsilon_line.dat')
  do i2=1,ndims(2)
   write(unt,'(1x,I8,1(1x,e22.15))')i2,eps(ndims(1)/2,i2,ndims(3)/2)
  end do
  call f_close(unt)

  contains
    pure function epsl(r,rc,delta)
      implicit none
      real(kind=8), intent(in) :: r,rc,delta
      real(kind=8) :: epsl
      !local variables
      real(kind=8) :: d

      d=(r-rc)/delta
      epsl=0.5d0*(erf(d)+1.d0)
    end function epsl

    pure function d1eps(r,rc,delta)
      use numerics, only: safe_exp
      implicit none
      real(kind=8), intent(in) :: r,rc,delta
      real(kind=8) :: d1eps
      !local variables
      real(kind=8) :: d

      d=(r-rc)/delta
      d1eps=(1.d0/(delta*sqrt(pi)))*max(safe_exp(-d**2),1.0d-24)
    end function d1eps

end subroutine epsilon_rigid_cavity_error_multiatoms

!> calculates the value of the dielectric function for a smoothed cavity 
!! given a set of centres and radii. Based on the Andreussi epsilon function
!! with a gaussian \rho^{elec}.
!! Need the epsilon0 as well as the radius of the cavit and its smoothness
subroutine epsilon_rigid_cavity_new_multiatoms(geocode,ndims,hgrids,nat,rxyz,radii,epsilon0,delta,&
     eps,dlogeps,oneoeps,oneosqrteps,corr,IntSur,IntVol)
  use f_utils
  implicit none
  character(len=1), intent(in) :: geocode !< @copydoc poisson_solver::doc::geocode
  integer, intent(in) :: nat !< number of centres defining the cavity
  real(kind=8), intent(in) :: epsilon0 !< dielectric constant of the solvent
  real(kind=8), intent(in) :: delta !< smoothness factor of the cavity
  integer, dimension(3), intent(in) :: ndims   !< dimensions of the simulation box
  real(kind=8), dimension(3), intent(in) :: hgrids !< grid spacings
  real(kind=8), dimension(nat), intent(in) :: radii !< radii of each of the atoms
  !> position of all the atoms in the grid coordinates
  real(kind=8), dimension(3,nat), intent(in) :: rxyz
  real(kind=8), dimension(ndims(1),ndims(2),ndims(3)), intent(out) :: eps !< dielectric function
  real(kind=8), dimension(3,ndims(1),ndims(2),ndims(3)), intent(out) :: dlogeps !< dlogeps
  real(kind=8), dimension(ndims(1),ndims(2),ndims(3)), intent(out) :: oneoeps !< inverse of epsilon. Needed for PI method.
  !> inverse square root of epsilon. Needed for PCG method.
  real(kind=8), dimension(ndims(1),ndims(2),ndims(3)), intent(out) :: oneosqrteps 
  real(kind=8), dimension(ndims(1),ndims(2),ndims(3)), intent(out) :: corr !< correction term of the Generalized Laplacian.
  !> Surface and volume integral needed for non-electrostatic contributions to the energy.
  real(kind=8), intent(out) :: IntSur,IntVol 

  !local variables
  logical :: perx,pery,perz
  integer :: i,i1,i2,i3,iat,jat,ii,nbl1,nbl2,nbl3,nbr1,nbr2,nbr3,unt
  real(kind=8) :: r2,x,y2,z2,d,d2,d12,y,z,eps_min,eps1,pi,de2,dde,d1,oneod,h,coeff,dmin,dmax,oneoeps0,oneosqrteps0
  real(kind=8) :: r,t,fact1,fact2,fact3,dd,dtx,curr
  real(kind=8), dimension(3) :: deps,ddeps,v
  real(kind=8), dimension(nat) :: ep,ddep
  real(kind=8), dimension(3,nat) :: dep

  !buffers associated to the geocode
  !conditions for periodicity in the three directions
  perx=(geocode /= 'F')
  pery=(geocode == 'P')
  perz=(geocode /= 'F')

  call ext_buffers(perx,nbl1,nbr1)
  call ext_buffers(pery,nbl2,nbr2)
  call ext_buffers(perz,nbl3,nbr3)

  IntSur=0.d0
  IntVol=0.d0
  pi = 4.d0*datan(1.d0)
  r=0.d0
  t=0.d0
  oneoeps0=1.d0/epsilon0
  oneosqrteps0=1.d0/dsqrt(epsilon0)

  do i3=1,ndims(3)
   z=hgrids(3)*(i3-1-nbl3)
   v(3)=z
   do i2=1,ndims(2)
    y=hgrids(2)*(i2-1-nbl2)
    v(2)=y
    do i1=1,ndims(1)
     x=hgrids(1)*(i1-1-nbl1)
     v(1)=x
     do iat=1,nat
      dmax = radii(iat) - 2.40d0*delta
      dmin = radii(iat) + 1.60d0*delta
      fact1=2.d0*pi/(-(dmax**2) + dmin**2)
      fact2=(dlog(2.d0))/(2.d0*pi)
      fact3=(dlog(2.d0))/(-(dmax**2) + dmin**2)
      d2=(x-rxyz(1,iat))**2+(y-rxyz(2,iat))**2+(z-rxyz(3,iat))**2
      if (d2.eq.0.d0) d2=1.0d-30
      d=dsqrt(d2)
      if (d.lt.dmax) then
       ep(iat)=0.d0
       do i=1,3
        dep(i,iat)=0.d0
       end do
       ddep(iat)=0.d0
      else if (d.gt.dmin) then
       ep(iat)=1.d0
       do i=1,3
        dep(i,iat)=0.d0
       end do
       ddep(iat)=0.d0
      else
       r=fact1*(-(dmax**2) + d2)
       t=fact2*(r-dsin(r)) 
       ep(iat)=dexp(t)-1.d0
       dtx=fact3*(1.d0-dcos(r))
       do i=1,3
        dep(i,iat)=dexp(t)*dtx*2.d0*(v(i)-rxyz(i,iat))
       end do
       ddep(iat) = dexp(t)*(4.d0*(dtx**2)*d2 + 4.d0*fact1*fact3*dsin(r)*d2 + 6.d0*dtx)
      end if
     end do

     IntVol = IntVol + (1.d0-product(ep))

     eps(i1,i2,i3)=(epsilon0-1.d0)*product(ep)+1.d0
     oneoeps(i1,i2,i3)=1.d0/eps(i1,i2,i3)
     oneosqrteps(i1,i2,i3)=1.d0/dsqrt(eps(i1,i2,i3))

     do i=1,3
      deps(i)=0.d0
      do jat=0,nat-1
       curr=dep(i,jat+1)
       do iat=1,nat-1
        curr=curr*ep(modulo(iat+jat,nat)+1)
       end do
        deps(i) = deps(i) + curr
      end do
      deps(i) = deps(i)*(epsilon0-1.d0)
     end do

     d12=0.d0
     do i=1,3
      dlogeps(i,i1,i2,i3)=deps(i)/eps(i1,i2,i3)
      d12 = d12 + deps(i)**2
     end do

     IntSur = IntSur + dsqrt(d12)

     dd=0.d0
     do jat=1,nat
      curr=ddep(jat)
      do iat=1,nat-1
       curr=curr*ep(modulo(iat+jat-1,nat)+1)
      end do
      dd = dd + curr
     end do

      do i=1,3
       do iat=1,nat-1
        do jat=iat+1,nat
         curr=dep(i,iat)*dep(i,jat)
         do ii=1,nat
          if ((ii.eq.iat).or.(ii.eq.jat)) then
          else
           curr=curr*ep(ii)
          end if
         end do
         curr=curr*2.d0
         dd = dd + curr
        end do
       end do
      end do

     dd=dd*(epsilon0-1.d0)
     corr(i1,i2,i3)=(-0.125d0/pi)*(0.5d0*d12/eps(i1,i2,i3)-dd)

    end do
   end do
  end do

  IntSur=IntSur*hgrids(1)*hgrids(2)*hgrids(3)/(epsilon0-1.d0)
  IntVol=IntVol*hgrids(1)*hgrids(2)*hgrids(3)

  unt=f_get_free_unit(21)
  call f_open_file(unt,file='epsilon.dat')
  i1=1!n03/2
  do i2=1,ndims(2)
     do i3=1,ndims(3)
        write(unt,'(2(1x,I4),2(1x,e14.7))')i2,i3,eps(i1,i2,i3),eps(ndims(1)/2,i2,i3)
     end do
  end do
  call f_close(unt)

  unt=f_get_free_unit(22)
  call f_open_file(unt,file='epsilon_line.dat')
  do i2=1,ndims(2)
   write(unt,'(1x,I8,1(1x,e22.15))')i2,eps(ndims(1)/2,i2,ndims(3)/2)
  end do
  call f_close(unt)

end subroutine epsilon_rigid_cavity_new_multiatoms


subroutine createEffectiveIonicPotential(iproc, nproc, verb, in, atoms, rxyz, shift, &
     Glr, hxh, hyh, hzh, rhopotd, pkernel, pot_ion, rho_ion, elecfield, psoffset)
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
  type(coulomb_operator), intent(inout) :: pkernel
  real(wp), dimension(*), intent(inout) :: pot_ion
  real(wp), dimension(*), intent(inout) :: rho_ion

  character(len = *), parameter :: subname = "createEffectiveIonicPotential"
  logical :: counterions
  real(dp), dimension(:), allocatable :: counter_ions

  ! Compute the main ionic potential.
  call createIonicPotential(atoms%astruct%geocode, iproc, nproc, verb, atoms, rxyz, hxh, hyh, hzh, &
       & elecfield, Glr%d%n1, Glr%d%n2, Glr%d%n3, rhopotd%n3pi, rhopotd%i3s + rhopotd%i3xcsh, &
       & Glr%d%n1i, Glr%d%n2i, Glr%d%n3i, pkernel, pot_ion, rho_ion, psoffset)

  !inquire for the counter_ion potential calculation (for the moment only xyz format)
  inquire(file='posinp_ci.xyz',exist=counterions)
  if (counterions) then
     if (rhopotd%n3pi > 0) then
        counter_ions = f_malloc(Glr%d%n1i*Glr%d%n2i*rhopotd%n3pi,id='counter_ions')
     else
        counter_ions = f_malloc(1,id='counter_ions')
     end if

     call CounterIonPotential(atoms%astruct%geocode,iproc,nproc,in,shift,&
          &   hxh,hyh,hzh,Glr%d,rhopotd%n3pi,rhopotd%i3s + rhopotd%i3xcsh,pkernel,counter_ions)

     !sum that to the ionic potential
     call axpy(Glr%d%n1i*Glr%d%n2i*rhopotd%n3pi,1.0_dp,counter_ions(1),1,&
          &   pot_ion(1),1)

     call f_free(counter_ions)
  end if


END SUBROUTINE createEffectiveIonicPotential


!> Create the ionic potential
subroutine createIonicPotential(geocode,iproc,nproc,verb,at,rxyz,&
     hxh,hyh,hzh,elecfield,n1,n2,n3,n3pi,i3s,n1i,n2i,n3i,pkernel,pot_ion,rho_ion,psoffset)
  use module_base, pi => pi_param
  use m_splines, only: splint
  use module_types
  use yaml_output
  use gaussians, only: initialize_real_space_conversion, finalize_real_space_conversion,mp_exp
!  use module_interfaces, except_this_one => createIonicPotential
  use Poisson_Solver, except_dp => dp, except_gp => gp, except_wp => wp
  use public_enums, only: PSPCODE_PAW
  implicit none
  character(len=1), intent(in) :: geocode !< @copydoc poisson_solver::doc::geocode
  integer, intent(in) :: iproc,nproc,n1,n2,n3,n3pi,i3s,n1i,n2i,n3i
  logical, intent(in) :: verb
  real(gp), intent(in) :: hxh,hyh,hzh,psoffset
  type(atoms_data), intent(in) :: at
  real(gp), dimension(3), intent(in) :: elecfield
  real(gp), dimension(3,at%astruct%nat), intent(in) :: rxyz
  type(coulomb_operator), intent(inout) :: pkernel
  real(wp), dimension(*), intent(inout) :: pot_ion
  real(dp), dimension(*), intent(out) :: rho_ion
  !local variables
  character(len=*), parameter :: subname='createIonicPotential'
  character(len = 3) :: quiet
  logical :: perx,pery,perz,gox,goy,goz,htoobig=.false.,efwrite,check_potion=.false.
  integer :: iat,i1,i2,i3,j1,j2,j3,isx,isy,isz,iex,iey,iez,ierr,ityp !n(c) nspin
  integer :: ind,nbl1,nbr1,nbl2,nbr2,nbl3,nbr3,nloc,iloc,indj3,indj23
  real(kind=8) :: rholeaked,rloc,charge,cutoff,x,y,z,r2,arg,xp,tt,rx,ry,rz
  real(kind=8) :: tt_tot,rholeaked_tot,potxyz
  real(kind=8) :: raux2,r2paw,rlocsq,zsq,yzsq
  real(kind=8) :: raux(1),rr(1)
  real(wp) :: maxdiff
  real(gp) :: ehart
  real(dp), dimension(2) :: charges_mpi
  real(dp), dimension(:), allocatable :: potion_corr
  !real(dp), dimension(:), allocatable :: den_aux

  call timing(iproc,'CrtLocPot     ','ON')

  !initialize the work arrays needed to integrate with isf
  if (at%multipole_preserving) call initialize_real_space_conversion(isf_m=at%mp_isf)

  ! Ionic charge (must be calculated for the PS active processes)
  rholeaked=0.d0
  ! Ionic energy (can be calculated for all the processors)

  !Creates charge density arising from the ionic PSP cores
  call f_zero(n1i*n2i*n3pi,pot_ion(1))

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
        rlocsq=rloc**2
        charge=real(at%nelpsp(ityp),kind=8)/(2.d0*pi*sqrt(2.d0*pi)*rloc**3)
        !cutoff of the range

        cutoff=10.d0*rloc
        if (at%multipole_preserving) then
           !We want to have a good accuracy of the last point rloc*10
           cutoff=cutoff+max(hxh,hyh,hzh)*real(at%mp_isf,kind=gp)
        end if

        isx=floor((rx-cutoff)/hxh)
        isy=floor((ry-cutoff)/hyh)
        isz=floor((rz-cutoff)/hzh)

        iex=ceiling((rx+cutoff)/hxh)
        iey=ceiling((ry+cutoff)/hyh)
        iez=ceiling((rz+cutoff)/hzh)

!       Calculate Ionic Density
!       using HGH parameters.
!       Eq. 1.104, T. Deutsch and L. Genovese, JDN. 12, 2011
        if( .not. any(at%npspcode == PSPCODE_PAW) ) then

           do i3=isz,iez
              z=real(i3,kind=8)*hzh-rz
              call ind_positions(perz,i3,n3,j3,goz) 
              j3=j3+nbl3+1
              if (goz .and. (j3<i3s.or.j3>i3s+n3pi-1)) cycle
              indj3=(j3-i3s)*n1i*n2i
              zsq=z**2
              do i2=isy,iey
                 y=real(i2,kind=8)*hyh-ry
                 call ind_positions(pery,i2,n2,j2,goy)
                 if (goz.and.(.not.goy)) cycle
                 indj23=1+nbl1+(j2+nbl2)*n1i+indj3
                 yzsq=y**2+zsq
                 do i1=isx,iex
                    x=real(i1,kind=8)*hxh-rx
                    call ind_positions(perx,i1,n1,j1,gox)
                    r2=x**2+yzsq
                    arg=r2/rlocsq
                    if (at%multipole_preserving) then
                       !use multipole-preserving function
                       xp=mp_exp(hxh,rx,0.5_gp/(rlocsq),i1,0,.true.)*&
                          mp_exp(hyh,ry,0.5_gp/(rlocsq),i2,0,.true.)*&
                          mp_exp(hzh,rz,0.5_gp/(rlocsq),i3,0,.true.)
                    else
                       xp=exp(-.5d0*arg)
                    end if
                    if (j3 >= i3s .and. j3 <= i3s+n3pi-1 .and. goy .and. gox) then
                       ind=j1+indj23
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
           r2paw=at%pawtab(ityp)%rpaw**2
           do i3=isz,iez
              z=real(i3,kind=8)*hzh-rz
              call ind_positions(perz,i3,n3,j3,goz)
              j3=j3+nbl3+1
              indj3=(j3-i3s)*n1i*n2i
              zsq=z**2
              do i2=isy,iey
                 y=real(i2,kind=8)*hyh-ry
                 call ind_positions(pery,i2,n2,j2,goy)
                 indj23=1+nbl1+(j2+nbl2)*n1i+indj3
                 yzsq=y**2+zsq
                 do i1=isx,iex
                    x=real(i1,kind=8)*hxh-rx
                    call ind_positions(perx,i1,n1,j1,gox)
                    r2=x**2+yzsq
                    !if(r2>r2paw) cycle
                    rr=sqrt(r2)
                    if(1==2) then
                      !This converges very slow                
                      call splint(at%pawtab(ityp)%wvl%rholoc%msz, &
                           & at%pawtab(ityp)%wvl%rholoc%rad, &
                           & at%pawtab(ityp)%wvl%rholoc%d(:,1), &
                           & at%pawtab(ityp)%wvl%rholoc%d(:,2), &
                           & 1,rr,raux,ierr)
                    else
                      !Take the HGH form for rho_L (long range)
                      arg=r2/rlocsq
                      if (at%multipole_preserving) then
                         !use multipole-preserving function
                         xp=mp_exp(hxh,rx,0.5_gp/(rlocsq),i1,0,.true.)*&
                            mp_exp(hyh,ry,0.5_gp/(rlocsq),i2,0,.true.)*&
                            mp_exp(hzh,rz,0.5_gp/(rlocsq),i3,0,.true.)
                      else
                         xp=exp(-.5d0*arg)
                      end if
                      raux=-xp*charge
                    end if
                    !raux=-4.d0**(3.0d0/2.0d0)*exp(-4.d0*pi*r2)

                    if (j3 >= i3s .and. j3 <= i3s+n3pi-1  .and. goy  .and. gox ) then
                       ind=j1+indj23
                       pot_ion(ind)=pot_ion(ind)+raux(1)
                    else if (.not. goz ) then
                       rholeaked=rholeaked-raux(1)
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
     indj3=(j3-1)*n1i*n2i
     do i2= -nbl2,2*n2+1+nbr2
        indj23=1+nbl1+(i2+nbl2)*n1i+indj3
        do i1= -nbl1,2*n1+1+nbr1
           ind=i1+indj23
           tt=tt+pot_ion(ind)
        enddo
     enddo
  enddo

  tt=tt*hxh*hyh*hzh
  rholeaked=rholeaked*hxh*hyh*hzh

  !print *,'test case input_rho_ion',iproc,i3start,i3end,n3pi,2*n3+16,tt
  !if rho_ion is needed for the SCF cycle copy in the array
  if (pkernel%method /= 'VAC') then
     call f_memcpy(n=n1i*n2i*n3pi,src=pot_ion(1),dest=rho_ion(1))
  end if

  if (pkernel%mpi_env%nproc > 1) then
     charges_mpi(1)=tt
     charges_mpi(2)=rholeaked

     call mpiallred(charges_mpi,MPI_SUM,comm=pkernel%mpi_env%mpi_comm)

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

     !in the case of vacuum the pot_ion treatment is as usual
     !otherwise the pot_ion array is set to zero and can be filled with external potentials
     !like the gaussian part of the PSP ad/or external electric fields
     if (pkernel%method /= 'VAC') then
        call f_zero(n1i*n2i*n3pi,pot_ion(1))
     else
        call H_potential('D',pkernel,pot_ion,pot_ion,ehart,-psoffset,.false.,quiet=quiet)
     end if

     call timing(iproc,'CrtLocPot     ','ON')
     
     if (check_potion) then
        !if (iproc == 0) write(*,'(1x,a)',advance='no') 'Check the ionic potential...'
          
        potion_corr = f_malloc0(n1i*n2i*n3pi,id='potion_corr')

        !call to_zero(n1i*n2i*n3pi,potion_corr)

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

        if (pkernel%mpi_env%nproc > 1) then
           call mpiallred(maxdiff,1,MPI_MAX,comm=pkernel%mpi_env%mpi_comm)
        end if

        if (iproc == 0) call yaml_map('Check the ionic potential',maxdiff,fmt='(1pe24.17)')
        !if (iproc == 0) write(*,'(1x,a,1pe24.17)')'...done. MaxDiff=',maxdiff

        stop

        call f_free(potion_corr)

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
        rlocsq=rloc**2
        cutoff=10.d0*rloc
        if (at%multipole_preserving) then
           !We want to have a good accuracy of the last point rloc*10
           cutoff=cutoff+max(hxh,hyh,hzh)*real(16,kind=gp)
        end if

        isx=floor((rx-cutoff)/hxh)
        isy=floor((ry-cutoff)/hyh)
        isz=floor((rz-cutoff)/hzh)

        iex=ceiling((rx+cutoff)/hxh)
        iey=ceiling((ry+cutoff)/hyh)
        iez=ceiling((rz+cutoff)/hzh)
        
        if( at%npspcode(1) /= PSPCODE_PAW) then

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
                 indj3=(j3-i3s)*n1i*n2i
                 zsq=z**2
                 if (goz .and. j3 >= i3s .and. j3 <=  i3s+n3pi-1) then
                    do i2=isy,iey
                       y=real(i2,kind=8)*hyh-ry
                       call ind_positions(pery,i2,n2,j2,goy)
                       indj23=1+nbl1+(j2+nbl2)*n1i+indj3
                       yzsq=y**2+zsq
                       if (goy) then
                          do i1=isx,iex
                             x=real(i1,kind=8)*hxh-rx
                             call ind_positions(perx,i1,n1,j1,gox)
                             if (gox) then
                                r2=x**2+yzsq
                                arg=r2/rlocsq

                                if (at%multipole_preserving) then
                                   !use multipole-preserving function
                                   xp=mp_exp(hxh,rx,0.5_gp/(rlocsq),i1,0,.true.)*&
                                      mp_exp(hyh,ry,0.5_gp/(rlocsq),i2,0,.true.)*&
                                      mp_exp(hzh,rz,0.5_gp/(rlocsq),i3,0,.true.)
                                 else
                                    xp=exp(-.5d0*arg)
                                 end if

                                tt=at%psppar(0,nloc,ityp)
                                do iloc=nloc-1,1,-1
                                   tt=arg*tt+at%psppar(0,iloc,ityp)
                                enddo
                                ind=j1+indj23
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
              indj3=(j3-i3s)*n1i*n2i
              zsq=z**2
              if (goz .and. j3 >= i3s .and. j3 <=  i3s+n3pi-1) then
                 do i2=isy,iey
                    y=real(i2,kind=8)*hyh-ry
                    call ind_positions(pery,i2,n2,j2,goy)
                    indj23=1+nbl1+(j2+nbl2)*n1i+indj3
                    yzsq=y**2+zsq
                    if (goy) then
                       do i1=isx,iex
                          x=real(i1,kind=8)*hxh-rx
                          call ind_positions(perx,i1,n1,j1,gox)
                          if (gox) then
                             r2=x**2+yzsq
                             rr(1)=sqrt(r2)
                             !1) V_L^HGH
                             if(rr(1)>0.01d0) then
                               arg=rr(1)/(sqrt(2.0)*rloc)
                               call derf_ab(tt,arg)
                               raux2=-charge/rr(1)*tt  
                             else
                               !In this case we deduce the values
                               !from a quadratic interpolation (due to 1/rr factor)
                               call interpol_vloc(rr(1),rloc,charge,raux2)
                             end if
                             !2) V^PAW from splines
                             call splint(at%pawtab(ityp)%wvl%rholoc%msz, &
                                  & at%pawtab(ityp)%wvl%rholoc%rad, &
                                  & at%pawtab(ityp)%wvl%rholoc%d(:,3), &
                                  & at%pawtab(ityp)%wvl%rholoc%d(:,4), &
                                  & 1,rr,raux,ierr)
                             
                             ind=j1+indj23
                             pot_ion(ind)=pot_ion(ind)+raux(1)-raux2
                          end if
                       enddo
                    end if
                 enddo
              end if
           end do
        end if ! at%npspcode(iat) /= PSPCODE_PAW
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

  !use rhopotential to calculate the potential from a constant electric field along y direction
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

  if (at%multipole_preserving) call finalize_real_space_conversion()

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
  use module_base, pi => pi_param
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
  real(wp) :: charge
  real(gp) :: r,sq2rl,rx,ry,rz,derf_val
  
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
  use module_base, pi => pi_param
  use module_types
  !use module_interfaces, except_this_one => CounterIonPotential
  use Poisson_Solver, except_dp => dp, except_gp => gp, except_wp => wp
  use module_input_dicts
  use public_keys, only: IG_OCCUPATION
  use dictionaries
  use yaml_output
  use module_atoms
  use gaussians, only: initialize_real_space_conversion, finalize_real_space_conversion,mp_exp
  implicit none
  character(len=1), intent(in) :: geocode !< @copydoc poisson_solver::doc::geocode
  integer, intent(in) :: iproc,nproc,n3pi,i3s
  real(gp), intent(in) :: hxh,hyh,hzh
  real(gp), dimension(3), intent(in) :: shift
  type(input_variables), intent(in) :: in
  type(grid_dimensions), intent(in) :: grid
  type(coulomb_operator), intent(inout) :: pkernel
  real(wp), dimension(*), intent(inout) :: pot_ion
  !local variables
  character(len=*), parameter :: subname='CounterIonPotential'
  logical :: htoobig=.false.,check_potion=.false.
  logical :: perx,pery,perz,gox,goy,goz
  integer :: iat,i1,i2,i3,j1,j2,j3,isx,isy,isz,iex,iey,iez,ityp,nspin
  integer :: ind,nbl1,nbr1,nbl2,nbr2,nbl3,nbr3
  real(kind=8) :: rholeaked,rloc,charge,cutoff,x,y,z,r2,arg,xp,tt,rx,ry,rz
  real(kind=8) :: tt_tot,rholeaked_tot,potxyz
  real(wp) :: maxdiff
  real(gp) :: ehart
  type(atoms_data) :: at
  type(dictionary), pointer :: dict
  real(dp), dimension(2) :: charges_mpi
  real(dp), dimension(:), allocatable :: potion_corr
!  real(gp), dimension(:,:), allocatable :: radii_cf

  call timing(iproc,'CrtLocPot     ','ON')
  
  !initialize the work arrays needed to integrate with isf
  if (at%multipole_preserving) call initialize_real_space_conversion(isf_m=at%mp_isf)

  if (iproc.eq.0) then
     write(*,'(1x,a)')&
          '--------------------------------------------------- Counter Ionic Potential Creation'
  end if

  at = atoms_data_null()
  !read the positions of the counter ions from file
  call dict_init(dict)
  call astruct_file_merge_to_dict(dict, "posinp", 'posinp_ci')
  call astruct_set_from_dict(dict // "posinp", at%astruct)

  call atoms_file_merge_to_dict(dict)
  do ityp = 1, at%astruct%ntypes, 1
     call psp_dict_fill_all(dict, at%astruct%atomnames(ityp), in%ixc, in%projrad, in%crmult, in%frmult)
  end do
  call psp_dict_analyse(dict, at)
  ! Read associated pseudo files.
  call atomic_data_set_from_dict(dict,IG_OCCUPATION, at, in%nspin)
  call dict_free(dict)

  !read the specifications of the counter ions from pseudopotentials
!  radii_cf = f_malloc((/ at%astruct%ntypes, 3 /),id='radii_cf')
!  radii_cf = at%radii_cf
  if (iproc == 0) call print_atomic_variables(at, max(in%hx,in%hy,in%hz), in%ixc, in%dispersion)

  ! Ionic charge (must be calculated for the PS active processes)
  rholeaked=0.d0
  ! Ionic energy (can be calculated for all the processors)

  !Creates charge density arising from the ionic PSP cores
  call f_zero(grid%n1i*grid%n2i*n3pi,pot_ion(1))


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
        if (at%multipole_preserving) then
           !We want to have a good accuracy of the last point rloc*10
           cutoff=cutoff+max(hxh,hyh,hzh)*real(16,kind=gp)
        end if

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

                 if (at%multipole_preserving) then
                    !multipole-preserving approach
                    xp=mp_exp(hxh,rx,0.5_gp/(rloc**2),i1,0,.true.)*&
                       mp_exp(hyh,ry,0.5_gp/(rloc**2),i2,0,.true.)*&
                       mp_exp(hzh,rz,0.5_gp/(rloc**2),i3,0,.true.)
                 else
                    xp=exp(-.5d0*arg)
                 end if

                 if (j3 >= i3s .and. j3 <= i3s+n3pi-1  .and. goy  .and. gox ) then
                    ind=j1+1+nbl1+(j2+nbl2)*grid%n1i+(j3-i3s)*grid%n1i*grid%n2i
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

     call mpiallred(charges_mpi,MPI_SUM,comm=pkernel%mpi_env%mpi_comm)

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
          
        potion_corr = f_malloc0(grid%n1i*grid%n2i*n3pi,id='potion_corr')

        !call to_zero(grid%n1i*grid%n2i*n3pi,potion_corr)

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

        if (pkernel%mpi_env%nproc > 1) then
           call mpiallred(maxdiff,1,MPI_MAX,comm=pkernel%mpi_env%mpi_comm)
        end if

        if (iproc == 0) call yaml_map('Check the ionic potential',maxdiff,fmt='(1pe24.17)')
        !if (iproc == 0) write(*,'(1x,a,1pe24.17)')'...done. MaxDiff=',maxdiff

        stop

        call f_free(potion_corr)

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
  call deallocate_atoms_data(at) 

!  call f_free(radii_cf)

  call f_free_ptr(at%astruct%rxyz)


  if (at%multipole_preserving) call finalize_real_space_conversion()

  call timing(iproc,'CrtLocPot     ','OF')

END SUBROUTINE CounterIonPotential
