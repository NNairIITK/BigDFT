!{\src2tex{textfont=tt}}
!!****f* ABINIT/isotemp
!! NAME
!! isotemp
!!
!! FUNCTION
!! performs one half step on isotemp parameters according to Martyna et al.
!!
!! COPYRIGHT
!! Copyright (C) 1998-2009 ABINIT group (DCA, XG, GMR, JCC, JYR, SE)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  amass(natom)=mass of each atom, in unit of electronic mass (=amu*1822...)
!!  dtion=
!!  dtset <type(dataset_type)>=all input variables for this dataset
!!   | natom=number of atoms in unit cell
!!  isotemp_data
!!  ktemp
!!  vel
!!
!! OUTPUT
!!  Only updates variables
!!
!! SIDE EFFECTS
!!  isotemp_data: updates the thermostat parameters
!!  vel=update the velocities
!!
!! NOTES
!! This program is decribed in the following paper 
!! Explicit integrators for extended systems dynamics
!! Glenn J Martyna et al.
!! Mol. Phys., 1996, Vol. 87, pp. 1117-1157
!!
!! SOURCE

subroutine isotemp(amass,dtion,ekin,iatfix,ktemp,mttk_vars,natom,nnos,qmass,vel)

 use defs_basis
 use m_abimover

 implicit none

!Arguments ------------------------------------
!scalars
 integer :: natom, nnos
 real(dp),intent(in) :: dtion,ktemp
 real(dp),intent(out) :: ekin
 type(mttk_type) :: mttk_vars
!arrays
 integer, intent(in) :: iatfix(3, natom)
 real(dp),intent(in) :: amass(natom), qmass(nnos)
 real(dp),intent(inout) :: vel(3,natom)

!Local variables ------------------------------
!scalars
 integer :: iatom,idir,inos
 real(dp) :: alocal,gnkt,nfree,scale
!arrays
 real(dp),allocatable :: glogs(:),vlogs(:),xlogs(:)

!***************************************************************************
!Beginning of executable session
!***************************************************************************
 allocate(glogs(nnos),vlogs(nnos),xlogs(nnos))
 glogs(:)=mttk_vars%glogs(:)
 vlogs(:)=mttk_vars%vlogs(:)
 xlogs(:)=mttk_vars%xlogs(:)
 scale=one
!Compute the ionic kinetic energy 
 nfree=zero
 ekin=zero
 do iatom=1,natom
  do idir=1,3
!  Warning : the fixing of atomis is implemented in reduced
!  coordinates, so that this expression is wrong
   if (iatfix(idir,iatom) == 0) then
    ekin=ekin+0.5d0*amass(iatom)*vel(idir,iatom)**2
!   Counts the degrees of freedom
    nfree=nfree+one
   end if
  end do
 end do
 gnkt=nfree*ktemp
!Update the forces
 glogs(1)=(two*ekin-gnkt)/qmass(1)
 vlogs(nnos)=vlogs(nnos)+glogs(nnos)*dtion/four
 do inos=1,nnos-1
  alocal=exp(-dtion/eight*vlogs(nnos+1-inos))
  vlogs(nnos-inos)=vlogs(nnos-inos)*alocal*alocal+&
&  dtion/four*glogs(nnos-inos)*alocal
 end do
!Update the particle velocities
 alocal=exp(-dtion/two*vlogs(1))
 scale=scale*alocal
!Update the forces
 glogs(1)=(scale*scale*two*ekin-gnkt)/qmass(1)
!Update the thermostat positions
 do inos=1,nnos
  xlogs(inos)=xlogs(inos)+vlogs(inos)*dtion/two
 end do
!Update the thermostat velocities
 do inos=1,nnos-1
  alocal=exp(-dtion/eight*vlogs(inos+1))
  vlogs(inos)=vlogs(inos)*alocal*alocal+dtion/four*glogs(inos)*alocal
  glogs(inos+1)=(qmass(inos)*vlogs(inos)*vlogs(inos)-ktemp)/qmass(inos+1)
 end do
 vlogs(nnos)=vlogs(nnos)+glogs(nnos)*dtion/four
 vel(:,:)=vel(:,:)*scale
!Compute the ionic kinetic energy 
 ekin=zero
 do iatom=1,natom
  do idir=1,3
!  Warning : the fixing of atomis is implemented in reduced
!  coordinates, so that this expression is wrong
   if (iatfix(idir,iatom) == 0) then
    ekin=ekin+half*amass(iatom)*vel(idir,iatom)**2
   end if
  end do
 end do
!Compute the thermostat kinetic energy and add it to the ionic one
 ekin=ekin+half*qmass(1)*vlogs(1)**2+xlogs(1)*nfree*ktemp
 do inos=2,nnos
  ekin=ekin+half*qmass(inos)*vlogs(inos)**2+xlogs(inos)*ktemp
 end do
 mttk_vars%glogs(:)=glogs(:)
 mttk_vars%vlogs(:)=vlogs(:)
 mttk_vars%xlogs(:)=xlogs(:)
 deallocate(glogs,vlogs,xlogs)

end subroutine isotemp

!!***
!!****f* ABINIT/isopress
!! NAME
!! isopress
!!
!! FUNCTION
!! performs one half step on isopress parameters according to Martyna et al.
!!
!! COPYRIGHT
!! Copyright (C) 1998-2009 ABINIT group (DCA, XG, GMR, JCC, JYR, SE)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  amass(natom)=mass of each atom, in unit of electronic mass (=amu*1822...)
!!  dtion= ionic time step
!!  dtset <type(dataset_type)>=all input variables for this dataset
!!   | natom=number of atoms in unit cell
!!  isotemp_data
!!  ktemp
!!  press= current pressure of the system
!!  prtarget= target pressure
!!  ucvol= unit cell volume
!!  vel= current velocity
!! OUTPUT
!!  Only updates variables
!! SIDE EFFECTS
!!  isotemp_data: updates the thermostat parameters (saved variables: bouh !)
!!  vel=update the velocities
!! NOTES
!! This program is decribed in the following paper 
!! Explicit integrators for extended systems dynamics
!! Glenn J Martyna et al.
!! Mol. Phys., 1996, Vol. 87, pp. 1117-1157
!!
!! SOURCE

subroutine isopress(amass,dtion,ekin,iatfix,ktemp,natom,nnos,qmass,strten,strtarget, &
&                   ucvol,mttk_vars,vel,vlogv,vmass)

 use defs_basis
 use m_abimover

 implicit none

!Arguments ------------------------------------
!scalars
 integer, intent(in) :: natom, nnos
 real(dp),intent(in) :: dtion,ktemp,ucvol,vmass
 real(dp),intent(inout) :: vlogv
 real(dp),intent(out) :: ekin
 type(mttk_type) :: mttk_vars
!arrays
 integer, intent(in) :: iatfix(3, natom)
 real(dp),intent(in) :: amass(natom), qmass(nnos),strtarget(6),strten(6)
 real(dp),intent(inout) :: vel(3,natom)

!Local variables ------------------------------
!scalars
 integer :: iatom,idir,inos
 real(dp) :: alocal,glogv,gn1kt,gnkt,nfree,odnf,press,prtarget,scale
!arrays
 real(dp),allocatable :: glogs(:),vlogs(:),xlogs(:)

!***************************************************************************
!Beginning of executable session
!***************************************************************************
 allocate(glogs(nnos),vlogs(nnos),xlogs(nnos))
 glogs(:)=mttk_vars%glogs(:)
 vlogs(:)=mttk_vars%vlogs(:)
 xlogs(:)=mttk_vars%xlogs(:)
 glogv   =mttk_vars%glogv
 scale=one
!Compute the ionic kinetic energy 
 nfree=zero
 ekin=zero
 do iatom=1,natom
  do idir=1,3
!  Warning : the fixing of atomis is implemented in reduced
!  coordinates, so that this expression is wrong
   if (iatfix(idir,iatom) == 0) then
    ekin=ekin+0.5d0*amass(iatom)*vel(idir,iatom)**2
!   Counts the degrees of freedom
    nfree=nfree+one
   end if
  end do
 end do
 prtarget=-(strtarget(1)+strtarget(2)+strtarget(3))/three
 press=-(strten(1)+strten(2)+strten(3))/three
 gnkt=nfree*ktemp
 gn1kt=(nfree+one)*ktemp
 odnf=one+three/nfree
!Update the forces
 glogs(1)=(two*ekin+vmass*vlogv*vlogv-gn1kt)/qmass(1)
 glogv=(odnf*two*ekin+three*(press-prtarget)*ucvol)/vmass
!Update thermostat velocity
 vlogs(nnos)=vlogs(nnos)+glogs(nnos)*dtion/four
 do inos=1,nnos-1
  alocal=exp(-dtion/eight*vlogs(nnos+1-inos))
  vlogs(nnos-inos)=vlogs(nnos-inos)*alocal*alocal+&
&  dtion/four*glogs(nnos-inos)*alocal
 end do
!Update dLog(V)/dt
 alocal=exp(-dtion/eight*vlogs(1))
 vlogv=vlogv*alocal**2+dtion/four*glogv*alocal
!Update the particle velocities
 alocal=exp(-dtion/two*(vlogs(1)+odnf*vlogv))
 scale=scale*alocal
 ekin=ekin*alocal**2
 glogv=(odnf*two*ekin+three*(press-prtarget)*ucvol)/vmass
!Update the thermostat positions
 do inos=1,nnos
  xlogs(inos)=xlogs(inos)+vlogs(inos)*dtion/two
 end do
!Update dLog(V)/dt
 alocal=exp(-dtion/eight*vlogs(1))
 vlogv=vlogv*alocal**2+dtion/four*glogv*alocal 
!Update the forces
 glogs(1)=(two*ekin+vmass*vlogv*vlogv-gn1kt)/qmass(1)
!Update the thermostat velocities
 do inos=1,nnos-1
  alocal=exp(-dtion/eight*vlogs(inos+1))
  vlogs(inos)=vlogs(inos)*alocal*alocal+dtion/four*glogs(inos)*alocal
  glogs(inos+1)=(qmass(inos)*vlogs(inos)*vlogs(inos)-ktemp)/qmass(inos+1)
 end do
 vlogs(nnos)=vlogs(nnos)+glogs(nnos)*dtion/four
 vel(:,:)=vel(:,:)*scale
!Compute the ionic kinetic energy 
 ekin=zero
 do iatom=1,natom
  do idir=1,3
!  Warning : the fixing of atomis is implemented in reduced
!  coordinates, so that this expression is wrong
   if (iatfix(idir,iatom) == 0) then
    ekin=ekin+half*amass(iatom)*vel(idir,iatom)**2
   end if
  end do
 end do
!Compute the thermostat kinetic energy and add it to the ionic one
!First thermostat
 ekin=ekin+half*qmass(1)*vlogs(1)**2+xlogs(1)*(nfree+one)*ktemp
!Other thermostats
 do inos=2,nnos
  ekin=ekin+half*qmass(inos)*vlogs(inos)**2+xlogs(inos)*ktemp
 end do
!Barostat
 ekin=ekin+half*vmass*vlogv**2+prtarget*ucvol
 
 deallocate(glogs,vlogs,xlogs)

end subroutine isopress

!!***
!!****f* ABINIT/isostress
!! NAME
!! isostress
!!
!! FUNCTION
!! performs one half step on isostress parameters according to Martyna et al.
!!
!! COPYRIGHT
!! Copyright (C) 1998-2009 ABINIT group (DCA, XG, GMR, JCC, JYR, SE)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  amass(natom)=mass of each atom, in unit of electronic mass (=amu*1822...)
!!  dtion= ionic time step
!!  dtset <type(dataset_type)>=all input variables for this dataset
!!   | natom=number of atoms in unit cell
!!  isotemp_data
!!  ktemp
!!  press= current pressure of the system
!!  prtarget= target pressure
!!  ucvol= unit cell volume
!!  vel= current velocity
!!
!! OUTPUT
!!  Only updates variables
!!
!! SIDE EFFECTS
!!  isotemp_data: updates the thermostat parameters (saved variables: bouh !)
!!  vel=update the velocities
!!
!! NOTES
!! This program is decribed in the following paper 
!! Explicit integrators for extended systems dynamics
!! Glenn J Martyna et al.
!! Mol. Phys., 1996, Vol. 87, pp. 1117-1157
!!
!! SOURCE

subroutine isostress(amass,bmass,dtion,ekin,iatfix,ktemp,natom,nnos,qmass,strten, &
&                    strtarget,ucvol,vel,mttk_vars)

 use defs_basis
 use m_abimover

 implicit none

!Arguments ------------------------------------
!scalars
 integer, intent(in) :: natom,nnos
 real(dp),intent(in) :: dtion,ktemp,ucvol,bmass
 real(dp),intent(out) :: ekin
 type(mttk_type) :: mttk_vars
!arrays
 integer, intent(in) :: iatfix(3, natom)
 real(dp),intent(in) :: amass(natom),qmass(nnos),strtarget(6),strten(6)
 real(dp),intent(inout) :: vel(3,natom)

!Local variables ------------------------------
!scalars
 integer,parameter :: lwork=8
 integer :: iatom,idir,info,inos,jdir
 real(dp) :: akinb,alocal,gn1kt,gnd2kt,nfree,odnf,trvg
!arrays
 real(dp) :: akin(3,3),expdiag(3),gboxg(3,3),identity(3,3),press(3,3)
 real(dp) :: prtarget(3,3),tvtemp(3,3),uv(3),vboxg(3,3),veig(3),vtemp(3,3)
 real(dp) :: work(lwork)
 real(dp),allocatable :: glogs(:),vlogs(:),xlogs(:)

!***************************************************************************
!Beginning of executable session
!***************************************************************************

 allocate(glogs(nnos),vlogs(nnos),xlogs(nnos))
 glogs(:)=mttk_vars%glogs(:)
 vlogs(:)=mttk_vars%vlogs(:)
 xlogs(:)=mttk_vars%xlogs(:)
 vboxg(:,:)=mttk_vars%vboxg(:,:)
 identity(:,:)=zero
 do idir=1,3
  identity(idir,idir)=one
 end do
!Compute the ionic kinetic energy 
 nfree=zero
 ekin=zero
 do iatom=1,natom
  do idir=1,3
!  Warning : the fixing of atomis is implemented in reduced
!  coordinates, so that this exprtargetion is wrong
   if (iatfix(idir,iatom) == 0) then
    ekin=ekin+0.5d0*amass(iatom)*vel(idir,iatom)**2
!   Counts the degrees of freedom
    nfree=nfree+one
   end if
  end do
 end do
 gn1kt=(nfree+one)*ktemp
 gnd2kt=(nfree+9)*ktemp
 odnf=one+three/nfree
 akin(:,:)=zero
 do iatom=1,natom
  do idir=1,3
   do jdir=1,3
!   Warning : the fixing of atomis is implemented in reduced
!   coordinates, so that this expression is wrong
    akin(idir,jdir)=akin(idir,jdir)+0.5d0*amass(iatom)*vel(idir,iatom)*vel(jdir,iatom)
   end do
  end do
 end do
 akinb=zero
 do idir=1,3
  do jdir=1,3
   akinb=akinb+0.5d0*bmass*vboxg(idir,jdir)**2
  end do
 end do
!Compute the pressure: from Voigt to tensor notation+kinetic energy
 do idir=1,3
  press(idir,idir)=-strten(idir)
  prtarget(idir,idir)=-strtarget(idir)
 end do
 press(3,2)=-strten(4); press(1,3)=-strten(5); press(2,1)=-strten(6)
 prtarget(3,2)=-strtarget(4); prtarget(1,3)=-strtarget(5); prtarget(2,1)=-strtarget(6)
 press(2,3)=press(3,2); press(3,1)=press(1,3); press(1,2)=press(2,1)
 prtarget(2,3)=prtarget(3,2); prtarget(3,1)=prtarget(1,3); prtarget(1,2)=prtarget(2,1)
!Update the forces
 glogs(1)=(two*ekin+two*akinb-gnd2kt)/qmass(1)
 gboxg(:,:)=(two*ekin/nfree*identity(:,:)+two*akin(:,:)+(press(:,:)-prtarget(:,:))*ucvol)/bmass
!Update thermostat velocity
 vlogs(nnos)=vlogs(nnos)+glogs(nnos)*dtion/four
 do inos=1,nnos-1
  alocal=exp(-dtion/eight*vlogs(nnos+1-inos))
  vlogs(nnos-inos)=vlogs(nnos-inos)*alocal*alocal+&
&  dtion/four*glogs(nnos-inos)*alocal
 end do
!Update box velocity
 alocal=exp(-dtion/eight*vlogs(1))
 vboxg(:,:)=vboxg(:,:)*alocal**2+dtion/four*gboxg(:,:)*alocal
!Update the thermostat positions
 do inos=1,nnos
  xlogs(inos)=xlogs(inos)+vlogs(inos)*dtion/two
 end do 
!Update the particle velocities 
 trvg=(vboxg(1,1)+vboxg(2,2)+vboxg(3,3))/nfree
 vtemp(:,:)=vboxg(:,:)+(trvg+vlogs(1))*identity(:,:)
 call dsyev('V','U',3,vtemp,3,veig,work,lwork,info)
!On exit, we have vtemp=U such that tU vtemp U = veig
 tvtemp(:,:)=transpose(vtemp)
 expdiag(1)=exp(-veig(1)*dtion/two)
 expdiag(2)=exp(-veig(2)*dtion/two)
 expdiag(3)=exp(-veig(3)*dtion/two)
 write(6,*)' isostress : expdiag(:)=',expdiag(:)  ! Do not remove this line : seems to be needed for g95 compilo
 do iatom=1,natom
  uv(:)=matmul(tvtemp,vel(:,iatom))
  uv(:)=uv(:)*expdiag(:)
  vel(:,iatom)=matmul(vtemp,uv)
 end do
!Compute the ionic kinetic energy 
 nfree=zero
 ekin=zero
 do iatom=1,natom
  do idir=1,3
!  Warning : the fixing of atomis is implemented in reduced
!  coordinates, so that this expression is wrong
   if (iatfix(idir,iatom) == 0) then
    ekin=ekin+0.5d0*amass(iatom)*vel(idir,iatom)**2
    write(6,*)'kin',iatom,ekin,vel(idir,iatom)
!   Counts the degrees of freedom
    nfree=nfree+one
   end if
  end do
 end do
 gn1kt=(nfree+one)*ktemp
 gnd2kt=(nfree+9)*ktemp
 odnf=one+three/nfree
 akin(:,:)=zero
 do iatom=1,natom
  do idir=1,3
   do jdir=1,3
!   Warning : the fixing of atomis is implemented in reduced
!   coordinates, so that this expression is wrong
    akin(idir,jdir)=akin(idir,jdir)+0.5d0*amass(iatom)*vel(idir,iatom)*vel(jdir,iatom)
   end do
  end do
 end do
 gboxg(:,:)=(two*ekin/nfree*identity(:,:)+two*akin(:,:)+(press(:,:)-prtarget(:,:))*ucvol)/bmass
!Update box velocity
 alocal=exp(-dtion/eight*vlogs(1))
 vboxg(:,:)=vboxg(:,:)*alocal**2+dtion/four*gboxg(:,:)*alocal
!Compute the box kinetic energy
 akinb=zero
 do idir=1,3
  do jdir=1,3
   akinb=akinb+0.5d0*bmass*vboxg(idir,jdir)**2
  end do
 end do
 glogs(1)=(two*ekin+two*akinb-gnd2kt)/qmass(1)
!Update the thermostat velocities
 do inos=1,nnos-1
  alocal=exp(-dtion/eight*vlogs(inos+1))
  vlogs(inos)=vlogs(inos)*alocal*alocal+dtion/four*glogs(inos)*alocal
  glogs(inos+1)=(qmass(inos)*vlogs(inos)*vlogs(inos)-ktemp)/qmass(inos+1)
 end do
 vlogs(nnos)=vlogs(nnos)+glogs(nnos)*dtion/four
!Compute the ionic kinetic energy 
 ekin=zero
 do iatom=1,natom
  do idir=1,3
!  Warning : the fixing of atomis is implemented in reduced
!  coordinates, so that this expression is wrong
   if (iatfix(idir,iatom) == 0) then
    ekin=ekin+half*amass(iatom)*vel(idir,iatom)**2
   end if
  end do
 end do
!Compute the thermostat kinetic energy and add it to the ionic one
!First thermostat
 ekin=ekin+half*qmass(1)*vlogs(1)**2+xlogs(1)*(nfree+nine)*ktemp
!Other thermostats
 do inos=2,nnos
  ekin=ekin+half*qmass(inos)*vlogs(inos)**2+xlogs(inos)*ktemp
 end do
!Barostat kinetic energy
 akinb=zero
 do idir=1,3
  do jdir=1,3
   akinb=akinb+0.5d0*bmass*vboxg(idir,jdir)**2
  end do
 end do
!ekin is the invariant minus the potential energy
 ekin=ekin+akinb+prtarget(1,1)*ucvol
 
 mttk_vars%vboxg(:,:)=vboxg(:,:)
 mttk_vars%glogs(:)=glogs(:)
 mttk_vars%vlogs(:)=vlogs(:)
 mttk_vars%xlogs(:)=xlogs(:)
 deallocate(glogs,vlogs,xlogs)

end subroutine isostress
!!***
