subroutine md_langevin(amass, dtion, fcart, fcart_mold, friction, itime, ktemp, &
     & mditemp, mdwall, natom, rprimd, vel, xcart, xcart_next, xred_next)
  
  use defs_basis

  implicit none

  integer, intent(in) :: natom, itime
  real(dp),intent(in) :: amass(natom), fcart(3, natom)
  real(dp),intent(inout) :: vel(3,natom)
  real(dp),intent(in) :: ktemp, mditemp, friction, dtion, mdwall
  real(dp),intent(in) :: rprimd(3,3)
  real(dp),intent(in) :: xcart(3, natom)
  real(dp),intent(out) :: xred_next(3,natom), xcart_next(3,natom)
  real(dp), intent(inout) :: fcart_mold(3, natom)
  
  interface
     function uniformrandom(seed) 
       implicit none
       integer :: seed
       double precision :: uniformrandom
     end function uniformrandom
  end interface


  real(dp),parameter :: v2tol=tol8
  integer :: iatom, idim, idum=-5
  real(dp) :: delxi, ran_num1, ran_num2, sig_gauss
  real(dp), allocatable :: ran_force(:,:), lang_force(:,:)
 
  if(itime==0)then
     call md_nose_init(amass, natom, mditemp, vel)
  end if

  !  This section is devoted to the optional atom permutation (JYR 001114)
  !  Two input variables are needed
  !  dtset%delayperm : is the interval (in time steps) at which
  !  atoms are tentatively permuted
  !  default value could be 0
  !  dtset%signperm  : is the type of bias for the permutation
  !  +1  to favor alternation of species
  !  -1  to favor segregation

!!$!  Force no permutation at initial step
!!$   if (itime/=0 .and. dtset%delayperm/=0 ) then
!!$    if (mod(itime,dtset%delayperm)==0) then
!!$!    Try commutation of atoms.
!!$     write(message, '(a)')' Attempt of commutation '
!!$     call wrtout(ab_out,message,'COLL')
!!$     call wrtout(std_out,message,'COLL')
!!$!    Compute a 'permutation potential'
!!$     do iatom=1,natom
!!$      pot_perm(iatom)=0.0_dp
!!$      do iatom1=1,natom
!!$       if (iatom1.ne.iatom) then
!!$        distx=xcart(1,iatom)-xcart(1,iatom1)
!!$        distx=distx-acell(1)*nint(distx/acell(1))
!!$        disty=xcart(2,iatom)-xcart(2,iatom1)
!!$        disty=disty-acell(2)*nint(disty/acell(2))
!!$        distz=xcart(3,iatom)-xcart(3,iatom1)
!!$        distz=distz-acell(3)*nint(distz/acell(3))
!!$!       Here we count each atom below 2 angstr as 1, could be customized
!!$        dist=sqrt(distx*distx+disty*disty+distz*distz)/3.7807
!!$        if (typat(iatom).ne.typat(iatom1)) then
!!$         mcfac=-1
!!$        else
!!$         mcfac=1
!!$        end if
!!$        if (dist<1.0_dp)  dist=1.0_dp
!!$        pot_perm(iatom)=pot_perm(iatom)+mcfac*(dtset%signperm)*1.0_dp&
!!$&        /exp(log(dist)*6.0_dp)
!!$       end if
!!$      end do
!!$     end do
!!$     write(message, '(a,10f12.5)' )' Perm_pot ',&
!!$&     (pot_perm(iatom1),iatom1=1,natom)
!!$     call wrtout(ab_out,message,'COLL')
!!$     call wrtout(std_out,message,'COLL')
!!$
!!$!    Find the two atoms, of different types, with the highest perm_pot
!!$     max_perm(:)=-1.0d9
!!$     do iatom=1,natom
!!$      if (pot_perm(iatom) > max_perm(typat(iatom))) then
!!$       max_perm(typat(iatom))=pot_perm(iatom)
!!$       imax_perm(typat(iatom))=iatom
!!$      end if
!!$     end do
!!$!    DEBUG
!!$!    write(message, '(a,10f12.5)' )' max_Perm ',&
!!$!    &      (max_perm(itypat),itypat=1,ntypat)
!!$!    call wrtout(ab_out,message,'COLL')
!!$!    call wrtout(std_out,message,'COLL')
!!$!    write(message, '(a,10i12)' )' imax_Perm ',&
!!$!    &      (imax_perm(itypat),itypat=1,ntypat)
!!$!    call wrtout(ab_out,message,'COLL')
!!$!    call wrtout(std_out,message,'COLL')
!!$!    ENDDEBUG
!!$
!!$!    Loop and keep the 2 largest values
!!$     if (max_perm(1)>max_perm(2)) then
!!$      maxp1=max_perm(1)
!!$      maxp2=max_perm(2)
!!$      iatom1=imax_perm(1)
!!$      iatom2=imax_perm(2)
!!$     else
!!$      maxp1=max_perm(2)
!!$      maxp2=max_perm(1)
!!$      iatom1=imax_perm(2)
!!$      iatom2=imax_perm(1)
!!$     end if
!!$
!!$     do itypat=3,ntypat
!!$      if (max_perm(itypat)>maxp1) then
!!$       maxp2=maxp1
!!$       iatom2=iatom1
!!$       maxp1=max_perm(itypat)
!!$       iatom1=imax_perm(itypat)
!!$      else if (max_perm(itypat)>maxp2) then
!!$       maxp2=max_perm(itypat)
!!$       iatom2=imax_perm(itypat)
!!$      end if
!!$     end do
!!$     write(message, '(2(a,i5))' )' Will commute atom...',iatom1,'...of type ',&
!!$&     typat(iatom1)
!!$     call wrtout(ab_out,message,'COLL')
!!$     call wrtout(std_out,message,'COLL')
!!$     write(message, '(2(a,i5))' )'         with atom...',iatom2,'...of type ',&
!!$&     typat(iatom2)
!!$     call wrtout(ab_out,message,'COLL')
!!$     call wrtout(std_out,message,'COLL')
!!$
!!$!    Commute the atoms positions
!!$     distx=xcart(1,iatom1)
!!$     disty=xcart(2,iatom1)
!!$     distz=xcart(3,iatom1)
!!$     xcart(1,iatom1)=xcart(1,iatom2)
!!$     xcart(2,iatom1)=xcart(2,iatom2)
!!$     xcart(3,iatom1)=xcart(3,iatom2)
!!$     xcart(1,iatom2)=distx
!!$     xcart(2,iatom2)=disty
!!$     xcart(3,iatom2)=distz
!!$!    Convert back to xred (reduced coordinates)
!!$     call xredxcart(natom,-1,rprimd,xcart,xred)
!!$
!!$!    Store current total energy
!!$     etotal_temp=etotal
!!$
!!$!    Compute LDA forces (big loop)
!!$     iapp=-1
!!$     if(itime>0)iapp=itime
!!$     call scfcv(acell,atindx,atindx1,cg,cpus,densymop_gs,dtefield,dtfil,&
!!$&     dtset,ecore,eigen,hdr,iapp,indsym,initialized,&
!!$&     irrzon,kg,mpi_enreg,&
!!$&     nattyp,nfftf,npwarr,nspinor,occ,&
!!$&     pawang,pawfgr,pawrad,pawrhoij,pawtab,phnons,psps,&
!!$&     pwind,pwind_alloc,pwnsfac,rec_set,resid,results_gs,rhog,rhor,rprimd,&
!!$&     scf_history,symrec,wffnew,wffnow,wvl,xred,xred_old,ylm,ylmgr)
!!$
!!$     if (etotal > etotal_temp) then
!!$
!!$!     Discard the changes
!!$      distx=xcart(1,iatom1)
!!$      disty=xcart(2,iatom1)
!!$      distz=xcart(3,iatom1)
!!$      xcart(1,iatom1)=xcart(1,iatom2)
!!$      xcart(2,iatom1)=xcart(2,iatom2)
!!$      xcart(3,iatom1)=xcart(3,iatom2)
!!$      xcart(1,iatom2)=distx
!!$      xcart(2,iatom2)=disty
!!$      xcart(3,iatom2)=distz
!!$
!!$!     Convert back to xred (reduced coordinates)
!!$      call xredxcart(natom,-1,rprimd,xcart,xred)
!!$      write(message, '(a)' )' Commutation unsuccessful, recomputing the forces'
!!$      call wrtout(ab_out,message,'COLL')
!!$      call wrtout(std_out,message,'COLL')
!!$
!!$!     And recompute the forces
!!$      iapp=-1
!!$      if(itime>0)iapp=itime
!!$      call scfcv(acell,atindx,atindx1,cg,cpus,densymop_gs,dtefield,dtfil,&
!!$&      dtset,ecore,eigen,hdr,iapp,indsym,initialized,&
!!$&      irrzon,kg,mpi_enreg,&
!!$&      nattyp,nfftf,npwarr,nspinor,occ,&
!!$&      pawang,pawfgr,pawrad,pawrhoij,pawtab,phnons,psps,&
!!$&      pwind,pwind_alloc,pwnsfac,rec_set,resid,results_gs,rhog,rhor,rprimd,&
!!$&      scf_history,symrec,wffnew,wffnow,wvl,xred,xred_old,ylm,ylmgr)
!!$
!!$     else
!!$
!!$      write(message, '(a)')' Commutation successful ! Going on'
!!$      call wrtout(ab_out,message,'COLL')
!!$      call wrtout(std_out,message,'COLL')
!!$
!!$!     Get rid of mean force on whole unit cell, but only if no generalized
!!$!     constraints are in effect
!!$      if(dtset%nconeq==0)then
!!$       do idir=1,3
!!$        favg=sum(fred(idir,:))/dble(natom)
!!$        fred_corrected(idir,:)=fred(idir,:)-favg
!!$        if(dtset%jellslab/=0.and.idir==3) fred_corrected(idir,:)=fred(idir,:)
!!$       end do
!!$      else
!!$       fred_corrected(:,:)=fred(:,:)
!!$      end if
!!$
!!$!     Update xfhist
!!$      xfhist(:,1:natom,1,nxfh)=xred(:,:)
!!$      xfhist(:,natom+1,1,nxfh)=acell(:)
!!$      xfhist(:,natom+2:natom+4,1,nxfh)=rprim(1:3,1:3)
!!$      xfhist(:,1:natom,2,nxfh)=fred_corrected(:,:)
!!$      xfhist(:,natom+2,2,nxfh)=strten(1:3)
!!$      xfhist(:,natom+3,2,nxfh)=strten(4:6)
!!$
!!$!     Store computed gradient in vout
!!$      option=3
!!$      call xfpack(acell,acell0,fred_corrected,&
!!$&      natom,ndim,nsym,optcell,option,rprim,rprimd0,&
!!$&      strtarget,strten,symrel,ucvol,ucvol0,vin,vout,xred)
!!$
!!$     end if
!!$
!!$
!!$    end if ! if(mod(itime,dtset%delayperm)==0)
!!$   end if ! if(itime/=0 .or. dtset%delayperm/=0)
!!$!  End of the commutation section

  allocate(ran_force(3,natom))

  !  Specific to Langevin dynamics
  !  Initialize an array of random forces
  !  No random force at itime=0
  !  if (itime==0) then
  if (itime<0) then

     ran_force(:,:)=0.0_dp

  else

     do iatom=1,natom
        !    sig_gauss is the std deviation of the random distribution
        sig_gauss=sqrt(2.0_dp*friction*amass(iatom)*ktemp)
        do idim=1,3
           delxi=2.0_dp
           do while (delxi >= 1.0_dp)
              ran_num1=2.0_dp*uniformrandom(idum)-1.0_dp
              ran_num2=2.0_dp*uniformrandom(idum)-1.0_dp
              delxi=ran_num1*ran_num1+ran_num2*ran_num2
           end do
           ran_force(idim,iatom)=ran_num1*sqrt(-2.0_dp*log(delxi)/delxi)&
                &      *sig_gauss/sqrt(dtion)

        end do
     end do
     !   DEBUG
     !   The distribution should be gaussian
     !   delxi=0.0_dp
     !   do iatom=1,natom
     !   do idim=1,3
     !   delxi=delxi+(ran_force(idim,iatom)*dtion)**2
     !   end do
     !   end do
     !   delxi=delxi/(3.0_dp*natom)
     !   write(message, '(2(a,es22.14))' )' variance =',delxi,'  asked =',&
     !   &    2.0_dp*(friction)*amass(2)*ktemp*dtion
     !   call wrtout(ab_out,message,'COLL')
     !   call wrtout(std_out,message,'COLL')
     !   ENDDEBUG
     !   end if itime\=0

  end if

  !  DEBUG
  !  write(message, '(a)' )' after initializing ran_force'
  !  call wrtout(ab_out,message,'COLL')
  !  call wrtout(std_out,message,'COLL')
  !  ENDDEBUG

  allocate(lang_force(3, natom))
  do iatom=1,natom
     do idim=1,3
        lang_force(idim,iatom)=fcart(idim,iatom)/amass(iatom)
        ran_force(idim,iatom)=ran_force(idim,iatom)/amass(iatom)
     end do
  end do
  lang_force(:,:)=ran_force(:,:)-(friction)*vel(:,:)+lang_force(:,:)

  deallocate(ran_force)

  !  DEBUG
  !  write(message, '(a)' )'before verlet'
  !  call wrtout(ab_out,message,'COLL')
  !  call wrtout(std_out,message,'COLL')
  !  ENDDEBUG

  !  Compute next atomic coordinates using Verlet algorithm

  !  Uses the velocity
  !  
  !  If an atom wants to cross the walls, velocity is reversed.
  !  
  do iatom=1,natom
     do idim=1,3
        delxi=xcart(idim,iatom)+dtion*vel(idim,iatom)+ &
             &     0.5_dp*dtion*dtion*lang_force(idim,iatom)
        if ( (delxi > (rprimd(idim,idim)+(mdwall)) ) .or. &
             &     (delxi < - (mdwall)                   )       ) then
           vel(idim,iatom)=-vel(idim,iatom)
           delxi=xcart(idim,iatom)+dtion*vel(idim,iatom)+ &
                &      0.5_dp*dtion*dtion*lang_force(idim,iatom)
        end if
        xcart_next(idim,iatom)=delxi
     end do
  end do

  !  Convert back to xred_next (reduced coordinates)
  call xredxcart(natom,-1,rprimd,xcart_next,xred_next)

  if (itime==0) then
     !   no old forces are available at first step
     !   Simple update of the velocity
     !   first compute vel_nexthalf for next steps
     vel(:,:)=vel(:,:)+dtion*lang_force(:,:)
  else
     !   case itime /= 0 normal verlet integration
     vel(:,:)=vel(:,:)+0.5_dp*dtion*(fcart_mold(:,:)+lang_force(:,:))
  end if

  !  Store 'current force' as 'old force'
  fcart_mold(:,:)=lang_force(:,:)

  deallocate(lang_force)

  !  End of case ionmov =9
end subroutine md_langevin
