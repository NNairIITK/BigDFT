subroutine abi_md_langevin(amass, dtion, fcart, fcart_mold, friction, itime, ktemp, &
     & mditemp, mdwall, natom, rprimd, vel, xcart, xcart_next, xred_next)
  
  use abi_defs_basis
  use abi_interfaces_geometry

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
     function abi_uniformrandom(seed) 
       implicit none
       integer :: seed
       double precision :: abi_uniformrandom
     end function abi_uniformrandom
  end interface


  real(dp),parameter :: v2tol=tol8
  integer :: iatom, idim, idum=-5
  real(dp) :: delxi, ran_num1, ran_num2, sig_gauss
  real(dp), allocatable :: ran_force(:,:), lang_force(:,:)
 
  if(itime==0)then
     call abi_md_nose_init(amass, natom, mditemp, vel)
  end if

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
              ran_num1=2.0_dp*abi_uniformrandom(idum)-1.0_dp
              ran_num2=2.0_dp*abi_uniformrandom(idum)-1.0_dp
              delxi=ran_num1*ran_num1+ran_num2*ran_num2
           end do
           ran_force(idim,iatom)=ran_num1*sqrt(-2.0_dp*log(delxi)/delxi)&
                &      *sig_gauss/sqrt(dtion)

        end do
     end do

  end if

  allocate(lang_force(3, natom))
  do iatom=1,natom
     do idim=1,3
        lang_force(idim,iatom)=fcart(idim,iatom)/amass(iatom)
        ran_force(idim,iatom)=ran_force(idim,iatom)/amass(iatom)
     end do
  end do
  lang_force(:,:)=ran_force(:,:)-(friction)*vel(:,:)+lang_force(:,:)

  deallocate(ran_force)

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
  call abi_xredxcart(natom,-1,rprimd,xcart_next,xred_next)

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
end subroutine abi_md_langevin
