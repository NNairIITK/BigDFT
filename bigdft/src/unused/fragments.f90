  ! copy a given full orbs structure to minimal orbs structure
  subroutine orbs_to_min_orbs_copy(orbs,forbs)
    use module_types
    implicit none
    ! Calling arguments
    type(orbitals_data),intent(in):: orbs
    type(minimal_orbitals_data),intent(inout):: forbs

    ! Local variables
    integer:: iis1, iie1, iis2, iie2, i1, i2, istat, iall
    character(len=200) :: subname

    subname='orbs_to_min_orbs_copy'

    forbs%norb = orbs%norb
    forbs%norbp = orbs%norbp
    forbs%isorb = orbs%isorb

    if(associated(forbs%norb_par)) then
        iall=-product(shape(forbs%norb_par))*kind(forbs%norb_par)
        deallocate(forbs%norb_par, stat=istat)
        call memocc(istat, iall, 'forbs%norb_par', subname)
    end if
    if(associated(orbs%norb_par)) then
        iis1=lbound(orbs%norb_par,1)
        iie1=ubound(orbs%norb_par,1)
        iis2=lbound(orbs%norb_par,2)
        iie2=ubound(orbs%norb_par,2)
        allocate(forbs%norb_par(iis1:iie1,iis2:iie2), stat=istat)
        call memocc(istat, forbs%norb_par, 'forbs%norb_par', subname)
        do i1=iis1,iie1
           do i2 = iis2,iie2
            forbs%norb_par(i1,i2) = orbs%norb_par(i1,i2)
           end do
        end do
    end if

    if(associated(forbs%inwhichlocreg)) then
        iall=-product(shape(forbs%inwhichlocreg))*kind(forbs%inwhichlocreg)
        deallocate(forbs%inwhichlocreg, stat=istat)
        call memocc(istat, iall, 'forbs%inwhichlocreg', subname)
    end if
    if(associated(orbs%inwhichlocreg)) then
        iis1=lbound(orbs%inwhichlocreg,1)
        iie1=ubound(orbs%inwhichlocreg,1)
        allocate(forbs%inwhichlocreg(iis1:iie1), stat=istat)
        call memocc(istat, forbs%inwhichlocreg, 'forbs%inwhichlocreg', subname)
        do i1=iis1,iie1
            forbs%inwhichlocreg(i1) = orbs%inwhichlocreg(i1)
        end do
    end if

    if(associated(forbs%onwhichatom)) then
        iall=-product(shape(forbs%onwhichatom))*kind(forbs%onwhichatom)
        deallocate(forbs%onwhichatom, stat=istat)
        call memocc(istat, iall, 'forbs%onwhichatom', subname)
    end if
    if(associated(orbs%onwhichatom)) then
        iis1=lbound(orbs%onwhichatom,1)
        iie1=ubound(orbs%onwhichatom,1)
        allocate(forbs%onwhichatom(iis1:iie1), stat=istat)
        call memocc(istat, forbs%onwhichatom, 'forbs%onwhichatom', subname)
        do i1=iis1,iie1
            forbs%onwhichatom(i1) = orbs%onwhichatom(i1)
        end do
    end if

    if(associated(forbs%isorb_par)) then
        iall=-product(shape(forbs%isorb_par))*kind(forbs%isorb_par)
        deallocate(forbs%isorb_par, stat=istat)
        call memocc(istat, iall, 'forbs%isorb_par', subname)
    end if
    if(associated(orbs%isorb_par)) then
        iis1=lbound(orbs%isorb_par,1)
        iie1=ubound(orbs%isorb_par,1)
        allocate(forbs%isorb_par(iis1:iie1), stat=istat)
        call memocc(istat, forbs%isorb_par, 'forbs%isorb_par', subname)
        do i1=iis1,iie1
            forbs%isorb_par(i1) = orbs%isorb_par(i1)
        end do
    end if

    if(associated(forbs%ispot)) then
        iall=-product(shape(forbs%ispot))*kind(forbs%ispot)
        deallocate(forbs%ispot, stat=istat)
        call memocc(istat, iall, 'forbs%ispot', subname)
    end if
    if(associated(orbs%ispot)) then
        iis1=lbound(orbs%ispot,1)
        iie1=ubound(orbs%ispot,1)
        allocate(forbs%ispot(iis1:iie1), stat=istat)
        call memocc(istat, forbs%ispot, 'forbs%ispot', subname)
        do i1=iis1,iie1
            forbs%ispot(i1) = orbs%ispot(i1)
        end do
    end if

  end subroutine orbs_to_min_orbs_copy



  pure function frg_center(frag)
    implicit none
    type(system_fragment), intent(in) :: frag
    real(gp), dimension(3) :: frg_center
    !local variables
    integer :: iat

    frg_center=0.0_gp
    do iat=1,frag%astruct_frg%nat
       frg_center=frg_center+frag%astruct_frg%rxyz(:,iat)
    end do
    frg_center=frg_center/real(frag%astruct_frg%nat,gp)

  end function frg_center
