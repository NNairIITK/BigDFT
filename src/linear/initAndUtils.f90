!> @file 
!!   Initializations
!! @author
!!   Copyright (C) 2011-2012 BigDFT group 
!!   This file is distributed under the terms of the
!!   GNU General Public License, see ~/COPYING file
!!   or http://www.gnu.org/copyleft/gpl.txt .
!!   For the list of contributors, see ~/AUTHORS 

subroutine allocateBasicArraysInputLin(lin, ntypes)
  use module_base
  use module_types
  implicit none
  
  ! Calling arguments
  type(linearInputParameters),intent(inout) :: lin
  integer, intent(in) :: ntypes
  
  ! Local variables
  integer :: istat
  character(len=*),parameter :: subname='allocateBasicArrays'
  
  allocate(lin%norbsPerType(ntypes), stat=istat)
  call memocc(istat, lin%norbsPerType, 'lin%norbsPerType', subname)

  allocate(lin%potentialPrefac_ao(ntypes), stat=istat)
  call memocc(istat, lin%potentialPrefac_ao, 'lin%potentialPrefac_ao', subname)

  allocate(lin%potentialPrefac_lowaccuracy(ntypes), stat=istat)
  call memocc(istat, lin%potentialPrefac_lowaccuracy, 'lin%potentialPrefac_lowaccuracy', subname)

  allocate(lin%potentialPrefac_highaccuracy(ntypes), stat=istat)
  call memocc(istat, lin%potentialPrefac_highaccuracy, 'lin%potentialPrefac_highaccuracy', subname)
  
  !added a second dimension to include the low and high accuracy values
  allocate(lin%locrad_type(ntypes,2),stat=istat)
  call memocc(istat,lin%locrad_type,'lin%locrad_type',subname)

  allocate(lin%kernel_cutoff_FOE(ntypes), stat=istat)
  call memocc(istat, lin%kernel_cutoff_FOE, 'lin%kernel_cutoff_FOE', subname)

  allocate(lin%kernel_cutoff(ntypes), stat=istat)
  call memocc(istat, lin%kernel_cutoff, 'lin%kernel_cutoff', subname)

end subroutine allocateBasicArraysInputLin

subroutine deallocateBasicArraysInput(lin)
  use module_base
  use module_types
  implicit none
  
  ! Calling arguments
  type(linearinputParameters),intent(inout) :: lin
  
  ! Local variables
  integer :: i_stat,i_all
  character(len=*),parameter :: subname='deallocateBasicArrays'
 
  if(associated(lin%potentialPrefac_ao)) then
    i_all = -product(shape(lin%potentialPrefac_ao))*kind(lin%potentialPrefac_ao)
    deallocate(lin%potentialPrefac_ao,stat=i_stat)
    call memocc(i_stat,i_all,'lin%potentialPrefac_ao',subname)
    nullify(lin%potentialPrefac_ao)
  end if 
  if(associated(lin%potentialPrefac_lowaccuracy)) then
    i_all = -product(shape(lin%potentialPrefac_lowaccuracy))*kind(lin%potentialPrefac_lowaccuracy)
    deallocate(lin%potentialPrefac_lowaccuracy,stat=i_stat)
    call memocc(i_stat,i_all,'lin%potentialPrefac_lowaccuracy',subname)
    nullify(lin%potentialPrefac_lowaccuracy)
  end if 
  if(associated(lin%potentialPrefac_highaccuracy)) then
    i_all = -product(shape(lin%potentialPrefac_highaccuracy))*kind(lin%potentialPrefac_highaccuracy)
    deallocate(lin%potentialPrefac_highaccuracy,stat=i_stat)
    call memocc(i_stat,i_all,'lin%potentialPrefac_highaccuracy',subname)
    nullify(lin%potentialPrefac_highaccuracy)
  end if 

  if(associated(lin%norbsPerType)) then
    i_all = -product(shape(lin%norbsPerType))*kind(lin%norbsPerType)
    deallocate(lin%norbsPerType,stat=i_stat)
    call memocc(i_stat,i_all,'lin%norbsPerType',subname)
    nullify(lin%norbsPerType)
  end if 

  if(associated(lin%locrad)) then
    i_all = -product(shape(lin%locrad))*kind(lin%locrad)
    deallocate(lin%locrad,stat=i_stat)
    call memocc(i_stat,i_all,'lin%locrad',subname)
    nullify(lin%locrad)
  end if 

  if(associated(lin%locrad_kernel)) then
    i_all = -product(shape(lin%locrad_kernel))*kind(lin%locrad_kernel)
    deallocate(lin%locrad_kernel,stat=i_stat)
    call memocc(i_stat,i_all,'lin%locrad_kernel',subname)
    nullify(lin%locrad_kernel)
  end if 

  if(associated(lin%locrad_mult)) then
    i_all = -product(shape(lin%locrad_mult))*kind(lin%locrad_mult)
    deallocate(lin%locrad_mult,stat=i_stat)
    call memocc(i_stat,i_all,'lin%locrad_mult',subname)
    nullify(lin%locrad_mult)
  end if 

  if(associated(lin%locrad_lowaccuracy)) then
    i_all = -product(shape(lin%locrad_lowaccuracy))*kind(lin%locrad_lowaccuracy)
    deallocate(lin%locrad_lowaccuracy,stat=i_stat)
    call memocc(i_stat,i_all,'lin%locrad_lowaccuracy',subname)
    nullify(lin%locrad_lowaccuracy)
  end if 

  if(associated(lin%locrad_highaccuracy)) then
    i_all = -product(shape(lin%locrad_highaccuracy))*kind(lin%locrad_highaccuracy)
    deallocate(lin%locrad_highaccuracy,stat=i_stat)
    call memocc(i_stat,i_all,'lin%locrad_highaccuracy',subname)
    nullify(lin%locrad_highaccuracy)
  end if 

  if(associated(lin%locrad_type)) then
    i_all = -product(shape(lin%locrad_type))*kind(lin%locrad_type)
    deallocate(lin%locrad_type,stat=i_stat)
    call memocc(i_stat,i_all,'lin%locrad_type',subname)
    nullify(lin%locrad_type)
  end if 

  if(associated(lin%kernel_cutoff_FOE)) then
    i_all = -product(shape(lin%kernel_cutoff_FOE))*kind(lin%kernel_cutoff_FOE)
    deallocate(lin%kernel_cutoff_FOE,stat=i_stat)
    call memocc(i_stat,i_all,'lin%kernel_cutoff_FOE',subname)
    nullify(lin%kernel_cutoff_FOE)
  end if 

  if(associated(lin%kernel_cutoff)) then
    i_all = -product(shape(lin%kernel_cutoff))*kind(lin%kernel_cutoff)
    deallocate(lin%kernel_cutoff,stat=i_stat)
    call memocc(i_stat,i_all,'lin%kernel_cutoff',subname)
    nullify(lin%kernel_cutoff)
  end if 

end subroutine deallocateBasicArraysInput



! lzd%llr already allocated, locregcenter and locrad already filled - could tidy this!
subroutine initLocregs(iproc, nproc, lzd, hx, hy, hz, astruct, orbs, Glr, locregShape, lborbs)
  use module_base
  use module_types
  use module_atoms, only: atomic_structure
  use module_interfaces, exceptThisOne => initLocregs
  implicit none
  
  ! Calling arguments
  integer,intent(in) :: iproc, nproc
  type(local_zone_descriptors),intent(inout) :: lzd
  real(kind=8),intent(in) :: hx, hy, hz
  type(atomic_structure),intent(in) :: astruct
  type(orbitals_data),intent(in) :: orbs
  type(locreg_descriptors),intent(in) :: Glr
  character(len=1),intent(in) :: locregShape
  type(orbitals_data),optional,intent(in) :: lborbs
  
  ! Local variables
  integer :: istat, jorb, jjorb, jlr, iall
  character(len=*),parameter :: subname='initLocregs'
  logical,dimension(:),allocatable :: calculateBounds

  
  allocate(calculateBounds(lzd%nlr), stat=istat)
  call memocc(istat, calculateBounds, 'calculateBounds', subname)
  calculateBounds=.false.
  
  do jorb=1,orbs%norbp
     jjorb=orbs%isorb+jorb
     jlr=orbs%inWhichLocreg(jjorb)
     calculateBounds(jlr)=.true.
  end do

  if(present(lborbs)) then
     do jorb=1,lborbs%norbp
        jjorb=lborbs%isorb+jorb
        jlr=lborbs%inWhichLocreg(jjorb)
        calculateBounds(jlr)=.true.
     end do
  end if
  
  if(locregShape=='c') then
      stop 'locregShape c is deprecated'
  else if(locregShape=='s') then
      call determine_locregSphere_parallel(iproc, nproc, lzd%nlr, hx, hy, hz, &
           astruct, orbs, Glr, lzd%Llr, calculateBounds)
  end if
  
  iall=-product(shape(calculateBounds))*kind(calculateBounds)
  deallocate(calculateBounds, stat=istat)
  call memocc(istat, iall, 'calculateBounds', subname)
  
  !DEBUG
  !do ilr=1,lin%nlr
  !    if(iproc==0) write(*,'(1x,a,i0)') '>>>>>>> zone ', ilr
  !    if(iproc==0) write(*,'(3x,a,4i10)') 'nseg_c, nseg_f, nvctr_c, nvctr_f', lin%Llr(ilr)%wfd%nseg_c, lin%Llr(ilr)%wfd%nseg_f, lin%Llr(ilr)%wfd%nvctr_c, lin%Llr(ilr)%wfd%nvctr_f
  !    if(iproc==0) write(*,'(3x,a,3i8)') 'lin%Llr(ilr)%d%n1i, lin%Llr(ilr)%d%n2i, lin%Llr(ilr)%d%n3i', lin%Llr(ilr)%d%n1i, lin%Llr(ilr)%d%n2i, lin%Llr(ilr)%d%n3i
  !    if(iproc==0) write(*,'(a,6i8)') 'lin%Llr(ilr)%d%nfl1,lin%Llr(ilr)%d%nfu1,lin%Llr(ilr)%d%nfl2,lin%Llr(ilr)%d%nfu2,lin%Llr(ilr)%d%nfl3,lin%Llr(ilr)%d%nfu3',&
  !    lin%Llr(ilr)%d%nfl1,lin%Llr(ilr)%d%nfu1,lin%Llr(ilr)%d%nfl2,lin%Llr(ilr)%d%nfu2,lin%Llr(ilr)%d%nfl3,lin%Llr(ilr)%d%nfu3
  !end do
  !END DEBUG
  
  lzd%linear=.true.

end subroutine initLocregs



subroutine init_foe(iproc, nproc, lzd, astruct, input, orbs_KS, orbs, foe_obj, reset, &
           cutoff_incr)
  use module_base
  use module_atoms, only: atomic_structure
  use module_types
  implicit none
  
  ! Calling arguments
  integer,intent(in) :: iproc, nproc
  type(local_zone_descriptors),intent(in) :: lzd
  type(atomic_structure),intent(in) :: astruct
  type(input_variables),intent(in) :: input
  type(orbitals_data),intent(in) :: orbs_KS, orbs
  type(foe_data),intent(out) :: foe_obj
  logical, intent(in) :: reset
  real(kind=8),optional,intent(in) :: cutoff_incr
  
  ! Local variables
  integer :: iorb, iiorb, jjorb, istat, iseg, ilr, jlr
  integer :: iwa, jwa, itype, jtype, ierr, iall, isegstart
  logical :: seg_started
  real(kind=8) :: tt, cut, incr
  logical,dimension(:,:),allocatable :: kernel_locreg
  character(len=*),parameter :: subname='initMatrixCompression'
!  integer :: ii, iseg

  if (present(cutoff_incr)) then
      incr=cutoff_incr
  else
      incr=0.d0
  end if
  
  call timing(iproc,'init_matrCompr','ON')

  if (reset) then
     foe_obj%ef=0.d0
     foe_obj%evlow=input%lin%evlow
     foe_obj%evhigh=input%lin%evhigh
     foe_obj%bisection_shift=1.d-1
     foe_obj%fscale=input%lin%fscale
     foe_obj%ef_interpol_det=input%lin%ef_interpol_det
     foe_obj%ef_interpol_chargediff=input%lin%ef_interpol_chargediff
     foe_obj%charge=0.d0
     do iorb=1,orbs_KS%norb
          foe_obj%charge=foe_obj%charge+orbs_KS%occup(iorb)
     end do
     foe_obj%evbounds_isatur=0
     foe_obj%evboundsshrink_isatur=0
     foe_obj%evbounds_nsatur=input%evbounds_nsatur
     foe_obj%evboundsshrink_nsatur=input%evboundsshrink_nsatur
     foe_obj%fscale_lowerbound=input%fscale_lowerbound
     foe_obj%fscale_upperbound=input%fscale_upperbound
  end if

  call nullify_foe(foe_obj)

  ! Initialize kernel_locreg
  !if (input%lin%scf_mode==LINEAR_FOE) then ! otherwise don't need to allocate just nullify as above
     allocate(kernel_locreg(orbs%norbp,orbs%norb), stat=istat)
     call memocc(istat, kernel_locreg, 'kernel_locreg', subname)
     allocate(foe_obj%nsegline(orbs%norb), stat=istat)
     call memocc(istat, foe_obj%nsegline, 'foe_obj%nsegline', subname)
     call to_zero(orbs%norb, foe_obj%nsegline(1))
     allocate(foe_obj%istsegline(orbs%norb), stat=istat)
     call memocc(istat, foe_obj%istsegline, 'foe_obj%nsegline', subname)
     call to_zero(orbs%norb, foe_obj%istsegline(1))
     do iorb=1,orbs%norbp
        iiorb=orbs%isorb+iorb
        ilr=orbs%inwhichlocreg(iiorb)
        iwa=orbs%onwhichatom(iiorb)
        itype=astruct%iatype(iwa)
        foe_obj%nsegline(iiorb)=0
        seg_started=.false.
        do jjorb=1,orbs%norb
           jlr=orbs%inwhichlocreg(jjorb)
           jwa=orbs%onwhichatom(jjorb)
           jtype=astruct%iatype(jwa)
           tt = (lzd%llr(ilr)%locregcenter(1)-lzd%llr(jlr)%locregcenter(1))**2 + &
                (lzd%llr(ilr)%locregcenter(2)-lzd%llr(jlr)%locregcenter(2))**2 + &
                (lzd%llr(ilr)%locregcenter(3)-lzd%llr(jlr)%locregcenter(3))**2
           cut = input%lin%kernel_cutoff_FOE(itype)+input%lin%kernel_cutoff_FOE(jtype)+2.d0*incr
           tt=sqrt(tt)
           if (tt<=cut) then
              kernel_locreg(iorb,jjorb)=.true.
              if (.not.seg_started) then
                 foe_obj%nsegline(iiorb)=foe_obj%nsegline(iiorb)+1
              end if
              seg_started=.true.
           else
              kernel_locreg(iorb,jjorb)=.false.
              seg_started=.false.
           end if
        end do
     end do

     if (nproc > 1) then
         call mpiallred(foe_obj%nsegline(1), orbs%norb, mpi_sum, bigdft_mpi%mpi_comm)
     end if

     ! Total number of segments
     foe_obj%nseg = sum(foe_obj%nsegline)
     
     ! Initialize istsegline, which gives the first segment of each line
     foe_obj%istsegline(1)=1
     do iorb=2,orbs%norb
         foe_obj%istsegline(iorb) = foe_obj%istsegline(iorb-1) + foe_obj%nsegline(iorb-1)
     end do

     allocate(foe_obj%keyg(2,foe_obj%nseg),stat=istat)
     call memocc(istat, foe_obj%keyg, 'foe_obj%keyg', subname)
     call to_zero(2*foe_obj%nseg, foe_obj%keyg(1,1))

     do iorb=1,orbs%norbp
        iiorb=orbs%isorb+iorb
        iseg=0
        seg_started=.false.
        isegstart=foe_obj%istsegline(iiorb)-1
        do jjorb=1,orbs%norb
           if(kernel_locreg(iorb,jjorb)) then
              if (.not.seg_started) then
                 iseg=iseg+1
                 foe_obj%keyg(1,isegstart+iseg)=(iiorb-1)*orbs%norb+jjorb
              end if
              seg_started=.true.
           else
              if (seg_started) then
                 foe_obj%keyg(2,isegstart+iseg)=(iiorb-1)*orbs%norb+jjorb-1
              end if
              seg_started=.false.
           end if
        end do
        if (seg_started) then
           foe_obj%keyg(2,isegstart+iseg)=(iiorb-1)*orbs%norb+orbs%norb
        end if
     end do

     if (nproc > 1) then
         call mpiallred(foe_obj%keyg(1,1), 2*foe_obj%nseg, mpi_sum, bigdft_mpi%mpi_comm)
     end if

     iall = -product(shape(kernel_locreg))*kind(kernel_locreg) 
     deallocate(kernel_locreg,stat=istat)
     call memocc(istat,iall,'kernel_locreg',subname)
  !end if

  call timing(iproc,'init_matrCompr','OF')


end subroutine init_foe


subroutine check_linear_and_create_Lzd(iproc,nproc,linType,Lzd,atoms,orbs,nspin,rxyz)
  use module_base
  use module_types
  use module_xc
  use ao_inguess, only: atomic_info
  implicit none

  integer, intent(in) :: iproc,nproc,nspin
  type(local_zone_descriptors), intent(inout) :: Lzd
  type(atoms_data), intent(in) :: atoms
  type(orbitals_data),intent(inout) :: orbs
  real(gp), dimension(3,atoms%astruct%nat), intent(in) :: rxyz
  integer, intent(in) :: linType
!  real(gp), dimension(atoms%astruct%ntypes,3), intent(in) :: radii_cf
  !Local variables
  character(len=*), parameter :: subname='check_linear_and_create_Lzd'
  logical :: linear
  real(gp) :: rcov
  integer :: iat,ityp,nspin_ig,i_all,i_stat,ilr
  real(gp), dimension(:), allocatable :: locrad
  logical,dimension(:),allocatable :: calculateBounds

  !default variables
  Lzd%nlr = 1

  if (nspin == 4) then
     nspin_ig=1
  else
     nspin_ig=nspin
  end if

  linear  = .true.
  if (linType == INPUT_IG_FULL) then
     Lzd%nlr=atoms%astruct%nat
     allocate(locrad(Lzd%nlr+ndebug),stat=i_stat)
     call memocc(i_stat,locrad,'locrad',subname)
     ! locrad read from last line of  psppar
     do iat=1,atoms%astruct%nat
        ityp = atoms%astruct%iatype(iat)
        call atomic_info(atoms%nzatom(ityp),atoms%nelpsp(ityp),rcov=rcov)
        locrad(iat) =  rcov * 10.0_gp ! locrad(iat) = atoms%rloc(ityp,1)
     end do  
     call timing(iproc,'check_IG      ','ON')
     call check_linear_inputguess(iproc,Lzd%nlr,rxyz,locrad,&
          Lzd%hgrids(1),Lzd%hgrids(2),Lzd%hgrids(3),&
          Lzd%Glr,linear) 
     call timing(iproc,'check_IG      ','OF')
     if(nspin >= 4) linear = .false. 
  end if

  ! If we are using cubic code : by choice or because locregs are too big
  Lzd%linear = .true.
  if (linType == INPUT_IG_LIG .or. linType == INPUT_IG_OFF .or. .not. linear) then
     Lzd%linear = .false.
     Lzd%nlr = 1
  end if


  if(linType /= INPUT_IG_TMO) then
     allocate(Lzd%Llr(Lzd%nlr+ndebug))
     do ilr=1,Lzd%nlr
        Lzd%Llr(ilr)=locreg_null()
     end do
     !for now, always true because we want to calculate the hamiltonians for all locregs
     if(.not. Lzd%linear) then
        Lzd%lintyp = 0
        !copy Glr to Llr(1)
        call nullify_locreg_descriptors(Lzd%Llr(1))
        call copy_locreg_descriptors(Lzd%Glr,Lzd%Llr(1))
     else 
        Lzd%lintyp = 1
        ! Assign orbitals to locreg (for LCAO IG each orbitals corresponds to an atomic function. WILL NEED TO CHANGE THIS)
        call assignToLocreg(iproc,nproc,orbs%nspinor,nspin_ig,atoms,orbs,Lzd)

        ! determine the localization regions
        ! calculateBounds indicate whether the arrays with the bounds (for convolutions...) shall also
        ! be allocated and calculated. In principle this is only necessary if the current process has orbitals
        ! in this localization region.
        allocate(calculateBounds(lzd%nlr),stat=i_stat)
        call memocc(i_stat,calculateBounds,'calculateBounds',subname)
        calculateBounds=.true.
!        call determine_locreg_periodic(iproc,Lzd%nlr,rxyz,locrad,hx,hy,hz,Lzd%Glr,Lzd%Llr,calculateBounds)
        call determine_locreg_parallel(iproc,nproc,Lzd%nlr,rxyz,locrad,&
             Lzd%hgrids(1),Lzd%hgrids(2),Lzd%hgrids(3),Lzd%Glr,Lzd%Llr,&
             orbs,calculateBounds)  
        i_all = -product(shape(calculateBounds))*kind(calculateBounds) 
        deallocate(calculateBounds,stat=i_stat)
        call memocc(i_stat,i_all,'calculateBounds',subname)
        i_all = -product(shape(locrad))*kind(locrad)
        deallocate(locrad,stat=i_stat)
        call memocc(i_stat,i_all,'locrad',subname)

        ! determine the wavefunction dimension
        call wavefunction_dimension(Lzd,orbs)
     end if
  else
     Lzd%lintyp = 2
  end if
  
!DEBUG
!!if(iproc==0)then
!!print *,'###################################################'
!!print *,'##        General information:                   ##'
!!print *,'###################################################'
!!print *,'Lzd%nlr,linear, ndimpotisf :',Lzd%nlr,Lzd%linear,Lzd%ndimpotisf
!!print *,'###################################################'
!!print *,'##        Global box information:                ##'
!!print *,'###################################################'
!!write(*,'(a24,3i4)')'Global region n1,n2,n3:',Lzd%Glr%d%n1,Lzd%Glr%d%n2,Lzd%Glr%d%n3
!!write(*,*)'Global fine grid: nfl',Lzd%Glr%d%nfl1,Lzd%Glr%d%nfl2,Lzd%Glr%d%nfl3
!!write(*,*)'Global fine grid: nfu',Lzd%Glr%d%nfu1,Lzd%Glr%d%nfu2,Lzd%Glr%d%nfu3
!!write(*,*)'Global inter. grid: ni',Lzd%Glr%d%n1i,Lzd%Glr%d%n2i,Lzd%Glr%d%n3i
!!write(*,'(a27,f6.2,f6.2,f6.2)')'Global dimension (1x,y,z):',Lzd%Glr%d%n1*hx,Lzd%Glr%d%n2*hy,Lzd%Glr%d%n3*hz
!!write(*,'(a17,f12.2)')'Global volume: ',Lzd%Glr%d%n1*hx*Lzd%Glr%d%n2*hy*Lzd%Glr%d%n3*hz
!!print *,'Global wfd statistics:',Lzd%Glr%wfd%nseg_c,Lzd%Glr%wfd%nseg_f,Lzd%Glr%wfd%nvctr_c,Lzd%Glr%wfd%nvctr_f
!!print *,'###################################################'
!!print *,'##        Local boxes information:               ##'
!!print *,'###################################################'
!!do i_stat =1, Lzd%nlr
!!   write(*,*)'=====> Region:',i_stat
!!   write(*,'(a24,3i4)')'Local region n1,n2,n3:',Lzd%Llr(i_stat)%d%n1,Lzd%Llr(i_stat)%d%n2,Lzd%Llr(i_stat)%d%n3
!!   write(*,*)'Local fine grid: nfl',Lzd%Llr(i_stat)%d%nfl1,Lzd%Llr(i_stat)%d%nfl2,Lzd%Llr(i_stat)%d%nfl3
!!   write(*,*)'Local fine grid: nfu',Lzd%Llr(i_stat)%d%nfu1,Lzd%Llr(i_stat)%d%nfu2,Lzd%Llr(i_stat)%d%nfu3
!!   write(*,*)'Local inter. grid: ni',Lzd%Llr(i_stat)%d%n1i,Lzd%Llr(i_stat)%d%n2i,Lzd%Llr(i_stat)%d%n3i
!!   write(*,'(a27,f6.2,f6.2,f6.2)')'Local dimension (1x,y,z):',Lzd%Llr(i_stat)%d%n1*hx,Lzd%Llr(i_stat)%d%n2*hy,&
!!            Lzd%Llr(i_stat)%d%n3*hz
!!   write(*,'(a17,f12.2)')'Local volume: ',Lzd%Llr(i_stat)%d%n1*hx*Lzd%Llr(i_stat)%d%n2*hy*Lzd%Llr(i_stat)%d%n3*hz
!!   print *,'Local wfd statistics:',Lzd%Llr(i_stat)%wfd%nseg_c,Lzd%Llr(i_stat)%wfd%nseg_f,Lzd%Llr(i_stat)%wfd%nvctr_c,&
!!            Lzd%Llr(i_stat)%wfd%nvctr_f
!!end do
!!end if
!!call mpi_finalize(i_stat)
!!stop
!END DEBUG

end subroutine check_linear_and_create_Lzd

subroutine create_LzdLIG(iproc,nproc,nspin,linearmode,hx,hy,hz,Glr,atoms,orbs,rxyz,nl,Lzd)
  use module_base
  use module_types
  use module_xc
  use ao_inguess, only: atomic_info
  implicit none

  integer, intent(in) :: iproc,nproc,nspin
  real(gp), intent(in) :: hx,hy,hz
  type(locreg_descriptors), intent(in) :: Glr
  type(atoms_data), intent(in) :: atoms
  type(orbitals_data),intent(inout) :: orbs
  integer, intent(in) :: linearmode
  real(gp), dimension(3,atoms%astruct%nat), intent(in) :: rxyz
  type(local_zone_descriptors), intent(inout) :: Lzd
  type(DFT_PSP_projectors), intent(inout) :: nl
  !  real(gp), dimension(atoms%astruct%ntypes,3), intent(in) :: radii_cf
  !Local variables
  character(len=*), parameter :: subname='check_linear_and_create_Lzd'
  logical :: linear
  integer :: iat,ityp,nspin_ig,ilr
  real(gp) :: rcov
  real(gp), dimension(:), allocatable :: locrad
  logical,dimension(:),allocatable :: calculateBounds,lr_mask

  call f_routine(id=subname)
  !default variables
  Lzd%nlr = 1

  Lzd%hgrids(1)=hx
  Lzd%hgrids(2)=hy
  Lzd%hgrids(3)=hz

  if (nspin == 4) then
     nspin_ig=1
  else
     nspin_ig=nspin
  end if

  linear  = .true.
  if (linearmode == INPUT_IG_LIG .or. linearmode == INPUT_IG_FULL) then
     Lzd%nlr=atoms%astruct%nat
     locrad=f_malloc(Lzd%nlr,id='locrad')
     ! locrad read from last line of  psppar
     do iat=1,atoms%astruct%nat
        ityp = atoms%astruct%iatype(iat)
        call atomic_info(atoms%nzatom(ityp),atoms%nelpsp(ityp),rcov=rcov)
        locrad(iat) =  rcov * 10.0_gp ! atoms%rloc(ityp,1)
        !locrad(iat)=18.d0
     end do
     call timing(iproc,'check_IG      ','ON')
     call check_linear_inputguess(iproc,Lzd%nlr,rxyz,locrad,hx,hy,hz,&
          Glr,linear) 
     call timing(iproc,'check_IG      ','OF')
     if(nspin >= 4) linear = .false. 
  end if

  ! If we are using cubic code : by choice or because locregs are too big
  if (linearmode == INPUT_IG_OFF .or. .not. linear) then
     linear = .false.
     Lzd%nlr = 1
  end if

  Lzd%linear = .true.
  if (.not. linear)  Lzd%linear = .false.

  !  print *,'before Glr => Lzd%Glr'
  call nullify_locreg_descriptors(Lzd%Glr)
  call copy_locreg_descriptors(Glr,Lzd%Glr)

  if(linearmode /= INPUT_IG_TMO) then
     allocate(Lzd%Llr(Lzd%nlr))
     do ilr=1,Lzd%nlr
        Lzd%Llr(ilr)=locreg_null()
     end do
     !for now, always true because we want to calculate the hamiltonians for all locregs

     if(.not. Lzd%linear) then
        Lzd%lintyp = 0
        call nullify_locreg_descriptors(Lzd%Llr(1))
        call copy_locreg_descriptors(Glr,Lzd%Llr(1))
     else 
        Lzd%lintyp = 1
        ! Assign orbitals to locreg (for LCAO IG each orbitals corresponds to an atomic function. WILL NEED TO CHANGE THIS)
        call assignToLocreg(iproc,nproc,orbs%nspinor,nspin_ig,atoms,orbs,Lzd)

        ! determine the localization regions
        ! calculateBounds indicate whether the arrays with the bounds (for convolutions...) shall also
        ! be allocated and calculated. In principle this is only necessary if the current process has orbitals
        ! in this localization region.
        calculateBounds=f_malloc(lzd%nlr,id='calculateBounds')
        calculateBounds=.true.
        !        call determine_locreg_periodic(iproc,Lzd%nlr,rxyz,locrad,hx,hy,hz,Glr,Lzd%Llr,calculateBounds)
        call determine_locreg_parallel(iproc,nproc,Lzd%nlr,rxyz,locrad,&
             hx,hy,hz,Glr,Lzd%Llr,&
             orbs,calculateBounds)  
        call f_free(calculateBounds)
        call f_free(locrad)

        ! determine the wavefunction dimension
        call wavefunction_dimension(Lzd,orbs)
        !in this case update the projector descriptor to be compatible with the locregs
        lr_mask=f_malloc0(Lzd%nlr,id='lr_mask')
        call update_lrmask_array(Lzd%nlr,orbs,lr_mask)
        !when the new tmbs are created the projector descriptors can be updated
        call update_nlpsp(nl,Lzd%nlr,Lzd%llr,Lzd%Glr,lr_mask)
        if (iproc == 0) call print_nlpsp(nl)
       
        call f_free(lr_mask)
     end if
  else
     Lzd%lintyp = 2
  end if

  call f_release_routine()

  !DEBUG
  !!if(iproc==0)then
  !!print *,'###################################################'
  !!print *,'##        General information:                   ##'
  !!print *,'###################################################'
  !!print *,'Lzd%nlr,linear, Lpsidimtot, ndimpotisf, Lnprojel:',Lzd%nlr,Lzd%linear,Lzd%ndimpotisf
  !!print *,'###################################################'
  !!print *,'##        Global box information:                ##'
  !!print *,'###################################################'
  !!write(*,'(a24,3i4)')'Global region n1,n2,n3:',Lzd%Glr%d%n1,Lzd%Glr%d%n2,Lzd%Glr%d%n3
  !!write(*,*)'Global fine grid: nfl',Lzd%Glr%d%nfl1,Lzd%Glr%d%nfl2,Lzd%Glr%d%nfl3
  !!write(*,*)'Global fine grid: nfu',Lzd%Glr%d%nfu1,Lzd%Glr%d%nfu2,Lzd%Glr%d%nfu3
  !!write(*,*)'Global inter. grid: ni',Lzd%Glr%d%n1i,Lzd%Glr%d%n2i,Lzd%Glr%d%n3i
  !!write(*,'(a27,f6.2,f6.2,f6.2)')'Global dimension (1x,y,z):',Lzd%Glr%d%n1*hx,Lzd%Glr%d%n2*hy,Lzd%Glr%d%n3*hz
  !!write(*,'(a17,f12.2)')'Global volume: ',Lzd%Glr%d%n1*hx*Lzd%Glr%d%n2*hy*Lzd%Glr%d%n3*hz
  !!print *,'Global wfd statistics:',Lzd%Glr%wfd%nseg_c,Lzd%Glr%wfd%nseg_f,Lzd%Glr%wfd%nvctr_c,Lzd%Glr%wfd%nvctr_f
  !!write(*,'(a17,f12.2)')'Global volume: ',Lzd%Glr%d%n1*input%hx*Lzd%Glr%d%n2*input%hy*Lzd%Glr%d%n3*input%hz
  !!print *,'Global wfd statistics:',Lzd%Glr%wfd%nseg_c,Lzd%Glr%wfd%nseg_f,Lzd%Glr%wfd%nvctr_c,Lzd%Glr%wfd%nvctr_f
  !!print *,'###################################################'
  !!print *,'##        Local boxes information:               ##'
  !!print *,'###################################################'
  !!do i_stat =1, Lzd%nlr
  !!   write(*,*)'=====> Region:',i_stat
  !!   write(*,'(a24,3i4)')'Local region n1,n2,n3:',Lzd%Llr(i_stat)%d%n1,Lzd%Llr(i_stat)%d%n2,Lzd%Llr(i_stat)%d%n3
  !!   write(*,*)'Local fine grid: nfl',Lzd%Llr(i_stat)%d%nfl1,Lzd%Llr(i_stat)%d%nfl2,Lzd%Llr(i_stat)%d%nfl3
  !!   write(*,*)'Local fine grid: nfu',Lzd%Llr(i_stat)%d%nfu1,Lzd%Llr(i_stat)%d%nfu2,Lzd%Llr(i_stat)%d%nfu3
  !!   write(*,*)'Local inter. grid: ni',Lzd%Llr(i_stat)%d%n1i,Lzd%Llr(i_stat)%d%n2i,Lzd%Llr(i_stat)%d%n3i
  !!   write(*,'(a27,f6.2,f6.2,f6.2)')'Local dimension (1x,y,z):',Lzd%Llr(i_stat)%d%n1*hx,Lzd%Llr(i_stat)%d%n2*hy,&
  !!            Lzd%Llr(i_stat)%d%n3*hz
  !!   write(*,'(a17,f12.2)')'Local volume: ',Lzd%Llr(i_stat)%d%n1*hx*Lzd%Llr(i_stat)%d%n2*hy*Lzd%Llr(i_stat)%d%n3*hz
  !!   print *,'Local wfd statistics:',Lzd%Llr(i_stat)%wfd%nseg_c,Lzd%Llr(i_stat)%wfd%nseg_f,Lzd%Llr(i_stat)%wfd%nvctr_c,&
  !!            Lzd%Llr(i_stat)%wfd%nvctr_f
  !!end do
  !!end if
  !call mpi_finalize(i_stat)
  !stop
  !END DEBUG

end subroutine create_LzdLIG





subroutine init_orbitals_data_for_linear(iproc, nproc, nspinor, input, astruct, rxyz, lorbs)
  use module_base
  use module_types
  use module_interfaces, except_this_one => init_orbitals_data_for_linear
  implicit none
  
  ! Calling arguments
  integer,intent(in) :: iproc, nproc, nspinor
  type(input_variables),intent(in) :: input
  type(atomic_structure),intent(in) :: astruct
  real(kind=8),dimension(3,astruct%nat),intent(in) :: rxyz
  type(orbitals_data),intent(out) :: lorbs
  
  ! Local variables
  integer :: norb, norbu, norbd, ityp, iat, ilr, istat, iall, iorb, nlr
  integer,dimension(:),allocatable :: norbsPerLocreg, norbsPerAtom
  real(kind=8),dimension(:,:),allocatable :: locregCenter
  character(len=*),parameter :: subname='init_orbitals_data_for_linear'

  call timing(iproc,'init_orbs_lin ','ON')
  
  call nullify_orbitals_data(lorbs)
 
  ! Count the number of basis functions.
  allocate(norbsPerAtom(astruct%nat), stat=istat)
  call memocc(istat, norbsPerAtom, 'norbsPerAtom', subname)
  norb=0
  nlr=0
  do iat=1,astruct%nat
      ityp=astruct%iatype(iat)
      norbsPerAtom(iat)=input%lin%norbsPerType(ityp)
      norb=norb+input%lin%norbsPerType(ityp)
      nlr=nlr+input%lin%norbsPerType(ityp)
  end do

  ! Distribute the basis functions among the processors.
  norbu=norb
  norbd=0
!!$  call orbitals_descriptors_forLinear(iproc, nproc, norb, norbu, norbd, input%nspin, nspinor,&
!!$       input%nkpt, input%kpt, input%wkpt, lorbs)
!!$  call repartitionOrbitals(iproc, nproc, lorbs%norb, lorbs%norb_par,&
!!$       lorbs%norbp, lorbs%isorb_par, lorbs%isorb, lorbs%onWhichMPI)
 
  call orbitals_descriptors(iproc, nproc, norb, norbu, norbd, input%nspin, nspinor,&
       input%gen_nkpt, input%gen_kpt, input%gen_wkpt, lorbs,.true.) !simple repartition

  allocate(locregCenter(3,nlr), stat=istat)
  call memocc(istat, locregCenter, 'locregCenter', subname)
  
  ilr=0
  do iat=1,astruct%nat
      ityp=astruct%iatype(iat)
      do iorb=1,input%lin%norbsPerType(ityp)
          ilr=ilr+1
          locregCenter(:,ilr)=rxyz(:,iat)
          ! DEBUGLR write(10,*) iorb,locregCenter(:,ilr)
      end do
  end do
 
  allocate(norbsPerLocreg(nlr), stat=istat)
  call memocc(istat, norbsPerLocreg, 'norbsPerLocreg', subname)
  norbsPerLocreg=1 !should be norbsPerLocreg
    
  iall=-product(shape(lorbs%inWhichLocreg))*kind(lorbs%inWhichLocreg)
  deallocate(lorbs%inWhichLocreg, stat=istat)
  call memocc(istat, iall, 'lorbs%inWhichLocreg', subname)
  call assignToLocreg2(iproc, nproc, lorbs%norb, lorbs%norb_par, astruct%nat, nlr, &
       input%nspin, norbsPerLocreg, locregCenter, lorbs%inwhichlocreg)

  iall=-product(shape(lorbs%onwhichatom))*kind(lorbs%onwhichatom)
  deallocate(lorbs%onwhichatom, stat=istat)
  call memocc(istat, iall, 'lorbs%onwhichatom', subname)
  call assignToLocreg2(iproc, nproc, lorbs%norb, lorbs%norb_par, astruct%nat, astruct%nat, &
       input%nspin, norbsPerAtom, rxyz, lorbs%onwhichatom)
  
  allocate(lorbs%eval(lorbs%norb), stat=istat)
  call memocc(istat, lorbs%eval, 'lorbs%eval', subname)
  lorbs%eval=-.5d0
  
  iall=-product(shape(norbsPerLocreg))*kind(norbsPerLocreg)
  deallocate(norbsPerLocreg, stat=istat)
  call memocc(istat, iall, 'norbsPerLocreg', subname)
  
  iall=-product(shape(locregCenter))*kind(locregCenter)
  deallocate(locregCenter, stat=istat)
  call memocc(istat, iall, 'locregCenter', subname)

  iall=-product(shape(norbsPerAtom))*kind(norbsPerAtom)
  deallocate(norbsPerAtom, stat=istat)
  call memocc(istat, iall, 'norbsPerAtom', subname)


  call timing(iproc,'init_orbs_lin ','OF')

end subroutine init_orbitals_data_for_linear


! initializes locrad and locregcenter
subroutine lzd_init_llr(iproc, nproc, input, astruct, rxyz, orbs, lzd)
  use module_base
  use module_types
  use module_interfaces
  implicit none
  
  ! Calling arguments
  integer,intent(in) :: iproc, nproc
  type(input_variables),intent(in) :: input
  type(atomic_structure),intent(in) :: astruct
  real(kind=8),dimension(3,astruct%nat),intent(in) :: rxyz
  type(orbitals_data),intent(in) :: orbs
  type(local_zone_descriptors),intent(inout) :: lzd
  
  ! Local variables
  integer :: iat, ityp, ilr, istat, iorb, iall
  real(kind=8),dimension(:,:),allocatable :: locregCenter
  character(len=*),parameter :: subname='lzd_init_llr'
  real(8):: t1, t2

  call timing(iproc,'init_locregs  ','ON')
  t1=mpi_wtime()
  
  nullify(lzd%llr)

  lzd%nlr=orbs%norb

  allocate(locregCenter(3,lzd%nlr), stat=istat)
  call memocc(istat, locregCenter, 'locregCenter', subname)
  
  ilr=0
  do iat=1,astruct%nat
      ityp=astruct%iatype(iat)
      do iorb=1,input%lin%norbsPerType(ityp)
          ilr=ilr+1
          locregCenter(:,ilr)=rxyz(:,iat)
      end do
  end do

  ! Allocate the array of localisation regions
  allocate(lzd%Llr(lzd%nlr),stat=istat)
  do ilr=1,lzd%nlr
     lzd%Llr(ilr)=locreg_null()
  end do
  do ilr=1,lzd%nlr
      lzd%llr(ilr)%locrad=input%lin%locrad(ilr)
      lzd%llr(ilr)%locrad_kernel=input%lin%locrad_kernel(ilr)
      lzd%llr(ilr)%locrad_mult=input%lin%locrad_mult(ilr)
      lzd%llr(ilr)%locregCenter=locregCenter(:,ilr)
  end do

  iall=-product(shape(locregCenter))*kind(locregCenter)
  deallocate(locregCenter, stat=istat)
  call memocc(istat, iall, 'locregCenter', subname)
  
  t2=mpi_wtime()
  !if(iproc==0) write(*,*) 'in lzd_init_llr: time',t2-t1
  call timing(iproc,'init_locregs  ','OF')

end subroutine lzd_init_llr


subroutine update_locreg(iproc, nproc, nlr, locrad, locrad_kernel, locrad_mult, locregCenter, glr_tmp, &
           useDerivativeBasisFunctions, nscatterarr, hx, hy, hz, astruct, input, &
           orbs_KS, orbs, lzd, npsidim_orbs, npsidim_comp, lbcomgp, lbcollcom, lfoe, lbcollcom_sr)
  use module_base
  use module_types
  use module_interfaces, except_this_one => update_locreg
  use communications_base, only: comms_linear_null
  use communications_init, only: init_comms_linear, init_comms_linear_sumrho, &
                                 initialize_communication_potential
  implicit none
  
  ! Calling arguments
  integer,intent(in) :: iproc, nproc, nlr
  integer,intent(out) :: npsidim_orbs, npsidim_comp
  logical,intent(in) :: useDerivativeBasisFunctions
  integer,dimension(0:nproc-1,4),intent(in) :: nscatterarr !n3d,n3p,i3s+i3xcsh-1,i3xcsh
  real(kind=8),intent(in) :: hx, hy, hz
  type(atomic_structure),intent(in) :: astruct
  type(input_variables),intent(in) :: input
  real(kind=8),dimension(nlr),intent(in) :: locrad, locrad_kernel, locrad_mult
  type(orbitals_data),intent(in) :: orbs_KS, orbs
  real(kind=8),dimension(3,nlr),intent(in) :: locregCenter
  type(locreg_descriptors),intent(in) :: glr_tmp
  type(local_zone_descriptors),intent(inout) :: lzd
  type(p2pComms),intent(inout) :: lbcomgp
  type(foe_data),intent(inout),optional :: lfoe
  type(comms_linear),intent(inout) :: lbcollcom
  type(comms_linear),intent(inout),optional :: lbcollcom_sr

  
  ! Local variables
  integer :: iorb, ilr, npsidim, istat
  character(len=*),parameter :: subname='update_locreg'

  call timing(iproc,'updatelocreg1','ON') 
  if (present(lfoe)) call nullify_foe(lfoe)
  !call nullify_comms_linear(lbcollcom)
  lbcollcom=comms_linear_null()
  if (present(lbcollcom_sr)) then
      !call nullify_comms_linear(lbcollcom_sr)
      lbcollcom_sr=comms_linear_null()
  end if
  call nullify_p2pComms(lbcomgp)
  call nullify_local_zone_descriptors(lzd)
  !!tag=1

  ! Allocate the array of localisation regions
  lzd%nlr=nlr
  allocate(lzd%Llr(lzd%nlr),stat=istat)
  do ilr=1,lzd%nlr
     lzd%Llr(ilr)=locreg_null()
  end do
  do ilr=1,lzd%nlr
      lzd%llr(ilr)%locrad=locrad(ilr)
      lzd%llr(ilr)%locrad_kernel=locrad_kernel(ilr)
      lzd%llr(ilr)%locrad_mult=locrad_mult(ilr)
      lzd%llr(ilr)%locregCenter=locregCenter(:,ilr)
  end do
  call timing(iproc,'updatelocreg1','OF') 
  call initLocregs(iproc, nproc, lzd, hx, hy, hz, astruct, orbs, glr_tmp, 's')!, llborbs)
  call timing(iproc,'updatelocreg1','ON') 
  call nullify_locreg_descriptors(lzd%glr)
  call copy_locreg_descriptors(glr_tmp, lzd%glr)
  lzd%hgrids(1)=hx
  lzd%hgrids(2)=hy
  lzd%hgrids(3)=hz

  npsidim = 0
  do iorb=1,orbs%norbp
   ilr=orbs%inwhichlocreg(iorb+orbs%isorb)
   npsidim = npsidim + lzd%llr(ilr)%wfd%nvctr_c+7*lzd%llr(ilr)%wfd%nvctr_f
  end do

  npsidim_orbs=max(npsidim,1)
  ! set npsidim_comp here too?!
  npsidim_comp=1

!  ! don't really want to keep this unless we do so right from the start, but for now keep it to avoid updating refs
!  orbs%eval=-.5d0

  call timing(iproc,'updatelocreg1','OF') 

  if (present(lfoe)) call init_foe(iproc, nproc, lzd, astruct, input, orbs_KS, orbs, lfoe, .false.)

  call init_comms_linear(iproc, nproc, npsidim_orbs, orbs, lzd, lbcollcom)
  if (present(lbcollcom_sr)) then
      call init_comms_linear_sumrho(iproc, nproc, lzd, orbs, nscatterarr, lbcollcom_sr)
  end if

  call initialize_communication_potential(iproc, nproc, nscatterarr, orbs, lzd, lbcomgp)
  call allocateCommunicationsBuffersPotential(lbcomgp, subname)

end subroutine update_locreg


subroutine update_ldiis_arrays(tmb, subname, ldiis)
  use module_base
  use module_types
  implicit none

  ! Calling arguments
  type(DFT_wavefunction),intent(in) :: tmb
  character(len=*),intent(in) :: subname
  type(localizedDIISParameters),intent(inout) :: ldiis

  ! Local variables
  integer :: iall, istat, ii, iorb, ilr

  iall=-product(shape(ldiis%phiHist))*kind(ldiis%phiHist)
  deallocate(ldiis%phiHist, stat=istat)
  call memocc(istat, iall, 'ldiis%phiHist', subname)
  iall=-product(shape(ldiis%hphiHist))*kind(ldiis%hphiHist)
  deallocate(ldiis%hphiHist, stat=istat)
  call memocc(istat, iall, 'ldiis%hphiHist', subname)

  ii=0
  do iorb=1,tmb%orbs%norbp
      ilr=tmb%orbs%inwhichlocreg(tmb%orbs%isorb+iorb)
      ii=ii+ldiis%isx*(tmb%lzd%llr(ilr)%wfd%nvctr_c+7*tmb%lzd%llr(ilr)%wfd%nvctr_f)
  end do

  allocate(ldiis%phiHist(ii), stat=istat)
  call memocc(istat, ldiis%phiHist, 'ldiis%phiHist', subname)
  allocate(ldiis%hphiHist(ii), stat=istat)
  call memocc(istat, ldiis%hphiHist, 'ldiis%hphiHist', subname)

end subroutine update_ldiis_arrays


subroutine allocate_auxiliary_basis_function(npsidim, subname, lphi, lhphi)
  use module_base
  implicit none

  ! Calling arguments
  integer,intent(in) :: npsidim
  real(kind=8),dimension(:),pointer,intent(out) :: lphi, lhphi
  character(len=*),intent(in) :: subname

  ! Local variables
  integer :: istat

  allocate(lphi(npsidim), stat=istat)
  call memocc(istat, lphi, 'lphi', subname)
  allocate(lhphi(npsidim), stat=istat)
  call memocc(istat, lhphi, 'lhphi', subname)

  call to_zero(npsidim, lphi(1))
  call to_zero(npsidim, lhphi(1))

end subroutine allocate_auxiliary_basis_function


subroutine deallocate_auxiliary_basis_function(subname, lphi, lhphi)
  use module_base
  implicit none

  ! Calling arguments
  real(kind=8),dimension(:),pointer :: lphi, lhphi
  character(len=*),intent(in) :: subname

  ! Local variables
  integer :: istat, iall

  iall=-product(shape(lphi))*kind(lphi)
  deallocate(lphi, stat=istat)
  call memocc(istat, iall, 'lphi', subname)
  iall=-product(shape(lhphi))*kind(lhphi)
  deallocate(lhphi, stat=istat)
  call memocc(istat, iall, 'lhphi', subname)

end subroutine deallocate_auxiliary_basis_function



subroutine destroy_new_locregs(iproc, nproc, tmb)
  use module_base
  use module_types
  use module_interfaces, except_this_one => destroy_new_locregs
  use communications_base, only: deallocate_comms_linear
  use communications, only: synchronize_onesided_communication
  implicit none

  ! Calling arguments
  integer,intent(in) :: iproc, nproc
  type(DFT_wavefunction),intent(inout) :: tmb

  ! Local variables
  character(len=*),parameter :: subname='destroy_new_locregs'

  !!call wait_p2p_communication(iproc, nproc, tmb%comgp)
  call synchronize_onesided_communication(iproc, nproc, tmb%comgp)
 ! call deallocateCommunicationsBuffersPotential(tmb%comgp, subname)
  call deallocate_p2pComms(tmb%comgp, subname)

  call deallocate_local_zone_descriptors(tmb%lzd, subname)
  call deallocate_orbitals_data(tmb%orbs, subname)
  call deallocate_foe(tmb%foe_obj, subname)

  call deallocate_comms_linear(tmb%collcom)
  call deallocate_comms_linear(tmb%collcom_sr)

end subroutine destroy_new_locregs


subroutine destroy_DFT_wavefunction(wfn)
  use module_base
  use module_types
  use module_interfaces, except_this_one => destroy_DFT_wavefunction
  use deallocatePointers
  use communications_base, only: deallocate_comms_linear
  use sparsematrix_base, only: deallocate_sparse_matrix, allocate_matrices, deallocate_matrices
  implicit none
  
  ! Calling arguments
  type(DFT_wavefunction),intent(inout) :: wfn

  ! Local variables
  integer :: istat, iall
  character(len=*),parameter :: subname='destroy_DFT_wavefunction'

  iall=-product(shape(wfn%psi))*kind(wfn%psi)
  deallocate(wfn%psi, stat=istat)
  call memocc(istat, iall, 'wfn%psi', subname)

  call deallocate_p2pComms(wfn%comgp, subname)
  call deallocate_foe(wfn%foe_obj, subname)
  !call deallocate_sparse_matrix(wfn%linmat%ovrlp, subname)
  !!call deallocate_sparse_matrix(wfn%linmat%ham, subname)
  !call deallocate_sparse_matrix(wfn%linmat%denskern_large, subname)
  !call deallocate_sparse_matrix(wfn%linmat%inv_ovrlp_large, subname)

  call deallocate_sparse_matrix(wfn%linmat%s, subname)
  call deallocate_sparse_matrix(wfn%linmat%m, subname)
  call deallocate_sparse_matrix(wfn%linmat%l, subname)
  call deallocate_sparse_matrix(wfn%linmat%ks, subname)
  call deallocate_sparse_matrix(wfn%linmat%ks_e, subname)
  call deallocate_matrices(wfn%linmat%ovrlp_)
  call deallocate_matrices(wfn%linmat%ham_)
  call deallocate_matrices(wfn%linmat%kernel_)

  call deallocate_orbitals_data(wfn%orbs, subname)
  !call deallocate_comms_cubic(wfn%comms, subname)
  call deallocate_comms_linear(wfn%collcom)
  call deallocate_comms_linear(wfn%collcom_sr)
  call deallocate_local_zone_descriptors(wfn%lzd, subname)

  if (associated(wfn%coeff)) then
      iall=-product(shape(wfn%coeff))*kind(wfn%coeff)
      deallocate(wfn%coeff, stat=istat)
      call memocc(istat, iall, 'wfn%coeff', subname)
  end if

  !call deallocate_p2pComms(wfn%ham_descr%comgp, subname)
  !call deallocate_local_zone_descriptors(wfn%ham_descr%lzd, subname)
  !call deallocate_matrixDescriptors_foe(wfn%ham_descr%mad, subname)
  !call deallocate_comms_linear(wfn%ham_descr%collcom, subname)

end subroutine destroy_DFT_wavefunction


subroutine update_wavefunctions_size(lzd,npsidim_orbs,npsidim_comp,orbs,iproc,nproc)
  use module_base
  use module_types
  implicit none

  ! Calling arguments
  type(local_zone_descriptors),intent(in) :: lzd
  type(orbitals_data),intent(in) :: orbs
  integer, intent(in) :: iproc, nproc
  integer, intent(out) :: npsidim_orbs, npsidim_comp

  ! Local variables
  integer :: npsidim, ilr, iorb
  integer :: nvctr_tot,jproc,istat,ierr,iall
  integer, allocatable, dimension(:) :: ncntt 
  integer, allocatable, dimension(:,:) :: nvctr_par
  character(len = *), parameter :: subname = "update_wavefunctions_size"

  npsidim = 0
  do iorb=1,orbs%norbp
   ilr=orbs%inwhichlocreg(iorb+orbs%isorb)
   npsidim = npsidim + lzd%Llr(ilr)%wfd%nvctr_c+7*lzd%Llr(ilr)%wfd%nvctr_f
  end do
  npsidim_orbs=max(npsidim,1)


  nvctr_tot = 1
  do iorb=1,orbs%norbp
     ilr=orbs%inwhichlocreg(iorb+orbs%isorb)
     nvctr_tot = max(nvctr_tot,lzd%llr(ilr)%wfd%nvctr_c+7*lzd%llr(ilr)%wfd%nvctr_f)
  end do
  if (nproc > 1) then
     call mpiallred(nvctr_tot, 1, mpi_max, bigdft_mpi%mpi_comm)
  end if

  allocate(nvctr_par(0:nproc-1,1),stat=istat)
  call memocc(istat,nvctr_par,'nvctr_par',subname)

  call kpts_to_procs_via_obj(nproc,1,nvctr_tot,nvctr_par)

  allocate(ncntt(0:nproc-1+ndebug),stat=istat)
  call memocc(istat,ncntt,'ncntt',subname)

  ncntt(:) = 0
  do jproc=0,nproc-1
     ncntt(jproc)=ncntt(jproc)+&
          nvctr_par(jproc,1)*orbs%norbp*orbs%nspinor
  end do

  npsidim_comp=sum(ncntt(0:nproc-1))

  iall=-product(shape(nvctr_par))*kind(nvctr_par)
  deallocate(nvctr_par,stat=istat)
  call memocc(istat,iall,'nvctr_par',subname) 

  iall=-product(shape(ncntt))*kind(ncntt)
  deallocate(ncntt,stat=istat)
  call memocc(istat,iall,'ncntt',subname)  

end subroutine update_wavefunctions_size


subroutine create_large_tmbs(iproc, nproc, KSwfn, tmb, denspot,nlpsp,input, at, rxyz, lowaccur_converged)
  use module_base
  use module_types
  use module_interfaces, except_this_one => create_large_tmbs
  implicit none

  ! Calling arguments
  integer,intent(in):: iproc, nproc
  type(DFT_Wavefunction),intent(inout):: KSwfn, tmb
  type(DFT_local_fields),intent(in):: denspot
  type(DFT_PSP_projectors), intent(inout) :: nlpsp
  type(input_variables),intent(in):: input
  type(atoms_data),intent(in):: at
  real(8),dimension(3,at%astruct%nat),intent(in):: rxyz
  logical,intent(in):: lowaccur_converged

  ! Local variables
  integer:: iorb, ilr, istat
  logical, dimension(:), allocatable :: lr_mask
  real(8),dimension(:,:),allocatable:: locrad_tmp
  real(8),dimension(:,:),allocatable:: locregCenter
  character(len=*),parameter:: subname='create_large_tmbs'

  call f_routine(id=subname)

  locregCenter=f_malloc((/3,tmb%lzd%nlr/),id='locregCenter')
  locrad_tmp=f_malloc((/tmb%lzd%nlr,3/),id='locrad_tmp')
  lr_mask=f_malloc0(tmb%lzd%nlr,id='lr_mask')

  do iorb=1,tmb%orbs%norb
      ilr=tmb%orbs%inwhichlocreg(iorb)
      locregCenter(:,ilr)=tmb%lzd%llr(ilr)%locregCenter
  end do
  do ilr=1,tmb%lzd%nlr
      locrad_tmp(ilr,1)=tmb%lzd%llr(ilr)%locrad+8.d0*tmb%lzd%hgrids(1)
      locrad_tmp(ilr,2)=tmb%lzd%llr(ilr)%locrad_kernel
      locrad_tmp(ilr,3)=tmb%lzd%llr(ilr)%locrad_mult
  end do

  !temporary,  moved from update_locreg
  tmb%orbs%eval=-0.5_gp
  call update_locreg(iproc, nproc, tmb%lzd%nlr, locrad_tmp(:,1), locrad_tmp(:,2), locrad_tmp(:,3), locregCenter, tmb%lzd%glr, &
       .false., denspot%dpbox%nscatterarr, tmb%lzd%hgrids(1), tmb%lzd%hgrids(2), tmb%lzd%hgrids(3), &
       at%astruct, input, KSwfn%orbs, tmb%orbs, tmb%ham_descr%lzd, tmb%ham_descr%npsidim_orbs, tmb%ham_descr%npsidim_comp, &
       tmb%ham_descr%comgp, tmb%ham_descr%collcom)

  call allocate_auxiliary_basis_function(max(tmb%ham_descr%npsidim_comp,tmb%ham_descr%npsidim_orbs), subname, &
       tmb%ham_descr%psi, tmb%hpsi)

  tmb%ham_descr%can_use_transposed=.false.
  nullify(tmb%ham_descr%psit_c)
  nullify(tmb%ham_descr%psit_f)
  allocate(tmb%confdatarr(tmb%orbs%norbp), stat=istat)

  if(.not.lowaccur_converged) then
      call define_confinement_data(tmb%confdatarr,tmb%orbs,rxyz,at,&
           tmb%ham_descr%lzd%hgrids(1),tmb%ham_descr%lzd%hgrids(2),tmb%ham_descr%lzd%hgrids(3),&
           4,input%lin%potentialPrefac_lowaccuracy,tmb%ham_descr%lzd,tmb%orbs%onwhichatom)
  else
      call define_confinement_data(tmb%confdatarr,tmb%orbs,rxyz,at,&
           tmb%ham_descr%lzd%hgrids(1),tmb%ham_descr%lzd%hgrids(2),tmb%ham_descr%lzd%hgrids(3),&
           4,input%lin%potentialPrefac_highaccuracy,tmb%ham_descr%lzd,tmb%orbs%onwhichatom)
  end if

  call f_free(locregCenter)
  call f_free(locrad_tmp)

  call update_lrmask_array(tmb%lzd%nlr,tmb%orbs,lr_mask)

  !when the new tmbs are created the projector descriptors can be updated
  call update_nlpsp(nlpsp,tmb%ham_descr%lzd%nlr,tmb%ham_descr%lzd%llr,KSwfn%Lzd%Glr,lr_mask)
  if (iproc == 0) call print_nlpsp(nlpsp)
  call f_free(lr_mask)
  call f_release_routine()
end subroutine create_large_tmbs

!>create the masking array to determine which localization regions have to be 
!! calculated
subroutine update_lrmask_array(nlr,orbs,lr_mask)
  use module_types, only: orbitals_data
  implicit none
  integer, intent(in) :: nlr
  type(orbitals_data), intent(in) :: orbs
  !> array of the masking, prior initialized to .false.
  logical, dimension(nlr), intent(inout) :: lr_mask
  !local variables
  integer :: ikpt,isorb,ieorb,nspinor,ilr,iorb

  !create the masking array according to the locregs which are known by the task
  if (orbs%norbp>0) then
      ikpt=orbs%iokpt(1)
      loop_kpt: do
         call orbs_in_kpt(ikpt,orbs,isorb,ieorb,nspinor)

         !activate all the localization regions which are present in the orbitals
         do iorb=isorb,ieorb
            ilr=orbs%inwhichlocreg(iorb+orbs%isorb)
            lr_mask(ilr)=.true.
         end do
         !last k-point has been treated
         if (ieorb == orbs%norbp) exit loop_kpt
         ikpt=ikpt+1
      end do loop_kpt
  end if
 
end subroutine update_lrmask_array

subroutine set_optimization_variables(input, at, lorbs, nlr, onwhichatom, confdatarr, &
     convCritMix, lowaccur_converged, nit_scc, mix_hist, alpha_mix, locrad, target_function, nit_basis, &
     convcrit_dmin, nitdmin, conv_crit_TMB)
  use module_base
  use module_types
  implicit none
  
  ! Calling arguments
  integer,intent(in) :: nlr
  type(orbitals_data),intent(in) :: lorbs
  type(input_variables),intent(in) :: input
  type(atoms_data),intent(in) :: at
  integer,dimension(lorbs%norb),intent(in) :: onwhichatom
  type(confpot_data),dimension(lorbs%norbp),intent(inout) :: confdatarr
  real(kind=8), intent(out) :: convCritMix, alpha_mix, convcrit_dmin, conv_crit_TMB
  logical, intent(in) :: lowaccur_converged
  integer, intent(out) :: nit_scc, mix_hist, nitdmin
  real(kind=8), dimension(nlr), intent(out) :: locrad
  integer, intent(out) :: target_function, nit_basis

  ! Local variables
  integer :: iorb, ilr, iiat

  if(lowaccur_converged) then
      do iorb=1,lorbs%norbp
          iiat=onwhichatom(lorbs%isorb+iorb)
          confdatarr(iorb)%prefac=input%lin%potentialPrefac_highaccuracy(at%astruct%iatype(iiat))
      end do
      target_function=TARGET_FUNCTION_IS_ENERGY
      nit_basis=input%lin%nItBasis_highaccuracy
      nit_scc=input%lin%nitSCCWhenFixed_highaccuracy
      mix_hist=input%lin%mixHist_highaccuracy
      do ilr=1,nlr
          locrad(ilr)=input%lin%locrad_highaccuracy(ilr)
      end do
      alpha_mix=input%lin%alpha_mix_highaccuracy
      convCritMix=input%lin%convCritMix_highaccuracy
      convcrit_dmin=input%lin%convCritDmin_highaccuracy
      nitdmin=input%lin%nItdmin_highaccuracy
      conv_crit_TMB=input%lin%convCrit_lowaccuracy
  else
      do iorb=1,lorbs%norbp
          iiat=onwhichatom(lorbs%isorb+iorb)
          confdatarr(iorb)%prefac=input%lin%potentialPrefac_lowaccuracy(at%astruct%iatype(iiat))
      end do
      target_function=TARGET_FUNCTION_IS_TRACE
      nit_basis=input%lin%nItBasis_lowaccuracy
      nit_scc=input%lin%nitSCCWhenFixed_lowaccuracy
      mix_hist=input%lin%mixHist_lowaccuracy
      do ilr=1,nlr
          locrad(ilr)=input%lin%locrad_lowaccuracy(ilr)
      end do
      alpha_mix=input%lin%alpha_mix_lowaccuracy
      convCritMix=input%lin%convCritMix_lowaccuracy
      convcrit_dmin=input%lin%convCritDmin_lowaccuracy
      nitdmin=input%lin%nItdmin_lowaccuracy
      conv_crit_TMB=input%lin%convCrit_highaccuracy
  end if

  !!! new hybrid version... not the best place here
  !!if (input%lin%nit_highaccuracy==-1) then
  !!    do iorb=1,lorbs%norbp
  !!        ilr=lorbs%inwhichlocreg(lorbs%isorb+iorb)
  !!        iiat=onwhichatom(lorbs%isorb+iorb)
  !!        confdatarr(iorb)%prefac=input%lin%potentialPrefac_lowaccuracy(at%astruct%iatype(iiat))
  !!    end do
  !!    wfnmd%bs%target_function=TARGET_FUNCTION_IS_HYBRID
  !!    wfnmd%bs%nit_basis_optimization=input%lin%nItBasis_lowaccuracy
  !!    wfnmd%bs%conv_crit=input%lin%convCrit_lowaccuracy
  !!    nit_scc=input%lin%nitSCCWhenFixed_lowaccuracy
  !!    mix_hist=input%lin%mixHist_lowaccuracy
  !!    do ilr=1,nlr
  !!        locrad(ilr)=input%lin%locrad_lowaccuracy(ilr)
  !!    end do
  !!    alpha_mix=input%lin%alpha_mix_lowaccuracy
  !!end if

end subroutine set_optimization_variables



subroutine adjust_locregs_and_confinement(iproc, nproc, hx, hy, hz, at, input, &
           rxyz, KSwfn, tmb, denspot, nlpsp,ldiis, locreg_increased, lowaccur_converged, locrad)
  use module_base
  use module_types
  use module_interfaces, except_this_one => adjust_locregs_and_confinement
  use yaml_output
  use communications_base, only: deallocate_comms_linear
  use communications, only: synchronize_onesided_communication
  use sparsematrix_base, only: sparse_matrix_null, deallocate_sparse_matrix, allocate_matrices, deallocate_matrices
  use sparsematrix_init, only: init_sparse_matrix, check_kernel_cutoff!, init_sparsity_from_distance
  implicit none
  
  ! Calling argument
  integer,intent(in) :: iproc, nproc
  real(8),intent(in) :: hx, hy, hz
  type(atoms_data),intent(in) :: at
  type(input_variables),intent(in) :: input
  real(8),dimension(3,at%astruct%nat),intent(in):: rxyz
  type(DFT_wavefunction),intent(inout) :: KSwfn, tmb
  type(DFT_local_fields),intent(inout) :: denspot
  type(DFT_PSP_projectors), intent(inout) :: nlpsp
  type(localizedDIISParameters),intent(inout) :: ldiis
  logical, intent(out) :: locreg_increased
  logical, intent(in) :: lowaccur_converged
  real(8), dimension(tmb%lzd%nlr), intent(inout) :: locrad

  ! Local variables
  integer :: iall, istat, ilr, npsidim_orbs_tmp, npsidim_comp_tmp
  real(kind=8),dimension(:,:),allocatable :: locregCenter
  real(kind=8),dimension(:),allocatable :: lphilarge, locrad_kernel, locrad_mult
  type(local_zone_descriptors) :: lzd_tmp
  character(len=*), parameter :: subname='adjust_locregs_and_confinement'

  locreg_increased=.false.
  if(lowaccur_converged ) then
      do ilr = 1, tmb%lzd%nlr
         if(input%lin%locrad_highaccuracy(ilr) /= input%lin%locrad_lowaccuracy(ilr)) then
             !!if(iproc==0) write(*,'(1x,a)') 'Increasing the localization radius for the high accuracy part.'
             if (iproc==0) call yaml_map('Increasing the localization radius for the high accuracy part',.true.)
             locreg_increased=.true.
             exit
         end if
      end do
  end if
  if (iproc==0) then
      if (locreg_increased) then
          call yaml_map('Locreg increased',.true.)
      else
          call yaml_map('Locreg increased',.false.)
      end if
  end if

  if(locreg_increased) then
     !tag=1
     !call wait_p2p_communication(iproc, nproc, tmb%comgp)
     call synchronize_onesided_communication(iproc, nproc, tmb%comgp)
     call deallocate_p2pComms(tmb%comgp, subname)

     call deallocate_comms_linear(tmb%collcom)
     call deallocate_comms_linear(tmb%collcom_sr)

     call nullify_local_zone_descriptors(lzd_tmp)
     call copy_local_zone_descriptors(tmb%lzd, lzd_tmp, subname)
     call deallocate_local_zone_descriptors(tmb%lzd, subname)

     npsidim_orbs_tmp = tmb%npsidim_orbs
     npsidim_comp_tmp = tmb%npsidim_comp

     call deallocate_foe(tmb%foe_obj, subname)

     !call deallocate_sparse_matrix(tmb%linmat%denskern_large, subname)
     !call deallocate_sparse_matrix(tmb%linmat%inv_ovrlp_large, subname)
     !call deallocate_sparse_matrix(tmb%linmat%ovrlp, subname)
     !!call deallocate_sparse_matrix(tmb%linmat%ham, subname)

     call deallocate_sparse_matrix(tmb%linmat%s, subname)
     call deallocate_sparse_matrix(tmb%linmat%m, subname)
     call deallocate_sparse_matrix(tmb%linmat%l, subname)
     call deallocate_sparse_matrix(tmb%linmat%ks, subname)
     call deallocate_sparse_matrix(tmb%linmat%ks_e, subname)
     call deallocate_matrices(tmb%linmat%ovrlp_)
     call deallocate_matrices(tmb%linmat%ham_)
     call deallocate_matrices(tmb%linmat%kernel_)

     allocate(locregCenter(3,lzd_tmp%nlr), stat=istat)
     call memocc(istat, locregCenter, 'locregCenter', subname)
     allocate(locrad_kernel(lzd_tmp%nlr),stat=istat)
     call memocc(istat,locrad_kernel,'locrad_kernel',subname)
     allocate(locrad_mult(lzd_tmp%nlr),stat=istat)
     call memocc(istat,locrad_mult,'locrad_mult',subname)
     do ilr=1,lzd_tmp%nlr
        locregCenter(:,ilr)=lzd_tmp%llr(ilr)%locregCenter
        locrad_kernel(ilr)=lzd_tmp%llr(ilr)%locrad_kernel
        locrad_mult(ilr)=lzd_tmp%llr(ilr)%locrad_mult
     end do

     !temporary,  moved from update_locreg
     tmb%orbs%eval=-0.5_gp
     call update_locreg(iproc, nproc, lzd_tmp%nlr, locrad, locrad_kernel, locrad_mult, locregCenter, lzd_tmp%glr, .false., &
          denspot%dpbox%nscatterarr, hx, hy, hz, at%astruct, input, KSwfn%orbs, tmb%orbs, tmb%lzd, &
          tmb%npsidim_orbs, tmb%npsidim_comp, tmb%comgp, tmb%collcom, tmb%foe_obj, tmb%collcom_sr)

     iall=-product(shape(locregCenter))*kind(locregCenter)
     deallocate(locregCenter, stat=istat)
     call memocc(istat, iall, 'locregCenter', subname)

     iall=-product(shape(locrad_kernel))*kind(locrad_kernel)
     deallocate(locrad_kernel, stat=istat)
     call memocc(istat, iall, 'locrad_kernel', subname)

     iall=-product(shape(locrad_mult))*kind(locrad_mult)
     deallocate(locrad_mult, stat=istat)
     call memocc(istat, iall, 'locrad_mult', subname)

     ! calculate psi in new locreg
     allocate(lphilarge(tmb%npsidim_orbs), stat=istat)
     call memocc(istat, lphilarge, 'lphilarge', subname)
     call to_zero(tmb%npsidim_orbs, lphilarge(1))
     call small_to_large_locreg(iproc, npsidim_orbs_tmp, tmb%npsidim_orbs, lzd_tmp, tmb%lzd, &
          tmb%orbs, tmb%psi, lphilarge)

     call deallocate_local_zone_descriptors(lzd_tmp, subname)
     iall=-product(shape(tmb%psi))*kind(tmb%psi)
     deallocate(tmb%psi, stat=istat)
     call memocc(istat, iall, 'tmb%psi', subname)
     allocate(tmb%psi(tmb%npsidim_orbs), stat=istat)
     call memocc(istat, tmb%psi, 'tmb%psi', subname)
     call vcopy(tmb%npsidim_orbs, lphilarge(1), 1, tmb%psi(1), 1)
     iall=-product(shape(lphilarge))*kind(lphilarge)
     deallocate(lphilarge, stat=istat)
     call memocc(istat, iall, 'lphilarge', subname) 
     
     call update_ldiis_arrays(tmb, subname, ldiis)

     ! Emit that lzd has been changed.
     if (tmb%c_obj /= 0) then
        call kswfn_emit_lzd(tmb, iproc, nproc)
     end if

     ! Now update hamiltonian descriptors
     !call destroy_new_locregs(iproc, nproc, tmblarge)

     ! to eventually be better sorted - replace with e.g. destroy_hamiltonian_descriptors
     call synchronize_onesided_communication(iproc, nproc, tmb%ham_descr%comgp)
     call deallocate_p2pComms(tmb%ham_descr%comgp, subname)
     call deallocate_local_zone_descriptors(tmb%ham_descr%lzd, subname)
     call deallocate_comms_linear(tmb%ham_descr%collcom)

     call deallocate_auxiliary_basis_function(subname, tmb%ham_descr%psi, tmb%hpsi)
     if(tmb%ham_descr%can_use_transposed) then
        iall=-product(shape(tmb%ham_descr%psit_c))*kind(tmb%ham_descr%psit_c)
        deallocate(tmb%ham_descr%psit_c, stat=istat)
        call memocc(istat, iall, 'tmb%ham_descr%psit_c', subname)
        iall=-product(shape(tmb%ham_descr%psit_f))*kind(tmb%ham_descr%psit_f)
        deallocate(tmb%ham_descr%psit_f, stat=istat)
        call memocc(istat, iall, 'tmb%ham_descr%psit_f', subname)
        tmb%ham_descr%can_use_transposed=.false.
     end if
     
     deallocate(tmb%confdatarr, stat=istat)

     call create_large_tmbs(iproc, nproc, KSwfn, tmb, denspot,nlpsp, input, at, rxyz, lowaccur_converged)

     ! check the extent of the kernel cutoff (must be at least shamop radius)
     call check_kernel_cutoff(iproc, tmb%orbs, at, tmb%lzd)

     ! Update sparse matrices
     !!call init_sparse_matrix_wrapper(iproc, nproc, tmb%orbs, tmb%ham_descr%lzd, at%astruct, &
     !!     input%store_index, imode=1, smat=tmb%linmat%ham)
     call init_sparse_matrix_wrapper(iproc, nproc, tmb%orbs, tmb%ham_descr%lzd, at%astruct, &
          input%store_index, imode=1, smat=tmb%linmat%m)
     call allocate_matrices(tmb%linmat%m, allocate_full=.false., &
          matname='tmb%linmat%ham_', mat=tmb%linmat%ham_)
     !!call init_matrixindex_in_compressed_fortransposed(iproc, nproc, tmb%orbs, &
     !!     tmb%collcom, tmb%ham_descr%collcom, tmb%collcom_sr, tmb%linmat%ham)
     call init_matrixindex_in_compressed_fortransposed(iproc, nproc, tmb%orbs, &
          tmb%collcom, tmb%ham_descr%collcom, tmb%collcom_sr, tmb%linmat%m)

     !call init_sparse_matrix_wrapper(iproc, nproc, tmb%orbs, tmb%lzd, at%astruct, &
     !     input%store_index, imode=1, smat=tmb%linmat%ovrlp)
     call init_sparse_matrix_wrapper(iproc, nproc, tmb%orbs, tmb%lzd, at%astruct, &
          input%store_index, imode=1, smat=tmb%linmat%s)
     call allocate_matrices(tmb%linmat%s, allocate_full=.false., &
          matname='tmb%linmat%ovrlp_', mat=tmb%linmat%ovrlp_)
     !call init_matrixindex_in_compressed_fortransposed(iproc, nproc, tmb%orbs, &
     !     tmb%collcom, tmb%ham_descr%collcom, tmb%collcom_sr, tmb%linmat%ovrlp)
     call init_matrixindex_in_compressed_fortransposed(iproc, nproc, tmb%orbs, &
          tmb%collcom, tmb%ham_descr%collcom, tmb%collcom_sr, tmb%linmat%s)

     call check_kernel_cutoff(iproc, tmb%orbs, at, tmb%lzd)
     !!call init_sparse_matrix_wrapper(iproc, nproc, tmb%orbs, tmb%lzd, at%astruct, &
     !!     input%store_index, imode=2, smat=tmb%linmat%denskern_large)
     call init_sparse_matrix_wrapper(iproc, nproc, tmb%orbs, tmb%lzd, at%astruct, &
          input%store_index, imode=2, smat=tmb%linmat%l)
     call allocate_matrices(tmb%linmat%l, allocate_full=.false., &
          matname='tmb%linmat%kernel_', mat=tmb%linmat%kernel_)
     !!call init_matrixindex_in_compressed_fortransposed(iproc, nproc, tmb%orbs, &
     !!     tmb%collcom, tmb%ham_descr%collcom, tmb%collcom_sr, tmb%linmat%denskern_large)
     call init_matrixindex_in_compressed_fortransposed(iproc, nproc, tmb%orbs, &
          tmb%collcom, tmb%ham_descr%collcom, tmb%collcom_sr, tmb%linmat%l)

     !tmb%linmat%inv_ovrlp_large=sparse_matrix_null()
     !call sparse_copy_pattern(tmb%linmat%l, tmb%linmat%inv_ovrlp_large, iproc, subname)


     tmb%linmat%ks = sparse_matrix_null()
     tmb%linmat%ks_e = sparse_matrix_null()
     if (input%lin%scf_mode/=LINEAR_FOE .or. input%lin%pulay_correction .or.  input%lin%new_pulay_correction .or. &
         (input%lin%plotBasisFunctions /= WF_FORMAT_NONE) .or. input%lin%diag_end) then
         call init_sparse_matrix_for_KSorbs(iproc, nproc, KSwfn%orbs, input, input%lin%extra_states, &
              tmb%linmat%ks, tmb%linmat%ks_e)
     end if



  else ! no change in locrad, just confining potential that needs updating

     call define_confinement_data(tmb%confdatarr,tmb%orbs,rxyz,at,&
          tmb%ham_descr%lzd%hgrids(1),tmb%ham_descr%lzd%hgrids(2),tmb%ham_descr%lzd%hgrids(3),&
          4,input%lin%potentialPrefac_highaccuracy,tmb%ham_descr%lzd,tmb%orbs%onwhichatom)

  end if

end subroutine adjust_locregs_and_confinement



subroutine adjust_DIIS_for_high_accuracy(input, denspot, mixdiis, lowaccur_converged, ldiis_coeff_hist, ldiis_coeff_changed)
  use module_base
  use module_types
  use module_interfaces, except_this_one => adjust_DIIS_for_high_accuracy
  implicit none
  
  ! Calling arguments
  type(input_variables),intent(in) :: input
  type(DFT_local_fields),intent(inout) :: denspot
  type(mixrhopotDIISParameters),intent(inout) :: mixdiis
  logical, intent(in) :: lowaccur_converged
  integer, intent(inout) :: ldiis_coeff_hist
  logical, intent(out) :: ldiis_coeff_changed  

  if(lowaccur_converged) then
     if(input%lin%mixHist_lowaccuracy==0 .and. input%lin%mixHist_highaccuracy>0) then
        call initializeMixrhopotDIIS(input%lin%mixHist_highaccuracy, denspot%dpbox%ndimrhopot, mixdiis)
     else if(input%lin%mixHist_lowaccuracy>0 .and. input%lin%mixHist_highaccuracy==0) then
        call deallocateMixrhopotDIIS(mixdiis)
     end if
     if (input%lin%scf_mode==LINEAR_DIRECT_MINIMIZATION) then
        ! check whether ldiis_coeff_hist arrays will need reallocating due to change in history length
        if (ldiis_coeff_hist /= input%lin%dmin_hist_highaccuracy) then
           ldiis_coeff_changed=.true.
        else
           ldiis_coeff_changed=.false.
        end if
        ldiis_coeff_hist=input%lin%dmin_hist_highaccuracy
     end if
  else
     if (input%lin%scf_mode==LINEAR_DIRECT_MINIMIZATION) then
        ldiis_coeff_changed=.false.
     end if
  end if
  
end subroutine adjust_DIIS_for_high_accuracy


subroutine check_whether_lowaccuracy_converged(itout, nit_lowaccuracy, lowaccuracy_convcrit, &
     lowaccur_converged, pnrm_out)
  use module_base
  use module_types
  implicit none

  ! Calling arguments
  integer,intent(in) :: itout
  integer,intent(in) :: nit_lowaccuracy
  real(8),intent(in) :: lowaccuracy_convcrit
  logical, intent(inout) :: lowaccur_converged
  real(kind=8), intent(in) :: pnrm_out
  
  if(.not.lowaccur_converged .and. &
       (itout>=nit_lowaccuracy+1 .or. pnrm_out<lowaccuracy_convcrit)) then
     lowaccur_converged=.true.
     !cur_it_highaccuracy=0
  end if 

end subroutine check_whether_lowaccuracy_converged



subroutine set_variables_for_hybrid(nlr, input, at, orbs, lowaccur_converged, confdatarr, &
           target_function, nit_basis, nit_scc, mix_hist, locrad, alpha_mix, convCritMix, &
           conv_crit_TMB)
  use module_base
  use module_types
  implicit none

  ! Calling arguments
  integer,intent(in) :: nlr
  type(input_variables),intent(in) :: input
  type(atoms_data),intent(in) :: at
  type(orbitals_data),intent(in) :: orbs
  logical,intent(out) :: lowaccur_converged
  type(confpot_data),dimension(orbs%norbp),intent(inout) :: confdatarr
  integer,intent(out) :: target_function, nit_basis, nit_scc, mix_hist
  real(kind=8),dimension(nlr),intent(out) :: locrad
  real(kind=8),intent(out) :: alpha_mix, convCritMix, conv_crit_TMB

  ! Local variables
  integer :: iorb, ilr, iiat

  lowaccur_converged=.false.
  do iorb=1,orbs%norbp
      ilr=orbs%inwhichlocreg(orbs%isorb+iorb)
      iiat=orbs%onwhichatom(orbs%isorb+iorb)
      confdatarr(iorb)%prefac=input%lin%potentialPrefac_lowaccuracy(at%astruct%iatype(iiat))
  end do
  target_function=TARGET_FUNCTION_IS_HYBRID
  nit_basis=input%lin%nItBasis_lowaccuracy
  nit_scc=input%lin%nitSCCWhenFixed_lowaccuracy
  mix_hist=input%lin%mixHist_lowaccuracy
  do ilr=1,nlr
      locrad(ilr)=input%lin%locrad_lowaccuracy(ilr)
  end do
  alpha_mix=input%lin%alpha_mix_lowaccuracy
  convCritMix=input%lin%convCritMix_lowaccuracy
  conv_crit_TMB=input%lin%convCrit_lowaccuracy

end subroutine set_variables_for_hybrid




subroutine increase_FOE_cutoff(iproc, nproc, lzd, astruct, input, orbs_KS, orbs, foe_obj, init)
  use module_base
  use module_types
  use module_interfaces, except_this_one => increase_FOE_cutoff
  use yaml_output
  implicit none

  ! Calling arguments
  integer,intent(in) :: iproc, nproc
  type(local_zone_descriptors),intent(in) :: lzd
  type(atomic_structure),intent(in) :: astruct
  type(input_variables),intent(in) :: input
  type(orbitals_data),intent(in) :: orbs_KS, orbs
  type(foe_data),intent(out) :: foe_obj
  logical,intent(in) :: init
  ! Local variables
  real(kind=8),save :: cutoff_incr
  character(len=*),parameter :: subname='increase_FOE_cutoff'

  ! Just initialize the save variable
  if (init) then
      cutoff_incr=0.d0
      return
  end if

  ! Deallocate the pointers
  call deallocate_foe(foe_obj, subname)

  ! How much should the cutoff be increased
  cutoff_incr=cutoff_incr+1.d0

  if (iproc==0) then
      call yaml_newline()
      call yaml_map('Need to re-initialize FOE cutoff',.true.)
      call yaml_newline()
      call yaml_map('Total increase of FOE cutoff wrt input values',cutoff_incr,fmt='(f5.1)')
  end if

  ! Re-initialize the foe data
  call init_foe(iproc, nproc, lzd, astruct, input, orbs_KS, orbs, foe_obj, reset=.false., &
       cutoff_incr=cutoff_incr)

end subroutine increase_FOE_cutoff


!> Set negative entries to zero
subroutine clean_rho(iproc, nproc, npt, rho)
  use module_base
  use yaml_output
  implicit none

  ! Calling arguments
  integer,intent(in) :: iproc, nproc, npt
  real(kind=8),dimension(npt),intent(inout) :: rho

  ! Local variables
  integer :: ncorrection, ipt, ierr
  real(kind=8) :: charge_correction

  if (iproc==0) then
      call yaml_newline()
      call yaml_map('Need to correct charge density',.true.)
      call yaml_warning('set to 1.d-20 instead of 0.d0')
  end if

  ncorrection=0
  charge_correction=0.d0
  do ipt=1,npt
      if (rho(ipt)<0.d0) then
          if (rho(ipt)>=-1.d-9) then
              ! negative, but small, so simply set to zero
              charge_correction=charge_correction+rho(ipt)
              !rho(ipt)=0.d0
              rho(ipt)=1.d-20
              ncorrection=ncorrection+1
          else
              ! negative, but non-negligible, so issue a warning
              call yaml_warning('considerable negative rho, value: '//trim(yaml_toa(rho(ipt),fmt='(es12.4)'))) 
              charge_correction=charge_correction+rho(ipt)
              !rho(ipt)=0.d0
              rho(ipt)=1.d-20
              ncorrection=ncorrection+1
          end if
      end if
  end do

  if (nproc > 1) then
      call mpiallred(ncorrection, 1, mpi_sum, bigdft_mpi%mpi_comm)
      call mpiallred(charge_correction, 1, mpi_sum, bigdft_mpi%mpi_comm)
  end if

  if (iproc==0) then
      call yaml_newline()
      call yaml_map('number of corrected points',ncorrection)
      call yaml_newline()
      call yaml_map('total charge correction',abs(charge_correction),fmt='(es14.5)')
      call yaml_newline()
  end if
  
end subroutine clean_rho



subroutine corrections_for_negative_charge(iproc, nproc, KSwfn, at, input, tmb, denspot)
  use module_types
  use module_interfaces
  use yaml_output
  implicit none

  ! Calling arguments
  integer,intent(in) :: iproc, nproc
  type(DFT_wavefunction),intent(in) :: KSwfn
  type(atoms_data),intent(in) :: at
  type(input_variables),intent(in) :: input
  type(DFT_wavefunction),intent(inout) :: tmb
  type(DFT_local_fields), intent(inout) :: denspot

  !!if (iproc==0) then
  !!    !call yaml_open_sequence()
  !!    !call yaml_open_map()
  !!    call yaml_newline()
  !!    call yaml_warning('Charge density contains negative points, need to increase FOE cutoff')
  !!end if
  !!call increase_FOE_cutoff(iproc, nproc, tmb%lzd, at%astruct, input, KSwfn%orbs, tmb%orbs, tmb%foe_obj, init=.false.)
  if (iproc==0) call yaml_warning('No increase of FOE cutoff')
  call clean_rho(iproc, nproc, KSwfn%Lzd%Glr%d%n1i*KSwfn%Lzd%Glr%d%n2i*denspot%dpbox%n3d, denspot%rhov)
  if (iproc==0) then
      !call yaml_close_map()
      !call yaml_close_sequence()
  end if

end subroutine corrections_for_negative_charge




subroutine determine_sparsity_pattern(iproc, nproc, orbs, lzd, nnonzero, nonzero)
      use module_base
      use module_types
      use module_interfaces, except_this_one => determine_sparsity_pattern
      implicit none
    
      ! Calling arguments
      integer,intent(in) :: iproc, nproc
      type(orbitals_data),intent(in) :: orbs
      type(local_zone_descriptors),intent(in) :: lzd
      integer,intent(out) :: nnonzero
      integer,dimension(:),pointer,intent(out) :: nonzero
    
      ! Local variables
      integer :: jproc, iorb, jorb, ioverlapMPI, ioverlaporb, ilr, jlr, ilrold
      integer :: iiorb, istat, iall, noverlaps, ierr, ii
      logical :: isoverlap
      integer :: onseg
      logical,dimension(:,:),allocatable :: overlapMatrix
      integer,dimension(:),allocatable :: noverlapsarr, displs, recvcnts, op_noverlaps
      integer,dimension(:,:),allocatable :: overlaps_op, op_overlaps
      integer,dimension(:,:,:),allocatable :: overlaps_nseg
      !character(len=*),parameter :: subname='determine_overlap_from_descriptors'

      call f_routine('determine_sparsity_pattern')
    
      overlapMatrix = f_malloc((/orbs%norb,maxval(orbs%norb_par(:,0))/),id='overlapMatrix')
      noverlapsarr = f_malloc(orbs%norbp,id='noverlapsarr')
      !!allocate(overlapMatrix(orbs%norb,maxval(orbs%norb_par(:,0))), stat=istat)
      !!call memocc(istat, overlapMatrix, 'overlapMatrix', subname)
      !!allocate(noverlapsarr(orbs%norbp), stat=istat)
      !!call memocc(istat, noverlapsarr, 'noverlapsarr', subname)
    
      overlapMatrix=.false.
      do iorb=1,orbs%norbp
         ioverlaporb=0 ! counts the overlaps for the given orbital.
         iiorb=orbs%isorb+iorb
         ilr=orbs%inWhichLocreg(iiorb)
         do jorb=1,orbs%norb
            jlr=orbs%inWhichLocreg(jorb)
            call check_overlap_cubic_periodic(lzd%Glr,lzd%llr(ilr),lzd%llr(jlr),isoverlap)
            if(isoverlap) then
               ! From the viewpoint of the box boundaries, an overlap between ilr and jlr is possible.
               ! Now explicitely check whether there is an overlap by using the descriptors.
               call check_overlap_from_descriptors_periodic(lzd%llr(ilr)%wfd%nseg_c, lzd%llr(jlr)%wfd%nseg_c,&
                    lzd%llr(ilr)%wfd%keyglob, lzd%llr(jlr)%wfd%keyglob, &
                    isoverlap, onseg)
               if(isoverlap) then
                  ! There is really an overlap
                  overlapMatrix(jorb,iorb)=.true.
                  ioverlaporb=ioverlaporb+1
               else
                  overlapMatrix(jorb,iorb)=.false.
               end if
            else
               overlapMatrix(jorb,iorb)=.false.
            end if
         end do
         noverlapsarr(iorb)=ioverlaporb
      end do


      overlaps_op = f_malloc((/maxval(noverlapsarr),orbs%norbp/),id='overlaps_op')
      !allocate(overlaps_op(maxval(noverlapsarr),orbs%norbp), stat=istat)
      !call memocc(istat, overlaps_op, 'overlaps_op', subname)
    
      ! Now we know how many overlaps have to be calculated, so determine which orbital overlaps
      ! with which one. This is essentially the same loop as above, but we use the array 'overlapMatrix'
      ! which indicates the overlaps.
      iiorb=0
      ilrold=-1
      do iorb=1,orbs%norbp
         ioverlaporb=0 ! counts the overlaps for the given orbital.
         iiorb=orbs%isorb+iorb
         do jorb=1,orbs%norb
            if(overlapMatrix(jorb,iorb)) then
               ioverlaporb=ioverlaporb+1
               overlaps_op(ioverlaporb,iorb)=jorb
            end if
         end do 
      end do


      nnonzero=0
      do iorb=1,orbs%norbp
          nnonzero=nnonzero+noverlapsarr(iorb)
      end do
      nonzero = f_malloc_ptr(nnonzero,id='nonzero')
      ii=0
      do iorb=1,orbs%norbp
          iiorb=orbs%isorb+iorb
          do jorb=1,noverlapsarr(iorb)
              ii=ii+1
              nonzero(ii)=(iiorb-1)*orbs%norb+overlaps_op(jorb,iorb)
          end do
      end do

      call f_free(overlapMatrix)
      call f_free(noverlapsarr)
      call f_free(overlaps_op)
    
    
    !!  iall=-product(shape(overlapMatrix))*kind(overlapMatrix)
    !!  deallocate(overlapMatrix, stat=istat)
    !!  call memocc(istat, iall, 'overlapMatrix', subname)
    !!
    !!  iall=-product(shape(noverlapsarr))*kind(noverlapsarr)
    !!  deallocate(noverlapsarr, stat=istat)
    !!  call memocc(istat, iall, 'noverlapsarr', subname)
    !!
    !!  iall=-product(shape(overlaps_op))*kind(overlaps_op)
    !!  deallocate(overlaps_op, stat=istat)
    !!  call memocc(istat, iall, 'overlaps_op', subname)

end subroutine determine_sparsity_pattern



subroutine determine_sparsity_pattern_distance(orbs, lzd, astruct, cutoff, nnonzero, nonzero)
  use module_base
  use module_types
  implicit none

  ! Calling arguments
  type(orbitals_data),intent(in) :: orbs
  type(local_zone_descriptors),intent(in) :: lzd
  type(atomic_structure),intent(in) :: astruct
  real(kind=8),dimension(lzd%nlr),intent(in) :: cutoff
  integer,intent(out) :: nnonzero
  integer,dimension(:),pointer,intent(out) :: nonzero

  ! Local variables
  integer :: iorb, iiorb, ilr, iwa, itype, jjorb, jlr, jwa, jtype, ii
  real(kind=8) :: tt, cut

  call f_routine('determine_sparsity_pattern_distance')

      nnonzero=0
      do iorb=1,orbs%norbp
         iiorb=orbs%isorb+iorb
         ilr=orbs%inwhichlocreg(iiorb)
         iwa=orbs%onwhichatom(iiorb)
         itype=astruct%iatype(iwa)
         do jjorb=1,orbs%norb
            jlr=orbs%inwhichlocreg(jjorb)
            jwa=orbs%onwhichatom(jjorb)
            jtype=astruct%iatype(jwa)
            tt = (lzd%llr(ilr)%locregcenter(1)-lzd%llr(jlr)%locregcenter(1))**2 + &
                 (lzd%llr(ilr)%locregcenter(2)-lzd%llr(jlr)%locregcenter(2))**2 + &
                 (lzd%llr(ilr)%locregcenter(3)-lzd%llr(jlr)%locregcenter(3))**2
            cut = cutoff(ilr)+cutoff(jlr)!+2.d0*incr
            tt=sqrt(tt)
            if (tt<=cut) then
               nnonzero=nnonzero+1
            end if
         end do
      end do
      !call mpiallred(nnonzero, 1, mpi_sum, bigdft_mpi%mpi_comm, ierr)
      nonzero = f_malloc_ptr(nnonzero,id='nonzero')

      ii=0
      do iorb=1,orbs%norbp
         iiorb=orbs%isorb+iorb
         ilr=orbs%inwhichlocreg(iiorb)
         iwa=orbs%onwhichatom(iiorb)
         itype=astruct%iatype(iwa)
         do jjorb=1,orbs%norb
            jlr=orbs%inwhichlocreg(jjorb)
            jwa=orbs%onwhichatom(jjorb)
            jtype=astruct%iatype(jwa)
            tt = (lzd%llr(ilr)%locregcenter(1)-lzd%llr(jlr)%locregcenter(1))**2 + &
                 (lzd%llr(ilr)%locregcenter(2)-lzd%llr(jlr)%locregcenter(2))**2 + &
                 (lzd%llr(ilr)%locregcenter(3)-lzd%llr(jlr)%locregcenter(3))**2
            cut = cutoff(ilr)+cutoff(jlr)!+2.d0*incr
            tt=sqrt(tt)
            if (tt<=cut) then
               ii=ii+1
               nonzero(ii)=(iiorb-1)*orbs%norb+jjorb
            end if
         end do
      end do

  call f_release_routine()

end subroutine determine_sparsity_pattern_distance


subroutine init_sparse_matrix_wrapper(iproc, nproc, orbs, lzd, astruct, store_index, imode, smat)
  use module_base
  use module_types
  use sparsematrix_init, only: init_sparse_matrix
  use module_interfaces, except_this_one => init_sparse_matrix_wrapper
  implicit none

  ! Calling arguments
  integer,intent(in) :: iproc, nproc, imode
  type(orbitals_data),intent(in) :: orbs
  type(local_zone_descriptors),intent(in) :: lzd
  type(atomic_structure),intent(in) :: astruct
  logical,intent(in) :: store_index
  type(sparse_matrix), intent(out) :: smat
  
  ! Local variables
  integer :: nnonzero, nnonzero_mult, ilr
  integer,dimension(:),pointer :: nonzero, nonzero_mult
  real(kind=8),dimension(:),allocatable :: cutoff
  integer,parameter :: KEYS=1
  integer,parameter :: DISTANCE=2

  cutoff = f_malloc(lzd%nlr,id='cutoff')

  do ilr=1,lzd%nlr
      cutoff(ilr)=lzd%llr(ilr)%locrad_mult
  end do

  if (imode==KEYS) then
      call determine_sparsity_pattern(iproc, nproc, orbs, lzd, nnonzero, nonzero)
  else if (imode==DISTANCE) then
      call determine_sparsity_pattern_distance(orbs, lzd, astruct, lzd%llr(:)%locrad_kernel, nnonzero, nonzero)
  else
      stop 'wrong imode'
  end if
  call determine_sparsity_pattern_distance(orbs, lzd, astruct, lzd%llr(:)%locrad_mult, nnonzero_mult, nonzero_mult)
  call init_sparse_matrix(iproc, nproc, orbs%norb, orbs%norbp, orbs%isorb, store_index, &
       nnonzero, nonzero, nnonzero_mult, nonzero_mult, smat)
  call f_free_ptr(nonzero)
  call f_free_ptr(nonzero_mult)
  call f_free(cutoff)

end subroutine init_sparse_matrix_wrapper


!> Initializes a sparse matrix type compatible with the ditribution of the KS orbitals
subroutine init_sparse_matrix_for_KSorbs(iproc, nproc, orbs, input, nextra, smat, smat_extra)
  use module_types
  use module_interfaces
  use sparsematrix_base, only: sparse_matrix
  use sparsematrix_init, only: init_sparse_matrix
  implicit none

  ! Calling arguments
  integer,intent(in) :: iproc, nproc, nextra
  type(orbitals_data),intent(in) :: orbs
  type(input_variables),intent(in) :: input
  type(sparse_matrix),intent(out) :: smat, smat_extra

  ! Local variables
  integer :: i, iorb, iiorb, jorb, ind
  integer,dimension(:),allocatable :: nonzero
  type(orbitals_data) :: orbs_aux
  character(len=*),parameter :: subname='init_sparse_matrix_for_KSorbs'

  call f_routine('init_sparse_matrix_for_KSorbs')

  ! First the type for the normal KS orbitals distribution
  nonzero = f_malloc(orbs%norb*orbs%norbp, id='nonzero')
  i=0
  do iorb=1,orbs%norbp
      iiorb=orbs%isorb+iorb
      do jorb=1,orbs%norb
          i=i+1
          ind=(iiorb-1)*orbs%norb+jorb
          nonzero(i)=ind
      end do
  end do
  call init_sparse_matrix(iproc, nproc, orbs%norb, orbs%norbp, orbs%isorb, input%store_index, &
       orbs%norb*orbs%norbp, nonzero, orbs%norb, nonzero, smat, print_info_=.false.)
  call f_free(nonzero)


  ! Now the distribution for the KS orbitals including the extr states. Requires
  ! first to calculate a corresponding orbs type.
  call nullify_orbitals_data(orbs_aux)
  call orbitals_descriptors(iproc, nproc, orbs%norb+nextra, orbs%norb+nextra, 0, input%nspin, orbs%nspinor,&
       input%gen_nkpt, input%gen_kpt, input%gen_wkpt, orbs_aux, .false.)
  nonzero = f_malloc(orbs_aux%norb*orbs_aux%norbp, id='nonzero')
  i=0
  do iorb=1,orbs_aux%norbp
      iiorb=orbs_aux%isorb+iorb
      do jorb=1,orbs_aux%norb
          i=i+1
          ind=(iiorb-1)*orbs_aux%norb+jorb
          nonzero(i)=ind
      end do
  end do
  call init_sparse_matrix(iproc, nproc, orbs_aux%norb, orbs_aux%norbp, orbs_aux%isorb, input%store_index, &
       orbs_aux%norb*orbs_aux%norbp, nonzero, orbs_aux%norb, nonzero, smat_extra, print_info_=.false.)
  call f_free(nonzero)
  call deallocate_orbitals_data(orbs_aux, subname)

  call f_release_routine()

end subroutine init_sparse_matrix_for_KSorbs
