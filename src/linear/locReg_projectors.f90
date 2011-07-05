!#############################################################################################################################################
!!****f* BigDFT/nlpspd_to_locreg
!#############################################################################################################################################
!! FUNCTION: Transform projector descriptors between Global region and localisation region
!!
!! WARNING: 
!!         
!! SOURCE:
!!
subroutine nlpspd_to_locreg(input_parameters,iproc,Glr,Llr,rxyz,atoms,orbs,&
       radii_cf,cpmult,fpmult,hx,hy,hz,nlpspd,Lnlpspd,projflg)

  use module_base
  use module_types
 
 implicit none

  !#######################################
  ! Subroutine Scalar Arguments
  !#######################################
  type(input_variables),intent(in) :: input_parameters
  integer,intent(in) :: iproc
  type(locreg_descriptors),intent(in) :: Glr  ! Global grid descriptor
  type(locreg_descriptors),intent(in) :: Llr  ! Local grid descriptor
  type(atoms_data),intent(in) :: atoms        ! atom descriptors
  type(orbitals_data),intent(in) :: orbs      ! orbital descriptors
  real(gp), intent(in) :: cpmult,fpmult,hx,hy,hz  ! grid descriptions
  type(nonlocal_psp_descriptors),intent(in) :: nlpspd  ! global descriptors for the projectors
  type(nonlocal_psp_descriptors),intent(out) :: Lnlpspd  ! local descriptors for the projectors 
  !########################################
  !Subroutine Array Arguments
  !########################################
  integer,dimension(atoms%nat),intent(out) :: projflg  ! atoms contributing projectors inside the locreg
  real(gp), dimension(3,atoms%nat), intent(in) :: rxyz !atomic positions
  real(gp), dimension(atoms%ntypes,3), intent(in) :: radii_cf  ! radii of the different atom types
  !#############################################
  !local variables
  !############################################
   integer :: iatom,igrid
   integer :: ii,jj,i1,i2,i3      !integers for loops
   integer :: mproj,natp
   integer :: mseg_c,mvctr_c,mseg_f,mvctr_f
   integer :: mseg !total number of segments
   integer :: iat  ! index of the atoms
   integer :: nprojelat ! total number of elements
   integer :: isx,isy,isz,iex,iey,iez
   integer :: iseg,jseg,Gseg,Gvctr
   integer :: nl1,nl2,nl3,nu1,nu2,nu3 ! bounds of projectors around atom iatom
   integer,dimension(1:2,1:2,1:3) :: bounds
   logical,dimension(0:Glr%d%n1,0:Glr%d%n2,0:Glr%d%n3) :: logrid
   character(len=*),parameter :: subname='nlpspd_to_locreg'

!Determine the number of projectors with components in locreg
! and also which atoms have such projectors and number of atoms
   call  number_of_projectors_in_locreg(atoms,cpmult,fpmult,Glr,hx,hy,hz,Llr,nlpspd,&
&        mproj,projflg,natp,radii_cf,rxyz)

   Lnlpspd%nproj = mproj

!DEBUG
!  print *,'Llr check:',Llr%ns1,Llr%ns2,Llr%ns3,Llr%d%n1,Llr%d%n2,Llr%d%n3
!  print *,'Number of projectors:', mproj
!  print *,'Projflg', projflg
!END DEBUG

!Allocate the arrays of Lnlpspd, except keyg_p and keyv_p
 call allocate_Lnlpspd(natp,Lnlpspd,subname)

  iat = 0
  nprojelat = 0
  Lnlpspd%nprojel = 0
  mseg = 0
  do iatom = 1,atoms%nat
     if(projflg(iatom) == 0) cycle 
     iat = iat + 1  !iatom is the global numbering of atoms, while iat is the numbering only in the locreg. 

!    Determine the bounds of the projectors
     call projector_box_in_locreg(iatom,Glr,Llr,nlpspd,bounds)

     do ii = 1,2
        do jj = 1,3
           Lnlpspd%nboxp_c(ii,jj,iat) = bounds(1,ii,jj)
           Lnlpspd%nboxp_f(ii,jj,iat) = bounds(2,ii,jj)
        end do
     end do

!    Rename the variables
     nl1 = nlpspd%nboxp_c(1,1,iatom)
     nu1 = nlpspd%nboxp_c(2,1,iatom)
     nl2 = nlpspd%nboxp_c(1,2,iatom)
     nu2 = nlpspd%nboxp_c(2,2,iatom)
     nl3 = nlpspd%nboxp_c(1,3,iatom)
     nu3 = nlpspd%nboxp_c(2,3,iatom)

!    Now we can determine the number of segments and elements of coarse grid
     call fill_logrid(atoms%geocode,Glr%d%n1,Glr%d%n2,Glr%d%n3,nl1,nu1,nl2,nu2,nl3,nu3,0,1,atoms%ntypes,&
&                     atoms%iatype(iatom),rxyz(1,iatom),radii_cf(:,3),cpmult,hx,hy,hz,logrid)

     call number_of_projector_elements_in_locreg(iatom,1,atoms,Glr,Llr,logrid,nlpspd,mproj,mseg_c,mvctr_c)

     Lnlpspd%nseg_p(2*iat-1) = mseg_c
     Lnlpspd%nvctr_p(2*iat-1) = mvctr_c 

! Do the same for fine grid

!    Rename the variables
     nl1 = nlpspd%nboxp_f(1,1,iatom)
     nu1 = nlpspd%nboxp_f(2,1,iatom)
     nl2 = nlpspd%nboxp_f(1,2,iatom)
     nu2 = nlpspd%nboxp_f(2,2,iatom)
     nl3 = nlpspd%nboxp_f(1,3,iatom)
     nu3 = nlpspd%nboxp_f(2,3,iatom)

     call fill_logrid(atoms%geocode,Glr%d%n1,Glr%d%n2,Glr%d%n3,nl1,nu1,nl2,nu2,nl3,nu3,0,1,atoms%ntypes,&
&                     atoms%iatype(iatom),rxyz(1,iatom),radii_cf(:,2),fpmult,hx,hy,hz,logrid)

     call number_of_projector_elements_in_locreg(iatom,2,atoms,Glr,Llr,logrid,nlpspd,mproj,mseg_f,mvctr_f)

     Lnlpspd%nseg_p(2*iat) = mseg_f
     Lnlpspd%nvctr_p(2*iat) = mvctr_f

!    Should not be useful, because if projflg is > 0 there should be some elements
!     if(mvctr_c == 0 .and. mvctr_f == 0) then
!        projflg(iatom) = 0 
!     end if    

     nprojelat = mvctr_c*projflg(iatom) + 7*mvctr_f*projflg(iatom)
     Lnlpspd%nprojel = max(Lnlpspd%nprojel,nprojelat)
     mseg = mseg + mseg_c + mseg_f
  end do

! Now allocate keyg_p,keyv_p following the needs
  call allocate_projd(mseg,Lnlpspd,subname)

! Renaming some variables to simply calling of routines
  !starting point of locreg
  isx = Llr%ns1
  isy = Llr%ns2
  isz = Llr%ns3
  !ending point of locreg
  iex = Llr%ns1 + Llr%d%n1
  iey = Llr%ns2 + Llr%d%n2
  iez = Llr%ns3 + Llr%d%n3
  
! At last, fill the projector descriptors (keyg_p,keyv_p)
  iseg = 1
  iat = 0
  do iatom= 1,atoms%nat
     if(projflg(iatom) == 0) cycle
     iat = iat + 1

!    number of segments for coarse
     jseg = nlpspd%nseg_p(2*iatom-2)+1 ! index where to start in keyg for global region (nlpspd)
     Gseg = nlpspd%nseg_p(2*iatom-1)-nlpspd%nseg_p(2*iatom-2) ! number of segments for global region
     Gvctr = nlpspd%nvctr_p(2*iatom-1)-nlpspd%nvctr_p(2*iatom-2)!number of elements for global region

!    Coarse part 
     if (Lnlpspd%nseg_p(2*iat-1) > 0) then
        call segkeys_loc(Glr%d%n1,Glr%d%n2,Glr%d%n3,isx,iex,isy,iey,isz,iez,&
            Gseg,Gvctr,nlpspd%keyg_p(1,jseg),nlpspd%keyv_p(jseg),&
            Lnlpspd%nseg_p(2*iat-1),Lnlpspd%nvctr_p(2*iat-1),&
            Lnlpspd%keyg_p(1,iseg),Lnlpspd%keyv_p(iseg))
     end if

     iseg = iseg + Lnlpspd%nseg_p(2*iat-1)      
     if(Lnlpspd%nseg_p(2*iat) > 0) then  !only do fine grid if present
!    Number of segments for fine
        jseg = nlpspd%nseg_p(2*iatom-1)+1 ! index where to start in keyg for global region (nlpspd)
        Gseg = nlpspd%nseg_p(2*iatom)-nlpspd%nseg_p(2*iatom-1) ! number of segments for global region
        Gvctr = nlpspd%nvctr_p(2*iatom)-nlpspd%nvctr_p(2*iatom-1)!number of elements for global region

!    Fine part 
        call segkeys_loc(Glr%d%n1,Glr%d%n2,Glr%d%n3,isx,iex,isy,iey,isz,iez,&
             Gseg,Gvctr,nlpspd%keyg_p(1,jseg),nlpspd%keyv_p(jseg),&
             Lnlpspd%nseg_p(2*iat),Lnlpspd%nvctr_p(2*iat),&
             Lnlpspd%keyg_p(1,iseg),Lnlpspd%keyv_p(iseg))

        iseg = iseg + Lnlpspd%nseg_p(2*iat)
     end if 
  end do

END SUBROUTINE nlpspd_to_locreg
!%***

!#############################################################################################################################################
!!****f* BigDFT/number_of_projectors_in_locreg
!#############################################################################################################################################
!! FUNCTION: Calculates the number of projectors with components in the locreg
!!           It also returns a vector, projflg, which identifies the atoms with projectors inside the region
!!                      projflg = 0, no projectors in locreg
!!                      projflg = nproj, nproj projectors from atom iatom in locreg
!! WARNING: 
!!         
!! SOURCE:
!!
subroutine number_of_projectors_in_locreg(atoms,cpmult,fpmult,Glr,hx,hy,hz,Llr,nlpspd,&
&          mproj,projflg,natp,radii_cf,rxyz)

  use module_base
  use module_types
 
 implicit none

  !#######################################
  ! Subroutine Scalar Arguments
  !#######################################
  real(gp),intent(in) :: cpmult,fpmult,hx,hy,hz  ! grid descriptions
  type(atoms_data),intent(in) :: atoms        ! atoms descriptor
  type(locreg_descriptors),intent(in) :: Glr  ! Global grid descriptor
  type(locreg_descriptors),intent(in) :: Llr  ! Local grid descriptor
  type(nonlocal_psp_descriptors),intent(in) :: nlpspd  ! global descriptors for the projectors
  integer,intent(out) :: mproj  ! number of projectors
  integer,intent(out) :: natp   ! number of atoms having projectors in region
  !#######################################
  ! Subroutine Array Arguments
  !#######################################
  integer,dimension(atoms%nat),intent(out) :: projflg ! flag which is equal to the number of projectors with components inside locreg for each atom
  real(gp), dimension(3,atoms%nat), intent(in) :: rxyz !atomic positions
  real(gp), dimension(atoms%ntypes,3), intent(in) :: radii_cf  ! radii of the different atom types
  !#######################################
  ! Local Variables
  !#######################################
  integer :: iatom,ii,izone                 ! integer for loop
  integer :: bound(1:2,1:3)           ! bound of locreg
  integer :: nproj                    ! temporary number of projectors
  integer :: i_stat                   ! allocation error
  logical :: intersect                ! logical for intersect of projector with locreg
  character(len=*), parameter :: subname='number_of_projectors_in_locreg'

  projflg = 0
  mproj = 0
  natp = 0
  do iatom=1,atoms%nat
!       check if projector of atom iatom overlap the locreg (coarse grid)        
        call check_projector_intersect_with_locreg(atoms,cpmult,Glr,hx,hy,hz,iatom,Llr,&
&            radii_cf(atoms%iatype(iatom),3),rxyz,intersect)

        if(intersect) then
           call numb_proj(atoms%iatype(iatom),atoms%ntypes,atoms%psppar,atoms%npspcode,nproj)
           mproj = mproj + nproj
           if(nproj > 0) then
              projflg(iatom) = nproj
              natp = natp + 1
           end if
        end if        

! Only have to do it if the atom is not yet selected
     if (projflg(iatom) .eq. 0) then
!         check if projector of atom iatom overlap the locreg (fine grid)
          call check_projector_intersect_with_locreg(atoms,fpmult,Glr,hx,hy,hz,iatom,Llr,&
&              radii_cf(atoms%iatype(iatom),2),rxyz,intersect)

          if(intersect) then        
             call numb_proj(atoms%iatype(iatom),atoms%ntypes,atoms%psppar,atoms%npspcode,nproj)
             mproj = mproj + nproj
             if(nproj > 0) projflg(iatom) = nproj
          end if
     end if        
  end do

END SUBROUTINE number_of_projectors_in_locreg
!%***

!#############################################################################################################################################
!!****f* BigDFT/check_projector_intersect_with_locreg
!#############################################################################################################################################
!! FUNCTION: Returns the limits of the various folded projector zones.
!!           
!! WARNING: Works only for overlaps (i.e. boundaries must be inside simulation box) 
!!         
!! SOURCE:
!!
subroutine check_projector_intersect_with_locreg(atoms,pmult,Glr,hx,hy,hz,iatom,Llr,radii_cf,rxyz,intersect)

  use module_base
  use module_types
 
  implicit none

  !#######################################
  ! Subroutine Scalar Arguments
  !#######################################
  integer,intent(in) :: iatom     !number of atom we are treating
  real(gp),intent(in) :: hx,hy,hz   ! grid spacing
  real(gp),intent(in) :: pmult     ! factor for the radius of projector
  real(gp),intent(in) :: radii_cf  ! radii of the atom type
  type(atoms_data),intent(in) :: atoms        ! atoms descriptor
  type(locreg_descriptors),intent(in) :: Glr ! global region descriptor
  type(locreg_descriptors),intent(in) :: Llr ! local region descriptor
  logical,intent(out) :: intersect !.true. if projector intersects zone
  !#######################################
  ! Subroutine Array Arguments
  !#######################################
  real(gp), dimension(3,atoms%nat), intent(in) :: rxyz !atomic positions
  !#######################################
  ! Local Variables
  !#######################################
  integer :: i1,i2,i3  !integer for loops
  real(gp) :: dx1,dx2,dy1,dy2,dz1,dz2 !two distance in X,Y,Z
  real(gp) :: rad !radius of projectors

  intersect = .false.
  rad=radii_cf*pmult 

!Check if zone is within the radius of the projectors
! Must also check the images in other cells  
  do i3 = Llr%ns3,Llr%ns3+Llr%d%n3
     dz1 = (real(i3,gp)*hz-rxyz(3,iatom))**2
     if (Glr%geocode == 'S' .or. Glr%geocode =='P') then
        dz2 = (real(i3,gp)*hz-(rxyz(3,iatom)+ (Glr%d%n3+1)*hz))**2 !translating to positive
        dz1 = min(dz1,dz2) 
        dz2 = (real(i3,gp)*hz-(rxyz(3,iatom)- (Glr%d%n3+1)*hz))**2 !translating to negative
        dz1 = min(dz1,dz2)
     end if

     do i2 = Llr%ns2,Llr%ns2+Llr%d%n2
        dy1 = (real(i2,gp)*hy-rxyz(2,iatom))**2
        if (Glr%geocode == 'P') then
           dy2 = (real(i2,gp)*hy-(rxyz(2,iatom)+ (Glr%d%n2+1)*hy))**2 !translating to positive
           dy1 = min(dy1,dy2) 
           dy2 = (real(i2,gp)*hy-(rxyz(2,iatom)- (Glr%d%n2+1)*hy))**2 !translating to negative
           dy1 = min(dy1,dy2)
        end if

        do i1 = Llr%ns1,Llr%ns1+Llr%d%n1
           dx1 = (real(i1,gp)*hx-rxyz(1,iatom))**2
           if (Glr%geocode == 'S' .or. Glr%geocode =='P') then
              dx2 = (real(i1,gp)*hx-(rxyz(1,iatom)+ (Glr%d%n1+1)*hx))**2 !translating to positive
              dx1 = min(dx1,dx2) 
              dx2 = (real(i1,gp)*hx-(rxyz(1,iatom)- (Glr%d%n1+1)*hx))**2 !translating to negative
              dx1 = min(dx1,dx2)
           end if

           if(dx1+dy1+dz1 <= rad**2) then
             intersect = .true.
             exit
           end if
        end do
        if (intersect) exit
     end do
        if (intersect) exit
  end do

END SUBROUTINE check_projector_intersect_with_locreg
!%***

!#############################################################################################################################################
!!****f* BigDFT/number_of_projector_elements_in_locreg
!#############################################################################################################################################
!! FUNCTION: Calculates the number of segments (mseg) and elements (mvctr) of projectors in locreg
!!           
!! WARNING: 
!!         
!! SOURCE:
!!
subroutine number_of_projector_elements_in_locreg(iatom,igrid,atoms,Glr,Llr,logrid,nlpspd,mproj,mseg,mvctr)

  use module_base
  use module_types
 
 implicit none

  !#######################################
  ! Subroutine Scalar Arguments
  !#######################################
  integer,intent(in) :: iatom  ! current atom
  integer,intent(in) :: igrid  ! treat coarse (1) or fine (2) grid
  type(atoms_data),intent(in) :: atoms        ! atoms descriptor
  type(locreg_descriptors),intent(in) :: Glr  ! Global grid descriptor
  type(locreg_descriptors),intent(in) :: Llr  ! Local grid descriptor
  type(nonlocal_psp_descriptors),intent(in) :: nlpspd  ! global descriptors for the projectors
  integer,intent(in) :: mproj  ! number of projectors
  integer,intent(out):: mseg   ! number of segments
  integer,intent(out):: mvctr  ! number of elements
  !#######################################
  ! Subroutine Array Arguments
  !#######################################
  logical, dimension(0:Glr%d%n1,0:Glr%d%n2,0:Glr%d%n3), intent(in) :: logrid
  !#######################################
  ! Local Variables
  !#######################################
  integer :: i1,i2,i3  ! integers for loops
  integer :: nl1,nl2,nl3,nu1,nu2,nu3   ! rename the bounds of projectors coarse grid
  integer :: nend,nsrt,mvctri,nsrti,nendi
  integer,dimension(1:2,1:3) :: bound  ! rename the bounds of locreg
  logical :: plogrid   ! logical to check start of new segment


! Set boundaries of projectors (coarse)
   if(igrid == 1) then
      nl1 = nlpspd%nboxp_c(1,1,iatom)
      nl2 = nlpspd%nboxp_c(1,2,iatom)
      nl3 = nlpspd%nboxp_c(1,3,iatom)

      nu1 = nlpspd%nboxp_c(2,1,iatom)
      nu2 = nlpspd%nboxp_c(2,2,iatom)
      nu3 = nlpspd%nboxp_c(2,3,iatom)

   end if

! Set boundaries of projectors (fine)
   if(igrid == 2) then
      nl1 = nlpspd%nboxp_f(1,1,iatom)
      nl2 = nlpspd%nboxp_f(1,2,iatom)
      nl3 = nlpspd%nboxp_f(1,3,iatom)

      nu1 = nlpspd%nboxp_f(2,1,iatom)
      nu2 = nlpspd%nboxp_f(2,2,iatom)
      nu3 = nlpspd%nboxp_f(2,3,iatom)
   end if

! bounds of the localization region
  !lower bound
  bound(1,1) = Llr%ns1 - Glr%ns1
  bound(1,2) = Llr%ns2 - Glr%ns2
  bound(1,3) = Llr%ns3 - Glr%ns3

  !upper bound  (WILL NOT WORK FOR PERIODICITY)
  bound(2,1) = Llr%ns1 + Llr%d%n1 - Glr%ns1
  bound(2,2) = Llr%ns2 + Llr%d%n2 - Glr%ns2
  bound(2,3) = Llr%ns3 + Llr%d%n3 - Glr%ns3 

  if (igrid == 1) then
!Initialize counters
     mvctr=0
     nsrt=0
     nend=0
     mvctri=0
     nsrti=0
     nendi=0

! Do coarse grid
     do i3=nl3,nu3
        if(i3 > bound(2,3) .or. i3 < bound(1,3)) cycle
        do i2=nl2,nu2
           if(i2 > bound(2,2) .or. i2 < bound(1,2)) cycle
           plogrid=.false.
           do i1=nl1,nu1
              if(i1 > bound(2,1) .or. i1 < bound(1,1))cycle

              if(logrid(i1,i2,i3)) then
                 mvctri=mvctri+1
                 if (.not. plogrid) then
                    nsrti=nsrti+1
                 endif
              else
                 if(plogrid) then
                    nendi=nendi+1
                 endif
              endif
              plogrid=logrid(i1,i2,i3)
           enddo
           if (i2 .le. bound(2,2) .and. i2 .ge. bound(1,2) .and. &
&              i3 .le. bound(2,3) .and. i3 .ge. bound(1,3) .and. &
               plogrid .eqv. .true.) then
              nendi=nendi+1
           endif
        enddo
     enddo

     mvctr=mvctr+mvctri
     nsrt=nsrt+nsrti
     nend=nend+nendi

     if (nend /= nsrt) then
        write(*,*)' ERROR in number_of_projector_elements_in_locreg : nend <> nsrt',nend,nsrt
        stop
     endif
     mseg=nend
  end if

  if(igrid == 2) then

     !Initialize counters
     mvctr=0
     nsrt=0
     nend=0
     mvctri=0
     nsrti=0
     nendi=0

! Do fine grid
     do i3=nl3,nu3
        if(i3 > bound(2,3) .or. i3 < bound(1,3)) cycle
        do i2=nl2,nu2
           if(i2 > bound(2,2) .or. i2 < bound(1,2)) cycle
           plogrid=.false.
           do i1=nl1,nu1
              if(i1 > bound(2,1) .or. i1 < bound(1,1))cycle

              if(logrid(i1,i2,i3)) then
                 mvctri=mvctri+1
                 if (.not. plogrid) then
                    nsrti=nsrti+1
                 endif
              else
                 if(plogrid) then
                    nendi=nendi+1
                 endif
              endif
              plogrid=logrid(i1,i2,i3)
           enddo
           if (i2 .le. bound(2,2) .and. i2 .ge. bound(1,2) .and. &
&              i3 .le. bound(2,3) .and. i3 .ge. bound(1,3) .and. &
               plogrid .eqv. .true.) then
              nendi=nendi+1
           endif
        enddo
     enddo

     mvctr=mvctr+mvctri
     nsrt=nsrt+nsrti
     nend=nend+nendi

     if (nend /= nsrt) then
        write(*,*)' ERROR in number_of_projector_elements_in_locreg (fine) : nend <> nsrt',nend,nsrt
        stop
     endif
     mseg=nend
  end if

END SUBROUTINE number_of_projector_elements_in_locreg
!%***


!#############################################################################################################################################
!!****f* BigDFT/ projector_box_in_locreg
!#############################################################################################################################################
!! FUNCTION: Calculates the bounds of the box of the projector in locreg
!!           bounds(1,:,:) for coarse grid
!!           bounds(2,:,:) for fine grid
!! WARNING: 
!!         
!! SOURCE:
!!
subroutine projector_box_in_locreg(iatom,Glr,Llr,nlpspd,bounds)

  use module_base
  use module_types
 
 implicit none

  !#######################################
  ! Subroutine Scalar Arguments
  !#######################################
  integer,intent(in) :: iatom  ! current atom
  type(locreg_descriptors),intent(in) :: Glr  ! Global grid descriptor
  type(locreg_descriptors),intent(in) :: Llr  ! Local grid descriptor
  type(nonlocal_psp_descriptors),intent(in) :: nlpspd  ! global descriptors for the projectors
  integer,dimension(1:2,1:2,1:3),intent(out) :: bounds
  !#######################################
  ! Local Variables
  !#######################################
  integer :: ii
  integer,dimension(1:2,1:3) :: Cnl,Fnl,Lnl
  
! Set boundaries of projectors (coarse)
  !lower bounds
  Cnl(1,1) = nlpspd%nboxp_c(1,1,iatom)
  Cnl(1,2) = nlpspd%nboxp_c(1,2,iatom)
  Cnl(1,3) = nlpspd%nboxp_c(1,3,iatom)
  !upper bounds
  Cnl(2,1) = nlpspd%nboxp_c(2,1,iatom)
  Cnl(2,2) = nlpspd%nboxp_c(2,2,iatom)
  Cnl(2,3) = nlpspd%nboxp_c(2,3,iatom)

! Set boundaries of projectors (fine)
  !lower bounds
  Fnl(1,1) = nlpspd%nboxp_f(1,1,iatom)
  Fnl(1,2) = nlpspd%nboxp_f(1,2,iatom)
  Fnl(1,3) = nlpspd%nboxp_f(1,3,iatom)
  !upper bounds
  Fnl(2,1) = nlpspd%nboxp_f(2,1,iatom)
  Fnl(2,2) = nlpspd%nboxp_f(2,2,iatom)
  Fnl(2,3) = nlpspd%nboxp_f(2,3,iatom)

! bounds of the localization region
  !lower bounds
  Lnl(1,1) = Llr%ns1 - Glr%ns1
  Lnl(1,2) = Llr%ns2 - Glr%ns2
  Lnl(1,3) = Llr%ns3 - Glr%ns3
  !upper bounds
  Lnl(2,1) = Llr%ns1 + Llr%d%n1 - Glr%ns1
  Lnl(2,2) = Llr%ns2 + Llr%d%n2 - Glr%ns2
  Lnl(2,3) = Llr%ns3 + Llr%d%n3 - Glr%ns3

! Calculate the bounds for the grids
  do ii=1,3
     !lower bounds
     bounds(1,1,ii) = max(Lnl(1,ii),Cnl(1,ii))
     bounds(2,1,ii) = max(Lnl(1,ii),Fnl(1,ii))
     ! upper bounds
     bounds(1,2,ii) = min(Lnl(2,ii),Cnl(2,ii))
     bounds(2,2,ii) = min(Lnl(2,ii),Fnl(2,ii))
  end do

END SUBROUTINE projector_box_in_locreg
!%***

!#############################################################################################################################################
!!****f* BigDFT/allocate_Lnlpspd
!#############################################################################################################################################
!! FUNCTION:  Allocates most of the arrays in Lnlpspd 
!!
!! WARNING: 
!!         
!! SOURCE:
!!
subroutine allocate_Lnlpspd(natom,Lnlpspd,subname)

  use module_base
  use module_types
 
 implicit none

  !#######################################
  ! Subroutine Scalar Arguments
  !#######################################
  integer,intent(in) :: natom
  type(nonlocal_psp_descriptors),intent(inout) :: Lnlpspd  ! Local descriptors for the projectors
  character(len=*), intent(in) :: subname
  !#######################################
  ! Local Variables 
  !#######################################
  integer :: i_stat

  allocate(Lnlpspd%nvctr_p(2*natom+ndebug),stat=i_stat)
  call memocc(i_stat,Lnlpspd%nvctr_p,'nvctr_p',subname)
  allocate(Lnlpspd%nseg_p(2*natom+ndebug),stat=i_stat)
  call memocc(i_stat,Lnlpspd%nseg_p,'nseg_p',subname)
  allocate(Lnlpspd%nboxp_c(2,3,2*natom),stat=i_stat)
  call memocc(i_stat,Lnlpspd%nboxp_c,'nbox_c',subname)
  allocate(Lnlpspd%nboxp_f(2,3,2*natom),stat=i_stat)
  call memocc(i_stat,Lnlpspd%nboxp_f,'nbox_f',subname)

END SUBROUTINE allocate_Lnlpspd
!%***

!#############################################################################################################################################
!!****f* BigDFT/allocate_Lnlpspd
!#############################################################################################################################################
!! FUNCTION:  Deallocates most of the arrays in Lnlpspd 
!!
!! WARNING: 
!!         
!! SOURCE:
!!
subroutine deallocate_Lnlpspd(Lnlpspd,subname)

  use module_base
  use module_types

 implicit none

  !#######################################
  ! Subroutine Scalar Arguments
  !#######################################
  type(nonlocal_psp_descriptors),intent(inout) :: Lnlpspd  ! Local descriptors for the projectors
  character(len=*), intent(in) :: subname
  !#######################################
  ! Local Variables 
  !#######################################
  integer :: i_stat,i_all

  i_all=-product(shape(Lnlpspd%nvctr_p)*kind(Lnlpspd%nvctr_p))
  deallocate(Lnlpspd%nvctr_p,stat=i_stat)
  call memocc(i_stat,i_all,'nvctr_p',subname)
  i_all=-product(shape(Lnlpspd%nseg_p)*kind(Lnlpspd%nseg_p))
  deallocate(Lnlpspd%nseg_p,stat=i_stat)
  call memocc(i_stat,i_all,'nseg_p',subname)
  i_all=-product(shape(Lnlpspd%nboxp_c)*kind(Lnlpspd%nboxp_c))
  deallocate(Lnlpspd%nboxp_c,stat=i_stat)
  call memocc(i_stat,i_all,'nboxp_c',subname)
  i_all=-product(shape(Lnlpspd%nboxp_f)*kind(Lnlpspd%nboxp_f))
  deallocate(Lnlpspd%nboxp_f,stat=i_stat)
  call memocc(i_stat,i_all,'nboxp_f',subname)
  i_all=-product(shape(Lnlpspd%keyg_p)*kind(Lnlpspd%keyg_p))
  deallocate(Lnlpspd%keyg_p,stat=i_stat)
  call memocc(i_stat,i_all,'keyg_p',subname)
  i_all=-product(shape(Lnlpspd%keyv_p)*kind(Lnlpspd%keyv_p))
  deallocate(Lnlpspd%keyv_p,stat=i_stat)
  call memocc(i_stat,i_all,'keyv_p',subname)

END SUBROUTINE deallocate_Lnlpspd
!%***



!#############################################################################################################################################
!!****f* BigDFT/allocate_projd
!#############################################################################################################################################
!! FUNCTION: allocates the keyg_p and keyv_p descriptors for the projectors
!!          
!!           
!! WARNING: 
!!         
!! SOURCE:
!!
subroutine allocate_projd(mseg,Lnlpspd,subname)

  use module_base
  use module_types
 
 implicit none

  !#######################################
  ! Subroutine Scalar Arguments
  !#######################################
  integer,intent(in) :: mseg
  type(nonlocal_psp_descriptors),intent(inout) :: Lnlpspd  ! Local descriptors for the projectors
  character(len=*), intent(in) :: subname
  !#######################################
  ! Local Variables 
  !#######################################
  integer :: i_stat
  allocate(Lnlpspd%keyg_p(2,mseg),stat=i_stat)
  call memocc(i_stat,Lnlpspd%keyg_p,'keyg_p',subname)
  allocate(Lnlpspd%keyv_p(mseg),stat=i_stat)
  call memocc(i_stat,Lnlpspd%keyv_p,'keyv_p',subname)

END SUBROUTINE allocate_projd
!%***


!#############################################################################################################################################
!!****f* BigDFT/apply_local_projectors
!#############################################################################################################################################
!! FUNCTION: Fills the projector pointer and applies the projectors to the wavefunctions
!!           
!!           
!! WARNING: 
!!         
!! SOURCE:
!!
subroutine apply_local_projectors(ilr,nspin,atoms,hx,hy,hz,Llr,Lnlpspd,orbs,Gorbs,projflg,psi,rxyz,hpsi,eproj)


  use module_base
  use module_types
  !use module_interfaces, exceptThisOne => apply_local_projectors
 
  implicit none

  !#######################################
  ! Subroutine Scalar Arguments
  !#######################################
  integer, intent(in) :: ilr,nspin
  real(gp), intent(in) :: hx,hy,hz
  type(atoms_data),intent(in) :: atoms
  type(locreg_descriptors),intent(in) :: Llr
  type(nonlocal_psp_descriptors),intent(in) :: Lnlpspd  ! Local descriptors for the projectors
  type(orbitals_data),intent(in) :: orbs
  type(orbitals_data),intent(in) :: Gorbs
  real(gp), intent(inout) :: eproj
  !#######################################
  ! Subroutine Array Arguments
  !#######################################
  integer,dimension(atoms%nat),intent(in) :: projflg
  real(wp),dimension((Llr%wfd%nvctr_c+7*Llr%wfd%nvctr_f)*orbs%nspinor*orbs%norb),intent(in) :: psi  !local wavefunction
  real(wp),dimension((Llr%wfd%nvctr_c+7*Llr%wfd%nvctr_f)*orbs%nspinor*orbs%norb),intent(inout):: hpsi ! local |p><p|Psi>
  real(gp), dimension(3,atoms%nat), intent(in) :: rxyz
  !#######################################
  ! Local Variables 
  !#######################################
  integer :: ikpt,istart_c,ncplx,jseg_c,iproj,iat,ityp,l,i,nwarnings
  integer :: isorb,ieorb,nspinor,iorb,istart_o,ispinor
  integer :: nels,ipsi,ii,iatom,iel,i_all,i_stat
  integer :: jj,orbtot,ispin,ind
  integer,dimension(Llr%localnorb*nspin) :: inthisLocreg
  !integer,dimension(:),allocatable :: inthisLocreg
  real(gp) :: kx,ky,kz,eproj_spinor
  real(wp),allocatable,dimension(:,:,:) :: psi_tmp
  real(wp),allocatable,dimension(:,:,:) :: hpsi_tmp
  real(wp),allocatable,dimension(:):: Lproj  !local projectors
  character(len=*), parameter :: subname='apply_local_projectors'


!  First reshape the wavefunctions: psi_tmp(nels,norbs,nspinor)
   nels = Llr%wfd%nvctr_c+7*Llr%wfd%nvctr_f

! Allocate arrays
   allocate(psi_tmp(nels,orbs%nspinor,orbs%norb),stat=i_stat)
   call memocc(i_stat,psi_tmp,'psi_tmp',subname)
   allocate(hpsi_tmp(nels,orbs%nspinor,orbs%norb),stat=i_stat)
   call memocc(i_stat,hpsi_tmp,'hpsi_tmp',subname)
   allocate(Lproj(Lnlpspd%nprojel),stat=i_stat)
   call memocc(i_stat,Lproj,'Lproj',subname)

   ! reshape the wavefunction
   ii=0
   do iorb=1,orbs%norb
       do ispinor=1,orbs%nspinor
           do iel=1,nels
               ii=ii+1
               psi_tmp(iel,ispinor,iorb)=psi(ii)
               hpsi_tmp(iel,ispinor,iorb)=hpsi(ii)
           end do
       end do
   end do
   
   ieorb = orbs%norbp   ! give an initial value because could skip whole loop on atoms (i.e. Li+ test)
   ikpt=orbs%iokpt(1)
   loop_kpt: do
      !features of the k-point ikpt
      kx=orbs%kpts(1,ikpt)
      ky=orbs%kpts(2,ikpt)
      kz=orbs%kpts(3,ikpt)

      !evaluate the complexity of the k-point
      if (kx**2 + ky**2 + kz**2 == 0.0_gp) then
         ncplx=1
      else
         ncplx=2
      end if

      ieorb = Gorbs%norbp  !initialize value in case no atoms have projectors
      jseg_c = 1
      iproj = 0
      iatom = 0
      do iat = 1,atoms%nat
         if(projflg(iat) == 0) cycle
         iatom = iatom +1
         istart_c = 1
         ityp=atoms%iatype(iat)

         do l=1,4 !generic case, also for HGHs (for GTH it will stop at l=2)
            do i=1,3 !generic case, also for HGHs (for GTH it will stop at i=2)
               if (atoms%psppar(l,i,ityp) /= 0.0_gp) then

!                 Second fill the projectors
!                 NOTE : idir was set to 0 because we don't care for derivatives
                  call local_projector(atoms%geocode,atoms%atomnames(ityp),iat,0,l,i,&
                       atoms%psppar(l,0,ityp),rxyz(1,iat),Llr,&
                       hx,hy,hz,kx,ky,kz,ncplx,Lnlpspd%nvctr_p(2*iatom-1),&
                       Lnlpspd%nvctr_p(2*iatom),Lnlpspd%nseg_p(2*iatom-1),Lnlpspd%nseg_p(2*iatom),&
                       Lnlpspd%keyv_p(jseg_c),Lnlpspd%keyg_p(1,jseg_c),Lproj(istart_c),nwarnings)

                  iproj=iproj+2*l-1
                  istart_c=istart_c+(Lnlpspd%nvctr_p(2*iatom-1)+7*Lnlpspd%nvctr_p(2*iatom))*(2*l-1)*ncplx
                  !print *,'iproc,istart_c,nlpspd%nprojel',istart_c,Lnlpspd%nprojel,ncplx,nlpspd%nprojel
                  if (istart_c > Lnlpspd%nprojel+1) stop 'istart_c > nprojel+1'
                  if (iproj > Lnlpspd%nproj) stop 'iproj > nproj'
               endif
            enddo
         enddo


!        Apply them on the wavefunctions in the overlap region
!        hpsi contains the new wavefunctions
         call orbs_in_kpt(ikpt,Gorbs,isorb,ieorb,nspinor) 

         do iorb=isorb,ieorb
            do ii=1,orbs%norb
               if (orbs%inWhichLocreg(ii) == iorb) then   !using ii and iorb to identify the orbitals because in linear case, the ordering is different
                                                          !orbitals are now orderer by locreg. So, iorb is the old numbering (i.e. in Global region)
                                                          !while ii is it's numbering in the locreg.

                  istart_o=1
                  do ispinor=1,nspinor,ncplx
                     eproj_spinor = 0.0_gp
                     if (ispinor >= 2) istart_o=1

                     !GTH and HGH pseudopotentials
                     do l=1,4
                        do i=1,3
                           if (atoms%psppar(l,i,ityp) /= 0.0_gp) then
                              call applyprojector(ncplx,l,i,atoms%psppar(0,0,ityp),atoms%npspcode(ityp),&
                                   Llr%wfd%nvctr_c,Llr%wfd%nvctr_f,Llr%wfd%nseg_c,&
                                   Llr%wfd%nseg_f,Llr%wfd%keyv,Llr%wfd%keyg,&
                                   Lnlpspd%nvctr_p(2*iatom-1),Lnlpspd%nvctr_p(2*iatom),Lnlpspd%nseg_p(2*iatom-1),&
                                   Lnlpspd%nseg_p(2*iatom),Lnlpspd%keyv_p(jseg_c),Lnlpspd%keyg_p(1,jseg_c),&
                                   Lproj(istart_o),psi_tmp(1,ispinor,ii),hpsi_tmp(1,ispinor,ii),eproj_spinor)
                               
                               istart_o=istart_o+(Lnlpspd%nvctr_p(2*iatom-1)+7*Lnlpspd%nvctr_p(2*iatom))*(2*l-1)*ncplx
                           end if
                        enddo
                     enddo
                     eproj=eproj+&
                          orbs%kwgts(orbs%iokpt(ii))*orbs%occup(ii+orbs%isorb)*eproj_spinor      
                  end do
               end if
            end do
         end do
         jseg_c = jseg_c + Lnlpspd%nseg_p(2*iatom - 1)+ Lnlpspd%nseg_p(2*iatom) 
      end do  !on iat
      
     ! hpsi = reshape(hpsi_tmp,(/ (Llr%wfd%nvctr_c+7*Llr%wfd%nvctr_f)*orbs%nspinor*LLr%localnorb*nspin/))!,order=(/ 2, 3, 1 /))

      ind = 0
      hpsi = 0.0_wp
      do ispin = 1,nspin                 !is the order correct for spin and spinor?
         do ispinor=1,orbs%nspinor
            do ii=1,orbs%norb/orbs%nspin
               do jj=1,Llr%wfd%nvctr_c+7*Llr%wfd%nvctr_f
                  hpsi(ind+jj) = hpsi_tmp(jj,ispinor,ii+(ispin-1)*orbs%norb/orbs%nspin)
               end do
               ind = ind + Llr%wfd%nvctr_c+7*Llr%wfd%nvctr_f
            end do
         end do
      end do
      if (iproj /= Lnlpspd%nproj) stop 'incorrect number of projectors created'
      if (ieorb == Gorbs%norbp) exit loop_kpt
      ikpt=ikpt+1
   end do loop_kpt

   !deallocate arrays
    i_all = -product(shape(psi_tmp))*kind(psi_tmp)
    deallocate(psi_tmp,stat=i_stat)
    call memocc(i_stat,i_all,'psi_tmp',subname)
    i_all = -product(shape(hpsi_tmp))*kind(hpsi_tmp)
    deallocate(hpsi_tmp,stat=i_stat)
    call memocc(i_stat,i_all,'hpsi_tmp',subname)
    i_all = -product(shape(Lproj))*kind(Lproj)
    deallocate(Lproj,stat=i_stat)
    call memocc(i_stat,i_all,'Lproj',subname)

END SUBROUTINE apply_local_projectors
!%***

!> BigDFT/projector
!!
!!
subroutine local_projector(geocode,atomname,iat,idir,l,i,gau_a,rxyz,Llr,&
     hx,hy,hz,kx,ky,kz,ncplx,&
     mbvctr_c,mbvctr_f,mseg_c,mseg_f,keyv_p,keyg_p,proj,nwarnings)
  use module_base
  use module_types
  implicit none
  character(len=1), intent(in) :: geocode
  character(len=20), intent(in) :: atomname
  type(locreg_descriptors),intent(in) :: Llr
  integer, intent(in) :: iat,idir,l,i,mbvctr_c,mbvctr_f,mseg_c,mseg_f,ncplx
  real(gp), intent(in) :: hx,hy,hz,gau_a,kx,ky,kz
  !integer, dimension(2,3), intent(in) :: nboxp_c,nboxp_f
  integer, dimension(mseg_c+mseg_f), intent(in) :: keyv_p
  integer, dimension(2,mseg_c+mseg_f), intent(in) :: keyg_p
  real(gp), dimension(3), intent(in) :: rxyz
  integer, intent(inout) :: nwarnings
  real(wp), dimension((mbvctr_c+7*mbvctr_f)*(2*l-1)*ncplx), intent(out) :: proj
  !local variables
  integer, parameter :: nterm_max=20 !if GTH nterm_max=4
  integer :: m,iterm
  !integer :: nl1_c,nu1_c,nl2_c,nu2_c,nl3_c,nu3_c,nl1_f,nu1_f,nl2_f,nu2_f,nl3_f,nu3_f
  integer :: istart_c,nterm
  real(gp) :: fpi,factor,rx,ry,rz
  real(dp) :: scpr
  integer, dimension(3) :: nterm_arr
  integer, dimension(nterm_max) :: lx,ly,lz
  integer, dimension(3,nterm_max,3) :: lxyz_arr
  real(gp), dimension(nterm_max) :: factors
  real(gp), dimension(nterm_max,3) :: fac_arr

  !this value can also be inserted as a parameter
  fpi=(4.0_gp*atan(1.0_gp))**(-.75_gp)

  rx=rxyz(1)
  ry=rxyz(2)
  rz=rxyz(3)
  
  istart_c=1
  !start of the projectors expansion routine
  factor=sqrt(2.0_gp)*fpi/(sqrt(gau_a)**(2*(l-1)+4*i-1))
  do m=1,2*l-1

     if (idir==0) then !normal projector calculation case
        call calc_coeff_proj(l,i,m,nterm_max,nterm,lx,ly,lz,factors)

        factors(1:nterm)=factor*factors(1:nterm)
     else !calculation of projector derivative
        call calc_coeff_derproj(l,i,m,nterm_max,gau_a,nterm_arr,lxyz_arr,fac_arr)

        nterm=nterm_arr(idir)
        do iterm=1,nterm
           factors(iterm)=factor*fac_arr(iterm,idir)
           lx(iterm)=lxyz_arr(1,iterm,idir)
           ly(iterm)=lxyz_arr(2,iterm,idir)
           lz(iterm)=lxyz_arr(3,iterm,idir)
        end do
     end if
     
     call crtproj(geocode,nterm,Llr,hx,hy,hz,kx,ky,kz,ncplx,&
          gau_a,factors,rx,ry,rz,lx,ly,lz,&
          mbvctr_c,mbvctr_f,mseg_c,mseg_f,keyv_p,keyg_p,proj(istart_c))

     ! testing
     if (idir == 0) then
        !here the norm should be done with the complex components
        call wnrm_wrap(ncplx,mbvctr_c,mbvctr_f,proj(istart_c),scpr)
        if (abs(1.d0-scpr) > 1.d-2) then
           if (abs(1.d0-scpr) > 1.d-1) then
              !if (iproc == 0) then
                !!write(*,'(1x,a,i4,a,a6,a,i1,a,i1,a,f6.3)')&
                !!      'The norm of the nonlocal PSP for atom n=',iat,&
                !!      ' (',trim(atomname),&
                !!      ') labeled by l=',l,' m=',m,' is ',scpr
                !! write(*,'(1x,a)')&
                !!      'while it is supposed to be about 1.0.'
              !end if
           else
              nwarnings=nwarnings+1
           end if
        end if
     end if
     !end testing
     istart_c=istart_c+(mbvctr_c+7*mbvctr_f)*ncplx
  enddo
END SUBROUTINE local_projector
