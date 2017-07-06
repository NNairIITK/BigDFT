!> @file
!! Datatypes and associated methods relative to the localization regions (mesh grid)
!! @author
!!    Copyright (C) 2007-2015 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS


!> Datatypes for localization regions descriptors
module locregs
  use module_base
  use box, only: cell,box_iterator
  use compression
  use bounds, only: convolutions_bounds
  implicit none

  private 

  !> Grid dimensions in all different wavelet basis
  type, public :: grid_dimensions
     integer :: n1,n2,n3                      !< Coarse grid dimensions
     integer :: nfl1,nfu1,nfl2,nfu2,nfl3,nfu3 !< Lower and upper indices of fine grid in 3D 
     integer :: n1i,n2i,n3i                   !< ISF grid dimension (roughly 2*n+buffer)
  end type grid_dimensions


  !> Contains the information needed for describing completely a wavefunction localisation region
  type, public :: locreg_descriptors
     character(len=1) :: geocode            !< @copydoc poisson_solver::doc::geocode
     logical :: hybrid_on                   !< Interesting for global, periodic, localisation regions
     integer :: ns1,ns2,ns3                 !< Starting point of the localisation region in global coordinates
     integer :: nsi1,nsi2,nsi3              !< Starting point of locreg for interpolating grid
     integer :: Localnorb                   !< Number of orbitals contained in locreg
     integer, dimension(3) :: outofzone     !< Vector of points outside of the zone outside Glr for periodic systems
     real(gp), dimension(3) :: locregCenter !< Center of the locreg 
     real(gp) :: locrad                     !< Cutoff radius of the localization region
     real(gp) :: locrad_kernel              !< Cutoff radius of the localization region (kernel)
     real(gp) :: locrad_mult                !< Cutoff radius of the localization region for the sparse matrix multiplications
     type(grid_dimensions) :: d             !< Grid dimensions in old different wavelet basis
     type(wavefunctions_descriptors) :: wfd
     type(convolutions_bounds) :: bounds
     type(cell) :: mesh !<defines the cell of the system 
                        !! (should replace the other geometrical informations)
     !>iterator over the mesh degrees of freedom
     type(box_iterator) :: bit
  end type locreg_descriptors

  public :: nullify_locreg_descriptors,locreg_null
  public :: deallocate_locreg_descriptors,deallocate_wfd
  public :: allocate_wfd,copy_locreg_descriptors,copy_grid_dimensions
  public :: check_overlap,check_overlap_cubic_periodic,check_overlap_from_descriptors_periodic,lr_box
  public :: init_lr

contains

  pure function grid_null() result(g)
    type(grid_dimensions) :: g
    g%n1   =0
    g%n2   =0
    g%n3   =0
    g%nfl1 =0
    g%nfu1 =0
    g%nfl2 =0
    g%nfu2 =0
    g%nfl3 =0
    g%nfu3 =0
    g%n1i  =0
    g%n2i  =0
    g%n3i  =0
  end function grid_null

  pure function locreg_null() result(lr)
    implicit none
    type(locreg_descriptors) :: lr
    call nullify_locreg_descriptors(lr)
  end function locreg_null

  pure subroutine nullify_locreg_descriptors(lr)
    use box
    use bounds, only: nullify_convolutions_bounds
    implicit none
    type(locreg_descriptors), intent(out) :: lr
    lr%geocode='F'
    lr%hybrid_on=.false.   
    lr%ns1=0
    lr%ns2=0
    lr%ns3=0 
    lr%nsi1=0
    lr%nsi2=0
    lr%nsi3=0  
    lr%Localnorb=0  
    lr%outofzone=(/0,0,0/) 
    lr%d=grid_null()
    call nullify_wfd(lr%wfd)
    call nullify_convolutions_bounds(lr%bounds)
    lr%locregCenter=(/0.0_gp,0.0_gp,0.0_gp/) 
    lr%locrad=0
    lr%mesh=cell_null()
    call nullify_box_iterator(lr%bit)
  end subroutine nullify_locreg_descriptors

  !> Destructors
  subroutine deallocate_locreg_descriptors(lr)
    use bounds
    implicit none
    ! Calling arguments
    type(locreg_descriptors),intent(inout):: lr

    call deallocate_wfd(lr%wfd)
    call deallocate_convolutions_bounds(lr%bounds)
    
  end subroutine deallocate_locreg_descriptors

  !> Methods for copying the structures, can be needed to avoid recalculating them
  !! should be better by defining a f_malloc inheriting the shapes and the structure from other array
  !! of the type dest=f_malloc(src=source,id='dest')
  subroutine copy_locreg_descriptors(glrin, glrout)
    use bounds
    implicit none
    ! Calling arguments
    type(locreg_descriptors), intent(in) :: glrin !<input locreg. Unchanged on exit.
    type(locreg_descriptors), intent(out):: glrout !<output locreg. Must be freed on input.

    glrout%geocode = glrin%geocode
    glrout%hybrid_on = glrin%hybrid_on
    glrout%ns1 = glrin%ns1
    glrout%ns2 = glrin%ns2
    glrout%ns3 = glrin%ns3
    glrout%nsi1 = glrin%nsi1
    glrout%nsi2 = glrin%nsi2
    glrout%nsi3 = glrin%nsi3
    glrout%Localnorb = glrin%Localnorb
    glrout%locrad=glrin%locrad
    glrout%locrad_kernel=glrin%locrad_kernel
    glrout%locrad_mult=glrin%locrad_mult
    glrout%locregCenter=glrin%locregCenter
    glrout%outofzone= glrin%outofzone

    call copy_grid_dimensions(glrin%d, glrout%d)
    call copy_wavefunctions_descriptors(glrin%wfd, glrout%wfd)
    call copy_convolutions_bounds(glrin%bounds, glrout%bounds)

    glrout%mesh=glrin%mesh
    glrout%bit=glrin%bit
  end subroutine copy_locreg_descriptors
  pure subroutine copy_grid_dimensions(din, dout)
    implicit none
    ! Calling arguments
    type(grid_dimensions),intent(in):: din
    type(grid_dimensions),intent(out):: dout

    dout%n1 = din%n1
    dout%n2 = din%n2
    dout%n3 = din%n3
    dout%nfl1 = din%nfl1
    dout%nfu1 = din%nfu1
    dout%nfl2 = din%nfl2
    dout%nfu2 = din%nfu2
    dout%nfl3 = din%nfl3
    dout%nfu3 = din%nfu3
    dout%n1i = din%n1i
    dout%n2i = din%n2i
    dout%n3i = din%n3i

  end subroutine copy_grid_dimensions

  !> Almost degenerate with get_number_of_overlap_region
  !! should merge the two... prefering this one since argument list is better 
  subroutine check_overlap_cubic_periodic(Glr,Ilr,Jlr,isoverlap)
    use module_base
    use bounds, only: check_whether_bounds_overlap
    implicit none
    type(locreg_descriptors), intent(in) :: Glr
    type(locreg_descriptors), intent(in) :: Ilr
    type(locreg_descriptors), intent(in) :: Jlr
    logical, intent(out) :: isoverlap
    !Local variables
    integer :: is1, ie1, is2, ie2, is3, ie3, js1, je1, js2, je2, js3, je3
    logical :: overlap1, overlap2, overlap3
  !!  integer :: azones,bzones,ii,izones,jzones !, i_stat, i_all
  !!  logical :: go1, go2, go3
  !!  integer,dimension(3,8) :: astart,bstart,aend,bend
  
  !!  azones = 1
  !!  bzones = 1
  !!! Calculate the number of regions to cut alr and blr
  !!  do ii=1,3
  !!     if(Ilr%outofzone(ii) > 0) azones = azones * 2
  !!     if(Jlr%outofzone(ii) > 0) bzones = bzones * 2
  !!  end do
  !!
  !!!FRACTURE THE FIRST LOCALIZATION REGION
  !!  call fracture_periodic_zone(azones,Glr,Ilr,Ilr%outofzone,astart,aend)
  !!
  !!!FRACTURE SECOND LOCREG
  !!  call fracture_periodic_zone(bzones,Glr,Jlr,Jlr%outofzone,bstart,bend)
  !!
  !!! Now check if they overlap
  !!  isoverlap = .false.
  !!  loop_izones: do izones=1,azones
  !!    do jzones=1,bzones
  !!      go1 = (bstart(1,jzones) .le. aend(1,izones) .and. bend(1,jzones) .ge. astart(1,izones))
  !!      go2 = (bstart(2,jzones) .le. aend(2,izones) .and. bend(2,jzones) .ge. astart(2,izones))
  !!      go3 = (bstart(3,jzones) .le. aend(3,izones) .and. bend(3,jzones) .ge. astart(3,izones))
  !!      if(go1 .and. go2 .and. go3) then
  !!        isoverlap = .true.
  !!        exit loop_izones
  !!      end if
  !!    end do
  !!  end do loop_izones
  
  
    !@ NEW VERSION #########################################
    ! Shift all the indices into the periodic cell. This can result is starting
    ! indices being larger than ending indices
    isoverlap = .false.
    is3 = modulo(ilr%ns3,glr%d%n3+1)
    ie3 = modulo(ilr%ns3+ilr%d%n3,glr%d%n3+1)
    js3 = modulo(jlr%ns3,glr%d%n3+1)
    je3 = modulo(jlr%ns3+jlr%d%n3,glr%d%n3+1)
    overlap3 = check_whether_bounds_overlap(is3, ie3, js3, je3)
    if (overlap3) then
        is2 = modulo(ilr%ns2,glr%d%n2+1)
        ie2 = modulo(ilr%ns2+ilr%d%n2,glr%d%n2+1)
        js2 = modulo(jlr%ns2,glr%d%n2+1)
        je2 = modulo(jlr%ns2+jlr%d%n2,glr%d%n2+1)
        overlap2 = check_whether_bounds_overlap(is2, ie2, js2, je2)
        if (overlap2) then
            is1 = modulo(ilr%ns1,glr%d%n1+1)
            ie1 = modulo(ilr%ns1+ilr%d%n1,glr%d%n1+1)
            js1 = modulo(jlr%ns1,glr%d%n1+1)
            je1 = modulo(jlr%ns1+jlr%d%n1,glr%d%n1+1)
            overlap1 = check_whether_bounds_overlap(is1, ie1, js1, je1)
            if (overlap1) then
                ! If we are here, all three overlaps are true
                isoverlap = .true.
            end if
        end if
    end if
  
    !!if (overlap1 .and. overlap2 .and. overlap3) then
    !!    isoverlap = .true.
    !!else
    !!    isoverlap = .false.
    !!end if
        
    !@ END NEW VERSION #####################################
  
    !!!debug
    !!isoverlap=.true.
  
  end subroutine check_overlap_cubic_periodic

    subroutine check_overlap(Llr_i, Llr_j, Glr, overlap)
      implicit none

      ! Calling arguments
      type(locreg_descriptors),intent(in) :: Llr_i, Llr_j, Glr
      logical, intent(out) :: overlap

      ! Local variables
      integer :: onseg

      call check_overlap_cubic_periodic(Glr,Llr_i,Llr_j,overlap)
      if(overlap) then
         call check_overlap_from_descriptors_periodic(Llr_i%wfd%nseg_c, Llr_j%wfd%nseg_c,&
              Llr_i%wfd%keyglob, Llr_j%wfd%keyglob, overlap, onseg)
      end if

    end subroutine check_overlap

    ! check if Llrs overlap from there descriptors
    ! The periodicity is hidden in the fact that we are using the keyglobs
    ! which are correctly defined. 
    subroutine check_overlap_from_descriptors_periodic(nseg_i, nseg_j, keyg_i, keyg_j,  &
         isoverlap, onseg)
      implicit none
      ! Calling arguments
      integer :: nseg_i, nseg_j
      integer,dimension(2,nseg_i),intent(in) :: keyg_i
      integer,dimension(2,nseg_j),intent(in) :: keyg_j
      logical,intent(out) :: isoverlap
      integer, intent(out) :: onseg
      ! Local variables
      integer :: iseg, jseg, istart, jstart, kstartg
      integer :: iend, jend, kendg, nseg_k


      ! Initialize some counters
      iseg=1
      jseg=1
      nseg_k=0
      isoverlap = .false.
      onseg = 0  ! in case they don't overlap
      ! Check whether all segments of both localization regions have been processed.
      if ((iseg>=nseg_i .and. jseg>=nseg_j) .or. nseg_i==0 .or. nseg_j==0) return

      segment_loop: do

         ! Starting point already in global coordinates
         istart=keyg_i(1,iseg)
         jstart=keyg_j(1,jseg)

         ! Ending point already in global coordinates
         iend=keyg_i(2,iseg)
         jend=keyg_j(2,jseg)
         ! Determine starting and ending point of the common segment in global coordinates.
         kstartg=max(istart,jstart)
         kendg=min(iend,jend)

         ! Check whether this common segment has a non-zero length
         if(kendg-kstartg+1>0) then
            isoverlap = .true.
            nseg_k=nseg_k+1
         end if

         ! Check whether all segments of both localization regions have been processed.
         if(iseg>=nseg_i .and. jseg>=nseg_j) exit segment_loop

         ! Increase the segment index
         if((iend<=jend .and. iseg<nseg_i) .or. jseg==nseg_j) then
            iseg=iseg+1
         else if(jseg<nseg_j) then
            jseg=jseg+1
         end if

      end do segment_loop

      if(isoverlap) then
         onseg = nseg_k
      end if

    end subroutine check_overlap_from_descriptors_periodic

    pure function grid_init(peri,n1,n2,n3,nfl1,nfl2,nfl3,nfu1,nfu2,nfu3,&
         ns1,ns2,ns3) result(g)
      implicit none
      integer, intent(in) :: n1,n2,n3
      integer, intent(in) :: nfl1,nfl2,nfl3,nfu1,nfu2,nfu3
      integer, intent(in) :: ns1,ns2,ns3
      logical, dimension(3), intent(in) :: peri !<periodic dimensions
      type(grid_dimensions) :: g
      !local variables
      integer, parameter :: ISF_GROW_BUFFER=31
      integer :: isx,isy,isz
      
      g%n1=n1-ns1
      g%n2=n2-ns2
      g%n3=n3-ns3

      !dimensions of the fine grid inside the localisation region
      g%nfl1=max(ns1,nfl1)-ns1
      g%nfl2=max(ns2,nfl2)-ns2
      g%nfl3=max(ns3,nfl3)-ns3

      !NOTE: This will not work with symmetries (must change it)
      g%nfu1=min(n1,nfu1)-ns1
      g%nfu2=min(n2,nfu2)-ns2
      g%nfu3=min(n3,nfu3)-ns3

      if (peri(1)) then
         g%n1i=2*g%n1+2
      else
         g%n1i=2*g%n1+ISF_GROW_BUFFER
      end if
      if (peri(2)) then
         g%n2i=2*g%n2+2
      else
         g%n2i=2*g%n2+ISF_GROW_BUFFER
      end if
      if (peri(3)) then
         g%n3i=2*g%n3+2
      else
         g%n3i=2*g%n3+ISF_GROW_BUFFER
      end if

    end function grid_init

    !> Create the localisation region information for cubic code
    subroutine init_lr(lr,geocode,hgridsh,n1,n2,n3,nfl1,nfl2,nfl3,nfu1,nfu2,nfu3,&
         hybrid_flag,isx,isy,isz,global_geocode,wfd,bnds)
      use compression
      use bounds
      use box
      implicit none
      logical, intent(in) :: hybrid_flag
      character(len=1), intent(in) :: geocode !< @copydoc poisson_solver::doc::geocode
      integer, intent(in) :: n1,n2,n3,nfl1,nfl2,nfl3,nfu1,nfu2,nfu3
      real(gp), dimension(3), intent(in) :: hgridsh
      character(len=1), intent(in), optional :: global_geocode
      !>have to be present with global_geocode
      integer, intent(in), optional :: isx,isy,isz 
      type(wavefunctions_descriptors), intent(in), optional :: wfd
      type(convolutions_bounds), intent(in), optional :: bnds
      type(locreg_descriptors), intent(inout) :: lr
      !local variables
      integer, parameter :: S0_GROW_BUFFER=14
      integer :: Gnbl1,Gnbl2,Gnbl3,Gnbr1,Gnbr2,Gnbr3
      integer :: Lnbl1,Lnbl2,Lnbl3,Lnbr1,Lnbr2,Lnbr3
      logical, dimension(3) :: peri,peri_glob

      lr%geocode=geocode
      lr%ns1=0
      lr%ns2=0
      lr%ns3=0
      if (present(isx)) lr%ns1=isx
      if (present(isy)) lr%ns2=isy
      if (present(isz)) lr%ns3=isz

      peri(1)=geocode /= 'F'
      peri(2)=geocode == 'P'
      peri(3)=geocode /= 'F'

      lr%d=grid_init(peri,n1,n2,n3,nfl1,nfl2,nfl3,nfu1,nfu2,nfu3,&
         lr%ns1,lr%ns2,lr%ns3)

      lr%mesh=cell_new(geocode,[lr%d%n1i,lr%d%n2i,lr%d%n3i],hgridsh)

      Gnbl1=0
      Gnbl2=0
      Gnbl3=0
      if (present(global_geocode)) then
         peri_glob(1)=(global_geocode /= 'F')
         peri_glob(2)=(global_geocode == 'P')
         peri_glob(3)=(global_geocode /= 'F')
         call ext_buffers(peri_glob(1),Gnbl1,Gnbr1)
         call ext_buffers(peri_glob(2),Gnbl2,Gnbr2)
         call ext_buffers(peri_glob(3),Gnbl3,Gnbr3)
         call ext_buffers(peri(1),Lnbl1,Lnbr1)
         call ext_buffers(peri(2),Lnbl2,Lnbr2)
         call ext_buffers(peri(3),Lnbl3,Lnbr3)
      end if
      lr%nsi1=0
      lr%nsi2=0
      lr%nsi3=0
      if (present(isx)) lr%nsi1= 2 * lr%ns1 - (Lnbl1 - Gnbl1)
      if (present(isy)) lr%nsi2= 2 * lr%ns2 - (Lnbl2 - Gnbl2)
      if (present(isz)) lr%nsi3= 2 * lr%ns3 - (Lnbl3 - Gnbl3)
      

      lr%hybrid_on = hybrid_flag
      lr%hybrid_on=lr%hybrid_on .and. (nfu1-nfl1+S0_GROW_BUFFER < n1+1)
      lr%hybrid_on=lr%hybrid_on .and. (nfu2-nfl2+S0_GROW_BUFFER < n2+1)
      lr%hybrid_on=lr%hybrid_on .and. (nfu3-nfl3+S0_GROW_BUFFER < n3+1)

      if (present(wfd)) lr%wfd=wfd !it just associates the pointers
      if (geocode == 'F' .and. present(bnds)) lr%bounds=bnds

      lr%bit=box_iter(lr%mesh,origin=locreg_mesh_origin(lr%mesh))

    END SUBROUTINE init_lr

    !> initalize the box-related components of the localization regions
    subroutine lr_box(lr,Glr,hgrids,nbox,correct)
      use bounds, only: ext_buffers
      implicit none
      logical, intent(in) :: correct
      !> Sub-box to iterate over the points (ex. around atoms)
      !! start and end points for each direction
      integer, dimension(2,3), intent(in) :: nbox
      real(gp), dimension(3), intent(in) :: hgrids
      type(locreg_descriptors), intent(in) :: Glr
      type(locreg_descriptors), intent(inout) :: lr
      !local variables
      character(len=1) :: geocode
      logical :: Gperx,Gpery,Gperz,xperiodic,yperiodic,zperiodic
      integer :: isx,iex,isy,iey,isz,iez
      integer :: Gnbl1,Gnbl2,Gnbl3,Gnbr1,Gnbr2,Gnbr3
      integer :: Lnbl1,Lnbl2,Lnbl3,Lnbr1,Lnbr2,Lnbr3
      integer :: ln1,ln2,ln3
      logical, dimension(3) :: peri
      integer, dimension(3) :: outofzone

      !initialize out of zone
      outofzone (:) = 0

      ! Localization regions should have free boundary conditions by default
      isx=nbox(1,1)
      iex=nbox(2,1)
      isy=nbox(1,2)
      iey=nbox(2,2)
      isz=nbox(1,3)
      iez=nbox(2,3)

      ln1 = iex-isx
      ln2 = iey-isy
      ln3 = iez-isz

      geocode='F'

      xperiodic = .false.
      yperiodic = .false.
      zperiodic = .false. 

      !assign the starting/ending points and outofzone for the different
      ! geometries
      select case(Glr%geocode)
      case('F')
         isx=max(isx,Glr%ns1)
         isy=max(isy,Glr%ns2)
         isz=max(isz,Glr%ns3)

         iex=min(iex,Glr%ns1+Glr%d%n1)
         iey=min(iey,Glr%ns2+Glr%d%n2)
         iez=min(iez,Glr%ns3+Glr%d%n3)

      case('S')
         ! Get starting and ending for x direction     
         if (iex - isx >= Glr%d%n1) then       
            isx=Glr%ns1
            iex=Glr%ns1 + Glr%d%n1
            xperiodic = .true.
         else
            if (correct) then
               isx=modulo(isx,Glr%d%n1+1) + Glr%ns1
               iex= ln1 + isx
            end if
            if (iex > Glr%ns1+Glr%d%n1) then
               outofzone(1)=modulo(iex,Glr%d%n1+1)
            end if
         end if

         ! Get starting and ending for y direction (perpendicular to surface)
         isy=max(isy,Glr%ns2)
         iey=min(iey,Glr%ns2 + Glr%d%n2)
         outofzone(2) = 0

         !Get starting and ending for z direction
         if (iez - isz >= Glr%d%n3) then
            isz=Glr%ns3 
            iez=Glr%ns3 + Glr%d%n3
            zperiodic = .true.
         else
            if (correct) then
               isz=modulo(isz,Glr%d%n3+1) +  Glr%ns3
               iez= ln3 + isz
            end if
            if (iez > Glr%ns3+Glr%d%n3) then
               outofzone(3)=modulo(iez,Glr%d%n3+1)
            end if
         end if
         if(xperiodic .and. zperiodic) then
            geocode = 'S'
         end if

      case('P')
         ! Get starting and ending for x direction     
         if (iex - isx >= Glr%d%n1) then       
            isx=Glr%ns1
            iex=Glr%ns1 + Glr%d%n1
            xperiodic = .true.
         else
            if (correct) then
               isx=modulo(isx,Glr%d%n1+1) + Glr%ns1
               iex= ln1 + isx
            end if
            if (iex > Glr%ns1+Glr%d%n1) then
               outofzone(1)=modulo(iex,Glr%d%n1+1)
            end if
         end if

         ! Get starting and ending for y direction (perpendicular to surface)
         if (iey - isy >= Glr%d%n2) then       
            isy=Glr%ns2
            iey=Glr%ns2 + Glr%d%n2
            yperiodic = .true.
         else
            if (correct) then
               isy=modulo(isy,Glr%d%n2+1) + Glr%ns2
               iey= ln2 + isy
            end if
            if (iey > Glr%ns2+Glr%d%n2) then
               outofzone(2)=modulo(iey,Glr%d%n2+1)
            end if
         end if

         !Get starting and ending for z direction
         if (iez - isz >= Glr%d%n3) then
            isz=Glr%ns3 
            iez=Glr%ns3 + Glr%d%n3
            zperiodic = .true.
         else
            if (correct) then
               isz=modulo(isz,Glr%d%n3+1) +  Glr%ns3
               iez= ln3 + isz
            end if
            if (iez > Glr%ns3+Glr%d%n3) then
               outofzone(3)=modulo(iez,Glr%d%n3+1)
            end if
         end if
         if(xperiodic .and. yperiodic .and. zperiodic ) then
            geocode = 'P'
         end if
      end select

      ! Make sure that the localization regions are not periodic
      if (xperiodic .or. yperiodic .or. zperiodic) then
         call f_err_throw('The size of the localization region '&
              &//&!trim(yaml_toa(ilr,fmt='(i0)'))//&
              &' is larger than that of the global region.&
              & Reduce the localization radii or use the cubic version',&
              & err_name='BIGDFT_RUNTIME_ERROR')
      end if

      call init_lr(lr,geocode,0.5_gp*hgrids,iex,iey,iez,&
           Glr%d%nfl1,Glr%d%nfl2,Glr%d%nfl3,&
           Glr%d%nfu1,Glr%d%nfu2,Glr%d%nfu3,&
           .false.,isx,isy,isz,Glr%geocode)

      !assign outofzone
      lr%outofzone(:) = outofzone(:)

!!$      !values for the starting point of the cube for wavelet grid
!!$      lr%ns1=isx
!!$      lr%ns2=isy
!!$      lr%ns3=isz
!!$
!!$      ! Set the conditions for ext_buffers (conditions for buffer size)
!!$      Gperx=(Glr%geocode /= 'F')
!!$      Gpery=(Glr%geocode == 'P')
!!$      Gperz=(Glr%geocode /= 'F')
!!$      peri(1)=(lr%geocode /= 'F')
!!$      peri(2)=(lr%geocode == 'P')
!!$      peri(3)=(lr%geocode /= 'F')
!!$
!!$      !calculate the size of the buffers of interpolating function grid
!!$      call ext_buffers(Gperx,Gnbl1,Gnbr1)
!!$      call ext_buffers(Gpery,Gnbl2,Gnbr2)
!!$      call ext_buffers(Gperz,Gnbl3,Gnbr3)
!!$      call ext_buffers(peri(1),Lnbl1,Lnbr1)
!!$      call ext_buffers(peri(2),Lnbl2,Lnbr2)
!!$      call ext_buffers(peri(3),Lnbl3,Lnbr3)
!!$
!!$      !starting point of the region for interpolating functions grid
!!$      lr%nsi1= 2 * lr%ns1 - (Lnbl1 - Gnbl1)
!!$      lr%nsi2= 2 * lr%ns2 - (Lnbl2 - Gnbl2)
!!$      lr%nsi3= 2 * lr%ns3 - (Lnbl3 - Gnbl3)
!!$      !write(*,*) 'ilr, lr%nsi3',ilr, lr%nsi3
!!$
!!$      lr%d=grid_init(peri,iex,iey,iez,&
!!$           Glr%d%nfl1,Glr%d%nfl2,Glr%d%nfl3,&
!!$           Glr%d%nfu1,Glr%d%nfu2,Glr%d%nfu3,&
!!$           isx,isy,isz)

!!$      !dimensions of the localisation region
!!$      lr%d%n1=iex-isx
!!$      lr%d%n2=iey-isy
!!$      lr%d%n3=iez-isz
!!$
!!$      !dimensions of the fine grid inside the localisation region
!!$      lr%d%nfl1=max(isx,Glr%d%nfl1)-isx ! should we really substract isx (probably because the routines are coded with 0 as origin)?
!!$      lr%d%nfl2=max(isy,Glr%d%nfl2)-isy
!!$      lr%d%nfl3=max(isz,Glr%d%nfl3)-isz
!!$
!!$      !NOTE: This will not work with symmetries (must change it)
!!$      lr%d%nfu1=min(iex,Glr%d%nfu1)-isx
!!$      lr%d%nfu2=min(iey,Glr%d%nfu2)-isy
!!$      lr%d%nfu3=min(iez,Glr%d%nfu3)-isz
!!$
!!$      !dimensions of the interpolating scaling functions grid (reduce to +2 for periodic)
!!$      if(lr%geocode == 'F') then
!!$         lr%d%n1i=2*lr%d%n1+31
!!$         lr%d%n2i=2*lr%d%n2+31
!!$         lr%d%n3i=2*lr%d%n3+31
!!$      else if(lr%geocode == 'S') then
!!$         lr%d%n1i=2*lr%d%n1+2
!!$         lr%d%n2i=2*lr%d%n2+31
!!$         lr%d%n3i=2*lr%d%n3+2
!!$      else
!!$         lr%d%n1i=2*lr%d%n1+2
!!$         lr%d%n2i=2*lr%d%n2+2
!!$         lr%d%n3i=2*lr%d%n3+2
!!$      end if

      ! Make sure that the extent of the interpolating functions grid for the
      ! locreg is not larger than the that of the global box.
      if (lr%d%n1i>Glr%d%n1i) then
         call f_err_throw('The interpolating functions grid in x dimension for locreg '&
              &//&!trim(yaml_toa(ilr,fmt='(i0)'))//&
              '('//trim(yaml_toa(lr%d%n1i,fmt='(i0)'))//')&
              & is larger than that of the global region('//trim(yaml_toa(Glr%d%n1i,fmt='(i0)'))//').&
              & Reduce the localization radii or use the cubic version',&
              & err_name='BIGDFT_RUNTIME_ERROR')
      end if
      if (lr%d%n2i>Glr%d%n2i) then
         call f_err_throw('The interpolating functions grid in y dimension for locreg '&
              !&//trim(yaml_toa(ilr,fmt='(i0)'))&
              //'('//trim(yaml_toa(lr%d%n2i,fmt='(i0)'))//')&
              & is larger than that of the global region('//trim(yaml_toa(Glr%d%n2i,fmt='(i0)'))//').&
              & Reduce the localization radii or use the cubic version',&
              & err_name='BIGDFT_RUNTIME_ERROR')
      end if
      if (lr%d%n3i>Glr%d%n3i) then
         call f_err_throw('The interpolating functions grid in z dimension for locreg '&
              !&//trim(yaml_toa(ilr,fmt='(i0)'))&
              //'('//trim(yaml_toa(lr%d%n3i,fmt='(i0)'))//')&
              & is larger than that of the global region('//trim(yaml_toa(Glr%d%n3i,fmt='(i0)'))//').&
              & Reduce the localization radii or use the cubic version',&
              & err_name='BIGDFT_RUNTIME_ERROR')
      end if
      
    end subroutine lr_box


end module locregs
