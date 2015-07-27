module locreg_operations
  use module_base
  use locregs
  implicit none

  private

  public :: Lpsi_to_global2
  public :: global_to_local_parallel
  public :: get_boundary_weight
  public :: small_to_large_locreg
  public :: psi_to_locreg2
  public :: global_to_local

  contains

    !> Tranform wavefunction between localisation region and the global region
    !!!!!#######!> This routine only works if both locregs have free boundary conditions.
    !! @warning 
    !! WARNING: Make sure psi is set to zero where Glr does not collide with Llr (or everywhere)
    subroutine Lpsi_to_global2(iproc, ldim, gdim, norb, nspinor, nspin, Glr, Llr, lpsi, psi)
    
      use module_base
    
     implicit none
    
      ! Subroutine Scalar Arguments
      integer,intent(in):: iproc
      integer :: Gdim          ! dimension of psi 
      integer :: Ldim          ! dimension of lpsi
      integer :: norb          ! number of orbitals
      integer :: nspinor       ! number of spinors
      integer :: nspin         ! number of spins 
      type(locreg_descriptors),intent(in) :: Glr  ! Global grid descriptor
      type(locreg_descriptors), intent(in) :: Llr  ! Localization grid descriptors 
      
      !Subroutine Array Arguments
      real(wp),dimension(Gdim),intent(inout) :: psi       !Wavefunction (compressed format)
      real(wp),dimension(Ldim),intent(in) :: lpsi         !Wavefunction in localization region
      
      !local variables
      integer :: igrid,isegloc,isegG,ix!,iorbs
      integer :: lmin,lmax,Gmin,Gmax
      integer :: icheck      ! check to make sure the dimension of loc_psi does not overflow 
      integer :: offset      ! gives the difference between the starting point of Lseg and Gseg
      integer :: length      ! Length of the overlap between Lseg and Gseg
      integer :: lincrement  ! Increment for writing orbitals in loc_psi
      integer :: Gincrement  ! Increment for reading orbitals in psi
      integer :: nseg        ! total number of segments in Llr
      !integer, allocatable :: keymask(:,:)  ! shift for every segment of Llr (with respect to Glr)
      character(len=*), parameter :: subname='Lpsi_to_global'
      integer :: i_all
      integer :: start,Gstart,Lindex
      integer :: lfinc,Gfinc,spinshift,ispin,Gindex,isegstart
      integer :: istart
      !integer :: i_stat
    
      call f_routine(id=subname)
    
      !!! This routine is only intended for conversions between locregs with the same boundary conditions.
      !!if (glr%geocode/= 'F' .or. llr%geocode/='F') then
      !!    call f_err_throw('Lpsi_to_global2 can only be used for locregs with free boundary conditions', &
      !!         err_name='BIGDFT_RUNTIME_ERROR')
      !!end if
    
      if(nspin/=1) stop 'not fully implemented for nspin/=1!'
    
    ! Define integers
      nseg = Llr%wfd%nseg_c + Llr%wfd%nseg_f
      lincrement = Llr%wfd%nvctr_c + 7*Llr%wfd%nvctr_f
      Gincrement = Glr%wfd%nvctr_c + 7*Glr%wfd%nvctr_f
      icheck = 0
      spinshift = Gdim / nspin
     
    ! Get the keymask: shift for every segment of Llr (with respect to Glr)
    ! allocate(keymask(2,nseg),stat=i_stat)
      !keymask = f_malloc((/2,nseg/),id='keymask')
    
      !call shift_locreg_indexes(Glr,Llr,keymask,nseg)
      !call shift_locreg_indexes_global(Glr,Llr,keymask,nseg)
      !!keymask = llr%wfd%keyglob
    
    !####################################################
    ! Do coarse region
    !####################################################
      isegstart=1
    
     
      !$omp parallel default(private) &
      !$omp shared(Glr,Llr, lpsi,icheck,psi,norb) &
      !$omp firstprivate(isegstart,nseg,lincrement,Gincrement,spinshift,nspin) 
    
      !$omp do reduction(+:icheck)
      local_loop_c: do isegloc = 1,Llr%wfd%nseg_c
         lmin = llr%wfd%keyglob(1,isegloc)
         lmax = llr%wfd%keyglob(2,isegloc)
         istart = llr%wfd%keyvglob(isegloc)-1
    
         
         global_loop_c: do isegG = isegstart,Glr%wfd%nseg_c
            Gmin = Glr%wfd%keyglob(1,isegG)
            Gmax = Glr%wfd%keyglob(2,isegG)
    
            ! For each segment in Llr check if there is a collision with the segment in Glr
            !if not, cycle
            if(lmin > Gmax) then
                isegstart=isegG
            end if
            if(Gmin > lmax) exit global_loop_c
    
            !if((lmin > Gmax) .or. (lmax < Gmin))  cycle global_loop_c
            if(lmin > Gmax)  cycle global_loop_c
    
            ! Define the offset between the two segments
            offset = lmin - Gmin
            if(offset < 0) then
               offset = 0
            end if
    
            ! Define the length of the two segments
            length = min(lmax,Gmax)-max(lmin,Gmin)
    
            !Find the common elements and write them to the new global wavefunction
            icheck = icheck + (length + 1)
    
            ! WARNING: index goes from 0 to length because it is the offset of the element
    
            do ix = 0,length     
               istart = istart + 1
               do ispin=1,nspin
                  Gindex = Glr%wfd%keyvglob(isegG)+offset+ix+spinshift*(ispin-1)
                  Lindex = istart+lincrement*norb*(ispin-1)
                  psi(Gindex) = lpsi(Lindex) 
               end do
            end do
         end do global_loop_c
      end do local_loop_c
      !$omp end do
    
    
    !##############################################################
    ! Now do fine region
    !##############################################################
    
      start = Llr%wfd%nvctr_c
      Gstart = Glr%wfd%nvctr_c
      lfinc  = Llr%wfd%nvctr_f
      Gfinc = Glr%wfd%nvctr_f
    
      isegstart=Glr%wfd%nseg_c+1
    
      !$omp do reduction(+:icheck)
      local_loop_f: do isegloc = Llr%wfd%nseg_c+1,nseg
         lmin = llr%wfd%keyglob(1,isegloc)
         lmax = llr%wfd%keyglob(2,isegloc)
         istart = llr%wfd%keyvglob(isegloc)-1
    
         global_loop_f: do isegG = isegstart,Glr%wfd%nseg_c+Glr%wfd%nseg_f
    
            Gmin = Glr%wfd%keyglob(1,isegG)
            Gmax = Glr%wfd%keyglob(2,isegG)
    
            ! For each segment in Llr check if there is a collision with the segment in Glr
            ! if not, cycle
            if(lmin > Gmax) then
                isegstart=isegG
            end if
            if(Gmin > lmax)  exit global_loop_f
            !if((lmin > Gmax) .or. (lmax < Gmin))  cycle global_loop_f
            if(lmin > Gmax)  cycle global_loop_f
    
            offset = lmin - Gmin
            if(offset < 0) offset = 0
    
            length = min(lmax,Gmax)-max(lmin,Gmin)
    
            !Find the common elements and write them to the new global wavefunction
            ! First set to zero those elements which are not copied. WARNING: will not work for npsin>1!!
     
            icheck = icheck + (length + 1)
    
            ! WARNING: index goes from 0 to length because it is the offset of the element
            do ix = 0,length
            istart = istart + 1
               do igrid=1,7
                  do ispin = 1, nspin
                     Gindex = Gstart + (Glr%wfd%keyvglob(isegG)+offset+ix-1)*7+igrid + spinshift*(ispin-1)
                     Lindex = start+(istart-1)*7+igrid + lincrement*norb*(ispin-1) 
                     psi(Gindex) = lpsi(Lindex) 
                  end do
               end do
            end do
         end do global_loop_f
      end do local_loop_f
      !$omp end do
    
      !$omp end parallel
    
      !Check if the number of elements in loc_psi is valid
      if(icheck .ne. Llr%wfd%nvctr_f+Llr%wfd%nvctr_c) then
        write(*,*)'There is an error in Lpsi_to_global2: sum of fine and coarse points used',icheck
        write(*,*)'is not equal to the sum of fine and coarse points in the region',Llr%wfd%nvctr_f+Llr%wfd%nvctr_c
        stop
      end if
    
      !!call f_free(keymask)
    
      call f_release_routine()
    
    END SUBROUTINE Lpsi_to_global2


    !> Projects a quantity stored with the global indexes (i1,i2,i3) within the localisation region.
    !! @warning       
    !!    The quantity must not be stored in a compressed form.
    subroutine global_to_local_parallel(Glr,Llr,size_rho,size_Lrho,rho,Lrho,i1s,i1e,i2s,i2e,i3s,i3e,ni1,ni2, &
               i1shift, i2shift, i3shift, ise)
    
     use module_base
     
     implicit none
    
     ! Arguments
     type(locreg_descriptors),intent(in) :: Llr   ! Local localization region
     type(locreg_descriptors),intent(in) :: Glr   ! Global localization region
     integer, intent(in) :: size_rho  ! size of rho array
     integer, intent(in) :: size_Lrho ! size of Lrho array
     real(wp),dimension(size_rho),intent(in) :: rho  ! quantity in global region
     real(wp),dimension(size_Lrho),intent(out) :: Lrho ! piece of quantity in local region
     integer,intent(in) :: i1s, i1e, i2s, i2e
     integer,intent(in) :: i3s, i3e ! starting and ending indices on z direction (related to distribution of rho when parallel)
     integer,intent(in) :: ni1, ni2 ! x and y extent of rho
     integer,intent(in) :: i1shift, i2shift, i3shift
     integer,dimension(6) :: ise
    
    ! Local variable
     integer :: ispin,i1,i2,i3,ii1,ii2,ii3  !integer for loops
     integer :: indSmall, indSpin, indLarge ! indexes for the arrays
     integer :: ist2S,ist3S, ist2L, ist3L, istsa, ists, istl
     integer :: ii1shift, ii2shift, ii3shift, i1glob, i2glob, i3glob
     integer :: iii1, iii2, iii3
    
     !THIS ROUTINE NEEDS OPTIMIZING
    
     !!write(*,'(a,8i8)') 'in global_to_local_parallel: i1s, i1e, i2s, i2e, i3s, i3e, ni1, ni2', i1s, i1e, i2s, i2e, i3s, i3e, ni1, ni2
     
     ! Cut out a piece of the quantity (rho) from the global region (rho) and
     ! store it in a local region (Lrho).
     indSmall=0
     indSpin=0
     ! Deactivate the spin for the moment
     do ispin=1,1!nspin
         !$omp parallel default(none) &
         !$omp shared(Glr, Llr, Lrho, rho, indSpin, i1s, i1e, i2s, i2e, i3s, i3e) &
         !$omp shared(i1shift, i2shift, i3shift, ni1, ni2, ise) &
         !$omp private(ii1, ii2, ii3, i1glob, i2glob, i3glob, ii1shift, ii2shift, ii3shift) &
         !$omp private(ist3S, ist3L, istsa, ist2S, ist2L, ists, istl, indSmall, indLarge) &
         !$omp private(iii1, iii2, iii3)
         !$omp do
         do ii3=i3s,i3e
             i3glob = ii3+ise(5)-1
             !i3=modulo(i3glob-1,glr%d%n3i)+1
             if (modulo(ii3-1,glr%d%n3i)+1>modulo(i3e-1,glr%d%n3i)+1) then
                 !This is a line before the wrap around, i.e. one needs a shift since 
                 ii3shift = i3shift
             else
                 ii3shift = 0
             end if
             if (i3glob<=glr%d%n3i) then
                 iii3=ii3+i3shift
             else
                 iii3=modulo(i3glob-1,glr%d%n3i)+1
             end if
             ist3S = (ii3-i3s)*Llr%d%n2i*Llr%d%n1i
             ist3L = (iii3-1)*ni2*ni1
             istsa=ist3S-i1s+1
             do ii2=i2s,i2e
                 i2glob = ii2+ise(3)-1
                 !i2=modulo(i2glob-1,glr%d%n2i)+1
                 if (modulo(ii2-1,glr%d%n2i)+1>modulo(i2e-1,glr%d%n2i)+1) then
                     !This is a line before the wrap around, i.e. one needs a shift since 
                     !the potential in the global region starts with the wrapped around part
                     ii2shift = i2shift
                 else
                     ii2shift = 0
                 end if
                 if (i2glob<=glr%d%n2i) then
                     iii2=ii2+i2shift
                 else
                     iii2=modulo(i2glob-1,glr%d%n2i)+1
                 end if
                 ist2S = (ii2-i2s)*Llr%d%n1i 
                 ist2L = (iii2-1)*ni1
                 ists=istsa+ist2S
                 istl=ist3L+ist2L
                 do ii1=i1s,i1e
                     i1glob = ii1+ise(1)-1
                     !i1=modulo(i1glob-1,glr%d%n1i)+1
                     if (modulo(ii1-1,glr%d%n1i)+1>modulo(i1e-1,glr%d%n1i)+1) then
                         !This is a line before the wrap around, i.e. one needs a shift since 
                         !the potential in the global region starts with the wrapped around part
                         ii1shift = i1shift
                     else
                         ii1shift = 0
                     end if
                     if (i1glob<=glr%d%n1i) then
                         iii1=ii1+i1shift
                     else
                         iii1=modulo(i1glob-1,glr%d%n1i)+1
                     end if
                     ! indSmall is the index in the local localization region
                     indSmall=ists+ii1
                     ! indLarge is the index in the global localization region. 
                     indLarge= iii1+istl
                     Lrho(indSmall)=rho(indLarge+indSpin)
                     !write(600+bigdft_mpi%iproc,'(a,14i7,2es16.8)') 'i1glob, i2glob, i3glob, i1, i2, i3, iii1, iii2, iii3, i1shift, i2shift, i3shift, indsmall, indlarge, val, testval', &
                     !    i1glob, i2glob, i3glob, i1, i2, i3, iii1, iii2, iii3, i1shift, i2shift, i3shift, indsmall, indlarge, Lrho(indSmall), real((i1+(i2-1)*glr%d%n1i+(i3-1)*glr%d%n1i*glr%d%n2i),kind=8)
                     !if (abs(Lrho(indSmall)-real((i1+(i2-1)*glr%d%n1i+(i3-1)*glr%d%n1i*glr%d%n2i),kind=8))>1.d-3) then
                     !    write(700+bigdft_mpi%iproc,'(a,11i7,2es16.8)') 'i1glob, i2glob, i3glob, i1, i2, i3, iii1, iii2, iii3, indsmall, indlarge, val, testval', &
                     !        i1glob, i2glob, i3glob, i1, i2, i3, iii1, iii2, iii3, indsmall, indlarge, Lrho(indSmall), real((i1+(i2-1)*glr%d%n1i+(i3-1)*glr%d%n1i*glr%d%n2i),kind=8)
                     !end if
                 end do
             end do
         end do
         !$omp end do
         !$omp end parallel
         indSpin=indSpin+Glr%d%n1i*Glr%d%n2i*Glr%d%n3i
     end do
    
    END SUBROUTINE global_to_local_parallel


    !> Check the relative weight which the support functions have at the
    !! boundaries of the localization regions.
    subroutine get_boundary_weight(iproc, nproc, orbs, lzd, atoms, crmult, nsize_psi, psi, crit)
      use module_base
      use module_types, only: orbitals_data, local_zone_descriptors
      use module_atoms, only: atoms_data
      use yaml_output
      implicit none

      ! Calling arguments
      integer,intent(in) :: iproc, nproc
      type(orbitals_data),intent(in) :: orbs
      type(local_zone_descriptors),intent(in) :: lzd
      type(atoms_data),intent(in) :: atoms
      real(kind=8),intent(in) :: crmult
      integer,intent(in) :: nsize_psi
      real(kind=8),dimension(nsize_psi),intent(in) :: psi
      real(kind=8),intent(in) :: crit

      ! Local variables
      integer :: iorb, iiorb, ilr, iseg, jj, j0, j1, ii, i3, i2, i0, i1, i, ind, iat, iatype
      integer :: ij3, ij2, ij1, jj3, jj2, jj1, ijs3, ijs2, ijs1, ije3, ije2, ije1, nwarnings
      real(kind=8) :: h, x, y, z, d, weight_inside, weight_boundary, points_inside, points_boundary, ratio
      real(kind=8) :: atomrad, rad, boundary, weight_normalized, maxweight, meanweight
      real(kind=8),dimension(:),allocatable :: maxweight_types, meanweight_types
      integer,dimension(:),allocatable :: nwarnings_types, nsf_per_type
      logical :: perx, pery, perz, on_boundary

      call f_routine(id='get_boundary_weight')

      maxweight_types = f_malloc0(atoms%astruct%ntypes,id='maxweight_types')
      meanweight_types = f_malloc0(atoms%astruct%ntypes,id='maxweight_types')
      nwarnings_types = f_malloc0(atoms%astruct%ntypes,id='nwarnings_types')
      nsf_per_type = f_malloc0(atoms%astruct%ntypes,id='nsf_per_type')

      if (iproc==0) then
          call yaml_sequence(advance='no')
      end if

      ! mean value of the grid spacing
      h = sqrt(lzd%hgrids(1)**2+lzd%hgrids(2)**2+lzd%hgrids(3)**2)

      ! periodicity in the three directions
      perx=(lzd%glr%geocode /= 'F')
      pery=(lzd%glr%geocode == 'P')
      perz=(lzd%glr%geocode /= 'F')

      ! For perdiodic boundary conditions, one has to check also in the neighboring
      ! cells (see in the loop below)
      if (perx) then
          ijs1 = -1
          ije1 = 1
      else
          ijs1 = 0
          ije1 = 0
      end if
      if (pery) then
          ijs2 = -1
          ije2 = 1
      else
          ijs2 = 0
          ije2 = 0
      end if
      if (perz) then
          ijs3 = -1
          ije3 = 1
      else
          ijs3 = 0
          ije3 = 0
      end if

      nwarnings = 0
      maxweight = 0.d0
      meanweight = 0.d0
      if (orbs%norbp>0) then
          ind = 0
          do iorb=1,orbs%norbp
              iiorb = orbs%isorb + iorb
              ilr = orbs%inwhichlocreg(iiorb)

              iat = orbs%onwhichatom(iiorb)
              iatype = atoms%astruct%iatype(iat)
              atomrad = atoms%radii_cf(iatype,1)*crmult
              rad = atoms%radii_cf(atoms%astruct%iatype(iat),1)*crmult

              boundary = min(rad,lzd%llr(ilr)%locrad)
              !write(*,*) 'rad, locrad, boundary', rad, lzd%llr(ilr)%locrad, boundary

              nsf_per_type(iatype) = nsf_per_type(iatype ) + 1

              weight_boundary = 0.d0
              weight_inside = 0.d0
              points_inside = 0.d0
              points_boundary = 0.d0
              do iseg=1,lzd%llr(ilr)%wfd%nseg_c
                  jj=lzd%llr(ilr)%wfd%keyvglob(iseg)
                  j0=lzd%llr(ilr)%wfd%keyglob(1,iseg)
                  j1=lzd%llr(ilr)%wfd%keyglob(2,iseg)
                  ii=j0-1
                  i3=ii/((lzd%glr%d%n1+1)*(lzd%glr%d%n2+1))
                  ii=ii-i3*(lzd%glr%d%n1+1)*(lzd%glr%d%n2+1)
                  i2=ii/(lzd%glr%d%n1+1)
                  i0=ii-i2*(lzd%glr%d%n1+1)
                  i1=i0+j1-j0
                  do i=i0,i1
                      ind = ind + 1
                      on_boundary = .false.
                      do ij3=ijs3,ije3!-1,1
                          jj3=i3+ij3*(lzd%glr%d%n3+1)
                          z = real(jj3,kind=8)*lzd%hgrids(3)
                          do ij2=ijs2,ije2!-1,1
                              jj2=i2+ij2*(lzd%glr%d%n2+1)
                              y = real(jj2,kind=8)*lzd%hgrids(2)
                              do ij1=ijs1,ije1!-1,1
                                  jj1=i+ij1*(lzd%glr%d%n1+1)
                                  x = real(i,kind=8)*lzd%hgrids(1)
                                  d = sqrt((x-lzd%llr(ilr)%locregcenter(1))**2 + &
                                           (y-lzd%llr(ilr)%locregcenter(2))**2 + &
                                           (z-lzd%llr(ilr)%locregcenter(3))**2)
                                  if (abs(d-boundary)<h) then
                                      on_boundary=.true.
                                  end if
                              end do
                          end do
                      end do
                      if (on_boundary) then
                          ! This value is on the boundary
                          !write(*,'(a,2f9.2,3i8,3es16.8)') 'on boundary: boundary, d, i1, i2, i3, x, y, z', &
                          !    boundary, d, i, i2, i3, x, y, z
                          weight_boundary = weight_boundary + psi(ind)**2
                          points_boundary = points_boundary + 1.d0
                      else
                          weight_inside = weight_inside + psi(ind)**2
                          points_inside = points_inside + 1.d0
                      end if
                  end do
              end do
              ! fine part, to be done only if nseg_f is nonzero
              do iseg=lzd%llr(ilr)%wfd%nseg_c+1,lzd%llr(ilr)%wfd%nseg_c+lzd%llr(ilr)%wfd%nseg_f
                  jj=lzd%llr(ilr)%wfd%keyvglob(iseg)
                  j0=lzd%llr(ilr)%wfd%keyglob(1,iseg)
                  j1=lzd%llr(ilr)%wfd%keyglob(2,iseg)
                  ii=j0-1
                  i3=ii/((lzd%glr%d%n1+1)*(lzd%glr%d%n2+1))
                  ii=ii-i3*(lzd%glr%d%n1+1)*(lzd%glr%d%n2+1)
                  i2=ii/(lzd%glr%d%n1+1)
                  i0=ii-i2*(lzd%glr%d%n1+1)
                  i1=i0+j1-j0
                  do i=i0,i1
                      ind = ind + 7
                      on_boundary = .false.
                      do ij3=ijs3,ije3!-1,1
                          jj3=i3+ij3*(lzd%glr%d%n3+1)
                          z = real(jj3,kind=8)*lzd%hgrids(3)
                          do ij2=ijs2,ije2!-1,1
                              jj2=i2+ij2*(lzd%glr%d%n2+1)
                              y = real(jj2,kind=8)*lzd%hgrids(2)
                              do ij1=ijs1,ije1!-1,1
                                  jj1=i+ij1*(lzd%glr%d%n1+1)
                                  x = real(i,kind=8)*lzd%hgrids(1)
                                  d = sqrt((x-lzd%llr(ilr)%locregcenter(1))**2 + &
                                           (y-lzd%llr(ilr)%locregcenter(2))**2 + &
                                           (z-lzd%llr(ilr)%locregcenter(3))**2)
                                  if (abs(d-boundary)<h) then
                                      on_boundary=.true.
                                  end if
                              end do
                          end do
                      end do
                      if (on_boundary) then
                          ! This value is on the boundary
                          !write(*,'(a,f9.2,3i8,3es16.8)') 'on boundary: d, i1, i2, i3, x, y, z', d, i, i2, i3, x, y, z
                          weight_boundary = weight_boundary + psi(ind-6)**2
                          weight_boundary = weight_boundary + psi(ind-5)**2
                          weight_boundary = weight_boundary + psi(ind-4)**2
                          weight_boundary = weight_boundary + psi(ind-3)**2
                          weight_boundary = weight_boundary + psi(ind-2)**2
                          weight_boundary = weight_boundary + psi(ind-1)**2
                          weight_boundary = weight_boundary + psi(ind-0)**2
                          points_boundary = points_boundary + 7.d0
                      else
                          weight_inside = weight_inside + psi(ind-6)**2
                          weight_inside = weight_inside + psi(ind-5)**2
                          weight_inside = weight_inside + psi(ind-4)**2
                          weight_inside = weight_inside + psi(ind-3)**2
                          weight_inside = weight_inside + psi(ind-2)**2
                          weight_inside = weight_inside + psi(ind-1)**2
                          weight_inside = weight_inside + psi(ind-0)**2
                          points_inside = points_inside + 7.d0
                      end if
                  end do
              end do
              ! Ratio of the points on the boundary with resepct to the total number of points
              ratio = points_boundary/(points_boundary+points_inside)
              weight_normalized = weight_boundary/ratio
              meanweight = meanweight + weight_normalized
              maxweight = max(maxweight,weight_normalized)
              meanweight_types(iatype) = meanweight_types(iatype) + weight_normalized
              maxweight_types(iatype) = max(maxweight_types(iatype),weight_normalized)
              if (weight_normalized>crit) then
                  nwarnings = nwarnings + 1
                  nwarnings_types(iatype) = nwarnings_types(iatype) + 1
              end if
              !write(*,'(a,i7,2f9.1,4es16.6)') 'iiorb, pi, pb, weight_inside, weight_boundary, ratio, xi', &
              !    iiorb, points_inside, points_boundary, weight_inside, weight_boundary, &
              !    points_boundary/(points_boundary+points_inside), &
              !    weight_boundary/ratio
          end do
          if (ind/=nsize_psi) then
              call f_err_throw('ind/=nsize_psi ('//trim(yaml_toa(ind))//'/='//trim(yaml_toa(nsize_psi))//')', &
                   err_name='BIGDFT_RUNTIME_ERROR')
          end if
      end if

      ! Sum up among all tasks... could use workarrays
      if (nproc>1) then
          call mpiallred(nwarnings, 1, mpi_sum, comm=bigdft_mpi%mpi_comm)
          call mpiallred(meanweight, 1, mpi_sum, comm=bigdft_mpi%mpi_comm)
          call mpiallred(maxweight, 1, mpi_max, comm=bigdft_mpi%mpi_comm)
          call mpiallred(nwarnings_types, mpi_sum, comm=bigdft_mpi%mpi_comm)
          call mpiallred(meanweight_types, mpi_sum, comm=bigdft_mpi%mpi_comm)
          call mpiallred(maxweight_types, mpi_max, comm=bigdft_mpi%mpi_comm)
          call mpiallred(nsf_per_type, mpi_sum, comm=bigdft_mpi%mpi_comm)
      end if
      meanweight = meanweight/real(orbs%norb,kind=8)
      do iatype=1,atoms%astruct%ntypes
          meanweight_types(iatype) = meanweight_types(iatype)/real(nsf_per_type(iatype),kind=8)
      end do
      if (iproc==0) then
          call yaml_sequence_open('Check boundary values')
          call yaml_sequence(advance='no')
          call yaml_mapping_open(flow=.true.)
          call yaml_map('type','overall')
          call yaml_map('mean / max value',(/meanweight,maxweight/),fmt='(2es9.2)')
          call yaml_map('warnings',nwarnings)
          call yaml_mapping_close()
          do iatype=1,atoms%astruct%ntypes
              call yaml_sequence(advance='no')
              call yaml_mapping_open(flow=.true.)
              call yaml_map('type',trim(atoms%astruct%atomnames(iatype)))
              call yaml_map('mean / max value',(/meanweight_types(iatype),maxweight_types(iatype)/),fmt='(2es9.2)')
              call yaml_map('warnings',nwarnings_types(iatype))
              call yaml_mapping_close()
          end do
          call yaml_sequence_close()
      end if

      ! Print the warnings
      if (nwarnings>0) then
          if (iproc==0) then
              call yaml_warning('The support function localization radii might be too small, got'&
                  &//trim(yaml_toa(nwarnings))//' warnings')
          end if
      end if

      call f_free(maxweight_types)
      call f_free(meanweight_types)
      call f_free(nwarnings_types)
      call f_free(nsf_per_type)

      call f_release_routine()

    end subroutine get_boundary_weight


    subroutine small_to_large_locreg(iproc, npsidim_orbs_small, npsidim_orbs_large, lzdsmall, lzdlarge, &
           orbs, phismall, philarge, to_global)
      use module_base
      use module_types, only: orbitals_data, local_zone_descriptors
      implicit none
      
      ! Calling arguments
      integer,intent(in) :: iproc, npsidim_orbs_small, npsidim_orbs_large
      type(local_zone_descriptors),intent(in) :: lzdsmall, lzdlarge
      type(orbitals_data),intent(in) :: orbs
      real(kind=8),dimension(npsidim_orbs_small),intent(in) :: phismall
      real(kind=8),dimension(npsidim_orbs_large),intent(out) :: philarge
      logical,intent(in),optional :: to_global
      
      ! Local variables
      integer :: ists, istl, iorb, ilr, sdim, ldim, nspin
      logical :: global
    
      call f_routine(id='small_to_large_locreg')
    
      if (present(to_global)) then
          global=to_global
      else
          global=.false.
      end if
    
      call timing(iproc,'small2large','ON') ! lr408t 
      ! No need to put arrays to zero, Lpsi_to_global2 will handle this.
      call f_zero(philarge)
      ists=1
      istl=1
      do iorb=1,orbs%norbp
          ilr = orbs%inwhichLocreg(orbs%isorb+iorb)
          sdim=lzdsmall%llr(ilr)%wfd%nvctr_c+7*lzdsmall%llr(ilr)%wfd%nvctr_f
          if (global) then
              ldim=lzdsmall%glr%wfd%nvctr_c+7*lzdsmall%glr%wfd%nvctr_f
          else
              ldim=lzdlarge%llr(ilr)%wfd%nvctr_c+7*lzdlarge%llr(ilr)%wfd%nvctr_f
          end if
          nspin=1 !this must be modified later
          if (global) then
              call Lpsi_to_global2(iproc, sdim, ldim, orbs%norb, orbs%nspinor, nspin, lzdsmall%glr, &
                   lzdsmall%llr(ilr), phismall(ists), philarge(istl))
          else
              call Lpsi_to_global2(iproc, sdim, ldim, orbs%norb, orbs%nspinor, nspin, lzdlarge%llr(ilr), &
                   lzdsmall%llr(ilr), phismall(ists), philarge(istl))
          end if
          ists=ists+sdim
          istl=istl+ldim
      end do
      if(orbs%norbp>0 .and. ists/=npsidim_orbs_small+1) then
          write(*,'(3(a,i0))') 'ERROR on process ',iproc,': ',ists,'=ists /= npsidim_orbs_small+1=',npsidim_orbs_small+1
          stop
      end if
      if(orbs%norbp>0 .and. istl/=npsidim_orbs_large+1) then
          write(*,'(3(a,i0))') 'ERROR on process ',iproc,': ',istl,'=istl /= npsidim_orbs_large+1=',npsidim_orbs_large+1
          stop
      end if
           call timing(iproc,'small2large','OF') ! lr408t 
      call f_release_routine()
    end subroutine small_to_large_locreg


    !> Tranform one wavefunction between Global region and localisation region
    subroutine psi_to_locreg2(iproc, ldim, gdim, Llr, Glr, gpsi, lpsi)
    
      use module_base
      use module_types
     
     implicit none
    
      ! Subroutine Scalar Arguments
      integer,intent(in) :: iproc                  ! process ID
      integer,intent(in) :: ldim          ! dimension of lpsi 
      integer,intent(in) :: gdim          ! dimension of gpsi 
      type(locreg_descriptors),intent(in) :: Llr  ! Local grid descriptor
      type(locreg_descriptors),intent(in) :: Glr  ! Global grid descriptor
      
      !Subroutine Array Arguments
      real(wp),dimension(gdim),intent(in) :: gpsi       !Wavefunction (compressed format)
      real(wp),dimension(ldim),intent(out) :: lpsi   !Wavefunction in localization region
      
      !local variables
      integer :: igrid,isegloc,isegG,ix!,iorbs
      integer :: lmin,lmax,Gmin,Gmax
      integer :: icheck      ! check to make sure the dimension of loc_psi does not overflow 
      integer :: offset      ! gives the difference between the starting point of Lseg and Gseg
      integer :: length      ! Length of the overlap between Lseg and Gseg
      integer :: lincrement  ! Increment for writing orbitals in loc_psi
      integer :: Gincrement  ! Increment for reading orbitals in psi
      integer :: nseg        ! total number of segments in Llr
      integer, allocatable :: keymask(:,:)  ! shift for every segment of Llr (with respect to Glr)
      character(len=*), parameter :: subname='psi_to_locreg'
    !  integer :: i_stat,i_all
      integer :: start,Gstart
      integer :: isegstart,istart
    
      call f_routine(id=subname)
    
    ! Define integers
      nseg = Llr%wfd%nseg_c + Llr%wfd%nseg_f
      lincrement = Llr%wfd%nvctr_c + 7*Llr%wfd%nvctr_f
      Gincrement = Glr%wfd%nvctr_c + 7*Glr%wfd%nvctr_f
      icheck = 0
    
    ! Initialize loc_psi
      call f_zero(lpsi)
    
    ! Get the keymask: shift for every segment of Llr (with respect to Glr)
    ! allocate(keymask(2,nseg),stat=i_stat)
      keymask = f_malloc((/ 2, nseg /),id='keymask')
    
      call shift_locreg_indexes(Glr,Llr,keymask,nseg)
    
    
    !####################################################
    ! Do coarse region
    !####################################################
      isegstart=1
      icheck = 0
    
    
    !$omp parallel default(private) &
    !$omp shared(icheck,lpsi,gpsi,Glr,Llr,keymask,lincrement,Gincrement,Gstart) &
    !$omp firstprivate(isegstart,nseg)
    
      !$omp do reduction(+:icheck)
      local_loop_c: do isegloc = 1,Llr%wfd%nseg_c
         lmin = keymask(1,isegloc)
         lmax = keymask(2,isegloc)
         istart = llr%wfd%keyvloc(isegloc)-1
     
         global_loop_c: do isegG = isegstart,Glr%wfd%nseg_c
            Gmin = Glr%wfd%keygloc(1,isegG)
            Gmax = Glr%wfd%keygloc(2,isegG)
    
            ! For each segment in Llr check if there is a collision with the segment in Glr
            ! if not, cycle
            if(lmin > Gmax) then
                isegstart=isegG
            end if
            if(Gmin > lmax) exit global_loop_c
            if((lmin > Gmax) .or. (lmax < Gmin)) cycle global_loop_c
            
            ! Define the offset between the two segments
            offset = lmin - Gmin
            if(offset < 0) then
               offset = 0
            end if
        
            ! Define the length of the two segments
            length = min(lmax,Gmax)-max(lmin,Gmin)
    
            icheck = icheck + (length + 1)
     
            !Find the common elements and write them to the new localized wavefunction
            ! WARNING: index goes from 0 to length because it is the offset of the element
    
            do ix = 0,length
               istart = istart + 1
               lpsi(istart) = gpsi(Glr%wfd%keyvloc(isegG)+offset+ix)
            end do
         end do global_loop_c
      end do local_loop_c
      !$omp end do
    
    ! Check if the number of elements in loc_psi is valid
     ! if(icheck .ne. Llr%wfd%nvctr_c) then
       ! write(*,*)'There is an error in psi_to_locreg2: number of coarse points used',icheck
       ! write(*,*)'is not equal to the number of coarse points in the region',Llr%wfd%nvctr_c
     ! end if
    
    !##############################################################
    ! Now do fine region
    !##############################################################
    
      !icheck = 0
      start = Llr%wfd%nvctr_c
      Gstart = Glr%wfd%nvctr_c
    
      isegstart=Glr%wfd%nseg_c+1
    
      !$omp do reduction(+:icheck)
      local_loop_f: do isegloc = Llr%wfd%nseg_c+1,nseg
         lmin = keymask(1,isegloc)
         lmax = keymask(2,isegloc)
         istart = llr%wfd%keyvloc(isegloc)-1
     
         global_loop_f: do isegG = isegstart,Glr%wfd%nseg_c+Glr%wfd%nseg_f
    
            Gmin = Glr%wfd%keygloc(1,isegG)
            Gmax = Glr%wfd%keygloc(2,isegG)
    
            ! For each segment in Llr check if there is a collision with the segment in Glr
            ! if not, cycle
            if(lmin > Gmax) then
                isegstart=isegG
            end if
            if(Gmin > lmax)  exit global_loop_f
            if((lmin > Gmax) .or. (lmax < Gmin))  cycle global_loop_f
    
            offset = lmin - Gmin
            if(offset < 0) offset = 0
    
            length = min(lmax,Gmax)-max(lmin,Gmin)
    
            icheck = icheck + (length + 1)
    
            !Find the common elements and write them to the new localized wavefunction
            ! WARNING: index goes from 0 to length because it is the offset of the element
            do ix = 0,length
               istart = istart+1
               do igrid=1,7
                  lpsi(start+(istart-1)*7+igrid) = gpsi(Gstart+(Glr%wfd%keyvloc(isegG)+offset+ix-1)*7+igrid)
               end do
            end do
         end do global_loop_f
      end do local_loop_f
      !$omp end do
    
      !$omp end parallel
    
     !! Check if the number of elements in loc_psi is valid
      if(icheck .ne. Llr%wfd%nvctr_f+Llr%wfd%nvctr_c) then
        write(*,'(a,i0,a,i0)')'process ',iproc,': There is an error in psi_to_locreg: number of fine points used ',icheck
        write(*,'(a,i0)')'is not equal to the number of fine points in the region ',Llr%wfd%nvctr_f+Llr%wfd%nvctr_c
      end if
    
    
    
    !  i_all=-product(shape(keymask))*kind(keymask)
    ! deallocate(keymask,stat=i_stat)
      call f_free(keymask)
      call f_release_routine()
    
    END SUBROUTINE psi_to_locreg2



    !> Find the shift necessary for the indexes of every segment of Blr
    !!   to make them compatible with the indexes of Alr. These shifts are
    !!   returned in the array keymask(nseg), where nseg should be the number
    !!   of segments in Blr.
    !! @warning 
    !!   This routine supposes that the region Blr is contained in the region Alr.
    !!   This should always be the case, if we concentrate on the overlap between two regions.
    subroutine shift_locreg_indexes(Alr,Blr,keymask,nseg)
    
      use module_base
      use module_types
     
     implicit none
    
    ! Arguments
     type(locreg_descriptors),intent(in) :: Alr,Blr   ! The two localization regions
     integer,intent(in) :: nseg
     integer,intent(out) :: keymask(2,nseg)
    
    ! Local variable
     integer :: iseg      !integer for the loop
     integer :: Bindex    !starting index of segments in Blr
     integer :: x,y,z     !coordinates of start of segments in Blr 
     integer :: shift(3)  !shift between the beginning of the segment in Blr and the origin of Alr
     integer ::  tmp
    
    
     ! This routine is only intended for conversions between locregs with the same boundary conditions.
     if (blr%geocode/='F') then
         call f_err_throw('shift_locreg_indexes can only be used for locregs with free boundary conditions', &
              err_name='BIGDFT_RUNTIME_ERROR')
     end if
    
    !Big loop on all segments
    !$omp parallel do default(private) shared(Blr,nseg,Alr,keymask)
     do iseg=1,nseg
    
    !##########################################
    ! For the Starting index
        Bindex = Blr%wfd%keygloc(1,iseg)
        tmp = Bindex -1
        z   = tmp / ((Blr%d%n2+1)*(Blr%d%n1+1))
        tmp = tmp - z*((Blr%d%n2+1)*(Blr%d%n1+1))
        y   = tmp / (Blr%d%n1+1)
        x   = tmp - y * (Blr%d%n1+1)
     
    ! Shift between the beginning of the segment and the start of the Alr region
        shift(1) = x + Blr%ns1 - Alr%ns1
        shift(2) = y + Blr%ns2 - Alr%ns2
        shift(3) = z + Blr%ns3 - Alr%ns3
    
    ! Write the shift in index form
        keymask(1,iseg) = shift(3)*(Alr%d%n1+1)*(Alr%d%n2+1) + shift(2)*(Alr%d%n1+1) + shift(1) + 1
    
    !######################################
    ! For the ending index
    
        Bindex = Blr%wfd%keygloc(2,iseg)
        tmp = Bindex -1
        z   = tmp / ((Blr%d%n2+1)*(Blr%d%n1+1))
        tmp = tmp - z*((Blr%d%n2+1)*(Blr%d%n1+1))
        y   = tmp / (Blr%d%n1+1)
        x   = tmp - y * (Blr%d%n1+1)
    
    ! Shift between the beginning of the segment and the start of the Alr region
        shift(1) = x + Blr%ns1 - Alr%ns1
        shift(2) = y + Blr%ns2 - Alr%ns2
        shift(3) = z + Blr%ns3 - Alr%ns3
    
    ! Write the shift in index form
        keymask(2,iseg) = shift(3)*(Alr%d%n1+1)*(Alr%d%n2+1) + shift(2)*(Alr%d%n1+1) + shift(1) + 1
     end do
    !$omp end parallel do
    
    END SUBROUTINE shift_locreg_indexes


    !> Projects a quantity stored with the global indexes (i1,i2,i3) within the localisation region.
    !! @warning: The quantity must not be stored in a compressed form.
    subroutine global_to_local(Glr,Llr,nspin,size_rho,size_Lrho,rho,Lrho)
    
     use module_base
     use module_types
     
     implicit none
    
    ! Arguments
     type(locreg_descriptors),intent(in) :: Llr   ! Local localization region
     type(locreg_descriptors),intent(in) :: Glr   ! Global localization region
     integer, intent(in) :: size_rho  ! size of rho array
     integer, intent(in) :: size_Lrho ! size of Lrho array
     integer, intent(in) :: nspin  !number of spins
     real(wp),dimension(size_rho),intent(in) :: rho  ! quantity in global region
     real(wp),dimension(size_Lrho),intent(out) :: Lrho ! piece of quantity in local region
    
    ! Local variable
     integer :: ispin,i1,i2,i3,ii1,ii2,ii3  !integer for loops
     integer :: indSmall, indSpin, indLarge ! indexes for the arrays
     logical:: z_inside, y_inside, x_inside
     integer:: iz, iy, m
     
    ! Cut out a piece of the quantity (rho) from the global region (rho) and
    ! store it in a local region (Lrho).
    
     if(Glr%geocode == 'F') then
         ! Use loop unrolling here
         indSmall=0
         indSpin=0
         do ispin=1,nspin
             ! WARNING: I added the factors 2.
             do i3=Llr%nsi3+1,Llr%d%n3i+Llr%nsi3
                 iz=(i3-1)*Glr%d%n2i*Glr%d%n1i
                 do i2=Llr%nsi2+1,Llr%d%n2i+Llr%nsi2
                     iy=(i2-1)*Glr%d%n1i
                     m=mod(Llr%d%n1i+Llr%nsi1-Llr%nsi1,4)
                     if(m/=0) then
                         do i1=Llr%nsi1+1,Llr%nsi1+m
                            indSmall=indSmall+1
                            indLarge=iz+iy+i1
                            Lrho(indSmall)=rho(indLarge+indSpin)
                         end do
                      end if
                      do i1=Llr%nsi1+1+m,Llr%d%n1i+Llr%nsi1,4
                         Lrho(indSmall+1)=rho(iz+iy+i1+0+indSpin)
                         Lrho(indSmall+2)=rho(iz+iy+i1+1+indSpin)
                         Lrho(indSmall+3)=rho(iz+iy+i1+2+indSpin)
                         Lrho(indSmall+4)=rho(iz+iy+i1+3+indSpin)
                         indSmall=indSmall+4
                      end do
                 end do
             end do
             indSpin=indSpin+Glr%d%n1i*Glr%d%n2i*Glr%d%n3i
         end do
     else
         ! General case
         indSmall=0
         indSpin=0
         do ispin=1,nspin
             ! WARNING: I added the factors 2.
             do ii3=Llr%nsi3+1,Llr%d%n3i+Llr%nsi3
                 i3 = mod(ii3-1,Glr%d%n3i)+1
                 z_inside = (i3>0 .and. i3<=Glr%d%n3i+1)
                 iz=(i3-1)*Glr%d%n2i*Glr%d%n1i
                 do ii2=Llr%nsi2+1,Llr%d%n2i+Llr%nsi2
                     i2 = mod(ii2-1,Glr%d%n2i)+1
                     y_inside = (i2>0 .and. i2<=Glr%d%n2i+1)
                     iy=(i2-1)*Glr%d%n1i
                     do ii1=Llr%nsi1+1,Llr%d%n1i+Llr%nsi1
                         i1 = mod(ii1-1,Glr%d%n1i)+1 
                         x_inside = (i1 > 0 .and. i1 <= Glr%d%n1i+1)
                         ! indSmall is the index in the local localization region
                         indSmall=indSmall+1
                         !!if (i3 > 0 .and. i2 > 0 .and. i1 > 0 .and.&                                       !This initializes the buffers of locreg to zeros if outside the simulation box.
                         !!    i3 <= Glr%d%n3i+1 .and. i2 <= Glr%d%n2i+1 .and. i1 <= Glr%d%n1i+1) then       !Should use periodic image instead... MUST FIX THIS.
                         !!   ! indLarge is the index in the global localization region. 
                         !!   indLarge=(i3-1)*Glr%d%n2i*Glr%d%n1i + (i2-1)*Glr%d%n1i + i1
                         if(z_inside .and. y_inside .and. x_inside) then
                            indLarge=iz+iy+i1
                            Lrho(indSmall)=rho(indLarge+indSpin)
                         else
                            Lrho(indSmall)= 0.0_wp
                         end if
                     end do
                 end do
             end do
             indSpin=indSpin+Glr%d%n1i*Glr%d%n2i*Glr%d%n3i
         end do
     end if
    
    END SUBROUTINE global_to_local


end module locreg_operations
