!> @file 
!!   Miscellaneous routines for linear toolbox
!! @author
!!   Copyright (C) 2011-2013 BigDFT group 
!!   This file is distributed under the terms of the
!!   GNU General Public License, see ~/COPYING file
!!   or http://www.gnu.org/copyleft/gpl.txt .
!!   For the list of contributors, see ~/AUTHORS 


!> Write the square of the wave functions (i.e. the orbital densities).
!! This routine can also be used to print the "support functions densities".
subroutine write_orbital_density(iproc, transform_to_global, iformat, &
           filename, npsidim, psi, orbs, lzd_g, at, rxyz, dens, lzd_l, in_which_locreg)
  use module_base
  use module_types
  !use module_interface2, except_this_one => write_orbital_density
  use locreg_operations, only: lpsi_to_global2
  use public_enums
  use bounds, only: locreg_bounds
  implicit none

  ! Calling arguments
  logical,intent(in) :: transform_to_global
  character(len=*),intent(in) :: filename
  integer,intent(in) :: iproc, npsidim, iformat
  real(kind=8),dimension(npsidim),intent(in),target :: psi
  type(orbitals_data),intent(in) :: orbs !< orbitals descriptors
  type(local_zone_descriptors),intent(inout) :: lzd_g !< global descriptors
  type(atoms_data),intent(in) :: at
  real(kind=8),dimension(3,at%astruct%nat),intent(in) :: rxyz
  logical,intent(in) :: dens !< density of wavefunctions or just wavefunctions
  type(local_zone_descriptors),intent(in),optional :: lzd_l !< local descriptors
  integer,dimension(orbs%norb),intent(in),optional :: in_which_locreg

  ! Local variables
  logical :: binary
  real(kind=8),dimension(:),pointer :: psi_g
  integer :: iunit0, iunitx, iunity, iunitz, iorb, ispinor, ist, ncount, iiorb, ilr
  integer :: iorb_out0, iorb_outx, iorb_outy, iorb_outz, sdim, ldim
  character(len=500) :: filebase0, filebasex, filebasey, filebasez
  character(len=500) :: file0, filex, filey, filez
  integer,dimension(:),allocatable :: in_which_locreg_

  call f_routine(id='write_orbital_density')

  if (transform_to_global) then
      if (.not.present(lzd_l)) call f_err_throw('lzd_l not present',err_name='BIGDFT_RUNTIME_ERROR')
  end if

  in_which_locreg_ = f_malloc(orbs%norb,id='in_which_locreg_')
  if (present(in_which_locreg)) then
      !in_which_locreg_ = in_which_locreg
      call f_memcpy(src=in_which_locreg, dest=in_which_locreg_)
  else
      !in_which_locreg_ = orbs%inwhichlocreg
      call f_memcpy(src=orbs%inwhichlocreg, dest=in_which_locreg_)
  end if

  !!! Transform to the global region
  !!if (transform_to_global) then
  !!    psi_g = f_malloc_ptr(orbs%norbp*(lzd_g%glr%wfd%nvctr_c+7*lzd_g%glr%wfd%nvctr_f), id='psi_g')
  !!    call small_to_large_locreg(iproc, npsidim, &
  !!         orbs%norbp*(lzd_l%glr%wfd%nvctr_c+7*lzd_l%glr%wfd%nvctr_f), lzd_l, &
  !!         lzd_g, orbs, psi, psi_g, to_global=.true.)
  !!else
  !!    psi_g => psi
  !!end if

  psi_g = f_malloc_ptr(lzd_g%glr%wfd%nvctr_c+7*lzd_g%glr%wfd%nvctr_f, id='psi_g')


  binary = (iformat==WF_FORMAT_BINARY)

  ! Need to create the convolution bounds
  ! check first if already allocated - maybe should be more thorough than just checking one array?
  if (.not. associated(lzd_g%glr%bounds%kb%ibyz_c)) then
     write(*,*) 'before calling locreg_bounds'
     call locreg_bounds(lzd_g%glr%d%n1, lzd_g%glr%d%n2, lzd_g%glr%d%n3, &
          lzd_g%glr%d%nfl1, lzd_g%glr%d%nfu1, &
          lzd_g%glr%d%nfl2, lzd_g%glr%d%nfu2, &
          lzd_g%glr%d%nfl3, lzd_g%glr%d%nfu3, &
          lzd_g%glr%wfd, lzd_g%glr%bounds)
     write(*,*) 'after locreg_bounds'
  end if

  ist = 1
  ncount = lzd_g%glr%wfd%nvctr_c+7*lzd_g%glr%wfd%nvctr_f
  do iorb=1,orbs%norbp
      iiorb = orbs%isorb + iorb
      !ilr = orbs%inwhichlocreg(iiorb)
      ilr = in_which_locreg_(iiorb)
      if (transform_to_global) then
          sdim=lzd_l%llr(ilr)%wfd%nvctr_c+7*lzd_l%llr(ilr)%wfd%nvctr_f
          ldim=lzd_l%glr%wfd%nvctr_c+7*lzd_l%glr%wfd%nvctr_f
          call f_zero(psi_g)
          call lpsi_to_global2(iproc, sdim, ldim, orbs%norb, 1, lzd_l%glr, &
               lzd_l%llr(ilr), psi(ist:), psi_g)
      else
          !sdim=lzd_g%llr(ilr)%wfd%nvctr_c+7*lzd_g%llr(ilr)%wfd%nvctr_f
          sdim=lzd_g%glr%wfd%nvctr_c+7*lzd_g%glr%wfd%nvctr_f
          call vcopy(sdim, psi(ist), 1, psi_g(1), 1)
      end if
      ist = ist + sdim
      do ispinor=1,orbs%nspinor
          if (orbs%nspinor/=1) stop 'write_orbital_density not implemented for nspinor/=1'
          call plot_one_orbdens(lzd_g%glr, at, orbs, rxyz, lzd_g%hgrids, filename, &
               iorb, ispinor, binary, psi_g, dens)
          !iunit0 = 101
          !iunit0 = 102
          !iunit0 = 103
          !iunit0 = 104
          !!!call open_filename_of_iorb(iunit0, binary, filename, orbs, iorb, ispinor, iorb_out0)
          !!!call open_filename_of_iorb(iunitx, binary, filename, orbs, iorb, ispinor, iorb_outx)
          !!!call open_filename_of_iorb(iunity, binary, filename, orbs, iorb, ispinor, iorb_outy)
          !!!call open_filename_of_iorb(iunitz, binary, filename, orbs, iorb, ispinor, iorb_outz)
          !call filename_of_iorb(binary, trim(input%dir_output)//filename, orbs, iorb, ispinor, filebase0, iorb_out0)
          !call filename_of_iorb(binary, trim(input%dir_output)//filename, orbs, iorb, ispinor, filebasex, iorb_outx)
          !call filename_of_iorb(binary, trim(input%dir_output)//filename, orbs, iorb, ispinor, filebasey, iorb_outy)
          !call filename_of_iorb(binary, trim(input%dir_output)//filename, orbs, iorb, ispinor, filebasez, iorb_outz)
          !file0 = trim(filebase0)//'.cube'
          !filex = trim(filebasex)//'.cube'
          !filey = trim(filebasey)//'.cube'
          !filez = trim(filebasez)//'.cube'
          !write(*,*) 'file0',file0
          !call f_open_file(iunit0, file=file0, binary=binary)
          !call f_open_file(iunitx, file=filex, binary=binary)
          !call f_open_file(iunity, file=filey, binary=binary)
          !call f_open_file(iunitz, file=filez, binary=binary)
          !!write(*,'(a,6i9)') 'iproc, iorb, iunit0, iunitx, iunity, iunitz',iproc, iorb, iunit0, iunitx, iunity, iunitz
          !call plot_wf(.true.,'', 2, at, 1.d0, lzd_g%glr, &
          !     lzd_g%hgrids(1), lzd_g%hgrids(2), lzd_g%hgrids(3), &
          !     rxyz, psi_g, &
          !     iunit0, iunitx, iunity, iunitz)
          !call f_close(iunit0)
          !call f_close(iunitx)
          !call f_close(iunity)
          !call f_close(iunitz)
          !!ist = ist + ncount
      end do
  end do

  !call deallocate_bounds(lzd_g%glr%geocode, lzd_g%glr%hybrid_on, lzd_g%glr%bounds)

  !if (transform_to_global) then
      call f_free_ptr(psi_g)
  call f_free(in_which_locreg_)
  !end if

  call f_release_routine()

end subroutine write_orbital_density



subroutine plot_one_orbdens(lr, at, orbs, rxyz, hgrids, filename, iorb, ispinor, binary, psi_g, dens)
  use module_base
  use module_types
  use module_interfaces, only: filename_of_iorb, plot_wf
  implicit none

  ! Calling arguments
  type(locreg_descriptors),intent(in) :: lr
  type(atoms_data),intent(in) :: at
  type(orbitals_data),intent(in) :: orbs
  real(kind=8),dimension(3,at%astruct%nat),intent(in) :: rxyz
  real(kind=8),dimension(3),intent(in) :: hgrids
  character(len=*),intent(in) :: filename
  integer,intent(in) :: iorb, ispinor
  logical,intent(in) :: binary
  real(kind=8),dimension(lr%wfd%nvctr_c+7*lr%wfd%nvctr_f),intent(in) :: psi_g
  logical,intent(in) :: dens !< density of wavefunction or just wavefunction

  ! Local variables
  integer :: iunit0, iunitx, iunity, iunitz, iorb_out0, iorb_outx, iorb_outy, iorb_outz
  character(len=500) :: filebase0, filebasex, filebasey, filebasez
  character(len=500) :: file0, filex, filey, filez

  iunit0 = 101
  iunitx = 102
  iunity = 103
  iunitz = 104
  call filename_of_iorb(binary, filename, orbs, iorb, ispinor, filebase0, iorb_out0)
  call filename_of_iorb(binary, filename, orbs, iorb, ispinor, filebasex, iorb_outx)
  call filename_of_iorb(binary, filename, orbs, iorb, ispinor, filebasey, iorb_outy)
  call filename_of_iorb(binary, filename, orbs, iorb, ispinor, filebasez, iorb_outz)
  file0 = trim(filebase0)//'.cube'
  filex = trim(filebasex)//'_x'
  filey = trim(filebasey)//'_y'
  filez = trim(filebasez)//'_z'
  !write(*,*) 'file0',file0
  !call f_open_file(iunit0, file=file0, binary=binary)
  !call f_open_file(iunitx, file=filex, binary=binary)
  !call f_open_file(iunity, file=filey, binary=binary)
  !call f_open_file(iunitz, file=filez, binary=binary)
  if (dens) then
     call plot_wf(.true.,trim(filebase0), 2, at, 1.d0, lr, &
          hgrids(1), hgrids(2), hgrids(3), &
          rxyz, psi_g, &
          iunit0, iunitx, iunity, iunitz)
  else
     call plot_wf(.true.,trim(filebase0), 1, at, 1.d0, lr, &
          hgrids(1), hgrids(2), hgrids(3), &
          rxyz, psi_g, &
          iunit0, iunitx, iunity, iunitz)
  end if
  !call f_close(iunit0)
  !call f_close(iunitx)
  !call f_close(iunity)
  !call f_close(iunitz)

end subroutine plot_one_orbdens
 

!> Plots the orbitals
subroutine plotOrbitals(iproc, tmb, phi, nat, rxyz, hxh, hyh, hzh, it, basename)
   use module_base
   use module_types
   use locreg_operations
   implicit none

   ! Calling arguments
   integer :: iproc
   type(DFT_wavefunction),intent(in) :: tmb
   real(kind=8), dimension((tmb%lzd%glr%wfd%nvctr_c+7*tmb%lzd%glr%wfd%nvctr_f)*tmb%orbs%nspinor*tmb%orbs%norbp) :: phi
   integer :: nat
   real(kind=8), dimension(3,nat) :: rxyz
   real(kind=8) :: hxh, hyh, hzh
   integer :: it
   character(len=4),intent(in) :: basename

   integer :: ix, iy, iz, ix0, iy0, iz0, iiAt, jj, iorb, i1, i2, i3, istart, ii, iat
   integer :: unit1, unit2, unit3, unit4, unit5, unit6, unit7, unit8, unit9, unit10, unit11, unit12
   integer :: ixx, iyy, izz, maxid, i, ixmin, ixmax, iymin, iymax, izmin, izmax
   integer :: ilr
   real(kind=8) :: dixx, diyy, dizz, prevdiff, maxdiff, diff, dnrm2
   real(kind=8), dimension(:), allocatable :: phir
   !real(kind=8), dimension(:,:,:), allocatable :: psig_c
   real(kind=8),dimension(3) :: rxyzdiff
   real(kind=8),dimension(3,11) :: rxyzref
   integer,dimension(4) :: closeid
   type(workarr_sumrho) :: w
   character(len=10) :: c1, c2, c3
   character(len=50) :: file1, file2, file3, file4, file5, file6, file7, file8, file9, file10, file11, file12
   logical :: dowrite, plot_axis, plot_diagonals, plot_neighbors

   plot_axis=.true.
   plot_diagonals=.false.!.true.
   plot_neighbors=.false.

   phir = f_malloc(tmb%lzd%glr%d%n1i*tmb%lzd%glr%d%n2i*tmb%lzd%glr%d%n3i,id='phir')

   call initialize_work_arrays_sumrho(tmb%lzd%glr,.true.,w)
   rxyzref=-555.55d0

   istart=0

   unit1 =20*(iproc+1)+3
   unit2 =20*(iproc+1)+4
   unit3 =20*(iproc+1)+5
   unit4 =20*(iproc+1)+6
   unit5 =20*(iproc+1)+7
   unit6 =20*(iproc+1)+8
   unit7 =20*(iproc+1)+9
   unit8 =20*(iproc+1)+10
   unit9 =20*(iproc+1)+11
   unit10=20*(iproc+1)+12
   unit11=20*(iproc+1)+13
   unit12=20*(iproc+1)+14

   !write(*,*) 'write, tmb%orbs%nbasisp', tmb%orbs%norbp
    orbLoop: do iorb=1,tmb%orbs%norbp
        !!phir=0.d0
        call f_zero(phir)
        call daub_to_isf(tmb%lzd%glr,w,phi(istart+1),phir(1))
        ilr=tmb%orbs%inwhichlocreg(tmb%orbs%isorb+iorb)
        iiAt=tmb%orbs%onwhichatom(tmb%orbs%isorb+iorb)
        ix0=nint(rxyz(1,iiAt)/hxh)!-15
        iy0=nint(rxyz(2,iiAt)/hyh)!-15
        iz0=nint(rxyz(3,iiAt)/hzh)!-15

        if (plot_neighbors) then
            ! Search the four closest atoms
            prevdiff=1.d-5 ! the same atom
            do i=1,4
                maxdiff=1.d100
                do iat=1,nat
                    rxyzdiff(:)=rxyz(:,iat)-rxyz(:,iiat)
                    diff=dnrm2(3,rxyzdiff,1)
                    if (diff<maxdiff .and. diff>prevdiff) then
                        maxdiff=diff
                        maxid=iat
                    end if
                end do
                closeid(i)=maxid
                prevdiff=maxdiff*1.00001d0 !just to be sure that not twice the same is chosen
            end do
        end if

        write(c1,'(i5.5)') iproc
        write(c2,'(i5.5)') iorb
        write(c3,'(i5.5)') it
        file1=basename//'_'//trim(c1)//'_'//trim(c2)//'_'//trim(c3)//'_x'
        file2=basename//'_'//trim(c1)//'_'//trim(c2)//'_'//trim(c3)//'_y'
        file3=basename//'_'//trim(c1)//'_'//trim(c2)//'_'//trim(c3)//'_z'
        file4=basename//'_'//trim(c1)//'_'//trim(c2)//'_'//trim(c3)//'_pxpypz'
        file5=basename//'_'//trim(c1)//'_'//trim(c2)//'_'//trim(c3)//'_mxpypz'
        file6=basename//'_'//trim(c1)//'_'//trim(c2)//'_'//trim(c3)//'_mxmypz'
        file7=basename//'_'//trim(c1)//'_'//trim(c2)//'_'//trim(c3)//'_pxmypz'
        file8=basename//'_'//trim(c1)//'_'//trim(c2)//'_'//trim(c3)//'_1st'
        file9=basename//'_'//trim(c1)//'_'//trim(c2)//'_'//trim(c3)//'_2nd'
        file10=basename//'_'//trim(c1)//'_'//trim(c2)//'_'//trim(c3)//'_3rd'
        file11=basename//'_'//trim(c1)//'_'//trim(c2)//'_'//trim(c3)//'_4th'
        file12=basename//'_'//trim(c1)//'_'//trim(c2)//'_'//trim(c3)//'_info'
        if (plot_axis) then
            open(unit=unit1, file=trim(file1))
            open(unit=unit2, file=trim(file2))
            open(unit=unit3, file=trim(file3))
        end if
        if (plot_diagonals) then
            open(unit=unit4, file=trim(file4))
            open(unit=unit5, file=trim(file5))
            open(unit=unit6, file=trim(file6))
            open(unit=unit7, file=trim(file7))
        end if
        if (plot_neighbors) then
            open(unit=unit8, file=trim(file8))
            open(unit=unit9, file=trim(file9))
            open(unit=unit10, file=trim(file10))
            open(unit=unit11, file=trim(file11))
        end if
        open(unit=unit12, file=trim(file12))
        
        !write(unit1,'(a,3i8)') '# ix0, iy0, iz0 ',ix0,iy0,iz0
        !write(unit2,'(a,3i8)') '# ix0, iy0, iz0 ',ix0,iy0,iz0
        !write(unit3,'(a,3i8)') '# ix0, iy0, iz0 ',ix0,iy0,iz0
        !write(unit4,'(a,3i8)') '# ix0, iy0, iz0 ',ix0,iy0,iz0
        !write(unit5,'(a,3i8)') '# ix0, iy0, iz0 ',ix0,iy0,iz0
        !write(unit6,'(a,3i8)') '# ix0, iy0, iz0 ',ix0,iy0,iz0
        !write(unit7,'(a,3i8)') '# ix0, iy0, iz0 ',ix0,iy0,iz0



        !!!! TEMPORARY ######################################
        !!!allocate(psig_c(0:tmb%lzd%glr%d%n1,0:tmb%lzd%glr%d%n2,0:tmb%lzd%glr%d%n3))
        !!!psig_c=0.d0
        !!!do iseg=1,tmb%lzd%glr%wfd%nseg_c
        !!!   jj=tmb%lzd%glr%wfd%keyvloc(iseg)
        !!!   j0=tmb%lzd%glr%wfd%keygloc(1,iseg)
        !!!   j1=tmb%lzd%glr%wfd%keygloc(2,iseg)
        !!!   ii=j0-1
        !!!   i3=ii/((tmb%lzd%glr%d%n1+1)*(tmb%lzd%glr%d%n2+1))
        !!!   ii=ii-i3*(tmb%lzd%glr%d%n1+1)*(tmb%lzd%glr%d%n2+1)
        !!!   i2=ii/(tmb%lzd%glr%d%n1+1)
        !!!   i0=ii-i2*(tmb%lzd%glr%d%n1+1)
        !!!   i1=i0+j1-j0
        !!!   do i=i0,i1
        !!!      psig_c(i,i2,i3)=phi(istart+i-i0+jj)
        !!!   enddo
        !!!enddo
        !!!ixmax=-1000000
        !!!ixmin= 1000000
        !!!iymax=-1000000
        !!!iymin= 1000000
        !!!izmax=-1000000
        !!!izmin= 1000000
        !!!do i3=0,tmb%lzd%glr%d%n3
        !!!   do i2=0,tmb%lzd%glr%d%n2
        !!!      do i1=0,tmb%lzd%glr%d%n1
        !!!         if (psig_c(i1,i2,i3)/=0.d0) then
        !!!             if (i1>ixmax) ixmax=i1
        !!!             if (i1<ixmin) ixmin=i1
        !!!             if (i2>iymax) iymax=i2
        !!!             if (i2<iymin) iymin=i2
        !!!             if (i3>izmax) izmax=i3
        !!!             if (i3<izmin) izmin=i3
        !!!         end if
        !!!      end do
        !!!   end do
        !!!end do

        !!!!!ixmax=ixmax-tmb%lzd%llr(ilr)%ns1
        !!!!!ixmin=ixmin-tmb%lzd%llr(ilr)%ns1
        !!!!!iymax=iymax-tmb%lzd%llr(ilr)%ns2
        !!!!!iymin=iymin-tmb%lzd%llr(ilr)%ns2
        !!!!!izmax=izmax-tmb%lzd%llr(ilr)%ns3
        !!!!!izmin=izmin-tmb%lzd%llr(ilr)%ns3

        !!!deallocate(psig_c)
        !!!write(unit12,'(a,2i6,3x,6i9)') '# id, ilr, ixconf, ix0, ixmin, ixmax, ns1, n1', & 
        !!!    tmb%orbs%isorb+iorb, ilr, nint(tmb%confdatarr(iorb)%rxyzconf(1)/(2.d0*hxh)), nint(rxyz(1,iiat)/(2.d0*hxh)), ixmin, ixmax, tmb%lzd%llr(ilr)%ns1, tmb%lzd%llr(ilr)%d%n1
        !!!write(unit12,'(a,2i6,3x,6i9)') '# id, ilr, iyconf, iy0, iymin, iymax, ns2, n2', &
        !!!    tmb%orbs%isorb+iorb, ilr, nint(tmb%confdatarr(iorb)%rxyzconf(2)/(2.d0*hyh)), nint(rxyz(2,iiat)/(2.d0*hyh)), iymin, iymax, tmb%lzd%llr(ilr)%ns2, tmb%lzd%llr(ilr)%d%n2
        !!!write(unit12,'(a,2i6,3x,6i9)') '# id, ilr, izconf, iz0, izmin, izmax, ns3, n3', &
        !!!    tmb%orbs%isorb+iorb, ilr, nint(tmb%confdatarr(iorb)%rxyzconf(3)/(2.d0*hzh)), nint(rxyz(3,iiat)/(2.d0*hzh)), izmin, izmax, tmb%lzd%llr(ilr)%ns3, tmb%lzd%llr(ilr)%d%n3

        !!!! END TEMPORARY ##################################





        ixmax=-1000000
        ixmin= 1000000
        iymax=-1000000
        iymin= 1000000
        izmax=-1000000
        izmin= 1000000
        jj=0
        do i3=1,tmb%lzd%glr%d%n3i
            do i2=1,tmb%lzd%glr%d%n2i
                do i1=1,tmb%lzd%glr%d%n1i
                   jj=jj+1
                   ! z component of point jj
                   iz=jj/(tmb%lzd%glr%d%n2i*tmb%lzd%glr%d%n1i)
                   ! Subtract the 'lower' xy layers
                   ii=jj-iz*(tmb%lzd%glr%d%n2i*tmb%lzd%glr%d%n1i)
                   ! y component of point jj
                   iy=ii/tmb%lzd%glr%d%n1i
                   ! Subtract the 'lower' y rows
                   ii=ii-iy*tmb%lzd%glr%d%n1i
                   ! x component
                   ix=ii
!if(phir(jj)>1.d0) write(*,'(a,3i7,es15.6)') 'WARNING: ix, iy, iz, phir(jj)', ix, iy, iz, phir(jj)


                   ! Shift the values due to the convolutions bounds
                   ix=ix-14
                   iy=iy-14
                   iz=iz-14
                   
                   if (phir(jj)/=0.d0) then
                       if (ix>ixmax) ixmax=ix
                       if (ix<ixmin) ixmin=ix
                       if (iy>iymax) iymax=iy
                       if (iy<iymin) iymin=iy
                       if (iz>izmax) izmax=iz
                       if (iz<izmin) izmin=iz
                   end if

                   ixx=ix-ix0
                   iyy=iy-iy0
                   izz=iz-iz0
                   dixx=hxh*dble(ixx)
                   diyy=hyh*dble(iyy)
                   dizz=hzh*dble(izz)

                   if (plot_axis) then
                       ! Write along x-axis
                       if(iyy==0 .and. izz==0) write(unit1,'(2es18.10)') dixx, phir(jj)

                       ! Write along y-axis
                       if(ixx==0 .and. izz==0) write(unit2,'(2es18.10)') diyy, phir(jj)

                       ! Write along z-axis
                       if(ixx==0 .and. iyy==0) write(unit3,'(2es18.10)') dizz, phir(jj)
                   end if

                   if (plot_diagonals) then
                       ! Write diagonal in octant +x,+y,+z
                       if (ixx==iyy .and. ixx==izz .and. iyy==izz) then
                           write(unit4,'(2es18.10)') sqrt(dixx**2+diyy**2+dizz**2)*dsign(1.d0,dizz), phir(jj)
                       end if

                       ! Write diagonal in octant -x,+y,+z
                       if (-ixx==iyy .and. -ixx==izz .and. iyy==izz) then
                           write(unit5,'(2es18.10)') sqrt(dixx**2+diyy**2+dizz**2)*dsign(1.d0,dizz), phir(jj)
                       end if

                       ! Write diagonal in octant -x,-y,+z
                       if (-ixx==-iyy .and. -ixx==izz .and. -iyy==izz) then
                           write(unit6,'(2es18.10)') sqrt(dixx**2+diyy**2+dizz**2)*dsign(1.d0,dizz), phir(jj)
                       end if

                       ! Write diagonal in octant +x,-y,+z
                       if (ixx==-iyy .and. ixx==izz .and. -iyy==izz) then
                           write(unit7,'(2es18.10)') sqrt(dixx**2+diyy**2+dizz**2)*dsign(1.d0,dizz), phir(jj)
                       end if
                   end if

                   if (plot_neighbors) then
                       ! Write along line in direction of the closest atom
                       dowrite=gridpoint_close_to_straightline(ix, iy, iz, &
                           rxyz(1,iiat), rxyz(1,closeid(1)), hxh, hyh, hzh)
                       if (dowrite) then
                           write(unit8,'(2es18.10)') sqrt(dixx**2+diyy**2+dizz**2)*dsign(1.d0,dizz), phir(jj)
                       end if

                       ! Write along line in direction of the second closest atom
                       dowrite=gridpoint_close_to_straightline(ix, iy, iz, &
                           rxyz(1,iiat), rxyz(1,closeid(2)), hxh, hyh, hzh)
                       if (dowrite) then
                           write(unit9,'(2es18.10)') sqrt(dixx**2+diyy**2+dizz**2)*dsign(1.d0,dizz), phir(jj)
                       end if

                       ! Write along line in direction of the third closest atom
                       dowrite=gridpoint_close_to_straightline(ix, iy, iz, &
                           rxyz(1,iiat), rxyz(1,closeid(3)), hxh, hyh, hzh)
                       if (dowrite) then
                           write(unit10,'(2es18.10)') sqrt(dixx**2+diyy**2+dizz**2)*dsign(1.d0,dizz), phir(jj)
                       end if

                       ! Write along line in direction of the fourth closest atom
                       dowrite=gridpoint_close_to_straightline(ix, iy, iz, &
                           rxyz(1,iiat), rxyz(1,closeid(4)), hxh, hyh, hzh)
                       if (dowrite) then
                           write(unit11,'(2es18.10)') sqrt(dixx**2+diyy**2+dizz**2)*dsign(1.d0,dizz), phir(jj)
                       end if
                   end if

                end do
            end do
        end do

        ! Write the positions of the atoms, following the same order as above.
        ! For each grid point, write those atoms which lie in the plane
        ! perpendicular to the axis under consideration.

        if (plot_axis) then
            ! Along the x axis
            rxyzref(1,1)=rxyz(1,iiat)+1.d0 ; rxyzref(2,1)=rxyz(2,iiat) ; rxyzref(3,1)=rxyz(3,iiat)

            ! Along the y axis
            rxyzref(1,2)=rxyz(1,iiat) ; rxyzref(2,2)=rxyz(2,iiat)+1.d0 ; rxyzref(3,2)=rxyz(3,iiat)

            ! Along the z axis
            rxyzref(1,3)=rxyz(1,iiat) ; rxyzref(2,3)=rxyz(2,iiat) ; rxyzref(3,3)=rxyz(3,iiat)+1.d0
        end if

        if (plot_diagonals) then
            ! Along the diagonal in the octant +x,+y,+z
            rxyzref(1,4)=rxyz(1,iiat)+1.d0 ; rxyzref(2,4)=rxyz(2,iiat)+1.d0 ; rxyzref(3,4)=rxyz(3,iiat)+1.d0

            ! Along the diagonal in the octant -x,+y,+z
            rxyzref(1,5)=rxyz(1,iiat)-1.d0 ; rxyzref(2,5)=rxyz(2,iiat)+1.d0 ; rxyzref(3,5)=rxyz(3,iiat)+1.d0

            ! Along the diagonal in the octant -x,-y,+z
            rxyzref(1,6)=rxyz(1,iiat)-1.d0 ; rxyzref(2,6)=rxyz(2,iiat)-1.d0 ; rxyzref(3,6)=rxyz(3,iiat)+1.d0

            ! Along the diagonal in the octant +x,-y,+z
            rxyzref(1,7)=rxyz(1,iiat)+1.d0 ; rxyzref(2,7)=rxyz(2,iiat)-1.d0 ; rxyzref(3,7)=rxyz(3,iiat)+1.d0
        end if

        if (plot_neighbors) then
            ! Along the line in direction of the closest atom
            rxyzref(:,8)=rxyz(1,closeid(1))

            ! Along the line in direction of the second closest atom
            rxyzref(:,9)=rxyz(1,closeid(2))

            ! Along the line in direction of the third closest atom
            rxyzref(:,10)=rxyz(1,closeid(3))

            ! Along the line in direction of the fourth closest atom
            rxyzref(:,11)=rxyz(1,closeid(4))
        end if

        !!!write(unit12,'(a,es16.8)') '# sum(phir)', sum(phir)
        !!!write(unit12,'(a,2i7)') '# inwhichlocreg, onwhichatom', &
        !!!    tmb%orbs%inwhichlocreg(tmb%orbs%isorb+iorb), tmb%orbs%onwhichatom(tmb%orbs%isorb+iorb)
        !!!write(unit12,'(a,3i9,2x,2i9)') '# ix0, ixmin, ixmax, nsi1, n1i, ', &
        !!!    ix0, ixmin, ixmax, tmb%lzd%llr(ilr)%nsi1, tmb%lzd%llr(ilr)%d%n1i
        !!!write(unit12,'(a,3i9,2x,2i9)') '# iy0, iymin, iymax, nsi2, n2i, ', &
        !!!    iy0, iymin, iymax, tmb%lzd%llr(ilr)%nsi2, tmb%lzd%llr(ilr)%d%n2i
        !!!write(unit12,'(a,3i9,2x,2i9)') '# iz0, izmin, izmax, nsi3, n3i, ', &
        !!!    iz0, izmin, izmax, tmb%lzd%llr(ilr)%nsi3, tmb%lzd%llr(ilr)%d%n3i

        !!!do iat=1,11
        !!!    write(unit12,'(a,5(3es12.4,4x))') '#  ', &
        !!!        rxyz(:,iiat), rxyzref(:,iat), rxyz(:,iat), rxyzref(:,iat)-rxyz(:,iiat), rxyz(:,iat)-rxyz(:,iiat)
        !!!end do

        do iat=1,nat
            if (iat/=iiat) then
                write(unit12,'(12es12.3)') 0.d0, &
                                            base_point_distance(rxyz(:,iiat), rxyzref(:,1), rxyz(:,iat)), &
                                            base_point_distance(rxyz(:,iiat), rxyzref(:,2), rxyz(:,iat)), &
                                            base_point_distance(rxyz(:,iiat), rxyzref(:,3), rxyz(:,iat)), &
                                            base_point_distance(rxyz(:,iiat), rxyzref(:,4), rxyz(:,iat)), &
                                            base_point_distance(rxyz(:,iiat), rxyzref(:,5), rxyz(:,iat)), &
                                            base_point_distance(rxyz(:,iiat), rxyzref(:,6), rxyz(:,iat)), &
                                            base_point_distance(rxyz(:,iiat), rxyzref(:,7), rxyz(:,iat)), &
                                            base_point_distance(rxyz(:,iiat), rxyzref(:,8), rxyz(:,iat)), &
                                            base_point_distance(rxyz(:,iiat), rxyzref(:,9), rxyz(:,iat)), &
                                            base_point_distance(rxyz(:,iiat), rxyzref(:,10), rxyz(:,iat)), &
                                            base_point_distance(rxyz(:,iiat), rxyzref(:,11), rxyz(:,iat))
            else
                write(unit12,'(12es12.3)') 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0 , 0.d0
            end if
        end do


        if (plot_axis) then
            close(unit=unit1)
            close(unit=unit2)
            close(unit=unit3)
        end if
        if (plot_diagonals) then
            close(unit=unit4)
            close(unit=unit5)
            close(unit=unit6)
            close(unit=unit7)
        end if
        if (plot_neighbors) then
            close(unit=unit8)
            close(unit=unit9)
            close(unit=unit10)
            close(unit=unit11)
        end if
        close(unit=unit12)

        istart=istart+(tmb%lzd%glr%wfd%nvctr_c+7*tmb%lzd%glr%wfd%nvctr_f)*tmb%orbs%nspinor

    end do orbLoop

call deallocate_work_arrays_sumrho(w)
call f_free(phir)


contains

  function gridpoint_close_to_straightline(ix, iy, iz, a, b, hxh, hyh, hzh)
    ! Checks whether the grid point (ix,iy,iz) is close to the straight line
    ! going through the points a and b.
    !! IGNORE THIS "Close" means that the point is closest to the line in that plane which is
    !! "most orthogonal" (the angle between the line and the plane normal is
    !! minimal) to the line.
    
    ! Calling arguments
    integer,intent(in) :: ix, iy, iz
    real(kind=8),dimension(3),intent(in) :: a, b
    real(kind=8),intent(in) :: hxh, hyh, hzh
    logical :: gridpoint_close_to_straightline

    ! Local variables
    real(kind=8),dimension(3) :: rxyz
    real(kind=8) :: dist, threshold
    
    !!! Determine which plane is "most orthogonal" to the straight line
    !!xx(1)=1.d0 ; xx(2)=0.d0 ; xx(3)=0.d0
    !!yy(1)=0.d0 ; yy(2)=1.d0 ; yy(3)=0.d0
    !!zz(1)=0.d0 ; zz(2)=0.d0 ; zz(3)=1.d0
    !!bma=b-a
    !!abs_bma=dnrm2(3,bma,1)
    !!! angle between line and xy plane
    !!cosangle(1)=ddot(3,bma,1,zz,1)/dnrm2(bma)
    !!! angle between line and xz plane
    !!cosangle(2)=ddot(3,bma,1,yy,1)/dnrm2(bma)
    !!! angle between line and yz plane
    !!cosangle(3)=ddot(3,bma,1,xx,1)/dnrm2(bma)
    !!plane=minloc(cosangle)

    ! Calculate the shortest distance between the grid point (ix,iy,iz) and the
    ! straight line through the points a and b.
    rxyz = ix*hxh + iy*hyh + iz*hzh
    dist=get_distance(a, b, rxyz)

    ! Calculate the threshold
    threshold = sqrt(0.5d0*hxh**2 + 0.5d0*hyh**2 + 0.5d0*hzh**2)

    ! Check whether the point is close
    if (dist<threshold) then
        gridpoint_close_to_straightline=.true.
    else
        gridpoint_close_to_straightline=.false.
    end if

  end function gridpoint_close_to_straightline

  function get_distance(a, b, c)
    ! Calculate the shortest distance between point C and the 
    ! straight line trough the points A and B.

    ! Calling arguments
    real(kind=8),dimension(3),intent(in) :: a, b, c
    real(kind=8) :: get_distance

    ! Local variables
    real(kind=8),dimension(3) :: cma, bma, f, distvec
    real(kind=8) :: lambda, ddot, dnrm2

    cma=c-a 
    bma=b-a
    lambda=ddot(3,bma,1,cma,1)/ddot(3,bma,1,bma,1)
    
    ! The point on the straight line which is closest to c
    f=a+lambda*bma

    ! Get the distance between c and f
    distvec=c-f
    get_distance=dnrm2(3,distvec,1)

  end function get_distance


  function cross_product(a,b)
    ! Calculates the crosss product of the two vectors a and b.

    ! Calling arguments
    real(kind=8),dimension(3),intent(in) :: a, b
    real(kind=8),dimension(3) :: cross_product

    cross_product(1) = a(2)*b(3) - a(3)*b(2)
    cross_product(2) = a(3)*b(1) - a(1)*b(3)
    cross_product(3) = a(1)*b(2) - a(2)*b(1)

  end function cross_product



  function base_point_distance(a, b, c)
    ! Determine the base point of the perpendicular of the point C with respect
    ! to the vector going through the points A and B.

    ! Calling arguments
    real(kind=8),dimension(3),intent(in) :: a, b, c
    real(kind=8) :: base_point_distance

    ! Local variables
    real(kind=8),dimension(3) :: base_point, distance_vector, ab, ac
    real(kind=8) :: diffp1, diffm1, ddot, dnrm2, cosangle, lambda

    ! Vectors from A to B and from A to C
    ab = b - a
    ac = c - a

    !lambda = (ab(1)*ac(1) + ab(2)*ac(2) + ab(3)*ac(3)) / (ab(1)**2 + ab(2)**2 + ab(3)**2)
    lambda = ddot(3,ab,1,ac,1)/ddot(3,ab,1,ab,1)


    ! Base point of the perpendicular
    base_point = a + lambda*ab
    !base_point(1) = (a(1)*c(1)-b(1)*c(1))/(a(1)-b(1))
    !base_point(2) = (a(2)*c(2)-b(2)*c(2))/(a(2)-b(2))
    !base_point(3) = (a(3)*c(3)-b(3)*c(3))/(a(3)-b(3))


    ! Vector from the point A to the base point
    distance_vector = base_point - a

    ! Angle between the distance vector and vector from A to B.
    ! A cosine of 1 means that they are parallel, -1 means that they are anti-parallel.
    !if (dnrm2(3,distance_vector,1)*dnrm2(3,ab,1)==0.0d0) print*,'Error in plot orbitals',&
    !     dnrm2(3,distance_vector,1),dnrm2(3,ab,1),ddot(3,distance_vector,1,ab,1)
    if (ddot(3,distance_vector,1,ab,1)==0.0d0) then
       cosangle=0.0d0
    else
       cosangle = ddot(3,distance_vector,1,ab,1)/(dnrm2(3,distance_vector,1)*dnrm2(3,ab,1))
    end if
    diffp1=abs(cosangle-1)
    diffm1=abs(cosangle+1)
    if (diffp1<diffm1) then
        ! parallel
        base_point_distance = dnrm2(3,distance_vector,1)
    else
        ! antiparallel
        base_point_distance = -dnrm2(3,distance_vector,1)
    end if

  end function base_point_distance


end subroutine plotOrbitals


subroutine local_potential_dimensions(iproc,Lzd,orbs,xc,ndimfirstproc)
  use module_base
  use module_types
  use module_xc
  implicit none
  integer, intent(in) :: iproc, ndimfirstproc
  type(local_zone_descriptors), intent(inout) :: Lzd
  type(orbitals_data), intent(inout) :: orbs
  type(xc_info), intent(in) :: xc
  !local variables
  character(len=*), parameter :: subname='local_potential_dimensions'
  logical :: newvalue
  integer :: ii,iilr,ilr,iorb,iorb2,nilr,ispin
  integer, dimension(:,:), allocatable :: ilrtable

  call timing(iproc, 'calc_bounds   ', 'ON')
  call f_routine(id='local_potential_dimensions')
  
  if(Lzd%nlr > 1) then
     ilrtable = f_malloc((/ orbs%norbp, 2 /),id='ilrtable')
     ilrtable=0
     ii=0
     do iorb=1,orbs%norbp
        newvalue=.true.
        !localization region to which the orbital belongs
        ilr = orbs%inwhichlocreg(iorb+orbs%isorb)
        !spin state of the orbital
        if (orbs%spinsgn(orbs%isorb+iorb) > 0.0_gp) then
           ispin = 1       
        else
           ispin=2
        end if
        !check if the orbitals already visited have the same conditions
        loop_iorb2: do iorb2=1,orbs%norbp
           if(ilrtable(iorb2,1) == ilr .and. ilrtable(iorb2,2)==ispin) then
              newvalue=.false.
              exit loop_iorb2
           end if
        end do loop_iorb2
        if (newvalue) then
           ii = ii + 1
           ilrtable(ii,1)=ilr
           ilrtable(ii,2)=ispin    !SOMETHING IS NOT WORKING IN THE CONCEPT HERE... ispin is not a property of the locregs, but of the orbitals
        end if
     end do
     !number of inequivalent potential regions
     nilr = ii


     !calculate the dimension of the potential in the gathered form
     lzd%ndimpotisf=0
     do iilr=1,nilr
        ilr=ilrtable(iilr,1)
        do iorb=1,orbs%norbp
           !put the starting point
           if (orbs%inWhichLocreg(iorb+orbs%isorb) == ilr) then
              !assignment of ispot array to the value of the starting address of inequivalent
              orbs%ispot(iorb)=lzd%ndimpotisf + 1

           end if
        end do
        lzd%ndimpotisf = lzd%ndimpotisf + &
             lzd%llr(ilr)%d%n1i*lzd%llr(ilr)%d%n2i*lzd%llr(ilr)%d%n3i!*orbs%nspin
     end do
     !part which refers to exact exchange (only meaningful for one region)
     if (xc_exctXfac(xc) /= 0.0_gp) then
        lzd%ndimpotisf = lzd%ndimpotisf + &
             max(max(lzd%llr(ilr)%d%n1i*lzd%llr(ilr)%d%n2i*lzd%llr(ilr)%d%n3i*orbs%norbp,ndimfirstproc*orbs%norb),1)
     end if

  else 
     ilrtable = f_malloc((/ 1, 2 /),id='ilrtable')
     nilr = 1
     ilrtable=1

     !calculate the dimension of the potential in the gathered form
     lzd%ndimpotisf=0
     do iorb=1,orbs%norbp
        !assignment of ispot array to the value of the starting address of inequivalent
        orbs%ispot(iorb)=lzd%ndimpotisf + 1
        if(orbs%spinsgn(orbs%isorb+iorb) <= 0.0_gp) then
           orbs%ispot(iorb)=lzd%ndimpotisf + &
                1 + lzd%Glr%d%n1i*lzd%Glr%d%n2i*lzd%Glr%d%n3i
        end if
     end do
     lzd%ndimpotisf = lzd%ndimpotisf + &
          lzd%Glr%d%n1i*lzd%Glr%d%n2i*lzd%Glr%d%n3i*orbs%nspin
          
     !part which refers to exact exchange (only meaningful for one region)
     if (xc_exctXfac(xc) /= 0.0_gp) then
        lzd%ndimpotisf = lzd%ndimpotisf + &
             max(max(lzd%Glr%d%n1i*lzd%Glr%d%n2i*lzd%Glr%d%n3i*orbs%norbp,ndimfirstproc*orbs%norb),1)
     end if


  end if


  call f_free(ilrtable)

  call f_release_routine()
  call timing(iproc, 'calc_bounds   ', 'OF')

end subroutine local_potential_dimensions



subroutine print_orbital_distribution(iproc, nproc, orbs)
use module_base
use module_types
use yaml_output
implicit none

integer, intent(in) :: iproc, nproc
type(orbitals_data), intent(in) :: orbs

! Local variables
integer :: jproc, jpst, norb0,  norb1
!integer :: space1, space2, len1, len2, 
!logical :: written
logical, parameter :: print_all=.false.
integer :: minn, maxn
real(kind=8) :: avn

!!write(*,'(1x,a)') '------------------------------------------------------------------------------------'
!!written=.false.
!!write(*,'(1x,a)') '>>>> Partition of the basis functions among the processes.'
call yaml_map('Total No. Support Functions',orbs%norb,fmt='(i6)')
call yaml_mapping_open('Support Function Repartition')
jpst=0
avn=0.0d0
minn=orbs%norb
maxn=0
do jproc=0,nproc-1
    norb0=orbs%norb_par(jproc,0)
    norb1=orbs%norb_par(min(jproc+1,nproc-1),0)
    if (norb0/=norb1 .or. jproc==nproc-1) then
        if (print_all) then
            call yaml_map('MPI tasks '//trim(yaml_toa(jpst,fmt='(i0)'))//'-'//&
              &trim(yaml_toa(jproc,fmt='(i0)')),norb0,fmt='(i0)')
        end if
        minn=min(minn,norb0)
        maxn=max(maxn,norb0)
        jpst=jproc+1
    end if
    avn=avn+real(norb0,kind=8)
    !!if(orbs%norb_par(jproc,0)<orbs%norb_par(jproc-1,0)) then
    !!    len1=1+ceiling(log10(dble(jproc-1)+1.d-5))+ceiling(log10(dble(orbs%norb_par(jproc-1,0)+1.d-5)))
    !!    len2=ceiling(log10(dble(jproc)+1.d-5))+ceiling(log10(dble(nproc-1)+1.d-5))+&
    !!         ceiling(log10(dble(orbs%norb_par(jproc,0)+1.d-5)))
    !!    if(len1>=len2) then
    !!        space1=1
    !!        space2=1+len1-len2
    !!    else
    !!        space1=1+len2-len1
    !!        space2=1
    !!    end if
    !!    write(*,'(4x,a,2(i0,a),a,a)') '| Processes from 0 to ',jproc-1,' treat ',&
    !!        orbs%norb_par(jproc-1,0), ' orbitals,', repeat(' ', space1), '|'
    !!    write(*,'(4x,a,3(i0,a),a,a)')  '| processes from ',jproc,' to ',nproc-1,' treat ', &
    !!        orbs%norb_par(jproc,0),' orbitals.', repeat(' ', space2), '|'
    !!    written=.true.
    !!    exit
    !!end if
end do
avn=avn/real(nproc,kind=8)
call yaml_map('Minimum ',minn,fmt='(i0)')
call yaml_map('Maximum ',maxn,fmt='(i0)')
call yaml_map('Average ',avn,fmt='(f8.1)')
call yaml_mapping_close()
!!if(.not.written) then
!!    write(*,'(4x,a,2(i0,a),a,a)') '| Processes from 0 to ',nproc-1, &
!!        ' treat ',orbs%norbp,' orbitals. |'!, &
!!end if
!!write(*,'(1x,a)') '-----------------------------------------------'

!!written=.false.
!!write(*,'(1x,a)') '>>>> Partition of the basis functions including the derivatives among the processes.'
!!do jproc=1,nproc-1
!!    if(derorbs%norb_par(jproc,0)<derorbs%norb_par(jproc-1,0)) then
!!        len1=1+ceiling(log10(dble(jproc-1)+1.d-5))+ceiling(log10(dble(derorbs%norb_par(jproc-1,0)+1.d-5)))
!!        len2=ceiling(log10(dble(jproc)+1.d-5))+ceiling(log10(dble(nproc-1)+1.d-5))+&
!!             ceiling(log10(dble(derorbs%norb_par(jproc,0)+1.d-5)))
!!        if(len1>=len2) then
!!            space1=1
!!            space2=1+len1-len2
!!        else
!!            space1=1+len2-len1
!!            space2=1
!!        end if
!!        write(*,'(4x,a,2(i0,a),a,a)') '| Processes from 0 to ',jproc-1,' treat ',&
!!            derorbs%norb_par(jproc-1,0), ' orbitals,', repeat(' ', space1), '|'
!!        write(*,'(4x,a,3(i0,a),a,a)')  '| processes from ',jproc,' to ',nproc-1,' treat ', &
!!            derorbs%norb_par(jproc,0),' orbitals.', repeat(' ', space2), '|'
!!        written=.true.
!!        exit
!!    end if
!!end do
!!if(.not.written) then
!!    write(*,'(4x,a,2(i0,a),a,a)') '| Processes from 0 to ',nproc-1, &
!!        ' treat ',derorbs%norbp,' orbitals. |'
!!end if
!!write(*,'(1x,a)') '------------------------------------------------------------------------------------'


end subroutine print_orbital_distribution




subroutine cut_at_boundaries(cut, tmb)
  use module_base
  use module_types
  implicit none

  ! Calling arguments
  real(kind=8),intent(in)  :: cut
  type(DFT_wavefunction),intent(inout) :: tmb

  ! Local variables
  integer :: iorb, iiorb, ilr, icount, iseg, jj, j0, j1, ii, i3, i2, i0, i1, i, ishift, istart, iend
  real(kind=8) :: dist, cut2

  ! square of the cutoff radius
  cut2=cut**2

  ishift=0
  do iorb=1,tmb%orbs%norbp
      iiorb=tmb%orbs%isorb+iorb
      ilr=tmb%orbs%inwhichlocreg(iiorb)

      icount=0
      do iseg=1,tmb%lzd%llr(ilr)%wfd%nseg_c
         jj=tmb%lzd%llr(ilr)%wfd%keyvloc(iseg)
         j0=tmb%lzd%llr(ilr)%wfd%keygloc(1,iseg)
         j1=tmb%lzd%llr(ilr)%wfd%keygloc(2,iseg)
         ii=j0-1
         i3=ii/((tmb%lzd%llr(ilr)%d%n1+1)*(tmb%lzd%llr(ilr)%d%n2+1))
         ii=ii-i3*(tmb%lzd%llr(ilr)%d%n1+1)*(tmb%lzd%llr(ilr)%d%n2+1)
         i2=ii/(tmb%lzd%llr(ilr)%d%n1+1)
         i0=ii-i2*(tmb%lzd%llr(ilr)%d%n1+1)
         i1=i0+j1-j0
         do i=i0,i1
            dist = ((tmb%lzd%llr(ilr)%ns1+i )*tmb%lzd%hgrids(1)-tmb%lzd%llr(ilr)%locregcenter(1))**2 &
                 + ((tmb%lzd%llr(ilr)%ns2+i2)*tmb%lzd%hgrids(2)-tmb%lzd%llr(ilr)%locregcenter(2))**2 &
                 + ((tmb%lzd%llr(ilr)%ns3+i3)*tmb%lzd%hgrids(3)-tmb%lzd%llr(ilr)%locregcenter(3))**2
            if (dist>=cut2) then
                icount=icount+1
                tmb%psi(ishift+icount)=0.d0
            else
                icount=icount+1
            end if
         end do
      end do
      if (icount/=tmb%lzd%llr(ilr)%wfd%nvctr_c) then
          write(*,*) 'ERROR: icount /= tmb%lzd%llr(ilr)%wfd%nvctr_c', icount, tmb%lzd%llr(ilr)%wfd%nvctr_c
          stop
      end if
      ishift=ishift+tmb%lzd%llr(ilr)%wfd%nvctr_c

      ! fine part
      istart=tmb%lzd%llr(ilr)%wfd%nseg_c+(min(1,tmb%lzd%llr(ilr)%wfd%nseg_f))
      iend=tmb%lzd%llr(ilr)%wfd%nseg_c+tmb%lzd%llr(ilr)%wfd%nseg_f
      icount=0
      do iseg=istart,iend
         jj=tmb%lzd%llr(ilr)%wfd%keyvloc(iseg)
         j0=tmb%lzd%llr(ilr)%wfd%keygloc(1,iseg)
         j1=tmb%lzd%llr(ilr)%wfd%keygloc(2,iseg)
         ii=j0-1
         i3=ii/((tmb%lzd%llr(ilr)%d%n1+1)*(tmb%lzd%llr(ilr)%d%n2+1))
         ii=ii-i3*(tmb%lzd%llr(ilr)%d%n1+1)*(tmb%lzd%llr(ilr)%d%n2+1)
         i2=ii/(tmb%lzd%llr(ilr)%d%n1+1)
         i0=ii-i2*(tmb%lzd%llr(ilr)%d%n1+1)
         i1=i0+j1-j0
         do i=i0,i1
            dist = ((tmb%lzd%llr(ilr)%ns1+i )*tmb%lzd%hgrids(1)-tmb%lzd%llr(ilr)%locregcenter(1))**2 &
                 + ((tmb%lzd%llr(ilr)%ns2+i2)*tmb%lzd%hgrids(2)-tmb%lzd%llr(ilr)%locregcenter(2))**2 &
                 + ((tmb%lzd%llr(ilr)%ns3+i3)*tmb%lzd%hgrids(3)-tmb%lzd%llr(ilr)%locregcenter(3))**2
            if (dist>=cut2) then
                icount=icount+1 ; tmb%psi(ishift+icount)=0.d0
                icount=icount+1 ; tmb%psi(ishift+icount)=0.d0
                icount=icount+1 ; tmb%psi(ishift+icount)=0.d0
                icount=icount+1 ; tmb%psi(ishift+icount)=0.d0
                icount=icount+1 ; tmb%psi(ishift+icount)=0.d0
                icount=icount+1 ; tmb%psi(ishift+icount)=0.d0
                icount=icount+1 ; tmb%psi(ishift+icount)=0.d0
            else
                icount=icount+7
            end if
         end do
      end do
      if (icount/=7*tmb%lzd%llr(ilr)%wfd%nvctr_f) then
          write(*,*) 'ERROR: icount /= 7*tmb%lzd%llr(ilr)%wfd%nvctr_f', icount, 7*tmb%lzd%llr(ilr)%wfd%nvctr_f
          stop
      end if
      ishift=ishift+7*tmb%lzd%llr(ilr)%wfd%nvctr_f

  end do

  if (ishift/=tmb%npsidim_orbs) then
      write(*,*) 'ERROR: ishift /= tmb%npsidim_orbs', ishift, tmb%npsidim_orbs
      stop
  end if

end subroutine cut_at_boundaries





subroutine charge_center(n1i, n2i, n3i, hgrids, phir, charge_center_elec)
  implicit none
  ! Calling arguments
  integer,intent(in) :: n1i, n2i, n3i
  real(kind=8),dimension(3),intent(in) :: hgrids
  real(kind=8),dimension(n1i*n2i*n3i),intent(in) :: phir
  real(kind=8),dimension(3),intent(out) :: charge_center_elec

  integer :: i1, i2, i3, jj, iz, iy, ix, ii
  real(kind=8) :: q, x, y, z, qtot

  charge_center_elec=0.d0


  qtot=0.d0
  jj=0
  do i3=1,n3i
      do i2=1,n2i
          do i1=1,n1i
              jj=jj+1
              ! z component of point jj
              iz=jj/(n2i*n1i)
              ! Subtract the 'lower' xy layers
              ii=jj-iz*(n2i*n1i)
              ! y component of point jj
              iy=ii/n1i
              ! Subtract the 'lower' y rows
              ii=ii-iy*n1i
              ! x component
              ix=ii

              ! Shift the values due to the convolutions bounds
              ix=ix-14
              iy=iy-14
              iz=iz-14

              q= -phir(jj)**2 * product(hgrids)
              x=ix*hgrids(1)
              y=iy*hgrids(2)
              z=iz*hgrids(3)
              charge_center_elec(1) = charge_center_elec(1) + q*x
              charge_center_elec(2) = charge_center_elec(2) + q*y
              charge_center_elec(3) = charge_center_elec(3) + q*z
              qtot=qtot+q
          end do
      end do
  end do
  charge_center_elec(1:3)=charge_center_elec(1:3)/qtot


end subroutine charge_center









subroutine analyze_wavefunctions(output, region, lzd, orbs, npsidim, psi, ioffset)
  use module_base
  use module_types
  use yaml_output
  use bounds, only: locreg_bounds
  implicit none

  ! Calling arguments
  character(len=*),intent(in) :: output, region
  type(local_zone_descriptors),intent(inout) :: lzd
  type(orbitals_data),intent(in) :: orbs
  integer,intent(in) :: npsidim
  real(kind=8),dimension(npsidim),intent(in) :: psi
  integer,dimension(3,orbs%norbp),intent(in) :: ioffset
  
  ! Local variables
  integer :: ist, iorb, iiorb, ilr, ncount
  real(kind=8),dimension(3) :: center, sigma
  real(kind=8),dimension(:),allocatable :: sigma_arr
  real(kind=8) :: dnrm2


  if (trim(region)=='global') then
      ! Need to create the convolution bounds
      call locreg_bounds(lzd%glr%d%n1, lzd%glr%d%n2, lzd%glr%d%n3, &
           lzd%glr%d%nfl1, lzd%glr%d%nfu1, &
           lzd%glr%d%nfl2, lzd%glr%d%nfu2, &
           lzd%glr%d%nfl3, lzd%glr%d%nfu3, &
           lzd%glr%wfd, lzd%glr%bounds)
  end if

  sigma_arr = f_malloc0(orbs%norb,id='sigma_arr')

  ist = 1
  do iorb=1,orbs%norbp
      iiorb = orbs%isorb + iorb
      ilr = orbs%inwhichlocreg(iiorb)
      if (trim(region)=='local') then
          ncount = lzd%llr(ilr)%wfd%nvctr_c + 7*lzd%llr(ilr)%wfd%nvctr_f
          call analyze_one_wavefunction(lzd%llr(ilr), lzd%hgrids, ncount, psi(ist), ioffset(1,iorb), center, sigma)
      else if (trim(region)=='global') then
          ncount = lzd%glr%wfd%nvctr_c + 7*lzd%glr%wfd%nvctr_f
          call analyze_one_wavefunction(lzd%glr, lzd%hgrids, ncount, psi(ist), ioffset(1,iorb), center, sigma)
      else
          call f_err_throw('wrong value of region',err_name='BIGDFT_RUNTIME_ERROR')
      end if
      sigma_arr(iiorb) = dnrm2(3, sigma(1), 1)
      ist = ist + ncount
  end do

  call mpiallred(sigma_arr, mpi_sum, comm=bigdft_mpi%mpi_comm)

  if (trim(region)=='global') then
      !call deallocate_bounds(lzd%glr%geocode, lzd%glr%hybrid_on, lzd%glr%bounds)
  end if

  if (bigdft_mpi%iproc==0) then
      call yaml_sequence_open(trim(output),flow=.true.)
      call yaml_newline()
      do iorb=1,orbs%norb
          call yaml_sequence()
          call yaml_mapping_open(flow=.true.)
          call yaml_map('eval',orbs%eval(iorb),fmt='(es19.12)')
          call yaml_map('sigma',sigma_arr(iorb),fmt='(es11.4)')
          call yaml_mapping_close()
          call yaml_comment(yaml_toa(iorb))
          call yaml_newline()
      end do
      call yaml_sequence_close()
  end if

  call f_free(sigma_arr)

end subroutine analyze_wavefunctions


subroutine analyze_one_wavefunction(lr, hgrids, npsidim, psi, ioffset, center, sigma)
  use module_base
  use module_types
  use locreg_operations
  implicit none

  ! Calling arguments
  integer,intent(in) :: npsidim
  type(locreg_descriptors),intent(in) :: lr
  real(kind=8),dimension(3),intent(in) :: hgrids
  real(kind=8),dimension(npsidim),intent(in) :: psi
  integer,dimension(3),intent(in) :: ioffset
  real(kind=8),dimension(3),intent(out) :: center, sigma

  ! Local variables
  type(workarr_sumrho) :: w
  real(kind=8),dimension(:),allocatable :: psir
  integer :: ind, i1, i2, i3
  real(kind=8) :: x, y, z, q
  real(kind=8),dimension(3) :: hhgrids, var

  call initialize_work_arrays_sumrho(lr, .true., w)

  psir = f_malloc(lr%d%n1i*lr%d%n2i*lr%d%n3i,id='psir')
  ! Initialisation
  if (lr%geocode == 'F') call f_zero(psir)

  call daub_to_isf(lr, w, psi, psir)

  hhgrids(1:3) = 0.5d0*hgrids(1:3)

  ind = 0
  center(1:3) = 0.d0
  q = 0.d0
  do i3=1,lr%d%n3i
      z = real(i3+ioffset(3),wp)*hhgrids(3)
      do i2=1,lr%d%n2i
          y = real(i3+ioffset(2),wp)*hhgrids(2)
          do i1=1,lr%d%n1i
              x = real(i3+ioffset(1),wp)*hhgrids(1)
              ind = ind + 1
              center(1) = center(1) + psir(ind)**2*x
              center(2) = center(2) + psir(ind)**2*y
              center(3) = center(3) + psir(ind)**2*z
              q = q + psir(ind)**2
          end do
      end do
  end do
  ! Normalize
  center(1:3) = center(1:3)/q

  !Calculate variance
  ind = 0
  var(1:3) = 0.d0
  do i3=1,lr%d%n3i
      z = real(i3+ioffset(3),wp)*hhgrids(3)
      do i2=1,lr%d%n2i
          y = real(i3+ioffset(2),wp)*hhgrids(2)
          do i1=1,lr%d%n1i
              x = real(i3+ioffset(1),wp)*hhgrids(1)
              ind = ind + 1
              var(1) = var(1) + psir(ind)**2*(x-center(1))**2
              var(2) = var(2) + psir(ind)**2*(y-center(2))**2
              var(3) = var(3) + psir(ind)**2*(z-center(3))**2
          end do
      end do
  end do
  !Normalize
  var(1:3) = var(1:3)/q
  ! Take square root
  sigma(1) = sqrt(var(1))
  sigma(2) = sqrt(var(2))
  sigma(3) = sqrt(var(3))

  call f_free(psir)
  call deallocate_work_arrays_sumrho(w)
  !write(*,'(a,es16.8,5x,3(3es16.8,3x))') 'q, center(1:3), var(1:3), locregcenter',q, center(1:3), var(1:3), lr%locregcenter

end subroutine analyze_one_wavefunction


!> Use the (non-sparse) coefficients to calculate a non-sparse kernel, then
!! analyze the magnitude of the elements.
!! WARNING: This routine must not be called in parallel
!! WARNING: This routine is not tested with spin polarization
subroutine analyze_kernel(ntmb, norb, nat, coeff, kernel, rxyz, on_which_atom)
  use module_base
  use module_types
  use sparsematrix_base, only: matrices, matrices_null, allocate_matrices, &
                               deallocate_matrices
  use yaml_output
  implicit none
  ! Calling arguments
  integer,intent(in) :: ntmb, norb, nat
  real(kind=8),dimension(ntmb,norb),intent(in) :: coeff
  real(kind=8),dimension(ntmb,ntmb),intent(in) :: kernel
  real(kind=8),dimension(3,nat),intent(in) :: rxyz
  integer,dimension(ntmb),intent(in) :: on_which_atom

  ! Local variables
  integer :: iorb, itmb, jtmb, iat, jat, iunit, nproc, ierr
  logical :: mpi_initd
  real(kind=8) :: d, asymm, maxdiff, diff
  real(kind=8),dimension(3) :: dist
  real(kind=8),dimension(:,:),allocatable :: kernel_full
  character(len=*),parameter :: filename='kernel_analysis.dat'
  real(kind=8),parameter :: print_limit=1.d-8

  call f_routine(id='analyze_kernel')

  ! Check that this is a monoproc run
  call mpi_initialized(mpi_initd, ierr)
  if (mpi_initd) then
      call mpi_comm_size(mpi_comm_world, nproc, ierr)
  else
      nproc = 1
  end if

  write(*,*) 'nproc', nproc
  if (nproc/=1) then
      call f_err_throw('analyze_kernel should only be called using 1 MPI task',err_name='BIGDT_RUNTIME_ERROR')
  end if

  kernel_full = f_malloc0((/ntmb,ntmb/),id='kernel_full')

  do iorb=1,norb
      call yaml_map('orbital being processed',iorb)
      call gemm('n', 't', ntmb, ntmb, 1, 1.d0, coeff(1,iorb), ntmb, &
           coeff(1,iorb), ntmb, 1.d0, kernel_full(1,1), ntmb)
  end do
  !call mpiallred(kernel%matrix, mpi_sum, comm=bigdft_mpi%mpi_comm)

  call yaml_map('Output file for kernel analysis',trim(filename))
  call f_open_file(iunit, file=trim(filename), binary=.false.)
  write(iunit,'(a)') '#     itmb,   jtmb,                d,              val'
  maxdiff = 0.d0
  do itmb=1,ntmb
      call yaml_map('basis function being processed',itmb)
      iat = on_which_atom(itmb)
      do jtmb=1,itmb
          jat = on_which_atom(jtmb)
          dist(1:3) = rxyz(1:3,iat)-rxyz(1:3,jat)
          d = nrm2(3, dist(1), 1)
          asymm = abs(kernel_full(jtmb,itmb)-kernel_full(itmb,jtmb))
          if (asymm>1.d-15) then
              call yaml_warning('kernel not symmetric, diff='//yaml_toa(asymm,fmt='(es9.2)'))
          end if
          diff = abs(kernel_full(jtmb,itmb)-kernel(jtmb,itmb))
          maxdiff = max(diff,maxdiff)
          if (abs(kernel_full(jtmb,itmb))>print_limit) then
              write(iunit,'(2x,2i8,2es18.10)') itmb, jtmb, d, kernel_full(jtmb,itmb)
          end if
      end do
  end do
  call yaml_map('maxdiff of sparse and full kernel',maxdiff)
  call f_close(iunit)

  call f_free(kernel_full)

  call f_release_routine()

end subroutine analyze_kernel
