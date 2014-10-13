!> @file 
!!   Miscellaneous routines for linear toolbox
!! @author
!!   Copyright (C) 2011-2013 BigDFT group 
!!   This file is distributed under the terms of the
!!   GNU General Public License, see ~/COPYING file
!!   or http://www.gnu.org/copyleft/gpl.txt .
!!   For the list of contributors, see ~/AUTHORS 
 

!> Plots the orbitals
subroutine plotOrbitals(iproc, tmb, phi, nat, rxyz, hxh, hyh, hzh, it, basename)
   use module_base
   use module_types
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

   call initialize_work_arrays_sumrho(1,tmb%lzd%glr,.true.,w)
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
        call to_zero(tmb%lzd%glr%d%n1i*tmb%lzd%glr%d%n2i*tmb%lzd%glr%d%n3i, phir(1))
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



subroutine plotGrid(iproc, norb, nspinor, nspin, orbitalNumber, llr, glr, atoms, rxyz, hx, hy, hz)
  use module_base
  use module_types
  implicit none
  
  ! Calling arguments
  integer, intent(in) :: iproc, norb, nspinor, nspin, orbitalNumber
  type(locreg_descriptors), intent(in) :: llr, glr
  type(atoms_data), intent(in) ::atoms
  real(kind=8), dimension(3,atoms%astruct%nat), intent(in) :: rxyz
  real(kind=8), intent(in) :: hx, hy, hz
  
  ! Local variables
  integer :: iseg, jj, j0, j1, ii, i3, i2, i0, i1, i, ishift, iat, ldim, gdim, jjj
  character(len=10) :: num
  character(len=20) :: filename
  real(kind=8), dimension(:), allocatable :: lphi, phi


    ldim=llr%wfd%nvctr_c+7*llr%wfd%nvctr_f
    gdim=glr%wfd%nvctr_c+7*glr%wfd%nvctr_f
    lphi = f_malloc(ldim,id='lphi')
    phi = f_malloc(gdim,id='phi')
    lphi=1.d0
    !!phi=0.d0
    call to_zero(gdim, phi(1))
    call Lpsi_to_global2(iproc, ldim, gdim, norb, nspinor, nspin, glr, llr, lphi, phi)
  
    write(num,'(i0)') orbitalNumber
    filename='orbital_'//trim(num)
  
    open(unit=2000+iproc,file=trim(filename)//'.xyz',status='unknown')
    !write(2000+iproc,*) llr%wfd%nvctr_c+llr%wfd%nvctr_f+atoms%astruct%nat,' atomic'
    write(2000+iproc,*) glr%wfd%nvctr_c+glr%wfd%nvctr_f+llr%wfd%nvctr_c+llr%wfd%nvctr_f+atoms%astruct%nat,' atomic'
    if (atoms%astruct%geocode=='F') then
       write(2000+iproc,*)'complete simulation grid with low and high resolution points'
    else if (atoms%astruct%geocode =='S') then
       write(2000+iproc,'(a,2x,3(1x,1pe24.17))')'surface',atoms%astruct%cell_dim(1),atoms%astruct%cell_dim(2),&
            atoms%astruct%cell_dim(3)
    else if (atoms%astruct%geocode =='P') then
       write(2000+iproc,'(a,2x,3(1x,1pe24.17))')'periodic',atoms%astruct%cell_dim(1),atoms%astruct%cell_dim(2),&
            atoms%astruct%cell_dim(3)
    end if

   do iat=1,atoms%astruct%nat
      write(2000+iproc,'(a6,2x,3(1x,e12.5),3x)') trim(atoms%astruct%atomnames(atoms%astruct%iatype(iat))),&
           rxyz(1,iat),rxyz(2,iat),rxyz(3,iat)
   end do

  
    jjj=0
    do iseg=1,glr%wfd%nseg_c
       jj=glr%wfd%keyvloc(iseg)
       j0=glr%wfd%keygloc(1,iseg)
       j1=glr%wfd%keygloc(2,iseg)
       ii=j0-1
       i3=ii/((glr%d%n1+1)*(glr%d%n2+1))
       ii=ii-i3*(glr%d%n1+1)*(glr%d%n2+1)
       i2=ii/(glr%d%n1+1)
       i0=ii-i2*(glr%d%n1+1)
       i1=i0+j1-j0
       do i=i0,i1
           jjj=jjj+1
           if(phi(jjj)==1.d0) write(2000+iproc,'(a4,2x,3(1x,e10.3))') '  lg ',&
                real(i,kind=8)*hx,real(i2,kind=8)*hy,real(i3,kind=8)*hz
           write(2000+iproc,'(a4,2x,3(1x,e10.3))') '  g ',real(i,kind=8)*hx,&
                real(i2,kind=8)*hy,real(i3,kind=8)*hz
       enddo
    enddo

    ishift=glr%wfd%nseg_c  
    ! fine part
    do iseg=1,glr%wfd%nseg_f
       jj=glr%wfd%keyvloc(ishift+iseg)
       j0=glr%wfd%keygloc(1,ishift+iseg)
       j1=glr%wfd%keygloc(2,ishift+iseg)
       ii=j0-1
       i3=ii/((glr%d%n1+1)*(glr%d%n2+1))
       ii=ii-i3*(glr%d%n1+1)*(glr%d%n2+1)
       i2=ii/(glr%d%n1+1)
       i0=ii-i2*(glr%d%n1+1)
       i1=i0+j1-j0
       do i=i0,i1
          jjj=jjj+1
          if(phi(jjj)==1.d0) write(2000+iproc,'(a4,2x,3(1x,e10.3))') '  lG ',real(i,kind=8)*hx,real(i2,kind=8)*hy,real(i3,kind=8)*hz
          write(2000+iproc,'(a4,2x,3(1x,e10.3))') '  G ',real(i,kind=8)*hx,real(i2,kind=8)*hy,real(i3,kind=8)*hz
          jjj=jjj+6
       enddo
    enddo
  
    close(unit=2000+iproc)

end subroutine plotGrid


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

!!write(*,'(1x,a)') '------------------------------------------------------------------------------------'
!!written=.false.
!!write(*,'(1x,a)') '>>>> Partition of the basis functions among the processes.'
call yaml_mapping_open('Support function repartition')
jpst=0
do jproc=0,nproc-1
    norb0=orbs%norb_par(jproc,0)
    norb1=orbs%norb_par(min(jproc+1,nproc-1),0)
    if (norb0/=norb1 .or. jproc==nproc-1) then
        call yaml_map('MPI tasks '//trim(yaml_toa(jpst,fmt='(i0)'))//'-'//trim(yaml_toa(jproc,fmt='(i0)')),norb0,fmt='(i0)')
        jpst=jproc+1
    end if

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



subroutine build_ks_orbitals(iproc, nproc, tmb, KSwfn, at, rxyz, denspot, GPU, &
           energs, nlpsp, input, order_taylor, &
           energy, energyDiff, energyold)
  use module_base
  use module_types
  use module_interfaces, except_this_one => build_ks_orbitals
  use communications_base, only: comms_cubic
  use communications_init, only: orbitals_communicators
  use communications, only: transpose_v, untranspose_v
  use sparsematrix_base, only: sparse_matrix
  use yaml_output
  implicit none
  
  ! Calling arguments
  integer:: iproc, nproc
  type(DFT_wavefunction),intent(inout) :: tmb, KSwfn
  type(atoms_data), intent(in) :: at
  real(gp), dimension(3,at%astruct%nat), intent(in) :: rxyz
  type(DFT_local_fields), intent(inout) :: denspot
  type(GPU_pointers), intent(inout) :: GPU
  type(energy_terms),intent(inout) :: energs
  type(DFT_PSP_projectors), intent(inout) :: nlpsp
  type(input_variables),intent(in) :: input
  integer,intent(inout) :: order_taylor
  real(kind=8),intent(out) :: energy, energyDiff
  real(kind=8), intent(inout) :: energyold

  ! Local variables
  type(orbitals_data) :: orbs
  type(comms_cubic) :: comms
  real(gp) :: fnrm
  logical :: rho_negative
  integer :: infoCoeff, nvctrp, npsidim_global
  real(kind=8),dimension(:),pointer :: phi_global, phiwork_global
  character(len=*),parameter :: subname='build_ks_orbitals'
  real(wp), dimension(:,:,:), pointer :: mom_vec_fake


  nullify(mom_vec_fake)

  !debug
  !integer :: iorb, jorb, ist, jst, ierr, i
  !real(kind=8) :: ddot, tt


  ! Get the expansion coefficients
  ! Start with a "clean" density, i.e. without legacy from previous mixing steps
  call communicate_basis_for_density_collective(iproc, nproc, tmb%lzd, &
       max(tmb%npsidim_orbs,tmb%npsidim_comp), tmb%orbs, tmb%psi, tmb%collcom_sr)
  call sumrho_for_TMBs(iproc, nproc, KSwfn%Lzd%hgrids(1), KSwfn%Lzd%hgrids(2), KSwfn%Lzd%hgrids(3), &
       tmb%collcom_sr, tmb%linmat%l, tmb%linmat%kernel_, denspot%dpbox%ndimrhopot, &
       denspot%rhov, rho_negative)
  if (rho_negative) then
      call corrections_for_negative_charge(iproc, nproc, KSwfn, at, input, tmb, denspot)
      !!if (iproc==0) call yaml_warning('Charge density contains negative points, need to increase FOE cutoff')
      !!call increase_FOE_cutoff(iproc, nproc, tmb%lzd, at%astruct, input, KSwfn%orbs, tmb%orbs, tmb%foe_obj, init=.false.)
      !!call clean_rho(iproc, nproc, KSwfn%Lzd%Glr%d%n1i*KSwfn%Lzd%Glr%d%n2i*denspot%dpbox%n3d, denspot%rhov)
  end if

  call updatePotential(input%nspin,denspot,energs%eh,energs%exc,energs%evxc)

  tmb%can_use_transposed=.false.
  call get_coeff(iproc, nproc, LINEAR_MIXDENS_SIMPLE, KSwfn%orbs, at, rxyz, denspot, GPU, infoCoeff, &
       energs, nlpsp, input%SIC, tmb, fnrm, .true., .false., .true., 0, 0, 0, 0, &
       order_taylor,input%lin%max_inversion_error,input%purification_quickreturn,&
       input%calculate_KS_residue,input%calculate_gap)

  if (bigdft_mpi%iproc ==0) then
     call write_eigenvalues_data(0.1d0,KSwfn%orbs,mom_vec_fake)
  end if


  !call communicate_basis_for_density_collective(iproc, nproc, tmb%lzd, &
  !     max(tmb%npsidim_orbs,tmb%npsidim_comp), tmb%orbs, tmb%psi, tmb%collcom_sr)
  !call sumrho_for_TMBs(iproc, nproc, KSwfn%Lzd%hgrids(1), KSwfn%Lzd%hgrids(2), KSwfn%Lzd%hgrids(3), &
  !     tmb%collcom_sr, tmb%linmat%denskern, KSwfn%Lzd%Glr%d%n1i*KSwfn%Lzd%Glr%d%n2i*denspot%dpbox%n3d, denspot%rhov)
  !call updatePotential(input%ixc,input%nspin,denspot,energs%eh,energs%exc,energs%evxc)
  !tmb%can_use_transposed=.false.
  !call get_coeff(iproc, nproc, LINEAR_MIXDENS_SIMPLE, KSwfn%orbs, at, rxyz, denspot, GPU, infoCoeff, &
  !     energs, nlpspd, proj, input%SIC, tmb, fnrm, .true., .false., .true., ham_small, 0, 0, 0, 0)
  !energy=energs%ebs-energs%eh+energs%exc-energs%evxc-energs%eexctX+energs%eion+energs%edisp
  !energyDiff=energy-energyold
  !energyold=energy

  !!if(tmb%can_use_transposed) then
  !!    call f_free_ptr(tmb%psit_c)
  !!    call f_free_ptr(tmb%psit_f)
  !!end if

  ! Create communication arrays for support functions in the global box
  
  call nullify_orbitals_data(orbs)
  call copy_orbitals_data(tmb%orbs, orbs, subname)
  call orbitals_communicators(iproc, nproc, tmb%lzd%glr, orbs, comms)


  ! Transform the support functions to the global box
  ! WARNING: WILL NOT WORK WITH K-POINTS, CHECK THIS
  npsidim_global=max(tmb%orbs%norbp*(tmb%lzd%glr%wfd%nvctr_c+7*tmb%lzd%glr%wfd%nvctr_f), &
                     tmb%orbs%norb*comms%nvctr_par(iproc,0)*orbs%nspinor)
  phi_global = f_malloc_ptr(npsidim_global,id='phi_global')
  phiwork_global = f_malloc_ptr(npsidim_global,id='phiwork_global')
  call small_to_large_locreg(iproc, tmb%npsidim_orbs, &
       tmb%orbs%norbp*(tmb%lzd%glr%wfd%nvctr_c+7*tmb%lzd%glr%wfd%nvctr_f), tmb%lzd, &
       KSwfn%lzd, tmb%orbs, tmb%psi, phi_global, to_global=.true.)
  call transpose_v(iproc, nproc, orbs, tmb%lzd%glr%wfd, comms, phi_global(1), phiwork_global(1))


  ! WARNING: WILL NOT WORK WITH K-POINTS, CHECK THIS
  nvctrp=comms%nvctr_par(iproc,0)*orbs%nspinor
  call dgemm('n', 'n', nvctrp, KSwfn%orbs%norb, tmb%orbs%norb, 1.d0, phi_global, nvctrp, tmb%coeff(1,1), &
             tmb%orbs%norb, 0.d0, phiwork_global, nvctrp)
  
  call untranspose_v(iproc, nproc, KSwfn%orbs, tmb%lzd%glr%wfd, KSwfn%comms, phiwork_global(1), phi_global(1))  

  call f_free_ptr(phi_global)

  !!ist=1
  !!do iorb=1,KSwfn%orbs%norbp
  !!    do i=1,tmb%lzd%glr%wfd%nvctr_c+7*tmb%lzd%glr%wfd%nvctr_f
  !!        write(800+iproc,*) iorb, i, phiwork_global(ist)
  !!        ist=ist+1
  !!    end do
  !!end do


  !!ierr=tmb%lzd%glr%wfd%nvctr_c+7*tmb%lzd%glr%wfd%nvctr_f
  !!do i=1,KSwfn%orbs%norb*ierr
  !!    write(401,*) i, phiwork_global(i)
  !!end do
  !!write(*,*) 'GLOBAL DDOT',ddot(KSwfn%orbs%norb*ierr, phi_global, 1, phi_global, 1)

  !!do i=1,KSwfn%orbs%norb*ierr
  !!     tmb%psi(i)=phi_global(i)
  !!end do
  !!call get_coeff(iproc, nproc, LINEAR_MIXDENS_SIMPLE, tmb%orbs, at, rxyz, denspot, GPU, infoCoeff, &
  !!     energs, nlpspd, proj, input%SIC, tmb, fnrm, .true., .false., .true., ham_small, 0, 0, 0, 0)

  !!do i=1,KSwfn%orbs%norb*(tmb%lzd%glr%wfd%nvctr_c+7*tmb%lzd%glr%wfd%nvctr_f)
  !!    write(600,'(i10,es16.7)') i, tmb%psi(i)
  !!end do


  !!write(*,*) 'iproc, input%output_wf_format',iproc, WF_FORMAT_PLAIN
  call writemywaves(iproc,trim(input%dir_output)//"wavefunction", WF_FORMAT_PLAIN, &
       KSwfn%orbs, KSwfn%Lzd%Glr%d%n1, KSwfn%Lzd%Glr%d%n2, KSwfn%Lzd%Glr%d%n3, &
       KSwfn%Lzd%hgrids(1), KSwfn%Lzd%hgrids(2), KSwfn%Lzd%hgrids(3), &
       at, rxyz, KSwfn%Lzd%Glr%wfd, phiwork_global)

   call f_free_ptr(phiwork_global)
   call deallocate_orbitals_data(orbs)
   call deallocate_comms_cubic(comms)

  ! To get consistent values of the energy and the Kohn-Sham residue with those
  ! which will be calculated by the cubic restart.
  call communicate_basis_for_density_collective(iproc, nproc, tmb%lzd, &
       max(tmb%npsidim_orbs,tmb%npsidim_comp), tmb%orbs, tmb%psi, tmb%collcom_sr)
  call sumrho_for_TMBs(iproc, nproc, KSwfn%Lzd%hgrids(1), KSwfn%Lzd%hgrids(2), KSwfn%Lzd%hgrids(3), &
       tmb%collcom_sr, tmb%linmat%l, tmb%linmat%kernel_, denspot%dpbox%ndimrhopot, &
       denspot%rhov, rho_negative)
  if (rho_negative) then
      call corrections_for_negative_charge(iproc, nproc, KSwfn, at, input, tmb, denspot)
      !!if (iproc==0) call yaml_warning('Charge density contains negative points, need to increase FOE cutoff')
      !!call increase_FOE_cutoff(iproc, nproc, tmb%lzd, at%astruct, input, KSwfn%orbs, tmb%orbs, tmb%foe_obj, init=.false.)
      !!call clean_rho(iproc, nproc, KSwfn%Lzd%Glr%d%n1i*KSwfn%Lzd%Glr%d%n2i*denspot%dpbox%n3d, denspot%rhov)
  end if
  call updatePotential(input%nspin,denspot,energs%eh,energs%exc,energs%evxc)
  tmb%can_use_transposed=.false.
  call get_coeff(iproc, nproc, LINEAR_MIXDENS_SIMPLE, KSwfn%orbs, at, rxyz, denspot, GPU, infoCoeff, &
       energs, nlpsp, input%SIC, tmb, fnrm, .true., .false., .true., 0, 0, 0, 0, &
       order_taylor, input%lin%max_inversion_error, input%purification_quickreturn, &
       input%calculate_KS_residue, input%calculate_gap, updatekernel=.false.)
  energy=energs%ebs-energs%eh+energs%exc-energs%evxc-energs%eexctX+energs%eion+energs%edisp
  energyDiff=energy-energyold
  energyold=energy

  !!if(tmb%can_use_transposed) then
  !!    call f_free_ptr(tmb%psit_c)
  !!    call f_free_ptr(tmb%psit_f)
  !!end if

end subroutine build_ks_orbitals



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



subroutine loewdin_charge_analysis(iproc,tmb,atoms,denspot,&
           calculate_overlap_matrix,calculate_ovrlp_half,meth_overlap)
  use module_base
  use module_types
  use module_interfaces, except_this_one => loewdin_charge_analysis
  use communications, only: transpose_localized
  use sparsematrix_base, only: sparse_matrix, sparsematrix_malloc, sparsematrix_malloc0, sparsematrix_malloc_ptr, &
                               DENSE_FULL, assignment(=), &
                               matrices_null, allocate_matrices, deallocate_matrices
  use sparsematrix, only: compress_matrix, uncompress_matrix
  use yaml_output
  implicit none
  integer,intent(in) :: iproc
  type(dft_wavefunction),intent(inout) :: tmb
  type(atoms_data),intent(in) :: atoms
  type(DFT_local_fields), intent(inout) :: denspot
  logical,intent(in) :: calculate_overlap_matrix, calculate_ovrlp_half
  integer,intent(in) :: meth_overlap

  !local variables
  !integer :: ifrag,ifrag_ref,isforb,jorb
  integer :: iorb,ierr
  real(kind=gp), allocatable, dimension(:,:,:) :: proj_mat
  real(kind=gp), allocatable, dimension(:,:) :: proj_ovrlp_half, weight_matrixp
  character(len=*),parameter :: subname='calculate_weight_matrix_lowdin'
  real(kind=gp) :: max_error, mean_error
  type(matrices),dimension(1) :: inv_ovrlp

  ! new variables
  integer :: iat
  real(kind=8),dimension(:,:),allocatable :: weight_matrix
  !real(kind=gp),dimension(:,:),pointer :: ovrlp
  real(kind=8) :: total_charge, total_net_charge
  real(kind=8),dimension(:),allocatable :: charge_per_atom
  !logical :: psit_c_associated, psit_f_associated


  ! needs parallelizing/converting to sparse
  ! re-use overlap matrix if possible either before or after

  call f_routine(id='loewdin_charge_analysis')

  inv_ovrlp(1) = matrices_null()
  call allocate_matrices(tmb%linmat%l, allocate_full=.true., matname='inv_ovrlp', mat=inv_ovrlp(1))



  if (calculate_overlap_matrix) then
     if(.not.tmb%can_use_transposed) then
         !!if(.not.associated(tmb%psit_c)) then
         !!    tmb%psit_c = f_malloc_ptr(sum(tmb%collcom%nrecvcounts_c),id='tmb%psit_c')
         !!    psit_c_associated=.false.
         !!else
         !!    psit_c_associated=.true.
         !!end if
         !!if(.not.associated(tmb%psit_f)) then
         !!    tmb%psit_f = f_malloc_ptr(7*sum(tmb%collcom%nrecvcounts_f),id='tmb%psit_f')
         !!    psit_f_associated=.false.
         !!else
         !!    psit_f_associated=.true.
         !!end if
         call transpose_localized(bigdft_mpi%iproc, bigdft_mpi%nproc, tmb%npsidim_orbs, tmb%orbs, tmb%collcom, &
              tmb%psi, tmb%psit_c, tmb%psit_f, tmb%lzd)
         tmb%can_use_transposed=.true.
     end if

     call calculate_overlap_transposed(bigdft_mpi%iproc, bigdft_mpi%nproc, tmb%orbs, tmb%collcom, tmb%psit_c, &
          tmb%psit_c, tmb%psit_f, tmb%psit_f, tmb%linmat%s, tmb%linmat%ovrlp_)
     ! This can then be deleted if the transition to the new type has been completed.
     !tmb%linmat%ovrlp%matrix_compr=tmb%linmat%ovrlp_%matrix_compr


     !!if (.not.psit_c_associated) then
     !!   call f_free_ptr(tmb%psit_c)
     !!   tmb%can_use_transposed=.false.
     !!end if
     !!if (.not.psit_f_associated) then
     !!   call f_free_ptr(tmb%psit_f)
     !!   tmb%can_use_transposed=.false.
     !!end if
  end if

  if (calculate_ovrlp_half) then
     tmb%linmat%ovrlp_%matrix = sparsematrix_malloc_ptr(tmb%linmat%s, iaction=DENSE_FULL, id='tmb%linmat%ovrlp_%matrix')
     call uncompress_matrix(bigdft_mpi%iproc, tmb%linmat%s, &
          inmat=tmb%linmat%ovrlp_%matrix_compr, outmat=tmb%linmat%ovrlp_%matrix)
     call overlapPowerGeneral(bigdft_mpi%iproc, bigdft_mpi%nproc, meth_overlap, 1, (/2/), &
          tmb%orthpar%blocksize_pdsyev, &
          imode=2, ovrlp_smat=tmb%linmat%s, inv_ovrlp_smat=tmb%linmat%l, &
          ovrlp_mat=tmb%linmat%ovrlp_, inv_ovrlp_mat=inv_ovrlp, check_accur=.true., &
          max_error=max_error, mean_error=mean_error)
     !!ovrlp_half=tmb%linmat%ovrlp%matrix
     call f_free_ptr(tmb%linmat%ovrlp_%matrix)
  end if

  ! optimize this to just change the matrix multiplication?
  proj_mat = sparsematrix_malloc0(tmb%linmat%l,iaction=DENSE_FULL,id='proj_mat')

  call uncompress_matrix(iproc, tmb%linmat%l, inmat=tmb%linmat%kernel_%matrix_compr, outmat=proj_mat)
  !!isforb=0
  !!do ifrag=1,input%frag%nfrag
  !!   ifrag_ref=input%frag%frag_index(ifrag)
  !!   if (ifrag==ifrag_charged(1)) then
  !!      do iorb=1,ref_frags(ifrag_ref)%fbasis%forbs%norb
  !!         proj_mat(iorb+isforb,iorb+isforb)=1.0_gp
  !!      end do
  !!   end if
  !!   !!if (nfrag_charged==2) then
  !!   !!   if (ifrag==ifrag_charged(2)) then
  !!   !!      do iorb=1,ref_frags(ifrag_ref)%fbasis%forbs%norb
  !!   !!         proj_mat(iorb+isforb,iorb+isforb)=-1.0_gp
  !!   !!      end do
  !!   !!   end if
  !!   !!end if
  !!   isforb=isforb+ref_frags(ifrag_ref)%fbasis%forbs%norb
  !!end do

  proj_ovrlp_half=f_malloc((/tmb%orbs%norb,tmb%orbs%norbp/),id='proj_ovrlp_half')
  if (tmb%orbs%norbp>0) then
     call dgemm('n', 'n', tmb%orbs%norb, tmb%orbs%norbp, &
            tmb%orbs%norb, 1.d0, &
            proj_mat(1,1,1), tmb%orbs%norb, &
            inv_ovrlp(1)%matrix(1,tmb%orbs%isorb+1,1), tmb%orbs%norb, 0.d0, &
            proj_ovrlp_half(1,1), tmb%orbs%norb)
  end if
  call f_free(proj_mat)
  weight_matrixp=f_malloc((/tmb%orbs%norb,tmb%orbs%norbp/), id='weight_matrixp')
  if (tmb%orbs%norbp>0) then
     call dgemm('n', 'n', tmb%orbs%norb, tmb%orbs%norbp, &
          tmb%orbs%norb, 1.d0, &
          inv_ovrlp(1)%matrix(1,1,1), tmb%orbs%norb, &
          proj_ovrlp_half(1,1), tmb%orbs%norb, 0.d0, &
          weight_matrixp(1,1), tmb%orbs%norb)
  end if
  !call f_free_ptr(ovrlp_half)
  call f_free(proj_ovrlp_half)
  weight_matrix=f_malloc((/tmb%orbs%norb,tmb%orbs%norb/), id='weight_matrix')
  if (bigdft_mpi%nproc>1) then
     call mpi_allgatherv(weight_matrixp, tmb%orbs%norb*tmb%orbs%norbp, mpi_double_precision, weight_matrix, &
          tmb%orbs%norb*tmb%orbs%norb_par(:,0), tmb%orbs%norb*tmb%orbs%isorb_par, &
          mpi_double_precision, bigdft_mpi%mpi_comm, ierr)
  else
     call vcopy(tmb%orbs%norb*tmb%orbs%norb,weight_matrixp(1,1),1,weight_matrix(1,1),1)
  end if
  call f_free(weight_matrixp)
  !call compress_matrix(bigdft_mpi%iproc,weight_matrix)

  charge_per_atom = f_malloc0(atoms%astruct%nat,id='charge_per_atom')
  !!do iorb=1,tmb%orbs%norb
  !!    do jorb=1,tmb%orbs%norb
  !!        if (iproc==0) write(*,'(a,2i7,es16.7)') 'iorb,jorb,weight_matrix(jorb,iorb)', iorb,jorb,weight_matrix(jorb,iorb)
  !!        if (iorb==jorb) then
  !!            total_charge = total_charge + weight_matrix(jorb,iorb)
  !!            iat=tmb%orbs%onwhichatom(iorb)
  !!            charge_per_atom(iat) = charge_per_atom(iat) + weight_matrix(jorb,iorb)
  !!        end if
  !!    end do
  !!end do
  !!if (iproc==0) then
  !!    do iat=1,atoms%astruct%nat
  !!        write(*,*) 'iat, partial total_charge', iat, charge_per_atom(iat)
  !!    end do
  !!    write(*,*) 'total total_charge',total_charge
  !!    if (iproc==0) call write_partial_charges()
  !!end if

  do iorb=1,tmb%orbs%norb
      iat=tmb%orbs%onwhichatom(iorb)
      charge_per_atom(iat) = charge_per_atom(iat) + weight_matrix(iorb,iorb)
  end do
  if (iproc==0) then
      call write_partial_charges()
      call yaml_sequence_open('Multipole analysis (based on the Loewdin charges)')
      call calculate_dipole()
      call calculate_quadropole()
      call yaml_sequence_close()
  end if
  !!call support_function_multipoles()

  call deallocate_matrices(inv_ovrlp(1))

  call f_free(charge_per_atom)
  call f_free(weight_matrix)
  call f_release_routine()



  contains

    subroutine write_partial_charges
      use yaml_output
      character(len=20) :: atomname
      real(kind=8),dimension(2) :: charges
      call yaml_sequence_open('Loewdin charge analysis (charge / net charge)')
      total_charge=0.d0
      total_net_charge=0.d0
      do iat=1,atoms%astruct%nat
          call yaml_sequence(advance='no')
          call yaml_mapping_open(flow=.true.)
          atomname=atoms%astruct%atomnames(atoms%astruct%iatype(iat))
          charges(1)=-charge_per_atom(iat)
          charges(2)=-(charge_per_atom(iat)-real(atoms%nelpsp(atoms%astruct%iatype(iat)),kind=8))
          total_charge = total_charge + charges(1)
          total_net_charge = total_net_charge + charges(2)
          call yaml_map(trim(atomname),charges,fmt='(1es20.12)')
          call yaml_mapping_close(advance='no')
          call yaml_comment(trim(yaml_toa(iat,fmt='(i4.4)')))
      end do
      call yaml_sequence(advance='no')
      call yaml_map('total charge',total_charge,fmt='(es16.8)')
      call yaml_sequence(advance='no')
      call yaml_map('total net charge',total_net_charge,fmt='(es16.8)')
      call yaml_sequence_close()
    end subroutine write_partial_charges


    subroutine calculate_dipole()
      use yaml_output
      real(kind=8),dimension(3) :: dipole_elec, dipole_cores, dipole_net

      dipole_cores(1:3)=0._gp
      do iat=1,atoms%astruct%nat
         dipole_cores(1:3)=dipole_cores(1:3)+atoms%nelpsp(atoms%astruct%iatype(iat))*atoms%astruct%rxyz(1:3,iat)
      end do

      dipole_elec=0.d0
      do iat=1,atoms%astruct%nat
          dipole_elec(1:3) = dipole_elec(1:3) -charge_per_atom(iat)*atoms%astruct%rxyz(1:3,iat)
      end do

      dipole_net=dipole_cores+dipole_elec

      if (iproc==0) then
          !!call yaml_map('core dipole', dipole_cores)
          !!call yaml_map('electronic dipole', dipole_elec)
          call yaml_map('net dipole', dipole_net,fmt='(es12.5)')
      end if

    end subroutine calculate_dipole


    subroutine calculate_quadropole()
      use yaml_output
      real(kind=8),dimension(3,3) :: quadropole_elec, quadropole_cores, quadropole_net
      real(kind=8),dimension(3) :: charge_center_cores, charge_center_charge
      integer :: i, j
      real(kind=8) :: delta_term, rj, ri, q, qtot


      ! charge center of the cores
      charge_center_cores(1:3)=0.d0
      qtot=0.d0
      do iat=1,atoms%astruct%nat
          q=atoms%nelpsp(atoms%astruct%iatype(iat))
          charge_center_cores(1:3) = charge_center_cores(1:3) + q*atoms%astruct%rxyz(1:3,iat)
          qtot=qtot+q
      end do
      charge_center_cores=charge_center_cores/qtot


      ! charge center of the charge
      charge_center_charge(1:3)=0.d0
      qtot=0.d0
      do iat=1,atoms%astruct%nat
          q=-charge_per_atom(iat)
          charge_center_charge(1:3) = charge_center_charge(1:3) + q*atoms%astruct%rxyz(1:3,iat)
          qtot=qtot+q
      end do
      charge_center_charge=charge_center_charge/qtot


      quadropole_cores(1:3,1:3)=0._gp
      do iat=1,atoms%astruct%nat
         q=atoms%nelpsp(atoms%astruct%iatype(iat))
         do i=1,3
             do j=1,3
                 if (i==j) then
                     delta_term = atoms%astruct%rxyz(1,iat)**2 + atoms%astruct%rxyz(2,iat)**2 + atoms%astruct%rxyz(3,iat)**2
                 else
                     delta_term=0.d0
                 end if
                 rj=atoms%astruct%rxyz(j,iat)
                 ri=atoms%astruct%rxyz(i,iat)
                 quadropole_cores(j,i) = quadropole_cores(j,i) + q*(3.d0*rj*ri-delta_term)
                 !!quadropole_cores(j,i) = quadropole_cores(j,i) + &
                 !!                        atoms%nelpsp(atoms%astruct%iatype(iat))* &
                 !!                          (3.d0*atoms%astruct%rxyz(j,iat)*atoms%astruct%rxyz(i,iat)-delta_term)
             end do
         end do
      end do


      quadropole_elec(1:3,1:3)=0._gp
      do iat=1,atoms%astruct%nat
         q=-charge_per_atom(iat)
         do i=1,3
             do j=1,3
                 if (i==j) then
                     delta_term = (atoms%astruct%rxyz(1,iat)+(charge_center_cores(1)-charge_center_charge(1)))**2 + &
                                  (atoms%astruct%rxyz(2,iat)+(charge_center_cores(2)-charge_center_charge(2)))**2 + &
                                  (atoms%astruct%rxyz(3,iat)+(charge_center_cores(3)-charge_center_charge(3)))**2
                 else
                     delta_term=0.d0
                 end if
                 rj=atoms%astruct%rxyz(j,iat)+(charge_center_cores(j)-charge_center_charge(j))
                 ri=atoms%astruct%rxyz(i,iat)+(charge_center_cores(i)-charge_center_charge(i))
                 quadropole_elec(j,i) = quadropole_elec(j,i) + q*(3.d0*rj*ri-delta_term)
                 !!quadropole_elec(j,i) = quadropole_elec(j,i) + &
                 !!                       -charge_per_atom(iat)* &
                 !!                         (3.d0*atoms%astruct%rxyz(j,iat)*atoms%astruct%rxyz(i,iat)-delta_term)
             end do
         end do
      end do

      quadropole_net=quadropole_cores+quadropole_elec

      if (iproc==0) then
          !!call yaml_sequence_open('core quadropole')
          !!do i=1,3
          !!   call yaml_sequence(trim(yaml_toa(quadropole_cores(i,1:3),fmt='(es12.5)')))
          !!end do
          !!call yaml_sequence_close()

          !!call yaml_sequence_open('electronic quadropole')
          !!do i=1,3
          !!   call yaml_sequence(trim(yaml_toa(quadropole_elec(i,1:3),fmt='(es12.5)')))
          !!end do
          !!call yaml_sequence_close()

          call yaml_sequence_open('net quadropole')
          do i=1,3
             call yaml_sequence(trim(yaml_toa(quadropole_net(i,1:3),fmt='(es12.5)')))
          end do
          call yaml_sequence(advance='no')
          call yaml_map('trace of quadropole matrix',&
               quadropole_net(1,1)+quadropole_net(2,2)+quadropole_net(3,3),fmt='(es12.2)')
          call yaml_sequence_close()
      end if

    end subroutine calculate_quadropole

end subroutine loewdin_charge_analysis


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


subroutine calculate_multipoles(n1i, n2i, n3i, hgrids, phir, charge_center_elec, rxyz_center, &
           dipole_net, quadropole_net)
  use yaml_output
  implicit none
  ! Calling arguments
  integer,intent(in) :: n1i, n2i, n3i
  real(kind=8),dimension(3),intent(in) :: hgrids
  real(kind=8),dimension(n1i*n2i*n3i),intent(in) :: phir
  real(kind=8),dimension(3),intent(in) :: charge_center_elec, rxyz_center
  real(kind=8),dimension(3),intent(out) :: dipole_net
  real(kind=8),dimension(3,3),intent(out) :: quadropole_net

  integer :: i1, i2, i3, jj, iz, iy, ix, ii, i, j
  real(kind=8) :: q, x, y, z, qtot, ri, rj, delta_term
  real(kind=8),dimension(3) :: dipole_center, dipole_el
  real(kind=8),dimension(3,3) :: quadropole_center, quadropole_el


  !!call yaml_map('rxyz_center',rxyz_center,fmt='(es16.6)')
  !!call yaml_map('charge_center_elec',charge_center_elec,fmt='(es16.6)')
  !!call yaml_map('sum phir',sum(phir),fmt='(es16.6)')

  ! Dipole and quadropole of the support function
  dipole_el=0.d0
  quadropole_el=0.d0
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

              q = phir(jj)**2 * product(hgrids)
              x = ix*hgrids(1) + (rxyz_center(1)-charge_center_elec(1))
              y = iy*hgrids(2) + (rxyz_center(2)-charge_center_elec(2))
              z = iz*hgrids(3) + (rxyz_center(3)-charge_center_elec(3))

              ! Dipole part
              dipole_el(1) = dipole_el(1) + q*x
              dipole_el(2) = dipole_el(2) + q*y
              dipole_el(3) = dipole_el(3) + q*z
              qtot=qtot+q

              ! Quadrupole part
              do i=1,3
                  ri=get_r(i, x, y, z)
                  do j=1,3
                      rj=get_r(j, x, y, z)
                      if (i==j) then
                          delta_term = x**2 + y**2 + z**2
                      else
                          delta_term=0.d0
                      end if
                      quadropole_el(j,i) = quadropole_el(j,i) + q*(3.d0*rj*ri-delta_term)
                  end do
              end do
          end do
      end do
  end do

  ! Dipole of the center
  dipole_center(1) = -qtot*rxyz_center(1)
  dipole_center(2) = -qtot*rxyz_center(2)
  dipole_center(3) = -qtot*rxyz_center(3)

  ! Quadropole of the center
  quadropole_center=0.d0
  do i=1,3
      ri=rxyz_center(i)
      do j=1,3
          rj=rxyz_center(j)
          if (i==j) then
              delta_term = rxyz_center(1)**2 + rxyz_center(2)**2 + rxyz_center(3)**2
          else
              delta_term=0.d0
          end if
          quadropole_center(j,i) = quadropole_center(j,i) -qtot*(3.d0*rj*ri-delta_term)
      end do
  end do

  ! Net dipole and quadropole
  dipole_net = dipole_el + dipole_center
  quadropole_net = quadropole_el + quadropole_center


!  call yaml_sequence_open(trim(yaml_toa(it))//'('//trim(atomname)//')')
!  !call yaml_map('qtot',qtot)
!  call yaml_sequence(advance='no')
!  !call yaml_map('center dipole',dipole_center,fmt='(es16.6)')
!  !call yaml_map('electronic dipole',dipole_el,fmt='(es18.10)')
!  call yaml_map('net dipole',dipole_net,fmt='(es18.10)')
!  call yaml_sequence(advance='no')
!  !call yaml_sequence_open('center quadropole')
!  !do i=1,3
!  !   call yaml_sequence(trim(yaml_toa(quadropole_center(i,1:3),fmt='(es15.8)')))
!  !end do
!  !call yaml_sequence_close()
!  !call yaml_sequence_open('electronic quadropole')
!  !do i=1,3
!  !   call yaml_sequence(trim(yaml_toa(quadropole_el(i,1:3),fmt='(es15.8)')))
!  !end do
!  !call yaml_sequence_close()
!  call yaml_sequence_open('net quadropole')
!  do i=1,3
!     call yaml_sequence(trim(yaml_toa(quadropole_net(i,1:3),fmt='(es15.8)')))
!  end do
!  call yaml_sequence_close()
!  call yaml_sequence_close()

  contains

    function get_r(i, x, y, z)
      integer,intent(in) :: i
      real(kind=8),intent(in) :: x, y, z
      real(kind=8) :: get_r

      select case (i)
      case (1)
          get_r=x
      case (2)
          get_r=y
      case (3)
          get_r=z
      case default
          stop 'wrong value of i'
      end select
    end function get_r

end subroutine calculate_multipoles



subroutine support_function_multipoles(iproc, tmb, atoms, denspot)
  use module_base
  use module_types
  use yaml_output
  
  ! Calling arguments
  integer,intent(in) :: iproc
  type(DFT_wavefunction),intent(in) :: tmb
  type(atoms_data),intent(in) :: atoms
  type(DFT_local_fields), intent(inout) :: denspot

  integer :: ist, istr, iorb, iiorb, ilr, i
  real(kind=8),dimension(3) :: charge_center_elec
  real(kind=8),dimension(:),allocatable :: phir
  type(workarr_sumrho) :: w
  character(len=20) :: atomname
  real(kind=8),dimension(:,:),allocatable :: dipole_net
  real(kind=8),dimension(:,:,:),allocatable :: quadropole_net

  call f_routine(id='support_function_multipoles')

  phir = f_malloc(tmb%collcom_sr%ndimpsi_c,id='phir')
  dipole_net = f_malloc((/3,tmb%orbs%norb/),id='dipole_net')
  quadropole_net = f_malloc((/3,3,tmb%orbs%norb/),id='quadropole_net')

  call to_zero(3*tmb%orbs%norb, dipole_net(1,1))
  call to_zero(9*tmb%orbs%norb, quadropole_net(1,1,1))

  ist=1
  istr=1
  do iorb=1,tmb%orbs%norbp
      iiorb=tmb%orbs%isorb+iorb
      ilr=tmb%orbs%inwhichlocreg(iiorb)
      iat=tmb%orbs%onwhichatom(iiorb)
      call initialize_work_arrays_sumrho(1,tmb%lzd%Llr(ilr),.true.,w)
      ! Transform the support function to real space
      call daub_to_isf(tmb%lzd%llr(ilr), w, tmb%psi(ist), phir(istr))
      call deallocate_work_arrays_sumrho(w)
      ! Calculate the charge center
      call charge_center(tmb%lzd%llr(ilr)%d%n1i, tmb%lzd%llr(ilr)%d%n2i, tmb%lzd%llr(ilr)%d%n3i, &
           denspot%dpbox%hgrids, phir(istr), charge_center_elec)
      !write(*,*) 'ilr, tmb%lzd%llr(ilr)%locregcenter', iat, tmb%lzd%llr(ilr)%locregcenter
      atomname=trim(atoms%astruct%atomnames(atoms%astruct%iatype(iat)))
      call calculate_multipoles(tmb%lzd%llr(ilr)%d%n1i, tmb%lzd%llr(ilr)%d%n2i, tmb%lzd%llr(ilr)%d%n3i, &
           denspot%dpbox%hgrids, phir(istr), charge_center_elec, tmb%lzd%llr(ilr)%locregcenter, &
           dipole_net(:,iiorb), quadropole_net(:,:,iiorb))
      !write(*,*) 'charge_center', charge_center_elec
      ist = ist + tmb%lzd%Llr(ilr)%wfd%nvctr_c + 7*tmb%lzd%Llr(ilr)%wfd%nvctr_f
      istr = istr + tmb%lzd%Llr(ilr)%d%n1i*tmb%lzd%Llr(ilr)%d%n2i*tmb%lzd%Llr(ilr)%d%n3i
  end do
  if(istr/=tmb%collcom_sr%ndimpsi_c+1) then
      write(*,'(a,i0,a)') 'ERROR on process ',iproc,' : istr/=tmb%collcom_sr%ndimpsi_c+1'
      stop
  end if


  if (bigdft_mpi%nproc>1) then
      call mpiallred(dipole_net(1,1), 3*tmb%orbs%norb, mpi_sum, bigdft_mpi%mpi_comm)
      call mpiallred(quadropole_net(1,1,1), 9*tmb%orbs%norb, mpi_sum, bigdft_mpi%mpi_comm)
  end if

  if (iproc==0) then
      call yaml_sequence_open('Support functions moments')
      do iorb=1,tmb%orbs%norb
          iat=tmb%orbs%onwhichatom(iorb)
          atomname=trim(atoms%astruct%atomnames(atoms%astruct%iatype(iat)))
          call yaml_sequence_open('number'//trim(yaml_toa(iorb))// &
               ' (atom number ='//trim(yaml_toa(iat))//', type = '//trim(atomname)//')')
          call yaml_sequence(advance='no')
          call yaml_map('net dipole',dipole_net(:,iorb),fmt='(es18.10)')
          call yaml_sequence(advance='no')
          call yaml_sequence_open('net quadropole')
          do i=1,3
             call yaml_sequence(trim(yaml_toa(quadropole_net(i,1:3,iorb),fmt='(es15.8)')))
          end do
          call yaml_sequence_close()
          call yaml_sequence_close()
      end do
      call yaml_sequence_close()
  end if

  call f_free(phir)
  call f_free(dipole_net)
  call f_free(quadropole_net)
  call f_release_routine()

  

end subroutine support_function_multipoles
