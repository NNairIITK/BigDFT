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

integer :: ix, iy, iz, ix0, iy0, iz0, iiAt, jj, iorb, i1, i2, i3, istart, ii, istat, iat
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

allocate(phir(tmb%lzd%glr%d%n1i*tmb%lzd%glr%d%n2i*tmb%lzd%glr%d%n3i), stat=istat)

call initialize_work_arrays_sumrho(tmb%lzd%glr,w)
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
deallocate(phir, stat=istat)


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
  integer :: iseg, jj, j0, j1, ii, i3, i2, i0, i1, i, ishift, iat, ldim, gdim, jjj, istat
  character(len=10) :: num
  character(len=20) :: filename
  real(kind=8), dimension(:), allocatable :: lphi, phi


    ldim=llr%wfd%nvctr_c+7*llr%wfd%nvctr_f
    gdim=glr%wfd%nvctr_c+7*glr%wfd%nvctr_f
    allocate(lphi(ldim), stat=istat)
    allocate(phi(gdim), stat=istat)
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



subroutine local_potential_dimensions(Lzd,orbs,ndimfirstproc)
  use module_base
  use module_types
  use module_xc
  implicit none
  integer, intent(in) :: ndimfirstproc
  type(local_zone_descriptors), intent(inout) :: Lzd
  type(orbitals_data), intent(inout) :: orbs
  !local variables
  character(len=*), parameter :: subname='local_potential_dimensions'
  logical :: newvalue
  integer :: i_all,i_stat,ii,iilr,ilr,iorb,iorb2,nilr,ispin
  integer, dimension(:,:), allocatable :: ilrtable
  
  if(Lzd%nlr > 1) then
     allocate(ilrtable(orbs%norbp,2),stat=i_stat)
     call memocc(i_stat,ilrtable,'ilrtable',subname)
     !call to_zero(orbs%norbp*2,ilrtable(1,1))
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
              if(orbs%spinsgn(orbs%isorb+iorb) <= 0.0_gp) then
                 orbs%ispot(iorb)=lzd%ndimpotisf + &
                      1 + lzd%llr(ilr)%d%n1i*lzd%llr(ilr)%d%n2i*lzd%llr(ilr)%d%n3i
              end if
           end if
        end do
        lzd%ndimpotisf = lzd%ndimpotisf + &
             lzd%llr(ilr)%d%n1i*lzd%llr(ilr)%d%n2i*lzd%llr(ilr)%d%n3i*orbs%nspin
     end do
     !part which refers to exact exchange (only meaningful for one region)
     if (xc_exctXfac() /= 0.0_gp) then
        lzd%ndimpotisf = lzd%ndimpotisf + &
             max(max(lzd%llr(ilr)%d%n1i*lzd%llr(ilr)%d%n2i*lzd%llr(ilr)%d%n3i*orbs%norbp,ndimfirstproc*orbs%norb),1)
     end if

  else 
     allocate(ilrtable(1,2),stat=i_stat)
     call memocc(i_stat,ilrtable,'ilrtable',subname)
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
     if (xc_exctXfac() /= 0.0_gp) then
        lzd%ndimpotisf = lzd%ndimpotisf + &
             max(max(lzd%Glr%d%n1i*lzd%Glr%d%n2i*lzd%Glr%d%n3i*orbs%norbp,ndimfirstproc*orbs%norb),1)
     end if


  end if


  i_all=-product(shape(ilrtable))*kind(ilrtable)
  deallocate(ilrtable,stat=i_stat)
  call memocc(i_stat,i_all,'ilrtable',subname)

end subroutine local_potential_dimensions



subroutine print_orbital_distribution(iproc, nproc, orbs)
use module_base
use module_types
use yaml_output
implicit none

integer, intent(in) :: iproc, nproc
type(orbitals_data), intent(in) :: orbs

! Local variables
integer :: jproc, len1, len2, space1, space2, jpst, norb0,  norb1
logical :: written

!!write(*,'(1x,a)') '------------------------------------------------------------------------------------'
!!written=.false.
!!write(*,'(1x,a)') '>>>> Partition of the basis functions among the processes.'
call yaml_open_map('Support function repartition')
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
call yaml_close_map()
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
           energs, nlpsp, input, &
           energy, energyDiff, energyold)
  use module_base
  use module_types
  use module_interfaces, except_this_one => build_ks_orbitals
  implicit none
  
  ! Calling arguments
  integer:: iproc, nproc
  type(DFT_wavefunction),intent(inout) :: tmb, KSwfn
  type(atoms_data), intent(inout) :: at
  real(gp), dimension(3,at%astruct%nat), intent(in) :: rxyz
  type(DFT_local_fields), intent(inout) :: denspot
  type(GPU_pointers), intent(inout) :: GPU
  type(energy_terms),intent(inout) :: energs
  type(DFT_PSP_projectors), intent(inout) :: nlpsp
  type(input_variables),intent(in) :: input
  real(kind=8),intent(out) :: energy, energyDiff, energyold

  ! Local variables
  type(orbitals_data) :: orbs
  type(communications_arrays) :: comms
  type(sparseMatrix) :: ham_small ! for FOE
  real(gp) :: fnrm
  integer :: infoCoeff, nvctrp, npsidim_global
  real(kind=8),dimension(:),pointer :: phi_global, phiwork_global
  character(len=*),parameter :: subname='build_ks_orbitals'

  !debug
  integer :: iorb, jorb, ist, jst, ierr, i, istat, iall
  real(kind=8) :: ddot, tt


  ! Get the expansion coefficients
  ! Start with a "clean" density, i.e. without legacy from previous mixing steps
  call communicate_basis_for_density_collective(iproc, nproc, tmb%lzd, &
       max(tmb%npsidim_orbs,tmb%npsidim_comp), tmb%orbs, tmb%psi, tmb%collcom_sr)
  call sumrho_for_TMBs(iproc, nproc, KSwfn%Lzd%hgrids(1), KSwfn%Lzd%hgrids(2), KSwfn%Lzd%hgrids(3), &
       tmb%collcom_sr, tmb%linmat%denskern, KSwfn%Lzd%Glr%d%n1i*KSwfn%Lzd%Glr%d%n2i*denspot%dpbox%n3d, denspot%rhov)
  call updatePotential(input%ixc,input%nspin,denspot,energs%eh,energs%exc,energs%evxc)

  tmb%can_use_transposed=.false.
  call get_coeff(iproc, nproc, LINEAR_MIXDENS_SIMPLE, KSwfn%orbs, at, rxyz, denspot, GPU, infoCoeff, &
       energs, nlpsp, input%SIC, tmb, fnrm, .true., .false., .true., ham_small, 0, 0, 0, 0, &
       input%lin%order_taylor,input%calculate_KS_residue)


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

  if(tmb%can_use_transposed) then
      iall=-product(shape(tmb%psit_c))*kind(tmb%psit_c)
      deallocate(tmb%psit_c, stat=istat)
      call memocc(istat, iall, 'tmb%psit_c', subname)
      iall=-product(shape(tmb%psit_f))*kind(tmb%psit_f)
      deallocate(tmb%psit_f, stat=istat)
      call memocc(istat, iall, 'tmb%psit_f', subname)
  end if

  ! Create communication arrays for support functions in the global box
  
  call nullify_orbitals_data(orbs)
  call copy_orbitals_data(tmb%orbs, orbs, subname)
  call orbitals_communicators(iproc, nproc, tmb%lzd%glr, orbs, comms)


  ! Transform the support functions to the global box
  ! WARNING: WILL NOT WORK WITH K-POINTS, CHECK THIS
  npsidim_global=max(tmb%orbs%norbp*(tmb%lzd%glr%wfd%nvctr_c+7*tmb%lzd%glr%wfd%nvctr_f), &
                     tmb%orbs%norb*comms%nvctr_par(iproc,0)*orbs%nspinor)
  allocate(phi_global(npsidim_global))
  allocate(phiwork_global(npsidim_global))
  call small_to_large_locreg(iproc, tmb%npsidim_orbs, &
       tmb%orbs%norbp*(tmb%lzd%glr%wfd%nvctr_c+7*tmb%lzd%glr%wfd%nvctr_f), tmb%lzd, &
       KSwfn%lzd, tmb%orbs, tmb%psi, phi_global, to_global=.true.)
  call transpose_v(iproc, nproc, orbs, tmb%lzd%glr%wfd, comms, phi_global, work=phiwork_global)


  ! WARNING: WILL NOT WORK WITH K-POINTS, CHECK THIS
  nvctrp=comms%nvctr_par(iproc,0)*orbs%nspinor
  call dgemm('n', 'n', nvctrp, KSwfn%orbs%norb, tmb%orbs%norb, 1.d0, phi_global, nvctrp, tmb%coeff(1,1), &
             tmb%orbs%norb, 0.d0, phiwork_global, nvctrp)
  
  call untranspose_v(iproc, nproc, KSwfn%orbs, tmb%lzd%glr%wfd, KSwfn%comms, phiwork_global, work=phi_global)  

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

   call deallocate_orbitals_data(orbs, subname)
   call deallocate_communications_arrays(comms, subname)

  ! To get consistent values of the energy and the Kohn-Sham residue with those
  ! which will be calculated by the cubic restart.
  call communicate_basis_for_density_collective(iproc, nproc, tmb%lzd, &
       max(tmb%npsidim_orbs,tmb%npsidim_comp), tmb%orbs, tmb%psi, tmb%collcom_sr)
  call sumrho_for_TMBs(iproc, nproc, KSwfn%Lzd%hgrids(1), KSwfn%Lzd%hgrids(2), KSwfn%Lzd%hgrids(3), &
       tmb%collcom_sr, tmb%linmat%denskern, KSwfn%Lzd%Glr%d%n1i*KSwfn%Lzd%Glr%d%n2i*denspot%dpbox%n3d, denspot%rhov)
  call updatePotential(input%ixc,input%nspin,denspot,energs%eh,energs%exc,energs%evxc)
  tmb%can_use_transposed=.false.
  call get_coeff(iproc, nproc, LINEAR_MIXDENS_SIMPLE, KSwfn%orbs, at, rxyz, denspot, GPU, infoCoeff, &
       energs, nlpsp, input%SIC, tmb, fnrm, .true., .false., .true., ham_small, 0, 0, 0, 0, &
       input%lin%order_taylor, input%calculate_KS_residue, updatekernel=.false.)
  energy=energs%ebs-energs%eh+energs%exc-energs%evxc-energs%eexctX+energs%eion+energs%edisp
  energyDiff=energy-energyold
  energyold=energy

  if(tmb%can_use_transposed) then
      iall=-product(shape(tmb%psit_c))*kind(tmb%psit_c)
      deallocate(tmb%psit_c, stat=istat)
      call memocc(istat, iall, 'tmb%psit_c', subname)
      iall=-product(shape(tmb%psit_f))*kind(tmb%psit_f)
      deallocate(tmb%psit_f, stat=istat)
      call memocc(istat, iall, 'tmb%psit_f', subname)
  end if

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
