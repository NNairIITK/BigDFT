!> @file
!!  Routines to do atomic analysis configuration
!! @author
!!    Copyright (C) 2009-2011 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 
!! LAST CHANGE
!!    14 Sep 2009


!>    Analyse atomic configurations
program find_angles
 implicit none
 integer, parameter :: ntypes=4,nnmax=15,nseg=1800,nsegr=10000
 real(kind=8), parameter :: rstep=0.005d0!factor=12.82,rad=3.2d0/factor
 real(kind=8) :: factor,rad,rij,distance
 integer :: i,j,istep,k,j1,nat,istart,iat,nrep,nata,natc
 real(kind=8) theta,anglemin,anglemax,th,facnorm
 integer, dimension(nnmax) :: nrstn
 integer, dimension(nseg,nnmax+1) :: isto
 integer, dimension(nsegr,0:nnmax) :: istor
 real(kind=8), dimension(nseg,nnmax+1) :: nisto
 real(kind=8), dimension(nsegr) :: nistor
 real(kind=8), dimension(nnmax+1) :: integrals
 character(len=1) :: whichone
 integer :: icount,ncount,ncountmax,tot,atcenter,atangles,nstep,jr,posout,iunit
 logical :: exists
 character(len=40) :: contcar
 character(len=40) :: xdatcar
 integer, dimension(:), pointer :: iatype
 real(kind=8), dimension(:,:), pointer :: pos
 ! Debug variables
 integer :: timecount, timeread_start, timeread_stop

 call system_clock(count_rate = timecount)

 inquire(file='input',exist=exists)
 if (exists) then
    print *,'Reading inputs from input file:'
    open(1,file='input',status='unknown')
    read(1,*)contcar
    read(1,*)nstep
    read(1,*)atcenter
    read(1,*)atangles
    read(1,*)rad
    read(1,*)xdatcar
    read(1,*)whichone
    read(1,*)nrep
    close(1)
    print *,'number of simulation steps:',nstep
    print *,'Atom type to be put in the center:',atcenter
    print *,'Atom type for calculating the angles:',atangles
    print *,'radius:',rad
    print *,'Xdatcar file: ',xdatcar
    print *,'Number of replica in each direction:',nrep
    if (whichone =='L') then
       print *,'LAMMPS type file'
    else if (whichone =='V') then
       print *,'VASP type file'
    else if (whichone =='B') then
       print *,'BigDFT type files'
    else
       print *,'ERROR: type of trajectory not recognized'
       stop
    end if
    if (whichone =='B') then
      !initial file to read the box features
      read(contcar,'(i5)') posout
      !write(fn4,'(i5.5)') posout
      !contcar='posmd_'//fn4
       contcar='posinp'
     end if

 else
    print *,"Choose trajectory type ('V' for VASP, 'L' for LAMMPS, 'B' for BigDFT): "
    read(5,*)whichone
    if (whichone =='V') then
       print *,'Contcar file: '
       read(5,*)contcar
    else if (whichone =='L') then
       print *,'Dump file: '
       read(5,*)contcar
    else if (whichone =='B') then
       print *,'posout starting file(value): '
       read(5,*)posout
    else
       print *,'ERROR: type of trajectory not recognized'
       stop
    end if
    print *,'Enter the number of simulation steps:'
    read(5,*)nstep
    print *,'Enter the atom type to be put in the center:'
    read(5,*)atcenter
    print *,'Enter the atom type for calculating the angles:'
    read(5,*)atangles
    print *,'Enter the radius:'
    read(5,*)rad
    print *,'Enter the number of replica in each direction:'
    read(5,*)nrep

    if (whichone =='V') then
       print *,'Xdatcar file: '
       read(5,*)xdatcar
    else if (whichone =='L') then
       xdatcar=contcar
    else if (whichone =='B') then
       !initial file to read the box features
       !write(fn4,'(i5.5)') posout
       !contcar='posmd_'//fn4
       contcar='posinp'
    end if
 end if
 call box_features(whichone,contcar,nrep,nat,ntypes,iatype,pos,factor)

 print *,'The factor (box size - ONLY CUBIC) is=',factor*real(nrep,kind=8)

 rad=rad/(factor*real(nrep,kind=8))

 istep=0
 isto=0
 istor=0
 anglemin=180.d0
 anglemax=0.d0
 ncountmax=0

 !for non-BigDFT cases the output is on one file only
 if (whichone /= 'B') then
    iunit=12
    open(iunit,file=xdatcar,status='unknown')
 else
    iunit=posout
 end if

 if (whichone == 'V') then
    !intialise the reading, unsed only for vasp
    do i=1,6
       read(iunit,*)
    end do
 end if

 call system_clock(timeread_start)
 !this loop is over the frames
 loop_step: do j1=1,nstep
    istep=istep+1
    if (mod(istep,100)==0) write(*,'(1x,i5)',advance='no')istep
    !now read atomic positions
    call read_pos(iunit,whichone,nat,pos,nrep)
    !increment the posout for BigDFT results
    if (whichone == 'B') iunit=iunit+1
    if (whichone == 'V' .and. (j1 /= nstep)) read(12,*)

    !process the atom species
    loop_center: do i=1,nat*nrep**3
       if (iatype(i)==atcenter) then
          ncount=0
          do j=1,nrep**3*nat
             
             rij=distance(pos(1,i),pos(1,j))
             if (iatype(j)==atangles .and. rij <= rad .and. i /= j) then
                ncount=ncount+1
                !positions of nearest neighbors
                if(ncount > nnmax) then
                   ncount=ncount-1
                else
                   nrstn(ncount)=j
                end if
             end if
             if (iatype(j)==atangles .and. i /= j) then
                !calculate the index of the istogram for the partial g(r)
                jr=int(real(nrep,kind=8)*rij/rstep)
                if (jr > nsegr) then
                   print *,'ERROR: jr>nsegr',jr,nsegr,'choose more points or increase rstep'
                   stop
                end if
                istor(jr,0)=istor(jr,0)+1
             end if
          end do
          ncountmax=max(ncount,ncountmax)
          !now evaluate the angles
          !and form the istogram with the angular distribution
          do icount=1,ncount
             rij=distance(pos(1,i),pos(1,nrstn(icount)))
             jr=int(real(nrep,kind=8)*rij/rstep)
             istor(jr,ncount)=istor(jr,ncount)+1
             do k=icount+1,ncount
                th=theta(pos(1,nrstn(icount)),pos(1,i),pos(1,nrstn(k)))
                anglemin=min(anglemin,th)
                anglemax=max(anglemax,th)
                j=int(th*real(nseg,kind=8)/180.d0) + 1
                isto(j,ncount)=isto(j,ncount)+1
                isto(j,nnmax+1)=isto(j,nnmax+1)+1
             end do
             !print *,'AAAAAAAAAAA',isto(j,ncount),ncount
          end do
       end if
    end do loop_center
    !print *,'ciao'
    !stop
 end do loop_step
 !close the reading file in the non-BigDFT case
 if (whichone /= 'B') then
    close(iunit)
 end if
 print *,''
 print *,'anglemin,anglemax,ncountmax=',anglemin,anglemax,ncountmax
 call system_clock(timeread_stop)
 write(0, "(A,F20.8,A)") "Global file read: ", &
      & real(timeread_stop - timeread_start) / real(timecount) , "s"
 call finaliseExtract()

 !values of the normalisations
 integrals(:)=0.d0

!Normalise histogram

 tot=0
 do i=1,nseg
    tot=tot+isto(i,nnmax+1)
 end do
 do j=1,nnmax+1
    if (tot/=0) then
       do i=1,nseg
          nisto(i,j)=real(isto(i,j),kind=8)/real(tot,kind=8)
          integrals(j)=integrals(j)+nisto(i,j)
       end do
    else
       do i=1,nseg
          nisto(i,j)=0.d0
       end do
    end if
 end do

 !print istogram
 do i=1,nseg
    write(13,'(f10.4,20(1pe12.4))')180.d0*real(i,kind=8)/real(nseg,kind=8),(smearing(nseg,nisto(:,j),0.0d0,i),j=1,nnmax+1)!(nisto(i,j),j=1,nnmax+1)
 end do

 print *,'Normalisations:'
 do j=1,nnmax+1
    print *,'# ',j,'int=',integrals(j)
 end do
!Normalise g(r)
 natc=0
 nata=0
 do i=1,nat
    if (iatype(i)==atcenter) then
       natc=natc+1
    else if (iatype(i)==atangles) then
       nata=nata+1
    end if
 end do
 if (atcenter == atangles) nata=natc
 facnorm=16.d0*atan(1.d0)*rstep*real(natc*nata*nrep**3*nstep,kind=8)/factor**2
 
!!$ tot=0
!!$ do i=1,nsegr
!!$    tot=tot+istor(i,0)*facnorm/(rstep*factor*real(i,kind=8))**2
!!$ end do
 do i=1,nsegr
!    nistor(i)=real(istor(i,0)*nrep**3,kind=8)*8.d0*factor**3*atan(1.d0)/(real(tot,kind=8)*1.d0)
    nistor(i)=real(istor(i,0),kind=8)/facnorm/(rstep*factor*real(i,kind=8))**2
 end do
 print *,'tot =', tot
 !print g(r)
 do i=1,nsegr
    write(14,*)rstep*i*factor,nistor(i)
 end do

 deallocate(iatype,pos)

contains

  !> calculate the smearing of the istogram
  pure function smearing(nseg,nisto,sigma,i)
    implicit none
    integer, intent(in) :: nseg
    real(kind=8), dimension(nseg), intent(in) :: nisto
    real(kind=8), intent(in) :: sigma !<gaussian spread in units of tenths of a degree
    integer, intent(in) :: i !< output point
    real(kind=8) :: smearing
    !local variables
    integer :: j
    real(kind=8) :: exponent

    if (sigma==0.0d0) then
       smearing=nisto(i)
       return
    end if
    
    smearing=0.0d0
    do j=1,nseg
       exponent=real(i-nisto(j),kind=8)
       exponent=exponent/sigma
       exponent=0.5d0*exponent**2
       smearing=smearing+exp(-exponent)
    end do

  end function smearing

  subroutine box_features(whichone,contcar,nrep,nat,ntypes,iatype,pos,factor)
    use BigDFT_API
    use module_interfaces
    use m_ab6_symmetry
    implicit none
    character(len=1), intent(in) :: whichone
    character(len=40), intent(in) :: contcar
    integer, intent(in) :: ntypes,nrep
    integer, intent(out) :: nat
    real(kind=8), intent(out) :: factor
    integer, dimension(:), pointer :: iatype
    real(kind=8), dimension(:,:), pointer :: pos
    !local variables
    integer :: i,ityp
    real(kind=8) :: xlo,xhi
    integer, dimension(ntypes) :: natoms
    type(atoms_data) :: atoms
    character(len=*), parameter :: subname='box_features'

    !!allocate arrays and set atom types
    if (whichone /= 'B') open(11,file=contcar,status='unknown')
    if (whichone == 'V') then
       print *,'Contcar file: ',contcar
       read(11,*)
       read(11,*)factor
       read(11,*)
       read(11,*)
       read(11,*)
       read(11,*)natoms
       close(11)
       nat=0
       do i=1,ntypes
          nat=nat+natoms(i)
       end do
       print *,'You have ',nat,' atoms, divided in groups of', natoms
       allocate(iatype(nrep**3*nat),pos(3,nrep**3*nat))
       !set atom types
       iatype(1:natoms(1))=1
       istart=natoms(1)
       do i=2,ntypes
          do j=istart+1,istart+natoms(i)
             iatype(j)=i
          end do
          istart=istart+natoms(i)
       end do

    else if (whichone == 'L') then
       print *,'Dump file: ',contcar
       read(11,*)
       read(11,*)
       read(11,*)
       read(11,*)nat
       read(11,*)
       read(11,*)xlo,xhi
       read(11,*)
       read(11,*)
       read(11,*)
       !read atom
       allocate(iatype(nrep**3*nat),pos(3,nrep**3*nat))
       do i=1,nat
          read(11,*)  iat,ityp
          if (ityp > ntypes) then
             print *,'ERROR: change ntypes'
             stop
          end if
          iatype(iat)=ityp
       enddo
       print *,'You have ',nat,' atoms'
       close(11)
       factor=xhi-xlo
    else if (whichone == 'B') then
       !open the first file to check box features
!print *,'here'
       call read_atomic_file(trim(contcar),0,atoms%astruct)
       call allocate_atoms_nat(atoms, subname)
       call allocate_atoms_ntypes(atoms, subname)
       nat=atoms%astruct%nat
!print *,'nat',nat
       allocate(iatype(nrep**3*nat),pos(3,nrep**3*nat))
       do i=1,nat
          iatype(i)=atoms%astruct%iatype(i)
       enddo
       factor=atoms%astruct%cell_dim(1) * Bohr_Ang

       deallocate(atoms%astruct%rxyz)
       call deallocate_atomic_structure(atoms%astruct,"box_features") 
!       call deallocate_atoms(atoms, "box_features")

    end if

    !replica of the atom types
    do i=2,nrep**3
       do iat=1,nat
          iatype((i-1)*nat+iat)=iatype(iat)
       end do
    end do

  END SUBROUTINE box_features


end program find_angles


subroutine read_pos(iunit,whichone,nat,pos,nrep)
  use BigDFT_API
  use module_interfaces
  use m_ab6_symmetry
  implicit none
  character(len=1), intent(in) :: whichone
  integer, intent(in) :: iunit,nat,nrep
  real(kind=8), dimension(3,nrep**3*nat), intent(out) :: pos
  !local variables
  character(len=5) :: fn4
  integer :: i,iat,ityp,i1,i2,i3
  real(kind=8) :: x,y,z,vx,vy,vz,xlo,xhi,ylo,yhi,zlo,zhi,alat(3)
  type(atoms_data) :: atoms
  character(len=*), parameter :: subname='read_pos'


  if (whichone == 'V') then
     do i=1,nat
        read(iunit,*)pos(1,i),pos(2,i),pos(3,i)
     end do
  else if (whichone == 'L') then
     do i=1,5
        read(iunit,*)
     end do
     read(iunit,*)xlo,xhi
     read(iunit,*)ylo,yhi
     read(iunit,*)zlo,zhi
     read(iunit,*)
     do i=1,nat
        read(iunit,*)  iat,ityp,x,y,z,vx,vy,vz
        pos(1,iat)=(x-xlo)/(xhi-xlo)
        pos(2,iat)=(y-ylo)/(yhi-ylo)
        pos(3,iat)=(z-zlo)/(zhi-zlo)
     enddo
  else if (whichone == 'B') then
     !use the BigDFT call with iunit to control the posout
     write(fn4,'(i5.5)') iunit
     call read_atomic_file('posmd_'//fn4,0,atoms%astruct)
     call allocate_atoms_nat(atoms, subname)
     call allocate_atoms_ntypes(atoms, subname)
     !transform the positions in reduced coordinates
     alat(1) = atoms%astruct%cell_dim(1)
     if (atoms%astruct%geocode == 'F') alat(1) = 1.d0
     alat(2) = atoms%astruct%cell_dim(2)
     if (atoms%astruct%geocode == 'F' .or. atoms%astruct%geocode == 'S') alat(2) = 1.d0
     alat(3) = atoms%astruct%cell_dim(3)
     if (atoms%astruct%geocode == 'F') alat(3) = 1.d0
     do iat=1,nat
        pos(1,iat)=atoms%astruct%rxyz(1,iat)/alat(1)
        pos(2,iat)=atoms%astruct%rxyz(2,iat)/alat(2)
        pos(3,iat)=atoms%astruct%rxyz(3,iat)/alat(3)
     enddo
     call deallocate_atomic_structure(atoms%astruct, 'distance')
  end if

  !replica of the atom positions
  do i1=1,nrep
     do i2=1,nrep
        do i3=1,nrep
           if (i1 + i2 + i3 /= 3) then
              i=i1+(i2-1)*nrep+(i3-1)*nrep**2
              do iat=1,nat
                 pos(1,(i-1)*nat+iat)=(pos(1,iat)+real(i1-1,kind=8))/real(nrep,kind=8)
                 pos(2,(i-1)*nat+iat)=(pos(2,iat)+real(i2-1,kind=8))/real(nrep,kind=8)
                 pos(3,(i-1)*nat+iat)=(pos(3,iat)+real(i3-1,kind=8))/real(nrep,kind=8)
              end do
           end if
        end do
     end do
  end do

  if (nrep > 1) then
     do iat=1,nat
        pos(1,iat)=pos(1,iat)/real(nrep,kind=8)
        pos(2,iat)=pos(2,iat)/real(nrep,kind=8)
        pos(3,iat)=pos(3,iat)/real(nrep,kind=8)
     end do
  end if

END SUBROUTINE read_pos


function theta(A,O,B)
 implicit none
 real(kind=8), dimension(3), intent(in) :: A,O,B
 real(kind=8) :: theta
 !local variables
 integer :: i
 real(kind=8), dimension(3) :: u,v
 real(kind=8) :: modu,modv,udv


 do i=1,3
    u(i)=A(i)-O(i)
    v(i)=B(i)-O(i)
 end do

 udv=u(1)*v(1)+u(2)*v(2)+u(3)*v(3)
 modu=sqrt(u(1)**2+u(2)**2+u(3)**2)
 modv=sqrt(v(1)**2+v(2)**2+v(3)**2)
 theta=45.d0/datan(1.d0)*dacos(udv/(modu*modv))
end function theta

function distance(A,B)
 implicit none
 real(kind=8), dimension(3), intent(in) :: A
 real(kind=8), dimension(3), intent(inout) :: B
 real(kind=8) :: distance
 !local variables
 integer :: i
 real(kind=8) :: mindist,shift,di
 !real(kind=8), dimension(3) :: v

 !valid for reduced coordinates
 !shift the B vector on its periodic image
 mindist=0.d0
 do i=1,3
    di=A(i)-B(i)
    !periodic image, if distance is bigger than half of the box
    shift=real(floor(di+0.5d0),kind=8)
    mindist=mindist+(di-shift)**2
    B(i)=B(i)+shift
 end do

!!$ do i=1,3
!!$    mini=10.d0
!!$    do ei=-1,1
!!$       distance=(A(i)-real(ei,kind=8)-B(i))**2
!!$       if (distance < mini) then
!!$          mini=min(mini,distance)
!!$          v(i)=real(ei,kind=8)+B(i)
!!$       end if
!!$    end do
!!$    mindist=mindist+mini
!!$ end do
 distance=sqrt(mindist)
 !print *,A,B,v,distance
 !new vector
 !B=v

end function distance


