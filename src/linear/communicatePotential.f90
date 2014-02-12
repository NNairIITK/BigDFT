!> @file
!! Routine to communicate the potential in linear version
!! @author
!!    Copyright (C) 2011-2012 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS


subroutine initialize_communication_potential(iproc, nproc, nscatterarr, orbs, lzd, comgp)
  use module_base
  use module_types
  implicit none
  
  ! Calling arguments
  integer,intent(in):: iproc, nproc
  integer,dimension(0:nproc-1,4),intent(in):: nscatterarr !n3d,n3p,i3s+i3xcsh-1,i3xcsh
  type(orbitals_data),intent(in):: orbs
  type(local_zone_descriptors),intent(in):: lzd
  type(p2pComms),intent(out):: comgp
  
  ! Local variables
  integer:: is1, ie1, is2, ie2, is3, ie3, ilr, ii, iorb, iiorb, jproc, kproc, istat, istsource
  integer:: ioverlap, is3j, ie3j, is3k, ie3k, mpidest, istdest, ioffsetx, ioffsety, ioffsetz
  integer :: is3min, ie3max, tag, ncount, ierr, nmaxoverlap
  logical :: datatype_defined
  character(len=*),parameter:: subname='initialize_communication_potential'

  call timing(iproc,'init_commPot  ','ON')
  
  call nullify_p2pComms(comgp)

  allocate(comgp%ise(6,0:nproc-1), stat=istat)
  call memocc(istat, comgp%ise, 'comgp%ise', subname)
  
  ! Determine the bounds of the potential that we need for
  ! the orbitals on this process.
  iiorb=0
  do jproc=0,nproc-1
      is1=1000000000
      ie1=-1000000000
      is2=1000000000
      ie2=-1000000000
      is3=1000000000
      ie3=-1000000000
      do iorb=1,orbs%norb_par(jproc,0)
          
          iiorb=iiorb+1 
          ilr=orbs%inwhichlocreg(iiorb)
      
          ii=1+lzd%Llr(ilr)%nsi1
          if(ii < is1 .or. iorb==1) then
              is1=ii
          end if
          ii=lzd%Llr(ilr)%nsi1+lzd%Llr(ilr)%d%n1i
          if(ii > ie1 .or. iorb==1) then
              ie1=ii
          end if
      
          ii=1+lzd%Llr(ilr)%nsi2
          if(ii < is2 .or. iorb==1) then
              is2=ii
          end if
          ii=lzd%Llr(ilr)%nsi2+lzd%Llr(ilr)%d%n2i
          if(ii > ie2 .or. iorb==1) then
              ie2=ii
          end if
      
          ii=1+lzd%Llr(ilr)%nsi3
          if(ii < is3 .or. iorb==1) then
              is3=ii
          end if
          ii=lzd%Llr(ilr)%nsi3+lzd%Llr(ilr)%d%n3i
          if(ii > ie3 .or. iorb==1) then
              ie3=ii
          end if
      
      end do
      comgp%ise(1,jproc)=is1
      comgp%ise(2,jproc)=ie1
      comgp%ise(3,jproc)=is2
      comgp%ise(4,jproc)=ie2
      comgp%ise(5,jproc)=is3
      comgp%ise(6,jproc)=ie3
  end do


  
  ! Determine how many slices each process receives.
  allocate(comgp%noverlaps(0:nproc-1), stat=istat)
  call memocc(istat, comgp%noverlaps, 'comgp%noverlaps', subname)
  nmaxoverlap=0
  do jproc=0,nproc-1
      is3j=comgp%ise(5,jproc)
      ie3j=comgp%ise(6,jproc)
      mpidest=jproc
      ioverlap=0
      do kproc=0,nproc-1
          is3k=nscatterarr(kproc,3)+1
          ie3k=is3k+nscatterarr(kproc,2)-1
          if(is3j<=ie3k .and. ie3j>=is3k) then
              ioverlap=ioverlap+1
              !if(iproc==0) write(*,'(2(a,i0),a)') 'process ',jproc,' gets potential from process ',kproc,'.' 
          !TAKE INTO ACCOUNT THE PERIODICITY HERE
          else if(ie3j > lzd%Glr%d%n3i .and. lzd%Glr%geocode /= 'F') then
              ie3j = comgp%ise(6,jproc) - lzd%Glr%d%n3i
              if(ie3j>=is3k) then
                 ioverlap=ioverlap+1
              end if
              if(is3j <= ie3k)then
                 ioverlap=ioverlap+1
              end if
          end if
      end do
      if (ioverlap>nmaxoverlap) nmaxoverlap=ioverlap
      comgp%noverlaps(jproc)=ioverlap
      !!if(iproc==0) write(*,'(2(a,i0),a)') 'Process ',jproc,' gets ',ioverlap,' potential slices.'
  end do
  
  ! Determine the parameters for the communications.
  allocate(comgp%comarr(6,nmaxoverlap,0:nproc-1))
  call memocc(istat, comgp%comarr, 'comgp%comarr', subname)
  call to_zero(6*nmaxoverlap*nproc, comgp%comarr(1,1,0))
  allocate(comgp%mpi_datatypes(0:nmaxoverlap,0:nproc-1), stat=istat)
  call memocc(istat, comgp%mpi_datatypes, 'comgp%mpi_datatypes', subname)
  call to_zero((nmaxoverlap+1)*nproc, comgp%mpi_datatypes(1,0))
  comgp%nrecvBuf = 0
  is3min=0
  ie3max=0

  ! Only do this if we have more than one MPI task
  nproc_if: if (nproc>1) then
      do jproc=0,nproc-1
          is3j=comgp%ise(5,jproc)
          ie3j=comgp%ise(6,jproc)
          mpidest=jproc
          ioverlap=0
          istdest=1
          datatype_defined =.false.
          do kproc=0,nproc-1
              is3k=nscatterarr(kproc,3)+1
              ie3k=is3k+nscatterarr(kproc,2)-1
      !SHOULD TAKE INTO ACCOUNT THE PERIODICITY HERE
      !Need to split the region
              if(is3j<=ie3k .and. ie3j>=is3k) then
                  is3=max(is3j,is3k) ! starting index in z dimension for data to be sent
                  ie3=min(ie3j,ie3k) ! ending index in z dimension for data to be sent
                  ioffsetz=is3-is3k ! starting index (in z direction) of data to be sent (actually it is the index -1)
                  ioffsety=comgp%ise(3,jproc)-1
                  ioffsetx=comgp%ise(1,jproc)
                  ioverlap=ioverlap+1
                  if(is3<is3min .or. ioverlap==1) then
                      is3min=is3
                  end if
                  if(ie3>ie3max .or. ioverlap==1) then
                      ie3max=ie3
                  end if
                  istsource = ioffsetz*lzd%glr%d%n1i*lzd%glr%d%n2i + ioffsety*lzd%glr%d%n1i + ioffsetx
                  ncount = 1
                  comgp%comarr(1,ioverlap,jproc)=kproc
                  comgp%comarr(2,ioverlap,jproc)=istsource
                  comgp%comarr(3,ioverlap,jproc)=jproc
                  comgp%comarr(4,ioverlap,jproc)=istdest
                  comgp%comarr(5,ioverlap,jproc)=ie3-is3+1
                  comgp%comarr(6,ioverlap,jproc)=lzd%glr%d%n1i*lzd%glr%d%n2i
                  if (.not. datatype_defined) then
                      call mpi_type_vector(comgp%ise(4,jproc)-comgp%ise(3,jproc)+1, comgp%ise(2,jproc)-comgp%ise(1,jproc)+1, &
                           lzd%glr%d%n1i, mpi_double_precision, comgp%mpi_datatypes(0,jproc), ierr)
                      call mpi_type_commit(comgp%mpi_datatypes(0,jproc), ierr)
                      !comgp%mpi_datatypes(2,jproc)=1
                      datatype_defined=.true.
                  end if
    
                  istdest = istdest + &
                            (ie3-is3+1)*(comgp%ise(2,jproc)-comgp%ise(1,jproc)+1)*(comgp%ise(4,jproc)-comgp%ise(3,jproc)+1)
                  if(iproc==jproc) then
                      comgp%nrecvBuf = comgp%nrecvBuf + &
                            (ie3-is3+1)*(comgp%ise(2,jproc)-comgp%ise(1,jproc)+1)*(comgp%ise(4,jproc)-comgp%ise(3,jproc)+1)
                  end if
              else if(ie3j > lzd%Glr%d%n3i .and. lzd%Glr%geocode /= 'F')then
                   stop 'WILL PROBABLY NOT WORK!'
                   ie3j = comgp%ise(6,jproc) - lzd%Glr%d%n3i
                   if(ie3j>=is3k) then
                       is3=max(0,is3k) ! starting index in z dimension for data to be sent
                       ie3=min(ie3j,ie3k) ! ending index in z dimension for data to be sent
                       ioffsetz=is3-0 ! starting index (in z direction) of data to be sent (actually it is the index -1)
                       ioverlap=ioverlap+1
                       !tag=tag+1
                       !!tag=p2p_tag(jproc)
                       if(is3<is3min .or. ioverlap==1) then
                           is3min=is3
                       end if
                       if(ie3>ie3max .or. ioverlap==1) then
                           ie3max=ie3
                       end if
                       !!call setCommunicationPotential(kproc, is3, ie3, ioffsetz, lzd%Glr%d%n1i, lzd%Glr%d%n2i, jproc,&
                       !!     istdest, tag, comgp%comarr(1,ioverlap,jproc))
                       istsource=ioffsetz*lzd%glr%d%n1i*lzd%glr%d%n2i+1
                       !ncount=(ie3-is3+1)*lzd%glr%d%n1i*lzd%glr%d%n2i
                       ncount=lzd%glr%d%n1i*lzd%glr%d%n2i
                       call setCommsParameters(kproc, jproc, istsource, istdest, ncount, tag, comgp%comarr(1,ioverlap,jproc))
                       comgp%comarr(7,ioverlap,jproc)=(ie3-is3+1)
                       comgp%comarr(8,ioverlap,jproc)=lzd%glr%d%n1i*lzd%glr%d%n2i
                       istdest = istdest + (ie3-is3+1)*ncount
                       if(iproc==jproc) then
                           comgp%nrecvBuf = comgp%nrecvBuf + (ie3-is3+1)*lzd%Glr%d%n1i*lzd%Glr%d%n2i
                       end if
                   end if
                   if(is3j <= ie3k)then
                       is3=max(is3j,is3k) ! starting index in z dimension for data to be sent
                       ie3=min(lzd%Glr%d%n3i,ie3k) ! ending index in z dimension for data to be sent
                       ioffsetz=is3-is3k ! starting index (in z direction) of data to be sent (actually it is the index -1)
                       ioverlap=ioverlap+1
                       !tag=tag+1
                       !!tag=p2p_tag(jproc)
                       if(is3<is3min .or. ioverlap==1) then
                           is3min=is3
                       end if
                       if(ie3>ie3max .or. ioverlap==1) then
                           ie3max=ie3
                       end if
                       !!call setCommunicationPotential(kproc, is3, ie3, ioffsetz, lzd%Glr%d%n1i, lzd%Glr%d%n2i, jproc,&
                       !!     istdest, tag, comgp%comarr(1,ioverlap,jproc))
                       istsource=ioffsetz*lzd%glr%d%n1i*lzd%glr%d%n2i+1
                       !ncount=(ie3-is3+1)*lzd%glr%d%n1i*lzd%glr%d%n2i
                       ncount=lzd%glr%d%n1i*lzd%glr%d%n2i
                       call setCommsParameters(kproc, jproc, istsource, istdest, ncount, tag, comgp%comarr(1,ioverlap,jproc))
                       comgp%comarr(7,ioverlap,jproc)=ie3-is3+1
                       comgp%comarr(8,ioverlap,jproc)=lzd%glr%d%n1i*lzd%glr%d%n2i
                       istdest = istdest + (ie3-is3+1)*ncount
                       if(iproc==jproc) then
                           comgp%nrecvBuf = comgp%nrecvBuf + (ie3-is3+1)*lzd%Glr%d%n1i*lzd%Glr%d%n2i
                       end if
                   end if
              end if
          end do
          !!comgp%ise3(1,jproc)=is3min
          !!comgp%ise3(2,jproc)=ie3max
          !!if (iproc==0) write(*,*) 'is3min,comgp%ise(5,jproc)', is3min,comgp%ise(5,jproc)
          !!if (iproc==0) write(*,*) 'ie3max,comgp%ise(6,jproc)', ie3max,comgp%ise(6,jproc)
          !if (comgp%ise(5,jproc)/=is3min) stop 'ERROR 1'
          !if (comgp%ise(6,jproc)/=ie3max) stop 'ERROR 2'
          if(ioverlap/=comgp%noverlaps(jproc)) stop 'ioverlap/=comgp%noverlaps(jproc)'
      end do

  else nproc_if ! monoproc

      comgp%nrecvbuf = (comgp%ise(2,iproc)-comgp%ise(1,iproc)+1)*(comgp%ise(4,iproc)-comgp%ise(3,iproc)+1)*&
                       (comgp%ise(6,iproc)-comgp%ise(5,iproc)+1)
  
  end if nproc_if

  comgp%nrecvbuf=max(comgp%nrecvbuf,1)
  
  ! To indicate that no communication is going on.
  comgp%communication_complete=.true.
  !!comgp%messages_posted=.false.

  call timing(iproc,'init_commPot  ','OF')

end subroutine initialize_communication_potential



