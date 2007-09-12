
! Fill the arrays occup and spinar
! if iunit /=0 this means that the file occup.dat does exist and it opens
subroutine input_occup(iproc,iunit,nelec,norb,norbu,norbd,nspin,occup,spinar)
  implicit none
! Arguments
  integer, intent(in) :: nelec,nspin,iproc,norb,norbu,norbd,iunit
  real(kind=8), intent(out) :: occup(norb),spinar(norb)
! Local variables
  integer :: iorb,nt,ne,it,ierror,iorb1
  real(kind=8) :: rocc

  do iorb=1,norb
     spinar(iorb)=1.0d0
  end do
  if (nspin/=1) then
     do iorb=1,norbu
        spinar(iorb)=1.0d0
     end do
     do iorb=norbu+1,norb
        spinar(iorb)=-1.0d0
     end do
  end if
!  write(*,'(1x,a,5i4,30f6.2)')'Spins: ',norb,norbu,norbd,norbup,norbdp,(spinar(iorb),iorb=1,norb)

! First fill the occupation numbers by default
  nt=0
  if (nspin==1) then
     ne=(nelec+1)/2
     do iorb=1,ne
        it=min(2,nelec-nt)
        occup(iorb)=real(it,kind=8)
        nt=nt+it
     enddo
     do iorb=ne+1,norb
        occup(iorb)=0.d0
     end do
  else
     do iorb=1,norb
        it=min(1,nelec-nt)
        occup(iorb)=real(it,kind=8)
        nt=nt+it
     enddo
  end if

! Then read the file "occup.dat" if does exist
  if (iunit /= 0) then
     nt=0
     do
        read(unit=iunit,fmt=*,iostat=ierror) iorb,rocc
        if (ierror/=0) then
           exit
        else
           nt=nt+1
           if (iorb<0 .or. iorb>norb) then
              if (iproc==0) then
                 write(*,'(1x,a,i0,a)') 'ERROR in line ',nt+1,' of the file "occup.dat"'
                 write(*,'(10x,a,i0,a)')     'The orbital index ',iorb,' is incorrect'
              end if
              stop
           elseif (rocc<0.d0 .or. rocc>2.d0) then
              if (iproc==0) then
                 write(*,'(1x,a,i0,a)') 'ERROR in line ',nt+1,' of the file "occup.dat"'
                 write(*,'(10x,a,f5.2,a)')     'The occupation number ',rocc,' is not between 0. and 2.'
              end if
              stop
           else
              occup(iorb)=rocc
           end if
        end if
     end do
     if (iproc==0) then
        write(*,'(1x,a,i0,a)') &
             'The occupation numbers are read from the file "occup.dat" (',nt,' lines read)'
     end if
     close(unit=iunit)
     !Check if sum(occup)=nelec
     rocc=sum(occup)
     if (abs(rocc-real(nelec,kind=8))>1.d-6) then
        if (iproc==0) then
           write(*,'(1x,a,f13.6,a,i0)') 'From the file "occup.dat", the total number of electrons ',rocc,&
                          ' is not equal to ',nelec
        end if
        stop
     end if
  end if
  if (iproc.eq.0) then 
     write(*,'(1x,a,i8)') &
          'Total Number of  Orbitals ',norb
     iorb1=1
     rocc=occup(1)
     do iorb=1,norb
        if (occup(iorb) /= rocc) then
           if (iorb1 == iorb-1) then
              write(*,'(1x,a,i0,a,f6.4)') 'occup(',iorb1,')= ',rocc
           else
              write(*,'(1x,a,i0,a,i0,a,f6.4)') 'occup(',iorb1,':',iorb-1,')= ',rocc
           end if
           rocc=occup(iorb)
           iorb1=iorb
        end if
     enddo
     if (iorb1 == norb) then
        write(*,'(1x,a,i0,a,f6.4)') 'occup(',norb,')= ',occup(norb)
     else
        write(*,'(1x,a,i0,a,i0,a,f6.4)') 'occup(',iorb1,':',norb,')= ',occup(norb)
     end if
  endif

end subroutine input_occup
