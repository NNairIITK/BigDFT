
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

subroutine read_system_variables(iproc,ntypes,atomnames,psppar,radii_cf,&
     npspcode,iasctype,nelpsp,nzatom)
  implicit none
  integer, intent(in) :: iproc,ntypes
  character(len=20), dimension(ntypes), intent(in) :: atomnames
  integer, dimension(ntypes), intent(out) :: npspcode,iasctype,nelpsp,nzatom
  real(kind=8), dimension(ntypes,2), intent(out) :: radii_cf
  real(kind=8), dimension(0:4,0:6,ntypes), intent(out) :: psppar
  !local variables
  character(len=2) :: symbol
  character(len=27) :: filename
  real(kind=8) :: rcov,rprb,ehomo,radfine
  integer :: i,j,l,nlterms,nprl,nn,ityp,ierror,i_stat,i_all
  integer, dimension(:,:), allocatable :: neleconf

  if (iproc == 0) then
     write(*,'(1x,a)')&
          'Atom Name   Ext.Electrons  PSP Code  Radii: Coarse     Fine   Calculated   From File'
  end if

  allocate(neleconf(6,0:3),stat=i_stat)
  call memocc(i_stat,product(shape(neleconf))*kind(neleconf),'neleconf','read_PSP_variables')


  do ityp=1,ntypes
     filename = 'psppar.'//atomnames(ityp)
     ! if (iproc.eq.0) write(*,*) 'opening PSP file ',filename
     open(unit=11,file=filename,status='old',iostat=ierror)
     !Check the open statement
     if (ierror /= 0) then
        write(*,*) 'iproc=',iproc,': Failed to open the file (it must be in ABINIT format!) "',&
             trim(filename),'"'
        stop
     end if
     read(11,*)
     read(11,*) nzatom(ityp),nelpsp(ityp)
     read(11,*) npspcode(ityp)
     psppar(:,:,ityp)=0.d0
     if (npspcode(ityp) == 2) then !GTH case
        read(11,*) (psppar(0,j,ityp),j=0,4)
        do i=1,2
           read(11,*) (psppar(i,j,ityp),j=0,3-i)
        enddo
     else if (npspcode(ityp) == 3) then !HGH case
        read(11,*) (psppar(0,j,ityp),j=0,4)
        read(11,*) (psppar(1,j,ityp),j=0,3)
        do i=2,4
           read(11,*) (psppar(i,j,ityp),j=0,3)
           read(11,*) !k coefficients, not used for the moment (no spin-orbit coupling)
        enddo
     else if (npspcode(ityp) == 10) then !HGH-K case
        read(11,*) psppar(0,0,ityp),nn,(psppar(0,j,ityp),j=1,nn) !local PSP parameters
        read(11,*) nlterms !number of channels of the pseudo
        prjloop: do l=1,nlterms
           read(11,*) psppar(l,0,ityp),nprl,psppar(l,1,ityp),&
                (psppar(l,j+2,ityp),j=2,nprl) !h_ij terms
           do i=2,nprl
              read(11,*) psppar(l,i,ityp),(psppar(l,i+j+1,ityp),j=i+1,nprl) !h_ij terms
           end do
           if (l==1) cycle
           do i=1,nprl
              read(11,*) !k coefficients, not used
           end do
        end do prjloop
     else
        if (iproc == 0) then
           write(*,'(1x,a,a)')trim(atomnames(ityp)),&
                'unrecognized pspcode: only GTH, HGH & HGH-K pseudos (ABINIT format)'
        end if
        stop
       end if
       !see whether the atom is semicore or not
       call eleconf(nzatom(ityp),nelpsp(ityp),symbol,rcov,rprb,ehomo,neleconf,iasctype(ityp))
       !if you want no semicore electrons, uncomment the following line
       !iasctype(ityp)=0

       !old way of calculating the radii, requires modification of the PSP files
       read(11,*,iostat=ierror) radii_cf(ityp,1),radii_cf(ityp,2)
       if (ierror.eq.0) then
          if (iproc==0) write(*,'(3x,a6,13x,i3,5x,i3,10x,2(1x,f8.5),a)')&
               trim(atomnames(ityp)),nelpsp(ityp),npspcode(ityp),&
               radii_cf(ityp,1),radii_cf(ityp,2),&
               '                   X    '
       else
          !assigning the radii by calculating physical parameters
          radii_cf(ityp,1)=1.d0/sqrt(abs(2.d0*ehomo))
          radfine=100.d0
          do i=0,4
             if (psppar(i,0,ityp)/=0.d0) then
                radfine=min(radfine,psppar(i,0,ityp))
             end if
          end do
          radii_cf(ityp,2)=radfine
          if (iproc==0) write(*,'(3x,a6,13x,i3,5x,i3,10x,2(1x,f8.5),a)')&
               trim(atomnames(ityp)),nelpsp(ityp),npspcode(ityp),&
               radii_cf(ityp,1),radii_cf(ityp,2),&
               '       X                '
       end if
       close(11)
    enddo

    !deallocation
    i_all=-product(shape(neleconf))*kind(neleconf)
    deallocate(neleconf,stat=i_stat)
    call memocc(i_stat,i_all,'neleconf','read_PSP_variables')


    !calculate number of electrons and orbitals


  end subroutine read_system_variables
