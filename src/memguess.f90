program memguess

  !implicit real(kind=8) (a-h,o-z)
  implicit none
  logical :: calc_tail,output_grid
  character(len=20) :: tatonam,units
  character(len=30) :: filename
  character(len=2) :: symbol
  integer :: ierror,nat,ntypes,iat,jat,ityp,nproc,n1,n2,n3,ixc,ncharge,itermax,i_stat,i_all,i,j
  integer :: ncong,ncongt,idsx,nzatom,npspcode,iasctype,nelec,norb,nateq,nt,it,iorb,iorb1,ne
  logical :: exists
  real(kind=8) :: hgrid,crmult,frmult,cpmult,fpmult,gnrm_cv,rbuf,elecfield
  real(kind=8) :: alat1,alat2,alat3,rcov,rprb,ehomo,radfine,rocc
  real(kind=8) :: cxmin,cxmax,cymin,cymax,czmin,czmax
  character(len=20), dimension(:), allocatable :: atomnames
  integer, dimension(:), allocatable :: iatype,nelpsp
  integer, dimension(:,:), allocatable :: neleconf
  real(kind=8), dimension(:,:), allocatable :: rxyz,radii_cf
  real(kind=8), dimension(:,:,:), allocatable :: psppar
  real(kind=8), dimension(:), allocatable :: occup 

  call getarg(1,tatonam)

  if(trim(tatonam)=='') then
     write(*,'(1x,a)')&
          'Usage: ./memguess <nproc> [y]'
     write(*,'(1x,a)')&
          'Indicate the number of processes after the executable'
     write(*,'(1x,a)')&
          'You can put a "y" in the second argument (optional) if you want the '
     write(*,'(1x,a)')&
          '  grid to be plotted with V_Sim'
     stop
  else
     read(unit=tatonam,fmt=*) nproc
     call getarg(2,tatonam)
     if(trim(tatonam)=='') then
        output_grid=.false.
     else if (trim(tatonam)=='y') then
        output_grid=.true.
        write(*,'(1x,a)')&
             'The system grid will be displayed in the "grid.ascii" file'
     else
        write(*,'(1x,a)')&
             'Usage: ./memguess <nproc> [y]'
        write(*,'(1x,a)')&
             'Indicate the number of processes after the executable'
        write(*,'(1x,a)')&
             'ERROR: The only second argument which is accepted is "y"'
        stop
     end if
  end if

  !initialize memory counting
  call memocc(0,0,'count','start')

  ! read atomic positions
  open(unit=9,file='posinp',status='old')
  read(9,*) nat,units
  allocate(rxyz(3,nat),stat=i_stat)
  call memocc(i_stat,product(shape(rxyz))*kind(rxyz),'rxyz','memguess')
  allocate(iatype(nat),stat=i_stat)
  call memocc(i_stat,product(shape(iatype))*kind(iatype),'iatype','memguess')
  allocate(atomnames(100),stat=i_stat)
  call memocc(i_stat,product(shape(atomnames))*kind(atomnames),'atomnames','memguess')

  ntypes=0
  do iat=1,nat
     read(9,*) rxyz(1,iat),rxyz(2,iat),rxyz(3,iat),tatonam
     do ityp=1,ntypes
        if (tatonam.eq.atomnames(ityp)) then
           iatype(iat)=ityp
           goto 200
        endif
     enddo
     ntypes=ntypes+1
     if (ntypes.gt.100) stop 'more than 100 atomnames not permitted'
     atomnames(ityp)=tatonam
     iatype(iat)=ntypes
200  continue
     if (units.eq.'angstroem') then
        ! if Angstroem convert to Bohr
        do i=1,3 
           rxyz(i,iat)=rxyz(i,iat)/.529177d0  
        enddo
     else if  (units.eq.'atomic' .or. units.eq.'bohr') then
     else
        write(*,*) 'length units in input file unrecognized'
        write(*,*) 'recognized units are angstroem or atomic = bohr'
        stop 
     endif
  enddo
  close(9)
  do ityp=1,ntypes
     write(*,'(1x,a,i0,a,a)') &
          'atoms of type ',ityp,' are ',trim(atomnames(ityp))
  enddo

  !control atom positions
  nateq=0
  do iat=1,nat
     do jat=iat+1,nat
        if ((rxyz(1,iat)-rxyz(1,jat))**2+(rxyz(2,iat)-rxyz(2,jat))**2+&
             (rxyz(3,iat)-rxyz(3,jat))**2 ==0.d0) then
           nateq=nateq+1
           write(*,'(1x,a,2(i0,a,a6,a))')'ERROR: atoms ',iat,&
                ' (',trim(atomnames(iatype(iat))),') and ',&
                jat,' (',trim(atomnames(iatype(jat))),&
                ') have the same positions'
        end if
     end do
  end do
  if (nateq /= 0) then
     write(*,'(1x,a)')'Control your posinp file, cannot proceed'
     stop
  end if

  ! Read the input variables.
  open(unit=1,file='input.dat',status='old')
  !First line for the main routine (the program)
  read(1,*) 
  !Parameters 
  read(1,*) hgrid
  read(1,*) crmult
  read(1,*) frmult
  read(1,*) cpmult
  read(1,*) fpmult
  if (fpmult.gt.frmult) write(*,*) 'NONSENSE: fpmult > frmult'
  read(1,*) ixc
  read(1,*) ncharge,elecfield
  read(1,*) gnrm_cv
  read(1,*) itermax
  read(1,*) ncong
  read(1,*) idsx
  read(1,*) calc_tail
  read(1,*) rbuf
  read(1,*) ncongt
  close(1)
 
  write(*,'(1x,a)')&
       '------------------------------------------------------------------- Input Parameters'
  write(*,'(1x,a)')&
       '    System Choice       Resolution Radii        SCF Iteration      Finite Size Corr.'
  write(*,'(1x,a,f7.3,1x,a,f5.2,1x,a,1pe8.1,1x,a,l4)')&
       'Grid spacing=',hgrid,    '|  Coarse Wfs.=',crmult,'| Wavefns Conv.=',gnrm_cv,&
       '| Calculate=',calc_tail
  write(*,'(1x,a,i7,1x,a,f5.2,1x,a,i8,1x,a,f4.1)')&
       '       XC id=',ixc,      '|    Fine Wfs.=',frmult,'| Max. N. Iter.=',itermax,&
       '| Extension=',rbuf
  write(*,'(1x,a,i7,1x,a,f5.2,1x,a,i8,1x,a,i4)')&
       'total charge=',ncharge,  '| Coarse Proj.=',cpmult,'| CG Prec.Steps=',ncong,&
       '|  CG Steps=',ncongt
  write(*,'(1x,a,1pe7.1,1x,a,0pf5.2,1x,a,i8)')&
       ' elec. field=',elecfield,'|   Fine Proj.=',fpmult,'| DIIS Hist. N.=',idsx

  allocate(psppar(0:4,0:4,ntypes),stat=i_stat)
  call memocc(i_stat,product(shape(psppar))*kind(psppar),'psppar','memguess')
  allocate(nelpsp(ntypes),stat=i_stat)
  call memocc(i_stat,product(shape(nelpsp))*kind(nelpsp),'nelpsp','memguess')
  allocate(radii_cf(ntypes,2),stat=i_stat)
  call memocc(i_stat,product(shape(radii_cf))*kind(radii_cf),'radii_cf','memguess')
  allocate(neleconf(6,0:3),stat=i_stat)
  call memocc(i_stat,product(shape(neleconf))*kind(neleconf),'neleconf','memguess')
  
  write(*,'(1x,a)')&
       '------------------------------------------------------------------ System Properties'
  write(*,'(1x,a)')&
       'Atom Name   Ext.Electrons  PSP Code  Radii: Coarse     Fine   Calculated   From File'

  do ityp=1,ntypes
     filename = 'psppar.'//atomnames(ityp)
     open(unit=11,file=filename,status='old',iostat=ierror)
     !Check the open statement
     if (ierror /= 0) then
        write(*,*)': Failed to open the file (it must be in ABINIT format!) "',&
             trim(filename),'"'
        stop
     end if
     read(11,*)
     read(11,*) nzatom,nelpsp(ityp)
     read(11,*) npspcode
     psppar(:,:,ityp)=0.d0
     read(11,*) (psppar(0,j,ityp),j=0,4)
     if (npspcode == 2) then !GTH case
        do i=1,2
           read(11,*) (psppar(i,j,ityp),j=0,3-i)
        enddo
     else if (npspcode == 3) then !HGH case
        read(11,*) (psppar(1,j,ityp),j=0,3)
        do i=2,4
           read(11,*) (psppar(i,j,ityp),j=0,3)
           read(11,*) !k coefficients, not used (no spin-orbit coupling)
        enddo
     else
        write(*,'(1x,a,a)')trim(atomnames(ityp)),&
             'unrecognized pspcode (accepts only GTH & HGH pseudopotentials in ABINIT format)'
        stop
     end if
     !see whether the atom is semicore or not
     call eleconf(nzatom,nelpsp(ityp),symbol,rcov,rprb,ehomo,neleconf,iasctype)

     !old way of calculating the radii, requires modification of the PSP files
     read(11,*,iostat=ierror) radii_cf(ityp,1),radii_cf(ityp,2)
     if (ierror.eq.0) then
        write(*,'(3x,a6,13x,i3,5x,i3,10x,2(1x,f8.5),a)')&
             trim(atomnames(ityp)),nelpsp(ityp),npspcode,&
             radii_cf(ityp,1),radii_cf(ityp,2),&
             '                   X    '
     else
        !new method for assigning the radii
        radii_cf(ityp,1)=1.d0/sqrt(abs(2.d0*ehomo))
        radfine=100.d0
        do i=0,4
           if (psppar(i,0,ityp)/=0.d0) then
              radfine=min(radfine,psppar(i,0,ityp))
           end if
        end do
        radii_cf(ityp,2)=radfine
        write(*,'(3x,a6,13x,i3,5x,i3,10x,2(1x,f8.5),a)')&
             trim(atomnames(ityp)),nelpsp(ityp),npspcode,&
             radii_cf(ityp,1),radii_cf(ityp,2),&
             '       X                '
     end if
     close(11)
  enddo

  !deallocation
  i_all=-product(shape(neleconf))*kind(neleconf)
  deallocate(neleconf,stat=i_stat)
  call memocc(i_stat,i_all,'neleconf','memguess')
  i_all=-product(shape(psppar))*kind(psppar)
  deallocate(psppar,stat=i_stat)
  call memocc(i_stat,i_all,'psppar','memguess')

! Number of electrons and number of semicore atoms
  nelec=0
  do iat=1,nat
     ityp=iatype(iat)
     nelec=nelec+nelpsp(ityp)
  enddo

  i_all=-product(shape(nelpsp))*kind(nelpsp)
  deallocate(nelpsp,stat=i_stat)
  call memocc(i_stat,i_all,'nelpsp','memguess')

  nelec=nelec-ncharge
  write(*,'(1x,a,i8)') &
       'Total Number of Electrons ',nelec
  if (mod(nelec,2).ne.0) write(*,*) &
       'WARNING: odd number of electrons, no closed shell system'

! Number of orbitals
  norb=(nelec+1)/2

! Test if the file 'occup.dat exists
  inquire(file='occup.dat',exist=exists)
  if (exists) then
     open(unit=24,file='occup.dat',form='formatted',action='read',status='old')
     !The first line gives the number of orbitals
     read(24,*,iostat=ierror) nt
     if (ierror /=0) then
         write(*,'(1x,a)') 'ERROR reading the number of orbitals in the file "occup.dat"'
        stop
     end if
     if (nt<=norb) then
        write(*,'(1x,a,i0,a,i0)') &
                'ERROR: In the file "occup.dat", the number of orbitals norb=',nt,&
                ' should be strictly greater than (nelec+1)/2=',norb
        stop
     else
        norb=nt
     end if
  end if

  if (exists) then
     write(*,'(1x,a,i8,a)') &
          'Total Number of  Orbitals ',norb,' (read from the file "occup.dat")'
  else
     write(*,'(1x,a,i8)') &
          'Total Number of  Orbitals ',norb
  end if

! Read the file "occup.dat" and test it
  if (exists) then
     allocate(occup(norb),stat=i_stat)
     call memocc(i_stat,product(shape(occup))*kind(occup),'occup','memguess')
     ! First fill the occupation numbers by default
     ne=(nelec+1)/2
     nt=0
     do iorb=1,ne
        it=min(2,nelec-nt)
        occup(iorb)=real(it,kind=8)
        nt=nt+it
     enddo
     do iorb=ne+1,norb
        occup(iorb)=0.d0
     end do
     ! Then read the file "occup.dat" if does exist
     nt=0
     do
        read(24,*,iostat=ierror) iorb,rocc
        if (ierror/=0) then
           exit
        else
           nt=nt+1
           if (iorb<0 .or. iorb>norb) then
              write(*,'(1x,a,i0,a)') 'ERROR in line ',nt+1,' of the file "occup.dat"'
              write(*,'(10x,a,i0,a)')     'The orbital index ',iorb,' is incorrect'
              stop
           elseif (rocc<0.d0 .or. rocc>2.d0) then
              write(*,'(1x,a,i0,a)') 'ERROR in line ',nt+1,' of the file "occup.dat"'
              write(*,'(10x,a,f5.2,a)')     'The occupation number ',rocc,' is not between 0. and 2.'
              stop
           else
              occup(iorb)=rocc
           end if
        end if
     end do
     write(*,'(1x,a,i0,a)') &
             'The occupation numbers are read from the file "occup.dat" (',nt,' lines read)'
     close(unit=24)
  end if

  !Check if sum(occup)=nelec
  rocc=sum(occup)
  if (abs(rocc-real(nelec,kind=8))>1.d-6) then
     write(*,'(1x,a,f13.6,a,i0)') 'From the file "occup.dat", the total number of electrons ',rocc,&
             ' is not equal to ',nelec
     stop
  end if
  write(*,'(1x,a,i8)') &
       'Total Number of  Orbitals ',norb
  iorb1=1
  rocc=occup(1)
  do iorb=1,norb
     if (occup(iorb) /= rocc) then
        if (iorb1 == iorb-1) then
           write(*,'(4x,a,i0,a,f6.4)') 'occup(',iorb1,')= ',rocc
        else
           write(*,'(4x,a,i0,a,i0,a,f6.4)') 'occup(',iorb1,':',iorb-1,')= ',rocc
        end if
        rocc=occup(iorb)
        iorb1=iorb
     end if
  enddo
  if (iorb1 == norb) then
     write(*,'(4x,a,i0,a,f6.4)') 'occup(',norb,')= ',occup(norb)
  else
     write(*,'(4x,a,i0,a,i0,a,f6.4)') 'occup(',iorb1,':',norb,')= ',occup(norb)
  end if
! De-allocation of occup
  deallocate(occup,stat=i_stat)
  call memocc(i_stat,i_all,'occup','cluster')

! determine size alat of overall simulation cell
  call system_size(nat,rxyz,radii_cf(1,1),crmult,iatype,ntypes, &
       cxmin,cxmax,cymin,cymax,czmin,czmax)
  alat1=(cxmax-cxmin)
  alat2=(cymax-cymin)
  alat3=(czmax-czmin)

  do iat=1,nat
     rxyz(1,iat)=rxyz(1,iat)-cxmin
     rxyz(2,iat)=rxyz(2,iat)-cymin
     rxyz(3,iat)=rxyz(3,iat)-czmin
  enddo

! grid sizes n1,n2,n3
  n1=int(alat1/hgrid)
  n2=int(alat2/hgrid)
  n3=int(alat3/hgrid)
  alat1=real(n1,kind=8)*hgrid 
  alat2=real(n2,kind=8)*hgrid 
  alat3=real(n3,kind=8)*hgrid
  write(*,'(1x,a,19x,a)') &
       '                                 Atomic Units:','grid spacing units:'
  write(*,'(1x,a,3(1x,1pe12.5),3x,3(1x,i9))')&
       '  Box Sizes=',alat1,alat2,alat3,n1,n2,n3

  call MemoryEstimator(nproc,idsx,n1,n2,n3,alat1,alat2,alat3,hgrid,nat,ntypes,iatype,&
          rxyz,radii_cf,crmult,frmult,norb,atomnames,output_grid)


  i_all=-product(shape(atomnames))*kind(atomnames)
  deallocate(atomnames,stat=i_stat)
  call memocc(i_stat,i_all,'atomnames','memguess')
  i_all=-product(shape(radii_cf))*kind(radii_cf)
  deallocate(radii_cf,stat=i_stat)
  call memocc(i_stat,i_all,'radii_cf','memguess')
  i_all=-product(shape(rxyz))*kind(rxyz)
  deallocate(rxyz,stat=i_stat)
  call memocc(i_stat,i_all,'rxyz','memguess')
  i_all=-product(shape(iatype))*kind(iatype)
  deallocate(iatype,stat=i_stat)
  call memocc(i_stat,i_all,'iatype','memguess')

  !finalize memory counting
  call memocc(0,0,'count','stop')

end program memguess
