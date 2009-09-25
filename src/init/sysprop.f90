!!****f* BigDFT/system_properties
!! FUNCTION
!!  Calculate the important objects related to the physical properties of the system
!!
!! SOURCE
!!
subroutine system_properties(iproc,nproc,in,atoms,orbs,radii_cf,nelec)
  use module_base
  use module_types
  implicit none
  integer, intent(in) :: iproc,nproc
  integer, intent(out) :: nelec
  type(input_variables), intent(in) :: in
  type(atoms_data), intent(inout) :: atoms
  type(orbitals_data), intent(out) :: orbs
  real(gp), dimension(atoms%ntypes,3), intent(out) :: radii_cf
  !local variables
  character(len=*), parameter :: subname='orbitals_descriptors'
  integer :: iunit,norb,norbu,norbd,nspinor,jpst,norbme,norbyou,i_all,i_stat,jproc,ikpts

  call read_system_variables(iproc,nproc,in,atoms,radii_cf,nelec,&
       norb,norbu,norbd,iunit)

  if(in%nspin==4) then
     nspinor=4
  else
     nspinor=1
  end if

  !temporary changement, to be controlled
  !nspinor=2

  call orbitals_descriptors(iproc,nproc,norb,norbu,norbd,nspinor,orbs)

  !distribution of wavefunction arrays between processors
  !tuned for the moment only on the cubic distribution
  if (iproc == 0 .and. nproc > 1) then
     jpst=0
     do jproc=0,nproc-1
        norbme=orbs%norb_par(jproc)
        norbyou=orbs%norb_par(min(jproc+1,nproc-1))
        if (norbme /= norbyou .or. jproc == nproc-1) then
           !this is a screen output that must be modified
           write(*,'(3(a,i0),a)')&
                ' Processes from ',jpst,' to ',jproc,' treat ',norbme,' orbitals '
           jpst=jproc+1
        end if
     end do
     !write(*,'(3(a,i0),a)')&
     !     ' Processes from ',jpst,' to ',nproc-1,' treat ',norbyou,' orbitals '
  end if


  !assign to each k-point the same occupation number
  do ikpts=1,orbs%nkpts
     call input_occup(iproc,iunit,nelec,norb,norbu,norbd,in%nspin,in%mpol,&
          orbs%occup(1+(ikpts-1)*orbs%norb),orbs%spinsgn(1+(ikpts-1)*orbs%norb))
  end do

end subroutine system_properties
!!***


!!****f* BigDFT/read_system_variables
!! FUNCTION
!!   Assign some of the physical system variables
!!   Performs also some cross-checks with other variables
!! DESCRIPTION
!!   The pointer in atoms structure have to be associated or nullify.
!! SOURCE
!!
subroutine read_system_variables(iproc,nproc,in,atoms,radii_cf,&
     nelec,norb,norbu,norbd,iunit)
  use module_base
  use module_types
  implicit none
  type(input_variables), intent(in) :: in
  integer, intent(in) :: iproc,nproc
  type(atoms_data), intent(inout) :: atoms
  integer, intent(out) :: nelec,norb,norbu,norbd,iunit
  real(gp), dimension(atoms%ntypes,3), intent(out) :: radii_cf
  !local variables
  character(len=*), parameter :: subname='read_system_variables'
  !real(kind=8), parameter :: eps_mach=1.d-12
  logical :: exists
  character(len=2) :: symbol
  character(len=24) :: message
  character(len=27) :: filename
  character(len=50) :: format
  character(len=100) :: line
  integer :: i,j,k,l,iat,nlterms,nprl,nn,nt,ntu,ntd,ityp,ierror,i_stat,i_all,ixcpsp,ispinsum,mxpl
  integer :: ispol,mxchg,ichg,natpol,ichgsum,nsccode,ierror1
  real(gp) :: rcov,rprb,ehomo,radfine,minrad,maxrad
  real(gp), dimension(3,3) :: hij
  real(gp), dimension(2,2,3) :: offdiagarr
  integer, dimension(6,0:3) :: neleconf

  !allocate atoms data variables
  ! store PSP parameters
  ! modified to accept both GTH and HGHs pseudopotential types
  allocate(atoms%psppar(0:4,0:6,atoms%ntypes+ndebug),stat=i_stat)
  call memocc(i_stat,atoms%psppar,'atoms%psppar',subname)
  allocate(atoms%nelpsp(atoms%ntypes+ndebug),stat=i_stat)
  call memocc(i_stat,atoms%nelpsp,'atoms%nelpsp',subname)
  allocate(atoms%npspcode(atoms%ntypes+ndebug),stat=i_stat)
  call memocc(i_stat,atoms%npspcode,'atoms%npspcode',subname)
  allocate(atoms%nzatom(atoms%ntypes+ndebug),stat=i_stat)
  call memocc(i_stat,atoms%nzatom,'atoms%nzatom',subname)
  allocate(atoms%iasctype(atoms%ntypes+ndebug),stat=i_stat)
  call memocc(i_stat,atoms%iasctype,'atoms%iasctype',subname)

  if (iproc == 0) then
     write(*,'(1x,a)')&
          ' Atom    N.Electr.  PSP Code  Radii: Coarse     Fine  CoarsePSP    Calculated   File'
  end if

  do ityp=1,atoms%ntypes
     filename = 'psppar.'//atoms%atomnames(ityp)

     inquire(file=filename,exist=exists)
     if (.not. exists) then
        !if (iproc == 0) 
            write(*,'(1x,3a)')&
             'ERROR: The pseudopotential parameter file "',trim(filename),&
             '" is lacking, exiting...'
        stop
     end if
     ! if (iproc.eq.0) write(*,*) 'opening PSP file ',filename
     open(unit=11,file=filename,status='old',iostat=ierror)
     !Check the open statement
     if (ierror /= 0) then
        write(*,*) 'iproc=',iproc,&
             ': Failed to open the file (it must be in ABINIT format!): "',&
             trim(filename),'"'
        stop
     end if
     read(11,*)
     read(11,*) atoms%nzatom(ityp),atoms%nelpsp(ityp)
     read(11,*) atoms%npspcode(ityp),ixcpsp
     !control if the PSP is calculated with the same XC value
     if (ixcpsp /= in%ixc .and. iproc==0) then
        write(*,'(1x,a)')&
             'WARNING: The pseudopotential file "'//trim(filename)//'"'
        write(*,'(1x,a,i0,a,i0)')&
             '         contains a PSP generated with an XC id=',&
             ixcpsp,' while for this run ixc=',in%ixc
     end if
     atoms%psppar(:,:,ityp)=0._gp
     if (atoms%npspcode(ityp) == 2) then !GTH case
        read(11,*) (atoms%psppar(0,j,ityp),j=0,4)
        do i=1,2
           read(11,*) (atoms%psppar(i,j,ityp),j=0,3-i)
        enddo
     else if (atoms%npspcode(ityp) == 3) then !HGH case
        read(11,*) (atoms%psppar(0,j,ityp),j=0,4)
        read(11,*) (atoms%psppar(1,j,ityp),j=0,3)
        do i=2,4
           read(11,*) (atoms%psppar(i,j,ityp),j=0,3)
           read(11,*) !k coefficients, not used for the moment (no spin-orbit coupling)
        enddo
     else if (atoms%npspcode(ityp) == 10) then !HGH-K case
        read(11,*) atoms%psppar(0,0,ityp),nn,(atoms%psppar(0,j,ityp),j=1,nn) !local PSP parameters
        read(11,*) nlterms !number of channels of the pseudo
        prjloop: do l=1,nlterms
           read(11,*) atoms%psppar(l,0,ityp),nprl,atoms%psppar(l,1,ityp),&
                (atoms%psppar(l,j+2,ityp),j=2,nprl) !h_ij terms
           do i=2,nprl
              read(11,*) atoms%psppar(l,i,ityp),(atoms%psppar(l,i+j+1,ityp),j=i+1,nprl) !h_ij 
           end do
           if (l==1) cycle
           do i=1,nprl
              read(11,*) !k coefficients, not used
           end do
        end do prjloop
     else
        !if (iproc == 0) then
           write(*,'(1x,a,a)')trim(atoms%atomnames(ityp)),&
                'unrecognized pspcode: only GTH, HGH & HGH-K pseudos (ABINIT format)'
        !end if
        stop
     end if
     !see whether the atom is semicore or not
     call eleconf(atoms%nzatom(ityp),atoms%nelpsp(ityp),symbol,rcov,rprb,ehomo,&
          neleconf,atoms%iasctype(ityp),mxpl,mxchg,atoms%amu(ityp))
     !if you want no semicore input guess electrons, uncomment the following line
     !atoms%iasctype(ityp)=0

     !here we must check of the input guess polarisation
     !control if the values are compatible with the atom configuration
     !do this for all atoms belonging to a given type
     !control the maximum polarisation allowed: consider only non-closed shells   

     do iat=1,atoms%nat
        if (atoms%iatype(iat) == ityp) then
           call charge_and_spol(atoms%natpol(iat),ichg,ispol)
           if (abs(ispol) > mxpl) then
              !if (iproc ==0) 
                      write(*,'(1x,a,i0,a,a,2(a,i0))')&
                   'ERROR: Input polarisation of atom No.',iat,&
                   ' (',trim(atoms%atomnames(ityp)),') must be <=',mxpl,&
                   ', while found ',ispol
              stop
           end if
           if (abs(ichg) > mxchg) then
              !if (iproc ==0) 
                   write(*,'(1x,a,i0,a,a,2(a,i0))')&
                   'ERROR: Input charge of atom No.',iat,&
                   ' (',trim(atoms%atomnames(ityp)),') must be <=',mxchg,&
                   ', while found ',ichg
              stop
           end if
        end if
     end do

     !control the hardest and the softest gaussian
     minrad=1.e10_gp
     maxrad=0.e0_gp ! This line added by Alexey, 03.10.08, to be able to compile with -g -C
     do i=0,4
        !the maximum radii is useful only for projectors
        if (i==1) maxrad=0.0_gp
        if (atoms%psppar(i,0,ityp)/=0._gp) then
           minrad=min(minrad,atoms%psppar(i,0,ityp))
           maxrad=max(maxrad,atoms%psppar(i,0,ityp))
        end if
     end do

     !old way of calculating the radii, requires modification of the PSP files
     read(11,'(a100)',iostat=ierror)line
     if (ierror /=0) then
        !if (iproc ==0) write(*,*)&
        !     ' WARNING: last line of pseudopotential missing, put an empty line'
        line=''
     end if
     read(line,*,iostat=ierror1) radii_cf(ityp,1),radii_cf(ityp,2),radii_cf(ityp,3)
     if (ierror1 /= 0 ) then
        read(line,*,iostat=ierror) radii_cf(ityp,1),radii_cf(ityp,2)
        radii_cf(ityp,3)=radii_cf(ityp,2)
     end if
     message='                   X ' 

     if (ierror /= 0) then
        !assigning the radii by calculating physical parameters
        radii_cf(ityp,1)=1._gp/sqrt(abs(2._gp*ehomo))
        radfine=100._gp
        do i=0,4
           if (atoms%psppar(i,0,ityp)/=0._gp) then
              radfine=min(radfine,atoms%psppar(i,0,ityp))
           end if
        end do
        radii_cf(ityp,2)=radfine
        message='         X              '
        radii_cf(ityp,3)=radfine
     end if
     close(11)

     !correct the coarse and the fine radius for projectors
     radii_cf(ityp,3)=max(min(in%crmult*radii_cf(ityp,1),15.0_gp*maxrad)/in%frmult,radii_cf(ityp,2))

     if (maxrad == 0.0_gp) then
        radii_cf(ityp,3)=0.0_gp
     end if

     if (iproc==0) write(*,'(1x,a6,8x,i3,5x,i3,10x,3(1x,f8.5),a)')&
          trim(atoms%atomnames(ityp)),atoms%nelpsp(ityp),atoms%npspcode(ityp),&
          radii_cf(ityp,1),radii_cf(ityp,2),radii_cf(ityp,3),message

     !control whether the grid spacing is too high
     if (iproc == 0 .and. max(in%hx,in%hy,in%hz) > 2.5_gp*minrad) then
        write(*,'(1x,a)')&
             'WARNING: The grid spacing value may be too high to treat correctly the above pseudo.' 
        write(*,'(1x,a,f5.2,a)')&
             '         Results can be meaningless if hgrid is bigger than',2.5_gp*minrad,&
             '. At your own risk!'
     end if

  enddo

  !print the pseudopotential matrices
  if (iproc == 0) then
     do l=1,3
        do i=1,2
           do j=i+1,3
              offdiagarr(i,j-i,l)=0._gp
              if (l==1) then
                 if (i==1) then
                    if (j==2)   offdiagarr(i,j-i,l)=-0.5_gp*sqrt(3._gp/5._gp)
                    if (j==3)   offdiagarr(i,j-i,l)=0.5_gp*sqrt(5._gp/21._gp)
                 else
                    offdiagarr(i,j-i,l)=-0.5_gp*sqrt(100._gp/63._gp)
                 end if
              else if (l==2) then
                 if (i==1) then
                    if (j==2)   offdiagarr(i,j-i,l)=-0.5_gp*sqrt(5._gp/7._gp)
                    if (j==3)   offdiagarr(i,j-i,l)=1._gp/6._gp*sqrt(35._gp/11._gp)
                 else
                    offdiagarr(i,j-i,l)=-7._gp/3._gp*sqrt(1._gp/11._gp)
                 end if
              else if (l==3) then
                 if (i==1) then
                    if (j==2)   offdiagarr(i,j-i,l)=-0.5_gp*sqrt(7._gp/9._gp)
                    if (j==3)   offdiagarr(i,j-i,l)=0.5_gp*sqrt(63._gp/143._gp)
                 else
                    offdiagarr(i,j-i,l)=-9._gp*sqrt(1._gp/143._gp)
                 end if
              end if
           end do
        end do
     end do

     write(*,'(1x,a)')&
          '------------------------------------ Pseudopotential coefficients (Upper Triangular)'
     do ityp=1,atoms%ntypes
        write(*,'(1x,a)')&
             'Atom Name    rloc      C1        C2        C3        C4  '
        do l=0,4
           if (l==0) then
              do i=4,0,-1
                 j=i
                 if (atoms%psppar(l,i,ityp) /= 0._gp) exit
              end do
              write(*,'(3x,a6,5(1x,f9.5))')&
                   trim(atoms%atomnames(ityp)),(atoms%psppar(l,i,ityp),i=0,j)
           else
              do i=3,0,-1
                 j=i
                 if (atoms%psppar(l,i,ityp) /= 0._gp) exit
              end do
              if (j /=0) then
                 write(*,'(1x,a,i0,a)')&
                      '    l=',l-1,' '//'     rl        h1j       h2j       h3j '
                 hij=0._gp
                 do i=1,j
                    hij(i,i)=atoms%psppar(l,i,ityp)
                 end do
                 if (atoms%npspcode(ityp) == 3) then !traditional HGH convention
                    hij(1,2)=offdiagarr(1,1,l)*atoms%psppar(l,2,ityp)
                    hij(1,3)=offdiagarr(1,2,l)*atoms%psppar(l,3,ityp)
                    hij(2,3)=offdiagarr(2,1,l)*atoms%psppar(l,3,ityp)
                 else if (atoms%npspcode(ityp) == 10) then !HGH-K convention
                    hij(1,2)=atoms%psppar(l,4,ityp)
                    hij(1,3)=atoms%psppar(l,5,ityp)
                    hij(2,3)=atoms%psppar(l,6,ityp)
                 end if
                 do i=1,j
                    if (i==1) then
                       write(format,'(a,2(i0,a))')"(9x,(1x,f9.5),",j,"(1x,f9.5))"
                       write(*,format)atoms%psppar(l,0,ityp),(hij(i,k),k=i,j)
                    else
                       write(format,'(a,2(i0,a))')"(19x,",i-1,"(10x),",j-i+1,"(1x,f9.5))"
                       write(*,format)(hij(i,k),k=i,j)
                    end if

                 end do
              end if
           end if
        end do
     end do
  end if

  !calculate number of electrons and orbitals
  ! Number of electrons and number of semicore atoms
  nelec=0
  atoms%natsc=0
  do iat=1,atoms%nat
     ityp=atoms%iatype(iat)
     nelec=nelec+atoms%nelpsp(ityp)
     nsccode=atoms%iasctype(ityp)
     call charge_and_spol(atoms%natpol(iat),ichg,ispol)
     if (ichg /=0) then
        call eleconf(atoms%nzatom(ityp),atoms%nelpsp(ityp),symbol,rcov,rprb,ehomo,&
             neleconf,atoms%iasctype(ityp),mxpl,mxchg,atoms%amu(ityp))
        call correct_semicore(atoms%atomnames(ityp),6,3,ichg,neleconf,nsccode)
     end if
     if (nsccode/= 0) atoms%natsc=atoms%natsc+1
  enddo
  nelec=nelec-in%ncharge
  if (iproc == 0) then
     write(*,'(1x,a,t28,i8)') 'Total Number of Electrons',nelec
  end if

  ! Number of orbitals
  if (in%nspin==1) then
     norb=(nelec+1)/2
     norbu=norb
     norbd=0
     if (mod(nelec,2).ne.0 .and. iproc==0) then
        write(*,'(1x,a)') 'WARNING: odd number of electrons, no closed shell system'
     end if
  else if(in%nspin==4) then
     if (iproc==0) write(*,'(1x,a)') 'Spin-polarized non-collinear calculation'
     norb=nelec
     norbu=norb
     norbd=0
  else 
     if (iproc==0) write(*,'(1x,a)') 'Spin-polarized calculation'
     norb=nelec
     norbu=min(norb/2+in%mpol,norb)
     norbd=norb-norbu

     !test if the spin is compatible with the input guess polarisations
     ispinsum=0
     ichgsum=0
     do iat=1,atoms%nat
        call charge_and_spol(atoms%natpol(iat),ichg,ispol)
        ispinsum=ispinsum+ispol
        ichgsum=ichgsum+ichg
     end do

     if (in%nspin == 2 .and. ispinsum /= norbu-norbd) then
        !if (iproc==0) then 
           write(*,'(1x,a,i0,a)')&
                'ERROR: Total input polarisation (found ',ispinsum,&
                ') must be equal to norbu-norbd.'
           write(*,'(1x,3(a,i0))')&
                'With norb=',norb,' and mpol=',in%mpol,' norbu-norbd=',norbu-norbd
        !end if
        stop
     end if

     if (ichgsum /= in%ncharge .and. ichgsum /= 0) then
        !if (iproc==0) then 
           write(*,'(1x,a,i0,a)')&
                'ERROR: Total input charge (found ',ichgsum,&
                ') cannot be different than charge.'
           write(*,'(1x,2(a,i0))')&
                'The charge is=',in%ncharge,' input charge=',ichgsum
        !end if
        stop
     end if

     !now warn if there is no input guess spin polarisation
     ispinsum=0
     do iat=1,atoms%nat
        call charge_and_spol(atoms%natpol(iat),ichg,ispol)
        ispinsum=ispinsum+abs(ispol)
     end do
     if (ispinsum == 0 .and. in%nspin==2) then
        !if (iproc==0) 
            write(*,'(1x,a)')&
             'WARNING: Found no input polarisation, add it for a correct input guess'
        !stop
     end if

  end if


  ! Test if the file 'occup.dat exists
  inquire(file='occup.dat',exist=exists)
  iunit=0
  if (exists) then
     iunit=25
     open(unit=iunit,file='occup.dat',form='formatted',action='read',status='old')
     if (in%nspin==1) then
        !The first line gives the number of orbitals
        read(unit=iunit,fmt=*,iostat=ierror) nt
     else
        !The first line gives the number of orbitals
        read(unit=iunit,fmt=*,iostat=ierror) ntu,ntd
     end if
     if (ierror /=0) then
        !if (iproc==0) 
          write(*,'(1x,a)') &
             'ERROR: reading the number of orbitals in the file "occup.dat"'
        stop
     end if
     !Check
     if (in%nspin==1) then
        if (nt<norb) then
           !if (iproc==0) 
               write(*,'(1x,a,i0,a,i0)') &
                'ERROR: In the file "occup.dat", the number of orbitals norb=',nt,&
                ' should be greater or equal than (nelec+1)/2=',norb
           stop
        else
           norb=nt
           norbu=norb
           norbd=0
        end if
     else
        nt=ntu+ntd
        if (nt<norb) then
           !if (iproc==0) 
               write(*,'(1x,a,i0,a,i0)') &
                'ERROR: In the file "occup.dat", the number of orbitals norb=',nt,&
                ' should be greater or equal than nelec=',norb
           stop
        else
           norb=nt
        end if
        if (ntu<norbu) then
           !if (iproc==0) 
                write(*,'(1x,a,i0,a,i0)') &
                'ERROR: In the file "occup.dat", the number of orbitals up norbu=',ntu,&
                ' should be greater or equal than min(nelec/2+mpol,nelec)=',norbu
           stop
        else
           norbu=ntu
        end if
        if (ntd<norbd) then
           !if (iproc==0) 
                  write(*,'(1x,a,i0,a,i0)') &
                'ERROR: In the file "occup.dat", the number of orbitals down norbd=',ntd,&
                ' should be greater or equal than min(nelec/2-mpol,0)=',norbd
           stop
        else
           norbd=ntd
        end if
     end if
  end if

!!$  tt=dble(norb)/dble(nproc)
!!$  norbp=int((1.d0-eps_mach*tt) + tt)
!!$  !if (iproc.eq.0) write(*,'(1x,a,1x,i0)') 'norbp=',norbp

end subroutine read_system_variables
!!***


!!****f* BigDFT/orbitals_descriptors
!! FUNCTION
!!    Define the descriptors of the orbitals from a given norb
!!    It uses the cubic strategy for partitioning the orbitals
!! SOURCE
!!
subroutine orbitals_descriptors(iproc,nproc,norb,norbu,norbd,nspinor,orbs)
  use module_base
  use module_types
  implicit none
  integer, intent(in) :: iproc,nproc,nspinor,norb,norbu,norbd
  type(orbitals_data), intent(out) :: orbs
  !local variables
  character(len=*), parameter :: subname='orbitals_descriptors'
  integer :: iorb,jproc,norb_tot,ikpt,i_stat,jorb,ierr,i_all
  real(gp) :: kx,alat
  logical, dimension(:), allocatable :: GPU_for_orbs

  allocate(orbs%norb_par(0:nproc-1+ndebug),stat=i_stat)
  call memocc(i_stat,orbs%norb_par,'orbs%norb_par',subname)

  !assign the value of the k-points
  orbs%nkpts=1!3
  !allocate vectors related to k-points
  allocate(orbs%kpts(3,orbs%nkpts+ndebug),stat=i_stat)
  call memocc(i_stat,orbs%kpts,'orbs%kpts',subname)
  allocate(orbs%kwgts(orbs%nkpts+ndebug),stat=i_stat)
  call memocc(i_stat,orbs%kwgts,'orbs%kwgts',subname)
  !only the gamma point for the moment

!!$  open(55)
!!$  read(55,*)kx,alat
!!$  close(55)

!!$  !toto point
!!$  alat=10.1901d0
!!$
!!$  kx=1.23456000d-01
!!$  ky=8.52147000d-01
!!$  kz=9.87452000d-01

  do ikpt=1,orbs%nkpts
     orbs%kpts(1,ikpt)=0.0_gp
     orbs%kpts(2,ikpt)=0.0_gp
     orbs%kpts(3,ikpt)=0.0_gp
!!$     orbs%kpts(1,ikpt)=kx*4.0*datan(1.d0)/(alat)
!!$     orbs%kpts(2,ikpt)=0.d0
!!$     orbs%kpts(3,ikpt)=0.d0
     orbs%kwgts(ikpt)=1.0_gp
  end do

!!$  orbs%kwgts(1)=0.2_gp
!!$  orbs%kwgts(2)=0.5_gp
!!$  orbs%kwgts(3)=0.3_gp

  !initialise the array
  do jproc=0,nproc-1
     orbs%norb_par(jproc)=0 !size 0 nproc-1
  end do


  !create an array which indicate which processor has a GPU associated 
  !from the viewpoint of the BLAS routines
  allocate(GPU_for_orbs(0:nproc-1+ndebug),stat=i_stat)
  call memocc(i_stat,GPU_for_orbs,'GPU_for_orbs',subname)

  if (nproc > 1 .and. .not. GPUshare) then
     call MPI_ALLGATHER(GPUconv,1,MPI_LOGICAL,GPU_for_orbs(0),1,MPI_LOGICAL,&
          MPI_COMM_WORLD,ierr)
  else
     GPU_for_orbs(0)=GPUconv
  end if

  !cubic-code strategy: balance the orbitals between processors
  !in the most symmetric way
  do iorb=1,norb*orbs%nkpts
     jproc=mod(iorb-1,nproc)
     orbs%norb_par(jproc)=orbs%norb_par(jproc)+1
  end do

  i_all=-product(shape(GPU_for_orbs))*kind(GPU_for_orbs)
  deallocate(GPU_for_orbs,stat=i_stat)
  call memocc(i_stat,i_all,'GPU_for_orbs',subname)


  !check the distribution
  norb_tot=0
  do jproc=0,iproc-1
     norb_tot=norb_tot+orbs%norb_par(jproc)
  end do
  !reference orbital for process
  orbs%isorb=norb_tot
  do jproc=iproc,nproc-1
     norb_tot=norb_tot+orbs%norb_par(jproc)
  end do

  if(norb_tot /= norb*orbs%nkpts) then
     write(*,*)'ERROR: partition of orbitals incorrect'
     stop
  end if

  !assign the values of the orbitals data
  orbs%norb=norb
  orbs%norbp=orbs%norb_par(iproc)
  orbs%nspinor=nspinor
  orbs%norbu=norbu
  orbs%norbd=norbd

  allocate(orbs%iokpt(orbs%norbp+ndebug),stat=i_stat)
  call memocc(i_stat,orbs%iokpt,'orbs%iokpt',subname)

  !assign the k-point to the given orbital, counting one orbital after each other
  jorb=0
  do ikpt=1,orbs%nkpts
     do iorb=1,orbs%norb
        jorb=jorb+1 !this runs over norb*nkpts values
        if (jorb > orbs%isorb .and. jorb <= orbs%isorb+orbs%norbp) then
           orbs%iokpt(jorb-orbs%isorb)=ikpt
        end if
     end do
  end do

  !allocate occupation number and spinsign
  !fill them in normal way
  allocate(orbs%occup(orbs%norb*orbs%nkpts+ndebug),stat=i_stat)
  call memocc(i_stat,orbs%occup,'orbs%occup',subname)
  allocate(orbs%spinsgn(orbs%norb*orbs%nkpts+ndebug),stat=i_stat)
  call memocc(i_stat,orbs%spinsgn,'orbs%spinsgn',subname)
  orbs%occup(1:orbs%norb*orbs%nkpts)=1.0_gp 
  orbs%spinsgn(1:orbs%norb*orbs%nkpts)=1.0_gp

  !allocate the array which assign the k-point to processor in transposed version
  allocate(orbs%ikptproc(orbs%nkpts+ndebug),stat=i_stat)
  call memocc(i_stat,orbs%ikptproc,'orbs%ikptproc',subname)

end subroutine orbitals_descriptors
!!***


!!****f* BigDFT/input_occup
!! FUNCTION
!!    Fill the arrays occup and spinsgn
!!    if iunit /=0 this means that the file 'occup.dat' does exist and it opens
!! SOURCE
!!
subroutine input_occup(iproc,iunit,nelec,norb,norbu,norbd,nspin,mpol,occup,spinsgn)
  use module_base
  implicit none
! Arguments
  integer, intent(in) :: nelec,nspin,mpol,iproc,norb,norbu,norbd,iunit
  real(gp), dimension(norb), intent(out) :: occup,spinsgn
! Local variables
  integer :: iorb,nt,ne,it,ierror,iorb1
  real(gp) :: rocc,rup,rdown
  character(len=100) :: line


  do iorb=1,norb
     spinsgn(iorb)=1.0_gp
  end do
  if (nspin/=1) then
     do iorb=1,norbu
        spinsgn(iorb)=1.0_gp
     end do
     do iorb=norbu+1,norb
        spinsgn(iorb)=-1.0_gp
     end do
  end if
! write(*,'(1x,a,5i4,30f6.2)')'Spins: ',norb,norbu,norbd,norbup,norbdp,(spinsgn(iorb),iorb=1,norb)

! First fill the occupation numbers by default
  nt=0
  if (nspin==1) then
     ne=(nelec+1)/2
     do iorb=1,ne
        it=min(2,nelec-nt)
        occup(iorb)=real(it,gp)
        nt=nt+it
     enddo
     do iorb=ne+1,norb
        occup(iorb)=0._gp
     end do
  else
     do iorb=1,norb
        it=min(1,nelec-nt)
        occup(iorb)=real(it,gp)
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
              !if (iproc==0) then
                 write(*,'(1x,a,i0,a)') 'ERROR in line ',nt+1,' of the file "occup.dat"'
                 write(*,'(10x,a,i0,a)') 'The orbital index ',iorb,' is incorrect'
              !end if
              stop
           elseif (rocc<0._gp .or. rocc>2._gp) then
              !if (iproc==0) then
                 write(*,'(1x,a,i0,a)') 'ERROR in line ',nt+1,' of the file "occup.dat"'
                 write(*,'(10x,a,f5.2,a)') 'The occupation number ',rocc,' is not between 0. and 2.'
              !end if
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
     if (abs(rocc-real(nelec,gp))>1.e-6_gp) then
        !if (iproc==0) then
           write(*,'(1x,a,f13.6,a,i0)') 'From the file "occup.dat", the total number of electrons ',rocc,&
                          ' is not equal to ',nelec
        !end if
        stop
     end if
     if (nspin/=1) then
!!$        !Check if the polarisation is respected (mpol)
!!$        rup=sum(occup(1:norbu))
!!$        rdown=sum(occup(norbu+1:norb))
!!$        if (abs(rup-rdown-real(norbu-norbd,gp))>1.e-6_gp) then
!!$           if (iproc==0) then
!!$              write(*,'(1x,a,f13.6,a,i0)') 'From the file "occup.dat", the polarization ',rup-rdown,&
!!$                             ' is not equal to ',norbu-norbd
!!$           end if
!!$           stop
!!$        end if
        !Fill spinsgn
        do iorb=1,norbu
           spinsgn(iorb)=1.0_gp
        end do
        do iorb=norbu+1,norb
           spinsgn(iorb)=-1.0_gp
        end do
      end if
  end if
  if (iproc==0) then 
     write(*,'(1x,a,t28,i8)') 'Total Number of Orbitals',norb
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
!!***
