!> @file
!!  Routines related to system properties
!! @author
!!    Copyright (C) 2010-2011 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 


!>  Calculate the important objects related to the physical properties of the system
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
  character(len=*), parameter :: subname='system_properties'
  integer :: iunit,norb,norbu,norbd,nspinor,jpst,norbme,norbyou,jproc,ikpts
  integer :: norbuempty,norbdempty

  call read_system_variables('input.occup',iproc,in,atoms,radii_cf,nelec,&
       norb,norbu,norbd,norbuempty,norbdempty,iunit)

  if(in%nspin==4) then
     nspinor=4
  else
     nspinor=1
  end if

  call orbitals_descriptors(iproc, nproc,norb,norbu,norbd,in%nspin,nspinor, &
       & in%nkpt,in%kpt,in%wkpt,orbs)

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
     call input_occup(iproc,iunit,nelec,norb,norbu,norbuempty,norbdempty,in%nspin,&
          orbs%occup(1+(ikpts-1)*orbs%norb),orbs%spinsgn(1+(ikpts-1)*orbs%norb))
  end do
END SUBROUTINE system_properties


!>  Check for the need of a core density and fill the rhocore array which
!!  should be passed at the rhocore pointer
subroutine calculate_rhocore(iproc,at,d,rxyz,hxh,hyh,hzh,i3s,i3xcsh,n3d,n3p,rhocore)
  use module_base
  use module_types
  implicit none
  integer, intent(in) :: iproc,i3s,n3d,i3xcsh,n3p
  real(gp), intent(in) :: hxh,hyh,hzh
  type(atoms_data), intent(in) :: at
  type(grid_dimensions), intent(in) :: d
  real(gp), dimension(3,at%nat), intent(in) :: rxyz
  real(wp), dimension(:), pointer :: rhocore
  !local variables
  character(len=*), parameter :: subname='calculate_rhocore'
  integer :: ityp,iat,i_stat,j3,i1,i2,ind,ierr
  real(wp) :: tt
  real(gp) :: rx,ry,rz,rloc,cutoff
  

  !check for the need of a nonlinear core correction
!!$  donlcc=.false.
!!$  chk_nlcc: do ityp=1,at%ntypes
!!$     filename = 'nlcc.'//at%atomnames(ityp)
!!$
!!$     inquire(file=filename,exist=exists)
!!$     if (exists) then
!!$        donlcc=.true.
!!$        exit chk_nlcc
!!$     end if
!!$  end do chk_nlcc

  if (at%donlcc) then
     !allocate pointer rhocore
     allocate(rhocore(d%n1i*d%n2i*n3d+ndebug),stat=i_stat)
     call memocc(i_stat,rhocore,'rhocore',subname)
     !initalise it 
     call razero(d%n1i*d%n2i*n3d,rhocore)
     !perform the loop on any of the atoms which have this feature
     do iat=1,at%nat
        ityp=at%iatype(iat)
!!$        filename = 'nlcc.'//at%atomnames(ityp)
!!$        inquire(file=filename,exist=exists)
!!$        if (exists) then
        if (at%nlcc_ngv(ityp)/=UNINITIALIZED(1) .or. at%nlcc_ngc(ityp)/=UNINITIALIZED(1) ) then
           if (iproc == 0) write(*,'(1x,a)',advance='no')&
                'NLCC: calculate core density for atom: '//&
                trim(at%atomnames(ityp))//';'
           rx=rxyz(1,iat) 
           ry=rxyz(2,iat)
           rz=rxyz(3,iat)

           rloc=at%psppar(0,0,ityp)
           cutoff=10.d0*rloc

           call calc_rhocore_iat(iproc,at,ityp,rx,ry,rz,cutoff,hxh,hyh,hzh,&
                d%n1,d%n2,d%n3,d%n1i,d%n2i,d%n3i,&
                i3s,n3d,rhocore)

           if (iproc == 0) write(*,'(1x,a)')'done.'
        end if
     end do

     !calculate total core charge in the grid
     !In general this should be really bad
     tt=0.0_wp
     do j3=1,n3p
        do i2=1,d%n2i
           do i1=1,d%n1i
              ind=i1+(i2-1)*d%n1i+(j3+i3xcsh-1)*d%n1i*d%n2i
              tt=tt+rhocore(ind)
           enddo
        enddo
     enddo
     call mpiallred(tt,1,MPI_SUM,MPI_COMM_WORLD,ierr)
     tt=tt*hxh*hyh*hzh
     if (iproc == 0) write(*,'(1x,a,f15.7)') &
       'Total core charge on the grid (generally bad, overestimated approx.): ',tt

  else
     !No NLCC needed, nullify the pointer 
     nullify(rhocore)
  end if

END SUBROUTINE calculate_rhocore


!>   Assign some of the physical system variables
!!   Performs also some cross-checks with other variables
!!   The pointer in atoms structure have to be associated or nullified.
subroutine read_system_variables(fileocc,iproc,in,atoms,radii_cf,&
     nelec,norb,norbu,norbd,norbuempty,norbdempty,iunit)
  use module_base
  use module_types
  use module_xc
  use ab6_symmetry
  implicit none
  character (len=*), intent(in) :: fileocc
  type(input_variables), intent(in) :: in
  integer, intent(in) :: iproc
  type(atoms_data), intent(inout) :: atoms
  integer, intent(out) :: nelec,norb,norbu,norbd,iunit,norbuempty,norbdempty
  real(gp), dimension(atoms%ntypes,3), intent(out) :: radii_cf
  !local variables
  character(len=*), parameter :: subname='read_system_variables'
  integer, parameter :: nelecmax=32,nmax=6,lmax=4,noccmax=2
  logical :: exists
  character(len=2) :: symbol
  character(len=24) :: message
  character(len=27) :: filename
  character(len=50) :: format
  character(len=100) :: line
  integer :: i,j,k,l,iat,nlterms,nprl,nn,nt,ntu,ntd,ityp,ierror,i_stat,ixcpsp,ispinsum,mxpl,ig
  integer :: ispol,mxchg,ichg,ichgsum,nsccode,ierror1,norbe,norbat,nspinor,nspin,nlcc_dim,ngv,ngc
  real(gp) :: rcov,rprb,ehomo,radfine,minrad,maxrad
  real(gp), dimension(3,3) :: hij
  real(gp), dimension(2,2,3) :: offdiagarr
  !integer, dimension(nmax,0:lmax-1) :: neleconf
  real(kind=8), dimension(nmax,0:lmax-1) :: neleconf
  integer, dimension(lmax) :: nl
  real(gp), dimension(0:4) :: fake_nlcc
  real(gp), dimension(noccmax,lmax) :: occup
  character(len=500) :: name_xc1, name_xc2


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
  allocate(atoms%iasctype(atoms%nat+ndebug),stat=i_stat)
  call memocc(i_stat,atoms%iasctype,'atoms%iasctype',subname)
  allocate(atoms%aocc(nelecmax,atoms%nat+ndebug),stat=i_stat)
  call memocc(i_stat,atoms%aocc,'atoms%aocc',subname)

  !parameters for NLCC
  allocate(atoms%nlcc_ngv(atoms%ntypes+ndebug),stat=i_stat)
  call memocc(i_stat,atoms%nlcc_ngv,'atoms%nlcc_ngv',subname)
  allocate(atoms%nlcc_ngc(atoms%ntypes+ndebug),stat=i_stat)
  call memocc(i_stat,atoms%nlcc_ngc,'atoms%nlcc_ngc',subname)

  if (iproc == 0) then
     write(*,'(1x,a)')&
          ' Atom    N.Electr.  PSP Code  Radii: Coarse     Fine  CoarsePSP    Calculated   File'
  end if

  !logical variable for non-linear core correction
  atoms%donlcc=.false.
  nlcc_dim=0
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
     if (ixcpsp < 0) then
        call xc_get_name(name_xc1, ixcpsp, XC_MIXED)
     else
        call xc_get_name(name_xc1, ixcpsp, XC_ABINIT)
     end if
     if (in%ixc < 0) then
        call xc_get_name(name_xc2, in%ixc, XC_MIXED)
     else
        call xc_get_name(name_xc2, in%ixc, XC_ABINIT)
     end if
     if (trim(name_xc1) /= trim(name_xc2) .and. iproc==0) then
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
     !and consider the ground state electronic configuration
     call eleconf(atoms%nzatom(ityp),atoms%nelpsp(ityp),symbol,rcov,rprb,ehomo,&
          neleconf,nsccode,mxpl,mxchg,atoms%amu(ityp))

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

     !correct the coarse radius for projectors
     !it is always multiplied by frmult
     !NOTE this radius is chosen such as to make the projector be defined always on the same sphere
     !     of the atom. This is clearly too much since such sphere is built to the exp decay of the wavefunction
     !     and not for the gaussian decaying of the pseudopotential projector
     !     add a proper varialbe in input.perf
     radii_cf(ityp,3)=max(min(in%crmult*radii_cf(ityp,1),in%projrad*maxrad)/in%frmult,radii_cf(ityp,2))

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

     call atomic_occupation_numbers(fileocc,ityp,in%nspin,atoms,nmax,lmax,nelecmax,&
          neleconf,nsccode,mxpl,mxchg)

     !inquire for non-linear core correction and count the components
     filename ='nlcc.'//atoms%atomnames(ityp)

     inquire(file=filename,exist=exists)
     if (exists) then
        atoms%donlcc=.true. !just one atom would be enough
        !associate the number of gaussians
        open(unit=79,file=filename,status='unknown')
        read(79,*)ngv
        if (ngv==0) then
           atoms%nlcc_ngv(ityp)=UNINITIALIZED(1)
        else
           atoms%nlcc_ngv(ityp)=ngv
        end if
        nlcc_dim=nlcc_dim+(ngv*(ngv+1)/2)
        do ig=1,(ngv*(ngv+1))/2
           read(79,*) (fake_nlcc(j),j=0,4)!jump the suitable lines (the file is organised with one element per line)
        end do
        read(79,*)ngc
        if (ngc==0) then
           atoms%nlcc_ngc(ityp)=UNINITIALIZED(1)
        else
           atoms%nlcc_ngc(ityp)=ngc
        end if
        nlcc_dim=nlcc_dim+(ngc*(ngc+1))/2
        !better to read values in a fake array
        do ig=1,(ngc*(ngc+1))/2
           read(79,*) (fake_nlcc(j),j=0,4)!jump the suitable lines (the file is organised with one element per line)
        end do
        !no need to go further for the moment
        close(unit=79)
     else
        atoms%nlcc_ngv(ityp)=UNINITIALIZED(1)
        atoms%nlcc_ngc(ityp)=UNINITIALIZED(1)
     end if
  enddo

  !process the nlcc parameters if present 
  !(allocation is performed also with zero size)
  allocate(atoms%nlccpar(0:4,max(nlcc_dim,1)+ndebug),stat=i_stat)
  call memocc(i_stat,atoms%nlccpar,'atoms%nlccpar',subname)
  !start again the file inspection to fill nlcc parameters
  if (atoms%donlcc) then
     nlcc_dim=0
     fill_nlcc: do ityp=1,atoms%ntypes
        filename = 'nlcc.'//atoms%atomnames(ityp)
        inquire(file=filename,exist=exists)
        if (exists) then
           !read the values of the gaussian for valence and core densities
           open(unit=79,file=filename,status='unknown')
           read(79,*)ngv
           do ig=1,(ngv*(ngv+1))/2
              nlcc_dim=nlcc_dim+1
              read(79,*)(atoms%nlccpar(j,nlcc_dim),j=0,4)!rhovxp(ig),(rhovc(ig,j),j=1,4)
           end do
           read(79,*)ngc
           do ig=1,(ngc*(ngc+1))/2
              nlcc_dim=nlcc_dim+1
              read(79,*)(atoms%nlccpar(j,nlcc_dim),j=0,4)!rhocxp(ig),(rhocc(ig,j),j=1,4)
           end do
           close(unit=79)
        end if
     end do fill_nlcc
  end if

  !print *,'iatsctype',atOMS%iasctype(:)

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
     nsccode=atoms%iasctype(iat)

     !print *,'nsccode,iat2',nsccode
     !this part should be removed one the occupation number has been passed
     !call charge_and_spol(atoms%natpol(iat),ichg,ispol)
     !if (ichg /=0) then
     !   call eleconf(atoms%nzatom(ityp),atoms%nelpsp(ityp),symbol,rcov,rprb,ehomo,&
     !        neleconf,atoms%iasctype(ityp),mxpl,mxchg,atoms%amu(ityp))
     !   call correct_semicore(6,3,ichg,neleconf,nsccode)
     !end if
     !end of part to be removed
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
     if (mod(norb+in%mpol,2) /=0) then
        write(*,*)'ERROR: the input polarization should have the same parity of the number of electrons'
        stop
     end if
     norbu=min((norb+in%mpol)/2,norb)
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
        if (iproc==0 .and. in%norbsempty == 0) &
             write(*,'(1x,a)')&
             'WARNING: Found no input polarisation, add it for a correct input guess'
        !stop
     end if

  end if

  !initialise the values for the empty orbitals
  norbuempty=0
  norbdempty=0

  ! Test if the file 'input.occ exists
  inquire(file='input.occ',exist=exists)
  iunit=0
  if (exists) then
     iunit=25
     open(unit=iunit,file='input.occ',form='formatted',action='read',status='old')
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
             'ERROR: reading the number of orbitals in the file "input.occ"'
        stop
     end if
     !Check
     if (in%nspin==1) then
        if (nt<norb) then
           !if (iproc==0) 
               write(*,'(1x,a,i0,a,i0)') &
                'ERROR: In the file "input.occ", the number of orbitals norb=',nt,&
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
                'ERROR: In the file "input.occ", the number of orbitals norb=',nt,&
                ' should be greater or equal than nelec=',norb
           stop
        else
           norb=nt
        end if
        if (ntu<norbu) then
           !if (iproc==0) 
                write(*,'(1x,a,i0,a,i0)') &
                'ERROR: In the file "input.occ", the number of orbitals up norbu=',ntu,&
                ' should be greater or equal than min((nelec+mpol)/2,nelec)=',norbu
           stop
        else
           norbu=ntu
        end if
        if (ntd<norbd) then
           !if (iproc==0) 
                  write(*,'(1x,a,i0,a,i0)') &
                'ERROR: In the file "input.occ", the number of orbitals down norbd=',ntd,&
                ' should be greater or equal than min((nelec-mpol/2),0)=',norbd
           stop
        else
           norbd=ntd
        end if
     end if
  else if (in%norbsempty > 0) then
     !total number of orbitals
     norbe=0
     if(in%nspin==4) then
        nspin=2
        nspinor=4
     else
        nspin=in%nspin
        nspinor=1
     end if

     do iat=1,atoms%nat
        ityp=atoms%iatype(iat)
        call count_atomic_shells(lmax,noccmax,nelecmax,nspin,nspinor,atoms%aocc(1,iat),occup,nl)
        norbat=(nl(1)+3*nl(2)+5*nl(3)+7*nl(4))
        norbe=norbe+norbat
     end do

     !value of empty orbitals up and down, needed to fill occupation numbers
     norbuempty=min(in%norbsempty,norbe-norbu)
     norbdempty=min(in%norbsempty,norbe-norbd)

     if (in%nspin == 4 .or. in%nspin==1) then
        norb=norb+norbuempty
        norbu=norbu+norbuempty
     else if (in%nspin ==2) then
        norbu=norbu+norbuempty
        norbd=norbd+norbdempty
        norb=norbu+norbd
     end if
  end if

  ! We modify the symmetry object with respect to the spin.
  if (atoms%symObj >= 0) then
     if (in%nspin == 2) then
        call ab6_symmetry_set_collinear_spin(atoms%symObj, atoms%nat, &
             & atoms%natpol, ierror)
!!$     else if (in%nspin == 4) then
!!$        call ab6_symmetry_set_spin(atoms%symObj, atoms%nat, &
!!$             & atoms%natpol, ierror)
     end if
  end if

!!!  tt=dble(norb)/dble(nproc)
!!!  norbp=int((1.d0-eps_mach*tt) + tt)
!!!  !if (iproc.eq.0) write(*,'(1x,a,1x,i0)') 'norbp=',norbp

END SUBROUTINE read_system_variables

!>find the correct position of the nlcc parameters
subroutine nlcc_start_position(ityp,atoms,ngv,ngc,islcc)
  use module_base
  use module_types
  implicit none
  integer, intent(in) :: ityp
  type(atoms_data), intent(in) :: atoms
  integer, intent(out) :: ngv,ngc,islcc
  !local variables
  integer :: ilcc,jtyp

  ilcc=0
  do jtyp=1,ityp-1
     ngv=atoms%nlcc_ngv(jtyp)
     if (ngv /= UNINITIALIZED(ngv)) ilcc=ilcc+(ngv*(ngv+1)/2)
     ngc=atoms%nlcc_ngc(jtyp)
     if (ngc /= UNINITIALIZED(ngc)) ilcc=ilcc+(ngc*(ngc+1))/2
  end do
  islcc=ilcc

  ngv=atoms%nlcc_ngv(ityp)
  if (ngv==UNINITIALIZED(1)) ngv=0
  ngc=atoms%nlcc_ngc(ityp)
  if (ngc==UNINITIALIZED(1)) ngc=0
END SUBROUTINE nlcc_start_position

!>   Fix all the atomic occupation numbers of the atoms which has the same type
!!   look also at the input polarisation and spin
!!   look at the file of the input occupation numbers and, if exists, modify the 
!!   occupations accordingly
subroutine atomic_occupation_numbers(filename,ityp,nspin,at,nmax,lmax,nelecmax,neleconf,nsccode,mxpl,mxchg)
  use module_base
  use module_types
  implicit none
  character(len=*), intent(in) :: filename
  integer, intent(in) :: ityp,mxpl,mxchg,nspin,nmax,lmax,nelecmax,nsccode
  type(atoms_data), intent(inout) :: at
  !integer, dimension(nmax,lmax), intent(in) :: neleconf
  real(gp), dimension(nmax,lmax), intent(in) :: neleconf
  !local variables
  integer, parameter :: noccmax=2
  character(len=100) :: string
  logical :: exists,found
  integer :: iat,ichg,ispol,nsp,nspinor,ierror,jat,l,iocc,icoll,noncoll,ispin,nl,inl,m
  real(gp) :: elec
  real(gp), dimension(nmax,lmax) :: eleconf

  !control the spin
  select case(nspin)
     case(1)
        nsp=1
        nspinor=1
        noncoll=1
     case(2)
        nsp=2
        nspinor=1
        noncoll=1
     case(4)
        nsp=1
        nspinor=4
        noncoll=2
     case default
        write(*,*)' ERROR: nspin not valid:',nspin
        stop
  end select

  inquire(file=filename,exist=exists)

  !search the corresponding atom
  if (exists) then
     open(unit=91,file=filename,status='old',iostat=ierror)
     !Check the open statement
     if (ierror /= 0) then
        write(*,*)'Failed to open the existing  file: '//filename
        stop
     end if
  end if

  !here we must check of the input guess polarisation
  !control if the values are compatible with the atom configuration
  !do this for all atoms belonging to a given type
  !control the maximum polarisation allowed: consider only non-closed shells   
  do iat=1,at%nat
     !control the atomic input guess occupation number
     !if you want no semicore input guess electrons, uncomment the following line
     !at%iasctype(iat)=0
     if (at%iatype(iat) == ityp) then
        !see whether the input.occup file contains the given atom
        found=.false.
        if (exists) then
           rewind(unit=91)
           parse_inocc: do
              read(91,'(a100)',iostat=ierror)string
              if (ierror /= 0) exit parse_inocc !file ends
              read(string,*,iostat=ierror)jat
              if (ierror /=0) stop 'Error reading line'
              if (jat==iat ) then
                 found=.true.
                 exit parse_inocc
              end if
           end do parse_inocc
        end if
        call charge_and_spol(at%natpol(iat),ichg,ispol)
        if (found) then
           call read_eleconf(string,nsp,nspinor,noccmax,nelecmax,lmax,&
                at%aocc(1,iat),at%iasctype(iat))
        else
           at%iasctype(iat)=nsccode
           if (abs(ispol) > mxpl+abs(ichg)) then
              !if (iproc ==0) 
              write(*,'(1x,a,i0,a,a,2(a,i0))')&
                   'ERROR: Input polarisation of atom No.',iat,&
                   ' (',trim(at%atomnames(ityp)),') must be <=',mxpl,&
                   ', while found ',ispol
              stop 
           end if
           if (abs(ichg) > mxchg) then
              !if (iproc ==0) 
              write(*,'(1x,a,i0,a,a,2(a,i0))')&
                   'ERROR: Input charge of atom No.',iat,&
                   ' (',trim(at%atomnames(ityp)),') must be <=',mxchg,&
                   ', while found ',ichg
              stop
           end if
           !correct the electronic configuration in case there is a charge
           !if (ichg /=0) then
           call correct_semicore(nmax,lmax-1,ichg,&
                neleconf,eleconf,at%iasctype(iat))
           !end if

           call at_occnums(ispol,nsp,nspinor,nmax,lmax,nelecmax,&
                eleconf,at%aocc(1,iat))
        end if

        !check the total number of electrons
        elec=0.0_gp
        iocc=0
        do l=1,lmax
           iocc=iocc+1
           nl=nint(at%aocc(iocc,iat))
           do inl=1,nl
              do ispin=1,nsp
                 do m=1,2*l-1
                    do icoll=1,noncoll !non-trivial only for nspinor=4
                       iocc=iocc+1
                       elec=elec+at%aocc(iocc,iat)
                    end do
                 end do
              end do
           end do
        end do
        if (nint(elec) /= at%nelpsp(ityp) - ichg) then
           write(*,*)'ERROR: the total atomic charge ',elec,&
                ' is different from the PSP charge ',at%nelpsp(ityp),&
                ' plus the charge ',-ichg
           stop
        end if
     end if
  end do

  if (exists) close(unit=91)

END SUBROUTINE atomic_occupation_numbers


!> Define the descriptors of the orbitals from a given norb
!! It uses the cubic strategy for partitioning the orbitals
subroutine orbitals_descriptors(iproc,nproc,norb,norbu,norbd,nspin,nspinor,nkpt,kpt,wkpt,orbs)
  use module_base
  use module_types
  implicit none
  integer, intent(in) :: iproc,nproc,norb,norbu,norbd,nkpt,nspin
  integer, intent(in) :: nspinor
  type(orbitals_data), intent(out) :: orbs
  real(gp), dimension(nkpt), intent(in) :: wkpt
  real(gp), dimension(3,nkpt), intent(in) :: kpt
  !local variables
  character(len=*), parameter :: subname='orbitals_descriptors'
  integer :: iorb,jproc,norb_tot,ikpt,i_stat,jorb,ierr,i_all
  logical, dimension(:), allocatable :: GPU_for_orbs
  integer, dimension(:), allocatable :: mykpts
  integer, dimension(:,:), allocatable :: norb_par !(with k-pts)

  allocate(orbs%norb_par(0:nproc-1+ndebug),stat=i_stat)
  call memocc(i_stat,orbs%norb_par,'orbs%norb_par',subname)

  !assign the value of the k-points
  orbs%nkpts=nkpt
  !allocate vectors related to k-points
  allocate(orbs%kpts(3,orbs%nkpts+ndebug),stat=i_stat)
  call memocc(i_stat,orbs%kpts,'orbs%kpts',subname)
  allocate(orbs%kwgts(orbs%nkpts+ndebug),stat=i_stat)
  call memocc(i_stat,orbs%kwgts,'orbs%kwgts',subname)
  orbs%kpts(:,1:nkpt) = kpt(:,:)
  orbs%kwgts(1:nkpt) = wkpt(:)

  ! Change the wavefunctions to complex if k-points are used (except gamma).
  orbs%nspinor=nspinor
  if (nspinor == 1) then
     if (maxval(abs(orbs%kpts)) > 0._gp) orbs%nspinor=2
     !nspinor=2 !fake, used for testing with gamma
  end if
  orbs%nspin = nspin

  !initialise the array
  do jproc=0,nproc-1
     orbs%norb_par(jproc)=0 !size 0 nproc-1
  end do

  !create an array which indicate which processor has a GPU associated 
  !from the viewpoint of the BLAS routines (deprecated, not used anymore)
  if (.not. GPUshare) then
     allocate(GPU_for_orbs(0:nproc-1+ndebug),stat=i_stat)
     call memocc(i_stat,GPU_for_orbs,'GPU_for_orbs',subname)
     
     if (nproc > 1) then
        call MPI_ALLGATHER(GPUconv,1,MPI_LOGICAL,GPU_for_orbs(0),1,MPI_LOGICAL,&
             MPI_COMM_WORLD,ierr)
     else
        GPU_for_orbs(0)=GPUconv
     end if
     
     i_all=-product(shape(GPU_for_orbs))*kind(GPU_for_orbs)
     deallocate(GPU_for_orbs,stat=i_stat)
     call memocc(i_stat,i_all,'GPU_for_orbs',subname)
  end if

  allocate(norb_par(0:nproc-1,orbs%nkpts+ndebug),stat=i_stat)
  call memocc(i_stat,norb_par,'norb_par',subname)

  !old system for calculating k-point repartition
!!$  call parallel_repartition_with_kpoints(nproc,orbs%nkpts,norb,orbs%norb_par)
!!$
!!$  !check the distribution
!!$  norb_tot=0
!!$  do jproc=0,iproc-1
!!$     norb_tot=norb_tot+orbs%norb_par(jproc)
!!$  end do
!!$  !reference orbital for process
!!$  orbs%isorb=norb_tot
!!$  do jproc=iproc,nproc-1
!!$     norb_tot=norb_tot+orbs%norb_par(jproc)
!!$  end do
!!$
!!$  if(norb_tot /= norb*orbs%nkpts) then
!!$     write(*,*)'ERROR: partition of orbitals incorrect, report bug.'
!!$     write(*,*)orbs%norb_par(:),norb*orbs%nkpts
!!$     stop
!!$  end if
!!$
!!$  !calculate the k-points related quantities
!!$  allocate(mykpts(orbs%nkpts+ndebug),stat=i_stat)
!!$  call memocc(i_stat,mykpts,'mykpts',subname)
!!$
!!$  call parallel_repartition_per_kpoints(iproc,nproc,orbs%nkpts,norb,orbs%norb_par,&
!!$       orbs%nkptsp,mykpts,norb_par)
!!$  if (orbs%norb_par(iproc) >0) then
!!$     orbs%iskpts=mykpts(1)-1
!!$  else
!!$     orbs%iskpts=0
!!$  end if
!!$  i_all=-product(shape(mykpts))*kind(mykpts)
!!$  deallocate(mykpts,stat=i_stat)
!!$  call memocc(i_stat,i_all,'mykpts',subname)

  !new system for k-point repartition
  call kpts_to_procs_via_obj(nproc,orbs%nkpts,norb,norb_par)
  !assign the values for norb_par and check the distribution
  norb_tot=0
  do jproc=0,nproc-1
     if (jproc==iproc) orbs%isorb=norb_tot
     do ikpt=1,orbs%nkpts
        orbs%norb_par(jproc)=orbs%norb_par(jproc)+norb_par(jproc,ikpt)
     end do
     norb_tot=norb_tot+orbs%norb_par(jproc)
  end do

  if(norb_tot /= norb*orbs%nkpts) then
     write(*,*)'ERROR: partition of orbitals incorrect, report bug.'
     write(*,*)orbs%norb_par(:),norb*orbs%nkpts
     stop
  end if


  !allocate(orbs%ikptsp(orbs%nkptsp+ndebug),stat=i_stat)
  !call memocc(i_stat,orbs%ikptsp,'orbs%ikptsp',subname)
  !orbs%ikptsp(1:orbs%nkptsp)=mykpts(1:orbs%nkptsp)

  !this array will be reconstructed in the orbitals_communicators routine
  i_all=-product(shape(norb_par))*kind(norb_par)
  deallocate(norb_par,stat=i_stat)
  call memocc(i_stat,i_all,'norb_par',subname)

  !assign the values of the orbitals data
  orbs%norb=norb
  orbs%norbp=orbs%norb_par(iproc)
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
  do ikpt=1,orbs%nkpts
     do iorb=1,orbs%norbu
        orbs%spinsgn(iorb+(ikpt-1)*orbs%norb)=1.0_gp
     end do
     do iorb=1,orbs%norbd
        orbs%spinsgn(iorb+orbs%norbu+(ikpt-1)*orbs%norb)=-1.0_gp
     end do
  end do

  !put a default value for the fermi energy
  orbs%efermi = UNINITIALIZED(orbs%efermi)
  !and also for the gap
  orbs%HLgap = UNINITIALIZED(orbs%HLgap)

  !allocate the array which assign the k-point to processor in transposed version
  allocate(orbs%ikptproc(orbs%nkpts+ndebug),stat=i_stat)
  call memocc(i_stat,orbs%ikptproc,'orbs%ikptproc',subname)

END SUBROUTINE orbitals_descriptors


!>    Fill the arrays occup and spinsgn
!!    if iunit /=0 this means that the file 'input.occ' does exist and it opens
subroutine input_occup(iproc,iunit,nelec,norb,norbu,norbuempty,norbdempty,nspin,occup,spinsgn)
  use module_base
  implicit none
  ! Arguments
  integer, intent(in) :: nelec,nspin,iproc,norb,norbu,iunit,norbuempty,norbdempty
  real(gp), dimension(norb), intent(out) :: occup,spinsgn
  ! Local variables
  integer :: iorb,nt,ne,it,ierror,iorb1,i
  real(gp) :: rocc
  character(len=20) :: string
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
     if (norbuempty+norbdempty == 0) then
        if (norb > nelec) then
           do iorb=1,min(norbu,norb/2+1)
              it=min(1,nelec-nt)
              occup(iorb)=real(it,gp)
              nt=nt+it
           enddo
           do iorb=min(norbu,norb/2+1)+1,norbu
              occup(iorb)=0.0_gp
           end do
           do iorb=norbu+1,norbu+min(norb-norbu,norb/2+1)
              it=min(1,nelec-nt)
              occup(iorb)=real(it,gp)
              nt=nt+it
           enddo
           do iorb=norbu+min(norb-norbu,norb/2+1)+1,norb
              occup(iorb)=0.0_gp
           end do
        else
           do iorb=1,norb
              occup(iorb)=1.0_gp
           end do
        end if
     else
        do iorb=1,norbu-norbuempty
           occup(iorb)=1.0_gp
        end do
        do iorb=norbu-norbuempty+1,norbu
           occup(iorb)=0.0_gp
        end do
        do iorb=1,norb-norbu-norbdempty
           occup(norbu+iorb)=1.0_gp
        end do
        do iorb=norb-norbu-norbdempty+1,norb-norbu
           occup(norbu+iorb)=0.0_gp
        end do
     end if
  end if
  ! Then read the file "input.occ" if does exist
  if (iunit /= 0) then
     nt=0
     do
        read(unit=iunit,fmt='(a100)',iostat=ierror) line
        if (ierror /= 0) then
           exit
        end if
        !Transform the line in case there are slashes (to ease the parsing)
        do i=1,len(line)
           if (line(i:i) == '/') then
              line(i:i) = ':'
           end if
        end do
        read(line,*,iostat=ierror) iorb,string
        call read_fraction_string(string,rocc,ierror) 
        if (ierror /= 0) then
           exit
        end if

        if (ierror/=0) then
           exit
        else
           nt=nt+1
           if (iorb<0 .or. iorb>norb) then
              !if (iproc==0) then
              write(*,'(1x,a,i0,a)') 'ERROR in line ',nt+1,' of the file "input.occ"'
              write(*,'(10x,a,i0,a)') 'The orbital index ',iorb,' is incorrect'
              !end if
              stop
           elseif (rocc<0._gp .or. rocc>2._gp) then
              !if (iproc==0) then
              write(*,'(1x,a,i0,a)') 'ERROR in line ',nt+1,' of the file "input.occ"'
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
             'The occupation numbers are read from the file "input.occ" (',nt,' lines read)'
     end if
     close(unit=iunit)

     if (nspin/=1) then
!!!        !Check if the polarisation is respected (mpol)
!!!        rup=sum(occup(1:norbu))
!!!        rdown=sum(occup(norbu+1:norb))
!!!        if (abs(rup-rdown-real(norbu-norbd,gp))>1.e-6_gp) then
!!!           if (iproc==0) then
!!!              write(*,'(1x,a,f13.6,a,i0)') 'From the file "input.occ", the polarization ',rup-rdown,&
!!!                             ' is not equal to ',norbu-norbd
!!!           end if
!!!           stop
!!!        end if
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

  !Check if sum(occup)=nelec
  rocc=sum(occup)
  if (abs(rocc-real(nelec,gp))>1.e-6_gp) then
     !if (iproc==0) then
     write(*,'(1x,a,f13.6,a,i0)') 'ERROR in determining the occupation numbers: the total number of electrons ',rocc,&
          ' is not equal to ',nelec
     !end if
     stop
  end if


END SUBROUTINE input_occup

!> Routine which assign to each processor the repartition of nobj*nkpts objects
subroutine kpts_to_procs_via_obj(nproc,nkpts,nobj,nobj_par)
  use module_base
  implicit none
  integer, intent(in) :: nproc,nkpts,nobj
  integer, dimension(0:nproc-1,nkpts), intent(out) :: nobj_par
  !local varaibles
  logical :: intrep
  integer :: jproc,ikpt,iobj,nobjp_max_kpt,nprocs_with_floor,jobj,nobjp,nprocs_with_ceiling
  integer :: jkpt,nproc_per_kpt,nproc_left,kproc,nkpt_per_proc,nkpts_left
  real(gp) :: robjp,rounding_ratio

  !decide the naive number of objects which should go to each processor.
  robjp=real(nobj,gp)*real(nkpts,gp)/real(nproc,gp)

  !maximum number of objects which has to go to each processor per k-point
  nobjp_max_kpt=ceiling(modulo(robjp-epsilon(1.0_gp),real(nobj,gp)))
  !print *,'eccoqui',nobjp_max_kpt,robjp
  !see the conditions for the integer repartition of k-points
  if (nobjp_max_kpt == nobj .or. (nobjp_max_kpt==1 .and. robjp < 1.0_gp)) then
     intrep=.true.
     rounding_ratio=0.0_gp
     nprocs_with_floor=0
  else
     intrep=.false.
     !the repartition is not obvious, some processors take nobj_max_kpt objects, others take the previous integer.
     !to understand how many, we round the percentage of processors which is given by
     rounding_ratio=(robjp-real(floor(robjp)))
     !then this is the number of processors which will take the floor
     nprocs_with_floor=ceiling((1.0_gp-rounding_ratio)*real(nproc,gp))!nproc-(nobj*nkpts-floor(robjp)*nproc)
     !print *,'rounding_ratio,nprocs_with_floor',rounding_ratio,nprocs_with_floor
     if (nprocs_with_floor > nproc) stop 'ERROR: should not happen'
     !if (nprocs_with_floor == nproc) nprocs_with_floor=nproc-1
  end if

  !start separating the objects for the repartition which is suggested by rounding_ratio and nprocs_with_floor
  nobj_par(0:nproc-1,1:nkpts)=0
  !integer repartition
  if (intrep) then 
     !strategy for the repartition
     if (nproc >= nkpts) then
        !decide in how many processors a single k-point can be partitioned
        nproc_per_kpt=max((nproc-1),1)/nkpts !this is the minimum
        !count how many processors are left that way
        !distribute the k-point among these
        nproc_left=nproc-nproc_per_kpt*nkpts
        ikpt=0
        jproc=0
        !print *,'qui',nproc_left,nproc_per_kpt
        do kproc=0,nproc_left-1
           ikpt=ikpt+1
           if (ikpt > nkpts) stop 'ERROR: also this should not happen3'
           do iobj=0,nobj-1
              nobj_par(jproc+modulo(iobj,nproc_per_kpt+1),ikpt)=nobj_par(jproc+modulo(iobj,nproc_per_kpt+1),ikpt)+1
           end do
           jproc=jproc+nproc_per_kpt+1
        end do
        !print *,'ciao'
        if ((nproc_per_kpt+1)*nproc_left < nproc) then
           do jproc=(nproc_per_kpt+1)*nproc_left,nproc-1,nproc_per_kpt
              ikpt=ikpt+1
              !print *,'passed through here',modulo(nproc,nkpts),nkpts,ikpt,nproc_per_kpt,nproc,jproc,nproc_left
              if (ikpt > nkpts .or. jproc > nproc-1) stop 'ERROR: also this should not happen3b'
              do iobj=0,nobj-1
                 nobj_par(jproc+modulo(iobj,nproc_per_kpt),ikpt)=nobj_par(jproc+modulo(iobj,nproc_per_kpt),ikpt)+1
              end do
           end do
        end if
        !print *,'passed through here',modulo(nproc,nkpts),nkpts,ikpt,nproc_per_kpt,nproc
     else
        !decide in how many kpoints single processor can be partitioned
        nkpt_per_proc=max((nkpts-1),1)/nproc !this is the minimum
        !count how many k-points are left that way
        !distribute the processors among these
        nkpts_left=nkpts-nkpt_per_proc*nproc
        ikpt=1
        jproc=-1
        !print *,'qui',nkpts_left,nkpts_per_proc
        do jkpt=1,nkpts_left
           jproc=jproc+1
           if (jproc > nproc-1) stop 'ERROR: also this should not happen4'
           do iobj=0,(nobj)*(nkpt_per_proc+1)-1
              nobj_par(jproc,ikpt+modulo(iobj,nkpt_per_proc+1))=nobj_par(jproc,ikpt+modulo(iobj,nkpt_per_proc+1))+1
           end do
           ikpt=ikpt+nkpt_per_proc+1
        end do
        !print *,'ciao'
        if ((nkpt_per_proc+1)*nkpts_left < nkpts) then
           do ikpt=(nkpt_per_proc+1)*nkpts_left+1,nkpts,nkpt_per_proc
              jproc=jproc+1
              !print *,'passed through here',modulo(nproc,nkpts),nkpts,ikpt,nproc_per_kpt,nproc,jproc,nproc_left
              if (ikpt > nkpts .or. jproc > nproc-1) stop 'ERROR: also this should not happen4b'
              do iobj=0,(nobj)*(nkpt_per_proc)-1
                 nobj_par(jproc,ikpt+modulo(iobj,nkpt_per_proc))=nobj_par(jproc,ikpt+modulo(iobj,nkpt_per_proc))+1
              end do
           end do
        end if
           !print *,'passed through here',modulo(nproc,nkpts),nkpts,ikpt,nproc_per_kpt,nproc
     end if
  else
     !non-integer repartition
     iobj=0
     ikpt=0
     do jproc=0,nproc-2 !leave the last processor at the end
        nobjp=floor(robjp)
        !respect the rounding ratio
        if (nproc-jproc > nprocs_with_floor) nobjp=nobjp+1
        !print *,'jproc,nobjp',jproc,nobjp,nkpts,nobj,nkpts*nobj,iobj,nprocs_with_floor
        do jobj=1,nobjp
           if (modulo(iobj,nobj) ==0) ikpt=ikpt+1
           iobj=iobj+1
           if (iobj > nobj*nkpts) stop 'ERROR: also this should not happen'
           nobj_par(jproc,ikpt)=nobj_par(jproc,ikpt)+1
        end do
     end do
     !in the last processor we put the objects which are lacking
     nobjp=nobj*nkpts-iobj
     do jobj=1,nobjp
        if (modulo(iobj,nobj) ==0) ikpt=ikpt+1
        iobj=iobj+1
        !print *,'finished',jobj,nobjp,iobj,nobj*nkpts,jproc,ikpt
        if (iobj > nobj*nkpts) stop 'ERROR: also this should not happen2'
        nobj_par(nproc-1,ikpt)=nobj_par(nproc-1,ikpt)+1
     end do
  end if
end subroutine kpts_to_procs_via_obj

subroutine components_kpt_distribution(nproc,nkpts,norb,nvctr,norb_par,nvctr_par)
  use module_base
  implicit none
  integer, intent(in) :: nproc,nkpts,nvctr,norb
  integer, dimension(0:nproc-1,nkpts), intent(in) :: norb_par
  integer, dimension(0:nproc-1,nkpts), intent(out) :: nvctr_par
  !local variables
  integer :: ikpt,jsproc,jeproc,kproc,icount,ivctr,jproc,numproc
  real(gp) :: strprc,endprc

  !for any of the k-points find the processors which have such k-point associated
  call to_zero(nproc*nkpts,nvctr_par(0,1))

  do ikpt=1,nkpts
     jsproc=UNINITIALIZED(1)
     jeproc=UNINITIALIZED(1)
     find_start: do jproc=0,nproc-1
        if(norb_par(jproc,ikpt) > 0) then 
           jsproc=jproc
           exit find_start
        end if
     end do find_start
     if (jsproc == UNINITIALIZED(1)) stop 'ERROR in kpt assignments'
     if(norb_par(jsproc,ikpt) /= norb) then
        strprc=real(norb_par(jsproc,ikpt),gp)/real(norb,gp)     
     else
        strprc=1.0_gp
     end if
     if (ikpt < nkpts) then
        find_end: do jproc=jsproc,nproc-1
           if(norb_par(jproc,ikpt+1) > 0) then
              if (norb_par(jproc,ikpt)==0) then
                 jeproc=jproc-1
              else
                 jeproc=jproc
              end if
              exit find_end
           end if
        end do find_end
        if (jeproc == UNINITIALIZED(1)) stop 'ERROR in kpt assignments'
     else
        jeproc=nproc-1
     end if
     if (jeproc /= jsproc) then
        endprc=real(norb_par(jeproc,ikpt),gp)/real(norb,gp)     
     else
        endprc=0.0_gp
     end if
     !if the number of processors is bigger than the number of orbitals this means 
     !that strprc and endprc are not correctly evaluated
     !evaluate the percentace on the number of components
     if (jeproc-jsproc+1 > norb) then
        strprc=1.0_gp/real(jeproc-jsproc+1,gp)
        endprc=strprc
     end if
     !assign the number of components which corresponds to the same orbital distribution
     numproc=jeproc-jsproc+1
     ivctr=0

     !print *,'kpoint',ikpt,jsproc,jeproc,strprc,endprc,ceiling(strprc*real(nvctr,gp)),nvctr
     !start filling the first processor
     nvctr_par(jsproc,ikpt)=min(ceiling(strprc*real(nvctr,gp)),nvctr)
     ivctr=min(ceiling(strprc*real(nvctr,gp)),nvctr)
     fill_array: do 
        if (ivctr==nvctr) exit fill_array
        icount=icount+1
        kproc=jsproc+modulo(icount,numproc)
        !put the floor of the components to the first processor
        if (strprc /= 1.0_gp .and. kproc==jsproc .and. nvctr_par(kproc,ikpt)==ceiling(strprc*real(nvctr,gp))) then
           !do nothing, skip away
        else
           nvctr_par(kproc,ikpt)=&
                nvctr_par(kproc,ikpt)+1
           ivctr=ivctr+1
        end if
     end do fill_array
  end do

end subroutine components_kpt_distribution

subroutine check_kpt_distributions(nproc,nkpts,norb,ncomp,norb_par,ncomp_par,info,lub_orbs,lub_comps)
  use module_base
  implicit none
  integer, intent(in) :: nproc,nkpts,norb,ncomp
  integer, dimension(0:nproc-1,nkpts), intent(in) :: norb_par
  integer, dimension(0:nproc-1,nkpts), intent(in) :: ncomp_par
  integer, intent(inout) :: info
  integer, intent(out) :: lub_orbs,lub_comps
  !local variables
  character(len=*), parameter :: subname='check_kpt_distributions'
  logical :: yesorb,yescomp,notcompatible
  integer :: ikpt,jorb,jproc,ierr,norbs,ncomps,i_all,i_stat,kproc,ieproc,isproc
  integer, dimension(:,:), allocatable :: load_unbalancing
  !before printing the distribution schemes, check that the two distributions contain
  !the same k-points
  if (info == 0) call print_distribution_schemes(6,nproc,nkpts,norb_par,ncomp_par)

  allocate(load_unbalancing(0:nproc-1,2+ndebug),stat=i_stat)
  call memocc(i_stat,load_unbalancing,'load_unbalancing',subname)

  do ikpt=1,nkpts
     isproc=UNINITIALIZED(1)
     find_isproc : do kproc=0,nproc-1
        if (ncomp_par(kproc,ikpt) > 0) then
           isproc=kproc
           exit find_isproc
        end if
     end do find_isproc
     if (isproc == UNINITIALIZED(1)) stop 'ERROR(check_kpt_distributions): isproc cannot be found'
     ieproc=UNINITIALIZED(1)
     find_ieproc : do kproc=nproc-1,0,-1
        if (ncomp_par(kproc,ikpt) > 0) then
           ieproc=kproc
           exit find_ieproc
        end if
     end do find_ieproc
     if (ieproc == UNINITIALIZED(1)) stop 'ERROR(check_kpt_distributions): ieproc cannot be found'

     norbs=0
     ncomps=0
     do jproc=0,nproc-1
        !count the total number of components
        norbs=norbs+norb_par(jproc,ikpt)
        ncomps=ncomps+ncomp_par(jproc,ikpt)
        notcompatible=(ncomp_par(jproc,ikpt) == 0 .neqv. norb_par(jproc,ikpt) == 0) 
        !check whether there are only 0 orbitals
        if (notcompatible .and. norb_par(jproc,ikpt)==0) then
           if (isproc < jproc .and. jproc <= ieproc) notcompatible=.false.
        end if
        if (notcompatible) then     
           if (info == 0) write(*,*)' ERROR: processor ', jproc,' kpt,',ikpt,&
                'have components and orbital distributions not compatible'
           info=1
           return
           !call MPI_ABORT(MPI_COMM_WORLD, ierr)
        end if
     end do
     if (norb/=norbs .or. ncomps /= ncomp) then
        if (info == 0) write(*,*)' ERROR: kpt,',ikpt,&
             'has components or orbital distributions not correct'
        info=2
        return
        !call MPI_ABORT(MPI_COMM_WORLD, ierr)
     end if
  end do
  do jproc=0,nproc-1
     load_unbalancing(jproc,:)=0
     do ikpt=1,nkpts
        load_unbalancing(jproc,1)=load_unbalancing(jproc,1)+norb_par(jproc,ikpt)
        load_unbalancing(jproc,2)=load_unbalancing(jproc,2)+ncomp_par(jproc,ikpt)
     end do
  end do

  !calculate the maximum load_unbalancing
  lub_orbs=0
  lub_comps=0
  do jproc=0,nproc-1
     do kproc=0,nproc-1
        lub_orbs=max(lub_orbs,load_unbalancing(jproc,1)-load_unbalancing(kproc,1))
        lub_comps=max(lub_comps,load_unbalancing(jproc,2)-load_unbalancing(kproc,2))
     end do
  end do

  if (info==0) write(*,*)' Kpoints Distribuitions are compatible, load unbalancings, orbs,comps:',lub_orbs,&
       '/',max(minval(load_unbalancing(:,1)),1),lub_comps,'/',minval(load_unbalancing(:,2))
  info=0
  i_all=-product(shape(load_unbalancing))*kind(load_unbalancing)
  deallocate(load_unbalancing,stat=i_stat)
  call memocc(i_stat,i_all,'load_unbalancing',subname)


end subroutine check_kpt_distributions

!>routine which associates to any of the processor a given number of objects
!! depending of the number of processors and k-points
subroutine parallel_repartition_with_kpoints(nproc,nkpts,nobj,nobj_par)
  use module_base
  implicit none
  integer, intent(in) :: nkpts,nobj,nproc
  integer, dimension(0:nproc-1), intent(out) :: nobj_par
  !local variables
  integer :: n_i,n_ip,rs_i,N_a,N_b,N_c,ikpt,jproc,i,ntmp
!!$  real(gp) :: rtmp

  ! Strategy to divide between k points.
  ! There is an nproc length to divide into orbs%nkpts segments.
  ! Segment (ikpt - 1) expand in 0 <= r_i < r_ip <= nproc.
  ! where r_i and r_ip are real values. There are two possibilities:
  !  - We can write r_i <= n_i <= n_ip <= r_ip with n_i and n_ip integers ;
  !  - or r_i <= n_i and n_ip <= r_ip and n_i = n_ip + 1.
  ! For both cases, we can divide nobj into the partition (real values):
  !  - N_a = (n_i - r_i)*nobj*nkpts/nproc (the initial part);
  !  - N_b = max((n_ip - n_i)*nobj*nkpts / nproc, 0) (the naive part, the only one if nkpts is a multiple of nproc);
  !  - N_c = (r_ip - n_ip) * nobj * orbs%nkpts / nproc (the final part);
  ! Before going to integer values, we have r_i = (ikpt - 1) * nproc / orbs%nkpts (the naive division)
  ! and r_ip = (ikpt) * nproc / orbs%nkpts (the segment endpoint)
  ! So N_a and N_b can be simplified and written instead:
  !  - N_a = int(nobj * (n_i * orbs%nkpts - (ikpt - 1) * nproc) / nproc);
  !  - N_c = int(nobj * ((ikpt) * nproc - n_ip * orbs%nkpts) / nproc)
  !  - N_b = nobj - N_a - N_c 
  ! After, if N_a > 0, we put this quantity to proc n_i - 1, if N_c > 0
  ! we put its quantity to proc n_ip ; and finally N_b is distributed
  ! among [n_i;n_ip[ procs.

  nobj_par(:)=0
  do ikpt=1,nkpts
     ! Calculation of n_i and n_ip, rs_i = r_i * orbs%nkpts to avoid rounding.
     rs_i=(ikpt-1)*nproc !integer variable for rounding purposes

     if (mod(rs_i,nkpts) == 0) then
        n_i=rs_i/nkpts 
     else
        n_i=rs_i/nkpts+1
     end if

     rs_i=ikpt*nproc
     n_ip=rs_i/nkpts
!!$     print *,'ikpt,ni,nip',ikpt,n_i,n_ip
     ! Calculation of N_a, N_b and N_c from given n_i and n_ip.
     if (n_ip >= n_i) then
        ntmp = (n_i*nkpts-(ikpt-1)*nproc) * nobj
        if (modulo(ntmp, nproc) == 0) then
           N_a = ntmp / nproc
        else
           N_a = (ntmp - modulo(ntmp, nproc) + nproc) / nproc
        end if
!!$        ntmp=n_i*nkpts-(ikpt-1)*nproc
!!$        rtmp=real(nobj,gp)/real(nproc,gp)
!!$        rtmp=rtmp*real(ntmp,gp)
!!$        N_a=floor(rtmp)
!!$        if (iproc == 0) print *,'ikpts,rtmp',ikpt,rtmp
        ntmp = (ikpt*nproc-n_ip*nkpts) * nobj
        if (modulo(ntmp, nproc) == 0) then
           N_c = ntmp / nproc
        else
           N_c = (ntmp - modulo(ntmp, nproc) + nproc) / nproc
        end if

!!$        ntmp=ikpt*nproc-n_ip*nkpts
!!$        rtmp=real(nobj,gp)/real(nproc,gp)
!!$        rtmp=rtmp*real(ntmp,gp)
!!$        N_c=ceiling(rtmp)
!!$        if (iproc == 0) print *,'ikpts,rtmp2',ikpt,rtmp,N_a,N_c
        !the corrections above are to avoid the 32 bit integer overflow
        !N_a=nint(real(nobj*(n_i*nkpts-(ikpt-1)*nproc),gp)/real(nproc,gp))
        !N_c=nint(real(nobj*(ikpt*nproc-n_ip*nkpts),gp)/real(nproc,gp))
     else
        N_c=nobj/2
        N_a=nobj-N_c
     end if
     N_b=nobj-N_a-N_c
     if (N_b == -1) then
        N_c = N_c - 1
        N_b = 0
     end if
!!$     write(*,*) ikpt, N_a, N_b, N_c
     if (nkpts > 1 .and. N_b < n_ip - n_i) stop 'ERROR:parallel_repartion_with_kpoints'
     !assign to procs the objects.
     if (N_a>0) nobj_par(n_i-1)=nobj_par(n_i-1)+N_a
     if (N_b>0) then
        do i=0,N_b-1
           jproc=n_i+mod(i,n_ip-n_i)
           nobj_par(jproc)=nobj_par(jproc)+1
        end do
     end if
     if (N_c>0) nobj_par(n_ip)=nobj_par(n_ip)+N_c
  end do
END SUBROUTINE parallel_repartition_with_kpoints


subroutine parallel_repartition_per_kpoints(iproc,nproc,nkpts,nobj,nobj_par,&
     nkptsp,mykpts,nobj_pkpt)
  implicit none
  integer, intent(in) :: iproc,nproc,nkpts,nobj
  integer, dimension(0:nproc-1), intent(in) :: nobj_par
  integer, intent(out) :: nkptsp
  integer, dimension(nkpts), intent(out) :: mykpts
  integer, dimension(0:nproc-1,nkpts), intent(out) :: nobj_pkpt
  !local variables
  integer :: ikpts,jproc,jobj,norb_tot,iorbp

  !initialise the array
  do ikpts=1,nkpts
     do jproc=0,nproc-1
        nobj_pkpt(jproc,ikpts)=0 
     end do
  end do

  !assign the k-point, counting one object after each other
  jobj=1
  ikpts=1
  !print *,'here',nobj_par(:)
  do jproc=0,nproc-1
     do iorbp=1,nobj_par(jproc)
        nobj_pkpt(jproc,ikpts)=nobj_pkpt(jproc,ikpts)+1
        if (mod(jobj,nobj)==0) then
           ikpts=ikpts+1
        end if
        jobj=jobj+1
     end do
  end do
  !some checks
  if (nobj /= 0) then
     !check the distribution
     do ikpts=1,nkpts
        !print *,'partition',ikpts,orbs%nkpts,'ikpts',nobj_pkpt(:,ikpts)
        norb_tot=0
        do jproc=0,nproc-1
           norb_tot=norb_tot+nobj_pkpt(jproc,ikpts)
        end do
        if(norb_tot /= nobj) then
           write(*,*)'ERROR: partition of objects incorrect, kpoint:',ikpts
           stop
        end if
     end do
  end if

  !calculate the number of k-points treated by each processor in both
  ! the component distribution and the orbital distribution.
  nkptsp=0
  do ikpts=1,nkpts
     if (nobj_pkpt(iproc,ikpts) /= 0) then
        nkptsp=nkptsp+1
        mykpts(nkptsp) = ikpts
     end if
  end do

END SUBROUTINE parallel_repartition_per_kpoints
