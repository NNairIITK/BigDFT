subroutine print_logo()
  implicit none
  write(*,'(23x,a)')'      ****         *       *****    '
  write(*,'(23x,a)')'     *    *               *         '
  write(*,'(23x,a)')'    *     *        *     *          '
  write(*,'(23x,a)')'    *    *         *     *        * '
  write(*,'(23x,a)')'    *****          *     *         *'
  write(*,'(23x,a)')'    *    *         *     *         *'
  write(*,'(23x,a)')'    *     *        *     *         *'
  write(*,'(23x,a)')'    *      *       *     *         *'
  write(*,'(23x,a)')'    *     *     ****     *         *'
  write(*,'(23x,a)')'    * ****         *      *        *'
  write(*,'(23x,a)')'    *             *        *      * '
  write(*,'(23x,a)')'*********    *****          ******  ' 
  !write(*,'(23x,a)')'---------------------------------------'
  write(*,'(23x,a)')'  ******          *****    *********'
  write(*,'(23x,a)')' *      *        *             *    '
  write(*,'(23x,a)')'*        *      *         **** *    '
  write(*,'(23x,a)')'*         *     ****     *     *    '
  write(*,'(23x,a)')'*         *     *       *      *    '
  write(*,'(23x,a)')'*         *     *        *     *    '
  write(*,'(23x,a)')'*         *     *         *    *    '
  write(*,'(23x,a)')'*         *     *          *****    '
  write(*,'(23x,a)')' *        *     *         *    *    '  
  write(*,'(23x,a)')'          *     *        *     *    ' 
  write(*,'(23x,a)')'         *               *    *     '
  write(*,'(23x,a)')'    *****       *         ****                       (Ver 1.0)'
  write(*,'(1x,a)')&
       '------------------------------------------------------------------------------------'
  write(*,'(1x,a)')&
       '|              Daubechies Wavelets for DFT Pseudopotential Calculations            |'
  write(*,'(1x,a)')&
          '------------------------------------------------------------------------------------'
end subroutine print_logo

subroutine read_input_variables(iproc,in)

  use module_types

  implicit none
  integer, intent(in) :: iproc
  type(input_variables), intent(out) :: in
  !local variables
  real(kind=4) :: hgrid,crmult,frmult,cpmult,fpmult
  integer :: ierror

  ! Read the input variables.
  open(unit=1,file='input.dat',status='old')
  read(1,*,iostat=ierror) in%ncount_cluster_x
  read(1,*,iostat=ierror) in%frac_fluct
  read(1,*,iostat=ierror) in%randdis
  read(1,*,iostat=ierror) in%betax
  read(1,*,iostat=ierror) hgrid
  read(1,*,iostat=ierror) crmult
  read(1,*,iostat=ierror) frmult
  read(1,*,iostat=ierror) cpmult
  read(1,*,iostat=ierror) fpmult
  in%hgrid  = real(hgrid,kind=8)
  in%crmult = real(crmult,kind=8)
  in%frmult = real(frmult,kind=8)
  in%cpmult = real(cpmult,kind=8)
  in%fpmult = real(fpmult,kind=8)
  if (in%fpmult.gt.in%frmult .and. iproc==0) write(*,*) ' NONSENSE: fpmult > frmult'
  read(1,*,iostat=ierror) in%ixc
  read(1,*,iostat=ierror) in%ncharge,in%elecfield
  read(1,*,iostat=ierror) in%gnrm_cv
  read(1,*,iostat=ierror) in%itermax
  read(1,*,iostat=ierror) in%ncong
  read(1,*,iostat=ierror) in%idsx
  read(1,*,iostat=ierror) in%calc_tail
  read(1,*,iostat=ierror) in%rbuf
  read(1,*,iostat=ierror) in%ncongt
  read(1,*,iostat=ierror) in%nspin,in%mpol

  if (ierror/=0) then
     write(*,'(1x,a)') 'Error while reading the file "input.dat"'
     stop
  end if

  !add optional line for the restart option
  read(1,*,iostat=ierror) in%inputPsiId,in%output_wf,in%output_grid

  close(1,iostat=ierror)

!!$  !these values are hard-coded for the moment but they can be entered in the input file
!!$  !in case of need also other variables can be entered without any changements
!!$  in%output_grid=.false. 
!!$  in%inputPsiId=0
!!$  in%output_wf=.false. 

  if (iproc == 0) then
     write(*,'(1x,a,i0)') 'Max. number of wavefnctn optim ',in%ncount_cluster_x
     write(*,'(1x,a,1pe10.2)') 'Convergence criterion for forces: fraction of noise ',&
          in%frac_fluct
     write(*,'(1x,a,1pe10.2)') 'Random displacement amplitude ',in%randdis
     write(*,'(1x,a,1pe10.2)') 'Steepest descent step ',in%betax
  end if

     if (in%nspin==2) then
        if (iproc == 0) write(*,'(1x,a)') 'Spin-polarised calculation: YES (Collinear)'
     else if (in%nspin==1) then
        if (iproc == 0) write(*,'(1x,a)') 'Spin-polarised calculation:  NO '
     else
        if (iproc == 0) write(*,'(1x,a,i0)')'Wrong spin polarisation id: ',in%nspin
        stop
     end if

end subroutine read_input_variables

subroutine print_input_parameters(in)
  use module_types
  implicit none
  type(input_variables), intent(in) :: in

  write(*,'(1x,a)')&
       '------------------------------------------------------------------- Input Parameters'
  write(*,'(1x,a)')&
       '    System Choice       Resolution Radii        SCF Iteration      Finite Size Corr.'
  write(*,'(1x,a,f7.3,1x,a,f5.2,1x,a,1pe8.1,1x,a,l4)')&
       'Grid spacing=',in%hgrid,   '|  Coarse Wfs.=',in%crmult,'| Wavefns Conv.=',in%gnrm_cv,&
       '| Calculate=',in%calc_tail
  write(*,'(1x,a,i7,1x,a,f5.2,1x,a,i8,1x,a,f4.1)')&
       '       XC id=',in%ixc,     '|    Fine Wfs.=',in%frmult,'| Max. N. Iter.=',in%itermax,&
       '| Extension=',in%rbuf
  write(*,'(1x,a,i7,1x,a,f5.2,1x,a,i8,1x,a,i4)')&
       'total charge=',in%ncharge, '| Coarse Proj.=',in%cpmult,'| CG Prec.Steps=',in%ncong,&
       '|  CG Steps=',in%ncongt
  write(*,'(1x,a,1pe7.1,1x,a,0pf5.2,1x,a,i8)')&
       ' elec. field=',in%elecfield,'|   Fine Proj.=',in%fpmult,'| DIIS Hist. N.=',in%idsx
  if (in%nspin==2) then
     write(*,'(1x,a,i7,1x,a)')&
          'Polarisation=',2*in%mpol, '|'
  end if
end subroutine print_input_parameters

! read atomic positions
subroutine read_atomic_positions(iproc,ifile,units,nat,ntypes,iatype,atomnames,rxyz)
  implicit none
  character(len=20), intent(in) :: units
  integer, intent(in) :: iproc,ifile,nat
  integer, intent(out) :: ntypes
  character(len=20), dimension(100), intent(out) :: atomnames
  integer, dimension(nat), intent(out) :: iatype
  real(kind=8), dimension(3,nat), intent(out) :: rxyz
  !local variables
  character(len=2) :: symbol
  character(len=3) :: suffix
  character(len=20) :: tatonam
  character(len=100) :: line
  integer :: nateq,iat,jat,ityp,i,ierror,ierrsfx
  real(kind=8), parameter :: bohr=0.5291772108d0 !1 AU in angstroem
! To read the file posinp (avoid differences between compilers)
  real(kind=4) :: rx,ry,rz

  if (iproc.eq.0) write(*,'(1x,a,i0)') 'Number of atoms     = ',nat

  !read from positions of .xyz format, but accepts also the old .ascii format
  read(ifile,'(a100)')line

  read(line,*,iostat=ierror) rx,ry,rz,tatonam

  ntypes=0
  do iat=1,nat
     if (ierror == 0) then
        if (iat /= 1) read(ifile,*) rx,ry,rz,tatonam
     else
        read(ifile,'(a100)')line 
        read(line,*,iostat=ierrsfx)symbol,rx,ry,rz,suffix
        if (ierrsfx ==0) then
           tatonam=trim(symbol)//'_'//trim(suffix)
        else
           read(line,*)symbol,rx,ry,rz
           tatonam=trim(symbol)
        end if
     end if
     rxyz(1,iat)=real(rx,kind=8)
     rxyz(2,iat)=real(ry,kind=8)
     rxyz(3,iat)=real(rz,kind=8)
     !! For reading in saddle points
     !!        open(unit=83,file='step',status='old')
     !!        read(83,*) step
     !!        close(83)
     !        step= .5d0
     !        if (iproc.eq.0) write(*,*) 'step=',step
     !        read(99,*) rxyz(1,iat),rxyz(2,iat),rxyz(3,iat),tatonam,t1,t2,t3
     !        rxyz(1,iat)=rxyz(1,iat)+step*t1
     !        rxyz(2,iat)=rxyz(2,iat)+step*t2
     !        rxyz(3,iat)=rxyz(3,iat)+step*t3

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
           rxyz(i,iat)=rxyz(i,iat)/bohr
        enddo
     else if  (units.eq.'atomic' .or. units.eq.'bohr') then
     else
        write(*,*) 'length units in input file unrecognized'
        write(*,*) 'recognized units are angstroem or atomic = bohr'
        stop 
     endif
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

  if (iproc.eq.0) write(*,'(1x,a,i0)') 'Number of atom types= ',ntypes

  do ityp=1,ntypes
     if (iproc.eq.0) &
          write(*,'(1x,a,i0,a,a)') 'Atoms of type ',ityp,' are ',trim(atomnames(ityp))
  enddo

end subroutine read_atomic_positions


! Fill the arrays occup and spinar
! if iunit /=0 this means that the file occup.dat does exist and it opens
subroutine input_occup(iproc,iunit,nelec,norb,norbu,norbd,nspin,mpol,occup,spinar)
  implicit none
! Arguments
  integer, intent(in) :: nelec,nspin,mpol,iproc,norb,norbu,norbd,iunit
  real(kind=8), intent(out) :: occup(norb),spinar(norb)
! Local variables
  integer :: iorb,nt,ne,it,ierror,iorb1
  real(kind=8) :: rocc,rup,rdown
  character(len=100) :: line

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
                 write(*,'(10x,a,i0,a)') 'The orbital index ',iorb,' is incorrect'
              end if
              stop
           elseif (rocc<0.d0 .or. rocc>2.d0) then
              if (iproc==0) then
                 write(*,'(1x,a,i0,a)') 'ERROR in line ',nt+1,' of the file "occup.dat"'
                 write(*,'(10x,a,f5.2,a)') 'The occupation number ',rocc,' is not between 0. and 2.'
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
     if (nspin/=1) then
        !Check if the polarisation is respected (mpol)
        rup=sum(occup(1:norbu))
        rdown=sum(occup(norbu+1:norb))
        if (abs(rup-rdown-real(mpol,kind=8))>1.d-6) then
           if (iproc==0) then
              write(*,'(1x,a,f13.6,a,i0)') 'From the file "occup.dat", the polarization ',rup-rdown,&
                             ' is not equal to ',mpol
           end if
           stop
        end if
        !Fill spinar
        do iorb=1,norbu
           spinar(iorb)=1.0d0
        end do
        do iorb=norbu+1,norb
           spinar(iorb)=-1.0d0
        end do
      end if
  end if
  if (iproc==0) then 
     write(*,'(1x,a,i8)') &
          'Total Number of Orbitals ',norb
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

subroutine read_system_variables(iproc,nproc,nat,ntypes,nspin,ncharge,mpol,hgrid,atomnames,iatype,&
     psppar,radii_cf,npspcode,iasctype,nelpsp,nzatom,nelec,natsc,norb,norbu,norbd,norbp,iunit)
  implicit none
  integer, intent(in) :: iproc,nproc,nat,ntypes,nspin,ncharge,mpol
  real(kind=8), intent(in) :: hgrid
  character(len=20), dimension(ntypes), intent(in) :: atomnames
  integer, dimension(nat), intent(in) :: iatype
  integer, intent(out) :: nelec,natsc,norb,norbu,norbd,norbp,iunit
  integer, dimension(ntypes), intent(out) :: npspcode,iasctype,nelpsp,nzatom
  real(kind=8), dimension(ntypes,2), intent(out) :: radii_cf
  real(kind=8), dimension(0:4,0:6,ntypes), intent(out) :: psppar
  !local variables
  real(kind=8), parameter :: eps_mach=1.d-12
  logical :: exists
  character(len=2) :: symbol
  character(len=27) :: filename
  integer :: i,j,l,iat,nlterms,nprl,nn,nt,ntu,ntd,ityp,ierror,i_stat,i_all
  real(kind=8) :: rcov,rprb,ehomo,radfine,tt,minrad
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
       !control the hardest gaussian
       minrad=1.d10
       do i=0,4
          if (psppar(i,0,ityp)/=0.d0) then
             minrad=min(minrad,psppar(i,0,ityp))
          end if
       end do
       !control whether the grid spacing is too high or not
       if (iproc == 0 .and. hgrid > 2.5d0*minrad) then
          write(*,'(1x,a)')'WARNING: The grid spacing value may be too high to treat correctly the above pseudo.' 
          write(*,'(1x,a,f5.2,a)')'         Results can be meaningless if hgrid is bigger than',2.5d0*minrad,'. At your own risk!'
       end if

    enddo



    !deallocation
    i_all=-product(shape(neleconf))*kind(neleconf)
    deallocate(neleconf,stat=i_stat)
    call memocc(i_stat,i_all,'neleconf','read_PSP_variables')


    !calculate number of electrons and orbitals
    ! Number of electrons and number of semicore atoms
    nelec=0
    natsc=0
    do iat=1,nat
       ityp=iatype(iat)
       nelec=nelec+nelpsp(ityp)
       if (iasctype(ityp) /= 0) natsc=natsc+1
    enddo
    nelec=nelec-ncharge
    if (iproc.eq.0) then
       write(*,'(1x,a,i8)') &
            'Total Number of Electrons ',nelec
    end if

    ! Number of orbitals
    if (nspin==1) then
       norb=(nelec+1)/2
       norbu=norb
       norbd=0
       if (mod(nelec,2).ne.0 .and. iproc==0) then
          write(*,'(1x,a)') 'WARNING: odd number of electrons, no closed shell system'
       end if
    else
       if (iproc==0) write(*,'(1x,a)') 'Spin-polarized calculation'
       norb=nelec
       norbu=min(norb/2+mpol,norb)
       norbd=norb-norbu
    end if

    ! Test if the file 'occup.dat exists
    inquire(file='occup.dat',exist=exists)
    iunit=0
    if (exists) then
       iunit=25
       open(unit=iunit,file='occup.dat',form='formatted',action='read',status='old')
       if (nspin==1) then
           !The first line gives the number of orbitals
           read(unit=iunit,fmt=*,iostat=ierror) nt
       else
           !The first line gives the number of orbitals
           read(unit=iunit,fmt=*,iostat=ierror) ntu,ntd
       end if
       if (ierror /=0) then
          if (iproc==0) write(*,'(1x,a)') &
               'ERROR reading the number of orbitals in the file "occup.dat"'
          stop
       end if
       !Check
       if (nspin==1) then
          if (nt<norb) then
              if (iproc==0) write(*,'(1x,a,i0,a,i0)') &
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
              if (iproc==0) write(*,'(1x,a,i0,a,i0)') &
                  'ERROR: In the file "occup.dat", the number of orbitals norb=',nt,&
                  ' should be greater or equal than nelec=',norb
             stop
          else
             norb=nt
          end if
          if (ntu<norbu) then
              if (iproc==0) write(*,'(1x,a,i0,a,i0)') &
                  'ERROR: In the file "occup.dat", the number of orbitals up norbu=',ntu,&
                  ' should be greater or equal than min(nelec/2+mpol,nelec)=',norbu
             stop
          else
             norbu=ntu
          end if
          if (ntd<norbd) then
              if (iproc==0) write(*,'(1x,a,i0,a,i0)') &
              'ERROR: In the file "occup.dat", the number of orbitals down norbd=',ntd,&
                  ' should be greater or equal than min(nelec/2-mpol,0)=',norbd
             stop
          else
             norbd=ntd
          end if
       end if
    end if

    tt=dble(norb)/dble(nproc)
    norbp=int((1.d0-eps_mach*tt) + tt)
    !if (iproc.eq.0) write(*,'(1x,a,1x,i0)') 'norbp=',norbp

end subroutine read_system_variables


subroutine system_size(iproc,nat,ntypes,rxyz,radii_cf,crmult,frmult,hgrid,iatype,atomnames, &
       alat1,alat2,alat3,n1,n2,n3,nfl1,nfl2,nfl3,nfu1,nfu2,nfu3)
  ! calculates the overall size of the simulation cell (cxmin,cxmax,cymin,cymax,czmin,czmax)
  !and shifts the atoms such that their position is the most symmetric possible
  implicit none
  integer, intent(in) :: iproc,nat,ntypes
  real(kind=8), intent(in) :: hgrid,crmult,frmult
  character(len=20), dimension(ntypes), intent(in) :: atomnames
  integer, dimension(nat), intent(in) :: iatype
  real(kind=8), dimension(3,nat), intent(inout) :: rxyz
  real(kind=8), dimension(ntypes,2), intent(in) :: radii_cf
  integer, intent(out) :: n1,n2,n3,nfl1,nfl2,nfl3,nfu1,nfu2,nfu3
  real(kind=8), intent(out) :: alat1,alat2,alat3
  !local variables
  real(kind=8), parameter ::eps_mach=1.d-12,onem=1.d0-eps_mach
  integer :: iat,j
  real(kind=8) :: rad,cxmin,cxmax,cymin,cymax,czmin,czmax,alatrue1,alatrue2,alatrue3

  cxmax=-1.d10 ; cxmin=1.d10
  cymax=-1.d10 ; cymin=1.d10
  czmax=-1.d10 ; czmin=1.d10
  do iat=1,nat
     rad=radii_cf(iatype(iat),1)*crmult
     cxmax=max(cxmax,rxyz(1,iat)+rad) ; cxmin=min(cxmin,rxyz(1,iat)-rad)
     cymax=max(cymax,rxyz(2,iat)+rad) ; cymin=min(cymin,rxyz(2,iat)-rad)
     czmax=max(czmax,rxyz(3,iat)+rad) ; czmin=min(czmin,rxyz(3,iat)-rad)
  enddo

  cxmax=cxmax+eps_mach 
  cymax=cymax+eps_mach  
  czmax=czmax+eps_mach  

  cxmin=cxmin-eps_mach
  cymin=cymin-eps_mach
  czmin=czmin-eps_mach

  alat1=(cxmax-cxmin)
  alat2=(cymax-cymin)
  alat3=(czmax-czmin)

  ! grid sizes n1,n2,n3
  n1=int(alat1/hgrid)
  !if (mod(n1+1,4).eq.0) n1=n1+1
  n2=int(alat2/hgrid)
  !if (mod(n2+1,8).eq.0) n2=n2+1
  n3=int(alat3/hgrid)
  alatrue1=real(n1,kind=8)*hgrid 
  alatrue2=real(n2,kind=8)*hgrid 
  alatrue3=real(n3,kind=8)*hgrid

  !balanced shift taking into account the missing space
  cxmin=cxmin+0.5d0*(alat1-alatrue1)
  cymin=cymin+0.5d0*(alat2-alatrue2)
  czmin=czmin+0.5d0*(alat3-alatrue3)

  alat1=alatrue1
  alat2=alatrue2
  alat3=alatrue3

  do iat=1,nat
     rxyz(1,iat)=rxyz(1,iat)-cxmin
     rxyz(2,iat)=rxyz(2,iat)-cymin
     rxyz(3,iat)=rxyz(3,iat)-czmin
  enddo

!!$! grid sizes n1,n2,n3 !added for testing purposes
!!$  n1=int(alat1/hgrid)
!!$  if (mod(n1,2).eq.0) n1=n1+1
!!$  n2=int(alat2/hgrid)
!!$  if (mod(n2,2).eq.0) n2=n2+1
!!$  n3=int(alat3/hgrid)
!!$  if (mod(n3,2).eq.0) n3=n3+1
!!$  alat1=n1*hgrid 
!!$  alat2=n2*hgrid 
!!$  alat3=n3*hgrid
!!$  do iat=1,nat
!!$     rxyz(1,iat)=(real(n1/2,kind=8)+0.5)*hgrid 
!!$     rxyz(2,iat)=(real(n1/2,kind=8)+0.5)*hgrid 
!!$     rxyz(3,iat)=(real(n1/2,kind=8)+0.5)*hgrid 
!!$  enddo

  ! fine grid size (needed for creation of input wavefunction, preconditioning)
  nfl1=n1 ; nfl2=n2 ; nfl3=n3
  nfu1=0 ; nfu2=0 ; nfu3=0
  do iat=1,nat
     rad=radii_cf(iatype(iat),2)*frmult
     nfl1=min(nfl1,ceiling((rxyz(1,iat)-rad)/hgrid - eps_mach))
     nfu1=max(nfu1,floor((rxyz(1,iat)+rad)/hgrid + eps_mach))

     nfl2=min(nfl2,ceiling((rxyz(2,iat)-rad)/hgrid - eps_mach))
     nfu2=max(nfu2,floor((rxyz(2,iat)+rad)/hgrid + eps_mach))

     nfl3=min(nfl3,ceiling((rxyz(3,iat)-rad)/hgrid - eps_mach)) 
     nfu3=max(nfu3,floor((rxyz(3,iat)+rad)/hgrid + eps_mach))

!!$     nfl1=min(nfl1,int(onem+(rxyz(1,iat)-rad)/hgrid))
!!$     nfu1=max(nfu1,int((rxyz(1,iat)+rad)/hgrid))
!!$
!!$     nfl2=min(nfl2,int(onem+(rxyz(2,iat)-rad)/hgrid))
!!$     nfu2=max(nfu2,int((rxyz(2,iat)+rad)/hgrid))
!!$
!!$     nfl3=min(nfl3,int(onem+(rxyz(3,iat)-rad)/hgrid)) 
!!$     nfu3=max(nfu3,int((rxyz(3,iat)+rad)/hgrid))

  enddo

  if (iproc.eq.0) then
     write(*,'(1x,a,19x,a)') 'Shifted atomic positions, Atomic Units:','grid spacing units:'
     do iat=1,nat
        write(*,'(1x,i5,1x,a6,3(1x,1pe12.5),3x,3(1x,0pf9.3))') &
             iat,trim(atomnames(iatype(iat))),&
             (rxyz(j,iat),j=1,3),rxyz(1,iat)/hgrid,rxyz(2,iat)/hgrid,rxyz(3,iat)/hgrid
     enddo
     write(*,'(1x,a,3(1x,1pe12.5))') &
          '   Shift of=',-cxmin,-cymin,-czmin
     write(*,'(1x,a,3(1x,1pe12.5),3x,3(1x,i9))')&
          '  Box Sizes=',alat1,alat2,alat3,n1,n2,n3
     write(*,'(1x,a,3x,3(3x,i4,a1,i0))')&
          '      Extremes for the high resolution grid points:',&
          nfl1,'<',nfu1,nfl2,'<',nfu2,nfl3,'<',nfu3
  endif

end subroutine system_size
