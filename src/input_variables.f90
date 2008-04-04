subroutine print_logo()
  implicit none
  write(*,'(23x,a)')'      BBBB         i       ggggg    '
  write(*,'(23x,a)')'     B    B               g         '
  write(*,'(23x,a)')'    B     B        i     g          '
  write(*,'(23x,a)')'    B    B         i     g        g '
  write(*,'(23x,a)')'    BBBBB          i     g         g'
  write(*,'(23x,a)')'    B    B         i     g         g'
  write(*,'(23x,a)')'    B     B        i     g         g'
  write(*,'(23x,a)')'    B      B       i     g         g'
  write(*,'(23x,a)')'    B     B     iiii     g         g'
  write(*,'(23x,a)')'    B BBBB         i      g        g'
  write(*,'(23x,a)')'    B             i        g      g '
  write(*,'(23x,a)')'BBBBBBBBB    iiiii          gggggg  ' 
  !write(*,'(23x,a)')'---------------------------------------'
  write(*,'(23x,a)')'  DDDDDD          FFFFF    TTTTTTTTT'
  write(*,'(23x,a)')' D      D        F             T    '
  write(*,'(23x,a)')'D        D      F         TTTT T    '
  write(*,'(23x,a)')'D         D     FFFF     T     T    '
  write(*,'(23x,a)')'D         D     F       T      T    '
  write(*,'(23x,a)')'D         D     F        T     T    '
  write(*,'(23x,a)')'D         D     F         T    T    '
  write(*,'(23x,a)')'D         D     F          TTTTT    '
  write(*,'(23x,a)')' D        D     F         T    T    '  
  write(*,'(23x,a)')'          D     F        T     T    ' 
  write(*,'(23x,a)')'         D               T    T     '
  write(*,'(23x,a)')'    DDDDD       F         TTTT                     (Ver 1.0.1)'
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

     if (in%nspin==4) then
        if (iproc == 0) write(*,'(1x,a)') 'Spin-polarised calculation: YES (Non-collinear)'
     else if (in%nspin==2) then
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
  if (in%nspin>=2) then
     write(*,'(1x,a,i7,1x,a)')&
          'Polarisation=',2*in%mpol, '|'
  end if
end subroutine print_input_parameters

! read atomic positions
subroutine read_atomic_positions(iproc,ifile,units,at,rxyz)
  use module_types
  implicit none
  character(len=20), intent(in) :: units
  integer, intent(in) :: iproc,ifile
  type(atoms_data), intent(inout) :: at
  real(kind=8), dimension(3,at%nat), intent(out) :: rxyz
  !local variables
  real(kind=8), parameter :: bohr=0.5291772108d0 !1 AU in angstroem
  character(len=3) :: suffix
  character(len=2) :: symbol
  character(len=20) :: tatonam
  character(len=100) :: line
  logical :: lpsdbl 
  integer :: nateq,iat,jat,ityp,i,ierror,ierrsfx,i_stat,natpol
! To read the file posinp (avoid differences between compilers)
  real(kind=4) :: rx,ry,rz
! case for which the atomic positions are given whithin machine precision
  real(kind=8) :: rxd0,ryd0,rzd0
  character(len=20), dimension(100) :: atomnames

  if (iproc.eq.0) write(*,'(1x,a,i0)') 'Number of atoms     = ',at%nat

  allocate(at%iatype(at%nat),stat=i_stat)
  call memocc(i_stat,product(shape(at%iatype))*kind(at%iatype),'iatype','read_atomic_positions')
  allocate(at%lfrztyp(at%nat),stat=i_stat)
  call memocc(i_stat,product(shape(at%lfrztyp))*kind(at%lfrztyp),'lfrztyp','read_atomic_positions')
  allocate(at%nspinat(at%nat),stat=i_stat)
  call memocc(i_stat,product(shape(at%nspinat))*kind(at%nspinat),'nspinat','read_atomic_positions')

  !controls if the positions are provided with machine precision
  if (units == 'angstroemd0' .or. units== 'atomicd0' .or. units== 'bohrd0') then
     lpsdbl=.true.
  else
     lpsdbl=.false.
  end if

  !this array is useful for frozen atoms
  !no atom is frozen by default
  at%lfrztyp(:)=.false.
  !also the spin polarisation is fixed to zero by default
  at%nspinat(:)=0

  !read from positions of .xyz format, but accepts also the old .ascii format
  read(ifile,'(a100)')line

  !old format, still here for backward compatibility
  !admits only simple precision calculation
  read(line,*,iostat=ierror) rx,ry,rz,tatonam

  at%ntypes=0
  do iat=1,at%nat
     if (ierror == 0) then
        if (iat /= 1) read(ifile,*) rx,ry,rz,tatonam
     else
        read(ifile,'(a100)')line 
        if (lpsdbl) then
           read(line,*,iostat=ierrsfx)symbol,rxd0,ryd0,rzd0,natpol,suffix
        else
           read(line,*,iostat=ierrsfx)symbol,rx,ry,rz,natpol,suffix
        end if
        if (ierrsfx ==0) then
           at%nspinat(iat)=natpol
           !tatonam=trim(symbol)//'_'//trim(suffix)
           if (suffix == 'f') then
              !the atom is considered as blocked
              at%lfrztyp(iat)=.true.
           else
              if (iproc == 0) then
                 print *,suffix
                 write(*,'(1x,a,i0,a)')'ERROR in input file for atom number ',&
                      iat,': the only value accepted in 5th column is "f"'
                 stop
              end if
           end if
        else
           if (lpsdbl) then
              read(line,*,iostat=ierrsfx)symbol,rxd0,ryd0,rzd0,suffix
           else
              read(line,*,iostat=ierrsfx)symbol,rx,ry,rz,suffix
           end if
           if (ierrsfx ==0) then
              if (suffix == 'f') then
                 !the atom is considered as blocked
                 at%lfrztyp(iat)=.true.
              else
                 read(suffix,*,iostat=i_stat)at%nspinat(iat)
                 if(i_stat /=0) then
                    if (iproc == 0) then
                       print *,suffix
                       write(*,'(1x,a,i0,a)')'ERROR in input file for atom number ',&
                            iat,': in 4th column you can put the input polarisation or "f"'
                       stop
                    end if
                 end if
              end if
           else
              if (lpsdbl) then
                 read(line,*)symbol,rx,ry,rz
              else
                 read(line,*)symbol,rxd0,ryd0,rzd0
              end if
           end if
        end if
        tatonam=trim(symbol)
     end if
     if (lpsdbl) then
        rxyz(1,iat)=rxd0
        rxyz(2,iat)=ryd0
        rxyz(3,iat)=rzd0
     else
        rxyz(1,iat)=real(rx,kind=8)
        rxyz(2,iat)=real(ry,kind=8)
        rxyz(3,iat)=real(rz,kind=8)
     end if
 
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

     do ityp=1,at%ntypes
        if (tatonam.eq.atomnames(ityp)) then
           at%iatype(iat)=ityp
           goto 200
        endif
     enddo
     at%ntypes=at%ntypes+1
     if (at%ntypes.gt.100) stop 'more than 100 atomnames not permitted'
     atomnames(ityp)=tatonam
     at%iatype(iat)=at%ntypes
200  continue
     if (units=='angstroem' .or. units=='angstroemd0') then
        ! if Angstroem convert to Bohr
        do i=1,3 
           rxyz(i,iat)=rxyz(i,iat)/bohr
        enddo
     else if  (units=='atomic' .or. units=='bohr'  .or.&
          units== 'atomicd0' .or. units== 'bohrd0') then
     else
        write(*,*) 'length units in input file unrecognized'
        write(*,*) 'recognized units are angstroem or atomic = bohr'
        stop 
     endif
  enddo

  !now that ntypes is determined allocate at%atomnames and copy the values
  allocate(at%atomnames(at%ntypes),stat=i_stat) 
  call memocc(i_stat,product(shape(at%atomnames))*kind(at%atomnames),'atomnames','read_atomic_positions')
  at%atomnames(1:at%ntypes)=atomnames(1:at%ntypes)

  !control atom positions
  nateq=0
  do iat=1,at%nat
     do jat=iat+1,at%nat
        if ((rxyz(1,iat)-rxyz(1,jat))**2+(rxyz(2,iat)-rxyz(2,jat))**2+&
             (rxyz(3,iat)-rxyz(3,jat))**2 ==0.d0) then
           nateq=nateq+1
           write(*,'(1x,a,2(i0,a,a6,a))')'ERROR: atoms ',iat,&
                ' (',trim(at%atomnames(at%iatype(iat))),') and ',&
                jat,' (',trim(at%atomnames(at%iatype(jat))),&
                ') have the same positions'
        end if
     end do
  end do
  if (nateq /= 0) then
     write(*,'(1x,a)')'Control your posinp file, cannot proceed'
     stop
  end if

  if (iproc.eq.0) write(*,'(1x,a,i0)') 'Number of atom types= ',at%ntypes

  do ityp=1,at%ntypes
     if (iproc.eq.0) &
          write(*,'(1x,a,i0,a,a)') 'Atoms of type ',ityp,' are ',trim(at%atomnames(ityp))
  enddo

  do iat=1,at%nat
     if (iproc.eq.0 .and. at%lfrztyp(iat)) &
          write(*,'(1x,a,i0,a,a)') 'FIXED Atom N.:',iat,', Name: ',trim(at%atomnames(at%iatype(iat)))
  enddo

  !print 

end subroutine read_atomic_positions


! Fill the arrays occup and spinsgn
! if iunit /=0 this means that the file occup.dat does exist and it opens
subroutine input_occup(iproc,iunit,nelec,norb,norbu,norbd,nspin,mpol,occup,spinsgn)
  implicit none
! Arguments
  integer, intent(in) :: nelec,nspin,mpol,iproc,norb,norbu,norbd,iunit
  real(kind=8), intent(out) :: occup(norb),spinsgn(norb)
! Local variables
  integer :: iorb,nt,ne,it,ierror,iorb1
  real(kind=8) :: rocc,rup,rdown
  character(len=100) :: line


  do iorb=1,norb
     spinsgn(iorb)=1.0d0
  end do
  if (nspin/=1) then
     do iorb=1,norbu
        spinsgn(iorb)=1.0d0
     end do
     do iorb=norbu+1,norb
        spinsgn(iorb)=-1.0d0
     end do
  end if
! write(*,'(1x,a,5i4,30f6.2)')'Spins: ',norb,norbu,norbd,norbup,norbdp,(spinsgn(iorb),iorb=1,norb)

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
!!$        !Check if the polarisation is respected (mpol)
!!$        rup=sum(occup(1:norbu))
!!$        rdown=sum(occup(norbu+1:norb))
!!$        if (abs(rup-rdown-real(norbu-norbd,kind=8))>1.d-6) then
!!$           if (iproc==0) then
!!$              write(*,'(1x,a,f13.6,a,i0)') 'From the file "occup.dat", the polarization ',rup-rdown,&
!!$                             ' is not equal to ',norbu-norbd
!!$           end if
!!$           stop
!!$        end if
        !Fill spinsgn
        do iorb=1,norbu
           spinsgn(iorb)=1.0d0
        end do
        do iorb=norbu+1,norb
           spinsgn(iorb)=-1.0d0
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

!this routine performs also some cross-checks with other variables
subroutine read_system_variables(iproc,nproc,in,at,radii_cf,nelec,norb,norbu,norbd,norbp,iunit)
  use module_types
  implicit none
  type(input_variables), intent(in) :: in
  integer, intent(in) :: iproc,nproc
  type(atoms_data), intent(inout) :: at
  integer, intent(out) :: nelec,norb,norbu,norbd,norbp,iunit
  real(kind=8), dimension(at%ntypes,2), intent(out) :: radii_cf
  !local variables
  real(kind=8), parameter :: eps_mach=1.d-12
  logical :: exists
  character(len=2) :: symbol
  character(len=27) :: filename
  character(len=50) :: format
  integer :: i,j,k,l,iat,nlterms,nprl,nn,nt,ntu,ntd,ityp,ierror,i_stat,i_all,ixcpsp,ispinsum,mxpl
  real(kind=8) :: rcov,rprb,ehomo,radfine,tt,minrad
  real(kind=8), dimension(3,3) :: hij
  real(kind=8), dimension(2,2,3) :: offdiagarr
  integer, dimension(6,0:3) :: neleconf

  !allocate atoms data variables
  ! store PSP parameters
  ! modified to accept both GTH and HGHs pseudopotential types
  allocate(at%psppar(0:4,0:6,at%ntypes),stat=i_stat)
  call memocc(i_stat,product(shape(at%psppar))*kind(at%psppar),'psppar','read_system_variables')
  allocate(at%nelpsp(at%ntypes),stat=i_stat)
  call memocc(i_stat,product(shape(at%nelpsp))*kind(at%nelpsp),'nelpsp','read_system_variables')
  allocate(at%npspcode(at%ntypes),stat=i_stat)
  call memocc(i_stat,product(shape(at%npspcode))*kind(at%npspcode),'npspcode','read_system_variables')
  allocate(at%nzatom(at%ntypes),stat=i_stat)
  call memocc(i_stat,product(shape(at%nzatom))*kind(at%nzatom),'nzatom','read_system_variables')
  allocate(at%iasctype(at%ntypes),stat=i_stat)
  call memocc(i_stat,product(shape(at%iasctype))*kind(at%iasctype),'iasctype','read_system_variables')

  if (iproc == 0) then
     write(*,'(1x,a)')&
          'Atom Name   Ext.Electrons  PSP Code  Radii: Coarse     Fine   Calculated   From File'
  end if

  do ityp=1,at%ntypes
     filename = 'psppar.'//at%atomnames(ityp)

     inquire(file=filename,exist=exists)
     if (.not. exists) then
        if (iproc == 0) write(*,'(1x,a)')&
             'ERROR: The pseudopotential parameter file "'//filename//'" is lacking, exiting...'
        stop
     end if
     ! if (iproc.eq.0) write(*,*) 'opening PSP file ',filename
     open(unit=11,file=filename,status='old',iostat=ierror)
     !Check the open statement
     if (ierror /= 0) then
        write(*,*) 'iproc=',iproc,': Failed to open the file (it must be in ABINIT format!): "',&
             trim(filename),'"'
        stop
     end if
     read(11,*)
     read(11,*) at%nzatom(ityp),at%nelpsp(ityp)
     read(11,*) at%npspcode(ityp),ixcpsp
     !control if the PSP is calculated with the same XC value
     if (ixcpsp /= in%ixc .and. iproc==0) then
        write(*,'(1x,a)')        'WARNING: The pseudopotential file "'//trim(filename)//'"'
        write(*,'(1x,a,i0,a,i0)')'         contains a PSP generated with an XC id=',&
             ixcpsp,' while for this run ixc=',in%ixc
     end if
     at%psppar(:,:,ityp)=0.d0
     if (at%npspcode(ityp) == 2) then !GTH case
        read(11,*) (at%psppar(0,j,ityp),j=0,4)
        do i=1,2
           read(11,*) (at%psppar(i,j,ityp),j=0,3-i)
        enddo
     else if (at%npspcode(ityp) == 3) then !HGH case
        read(11,*) (at%psppar(0,j,ityp),j=0,4)
        read(11,*) (at%psppar(1,j,ityp),j=0,3)
        do i=2,4
           read(11,*) (at%psppar(i,j,ityp),j=0,3)
           read(11,*) !k coefficients, not used for the moment (no spin-orbit coupling)
        enddo
     else if (at%npspcode(ityp) == 10) then !HGH-K case
        read(11,*) at%psppar(0,0,ityp),nn,(at%psppar(0,j,ityp),j=1,nn) !local PSP parameters
        read(11,*) nlterms !number of channels of the pseudo
        prjloop: do l=1,nlterms
           read(11,*) at%psppar(l,0,ityp),nprl,at%psppar(l,1,ityp),&
                (at%psppar(l,j+2,ityp),j=2,nprl) !h_ij terms
           do i=2,nprl
              read(11,*) at%psppar(l,i,ityp),(at%psppar(l,i+j+1,ityp),j=i+1,nprl) !h_ij terms
           end do
           if (l==1) cycle
           do i=1,nprl
              read(11,*) !k coefficients, not used
           end do
        end do prjloop
     else
        if (iproc == 0) then
           write(*,'(1x,a,a)')trim(at%atomnames(ityp)),&
                'unrecognized pspcode: only GTH, HGH & HGH-K pseudos (ABINIT format)'
        end if
        stop
       end if
       !see whether the atom is semicore or not
       call eleconf(at%nzatom(ityp),at%nelpsp(ityp),symbol,rcov,rprb,ehomo,&
            neleconf,at%iasctype(ityp))
       !if you want no semicore input guess electrons, uncomment the following line
       !at%iasctype(ityp)=0

       !here we must check of the input guess polarisation
       !first of all, control if the values are compatible with the atom configuration
       !do this for all atoms belonging to a given type
       !control the maximum polarisation allowed: consider only non-closed shells
       do iat=1,at%nat
          if (at%iatype(iat) == ityp) then
             mxpl=0
             do l=0,3
                do i=1,6
                   if (neleconf(i,l) /= 0 .and. neleconf(i,l) /= 2*(2*l+1)) then
                      mxpl=mxpl+neleconf(i,l)
                   end if
                end do
             end do
             if (abs(at%nspinat(iat)) > mxpl) then
                if (iproc ==0) write(*,'(1x,a,i0,a,a,2(a,i0))')&
                     'ERROR: Input polarisation of atom No.',iat,&
                     ' (',trim(at%atomnames(ityp)),') must be <=',mxpl,&
                     ', while found ',at%nspinat(iat)
!                stop
             end if
          end if
       end do


       !old way of calculating the radii, requires modification of the PSP files
       read(11,*,iostat=ierror) radii_cf(ityp,1),radii_cf(ityp,2)
       if (ierror.eq.0) then
          if (iproc==0) write(*,'(3x,a6,13x,i3,5x,i3,10x,2(1x,f8.5),a)')&
               trim(at%atomnames(ityp)),at%nelpsp(ityp),at%npspcode(ityp),&
               radii_cf(ityp,1),radii_cf(ityp,2),&
               '                   X    '
       else
          !assigning the radii by calculating physical parameters
          radii_cf(ityp,1)=1.d0/sqrt(abs(2.d0*ehomo))
          radfine=100.d0
          do i=0,4
             if (at%psppar(i,0,ityp)/=0.d0) then
                radfine=min(radfine,at%psppar(i,0,ityp))
             end if
          end do
          radii_cf(ityp,2)=radfine
          if (iproc==0) write(*,'(3x,a6,13x,i3,5x,i3,10x,2(1x,f8.5),a)')&
               trim(at%atomnames(ityp)),at%nelpsp(ityp),at%npspcode(ityp),&
               radii_cf(ityp,1),radii_cf(ityp,2),&
               '       X                '
       end if
       close(11)
       !control the hardest gaussian
       minrad=1.d10
       do i=0,4
          if (at%psppar(i,0,ityp)/=0.d0) then
             minrad=min(minrad,at%psppar(i,0,ityp))
          end if
       end do
       !control whether the grid spacing is too high or not
       if (iproc == 0 .and. in%hgrid > 2.5d0*minrad) then
          write(*,'(1x,a)')&
               'WARNING: The grid spacing value may be too high to treat correctly the above pseudo.' 
          write(*,'(1x,a,f5.2,a)')&
               '         Results can be meaningless if hgrid is bigger than',2.5d0*minrad,&
               '. At your own risk!'
       end if

    enddo

    !print the pseudopotential matrices
    if (iproc == 0) then
       do l=1,3
          do i=1,2
             do j=i+1,3
                  offdiagarr(i,j-i,l)=0.d0
                if (l==1) then
                   if (i==1) then
                      if (j==2)   offdiagarr(i,j-i,l)=-0.5d0*sqrt(3.d0/5.d0)
                      if (j==3)   offdiagarr(i,j-i,l)=0.5d0*sqrt(5.d0/21.d0)
                   else
                        offdiagarr(i,j-i,l)=-0.5d0*sqrt(100.d0/63.d0)
                   end if
                else if (l==2) then
                   if (i==1) then
                      if (j==2)   offdiagarr(i,j-i,l)=-0.5d0*sqrt(5.d0/7.d0)
                      if (j==3)   offdiagarr(i,j-i,l)=1.d0/6.d0*sqrt(35.d0/11.d0)
                   else
                        offdiagarr(i,j-i,l)=-7.d0/3.d0*sqrt(1.d0/11.d0)
                   end if
                else if (l==3) then
                   if (i==1) then
                      if (j==2)   offdiagarr(i,j-i,l)=-0.5d0*sqrt(7.d0/9.d0)
                      if (j==3)   offdiagarr(i,j-i,l)=0.5d0*sqrt(63.d0/143.d0)
                   else
                        offdiagarr(i,j-i,l)=-9.d0*sqrt(1.d0/143.d0)
                   end if
                end if
             end do
          end do
       end do

       write(*,'(1x,a)')&
            '------------------------------------ Pseudopotential coefficients (Upper Triangular)'
        do ityp=1,at%ntypes
        write(*,'(1x,a)')&
          'Atom Name    rloc      C1        C2        C3        C4  '
           do l=0,3
              if (l==0) then
                 do i=4,0,-1
                    j=i
                    if (at%psppar(l,i,ityp) /= 0.d0) exit
                 end do
                 write(*,'(3x,a6,5(1x,f9.5))')&
                      trim(at%atomnames(ityp)),(at%psppar(l,i,ityp),i=0,j)
              else
                 do i=3,0,-1
                    j=i
                    if (at%psppar(l,i,ityp) /= 0.d0) exit
                 end do
                 if (j /=0) then
                    write(*,'(1x,a,i0,a)')&
                         '    l=',l-1,' '//'     rl        h1j       h2j       h3j '
                    hij=0.d0
                    do i=1,j
                       hij(i,i)=at%psppar(l,i,ityp)
                    end do
                    if (at%npspcode(ityp) == 3) then !traditional HGH convention
                       hij(1,2)=offdiagarr(1,1,l)*at%psppar(l,2,ityp)
                       hij(1,3)=offdiagarr(1,2,l)*at%psppar(l,3,ityp)
                       hij(2,3)=offdiagarr(2,1,l)*at%psppar(l,3,ityp)
                    else if (at%npspcode(ityp) == 10) then !HGH-K convention
                       hij(1,2)=at%psppar(l,4,ityp)
                       hij(1,3)=at%psppar(l,5,ityp)
                       hij(2,3)=at%psppar(l,6,ityp)
                    end if
                    do i=1,j
                       if (i==1) then
                          write(format,'(a,2(i0,a))')"(9x,(1x,f9.5),",j,"(1x,f9.5))"
                          write(*,format)at%psppar(l,0,ityp),(hij(i,k),k=i,j)
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
    at%natsc=0
    do iat=1,at%nat
       ityp=at%iatype(iat)
       nelec=nelec+at%nelpsp(ityp)
       if (at%iasctype(ityp) /= 0) at%natsc=at%natsc+1
    enddo
    nelec=nelec-in%ncharge
    if (iproc.eq.0) then
       write(*,'(1x,a,i8)') &
            'Total Number of Electrons ',nelec
    end if

    ! Number of orbitals
    if (in%nspin==1) then
       norb=(nelec+1)/2
       norbu=norb
       norbd=0
       if (mod(nelec,2).ne.0 .and. iproc==0) then
          write(*,'(1x,a)') 'WARNING: odd number of electrons, no closed shell system'
       end if
!    else if(in%nspin==4) then
!       if (iproc==0) write(*,'(1x,a)') 'Spin-polarized non-collinear calculation'
!       norb=nelec
!       norbu=norb
!       norbd=0
    else 
       if (iproc==0) write(*,'(1x,a)') 'Spin-polarized calculation'
       norb=nelec
       norbu=min(norb/2+in%mpol,norb)
       norbd=norb-norbu

       !test if the spin is compatible with the input guess polarisations
       ispinsum=0
       do iat=1,at%nat
          ispinsum=ispinsum+at%nspinat(iat)
       end do
       if (in%nspin == 2 .and. ispinsum /= norbu-norbd) then
          if (iproc==0) then 
             write(*,'(1x,a,i0,a)')&
                  'ERROR: Total input polarisation (found ',ispinsum,&
                  ') must be equal to norbu-norbd.'
             write(*,'(1x,3(a,i0))')&
                  'With norb=',norb,' and mpol=',in%mpol,' norbu-norbd=',norbu-norbd
             stop
          end if
       end if

       !now warn if there is no input guess spin polarisation
       ispinsum=0
       do iat=1,at%nat
          ispinsum=ispinsum+abs(at%nspinat(iat))
       end do
       if (ispinsum == 0 .and. in%nspin==2) then
          if (iproc==0) write(*,'(1x,a)')&
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
          if (iproc==0) write(*,'(1x,a)') &
               'ERROR: reading the number of orbitals in the file "occup.dat"'
          stop
       end if
       !Check
       if (in%nspin==1) then
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


subroutine system_size(iproc,geocode,atoms,rxyz,radii_cf,crmult,frmult,hx,hy,hz,&
     alat1,alat2,alat3,n1,n2,n3,nfl1,nfl2,nfl3,nfu1,nfu2,nfu3,n1i,n2i,n3i)
  !calculates the overall size of the simulation cell (cxmin,cxmax,cymin,cymax,czmin,czmax)
  !and shifts the atoms such that their position is the most symmetric possible
  use module_types
  implicit none
  type(atoms_data), intent(in) :: atoms
  character(len=1), intent(in) :: geocode
  integer, intent(in) :: iproc
  real(kind=8), intent(in) :: crmult,frmult
  real(kind=8), dimension(3,atoms%nat), intent(inout) :: rxyz
  real(kind=8), dimension(atoms%ntypes,2), intent(in) :: radii_cf
  integer, intent(out) :: n1,n2,n3,nfl1,nfl2,nfl3,nfu1,nfu2,nfu3,n1i,n2i,n3i
  real(kind=8), intent(inout) :: hx,hy,hz,alat1,alat2,alat3
  !local variables
  real(kind=8), parameter ::eps_mach=1.d-12,onem=1.d0-eps_mach
  integer :: iat,j
  real(kind=8) :: rad,cxmin,cxmax,cymin,cymax,czmin,czmax,alatrue1,alatrue2,alatrue3

  !check the geometry code with the grid spacings
  if (geocode == 'F' .and. (hx/=hy .or. hx/=hz .or. hy/=hz)) then
     write(*,'(1x,a)')'ERROR: The values of the grid spacings must be equal in the Free BC case'
     stop
  end if


  !calculate the extremes of the boxes taking into account the spheres arount the atoms
  cxmax=-1.d10 ; cxmin=1.d10
  cymax=-1.d10 ; cymin=1.d10
  czmax=-1.d10 ; czmin=1.d10
  do iat=1,atoms%nat
     rad=radii_cf(atoms%iatype(iat),1)*crmult
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


  !define the box sizes for free BC, and calculate dimensions for the fine grid with ISF
  if (geocode == 'F') then
     alat1=(cxmax-cxmin)
     alat2=(cymax-cymin)
     alat3=(czmax-czmin)

     ! grid sizes n1,n2,n3
     n1=int(alat1/hx)
     n2=int(alat2/hy)
     n3=int(alat3/hz)
     alatrue1=real(n1,kind=8)*hx
     alatrue2=real(n2,kind=8)*hy
     alatrue3=real(n3,kind=8)*hz

     n1i=2*n1+31
     n2i=2*n2+31
     n3i=2*n3+31

  else if (geocode == 'P') then !define the grid spacings, controlling the FFT compatibility
     call correct_grid(alat1,hx,n1)
     call correct_grid(alat2,hy,n2)
     call correct_grid(alat3,hz,n3)
     alatrue1=(cxmax-cxmin)
     alatrue2=(cymax-cymin)
     alatrue3=(czmax-czmin)

     n1i=2*n1+2
     n2i=2*n2+2
     n3i=2*n3+2

  else if (geocode == 'S') then
     call correct_grid(alat1,hx,n1)
     alat2=(cymax-cymin)
     call correct_grid(alat3,hz,n3)

     alatrue1=(cxmax-cxmin)
     alat2=(cymax-cymin)
     n2=int(alat2/hy)
     alatrue2=real(n2,kind=8)*hy
     alatrue3=(czmax-czmin)

     n1i=2*n1+2
     n2i=2*n2+31
     n3i=2*n3+2

  end if

  !balanced shift taking into account the missing space
  cxmin=cxmin+0.5d0*(alat1-alatrue1)
  cymin=cymin+0.5d0*(alat2-alatrue2)
  czmin=czmin+0.5d0*(alat3-alatrue3)

  !correct the box sizes for the isolated case
  if (geocode == 'F') then
     alat1=alatrue1
     alat2=alatrue2
     alat3=alatrue3
  else if (geocode == 'S') then
     alat2=alatrue2
  end if

  do iat=1,atoms%nat
     rxyz(1,iat)=rxyz(1,iat)-cxmin
     rxyz(2,iat)=rxyz(2,iat)-cymin
     rxyz(3,iat)=rxyz(3,iat)-czmin
  enddo

  ! fine grid size (needed for creation of input wavefunction, preconditioning)
  nfl1=n1 ; nfl2=n2 ; nfl3=n3
  nfu1=0 ; nfu2=0 ; nfu3=0
  do iat=1,atoms%nat
     rad=radii_cf(atoms%iatype(iat),2)*frmult
     nfl1=min(nfl1,ceiling((rxyz(1,iat)-rad)/hx - eps_mach))
     nfu1=max(nfu1,floor((rxyz(1,iat)+rad)/hx + eps_mach))

     nfl2=min(nfl2,ceiling((rxyz(2,iat)-rad)/hy - eps_mach))
     nfu2=max(nfu2,floor((rxyz(2,iat)+rad)/hy + eps_mach))

     nfl3=min(nfl3,ceiling((rxyz(3,iat)-rad)/hz - eps_mach)) 
     nfu3=max(nfu3,floor((rxyz(3,iat)+rad)/hz + eps_mach))
  enddo

  !correct the values of the delimiter if they go outside the box
  if (nfl1 < 0) nfl1=0
  if (nfl2 < 0) nfl2=0
  if (nfl3 < 0) nfl3=0

  if (nfu1 > n1) nfl1=n1
  if (nfu2 > n2) nfl2=n2
  if (nfu3 > n3) nfl3=n3


  if (iproc.eq.0) then
     write(*,'(1x,a,19x,a)') 'Shifted atomic positions, Atomic Units:','grid spacing units:'
     do iat=1,atoms%nat
        write(*,'(1x,i5,1x,a6,3(1x,1pe12.5),3x,3(1x,0pf9.3))') &
             iat,trim(atoms%atomnames(atoms%iatype(iat))),&
             (rxyz(j,iat),j=1,3),rxyz(1,iat)/hx,rxyz(2,iat)/hy,rxyz(3,iat)/hz
     enddo
     write(*,'(1x,a,3(1x,1pe12.5),a,3(1x,0pf5.2))') &
          '   Shift of=',-cxmin,-cymin,-czmin,' Grid Spacings=',hx,hy,hz
     write(*,'(1x,a,3(1x,1pe12.5),3x,3(1x,i9))')&
          '  Box Sizes=',alat1,alat2,alat3,n1,n2,n3
     write(*,'(1x,a,3x,3(3x,i4,a1,i0))')&
          '      Extremes for the high resolution grid points:',&
          nfl1,'<',nfu1,nfl2,'<',nfu2,nfl3,'<',nfu3
  endif

end subroutine system_size

subroutine correct_grid(a,h,n)
  use Poisson_Solver
  implicit none
  real(kind=8), intent(in) :: a
  integer, intent(inout) :: n
  real(kind=8), intent(inout) :: h
  !local variables
  integer :: m

  n=int(a/h)
  m=2*n+2
  do 
     call fourier_dim(m,m)
     if ((m/2)*2==m) then
        n=(m-2)/2
        exit
     else
        m=m+1
     end if
  end do
  h=a/real(n,kind=8)
  
end subroutine correct_grid
