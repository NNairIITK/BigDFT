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
  use module_base
  use module_types
  implicit none
  integer, intent(in) :: iproc
  type(input_variables), intent(out) :: in
  !local variables
  character(len=7) :: cudagpu
  character(len=100) :: line
  real(kind=4) :: hgrid,crmult,frmult,cpmult,fpmult
  integer :: ierror,ierrfrc,iconv,iblas

  ! Read the input variables.
  open(unit=1,file='input.dat',status='old')

  !read the line for force the CUDA GPU calculation for all processors
  read(1,'(a100)')line
  read(line,*,iostat=ierrfrc) cudagpu
  if (ierrfrc == 0 .and. cudagpu=='CUDAGPU') then
     call set_cpu_gpu_aff(iproc,iconv,iblas)
     if (iconv == 0) then
        !change the value of the GPU convolution flag defined in the module_base
        GPUconv=.true.
     end if
     if (iblas == 0) then
        !change the value of the GPU convolution flag defined in the module_base
        GPUblas=.true.
     end if
     read(1,*,iostat=ierror) in%ncount_cluster_x
  else
     read(line,*,iostat=ierror) in%ncount_cluster_x
  end if

  !read(1,*,iostat=ierror) in%ncount_cluster_x
  read(1,'(a100)')line
  read(line,*,iostat=ierrfrc) in%frac_fluct,in%forcemax
  if (ierrfrc /= 0) then
     read(line,*,iostat=ierror) in%frac_fluct
     in%forcemax=0.0_gp
  end if
  read(1,*,iostat=ierror) in%randdis
  read(1,*,iostat=ierror) in%betax
  read(1,*,iostat=ierror) hgrid
  read(1,*,iostat=ierror) crmult
  read(1,*,iostat=ierror) frmult
  read(1,*,iostat=ierror) cpmult
  read(1,*,iostat=ierror) fpmult
  !put the value at the max, such that to coincide with the maximum possible extension
  fpmult=100.e0
  in%hgrid  = real(hgrid,gp)
  in%crmult = real(crmult,gp)
  in%frmult = real(frmult,gp)
  in%cpmult = real(cpmult,gp)
  in%fpmult = real(fpmult,gp)
  in%cpmult = in%frmult
  if (in%fpmult > in%frmult) then
     if (iproc == 0) write(*,*) ' NONSENSE: fpmult > frmult, putting them equal'
     in%fpmult=in%frmult
  end if
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
  read(1,*,iostat=ierror) in%inputPsiId,in%output_wf,in%output_grid

  if (ierror/=0) then
     if (iproc == 0) write(*,'(1x,a)') 'Error while reading the file "input.dat"'
     stop
  end if

  !add reading lines for Davidson treatment (optional for backward compatibility)
  read(1,*,iostat=ierror) in%nvirt, in%nplot
  
  if (ierror/=0) then
     in%nvirt=0
     in%nplot=0
  else
     !performs some check: for the moment Davidson treatment is allowed only for spin-unpolarised
     !systems
     if (in%nspin/=1 .and. in%nvirt/=0) then
        if (iproc==0) then
           write(*,'(1x,a)')'ERROR: Davidson treeatment allowed on fon non spin-polarised systems'
        end if
        stop
     end if
  end if
 
  close(1,iostat=ierror)

  if (iproc == 0) then
     write(*,'(1x,a,i0)') 'Max. number of wavefnctn optim ',in%ncount_cluster_x
     write(*,'(1x,a,1pe10.2)') 'Convergence criterion for forces: fraction of noise ',&
          in%frac_fluct
     write(*,'(1x,a,1pe10.2)') '                                : maximal component ',&
          in%forcemax
     write(*,'(1x,a,1pe10.2)') 'Random displacement amplitude ',in%randdis
     write(*,'(1x,a,1pe10.2)') 'Steepest descent step ',in%betax
     if (in%nvirt > 0) then
        !read virtual orbital and plotting request
        write(*,'(1x,a,i0)')'Virtual orbitals ',in%nvirt
        write(*,'(1x,a,i0,a)')'Output for density plots is requested for ',in%nplot,' orbitals'
     end if
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
  if (in%geocode /= 'F') then
     write(*,'(1x,a,1x,a,3(1x,1pe12.5))')&
          '  Geom. Code=    '//in%geocode//'   |',&
          '  Box Sizes (Bohr) =',in%alat1,in%alat2,in%alat3

  end if
end subroutine print_input_parameters

! read atomic positions
subroutine read_atomic_positions(iproc,ifile,units,in,at,rxyz)
  use module_base
  use module_types
  implicit none
  character(len=20), intent(in) :: units
  integer, intent(in) :: iproc,ifile
  type(atoms_data), intent(inout) :: at
  type(input_variables), intent(out) :: in
  real(gp), dimension(3,at%nat), intent(out) :: rxyz
  !local variables
  character(len=*), parameter :: subname='read_atomic_positions'
  real(gp), parameter :: bohr=0.5291772108_gp !1 AU in angstroem
  character(len=3) :: suffix
  character(len=2) :: symbol
  character(len=20) :: tatonam
  character(len=50) :: extra
  character(len=100) :: line
  logical :: lpsdbl 
  integer :: nateq,iat,jat,ityp,i,ierror,ierrsfx,i_stat,natpol
! To read the file posinp (avoid differences between compilers)
  real(kind=4) :: rx,ry,rz,alat1,alat2,alat3
! case for which the atomic positions are given whithin general precision
  real(gp) :: rxd0,ryd0,rzd0,alat1d0,alat2d0,alat3d0
  character(len=20), dimension(100) :: atomnames

  if (iproc.eq.0) write(*,'(1x,a,i0)') 'Number of atoms     = ',at%nat

  allocate(at%iatype(at%nat+ndebug),stat=i_stat)
  call memocc(i_stat,at%iatype,'at%iatype',subname)
  allocate(at%lfrztyp(at%nat+ndebug),stat=i_stat)
  call memocc(i_stat,at%lfrztyp,'at%lfrztyp',subname)
  allocate(at%natpol(at%nat+ndebug),stat=i_stat)
  call memocc(i_stat,at%natpol,'at%natpol',subname)

  !controls if the positions are provided with machine precision
  if (units == 'angstroemd0' .or. units== 'atomicd0' .or. units== 'bohrd0') then
     lpsdbl=.true.
  else
     lpsdbl=.false.
  end if

  !this array is useful for frozen atoms
  !no atom is frozen by default
  at%lfrztyp(:)=.false.
  !also the spin polarisation and the charge are is fixed to zero by default
  !this corresponds to the value of 100
  !RULE natpol=charge*1000 + 100 + spinpol
  at%natpol(:)=100

  !read from positions of .xyz format, but accepts also the old .ascii format
  read(ifile,'(a100)')line

  !old format, still here for backward compatibility
  !admits only simple precision calculation
  read(line,*,iostat=ierror) rx,ry,rz,tatonam

  !in case of old format, put geocode to F and alat to 0.
  if (ierror == 0) then
     in%geocode='F'
     in%alat1=0.0_gp
     in%alat2=0.0_gp
     in%alat3=0.0_gp
  else
     if (lpsdbl) then
        read(line,*,iostat=ierrsfx) tatonam,alat1d0,alat2d0,alat3d0
     else
        read(line,*,iostat=ierrsfx) tatonam,alat1,alat2,alat3
     end if
     if (ierrsfx == 0) then
        if (trim(tatonam)=='periodic') then
           in%geocode='P'
        else if (trim(tatonam)=='surface') then 
           in%geocode='S'
        else !otherwise free bc
           in%geocode='F'
           in%alat1=0.0_gp
           in%alat2=0.0_gp
           in%alat3=0.0_gp
        end if
        if (.not. lpsdbl) then
           alat1d0=real(alat1,gp)
           alat2d0=real(alat2,gp)
           alat3d0=real(alat3,gp)
        end if
     else
        in%geocode='F'
        in%alat1=0.0_gp
        in%alat2=0.0_gp
        in%alat3=0.0_gp
     end if
  end if

  at%ntypes=0
  do iat=1,at%nat
     if (ierror == 0) then
        !old case of ascii file, added for backward compatibility
        if (iat /= 1) read(ifile,*) rx,ry,rz,tatonam
     else
        !xyz input file, allow extra information
        read(ifile,'(a100)')line 
        if (lpsdbl) then
           read(line,*,iostat=ierrsfx)symbol,rxd0,ryd0,rzd0,extra
        else
           read(line,*,iostat=ierrsfx)symbol,rx,ry,rz,extra
        end if
        call find_extra_info(line,extra)
        call parse_extra_info(iproc,iat,extra,at)
!!$        if (ierrsfx ==0) then
!!$           call parse_extra_info(iproc,iat,extra,at)
!!$        else
!!$           if (lpsdbl) then
!!$              read(line,*)symbol,rx,ry,rz
!!$           else
!!$              read(line,*)symbol,rxd0,ryd0,rzd0
!!$           end if
!!$        end if
        
!!$        if (lpsdbl) then
!!$           read(line,*,iostat=ierrsfx)symbol,rxd0,ryd0,rzd0,natpol,suffix
!!$        else
!!$           read(line,*,iostat=ierrsfx)symbol,rx,ry,rz,natpol,suffix
!!$        end if
!!$        if (ierrsfx ==0) then
!!$           at%natpol(iat)=natpol
!!$           !tatonam=trim(symbol)//'_'//trim(suffix)
!!$           if (suffix == 'f') then
!!$              !the atom is considered as blocked
!!$              at%lfrztyp(iat)=.true.
!!$           else
!!$              if (iproc == 0) then
!!$                 print *,suffix
!!$                 write(*,'(1x,a,i0,a)')'ERROR in input file for atom number ',&
!!$                      iat,': the only value accepted in 5th column is "f"'
!!$                 stop
!!$              end if
!!$           end if
!!$        else
!!$           if (lpsdbl) then
!!$              read(line,*,iostat=ierrsfx)symbol,rxd0,ryd0,rzd0,suffix
!!$           else
!!$              read(line,*,iostat=ierrsfx)symbol,rx,ry,rz,suffix
!!$           end if
!!$           if (ierrsfx ==0) then
!!$              if (suffix == 'f') then
!!$                 !the atom is considered as blocked
!!$                 at%lfrztyp(iat)=.true.
!!$              else
!!$                 read(suffix,*,iostat=i_stat)at%natpol(iat)
!!$                 if(i_stat /=0) then
!!$                    if (iproc == 0) then
!!$                       print *,suffix
!!$                       write(*,'(1x,a,i0,a)')'ERROR in input file for atom number ',&
!!$                            iat,': in 4th column you can put the input polarisation or "f"'
!!$                       stop
!!$                    end if
!!$                 end if
!!$              end if
!!$           else
!!$              if (lpsdbl) then
!!$                 read(line,*)symbol,rx,ry,rz
!!$              else
!!$                 read(line,*)symbol,rxd0,ryd0,rzd0
!!$              end if
!!$           end if
!!$        end if
        tatonam=trim(symbol)
     end if
     if (lpsdbl) then
        rxyz(1,iat)=rxd0
        rxyz(2,iat)=ryd0
        rxyz(3,iat)=rzd0
     else
        rxyz(1,iat)=real(rx,gp)
        rxyz(2,iat)=real(ry,gp)
        rxyz(3,iat)=real(rz,gp)
     end if

     if (in%geocode == 'P') then
        rxyz(1,iat)=modulo(rxyz(1,iat),alat1d0)
        rxyz(2,iat)=modulo(rxyz(2,iat),alat2d0)
        rxyz(3,iat)=modulo(rxyz(3,iat),alat3d0)
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
        in%alat1=alat1d0/bohr
        in%alat2=alat2d0/bohr
        in%alat3=alat3d0/bohr
        do i=1,3 
           rxyz(i,iat)=rxyz(i,iat)/bohr
        enddo
     else if  (units=='atomic' .or. units=='bohr'  .or.&
          units== 'atomicd0' .or. units== 'bohrd0') then
        in%alat1=alat1d0
        in%alat2=alat2d0
        in%alat3=alat3d0
     else
        write(*,*) 'length units in input file unrecognized'
        write(*,*) 'recognized units are angstroem or atomic = bohr'
        stop 
     endif
  enddo

  !now that ntypes is determined allocate at%atomnames and copy the values
  allocate(at%atomnames(at%ntypes+ndebug),stat=i_stat)
  call memocc(i_stat,at%atomnames,'at%atomnames',subname)
  at%atomnames(1:at%ntypes)=atomnames(1:at%ntypes)

  !control atom positions
  nateq=0
  do iat=1,at%nat
     do jat=iat+1,at%nat
        if ((rxyz(1,iat)-rxyz(1,jat))**2+(rxyz(2,iat)-rxyz(2,jat))**2+&
             (rxyz(3,iat)-rxyz(3,jat))**2 ==0.0_gp) then
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

end subroutine read_atomic_positions

subroutine find_extra_info(line,extra)
  implicit none
  character(len=100), intent(in) :: line
  character(len=50), intent(out) :: extra
  !local variables
  logical :: space
  integer :: i,nspace
  i=1
  space=.true.
  nspace=-1
  find_space : do
     !toggle the space value for each time
     if (line(i:i) == ' ' .neqv. space) then
        nspace=nspace+1
        space=.not. space
     end if
     !print *,line(i:i),nspace
     if (nspace==7) then
        extra=line(i:min(100,i+49))
        exit find_space
     end if
     i=i+1
  end do find_space

end subroutine find_extra_info

subroutine parse_extra_info(iproc,iat,extra,at)
  use module_types
  implicit none
  character(len=50), intent(in) :: extra
  integer, intent(in) :: iat,iproc
  type(atoms_data), intent(out) :: at
  !local variables
  character(len=3) :: suffix
  integer :: ierr,ierr1,ierr2,nspol,nchrg,nsgn
  !case with all the information
  read(extra,*,iostat=ierr)nspol,nchrg,suffix
  if (extra == '') then !case with empty information
     nspol=0
     nchrg=0
     suffix='   '
  else if (ierr /= 0) then !case with partial information
     read(extra,*,iostat=ierr1)nspol,suffix
     if (ierr1 /=0) then
        if (trim(extra) == 'f') then
           suffix=trim(extra)
           nspol=0
           nchrg=0
        else
           read(extra,*,iostat=ierr2)nspol
           if (ierr2 /=0) then
              call error
           end if
           suffix='   '
           nchrg=0
        end if
     else
        nchrg=0
        if (trim(suffix) /= 'f') then
           read(suffix,*,iostat=ierr2)nchrg
           if (ierr2 /= 0) then
              call error
           else
              suffix='   '
           end if
        end if
     end if
  end if

  !now assign the array, following the rule
  if(nchrg>=0) then
     nsgn=1
  else
     nsgn=-1
  end if
  at%natpol(iat)=1000*nchrg+nsgn*100+nspol

  !print *,'natpol atomic',iat,at%natpol(iat)
  
  if (trim(suffix) == 'f') then
     !the atom is considered as blocked
     at%lfrztyp(iat)=.true.
  end if

contains

 subroutine error
   if (iproc == 0) then
      print *,extra
      write(*,'(1x,a,i0,a)')&
           'ERROR in input file for atom number ',iat,&
           ': after 4th column you can put the input polarisation(s) or "f"'
   end if
   stop
 end subroutine error
  
end subroutine parse_extra_info


! Fill the arrays occup and spinsgn
! if iunit /=0 this means that the file occup.dat does exist and it opens
subroutine input_occup(iproc,iunit,nelec,norb,norbu,norbd,nspin,mpol,occup,spinsgn)
  use module_base
  implicit none
! Arguments
  integer, intent(in) :: nelec,nspin,mpol,iproc,norb,norbu,norbd,iunit
  real(gp), intent(out) :: occup(norb),spinsgn(norb)
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
              if (iproc==0) then
                 write(*,'(1x,a,i0,a)') 'ERROR in line ',nt+1,' of the file "occup.dat"'
                 write(*,'(10x,a,i0,a)') 'The orbital index ',iorb,' is incorrect'
              end if
              stop
           elseif (rocc<0._gp .or. rocc>2._gp) then
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
     if (abs(rocc-real(nelec,gp))>1.e-6_gp) then
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

!calculate the charge and the spin polarisation to be placed on a given atom
!RULE: natpol = c*1000 + sgn(c)*100 + s: charged and polarised atom (charge c, polarisation s)
subroutine charge_and_spol(natpol,nchrg,nspol)
  implicit none
  integer, intent(in) :: natpol
  integer, intent(out) :: nchrg,nspol
  !local variables
  integer :: nsgn

  nchrg=natpol/1000
  if (nchrg>=0) then
     nsgn=1
  else
     nsgn=-1
  end if

  nspol=natpol-1000*nchrg-nsgn*100

end subroutine charge_and_spol

!this routine performs also some cross-checks with other variables
subroutine read_system_variables(iproc,nproc,in,at,radii_cf,nelec,norb,norbu,norbd,norbp,iunit)
  use module_base
  use module_types
  implicit none
  type(input_variables), intent(in) :: in
  integer, intent(in) :: iproc,nproc
  type(atoms_data), intent(inout) :: at
  integer, intent(out) :: nelec,norb,norbu,norbd,norbp,iunit
  real(gp), dimension(at%ntypes,3), intent(out) :: radii_cf
  !local variables
  character(len=*), parameter :: subname='read_system_variables'
  real(kind=8), parameter :: eps_mach=1.d-12
  logical :: exists
  character(len=2) :: symbol
  character(len=24) :: message
  character(len=27) :: filename
  character(len=50) :: format
  character(len=100) :: line
  integer :: i,j,k,l,iat,nlterms,nprl,nn,nt,ntu,ntd,ityp,ierror,i_stat,i_all,ixcpsp,ispinsum,mxpl
  integer :: ispol,mxchg,ichg,natpol,ichgsum,nsccode,ierror1
  real(gp) :: rcov,rprb,ehomo,radfine,tt,minrad,maxrad
  real(gp), dimension(3,3) :: hij
  real(gp), dimension(2,2,3) :: offdiagarr
  integer, dimension(6,0:3) :: neleconf

  !allocate atoms data variables
  ! store PSP parameters
  ! modified to accept both GTH and HGHs pseudopotential types
  allocate(at%psppar(0:4,0:6,at%ntypes+ndebug),stat=i_stat)
  call memocc(i_stat,at%psppar,'at%psppar',subname)
  allocate(at%nelpsp(at%ntypes+ndebug),stat=i_stat)
  call memocc(i_stat,at%nelpsp,'at%nelpsp',subname)
  allocate(at%npspcode(at%ntypes+ndebug),stat=i_stat)
  call memocc(i_stat,at%npspcode,'at%npspcode',subname)
  allocate(at%nzatom(at%ntypes+ndebug),stat=i_stat)
  call memocc(i_stat,at%nzatom,'at%nzatom',subname)
  allocate(at%iasctype(at%ntypes+ndebug),stat=i_stat)
  call memocc(i_stat,at%iasctype,'at%iasctype',subname)

  if (iproc == 0) then
     write(*,'(1x,a)')&
          ' Atom    N.Electr.  PSP Code  Radii: Coarse     Fine  CoarsePSP    Calculated   File'
  end if

  do ityp=1,at%ntypes
     filename = 'psppar.'//at%atomnames(ityp)

     inquire(file=filename,exist=exists)
     if (.not. exists) then
        if (iproc == 0) write(*,'(1x,3a)')&
             'ERROR: The pseudopotential parameter file "',trim(filename),'" is lacking, exiting...'
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
     at%psppar(:,:,ityp)=0._gp
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
          neleconf,at%iasctype(ityp),mxpl,mxchg)
     !if you want no semicore input guess electrons, uncomment the following line
     !at%iasctype(ityp)=0

     !here we must check of the input guess polarisation
     !control if the values are compatible with the atom configuration
     !do this for all atoms belonging to a given type
     !control the maximum polarisation allowed: consider only non-closed shells   

     do iat=1,at%nat
        if (at%iatype(iat) == ityp) then
           call charge_and_spol(at%natpol(iat),ichg,ispol)
           if (abs(ispol) > mxpl) then
              if (iproc ==0) write(*,'(1x,a,i0,a,a,2(a,i0))')&
                   'ERROR: Input polarisation of atom No.',iat,&
                   ' (',trim(at%atomnames(ityp)),') must be <=',mxpl,&
                   ', while found ',ispol
              stop
           end if
           if (abs(ichg) > mxchg) then
              if (iproc ==0) write(*,'(1x,a,i0,a,a,2(a,i0))')&
                   'ERROR: Input charge of atom No.',iat,&
                   ' (',trim(at%atomnames(ityp)),') must be <=',mxchg,&
                   ', while found ',ichg
              stop
           end if
        end if
     end do

     !control the hardest and the softest gaussian
     minrad=1.e10_gp
     do i=0,4
        !the maximum radii is useful only for projectors
        if (i==1) maxrad=0.0_gp
        if (at%psppar(i,0,ityp)/=0._gp) then
           minrad=min(minrad,at%psppar(i,0,ityp))
           maxrad=max(maxrad,at%psppar(i,0,ityp))
        end if
     end do

     !old way of calculating the radii, requires modification of the PSP files
     read(11,'(a100)',iostat=ierror)line
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
           if (at%psppar(i,0,ityp)/=0._gp) then
              radfine=min(radfine,at%psppar(i,0,ityp))
           end if
        end do
        radii_cf(ityp,2)=radfine
        message='         X              '
        radii_cf(ityp,3)=radfine
     end if
     close(11)

     !correct the coarse and the fine radius for projectors
     radii_cf(ityp,3)=min(in%crmult*radii_cf(ityp,1),15.0_gp*maxrad)/in%frmult

     if (iproc==0) write(*,'(1x,a6,8x,i3,5x,i3,10x,3(1x,f8.5),a)')&
          trim(at%atomnames(ityp)),at%nelpsp(ityp),at%npspcode(ityp),&
          radii_cf(ityp,1),radii_cf(ityp,2),radii_cf(ityp,3),message

     !control whether the grid spacing is too high
     if (iproc == 0 .and. in%hgrid > 2.5_gp*minrad) then
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
     do ityp=1,at%ntypes
        write(*,'(1x,a)')&
             'Atom Name    rloc      C1        C2        C3        C4  '
        do l=0,3
           if (l==0) then
              do i=4,0,-1
                 j=i
                 if (at%psppar(l,i,ityp) /= 0._gp) exit
              end do
              write(*,'(3x,a6,5(1x,f9.5))')&
                   trim(at%atomnames(ityp)),(at%psppar(l,i,ityp),i=0,j)
           else
              do i=3,0,-1
                 j=i
                 if (at%psppar(l,i,ityp) /= 0._gp) exit
              end do
              if (j /=0) then
                 write(*,'(1x,a,i0,a)')&
                      '    l=',l-1,' '//'     rl        h1j       h2j       h3j '
                 hij=0._gp
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
     nsccode=at%iasctype(ityp)
     call charge_and_spol(at%natpol(iat),ichg,ispol)
     if (ichg /=0) then
        call eleconf(at%nzatom(ityp),at%nelpsp(ityp),symbol,rcov,rprb,ehomo,&
             neleconf,at%iasctype(ityp),mxpl,mxchg)
        call correct_semicore(at%atomnames(ityp),6,3,ichg,neleconf,nsccode)
     end if
     if (nsccode/= 0) at%natsc=at%natsc+1
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
     ichgsum=0
     do iat=1,at%nat
        call charge_and_spol(at%natpol(iat),ichg,ispol)
        ispinsum=ispinsum+ispol
        ichgsum=ichgsum+ichg
     end do

     if (in%nspin == 2 .and. ispinsum /= norbu-norbd) then
        if (iproc==0) then 
           write(*,'(1x,a,i0,a)')&
                'ERROR: Total input polarisation (found ',ispinsum,&
                ') must be equal to norbu-norbd.'
           write(*,'(1x,3(a,i0))')&
                'With norb=',norb,' and mpol=',in%mpol,' norbu-norbd=',norbu-norbd
        end if
        stop
     end if

     if (ichgsum /= in%ncharge .and. ichgsum /= 0) then
        if (iproc==0) then 
           write(*,'(1x,a,i0,a)')&
                'ERROR: Total input charge (found ',ichgsum,&
                ') cannot be different than charge.'
           write(*,'(1x,2(a,i0))')&
                'The charge is=',in%ncharge,' input charge=',ichgsum
        end if
        stop
     end if

     !now warn if there is no input guess spin polarisation
     ispinsum=0
     do iat=1,at%nat
        call charge_and_spol(at%natpol(iat),ichg,ispol)
        ispinsum=ispinsum+abs(ispol)
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
  use module_base
  use module_types
  implicit none
  type(atoms_data), intent(in) :: atoms
  character(len=1), intent(in) :: geocode
  integer, intent(in) :: iproc
  real(gp), intent(in) :: crmult,frmult
  real(gp), dimension(3,atoms%nat), intent(inout) :: rxyz
  real(gp), dimension(atoms%ntypes,3), intent(in) :: radii_cf
  integer, intent(out) :: n1,n2,n3,nfl1,nfl2,nfl3,nfu1,nfu2,nfu3,n1i,n2i,n3i
  real(gp), intent(inout) :: hx,hy,hz,alat1,alat2,alat3
  !local variables
  real(gp), parameter ::eps_mach=1.e-12_gp,onem=1.0_gp-eps_mach
  integer :: iat,j
  real(gp) :: rad,cxmin,cxmax,cymin,cymax,czmin,czmax,alatrue1,alatrue2,alatrue3

  !check the geometry code with the grid spacings
  if (geocode == 'F' .and. (hx/=hy .or. hx/=hz .or. hy/=hz)) then
     write(*,'(1x,a)')'ERROR: The values of the grid spacings must be equal in the Free BC case'
     stop
  end if


  !calculate the extremes of the boxes taking into account the spheres around the atoms
  cxmax=-1.e10_gp 
  cxmin=1.e10_gp

  cymax=-1.e10_gp 
  cymin=1.e10_gp

  czmax=-1.e10_gp 
  czmin=1.e10_gp

  do iat=1,atoms%nat

     rad=radii_cf(atoms%iatype(iat),1)*crmult

     cxmax=max(cxmax,rxyz(1,iat)+rad) 
     cxmin=min(cxmin,rxyz(1,iat)-rad)

     cymax=max(cymax,rxyz(2,iat)+rad) 
     cymin=min(cymin,rxyz(2,iat)-rad)
     
     czmax=max(czmax,rxyz(3,iat)+rad) 
     czmin=min(czmin,rxyz(3,iat)-rad)
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
     alatrue1=real(n1,gp)*hx
     alatrue2=real(n2,gp)*hy
     alatrue3=real(n3,gp)*hz

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
     alatrue2=real(n2,gp)*hy
     alatrue3=(czmax-czmin)

     n1i=2*n1+2
     n2i=2*n2+31
     n3i=2*n3+2

  end if

  !balanced shift taking into account the missing space
  cxmin=cxmin+0.5_gp*(alat1-alatrue1)
  cymin=cymin+0.5_gp*(alat2-alatrue2)
  czmin=czmin+0.5_gp*(alat3-alatrue3)

  !correct the box sizes for the isolated case
  if (geocode == 'F') then
     alat1=alatrue1
     alat2=alatrue2
     alat3=alatrue3
  else if (geocode == 'S') then
     alat2=alatrue2
  else if (geocode == 'P') then
     !for the moment we do not put the shift, at the end it will be tested
     !here we should put the center of mass
     cxmin=0.0_gp
     cymin=0.0_gp
     czmin=0.0_gp
  end if

  !if (geocode /= 'P') then
     do iat=1,atoms%nat
        rxyz(1,iat)=rxyz(1,iat)-cxmin
        rxyz(2,iat)=rxyz(2,iat)-cymin
        rxyz(3,iat)=rxyz(3,iat)-czmin
     enddo
  !else !place the atoms inside the box
  !   do iat=1,atoms%nat
  !      rxyz(1,iat)=modulo(rxyz(1,iat),alat1)
  !      rxyz(2,iat)=modulo(rxyz(2,iat),alat2)
  !      rxyz(3,iat)=modulo(rxyz(3,iat),alat3)
  !   enddo
  !end if

  ! fine grid size (needed for creation of input wavefunction, preconditioning)
  nfl1=n1 
  nfl2=n2 
  nfl3=n3

  nfu1=0 
  nfu2=0 
  nfu3=0

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
  if (nfl1 < 0 .or. nfu1 > n1) then
     nfl1=0
     nfu1=n1
  end if
  if (nfl2 < 0 .or. nfu2 > n2) then
     nfl2=0
     nfu2=n2
  end if
  if (nfl3 < 0 .or. nfu3 > n3) then
     nfl3=0
     nfu3=n3
  end if

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
  use module_base
  use Poisson_Solver
  implicit none
  real(gp), intent(in) :: a
  integer, intent(inout) :: n
  real(gp), intent(inout) :: h
  !local variables
  integer :: m

  !here the dimensions should be corrected in order to 
  !allow the fft for the preconditioner

  n=int(a/h)-1
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
  h=a/real(n+1,gp)
  
end subroutine correct_grid
