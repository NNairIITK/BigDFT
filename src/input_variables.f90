!!****f* BigDFT/print_logo
!! FUNCTION
!!    Display the logo of BigDFT 
!! COPYRIGHT
!!    Copyright (C) 2007-2008 BigDFT group 
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 
!! SOURCE
!!
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
  write(*,'(23x,a)')'    DDDDD       F         TTTT                     (Ver 1.2.2)'
  write(*,'(1x,a)')&
       '------------------------------------------------------------------------------------'
  write(*,'(1x,a)')&
       '|              Daubechies Wavelets for DFT Pseudopotential Calculations            |'
  write(*,'(1x,a)')&
       '------------------------------------------------------------------------------------'
end subroutine print_logo
!!***


!!****f* BigDFT/read_input_variables
!! FUNCTION
!!    Read the input variables in the file 'input.dat'
!! SOURCE
!!
subroutine read_input_variables(iproc,filename,in)
  use module_base
  use module_types
  implicit none
  character(len=*), intent(in) :: filename
  integer, intent(in) :: iproc
  type(input_variables), intent(out) :: in
  !local variables
  character(len=7) :: cudagpu
  character(len=100) :: line
  real(kind=4) :: hgrid,crmult,frmult,cpmult,fpmult
  integer :: ierror,ierrfrc,iconv,iblas,iline

  ! Read the input variables.
  open(unit=1,file=filename,status='old')

  iline=0
  !read the line for force the CUDA GPU calculation for all processors
  read(1,'(a100)')line
  read(line,*,iostat=ierrfrc) cudagpu
  if (ierrfrc == 0 .and. cudagpu=='CUDAGPU') then
    ! call set_cpu_gpu_aff(iproc,iconv,iblas)
     !GPUshare=.false.
     call init_gpu_sharing(8,2)
     GPUshare=.true.
     if (iconv == 0) then
        !change the value of the GPU convolution flag defined in the module_base
        GPUconv=.true.
     end if
     if (iblas == 0) then
        !change the value of the GPU convolution flag defined in the module_base
        GPUblas=.true.
     end if
     read(1,*,iostat=ierror) in%ncount_cluster_x
     call check()
  else
     read(line,*,iostat=ierror) in%ncount_cluster_x
     call check()
  end if

  read(1,'(a100)')line
  read(line,*,iostat=ierrfrc) in%frac_fluct,in%forcemax
  if (ierrfrc /= 0) then
     read(line,*,iostat=ierror) in%frac_fluct
     in%forcemax=0.0_gp
  end if
  call check()
  read(1,*,iostat=ierror) in%randdis
  call check()
  read(1,*,iostat=ierror) in%betax
  call check()
  read(1,*,iostat=ierror) hgrid
  call check()
  read(1,*,iostat=ierror) crmult
  call check()
  read(1,*,iostat=ierror) frmult
  call check()
  !read(1,*,iostat=ierror) cpmult !this value can be removed from the input files
  !read(1,*,iostat=ierror) fpmult !this value can be removed from the input files
  !put the value at the max, such that to coincide with the maximum possible extension
  in%hgrid  = real(hgrid,gp)
  in%crmult = real(crmult,gp)
  in%frmult = real(frmult,gp)

  !in%cpmult = in%frmult
  !in%fpmult=in%frmult

  read(1,*,iostat=ierror) in%ixc
  call check()
  read(1,'(a100)')line
  read(line,*,iostat=ierror) in%ncharge,in%ef(1)
  if (ierror == 0 .and. in%ef(1) /= 0.0_gp) then
     read(line,*,iostat=ierror) in%ncharge,in%ef(1),in%ef(2),in%ef(3)
  else
     in%ef(2)=0.0_gp
     in%ef(3)=0.0_gp
  end if
  call check()
  read(1,*,iostat=ierror) in%gnrm_cv
  call check()
  read(1,'(a100)')line
  read(line,*,iostat=ierror) in%itermax,in%nrepmax
  if (ierror == 0) then
     !read(line,*,iostat=ierror) in%ncharge,in%ef(1),in%ef(2),in%ef(3)
  else
     read(line,*,iostat=ierror)in%itermax
     in%nrepmax=10
  end if
  call check()
  read(1,*,iostat=ierror) in%ncong
  call check()
  read(1,*,iostat=ierror) in%idsx
  call check()
  read(1,*,iostat=ierror) in%calc_tail
  call check()
  read(1,*,iostat=ierror) in%rbuf
  call check()
  read(1,*,iostat=ierror) in%ncongt
  call check()
  read(1,*,iostat=ierror) in%nspin,in%mpol
  call check()
  read(1,*,iostat=ierror) in%inputPsiId,in%output_wf,in%output_grid
  call check()

  !project however the wavefunction on gaussians if asking to write them on disk
  in%gaussian_help=(in%inputPsiId >= 10) .or. in%output_wf 
  !switch on the gaussian auxiliary treatment 
  !and the zero of the forces
  if (in%inputPsiId == 10) then
     in%inputPsiId=0
  end if

  !add reading lines for Davidson treatment (optional for backward compatibility)
  read(1,*,iostat=ierror) in%nvirt, in%nplot
  !call check()
  
  if (ierror/=0) then
     in%nvirt=0
     in%nplot=0
  else
     !performs some check: for the moment Davidson treatment is allowed only for spin-unpolarised
     !systems
     if (in%nspin/=1 .and. in%nvirt/=0) then
        if (iproc==0) then
           write(*,'(1x,a)')'ERROR: Davidson treatment allowed only for non spin-polarised systems'
        end if
        stop
     end if
  end if
 
  close(unit=1,iostat=ierror)

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

contains

  subroutine check()
    iline=iline+1
    if (ierror/=0) then
       if (iproc == 0) write(*,'(1x,a,a,a,i3)') &
            'Error while reading the file "',trim(filename),'", line=',iline
       stop
    end if
  end subroutine check

end subroutine read_input_variables
!!***


!!****f* BigDFT/print_input_parameters
!! FUNCTION
!!    Print all input parameters
!! SOURCE
!!
subroutine print_input_parameters(in,atoms)
  use module_types
  implicit none
  type(input_variables), intent(in) :: in
  type(atoms_data), intent(in) :: atoms

  write(*,'(1x,a)')&
       '------------------------------------------------------------------- Input Parameters'
  write(*,'(1x,a)')&
       '    System Choice       Resolution Radii        SCF Iteration      Finite Size Corr.'
  write(*,'(1x,a,f7.3,1x,a,f5.2,1x,a,1pe8.1,1x,a,l4)')&
       'Grid spacing=',in%hgrid,   '|  Coarse Wfs.=',in%crmult,'| Wavefns Conv.=',in%gnrm_cv,&
       '| Calculate=',in%calc_tail
  write(*,'(1x,a,i7,1x,a,f5.2,1x,a,i5,a,i2,1x,a,f4.1)')&
       '       XC id=',in%ixc,     '|    Fine Wfs.=',in%frmult,'| Max. N. Iter.=',in%itermax,&
       'x',in%nrepmax,'| Extension=',in%rbuf
  write(*,'(1x,a,i7,1x,a,1x,a,i8,1x,a,i4)')&
       'total charge=',in%ncharge, '|                   ','| CG Prec.Steps=',in%ncong,&
       '|  CG Steps=',in%ncongt
  write(*,'(1x,a,1pe7.1,1x,a,1x,a,i8)')&
       ' elec. field=',in%ef(1),'|                   ','| DIIS Hist. N.=',in%idsx
  if (in%nspin>=2) then
     write(*,'(1x,a,i7,1x,a)')&
          'Polarisation=',2*in%mpol, '|'
  end if
  if (atoms%geocode /= 'F') then
     write(*,'(1x,a,1x,a,3(1x,1pe12.5))')&
          '  Geom. Code=    '//atoms%geocode//'   |',&
          '  Box Sizes (Bohr) =',atoms%alat1,atoms%alat2,atoms%alat3

  end if
end subroutine print_input_parameters
!!***


!!****f* BigDFT/read_atomic_positions
!! FUNCTION
!!    Read atomic positions
!! SOURCE
!!
subroutine read_atomic_positions(iproc,ifile,at,rxyz)
  use module_base
  use module_types
  implicit none
  integer, intent(in) :: iproc,ifile
  type(atoms_data), intent(inout) :: at
  real(gp), dimension(3,at%nat), intent(out) :: rxyz
  !local variables
  character(len=*), parameter :: subname='read_atomic_positions'
  real(gp), parameter :: bohr=0.5291772108_gp !1 AU in angstroem
  character(len=3) :: suffix
  character(len=2) :: symbol
  character(len=20) :: tatonam
  character(len=50) :: extra
  character(len=150) :: line
  logical :: lpsdbl,dowrite
  integer :: nateq,iat,jat,ityp,i,ierror,ierrsfx,i_stat,natpol,j
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
  if (at%units == 'angstroemd0' .or. at%units== 'atomicd0' .or. at%units== 'bohrd0') then
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
  read(ifile,'(a150)')line

  !old format, still here for backward compatibility
  !admits only simple precision calculation
  read(line,*,iostat=ierror) rx,ry,rz,tatonam

  !in case of old format, put geocode to F and alat to 0.
  if (ierror == 0) then
     at%geocode='F'
     at%alat1=0.0_gp
     at%alat2=0.0_gp
     at%alat3=0.0_gp
  else
     if (lpsdbl) then
        read(line,*,iostat=ierrsfx) tatonam,alat1d0,alat2d0,alat3d0
     else
        read(line,*,iostat=ierrsfx) tatonam,alat1,alat2,alat3
     end if
     if (ierrsfx == 0) then
        if (trim(tatonam)=='periodic') then
           at%geocode='P'
        else if (trim(tatonam)=='surface') then 
           at%geocode='S'
           at%alat2=0.0_gp
        else !otherwise free bc
           at%geocode='F'
           at%alat1=0.0_gp
           at%alat2=0.0_gp
           at%alat3=0.0_gp
        end if
        if (.not. lpsdbl) then
           alat1d0=real(alat1,gp)
           alat2d0=real(alat2,gp)
           alat3d0=real(alat3,gp)
        end if
     else
        at%geocode='F'
        at%alat1=0.0_gp
        at%alat2=0.0_gp
        at%alat3=0.0_gp
     end if
  end if

  !reduced coordinates are possible only with periodic units
  if (at%units == 'reduced' .and. at%geocode == 'F') then
     if (iproc==0) write(*,'(1x,a)')&
          'ERROR: Reduced coordinates are not allowed with isolated BC'
  end if

  !convert the values of the cell sizes in bohr
  if (at%units=='angstroem' .or. at%units=='angstroemd0') then
     ! if Angstroem convert to Bohr
     at%alat1=alat1d0/bohr
     at%alat2=alat2d0/bohr
     at%alat3=alat3d0/bohr
  else if  (at%units=='atomic' .or. at%units=='bohr'  .or.&
       at%units== 'atomicd0' .or. at%units== 'bohrd0') then
     at%alat1=alat1d0
     at%alat2=alat2d0
     at%alat3=alat3d0
  else if (at%units == 'reduced') then
     !assume that for reduced coordinates cell size is in bohr
     at%alat1=real(alat1,gp)
     at%alat2=real(alat2,gp)
     at%alat3=real(alat3,gp)
  else
     write(*,*) 'length units in input file unrecognized'
     write(*,*) 'recognized units are angstroem or atomic = bohr'
     stop 
  endif

  at%ntypes=0
  do iat=1,at%nat
     if (ierror == 0) then
        !old case of ascii file, added for backward compatibility
        if (iat /= 1) read(ifile,*) rx,ry,rz,tatonam
     else
        !xyz input file, allow extra information
        read(ifile,'(a150)')line 
        if (lpsdbl) then
           read(line,*,iostat=ierrsfx)symbol,rxd0,ryd0,rzd0,extra
        else
           read(line,*,iostat=ierrsfx)symbol,rx,ry,rz,extra
        end if
        !print *,line
        call find_extra_info(line,extra)
        call parse_extra_info(iproc,iat,extra,at)

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
     if (at%units == 'reduced') then !add treatment for reduced coordinates
        rxyz(1,iat)=modulo(rxyz(1,iat),1.0_gp)
        if (at%geocode == 'P') rxyz(2,iat)=modulo(rxyz(2,iat),1.0_gp)
        rxyz(3,iat)=modulo(rxyz(3,iat),1.0_gp)
     else if (at%geocode == 'P') then
        rxyz(1,iat)=modulo(rxyz(1,iat),alat1d0)
        rxyz(2,iat)=modulo(rxyz(2,iat),alat2d0)
        rxyz(3,iat)=modulo(rxyz(3,iat),alat3d0)
     else if (at%geocode == 'S') then
        rxyz(1,iat)=modulo(rxyz(1,iat),alat1d0)
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
        if (tatonam == atomnames(ityp)) then
           at%iatype(iat)=ityp
           goto 200
        endif
     enddo
     at%ntypes=at%ntypes+1
     if (at%ntypes > 100) stop 'more than 100 atomnames not permitted'
     atomnames(ityp)=tatonam
     at%iatype(iat)=at%ntypes
200  continue

     if (at%units=='angstroem' .or. at%units=='angstroemd0') then
        ! if Angstroem convert to Bohr
        do i=1,3 
           rxyz(i,iat)=rxyz(i,iat)/bohr
        enddo
     else if (at%units == 'reduced') then 
        rxyz(1,iat)=rxyz(1,iat)*at%alat1
        if (at%geocode == 'P') rxyz(2,iat)=rxyz(2,iat)*at%alat2
        rxyz(3,iat)=rxyz(3,iat)*at%alat3
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
     if (iproc == 0) then
        write(*,'(1x,a)')'Control your posinp file, cannot proceed'
        write(*,'(1x,a)',advance='no')&
             'Writing tentative alternative positions in the file posinp_alt...'
        open(unit=9,file='posinp_alt')
        write(9,'(1x,a)')' ??? atomicd0'
        write(9,*)
        do iat=1,at%nat
           dowrite=.true.
           do jat=iat+1,at%nat
              if ((rxyz(1,iat)-rxyz(1,jat))**2+(rxyz(2,iat)-rxyz(2,jat))**2+&
                   (rxyz(3,iat)-rxyz(3,jat))**2 ==0.0_gp) then
                 dowrite=.false.
              end if
           end do
           if (dowrite) & 
                write(9,'(a2,4x,3(1x,1pe21.14))')trim(at%atomnames(at%iatype(iat))),&
                (rxyz(j,iat),j=1,3)
        end do
        close(9)
        write(*,'(1x,a)')' done.'
        write(*,'(1x,a)')' Replace ??? in the file heading with the actual atoms number'               
     end if
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
!!***


!!****f* BigDFT/find_extra_info
!! FUNCTION
!!    Find extra information
!!
!! SOURCE
!!
subroutine find_extra_info(line,extra)
  implicit none
  character(len=150), intent(in) :: line
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
     if (nspace==8) then
        extra=line(i:min(150,i+49))
        exit find_space
     end if
     if (i==150) then
        extra='nothing'
        exit find_space
     end if
     i=i+1
  end do find_space

end subroutine find_extra_info
!!***


!!****f* BigDFT/parse_extra_info
!! FUNCTION
!!    Parse extra information
!!
!! SOURCE
!!
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
  !print *,iat,'ex'//trim(extra)//'ex'
  read(extra,*,iostat=ierr)nspol,nchrg,suffix
  if (extra == 'nothing') then !case with empty information
     nspol=0
     nchrg=0
     suffix='   '
  else if (ierr /= 0) then !case with partial information
     read(extra,*,iostat=ierr1)nspol,suffix
     if (ierr1 /=0) then
        if (trim(extra) == 'f') then
           suffix='f'
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
!!***


!!****f* BigDFT/charge_and_spol
!! FUNCTION
!!   Calculate the charge and the spin polarisation to be placed on a given atom
!!   RULE: natpol = c*1000 + sgn(c)*100 + s: charged and polarised atom (charge c, polarisation s)
!! SOURCE
!!
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
!!***

