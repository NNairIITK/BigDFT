      program pseudo
      use libxcModule


c     Parallel fitting of HGH pseudopotentials using a simplex downhill method.
c     Takes into account multiple atomic references, excitation energies and
c     softness monitored by errors in 3d wavelet transformations.
c     Uses MPI and libXC, supoorts collinear spin polarization as well as 
c     nonlinear core corrections and has a GPU accelerated version of the
c     wavelet part.
c
c     This program is based on the sources available at
c     http://cvs.berlios.de/cgi-bin/viewcvs.cgi/cp2k/potentials/Goedecker/pseudo/v2.2/
c
c     Authors:                  Raffael.Widmer (Cuda Acceleration)
c                             and Alex Willand (all other modifications)  
c                      under the supevision of  
c              Stefan Goedecker, December 2010   
c             Universitaet  Basel, Switzerland 
c                                                          


      implicit real*8 (a-h,o-z)
      logical fullac, avgl1,avgl2,avgl3,plotwf,denbas,mixref,pol,
     :    ortprj, litprj, energ,igrad,verbose, info,exists,ldump

      parameter ( norbmx=40, nrmax=10000, maxdim=30 )
      parameter ( lmx=5, lpmx= 4, noccmx=7, nsmx=2 )
      parameter ( ngmx=32, nintmx=5*(ngmx+14) )

      dimension aeval(noccmx,lmx,nsmx),chrg(noccmx,lmx,nsmx),
     1     dhrg(noccmx,lmx,nsmx),ehrg(noccmx,lmx,nsmx),
     1     res(noccmx,lmx,nsmx),
     1     wfnode(noccmx,lmx,nsmx,3),
     1     psi(0:ngmx,noccmx,lmx,nsmx),wght(noccmx,lmx,nsmx,8),
     1     gpot(4),hsep(6,lpmx,nsmx),r_l(lpmx),
     1     occup(noccmx,lmx,nsmx),
     1     vh((lmx+1)*((ngmx+1)*(ngmx+2))/2,
     1     (lmx+1)*((ngmx+1)*(ngmx+2))/2),xp(0:ngmx),
     1     rmt(nintmx,((ngmx+1)*(ngmx+2))/2,lmx+1),
     1     rmtg(nintmx,((ngmx+1)*(ngmx+2))/2,lmx+1),
     1     wghtcf(20),wghtex(20),
     1     ud(nintmx,((ngmx+1)*(ngmx+2))/2,lmx+1)

      dimension pp(maxdim*(maxdim+1)),yp(maxdim+1),rae(nrmax)
      dimension psiold(nrmax,noccmx,lmx+1,nsmx)
      dimension rr(nintmx),rw(nintmx),rd(nintmx)
      dimension ofdcoef(3,4), psppar(6,4,2)
      dimension excit(200)! just a temporary array for excitations

      dimension no(norbmx),noae(norbmx),lo(norbmx),so(norbmx),
     :     zo(norbmx),ev(norbmx),crcov(norbmx),dcrcov(norbmx),
     :     ddcrcov(norbmx),
     :     gf(nrmax,norbmx,nsmx),
     :     nomax(0:lmx),nomin(0:lmx),hso(6),havg(6),hhsep(6),
     :     time(3)

      character(len=1) :: il(5),ispp,isppp
      character(len=30) :: plotfile
      character(len=80) :: errmsg
      character(len=7) :: is(2)
      character(len=35) :: fname
      character(len=10) :: tname
      character(len=525) :: string
      character(len=10) :: strcyc
      character(len=25) :: form
      character(len=80) :: label
      integer :: namoeb,nsmplx
      integer ::  nconfpaw, nchannelspaw, npawl, pawstN, pawstL, pawstP
      real(8) pawrcovfact
      character(len=125) :: pawstatom
      character(len=8) :: dateYMD


c     INCLUDE 'func.inc'
      include 'mpif.h'
c     initialize some counters
      ntime=0
      itertot=0

c     and time profiling
      time=0d0
      call cpu_time(t0)


c     this is a large number passed as reference when the penalty
c     contributions shall be written out by the penalty subroutine.

      penref=1d100


c     MPI initialization
      iproc=0
      nproc=1
      call MPI_INIT(ierr)
      call MPI_COMM_RANK(MPI_COMM_WORLD,iproc,ierr)
      call MPI_COMM_SIZE(MPI_COMM_WORLD,nproc,ierr)


c     This might not work with all compilers:
c     Generate additional output files for each
c     process iproc > 0 connected to unit 6,
c     which serves as standard output for process 0.
      if(iproc>0)then
         write(label,'(a,i2.2,a)')'proc.',iproc,'.out'
         open(6,file=trim(label))
      end if

c     DEBUG FILES: Intended to split part of the output per process       
c         write(label,'(a,i2.2,a)')'dbg.',iproc,'.out'
c         open(66,file=trim(label))
c
      write(6,*) '*********************************************'
      write(6,*) '***              pseudo_2.5               ***'
      write(6,*) '***              fitting of               ***'
      write(6,*) '***   goedecker type pseudopotentials     ***'
      write(6,*) '***   last changes:     december 2010     ***'

      if(iproc>0)then 
         write(6,'(a,i2,a,i2,a)')          
     :          ' ***   output file for process ',iproc, ' of ',
     :                                           nproc,'    ***'
      elseif(nproc>1)then
         write(6,'(a,i2,a)') 
     :          ' ***   parallel run with ', nproc,
     :                                     ' processes      ***'
      end if
      write(6,*) '*********************************************'
      write(6,*)

c     Optional info when call system is allowed
      call system('hostname >> my.hostnames')


c
c     default values:
c
      NAMOEB=0
      nsmplx=2000
      ng=0
      rij=0
      FULLAC=.false.
      denbas=.false.
      avgl1=.false.
      avgl2=.false.
      avgl3=.false.
      plotwf=.false.
      ortprj=.false.
      litprj=.false.
      energ=.false.
      igrad=.false.
      verbose=.false.
      info=.false.
      mixref=.false.
      ldump=.false.

      npawconf=-1
      npawl   = 3
      nchannelspaw = 4
      pawstatom=''
      pawrcovfact=1.0_8
      
c
c     DO NOT read command line options as in previous versions.
c     Instead, read those options from the first line of a file input.dat
c     This modification is needed for a parallel run.

       write(6,*)'Reading the file input.dat'
       write(6,*)'__________________________'
       write(6,*)


       open(11,file='input.dat')
       read(11,'(a)',iostat=ierr)string
       if(ierr/=0)then
          write(6,'(a,i2)')'Could not read the file input.dat'
          string=''
       end if

       if(trim(string)=='')then  
            write(6,*) 'Could not read any options from the'//
     :                 ' first line of the file input.dat.'
            write(6,*) 'possible options: -cN       ',
     :           'number of fitting cycles (-1< N< 10000)'
            write(6,*) '                  -nN       ',
     :           'number of iterations per simplex'
            write(6,*) '                  -gN       ',
     :           'number gaussians'
            write(6,*) '                  -rN       ',
     :           'set rij=N/10 '
            write(6,*) '                  -fullacc  ',
     :           'use max. number of gaussians'
            write(6,*) '                  -orth     ',
     :           'orthogonalisation of the projectors'
            write(6,*) '                  -lith     ',
     :           'transformation of the projectors as in lit.'
            write(6,*) '                  -denbas   ',
     :           'use dense gaussian basis    '
            write(6,*) '                  -mixref   ',
     :           'allow contradictions in AE ref data'
            write(6,*) '                  -info     ',
     :           'write gaussian coefficients of final'
            write(6,*) '                            ',
     :           'wavefunctions to gauss.par'
            write(6,*) '                  -plot     ',
     :           'plot wfs after each iteration'
            write(6,*) '                  -lNso     ',
     :           'average nonlocal potential = 0 for l=N'
            write(6,*) '                            ',
     :           '(only for the highest projector)'
            write(6,*) '                  -pawN       ',
     :           ' no opt.,calculate  pawpatch projectors for'//
     :           ' the Nth configuration '
            write(6,*) '                  -noflpawN       ',
     :           ' pawpatch patches for the first N Ls (defaults to 3)'
            write(6,*) '                  -nchannelspawN       ',
     :           ' set number of paw projectors to N (defaults to 4)'

            write(6,*) '                  -pawstatomName       ',
     :           '  file named Name will be read for initial ' , 
     :           '  wave funct. '
            write(6,*) '                  -pawstnN       ',
     :           '  initial wave function has first quantum number N ' 
            write(6,*) '                  -pawstlL       ',
     :           '  initial wave function has angular momentum L' 
            write(6,*) '                  -pawstpP       ',
     :           '  initial wave function is multiplied by r**P' 
            write(6,*) '                  -pawrcovfactF     ',
     :      'Rbox for paw is equal to rcov*pawrcovfact. Defaults to 1' 



        goto 11
        end if

c      No loop over command line arguments needed.

            ii=index(string,'-c')
            if (ii.ne.0) then
               label=string(ii+2:min(ii+10,520))
               read(label,*)namoeb
               write(6,*) namoeb, 'fit cycles'
            endif
            ii=index(string,'-orth')
            if (ii.ne.0) then
               ortprj=.true.
               write(6,*) 'orthogonalize the projectors'
            endif
            ii=index(string,'-lith')
            if (ii.ne.0) then
               litprj=.true.
               write(6,*) 'transform the projectors as in ',
     :              'literature'
            endif

            ii=index(string,'-n')
            if (ii.ne.0 .and.  string(1:2).eq."-n"  .and. string(3:3) 
     :   .ne. 'o' .and.
     :            string(3:3) .ne.'c') then
               print *, label
               label=string(ii+2:min(ii+10,520))
               read(label,*)  nsmplx
               write(6,*) nsmplx, 'max. simplex iterations'
            endif

            ii=index(string,'-g')
            if (ii.ne.0) then
               label=string(ii+2:min(ii+10,520))
               read(label,*)ng
               write(6,*)ng,"gaussians (don't use value from atom..ae)"
            endif
            ii=index(string,'-r')
            if (ii.ne.0) then
               label=string(ii+2:min(ii+10,520))
               read(label,*)rij
               write(6,*)rij,' rij (don''t use value from atom..ae)'
            endif
            ii=index(string,'-fullacc')
            if (ii.ne.0) then
               fullac=.true.
               write(6,*) 'Use max. number of gaussians'
               
            endif
            ii=index(string,'-plot')
            if (ii.ne.0) then
               plotwf=.true.
               write(6,*) 'plot wfs after each iteration'
               
            endif
            ii=index(string,'-denbas')
            if (ii.ne.0) then
               denbas=.true.
               write(6,*) 'use dense gaussian basis'
               
            endif
            ii=index(string,'-mixref')
            if (ii.ne.0) then
               mixref=.true.
               write(6,*) 'allow mixed reference data'
               
            endif
            ii=index(string,'-info')
            if (ii.ne.0) then
               info=.true.
               write(6,*) 'write gaussian coefficients of final'
               write(6,*) 'wavefunctions to gauss.par'
              
            endif
            ii=index(string,'-l1so')
            if (ii.ne.0) then
               avgl1=.true.
               write(6,*) 'average nonlocal potential zero for l=1'
               write(6,*) '(only for the highest projector )'
               
            endif
            ii=index(string,'-l2so')
            if (ii.ne.0) then
               avgl2=.true.
               write(6,*) 'average nonlocal potential zero for l=2'
               write(6,*) '(only for the highest projector)'
               
            endif
            ii=index(string,'-l3so')
            if (ii.ne.0) then
               avgl3=.true.
               write(6,*) 'average nonlocal potential zero for l=3'
               write(6,*) '(only for the highest projector)'
               
            endif
            ii=index(string,'-dump')
            if (ii.ne.0) then
               ldump=.true.
               write(6,*) 'dumpfile for ekin.test.f requested'

            endif
            ii=index(string,'-paw')
            if (ii.ne.0) then
               label=string(ii+4:min(ii+12,520))
               read(label,*) nconfpaw
               write(6,*)  'will calculate pawpatches for conf No ',
     :                 nconfpaw
               

            endif
            ii=index(string,'-noflpaw')
            if (ii.ne.0) then
               label=string(ii+8:min(ii+16,520))
               read(label,*) npawl
               write(6,*)  'will calculate paw patches for the',
     :          npawl, ' first Ls '

            endif

            ii=index(string,'-nchannelspaw')
            if (ii.ne.0) then
               label=string(ii+13:min(ii+21,520))
               read(label,*)  nchannelspaw
               write(6,*)  'will consider',
     :         nchannelspaw , ' paw channels  '

            endif

            ii=index(string,'-pawstatom')
            if (ii.ne.0) then
               pawstatom=trim(string(ii+10:min(ii+130,520)))
               ii=index(pawstatom,' ')
               if(ii.ne.0) then
                  pawstatom=trim(pawstatom(:ii-1))
               endif
               write(6,*)  'will consider  ',
     :        trim(pawstatom) ,'file for reading the initial potential '
            endif


            ii=index(string,'-pawstn')
            if (ii.ne.0) then
               label=string(ii+7:min(ii+15,520))
               read(label,*)  pawstN
               write(6,*)  ' N of st. wf. is ', pawstN
            endif

            ii=index(string,'-pawstl')
            if (ii.ne.0) then
               label=string(ii+7:min(ii+15,520))
               read(label,*)  pawstL
               write(6,*)  ' L of st. wf. is ' , pawstL
            endif

            ii=index(string,'-pawstp')
            if (ii.ne.0) then
               label=string(ii+7:min(ii+15,520))
               read(label,*)  pawstP
               write(6,*)  ' initial wf radial part is ',
     :          ' multiplied by r**' , pawstP
            endif

            ii=index(string,'-pawrcovfact')
            if (ii.ne.0) then
               label=string(ii+12:min(ii+20,520))
               print *, label
               read(label,*)  pawrcovfact
               write(6,*)  ' rbox is rcov   ',
     :          ' multiplied by ' , pawrcovfact
            endif

c           End of rading options from input.dat's first line 
 11         continue

      if (nconfpaw /= -1) then
         namoeb=0
         write(6,*) ' fitting disactivated because paw option is active'
      endif

      if (namoeb.eq.0) then
         write(6,*) 'Do one pseudopotential calculation.'
         write(6,*) 'No fitting.'
      endif

      if (ortprj .and. litprj ) then
         write(6,*) 'use only one option -orth or -lith!'
c        stop
      endif


     
c       Read additional parameters for wavelet tranforms
c       from the following lines of the file input.dat

c       grid spacing range: samples, min, max, power law
        read(11,*,iostat=ierr) nhgrid, hgridmin,hgridmax, nhpow 
        read(11,*,iostat=jerr) ampl,crmult,frmult  
        if(ierr.ne.0.or.jerr.ne.0)then
                write(6,*)'input.dat does not contain input '//
     :         'for wavelets. No fitting of softness possible.'
                nhgrid=0
        else
         write(6,'(a,i3,2(1pe12.4),i3)')'hgrid parameters ',
     :                    nhgrid, hgridmin,hgridmax, nhpow
         write(6,'(a,3(1pe12.4))')'offset and radii ', ampl,crmult,frmult
        end if
c       random shift, coarse and fine grid radii
        close(11)
     


c
c     ----------------- read data from AE calculation -----------------
c
      write(6,*)
      write(label,'(a,i2.2,a)')    'atom.',iproc,'.ae'
      write(6,*)'Reading data from '//trim(label)
      write(6,*)'____________________________'
      write(6,*)
      inquire(file=label, exist=exists)
      ierr=0
      if(.not.exists)then
          ierr=2
          write(6,*)
          write(6,*)'No such file! Attempt to read atom.00.ae instead.' 
          write(6,*)'This is a waste of resources, as you could treat'
          write(6,*)'more excitations or different parameters (rprb,'
          write(6,*)'rloc, or even iXC) with this number of processes'
          label='atom.00.ae'
          inquire(file=label, exist=exists)
          if(.not.exists)then
             write(6,*)'Cannot proceed. File atom.00.ae not found. '
             ierr=3
          end if
      end if
      write(errmsg,*)'atomic reference file missing!'
      call errorhandler(ierr,iproc,nproc,errmsg)

      open(unit=40,file=trim(label),form='formatted',status='unknown')
      read(40,*,iostat=ierr) norb, wghtconf
         if(ierr/=0)ierr=3
         write(errmsg,*)'error: 1st line of AE ref data'
         call errorhandler(ierr,iproc,nproc,errmsg)
      read(40,*,iostat=ierr) znucp,zionp,rcovp,rprbp
         if(ierr/=0)ierr=3
         write(errmsg,*)'error: 2nd line of AE ref data'
         call errorhandler(ierr,iproc,nproc,errmsg)
      read(40,'(a)',iostat=ierr) label
         if(ierr/=0)ierr=3
         write(errmsg,*)'error: 3rd line of AE ref data'
         call errorhandler(ierr,iproc,nproc,errmsg)
      j=1

      do i=len(label),1,-1
         if (label(i:i).ne.' ') j=i
      enddo
      isppp=label(j:j)
      read(40,'(a)',iostat=ierr) label
      j1=1
      j2=2
      do i=len(label),1,-1
         if (label(i:i).ne.' ') j1=i
      enddo
      do i=len(label),j1,-1
         if (label(i:i).eq.' ') j2=i
      enddo
      j2=j2-1

c     READING of iXC
c     for now, only keep the two most commonly used functionals
c     backwards compatible. Otherwise, require ABINITs iXC < 0

c     icorrp=label(j1:j2)
      if    (label(j1:j2)=='PADE')then
         iXC=-20
      elseif(label(j1:j2)=='PBE')then   
         iXC=-101130
      else
         read(label(j1:j2),*,iostat=ierr)iXC
         if(ierr/=0)then
          write(6,*)'Could not read the XC input in psppar'
          stop
         end if
      end if

      write(6,*)
c      read(40,'(t2,a)',iostat=ierr) isppp
c      read(40,'(t2,a)',iostat=ierr) icorrp
      read(40,*,iostat=ierr) ngrid
         if(ierr/=0.or.ngrid .gt. nrmax )ierr=3
         if(ngrid .gt. nrmax )
     :   write(6,*)'Input number value is to large.'
         write(errmsg,*)'error: nr gridpoints in AE ref'
         call errorhandler(ierr,iproc,nproc,errmsg)
      write(6,*)
      write(6,*)'pseudo states = ', norb
      write(6,*)'znuc          = ', znucp
      write(6,*)'zpseudo       = ', zionp
      write(6,*)'r_covalent    = ', rcovp
      write(6,*)'r_confining   = ', rprbp
      write(6,*)'ispp          = ', isppp
      write(6,*)'gridpoints    = ', ngrid
      il(1) = 's'
      il(2) = 'p'
      il(3) = 'd'
      il(4) = 'f'
      il(5) = 'g'
      nspin=1
      is(1) = '  so=0'
      if (isppp.eq.'r'.or.isppp.eq.'s') then
         nspin=2
         is(1)= 'so=+0.5'
         is(2)= 'so=-0.5'
      endif
c     for now, use this convention
      pol=.false.
      nspol=1
      if (isppp.eq.'s') then
         pol=.true.
         nspol=2
      end if

      write(6,*)
      write(6,'(a,i7,a,i5)')
     :     ' Initializing libXC with iXC =', iXC,'; nspol =',nspol
c      the barriers are here because the init routine writes to stout
c      with each process. Need to find a way to fix this.
       if (nproc > 1) call MPI_BARRIER(MPI_COMM_WORLD,ierr)  
          call libxc_functionals_init(iXC,nspol)
       if (nproc > 1) call MPI_BARRIER(MPI_COMM_WORLD,ierr)  
         
        
      write(6,*)
      read(40,*) !some comment line
      write(6,*)' nl    s      occ        ',
     :                 'eigenvalue     charge(rcov)    '
c      write(6,*)' nl    s      occ        ',
c     :     'eigenvalue     charge(rcov)    ',
c     :     'dcharge         ddcharge'


c     this file should hold all ae plots
      if(plotwf.and.iproc==0) open(unit=41,file='ae.orbitals.plt')
c     read the AE data
      do iorb=1,norb
         read(40,*,iostat=ierr) no(iorb),lo(iorb),so(iorb),zo(iorb),
     :        ev(iorb),crcov(iorb),dcrcov(iorb),ddcrcov(iorb)!,plotfile
c        write(6,*,iostat=ierr) no(iorb),lo(iorb),so(iorb),zo(iorb),
c    :        ev(iorb),crcov(iorb),dcrcov(iorb),ddcrcov(iorb),plotfile
         if(ierr/=0)ierr=3
         write(errmsg,*)'reading error in AE ref data'
         call errorhandler(ierr,iproc,nproc,errmsg)
         write(6,30) no(iorb),il(lo(iorb)+1),so(iorb),zo(iorb),
     :        ev(iorb),crcov(iorb)
c    :        ev(iorb),crcov(iorb),dcrcov(iorb),ddcrcov(iorb)
 30      format(1x,i1,a1,f6.1,f10.4,2f16.10,2f16.7,a)
         if(plotwf.and.iproc==0)then
c           use this statement if atom is compiled to write one file
c           per configuration and orbital
c           open(41,file=trim(plotfile))
            read(41,*)
            read(41,*)
            do igrid=1,ngrid
              read(41,*,iostat=ierr) rae(igrid),(gf(igrid,iorb,igf),
     :          igf=1,nspin)
c             error handling in the loop is slow, but better give detailed feedback
c             for now, we only plot the ground state, i.e. atom.00.ae
              if(ierr/=0)ierr=2
c              write(errmsg,*)'error reading AE plots',
c    :                    trim(plotfile),
c     :                   'orb',iorb,'pt',igrid
c              call errorhandler(ierr,iproc,nproc,errmsg)
            end do
            read(41,*)
c           do not close this unit when reading from one single file
c           close(unit=41)
         end if
      enddo
      if(plotwf) close(unit=41)
      goto 457
c 456  write(6,*) 'error during reading atom.ae'
c      stop
 457  continue
      lmax=0
      lcx=0
      do iorb=1,norb
         lmax=max(lo(iorb),lmax)
         if (zo(iorb).gt.1.d-10)lcx=max(lo(iorb),lcx)
      enddo
c     print*,'lmax=',lmax
c     print*,'lcx=',lcx, '( charge > 1.0d-10)'
      if(lmax.gt.lmx+1)ierr=3
      write(errmsg,*)'array dimension problem:lmax'
      call errorhandler(ierr,iproc,nproc,errmsg)
c     compute corresponding n-quantum numbers of the pseudostates
c     no()   will contain n quantum numbers of the pseudostates afterwards
c     noae() will contain n quantum numbers of the ae-states afterwards
c     no() starts from n=1 for each(!) l
      noccmax=0
      do l=0,lmax
         nomin(l)=100
         nomax(l)=0
      enddo
      do iorb=1,norb
         nomin(lo(iorb))=min(no(iorb),nomin(lo(iorb)))
         nomax(lo(iorb))=max(no(iorb),nomax(lo(iorb)))
      enddo
      do iorb=1,norb
         noae(iorb)=no(iorb)
         no(iorb)=no(iorb)-nomin(lo(iorb))+ 1
      enddo
      do l=0,lmax
         noccmax= max(noccmax,nomax(l)-nomin(l)+1)
      enddo
      write(6,*) 'noccmax ', noccmax
      if (noccmax.gt.noccmx)ierr=3 
      write(errmsg,*) 'array dimension problem: noccmax'
      call errorhandler(ierr,iproc,nproc,errmsg)
c      print*,'noccmax=',noccmax
      do nocc=1,noccmax
         do l=0,lmax
            do ispin=1,nspin
               occup(nocc,l+1,ispin)=0.0d0
               aeval(nocc,l+1,ispin)=0.0d0
               chrg (nocc,l+1,ispin)=0.0d0
               dhrg (nocc,l+1,ispin)=0.0d0
               ehrg (nocc,l+1,ispin)=0.0d0
            enddo
         enddo
      enddo
      !why not initialize ALL elements of occup with zeroes?
                              occup=0d0
      do iorb=1,norb
         nocc=no(iorb)
         l=lo(iorb)
         ispin=1
         if (so(iorb).lt.0) ispin=2
         occup(nocc,l+1,ispin)=zo(iorb)
c        TEST  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!|||||||||||||||||||||||||||||||||
         ispin=2
         if (so(iorb).le.0) ispin=1


         do j=1,ngrid
            psiold(j,nocc,l+1,ispin) = 0.0d0
c     use major comp. as reference
            if (rae(j).ne.0.0)
     :           psiold(j,nocc,l+1,ispin)=gf(j,iorb,1)/rae(j)
         enddo
      enddo

      write(6,*)
      write(6,*) 'All electron and pseudo-wfn quantum numbers'
      write(6,*) '        n(AE)   l   inl(PS)   '
      do iorb=1,norb
         write(6,'(6x,3i6)')  noae(iorb),lo(iorb), no(iorb)
      enddo



c
c     weights will be read from weights.par
c

c     pseudo 2.4 was backwards compatible with this files
c     reading conventions from pseudo2.2 and 2.3.
c     For this version, this is not the case anymore!

         write(6,*)
         write(6,*) 'Reading actual weights from file weights.par'
         write(6,*) '____________________________________________'
         write(6,*)
         open(unit=24,file='weights.par',form='formatted')
         read(24,*)
         wghtp0   = 0d0
         wghtsoft = 0d0
         wghtrad  = 0d0
         wghthij  = 0d0
         wghtexci = 0d0
c        what about a different order for reading these?
         read(24,*,iostat=ierr)wghtp0,wghtsoft,wghtrad,wghthij,wghtexci
c        no need to call error handler here, shared input file
         if(ierr/=0)
     :   write(6,*)'Reading error for weights of psi(r=0) and softness.'
         write(6,'(a,1pe20.10)') 
     :        ' Weight for psi(r=0)=0 is     ',wghtp0
         write(6,'(a,1pe20.10)') 
     :        ' for Ekin(gauss-wavelet)      ',wghtsoft
         write(6,'(a,1pe20.10)') 
     :        ' for small radii (norm prj)   ',wghtrad
         write(6,'(a,1pe20.10)') 
     :        ' for large hij  (instable)    ',wghthij
         write(6,'(a,1pe20.10)') 
     :        ' and for excitation energies  ',wghtexci
            
         read(24,*)!comment line
        
c        read the weights for eigenvalues and integrals into the array wght

         do iorb=1,norb
            nocc=no(iorb)
            l=lo(iorb)
            ispin=1
            ss = so(iorb)
            if (ss.lt.0) ispin=2
            read(24,*) nw,lw,sw,(wght(nocc,l+1,ispin,i),i=1,8)
            if (noae(iorb).ne.nw .or. l.ne.lw .or. ss.ne.sw) then
               write(6,*) 'error in file ''weights.par'' '
               write(6,*) 'need weights for n,l,s:',
     :              noae(iorb),l,so(iorb)
               write(6,*) 'found            n,l,s:',nw,lw,sw
               if(nproc>1) call MPI_FINALIZE(ierr)
               stop
            endif
         enddo
         close(24)

c     Read compute excitation energies from the last line of atom.??.ae

      if(nproc>1)then
         write(6,*)
c        read etotal and exchange data with all processes
c        it would be enough to broadcast etot of system 0,
c        but this will not take much time and allow some 
c        output and proofreading.
         excit=0d0
         read(40,*,iostat=ierr)excit(1)
         write(6,*)'all electron energy (Ryd):',excit(1)
         if(ierr/=0)then
           write(6,*)
           write(6,*)'               WARNING'
           write(6,*)'The last line of the atomic reference file'
           write(6,*)'must specify the total energy.'
         end if
         excit(nproc+1:2*nproc)=excit(1)
         call MPI_ALLTOALL(
     :          excit(nproc+1:2*nproc),1,MPI_DOUBLE_PRECISION, 
     :          excit(1:nproc),1,MPI_DOUBLE_PRECISION,
     :          MPI_COMM_WORLD,ierr)
         
         if(ier/=0)write(6,*)'           WARNING: MPI error'        
c        total energies to excitation energies
c        different conventions of etot in gatom and atom.f
         excit=0.5d0*(excit-excit(1))
         excitAE=excit(iproc+1)
         if(ierr/=0)excitAE=0d0
         write(6,*)
         write(6,*)'excitation energies (AE)'
         do i=0,nproc-1
             write(6,'(1x,i3,3(1x,1pe20.10))')
     :       i, excit(i+1)
         end do
         ierr=0
         if(excitAE==0d0.and.iproc/=0)ierr=1
         write(errmsg,*)'Excitation energy is zero and thus ignored'
         call errorhandler(ierr,iproc,nproc,errmsg)
      else
         excitAE = 0d0
      end if
c     done reading atomic ref data
      close(40)




c     ---------------------------------------------------------------
c     main loop begins here
c     ---------------------------------------------------------------

      do iiter=1,max(namoeb,1)

c        If the excitation energy has a nonzero weight
c        we always need to compute the total energy. 
c        The ground state energy is needed as a reference
c        for the excitation  energies, thus iproc==zero
         energ=(wghtexci>0d0.or.(iproc==0.and.nproc/=1))



c     read initial pseudopotential parameter from psppar
c     test if parameter from atom.ae and psppar are consistent
c     notice that fitting to inconsistent AE data MIGHT be desired
c     in a parellel run, for instance several values for rprb.



         if (iiter.eq.1) then
            write(6,*)
            write(6,*) 'Reading data from psppar'
            write(6,*) '________________________'
            write(6,*) 

           open(99,file='vertex.dump')

c     the format has been changed from previous version
c     to allow direct compatibility with ABINIT pseudos!

c     Output will always be in HGH-k format (pspcod=10)
c     Input formats are supported as GTH (2) and HGH (3) and HGH-k (10)
c     Notice some additional data needs to be read from this file.

            inquire(file='psppar',exist=exists) 
            if(.not.exists)then
c           no need to call errorhandler, shared file
                write(6,*)'The file psppar is lacking.'
                write(6,*)'Pseudo potentials are available from'
                write(6,*)'http://www.abinit.org/downloads/psp-links'
                if(nproc>1) call MPI_FINALIZE(ierr)
                stop
            end if


            open(unit=11,file='psppar',form='formatted',
     1           status='unknown')


c           the very first character (when reading in list format)
c           defines the method:
c             n   non relativistc
c             r   relativistc
c             s   spin polarized (and relativistic)
c           

c           read(23,'(t2,a,t10,a)') ispp,icorr
            read(11,*,iostat=ierr) label, ng, tmprij, rcov, rprb
            if(ierr/=0)then
c               no need to call error handler, shared input file
c               thus some stop statements have not been eliminated here
                write(6,*)
                write(6,*)'             WARNING'
                write(6,*)
                write(6,*)'The first line of psppar MUST specify'
                write(6,*)'the spin treatment and some more params'
                write(6,*)'for the fit: ispp  ng rij rcov and rprb.'
                write(6,*)'The leading character of the first word'
                write(6,*)'ispp must be n, r or p. ng and rij define'
                write(6,*)'the gaussian basis, rcov the charge inte-'
                write(6,*)'grals and rprb the confining potential.'
                write(6,*)
                write(6,*)'          the first line was:'
                write(6,*)trim(label), ng, tmprij, rcov, rprb
                if(nproc>1) call MPI_FINALIZE(ierr)
                stop
            end if
            j=1
            do i=len(label),1,-1
               if (label(i:i).ne.' ') j=i
            enddo
            ispp=label(j:j)
            ierr=0
            if(ispp.ne.'r'.and.ispp.ne.'n'.and.ispp.ne.'s')then
               write(6,*)'The first non-blank character of psppar'
               write(6,*)'must be one of' 
               write(6,*)'n: for non relativistic calculations'
               write(6,*)'r: for relativistic unpolarized calculations'
               write(6,*)'s: for spin polarized calculations'
               write(6,*)
               write(6,*)'Character found:',ispp
               write(6,*)'Exiting.'
               call MPI_FINALIZE(ierr)
               stop
            end if
            if (isppp.ne.ispp) ierr=3
            write(errmsg,*)'inconsistent spin treatment!'
            call errorhandler(ierr,iproc,nproc,errmsg)
c              below option does not really make sense.
c              it could actually be useful, but needs testing for
c              allocation errors and other causes of trouble!

c              if (mixref) then
c                 write(6,*) 'Warning! continue program using ispp',
c    :                 ' from psp.par'
c              else
c                 write(6,*) 'option ''-mixref'' allows such settings'
c                 stop
c              endif
 
            read(11,*,iostat=ierr) znuc, zion
            if(ierr/=0)then
c               no need to call error handler, shared input file
c               thus some stop statements have not been eliminated here
                write(6,*)
                write(6,*)'             WARNING'
                write(6,*)'Could not read nuclear and valence charges'
                write(6,*)'on the second line of psppar.'
                if(nproc>1) call MPI_FINALIZE(ierr)
                stop
            end if


 
            read(11,*,iostat=ierr) ipspcod, ixcpp
            if(ierr/=0)then
c               no need to call error handler, shared input file
                write(6,*)
                write(6,*)'             WARNING'
                write(6,*)'Could not read PSP format and iXC from'
                write(6,*)'the third line of psppar.'
                if(nproc>1) call MPI_FINALIZE(ierr)
                stop
            end if

c           already read from psppar
c           read(23,*) nng,tmprij
            if (ng.eq.0) ng=nng
            if (rij.eq.0) rij=tmprij
            if (fullac) then
               ng=ngmx
            elseif (ng .gt. ngmx )then 
               write(6,*) 'gaussians:',ng
               write(6,*) 'max is   :',ngmx
                if(nproc>1) call MPI_FINALIZE(ierr)
               stop
            endif
            if (noccmax.gt.ng+1)then
               write(6,*) 'noccmax>ng+1'
               write(6,*) 'ng+1,rij ',ng+1,rij
               if(nproc>1) call MPI_FINALIZE(ierr)
               stop 
            end if


c           already read from psppar
c           read(23,*) rcov,rprb
            ierr=0
            if (rcov.ne.rcovp) then
               ierr=1
               write(6,*)'rcov from atom.ae and psppar not identical'
               write(6,*) 'atom.ae   rcov=',rcovp
               write(6,*) 'psppar    rcov=',rcov
               if (mixref) then
                  write(6,*)'No excitation energies for mixed rcov'
                  excitAE=0d0
                  rcov=rcovp
                  ierr=1
               else
                  ierr=3
                  write(6,*) 'option ''-mixref'' allows mixed rcov'
               endif
            endif
            write(errmsg,*)'different integration radii found'
            call errorhandler(ierr,iproc,nproc,errmsg)
            
            if (3*rcov > max(crmult,frmult) ) then
               ierr=1
               write(6,*)'                NOTE'
               write(6,*)
               write(6,*)'The localization radii for the wavelet basis'
               write(6,*)'are smaller than 3 times rcov'
               write(6,*)'psppar       rcov=',rcov
               write(6,*)'input.dat  crmult=',crmult
               write(6,*)'input.dat  frmult=',frmult
               write(6,*)'If the curve dEkin.orb.dat  gets bumpy,'
               write(6,*)'try raising these two radii'
            endif
c           no error handler needed, same for all processes.

            ierr=0
            if (rprb.ne.rprbp) then
               write(6,*)'rprb in atomic reference differs', 
     :                   ' from the value in psppar.'
               write(6,*) 'atom.ae   rprb=',rprbp
               write(6,*) 'psppar    rprb=',rprb

               if (mixref) then
                  write(6,*)'No excitation energies for mixed rprb'
                  excitAE=0d0
                  rprb=rprbp
                  ierr=1
               else
                  write(6,*) 'option ''-mixref'' allows mixed rprb'
                  ierr=3
               endif
            endif
            write(errmsg,*)'different confining potential found' 
            call errorhandler(ierr,iproc,nproc,errmsg)

            is(1) = '  so=0'
            if (ispp.eq.'r'.or.ispp.eq.'s') then
               nspin=2
               is(1)= 'so=+0.5'
               is(2)= 'so=-0.5'
            endif
c           read(23,'(a)') label
c           j1=1
c           j2=2
c           do i=len(label),1,-1
c              if (label(i:i).ne.' ') j1=i
c           enddo
c           do i=len(label),j1,-1
c              if (label(i:i).eq.' ') j2=i
c           enddo
c           j2=j2-1
c           icorr=label(j1:j2)
c           if (icorr.ne.icorrp) then
            ierr=0
            if (iXC.ne.iXCpp) then
                write(6,*) 'contradiction in exchange correlation'
                write(6,*) 'atom.ae   iXC=',iXC
                write(6,*) 'psppar    iXC=',iXCpp
                if (mixref) then
                   write(errmsg,*)'no excitation energy for this system'
                   ExcitAE=0d0
                   ierr=1
                else
                   write(errmsg,*) 'option ''-mixref'' allows mixed iXC'
                   ierr=3
                endif
            endif
            write(errmsg,*) 'Mixed iXC in atomic reference'
            call errorhandler(ierr,iproc,nproc,errmsg)

c     here we previously set the parameter/variables for the xc-functional(s)
c     

c     NOTE: znuc MUST be consistent, mixed references make no sense here
c           zion WILL be different as soon as we have ionic configurations. 

c           read(11,*) znuc, zion, rloc, gpot(1),gpot(2),gpot(3),gpot(4)
            ierr=0
            if (znucp.ne.znuc) then
               write(6,*) 'znuc from atom.ae and psppar not identical'
               write(6,*) 'atom.ae   znuc=',znucp
               write(6,*) 'psppar    znuc=',znuc
               ierr=3
            endif
            write(errmsg,*) 'nucleonic charge differs from AE data'
            call errorhandler(ierr,iproc,nproc,errmsg)
            ierr=0
            if (zionp.ne.zion) then
               write(6,*) 'zion from atom.ae and psppar not identical'
               write(6,*) 'atom.ae  zion=',zionp
               write(6,*) 'psppar   zion=',zion
               ierr=1
c              zion=zionp
            endif
            write(errmsg,*) 'valence charge differs from AE data'
            call errorhandler(ierr,iproc,nproc,errmsg)

c           read the pseudopotential the way BigDFT does

c           be sure to have zeroes for undefined entries.
            psppar = 0d0
            rcore  = 0d0
            zcore  = 0d0

            if(ipspcod==10.or.ipspcod==11)then
                 write(6,*)'HGH matrix format'
c              ! HGH-K format: all projector elements given.
               ! dimensions explicitly defined for nonzero output.

c              ! local part
               read(11,*) rloc,nloc,(gpot(j),j=1,nloc) 
c              ! separable part
c              ! relativistic; hij are averages, kij separatons of hsep
               read(11,*,iostat=ierr) lpx,rcore,zcore !number of channels of the pseudo
               lpx=lpx-1      !max angular momentum 
               if (lpx .gt. lmx ) then
                 write(6,*) 'array dimension problem: lpx,lpmx',lpx,lpmx
                 if (nproc > 1) call MPI_FINALIZE(ierr)
                 stop
               end if
               do l=1,lpx+1 !l channels
                  read(11,*) r_l(l),nprl,
     :                   psppar(1,l,1),(psppar(j+2,l,1),
     :                                              j=2,nprl)  !h_ij 1st line
                  do i=2,nprl
c                    spin up
                     read(11,*) psppar(i,l,1),(psppar(i+j+1,l,1),
     :                                              j=i+1,nprl)!remaining h_ij 
                  end do
c                 This concept has been discarded:
c                 NEW convention: pspcod = 11 has kij elements for s,  too
c                 if (l==1.and.ipspcod==10) cycle
                  if (l==1) cycle
                   do i=1,nprl
                      read(11,*) psppar(i,l,2),(psppar(i+j+1,l,2),
     :                                       j=i+1,nprl)!all k_ij
                   end do
               end do ! loop over l 
            elseif(ipspcod==3)then
                 write(6,*)'HGH diagonal format'
c                HGH diagonal part case
c                technically, lpx is fixed at
                 lpx=3
                 read(11,*) rloc,(gpot(j),j=1,4)
                 read(11,*) r_l(1),psppar(1:3,1,1) 
                 do l=2,4
                    read(11,*) r_l(l),psppar(1:3,l,1) 
                    read(11,*)        psppar(1:3,l,2) 
                 end do
            elseif(ipspcod==2)then
                 write(6,*)'GTH format'
c                ! GTH case
c                technically, lpx is fixed at
                 lpx=1
                 read(11,*) rloc,(gpot(j),j=1,4)
                 read(11,*) r_l(1),psppar(1:2,1,1)
                 read(11,*) r_l(2),psppar(1  ,2,1)
            else
c                no need to call error handler, shared input file
                 write(6,*)'               WARNING'
                 write(6,*)'pspcod (1st number of 3rd line) read from' 
                 write(6,*)'psppar is unknown or not supported.' 
                 write(6,*)'supported are 2,3, or 10, not ',ipspcod
                 if (nproc > 1) call MPI_FINALIZE(ierr)
                 stop
            end if
c           done with reading psppar
            close(11) 

c           avoid radii equal zero, even for unused projectors. 
            do l=1,lpx+1
              if(r_l(l)==0d0)then
                write(6,*)'all r_l should be nonzero.'
                write(6,*)'The r_l of the ',il(l),'-projector has been'
                write(6,*)'adjusted from 0 to 1. Check your psppar.'
                r_l(l)=1d0
              end if
            end do
c           and to be very sure
c           in case lpx increses later
            r_l(lpx+2:lpmx)=1d0         


c           Then rearrange pspcod into hsep accordingly and 
c           convert hij,kij to hij(up,down)  

c           pspcod  as read are packed upper diagonal ROW elements of

c               h =  ((l+1) hup + l hdn)/(2l+1)
c               k =      2( hup -   hdn)/(2l+1)

c           we want hsep,  packed upper diagonal COL elements of

c             hup = h +   l/2   k
c             hdn = h - (l+1)/2 k


c           1) get offdiagonal elements where needed
c           2) permute from row to col ordering
c           3) transform from h/k to up/down


c           just to be sure no rubbish will be added 
            hsep=0d0
            
            if (ipspcod == 2) then !GTH case
c              offdiagonal elements are zero per definition.
c              simply rearrange hij and fill zero elements
               do l=1,lpx+1
                  hsep(1,l,1)=psppar(1,l,1)
                  hsep(2,l,1)=0.0d0
                  hsep(3,l,1)=psppar(2,l,1)
                  hsep(4,l,1)=0.0d0
                  hsep(5,l,1)=0.0d0
                  hsep(6,l,1)=psppar(3,l,1)
c                 in the polarized or relativistic case,
c                 we assume all kij to be zero,
c                 i.e. same up and down projectors
                  if(nspin==2) hsep(:,l,2)=hsep(:,l,1)
               end do
            elseif (ipspcod == 3) then !HGH diagonal case
c              we need to compute the offdiagonal elements with the following coeffs
               ofdcoef(1,1)=-0.5d0*sqrt(3.d0/5.d0) !h2
               ofdcoef(2,1)=0.5d0*sqrt(5.d0/21.d0) !h4
               ofdcoef(3,1)=-0.5d0*sqrt(100.0d0/63.d0) !h5
            
               ofdcoef(1,2)=-0.5d0*sqrt(5.d0/7.d0) !h2
               ofdcoef(2,2)=1.d0/6.d0*sqrt(35.d0/11.d0) !h4
               ofdcoef(3,2)=-7.d0/3.d0*sqrt(1.d0/11.d0) !h5
            
               ofdcoef(1,3)=-0.5d0*sqrt(7.d0/9.d0) !h2
               ofdcoef(2,3)=0.5d0*sqrt(63.d0/143.d0) !h4
               ofdcoef(3,3)=-9.d0*sqrt(1.d0/143.d0) !h5
            
               ofdcoef(1,4)=0.0d0 !h2
               ofdcoef(2,4)=0.0d0 !h4
               ofdcoef(3,4)=0.0d0 !h5
            
c              this could possibly be done in a simpler way ...

               do l=1,lpx+1; do ispin=1,nspin
                  hsep(1,l,ispin)=psppar(1,l,ispin)
                  hsep(2,l,ispin)=psppar(2,l,ispin)*ofdcoef(1,l)
                  hsep(3,l,ispin)=psppar(2,l,ispin)
                  hsep(4,l,ispin)=psppar(3,l,ispin)*ofdcoef(2,l)
                  hsep(5,l,ispin)=psppar(3,l,ispin)*ofdcoef(3,l)
                  hsep(6,l,ispin)=psppar(3,l,ispin)
               end do; end do


c              in the nonrelativistic case, we are done.
               if(nspin==2)then
c                in the polarized case, copy the missing s projector
                 if(pol)  hsep(:,1,2)=hsep(:,1,1)
c                use psppar as a temporary array 
                 psppar=hsep
c                and then convert hij/kij to hij up/down 
                 do l=2,lpx+1
c                 l is index, angular momentum +one
c                 up
                  hsep(1,l,1)=psppar(1,l,1) +.5d0*(l-1)* psppar(1,l,2)!h11
                  hsep(2,l,1)=psppar(2,l,1) +.5d0*(l-1)* psppar(2,l,2)!h12
                  hsep(3,l,1)=psppar(3,l,1) +.5d0*(l-1)* psppar(3,l,2)!h22
                  hsep(4,l,1)=psppar(4,l,1) +.5d0*(l-1)* psppar(4,l,2)!h13
                  hsep(5,l,1)=psppar(5,l,1) +.5d0*(l-1)* psppar(5,l,2)!h23
                  hsep(6,l,1)=psppar(6,l,1) +.5d0*(l-1)* psppar(6,l,2)!h33
c                 down
                  hsep(1,l,2)=psppar(1,l,1) -.5d0* l   * psppar(1,l,2)!h11
                  hsep(2,l,2)=psppar(2,l,1) -.5d0* l   * psppar(2,l,2)!h12
                  hsep(3,l,2)=psppar(3,l,1) -.5d0* l   * psppar(3,l,2)!h22
                  hsep(4,l,2)=psppar(4,l,1) -.5d0* l   * psppar(4,l,2)!h13
                  hsep(5,l,2)=psppar(5,l,1) -.5d0* l   * psppar(5,l,2)!h23
                  hsep(6,l,2)=psppar(6,l,1) -.5d0* l   * psppar(6,l,2)!h33
                 end do
               end if
            end if

            if (ipspcod>9) then !HGH-K or HGH case,
c              psppar holds hij and kij in HGHK convention
c              fill hsep(up,dn) upper diagonal col by col, as needed for the fit

c             for a nonrelativistic calculation, discard the kij elements
              if(nspin==1)then
               do l=1,lpx+1
                  hsep(1,l,1)=psppar(1,l,1) !h11
                  hsep(2,l,1)=psppar(4,l,1) !h12
                  hsep(3,l,1)=psppar(2,l,1) !h22
                  hsep(4,l,1)=psppar(5,l,1) !h13
                  hsep(5,l,1)=psppar(6,l,1) !h23
                  hsep(6,l,1)=psppar(3,l,1) !h33
               end do
              else
c              relativistic or polarized calculation
               do l=1,lpx+1
c                 l is the index, angular momentum +one
c                 up
                  hsep(1,l,1)=psppar(1,l,1) +.5d0*(l-1)* psppar(1,l,2)!h11
                  hsep(2,l,1)=psppar(4,l,1) +.5d0*(l-1)* psppar(4,l,2)!h12
                  hsep(3,l,1)=psppar(2,l,1) +.5d0*(l-1)* psppar(2,l,2)!h22
                  hsep(4,l,1)=psppar(5,l,1) +.5d0*(l-1)* psppar(5,l,2)!h13
                  hsep(5,l,1)=psppar(6,l,1) +.5d0*(l-1)* psppar(6,l,2)!h23
                  hsep(6,l,1)=psppar(3,l,1) +.5d0*(l-1)* psppar(3,l,2)!h33
c                 if nspol==1 and l==1
c                 if(l==2-nspol)cycle
c                 down
                  hsep(1,l,2)=psppar(1,l,1) -.5d0* l   * psppar(1,l,2)!h11
                  hsep(2,l,2)=psppar(4,l,1) -.5d0* l   * psppar(4,l,2)!h12
                  hsep(3,l,2)=psppar(2,l,1) -.5d0* l   * psppar(2,l,2)!h22
                  hsep(4,l,2)=psppar(5,l,1) -.5d0* l   * psppar(5,l,2)!h13
                  hsep(5,l,2)=psppar(6,l,1) -.5d0* l   * psppar(6,l,2)!h23
                  hsep(6,l,2)=psppar(3,l,1) -.5d0* l   * psppar(3,l,2)!h33
               end do
              end if
            end if





c           output
            write(6,*)
            write(6,'(4f10.3,t30,a)') rcov,rprb,rcore,zcore,
     :           'rcov, rprb, rcore,zcore (NLCC)'
            if (ispp.eq.'r') then
               write(6,'(t30,a)')'relativistic calculation'
            elseif(ispp.eq.'s')then
               write(6,'(t30,a)')'spin polarized (dev) calculation'
            else
               write(6,'(t30,a)')'non relativistic calculation'
            endif
            write(6,'(t2,i10,t30,a)')iXC ,'XC-functional'
            write(6,*) 
            write(6,*) 'local part'
            write(6,'(f5.0,f7.2,f7.3,4(1pe11.3),t65,a)')
     :           znuc,zion,rloc,gpot(1),gpot(2),gpot(3),gpot(4),
     :           'znuc,zion, rloc, gpot() '
            if (lpx.ge.0) then
               write(6,*) 'nonlocal part in internal format'
               write(6,'(i4,t60,a)') lpx ,
     :              'lpx, (Projectors for l=0..lpx)'
               do l=0,lpx
                  write(6,*) 'l=',l
                  write(6,'(f7.3,t8,6(1pe11.3),t76,a)') r_l(l+1),
     :                 (hsep(i,l+1,1),i=1,6),'r_l(),hsep(), '//is(1)
                  if (l.gt.0 .and. nspin.eq.2)
     :                 write(6,'(t8,6(1pe11.3),t76,a)')
     :                 (hsep(i,l+1,2),i=1,6),'      hsep(), '//is(2)
               enddo
            endif
            if(zcore>0)then
               write(6,*)'This pseudopotential inplies a'
               write(6,*)'nonlinear core correction.'
               if(rcore>0d0)then
                 write(6,*)'rhocore is a guassian.'
               else
                 write(6,*)'rhocore is an exponential form.'
               end if
               write(6,*)'Parameters rcore and zcore:',rcore,zcore
            end if
c           !done 
         endif ! first iteration

c
c     previous versions have read weights.par at each iteration step.
c

c     let us sacrifice this feature for the sake of stability in
c     parallel runs.


c

c     calc. exponents of gaussians
         a0=rloc/rij
c     take this for an crude initial fit:
c     tt=2.d0**.4d0
         if (denbas) then
c     fine fit:
            tt=sqrt(sqrt(2.d0))
         else
c     normal fit:
            tt=2.d0**.3d0
         endif
         a=a0
         do i=0,ng
c           
c           a=a0*tt**i
            xp(i)=.5d0/a**2
            a=a*tt
         enddo
         write(6,*)
         write(6,*)'Gaussian basis'
         write(6,*)'______________'
         write(6,*)
         write(6,'(a,4(1pe11.4))') ' amin,amax',a0,a
         write(6,'(a,t10,3(1pe11.4),a,2(1pe11.4))') ' gaussians ',
     &        xp(1),xp(2),xp(3),' .... ',xp(ng-1),xp(ng)
         write(6,*)'gaussians:',ng
c     set up radial grid
         nint=5*(ng+14)
         rmax=min(15.d0*rprb,120.d0)
         a_grd=a0/400.d0
         b_grd=log(rmax/a_grd)/(nint-2)
         call radgrid(nint,rr,rw,rd,a_grd,b_grd,rmax)
         write(6,'(a,t10,3(1pe11.4),a,2(1pe11.4))') ' r-grid: ',
     &        rr(1),rr(2),rr(3),' .... ',rr(nint-1),rr(nint)
         write(6,*)'gridpoints:',nint
         write(6,*)
         call crtvh(ng,lcx,lmax,xp,vh,nint,rmt,rmtg,ud,rr)
         write(6,*)


         if (namoeb.gt.0) then

c     some OLD, commmented out feature
c     refine simplex only every 10.th step
c     start of if-block
c          if (mod(iter,10).eq.0) then



         write(6,*)'Reading fitting parameters from FITPAR'
         write(6,*)'______________________________________'
         write(6,*)

         open(11,file='FITPAR')
c        convention: a comment line first
         read(11,*)
c        New optional input line from FITPAR's header:

c        initial simplex width (dh), max successive steps
c        with no improvement (ntrymax) and convergence
c        criterion for the simplex (FTOL)
         read(11,*,iostat=ierr)dh,ntrymax,FTOL
         if(ierr.eq.0)then
c            write(6,*)'Optional fitting paramaters from'//
c     :                ' the first line of FITPAR:'
            write(6,'(1x,a,1pe20.10,1x,i3,1x,g15.5)')
     :           'initial width, tries and spread for convergence ',
     :                          dh,ntrymax,FTOL
         else
            dh=0.2d0
            ntrymax=200
            ftol=1d-7
         end if
         close(11)
         

c
c     pack initial guess
c
            write(6,*)
            write(6,*)'Free parameters of the fit:'
            write(6,*)
c           if(ipspcod.le.9)lpx=-2
            verbose=.true.
            call  ppack (verbose,rloc,gpot,hsep,r_l,pp(1),
     :           lpx,lpmx,nspin,pol,nsmx,maxdim,nfit,'init',
     :           avgl1,avgl2,avgl3,ortprj,litprj,
     :           rcore,zcore,znuc,zion)
            verbose=.false.
c
c     initial simplex
c

c           copy the packed initial vertex to the other vertices
c           does not make much sense, though, as we copy zeroes only.
c           keep it in case ppack.f is changed
            do i=1,nfit
               do j=1,nfit
                  pp(j+i*nfit)=pp(j)
               enddo
            enddo

c           This width is not hard coded anymore
c           dh=0.2d0

c           shift vertices number 2 to nfit+1 by random numbers in [-dh/2:dh/2]
c           in this convention, only move the k-th component of vertex k+1

c           the random numbers generated SHOULD be the same for all processes.
c           Though, let us enforce equality with MPI broadcast from process 0.
            if(irpoc==0)then
              do i=1,nfit
c     f90 intrinsic
                call random_number(randNr)
                pp(i+i*nfit)=pp(i)+dh*1.0d0*(randNr-.5d0)
c     Intel (ifc)
c              pp(i+i*nfit)=pp(i)+dh*1.0d0*(dble(rand(0.0d0))-.5d0)
c     IBM/DEC/PGI
c              pp(i+i*nfit)=pp(i)+dh*1.0d0*(dble(rand())-.5d0)
c     CRAY
c             pp(i+i*nfit)=pp(i)+dh*1.0d0*(ranf()-.5d0)
              enddo
            end if
            call MPI_BCAST(pp,nfit*(nfit+1),MPI_DOUBLE_PRECISION,0,
     :                                      MPI_COMM_WORLD,ierr)

            write(6,*)
            write(6,'(1x,a,i4)')'starting amoeba cycle',iiter
            write(6,*)          '_________________________'
            write(6,*)
             write(6,'(2a)')' Penalty contributions of the initial ',
     :                      'parameters and random simplex vertices'
             write(6,'(2a)')' _______________________________________',
     :                      '____________________________________'
             write(6,*)
             if(nproc>1)then
                write(6,'(a)')'  amoeba   vertex    overall penalty'//
c    :               '      softness             narrow radii      '//
     :               '      softness             psp empirical     '//
     :               '   all configurations   all excitation       '//
     :               'this configuration'
                write(6,'(a)')'    step   number    value          '//
     :               '      gauss-wvlt           (hgridmax/r)**12  '//
     :               '   E_KS, integrals      energies             '//
     :               'and its excitation'
             else
                write(6,'(a)')'    step    vertex   penalty,       '//
c    :               '      gaus-wvlt            narrow radii      '//
     :               '      gaus-wvlt            psp empirical     '//
     :               '   eigenvalues etc'
             end if
c           do i=2,nfit+1  would be the random vertices, 1 has no shift
            iter=0
            do i=1,nfit+1
               write(99,*)'vertex',i
               call penalty(energ,verbose,nfit,pp(1+(i-1)*nfit),yp(i),
     :           noccmax,noccmx,lmax,lmx,lpx,lpmx,lcx,nspin,pol,nsmx,
     :           no,lo,so,ev,crcov,dcrcov,ddcrcov,norb,
     :           occup,aeval,chrg,dhrg,ehrg,res,wght,
     :           wfnode,psir0,wghtp0,
     :           rcov,rprb,rcore,zcore,znuc,zion,rloc,gpot,r_l,hsep,
     :           vh,xp,rmt,rmtg,ud,nint,ng,ngmx,psi,
     :           avgl1,avgl2,avgl3,ortprj,litprj,igrad,rr,rw,rd,
c                the following lines differ from pseudo2.2
     :           iproc,nproc,wghtconf,wghtexci,wghtsoft,wghtrad,wghthij,
     :           nhgrid, hgridmin,hgridmax, nhpow,ampl,crmult,frmult,
     :           excitAE,ntime,iter,itertot,penref,time)
c              write(6,'(4x,i3,4x,1pe20.10)')i,yp(i)
CMK            CALL FLUSH(6)
            enddo
            write(99,*)'history of lowest vertices'


             write(6,*)
             write(6,'(2a)')' Penalty contributions of the',
     :                      ' currently lowest vertex'
             write(6,'(2a)')' _______________________________________',
     :                      '_____________'
             write(6,*)
             if(nproc>1)then
                write(6,'(a)')'  amoeba   gatom     overall penalty'//
c    :               '      softness             narrow radii      '//
     :               '      softness             psp empirical     '//
     :               '   all configurations   all excitation       '//
     :               'this configuration'
                write(6,'(a)')'    step   calls     value          '//
     :               '      gauss-wvlt           (hgridmax/r)**12  '//
     :               '   E_KS, integrals      energies             '//
     :               'and its excitation'
             else
                write(6,'(a)')'    step    gatom    penalty,       '//
c    :               '      gaus-wvlt            narrow radii      '//
     :               '      gaus-wvlt            psp empirical     '//
     :               '   eigenvalues etc'
             end if



c            The vertex number 1 holds the initial psp

c            this call is not be needed anymore, because
c            we already wrote details of all vertices

c            iter=0
c            write(*,*)'test line, vertex 1 again'
c            call penalty(energ,verbose,nfit,pp(1),yp(1),
c    :          noccmax,noccmx,lmax,lmx,lpx,lpmx,lcx,nspin,pol,nsmx,
c    :          no,lo,so,ev,crcov,dcrcov,ddcrcov,norb,
c    :          occup,aeval,chrg,dhrg,ehrg,res,wght,
c    :          wfnode,psir0,wghtp0,
c    :          rcov,rprb,rcore,zcore,znuc,zion,rloc,gpot,r_l,hsep,
c    :          vh,xp,rmt,rmtg,ud,nint,ng,ngmx,psi,
c    :          avgl1,avgl2,avgl3,ortprj,litprj,igrad,rr,rw,rd,
c               the following lines differ from pseudo2.2
c    :          iproc,nproc,wghtconf,wghtexci,wghtsoft,wghtrad,wghthij,
c    :          nhgrid, hgridmin,hgridmax, nhpow,ampl,crmult,frmult,
c    :          excitAE,ntime,iter,itertot,penref,time)


c     refine simplex only every 10.th step
c     end of if-block
c         endif
c
c     starting amoeba
c
c           ftol=1.d-7 is not hardcoded anymore
            call AMOEBA(pp,yp,nfit,FTOL,ITER,nsmplx,namoeb,
     :           noccmax,noccmx,lmax,lmx,lpx,lpmx,lcx,nspin,pol,nsmx,
     :           no,lo,so,ev,crcov,dcrcov,ddcrcov,norb,
     :           occup,aeval,chrg,dhrg,ehrg,res,wght,
     :           wfnode,psir0,wghtp0,
     :           rcov,rprb,rcore,zcore,znuc,zion,rloc,gpot,r_l,hsep,
     :           vh,xp,rmt,rmtg,ud,nint,ng,ngmx,psi,
     :           avgl1,avgl2,avgl3,ortprj,litprj,igrad,rr,rw,rd,
c                the following line differs from pseudo2.2
     :           iproc,nproc,wghtconf,wghtexci,wghtsoft,wghtrad,wghthij,
     :           nhgrid,hgridmin,hgridmax,nhpow,ampl,crmult,frmult,
     :           ntrymax,excitAE,ntime,itertot,energ,verbose,time)
c           write(6,*) 'Finished amoeba with ',iter,'iterations'
         else
            verbose=.false.
            energ=.true.
c           if(ipspcod.le.9)lpx=-2
               call ppack (verbose,rloc,gpot,hsep,r_l,pp(1),
     :              lpx,lpmx,nspin,pol,nsmx,maxdim,nfit,'init',
     :              avgl1,avgl2,avgl3,ortprj,litprj,
     :              rcore,zcore,znuc,zion)

            if (nconfpaw/=-1) then
               verbose=.true.

               call pawpatch(energ,verbose,nfit,pp(1),yp(1),
     :              noccmax,noccmx,lmax,lmx,lpx,lpmx,lcx,nspin,pol,nsmx,
     :              no,lo,so,ev,crcov,dcrcov,ddcrcov,norb,
     :              occup,aeval,chrg,dhrg,ehrg,res,wght,
     :              wfnode,psir0,wghtp0,
     :              rcov,rprb,rcore,zcore,znuc,zion,rloc,gpot,r_l,hsep,
     :              vh,xp,rmt,rmtg,ud,nint,ng,ngmx,psi,
     :              avgl1,avgl2,avgl3,ortprj,litprj,igrad,rae,
c                the following lines differ from pseudo2.2
     :              iproc,nproc,wghtconf,wghtexci,wghtsoft,wghtrad,
     :          wghthij,
     :              nhgrid,hgridmin,hgridmax,nhpow,ampl,crmult,frmult,
     :              excitAE,ntime,iter,itertot,penref,time,ngrid, 
     :               nconfpaw, npawl, nchannelspaw , ispp, pawstatom,
     :              pawstN, pawstL  , pawstP,     pawrcovfact    )


            else
c              No fitting was requested, evaluate penalty contributions once         
c              call penalty with the verbose flag to print the details
               verbose=.true.
               call penalty(energ,verbose,nfit,pp(1),yp(1),
     :              noccmax,noccmx,lmax,lmx,lpx,lpmx,lcx,nspin,pol,nsmx,
     :              no,lo,so,ev,crcov,dcrcov,ddcrcov,norb,
     :              occup,aeval,chrg,dhrg,ehrg,res,wght,
     :              wfnode,psir0,wghtp0,
     :              rcov,rprb,rcore,zcore,znuc,zion,rloc,gpot,r_l,hsep,
     :              vh,xp,rmt,rmtg,ud,nint,ng,ngmx,psi,
     :              avgl1,avgl2,avgl3,ortprj,litprj,igrad,rr,rw,rd,
c                the following lines differ from pseudo2.2
     :              iproc,nproc,wghtconf,wghtexci,wghtsoft,wghtrad,wghthij,
     :              nhgrid,hgridmin,hgridmax,nhpow,ampl,crmult,frmult,
     :              excitAE,ntime,iter,itertot,penref,time)
               write(6,*)'Overall penalty value would be',yp(1)

c          IMPORTANT TEST: Calling gatom once more should not make any difference 
c     call gatom(nspol,energ,verbose,
c    :     noccmax,noccmx,lmax,lmx,lpx,lpmx,lcx,nspin,nsmx,
c    :     occup,aeval,chrg,dhrg,ehrg,res,wght,wfnode,psir0,
c    :     rcov,rprb,rcore,zcore,znuc,zion,rloc,gpot,r_l,hsep,
c    :     vh,xp,rmt,rmtg,ud,nint,ng,ngmx,psi,igrad,
c    :     rr,rw,rd,ntime,itertot,etotal)

c          ...true

            endif

          endif
c         fit or evaluate once


c
c     print results
c
         write(6,*)
         write(6,*)'Penalty contributions from this configuration'
         write(6,*)'_____________________________________________'
         write(6,*)
         write(6,'(2(tr10,a,1pe12.4))')
     :        'psi(r=0) =',psir0,'; psi(0)*wght=',abs(psir0*wghtp0)
         write(6,*)
         write(6,'(a,t32,a,t42,a,t55,a,t64,a)')
     :        ' nl    s      occ','ae','pseudo','diff','diff*weight'

         do iorb=1,norb
            write(6,31) noae(iorb),il(lo(iorb)+1),so(iorb),zo(iorb)
 31         format(1x,i1,a1,f6.1,f10.4)
            nocc=no(iorb)
            l=lo(iorb)
            ispin=1
            if (so(iorb).lt.0) ispin=2
            write(6,32) 'eigenvalue ',
     :           ev(iorb),aeval(nocc,l+1,ispin),
     :           aeval(nocc,l+1,ispin)-ev(iorb),
     :           abs(wght(nocc,l+1,ispin,1)*
     :           (aeval(nocc,l+1,ispin)-ev(iorb)))
            write(6,32) 'charge    ',
     :           crcov(iorb),chrg(nocc,l+1,ispin),
     :           chrg(nocc,l+1,ispin)-crcov(iorb),
     :           abs(wght(nocc,l+1,ispin,2)*
     :           (chrg(nocc,l+1,ispin)-crcov(iorb)))
            if (wght(nocc,l+1,ispin,3).ne.0.0d0)
     :           write(6,32) 'dcharge   ',
     :           dcrcov(iorb),dhrg(nocc,l+1,ispin),
     :           100.d0*abs(1.d0-dhrg(nocc,l+1,ispin)/dcrcov(iorb)),
     :           abs(wght(nocc,l+1,ispin,3))*
     :           100.d0*abs(1.d0-dhrg(nocc,l+1,ispin)/dcrcov(iorb))
            if (wght(nocc,l+1,ispin,4).ne.0.0d0)
     :           write(6,32) 'echarge   ',
     :           ddcrcov(iorb),ehrg(nocc,l+1,ispin),
     :           100.d0*abs(1.d0-ehrg(nocc,l+1,ispin)/ddcrcov(iorb)),
     :           abs(wght(nocc,l+1,ispin,4))*
     :           100.d0*abs(1.d0-ehrg(nocc,l+1,ispin)/ddcrcov(iorb))
            write(6,33) 'residue   ',
     :           res(nocc,l+1,ispin),
     :           abs(wght(nocc,l+1,ispin,5)*res(nocc,l+1,ispin))
            if (wght(nocc,l+1,ispin,6).ne.0.0d0)
     :           write(6,33) 'rnode     ',
     :           wfnode(nocc,l+1,ispin,1),
     :           abs(wght(nocc,l+1,ispin,6)*wfnode(nocc,l+1,ispin,1))
            if (wght(nocc,l+1,ispin,7).ne.0.0d0)
     :           write(6,33) 'dnode     ',
     :           wfnode(nocc,l+1,ispin,2),
     :           abs(wght(nocc,l+1,ispin,7)*wfnode(nocc,l+1,ispin,2))
            if (wght(nocc,l+1,ispin,8).ne.0.0d0)
     :           write(6,33) 'ddnode    ',
     :           wfnode(nocc,l+1,ispin,3),
     :           abs(wght(nocc,l+1,ispin,8)*wfnode(nocc,l+1,ispin,3))
            write(6,*)
         enddo
 32      format (t10,a,t25,4(1pe12.4))
 33      format (t10,a,t25,2(1pe24.4))
         write(6,*) 'diff for dcharg and echarge is given in (%)'

         if(namoeb.ne.0)then
           write(6,*)
           write(6,*) 'Resulting Pseudpotential Parameters'
           write(6,*) '___________________________________'
           write(6,*)

c     At this point, overwrite  the psppar in ABINIT pspcod=10 format

c     dismissed:
c          Spin polarized runs will use an experimental format pspcod=11
c          with two spin channels also for s projectors.

          if(iproc==0)then
          open(unit=13,file='psppar')!.out'),position='append')

          if(ispp=='r')then
            write(13,'(a)',advance='no')
     :      'relativistic '
          elseif(ispp=='s')then
            write(13,'(a)',advance='no')
     :      'spin-polarized '
          else
            write(13,'(a)',advance='no')
     :      'nonrelativistic '
          end if
          write(13,'(i3,3(1x,g9.3),a)')ng,rij,rcov,rprb,
     :                        ' spin treatment, ng rij rcov and rprb'

c         this intrinsic may have different conventions on its arguments
          call date_and_time(dateYMD)
          write(13,'(2i5,2x,2a)') int(znuc+.1),int(zion+.1),dateYMD,
     :                     ' zatom, zion, date (yymmdd)'
          write( 6,'(2i5,2x,2a)') int(znuc+.1),int(zion+.1),dateYMD,
     :                     ' zatom, zion, date (yymmdd)'
c         use pspcod=10 per default, always. When we tried to write hij
c         and kij also for s projectors in a polarized fit, we used a
c         special format pspcod 11= 9+nspol, which does not make much
c         sense and has been discarded.
          write(13,'(3x,i3,i8,i2,a)') 10,ixc,lpx,
     :               ' 0 2002 0     pspcod,IXC,lmax,lloc,mmax,r2well'
          write( 6,'(3x,i3,i7,i2,a)') 10,ixc,lpx,
     :               ' 0 2002 0     pspcod,IXC,lmax,lloc,mmax,r2well'
          ngpot=0
          do j=1,4
             if (gpot(j).ne.0.d0) ngpot=j
          enddo
          if(ngpot==0)then
            write(13,'(1pe18.8,a)') rloc,' 0 rloc nloc ck (none)'
            write( 6,'(1pe18.8,a)') rloc,' 0 rloc nloc ck (none)'
          else
            write(label,'(a,i1.1,a)')'(1pe18.8,i2,',ngpot,'(1pe18.8),a)'
            write(13,label) rloc,ngpot,gpot(1:ngpot),
     :             ' rloc nloc c1 .. cnloc'
            write( 6,label) rloc,ngpot,gpot(1:ngpot),
     :             ' rloc nloc c1 .. cnloc'
          end if
          write(13,'(i3,2g16.8,a)')lpx+1, rcore,zcore,
     :                             'nsep, rcore&zcore (NLCC)'
          write( 6,'(i3,2g16.8,a)')lpx+1, rcore,zcore,
     :                             'nsep, rcore&zcore (NLCC)'
          end if

           if (lpx.ge.0) then
              if((pol.or.nspin==1).and.iproc==0)then
                 do l=0,lpx
                    lpj=1
                    if(abs(hsep(3,l+1,1))>1d-8)lpj=2
                    if(abs(hsep(6,l+1,1))>1d-8)lpj=3
                    if(lpj==1)then
                        write(13,'(2x,g16.8,i3,1pe16.8,6x,a)')
     :                            r_l(l+1),lpj, hsep(1,l+1,1)
     :                            ,il(l+1)//'-projector'
                        write( 6,'(2x,g16.8,i3,1pe16.8,6x,a)')
     :                            r_l(l+1),lpj, hsep(1,l+1,1)
     :                            ,il(l+1)//'-projector'
                    elseif(lpj==2)then
                        write(13,'(2x,g16.8,i3,2g16.8,6x,a)')
     :                            r_l(l+1),lpj, hsep(1:2,l+1,1)
     :                            ,il(l+1)//'-projector'
                        write( 6,'(2x,g16.8,i3,2g16.8,6x,a)')
     :                            r_l(l+1),lpj, hsep(1:2,l+1,1)
     :                            ,il(l+1)//'-projector'
                        write(13,'(37x,g16.8)') hsep(3,l+1,1)
                        write( 6,'(37x,g16.8)') hsep(3,l+1,1)
                    elseif(lpj==3)then
                        write(13,'(2x,g16.8,i3,3g16.8,6x,a)')
     :                      r_l(l+1),lpj, hsep(1:2,l+1,1), hsep(4,l+1,1)
     :                            ,il(l+1)//'-projector'
                        write( 6,'(2x,g16.8,i3,3g16.8,6x,a)')
     :                      r_l(l+1),lpj, hsep(1:2,l+1,1), hsep(4,l+1,1)
     :                            ,il(l+1)//'-projector'
                        write(13,'(37x,2g16.8)')
     :                      hsep(3,l+1,1),hsep(5,l+1,1)
                        write( 6,'(37x,2g16.8)')
     :                      hsep(3,l+1,1),hsep(5,l+1,1)
                        write(13,'(53x,g16.8)') hsep(6,l+1,1)
                        write( 6,'(53x,g16.8)') hsep(6,l+1,1)
                    end if
                    if(l>0.and.lpj==1) write(13,'(23x,a)')'0 (no kij)'
                    if(l>0.and.lpj==2) then 
                         write(13,'(23x,a)')'0  0 (no kij)'
                         write( 6,'(23x,a)')'0  0 (no kij)'
                         write(13,'(23x,a)')'   0'
                         write( 6,'(23x,a)')'   0'
                    end if
                    if(l>0.and.lpj==3) then 
                         write(13,'(23x,a)')'0  0  0 (no kij)'
                         write( 6,'(23x,a)')'0  0  0 (no kij)'
                         write(13,'(23x,a)')'   0  0'
                         write( 6,'(23x,a)')'   0  0'
                    end if
                 enddo
              else
c                output in the relativistic case:
c                kij terms need to be computed
                 do l=0,lpx
                    do i=1,6
                       if (l.gt.0) then
                          havg(i)=((l+1)*hsep(i,l+1,1)+l*hsep(i,l+1,2))
     :                         /(2*l+1)
                          hso(i)=2*(hsep(i,l+1,1)-hsep(i,l+1,2))
     :                         /(2*l+1)
                       else
                          havg(i)=hsep(i,l+1,1)
                       endif
                    enddo
c                   only process 0 writes out psppar.out
                    if(iproc>0) cycle
c                   get the matrix dimension lpj = 1,2 or 3
                    lpj=1
                    if(max(abs(havg(3)),abs(hso(3)))>1d-8)lpj=2
                    if(max(abs(havg(6)),abs(hso(6)))>1d-8)lpj=3
                    if(lpj==1)then
                        write(13,'(2x,1pe16.8,i3,1pe16.8,6x,a)')
     :                             r_l(l+1),lpj, havg(1)
     :                            ,il(l+1)//'-projector'
                        write( 6,'(2x,1pe16.8,i3,1pe16.8,6x,a)')
     :                             r_l(l+1),lpj, havg(1)
     :                            ,il(l+1)//'-projector'
                        if (l.gt.0)write(13,'(21x,1pe16.8)')hso(1)
                        if (l.gt.0)write( 6,'(21x,1pe16.8)')hso(1)
                    elseif(lpj==2)then
                        write(13,'(2x,1pe16.8,i3,2(1pe16.8),6x,a)')
     :                             r_l(l+1),lpj, havg(1:2)
     :                            ,il(l+1)//'-projector'
                        write( 6,'(2x,1pe16.8,i3,2(1pe16.8),6x,a)')
     :                             r_l(l+1),lpj, havg(1:2)
     :                            ,il(l+1)//'-projector'
                        write(13,'(37x,1pe16.8)') havg(3)
                        write( 6,'(37x,1pe16.8)') havg(3)
                        if (l.gt.0)then
                          write(13,'(21x,2(1pe16.8))') hso(1:2)
                          write( 6,'(21x,2(1pe16.8))') hso(1:2)
                          write(13,'(37x,1pe16.8)') hso(3)
                          write( 6,'(37x, 1pe16.8)') hso(3)
                        end if
                    elseif(lpj==3)then
                        write(13,'(2x,1pe16.8,i3,3(1pe16.8),6x,a)')
     :                             r_l(l+1),lpj,havg(1:2), havg(4)
     :                            ,il(l+1)//'-projector'
                        write( 6,'(2x,1pe16.8,i3,3(1pe16.8),6x,a)')
     :                             r_l(l+1),lpj,havg(1:2), havg(4)
     :                            ,il(l+1)//'-projector'
                        write(13,'(37x,2(1pe16.8))') havg(3),havg(5)
                        write( 6,'(37x,2(1pe16.8))') havg(3),havg(5)
                        write(13,'(53x,1pe16.8)') havg(6)
                        write( 6,'(53x,1pe16.8)') havg(6)
                        if (l.gt.0)then
                          write(13,'(21x,3(1pe16.8))') hso(1:2),hso(4)
                          write( 6,'(21x,3(1pe16.8))') hso(1:2),hso(4)
                          write(13,'(37x,2(1pe16.8))') hso(3),hso(5)
                          write( 6,'(37x,2(1pe16.8))') hso(3),hso(5)
                          write(13,'(53x, 1pe16.8)') hso(6)
                          write( 6,'(53x, 1pe16.8)') hso(6)
                        end if
                    end if
                 enddo
              endif
           endif ! there is a nonlocal part
         end if
         if(iproc==0) close(13)
c        end of writing the psppar
         
c        nonlinear core correction file for BigDFT
         if(zcore>0d0.and.iproc==0)then
            open(13,file='nlcc')
            write(13,*)'0  no valence distribution to subtract'
            write(13,*)'1  one single gaussian for the core'
c           problem: analytic volume integral of a guassian should give
c                 zcore = c * (sqrt(2pi)*sigma)**3
c                 where g(x) = c * exp(-0.5*(|x|/sigma)**2)
c           but BigDFT includes a factor of 4pi for rhocore, thus 
c                 zcore = c *sigma**3 * sqrt(pi/2)
            tt=sqrt(2d0*atan(1d0))
            write(13,*)rcore, zcore/(tt*rcore**3),' 0 0 0  no polyn'
c           write(13,*)c(1,1) = zcore/(sqrt(2pi)*rcore)**3 else c = 0  
            close(13)
c           and a plot to compare with the full AE core charge
            open(13,file='nlcc.gnuplt')
            write(13,*)'rcore=',abs(rcore)
            write(13,*)'zcore=',zcore
            if(rcore>0d0)
     :      write(13,*)'rho(x)=zcore/(sqrt(2*pi)*rcore)**3'
     :               //'*exp(-0.5*(x/rcore)**2)'
            if(rcore<0d0)
     :      write(13,*)'rho(x)=zcore/(4*pi*rcore**3)'
     :               //'*(1+2*x/rcore)*exp(-2*(x/rcore))'
 
            write(13,*)"dS(x)=4*pi*x*x"
            write(13,*)"  p rho(x)*dS(x)"
            write(13,*)"rep 'ae.core.dens.plt'"
            write(13,*)"rep 'ae.val.dens.plt'"
            write(13,*)"show function"
            close(13)
         end if





c write out the local potential
        if(iproc==0)then
        open(17,file='local.pot')
        write(17,'(a)')'# r,   Vloc+z/r,  gpot term, Vprb(r)'
           do k=1,nint
              r=rr(k)
              dd=exp(-.5d0*(r/rloc)**2)*
     :                  (gpot(1) + gpot(2)*(r/rloc)**2+  
     :                   gpot(3)*(r/rloc)**4 +
     :                   gpot(4)*(r/rloc)**6 )

              write(17,'(4(1pe18.8))') r,-zion*Derf(r/(sqrt(2.d0)*rloc))/r
     1                    +dd+zion/r, dd,.5d0*(r/rprb**2)**2
           enddo
        end if


         if(ldump.and.iproc==0)then
            write(6,*)'Writing out a dumpfile of',
     :        8*4+size(xp)+size(psi)+size(occup),'byte'           
            open(13,file='dumpfile.bin',form='unformatted')
c           open(13,file='dumpfile')
            write(13)ng,noccmax,lmax,lpx,
     :               lcx,nspin,xp,psi,occup
            close(13)
         end if
c
c     here we used to overwrite old values of 'psp.par' with the current ones
c

c     for now, only used with the -info flag

         if (iproc==0.and.namoeb.ne.0) then
            if (info) then
               open(unit=23,file='gauss.par',form='formatted')
               write(23,*) 'Additional information (last calculation):'
               write(23,*) 'gaussian exponents:'
               write(23,*) 'xp():',(xp(i),i=1,ng)
               write(23,*) 'orbital coefficients ',
     :              '(one orbital per column)'
               do l=0,lmax
c        no spin down s orbitals in the r case: nspin=2, nspol=1, 2l+1=1
                  do ispin=1,max(nspol,min(2*l+1,nspin))
                     write(23,*) 'l,ispin:',l,ispin
                     write(23,*) 'psi n=1,noccmax(l):'
                     do i=0,ng
                        write(23,'(10(1pe20.10))')
     :                       (psi(i,nocc,l+1,ispin),
     :                       nocc=1,nomax(l))
                     enddo
                  enddo
               enddo
            endif
            close(unit=23)
         endif


c
c     PLOT WAVEFUNCTIONS (up to 5*rcov)   c.hartwigsen: version for gnuplot
c

c     NOTE: For now, only do this with config0 from atom.00.ae
c           should be straight forward to write out all plots in parallel


         if (plotwf.and.iproc==0) then
            call detnp(ngrid,rae,5*rcov,np)
            open(32,file='pswf.gnu',form='formatted',status='unknown')
            write (32,*) 'set data style lines'
            do iorb=1,norb
               nocc=no(iorb)
               l=lo(iorb)
               ispin=1
               if(so(iorb).lt.0) ispin=2
               if (ispp.eq.'r')then
                  fname= 'ps.'//char(ichar('0')+noae(iorb))
     :                 //il(lo(iorb)+1)
     :                 //char(ichar('0')+int(2*(lo(iorb)+so(iorb))))
     :                 //'by2.dat'
                  tname=char(ichar('0')+noae(iorb))//il(lo(iorb)+1)
     :                 //char(ichar('0')+int(2*(lo(iorb)+so(iorb))))
     :                 //'/2'
               else
c                 either nonrel or spin pol
                  fname= 'ps.'//char(ichar('0')+noae(iorb))
     :                 //il(lo(iorb)+1)
                  tname=char(ichar('0')+noae(iorb))//il(lo(iorb)+1)
                  if(so(iorb)<0)then
                      fname=trim(fname)//'.down'
                      tname=trim(tname)//'.down'
                  elseif(so(iorb)>0)then
                      fname=trim(fname)//'.up'
                      tname=trim(tname)//'.up'
                  end if
                  fname=trim(fname)//'.dat'
               endif
               open(unit=2,file=trim(fname),
     :                 form='formatted',status='unknown')
c     find outer max of psi (approx), search from 10 bohr down
                  ra=10.d0
                  ttrmax=ra
                  ttmax= dabs(wave(ng,l,xp,psi(0,nocc,l+1,ispin),ra))
                  do i=100,0, -1
                     ra= 0.1d0 * i
                     ttpsi=dabs(wave(ng,l,xp,psi(0,nocc,l+1,ispin),ra))
c                     print*,ra,ttpsi
                     if ( ttpsi .gt. ttmax
     :                    .and. ttpsi .gt. 1.0d-4 ) then
                        ttmax=ttpsi
                        ttrmax=ra
                     endif
                  if (ttpsi.lt.ttmax .and. ttpsi.gt.1.0d-4) goto 3456
                  enddo
 3456             continue
c     ae/pseudo wfs should have the same sign for large r when plotted
               call detnp(ngrid,rae,ttrmax,nsign)
c     never use first gridpoint! (only relevant for H and He)
               if (nsign.eq.1) nsign=nsign+1
               tt=psiold(nsign,nocc,l+1,ispin)
               sign1=tt/abs(tt)
               tt= wave(ng,l,xp,psi(0,nocc,l+1,ispin),rae(nsign))
               sign2=tt/abs(tt)
               do i=2,np
                  ttold=psiold(i,nocc,l+1,ispin)*sign1*rae(i)
                  ttold=max(min(3.d0,ttold),-3.d0)
                  tt=wave(ng,l,xp,psi(0,nocc,l+1,ispin),rae(i))
     :                 *sign2*rae(i)
                  tt=max(min(3.d0,tt),-3.d0)
                  ttdiff=psiold(i,nocc,l+1,ispin)*sign1-
     :                 wave(ng,l,xp,psi(0,nocc,l+1,ispin),rae(i))*sign2
                  ttdiff= ttdiff*rae(i)
                  ttdiff=log(max(abs(ttdiff),1.d-8))/log(10.d0)
                  write(2,'(7e20.10)') rae(i),ttold,tt,ttdiff
c     plot of the wavefunction and the higher derivatives
c                  write(2,'(7e20.10)') rae(i),
c     :                 psiold(i,nocc,l+1,ispin),
c     :                 wave(ng,l,xp,psi(0,nocc,l+1,ispin),rae(i)),
c     :                 dwave(ng,l,xp,psi(0,nocc,l+1,ispin),rae(i)),
c     :                 ddwave(ng,l,xp,psi(0,nocc,l+1,ispin),rae(i))
               enddo
               close(2)
               write (32,*)'  plot "'//fname(1:index(fname,' ')-1)
     :              //'"     title "'//tname(1:index(tname,' ')-1)
     :              //'"'
               write (32,*)'replot "'//fname(1:index(fname,' ')-1)
     :              //'" using 1:3 title "pseudo"'
               write(32,*) 'pause -10 "Hit return to continue"'
               write (32,*)'replot "'//fname(1:index(fname,' ')-1)
     :              //'" using 1:4 title "diff"'
               write(32,*) 'pause -10 "Hit return to continue"'
            enddo
            write (32,*) 'set nokey'
            close(unit=32)
            write(6,*)
            write(6,*)'Plot files for gnuplot written.'
            write(6,*)'To plot wfs type ''gnuplot pswf.gnu'' '
            write(6,*)
         endif



c    -----------------------------------------------
c                     MAIN LOOP END
      if(iiter>namoeb)exit
      enddo
 1000 continue


c     test orthogonality of the projectors
c
      if (lpx.ge.0.and.ortprj)
     :     call pj2test(hsep,lpx,lpmx,lmx,nspin,nsmx,r_l,is)
c
      write(6,*)
      write(6,*) 'Total SCF-cycles:',itertot
      write(6,*) 'Pseudoatom calculations:',ntime
      call MPI_FINALIZE(ierr)
      call libxc_functionals_end()

      call cpu_time(t)
      write(6,*)
      write(6,*)'Simple Time Profiling'
      write(6,*)'_____________________'
      write(6,*)
      write(6,'(a,f9.3,a)')' gatom  overall runtime ',time(1),' seconds'
      write(6,'(a,f9.3,a)')' MPI comminactions time ',time(2),' seconds'
      write(6,'(a,f9.3,a)')' wavelet libraries time ',time(3),' seconds'
      write(6,'(a)')' ________________________________________'
      write(6,'(a,f9.3,a)')'               CPU time ',t-t0,' seconds'
      write(6,*)       
      write(6,*)'______________________________________________________'
      write(6,*)'                                              finished'
      write(6,*)       
c     end


c     optional info when call to bash is available
      call system('echo -n " process finished " >> my.hostnames')
      end program

c
c     CRAY: no derf() -> user erf()
c
c      real*8 function derf(x)
c      REAL*8 X
c      DERF=ERF(X)
c      RETURN
c      END


