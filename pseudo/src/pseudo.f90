!> @file
!! Program to genereate the pseudopotentials used by BigDFT
!! This program is based on the sources available at
!! http://cvs.berlios.de/cgi-bin/viewcvs.cgi/cp2k/potentials/Goedecker/pseudo/v2.2/
!!
!! @author
!!     Authors:                  Raffael.Widmer (Cuda Acceleration)
!!                             and Alex Willand (all other modifications)  
!!                      under the supevision of  
!!                 Stefan Goedecker, July 2011   
!!             Universitaet  Basel, Switzerland 


!> Parallel fitting of HGH pseudopotentials using a simplex downhill method.
!! Takes into account multiple atomic references, excitation energies and
!! softness monitored by errors in 3d wavelet transformations.
!! Uses MPI and libXC, supoorts collinear spin polarization as well as 
!! nonlinear core corrections and has a GPU accelerated version of the
!! wavelet part.
      program pseudo
      use libxcModule


      implicit real*8 (a-h,o-z)
      logical fullac, avgl1,avgl2,avgl3,plotwf,denbas,mixref,pol,  &
          ortprj, litprj, energ,igrad,verbose, info,exists,ldump

      parameter ( norbmx=40, nrmax=10000, maxdim=30 )
      parameter ( lmx=5, lpmx= 4, noccmx=7, nsmx=2 )
      parameter ( ngmx=32, nintmx=5*(ngmx+14) )

      dimension aeval(noccmx,lmx,nsmx),chrg(noccmx,lmx,nsmx),  &
           dhrg(noccmx,lmx,nsmx),ehrg(noccmx,lmx,nsmx),  &
           res(noccmx,lmx,nsmx),  &
           wfnode(noccmx,lmx,nsmx,3),  &
           psi(0:ngmx,noccmx,lmx,nsmx),wght(noccmx,lmx,nsmx,8),  &
           gpot(4),hsep(6,lpmx,nsmx),r_l(lpmx),r_l2(lpmx),  &
           gcore(4),gcorepp(4),  &
           occup(noccmx,lmx,nsmx),  &
           vh((lmx+1)*((ngmx+1)*(ngmx+2))/2,  &
           (lmx+1)*((ngmx+1)*(ngmx+2))/2),xp(0:ngmx),  &
           rmt(nintmx,((ngmx+1)*(ngmx+2))/2,lmx+1),  &
           rmtg(nintmx,((ngmx+1)*(ngmx+2))/2,lmx+1),  &
           wghtcf(20),wghtex(20),  &
           ud(nintmx,((ngmx+1)*(ngmx+2))/2,lmx+1)

      dimension pp(maxdim*(maxdim+1)),yp(maxdim+1),rae(nrmax)
      dimension psiold(nrmax,noccmx,lmx+1,nsmx)
      dimension rr(nintmx),rw(nintmx),rd(nintmx)
      dimension ofdcoef(3,4), psppar(6,4,2)
      dimension excit(200)! just a temporary array for excitations

      dimension no(norbmx),noae(norbmx),lo(norbmx),so(norbmx),  &
           zo(norbmx),ev(norbmx),crcov(norbmx),dcrcov(norbmx),  &
           ddcrcov(norbmx),  &
           gf(nrmax,norbmx,nsmx),  &
           nomax(0:lmx),nomin(0:lmx),hso(6),havg(6),hhsep(6),  &
           time(3)

      character*1 il(5),ispp,isppp
      character*30 plotfile
      character(80):: errmsg
      character*7 is(2)
      character fname*35,tname*10,string*520,strcyc*10,form*25
      character*80 label
      integer namoeb,nsmplx
!      integer*4:: dateDMY(3)
      character(len=8) :: dateYMD
 
      integer ::  nconfpaw, nchannelspaw, npawl, pawstN, pawstL, pawstP
      real(8) pawrcovfact
      character(len=125) :: pawstatom


!     INCLUDE 'func.inc'
      include 'mpif.h'
!     initialize some counters
      ntime=0
      itertot=0

!     and time profiling
      time=0d0
      call cpu_time(t0)


!     this is a large number passed as reference when the penalty
!     contributions shall be written out by the penalty subroutine.

      penref=1d100


!     MPI initialization
      iproc=0
      nproc=1
      call MPI_INIT(ierr)
      call MPI_COMM_RANK(MPI_COMM_WORLD,iproc,ierr)
      call MPI_COMM_SIZE(MPI_COMM_WORLD,nproc,ierr)


!     This might not work with all compilers:
!     Generate additional output files for each
!     process iproc > 0 connected to unit 6,
!     which serves as standard output for process 0.
      if(iproc>0)then
         write(label,'(a,i2.2,a)')'proc.',iproc,'.out'
         open(6,file=trim(label))
      end if

!     DEBUG FILES: Intended to split part of the output per process       
!         write(label,'(a,i2.2,a)')'dbg.',iproc,'.out'
!         open(66,file=trim(label))
!
      write(6,*) '*********************************************'
      write(6,*) '***              pseudo_2.5               ***'
      write(6,*) '***              fitting of               ***'
      write(6,*) '***   goedecker type pseudopotentials     ***'
      write(6,*) '***   last changes:         July 2011     ***'

      if(iproc>0)then  
         write(6,'(a,i2,a,i2,a)')            &
                ' ***   output file for process ',iproc, ' of ',  &
                                                 nproc,'    ***'  
      elseif(nproc>1)then
         write(6,'(a,i2,a)')   &
                ' ***   parallel run with ', nproc,  &
                                           ' processes      ***'  
      end if
      write(6,*) '*********************************************'
      write(6,*)

!     optional info when call to the shell is ava-1ilable
      write(label,'(4a,i3.3)')'echo -n  ',  &
                              ' started  on `hostname`',  &
                              ' at `date` "... "',  &
                              ' >> hostnames.',iproc
      call system(label)


!
!     default values:
!
      NAMOEB=0
      nsmplx=2000
      ng=20
      rij=2.0
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
 

!     some defaults for new input related to the simplex
      dh=0.2d0
      ntrymax=200
      ftol=1d-7

!     some defaults related to wavelet transformations:
      nhpow=14
      nhgrid=0
      hgridmin=0.3d0
      hgridmax=0.5d0
      ampl=0d0
      crmult=8d0
      frmult=8d0

!
!     DO NOT read command line options as in previous versions.
!     This modification is needed for a parallel run.
!     Instead, read those options from a file input.dat
!     The file will be parsed for keywords line by line.
!     Some keywords are to be followed with input values.
!     For backwards compatibility, it is allowed to ommit
!     blank characters between keyword and value.



       write(6,*)'Reading the file input.dat'
       write(6,*)'__________________________'
       write(6,*)


       open(11,file='input.dat')
       do iline=1,999
            read(11,'(a)',iostat=ierr)string
            if(ierr<0)then
               write(6,'(i3,a)')iline-1,' lines have been read.'

               if(iline==1)then
                  write(6,'(a,i2)')'   Could not read the file input.dat'
                  write(6,*)
                  write(6,*) 'possible options (with input N/Float)'
                  write(6,*) '                  -cy      N ',  &
                       'number of fitting cycles (-1< N< 10000)'
                  write(6,*) '                  -fit       ',  &
                       'equivalent to -cy 1'
                  write(6,*) '                  -maxiter N ',  &
                       'number of iterations per simplex'
                  write(6,*) '                  -ng      N ',  &
                       'number of Gaussians'
                  write(6,*) '                  -rij     F ',  &
                       'max length scale of Gaussians/rloc '
                  write(6,*) '                  -fullacc   ',  &
                       'use max. number of gaussians'
                  write(6,*) '                  -orth      ',  &
                       'orthogonalisation of the projectors'
                  write(6,*) '                  -lith      ',  &
                       'transformation of the projectors as in lit.'
                  write(6,*) '                  -denbas    ',  &
                       'use dense gaussian basis    '
                  write(6,*) '                  -mixref    ',  &
                       'allow contradictions in AE ref data'
                  write(6,*) '                  -info      ',  &
                       'write gaussian coefficients of final'
                  write(6,*) '                             ',  &
                       'wavefunctions to gauss.par'
                  write(6,*) '                  -plot      ',  &
                       'plot wfs after each iteration'
                  write(6,*) '                  -lNso      ',  &
                       'avg nonlocal potential =0 for l=N  '
                  write(6,*) '                            ',  &
                       '(only for the highest projector)'
                  write(6,*) '                  -inwidth F ',  &
                       'initial width of the simplex'
                  write(6,*) '                  -stuck   N ',  &
                       'max steps to consider simplex stuck'
                  write(6,*) '                  -cvg     F ',  &
                       'convergence crit for simplex spread'
                  write(6,*)   
                  write(6,*)' Keywords related to wavelets:'
                  write(6,*) '                  -nh      N ',  &
                       'no of grids with different spacings'
                  write(6,*) '                  -hmin    F ',  &
                       'minimum value for the grid spacing '
                  write(6,*) '                  -hmax    F ',  &
                       'maximum value for the grid spacing'
                  write(6,*) '                  -nhpow   N ',  &
                       'power for weighting softness samples'
                  write(6,*) '                  -crmult  F ',  &
                       'localization radius for scaling functions'
                  write(6,*) '                  -frmult  F ',  &
                       'localization radius for wavelets'
                  write(6,*) '                  -offset  F ',  &
                       'offset between grid and origin'
                  write(6,*)
                  
                  write(6,*) '                  -pawN       ',&
                       ' no opt.,calculate  pawpatch projectors for',&
                       ' the Nth configuration '
                  write(6,*) '                  -noflpawN       ',&
                       ' pawpatch patches for the first N Ls (defaults to 3)'
                  write(6,*) '                  -nchannelspawN       ',&
                       ' set number of paw projectors to N (defaults to 4)'
                  
                  write(6,*) '                  -pawstatomName       ',&
                       '  file named Name will be read for initial ' , &
                       '  wave funct. '
                  write(6,*) '                  -pawstnN       ',&
                       '  initial wave function has first quantum number N ' 
                  write(6,*) '                  -pawstlL       ',&
                       '  initial wave function has angular momentum L' 
                  write(6,*) '                  -pawstpP       ',&
                       '  initial wave function is multiplied by r**P' 
                  write(6,*) '                  -pawrcovfactF     ',&
                       'Rbox for paw is equal to rcov*pawrcovfact. Defaults to 1' 
               end if
               exit
            elseif(ierr>0)then
               write(6,'(a,i3,a)')'Line ', iline,  &
                    ' skipped, reading error.'
               cycle
            end if
            ii=index(string,'-cy')
            if (ii.ne.0) then
               label=string(ii+3:min(ii+13,120))
               read(label,*,iostat=ierr)namoeb
               write(6,'(i4,a)') namoeb, 'fit cycles'
               if(ierr/=0)write(6,'(a,a,i3)')'Above value was set to ',   &
                    'its default due to a READING ERROR. Check line',iline
            endif
            ii=index(string,'-fit')
            if (ii.ne.0) then
               namoeb=1
               write(6,*) 'Do one fit cycle'
            endif
            ii=index(string,'-orth')
            if (ii.ne.0) then
               ortprj=.true.
               write(6,*) 'orthogonalize the projectors'
            endif
            ii=index(string,'-lith')
            if (ii.ne.0) then
               litprj=.true.
               write(6,*) 'transform the projectors as in ',  &
                    'literature'
            endif
            ii=index(string,'-maxiter')
            if (ii.ne.0) then
               label=string(ii+8:min(ii+18,120))
               read(label,*,iostat=ierr)nsmplx
               write(6,'(i4,a)') nsmplx, 'max. simplex iterations'
               if(ierr/=0)write(6,'(a,a,i3)')'Above value was set to ',   &
               'its default due to a READING ERROR. Check line',iline
            endif
            ii=index(string,'-ng')
            if (ii.ne.0) then
               label=string(ii+3:min(ii+13,120))
               read(label,*,iostat=ierr)ng
               write(6,'(i4,a)')ng,"gaussians (don't use value from atom..ae)"
               if(ierr/=0)write(6,'(a,a,i3)')'Above value was set to ',   &
               'its default due to a READING ERROR. Check line',iline
            endif
            ii=index(string,'-rij')
            if (ii.ne.0) then
               label=string(ii+4:min(ii+14,120))
               read(label,*,iostat=ierr)rij
               write(6,'(f12.4,a)')rij,' rij (overrides optional value in psppar)'
               if(ierr/=0)write(6,'(a,a,i3)')'Above value was set to ',   &
               'its default due to a READING ERROR. Check line',iline
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

           ii=index(string,'-inwidth')
            if (ii.ne.0) then
               label=string(ii+8:min(ii+18,120))
               read(label,*,iostat=ierr)dh
               write(6,'(f12.4,a)')dh,' initial width of the simplex'
               if(ierr/=0)write(6,'(a,a,i3)')'Above value was set to ',  &
               'its default due to a READING ERROR. Check line',iline
            endif
           ii=index(string,'-stuck')
            if (ii.ne.0) then
               label=string(ii+6:min(ii+16,120))
               read(label,*,iostat=ierr)ntrymax
               write(6,'(i4,a)')ntrymax,' max simplex steps when stuck'
               if(ierr/=0)write(6,'(a,a,i3)')'Above value was set to ',  &
               'its default due to a READING ERROR. Check line',iline
            endif
           ii=index(string,'-cvg')
            if (ii.ne.0) then
               label=string(ii+4:min(ii+14,120))
               read(label,*,iostat=ierr)ftol
               write(6,'(f12.4,a)')ftol,' convergence crit for the simplex'
               if(ierr/=0)write(6,'(a,a,i3)')'Above value was set to ',  &
               'its default due to a READING ERROR. Check line',iline
            endif


!          More options related to and only needed with wavelet
!          transformations, i.e. for positive weight on softness

            ii=index(string,'-nh')
            if (ii.ne.0) then
               label=string(ii+3:min(ii+13,120))
               read(label,*,iostat=ierr)nhgrid
               write(6,'(i4,a)') nhgrid, ' samples of grid spacings for WVLT'
               if(ierr/=0)write(6,'(a,a,i3)')'Above value was set to ',   &
               'its default due to a READING ERROR. Check line',iline
            endif

            ii=index(string,'-hmin')
            if (ii.ne.0) then
               label=string(ii+5:min(ii+15,120))
               read(label,*,iostat=ierr)hgridmin
               write(6,'(f12.4,a)') hgridmin, 'min grid spacing for WVLT'
               if(ierr/=0)write(6,'(a,a,i3)')'Above value was set to ',   &
               'its default due to a READING ERROR. Check line',iline
            endif

            ii=index(string,'-hmax')
            if (ii.ne.0) then
               label=string(ii+5:min(ii+15,120))
               read(label,*,iostat=ierr)hgridmax
               write(6,'(f12.4,a)') hgridmax, 'max grid spacing for WVLT'
               if(ierr/=0)write(6,'(a,a,i3)')'Above value was set to ',   &
               'its default due to a READING ERROR. Check line',iline
            endif


            ii=index(string,'-nhpow')
            if (ii.ne.0) then
               label=string(ii+6:min(ii+16,120))
               read(label,*,iostat=ierr)nhpow
               write(6,'(i4,a)') nhpow, ' power weighting for WVLT'
               if(ierr/=0)write(6,'(a,a,i3)')'Above value was set to ',   &
               'its default due to a READING ERROR. Check line',iline
            endif

            ii=index(string,'-offset')
            if (ii.ne.0) then
               label=string(ii+6:min(ii+16,120))
               read(label,*,iostat=ierr)ampl
               write(6,'(f12.4,a)') ampl, ' offset for grid for WVLT'
               if(ierr/=0)write(6,'(a,a,i3)')'Above value was set to ',   &
               'its default due to a READING ERROR. Check line',iline
            endif

            ii=index(string,'-crmult')
            if (ii.ne.0) then
               label=string(ii+7:min(ii+17,120))
               read(label,*,iostat=ierr)crmult
               write(6,'(f12.4,a)') crmult, ' coarse cutoff for WVLT'
               if(ierr/=0)write(6,'(a,a,i3)')'Above value was set to ',   &
               'its default due to a READING ERROR. Check line',iline
            endif

            ii=index(string,'-frmult')
            if (ii.ne.0) then
               label=string(ii+7:min(ii+17,120))
               read(label,*,iostat=ierr)frmult
               write(6,'(f12.4,a)') frmult, ' fine cutoff for WVLT'
               if(ierr/=0)write(6,'(a,a,i3)')'Above value was set to ',   &
               'its default due to a READING ERROR. Check line',iline
            endif


            ii=index(string,'-paw')
            if (ii.ne.0) then
               label=string(ii+4:min(ii+12,520))
               read(label,*) nconfpaw
               write(6,*)  'will calculate pawpatches for conf No ',&
                      nconfpaw
            endif
            ii=index(string,'-noflpaw')
            if (ii.ne.0) then
               label=string(ii+8:min(ii+16,520))
               read(label,*) npawl
               write(6,*)  'will calculate paw patches for the',&
               npawl, ' first Ls '

            endif

            ii=index(string,'-nchannelspaw')
            if (ii.ne.0) then
               label=string(ii+13:min(ii+21,520))
               read(label,*)  nchannelspaw
               write(6,*)  'will consider',&
              nchannelspaw , ' paw channels  '

            endif

            ii=index(string,'-pawstatom')
            if (ii.ne.0) then
               pawstatom=trim(string(ii+10:min(ii+130,520)))
               ii=index(pawstatom,' ')
               if(ii.ne.0) then
                  pawstatom=trim(pawstatom(:ii-1))
               endif
               write(6,*)  'will consider  ',&
             trim(pawstatom) ,'file for reading the initial potential '
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
               write(6,*)  ' initial wf radial part is ',&
               ' multiplied by r**' , pawstP
            endif

            ii=index(string,'-pawrcovfact')
            if (ii.ne.0) then
               label=string(ii+12:min(ii+20,520))
               print *, label
               read(label,*)  pawrcovfact
               write(6,*)  ' rbox is rcov   ',&
               ' multiplied by ' , pawrcovfact
            endif

!           loop over input lines from input.dat ends here
       end do
       close(11)


!      further output about input variables

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
!        stop
      endif


!     A little estimation on memory reqs with WVLT
         if(nhgrid>0) write(6,'(a,f4.2,a)')  &
                      'The wavelet coeffs will use about ',  &
                       8.0/2**20*int(2*crmult/hgridmin)**3,  &
                      ' MB per orbital.'
     


!
!     ----------------- read data from AE calculation -----------------
!
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
!     open(unit=40,file='atom.ae',form='formatted',status='unknown')
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

!     READING of iXC
!     for now, only keep the two most commonly used functionals
!     backwards compatible. Otherwise, require ABINITs iXC < 0

!     icorrp=label(j1:j2)
      if    (label(j1:j2)=='PADE')then
         iXC=-20
      elseif(label(j1:j2)=='PBE')then   
         iXC=-101130
      else
         read(label(j1:j2),*,iostat=ierr)iXC
         if(ierr/=0)then
          write(6,*)'Could not read the XC input in atom..ae'
          stop
         end if
      end if

      write(6,*)
!      read(40,'(t2,a)',iostat=ierr) isppp
!      read(40,'(t2,a)',iostat=ierr) icorrp
      read(40,*,iostat=ierr) ngrid
         if(ierr/=0.or.ngrid .gt. nrmax )ierr=3
         if(ngrid .gt. nrmax )  &
         write(6,*)'Input number value is to large.'
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
!     for now, use this convention
      pol=.false.
      nspol=1
      if (isppp.eq.'s') then
         pol=.true.
         nspol=2
      end if

      write(6,*)
      write(6,'(a,i7,a,i5)')  &
           ' Initializing libXC with iXC =', iXC,'; nspol =',nspol
!      the barriers are here because the init routine writes to stout
!      with each process. Need to find a way to fix this.
       if (nproc > 1) call MPI_BARRIER(MPI_COMM_WORLD,ierr)  
          call libxc_functionals_init(iXC,nspol)
       if (nproc > 1) call MPI_BARRIER(MPI_COMM_WORLD,ierr)  
         
        
      write(6,*)
      read(40,*) !some comment line
      write(6,*)' nl    s      occ        ',  &
                       'eigenvalue     charge(rcov)    '
!      write(6,*)' nl    s      occ        ',
!     :     'eigenvalue     charge(rcov)    ',
!     :     'dcharge         ddcharge'


!     this file should hold all ae plots
      if(plotwf.and.iproc==0) open(41,file='ae.orbitals.plt')
!     read the AE data
      do iorb=1,norb
         read(40,*,iostat=ierr) no(iorb),lo(iorb),so(iorb),zo(iorb),  &
              ev(iorb),crcov(iorb),dcrcov(iorb),ddcrcov(iorb)!,plotfile
!        write(6,*,iostat=ierr) no(iorb),lo(iorb),so(iorb),zo(iorb),
!    :        ev(iorb),crcov(iorb),dcrcov(iorb),ddcrcov(iorb),plotfile
         if(ierr/=0)ierr=3
         write(errmsg,*)'reading error in AE ref data'
         call errorhandler(ierr,iproc,nproc,errmsg)
       write(6,30) no(iorb),il(lo(iorb)+1),so(iorb),zo(iorb),  &
              ev(iorb),crcov(iorb)
!    :        ev(iorb),crcov(iorb),dcrcov(iorb),ddcrcov(iorb)
 30      format(1x,i1,a1,f6.1,f10.4,2f16.10,2f16.7,a)
         if(plotwf.and.iproc==0)then
!         use this statement if atom is compiled to write one file
!         per configuration and orbital
!         open(41,file=trim(plotfile))
          read(41,*)
          read(41,*)
          do igrid=1,ngrid
            read(41,*,iostat=ierr) rae(igrid),(gf(igrid,iorb,igf),  &
              igf=1,nspin)  
!           error handling in the loop is slow, but better give detailed feedback
!           for now, we only plot the ground state, i.e. atom.00.ae
            if(ierr/=0)ierr=2
!            write(errmsg,*)'error reading AE plots',
!    :                  trim(plotfile),
!     :                 'orb',iorb,'pt',igrid
!            call errorhandler(ierr,iproc,nproc,errmsg)
          end do
          read(41,*)
!         don't close this unit when reading from one single file
!         close(41)
         end if
      enddo
      if(plotwf) close(41)
      goto 457
! 456  write(6,*) 'error during reading atom.ae'
!      stop
 457  continue
      lmax=0
      lcx=0
      do iorb=1,norb
         lmax=max(lo(iorb),lmax)
         if (zo(iorb).gt.1.d-10)lcx=max(lo(iorb),lcx)
      enddo
!     print*,'lmax=',lmax
!     print*,'lcx=',lcx, '( charge > 1.0d-10)'
      if(lmax.gt.lmx+1)ierr=3
      write(errmsg,*)'array dimension problem:lmax'
      call errorhandler(ierr,iproc,nproc,errmsg)
!     compute corresponding n-quantum numbers of the pseudostates
!     no()   will contain n quantum numbers of the pseudostates afterwards
!     noae() will contain n quantum numbers of the ae-states afterwards
!     no() starts from n=1 for each(!) l
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
!      print*,'noccmax=',noccmax
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
         do j=1,ngrid
            psiold(j,nocc,l+1,ispin) = 0.0d0
!           use major comp. as reference
            if (rae(j).ne.0.0)  &
                 psiold(j,nocc,l+1,ispin)=gf(j,iorb,1)/rae(j)
         enddo
      enddo

      write(6,*)
      write(6,*) 'All electron and pseudo-wfn quantum numbers'
      write(6,*) '        n(AE)   l   inl(PS)   '
      do iorb=1,norb
         write(6,'(6x,3i6)')  noae(iorb),lo(iorb), no(iorb)
      enddo



!
!     weights will be read from weights.par
!

!     pseudo 2.4 was backwards compatible with this files
!     reading conventions from pseudo2.2 and 2.3.
!     For this version, this is not the case anymore!

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
!        what about a different order for reading these?
         read(24,*,iostat=ierr)wghtp0,wghtsoft,wghtrad,wghthij,wghtloc,wghtexci
!        no need to call error handler here, shared input file
         if(ierr/=0)  &
         write(6,*)'Reading error for weights of psi(r=0) and softness.'
         write(6,'(a,e10.3)') ' Weight for psi(r=0)=0 is     ',wghtp0
         write(6,'(a,e10.3)') ' for Ekin(gauss-wavelet)      ',wghtsoft
         write(6,'(a,e10.3)') ' for keeping radii wide       ',wghtrad
         write(6,'(a,e10.3)') ' for keeping hij small        ',wghthij
         write(6,'(a,e10.3)') ' for keeping Vloc local       ',wghtloc
         write(6,'(a,e10.3)') ' and for excitation energies  ',wghtexci
            
         read(24,*)!comment line
        
!        read the weights for eigenvalues and integrals into the array wght

         do iorb=1,norb
            nocc=no(iorb)
            l=lo(iorb)
            ispin=1
            ss = so(iorb)
            if (ss.lt.0) ispin=2
            read(24,*) nw,lw,sw,(wght(nocc,l+1,ispin,i),i=1,8)
            if (noae(iorb).ne.nw .or. l.ne.lw .or. ss.ne.sw) then
               write(6,*) 'error in file ''weights.par'' '
               write(6,*) 'need weights for n,l,s:',  &
                    noae(iorb),l,so(iorb)
               write(6,*) 'found            n,l,s:',nw,lw,sw
               if(nproc>1) call MPI_FINALIZE(ierr)
               stop
            endif
         enddo
         close(24)

!     Read excitation energies from the last line of atom.??.ae

      if(nproc>1)then
         write(6,*)
!        read etotal and exchange data with all processes
!        it would be enough to broadcast etot of system 0,
!        but this will not take much time and allow some 
!        output and proofreading.
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
         call MPI_ALLTOALL(  &
                excit(nproc+1:2*nproc),1,MPI_DOUBLE_PRECISION,   &
                excit(1:nproc),1,MPI_DOUBLE_PRECISION,  &
                MPI_COMM_WORLD,ierr)
         
         if(ierr/=0)write(6,*)'           WARNING: MPI error'        
!        total energies to excitation energies
!        different conventions of etot in gatom and atom.f
         excit=0.5d0*(excit-excit(1))
         !write(6,*)'DEBUG: No factor two in EXCITATION energies'
         !excit=(excit-excit(1))
         excitAE=excit(iproc+1)
         if(ierr/=0)excitAE=0d0
         write(6,*)
         write(6,*)'excitation energies (AE)'
         do i=0,nproc-1
             write(6,'(1x,i3,3(1x,e20.10))')  &
             i, excit(i+1)
         end do
         ierr=0
         if(excitAE==0d0.and.iproc/=0)ierr=1
         write(errmsg,*)'Excitation energy is zero and thus ignored'
         call errorhandler(ierr,iproc,nproc,errmsg)
      else
         excitAE = 0d0
      end if
!     done reading atomic ref data
      close(40)




!     ---------------------------------------------------------------
!     main loop begins here
!     ---------------------------------------------------------------

      do iiter=1,max(namoeb,1)

!        If the excitation energy has a nonzero weight
!        we always need to compute the total energy. 
!        The ground state energy is needed as a reference
!        for the excitation  energies, thus iproc==zero
         energ=(wghtexci>0d0.or.(iproc==0.and.nproc/=1))



!     read initial pseudopotential parameter from psppar
!     test if parameter from atom.ae and psppar are consistent
!     notice that fitting to inconsistent AE data MIGHT be desired
!     in a parellel run, for instance several values for rprb.



         if (iiter.eq.1) then



!           NLCC coeffs are read from both, psppar and a file 'nlcc' as used by BigDFT.
!           Both inputs are optional and one is sufficient. The input in psppar is 
!           supported in HGH-K format only, after the value nsep and on the same line. 
!           The conventions for the polynomials in BigDFT and pseudo
!           are slightly different. Values are checked for consistency.

!           initial values - negative radius means no NLCC is used 

            rcore=-1d0
            gcore(1:4)=0d0
            rcorepp=-1d0
            gcorepp(1:4)=0d0

!           read optional input file as used by BigDFT
            open(12,file='nlcc')
            read(12,*,iostat=ierr)
            read(12,*,iostat=ierr)
            read(12,*,iostat=ierrnlcc)rcore,gcore
            close(12)
            if(ierrnlcc==0)then
            write(6,*)
            write(6,*)'Reading core charge coefficients from nlcc'
            write(6,*)'__________________________________________'
            write(6,*)
!           pseudos convention is to scale the polynomial by rcore for fitting
!           convert coefficients to this convention here
            write(6,'(a,2f16.8)')' rcore and c0      ',rcore, gcore(1)
            write(6,'(a,3f16.8)')' c2 to c4          ', gcore(2:4)
            do i=2,4
               gcore(i)=gcore(i)*rcore**(2*i-2) 
            end do
            write(6,'(a,3f16.8)')' scaled by rcore:  ', gcore(2:4)
            write(6,*)
            end if

            write(6,*)
            write(6,*) 'Reading data from psppar'
            write(6,*) '________________________'
            write(6,*) 

           open(99,file='vertex.dump')

!     the format has been changed from previous version
!     to allow direct compatibility with ABINIT pseudos!

!     Output will always be in HGH-k format (pspcod=10)
!     Input formats are supported as GTH (2) and HGH (3) and HGH-k (10)
!     Notice some additional data needs to be read from this file.

            inquire(file='psppar',exist=exists) 
            if(.not.exists)then
!           no need to call errorhandler, shared file
                write(6,*)'The file psppar is lacking.'
                write(6,*)'Pseudo potentials are available from'
                write(6,*)'http://www.abinit.org/downloads/psp-links'
                if(nproc>1) call MPI_FINALIZE(ierr)
                stop
            end if


            open(unit=11,file='psppar',form='formatted',  &
                 status='unknown')

!           The first line is usually for comments only.
!           Pseudo uses it to write, read and compare some additional data


!           

!           Read 1st line into a string
            read(11,'(a)',iostat=ierr) label
!           then get optional data
            read(label,*,iostat=ierr)ispp , ng, rij, rcov, rprb
!           the very first character (when reading in list format)
!           defines the method:
!             n   non relativistc
!             r   relativistc
!             s   spin polarized (and relativistic)

!           ng and rij are the number and relative max width of the gaussians
            if(ierr/=0)then
                write(6,*)
                write(6,*)'                NOTICE'
                write(6,*)
                write(6,*)'The first line of psppar does not specify'
                write(6,*)'optional information about the calculation'
                write(6,*)'type (polarized, non-/relativistic),'
                write(6,*)'the confinement, integration radius and'
                write(6,*)'the gaussian basis. Values are taken from'
                write(6,*)'atom ae files without proofreading.'
                write(6,*)
!               This is not a reason to stop.
!               if(nproc>1) call MPI_FINALIZE(ierr)
!               stop
!               simply take the value from atom.ae, don't worry.
!               This could cause serious trouble with wavelets
!               when process 0 has a different no of spin channels
!               than some of the other processes, though.
                ispp=isppp
                rprb=rprbp
                rcov=rcovp
!               ng and rij are taken from input.dat 
!               if not specified, 20 and 2.0 is the default

                ierr=0
            else
!                read the calculation type from the label
                 j=1
                 do i=len(label),1,-1
                    if (label(i:i).ne.' ') j=i
                 enddo
                 ispp=label(j:j)
                 ierr=0
                 if(ispp.ne.'r'.and.ispp.ne.'n'.and.ispp.ne.'s')then
                    write(6,*)'The first non-blank character of psppar'
                    write(6,*)'must be one of' 
                    write(6,*)'n: for nonrelativistic calculations'
                    write(6,*)'r: for relativisitc    calculations'
                    write(6,*)'s: for spin polarized  calculations'
                    write(6,*)
                    write(6,*)'Character found:',ispp
                    write(6,*)'Exiting.'
                    call MPI_FINALIZE(ierr)
                    stop
                 end if
            end if
            if (isppp.ne.ispp) ierr=3
            write(errmsg,*)'Inconsistent spin treatment.'
            call errorhandler(ierr,iproc,nproc,errmsg)
!              below option does not really make sense.
!              it could actually be useful, but needs testing for
!              allocation errors and other causes of trouble.

!              if (mixref) then
!                 write(6,*) 'Warning! continue program using ispp',
!    :                 ' from psp.par'
!              else
!                 write(6,*) 'option ''-mixref'' allows such settings'
!                 stop
!              endif
 
            read(11,*,iostat=ierr) znuc, zion
            if(ierr/=0)then
!               no need to call error handler, shared input file
!               thus some stop statements have not been eliminated here
                write(6,*)
                write(6,*)'             WARNING'
                write(6,*)'Could not read nuclear and valence charges'
                write(6,*)'on the second line of psppar.'
                if(nproc>1) call MPI_FINALIZE(ierr)
                stop
            end if


 
            read(11,*,iostat=ierr) ipspcod, ixcpp
            if(ierr/=0)then
!               no need to call error handler, shared input file
                write(6,*)
                write(6,*)'             WARNING'
                write(6,*)'Could not read PSP format and iXC from'
                write(6,*)'the third line of psppar.'
                if(nproc>1) call MPI_FINALIZE(ierr)
                stop
            end if

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


!           already read from psppar
!           read(23,*) rcov,rprb
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
            
!           THIS DEFINITELY NEEDS MORE TESTING. GPU FAILS FOR SOME CHOICES OF RMULT
!           i.e. crmult = 2 frmult looks savest so far
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
!           no error handler needed, same for all processes.




            ierr=0
            if (rprb.ne.rprbp) then
               write(6,*)'rprb in atomic reference differs',   &
                         ' from the value in psppar.'
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
!           read(23,'(a)') label
!           j1=1
!           j2=2
!           do i=len(label),1,-1
!              if (label(i:i).ne.' ') j1=i
!           enddo
!           do i=len(label),j1,-1
!              if (label(i:i).eq.' ') j2=i
!           enddo
!           j2=j2-1
!           icorr=label(j1:j2)
!           if (icorr.ne.icorrp) then
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

!     here we previously set the parameter/variables for the xc-functional(s)
!     

!     NOTE: znuc MUST be consistent, mixed references make no sense here
!           zion WILL be different as soon as we have ionic configurations. 

!           read(11,*) znuc, zion, rloc, gpot(1),gpot(2),gpot(3),gpot(4)
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
!              zion=zionp
            endif
            write(errmsg,*) 'valence charge differs from AE data'
            call errorhandler(ierr,iproc,nproc,errmsg)

!           read the pseudopotential the way BigDFT does

!           be sure to have zeroes for undefined entries.
            psppar = 0d0

            if(ipspcod==10.or.ipspcod==11)then
                 write(6,*)'HGH matrix format'
!              ! HGH-K format: all projector elements given.
               ! dimensions explicitly defined for nonzero output.

!              ! local part
               read(11,*) rloc,nloc,(gpot(j),j=1,nloc) 

!              read line with optional nlcc data via some string
!              to prevent jumping to the next line when the input is absent
               read(11,'(a)')string
               read(string,*,iostat=ierrpp) lpx,rcorepp,gcorepp 


!              report errors or contradictions
               if(ierrpp/=0.and.ierrnlcc==0)then
                 write(6,*)'                NOTICE'
                 write(6,*)  &
              'nlcc file present and used, but no NLCC in psppar found.'
                 write(6,*)
                 rcorepp=rcore
                 gcorepp=gcore
               end if
               if(ierrpp==0.and.ierrnlcc/=0 &
                           .and.sum(gcorepp**2)>1d-12)then
                 write(6,*)'                NOTICE'
                 write(6,*)  &
              'Using NLCC from psppar, but no nlcc file could be read.'
                 write(6,*)
                 rcore=rcorepp
                 gcore=gcorepp
               end if
               tt=maxval(abs(gcore-gcorepp))
               tt=max(tt,abs(rcore-rcorepp))
               if(tt>1d-4)then
                 write(6,*)'                WARNING'
                 write(6,*)  &
              'NLCC data in "nlcc" (used) and "psppar" (ignored) differ'
                 write(6,'(a,5f12.6)')'   nlcc',rcore,gcore
                 write(6,'(a,5f12.6)')' psppar',rcorepp,gcorepp
                 write(6,*)
               end if

!              lpx is here equivalent to nsep. Previous versions used
!              it as the max l quantum number, subracting one.
!              Here, 0 does not mean s only, but a purely local psppar.
!              Negative is no longer used for local, but for r_l2.
               if (lpx-1 .gt. lmx ) then
                 write(6,*) 'array dimension problem: lpx,lpmx',lpx,lpmx
                 if (nproc > 1) call MPI_FINALIZE(ierr)
                 stop
               end if
!              ! separable part
!              ! relativistic; hij are averages, kij separatons of hsep
               do l=1,lpx !l channels
                  ! add some line to read r_l2 if nprl < 0
                  read(11,'(a)') string
                  read(string,*) r_l(l),nprl
                  if(nprl>0)then
                     read(string,*) r_l(l),nprl,  &
                         psppar(1,l,1),(psppar(j+2,l,1),  &
                                               j=2,nprl)  !h_ij 1st line
                     do i=2,nprl
!                       spin up
                        read(11,*) psppar(i,l,1),(psppar(i+j+1,l,1),  &
                                                       j=i+1,nprl)!remaining h_ij 
                     end do
!                    disable r_l2, i.e set it equal to r_l
                     r_l2(l) = r_l(l)
                  else
!                    if nprl is negative, read r_l2 from the 2nd line of hij
                     nprl=-nprl
                     read(string,*) r_l(l),i,  &
                         psppar(1,l,1),(psppar(j+2,l,1),  &
                                               j=2,nprl)  !h_ij 1st line
                     read(11,*)r_l2(l), psppar(2,l,1),(psppar(2+j+1,l,1),  &
                                                       j=2+1,nprl)!2nd line
                     if(nprl==3) read(11,*) psppar(3,l,1)! thid line
                  end if

!                 there are no kij the s-projector
                  if (l==1) cycle
                   do i=1,nprl
                      read(11,*) psppar(i,l,2),(psppar(i+j+1,l,2),  &
                                             j=i+1,nprl)!all k_ij
                   end do
               end do ! loop over l 
            elseif(ipspcod==3)then
                 write(6,*)'HGH diagonal format'
!                HGH diagonal part case
!                technically, lpx is fixed at the max value of
                 lpx=4
                 read(11,*) rloc,(gpot(j),j=1,4)
                 read(11,*) r_l(1),psppar(1:3,1,1) 
                 do l=2,4
                    read(11,*) r_l(l),psppar(1:3,l,1) 
                    read(11,*)        psppar(1:3,l,2) 
                 end do
            elseif(ipspcod==2)then
                 write(6,*)'GTH format'
!                ! GTH case
!                technically, lpx is fixed at s and p
                 lpx=2
                 read(11,*) rloc,(gpot(j),j=1,4)
                 read(11,*) r_l(1),psppar(1:2,1,1)
                 read(11,*) r_l(2),psppar(1  ,2,1)
            else
!                no need to call error handler, shared input file
                 write(6,*)'               WARNING'
                 write(6,*)'pspcod (1st number of 3rd line) read from' 
                 write(6,*)'psppar is unknown or not supported.' 
                 write(6,*)'supported are 2,3, or 10, not ',ipspcod
                 if (nproc > 1) call MPI_FINALIZE(ierr)
                 stop
            end if

!           done with reading psppar
            close(11) 

!           avoid radii equal zero, even for unused projectors. 
            do l=1,lpx
              if(r_l(l)==0d0)then
                write(6,*)'all r_l should be nonzero.'
                write(6,*)'The r_l of the ',il(l),'-projector has been'
                write(6,*)'adjusted from 0 to 1. Check your psppar.'
                r_l(l)=1d0
                r_l2(l)=1d0
              end if
              if(r_l2(l)==0d0) r_l2(l) = r_l(l)
            end do
            if( sum((r_l2-r_l)**2) > 1d-6) then
                  write(6,*)
                  write(6,*)'               NOTE'
                  write(6,*)'The separable part will use two length'
                  write(6,*)'scales as read from psppar.'
                  write(6,*)'  l             r_l            r_l2'
                  do l=1,lpx
                      write(6,'(i4,5x,2f15.5)')l-1,r_l(l), r_l2(l)
                  end do
                  write(6,*)
            end if
 


!           and to be very sure
!           in case lpx increses later
            r_l(lpx+1:lpmx)=1d0         
            r_l2(lpx+1:lpmx)=1d0         


!           Then rearrange pspcod into hsep accordingly and 
!           convert hij,kij to hij(up,down)  

!           pspcod  as read are packed upper diagonal ROW elements of

!               h =  ((l+1) hup + l hdn)/(2l+1)
!               k =      2( hup -   hdn)/(2l+1)

!           we want hsep,  packed upper diagonal COL elements of

!             hup = h +   l/2   k
!             hdn = h - (l+1)/2 k


!           1) get offdiagonal elements where needed
!           2) permute from row to col ordering
!           3) transform from h/k to up/down


!           just to be sure no rubbish will be added 
            hsep=0d0
            
            if (ipspcod == 2) then !GTH case
!              offdiagonal elements are zero per definition.
!              simply rearrange hij and fill zero elements
               do l=1,lpx
                  hsep(1,l,1)=psppar(1,l,1)
                  hsep(2,l,1)=0.0d0
                  hsep(3,l,1)=psppar(2,l,1)
                  hsep(4,l,1)=0.0d0
                  hsep(5,l,1)=0.0d0
                  hsep(6,l,1)=psppar(3,l,1)
!                 in the polarized or relativistic case,
!                 we assume all kij to be zero,
!                 i.e. same up and down projectors
                  if(nspin==2) hsep(:,l,2)=hsep(:,l,1)
               end do
            elseif (ipspcod == 3) then !HGH diagonal case
!              we need to compute the offdiagonal elements with the following coeffs
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
            
!              this could possibly be done in a simpler way ...

               do l=1,lpx; do ispin=1,nspin
                  hsep(1,l,ispin)=psppar(1,l,ispin)
                  hsep(2,l,ispin)=psppar(2,l,ispin)*ofdcoef(1,l)
                  hsep(3,l,ispin)=psppar(2,l,ispin)
                  hsep(4,l,ispin)=psppar(3,l,ispin)*ofdcoef(2,l)
                  hsep(5,l,ispin)=psppar(3,l,ispin)*ofdcoef(3,l)
                  hsep(6,l,ispin)=psppar(3,l,ispin)
               end do; end do


!              in the nonrelativistic case, we are done.
               if(nspin==2)then
!                in the polarized case, copy the missing s projector
                 if(pol)  hsep(:,1,2)=hsep(:,1,1)
!                use psppar as a temporary array 
                 psppar=hsep
!                and then convert hij/kij to hij up/down 
                 do l=2,lpx
!                 l is index, angular momentum +one
!                 up
                  hsep(1,l,1)=psppar(1,l,1) +.5d0*(l-1)* psppar(1,l,2)!h11
                  hsep(2,l,1)=psppar(2,l,1) +.5d0*(l-1)* psppar(2,l,2)!h12
                  hsep(3,l,1)=psppar(3,l,1) +.5d0*(l-1)* psppar(3,l,2)!h22
                  hsep(4,l,1)=psppar(4,l,1) +.5d0*(l-1)* psppar(4,l,2)!h13
                  hsep(5,l,1)=psppar(5,l,1) +.5d0*(l-1)* psppar(5,l,2)!h23
                  hsep(6,l,1)=psppar(6,l,1) +.5d0*(l-1)* psppar(6,l,2)!h33
!                 down
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
!              psppar holds hij and kij in HGHK convention
!              fill hsep(up,dn) upper diagonal col by col, as needed for the fit

!             for a nonrelativistic calculation, discard the kij elements
              if(nspin==1)then
               do l=1,lpx
                  hsep(1,l,1)=psppar(1,l,1) !h11
                  hsep(2,l,1)=psppar(4,l,1) !h12
                  hsep(3,l,1)=psppar(2,l,1) !h22
                  hsep(4,l,1)=psppar(5,l,1) !h13
                  hsep(5,l,1)=psppar(6,l,1) !h23
                  hsep(6,l,1)=psppar(3,l,1) !h33
               end do
              else
!              relativistic or polarized calculation
               do l=1,lpx
!                 l is the index, angular momentum +one
!                 up
                  hsep(1,l,1)=psppar(1,l,1) +.5d0*(l-1)* psppar(1,l,2)!h11
                  hsep(2,l,1)=psppar(4,l,1) +.5d0*(l-1)* psppar(4,l,2)!h12
                  hsep(3,l,1)=psppar(2,l,1) +.5d0*(l-1)* psppar(2,l,2)!h22
                  hsep(4,l,1)=psppar(5,l,1) +.5d0*(l-1)* psppar(5,l,2)!h13
                  hsep(5,l,1)=psppar(6,l,1) +.5d0*(l-1)* psppar(6,l,2)!h23
                  hsep(6,l,1)=psppar(3,l,1) +.5d0*(l-1)* psppar(3,l,2)!h33
!                 if nspol==1 and l==1
!                 if(l==2-nspol)cycle
!                 down
                  hsep(1,l,2)=psppar(1,l,1) -.5d0* l   * psppar(1,l,2)!h11
                  hsep(2,l,2)=psppar(4,l,1) -.5d0* l   * psppar(4,l,2)!h12
                  hsep(3,l,2)=psppar(2,l,1) -.5d0* l   * psppar(2,l,2)!h22
                  hsep(4,l,2)=psppar(5,l,1) -.5d0* l   * psppar(5,l,2)!h13
                  hsep(5,l,2)=psppar(6,l,1) -.5d0* l   * psppar(6,l,2)!h23
                  hsep(6,l,2)=psppar(3,l,1) -.5d0* l   * psppar(3,l,2)!h33
               end do
              end if
            end if





!           output
            write(6,*)
            write(6,'(2f10.3,3x,a)') rcov,rprb,  &
                 ' rcov and rprb (charge integr and confinement)'
            if (ispp.eq.'r') then
               write(6,'(t30,a)')'relativistic calculation'
            elseif(ispp.eq.'s')then
               write(6,'(t30,a)')'spin polarized calculation'
            else
               write(6,'(t30,a)')'non relativistic calculation'
            endif
            write(6,'(t2,i10,t30,a)')iXC ,'ixc for XC-functional'
            write(6,*) 
            write(6,*) 'local part'
            write(6,'(f5.0,f7.2,f7.3,4e11.3,t65,a)')  &
                 znuc,zion,rloc,gpot(1),gpot(2),gpot(3),gpot(4),  &
                 'znuc,zion, rloc, gpot() '
            if (lpx.ge.1) then
               write(6,*) 'nonlocal part in internal format'
               write(6,'(i4,t60,a)') lpx ,  &
                    'lpx, (Projectors for l=0..lpx)'
               do l=1,lpx
                  write(6,*) il(l)//'-projector'
                  write(6,'(f7.3,t8,6e11.3,t76,a)') r_l(l),  &
                       (hsep(i,l,1),i=1,6),'r_l(),hsep(), '//is(1)
                  if (l.gt.1 .and. nspin.eq.2)  &
                       write(6,'(t8,6e11.3,t76,a)')  &
                       (hsep(i,l,2),i=1,6),'      hsep(), '//is(2)
               enddo
            endif
            if(rcore>0)then
               write(6,*)
               write(6,*)'This calculation will use NLCC'
               write(6,*)'(nonlinear core corrections).'
               write(6,*)'Parameters rcore and gcore(1:4) are'
               write(6,'(1x,5e18.6)')rcore,gcore
               write(6,*)
            end if
         endif ! first iteration

!
!     previous versions have read weights.par at each iteration step.
!

!     let us sacrifice this feature for the sake of stability in
!     parallel runs.


!

!     calc. exponents of gaussians
         a0=rloc/rij
!     take this for an crude initial fit:
!     tt=2.d0**.4d0
         if (denbas) then
!     fine fit:
            tt=sqrt(sqrt(2.d0))
         else
!     normal fit:
            tt=2.d0**.3d0
         endif
         a=a0
         do i=0,ng
!           
!           a=a0*tt**i
            xp(i)=.5d0/a**2
            a=a*tt
         enddo
         write(6,*)
         write(6,*)'Gaussian basis'
         write(6,*)'______________'
         write(6,*)
         write(6,'(a,4e11.4)') ' amin,amax',a0,a
         write(6,'(a,t10,3(e11.4),a,2(e11.4))') ' gaussians ',  &
              xp(1),xp(2),xp(3),' .... ',xp(ng-1),xp(ng)  
         write(6,*)'gaussians:',ng
!     set up radial grid
         nint=5*(ng+14)
         rmax=min(15.d0*rprb,120.d0)
         a_grd=a0/400.d0
         b_grd=log(rmax/a_grd)/(nint-2)
         call radgrid(nint,rr,rw,rd,a_grd,b_grd,rmax)
         write(6,'(a,t10,3(e11.4),a,2(e11.4))') ' r-grid: ',  &
              rr(1),rr(2),rr(3),' .... ',rr(nint-1),rr(nint)  
         write(6,*)'gridpoints:',nint
         write(6,*)
         call crtvh(ng,lcx,lmax,xp,vh,nint,rmt,rmtg,ud,rr)
         write(6,*)


         if (namoeb.gt.0) then

!     some OLD, commmented out feature
!     refine simplex only every 10.th step
!     start of if-block
!          if (mod(iter,10).eq.0) then


!
!     pack initial guess
!

         write(6,*)'Reading fitting parameters from FITPAR'
         write(6,*)'______________________________________'
         write(6,*)

            verbose=.true.
            call  ppack (verbose,rloc,gpot,hsep,r_l,r_l2,pp(1),  &
                 lpx,lpmx,nspin,pol,nsmx,maxdim,nfit,'init',  &
                 avgl1,avgl2,avgl3,ortprj,litprj,  &
                 rcore,gcore,znuc,zion)
            verbose=.false.
!
!     initial simplex
!

!           copy the packed initial vertex to the other vertices
!           does not make much sense, though, as we copy zeroes only.
!           keep it in case ppack.f is changed
            do i=1,nfit
               do j=1,nfit
                  pp(j+i*nfit)=pp(j)
               enddo
            enddo

!           This width is not hard coded anymore
!           dh=0.2d0

!           shift vertices number 2 to nfit+1 by random numbers in [-dh/2:dh/2]
!           in this convention, only move the k-th component of vertex k+1

!           the random numbers generated SHOULD be the same for all processes.
!           Though, let us enforce equality with MPI broadcast from process 0.
            if(iproc==0)then
              do i=1,nfit
!     f90 intrinsic
                call random_number(randNr)
                pp(i+i*nfit)=pp(i)+dh*1.0d0*(randNr-.5d0)
!     Intel (ifc)
!              pp(i+i*nfit)=pp(i)+dh*1.0d0*(dble(rand(0.0d0))-.5d0)
!     IBM/DEC/PGI
!              pp(i+i*nfit)=pp(i)+dh*1.0d0*(dble(rand())-.5d0)
!     CRAY
!             pp(i+i*nfit)=pp(i)+dh*1.0d0*(ranf()-.5d0)
              enddo
            end if
            call MPI_BCAST(pp,nfit*(nfit+1),MPI_DOUBLE_PRECISION,0,  &
                                            MPI_COMM_WORLD,ierr)

            write(6,*)
            write(6,'(1x,a,i4)')'starting amoeba cycle',iiter
            write(6,*)          '_________________________'
            write(6,*)
             write(6,'(2a)')' Penalty contributions of the initial ',  &
                            'parameters and random simplex vertices'
             write(6,'(2a)')' _______________________________________',  &
                            '____________________________________'
             write(6,*)
             if(nproc>1)then
                write(6,'(a)')'  amoeba   vertex    overall penalty'//  &
!    :               '      softness             narrow radii      '//  &
                     '      softness             psp empirical     '//  &
                     '   all configurations   all excitation       '//  &
                     'this configuration'
                write(6,'(a)')'    step   number    value          '//  &
                     '      gauss-wvlt           (hgridmax/r)**12  '//  &
                     '   E_KS, integrals      energies             '//  &
                     'and its excitation'
             else
                write(6,'(a)')'    step    vertex   penalty,       '//  &
!    :               '      gaus-wvlt            narrow radii      '//  &
                     '      gaus-wvlt            psp empirical     '//  &
                     '   eigenvalues etc'
             end if
!           do i=2,nfit+1  would be the random vertices, 1 has no shift
            iter=0
            do i=1,nfit+1
               call penalty(energ,verbose,nfit,pp(1+(i-1)*nfit),yp(i),  &
                 noccmax,noccmx,lmax,lmx,lpx,lpmx,lcx,nspin,pol,nsmx,  &
                 no,lo,so,ev,crcov,dcrcov,ddcrcov,norb,  &
                 occup,aeval,chrg,dhrg,ehrg,res,wght,  &
                 wfnode,psir0,wghtp0,  &
                 rcov,rprb,rcore,gcore,znuc,zion,rloc,gpot,r_l,r_l2,hsep,  &
                 vh,xp,rmt,rmtg,ud,nint,ng,ngmx,psi,  &
                 avgl1,avgl2,avgl3,ortprj,litprj,igrad,rr,rw,rd,  &
!                the following lines differ from pseudo2.2  &
                 iproc,nproc,wghtconf,wghtexci,wghtsoft,wghtrad,wghthij,wghtloc,  &
                 nhgrid, hgridmin,hgridmax, nhpow,ampl,crmult,frmult,  &
                 excitAE,ntime,iter,itertot,penref,time)
!              write(6,'(4x,i3,4x,e20.10)')i,yp(i)
!MK            CALL FLUSH(6)
            enddo
            write(99,*)'history of lowest vertices'


             write(6,*)
             write(6,'(2a)')' Penalty contributions of the',  &
                            ' currently lowest vertex'
             write(6,'(2a)')' _______________________________________',  &
                            '_____________'
             write(6,*)
             if(nproc>1)then
                write(6,'(a)')'  amoeba   gatom     overall penalty'//  &
!    :               '      softness             narrow radii      '//  &
                     '      softness             psp empirical     '//  &
                     '   all configurations   all excitation       '//  &
                     'this configuration'
                write(6,'(a)')'    step   calls     value          '//  &
                     '      gauss-wvlt           (hgridmax/r)**12  '//  &
                     '   E_KS, integrals      energies             '//  &
                     'and its excitation'
             else
                write(6,'(a)')'    step    gatom    penalty,       '//  &
!    :               '      gaus-wvlt            narrow radii      '//  &
                     '      gaus-wvlt            psp empirical     '//  &
                     '   eigenvalues etc'
             end if



!            The vertex number 1 holds the initial psp

!            this call is not be needed anymore, because
!            we already wrote details of all vertices

!            iter=0
!            write(*,*)'test line, vertex 1 again'
!            call penalty(energ,verbose,nfit,pp(1),yp(1),
!    :          noccmax,noccmx,lmax,lmx,lpx,lpmx,lcx,nspin,pol,nsmx,
!    :          no,lo,so,ev,crcov,dcrcov,ddcrcov,norb,
!    :          occup,aeval,chrg,dhrg,ehrg,res,wght,
!    :          wfnode,psir0,wghtp0,
!    :          rcov,rprb,rcore,gcore,znuc,zion,rloc,gpot,r_l,r_l2,hsep,
!    :          vh,xp,rmt,rmtg,ud,nint,ng,ngmx,psi,
!    :          avgl1,avgl2,avgl3,ortprj,litprj,igrad,rr,rw,rd,
!               the following lines differ from pseudo2.2
!    :          iproc,nproc,wghtconf,wghtexci,wghtsoft,wghtrad,wghthij,wghtloc,
!    :          nhgrid, hgridmin,hgridmax, nhpow,ampl,crmult,frmult,
!    :          excitAE,ntime,iter,itertot,penref,time)


!     refine simplex only every 10.th step
!     end of if-block
!         endif
!
!     starting amoeba
!
!           ftol=1.d-7 is not hardcoded anymore
            call AMOEBA(pp,yp,nfit,FTOL,ITER,nsmplx,namoeb,  &
                 noccmax,noccmx,lmax,lmx,lpx,lpmx,lcx,nspin,pol,nsmx,  &
                 no,lo,so,ev,crcov,dcrcov,ddcrcov,norb,  &
                 occup,aeval,chrg,dhrg,ehrg,res,wght,  &
                 wfnode,psir0,wghtp0,  &
                 rcov,rprb,rcore,gcore,znuc,zion,rloc,gpot,r_l,r_l2,hsep,  &
                 vh,xp,rmt,rmtg,ud,nint,ng,ngmx,psi,  &
                 avgl1,avgl2,avgl3,ortprj,litprj,igrad,rr,rw,rd,  &
!                the following line differs from pseudo2.2  &
                 iproc,nproc,wghtconf,wghtexci,wghtsoft,wghtrad,wghthij,wghtloc,  &
                 nhgrid,hgridmin,hgridmax,nhpow,ampl,crmult,frmult,  &
                 ntrymax,excitAE,ntime,itertot,energ,verbose,time)
!           write(6,*) 'Finished amoeba with ',iter,'iterations'
         else

!            No fitting was requested, evaluate penalty contributions once

         
!        call penalty with the verbose flag to print the details
            verbose=.false.
            energ=.true.
            penref=-1d9
  
            call ppack (verbose,rloc,gpot,hsep,r_l,r_l2,pp(1),  &
                 lpx,lpmx,nspin,pol,nsmx,maxdim,nfit,'init',  &
                 avgl1,avgl2,avgl3,ortprj,litprj,  &
                 rcore,gcore,znuc,zion)
            verbose=.true.

            if (nconfpaw/=-1) then
               call pawpatch(energ,verbose,nfit,pp(1),yp(1), &
                             noccmax,noccmx,lmax,lmx,lpx,lpmx,lcx,nspin,pol,nsmx,&
                             no,lo,so,ev,crcov,dcrcov,ddcrcov,norb,&
                             occup,aeval,chrg,dhrg,ehrg,res,wght,&
                             wfnode,psir0,wghtp0,&
                             rcov,rprb,rcore,zcore,znuc,zion,rloc,gpot,r_l,hsep,&
                             vh,xp,rmt,rmtg,ud,nint,ng,ngmx,psi,&
                             avgl1,avgl2,avgl3,ortprj,litprj,igrad,rae,&
!                               the following lines differ from pseudo2.2
                             iproc,nproc,wghtconf,wghtexci,wghtsoft,wghtrad,&
                         wghthij,&
                             nhgrid,hgridmin,hgridmax,nhpow,ampl,crmult,frmult,&
                             excitAE,ntime,iter,itertot,penref,time,ngrid, &
                              nconfpaw, npawl, nchannelspaw , ispp, pawstatom,&
                             pawstN, pawstL  , pawstP,     pawrcovfact    )
               
               stop " pawpatch normal stop"
            endif

            call penalty(energ,verbose,nfit,pp(1),yp(1),  &
                 noccmax,noccmx,lmax,lmx,lpx,lpmx,lcx,nspin,pol,nsmx,  &
                 no,lo,so,ev,crcov,dcrcov,ddcrcov,norb,  &
                 occup,aeval,chrg,dhrg,ehrg,res,wght,  &
                 wfnode,psir0,wghtp0,  &
                 rcov,rprb,rcore,gcore,znuc,zion,rloc,gpot,r_l,r_l2,hsep,  &
                 vh,xp,rmt,rmtg,ud,nint,ng,ngmx,psi,  &
                 avgl1,avgl2,avgl3,ortprj,litprj,igrad,rr,rw,rd,  &
!                the following lines differ from pseudo2.2  &
                 iproc,nproc,wghtconf,wghtexci,wghtsoft,wghtrad,wghthij,wghtloc,  &
                 nhgrid,hgridmin,hgridmax,nhpow,ampl,crmult,frmult,  &
                 excitAE,ntime,iter,itertot,penref,time)
                 write(6,*)'Overall penalty value would be',yp(1)


!          IMPORTANT TEST: Calling gatom once more should not make any difference 
!     call gatom(nspol,energ,verbose,
!    :     noccmax,noccmx,lmax,lmx,lpx,lpmx,lcx,nspin,nsmx,
!    :     occup,aeval,chrg,dhrg,ehrg,res,wght,wfnode,psir0,
!    :     rcov,rprb,rcore,gcore,znuc,zion,rloc,gpot,r_l,r_l2,hsep,
!    :     vh,xp,rmt,rmtg,ud,nint,ng,ngmx,psi,igrad,
!    :     rr,rw,rd,ntime,itertot,etotal)

!          ...true
          endif
!         fit or evaluate once


!
!     print results
!
         write(6,*)
         write(6,*)'Penalty contribtions from this configuration'
         write(6,*)'____________________________________________'
         write(6,*)
         write(6,'(2(tr10,a,e12.4))')  &
              'psi(r=0) =',psir0,'; psi(0)*wght=',abs(psir0*wghtp0)
         write(6,*)
         write(6,'(a,t32,a,t42,a,t55,a,t64,a)')  &
              ' nl    s      occ','ae','pseudo','diff','diff*weight'

         do iorb=1,norb
            write(6,31) noae(iorb),il(lo(iorb)+1),so(iorb),zo(iorb)
 31         format(1x,i1,a1,f6.1,f10.4)
            nocc=no(iorb)
            l=lo(iorb)
            ispin=1
            if (so(iorb).lt.0) ispin=2
            write(6,32) 'eigenvalue ',  &
                 ev(iorb),aeval(nocc,l+1,ispin),  &
                 aeval(nocc,l+1,ispin)-ev(iorb),  &
                 abs(wght(nocc,l+1,ispin,1)*  &
                 (aeval(nocc,l+1,ispin)-ev(iorb)))
            write(6,32) 'charge    ',  &
                 crcov(iorb),chrg(nocc,l+1,ispin),  &
                 chrg(nocc,l+1,ispin)-crcov(iorb),  &
                 abs(wght(nocc,l+1,ispin,2)*  &
                 (chrg(nocc,l+1,ispin)-crcov(iorb)))
            if (wght(nocc,l+1,ispin,3).ne.0.0d0)  &
                 write(6,32) 'dcharge   ',  &
                 dcrcov(iorb),dhrg(nocc,l+1,ispin),  &
                 100.d0*abs(1.d0-dhrg(nocc,l+1,ispin)/dcrcov(iorb)),  &
                 abs(wght(nocc,l+1,ispin,3))*  &
                 100.d0*abs(1.d0-dhrg(nocc,l+1,ispin)/dcrcov(iorb))
            if (wght(nocc,l+1,ispin,4).ne.0.0d0)  &
                 write(6,32) 'echarge   ',  &
                 ddcrcov(iorb),ehrg(nocc,l+1,ispin),  &
                 100.d0*abs(1.d0-ehrg(nocc,l+1,ispin)/ddcrcov(iorb)),  &
                 abs(wght(nocc,l+1,ispin,4))*  &
                 100.d0*abs(1.d0-ehrg(nocc,l+1,ispin)/ddcrcov(iorb))
            write(6,33) 'residue   ',  &
                 res(nocc,l+1,ispin),  &
                 abs(wght(nocc,l+1,ispin,5)*res(nocc,l+1,ispin))
            if (wght(nocc,l+1,ispin,6).ne.0.0d0)  &
                 write(6,33) 'rnode     ',  &
                 wfnode(nocc,l+1,ispin,1),  &
                 abs(wght(nocc,l+1,ispin,6)*wfnode(nocc,l+1,ispin,1))
            if (wght(nocc,l+1,ispin,7).ne.0.0d0)  &
                 write(6,33) 'dnode     ',  &
                 wfnode(nocc,l+1,ispin,2),  &
                 abs(wght(nocc,l+1,ispin,7)*wfnode(nocc,l+1,ispin,2))
            if (wght(nocc,l+1,ispin,8).ne.0.0d0)  &
                 write(6,33) 'ddnode    ',  &
                 wfnode(nocc,l+1,ispin,3),  &
                 abs(wght(nocc,l+1,ispin,8)*wfnode(nocc,l+1,ispin,3))
            write(6,*)
         enddo
 32      format (t10,a,t25,4e12.4)
 33      format (t10,a,t25,2e24.4)
         write(6,*) 'diff for dcharg and echarge is given in (%)'

!        always print out the psppar, also if no fitting was done.
!        This is useful to convert to pspcod 10 without fitting,
!        and to include some extra information.
         write(6,*)
         write(6,*) 'Resulting Pseudpotential Parameters'
         write(6,*) '___________________________________'
         write(6,*) 

!     At this point, overwrite  the psppar in ABINIT pspcod=10 format
!     Do the output twice, to psppar and to the logfile(s)
!     Additional output beyond the standad format:
!        some info about method, basis set, rcov and rprb on line 1
!        if NLCC is used, coeffs are written after nsep on the same line
!        If r_l2 are used, put them in front of the h12
!             and give lpj a negative sign, usually -2.


         if(iproc==0)then
            open(unit=13,file='psppar')!,position='append')

            if(ispp=='r')then
              write(13,'(a)',advance='no')  &
              'relativistic '
            elseif(ispp=='s')then
              write(13,'(a)',advance='no')  &
              'spin-polarized '
            else
              write(13,'(a)',advance='no')  &
              'nonrelativistic '
            end if
            write(13,'(i3,3(1x,g9.3),a)')ng,rij,rcov,rprb,  &
                                ' spin treatment, ng rij rcov and rprb'

!           this intrinsic may have different conventions on its arguments
!            call idate(dateDMY)
!            call idate(dateDMY(1),dateDMY(2),dateDMY(3))
!            write(13,'(2i5,2x,3i2.2,a)')  &
!                              int(znuc+.1),int(zion+.1),dateDMY,  &
!                              ' zatom, zion, date (ddmmyy)'
!            write( 6,'(2i5,2x,3i2.2,a)')  &
!                              int(znuc+.1),int(zion+.1),dateDMY,  &
!                             ' zatom, zion, date (ddmmyy)'
            call date_and_time(dateYMD)
            write(13,'(2i5,2x,2a)') int(znuc+.1),int(zion+.1),dateYMD,' zatom, zion, date (yymmdd)'
             write( 6,'(2i5,2x,2a)') int(znuc+.1),int(zion+.1),dateYMD,' zatom, zion, date (yymmdd)'
 
            write(13,'(2x,i3,i8,i2,a)')10,ixc,lpx-1,  &
                   ' 0 2002 0     pspcod, IXC, lmax, lloc, mmax, r2well'
            write( 6,'(2x,i4,i7,i2,a)')10,ixc,lpx-1,  &
                   ' 0 2002 0     pspcod, IXC, lmax, lloc, mmax, r2well'

!           determine the number of nonzero terms in the local potential
            ngpot=0
            do j=1,4
               if (gpot(j).ne.0.d0) ngpot=j
            enddo
            if(ngpot==0)then
              write(13,'(2x,f16.8,a)')rloc,' 0 rloc nloc ck (none)'
              write( 6,'(2x,f16.8,a)')rloc,' 0 rloc nloc ck (none)'
            else
              write(label,'(a,i1.1,a)')'(2x,f16.8,i3,',ngpot,'f16.8,a)'
              write(13,label)rloc,ngpot,gpot(1:ngpot),  &
                     ' rloc nloc c1 .. cnloc'
              write( 6,label)rloc,ngpot,gpot(1:ngpot),  &
                     ' rloc nloc c1 .. cnloc'
            end if

!           write out NLCC coeffs only if used.
            if(rcore>0d0)then
               write(13,'(i4,1x,5f16.8,a)')lpx, rcore,gcore,  &
                                     ' nsep, rcore, gcore (NLCC)'
               write( 6,'(i4,1x,5f16.8,a)')lpx, rcore,gcore,  &
                                     ' nsep, rcore, gcore (NLCC)'
            else
               write(13,'(i4,55x,a)')lpx,'nsep' 
               write( 6,'(i4,55x,a)')lpx,'nsep' 
            end if

!           Write out the separable part
            do l=1,lpx
!              l is a positive index, lq the l quantum number
               lq=l-1
               do i=1,6
!                 In the relativistic case and for l > s,
!                 we need to compute havg and hso from hup and hdown
!                 needed if: l > 1, nspin =2 and nspol = 1
                  if ( l>1 .and. nspin > nspol ) then
                     
                     havg(i)=((lq+1)*hsep(i,l,1)+lq*hsep(i,l,2))  &
                          /(2*lq+1)
                     hso(i)=2*(hsep(i,l,1)-hsep(i,l,2))  &
                          /(2*lq+1)
                  else
                     havg(i)=hsep(i,l,1)
                     hso(i)=0d0
                  endif
               enddo
!              get the matrix dimension lpj = 1, 2 or 3
               lpj=1
               if(max(abs(havg(3)),abs(hso(3)))>1d-8)lpj=2
               if(max(abs(havg(6)),abs(hso(6)))>1d-8)lpj=3

!              Then see if the r_l2 feature was used, compare with r_l
               if((r_l2(l)-r_l(l))**2>1d-8) then 
!                  this will be written in front of h12 
                   write(label,'(f16.8)') r_l2(l)
!                  negative sign for the dimension integer 
                   write(tname,'(i3)') -lpj
               else
!                  No r_l2 used -> 16 spaces 
                   label='                ' 
                   write(tname,'(i3)') lpj
               end if
!              formatted output by case of lpj
               select case(lpj)
               case(1)
!                  Note: r_l2 has no meaning when lpj is one.
!                  Keep it, though, as we may want to add a convention later. 
                   write(13,'(2x,f16.8,a,f16.8,6x,a)')  &
                              r_l(l),trim(tname), havg(1)  &
                             ,il(l)//'-projector'
                   write( 6,'(2x,f16.8,i3,f16.8,6x,a)')  &
                              r_l(l),lpj, havg(1)  &
                             ,il(l)//'-projector'
                   if (l.gt.1)write(13,'(21x,f16.8)')hso(1)
                   if (l.gt.1)write( 6,'(21x,f16.8)')hso(1)
               case(2)
                   write(13,'(2x,f16.8,a,2f16.8,6x,a)')  &
                              r_l(l),trim(tname), havg(1:2)  &
                             ,il(l)//'-projector'
                   write( 6,'(2x,f16.8,a,2f16.8,6x,a)')  &
                              r_l(l),trim(tname), havg(1:2)  &
                             ,il(l)//'-projector'
                   write(13,'(2x,a,19x,f16.8)')trim(label),havg(3)
                   write( 6,'(2x,a,19x,f16.8)')trim(label),havg(3)
                   if (l.gt.1)then
                     write(13,'(21x,2f16.8)') hso(1:2)
                     write( 6,'(21x,2f16.8)') hso(1:2)
                     write(13,'(37x, f16.8)') hso(3)
                     write( 6,'(37x, f16.8)') hso(3)
                   end if
               case(3)
                   write(13,'(2x,f16.8,a,3f16.8,6x,a)')  &
                              r_l(l),trim(tname),havg(1:2), havg(4)  &
                             ,il(l)//'-projector'
                   write( 6,'(2x,f16.8,a,3f16.8,6x,a)')  &
                              r_l(l),trim(tname),havg(1:2), havg(4)  &
                             ,il(l)//'-projector'
                   write(13,'(2x,a,19x,2f16.8)')  &
                             trim(label),havg(3),havg(5)
                   write( 6,'(2x,a,19x,2f16.8)')  &
                             trim(label),havg(3),havg(5)
                   write(13,'(53x,f16.8)') havg(6)
                   write( 6,'(53x,f16.8)') havg(6)
                   if (l.gt.1)then
                     write(13,'(21x,3f16.8)') hso(1:2),hso(4)
                     write( 6,'(21x,3f16.8)') hso(1:2),hso(4)
                     write(13,'(37x,2f16.8)') hso(3),hso(5)
                     write( 6,'(37x,2f16.8)') hso(3),hso(5)
                     write(13,'(53x, f16.8)') hso(6)
                     write( 6,'(53x, f16.8)') hso(6)
                   end if
               end select
!              dimension of hij
            end do
!           loop over l
            close(13)
!        iproc is zero
         end if
!        end of writing the psppar by process zero
         
!        nonlinear core correction file for BigDFT
         if(rcore>0d0.and.iproc==0)then
            open(13,file='nlcc')
            write(13,*)'0  no valence distribution to subtract'
            write(13,*)'1  one single gaussian for the core'
!           BigDFT has a different convention here,
!           the polynomial is not scaled by rcore.
            do i=1,4
               gcorepp(i)=gcore(i)*rcore**(2-2*i)
            end do
            write(13,*)rcore, gcorepp,  &
                       'sigma c0 c2 c4 c6'
            write(13,*)'end of NLCC input for BigDFT'
            write(13,*)'analytic form: gaussian(r,sigma)*sum(ck*x**k)'

            close(13)
!           print out the analytic charge of the resulting NLCC
            write(6,*)
            write(6,*)'Analytic core charge of the NLCC:',  &
            sqrt( 2.0d0*atan(1.0) )*(  &
                     gcorepp(1)*rcore**3  &
            +  3.0d0*gcorepp(2)*rcore**5  &
            + 15.0d0*gcorepp(3)*rcore**7  &
            +105.0d0*gcorepp(4)*rcore**9)
            fourpi = 16.d0*atan(1.d0)
            tt=0d0
            do k= 1,nint
                r2=(rr(k)/rcore)**2
                tt=tt+  &
                exp(-.5d0*r2)/fourpi *rw(k)  *(  &
                gcore(1)  &
              + gcore(2)*r2   &
              + gcore(3)*r2**2   &
              + gcore(4)*r2**3 )
            end do
            write(6,*)'Value for the radial grid used:  ',tt
            write(6,*)



!           and a plot to compare with the full AE core charge
            open(13,file='nlcc.gnuplt')
            write(13,*)'rcore=',rcore
            write(13,*)'c0=',gcore(1)
            write(13,*)'c2=',gcore(2)
            write(13,*)'c4=',gcore(3)
            write(13,*)'c6=',gcore(4)
            write(13,*)'r(x)=x/rcore'
            write(13,*)'p(x)=(c0 +c2*x*x +c4*x**4 +c6*x**6)'
            write(13,*)'g(x)=exp(-0.5*x**2)'
            write(13,*)'rho(x)=p(r(x))*g(r(x))'
            write(13,*)"set xrange [0.2*rcore:5.0*rcore]"
            write(13,*)"  p rho(x)*x*x"
            write(13,*)"rep 'ae.core.dens.plt','ae.core.dens.plt' u 1:3"
            write(13,*)"show function"
            close(13)
         end if




!  dumpfile for testing with another program
         if(ldump.and.iproc==0)then
            write(6,*)'Writing out a dumpfile of',  &
              8*4+size(xp)+size(psi)+size(occup),'byte'           
            open(13,file='dumpfile.bin',form='unformatted')
!           open(13,file='dumpfile')
            write(13)ng,noccmax,lmax,lpx,  &
                     lcx,nspin,xp,psi,occup
            close(13)
         end if
!
!     here we used to overwrite old values of 'psp.par' with the current ones
!     There was an info flag used to append some gaussian coeffs to psp.par
         if (iproc==0.and.namoeb.ne.0) then
            if (info) then
               open(unit=23,file='gauss.par',form='formatted')
               write(23,*) 'Additional information (last calculation):'
               write(23,*) 'gaussian exponents:'
               write(23,*) 'xp():',(xp(i),i=1,ng)
               write(23,*) 'orbital coefficients ',  &
                    '(one orbital per column)'
               do l=0,lmax
!        no spin down s orbitals in the r case: nspin=2, nspol=1, 2l+1=1
                  do ispin=1,max(nspol,min(2*l+1,nspin))
                     write(23,*) 'l,ispin:',l,ispin
                     write(23,*) 'psi n=1,noccmax(l):'
                     do i=0,ng
                        write(23,'(10e20.10)')  &
                             (psi(i,nocc,l+1,ispin),  &
                             nocc=1,nomax(l))
                     enddo
                  enddo
               enddo
            endif
            close(unit=23)
         endif


!
!     PLOT WAVEFUNCTIONS (up to 5*rcov)   c.hartwigsen: version for gnuplot
!

!     NOTE: For now, only do this with config0 from atom.00.ae
!           should be straight forward to write out all plots in parallel


      if (plotwf.and.iproc==0) then
!     generate various files for plotting

! write out the local potential in real space
        open(17,file='local.pot')
        write(17,'(a)')'# r, Vloc-Vnuc, Vnuc=-z/r, erf term, gpot term'
        do k=1,nint
              r=rr(k)
              gt=exp(-.5d0*(r/rloc)**2)*  &
                        (gpot(1) + gpot(2)*(r/rloc)**2+    &
                         gpot(3)*(r/rloc)**4 +  &
                         gpot(4)*(r/rloc)**6 )
              et= -zion*Derf(r/(sqrt(2.d0)*rloc))/r

              write(17,'(5e18.8)') r,et+gt+zion/r,zion/r, et, gt
                          !.5d0*(r/rprb**2)**2
        enddo
        close(17)

!       write out some plotable kernel for the separable part.
!       For now let us try the diagonal part in real space
!       and the contour at r'=rcov for each l-component.
!       ignore the case kij /= 0 
        do l=1,lpx
           if(iproc>1) exit
           lq=l-1
           rnrm1=1.d0/sqrt(.5d0*gamma(lq+1.5d0)*r_l(l)**(2*lq+3))
           rnrm2=1.d0/sqrt(.5d0*gamma(lq+3.5d0)*r_l2(l)**(2*lq+7))
           rnrm3=1.d0/sqrt(.5d0*gamma(lq+5.5d0)*r_l2(l)**(2*lq+11))
           open(17,file=trim(il(l))//'.kernel.pot')
           write(17,'(3(9x,a))')'#   r          ',&         
                                 'V_'//il(l)//'(r,r)   ',&
                                 'V_'//il(l)//'(r,rcov)'
           do k=1,nint
              r=rr(k)
              ppr1=rnrm1*r**lq    *exp(-.5d0*(r/r_l(l))**2)
              ppr2=rnrm2*r**(lq+2)*exp(-.5d0*(r/r_l2(l))**2)
              ppr3=rnrm3*r**(lq+4)*exp(-.5d0*(r/r_l2(l))**2)
              ppc1=rnrm1*rcov**lq    *exp(-.5d0*(r/r_l(l))**2)
              ppc2=rnrm2*rcov**(lq+2)*exp(-.5d0*(rcov/r_l2(l))**2)
              ppc3=rnrm3*rcov**(lq+4)*exp(-.5d0*(rcov/r_l2(l))**2)
              write(17,'(3f20.5)')r, &
                            ppr1*hsep(1,l,1)*ppr1   + &
                            ppr1*hsep(2,l,1)*ppr2*2 + &  
                            ppr2*hsep(3,l,1)*ppr2   + &
                            ppr1*hsep(4,l,1)*ppr3*2 + &
                            ppr2*hsep(5,l,1)*ppr3*2 + &  
                            ppr3*hsep(6,l,1)*ppr3 , &
!                           ---------------------------
                            ppc1*hsep(1,l,1)*ppr1   + &
                            ppc1*hsep(2,l,1)*ppr2*2 + &  
                            ppc2*hsep(3,l,1)*ppr2   + &
                            ppc1*hsep(4,l,1)*ppr3*2 + &
                            ppc2*hsep(5,l,1)*ppr3*2 + &  
                            ppc3*hsep(6,l,1)*ppr3 
           end do
           close(17)
        end do

!       plotting of the orbitals
            call detnp(ngrid,rae,5*rcov,np)
            open(32,file='pswf.gnu',form='formatted',status='unknown')
            write (32,*) 'set data style lines'
            do iorb=1,norb
               nocc=no(iorb)
               l=lo(iorb)
               ispin=1
               if(so(iorb).lt.0) ispin=2
               if (ispp.eq.'r')then
                  fname= 'ps.'//char(ichar('0')+noae(iorb))  &
                       //il(lo(iorb)+1)  &
                       //char(ichar('0')+int(2*(lo(iorb)+so(iorb))))  &
                       //'by2.dat'
                  tname=char(ichar('0')+noae(iorb))//il(lo(iorb)+1)  &
                       //char(ichar('0')+int(2*(lo(iorb)+so(iorb))))  &
                       //'/2'
               else
!                 either nonrel or spin pol
                  fname= 'ps.'//char(ichar('0')+noae(iorb))  &
                       //il(lo(iorb)+1)
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
               open(unit=2,file=trim(fname),  &
                       form='formatted',status='unknown')
!     find outer max of psi (approx), search from 10 bohr down
                  ra=10.d0
                  ttrmax=ra
                  ttmax= dabs(wave(ng,l,xp,psi(0,nocc,l+1,ispin),ra))
                  do i=100,0, -1
                     ra= 0.1d0 * i
                     ttpsi=dabs(wave(ng,l,xp,psi(0,nocc,l+1,ispin),ra))
!                     print*,ra,ttpsi
                     if ( ttpsi .gt. ttmax  &
                          .and. ttpsi .gt. 1.0d-4 ) then
                        ttmax=ttpsi
                        ttrmax=ra
                     endif
                  if (ttpsi.lt.ttmax .and. ttpsi.gt.1.0d-4) goto 3456
                  enddo
 3456             continue
!     ae/pseudo wfs should have the same sign for large r when plotted
               call detnp(ngrid,rae,ttrmax,nsign)
!     never use first gridpoint! (only relevant for H and He)
               if (nsign.eq.1) nsign=nsign+1
               tt=psiold(nsign,nocc,l+1,ispin)
               sign1=tt/abs(tt)
               tt= wave(ng,l,xp,psi(0,nocc,l+1,ispin),rae(nsign))
               sign2=tt/abs(tt)
               do i=2,np
                  ttold=psiold(i,nocc,l+1,ispin)*sign1*rae(i)
                  ttold=max(min(3.d0,ttold),-3.d0)
                  tt=wave(ng,l,xp,psi(0,nocc,l+1,ispin),rae(i))  &
                       *sign2*rae(i)
                  tt=max(min(3.d0,tt),-3.d0)
                  ttdiff=psiold(i,nocc,l+1,ispin)*sign1-  &
                       wave(ng,l,xp,psi(0,nocc,l+1,ispin),rae(i))*sign2
                  ttdiff= ttdiff*rae(i)
                  ttdiff=log(max(abs(ttdiff),1.d-8))/log(10.d0)
                  write(2,'(7e20.10)') rae(i),ttold,tt,ttdiff
!     plot of the wavefunction and the higher derivatives
!                  write(2,'(7e20.10)') rae(i),
!     :                 psiold(i,nocc,l+1,ispin),
!     :                 wave(ng,l,xp,psi(0,nocc,l+1,ispin),rae(i)),
!     :                 dwave(ng,l,xp,psi(0,nocc,l+1,ispin),rae(i)),
!     :                 ddwave(ng,l,xp,psi(0,nocc,l+1,ispin),rae(i))
               enddo
               close(2)
               write (32,*)'  plot "'//fname(1:index(fname,' ')-1)  &
                    //'"     title "'//tname(1:index(tname,' ')-1)  &
                    //'"'
               write (32,*)'replot "'//fname(1:index(fname,' ')-1)  &
                    //'" using 1:3 title "pseudo"'
               write(32,*) 'pause -10 "Hit return to continue"'
               write (32,*)'replot "'//fname(1:index(fname,' ')-1)  &
                    //'" using 1:4 title "diff"'
               write(32,*) 'pause -10 "Hit return to continue"'
            enddo
            write (32,*) 'set nokey'
            close(unit=32)
            write(6,*)
            write(6,*)'Plot files for gnuplot written.'
            write(6,*)'To plot wfs type ''gnuplot pswf.gnu'' '
            write(6,*)
         endif



!    -----------------------------------------------
!                     MAIN LOOP END
      if(iiter>namoeb)exit
      enddo
 1000 continue


!     test orthogonality of the projectors
!
      if (lpx.ge.1.and.ortprj)  &
           call pj2test(hsep,lpx-1,lpmx,lmx,nspin,nsmx,r_l,is)
!
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
!     end


!     optional info when call to the shell is available
      write(label,'(2a,i3.3,a)')'echo finished at `date`',  &
                                '>> hostnames.',iproc
      call system(label)
      end program

!
!     CRAY: no derf() -> user erf()
!
!      real*8 function derf(x)
!      REAL*8 X
!      DERF=ERF(X)
!      RETURN
!      END

