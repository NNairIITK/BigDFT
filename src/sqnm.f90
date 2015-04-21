!> @file
!!  Routines for Stefan's new minimization method
!! @author Stefan Goedecker and Bastian Schaefer
!! @section LICENCE
!!    Copyright (C) 2014 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS
!subroutine geopt(nat,wpos,etot,fout,fnrmtol,count,count_sd,displr)
subroutine sqnm(runObj,outsIO,nproc,iproc,verbosity,ncount_bigdft,fail)
!bigdft_state has to be run once on runObj and outs
!before calling this routine. sqnm will return to the caller the energies
!and coordinates used/obtained at the last accepted iteration step
   use module_base
   use bigdft_run!module_types
   use yaml_output
   use module_sqn, only: modify_gradient, getSubSpaceEvecEval, findbonds
   implicit none
   !parameter
   integer, intent(in)                    :: nproc
   integer, intent(in)                    :: iproc
   integer, intent(in)                    :: verbosity
   type(run_objects), intent(inout)       :: runObj
   type(state_properties), intent(inout) :: outsIO
   integer, intent(inout)                 :: ncount_bigdft
   logical, intent(out)                   :: fail
   !local variables
   character(len=*), parameter :: subname='sqnm'
   integer :: infocode,info !< variables containing state codes
   integer :: nhistx !< maximum history length
   integer :: nhist  !< actual history length
   integer :: ndim   !< dimension of significant subspace
   integer :: nit    !< maximum number of iterations
   integer :: nat    !< number of atoms
   integer :: istat,iall
   integer :: lwork
   integer :: it,i,iat,l,j,idim,jdim,ihist,icheck !<counter variables
   integer :: itswitch
   integer :: imode
   integer :: nbond
   type(state_properties) :: outs
   logical :: success=.false.
   logical :: debug !< set .true. for debug output to fort.100
   logical :: steep !< steepest descent flag
   real(gp) :: displr !< (non-physical) integrated path length,
                      !< includes displacements from rejctions
                      !< (denoted as dsplr in geopt.mon)
   real(gp) :: displp !< (physical) integrated path length,
                      !< includes NO displacements from rejections
                      !< (denoted as dsplp in geopt.mon)
   real(gp) :: etot 
   real(gp) :: fnrm
   real(gp) :: fmax
   real(gp) :: fluct
   real(gp) :: fnoise
   real(gp) :: betax !< initial step size (gets not changed)
   real(gp) :: beta_stretchx
   real(gp) :: beta  !< current step size
   real(gp) :: beta_stretch  !< current step size in bond-stretching directions
   real(gp) :: cosangle
   real(gp) :: tt
   real(gp) :: maxd !< maximum displacement of single atom
   real(gp) :: scl
   real(gp) :: etotold
   real(gp) :: detot
   real(gp) :: etotp
   real(gp) :: dt
   real(gp) :: cutoffRatio !< if fraction of eigenvalues of overlapmatrix of
                           !< displacements is smaller that cutoffRatio,
                           !< then those displacements are regarded
                           !< as linear dependent and are not taken into account
                           !< for building the significant subspace
   real(gp) :: maxrise !< energy ist not allowed to rise more than maxrise in single step
   real(gp) :: steepthresh !< if fnrm is larger that steepthresh, steepest descent is used
   real(gp) :: trustr !< a single atoms is not allowed to be dsiplaced more than by trustr
   integer :: idxtmp
   integer,  allocatable, dimension(:,:)   :: iconnect
   integer,  allocatable, dimension(:)     :: idx!index array for keeping track of history
   real(gp), allocatable, dimension(:,:,:) :: rxyz
   real(gp), allocatable, dimension(:,:,:) :: fxyz
   real(gp), allocatable, dimension(:,:,:) :: fxyzraw
   real(gp), allocatable, dimension(:,:,:) :: fstretch
   real(gp), allocatable, dimension(:,:)   :: rxyzOld
   real(gp), allocatable, dimension(:,:)   :: delta
   real(gp), allocatable, dimension(:,:,:) :: ff
   real(gp), allocatable, dimension(:,:,:) :: rr
   real(gp), allocatable, dimension(:,:,:) :: rrr
   real(gp), allocatable, dimension(:,:,:) :: fff
   real(gp), allocatable, dimension(:,:)   :: aa
   real(gp), allocatable, dimension(:,:)   :: dd
   real(gp), allocatable, dimension(:)     :: eval
   real(gp), allocatable, dimension(:)     :: work
   real(gp), allocatable, dimension(:)     :: res
   real(gp), allocatable, dimension(:)     :: scpr
   real(gp), allocatable, dimension(:)     :: rnorm
   real(gp), allocatable, dimension(:)     :: wold
   real(gp), allocatable, dimension(:)     :: rcov
   character(len=4)                        :: fn4
   character(len=40)                       :: comment
   character(len=12)                        :: cdmy12_1
   character(len=12)                        :: cdmy12_2
   character(len=9)                        :: cdmy9
   character(len=8)                        :: cdmy8
   integer :: ifail
   !functions
   real(gp) :: dnrm2
   type(f_tree) :: f_info

   f_info=f_tree_new()

   !set parameters
   nit=runObj%inputs%ncount_cluster_x
   nat=runObj%atoms%astruct%nat
   betax=runObj%inputs%betax
   nhistx=runObj%inputs%nhistx
   maxrise=runObj%inputs%maxrise
   cutoffRatio=runObj%inputs%cutoffratio
   steepthresh=runObj%inputs%steepthresh
   trustr=runObj%inputs%trustr
   nbond=1
   ifail=0
   if(runObj%inputs%biomode)imode=2

   if (iproc==0.and.verbosity > 0) then
      call yaml_mapping_open('Geometry parameters')
         call yaml_map('Geometry Method','GEOPT_SQNM')
         call yaml_map('nhistx',nhistx)
         call yaml_map('biomode',runObj%inputs%biomode)
         call yaml_map('betax', betax,fmt='(1pe21.14)')
         call yaml_map('beta_stretchx', runObj%inputs%beta_stretchx,fmt='(1pe21.14)')
         call yaml_map('maxrise', maxrise,fmt='(1pe21.14)')
         call yaml_map('cutoffRatio', cutoffRatio,fmt='(1pe21.14)')
         call yaml_map('steepthresh', steepthresh,fmt='(1pe21.14)')
         call yaml_map('trustr', trustr,fmt='(1pe21.14)')
      call yaml_mapping_close()
   end if

   !init variables
   debug=.false.
   fail=.true.
   displr=0.0_gp
   displp=0.0_gp
   fluct=0.0_gp
   icheck=0
   detot=0.0_gp
   itswitch=0
   ndim=0
   nhist=0
   beta=betax
   beta_stretch=runObj%inputs%beta_stretchx
   maxd=0.0_gp

   ! allocate arrays
   lwork=1000+10*nat**2
   rxyz = f_malloc((/ 1.to.3, 1.to.nat, 0.to.nhistx /),id='rxyz')
   idx = f_malloc((/ 0.to.nhistx /),id='ixyz')
   fxyz = f_malloc((/ 1.to.3, 1.to.nat, 0.to.nhistx /),id='fxyz')
   fxyzraw = f_malloc((/ 1.to.3, 1.to.nat, 0.to.nhistx /),id='fxyzraw')
   fstretch = f_malloc((/ 1.to.3, 1.to.nat, 0.to.nhistx /),id='fstretch')
   rxyzOld = f_malloc((/ 1.to.3, 1.to.nat/),id='rxyzOld')
   delta = f_malloc((/ 1.to.3, 1.to.nat/),id='delta')
   aa = f_malloc((/ nhistx, nhistx /),id='aa')
   eval = f_malloc(nhistx,id='eval')
   res = f_malloc(nhistx,id='res')
   rnorm = f_malloc(nhistx,id='rnorm')
   work = f_malloc(lwork,id='work')
   ff = f_malloc((/ 1.to.3, 1.to.nat, 0.to.nhistx /),id='ff')
   rr = f_malloc((/ 1.to.3, 1.to.nat, 0.to.nhistx /),id='rr')
   dd = f_malloc((/ 3, nat /),id='dd')
   fff = f_malloc((/ 1.to.3, 1.to.nat, 0.to.nhistx /),id='fff')
   rrr = f_malloc((/ 1.to.3, 1.to.nat, 0.to.nhistx /),id='rrr')
   scpr = f_malloc(nhistx,id='scpr')
   rcov     = f_malloc((/ 1.to.nat/),id='rcov')
   iconnect = f_malloc((/ 1.to.2, 1.to.1000/),id='iconnect')
   if(runObj%inputs%biomode)then
        call give_rcov_sqnm(iproc,runObj%atoms,runObj%atoms%astruct%nat,rcov)
        call findbonds('(SQNM)',iproc,10,runObj%atoms%astruct%nat,&
             rcov,runObj%atoms%astruct%rxyz,nbond,iconnect)
   endif 
   wold = f_malloc((/ 1.to.nbond/),id='wold')
   wold =0.0_gp
   do i=0,nhistx
    idx(i)=i
   enddo


   call init_state_properties(outs, runObj%atoms%astruct%nat)

   !copy outs_datatype
   call copy_state_properties(outsIO,outs)



!!!!!!   call energyandforces(nat,rxyz(1,1,0),fxyz(1,1,0),etot)
!!  not necessary, bigdft_state allready called outside
!   call bigdft_state(runObj,outs,nproc,iproc,infocode)
!   ncount_bigdft=ncount_bigdft+1

!! copy to internal variables
   call vcopy(3*runObj%atoms%astruct%nat, runObj%atoms%astruct%rxyz(1,1), 1,rxyz(1,1,idx(0)), 1)
   call vcopy(3*runObj%atoms%astruct%nat, runObj%atoms%astruct%rxyz(1,1), 1,rxyzOld(1,1), 1)
   call vcopy(3*outs%fdim, outs%fxyz(1,1), 1, fxyz(1,1,idx(0)), 1)
   etot=outs%energy

   call minenergyandforces(iproc,nproc,.false.,imode,runObj,outs,nat,rxyz(1,1,idx(0)),&
       fxyz(1,1,idx(0)),fstretch(1,1,idx(0)),fxyzraw(1,1,idx(0)),&
       etot,iconnect,nbond,wold,beta_stretchx,beta_stretch,infocode)
   if(imode==2)rxyz(:,:,idx(0))=rxyz(:,:,idx(0))+beta_stretch*fstretch(:,:,idx(0))

   call fnrmandforcemax(fxyzraw(1,1,idx(0)),fnrm,fmax,nat)
   fnrm=sqrt(fnrm)
   if (fmax < 3.e-1_gp) call updatefluctsum(outs%fnoise,fluct)

   etotold=etot
   etotp=etot

   if (iproc==0.and.verbosity > 0) then
       !avoid space for leading sign (numbers are positive, anyway)
       write(cdmy8,'(es8.1)')abs(maxd)
       write(cdmy12_1,'(es12.5)')abs(displr)
       write(cdmy12_2,'(es12.5)')abs(displp)
       write(cdmy9,'(es9.2)')abs(beta)


       write(16,'(i5,1x,i5,2x,a10,2x,1es21.14,2x,es9.2,es11.3,3es10.2,2x,a6,a8,1x,a4,i3.3,1x,a5,a7,2(1x,a6,a11))') &
       ncount_bigdft,0,'GEOPT_SQNM',etotp,detot,fmax,fnrm,fluct*runObj%inputs%frac_fluct,fluct, &
       'beta=',trim(adjustl(cdmy9)),'dim=',ndim,'maxd=',&
       trim(adjustl(cdmy8)),'dsplr=',trim(adjustl(cdmy12_1)),&
       'dsplp=',trim(adjustl(cdmy12_2))
       call f_utils_flush(16)
   endif

   do it=1,nit!start main loop
!  do it=1,nit-1!start main loop (nit-1 if first bigdft call is NOT done outside, but inside this subroutine)
      if (debug.and.iproc==0) write(100,*) 'it:',it,etot,fnrm,itswitch
      nhist=nhist+1

      if (fnrm.gt.steepthresh .or. it.le.itswitch ) then
         ndim=0
         steep=.true.
         if (it.gt.itswitch) itswitch=it+nhistx
         if (debug.and.iproc==0) write(100,*) "STEEP"
      else
         steep=.false.
      endif

      ! cycle the index array
      if (nhist.gt.nhistx) then
         nhist=nhistx
         idxtmp=idx(0)
         do ihist=0,nhist-1
            idx(ihist)=idx(ihist+1)
         enddo
         idx(nhist)=idxtmp
      endif
   
      ! decompose gradient
500 continue
    call modify_gradient(nat,ndim,rrr(1,1,1),eval(1),res(1),fxyz(1,1,idx(nhist-1)),beta,dd(1,1))
   
      tt=0.0_gp
      dt=0.0_gp
      maxd=-huge(1.0_gp)
      do iat=1,nat
         dt=dd(1,iat)**2+dd(2,iat)**2+dd(3,iat)**2
         tt=tt+dt
         maxd=max(maxd,dt)
      enddo
      tt=sqrt(tt)
      maxd=sqrt(maxd)
   
      !trust radius approach: avoids too large steps due to large forces
      !only used when in steepest decent mode
      if(maxd>trustr .and. steep)then
         if(debug.and.iproc==0)write(100,'(a,1x,es24.17,1x,i0)')'step too large',maxd,it
         if(iproc==0)then
            write(16,'(a,2(1x,es9.2))')'WARNING GEOPT_SQNM: step too large: maxd, trustradius ',maxd,trustr
            call f_utils_flush(16)
         endif
         scl=0.50_gp*trustr/maxd
         dd=dd*scl
         tt=tt*scl
         maxd=maxd*scl
      endif
!      displr=displr+tt
   
      !update positions
      do iat=1,nat
         rxyz(1,iat,idx(nhist))=rxyz(1,iat,idx(nhist-1))-dd(1,iat)
         rxyz(2,iat,idx(nhist))=rxyz(2,iat,idx(nhist-1))-dd(2,iat)
         rxyz(3,iat,idx(nhist))=rxyz(3,iat,idx(nhist-1))-dd(3,iat)
      enddo
   
      delta=rxyz(:,:,idx(nhist))-rxyzOld
      displr=displr+dnrm2(3*nat,delta(1,1),1)
      runObj%inputs%inputPsiId=1
      call minenergyandforces(iproc,nproc,.true.,imode,runObj,outs,nat,rxyz(1,1,idx(nhist)),&
                             fxyz(1,1,idx(nhist)),fstretch(1,1,idx(nhist)),fxyzraw(1,1,idx(nhist)),&
                             etotp,iconnect,nbond,wold,beta_stretchx,beta_stretch,infocode)
      detot=etotp-etotold
      ncount_bigdft=ncount_bigdft+1


      call fnrmandforcemax(fxyzraw(1,1,idx(nhist)),fnrm,fmax,nat)
      fnrm=sqrt(fnrm)

      if (iproc == 0) then
         write(fn4,'(i4.4)') ncount_bigdft
         write(comment,'(a,1pe10.3)')'SQNM:fnrm= ',fnrm
         call bigdft_write_atomic_file(runObj,outs,'posout_'//fn4,&
              trim(comment))
      endif
      if(infocode==0)then
        ifail=0
      else
        ifail=ifail+1
      endif
      if ((infocode==0 .or. ifail>20).and.(fmax < 3.e-1_gp)) call updatefluctsum(outs%fnoise,fluct)
      cosangle=-dot_double(3*nat,fxyz(1,1,idx(nhist)),1,dd(1,1),1)/&
              sqrt(dot_double(3*nat,fxyz(1,1,idx(nhist)),1,fxyz(1,1,idx(nhist)),1)*&
              dot_double(3*nat,dd(1,1),1,dd(1,1),1))

      if (detot.gt.maxrise .and. beta > 1.e-1_gp*betax) then !
         if (debug.and.iproc==0) write(100,'(a,i0,1x,e9.2)') "WARN: it,detot", it,detot
         if (debug.and.iproc==0) then
             write(16,'(a,i0,4(1x,e9.2))') &
             "WARNING GEOPT_SQNM: Prevent energy to rise by more than maxrise: it,maxrise,detot,beta,1.e-1*betax ",&
             it,maxrise,detot,beta,1.e-1_gp*betax
             call f_utils_flush(16)
         endif
         if (iproc==0.and.verbosity > 0) then
            !avoid space for leading sign (numbers are positive, anyway)
            write(cdmy8,'(es8.1)')abs(maxd)
            write(cdmy12_1,'(es12.5)')abs(displr)
            write(cdmy12_2,'(es12.5)')abs(displp)
            write(cdmy9,'(es9.2)')abs(beta)

   
            write(16,'(i5,1x,i5,2x,a10,2x,1es21.14,2x,es9.2,es11.3,3es10.2,2x,a6,a8,1x,a4,i3.3,1x,a5,a7,2(1x,a6,a11))') &
             ncount_bigdft,it,'GEOPT_SQNM',etotp,detot,fmax,fnrm,fluct*runObj%inputs%frac_fluct,fluct, &
             'beta=',trim(adjustl(cdmy9)),'dim=',ndim,&
             'maxd=',trim(adjustl(cdmy8)),'dsplr=',trim(adjustl(cdmy12_1)),&
             'dsplp=',trim(adjustl(cdmy12_2))
            call f_utils_flush(16)

            call f_tree_push(f_info//'ndim'        ,yaml_toa(ndim))
            call f_tree_push(f_info//'etot'        ,yaml_toa(etotp,fmt='(1pe21.14)'))
            call f_tree_push(f_info//'detot'       ,yaml_toa(detot,fmt='(1pe21.14)'))
            call f_tree_push(f_info//'beta'        ,yaml_toa(beta,fmt='(1pe21.14)'))
            call f_tree_push(f_info//'beta_stretch',yaml_toa(beta_stretch,fmt='(1pe21.14)'))
            call geometry_output('GEOPT_SQNM',ncount_bigdft,it,fmax,fnrm,fluct,f_info)

!!$            call yaml_mapping_open('Geometry')
!!$               call yaml_map('Ncount_BigDFT',ncount_bigdft) !universal
!!$               call yaml_map('Geometry step',it)
!!$               call yaml_map('Geometry Method','GEOPT_SQNM')
!!$               call yaml_map('ndim',ndim)
!!$               call yaml_map('etot', etotp,fmt='(1pe21.14)')
!!$               call yaml_map('detot',detot,fmt='(1pe21.14)')
!!$               call yaml_map('fmax',fmax,fmt='(1pe21.14)')
!!$               call yaml_map('fnrm',fnrm,fmt='(1pe21.14)')
!!$               call yaml_map('beta',beta,fmt='(1pe21.14)')
!!$               call yaml_map('beta_stretch',beta_stretch,fmt='(1pe21.14)')
!!$               call geometry_output(fmax,fnrm,fluct)
!!$            call yaml_mapping_close()
         end if
    
         if(ncount_bigdft >= nit)then!no convergence within ncount_cluster_x energy evaluations
            !following copy of rxyz(1,1,nhist-1) to runObj is necessary for returning to the caller
            !the energies and coordinates used/obtained from/in the last ACCEPTED iteration step
            !(otherwise coordinates of last call to bigdft_state would be returned)
            call vcopy(3 * runObj%atoms%astruct%nat, rxyz(1,1,idx(nhist-1)), 1,runObj%atoms%astruct%rxyz(1,1), 1)
            goto 900  !sqnm will return to caller the energies and coordinates used/obtained from the last ACCEPTED iteration step
         endif

         beta=.50_gp*beta
         if (debug.and.iproc==0) write(100,'(a,1x,e9.2)') 'WARNING GEOPT_SQNM: beta reset ',beta
         ndim=0
         wold=0.0_gp
         if(.not.steep)then
            do iat=1,nat
               rxyz(1,iat,0)=rxyz(1,iat,idx(nhist-1))
               rxyz(2,iat,0)=rxyz(2,iat,idx(nhist-1))
               rxyz(3,iat,0)=rxyz(3,iat,idx(nhist-1))
   
               fxyz(1,iat,0)=fxyz(1,iat,idx(nhist-1))
               fxyz(2,iat,0)=fxyz(2,iat,idx(nhist-1))
               fxyz(3,iat,0)=fxyz(3,iat,idx(nhist-1))
               fxyzraw(1,iat,0)=fxyzraw(1,iat,idx(nhist-1))
               fxyzraw(2,iat,0)=fxyzraw(2,iat,idx(nhist-1))
               fxyzraw(3,iat,0)=fxyzraw(3,iat,idx(nhist-1))
            enddo
            nhist=1
            do i=0,nhistx
                idx(i)=i
            enddo
         endif
         goto  500
      endif

      delta=rxyz(:,:,idx(nhist))-rxyzOld
      displp=displp+dnrm2(3*nat,delta(1,1),1)
      rxyzOld=rxyz(:,:,idx(nhist))
!      displp=displp+tt
      if (iproc==0.and.verbosity > 0) then
         !avoid space for leading sign (numbers are positive, anyway)
         write(cdmy8,'(es8.1)')abs(maxd)
         write(cdmy12_1,'(es12.5)')abs(displr)
         write(cdmy12_2,'(es12.5)')abs(displp)
         write(cdmy9,'(es9.2)')abs(beta)


         write(16,'(i5,1x,i5,2x,a10,2x,1es21.14,2x,es9.2,es11.3,3es10.2,2x,a6,a8,1x,a4,i3.3,1x,a5,a7,2(1x,a6,a11))') &
          ncount_bigdft,it,'GEOPT_SQNM',etotp,detot,fmax,fnrm,fluct*runObj%inputs%frac_fluct,fluct, &
          'beta=',trim(adjustl(cdmy9)),'dim=',ndim,'maxd=',&
          trim(adjustl(cdmy8)),'dsplr=',trim(adjustl(cdmy12_1)),&
          'dsplp=',trim(adjustl(cdmy12_2))
         call f_utils_flush(16)

         call f_tree_push(f_info//'ndim'        ,yaml_toa(ndim))
         call f_tree_push(f_info//'etot'        ,yaml_toa(etotp,fmt='(1pe21.14)'))
         call f_tree_push(f_info//'detot'       ,yaml_toa(detot,fmt='(1pe21.14)'))
         call f_tree_push(f_info//'beta'        ,yaml_toa(beta,fmt='(1pe21.14)'))
         call f_tree_push(f_info//'beta_stretch',yaml_toa(beta_stretch,fmt='(1pe21.14)'))
         call geometry_output('GEOPT_SQNM',ncount_bigdft,it,fmax,fnrm,fluct,f_info)

!!$         call yaml_mapping_open('Geometry')
!!$            call yaml_map('Ncount_BigDFT',ncount_bigdft)
!!$            call yaml_map('Geometry step',it)
!!$            call yaml_map('Geometry Method','GEOPT_SQNM')
!!$            call yaml_map('ndim',ndim)
!!$            call yaml_map('etot', etotp,fmt='(1pe21.14)')
!!$            call yaml_map('detot',detot,fmt='(1pe21.14)')
!!$            call yaml_map('fmax',fmax,fmt='(1pe21.14)')
!!$            call yaml_map('fnrm',fnrm,fmt='(1pe21.14)')
!!$            call yaml_map('beta',beta,fmt='(1pe21.14)')
!!$            call yaml_map('beta_stretch',beta_stretch,fmt='(1pe21.14)')
!!$            call geometry_output(fmax,fnrm,fluct)
!!$         call yaml_mapping_close()
      end if

      etot    = etotp
      etotold = etot
      !copy outs_datatype
      call copy_state_properties(outs,outsIO)

      if(detot .gt. maxrise)then
         if (iproc==0) then
            write(16,'(a,i0,4(1x,e9.2))') &
             "WARNING GEOPT_SQNM: Allowed energy to rise by more than maxrise: it,maxrise,detot,beta,1.d-1*betax ",&
             it,maxrise,detot,beta,1.e-1_gp*betax
            call f_utils_flush(16)
         endif
      endif


!      if (fnrm.le.fnrmtol) goto 1000
      call convcheck(fmax,fluct*runObj%inputs%frac_fluct,runObj%inputs%forcemax,icheck)
      if(icheck>5)then
         goto 1000
      endif
     if(imode==2)rxyz(:,:,idx(nhist))=rxyz(:,:,idx(nhist))+&
                 beta_stretch*fstretch(:,:,idx(nhist)) !has to be after convergence check,
                                                       !otherwise energy will not match
                                                       !the true energy of rxyz(:,:,idx(nhist))

      if(ncount_bigdft >= nit)then!no convergence within ncount_cluster_x energy evaluations
            goto 900  !sqnm will return to caller the energies and coordinates used/obtained from the last accepted iteration step
      endif
   
      if (cosangle.gt..200_gp) then
         beta=beta*1.100_gp
      else
         beta=max(beta*.850_gp,betax)
      endif
   
      if (debug.and.iproc==0) write(100,*) 'cosangle ',cosangle,beta

      call getSubSpaceEvecEval('(SQNM)',iproc,verbosity,nat,nhist,&
           nhistx,ndim,cutoffratio,lwork,work,idx,rxyz,fxyz,aa,rr,ff,&
           rrr,fff,eval,res,success)
      if(.not.success)stop 'subroutine minimizer_sqnm: no success in getSubSpaceEvecEval.'

   enddo!end main loop

900 continue

   !if code gets here, it failed
   if(debug.and.iproc==0) write(100,*) it,etot,fnrm
   if(iproc==0) then
        write(16,'(a,3(1x,i0))') &
       "WARNING GEOPT_SQNM: SQNM not converged: it,ncount_bigdft,ncount_cluster_x: ", &
       it,ncount_bigdft,runObj%inputs%ncount_cluster_x
       call f_utils_flush(16)
    endif
!   stop "No convergence "
   fail=.true.
   goto 2000

1000 continue!converged successfully
   
   if(iproc==0)then
         write(16,'(2(a,1x,i0))') "SQNM converged at iteration ",it,". Needed bigdft calls: ",ncount_bigdft
         call f_utils_flush(16)
   endif
   if(iproc==0)  call yaml_map('Iterations when SQNM converged',it)
   fail=.false.
   
2000 continue
!deallocations
   call f_free(rxyz)
   call f_free(idx)
   call f_free(rxyzOld)
   call f_free(delta)
   call f_free(fxyz)
   call f_free(fxyzraw)
   call f_free(fstretch)
   call f_free(aa)
   call f_free(eval)
   call f_free(res)
   call f_free(rnorm)
   call f_free(work)
   call f_free(ff)
   call f_free(rr)
   call f_free(dd)
   call f_free(fff)
   call f_free(rrr)
   call f_free(scpr)
   call f_free(wold)
   call f_free(rcov )   
   call f_free(iconnect)
   call f_tree_free(f_info)
   call deallocate_state_properties(outs)
end subroutine
subroutine minenergyandforces(iproc,nproc,eeval,imode,runObj,outs,nat,rat,fat,fstretch,&
           fxyzraw,epot,iconnect,nbond_,wold,alpha_stretch0,alpha_stretch,infocode)
    use module_base
    use bigdft_run!module_types
    use module_sqn
    use module_interfaces
    implicit none
    !parameter
    integer, intent(in)           :: iproc,nproc,imode
    integer, intent(in)           :: nat
    type(run_objects), intent(inout)       :: runObj
    type(state_properties), intent(inout) :: outs
    integer, intent(in)           :: nbond_
    integer, intent(in)           :: iconnect(2,nbond_)
    real(gp),intent(inout)        :: rat(3,nat)
    real(gp),intent(out)          :: fxyzraw(3,nat)
    real(gp),intent(inout)        :: fat(3,nat)
    real(gp),intent(out)          :: fstretch(3,nat)
    real(gp), intent(inout)       :: wold(nbond_)
    real(gp), intent(in)          :: alpha_stretch0
    real(gp), intent(inout)       :: alpha_stretch
    real(gp), intent(inout)       :: epot
    logical, intent(in)           :: eeval
    integer,intent(out) :: infocode
    !internal

    infocode=0
    if(eeval)then
        call vcopy(3 * runObj%atoms%astruct%nat, rat(1,1), 1,runObj%atoms%astruct%rxyz(1,1), 1)
        runObj%inputs%inputPsiId=1
        call bigdft_state(runObj,outs,infocode)
    endif
    call vcopy(3 * outs%fdim, outs%fxyz(1,1), 1, fat(1,1), 1)
    call vcopy(3 * outs%fdim, fat(1,1), 1,fxyzraw(1,1), 1)
    epot=outs%energy
    fstretch=0.0_gp
    if(imode==2)then
        call projectbond(nat,nbond_,rat,fat,fstretch,iconnect,&
             wold,alpha_stretch0,alpha_stretch)
    endif

end subroutine minenergyandforces
subroutine give_rcov_sqnm(iproc,atoms,nat,rcov)
  use module_base, only: gp
  use module_types
  use yaml_output
  implicit none
  !Arguments
  integer, intent(in) :: iproc,nat
  type(atoms_data), intent(in) :: atoms
  real(gp), intent(out) :: rcov(nat)
  !Local variables
  integer :: iat

  do iat=1,nat
     if (trim(atoms%astruct%atomnames(atoms%astruct%iatype(iat)))=='H') then
        rcov(iat)=0.75d0
!        rcov(iat)=0.75d0*0.529177211d0
     else if (trim(atoms%astruct%atomnames(atoms%astruct%iatype(iat)))=='LJ')then
        rcov(iat)=0.56d0
     else if (trim(atoms%astruct%atomnames(atoms%astruct%iatype(iat)))=='He')then
        rcov(iat)=0.75d0
     else if (trim(atoms%astruct%atomnames(atoms%astruct%iatype(iat)))=='Li')then
        rcov(iat)=3.40d0
     else if (trim(atoms%astruct%atomnames(atoms%astruct%iatype(iat)))=='Be')then
        rcov(iat)=2.30d0
     else if (trim(atoms%astruct%atomnames(atoms%astruct%iatype(iat)))=='B' )then
        rcov(iat)=1.55d0
     else if (trim(atoms%astruct%atomnames(atoms%astruct%iatype(iat)))=='C' )then
        rcov(iat)=1.45d0
!        rcov(iat)=1.45d0*0.529177211d0
     else if (trim(atoms%astruct%atomnames(atoms%astruct%iatype(iat)))=='N' )then
        rcov(iat)=1.42d0
!        rcov(iat)=1.42d0*0.529177211d0
     else if (trim(atoms%astruct%atomnames(atoms%astruct%iatype(iat)))=='O' )then
        rcov(iat)=1.38d0
!        rcov(iat)=1.38d0*0.529177211d0
     else if (trim(atoms%astruct%atomnames(atoms%astruct%iatype(iat)))=='F' )then
        rcov(iat)=1.35d0
     else if (trim(atoms%astruct%atomnames(atoms%astruct%iatype(iat)))=='Ne')then
        rcov(iat)=1.35d0
     else if (trim(atoms%astruct%atomnames(atoms%astruct%iatype(iat)))=='Na')then
        rcov(iat)=3.40d0
     else if (trim(atoms%astruct%atomnames(atoms%astruct%iatype(iat)))=='Mg')then
        rcov(iat)=2.65d0
     else if (trim(atoms%astruct%atomnames(atoms%astruct%iatype(iat)))=='Al')then
        rcov(iat)=2.23d0
     else if (trim(atoms%astruct%atomnames(atoms%astruct%iatype(iat)))=='Si')then
        rcov(iat)=2.09d0
     else if (trim(atoms%astruct%atomnames(atoms%astruct%iatype(iat)))=='P' )then
        rcov(iat)=2.00d0
     else if (trim(atoms%astruct%atomnames(atoms%astruct%iatype(iat)))=='S' )then
        rcov(iat)=1.92d0
     else if (trim(atoms%astruct%atomnames(atoms%astruct%iatype(iat)))=='Cl')then
        rcov(iat)=1.87d0
     else if (trim(atoms%astruct%atomnames(atoms%astruct%iatype(iat)))=='Ar')then
        rcov(iat)=1.80d0
     else if (trim(atoms%astruct%atomnames(atoms%astruct%iatype(iat)))=='K' )then
        rcov(iat)=4.00d0
     else if (trim(atoms%astruct%atomnames(atoms%astruct%iatype(iat)))=='Ca')then
        rcov(iat)=3.00d0
     else if (trim(atoms%astruct%atomnames(atoms%astruct%iatype(iat)))=='Sc')then
        rcov(iat)=2.70d0
     else if (trim(atoms%astruct%atomnames(atoms%astruct%iatype(iat)))=='Ti')then
        rcov(iat)=2.70d0
     else if (trim(atoms%astruct%atomnames(atoms%astruct%iatype(iat)))=='V' )then
        rcov(iat)=2.60d0
     else if (trim(atoms%astruct%atomnames(atoms%astruct%iatype(iat)))=='Cr')then
        rcov(iat)=2.60d0
     else if (trim(atoms%astruct%atomnames(atoms%astruct%iatype(iat)))=='Mn')then
        rcov(iat)=2.50d0
     else if (trim(atoms%astruct%atomnames(atoms%astruct%iatype(iat)))=='Fe')then
        rcov(iat)=2.50d0
     else if (trim(atoms%astruct%atomnames(atoms%astruct%iatype(iat)))=='Co')then
        rcov(iat)=2.40d0
     else if(trim(atoms%astruct%atomnames(atoms%astruct%iatype(iat)))=='Ni')then
        rcov(iat)=2.30d0
     else if (trim(atoms%astruct%atomnames(atoms%astruct%iatype(iat)))=='Cu')then
        rcov(iat)=2.30d0
     else if (trim(atoms%astruct%atomnames(atoms%astruct%iatype(iat)))=='Zn')then
        rcov(iat)=2.30d0
     else if (trim(atoms%astruct%atomnames(atoms%astruct%iatype(iat)))=='Ga')then
        rcov(iat)=2.10d0
     else if (trim(atoms%astruct%atomnames(atoms%astruct%iatype(iat)))=='Ge')then
        rcov(iat)=2.40d0
     else if (trim(atoms%astruct%atomnames(atoms%astruct%iatype(iat)))=='As')then
        rcov(iat)=2.30d0
     else if (trim(atoms%astruct%atomnames(atoms%astruct%iatype(iat)))=='Se')then
        rcov(iat)=2.30d0
     else if (trim(atoms%astruct%atomnames(atoms%astruct%iatype(iat)))=='Br')then
        rcov(iat)=2.20d0
     else if (trim(atoms%astruct%atomnames(atoms%astruct%iatype(iat)))=='Kr')then
        rcov(iat)=2.20d0
     else if (trim(atoms%astruct%atomnames(atoms%astruct%iatype(iat)))=='Rb')then
        rcov(iat)=4.50d0
     else if (trim(atoms%astruct%atomnames(atoms%astruct%iatype(iat)))=='Sr')then
        rcov(iat)=3.30d0
     else if (trim(atoms%astruct%atomnames(atoms%astruct%iatype(iat)))=='Y' )then
        rcov(iat)=3.30d0
     else if (trim(atoms%astruct%atomnames(atoms%astruct%iatype(iat)))=='Zr')then
        rcov(iat)=3.00d0
     else if (trim(atoms%astruct%atomnames(atoms%astruct%iatype(iat)))=='Nb')then
        rcov(iat)=2.92d0
     else if (trim(atoms%astruct%atomnames(atoms%astruct%iatype(iat)))=='Mo')then
        rcov(iat)=2.83d0
     else if (trim(atoms%astruct%atomnames(atoms%astruct%iatype(iat)))=='Tc')then
        rcov(iat)=2.75d0
     else if (trim(atoms%astruct%atomnames(atoms%astruct%iatype(iat)))=='Ru')then
        rcov(iat)=2.67d0
     else if (trim(atoms%astruct%atomnames(atoms%astruct%iatype(iat)))=='Rh')then
        rcov(iat)=2.58d0
     else if (trim(atoms%astruct%atomnames(atoms%astruct%iatype(iat)))=='Pd')then
        rcov(iat)=2.50d0
     else if (trim(atoms%astruct%atomnames(atoms%astruct%iatype(iat)))=='Ag')then
        rcov(iat)=2.50d0
     else if (trim(atoms%astruct%atomnames(atoms%astruct%iatype(iat)))=='Cd')then
        rcov(iat)=2.50d0
     else if (trim(atoms%astruct%atomnames(atoms%astruct%iatype(iat)))=='In')then
        rcov(iat)=2.30d0
     else if (trim(atoms%astruct%atomnames(atoms%astruct%iatype(iat)))=='Sn')then
        rcov(iat)=2.66d0
     else if (trim(atoms%astruct%atomnames(atoms%astruct%iatype(iat)))=='Sb')then
        rcov(iat)=2.66d0
     else if (trim(atoms%astruct%atomnames(atoms%astruct%iatype(iat)))=='Te')then
        rcov(iat)=2.53d0
     else if (trim(atoms%astruct%atomnames(atoms%astruct%iatype(iat)))=='I' )then
        rcov(iat)=2.50d0
     else if (trim(atoms%astruct%atomnames(atoms%astruct%iatype(iat)))=='Xe')then
        rcov(iat)=2.50d0
     else if (trim(atoms%astruct%atomnames(atoms%astruct%iatype(iat)))=='Cs')then
        rcov(iat)=4.50d0
     else if (trim(atoms%astruct%atomnames(atoms%astruct%iatype(iat)))=='Ba')then
        rcov(iat)=4.00d0
     else if (trim(atoms%astruct%atomnames(atoms%astruct%iatype(iat)))=='La')then
        rcov(iat)=3.50d0
     else if (trim(atoms%astruct%atomnames(atoms%astruct%iatype(iat)))=='Ce')then
        rcov(iat)=3.50d0
     else if (trim(atoms%astruct%atomnames(atoms%astruct%iatype(iat)))=='Pr')then
        rcov(iat)=3.44d0
     else if (trim(atoms%astruct%atomnames(atoms%astruct%iatype(iat)))=='Nd')then
        rcov(iat)=3.38d0
     else if (trim(atoms%astruct%atomnames(atoms%astruct%iatype(iat)))=='Pm')then
        rcov(iat)=3.33d0
     else if (trim(atoms%astruct%atomnames(atoms%astruct%iatype(iat)))=='Sm')then
        rcov(iat)=3.27d0
     else if (trim(atoms%astruct%atomnames(atoms%astruct%iatype(iat)))=='Eu')then
        rcov(iat)=3.21d0
     else if (trim(atoms%astruct%atomnames(atoms%astruct%iatype(iat)))=='Gd')then
        rcov(iat)=3.15d0
     else if (trim(atoms%astruct%atomnames(atoms%astruct%iatype(iat)))=='Td')then
        rcov(iat)=3.09d0
     else if (trim(atoms%astruct%atomnames(atoms%astruct%iatype(iat)))=='Dy')then
        rcov(iat)=3.03d0
     else if (trim(atoms%astruct%atomnames(atoms%astruct%iatype(iat)))=='Ho')then
        rcov(iat)=2.97d0
     else if (trim(atoms%astruct%atomnames(atoms%astruct%iatype(iat)))=='Er')then
        rcov(iat)=2.92d0
     else if (trim(atoms%astruct%atomnames(atoms%astruct%iatype(iat)))=='Tm')then
        rcov(iat)=2.92d0
     else if (trim(atoms%astruct%atomnames(atoms%astruct%iatype(iat)))=='Yb')then
        rcov(iat)=2.80d0
     else if (trim(atoms%astruct%atomnames(atoms%astruct%iatype(iat)))=='Lu')then
        rcov(iat)=2.80d0
     else if (trim(atoms%astruct%atomnames(atoms%astruct%iatype(iat)))=='Hf')then
        rcov(iat)=2.90d0
     else if (trim(atoms%astruct%atomnames(atoms%astruct%iatype(iat)))=='Ta')then
        rcov(iat)=2.70d0
     else if (trim(atoms%astruct%atomnames(atoms%astruct%iatype(iat)))=='W' )then
        rcov(iat)=2.60d0
     else if (trim(atoms%astruct%atomnames(atoms%astruct%iatype(iat)))=='Re')then
        rcov(iat)=2.60d0
     else if (trim(atoms%astruct%atomnames(atoms%astruct%iatype(iat)))=='Os')then
        rcov(iat)=2.50d0
     else if (trim(atoms%astruct%atomnames(atoms%astruct%iatype(iat)))=='Ir')then
        rcov(iat)=2.50d0
     else if (trim(atoms%astruct%atomnames(atoms%astruct%iatype(iat)))=='Pt')then
        rcov(iat)=2.60d0
     else if (trim(atoms%astruct%atomnames(atoms%astruct%iatype(iat)))=='Au')then
        rcov(iat)=2.70d0
     else if (trim(atoms%astruct%atomnames(atoms%astruct%iatype(iat)))=='Hg')then
        rcov(iat)=2.80d0
     else if (trim(atoms%astruct%atomnames(atoms%astruct%iatype(iat)))=='Tl')then
        rcov(iat)=2.50d0
     else if (trim(atoms%astruct%atomnames(atoms%astruct%iatype(iat)))=='Pb')then
        rcov(iat)=3.30d0
     else if (trim(atoms%astruct%atomnames(atoms%astruct%iatype(iat)))=='Bi')then
        rcov(iat)=2.90d0
     else if (trim(atoms%astruct%atomnames(atoms%astruct%iatype(iat)))=='Po')then
        rcov(iat)=2.80d0
     else if (trim(atoms%astruct%atomnames(atoms%astruct%iatype(iat)))=='At')then
        rcov(iat)=2.60d0
     else if (trim(atoms%astruct%atomnames(atoms%astruct%iatype(iat)))=='Rn')then
        rcov(iat)=2.60d0
     else
        call yaml_comment('(SQNM) no covalent radius stored for this atomtype '&
             //trim(atoms%astruct%atomnames(atoms%astruct%iatype(iat))))
        stop
     endif
     if (iproc == 0) then
        call yaml_map('(SQNM) RCOV:'//trim(atoms%astruct%atomnames&
                 (atoms%astruct%iatype(iat))),rcov(iat))
     endif
  enddo
end subroutine give_rcov_sqnm
