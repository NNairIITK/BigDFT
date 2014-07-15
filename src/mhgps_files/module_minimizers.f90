module module_minimizers
implicit none

contains

!> @file
!!  Routines for Stefan's new minimization method
!! @author
!!    Copyright (C) 2014 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS
!subroutine geopt(nat,wpos,etot,fout,fnrmtol,count,count_sd,displr)
subroutine minimizer_sbfgs(imode,nat,alat,atomnames,nbond,iconnect,rxyzio,fxyzio,fnoiseio,energyio,energycounter,converged)
!call_bigdft has to be run once on runObj and outs !before calling this routine
!sbfgs will return to caller the energies and coordinates used/obtained from the last accepted iteration step
   use module_base
   use module_types
   use module_interfaces
   use yaml_output
   use module_sbfgs
   use module_energyandforces
   use module_global_variables, only: iproc,&
                                      inputPsiId,&
                                      mhgps_verbosity,&
                                      ixyz_int,&
                                      atoms,&
                                      mini_frac_fluct,&
                                      mini_ncluster_x,&
                                      mini_betax,&
                                      mini_beta_stretchx,&
                                      mini_nhistx,&
                                      mini_maxrise,&
                                      mini_cutoffratio,&
                                      mini_steepthresh,&
                                      mini_trustr,&
                                      mini_forcemax,&
                                      fdim
   implicit none
   !parameter
   integer, intent(in)                    :: nat, nbond,imode
   real(gp), intent(inout)                :: energycounter
   logical, intent(out)                   :: converged
   real(gp), intent(inout)                :: rxyzio(3,nat),fxyzio(3,nat),alat(3,nat)
   real(gp), intent(inout)                :: energyio,fnoiseio
   character(20), intent(in)              :: atomnames(nat)
   integer, intent(in)                    :: iconnect(2,nbond)
   !local variables
   character(len=*), parameter :: subname='sbfgs'
   integer :: infocode,info !< variables containing state codes
   integer :: nhistx !< maximum history length
   integer :: nhist  !< actual history length
   integer :: ndim   !< dimension of significant subspace
   integer :: nit    !< maximum number of iterations
   integer :: istat,iall
   integer :: lwork
   integer :: it,i,iat,l,j,idim,jdim,ihist,icheck !<counter variables
   integer :: itswitch
   logical :: debug !< set .true. for debug output to fort.100
   logical :: steep !< steepest descent flag
   logical :: success
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
   real(gp) :: beta_stretch  !< current step size
   real(gp) :: cosangle
   real(gp) :: ts
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
   real(gp), allocatable, dimension(:,:,:) :: rxyz
   real(gp), allocatable, dimension(:,:,:) :: rxyzraw
   real(gp), allocatable, dimension(:,:,:) :: fxyz
   real(gp), allocatable, dimension(:,:,:) :: fxyzraw
   real(gp), allocatable, dimension(:,:,:) :: fstretch
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
   character(len=4)                        :: fn4
   character(len=40)                       :: comment
   character(len=9)                        :: cdmy9_1
   character(len=9)                        :: cdmy9_2
   character(len=9)                        :: cdmy9_3
   character(len=8)                        :: cdmy8


   !set parameters
!   nit=runObj%inputs%ncount_cluster_x
!   nat=runObj%atoms%astruct%nat
!   betax=runObj%inputs%betax
!   nhistx=runObj%inputs%nhistx
!   maxrise=runObj%inputs%maxrise
!   cutoffRatio=runObj%inputs%cutoffratio
!   steepthresh=runObj%inputs%steepthresh
!   trustr=runObj%inputs%trustr
   nit           =mini_ncluster_x
   betax         =mini_betax
   beta_stretchx =mini_beta_stretchx
   nhistx        =mini_nhistx
   maxrise       =mini_maxrise
   cutoffRatio   =mini_cutoffratio
   steepthresh   =mini_steepthresh
   trustr        =mini_trustr

   if (iproc==0.and.mhgps_verbosity > 0) then
      call yaml_mapping_open('Geometry parameters')
         call yaml_map('Geometry Method','GEOPT_SBFGS')
         call yaml_map('nhistx',nhistx)
         call yaml_map('betax', betax,fmt='(1pe21.14)')
         call yaml_map('maxrise', maxrise,fmt='(1pe21.14)')
         call yaml_map('cutoffRatio', cutoffRatio,fmt='(1pe21.14)')
         call yaml_map('steepthresh', steepthresh,fmt='(1pe21.14)')
         call yaml_map('trustr', trustr,fmt='(1pe21.14)')
      call yaml_mapping_close()
   end if

   !init varaibles
   debug=.false.
   converged=.false.
   displr=0.0_gp
   displp=0.0_gp
   fluct=0.0_gp
   icheck=0
   detot=0.0_gp
   itswitch=0
   ndim=0
   nhist=0
   beta=betax
   beta_stretch=beta_stretchx
   maxd=0.0_gp

   ! allocate arrays
   lwork=1000+10*nat**2
   rxyz = f_malloc((/ 1.to.3, 1.to.nat, 0.to.nhistx /),id='rxyz')
   rxyzraw = f_malloc((/ 1.to.3, 1.to.nat, 0.to.nhistx /),id='rxyzraw')
   fxyz = f_malloc((/ 1.to.3, 1.to.nat, 0.to.nhistx /),id='fxyz')
   fxyzraw = f_malloc((/ 1.to.3, 1.to.nat, 0.to.nhistx /),id='fxyzraw')
   fstretch = f_malloc((/ 1.to.3, 1.to.nat, 0.to.nhistx /),id='fstretch')
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
   wold = f_malloc((/ 1.to.nbond/),id='wold')
   



!!!!!!   call energyandforces(nat,rxyz(1,1,0),fxyz(1,1,0),etot)
!!  not necessary, call_bigdft allready called outside
!   call call_bigdft(runObj,outs,nproc,iproc,infocode)
!   energycounter=energycounter+1

!! copy to internal variables
   call vcopy(3*nat, rxyzio(1,1), 1,rxyz(1,1,0), 1)
   call vcopy(3*fdim, fxyzio(1,1), 1, fxyz(1,1,0), 1)
   etot=energyio
   fnoise=fnoiseio
   call minenergyandforces(.false.,imode,nat,alat,rxyz(1,1,0),&
       rxyzraw(1,1,0),fxyz(1,1,0),fstretch(1,1,0),fxyzraw(1,1,0),&
       etot,iconnect,nbond,atomnames,wold,beta_stretchx,beta_stretch)
   if(imode==2)rxyz(:,:,0)=rxyz(:,:,0)+beta_stretch*fstretch(:,:,0)

   call fnrmandforcemax(fxyzraw(1,1,0),fnrm,fmax,nat)
   fnrm=sqrt(fnrm)
   if (fmax < 3.e-1_gp) call updatefluctsum(fnoise,fluct)

   etotold=etot
   etotp=etot

   if (iproc==0.and.mhgps_verbosity > 0) then
       !avoid space for leading sign (numbers are positive, anyway)
       write(cdmy8,'(es8.1)')abs(maxd)
       write(cdmy9_1,'(es9.2)')abs(displr)
       write(cdmy9_2,'(es9.2)')abs(displp)
       write(cdmy9_3,'(es9.2)')abs(beta)

       write(16,'(i5,1x,i5,2x,a10,2x,1es21.14,2x,es9.2,es11.3,3es10.2,2x,a6,a8,1x,a4,i3.3,1x,a5,a7,2(1x,a6,a8))') &
       int(energycounter),0,'GEOPT_SBFGS',etotp,detot,fmax,fnrm,fluct*mini_frac_fluct,fluct, &
       'beta=',trim(adjustl(cdmy9_3)),'dim=',ndim,'maxd=',trim(adjustl(cdmy8)),'dsplr=',trim(adjustl(cdmy9_1)),'dsplp=',trim(adjustl(cdmy9_2))
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

      ! make space in the history list
      if (nhist.gt.nhistx) then
         nhist=nhistx
         do ihist=0,nhist-1
            do iat=1,nat
               do l=1,3
                  rxyz(l,iat,ihist)=rxyz(l,iat,ihist+1)
                  rxyzraw(l,iat,ihist)=rxyzraw(l,iat,ihist+1)
                  fxyz(l,iat,ihist)=fxyz(l,iat,ihist+1)
                  fxyzraw(l,iat,ihist)=fxyzraw(l,iat,ihist+1)
                  fstretch(l,iat,ihist)=fstretch(l,iat,ihist+1)
               enddo
            enddo
         enddo
      endif
   
      ! decompose gradient
500 continue
    call modify_gradient(nat,ndim,rrr(1,1,1),eval(1),res(1),fxyz(1,1,nhist-1),beta,dd(1,1))
   
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
         if(iproc==0)write(16,'(a,2(1x,es9.2))')'WARNING GEOPT_SBFGS: step too large: maxd, trustradius ',maxd,trustr
         scl=0.50_gp*trustr/maxd
         dd=dd*scl
         tt=tt*scl
         maxd=maxd*scl
      endif
      displr=displr+tt
   
      !update positions
      do iat=1,nat
         rxyz(1,iat,nhist)=rxyz(1,iat,nhist-1)-dd(1,iat)
         rxyz(2,iat,nhist)=rxyz(2,iat,nhist-1)-dd(2,iat)
         rxyz(3,iat,nhist)=rxyz(3,iat,nhist-1)-dd(3,iat)
      enddo
   
!      call energyandforces(nat,rxyz(1,1,nhist),fxyz(1,1,nhist),etotp)
!      call vcopy(3 * runObj%atoms%astruct%nat, rxyz(1,1,nhist), 1,runObj%atoms%astruct%rxyz(1,1), 1)
!      runObj%inputs%inputPsiId=1
!      call call_bigdft(runObj,outs,nproc,iproc,infocode)
!      energycounter=energycounter+1
!      call vcopy(3 * outs%fdim, outs%fxyz(1,1), 1, fxyz(1,1,nhist), 1)
!      etotp=outs%energy
!      detot=etotp-etotold

      inputPsiId=1
!      call energyandforces(nat,alat,rxyz(1,1,nhist),fxyz(1,1,nhist),fnoise,etotp)
      call minenergyandforces(.true.,imode,nat,alat,rxyz(1,1,nhist),rxyzraw(1,1,nhist),&
                             fxyz(1,1,nhist),fstretch(1,1,nhist),fxyzraw(1,1,nhist),&
                             etotp,iconnect,nbond,atomnames,wold,beta_stretchx,beta_stretch)
      detot=etotp-etotold
      energycounter=energycounter+1.0_gp

      if(debug.and.iproc==0)write(100,'(a,i6,2(1x,e21.14),1x,5(1x,e10.3),1x,i0)')&
            'SBFGS it,etot,etotold,Detot,fnrm,fnrmp/fnrm,dnrm/fnrm,beta,ndim',&
             it-1,etotp,etotold,detot,fnrm,sqrt(ts)/fnrm,sqrt(tt)/fnrm,beta,ndim


      call fnrmandforcemax(fxyzraw(1,1,nhist),fnrm,fmax,nat)
      fnrm=sqrt(fnrm)

      if (iproc == 0 .and. mhgps_verbosity >=4) then
         write(fn4,'(i4.4)') int(energycounter)
         write(comment,'(a,1pe10.3)')'SBFGS:fnrm= ',fnrm
         call write_atomic_file('posmini_'//fn4, &
              etotp,rxyz(1,1,nhist),ixyz_int,&
              atoms,trim(comment),forces=fxyz(1,1,nhist))
      endif

      if (fmax < 3.e-1_gp) call updatefluctsum(fnoise,fluct)
      cosangle=-dot_double(3*nat,fxyz(1,1,nhist),1,dd(1,1),1)/&
              sqrt(dot_double(3*nat,fxyz(1,1,nhist),1,fxyz(1,1,nhist),1)*&
              dot_double(3*nat,dd(1,1),1,dd(1,1),1))

      if (detot.gt.maxrise .and. beta > 1.e-1_gp*betax) then !
         if (debug.and.iproc==0) write(100,'(a,i0,1x,e9.2)') "WARN: it,detot", it,detot
         if (debug.and.iproc==0) write(16,'(a,i0,4(xe9.2))') &
             "WARNING GEOPT_SBFGS: Prevent energy to rise by more than maxrise: it,maxrise,detot,beta,1.e-1*betax ",&
             it,maxrise,detot,beta,1.e-1_gp*betax
         if (iproc==0.and.mhgps_verbosity > 0) then
            !avoid space for leading sign (numbers are positive, anyway)
            write(cdmy8,'(es8.1)')abs(maxd)
            write(cdmy9_1,'(es9.2)')abs(displr)
            write(cdmy9_2,'(es9.2)')abs(displp)
            write(cdmy9_3,'(es9.2)')abs(beta)
   
            write(16,'(i5,1x,i5,2x,a10,2x,1es21.14,2x,es9.2,es11.3,3es10.2,2x,a6,a8,1x,a4,i3.3,1x,a5,a7,2(1x,a6,a8))') &
             int(energycounter),it,'GEOPT_SBFGS',etotp,detot,fmax,fnrm,fluct*mini_frac_fluct,fluct, &
             'beta=',trim(adjustl(cdmy9_3)),'dim=',ndim,'maxd=',trim(adjustl(cdmy8)),'dsplr=',trim(adjustl(cdmy9_1)),'dsplp=',trim(adjustl(cdmy9_2))
            call yaml_mapping_open('Geometry')
               call yaml_map('Ncount_BigDFT',int(energycounter))
               call yaml_map('Geometry step',it)
               call yaml_map('Geometry Method','GEOPT_SBFGS')
               call yaml_map('ndim',ndim)
               call yaml_map('etot', etotp,fmt='(1pe21.14)')
               call yaml_map('detot',detot,fmt='(1pe21.14)')
               call yaml_map('fmax',fmax,fmt='(1pe21.14)')
               call yaml_map('fnrm',fnrm,fmt='(1pe21.14)')
               call yaml_map('beta',beta,fmt='(1pe21.14)')
               call geometry_output(fmax,fnrm,fluct)
            call yaml_mapping_close()
         end if
    
         if(int(energycounter) >= nit)then!no convergence within ncount_cluster_x energy evaluations
            !following copy of rxyz(1,1,nhist-1) to runObj is necessary for returning to the caller
            !the energies and coordinates used/obtained from/in the last ACCEPTED iteration step
            !(otherwise coordinates of last call to call_bigdft would be returned)
!            call vcopy(3 * nat, rxyz(1,1,nhist-1), 1,runObj%atoms%astruct%rxyz(1,1), 1)
            energyio=etotold
            do iat=1,nat
              do l=1,3
                 rxyzio(l,iat)= rxyz(l,iat,nhist-1)
                 fxyzio(l,iat)= fxyzraw(l,iat,nhist-1)
              enddo
           enddo
           goto 900  !sbfgs will return to caller the energies and coordinates used/obtained from the last ACCEPTED iteration step
         endif

         !beta=min(.50_gp*beta,betax)
         beta=.50_gp*beta
         if (debug.and.iproc==0) write(100,'(a,1x,e9.2)') 'WARNING GEOPT_SBFGS: beta reset ',beta
         ndim=0
         if(.not.steep)then
            do iat=1,nat
               rxyz(1,iat,0)=rxyz(1,iat,nhist-1)
               rxyz(2,iat,0)=rxyz(2,iat,nhist-1)
               rxyz(3,iat,0)=rxyz(3,iat,nhist-1)
               rxyzraw(1,iat,0)=rxyzraw(1,iat,nhist-1)
               rxyzraw(2,iat,0)=rxyzraw(2,iat,nhist-1)
               rxyzraw(3,iat,0)=rxyzraw(3,iat,nhist-1)
   
               fxyz(1,iat,0)=fxyz(1,iat,nhist-1)
               fxyz(2,iat,0)=fxyz(2,iat,nhist-1)
               fxyz(3,iat,0)=fxyz(3,iat,nhist-1)
               fxyzraw(1,iat,0)=fxyzraw(1,iat,nhist-1)
               fxyzraw(2,iat,0)=fxyzraw(2,iat,nhist-1)
               fxyzraw(3,iat,0)=fxyzraw(3,iat,nhist-1)
            enddo
            nhist=1
         endif
         goto  500
      endif

      displp=displp+tt
      if (iproc==0.and.mhgps_verbosity > 0) then
         !avoid space for leading sign (numbers are positive, anyway)
         write(cdmy8,'(es8.1)')abs(maxd)
         write(cdmy9_1,'(es9.2)')abs(displr)
         write(cdmy9_2,'(es9.2)')abs(displp)
         write(cdmy9_3,'(es9.2)')abs(beta)

         write(16,'(i5,1x,i5,2x,a10,2x,1es21.14,2x,es9.2,es11.3,3es10.2,2x,a6,a8,1x,a4,i3.3,1x,a5,a7,2(1x,a6,a8))') &
          int(energycounter),it,'GEOPT_SBFGS',etotp,detot,fmax,fnrm,fluct*mini_frac_fluct,fluct, &
          'beta=',trim(adjustl(cdmy9_3)),'dim=',ndim,'maxd=',trim(adjustl(cdmy8)),'dsplr=',trim(adjustl(cdmy9_1)),'dsplp=',trim(adjustl(cdmy9_2))
         call yaml_mapping_open('Geometry')
            call yaml_map('Ncount_BigDFT',int(energycounter))
            call yaml_map('Geometry step',it)
            call yaml_map('Geometry Method','GEOPT_SBFGS')
            call yaml_map('ndim',ndim)
            call yaml_map('etot', etotp,fmt='(1pe21.14)')
            call yaml_map('detot',detot,fmt='(1pe21.14)')
            call yaml_map('fmax',fmax,fmt='(1pe21.14)')
            call yaml_map('fnrm',fnrm,fmt='(1pe21.14)')
            call yaml_map('beta',beta,fmt='(1pe21.14)')
            call geometry_output(fmax,fnrm,fluct)
         call yaml_mapping_close()
      end if

      etot    = etotp
      etotold = etot

      if(detot .gt. maxrise)then
         if (iproc==0) write(16,'(a,i0,4(xe9.2))') &
             "WARNING GEOPT_SBFGS: Allowed energy to rise by more than maxrise: it,maxrise,detot,beta,1.e-1*betax ",&
             it,maxrise,detot,beta,1.e-1_gp*betax
      endif


!      if (fnrm.le.fnrmtol) goto 1000
      call convcheck(fmax,fluct*mini_frac_fluct,mini_forcemax,icheck)
      if(icheck>5)then
         goto 1000
      endif
     if(imode==2)rxyz(:,:,nhist)=rxyz(:,:,nhist)+beta_stretch*fstretch(:,:,nhist) !has to be after convergence check,
                                                                       !otherwise energy will not match
                                                                       !the true energy of rxyz(:,:,nhist)

      if(int(energycounter) >= nit)then!no convergence within ncount_cluster_x energy evaluations
            energyio=etot
            do iat=1,nat
              do l=1,3
                 rxyzio(l,iat)= rxyz(l,iat,nhist)
                 fxyzio(l,iat)= fxyzraw(l,iat,nhist)
              enddo
           enddo
            goto 900  !sbfgs will return to caller the energies and coordinates used/obtained from the last accepted iteration step
      endif
   
      if (cosangle.gt..200_gp) then
         beta=beta*1.100_gp
      else
         beta=max(beta*.850_gp,betax)
      endif
   
      if (debug.and.iproc==0) write(100,*) 'cosangle ',cosangle,beta

      call getSubSpaceEvecEval(nat,nhist,nhistx,ndim,cutoffratio,lwork,work,rxyz,&
                   &fxyz,aa,rr,ff,rrr,fff,eval,res,success)
      if(.not.success)stop 'subroutine minimizer_sbfgs: no success in getSubSpaceEvecEval.'

   enddo!end main loop

900 continue

   !if code gets here, it failed
   if(debug.and.iproc==0) write(100,*) it,etot,fnrm
   if(iproc==0) write(16,'(a,3(1x,i0))') &
       "WARNING GEOPT_SBFGS: SBFGS not converged: it,energycounter,ncount_cluster_x: ", &
       it,int(energycounter),mini_ncluster_x
!   stop "No convergence "
   converged=.false.
   goto 2000

1000 continue!converged successfully
   
   if(iproc==0) write(16,'(2(a,1x,i0))') "SBFGS converged at iteration ",it,". Needed energy calls: ",int(energycounter)
   if(iproc==0)  call yaml_map('Iterations when SBFGS converged',it)
   converged=.true.
   
   energyio=etotp
   do iat=1,nat
      do l=1,3
         rxyzio(l,iat)= rxyz(l,iat,nhist)
         fxyzio(l,iat)= fxyzraw(l,iat,nhist)
      enddo
   enddo
2000 continue
!deallocations
   call f_free(rxyz)
   call f_free(fxyz)
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
end subroutine

end module
