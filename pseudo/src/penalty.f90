!> @file
!! penalty routine for amoeba
!! @author
!!    Alex Willand, under the supervision of Stefan Goedecker
!!    gpu accelerated routines by Raffael Widmer
!!    parts of this program were based on the fitting program by matthias krack
!!    http://cvs.berlios.de/cgi-bin/viewcvs.cgi/cp2k/potentials/goedecker/pseudo/v2.2/


!> penalty routine for amoeba
subroutine penalty(energ,verbose,pp,penal,  &
        iproc,nproc,iter,penref)
   
   use pseudovars
   use penaltyvars
   use gatomvars

   implicit none

   !Arguments
   logical, intent(in) :: energ, verbose
   integer, intent(in) :: iproc,nproc, iter
   real(kind=8), dimension(maxdim) :: pp
   real(kind=8), dimension(7) :: pen_cont
   real(kind=8), dimension(4*nproc) :: exverbose
   real(kind=8), intent(out) :: penal
   real(kind=8), intent(in) :: penref
   !Local variables
   integer :: iorb,nocc,l,ispin,i,ierr,k
   real(kind=8) :: dd,ekin_pen,eref,excit
   real(kind=8) :: pen_h, pen_k, pen_loc,pen_r,penalorb,r,t,tk
   
   include 'mpif.h'
   
   ! set res() =-1 for orbitals with zero charge and
   !  wght(nocc,l+1,ispin,5) set to zero
   ! this avoids unneccessay computations of res() in the
   ! routine resid()
   
   ! when we print out detailed info with verbose,
   ! we actually do not want to miss any residues
   
   do iorb=1,norb
      nocc=no(iorb)
      l=lo(iorb)
      ispin=1
      if (so(iorb).lt.0) ispin=2
      if ( (wght(nocc,l+1,ispin,5).eq.0.0d0) .and.  &
           (occup(nocc,l+1,ispin).lt.1.0d-8) .and.  &
           ( .not. verbose) )  then
         res(nocc,l+1,ispin) = -1.d0
      else
         res(nocc,l+1,ispin) = 0.d0
      endif
   enddo
   ! unpack variables
   ! print*,'penalty: maxdim=',maxdim
   ! convention for spin treatment with libxc using the pol variable
   call ppack (verbose, pp(1), 'unpack')
   
   call cpu_time(t)
   time(1)=time(1)-t
   call gatom(energ,verbose)
   call cpu_time(t)
   time(1)=time(1)+t
   penal=0d0
   !     diff for dcharg and echarge is calc. in (%)
   do iorb=1,norb
      nocc=no(iorb)
      l=lo(iorb)
      ispin=1
      if (so(iorb).lt.0) ispin=2
      
      penalorb= &
           ((aeval(nocc,l+1,ispin)-ev(iorb))  &
           *wght(nocc,l+1,ispin,1))**2 +  &
           ((chrg(nocc,l+1,ispin)-crcov(iorb))  &
           *wght(nocc,l+1,ispin,2))**2 +  &
           (100.d0*(1.d0-dhrg(nocc,l+1,ispin)/dcrcov(iorb))  &
           *wght(nocc,l+1,ispin,3))**2 +  &
           (100.d0*(1.d0-ehrg(nocc,l+1,ispin)/ddcrcov(iorb))  &
           *wght(nocc,l+1,ispin,4))**2 +  & 
           (res(nocc,l+1,ispin)*wght(nocc,l+1,ispin,5))**2 +  &
           (wfnode(nocc,l+1,ispin,1)*wght(nocc,l+1,ispin,6))**2 +  &
           (wfnode(nocc,l+1,ispin,2)*wght(nocc,l+1,ispin,7))**2 +  &
           (wfnode(nocc,l+1,ispin,3)*wght(nocc,l+1,ispin,8))**2  
      
      if (penalorb**2>= 0d0) then
         !                some isnan test... 
         !                we dont want to add nan to the penalty,
         !                but rather give info for debugging
         penal=penal+penalorb
      else
         write(6,*) 'warning, nan penalty for orbital',iorb, penalorb
         
         write(6,*) 'debug: nocc,l,ispin',nocc,l,ispin
         write(6,*) 'ae/ps evals',aeval(nocc,l+1,ispin),ev(iorb)
         write(6,*) 'node integrals ', wfnode(nocc,l+1,ispin,:)
         write(6,*) 'all wghts',wght(nocc,l+1,ispin,1:8)
         write(6,*)
      end if
   enddo
   !     for now, keep this contribution here, no special output
   penal=penal + (psir0*wghtp0)**2
   
   !     add weight for narrow radii
   !     use: error in norm conversation of projectors
   !          is approximately a power law of the radius
   !
   !         1-|projector|  ~ 1d-6 (hgrid/r)**12
   
   !     write(16,*)'debug: sqpenal without rad',penal
   pen_r=0d0
   do l=1,lpx
      pen_r=pen_r+(hgridmax/r_l(l))**24
   end do
   !     just use the same convention for rloc
   pen_r=pen_r+(hgridmax/rloc)**24
   !     prefactors
   pen_r=pen_r*1d-12*wghtrad**2
   !     write(16,*)'debug: sqpenal with rad',penal
   
   !     add an empirical term to penalize all hsep larger than ~10
   !     also calculate kij and penalize values larger than ~ 2,
   !     otherwise we might favor big splittings.
   !     use a power law such that h~10 and k~2 contribute each about
   !     1% of wghthij to the penalty, but larger values are strictly
   !     penalized.
   pen_h=0d0
   pen_k=0d0
   do l=1,lpx
      do i=1,6
         if (hsep(i,l,1).eq.0)cycle
         pen_h= pen_h+(hsep(i,l,1)*0.1d0)**12
         if (hsep(i,l,2).eq.0)cycle
         pen_h= pen_h+(hsep(i,l,2)*0.1d0)**12
         tk=2d0*(hsep(i,l,1)-hsep(i,l,2))/(2*l-1)
         if (tk.eq.0)cycle
         !          we may want to try different weights on kij later
         pen_k= pen_k+(tk*0.5d0)**12
      end do
   end do
   !     the prefactor of 100 for the weight can be dismissed
   pen_h=pen_h*wghthij**2*1d-4
   pen_k=pen_k*wghthij**2*1d-4
   
   !     let us add one more empirical term to keep the local part "local"
   pen_loc=0d0
   do k=1,nint
      r=rr(k)
      if (r<rcov) cycle
      if (r>1.5*rcov) exit
      !        local potential vloc(r)
      dd=exp(-.5d0*(r/rloc)**2)*  &
           (gpot(1) + gpot(2)*(r/rloc)**2+    &
           gpot(3)*(r/rloc)**4 +  &
           gpot(4)*(r/rloc)**6 )
      pen_loc=pen_loc+dd**2
   end do
   pen_loc=pen_loc*wghtloc**2
   
   !     if we're not in the middle of a fit,
   !     give some more information about these penalty terms
   if (verbose) then
      write(6,*)
      write(6,*) 'empirical penalty terms from psppar'
      write(6,*) '___________________________________'
      write(6,*) 
      write(6,*) 'narrow radii   ',sqrt(pen_r)
      write(6,*) 'large  hij     ',sqrt(pen_h)
      write(6,*) 'large  kij     ',sqrt(pen_k)
      write(6,*) 'vloc(r>rcov)   ',sqrt(pen_loc)
      write(6,*)
      
   end if
   
   pen_cont(7)=pen_r+ pen_h+pen_k + pen_loc
   
   
   !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   !---------------- parallel computation of excitation energies and softness ----------------
   !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   
   !     first, make sure that excitation energies of value zero are ignored
   if (excitae.eq.0d0)wghtexci=0d0
   
   !     if the weight for softness is nonzero,
   !     calculate the change in kinetic energy when
   !     transforming psi into a debauchie wavelet basis
   
   if (wghtsoft .gt. 0d0 .and. nhgrid .gt. 0 ) then ! .and. mod(iter,nskip).eq.0)
      call ekin_wvlt(verbose,iproc,nproc,ng,ngmx,  &
           noccmx,lmx,nspin,nsmx,  &
           nhgrid, hgridmin,hgridmax, nhpow,ampl,crmult,  &
           frmult, xp,psi,occup,ekin_pen,time)
      !        important: here we scale the radii for the coarse
      !                   and fine grid with the parameter rcov.
   else
      ekin_pen=0d0
   end if
   
   if (nproc.eq.1) then
      !        serial case: one configuration, no excitation energies
      pen_cont(4)=wghtconf**2*penal
      pen_cont(6)=wghtsoft**2*ekin_pen
      penal=sqrt(penal+sum(pen_cont(6:7)))
      !        write(16,*)'debug:ekin_pen,penal',ekin_pen,penal,wghtsoft
   else 
      !        parallel case
      !        set up penalty contributions to be summed over all processes
      call cpu_time(t)
      time(2)=time(2)-t
      
      !        compute excitation energy for this configuration. 
      !        for this, get total energy of configuration 1.
      if (iproc.eq.0) eref=etotal
      call mpi_bcast(eref,1,mpi_double_precision,0,  &
           mpi_comm_world,ierr)
      if (ierr.ne.0) write(*,*)'mpi_bcast ierr',ierr
      excit=etotal-eref
      
      !        sum up penalty contributions from this mpi process
      pen_cont(1)=wghtconf**2 *penal
      pen_cont(2)=wghtexci**2 *(excit-excitae)**2
      pen_cont(3)=wghtsoft**2 *ekin_pen
      
      !        write(6,*)'debug lines for penalty over processes'
      !        write(6,*)'wghtconf,penal',wghtconf ,penal
      !        write(6,*)'wghtexci,(excit-excitae)**2',&
           !    &       wghtexci,(excit-excitae)**2
      !        write(6,*)'wghtsoft,ekin_pen',wghtsoft,ekin_pen
      !        write(6,*)'products',pen_cont(1:3)
      
      !        sum up over all processes, that is over all configurations
      !        and over the orbitals and hgrid samples of the wavelet transform
      call mpi_allreduce(pen_cont(1:3),pen_cont(4:6),3,  &
           mpi_double_precision,mpi_sum,mpi_comm_world,ierr)
      if (ierr.ne.0) write(*,*)'mpi_allreduce ierr',ierr
      call cpu_time(t)
      time(2)=time(2)+t
      penal=sqrt(sum(pen_cont(4:7)))
      !        write(6,*)'mpi reduced:',pen_cont(4:6)
      !        write(6,*)'pen7=', pen_cont(7),'; sqrtsum=penal',penal
      
   end if
   
   !     if the penalty is below the reference - usually the currently
   !     lowest vertex of the simplex -  write out the major components:
   
   if (penal<penref) then
      if (nproc.eq.1) then
         write(6,'(2i8,4(1x,1pe20.10))')  iter, ntime, penal, sqrt(pen_cont(6)), sqrt(pen_cont(7)),&
              sqrt(pen_cont(4))
      else
         pen_cont(4:7)=sqrt(pen_cont(4:7))
         write(6, '(2i8,7(1x,1pe20.10))') iter, ntime, penal, pen_cont(6),pen_cont(7), pen_cont(4:5), &
              sqrt(sum(pen_cont(1:2)))
      end if
      !         write the vertex to the dumpfile without transforming  psppar
      !         this gives the psppar of all vertices ever found to be the lowest
      
      !         disabling this should speed up things a little
      
      !         backspace(99)
      if (iproc.eq.0) then
         write(99,'(1pe11.3,a)') penal,' penalty' 
         write(99,'(5(1pe11.3),t65,a)') rloc,gpot,  ' rloc, gpot'
         write(99,'(5(1pe11.3),t65,a)') rcore,gcore,'rcore, gcore'
         do l=1,lpx
            write(99,'(f7.3,t8,6e11.3,t76,a)') r_l(l), (hsep(i,l,1),i=1,6),'r_l,hij(up)'
            if (l.gt.1-nspol .and. nspin.eq.2)  &
                 write(99,'(t8,6e11.3,t76,a)') (hsep(i,l,2),i=1,6),'      hij(dn)'
         enddo
         write(99,*)  
      end if
   end if
   
   
   if (nproc>1.and.verbose) then
      !         write out the excitation energies
      !         get excitation energies from all processes  
      !         using mpiallreduce may be a clumsy way of doing this
      exverbose=0d0
      exverbose(2*iproc+1)=excit
      exverbose(2*iproc+2)=sqrt(pen_cont(2))
      !         note: write out energies in the convention they have been read
      call mpi_allreduce(exverbose(1),  &
           exverbose(2*nproc+1),2*nproc,  &
           mpi_double_precision,mpi_sum,mpi_comm_world,ierr)
      write(6,'(/,1x,a)') 'excitation energies'
      write(6,'(1x,a,/)') '___________________'
      write(6,'(1x,a)')    'configuration,    de=etot-etot1,    (de-de_ae)*weight'
      do i=1,nproc
         write(6,'(10x,i4,3e20.12)')  &
              i-1,exverbose(2*nproc+2*i-1),  &
              exverbose(2*nproc+2*i)
      end do
      write(6,*)
   end if
   ierr=0
   if (nproc>1) call mpi_barrier(mpi_comm_world,ierr)  
   if (ierr.ne.0) write(*,*)'mpi_barrier ierr',ierr
   
end subroutine penalty
