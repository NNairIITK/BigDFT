!> @file
!! File for the pseudopotential generator
!! @author
!!    Copyright (C) 2010-2013 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS


!> pack the fittingparameters into one array for amoeba
!! and vice versa
!! if nstring = 'init'   -> initialisation + packing
!! if nstring = 'pack'   -> packing
!! if nstring = 'unpack' -> unpacking

!! modifications:
!!               fitting of frozen core charge,
!!               spin polarization
!!               moved save block of old psppar to a module

!! let us define a new integer nso for spin orbit,
!! the psp spin components. if nso=1, all kij=0.

!! it widely replaces npsin here, which is the number of
!! spin channels for the wavefunctions.

!! the number of spin channels for charge and xc is
!! nspol, which equals 2 for polarized dft calculations
!! and 1 otherwise. thus, nspin = max(nspol,nso).

!! for now, relativistic polarized dft with
!! nspol=nso=nspin=2 is not supported.

!! the list format reading of booleans from input.fitpar
!! in fixed order is discarded in this version.
!! the user rather gives a list of keywords in input.fitpar
!! that name the psppar that are to be fitted.
!! if input.fitpar is empty and the verbose flag is on,
!! some help text will be printed.

!! @note
!!   many lines of this routine are here only
!!   to perform rather meaningless transformations
!!   of the hij and kij for given l. these can
!!   be eliminated, but are still here for testing.
subroutine ppack (verbose, pp, nstring)
   
   use pseudovars
   use ppackvars
   implicit none
   
   logical verbose
   real(8):: pp(maxdim) ! this time, it is just one vertex!
   character(*):: nstring
   
   !local variables
   integer increm,lpxfit,ierr
   real(8):: h11,h12,h13,h22,h23,h33,hh11,hh12,hh13,hh22,hh23,hh33
   
   integer maxpar,iline
   parameter (maxpar = 42) !up to 5 local,  4*7 nonlocal, 5 nlcc, 4 r_l2
   character*10 spack(maxpar)
   character*100 string
   logical lpack(maxpar),match
   integer i,j,ncount,nhsep1,nhsep2,nhsep3,ll,nn,nnmin,nso
   
   real(8):: pih,xshift
   
   !> The spack array gives the names of all possible input.fitpar
   !! the meaning of the index is the same as for lpack
   data spack /'rloc', 'gpot1', 'gpot2', 'gpot3', 'gpot4',&
        'rcore','gcore1','gcore2','gcore3','gcore4',&
        'rs','hs11','hs12','hs22','hs13','hs23','hs33','rs2',&
        'rp','hp11','hp12','hp22','hp13','hp23','hp33','rp2',&
        'rd','hd11','hd12','hd22','hd13','hd23','hd33','rd2',&
        'rf','hf11','hf12','hf22','hf13','hf23','hf33','rf2'/ 
   
   !> @note
   !!   hsep index convention:
   !!   hsep 1   -    h11
   !!   hsep 2   -    h12
   !!   hsep 3   -    h22
   !!   hsep 4   -    h13
   !!   hsep 5   -    h23
   !!   hsep 6   -    h33
   
   nso=1
   if(nspin>nspol)nso=2
   
   !
   !------------------------------------------------------------------
   !     initialisation: set up the list of parameters to be (un)packed
   !     lpack() = 1 -> (un)pack, lpack() = 0 -> keep aktual value
   !     if a file "input.fitpar" exist in the current directory lpack() will be
   !     read from there
   
   !      print*,'entered ppack, nstring=',nstring
   !      print*,'nfit=',nfit
   
   pih=2.d0*atan(1.d0)
   if (nstring.eq.'init') then
      do i=1,maxdim
         lpack(i)=.false.
         pp(i) = 0.0d0
      enddo
      
      open(20,file='input.fitpar')
      !        let us parse the input.fitpar file line by line
      do iline=1,999
         read(20,'(a)',iostat=ierr)string
         
         if(ierr<0)then
            !                end of file, exit from the reading loop
            if(verbose) write(6,'(i3,a)')&
                 iline-1,' lines have been read. free params are:'
            exit
         elseif(ierr>0) then
            !                errors for reading strings should not happen
            write(6,*)'reading error for file input.fitpar'
            exit
         end if
         
         !             "auto", automatic assignment requested:
         if( index(string,'auto')>0)then
            if(verbose) write(6,*)&
                 'auto: free all params except for rloc and rcore'
            !                all nonzero psppar except for rloc are freed.
            !                unused r_l are assigned a default value of 1d0,
            !                so free them if they differ from that.
            !                all nlcc params are disabled per default.
            !                experimental: r_l2 are freed if they differ from r_l.
            !                experimantal: we still have gcore(1:4) for nlcc.
            !                lpack(6) = ( rcore > 0d0)
            do i=1,4
               lpack(1+i)=( gpot(i)**2> 1d-8)
               !lpack(6+i)=( gcore(i)**2> 1d-8)
               lpack(3+8*i)= (r_l(i)/=1d0)
               lpack(10+8*i)= ( r_l2(i)/=r_l(i) )
               do j=1,6
                  lpack(3+8*i+j)=  ((hsep(j,i,1)**2+&
                       hsep(j,i,2)**2) > 1d-8)
               end do
            end do
         end if    
         
         !             to avoid confusion, use qcore as an alias for gcore(1)
         if( index(string,'qcore')>0)then
            lpack(6+1)=.true.
            if(verbose) write(6,*) 'qcore corresponds to gcore1'
         end if
         
         
         !             check each line for occurences of any fitting
         !             parameter name as specified in the spack  array.
         !             if the name occurs at least once,
         !             the corresponding lpack is true
         do i=1,maxpar
            match=( index(string, trim(spack(i)) ) >0)
            lpack(i) = (lpack(i) .or. match ) 
            !                 if(match.and.verbose) write(6,*) spack(i)
         end do
      end do
      close(20)
      
      !        print some help if nothing matches
      if(.not.any(lpack) .and. verbose)then
         write(6,*)'could not recognize any keyword for free'
         write(6,*)'fitting variables in the file input.fitpar.'
         write(6,*)
         write(6,*)'valid keywords are the following:'
         write(6,*)'auto'
         write(6,'(5(2x,a))') spack(1:5)
         write(6,'(5(2x,a))') spack(6:10)
         write(6,'(8(2x,a))') spack(11:18)
         write(6,'(8(2x,a))') spack(19:26)
         write(6,'(8(2x,a))') spack(27:34)
         write(6,'(8(2x,a))') spack(35:42)
         write(6,*)
         write(6,*)'the keyword "auto" will free all parameters read from psppar and'
         write(6,*)'override previous keywords for params that are not present there.'
         write(6,*)'exception: does not free rloc.'
         write(6,*)
         
      end if
      
      !        test the highest angular momentum of the separable part and
      !        compare it with the lpx read from psppar
      lpxfit=0
      do i=1,4
         if( any(lpack(3+8*i: 10+8*i)) ) lpxfit= i
      end do
      if ( lpxfit .gt. lpx ) then
         !           this warning message does not need the error handler
         write(6,*)
         write(6,*)'                warning' 
         write(6,*)'in input.fitpar there are free parameters for the'
         write(6,*)'separable part for which no value has been'
         write(6,*)'assigned when reading psppar.'
         write(6,*)
         write(6,*)'the number of l components differ:'
         write(6,'(2(a,i3))')'psppar:',lpx,', input.fitpar:',lpxfit
         write(6,*)
         lpx = lpxfit
      endif
      
      !     pri,'lpack',lpack
      !     wri6,*)'lpack',lpack
      !
      !     if  projectors have to be orthogonalized do not count
      !     h(l2),h(l+1,4),h(l+1,5)  as fitting parameter 
      
      !
      if(ortprj.or.litprj) then
         !        due to gcore&gcore, indices of rl&hsep are shifted by +5
         do i=0,lpx-1
            if (lpack( (i*7+6) +7) ) then
               write(6,*) ' l=',i
               write(6,*) 'warning! transformation of hsep(): ',&
                    'no fit of h(1,2)'
               lpack((i*7+6) +4)=.false.
            endif
            if (lpack( (i*7+6) +9) ) then
               write(6,*) ' l=',i
               write(6,*) 'warning! transformation of hsep(): ',&
                    'no fit of h(1,3)'
               lpack((i*7+6) +6)=.false.
            endif
            if (lpack( (i*7+6) +10) ) then
               write(6,*) ' l=',i
               write(6,*) 'warning! transformation of hsep(): ',&
                    'no fit of h(2,3)'
               lpack((i*7+6) +10)=.false.
            endif
         enddo
      end if
      !        count fitting params --> simplex dimension
      nfit = 0
      do i=1,maxpar
         if (lpack(i)) then 
            if(verbose) write(6,*) spack(i)
            !              write(99,*) spack(i)
            !              also write to the dumpfile
            nfit = nfit +1
         endif
      enddo
      !
      !     double count projectors for nso=2 and l>0
      !
      !     index offset from new input.fitpar for core (5) and r_l2 (one per l)
      
      !         print*,'nfit before double count:',nfit
      if (nso.eq.2) then
         !     l=1   hsep(1,1)  at position 14+5+1
         do i=20,25
            if (lpack(i)) nfit =nfit + 1
         enddo
         !     l=2   hsep(1,1)  at position 21+5+2
         do i=28,33
            if (lpack(i)) nfit =nfit + 1
         enddo
         !     l=3   hsep(1,1)  at position 28+5+3
         do i=36,41
            if (lpack(i)) nfit =nfit + 1
         enddo
      endif
      !         print*,'nfit after double count:',nfit
      if (nfit.gt.maxdim)then
         write(6,*) 'aborting: maxdim,  nfit:',maxdim,nfit
         stop 'nfit > maxdim'
      end if
      !        if (maxdim.lt.nfit) stop
      !
      !     save initial parameter 
      !
      orloc = rloc
      do i=1,4
         ogpot(i) = gpot(i)
         ogcore(i)=gcore(i)
      enddo
      orcore=rcore
      do ll=1,lpx
         orl(ll)=r_l(ll)
         orl2(ll)=r_l2(ll)
         do i=1,min(2*ll,nso)
            do j=1,6
               ohsep(j,ll,i)= hsep(j,ll,i)
            enddo
            !     if the projectors are othonormalized we transform ohsep to the orthognomal basis
            !     and save the diagonal terms of the new ohsep() matrix
            if (ortprj) then
               h11=ohsep(1,ll,i)
               h12=ohsep(2,ll,i)
               h22=ohsep(3,ll,i)
               h13=ohsep(4,ll,i)
               h23=ohsep(5,ll,i)
               h33=ohsep(6,ll,i)
               if (ll.eq.1) then
                  
                  hh11=h11 + 1.549193338482967d0*h12 + 0.975900072948533d0*h13 &
                       + 0.6d0*h22 + 0.7559289460184545d0*h23 &
                       + 0.2380952380952381d0*h33
                  hh12=0.6324555320336759d0*h12 + 0.7968190728895957d0*h13 + &
                       0.4898979485566356d0*h22 + 0.925820099772551d0*h23 &
                       + 0.3888078956798695d0*h33
                  hh13=0.3563483225498992d0*h13 + 0.2760262237369417d0*h23 + &
                       0.173880176985767d0*h33
                  hh22=0.4d0*h22+1.007905261357939d0*h23 + 0.6349206349206349d0*h33
                  hh23=0.02839451399999733d0*(7.937253933193773d0*h23 + 10.d0*h33)
                  hh33=0.126984126984127d0*h33
                  
               elseif (ll.eq.2) then
                  
                  hh11=h11 + 1.690308509457033d0*h12 + 1.189176780021126d0*h13&
                       + 0.7142857142857144d0*h22 + 1.005037815259212d0*h23 &
                       + 0.3535353535353535d0*h33
                  hh12=0.5345224838248488d0*h12 + 0.7521014330903549d0*h13 &
                       + 0.4517539514526256d0*h22 + 0.953462589245592d0*h23 &
                       + 0.4471907802258314d0*h33
                  hh13=0.2842676218074805d0*h13 + 0.240249990052149d0*h23 &
                       + 0.1690222275826415d0*h33
                  hh22=0.2857142857142857d0*h22 + 0.80403025220737d0*h23 + &
                       0.5656565656565657d0*h33
                  hh23=0.01527129183875666d0*(9.9498743710662d0*h23 + 14.d0*h33)
                  hh33=0.0808080808080808d0*h33
                  
               elseif (ll.eq.3) then
                  
                  
                  hh11=h11 + 1.763834207376394d0*h12 + 1.327493036606129d0*h13 + &
                       0.7777777777777778d0*h22 + 1.170738814009927d0*h23 &
                       + 0.4405594405594406d0*h33
                  hh12=0.4714045207910317d0*h12 + 0.7095748751868991d0*h13 + &
                       0.4157397096415491d0*h22 + 0.938679328162116d0*h23 + &
                       0.4709778528806361d0*h33
                  hh13=0.236524958395633d0*h13 + 0.2085954062582479d0*h23 + &
                       0.1569926176268787d0*h33
                  hh12=0.4714045207910317d0*h12 + 0.7095748751868991d0*h13 + &
                       0.4157397096415491d0*h22 + 0.938679328162116d0*h23 + &
                       0.4709778528806361d0*h33
                  hh22=0.2222222222222222d0*h22 + 0.6689936080056727d0*h23 + &
                       0.5034965034965035d0*h33
                  hh23=0.00932400932400932d0*(11.9582607431014d0*h23 + 18.d0*h33)
                  hh33=0.05594405594405595d0*h33
                  
               elseif (ll.eq.4) then
                  
                  hh11=h11 + 1.809068067466582d0*h12 + 1.425050606388851d0*h13 + &
                       0.818181818181818d0*h22 + 1.289006773270979d0*h23 + &
                       0.5076923076923077d0*h33
                  hh12=0.0006593070220853591d0*(646.741834119303d0*h12 + &
                       1018.911183568028d0*h13 + 585.d0*h22 + &
                       1382.459764333125d0*h23 + 726.d0*h33)
                  hh13=0.2025478734167333d0*h13 + 0.1832114449657378d0*h23 + &
                       0.144320484917644d0*h33
                  hh22=0.1818181818181818d0*h22 + 0.5728918992315464d0*h23 + &
                       0.4512820512820513d0*h33
                  hh23=0.006184848093902844d0*(13.96424004376894d0*h23+22.d0*h33)
                  hh33=0.04102564102564103d0*h33
                  
               endif
               ohsep(1,ll,i)=hh11
               ohsep(2,ll,i)=hh12
               ohsep(3,ll,i)=hh22
               ohsep(4,ll,i)=hh13
               ohsep(5,ll,i)=hh23
               ohsep(6,ll,i)=hh33
            endif
         enddo
      enddo
   endif
   !
   !------------------------------------------------------------------
   !
   !     pack the parameters into the array pp()
   !
   if (nstring.eq.'pack'.or.nstring.eq.'init') then
      !         print*,'nfit=',nfit
      do i=1,nfit
         pp(i)=0.0d0
      enddo
      !         write(*,*) 'packed array pp:'
      !         write(*,'(20e8.2)') (pp(i),i=1,nfit)
   endif
   !
   !------------------------------------------------------------------
   !
   !     unpack array pp()
   !
   
   !     frequently used weighting scheme: factor bound to +/- 10% 
   !                                       using  arctan(pp)/5pi 
   
   if (nstring.eq.'unpack') then
      !        write(*,*) 'unpacking array pp:'
      !        write(*,'(20e8.2)') (pp(i),i=1,nfit)
      ncount = 1
      !     rloc
      if (lpack(1)) then 
         rloc = orloc+.10d0*orloc*atan(pp(ncount))/pih 
         ncount = ncount + 1
      endif  
      !     gpot(1-4)
      do i=1,4
         if (lpack(1+i)) then 
            gpot(i) = ogpot(i)+pp(ncount)
            ncount = ncount + 1
         endif  
      enddo
      !        packing of gcore and rcore may need other conventions;
      !        let us handle rcore like other radii 
      if (lpack(6)) then 
         rcore = orcore+.10d0*orcore*atan(pp(ncount))/pih 
         ncount=ncount+1
      end if
      !        gcore(1-4)
      !        we need to be careful with the moves for the nlcc coeffs...
      !        a slight deviation can add a lot of charge
      !        let us consider pp a percentage to alter the ogcore
      !        downside: we cannot change the sign of a gcore curing the fit.
      do i=1,4
         if (lpack(6+i)) then 
            gcore(i) = ogcore(i) *( 1d0 + 0.01d0 * pp(ncount) )
            ncount = ncount + 1
         endif  
      enddo
      
      
      
      !     projectors for l=0: r_l(1),hsep(1,1-6)
      !        in the polarized case, be sure to keep the down
      !        component equal to the up component.
      !        this is done after the following section with
      !        explicit treatment for each l-component.
      increm=6+5
      if ( lpack(increm) ) then 
         r_l(1) = orl(1)+.10d0*orl(1)*atan(pp(ncount))/pih 
         ncount = ncount + 1
      endif  
      do i=1,6
         if (lpack(increm+i)) then 
            hsep(i,1,1)= ohsep(i,1,1)+pp(ncount)
            !              ncount = ncount + nspol
            ncount = ncount + 1
         else
            hsep(i,1,1)= ohsep(i,1,1)
         endif  
      enddo
      if ( lpack(increm+7) ) then 
         r_l2(1) = orl2(1)+.10d0*orl2(1)*atan(pp(ncount))/pih 
         ncount = ncount + 1
      else
         !           convention: if r_l2 is not free but r_l is, 
         !                       then r_l2 must be disabled 
         if(lpack(increm)) r_l2(1) = r_l(1) 
      endif  
      if (ortprj) then
         !     do back transformation from orthonormal projectors to 
         !     unnormalized projectors:
         h11=hsep(1,1,1)
         h22=hsep(3,1,1)
         h33=hsep(6,1,1)
         
         hh11=h11 + 1.5d0*h22 + 1.875d0*h33
         hh12=-1.936491673103709d0*h22 - 4.841229182759272d0*h33
         hh13=3.842606537234849d0*h33
         hh22=2.5d0*h22 + 12.5d0*h33
         hh23=-9.92156741649221d0*h33
         hh33=7.875d0*h33
         
         hsep(1,1,1)=hh11
         hsep(2,1,1)=hh12
         hsep(3,1,1)=hh22
         hsep(4,1,1)=hh13
         hsep(5,1,1)=hh23
         hsep(6,1,1)=hh33
      else if (litprj) then
         hsep(2,1,1)= -0.5d0 * sqrt(3.d0/5.d0)  *hsep(3,1,1)
         hsep(4,1,1)= +0.5d0 * sqrt(5.d0/21.d0) *hsep(6,1,1)
         hsep(5,1,1)= -0.5d0 * sqrt(100.d0/63.0d0) *hsep(6,1,1)
      endif
      !     projectors for l=1: r_l(1),hsep(1,1-6,1-2)
      increm=13+5+1
      if (lpack(increm)) then 
         r_l(2) = orl(2)+.10d0*orl(2)*atan(pp(ncount))/pih 
         ncount = ncount + 1
      endif  
      do i=1,6
         if (lpack(increm+i)) then 
            hsep(i,2,1)= ohsep(i,2,1) + pp(ncount)
            if (nso.eq.2) hsep(i,2,2)= ohsep(i,2,2)+pp(ncount+1)
            ncount = ncount + nso
         else
            hsep(i,2,1)= ohsep(i,2,1)
            !              relativistic : copy back fixed pars
            if (nso.eq.2) hsep(i,2,2)= ohsep(i,2,2)
         endif  
      enddo
      if ( lpack(increm+7) ) then 
         r_l2(2) = orl2(2)+.10d0*orl2(2)*atan(pp(ncount))/pih 
         ncount = ncount + 1
      else
         !           convention: if r_l2 is not free but r_l is, 
         !                       then r_l2 must be disabled 
         if(lpack(increm)) r_l2(2) = r_l(2) 
      endif  
      if (ortprj) then
         do i=1,nso
            !     
            !     do back transformation from orthonormal projectors to 
            !     unnormalized projectors:
            h11=hsep(1,2,i)
            h22=hsep(3,2,i)
            h33=hsep(6,2,i)
            
            hh11=h11 + 2.5d0*h22 + 4.375d0*h33
            hh12=-2.958039891549808d0*h22 - 10.35313962042433d0*h33
            hh13=7.358031326380719d0*h33
            hh22=3.5d0*h22 + 24.5d0*h33
            hh23=-17.41228014936585d0*h33
            hh33=12.375d0*h33
            
            hsep(1,2,i)=hh11
            hsep(2,2,i)=hh12
            hsep(3,2,i)=hh22
            hsep(4,2,i)=hh13
            hsep(5,2,i)=hh23
            hsep(6,2,i)=hh33
            
         enddo
      endif
      if (litprj) then
         do i=1,nso
            hsep(2,2,i)= -0.5d0*sqrt(5.d0/7.d0)       *hsep(3,2,i)
            hsep(4,2,i)=  0.5d0*sqrt(35.d0/11.d0)/3.d0*hsep(6,2,i)
            hsep(5,2,i)= -0.5d0*14.d0/sqrt(11.d0)/3.d0*hsep(6,2,i)
         enddo
      endif
      !     projectors for l=2: r_l(1),hsep(1,1-6,1-2)
      increm=20+5+2
      if (lpack(increm)) then 
         r_l(3) = orl(3)+.10d0*orl(3)*atan(pp(ncount))/pih 
         ncount = ncount + 1
      endif  
      do i=1,6
         if (lpack(increm+i)) then 
            hsep(i,3,1)= ohsep(i,3,1)+pp(ncount)
            if (nso.eq.2) hsep(i,3,2)= ohsep(i,3,2)+pp(ncount+1)
            ncount = ncount + nso
         else
            hsep(i,3,1)= ohsep(i,3,1)
            !              relativistic : copy back fixed pars
            if (nso.eq.2) hsep(i,3,2)= ohsep(i,3,2)
         endif  
      enddo
      if ( lpack(increm+7) ) then 
         r_l2(3) = orl2(3)+.10d0*orl2(3)*atan(pp(ncount))/pih 
         ncount = ncount + 1
      else
         !           convention: if r_l2 is not free but r_l is, 
         !                       then r_l2 must be disabled 
         if(lpack(increm)) r_l2(3) = r_l(3) 
      endif  
      if (ortprj) then
         do i=1,nso
            !     
            !     do back transformation from orthonormal projectors to 
            !     unnormalized projectors:
            h11=hsep(1,3,i)
            h22=hsep(3,3,i)
            h33=hsep(6,3,i)
            
            
            hh11=h11 + 3.5d0*h22 + 7.875d0*h33
            hh12=-3.968626966596886d0*h22-17.85882134968598d0*h33
            hh13=11.86446901466728d0*h33
            hh22=4.5d0*h22 + 40.5d0*h33
            hh23=-26.90608667197814d0*h33
            hh33=17.875d0*h33
            
            hsep(1,3,i)=hh11
            hsep(2,3,i)=hh12
            hsep(3,3,i)=hh22
            hsep(4,3,i)=hh13
            hsep(5,3,i)=hh23
            hsep(6,3,i)=hh33
            
         enddo
      endif
      if (litprj) then
         do i=1,nso
            hsep(2,3,i)= -0.5d0*sqrt(7.d0/9.d0)   *hsep(3,3,i)
            hsep(4,3,i)=  0.5d0*3.d0*sqrt(7.0d0/143.d0)*hsep(6,3,i)
            hsep(5,3,i)= -0.5d0*18.d0*sqrt(1/143.0d0)*hsep(6,3,i)
         enddo
      endif
      !     projectors for l=3: r_l(1),hsep(1,1-6,1-2)
      increm=27+5+3
      if (lpack(increm)) then 
         r_l(4) = orl(4)+.10d0*orl(4)*atan(pp(ncount))/pih 
         ncount = ncount + 1
      endif  
      do i=1,6
         if (lpack(increm+i)) then 
            hsep(i,4,1)= ohsep(i,4,1)+pp(ncount)
            if (nso.eq.2) hsep(i,4,2)= ohsep(i,4,2)+pp(ncount+1)
            ncount = ncount + nso
         else
            hsep(i,4,1)= ohsep(i,4,1)
            !              relativistic : copy back fixed pars
            if (nso.eq.2) hsep(i,4,2)= ohsep(i,4,2)
         endif  
      enddo
      if ( lpack(increm+7) ) then 
         r_l2(4) = orl2(4)+.10d0*orl2(4)*atan(pp(ncount))/pih 
         ncount = ncount + 1
      else
         !           convention: if r_l2 is not free but r_l is, 
         !                       then r_l2 must be disabled 
         if(lpack(increm)) r_l2(4) = r_l(4) 
      endif  
      if (ortprj) then
         do i=1,nso
            !     
            !     do back transformation from orthonormal projectors to 
            !     unnormalized projectors:
            h11=hsep(1,4,i)
            h22=hsep(3,4,i)
            h33=hsep(6,4,i)
            
            hh11=h11 + 4.5d0*h22 + 12.375d0*h33
            hh12=-4.9749371855331d0*h22 - 27.36215452043205d0*h33
            hh13=17.36780426536412d0*h33
            hh22=5.5d0*h22 + 60.5d0*h33
            hh23=-38.40166012036458d0*h33
            hh33=24.375d0*h33
            
            hsep(2,4,i)=hh12
            hsep(3,4,i)=hh22
            hsep(4,4,i)=hh13
            hsep(5,4,i)=hh23
            hsep(6,4,i)=hh33
            
         enddo
      endif
      !        important: use polarized psps only in the relativistic case!
      !        in the polarized, nonrelativistic case, the projectors
      !        must be the same for up and down states!
      !        below line is to be on the safe side.
      if(nspol>nso)  hsep(:,:,2)= hsep(:,:,1)
      
      if (litprj) then
         do i=1,nso
            hsep(2,4,i)= -0.5d0*sqrt(9.d0/11.d0)  *hsep(3,4,i)
            hsep(4,4,i)= +0.5d0*sqrt(33.d0/65.d0) *hsep(6,4,i)
            hsep(5,4,i)= -0.5d0*22.d0/sqrt(195.d0)*hsep(6,4,i)
         enddo
      endif
      !
      !     if avgl is set: modify hsep() so that average potential for the highest
      !     projector is zero; if the projectors have to be orthogonalized modify
      !     also the corcesponding offdialgonal elements
      !     (only for relativistic calculations)
      if (nso.eq.2) then
         nhsep1=0
         nhsep2=0
         nhsep3=0
         do i=1,6
            if (hsep(i,2,1).ne.0) nhsep1=i
            if (hsep(i,3,1).ne.0) nhsep2=i
            if (hsep(i,4,1).ne.0) nhsep3=i
         enddo
         if (avgl1) then
            nnmin=1
            if (nhsep1.eq.3 .and. litprj) nnmin=2
            if (nhsep1.eq.6 .and. litprj) nnmin=4
            if (.not.litprj) nnmin=nhsep1
            if (ortprj) then
               !     for orthogonal projectors remove average-part form the highest 
               !     orthonormal projector 
               do i=1,nso
                  h11=hsep(1,2,i)
                  h12=hsep(2,2,i)
                  h22=hsep(3,2,i)
                  h13=hsep(4,2,i)
                  h23=hsep(5,2,i)
                  h33=hsep(6,2,i)
                  hh11=h11 + 1.690308509457033d0*h12 + 1.189176780021126d0*h13&
                       + 0.7142857142857144d0*h22 + 1.005037815259212d0*h23 &
                       + 0.3535353535353535d0*h33
                  hh22=0.2857142857142857d0*h22 + 0.80403025220737d0*h23 + &
                       0.5656565656565657d0*h33
                  hh33=0.0808080808080808d0*h33
                  hsep(1,2,i)=hh11
                  hsep(3,2,i)=hh22
                  hsep(6,2,i)=hh33
               enddo
            endif
            do nn=nnmin,nhsep1
               xshift= hsep(nn,2,2)+2.0d0*hsep(nn,2,1)
               hsep(nn,2,1) =  hsep(nn,2,1) - xshift*1.d0/3.d0
               hsep(nn,2,2) =  hsep(nn,2,2) - xshift*1.d0/3.d0
            enddo
            if (ortprj) then
               do i=1,nso
                  !     do back transformation of hsep()
                  h11=hsep(1,2,i)
                  h22=hsep(3,2,i)
                  h33=hsep(6,2,i)
                  
                  hh11=h11 + 2.5d0*h22 + 4.375d0*h33
                  hh12=-2.958039891549808d0*h22 - 10.35313962042433d0*h33
                  hh13=7.358031326380719d0*h33
                  hh22=3.5d0*h22 + 24.5d0*h33
                  hh23=-17.41228014936585d0*h33
                  hh33=12.375d0*h33
                  
                  hsep(1,2,i)=hh11
                  hsep(2,2,i)=hh12
                  hsep(3,2,i)=hh22
                  hsep(4,2,i)=hh13
                  hsep(5,2,i)=hh23
                  hsep(6,2,i)=hh33
                  
               enddo
            endif
         endif
         if (avgl2) then
            nnmin=1
            if (nhsep1.eq.3 .and. litprj) nnmin=2
            if (nhsep1.eq.6 .and. litprj) nnmin=4
            if (.not.litprj) nnmin=nhsep2
            if (ortprj) then
               !     for orthogonal projectors remove average-part form the highest 
               !     orthonormal projector 
               do i=1,nso
                  h11=hsep(1,3,i)
                  h12=hsep(2,3,i)
                  h22=hsep(3,3,i)
                  h13=hsep(4,3,i)
                  h23=hsep(5,3,i)
                  h33=hsep(6,3,i)
                  hh11=h11 + 1.763834207376394d0*h12 + 1.327493036606129d0*h13 + &
                       0.7777777777777778d0*h22 + 1.170738814009927d0*h23 &
                       + 0.4405594405594406d0*h33
                  hh12=0.4714045207910317d0*h12 + 0.7095748751868991d0*h13 + &
                       0.4157397096415491d0*h22 + 0.938679328162116d0*h23 + &
                       0.4709778528806361d0*h33
                  hh13=0.236524958395633d0*h13 + 0.2085954062582479d0*h23 + &
                       0.1569926176268787d0*h33
                  hh12=0.4714045207910317d0*h12 + 0.7095748751868991d0*h13 + &
                       0.4157397096415491d0*h22 + 0.938679328162116d0*h23 + &
                       0.4709778528806361d0*h33
                  hh22=0.2222222222222222d0*h22 + 0.6689936080056727d0*h23 + &
                       0.5034965034965035d0*h33
                  hh23=0.00932400932400932d0*(11.9582607431014d0*h23 + 18.d0*h33)
                  hh33=0.05594405594405595d0*h33
                  
                  
                  hsep(1,3,i)=hh11
                  hsep(3,3,i)=hh22
                  hsep(6,3,i)=hh33
               enddo
            endif
            
            do nn=nnmin,nhsep2
               xshift= 2*hsep(nn,3,2)+3.0d0*hsep(nn,3,1)
               hsep(nn,3,1) =  hsep(nn,3,1) - xshift*1.d0/5.d0
               hsep(nn,3,2) =  hsep(nn,3,2) - xshift*1.d0/5.d0
            enddo
            if (ortprj) then
               do i=1,nso
                  !     do back transformation of hsep()
                  h11=hsep(1,3,i)
                  h22=hsep(3,3,i)
                  h33=hsep(6,3,i)
                  
                  hh11=h11 + 3.5d0*h22 + 7.875d0*h33
                  hh12=-3.968626966596886d0*h22-17.85882134968598d0*h33
                  hh13=11.86446901466728d0*h33
                  hh22=4.5d0*h22 + 40.5d0*h33
                  hh23=-26.90608667197814d0*h33
                  hh33=17.875d0*h33
                  
                  hsep(1,3,i)=hh11
                  hsep(2,3,i)=hh12
                  hsep(3,3,i)=hh22
                  hsep(4,3,i)=hh13
                  hsep(5,3,i)=hh23
                  hsep(6,3,i)=hh33
               enddo
            endif
         endif
         if (avgl3) then
            nnmin=1
            if (nhsep1.eq.3 .and. litprj) nnmin=2
            if (nhsep1.eq.6 .and. litprj) nnmin=4
            if (.not.litprj) nnmin=nhsep3
            if (ortprj) then
               !     for orthogonal projectors remove average-part form the highest 
               !     orthonormal projector 
               do i=1,nso
                  h11=hsep(1,4,i)
                  h12=hsep(2,4,i)
                  h22=hsep(3,4,i)
                  h13=hsep(4,4,i)
                  h23=hsep(5,4,i)
                  h33=hsep(6,4,i)
                  
                  hh11=h11 + 1.809068067466582d0*h12 + 1.425050606388851d0*h13 + &
                       0.818181818181818d0*h22 + 1.289006773270979d0*h23 + &
                       0.5076923076923077d0*h33
                  hh12=0.0006593070220853591d0*(646.741834119303d0*h12 + &
                       1018.911183568028d0*h13 + 585.d0*h22 + &
                       1382.459764333125d0*h23 + 726.d0*h33)
                  hh13=0.2025478734167333d0*h13 + 0.1832114449657378d0*h23 + &
                       0.144320484917644d0*h33
                  hh22=0.1818181818181818d0*h22 + 0.5728918992315464d0*h23 + &
                       0.4512820512820513d0*h33
                  hh23=0.006184848093902844d0*(13.96424004376894d0*h23+22.d0*h33)
                  hh33=0.04102564102564103d0*h33
                  
                  hsep(1,4,i)=hh11
                  hsep(3,4,i)=hh22
                  hsep(6,4,i)=hh33
               enddo
            endif
            do nn=nnmin,nn
               xshift= 5*hsep(nn,4,2)+7.0d0*hsep(nn,4,1)
               hsep(nn,4,1) =  hsep(nn,4,1) - xshift*1.d0/7.d0
               hsep(nn,4,2) =  hsep(nn,4,2) - xshift*1.d0/7.d0
            enddo
            if (ortprj) then
               do i=1,nso
                  !     do back transformation of hsep()
                  h11=hsep(1,4,i)
                  h22=hsep(3,4,i)
                  h33=hsep(6,4,i)
                  
                  hh11=h11 + 4.5d0*h22 + 12.375d0*h33
                  hh12=-4.9749371855331d0*h22 - 27.36215452043205d0*h33
                  hh13=17.36780426536412d0*h33
                  hh22=5.5d0*h22 + 60.5d0*h33
                  hh23=-38.40166012036458d0*h33
                  hh33=24.375d0*h33
                  
                  hsep(1,4,i)=hh11
                  hsep(2,4,i)=hh12
                  hsep(3,4,i)=hh22
                  hsep(4,4,i)=hh13
                  hsep(5,4,i)=hh23
                  hsep(6,4,i)=hh33
               enddo
            endif
         endif
         
      endif
   endif  
   !------------------------------------------------------------------
   !      print*,'leave ppack with nfit=',nfit
   return
end subroutine ppack
