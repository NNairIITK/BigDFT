!{\src2tex{textfont=tt}}
!!****f* ABINIT/ewald2
!!
!! NAME
!! ewald2
!!
!! FUNCTION
!! Compute the part of the stress tensor coming from the Ewald energy
!! which is calculated by derivating the Ewald energy with respect to
!! strain.
!! See Nielsen and Martin, Phys. Rev. B 32, 3792 (1985).
!! Definition of stress tensor is $(1/ucvol)*d(Etot)/d(strain(a,b))$.
!!
!! COPYRIGHT
!! Copyright (C) 1998-2010 ABINIT group (JCC, DCA, XG, GMR)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! gmet(3,3)=metric tensor in reciprocal space (bohr^-2)
!! natom=number of atoms in umit cell
!! ntypat=number of type of atoms
!! rmet(3,3)=metric tensor in real space (bohr^2) (inverse transpose of gmet)
!! rprimd(3,3)=dimensional primitive translations in real space (bohr)
!! typat(natom)=integer label of each type of atom (1,2,...)
!! ucvol=unit cell volume (bohr^3)
!! xred(3,natom)=relative coords of atoms in unit cell (dimensionless)
!! zion(ntypat)=charge on each type of atom (real number)
!!
!! OUTPUT
!! $stress(6)=(1/ucvol)*gradient$ of Ewald energy with respect to strain,
!!      in hartrees/bohr^3
!! Cartesian components of stress are provided for this symmetric
!! tensor in the order 11 22 33 32 31 21.
!!
!! PARENTS
!!      stress
!!
!! CHILDREN
!!      derfc,matr3inv
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.inc"
#endif

subroutine ewald2(iproc,nproc,gmet,natom,ntypat,rmet,rprimd,stress,&
&  typat,ucvol,xred,zion)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
 use interfaces_32_util
!End of the abilint section

 implicit none
 !SM there are probably better ways than this...
 include 'mpif.h'

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: iproc,nproc,natom,ntypat
 real(dp),intent(in) :: ucvol
!arrays
 integer,intent(in) :: typat(natom)
 real(dp),intent(in) :: gmet(3,3),rmet(3,3),rprimd(3,3),xred(3,natom)
 real(dp),intent(in) :: zion(ntypat)
 real(dp),intent(out) :: stress(6)

!Local variables-------------------------------
!scalars
 integer :: ia,ib,ig1,ig2,ig3,ir1,ir2,ir3,newg,newr,ng,nr
 real(dp) :: arg1,arg2,arg3,ch,dderfc,derfc_arg,direct,eta,fac,fraca1
 real(dp) :: fraca2,fraca3,fracb1,fracb2,fracb3,g1,g2,g3,gsq,r1,r1c,r2,r2c
 real(dp) :: r3,r3c,recip,reta,rmagn,rsq,summi,summr,t1,t2,t3,t4,t5,t6,term1
 real(dp) :: term2,term3,term4
!arrays
 real(dp) :: gprimd(3,3),strg(6),strr(6)
 real(dp) :: tt
 integer :: natp, isat, ii, iia, ierr
 real(dp),dimension(6) :: strr_tmp

! *************************************************************************


!SM: MPI parallelization over the atoms
 tt = real(natom,kind=dp)/nproc
 natp = floor(tt) !number of atoms per proc
 isat = iproc*natp !offset for each proc
 ii = natom-nproc*natp !remaining atoms
 if (iproc<ii) then
     natp = natp+1 !one more atom for this proc
     isat = isat+iproc !offset increases by the number of additional atoms up to iproc
 else
     isat = isat+ii !offset increases by the number of additional atoms
 end if
 ! Check
 ii = natp
 if (nproc>1) then
     !call mpiallred(ii, 1, mpi_sum, comm=bigdft_mpi%mpi_comm)
     iia=0
     call mpi_allreduce(ii, iia, 1, mpi_integer, mpi_sum, mpi_comm_world, ierr)
     ii = iia
 end if
 !if (ii/=natom) call f_err_throw('ii/=natom',err_name='BIGDFT_RUNTIME_ERROR')
 if (ii/=natom) stop 'ii/=natom'

!Define dimensional reciprocal space primitive translations gprimd
!(inverse transpose of rprimd)
 call matr3inv(rprimd,gprimd)

!Add up total charge and sum of charge^2 in cell
 ch=0._dp
 !$omp parallel default(none) shared(natom,zion,typat,ch) private(ia)
 !$omp do reduction(+:ch) schedule(static)
 do ia=1,natom
   ch=ch+zion(typat(ia))
 end do
 !$omp end do
 !$omp end parallel

!Compute eta, the Ewald summation convergence parameter,
!for approximately optimized summations:
 direct=rmet(1,1)+rmet(1,2)+rmet(1,3)+rmet(2,1)+&
& rmet(2,2)+rmet(2,3)+rmet(3,1)+rmet(3,2)+rmet(3,3)
 recip=gmet(1,1)+gmet(1,2)+gmet(1,3)+gmet(2,1)+&
& gmet(2,2)+gmet(2,3)+gmet(3,1)+gmet(3,2)+gmet(3,3)
!Here, a bias is introduced, because G-space summation scales
!better than r space summation !
 eta=pi*200.0_dp/33.0_dp*sqrt(1.69_dp*recip/direct)

 fac=pi**2/eta

!Conduct reciprocal space summations
 strg(1:6)=0.0_dp

!Sum over G space, done shell after shell until all
!contributions are too small
 ng=0
 do
   ng=ng+1
   newg=0

   do ig3=-ng,ng
     do ig2=-ng,ng
       do ig1=-ng,ng

!        Exclude shells previously summed over
         if(abs(ig1)==ng .or. abs(ig2)==ng .or. abs(ig3)==ng&
&         .or. ng==1 ) then

!          Compute Cartesian components of each G
           g1=gprimd(1,1)*ig1+gprimd(1,2)*ig2+gprimd(1,3)*ig3
           g2=gprimd(2,1)*ig1+gprimd(2,2)*ig2+gprimd(2,3)*ig3
           g3=gprimd(3,1)*ig1+gprimd(3,2)*ig2+gprimd(3,3)*ig3
!          Compute |G|^2 (no pi factors)
           gsq=(g1**2+g2**2+g3**2)

!          skip g=0:
           if (gsq>1.0d-20) then
             arg1=fac*gsq

!            larger arg1 gives 0 contribution because of exp(-arg1)
             if (arg1<=80._dp) then
!              When any term contributes then include next shell
               newg=1
               term1=exp(-arg1)/arg1
               summr = 0.0_dp
               summi = 0.0_dp
               !$omp parallel default(none) &
               !$omp shared(natom,ig1,ig2,ig3,xred,zion,typat,summr,summi)private(ia,arg2)
               !$omp do reduction(+:summr,summi) schedule(static)
               do ia=1,natom
                 arg2=two_pi*(ig1*xred(1,ia)+ig2*xred(2,ia)+ig3*xred(3,ia))
!                Sum real and imaginary parts (avoid complex variables)
                 summr=summr+zion(typat(ia))*cos(arg2)
                 summi=summi+zion(typat(ia))*sin(arg2)
               end do
               !$omp end do
               !$omp end parallel

!              Avoid underflow error messages
               if (abs(summr)<1.d-16) summr=0.0_dp
               if (abs(summi)<1.d-16) summi=0.0_dp

               term2=(2._dp/gsq)*(1._dp+arg1)
               t1=term2*g1*g1-1._dp
               t2=term2*g2*g2-1._dp
               t3=term2*g3*g3-1._dp
               t4=term2*g2*g3
               t5=term2*g1*g3
               t6=term2*g1*g2
               term3=term1*(summr*summr+summi*summi)
               strg(1)=strg(1)+t1*term3
               strg(2)=strg(2)+t2*term3
               strg(3)=strg(3)+t3*term3
               strg(4)=strg(4)+t4*term3
               strg(5)=strg(5)+t5*term3
               strg(6)=strg(6)+t6*term3

!              End condition not being larger than 80.0
             end if

!            End skip g=0
           end if

!          End triple loop and condition of new shell
         end if
       end do
     end do
   end do

!  Check if new shell must be calculated
   if (newg==0) exit

!  End loop on new shell. Note that there is an "exit" instruction within the loop
 end do


!Conduct real space summations
 reta=sqrt(eta)
 strr(1:6)=0.0_dp
 strr_tmp(1:6)=0.0_dp

!Loop on shells in r-space as was done in g-space
 nr=0
 do
   nr=nr+1
   newr=0

   do ir3=-nr,nr
     do ir2=-nr,nr
       do ir1=-nr,nr
         if( abs(ir3)==nr .or. abs(ir2)==nr .or. abs(ir1)==nr&
&         .or. nr==1 )then

           do ia=1,natp!natom
             iia=isat+ia
!            Convert reduced atomic coordinates to [0,1)
             fraca1=xred(1,iia)-aint(xred(1,iia))+0.5_dp-sign(0.5_dp,xred(1,iia))
             fraca2=xred(2,iia)-aint(xred(2,iia))+0.5_dp-sign(0.5_dp,xred(2,iia))
             fraca3=xred(3,iia)-aint(xred(3,iia))+0.5_dp-sign(0.5_dp,xred(3,iia))
             !$omp parallel default(none) &
             !$omp shared(natom,xred,fraca1,fraca2,fraca3,ir1,ir2,ir3,rprimd,reta,strr_tmp,newr,eta,zion,typat,iia) &
             !$omp private(ib,fracb1,fracb2,fracb3,r1,r2,r3,r1c,r2c,r3c,rsq,rmagn,arg3,term3,term4,dderfc,derfc_arg)
             !$omp do reduction(+:strr_tmp,newr)
             do ib=1,natom
               fracb1=xred(1,ib)-aint(xred(1,ib))+0.5_dp-sign(0.5_dp,xred(1,ib))
               fracb2=xred(2,ib)-aint(xred(2,ib))+0.5_dp-sign(0.5_dp,xred(2,ib))
               fracb3=xred(3,ib)-aint(xred(3,ib))+0.5_dp-sign(0.5_dp,xred(3,ib))
               r1=ir1+fracb1-fraca1
               r2=ir2+fracb2-fraca2
               r3=ir3+fracb3-fraca3
!              Convert from reduced to cartesian coordinates
               r1c=rprimd(1,1)*r1+rprimd(1,2)*r2+rprimd(1,3)*r3
               r2c=rprimd(2,1)*r1+rprimd(2,2)*r2+rprimd(2,3)*r3
               r3c=rprimd(3,1)*r1+rprimd(3,2)*r2+rprimd(3,3)*r3
!              Compute |r|^2
               rsq=r1c**2+r2c**2+r3c**2
               rmagn=sqrt(rsq)

!              Avoid zero denominators in 'term':
               if (rmagn>=1.0d-12) then

!                Note: erfc(8) is about 1.1e-29,
!                so do not bother with larger arg.
!                Also: exp(-64) is about 1.6e-28,
!                so do not bother with larger arg**2 in exp.
                 arg3=reta*rmagn
                 if (arg3<8.0_dp) then
                   newr=newr+1
!                  derfc computes the complementary error function
!                  dderfc is the derivative of the complementary error function
                   dderfc=(-2/sqrt(pi))*exp(-eta*rsq)
                   call derfcf(derfc_arg,arg3)
                   term3=dderfc-derfc_arg/arg3
                   term4=zion(typat(iia))*zion(typat(ib))*term3
                   strr_tmp(1)=strr_tmp(1)+term4*r1c*r1c/rsq
                   strr_tmp(2)=strr_tmp(2)+term4*r2c*r2c/rsq
                   strr_tmp(3)=strr_tmp(3)+term4*r3c*r3c/rsq
                   strr_tmp(4)=strr_tmp(4)+term4*r2c*r3c/rsq
                   strr_tmp(5)=strr_tmp(5)+term4*r1c*r3c/rsq
                   strr_tmp(6)=strr_tmp(6)+term4*r1c*r2c/rsq
!                  End the condition of not being to large
                 end if

!                End avoid zero denominator
               end if

!              End loop over ib:
             end do
             !$omp end do
             !$omp end parallel

!            End loop over ia:
           end do

!          End triple loop overs real space points, and associated new shell condition
         end if
       end do
     end do
   end do

   if (nproc>1) then
     !call mpiallred(newr, mpi_sum, bigdft_mpi%mpi_comm)
     ii = 0
     call mpi_allreduce(newr, ii, 1, &
          mpi_integer, mpi_sum, mpi_comm_world, ierr)
     newr=ii
   end if

!  Check if new shell must be calculated
   if(newr==0) exit

!  End loop on new shells
 end do

 if (nproc>1) then
   !call mpiallred(stress_tmp, mpi_sum, bigdft_mpi%mpi_comm)
   strr=0.0_dp
   call mpi_allreduce(strr_tmp, strr, 6, &
        mpi_double_precision, mpi_sum, mpi_comm_world, ierr)
 else
   strr=strr_tmp
 end if

!Finally assemble stress tensor coming from Ewald energy, stress
!(note division by unit cell volume in accordance with definition
!found in Nielsen and Martin, Phys. Rev. B 32, 3792 (1985).)

 fac = pi/(2._dp*ucvol*eta)
 stress(1)=(0.5_dp*reta*strr(1)+fac*(strg(1)+(ch**2)))/ucvol
 stress(2)=(0.5_dp*reta*strr(2)+fac*(strg(2)+(ch**2)))/ucvol
 stress(3)=(0.5_dp*reta*strr(3)+fac*(strg(3)+(ch**2)))/ucvol
 stress(4)=(0.5_dp*reta*strr(4)+fac*strg(4))/ucvol
 stress(5)=(0.5_dp*reta*strr(5)+fac*strg(5))/ucvol
 stress(6)=(0.5_dp*reta*strr(6)+fac*strg(6))/ucvol

end subroutine ewald2
!!***
