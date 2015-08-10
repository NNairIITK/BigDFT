!> @file
!!    Main routine to perform Poisson solver calculation
!! @author
!!    Creation date: February 2007
!!    Luigi Genovese
!!    Copyright (C) 2002-2013 BigDFT group 
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS


!> Calculate the Hartree potential by solving the Poisson equation 
!! @f$\nabla^2 V(x,y,z)=-4 \pi \rho(x,y,z)@f$
!! from a given @f$\rho@f$, 
!! for different boundary conditions an for different data distributions.
!! Following the boundary conditions, it applies the Poisson Kernel previously calculated.
!!    
!! @warning
!!    The dimensions of the arrays must be compatible with geocode, datacode, nproc, 
!!    ixc and iproc. Since the arguments of these routines are indicated with the *, it
!!    is IMPERATIVE to use the PS_dim4allocation routine for calculation arrays sizes.
!!    Moreover, for the cases with the exchange and correlation the density must be initialised
!!    to 10^-20 and not to zero.
!!
!! @todo
!!    Wire boundary condition is missing
subroutine H_potential(datacode,kernel,rhopot,pot_ion,eh,offset,sumpion,&
      quiet,rho_ion,stress_tensor) !optional argument
   implicit none

   !> kernel of the Poisson equation. It is provided in distributed case, with
   !! dimensions that are related to the output of the PS_dim4allocation routine
   !! it MUST be created by following the same geocode as the Poisson Solver.
   !! it is declared as inout as in the cavity treatment it might be modified
   type(coulomb_operator), intent(inout) :: kernel

   !> @copydoc poisson_solver::doc::datacode
   !! To be used only in the periodic case, ignored for other boundary conditions.
   character(len=1), intent(in) :: datacode

   !> Logical value which states whether to sum pot_ion to the final result or not
   !!   .true.  rhopot will be the Hartree potential + pot_ion
   !!           pot_ion will be untouched
   !!   .false. rhopot will be only the Hartree potential
   !!           pot_ion will be ignored
   logical, intent(in) :: sumpion

   !> Total integral on the supercell of the final potential on output
   real(dp), intent(in) :: offset

   !> Hartree Energy (Hartree)
   real(gp), intent(out) :: eh

   !> On input, it represents the density values on the grid points
   !! On output, it is the Hartree potential (and maybe also pot_ion)
   real(dp), dimension(*), intent(inout) :: rhopot

   !> Additional external potential that is added to the output, 
   !! when the XC parameter ixc/=0 and sumpion=.true.
   !! When sumpion=.true., it is always provided in the distributed form,
   !! clearly without the overlapping terms which are needed only for the XC part
   real(wp), dimension(*), intent(inout) :: pot_ion

   !> Optional argument to avoid output writings
   character(len=3), intent(in), optional :: quiet
   
   !> Additional external density added to the input, regardless of the value of sumpion,
   !! without the overlapping terms which are needed only for the XC part
   real(wp), dimension(kernel%grid%m1,kernel%grid%m3*kernel%grid%n3p), intent(inout), optional :: rho_ion


   !> Stress tensor: Add the stress tensor part from the Hartree potential
   real(dp), dimension(6), intent(out), optional :: stress_tensor

   !local variables
   character(len=*), parameter :: subname='H_potential'
   real(dp), parameter :: max_ratioex = 1.0e10_dp,eps0=78.36d0 !to be inserted in pkernel
   logical :: wrtmsg,cudasolver
   integer :: m1,m2,m3,md1,md2,md3,n1,n2,n3,nd1,nd2,nd3
   integer :: ierr,ind,ind2,ind3,indp,ind2p,ind3p,i
   integer :: i1,i2,i3,j2,istart,iend,i3start,jend,jproc,i3xcsh
   integer :: nxc,istden,istglo,ip,irho,i3s,i23,i23s,n23
   real(dp) :: ehartreeLOC,pot,rhores2,beta
   real(dp) :: alpha,ratio,normb,normr,norm_nonvac,e_static,rpoints
   !real(dp) :: scal
   real(dp), dimension(6) :: strten
   real(dp), dimension(:,:), allocatable :: rho,rhopol,x,q,p,r,z,depsdrho
   real(dp), dimension(:,:,:), allocatable :: zf,work_full
   !integer, dimension(:,:), allocatable :: gather_arr

   call f_routine(id='H_potential')
   
   cudasolver=.false.
   
   !do not write anything on screen if quiet is set to yes
   if (present(quiet)) then
      if(quiet == 'yes' .or. quiet == 'YES') then
         wrtmsg=.false.
      else if(trim(quiet) == 'no' .or. trim(quiet) == 'NO') then
         wrtmsg=.true.
      else
         call f_err_throw('ERROR: Unrecognised value for "quiet" option: ' // trim(quiet))
         !write(*,*)'ERROR: Unrecognised value for "quiet" option:',quiet
      end if
   else
      wrtmsg=.true.
   end if
   wrtmsg=wrtmsg .and. kernel%mpi_env%iproc==0 .and. kernel%mpi_env%igroup==0
   ! rewrite

   if (wrtmsg) call yaml_mapping_open('Poisson Solver')
   
   !call timing(kernel%mpi_env%iproc,'PSolv_comput  ','ON')
   call f_timing(TCAT_PSOLV_COMPUT,'ON')
   !calculate the dimensions wrt the geocode
   select case(kernel%geocode)
   case('P')
      if (wrtmsg) &
           call yaml_map('BC','Periodic')
      call P_FFT_dimensions(kernel%ndims(1),kernel%ndims(2),kernel%ndims(3),m1,m2,m3,n1,n2,n3,&
           md1,md2,md3,nd1,nd2,nd3,kernel%mpi_env%nproc,.false.)
   case('S')
      if (wrtmsg) &
           call yaml_map('BC','Surface')
      call S_FFT_dimensions(kernel%ndims(1),kernel%ndims(2),kernel%ndims(3),m1,m2,m3,n1,n2,n3,&
           md1,md2,md3,nd1,nd2,nd3,&
           kernel%mpi_env%nproc,kernel%igpu,.false.)
   case('F')
      if (wrtmsg) &
           call yaml_map('BC','Free')
      call F_FFT_dimensions(kernel%ndims(1),kernel%ndims(2),kernel%ndims(3),m1,m2,m3,n1,n2,n3,&
           md1,md2,md3,nd1,nd2,nd3,&
           kernel%mpi_env%nproc,kernel%igpu,.false.)
   case('W')
      if (wrtmsg) &
           call yaml_map('BC','Wires')
      call W_FFT_dimensions(kernel%ndims(1),kernel%ndims(2),kernel%ndims(3),m1,m2,m3,n1,n2,n3,&
           md1,md2,md3,nd1,nd2,nd3,kernel%mpi_env%nproc,kernel%igpu,.false.)
   case default
      call f_err_throw('geometry code ("'//kernel%geocode//'" not admitted in PSolver')
      !stop 'PSolver: geometry code not admitted'
   end select

   if (kernel%method /= PS_VAC_ENUM .and. .not. present(rho_ion) .and. sumpion) then
      call f_err_throw('Error in H_potential, for a cavity treatment the array rho_ion is needed')
   end if
   
   cudasolver= (kernel%igpu==1 .and. .not. present(stress_tensor))
   
   if (wrtmsg) then
      call yaml_map('Box',kernel%ndims,fmt='(i5)')
      call yaml_map('MPI tasks',kernel%mpi_env%nproc,fmt='(i5)')
      if (cudasolver) call yaml_map('GPU acceleration',.true.)
   end if
   
   !array allocations
   !initalise to zero the zf array
   zf = f_malloc([md1, md3, 2*md2/kernel%mpi_env%nproc],id='zf')

   select case(datacode)
   case('G')
      !starting address of rhopot in the case of global i/o
      i3start=kernel%grid%istart*kernel%ndims(1)*kernel%ndims(2)+1
      if (kernel%grid%n3p == 0) i3start=1
   case('D')
      !distributed i/o
      i3start=1
   case default
      call f_err_throw('datacode ("'//datacode//'") not admitted in PSolver')
   end select

   rpoints=product(real(kernel%ndims,dp))

   n23=kernel%grid%m3*kernel%grid%n3p
   n1=kernel%grid%m1

   !in the case of SC cavity, gather the full density and determine the depsdrho
   !here the if statement for the SC cavity should be put
   !print *,'method',trim(char(kernel%method)),associated(kernel%method%family),trim(char(kernel%method%family))
   if (kernel%method .hasattr. PS_SCCS_ENUM) then
      work_full=f_malloc(kernel%ndims,id='work_full')
      depsdrho=f_malloc([n1,n23],id='depsdrho')
      if (kernel%mpi_env%nproc > 1) then
         call mpiallgather(rhopot(i3start),recvbuf=work_full(1,1,1),recvcounts=kernel%counts,&
              displs=kernel%displs,comm=kernel%mpi_env%mpi_comm)
      else
         call f_memcpy(n=product(kernel%ndims),src=rhopot(i3start),dest=work_full(1,1,1))
      end if
      call pkernel_build_epsilon(kernel,work_full,eps0,depsdrho)
   end if

   !add the ionic density to the potential
   e_static=0.0_dp
   if (kernel%method /= PS_VAC_ENUM .and. kernel%grid%n3p>0 .and. sumpion) then
      !call yaml_map('Rho_ion monopole',sum(rho_ion)*product(kernel%hgrids))
      call finalize_hartree_results(.true.,rho_ion,kernel%grid%m1,kernel%grid%m3,kernel%grid%n3p,&
           kernel%grid%m1,kernel%grid%m3,kernel%grid%n3p,&
           pot_ion,rhopot(i3start),rhopot(i3start),e_static)
      if (wrtmsg) then
         call yaml_map('Rho_ion monopole',sum(rho_ion)*product(kernel%hgrids))
         call yaml_map('Summing up ionic density',.true.)
         call yaml_map('Estatic',e_static*product(kernel%hgrids))
         !call yaml_map('Total monopole',sum(rhopot(i3start:i3start+n1*n23-1)*product(kernel%hgrids)))
      end if
   end if

   !now switch the treatment according to the method used
   select case(trim(str(kernel%method)))
   case('VAC')
      !core psolver routine
      call apply_kernel(cudasolver,kernel,rhopot(i3start),offset,strten,zf,.false.)

      call finalize_hartree_results(sumpion,pot_ion,&
           kernel%grid%m1,kernel%grid%m3,kernel%grid%n3p,&
           kernel%grid%md1,kernel%grid%md3,2*(kernel%grid%md2/kernel%mpi_env%nproc),&
           rhopot(i3start),zf,rhopot(i3start),ehartreeLOC)
      !gathering the data to obtain the distribution array
      if (datacode == 'G' .and. kernel%mpi_env%nproc > 1) then
         call mpiallgather(rhopot(1),recvcounts=kernel%counts,&
              displs=kernel%displs,comm=kernel%mpi_env%mpi_comm)
      end if
   case('PI')
      if (.not. associated(kernel%oneoeps) .or. .not. associated(kernel%dlogeps))&
           call f_err_throw('The PI method needs the arrays of dlogeps and oneoeps,'//&
           ' use pkernel_set_epsilon routine to set them')
      
      if (datacode /= 'G' .and. kernel%mpi_env%nproc > 1) &
           call f_err_throw('Error in H_potential, PI method only works with datacode=G')
      !allocate arrays
      rhopol=f_malloc0([kernel%ndims(1),kernel%ndims(2)*kernel%ndims(3)],&
           id='rhopol')
      rho=f_malloc([n1,n23],id='rho')
      call f_memcpy(n=size(rho),src=rhopot(i3start),dest=rho)
      !check the integrity of rhopot with respect to the vacuum region
      call nonvacuum_projection(n1,n23,rho,kernel%oneoeps,norm_nonvac)
      norm_nonvac=norm_nonvac*product(kernel%hgrids)
      if (kernel%mpi_env%nproc > 1) &
           call mpiallred(norm_nonvac,1,MPI_SUM,comm=kernel%mpi_env%mpi_comm)
      if (wrtmsg) call yaml_map('Integral of the density in the nonvacuum region',norm_nonvac)

      if (wrtmsg) &
           call yaml_sequence_open('Embedded PSolver Polarization Iteration Method')
      pi_loop: do ip=1,kernel%max_iter

         !update the needed part of rhopot array
         irho=i3start
         !i3s=kernel%grid%istart
         i23s=kernel%grid%istart*kernel%grid%m3
         do i23=1,n23
            do i1=1,n1
               rhopot(irho)=&
                    kernel%oneoeps(i1,i23)*rho(i1,i23)+&
                    rhopol(i1,i23+i23s)
!!$                       kernel%oneoeps(i1,i2,i3+i3s)*rho(i1,i2,i3)+&
!!$                       rhopol(i1,i2,i3+i3s)
               irho=irho+1
            end do
         end do
         
         call apply_kernel(cudasolver,kernel,rhopot(i3start),offset,strten,zf,.true.)
         !gathering the data to obtain the distribution array
         !this method only works with datacode == 'G'
         if (kernel%mpi_env%nproc > 1) then
            call mpiallgather(rhopot(1),recvcounts=kernel%counts,&
                 displs=kernel%displs,comm=kernel%mpi_env%mpi_comm)
         end if
         
         !update rhopol and calculate residue
         call fssnord3DmatNabla_LG(kernel%geocode,kernel%ndims(1),kernel%ndims(2),&
              kernel%ndims(3),&
              rhopot,kernel%nord,kernel%hgrids,kernel%PI_eta,kernel%dlogeps,rhopol,rhores2)

         !reduction of the residue not necessary as datacode==G
         !if (kernel%mpi_env%nproc > 1) &
         !     call mpiallred(rhores2,1,MPI_SUM,comm=kernel%mpi_env%mpi_comm)
         rhores2=sqrt(rhores2/rpoints)

         if (wrtmsg) then 
            call yaml_newline()
            call yaml_sequence(advance='no')
            call EPS_iter_output(ip,0.0_dp,rhores2,0.0_dp,0.0_dp,0.0_dp)
         end if

         if (rhores2 < kernel%minres) exit pi_loop

      end do pi_loop
      if (wrtmsg) call yaml_sequence_close()

      !if statement for SC cavity
      if (kernel%method .hasattr. PS_SCCS_ENUM) &
           call extra_sccs_potential(kernel,work_full,depsdrho,rhopot(i3start))

      !here the harteee energy can be calculated and the ionic potential
      !added
      call finalize_hartree_results(sumpion,pot_ion,&
           kernel%grid%m1,kernel%grid%m3,kernel%grid%n3p,&
           kernel%grid%m1,kernel%grid%m3,kernel%grid%n3p,&
           rho,rhopot(i3start),rhopot(i3start),ehartreeLOC)
      
      if (kernel%method .hasattr. PS_SCCS_ENUM) then
           !then add the extra potential after having calculated the hartree energy
         call axpy(n1*n23,1.0_dp,depsdrho(1,1),1,rhopot(i3start),1)
      end if

      if (kernel%mpi_env%nproc > 1 .and. sumpion) then
         call mpiallgather(rhopot(1),recvcounts=kernel%counts,&
              displs=kernel%displs,comm=kernel%mpi_env%mpi_comm)
      end if

      call f_free(rho)
      call f_free(rhopol)

   case('PCG')
      if (.not. associated(kernel%oneoeps) .or. .not. associated(kernel%corr))&
           call f_err_throw('The PCG method needs the arrays corr and oneosqrteps,'//&
           ' use pkernel_set_epsilon routine to set them')

      call nonvacuum_projection(n1,n23,rhopot(i3start),kernel%oneoeps,norm_nonvac)
      norm_nonvac=norm_nonvac*product(kernel%hgrids)
      if (kernel%mpi_env%nproc > 1) &
           call mpiallred(norm_nonvac,1,MPI_SUM,comm=kernel%mpi_env%mpi_comm)
      if (wrtmsg) call yaml_map('Integral of the density in the nonvacuum region',norm_nonvac)


      r=f_malloc([n1,n23],id='r')
      x=f_malloc0([n1,n23],id='x')
      q=f_malloc0([n1,n23],id='q')
      p=f_malloc0([n1,n23],id='p')
      z=f_malloc([n1,n23],id='z')

      if (wrtmsg) &
           call yaml_sequence_open('Embedded PSolver, Preconditioned Conjugate Gradient Method')
      beta=1.d0
      ratio=1.d0
      normr=1.d0

      normb=dot(n1*n23,rhopot(i3start),1,rhopot(i3start),1)
      if (kernel%mpi_env%nproc > 1) &
           call mpiallred(normb,1,MPI_SUM,comm=kernel%mpi_env%mpi_comm)
      normb=sqrt(normb/rpoints)

      call f_memcpy(n=kernel%grid%m1*kernel%grid%m3*kernel%grid%n3p,&
           src=rhopot(i3start),dest=r)

      !$omp parallel do default(shared) private(i1,i23)
      do i23=1,n23
         do i1=1,n1
             z(i1,i23)=r(i1,i23)*kernel%oneoeps(i1,i23)
         end do
      end do
      !$omp end parallel do

      PCG_loop: do ip=1,kernel%max_iter

         if (normr < kernel%minres .or. normr > max_ratioex) exit PCG_loop

         !  Apply the Preconditioner
         call apply_kernel(cudasolver,kernel,z,offset,strten,zf,.true.)

         call apply_reductions(ip, cudasolver, kernel, r, x, p, q, z, alpha, beta, normr)

         normr=sqrt(normr/rpoints)

         ratio=normr/normb
         if (wrtmsg) then
            call yaml_newline()
            call yaml_sequence(advance='no')
            call EPS_iter_output(ip,normb,normr,ratio,alpha,beta)
         end if
      end do PCG_loop
      if (wrtmsg) call yaml_sequence_close()
      
      !only useful for gpu, bring back the x array
      call update_pot_from_device(cudasolver, kernel, x)
      !if statement for SC cavity
      if (kernel%method .hasattr. PS_SCCS_ENUM)&
           call extra_sccs_potential(kernel,work_full,depsdrho,x)

      !here the harteee energy can be calculated and the ionic potential
      !added
      call finalize_hartree_results(sumpion,pot_ion,&
           kernel%grid%m1,kernel%grid%m3,kernel%grid%n3p,&
           kernel%grid%m1,kernel%grid%m3,kernel%grid%n3p,&
           rhopot(i3start),x,rhopot(i3start),ehartreeLOC)

      if (kernel%method .hasattr. PS_SCCS_ENUM) then
         !then add the extra potential after having calculated the hartree energy
         call axpy(n1*n23,1.0_dp,depsdrho(1,1),1,rhopot(i3start),1)
      end if

      if (kernel%mpi_env%nproc > 1 .and. datacode =='G') then
         call mpiallgather(rhopot(1),recvcounts=kernel%counts,&
              displs=kernel%displs,comm=kernel%mpi_env%mpi_comm)
      end if
      call f_free(x)
      call f_free(z)
      call f_free(p)
      call f_free(q)
      call f_free(r)

   end select

   call f_free(zf)

   !if statement for SC cavity to be added
   if (kernel%method .hasattr. PS_SCCS_ENUM) then
      call f_free(work_full)
      call f_free(depsdrho)
   end if

   !check for the presence of the stress tensor
   if (present(stress_tensor)) call f_memcpy(src=strten,dest=stress_tensor)

   !evaluating the total ehartree + e_static if needed
   eh=(ehartreeLOC*0.5_dp+e_static)*product(kernel%hgrids)
   if (kernel%mpi_env%nproc > 1) then

      call mpiallred(eh,1,MPI_SUM,comm=kernel%mpi_env%mpi_comm)
      !reduce also the value of the stress tensor
   
      if (present(stress_tensor)) then
         call mpiallred(stress_tensor,MPI_SUM,comm=kernel%mpi_env%mpi_comm)
      end if
   end if
   
      !call timing(kernel%mpi_env%iproc,'PSolv_commun  ','OF')
      !call f_timing(TCAT_PSOLV_COMPUT,'OF')

      if (datacode == 'G') then
         !building the array of the data to be sent from each process
         !and the array of the displacement
   
!!$         gather_arr = f_malloc((/ 0.to.kernel%mpi_env%nproc-1, 1.to.2 /),id='gather_arr')
!!$         do jproc=0,kernel%mpi_env%nproc-1
!!$            istart=min(jproc*(md2/kernel%mpi_env%nproc),m2-1)
!!$            jend=max(min(md2/kernel%mpi_env%nproc,m2-md2/kernel%mpi_env%nproc*jproc),0)
!!$            gather_arr(jproc,1)=m1*m3*jend
!!$            gather_arr(jproc,2)=m1*m3*istart
!!$         end do
!!$         !gather all the results in the same rhopot array
!!$         istart=min(kernel%mpi_env%iproc*(md2/kernel%mpi_env%nproc),m2-1)
!!$   
!!$         istden=1+kernel%ndims(1)*kernel%ndims(2)*istart
!!$         istglo=1
!!$         !call mpiallgatherv(rhopot(istglo), gather_arr(:,1), gather_arr(:,2), &
         !call mpiallgatherv(rhopot, kernel%counts,kernel%displs, &
         !     kernel%mpi_env%iproc, kernel%mpi_env%mpi_comm)
!!$         call f_free(gather_arr)
         !call timing(kernel%mpi_env%iproc,'PSolv_comput  ','OF')
      end if
!!$   end if

   if (wrtmsg) call yaml_mapping_close()
   call f_timing(TCAT_PSOLV_COMPUT,'OF')
   call f_release_routine()

END SUBROUTINE H_potential

subroutine extra_sccs_potential(kernel,work_full,depsdrho,pot)
  implicit none
  type(coulomb_operator), intent(in) :: kernel
  real(dp), dimension(kernel%ndims(1),kernel%ndims(2),kernel%ndims(3)), intent(out) :: work_full
  real(dp), dimension(kernel%ndims(1),kernel%ndims(2)*kernel%grid%n3p), intent(inout) :: depsdrho
  real(dp), dimension(kernel%ndims(1),kernel%ndims(2)*kernel%grid%n3p) :: pot !intent in

  !first gather the potential to calculate the derivative
  if (kernel%mpi_env%nproc > 1) then
     call mpiallgather(pot,recvbuf=work_full,recvcounts=kernel%counts,&
          displs=kernel%displs,comm=kernel%mpi_env%mpi_comm)
  else
     call f_memcpy(src=pot,dest=work_full)
  end if

  !then calculate the extra potential and add it to pot
  call sccs_extra_potential(kernel,work_full,depsdrho)
  
end subroutine extra_sccs_potential

!> verify that the density is considerably zero in the region where epsilon is different from one
subroutine nonvacuum_projection(n1,n23,rho,oneoeps,norm)
  implicit none
  integer, intent(in) :: n1,n23 !< parallelized box dimensions
  real(dp), dimension(n1,n23), intent(in) :: rho !<charge density
  !>inverse of epsilon (might also be inverse of sqrt(eps))
  real(dp), dimension(n1,n23), intent(in) :: oneoeps 
  real(dp), intent(out) :: norm !< \int of rho where epsilon /=1
  !local variables
  integer :: i1,i23
  real(dp), parameter :: tol= 5.d-1
 
  norm=0.0_dp
  !$omp parallel do default(shared) private(i1,i23)&
  !$omp reduction(+:norm)
  do i23=1,n23
     do i1=1,n1
        if (abs(oneoeps(i1,i23) - 1.0_dp) > tol) norm=norm+rho(i1,i23)
     end do
  end do
  !$omp end parallel do

end subroutine nonvacuum_projection

!regroup the psolver from here -------------
subroutine apply_kernel(gpu,kernel,rho,offset,strten,zf,updaterho)
  use f_utils, only: f_zero
  use time_profiling, only: f_timing
  implicit none
  !>when true, the density is updated with the value of zf
  logical, intent(in) :: updaterho
  logical, intent(in) :: gpu !< logical variable controlling the gpu acceleration
  type(coulomb_operator), intent(in) :: kernel 
  !> Total integral on the supercell of the final potential on output
  real(dp), intent(in) :: offset
  real(dp), dimension(kernel%grid%m1,kernel%grid%m3*kernel%grid%n3p), intent(inout) :: rho
  !>work array for the usage of the main routine
  real(dp), dimension(kernel%grid%md1,kernel%grid%md3*2*(kernel%grid%md2/kernel%mpi_env%nproc)), intent(inout) :: zf
  real(dp), dimension(6), intent(out) :: strten !< stress tensor associated to the cell
  !local variables
  integer, dimension(3) :: n
  integer :: size1,size2,switch_alg,i_stat,ierr,i23,j23,j3,i1,n3delta
  real(dp) :: pt,rh
  real(dp), dimension(:), allocatable :: zf1

  call f_zero(zf)
  !this routine builds the values for each process of the potential (zf), multiplying by scal   
  !fill the array with the values of the charge density
  !no more overlap between planes
  !still the complex case should be defined
  n3delta=kernel%grid%md3-kernel%grid%m3
  !$omp parallel do default(shared) private(i1,i23,j23,j3)
  do i23=1,kernel%grid%n3p*kernel%grid%m3
     j3=(i23-1)/kernel%grid%m3
     !j2=i23-kernel%grid%m3*j3
     j23=i23+n3delta*j3!=j2+kernel%grid%md3*j3
     do i1=1,kernel%grid%m1
        zf(i1,j23)=rho(i1,i23)
     end do
  end do
  !$omp end parallel do

  call f_zero(strten)
  if (.not. gpu) then !CPU case
     call f_timing(TCAT_PSOLV_COMPUT,'OF')
     call G_PoissonSolver(kernel%mpi_env%iproc,kernel%mpi_env%nproc,&
          kernel%part_mpi%mpi_comm,kernel%inplane_mpi%iproc,&
          kernel%inplane_mpi%mpi_comm,kernel%geocode,1,&
          kernel%grid%n1,kernel%grid%n2,kernel%grid%n3,&
          kernel%grid%nd1,kernel%grid%nd2,kernel%grid%nd3,&
          kernel%grid%md1,kernel%grid%md2,kernel%grid%md3,&
          kernel%kernel,zf,&
          kernel%grid%scal,kernel%hgrids(1),kernel%hgrids(2),kernel%hgrids(3),offset,strten)
     call f_timing(TCAT_PSOLV_COMPUT,'ON')

  else !GPU case

     n(1)=kernel%grid%n1!kernel%ndims(1)*(2-kernel%geo(1))
     n(2)=kernel%grid%n3!kernel%ndims(2)*(2-kernel%geo(2))
     n(3)=kernel%grid%n2!kernel%ndims(3)*(2-kernel%geo(3))

     size1=kernel%grid%md1*kernel%grid%md2*kernel%grid%md3! nproc always 1 kernel%ndims(1)*kernel%ndims(2)*kernel%ndims(3)

     if (kernel%keepGPUmemory == 0) then
        size2=2*kernel%grid%n1*kernel%grid%n2*kernel%grid%n3
        call cudamalloc(size2,kernel%work1_GPU,i_stat)
        if (i_stat /= 0) print *,'error cudamalloc',i_stat
        call cudamalloc(size2,kernel%work2_GPU,i_stat)
        if (i_stat /= 0) print *,'error cudamalloc',i_stat
     endif

     if (kernel%mpi_env%nproc > 1) then
        zf1 = f_malloc(size1,id='zf1')

        call mpi_gather(zf,size1/kernel%mpi_env%nproc,mpidtypd,zf1,size1/kernel%mpi_env%nproc, &
             mpidtypd,0,kernel%mpi_env%mpi_comm,ierr)

        if (kernel%mpi_env%iproc == 0) then
           !fill the GPU memory
           call reset_gpu_data(size1,zf1,kernel%work1_GPU)

           switch_alg=0
           if (kernel%initCufftPlan == 1) then
              call cuda_3d_psolver_general(n,kernel%plan,kernel%work1_GPU,kernel%work2_GPU, &
                   kernel%k_GPU,switch_alg,kernel%geo,kernel%grid%scal)
           else
              call cuda_3d_psolver_plangeneral(n,kernel%work1_GPU,kernel%work2_GPU, &
                   kernel%k_GPU,kernel%geo,kernel%grid%scal)
           endif

           !take data from GPU
           call get_gpu_data(size1,zf1,kernel%work1_GPU)
        endif

        call MPI_Scatter(zf1,size1/kernel%mpi_env%nproc,mpidtypd,zf,size1/kernel%mpi_env%nproc, &
             mpidtypd,0,kernel%mpi_env%mpi_comm,ierr)

        call f_free(zf1)

     else

        !fill the GPU memory
        call reset_gpu_data(size1,zf,kernel%work1_GPU)

        switch_alg=0

        if (kernel%initCufftPlan == 1) then
           call cuda_3d_psolver_general(n,kernel%plan,kernel%work1_GPU,kernel%work2_GPU, &
                kernel%k_GPU,switch_alg,kernel%geo,kernel%grid%scal)
        else
           call cuda_3d_psolver_plangeneral(n,kernel%work1_GPU,kernel%work2_GPU, &
                kernel%k_GPU,kernel%geo,kernel%grid%scal)
        endif


        !take data from GPU
        call get_gpu_data(size1,zf,kernel%work1_GPU)


     endif

     if (kernel%keepGPUmemory == 0) then
        call cudafree(kernel%work1_GPU)
        call cudafree(kernel%work2_GPU)
     endif

  endif

  if (updaterho) then
     !$omp parallel do default(shared) private(i1,i23,j23,j3)
     do i23=1,kernel%grid%n3p*kernel%grid%m3
        j3=(i23-1)/kernel%grid%m3
        j23=i23+n3delta*j3
        do i1=1,kernel%grid%m1
           rho(i1,i23)=zf(i1,j23)
        end do
     end do
     !$omp end parallel do
  end if

end subroutine apply_kernel



subroutine apply_reductions(ip, gpu, kernel, r, x, p, q, z, alpha, beta, normr)
  use f_utils, only: f_zero
  use time_profiling, only: f_timing
  use yaml_output
  implicit none
  logical, intent(in) :: gpu !< logical variable controlling the gpu acceleration
  type(coulomb_operator), intent(in) :: kernel 
  real(dp), dimension(kernel%grid%m1,kernel%grid%m3*kernel%grid%n3p), intent(inout) :: r,x,p,q,z
  !local variables
  integer :: n1,n23,i_stat,ierr,i23,i1,size1, ip
  real(dp) :: beta, zeta, rval, alpha, normr, beta0, epsc, kappa, pval, qval



  n23=kernel%grid%m3*kernel%grid%n3p
  n1=kernel%grid%m1
  beta0 = beta
  beta=0.d0
!  if (.true.) then !CPU case
  if (.not. gpu) then !CPU case
         !$omp parallel do default(shared) private(i1,i23,rval,zeta) &
         !$omp reduction(+:beta)
         do i23=1,n23
            do i1=1,n1
               zeta=z(i1,i23)
               zeta=zeta*kernel%oneoeps(i1,i23)
               rval=r(i1,i23)
               rval=rval*zeta
               beta=beta+rval
               z(i1,i23)=zeta
            end do
         end do
         !$omp end parallel do


         if (kernel%mpi_env%nproc > 1) &
              call mpiallred(beta,1,MPI_SUM,comm=kernel%mpi_env%mpi_comm)
         kappa=0.d0
         !$omp parallel do default(shared) private(i1,i23,epsc,zeta)&
         !$omp private(pval,qval,rval) reduction(+:kappa)
         do i23=1,n23
            do i1=1,n1
               zeta=z(i1,i23)
               epsc=kernel%corr(i1,i23)
               pval=p(i1,i23)
               qval=q(i1,i23)
               rval=r(i1,i23)
               pval = zeta+(beta/beta0)*pval
               qval = zeta*epsc+rval+(beta/beta0)*qval
               p(i1,i23) = pval
               q(i1,i23) = qval
               rval=pval*qval
               kappa = kappa+rval 
            end do
         end do
         !$omp end parallel do
         if (kernel%mpi_env%nproc > 1) &
              call mpiallred(kappa,1,MPI_SUM,comm=kernel%mpi_env%mpi_comm)
         alpha = beta/kappa

         normr=0.d0
         !$omp parallel do default(shared) private(i1,i23,rval) &
         !$omp reduction(+:normr)
         do i23=1,n23
            do i1=1,n1
               x(i1,i23) = x(i1,i23) + alpha*p(i1,i23)
               r(i1,i23) = r(i1,i23) - alpha*q(i1,i23)
               z(i1,i23) = r(i1,i23)*kernel%oneoeps(i1,i23)
               rval=r(i1,i23)*r(i1,i23)
               normr=normr+rval
            end do
         end do
         !$omp end parallel do
         if (kernel%mpi_env%nproc > 1) &
              call mpiallred(normr,1,MPI_SUM,comm=kernel%mpi_env%mpi_comm)

  else
!naive method, with allocations/free at each time .. 
!may need to store more pointers inside kernel
  size1=n1*n23
  if (kernel%keepGPUmemory == 0) then
    call cudamalloc(size1,kernel%z_GPU,i_stat)
    if (i_stat /= 0) print *,'error cudamalloc',i_stat
    call cudamalloc(size1,kernel%r_GPU,i_stat)
    if (i_stat /= 0) print *,'error cudamalloc',i_stat
    call cudamalloc(size1,kernel%oneoeps_GPU,i_stat)
    if (i_stat /= 0) print *,'error cudamalloc',i_stat
    call cudamalloc(size1,kernel%p_GPU,i_stat)
    if (i_stat /= 0) print *,'error cudamalloc',i_stat
    call cudamalloc(size1,kernel%q_GPU,i_stat)
    if (i_stat /= 0) print *,'error cudamalloc',i_stat
    call cudamalloc(size1,kernel%x_GPU,i_stat)
    if (i_stat /= 0) print *,'error cudamalloc',i_stat
    call cudamalloc(size1,kernel%corr_GPU,i_stat)
    if (i_stat /= 0) print *,'error cudamalloc',i_stat


    call cudamalloc(sizeof(alpha),kernel%alpha_GPU,i_stat)
    if (i_stat /= 0) print *,'error cudamalloc',i_stat
    call cudamalloc(sizeof(beta),kernel%beta_GPU,i_stat)
    if (i_stat /= 0) print *,'error cudamalloc',i_stat
    call cudamalloc(sizeof(beta0),kernel%beta0_GPU,i_stat)
    if (i_stat /= 0) print *,'error cudamalloc',i_stat
    call cudamalloc(sizeof(kappa),kernel%kappa_GPU,i_stat)
    if (i_stat /= 0) print *,'error cudamalloc',i_stat
  end if

  if (ip == 1) then 
    call reset_gpu_data(size1,r,kernel%r_GPU)
    call reset_gpu_data(size1,kernel%oneoeps,kernel%oneoeps_GPU)
    call cudamemset(kernel%p_GPU, 0, size1,i_stat)
    if (i_stat /= 0) print *,'error cudamemset p',i_stat
    call cudamemset(kernel%q_GPU, 0, size1,i_stat)
    if (i_stat /= 0) print *,'error cudamemset q',i_stat
    call cudamemset(kernel%x_GPU, 0, size1,i_stat)
    if (i_stat /= 0) print *,'error cudamemset x',i_stat
  end if

  call reset_gpu_data(size1,z,kernel%z_GPU)
  call cudamemset(kernel%beta_GPU, 0, sizeof(beta),i_stat)
  call first_reduction_kernel(n1,n23,kernel%p_GPU,kernel%q_GPU,kernel%r_GPU,&
kernel%x_GPU,kernel%z_GPU,kernel%corr_GPU, kernel%oneoeps_GPU, kernel%alpha_GPU,&
 kernel%beta_GPU, kernel%beta0_GPU, kernel%kappa_GPU, beta)

  call cudamemset(kernel%kappa_GPU, 0, sizeof(kappa),i_stat)

!TODO : gpudirect.
         if (kernel%mpi_env%nproc > 1) &
              call mpiallred(beta,1,MPI_SUM,comm=kernel%mpi_env%mpi_comm)
         kappa=0.d0

  call reset_gpu_data(sizeof(beta),beta, kernel%beta_GPU)
  call reset_gpu_data(sizeof(beta0),beta0, kernel%beta0_GPU)

  if (ip == 1) then 
    call reset_gpu_data(size1,kernel%corr,kernel%corr_GPU)
  end if

  call second_reduction_kernel(n1,n23,kernel%p_GPU,kernel%q_GPU,kernel%r_GPU,&
kernel%x_GPU,kernel%z_GPU,kernel%corr_GPU, kernel%oneoeps_GPU, kernel%alpha_GPU,&
 kernel%beta_GPU, kernel%beta0_GPU, kernel%kappa_GPU, kappa)

!  call get_gpu_data(size1,p,kernel%p_GPU)
!  call get_gpu_data(size1,q,kernel%q_GPU)

  if (kernel%mpi_env%nproc > 1) &
    call mpiallred(kappa,1,MPI_SUM,comm=kernel%mpi_env%mpi_comm)

  alpha = beta/kappa
  call reset_gpu_data(sizeof(alpha),alpha,kernel%alpha_GPU)
  normr=0.d0

  call third_reduction_kernel(n1,n23,kernel%p_GPU,kernel%q_GPU,kernel%r_GPU,&
kernel%x_GPU,kernel%z_GPU,kernel%corr_GPU, kernel%oneoeps_GPU, kernel%alpha_GPU,&
 kernel%beta_GPU, kernel%beta0_GPU, kernel%kappa_GPU, normr)

         if (kernel%mpi_env%nproc > 1) &
              call mpiallred(normr,1,MPI_SUM,comm=kernel%mpi_env%mpi_comm)
  
  call get_gpu_data(size1,z,kernel%z_GPU)


  if (kernel%keepGPUmemory == 0) then
    call cudafree(kernel%z_GPU)
    call cudafree(kernel%r_GPU)
    call cudafree(kernel%oneoeps_GPU)
    call cudafree(kernel%p_GPU)
    call cudafree(kernel%q_GPU)
    call cudafree(kernel%x_GPU)
    call cudafree(kernel%corr_GPU)
    call cudafree(kernel%alpha_GPU)
    call cudafree(kernel%beta_GPU)
    call cudafree(kernel%beta0_GPU)
    call cudafree(kernel%kappa_GPU)
  end if

  end if

end subroutine apply_reductions

!at the end of the loop, we have to synchronize potential from gpu
subroutine update_pot_from_device(gpu, kernel,x)
  type(coulomb_operator), intent(in) :: kernel 
  real(dp), dimension(kernel%grid%m1,kernel%grid%m3*kernel%grid%n3p), intent(inout) :: x
  logical, intent(in) :: gpu !< logical variable controlling the gpu acceleration
  integer size1
  if (gpu) then !CPU case
    size1=kernel%grid%m3*kernel%grid%n3p*kernel%grid%m1
    call get_gpu_data(size1,x,kernel%x_GPU)
  end if 
end subroutine update_pot_from_device

!> calculate the integral between the array in zf and the array in rho
!! copy the zf array in pot and sum with the array pot_ion if needed
subroutine finalize_hartree_results(sumpion,pot_ion,m1,m2,m3p,md1,md2,md3p,&
     rho,zf,pot,eh)
  implicit none
  !if .true. the array pot is zf+pot_ion
  !if .false., pot is only zf
  logical, intent(in) :: sumpion
  integer, intent(in) :: m1,m2,m3p !< dimension of the grid
  integer, intent(in) :: md1,md2,md3p !< dimension of the zf array
  !> original density and final potential (can point to the same array)
  real(dp), dimension(m1,m2*m3p) :: rho,pot 
  !> ionic potential, to be added to the potential
  real(dp), dimension(m1,m2*m3p), intent(in) :: pot_ion
  !> work array for the usage of the main routine (might have the same dimension of rho)
  real(dp), dimension(md1,md2*md3p), intent(inout) :: zf
  !>hartree energy, being \int d^3 x \rho(x) V_H(x)
  real(dp), intent(out) :: eh
  !local variables
  integer :: i1,i23,j23,j3,n3delta
  real(dp) :: pt,rh

  eh=0.0_dp
  !recollect the final data
  n3delta=md2-m2 !this is the y dimension
  if (sumpion) then
     !$omp parallel do default(shared) private(i1,i23,j23,j3,pt,rh) &
     !$omp reduction(+:eh)
     do i23=1,m2*m3p
        j3=(i23-1)/m2
        j23=i23+n3delta*j3
        do i1=1,m1
           pt=zf(i1,j23)
           rh=rho(i1,i23)
           rh=rh*pt
           eh=eh+rh
           pot(i1,i23)=pt+pot_ion(i1,i23)
        end do
     end do
     !$omp end parallel do
  else
     !$omp parallel do default(shared) private(i1,i23,j23,j3,pt,rh) &
     !$omp reduction(+:eh)
     do i23=1,m2*m3p
        j3=(i23-1)/m2
        j23=i23+n3delta*j3
        do i1=1,m1
           pt=zf(i1,j23)
           rh=rho(i1,i23)
           rh=rh*pt      
           eh=eh+rh
           pot(i1,i23)=pt
        end do
     end do
     !$omp end parallel do
  end if
end subroutine finalize_hartree_results

subroutine EPS_iter_output(iter,normb,normr,ratio,alpha,beta)
  implicit none
  integer, intent(in) :: iter
  real(dp), intent(in) :: normb,normr,ratio,beta,alpha
  !local variables
  character(len=*), parameter :: vrb='(1pe25.17)'!'(1pe16.4)'

  call yaml_mapping_open('Iteration quality',flow=.true.)
  if (beta /= 0.0_dp) call yaml_comment('Iteration '//iter,hfill='_')
  !write the PCG iteration
  call yaml_map('iter',iter,fmt='(i4)')
  !call yaml_map('rho_norm',normb)
  if (normr/=0.0_dp) call yaml_map('res',normr,fmt=vrb)
  if (ratio /= 0.0_dp) call yaml_map('ratio',ratio,fmt=vrb)
  if (alpha /= 0.0_dp) call yaml_map('alpha',alpha,fmt=vrb)
  if (beta /= 0.0_dp) call yaml_map('beta',beta,fmt=vrb)

  call yaml_mapping_close()
end subroutine EPS_iter_output


!> Calculate the dimensions needed for the allocation of the arrays 
!! related to the Poisson Solver
!!
!! @warning
!!    The XC enlarging due to GGA part is not present for surfaces and 
!!    periodic boundary condition. This is related to the fact that the calculation of the
!!    gradient and the White-Bird correction are not yet implemented for non-isolated systems
!! @author Luigi Genovese
!! @date February 2007
subroutine PS_dim4allocation(geocode,datacode,iproc,nproc,n01,n02,n03,use_gradient,use_wb_corr,&
      n3d,n3p,n3pi,i3xcsh,i3s)
   implicit none
   character(len=1), intent(in) :: geocode  !< @copydoc poisson_solver::doc::geocode
   character(len=1), intent(in) :: datacode !< @copydoc poisson_solver::doc::datacode
   integer, intent(in) :: iproc        !< Process Id
   integer, intent(in) :: nproc        !< Number of processes
   integer, intent(in) :: n01,n02,n03  !< Dimensions of the real space grid to be hit with the Poisson Solver
   logical, intent(in) :: use_gradient !< .true. if functional is using the gradient.
   logical, intent(in) :: use_wb_corr  !< .true. if functional is using WB corrections.
   !> Third dimension of the density. For distributed data, it takes into account 
   !! the enlarging needed for calculating the XC functionals.
   !! For global data it is simply equal to n03. 
   !! When there are too many processes and there is no room for the density n3d=0
   integer, intent(out) :: n3d
   !> Third dimension for the potential. The same as n3d, but without 
   !! taking into account the enlargment for the XC part. For non-GGA XC, n3p=n3d.
   integer, intent(out) :: n3p
   !> Dimension of the pot_ion array, always with distributed data. 
   !! For distributed data n3pi=n3p
   integer, intent(out) :: n3pi
   !> Shift of the density that must be performed to enter in the 
   !! non-overlapping region. Useful for recovering the values of the potential
   !! when using GGA XC functionals. If the density starts from rhopot(1,1,1),
   !! the potential starts from rhopot(1,1,i3xcsh+1). 
   !! For non-GGA XCs and for global distribution data i3xcsh=0
   integer, intent(out) :: i3xcsh
   !> Starting point of the density effectively treated by each processor 
   !! in the third direction.
   !! It takes into account also the XC enlarging. The array rhopot will correspond
   !! To the planes of third coordinate from i3s to i3s+n3d-1. 
   !! The potential to the planes from i3s+i3xcsh to i3s+i3xcsh+n3p-1
   !! The array pot_ion to the planes from i3s+i3xcsh to i3s+i3xcsh+n3pi-1
   !! For global disposition i3s is equal to distributed case with i3xcsh=0.
   integer, intent(out) :: i3s
   !local variables
   !n(c) integer, parameter :: nordgr=4
   integer :: m1,m2,m3,md1,md2,md3,n1,n2,n3,nd1,nd2,nd3
   integer :: istart,iend,nxc,nwb,nxt,nxcl,nxcr,nwbl,nwbr
   integer :: n3pr1,n3pr2
   
   
   !calculate the dimensions wrt the geocode
   select case(geocode)
   case('P')
      call P_FFT_dimensions(n01,n02,n03,m1,m2,m3,n1,n2,n3,md1,md2,md3,nd1,nd2,nd3,nproc,.false.)
       if (nproc>2*(n3/2+1)-1 .and. .false.) then
        n3pr1=nproc/(n3/2+1)
        n3pr2=n3/2+1
        if ((md2/nproc)*n3pr1*n3pr2 < n2) then
           md2=(md2/nproc+1)*nproc
        endif
      endif
   case('S')
      call S_FFT_dimensions(n01,n02,n03,m1,m2,m3,n1,n2,n3,md1,md2,md3,nd1,nd2,nd3,nproc,0,.false.)
      if (nproc>2*(n3/2+1)-1 .and. .false.) then
        n3pr1=nproc/(n3/2+1)
        n3pr2=n3/2+1
        if ((md2/nproc)*n3pr1*n3pr2 < n2) then
           md2=(md2/nproc+1)*nproc
        endif
      endif
   case('F','H')
      call F_FFT_dimensions(n01,n02,n03,m1,m2,m3,n1,n2,n3,md1,md2,md3,nd1,nd2,nd3,nproc,0,.false.)
      if (nproc>2*(n3/2+1)-1 .and. .false.) then
        n3pr1=nproc/(n3/2+1)
        n3pr2=n3/2+1
        if ((md2/nproc)*n3pr1*n3pr2 < n2/2) then
           md2=(md2/nproc+1)*nproc 
        endif
      endif
   case('W')
      call W_FFT_dimensions(n01,n02,n03,m1,m2,m3,n1,n2,n3,md1,md2,md3,nd1,nd2,nd3,nproc,0,.false.)
      if (nproc>2*(n3/2+1)-1 .and. .false.) then
        n3pr1=nproc/(n3/2+1)
        n3pr2=n3/2+1
        if ((md2/nproc)*n3pr1*n3pr2 < n2) then
           md2=(md2/nproc+1)*nproc
        endif
      endif
   case default
      !write(*,*) geocode
      !stop 'PS_dim4allocation: geometry code not admitted'
      call f_err_throw('Geometry code "'//geocode//'" not admitted in PS_dim4allocation')
   end select
   
   !formal start and end of the slice
   istart=iproc*(md2/nproc)
   iend=min((iproc+1)*md2/nproc,m2)
   
   select case(datacode)
   case('D')
      call xc_dimensions(geocode,use_gradient,use_wb_corr,istart,iend,m2,nxc,nxcl,nxcr,nwbl,nwbr,i3s,i3xcsh)
   
      nwb=nxcl+nxc+nxcr-2
      nxt=nwbr+nwb+nwbl
   
      n3p=nxc
      n3d=nxt
      n3pi=n3p
   case('G')
      n3d=n03
      n3p=n03
      i3xcsh=0
      i3s=min(istart,m2-1)+1
      n3pi=max(iend-istart,0)
   case default
      !print *,datacode
      !stop 'PS_dim4allocation: data code not admitted'
      call f_err_throw('data code "'//datacode//'" not admitted in PS_dim4allocation')
   end select

!!  print *,'P4,iproc',iproc,'nxc,ncxl,ncxr,nwbl,nwbr',nxc,nxcl,nxcr,nwbl,nwbr,&
!!       'ixc,n3d,n3p,i3xcsh,i3s',ixc,n3d,n3p,i3xcsh,i3s

END SUBROUTINE PS_dim4allocation


!> Calculate the dimensions to be used for the XC part, taking into account also
!! the White-bird correction which should be made for some GGA functionals
!!
!! SYNOPSIS
!!    @param use_gradient .true. if functional is using the gradient.
!!    @param use_wb_corr  .true. if functional is using WB corrections.
!!    @param m2        dimension to be parallelised
!!    @param nxc       size of the parallelised XC potential
!!    @param ncxl,ncxr left and right buffers for calculating the WB correction after call drivexc
!!    @param nwbl,nwbr left and right buffers for calculating the gradient to pass to drivexc    
!!    @param i3s       starting addres of the distributed dimension
!!    @param i3xcsh    shift to be applied to i3s for having the striting address of the potential
!!
!! @warning It is imperative that iend <=m2
subroutine xc_dimensions(geocode,use_gradient,use_wb_corr,&
     & istart,iend,m2,nxc,nxcl,nxcr,nwbl,nwbr,i3s,i3xcsh)
  implicit none

  character(len=1), intent(in) :: geocode !< @copydoc poisson_solver::doc::geocode
  logical, intent(in) :: use_gradient, use_wb_corr
  integer, intent(in) :: istart,iend,m2
  integer, intent(out) :: nxc,nxcl,nxcr,nwbl,nwbr,i3s,i3xcsh
  !local variables
  integer, parameter :: nordgr=4 !the order of the finite-difference gradient (fixed)

  if (istart <= m2-1) then
     nxc=iend-istart
     if (use_gradient .and. geocode == 'F') then
        if (.not. use_wb_corr) then
           !now the dimension of the part required for the gradient
           nwbl=min(istart,nordgr)
           nwbr=min(m2-iend,nordgr) !always m2 < iend
           nxcl=1
           nxcr=1
        else
           !now the dimension of the part required for the gradient
           if(istart<=nordgr) then
              nxcl=istart+1
              nwbl=0
           else
              nxcl=nordgr+1
              nwbl=min(nordgr,istart-nordgr)
           end if
           if(iend>=m2-nordgr+1) then
              nxcr=m2-iend+1
              nwbr=0
           else
              nxcr=nordgr+1
              nwbr=min(nordgr,m2-nordgr-iend)
           end if
        end if
     else if (geocode /= 'F' .and. use_gradient .and. nxc /= m2) then
        if (.not. use_wb_corr) then
           !now the dimension of the part required for the gradient
           nwbl=nordgr
           nwbr=nordgr
           nxcl=1
           nxcr=1
        else
           nxcl=nordgr+1
           nwbl=nordgr
           nxcr=nordgr+1
           nwbr=nordgr
        end if
     !this case is also considered below
     !else if (geocode /= 'F' .and. use_gradient .and. nxc == m2) then
     else 
        nwbl=0
        nwbr=0
        nxcl=1
        nxcr=1
     end if
     i3xcsh=nxcl+nwbl-1
     i3s=istart+1-i3xcsh
  else
     nwbl=0
     nwbr=0
     nxcl=1
     nxcr=1
     nxc=0
     i3xcsh=0
     i3s=m2
  end if
END SUBROUTINE xc_dimensions


!> Calculate four sets of dimension needed for the calculation of the
!! convolution for the periodic system
!!
!! @warning
!!    This four sets of dimensions are actually redundant (mi=n0i), 
!!    due to the backward-compatibility
!!    with the other geometries of the Poisson Solver.
!!    The dimensions 2 and 3 are exchanged.
!! @author Luigi Genovese
!! @date October 2006
subroutine P_FFT_dimensions(n01,n02,n03,m1,m2,m3,n1,n2,n3,md1,md2,md3,nd1,nd2,nd3,nproc,enlarge_md2)
 implicit none
 !Arguments
 logical, intent(in) :: enlarge_md2
 integer, intent(in) :: n01,n02,n03   !< original real dimensions (input)
 integer, intent(in) :: nproc
 integer, intent(out) :: m1,m2,m3     !< original real dimension, with m2 and m3 exchanged
 integer, intent(out) :: n1,n2,n3     !< the first FFT dimensions (even for the moment - the medium point being n/2+1)
 integer, intent(out) :: md1,md2,md3  !< the n1,n2,n3 dimensions. They contain the real unpadded space.
                                      !! md2 is further enlarged to be a multiple of nproc
 integer, intent(out) :: nd1,nd2,nd3  !< Fourier dimensions for which the kernel is injective,
                                      !! formally 1/8 of the fourier grid. Here the dimension nd3 is
                                      !! enlarged to be a multiple of nproc
 !Local variables
 integer :: l1,l2,l3
 
 !dimensions of the density in the real space
 m1=n01
 m2=n03
 m3=n02

 ! real space grid dimension (suitable for number of processors)
 l1=m1
 l2=m2
 l3=m3 

 !initialize the n dimension to solve Cray compiler bug
 n1=l1
 n2=l2
 n3=l3

 call fourier_dim(l1,n1)
 if (n1 /= m1) then
    !print *,'the FFT in the x direction is not allowed'
    !print *,'n01 dimension',n01
    !stop
    call f_err_throw('The FFT in the x direction is not allowed, n01 dimension '//n01)
 end if

 call fourier_dim(l2,n2)
 if (n2 /= m2) then
    !print *,'the FFT in the z direction is not allowed'
    !print *,'n03 dimension',n03
    !stop
    call f_err_throw('The FFT in the z direction is not allowed, n03 dimension '//n03)
 end if
 
 call fourier_dim(l3,n3)
 if (n3 /= m3) then
    !print *,'the FFT in the y direction is not allowed'
    !print *,'n02 dimension',n02
    !stop
    call f_err_throw('The FFT in the y direction is not allowed, n02 dimension '//n02)
 end if

 !dimensions that contain the unpadded real space,
 ! compatible with the number of processes
 md1=n1
 md2=n2
 md3=n3

 !enlarge the md2 dimension to be compatible with MPI_ALLTOALL communication
 do while(nproc*(md2/nproc) < n2)
    !151 if (nproc*(md2/nproc) < n2) then
    md2=md2+1
 end do
!    goto 151
 !endif
 
 if (enlarge_md2) md2=(md2/nproc+1)*nproc

 !dimensions of the kernel, 1/8 of the total volume,
 !compatible with nproc
 nd1=n1/2+1 
 nd2=n2/2+1
 nd3=n3/2+1
 
 !enlarge the md2 dimension to be compatible with MPI_ALLTOALL communication
 do while(modulo(nd3,nproc) /= 0)
!250 if (modulo(nd3,nproc) /= 0) then
    nd3=nd3+1
!    goto 250
! endif
 end do

END SUBROUTINE P_FFT_dimensions


!> Calculate four sets of dimension needed for the calculation of the
!! convolution for the surface system
!!
!! SYNOPSIS
!!    n01,n02,n03 original real dimensions (input)
!!
!!    m1,m2,m3 original real dimension, with 2 and 3 exchanged
!!
!!    n1,n2 the first FFT dimensions, for the moment supposed to be even
!!    n3    the double of the first FFT even dimension greater than m3
!!          (improved for the HalFFT procedure)
!!
!!    md1,md2     the n1,n2 dimensions. 
!!    md3         half of n3 dimension. They contain the real unpadded space,
!!                which has been properly enlarged to be compatible with the FFT dimensions n_i.
!!                md2 is further enlarged to be a multiple of nproc
!!
!!    nd1,nd2,nd3 fourier dimensions for which the kernel FFT is injective,
!!                formally 1/8 of the fourier grid. Here the dimension nd3 is
!!                enlarged to be a multiple of nproc
!!
!! @warning
!!    This four sets of dimensions are actually redundant (mi=n0i), 
!!    due to the backward-compatibility
!!    with the Poisson Solver with other geometries.
!!    Dimensions n02 and n03 were exchanged
!! @author Luigi Genovese
!! @date October 2006
subroutine S_FFT_dimensions(n01,n02,n03,m1,m2,m3,n1,n2,n3,md1,md2,md3,nd1,nd2,nd3,nproc,gpu,enlarge_md2)
 implicit none
 logical, intent(in) :: enlarge_md2
 integer, intent(in) :: n01,n02,n03,nproc,gpu
 integer, intent(out) :: m1,m2,m3,n1,n2,n3,md1,md2,md3,nd1,nd2,nd3
 integer :: l1,l2,l3

 !dimensions of the density in the real space
 m1=n01
 m2=n03
 m3=n02

 ! real space grid dimension (suitable for number of processors)
 l1=m1
 l2=m2
 if (gpu.eq.0) then
  l3=m3 !beware of the half dimension
 else
  l3=2*m3
 endif

 !initialize the n dimension to solve Cray compiler bug
 n1=l1
 n2=l2
 n3=l3

 call fourier_dim(l1,n1)
 if (n1 /= m1) then
    print *,'the FFT in the x direction is not allowed'
    print *,'n01 dimension',n01
    stop
 end if
 
 call fourier_dim(l2,n2)
 if (n2 /= m2) then
    print *,'the FFT in the z direction is not allowed'
    print *,'n03 dimension',n03
    stop
 end if

 do
    call fourier_dim(l3,n3)
    if (modulo(n3,2) == 0) then
       exit
    end if
    l3=l3+1
 end do
 if (gpu.eq.0) n3=2*n3

 !dimensions that contain the unpadded real space,
 ! compatible with the number of processes
 md1=n1
 md2=n2
 md3=n3/2
 do while(nproc*(md2/nproc) < n2)
    !151 if (nproc*(md2/nproc).lt.n2) then
    md2=md2+1
    !goto 151
    !endif
 end do

 if (enlarge_md2) md2=(md2/nproc+1)*nproc

 !dimensions of the kernel, 1/8 of the total volume,
 !compatible with nproc

 !these two dimensions are like that since they are even
 nd1=n1/2+1
 nd2=n2/2+1
 nd3=n3/2+1
 do while(modulo(nd3,nproc) /= 0)
    !250 if (modulo(nd3,nproc) .ne. 0) then
    nd3=nd3+1
    !goto 250
    !endif
 end do

END SUBROUTINE S_FFT_dimensions


!> Calculate four sets of dimension needed for the calculation of the
!! convolution for the Wires BC system
!!
!! SYNOPSIS
!!    n01,n02,n03 original real dimensions (input)
!!
!!    m1,m2,m3 original real dimension, with 2 and 3 exchanged
!!
!!    n1,n2 the first FFT dimensions, for the moment supposed to be even
!!    n3    the double of the first FFT even dimension greater than m3
!!          (improved for the HalFFT procedure)
!!
!!    md1,md2     the n1,n2 dimensions. 
!!    md3         half of n3 dimension. They contain the real unpadded space,
!!                which has been properly enlarged to be compatible with the FFT dimensions n_i.
!!                md2 is further enlarged to be a multiple of nproc
!!
!!    nd1,nd2,nd3 fourier dimensions for which the kernel FFT is injective,
!!                formally 1/8 of the fourier grid. Here the dimension nd3 is
!!                enlarged to be a multiple of nproc
!!
!! @warning
!!    This four sets of dimensions are actually redundant (mi=n0i), 
!!    due to the backward-compatibility
!!    with the Poisson Solver with other geometries.
!!    Dimensions n02 and n03 were exchanged
!! @author Luigi Genovese
!! @date October 2006
subroutine W_FFT_dimensions(n01,n02,n03,m1,m2,m3,n1,n2,n3,md1,md2,md3,nd1,nd2,nd3,nproc,gpu,enlarge_md2)
 implicit none
 logical, intent(in) :: enlarge_md2
 integer, intent(in) :: n01,n02,n03,nproc,gpu
 integer, intent(out) :: m1,m2,m3,n1,n2,n3,md1,md2,md3,nd1,nd2,nd3
 integer :: l1,l2,l3

 !dimensions of the density in the real space
 m1=n01
 m2=n03
 m3=n02

 ! real space grid dimension (suitable for number of processors)
 l1=2*m1
 l2=m2
 if (gpu.eq.0) then
  l3=m3 !beware of the half dimension
 else
  l3=2*m3
 endif

 do
    call fourier_dim(l1,n1)
    if (modulo(n1,2) == 0) then
       exit
    end if
    l1=l1+1
 end do

 
 call fourier_dim(l2,n2)
 if (n2 /= m2) then
    print *,'the FFT in the z direction is not allowed'
    print *,'n03 dimension',n03
    stop
 end if

 do
    call fourier_dim(l3,n3)
    if (modulo(n3,2) == 0) then
       exit
    end if
    l3=l3+1
 end do

 if (gpu.eq.0) n3=2*n3

 !dimensions that contain the unpadded real space,
 ! compatible with the number of processes
 md1=n1/2
 md2=n2
 md3=n3/2
 do while(nproc*(md2/nproc) < n2)
    !151 if (nproc*(md2/nproc).lt.n2) then
    md2=md2+1
    !goto 151
    !endif
 end do

 if (enlarge_md2) md2=(md2/nproc+1)*nproc

 !dimensions of the kernel, 1/8 of the total volume,
 !compatible with nproc

 !these two dimensions are like that since they are even
 nd1=n1/2+1
 nd2=n2/2+1
 nd3=n3/2+1
 do while(modulo(nd3,nproc) /= 0)
    !250 if (modulo(nd3,nproc) .ne. 0) then
    nd3=nd3+1
    !goto 250
    !endif
 end do

END SUBROUTINE W_FFT_dimensions


!> Calculate four sets of dimension needed for the calculation of the
!! zero-padded convolution
!!
!!    @param n01,n02,n03 original real dimensions (input)
!!
!!    @param m1,m2,m3 original real dimension with the dimension 2 and 3 exchanged
!!
!!    @param n1,n2 the first FFT even dimensions greater that 2*m1, 2*m2
!!    @param n3    the double of the first FFT even dimension greater than m3
!!          (improved for the HalFFT procedure)
!!
!!    @param md1,md2,md3 half of n1,n2,n3 dimension. They contain the real unpadded space,
!!                which has been properly enlarged to be compatible with the FFT dimensions n_i.
!!                md2 is further enlarged to be a multiple of nproc
!!
!!    @param nd1,nd2,nd3 fourier dimensions for which the kernel FFT is injective,
!!                formally 1/8 of the fourier grid. Here the dimension nd3 is
!!                enlarged to be a multiple of nproc
!!
!! @warning
!!    The dimension m2 and m3 correspond to n03 and n02 respectively
!!    this is needed since the convolution routine manage arrays of dimension
!!    (md1,md3,md2/nproc)
!! @author Luigi Genovese
!! @date February 2006
subroutine F_FFT_dimensions(n01,n02,n03,m1,m2,m3,n1,n2,n3,md1,md2,md3,nd1,nd2,nd3,nproc,gpu,enlarge_md2)
 implicit none
 logical, intent(in) :: enlarge_md2
 integer, intent(in) :: n01,n02,n03,nproc,gpu
 integer, intent(out) :: m1,m2,m3,n1,n2,n3,md1,md2,md3,nd1,nd2,nd3
 integer :: l1,l2,l3

 !dimensions of the density in the real space, inverted for convenience
 m1=n01
 m2=n03
 m3=n02
 ! real space grid dimension (suitable for number of processors)
 l1=2*m1
 l2=2*m2
 if (gpu.eq.0) then
  l3=m3 !beware of the half dimension
 else
  l3=2*m3
 endif
 !initialize the n dimension to solve Cray compiler bug
 n1=l1
 n2=l2
 n3=l3
 do
    call fourier_dim(l1,n1)
    if (modulo(n1,2) == 0) then
       exit
    end if
    l1=l1+1
 end do
 do
    call fourier_dim(l2,n2)
    if (modulo(n2,2) == 0) then
       exit
    end if
    l2=l2+1
 end do
 do
    call fourier_dim(l3,n3)
    if (modulo(n3,2) == 0) then
       exit
    end if
    l3=l3+1
 end do
 if (gpu.eq.0) n3=2*n3

 !dimensions that contain the unpadded real space,
 ! compatible with the number of processes
 md1=n1/2
 md2=n2/2
 md3=n3/2
 do while(nproc*(md2/nproc) < n2/2)
   !151 if (nproc*(md2/nproc).lt.n2/2) then
    md2=md2+1
   !goto 151
   !endif
 end do

 if (enlarge_md2) md2=(md2/nproc+1)*nproc

 !dimensions of the kernel, 1/8 of the total volume,
 !compatible with nproc
 nd1=n1/2+1
 nd2=n2/2+1
 nd3=n3/2+1

 do while(modulo(nd3,nproc) /= 0)
    !250 if (modulo(nd3,nproc) .ne. 0) then
    nd3=nd3+1
    !    goto 250
    ! endif
 end do

END SUBROUTINE F_FFT_dimensions
