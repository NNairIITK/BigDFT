!> @file
!! Copy the different type used by linear version
!! @author
!!    Copyright (C) 2011-2013 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS


!> Currently incomplete - need to add comms arrays etc
subroutine copy_tmbs(iproc, tmbin, tmbout, subname)
  use module_base
  use module_types
  use module_interfaces
  implicit none

  integer,intent(in) :: iproc
  type(DFT_wavefunction), intent(in) :: tmbin
  type(DFT_wavefunction), intent(out) :: tmbout
  character(len=*),intent(in):: subname

  call f_routine(id='copy_tmbs')

  call nullify_orbitals_data(tmbout%orbs)
  call copy_orbitals_data(tmbin%orbs, tmbout%orbs, subname)
  call nullify_local_zone_descriptors(tmbout%lzd)
  call copy_old_supportfunctions(iproc,tmbin%orbs,tmbin%lzd,tmbin%psi,tmbout%lzd,tmbout%psi)

  tmbout%npsidim_orbs = tmbin%npsidim_orbs

  if (associated(tmbin%coeff)) then !(in%lin%scf_mode/=LINEAR_FOE) then ! should move this check to copy_old_coeffs
      call copy_old_coefficients(tmbin%orbs%norb, tmbin%coeff, tmbout%coeff)
  else
      nullify(tmbout%coeff)
  end if

  ! Parts of tmbout%lzd have been allocated in copy_old_supportfunctions, so deallocate everything and reallocate everything
  ! properly. Of course this is a very bad solution.
  call deallocate_local_zone_descriptors(tmbout%lzd, subname)
  call copy_local_zone_descriptors(tmbin%lzd, tmbout%lzd, subname)

  call copy_linear_matrices(tmbin%linmat, tmbout%linmat)

  call copy_comms_linear(tmbin%collcom, tmbout%collcom)

  ! should technically copy these across as well but not needed for restart and will eventually be removing wfnmd as a type
  !nullify(tmbout%linmat%denskern%matrix_compr)
  !nullify(tmbout%linmat%denskern_large%matrix_compr)

  ! should also copy/nullify p2pcomms etc

  !call copy_old_inwhichlocreg(tmbin%orbs%norb, tmbin%orbs%inwhichlocreg, tmbout%orbs%inwhichlocreg, &
  !     tmbin%orbs%onwhichatom, tmbout%orbs%onwhichatom)

  call f_release_routine()

end subroutine copy_tmbs

subroutine copy_convolutions_bounds(geocode,boundsin, boundsout, subname)
  use module_base
  use module_types
  implicit none
  
  ! Calling arguments
  character(len=1),intent(in) :: geocode !< @copydoc poisson_solver::doc::geocode
  type(convolutions_bounds),intent(in):: boundsin
  type(convolutions_bounds),intent(inout):: boundsout
  character(len=*),intent(in):: subname
  
  ! Local variables
  integer:: iis1, iie1, iis2, iie2, iis3, iie3, i1, i2, i3, istat, iall
  
  call copy_kinetic_bounds(geocode, boundsin%kb, boundsout%kb, subname)
  call copy_shrink_bounds(geocode, boundsin%sb, boundsout%sb, subname)
  call copy_grow_bounds(geocode, boundsin%gb, boundsout%gb, subname)
  
  if(geocode == 'F') then
     if(associated(boundsout%ibyyzz_r)) then
         call f_free_ptr(boundsout%ibyyzz_r)
     end if
  
     if(associated(boundsin%ibyyzz_r)) then
         iis1=lbound(boundsin%ibyyzz_r,1)
         iie1=ubound(boundsin%ibyyzz_r,1)
         iis2=lbound(boundsin%ibyyzz_r,2)
         iie2=ubound(boundsin%ibyyzz_r,2)
         iis3=lbound(boundsin%ibyyzz_r,3)
         iie3=ubound(boundsin%ibyyzz_r,3)
         boundsout%ibyyzz_r = f_malloc_ptr((/ iis1.to.iie1,iis2.to.iie2,iis3.to.iie3 /),id='boundsout%ibyyzz_r')
         do i3=iis3,iie3
             do i2=iis2,iie2
                 do i1=iis1,iie1
                     boundsout%ibyyzz_r(i1,i2,i3) = boundsin%ibyyzz_r(i1,i2,i3)
                 end do
             end do
         end do
     end if
  end if
end subroutine copy_convolutions_bounds

subroutine copy_kinetic_bounds(geocode,kbin, kbout, subname)
use module_base
use module_types
implicit none

! Calling arguments
character(len=1),intent(in) :: geocode !< @copydoc poisson_solver::doc::geocode 
type(kinetic_bounds),intent(in):: kbin
type(kinetic_bounds),intent(inout):: kbout
character(len=*),intent(in):: subname

! Local variables
integer:: iis1, iie1, iis2, iie2, iis3, iie3, i1, i2, i3, istat, iall

if(geocode == 'F') then
   if(associated(kbout%ibyz_c)) then
       call f_free_ptr(kbout%ibyz_c)
   end if
   if(associated(kbin%ibyz_c)) then
       iis1=lbound(kbin%ibyz_c,1)
       iie1=ubound(kbin%ibyz_c,1)
       iis2=lbound(kbin%ibyz_c,2)
       iie2=ubound(kbin%ibyz_c,2)
       iis3=lbound(kbin%ibyz_c,3)
       iie3=ubound(kbin%ibyz_c,3)
       kbout%ibyz_c = f_malloc_ptr((/ iis1.to.iie1 , iis2.to.iie2 , iis3.to.iie3 /),id='kbout%ibyz_c')
       do i3=iis3,iie3
           do i2=iis2,iie2
               do i1=iis1,iie1
                   kbout%ibyz_c(i1,i2,i3) = kbin%ibyz_c(i1,i2,i3)
               end do
           end do
       end do
   end if
   
   
   if(associated(kbout%ibxz_c)) then
       call f_free_ptr(kbout%ibxz_c)
   end if
   if(associated(kbin%ibxz_c)) then
       iis1=lbound(kbin%ibxz_c,1)
       iie1=ubound(kbin%ibxz_c,1)
       iis2=lbound(kbin%ibxz_c,2)
       iie2=ubound(kbin%ibxz_c,2)
       iis3=lbound(kbin%ibxz_c,3)
       iie3=ubound(kbin%ibxz_c,3)
       kbout%ibxz_c = f_malloc_ptr((/ iis1.to.iie1 , iis2.to.iie2 , iis3.to.iie3 /),id='kbout%ibxz_c')
       do i3=iis3,iie3
           do i2=iis2,iie2
               do i1=iis1,iie1
                   kbout%ibxz_c(i1,i2,i3) = kbin%ibxz_c(i1,i2,i3)
               end do
           end do
       end do
   end if
   
   
   if(associated(kbout%ibxy_c)) then
       call f_free_ptr(kbout%ibxy_c)
   end if
   if(associated(kbin%ibxy_c)) then
       iis1=lbound(kbin%ibxy_c,1)
       iie1=ubound(kbin%ibxy_c,1)
       iis2=lbound(kbin%ibxy_c,2)
       iie2=ubound(kbin%ibxy_c,2)
       iis3=lbound(kbin%ibxy_c,3)
       iie3=ubound(kbin%ibxy_c,3)
       kbout%ibxy_c = f_malloc_ptr((/ iis1.to.iie1 , iis2.to.iie2 , iis3.to.iie3 /),id='kbout%ibxy_c')
       do i3=iis3,iie3
           do i2=iis2,iie2
               do i1=iis1,iie1
                   kbout%ibxy_c(i1,i2,i3) = kbin%ibxy_c(i1,i2,i3)
               end do
           end do
       end do
   end if
end if

if(associated(kbout%ibyz_f)) then
    call f_free_ptr(kbout%ibyz_f)
end if
if(associated(kbin%ibyz_f)) then
    iis1=lbound(kbin%ibyz_f,1)
    iie1=ubound(kbin%ibyz_f,1)
    iis2=lbound(kbin%ibyz_f,2)
    iie2=ubound(kbin%ibyz_f,2)
    iis3=lbound(kbin%ibyz_f,3)
    iie3=ubound(kbin%ibyz_f,3)
    kbout%ibyz_f = f_malloc_ptr((/ iis1.to.iie1, iis2.to.iie2, iis3.to.iie3 /),id='kbout%ibyz_f')
    do i3=iis3,iie3
        do i2=iis2,iie2
            do i1=iis1,iie1
                kbout%ibyz_f(i1,i2,i3) = kbin%ibyz_f(i1,i2,i3)
            end do
        end do
    end do
end if


if(associated(kbout%ibxz_f)) then
    call f_free_ptr(kbout%ibxz_f)
end if
if(associated(kbin%ibxz_f)) then
    iis1=lbound(kbin%ibxz_f,1)
    iie1=ubound(kbin%ibxz_f,1)
    iis2=lbound(kbin%ibxz_f,2)
    iie2=ubound(kbin%ibxz_f,2)
    iis3=lbound(kbin%ibxz_f,3)
    iie3=ubound(kbin%ibxz_f,3)
    kbout%ibxz_f = f_malloc_ptr((/ iis1.to.iie1, iis2.to.iie2, iis3.to.iie3 /),id='kbout%ibxz_f')
    do i3=iis3,iie3
        do i2=iis2,iie2
            do i1=iis1,iie1
                kbout%ibxz_f(i1,i2,i3) = kbin%ibxz_f(i1,i2,i3)
            end do
        end do
    end do
end if


if(associated(kbout%ibxy_f)) then
    call f_free_ptr(kbout%ibxy_f)
end if
if(associated(kbin%ibxy_f)) then
    iis1=lbound(kbin%ibxy_f,1)
    iie1=ubound(kbin%ibxy_f,1)
    iis2=lbound(kbin%ibxy_f,2)
    iie2=ubound(kbin%ibxy_f,2)
    iis3=lbound(kbin%ibxy_f,3)
    iie3=ubound(kbin%ibxy_f,3)
    kbout%ibxy_f = f_malloc_ptr((/ iis1.to.iie1, iis2.to.iie2, iis3.to.iie3 /),id='kbout%ibxy_f')
    do i3=iis3,iie3
        do i2=iis2,iie2
            do i1=iis1,iie1
                kbout%ibxy_f(i1,i2,i3) = kbin%ibxy_f(i1,i2,i3)
            end do
        end do
    end do
end if


end subroutine copy_kinetic_bounds




subroutine copy_shrink_bounds(geocode, sbin, sbout, subname)
use module_base
use module_types
implicit none

! Calling arguments
character(len=1), intent(in) :: geocode !< @copydoc poisson_solver::doc::geocode
type(shrink_bounds),intent(in):: sbin
type(shrink_bounds),intent(inout):: sbout
character(len=*),intent(in):: subname

! Local variables
integer:: iis1, iie1, iis2, iie2, iis3, iie3, i1, i2, i3, istat, iall

if(geocode == 'F') then
   if(associated(sbout%ibzzx_c)) then
       call f_free_ptr(sbout%ibzzx_c)
   end if
   if(associated(sbin%ibzzx_c)) then
       iis1=lbound(sbin%ibzzx_c,1)
       iie1=ubound(sbin%ibzzx_c,1)
       iis2=lbound(sbin%ibzzx_c,2)
       iie2=ubound(sbin%ibzzx_c,2)
       iis3=lbound(sbin%ibzzx_c,3)
       iie3=ubound(sbin%ibzzx_c,3)
       sbout%ibzzx_c = f_malloc_ptr((/ iis1.to.iie1, iis2.to.iie2, iis3.to.iie3 /),id='sbout%ibzzx_c')
       do i3=iis3,iie3
           do i2=iis2,iie2
               do i1=iis1,iie1
                   sbout%ibzzx_c(i1,i2,i3) = sbin%ibzzx_c(i1,i2,i3)
               end do
           end do
       end do
   end if
   
   
   if(associated(sbout%ibyyzz_c)) then
       call f_free_ptr(sbout%ibyyzz_c)
   end if
   if(associated(sbin%ibyyzz_c)) then
       iis1=lbound(sbin%ibyyzz_c,1)
       iie1=ubound(sbin%ibyyzz_c,1)
       iis2=lbound(sbin%ibyyzz_c,2)
       iie2=ubound(sbin%ibyyzz_c,2)
       iis3=lbound(sbin%ibyyzz_c,3)
       iie3=ubound(sbin%ibyyzz_c,3)
       sbout%ibyyzz_c = f_malloc_ptr((/ iis1.to.iie1, iis2.to.iie2, iis3.to.iie3 /),id='sbout%ibyyzz_c')
       do i3=iis3,iie3
           do i2=iis2,iie2
               do i1=iis1,iie1
                   sbout%ibyyzz_c(i1,i2,i3) = sbin%ibyyzz_c(i1,i2,i3)
               end do
           end do
       end do
   end if
end if

if(associated(sbout%ibxy_ff)) then
    call f_free_ptr(sbout%ibxy_ff)
end if
if(associated(sbin%ibxy_ff)) then
    iis1=lbound(sbin%ibxy_ff,1)
    iie1=ubound(sbin%ibxy_ff,1)
    iis2=lbound(sbin%ibxy_ff,2)
    iie2=ubound(sbin%ibxy_ff,2)
    iis3=lbound(sbin%ibxy_ff,3)
    iie3=ubound(sbin%ibxy_ff,3)
    sbout%ibxy_ff = f_malloc_ptr((/ iis1.to.iie1, iis2.to.iie2, iis3.to.iie3 /),id='sbout%ibxy_ff')
    do i3=iis3,iie3
        do i2=iis2,iie2
            do i1=iis1,iie1
                sbout%ibxy_ff(i1,i2,i3) = sbin%ibxy_ff(i1,i2,i3)
            end do
        end do
    end do
end if


if(associated(sbout%ibzzx_f)) then
    call f_free_ptr(sbout%ibzzx_f)
end if
if(associated(sbin%ibzzx_f)) then
    iis1=lbound(sbin%ibzzx_f,1)
    iie1=ubound(sbin%ibzzx_f,1)
    iis2=lbound(sbin%ibzzx_f,2)
    iie2=ubound(sbin%ibzzx_f,2)
    iis3=lbound(sbin%ibzzx_f,3)
    iie3=ubound(sbin%ibzzx_f,3)
    sbout%ibzzx_f = f_malloc_ptr((/ iis1.to.iie1, iis2.to.iie2, iis3.to.iie3 /),id='sbout%ibzzx_f')
    do i3=iis3,iie3
        do i2=iis2,iie2
            do i1=iis1,iie1
                sbout%ibzzx_f(i1,i2,i3) = sbin%ibzzx_f(i1,i2,i3)
            end do
        end do
    end do
end if


if(associated(sbout%ibyyzz_f)) then
    call f_free_ptr(sbout%ibyyzz_f)
end if
if(associated(sbin%ibyyzz_f)) then
    iis1=lbound(sbin%ibyyzz_f,1)
    iie1=ubound(sbin%ibyyzz_f,1)
    iis2=lbound(sbin%ibyyzz_f,2)
    iie2=ubound(sbin%ibyyzz_f,2)
    iis3=lbound(sbin%ibyyzz_f,3)
    iie3=ubound(sbin%ibyyzz_f,3)
    sbout%ibyyzz_f = f_malloc_ptr((/ iis1.to.iie1, iis2.to.iie2, iis3.to.iie3 /),id='sbout%ibyyzz_f')
    do i3=iis3,iie3
        do i2=iis2,iie2
            do i1=iis1,iie1
                sbout%ibyyzz_f(i1,i2,i3) = sbin%ibyyzz_f(i1,i2,i3)
            end do
        end do
    end do
end if



end subroutine copy_shrink_bounds




subroutine copy_grow_bounds(geocode, gbin, gbout, subname)
use module_base
use module_types
implicit none

! Calling arguments
character(len=1),intent(in) :: geocode !< @copydoc poisson_solver::doc::geocode
type(grow_bounds),intent(in):: gbin
type(grow_bounds),intent(inout):: gbout
character(len=*),intent(in):: subname

! Local variables
integer:: iis1, iie1, iis2, iie2, iis3, iie3, i1, i2, i3, istat, iall

if(geocode == 'F')then
   if(associated(gbout%ibzxx_c)) then
       call f_free_ptr(gbout%ibzxx_c)
   end if
   if(associated(gbin%ibzxx_c)) then
       iis1=lbound(gbin%ibzxx_c,1)
       iie1=ubound(gbin%ibzxx_c,1)
       iis2=lbound(gbin%ibzxx_c,2)
       iie2=ubound(gbin%ibzxx_c,2)
       iis3=lbound(gbin%ibzxx_c,3)
       iie3=ubound(gbin%ibzxx_c,3)
       gbout%ibzxx_c = f_malloc_ptr((/ iis1.to.iie1, iis2.to.iie2, iis3.to.iie3 /),id='gbout%ibzxx_c')
       do i3=iis3,iie3
           do i2=iis2,iie2
               do i1=iis1,iie1
                   gbout%ibzxx_c(i1,i2,i3) = gbin%ibzxx_c(i1,i2,i3)
               end do
           end do
       end do
   end if
       
   
   if(associated(gbout%ibxxyy_c)) then
       call f_free_ptr(gbout%ibxxyy_c)
   end if
   if(associated(gbin%ibxxyy_c)) then
       iis1=lbound(gbin%ibxxyy_c,1)
       iie1=ubound(gbin%ibxxyy_c,1)
       iis2=lbound(gbin%ibxxyy_c,2)
       iie2=ubound(gbin%ibxxyy_c,2)
       iis3=lbound(gbin%ibxxyy_c,3)
       iie3=ubound(gbin%ibxxyy_c,3)
       gbout%ibxxyy_c = f_malloc_ptr((/ iis1.to.iie1, iis2.to.iie2, iis3.to.iie3 /),id='gbout%ibxxyy_c')
       do i3=iis3,iie3
           do i2=iis2,iie2
               do i1=iis1,iie1
                   gbout%ibxxyy_c(i1,i2,i3) = gbin%ibxxyy_c(i1,i2,i3)
               end do
           end do
       end do
   end if
end if

if(associated(gbout%ibyz_ff)) then
    call f_free_ptr(gbout%ibyz_ff)
end if
if(associated(gbin%ibyz_ff)) then
    iis1=lbound(gbin%ibyz_ff,1)
    iie1=ubound(gbin%ibyz_ff,1)
    iis2=lbound(gbin%ibyz_ff,2)
    iie2=ubound(gbin%ibyz_ff,2)
    iis3=lbound(gbin%ibyz_ff,3)
    iie3=ubound(gbin%ibyz_ff,3)
    gbout%ibyz_ff = f_malloc_ptr((/ iis1.to.iie1, iis2.to.iie2, iis3.to.iie3 /),id='gbout%ibyz_ff')
    do i3=iis3,iie3
        do i2=iis2,iie2
            do i1=iis1,iie1
                gbout%ibyz_ff(i1,i2,i3) = gbin%ibyz_ff(i1,i2,i3)
            end do
        end do
    end do
end if

if(associated(gbout%ibzxx_f)) then
    call f_free_ptr(gbout%ibzxx_f)
end if
if(associated(gbin%ibzxx_f)) then
    iis1=lbound(gbin%ibzxx_f,1)
    iie1=ubound(gbin%ibzxx_f,1)
    iis2=lbound(gbin%ibzxx_f,2)
    iie2=ubound(gbin%ibzxx_f,2)
    iis3=lbound(gbin%ibzxx_f,3)
    iie3=ubound(gbin%ibzxx_f,3)
    gbout%ibzxx_f = f_malloc_ptr((/ iis1.to.iie1, iis2.to.iie2, iis3.to.iie3 /),id='gbout%ibzxx_f')
    do i3=iis3,iie3
        do i2=iis2,iie2
            do i1=iis1,iie1
                gbout%ibzxx_f(i1,i2,i3) = gbin%ibzxx_f(i1,i2,i3)
            end do
        end do
    end do
end if


if(associated(gbout%ibxxyy_f)) then
    call f_free_ptr(gbout%ibxxyy_f)
end if
if(associated(gbin%ibxxyy_f)) then
    iis1=lbound(gbin%ibxxyy_f,1)
    iie1=ubound(gbin%ibxxyy_f,1)
    iis2=lbound(gbin%ibxxyy_f,2)
    iie2=ubound(gbin%ibxxyy_f,2)
    iis3=lbound(gbin%ibxxyy_f,3)
    iie3=ubound(gbin%ibxxyy_f,3)
    gbout%ibxxyy_f = f_malloc_ptr((/ iis1.to.iie1, iis2.to.iie2, iis3.to.iie3 /),id='gbout%ibxxyy_f')
    do i3=iis3,iie3
        do i2=iis2,iie2
            do i1=iis1,iie1
                gbout%ibxxyy_f(i1,i2,i3) = gbin%ibxxyy_f(i1,i2,i3)
            end do
        end do
    end do
end if


end subroutine copy_grow_bounds

subroutine copy_orbitals_data(orbsin, orbsout, subname)
use module_base
use module_types
implicit none

! Calling arguments
type(orbitals_data),intent(in):: orbsin
type(orbitals_data),intent(inout):: orbsout
character(len=*),intent(in):: subname

! Local variables
integer:: iis1, iie1, iis2, iie2, i1, i2, istat, iall

orbsout%norb = orbsin%norb
orbsout%norbp = orbsin%norbp
orbsout%norbu = orbsin%norbu
orbsout%norbd = orbsin%norbd
orbsout%nspin = orbsin%nspin
orbsout%nspinor = orbsin%nspinor
orbsout%isorb = orbsin%isorb
orbsout%npsidim_orbs = orbsin%npsidim_orbs
orbsout%npsidim_comp = orbsin%npsidim_comp
orbsout%nkpts = orbsin%nkpts
orbsout%nkptsp = orbsin%nkptsp
orbsout%iskpts = orbsin%iskpts
orbsout%efermi = orbsin%efermi

call f_free_ptr(orbsout%norb_par)
if(associated(orbsin%norb_par)) then
   orbsout%norb_par = &
        f_malloc_ptr(src=orbsin%norb_par,lbounds=lbound(orbsin%norb_par),id='orbsout%norb_par')
!!$
!!$    iis1=lbound(orbsin%norb_par,1)
!!$    iie1=ubound(orbsin%norb_par,1)
!!$    iis2=lbound(orbsin%norb_par,2)
!!$    iie2=ubound(orbsin%norb_par,2)
!!$    orbsout%norb_par = f_malloc_ptr((/ iis1.to.iie1 , iis2.to.iie2 /),id='orbsout%norb_par')
!!$    do i1=iis1,iie1
!!$       do i2 = iis2,iie2
!!$        orbsout%norb_par(i1,i2) = orbsin%norb_par(i1,i2)
!!$       end do
!!$    end do
end if

call f_free_ptr(orbsout%iokpt)
if(associated(orbsin%iokpt)) then
   orbsout%iokpt = f_malloc_ptr(src=orbsin%iokpt,lbounds=lbound(orbsin%iokpt),id='orbsout%iokpt')
!!$    iis1=lbound(orbsin%iokpt,1)
!!$    iie1=ubound(orbsin%iokpt,1)
!!$    orbsout%iokpt = f_malloc_ptr(iis1.to.iie1,id='orbsout%iokpt')
!!$    do i1=iis1,iie1
!!$        orbsout%iokpt(i1) = orbsin%iokpt(i1)
!!$    end do
end if

call f_free_ptr(orbsout%ikptproc)
if(associated(orbsin%ikptproc)) then
   orbsout%ikptproc = f_malloc_ptr(src=orbsin%ikptproc,&
        lbounds=lbound(orbsin%ikptproc),id='orbsout%ikptproc')
!!$    iis1=lbound(orbsin%ikptproc,1)
!!$    iie1=ubound(orbsin%ikptproc,1)
!!$    orbsout%ikptproc = f_malloc_ptr(iis1.to.iie1,id='orbsout%ikptproc')
!!$    do i1=iis1,iie1
!!$        orbsout%ikptproc(i1) = orbsin%ikptproc(i1)
!!$    end do
end if

call f_free_ptr(orbsout%inwhichlocreg)
    
if(associated(orbsin%inwhichlocreg)) then
   orbsout%inwhichlocreg = &
        f_malloc_ptr(src=orbsin%inwhichlocreg,&
        lbounds=lbound(orbsin%inwhichlocreg),id='orbsout%inwhichlocreg')
       
!!$    iis1=lbound(orbsin%inwhichlocreg,1)
!!$    iie1=ubound(orbsin%inwhichlocreg,1)
!!$    orbsout%inwhichlocreg = f_malloc_ptr(iis1.to.iie1,id='orbsout%inwhichlocreg')
!!$    do i1=iis1,iie1
!!$        orbsout%inwhichlocreg(i1) = orbsin%inwhichlocreg(i1)
!!$    end do
end if

call f_free_ptr(orbsout%onwhichatom)
if(associated(orbsin%onwhichatom)) then
    orbsout%onwhichatom = &
         f_malloc_ptr(src=orbsin%onwhichatom,lbounds=lbound(orbsin%onwhichatom),&
         id='orbsout%onwhichatom')
!!$    iis1=lbound(orbsin%onwhichatom,1)
!!$    iie1=ubound(orbsin%onwhichatom,1)
!!$    orbsout%onwhichatom = f_malloc_ptr(iis1.to.iie1,id='orbsout%onwhichatom')
!!$    do i1=iis1,iie1
!!$        orbsout%onwhichatom(i1) = orbsin%onwhichatom(i1)
!!$    end do
end if

call f_free_ptr(orbsout%isorb_par)
if(associated(orbsin%isorb_par)) then
   orbsout%isorb_par = f_malloc_ptr(src=orbsin%isorb_par,id='orbsout%isorb_par')
!!$    iis1=lbound(orbsin%isorb_par,1)
!!$    iie1=ubound(orbsin%isorb_par,1)
!!$    orbsout%isorb_par = f_malloc_ptr(iis1.to.iie1,id='orbsout%isorb_par')
!!$    do i1=iis1,iie1
!!$        orbsout%isorb_par(i1) = orbsin%isorb_par(i1)
!!$    end do
end if

    call f_free_ptr(orbsout%eval)
    if(associated(orbsin%eval)) then
       orbsout%eval = f_malloc_ptr(src=orbsin%eval,id='orbsout%eval')
!!$    iis1=lbound(orbsin%eval,1)
!!$    iie1=ubound(orbsin%eval,1)
!!$    if(iie1 /= iis1 ) then
!!$       orbsout%eval = f_malloc_ptr(iis1.to.iie1,id='orbsout%eval')
!!$       do i1=iis1,iie1
!!$           orbsout%eval(i1) = orbsin%eval(i1)
!!$       end do
!!$    end if
end if

    call f_free_ptr(orbsout%occup)
if(associated(orbsin%occup)) then
   orbsout%occup = f_malloc_ptr(src=orbsin%occup,id='orbsout%occup')
!!$    iis1=lbound(orbsin%occup,1)
!!$    iie1=ubound(orbsin%occup,1)
!!$    orbsout%occup = f_malloc_ptr(iis1.to.iie1,id='orbsout%occup')
!!$    do i1=iis1,iie1
!!$        orbsout%occup(i1) = orbsin%occup(i1)
!!$    end do
end if


    call f_free_ptr(orbsout%spinsgn)

if(associated(orbsin%spinsgn)) then
   orbsout%spinsgn = f_malloc_ptr(src=orbsin%spinsgn,id='orbsout%spinsgn')
!!$    iis1=lbound(orbsin%spinsgn,1)
!!$    iie1=ubound(orbsin%spinsgn,1)
!!$    orbsout%spinsgn = f_malloc_ptr(iis1.to.iie1,id='orbsout%spinsgn')
!!$    do i1=iis1,iie1
!!$        orbsout%spinsgn(i1) = orbsin%spinsgn(i1)
!!$    end do
end if


   call f_free_ptr(orbsout%kwgts)
if(associated(orbsin%kwgts)) then
   orbsout%kwgts = f_malloc_ptr(src=orbsin%kwgts,id='orbsout%kwgts')
!!$    iis1=lbound(orbsin%kwgts,1)
!!$    iie1=ubound(orbsin%kwgts,1)
!!$    orbsout%kwgts = f_malloc_ptr(iis1.to.iie1,id='orbsout%kwgts')
!!$    do i1=iis1,iie1
!!$        orbsout%kwgts(i1) = orbsin%kwgts(i1)
!!$    end do
end if

call f_free_ptr(orbsout%kpts)
if(associated(orbsin%kpts)) then
   orbsout%kpts = f_malloc_ptr(src=orbsin%kpts,id='orbsout%kpts')
!!$    iis1=lbound(orbsin%kpts,1)
!!$    iie1=ubound(orbsin%kpts,1)
!!$    iis2=lbound(orbsin%kpts,2)
!!$    iie2=ubound(orbsin%kpts,2)
!!$    orbsout%kpts = f_malloc_ptr((/ iis1.to.iie1 , iis2.to.iie2 /),id='orbsout%kpts')
!!$    do i2=iis2,iie2
!!$        do i1=iis1,iie1
!!$            orbsout%kpts(i1,i2) = orbsin%kpts(i1,i2)
!!$        end do
!!$    end do
end if

call f_free_ptr(orbsout%ispot)
if(associated(orbsin%ispot)) then
   orbsout%ispot = f_malloc_ptr(src=orbsin%ispot,id='orbsout%ispot')
!!$    iis1=lbound(orbsin%ispot,1)
!!$    iie1=ubound(orbsin%ispot,1)
!!$    orbsout%ispot = f_malloc_ptr(iis1.to.iie1,id='orbsout%ispot')
!!$    do i1=iis1,iie1
!!$        orbsout%ispot(i1) = orbsin%ispot(i1)
!!$    end do
end if


end subroutine copy_orbitals_data


subroutine copy_local_zone_descriptors(lzd_in, lzd_out, subname)
  use locregs
  use module_types, only: local_zone_descriptors
  implicit none

  ! Calling arguments
  type(local_zone_descriptors),intent(in):: lzd_in
  type(local_zone_descriptors),intent(inout):: lzd_out
  character(len=*),intent(in):: subname

  ! Local variables
  integer:: istat, i1, iis1, iie1

  lzd_out%linear=lzd_in%linear
  lzd_out%nlr=lzd_in%nlr
  lzd_out%lintyp=lzd_in%lintyp
  lzd_out%ndimpotisf=lzd_in%ndimpotisf
  lzd_out%hgrids=lzd_in%hgrids

  call nullify_locreg_descriptors(lzd_out%glr)
  call copy_locreg_descriptors(lzd_in%glr, lzd_out%glr)

  if(associated(lzd_out%llr)) then
      deallocate(lzd_out%llr, stat=istat)
  end if
  if(associated(lzd_in%llr)) then
      iis1=lbound(lzd_in%llr,1)
      iie1=ubound(lzd_in%llr,1)
      allocate(lzd_out%llr(iis1:iie1), stat=istat)
      do i1=iis1,iie1
          call nullify_locreg_descriptors(lzd_out%llr(i1))
          call copy_locreg_descriptors(lzd_in%llr(i1), lzd_out%llr(i1))
      end do
  end if

end subroutine copy_local_zone_descriptors




!only copying sparsity pattern here, not copying whole matrix (assuming matrices not allocated)
subroutine sparse_copy_pattern(sparseMat_in, sparseMat_out, iproc, subname)
  use module_base
  use module_types
  use sparsematrix_base, only: sparse_matrix
  implicit none

  ! Calling arguments
  type(sparse_matrix),intent(in):: sparseMat_in
  type(sparse_matrix),intent(inout):: sparseMat_out
  integer, intent(in) :: iproc
  character(len=*),intent(in):: subname

  ! Local variables
  integer:: iis1, iie1, iis2, iie2, i1, i2, istat, iall

  call timing(iproc,'sparse_copy','ON')

  sparsemat_out%nseg = sparsemat_in%nseg
  sparsemat_out%store_index = sparsemat_in%store_index
  sparsemat_out%nvctr = sparsemat_in%nvctr
  sparsemat_out%nvctrp = sparsemat_in%nvctrp
  sparsemat_out%isvctr = sparsemat_in%isvctr
  sparsemat_out%nfvctr = sparsemat_in%nfvctr
  sparsemat_out%nfvctrp = sparsemat_in%nfvctrp
  sparsemat_out%isfvctr = sparsemat_in%isfvctr
  sparsemat_out%parallel_compression = sparsemat_in%parallel_compression

  if(associated(sparsemat_out%nvctr_par)) then
     call f_free_ptr(sparsemat_out%nvctr_par)
  end if
  if(associated(sparsemat_in%nvctr_par)) then
     iis1=lbound(sparsemat_in%nvctr_par,1)
     iie1=ubound(sparsemat_in%nvctr_par,1)
     sparsemat_out%nvctr_par=f_malloc_ptr((/iis1.to.iie1/),id='sparsemat_out%nvctr_par')
     do i1=iis1,iie1
        sparsemat_out%nvctr_par(i1) = sparsemat_in%nvctr_par(i1)
     end do
  end if
  if(associated(sparsemat_out%isvctr_par)) then
     call f_free_ptr(sparsemat_out%isvctr_par)
  end if
  if(associated(sparsemat_in%isvctr_par)) then
     iis1=lbound(sparsemat_in%isvctr_par,1)
     iie1=ubound(sparsemat_in%isvctr_par,1)
     sparsemat_out%isvctr_par=f_malloc_ptr((/iis1.to.iie1/),id='sparsemat_out%isvctr_par')
     do i1=iis1,iie1
        sparsemat_out%isvctr_par(i1) = sparsemat_in%isvctr_par(i1)
     end do
  end if
  if(associated(sparsemat_out%nfvctr_par)) then
     call f_free_ptr(sparsemat_out%nfvctr_par)
  end if
  if(associated(sparsemat_in%nfvctr_par)) then
     iis1=lbound(sparsemat_in%nfvctr_par,1)
     iie1=ubound(sparsemat_in%nfvctr_par,1)
     sparsemat_out%nfvctr_par=f_malloc_ptr((/iis1.to.iie1/),id='sparsemat_out%nfvctr_par')
     do i1=iis1,iie1
        sparsemat_out%nfvctr_par(i1) = sparsemat_in%nfvctr_par(i1)
     end do
  end if
  if(associated(sparsemat_out%isfvctr_par)) then
     call f_free_ptr(sparsemat_out%isfvctr_par)
  end if
  if(associated(sparsemat_in%isfvctr_par)) then
     iis1=lbound(sparsemat_in%isfvctr_par,1)
     iie1=ubound(sparsemat_in%isfvctr_par,1)
     sparsemat_out%isfvctr_par=f_malloc_ptr((/iis1.to.iie1/),id='sparsemat_out%isfvctr_par')
     do i1=iis1,iie1
        sparsemat_out%isfvctr_par(i1) = sparsemat_in%isfvctr_par(i1)
     end do
  end if

  nullify(sparsemat_out%matrix)
  nullify(sparsemat_out%matrix_compr)
  nullify(sparsemat_out%matrixp)
  nullify(sparsemat_out%matrix_comprp)

  if(associated(sparsemat_out%keyv)) then
     call f_free_ptr(sparsemat_out%keyv)
  end if
  if(associated(sparsemat_in%keyv)) then
     iis1=lbound(sparsemat_in%keyv,1)
     iie1=ubound(sparsemat_in%keyv,1)
     sparsemat_out%keyv=f_malloc_ptr(iis1.to.iie1,id='sparsemat_out%keyv')
     do i1=iis1,iie1
        sparsemat_out%keyv(i1) = sparsemat_in%keyv(i1)
     end do
  end if

  if(associated(sparsemat_out%nsegline)) then
     call f_free_ptr(sparsemat_out%nsegline)
  end if
  if(associated(sparsemat_in%nsegline)) then
     iis1=lbound(sparsemat_in%nsegline,1)
     iie1=ubound(sparsemat_in%nsegline,1)
     sparsemat_out%nsegline=f_malloc_ptr(iis1.to.iie1,id='sparsemat_out%nsegline')
     do i1=iis1,iie1
        sparsemat_out%nsegline(i1) = sparsemat_in%nsegline(i1)
     end do
  end if

  if(associated(sparsemat_out%istsegline)) then
     call f_free_ptr(sparsemat_out%istsegline)
  end if
  if(associated(sparsemat_in%istsegline)) then
     iis1=lbound(sparsemat_in%istsegline,1)
     iie1=ubound(sparsemat_in%istsegline,1)
     sparsemat_out%istsegline=f_malloc_ptr(iis1.to.iie1,id='sparsemat_out%istsegline')
     do i1=iis1,iie1
        sparsemat_out%istsegline(i1) = sparsemat_in%istsegline(i1)
     end do
  end if

  if(associated(sparsemat_out%keyg)) then
     call f_free_ptr(sparsemat_out%keyg)
  end if
  if(associated(sparsemat_in%keyg)) then
     iis1=lbound(sparsemat_in%keyg,1)
     iie1=ubound(sparsemat_in%keyg,1)
     iis2=lbound(sparsemat_in%keyg,2)
     iie2=ubound(sparsemat_in%keyg,2)
     sparsemat_out%keyg=f_malloc_ptr((/iis1.to.iie1,iis2.to.iie2/),id='sparsemat_out%keyg')
     do i1=iis1,iie1
        do i2 = iis2,iie2
           sparsemat_out%keyg(i1,i2) = sparsemat_in%keyg(i1,i2)
        end do
     end do
  end if

  if(associated(sparsemat_out%matrixindex_in_compressed_arr)) then
     call f_free_ptr(sparsemat_out%matrixindex_in_compressed_arr)
  end if
  if(associated(sparsemat_in%matrixindex_in_compressed_arr)) then
     iis1=lbound(sparsemat_in%matrixindex_in_compressed_arr,1)
     iie1=ubound(sparsemat_in%matrixindex_in_compressed_arr,1)
     iis2=lbound(sparsemat_in%matrixindex_in_compressed_arr,2)
     iie2=ubound(sparsemat_in%matrixindex_in_compressed_arr,2)
     sparsemat_out%matrixindex_in_compressed_arr=f_malloc_ptr((/iis1.to.iie1,iis2.to.iie2/),&
         id='sparsemat_out%matrixindex_in_compressed_ar')
     do i1=iis1,iie1
        do i2 = iis2,iie2
           sparsemat_out%matrixindex_in_compressed_arr(i1,i2) = sparsemat_in%matrixindex_in_compressed_arr(i1,i2)
        end do
     end do
  end if

  if(associated(sparsemat_out%orb_from_index)) then
     call f_free_ptr(sparsemat_out%orb_from_index)
  end if
  if(associated(sparsemat_in%orb_from_index)) then
     iis1=lbound(sparsemat_in%orb_from_index,1)
     iie1=ubound(sparsemat_in%orb_from_index,1)
     iis2=lbound(sparsemat_in%orb_from_index,2)
     iie2=ubound(sparsemat_in%orb_from_index,2)
     sparsemat_out%orb_from_index=f_malloc_ptr((/iis1.to.iie1,iis2.to.iie2/),id='sparsemat_out%orb_from_index')
     do i1=iis1,iie1
        do i2 = iis2,iie2
           sparsemat_out%orb_from_index(i1,i2) = sparsemat_in%orb_from_index(i1,i2)
        end do
     end do
  end if

  if(associated(sparsemat_out%matrixindex_in_compressed_fortransposed)) then
     call f_free_ptr(sparsemat_out%matrixindex_in_compressed_fortransposed)
  end if
  if(associated(sparsemat_in%matrixindex_in_compressed_fortransposed)) then
     iis1=lbound(sparsemat_in%matrixindex_in_compressed_fortransposed,1)
     iie1=ubound(sparsemat_in%matrixindex_in_compressed_fortransposed,1)
     iis2=lbound(sparsemat_in%matrixindex_in_compressed_fortransposed,2)
     iie2=ubound(sparsemat_in%matrixindex_in_compressed_fortransposed,2)
     sparsemat_out%matrixindex_in_compressed_fortransposed=f_malloc_ptr((/iis1.to.iie1,iis2.to.iie2/),&
         id='sparsemat_out%matrixindex_in_compressed_fortransposed')
     do i1=iis1,iie1
        do i2 = iis2,iie2
           sparsemat_out%matrixindex_in_compressed_fortransposed(i1,i2) = &
                sparsemat_in%matrixindex_in_compressed_fortransposed(i1,i2)
        end do
     end do
  end if


  sparsemat_out%smmm%nout = sparsemat_in%smmm%nout
  sparsemat_out%smmm%nseq = sparsemat_in%smmm%nseq
  sparsemat_out%smmm%nmaxsegk = sparsemat_in%smmm%nmaxsegk
  sparsemat_out%smmm%nmaxvalk = sparsemat_in%smmm%nmaxvalk
  sparsemat_out%smmm%nseg = sparsemat_in%smmm%nseg
  
  if(associated(sparsemat_out%smmm%ivectorindex)) then
     call f_free_ptr(sparsemat_out%smmm%ivectorindex)
  end if
  if(associated(sparsemat_in%smmm%ivectorindex)) then
     iis1=lbound(sparsemat_in%smmm%ivectorindex,1)
     iie1=ubound(sparsemat_in%smmm%ivectorindex,1)
     sparsemat_out%smmm%ivectorindex=f_malloc_ptr(iis1.to.iie1,id='sparsemat_out%smmm%ivectorindex')
     do i1=iis1,iie1
        sparsemat_out%smmm%ivectorindex(i1) = sparsemat_in%smmm%ivectorindex(i1)
     end do
  end if

  if(associated(sparsemat_out%smmm%nsegline)) then
     call f_free_ptr(sparsemat_out%smmm%nsegline)
  end if
  if(associated(sparsemat_in%smmm%nsegline)) then
     iis1=lbound(sparsemat_in%smmm%nsegline,1)
     iie1=ubound(sparsemat_in%smmm%nsegline,1)
     sparsemat_out%smmm%nsegline=f_malloc_ptr(iis1.to.iie1,id='sparsemat_out%smmm%nsegline')
     do i1=iis1,iie1
        sparsemat_out%smmm%nsegline(i1) = sparsemat_in%smmm%nsegline(i1)
     end do
  end if

  if(associated(sparsemat_out%smmm%istsegline)) then
     call f_free_ptr(sparsemat_out%smmm%istsegline)
  end if
  if(associated(sparsemat_in%smmm%istsegline)) then
     iis1=lbound(sparsemat_in%smmm%istsegline,1)
     iie1=ubound(sparsemat_in%smmm%istsegline,1)
     sparsemat_out%smmm%istsegline=f_malloc_ptr(iis1.to.iie1,id='sparsemat_out%smmm%istsegline')
     do i1=iis1,iie1
        sparsemat_out%smmm%istsegline(i1) = sparsemat_in%smmm%istsegline(i1)
     end do
  end if

  if(associated(sparsemat_out%smmm%indices_extract_sequential)) then
     call f_free_ptr(sparsemat_out%smmm%indices_extract_sequential)
  end if
  if(associated(sparsemat_in%smmm%indices_extract_sequential)) then
     iis1=lbound(sparsemat_in%smmm%indices_extract_sequential,1)
     iie1=ubound(sparsemat_in%smmm%indices_extract_sequential,1)
     sparsemat_out%smmm%indices_extract_sequential=f_malloc_ptr(iis1.to.iie1,id='sparsemat_out%smmm%indices_extract_sequential')
     do i1=iis1,iie1
        sparsemat_out%smmm%indices_extract_sequential(i1) = sparsemat_in%smmm%indices_extract_sequential(i1)
     end do
  end if

  if(associated(sparsemat_out%smmm%onedimindices)) then
     call f_free_ptr(sparsemat_out%smmm%onedimindices)
  end if
  if(associated(sparsemat_in%smmm%onedimindices)) then
     iis1=lbound(sparsemat_in%smmm%onedimindices,1)
     iie1=ubound(sparsemat_in%smmm%onedimindices,1)
     iis2=lbound(sparsemat_in%smmm%onedimindices,2)
     iie2=ubound(sparsemat_in%smmm%onedimindices,2)
     sparsemat_out%smmm%onedimindices=f_malloc_ptr((/iis1.to.iie1,iis2.to.iie2/),&
         id='sparsemat_out%smmm%onedimindices')
     do i1=iis1,iie1
        do i2 = iis2,iie2
           sparsemat_out%smmm%onedimindices(i1,i2) = &
                sparsemat_in%smmm%onedimindices(i1,i2)
        end do
     end do
  end if

  if(associated(sparsemat_out%smmm%keyg)) then
     call f_free_ptr(sparsemat_out%smmm%keyg)
  end if
  if(associated(sparsemat_in%smmm%keyg)) then
     iis1=lbound(sparsemat_in%smmm%keyg,1)
     iie1=ubound(sparsemat_in%smmm%keyg,1)
     iis2=lbound(sparsemat_in%smmm%keyg,2)
     iie2=ubound(sparsemat_in%smmm%keyg,2)
     sparsemat_out%smmm%keyg=f_malloc_ptr((/iis1.to.iie1,iis2.to.iie2/),&
         id='sparsemat_out%smmm%keyg')
     do i1=iis1,iie1
        do i2 = iis2,iie2
           sparsemat_out%smmm%keyg(i1,i2) = &
                sparsemat_in%smmm%keyg(i1,i2)
        end do
     end do
  end if


  call timing(iproc,'sparse_copy','OF')

end subroutine sparse_copy_pattern



subroutine copy_linear_matrices(linmat_in, linmat_out)
  use module_types
  implicit none

  ! Calling arguments
  type(linear_matrices),intent(in) :: linmat_in
  type(linear_matrices),intent(out) :: linmat_out

  call copy_sparse_matrix(linmat_in%s, linmat_out%s)
  call copy_sparse_matrix(linmat_in%m, linmat_out%m)
  call copy_sparse_matrix(linmat_in%l, linmat_out%l)
  call copy_sparse_matrix(linmat_in%ks, linmat_out%ks)
  call copy_sparse_matrix(linmat_in%ks_e, linmat_out%ks_e)
  call copy_matrices(linmat_in%ham_, linmat_out%ham_)
  call copy_matrices(linmat_in%ovrlp_, linmat_out%ovrlp_)
  call copy_matrices(linmat_in%kernel_, linmat_out%kernel_)

end subroutine copy_linear_matrices


subroutine copy_sparse_matrix(smat_in, smat_out)
  use sparsematrix_base, only: sparse_matrix
  use copy_utils, only: allocate_and_copy
  implicit none

  ! Calling arguments
  type(sparse_matrix),intent(in) :: smat_in
  type(sparse_matrix),intent(out) :: smat_out



  smat_out%nvctr = smat_in%nvctr
  smat_out%nseg = smat_in%nseg
  smat_out%nvctrp = smat_in%nvctrp
  smat_out%isvctr = smat_in%isvctr
  smat_out%parallel_compression = smat_in%parallel_compression
  smat_out%nfvctr = smat_in%nfvctr
  smat_out%nfvctrp = smat_in%nfvctrp
  smat_out%isfvctr = smat_in%isfvctr
  smat_out%store_index = smat_in%store_index
  smat_out%can_use_dense = smat_in%store_index

  call allocate_and_copy(smat_in%keyv, smat_out%keyv, id='smat_out%')
  call allocate_and_copy(smat_in%nsegline, smat_out%nsegline, id='smat_out%nsegline')
  call allocate_and_copy(smat_in%istsegline, smat_out%istsegline, id='smat_out%istsegline')
  call allocate_and_copy(smat_in%isvctr_par, smat_out%isvctr_par, id='smat_out%isvctr_par')
  call allocate_and_copy(smat_in%nvctr_par, smat_out%nvctr_par, id='smat_out%nvctr_par')
  call allocate_and_copy(smat_in%isfvctr_par, smat_out%isfvctr_par, id='smat_out%isfvctr_par')
  call allocate_and_copy(smat_in%nfvctr_par, smat_out%nfvctr_par, id='smat_out%nfvctr_par')

  call allocate_and_copy(smat_in%keyg, smat_out%keyg, id='smat_out%keyg')
  call allocate_and_copy(smat_in%matrixindex_in_compressed_arr, smat_out%matrixindex_in_compressed_arr, &
                         id='smat_out%matrixindex_in_compressed_arr')
  call allocate_and_copy(smat_in%orb_from_index, smat_out%orb_from_index, id='smat_out%orb_from_index')
  call allocate_and_copy(smat_in%matrixindex_in_compressed_fortransposed, smat_out%matrixindex_in_compressed_fortransposed, &
                         id='smat_out%matrixindex_in_compressed_fortransposed')

  call allocate_and_copy(smat_in%matrix_compr, smat_out%matrix_compr, id='smat_out%matrix_compr')
  call allocate_and_copy(smat_in%matrix_comprp, smat_out%matrix_comprp, id='smat_out%matrix_comprp')

  call allocate_and_copy(smat_in%matrix, smat_out%matrix, id='smat_out%matrix')
  call allocate_and_copy(smat_in%matrixp, smat_out%matrixp, id='smat_out%matrixp')


end subroutine copy_sparse_matrix


subroutine copy_sparse_matrix_matrix_multiplication(smmm_in, smmm_out)
  use sparsematrix_base, only: sparse_matrix_matrix_multiplication
  use copy_utils, only: allocate_and_copy
  implicit none

  ! Calling arguments
  type(sparse_matrix_matrix_multiplication),intent(in) :: smmm_in
  type(sparse_matrix_matrix_multiplication),intent(out) :: smmm_out
  smmm_out%nout = smmm_in%nout
  smmm_out%nseq = smmm_in%nseq
  smmm_out%nmaxsegk = smmm_in%nmaxsegk
  smmm_out%nmaxvalk = smmm_in%nmaxvalk
  smmm_out%nseg = smmm_in%nseg

  call allocate_and_copy(smmm_in%ivectorindex, smmm_out%ivectorindex, id='ivectorindex')
  call allocate_and_copy(smmm_in%nsegline, smmm_out%nsegline, id='nsegline')
  call allocate_and_copy(smmm_in%istsegline, smmm_out%istsegline, id='istsegline')
  call allocate_and_copy(smmm_in%indices_extract_sequential, smmm_out%indices_extract_sequential, id='indices_extract_sequential')
end subroutine copy_sparse_matrix_matrix_multiplication


subroutine copy_matrices(mat_in, mat_out)
  use sparsematrix_base, only: matrices
  use copy_utils, only: allocate_and_copy
  implicit none

  ! Calling arguments
  type(matrices),intent(in) :: mat_in
  type(matrices),intent(out) :: mat_out


  call allocate_and_copy(mat_in%matrix_compr, mat_out%matrix_compr, id='mat_out%matrix_compr')
  call allocate_and_copy(mat_in%matrix_comprp, mat_out%matrix_comprp, id='mat_out%matrix_comprp')

  call allocate_and_copy(mat_in%matrix, mat_out%matrix, id='mat_out%matrix')
  call allocate_and_copy(mat_in%matrixp, mat_out%matrixp, id='mat_out%matrixp')

end subroutine copy_matrices


subroutine copy_comms_linear(comms_in, comms_out)
  use communications_base, only: comms_linear
  use copy_utils, only: allocate_and_copy
  implicit none

  ! Calling arguments
  type(comms_linear),intent(in) :: comms_in
  type(comms_linear),intent(out) :: comms_out


    comms_out%nptsp_c = comms_in%nptsp_c
    comms_out%ndimpsi_c = comms_in%ndimpsi_c
    comms_out%ndimind_c = comms_in%ndimind_c
    comms_out%ndimind_f = comms_in%ndimind_f
    comms_out%nptsp_f = comms_in%nptsp_f
    comms_out%ndimpsi_f = comms_in%ndimpsi_f
    comms_out%ncomms_repartitionrho = comms_in%ncomms_repartitionrho
    comms_out%window = comms_in%window

    call allocate_and_copy(comms_in%nsendcounts_c, comms_out%nsendcounts_c, id='comms_out%nsendcounts_c')
    call allocate_and_copy(comms_in%nsenddspls_c, comms_out%nsenddspls_c, id='comms_out%nsenddspls_c')
    call allocate_and_copy(comms_in%nrecvcounts_c, comms_out%nrecvcounts_c, id='comms_out%nrecvcounts_c')
    call allocate_and_copy(comms_in%nrecvdspls_c, comms_out%nrecvdspls_c, id='comms_out%nrecvdspls_c')
    call allocate_and_copy(comms_in%isendbuf_c, comms_out%isendbuf_c, id='comms_out%isendbuf_c')
    call allocate_and_copy(comms_in%iextract_c, comms_out%iextract_c, id='comms_out%iextract_c')
    call allocate_and_copy(comms_in%iexpand_c, comms_out%iexpand_c, id='comms_out%iexpand_c')
    call allocate_and_copy(comms_in%irecvbuf_c, comms_out%irecvbuf_c, id='comms_out%irecvbuf_c')
    call allocate_and_copy(comms_in%norb_per_gridpoint_c, comms_out%norb_per_gridpoint_c, id='comms_out%norb_per_gridpoint_c')
    call allocate_and_copy(comms_in%indexrecvorbital_c, comms_out%indexrecvorbital_c, id='comms_out%indexrecvorbital_c')
    call allocate_and_copy(comms_in%nsendcounts_f, comms_out%nsendcounts_f, id='comms_out%nsendcounts_f')
    call allocate_and_copy(comms_in%nsenddspls_f, comms_out%nsenddspls_f, id='comms_out%nsenddspls_f')
    call allocate_and_copy(comms_in%nrecvcounts_f, comms_out%nrecvcounts_f, id='comms_out%nrecvcounts_f')
    call allocate_and_copy(comms_in%nrecvdspls_f, comms_out%nrecvdspls_f, id='comms_out%nrecvdspls_f')
    call allocate_and_copy(comms_in%isendbuf_f, comms_out%isendbuf_f, id='comms_out%isendbuf_f')
    call allocate_and_copy(comms_in%iextract_f, comms_out%iextract_f, id='comms_out%iextract_f')
    call allocate_and_copy(comms_in%iexpand_f, comms_out%iexpand_f, id='comms_out%iexpand_f')
    call allocate_and_copy(comms_in%irecvbuf_f, comms_out%irecvbuf_f, id='comms_out%irecvbuf_f')
    call allocate_and_copy(comms_in%norb_per_gridpoint_f, comms_out%norb_per_gridpoint_f, id='ncomms_out%orb_per_gridpoint_f')
    call allocate_and_copy(comms_in%indexrecvorbital_f, comms_out%indexrecvorbital_f, id='comms_out%indexrecvorbital_f')
    call allocate_and_copy(comms_in%isptsp_c, comms_out%isptsp_c, id='comms_out%isptsp_c')
    call allocate_and_copy(comms_in%isptsp_f, comms_out%isptsp_f, id='comms_out%isptsp_f')
    call allocate_and_copy(comms_in%nsendcounts_repartitionrho, comms_out%nsendcounts_repartitionrho, &
                           id='comms_out%nsendcounts_repartitionrho')
    call allocate_and_copy(comms_in%nrecvcounts_repartitionrho, comms_out%nrecvcounts_repartitionrho, &
                           id='comms_out%nrecvcounts_repartitionrho')
    call allocate_and_copy(comms_in%nsenddspls_repartitionrho, comms_out%nsenddspls_repartitionrho, &
                           id='comms_out%nsenddspls_repartitionrho')
    call allocate_and_copy(comms_in%nrecvdspls_repartitionrho, comms_out%nrecvdspls_repartitionrho, &
                           id='comms_out%nrecvdspls_repartitionrho')

    call allocate_and_copy(comms_in%commarr_repartitionrho, comms_out%commarr_repartitionrho, id='comms_in%commarr_repartitionrho')

    call allocate_and_copy(comms_in%psit_c, comms_out%psit_c, id='comms_out%psit_c')
    call allocate_and_copy(comms_in%psit_f, comms_out%psit_f, id='comms_out%psit_f')

end subroutine copy_comms_linear
