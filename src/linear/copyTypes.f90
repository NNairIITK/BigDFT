subroutine copy_locreg_descriptors(glrin, glrout, subname)
use module_base
use module_types
use module_interfaces, exceptThisOne => copy_locreg_descriptors
implicit none

! Calling arguments
type(locreg_descriptors),intent(in):: glrin
type(locreg_descriptors),intent(out):: glrout
character(len=*),intent(in):: subname

! Local variables
integer:: iis, iie, istat, i, iall
  
glrout%geocode = glrin%geocode
glrout%hybrid_on = glrin%hybrid_on
glrout%ns1 = glrin%ns1
glrout%ns2 = glrin%ns2
glrout%ns3 = glrin%ns3
glrout%nsi1 = glrin%nsi1
glrout%nsi2 = glrin%nsi2
glrout%nsi3 = glrin%nsi3
glrout%Localnorb = glrin%Localnorb

glrout%outofzone(1) = glrin%outofzone(1)
glrout%outofzone(2) = glrin%outofzone(2)
glrout%outofzone(3) = glrin%outofzone(3)

if(associated(glrout%projflg)) then
    iall=-product(shape(glrout%projflg))*kind(glrout%projflg)
    deallocate(glrout%projflg, stat=istat)
    call memocc(istat, iall, 'glrout%projflg', subname)
end if
   iis=lbound(glrin%projflg,1)
   iie=ubound(glrin%projflg,1)
   allocate(glrout%projflg(iis:iie), stat=istat)
   call memocc(istat, glrout%projflg, 'glrout%projflg', subname)
   do i=iis,iie
       glrout%projflg(i) = glrin%projflg(i)
   end do

call copy_grid_dimensions(glrin%d, glrout%d)
call copy_wavefunctions_descriptors(glrin%wfd, glrout%wfd, subname)
call copy_convolutions_bounds(glrin%bounds, glrout%bounds, subname)


end subroutine copy_locreg_descriptors



subroutine copy_grid_dimensions(din, dout)
use module_base
use module_types
implicit none

! Calling arguments
type(grid_dimensions),intent(in):: din
type(grid_dimensions),intent(out):: dout

dout%n1 = din%n1
dout%n2 = din%n2
dout%n3 = din%n3
dout%nfl1 = din%nfl1
dout%nfu1 = din%nfu1
dout%nfl2 = din%nfl2
dout%nfu2 = din%nfu2
dout%nfl3 = din%nfl3
dout%nfu3 = din%nfu3
dout%n1i = din%n1i
dout%n2i = din%n2i
dout%n3i = din%n3i


end subroutine copy_grid_dimensions



subroutine copy_wavefunctions_descriptors(wfdin, wfdout, subname)
use module_base
use module_types
implicit none

! Calling arguments
type(wavefunctions_descriptors),intent(in):: wfdin
type(wavefunctions_descriptors),intent(out):: wfdout
character(len=*),intent(in):: subname

! Local variables
integer:: i1, i2, iis1, iie1, iis2, iie2, istat, iall


wfdout%nvctr_c = wfdin%nvctr_c
wfdout%nvctr_f = wfdin%nvctr_f
wfdout%nseg_c = wfdin%nseg_c
wfdout%nseg_f = wfdin%nseg_f

if(associated(wfdout%keyg)) then
    iall=-product(shape(wfdout%keyg))*kind(wfdout%keyg)
    deallocate(wfdout%keyg, stat=istat)
    call memocc(istat, iall, 'wfdout%keyg', subname)
end if
iis1=lbound(wfdin%keyg,1)
iie1=ubound(wfdin%keyg,1)
iis2=lbound(wfdin%keyg,2)
iie2=ubound(wfdin%keyg,2)
allocate(wfdout%keyg(iis1:iie1,iis2:iie2), stat=istat)
call memocc(istat, wfdout%keyg, 'wfdout%keyg', subname)
do i2=iis2,iie2
    do i1=iis1,iie1
        wfdout%keyg(i1,i2) = wfdin%keyg(i1,i2)
    end do
end do
    
if(associated(wfdout%keyv)) then
    iall=-product(shape(wfdout%keyv))*kind(wfdout%keyv)
    deallocate(wfdout%keyv, stat=istat)
    call memocc(istat, iall, 'wfdout%keyv', subname)
end if
iis1=lbound(wfdin%keyv,1)
iie1=ubound(wfdin%keyv,1)
allocate(wfdout%keyv(iis1:iie1), stat=istat)
call memocc(istat, wfdout%keyv, 'wfdout%keyv', subname)
do i1=iis1,iie1
    wfdout%keyv(i1) = wfdin%keyv(i1)
end do


end subroutine copy_wavefunctions_descriptors





subroutine copy_convolutions_bounds(boundsin, boundsout, subname)
use module_base
use module_types
use module_interfaces, expectThisOne => copy_convolutions_bounds
implicit none

! Calling arguments
type(convolutions_bounds),intent(in):: boundsin
type(convolutions_bounds),intent(out):: boundsout
character(len=*),intent(in):: subname

! Local variables
integer:: iis1, iie1, iis2, iie2, iis3, iie3, i1, i2, i3, istat, iall

call copy_kinetic_bounds(boundsin%kb, boundsout%kb, subname)
call copy_shrink_bounds(boundsin%sb, boundsout%sb, subname)
call copy_grow_bounds(boundsin%gb, boundsout%gb, subname)

iis1=lbound(boundsin%ibyyzz_r,1)
iie1=ubound(boundsin%ibyyzz_r,1)
iis2=lbound(boundsin%ibyyzz_r,2)
iie2=ubound(boundsin%ibyyzz_r,2)
iis3=lbound(boundsin%ibyyzz_r,3)
iie3=ubound(boundsin%ibyyzz_r,3)

if(associated(boundsout%ibyyzz_r)) then
    iall=-product(shape(boundsout%ibyyzz_r))*kind(boundsout%ibyyzz_r)
    deallocate(boundsout%ibyyzz_r, stat=istat)
    call memocc(istat, iall, 'boundsout%ibyyzz_r', subname)
end if
allocate(boundsout%ibyyzz_r(iis1:iie1,iis2:iie2,iis3:iie3), stat=istat)
call memocc(istat, boundsout%ibyyzz_r, 'boundsout%ibyyzz_r', subname)
do i3=iis3,iie3
    do i2=iis2,iie2
        do i1=iis1,iie1
            boundsout%ibyyzz_r(i1,i2,i3) = boundsin%ibyyzz_r(i1,i2,i3)
        end do
    end do
end do

end subroutine copy_convolutions_bounds



subroutine copy_kinetic_bounds(kbin, kbout, subname)
use module_base
use module_types
implicit none

! Calling arguments
type(kinetic_bounds),intent(in):: kbin
type(kinetic_bounds),intent(out):: kbout
character(len=*),intent(in):: subname

! Local variables
integer:: iis1, iie1, iis2, iie2, iis3, iie3, i1, i2, i3, istat, iall

if(associated(kbout%ibyz_c)) then
    iall=-product(shape(kbout%ibyz_c))*kind(kbout%ibyz_c)
    deallocate(kbout%ibyz_c, stat=istat)
    call memocc(istat, iall, 'kbout%ibyz_c', subname)
end if
iis1=lbound(kbin%ibyz_c,1)
iie1=ubound(kbin%ibyz_c,1)
iis2=lbound(kbin%ibyz_c,2)
iie2=ubound(kbin%ibyz_c,2)
iis3=lbound(kbin%ibyz_c,3)
iie3=ubound(kbin%ibyz_c,3)
allocate(kbout%ibyz_c(iis1:iie1,iis2:iie2,iis3:iie3), stat=istat)
call memocc(istat, kbout%ibyz_c, 'kbout%ibyz_c', subname)
do i3=iis3,iie3
    do i2=iis2,iie2
        do i1=iis1,iie1
            kbout%ibyz_c(i1,i2,i3) = kbin%ibyz_c(i1,i2,i3)
        end do
    end do
end do


if(associated(kbout%ibxz_c)) then
    iall=-product(shape(kbout%ibxz_c))*kind(kbout%ibxz_c)
    deallocate(kbout%ibxz_c, stat=istat)
    call memocc(istat, iall, 'kbout%ibxz_c', subname)
end if
iis1=lbound(kbin%ibxz_c,1)
iie1=ubound(kbin%ibxz_c,1)
iis2=lbound(kbin%ibxz_c,2)
iie2=ubound(kbin%ibxz_c,2)
iis3=lbound(kbin%ibxz_c,3)
iie3=ubound(kbin%ibxz_c,3)
allocate(kbout%ibxz_c(iis1:iie1,iis2:iie2,iis3:iie3), stat=istat)
call memocc(istat, kbout%ibxz_c, 'kbout%ibxz_c', subname)
do i3=iis3,iie3
    do i2=iis2,iie2
        do i1=iis1,iie1
            kbout%ibxz_c(i1,i2,i3) = kbin%ibxz_c(i1,i2,i3)
        end do
    end do
end do


if(associated(kbout%ibxy_c)) then
    iall=-product(shape(kbout%ibxy_c))*kind(kbout%ibxy_c)
    deallocate(kbout%ibxy_c, stat=istat)
    call memocc(istat, iall, 'kbout%ibxy_c', subname)
end if
iis1=lbound(kbin%ibxy_c,1)
iie1=ubound(kbin%ibxy_c,1)
iis2=lbound(kbin%ibxy_c,2)
iie2=ubound(kbin%ibxy_c,2)
iis3=lbound(kbin%ibxy_c,3)
iie3=ubound(kbin%ibxy_c,3)
allocate(kbout%ibxy_c(iis1:iie1,iis2:iie2,iis3:iie3), stat=istat)
call memocc(istat, kbout%ibxy_c, 'kbout%ibxy_c', subname)
do i3=iis3,iie3
    do i2=iis2,iie2
        do i1=iis1,iie1
            kbout%ibxy_c(i1,i2,i3) = kbin%ibxy_c(i1,i2,i3)
        end do
    end do
end do


if(associated(kbout%ibyz_f)) then
    iall=-product(shape(kbout%ibyz_f))*kind(kbout%ibyz_f)
    deallocate(kbout%ibyz_f, stat=istat)
    call memocc(istat, iall, 'kbout%ibyz_f', subname)
end if
iis1=lbound(kbin%ibyz_f,1)
iie1=ubound(kbin%ibyz_f,1)
iis2=lbound(kbin%ibyz_f,2)
iie2=ubound(kbin%ibyz_f,2)
iis3=lbound(kbin%ibyz_f,3)
iie3=ubound(kbin%ibyz_f,3)
allocate(kbout%ibyz_f(iis1:iie1,iis2:iie2,iis3:iie3), stat=istat)
call memocc(istat, kbout%ibyz_f, 'kbout%ibyz_f', subname)
do i3=iis3,iie3
    do i2=iis2,iie2
        do i1=iis1,iie1
            kbout%ibyz_f(i1,i2,i3) = kbin%ibyz_f(i1,i2,i3)
        end do
    end do
end do


if(associated(kbout%ibxz_f)) then
    iall=-product(shape(kbout%ibxz_f))*kind(kbout%ibxz_f)
    deallocate(kbout%ibxz_f, stat=istat)
    call memocc(istat, iall, 'kbout%ibxz_f', subname)
end if
iis1=lbound(kbin%ibxz_f,1)
iie1=ubound(kbin%ibxz_f,1)
iis2=lbound(kbin%ibxz_f,2)
iie2=ubound(kbin%ibxz_f,2)
iis3=lbound(kbin%ibxz_f,3)
iie3=ubound(kbin%ibxz_f,3)
allocate(kbout%ibxz_f(iis1:iie1,iis2:iie2,iis3:iie3), stat=istat)
call memocc(istat, kbout%ibxz_f, 'kbout%ibxz_f', subname)
do i3=iis3,iie3
    do i2=iis2,iie2
        do i1=iis1,iie1
            kbout%ibxz_f(i1,i2,i3) = kbin%ibxz_f(i1,i2,i3)
        end do
    end do
end do


if(associated(kbout%ibxy_f)) then
    iall=-product(shape(kbout%ibxy_f))*kind(kbout%ibxy_f)
    deallocate(kbout%ibxy_f, stat=istat)
    call memocc(istat, iall, 'kbout%ibxy_f', subname)
end if
iis1=lbound(kbin%ibxy_f,1)
iie1=ubound(kbin%ibxy_f,1)
iis2=lbound(kbin%ibxy_f,2)
iie2=ubound(kbin%ibxy_f,2)
iis3=lbound(kbin%ibxy_f,3)
iie3=ubound(kbin%ibxy_f,3)
allocate(kbout%ibxy_f(iis1:iie1,iis2:iie2,iis3:iie3), stat=istat)
call memocc(istat, kbout%ibxy_f, 'kbout%ibxy_f', subname)
do i3=iis3,iie3
    do i2=iis2,iie2
        do i1=iis1,iie1
            kbout%ibxy_f(i1,i2,i3) = kbin%ibxy_f(i1,i2,i3)
        end do
    end do
end do


end subroutine copy_kinetic_bounds




subroutine copy_shrink_bounds(sbin, sbout, subname)
use module_base
use module_types
implicit none

! Calling arguments
type(shrink_bounds),intent(in):: sbin
type(shrink_bounds),intent(out):: sbout
character(len=*),intent(in):: subname

! Local variables
integer:: iis1, iie1, iis2, iie2, iis3, iie3, i1, i2, i3, istat, iall


if(associated(sbout%ibzzx_c)) then
    iall=-product(shape(sbout%ibzzx_c))*kind(sbout%ibzzx_c)
    deallocate(sbout%ibzzx_c, stat=istat)
    call memocc(istat, iall, 'sbout%ibzzx_c', subname)
end if
iis1=lbound(sbin%ibzzx_c,1)
iie1=ubound(sbin%ibzzx_c,1)
iis2=lbound(sbin%ibzzx_c,2)
iie2=ubound(sbin%ibzzx_c,2)
iis3=lbound(sbin%ibzzx_c,3)
iie3=ubound(sbin%ibzzx_c,3)
allocate(sbout%ibzzx_c(iis1:iie1,iis2:iie2,iis3:iie3), stat=istat)
call memocc(istat, sbout%ibzzx_c, 'sbout%ibzzx_c', subname)
do i3=iis3,iie3
    do i2=iis2,iie2
        do i1=iis1,iie1
            sbout%ibzzx_c(i1,i2,i3) = sbin%ibzzx_c(i1,i2,i3)
        end do
    end do
end do


if(associated(sbout%ibyyzz_c)) then
    iall=-product(shape(sbout%ibyyzz_c))*kind(sbout%ibyyzz_c)
    deallocate(sbout%ibyyzz_c, stat=istat)
    call memocc(istat, iall, 'sbout%ibyyzz_c', subname)
end if
iis1=lbound(sbin%ibyyzz_c,1)
iie1=ubound(sbin%ibyyzz_c,1)
iis2=lbound(sbin%ibyyzz_c,2)
iie2=ubound(sbin%ibyyzz_c,2)
iis3=lbound(sbin%ibyyzz_c,3)
iie3=ubound(sbin%ibyyzz_c,3)
allocate(sbout%ibyyzz_c(iis1:iie1,iis2:iie2,iis3:iie3), stat=istat)
call memocc(istat, sbout%ibyyzz_c, 'sbout%ibyyzz_c', subname)
do i3=iis3,iie3
    do i2=iis2,iie2
        do i1=iis1,iie1
            sbout%ibyyzz_c(i1,i2,i3) = sbin%ibyyzz_c(i1,i2,i3)
        end do
    end do
end do


if(associated(sbout%ibxy_ff)) then
    iall=-product(shape(sbout%ibxy_ff))*kind(sbout%ibxy_ff)
    deallocate(sbout%ibxy_ff, stat=istat)
    call memocc(istat, iall, 'sbout%ibxy_ff', subname)
end if
iis1=lbound(sbin%ibxy_ff,1)
iie1=ubound(sbin%ibxy_ff,1)
iis2=lbound(sbin%ibxy_ff,2)
iie2=ubound(sbin%ibxy_ff,2)
iis3=lbound(sbin%ibxy_ff,3)
iie3=ubound(sbin%ibxy_ff,3)
allocate(sbout%ibxy_ff(iis1:iie1,iis2:iie2,iis3:iie3), stat=istat)
call memocc(istat, sbout%ibxy_ff, 'sbout%ibxy_ff', subname)
do i3=iis3,iie3
    do i2=iis2,iie2
        do i1=iis1,iie1
            sbout%ibxy_ff(i1,i2,i3) = sbin%ibxy_ff(i1,i2,i3)
        end do
    end do
end do


if(associated(sbout%ibzzx_f)) then
    iall=-product(shape(sbout%ibzzx_f))*kind(sbout%ibzzx_f)
    deallocate(sbout%ibzzx_f, stat=istat)
    call memocc(istat, iall, 'sbout%ibzzx_f', subname)
end if
iis1=lbound(sbin%ibzzx_f,1)
iie1=ubound(sbin%ibzzx_f,1)
iis2=lbound(sbin%ibzzx_f,2)
iie2=ubound(sbin%ibzzx_f,2)
iis3=lbound(sbin%ibzzx_f,3)
iie3=ubound(sbin%ibzzx_f,3)
allocate(sbout%ibzzx_f(iis1:iie1,iis2:iie2,iis3:iie3), stat=istat)
call memocc(istat, sbout%ibzzx_f, 'sbout%ibzzx_f', subname)
do i3=iis3,iie3
    do i2=iis2,iie2
        do i1=iis1,iie1
            sbout%ibzzx_f(i1,i2,i3) = sbin%ibzzx_f(i1,i2,i3)
        end do
    end do
end do


if(associated(sbout%ibyyzz_f)) then
    iall=-product(shape(sbout%ibyyzz_f))*kind(sbout%ibyyzz_f)
    deallocate(sbout%ibyyzz_f, stat=istat)
    call memocc(istat, iall, 'sbout%ibyyzz_f', subname)
end if
iis1=lbound(sbin%ibyyzz_f,1)
iie1=ubound(sbin%ibyyzz_f,1)
iis2=lbound(sbin%ibyyzz_f,2)
iie2=ubound(sbin%ibyyzz_f,2)
iis3=lbound(sbin%ibyyzz_f,3)
iie3=ubound(sbin%ibyyzz_f,3)
allocate(sbout%ibyyzz_f(iis1:iie1,iis2:iie2,iis3:iie3), stat=istat)
call memocc(istat, sbout%ibyyzz_f, 'sbout%ibyyzz_f', subname)
do i3=iis3,iie3
    do i2=iis2,iie2
        do i1=iis1,iie1
            sbout%ibyyzz_f(i1,i2,i3) = sbin%ibyyzz_f(i1,i2,i3)
        end do
    end do
end do



end subroutine copy_shrink_bounds




subroutine copy_grow_bounds(gbin, gbout, subname)
use module_base
use module_types
implicit none

! Calling arguments
type(grow_bounds),intent(in):: gbin
type(grow_bounds),intent(out):: gbout
character(len=*),intent(in):: subname

! Local variables
integer:: iis1, iie1, iis2, iie2, iis3, iie3, i1, i2, i3, istat, iall


if(associated(gbout%ibzxx_c)) then
    iall=-product(shape(gbout%ibzxx_c))*kind(gbout%ibzxx_c)
    deallocate(gbout%ibzxx_c, stat=istat)
    call memocc(istat, iall, 'gbout%ibzxx_c', subname)
end if
iis1=lbound(gbin%ibzxx_c,1)
iie1=ubound(gbin%ibzxx_c,1)
iis2=lbound(gbin%ibzxx_c,2)
iie2=ubound(gbin%ibzxx_c,2)
iis3=lbound(gbin%ibzxx_c,3)
iie3=ubound(gbin%ibzxx_c,3)
allocate(gbout%ibzxx_c(iis1:iie1,iis2:iie2,iis3:iie3), stat=istat)
call memocc(istat, gbout%ibzxx_c, 'gbout%ibzxx_c', subname)
do i3=iis3,iie3
    do i2=iis2,iie2
        do i1=iis1,iie1
            gbout%ibzxx_c(i1,i2,i3) = gbin%ibzxx_c(i1,i2,i3)
        end do
    end do
end do


if(associated(gbout%ibxxyy_c)) then
    iall=-product(shape(gbout%ibxxyy_c))*kind(gbout%ibxxyy_c)
    deallocate(gbout%ibxxyy_c, stat=istat)
    call memocc(istat, iall, 'gbout%ibxxyy_c', subname)
end if
iis1=lbound(gbin%ibxxyy_c,1)
iie1=ubound(gbin%ibxxyy_c,1)
iis2=lbound(gbin%ibxxyy_c,2)
iie2=ubound(gbin%ibxxyy_c,2)
iis3=lbound(gbin%ibxxyy_c,3)
iie3=ubound(gbin%ibxxyy_c,3)
allocate(gbout%ibxxyy_c(iis1:iie1,iis2:iie2,iis3:iie3), stat=istat)
call memocc(istat, gbout%ibxxyy_c, 'gbout%ibxxyy_c', subname)
do i3=iis3,iie3
    do i2=iis2,iie2
        do i1=iis1,iie1
            gbout%ibxxyy_c(i1,i2,i3) = gbin%ibxxyy_c(i1,i2,i3)
        end do
    end do
end do


if(associated(gbout%ibyz_ff)) then
    iall=-product(shape(gbout%ibyz_ff))*kind(gbout%ibyz_ff)
    deallocate(gbout%ibyz_ff, stat=istat)
    call memocc(istat, iall, 'gbout%ibyz_ff', subname)
end if
iis1=lbound(gbin%ibyz_ff,1)
iie1=ubound(gbin%ibyz_ff,1)
iis2=lbound(gbin%ibyz_ff,2)
iie2=ubound(gbin%ibyz_ff,2)
iis3=lbound(gbin%ibyz_ff,3)
iie3=ubound(gbin%ibyz_ff,3)
allocate(gbout%ibyz_ff(iis1:iie1,iis2:iie2,iis3:iie3), stat=istat)
call memocc(istat, gbout%ibyz_ff, 'gbout%ibyz_ff', subname)
do i3=iis3,iie3
    do i2=iis2,iie2
        do i1=iis1,iie1
            gbout%ibyz_ff(i1,i2,i3) = gbin%ibyz_ff(i1,i2,i3)
        end do
    end do
end do



if(associated(gbout%ibzxx_f)) then
    iall=-product(shape(gbout%ibzxx_f))*kind(gbout%ibzxx_f)
    deallocate(gbout%ibzxx_f, stat=istat)
    call memocc(istat, iall, 'gbout%ibzxx_f', subname)
end if
iis1=lbound(gbin%ibzxx_f,1)
iie1=ubound(gbin%ibzxx_f,1)
iis2=lbound(gbin%ibzxx_f,2)
iie2=ubound(gbin%ibzxx_f,2)
iis3=lbound(gbin%ibzxx_f,3)
iie3=ubound(gbin%ibzxx_f,3)
allocate(gbout%ibzxx_f(iis1:iie1,iis2:iie2,iis3:iie3), stat=istat)
call memocc(istat, gbout%ibzxx_f, 'gbout%ibzxx_f', subname)
do i3=iis3,iie3
    do i2=iis2,iie2
        do i1=iis1,iie1
            gbout%ibzxx_f(i1,i2,i3) = gbin%ibzxx_f(i1,i2,i3)
        end do
    end do
end do


if(associated(gbout%ibxxyy_f)) then
    iall=-product(shape(gbout%ibxxyy_f))*kind(gbout%ibxxyy_f)
    deallocate(gbout%ibxxyy_f, stat=istat)
    call memocc(istat, iall, 'gbout%ibxxyy_f', subname)
end if
iis1=lbound(gbin%ibxxyy_f,1)
iie1=ubound(gbin%ibxxyy_f,1)
iis2=lbound(gbin%ibxxyy_f,2)
iie2=ubound(gbin%ibxxyy_f,2)
iis3=lbound(gbin%ibxxyy_f,3)
iie3=ubound(gbin%ibxxyy_f,3)
allocate(gbout%ibxxyy_f(iis1:iie1,iis2:iie2,iis3:iie3), stat=istat)
call memocc(istat, gbout%ibxxyy_f, 'gbout%ibxxyy_f', subname)
do i3=iis3,iie3
    do i2=iis2,iie2
        do i1=iis1,iie1
            gbout%ibxxyy_f(i1,i2,i3) = gbin%ibxxyy_f(i1,i2,i3)
        end do
    end do
end do


end subroutine copy_grow_bounds



subroutine copy_nonlocal_psp_descriptors(nlpspin, nlpspout, subname)
use module_base
use module_types
implicit none

! Calling arguments
type(nonlocal_psp_descriptors),intent(in):: nlpspin
type(nonlocal_psp_descriptors),intent(out):: nlpspout
character(len=*),intent(in):: subname

! Local variables
integer:: iis1, iie1, iis2, iie2, iis3, iie3, i1, i2, i3, istat, iall


nlpspout%nproj = nlpspin%nproj
nlpspout%nprojel = nlpspin%nprojel


if(associated(nlpspout%nvctr_p)) then
    iall=-product(shape(nlpspout%nvctr_p))*kind(nlpspout%nvctr_p)
    deallocate(nlpspout%nvctr_p, stat=istat)
    call memocc(istat, iall, 'nlpspout%nvctr_p', subname)
end if
iis1=lbound(nlpspin%nvctr_p,1)
iie1=ubound(nlpspin%nvctr_p,1)
allocate(nlpspout%nvctr_p(iis1:iie1), stat=istat)
call memocc(istat, nlpspout%nvctr_p, 'nlpspout%nvctr_p', subname)
do i1=iis1,iie1
    nlpspout%nvctr_p(i1) = nlpspin%nvctr_p(i1)
end do


if(associated(nlpspout%nseg_p)) then
    iall=-product(shape(nlpspout%nseg_p))*kind(nlpspout%nseg_p)
    deallocate(nlpspout%nseg_p, stat=istat)
    call memocc(istat, iall, 'nlpspout%nseg_p', subname)
end if
iis1=lbound(nlpspin%nseg_p,1)
iie1=ubound(nlpspin%nseg_p,1)
allocate(nlpspout%nseg_p(iis1:iie1), stat=istat)
call memocc(istat, nlpspout%nseg_p, 'nlpspout%nseg_p', subname)
do i1=iis1,iie1
    nlpspout%nseg_p(i1) = nlpspin%nseg_p(i1)
end do


if(associated(nlpspout%keyv_p)) then
    iall=-product(shape(nlpspout%keyv_p))*kind(nlpspout%keyv_p)
    deallocate(nlpspout%keyv_p, stat=istat)
    call memocc(istat, iall, 'nlpspout%keyv_p', subname)
end if
iis1=lbound(nlpspin%keyv_p,1)
iie1=ubound(nlpspin%keyv_p,1)
allocate(nlpspout%keyv_p(iis1:iie1), stat=istat)
call memocc(istat, nlpspout%keyv_p, 'nlpspout%keyv_p', subname)
do i1=iis1,iie1
    nlpspout%keyv_p(i1) = nlpspin%keyv_p(i1)
end do


if(associated(nlpspout%keyg_p)) then
    iall=-product(shape(nlpspout%keyg_p))*kind(nlpspout%keyg_p)
    deallocate(nlpspout%keyg_p, stat=istat)
    call memocc(istat, iall, 'nlpspout%keyg_p', subname)
end if
iis1=lbound(nlpspin%keyg_p,1)
iie1=ubound(nlpspin%keyg_p,1)
iis2=lbound(nlpspin%keyg_p,2)
iie2=ubound(nlpspin%keyg_p,2)
allocate(nlpspout%keyg_p(iis1:iie1,iis2:iie2), stat=istat)
call memocc(istat, nlpspin%keyg_p, 'nlpspin%keyg_p', subname)
do i2=iis2,iie2
    do i1=iis1,iie1
        nlpspout%keyg_p(i1,i2) = nlpspin%keyg_p(i1,i2)
    end do
end do


if(associated(nlpspout%nboxp_c)) then
    iall=-product(shape(nlpspout%nboxp_c))*kind(nlpspout%nboxp_c)
    deallocate(nlpspout%nboxp_c, stat=istat)
    call memocc(istat, iall, 'nlpspout%nboxp_c', subname)
end if
iis1=lbound(nlpspin%nboxp_c,1)
iie1=ubound(nlpspin%nboxp_c,1)
iis2=lbound(nlpspin%nboxp_c,2)
iie2=ubound(nlpspin%nboxp_c,2)
iis3=lbound(nlpspin%nboxp_c,3)
iie3=ubound(nlpspin%nboxp_c,3)
allocate(nlpspout%nboxp_c(iis1:iie1,iis2:iie2,iis3:iie3), stat=istat)
call memocc(istat, nlpspout%nboxp_c, 'nlpspout%nboxp_c', subname)
do i3=iis3,iie3
    do i2=iis2,iie2
        do i1=iis1,iie1
            nlpspout%nboxp_c(i1,i2,i3) = nlpspin%nboxp_c(i1,i2,i3)
        end do
    end do
end do


if(associated(nlpspout%nboxp_f)) then
    iall=-product(shape(nlpspout%nboxp_f))*kind(nlpspout%nboxp_f)
    deallocate(nlpspout%nboxp_f, stat=istat)
    call memocc(istat, iall, 'nlpspout%nboxp_f', subname)
end if
iis1=lbound(nlpspin%nboxp_f,1)
iie1=ubound(nlpspin%nboxp_f,1)
iis2=lbound(nlpspin%nboxp_f,2)
iie2=ubound(nlpspin%nboxp_f,2)
iis3=lbound(nlpspin%nboxp_f,3)
iie3=ubound(nlpspin%nboxp_f,3)
allocate(nlpspout%nboxp_f(iis1:iie1,iis2:iie2,iis3:iie3), stat=istat)
call memocc(istat, nlpspout%nboxp_f, 'nlpspout%nboxp_f', subname)
do i3=iis3,iie3
    do i2=iis2,iie2
        do i1=iis1,iie1
            nlpspout%nboxp_f(i1,i2,i3) = nlpspin%nboxp_f(i1,i2,i3)
        end do
    end do
end do


end subroutine copy_nonlocal_psp_descriptors



subroutine copy_orbitals_data(orbsin, orbsout, subname)
use module_base
use module_types
implicit none

! Calling arguments
type(orbitals_data),intent(in):: orbsin
type(orbitals_data),intent(out):: orbsout
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
orbsout%npsidim = orbsin%npsidim
orbsout%nkpts = orbsin%nkpts
orbsout%nkptsp = orbsin%nkptsp
orbsout%iskpts = orbsin%iskpts
orbsout%efermi = orbsin%efermi

if(associated(orbsout%norb_par)) then
    iall=-product(shape(orbsout%norb_par))*kind(orbsout%norb_par)
    deallocate(orbsout%norb_par, stat=istat)
    call memocc(istat, iall, 'orbsout%norb_par', subname)
end if
iis1=lbound(orbsin%norb_par,1)
iie1=ubound(orbsin%norb_par,1)
allocate(orbsout%norb_par(iis1:iie1), stat=istat)
call memocc(istat, orbsout%norb_par, 'orbsout%norb_par', subname)
do i1=iis1,iie1
    orbsout%norb_par(i1) = orbsin%norb_par(i1)
end do


if(associated(orbsout%iokpt)) then
    iall=-product(shape(orbsout%iokpt))*kind(orbsout%iokpt)
    deallocate(orbsout%iokpt, stat=istat)
    call memocc(istat, iall, 'orbsout%iokpt', subname)
end if
iis1=lbound(orbsin%iokpt,1)
iie1=ubound(orbsin%iokpt,1)
allocate(orbsout%iokpt(iis1:iie1), stat=istat)
call memocc(istat, orbsout%iokpt, 'orbsout%iokpt', subname)
do i1=iis1,iie1
    orbsout%iokpt(i1) = orbsin%iokpt(i1)
end do


if(associated(orbsout%ikptproc)) then
    iall=-product(shape(orbsout%ikptproc))*kind(orbsout%ikptproc)
    deallocate(orbsout%ikptproc, stat=istat)
    call memocc(istat, iall, 'orbsout%ikptproc', subname)
end if
iis1=lbound(orbsin%ikptproc,1)
iie1=ubound(orbsin%ikptproc,1)
allocate(orbsout%ikptproc(iis1:iie1), stat=istat)
call memocc(istat, orbsout%ikptproc, 'orbsout%ikptproc', subname)
do i1=iis1,iie1
    orbsout%ikptproc(i1) = orbsin%ikptproc(i1)
end do

if(associated(orbsout%inwhichlocreg)) then
    iall=-product(shape(orbsout%inwhichlocreg))*kind(orbsout%inwhichlocreg)
    deallocate(orbsout%inwhichlocreg, stat=istat)
    call memocc(istat, iall, 'orbsout%inwhichlocreg', subname)
end if
iis1=lbound(orbsin%inwhichlocreg,1)
iie1=ubound(orbsin%inwhichlocreg,1)
allocate(orbsout%inwhichlocreg(iis1:iie1), stat=istat)
call memocc(istat, orbsout%inwhichlocreg, 'orbsout%inwhichlocreg', subname)
do i1=iis1,iie1
    orbsout%inwhichlocreg(i1) = orbsin%inwhichlocreg(i1)
end do

if(associated(orbsout%inWhichLocregP)) then
    iall=-product(shape(orbsout%inWhichLocregP))*kind(orbsout%inWhichLocregP)
    deallocate(orbsout%inWhichLocregP, stat=istat)
    call memocc(istat, iall, 'orbsout%inWhichLocregP', subname)
end if
iis1=lbound(orbsin%inWhichLocregP,1)
iie1=ubound(orbsin%inWhichLocregP,1)
allocate(orbsout%inWhichLocregP(iis1:iie1), stat=istat)
call memocc(istat, orbsout%inWhichLocregP, 'orbsout%inWhichLocregP', subname)
do i1=iis1,iie1
    orbsout%inWhichLocregP(i1) = orbsin%inWhichLocregP(i1)
end do

if(associated(orbsout%onWhichMPI)) then
    iall=-product(shape(orbsout%onWhichMPI))*kind(orbsout%onWhichMPI)
    deallocate(orbsout%onWhichMPI, stat=istat)
    call memocc(istat, iall, 'orbsout%onWhichMPI', subname)
end if
iis1=lbound(orbsin%onWhichMPI,1)
iie1=ubound(orbsin%onWhichMPI,1)
allocate(orbsout%onWhichMPI(iis1:iie1), stat=istat)
call memocc(istat, orbsout%onWhichMPI, 'orbsout%onWhichMPI', subname)
do i1=iis1,iie1
    orbsout%onWhichMPI(i1) = orbsin%onWhichMPI(i1)
end do


if(associated(orbsout%isorb_par)) then
    iall=-product(shape(orbsout%isorb_par))*kind(orbsout%isorb_par)
    deallocate(orbsout%isorb_par, stat=istat)
    call memocc(istat, iall, 'orbsout%isorb_par', subname)
end if
   iis1=lbound(orbsin%isorb_par,1)
   iie1=ubound(orbsin%isorb_par,1)
   allocate(orbsout%isorb_par(iis1:iie1), stat=istat)
   call memocc(istat, orbsout%isorb_par, 'orbsout%isorb_par', subname)
   do i1=iis1,iie1
       orbsout%isorb_par(i1) = orbsin%isorb_par(i1)
   end do


if(associated(orbsout%eval)) then
    iall=-product(shape(orbsout%eval))*kind(orbsout%eval)
    deallocate(orbsout%eval, stat=istat)
    call memocc(istat, iall, 'orbsout%eval', subname)
end if
iis1=lbound(orbsin%eval,1)
iie1=ubound(orbsin%eval,1)
allocate(orbsout%eval(iis1:iie1), stat=istat)
call memocc(istat, orbsout%eval, 'orbsout%eval', subname)
do i1=iis1,iie1
    orbsout%eval(i1) = orbsin%eval(i1)
end do


if(associated(orbsout%occup)) then
    iall=-product(shape(orbsout%occup))*kind(orbsout%occup)
    deallocate(orbsout%occup, stat=istat)
    call memocc(istat, iall, 'orbsout%occup', subname)
end if
iis1=lbound(orbsin%occup,1)
iie1=ubound(orbsin%occup,1)
allocate(orbsout%occup(iis1:iie1), stat=istat)
call memocc(istat, orbsout%occup, 'orbsout%occup', subname)
do i1=iis1,iie1
    orbsout%occup(i1) = orbsin%occup(i1)
end do


if(associated(orbsout%spinsgn)) then
    iall=-product(shape(orbsout%spinsgn))*kind(orbsout%spinsgn)
    deallocate(orbsout%spinsgn, stat=istat)
    call memocc(istat, iall, 'orbsout%spinsgn', subname)
end if
iis1=lbound(orbsin%spinsgn,1)
iie1=ubound(orbsin%spinsgn,1)
allocate(orbsout%spinsgn(iis1:iie1), stat=istat)
call memocc(istat, orbsout%spinsgn, 'orbsout%spinsgn', subname)
do i1=iis1,iie1
    orbsout%spinsgn(i1) = orbsin%spinsgn(i1)
end do


if(associated(orbsout%kwgts)) then
    iall=-product(shape(orbsout%kwgts))*kind(orbsout%kwgts)
    deallocate(orbsout%kwgts, stat=istat)
    call memocc(istat, iall, 'orbsout%kwgts', subname)
end if
iis1=lbound(orbsin%kwgts,1)
iie1=ubound(orbsin%kwgts,1)
allocate(orbsout%kwgts(iis1:iie1), stat=istat)
call memocc(istat, orbsout%kwgts, 'orbsout%kwgts', subname)
do i1=iis1,iie1
    orbsout%kwgts(i1) = orbsin%kwgts(i1)
end do


if(associated(orbsout%kpts)) then
    iall=-product(shape(orbsout%kpts))*kind(orbsout%kpts)
    deallocate(orbsout%kpts, stat=istat)
    call memocc(istat, iall, 'orbsout%kpts', subname)
end if
iis1=lbound(orbsin%kpts,1)
iie1=ubound(orbsin%kpts,1)
iis2=lbound(orbsin%kpts,2)
iie2=ubound(orbsin%kpts,2)
allocate(orbsout%kpts(iis1:iie1,iis2:iie2), stat=istat)
call memocc(istat, orbsout%kpts, 'orbsout%kpts', subname)
do i2=iis2,iie2
    do i1=iis1,iie1
        orbsout%kpts(i1,i2) = orbsin%kpts(i1,i2)
    end do
end do


end subroutine copy_orbitals_data

