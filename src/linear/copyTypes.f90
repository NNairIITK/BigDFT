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
integer:: iis, iie, istat, i
  
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
integer:: i1, i2, iis1, iie1, iis2, iie2, istat


wfdout%nvctr_c = wfdin%nvctr_c
wfdout%nvctr_f = wfdin%nvctr_f
wfdout%nseg_c = wfdin%nseg_c
wfdout%nseg_f = wfdin%nseg_f

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
integer:: iis1, iie1, iis2, iie2, iis3, iie3, i1, i2, i3, istat

call copy_kinetic_bounds(boundsin%kb, boundsout%kb, subname)
call copy_shrink_bounds(boundsin%sb, boundsout%sb, subname)
call copy_grow_bounds(boundsin%gb, boundsout%gb, subname)

iis1=lbound(boundsin%ibyyzz_r,1)
iie1=ubound(boundsin%ibyyzz_r,1)
iis2=lbound(boundsin%ibyyzz_r,2)
iie2=ubound(boundsin%ibyyzz_r,2)
iis3=lbound(boundsin%ibyyzz_r,3)
iie3=ubound(boundsin%ibyyzz_r,3)

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
integer:: iis1, iie1, iis2, iie2, iis3, iie3, i1, i2, i3, istat

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
integer:: iis1, iie1, iis2, iie2, iis3, iie3, i1, i2, i3, istat


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
integer:: iis1, iie1, iis2, iie2, iis3, iie3, i1, i2, i3, istat


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
integer:: iis1, iie1, iis2, iie2, iis3, iie3, i1, i2, i3, istat


nlpspout%nproj = nlpspin%nproj
nlpspout%nprojel = nlpspin%nprojel


iis1=lbound(nlpspin%nvctr_p,1)
iie1=ubound(nlpspin%nvctr_p,1)
allocate(nlpspout%nvctr_p(iis1:iie1), stat=istat)
call memocc(istat, nlpspout%nvctr_p, 'nlpspout%nvctr_p', subname)
do i1=iis1,iie1
    nlpspout%nvctr_p(i1) = nlpspin%nvctr_p(i1)
end do


iis1=lbound(nlpspin%nseg_p,1)
iie1=ubound(nlpspin%nseg_p,1)
allocate(nlpspout%nseg_p(iis1:iie1), stat=istat)
call memocc(istat, nlpspout%nseg_p, 'nlpspout%nseg_p', subname)
do i1=iis1,iie1
    nlpspout%nseg_p(i1) = nlpspin%nseg_p(i1)
end do


iis1=lbound(nlpspin%keyv_p,1)
iie1=ubound(nlpspin%keyv_p,1)
allocate(nlpspout%keyv_p(iis1:iie1), stat=istat)
call memocc(istat, nlpspout%keyv_p, 'nlpspout%keyv_p', subname)
do i1=iis1,iie1
    nlpspout%keyv_p(i1) = nlpspin%keyv_p(i1)
end do


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
integer:: iis1, iie1, iis2, iie2, i1, i2, istat

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

iis1=lbound(orbsin%norb_par,1)
iie1=ubound(orbsin%norb_par,1)
allocate(orbsout%norb_par(iis1:iie1), stat=istat)
call memocc(istat, orbsout%norb_par, 'orbsout%norb_par', subname)
do i1=iis1,iie1
    orbsout%norb_par(i1) = orbsin%norb_par(i1)
end do


iis1=lbound(orbsin%iokpt,1)
iie1=ubound(orbsin%iokpt,1)
allocate(orbsout%iokpt(iis1:iie1), stat=istat)
call memocc(istat, orbsout%iokpt, 'orbsout%iokpt', subname)
do i1=iis1,iie1
    orbsout%iokpt(i1) = orbsin%iokpt(i1)
end do


iis1=lbound(orbsin%ikptproc,1)
iie1=ubound(orbsin%ikptproc,1)
allocate(orbsout%ikptproc(iis1:iie1), stat=istat)
call memocc(istat, orbsout%ikptproc, 'orbsout%ikptproc', subname)
do i1=iis1,iie1
    orbsout%ikptproc(i1) = orbsin%ikptproc(i1)
end do


iis1=lbound(orbsin%inwhichlocreg,1)
iie1=ubound(orbsin%inwhichlocreg,1)
allocate(orbsout%inwhichlocreg(iis1:iie1), stat=istat)
call memocc(istat, orbsout%inwhichlocreg, 'orbsout%inwhichlocreg', subname)
do i1=iis1,iie1
    orbsout%inwhichlocreg(i1) = orbsin%inwhichlocreg(i1)
end do


iis1=lbound(orbsin%inWhichLocregP,1)
iie1=ubound(orbsin%inWhichLocregP,1)
allocate(orbsout%inWhichLocregP(iis1:iie1), stat=istat)
call memocc(istat, orbsout%inWhichLocregP, 'orbsout%inWhichLocregP', subname)
do i1=iis1,iie1
    orbsout%inWhichLocregP(i1) = orbsin%inWhichLocregP(i1)
end do


iis1=lbound(orbsin%onWhichMPI,1)
iie1=ubound(orbsin%onWhichMPI,1)
allocate(orbsout%onWhichMPI(iis1:iie1), stat=istat)
call memocc(istat, orbsout%onWhichMPI, 'orbsout%onWhichMPI', subname)
do i1=iis1,iie1
    orbsout%onWhichMPI(i1) = orbsin%onWhichMPI(i1)
end do


iis1=lbound(orbsin%isorb_par,1)
iie1=ubound(orbsin%isorb_par,1)
allocate(orbsout%isorb_par(iis1:iie1), stat=istat)
call memocc(istat, orbsout%isorb_par, 'orbsout%isorb_par', subname)
do i1=iis1,iie1
    orbsout%isorb_par(i1) = orbsin%isorb_par(i1)
end do


iis1=lbound(orbsin%eval,1)
iie1=ubound(orbsin%eval,1)
allocate(orbsout%eval(iis1:iie1), stat=istat)
call memocc(istat, orbsout%eval, 'orbsout%eval', subname)
do i1=iis1,iie1
    orbsout%eval(i1) = orbsin%eval(i1)
end do


iis1=lbound(orbsin%occup,1)
iie1=ubound(orbsin%occup,1)
allocate(orbsout%occup(iis1:iie1), stat=istat)
call memocc(istat, orbsout%occup, 'orbsout%occup', subname)
do i1=iis1,iie1
    orbsout%occup(i1) = orbsin%occup(i1)
end do


iis1=lbound(orbsin%spinsgn,1)
iie1=ubound(orbsin%spinsgn,1)
allocate(orbsout%spinsgn(iis1:iie1), stat=istat)
call memocc(istat, orbsout%spinsgn, 'orbsout%spinsgn', subname)
do i1=iis1,iie1
    orbsout%spinsgn(i1) = orbsin%spinsgn(i1)
end do


iis1=lbound(orbsin%kwgts,1)
iie1=ubound(orbsin%kwgts,1)
allocate(orbsout%kwgts(iis1:iie1), stat=istat)
call memocc(istat, orbsout%kwgts, 'orbsout%kwgts', subname)
do i1=iis1,iie1
    orbsout%kwgts(i1) = orbsin%kwgts(i1)
end do


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

