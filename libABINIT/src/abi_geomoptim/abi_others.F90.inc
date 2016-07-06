! ######## Begin case ionmov==14 ######
! 
    a1 =  0.0378593198406116_dp;
    a2 =  0.102635633102435_dp;
    a3 = -0.0258678882665587_dp;
    a4 =  0.314241403071447_dp;
    a5 = -0.130144459517415_dp; 
    a6 =  0.106417700369543_dp;
    a7 = -0.00879424312851058_dp;
    a8 =  1._dp-2._dp*(a1+a2+a3+a4+a5+a6+a7);
    a9 =  a7;
    a10=  a6;
    a11=  a5;
    a12=  a4;
    a13=  a3;
    a14=  a2;
    a15=  a1;

    b1 =  0.09171915262446165_dp;
    b2 =  0.183983170005006_dp;
    b3 = -0.05653436583288827_dp; 
    b4 =  0.004914688774712854_dp;
    b5 =  0.143761127168358_dp;
    b6 =  0.328567693746804_dp;
    b7 =  0.5_dp - (b1+b2+b3+b4+b5+b6);
    b8 =  b6;
    b9 =  b5;
    b10=  b4;
    b11=  b3;
    b12=  b2;
    b13=  b1;

   acell_next(:)=acell(:)
   ucvol_next=ucvol
   rprim_next(:,:)=rprim(:,:)
   rprimd_next(:,:)=rprimd(:,:)

! step 1 of 15
!  Convert input xred (reduced coordinates) to xcart (cartesian)
   call abi_xredxcart(natom,1,rprimd,xcart,xred)
   xcart_next(:,:) = xcart(:,:) + dtion * a1 * vel(:,:);
!   Convert back to xred (reduced coordinates)
    call abi_xredxcart(natom,-1,rprimd,xcart_next,xred_next)

!step 2 of 15
!   Computation of the forces for the new positions
    iapp=-1
    if(itime>0)iapp=itime
!!$    call scfcv(acell,atindx,atindx1,cg,cpus,densymop_gs,dtefield,dtfil,dtset,ecore,&
!!$&    eigen,hdr,iapp,indsym,initialized,&
!!$&    irrzon,kg,mpi_enreg,&
!!$&    nattyp,nfftf,npwarr,nspinor,occ,&
!!$&    pawang,pawfgr,pawrad,pawrhoij,pawtab,phnons,psps,&
!!$&    pwind,pwind_alloc,pwnsfac,rec_set,resid,results_gs,rhog,rhor,rprimd,&
!!$&    scf_history,symrec,wffnew,wffnow,wvl,xred_next,xred,ylm,ylmgr)
    do iatom=1,natom
     do idim=1,3
      fcart_m(idim,iatom)=fcart(idim,iatom)/amass(iatom)
     end do
    end do
!   update of velocities and positions
    vel(:,:) = vel(:,:) + b1 * dtion * fcart_m(:,:)
    xred(:,:) = xred_next(:,:)
    xcart_next(:,:) = xcart_next(:,:) + a2 * dtion * vel(:,:)
!   Convert xcart_next to xred_next (reduced coordinates) for scfcv
    call abi_xredxcart(natom,-1,rprimd,xcart_next,xred_next)

!step 3 of 15
!   Computation of the forces for the new positions
    iapp=-1
    if(itime>0)iapp=itime
!!$    call scfcv(acell,atindx,atindx1,cg,cpus,densymop_gs,dtefield,dtfil,dtset,ecore,&
!!$&    eigen,hdr,iapp,indsym,initialized,&
!!$&    irrzon,kg,mpi_enreg,&
!!$&    nattyp,nfftf,npwarr,nspinor,occ,&
!!$&    pawang,pawfgr,pawrad,pawrhoij,pawtab,phnons,psps,&
!!$&    pwind,pwind_alloc,pwnsfac,rec_set,resid,results_gs,rhog,rhor,rprimd,&
!!$&    scf_history,symrec,wffnew,wffnow,wvl,xred_next,xred,ylm,ylmgr)
    do iatom=1,natom
     do idim=1,3
      fcart_m(idim,iatom)=fcart(idim,iatom)/amass(iatom)
     end do
    end do
!   update of velocities and positions
    vel(:,:) = vel(:,:) + b2 * dtion * fcart_m(:,:)
    xred(:,:) = xred_next(:,:)
    xcart_next(:,:) = xcart_next(:,:) + a3 * dtion * vel(:,:)
!   Convert xcart_next to xred_next (reduced coordinates) for scfcv
    call abi_xredxcart(natom,-1,rprimd,xcart_next,xred_next)    

!step 4 of 15
!   Computation of the forces for the new positions
    iapp=-1
    if(itime>0)iapp=itime
!!$    call scfcv(acell,atindx,atindx1,cg,cpus,densymop_gs,dtefield,dtfil,dtset,ecore,&
!!$&    eigen,hdr,iapp,indsym,initialized,&
!!$&    irrzon,kg,mpi_enreg,&
!!$&    nattyp,nfftf,npwarr,nspinor,occ,&
!!$&    pawang,pawfgr,pawrad,pawrhoij,pawtab,phnons,psps,&
!!$&    pwind,pwind_alloc,pwnsfac,rec_set,resid,results_gs,rhog,rhor,rprimd,&
!!$&    scf_history,symrec,wffnew,wffnow,wvl,xred_next,xred,ylm,ylmgr)
    do iatom=1,natom
     do idim=1,3
      fcart_m(idim,iatom)=fcart(idim,iatom)/amass(iatom)
     end do
    end do
!   update of velocities and positions
    vel(:,:) = vel(:,:) + b3 * dtion * fcart_m(:,:)
    xred(:,:) = xred_next(:,:)
    xcart_next(:,:) = xcart_next(:,:) + a4 * dtion * vel(:,:)
!   Convert xcart_next to xred_next (reduced coordinates) for scfcv
    call abi_xredxcart(natom,-1,rprimd,xcart_next,xred_next)    

!step 5 of 15
!   Computation of the forces for the new positions
    iapp=-1
    if(itime>0)iapp=itime
!!$    call scfcv(acell,atindx,atindx1,cg,cpus,densymop_gs,dtefield,dtfil,dtset,ecore,&
!!$&    eigen,hdr,iapp,indsym,initialized,&
!!$&    irrzon,kg,mpi_enreg,&
!!$&    nattyp,nfftf,npwarr,nspinor,occ,&
!!$&    pawang,pawfgr,pawrad,pawrhoij,pawtab,phnons,psps,&
!!$&    pwind,pwind_alloc,pwnsfac,rec_set,resid,results_gs,rhog,rhor,rprimd,&
!!$&    scf_history,symrec,wffnew,wffnow,wvl,xred_next,xred,ylm,ylmgr)
    do iatom=1,natom
     do idim=1,3
      fcart_m(idim,iatom)=fcart(idim,iatom)/amass(iatom)
     end do
    end do
!   update of velocities and positions
    vel(:,:) = vel(:,:) + b4 * dtion * fcart_m(:,:)
    xred(:,:) = xred_next(:,:)
    xcart_next(:,:) = xcart_next(:,:) + a5 * dtion * vel(:,:)
!   Convert xcart_next to xred_next (reduced coordinates) for scfcv
    call abi_xredxcart(natom,-1,rprimd,xcart_next,xred_next)
    
!step 6 of 15
!   Computation of the forces for the new positions
    iapp=-1
    if(itime>0)iapp=itime
!!$    call scfcv(acell,atindx,atindx1,cg,cpus,densymop_gs,dtefield,dtfil,dtset,ecore,&
!!$&    eigen,hdr,iapp,indsym,initialized,&
!!$&    irrzon,kg,mpi_enreg,&
!!$&    nattyp,nfftf,npwarr,nspinor,occ,&
!!$&    pawang,pawfgr,pawrad,pawrhoij,pawtab,phnons,psps,&
!!$&    pwind,pwind_alloc,pwnsfac,rec_set,resid,results_gs,rhog,rhor,rprimd,&
!!$&    scf_history,symrec,wffnew,wffnow,wvl,xred_next,xred,ylm,ylmgr)
    do iatom=1,natom
     do idim=1,3
      fcart_m(idim,iatom)=fcart(idim,iatom)/amass(iatom)
     end do
    end do
!   update of velocities and positions
    vel(:,:) = vel(:,:) + b5 * dtion * fcart_m(:,:)
    xred(:,:) = xred_next(:,:)
    xcart_next(:,:) = xcart_next(:,:) + a6 * dtion * vel(:,:)
!   Convert xcart_next to xred_next (reduced coordinates) for scfcv
    call abi_xredxcart(natom,-1,rprimd,xcart_next,xred_next)    

!step 7 of 15
!   Computation of the forces for the new positions
    iapp=-1
    if(itime>0)iapp=itime
!!$    call scfcv(acell,atindx,atindx1,cg,cpus,densymop_gs,dtefield,dtfil,dtset,ecore,&
!!$&    eigen,hdr,iapp,indsym,initialized,&
!!$&    irrzon,kg,mpi_enreg,&
!!$&    nattyp,nfftf,npwarr,nspinor,occ,&
!!$&    pawang,pawfgr,pawrad,pawrhoij,pawtab,phnons,psps,&
!!$&    pwind,pwind_alloc,pwnsfac,rec_set,resid,results_gs,rhog,rhor,rprimd,&
!!$&    scf_history,symrec,wffnew,wffnow,wvl,xred_next,xred,ylm,ylmgr)
    do iatom=1,natom
     do idim=1,3
      fcart_m(idim,iatom)=fcart(idim,iatom)/amass(iatom)
     end do
    end do
!   update of velocities and positions
    vel(:,:) = vel(:,:) + b6 * dtion * fcart_m(:,:)
    xred(:,:) = xred_next(:,:)
    xcart_next(:,:) = xcart_next(:,:) + a7 * dtion * vel(:,:)
!   Convert xcart_next to xred_next (reduced coordinates) for scfcv
    call abi_xredxcart(natom,-1,rprimd,xcart_next,xred_next)    

!step 8 of 15
!   Computation of the forces for the new positions
    iapp=-1
    if(itime>0)iapp=itime
!!$    call scfcv(acell,atindx,atindx1,cg,cpus,densymop_gs,dtefield,dtfil,dtset,ecore,&
!!$&    eigen,hdr,iapp,indsym,initialized,&
!!$&    irrzon,kg,mpi_enreg,&
!!$&    nattyp,nfftf,npwarr,nspinor,occ,&
!!$&    pawang,pawfgr,pawrad,pawrhoij,pawtab,phnons,psps,&
!!$&    pwind,pwind_alloc,pwnsfac,rec_set,resid,results_gs,rhog,rhor,rprimd,&
!!$&    scf_history,symrec,wffnew,wffnow,wvl,xred_next,xred,ylm,ylmgr)
    do iatom=1,natom
     do idim=1,3
      fcart_m(idim,iatom)=fcart(idim,iatom)/amass(iatom)
     end do
    end do
!   update of velocities and positions
    vel(:,:) = vel(:,:) + b7 * dtion * fcart_m(:,:)
    xred(:,:) = xred_next(:,:)
    xcart_next(:,:) = xcart_next(:,:) + a8 * dtion * vel(:,:)
!   Convert xcart_next to xred_next (reduced coordinates) for scfcv
    call abi_xredxcart(natom,-1,rprimd,xcart_next,xred_next)    

!step 9 of 15
!   Computation of the forces for the new positions
    iapp=-1
    if(itime>0)iapp=itime
!!$    call scfcv(acell,atindx,atindx1,cg,cpus,densymop_gs,dtefield,dtfil,dtset,ecore,&
!!$&    eigen,hdr,iapp,indsym,initialized,&
!!$&    irrzon,kg,mpi_enreg,&
!!$&    nattyp,nfftf,npwarr,nspinor,occ,&
!!$&    pawang,pawfgr,pawrad,pawrhoij,pawtab,phnons,psps,&
!!$&    pwind,pwind_alloc,pwnsfac,rec_set,resid,results_gs,rhog,rhor,rprimd,&
!!$&    scf_history,symrec,wffnew,wffnow,wvl,xred_next,xred,ylm,ylmgr)
    do iatom=1,natom
     do idim=1,3
      fcart_m(idim,iatom)=fcart(idim,iatom)/amass(iatom)
     end do
    end do
!   update of velocities and positions
    vel(:,:) = vel(:,:) + b7 * dtion * fcart_m(:,:)
    xred(:,:) = xred_next(:,:)
    xcart_next(:,:) = xcart_next(:,:) + a9 * dtion * vel(:,:)
!   Convert xcart_next to xred_next (reduced coordinates) for scfcv
    call abi_xredxcart(natom,-1,rprimd,xcart_next,xred_next)    

!step 10 of 15
!   Computation of the forces for the new positions
    iapp=-1
    if(itime>0)iapp=itime
!!$    call scfcv(acell,atindx,atindx1,cg,cpus,densymop_gs,dtefield,dtfil,dtset,ecore,&
!!$&    eigen,hdr,iapp,indsym,initialized,&
!!$&    irrzon,kg,mpi_enreg,&
!!$&    nattyp,nfftf,npwarr,nspinor,occ,&
!!$&    pawang,pawfgr,pawrad,pawrhoij,pawtab,phnons,psps,&
!!$&    pwind,pwind_alloc,pwnsfac,rec_set,resid,results_gs,rhog,rhor,rprimd,&
!!$&    scf_history,symrec,wffnew,wffnow,wvl,xred_next,xred,ylm,ylmgr)
    do iatom=1,natom
     do idim=1,3
      fcart_m(idim,iatom)=fcart(idim,iatom)/amass(iatom)
     end do
    end do
!   update of velocities and positions
    vel(:,:) = vel(:,:) + b8 * dtion * fcart_m(:,:)
    xred(:,:) = xred_next(:,:)
    xcart_next(:,:) = xcart_next(:,:) + a10 * dtion * vel(:,:)
!   Convert xcart_next to xred_next (reduced coordinates) for scfcv
    call abi_xredxcart(natom,-1,rprimd,xcart_next,xred_next)    

!step 11 of 15
!   Computation of the forces for the new positions
    iapp=-1
    if(itime>0)iapp=itime
!!$    call scfcv(acell,atindx,atindx1,cg,cpus,densymop_gs,dtefield,dtfil,dtset,ecore,&
!!$&    eigen,hdr,iapp,indsym,initialized,&
!!$&    irrzon,kg,mpi_enreg,&
!!$&    nattyp,nfftf,npwarr,nspinor,occ,&
!!$&    pawang,pawfgr,pawrad,pawrhoij,pawtab,phnons,psps,&
!!$&    pwind,pwind_alloc,pwnsfac,rec_set,resid,results_gs,rhog,rhor,rprimd,&
!!$&    scf_history,symrec,wffnew,wffnow,wvl,xred_next,xred,ylm,ylmgr)
    do iatom=1,natom
     do idim=1,3
      fcart_m(idim,iatom)=fcart(idim,iatom)/amass(iatom)
     end do
    end do
!   update of velocities and positions
    vel(:,:) = vel(:,:) + b9 * dtion * fcart_m(:,:)
    xred(:,:) = xred_next(:,:)
    xcart_next(:,:) = xcart_next(:,:) + a11 * dtion * vel(:,:)
!   Convert xcart_next to xred_next (reduced coordinates) for scfcv
    call abi_xredxcart(natom,-1,rprimd,xcart_next,xred_next)    

!step 12 of 15
!   Computation of the forces for the new positions
    iapp=-1
    if(itime>0)iapp=itime
!!$    call scfcv(acell,atindx,atindx1,cg,cpus,densymop_gs,dtefield,dtfil,dtset,ecore,&
!!$&    eigen,hdr,iapp,indsym,initialized,&
!!$&    irrzon,kg,mpi_enreg,&
!!$&    nattyp,nfftf,npwarr,nspinor,occ,&
!!$&    pawang,pawfgr,pawrad,pawrhoij,pawtab,phnons,psps,&
!!$&    pwind,pwind_alloc,pwnsfac,rec_set,resid,results_gs,rhog,rhor,rprimd,&
!!$&    scf_history,symrec,wffnew,wffnow,wvl,xred_next,xred,ylm,ylmgr)
    do iatom=1,natom
     do idim=1,3
      fcart_m(idim,iatom)=fcart(idim,iatom)/amass(iatom)
     end do
    end do
!   update of velocities and positions
    vel(:,:) = vel(:,:) + b10 * dtion * fcart_m(:,:)
    xred(:,:) = xred_next(:,:)
    xcart_next(:,:) = xcart_next(:,:) + a12 * dtion * vel(:,:)
!   Convert xcart_next to xred_next (reduced coordinates) for scfcv
    call abi_xredxcart(natom,-1,rprimd,xcart_next,xred_next)    

!step 13 of 15
!   Computation of the forces for the new positions
    iapp=-1
    if(itime>0)iapp=itime
!!$    call scfcv(acell,atindx,atindx1,cg,cpus,densymop_gs,dtefield,dtfil,dtset,ecore,&
!!$&    eigen,hdr,iapp,indsym,initialized,&
!!$&    irrzon,kg,mpi_enreg,&
!!$&    nattyp,nfftf,npwarr,nspinor,occ,&
!!$&    pawang,pawfgr,pawrad,pawrhoij,pawtab,phnons,psps,&
!!$&    pwind,pwind_alloc,pwnsfac,rec_set,resid,results_gs,rhog,rhor,rprimd,&
!!$&    scf_history,symrec,wffnew,wffnow,wvl,xred_next,xred,ylm,ylmgr)
    do iatom=1,natom
     do idim=1,3
      fcart_m(idim,iatom)=fcart(idim,iatom)/amass(iatom)
     end do
    end do
!   update of velocities and positions
    vel(:,:) = vel(:,:) + b11 * dtion * fcart_m(:,:)
    xred(:,:) = xred_next(:,:)
    xcart_next(:,:) = xcart_next(:,:) + a13 * dtion * vel(:,:)
!   Convert xcart_next to xred_next (reduced coordinates) for scfcv
    call abi_xredxcart(natom,-1,rprimd,xcart_next,xred_next)    

!step 14 of 15
!   Computation of the forces for the new positions
    iapp=-1
    if(itime>0)iapp=itime
!!$    call scfcv(acell,atindx,atindx1,cg,cpus,densymop_gs,dtefield,dtfil,dtset,ecore,&
!!$&    eigen,hdr,iapp,indsym,initialized,&
!!$&    irrzon,kg,mpi_enreg,&
!!$&    nattyp,nfftf,npwarr,nspinor,occ,&
!!$&    pawang,pawfgr,pawrad,pawrhoij,pawtab,phnons,psps,&
!!$&    pwind,pwind_alloc,pwnsfac,rec_set,resid,results_gs,rhog,rhor,rprimd,&
!!$&    scf_history,symrec,wffnew,wffnow,wvl,xred_next,xred,ylm,ylmgr)
    do iatom=1,natom
     do idim=1,3
      fcart_m(idim,iatom)=fcart(idim,iatom)/amass(iatom)
     end do
    end do
!   update of velocities and positions
    vel(:,:) = vel(:,:) + b12 * dtion * fcart_m(:,:)
    xred(:,:) = xred_next(:,:)
    xcart_next(:,:) = xcart_next(:,:) + a14 * dtion * vel(:,:)
!   Convert xcart_next to xred_next (reduced coordinates) for scfcv
    call abi_xredxcart(natom,-1,rprimd,xcart_next,xred_next)    

!step 15 of 15
!   Computation of the forces for the new positions
    iapp=-1
    if(itime>0)iapp=itime
!!$    call scfcv(acell,atindx,atindx1,cg,cpus,densymop_gs,dtefield,dtfil,dtset,ecore,&
!!$&    eigen,hdr,iapp,indsym,initialized,&
!!$&    irrzon,kg,mpi_enreg,&
!!$&    nattyp,nfftf,npwarr,nspinor,occ,&
!!$&    pawang,pawfgr,pawrad,pawrhoij,pawtab,phnons,psps,&
!!$&    pwind,pwind_alloc,pwnsfac,rec_set,resid,results_gs,rhog,rhor,rprimd,&
!!$&    scf_history,symrec,wffnew,wffnow,wvl,xred_next,xred,ylm,ylmgr)
    do iatom=1,natom
     do idim=1,3
      fcart_m(idim,iatom)=fcart(idim,iatom)/amass(iatom)
     end do
    end do
!   update of velocities and positions
    vel(:,:) = vel(:,:) + b13 * dtion * fcart_m(:,:)
    xred(:,:) = xred_next(:,:)
    xcart_next(:,:) = xcart_next(:,:) + a15 * dtion * vel(:,:)
!   Convert xcart_next to xred_next (reduced coordinates) for scfcv
    call abi_xredxcart(natom,-1,rprimd,xcart_next,xred_next)    
! end of step 15

! step 16 : add a new call to scfcv to get the total energy at the final position
    iapp=-1
    if(itime>0)iapp=itime
!!$    call scfcv(acell,atindx,atindx1,cg,cpus,densymop_gs,dtefield,dtfil,dtset,ecore,&
!!$&    eigen,hdr,iapp,indsym,initialized,&
!!$&    irrzon,kg,mpi_enreg,&
!!$&    nattyp,nfftf,npwarr,nspinor,occ,&
!!$&    pawang,pawfgr,pawrad,pawrhoij,pawtab,phnons,psps,&
!!$&    pwind,pwind_alloc,pwnsfac,rec_set,resid,results_gs,rhog,rhor,rprimd,&
!!$&    scf_history,symrec,wffnew,wffnow,wvl,xred_next,xred,ylm,ylmgr)
    do iatom=1,natom
     do idim=1,3
      fcart_m(idim,iatom)=fcart(idim,iatom)/amass(iatom)
     end do
    end do
