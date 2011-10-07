find . -iname "*.f90"  | grep -v unused | xargs ctags -e


-----------------------------
rxyz(ix, iat)
        ix da 1 a 3
        iat da 1 a nat ; nat contenuto in atoms%nat

	 Le posizioni sono shiftate da system_size in modo 
        che le sfere coarse centrate sugli atomi non 
        vadano al di sotto di zero in alcuna direzione Free.
         Nessuno shift viene applicato da system_size nelle direzioni
        periodiche.

-----------------------------------------------------
subroutine	explain_pspreading()
end

psppar
   the are three parametrisations: 

 GTH Goedecker Teter Hutter     PHYSICAL REVIEW B VOLUME 54, NUMBER 3 15 JULY 1996
   
 HGH Hartwigsen C, Goedecker S, Hutter J (1998)   PHYSICAL REVIEW B VOLUME 58, NUMBER 7

 HGH-K :KRACK Theor Chem Acc (2005) 114: 145

 The first row in the file it is always a comment that one skips
 From the second row one reads (ityp is the atom type)
	atoms%nzatom(ityp)  which is Z (real one )  of the atom and
        atoms%nelpsp(ityp)  which is the number of electron in the valence
        or in other terms the charge  Zion of the core
 from the third row one read
      	atoms%npspcode(ityp)  which is the numcode of the  pesudo : 2=GTH , 3=HGH, 10=HGH-K
        ixcpsp  is the numcode  of the functional  of density  with which the pseudo has been calculated .

 When  npspcode = 2 oppure 3 :
    atoms%psppar(0,j,ityp),j=0,4
       that' s  rloc e c1,c2,c3,c4 which appear in the formul for the local part of the potential
              Vloc(r)=-Zion/r erf(r/(sqrt2 rloc)) + exp(-(r/rloc)**2/2) *Sum_i( Ci *(r/rloc)**(2*(i-1)))
       where zion is  atoms%nelpsp(ityp) 

    ---npspcode=2 -----

     the non local part for  npspcode=2  is    diagonally written in the radial projectors basis.
     There are two projectors for the s waves and one for p 
       Thus  psppar(1,j,ityp) is for  s waves con j=0 for the tipical radius  r_l and j=1,2 for the  two coefficients  h.
    
	The form of the nonlocal part  est, per npscode==2

               Somma_i somma_lm Y_lm(r) Y*_lm(r')   h^l_i  p^l_i(r)  p^l_i(r')

	where the projectors  p  are, besides   factors which renormalise them to 1,

        p^l_i (r) = r^(l+2(i-1)) exp(-(r/r_l)**2 ) 

    ---npspcode=3 -----

	This case is very similar to npspcode=2 the difference is that the pesudo potentials are read 
        for the waves  s,p,d,f. Then one read cioe' psppar(l,j,ityp) with l going from 1 a 4 
        ( missing data from the psp file are continued with zeros )
        The other difference is that  h coeffs are read for  j=1,2,3 (  missing are zero filled).


	Another difference is that the h coeffs  h_i that one reads, become in reality the h_{ii} diagonal
	of a matrix that is filled with the formulas below.
	The last difference is that for  l>1  (l=1 is s )  one  skips a line containing  k coeffs 
	who are there for spin-orbit but are not yet treated  in Big-dft

	Out-of-diagonal terms are       

               L=0
                       (1,2) = -1/2 sqrt(3/5) (2,2)
                       (1,3) =  1/2  sqrt(5/21) (3,3)
                       (2,3) = -1/2 sqrt(100/63) (3,3)

               L = 1  
                       (1,2) = -1/2 sqrt(5/7) (2,2)
                       (1,3) =  1/6  sqrt(35/11) (3,3)
                       (2,3) = -1/6 14/sqrt(11) (3,3)

	       L = 2
                       (1,2) = -1/2 sqrt(7/9) (2,2)
                       (1,3) =  1/2  sqrt(63/143) (3,3)
                       (2,3) = -1/2 18/sqrt(143) (3,3)

	                       
    -- npspcode=10 ------------------

	very similar to the case above but the file has a variable number of terms.
	The out-of-diagonal part is not completed but read from the file 

	The row for the local potentail has a number nn indicating the number of Ci terms (i=1,nn)
	 contained in the row
        Before the eventual rows for s p d f non local parts the is a line with the number nlterms
	which is the number of channel ( rows for s p d f ) to be read 


       read(11,*) atoms%psppar(0,0,ityp),nn,(atoms%psppar(0,j,ityp),j=1,nn) !local PSP parameters
        read(11,*) nlterms !number of channels of the pseudo
        prjloop: do l=1,nlterms
           read(11,*) atoms%psppar(l,0,ityp),nprl,atoms%psppar(l,1,ityp),&  !  H11
                (atoms%psppar(l,j+2,ityp),j=2,nprl) !h_ij terms  j+2 parte da 4 :  eventually H12  eventually H13
           do i=2,nprl
              read(11,*) atoms%psppar(l,i,ityp),&  ! ev H22 and after ev H33
	               (atoms%psppar(l,i+j+1,ityp),j=i+1,nprl) ! h_ij  if nprl=3  this will be  H23
           end do
           if (l==1) cycle
           do i=1,nprl
              read(11,*) !k coefficients, not used
           end do
        end do prjloop

	ATTENTION : nprl non e'mai superiore a 3
	0        1         2         3           4        5         6          7
        rc       h11       h22       h23         h12      h13       h23

-----------------------------------------------
subroutine explain_nsccode
end

nsccode:
  based on powers of 4 in BigDft. When initialised from subroutine eleconf
  a conversion is applied cause in eleconf it is coded with powers of 10.
  Used to tag some orbitals n l( for each l the lowest in n are taken)
  This orbitals will be taken anyway by input_wf_diag regardeles of the energy
  given at the first diag step, cause this energy could bring to loss of some 
  important orbitals , like Zn 3d shell close to Ef in porfirine.

  nsccode   Semicore orbitals, indicated as an integer.
!!             The integer is the n_s + 4*n_p + 16* n_d + 64* n_f
!!             where n_l are the number of semicore orbitals for a given angular momentum
!!             starting from the lower level of course

  Questo in read_eleconf:

  !then calculate the nsccode
  nsccode=0
  do lsc=1,4
     do i=1,nlsc(lsc)
        nsccode=nsccode+4**(lsc-1)
     end do
  end do


nsccode  is reset to zero in routine correct_semicore
if one impose by in input a modification of the charge.

------------------------------------------
atoms%alat1 , atoms%alat2, atoms%alat3 
   le tre dimensioni della scatola.
   In realta in input si leggono 1,2,3,4,5,6
   cioe 1 (2) 3 (4,5) 6
   Dove 2,4,5 sono nulli nel caso ortorombico
   e 1->1 3->2 6->3


----------------------------------------------

subroutine explain_radii_cf
end

radii_cf(ntypes,3)

  The first two are coarse and fine radius for functions,
the third one is the coarse for the projectors.
  They are read from the psp file. If on this file
only first two are given, the proj-coarse is set equal to the fine radius for functions.
  The fine radius for projectors is always set equal to the fine for functions.

  These radius  are expanded by factors to define the radius of the 
containing spheres.
  To do this one aplies crmult, frmult for functions, while for projectos
one applies cpmult=fpmult=crmult ( da riverificare )

  If the radius are completely missing form pspfile then they are obtained
the coarse as the tipical radius of an evanesceting exponential wave having 
energy=ehomo  defined in eleconf.
  The fine radius is taken as the maximum ot te tipical radius of the nonlocal part


----------------------------------------------------

subroutine explain_nspinor
end

nspinor can be retouched later, for example by function orbitals_descriptors:
if it comes out that the calculation is complex ( not a gamma point ) 
one sets nspinor=2 for the two components : real and imaginary

  select case(nspin)
     case(1)
        nsp=1
        nspinor=1
        noncoll=1
     case(2)
        nsp=2
        nspinor=1
        noncoll=1
     case(4)
        nsp=1
        nspinor=4
        noncoll=2
     case default
        write(*,*)' ERROR: nspin not valid:',nspin
        stop
  end select

-------------------------------------------------------------------
subroutine explain_natpol
end

natpol is a code for  nchrg e nspol  and sign of the charge
  atoms%natpol(iat)=1000*nchrg+nsgn*100+nspol

natpol is set at the beginning, when reading atomic positions ( routine read_atomic_positions )
by routine parse_extra_info if there are extra info on the position line

-------------------------------------------------------------------
subroutine explain_aocc
end
   Scheme of building aocc:
          when atom is not in input.occup all this in atomic_occupation_numbers
             neleconf ---> eleconf( changed to satisfie ichg charging charge in correct_semicore) -->  aocc by at_occnums
          if iat  in input.occup  then read_eleconf  replace the above reading the input string


   at%aocc(nelecmax,atoms%nat+ndebug) :
	compressed form :
            for every  l=L+1 channel : 
                     a number equal to the number of  quantum numbers n :
                           for every n:  
                                nsp*noncoll*(2*L+1)  occupation   numbers
   all this on the linear array at%aocc with integers being converted to floats
                              

     n the inner nsp*noncoll*(2*L+1)  long block, 
          the fastest dimension runs over  ncoll ( when it is 2 for noncollinear)
      second fastest is over m
     slowest is over  spin ( if >1)

    here the extract from read_eleconf

     aocc(iocc)=real(nl(l),gp)
     do inl=1,nl(l)
        do ispin=1,nspin
           do m=1,2*l-1
              do icoll=1,noncoll !non-trivial only for nspinor=4
                 iocc=iocc+1
                 aocc(iocc)=allocc(icoll+(m-1)*noncoll+(ispin-1)*(2*l-1)*noncoll,inl,l)
              end do
           end do
        end do
     end do
  end do

------------------------------------------------------------------------------------------------
subroutine explain_orbital_data


- orbitals_data

   - orbs%nkpt number of  points
   - orbs%kpts(3,orbs%nkpts+ndebug)  real(gp)    : the k points
   - orbs%norb_par(0:npar-1)  for each process the number of orbitals trated by the process 
     Their total number is  norb_tot = norb*orbs%nkpts.

     The fractioning is made by  parallel_repartition_with_kpoints.
     One considers the fractioning for [0:nproc-1] in  nkpts segments.
     For each segment one considers if it contains full segments : processes that 
      will be entirely dedicated to a single  ikpt .
      On the wings  points are distributed pro-rata 
      !!  Attention : 
      when a segment is on the saddle of two processesors 
      one distributes nobj/2 on the firs and the rest to the other . 
      This distribution half/half might be non equilibrated 

    - orbs%iskpts +1  is the first ikpt treated paritally or entirely by iproc 
      !!! se l' indicizzazione parte da 1 !!!
   
   - orbs%nkptsp  is the number of  kpoints  touched by  iproc

   - orbs%isorb is the number of orbitals treated by all process having rank 
   less than iproc

   -  orbs%norb=norb is the number of orbital per  kpoint

   - orbs%norbp=orbs%norb_par(iproc)  is the number of orbitals treated by iproc 

   - orbs%iokpt is a vector having lenght  norbp and for each orbitals treated by iproc 
   tells the index of the corresponding kpoint

--------------------------------------------------------------------------
subroutine explain_grid_dimensions()
end

n1,n2,n3,n1i,n2i,n3i
   ~~~~~~~~~~`
   These variables are containe in the member d of the structur of type locreg_descriptors
   The memebr  d is of the type 
     type, public :: grid_dimensions
        integer :: n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,n1i,n2i,n3i
     end type grid_dimensions
           The nfu_s  e nfl_s  are the extrema  (upper lower)
           of the box containing fine regions.
   ~~~~~~~~~~~~~~~~~~~


   NOTA:  important difference : n  is the maximum index starting from 0
         for the coarse grid while ni is the number of points
         of the fine  grid

 n1 n2 , n3  are the maxium value of the grid starting from zero.

-F-  For this reason in the F case 
   alat=n*h
  The point n being at distance  alat  from point  0.

-P-  In the case  P instead
   alat=(n+1)*h 
 the point n being at distance h from the following point n+1 ( not comprised in the grid)
 which is equivalent to the point 0. The distance from 0 and n+1 being alat.
 This is made sure by  correct_grid routine.

n1i,n2i,n3i


  -F-
    
    If the coarse point are in the range 0:n
    the points of the fine grid will be in number

        2*n+31

    this because the filter h and g go from -7 to +8
    and the magic filters go from -7 to +8
   When doing wavelet synthesys ome starts from a range  0 : n
   and one gets a range goinf for  -7 to 2n+8.
   The  synthesys is:

	S(i) = Sum_j  H(i-2*j) *s_j +  G(i-2*j) *d_j

   where  S are the coeff on the  fine grid.

   Applying further the filter one obtains a range going from  -14 to 2n+16

   The Daubechies in  BigDft are maximally  symmetric. 

 - P -
    
   In this case on simply doubles the number of points.
   If the  range e'0:n  the number of coarse points is , 
    = n+1 
   therefore ni is
    ni=2n+2




-----------------------------------------
hybrid_on
      ~~~~~~~~~~~~~~~~~~~
       memorizzato in 
           Glr%hybrid_on
       dove Glr e di tipo locreg_descriptor
      ~~~~~~~~~~~~~~~~~~~~

 Se la scatola contenente le fini, estesa di un margine =14
che corriponde all'espanzione per sintesi, e'contenuta 
nella griglia 0:n1, 0:n2, 0:n2 allora si attiva l'ibrido.

!!!!!!!!!Da Controllare !!!!!!!!1
Questo vuole dire che nel caso periodico la parte esterna alla zona fine
la si manterra  in coarse.
Nel caso Free questo avviene sempre perche gli alat si adattano al raggio coarse.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


--------------------------------------------------------------------------------------
subroutine explain_kinetic_bounds()
end

these kinetic-bounds are active in F case or in the case hybrid_on ( make_bounds_per ).
They are used to calculate the kinetic operator only for the elements
of the box which will be conserved in the compressed form
In the calculation one reads element from the comprexed form

The loop on the result will go ( in x ) from  ibyz_c(1,i2,i3) to ibyz_c(2,i2,i3)
for i2 e i3 which run  on y and z of the box when one derives along x 
and so on  for the  other  directions of  derivation.

Same thing for fine quantities.


Glr%bounds%kb%ibyz_c, Glr%bounds%kb%ibxz_c, Glr%bounds%kb%ibxz_c


Glr%bounds%kb%ibyz_c(2,0:n2,0:n3+ndebug)
Glr%bounds%kb%ibxz_c(2,0:n1,0:n3+ndebug)
Glr%bounds%kb%ibxy_c(2,0:n1,0:n2+ndebug)

Glr%bounds%kb%ibyz_f(2,0:n2,0:n3+ndebug)
Glr%bounds%kb%ibxz_f(2,0:n1,0:n3+ndebug)
Glr%bounds%kb%ibxy_f(2,0:n1,0:n2+ndebug)


These quantities are calculated on the base of a file logrid
using spheres ( coarse and fine ).

This would complicated in the non-ortorombic case but
would be not-necessary because the non-ortorombic case is always
compact-periodic


--------------------------------------------------------------------

subroutine explain_segments()
end

mseg,mvctr

are the number of segments for the compression along x ( first index )
and the total number of points contained in all the segments.
The points are all the points contained in all coarse sphere for _c postfix
and the fine spheres for _f postfix.

Glr%wfd%nseg_c,Glr%wfd%nvctr_c
Glr%wfd%nseg_f,Glr%wfd%nvctr_f

num_segkeys(n1,n2,n3,nl1,nu1,nl2,nu2,nl3,nu3,logrid,mseg,mvctr)


The informations about segments will be contained in  keyg ( g as griglia )
and keyv ( v as vettore ).
keyv has dimension mseg and gives the  position in the --compressed--  vector
for the segment beginning.

keyg has dimensions (2,mseg) and gives for every segment, the beginning  (1) and end (2) 
of the  segment in the volume taken as array  1D ( total size (n1+1)*(n2+1)*(n3+1) )
( the x,y,z coordinates are obtainable from the sequential numbers) 

The information for the fine  quantities are storead in the same vector. They begin 
at  iseg=Glr%wfd%nseg_c+1


------------------------------------------------------------------------------
subroutine  explain_grow_shrink()
end

          Glr%bounds%kb%ibyz_c,Glr%bounds%gb%ibzxx_c,Glr%bounds%gb%ibxxyy_c,&
          Glr%bounds%kb%ibyz_f,Glr%bounds%gb%ibyz_ff,Glr%bounds%gb%ibzxx_f,Glr%bounds%gb%ibxxyy_f,&

These bounds are used for the magic filters expansion ( grow)
and the other way round ( shrink)  .. ( bra V ket ) 
in the case F and in the case P-hybridon.

In these case one can save on the loop extensions because the areas
dont cover completely the cube.


In the growin case one loops with magic filters on the fastest variable
using first  kb%ibyz_c limits taken from kinetic bounds.
The result are store permutating ciclically x,y,z. So that 
at the following steps ( y, z ) one always loops on the fastest variable 
One used a filter already convoluted with synthesys

For  shrinking 

  lr%bounds%kb%ibxy_c,lr%bounds%sb%ibzzx_c,lr%bounds%sb%ibyyzz_c,&
  lr%bounds%sb%ibxy_ff,lr%bounds%sb%ibzzx_f,lr%bounds%sb%ibyyzz_f,&

the loop is organized in the same way but one proceeds in the other directions
 fine grid -> coarse-grid-scaling-functions+wavelet.
Therefore first   yyzz_f poi zzx_f e xy_ff. Permutando alla stessa maniera


-----------------------------------------------------------------------------------------------------------------

nlpspd : non local psp descriptor

	nlpspd%nboxp_c(low=1 up=2, 1=x, y=2, z=3 , iat) : box contenente projector coarse

        --------------------------------------
        nlpspd%nseg_p(2*iat-1)
        nlpspd%nvctr_p(2*iat-1)
        nlpspd%nseg_p(2*iat)
        nlpspd%nvctr_p(2*iat)          
		questi array sono comulativi : 
	              (2*iat-1) numero dei segmenti del coarse di iat + tutti i precedenti fini e grossi 
	              (2*iat) numero dei segmenti del fini di iat + tutti i precedenti grossi  e fini 


	nlpspd%nproj  :  numero di proiettori ( per un solo k ) 
	nlpspd%nprojel 
	Se DistProjApply est true
		e' il massimo dei 
		( nlpspd%nvctr_p(2*iat-1)-nlpspd%nvctr_p(2*iat-2)+(nlpspd%nvctr_p(2*iat)-nlpspd%nvctr_p(2*iat-1))*7)*mproj
		( mproj e' il numero totale di projettori per l' atomo in questione. Cioe' per tutti i canali
	          e per tuttii proiettori del canale si somma 2l+1 ) 
	se e' false
	        non e' il massimo ma la somma di tutti quanti

	C'e ' infine un ulteriore fattore che lo moltiplica : questo e' 1 o 2 nel caso DistProjApply ( 2 se c' e' un k complesso )
        se no e' la somma, per tutti i k toccati dal processore (nkptsp a partire da iskpts+1), di 1 o 2 ( se c' e' il complesso ) 


	Per i proiettori si keyg e keyv sono messi tutti negli stessi array

	nlpspd%keyg_p(1,iseg),nlpspd%keyv_p(iseg)
	
	Dove quelli per la coarse iniziano in    iseg=nlpspd%nseg_p(2*iat-2)+1 e sono in numero 
        mseg=nlpspd%nseg_p(2*iat-1)-nlpspd%nseg_p(2*iat-2)
	Similmente per i fini.


	L'ordinamento degli m, quando si somma su uno spazio angolare dei proiettori, e'

        l=p     x,y,z
        l=d     yz, xz, zy, x2-y2, 2z2-x2-y2
        l=f     x( -4z2 +x2+y2 )    ,  y( -4z2 +x2+y2 ),     z(-2z2+ 3x2+3y2    ),
                x( x2-3y2) , y(-3x2+y2  ) ,    z(x2-y2) , xyz



	Dentro all'array proj ( che e'tenuto  ad alto livello i proiettori sono ( loop esterno
           verso interno ) messi dentro per k (quelli toccati dal processore ),
             per atomo, per l , per indice i (massimo 3),
              eventualmente per parte reale/ immaginaria se e'il caso,
            infine per compresso coarse e compresso fine ( wavelet 7 a 7)

----------------------------------------------------------------------------

  character(len=*), parameter :: subname='orbitals_communicators'
  logical :: yesorb,yescomp
  integer :: jproc,nvctr_tot,ikpts,iorbp,jorb,norb_tot,ikpt,i_stat,i_all
  integer :: ncomp_res,nkptsp,ierr
  integer, dimension(:), allocatable :: mykpts
  logical, dimension(:), allocatable :: GPU_for_comp
  integer, dimension(:,:), allocatable ::


  -- communication_arrays

   -  nvctr_par che viene poi messo in comm%nvctr_par di comm ( communication_array)
       contiene dapprima in 
                 nvctr(jproc,0)
       il numero di gradi di liberta ospitati dal processo jproc su 
       un totale di nkpts*(nvctr_c+7*nvctr_f) sequenzialmente, per la rappresentazione
       trasposta.
       Questo frazionemento e' stabilito da parallel_repartition_with_k_points
       
       Dopodiche il numero nvctr(jproc,0) viene distribuito(e lasciato locale ) su 
       comms%nvctr_par(jproc,ikpts) dove ikpts parte da 1 e arriva a orbs%nkptsp
       Questo da il numero di componenti trattate da jproc che appartengono 
       al punto k :  orbs%iskpts+ikpts.

       Questo viene fatto all interno della routine  orbitals_communicators.

       ATTENZione : questa routine ricalcola alcune quantita come 
       orbs%iskpts che gia erano calcolate da orbitals_descriptors.
       Cio' va bene perche tali quantita coincidono nella distribuzione
       degli orbitali con quelle nella distribuzione delle componenti.
        
       --  orbs%ikptproc(1:nkpts)
         !this function which associates a given k-point to a processor in the component distribution
  	 !the association is chosen such that each k-point is associated to only
  	 !one processor
  	 !if two processors treat the same k-point the processor which highest rank is chosen
                 


  -- comms%ncntd(0:nproc-1+ndebug)
    la lunghezza del receive buffer in provenienza da jproc per passare dalla 
    distribuzione delle componenti ( t come trasposta )
    alla distribuzione degli orbitali ( d come diretta)
    Se iproc deve tenere norb_par(iproc,ikpts) per il punto ikpts
    nella distribuzione diretta , ogni processo jproc che ha nella distribuzione
    trasposta  nvctr_par(jproc,ikpts)  componenti, le deve mandare
    per tutti gli orbitali norb_par(iproc,ikpts) visto che li ha tutti.
    Il tutto moltiplicato per orbs%nspinor

  -- comms%ndspld(jproc)
    e' il displacement per l' all gather, comms%ndspld(jproc)=comms%ndspld(jproc-1)+comms%ncntd(jproc-1)
    Il buffer dovra poi essere smistato

  -- comms%ncntt(jproc)  

    lunghezza del receive buffer per passare alla forma trasposta.
    Ogni jproc  contiene norb_par(jproc,ikpts) orbitali per un ikpts
    e per ognuno di questi deve trasmettere  nvctr_par(iproc,ikpts)
    componenti ( le ha tutte ).

  -- comms%ndsplt(jproc)
    displacement per l' all gather
       comms%ndsplt(jproc)=comms%ndsplt(jproc-1)+comms%ncntt(jproc-1)




-----------------------------------------------------------------------------------

!!    n3d         third dimension of the density. For distributed data, it takes into account 
!!                the enlarging needed for calculating the XC functionals.
!!                For global data it is simply equal to n03. 
!!                When there are too many processes and there is no room for the density n3d=0
!!    n3p         third dimension for the potential. The same as n3d, but without 
!!                taking into account the enlargment for the XC part. For non-GGA XC, n3p=n3d.
!!    n3pi        Dimension of the pot_ion array, always with distributed data. 
!!                For distributed data n3pi=n3p
!!    i3xcsh      Shift of the density that must be performed to enter in the 
!!                non-overlapping region. Useful for recovering the values of the potential
!!                when using GGA XC functionals. If the density starts from rhopot(1,1,1),
!!                the potential starts from rhopot(1,1,i3xcsh+1). 
!!                For non-GGA XCs and for global distribution data i3xcsh=0
!!    i3s         Starting point of the density effectively treated by each processor 
!!                in the third direction.
!!                It takes into account also the XC enlarging. The array rhopot will correspond
!!                To the planes of third coordinate from i3s to i3s+n3d-1. 
!!                The potential to the planes from i3s+i3xcsh to i3s+i3xcsh+n3p-1
!!                The array pot_ion to the planes from i3s+i3xcsh to i3s+i3xcsh+n3pi-1
!!                For global disposition i3s is equal to distributed case with i3xcsh=0.
!!


  nscatterarr(iproc,1)=n3d
  nscatterarr(iproc,2)=n3p
  nscatterarr(iproc,3)=i3s+i3xcsh-1
  nscatterarr(iproc,4)=i3xcsh

  Questi sopra calcolati per ogni jproc pure.

  ngatherarr(:,1)=n1i*n2i*nscatterarr(:,2)
  ngatherarr(:,2)=n1i*n2i*nscatterarr(:,3)

-------------------------------------------------------------------------------

domanda : 
        
         perche ci vuole psoffset nel caso periodico? 

------------------------------------------------------------------

subroutine	explain_paw_psp_terms_in_atom_data()
end

      -- for each  itype, paw_NofL(itype) is the number of different L having paw projectors
      These paw_NofL are stored sequentially in 
      -- paw_l(:),  in ascending order of itype and      ascending order of l
      --paw_nofchannels : The array runs parallel to paw_l, its entries are, for each l, the number of channels
      -- The array paw_nofgaussians contains, for each l, the number of gaussian
               concurring to it, they are put in sequential order : loopover_itype, over_l.
      -- paw_Greal, in sequeantial order  for each (loopover_itype, over_l) ,
        is the real part of the coefficient of r**2 in the exponent it is constant for all gaussians
        there is one element in paw_Greal for each l   
      -- paw_Gimag is the imaginary part ( cos/sin). There are paw_nofgaussians elements for each
          (loopover_itype, over_l). The first exponent coefficient for cos is 0 ( instead of startign from
      the first   paw_Gimag 
      __ paw_Gcoeffs : there are 2*paw_nofgaussians elements for each (loopover_itype, over_l, over_channel) 
         the 2* is there because it comprises the loop over cos/sin  coefficients : cos sin cos sin...

      __ paw_H_matrices, paw_S_matrices, in sequential order, for each l, the matrix connecting the channels


           explain_paw_psp_terms_in_PAWproj_data_type



subroutine explain_paw_psp_terms_in_PAWprojdatatype()
end
     
     !! all data are in sequential order, projector after projector
     !! Each iproj belongs to a block of iproj_to_paw_nchannels(iproj)*(2*l-1), 
     !! corresponding to a give l = iproj_to_l(l) for a given atom
     !! The projectors of a given block are linked by matrices contained in 
     !! paw_matrices which gives the patch
     !! The corrections to overlap is contained in paw_matrices_S.
     !!  The loop over m is the innest one.

end
