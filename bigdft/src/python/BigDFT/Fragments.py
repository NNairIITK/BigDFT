AU_to_A=0.52917721092

class XYZfile():
    def __init__(self,filename=None,units='atomic'):
        self.filename=filename
        self.lines=[]
        self.units=units
        self.fac=1.0
        if units == 'angstroem': self.fac=AU_to_A
    def append(self,array,basename='',names=None):
        "Add lines to the file position list"
        nm=basename
        for i,r in enumerate(array):
            if names is not None: nm=basename+names[i]
            line=str(nm)
            for t in r:
                line+=' '+str(self.fac*t)
            self.lines.append(line+'\n')
    def dump(self,position='w'):
        "Dump the file as it is now ready for writing"
        import sys
        f=sys.stdout
        if self.filename is not None: f=open(self.filename,position)
        f.write(str(len(self.lines))+' '+str(self.units)+'\n')
        f.write('# xyz dump \n')
        #then the positions
        for l in self.lines:
            f.write(l)
        if self.filename is not None: f.close()

def open_xyz(filename,nat,unit,comment,position='a'):
    import sys
    f=sys.stdout
    if filename is not None: f=open(filename,position)
    if (position != 'a'):
        f.write(str(nat)+' '+str(unit)+'\n')
        f.write(comment+'\n')
    return f

def close_xyz(f,filename):
    if filename is not None: f.close()

def dump_xyz_positions(f,array,basename='',names=None):
    nm=basename
    for i,r in enumerate(array):
        if names is not None: nm=basename+names[i]
        f.write(str(nm)+' '+str(r[0])+' '+str(r[1])+' '+str(r[2])+'\n')

def dump_xyz(array,basename='',units='atomic',names=None,filename=None,position='a'):
    f=open_xyz(filename,len(array),units,'# xyz dump with basename "'+basename+'"',position)
    dump_xyz_positions(f,array,basename=basename,names=names)
    close_xyz(f,filename)

class Lattice():
    "Defines the fundamental objects to deal with periodic systems"
    def __init__(self,vectors):
        self.vectors=vectors
    def grid(self,origin=[0.0,0.0,0.0],extremes=None,radius=None):
        "produces a set of translation vectors from a given origin"
        import numpy as np
        transl=[]
        g=[[],[],[]] # the grid of discrete translations
        if extremes is not None:
            #print extremes
            for i,d in enumerate(extremes):
                for k in range(d[0],d[1]+1):
                   g[i].append(k)
            #print g
            for i in g[0]:
                arri=np.array(self.vectors[0])*i
                for j in g[1]:
                    arrj=np.array(self.vectors[1])*j+arri
                    for k in g[2]:
                        arrk=np.array(self.vectors[2])*k+arrj
                        vect=np.array(origin)+arrk
                        app=True
                        if radius is not None: app=np.linalg.norm(arrk) < radius
                        if app: transl.append(vect)
        return transl

class Fragment():
    protected_keys=['q0','q1','q2','sigma']
    def __init__(self,atomlist=None,id='Unknown',units='AU'):
        self.atoms=[]
        self.id=id
        self.to_AU=1.0
        if units == 'A': self.to_AU=1.0/AU_to_A
        if atomlist is not None:
            for i in atomlist:
                self.append(i)
            self.positions=self.__positions()  #update positions
    def __len__(self):
        return len(self.atoms)
    #def __str__(self):
    #    import yaml
    #    return yaml.dump({'Positions': self.atoms,'Properties': {'name': self.id}})
    def xyz(self,filename=None,units='atomic'):
        "Write the fragment positions in a xyz file"
        import numpy as np
        f=XYZfile(filename,units=units)
        names=[ self.element(at) for at in self.atoms]
        posarr=[ np.ravel(r) for r in self.positions]
        f.append(posarr,names=names)
        f.dump()
    def dict(self):
        "Transform the fragment information into a dictionary ready to be put as external potential"
        lat=[]
        for at in self.atoms:
            dat={}
            dat['r']=list(at[self.__torxyz(at)])
            dat['sym']=self.element(at)
            for k in self.protected_keys:
                if k in at: dat[k]=at[k].tolist()
            lat.append(dat)
        return lat
    def append(self,atom=None,sym=None,positions=None):
        if atom is not None:
            self.atoms.append(atom)
        elif sym is not None:
            self.atoms.append({sym: positions})
        self.positions=self.__positions()  #update positions
    def element(self,atom):
        "Provides the name of the element"
        el=self.__torxyz(atom)
        if el == 'r': el=atom['sym']
        return el
    def __torxyz(self,atom):
        "provide the key which contains the positions"
        ks=atom.keys()
        for k in ks:
            if k not in self.protected_keys and type(atom[k])==type([]):
                if len(atom[k])==3: return k
        print 'atom',atom
        raise ValueError
    def rxyz(self,atom):
        import numpy as np
        k=self.__torxyz(atom)
        return self.to_AU*np.array(atom[k])
    def __positions(self):
        import numpy
        return numpy.mat([ self.rxyz(at) for at in self.atoms ])
    def centroid(self):
        import numpy
        return numpy.ravel(numpy.mean(self.positions, axis=0))
    def transform(self,R=None,t=None):
        "Apply a rototranslation of the fragment positions"
        import wahba as w,numpy as np
        if t is None:
            self.positions=w.apply_R(R,self.positions)
        elif R is None:
            self.positions=w.apply_t(t,self.positions)
        else:
            self.positions=w.apply_Rt(R,t,self.positions)
        #then replace the correct positions at the atoms
        for at,r in zip(self.atoms,self.positions):
            k=self.__torxyz(at)
            at[k]=np.ravel(r).tolist()
        #further treatments have to be added for the atomic multipoles
    def q0(self,atom):
        "Provides the charge of the atom"
        charge=atom.get('q0')
        if charge is not None: charge=charge[0]
        return charge
    def q1(self,atom):
        "Provides the dipole of the atom"
        import numpy as np
        dipole=atom.get('q1') #they are (so far) always given in AU
        if dipole is not None: dipole=np.array([dipole[2],dipole[0],dipole[1]])
        return dipole
    def Q(self,atom=None):
        "Fragment Monopole, when the information is available. Select the element if possible"
        charge=0
        found=False
        for at in self.atoms:
            q0=self.q0(at)
            if q0 is not None and (atom is None or self.element(at)==atom):
                found=True
                charge+=q0
        if found:
            return charge
        else:
            return None
    def d0(self,center=None):
        "Fragment dipole, calculated only from the atomic charges"
        #one might added a treatment for non-neutral fragments
        #but if the center of charge is used the d0 value is zero
        import numpy as np
        if center is not None:
            cxyz=center
        else:
            cxyz=np.ravel(self.centroid())
        d0=np.zeros(3)
        found=False
        for at in self.atoms:
            if self.q0(at) is not None:
                found=True
                d0+=at['q0'][0]*(self.rxyz(at)-cxyz)
        if found:
            return d0
        else:
            return None
    def d1(self,center=None):
        "Fragment dipole including the atomic dipoles"
        import numpy as np
        d1=np.zeros(3)
        dtot=self.d0(center)
        if dtot is None: return dtot
        found=False
        for at in self.atoms:
            q1=self.q1(at)
            if q1 is not None:
                found=True
                d1+=q1
        if found:
            return d1+dtot
        else:
            return None

                        
class System():
    "A system is defined by a collection of Fragments. It might be given by one single fragment"
    def __init__(self,mp_dict=None,xyz=None,nat_reference=None,units='AU',transformations=None,reference_fragments=None):
        self.fragments=[]
        self.CMs=[]
        self.units=units
        if xyz is not None: self.fill_from_xyz(xyz,nat_reference)
        if mp_dict is not None: self.fill_from_mp_dict(mp_dict,nat_reference)
        if transformations is not None: self.recompose(transformations,reference_fragments)
    def __len__(self):
        return sum([len(frag) for frag in self.fragments])
    def fill_from_xyz(self,file,nat_reference):
        "Import the fragment information from a xyz file"
        fil=open(file,'r')
        nat=0
        iat=0
        frag=None
        for l in fil:
            try:
                pos=l.split()
                if len(pos) <=2: #these are the number of atoms
                    nt=int(pos[0])
                    nat-=nt
                    if len(pos)==2:
                        unt=pos[1]
                        if unt=='angstroem':
                            self.units='A'
                        elif unt=='atomic' or unt=='bohr':
                            self.units='AU'
                    if frag is not None: self.append(frag)
                    frag=Fragment(units=self.units)
                    iat=0
                elif len(pos)>0:
                    if nat_reference is not None and iat == nat_reference: #we should break the fragment, alternative strategy
                        if frag is not None: self.append(frag)
                        frag=Fragment(units=self.units)
                        iat=0  
                    frag.append({pos[0]: map(float,pos[1:])})
                    nat+=1
                    iat+=1
            except Exception,e:
                print 'Warning, line not parsed: "',l,e,'"'
        if iat != 0: self.append(frag) #append the remaining fragment
    def fill_from_mp_dict(self,mpd,nat_reference=None):
        "Fill the System from a dictionary of multipole coefficients"
        frag=Fragment(units=self.units)
        iat=0
        for at in mpd:
            #frag.append(sym=at['sym'],positions=at['r'])
            frag.append(at)
            iat+=1
            if nat_reference is not None and iat == nat_reference: 
                if len(frag) !=0: self.append(frag)
                frag=Fragment(units=self.units)
                iat=0
    def xyz(self,filename=None,units='atomic'):
        import numpy as np
        f=XYZfile(filename,units)
        for frag in self.fragments:
            names=[ frag.element(at) for at in frag.atoms]
            posarr=[ np.ravel(r) for r in frag.positions]
            f.append(posarr,names=names)
        f.dump()
    def dict(self,filename=None):
        atoms=[]
        for f in self.fragments:
            atoms+=f.dict()
        if self.units != 'A': 
            print 'Dictionary version not available if the system is given in AU'
            raise Exception
        dc={'units': 'angstroem','global monopole': float(self.Q()), 'values': atoms}
        return dc
    def append(self,frag):
        "Append a fragment to the System class"
        assert isinstance(frag,Fragment)
        self.fragments.append(frag)
        self.CMs.append(frag.centroid())  #update center of mass
    def pop(self,ifrag):
        "Pop the fragment ifrag from the list of fragments"
        self.CMs.pop(ifrag)
        return self.fragments.pop(ifrag)
    def centroid(self):
        "Center of mass of the system"
        import numpy
        return numpy.mean(self.CMs, axis=0)
    def central_fragment(self):
        "Returns the fragment whose center of mass is closest to the centroid"
        import numpy as np
        return np.argmin([ np.dot(dd,dd.T) for dd in (self.CMs - self.centroid())])
        #return self.fragments[imin]
    def fragment_transformation(self,frag1,frag2):
        "returns the transformation among fragments if exists"
        try:
            import wahba
            roto,translation,J=wahba.rigid_transform_3D(frag1.positions,frag2.positions)
        except:
            roto,translation,J=(None,None,1.0e10)
        return roto,translation,J
    def decompose(self,reference_fragments):
        "Decompose the system into reference fragments"
        assert type(reference_fragments) == type([])
        self.decomposition=[]
        for frag in self.fragments:
            transf=[]
            Js=[]
            for ref in reference_fragments:
                r,t,j=self.fragment_transformation(ref,frag)
                transf.append({'R':r,'t':t})
                Js.append(j)
            #choose the minimal one
            import numpy
            Jchosen=numpy.argmin(Js)
            ref=transf[Jchosen]
            ref['ref']=reference_fragments[Jchosen]
            ref['J']=Js[Jchosen]
            ref['id']=Jchosen
            self.decomposition.append(ref)
    def recompose(self,transformations=None,reference_fragments=None):
        "Rebuild the system from a set of transformations"
        import copy,numpy as np
        if transformations is not None: 
            RT=transformations
            self.decomposition=RT
        else:
            RT=self.decomposition
        self.fragments=[]
        self.CMs=[]
        for item in RT:
            if reference_fragments:
                idf=item['id']
                frag=copy.deepcopy(reference_fragments[idf])
            else:
                frag=copy.deepcopy(item['ref'])
            frag.transform(item['R'],item['t'])
            self.append(frag)
    def Q(self):
            "Provides the global monopole of the system given as a sum of the monopoles of the atoms"
            return sum([ f.Q() for f in self.fragments])

def frag_average(ref,flist,clean_monopole=True):
    "Perform the average in attributes of the fragments provided by the list, with the position of the first fragment"
    import copy,numpy as np
    keys=['q0','q1','q2']
    navg=len(flist)
    if navg==0: return ref
    favg=copy.deepcopy(ref)
    qtot=0.0
    for i,at in enumerate(favg.atoms):
        #form a fragment which has the positions of the references and 
        #neutral total monopole if asked for.
        for k in keys:
            population=[ f.atoms[i][k] for f in flist ] 
            vals=np.mean(population,axis=0)
            st=np.std(population,axis=0)
            at[k]=vals
            at['sigma'+k]=st
            #print 'test',k,np.max(abs(at['sigma'+k]/at[k]))
        qtot+=at['q0'][0]
    qtot/=float(len(favg.atoms))
    for at in favg.atoms:
        at['q0'][0]-=qtot
        #print 'retest',i,at
    return favg
                        
def distance(i,j):
    "Distance between fragments, defined as distance between center of mass"
    import numpy
    vec=i.centroid()-j.centroid()
    return numpy.sqrt(numpy.dot(vec,vec.T))

def wahba_fragment(frag1,frag2):
    "Solve the wahba's problem among fragments"
    import wahba #should be cleaned
    #For each of the fragment build the list of the coordinated
    roto,translation,J=wahba.rigid_transform_3D(frag1.positions,frag2.positions)
    return roto,translation,J

def rotot_collection(ref_frag,lookup,fragments):
    W=[]
    for f in lookup:
        refF=lookup[ref_frag]
        roto,translation,J=wahba_fragment(fragments[f],fragments[refF])
        if (J > 1.e-12):
            print 'Error',f,J,refF
            #try with the second ref
            refF2=lookup[ref_frag+1]
            roto2,translation2,J2=wahba_fragment(fragments[f],fragments[refF2])
            #print 'Error Now',f,J2,refF2
            if (J2< J):
                roto=roto2
                translation=translation2
                J=J2
                refF=refF2
        W.append({'R':roto,'t':translation,'J':J,'ref':refF})
    return W

if __name__ == '__main__': 
    #extract fragments
    import sys,numpy
    one1=System(xyz='one-1.xyz')
    print 'Parsed',len(one1.fragments)
    print one1.xyz()
    two=System(xyz='two.xyz',nat_reference=len(one1),units='A')
    two.decompose(one1.fragments)
    trans=two.decomposition
    PC1=one1.fragments[0]
    print 'one',PC1.centroid()
    for frag,t in zip(two.fragments,trans):
        print 'ff',frag.centroid()
        print 't',t['t']
        print 'R',t['R']
    print two.CMs
    print trans
    two2=System(transformations=trans)
    two2.xyz('two-2.xyz',units='angstroem')
    #now rigidify the big case scenario
    #read lattice coordinates
    fil=open('lattice.txt','r')
    acell=[eval(l.strip('\r\n')) for l in fil]
    latt=Lattice(acell)
    print latt.vectors
    print 0.5*numpy.array(latt.vectors)
    print two.CMs[1]-two.CMs[0]
    #find the positions of the center of mass of the big system
    bigold=System(xyz='BigCase.xyz',nat_reference=36)
    icen=bigold.central_fragment()
    print icen,bigold.CMs[icen],bigold.CMs[icen-1]
    samples=[bigold.CMs[icen],bigold.CMs[icen-1]]
    extremes=[[-5,5],[-5,5],[-1,1]]
    grid=[]
    for oxyz in samples: #two.CMs:
      grid += latt.grid(extremes=extremes,origin=oxyz)
      #grid centroid
    cent=numpy.ravel(numpy.mean(numpy.mat(grid), axis=0))
    trans=[]
    limit=len(grid)/len(two.CMs)
    for i,d in enumerate(grid):
        ref=two2.fragments[i/limit]
        if numpy.linalg.norm(d-cent) < 25.0:
            trans.append({'t':numpy.mat(d).T/AU_to_A,'ref':ref,'R':None})

    big=System(transformations=trans)
    big.xyz('Bigcase-2.xyz',units='angstroem')
    cents=XYZfile('Bigcase2-centroids.xyz',units='angstroem')
    cents.append(big.CMs,basename='Cen')
    cents.dump()
    icen = big.central_fragment()
    print 'the central fragment is',icen
    #find the atoms
    iat=0
    for i,f in enumerate(big.fragments):
        if i==icen: print 'from',iat+1
        iat+=len(f)
        if i==icen:
            print 'to',iat+1
            break
    exit(0)
    filename=sys.argv[1]
    limit=36 #maximum value of each fragment
    fragments=[]
    #try to initialize with the class
    
    fil=open(filename,'r')
    count=0
    nat=0
    iat=0
    frag=None
    for l in fil:
        count+=1
        try:
          pos=l.split()
          if len(pos) ==1: #these are the number of atoms
            nt=int(pos[0])
            nat-=nt
            if frag is not None: fragments.append(frag)
            frag=Fragment()
            iat=0
          elif len(pos)>0:
            if iat == limit: #we should break the fragment, alternative strategy
                if frag is not None: fragments.append(frag)
                frag=Fragment()
                iat=0  
            frag.append({pos[0]: map(float,pos[1:])})
            nat+=1
            iat+=1
        except Exception,e:
            print 'error',l,e
            break
    
    print 'calculation finished',len(fragments),'balance',nat
    
    #find the F4TCNQ
    F4TCNQs=[]
    PCs=[]
    CMs=[] #ordered center of mass
    for f,frag in enumerate(fragments):
        CMs.append(frag.centroid())
        if len(frag) < 36:
            F4TCNQs.append(f)
        else:
            PCs.append(f)
    #find the central molecule
    centroid= numpy.mean(CMs, axis=0)
            
    print 'species identified:',len(F4TCNQs),'F4TCNQ and',len(PCs),' pentacenes, tot',len(fragments),len(CMs)
    
    #now append the fragments to the System class
    stm=System(xyz=filename,nat_reference=36)
    print len(stm.fragments),'before'
    for frag in fragments:
        stm.append(frag)
    print len(stm.fragments),'after'
    print centroid,[i for i in CMs[0]],[i for i in centroid],stm.centroid(),stm.central_fragment()
    
    refF4=0
    refPC=0
    #try:
    #    icen=F4TCNQs.index(imin)
    #    refF4=icen
    #except:
    #    icen=PCs.index(imin)
    #    refPC=icen
        
    #check if now all the atoms are the rototranslation of the same fragment and find the transformation
    W_F4=rotot_collection(refF4,F4TCNQs,fragments)
    W_PEN=rotot_collection(refPC,PCs,fragments)
    
        
    #print CMs
    #search the NN of each of the F4TCNQs
    DFP=[]
    DFF=[]
    for f in F4TCNQs:
        DFP.append(numpy.array([distance(f,p) for p in PCs]))
        DFF.append(numpy.array([distance(f,p) for p in F4TCNQs]))
    
        
    #Then we can classify the attributes of each pentacene accordingly to the limiting distance
    
    threshold=10.0
    
    for i,f in enumerate(F4TCNQs):
        import yaml
        if (i==0):
            print yaml.dump(fragments[f])
            print 'test wahba'
            roto,translation,J=wahba_fragment(fragments[f],fragments[f])
            print roto
            print translation
            print 'Rototranslation Error',J
        iPC=0
        for dist in DFP[i]:
            if dist < threshold:
                iPC+=1
        iFF=0
        for dist in DFF[i]:
            if dist < threshold and dist !=0:
                iFF+=1
        print i,iFF,iPC

    
