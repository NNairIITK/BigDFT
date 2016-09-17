import numpy

AU_eV=27.31138386

def get_ev(ev,keys=None,ikpt=1):
    "Get the correct list of the energies for this eigenvalue"
    res=False
    if keys is None:
        ener=ev.get('e')
        spin=ev.get('s')
        kpt=ev.get('k')
        if not kpt and ikpt==1: 
            kpt=True
        elif kpt and kpt != ikpt: 
            kpt=False
        if ener and (spin==1 or not spin):
            if kpt: res=[ener]
        elif ener and spin==-1:
            if kpt: res=[None,ener]
    else:
        for k in keys:
            if k in ev:
                res=ev[k]
                if type(res) != type([]): res=[res]
                break
    return res

class BandArray(numpy.ndarray):
    "Defines the array of data for one band. It is a dictionary which contains a numpy array per spin channel"
    def __new__(cls,*args,**kwargs): #logdata,ikpt=0,kpt=(0.0,0.0,0.0):
        "Takes the data from the logfile and convert it"
        evs=[[],[]]
        logdata=kwargs.get("logdata", args[0])
        ikpt=kwargs.get("ikpt", 1 if len(args) < 2 else args[1])
        for ev in logdata:
            occ=get_ev(ev,['e_occ','e_occupied'],ikpt=ikpt)
            vrt=get_ev(ev,['e_vrt','e_virt'],ikpt=ikpt)
            eigen=occ or vrt
            if not eigen: eigen=get_ev(ev,ikpt=ikpt)
            if not eigen: continue
            for i,e in enumerate(eigen):
                if e: evs[i].append(e)
        #here we can create the numpy array
        #data=numpy.full((2 if len(evs[1])>0 else 1, max(map(len,evs))), numpy.nan)
        data=numpy.ndarray.__new__(cls,shape=(2 if len(evs[1])>0 else 1, max(map(len,evs))),dtype=numpy.float)
        data.fill(numpy.nan)
        data[0,:len(evs[0])]=evs[0]
        if len(evs[1])>0: data[1,:len(evs[1])]=evs[1]
        data.info=(map(len,evs))
        return data
    def __init__(self,*args,**kwargs):
        ikpt=kwargs.get("ikpt", 1 if len(args) < 2 else args[1])
        kpt=kwargs.get("kpt", (0.,0.,0.) if len(args) < 3 else args[2])
        kwgt=kwargs.get('kwgt',1.0)
        self.set_kpt(ikpt,kpt,kwgt)
    def set_kpt(self,ikpt,kpt,kwgt=1.0):      
        if not isinstance(ikpt, int ): raise TypeError('ikpt should be a integer')
        if len(kpt) !=3: raise TypeError('kpt should be a object of len 3')
        self.ikpt=ikpt
        self.kpt=kpt
        self.kwgt=kwgt

#definition of the density of state class
class DoS():
    def __init__(self,bandarrays=None,energies=None,evals=None,units='eV',label='1',sigma=0.2/AU_eV,npts=2500,fermi_level=None,norm=1.0):
        "Extract a quantity which is associated to the DoS, that can be plotted"
        import numpy as np
        self.ens=[]
        self.labels=[]
        self.norms=[]
        if bandarrays is not None: self.append_from_bandarray(bandarrays,label)
        if evals is not None: self.append_from_dict(evals,label)
        if energies is not None: self.append(energies,label=label,units=units,norm=norm)
        self.sigma=self.conversion_factor(units)*sigma
        self.npts=npts
        self.fermi_level(fermi_level,units=units)
    def append_from_bandarray(self,bandarrays,label):
        import numpy as np
        kptlists=[[],[]]
        for orbs in bandarrays:
            for ispin,norbs in enumerate(orbs.info):
                if norbs==0: continue
                kptlists[ispin]+=list(orbs[ispin,:norbs]*orbs.kwgt)
        #print 'kpt',kptlists
        for ispin,norbs in enumerate(kptlists):
            if len(norbs)==0: continue
            self.append(np.array(norbs),label=label,units='AU',norm=(1.0-2*ispin))
    def append_from_dict(self,evals,label):
        import numpy as np
        "Get the energies from the different flavours given by the dict"
        evs=[[],[]]
        ef=None
        for ev in evals:
            occ=self.get_ev(ev,['e_occ','e_occupied'])
            if occ: ef=max(occ)
            vrt=self.get_ev(ev,['e_vrt','e_virt'])
            eigen=False
            if occ: eigen=occ
            if vrt: eigen=vrt
            if not eigen: eigen=self.get_ev(ev)
            if not occ and not vrt and eigen: ef=max(eigen)
            if not eigen: continue
            for i,e in enumerate(eigen):
                if e: evs[i].append(e)
        for i,energs in enumerate(evs):
            if len(energs)==0: continue
            self.append(np.array(energs),label=label,units='AU',norm=1.0-2.0*i)
        if ef: self.fermi_level(ef,units='AU')
    def get_ev(self,ev,keys=None):
        "Get the correct list of the energies for this eigenvalue"
        res=False
        if keys is None:
            ener=ev.get('e')
            spin=ev.get('s')
            if ener and spin==1:
                res=[ener]
            elif ener and spin==-1:
                res=[None,ener]
        else:
            for k in keys:
                if k in ev:
                    res=ev[k]
                    if type(res) != type([]): res=[res]
                    break
        return res
    def append(self,energies,label=None,units='eV',norm=1.0):
        self.ens.append(self.conversion_factor(units)*energies)
        if label is not None:
            self.labels.append(label)
        else:
            self.labels.append(str(len(self.labels)+1))
        self.norms.append(norm)
    def conversion_factor(self,units):
        if units == 'AU':
            fac = AU_eV
        elif units == 'eV':
            fac=1.0
        else:
            raise 'Unrecognized units (',unit,')'
        return fac
    def fermi_level(self,fermi_level,units='eV'):
        if fermi_level is not None:
            self.ef=fermi_level*self.conversion_factor(units)
    def range(self,npts=None):
        import numpy as np
        if npts is None: npts=self.npts
        e_min=1.e100
        e_max=-1.e100
        for dos in self.ens:
            e_min = min(e_min,np.min(dos) - 0.05*(np.max(dos) - np.min(dos)))
            e_max = max(e_max,np.max(dos) + 0.05*(np.max(dos) - np.min(dos)))
        return np.arange( e_min, e_max, (e_max-e_min)/npts )
    def curve(self,dos,norm,sigma=None):
        import numpy as np
        if sigma is None: sigma=self.sigma
        nrm=np.sqrt(2.0*np.pi)*sigma/norm
        dos_g = []
        for e_i in self.range():
            dos_g.append(np.sum(np.exp( - (e_i - dos[:])**2 / (2.0 * sigma**2))/nrm)) #Append data corresponding to each energy grid
        return np.array(dos_g)
    def dump(self,sigma=None):
        "For Gnuplot"
        if sigma is None: sigma=self.sigma
        data=[self.curve(dos,norm=self.norms[i],sigma=sigma) for i,dos in enumerate(self.ens)]
        for i,e in enumerate(self.range()):
            print e,' '.join(map(str,[d[i] for d in data]))
    def plot(self,sigma=None):
        import matplotlib.pyplot as plt
        from matplotlib.widgets import Slider#, Button, RadioButtons
        if sigma is None: sigma=self.sigma
        self.fig, self.ax1 = plt.subplots()
        self.plotl=[]
        for i,dos in enumerate(self.ens):
            self.plotl.append(self.ax1.plot(self.range(),self.curve(dos,norm=self.norms[i],sigma=sigma),label=self.labels[i]))
        plt.xlabel('Energy [eV]', fontsize=18)
        plt.ylabel('DoS', fontsize=18)
        if self.ef is not None:
            self.ax1.annotate('Fermi level', xy=(self.ef,2),
                    xytext=(self.ef, 10),
                arrowprops=dict(facecolor='white', shrink=0.05),
            )
        if len(self.labels) > 1: plt.legend()
        axcolor = 'lightgoldenrodyellow'
        axsigma = plt.axes([0.2, 0.93, 0.65, 0.03], axisbg=axcolor)
        self.ssig = Slider(axsigma, 'Smearing', 0.0, 0.4, valinit=sigma)
        self.ssig.on_changed(self.update)
        plt.show()
    def update(self,val):
        sig = self.ssig.val
        for i,dos in enumerate(self.ens):
            self.plotl[i][0].set_ydata(self.curve(dos,norm=self.norms[i],sigma=sig))
        self.ax1.relim()
        self.ax1.autoscale_view()
        self.fig.canvas.draw_idle()

class BrillouinZone():
    def __init__(self,astruct,mesh,evals):
        import spglib,numpy
	self.lattice=numpy.diag(astruct['Cell'])
	pos=[ [ a/b for a,b in zip(at.values()[0], astruct['Cell'])]for at in astruct['Positions']]
        atoms=[ at.keys()[0] for at in astruct['Positions']]
	ianames,iatype=numpy.unique(atoms,return_inverse=True) #[1,]*4+[2,]*4 #we should write a function for the iatype
        #print 'iatype', iatype
	cell=(self.lattice,pos,iatype)
	print 'spacegroup',spglib.get_spacegroup(cell, symprec=1e-5)
        #then define the pathes and special points
        import ase.dft.kpoints as ase
        #we should adapt the 'cubic'
        cell_tmp=astruct['Cell']
        #print 'cell',cell_tmp,numpy.allclose(cell_tmp,[cell_tmp[0],]*len(cell_tmp))
        if numpy.allclose(cell_tmp,[cell_tmp[0],]*len(cell_tmp)):
            lattice_string='cubic'
        else:
            lattice_string='orthorhombic'
        print 'Lattice found:',lattice_string
        self.special_points=ase.get_special_points(lattice_string, self.lattice, eps=0.0001)
        self.special_paths=ase.special_paths[lattice_string]
	#dataset = spglib.get_symmetry_dataset(cell, symprec=1e-3)
	#print dataset
	#the shift has also to be put if present
	mapping, grid = spglib.get_ir_reciprocal_mesh(mesh, cell, is_shift=[0, 0, 0])
	lookup=[]
	for ikpt in numpy.unique(mapping):
	    ltmp=[]
	    for m, g in zip(mapping, grid):
	        if m==ikpt:
	            ltmp.append(g)
	    lookup.append(ltmp)
        #print 'irreductible k-points',len(lookup)
	coords=numpy.array(grid, dtype = numpy.float)/mesh
	#print grid #[ x+mesh[0]*y+mesh[0]*mesh[1]*z for x,y,z in grid]
	#brillouin zone
        kp=numpy.array([ k.kpt for k in evals])
	ourkpt=numpy.rint(kp*(numpy.array(mesh))).astype(int)
	bz=numpy.ndarray(tuple(mesh) + (evals[0].size, ),dtype=float)
	#print bz
	shift=(numpy.array(mesh)-1)/2
	for ik in lookup:
	    irrk=None
	    for orbs,bzk in zip(evals,ourkpt):
	        for kt in ik:
	            if (bzk==kt).all():
	                irrk=orbs
	                #print 'hello',orbs.kpt,kt
	                break 
	        if irrk is not None: break
	    if irrk is None:
	        print 'error in ik',ik
	        print 'our',ourkpt
	        print 'spglib',grid
	        print 'mapping',mapping
	    for kt in ik:
	        r=kt+shift
	        #print irrk.shape, bz.shape
	        bz[r[0],r[1],r[2],:]=irrk.reshape(irrk.size)
        #duplicate coordinates for the interpolation
        bztmp=bz.reshape((mesh[0]*mesh[1]*mesh[2], -1))
        #print bztmp
        bztot=numpy.ndarray((7,bztmp.shape[0],bztmp.shape[1]))
        bztot[0,:,:]=bztmp
        ctot=numpy.ndarray((7,coords.shape[0],coords.shape[1]))
        ctot[0,:,:]=coords
        for i,sh in enumerate([[-1,0,0],[1,0,0],[0,-1,0],[0,1,0],[0,0,-1],[0,0,1]]):
            bztot[i,:,:]=bztmp
            ctot[i,:,:]=coords+sh
            #print 'coors',coords,coords+[1.0,0,0]
        bztot=bztot.reshape((7*bztmp.shape[0], -1))
        ctot=ctot.reshape((7*coords.shape[0], -1))
        import scipy.interpolate.interpnd as interpnd
        self.interpolator=interpnd.LinearNDInterpolator(ctot,bztot)
    def plot(self,path=None,npts=50):
        import ase.dft.kpoints as ase
        path_tmp=path
        if path is None: path_tmp=self.special_paths[0]+['G',]
        ppath,xaxis,xlabel=ase.get_bandpath([self.special_points[p] for p in path_tmp], self.lattice,npts)
        symbols = [t.replace('G', '$\Gamma$') for t in path_tmp]
        toto=self.interpolator(ppath)
        import matplotlib.pyplot as plt
        #print toto.min(),toto.max()
        for b in toto.transpose():
            plt.plot(xaxis, b)#, 'b-')
        plt.axhline(0,color='k',linestyle='--')
        for p in xlabel:
            plt.axvline(p,color='k',linestyle='-')
            plt.xticks(xlabel, symbols)
        plt.show()



if __name__ == "__main__":
    import numpy as np
    energies=np.array([-0.815924953235059, -0.803163374736654, -0.780540200987971, -0.7508806541364, -0.723626807289917, -0.714924448617026, -0.710448085701742, -0.68799028016451, -0.67247569974853, -0.659038909236607, -0.625396293324399, -0.608009041659988, -0.565337910777367, -0.561250536074343, -0.551767438323268, -0.541295070404525, -0.532326667587434, -0.515961980147107, -0.474601108285518, -0.473408476151224, -0.46509070541069, -0.445709086452906, -0.433874403837837, -0.416121660651406, -0.407871082254237, -0.406123490618786, -0.403004188319382, -0.38974739285104, -0.380837488456638, -0.375163102271681, -0.375007771592681, -0.367898783582561, -0.367518948507212, -0.359401585874402, -0.358189406008502, -0.354517727598174, -0.334286389724978, -0.332921810616845, -0.315466259109401, -0.308028853904577, -0.29864142362141, -0.294024743731349, -0.292104129933301, -0.285165738729842, -0.28419932605141, -0.267399999874122, -0.259487769142101, -0.239899780812716, -0.224858003804207, -0.20448050758473, -0.164155133452971, -0.117617164459898, -0.0717938081884113, -0.0526986239898579, -0.0346031190163735, -0.0167949342608791, -0.0135168064347152, -0.0102971895842409, 0.00759271179427191, 0.00974950976249545, 0.010176021051287, 0.0217652761059223, 0.0239924727094222, 0.0413057846713024, 0.0422334333464529, 0.0459150454793617, 0.0517637894860314])
    dos=DoS(energies,fermi_level=-0.1)
    dos.append(0.2+energies)
    dos.dump(sigma=0.01)
    dos.plot()
