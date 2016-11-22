import numpy

AU_eV=27.31138386

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
        "Important for kpoints DOS"
        import numpy as np
        for jspin in range(2):
            lbl= 'up' if jspin==0 else 'dw'
            kptlists=[[],[]]
            for orbs in bandarrays:
                for ispin,norbs in enumerate(orbs.info):
                    if norbs==0 or ispin !=jspin: continue
                    #energy values
                    kptlists[0].append(orbs[ispin,:norbs])
                    #normalization
                    kptlists[1].append(orbs.kwgt*(1.0-2*ispin))
                #print 'kpt',kptlists
            self.append(np.array(kptlists[0]),label=label+lbl,units='AU',
                        norm=np.array(kptlists[1]))
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
            ddos=np.ravel(dos)
            if len(ddos)==0: continue
            e_min = min(e_min,np.min(ddos) - 0.05*(np.max(ddos) - np.min(ddos)))
            e_max = max(e_max,np.max(ddos) + 0.05*(np.max(ddos) - np.min(ddos)))
        return np.arange( e_min, e_max, (e_max-e_min)/npts )
    def curve(self,dos,norm,sigma=None):
        import numpy as np
        if sigma is None: sigma=self.sigma
        nrm=np.sqrt(2.0*np.pi)*sigma/norm
        dos_g = []
        for e_i in self.range():
            if len(dos.shape)==2:
                nkpt=dos.shape[0]
                value=0.0
                for ikpt in range(nkpt):
                    value+=np.sum(np.exp( - (e_i - dos[ikpt,:])**2 / (2.0 * sigma**2))/nrm[ikpt])
            else:
                value=np.sum(np.exp( - (e_i - dos[:])**2 / (2.0 * sigma**2))/nrm)
            dos_g.append(value) #Append data corresponding to each energy grid
        return np.array(dos_g)
    def dump(self,sigma=None):
        "For Gnuplot"
        if sigma is None: sigma=self.sigma
        data=[self.curve(dos,norm=self.norms[i],sigma=sigma) for i,dos in enumerate(self.ens)]
        for i,e in enumerate(self.range()):
            print e,' '.join(map(str,[d[i] for d in data]))
    def plot(self,sigma=None,legend=False):
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
            plt.axvline(self.ef,color='k',linestyle='--')
            #self.ax1.annotate('Fermi level', xy=(self.ef,2),
            #        xytext=(self.ef, 10),
            #    arrowprops=dict(facecolor='white', shrink=0.05),
            #)
        if len(self.labels) > 1 or legend: plt.legend()
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

if __name__ == "__main__":
    import numpy as np
    energies=np.array([-0.815924953235059, -0.803163374736654, -0.780540200987971, -0.7508806541364, -0.723626807289917, -0.714924448617026, -0.710448085701742, -0.68799028016451, -0.67247569974853, -0.659038909236607, -0.625396293324399, -0.608009041659988, -0.565337910777367, -0.561250536074343, -0.551767438323268, -0.541295070404525, -0.532326667587434, -0.515961980147107, -0.474601108285518, -0.473408476151224, -0.46509070541069, -0.445709086452906, -0.433874403837837, -0.416121660651406, -0.407871082254237, -0.406123490618786, -0.403004188319382, -0.38974739285104, -0.380837488456638, -0.375163102271681, -0.375007771592681, -0.367898783582561, -0.367518948507212, -0.359401585874402, -0.358189406008502, -0.354517727598174, -0.334286389724978, -0.332921810616845, -0.315466259109401, -0.308028853904577, -0.29864142362141, -0.294024743731349, -0.292104129933301, -0.285165738729842, -0.28419932605141, -0.267399999874122, -0.259487769142101, -0.239899780812716, -0.224858003804207, -0.20448050758473, -0.164155133452971, -0.117617164459898, -0.0717938081884113, -0.0526986239898579, -0.0346031190163735, -0.0167949342608791, -0.0135168064347152, -0.0102971895842409, 0.00759271179427191, 0.00974950976249545, 0.010176021051287, 0.0217652761059223, 0.0239924727094222, 0.0413057846713024, 0.0422334333464529, 0.0459150454793617, 0.0517637894860314])
    dos=DoS(energies,fermi_level=-0.1)
    dos.append(0.2+energies)
    dos.dump(sigma=0.01)
    dos.plot()
