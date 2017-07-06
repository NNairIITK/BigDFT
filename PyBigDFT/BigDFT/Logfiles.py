#This module needs: yaml, futile, matplotlib, numpy, BZ, DoS
import yaml
from futile.Utils import write

EVAL = "eval"
SETUP = "let"
INITIALIZATION = "globals"

PATH='path'
PRINT='print'
GLOBAL='global'
FLOAT_SCALAR='scalar'

PRE_POST = [EVAL, SETUP, INITIALIZATION]

#Builtin pathes to define the search paths
BUILTIN={
    'nat': {PATH: [ ['Atomic System Properties','Number of atoms']],
                 PRINT: "Number of Atoms", GLOBAL: True},
    "energy": {PATH: [["Last Iteration", "FKS"],["Last Iteration", "EKS"], ["Energy (Hartree)"]], 
                    PRINT: "Energy", GLOBAL: False},
    "fermi_level": {PATH:[["Ground State Optimization", -1, "Fermi Energy"], ["Ground State Optimization", -1, "Hamiltonian Optimization", -1, "Subspace Optimization", "Fermi Energy"]],
                         PRINT: True, GLOBAL: False},
    'astruct': {PATH: [ ['Atomic structure']]},
    'evals': {PATH: [ ["Complete list of energy eigenvalues"], [ "Ground State Optimization", -1, "Orbitals"],
                    ["Ground State Optimization",-1,"Hamiltonian Optimization",-1,"Subspace Optimization","Orbitals"] ]},
    'kpts': {PATH: [["K points"]],
             PRINT: False, GLOBAL: True},
    'gnrm_cv': {PATH: [["dft","gnrm_cv"]], PRINT: "Convergence criterion on Wfn. Residue", GLOBAL: True},
    'kpt_mesh': {PATH:[['kpt','ngkpt']], PRINT: True, GLOBAL: True},
    'forcemax': {PATH: [["Geometry","FORCES norm(Ha/Bohr)","maxval"],['Clean forces norm (Ha/Bohr)','maxval']],
                 PRINT: "Max val of Forces"},
    'pressure': {PATH: [['Pressure','GPa']], PRINT: True},
    'forces': {PATH: [['Atomic Forces (Ha/Bohr)']]},
    'forcemax_cv': {PATH: [['geopt','forcemax']], PRINT: 'Convergence criterion on forces', GLOBAL: True, FLOAT_SCALAR: True},
    'force_fluct': {PATH:[["Geometry","FORCES norm(Ha/Bohr)","fluct"]], PRINT: "Threshold fluctuation of Forces"},
    'magnetization': {PATH:[[ "Ground State Optimization", -1, "Total magnetization"],
                            ["Ground State Optimization",-1,"Hamiltonian Optimization",-1,"Subspace Optimization","Total magnetization"]],
                      PRINT: "Total magnetization of the system"},
    'support_functions': {PATH: [["Gross support functions moments",'Multipole coefficients','values']]},
    'electrostatic_multipoles': {PATH: [['Multipole coefficients','values']]},
    'sdos': {PATH: [['SDos files']], GLOBAL: True},
    'symmetry': {PATH: [ ['Atomic System Properties','Space group']], 
                 PRINT: "Symmetry group", GLOBAL: True}}

def get_logs(files,select_document=None):
    """
    Return a list of loaded logfiles from files, which is a list
    of paths leading to logfiles.
    
    Optional arguments:
    - select_document:
    """
    from futile import Yaml
    logs=[]
    for filename in files:
        logs+=Yaml.load(filename,doc_lists=True,safe_mode=True)
    return logs

##
##def get_logs(files,safe_mode=False,select_document=None):
##    """
##    Return a list of loaded logfiles from files, which is a list
##    of paths leading to logfiles.
##    
##    Optional arguments:
##    - safe_mode:
##    - select_document:
##    """
##    logs=[]
##    for filename in files:
##      rawfile=open(filename, "r").read()
##      try:
##         logs+=[yaml.load(rawfile, Loader = yaml.CLoader)]
##      except Exception,e:
##         print 'WARNING: More than one document are present',e
##         if safe_mode or select_document is not None:
##             documents=rawfile.split('---\n')
##             print 'Safe mode, Found',len(documents),'documents,try loading them separately'
##             actual_doc=-1
##             for i,raw_doc in enumerate(documents):
##                 if len(raw_doc)==0: continue
##                 actual_doc+=1
##                 if select_document is not None and actual_doc not in select_document: continue
##                 try:
##                     logs.append(yaml.load(raw_doc,Loader=yaml.CLoader))
##                     print 'Document',i,'...loaded.'
##                 except Exception,f:
##                     print 'Document',i,'...NOT loaded.'
##                     print f
##                     #logs+=[None]
##                     #print "warning, skipping logfile",filename
##         else:
##             try: 
##                 from futile import Yaml
##                 test=Yaml.YamlDB(rawfile)
##                 for a in range(len(test)):
##                     # we should use another representation for the logfile, to be changed
##                     lg=dict(test[a]) 
##                     if lg is not None: logs+=[lg]
##                 #logs+=yaml.load_all(rawfile, Loader = yaml.CLoader)
##             except Exception,e:
##                 print e
##                 print 'WARNING: Usual loading of the document have some errors, some documents might not be there'
##                 print 'Consider to put safe_mode=True'
##    return logs
##

def floatify(scalar):
    """Useful to make float from strings compatible from fortran"""
    import numpy
    if isinstance(scalar,str):
        return float(scalar.replace('d','e').replace('D','E'))
    else:
        return scalar


# This is a tentative function written to extract information from the runs
def document_quantities(doc,to_extract):
    """Extract information from the runs."""
    analysis={}
    for quantity in to_extract:
      if quantity in PRE_POST: continue
      #follow the levels indicated to find the quantity
      field=to_extract[quantity]
      if type(field) is not type([]) is not type({}) and field in BUILTIN:
          paths=BUILTIN[field][PATH]
      else:
          paths=[field]
      #now try to find the first of the different alternatives
      for path in paths:
        #print path,BUILTIN,BUILTIN.keys(),field in BUILTIN,field
        value=doc
        for key in path:
          #as soon as there is a problem the quantity is null
          try:
            value=value[key]
          except:
            value=None
            break
        if value is not None: break
      analysis[quantity]=value
    return analysis


def perform_operations(variables,ops,debug=False):
    """Perform operations given by 'ops'.
    'variables' is a dictionary of variables i.e. key=value."""
##    glstr=''
##    if globs is not None:
##        for var in globs:
##            glstr+= "global "+var+"\n"
##        if debug: print '###Global Strings: \n',glstr
##    #first evaluate the given variables
    for key in variables:
        command=key+"="+str(variables[key])
        if debug: print(command)
        exec(command)
        #then evaluate the given expression
    if debug: print(ops)
    #exec(glstr+ops, globals(), locals())
    exec(ops, globals(), locals())


def process_logfiles(files,instructions,debug=False):
    """Process the logfiles in files with the dictionary 'instructions'."""
    import sys
    glstr='global __LAST_FILE__ \n'
    glstr+='__LAST_FILE__='+str(len(files))+'\n'
    if INITIALIZATION in instructions:
        for var in instructions[INITIALIZATION]:
            glstr+= "global "+var+"\n"
            glstr+= var +" = "+ str(instructions[INITIALIZATION][var])+"\n"
            #exec var +" = "+ str(instructions[INITIALIZATION][var])
    exec(glstr, globals(), locals())
    for f in files:
        sys.stderr.write("#########processing "+f+"\n")
        datas=get_logs([f])
        for doc in datas:
            doc_res=document_quantities(doc,instructions)
            #print doc_res,instructions
            if EVAL in instructions: perform_operations(doc_res,instructions[EVAL],debug=debug)


def find_iterations(log):
    """Identify the different block of the iterations of the wavefunctions optimization."""
    import numpy
    for itrp in log['Ground State Optimization']:
        rpnrm=[]
        for itsp in itrp['Hamiltonian Optimization']:
            gnrm_sp=[]
            for it in itsp['Subspace Optimization']['Wavefunctions Iterations']:
                if 'gnrm' in it: gnrm_sp.append(it['gnrm'])
            rpnrm.append(numpy.array(gnrm_sp))
    rpnrm=numpy.array(rpnrm)
    return rpnrm


def plot_wfn_convergence(wfn_it,gnrm_cv):
    """Plot the convergence of the wavefunction coming from the find_iterations function.""" 
    import matplotlib.pyplot as plt
    import numpy 
    plt.semilogy(numpy.ravel(wfn_it))
    plt.axhline(gnrm_cv,color='k',linestyle='--')
    it=0
    for itrp in wfn_it:
        it+=len(itrp)
        plt.axvline(it,color='k',linestyle='--')



class Logfile():
    """Import a Logfile from a filename in yaml format, a list of filenames,
        an archive (compressed tar file), aidctionaory or a list od dictionaries:
        l = Logfile('one.yaml','two.yaml')
        l = Logfile(archive='calc.tgz')
        l = Logfile(archive='calc.tgz',member='one.yaml')
        l = Logfile(dictionary=dict)
        l = Logfile(dictionary=[dict1, dict2])"""
    def __init__(self,*args,**kwargs):
        dicts = []
        #Read the dictionary kwargs
        arch = kwargs.get("archive")
        member = kwargs.get("member")
        label = kwargs.get("label")
        dictionary=kwargs.get("dictionary")
        if arch:
            #An archive is detected
            import tarfile,os
            from futile import Yaml
            tar = tarfile.open(arch)
            members = [ tar.getmember(member) ] if member else tar.getmembers()
            #print members
            for memb in members:
                f = tar.extractfile(memb)
                dicts.append(Yaml.load(stream=f.read()))
                #dicts[-1]['label'] = memb.name #Add the label (name of the file)
            srcdir=os.path.dirname(arch)
        elif dictionary:
            #Read the dictionary or a list of dictionaries or from a generator
            dicts = dictionary if isinstance(dictionary,list) else [d for d in dictionary]
            srcdir=''
        elif args:
            import os
            #Read the list of files (member replaces load_only...)
            dicts=get_logs(args,select_document=member)
            label = label if label else args[0]
            srcdir=os.path.dirname(args[0])
        #Set the label
        self.label=label
        self.srcdir=os.path.abspath('.' if srcdir == '' else srcdir)
        if not dicts:
            raise ValueError("No log information provided.")
        #Initialize the logfile with the first document
        self._initialize_class(dicts[0])
        if len(dicts)>1:
            #first initialize the instances with the previous logfile such as to provide the
            #correct information
            #(we should however decide what to do if some run did not converged)
            self._instances=[]
            for i,d in enumerate(dicts):
                #label=d.get('label','log'+str(i))
                label = 'log'+str(i)
                dtmp=dicts[0]
                instance=Logfile(dictionary=dtmp,label=label)
                #now update the instance with the other value
                instance._initialize_class(d)
                self._instances.append(instance)
            #then we should find the best values for the dictionary
            print('Found',len(self._instances),'different runs')	
            import numpy
            #Initialize the class with the dictionary corresponding to the lower value of the energy
            ens=[(l.energy if hasattr(l,'energy') else 1.e100) for l in self._instances] 
            self.reference_log=numpy.argmin(ens)
            #print 'Energies',ens
            self._initialize_class(dicts[self.reference_log])
    #
    def __getitem__(self,index):
        if hasattr(self,'_instances'):
            return self._instances[index]
        else:
            print('index not available')
            raise
    #
    def __str__(self):
        """Display short information about the logfile"""
        return self._print_information()
    #
    def __len__(self):
        if hasattr(self,'_instances'):
            return len(self._instances)
        else:
            return 0 #single point run
    def _initialize_class(self,d):
        import numpy
        self.log=d
        #here we should initialize different instances of the logfile class again
        sublog=document_quantities(self.log,{val: val for val in BUILTIN})
        for att in sublog:
            val=sublog[att]
            if val is not None: 
                val_tmp=floatify(val) if BUILTIN[att].get(FLOAT_SCALAR) else val
                setattr(self,att,val_tmp)
            elif hasattr(self,att) and not BUILTIN[att].get(GLOBAL):
                delattr(self,att)
        #then postprocess the particular cases
        if not hasattr(self,'fermi_level') and hasattr(self,'evals'):
            self._fermi_level_from_evals(self.evals)

        if hasattr(self,'kpts'):
            self.nkpt=len(self.kpts)
            if hasattr(self,'evals'): self.evals=self._get_bz(self.evals,self.kpts)
            if hasattr(self,'forces') and hasattr(self,'astruct'): 
                self.astruct.update({'forces': self.forces})
                delattr(self,'forces')
        elif hasattr(self,'evals'):
            import BZ
            self.evals=[BZ.BandArray(self.evals),]
        if hasattr(self,'sdos'):
            import os
            #load the different sdos files
            sd=[]
            for f in self.sdos:
                try:
                    data=numpy.loadtxt(os.path.join(self.srcdir,f))
                except:
                    data=None
                if data is not None:
                    xs=[]
                    ba=[[],[]]
                    for line in data:
                        xs.append(line[0])
                        ss=self._sdos_line_to_orbitals(line)
                        for ispin in [0,1]:
                            ba[ispin].append(ss[ispin])
                    sd.append({'coord':xs,'dos':ba})
                else:
                    sd.append(None)
            self.sdos=sd
    #
    def _fermi_level_from_evals(self,evals):
        import numpy
        #this works when the representation of the evals is only with occupied states
        #write('evals',self.evals)
        fl=None
        for iorb,ev in enumerate(evals):
            e=ev.get('e')
            if e is not None:
                fref=ev['f'] if iorb==0 else fref
                fl=e
                if ev['f']<0.5*fref: break
            e=ev.get('e_occ',ev.get('e_occupied'))
            if e is not None: fl=e
            e=ev.get('e_vrt',ev.get('e_virt'))
            if e is not None: break
        self.fermi_level=fl
    #
    def _sdos_line_to_orbitals_old(self,sorbs):
        import BZ
        evals=[]
        iorb=1
        #renorm=len(xs)
        #iterate on k-points
        kpts=self.kpts if hasattr(self,'kpts') else [{'Rc':[0.0,0.0,0.0],'Wgt':1.0}]
        for i,kp in enumerate(kpts):
            ev=[]
            #iterate on the subspaces of the kpoint
            for ispin,norb in enumerate(self.evals[0].info):
                for iorbk in range(norb):
                    #renorm postponed
                    ev.append({'e':sorbs[iorb+iorbk],'s':1-2*ispin,'k':i+1})
                    #ev.append({'e':np.sum([ so[iorb+iorbk] for so in sd]),'s':1-2*ispin,'k':i+1})
                iorb+=norb
            evals.append(BZ.BandArray(ev,ikpt=i+1,kpt=kp['Rc'],kwgt=kp['Wgt']))
        return evals
    #
    def _sdos_line_to_orbitals(self,sorbs):
        import BZ,numpy as np
        evals=[]
        iorb=1
        sdos=[[],[]]
        for ikpt,band in enumerate(self.evals):
            sdoskpt=[[],[]]
            for ispin,norb in enumerate(band.info):
                if norb==0: continue
                bands=band[ispin]
                for i in range(norb):
                    val=sorbs[iorb]
                    e=bands[i]
                    #val/=e #not needed anymore
                    iorb+=1
                    sdoskpt[ispin].append(val)
                sdos[ispin].append(np.array(sdoskpt[ispin]))
        return sdos
    #
    def _get_bz(self,ev,kpts):
        """Get the Brillouin Zone."""
        evals=[]
        import BZ
        for i,kp in enumerate(kpts):
            evals.append(BZ.BandArray(ev,ikpt=i+1,kpt=kp['Rc'],kwgt=kp['Wgt']))
        return evals
    #
    def get_dos(self,label=None,npts=2500):
        """Get the density of states from the logfile."""
        import DoS
        #reload(DoS)
        lbl=self.label if label is None else label
        sdos=self.sdos if hasattr(self,'sdos') else None
        return DoS.DoS(bandarrays=self.evals,label=lbl,units='AU',fermi_level=self.fermi_level,npts=npts,sdos=sdos)
    #
    def get_brillouin_zone(self):
        """Returns an instance of the BrillouinZone class, useful for band structure."""
        import BZ
        if self.nkpt==1: 
            print('WARNING: Brillouin Zone plot cannot be defined properly with only one k-point')
            #raise
        mesh=self.kpt_mesh
        if isinstance(mesh,int): mesh=[mesh,]*3
        if self.astruct['Cell'][1]==float('inf'): mesh[1]=1
        return BZ.BrillouinZone(self.astruct,mesh,self.evals,self.fermi_level)
    #
    def wfn_plot(self):
        """Plot the wavefunction convergence."""
        wfn_it=find_iterations(self.log)
        plot_wfn_convergence(wfn_it,self.gnrm_cv)
    #
    def geopt_plot(self):
        """For a set of logfiles construct the convergence plot if available."""
        import numpy
        energies=[]
        forces=[]
        ferr=[]
        if not hasattr(self,'_instances'): 
            print('ERROR: No geopt plot possible, single point run')
            return
        for l in self._instances:
            if hasattr(l,'forcemax') and hasattr(l,'energy'):
                forces.append(l.forcemax)
                energies.append(l.energy-self.energy)
                ferr.append(0.0 if not hasattr(l,'force_fluct') else self.force_fluct)
        if len(forces) > 1:
            import matplotlib.pyplot as plt
            plt.errorbar(energies, forces,yerr=ferr, fmt='.-',label=self.label)
            plt.legend(loc='upper right')
            plt.loglog()
            plt.xlabel('Energy - min(Energy)')
            plt.ylabel('Forcemax')
            if hasattr(self,'forcemax_cv'): plt.axhline(self.forcemax_cv,color='k',linestyle='--')
            plt.show()
        else:
            print('No plot necessary, less than two points found')
    #
    #
    def _print_information(self):
        """Display short information about the logfile (used by str)."""
        import yaml,numpy
        summary=[{'Atom types': 
                  numpy.unique([ at.keys()[0] for at in self.astruct['Positions']]).tolist()},
                 {'Cell': 
                  self.astruct.get('Cell','Free BC')}]
        #normal printouts in the document, according to definition
        for field in BUILTIN:
            name=BUILTIN[field].get(PRINT)
            if name==True: name=field
            if not name or not hasattr(self,field): continue
            summary.append({name: getattr(self,field)})
        if hasattr(self,'evals'): summary.append({'No. of KS orbitals per k-point': self.evals[0].info})
        return yaml.dump(summary,default_flow_style=False)

if __name__ == "__main__":
    #Create a logfile: should give an error (ValueError: No log information provided.)
    l = Logfile()
