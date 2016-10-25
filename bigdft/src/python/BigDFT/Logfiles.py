import yaml

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
    "fermi_level": {PATH:[["Ground State Optimization", -1, "Fermi Energy"]], 
                         PRINT: True, GLOBAL: False},
    'astruct': {PATH: [ ['Atomic structure']]},
    'evals': {PATH: [ ["Complete list of energy eigenvalues"], [ "Ground State Optimization", -1, "Orbitals"],
                    ["Ground State Optimization",-1,"Hamiltonian Optimization",-1,"Subspace Optimization","Orbitals"] ]},
    'kpts': {PATH: [["K points"]],
             PRINT: False, GLOBAL: True},
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
    'symmetry': {PATH: [ ['Atomic System Properties','Space group']], 
                 PRINT: "Symmetry group", GLOBAL: True}}

def get_log(f):
    "Transform a logfile into a python dictionary"
    return yaml.load(open(f, "r").read(), Loader = yaml.CLoader)

def get_logs(files,safe_mode=False,select_document=None):
   logs=[]
   for filename in files:
     rawfile=open(filename, "r").read()
     try:
        logs+=[yaml.load(rawfile, Loader = yaml.CLoader)]
     except:
        if safe_mode or select_document is not None:
            documents=rawfile.split('---\n')
            print 'Safe mode, Found',len(documents),'documents,try loading them separately'
            #print 'first',documents[0]
            #print 50*'X'
            #print 'last',documents[-1]
            actual_doc=-1
            for i,raw_doc in enumerate(documents):
                if len(raw_doc)==0: continue
                actual_doc+=1
                if select_document is not None and actual_doc not in select_document: continue
                try:
                    logs.append(yaml.load(raw_doc,Loader=yaml.CLoader))
                    print 'Document',i,'...loaded.'
                except Exception,f:
                    print 'Document',i,'...NOT loaded.'
                    print f
                    #logs+=[None]
                    #print "warning, skipping logfile",filename
        else:
            try: 
                logs+=yaml.load_all(rawfile, Loader = yaml.CLoader)
            except Exception,e:
                print e
                print 'WARNING: Usual loading of the document have some errors, some documents might not be there'
                print 'Consider to put safe_mode=True'
   return logs

def floatify(scalar):
    "Useful to make float from strings compatible from fortran"
    import numpy
    if isinstance(scalar,str):
        return float(scalar.replace('d','e').replace('D','E'))
    else:
        return scalar

# this is a tentative function written to extract information from the runs
def document_quantities(doc,to_extract):
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
##    glstr=''
##    if globs is not None:
##        for var in globs:
##            glstr+= "global "+var+"\n"
##        if debug: print '###Global Strings: \n',glstr
##    #first evaluate the given variables
    for key in variables:
        command=key+"="+str(variables[key])
        if debug: print command
        exec(command)
        #then evaluate the given expression
    if debug: print ops
    #exec(glstr+ops, globals(), locals())
    exec(ops, globals(), locals())

def process_logfiles(files,instructions,debug=False):
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


class Logfile():
    def __init__(self,filename=None,dictionary=None,filename_list=None,label=None,load_only=None):
        "Import a Logfile from a filename in yaml format"
        filelist=None
        self.label=label
        if filename is not None: 
            if self.label is None: self.label=filename
            filelist=[filename]
        elif filename_list is not None:
            if self.label is None: self.label=filename_list[0]
            filelist=filename_list
        if filelist:
            #print 'here',label
            dicts=get_logs(filelist,select_document=load_only)
        elif dictionary:
            #print 'there',label
            dicts=[dictionary]
        #print 'dicts',len(dicts)
        #initialize the logfile with the first document
        self._initialize_class(dicts[0])
        if len(dicts)>1:
            #first initialize the instances with the previous logfile such as to provide the
            #correct information (we should however decide what to do if some run did not converged)
            self._instances=[]
            for i,d in enumerate(dicts):
                dtmp=dicts[0]
                instance=Logfile(dictionary=dtmp,label='log'+str(i))
                #now update the instance with the other value
                instance._initialize_class(d)
                self._instances.append(instance)
            #then we should find the best values for the dictionary
            print 'Found',len(self._instances),'different runs'
            import numpy
            #initalize the class with the dictionary corresponding to the lower value of the energy
            ens=[(l.energy if hasattr(l,'energy') else 1.e100) for l in self._instances] 
            self.reference_log=numpy.argmin(ens)
            #print 'Energies',ens
            self._initialize_class(dicts[self.reference_log])
    def __getitem__(self,index):
        if hasattr(self,'_instances'):
            return self._instances[index]
        else:
            print 'index not available'
            raise 
    def __str__(self):
        return self._print_information()
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
        if hasattr(self,'kpts'):
            self.nkpt=len(self.kpts)
            if hasattr(self,'evals'): self.evals=self._get_bz(self.evals,self.kpts)
            if hasattr(self,'forces') and hasattr(self,'astruct'): 
                self.astruct.update({'forces': self.forces})
                delattr(self,'forces')
        elif hasattr(self,'evals'):
            import BZ
            self.evals=[BZ.BandArray(self.evals),]
        if not hasattr(self,'fermi_level') and hasattr(self,'evals'):
            import numpy
            self.fermi_level=float(max(numpy.ravel(self.evals)))
    def _get_bz(self,ev,kpts):
        evals=[]
        import BZ
        for i,kp in enumerate(kpts):
            evals.append(BZ.BandArray(ev,ikpt=i+1,kpt=kp['Rc'],kwgt=kp['Wgt']))
        return evals
    def get_dos(self,label=None,npts=2500):
        "Get the density of states from the logfile"
        import DoS
        reload(DoS)
        lbl=self.label if label is None else label
        return DoS.DoS(bandarrays=self.evals,label=lbl,units='AU',fermi_level=self.fermi_level,npts=npts)
    def get_brillouin_zone(self):
        "Returns an instance of the BrillouinZone class, useful for band strucure"
        import BZ
        if self.nkpt==1: 
            print 'WARNING: Brillouin Zone plot cannot be defined properly with only one k-point'
            #raise
        mesh=self.kpt_mesh
	if isinstance(mesh,int): mesh=[mesh,]*3
        if self.astruct['Cell'][1]==float('inf'): mesh[1]=1
        return BZ.BrillouinZone(self.astruct,mesh,self.evals,self.fermi_level)
    def geopt_plot(self):
        import numpy
        #for a set of logfiles construct the convergence plot if available
        energies=[]
        forces=[]
        ferr=[]
        if not hasattr(self,'_instances'): 
            print 'ERROR: No geopt plot possible, single point run'
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
            print 'No plot necessary, less than two points found'

    def _print_information(self):
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
