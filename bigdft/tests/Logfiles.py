import yaml

EVAL = "eval"
SETUP = "let"
INITIALIZATION = "globals"

PRE_POST = [EVAL, SETUP, INITIALIZATION]

ENERGY = "BigDFT.energy"
FERMI_LEVEL= "__FERMI_LEVEL__"
NUMBER_OF_ATOMS = 'BigDFT.nat'
EIGENVALUES = 'BigDFT.evals'

#Builtin pathes to define the search paths
BUILTIN={ENERGY: [["Last Iteration", "FKS"],["Last Iteration", "EKS"], ["Energy (Hartree)"]],
         FERMI_LEVEL: [["Ground State Optimization", -1, "Fermi Energy"]],
         NUMBER_OF_ATOMS: [ ['Atomic System Properties','Number of atoms']],
         EIGENVALUES: [ ["Complete list of energy eigenvalues"], [ "Ground State Optimization", -1, "Orbitals"],
                        ["Ground State Optimization",-1,"Hamiltonian Optimization",-1,"Subspace Optimization","Orbitals"] ]}

def get_log(f):
    "Transform a logfile into a python dictionary"
    return yaml.load(open(f, "r").read(), Loader = yaml.CLoader)

def get_logs(files):
   logs=[]
   for filename in files:
     try:
        logs+=[yaml.load(open(filename, "r").read(), Loader = yaml.CLoader)]
     except:
        try: 
            logs+=yaml.load_all(open(filename, "r").read(), Loader = yaml.CLoader)
        except:
            logs+=[None]
            print "warning, skipping logfile",filename
   return logs

# this is a tentative function written to extract information from the runs
def document_quantities(doc,to_extract):
  analysis={}
  for quantity in to_extract:
    if quantity in PRE_POST: continue
    #follow the levels indicated to find the quantity
    field=to_extract[quantity]
    if type(field) is not type([]) is not type({}) and field in BUILTIN:
        paths=BUILTIN[field]
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
