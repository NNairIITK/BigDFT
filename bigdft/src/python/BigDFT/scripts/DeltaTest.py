cpickle_file = "bigdft-results.cpickle"

def process_element(run,var,dico,force_run,Elements,strains):
    import math,cPickle
    name= dico["name"]
    hgrids=var['dft']['hgrids']
    ngrids = (math.ceil(strains[0] * dico["a"] / hgrids[0]),
              math.ceil(strains[0] * dico["b"] / hgrids[1]),
              math.ceil(strains[0] * dico["c"] / hgrids[2]))
    for strain in strains:
        if run.iproc == 0:
            print "(%s)" % strain
        if dico["eKS"].has_key(strain) and not force_run:
            #The calculation is already done
            if run.iproc == 0:
                print (name,strain),": ", dico["eKS"][strain]
            continue
        var.update(set_strain(strain,ngrids,dico))

        run.update(var)

        #Run and store the results
        out = run.run()

        var.update(set_restart())

        run.update(var)
        #And restart
        out = run.run()

        Elements[name] = dico
        dico["eKS"][strain] = out.eKS
        dico["FKS"][strain] = out.etot
        dico["pressure"][strain] = out.pressure
        if run.iproc == 0: cPickle.dump(Elements,open(cpickle_file,"w"))
        # Parse the generated log to get the number of grid points.
        for l in open("log-%s.yaml" % var["radical"], "r").xreadlines():
            if "Grid Spacing Units" in l:
                grid = l.split()
        ngrids = (int(grid[5][:-1]), int(grid[6][:-1]), int(grid[7]))
        ngrids = [a + 1 for a in ngrids]
    #Freeing memory
    run = None
    out = None       
    #Build the file
    if run.iproc == 0:
        fd = open("%s.dat" % name,'w')
        for strain in strains:
            dico = Elements[name]
            #Volume in A^3/atom
            volume = dico["volume"]*strain**3/dico["nat"]
            HatoeV = 27.21138386
            eKS = float(dico["eKS"][strain])*HatoeV/dico["nat"]
            fd.write("%16.9f %16.9f\n" % (volume,eKS))
            print strain,dico["eKS"][strain],volume,eKS
        fd.close()


def set_inputfile(hgrid,dico):
    basicinput="""
logfile: Yes
dft:
    ixc: PBE
    ncong: 2
    rmult: [10, 8]
    itermax: 3
    idsx: 0
    gnrm_cv: 1e-8
#Control the diagonalisation scheme
mix:
    iscf: 7
    itrpmax: 200
    rpnrm_cv: 1e-12
    tel: 1e-3
    alphamix: 0.5
    norbsempty: 1000
    alphadiis: 1.0

#perf:
#    accel: OCLGPU
#    ocl_devices: Tesla K40c
#    blas: Yes

"""
    import yaml,os
    var=yaml.load(basicinput)
    #Spin parameters
    var["dft"].update(set_spin(dico["name"], dico["nat"]))
    #K point parameters
    var["kpt"]=set_kpoints(dico["nat"])
    var["dft"]["hgrids"] = (hgrid, hgrid, hgrid)
    #Control the diagonalisation scheme
    if dico["name"] in ("Cr", ):
        var["mix"]["iscf"] = 3
	var["mix"]["alphamix"] = 0.9
    if dico["name"] in ("Ba", "Ca"):
        var["mix"]["norbsempty"] = 8
    var["ig_occupation"] = {dico["name"]: {"empty_shells": ("s", "p", "d")}}
    #Atoms
    pspfile="psppar."+dico["name"]
    if not os.path.isfile(pspfile):
        print "WARNING: Using default PSP for atom",dico["name"]
    else:
        var[pspfile]=open(pspfile,'r').read()
    #var["posinp"] = {"positions": [{dico["name"]: map(float, dico[i + 1].split())} for i in range(dico["nat"])], "units": "reduced", "cell": (dico["a"], dico["b"], dico["c"])}
    var["posinp"] = {"positions": [{dico["name"]: dico[i + 1]} for i in range(dico["nat"])], "units": "reduced", "cell": (dico["a"], dico["b"], dico["c"])}
    # We round robin the igspins.
    if "mpol" in var["dft"]:
        mpol = 0
        while mpol < var["dft"]["mpol"]:
            for at in var["posinp"]["positions"]:
                if mpol < var["dft"]["mpol"]:
                    if "IGSpin" in at:
                        at["IGSpin"] += 1
                    else:
                        at["IGSpin"]  = 1
                    mpol += 1
    elif "nspin" in var["dft"] and var["dft"]["nspin"] == 2:
        for (i, at) in enumerate(var["posinp"]["positions"]):
            at["IGSpin"] = 1 - 2 * (i % 2)
    return var

def set_spin(name,nat):
    "Define the spin in function of the nature of the atoms"
    dspin={}
    if name == 'O':
        dspin["nspin"] = 2
    elif name == 'Cr' or name=='Mn':
        dspin["nspin"] = 2
    elif name == 'Fe' or name == 'Co' or name == 'Ni':
        dspin["nspin"] = 2
        mpol = {"Fe": 2.22, "Co": 1.72, "Ni": 0.60}
        dspin["mpol"] = int(mpol[name] * nat)
        if dspin["mpol"] % 2 == 1:
            dspin["mpol"] += 1
    else:
        dspin["nspin"] = 1
    return dspin

def set_kpoints(nat):
    "Define the k point mesh"
    dkpt={}
    dkpt["method"] = "mpgrid"
    #This line is useless
    dkpt["shiftk"] = ((0., 0., 0.), )
    if nat == 1:
        dkpt["ngkpt"] = (19, 19, 19)
    elif nat == 2:
        dkpt["ngkpt"] = (15, 15, 15)
    elif nat == 3:
        dkpt["ngkpt"] = (14, 14, 14)
    elif nat == 4:
        dkpt["ngkpt"] = (12, 12, 12)
    elif nat == 6:
        dkpt["ngkpt"] = (11, 11, 11)
    elif nat == 8:
        dkpt["ngkpt"] = (10, 10, 10)
    elif nat == 12:
        dkpt["ngkpt"] = (9, 9, 9)
    return dkpt

def set_strain(strain,ngrids,dico):
    return {"radical": dico["name"] + "-%s" % strain,
            "posinp": {"cell": (strain * dico["a"],
                                strain * dico["b"],
                                strain * dico["c"])},
            "dft": {"hgrids": ((strain * dico["a"] + 1e-4) / ngrids[0],
                               (strain * dico["b"] + 1e-4) / ngrids[1],
                               (strain * dico["c"] + 1e-4) / ngrids[2])}}

def set_restart():
    return {"mix": {"tel": 1.e-5}, "dft": {"inputpsiid": 1}}

#Build a dictionary
def elements_from_cif(files):
    #filename match as Unix shell
    import fnmatch,os,math
    Elements = dict()
    nonortho = list()
    ortho = list()
    for file in files:
        if not fnmatch.fnmatch(file,"*.cif"): continue
        dico = dict()
        dico["file"] = file
        for line in open(os.path.abspath(file)).readlines(): #open("%s/%s" %(dirCIF,file)).readlines():
            if "_" in line:
                #Remove quote '
                items = line.replace("'","").split()
                if len(items) == 2:
                    dico[items[0]] = items[1]
            elif "Biso" in line:
                items = line.split()
                #We know that the maximum of atoms is 8 in the cell
                n = int(items[0][-1])
                dico['nat'] = n
                #dico[n] = "%s  %s  %s " % (items[2], items[3], items[4])
                dico[n] = map(float, items[2:5])
        #Convert name
        dico["name"] = dico["_pd_phase_name"]
        #We use bohrs
        atob = 1.0/0.5291772108
        dico["a"] = float(dico["_cell_length_a"])*atob
        dico["b"] = float(dico["_cell_length_b"])*atob
        dico["c"] = float(dico["_cell_length_c"])*atob
        dico["alpha"] = dico["_cell_angle_alpha"]
        dico["beta"] = dico["_cell_angle_beta"]
        dico["gamma"] = dico["_cell_angle_gamma"]
        #Only for non-orthorhombic and in angstroem^3
        dico["volume"] = dico["a"]*dico["b"]*dico["c"]/atob**3
        #Create a key results
        dico["eKS"] = dict()
        dico["FKS"] = dict()
        dico["pressure"] = dict()
        btype = None
        # Look for possible orthorombic transformation :
        if (dico["alpha"] != "90" and dico["b"] == dico["c"] and dico["beta"] == "90" and dico["gamma"] == "90"):
          btype = "hexagonal"
          la = "b"
          lb = "c"
          ang = "alpha"
          P = ((1.,  0.,  0.),
               (0.,  0.5, 0.5),
               (0., -0.5, 0.5))
          du = 0.
          dv = 0.
          dw = 1.
        elif (dico["beta"] != "90" and dico["a"] == dico["c"] and dico["alpha"] == "90" and dico["gamma"] == "90"):
          la = "a"
          lb = "c"
          btype = "hexagonal"
          ang = "beta"
          P = (( 0.5, 0., 0.5),
               ( 0.,  1., 0. ),
               (-0.5, 0., 0.5))
          du = 0.
          dv = 0.
          dw = 1.
        elif (dico["gamma"] != "90" and dico["a"] == dico["b"] and dico["alpha"] == "90" and dico["beta"] == "90"):
          la = "a"
          lb = "b"
          btype = "hexagonal"
          ang = "gamma"
          P = (( 0.5, 0.5, 0.),
               (-0.5, 0.5, 0.),
               ( 0.,  0.,  1.))
          du = 0.
          dv = 1.
          dw = 0.
        elif (dico["gamma"] == "90" and dico["alpha"] == "90" and dico["beta"] == "90"):
          btype = "orthorombic"
        elif dico["a"] == dico["b"] and dico["a"] == dico["c"] and dico["alpha"] == dico["beta"]and dico["alpha"] == dico["gamma"]:
          btype = "rhombohedral"
    
        # Transform to orthorombic when possible.
        if btype == "hexagonal":
          a = dico[la]
          b = dico[lb]
          # - hexagonal case.
          alpha = float(dico[ang]) * math.pi / 180.
          dico[la] = math.sqrt((a+b*math.cos(alpha))**2 + (b*math.sin(alpha))**2)
          dico[lb] = math.sqrt((a-b*math.cos(alpha))**2 + (b*math.sin(alpha))**2)
          for i in range(dico["nat"]):
            u, v, w = dico[i+1]
            a = P[0][0] * u + P[0][1] * v + P[0][2] * w
            b = P[1][0] * u + P[1][1] * v + P[1][2] * w
            c = P[2][0] * u + P[2][1] * v + P[2][2] * w
            dico[i + 1] = (a, b, c)
            u += du
            v += dv
            w += dw
            a = P[0][0] * u + P[0][1] * v + P[0][2] * w
            b = P[1][0] * u + P[1][1] * v + P[1][2] * w
            c = P[2][0] * u + P[2][1] * v + P[2][2] * w
            dico[dico["nat"] + i + 1] = (a, b, c)
          dico["nat"] *= 2
          dico[ang] = "90"
        elif btype == "rhombohedral" and float(dico["alpha"]) == 60.:
          a = dico["a"]
          dico["a"] = a * math.sqrt(2.)
          dico["b"] = a * math.sqrt(2.)
          dico["c"] = a * math.sqrt(2.)
          P = ((0.,  0.5, 0.5),
               (0.5, 0.,  0.5),
               (0.5, 0.5, 0. ))
          dd = ((0., 0., 0.),
                (1., 0., 0.),
                (0., 1., 0.),
                (0., 0., 1.))
          for i in range(dico["nat"]):
            for (j, (du, dv, dw)) in enumerate(dd):
              u, v, w = dico[i+1]
              u += du
              v += dv
              w += dw
              a = P[0][0] * u + P[0][1] * v + P[0][2] * w
              b = P[1][0] * u + P[1][1] * v + P[1][2] * w
              c = P[2][0] * u + P[2][1] * v + P[2][2] * w
              dico[j * dico["nat"] + i + 1] = (a, b, c)
    	      dico["nat"] *= len(dd)
    	      dico["alpha"] = "90"
    	      dico["beta"] = "90"
    	      dico["gamma"] = "90"
        elif btype != "orthorombic":
            print "to be treated", dico["name"], dico["alpha"], dico["beta"], dico["gamma"], dico["a"], dico["b"], dico["c"]
    	    # Update volume after orthomrombic tranformation and in angstroem^3
    	    dico["volume"] = dico["a"]*dico["b"]*dico["c"]/atob**3
        name = dico['name']
        if dico["alpha"] != "90" or dico["beta"] != "90" or dico["gamma"] != "90":
            nonortho.append(name)
            continue
        else:
            ortho.append(name)
        #Add in the large dictionary
        Elements[dico["_pd_phase_name"]] = dico
    ortho.sort()
    nonortho.sort()
    return Elements,ortho,nonortho

def xyz_from_elements(dirXYZ,Elements,ortho,nonortho):
    import shutil,os,time
    format_xyz = """{0[nat]} reduced
 periodic {0[a]}   {0[b]}   {0[c]} 
"""
    #print "Remove the directory '%s'." % dirXYZ
    shutil.rmtree(dirXYZ,ignore_errors=True)
    #Create it
    os.mkdir(dirXYZ)
    print "---"
    #Start and create the xyz files
    print "Delta-test timestamp:", time.strftime('%X %x %Z')
    print "Delta-test code: BigDFT"
    print "Number of elements: ", len(Elements)
    print "List of elements:", Elements.keys()
    print "Number of orthorhombic elements: ", len(ortho)
    print "Orthorhombic elements: ", ortho
    print "Number of non-orthorhombic elements: ", len(nonortho)
    print "Non-orthorhombic elements: ", nonortho
    for dico in Elements.values():
        name = dico['name']
        #We have all the specification
        fnew = os.path.join(dirXYZ,name+".xyz") #"%s/%s.xyz" % (dirXYZ,name)
        fd = open(fnew,"w")
        fd.write(format_xyz.format(dico))
        for i in range(dico['nat']):
            fd.write("%s " % name)
            fd.write("%f %f %f\n" % tuple(dico[i+1]))
            #fd.write("%s %s\n" % (name,dico[i+1]))
        fd.close()
        #print "#Creation of the file '{0:s}' from '{1:s}/{2:s}'".format(fnew,dirCIF,dico['file'])


class Benchmark():
    def __init__(self,calculator):
        import os
        #import the positions of the periodic table from the CIFs directory
        dirCIF = os.path.join(os.path.dirname(__file__),'CIFs')
        files = [os.path.abspath(os.path.join(dirCIF,f)) for f in os.listdir(dirCIF)]
        self.Elements,self.ortho,self.nonortho=elements_from_cif(files)
        #Strains (The original)
        self.strains = [ 0.94, 0.96, 0.98, 1.00, 1.02, 1.04, 1.06 ]
        #From Francois Jollet
        #strains = [ 0.9795861087155615, 0.986484829732188, 0.9932883883792687,  1.00, 1.006622709560113, 1.0131594038201772, 1.0196128224222163 ]
        #ortho = [ "Ag"]
        self.strains.reverse()
        #strains = [ 1.00 ]
        self.calculator=calculator
        #Now build "*.xyz" for each elements
        if self.calculator.iproc == 0 and False: 
            xyz_from_elements("XYZs",self.Elements,self.ortho,self.nonortho)
        
    def run(self,atoms=None,force=False):
        import os
        #To dump python object in a file
        import cPickle
        if atoms is None: 
            atomlist=self.ortho
        else:
            atomlist=atoms
        #Loop over the elements 
        #Test if the file bigdft-results.cpickle exists and use it instead
        if os.path.exists(cpickle_file) and not force:
            self.Elements = cPickle.load(open(cpickle_file,"r"))
        for name in atomlist:
            dico=self.Elements[name]
            if self.calculator.iproc == 0: print "Start calculation of %s" % dico["name"]
            hgrid = 0.30
            var = set_inputfile(hgrid,dico)
            self.calculator.set(var)
            process_element(self.calculator,var,dico,force,self.Elements,self.strains)
