#!/usr/bin/env python



BIGDFT_CFG='BIGDFT_CONFIGURE_FLAGS'
CLEAN=' clean '
CLEANONE=' cleanone '
UNINSTALL=' uninstall '
LIST=' list '
BUILD=' build '
TINDERBOX=' tinderbox -o build '
DOT=' dot | dot -Tpng > buildprocedure.png '
DIST=' distone bigdft-suite '

CHECKMODULES= ['flib','bigdft']
MAKEMODULES= ['flib','libABINIT','bigdft']

#allowed actions and corresponfing description
ACTIONS={'build':
         'Compiles and install the code with the given configuration',
         'make':
         'Recompiles the bigdft internal branches, skip configuring step',
         'clean':
         'Clean the branches for a fresh reinstall',
         'autogen':
         'Perform the autogen in the modules which need that. For developers only',
         'dist':
         'Creates a tarfile for the bigdft-suite tailored to reproduce the compilation options specified',
         'check':
         'Perform check in the bigdft branches, skip external libraries',
         'dry_run':
         'Visualize the list of modules that will be compiled with the provided configuration in the buildprocedure.png file'}


class BigDFTInstaller():
    def __init__(self,action,rcfile,verbose):
        import os
        self.action=action
        self.verbose=verbose
        #look where we are
        self.srcdir = os.path.dirname(__file__)
        #look the builddir
        self.builddir=os.getcwd()
        if os.path.abspath(self.srcdir) == os.path.abspath(self.builddir):
            print 50*'-'
            print "ERROR: BigDFT Installer works better with a build directory different from the source directory, install from another directory"
            print "SOLUTION: Create a separate directory and invoke this script from it"
            print 50*'-'
            exit(1)
        #hostname
        hostname=os.uname()[1]
        #determine the rcfile
        if rcfile is not None:
            self.rcfile=rcfile
        else:
            self.rcfile='jhbuildrc'
            if not os.path.exists(self.rcfile): self.rcfile=''
            #check for environment variable or for configuration file
            if self.rcfile == '' and (BIGDFT_CFG not in os.environ.keys()):
                #here we might explore the list of present rc files and prompt for their usage
                #print hostname
                #rcfile might be changed
                #rcfile=''
                raise 'ERROR: no rcfile provided and '+BIGDFT_CFG+' variable not present, exiting...'

        #jhbuild script
        self.jhb=os.path.join(self.srcdir,'jhbuild.py ')
        if self.rcfile != '': self.jhb += '-f '+self.rcfile
                
        #now get the list of modules that has to be treated with the given command
        self.modulelist=self.get_output(self.jhb + LIST).split('\n')
        self.__dump("List of modules to be treated",self.modulelist)

        #then choose the actions to be taken
        getattr(self,action)()
                        
    def __dump(self,*msg):
        if self.verbose:
            for m in msg:
                print m
                
    def selected(self,l):
        return [val for val in l if val in self.modulelist]

    def shellaction(self,path,modules,action):
        import os
        for mod in self.selected(modules):
            directory=os.path.join(path,mod)
            here = os.getcwd()
            if os.path.isdir(directory):
                self.__dump('Treating directory '+directory)
                os.chdir(directory)
                os.system(action)
                os.chdir(here)
                self.__dump('done.')
            else:
                print 'Cannot perform action "',action,'" on module "',mod,'" directory not present in the build'
    
    def get_output(self,cmd):
        import subprocess
        self.__dump('executing:',cmd)
        proc=subprocess.Popen(cmd,stdout=subprocess.PIPE,shell=True)
        (out, err) = proc.communicate()
        self.__dump("program output:", out)
        return out.rstrip('\n')

    def removefile(self,pattern,dirname,names):
        import os,fnmatch
        "Return the files given by the pattern"
        for name in names:
            if fnmatch.fnmatch(name,pattern):
                self.__dump('removing',os.path.join(dirname,name))
                os.remove(os.path.join(dirname,name))

    def autogen(self):
        self.shellaction(self.srcdir,self.modulelist,'autoreconf -fi')

    def check(self):
        self.shellaction('.',CHECKMODULES,'make check')

    def make(self):
        self.shellaction('.',MAKEMODULES,'make -j6 && make install')
        
    def dist(self):
        self.shellaction('.',self.modulelist,'make dist')
        self.get_output(self.jhb+DIST)
                                
    def build(self):
        "Build the bigdft module with the options provided by the rcfile"
        import os
        if (self.verbose): 
            os.system(self.jhb+BUILD)
        else:
            os.system(self.jhb+TINDERBOX)

    def clean(self):#clean files
        import os
        for mod in self.selected(MAKEMODULES):
            self.get_output(self.jhb+UNINSTALL+mod)
            self.get_output(self.jhb+CLEANONE+mod)
            #here we should eliminate residual .mod files
            os.path.walk(mod,self.removefile,"*.mod")
            os.path.walk(mod,self.removefile,"*.MOD")
        #self.get_output(self.jhb+CLEAN)

    def dry_run(self):
        self.get_output(self.jhb+DOT)

    def __del__(self):
        print 'Thank you for using the Installer of BigDFT suite.'
        print 'Your used configuration options have been saved in the file jhbuildrc'
        print 'The action taken was:',self.action


#now follows the available actions, argparse might be called
import argparse
parser = argparse.ArgumentParser(description='BigDFT suite Installer',
                                 epilog='For more information, visit www.bigdft.org')
parser.add_argument('-f','--file',
                   help='Use an alternative configuration file instead of the default given by the environment variable BIGDFT_CONFIGURE_FLAGS')
parser.add_argument('-d','--verbose',action='store_true',
                   help='Verbose output')

parser.add_argument('action',choices=[act for act in ACTIONS],
                   help='Define the installer action to be taken')
args = parser.parse_args()

BigDFTInstaller(args.action,args.file,args.verbose)
