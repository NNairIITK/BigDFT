#!/usr/bin/env python



BIGDFT_CFG='BIGDFT_CONFIGURE_FLAGS'
CLEAN=' clean '
CLEANONE=' cleanone '
UNINSTALL=' uninstall '
LIST=' list '
BUILD=' build '
TINDERBOX=' tinderbox -o build '
DOT=' dot | dot -Tpng > buildprocedure.png '

CHECKMODULES= ['flib','bigdft']

#allowed actions and corresponfing description
ACTIONS={'build':
         'compiles the code with the given configuration',
         'clean':
         'clean the branches for a fresh reinstall',
         'autogen':
         'perform the autogen in the modules which need that. For developers only',
         'dist':
         'creates a tarfile for the bigdft-suite tailored to reproduce the compilation options specified',
         'check':
         'perform check in the bigdft branches, skip external libraries',
         'dry_run':
         'visualize the list of modules that will be compiled with the provided configuration in the buildprocedure.png file'}


class BigDFTInstaller():
    def __init__(self,action,rcfile,verbose):
        import os
        self.action=action
        self.verbose=verbose
        #look where we are
        self.srcdir = os.path.dirname(__file__)
        #hostname
        hostname=os.uname()[1]
        #determine the rcfile
        if rcfile is not None:
            self.rcfile=rcfile
        else:
            rcfile='jhbuildrc'
            if not os.path.exists(rcfile): self.rcfile=''
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

    def get_output(self,cmd):
        import subprocess
        self.__dump('executing:',cmd)
        proc=subprocess.Popen(cmd,stdout=subprocess.PIPE,shell=True)
        (out, err) = proc.communicate()
        self.__dump("program output:", out)
        return out.rstrip('\n')

    def removefile(self,pattern,dirname,names):
        import fnmatch
        "Return the files given by the pattern"
        for name in names:
            if fnmatch.fnmatch(name,pattern):
                self.__dump('removing',os.path.join(dirname,name))
                os.remove(os.path.join(dirname,name))

    def autogen(self):
        import os
        for mod in self.modulelist:
            directory=os.path.join(self.srcdir,mod)
            here = os.getcwd()
            if os.path.isdir(directory):
                self.__dump('Treating directory '+directory)
                os.chdir(directory)
                self.get_output('autoreconf -fi')
                os.chdir(here)
                self.__dump('done.')

    def check(self):
        import os
        for module in [mod for mod in CHECKMODULES if mod in self.modulelist]:
            here = os.getcwd()
            directory=os.path.join(here,module)
            if os.path.isdir(directory):
                self.__dump('Treating directory '+directory)
                os.chdir(directory)
                os.system('make check')
                os.chdir(here)
                self.__dump('done.')
            else:
                print 'Cannot check module "',module,'", directory not present in the build'

    def dist(self):
        print 'dist to be defined'
                        
    def build(self):
        "Build the bigdft module with the options provided by the rcfile"
        #also the force-checkout might be interesting
        import os
        if (self.verbose): 
            os.system(self.jhb+BUILD)
        else:
            os.system(self.jhb+TINDERBOX)
                
    def clean(self):#clean files
        print 'clean'
        for mod in self.modulelist:
            self.get_output(self.jhb+UNINSTALL+mod)
            self.get_output(self.jhb+CLEANONE+mod)
        self.get_output(self.jhb+CLEAN)
        #here we should eliminate residual .mod files
        os.path.walk(".",self.removefile,"*.mod")
        os.path.walk(".",self.removefile,"*.MOD")

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
