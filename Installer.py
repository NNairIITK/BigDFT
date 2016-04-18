#!/usr/bin/env python
# -*- coding: us-ascii -*-

#--------------------------------------------------------------------------------
# Copyright (C) 2015-2016 BigDFT group
# This file is distributed under the terms of the
# GNU General Public License, see ~abinit/COPYING
# or http://www.gnu.org/copyleft/gpl.txt .
#--------------------------------------------------------------------------------

BIGDFT_CFG='BIGDFT_CONFIGURE_FLAGS'
CLEAN=' clean '
CLEANONE=' cleanone '
UNINSTALL=' uninstall '
LIST=' list '
BUILD=' build '
BUILDONE=' buildone '
TINDERBOX=' tinderbox -o build '
DOT=' dot '
DOTCMD=' | dot -Edir=back -Tpng > buildprocedure.png '
DIST='  dist --dist-only bigdft-suite '
RCFILE='buildrc'
SETUP=' setup '

CHECKMODULES= ['futile','psolver','bigdft','spred']
MAKEMODULES= ['futile','psolver','libABINIT','bigdft','spred']

#allowed actions and corresponding description
ACTIONS={'build':
         'Compile and install the code with the given configuration.',
         'make':
         'Recompile the bigdft internal branches, skip configuring step.',
         'clean':
         'Clean the branches for a fresh reinstall.',
         'startover':
         'Wipe out all the build directories and recompile the important parts',
         'autogen':
         'Perform the autogen in the modules which need that. For developers only.',
         'dist':
         'Creates a tarfile for the bigdft-suite tailored to reproduce the compilation options specified.',
         'check':
         'Perform check in the bigdft branches, skip external libraries.',
         'dry_run':
         "Visualize the list of modules that will be compiled with the provided configuration in the 'buildprocedure.png' file."}


class BigDFTInstaller():
    def __init__(self,action,package,rcfile,verbose,quiet):
        import os
        self.action=action    #Action to be performed
        self.package=package  #Package
        self.verbose=verbose  #verbose option
        self.quiet=quiet      #Ask a question
        #look where we are
        self.srcdir = os.path.dirname(__file__)
        #look the builddir
        self.builddir=os.getcwd()
        #look if we are building from a branch
        bigdftdir=os.path.join(self.srcdir,'bigdft')
        self.branch=os.path.isfile(os.path.join(bigdftdir,'branchfile'))

        #To be done BEFORE any exit instruction in __init__ (get_rcfile)
        self.time0 = None

        if os.path.abspath(self.srcdir) == os.path.abspath(self.builddir) and self.action != ['autogen','dry_run']:
            print 50*'-'
            print "ERROR: BigDFT Installer works better with a build directory different from the source directory, install from another directory"
            print "SOLUTION: Create a separate directory and invoke this script from it"
            exit(1)
        #hostname
        self.hostname=os.uname()[1]

        #rcfile
        self.get_rcfile(rcfile)

        #jhbuild script
        self.jhb=os.path.join(self.srcdir,'jhbuild.py ')
        if self.rcfile != '': self.jhb += '-f '+self.rcfile

        #date of bigdft executable if present
        self.time0=self.bigdft_time()

        self.print_present_configuration()

        #now get the list of modules that has to be treated with the given command
        self.modulelist=self.get_output(self.jhb + LIST+self.package).split('\n')
        print " List of modules to be treated:",self.modulelist

        #then choose the actions to be taken
        getattr(self,action)()

    def bigdft_time(self):
        import os
        return self.filename_time(os.path.join(self.builddir,'install','bin','bigdft'))

    def filename_time(self,filename):
        import os
        if os.path.isfile(filename):
            return os.path.getmtime(filename)
        else:
            return 0

    def get_rcfile(self,rcfile):
        "Determine the rcfile"
        import os
        if rcfile is not None:
            self.rcfile=rcfile
        else:
            self.rcfile=RCFILE
        #see if it exists where specified
        if os.path.exists(self.rcfile): return
        #otherwise search again in the rcfiles
        rcdir=os.path.join(self.srcdir,'rcfiles')
        self.rcfile=os.path.join(rcdir,self.rcfile)
        if os.path.exists(self.rcfile): return
        #see if the environment variables BIGDFT_CFG is present
        self.rcfile = ''
        if BIGDFT_CFG in os.environ.keys(): return
        #otherwise search for rcfiles similar to hostname and propose a choice
        rcs=[]
        for file in os.listdir(rcdir):
            testname=os.path.basename(file)
            base=os.path.splitext(testname)[0]
            if base in self.hostname or self.hostname in base: rcs.append(file)
        print "Search in the configuration directory '%s'" % rcdir
        if len(rcs)==1:
            self.rcfile=os.path.join(rcdir,rcs[0])
        elif len(rcs) > 0:
            print "No valid configuration file specified, found various that matches the hostname '%s'" % self.hostname
            print 'In the directory "'+rcdir+'"'
            print 'Choose among the following options'
            for i,rc in enumerate(rcs):
                print str(i+1)+'. '+rc
            while True:
                choice=raw_input('Pick your choice (q to quit) ')
                if choice == 'q': exit(0)
                try:
                    ival=int(choice)
                    if (ival <= 0): raise
                    ch=rcs[ival-1]
                    break
                except:
                    print 'The choice must be a valid integer among the above'
            self.rcfile=os.path.join(rcdir,ch)
        elif len(rcs) == 0:
            print 'No valid configuration file provided and '+BIGDFT_CFG+' variable not present, exiting...'
            exit(1)

    def __dump(self,*msg):
        if self.verbose:
            for m in msg:
                print m

    def print_present_configuration(self):
        import  os
        indent = ' '*2
        print 'Configuration chosen for the Installer:'
        print indent + 'Hostname:',self.hostname
        print indent + 'Source directory:',os.path.abspath(self.srcdir)
        print indent + 'Compiling from a branch:',self.branch
        print indent + 'Build directory:',os.path.abspath(self.builddir)
        print indent + 'Action chosen:',self.action
        print indent + 'Verbose:',self.verbose
        print indent + 'Configuration options:'
        if self.rcfile=='':
            print indent*2 + "Source: Environment variable '%s'" % BIGDFT_CFG
	    print indent*2 + "Value: '%s'" % os.environ[BIGDFT_CFG]
        else:
            print indent*2 + "Source: Configuration file '%s'" % os.path.abspath(self.rcfile)
        while not self.quiet:
            ok = raw_input('Do you want to continue (Y/n)? ')
            if ok == 'n' or ok=='N':
                exit(0)
            elif ok != 'y' and ok != 'Y' and repr(ok) != repr(''):
                print 'Please answer y or n'
            else:
                break

    def selected(self,l):
        return [val for val in l if val in self.modulelist]

    def shellaction(self,path,modules,action,hidden=False):
        "Perform a shell action, dump also the result if verbose is True."
        import os
        import sys
        for mod in self.selected(modules):
            directory=os.path.join(path,mod)
            here = os.getcwd()
            if os.path.isdir(directory):
                #self.__dump('Treating directory '+directory)
                sys.stdout.write('Module '+mod+' ['+directory+']: '+action)
                sys.stdout.flush()
                os.chdir(directory)
                if hidden:
                    self.get_output(action)
                else:
                    os.system(action)
                os.chdir(here)
                #self.__dump('done.')
                sys.stdout.write(' (done)\n')
            else:
                sys.stdout.write('Cannot perform action "'+action+'" on module "'+mod+'" directory not present in the build.\n')
            sys.stdout.flush()

    def get_output(self,cmd):
        import subprocess
        self.__dump('executing:',cmd)
        proc=subprocess.Popen(cmd,stdout=subprocess.PIPE,shell=True)
        (out, err) = proc.communicate()
        self.__dump("program output:", out)
        return out.rstrip('\n')

    def removefile(self,pattern,dirname,names):
        "Return the files given by the pattern"
        import os,fnmatch
        for name in names:
            if fnmatch.fnmatch(name,pattern):
                self.__dump('removing',os.path.join(dirname,name))
                os.remove(os.path.join(dirname,name))

    def autogen(self):
        "Perform the autogen action"
        self.shellaction(self.srcdir,self.modulelist,'autoreconf -fi')

    def check(self):
        "Perform the check action"
        self.shellaction('.',CHECKMODULES,'make check',hidden=not self.verbose)

    def make(self):
        "Perform the simple make action"
        self.shellaction('.',MAKEMODULES,'make -j6 && make install',hidden=not self.verbose)

    def dist(self):
        "Perform make dist action"
        import os
        tarfile=os.path.join(self.builddir,'bigdft-suite.tar.gz')
        disttime0=self.filename_time(tarfile)
        os.system(self.jhb+DIST)
        disttime1=self.filename_time(tarfile)
        if not (disttime1 == disttime0):
            print 'SUCCESS: distribution file "bigdft-suite.tar.gz" generated correctly'
        else:
            print 'WARNING: the dist file seems not have been updated or generated correctly'

    def build(self):
        "Build the bigdft module with the options provided by the rcfile"
        import os
        #in the case of a nonbranch case, like a dist build, force checkout
        #should the make would not work
        if self.branch:
            co=''
        else:
            co=' -C'
        if (self.verbose):
            os.system(self.jhb+BUILD+self.package+co)
        else:
            os.system(self.jhb+TINDERBOX+self.package+co)

    def clean(self):
        "Clean files in the build directory"
        import os
        for mod in self.selected(MAKEMODULES):
            self.get_output(self.jhb+UNINSTALL+mod)
            self.get_output(self.jhb+CLEANONE+mod)
            #here we should eliminate residual .mod files
            os.path.walk(mod,self.removefile,"*.mod")
            os.path.walk(mod,self.removefile,"*.MOD")
        #self.get_output(self.jhb+CLEAN)

    def startover(self):
        "Wipe files in the makemodules directory"
        if not self.branch:
            print 'ERROR: The action "startover" is allowed only from a developer branch'
            exit(1)
        import shutil
        import os
        for mod in self.selected(MAKEMODULES):
            self.get_output(self.jhb+UNINSTALL+mod)
            print 'Wipe directory: ',mod
            shutil.rmtree(mod, ignore_errors=True)
        print 'Building again...'
        startat=' -t '
        for mod in self.selected(MAKEMODULES):
            print 'Resetting: ',mod
            self.get_output(self.jhb+SETUP+mod+startat+mod)
            print 'Building: ',mod
            self.get_output(self.jhb+BUILDONE+mod)
        self.build()
        
    def dry_run(self):
        "Do dry build"
        self.get_output(self.jhb+DOT+self.package+DOTCMD)

    def rcfile_from_env(self):
        "Build the rcfile information from the chosen "+BIGDFT_CFG+" environment variable"
        import os
        if os.path.isfile(self.rcfile) and not os.path.isfile(RCFILE):
            from shutil import copyfile
            copyfile(self.rcfile,RCFILE)
            return
        if BIGDFT_CFG not in os.environ.keys() or os.path.isfile(RCFILE): return
        print 'The suite has been built without configuration file.'
        rclist=[]
        rclist.append("""#This is the configuration file for the BigDFT installer""")
        rclist.append("""#This is a python script which is executed by the build suite """)
        rclist.append(" ")
        rclist.append("""#Add the condition testing to run tests and includes PyYaml""")
        rclist.append("""conditions.add("testing")""")
        rclist.append("""#List the module the this rcfile will build""")
        rclist.append("modules = ['"+self.package+"',]")
        sep='"""'
        confline=sep+os.environ[BIGDFT_CFG]+sep
        rclist.append("#example of the potentialities of the python syntax in this file")
        rclist.append("def env_configuration():")
        rclist.append("    return "+confline)
        rclist.append("""#here follow the configuration instructions for the modules built""")
        rclist.append("module_autogenargs.update({")
        rclist.append("   ")
        for mod in self.modulelist:
            rclist.append("'"+mod+"': env_configuration(),")
            rclist.append("   ")
        rclist.append("})")
        #then write the file
        rcfile=open(RCFILE,'w')
        for item in rclist:
            rcfile.write("%s\n" % item)
            #rcfile.write("\n")
        rcfile.close()
        print "Your used configuration options have been saved in the file '%s'." % RCFILE
        print "Such file will be used for next builds, you might also save it in the 'rcfiles/'."
        print "Directory of the source for future use. The name might contain the hostname."

    def __del__(self):
        import os
        print 50*'-'
        print 'Thank you for using the Installer of BigDFT suite.'
        print 'The action considered was:',self.action
        try:
           if self.time0 is not None:
               if self.action in ['build','dry_run']: self.rcfile_from_env()
               if not (self.time0==self.bigdft_time()):
                   print 'SUCCESS: The Installer seems to have built correctly bigdft bundle'
                   print 'All the available executables and scripts can be found in the directory'
                   print '"'+os.path.join(os.path.abspath(self.builddir),'install','bin')+'"'
               elif (self.action == 'build' or self.action == 'make'):
                   print 'WARNING: The Installer seems NOT have created or updated bigdft executable'
                   print '        (maybe everything was already compiled?)'
                   print 'ACTION: check the compiling procedure.'
                   if self.branch:
                       print 'HINT: It appears you are compiling from a branch source tree. Did you perform the action "autogen"?'
                   if not self.verbose and self.action == 'build':
                      print '  HINT: Have a look at the file index.html of the build/ directory to find the reason'
        except:
            print 'Goodbye...'

#Now follows the available actions, argparse might be called
import argparse

#Redefine ArgumentParser to have the help message if no arguments
class Installer_Parser(argparse.ArgumentParser):
    def error(self, message):
        import sys
        sys.stderr.write('error: %s\n' % message)
        self.print_help()
        self.exit()

parser = Installer_Parser(description='BigDFT suite Installer',
                            epilog='''
If you want more help type "%(prog)s help"
------------------------------------------------
For more information, visit www.bigdft.org''',
                            formatter_class=argparse.RawDescriptionHelpFormatter)

parser.add_argument('action',nargs='?',default='help',
                    help='Action to be performed by the Installer.'
                    ' (default: %(default)s)',choices=['help']+[a for a in ACTIONS])
parser.add_argument('package',nargs='?',default='spred',
                    help='Package to be built by the installer. (default: %(default)s)',
                    choices=CHECKMODULES)
parser.add_argument('-f','--file',
                   help='Use an alternative configuration file instead of the default configuration '
                    + 'given by the environment variable %s' % BIGDFT_CFG)
parser.add_argument('-d','--verbose',action='store_true',
                   help='Verbose output')
parser.add_argument('-q','--quiet',action='store_true',
                   help='Skip dialog after setup')
parser.add_argument('-c','--configure-line',nargs=argparse.REMAINDER,
                   help='Specify the configure line to be passed (set BIGDFT_CONFIGURE_FLAGS variable)')



###Define the possible actions
##subparsers = parser.add_subparsers(title='The following actions are available',
##                    dest='action',
##                    help='Action to be performed by the Installer.')
##for (k,v) in ACTIONS.items():
##    subparsers.add_parser(k,help=v)
##
args = parser.parse_args()


if args.configure_line is not None:
  cfg=''
  for i in args.configure_line:
      cfg+='"'+i+'" '
  #scratch the BIGDFT_CFG environment variable
  import os
  os.environ[BIGDFT_CFG]=cfg

if args.action=='help':
    print "Quick overview of the BigDFT suite Installer program"
    print 50*'-'
    print "USAGE: Installer.py <action> <package>"
    print 50*'-'+'Available actions'
    for a in ACTIONS:
        print a,':'
        print '     ',ACTIONS[a]
    print 50*'-'
    print 'Available packages:',CHECKMODULES
    print 50*'-'
    print 10*"QIFI-"+' (Quick Instructions For the Impatient)'
    print 'Ideally, there are two different policies:'
    print 'Developer: From a development branch, start by "startover", then "build"'
    print '     User: From a tarball, start by "build"'
    print 'Perform the "dry_run" command to have a graphical overview of the building procedure'
else:
    BigDFTInstaller(args.action,args.package,args.file,args.verbose,args.quiet)
