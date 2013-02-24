#!/usr/bin/env python
# -*- coding: us-ascii -*-
#--------------------------------------------------------------------------------
# Copyright (C) 2008-2012 BigDFT group (TD)
# This file is distributed under the terms of the
# GNU General Public License, see ~abinit/COPYING
# or http://www.gnu.org/copyleft/gpl.txt .
# For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
#--------------------------------------------------------------------------------
# VERSION for BigDFT
#
# Try to have a common definition of classes with abilint (ABINIT)
#
# Date: 04/12/2012
#--------------------------------------------------------------------------------
#i# Lines commented: before used for #ifdef interfaces

#Principles
#----------
# 1 - Analyze without changing any part of the code
# 2 - Finally change if it is needed


#Plan
#----
# 1 - List all routines and store all in memory.
# 2 - Search all routines belonging to a directory.
# 3 - Remove all interfaces contained in the source files.
# 4 - Create a module containing the interface of all routines
#     for each directory.
# 5 - Add this file in the directory 'src/defs' (this is a choice which is not the best I think)
# 6 - Beautify the code if required
# 7 - Build some graphs concerning the calls of routines
# 8 - Copy the new version


#TODO
#----
# Add include files for dependencies
# Pb with pre-processing (?)
# The module ftfvw2 needs interfaces: Check if everything is correct
# Add correct children to _xxx_ files.
# Add children and parents to the files interfaces_xxx.F90.
# Add a true functionality to remove variables
# Call after a ';' is not detected.
# Handle the preprocessing directives (?)
# Unify the detection and the creation of the interfaces by means of the class Fortran_Interface
# Improve the experimental edition


#Description
#-----------
# The script abilint uses class to implement the notion of project, files, routines and so on.
# The class Project (maybe Package is more appropriate) 
# contains all the information about directories, files and routines.
# There are lists of:
# - directories (abilint.dirs) which are only names of directories,
# - files (abilint.files) which are the File objects (see File class) and
# - routines (abilint.routines) which are the Routine objects (see Routine class).
# So it is easy to access directly to an object in the project.
# The only restriction is that the names for the routines and the files have to be unique.
# The class Structure is a general class which implements the method of read or write which is
# propagated by inheritance for File, File_F90, File_F77, Module, Routine, Function, Program, Code,
# Comment, Declaration, Execution, Fortran_Type, ...

#Hierarchy of classes
#--------------------
# File_F77 is a special class of the class File_F90 with a filter to parse the code Fortran77.

# File_F90 -> Comment (Robodoc_Header, Robodoc_Footer)
#          -> Module   ----> Comment, Declaration, Routine, Function
#          -> Routine  --|
#          -> Function   |
#                        --> Comment, Header_Routine, Use, Implicit, Declaration, Execution
#                                                                    |
#                                                                    |
#                                Fortran_Type, Fortran_Interface <----
#
# The class Code is a generic class for all program structures.
#
# Each part of a code is handled by a class:
# Comment   lines of comments, specially robodoc
# Module    module fortran
# Routine   subroutine fortran
# Function  function fortran

# The classes 'Module', 'Routine', 'Function' 
#              contain methods to handle the code and specially to analyze it into:
# Comment      lines of comments
# Header       header of the routine (ex. subroutine foo(one,two))
# Use          use statements
# Implicit     implicit statements
# Declaration  declaration statements of all variables
# Execution    execution statements
#
# The class Include is a special class to handle the include statement.
# abilint gives an ersatz of 'mpif.f' file for commodity.
#
# Two special classes are not added as children because they duplicate the code:
# Fortran_Type      inherited of the class Declaration to handle fortran types.
# Fortran_Interface inherited of the class Fortran_Type to handle declaration of interfaces.
#
# The class Include is a special class to handle the include statement.
# abilint gives an ersatz of 'mpif.f' file for commodity.

# The class Variable is created to have all information for a given variable.
# A structure has a 'dict_vars' dictionary.
#
#Edition of the files
#--------------------
# Vim is used to help the edition


#Usage (and help)
def usage(error=0):
    "Display the usage of abilint"
    sys.stderr.write("""
abilint [options] <source> <dest>
    where:
      <source> is the directory that contains the bigdft project
      <dest>   is the directory of destination
    options:
      --beautify  beautify the code (experimental, use vim)
      --graph=<routine1,routine2,...> or --graph=all (experimental)
                  build the graph of calls for the <routine> in the file 'routine.dot'
      --graph=directories
                  build the graph of interdependences between directories
                  (need the package graphviz)
      --help      display this message
      --lint      complete analysis (experimental)
      --verbose   display more information
""")
    sys.exit(error)


#Operating System (used for os.path and os.system)
import os
#filename match as Unix shell
import fnmatch
#stderr and stdout
import sys
#Options
import getopt
#Regular expressions
import re
#To dump python object in a file
import cPickle


#Check the version of python
version_info = map(int,sys.version_info[0:3])
if  version_info < [2,3,0]:
    sys.stderr.write("Detected version %d.%d.%d\n" % tuple(version_info))
    sys.stderr.write("Minimal required version is python 2.3.0")
    sys.exit(1)

#If the class 'set' does not exist use instead Set
if 'set' not in dir(__builtins__):
    from sets import Set as set

#Indentation (one blank character)
indent = " "*2

#Prefix for the modules containing the interfaces
prefix = "interfaces_"

#To add before and after statements added by abilint
abilint_start = "!This section has been created automatically by the script Abilint (TD).\n" \
                + "!Do not modify the following lines by hand.\n"
abilint_end =  "!End of the abilint section\n"

use_before = abilint_start
notuse_before = abilint_start

#i# use_after = "#endif\n" + abilint_end
#i# notuse_after = use_after
use_after = abilint_end
notuse_after = abilint_end

#End of local variables (for abirules)
re_end_variables = re.compile("[\n]*![ ]*[*]{6,}.*\n",re.MULTILINE)
end_variables = "\n!"+72*"*"+"\n"

#Intrinsic routines of the Fortran
intrinsic_routines = [ "cpu_time",
                       "date_and_time",
                       "mvbits",
                       "random_number",
                       "random_seed",
                       "system_clock" ]
#Intrinsic functions of the Fortran
intrinsic_functions = [ "abs", "achar", "acos", "adjustl", "adjustr", "aimag", "aint", "all",
                       "allocated", "anint", "any", "asin", "associated", "atan", "atan2",
                       "bit_size", "btest", "ceiling", "char", "cmplx", "conjg", "cos", "cosh",
                       "count", 
                       "dabs", "dble", "dcmplx", "dfloat", "digits", "dim", "dot_product", "dprod",
                       "dreal", "eoshift", "epsilon", "exp", "exponent", 
                       "float", "floor", "fraction", "huge", "iachar", "iand",
                       "ibclr", "ibits", "ibset", "ichar", "ieor", "index", "int", "ior", "ishft",
                       "ishftc", "kind", "lbound", "len", "len_trim", "lge", "lgt", "ll", "llt",
                       "log", "logical", "log10", "matmul", "max", "maxexponent", "maxloc",
                       "maxval", "merge", "min", "minexponent", "minloc", "minval", "mod", "modulo",
                       "nearest", "nint", "not", "null", "pack", "pad", "precision", "present", "product",
                       "radix", "range", "real", "repeat", "reshape", "rrspacing", "scale", "scan",
                       "selected_int_kind", "selected_real_kind", "set_exponent", "shape", "sign",
                       "sin", "sinh", "size", "spacing", "spread", "sqrt", "tan", "tanh",
                       "tiny", "transfer", "transpose", "trim", "ubound", "unpack", "verify",
                       "datan", "dcos", "dexp", "dlog", "dsin", "dsqrt", "dtan"
                      ]

#Comment "!Arguments"
comment_arguments = "!Arguments " + "\n"
#Comment "!Local variables"
comment_local_variables = "!Local variables " + "\n"
#Last line of declaration
comment_declaration = "!" + "*"*60

#Regular expression to detect this section
re_section_abilint = re.compile("[ \n]*" + abilint_start[:11] + ".+?" \
                        + abilint_end.rstrip() + "[^\n]*[\n]*", \
                        re.MULTILINE+re.IGNORECASE+re.DOTALL)
re_section_abilint_start = re.compile("[ \n]*" + abilint_start[:11] + ".+?" \
                        + "#if" + "[^\n]*[\n]*", \
                        re.MULTILINE+re.IGNORECASE+re.DOTALL)
re_section_abilint_end = re.compile("[ \n]*" + "#endif" + ".+?" \
                        + abilint_end.rstrip() + "[^\n]*[\n]*", \
                        re.MULTILINE+re.IGNORECASE+re.DOTALL)

#Regular Expressions
#-------------------
#Detect !Arguments
re_arguments = re.compile("[\n ]*!Argument[^\n]+\n")
#Detect !Local
re_local_variables = re.compile("[\n ]*!Local[^\n]+\n")
#Remove use replaced by the system of interfaces
re_use_replaced = re.compile('^[ ]+use[ ]*(defs_berry|defs_dyson|defs_interfaces|defs_xc|defs_xderive|defs_xfuncmpi)[ ]*?\n',re.MULTILINE+re.IGNORECASE)


#Function to split variable in declaration statements (critics and robust but ugly)
def split_variables(string):
    """Split declaration into variables: [ (varname,declaration), ...]
    Accept comments and \n (!)"""
    a = list()
    varname = ""
    declaration = ""
    it = iter(string)
    for c in it:
        if c == "=":
            #Parameters: iterate up to a comma or the end
            declaration += c
            for c in it:
                #Check if inside an array (/
                if c == "(":
                    par = 1
                    while par:
                        declaration += c
                        c = it.next()
                        if c == "(":
                            par += 1
                        elif c == ")":
                            par -= 1
                    declaration += c
                    continue
                if c == '"' or c == "'":
                    #Iterate up to the last quote
                    quote = c
                    declaration += c
                    c = it.next()
                    while c != quote:
                        declaration += c
                        c = it.next()
                    declaration += c
                elif c == " ":
                    #Remove
                    c = ""
                elif c == ",":
                    #Remove "," and new variable
                    a.append( (varname.lower(),varname+declaration) )
                    varname = ""
                    declaration = ""
                    c = ""
                    #Finish this loop
                    break
                elif c ==  "!":
                    #Comment: Finish all
                    break
                else:
                    declaration += c
            #Continue except for a comment
            if c != "!":
                continue
        #Inside a declaration
        if c == "(":
            par = 1
            while par:
                declaration += c
                c = it.next()
                if c == "(":
                    par += 1
                elif c == ")":
                    par -= 1
            declaration += c
            continue
        if c == " ":
            #Remove
            c = ""
        if c == ",":
            #Remove "," and new variable
            a.append( (varname.lower(),varname+declaration) )
            varname = ""
            declaration = ""
            c = ""
        if c == "!":
            #Comment: finished
            break
        varname += c
    #Append and return
    a.append( (varname.lower(),varname+declaration) )
    return a


#Functions to split the code (experimental but quite robust)
def split_code(text):
    "Split a statement of code, do not keep floats and structure components"
    logical = [ "and", "eq", "eqv", "false", "ge", "gt", "le", "lt", "ne", "neqv", "not", "or", "true" ]
    re_eq = re.compile("[0-9]*[.]e")
    statements = list()
    statement = list()
    word = ''
    #True if continuation line
    continuation = False
    iter_text = iter(text.lower())
    for char in iter_text:
        #Numbers
        while char.isdigit() or char == '.':
            word += char
            char = iter_text.next()
            if char == 'd' or char == 'e':
                word += char
                char = iter_text.next()
                if char == '-' or char == '+':
                    word += char
                    char = iter_text.next()
                elif not char.isdigit():
                    #Assume .e q
                    if re_eq.match(word) and char == 'q':
                        if len(word) > 2:
                            statement.append(word[:-2])
                        word = '.'
                        char = 'eq'
                    elif char == 's':
                        #We assume that this is the format es (as 3es)
                        char += word
                    else:
                        sys.stdout.write('text=%s\n' % text)
                        sys.stdout.write('statement=%s\n' % statement)
                        sys.stdout.write('word=%s, char=%s\n' % (word,char))
                        sys.stdout.write('Error in split_code\n')
                        sys.exit(1)
            elif char in 'aglno' and word[-1] == '.':
                statement.append(word[:-1])
                word = '.'
            elif char == "_":
                #1._
                statement.append(word)
                word = char
                char = ''
        if word:
            statement.append(word)
            continuation = False
            word = ''
        #Alpha-Alphanumeric
        while char.isalnum() or char == '_':
            word += char
            char = iter_text.next()
        if word:
            statement.append(word)
            continuation = False
            word = ''
        if char == '%':
            #Remove the next word
            char = iter_text.next()
            while char.isalnum() or char == '_' or char == '%':
                char = iter_text.next()
        if char == ';' or char == '\n' or char == '!' or char == '#':
            #End of statement
            if (not continuation) and statement:
                statements.append(statement)
                statement = list()
            #Remove after ! and #
            if char == '!' or char == '#':
                while char != '\n':
                    char = iter_text.next()
            #New line
            continuation = False
        elif char == "'" or char == '"':
            #Remove boolean operator or characters
            end = char
            char = iter_text.next()
            while char != end:
                char = iter_text.next()
        elif char == '&':
            #Ignore it and flag continuation
            continuation = True
        elif char == ' ' or char == '\t':
            #Remove
            pass
        else:
            #Check if continuation and remove it because it is & at the beginning of a line
            continuation = False
            if char == '.':
                ll = statement[-1]
                if ll in logical and statement[-2] == '.':
                    #Detect logical operator
                    char = '.%s.' % ll
                    statement.pop()
                    statement.pop()
            statement.append(char)
    return statements


#Indent the Fortran code
def indent_code(text,n=0):
    "Indent a Fortran code"
    re_indented = re.compile("^[ \t]*(interface|function|((pure)?[ ]*subroutine))",re.IGNORECASE)
    re_end = re.compile("^[ \t]*end",re.IGNORECASE)
    new = ""
    for line in text.splitlines():
        if re_end.match(line):
            n -= 1
        #sys.stdout("[%d][%s]\n" % (n,n*indent+line.lstrip()[:-1]))
        new += n*indent+line.lstrip()+"\n"
        if re_indented.match(line):
            n += 1
    return new


#Build the declaration
def build_declaration(dict_vars):
    "Build the declaration from a dictionary"
    arguments = list()
    locals = list()
    declaration = ""
    for var in dict_vars.values():
        if var.is_argument:
            arguments.append( (var.order(),var.lower) )
        else:
            locals.append( (var.order(),var.lower) )
    arguments.sort()
    locals.sort()
    if arguments:
        declaration = comment_arguments
    line = ""
    type = ""
    for (order,name) in arguments:
        arg = dict_vars[name]
        (decl,var) = dict_vars[name].build_declaration()
        if len(line) > 80 or (decl != "type" and decl != "interface"):
            if line:
                declaration += line+"\n"
            if decl == "interface":
                line = var
                type = ""
            else:
                type = decl 
                line = indent+"%s :: %s" % (decl,var)
        else:
            if line:
                line += "," + var
    if line:
        declaration += line+"\n"
    if locals:
        declaration += comment_local_variables
    line = ""
    type = ""
    for (order,name) in locals:
        (decl,var) = dict_vars[name].build_declaration()
        if len(line) > 80 or  decl != type:
            if line:
                declaration += line+"\n"
            type = decl 
            line = indent+"%s :: %s" % (decl,var)
        else:
            if line:
                line += "," + var
    if line:
        declaration += line+"\n"
    return declaration


#Class for a project
class Project:
    "Class to define a project which is a set of directories and files."
    def __init__(self,dir,pat_dir=['*'],pat_file=["*"],exclude=[],name="",logfile="project.log",\
                 given_include=dict(),File_Class=dict(),style_comment="robodoc"):
        "Initialisation"
        self.ROOT = dir
        self.re_ROOT = re.compile("^"+self.ROOT+"/")
        self.name = name
        #Instantiation of the class message
        self.message = Message(verbose=verbose,logfile=logfile)
        #List of directories
        self.dirs = list()
        #List of files in each directory
        self.dirsfiles = dict()
        self.dirsfiles[dir] = list()
        #Dictionary of files referenced by file (the name should be unique: file.dir/file.name)
        self.files = dict()
        #Dictionary of routines (the name should be unique)
        self.routines = dict()
        #Dictionary of modules
        self.modules = dict()
        #Includes given by the script if missing
        self.given_include = given_include
        #Dictionary of class of files for each pattern
        self.File_Class = File_Class
        #Add directories and files
        self.add(self.ROOT,pat_dir,pat_file,exclude)
        #List of functions
        self.functions = list()
        #Dictionary of fortran type
        self.ftypes = dict()
        #Number of copied files
        self.copied = 0
        #Style of the comments for the documentation
        self.style_comment = style_comment
        #Data which are saved for a new process
        self.cached = dict()
        #Message in the case of routines with the same name and building of interfaces
        self.message_interfaces = ""
        #Statistics
        self.message.write("Project %s (%d directories, %d files)\n" % \
                         (self.name,len(self.dirs),len(self.files)),verbose=-10)
    #
    def add(self,dir,pat_dir=['*'],pat_file=['*'],exclude=[],read_only=False,File_Class=None):
        "Add directories and files with a pattern"
        try:
            files = os.listdir(dir)
        except OSError:
            self.message.error("The directory '%s' does not exist!\n" % dir)
            return
        #Check if files are in pattern and not in exclude
        #Use [:] because we modify files in the block
        for file in files[:]:
            #Check also as */*
            d1 = "%s/%s" % (os.path.basename(dir),file)
            OK = False
            if os.path.isdir("%s/%s" % (dir,file)):
                pattern = pat_dir
            else:
                pattern = pat_file
            for pat in pattern:
                OK = OK or fnmatch.fnmatch(file,pat) or fnmatch.fnmatch(d1,pat)
            if not OK:
                files.remove(file)
                continue
            for excl in exclude:
                if fnmatch.fnmatch(file,excl) or fnmatch.fnmatch(d1,excl):
                    files.remove(file)
        for file in files:
            dd = "%s/%s" % (dir,file)
            if os.path.isdir(dd):
                self.message.write("[add dir %s ---> " % dd)
                if not self.dirsfiles.has_key(dd):
                    self.dirs.append(dd)
                    self.dirsfiles[dd] = list()
                #Add all files in this sub-directory
                self.add(dd,pat_dir,pat_file,exclude,read_only,File_Class)
                self.message.write(" <--- add dir %s]\n" % dd)
            elif os.path.isfile(dd):
                self.add_file(dir,file,read_only=read_only,File_Class=File_Class)
    #
    def add_file(self,dir,file,create=False,read_only=False,File_Class=None):
        "Add a file in the project and return the file."
        fullname = "%s/%s" % (dir,file)
        if not self.dirsfiles.has_key(dir):
            self.dirs.append(dir)
            self.dirsfiles[dir] = list()
        if not file in self.dirsfiles[dir]:
            self.dirsfiles[dir].append(file)
        else:
            #Already created: Do nothing
            return self.files[fullname]
        tt = file.split(".")
        if len(tt) > 1:
            suffix = tt[-1]
        else:
            suffix = ""
        self.message.write("(%s)" % file)
        if not File_Class:
            File_Class=File
            for (pattern,file_class) in self.File_Class:
                if fnmatch.fnmatch(file,pattern):
                    File_Class=file_class
                    break
        self.files[fullname] = File_Class(file,dir=dir,read_only=read_only,message=self.message)
        if not create:
            self.files[fullname].read()
        #Return the File instance
        return self.files[fullname]
    #
    def add_routine(self,routine):
        "Add a routine to the project (routine.lower have to be defined)"
        if routine.lower in self.routines:
            old = self.routines[routine.lower]
            #Then add a fatal message in the case of building interfaces
            self.message_interfaces = self.message_interfaces \
                    + "Two subroutines have the same name:\n" \
                    + "   %s in %s/%s\n" % (old.name,old.dir,old.file) \
                    + "   %s in %s/%s\n" % (routine.name,routine.dir,routine.file)
            #Change the name of the routine (to be improved)
            name = routine.lower
            i = 1
            while "%s_%d" % (name,i) in self.routines:
                i += 1
            self.routines["%s_%d" % (name,i)] = routine
        else:
            self.routines[routine.lower] = routine
    #
    def add_use_interface(self):
        "Add 'use interface_' at each file"
        self.message.section("Add 'use interface_xxx' in each routine...")
        #First we add a module for each routine
        for routine in self.routines.values():
            if not routine.module:
                routine.module = "%s%s" % (prefix,os.path.basename(routine.dir))
        #Do for each routine stored in self.routines
        for (name,routine) in self.routines.items():
            self.message.write("(%s/%s <%s>)\n" % (routine.dir,routine.file,name))
            routine.add_use_interface(self)
        self.message.done()
    #
    def analyze_all(self,exclude="a"*40):
        "Analyze all the project except the files which contain exclude"
        self.message.section("Analyze all the project...")
        for file in self.files.values():
            #print file.name
            if exclude not in file.name:
                file.analyze(self)
        self.message.done()
        #Analyze all variables of all routines
        self.analyze_variables()
        #Build the list of arguments of all routines
        self.get_arguments()
        #Determine the list of functions
        self.set_functions()
    #
    def analyze_execution(self,edition=False):
        "Analyze the execution statements of all routines"
        self.message.section("Analyze the execution statements of routines...")
        for (name,routine) in self.routines.items():
            self.message.write("(%s)" % name)
            routine.analyze_execution(self,edition=edition)
        self.message.done()
    #
    def analyze_comments(self,exclude="a"*40,edition=False):
        "Analyze the comments of all the files except the files which contain exclude"
        self.message.section("Analyze the comments...")
        for file in self.files.values():
            if isinstance(file,File_F90) and exclude not in file.name:
                file.analyze_comments(self,edition=edition)
        self.message.done()
    #
    def analyze_directories(self):
        "Analyze the hierarchy of directories given by 'rank_dir'"
        #First establish the hierarchy
        self.hierarchy()
        for routine in self.routines.values():
            rank_routine = routine.rank
            for called in routine.called:
                try:
                    if rank_routine < self.routines[called].rank:
                        self.message.error("Misplacing <%s> (%s) calls <%s> (%s)" \
                                % (routine.name,routine.dir,called,self.routines[called].dir))
                except KeyError:
                    pass
    #
    def analyze_variables(self):
        "Analyze the variables of routines and modules"
        self.message.section("Analyze the variables of all modules and routines...")
        #Do for each module stored in self.modules
        for (name,module) in self.modules.items():
            self.message.write("(%s/%s <%s>)" % (module.dir,module.file,name))
            module.analyze_variables(self)
        #Do for each routine stored in self.routines
        for (name,routine) in self.routines.items():
            if not isinstance(routine,Generic_Routine):
                #Analyze only for not Generic_Routine
                self.message.write("(%s/%s <%s>)" % (routine.dir,routine.file,name))
                routine.analyze_variables(self)
        self.message.done()
    #
    def backup(self):
        "Backup the code"
        self.message.section("Backup all the code before changing it...")
        for file in self.files.values():
            file.backup()
        self.message.done()
    #
    def beautify(self):
        "Improve the appearance of the code"
        self.message.section("Beautify the code...\n")
        for file in self.files.values():
            file.beautify()
        self.message.done()
    #
    def copy_all(self,NEW,only_if_changed=True):
        "Copy all the project in a newdir."
        self.message.section("Copy all the project in the directory '%s'..." % NEW)
        for file in self.files.values():
            newdir = file.dir.replace(self.ROOT+"/",NEW+"/")
            if not os.path.exists(newdir):
                os.makedirs(newdir)
            file.write(newdir,file.name,fd=None,only_if_changed=only_if_changed)
            self.copied += file.copied
        if self.copied <= 1:
            self.message.write("(%d copied file)" % self.copied)
        else:
            self.message.write("%d copied files)" % self.copied)
        self.message.done()
        #Statistics
        self.message.final()
    #
    def cache_load(self,NEW,OLD=None,nocache=False):
        "Load a cached file for the project in order to decrease the execution of abilint next time"
        cached_name = "%s/.abilint" % NEW
        if os.path.exists(cached_name):
            if nocache:
                #Remove the cache file
                self.message.section("Remove the cached file '%s'..." % cached_name)
                os.remove(cached_name)
                self.message.done()
                return
            try:
                self.message.section("Read the cached file '%s'..." % cached_name)
                self.cached = cPickle.load(open(cached_name,"r"))
                self.message.done()
            except EOFError, cPickle.UnplickingError:
                self.cached = dict()
            #Check if
            # - the ROOT directory is the same to use safely the cache
            # - the python version is the same (bug from MT)
        version_info = map(int,sys.version_info[0:3])
        if self.cached.has_key("ROOT") and self.cached["ROOT"] == OLD and \
           self.cached.has_key("python_version") and self.cached["python_version"] == version_info:
            #We can use the cache
            pass
        else:
            #Do not use the cache
            self.cached = dict()
    #
    def cache_save(self,NEW):
        "Save a cached file for the project in order to decrease the execution of abilint next time"
        #Add the ROOT directory
        self.cached["ROOT"] = self.ROOT
        #Add the version of python
        self.cached["python_version"] = map(int,sys.version_info[0:3])
        for (name,routine) in self.routines.items():
            self.cached[name] = routine.cache_save()
        cached_name = "%s/.abilint" % NEW
        self.message.section("Write the cached file '%s'..." % cached_name)
        cPickle.dump(self.cached,open(cached_name,"w"))
        self.message.done()
    #
    def cache_use(self):
        "Use the information in the cache"
        for routine in self.routines.values():
            routine.cache_use(self)
    #
    def dependencies(self,inside_name,outside_name):
        "Build the list of dependencies for all routines inside and outside the directory"
        self.message.section("Build the list of dependencies for each file...")
        #First create a list of all modules per directories (dir_modules)
        dir_modules = dict()
        for module in self.modules.values():
            if not dir_modules.has_key(module.dir):
                dir_modules[module.dir] = [ module.name ]
            else:
                dir_modules[module.dir].append(module.name)
        #Then build the dependencies in the same directory (dependencies) and outside
        for (dir,files) in self.dirsfiles.items():
            dep = dict()
            dirs = set()
            dependencies = None
            for filename in files:
                fullname = "%s/%s" %(dir,filename)
                file = self.files[fullname]
                if isinstance(file,File_F90):
                    (dependencies,dirs_child) = file.dependencies(self)
                    dirs.update(dirs_child)
                    if dependencies:
                        dep[file.name] = map(lambda x: x.file,dependencies)
                        dep[file.name].sort()
            #Create the file 'inside_name' only if modules are inside the directories
            if dir_modules.has_key(dir) or dirs:
                self.message.write("[%s:" % dir,verbose=-10)
            if dir_modules.has_key(dir):
                #Build the file inside_name
                self.message.write("(%s)" % inside_name,verbose=-10)
                depfile = self.add_file(dir,inside_name,create=True)
                #Add header of the file
                depfile.add_code(head_dependencies % {"dir": dir, 
                                                      "message": "(inside the directory)"})
                #Add to clean files .mod
                depfile.add_code("CLEANFILES += ")
                linmod = dir_modules[dir]
                linmod.sort()
                for module in linmod:
                    depfile.add_code("\\\n\t%s.$(MODEXT) " % module)
                depfile.add_code("\n")
            if dep:
                files = dep.keys()
                files.sort()
                for file in files:
                    #Remove all suffixes (ex. .F90 and .F90.in)
                    depfile.add_code("\n%s.%s: " % (file.split(".")[0],"$(OBJEXT)"))
                    for module in dep[file]:
                        depfile.add_code("%s.%s " % (module.split(".")[0],"$(OBJEXT)"))
                    depfile.add_code("\n")
            if dirs:
                #Build the file outside_name
                self.message.write("(%s)" % outside_name,verbose=-10)
                depfile = self.add_file(dir,outside_name,create=True)
                #Add header of the file
                depfile.add_code(head_dependencies % {"dir": dir,
                                                      "message": "(outside the directory)"})
                #Add the current directory
                depfile.add_code("include_dirs = [")
                comma = ""
                for d in dirs:
                    depfile.add_code('%s\n\t"%s"' % (comma,d.replace("%s/src/" % OLD,"")))
                    comma = ","
                depfile.add_code("]\n")
            if dir_modules.has_key(dir) or dirs:
                self.message.write("]",verbose=-10)
        self.message.done()
    #
    def dump_ftype(self,name,filename):
        "Dump with cPickle a Fortran type in a file"
        if name in self.ftypes:
            self.message.section("Write the fortran type '%s' into the file '%s'..." % (name,filename))
            cPickle.dump(self.ftypes[name].dict_vars,open(filename,"w"))
            self.message.done()
        else:
            self.message.fatal("The Fortran type '%s' is not inside the project!\n" % name)
    #
    def find_file(self,filename):
        "Find the instance File given by its name"
        for file in self.files.values():
            if file.name == filename:
                return file
        return None
    #
    def generate(self,format,pattern_dir='.*'):
        "Generate a file which contains the interfaces of all routines in dir"
        self.message.section("Generate all signature of directories '%s'..." % pattern_dir)
        re_dir = re.compile(pattern_dir)
        for dir in self.dirs:
            if re_dir.match(dir):
                signature = ""
                generated = format % os.path.basename(dir)
                fullgenerated = "%s/%s" % (dir,generated)
                file=self.files[fullgenerated]
                for file in self.dirsfiles[dir]:
                    if file != generated:
                        fullname = "%s/%s" % (dir,filename)
                        signature += self.files[fullname].signature()
                #Create the file if it doesn't exist
                if not self.files.has_key(fullgenerated):
                    file = self.add_file(dir,generated,create=True,File_Class=File_F90_Library)
                else:
                    file = self.files[fullgenerated]
                #Add code to the file
                file.add_code(signature)
        self.message.done()
    #
    def get_arguments(self):
        "Build the list of arguments of all routines"
        self.message.section("Build the list of arguments in each routine...")
        #Do for each routine stored in self.routines
        for (name,routine) in self.routines.items():
            self.message.write("(%s/%s <%s>)" % (routine.dir,routine.file,name))
            routine.get_arguments(self)
        self.message.done()
    #
    def graph(self,names,graph_excluded=list()):
        "Build some graphs of calls"
        #Build a graph dictionary
        if names == "directories":
            #Determine the rank of all routines
            self.analyze_directories()
            filename = "directories.dot"
        else:
            filename = "bigdft.dot"
        fd = open(filename,"w")
        if names == "directories":
            #Build the graph of module interdependences
            self.message.section("Build the graphs of all directories...")
            dict_graph = dict()
            #Dictionary of the color of nodes (red, wrong calls, green OK)
            dict_color = dict()
            for routine in self.routines.values():
                dir = os.path.basename(routine.dir)
                if dir not in dict_graph:
                    dict_graph[dir] = set()
                    dict_color[dir] = "green"
                for called in routine.called:
                    try:
                        called_dir = self.routines[called].dir.split('/')[-1]
                        num = rank_dir(called_dir)
                        dict_graph[dir].add(called_dir)
                        if routine.rank < num:
                            dict_color[dir] = "red"
                            self.message.warning("Misplacing <%s> (%s) calls <%s> (%s)" \
                                    % (routine.name,dir,called,called_dir),verbose=-10)
                    except KeyError:
                        pass
            self.graph_build("directories",fd,dict_graph,dict_color,graph_excluded=graph_excluded)
        elif names == "all":
            #Build the graph of all routines
            self.message.section("Build the graphs of all routines (parents and children)...")
            names = self.routines.keys()
            names.sort()
            for name in names:
                dict_graph = dict()
                self.message.write("(%s)" % name,verbose=-10,flush=True)
                self.routines[name].graph_children(dict_graph,routines=self.routines)
                self.routines[name].graph_parents(dict_graph,routines=self.routines)
                self.graph_build(name,fd,dict_graph,graph_excluded=graph_excluded)
        else:
            self.message.section("Build the graph for the routines...")
            for name in names.split(','):
                dict_graph = dict()
                self.message.write("(%s)" % name,verbose=-10,flush=True)
                if name not in self.routines.keys():
                    self.message.error("The routine %s does not exist: no graph\n" % name)
                    return
                self.routines[name].graph_children(dict_graph,routines=self.routines,recursive=True)
                self.routines[name].graph_parents(dict_graph,routines=self.routines,recursive=True)
                self.graph_build(name,fd,dict_graph,graph_excluded=graph_excluded)
        fd.close()
        self.message.write("(%s)" % filename,verbose=-10)
        self.message.done()
    #
    #Build the tree and give an identifier at each node
    def graph_build(self,name,fd,dict_graph,dict_color=None,graph_excluded=list()):
        "Build the graph"
        fd.write("digraph routine_%s {\n" % name)
        fd.write("   node [shape = record,height=.1,style=filled];\n")
        #Build the dict of colors
        if not dict_color:
            dict_color = dict()
            for key in dict_graph.keys():
                dict_color[key] = "white"
                for called in dict_graph[key]:
                    dict_color[called] = "white"
        id = 0
        node_graph = dict()
        for (calling,called_set) in dict_graph.iteritems():
            if calling not in node_graph and calling not in graph_excluded:
                id += 1
                node_graph[calling] = id
                fd.write(3*' '+'n%d [label="%s",color=%s];\n' \
                    % (id,calling,dict_color[calling]))
            for called in called_set:
                if called not in node_graph and called not in graph_excluded:
                    id += 1
                    node_graph[called] = id
                    fd.write(3*' '+'n%d [label="%s",color=%s];\n' \
                        % (id,called,dict_color[called]))
        for (calling,called_set) in dict_graph.iteritems():
            for called in called_set:
                if called not in graph_excluded:
                    fd.write(3*' '+'n%d -> n%d;\n' % (node_graph[calling],node_graph[called]))
        fd.write("}\n")
    #
    def has_module(self,module):
        "Check if the module exists and check if also in the given modules."
        return self.modules.has_key(module)
    #
    def hierarchy(self):
        "Establish the hierarchy of routines: a routine rank_dir should be defined."
        for routine in self.routines.values():
            dir = os.path.basename(routine.dir)
            routine.rank = rank_dir(dir)
    #
    def interfaces_all(self,dir_mod,exclude,warning=""):
        "Build the interfaces of all files"
        self.message.section("Build the interfaces...")
        #First test if there are routines with the same name
        if self.message_interfaces:
            self.message.fatal("\n"+self.message_interfaces)
        #Do for each directory
        for dir in self.dirs:
            main = os.path.basename(dir)
            #exclude contains the excluded directories
            if main in exclude:
                continue
            #Directory where the interface files are stored
            dir_module = "%s/src/%s" % (self.ROOT,main)
            if ( not os.path.exists(dir_module) ):
                dir_module = "%s/%s" % (self.ROOT,dir_mod)
            module = "%s%s" % (prefix,main)
            file_module = "%s.F90" % module
            pfile = self.add_file(dir_module,file_module,create=True)
            directory = self.re_ROOT.sub('',dir)
            pfile.code = head_interface % { \
                    'name': module,
                    'description': "in the directory %s" % directory,
                    'warning': warning }
            #Add interfaces in the file
            flist = self.dirsfiles[dir]
            #Sort alphabetically
            flist.sort()
            for file in flist:
                fullname = "%s/%s" % (dir,file)
                pfile.code += self.files[fullname].interface()
            pfile.code += "end module %s\n" % module \
                          + "!!***\n"
            #Analyze this file
            pfile.analyze(self)
        self.message.done()
    #
    def no_fatal(self):
        "Routine which prepares the project to have no fatal error in order to generate interfaces"
        #Remove possible message for duplicated routines
        self.message_interfaces = ""
        #Remove all children
        self.remove_children()
    #
    def remove_children(self):
        "Remove children of all files"
        for file in self.files.values():
            file.remove_children()
    #
    def remove_code(self,pattern_file):
        "Remove code except comments of files"
        re_file = re.compile(pattern_file)
        for (dir,files) in self.dirsfiles.items():
            for file in files: 
                if re_file.match(file):
                    fullname = "%s/%s" % (dir,file)
                    self.files[fullname].remove_code()
    #
    def remove_files(self,pattern):
        "Remove files which correspond to the given pattern (do not clean properly!)"
        self.message.section("Remove the files with pattern '%s'..." % pattern)
        for (dir,files) in self.dirsfiles.items():
            for file in files: 
                if pattern in file:
                    fullname = "%s/%s" % (dir,file)
                    self.message.write("[%s]" % (fullname))
                    del self.files[fullname]
                    self.dirsfiles[dir].remove(file)
        self.message.done()
    #
    def set_children_parents(self):
        "Build the list of children and parents routines"
        self.message.section("Determine children and parents of all routines...")
        for routine in self.routines.values():
            routine.set_children_parents(self)
        self.message.done()
    #
    def set_functions(self):
        "Build the list of functions"
        self.message.section("Build the list of functions...")
        self.message.write("\n")
        for routine in self.routines.values():
            if isinstance(routine,Function):
                self.message.write("Function %s -- (%s)\n" % (routine.name,routine.type))
                self.functions.append(routine.name.lower())
        #We sort the functions
        self.functions.sort()
        functions = ""
        for function in self.functions:
            #We do not detect function with only one character
            if len(function) > 1:
                functions = "%s|%s" % (functions,function)
            else:
                self.message.warning( \
                     "The function %s (%s/%s) has a too short name for the detection" \
                     % (function,self.routines[function].dir,self.routines[function].file))
        nf = len(self.functions)
        self.message.write("Number of functions: %d\n" % nf)
        #Build the regular expression to detect functions ('xx(')
        if nf  == 0:
            self.re_functions = None
        else:
            self.re_functions = re.compile("(?<![a-z])(%s)[ ]*[(]" % functions[1:],re.IGNORECASE)
        #Warning: This detection is not fully robust but works
        self.message.done()
    #
    def special(self,special_modif):
        "Do some special modifications"
        if len(special_modif) == 0:
            #Nothing to be done
            return
        self.message.write("\nSpecial transformations:\n")
        for name in special_modif:
            self.message.write("[special:%s]" % name)
            eval("self.routines[name].%s" % special_modif[name])
        self.message.write("\n")
    #
    def statistics(self):
        "Display some statistics about the project"
        self.message.write("Statistics of the project %s:\n" % self.name \
            + "%d directories, " % len(self.dirs) \
            + "%d files, " % len(self.files) \
            + "%d modules, " % len(self.modules) \
            + "%d routines (from which %d functions)" % (len(self.routines),len(self.functions)) \
            + " -- [%d files copied]\n" % self.copied,verbose=-10)
    #
    def unused_routines(self):
        "Give list of unused routines and if a routine inside a generic routine is not called"
        #Do for each routine stored in self.routines
        unused = list()
        for (name,routine) in self.routines.items():
            #Only routine (not program) and not in lib directory 
            if len(routine.calling) == 0 and routine.nogeneric and routine.rank < 100 and routine.rank >= 0:
                unused.append((routine.dir,routine.file,name))
            if len(routine.calling) != 0 and (not routine.nogeneric) and routine.rank < 100 and routine.rank >= 0:
                #A routine part of a generic routine is called directly
                self.message.error("<%s> (%s/%s), part of a generic routine, is called directly by %s" % \
                        (name,routine.dir,routine.file,list(routine.calling)))
        unused.sort()
        self.message.section("List of unused routines:")
        for lists in unused:
                self.message.write("(%s/%s <%s>)" % lists,verbose=-10)
        self.message.done()


#General class handling File and the other ones.
class Structure:
    "Define a structure (a routine, header, etc) which can have a portion of code"
    #Comments and also preprocessing commands
    re_comment_match = re.compile("^([ \t]*(([!#].*)|[ ]*)\n)+")
    #Detect the beginning and the end of a block
    re_sub_start = re.compile('^[ \t]*(module|program|(((pure|recursive)[ ]*)*subroutine)|(([^!\'"\n]*?)function))',re.IGNORECASE)
    re_sub_end   = re.compile('^[ \t]*end[ \t]*(function|module|program|subroutine|\n)',re.IGNORECASE)
    #
    def __init__(self,name=None,parent=None,message=None):
        "Initialisation"
        self.name = name
        if name:
            self.lower = self.name.lower()
        self.parent = parent
        if parent:
            #Add self in child parent
            self.message = self.parent.message
            #Add to parent's children (if it exists)
            parent.children.append(self)
        elif message:
            self.message = message
        else:
            sys.stderr.write("Structure %s (%s) has no message class!\n" \
                    % (str(self.__class__),self.name))
            sys.exit(1)
        self.children = list()
        #Code of the structure
        self.code = ""
        self.code_backup = ""
        #Is analyzed
        self.is_analyzed = False
    #Add code
    def add_code(self,code):
        "Add code"
        self.code += code
    #Analyze
    def analyze(self):
        "Check only if it is already analyzed"
        if self.is_analyzed:
            self.message.fatal("The structure %s is already analyzed" % self.name) 
    #Ancestry (debugging)
    def ancestry(self,n):
        if self.parent:
            self.parent.ancestry(n+1)
    #Backup
    def backup(self):
        "Do a copy of the code"
        self.code_backup = self.code
        for child in self.children:
            child.backup()
    #Beautify
    def beautify(self):
        "Beautify the children  "
        for child in self.children:
            child.beautify()
    #Compare
    def cmp(st1,st2):
        "Compare two structures"
        if st1.lower < st2.lower:
            return -1
        else:
            return 1
    #Check if the code and the backup are identical
    def is_changed(self):
        changed = self.code != self.code_backup
        if changed:
            return True
        #Test the children
        for child in self.children:
            changed = child.is_changed()
            if changed:
                return True
        return False
    #String representation
    def __str__(self):
        "Display all the code with children"
        #The code is inserted before the child codes
        code = self.code
        for child in self.children:
            code +=str(child)
        return code
    #Read
    def read(self,fd):
        "Read from a file descriptor"
        #We ban \r !!
        self.code = fd.read().replace('\r','')
    #Remove children
    def remove_children(self):
        "Remove all children"
        self.children = list()
    #Replace
    def replace(self,old,new):
        "Replace the regular expressions 'old' by the string 'new'"
        self.code = re.sub(old,new,self.code)
    #Write
    def write(self,fd,only_if_changed):
        "write into a file descriptor"
        if only_if_changed:
            #Check if the code is changed
            if not self.is_changed():
                #Do nothing
                return
        #The code is inserted before the child codes
        fd.write(self.code)
        for child in self.children:
            child.write(fd,only_if_changed)


class File(Structure):
    "Class to define a file (for all files) which is a code with a structure"
    def __init__(self,name,dir=".",read_only=False,message=None):
        "Initialisation"
        Structure.__init__(self,name=name,parent=None,message=message)
        self.dir = dir
        self.copied = 0
        self.read_only = read_only
        #File which contains the structure
        self.file = self.name
        #Timestamp for cache (default None)
        self.timestamp = None
    #Analyze
    def analyze(self,project=None):
        "Do nothing"
        if self.is_analyzed:
            self.message.fatal("The file %s/%s is already analyzed" % (self.dir,self.file)) 
    #Interface
    def interface(self):
        "Do nothing"
        return ""
    #
    def read(self,fd=None):
        "Read the file."
        if fd == None:
            fd = open("%s/%s" % (self.dir,self.name))
        #Time stamp i.e. the time of the last modification
        self.timestamp = os.path.getmtime(fd.name)
        #We ban \r !!
        self.code = fd.read().replace('\r','')
        fd.close()
    #
    def write(self,dir,file,fd=None,only_if_changed=True):
        "Write the file as '%s/%s' % (dir,file)."
        #If read_only, return
        if self.read_only:
            return
        #Check if the code is changed
        if only_if_changed:
            if not self.is_changed():
                #Do nothing
                return
        self.message.write("[%s/%s]" % (dir,file))
        if fd == None:
            fd = open("%s/%s" % (dir,file),"w")
        #The code is inserted before the child codes
        fd.write(self.code)
        for child in self.children:
            child.write(fd,only_if_changed=False)
        fd.close()
        self.copied += 1


class File_F90(File):
    "Class to define a fortran90 file which is a code with a structure"
    #
    def __init__(self,name,dir=".",read_only=False,message=None):
        "Initialisation"
        File.__init__(self,name,dir,read_only,message=message)
        self.is_module = False
        self.module = None
    #
    def analyze(self,project):
        "Analyze the code of the file"
        File.analyze(self,project)
        self.message.write("[%s/%s:" % (self.dir,self.name),verbose=10)
        #Create an iterator
        iter_code = iter(self.code.splitlines(1))
        #All the code will be inside child structures
        self.code = ""
        struct = None
        for line in iter_code:
            if self.re_comment_match.match(line):
                if not isinstance(struct,Comment):
                    struct = Comment(parent=self)
                struct.add_code(line)
            elif self.re_sub_start.match(line):
                line_lower = line.lower()
                if "subroutine" in line_lower:
                    struct = Routine(parent=self)
                elif "module" in line_lower:
                    struct = Module(parent=self)
                    self.is_module = True
                elif "program" in line_lower:
                    struct = Program(parent=self)
                else:
                    struct = Function(parent=self)
                struct.analyze(line,iter_code,project)
            else:
                self.message.fatal("\n'%s'\n--> No detected routine (File F90 analysis)!\n" % line \
                        + "This part of code can not be parsed as Fortran file:\n" \
                        + "Analysis Error in %s/%s\n" % (self.dir,self.file))
        self.message.write("]\n",verbose=10)
    #
    def analyze_comments(self,project,edition=False):
        "Analyze comments to detect headers and footers (robodoc or doxygen style)"
        temp_structs = self.children
        self.children = []
        for struct in temp_structs:
            if isinstance(struct,Comment):
                #Detect robodc comments, split comments if necessary and add to children of self
                if project.style_comment == "robodoc":
                    struct.detect_robodoc_comment()
            else:
                #Add the structure
                self.children.append(struct)
        if project.style_comment == "robodoc":
            #We have splitted the comments. Now we build robodoc structure
            #(Robodoc_Header -- [comment] -- Routine -- Robodoc_Footer)
            self.build_robodoc_structure()
        #Finally, analyze comment
        for struct in self.children:
            if isinstance(struct,Comment):
                struct.analyze_comment(edition=edition)
    #
    def build_robodoc_structure(self):
        "Build the robodoc structure around each routine"
        robodoc_header = None
        robodoc_footer = None
        routine = None
        for struct in self.children:
            if isinstance(struct,Robodoc_Header):
                if robodoc_header != None:
                    self.message.fatal("\n%s%s" % (robodoc_header.code,struct.code) \
                            + "Two robodoc headers for a subroutine or a robodoc header without footer!\n" \
                            + "Analysis Error in %s/%s" % (self.dir,self.name))
                else:
                    robodoc_header = struct
                robodoc_footer = None
                routine = None
            elif isinstance(struct,Robodoc_Footer):
                if routine == None:
                    self.message.error("Robodoc footer without routine in %s/%s" % (self.dir,self.name))
                    if robodoc_header:
                        robodoc_header.routine = None
                    else:
                        self.message.fatal("Robodoc footer without header!\n" \
                                + "Analysis error in %s/%s" % (self.dir,self.name))
                robodoc_footer = struct
                routine = None
                robodoc_header = None
            elif isinstance(struct,Module):
                if robodoc_header:
                    if routine:
                        #There is already a routine
                        self.message.fatal("\n%s%s" % (robodoc_header.code,struct.code) \
                                + "Robodoc header for more than one subroutine or a robodoc header without footer!\n" \
                                + "Analysis Error in %s/%s" % (self.dir,self.name))
                    robodoc_header.routine = struct
                else:
                    self.message.error("Routine %s without robodoc header in %s/%s" % (struct.name,self.dir,self.name))
                routine = struct
    #
    def build_robodoc_structure_generic(self):
        "Build the robodoc structure for a library or a generic routine (Robodoc_Header [Routine] Robodoc_Footer"
        robodoc_header = None
        robodoc_footer = None
        routine = None
        for struct in self.children:
            if isinstance(struct,Robodoc_Header):
                if robodoc_header != None:
                    self.message.fatal("\n%s%s" % (robodoc_header.code,struct.code) \
                            + "Two robodoc headers for a generic subroutine!\n" \
                            + "Analysis Error in %s/%s" % (self.dir,self.name))
                elif routine:
                    self.message.fatal("\n%s%s" % (robodoc_header.code,struct.code) \
                            + "Only one robodoc header for a generic subroutine on top of the file!\n" \
                            + "Analysis Error in %s/%s" % (self.dir,self.name))
                else:
                    robodoc_header = struct
                    robodoc_header.routine = self
                    self.robodoc_header = robodoc_header
            elif isinstance(struct,Robodoc_Footer):
                if routine == None:
                    self.message.error("Robodoc footer without routine in %s/%s" % (self.dir,self.name))
                    if not robodoc_header:
                        self.message.fatal("Robodoc footer without header!\n" \
                                + "Analysis error in %s/%s\n" % (self.dir,self.name))
                if robodoc_footer != None:
                    self.message.fatal("\n%s%s" % (robodoc_footer.code,struct.code) \
                            + "Two robodoc footers for a generic subroutine!\n" \
                            + "Analysis Error in %s/%s" % (self.dir,self.name))
                robodoc_footer = struct
                robodoc_footer.routine = self
                self.robodoc_footer = robodoc_footer
            elif isinstance(struct,Routine):
                if not robodoc_header:
                    self.message.error("Generic routine without robodoc header in %s/%s" % (self.dir,self.name))
                if robodoc_footer:
                    self.message.fatal("\n%s%s" % (robodoc_header.code,struct.code) \
                            + "Robodoc footer for a generic subroutine before a routine!\n" \
                            + "Analysis Error in %s/%s" % (self.dir,self.name))
                routine = struct
    #
    def dependencies(self,project):
        "Build the dependencies"
        dependencies = set()
        dirs = set()
        for child in self.children:
            if isinstance(child,Routine) or isinstance(child,Function) or isinstance(child,Module):
                (dependencies_child,dirs_child) = child.dependencies(project)
                dependencies.update(dependencies_child)
                dirs.update(dirs_child)
        return (dependencies,dirs)
    #
    def interface(self):
        "Build the interface"
        code = ""
        for child in self.children:
            if isinstance(child,Routine) or isinstance(child,Function):
                code += child.interface() + "\n"
        return code
    #
    def remove_code(self):
        "Remove code except first comments"
        code = ""
        for line in self.code.splitlines(1):
            if self.re_comment_match.match(line):
                code += line
            else:
                #We have finished
                break
        self.code = code
    #
    def signature(self):
        "Return all signatures of routines"
        code = ""
        for child in self.children:
            if isinstance(child,Routine) or isinstance(child,Function):
                code += child.signature + "\n"
        return code
    #
    def special(self,modif):
        "Special modifications"
        for child in self.children:
            if isinstance(child,Routine):
                eval("child.%s" % modif)


class File_F90_Library(File_F90):
    "Special class to generate the interfaces from a library"
    def __init__(self,name,dir=".",read_only=False,message=None):
        "Initialisation"
        File_F90.__init__(self,name,dir,read_only,message=message)
        self.robodoc_header = None
        self.robodoc_footer = None
    #
    def build_robodoc_structure(self):
        "Build the robodoc structure for a library (Robodoc_Header [Routine] Robodoc_Footer"
        File_F90.build_robodoc_structure_generic(self)

class File_F90_Generic(File_F90):
    "Special class to generate a generic routine from Fortran90 files"
    def __init__(self,name,dir=".",read_only=False,message=None):
        "Initialisation"
        File_F90.__init__(self,name,dir,read_only,message=message)
        #Generic_Routine will be associated to the generic routine (see analyze)
        self.robodoc_header = None
        self.robodoc_footer = None
        self.generic_routine = None
    #
    def analyze(self,project=None):
        "Analyze the code of the file"
        #First analyze normally the file
        File_F90.analyze(self,project)
        #Add a routine which has the name of the file
        name = self.name.replace(".F90","")
        #Do not add to the list of children of the file
        struct = Generic_Routine(name=name,file=self)
        self.generic_routine = struct
        if project != None:
            #Add to the list of routines
            project.routines[struct.lower] = struct
        #Add the cached information
        struct.cache_use(project)
        #All routines are a part of the generic routine (used of unused_routines)
        for child in self.children:
            if isinstance(child,Routine):
                child.nogeneric = False
    #
    def build_robodoc_structure(self):
        "Build the robodoc structure for a generic interface (Robodoc_Header [Routine] Robodoc_Footer"
        File_F90.build_robodoc_structure_generic(self)
        if self.robodoc_header:
            self.robodoc_header.routine = self.generic_routine
        if self.robodoc_footer:
            self.robodoc_footer.routine = self.generic_routine
    #
    def interface(self):
        "Build a generic interface"
        name = self.name.replace(".F90","")
        code = "\n!Generic interface of the routines %s\n" % name \
                + "interface %s\n" % name
        for child in self.children:
            if isinstance(child,Routine):
                child.interface()
                code += child.signature
        code += "end interface\n" \
             + "!End of the generic interface of %s\n\n" % name
        self.generic_routine.code_interface = indent_code(code)
        return self.generic_routine.code_interface


class File_F77(File_F90):
    "Class to define a fortran 77 file (always read only)"
    def analyze(self,project):
        "Analyze the code of the file"
        File.analyze(self,project)
        #First, we transform the file into Fortran 90
        #Replace tabulation as blank characters
        self.code = self.code.replace("\t"," "*6)
        code = ""
        for line in self.code.splitlines(1):
            if len(line) < 6:
                continue
            first_char = line[0].lower()
            if first_char == "c" or first_char == "*" or first_char == "!":
                code += "!" + line[1:]
            elif first_char == "#" or first_char == "$":
                code += "#" + line[1:]
            elif line[5] != " ":
                code = code[:-1] + "&\n" + line[6:]
            else:
                code += line
        self.code = code
        #Then we analyze as a F90 File
        File_F90.analyze(self,project)


class Code(Structure):
    "Generic Class for fortran code inside the routines or the functions"
    #Comments and also preprocessing commands
    re_allcomment = re.compile("([ \t]*(([!#].*)|[ ]*)\n)+")
    #Commment inside a line
    re_comment = re.compile("[!].*")
    #Contains statement
    re_contains = re.compile('^[ \t]*contains',re.IGNORECASE)
    #Continuation
    re_continuation = re.compile("&[ \t]*(![^\n]*)?\n")
    #Include command
    re_include = re.compile('^[ ]*include.*?\n',re.MULTILINE+re.IGNORECASE)
    #
    def __init__(self,name=None,parent=None,message=None):
        "Initialisation"
        Structure.__init__(self,name=name,parent=parent,message=message)
        if parent:
            self.file = parent.file
            self.dir = parent.dir
        else:
            self.file = None
            self.dir = None


class Comment(Structure):
    "Class for comments inside Fortran files"
    #Detect robodoc header (!***)
    re_robodoc_header = re.compile("!!\*\*\*\*[a-z]\*")
    #
    def __init__(self,parent):
        "Initialisation"
        Structure.__init__(self,name="comment",parent=parent)
        self.file = self.parent.name
        self.dir = self.parent.dir
    #
    def analyze_comment(self,edition=False):
        "Analyze the comment"
        pass
    #
    def detect_robodoc_comment(self):
        "Split the comment (self) to have robodoc header or footer"
        struct = None
        for line in self.code.splitlines(1):
            if "!!" == line[0:2]:
                if "!!***" == line[:5]:
                    if self.re_robodoc_header.match(line):
                        #Add struct to the children of self.parent
                        struct = Robodoc_Header(parent=self.parent)
                    else:
                        struct = Robodoc_Footer(parent=self.parent)
                else:
                    if struct == None:
                        struct = Comment(parent=self.parent)
            else:
                if struct == None or isinstance(struct,Robodoc_Header) or isinstance(struct,Robodoc_Footer):
                    struct = Comment(parent=self.parent)
            struct.add_code(line)


class Robodoc_Header(Comment):
    "Class for the robodoc header (comment)"
    #Sections of robodoc
    robodoc_sections = ["NAME", "FUNCTION", "NOTE", "COPYRIGHT", "INPUTS", "OUTPUT", "PARENTS", "CHILDREN", "SOURCE"]
    #
    def __init__(self,parent=None):
        "Initialisation"
        Comment.__init__(self,parent)
        self.sections = dict()
        self.categories = []
    #
    def analyze_comment(self,edition=False):
        "Analyze the robodoc header"
        iter_code = iter(self.code.splitlines())
        #First line
        line = iter_code.next()
        self.type = line[6]
        self.title = line[9:-1]
        category = None
        for line in iter_code:
            sec = line[2:].strip()
            if sec != "" and sec == sec.upper():
                category = sec.upper()
                self.sections[category] = ""
                self.categories.append(category)
            elif category:
                self.sections[category] += "%s\n" % line[2:]
        #Check the coherency between the robodoc header and the routine
        if isinstance(self.parent,File_F90_Library) and self.type != "d":
            self.message.error("The type for the robodoc header '%s' is not for the file '%s'!\n" \
                    % (self.type,self.routine.file) \
                + "Analysis error in %s/%s\n" % (self.dir,self.file))
            return
        if self.type != "m" and isinstance(self.routine,Module) and not isinstance(self.routine,Routine):
            self.message.error("The type for the robodoc header '%s' is not for the module '%s'!\n" \
                    % (self.type,self.routine.name) \
                + "Analysis error in %s/%s\n" % (self.dir,self.file))
        if self.type != "f" and isinstance(self.routine,Routine) and not isinstance(self.routine,Program):
            self.message.error("The type for the robodoc header '%s' is not for the routine '%s'!\n" \
                    % (self.type,self.routine.name)\
                + "Analysis error in %s/%s\n" % (self.dir,self.file))
        if self.type != "p" and isinstance(self.routine,Program):
            self.message.error("The type for the robodoc header '%s' is not for the program '%s'!\n" \
                    % (self.type,self.routine.name) \
                + "Analysis error in %s/%s\n" % (self.dir,self.file))
    #
    def beautify(self):
        "Display Robodoc header (experimental)"
        #Define the type
        if self.routine == None:
            #Do nothing
            return
        elif isinstance(self.routine,Program):
            self.type = "p"
        elif isinstance(self.routine,Routine):
            self.type = "f"
        elif isinstance(self.routine,Module):
            self.type = "m"
        else:
            #Keep the type.
            pass
        #Define the title
        self.title = "bigdft/%s" % self.routine.name
        if isinstance(self.routine,Routine):
            #Add parents
            line = 2*indent
            routines = list(self.routine.calling)
            routines.sort()
            for name in routines:
                line += name+","
            self.sections["PARENTS"] = line[:-1]+"\n"
            #Add children
            line = 2*indent
            routines = list(self.routine.called)
            routines.sort()
            for name in routines:
                line += name+","
            self.sections["CHILDREN"] = line[:-1]+"\n"
        #Header
        self.code = "!!****%s* %s\n" % (self.type, self.title)
        for section in self.categories:
            self.code += "!! %s\n" % section
            for line in self.sections[section].splitlines(1):
                self.code += "!!%s" % line


class Robodoc_Footer(Comment):
    "Class for the robodoc footer (comment)"


class Module(Code):
    "Fortran module"
    #Detect module
    re_module = re.compile("^[ ]*module[ ]*(?P<name>\w+)",re.MULTILINE+re.IGNORECASE)
    #
    def __init__(self,name=None,parent=None,message=None):
        "Initialisation"
        Structure.__init__(self,name,parent=parent,message=message)
        #No implicit from parent (no module inside a module)
        self.implicit = None
        if parent:
            self.dir = self.parent.dir
            self.file = self.parent.file
            #None if not in a module
            self.module = self.parent.module
        else:
            self.file = None
            self.dir = None
            self.module = None
        #List of used modules
        self.use_modules = set()
        #There is a contains statement
        self.contains = True
        #True if already analyzed
        self.analyzed_variables = False
        #Gestion of the cache
        #1 - Time stamp
        if self.parent:
            self.timestamp = self.parent.timestamp
        else:
            self.timestamp = None
    #
    def analyze(self,line,iter_code,project):
        "Analyze the module"
        Code.analyze(self)
        self.Header = Code(parent=self)
        self.Header.add_code(line)
        #No argument
        self.Header.arguments = list()
        #First line contains the name of the module
        self.name = self.re_module.match(line).groupdict()['name']
        self.message.write("{%s:" % self.name)
        self.lower = self.name.lower()
        #The module belongs to him!
        self.module = self.lower
        #Add to the list of modules
        if project.modules.has_key(self.lower):
            previous = project.modules[self.lower]
            self.message.fatal("\n%s/%s: There does exist an already module <%s>\n" \
                               % (self.dir,self.file,self.name) + \
                               "in %s/%s: <%s>\n" % (previous.dir,previous.file,previous.name))
        else:
            project.modules[self.lower] = self
        #Analyze the use statement
        self.Use = Use(parent=self)
        (final_comments,line) = self.Use.analyze(iter_code)
        #Add the modules
        self.use_modules.update(self.Use.modules)
        #Test if implicit statement
        self.Implicit = Implicit(parent=self)
        #More than one line are analyzed
        line = self.Implicit.analyze(line,iter_code,comments=final_comments)
        #Analyze the declarations
        self.Declaration = Declaration(parent=self,comments=final_comments)
        line = self.Declaration.analyze(line,iter_code,self.Implicit.dict,project)
        #Now we have contains or end module
        self.contains = self.analyze_contains(line,iter_code,project)
        self.message.write("}")
    #
    def analyze_contains(self,line,iter_code,project):
        "Analyze a contains statement or end module"
        n_struct = 0
        if line:
            if self.re_contains.match(line):
                Code(parent=self).add_code(line)
            elif self.re_sub_end.match(line):
                #Detect the end of a subroutine or module: We have finished
                Code(parent=self).add_code(line)
                return False
            else:
                self.message.fatal("\n%s\n--> No detected 'contains' or 'end' statements!\n" % line \
                        + "This part of code can not be parsed as Fortran file:\n" \
                        + "Analysis Error in %s/%s:<%s>\n" % (self.dir,self.file,self.name))
        #Start the analyze of the 'contains' section
        struct = None
        for line in iter_code:
            if self.re_comment_match.match(line):
                if not isinstance(struct,Comment):
                    struct = Comment(parent=self)
                struct.add_code(line)
            elif self.re_sub_end.match(line):
                #Detect the end of a subroutine or module (before re_sub_start ortherwise trouble)
                #We have finished
                Code(parent=self).add_code(line)
                return (n_struct != 0)
            elif self.re_sub_start.match(line):
                line_lower = line.lower()
                if "subroutine" in line_lower:
                    struct = Routine(parent=self,implicit=self.Implicit.dict)
                else:
                    struct = Function(parent=self,implicit=self.Implicit.dict)
                n_struct += 1
                #Analyze the code
                struct.analyze(line,iter_code,project)
                #Add the modules used in the module to the used modules by the routine
                struct.use_modules.update(self.use_modules)
                #Add also the module used by the routine
                self.use_modules.update(struct.use_modules)
            elif self.re_include.match(line):
                #Include directive
                struct = Include(project,parent=self,line=line)
            else:
                self.message.fatal("\n%s\n--> No detected routine (Structure analysis)!\n" % line \
                        + "This part of code can not be parsed as Fortran file:\n" \
                        + "Analysis Error in %s/%s:<%s>\n" % (self.dir,self.file,self.name))
    #
    def analyze_variables(self,project):
        "Analyze the variables of the routine"
        if self.analyzed_variables:
            #Already done
            return
        else:
            self.analyzed_variables = True
        self.Declaration.analyze_variables(self.Implicit.dict,project)
        #If it is inside a module i.e. the parent is a module, then
        if isinstance(self.parent,Module):
            #Add the variables declared in the parent (module)
            if self.parent.Declaration.dict_vars == None:
                self.parent.analyze_variables(project)
            self.Declaration.dict_vars.update(self.parent.Declaration.dict_vars)
            #Add also the private variables
            self.Declaration.dict_vars.update(self.parent.Declaration.dict_vars_private)
        #Add the variables of the modules as its own variables
        for module in self.use_modules:
            a = project.modules.get(module)
            if a:
                a.analyze_variables(project)
                self.Declaration.dict_vars.update(a.Declaration.dict_vars)
            else:
                self.message.error("[%s/%s:%s] The module '%s' is missing!" \
                    % (self.dir,self.file,self.name,module))
    #
    def dependencies(self,project):
        "Build the dependencies from the list of use"
        dependencies = set()
        dirs = set()
        for module in self.use_modules:
            a = project.modules.get(module)
            if a:
                if a.dir == self.dir:
                    #Only dependencies of modules with the same directories
                    dependencies.add(a)
                else:
                    #Add dir which contains a module needed for the compilation of self
                    dirs.add(a.dir)
        return (dependencies,dirs)
    #
    def update_declaration(self,variable):
        "Add a variable into the declaration and also in each children"
        if not self.Declaration.dict_vars.has_key(variable.name):
            self.Declaration.dict_vars[variable.name] = variable
        for child in self.children:
            if isinstance(child,Module):
                child.update_declaration(variable)


class Routine(Module):
    "Class to handle subroutines or function"
    #Ampersand (continuation line)
    re_amp = re.compile('[ ]*&[ ]*')
    #Detect call (not after an if or ;)
    re_call = re.compile('(^[ \t]*call[ ]*)(\w+)', re.MULTILINE+re.IGNORECASE)
    #Detect optional arguments
    re_optional = re.compile("optional.*::", re.IGNORECASE)
    #
    def __init__(self,name=None,parent=None,implicit=None,message=None):
        "Initialisation"
        Module.__init__(self,name,parent=parent,message=message)
        self.nature = "subroutine"
        self.called = set()
        self.calling = set()
        #No interface defined for this routine (call self.interface to do that)
        self.has_interface = False
        self.optional = False
        #Signature (header+arguments) of the routine
        self.signature = ""
        #Use implicit (if inside module)
        self.implicit = implicit
        #If not a part of a generic routine
        self.nogeneric = True
        #Private
        self.private = False
        #Public
        self.public = False
        self.from_implicit = False
        self.is_argument = False
        #Gestion of the cache
        #1 - Time stamp (done in Module class)
        #2 - The information of the cache is updated
        self.cached_updated = False
        #3 - called routines
        self.cached_called = None
        #4 - arguments
        self.cached_arguments = None
        #5 - calling
        self.cached_calling = None
    #
    def add_use_interface(self,project):
        "Add 'use self.module'"
        modules = set()
        #Add routines and functions
        functions = list()
        else_modules = set()
        for called in self.called:
            try:
                routine = project.routines[called]
            except KeyError:
                if called in intrinsic_routines:
                    #Next routine
                    continue
                else:
                    #Next routine
                    self.message.warning_no_interface(self.dir,self.file,self.name,called)
                    continue
            modules.add(routine.module)
            if isinstance(routine,Function):
                functions.append(routine)
        #Remove main (for the program "optic" which has some subroutines !!!)
        modules.discard("main")
        #Add all modules for routine
        self.use_modules.update(modules)
        #Build a list in order to sort
        list_modules = list(modules)
        list_modules.sort()
        list_else_modules = list(else_modules)
        list_else_modules.sort()
        #Modify the interfaces
        self.Use.add_use_interface(list_modules,list_else_modules)
        #Add declarations of functions
        if functions:
            #Sort functions (to have unique ordering)
            functions.sort(Structure.cmp)
            self.Declaration.add_functions(functions)
    #
    def analyze(self,line,iter_code,project):
        "Iterate the iterator up to the end of the routine"
        Code.analyze(self)
        #We are already inside a routine
        #We have the first line of header
        self.Header = Header_Routine(parent=self)
        #One line is not enough: Add line if continuation line
        self.Header.analyze(line,iter_code)
        #Add the routine to the project
        project.add_routine(self)
        #Add the cached information
        self.cache_use(project)
        #Analyze the use statement
        self.Use = Use(parent=self)
        (final_comments,line) = self.Use.analyze(iter_code)
        #Add the modules
        self.use_modules.update(self.Use.modules)
        #Test if implicit statement
        self.Implicit = Implicit(parent=self)
        #More than one line are analyzed
        line = self.Implicit.analyze(line,iter_code,comments=final_comments)
        #Analyze the declarations
        self.Declaration = Declaration(parent=self,comments=final_comments)
        line = self.Declaration.analyze(line,iter_code,self.Implicit.dict,project)
        #Analyze execution statements
        self.Execution = Execution(parent=self)
        self.Execution.analyze(line,iter_code,project)
    #
    def analyze_execution(self,project,edition=False):
        "Analyze the execution statements of the routine"
        #Analyze the execution statements
        self.Execution.analyze_execution(self.Declaration.dict_vars,self.Implicit.dict,project,edition=edition)
    #
    def beautify(self):
        "Beautify the code of the routines (experimental)"
        #Remove tabulations
        self.code = self.code.replace("\t",indent)
        #Format the comment "!Arguments"
        self.code = re_arguments.sub(self.code,"\n\n" + comment_arguments)
        #Format the comment "!Local variables"
        self.code = re_local_variables.sub(self.code,"\n\n" + comment_local_variables)
    #
    def build_declaration(self):
        "Build declaration for routine declared in an interface"
        if isinstance(self.parent,Fortran_Interface):
            self.ancestry(0)
            return ("interface",self.__str__())
    #
    def cache_save(self):
         "Save some information for the cache"
         return { "timestamp": self.timestamp,
                  "called": self.called,
                  "arguments": self.arguments,
                  "calling": self.calling }
    #
    def cache_use(self,project):
        "Use the cached information"
        #Test time stamp to save time
        if project.cached.has_key(self.lower):
            cached = project.cached[self.lower]
            if cached["timestamp"] == self.timestamp:
                self.cached_updated = True
                self.cached_called = cached["called"]
                self.cached_arguments = cached["arguments"]
                self.cached_calling = cached["calling"]
            else:
                #The cache is not updated. Nevertheless, we use some information
                self.cached_updated = False
                self.cached_called = cached["called"]
                self.cached_calling = cached["calling"]
    #
    def display_information(self):
        "Display information about the routine"
        message = "%s %s\n" % (self.nature,self.name)
        if self.public:
            message += 3*" " + "is public\n"
        elif self.private:
            message += 3*" " + "is private\n"
        if self.nogeneric:
            message += 3*" " + "is a generic routine\n"
        self.message.write(message,verbose=-10)
    #
    def get_arguments(self,project):
        "Build the list of arguments of a routine"
        #Arguments
        if self.cached_arguments != None:
            self.arguments = self.cached_arguments
        else:
            self.arguments = self.Declaration.arguments(self.Header.arguments,self.Implicit.dict,project)
    #
    def get_dependencies(self):
        "Give a list of variables on which depends the function"
        self.dependencies = set()
        return self.dependencies
    #
    def graph_children(self,dict_graph,routines=list(),recursive=False):
        "Add children in the graph"
        calling = self.name
        if calling in dict_graph:
            #Already done: To avoid Recursive call, return
            return
        else:
            dict_graph[calling] = set()
        for name in self.called:
            try:
                struct = routines[name]
                dict_graph[calling].add(struct.name)
                if recursive:
                    struct.graph_children(dict_graph,routines=routines,recursive=True)
            except KeyError:
                dict_graph[calling].add(name)
    #
    def graph_parents(self,dict_graph,routines=list(),recursive=False):
        "Add parents in the graph"
        for name in self.calling:
            #Add the parents routines
            try:
                struct = routines[name]
                calling = struct.name
                if calling not in dict_graph:
                    dict_graph[calling] = set()
                called = self.name
                if called in dict_graph[calling]:
                    #Already done: To avoid recursive call, continue
                    continue
                dict_graph[calling].add(called)
                if recursive:
                    struct.graph_parents(dict_graph,routines=routines,recursive=True)
            except KeyError:
                pass
    #
    def has_type(self):
        "A subroutine has always a type"
        return True
    #
    def interface(self):
        "Build the interface"
        if self.arguments == "":
            #Nothing to be detected
            self.has_interface = False
            return ""
        #Header (de-indent inside a continuation line)
        header = self.re_amp.sub('&'+indent*2,self.Header.code.lstrip())
        #Check if an optional argument is present
        if self.re_optional.search(self.arguments):
            self.optional = True
        code = self.arguments
        #Store the interface without interface...end interface
        self.signature = header \
                        + code \
                        + "end %s %s\n" % (self.nature,self.name)
        #For the interface, no implicit none
        code = indent_code("interface\n" \
                           + header \
                           + code \
                           + "end %s %s\n" % (self.nature,self.name) \
                           + "end interface\n")
        #Store the interface for this routine
        self.code_interface = code
        self.has_interface = True
        return code
    #
    def order(self):
        "Determine the order"
        return 10000
    #
    def set_called_routines(self,project):
        "Determine routines and functions which are called by the routine (only this part is different for a generic routine"
        #Define the called routines
        list = self.re_call.findall(self.Execution.code)
        for tt in list:
            called = tt[1].lower()
            self.called.add(called)
        #Add functions
        if project.re_functions:
            #We detect function in comments!!
            list = project.re_functions.findall(self.Execution.code)
            functions = set()
            for tt in list:
                called = tt.lower()
                if called not in intrinsic_routines and called != self.name:
                    #Add the functions
                    functions.add(called)
            self.called.update(functions)
    #
    def set_children_parents(self,project):
        "Create the list of children and parents of the routines"
        #Save time if called is already done
        if self.cached_calling != None:
            #Even it is not updated, we use this information
            self.calling = self.cached_calling
        if self.cached_updated:
            self.called = self.cached_called
        else:
            #In this case, read the code and search called routines (time consuming)
            self.set_called_routines(project)
            #Now change the calling routines of other routines.
            if self.cached_called:
                called_to_remove = self.cached_called.difference(self.called)
                #Remove self into calling routines which are not already called
                for called in called_to_remove:
                    if project.routines.has_key(called):
                        #Remove itself from the call
                        project.routines[called].calling.discard(self.lower)
                #New called routine to add
                called_to_add = self.called.difference(self.cached_called)
            else:
                called_to_add = self.called
            #Add self into calling routines for called
            for called in called_to_add:
                if project.routines.has_key(called):
                    project.routines[called].calling.add(self.lower)
                #elif called in intrinsic_routines:
                #    #Do nothing
                #    pass
                #else:
                #    self.message.warning_no_reference(self.dir,self.file,self.name,called)
        #Verbose message
        #text = "(%s:" % self.name
        #for called in self.called:
        #    text += "{%s}" % called
        #self.message.write(text+")")
        self.message.write("<%s>" % self.name)


class Generic_Routine(Routine):
    "Class to handle generic routine"
    def __init__(self,name,file):
        "Initialisation"
        Routine.__init__(self,name,parent=None,message=file.message)
        #A generic routine is not the child of a file but it is related to a file.
        #Initialisation done
        self.generic_file = file
        self.file = file.name
        self.dir = file.dir
        self.module = file.module
    #
    def add_use_interface(self,project):
        "Do nothing"
        return
    #
    def analyze_execution(self,project,edition=False):
        "Do nothing"
        return
    #
    def get_arguments(self,project):
        "Do nothing"
        self.arguments = ""
        return
    #
    def interface(self):
        "Do nothing"
        return ""
    #
    def set_called_routines(self,project):
        "Determine routines and functions which are called by the routine: Use all routines from the other routines"
        for routine in self.generic_file.children:
            if isinstance(routine,Routine):
                self.called.update(routine.called)


class Function(Routine):
    "Class to handle Function"
    def __init__(self,name=None,parent=None,implicit=None,message=None):
        "Initialisation"
        Routine.__init__(self,name,parent,implicit,message)
        self.nature = "function"
        #Type of the function
        self.type = None
    #
    def analyze(self,line,iter_code,project):
        "Analyze the function"
        Code.analyze(self)
        self.Header = Header_Function(parent=self)
        self.Header.analyze(line,iter_code)
        #Add in the dictionary of all routines
        project.add_routine(self)
        #Add the cached information
        self.cache_use(project)
        #Analyze the use statement
        self.Use = Use(parent=self)
        (final_comments,line) = self.Use.analyze(iter_code)
        #Add the modules
        self.use_modules.update(self.Use.modules)
        #Test if implicit statement
        self.Implicit = Implicit(parent=self)
        #More than one line are analyzed
        line = self.Implicit.analyze(line,iter_code,comments=final_comments)
        #Analyze the declarations
        self.Declaration = Declaration(parent=self,comments=final_comments)
        line = self.Declaration.analyze(line,iter_code,self.Implicit.dict,project)
        #Analyze the execution
        self.Execution = Execution(parent=self)
        self.Execution.analyze(line,iter_code,project)
    #
    def get_arguments(self,project):
        "Build the list of arguments of a function"
        #Arguments (use the method of Routine class)
        Routine.get_arguments(self,project)
        #Check if the function has already a type
        if not self.type:
            #Check if dict_vars is not equal to None (case if the cache is used)
            if self.Declaration.dict_vars == None:
                #Dictionary of variables
                self.Declaration.analyze_variables(self.Implicit.dict,project)
            name = self.func_name.lower()
            if not self.Declaration.dict_vars.has_key(name):
                if self.Implicit.dict:
                    #There is an implicit dictionary
                    Variable.type_from_implicit(self,self.Implicit.dict)
                else:
                    self.message.fatal("\n%s/%s: Check if the arguments are well formatted:\n" \
                        % (self.dir,self.file) \
                        + "Type not declared for the function '%s'\n" % self.name)
            else:
                #We are looking for self.name in dict_vars
                self.type = self.Declaration.dict_vars[name].type
        #Add the type of the function
        name = self.name.lower()
        if not self.Declaration.dict_vars.has_key(name):
            self.Declaration.dict_vars[name] = Variable(self.name,parent=self,type=self.type)
        #If in a module or a subroutine, add the type of the function in parent.dict_vars
        if isinstance(self.parent,Module):
            #Need to update dict_vars for dict_vars of the module and other functions
            self.parent.update_declaration(self.Declaration.dict_vars[name])
    #
    def has_type(self):
        return self.type != None


class Program(Routine):
    "Fortran program (as a routine)"
    def interface(self):
        "No interface for program"
        self.code_interface = ""
        self.has_interface = False
        return ""

class Header_Routine(Code):
    "Header of a routine"
    #Detect the beginning of a block (program, subroutine, function)
    re_startblock = re.compile('(?P<header>[ \t]*' \
        + '(?P<type>(program|((pure|recursive)[ ]*)*subroutine|(([^!\n]*?)function)))' \
        + '[ ]*(?P<name>\w+)[ ]*(?P<arguments>[(][^)]*[)])?[^\n]*)\n', \
           re.MULTILINE+re.IGNORECASE)
    #In arguments, remove some characters
    re_inarg = re.compile('[ \t\n&()]+')
    #
    #Analyze the code
    def analyze(self,line,iter_code):
        "Analyze the header"
        Code.analyze(self)
        self.add_code(line)
        #Remove comments at the end of the line only for analysis
        code_an = self.re_comment.sub('',line)
        while self.re_continuation.search(line) or self.re_comment_match.match(line):
            line = iter_code.next()
            self.add_code(line)
            #Remove comments
            line = self.re_comment.sub('',line)
            code_an += line
        #Analyze the header
        search = self.re_startblock.match(code_an).groupdict()
        self.parent.type = search['type']
        self.parent.name = search['name']
        args = search['arguments']
        if args:
            #Remove blank characters
            args = args.replace(' ','')
            if args == '()':
                args = None
        if args:
            self.arguments = self.re_inarg.sub('',args)
            self.arguments = self.arguments.split(',')
        else:
            self.arguments = list()
        self.parent.lower = self.parent.name.lower()
        self.message.write("<%s" % self.parent.name)


class Header_Function(Header_Routine):
    "Header of a function"
    #Use to determine the variable in result statement for a function
    re_result = re.compile('result[ ]*[(](?P<result>[^)]+)[)]')
    def analyze(self,line,iter_code):
        "Analyze the header"
        Code.analyze(self)
        Header_Routine.analyze(self,line,iter_code)
        #Determine the type of the function
        if (len(self.parent.type) > 9):
            #'xxx function'
            self.parent.type = self.parent.type[:-9].strip()
        else:
            #Determine inside the code
            #First check if result and use the corresponding name
            search = self.re_result.search(self.code)
            if search:
                self.parent.func_name = search.groupdict()['result']
            else:
                self.parent.func_name = self.parent.name
            #Add at the beginning of the list of arguments
            self.arguments.insert(0,self.parent.func_name)


class Use(Code):
    "Use statement"
    #Empty preprocessing directives
    re_empty_preproc = re.compile(\
                "^[ ]*#[ \t]*if[^\n]+\n([ ]*#[ \t]*else[ ]*\n)?[ ]*#[ \t]*endif[^\n]*\n?",re.MULTILINE)
    #Multiple (more than 3 \n): bug or unavailable
    re_multi_n = re.compile("[\n]{3,}")
    #Use statement
    re_use = re.compile('^[ \t]*use[ ]*(?P<name>\w+)', re.MULTILINE+re.IGNORECASE)
    re_use_prefix = re.compile('^[ \t]*use[ ]*'+prefix+'.*?\n', re.MULTILINE+re.IGNORECASE)
    #
    #Add use
    def add_use_interface(self,modules,else_modules):
        "Add 'use interfaces_xxx'"
        #Remove section created by abilint (replace by \n if a use after)
        self.code = re_section_abilint.sub('\n',self.code)
        if re_section_abilint_start.search(self.code):
            self.message.error("Routine %s: alone header of an abilint section" % self.parent.name)
            #Remove the beginning of a section created by abilint (if separated)
            self.code = re_section_abilint_start.sub('',self.code)
        if re_section_abilint_end.search(self.code):
            self.message.error("Routine %s: alone footer of an abilint section" % self.parent.name)
            #Remove the end of a section created by abilint (if separated)
            self.code = re_section_abilint_end.sub('\n',self.code)
        #Remove line with use prefix
        self.code = self.re_use_prefix.sub('',self.code)
        #Remove some replaced use
        self.code = re_use_replaced.sub('',self.code)
        #Remove empty preprocessing directives
        self.code = self.re_empty_preproc.sub('',self.code)
        #Add preprocessing commands, comments and implicit none
        if modules or else_modules:
            text_use = ""
            #Add modules
            for module in modules:
                #Special and ugly
                if "contract" in module:
                    text_use += "#if defined DEBUG_CONTRACT\n"
                if "cuda" in module:
                    text_use += "#if defined HAVE_GPU_CUDA\n"
                if prefix in module:
                    if module == self.parent.module and self.parent.has_interface:
                        #We except the given subroutine
                        text_use += indent + "use %s, except_this_one => %s\n" \
                               % (self.parent.module,self.parent.name)
                    else:
                        text_use += indent + "use %s\n" % module
                if "contract" in module:
                    text_use += "#endif\n"
                if "cuda" in module:
                    text_use += "#endif\n"
            if else_modules:
                text_use += "#else\n"
                for name in else_modules:
                    text_use += indent + "use %s\n" % name
            if text_use != "":
                text_use = "\n\n" + use_before + text_use + use_after + "\n"
            #Add text_use inside use statements
            self.code += text_use
        else:
            #Be sure to have 2 \n
            self.code += (2-self.code[-2:].count("\n"))*"\n"
        #
        # MG Add CPP variable with the name of the procedure.
        #self.code += "\n #undef ABI_FUNC \n #def ABI_FUNC " + str(self.parent.name) 
        #print self.code + "\n #undef ABI_FUNC \n #define ABI_FUNC " + str(self.parent.name) 
        #
        #Remove multiple \n
        self.code = self.re_multi_n.sub('\n\n',self.code)
    #
    #Analyze the code
    def analyze(self,iter_code):
        """Analyze use statements 
           (special treatment for preprocessing commands which needs to be clarify)
           only is ignored."""
        Code.analyze(self)
        #List of used modules
        self.modules = set()
        comments = ""
        inside_if = False
        for line in iter_code:
            res = self.re_use.match(line)
            if res:
                #Add comments and the line
                self.add_code(comments+line)
                #Add this module
                self.modules.add(res.groupdict()["name"].lower())
                while self.re_continuation.search(line):
                    #print line,res.groupdict()["name"].lower()
                    line = iter_code.next()
                    self.add_code(line)
                comments=""
            elif self.re_allcomment.match(line):
                if "#if" in line:
                    inside_if = True
                elif "#endif" in line:
                    inside_if = False
                comments += line
            else:
                #We have finished
                final_comments = ""
                if inside_if:
                    #We remove the beginning of the preprocessing directives
                    a_comments = comments.split("\n")
                    #Last element is not a line
                    a_comments.pop()
                    a_comments.reverse()
                    comments = ""
                    l_final = len(a_comments)
                    for cline in a_comments:
                        if "#if" in cline:
                            l_final = 2
                        if l_final > 0:
                            final_comments = cline + "\n" + final_comments
                        else:
                            comments = cline + "\n" + comments
                        l_final -= 1
                if comments:
                    self.add_code(comments)
                return (final_comments,line)


class Implicit(Code):
    "Class for the statement 'implicit'"
    re_imp_stat = re.compile("(?P<type>.*)[ ]*[(](?P<range>[^)]+)[)]")
    default_dict = dict()
    #Default in Fortran
    for i in range(97,105):
        #real
        default_dict[chr(i)] = "real"
    for i in range(105,111):
        #integer
        default_dict[chr(i)] = "integer"
    for i in range(111,123):
        #real
        default_dict[chr(i)] = "real"
    #
    #Detect implicit statements
    def analyze(self,line,iter_code,comments=""):
        "Test if implicit statements are present"
        Code.analyze(self)
        if "implicit" in line.lower():
            if comments:
                #Fatal error
                self.message.fatal("%s/%s[%s]: " % (self.dir,self.file,self.parent.name) \
                    + "Preprocessing block not closed for use statements\n" )
            #Check if many implicit statements
            while True:
                search = self.re_continuation.search(line)
                self.code += line
                line = iter_code.next()
                if not ("implicit" in line.lower() or search):
                    #The last line is not an implicit statement or the previous one has not a continuation sign
                    break
            #Build a dictionary
            self.implicit()
            return line
        elif self.parent.implicit != None:
            self.dict = self.parent.implicit
            return line
        else:
            #Error (no implicit statement) and build a dictionary
            self.implicit()
            self.message.error("No implicit statement in %s (%s/%s)" \
                % (self.parent.name,self.dir,self.file))
            return line
    #    
    #Build a dictionary
    def implicit(self):
        "Give the type of variables given by implicit statement"
        lines = self.code.lower()
        if "none" in lines:
            self.dict = dict()
            return
        #Default implicit in Fortran (warning: copy the defaut dict)
        self.dict = self.default_dict.copy()
        if self.code == "":
            #Finished
            return
        #Split in lines and remove blank lines
        lines = lines.split("\n")
        lines.remove('')
        for line in lines:
            #We analyze the implicit statements
            line = line.replace('implicit','').strip()
            search = self.re_imp_stat.match(line)
            if search:
                search = self.re_imp_stat.match(line).groupdict()
            else:
                self.message.fatal("%s/%s[%s]: " % (self.dir,self.file,self.parent.name) \
                    + "%s" % line \
                    + "Analysis error of implicit statement\n" )
            type_imp = search["type"].strip()
            #a-h,...
            table = search["range"].split(",")
            for letters in table:
                (a1,a2) = letters.strip().split('-')
                a1 = ord(a1)
                a2 = ord(a2)+1
                for i in range(a1,a2):
                    self.dict[chr(i)] = type_imp
        #All lines are analyzed
        return


class Declaration(Code):
    "Declaration statement"
    #Detection of character as style f77 + f95
    re_character_f77 = re.compile('[ ]*character[ ]*([*]?[(][^)]+[)]|[*][0-9]+|)',re.IGNORECASE)
    #Detect declarations
    re_declaration = re.compile('^[ \t]*' \
        + '(allocatable|character|common|complex|data|dimension|double|end[ ]+type|equivalence|external|'\
        + 'integer|intrinsic|logical|optional|parameter|private|public|real|save|type)', re.IGNORECASE)
    #Detect "type " or "type," or "type::"
    re_def_type = re.compile('^[ \t]*type[ ]*(?![(])',re.IGNORECASE)
    #Detect digits only (1.2d0 or 1.3e-4 etc.)
    re_digits = re.compile("^\d+[.]?\d*[de]?[+-]?[\d]*$")
    #Detect digits at the beginning for 0.0_dp
    re_digits_start = re.compile("^\d+[.]?\d*[de]?[+-]?[\d]*_")
    #Remove elements of type
    re_element = re.compile("%\w+")
    #Detect a group after =
    re_equal = re.compile('[ ]*=.*')
    #Detect interface
    re_interface_start = re.compile('^[ \t]*interface',re.IGNORECASE)
    re_interface_end   = re.compile('^[ \t]*end[ ]+interface',re.IGNORECASE)
    #Multiple \n
    re_mn = re.compile('\n+',re.MULTILINE)
    #No ampersand
    re_noamp = re.compile("[ ]*[&]+[ ]*")
    #Detect only declarations without include, data, external and save
    re_only_declaration = re.compile('^[ \t]*' \
        + '(allocatable|character|complex|dimension|double|end[ ]+type|integer|interface|logical|parameter|private|public|real|type)',\
        re.IGNORECASE)
    #For character(len=xxx)
    re_character = re.compile('character[(](len[ =]+)?(?P<len>[^)]+)[)]',re.IGNORECASE)
    #For complex(kind=dp) or complex(kind(dp))
    re_complex = re.compile('complex[(](kind[=(])?(?P<kind>[^)]+)[)]+',re.IGNORECASE)
    #Detect !Local
    re_local = re.compile("!Local",re.IGNORECASE)
    #For parameter(xxx)
    re_parameter = re.compile('parameter[ ]*[(](?P<liste>[^)]+)[)]',re.IGNORECASE)
    #For real(kind=dp) or real(kind(dp))
    re_real = re.compile('real[(](kind[=(])?(?P<kind>[^)]+)[)]+',re.IGNORECASE)
    #For type(xxx)
    re_type = re.compile('type[(](?P<type>[^)]+)[)]',re.IGNORECASE)
    #Detect variable var(.*) and also xxx=
    #re_var = re.compile('(\w+([(][^)]+[)])?)')
    re_var = re.compile('(\w+[ ]*((=[^,]+)|([(][^)]+[)])?))')
    #
    def __init__(self,name=None,parent=None,comments=""):
        "Initialisation"
        Code.__init__(self,name,parent)
        self.code = comments
        #Dictionary of all variables
        self.dict_vars = None
        #Include directives
        self.includes = list()
        #True if all variables are private by default
        self.private = False
        #True if all variables are public by default for module
        if hasattr(self.parent,"is_module"):
            self.public = self.parent.is_module
        else:
            self.public = False
    #
    #Add functions in the declaration
    def add_functions(self,functions):
        "Add declaration of function with preprocessing"
        #Remove section created by abilint
        self.code = re_section_abilint.sub('\n',self.code)
        if re_section_abilint_start.search(self.code):
            self.message.error("Routine %s: alone header of an abilint section" % self.parent.name)
            #Remove the beginning of a section created by abilint (if separated)
            self.code = re_section_abilint_start.sub('',self.code)
        if re_section_abilint_end.search(self.code):
            self.message.error("Routine %s: alone footer of an abilint section" % self.parent.name)
            #Remove the end of a section created by abilint (if separated)
            self.code = re_section_abilint_end.sub('',self.code)
        #We remove the declaration of the function
        for struct in functions:
            #Match at least 4 characters before :: which are not !, ' and $ in order to keep preprocessing directives
            #Warning: does not detect '& xxxxx' or in a continuation line
            self.code = re.sub("(?i)([^$!']{4,}[ ]*::.*?),?[ ]*%s(?=[, \n])" % struct.name, r"\1",self.code)
        #We remove two commas ',,'
        self.code = re.sub(",,",",",self.code)
        #We remove a comma behind "::"
        self.code = re.sub( "::[ ]*,", ":: ",self.code)
        #We remove a line ending by "::" with a comma or not behind
        self.code = re.sub( ".*::[ ]*[,]*[ ]*\n", "",self.code)
        #i# text_functions = ""
        #i# for struct in functions:
        #i#     text_functions += indent + "%s :: %s\n" % (struct.type,struct.name)
        #i# #Add at the end of declaration
        #i# if re_end_variables.search(self.code):
        #i#     #Add before "!******" only once (count=1)
        #i#     self.code = re_end_variables.sub("\n" + notuse_before + text_functions \
        #i#             + notuse_after + end_variables,self.code,count=1)
        #i# else:
        #i#     self.code += notuse_before + text_functions + notuse_after + "\n"
    #
    #Analyze the code
    def analyze(self,line,iter_code,dict_implicit,project):
        "Analyze the declaration statements"
        Code.analyze(self)
        if not line:
            line = iter_code.next()
        while True:
            if self.re_declaration.search(line):
                #We can not avoid assignation as 'reali=2'
                self.add_code(line)
                #Is it a continuation line or a comment line inside continuations lines?
                while self.re_continuation.search(line) or self.re_comment_match.match(line):
                    line = iter_code.next()
                    self.add_code(line)
            elif self.re_include.match(line):
                #Include directive (not add as child)
                struct = Include(project,line=line,file=self.file,dir=self.dir)
                self.includes.append(struct)
                self.add_code(line)
            elif self.re_comment_match.match(line):
                self.add_code(line)
            elif self.re_interface_start.match(line):
                #We enter in an interface block
                self.add_code(line)
                self.analyze_interface(iter_code)
            else:
                #We have finished
                #Test if not implicit
                if "implicit" in line:
                    self.message.fatal("\n%s/%s: [%s]\n--> Implicit statement '%s' inside declaration of variables" \
                            % (self.parent.dir,self.parent.file,self.parent.name,line[:-1]) \
                            + "\nPossible reason: An include statement before implicit statement\n")
                else:
                    break
            try:
                line = iter_code.next()
            except StopIteration:
                self.message.fatal("\n%s/%s: [%s]\n--> End of the code inside declaration!!\n" \
                        % (self.parent.dir,self.parent.file,self.parent.name))
        #Get the last line
        return line
    #
    def analyze_interface(self,iter_code):
        "Iterate inside an interface"
        in_interface = 1
        while in_interface > 0:
            try:
                line = iter_code.next()
            except StopIteration:
                self.message.fatal("\n%s\n--> No detected end of the interface!\n" % line\
                                + "Analysis Error in %s/%s\n" % (self.parent.dir,self.parent.file))
            self.add_code(line)
            line_lower = line.lower()
            if "interface" in line_lower:
                if self.re_interface_end.match(line):
                    in_interface -= 1
                elif self.re_interface_start.match(line):
                    in_interface += 1
    #
    #Analyze all variables xxxx :: var -> dict[var] == xxx
    def analyze_variables(self,dict_implicit,project):
        "Analyze all variables xxxx :: var -> dict[var] == xxx"
        #First, we build a long line without continuation, comments and interface
        code = self.code
        #We add the include files at the end of declaration
        for struct in self.includes:
            code += struct.struct.code
        iter_code = iter(code.splitlines(-1))
        code = list()
        self.dict_vars = dict()
        for line in iter_code:
            if self.re_only_declaration.search(line):
                #Definition of a type?
                if self.re_def_type.match(line):
                    #Create a type
                    ftype = Fortran_Type(parent=self.parent,private=self.private,public=self.public)
                    #Analyze the type (and go to "end type")
                    ftype.analyze(line,iter_code,dict_implicit,project)
                    #Add the type and the name in the dictionary of variables
                    self.dict_vars[ftype.lower] = ftype
                    #Next line
                    continue
                elif self.re_interface_start.match(line):
                    #Definition of an interface : Create an interface
                    interface = Fortran_Interface(parent=self.parent)
                    #Analyze the interface (and go to "end interface" and add the subroutines or the interface)
                    self.dict_vars.update(interface.analyze(line,iter_code,dict_implicit,project))
                    continue
                #Remove '\n'
                line = self.re_comment.sub('',line).lstrip()
                code.append(line[:-1])
                while self.re_continuation.search(line):
                    null_line = True
                    #Detect if comment lines inside continuation lines
                    while null_line:
                        line = iter_code.next()
                        line = self.re_comment.sub('',line).lstrip()
                        null_line = (len(line) == 0)
                        code[-1] += line[:-1]
                #Remove ampersand
                code[-1] = self.re_noamp.sub('',code[-1])
        #First we assume to have ::
        iter_code = iter(code)
        #Build: declaration (decl) :: list of variables
        for line in iter_code:
            # xxx :: yyy
            i = line.find("::")
            if i > -1:
                #Easy
                decl = line[:i]
                liste = line[i+2:]
            else:
                line_lower = line.lower()
                #Keyword variables (not very robust)
                liste = line.strip().split()
                decl = liste[0]
                decl_lower = decl.lower()
                if decl_lower == "private":
                    #All variables are private (do not work if comment after private)
                    self.private = True
                    continue
                elif decl_lower == "public":
                    #All variables are public
                    self.public = False
                    continue
                elif decl_lower == "double":
                    #Add 'precision' or 'complex'
                    decl += " " + liste[1]
                    liste = reduce(lambda x,y: x + " " + y, liste[2:])
                elif decl_lower[0:4] == "real" or \
                        decl_lower[0:7] == "complex" or \
                        decl_lower[0:7] == "integer" or \
                        decl_lower[0:7] == "logical" or \
                        decl_lower == "dimension" or decl_lower == "allocatable":
                    #In this case: decl = liste[0] and liste = the rest
                    liste = reduce(lambda x,y: x + " " + y, liste[1:])
                elif "character" in decl_lower:
                    #Detect character without ::
                    decl = self.re_character_f77.match(line).group()
                    liste = line.replace(decl,"")
                    decl = decl.replace(" ","")
                elif decl_lower[0:9] == "parameter":
                    res = self.re_parameter.search(line.strip())
                    if res:
                        decl="parameter"
                        liste=res.groupdict()["liste"]
                    else:
                        self.message.fatal("\n%s\n--> Parameter statement not correct!\n" % line\
                            + "Analysis Error in %s/%s\n" % (self.parent.dir,self.parent.file))
                elif decl_lower[0:4] == "type":
                    #Detect the used type without ::
                    liste = liste[1]
                else:
                    self.message.fatal("\n%s\n--> Strange declaration!\n" % line\
                        + "Analysis Error in %s/%s\n" % (self.parent.dir,self.parent.file))
            #Declaration -- lists of variables
            decl0 = decl.lower().strip()
            liste = liste.strip()
            for (name,var) in split_variables(liste):
                if "character" in decl0 and '*' in name:
                    #The name is not correct
                    (var,l) = name.split('*')
                    (name,l) = name.split('*')
                    decl = decl0.replace("character","character(len="+l+")")
                else:
                    decl = decl0
                if self.dict_vars.has_key(name):
                    self.dict_vars[name].update(decl,truename=var)
                else:
                    self.dict_vars[name] = Variable(name,parent=self,decl=decl,truename=var,\
                                                    public=self.public,private=self.private)
        #Check if all variables have a type
        if not dict_implicit:
            #Implicit none
            for (name,variable) in self.dict_vars.items():
                if not variable.has_type():
                    self.message.fatal( "[%s/%s:%s] Variable {%s} has no type\n" \
                        % (self.parent.dir,self.parent.file,self.parent.name,name))
        else:
            #There is an implicit dictionary
            for (name,variable) in self.dict_vars.items():
                if variable.type == "":
                    #Add type of the implicit dictionary
                    variable.type_from_implicit(dict_implicit)
        #Put private variables in a private directory and remove it. Get arguments for routines
        self.dict_vars_private = dict()
        for (name,variable) in self.dict_vars.items():
            if variable.private:
                self.dict_vars_private[name] = variable
                del(self.dict_vars[name])
            #if isinstance(variable,Routine):
                #variable.get_arguments(project)
    #
    #Give the declaration of arguments
    def arguments(self,arguments,dict_implicit,project):
        "Give declarations of arguments"
        if not self.dict_vars:
            return ""
        #List of declarations
        declarations = list()
        #Set for parameters in argument (ex: selected_kind)
        parameters = dict()
        #Set for modules in parameters
        modules = set()
        arguments_lower = map(lambda x: x.lower(),arguments)
        unfound_modules = []
        for argument in arguments:
            names = set()
            argument_lower = argument.lower()
            if self.dict_vars.has_key(argument_lower):
                arg = self.dict_vars[argument_lower]
                #Is an argument
                arg.is_argument = True
            else:
                #We use the implicit dictionary
                arg = Variable(argument_lower,parent=self.parent,truename=argument,is_argument=True)
                if not arg.type_from_implicit(dict_implicit):
                    arg.display_information()
                    text = "\nArguments of the routine '%s':" % self.parent.name
                    for arg in arguments_lower:
                        text += " '%s'" % arg
                    text += "\n==Code"+"="*50+">\n%s" % self.code + "="*56+"<\n"
                    self.message.fatal( \
                        text + "[%s/%s:%s] Argument '%s' is not declared\n" \
                        % (self.parent.dir,self.parent.file,self.parent.name,argument))
                #Add in the dictionary of variables
                self.dict_vars[argument_lower] = arg
            #Determine the dependencies and the order
            names = arg.get_dependencies()
            order = arg.order()
            #Add the declaration of the argument with order
            declarations.append((order,arg))
            #For each variables which depend arguments, we add information
            for name in names:
                #Find where some information about 'name' is stored
                has_name = False
                if name in arguments_lower:
                    #We have already
                    has_name = True
                elif self.dict_vars.has_key(name):
                    #Declared in the variable for the subroutine
                    parameters[name] = self.dict_vars[name]
                    has_name = True
                elif name in intrinsic_functions:
                    #This is an intrinsic functions
                    has_name = True
                else:
                    #Check if some information is not stored in a module
                    unfound_module = None
                    for module in self.parent.use_modules:
                        if project.has_module(module):
                            dict_vars = project.modules[module].Declaration.dict_vars
                            if dict_vars.has_key(name):
                                modules.add(module)
                                has_name = True
                                break
                        else:
                            #A module is missing
                            unfound_module=module
                            if unfound_module not in unfound_modules:
                                self.message.warning("[%s/%s:%s] The module '%s' is missing!" \
                                    % (self.parent.dir,self.parent.file,self.parent.name,module))
                                unfound_modules.append(module)
                if not has_name:
                    texte = "%s/%s:%s\n" % (self.parent.dir,self.parent.file,self.parent.name)
                    texte += "--> Arguments of the routine '%s':" % self.parent.name
                    for argu in arguments:
                        texte += " '%s'" % argu
                    if len(self.parent.use_modules) > 0:
                        texte += "\n    Used modules:"
                        for module in self.parent.use_modules:
                            texte += " %s" % module
                    #If inside a module: don't care for interface
                    if self.parent.module:
                        message = self.message.error
                    else:
                        message = self.message.fatal
                    if unfound_module:
                        texte += "\n Unfound modules:"
                        for unfound_module in unfound_modules:
                            texte += " %s" % unfound_module
                        message("%s\n   " % texte \
                                + " This modules are not found and"  \
                                + " the argument '%s' depends on '%s' which could be in these modules." \
                                % (arg.name,name))
                    else:
                        message("%s\n[%s/%s:%s]:" % \
                                           (texte,self.parent.dir,self.parent.file,self.parent.name) \
                                + " The argument '%s' depends on '%s'" % (arg.name,name) \
                                + " which is not found in the declaration even in a module." )
        #Build the declaration of arguments
        decl_args = ""
        #Sort arguments (integer, scalar, arrays)
        declarations.sort()
        #Add modules
        for module in modules:
            decl_args += "use %s\n" % module
        #Add implicit none
        decl_args += "implicit none\n"
        #Add parameters
        decl_args += build_declaration(parameters)
        for (a,arg) in declarations:
            if isinstance(arg,Routine):
                #Special case
                line = "interface\n%send interface\n" % arg
            else:
                line = "%s :: %s\n" % arg.build_declaration()
                if len(line) > 100:
                    #Split the line
                    line = line[:60] + line[60:].replace(',',', &\n&'+9*' ',1)
            decl_args += line
        return decl_args


class Execution(Code):
    "The execution of the routine or a function"
    #Comment at the end of a line
    re_comment_line = re.compile('[!].*\n')
    #Character string
    re_string = re.compile('''['"][^'"]*['"]''')
    #Alphanumeric
    re_word = re.compile('[a-z]\w*')
    #List of reserved keywords (special case for call, end, format)
    reserved_keywords = [ "case", "close", "continue", "cycle",
                         "default", "do", "else", "elseif", "exit", "goto","go","to",
                         "if", "nullify",
                         "print", "return", "select", "stop", "then", "where", "while" ]
    #Add intrinsic functions
    reserved_keywords.extend(intrinsic_functions)
    #List of reserved keywords with variables
    dict_keywords = { "allocate": [ "stat" ],
                      "backspace": [ "unit", "iostat", "err" ],
                      "close": [ "unit", "iostat", "err", "status"],
                      "deallocate": [ "stat" ],
                      "endfile": [ "unit", "iostat", "err" ],
                      "inquire": [ "unit", "exist", "file", "opened", "number", "named", "access", "sequential", "form", "formatted",
                                   "recl", "nextrec", "blank", "position", "action", "read", "write", "readwrite", "delim", "pad"],
                      "maxval": [ "mask" ],
                      "minval": [ "mask" ],
                      "open": [ "access", "action", "blank", "delim", "file", "form", "iostat", "pad", "position", "recl", "status", "unit"],
                      "pack": [ "mask" ],
                      "read": [ "fmt", "iostat", "rec", "unit"],
                      "rewind": [ "unit", "iostat", "err" ],
                      "sum": [ "mask" ],
                      "write": [ "advance", "fmt", "iostat", "unit"]}
    #
    def analyze(self,line,iter_code,project):
        "Iterate the iterator up to the end of the routine."
        Code.analyze(self)
        self.code = ""
        while True:
            line_lower = line.lstrip()[:8].lower()
            if "contains" == line_lower:
                if self.re_contains.match(line):
                    self.message.write(":")
                    Module.analyze_contains(self.parent,line,iter_code,project)
                    #The end is detected
                    self.message.write(">")
                    return
            elif "end" == line_lower[:3]:
                if self.re_sub_end.match(line):
                    #Detect the end of the subroutine: We have finished
                    self.code += line
                    self.message.write(">")
                    return
            self.code += line
            try:
                line = iter_code.next()
            except StopIteration:
                self.message.fatal("\n%s\n--> No detected end of the routine!\n" % line\
                        + "Analysis Error in %s/%s\n" % (self.parent.dir,self.parent.file))
    #
    def analyze_execution(self,dict_vars,dict_implicit,project,edition=False):
        "Analyze the execution statements (initialization, undeclared variables)"
        if dict_vars == None:
            self.message.fatal("No dictionary of variables\n" \
                    + "Analysis Error in %s/%s\n" % (self.parent.dir,self.parent.file))
        #We have a dictionary
        code = split_code(self.code)
        variables = set()
        #Check if all variables are declared
        for line in code:
            iter_line = iter(line)
            reserved = list()
            #Use to detect inside a call and remove xxx=
            in_call = False
            for word in iter_line:
                if self.re_word.match(word):
                    if word in self.dict_keywords.keys():
                        #Is a keyword with variables
                        reserved = self.dict_keywords[word]
                    elif word in self.reserved_keywords:
                        #Is a keyword
                        pass
                    elif word in reserved:
                        #Reserved keywords
                        pass
                    elif word == "call":
                        #We iterate (not consider the next word)
                        in_call = True
                        iter_line.next()
                    elif word == "end" or word == "enddo" or word == "endif":
                        #We iterate at the end of the line
                        #if "end do", we detect if there is a block name
                        l=len(line)
                        if (word == "endif" and l == 2) or (l == 3 and line[1] == "if") or \
                           (word == "enddo" and l == 2) or (l == 3 and line[1] == "do"):
                            variables.remove(line[-1])
                        break
                    elif word == "format":
                        #We skip the line
                        break
                    else:
                        if in_call:
                            #Check if xxx= (should be better to have also reserved words)
                            try:
                                next = iter_line.next()
                            except:
                                print in_call
                                print code
                                print line
                                sys.exit(1)
                            if next == "=":
                                #Do not add
                                continue
                        #Add the variable
                        variables.add(word)
        declared_variables = set(dict_vars.keys())
        #Add the variables from modules
        for module in self.parent.use_modules:
            if project.has_module(module):
                declared_variables.update(project.modules[module].Declaration.dict_vars.keys())
            else:
                #A module is missing
                self.message.warning("[%s/%s:%s] The module '%s' is missing!" \
                    % (self.parent.dir,self.parent.file,self.parent.name,module))
        #Treat undeclared variables
        undeclared = variables.difference(declared_variables)
        #Add undeclared variables inside dict_vars
        for name in undeclared:
            dict_vars[name] = Variable(name,parent=self)
            dict_vars[name].type_from_implicit(dict_implicit)
        #Add arguments or variables not typed inside dict_vars
        undeclared= dict()
        for var in dict_vars.values():
            if var.from_implicit:
                undeclared[var.lower] = var
        if undeclared:
            #Build the declaration of all variables (to build only undeclared, use undeclared)
            declaration = build_declaration(dict_vars)
            if declaration:
                self.message.error("%s/%s: [%s]\n--> Undeclared variables\n%s" \
                    % (self.dir,self.file,self.parent.name,build_declaration(undeclared)))
                #Edition using vim
                if edition:
                    for var in dict_vars.values():
                        var.display_information()
                    #Edit inside vim the whole file: put in the register "a the declarations
                    vim_message = "Register v = %s" % declaration \
                            + "Register a = !Arguments\n" \
                            + "Register l = !Local variables\n" \
                            + "Edit the routine %s (use reg v, reg a or reg l)\n" % self.parent.name \
                            + "Create a file EXIT to stop all editions"
                    vim_commands = """let @a ="%s!Arguments\\n"\n""" % indent \
                            + """let @l = "%s!Local variables\\n"\n""" % indent \
                            + """let @v ="%s"\n""" % declaration.replace("\n","\\n") \
                            + """let @/ ="%s"\n""" % self.parent.name \
                            + """echo "%s"\n""" % vim_message.replace("\n","\\n")
                    open("temp.vim","w").write(vim_commands)
                    #subprocess.call(["vim", "%s/%s" % (self.parent.dir,self.parent.file)])
                    os.system("vim -S temp.vim %s/%s" % (self.parent.dir,self.parent.file))
                    os.remove("temp.vim")
                    if os.path.exists("EXIT"):
                        self.message.write("\nThe file 'EXIT' has been edited.\nSTOP EDITION and abilint.\n",verbose=-10)
                        sys.exit(1)
            self.message.done()


#Class for Fortran variables
class Variable:
    "Class for Fortran variables (type, dimension, ...)"
    #Built-in types
    builtin_types = [ "character", "complex", "double complex", "double precision", "integer",
                     "logical" ]
    #Reserved words for declaration
    reserved = [ "allocatable", "character", "complex", "dimension", "doublecomplex", "doubleprecision",
                 "integer","intent",
                 "in", "out", "inout", "kind", "len", "optional", "real", "target", "pointer",
                 "logical", "type", "parameter", "max", "min", "modulo", "merge",
                 "private", "public", "len_trim", "interface", "function" ]
    #No letter (add . to have 1.d0 and % for type)
    re_noletter = re.compile("[^a-zA-Z0-9_.%]+")
    #Detect digits only (1.2d0 or 1.3e-4 etc.)
    re_digits = re.compile("^\d+[.]?\d*[de]?[+-]?[\d]*$")
    #
    def __init__(self,name,parent=None,message=None,decl="",type="",truename="",\
             public=False,private=False,is_argument=False):
        "Initialization: decl=declaration (lower case), truename = exact name of the variable"
        self.name = name
        self.lower = name.lower()
        #Proper attributes
        self.type=type
        self.is_argument = is_argument
        self.allocatable = False
        self.dimension = ''
        self.intent = ''
        self.external = False
        self.optional = False
        self.parameter = False
        self.pointer = False
        self.private = private
        self.public = public
        self.save = False
        self.target = False
        self.value = None
        self.parent = parent
        self.dependencies = None
        if self.type:
            self.from_implicit = None
        else:
            self.from_implicit = False
        if parent:
            self.message = self.parent.message
        elif message:
            self.message = message
        else:
            sys.stderr.write("Structure %s (%s) has no message class!\n" \
                    % (str(self.__class__),self.name))
            sys.exit(1)
        #Truename is a, a(5) or A=5
        if truename:
            self.analyze_truename(truename)
        else:
            self.truename = name
        #Analyze the declaration
        self.analyze_decl(decl)
    #
    def analyze_decl(self,decl):
        "Analyze the declaration of a variable"
        for (head,all) in split_variables(decl):
            if head == '':
                return
            elif head[0:9] == "character" or head[0:7] == "complex" or \
                    head == "doublecomplex" or head == "doubleprecision" or \
                    head[0:7] == "integer" or head == "logical" or \
                    head[0:4] == "real" or head == "type":
                self.type = all
            elif head == "allocatable":
                self.allocatable = True
            elif head == "dimension":
                dim = all[9:]
                if dim:
                    #To avoid to remove dimension information as dimension :: a(5)
                    self.dimension = dim
            elif head == "intent":
                self.intent = head[7:-2]
            elif head == "external":
                self.external = True
            elif head == "optional":
                self.optional = True
            elif head == "parameter":
                self.parameter = True
            elif head == "pointer":
                self.pointer = True
            elif head == "private":
                self.private = True
                self.public = False
            elif head == "public":
                self.public = True
                self.private = False
            elif head == "save":
                self.save = True
            elif head == "target":
                self.target = True
            else:
                #Fatal error
                if self.parent:
                    message = "%s/%s: " % (self.parent.dir,self.parent.file)
                else:
                    message = ""
                message += "[%s]\n--> Declaration error: %s [%s][%s]" % (self.name,decl,head,all)
                self.message.fatal(message)
    #
    def analyze_truename(self,truename):
        "Analyze the variable a(5) or a=5"
        if "=" in truename:
            a = truename.split("=")
            self.truename = a[0]
            self.value = a[1]
        elif "(" in truename:
            a = truename.split("(")
            self.truename = a[0]
            self.dimension = truename.replace(a[0],'')
        else:
            self.truename = truename
    #
    def get_dependencies(self):
        "Give a list of variables on which depends the variable in dimension, kind or type"
        self.dependencies = set()
        if self.type:
            self.search_name(self.type)
        elif self.dimension:
            self.search_name(self.dimension)
        return self.dependencies
    #
    def build_declaration(self):
        "Build the declaration of the variable"
        decl = self.type
        #var = self.truename+self.dimension
        var = self.truename
        if self.parameter:
            decl += ", parameter"
            var += "=" + self.value
        if self.intent:
            decl += ", " + self.intent
        if self.public:
            decl += ", public"
        elif self.private:
            decl += ", private"
        if self.allocatable:
            decl += ", allocatable"
        if self.optional:
            decl += ", optional"
        elif self.external:
            decl += ", external"
        if self.pointer:
            decl += ", pointer"
        if self.save:
            decl += ", save"
        if self.target:
            decl += ", target"
        if self.dimension:
            decl += ", dimension%s" % self.dimension
        return (decl,var)
    #
    def display_information(self):
        "Display information about the variable"
        if self.is_argument:
            message = "Argument "
        else:
            message = "Variable "
        message += "%s [%s]\n" % (self.name,self.truename) \
                + 3*" " + "type: %s" % self.type
        if self.from_implicit:
            message += " (from the implicit statement)"
        message += "\n"
        if self.dimension:
            message += 3*" " + "dimension %s\n" % self.dimension
        if self.parameter:
            message += 3*" " + "is a parameter with a value of '%s'\n" % self.value
        if self.intent:
            message += 3*" " + "intent(%s)\n" % self.intent
        if self.public:
            message += 3*" " + "is public\n"
        if self.private:
            message += 3*" " + "is private\n"
        if self.allocatable:
            message += 3*" " + "is allocatable\n"
        if self.optional:
            message += 3*" " + "is optional\n"
        elif self.external:
            message += 3*" " + "is external\n"
        if self.pointer:
            message += 3*" " + "is a pointer\n"
        if self.save:
            message += 3*" " + "is saved\n"
        if self.target:
            message += 3*" " + "is a target\n"
        if self.dependencies:
            liste = ""
            for name in self.dependencies:
                liste += name + " "
            message += 3*" " + "depends on %s\n" % liste
        self.message.write(message,verbose=-10)
    #
    def has_type(self):
        "True if has a type (self.type) of public or private"
        return (self.type != None) or self.public or self.private
    #
    def order(self):
        "Determine the order"
        self.get_dependencies()
        if self.dimension:
            if self.dependencies:
                order = 100
            else:
                order = 10
        elif "integer" in self.type:
            #Insert first the integer
            order = 0
        else:
            #Scalar after the integer
            order = 1
        return order
    #
    def search_name(self,string):
        "Each non-alphanumeric character is converted into blank, split and remove keywords."
        #If we find a type, do not select built-in types
        if string in self.builtin_types:
            return
        liste = self.re_noletter.sub(" ",string).split()
        for name in liste:
            if name in self.reserved:
                pass
            elif not self.re_digits.match(name):
                #Remove %xxx in order to keep the name of the type
                i = name.find("%")
                if i > -1:
                    name = name[:i]
                self.dependencies.add(name)
    #
    def type_from_implicit(self,dict_implicit):
        "Determine the type from the implicit dictionary"
        if dict_implicit != {}:
            self.type = dict_implicit[self.name[0]]
            self.from_implicit = True
        return self.type
    #
    def update(self,decl,truename=None):
        "Add information about the variable"
        if truename:
            self.analyze_truename(truename)
        self.analyze_decl(decl)


#Class to handle the type of Fortran datastructure
class Fortran_Type(Declaration):
    "Class to handle the Fortran type of datastructure"
    #Name of the type ("type xxxx", "type :: xxxx" or "type, yyy :: xxxx"
    re_name_type = re.compile('^[ \t]*type[ ]*(|.*::)[ ]*(?P<name>\w+)',re.IGNORECASE)
    #Detect "end type"
    re_end_type = re.compile('^[ \t]*end[ ]+type',re.IGNORECASE)
    #Detect only declarations without include, data, external and save
    re_type_declaration = re.compile('^[ \t]*' \
        + '(character|complex|dimension|double|integer|logical|parameter|private|public|real|type)', re.IGNORECASE)
    #
    def __init__(self,parent,name=None,public=False,private=False):
        "Special initialisation: we do not add as child to the parent"
        self.name = name
        if name:
            self.lower = self.name.lower()
        self.parent = parent
        #From parent, pick up the implicit dictionary
        self.implicit = parent.implicit
        self.message = parent.message
        self.file = parent.file
        self.dir = parent.dir
        self.code = ""
        self.code_backup = ""
        #Below are defined for compatibility with the class Routine
        self.module = self.parent.module
        self.timestamp = None
        self.children = list()
        self.cached = dict()
        #List of used modules
        self.use_modules = set()
        #List of includes directives
        self.includes = list()
        self.is_analyzed = False
        #Type (for dict_vars)
        self.type = "type"
        self.dict_vars = dict()
        #Private
        self.private = private
        #Public
        self.public = public
        self.from_implicit = False
        self.is_argument = False
    #
    def analyze(self,line,iter_code,dict_implicit,project):
        "Analyze the type statements"
        Code.analyze(self)
        #Determine the name
        self.name = self.re_name_type.match(line).groupdict()["name"]
        self.lower = self.name.lower()
        #Add a new fortran type and check if it is already defined
        if project.ftypes.has_key(self.lower):
            struct = project.ftypes[self.lower]
            self.message.error("%s/%s: [%s]\n--> This fortran type is already defined in %s/%s!" \
                    % (self.dir,self.file,self.name,struct.dir,struct.file))
        else:
            project.ftypes[self.lower] = self
        #Header (type <name>)
        self.header = line
        while True:
            try:
                line = iter_code.next()
            except StopIteration:
                self.message.fatal("%s/%s: [%s]\n--> Incorrect end of a fortran type!!\n'%s'\n" \
                        % (self.dir,self.file,self.name,line[:-1]))
            if self.re_type_declaration.search(line):
                self.add_code(line)
                #Is it a continuation line or a comment line inside continuations lines?
                while self.re_continuation.search(line) or self.re_comment_match.match(line):
                    line = iter_code.next()
                    self.add_code(line)
            elif self.re_comment_match.match(line):
                self.add_code(line)
            elif self.re_end_type.search(line):
                #We have finished
                self.end = line
                break
            else:
                self.message.fatal("%s/%s: [%s]\n--> Incorrect line in a fortran type declaration!\n'%s'\n" \
                        % (self.dir,self.file,self.name,line[:-1]))
        #We analyze all variables
        self.analyze_variables(dict_implicit,project)
    #
    def build_declaration(self):
        "Need to code"
        decl = "type"
        if self.public:
            decl += ", public"
        elif self.private:
            decl += ", private"
        decl += ":: %s" % self.name
        decl += self.code
        decl += "end type %s" % self.name
        return ("type", decl)
    #
    def display_information(self):
        "Display information about the fortran type"
        message = "Type %s\n" % self.name
        if self.public:
            message += 3*" " + "is public\n"
        elif self.private:
            message += 3*" " + "is private\n"
        message += self.code
        self.message.write(message,verbose=-10)
    #
    def get_dependencies(self):
        "Give a list of variables on which depends the variable in dimension, kind or type"
        self.dependencies = set()
        for child in self.children:
            self.dependencies.update(child.get_dependencies())
        return self.dependencies
    #
    def has_type(self):
        return True
    #
    def order(self):
        "Determine the order"
        order = 100
        return order
    #
    def update(self,decl,truename=None):
        "Add information as public or private"
        if "private" in decl:
            self.private = True
            self.public = False
        elif "public" in decl:
            self.public = True
            self.private = False


#Class to handle the fortran interfaces inside declarations (interface itself is not a child but
#defined functions in the interface do
class Fortran_Interface(Fortran_Type):
    "Class to handle a fortran interface"
    #Detect the name of the interface
    re_interface_name = re.compile('^[ \t]*interface[ ]+(?P<name>\w+)?',re.IGNORECASE)
    re_module_procedure = re.compile('^[ \t]*module[ ]+procedure[ ]+(?P<name>\w+)?',re.IGNORECASE)
    #
    def __init__(self,parent,name=None):
        "__init__ function as Fortran_Type, only the type differs"
        Fortran_Type.__init__(self,parent,name)
        self.type = "interface"
        #Private
        self.private = False
        #Public
        self.public = False
        self.from_implicit = False
    #
    def add_routine(self,routine):
        "Do nothing"
        pass
    #
    def analyze(self,line,iter_code,dict_implicit,project):
        "Analyze an interface"
        Code.analyze(self)
        #Determine the name of the interface
        reg = self.re_interface_name.match(line)
        if reg:
            #This is a generic interface
            self.name = reg.groupdict()['name']
            self.lower = self.name.lower()
            self.generic = True
        else:
            #Declare only interfaces for the compiler
            self.generic = False
            #Defined children
            self.children = list()
        in_interface = 1
        #Detect subroutines, functions or 'module procedure'
        for line in iter_code:
            line_lower = self.re_allcomment.sub('',line.lower())
            if self.re_comment_match.match(line):
                continue
            elif "subroutine" in line_lower:
                struct = Routine(parent=self,implicit=dict_implicit)
                struct.analyze(line,iter_code,self)
                #Add the modules used in the interface
                self.use_modules.update(struct.use_modules)
            elif "function" in line_lower:
                struct = Function(parent=self,implicit=dict_implicit)
                struct.analyze(line,iter_code,self)
                #Add the modules used in the module to the used modules by the routine
                self.use_modules.update(struct.use_modules)
            elif "module" in line_lower and "procedure" in line_lower:
                name = self.re_module_procedure.match(line_lower).groupdict()['name']
                struct = Routine(parent=self,name=name,implicit=dict_implicit)
                #Add as child (generic routine)
                self.children.append(struct)
                #If continuation mark, iterate
                while self.re_continuation.search(line):
                    line = iter_code.next()
            elif self.re_interface_end.match(line_lower):
                in_interface -= 1
                self.add_code(line)
                #We have finished
                break
            else:
                self.message.fatal("\n'%s'\n--> No detected routine (Interface analysis)!\n" % line \
                        + "This part of code can not be parsed as Fortran file:\n" \
                                   + "Analysis Error in %s/%s:<%s>\n" % (self.dir,self.file,self.name))
        if in_interface:
            self.message.fatal("\n%s\n--> No detected end of the interface!\n" % line\
                               + "Analysis Error in %s/%s:<%s>\n" % (self.dir,self.file,self.name))
        #Give modules in parent
        self.parent.use_modules.update(self.use_modules)
        #Return a dictionary
        if self.generic:
            return { self.lower: self}
        else:
            children = dict()
            for child in self.children:
                children[child.lower] = child
            return children
    #
    def update(self,decl,truename=None):
        "We do not need to analyze generic interface"
        if self.generic:
            pass
        else:
            #Not implemented: error
                self.message.fatal("\n--> Not implemented for function interface!\n" \
                        + "This part of code can not be parsed as Fortran file:\n" \
                        + "Analysis Error in %s/%s\n" % (self.dir,self.file))


class Include(Code):
    "Class to handle include statement"
    #Extract file of the include file
    re_include_file = re.compile("""^[ ]*include[ ]+['"](?P<file>[\w./_-]+)['"]""",re.IGNORECASE)
    def __init__(self,project,parent=None,line=None,file=None,dir=None,):
        "Initialisation"
        Code.__init__(self,parent=parent,message=project.message)
        if not parent:
            if file:
                self.file = file
            if dir:
                self.dir = dir
        if line:
            try:
                self.includefile = self.re_include_file.match(line).groupdict()['file']
            except AttributeError:
                self.message.fatal("\n'%s'\n--> Wrong include statement!\n" % line[:-1] \
                            + "This part of code can not be parsed as Fortran file:\n" \
                            + "Analysis Error in %s/%s\n" % (self.dir,self.file))
            self.message.write("(include file %s)" % self.includefile)
            self.code = line
            self.struct = None
            if project:
                #Find the include file in the project
                self.struct = project.find_file(self.includefile)
                if not self.struct:
                    if self.includefile in project.given_include:
                        #We use own version of these include files as 'mpif.h'
                        self.message.error("Use own version of '%s'" % self.includefile)
                        file = project.add_file('.',self.includefile,create=True,read_only=True,File_Class=File_F90)
                        file.add_code(project.given_include[self.includefile])
                        self.struct = file
                    else:
                        self.message.error("[%s/%s]: The include file '%s' does not exist in the project!" % \
                                           (self.dir,self.file,self.includefile))
                        #We create one with no code
                        file = project.add_file('.',self.includefile,create=True,read_only=True,File_Class=File_F90)
                        self.struct = file


class FatalError(Exception):
    "Error raised when Fatal message"
    def __init__(self,value):
       self.value = value
    def __str__(self):
        text = ""
        for line in self.value.splitlines(1):
            text += "!!   %s" % line
        return text

class Message:
    "Message management"
    def __init__(self,verbose=0,logfile="message.log"):
        "Initialisation"
        #List of errors
        self.terrors = list()
        #List of warnings
        self.twarnings = list()
        #Text of the Fatal error (None if not fatal error)
        self.tfatal = None
        #Message to avoid fatal error
        self.mfatal = ""
        #Raise an exception (if nofatal==True) for a fatal error instead of exiting
        self.nofatal = False
        #Set of not referenced routines
        self.noreference = dict()
        #Verbosity
        self.verbose = verbose
        #Open the log file
        self.logfile = logfile
        self.log = open(logfile,"w")
        #Inform if the last character was a '\n'
        self.ret = True
    #
    def bug(self,text):
        "Display a message about a possible bug and stops"
        tt = "Bug: %s\n" % text
        sys.stderr.write(tt)
        sys.exit(1)
    #
    def close(self):
        "Close the log file"
        self.log.close()
    #
    def done(self):
        "Write done"
        self.write("done.\n",verbose=-10,flush=True)
    #
    def error(self,text):
        "Display an error message"
        tt = "Error[%d]: %s\n" % (len(self.terrors)+1,text)
        self.terrors.append(tt)
        if len(self.terrors) == 1 or not self.ret:
            sys.stdout.write("\n"+tt)
        else:
            sys.stdout.write(tt)
        self.ret = True
    #
    def flush(self):
        "Flush the buffer"
        sys.stdout.flush()
    #
    def fatal(self,text):
        "Display an error message and stops"
        #Flush before writing into stderr
        sys.stdout.flush()
        sys.stderr.write("\nFatal Error: %s\n" % text)
        sys.stderr.write(self.mfatal)
        self.tfatal = text
        if self.nofatal:
            #Raise an exception
            raise FatalError,text
        else:
            sys.exit(1)
    #
    def final(self):
        "Final statistics: Write in a log file and resume."
        #Write at the end of the log file all warning and errors
        for tt in self.twarnings:
            self.log.write(tt)
        #Not referenced routines
        items = self.noreference.items()
        items.sort()
        for item in items:
            if item == 1:
                self.log.write("The routine %s is not referenced once.\n")
            else:
                self.log.write("The routine %s is not referenced %d times.\n" % item)
        #Errors
        for tt in self.terrors:
            self.log.write(tt)
        if self.tfatal:
            self.log.write(self.tfatal)
            self.log.write("\n" + "*"*20 + " Fatal Error " + "*"*20 + "\n\n")
            sys.stdout.write("\n" + "*"*20 + " Fatal Error " + "*"*20 + "\n\n")
        else:
            self.log.write("%d warnings.\n" % len(self.twarnings))
            self.log.write("%d not referenced routines.\n" % len(self.noreference))
            if self.terrors <= 1:
                self.log.write("%d error.\n" % (len(self.terrors)))
            else:
                self.log.write("%d errors.\n" % (len(self.terrors)))
            sys.stdout.write("%d warnings (%d not referenced routines) -- %d errors.\n" \
                    % (len(self.twarnings),len(self.noreference),len(self.terrors)))
        sys.stdout.write("See the file '%s' for warnings.\n" % self.logfile)
    #
    def section(self,text):
        "Display a message and flush"
        if not self.ret:
            self.write('\n')
        self.write(text,verbose=-10,flush=True)
        self.ret = (text[-1] == '\n')
    #
    def warning(self,text,verbose=0):
        "Display a warning message"
        tt = "Warning[%d]: %s\n" % (len(self.twarnings)+1,text)
        self.twarnings.append(tt)
        if verbose <= self.verbose:
            #Display in stdout
            if not self.ret:
                tt = '\n'+tt
            sys.stdout.write(tt)
            self.ret = True
    #
    def warning_no_interface(self,dir,file,routine,called):
        "Display a warning message concerning not referenced routine"
        if called in self.noreference.keys():
            self.noreference[called] += 1
        else:
            self.noreference[called] = 1
        self.warning("[%s/%s:%s] No interface for the routine '%s'" \
               % (dir,file,routine,called))
    #
    def warning_no_reference(self,dir,file,routine,called):
        "Display a warning message concerning not referenced routine (for interface)"
        if called in self.noreference.keys():
            self.noreference[called] += 1
        else:
            self.noreference[called] = 1
        self.warning("[%s/%s:%s] No reference for the routine '%s'" \
               % (dir,file,routine,called))
    #
    def write(self,text,verbose=0,flush=False):
        "Display a message"
        if not text:
            return
        if verbose <= self.verbose:
            sys.stdout.write(text)
            #Copy to the log file
            self.log.write(text)
            self.ret = (text[-1] == '\n')
        if flush:
            #Flush stdout
            sys.stdout.flush()

#-----------------------
#End of class defintions
#-----------------------


#Header of the interface module.
head_interface = \
"""!!****m* BigDFT/%(name)s
!! NAME
!! %(name)s
!!
!! FUNCTION
!! This module contains the interfaces of the routines
!! %(description)s
!!
!! COPYRIGHT
!! Copyright (C) 2008-2011 BigDFT group
!! This file is distributed under the terms of the
!! GNU General Public License, see ~bigdft/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! NOTES
!! THIS FILE IS GENERATED AUTOMATICALLY BY abilint.
!! To do that: config/scripts/abilint . .
!! 
%(warning)s!!
!! SOURCE

module %(name)s

 implicit none

"""


#Header of the dependencies files (abinit.dep).
head_dependencies = \
"""#Dependencies of the directory %(dir)s
#
#COPYRIGHT
#Copyright (C) 2008-2011 BigDFT group
#This file is distributed under the terms of the
#GNU General Public License, see ~abinit/COPYING
#or http://www.gnu.org/copyleft/gpl.txt .
#
#THIS FILE IS GENERATED AUTOMATICALLY BY abilint.
#To do that: config/scripts/abilint --dependencies . .

"""


#Ersatz of 'mpif.h' file
mpif_file = """
integer, parameter :: mpi_comm_self = 1, mpi_comm_null = 2
integer, parameter :: mpi_comm_world = 0, mpi_max = 1, mpi_min = 2, mpi_sum = 3
integer, parameter :: mpi_character = 5
integer, parameter :: mpi_integer = 7, mpi_integer8=11
integer, parameter :: mpi_real = 13, mpi_complex=18
integer, parameter :: mpi_double_precision = 17, mpi_double_complex=22
!Seems to have no value ?
integer, parameter :: mpi_in_place = 999
"""
#Ersatz of 'fftw3.f' file
fftw3_file = """
integer, parameter :: fftw_estimate = 64
integer, parameter :: fftw_forward = -1
"""

#Include files given by abilint
bigdft_include = { "mpif.h": mpif_file,
                   "fftw3.f": fftw3_file}


#Exclude files (*/*.f90)
bigdft_exclude = [ "unused", "tools", "PSolver/base.f90", \
                   "others.F90",\
                   "lib/lbfgs.f90", "wavelib/i-o-etsf.f90" ]


#Files to generate generic routines
generic_routines = []

#File class per pattern
bigdft_File_Class = [("*.F90",File_F90), ("*.f90",File_F90), ("*.f",File_F77)]
#Add generic routine
for file in generic_routines:
    bigdft_File_Class.insert(0,(file,File_F90_Generic))


#Special modification inside the code
special_modif =  {}


#Define the hierarchy of directories to check the correctness of calls
def rank_dir(dir):
    "Define the hierarchy of directories in the project"
    return 0


#Routines excluded in the graph
graph_excluded = [ "MPI_INIT", "MPI_COMM_RANK", "MPI_COMM_SIZE" ]
#Add intrinsic routines
graph_excluded.extend(intrinsic_routines)

#----------------
#Start the script
#----------------
if __name__ == "__main__":
    #Check arguments
    try:
        optlist, args = getopt.getopt(sys.argv[1:],\
            "bcdug:hiv",
                ["beautify","nocache","dependencies","dump_dtset=","graph=","help","lint","verbose"])
    except getopt.error:
        sys.stderr.write("Error in arguments\n")
        usage(error=1)
    #Help message
    verbose = -10
    beautify = False
    graph = False
    lint = False
    dump_dtset = False
    dependencies = False
    edition = False
    nocache = False
    for opt,arg in optlist:
        if   opt == "-b" or opt == "--beautify":
            beautify = True
            edition = True
            #Complete analysis
            lint = True
        elif opt == "-d" or opt == "--dependencies":
            dependencies = True
        elif opt == "-u" or opt == "--dump_dtset":
            dump_dtset = True
            dtset_arg = arg
        elif opt == "-g" or opt == "--graph":
            graph = True
            graph_arg = arg
        elif opt == "-h" or opt == "--help":
            usage(error=0)
        elif opt == "-i" or opt == "--lint":
            lint = True
        elif opt == "-c" or opt == "--nocache":
            nocache = True
        elif opt == "-v" or opt == "--verbose":
            verbose = 10
    #Two arguments are required
    if len(args) != 2:
        sys.stderr.write('You have to specify the bigdft source directory and the destination.\n')
        sys.stderr.write('    Ex: abilint . .\n')
        usage(error=1)
    #Define OLD and NEW directories
    OLD = args[0]
    NEW = args[1]
    #Create the project and read all files
    bigdft = Project(OLD,name="BigDFT",\
                     pat_dir=["src","src/*","libABINIT","libABINIT/src","libABINIT/src/*"],\
                     pat_file=["*.F90","*.f90", "*.inc"],\
                     logfile="abilint.log",\
                     exclude=bigdft_exclude,given_include=bigdft_include,\
                     File_Class=bigdft_File_Class,\
                     style_comment="doxygen")
    bigdft.message.write("(%d directories, %d files)" \
            % (len(bigdft.dirs),len(bigdft.files)),verbose=-10)
    bigdft.message.done()
    #We check for a cache file
    if NEW == OLD:
        bigdft.cache_load(NEW)
    #Analyze the project.
    #bigdft.analyze_all(exclude="interfaces_")
    bigdft.analyze_all()
    #Set the called routines
    bigdft.set_children_parents()
    if lint:
        #Analyze the interdependencies between directories
        bigdft.analyze_directories()
        #Unused routines
        #bigdft.unused_routines()
        #Analyze all the comments in order to detect the robodoc structure
        bigdft.analyze_comments(edition=edition)
        #Analyze the code (body)
        bigdft.analyze_execution(edition=edition)
    if NEW == OLD:
        #Backup before changing to optimize the minimal number of copies
        bigdft.backup()
    #Special modifications (to be done once).
    bigdft.special(special_modif)
    if beautify:
        #Edit and improve the appearance of the code 
        bigdft.beautify()
    if graph:
        #Build in the file bigdft.dot some graph
        bigdft.graph(graph_arg,graph_excluded=graph_excluded)
    if edition:
        bigdft.message.write("The files have been edited: copy nothing.\n",verbose=-10)
    else:
        #We copy everything
        bigdft.copy_all(NEW,only_if_changed=(NEW == OLD))
    #We copy a cache file
    bigdft.cache_save(NEW)
    #Display some statistics
    bigdft.statistics()
    #Close the log file
    bigdft.message.close()

