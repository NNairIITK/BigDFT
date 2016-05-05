#!/usr/bin/env python
# -*- coding: us-ascii -*-
#----------------------------------------------------------------------------

import math
import sys
import os,copy,optparse
import cProfile

#path of the script
path=os.path.dirname(sys.argv[0])
#yaml_folder = '/local/gigi/binaries/ifort-OMP-OCL-CUDA-gC/PyYAML-3.10/lib'
#if yaml_folder not in sys.path:
#  sys.path.insert(0,'/local/gigi/binaries/ifort-OMP-OCL-CUDA-gC/PyYAML-3.10/lib')

import yaml

EVAL = "eval"
SETUP = "let"
INITIALIZATION = "globals"

PRE_POST = [EVAL, SETUP, INITIALIZATION]

ENERGY = "BigDFT.energy"
FERMI_LEVEL= "__FERMI_LEVEL__"
NUMBER_OF_ATOMS = 'BigDFT.nat'

#Builtin pathes to define the search paths
BUILTIN={ENERGY: [["Last Iteration", "FKS"],["Last Iteration", "EKS"]],
         FERMI_LEVEL: [["Ground State Optimization", -1, "Fermi Energy"]],
         NUMBER_OF_ATOMS: [ ['Atomic System Properties','Number of atoms']]}

# print out a python dictionary in yaml syntax
def dict_dump(dict):
  sys.stdout.write(yaml.dump(dict,default_flow_style=False,explicit_start=True))

#function which removes from a set of lines the yaml_fields contained in the to_remove list
def clean_logfile(logfile_lines,to_remove):
  line_rev=logfile_lines #list of the lines of the logfile
  #loop in the reversed from (such as to parse by blocks)
  extra_lines=20 #internal variable to be customized
  line_rev.reverse()
  #clean the log
  cleaned_logfile=[]
  removed=[]
  #for line in line_rev: #line_iter:
  while len(line_rev) >0:
    line=line_rev.pop()
    to_print=line
    #check if the line contains interesting information
    for remove_it in to_remove :
      stream_list=[]
      #line without comments
      valid_line=line.split('#')[0]
      spaces='nospace'
      #control that the string between the key and the semicolon is only spaces
      if remove_it in valid_line and ":" in valid_line:
        #print "here",remove_it,remove_it in valid_line and ":" in valid_line,valid_line
        starting_point=valid_line.find(remove_it)
        tmp_buf=valid_line[:starting_point]
        #find the closest comma to the staring point, if exists
        tmp_buf=tmp_buf[::-1]
        starting_comma=tmp_buf.find(',')
        if starting_comma <0: st=0
        tmp_buf=tmp_buf[st:]
        tmp_buf=tmp_buf[::-1]
        tmp_buf=tmp_buf.strip(' ')
        #print "there",tmp_buf,'starting',starting_point,len(tmp_buf)
        valid_line= valid_line[starting_point+len(remove_it):]
        spaces= valid_line[1:valid_line.find(':')]
        #if remove_it+':' in line.split('#')[0]:
      if len(spaces.strip(' ')) == 0 and len(tmp_buf)==0: #this means that the key has been found
         #creates a new Yaml document starting from the line
         #treat the rest of the line following the key to be removed
         header=''.join(line.split(':')[1:])
         header=header.rstrip()+'\n'
         #eliminate the anchor
         header=header.lstrip(' ')
         header=header.lstrip('*') 
         if len(header) > 0 :
            stream_list.append(header)
         #part to be printed, updated
         to_print = line.split(':')[0] + ": <folded> \n"

         #then check when the mapping will end:
         while True:
            #create a stream with extra_lines block           
            for i in range(0,min(extra_lines,len(line_rev))):
               stream_list.append(line_rev.pop())
            #create a stream to be parsed
            stream=''.join(stream_list)
            #then parse the stream until the last valid position has been found
            try:
              for i in yaml.parse(stream,Loader=yaml.CLoader):
                endpos=i.end_mark.index
            except Exception, e:
              #  print 'error',str(e),stream
              #convert back the valid stream into a list
              #if needed the stream can be loaded into a document
              item_list=stream[:endpos].split('\n')
              #if lengths are different there is no need to add lines
              if len(item_list) != len(stream_list):
                #last line might be shorter, therefore treat it separately
                last_line=item_list.pop()
                #purge the stream
                for item in item_list:
                  stream_list.remove(item+'\n')
                #extract the remaining line which should be compared with the last one
                strip_size=len(last_line.rstrip())
                if strip_size > 0:
                  first_line=stream_list.pop(0)[strip_size:]
                  if '*' in first_line or '&' in first_line:
                    first_line='' #eliminate anchors
                else:
                  first_line=''
                #then put the rest in the line to be treated
                to_print.rstrip('\n')
                to_print += first_line+'\n'
                # the item has been found
                break
         stream_list.reverse()
         #put back the unused part in the document
         line_rev.extend(stream_list)
         # mark that the key has been removed
         if (remove_it not in removed):
           removed.append(remove_it)
           print 'removed: ',remove_it
  # then print out the line 
    cleaned_logfile.append(to_print)

  # check that everything has been removed, at least once
  if (set(removed) != set(to_remove)):
    print 'WARNING, not all the requested items have been removed!'
    print 'To_remove : ',to_remove
    print 'removed   : ',removed
    print 'Difference: ',list(set(to_remove) - set(removed) )
  return cleaned_logfile


# this is a tentative function written to extract information from the runs
def document_analysis(doc,to_extract):
  analysis=[]
  for quantity in to_extract:
    #print quantity
    #follow the levels indicated to find the quantity
    value=doc
    for key in quantity:
      #as soon as there is a problem the quantity is null
      try:
        value=value[key]
      except:
        value=None
        break
    if value is not None:
      analysis.append(value)
    else:
      #skip the document is one of the value is None
      return None
  return analysis    

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

def perform_operations(variables,ops,globs,debug=False):
    glstr=''
    if globs is not None:
        for var in globs:
            glstr+= "global "+var+"\n"
    #first evaluate the given variables
    for key in variables:
        command=key+"="+str(variables[key])
        if debug: print command
        exec(command)
        #then evaluate the given expression
    if debug: print ops
    exec(glstr+ops, globals(), locals())

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
    

def parse_arguments():
  parser = optparse.OptionParser("This script is used to extract some information from a logfile")
  parser.add_option('-d', '--data', dest='data',default=None,action="store_true", #sys.argv[2],
                    help="BigDFT logfile, yaml compliant (check this if only this option is given)", metavar='FILE')
  parser.add_option('-v', '--invert-match', dest='remove',default=None, #sys.argv[2],
                    help="File containing the keys which have to be excluded from the logfile", metavar='FILE')
  parser.add_option('-a', '--analyze', dest='analyze',default=None, #sys.argv[2],
                    help="File containing the keys which have to be analysed and extracted to build the quantities", metavar='FILE')
  parser.add_option('-e', '--extract', dest='extract',default=None, #sys.argv[2],
                    help="File containing the keys which have to be extracted to build the quantities", metavar='FILE')
  parser.add_option('-o', '--output', dest='output', default="/dev/null", #sys.argv[4],
                    help="set the output file (default: /dev/null)", metavar='FILE')
  parser.add_option('-t', '--timedata', dest='timedata',default=False,action="store_true",
                    help="BigDFT time.yaml file, a quick report is dumped on screen if this option is given", metavar='FILE')
  parser.add_option('-n', '--name', dest='name',default=None,
                    help="Give a name to the set of the plot represented", metavar='FILE')
  parser.add_option('-p', '--plot', dest='plottype',default='Seconds',
                    help="Decide the default yscale for the plotting", metavar='FILE')
  parser.add_option('-s', '--static', dest='static',default=False,action="store_true",
                    help="Show the plot statically for screenshot use", metavar='FILE')
  parser.add_option('-f', '--fontsize', dest='fontsize',default=15,
                    help="Determine fontsize of the bar chart plot", metavar='FILE')
  parser.add_option('-k', '--remove-key', dest='nokey',default=False,action="store_true",
                    help="Remove the visualisation of the key from the main plot", metavar='FILE')

  
  #Return the parsing
  return parser

def nopc(string):
  pc=string.rstrip("%")
  try:
    fl=float(pc)
  except:
    fl=-1.0
  return fl

def dump_timing_level(level,ilev=0,theta=0,
                      data={"time":[],"level":[],"theta":[],"names":[[]]}):
  """Inspect the first level of the given dictionary and dump the profile subroutines at this level"""
  import pylab
  subs=data["time"]
  tht=data["theta"]
  lev=data["level"]
  nms=data["names"]
  if ilev == len(nms):
    nms.append([])
  #tel=0
  tel=theta #entry point of the routine
  for routine in level:
    #first eliminate the subroutines from the level
    try:
      sublevel=routine.pop("Subroutines")
    except:
      sublevel=None
    for name in routine:
      t0=routine[name][0] #take the time
      subs.append(t0) #take the time
      tht.append(tel)
      lev.append(ilev)
      nms[ilev].append(name)
    if sublevel is not None:
      jlev=ilev+1
      dump_timing_level(sublevel,ilev=jlev,theta=tel,data=data)
    tel=tel+t0
  return data

#defines the class of the polar plot
class polar_axis():
  def __init__(self,data):
    import pylab
    self.fig = pylab.figure()
    self.ax = pylab.axes([0.025,0.025,0.95,0.95], polar=True)
    self.tot=data["time"][0]
    self.times=data["time"]
    self.N=len(self.times)
    self.step=5
    self.width=pylab.np.array(data["time"])/self.tot*2*pylab.np.pi
    self.theta=pylab.np.array(data["theta"])/self.tot*2*pylab.np.pi
    self.bot=pylab.np.array(data["level"])
    self.radii = pylab.np.array(self.N*[self.step])
    self.bars = pylab.bar(self.theta,self.radii,
                          width=self.width,
                          bottom=self.step*self.bot,picker=True)
    self.names=data["names"]
    
    ilev=0
    #maxlev=max(self.bot)
    for r,bar,ilev in zip(self.radii, self.bars,self.theta):
       #print ilev,'hello',float(ilev)/float(N),maxlev
       #bar.set_facecolor( pylab.cm.jet(float(ilev)/maxlev))
       bar.set_facecolor(pylab.cm.jet(float(ilev)/(2*pylab.np.pi)))
       bar.set_alpha(0.5)
       ilev+=1

    self.ax.set_xticklabels([])
    self.ax.set_yticklabels([])
    self.info=self.ax.text(0,0,'',fontsize=150)
    # savefig('../figures/polar_ex.png',dpi=48)
    self.fig.canvas.mpl_connect('pick_event',self.info_callback)

  # to be used what it is the moment to show the plot(s)
  def show(self):
    import pylab
    try:
      pylab.show()
    except KeyboardInterrupt:
      raise
        
  def find_name(self,th,level):
    import pylab
    levth=[]
    ipiv=[]
    for i in range(self.N):
      if self.bot[i]==level:
        levth.append(self.theta[i])
        ipiv.append(i)
    to_find_zero=pylab.np.array(levth)-th
    if min(abs(to_find_zero))==0.0:
      routine=pylab.np.argmin(abs(to_find_zero))
      return (th,self.names[level][routine],self.times[ipiv[routine]])
    else:
      to_find_zero[to_find_zero > 0]=-4*pylab.np.pi
      routine=pylab.np.argmin(-to_find_zero)
      return (self.theta[ipiv[routine]],self.names[level][routine],self.times[ipiv[routine]])

    
  def info_string(self,xdata,level):
    #create string to be plotted in plot
    (tt,name,time)=self.find_name(xdata,level)
    info=name+':\n'+str(time)+" s ("
    info+= " %6.2f %s" % (time/self.tot*100.0,"% ) \n")
    #find lower level
    it=range(level)
    it.reverse()
    parents=''
    for i in it:
      (tt,name,time)=self.find_name(tt,i)
      parents+="-"+name+"("+str(time)+"s,"
      parents+= " %6.2f %s" % (time/self.tot*100.0,"% ) \n")
    if len(parents) > 0:
      info+="Calling path:\n"+parents
    return info
        
    
  def info_callback(self,event):
    thisline = event.artist
    xdata, ydata = thisline.get_xy()
    level = ydata/self.step
    #once that the level has been found filter the list of theta
    (tt,name,time)=self.find_name(xdata,level)
    #print "======================="
    #print "Routine Picked:",name,",",time,"s (",time/self.tot*100.0,"% of total time)"
    #find lower level
    #it=range(level)
    #it.reverse()
    #for i in it:
    #  (tt,name,time)=self.find_name(tt,i)
    #  print "  Called by:",name

    #then plot the string
    self.ax.texts.remove(self.info)
    offset = 0.02
    #print 'nowinfo',info
    self.info= self.ax.text( offset, offset, self.info_string(xdata,level),
                             fontsize = 15,transform = self.ax.transAxes )
    self.fig.canvas.draw()

          
class BigDFTiming:
  def __init__(self,filenames,args):
    #here a try-catch section should be added for multiple documents
    #if (len(filename) > 1
    self.log=[]
    for filename in filenames:
      try:
        self.log+=[yaml.load(open(filename, "r").read(), Loader = yaml.CLoader)]
      except:
        self.log+=yaml.load_all(open(filename, "r").read(), Loader = yaml.CLoader)

    #cat=self.log["Dictionary of active categories"]
    #print 'categories',cat
    #icat=1
    #for i in cat:
    #  print icat,i["Category"]
    #  icat +=1
    #dssad
    #shallow copies of important parts
    self.routines=[]
    self.hostnames=[]
    self.scf=[]
    self.ids=[]
    self.barfig = None
    self.axbars = None
    self.newfigs =[]
    self.radio = None
    self.toggle_unbalancing = False
    self.quitButton = None
    self.plot_start=args.plottype
    self.static = args.static
    self.fontsize=args.fontsize
    self.nokey=args.nokey
    for doc in self.log:
        self.routines.append(doc.get("Routines timing and number of calls"))
        self.hostnames.append(doc.get("Hostnames"))
        scf=doc.get("WFN_OPT")
        if scf is not None:
            self.scf.append(scf)
            if "Run name" in doc:
                self.ids.append(doc["Run name"])
            else:
                mpit=doc.get("CPU parallelism")
                if mpit is not None:
                    self.ids.append(mpit["MPI tasks"])
                else:
                    self.ids.append("Unknown")
    self.classes=["Communications","Convolutions","BLAS-LAPACK","Linear Algebra",
            "Other","PS Computation","Potential",
            "Flib LowLevel","Initialization","Unknown"]

  def bars_data(self,vals=None,title='Time bar chart'):
    """Extract the data for plotting the different categories in bar chart"""
    import numpy as np
    import matplotlib.pyplot as plt
    from pylab import cm as cm
    from matplotlib.widgets import Button,RadioButtons
    if vals is None:
      self.vals=self.plot_start
    else:
      self.vals=vals
    if self.barfig is None:
      self.barfig, self.axbars = plt.subplots()
      if self.static: self.barfig.patch.set_facecolor("white")
    dict_list=self.scf
    self.plts=[]
    self.draw_barplot(self.axbars,self.collect_categories(dict_list,self.vals),self.vals,title=title,nokey=self.nokey)
    active=0
    if self.vals == 'Percent':
      self.axbars.set_yticks(np.arange(0,100,10))
      active=1
    if self.radio is None and not self.static:
      self.radio = RadioButtons(plt.axes([0.0, 0.75, 0.08, 0.11], axisbg='lightgoldenrodyellow'), ('Seconds', 'Percent'),active=1)
      self.radio.on_clicked(self.replot)

    if self.quitButton is None and not self.static:
      self.quitButton = Button(plt.axes([0.0, 0.0, 0.1, 0.075]), 'Quit')
      self.quitButton.on_clicked(self.onclick_quitButton)
      self.barfig.canvas.mpl_connect('pick_event',self.onclick_ev)
      self.barfig.canvas.mpl_connect('key_press_event',self.onkey_ev)

  def find_items(self,category,dict_list):
    """For a given category find the items which has them"""
    import numpy as np
    items={}
    for idoc in range( len(dict_list) ):
        for cat in dict_list[idoc]["Categories"]:
            dicat=dict_list[idoc]["Categories"][cat]
            if dicat["Class"] == category:
                if cat not in items:
                    items[cat]=np.zeros(len(dict_list))
                items[cat][idoc]=dicat["Data"][self.iprc]
    return [ (cat,items[cat]) for cat in items]
    
  def collect_categories(self,dict_list,vals):
    """ Collect all the categories which belong to all the dictionaries """
    import numpy as np
    #classes=[]
    #for doc in dict_list:
    #  for cat in doc["Classes"]:
    #    if cat not in classes and cat != 'Total': classes.append(cat)
    #print 'found',classes
    if vals == 'Percent':
      self.iprc=0
    elif vals == 'Seconds':
      self.iprc=1
    catsdats=[]
    self.values_legend=[]
    for cat in self.classes:
      try:
        if cat == "Unknown":
            print "here"
            data_unk=[]
            for doc in dict_list:
                percent_unk=100.0-doc["Classes"]["Total"][0]
                if self.iprc==0:
                    data_unk.append(percent_unk)
                elif self.iprc==1:
                    time_unk=(doc["Classes"]["Total"][1])*percent_unk/100.0
                    data_unk.append(time_unk) 
            dat=np.array(data_unk)
        else:
            dat=np.array([doc["Classes"][cat][self.iprc] for doc in dict_list])
        catsdats.append((cat,dat))
        self.values_legend.append(cat)
      except Exception,e:
        print 'EXCEPTION FOUND',e
        print "category",cat,"not present everywhere"
    return catsdats

  def onkey_ev(self,event):
    number=event.key
    if number == 'u' or number == 'U':
      self.toggle_unbalancing = not self.toggle_unbalancing
      print 'Unbalancing',self.toggle_unbalancing
    try:
      number=int(number)
      if not self.toggle_unbalancing:
        if self.routines[number] is not None:
          toplt=self.routines[number]
          data=dump_timing_level(toplt)
          plt=polar_axis(data)
          plt.show()
      else:
        if self.scf[number] is not None:
          self.load_unbalancing(self.scf[number]["Classes"],self.hostnames[number])
          plt.show()
    except:
      print 'not present or out of range'

  def onclick_ev(self,event):
    import matplotlib.pyplot as plt
    thisline = event.artist
    xdata, ydata = thisline.get_xy()
    #print 'data',xdata,ydata
    #find the category which has been identified
    y0data=0.0
    for cat in self.values_legend:
      y0data+=self.scf[xdata]["Classes"][cat][self.iprc]
      #print 'cat,y0data',cat,y0data,ydata
      if y0data > ydata:
        category=cat
        break
    print 'category',category
    print self.find_items(category,self.scf)
    #self.axbars.cla()
    #self.draw_barplot(self.axbars,self.find_items(category,self.scf),self.vals)
    #self.barfig.canvas.draw()

    newfig=plt.figure()
    newax=newfig.add_subplot(111)
    self.draw_barplot(newax,self.find_items(category,self.scf),self.vals,title=category)
    newfig.show()
    self.newfigs.append((newfig,newax))
  
    
  def draw_barplot(self,axbars,data,vals,title='Time bar chart',static=False,nokey=False):
    import numpy as np
    import matplotlib.pyplot as plt
    from pylab import cm as cm
    from matplotlib.widgets import Button,RadioButtons

    ndata=len(data[0][1])
    ind=np.arange(ndata)
    width=0.9#0.35
    bot=np.array(ndata*[0.0])
    icol=1.0
    for cat,dat in data:
      print 'cat',cat,dat
      for i in range(len(self.ids)):
          print self.ids[i],dat[i]
      plt=axbars.bar(ind,dat,width,bottom=bot,color=cm.jet(icol/len(self.classes)),picker=True,label=cat)
      self.plts.append(plt)
      bot+=dat
      icol+=1.0
    drawn_classes=np.array(self.values_legend)
    axbars.set_title(title,fontsize=self.fontsize*1.2)
    axbars.set_ylabel(vals,fontsize=self.fontsize)
    axbars.set_xticks(ind+width/2.)
    axbars.set_xticklabels(np.array(self.ids),size=self.fontsize)
    if not nokey:
      self.leg = axbars.legend(loc='upper right',fontsize=self.fontsize)
      self.leg.get_frame().set_alpha(0.4)  
          
  def onclick_quitButton(self,event):
    print "Good bye!"
    for figax in self.newfigs:
      pylab.close(figax[0])
    pylab.close(self.barfig)
    
  def replot(self,label):
    title=self.axbars.get_title()
    self.axbars.cla()
    self.bars_data(vals=label,title=title)
    self.barfig.canvas.draw()
    for figax in self.newfigs:
      ax=figax[1]
      fi=figax[0]
      category=ax.get_title()
      ax.cla()
      self.draw_barplot(ax,self.find_items(category,self.scf),self.vals,title=category)
      fi.canvas.draw()


  def func(self,label):
    print 'label,cid',label,self.cid
    self.barfig.draw()

  def unbalanced(self,val):
    """Criterion for unbalancing"""
    return val > 1.1 or val < 0.9

  def find_unbalanced(self,data):
    """Determine lookup array of unbalanced categories"""
    ipiv=[]
    for i in range(len(data)):
      if self.unbalanced(data[i]):
        ipiv.append(i)
    return ipiv

  def load_unbalancing(self,dict,hosts):
    """Extract the data for plotting the hostname balancings between different categories in bar chart"""
    import matplotlib.pyplot as plt
    import numpy as np
    width=0.50
    plts=[]
    key_legend=[]
    values_legend=[]
    icol=1.0
    print "dict",dict
    #open a new figure
    newfig=plt.figure()
    newax=newfig.add_subplot(111)
    for cat in self.classes:
      try:
        dat=np.array([dict[cat][0]])
        print 'data',dat
        unb=np.array(dict[cat][2:])
        print 'unbalancing',unb
        unb2=self.find_unbalanced(unb)
        print 'unbalanced objects',cat
        if hosts is not None and (cat=='Convolutions' or cat =='Communications'):
          print 'unb2',unb2,len(unb),len(hosts)
          print 'vals',[ [i,unb[i],hosts[i]] for i in unb2]
        ind=np.arange(len(unb))
        plt=pylab.bar(ind,unb,width,color=pylab.cm.jet(icol/len(self.classes)))
        plts.append(plt)
        key_legend.append(plt[0])
        values_legend.append(cat)
        icol+=1.0
        if (width > 0.05):
          width -= 0.05
      except Exception,e:
        print 'EXCEPTION FOUND',e
        print "cat",cat,"not found"

    if len(ind) > 2:
      if hosts is not None:
        tmp=np.array(hosts)
      else:
        tmp=None
    else:
      tmp=np.array(["max","min"])
    pylab.ylabel('Load Unbalancing wrt average')
    pylab.title('Work Load of different classes')
    if tmp is not None: pylab.xticks(ind+width/2., tmp,rotation=90,verticalalignment='bottom')
    pylab.yticks(pylab.np.arange(0,2,0.25))
    pylab.legend(pylab.np.array(key_legend), pylab.np.array(values_legend))
    newfig.show()
    self.newfigs.append((newfig,newax))


if __name__ == "__main__":
  parser = parse_arguments()
  (args, argcl) = parser.parse_args()


#logfile
#check if timedata is given
if args.timedata:
  import pylab
  print 'args of time',args.timedata,argcl
  if args.name is not None:
    title=args.name
  else:
    title='Time bar chart'
  #in the case of more than one file to analyse
  #or in the case of more than one yaml document per file
  #just load the bars data script
  
  #load the first yaml document
  bt=BigDFTiming(argcl,args)
  print "hosts",bt.hostnames
  if bt.scf is not None:
    bt.bars_data(title=title) #timing["WFN_OPT"]["Classes"])
    
  if bt.scf[0] is not None and False:
    bt.load_unbalancing(bt.scf[0]["Classes"]) #timing["WFN_OPT"]["Classes"])
      
  if bt.routines[0] is not None and False:
    data=dump_timing_level(bt.routines[0]) #dict_routines)
    plt=polar_axis(data)
    plt.show()
  else:
    pylab.show()

  #print allev
  #dump the loaded info
  

if args.data is None:
  print "No input file given, exiting..."
  exit(0)

if args.analyze is not None and args.data:
  instructions= yaml.load(open(args.analyze, "r").read(), Loader = yaml.CLoader)
  if INITIALIZATION in instructions:
      for var in instructions[INITIALIZATION]:
          exec var +" = "+ str(instructions[INITIALIZATION][var])
  print '#',args.data,argcl
  for f in argcl:
    sys.stderr.write("#########processing "+f+"\n")
    datas=get_logs([f])
    for doc in datas:
      doc_res=document_quantities(doc,instructions)
      #print doc_res,instructions
      if EVAL in instructions: perform_operations(doc_res,instructions[EVAL],instructions.get(INITIALIZATION))
  exit(0)

if args.data:
  with open(argcl[0], "r") as fp:
    logfile_lines = fp.readlines()

    
#output file
file_out=open(args.output, "w")
#to_remove list
if args.remove is not None:
  to_remove = yaml.load(open(args.remove, "r").read(), Loader = yaml.CLoader)
else:
  #standard list which removes long items from the logfile
  to_remove= ["Atomic positions within the cell (Atomic and Grid Units)",
              "Atomic Forces (Ha/Bohr)",
              #"Orbitals",
              #"Energies",
              "Properties of atoms in the system"]
  #to_remove=[]

    
#obtain the cleaned document
cleaned_logfile = clean_logfile(logfile_lines,to_remove)

print 'Logfile cleaned, writing output file...'
#dump it
for line in cleaned_logfile:
  file_out.write(line) 

file_out.close()
print 'Output file written'
#sys.exit(0)
#experimental: extract for any document of the cleaned logfile the final energy and positions
if args.extract is not None:
  to_extract = yaml.load(open(args.extract, "r").read(), Loader = yaml.CLoader)
else:
  to_extract=[ ["Geometry","FORCES norm(Ha/Bohr)","maxval"], ["Last Iteration","EKS"]]

datas=yaml.load_all(''.join(cleaned_logfile), Loader = yaml.CLoader)
extracted_result=[]
for doc in datas:
  doc_res=document_analysis(doc,to_extract)
  print doc_res,to_extract
  if doc_res is not None:
    extracted_result.append(doc_res)

print "Number of valid documents:",len(extracted_result)
for it in extracted_result:
  print it

exit(0)
    
iterations = range(len(extracted_result))
energies = [en for [f, en] in extracted_result]
energy_min=min(energies)
energies = [en-energy_min for en in energies]
forces = [f for [f, en] in extracted_result]

import matplotlib.pyplot as plt
#plt.plot(iterations, energies, '.-',label='E - min(E)')
plt.plot(energies, forces, '.-',label='Energy')
#plt.yscale('log')
plt.legend(loc='lower left')
plt.show()
  
sys.exit(0)

#do some tests
document = """---
block sequence:
 - BlockEntryToken
 block mapping:
   ? KeyToken
   : ValueToken
 flow sequence: [FlowEntryToken, FlowEntryToken]
 flow mapping: {KeyToken: ValueToken}
 anchors and tags:
 - &A !!int '5'
 - *A
...
 """

#for token in yaml.scan(document):
#     print token

#ddd
#print args.ref,args.data,args.output

datas    = [a for a in yaml.load_all(open(argcl[0], "r").read(), Loader = yaml.CLoader)]
#i=0
#for doc in yaml.load_all(open(args.data, "r").read(), Loader = yaml.CLoader):
#  i+=1
#  print 'doc No.',i#,"Last Iteration" in doc.keys()
  #try:
  #  print yaml.dump(doc["Last Iteration"])
  #except:
  #  print "Last Iteration not found"

#datas    = [a for a in yaml.safe_load_all(open(args.data, "r").read())]
#Profile.run('datas    = [a for a in yaml.load_all(open(args.data, "r").read(), Loader = yaml.CLoader)]')
#gyi
#runfile = open(args.data, "r").read()
sss
try:
    datas    = [a for a in yaml.load_all(runfile, Loader = yaml.Loader)]
    #Profile.run('datas    = [a for a in yaml.load_all(open(args.data, "r").read())]')
#for a in yaml.load_all(open(args.data, "r").read(), Loader = yaml.CLoader):
#    print 'New Document'
#    document_analysis(a)
  
except:
  datas = []
  reports = open(args.output, "w")
  fatal_error(args,reports)

ndocs=len(datas)

print 'No. of Documents:',ndocs

for run in datas:
  print 'New Document'
  document_analysis(run)



print 'end'

test = run['Atomic positions within the cell (Atomic and Grid Units)']

#value = test[1]["Ag"]["AU"][1]

#ex={}
#ex["value"]=complex(12,13)
#print 'example',value,3+4j,value,(3+4j),complex(7,8),ex["value"]+complex(value),test[1]["Ag"]["AU"][0]+test[1]["Ag"]["AU"][2]

sys.stdout.write(yaml.dump(run["Last Iteration"],default_flow_style=False,explicit_start=True))
#sys.stdout.write(yaml.dump(run["Direct and transposed data repartition"],default_flow_style=False,explicit_start=True))
#sys.stdout.write(yaml.dump(ex,default_flow_style=False,explicit_start=True))

sys.exit(0)


