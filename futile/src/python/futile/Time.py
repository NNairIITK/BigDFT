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

          
class BigDFTiming:
  def __init__(self,filenames,args):
    #here a try-catch section should be added for multiple documents
    #if (len(filename) > 1
    import yaml
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
    import pylab
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
    import pylab
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
