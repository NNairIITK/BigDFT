def axis_from_data(fig,ax,data):
  "Transform a data tuple into axis coordinates"
  trans=fig.transFigure.inverted()
  ll=trans.transform(ax.transData.transform(data))
  return list(ll)

class FigureSet():
  """Define a container for a plot with the possiblity to switch between simple and gnuplot plotting"""
  def __init__(self,**kwargs):
    import matplotlib.pyplot as plt
    from Yaml import kw_pop
    newkw,title=kw_pop('title','',**kwargs)    
    self.title=title
    self.figures=[]
    self.showing=False
    self.add(**newkw)

  def __getitem__(self,idx):
    if isinstance(idx,int):
      return self._get_figure(idx)
    else:
      return self._get_figure(figname=idx)

  def _locate(self,figname):
    for i,fig in enumerate(self.figures):
      if figname==fig['Title']: 
        return i,fig
    return None

  def exists(self,figname):
    "True if the Figure exists in the Set"
    return self._locate(figname) is not None

  def invoke(self,idx):
    import matplotlib.pyplot as plt
    "Invoke the Figure if it exists. Return None otherwise"
    if isinstance(idx,int):
      if idx >= len(self.figures): return None
      i=idx
    else:
      if not self.exists(idx): return None
      i,fig=self._locate(idx)
    figname=self.figures[i]['TrueTitle']
    return plt.figure(figname)

  def add(self,**kwargs):
    import matplotlib.pyplot as plt
    toadd=str(len(self.figures)+1)
    title=kwargs.get('title','Figure '+toadd if self.title is '' else self.title)
    if self.title != '' and len(self.figures)>0: 
      newtitle=title+' ('+self.title+')'
    else:
      newtitle=title
    quitbutton=kwargs.get('QuitButton',False)
    newfig=plt.figure(newtitle)
    axargs=kwargs.get('axes')
    #cfm=plt.get_current_fig_manager()
    if axargs is None:
      newax=newfig.add_subplot(111)
    else:
      rect=axargs.pop('rect')
      newax=plt.axes(rect,**axargs)
    #connect the home key execpt for the first figure
    from functools import partial
    newfig.canvas.mpl_connect('key_press_event',partial(self._onkey_home,len(self.figures)))
    figdict={'Figure':newfig,'Axes':newax,'Title':title,'TrueTitle':newtitle}
    if quitbutton: figdict.update(self._put_quitbutton())
    self.figures.append(figdict)
    return newfig,newax

  def _onkey_home(self,ifig,event):
    import matplotlib.pyplot as plt
    number=event.key
    #print "pressed",number,'end'
    #return to the main figure if the key "home" is pressed
    if str(number)=='home': 
      #print 'raising window'
      self._raisewindow(0)
    elif number == 'q' or number == 'Q':
      if ifig==0: 
        self._quitall()
      else:
        figax=self.figures[ifig]
        plt.close(figax['Figure'])

  def _put_quitbutton(self):
    import matplotlib.pyplot as plt
    from matplotlib.widgets import Button
    quitbutton=Button(plt.axes([0.0, 0.0, 0.1, 0.075]), 'Quit')
    quitbutton.on_clicked(self._onclick_quitButton)
    return {'QuitButton':quitbutton}

  def _raisewindow(self,index):
    import matplotlib.pyplot as plt
    if self.invoke(index) is None: return
    cfm=plt.get_current_fig_manager()
    try:
      cfm.window.activateWindow()
      cfm.window.raise_()
    except:
      cfm.window.attributes('-topmost', True)
      #cfm.window.attributes('-topmost', False)
    return cfm

  def _quitall(self):
    import matplotlib.pyplot as plt
    for figax in self.figures:
      plt.close(figax['Figure'])
    self.showing=False
      
  def _onclick_quitButton(self,event):
    print "Good bye!"
    self._quitall()

  def _get_figure(self,index=0,figname=None):
    if figname is not None:
      for fig in self.figures:
        if figname==fig['Title']: break
    else:
      fig=self.figures[index]
    fx=fig["Figure"]
    ax=fig["Axes"]
    return fx,ax

  def show(self,figname=None):
    import matplotlib.pyplot as plt
    if figname is not None and self.exists(figname): 
      fx,ax=self._get_figure(figname=figname)
      if not self.showing: self.show()
      fx.show()
    elif figname is None:
      self.showing=True
      plt.show()

#defines the class of the polar plot
class polar_axis():
  step=5
  def __init__(self,fig,ax,data):
    self.fig = fig #pylab.figure()
    self.ax = ax #pylab.axes([0.025,0.025,0.95,0.95], polar=True)
    self.data=data
    self.reset=True
    self.ax.cla()
    datatmp=self.data
    self.subdata=self.dump_timing_level(datatmp)
    self.draw_polarplot()
    self.fig.canvas.mpl_connect('pick_event',self._new_level)
    self.fig.canvas.mpl_connect('motion_notify_event',self._info_callback)

  def draw_polarplot(self):
    import pylab
    import numpy as np
    #self.subdata=data
    self.tot=self.subdata["time"][0]
    self.times=self.subdata["time"]
    self.N=len(self.times)

    self.width=np.array(self.subdata["time"])/self.tot*2*np.pi
    self.theta=np.array(self.subdata["theta"])/self.tot*2*np.pi
    self.bot=np.array(self.subdata["level"])
    self.names=self.subdata["names"]
    self.radii = np.array(self.N*[self.step])
    
    self.bars = self.ax.bar(self.theta,self.radii,
                          width=self.width,
                          bottom=self.step*self.bot,picker=True)
   
    ilev=0
    #maxlev=max(self.bot)
    for r,bar,ilev in zip(self.radii, self.bars,self.theta):
       #print ilev,'hello',float(ilev)/float(N),maxlev
       #bar.set_facecolor( pylab.cm.jet(float(ilev)/maxlev))
       bar.set_facecolor(pylab.cm.jet(float(ilev)/(2.0*np.pi)))
       bar.set_alpha(0.5)
       ilev+=1

    self.ax.set_xticklabels([])
    self.ax.set_yticklabels([])
    self._wrotesome=False
    self.info=self.ax.text(0,0,'',fontsize=150)

  ## to be used what it is the moment to show the plot(s)
  #def show(self):
  #  import pylab
  #  try:
  #    pylab.show()
  #  except KeyboardInterrupt:
  #    raise
        
  def _find_name(self,th,level):
    import numpy as np
    levth=[]
    ipiv=[]
    for i in range(self.N):
      if self.bot[i]==level:
        levth.append(self.theta[i])
        ipiv.append(i)
    to_find_zero=np.array(levth)-th
    to_find_zero[to_find_zero > 0]=-4*np.pi
    routine=np.argmin(-to_find_zero)
    if routine == len(to_find_zero)-1 and abs(routine)> 0.5:
      return None
    else:
      return (self.theta[ipiv[routine]],self.names[level][routine],self.times[ipiv[routine]])
    #if min(abs(to_find_zero))==0.0:
    #  routine=np.argmin(abs(to_find_zero))
    #  return (th,self.names[level][routine],self.times[ipiv[routine]])
    #else:
    #  to_find_zero[to_find_zero > 0]=-4*np.pi
    #  routine=np.argmin(-to_find_zero)
    #  return (self.theta[ipiv[routine]],self.names[level][routine],self.times[ipiv[routine]])

    
  def _info_string(self,xdata,level):
    #create string to be plotted in plot
    (tt,name,time)=self._find_name(xdata,level)
    info=name+':\n'+str(time)+" s ("
    info+= " %6.2f %s" % (time/self.tot*100.0,"% ) \n")
    #find lower level
    it=range(level)
    it.reverse()
    parents=''
    for i in it:
      (tt,name,time)=self._find_name(tt,i)
      parents+="-"+name+"("+str(time)+"s,"
      parents+= " %6.2f %s" % (time/self.tot*100.0,"% ) \n")
    if len(parents) > 0:
      info+="Calling path:\n"+parents
    return info

  def _new_level(self,event):
    click=event.mouseevent.button
    datatmp=self.data
    if click==3:
      if not self.reset:
        self.ax.texts.remove(self.comebackinfo)
        self.reset=True
      nameandlev=None
    elif click==1:
      thisline = event.artist
      xdata, ydata = thisline.get_xy()
      level = int(ydata)/self.step
      goodpoint=self._find_name(xdata,level)
      if goodpoint is None: return
      (tt,name,time)=goodpoint
      if not self.reset: return #no double level
      self.reset=False
      nameandlev=(name,level)
      #extract the data to be plotted for the internal level
    self.subdata=self.dump_timing_level(datatmp,starting_point=nameandlev)
    self.ax.cla()
    self.draw_polarplot()
    if not self.reset: 
      self.comebackinfo=self.ax.text(0.02,0.95,
                                     'Subroutines of '+name+' (Right click to reset)',
                                     fontsize=15,transform = self.ax.transAxes)
    self.fig.canvas.draw()
    
  def _info_callback(self,event):
    #print dir(event)
    #thisline = event.artist
    #xdata, ydata = thisline.get_xy()
    xdata=event.xdata#data
    ydata=event.ydata#data
    #print 'test',xdata,ydata,event.x,event.y
    if xdata is None or ydata is None: return
    level = int(ydata)/self.step
    goodpoint=self._find_name(xdata,level)
    #print 'level',level
    #once that the level has been found filter the list of theta
    #(tt,name,time)=self._find_name(xdata,level)
    #print "======================="
    #print "Routine Picked:",name,",",time,"s (",time/self.tot*100.0,"% of total time)"
    #find lower level
    #it=range(level)
    #it.reverse()
    #for i in it:
    #  (tt,name,time)=self.find_name(tt,i)
    #  print "  Called by:",name

    #then plot the string
    if self._wrotesome:
      self.ax.texts.remove(self.info)
      self.info=self.ax.text(0,0,'',fontsize=150)
      self._wrotesome=False
    if goodpoint is None: return
    offset = 0.02
    #print 'nowinfo',info
    self.info= self.ax.text( offset, offset, self._info_string(xdata,level),
                             fontsize = 15,transform = self.ax.transAxes )
    self._wrotesome=True
    self.fig.canvas.draw()

  def dump_timing_level(self,level,starting_point=None,ilev=0,theta=0,data=None): ##data={"time":[],"level":[],"theta":[],"names":[[]]},reset=False):
    """Inspect the first level of the given dictionary and dump the profile subroutines at this level"""
    import pylab
    if data is None: data={"time":[],"level":[],"theta":[],"names":[[]]}
    subs=data["time"]
    tht=data["theta"]
    lev=data["level"]
    nms=data["names"]
    if ilev == len(nms):
      nms.append([])
    #tel=0
    tel=theta #entry point of the routine
    t0=0
    start=starting_point
    import copy
    leveltmp=copy.deepcopy(level)
    brk=False
    for routine in leveltmp:
      #first eliminate the subroutines from the level
      try:
        sublevel=routine.pop("Subroutines")
      except:
        sublevel=None
      for name in routine:
        if starting_point is not None:
          if (name,ilev) != starting_point:
            continue
          else:
            brk=True
            start=None #from now on take the data
            ilev=0
        t0=routine[name][0] #take the time
        subs.append(t0) #take the time
        tht.append(tel)
        lev.append(ilev)
        nms[ilev].append(name)
      if sublevel is not None:
        jlev=ilev+1
        self.dump_timing_level(sublevel,starting_point=start,ilev=jlev,theta=tel,data=data)
      tel=tel+t0
      if brk: break
    return data
  
class TimeData:
  barwidth=0.9#0.35
  def __init__(self,*filenames,**kwargs):
    """
    Class to analyze the timing from a Futile run
    ------
    Arguments
    filenames: list of the files which have to be treated
    plottype: Decide the default yscale for the plotting
              May be "Percent" or "Seconds" 
    static: Show the plot statically for screenshot use
    fontsize: Determine fontsize of the bar chart plot 
    nokey: Remove the visualization of the key from the main plot
    """
    #here a try-catch section should be added for multiple documents
    #if (len(filename) > 1
    import yaml
    self.log=[]
    for filename in filenames:
      try:
        self.log+=[yaml.load(open(filename, "r").read(), Loader = yaml.CLoader)]
      except:
        self.log+=yaml.load_all(open(filename, "r").read(), Loader = yaml.CLoader)
    #create the figure environemnt
    self.figures=FigureSet(title='Profiling',QuitButton=True)
    #self.barfig = None
    #self.axbars = None
    #self.newfigs =[]
    self.radio = None
    self.toggle_unbalancing = False
    self.quitButton = None
    self.lined=[]
    self.plot_start=kwargs.get('plottype','Seconds')
    self.static = kwargs.get('static',False)
    self.fontsize=kwargs.get('fontsize',15)
    self.nokey=kwargs.get('nokey',False)
    counter=kwargs.get('counter','WFN_OPT') #the default value, to be customized
    self.inspect_counter(counter)

  def _refill_classes(self,main):
    self.classes=[]
    for doc in self.log:
        self.routines.append(doc.get("Routines timing and number of calls"))
        self.hostnames.append(doc.get("Hostnames"))
        scf=doc.get(main)
        if scf is not None:
            self.scf.append(scf)
            loclass=scf["Classes"].keys()
            for cs in loclass:
              if cs not in self.classes and cs != "Total": self.classes.append(cs)
            if "Run name" in doc:
                self.ids.append(doc["Run name"])
            else:
                mpit=doc.get("CPU parallelism")
                if mpit is not None:
                    self.ids.append(mpit["MPI tasks"])
                else:
                    self.ids.append("Unknown")
    #to be updated
    self.classes.append("Unknown") #an evergreen
    self.classes.sort()
    #self.classes=["Communications","Convolutions","BLAS-LAPACK","Linear Algebra",
    #        "Other","PS Computation","Potential",
    #        "Flib LowLevel","Initialization","Unknown"]
      
  def inspect_counter(self,counter,unit=None):
    self.routines=[]
    self.hostnames=[]
    self.scf=[]
    self.ids=[]
    self.vals=self.plot_start if unit is None else unit
    self._refill_classes(counter)
    #here we might add the policy to add new figures or delete previous ones
    data=self.collect_categories(self.scf,self.vals)
    fig,axis=self.figures[0]
    axis.cla()
    self.draw_barfigure(fig,axis,data,title="Counter "+counter)

  def inspect_category(self,cat):
    data=self.find_items(cat,self.scf)
    if self.figures.exists(cat):
      newfig,newax=self.figures[cat]
      #erase the plot in case it is already available
      newax.cla()
    else:
      newfig,newax=self.figures.add(title=cat)
    self.draw_barfigure(newfig,newax,data,title=cat)
    newfig.show()

  def show(self):
    self.figures.show()

  def draw_barfigure(self,fig,axis,data,title):
    import numpy as np
    from matplotlib.widgets import RadioButtons,Button
    import matplotlib.pyplot as plt
    if self.static: fig.patch.set_facecolor("white")
    self.draw_barplot(axis,data,self.vals,title=title,nokey=self.nokey)
    if self.vals == 'Percent':
      axis.set_yticks(np.arange(0,100,10))
    if self.radio is None and not self.static:
      self.radio = RadioButtons(plt.axes([0.0, 0.75, 0.08, 0.11], axisbg='lightgoldenrodyellow'), ('Seconds', 'Percent'),
                                active=1 if self.vals=='Percent' else 0)
      self.radio.on_clicked(self.replot)
      tt=axis.get_xticks()
      routics=[axis_from_data(fig,axis,(tic-0.45*self.barwidth,0.))[0] for tic in tt]
      unbtics=[axis_from_data(fig,axis,(tic+0.05*self.barwidth,0.))[0] for tic in tt]
      self.routine_buttons=[]
      from functools import partial
      for i,t in enumerate(routics):
        if self.routines[i] is None: continue
        but=Button(plt.axes([t, 0.0, 0.1, 0.05]), 'R')
        but.on_clicked(partial(self.routines_plot,i))
        self.routine_buttons.append(but)
      for i,t in enumerate(unbtics):
        if self.scf[i] is None: continue
        but=Button(plt.axes([t, 0.0, 0.1, 0.05]), 'U')
        but.on_clicked(partial(self.workload_plot,i))
        self.routine_buttons.append(but)
      fig.canvas.mpl_connect('pick_event',self.onclick_ev)
      #fig.canvas.mpl_connect('key_press_event',self.onkey_ev)
    fig.canvas.draw()

  def routines_plot(self,index,event=None):
    "Draw the plot of the routines level for the run identified by index"
    toplt=self.routines[index]
    title='Routines '+str(self.ids[index])+'-'+str(index)
    if self.figures.exists(title): return
    #data=dump_timing_level(toplt)#,starting_point='cluster')
    fig,ax=self.figures.add(title=title,
                            axes={'rect':[0.025,0.025,0.95,0.95],'polar':True})
    if not hasattr(self,'routine_plots'): self.routine_plots=[]
    plt=polar_axis(fig,ax,toplt)
    self.routine_plots.append(plt)
    fig.show()

  def workload_plot(self,index,event=None):
    "Draw the plot of the workload of different classes for the run identified by index"
    from functools import partial
    title='Unbalancing '+str(self.ids[index])+'-'+str(index)
    if self.figures.exists(title): return
    fig,ax=self.figures.add(title=title)
    leg,plts=self.load_unbalancing(ax,self.scf[index]["Classes"],self.hostnames[index])
    idx=len(self.lined)
    self.lined.append({})
    for legline, origline in zip(leg.get_patches(), plts):
          legline.set_picker(5)  # 5 pts tolerance
          self.lined[idx][legline] = origline
    fig.canvas.mpl_connect('pick_event', partial(self._onpick_curve,idx,fig))
    fig.show()

##  def bars_data(self,counter="WFN_OPT",vals=None,title='Time bar chart'):
##    """Extract the data for plotting the different categories in bar chart"""
##    import numpy as np
##    import matplotlib.pyplot as plt
##    from pylab import cm as cm
##    from matplotlib.widgets import Button,RadioButtons
##    if not hasattr(self,'counter'): 
##      self.counter=counter #only for the first time
##      self._refill_classes(self.counter)
##    self.vals=self.plot_start if vals is None else vals
##    if self.barfig is None:
##      self.barfig, self.axbars = plt.subplots()
##      if self.static: self.barfig.patch.set_facecolor("white")
##    dict_list=self.scf
##    self.plts=[]
##    self.draw_barplot(self.axbars,self.collect_categories(dict_list,self.vals),self.vals,title=title,nokey=self.nokey,units=self.vals)
##    active=0
##    if self.vals == 'Percent':
##      self.axbars.set_yticks(np.arange(0,100,10))
##      active=1
##    if self.radio is None and not self.static:
##      self.radio = RadioButtons(plt.axes([0.0, 0.75, 0.08, 0.11], axisbg='lightgoldenrodyellow'), ('Seconds', 'Percent'),
##                                active=1 if self.vals=='Percent' else 0)
##      self.radio.on_clicked(self.replot)
##
##    if self.quitButton is None and not self.static:
##      self.quitButton = Button(plt.axes([0.0, 0.0, 0.1, 0.075]), 'Quit')
##      self.quitButton.on_clicked(self.onclick_quitButton)
##      self.barfig.canvas.mpl_connect('pick_event',self.onclick_ev)
##      self.barfig.canvas.mpl_connect('key_press_event',self.onkey_ev)

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
          print self.routines[number]
          data=dump_timing_level(toplt)#,starting_point='cluster')
          print data
          plt=polar_axis(data)
          plt.show()
      else:
        if self.scf[number] is not None:
          self.load_unbalancing(self.scf[number]["Classes"],self.hostnames[number])
          plt.show()
    except Exception,e:
      print e
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
    if category != "Unknown": self.inspect_category(category)
    #data=self.find_items(category,self.scf)
    #print data
    ##self.axbars.cla()
    ##self.draw_barplot(self.axbars,self.find_items(category,self.scf),self.vals)
    ##self.barfig.canvas.draw()
    #newfig,newax=self.figures.add(title=category)
    #newfig=plt.figure()
    #newax=newfig.add_subplot(111)
    #here we should add the check not to replot the same thing if already done
    #self.draw_barfigure(newfig,newax,data,title=category)
    #self.draw_barplot(newax,self.find_items(category,self.scf),self.vals,title=category)
    #newfig.show()
    #self.newfigs.append((newfig,newax))
  
    
  def draw_barplot(self,axbars,data,vals,title='Time bar chart',static=False,nokey=False):
    import numpy as np
    from pylab import cm as cm
    ndata=len(data[0][1])
    ind=np.arange(ndata)
    bot=np.array(ndata*[0.0])
    width=self.barwidth
    icol=1.0
    for cat,dat in data:
      print 'cat',cat,dat
      for i in range(len(self.ids)):
        print self.ids[i],dat[i]
      plt=axbars.bar(ind,dat,width,bottom=bot,color=cm.jet(icol/len(self.classes)),picker=True,label=cat)
      #self.plts.append(plt)
      bot+=dat
      icol+=1.0
    #drawn_classes=np.array(self.values_legend)
    axbars.set_title(title,fontsize=self.fontsize)
    axbars.set_ylabel(vals,fontsize=self.fontsize)
    axbars.set_xticks(ind+width/2.)
    axbars.set_xticklabels(np.array(self.ids),size=self.fontsize)
    if not nokey:
      self.leg = axbars.legend(loc='upper right',fontsize=self.fontsize)
      self.leg.get_frame().set_alpha(0.4)  
          
##  def onclick_quitButton(self,event):
##    import pylab
##    print "Good bye!"
##    for figax in self.newfigs:
##      pylab.close(figax[0])
##    pylab.close(self.barfig)
##    
  def replot(self,label):
    self.vals=label
    #print self.vals
    for i,fig in enumerate(self.figures):
      fx,ax=fig
      cat=ax.get_title()
      ax.cla()
      if i==0:
        data=self.collect_categories(self.scf,label)
      else:
        data=self.find_items(cat,self.scf)
      if data != []: self.draw_barfigure(fx,ax,data,cat)
      #fx.canvas.draw()
    #fig,axis=self.figures.get_figure(0)
    #title=axis.get_title()
    ##title=self.axbars.get_title()
    #self.axbars.cla()
    #self.bars_data(vals=label,title=title)
    #self.barfig.canvas.draw()
    #for figax in self.newfigs:
    #  ax=figax[1]
    #  fi=figax[0]
    #  category=ax.get_title()
    #  ax.cla()
    #  self.draw_barplot(ax,self.find_items(category,self.scf),self.vals,title=category)
    #  fi.canvas.draw()


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

  def load_unbalancing(self,ax,dict,hosts):
    """Extract the data for plotting the hostname balancings between different categories in bar chart"""
    import numpy as np
    from pylab import cm
    width=0.50
    plts=[]
    key_legend=[]
    values_legend=[]
    icol=1.0
    ##print "dict",dict
    #open a new figure
    #newfig=plt.figure()
    #newax=newfig.add_subplot(111)
    for cat in self.classes:
      if cat=='Unknown': continue
      try:
        dat=np.array([dict[cat][0]])
        #print 'data',dat
        unb=np.array(dict[cat][2:])
        #print 'unbalancing',unb
        #unb2=self.find_unbalanced(unb)
        #print 'unbalanced objects',cat
        #if hosts is not None and (cat=='Convolutions' or cat =='Communications'):
          #print 'unb2',unb2,len(unb),len(hosts)
          #print 'vals',[ [i,unb[i],hosts[i]] for i in unb2]
        ind=np.arange(len(unb))
        pltmp=ax.bar(ind,unb,width,color=cm.jet(icol/len(self.classes)))
        plts.append(pltmp)
        key_legend.append(pltmp[0])
        values_legend.append(cat)
        icol+=1.0
        if (width > 0.05):
          width -= 0.05
      except Exception,e:
        print 'EXCEPTION FOUND',e
        print "cat",cat,"not found in workload data"

    if len(ind) > 2:
      tmp=np.array(hosts) if hosts is not None else None
    else:
      tmp=np.array(["max","min"])
    ax.set_ylabel('Load Unbalancing wrt average')
    ax.set_title('Work Load of different classes')
    if tmp is not None: 
      ax.set_xticks(ind+width/2.)
      ax.set_xticklabels(tmp,rotation=90,verticalalignment='bottom')
    ax.set_yticks(np.arange(0,2,0.25))
    ax.axhline(1.0,color='k',linestyle='--')
    leg=ax.legend(np.array(key_legend),np.array(values_legend))#,fancybox=True, shadow=True)
    #leg.get_frame().set_alpha(0.4)
    return leg,plts

  def _onpick_curve(self,index,fig,event):
        # on the pick event, find the orig line corresponding to the
    # legend proxy line, and toggle the visibility
    legline = event.artist
    origline = self.lined[index][legline]
    for val in origline:
      vis = not val.get_visible()
      val.set_visible(vis)
    if vis:
        legline.set_alpha(1.0)
    else:
        legline.set_alpha(0.2)
    fig.canvas.draw()
  
