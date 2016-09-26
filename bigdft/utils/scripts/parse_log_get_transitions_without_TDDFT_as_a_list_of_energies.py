#!/usr/bin/env python
# -*- coding: us-ascii -*-
#----------------------------------------------------------------------------

import math
import sys
import os,copy,optparse
import cProfile
from numpy import *
import numpy as np
from Tkinter import *
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg
from matplotlib.figure import Figure
import sys
if sys.version_info[0] < 3:
    import Tkinter as Tk
else:
    import tkinter as Tk



#path of the script
path=os.path.dirname(sys.argv[0])
#yaml_folder = '/local/gigi/binaries/ifort-OMP-OCL-CUDA-gC/PyYAML-3.10/lib'
#if yaml_folder not in sys.path:
#  sys.path.insert(0,'/local/gigi/binaries/ifort-OMP-OCL-CUDA-gC/PyYAML-3.10/lib')

import yaml
#from yaml_hl import *

#import yaml_hl

class Pouet:
  def __init__(self):
    self.input = "report"
    self.output = None
    self.style = "ascii"
    self.config = os.path.join(path,'yaml_hl.cfg')

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
      if remove_it+':' in line.split('#')[0]:
         #creates a new Yaml document starting from the line
         #treat the rest of the line following the key to be removed
         header=''.join(line.split(':')[1:])
         header=header.rstrip()+'\n'
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
            except:
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

def parse_arguments():
  parser = optparse.OptionParser("This script is used to extract some information from a logfile")
  parser.add_option('-d', '--data', dest='data',default=None, #sys.argv[2],
                    help="BigDFT logfile, yaml compliant (check this if only this option is given)", metavar='FILE')
  parser.add_option('-v', '--invert-match', dest='remove',default=None, #sys.argv[2],
                    help="File containing the keys which have to be excluded from the logfile", metavar='FILE')
  parser.add_option('-e', '--extract', dest='extract',default=None, #sys.argv[2],
                    help="File containing the keys which have to be extracted to build the quantities", metavar='FILE')
  parser.add_option('-o', '--output', dest='output', default="/dev/null", #sys.argv[4],
                    help="set the output file (default: /dev/null)", metavar='FILE')
  parser.add_option('-n', '--nvirt', dest='nvirt',default=None, #sys.argv[2],
                    help="Number of virtual states to be written in the output file dos.dat (default is the number of virtual states read)", metavar='INTEGER')
  parser.add_option('-u', '--unit', dest='unit',default='eV', #sys.argv[2],
                    help="Unit of the output energies: eV or Ha", metavar='STRING (eV or Ha)')
  #Return the parsing
  return parser

def conv(res,sigma,disc):

  energy_min=res[0,0]
  energy_max=res[-1][0]
  h=(energy_max-energy_min)/disc
  grid = linspace(energy_min,energy_max,disc)
  grid.shape=(1,disc)
  l=len(res)
  a = zeros((l,disc))
  out = zeros(disc)
  tmp = zeros(l)
  for i in range(l):
    for j in range(disc):
      a[i][j] = exp(-(grid[0][j]-res[i][0])*(grid[0][j]-res[i][0])/(2*sigma*sigma))*res[i][1]
      tmp[i] += a[i][j]*h/res[i][1]
   # tmp[i]=1/tmp[i]
  TMP = eye(l)*tmp
  A = mat(a)
  A = TMP*A
  concat = ones((1,l))
  A = concat*A
  A = A[0][:]
  A=A.transpose()
  grid = grid.transpose()
  return [A,grid]

def replot(bool_plot,var_chk_button, nb_curves,sigma,sigma_old,disc,disc_old,smooth,grid,RES,a,canvas,name):
    update_to_be_plotted_curves(bool_plot,var_chk_button,nb_curves)
    calculation_new_parameters(nb_curves,sigma, sigma_old,disc,disc_old,smooth,grid,RES,bool_plot)
    plot_spectra_checkbox(a,nb_curves,smooth,grid,bool_plot,canvas,name)

def bind_return(bool_plot,var_chk_button,nb_curves,sigma, sigma_old,disc,disc_old,smooth,grid,RES,a,canvas,name):
    update_to_be_plotted_curves(bool_plot,var_chk_button,nb_curves)
    calculation_new_parameters(nb_curves,sigma, sigma_old,disc,disc_old,smooth,grid,RES,bool_plot)
    plot_spectra_checkbox(a,nb_curves,smooth,grid,bool_plot,canvas,name)
    
  
def update_to_be_plotted_curves(bool_plot,var_chk_button,nb_curves):
    for i in range(nb_curves[0]):
        bool_plot[i]=var_chk_button[i].get()
            
def calculation_new_parameters(nb_curves,sigma, sigma_old,disc,disc_old,smooth,grid,RES,bool_plot):
    for i in range(nb_curves[0]):
        if (float(sigma[i].get())!=sigma_old[i] or int(disc[i].get())!=disc_old[i]) and bool_plot[i]==1:
            smooth[i]=zeros(int(disc[i].get()))
            grid[i]=zeros(int(disc[i].get()))
            [smooth[i],grid[i]]=conv(array([RES[2*i],RES[2*i+1]]).transpose(),float(sigma[i].get()),int(disc[i].get()))  
        sigma_old[i]=float(sigma[i].get())
        disc_old[i]=int(disc[i].get())
def plot_spectra_checkbox(a,nb_curves,smooth,grid,bool_plot,canvas,name):
        a.clear()
        y_max, y_min, x_max, x_min = 0.01, 0, 0.01, 0
        for i in range(nb_curves[0]):
            if y_max<np.amax(smooth[i]):
                y_max=np.amax(smooth[i])
            if x_max<np.amax(grid[i]):
                x_max=np.amax(grid[i])
            if y_min>np.amin(smooth[i]):
                y_min=np.amin(smooth[i])
            if x_min>np.amin(grid[i]):
                x_min=np.amin(grid[i])

        y_max=y_max*1.15
        a.axis([x_min, x_max, y_min, y_max])
        a.set_title('Energy spectra')
        a.set_xlabel('Energy')
        a.set_ylabel('Extinction coefficients')
        for i in range(nb_curves[0]):
            if bool_plot[i]==1:
                a.plot(grid[i],smooth[i],label=name[i])
                a.legend(loc=2,prop={'size':6})
        canvas.show() 


def save_visible(nb_curves,bool_plot,save_plot,grid,smooth):
    for i in range(nb_curves[0]):
        if bool_plot[i]==1:
            my_file = open(save_plot[i].get(), "w")          
            for j in range(len(grid[i])):
                my_file.write(str(float(grid[i][j]))+"  "+str(float(smooth[i][j]))+"\n")
            my_file.close()
                
def add_curve(new_curve, data_frame, curves_frame,var_chk_button, nb_curves,chk_button_list,sigma_label,sigma,sigma_old,sigma_entry,disc_label,disc,disc_entry,disc_old,smooth,grid,save_label,save_plot,save_entry,RES,a,canvas,bool_plot, name):

    
    if (new_curve.get())[-4:]!="yaml":
        data=array(loadtxt(str(new_curve.get())))
        data.transpose()
        RES.append(zeros(len(data[:,0])))
        RES[2*nb_curves[0]]=data[:,0]
        RES.append(zeros(len(data[:,1])))
        RES[2*nb_curves[0]+1]=data[:,1]
        norm=np.linalg.norm(array(RES[2*nb_curves[0]+1]))
        RES[2*nb_curves[0]+1]=RES[2*nb_curves[0]+1]/norm
 
        curves_frame.append(Frame(data_frame))
        curves_frame[nb_curves[0]].pack(side='top')
        
        var_chk_button.append(IntVar())
        name.append(new_curve.get())
        chk_button_list.append(Checkbutton(curves_frame[nb_curves[0]],text=new_curve.get(),variable=var_chk_button[nb_curves[0]]))
        chk_button_list[nb_curves[0]].pack(side='left')
        bool_plot.append(var_chk_button[nb_curves[0]].get())
        
        sigma_label.append(Label(curves_frame[nb_curves[0]], text ='Sigma : '))
        sigma_label[nb_curves[0]].pack(side='left')
        sigma.append(StringVar())
        sigma[nb_curves[0]].set('0.1')
        sigma_old.append(float(sigma[nb_curves[0]].get()))
        sigma_entry.append(Entry(curves_frame[nb_curves[0]],width=6,textvariable=sigma[nb_curves[0]]))
        sigma_entry[nb_curves[0]].pack(side='left')
        
        disc_label.append(Label(curves_frame[nb_curves[0]], text ='Disc : '))
        disc_label[nb_curves[0]].pack(side='left')
        disc.append(StringVar())
        disc[nb_curves[0]].set('1400')
        disc_old.append(int(disc[nb_curves[0]].get()))
        disc_entry.append(Entry(curves_frame[nb_curves[0]],width=6,textvariable=disc[nb_curves[0]]))
        disc_entry[nb_curves[0]].pack(side='left')
        
        sigma_entry[nb_curves[0]].bind('<Return>',lambda event : bind_return(bool_plot,var_chk_button,nb_curves,sigma, sigma_old,disc,disc_old,smooth,grid,RES,a,canvas, name))
        disc_entry[nb_curves[0]].bind('<Return>',lambda event : bind_return(bool_plot,var_chk_button,nb_curves,sigma, sigma_old,disc,disc_old,smooth,grid,RES,a,canvas, name))
        
        smooth.append(zeros(int(disc[nb_curves[0]].get())))
        grid.append(zeros(int(disc[nb_curves[0]].get())))
        
        save_label.append(Label(curves_frame[nb_curves[0]],text='Name with extension : '))
        save_label[nb_curves[0]].pack(side='left')
        save_plot.append(StringVar())
        save_plot[nb_curves[0]].set('default'+str(nb_curves[0]+1)+'.dat')
        save_entry.append(Entry(curves_frame[nb_curves[0]], width=20,textvariable=save_plot[nb_curves[0]]))
        save_entry[nb_curves[0]].pack(side='left')
        
        
        [smooth[nb_curves[0]],grid[nb_curves[0]]]=conv(array([RES[2*nb_curves[0]],RES[2*nb_curves[0]+1]]).transpose(),float(sigma[nb_curves[0]].get()),int(disc[nb_curves[0]].get()))
        
        sigma_old[nb_curves[0]]=float(sigma[nb_curves[0]].get())
        disc_old[nb_curves[0]]=int(disc[nb_curves[0]].get())
        
        nb_curves[0]=nb_curves[0]+1
        update_to_be_plotted_curves(bool_plot,var_chk_button,nb_curves)
        plot_spectra_checkbox(a,nb_curves,smooth,grid,bool_plot,canvas)

    else:
        with open(new_curve.get(), "r") as fp:
            logfile_lines = fp.readlines()
        to_remove = []
        #to_extract = [yaml.load('[Excitation Energy and Oscillator Strength]', Loader = yaml.CLoader)]
        #to_extract=[[['Excitation Energy and Oscillator Strength']]]
        #to_extract = [yaml.load('[Complete list of energy eigenvalues]', Loader = yaml.CLoader)]
        to_extract=[['Complete list of energy eigenvalues']]
        cleaned_logfile=clean_logfile(logfile_lines,to_remove)
       
        datas=yaml.load_all(''.join(cleaned_logfile), Loader = yaml.CLoader)
        extracted_result=[]
        #print datas
        #blahblah
        for doc in datas:
            doc_res=document_analysis(doc,to_extract)
            if doc_res is not None:
                extracted_result.append(doc_res)
 
        #for doc in datas:
        #    doc_res=document_analysis(doc,to_extract)
        #    print doc_res
        #    if doc_res is not None:
        #        print extract_dos(doc_res)
        #        extracted_result.append(extract_dos(doc_res))
            #blahblah

        tmp=array(extracted_result[0][0])
        RES.append(zeros(len(tmp[:,0])))
        RES.append(zeros(len(tmp[:,1])))
        RES[2*nb_curves[0]]=tmp[:,0]
        RES[2*nb_curves[0]+1]=tmp[:,1]

        norm=np.linalg.norm(array(RES[2*nb_curves[0]+1]))
        RES[2*nb_curves[0]+1]=RES[2*nb_curves[0]+1]/norm
 
        curves_frame.append(Frame(data_frame))
        curves_frame[nb_curves[0]].pack(side='top')
        
        var_chk_button.append(IntVar())
        name.append(new_curve.get())
        chk_button_list.append(Checkbutton(curves_frame[nb_curves[0]],text=new_curve.get(),variable=var_chk_button[nb_curves[0]]))
        chk_button_list[nb_curves[0]].pack(side='left')
        bool_plot.append(var_chk_button[nb_curves[0]].get())
        
        sigma_label.append(Label(curves_frame[nb_curves[0]], text ='Sigma : '))
        sigma_label[nb_curves[0]].pack(side='left')
        sigma.append(StringVar())
        sigma[nb_curves[0]].set('0.1')
        sigma_old.append(float(sigma[nb_curves[0]].get()))
        sigma_entry.append(Entry(curves_frame[nb_curves[0]],width=6,textvariable=sigma[nb_curves[0]]))
        sigma_entry[nb_curves[0]].pack(side='left')
        
        disc_label.append(Label(curves_frame[nb_curves[0]], text ='Disc : '))
        disc_label[nb_curves[0]].pack(side='left')
        disc.append(StringVar())
        disc[nb_curves[0]].set('1400')
        disc_old.append(int(disc[nb_curves[0]].get()))
        disc_entry.append(Entry(curves_frame[nb_curves[0]],width=6,textvariable=disc[nb_curves[0]]))
        disc_entry[nb_curves[0]].pack(side='left')
        
        sigma_entry[nb_curves[0]].bind('<Return>',lambda event : bind_return(bool_plot,var_chk_button,nb_curves,sigma, sigma_old,disc,disc_old,smooth,grid,RES,a,canvas, name))
        disc_entry[nb_curves[0]].bind('<Return>',lambda event : bind_return(bool_plot,var_chk_button,nb_curves,sigma, sigma_old,disc,disc_old,smooth,grid,RES,a,canvas, name))
        
        smooth.append(zeros(int(disc[nb_curves[0]].get())))
        grid.append(zeros(int(disc[nb_curves[0]].get())))
        
        save_label.append(Label(curves_frame[nb_curves[0]],text='Name with extension : '))
        save_label[nb_curves[0]].pack(side='left')
        save_plot.append(StringVar())
        save_plot[nb_curves[0]].set('default'+str(nb_curves[0]+1)+'.dat')
        save_entry.append(Entry(curves_frame[nb_curves[0]], width=20,textvariable=save_plot[nb_curves[0]]))
        save_entry[nb_curves[0]].pack(side='left')
        
        
        [smooth[nb_curves[0]],grid[nb_curves[0]]]=conv(array([RES[2*nb_curves[0]],RES[2*nb_curves[0]+1]]).transpose(),float(sigma[nb_curves[0]].get()),int(disc[nb_curves[0]].get()))
        
        sigma_old[nb_curves[0]]=float(sigma[nb_curves[0]].get())
        disc_old[nb_curves[0]]=int(disc[nb_curves[0]].get())
        
        nb_curves[0]=nb_curves[0]+1
        update_to_be_plotted_curves(bool_plot,var_chk_button,nb_curves)
        plot_spectra_checkbox(a,nb_curves,smooth,grid,bool_plot,canvas,name)
        
    

def window_viz(res):
    res=array(res)
    nb_curves = [int(min(shape(res))/2)]
    RES = []
    for i in range(nb_curves[0]):
        RES.append(zeros(len(res[:,2*i])))
        RES[2*i]=res[:,2*i]
        RES.append(zeros(len(res[:,2*i+1])))
        RES[2*i+1]=res[:,2*i+1]
        norm=np.linalg.norm(array(RES[2*i+1]))
        RES[2*i+1]=RES[2*i+1]/norm
 
    y_max, y_min, x_max, x_min = 0.01, 0, 0.01, 0
    curves_frame=[]
    sigma, sigma_old, disc,disc_old, save_plot = [], [], [], [], []
    sigma_label, disc_label, save_label=[], [], []
    sigma_entry, disc_entry, save_entry= [], [], []
    save_button =[]
    var_chk_button, chk_button_list=[], []
    smooth, grid = [], []
    bool_plot=[]
    name=[]

    root = Tk.Tk()  #main window
    root.wm_title("Energy spectra")
    top=Frame(root)
    top.pack(side='top')

    data_frame = Frame(root)
    data_frame.pack(side='top')

    for i in range(nb_curves[0]):
        curves_frame.append(Frame(data_frame))
        curves_frame[i].pack(side='top')
        
        var_chk_button.append(IntVar())
        name.append('original_curve_'+str(i+1))
        chk_button_list.append(Checkbutton(curves_frame[i],text=name[i],variable=var_chk_button[i],indicatoron=1))
        chk_button_list[i].pack(side='left')
        bool_plot.append(var_chk_button[i].get())
      
        sigma_label.append(Label(curves_frame[i], text ='Sigma : '))
        sigma_label[i].pack(side='left')
        sigma.append(StringVar())
        sigma[i].set('0.1')
        sigma_old.append(float(sigma[i].get()))
        sigma_entry.append(Entry(curves_frame[i],width=6,textvariable=sigma[i]))
        sigma_entry[i].pack(side='left')
        
        disc_label.append(Label(curves_frame[i], text ='Disc : '))
        disc_label[i].pack(side='left')
        disc.append(StringVar())
        disc[i].set('1400')
        disc_old.append(int(disc[i].get()))
        disc_entry.append(Entry(curves_frame[i],width=6,textvariable=disc[i]))
        disc_entry[i].pack(side='left')

        smooth.append(zeros(int(disc[i].get())))
        grid.append(zeros(int(disc[i].get())))

        save_label.append(Label(curves_frame[i],text='Name with extension : '))
        save_label[i].pack(side='left')
        save_plot.append(StringVar())
        save_plot[i].set('default'+str(i+1)+'.dat')
        save_entry.append(Entry(curves_frame[i], width=20,textvariable=save_plot[i]))
        save_entry[i].pack(side='left')    

        


    add_frame = Frame(root)
    add_frame.pack(side='top')
    add_label=Label(add_frame,text='Add the curve : ')
    add_label.pack(side='left')
    new_curve=StringVar()
    new_curve.set('log_crmult11.yaml')
    new_curve.set('log10.yaml')
    add_entry = Entry(add_frame,width = 30,textvariable=new_curve)
    add_entry.pack(side='left')

       
        
    save_frame = Frame(root)
    save_frame.pack(side='top')
    save_button = Button(save_frame,text='save visible curves',command=lambda :save_visible(nb_curves,bool_plot,save_plot,grid,smooth))
    save_button.pack(side='left')
        

    f = plt.figure(figsize=(5,5), dpi=150)
    a = f.add_subplot(111)                  
    # create a frame for the plot
    plot_frame=Frame(root)
    plot_frame.pack(side='top')
    # a tk.DrawingArea
    canvas = FigureCanvasTkAgg(f, master=plot_frame)
    canvas.get_tk_widget().pack(side='top')
    canvas._tkcanvas.pack(side ='top')

    add_button = Button(add_frame,text='add a curve',command=lambda :add_curve(new_curve, data_frame, curves_frame,var_chk_button, nb_curves,chk_button_list,sigma_label,sigma,sigma_old,sigma_entry,disc_label,disc,disc_entry,disc_old,smooth,grid,save_label,save_plot,save_entry,RES,a,canvas,bool_plot, name))
    add_button.pack(side='left')
    replot_button = Button(save_frame,text='replot',command= lambda :replot(bool_plot,var_chk_button, nb_curves,sigma,sigma_old,disc,disc_old,smooth,grid,RES,a,canvas,name))
    replot_button.pack(side='left')

    def stop():
        root.destroy()
        root.quit()
        
    root.bind('<q>',stop)
    for i in range(nb_curves[0]):
        sigma_entry[i].bind('<Return>',lambda event : bind_return(bool_plot,var_chk_button,nb_curves,sigma, sigma_old,disc,disc_old,smooth,grid,RES,a,canvas,name))
        disc_entry[i].bind('<Return>',lambda event : bind_return(bool_plot,var_chk_button,nb_curves,sigma, sigma_old,disc,disc_old,smooth,grid,RES,a,canvas,name))

    button = Tk.Button(master = save_frame, text='Quit',command =stop)
    button.pack(side='left')
               
    # first calculations of the convolution with default values

    for i in range(nb_curves[0]):   
        [smooth[i],grid[i]]=conv(array([RES[2*i],RES[2*i+1]]).transpose(),float(sigma[i].get()),int(disc[i].get()))
   
    for i in range(nb_curves[0]):
        if y_max<np.amax(smooth[i]): 
            y_max=np.amax(smooth[i])
        if x_max<np.amax(grid[i]):
            x_max=np.amax(grid[i])
        if y_min>np.amin(smooth[i]):
            y_min=np.amin(smooth[i])
        if x_min>np.amin(grid[i]):
            x_min=np.amin(grid[i])

    y_max=y_max*1.15
    a.axis([x_min, x_max, y_min, y_max])
    a.set_title('Energy spectra')
    a.set_xlabel('Energy')
    a.set_ylabel('Extinction coefficients')
    
    for i in range(nb_curves[0]):               
        a.plot(grid[i],smooth[i],label=name[i])
        a.legend(loc=2,prop={'size':6})


    canvas.show()
    canvas.get_tk_widget().pack(side='top')
    canvas._tkcanvas.pack(side ='top')

    
    Tk.mainloop()
    return


def extract_dos(doc2):
    ev=[]
    for i in doc2:
        if not i.has_key('HOMO LUMO gap (AU, eV)'):
            ev.append(i.values() + [1.0])
    return ev
    
    
if __name__ == "__main__":
  parser = parse_arguments()
  (args, argtmp) = parser.parse_args()


#args=parse_arguments()
#logfile
print args.data
with open(args.data, "r") as fp:
  logfile_lines = fp.readlines()
#output file
file_out=open(args.output, "w")
#to_remove list
if args.remove is not None:
  to_remove = yaml.load(open(args.remove, "r").read(), Loader = yaml.CLoader)
else:
  #standard list which removes long items from the logfile
  # to_remove=["Atomic positions within the cell (Atomic and Grid Units)",
  #            "Atomic Forces (Ha/Bohr)",
  #            "Orbitals",
  #            "Energies",
  #            "Properties of atoms in the system"]
  to_remove=[]

    
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
#if args.extract is not None:
#  to_extract = yaml.load(open(args.extract, "r").read(), Loader = yaml.CLoader)
#  print to_extract
#else:
# to_extract=[[['Excitation Energy and Oscillator Strength']]]
# to_extract.append([['Complete list of energy eigenvalues']])
# to_extract.append([['Transition energies (eV)']])

#Extract the eigenvalues of the calculation            
to_extract=[[['Complete list of energy eigenvalues']]]
datas=yaml.load_all(''.join(cleaned_logfile), Loader = yaml.CLoader)
extracted_result=[]
for doc in datas:
    for i in range(len(to_extract)):
        doc_res=document_analysis(doc,to_extract[i])
        if doc_res is not None:
            extracted_result.append(doc_res)

#print "Number of valid documents:",len(extracted_result)
allenergies=extracted_result[0][0] #store all the information concerning the eigenvalues of the BigDFT calculation
##print ''
##print 'to_extract',to_extract
###print 'it=', it
##print 'allenergies=', allenergies
##print 'len(allenergies)', len(allenergies)
#print 'range(len(it[0]))', range(len(it[0]))
enocc=[] #initialize the list of the energies of the occupied states
envirt=[] #initialize the list of the energies of unoccupied states
#Get the occupied and virtual energies
for energy in allenergies:
    #print ''
    #print 'ien:', ien, 'allenergies[ien]:', allenergies[ien], 'type:', type(allenergies[ien]), 'len(allenergies[ien])', len(allenergies[ien])
    if energy.has_key('e_occ'): #if it is an occupied state energy
        #print 'new enocc:', energy['e_occupied']
        enocc.append(energy['e_occ']) #then add it to the list of occupied states energies
        #print 'enocc:', enocc
    elif energy.has_key('e_vrt'): #if it is an unoccupied state energy
        #print 'new envirt:', energy['e_virtual']
        envirt.append(energy['e_vrt']) #then add it to the list of unoccupied states energies
        #print 'envirt:', envirt

#define the number of vitrual states to be written.
if args.nvirt == None: 
   args.nvirt = int(len(envirt))
else:
   print int(args.nvirt), len(envirt)
   if int(args.nvirt) > len(envirt):
       print 'ERROR: too many virtual states asked'
       exit()
   else:
       args.nvirt = int(args.nvirt)
       #take only the first args.nvirt unoccupied states energies
       envirt = envirt[:args.nvirt]

print '\noccupied states energies', enocc, 'length:', len(enocc),'\n'
print 'virtual states energies', envirt, 'length:', len(envirt),'\n'


#Compute all the energy differences between occupied and unoccupied states
e_diff = [] #initialize the variables that will contain all differences between occupied and unoccupied states
for e_o in enocc: #loop over all occupied states energies
    for e_v in envirt: #loop over all unoccupied states energies
        e_diff.append(e_v - e_o) #add the difference to e_diff
e_diff.sort() #sort the list
e_diff = np.asarray(e_diff) #transform it to a numpy array
print 'energy differences: \n', e_diff #print the numpy array

#define the unit of the output energies
if args.unit == 'eV': e_diff = 27.211396132 * e_diff #convert from Ha to eV

#write an output file
#This output file is made of a list of energy differences, from the lowest to the highest
write_output = True

if write_output:

    print('Writing output\n')

    #open the output file
    mydosfile_e = open('transitions_noTDDFT_'+str(args.nvirt).strip()+'v.dat','w')

    #write energy differences
    for e_d in e_diff:
        mydosfile_e.write(  ' %.16e \n' % ( e_d ) )
 
    #close file
    mydosfile_e.close()


##
###Extract peak energies and transitions associated to it
##print ''
##to_extract=[[['Transition energies (eV)']]]
##datas=yaml.load_all(''.join(cleaned_logfile), Loader = yaml.CLoader)
##extracted_result=[]
##for doc in datas:
##    for i in range(len(to_extract)):
##        doc_res=document_analysis(doc,to_extract[i])
##        if doc_res is not None:
##            extracted_result.append(doc_res)
##
###print "Number of valid documents:",len(extracted_result)
##print ''
##allpeaks=extracted_result[0][0] #store all the information concerning the peaks in TDDFT spectra and of transitions involved in these peaks
##print 'to_extract',to_extract
##print 'allpeaks=', allpeaks
##print 'len(allpeaks)=', len(allpeaks)
##peak_en=[] #initialize the peak energies
##j=0 #initialize the counter of the number of transitions
##sp="   " #space needed to write the lines in the output files
##
###Get the energy of the peaks and the transitions involved in these peaks.
###At the same time, define the lines to be written in the file 'transitions.txt'.
##for ipeak,peak in enumerate(allpeaks): #loop over all the peaks
##    #print ''
##    #print 'ipeak:', ipeak, 'peak:', peak, 'type:', type(peak)
##    if (ipeak==0): #if first peak then initialize the transition and coefficient list, as well as the list of lines to write in 'transitions.txt
##        #print 'initialize the transition and coefficient list, as well as the list of lines to write in 'transitions.txt'.
##        trans=[[]]
##        coeff=[[]]
##        linestowrite=[]
##    else: #else add a new transition and coefficient list
##        #print 'add a new transition and coefficient list'
##        trans.append([])
##        coeff.append([])
##    for itrans,alltrans in enumerate(peak): #loop over all the transitions involved in the studied peak
##        #print ''
##        #print 'itrans', itrans, '(len(peak)-1)/2', (len(peak)-1)/2
##        if (itrans>(len(peak)-1)/2): #at the moment, half of the transitions are redundant (due to TDA ?), so we can skip them
##            break
##        elif(alltrans.has_key('Energy')): #else, if data concern the energy,
##            #print 'new peak_en:', alltrans['Energy']
##            peak_en.append(alltrans['Energy']) #then memorize it
##            linestowrite.append(str(alltrans['Energy'])) #and write it the output line
##            #print ipeak, itrans, linestowrite
##            #print 'peak_en', peak_en
##        elif(alltrans.has_key('Transition')): #else if data concern the transitions,
##            #print 'ipeak', ipeak, 'itrans', itrans,'new transition', alltrans['Transition'], 'old transitions:', trans[ipeak]
##            trans[ipeak].append(alltrans['Transition']) #then memorize it
##            coeff[ipeak].append(alltrans['Coeff']) #as well as the coefficient associated to it.
##            iocc=alltrans['Transition'][0] #define the occupied state involved in the transition.
##            ivirt=alltrans['Transition'][1] #define the virtual state involved in the transition.
##            if (itrans==1): #if it is the first transition associated to the peak
##                linestowrite[j]=linestowrite[j]+sp+str(iocc)+sp+str(enocc[iocc-1])+sp+str(ivirt)+sp+str(envirt[ivirt-1])+sp+str(alltrans['Coeff'])+'\n' #add information to the line to write
##            elif (itrans>1): #else, if it is not the first transitions
##                linestowrite.append(str(peak_en[ipeak])+sp+str(iocc)+sp+str(enocc[iocc-1])+sp+str(ivirt)+sp+str(envirt[ivirt-1])+sp+str(alltrans['Coeff'])+'\n') #define a new line to write in
##            #print ipeak, itrans, linestowrite
##            j=j+1 #increment the counter of transitions
##
##
##print ''
##print 'peaks energy', peak_en, 'length:', len(peak_en)
##print ''
##print 'transitions', trans, 'length:', len(trans)
##print ''
##print 'coeff', coeff, 'length:', len(coeff)
##
##
##print ''
##print 'writing an output file containing transitions'
##print 'total number of transitions:', j
##
##myfile=open('transitions.txt','w') #open an output file to write
##myfile.write(str(j)+'\n') #write the number of transitions in the first line
##for line in linestowrite: #loop over the lines to write
##    myfile.write(line) #write the line in my file
##myfile.close() #close the file
##
##
##sys.exit(0)
##
###do some tests
##document = """---
##block sequence:
## - BlockEntryToken
## block mapping:
##   ? KeyToken
##   : ValueToken
## flow sequence: [FlowEntryToken, FlowEntryToken]
## flow mapping: {KeyToken: ValueToken}
## anchors and tags:
## - &A !!int '5'
## - *A
##...
## """
##
###for token in yaml.scan(document):
###     print token
##
###ddd
###print args.ref,args.data,args.output
##
##datas    = [a for a in yaml.load_all(open(args.data, "r").read(), Loader = yaml.CLoader)]
###i=0
###for doc in yaml.load_all(open(args.data, "r").read(), Loader = yaml.CLoader):
###  i+=1
###  print 'doc No.',i#,"Last Iteration" in doc.keys()
##  #try:
##  #  print yaml.dump(doc["Last Iteration"])
##  #except:
##  #  print "Last Iteration not found"
##
###datas    = [a for a in yaml.safe_load_all(open(args.data, "r").read())]
###Profile.run('datas    = [a for a in yaml.load_all(open(args.data, "r").read(), Loader = yaml.CLoader)]')
###gyi
###runfile = open(args.data, "r").read()
##sss
##try:
##    datas    = [a for a in yaml.load_all(runfile, Loader = yaml.Loader)]
##    #Profile.run('datas    = [a for a in yaml.load_all(open(args.data, "r").read())]')
###for a in yaml.load_all(open(args.data, "r").read(), Loader = yaml.CLoader):
###    print 'New Document'
###    document_analysis(a)
##  
##except:
##  datas = []
##  reports = open(args.output, "w")
##  fatal_error(args,reports)
##
##ndocs=len(datas)
##
##print 'No. of Documents:',ndocs
##
##for run in datas:
##  print 'New Document'
##  document_analysis(run)
##
##
##
##print 'end'
##
##test = run['Atomic positions within the cell (Atomic and Grid Units)']
##
###value = test[1]["Ag"]["AU"][1]
##
###ex={}
###ex["value"]=complex(12,13)
###print 'example',value,3+4j,value,(3+4j),complex(7,8),ex["value"]+complex(value),test[1]["Ag"]["AU"][0]+test[1]["Ag"]["AU"][2]
##
##sys.stdout.write(yaml.dump(run["Last Iteration"],default_flow_style=False,explicit_start=True))
###sys.stdout.write(yaml.dump(run["Direct and transposed data repartition"],default_flow_style=False,explicit_start=True))
###sys.stdout.write(yaml.dump(ex,default_flow_style=False,explicit_start=True))
##
##sys.exit(0)


