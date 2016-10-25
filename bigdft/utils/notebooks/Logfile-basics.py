
# coding: utf-8

# # Analyze a BigDFT run

# In this example we will inspect how the BigDFT python library can be used to retrieve the results from a run.
# 
# As a first operation we have to be sure that the path for the BigDFT scripts is present in the PYTHONPATH.
# If not so, this might be updated by the script in the following way:

# In[1]:

#write here the position of the BigDFT installation (usually <buildtree>/install/lib/pythonX.X/site-packages)
import os
BIGDFT_PYTHONDIR=os.path.abspath(
    os.path.join(os.pardir,os.pardir,'src','python')
    ) #refer to the sources, only for test 
#then update the path
import sys
if BIGDFT_PYTHONDIR not in sys.path: 
    sys.path+=[BIGDFT_PYTHONDIR]


# Now we can load the Logfiles module:

# In[2]:

from BigDFT import Logfiles as lf 


# Let us now load a file into a instance of a Logfile class. Imagine that our logfile corresponds to a single-point run, and it is present in a file named "log-HBDMI.yaml":
# 

# In[3]:

HBDMI=lf.Logfile('testfiles/log-HBDMI.yaml')


# The run is now loaded. To inspect its behaviour we might print it down to inspect the usual information:

# In[4]:

print HBDMI
HBDMI.geopt_plot()


# The above information can also be accessed separately, by having a look at the attributes of the HBDMI object:

# In[5]:

HBDMI.__dict__.keys() #you may also type: dir(HBDMI)


# For this reason, we might consider to postprocess some of the variables for study of the system. Here an example:

# In[6]:

print 'The average energy per atom is:',HBDMI.energy/HBDMI.nat,    '(',HBDMI.energy,' Ha),',HBDMI.nat,' atoms)'
print 'There are also,',HBDMI.evals[0].info,' (up,down) orbitals in the run'


# We might access the Density of States of this system (note that if the matplotlib inlining is deactivated, the RadioButton widget below allows to adjust the smearing):

# In[7]:

HBDMIDos=HBDMI.get_dos(label='HBDMI molecule')
HBDMIDos.plot(sigma=0.2)


# ## Case of a periodic system

# The above case was a Free BC molecule single point run. Let us now consider the case of a periodic calculation.
# We take as an example a logfile coming from one run of the DeltaTest benchmark (see [this page](https://molmod.ugent.be/deltacodesdft) to know what it is all about). In any case, let us load the "log-K-1.0.yaml" file:

# In[8]:

K=lf.Logfile('testfiles/log-K.yaml')
print K


# Here we can see that there are also other attributes available, like the k-points and the pressure (in GPa):

# In[9]:

dir(K)


# Here we might trace the density of states but also the band structure plot, in a similar fashion (sse also the Band Structure notebook):

# In[10]:

K.get_dos().plot()
K.get_brillouin_zone().plot()


# As another example for the system we might inspect the kpoints:

# In[11]:

print K.kpts[0]
print K.kpts[-1]


# ## Case of a Geometry optimization

# For a geometry optimization the situation is similar, with the extra point that the code automatically recognize  multiple runs inside the logfile. Let us see the example of the following logfile:

# In[12]:

geopt=lf.Logfile('testfiles/GEOPT-all_sqnmbiomode.out.ref.yaml')
print geopt


# The interesting point is that now the logfile can be iterated among the different geometry steps:

# In[13]:

en=[l.energy for l in geopt]
for i,e in enumerate(en):
    print i,e


# The geopt_plot() function allows to plot the relation beween energy and forces, where it can be also seen that the desired criterion is reached. Errorbars show the local fluctuation of the forces, an indication of the (cleaned) center of mass drift. See the example:

# In[14]:

geopt.geopt_plot()

