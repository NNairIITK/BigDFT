
# coding: utf-8

# # Band Structure Visualization

# This notebook presents some of the basics features that can be used to visualize Band Structure from a BigDFT logfile.
# The only input which is needed is the logfile of the system, wuth a k-point meshi that is dense enough to achieve convergence of the result. The bans structure is then obtained via interpolation over the Brillouin Zone of the mesh defined in this way.
# To use this notebook you must have installed in your system the [spglib](http://atztogo.github.io/spglib/python-spglib.html#python-spglib) and the [ASE](http://atztogo.github.io/spglib/python-spglib.html#python-spglib) python library, whose installation is really standard and can be done in few minutes.
# 
# Let us start by setting up the PYTHONPATH to the BigDFT installation directory.
# Here we point to the sources as we do not know where the code is installed, but the directory below should be <build>/install/lib/python2.7/site-packages. However:

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


# Next we load a logfile which is suitable for the calculation, for example an evergreen (this log has been shortened for practical purposes):

# In[2]:

import BigDFT.Logfiles as lf
#get_ipython().magic(u'matplotlib inline')
Graph=lf.Logfile('testfiles/log-Graphene.yaml')


# In[3]:

#inform
print Graph


# In[4]:

#How many k-points?
Graph.nkpt


# To start the analysis, it might be interesting to plot the Density of States (also see the DoS example notebook):

# In[5]:

Graph.get_dos(label='C').plot()


# Where we can see that the density of states at the fermi level is zero as it should.
# Let us now instance the Brillouin Zone class:

# In[6]:

GBZ=Graph.get_brillouin_zone()


# The inspection of the symmetry of the system provides us with quantities that we might compare with the log properties above. The BrillouinZone class has also different attributes. In particulat we might be interested in kwowing the coorinates of the special points or the special paths of the corresponding structure:

# In[7]:

dir(GBZ)


# In[8]:

#let us print the special points
import yaml
print yaml.dump(GBZ.special_points)


# In[9]:

#and the special paths
print yaml.dump(GBZ.special_paths)


# Let us now plot the band structure. Without specifying anything the plot is performed on the first special path, wrappend again at Gamma. We might also choose a fine discretization for the interpolation:

# In[10]:

GBZ.plot(npts=5000)


# This path does not really makes sense for graphene as it also involves points which have nonzero coordinate in the second index (which is a isolated direction according to BigDFT conventions)
# We might of course also choose other paths:

# In[11]:

import BigDFT.BZ as BZ
reload(BZ)
path1=BZ.BZPath(GBZ.lattice,GBZ.special_paths[2],GBZ.special_points,npts=1000)
GBZ.plot(path=path1)


# Or we might remove the points from the first path which have nonzero second coodinate:

# In[15]:

path0=[]
for p in GBZ.special_paths[0]+['G']:
    if GBZ.special_points[p][1]==0.0: path0.append(p)
print path0
path0=BZ.BZPath(GBZ.lattice,path0,GBZ.special_points,npts=1000)
GBZ.plot(path=path0)


# Where we can see the K and K' points with zero gap.
# 
# Let us now construct a Customized path that passes through the point of zero gap.
# But first, we should inspect the k-point eigenvalues to find the k-point coordinates which correspond to this:

# In[149]:

import numpy
print numpy.argmin([numpy.min(abs(e-Graph.fermi_level)) for e in Graph.evals])


# It seems that it is the k-point number 12 (or 13 in fortran index). Let us verify:

# In[16]:

Kprime=Graph.evals[12]
print Kprime,Graph.fermi_level


# Indeed it is so. Let us now check which is its coordinate. There are attributes of the BandArray instance Kprime which might be of interest to us:

# In[17]:

#brilloiun zone coordinate
print Kprime.kpt
#kpoint weight
print Kprime.kwgt


# Anyhow, let us create a path that passes through this point. As it is not among the list of special point (we provided Orthorhombic coordinates for the Graphene), we should create a path via ASE. This is automatically done by the BZ instance for special points.

# In[18]:

path2=BZ.BZPath(GBZ.lattice,path=['G',{'Kp':Kprime.kpt},'X','G'],special_points=GBZ.special_points,npts=1000)


# In[19]:

GBZ.plot(path=path2)


# ## Case of a 3d periodic (cubic) system

# The same thing can be done for a cubic system, like the one already presented in the Logfile analysis example.
# The following sections should be straightforward.

# In[20]:

K=lf.Logfile('testfiles/log-K.yaml')


# In[21]:

BZK=K.get_brillouin_zone()


# In[22]:

BZK.plot()

