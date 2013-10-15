#!python 

from numpy import *
import numpy
import string
import sys



#/*##########################################################################
# Copyright (C) 2004-2007 European Synchrotron Radiation Facility
#
# This file is part of the PyMCA X-ray Fluorescence Toolkit developed at
# the ESRF by the Beamline Instrumentation Software Support (BLISS) group.
#
# This toolkit is free software; you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the Free
# Software Foundation; either version 2 of the License, or (at your option)
# any later version.
#
# PyMCA is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
# details.
#
# You should have received a copy of the GNU General Public License along with
# PyMCA; if not, write to the Free Software Foundation, Inc., 59 Temple Place,
# Suite 330, Boston, MA 02111-1307, USA.
#
# PyMCA follows the dual licensing model of Trolltech's Qt and Riverbank's PyQt
# and cannot be used as a free plugin for a non-free program.
#
# Please contact the ESRF industrial unit (industry@esrf.fr) if this license
# is a problem for you.
#############################################################################*/

import random
import math
import numpy
from numpy.oldnumeric import *
from numpy.oldnumeric.linear_algebra import inverse
from numpy.linalg import svd
import time
__author__ = "V.A. Sole <sole@esrf.fr>"
__revision__ = "$Revision: 1.19 $"
# codes understood by the routine
CFREE       = 0
CPOSITIVE   = 1
CQUOTED     = 2
CFIXED      = 3
CFACTOR     = 4
CDELTA      = 5
CSUM        = 6
CIGNORED    = 7

ONED = 0

def LeastSquaresFit(model, parameters0, data=None, maxiter = 50,constrains=[],
                        weightflag = 0,model_deriv=None,deltachi=None,fulloutput=0,
                        xdata=None,ydata=None,sigmadata=None,linear=None):
    parameters = array(parameters0).astype(Float)
    if linear is None:linear=0
    if deltachi is None:
        deltachi = 0.01
    if ONED:
      data0 = array(data)
      x = data0[0:2,0]
    #import SimplePlot
    #SimplePlot.plot([data[:,0],data[:,1]],yname='Received data')
    else:
        if xdata is None:
            x=array(map(lambda y:y[0],data))
        else:
            x=xdata
    if linear:
           return LinearLeastSquaresFit(model,parameters0,
                                        data,maxiter,
                                        constrains,weightflag,model_deriv=model_deriv,
                                        deltachi=deltachi,
                                        fulloutput=fulloutput,
                                        xdata=xdata,
                                        ydata=ydata,
                                        sigmadata=sigmadata)
    elif len(constrains) == 0:
        try:
            model(parameters,x)
            constrains = [[],[],[]]
            for i in range(len(parameters0)):
                constrains[0].append(0)
                constrains[1].append(0)
                constrains[2].append(0)
            return RestreinedLeastSquaresFit(model,parameters0,
                                    data,maxiter,
                                    constrains,weightflag,
                                    model_deriv=model_deriv,
                                    deltachi=deltachi,
                                    fulloutput=fulloutput,
                                    xdata=xdata,
                                    ydata=ydata,
                                    sigmadata=sigmadata)
        except TypeError:
            print "You should reconsider how to write your function"
            raise TypeError
    else:
            return RestreinedLeastSquaresFit(model,parameters0,
                                    data,maxiter,
                                    constrains,weightflag,model_deriv=model_deriv,
                                    deltachi=deltachi,
                                    fulloutput=fulloutput,
                                    xdata=xdata,
                                    ydata=ydata,
                                    sigmadata=sigmadata)

def LinearLeastSquaresFit(model0,parameters0,data0,maxiter,
                                constrains0,weightflag,model_deriv=None,deltachi=0.01,fulloutput=0,
                                    xdata=None,
                                    ydata=None,
                                    sigmadata=None):
    #get the codes:
    # 0 = Free       1 = Positive     2 = Quoted
    # 3 = Fixed      4 = Factor       5 = Delta
    # 6 = Sum        7 = ignored
    constrains = [[],[],[]]
    if len(constrains0) == 0:
        for i in range(len(parameters0)):
            constrains[0].append(0)
            constrains[1].append(0)
            constrains[2].append(0)
    else:
        for i in range(len(parameters0)):
            constrains[0].append(constrains0[0][i])
            constrains[1].append(constrains0[1][i])
            constrains[2].append(constrains0[2][i])
    for i in range(len(parameters0)):
        if type(constrains[0][i]) == type('string'):
            #get the number
            if   constrains[0][i] == "FREE":
                 constrains[0][i] = CFREE
            elif constrains[0][i] == "POSITIVE":
                 constrains[0][i] = CPOSITIVE
            elif constrains[0][i] == "QUOTED":
                 constrains[0][i] = CQUOTED
            elif constrains[0][i] == "FIXED":
                 constrains[0][i] = CFIXED
            elif constrains[0][i] == "FACTOR":
                 constrains[0][i] = CFACTOR
                 constrains[1][i] = int(constrains[1][i])
            elif constrains[0][i] == "DELTA":
                 constrains[0][i] = CDELTA
                 constrains[1][i] = int(constrains[1][i])
            elif constrains[0][i] == "SUM":
                 constrains[0][i] = CSUM
                 constrains[1][i] = int(constrains[1][i])
            elif constrains[0][i] == "IGNORED":
                 constrains[0][i] = CIGNORED
            elif constrains[0][i] == "IGNORE":
                 constrains[0][i] = CIGNORED
            else:
               #I should raise an exception
                #constrains[0][i] = 0
                raise ValueError,"Unknown constraint %s" % constrains[0][i]
        if (constrains[0][i] == CQUOTED):
            raise ValueError, "Linear fit cannot handle quoted constraint"
    # make a local copy of the function for an easy speed up ...
    model = model0
    parameters = array(parameters0)
    if data0 is not None:
        selfx = array(map(lambda x:x[0],data0))
        selfy = array(map(lambda x:x[1],data0))
    else:
        selfx = xdata
        selfy = ydata
    selfweight = ones(selfy.shape,Float)
    nr0 = len(selfy)
    if data0 is not None:
        nc =  len(data0[0])
    else:
        if sigmadata is None:
            nc = 2
        else:
            nc = 3
    if weightflag == 1:
        if nc == 3:
            #dummy = abs(data[0:nr0:inc,2])
            if data0 is not None:
                dummy = abs(array(map(lambda x:x[2],data0)))
            else:
                dummy = abs(array(sigmadata))
            selfweight = 1.0 / (dummy + equal(dummy,0))
            selfweight = selfweight * selfweight
        else:
            selfweight = 1.0 / (abs(selfy) + equal(abs(selfy),0))
    n_param = len(parameters)
    #linear fit, use at own risk since there is no check for the
    #function being linear on its parameters.
    #Only the fixed constrains are handled properly
    x=selfx
    y=selfy
    weight = selfweight
    iter  = maxiter
    niter = 0
    newpar = parameters.__copy__()
    while (iter>0):
        niter+=1
        chisq0, alpha0, beta,\
        n_free, free_index, noigno, fitparam, derivfactor  =ChisqAlphaBeta(
                                                 model,newpar,
                                                 x,y,weight,constrains,model_deriv=model_deriv,
                                                 linear=1)
        print "A", chisq0

        nr, nc = alpha0.shape
        fittedpar = dot(beta, inverse(alpha0))
        #check respect of constraints (only positive is handled -force parameter to 0 and fix it-)
        error = 0
        for i in range(n_free):
            if constrains [0] [free_index[i]] == CPOSITIVE:
                if fittedpar[0,i] < 0:
                    #fix parameter to 0.0 and re-start the fit
                    newpar[free_index[i]] = 0.0
                    constrains[0][free_index[i]] = CFIXED
                    error = 1
        if error:
            continue
        for i in range(n_free):
            newpar[free_index[i]] = fittedpar[0,i]
        newpar=array(getparameters(newpar,constrains))
        iter=-1
    yfit = model(newpar,x)
    chisq = sum( weight * (y-yfit) * (y-yfit))
    sigma0 = sqrt(abs(diagonal(inverse(alpha0))))
    sigmapar = getsigmaparameters(newpar,sigma0,constrains)
    lastdeltachi = chisq
    if not fulloutput:
        return newpar.tolist(), chisq/(len(y)-len(sigma0)), sigmapar.tolist()
    else:
        return newpar.tolist(), chisq/(len(y)-len(sigma0)), sigmapar.tolist(),niter,lastdeltachi

def RestreinedLeastSquaresFit(model0,parameters0,data0,maxiter,
                constrains0,weightflag,model_deriv=None,deltachi=0.01,fulloutput=0,
                                    xdata=None,
                                    ydata=None,
                                    sigmadata=None):
    #get the codes:
    # 0 = Free       1 = Positive     2 = Quoted
    # 3 = Fixed      4 = Factor       5 = Delta
    # 6 = Sum        7 = ignored
    constrains=[[],[],[]]
    for i in range(len(parameters0)):
        constrains[0].append(constrains0[0][i])
        constrains[1].append(constrains0[1][i])
        constrains[2].append(constrains0[2][i])
    for i in range(len(parameters0)):
        if type(constrains[0][i]) == type('string'):
            #get the number
            if   constrains[0][i] == "FREE":
                 constrains[0][i] = CFREE
            elif constrains[0][i] == "POSITIVE":
                 constrains[0][i] = CPOSITIVE
            elif constrains[0][i] == "QUOTED":
                 constrains[0][i] = CQUOTED
            elif constrains[0][i] == "FIXED":
                 constrains[0][i] = CFIXED
            elif constrains[0][i] == "FACTOR":
                 constrains[0][i] = CFACTOR
                 constrains[1][i] = int(constrains[1][i])
            elif constrains[0][i] == "DELTA":
                 constrains[0][i] = CDELTA
                 constrains[1][i] = int(constrains[1][i])
            elif constrains[0][i] == "SUM":
                 constrains[0][i] = CSUM
                 constrains[1][i] = int(constrains[1][i])
            elif constrains[0][i] == "IGNORED":
                 constrains[0][i] = CIGNORED
            elif constrains[0][i] == "IGNORE":
                 constrains[0][i] = CIGNORED
            else:
               #I should raise an exception
                #constrains[0][i] = 0
                raise ValueError,"Unknown constraint %s" % constrains[0][i]
    # make a local copy of the function for an easy speed up ...
    model = model0
    parameters = array(parameters0)
    if ONED:
        data = array(data0)
        x = data[1:2,0]
    fittedpar = parameters.__copy__()
    flambda = 0.001
    iter = maxiter
    niter = 0
    if ONED:
        selfx = data [:,0]
        selfy = data [:,1]
    else:
        if data0 is not None:
            selfx = array(map(lambda x:x[0],data0))
            selfy = array(map(lambda x:x[1],data0))
        else:
            selfx = xdata
            selfy = ydata
    selfweight = ones(selfy.shape,Float)
    if ONED:
        nr0, nc = data.shape
    else:
        nr0 = len(selfy)
        if data0 is not None:
            nc =  len(data0[0])
        else:
            if sigmadata is None:
                nc = 2
            else:
                nc = 3

    if weightflag == 1:
            if nc == 3:
                #dummy = abs(data[0:nr0:inc,2])
                if ONED:
                    dummy = abs(data [:,2])
                else:
                    if data0 is not None:
                        dummy = abs(array(map(lambda x:x[2],data0)))
                    else:
                        dummy = abs(array(sigmadata))
                selfweight = 1.0 / (dummy + equal(dummy,0))
                selfweight = selfweight * selfweight
            else:
                selfweight = 1.0 / (abs(selfy) + equal(abs(selfy),0))
    n_param = len(parameters)
    selfalphazeros = zeros((n_param, n_param),Float)
    selfbetazeros = zeros((1,n_param),Float)
    index = arange(0,nr0,1)
    while (iter > 0):
        niter = niter + 1
        if (niter < 2) and (n_param*3 < nr0):
                x=take(selfx,index)
                y=take(selfy,index)
                weight=take(selfweight,index)
        else:
                x=selfx
                y=selfy
                weight = selfweight

        chisq0, alpha0, beta,\
        n_free, free_index, noigno, fitparam, derivfactor  =ChisqAlphaBeta(
                                                 model,fittedpar,
                                                 x,y,weight,constrains,model_deriv=model_deriv)
        print"B ",  chisq0


        nr, nc = alpha0.shape
        flag = 0
        lastdeltachi = chisq0
        while flag == 0:
            newpar = parameters.__copy__()
            if(1):
                alpha = alpha0 + flambda * identity(nr) * alpha0
                deltapar = dot(beta, inverse(alpha))
            else:
                #an attempt to increase accuracy
                #(it was unsuccessful)
                alphadiag=sqrt(diagonal(alpha0))
                npar = len(sqrt(diagonal(alpha0)))
                narray = zeros((npar,npar),Float)
                for i in range(npar):
                    for j in range(npar):
                        narray[i,j] = alpha0[i,j]/(alphadiag[i]*alphadiag[j])
                narray = inverse(narray + flambda * identity(nr))
                for i in range(npar):
                    for j in range(npar):
                        narray[i,j] = narray[i,j]/(alphadiag[i]*alphadiag[j])
                deltapar = dot(beta, narray)
            pwork = zeros(deltapar.shape, Float)
            for i in range(n_free):
                if constrains [0] [free_index[i]] == CFREE:
                    pwork [0] [i] = fitparam [i] + deltapar [0] [i]
                elif constrains [0] [free_index[i]] == CPOSITIVE:
                    #abs method
                    pwork [0] [i] = fitparam [i] + deltapar [0] [i]
                    #square method
                    #pwork [0] [i] = (sqrt(fitparam [i]) + deltapar [0] [i]) * \
                    #                (sqrt(fitparam [i]) + deltapar [0] [i])
                elif constrains [0] [free_index[i]] == CQUOTED:
                    pmax=max(constrains[1] [free_index[i]],
                            constrains[2] [free_index[i]])
                    pmin=min(constrains[1] [free_index[i]],
                            constrains[2] [free_index[i]])
                    A = 0.5 * (pmax + pmin)
                    B = 0.5 * (pmax - pmin)
                    if (B != 0):
                        pwork [0] [i] = A + \
                                    B * sin(arcsin((fitparam[i] - A)/B)+ \
                                    deltapar [0] [i])
                    else:
                        print "Error processing constrained fit"
                        print "Parameter limits are",pmin,' and ',pmax
                        print "A = ",A,"B = ",B
                newpar [free_index[i]] = pwork [0] [i]
            newpar=array(getparameters(newpar,constrains))
            workpar = take(newpar,noigno)
            #yfit = model(workpar.tolist(), x)
            yfit = model(workpar,x)
            chisq = sum( weight * (y-yfit) * (y-yfit))
            print "chisq ", chisq, "chisq0 ", chisq0
            if chisq > chisq0:
                flambda = flambda * 10.0
                if flambda > 1000:
                    flag = 1
                    iter = 0
            else:
                flag = 1
                fittedpar = newpar.__copy__()
                lastdeltachi = (chisq0-chisq)/(chisq0+(chisq0==0))
                if (lastdeltachi) < deltachi:
                    pass
                    # iter = 0
                chisq0 = chisq
                flambda = flambda / 10.0
                print "iter = ",iter,"chisq = ", chisq
            iter = iter -1
    sigma0 = sqrt(abs(diagonal(inverse(alpha0))))
    sigmapar = getsigmaparameters(fittedpar,sigma0,constrains)
    if not fulloutput:
        return fittedpar.tolist(), chisq/(len(yfit)-len(sigma0)), sigmapar.tolist()
    else:
        return fittedpar.tolist(), chisq/(len(yfit)-len(sigma0)), sigmapar.tolist(),niter,lastdeltachi

def ChisqAlphaBeta(model0, parameters, x,y,weight, constrains,model_deriv=None,linear=None):
    if linear is None:linear=0
    model = model0
    #nr0, nc = data.shape
    n_param = len(parameters)
    n_free = 0
    fitparam=[]
    free_index=[]
    noigno = []
    derivfactor = []
    for i in range(n_param):
        if constrains[0] [i] != CIGNORED:
            noigno.append(i)
        if constrains[0] [i] == CFREE:
            fitparam.append(parameters [i])
            derivfactor.append(1.0)
            free_index.append(i)
            n_free += 1
        elif constrains[0] [i] == CPOSITIVE:
            fitparam.append(abs(parameters[i]))
            derivfactor.append(1.0)
            #fitparam.append(sqrt(abs(parameters[i])))
            #derivfactor.append(2.0*sqrt(abs(parameters[i])))
            free_index.append(i)
            n_free += 1
        elif constrains[0] [i] == CQUOTED:
            pmax=max(constrains[1] [i],constrains[2] [i])
            pmin=min(constrains[1] [i],constrains[2] [i])
            if ((pmax-pmin) > 0) & \
               (parameters[i] <= pmax) & \
               (parameters[i] >= pmin):
                A = 0.5 * (pmax + pmin)
                B = 0.5 * (pmax - pmin)
                if 1:
                    fitparam.append(parameters[i])
                    derivfactor.append(B*cos(arcsin((parameters[i] - A)/B)))
                else:
                    help0 = arcsin((parameters[i] - A)/B)
                    fitparam.append(help0)
                    derivfactor.append(B*cos(help0))
                free_index.append(i)
                n_free += 1
    fitparam = array(fitparam, Float)
    alpha = zeros((n_free, n_free),Float)
    beta = zeros((1,n_free),Float)
    delta = (fitparam + equal(fitparam,0.0)) * 0.00001
    nr  = x.shape[0]
    ##############
    # Prior to each call to the function one has to re-calculate the
    # parameters
    pwork = parameters.__copy__()
    for i in range(n_free):
        pwork [free_index[i]] = fitparam [i]
    newpar = getparameters(pwork.tolist(),constrains)
    newpar = take(newpar,noigno)
    for i in range(n_free):
        if model_deriv is None:
            #pwork = parameters.__copy__()
            pwork [free_index[i]] = fitparam [i] + delta [i]
            newpar = getparameters(pwork.tolist(),constrains)
            newpar=take(newpar,noigno)
            f1 = model(newpar, x)
            pwork [free_index[i]] = fitparam [i] - delta [i]
            newpar = getparameters(pwork.tolist(),constrains)
            newpar=take(newpar,noigno)
            f2 = model(newpar, x)
            help0 = (f1-f2) / (2.0 * delta [i])
            help0 = help0 * derivfactor[i]
            pwork [free_index[i]] = fitparam [i]
            #removed I resize outside the loop:
            #help0 = resize(help0,(1,nr))
        else:
            newpar = getparameters(pwork.tolist(),constrains)
            help0=model_deriv(pwork,free_index[i],x)
            help0 = help0 * derivfactor[i]

        if i == 0 :
            deriv = help0
        else:
            deriv = concatenate ((deriv,help0), 0)
    #line added to resize outside the loop
    deriv=resize(deriv,(n_free,nr))
    if linear:
        pseudobetahelp = weight * y
    else:
        yfit = model(newpar, x)
        deltay = y - yfit
        help0 = weight * deltay
    for i in range(n_free):
        derivi = resize(deriv [i,:], (1,nr))
        if linear:
            if i==0:
                beta = resize(sum((pseudobetahelp * derivi),1),(1,1))
            else:
                beta = concatenate((beta, resize(sum((pseudobetahelp * derivi),1),(1,1))), 1)
        else:
            help1 = resize(sum((help0 * derivi),1),(1,1))
            if i == 0:
                beta = help1
            else:
                beta = concatenate ((beta, help1), 1)
        help1 = innerproduct(deriv,weight*derivi)
        if i == 0:
            alpha = help1
        else:
            alpha = concatenate((alpha, help1),1)
    if linear:
        #not used
        chisq = 0.0
    else:
        chisq = sum(help0 * deltay)
    return chisq, alpha, beta, \
           n_free, free_index, noigno, fitparam, derivfactor

def getparameters(parameters,constrains):
    # 0 = Free       1 = Positive     2 = Quoted
    # 3 = Fixed      4 = Factor       5 = Delta
    newparam=[]
    #first I make the free parameters
    #because the quoted ones put troubles
    for i in range(len(constrains [0])):
        if constrains[0][i] == CFREE:
            newparam.append(parameters[i])
        elif constrains[0][i] == CPOSITIVE:
            #newparam.append(parameters[i] * parameters[i])
            newparam.append(abs(parameters[i]))
        elif constrains[0][i] == CQUOTED:
            if 1:
                newparam.append(parameters[i])
            else:
                pmax=max(constrains[1] [i],constrains[2] [i])
                pmin=min(constrains[1] [i],constrains[2] [i])
                A = 0.5 * (pmax + pmin)
                B = 0.5 * (pmax - pmin)
                newparam.append(A + B * sin(parameters[i]))
        elif abs(constrains[0][i]) == CFIXED:
            newparam.append(parameters[i])
        else:
            newparam.append(parameters[i])
    for i in range(len(constrains [0])):
        if constrains[0][i] == CFACTOR:
            newparam[i] = constrains[2][i]*newparam[int(constrains[1][i])]
        elif constrains[0][i] == CDELTA:
            newparam[i] = constrains[2][i]+newparam[int(constrains[1][i])]
        elif constrains[0][i] == CIGNORED:
            newparam[i] = 0
        elif constrains[0][i] == CSUM:
            newparam[i] = constrains[2][i]-newparam[int(constrains[1][i])]
    return newparam

def getsigmaparameters(parameters,sigma0,constrains):
    # 0 = Free       1 = Positive     2 = Quoted
    # 3 = Fixed      4 = Factor       5 = Delta
    n_free = 0
    sigma_par = zeros(parameters.shape,Float)
    for i in range(len(constrains [0])):
        if constrains[0][i] == CFREE:
            sigma_par [i] = sigma0[n_free]
            n_free += 1
        elif constrains[0][i] == CPOSITIVE:
            #sigma_par [i] = 2.0 * sigma0[n_free]
            sigma_par [i] = sigma0[n_free]
            n_free += 1
        elif constrains[0][i] == CQUOTED:
            pmax = max(constrains [1] [i], constrains [2] [i])
            pmin = min(constrains [1] [i], constrains [2] [i])
            A = 0.5 * (pmax + pmin)
            B = 0.5 * (pmax - pmin)
            if (B > 0) & (parameters [i] < pmax) & (parameters [i] > pmin):
                sigma_par [i] = abs(B) * cos(parameters[i]) * sigma0[n_free]
                n_free += 1
            else:
                sigma_par [i] = parameters[i]
        elif abs(constrains[0][i]) == CFIXED:
            sigma_par[i] = parameters[i]
    for i in range(len(constrains [0])):
        if constrains[0][i] == CFACTOR:
            sigma_par [i] = constrains[2][i]*sigma_par[int(constrains[1][i])]
        elif constrains[0][i] == CDELTA:
            sigma_par [i] = sigma_par[int(constrains[1][i])]
        elif constrains[0][i] == CSUM:
            sigma_par [i] = sigma_par[int(constrains[1][i])]
    return sigma_par

def fitpar2par(fitpar,constrains,free_index):
    newparam = []
    for i in range(len(constrains [0])):
        if constrains[0][free_index[i]] == CFREE:
            newparam.append(fitpar[i])
        elif constrains[0][free_index[i]] == CPOSITIVE:
            newparam.append(fitpar[i] * fitpar [i])
        elif abs(constrains[0][free_index[i]]) == CQUOTED:
            pmax=max(constrains[1] [free_index[i]],constrains[2] [free_index[i]])
            pmin=min(constrains[1] [free_index[i]],constrains[2] [free_index[i]])
            A = 0.5 * (pmax + pmin)
            B = 0.5 * (pmax - pmin)
            newparam.append(A + B * sin(fitpar[i]))
    return newparam


class GaussFit:
    def __init__(self,ngauss, R, refinement, Labs):
        self.ngauss=ngauss
        self.exps=zeros(ngauss,"d")
        self.expsI=zeros(ngauss,"d")
        for i in range(ngauss):
            self.expsI[i]=(1+i)* math.pi/(R*R*refinement)
            self.exps[i] =3.0/R/R
            
        self.memcos = None
        self.memsin = None
        self.Labs=Labs
        pass

    def gaussCos(self,param,x, N=None):
        res=0.0
        if self.memcos is None:
            mask=less(x,3.0)
            self.memcos={}
            for i in range((len(self.exps) ) ):
                self.memcos[i] = myexp(-self.exps[i]*x*x)*cos( (self.expsI[i]- self.expsI[0] )*x*x )*(x**(self.Labs+1)) # *mask
        if N is not None:
            return self.memcos[N]
        for i in range(len(self.exps)):
            res=res+(param[i])*self.memcos[i]
        return res

    def gaussSin(self,param,x, N=None):
        res=0.0
        if self.memsin is None:
            mask=less(x,3.0)
            self.memsin={}
            for i in range((len(self.exps) ) ):
                self.memsin[i] = myexp(-self.exps[i]*x*x)*sin( self.expsI[i]*x*x )*(x**(self.Labs+1)) # *mask
        if N is not None:
            return self.memsin[N]
        for i in range(len(self.exps)):
            res=res+(param[i])*self.memsin[i]
        return res
    

    def gauss (self,param,x, N=None):
        # print param[ :self.ngauss ]
        if N is None:
            a=self.gaussCos (param[ :self.ngauss ], x)  
            b=self.gaussSin (param[ self.ngauss: ], x)
            return a+b
        else:
            if N< self.ngauss:
                return self.gaussCos (None, x,  N)
            else:
                return self.gaussSin (None, x , N=N - self.ngauss  )
                

##     def gauss (self,param,x, N=None):

##         if N is None:
##             a=self.gaussCos (param, x)  

##             return a
##         else:
##                 return self.gaussCos (None, x,  N)
                
    

    def deriv(self,param,ideriv, x):
        return  self.gauss( param,x, ideriv  )
    

## +param[i+self.ngauss]*x*x+param[i+self.ngauss*2]*x*x*x*x

def myexp(x):
    # put a (bad) filter to avoid over/underflows
    # with no python looping
    return exp(x*less(abs(x),250))-1.0*greater_equal(abs(x),250)


def fit( ngauss,  refinement, xdata, ydata, Labs ):
    
    if(1):

        xx=[xdata[0] ]
        yy=[ydata[0] ]
        step=0.002
        last=xx[0]
        
        for x,y in zip (xdata[:], ydata[:]) :
            if x>last+step or x==xdata[-1]:
                last=x
                xx.append(x)
                yy.append(y)

        for x in arange(xx[-1]+step,10.0, step):
            xx.append(x)
            yy.append(0.0)
            
        R= xdata[-1]
        npoints=len(xx)
        xx=resize(xx,(npoints,1))
        yy=resize(yy,(npoints,1))
        sy=yy*0 +1
        sy=resize(sy,(npoints,1))
    
    stime = time.time()


    
    gaussfit=GaussFit( ngauss,R, refinement   , Labs   )

    data = concatenate((xx, yy, sy),1)

    if(1):
        parameters = array([1.0]*2*ngauss)
        # parameters = array([1.0]*ngauss)
        gaussfit.mem=None
        fittedpar, chisq, sigmapar = LeastSquaresFit(gaussfit.gauss,parameters,data, deltachi=0.0000000001,model_deriv=gaussfit.deriv  )

##         fakepar=[0.0 for tok in fittedpar]
##         fakepar[0]=1.0

        cy=gaussfit.gauss(fittedpar, xx)

        print cy.shape
        print yy.shape
        print xx.shape
        
##         f=open("p"+"_fit.dat", "w")
##         for i in range(len(xx)):
##             f.write(" %e %e %e %e \n" %(  xx[i,0], yy[i,0], cy[i],  yy[i,0]-cy[i]  )  ) 
##         raise "OK "

##         f=open(nomefunc+"_parameters", "w")
##         for i in range(len( fittedpar)):
##             f.write(" %30.20e %30.20e \n"%(gaussfit.exps[i],  fittedpar[i] ))
            

    f=None
    
    
    etime = time.time()
    print "Took ",etime - stime, "seconds"
    print "chi square  = ",chisq
    print "Fitted pars = ",fittedpar
    print "Sigma pars  = ",sigmapar

    result={}
    result["fitobject"]=gaussfit
    result["xdy"]=[ xx[:,0], yy[:,0], cy[:] ]
    result["fittedpar"]=fittedpar
    return result

            
def AsciiMatrixReader(filename):
   f=open(filename,"r")
   datas=f.read()
   datalines=string.split(datas,"\n")

   Np=len(datalines)

   data=map(string.split,datalines)

   dim=len(data[0])
   for i in range(Np):
     data[i]=map(string.atof,data[i])

   newdata=[]
   for i in range(Np):
     if(len(data[i])==dim):
         newdata.append(data[i])
   data=numpy.array(newdata)
   return data


## if __name__ == "__main__":
##     test(200,4 , x,y)


#   import profile
#   profile.run('test(100,4 , "ptilde_0")',"test")
#   import pstats
#   p=pstats.Stats("test")
#   p.strip_dirs().sort_stats(-1).print_stats()












# pos=string.find(sys.argv[1] , "_in_coeffs")
# if pos==-1:
#     raise " filename expected to end with _in_coeffs" 
# if pos+len("_in_coeffs" )!=len(sys.argv[1]):
#     raise " filename expected to end with _in_coeffs" 



lstart = string.atoi(sys.argv[1])
lend   = string.atoi(sys.argv[2])
ngauss =  string.atoi(sys.argv[3])
ltot=0
l_ltot=[]

for lpartial in range(lstart, lend+1):
    nomefile="pawdata.L=%0d.dat" % lpartial
    s=open( nomefile ,"r").read()
    sl=string.split(s,"\n")
    lpar=0
    for l in sl:
        if "now writing the dual functions and PAWpatch"  in l:
            ltot=ltot+1
            lpar=lpar+1
    l_ltot.append(lpar)

lcount=-1
for lpartial in range(lstart, lend+1):
  lcount+=1
  nomefile="pawdata.L=%0d.dat" % lpartial

  s=open( nomefile ,"r").read()
  sl=string.split(s,"\n")
  for volta in range(l_ltot[lcount]):
    while "now writing the dual functions and PAWpatch" not in sl[0]:
        del sl[0]
    del sl[0]

    Npaw, Ngrid , Labs= map(string.atoi, string.split(sl[0]))


    if Labs<0:
        Labs=-Labs-1
        is_initial=True
    else:
        is_initial=False

    del sl[0]

    rgrid=[]

    for i in range(Ngrid):
        rgrid.append(string.atof(sl[i]))

    ptildes = numpy.zeros([Npaw,Ngrid],numpy.float64)

    for i in range(Npaw):
        sl = sl[Ngrid:]
        del sl[0]
        for igrid in range(Ngrid):
            ptildes[i, igrid]  = string.atof(sl[igrid])

    sl = sl[Ngrid:]
    del sl[0]

    PAWpatch=zeros( [Npaw, Npaw] , numpy.float64 )

    for i in range(Npaw):
        for j in range(Npaw):
            PAWpatch[i,j]=string.atof(sl[0])
            del sl[0]


    Sm1 = numpy.zeros([Npaw,Npaw],numpy.float64)

    for i in range(Npaw):
        for j in range(Npaw):
            y = ptildes[i]*ptildes[j]
            Sm1[i,j]=numpy.trapz(y, x=rgrid)

    S=numpy.linalg.inv(Sm1)



    ## print [ S[i,i] for i in range(len(Sm1)) ]
    ## print S
    ## raise "OK "


    SSm1=numpy.dot( S, S) - S


    SS=-S
    for i in range(len(SS)):
        SS[i,i]=1.0+SS[i,i]


    # ngauss=100
    refinement=2
    results=[]
    for i in range(Npaw):
        result = fit( ngauss,  refinement, rgrid , ptildes[i], Labs )
        results.append(result)

        if(not is_initial) :
            fres=open("fitresult_l_"+str(lpartial)+"_channel_"+str(i),"w")
        else:
            fres=open("fitresult_l_"+str(lpartial)+"_initial","w")

        xx,yy,cc = result["xdy"]
        for x,y,c in zip(xx,yy,cc):
            fres.write("%20.15e %20.15e %20.15e  %20.15e\n"%(x,y,c, c-y))
        fres=None



    gaussfit = results[0]["fitobject"]


    if lpartial==lstart and volta==0:
        f=open("pseudopaw","w")
        f.write("PAWPATCH\n")
        f.write("%d           !!  paw_NofL\n "%  ltot    )



    if is_initial:
        f.write("%d           !! paw_l \n"%  (-Labs-1)   )
    else:
        f.write("%d           !! paw_l \n"%  Labs   )

    f.write("%20.15e           !!  paw_Greal ==> exp( - paw_Greal*r**2) * cos, sin ( paw_Gimag*r**2) \n "%   gaussfit.exps[0]  )
    f.write("%d           !!  paw_nofgaussians \n "%   (gaussfit.ngauss +1) )
    f.write("%d           !! paw_nofchannels  \n"%  Npaw   )
    f.write("!!  follow %d PAW_Gimag factors\n"% (gaussfit.ngauss +1) )

    f.write("%20.15e     \n"%   0.0 )
    for i in range( gaussfit.ngauss):
        f.write("%20.15e     \n"%   gaussfit.expsI[i]  )

    f.write("!!  follow for each of the %d channels %d (cos_factor, sin_factor)  pairs\n"% (   Npaw,  gaussfit.ngauss)  )
    for j in range( Npaw):
        fittedpar = results[j]["fittedpar"]

        f.write("%20.15e %20.15e\n"%   ( fittedpar [ 0 ], 0.0 )  )

        for i in range( gaussfit.ngauss-1):
            f.write("%20.15e %20.15e\n"%   ( fittedpar [ i +1], fittedpar [  gaussfit.ngauss + i  ])  )

        f.write("%20.15e %20.15e\n"%   ( 0.0 , fittedpar [  gaussfit.ngauss +  gaussfit.ngauss-1 ])  )

    f.write("!! pawpatch matrix \n")
    for i in range( Npaw):
        for j in range( Npaw):
            f.write("%20.15e\n"% PAWpatch[i,j]) 


    f.write("!! S matrix \n")
    for i in range( Npaw):
        for j in range( Npaw):
            f.write("%20.15e\n"%SS[i,j]) 
    f.write("!! Sm1 matrix \n")
    for i in range( Npaw):
        for j in range( Npaw):
            f.write("%20.15e\n"%SSm1[i,j]) 
