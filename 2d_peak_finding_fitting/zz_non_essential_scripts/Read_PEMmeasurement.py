# -*- coding: utf-8 -*-
"""
Created on Thu Mar 28 16:44:49 2019

@author: rober
"""

import os, numpy as np, csv, matplotlib.pyplot as plt, scipy.optimize as opt, math, struct, binascii, gc, time, random
import multiprocessing
from operator import sub
from joblib import Parallel, delayed
#import scipy, lmfit
from scipy.optimize import minimize # used for implementation of maximum likelihood exponential fit
from matplotlib import gridspec
from math import factorial
from math import *
from scipy.stats import poisson
get_ipython().run_line_magic('matplotlib', 'auto')
import matplotlib as mpl
import pickle
import numba as nb
from joblib import Parallel, delayed
import multiprocessing
from tqdm import tqdm

#%matplotlib auto

# In[10]:


#I get errors unless the function that the cores execute in parallel is defined outside the class function 
def hist2(x,y,bins):
    store = np.zeros(len(bins)-1,dtype='float');
    for i in x:
        res = y.searchsorted(bins+i)
        store += res[1:]-res[:-1]
    return store
def load_obj(name, folder ):
    with open(folder + name + '.pkl', 'rb') as f:
        return pickle.load(f)

# In[23]:


def ImportT3(filename):
    with open(filename, "rb+") as f:
    
#        while True:
#            if f.read(1) == b'M':
#                if f.read(18) == b'easDesc_Resolution': # recognize time unit entry
#                    break
#        f.read(21) # rest of the tag
#        dtmicro = struct.unpack('d',f.read(8))[0]
#        #print('Microtime unit:', dtmicro)
#    
#        while True:
#            if f.read(1) == b'M':
#                if f.read(24) == b'easDesc_GlobalResolution': # recognize time unit entry
#                    break
#        f.read(15) # rest of the tag
#        dtmacro = struct.unpack('d',f.read(8))[0]
#        #print('Macrotime unit:', dtmacro)
#    
#        while True:
#            if f.read(1) == b'T':
#                if f.read(23) == b'TResult_NumberOfRecords': # recognize number of records entry
#                    break
#        f.read(16) # rest of the tag
#        nrrec = struct.unpack('q',f.read(8))[0] # extract number of records
#        #print('Number of records in file:', nrrec)
#
#        while True:
#            if f.read(1) == b'H':
#                if f.read(9) == b'eader_End':
#                    #print('Header_End found')
#                    break
#        f.read(38) # rest of Header_End
        nrrec=HHsettings["overallCounts"]
#        nrrec=int(1e5)
        dtmicro=HHsettings["resolution"]*1e-12    #in s
        dtmacro=1/HHsettings["syncRate"]    #in s
        macrotimes0 = np.zeros(nrrec,dtype='int64');
        microtimes0 = np.zeros(nrrec,dtype='int64');
        macrotimes1 = np.zeros(nrrec,dtype='int64');
        microtimes1 = np.zeros(nrrec,dtype='int64');
        macrotimesfireA = np.zeros(nrrec,dtype='int64');
        microtimesfireA = np.zeros(nrrec,dtype='int64');
        macrotimesfireB = np.zeros(nrrec,dtype='int64');
        microtimesfireB = np.zeros(nrrec,dtype='int64');
        overflows = 0
        nrphotons0 = 0
        nrphotons1 = 0
        nrfireA = 0
        nrfireB = 0
        prevchann = 0
        
        for i in range(nrrec):
            entry = f.read(4)
            channel = struct.unpack("I",entry)[0] >> 25 # read channel number, first 7 bits
            if channel == 0:
                macrotime = (struct.unpack("I",entry)[0] & 0x3FF)
                macrotimes0[nrphotons0] = macrotime + 1024*overflows
                microtime = ((struct.unpack("I",entry)[0] >> 10) & 0x7FFF)
                microtimes0[nrphotons0] = microtime
                nrphotons0 += 1
                prevchann = 0
            elif channel == 1:
                macrotime = (struct.unpack("I",entry)[0] & 0x3FF)
                macrotimes1[nrphotons1] = macrotime + 1024*overflows
                microtime = ((struct.unpack("I",entry)[0] >> 10) & 0x7FFF)
                microtimes1[nrphotons1] = microtime
                nrphotons1 += 1
                prevchann = 1
            elif channel == 127:
                nroverflows = (struct.unpack("I",entry)[0] & 0x3FF)
                overflows += nroverflows
                prevchann = 127
            elif channel == 65:
                macrotime = (struct.unpack("I",entry)[0] & 0x3FF)
                macrotimesfireA[nrfireA] = macrotime + 1024*overflows
                microtime = ((struct.unpack("I",entry)[0] >> 10) & 0x7FFF)
                microtimesfireA[nrfireA] = microtime
                nrfireA += 1
            elif channel == 66:
                macrotime = (struct.unpack("I",entry)[0] & 0x3FF)
                macrotimesfireB[nrfireB] = macrotime + 1024*overflows
                microtime = ((struct.unpack("I",entry)[0] >> 10) & 0x7FFF)
                microtimesfireB[nrfireB] = microtime
                nrfireB += 1                                
            else:
                print('bad channel:',channel)
                
    microtimes0 = microtimes0[:nrphotons0]
    macrotimes0 = macrotimes0[:nrphotons0]
    microtimes1 = microtimes1[:nrphotons1]
    macrotimes1 = macrotimes1[:nrphotons1]
    microtimesfireA = microtimesfireA[:nrfireA]
    macrotimesfireA = macrotimesfireA[:nrfireA]
    microtimesfireB = microtimesfireB[:nrfireB]
    macrotimesfireB = macrotimesfireB[:nrfireB]
    
    print('nrphotons0:',nrphotons0)
    print('nrphotons1:',nrphotons1)
    print('nrfireA:',nrfireA)
    print('nrfireB:',nrfireB)
    print('overflows:',overflows)
    
    return [dtmicro, dtmacro, microtimes0, macrotimes0, microtimes1, macrotimes1, nrphotons0,nrphotons1,overflows,microtimesfireA,macrotimesfireA,nrfireA,microtimesfireB,macrotimesfireB,nrfireB]

def ShiftPulsedData(microtimes0,microtimes1,macrotimes0,macrotimes1,dtmicro,dtmacro):
    dtmax = 8
    
    [ylist1,xlist1] = np.histogram(microtimes1,int(dtmacro/dtmicro),[0,int(dtmacro/dtmicro)])
    [ylist0,xlist0] = np.histogram(microtimes0,int(dtmacro/dtmicro),[0,int(dtmacro/dtmicro)])
    tlist = (xlist0[:-1]+0.5*(xlist0[1]-xlist0[0]))*dtmicro*1e9

    corrx = []; corry = [] # find shift for which the two decay curves overlap most
    for i in range(-dtmax,dtmax):
        corrx.append(i)
        corry.append(sum(ylist1[dtmax:-dtmax]*ylist0[dtmax+i:-dtmax+i]))
    xmax = corry.index(max(corry))
    shift = corrx[xmax]
    
    tlist0 = (microtimes0-shift) + macrotimes0*int(dtmacro/dtmicro)
    tlist1 = microtimes1 + macrotimes1*int(dtmacro/dtmicro) #in units of dtmicro
    
#    plt.xlabel('time (ns)')
#    plt.ylabel('counts (a.u.)')
#    p1, = plt.plot(tlist,ylist0+ylist1)
#    p2, = plt.plot(tlist,ylist0)
#    p3, = plt.plot(tlist,ylist1)
#    plt.legend([p1,p2,p3], ["APD0 + APD1","APD0","APD1"])
           
    return(microtimes0-shift,microtimes1,tlist0,tlist1,dtmicro,dtmacro,tlist,ylist0+ylist1)

def GetLifetime(microtimes,dtmicro,dtmacro,dtfit,tstart,binwidth=1,ybg=0): #dtfit is the time interval considered for t
    # binwidth is a multiplier. actual binwidth is given as binwidth*dtmicro[s]
    [ylist,xlist] = np.histogram(microtimes,int(dtmacro/(dtmicro*binwidth)),[0,int(dtmacro/dtmicro)])
    tlist = (xlist[:-1]+0.5*(xlist[1]-xlist[0]))*dtmicro*1e9
    
    istart = int(tstart/dtmicro) #find index of maximum element in ylist
    if istart < 0:
        istart = ylist.argmax()
    iend = istart + int(dtfit/(dtmicro*binwidth))
    if iend>len(tlist):
        iend = len(tlist) 
        
    # get background (by simply looking at last ten data points) and substract from intensity data.
    if ybg == 0:
        ybg = np.mean(ylist[-10:-1]) # mean background per histogram bin bin of length
            
    # exponential decay    
    #expfit = np.polyfit(tlist[istart:iend], np.log(ylist[istart:iend]), 1, w=np.sqrt(ylist[istart:iend])) #weighted exponential fit
    #tauexp = -1/expfit[0]
    
    # new exponential fit
    #[expfit,_] = scipy.optimize.curve_fit(lambda t,a,b: a*np.exp(b*t)+ybg,  tlist[istart:iend],  ylist[istart:iend],  p0=(1e3, -1/40))
    #tauexp = -1/expfit[1]
    #Aexp = expfit[0]
    
    #maximum likelihood exponential fit    
    [tauexp,Aexp] = MaxLikelihoodFit(tlist,ylist,istart,iend,ybg,False)
    
    # average lifetime
    tauave = np.sum((tlist[istart:iend]-tlist[istart])*(ylist[istart:iend]-ybg)/np.sum(ylist[istart:iend]-ybg))
      
    plt.xlabel('time (ns)')
    plt.ylabel('')
    #plt.semilogy(tlist,ylist,'.',tlist[istart:iend],np.exp(expfit[1])*np.exp(expfit[0]*tlist[istart:iend]))
    plt.plot(tlist,ylist/np.mean(ylist),'.',tlist[istart:iend],Aexp*np.exp(-(tlist[istart:iend]-tlist[istart])/tauexp)+ybg)
    #plt.semilogy([tlist[0],tlist[-1]],[ybg,ybg],'k--')
    #plt.show()

    # Amax is the maximum y-value
    Amax = np.max(ylist)
        
    print('single-exp lifetime:',tauexp,'ns; average lifetime:',tauave,'ns; Amax:',Amax)
    
    return(tauave,Amax,ybg) 

def HistPhotons(photontlist,binwidth,Texp): # finds the index of the first photon in each bin. photontlist is in [s]
    histmax = Texp # experiment duration [s]

    nrbins = int(Texp/binwidth)
    limits = np.full(nrbins,len(photontlist))
    counter,i = 0,0
    while counter < nrbins and i < len(photontlist):
        while photontlist[i] > counter*binwidth:
            limits[counter] = i
            counter += 1
        i += 1
    
    return(limits)

def MakeIntTrace(limits0,limits1,binwidth,Texp):
    nrbins = int(Texp/binwidth)
    inttrace = np.array([(limits0[binnr+1]-limits0[binnr]+limits1[binnr+1]-limits1[binnr]) for binnr in range(nrbins-1)])

    fig = plt.figure(figsize=(15,3))
    gs = gridspec.GridSpec(1,2,width_ratios=[4,1])
    ax0 = plt.subplot(gs[0])
    p0 = ax0.plot(np.arange(len(inttrace))*binwidth,inttrace,'-',linewidth=0.5)
    plt.xlabel('time (s)')
    plt.ylabel('counts / %i ms' %(binwidth*1e3))
    plt.xlim([0,Texp])
    plt.ylim([0,1.1*np.max(inttrace)])

    histogram = np.histogram(inttrace,max(inttrace),[0,max(inttrace)])
    
    ax1 = plt.subplot(gs[1])
    ax1.plot(histogram[0],0.5*(histogram[1][:-1]+histogram[1][1:]))
    plt.xlabel('occurrence')
    plt.ylabel('counts / %i ms' %(binwidth*1e3))
    plt.ylim([0,1.1*np.max(inttrace)])
    
    return(inttrace)

def MakeTauTrace(taubinlist,intbinlist,binwidth,Texp,taumin=0,taumax=100,intmin=0,intmax=100,col='k'):
    nrbins = int(Texp/binwidth)
  
    fig = plt.figure(figsize=(15,7))
    gs = gridspec.GridSpec(2,2,width_ratios=[4,1])
    ax0 = plt.subplot(gs[0])
    p0 = ax0.plot(np.arange(len(intbinlist))*binwidth,intbinlist,'-',linewidth=0.5,color=col)
    plt.xlabel('time (s)')
    plt.ylabel('counts / %i ms' %(binwidth*1e3))
    plt.xlim([0,Texp])
    plt.ylim([intmin,intmax])

    histogram = np.histogram(intbinlist,int(np.max(intbinlist)),[0,int(np.max(intbinlist))])
    
    ax1 = plt.subplot(gs[1])
    ax1.plot(histogram[0],0.5*(histogram[1][:-1]+histogram[1][1:]),color=col)
    plt.xlabel('occurrence')
    plt.ylabel('counts / %i ms' %(binwidth*1e3))
    plt.ylim([intmin,intmax])
    
    ax2 = plt.subplot(gs[2])
    p2 = ax2.plot(np.arange(len(intbinlist))*binwidth,taubinlist,'.',markersize=1.5,color=col)
    plt.xlabel('time (s)')
    plt.ylabel('lifetime (ns)')
    plt.xlim([0,Texp])
    plt.ylim([taumin,taumax])

    histogram = np.histogram(taubinlist,taumax-taumin,[taumin,taumax])
    
    ax3 = plt.subplot(gs[3])
    ax3.plot(histogram[0],0.5*(histogram[1][:-1]+histogram[1][1:]),color=col)
    plt.xlabel('occurrence')
    plt.ylabel('lifetime (ns)')
    plt.ylim([taumin,taumax])
    
    plt.hold(True)
    

def BinIntensity(microtimes0,times0,limits0,microtimes1,times1,limits1,dtmicro,dtmacro,onintlim,offintlim):
    ## select only data with high or low intensity

    plt.title('total decay')
    tauave = GetLifetime(np.append(microtimes0,microtimes1),dtmicro,dtmacro,200e-9,-1)

    nrbins = len(limits0)
    inttrace = np.array([limits0[binnr+1]-limits0[binnr]+limits1[binnr+1]-limits1[binnr] for binnr in range(nrbins-1)])

    # find photons in on period
    onphotonlist0 = np.array([np.arange(limits0[binnr],limits0[binnr+1]) for binnr in range(nrbins-1) if inttrace[binnr] >= onintlim])
    onphotonlist0 = np.concatenate(onphotonlist0).ravel()
    onphotonlist1 = np.array([np.arange(limits1[binnr],limits1[binnr+1]) for binnr in range(nrbins-1) if inttrace[binnr] >= onintlim])
    onphotonlist1 = np.concatenate(onphotonlist1).ravel()

    onmicrotimes0 = np.array([microtimes0[i] for i in onphotonlist0])
    onmicrotimes1 = np.array([microtimes1[i] for i in onphotonlist1])
    ontimes0 = np.array([times0[i] for i in onphotonlist0])
    ontimes1 = np.array([times1[i] for i in onphotonlist1])
    plt.title('on decay')
    ontauave = GetLifetime(np.append(onmicrotimes0,onmicrotimes1),dtmicro,dtmacro,200e-9,-1)

    # find photons in off period
    offphotonlist0 = np.array([np.arange(limits0[binnr],limits0[binnr+1]) for binnr in range(nrbins-1) if inttrace[binnr] < offintlim])
    offphotonlist0 = np.concatenate(offphotonlist0).ravel()
    offphotonlist1 = np.array([np.arange(limits1[binnr],limits1[binnr+1]) for binnr in range(nrbins-1) if inttrace[binnr] < offintlim])
    offphotonlist1 = np.concatenate(offphotonlist1).ravel()

    offmicrotimes0 = np.array([microtimes0[i] for i in offphotonlist0])
    offmicrotimes1 = np.array([microtimes1[i] for i in offphotonlist1])
    offtimes0 = np.array([times0[i] for i in offphotonlist0])
    offtimes1 = np.array([times1[i] for i in offphotonlist1])
    plt.title('off decay')
    offtauave = GetLifetime(np.append(offmicrotimes0,offmicrotimes1),dtmicro,dtmacro,10e-9,-1)
    
    return(onmicrotimes0,offmicrotimes0,ontimes0,offtimes0,onmicrotimes1,offmicrotimes1,ontimes1,offtimes1)

def SliceHistogram(microtimes0,times0,limits0,microtimes1,times1,limits1,dtmicro,dtmacro,Imin,Imax):
    ## select only data with intensity between Imin and Imax

    nrbins = len(limits0)
    inttrace = np.array([limits0[binnr+1]-limits0[binnr]+limits1[binnr+1]-limits1[binnr] for binnr in range(nrbins-1)])

    # find photons in bins with intensities in (Imin,Imax] range
    onphotonlist0 = np.array([np.arange(limits0[binnr],limits0[binnr+1]) for binnr in range(nrbins-1) if Imin < inttrace[binnr] <= Imax])
    onphotonlist0 = np.concatenate(onphotonlist0).ravel()
    onphotonlist1 = np.array([np.arange(limits1[binnr],limits1[binnr+1]) for binnr in range(nrbins-1) if Imin < inttrace[binnr] <= Imax])
    onphotonlist1 = np.concatenate(onphotonlist1).ravel()
      
    onmicrotimes0 = np.array([microtimes0[i] for i in onphotonlist0])
    onmicrotimes1 = np.array([microtimes1[i] for i in onphotonlist1])
    ontimes0 = np.array([times0[i] for i in onphotonlist0])
    ontimes1 = np.array([times1[i] for i in onphotonlist1])
    
    # count nr of time bins with intensities corresponding to slice intensity
    onbincount = 0
    for binnr in range(nrbins-1):
        if Imin < inttrace[binnr] <= Imax:
            onbincount +=1

    return(onmicrotimes0,ontimes0,onmicrotimes1,ontimes1,onbincount)



def MaxLikelihoodFit(tlist,ylist,istart,iend,bgcpb,plotbool=False):
    ### Maximum likelihood routine to fit single exponential. Pro: Works also for small amount of data (single bins of 10ms!)
    # tlist: x-axis values, here time in ns; ylist: y-axis values, here cts per tlist-bin; istart and iend: first and last element of tlist and ylist that are considered for the fit.

    # check if istart and iend are good numbers
    if istart<0 or istart>=len(ylist):
        istart = 0
        print('WARNING: adapted istart in MaxLikelihoodExpFit')
    if iend<=istart or iend>len(ylist):
        iend = len(ylist)
        print('WARNING: adapted iend in MaxLikelihoodExpFit')

    # shift t0 to t=0
    Seterr=1/1.16
    duration=dtmacro*1e9
    ydata = ylist[istart:iend]
    xdata = tlist[istart:iend]

    # do calculations
    initParams = [np.max(ydata), 0.8, 0] #initial guess for A and tau
    results = minimize(MaxLikelihoodFunction, initParams, args=(xdata,ydata,bgcpb),method='Nelder-Mead') # minimize the negative of the maxlikelihood function instead of maximimizing
    Aest = results.x[0] # get results of fit, A
#    Aest = 280
#    Alphaest = 0.5
#    Thetaest = 0
    Alphaest = results.x[1] # get results of fit, tau
    Thetaest = results.x[2]

    if plotbool == True:
#        yest = np.array([Aest*np.exp(-(xdata[i]-xdata[0])/tauest)+bgcpb for i in range(len(xdata))])
        yest = np.array([Aest*((1+Alphaest)*(1+Seterr**2+(-1+Seterr**2)*np.cos(np.pi*np.cos(2*np.pi*(xdata[i]-xdata[0])/duration))) + (-1+Alphaest)*np.cos(2*Thetaest)*(-1+Seterr**2+(1+Seterr**2)*np.cos(np.pi*np.cos(2*np.pi*(xdata[i]-xdata[0])/duration)))+2*Seterr*(-1+Alphaest)*np.sin(2*Thetaest)*np.sin(np.pi*np.cos(2*np.pi*(xdata[i]-xdata[0])/duration)))/4 for i in range(len(xdata))])
        plt.plot(tlist,ylist,'.',xdata,yest,[xdata[1],xdata[-1]],[bgcpb,bgcpb],'k--')
#        plt.plot(tlist,ylist)
        plt.show()        
    
    return(Aest,Alphaest,Thetaest)
    
def MaxLikelihoodFit_c(tlist,ylist,istart,iend,bgcpb,initParams,plotbool=False):
    ### Maximum likelihood routine to fit single exponential. Pro: Works also for small amount of data (single bins of 10ms!)
    # tlist: x-axis values, here time in ns; ylist: y-axis values, here cts per tlist-bin; istart and iend: first and last element of tlist and ylist that are considered for the fit.

    # check if istart and iend are good numbers
    if istart<0 or istart>=len(ylist):
        istart = 0
        print('WARNING: adapted istart in MaxLikelihoodExpFit')
    if iend<=istart or iend>len(ylist):
        iend = len(ylist)
        print('WARNING: adapted iend in MaxLikelihoodExpFit')

    # shift t0 to t=0
    Seterr=1/1.03
    duration=dtmacro*1e9
    ydata = ylist[istart:iend]
    xdata = tlist[istart:iend]
    initParams[0]=np.max(ylist)
    # do calculations
    results = minimize(MaxLikelihoodFunction_c, initParams, args=(xdata,ydata,bgcpb,duration,Seterr),method='Nelder-Mead',options={'disp':False,'adaptive':True}) # minimize the negative of the maxlikelihood function instead of maximimizing
    
    Aest = results.x[0] # get results of fit, A
#    Aest = 280
#    Alphaest = 0.5
#    Thetaest = 0
    Alphaest = results.x[1] # get results of fit, tau
    Thetaest = results.x[2]    
    
    if plotbool == True:
        yest = np.array([QWPFitfunction(Aest,Alphaest,Seterr,duration,Thetaest,xdata[i]-xdata[0],bgcpb) for i in range(len(xdata))])
        plt.figure()
        plt.plot(tlist,ylist,'.',xdata,yest,[xdata[1],xdata[-1]],[bgcpb,bgcpb],'k--')
#        plt.plot(xdata,yest,[xdata[1],xdata[-1]],[bgcpb,bgcpb],'k--')
        plt.xlim([xdata[1],xdata[-1]])
        plt.show()        
        
    return(Aest,Alphaest,Thetaest)
    
def QWPFitfunction(A,Alpha,Seterr,duration,Theta,t,const):
    QWPfunc=(A*((1+Alpha)*(1+Seterr**2+(-1+Seterr**2)*np.cos(np.pi*np.cos(2*np.pi*(t)/duration))) + (-1+Alpha)*np.cos(2*Theta)*(-1+Seterr**2+(1+Seterr**2)*np.cos(np.pi*np.cos(2*np.pi*(t)/duration)))+2*Seterr*(-1+Alpha)*np.sin(2*Theta)*np.sin(np.pi*np.cos(2*np.pi*(t)/duration)))/4+const)
    return QWPfunc
    
def MaxLikelihoodFunction(params,xdata,ydata,const,duration,Seterr):
    # max likelihood function for A*exp(-t/tau), needed in function MakLikelihoodFit
    # params = [A,tau]
    A = params[0]
    Alpha = params[1]
    Theta=params[2]
    E = 0.;
    for i in range(len(xdata)):
#        E = E + ydata[i]*np.log(A*np.exp(-(xdata[i]-xdata[0])/tau)+const)-(A*np.exp(-(xdata[i]-xdata[0])/tau)+const)
        E = E +ydata[i]*np.log((A*((1+Alpha)*(1+Seterr**2+(-1+Seterr**2)*np.cos(np.pi*np.cos(2*np.pi*(xdata[i]-xdata[0])/duration))) + (-1+Alpha)*np.cos(2*Theta)*(-1+Seterr**2+(1+Seterr**2)*np.cos(np.pi*np.cos(2*np.pi*(xdata[i]-xdata[0])/duration)))+2*Seterr*(-1+Alpha)*np.sin(2*Theta)*np.sin(np.pi*np.cos(2*np.pi*(xdata[i]-xdata[0])/duration)))/4+const))-(A*((1+Alpha)*(1+Seterr**2+(-1+Seterr**2)*np.cos(np.pi*np.cos(2*np.pi*(xdata[i]-xdata[0])/duration))) + (-1+Alpha)*np.cos(2*Theta)*(-1+Seterr**2+(1+Seterr**2)*np.cos(np.pi*np.cos(2*np.pi*(xdata[i]-xdata[0])/duration)))+2*Seterr*(-1+Alpha)*np.sin(2*Theta)*np.sin(np.pi*np.cos(2*np.pi*(xdata[i]-xdata[0])/duration)))/4+const)
        
    return(-E) # This function needs to be MINIMIZED (because of the minus sign) to have the maximum likelihood fit!

def processInput(tbinnr):
    microtimes = np.append(microtimes0[limits0[tbinnr]:limits0[tbinnr+1]],microtimes1[limits1[tbinnr]:limits1[tbinnr+1]])
    [ylist,xlist] = np.histogram(microtimes,int(dtmacro/(dtmicro*histbinmultiplier)),[0,int(dtmacro/dtmicro)])    
    tlist = (xlist[:-1]+0.5*(xlist[1]-xlist[0]))*dtmicro*1e9 # convert x-axis to time in ns
#    plt.clf()
#    plt11 = plt.figure(11)
    result=MaxLikelihoodFit_c(tlist,ylist,istart,iend,bgcpb,initParams,plotbool=False)
    Ampl[tbinnr]=result[0]
    Alpha[tbinnr]=result[1]
    Theta[tbinnr]=result[2]
    return(result)

def MakeG2(times0,times1,dtmicro,g2restime=8e-9,nrbins=200):
    i0=0
    i1=0
    lim1=0
    g2 = np.zeros(2*nrbins)
    #g2B = np.zeros(2*nrbins)
    #g2C = np.zeros(2*nrbins)
    #blindB = 2e-9
    #blindC = 5e-9

    g2res = g2restime/dtmicro #transform g2restime [s] to g2res [microtime units]
    #blindB = blindB/tmicro
    #blindC = blindC/tmicro

    # correlate det0 with det1 (positive time differences)
    for i0 in range(len(times0)):
        t0 = times0[i0]
        i1 = 0
        q = 0
        while q == 0: 
            if lim1 + i1 < len(times1): # check if we've already reached end of photon stream on det1
                dt = times1[lim1+i1]-t0 
                if dt < 0: # find index lim1 of first photon on det1 that came after photon i0 on det0
                    lim1 = lim1 + 1
                else:
                    binnr = int(dt/g2res) # calculate binnr that corresponds to dt
                    if binnr < nrbins: # check if time difference is already large enough to stop correlation
                        g2[nrbins + binnr] += 1 # increase counter in corresponding bin by one
                        #if microtimes0[i0] > blindB and microtimes1[lim1+i1] > blindB:poi
                            #g2B[nrbins + binnr] += 1
                        #if microtimes0[i0] > blindC and microtimes1[lim1+i1] > blindC:
                            #g2C[nrbins + binnr] += 1
                        i1 = i1 + 1 # look at next photon on det1
                    else:
                        q = 1 # dt larger than maximum correlation width. stop. 
            else:
                q = 1 # end of photon stream on det1 reached. stop.

    # correlate det1 with det0 (positive time differences)
    lim1=0
    for i0 in range(len(times1)):
        t0 = times1[i0]
        i1 = 0
        q = 0
        while q == 0:
            if lim1 + i1 < len(times0):
                dt = times0[lim1+i1]-t0
                if dt < 0:
                    lim1 = lim1 + 1
                else:
                    binnr = int(dt/g2res)
                    if binnr < nrbins:
                        g2[nrbins - 1 - binnr] += 1
                        #if microtimes0[lim1+i1] > blindB and microtimes1[i0] > blindB:
                        #    g2B[nrbins - 1 - binnr] += 1
                        #if microtimes0[lim1+i1] > blindC and microtimes1[i0] > blindC:
                        #    g2C[nrbins - 1 - binnr] += 1
                        i1 = i1 + 1
                    else:
                        q = 1
            else:
                q = 1
                                
    g2tlist = np.arange(-g2res*dtmicro*(nrbins-0.5),g2res*dtmicro*nrbins,g2restime)/dtmacro
    plt.plot(g2tlist,g2,'-.')
    # plt.plot(g2)
    #plt.plot(g2tlist,g2B)
    #plt.plot(g2tlist,g2C)
    plt.title('g(2) correlation')
    plt.xlabel('delay (pulses)')
    plt.ylabel('occurence (a.u.)')
    plt.ylim([0,max(g2)])
    plt.show()

    return(g2tlist,g2,g2restime,nrbins)
# In[25]:

### GENERAL EVALUATION OF TTTR data     

# parameters (use forward slashes!)
folder = 'C:/Users/rober/Documents/Doktorat/Projects/Linear_absorption_anisotropy/Andor Spectrometer/20190521_linearanisotropy/scans/measurements/'
folder = 'C:/Users/rober/Documents/Doktorat/GitLab/apd_scanning_imaging/demo_python_preliminary/64bit/output/'
folder = 'E:/LAB_DATA/Robert/20200527_cryo_HR63HR165/'
#folder = 'E:/LAB_DATA/Yannik_G/20190415_rods/Scans/measurements_HWP/'
namelist = ['HR63_film_RT_codetest_singledot_0mV_OD1']    #Paste the file that ends in _tttrmode.out in here without the .out
# settingsfile=namelist[0].rsplit('_',1)[0]+'_HHsettings'     #HHsettings file should have the same file name apart from this last identifier
settingsfile= namelist[0]+'_settings'
HHsettings=load_obj(settingsfile, folder )
#MCLdata=load_obj("Sample1_x174_y187_z11_range20_pts100_MCLdata", folder )

Texp = 500 # total time [s] of measurements
binwidth = 0.01 # width of time bins to create intensity trace [s]
nrbins = 5 # intensity bin to separate intensity trace in (nrbins) to separately evaluate different intensity levels
histbinmultiplier = 1 #(histbinmultiplier)*dtmicro = bin size for lifetime histogram
dttau = 50e-9 # end time of lifetime fit [s]. Note: this is not the width but the end time! NEED TO CHANGE LATER!

savebool = False;
# subtractbg = True; # not used?! 
singlebintaubool = False;
g2bool = False;

colmap = np.array([[0.0,0.0,0.0],
                   [0.6,0.6,1.0],
                   [0.3,0.3,1.0],
                   [0.0,0.0,1.0],
                   [1.0,0.6,0.6],
                   [1.0,0.3,0.3],
                   [1.0,0.0,0.0],
                  ])

# init some arrays and values
nrmeas = len(namelist) #nr of data files
nrtbins = int(Texp/binwidth) # nr of bins in intensity trace
Ilimslist = np.ndarray(shape=(nrmeas,nrbins),dtype=float)
Alphabinlist = np.full((nrmeas,nrbins),0,dtype=float) # init list to save Alpha of each intensity bin
Thetabinlist = np.full((nrmeas,nrbins),0,dtype=float) # init list to save Alpha of each intensity bin
Abinlist = np.full((nrmeas,nrbins),0,dtype=float) # init list to save maximum histogram amplitudes (at t=t0)
taubinlist = np.zeros(nrtbins-1,dtype=float) # init list to save lifetimes of each time bin (binwidth)
intbinlist = np.zeros(nrtbins-1,dtype=float) # init list to save intensity per bin

for measnr in range(nrmeas): #for each measurement file

    # load data
    data = ImportT3(folder + namelist[measnr] + '.out') # import data from .out file; structure of data can be seen in function definition of ImportT3
    Texp = round(data[8]*data[1]*1024,0) # macrotime values can store up to 2^10 = 1024 before overflow. overflows*1024*dtmacro gives experiment time [s]
    print('averaged cps on det0 and det1:',np.array(data[6:8])/Texp)
    print('experimental time in s:',Texp)

    # compensate detector offset; here we also make the microtimetrace
    [microtimes0,microtimes1,times0,times1,dtmicro,dtmacro,decaytlist,decayylist] = ShiftPulsedData(data[2],data[4],data[3],data[5],data[0],data[1]) #decaytlist and decayylist are the two variables you want to check for the modulation trace
    histbinmultiplier=4
    [ylist,xlist] = np.histogram(microtimes0,int(dtmacro/(dtmicro*histbinmultiplier)),[0,int(dtmacro/dtmicro)])    
    tlist = (xlist[:-1]+0.5*(xlist[1]-xlist[0]))*dtmicro*1e9 # convert x-axis to time in ns
    
    nrnewbins=int(dtmacro/(dtmicro*histbinmultiplier))
    tlist = (xlist[:-1]+0.5*(xlist[1]-xlist[0]))*dtmicro*1e9 # convert x-axis to time in ns
    istart=0
    iend=max(tlist)
    bgcounts=Texp*80
    bgcpb=bgcounts/nrnewbins
#    plt.plot(tlist,ylist)
    #    plt.clf()
#    MaxLikelihoodFunction_c = nb.jit(nopython=True)(MaxLikelihoodFunction) # implement c function
#    initParams = [np.max(ylist), 0.5, 1] #initial guess for A, Alpha and Theta; A will get overridden with the maximum value of ylist
#    result=MaxLikelihoodFit_c(tlist,ylist,istart,iend,bgcpb,initParams,plotbool=True)
#    #binning[i]=histbinmultiplier
#    Alpha_av=result[1]
#    Theta_av=result[2]
#    if result[1]>1:
#        Alpha_av=1/result[1]
#        Theta_av=result[2]+np.pi/4
#        print("Corrected fit")
#    print(Alpha_av,Theta_av)
#    
    # Potentially it makes sense to bin this for the fitting
#    plt.show()
    if savebool == True:        
        np.save(folder + namelist[measnr] + '_decay_tlist', decaytlist)
        np.save(folder + namelist[measnr] + '_decay_ylist', decayylist)   
    tmax = decaytlist[decayylist.argmax()] # find histogram time with max photon count [ns]
    
    # get number of first photon per bin for both detectors (APD1 = 0, APD2 = 1)
#    binwidth = 1
    limits0 = HistPhotons(times0*dtmicro,binwidth,Texp) #gives index of first photon of det0 in each bin
    limits1 = HistPhotons(times1*dtmicro,binwidth,Texp)

    # make an intensity trace and find Imax
    inttrace = MakeIntTrace(limits0,limits1,binwidth,Texp)  #this gives the trace of integrated counts vs time where we see the blinking. You might need plt.show() to see it
    Imax = np.max(inttrace) # search overall maximum intensity per intensity bin
    Ilims = np.arange(0,Imax*0.95+1,Imax*0.95/nrbins,dtype=float) # separate intensities into (nrbins) intensity bins
    Ilimslist[measnr,:] = np.array([Ilims[binnr]*0.5+Ilims[binnr+1]*0.5 for binnr in range(nrbins)]) # get centered intensity value per intensity bin

#%% Fit to intensity bins

for binnr in range(nrbins): # for each intensity bin
    # select photons in bins with intensity within Ilims(binnr) and Ilims(binnr+1)
    if binnr>-1:
        histbinmultiplier=16
        [binmicrotimes0,bintimes0,binmicrotimes1,bintimes1,onbincount] = SliceHistogram(microtimes0,times0,limits0,microtimes1,times1,limits1,dtmicro,dtmacro,Ilims[binnr],Ilims[binnr+1])  
        [ylist,xlist] = np.histogram(binmicrotimes0,int(dtmacro/(dtmicro*histbinmultiplier)),[0,int(dtmacro/dtmicro)])    
        tlist = (xlist[:-1]+0.5*(xlist[1]-xlist[0]))*dtmicro*1e9 # convert x-axis to time in ns
        istart=0
        iend=max(tlist)
        nrnewbins=int(dtmacro/(dtmicro*histbinmultiplier))
        bgcounts=binwidth*onbincount*80
        bgcpb=bgcounts/nrnewbins
        initParams = [np.max(ylist), Alpha_av, Theta_av]
        [Abinlist[measnr,binnr],Alphabinlist[measnr,binnr],Thetabinlist[measnr,binnr]]=MaxLikelihoodFit_c(tlist,ylist,istart,iend,bgcpb,initParams,plotbool=True)
    # calculate average lifetimes
#        plt.figure(15)
        plt.title(['binnr:', str(binnr),'Alpha: ', str(round(Alphabinlist[measnr,binnr],2)),'Theta: ', str(round(Thetabinlist[measnr,binnr]*180/np.pi))])
    
    
#        [tauavelist[measnr,binnr],Alist[measnr,binnr],_] = GetLifetime(np.append(binmicrotimes0,binmicrotimes1),dtmicro,dtmacro,dttau,-1,histbinmultiplier,ybg*onbincount)
        print('onbincount:',onbincount)
    #%% Fit polarization modulation to bins; compile for first time
tbinnr=0
maxbinning=10
Alpha=np.zeros(maxbinning)
Theta=np.zeros(maxbinning)
i=0
binning=np.zeros(maxbinning)
for i in range(maxbinning):
    histbinmultiplier=2**(i)
#    histbinmultiplier=20
    bgcounts=binwidth*80
    nrnewbins=int(dtmacro/(dtmicro*histbinmultiplier))
    #    print(20e-6/(dtmicro*histbinmultiplier))
    microtimes = np.append(microtimes0[limits0[tbinnr]:limits0[tbinnr+1]],microtimes1[limits1[tbinnr]:limits1[tbinnr+1]])
    [ylist,xlist] = np.histogram(microtimes,int(dtmacro/(dtmicro*histbinmultiplier)),[0,int(dtmacro/dtmicro)])    
    tlist = (xlist[:-1]+0.5*(xlist[1]-xlist[0]))*dtmicro*1e9 # convert x-axis to time in ns
    istart=0
    iend=max(tlist)
    bgcpb=bgcounts/nrnewbins
#    plt.plot(tlist,ylist)
    #    plt.clf()
    initParams = [np.max(ylist), Alpha_av, Theta_av]
#    MaxLikelihoodFunction_c = nb.jit(nopython=True)(MaxLikelihoodFunction) # implement c function
    result=MaxLikelihoodFit_c(tlist,ylist,istart,iend,bgcpb,initParams,plotbool=False)
    binning[i]=histbinmultiplier
    Alpha[i]=result[1]
    Theta[i]=result[2]*180/np.pi
    if result[1]>1:
        Alpha[i]=1/result[1]
        Theta[i]=result[2]*180/np.pi+90
plt9=plt.figure(9)
plt.plot(binning,Alpha)
plt.xlabel('Binning')
plt.ylabel('Estimated Alpha')
print(binning,Alpha,Theta)
#print(result)
#%% Loop
#Ampl=np.zeros(nrtbins);
#Alpha=np.zeros(nrtbins);
#Theta=np.zeros(nrtbins);
#histbinmultiplier=1
#bgcounts=binwidth*80 #assume 80 dark counts for now
#nrnewbins=int(dtmacro/(dtmicro*histbinmultiplier))
#istart=0
#iend=max(tlist)
#bgcpb=bgcounts/nrnewbins
#for tbinnr in range(nrtbins):
##    print(20e-6/(dtmicro*histbinmultiplier))
#    microtimes = np.append(microtimes0[limits0[tbinnr]:limits0[tbinnr+1]],microtimes1[limits1[tbinnr]:limits1[tbinnr+1]])
#    [ylist,xlist] = np.histogram(microtimes,int(dtmacro/(dtmicro*histbinmultiplier)),[0,int(dtmacro/dtmicro)])    
#    tlist = (xlist[:-1]+0.5*(xlist[1]-xlist[0]))*dtmicro*1e9 # convert x-axis to time in ns
##    plt.clf()
##    plt11 = plt.figure(11)
##    MaxLikelihoodFunction_c = nb.jit(nopython=True)(MaxLikelihoodFunction) # implement c function
#    result=MaxLikelihoodFit_c(tlist,ylist,istart,iend,bgcpb,plotbool=False)
#    Ampl[tbinnr]=result[0]
#    Alpha[tbinnr]=result[1]
#    Theta[tbinnr]=result[2]
#    print(str(round((tbinnr+1)/nrtbins*100))+"% done")
    
#%% Loop with parallel processing
Ampl=np.zeros(nrtbins-1);
Alpha=np.zeros(nrtbins-1);
Theta=np.zeros(nrtbins-1);
WeightedAlpha=np.zeros(nrtbins-1);
WeightedTheta=np.zeros(nrtbins-1);
Mactime=np.zeros(nrtbins-1);
histbinmultiplier=1
bgcounts=binwidth*80 #assume 80 dark counts for now
nrnewbins=int(dtmacro/(dtmicro*histbinmultiplier))
istart=0
iend=max(tlist)
bgcpb=bgcounts/nrnewbins
num_cores = multiprocessing.cpu_count()
initParams = [np.max(ylist), Alpha_av, Theta_av]

results = Parallel(n_jobs=-1, max_nbytes=None)(delayed(processInput)(tbinnr) for tbinnr in tqdm(range(nrtbins-1)))

# Plot results
for i in range(nrtbins-1):
    Mactime[i]=i*binwidth
    Ampl[i]=results[i][0]
    Alpha[i]=results[i][1]
    Theta[i]=results[i][2]*180/np.pi
    if results[i][1]>1:
        Alpha[i]=1/results[i][1]
        Theta[i]=results[i][2]*180/np.pi+90
    WeightedAlpha[i]=Alpha[i]*Ampl[i]
    WeightedTheta[i]=Theta[i]*Ampl[i]
plt.subplot(3,1,1)
plt.plot(Mactime,Ampl, 'b-')
plt.xlabel('time (s)')
plt.ylabel('Amplitude')

plt.subplot(3,1,2)
plt.plot(Mactime,Alpha, 'b-')
plt.xlabel('time (s)')
plt.ylabel('Anisotropy Alpha')

plt.subplot(3,1,3)
plt.plot(Mactime,Theta, 'b-')
plt.xlabel('time (s)')
plt.ylabel('Angle of max. Anisotropy (°)')

Correlation=np.corrcoef(Ampl,Alpha)
AlphaAverage=sum(WeightedAlpha)/sum(Ampl)
ThetaAverage=sum(Theta*Ampl)/sum(Ampl)
print(AlphaAverage)
#%% 
plt.plot(Ampl,Alpha,'.')