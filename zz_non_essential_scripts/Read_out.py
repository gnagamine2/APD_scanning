# -*- coding: utf-8 -*-
"""
Created on Thu Mar 28 16:44:49 2019

@author: rober
"""

import os, numpy as np, csv, matplotlib.pyplot as plt, scipy.optimize as opt, math, struct, binascii, gc, time, random
import multiprocessing
from operator import sub
from joblib import Parallel, delayed
import scipy, lmfit
from scipy.optimize import minimize # used for implementation of maximum likelihood exponential fit
from matplotlib import gridspec
from math import factorial
from scipy.stats import poisson
get_ipython().run_line_magic('matplotlib', 'auto')
import matplotlib as mpl
import pickle
import numba as nb
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
    
        while True:
            if f.read(1) == b'M':
                if f.read(18) == b'easDesc_Resolution': # recognize time unit entry
                    break
        f.read(21) # rest of the tag
        dtmicro = struct.unpack('d',f.read(8))[0]
        #print('Microtime unit:', dtmicro)
    
        while True:
            if f.read(1) == b'M':
                if f.read(24) == b'easDesc_GlobalResolution': # recognize time unit entry
                    break
        f.read(15) # rest of the tag
        dtmacro = struct.unpack('d',f.read(8))[0]
        #print('Macrotime unit:', dtmacro)
    
        while True:
            if f.read(1) == b'T':
                if f.read(23) == b'TResult_NumberOfRecords': # recognize number of records entry
                    break
        f.read(16) # rest of the tag
        nrrec = struct.unpack('q',f.read(8))[0] # extract number of records
        #print('Number of records in file:', nrrec)

        while True:
            if f.read(1) == b'H':
                if f.read(9) == b'eader_End':
                    #print('Header_End found')
                    break
        f.read(38) # rest of Header_End
#        nrrec=HHsettings["overallCounts"]
#        dtmicro=HHsettings["resolution"]*1e-12    #in s
#        dtmacro=1/HHsettings["syncRate"]    #in s
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


def importT3marked(filename):
    """Import ptu files from HydraHarp with Marker inputs for frame synchronization.

    Parameters
    ----------
    filename: String with the complete path to the file including file ending '.ptu'.

    Returns
    -------


    Notes
    -----
    Save T3 times directly in PicoHarp. Connect FIRE output of ANDOR EMCCD to Marker channel M1 and M2.
    Set active edge of channel M1 to "rising" and M2 to "falling".

    Discards all photon events before first frame and after last frame.

    References
    ----------


    Examples
    --------
    # >>> channel0, channel1, dt1, dt2 = importT3marked(filename) #FIXME: This doctest function does not work as intended

    """

    print('*** START READING PTU FILE ***')

    with open(filename, "rb+") as f:

        nrrec = HHsettings["overallCounts"]
        dtmicro = HHsettings["resolution"] * 1e-12  # in s
        dtmacro = 1 / HHsettings["syncRate"]  # in s
        macrotimes0 = np.zeros(nrrec, dtype='int64');
        microtimes0 = np.zeros(nrrec, dtype='int64');
        macrotimes1 = np.zeros(nrrec, dtype='int64');
        microtimes1 = np.zeros(nrrec, dtype='int64');
        macrotimesfireA = np.zeros(nrrec, dtype='int64');
        microtimesfireA = np.zeros(nrrec, dtype='int64');
        macrotimesfireB = np.zeros(nrrec, dtype='int64');
        microtimesfireB = np.zeros(nrrec, dtype='int64');
        overflows = 0
        nrphotons0 = 0
        nrphotons1 = 0
        nrfireA = 0
        nrfireB = 0
        prevchann = 0

        # initialize arrays
        macrotimes0 = np.zeros(nrrec, dtype='int64');
        microtimes0 = np.zeros(nrrec, dtype='int64');
        macrotimes1 = np.zeros(nrrec, dtype='int64');
        microtimes1 = np.zeros(nrrec, dtype='int64');
        macrotimesP = np.zeros(nrrec, dtype='int64');
        microtimesP = np.zeros(nrrec, dtype='int64');
        macrotimesL = np.zeros(nrrec, dtype='int64');
        microtimesL = np.zeros(nrrec, dtype='int64');

        overflows = 0  # overflow counter
        nrphotons0 = 0  # photon counter APD1
        nrphotons1 = 0  # photon counter APD2
        nrP = 0  # pixel tag counter
        nrL = 0  # line tag counter
        prevchann = 0  # remember last channel. Not sure what this is good for?

        for i in range(nrrec):
            entry = f.read(4)
            channel = struct.unpack("I", entry)[0] >> 25  # read channel number, first 7 bits
            if channel == 0:
                macrotime = (struct.unpack("I", entry)[0] & 0x3FF)
                macrotimes0[nrphotons0] = macrotime + 1024 * overflows
                microtime = ((struct.unpack("I", entry)[0] >> 10) & 0x7FFF)
                microtimes0[nrphotons0] = microtime
                nrphotons0 += 1
                prevchann = 0
            elif channel == 1:  # APD2: not implemented yet.
                macrotime = (struct.unpack("I", entry)[0] & 0x3FF)
                macrotimes1[nrphotons1] = macrotime + 1024 * overflows
                microtime = ((struct.unpack("I", entry)[0] >> 10) & 0x7FFF)
                microtimes1[nrphotons1] = microtime
                nrphotons1 += 1
                prevchann = 1
            elif channel == 127:
                nroverflows = (struct.unpack("I", entry)[0] & 0x3FF)
                overflows += nroverflows
                prevchann = 127
            elif channel == 65:
                macrotime = (struct.unpack("I", entry)[0] & 0x3FF)
                macrotimesP[nrP] = macrotime + 1024 * overflows
                microtime = ((struct.unpack("I", entry)[0] >> 10) & 0x7FFF)
                microtimesP[nrP] = microtime
                nrP += 1
            elif channel == 66:
                macrotime = (struct.unpack("I", entry)[0] & 0x3FF)
                macrotimesL[nrL] = macrotime + 1024 * overflows
                microtime = ((struct.unpack("I", entry)[0] >> 10) & 0x7FFF)
                microtimesL[nrL] = microtime
                nrL += 1
            else:
                print('bad channel:', channel)

    # cut emtpy entries at end of each storage array
    microtimes0 = microtimes0[:nrphotons0]
    macrotimes0 = macrotimes0[:nrphotons0]
    microtimes1 = microtimes1[:nrphotons1]
    macrotimes1 = macrotimes1[:nrphotons1]
    microtimesP = microtimesP[:nrP]
    macrotimesP = macrotimesP[:nrP]
    microtimesL = microtimesL[:nrL]
    macrotimesL = macrotimesL[:nrL]

    Texp = round(overflows * dtmacro * 1024,
                 0)  # macrotime values can store up to 2^10 = 1024 before overflow. overflows*1024*dtmacro gives experiment time [s]

    print('nrphotons0:', nrphotons0)
    print('nrphotons1:', nrphotons1)
    print('nrP:', nrP)
    print('nrL:', nrL)
    print('overflows:', overflows)

    print('averaged cps on det0 and det1:', np.array([nrphotons0, nrphotons1]) / Texp)
    print('experimental time in s:', Texp)
    print('*** FILE READING DONE ***')

    return [dtmicro, dtmacro, microtimes0, macrotimes0, microtimes1, macrotimes1, nrphotons0, nrphotons1, overflows,
            microtimesP, macrotimesP, nrP, microtimesL, macrotimesL, nrL]


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
    
    plt.xlabel('time (ns)')
    plt.ylabel('counts (a.u.)')
    p1, = plt.semilogy(tlist,ylist0+ylist1)
    p2, = plt.semilogy(tlist,ylist0)
    p3, = plt.semilogy(tlist,ylist1)
    plt.legend([p1,p2,p3], ["APD0 + APD1","APD0","APD1"])
           
    return(microtimes0-shift,microtimes1,tlist0,tlist1,dtmicro,dtmacro,tlist,ylist0+ylist1)

def CropData(limits0,limits1,microtimes0,microtimes1,macrotimes0,macrotimes1,tstart,tend,binwidth):
    binstart = int(tstart/binwidth)
    binend = int(tend/binwidth)

    limits0crop = np.array([limits0[i] for i in range(len(limits0)) if int(tstart/binwidth) <= i < int(tend/binwidth)])
    limits1crop = np.array([limits1[i] for i in range(len(limits1)) if int(tstart/binwidth) <= i < int(tend/binwidth)])
    microtimes0crop = np.array([microtimes0[i] for i in range(len(microtimes0)) if limits0crop[0] <= i <= limits0crop[-1]])
    microtimes1crop = np.array([microtimes1[i] for i in range(len(microtimes1)) if limits1crop[0] <= i <= limits1crop[-1]])
    macrotimes0crop = np.array([macrotimes0[i] for i in range(len(macrotimes0)) if limits0crop[0] <= i <= limits0crop[-1]])
    macrotimes1crop = np.array([macrotimes1[i] for i in range(len(macrotimes1)) if limits1crop[0] <= i <= limits1crop[-1]])

    Texp = tend-tstart

    return(microtimes0crop,microtimes1crop,macrotimes0crop,macrotimes1crop,limits0crop,limits1crop,Texp)

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

def GetVoltageLimits(times,macrotimesfireA): # finds the index of the first photon in each bin. times in [microtimeunits]. macrotimes in [macrotimeunits]
    
    nrVbins = len(macrotimesfireA)
    Vlimits = np.full(nrVbins,len(times))
    counter,i = 0,0
    while i < len(times):
        while (counter < nrVbins-1) & (times[i]*dtmicro/dtmacro > macrotimesfireA[counter]):
            Vlimits[counter] = i
            counter += 1
        i += 1
    
    return(Vlimits)

def BinVoltage(microtimes0,times0,microtimes1,times1,macrotimesfireA):
    nrVbins = len(macrotimesfireA)

    limitsfire0 = GetVoltageLimits(times0,macrotimesfireA)
    limitsfire1 = GetVoltageLimits(times1,macrotimesfireA)

    # find photons in V=on and V=off periods
    Vonphotonlist0 = np.array([np.arange(limitsfire0[binnr],limitsfire0[binnr+1]) for binnr in range(nrVbins-1) if np.mod(binnr,2)==0])
    Voffphotonlist0 = np.array([np.arange(limitsfire0[binnr],limitsfire0[binnr+1]) for binnr in range(nrVbins-1) if np.mod(binnr,2)!=0])
    Vonphotonlist1 = np.array([np.arange(limitsfire1[binnr],limitsfire1[binnr+1]) for binnr in range(nrVbins-1) if np.mod(binnr,2)==0])
    Voffphotonlist1 = np.array([np.arange(limitsfire1[binnr],limitsfire1[binnr+1]) for binnr in range(nrVbins-1) if np.mod(binnr,2)!=0])

    # reduce arrays
    Vonphotonlist0 = np.concatenate(Vonphotonlist0).ravel()
    Voffphotonlist0 = np.concatenate(Voffphotonlist0).ravel()
    Vonphotonlist1 = np.concatenate(Vonphotonlist1).ravel()
    Voffphotonlist1 = np.concatenate(Voffphotonlist1).ravel()

    # get V=on microtimes
    Vonmicrotimes0 = np.array([microtimes0[i] for i in Vonphotonlist0])
    Vonmicrotimes1 = np.array([microtimes1[i] for i in Vonphotonlist1])
    Vontimes0 = np.array([times0[i] for i in Vonphotonlist0])
    Vontimes1 = np.array([times1[i] for i in Vonphotonlist1])
    #plt.title('V on')
    #Vontauave = GetLifetime(np.append(Vonmicrotimes0,Vonmicrotimes1),dtmicro,dtmacro,200e-9,-1)

    # get V=off microtimest
    Voffmicrotimes0 = np.array([microtimes0[i] for i in Voffphotonlist0])
    Voffmicrotimes1 = np.array([microtimes1[i] for i in Voffphotonlist1])
    Vofftimes0 = np.array([times0[i] for i in Voffphotonlist0])
    Vofftimes1 = np.array([times1[i] for i in Voffphotonlist1])
    #plt.title('V off')
    #Vofftauave = GetLifetime(np.append(Voffmicrotimes0,Voffmicrotimes1),dtmicro,dtmacro,200e-9,-1)

    return(Vonmicrotimes0,Voffmicrotimes0,Vontimes0,Vofftimes0,Vonmicrotimes1,Voffmicrotimes1,Vontimes1,Vofftimes1)

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
                                
    g2tlist = np.arange(-g2res*dtmicro*(nrbins-0.5),g2res*dtmicro*nrbins,g2restime)*1e9
    plt.plot(g2tlist,g2)
    #plt.plot(g2tlist,g2B)
    #plt.plot(g2tlist,g2C)
    plt.title('g(2) correlation')
    plt.xlabel('delay (ns)')
    plt.ylabel('occurence (a.u.)')
    plt.ylim([0,max(g2)])
    plt.show()

    return(g2tlist,g2,g2restime,nrbins)

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
    ydata = ylist[istart:iend]
    xdata = tlist[istart:iend]

    # do calculations
    initParams = [np.max(ydata), 25] #initial guess for A and tau
    results = minimize(MaxLikelihoodFunction, initParams, args=(xdata,ydata,bgcpb),method='Nelder-Mead') # minimize the negative of the maxlikelihood function instead of maximimizing
    Aest = results.x[0] # get results of fit, A
    tauest = results.x[1] # get results of fit, tau

    if plotbool == True:
        yest = np.array([Aest*np.exp(-(xdata[i]-xdata[0])/tauest)+bgcpb for i in range(len(xdata))])
        plt.semilogy(tlist,ylist,'.',xdata,yest,[xdata[1],xdata[-1]],[bgcpb,bgcpb],'k--')
        plt.show()        
    
    return(tauest,Aest)

def MaxLikelihoodFunction(params,xdata,ydata,const):
    # max likelihood function for A*exp(-t/tau), needed in function MakLikelihoodFit
    # params = [A,tau]
    A = params[0]
    tau = params[1]   
    E = 0;
    for i in range(len(xdata)):
        E = E + ydata[i]*np.log(A*np.exp(-(xdata[i]-xdata[0])/tau)+const)-(A*np.exp(-(xdata[i]-xdata[0])/tau)+const)
        
    return(-E) # This function needs to be MINIMIZED (because of the minus sign) to have the maximum likelihood fit!

def MaxLikelihoodFit_c(tlist,ylist,istart,iend,bgcpb,plotbool=False):
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
    ydata = ylist[istart:iend]
    xdata = tlist[istart:iend]

    # do calculations
    initParams = [np.max(ydata), 25] #initial guess for A and tau
    results = minimize(MaxLikelihoodFunction_c, initParams, args=(xdata,ydata,bgcpb),method='Nelder-Mead') # minimize the negative of the maxlikelihood function instead of maximimizing
    Aest = results.x[0] # get results of fit, A
    tauest = results.x[1] # get results of fit, tau

#    if plotbool == True:
#        yest = np.array([Aest*np.exp(-(xdata[i]-xdata[0])/tauest)+bgcpb for i in range(len(xdata))])
#        plt.semilogy(tlist,ylist,'.',xdata,yest,[xdata[1],xdata[-1]],[bgcpb,bgcpb],'k--')
#        plt.show()        


    if plotbool == True:
        yest = np.array([Aest*np.exp(-(xdata[i]-xdata[0])/tauest)+bgcpb for i in range(len(xdata))])
        plt.figure()
        plt.plot(tlist,ylist,'.',xdata,yest,[xdata[1],xdata[-1]],[bgcpb,bgcpb],'k--')
        plt.xlim([xdata[1],xdata[-1]])
        plt.show()        
        
    return(tauest,Aest)

def CorrelateSlices(slicewidth,tlist0,tlist1,limits0,limits1,intlist,binwidth,tmicro,tmacro):
    base = 1.2
    
    intlimits = np.arange(slicewidth*int(min(intlist)/slicewidth),slicewidth*int(max(intlist)/slicewidth),slicewidth)
    slices = np.array([np.array([i for i,x in enumerate(intlist) if intlimits[j] < x <= intlimits[j+1]]) for j in range(len(intlimits)-1)])

    correlations = []; correlations2 = []; norm1list = []; norm2list = []; norm3list = []; norm4list = [];
    tbin = binwidth/tmicro; tbinp = int(tbin) # bin width in units of the fundamental microtime unit
    Tbin = binwidth/tmacro; Tbinp = int(Tbin) # bin width in units of the fundamental macrotime unit
    resratio = int(tmacro/tmicro)

    #logbins = np.append(np.insert(base**np.arange(0,int(np.log(tmax)/np.log(base)+1),dtype='int64'),0,0),tmax).astype(np.int64)
    logbins = np.insert(int(0.5*resratio)+(resratio*np.append(np.insert(base**np.arange(0,int(np.log(Tbin)/np.log(base)+1),dtype='int64'),0,0),Tbin).astype(np.int64)),0,0)
    logbins = np.unique(logbins)
    logbins = np.append(logbins,[int(x) for x in np.divide([2e-7,6e-7,2e-6,6e-6,2e-5,6e-5,2e-4,6e-4,2e-3],tmicro)])

    for slicenr in (np.arange(len(slices),0,-1)-1): # count down from the highest intensity
        lower0 = np.array([limits0[x] for x in slices[slicenr]]) # list of first photons on detector 0 in the bins belonging to the slice
        upper0 = np.array([limits0[x+1] for x in slices[slicenr]])
        lower1 = np.array([limits1[x] for x in slices[slicenr]])
        upper1 = np.array([limits1[x+1] for x in slices[slicenr]])

        total0 = 0; total1 = 0;
        histavg = np.zeros(len(logbins)-1,dtype='float');get
        histavg2 = np.zeros(len(logbins)-1,dtype='float');
        calctime0 = time.time()
    
        for i in range(len(lower0)):
            selection0 = tlist0[lower0[i]:upper0[i]]
            selection1 = tlist1[lower1[i]:upper1[i]]
            total0 += upper0[i]-lower0[i]#len(selection0)
            total1 += upper1[i]-lower1[i]#len(selection1)
                
            if len(selection0) == 0 or len(selection1) == 0:
                continue
               
            hist = np.zeros(len(logbins)-1,dtype='float');
            for counterA in range(len(selection1)):
                res = selection0.searchsorted(logbins+selection1[counterA])
                hist += (res[1:]-res[:-1])
            
            #hist = np.zeros(len(logbins)-1,dtype='float');
            for counterB in range(len(selection0)):
                res = selection1.searchsorted(logbins+selection0[counterB])
                hist += (res[1:]-res[:-1])
            
            histavg += hist/(len(selection1)+len(selection0))*(2*Tbinp)
            histavg2 += hist
        
        if slicenr % 5 == 0:
            print('Analysis time for slice #',slicenr,':', (time.time() - calctime0), 's')

        norm = 0.5*(logbins[:-1] - logbins[1:])*(logbins[:-1] + logbins[1:] - 2*tbin - 1)
        tlist = 0.5*(logbins[:-1] + logbins[1:])*data[0]
        norm1 = total1
        norm2 = len(slices[slicenr])
        norm3 = total0
        norm4 = tbinp
    
        correlations.append(np.divide(histavg,norm/norm[1]))
        correlations2.append(np.divide(histavg2/total1,norm/norm[1]))
        norm1list.append(norm1); norm2list.append(norm2); norm3list.append(norm3); norm4list.append(norm4)
    
    return(tlist,correlations,correlations2,norm1list,norm2list,norm3list,norm4list)

def DecaySlices(slicewidth,micro0,micro1,limits0,limits1,inttrace,dt,tmicro,tmacro):
    intlimits = np.arange(slicewidth*int(min(inttrace)/slicewidth),slicewidth*int(max(inttrace)/slicewidth),slicewidth)
    slices = np.array([np.array([i for i,x in enumerate(inttrace) if intlimits[j] < x <= intlimits[j+1]]) for j in range(len(intlimits)-1)])

    decaycurves = np.zeros((len(slices),int(tmacro/dt)+1))

    for slicenr in (np.arange(len(slices),0,-1)-1): # count down from the highest intensity
        lower0 = np.array([limits0[x] for x in slices[slicenr]]) # list of first photons on detector 0 in the bins belonging to the slice
        upper0 = np.array([limits0[x+1] for x in slices[slicenr]])
        lower1 = np.array([limits1[x] for x in slices[slicenr]])
        upper1 = np.array([limits1[x+1] for x in slices[slicenr]])

        lower0 = [int(x) for x in lower0]
        lower1 = [int(x) for x in lower1]
        upper0 = [int(x) for x in upper0]
        upper1 = [int(x) for x in upper1]

        ylist=np.zeros(int(tmacro/dt))
        for j in range(len(lower0)):
            ylist += (np.histogram(micro0[lower0[j]:upper0[j]],int(tmacro/dt),[0,int(tmacro/dt)*int(dt/tmicro)])[0])
            ylist += (np.histogram(micro1[lower1[j]:upper1[j]],int(tmacro/dt),[0,int(tmacro/dt)*int(dt/tmicro)])[0])
        
        ylist = np.append(ylist,len(slices[slicenr]))
        
        decaycurves[slicenr]=ylist
        
    return(decaycurves)

def doublepoisson(k,a,b,c,d):
    return a*poisson.pmf(k.astype(int),b)+c*poisson.pmf(k.astype(int),d)

def singlepoisson(k,a,b):
    return a*poisson.pmf(k.astype(int),b)

def singlepoissonfit(xdata,ydata,params0,plotbool=False):
    # do calculations
    [popt, _]= scipy.optimize.curve_fit(singlepoisson, xdata, ydata,p0 = params0)
    if plotbool == True:
        plt.plot(xdata,ydata,'k.')
        plt.plot(xdata,singlepoisson(xdata,*popt))
    return popt

def doublepoissonfit(xdata,ydata,params0,plotbool=False):
    # do calculations
    [popt, _]= scipy.optimize.curve_fit(doublepoisson, xdata, ydata,p0 = params0)
    if plotbool == True:
        plt.plot(xdata,ydata,'k.')
        plt.plot(xdata,doublepoisson(xdata,*popt))
    return popt

def readspec(filename,dy=0): # read andor's spectral data (ascii)
    #dy gives width in y-direction for spectral averaging
    
    data = np.genfromtxt(filename, delimiter=',')
    data[:,0]

    # find length of spectra
    i=0
    while data[i,0]<data[i+1,0]: i+=1 # find end of first spectrum by checking when maximum lambda is reached
    xlen = i+1
    ylen = data.shape[1]-1
    nrimgs = int(len(data)/xlen)

    # sum y
    ysums = np.array([np.sum(data[:,i]) for i in range(1,ylen)])
    ymaxind = np.argmax(ysums)+1 #+1 because zero are lambda values!
    yminind = 1 # BACKGROUND IS HARD-COPIED HERE!!!!

    # get lambda
    lambd = data[0:xlen,0]

    # get spec
    spec = np.zeros(shape=(xlen,nrimgs),dtype=float)

    for imgnr in range(nrimgs):
        tmpspec = np.mean(data[imgnr*xlen:xlen+imgnr*xlen,ymaxind-dy:ymaxind+dy],axis=1)
        bg = np.mean(data[0:xlen,yminind:yminind+3],axis=1)
        spec[:,imgnr] = tmpspec - bg

    return(lambd,spec,xlen,nrimgs)

def Gauss(x, a, x0, sigma): # gaussian function
    return a * np.exp(-(x - x0)**2 / (2 * sigma**2))

def GaussFit(x,y,plotbool=False): # fits gaussian to given curve
    # calculate starting values
    mean = sum(x * y) / sum(y)
    sigma = np.sqrt(sum(y * (x - mean)**2) / sum(y))

    # fitting
    try:
        popt,pcov = scipy.optimize.curve_fit(Gauss, x, y, p0=[max(y), mean, sigma],maxfev = 10000) 
    except RuntimeError:
        #print("Gaussian fit failed. Replaced data with NaN")
        popt= [float('nan'),float('nan'),float('nan')]

    # plot if plotbool == True
    if plotbool == True:
        plt.plot(x, y, 'b+:', label='data')
        plt.plot(x, Gauss(x, *popt), 'r-', label='fit')
        plt.legend()
        plt.show()

    return(popt) # [maximum intensity, maximum lambda, sigma]

# In[25]:

### GENERAL EVALUATION OF TTTR data     

# parameters (use forward slashes!)
#folder = 'C:/Users/omel-acton/Desktop/gitfolder/apd_scanning_imaging/output/'
#folder = 'C:/Users/rober/Documents/Doktorat\Projects/Linear_absorption_anisotropy/Andor Spectrometer/20190416_rods_polymer/scans/measurements/'
#folder = 'C:/Users/rober/Documents/Doktorat/Projects/SingleParticle_PLE/20190911_PLE_PM111/'
#namelist = ['Dot1_405']

#folder = "/Users/JulianSanten/GIT/OMEL_GIT/Data/22-09-01_ScanningStageTTTRCodeTest3/measurements/"
#namelist = ["Test_AP_4_17_10-6_RepRate2,5x10^-6_x_274_y60_z94_range30_pts150_Particle_0_x281y_70_tttrmode"]

#folder = "/Users/JulianSanten/GIT/OMEL_GIT/Data/22-09-01_ScanningStageTTTRCodeTest3/"
#namelist = ["Test1_x_266_y74_z93_range30_pts150_tttrmode"]

folder = "/Users/JulianSanten/GIT/OMEL_GIT/Data/22-09-01_ScanningStageTTTRCodeTest3/"
namelist = ["Test_AP_4_17_10-6_RepRate2,5x10^-6_x_274_y60_z94_range30_pts150_tttrmode"]

#settingsfile=namelist[0].rsplit('_',1)[0]+'_HHsettings'
#HHsettings=load_obj(settingsfile, folder )
#MCLdata=load_obj("Sample1_x174_y187_z11_range20_pts100_MCLdata", folder )

Texp = 240.0 # total time [s] of measurements
binwidth = 0.01 # width of time bins to create intensity trace [s]
nrbins = 4 # intensity bin to separate intensity trace in (nrbins) to separately evaluate different intensity levels
histbinmultiplier = 1 #(histbinmultiplier)*dtmicro = bin size for lifetime histogram
dttau = 50e-9 # end time of lifetime fit [s]. Note: this is not the width but the end time! TODO: NEED TO CHANGE LATER!

savebool = True;
# subtractbg = True; # not used?! 
singlebintaubool = True;
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
tauavelist = np.full((nrmeas,nrbins),0,dtype=float) # init list to save lifetimes of each intensity bin
Alist = np.full((nrmeas,nrbins),0,dtype=float) # init list to save maximum histogram amplitudes (at t=t0)
taubinlist = np.zeros(nrtbins-1,dtype=float) # init list to save lifetimes of each time bin (binwidth)
intbinlist = np.zeros(nrtbins-1,dtype=float) # init list to save intensity per bin

for measnr in range(nrmeas): #for each measurement file

    # load data
    #data = ImportT3(folder + namelist[measnr] + '.ptu') # import data from .ptu file
    #data = ImportT3(folder + namelist[measnr] + '.out') # import data from .ptu file

    #%%
    [dtmicro, dtmacro, microtimes0, macrotimes0, microtimes1, macrotimes1, nrphotons0,nrphotons1,overflows,microtimesP,macrotimesP,nrP,microtimesL,macrotimesL,nrL] = importT3marked(folder + namelist[measnr] + '.out') # import data from .ptu file




    #%%
    Texp = round(data[8]*data[1]*1024,0) # macrotime values can store up to 2^10 = 1024 before overflow. overflows*1024*dtmacro gives experiment time [s]
    print('averaged cps on det0 and det1:',np.array(data[6:8])/Texp)
    print('experimental time in s:',Texp)

    # compensate detector offset
    plt11 = plt.figure(11)
    plt.title(str(measnr))
    [microtimes0,microtimes1,times0,times1,dtmicro,dtmacro,decaytlist,decayylist] = ShiftPulsedData(data[2],data[4],data[3],data[5],data[0],data[1])
    plt.show()
    if savebool == True:        
        np.save(folder + 'output/' + namelist[measnr] + '_decay_tlist', decaytlist)
        np.save(folder + 'output/' + namelist[measnr] + '_decay_ylist', decayylist)   
    tmax = decaytlist[decayylist.argmax()] # find histogram time with max photon count [ns]
    
    # get number of first photon per bin for both detectors (APD1 = 0, APD2 = 1)
    limits0 = HistPhotons(times0*dtmicro,binwidth,Texp) #gives index of first photon of det0 in each bin
    limits1 = HistPhotons(times1*dtmicro,binwidth,Texp)

    # make an intensity trace and find Imax
    inttrace = MakeIntTrace(limits0,limits1,binwidth,Texp)
    Imax = np.max(inttrace) # search overall maximum intensity per intensity bin
    Ilims = np.arange(0,Imax*0.95+1,Imax*0.95/nrbins,dtype=float) # separate intensities into (nrbins) intensity bins
    Ilimslist[measnr,:] = np.array([Ilims[binnr]*0.5+Ilims[binnr+1]*0.5 for binnr in range(nrbins)]) # get centered intensity value per intensity bin
    
    # get G2
    if g2bool == True:
        [g2tlist,g2ylist,_,_] = MakeG2(times0,times1,dtmicro,8e-9,150)
        if savebool == True:        
            np.save(folder + 'output/' + namelist[measnr] + '_G2_tlist', g2tlist)
            np.save(folder + 'output/' + namelist[measnr] + '_G2_ylist', g2ylist)
            
    for binnr in range(nrbins): # for each intensity bin
        # select photons in bins with intensity within Ilims(binnr) and Ilims(binnr+1)
        [binmicrotimes0,bintimes0,binmicrotimes1,bintimes1,onbincount] = SliceHistogram(microtimes0,times0,limits0,microtimes1,times1,limits1,dtmicro,dtmacro,Ilims[binnr],Ilims[binnr+1])  

        # calculate average lifetimes
        plt.figure(15)
        plt.title(['binnr:', str(binnr)])

        if binnr == 0: # use bin with lowest intensity to get background counts (idea: lowest intensity usually means fastest decay. chance is high to avoid build-up of photons)
            [tauavelist[measnr,binnr],Alist[measnr,binnr],ybg] = GetLifetime(np.append(binmicrotimes0,binmicrotimes1),dtmicro,dtmacro,dttau,-1,histbinmultiplier,0)
            bgcts = ybg*int(dtmacro/(dtmicro*histbinmultiplier)) #ybg is the background counts per histogram bin of lifetime binning. multiply with number of hist bins (int(dtmacro/(dtmicro*histbinmultiplier))) to get total nr of bgcounts
            bgcps = bgcts/(onbincount*binwidth) # divide by nr of bins of this slice times the binwidth to get cps
            ybg=ybg/onbincount #normalize to nr of time bins within slice
            print('bgcts:',bgcts)
            print('bgcps:',bgcps)
            print('onbincount:',onbincount)

            # convert data for fitting with max likelihood
            bgcpsfit = bgcps
            bgcpbfit = bgcps*binwidth/int(dtmacro/(dtmicro*histbinmultiplier)); # counts per fitting bin (see below!)
            print('bgcpbfit',bgcpbfit)
        else: # for all other bins just excract lifetimes
            [tauavelist[measnr,binnr],Alist[measnr,binnr],_] = GetLifetime(np.append(binmicrotimes0,binmicrotimes1),dtmicro,dtmacro,dttau,-1,histbinmultiplier,ybg*onbincount)
            print('onbincount:',onbincount)
            
    # get lifetimes for each time bin
    if singlebintaubool == True:
        MaxLikelihoodFunction_c = nb.jit(nopython=True)(MaxLikelihoodFunction)
        for tbinnr in tqdm(range(nrtbins-1)):
            microtimes = np.append(microtimes0[limits0[tbinnr]:limits0[tbinnr+1]],microtimes1[limits1[tbinnr]:limits1[tbinnr+1]])
            [ylist,xlist] = np.histogram(microtimes,int(dtmacro/(dtmicro*histbinmultiplier)),[0,int(dtmacro/dtmicro)])
            tlist = (xlist[:-1]+0.5*(xlist[1]-xlist[0]))*dtmicro*1e9 # convert x-axis to time in ns
            idxstart, = np.where(tlist == next( (x for x in tlist if x>tmax), 1)) # find index of dttau
            istart = idxstart.item() # convert index to normal scalar
            idxend, = np.where(tlist == next( (x for x in tlist if x>dttau*1e9), len(tlist))) # find index of dttau
            iend = idxend.item()+istart # convert index to normal scalar
            bgcpbfit = bgcps*binwidth/int(dtmacro/(dtmicro*histbinmultiplier)); # counts per fitting bin
            [tau,A] = MaxLikelihoodFit_c(tlist,ylist,istart,iend,bgcpbfit,False)
            taubinlist[tbinnr]=tau
            intbinlist[tbinnr]=len(microtimes)
#            if np.mod(tbinnr,5)==0: #for every 1000th bin output a plot of the fit (just to visualize the process)
#                print('maxlikelihood tbinnr:',tbinnr)
#                [tau,A] = MaxLikelihoodFit(tlist,ylist,istart,iend,bgcpbfit,True) # do again just for plotting
        
    # write out lifetimes and intensities
        if savebool == True:        
            np.save(folder + 'output/' + namelist[measnr] + '_FLID_taubin', taubinlist)
            np.save(folder + 'output/' + namelist[measnr] + '_FLID_intbin', intbinlist)   

# In[34]:

### PLOT DATA
            
# input
savebool = False # safe figures to folder /output/
intmin = 0
intmax = 500
intbins = 50
taumin = 0
taumax = 60
taubins = 61
nrmeas = len(namelist)


# initialize FLID array
FLID = np.ndarray(shape=(intbins,taubins),dtype=int) # init list to save intensity per bin

fig2 = plt.figure(11,figsize=(8,4))
plt.xlabel('time (ns)')
plt.ylabel('normalized count rate')
for measnr in range(nrmeas):
    decaytlist = np.load(folder + 'output/' + namelist[measnr] + '_decay_tlist.npy')    
    decayylist = np.load(folder + 'output/' + namelist[measnr] + '_decay_ylist.npy')
    plt.semilogy(decaytlist,decayylist/np.max(decayylist))

    
for measnr in range(nrmeas):
    intbinlist = np.load(folder + 'output/' + namelist[measnr] + '_FLID_intbin.npy')
    taubinlist = np.load(folder + 'output/' + namelist[measnr] + '_FLID_taubin.npy')
    timelist = binwidth*np.arange(len(intbinlist))

    fig = plt.figure(10,figsize=(18,10))
    gs = gridspec.GridSpec(3,3,width_ratios=[1.5,3,1.5])
    
    # FLID MAP
    plt.subplot(gs[0])
    plt.title('',y=0.85)
    plt.ylabel('emission intensity (cts / %i ms)' %(binwidth*1e3))
    plt.xlabel('PL lifetime (ns)')
    [H,intedges,tauedges,img] = plt.hist2d(taubinlist,intbinlist,[taubins,intbins],[[taumin,taumax],[intmin,intmax]],cmin=0,cmap='RdPu',norm=mpl.colors.LogNorm())
    plt.colorbar(img)
        
    # Intensity trace
    plt.subplot(gs[1])
    plt.title(namelist[measnr])
    plt.plot(timelist,intbinlist,'-',linewidth=1.5,color='k')
    plt.xlabel('time (s)')
    plt.ylabel('counts / %i ms' %(binwidth*1e3))
    plt.xlim([0,Texp])
    plt.ylim([intmin,intmax])
    
    # Intensity histogram
    inthistogram = np.histogram(intbinlist,intbins,[0,int(np.max(intbinlist))])
    plt.subplot(gs[2])
    xdata = 0.5*(inthistogram[1][:-1]+inthistogram[1][1:])
    plt.plot(inthistogram[0],xdata,color='k',linewidth=1.5)
    plt.xlabel('occurrence')
    plt.ylabel('counts / %i ms' %(binwidth*1e3))
    plt.ylim([intmin,intmax])
    
    # poisson fits
#    popt = singlepoissonfit(xdata[int(intbins/2):-1],inthistogram[0][int(intbins/2):-1],[np.sqrt(0.9*np.max(intbinlist)),0.9*np.max(intbinlist)]);
#    plt.plot(singlepoisson(xdata,*popt),xdata)
#    popt = singlepoissonfit(xdata[0:int(intbins/2)],inthistogram[0][0:int(intbins/2)],[np.sqrt(0.12*np.max(intbinlist)),0.12*np.max(intbinlist)]);
#    plt.plot(singlepoisson(xdata,*popt),xdata)
    
    #g2 correlation
    if g2bool == True:
        g2tlist = np.load(folder + 'output/' + namelist[measnr] + '_G2_tlist.npy')
        g2ylist = np.load(folder + 'output/' + namelist[measnr] + '_G2_ylist.npy')
        plt.subplot(gs[3])
        plt.plot(g2tlist,g2ylist/np.max(g2ylist),color='k',linewidth=1.5)
        plt.xlabel('delay (ns)')
        plt.ylabel('occurence (a.u.)')
        plt.ylim([0,1])
    else:
        plt.subplot(gs[3])
        plt.xlabel('time (ns)')
        plt.ylabel('normalized count rate')
        for measnr in range(nrmeas):
            decaytlist = np.load(folder + 'output/' + namelist[measnr] + '_decay_tlist.npy')    
            decayylist = np.load(folder + 'output/' + namelist[measnr] + '_decay_ylist.npy')
        plt.semilogy(decaytlist,decayylist/np.max(decayylist),'k')
    
    # Lifetime trace
    plt.subplot(gs[4])
    plt.plot(timelist,taubinlist,'-',linewidth=1.5,color='k')
    plt.xlabel('time (s)')
    plt.ylabel('lifetime (ns)')
    plt.xlim([0,Texp])
    plt.ylim([taumin,taumax])
    
    # lifetime histogram    
    tauhistogram = np.histogram(taubinlist,taubins,[taumin,taumax])
    xdata = 0.5*(tauhistogram[1][:-1]+tauhistogram[1][1:])
    plt.subplot(gs[5])
    plt.plot(tauhistogram[0],0.5*(tauhistogram[1][:-1]+tauhistogram[1][1:]),color='k',linewidth=1.5)
    plt.xlabel('occurrence')
    plt.ylabel('lifetime (ns)')
    plt.ylim([taumin,taumax])
    
    # poisson fit
#    popt = singlepoissonfit(xdata,tauhistogram[0],[np.sqrt(0.95*np.max(taubinlist)),0.95*np.max(taubinlist)]);
#    plt.plot(singlepoisson(xdata,*popt),xdata)
        
    if savebool == True:
        plt.savefig(folder + 'output/' + namelist[measnr] + '.pdf', bbox_inches='tight',format='pdf')
        plt.savefig(folder + 'output/' + namelist[measnr], bbox_inches='tight')
    plt.show()