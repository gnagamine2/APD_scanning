# -*- coding: utf-8 -*-
"""
Created on Tue Jan 15 17:50:02 2019
TTTR analysis code with marker from Felipe
@author: robert
"""

import struct
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import scipy, lmfit # used for implementation of least-squares fit
from scipy.optimize import minimize # used for implementation of maximum likelihood exponential fit
import time # used to count run time of loops
import numba as nb
from matplotlib.colors import LogNorm

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
    >>> channel0, channel1, dt1, dt2 = importT3marked(filename)
    
    """
    
    print('*** START READING PTU FILE ***')

    with open(filename, "rb") as f:

        while True:
            if f.read(1) == b'M':
                if f.read(18) == b'easDesc_Resolution': # recognize time unit entry
                    break
        f.read(21) # rest of the tag
        dtmicro = struct.unpack('d',f.read(8))[0]
        print('Microtime unit:', dtmicro)

        while True:
            if f.read(1) == b'M':
                if f.read(24) == b'easDesc_GlobalResolution': # recognize time unit entry
                    break
        f.read(15) # rest of the tag
        dtmacro = struct.unpack('d',f.read(8))[0]
        print('Macrotime unit:', dtmacro)

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

        # initialize arrays
        macrotimes0 = np.zeros(nrrec,dtype='int64');
        microtimes0 = np.zeros(nrrec,dtype='int64');
        macrotimes1 = np.zeros(nrrec,dtype='int64');
        microtimes1 = np.zeros(nrrec,dtype='int64');
        macrotimesP = np.zeros(nrrec,dtype='int64');
        microtimesP = np.zeros(nrrec,dtype='int64');
        macrotimesL = np.zeros(nrrec,dtype='int64');
        microtimesL = np.zeros(nrrec,dtype='int64');

        overflows = 0 # overflow counter
        nrphotons0 = 0 # photon counter APD1
        nrphotons1 = 0 # photon counter APD2
        nrP = 0 # pixel tag counter
        nrL = 0 # line tag counter
        prevchann = 0 # remember last channel. Not sure what this is good for?

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
            elif channel == 1: # APD2: not implemented yet.
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
                macrotimesP[nrP] = macrotime + 1024*overflows
                microtime = ((struct.unpack("I",entry)[0] >> 10) & 0x7FFF)
                microtimesP[nrP] = microtime
                nrP += 1
            elif channel == 66:
                macrotime = (struct.unpack("I",entry)[0] & 0x3FF)
                macrotimesL[nrL] = macrotime + 1024*overflows
                microtime = ((struct.unpack("I",entry)[0] >> 10) & 0x7FFF)
                microtimesL[nrL] = microtime
                nrL += 1
            else:
                print('bad channel:',channel)

    # cut emtpy entries at end of each storage array
    microtimes0 = microtimes0[:nrphotons0]
    macrotimes0 = macrotimes0[:nrphotons0]
    microtimes1 = microtimes1[:nrphotons1]
    macrotimes1 = macrotimes1[:nrphotons1]
    microtimesP = microtimesP[:nrP]
    macrotimesP = macrotimesP[:nrP]
    microtimesL = microtimesL[:nrL]
    macrotimesL = macrotimesL[:nrL]

    Texp = round(overflows*dtmacro*1024,0) # macrotime values can store up to 2^10 = 1024 before overflow. overflows*1024*dtmacro gives experiment time [s]

    print('nrphotons0:',nrphotons0)
    print('nrphotons1:',nrphotons1)
    print('nrP:',nrP)
    print('nrL:',nrL)
    print('overflows:',overflows)
    
    print('averaged cps on det0 and det1:',np.array([nrphotons0,nrphotons1])/Texp)
    print('experimental time in s:',Texp)
    print('*** FILE READING DONE ***')
    
    return [dtmicro, dtmacro, microtimes0, macrotimes0, microtimes1, macrotimes1, nrphotons0,nrphotons1,overflows,microtimesP,macrotimesP,nrP,microtimesL,macrotimesL,nrL]

def FindFirstPhotonAfterTrigger(macrotimes,macrotimesT): 
    # finds the index of the first photon i of the array macrotimes after each event of macrotimesT (e.g. of trigger). 
    # gives out the index of the first event in array macrotimes that happens after each trigger event in array macrotimesT
    # usually macrotimes = macrotimes of photons
    # macrotimesT = macrotimes of trigger
    # output firstPhot gives array with indices of the first photons per macrotimesT bin. Length is the same as len(macrotimesT)
    
    nrbins = len(macrotimesT)
    nrphot = len(macrotimes)
    firstPhot = np.full(nrbins,nrphot) #initialize array and assign nr of last photon to all elements (for safety)
    counter,i = 0,0
    while counter < nrbins and i < nrphot: # for all trigger signals and photons
        while macrotimes[i] > macrotimesT[counter] : # as long as photon i arrives after trigger counter do:
            firstPhot[counter] = i # store
            counter += 1 # # increase trigger
            if counter == nrbins:
                break
        if counter == nrbins:
                break
        i += 1 # as long as photon arrived before trigger, go to next photon
    
    return(firstPhot)

def GetLifetime(microtimes,dtmicro,dtmacro,dtfit,tstart=-1,binwidth=1,ybg=0,plotbool=True,method='WLS'): 
    # microtimes = microtimes array with photon events
    # dtfit is the time interval considered for the fit [s], tstart [s] is the starting point of the fit within the histogram. If set to -1 it starts at the time with the highest intensity.
    # binwidth is a multiplier. actual binwidth is given as binwidth*dtmicro[s]
    # ybg is the background considered for the fit (CHECK UNITS!!). If set to -1 --> try to estimate background based on last bins. set to 0 --> no background subtraction
    # plotbool: plot histogram with fit
    
    [ylist,xlist] = np.histogram(microtimes,int(dtmacro/(dtmicro*binwidth)),[0,int(dtmacro/dtmicro)])
    tlist = (xlist[:-1]+0.5*(xlist[1]-xlist[0]))*dtmicro*1e9
    
    istart = int(tstart/dtmicro) #find index of maximum element in ylist
    if istart < 0:
        istart = ylist.argmax()
    iend = istart + int(dtfit/(dtmicro*binwidth))
    if iend>len(tlist):
        iend = len(tlist) 
        
    # get background (by simply looking at last ten data points) and substract from intensity data.
    if ybg < 0:
        ybg = np.mean(ylist[-10:]) # mean background per histogram bin bin of length
            
    if method == 'ML': #maximum likelihood exponential fit
        [taufit,Afit] = MaxLikelihoodFit(tlist,ylist,istart,iend,ybg,False)
    if method == 'ML_c': #maximum likelihood exponential fit
        [taufit,Afit] = MaxLikelihoodFit_c(tlist,ylist,istart,iend,ybg,False)
    elif method == 'WLS': # weighted least squares fit
        [taufit,Afit] = WeightedLeastSquareFit(tlist,ylist,istart,iend,ybg,plotbool=False)
    else:
        taufit = 0; Afit = 0;
        print('Error: invalid fit method')
    
    if plotbool == True:
        plt.xlabel('time (ns)')
        plt.ylabel('')
        plt.semilogy(tlist,ylist,'.',tlist[istart:iend],Afit*np.exp(-(tlist[istart:iend]-tlist[istart])/taufit)+ybg)
        plt.semilogy([tlist[0],tlist[-1]],[ybg,ybg],'k--')
        plt.show()
        print('Fitted lifetime:',taufit,'ns; Amax:',Afit)

    # Amax is the maximum y-value
    Amax = np.max(ylist)
        
    return(taufit,Afit,ybg)  
    
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

#    if plotbool == True:
#        yest = np.array([Aest*np.exp(-(xdata[i]-xdata[0])/tauest)+bgcpb for i in range(len(xdata))])
#        plt.semilogy(tlist,ylist,'.',xdata,yest,[xdata[1],xdata[-1]],[bgcpb,bgcpb],'k--')
#        plt.show()        
    
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

def MaxLikelihoodFunction(params,xdata,ydata,const): 
    # max likelihood function for A*exp(-t/tau), needed in function MakLikelihoodFit
    # params = [A,tau]
    A = params[0]
    tau = params[1]   
    E = 0;
    for i in range(len(xdata)):
        E = E + ydata[i]*np.log(A*np.exp(-(xdata[i]-xdata[0])/tau)+const)-(A*np.exp(-(xdata[i]-xdata[0])/tau)+const)
        
    return(-E) # This function needs to be MINIMIZED (because of the minus sign) to have the maximum likelihood fit!
    
def WeightedLeastSquareFit(tlist,ylist,istart,iend,bgcpb,plotbool=False):

    # check if istart and iend are good numbers
    if istart<0 or istart>=len(ylist):
        istart = 0
        print('WARNING: adapted istart in MaxLikelihoodExpFit')
    if iend<=istart or iend>len(ylist):
        iend = len(ylist)
        print('WARNING: adapted iend in MaxLikelihoodExpFit')

    # shift t0 to t=0
    ydata = ylist[istart:iend]-bgcpb #subtract constant bg
    ydatalog = np.log(ydata) #take log (weighting)
    xdata = tlist[istart:iend]-tlist[istart] #shift to 0

    # do calculations    
    mod = lmfit.models.LinearModel()
    pars = mod.guess(ydatalog,x=xdata)
    out = mod.fit(ydatalog,pars,x=xdata)
    Aest = np.exp(out.best_values['intercept'])
    tauest = -1.0/out.best_values['slope']

    if plotbool == True:
        yest = np.array([Aest*np.exp(-(xdata[i])/tauest)+bgcpb for i in range(len(xdata))])
        plt.semilogy(tlist,ylist,'.',xdata+tlist[istart],yest,[xdata[0],xdata[-1]],[bgcpb,bgcpb],'k--')
        plt.show()
        print(out.fit_report())         
    
    return(tauest,Aest)
# In[1]: Parameters
    
#folder = 'C:/Users/rapha/OMEL/02 Projects/02 MD Plasmons/Adi Drexhage experiment/Optical Data/20190116_TTTRscannin_QD_test/'
folder = 'C:/Users/rober/Documents/Doktorat/Projects/Scanning_imaging/20190116_TTTRscannin_QD_test/'
name = 'scan3'

xbinning = 2
ybinning = 2

savebool = False

# In[2]: Analysis

# read data
[dtmicro, dtmacro, microtimes0, macrotimes0, microtimes1, macrotimes1, nrphotons0,nrphotons1,overflows,microtimesP,macrotimesP,nrP,microtimesL,macrotimesL,nrL] = importT3marked(folder + name + '.ptu') # import data from .ptu file
Texp = round(overflows*dtmacro*1024,0) # macrotime values can store up to 2^10 = 1024 before overflow. overflows*1024*dtmacro gives experiment time [s]

# get number of first photon per trigger pulse
limits0P = FindFirstPhotonAfterTrigger(macrotimes0,macrotimesP)
limits0L= FindFirstPhotonAfterTrigger(macrotimes0,macrotimesL)
 
# Plot decay histogram for all photons
plt.figure(1)
fitOut = GetLifetime(microtimes0,dtmicro,dtmacro,20e-9,plotbool=True,method='WLS',ybg=-1) # WLS = Weighted least squares fit

# Numba booster (c compilation of maximum-likelihood function)
MaxLikelihoodFunction_c = nb.jit(nopython=True)(MaxLikelihoodFunction) # implement c function
fitOut = GetLifetime(microtimes0,dtmicro,dtmacro,20e-9,plotbool=False,method='ML_c',ybg=-1) # compile c function for the first time

# Get intensities per pixel
ylen = int(nrL)
xlen = int(nrP/nrL)
Imat = np.zeros((xlen,ylen))
for l in range(ylen):
    for p in range(xlen):
        if p>0: # this is a special case because there is no pixel trigger at start of each line, but instead a line trigger
            Imat[p,l] = limits0P[p+xlen*l]-limits0P[p+xlen*l-1]
        else:
            Imat[p,l] = limits0P[p+xlen*l]-limits0L[l]

# Get lifetimes and intensities (binned)
taumat = np.zeros((int(xlen/xbinning),int(ylen/ybinning)))
Imat_bin = np.zeros((int(xlen/xbinning),int(ylen/ybinning)))
Iclustmat_bin = np.zeros((int(xlen/xbinning),int(ylen/ybinning)))
xcnt = 0;
ycnt = 0;

start_time = time.time()
for l in range(0,ylen,ybinning):
    print('line',l,'out of',ylen-1) #start counting at 0
    for p in range(0,xlen,xbinning):
        if p>0: # this is a special case because there is no pixel trigger at start of each line, but instead a line trigger
            microtimesp = microtimes0[limits0P[p+xlen*l-1]:limits0P[p+(xbinning-1)+xlen*l]]
            Ip = limits0P[p+xlen*l]-limits0P[p+xlen*l-1]
            for dl in range(1,ybinning-1):
                np.append(microtimesp,microtimes0[limits0P[p+xlen*(l+dl)-1]:limits0P[p+(xbinning-1)+xlen*(l+dl)]])
                Ip += limits0P[p+(xbinning-1)+xlen*(l+dl)] - limits0P[p+xlen*(l+dl)-1]
        else:
            microtimesp = microtimes0[limits0L[l]:limits0P[p+(xbinning-1)+xlen*l]]
            Ip = limits0P[p+(xbinning-1)+xlen*l]-limits0L[l]
            for dl in range(1,ybinning-1):
                np.append(microtimesp,microtimes0[limits0L[l+dl]:limits0P[p+(xbinning-1)+xlen*(l+dl)]])
                Ip += limits0P[p+(xbinning-1)+xlen*(l+dl)] - limits0L[l+dl]
        #plt.figure()
        fitOut = GetLifetime(microtimesp,dtmicro,dtmacro,20e-9,plotbool=False,method='ML_c') # ML = Maximum likelihood fit
        taumat[xcnt,ycnt] = fitOut[0]
        Imat_bin[xcnt,ycnt] = Ip
        Iclustmat_bin[xcnt,ycnt] = Ip>500
        xcnt += 1
    xcnt = 0 # reset column counter when reaching new line
    ycnt+=1
print("--- %s s for execution of lifetime fitting ---" % (time.time() - start_time))
# %% Timegated PL images
Imat_bing = np.zeros((int(xlen/xbinning),int(ylen/ybinning)))
xcnt = 0;
ycnt = 0;
measstart=10*1e-9/dtmicro #start is at 10 ns
timemin=5*1e-9/dtmicro
timemax=10e-9/dtmicro
myphotons = np.array([i for i in range(len(microtimes0)) if (measstart+timemin)<microtimes0[i]<(measstart+timemax)])
for l in range(0,ylen,ybinning):
    print('line',l,'out of',ylen-1) #start counting at 0
    for p in range(0,xlen,xbinning):
        if p>0: # this is a special case because there is no pixel trigger at start of each line, but instead a line trigger
            microtimesp = myphotons[limits0P[p+xlen*l-1]:limits0P[p+(xbinning-1)+xlen*l]]
            Ipg = limits0P[p+xlen*l]-limits0P[p+xlen*l-1]
            for dl in range(1,ybinning-1):
                np.append(microtimesp,myphotons[limits0P[p+xlen*(l+dl)-1]:limits0P[p+(xbinning-1)+xlen*(l+dl)]])
                Ipg += limits0P[p+(xbinning-1)+xlen*(l+dl)] - limits0P[p+xlen*(l+dl)-1]
        else:
            microtimesp = myphotons[limits0L[l]:limits0P[p+(xbinning-1)+xlen*l]]
            Ipg = limits0P[p+(xbinning-1)+xlen*l]-limits0L[l]
            for dl in range(1,ybinning-1):
                np.append(microtimesp,myphotons[limits0L[l+dl]:limits0P[p+(xbinning-1)+xlen*(l+dl)]])
                Ipg += limits0P[p+(xbinning-1)+xlen*(l+dl)] - limits0L[l+dl]
        #plt.figure()
        Imat_bing[xcnt,ycnt] = Ipg
        xcnt += 1
    xcnt = 0 # reset column counter when reaching new line
    ycnt+=1
# In[3]: Plotting
scancenter=[0,0]

scanrange=20     #in micrometer
points=160
xtargetvals=list(scancenter[0]+((j+0.5-points/2)/points)*scanrange for j in range(points))
ytargetvals=list(scancenter[1]+((j+0.5-points/2)/points)*scanrange for j in range(points))
# Plot intensity and lifetime maps
fig = plt.figure(1)
plt.subplot(1,3,1)
ax1 = fig.add_subplot(1,3,1)
plt.imshow(Imat_bin,extent=[min(xtargetvals),max(xtargetvals),min(ytargetvals),max(ytargetvals)])
plt.colorbar()
plt.title('PL intensity (counts per bin)')
plt.xlabel('x (um)')
plt.ylabel('y (um)')
    
ax2 = fig.add_subplot(1,3,2)
plt.imshow(taumat,vmin=3.5,vmax=7,extent=[min(xtargetvals),max(xtargetvals),min(ytargetvals),max(ytargetvals)])
plt.colorbar()
plt.show()
plt.title('Lifetime (ns)')
plt.xlabel('x (um)')
plt.ylabel('y (um)')

# Plot FLID map

ax3 = fig.add_subplot(1,3,3)
plt.hist2d(Imat_bin.reshape(Imat_bin.size,1)[:,0],taumat.reshape(taumat.size,1)[:,0],(50,50))
plt.title('FLID map')
plt.xlabel('Counts per bin')
plt.ylabel('Lifetime (ns)')
plt.colorbar()
#plt.axes().set_aspect('equal', 'datalim')

# plot triggers
fig2 = plt.figure(2)
plt.plot(macrotimesP,np.ones(macrotimesP.shape),'r.')
plt.plot(macrotimesL,np.ones(macrotimesL.shape),'b.')
plt.show()
testarray=np.array([macrotimesP[i+1]-macrotimesP[i] for i in range(len(macrotimesP)-1)])
# In[4]: Save arrays and plot to hard drive
if savebool == True:
    np.save(folder +'output/' + name +'_taumat',taumat)
    np.save(folder +'output/' + name +'_Imat',Imat)
    np.save(folder +'output/' + name + '_Imat_bin',Imat_bin)
    plt.savefig(folder + 'output/' + name + '.pdf', bbox_inches='tight',format='pdf')
    plt.savefig(folder + 'output/' + name, bbox_inches='tight')