import os, numpy as np, csv, matplotlib.pyplot as plt, scipy.optimize as opt, math, struct, binascii, gc, time, random

import multiprocessing
from operator import sub
from joblib import Parallel, delayed
import scipy, lmfit
import scipy.optimize
import scipy.signal
from scipy.optimize import minimize  # used for implementation of maximum likelihood exponential fit
from matplotlib import gridspec
from matplotlib import rc

from math import factorial
from scipy.stats import poisson
from enum import Enum
import pandas as pd

# Custom Style jsanten
import seaborn as sns
sns.set()
sns.color_palette("bright")
sns.set_color_codes(palette="bright")
#rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
#rc('text', usetex=True)


import matplotlib
matplotlib.use("Qt5Agg")


import pickle
import numba as nb
from tqdm import tqdm
import warnings


# I get errors unless the function that the cores execute in parallel is defined outside the class function
def hist2(x, y, bins):
    store = np.zeros(len(bins) - 1, dtype='float')
    for i in x:
        res = y.searchsorted(bins + i)
        store += res[1:] - res[:-1]
    return store


def load_obj(name, folder):
    with open(folder + name + '.pkl', 'rb') as f:
        return pickle.load(f)


# In[23]:


def ImportT3(filename, HHsettings):

    #Get the last three characters of the filename to determine if it is .out or .ptu
    file_extension = filename[len(filename) - 3:]

    if file_extension == "out":
        file_extension_is_out = True
        print(".out recognized")
    elif file_extension == "ptu":
        file_extension_is_out = False
        print(".ptu recognized")
    else:
        file_extension_is_out = False
        warnings.warn("file extension not set correctly")


    with open(filename, "rb+") as f:

        if not file_extension_is_out: #file is a .ptu file
            while True:
                if f.read(1) == b'M':
                    if f.read(18) == b'easDesc_Resolution':  # recognize time unit entry
                        break
            f.read(21)  # rest of the tag
            dtmicro = struct.unpack('d', f.read(8))[0]
            print('Microtime unit:', dtmicro)

            while True:
                if f.read(1) == b'M':
                    if f.read(24) == b'easDesc_GlobalResolution':  # recognize time unit entry
                        break
            f.read(15)  # rest of the tag
            dtmacro = struct.unpack('d', f.read(8))[0]
            print('Macrotime unit:', dtmacro)

            while True:
                if f.read(1) == b'T':
                    if f.read(23) == b'TResult_NumberOfRecords':  # recognize number of records entry
                        break
            f.read(16)  # rest of the tag
            nrrec = struct.unpack('q', f.read(8))[0]  # extract number of records
            # print('Number of records in file:', nrrec)

            while True:
                if f.read(1) == b'H':
                    if f.read(9) == b'eader_End':
                        # print('Header_End found')
                        break
            f.read(38)  # rest of Header_End

        elif file_extension_is_out:
            nrrec = HHsettings["overallCounts"]
            dtmicro = HHsettings["resolution"] * 1e-12  #in s
            dtmacro = 1 / HHsettings["syncRate"]    #in s
            print('Microtime unit:', dtmicro)
            print('Macrotime unit:', dtmacro)

        macrotimes0 = np.zeros(nrrec, dtype='int64')
        microtimes0 = np.zeros(nrrec, dtype='int64')
        macrotimes1 = np.zeros(nrrec, dtype='int64')
        microtimes1 = np.zeros(nrrec, dtype='int64')
        macrotimesfireA = np.zeros(nrrec, dtype='int64')
        microtimesfireA = np.zeros(nrrec, dtype='int64')
        macrotimesfireB = np.zeros(nrrec, dtype='int64')
        microtimesfireB = np.zeros(nrrec, dtype='int64')
        overflows = 0
        nrphotons0 = 0
        nrphotons1 = 0
        nrfireA = 0
        nrfireB = 0
        prevchann = 0

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
            elif channel == 1:
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
                macrotimesfireA[nrfireA] = macrotime + 1024 * overflows
                microtime = ((struct.unpack("I", entry)[0] >> 10) & 0x7FFF)
                microtimesfireA[nrfireA] = microtime
                nrfireA += 1
            elif channel == 72:
                macrotime = (struct.unpack("I", entry)[0] & 0x3FF)
                macrotimesfireB[nrfireB] = macrotime + 1024 * overflows
                microtime = ((struct.unpack("I", entry)[0] >> 10) & 0x7FFF)
                microtimesfireB[nrfireB] = microtime
                nrfireB += 1
            else:
                print('bad channel:', channel)

    microtimes0 = microtimes0[:nrphotons0]
    macrotimes0 = macrotimes0[:nrphotons0]
    microtimes1 = microtimes1[:nrphotons1]
    macrotimes1 = macrotimes1[:nrphotons1]
    microtimesfireA = microtimesfireA[:nrfireA]
    macrotimesfireA = macrotimesfireA[:nrfireA]
    microtimesfireB = microtimesfireB[:nrfireB]
    macrotimesfireB = macrotimesfireB[:nrfireB]

    print("")
    print('nrphotons0:', nrphotons0)
    print('nrphotons1:', nrphotons1)
    print('nrfireA:', nrfireA)
    print('nrfireB:', nrfireB)
    print('overflows:', overflows)

    return [dtmicro, dtmacro, microtimes0, macrotimes0, microtimes1, macrotimes1, nrphotons0, nrphotons1, overflows,
            microtimesfireA, macrotimesfireA, nrfireA, microtimesfireB, macrotimesfireB, nrfireB]


def importT3marked(filename, HHsettings):
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


    print('*** START READING FILE ***')

    with open(filename, "rb+") as f:

        nrrec = HHsettings["overallCounts"]
        dtmicro = HHsettings["resolution"] * 1e-12  # in s
        dtmacro = 1 / HHsettings["syncRate"]  # in s
        macrotimes0 = np.zeros(nrrec, dtype='int64')
        microtimes0 = np.zeros(nrrec, dtype='int64')
        macrotimes1 = np.zeros(nrrec, dtype='int64')
        microtimes1 = np.zeros(nrrec, dtype='int64')
        macrotimesfireA = np.zeros(nrrec, dtype='int64')
        microtimesfireA = np.zeros(nrrec, dtype='int64')
        macrotimesfireB = np.zeros(nrrec, dtype='int64')
        microtimesfireB = np.zeros(nrrec, dtype='int64')
        overflows = 0
        nrphotons0 = 0
        nrphotons1 = 0
        nrfireA = 0
        nrfireB = 0
        prevchann = 0

        # initialize arrays
        macrotimes0 = np.zeros(nrrec, dtype='int64')
        microtimes0 = np.zeros(nrrec, dtype='int64')
        macrotimes1 = np.zeros(nrrec, dtype='int64')
        microtimes1 = np.zeros(nrrec, dtype='int64')
        macrotimesP = np.zeros(nrrec, dtype='int64')
        microtimesP = np.zeros(nrrec, dtype='int64')
        macrotimesL = np.zeros(nrrec, dtype='int64')
        microtimesL = np.zeros(nrrec, dtype='int64')

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
            elif channel == 1:  # TODO: APD2: not implemented yet.
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
            elif channel == 72:  # channel M4 of the hydraharp (EMCCD)
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


def ShiftPulsedData(microtimes0, microtimes1, macrotimes0, macrotimes1, dtmicro, dtmacro):
    dtmax = 8

    [ylist1, xlist1] = np.histogram(microtimes1, int(dtmacro / dtmicro), [0, int(dtmacro / dtmicro)])
    [ylist0, xlist0] = np.histogram(microtimes0, int(dtmacro / dtmicro), [0, int(dtmacro / dtmicro)])
    tlist = (xlist0[:-1] + 0.5 * (xlist0[1] - xlist0[0])) * dtmicro * 1e9

    corrx = []
    corry = []  # find shift for which the two decay curves overlap most
    for i in range(-dtmax, dtmax):
        corrx.append(i)
        corry.append(sum(ylist1[dtmax:-dtmax] * ylist0[dtmax + i:-dtmax + i]))
    xmax = corry.index(max(corry))
    shift = corrx[xmax]

    tlist0 = (microtimes0 - shift) + macrotimes0 * int(dtmacro / dtmicro)
    tlist1 = microtimes1 + macrotimes1 * int(dtmacro / dtmicro)  # in units of dtmicro
    
    fig = plt.figure()
    plt.xlabel('time [ns]')
    plt.ylabel('counts [a.u.]')
    plt.semilogy(tlist, ylist0 + ylist1, label = "APD0 + APD1", color="k")
    plt.semilogy(tlist, ylist0, label = "APD0", color="r")
    plt.semilogy(tlist, ylist1, label = "APD1", color="b" )
    plt.legend()


    return (microtimes0 - shift, microtimes1, tlist0, tlist1, dtmicro, dtmacro, tlist, ylist0 + ylist1, fig)


def CropData(limits0, limits1, microtimes0, microtimes1, macrotimes0, macrotimes1, tstart, tend, binwidth):
    binstart = int(tstart / binwidth)
    binend = int(tend / binwidth)

    limits0crop = np.array(
        [limits0[i] for i in range(len(limits0)) if int(tstart / binwidth) <= i < int(tend / binwidth)])
    limits1crop = np.array(
        [limits1[i] for i in range(len(limits1)) if int(tstart / binwidth) <= i < int(tend / binwidth)])
    microtimes0crop = np.array(
        [microtimes0[i] for i in range(len(microtimes0)) if limits0crop[0] <= i <= limits0crop[-1]])
    microtimes1crop = np.array(
        [microtimes1[i] for i in range(len(microtimes1)) if limits1crop[0] <= i <= limits1crop[-1]])
    macrotimes0crop = np.array(
        [macrotimes0[i] for i in range(len(macrotimes0)) if limits0crop[0] <= i <= limits0crop[-1]])
    macrotimes1crop = np.array(
        [macrotimes1[i] for i in range(len(macrotimes1)) if limits1crop[0] <= i <= limits1crop[-1]])

    Texp = tend - tstart

    return (microtimes0crop, microtimes1crop, macrotimes0crop, macrotimes1crop, limits0crop, limits1crop, Texp)


def GetLifetime(microtimes, dtmicro, dtmacro, dtfit, tstart, binwidth=1,
                ybg=0):  # dtfit is the time interval considered for t
    # binwidth is a multiplier. actual binwidth is given as binwidth*dtmicro[s]
    [ylist, xlist] = np.histogram(microtimes, bins = int(dtmacro / (dtmicro * binwidth)), range = [0, int(dtmacro / dtmicro)])
    tlist = (xlist[:-1] + 0.5 * (xlist[1] - xlist[0])) * dtmicro * 1e9

    istart = int(tstart / dtmicro)  # find index of maximum element in ylist #TODO: Why does this work?
    if istart < 0:
        istart = ylist.argmax()
    iend = istart + int(dtfit / (dtmicro * binwidth))
    if iend > len(tlist):
        iend = len(tlist)

        # get background (by simply looking at last ten data points) and substract from intensity data.
    if ybg == 0:
        ybg = np.mean(ylist[-10:-1])  # mean background per histogram bin of length

    # exponential decay
    # expfit = np.polyfit(tlist[istart:iend], np.log(ylist[istart:iend]), 1, w=np.sqrt(ylist[istart:iend])) #weighted exponential fit
    # tauexp = -1/expfit[0]

    # new exponential fit
    #[expfit,_] = scipy.optimize.curve_fit(lambda t,a,b: a*np.exp(b*t)+ybg,  tlist[istart:iend],  ylist[istart:iend],  p0=(1e3, -1/40))
    #tauexp = -1/expfit[1]
    #Aexp = expfit[0]

    # maximum likelihood exponential fit
    [tauexp, Aexp] = MaxLikelihoodFit(tlist, ylist, istart, iend, ybg, False)
    print("Aexp", Aexp, " , tauexp:", tauexp)
    # average lifetime
    tauave = np.sum(
        (tlist[istart:iend] - tlist[istart]) * (ylist[istart:iend] - ybg) / np.sum(ylist[istart:iend] - ybg))

    plt.xlabel('time (ns)')
    plt.ylabel('')
    plt.semilogy(tlist,ylist,'c.', alpha = 0.1)
    #plt.semilogy(tlist[istart:iend], np.exp(expfit[1])*np.exp(expfit[0]*tlist[istart:iend]), color = "y")
    #plt.plot(tlist, ylist / np.mean(ylist), 'r.')
    plt.plot(tlist[istart:iend], Aexp * np.exp(-(tlist[istart:iend] - tlist[istart]) / tauexp) + ybg, "r",
             label = "exponential fit: Axp = {:.2f}, tauexp = {:.2f}".format(Aexp, tauexp))
    plt.semilogy([tlist[0],tlist[-1]],[ybg,ybg],'k--', label = "background level = {:.2f}".format(ybg))
    plt.legend()
    #plt.show()

    # Amax is the maximum y-value
    Amax = np.max(ylist)

    print('single-exp lifetime:', tauexp, 'ns; average lifetime:', tauave, 'ns; Amax:', Amax)

    return (tauave, Amax, ybg)


def HistPhotons(photontlist, binwidth, Texp):  # finds the index of the first photon in each bin. photontlist is in [s]
    histmax = Texp  # experiment duration [s]

    nrbins = int(Texp / binwidth)
    limits = np.full(nrbins, len(photontlist))
    counter, i = 0, 0
    while counter < nrbins and i < len(photontlist):
        while photontlist[i] > counter * binwidth:
            limits[counter] = i
            counter += 1
        i += 1

    return (limits)


def MakeIntTrace(limits0, limits1, binwidth, Texp):
    nrbins = int(Texp / binwidth)
    inttrace = np.array(
        [(limits0[binnr + 1] - limits0[binnr] + limits1[binnr + 1] - limits1[binnr]) for binnr in range(nrbins - 1)])

    fig = plt.figure(figsize=(15*0.73, 5*0.73))#plt.figure(figsize=(15, 5))
    gs = gridspec.GridSpec(1, 2, width_ratios=[4, 1])
    ax0 = plt.subplot(gs[0])
    xarray = np.arange(len(inttrace)) * binwidth
    yarray = inttrace
    p0 = ax0.plot(xarray, yarray, '-', linewidth=0.5, color = "k")
    plt.xlabel('time (s)', fontsize="large")
    plt.ylabel('counts / %i ms' % (binwidth * 1e3),fontsize="large")
    plt.xlim([0, Texp])
    plt.ylim([0, 1.1 * np.max(inttrace)])

    histogram = np.histogram(inttrace, max(inttrace), [0, max(inttrace)])

    ax1 = plt.subplot(gs[1])
    ax1.plot(histogram[0], 0.5 * (histogram[1][:-1] + histogram[1][1:]), color="k")
    plt.xlabel('occurrence',fontsize="large")
    plt.ylabel('counts / %i ms' % (binwidth * 1e3),fontsize="large")
    plt.ylim([0, 1.1 * np.max(inttrace)])
    fig.tight_layout()

    return (inttrace, fig, xarray, yarray)


def MakeTauTrace(taubinlist, intbinlist, binwidth, Texp, taumin=0, taumax=100, intmin=0, intmax=100, col='k'):
    nrbins = int(Texp / binwidth)

    fig = plt.figure(figsize=(15, 7))
    gs = gridspec.GridSpec(2, 2, width_ratios=[4, 1])
    ax0 = plt.subplot(gs[0])
    p0 = ax0.plot(np.arange(len(intbinlist)) * binwidth, intbinlist, '-', linewidth=0.5, color=col)
    plt.xlabel('time (s)')
    plt.ylabel('counts / %i ms' % (binwidth * 1e3))
    plt.xlim([0, Texp])
    plt.ylim([intmin, intmax])

    histogram = np.histogram(intbinlist, int(np.max(intbinlist)), [0, int(np.max(intbinlist))])

    ax1 = plt.subplot(gs[1])
    ax1.plot(histogram[0], 0.5 * (histogram[1][:-1] + histogram[1][1:]), color=col)
    plt.xlabel('occurrence')
    plt.ylabel('counts / %i ms' % (binwidth * 1e3))
    plt.ylim([intmin, intmax])

    ax2 = plt.subplot(gs[2])
    p2 = ax2.plot(np.arange(len(intbinlist)) * binwidth, taubinlist, '.', markersize=1.5, color=col)
    plt.xlabel('time (s)')
    plt.ylabel('lifetime (ns)')
    plt.xlim([0, Texp])
    plt.ylim([taumin, taumax])

    histogram = np.histogram(taubinlist, taumax - taumin, [taumin, taumax])

    ax3 = plt.subplot(gs[3])
    ax3.plot(histogram[0], 0.5 * (histogram[1][:-1] + histogram[1][1:]), color=col)
    plt.xlabel('occurrence')
    plt.ylabel('lifetime (ns)')
    plt.ylim([taumin, taumax])

    plt.hold(True)


def BinIntensity(microtimes0, times0, limits0, microtimes1, times1, limits1, dtmicro, dtmacro, onintlim, offintlim):
    ## select only data with high or low intensity

    plt.title('total decay')
    tauave = GetLifetime(np.append(microtimes0, microtimes1), dtmicro, dtmacro, 200e-9, -1)

    nrbins = len(limits0)
    inttrace = np.array(
        [limits0[binnr + 1] - limits0[binnr] + limits1[binnr + 1] - limits1[binnr] for binnr in range(nrbins - 1)])

    # find photons in on period
    onphotonlist0 = np.array(
        [np.arange(limits0[binnr], limits0[binnr + 1]) for binnr in range(nrbins - 1) if inttrace[binnr] >= onintlim])
    onphotonlist0 = np.concatenate(onphotonlist0).ravel()
    onphotonlist1 = np.array(
        [np.arange(limits1[binnr], limits1[binnr + 1]) for binnr in range(nrbins - 1) if inttrace[binnr] >= onintlim])
    onphotonlist1 = np.concatenate(onphotonlist1).ravel()

    onmicrotimes0 = np.array([microtimes0[i] for i in onphotonlist0])
    onmicrotimes1 = np.array([microtimes1[i] for i in onphotonlist1])
    ontimes0 = np.array([times0[i] for i in onphotonlist0])
    ontimes1 = np.array([times1[i] for i in onphotonlist1])
    plt.title('on decay')
    ontauave = GetLifetime(np.append(onmicrotimes0, onmicrotimes1), dtmicro, dtmacro, 200e-9, -1)

    # find photons in off period
    offphotonlist0 = np.array(
        [np.arange(limits0[binnr], limits0[binnr + 1]) for binnr in range(nrbins - 1) if inttrace[binnr] < offintlim])
    offphotonlist0 = np.concatenate(offphotonlist0).ravel()
    offphotonlist1 = np.array(
        [np.arange(limits1[binnr], limits1[binnr + 1]) for binnr in range(nrbins - 1) if inttrace[binnr] < offintlim])
    offphotonlist1 = np.concatenate(offphotonlist1).ravel()

    offmicrotimes0 = np.array([microtimes0[i] for i in offphotonlist0])
    offmicrotimes1 = np.array([microtimes1[i] for i in offphotonlist1])
    offtimes0 = np.array([times0[i] for i in offphotonlist0])
    offtimes1 = np.array([times1[i] for i in offphotonlist1])
    plt.title('off decay')
    offtauave = GetLifetime(np.append(offmicrotimes0, offmicrotimes1), dtmicro, dtmacro, 10e-9, -1)

    return (onmicrotimes0, offmicrotimes0, ontimes0, offtimes0, onmicrotimes1, offmicrotimes1, ontimes1, offtimes1)


def SliceHistogram(microtimes0, times0, limits0, microtimes1, times1, limits1, dtmicro, dtmacro, Imin, Imax):
    ## select only data with intensity between Imin and Imax

    nrbins = len(limits0)
    inttrace = np.array(
        [limits0[binnr + 1] - limits0[binnr] + limits1[binnr + 1] - limits1[binnr] for binnr in range(nrbins - 1)])

    # find photons in bins with intensities in (Imin,Imax] range
    onphotonlist0 = np.array([np.arange(limits0[binnr], limits0[binnr + 1]) for binnr in range(nrbins - 1) if
                              Imin < inttrace[binnr] <= Imax], dtype=object)
    onphotonlist0 = np.concatenate(onphotonlist0).ravel()
    onphotonlist1 = np.array([np.arange(limits1[binnr], limits1[binnr + 1]) for binnr in range(nrbins - 1) if
                              Imin < inttrace[binnr] <= Imax], dtype=object)
    onphotonlist1 = np.concatenate(onphotonlist1).ravel()

    onmicrotimes0 = np.array([microtimes0[i] for i in onphotonlist0])
    onmicrotimes1 = np.array([microtimes1[i] for i in onphotonlist1])
    ontimes0 = np.array([times0[i] for i in onphotonlist0])
    ontimes1 = np.array([times1[i] for i in onphotonlist1])

    # count nr of time bins with intensities corresponding to slice intensity
    onbincount = 0
    for binnr in range(nrbins - 1):
        if Imin < inttrace[binnr] <= Imax:
            onbincount += 1

    return (onmicrotimes0, ontimes0, onmicrotimes1, ontimes1, onbincount)


def GetVoltageLimits(times,
                     macrotimesfireA):  # finds the index of the first photon in each bin. times in [microtimeunits]. macrotimes in [macrotimeunits]

    nrVbins = len(macrotimesfireA)
    Vlimits = np.full(nrVbins, len(times))
    counter, i = 0, 0
    while i < len(times):
        while (counter < nrVbins - 1) & (times[i] * dtmicro / dtmacro > macrotimesfireA[counter]):
            Vlimits[counter] = i
            counter += 1
        i += 1

    return (Vlimits)


def BinVoltage(microtimes0, times0, microtimes1, times1, macrotimesfireA):
    nrVbins = len(macrotimesfireA)

    limitsfire0 = GetVoltageLimits(times0, macrotimesfireA)
    limitsfire1 = GetVoltageLimits(times1, macrotimesfireA)

    # find photons in V=on and V=off periods
    Vonphotonlist0 = np.array([np.arange(limitsfire0[binnr], limitsfire0[binnr + 1]) for binnr in range(nrVbins - 1) if
                               np.mod(binnr, 2) == 0])
    Voffphotonlist0 = np.array([np.arange(limitsfire0[binnr], limitsfire0[binnr + 1]) for binnr in range(nrVbins - 1) if
                                np.mod(binnr, 2) != 0])
    Vonphotonlist1 = np.array([np.arange(limitsfire1[binnr], limitsfire1[binnr + 1]) for binnr in range(nrVbins - 1) if
                               np.mod(binnr, 2) == 0])
    Voffphotonlist1 = np.array([np.arange(limitsfire1[binnr], limitsfire1[binnr + 1]) for binnr in range(nrVbins - 1) if
                                np.mod(binnr, 2) != 0])

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
    # plt.title('V on')
    # Vontauave = GetLifetime(np.append(Vonmicrotimes0,Vonmicrotimes1),dtmicro,dtmacro,200e-9,-1)

    # get V=off microtimest
    Voffmicrotimes0 = np.array([microtimes0[i] for i in Voffphotonlist0])
    Voffmicrotimes1 = np.array([microtimes1[i] for i in Voffphotonlist1])
    Vofftimes0 = np.array([times0[i] for i in Voffphotonlist0])
    Vofftimes1 = np.array([times1[i] for i in Voffphotonlist1])
    # plt.title('V off')
    # Vofftauave = GetLifetime(np.append(Voffmicrotimes0,Voffmicrotimes1),dtmicro,dtmacro,200e-9,-1)

    return (
    Vonmicrotimes0, Voffmicrotimes0, Vontimes0, Vofftimes0, Vonmicrotimes1, Voffmicrotimes1, Vontimes1, Vofftimes1)




def MakeG2(times0, times1, dtmicro, dtmacro, microtimes0, microtimes1, g2restime=8e-9, nrbins=200, normalized=True):
    i0 = 0
    i1 = 0
    lim1 = 0
    g2 = np.zeros(2 * nrbins)
    g2B = np.zeros(2*nrbins)
    g2C = np.zeros(2*nrbins)
    blindB = 2e-9
    blindC = 5e-9

    g2res = g2restime / dtmicro  # transform g2restime [s] to g2res [microtime units]
    blindB = blindB / dtmicro
    blindC = blindC / dtmicro

    # correlate det0 with det1 (positive time differences)
    for i0 in range(len(times0)):
        t0 = times0[i0]
        i1 = 0
        while True:
            if lim1 + i1 < len(times1):  # check if we've already reached end of photon stream on det1
                dt = times1[lim1 + i1] - t0
                if dt < 0:  # find index lim1 of first photon on det1 that came after photon i0 on det0
                    lim1 = lim1 + 1
                else:
                    binnr = int(dt / g2res)  # calculate binnr that corresponds to dt
                    if binnr < nrbins:  # check if time difference is already large enough to stop correlation
                        g2[nrbins + binnr] += 1  # increase counter in corresponding bin by one
                        if microtimes0[i0] > blindB and microtimes1[lim1+i1] > blindB:
                            g2B[nrbins + binnr] += 1
                        if microtimes0[i0] > blindC and microtimes1[lim1+i1] > blindC:
                            g2C[nrbins + binnr] += 1
                        i1 = i1 + 1  # look at next photon on det1
                    else:
                        break  # dt larger than maximum correlation width. stop.
            else:
                break  # end of photon stream on det1 reached. stop.

    # correlate det1 with det0 (positive time differences)
    lim1 = 0
    for i0 in range(len(times1)):
        t0 = times1[i0]
        i1 = 0
        while True:
            if lim1 + i1 < len(times0):
                dt = times0[lim1 + i1] - t0
                if dt < 0:
                    lim1 = lim1 + 1
                else:
                    binnr = int(dt / g2res)
                    if binnr < nrbins:
                        g2[nrbins - 1 - binnr] += 1
                        if microtimes0[lim1+i1] > blindB and microtimes1[i0] > blindB:
                           g2B[nrbins - 1 - binnr] += 1
                        if microtimes0[lim1+i1] > blindC and microtimes1[i0] > blindC:
                           g2C[nrbins - 1 - binnr] += 1
                        i1 = i1 + 1
                    else:
                        break
            else:
                break

    g2tlist = np.arange(-g2res * dtmicro * (nrbins - 0.5), g2res * dtmicro * nrbins, g2restime) * 1e9

    if normalized:
        if np.max(g2) < 30: # if the emission is too low, do not fit the plot
            state = "discarded"
            title = r"non-normalized $g^{2}(\tau)$ correlation" + " | max peak height: {} counts".format(
                np.round(np.max(g2), 2)) + " | " + state
            print("np.max(g2): ", np.max(g2))
            g2_tau_zero = np.NaN
        else:  # emission is high enough
            # Calculate time difference based on binning
            time_diff = g2restime * 1e9  # [nano seconds]

            # prepare finding the center peak
            center_delta = 35  # [nano seconds] Assuming the center value is in an interval = [-center_delta,+center_delta]
            center_delta_idx = int(
                center_delta / time_diff)  # check how many values to include in the search for the center peak

            # Find Peaks that are at least 50% of the maximum value
            g2_no_center = g2.copy()  # temp copy of g2 that has the zero peak removed
            g2_no_center[nrbins - center_delta_idx: nrbins + center_delta_idx] = 0

            peaks, _ = scipy.signal.find_peaks(g2_no_center, height=0.5 * np.max(g2), prominence=0.8, distance=5)
            yvalues_peaks = g2[peaks]
            tvalues_peaks = g2tlist[peaks]

            # Normalize by mean of peak height
            mean_peaks = np.mean(yvalues_peaks)
            g2 = g2 / mean_peaks

            g2_center = g2[nrbins - center_delta_idx: nrbins + center_delta_idx]  # yvalues in center
            g2tlist_center = g2tlist[nrbins - center_delta_idx: nrbins + center_delta_idx]  # tvalues in center

            print("mean_peaks: ", mean_peaks)
            # find peak(s) in center
            peaks_center, _ = scipy.signal.find_peaks(g2_center, height=0.04 * np.max(g2), prominence=0.4)

            # check that the peaks array is non-empty. IF empty: decrease the prominence iteratively until
            # no longer empty
            prominence_decreasing = 0.4
            while np.size(peaks_center) == 0:
                prominence_decreasing *= 0.9
                peaks_center, _ = scipy.signal.find_peaks(g2_center, height=0.04 * np.max(g2), prominence=prominence_decreasing)
                warnings.warn("Adapted prominence for center peak. Will be removed in future versions")
            center_peak = g2_center[peaks_center][0]  # define a singular center peak


            # Determine the Boundaries for the Integral (i.e., each peak should be in one interval)
            array_containing_peaks = [] # this array contains the cut-off values for the peaks
            for i in range(-int(np.size(yvalues_peaks) / 2), int(np.size(yvalues_peaks) / 2) + 2):
                i -= 1/2
                cut_value = center_peak + i * dtmacro * 1e9

                if not ((cut_value > np.max(g2tlist)) or (cut_value < np.min(g2tlist))):
                    array_containing_peaks.append(cut_value)
                    #plt.axvline(x=cut_value, color="m", alpha=0.05)

            #print(array_containing_peaks)
            integral = []
            integral_idx = []
            #print("array_containing_peaks: ", array_containing_peaks)
            for i in range(0,np.size(array_containing_peaks) - 1):
                array_indices = [ind for ind, x in enumerate(g2tlist) if array_containing_peaks[i] <= x <= array_containing_peaks[i+1]]
                #print("array_indices: ", array_indices)
                if np.size(array_indices) == 0:
                    warnings.warn("Empty array. This might indicate that the wrong value is used for tau = 0")
                    continue
                integral.append(np.sum(g2[np.min(array_indices): np.max(array_indices)+1]))
                integral_idx.append((np.min(array_indices), np.max(array_indices)+1))

            #print("integral_idx: ", integral_idx)
            # normalize by outer outer_peaks_no outer peaks (outer_peaks_no on the negative / positive side each)
            outer_peaks_no = 5  # TODO: implement a safeguard, s.t. outer peaks does not include the center peak
            g2_tau_zero = integral[int(np.floor(np.size(integral) / 2))] / np.mean(
                np.append(np.array(integral[-outer_peaks_no:]), integral[:outer_peaks_no]))

            for i in range(0, outer_peaks_no):  # start of array
                plt.axvline(x=g2tlist[integral_idx[i][0]], linestyle=":", color="b", alpha=1)
                plt.axvline(x=g2tlist[integral_idx[i][1]], linestyle=":", color="b", alpha=1)

            for i in range(1, outer_peaks_no + 1):  # end of array
                plt.axvline(x=g2tlist[integral_idx[-i][0]], linestyle=":", color="b", alpha=1)
                plt.axvline(x=g2tlist[integral_idx[-i][1]], linestyle=":", color="b", alpha=1)
            print("g2_tau_zero: ", g2_tau_zero)

            # This is also an intersting plot that allows you to observe the envelope
            #plt.plot(integral)

            # Determine the state {singular, non-non_singular, discarded} of the dots
            # FIXME: This should be implemented using an enum / some type of class
            if mean_peaks < 30:
                state = "discarded"
                g2_tau_zero = np.NaN
            elif 0 <= g2_tau_zero < 0.5:
                state = "singular"
            elif 0.5 <= g2_tau_zero <= 1.8:
                state = "non_singular"
            else:
                state = "discarded"

            #plt.axhline(y=1, linestyle=":", alpha=0.5, color="k")
            plt.plot(g2tlist[peaks], g2[peaks], "x", color="r") # mark the peaks used for the fit
            plt.plot(g2tlist_center[peaks_center], g2_center[peaks_center], "*", color="b",
                     label=r"$g^{2}(\tau = 0) =$" + str(np.round(g2_center[peaks_center], 2)) + "(" + str(np.round(g2_tau_zero,2)) + ")")# mark the center peak
            title = r"normalized $g^{2}(\tau)$ correlation " + "| mean peak height: {}".format(
                np.round(mean_peaks, 2)) + r" $\pm$ " + str(np.round(np.std(yvalues_peaks),2)) + " counts | " + state
            plt.legend(loc=1)
            # TODO: Implement peak finding threshold based on known distance between peaks in x

    else:
        title = r"non-normalized $G^{2}(\tau)$ correlation"

    plt.plot(g2tlist, g2, color="k")
    #plt.plot(g2tlist,g2B) #TODO: What are these plots?
    #plt.plot(g2tlist,g2C)
    plt.title(title, fontsize="large")
    plt.xlabel('delay [ns]', fontsize="large")
    plt.ylabel('occurrence [a.u.]', fontsize="large")
    plt.ylim([0, 1.1 * max(g2)])
    #plt.show()

    return (g2tlist, g2, g2restime, nrbins, state, g2_tau_zero)


def MaxLikelihoodFit(tlist, ylist, istart, iend, bgcpb, plotbool=False):
    ### Maximum likelihood routine to fit single exponential. Pro: Works also for small amount of data (single bins of 10ms!)
    # tlist: x-axis values, here time in ns; ylist: y-axis values, here cts per tlist-bin; istart and iend: first and last element of tlist and ylist that are considered for the fit.

    # check if istart and iend are good numbers
    if istart < 0 or istart >= len(ylist):
        istart = 0
        print('WARNING: adapted istart in MaxLikelihoodExpFit')
    if iend <= istart or iend > len(ylist):
        iend = len(ylist)
        # FIXME: For g2_AP-3-120_SingleDot_2x10+7Hz_laserexct_25-ND03-ND20_dot9.ptu this gets printed every time
        #print('WARNING: adapted iend in MaxLikelihoodExpFit')

    # shift t0 to t=0
    ydata = ylist[istart:iend]
    xdata = tlist[istart:iend]

    # do calculations
    initParams = [np.max(ydata), 25]  # initial guess for A and tau
    results = minimize(MaxLikelihoodFunction, initParams, args=(xdata, ydata, bgcpb),
                       method='Nelder-Mead')  # minimize the negative of the maxlikelihood function instead of maximimizing
    Aest = results.x[0]  # get results of fit, A
    tauest = results.x[1]  # get results of fit, tau

    if plotbool == True:
        yest = np.array([Aest * np.exp(-(xdata[i] - xdata[0]) / tauest) + bgcpb for i in range(len(xdata))])
        plt.semilogy(tlist, ylist, '.', xdata, yest, [xdata[1], xdata[-1]], [bgcpb, bgcpb], 'k--')
        #plt.show()

    return (tauest, Aest)


def MaxLikelihoodFunction(params, xdata, ydata, const):
    # max likelihood function for A*exp(-t/tau), needed in function MakLikelihoodFit
    # params = [A,tau]
    A = params[0]
    tau = params[1]
    E = 0
    for i in range(len(xdata)):
        E = E + ydata[i] * np.log(A * np.exp(-(xdata[i] - xdata[0]) / tau) + const) - (
                    A * np.exp(-(xdata[i] - xdata[0]) / tau) + const)

    return (-E)  # This function needs to be MINIMIZED (because of the minus sign) to have the maximum likelihood fit!


def MaxLikelihoodFit_c(tlist, ylist, istart, iend, bgcpb, MaxLikelihoodFunction_c, plotbool=False):
    ### Maximum likelihood routine to fit single exponential. Pro: Works also for small amount of data (single bins of 10ms!)
    # tlist: x-axis values, here time in ns; ylist: y-axis values, here cts per tlist-bin; istart and iend: first and last element of tlist and ylist that are considered for the fit.

    # check if istart and iend are good numbers
    if istart < 0 or istart >= len(ylist):
        istart = 0
        print('WARNING: adapted istart in MaxLikelihoodExpFit')
    if iend <= istart or iend > len(ylist):
        iend = len(ylist)
        # FIXME: For g2_AP-3-120_SingleDot_2x10+7Hz_laserexct_25-ND03-ND20_dot9.ptu this gets printed every time
        #print('WARNING: adapted iend in MaxLikelihoodExpFit')

    # shift t0 to t=0
    ydata = ylist[istart:iend]
    xdata = tlist[istart:iend]

    # do calculations
    initParams = [np.max(ydata), 25]  # initial guess for A and tau
    results = minimize(MaxLikelihoodFunction_c, initParams, args=(xdata, ydata, bgcpb),
                       method='Nelder-Mead')  # minimize the negative of the maxlikelihood function instead of maximimizing
    Aest = results.x[0]  # get results of fit, A
    tauest = results.x[1]  # get results of fit, tau

    #    if plotbool == True:
    #        yest = np.array([Aest*np.exp(-(xdata[i]-xdata[0])/tauest)+bgcpb for i in range(len(xdata))])
    #        plt.semilogy(tlist,ylist,'.',xdata,yest,[xdata[1],xdata[-1]],[bgcpb,bgcpb],'k--')
    #        plt.show()

    if plotbool == True:
        yest = np.array([Aest * np.exp(-(xdata[i] - xdata[0]) / tauest) + bgcpb for i in range(len(xdata))])
        plt.plot(tlist, ylist, 'r.', xdata, yest, [xdata[1], xdata[-1]], [bgcpb, bgcpb], 'k--')
        plt.xlim([xdata[1], xdata[-1]])
        #plt.show()

    return (tauest, Aest)


def CorrelateSlices(slicewidth, tlist0, tlist1, limits0, limits1, intlist, binwidth, tmicro, tmacro):
    base = 1.2

    intlimits = np.arange(slicewidth * int(min(intlist) / slicewidth), slicewidth * int(max(intlist) / slicewidth),
                          slicewidth)
    slices = np.array([np.array([i for i, x in enumerate(intlist) if intlimits[j] < x <= intlimits[j + 1]]) for j in
                       range(len(intlimits) - 1)])

    correlations = []
    correlations2 = []
    norm1list = []
    norm2list = []
    norm3list = []
    norm4list = []
    tbin = binwidth / tmicro
    tbinp = int(tbin)  # bin width in units of the fundamental microtime unit
    Tbin = binwidth / tmacro
    Tbinp = int(Tbin)  # bin width in units of the fundamental macrotime unit
    resratio = int(tmacro / tmicro)

    # logbins = np.append(np.insert(base**np.arange(0,int(np.log(tmax)/np.log(base)+1),dtype='int64'),0,0),tmax).astype(np.int64)
    logbins = np.insert(int(0.5 * resratio) + (resratio * np.append(
        np.insert(base ** np.arange(0, int(np.log(Tbin) / np.log(base) + 1), dtype='int64'), 0, 0), Tbin).astype(
        np.int64)), 0, 0)
    logbins = np.unique(logbins)
    logbins = np.append(logbins,
                        [int(x) for x in np.divide([2e-7, 6e-7, 2e-6, 6e-6, 2e-5, 6e-5, 2e-4, 6e-4, 2e-3], tmicro)])

    for slicenr in (np.arange(len(slices), 0, -1) - 1):  # count down from the highest intensity
        lower0 = np.array([limits0[x] for x in
                           slices[slicenr]])  # list of first photons on detector 0 in the bins belonging to the slice
        upper0 = np.array([limits0[x + 1] for x in slices[slicenr]])
        lower1 = np.array([limits1[x] for x in slices[slicenr]])
        upper1 = np.array([limits1[x + 1] for x in slices[slicenr]])

        total0 = 0
        total1 = 0
        histavg = np.zeros(len(logbins) - 1, dtype='float')
        histavg2 = np.zeros(len(logbins) - 1, dtype='float')
        calctime0 = time.time()

        for i in range(len(lower0)):
            selection0 = tlist0[lower0[i]:upper0[i]]
            selection1 = tlist1[lower1[i]:upper1[i]]
            total0 += upper0[i] - lower0[i]  # len(selection0)
            total1 += upper1[i] - lower1[i]  # len(selection1)

            if len(selection0) == 0 or len(selection1) == 0:
                continue

            hist = np.zeros(len(logbins) - 1, dtype='float')
            for counterA in range(len(selection1)):
                res = selection0.searchsorted(logbins + selection1[counterA])
                hist += (res[1:] - res[:-1])

            # hist = np.zeros(len(logbins)-1,dtype='float');
            for counterB in range(len(selection0)):
                res = selection1.searchsorted(logbins + selection0[counterB])
                hist += (res[1:] - res[:-1])

            histavg += hist / (len(selection1) + len(selection0)) * (2 * Tbinp)
            histavg2 += hist

        if slicenr % 5 == 0:
            print('Analysis time for slice #', slicenr, ':', (time.time() - calctime0), 's')

        norm = 0.5 * (logbins[:-1] - logbins[1:]) * (logbins[:-1] + logbins[1:] - 2 * tbin - 1)
        tlist = 0.5 * (logbins[:-1] + logbins[1:]) * data[0]
        norm1 = total1
        norm2 = len(slices[slicenr])
        norm3 = total0
        norm4 = tbinp

        correlations.append(np.divide(histavg, norm / norm[1]))
        correlations2.append(np.divide(histavg2 / total1, norm / norm[1]))
        norm1list.append(norm1)
        norm2list.append(norm2)
        norm3list.append(norm3)
        norm4list.append(norm4)

    return (tlist, correlations, correlations2, norm1list, norm2list, norm3list, norm4list)


def DecaySlices(slicewidth, micro0, micro1, limits0, limits1, inttrace, dt, tmicro, tmacro):
    intlimits = np.arange(slicewidth * int(min(inttrace) / slicewidth), slicewidth * int(max(inttrace) / slicewidth),
                          slicewidth)
    slices = np.array([np.array([i for i, x in enumerate(inttrace) if intlimits[j] < x <= intlimits[j + 1]]) for j in
                       range(len(intlimits) - 1)])

    decaycurves = np.zeros((len(slices), int(tmacro / dt) + 1))

    for slicenr in (np.arange(len(slices), 0, -1) - 1):  # count down from the highest intensity
        lower0 = np.array([limits0[x] for x in
                           slices[slicenr]])  # list of first photons on detector 0 in the bins belonging to the slice
        upper0 = np.array([limits0[x + 1] for x in slices[slicenr]])
        lower1 = np.array([limits1[x] for x in slices[slicenr]])
        upper1 = np.array([limits1[x + 1] for x in slices[slicenr]])

        lower0 = [int(x) for x in lower0]
        lower1 = [int(x) for x in lower1]
        upper0 = [int(x) for x in upper0]
        upper1 = [int(x) for x in upper1]

        ylist = np.zeros(int(tmacro / dt))
        for j in range(len(lower0)):
            ylist += (
            np.histogram(micro0[lower0[j]:upper0[j]], int(tmacro / dt), [0, int(tmacro / dt) * int(dt / tmicro)])[0])
            ylist += (
            np.histogram(micro1[lower1[j]:upper1[j]], int(tmacro / dt), [0, int(tmacro / dt) * int(dt / tmicro)])[0])

        ylist = np.append(ylist, len(slices[slicenr]))

        decaycurves[slicenr] = ylist

    return (decaycurves)


def doublepoisson(k, a, b, c, d):
    return a * poisson.pmf(k.astype(int), b) + c * poisson.pmf(k.astype(int), d)


def singlepoisson(k, a, b):
    return a * poisson.pmf(k.astype(int), b)


def singlepoissonfit(xdata, ydata, params0, plotbool=False):
    # do calculations
    [popt, _] = scipy.optimize.curve_fit(singlepoisson, xdata, ydata, p0=params0)
    if plotbool == True:
        plt.plot(xdata, ydata, 'k.')
        plt.plot(xdata, singlepoisson(xdata, *popt))
    return popt


def doublepoissonfit(xdata, ydata, params0, plotbool=False):
    # do calculations
    [popt, _] = scipy.optimize.curve_fit(doublepoisson, xdata, ydata, p0=params0)
    if plotbool == True:
        plt.plot(xdata, ydata, 'k.')
        plt.plot(xdata, doublepoisson(xdata, *popt))
    return popt


def readspec(filename, dy=0):  # read andor's spectral data (ascii)
    # dy gives width in y-direction for spectral averaging

    data = np.genfromtxt(filename, delimiter=',')
    data[:, 0]

    # find length of spectra
    i = 0
    while data[i, 0] < data[i + 1, 0]: i += 1  # find end of first spectrum by checking when maximum lambda is reached
    xlen = i + 1
    ylen = data.shape[1] - 1
    nrimgs = int(len(data) / xlen)

    # sum y
    ysums = np.array([np.sum(data[:, i]) for i in range(1, ylen)])
    ymaxind = np.argmax(ysums) + 1  # +1 because zero are lambda values!
    yminind = 1  # BACKGROUND IS HARD-COPIED HERE!!!!

    # get lambda
    lambd = data[0:xlen, 0]

    # get spec
    spec = np.zeros(shape=(xlen, nrimgs), dtype=float)

    for imgnr in range(nrimgs):
        tmpspec = np.mean(data[imgnr * xlen:xlen + imgnr * xlen, ymaxind - dy:ymaxind + dy], axis=1)
        bg = np.mean(data[0:xlen, yminind:yminind + 3], axis=1)
        spec[:, imgnr] = tmpspec - bg

    return (lambd, spec, xlen, nrimgs)


def Gauss(x, a, x0, sigma):  # gaussian function
    return a * np.exp(-(x - x0) ** 2 / (2 * sigma ** 2))


def GaussFit(x, y, plotbool=False):  # fits gaussian to given curve
    # calculate starting values
    mean = sum(x * y) / sum(y)
    sigma = np.sqrt(sum(y * (x - mean) ** 2) / sum(y))

    # fitting
    try:
        popt, pcov = scipy.optimize.curve_fit(Gauss, x, y, p0=[max(y), mean, sigma], maxfev=10000)
    except RuntimeError:
        # print("Gaussian fit failed. Replaced data with NaN")
        popt = [float('nan'), float('nan'), float('nan')]

    # plot if plotbool == True
    if plotbool == True:
        plt.plot(x, y, 'b+:', label='data')
        plt.plot(x, Gauss(x, *popt), 'r-', label='fit')
        plt.legend()
        #plt.show()

    return (popt)  # [maximum intensity, maximum lambda, sigma]