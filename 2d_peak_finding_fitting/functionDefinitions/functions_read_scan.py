import struct
from math import sqrt
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Qt5Agg")
from matplotlib import pyplot as plt
import scipy, lmfit  # used for implementation of least-squares fit
from scipy.optimize import minimize  # used for implementation of maximum likelihood exponential fit
import time  # used to count run time of loops
import numba as nb
from matplotlib.colors import LogNorm
import pickle
import pylab as plb
from skimage.feature import blob_dog, blob_log, blob_doh
import os

# Custom Style jsanten
from matplotlib import rc
import seaborn as sns
sns.set()
sns.color_palette("bright")
sns.set_color_codes(palette="bright")
sns.set_style('ticks')
# Computer Modern Font (like in LaTeX)
#rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
#rc('text', usetex=True)


def load_obj(name, folder):
    with open(folder + name + '.pkl', 'rb') as f:
        return pickle.load(f)


def save_obj(obj, name, folder):
    with open(folder + name + '.pkl', 'wb') as f:
        pickle.dump(obj, f, pickle.HIGHEST_PROTOCOL)


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


def FindFirstPhotonAfterTrigger(macrotimes, macrotimesT):
    # finds the index of the first photon i of the array macrotimes after each event of macrotimesT (e.g. of trigger).
    # gives out the index of the first event in array macrotimes that happens after each trigger event in array macrotimesT
    # usually macrotimes = macrotimes of photons
    # macrotimesT = macrotimes of trigger
    # output firstPhot gives array with indices of the first photons per macrotimesT bin. Length is the same as len(macrotimesT)

    nrbins = len(macrotimesT)
    nrphot = len(macrotimes)
    firstPhot = np.full(nrbins, nrphot)  # initialize array and assign nr of last photon to all elements (for safety)
    counter, i = 0, 0
    while counter < nrbins and i < nrphot:  # for all trigger signals and photons
        while macrotimes[i] > macrotimesT[counter]:  # as long as photon i arrives after trigger counter do:
            firstPhot[counter] = i  # store
            counter += 1  # # increase trigger
            if counter == nrbins:
                break
        if counter == nrbins:
            break
        i += 1  # as long as photon arrived before trigger, go to next photon

    return (firstPhot)


def GetLifetime(microtimes, dtmicro, dtmacro, dtfit, tstart=-1, binwidth=1, ybg=0, plotbool=True, method='WLS'):
    # microtimes = microtimes array with photon events
    # dtfit is the time interval considered for the fit [s], tstart [s] is the starting point of the fit within the histogram. If set to -1 it starts at the time with the highest intensity.
    # binwidth is a multiplier. actual binwidth is given as binwidth*dtmicro[s]
    # ybg is the background considered for the fit (CHECK UNITS!!). If set to -1 --> try to estimate background based on last bins. set to 0 --> no background subtraction
    # plotbool: plot histogram with fit

    [ylist, xlist] = np.histogram(microtimes, int(dtmacro / (dtmicro * binwidth)), [0, int(dtmacro / dtmicro)])
    tlist = (xlist[:-1] + 0.5 * (xlist[1] - xlist[0])) * dtmicro * 1e9

    istart = int(tstart / dtmicro)  # find index of maximum element in ylist
    if istart < 0:
        istart = ylist.argmax()
    iend = istart + int(dtfit / (dtmicro * binwidth))
    if iend > len(tlist):
        iend = len(tlist)

        # get background (by simply looking at last ten data points) and substract from intensity data.
    if ybg < 0:
        ybg = np.mean(ylist[-10:])  # mean background per histogram bin bin of length

    if method == 'ML':  # maximum likelihood exponential fit
        [taufit, Afit] = MaxLikelihoodFit(tlist, ylist, istart, iend, ybg, False)
    if method == 'ML_c':  # maximum likelihood exponential fit
        [taufit, Afit] = MaxLikelihoodFit_c(tlist, ylist, istart, iend, ybg, False)
    elif method == 'WLS':  # weighted least squares fit
        [taufit, Afit] = WeightedLeastSquareFit(tlist, ylist, istart, iend, ybg, plotbool=False)
    else:
        taufit = 0;
        Afit = 0;
        print('Error: invalid fit method')

    if plotbool == True:
        fig = plt.figure()
        plt.xlabel('time (ns)')
        plt.ylabel('logy(??)')
        plt.semilogy(tlist, ylist, 'b.', alpha=0.2)
        plt.semilogy(tlist[istart:iend], Afit * np.exp(-(tlist[istart:iend] - tlist[istart]) / taufit) + ybg, "r")
        plt.semilogy([tlist[0], tlist[-1]], [ybg, ybg], 'k--')
        plt.tight_layout()
        plt.show()
        fig.savefig("png/" + name + "_lifetime_fit_" + method + ".png", dpi=300, transparent=False)
        fig.savefig(
            "pdf/" + name + "_lifetime_fit_" + method + ".pdf")  # TODO: change the plots to be unique somehow else
        print('Fitted lifetime:', taufit, 'ns; Amax:', Afit)

    # Amax is the maximum y-value
    Amax = np.max(ylist)

    return (taufit, Afit, ybg)


def MaxLikelihoodFit(tlist, ylist, istart, iend, bgcpb, plotbool=False):
    ### Maximum likelihood routine to fit single exponential. Pro: Works also for small amount of data (single bins of 10ms!)
    # tlist: x-axis values, here time in ns; ylist: y-axis values, here cts per tlist-bin; istart and iend: first and last element of tlist and ylist that are considered for the fit.

    # check if istart and iend are good numbers
    if istart < 0 or istart >= len(ylist):
        istart = 0
        print('WARNING: adapted istart in MaxLikelihoodExpFit')
    if iend <= istart or iend > len(ylist):
        iend = len(ylist)
        print('WARNING: adapted iend in MaxLikelihoodExpFit')

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
        # yest = np.array([Aest*np.exp(-(xdata[i]-xdata[0])/tauest)+bgcpb for i in range(len(xdata))])
        # plt.semilogy(tlist,ylist,'.',xdata,yest,[xdata[1],xdata[-1]],[bgcpb,bgcpb],'k--')
        # plt.show()
        # plt.savefig("TestSave.png")
        print("plotbool should have been executed")


def MaxLikelihoodFit_c(tlist, ylist, istart, iend, bgcpb, plotbool=False):
    ### Maximum likelihood routine to fit single exponential. Pro: Works also for small amount of data (single bins of 10ms!)
    # tlist: x-axis values, here time in ns; ylist: y-axis values, here cts per tlist-bin; istart and iend: first and last element of tlist and ylist that are considered for the fit.

    # check if istart and iend are good numbers
    if istart < 0 or istart >= len(ylist):
        istart = 0
        print('WARNING: adapted istart in MaxLikelihoodExpFit')
    if iend <= istart or iend > len(ylist):
        iend = len(ylist)
        print('WARNING: adapted iend in MaxLikelihoodExpFit')

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
        plt.figure()
        plt.plot(tlist, ylist, '.', xdata, yest, [xdata[1], xdata[-1]], [bgcpb, bgcpb], 'k--')
        plt.xlim([xdata[1], xdata[-1]])
        plt.show()

    return (tauest, Aest)


def MaxLikelihoodFunction(params, xdata, ydata, const):
    # max likelihood function for A*exp(-t/tau), needed in function MakLikelihoodFit
    # params = [A,tau]
    A = params[0]
    tau = params[1]
    E = 0;
    for i in range(len(xdata)):
        E = E + ydata[i] * np.log(A * np.exp(-(xdata[i] - xdata[0]) / tau) + const) - (
                    A * np.exp(-(xdata[i] - xdata[0]) / tau) + const)

    return (-E)  # This function needs to be MINIMIZED (because of the minus sign) to have the maximum likelihood fit!


def WeightedLeastSquareFit(tlist, ylist, istart, iend, bgcpb, plotbool=False):
    # check if istart and iend are good numbers
    if istart < 0 or istart >= len(ylist):
        istart = 0
        print('WARNING: adapted istart in MaxLikelihoodExpFit')
    if iend <= istart or iend > len(ylist):
        iend = len(ylist)
        print('WARNING: adapted iend in MaxLikelihoodExpFit')

    # shift t0 to t=0
    ydata = ylist[
            istart:iend] - bgcpb  # subtract constant bg %FIXME: jps. For uniform y, bgcpb is too high and therefore y data < 0

    ydatalog = np.log(np.abs(
        ydata))  # take log (weighting) #FIXME: jps. Why can ydata be negative here? (remove the abs value!!) or add the subtraction from the background
    xdata = tlist[istart:iend] - tlist[istart]  # shift to 0

    # do calculations
    mod = lmfit.models.LinearModel()
    pars = mod.guess(ydatalog, x=xdata)
    out = mod.fit(ydatalog, pars, x=xdata)
    Aest = np.exp(out.best_values['intercept'])
    tauest = -1.0 / out.best_values['slope']

    if plotbool == True:
        yest = np.array([Aest * np.exp(-(xdata[i]) / tauest) + bgcpb for i in range(len(xdata))])
        plt.semilogy(tlist, ylist, '.')
        plt.semilogy(xdata + tlist[istart], yest, [xdata[0], xdata[-1]], [bgcpb, bgcpb], 'k--')
        plt.show()
        print(out.fit_report())

    return (tauest, Aest)


def TimeGate(microtimes, macrotimes, dtmicro, tthreshold):
    # microtimes: input microtimes array
    # dtmicro: microtime unit [s]
    # tthreshold: threshold time [s]
    # returns array with all microtimes below threshold and array with all microtimes above threshold

    microtimesbelow = microtimes[microtimes < int(tthreshold / dtmicro)]
    microtimesabove = microtimes[microtimes >= int(tthreshold / dtmicro)]
    macrotimesbelow = macrotimes[microtimes < int(tthreshold / dtmicro)]
    macrotimesabove = macrotimes[microtimes >= int(tthreshold / dtmicro)]

    return (microtimesbelow, macrotimesbelow, microtimesabove, macrotimesabove)


def GetHistogram(microtimes, dtmicro, dtmacro, binning):
    [ylist, xlist] = np.histogram(microtimes, int(dtmacro / dtmicro / binning), [0, int(dtmacro / dtmicro)])
    tlist = (xlist[:-1] + 0.5 * (xlist[1] - xlist[0])) * dtmicro
    return (ylist, tlist)