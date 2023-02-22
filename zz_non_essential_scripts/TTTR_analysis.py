# -*- coding: utf-8 -*-
"""
Created on Tue Jan 15 17:50:02 2019
TTTR analysis code with marker from Felipe
@author: rober
"""

import struct
import numpy as np
import pandas as pd

def importT3marked(filename):
    """Import ptu files from HydraHarp with Marker inputs for frame synchronization.
    
    Parameters
    ----------
    filename: String with the complete path to the file.

    Returns
    -------
    raw_data : 2D numpy array containing the counts [NxN]
    x: horizontal pixel number (image) or corresponding wavelength (spectrum), depending on input file type [Nx1]
    y: vertical pixel number [Nx1]
    info: dictionary with additional information such as exposure time
    
    channel0: DataFrame [Nx6] ['frame', 'macrotime', 'macrotime_s', 'microtime', 'microtime_s', 'totaltime_s']
    channel1: DataFrame [Nx6] ['frame', 'macrotime', 'macrotime_s', 'microtime', 'microtime_s', 'totaltime_s']
    dt1: macrotime resolution in ns (sync channel repetition rate)
    dt2: microtime resolution in ns (photon arrival time after pulse)

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
    with open(filename, "rb") as f:

        while True:
            if f.read(1) == b'M':
                if f.read(18) == b'easDesc_Resolution': # recognize time unit entry
                    break
        f.read(21) # rest of the tag
        dt1 = struct.unpack('d',f.read(8))[0]
        #print('Microtime unit:', dt1)

        while True:
            if f.read(1) == b'M':
                if f.read(24) == b'easDesc_GlobalResolution': # recognize time unit entry
                    break
        f.read(15) # rest of the tag
        dt2 = struct.unpack('d',f.read(8))[0]
        #print('Macrotime unit:', dt2)

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

        macrotimes0 = np.zeros(nrrec,dtype='int64');
        microtimes0 = np.zeros(nrrec,dtype='int64');
        framenr0 = np.zeros(nrrec,dtype='int64');
        pixelnr0 = np.zeros(nrrec,dtype='int64');
        linenr0 = np.zeros(nrrec,dtype='int64');
        accphotons = np.zeros(nrrec,dtype='int64');
        photonsppx = np.zeros(nrrec,dtype='int64');
        # Optional as matrix: photonsppx = np.zeros((points,points),dtype='int64');

        macrotimes1 = np.zeros(nrrec,dtype='int64');
        microtimes1 = np.zeros(nrrec,dtype='int64');
        framenr1 = np.zeros(nrrec,dtype='int64');

        macrotimesM1 = np.zeros(nrrec,dtype='int64');
        macrotimesM2 = np.zeros(nrrec,dtype='int64');

        overflows = 0
        nrphotons0 = 0
        nrphotons1 = 0
        nrM1 = 0
        nrM2 = 0
        pxN1 = 0
        prevchann = 0

        for i in range(nrrec):
            entry = f.read(4)
            channel = struct.unpack("I",entry)[0] >> 25 # read channel number, first 7 bits
            if channel == 0:
                if nrM1>0:
                    macrotime = (struct.unpack("I",entry)[0] & 0x3FF)
                    macrotimes0[nrphotons0] = macrotime + 1024*overflows
                    microtime = ((struct.unpack("I",entry)[0] >> 10) & 0x7FFF)
                    microtimes0[nrphotons0] = microtime
                    framenr0[nrphotons0] = nrM1
                    pixelnr0[nrphotons0] = pxN1
                    linenr0[nrphotons0] = nrM2
                    nrphotons0 += 1
                    prevchann = 0
            elif channel == 65:
                accphotons[nrM1]=nrphotons0
                if nrM1>0:
                    photonsppx[nrM1]=accphotons[nrM1]-accphotons[nrM1-1]
                    #as matrix: photonsppx[nrM1,nrM2]=accphotons[nrM1,nrM2]-accphotons[nrM1-1,nrM2]
                else:
                    photonsppx[nrM1,nrM2]=accphotons[nrM1,nrM2]
                    
                nrM1 += 1
                pxN1 += 1
                macrotime = (struct.unpack("I",entry)[0] & 0x3FF)
                macrotimesM1[nrM1] = macrotime + 1024*overflows
                prevchann = 65
            elif channel == 66:
                nrM2 += 1
                macrotime = (struct.unpack("I",entry)[0] & 0x3FF)
                macrotimesM2[nrM2] = macrotime + 1024*overflows
                prevchann = 66
                pxN1 = 0
            elif channel == 127:
                nroverflows = (struct.unpack("I",entry)[0] & 0x3FF)
                overflows += nroverflows
                prevchann = 127
            else:
                print('bad')
                print(channel)

    microtimes0 = microtimes0[:nrphotons0]
    macrotimes0 = macrotimes0[:nrphotons0]
    framenr0 = framenr0[:nrphotons0]
    framenr1 = framenr1[:nrphotons1]

    microtimes1 = microtimes1[:nrphotons1]
    macrotimes1 = macrotimes1[:nrphotons1]

    macrotimesM1 = macrotimesM1[:nrM1]
    macrotimesM2 = macrotimesM2[:nrM2]
    
    #channel0 = pd.DataFrame({'macrotime': macrotimes0-macrotimesM1[0], 'macrotime_s': dt2*(macrotimes0-macrotimesM1[0]), 'microtime': microtimes0, 'microtime_s': dt1*microtimes0, 'totaltime_s' : dt2*(macrotimes0-macrotimesM1[0])+dt1*microtimes0, 'frame': framenr0, 'line': linenr0, 'pixel': pixelnr0, 'frame': framenr0})
    
    
    #channel0 = channel0.loc[channel0['frame']>=0]
    #channel0.index = np.arange(len(channel0))

    print(nrphotons0,nrM1,nrM2,overflows)
    
    print('Macrotime resolution: '+str(dt2*1e9)+' ns')
    print('Microtime resolution: '+str(int(dt1*1e12))+' ps')
    return [photonsppx]
    #return [channel0, channel1, dt1, dt2]
    
