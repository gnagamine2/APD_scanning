# -*- coding: utf-8 -*-
"""
Created on Mon Apr 15 13:31:33 2019
Takes the output from read scan and does a TTTR measurement on each particle
@author: omel-acton
"""

import struct
from math import sqrt
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import scipy, lmfit # used for implementation of least-squares fit
from scipy.optimize import minimize # used for implementation of maximum likelihood exponential fit
import time # used to count run time of loops
import numba as nb
from matplotlib.colors import LogNorm
import pickle
from multiprocessing import Process, Pipe,Value
from time import sleep
import time
import numpy as np
from matplotlib import pyplot as plt
from ctypes import * #C is needed to use DLL. madlib.dll has to be in same folder as python.exe
import timeit
import ctypes as ct
import sys
import struct
import pickle
import cv2

# Imports for the LabFlipper
import os

# Import for andor
import nidaqmx

########################################################
# Data Import Here --> Need to change paths
########################################################
folder = "C:/Users/omel-acton/Desktop/Gabriel_Local_EMCCD/22-11-18/"
Mname="AP-3-108_C1x10-6_SingleDot_1x10+6Hz_405exct_220nW_x_49_y78_z124_range40_pts200"
MEASUREMENT_TIME = 300  # [seconds] for intensity tracing
exposuretimeandor = 21  # [seconds] Define the same exposure time as in Andor

########################################################
# LEAVE UNCHANGED
########################################################

#%% Define functions
def load_obj(name, folder ):
    with open(folder + name + '.pkl', 'rb') as f:
        return pickle.load(f)
    
def save_obj(obj, name, folder ):
    with open(folder + name + '.pkl', 'wb') as f:
        pickle.dump(obj, f, pickle.HIGHEST_PROTOCOL)

def MCLReadpos(handle):
    xval=mcldll.MCL_SingleReadN(c_uint(1), handle)
    yval=mcldll.MCL_SingleReadN(c_uint(2), handle)
    zval=mcldll.MCL_SingleReadN(c_uint(3), handle)
    pos=[xval,yval,zval]
    return pos
def closeDevices():
    for i in range(0, MAXDEVNUM):
        hhlib.HH_CloseDevice(ct.c_int(i))
    #exit(0)

def stoptttr():
    retcode = hhlib.HH_StopMeas(ct.c_int(dev[0]))
    if retcode < 0:
        print("HH_StopMeas error %1d. Aborted." % retcode)
#    closeDevices()

def tryfunc(retcode, funcName, measRunning=False):
    if retcode < 0:
        hhlib.HH_GetErrorString(errorString, ct.c_int(retcode))
        print("HH_%s error %d (%s). Aborted." % (funcName, retcode,\
              errorString.value.decode("utf-8")))
        if measRunning:
            stoptttr()
        else:
            closeDevices()
#%%Initialize Stage
mcldll = cdll.madlib
# For all remaining calls do : mcldll.FunctionName(function parameters)
handle = mcldll.MCL_InitHandle()
if handle == 0:
    print("Error: Handle not intialized correctly\nExiting")
    sys.exit()

# find some basic information about the NanoDrive
class ProductInfo(Structure):
    _fields_ = [("axis_bitmap", c_ubyte),
                ("ADC_resolution", c_short),
                ("DAC_resolution", c_short),
                ("Product_id", c_short),
                ("FirmwareVersion", c_short),
                ("FirmwareProfile", c_short)]
    _pack_ = 1 # this is how it is packed in the Madlib dll

pi = ProductInfo()
ppi = pointer(pi)

err = mcldll.MCL_GetProductInfo(ppi, handle)
if err != 0:
    print("Error: NanoDrive could not get productInformation. Error Code:", err, "\nExiting")
    mcldll.MCL_ReleaseHandle(handle) # be sure to release handle anytime before returning
    sys.exit()
else:
    print("Information about the NanoDrive:")
    print("axis bitmap:", pi.axis_bitmap)
    print("ADC resolution:", pi.ADC_resolution)
    print("DAC resolution:", pi.DAC_resolution)
    print("Product ID:", pi.Product_id)
    print("Firmware Version:", pi.FirmwareVersion)
    print("Firmware Profile:", pi.FirmwareProfile)
cal = mcldll.MCL_GetCalibration
cal.restype = c_double
readpos = mcldll.MCL_SingleReadN
readpos.restype = c_double
xaxis = c_uint(1)
yaxis = c_uint(2)
zaxis = c_uint(3)
mcldll.MCL_IssBindClockToAxis(c_uint(1), c_uint(2), xaxis, handle)
xcalibration = mcldll.MCL_GetCalibration(xaxis, handle)
print("axis calibration range is ",xcalibration)
ycalibration = mcldll.MCL_GetCalibration(yaxis, handle)
zcalibration = mcldll.MCL_GetCalibration(zaxis, handle)
print("axis calibration range is ",xcalibration,ycalibration,zcalibration)
inpos=MCLReadpos(handle)
print("Current position is ",inpos[0],inpos[1],inpos[2])
print("Initialization Successfull")

#%% Initialize Hydraharp
LIB_VERSION = "3.0"
MAXDEVNUM = 8
MODE_T2 = 2
MODE_T3 = 3
MAXLENCODE = 6
HHMAXINPCHAN = 8
TTREADMAX = 131072
FLAG_OVERFLOW = 0x0001
FLAG_FIFOFULL = 0x0002

HHsettings={}
# Measurement parameters, these are hardcoded since this is just a demo
mode = MODE_T3 # set T2 or T3 here, observe suitable Syncdivider and Range!
HHsettings["binning"] = 10 # you can change this, meaningful only in T3 mode; 32768 ps*2^binning NEEDS TO BE AT LEAST THE TIME BETWEEN PULSES TO CAPTURE ALL THE DATA
HHsettings["offset"] = 0 # you can change this, meaningful only in T3 mode
HHsettings["tacq"] = MEASUREMENT_TIME * 1000 # Measurement time in millisec, you can change this
HHsettings["syncDivider"] = 2 # you can change this, observe mode! READ MANUAL!; Check in software!!! In software it is 8
HHsettings["syncCFDZeroCross"] = 10 # you can change this (in mV); Check in software!!! - correct
HHsettings["syncCFDLevel"] = 80 # you can change this (in mV); Check in software!!! - correct
HHsettings["syncChannelOffset"] = 48000 # you can change this (in ps, like a cable delay); Check in software!!!
HHsettings["inputCFDZeroCross"] = 19 # you can change this (in mV); Check in software!!!
HHsettings["inputCFDLevel"] = 200 # you can change this (in mV); Check in software!!!
HHsettings["inputChannelOffset"] = -500 # you can change this (in ps, like a cable delay); Check in software!!!

# Variables to store information read from DLLs
buffer = (ct.c_uint * TTREADMAX)()
dev = []
libVersion = ct.create_string_buffer(b"", 8)
hwSerial = ct.create_string_buffer(b"", 8)
hwPartno = ct.create_string_buffer(b"", 8)
hwVersion = ct.create_string_buffer(b"", 8)
hwModel = ct.create_string_buffer(b"", 16)
errorString = ct.create_string_buffer(b"", 40)
numChannels = ct.c_int()
resolution = ct.c_double()
syncRate = ct.c_int()
countRate = ct.c_int()
flags = ct.c_int()
nRecords = ct.c_int()
ctcstatus = ct.c_int()
warnings = ct.c_int()
warningstext = ct.create_string_buffer(b"", 16384)

hhlib = ct.CDLL("hhlib64.dll")
hhlib.HH_GetLibraryVersion(libVersion)
print("Library version is %s" % libVersion.value.decode("utf-8"))
if libVersion.value.decode("utf-8") != LIB_VERSION:
    print("Warning: The application was built for version %s" % LIB_VERSION)

print("\nMode             : %d" % mode)
print("Binning           : %d" % HHsettings["binning"])
print("Offset            : %d" % HHsettings["offset"])
print("AcquisitionTime   : %d" % HHsettings["tacq"])
print("SyncDivider       : %d" % HHsettings["syncDivider"])
print("SyncCFDZeroCross  : %d" % HHsettings["syncCFDZeroCross"])
print("SyncCFDLevel      : %d" % HHsettings["syncCFDLevel"])
print("InputCFDZeroCross : %d" % HHsettings["inputCFDZeroCross"])
print("InputCFDLevel     : %d" % HHsettings["inputCFDLevel"])

print("\nSearching for HydraHarp devices...")
print("Devidx     Status")

for i in range(0, MAXDEVNUM):
    retcode = hhlib.HH_OpenDevice(ct.c_int(i), hwSerial)
    if retcode == 0:
        print("  %1d        S/N %s" % (i, hwSerial.value.decode("utf-8")))
        dev.append(i)
    else:
        if retcode == -1: # HH_ERROR_DEVICE_OPEN_FAIL
            print("  %1d        no device" % i)
        else:
            hhlib.HH_GetErrorString(errorString, ct.c_int(retcode))
            print("  %1d        %s" % (i, errorString.value.decode("utf8")))

# In this demo we will use the first HydraHarp device we find, i.e. dev[0].
# You can also use multiple devices in parallel.
# You can also check for specific serial numbers, so that you always know 
# which physical device you are talking to.

if len(dev) < 1:
    print("No device available.")
    closeDevices()
print(dev)
print("Using device #%1d" % dev[0])
print("\nInitializing the device...")

# with internal clock; CHECK THAT THIS IS GOOD AND MAKES SENSE
tryfunc(hhlib.HH_Initialize(ct.c_int(dev[0]), ct.c_int(mode), ct.c_int(0)),\
        "Initialize")

# Only for information
tryfunc(hhlib.HH_GetHardwareInfo(dev[0], hwModel, hwPartno, hwVersion),\
        "GetHardwareInfo")
print("Found Model %s Part no %s Version %s" % (hwModel.value.decode("utf-8"),\
      hwPartno.value.decode("utf-8"), hwVersion.value.decode("utf-8")))

tryfunc(hhlib.HH_GetNumOfInputChannels(ct.c_int(dev[0]), byref(numChannels)),\
        "GetNumOfInputChannels")
print("Device has %i input channels." % numChannels.value)

print("\nCalibrating...")
tryfunc(hhlib.HH_Calibrate(ct.c_int(dev[0])), "Calibrate")
tryfunc(hhlib.HH_SetSyncDiv(ct.c_int(dev[0]), ct.c_int(HHsettings["syncDivider"])), "SetSyncDiv")

tryfunc(
    hhlib.HH_SetSyncCFD(ct.c_int(dev[0]), ct.c_int(HHsettings["syncCFDLevel"]),
                        ct.c_int(HHsettings["syncCFDZeroCross"])),\
    "SetSyncCFD"
)
tryfunc(hhlib.HH_SetSyncChannelOffset(ct.c_int(dev[0]), ct.c_int(HHsettings["syncChannelOffset"])),\
        "SetSyncChannelOffset")

# we use the same input settings for all channels, you can change this
for i in range(0, numChannels.value):
    tryfunc(
        hhlib.HH_SetInputCFD(ct.c_int(dev[0]), ct.c_int(i), ct.c_int(HHsettings["inputCFDLevel"]),\
                             ct.c_int(HHsettings["inputCFDZeroCross"])),\
        "SetInputCFD"
    )

    tryfunc(
        hhlib.HH_SetInputChannelOffset(ct.c_int(dev[0]), ct.c_int(i),\
                                       ct.c_int(HHsettings["inputChannelOffset"])),\
        "SetInputChannelOffset"
    )
# Meaningful only in T3 mode
if mode == MODE_T3:
    tryfunc(hhlib.HH_SetBinning(ct.c_int(dev[0]), ct.c_int(HHsettings["binning"])), "SetBinning")
    tryfunc(hhlib.HH_SetOffset(ct.c_int(dev[0]), ct.c_int(HHsettings["offset"])), "SetOffset")
    tryfunc(hhlib.HH_GetResolution(ct.c_int(dev[0]), byref(resolution)), "GetResolution")
    HHsettings["resolution"]=resolution.value
    print("Resolution is %1.1lfps" % resolution.value)

# Note: after Init or SetSyncDiv you must allow >100 ms for valid  count rate readings
time.sleep(0.2)

tryfunc(hhlib.HH_GetSyncRate(ct.c_int(dev[0]), byref(syncRate)), "GetSyncRate")
HHsettings["syncRate"]=syncRate.value
print("\nSyncrate=%1d/s" % syncRate.value)
for i in range(0, numChannels.value):
    tryfunc(hhlib.HH_GetCountRate(ct.c_int(dev[0]), ct.c_int(i), byref(countRate)),\
            "GetCountRate")
    print("Countrate[%1d]=%1d/s" % (i, countRate.value))


#%% Initialize MFF102 (Mirror flip)

os.chdir(r"C:\Program Files\Thorlabs\Kinesis")
lib = cdll.LoadLibrary("Thorlabs.MotionControl.FilterFlipper.dll")

# Build device list
# This function builds an internal collection of all devices found on the USB that are not currently open. <br />
# NOTE, if a device is open, it will not appear in the list until the device has been closed.
lib.TLI_BuildDeviceList()

# set up serial number variable
# manually replace serial number of specific device
serialNumber_MFF102 = c_char_p(b"37001165")
moveTimeout = 60.0

# set up device
lib.FF_Open(serialNumber_MFF102)
lib.FF_StartPolling(serialNumber_MFF102, c_int(200))

time.sleep(3)

lib.FF_ClearMessageQueue(serialNumber_MFF102)

# homing
lib.FF_Home(serialNumber_MFF102)

#%% Define Settings for Andor

# intialize camera (connect NI box)
tacq = exposuretimeandor;  # this is set at ANDOR
T = [True]  # True starts it
F = [False]  # True ends it

#%% Load settings

# Check if Necessary Folders are created (namely: /measurements)
measurements_folder = os.path.join(folder, 'measurements')

if not (os.path.exists(measurements_folder)): #if folder not found
    os.mkdir(measurements_folder)
    print("measurements folder created: ", measurements_folder)
else:
    print("measurements folder found: ", measurements_folder)

savefolder = folder + "measurements/"

name = Mname+'_tttrmode'
MCLdata=load_obj(Mname+'_MCLdata', folder)

#%% Measure each particle
with nidaqmx.Task() as task:
    task.do_channels.add_do_chan('Dev1/port0/line2')

    for j in range(MCLdata["ParticleLocationsIndices"].shape[0]):
        
        
        index0, index1, r=MCLdata["ParticleLocationsIndices"][j]
        index0=int(index0)
        index1=int(index1)
        mcldll.MCL_SingleWriteN(c_double(MCLdata["xposmatrix"][index0,index1]), xaxis, handle)    #TODO index0 and index1 might have to be exchanged
        mcldll.MCL_SingleWriteN(c_double(MCLdata["yposmatrix"][index0,index1]), yaxis, handle)
        sleep(200/1000)
        xposcurr=mcldll.MCL_SingleReadN(xaxis, handle)
        yposcurr=mcldll.MCL_SingleReadN(yaxis, handle)
        print("Starting measurement at x="+str(xposcurr)+", y="+str(yposcurr))
        print("Target position was x="+str(MCLdata["xposmatrix"][index0,index1])+", y="+str(MCLdata["yposmatrix"][index0,index1]))
        
        
        # JSANTEN: Put the mirror in the side position to measure the spectrum (TBD)
        lib.FF_MoveToPosition(serialNumber_MFF102, 1)  # move to position 2 (to the side); 1 (up)
        print('Moving Mirror to the Up Position')  # send move command
        sleep(2)

        # Measure the Spectrum
        task.write(T)
        time.sleep(tacq + 0.1)
        task.write(F)
        #sleep(0.5)
        
        lib.FF_MoveToPosition(serialNumber_MFF102, 2)
        print('Moving Mirror to the Side Position')  # send move command
        sleep(0.5)
        #JSANTEN END
        

       
        # Start measuring
        savename=Mname+"_Particle_"+str(j)+"_x"+str(round(MCLdata["xposmatrix"][index0,index1]))+"y_"+str(round(MCLdata["yposmatrix"][index0,index1]))
        outputfile = open(savefolder+savename+"_tttrmode.out", "wb+")
        progress=0
        tryfunc(hhlib.HH_GetCountRate(ct.c_int(dev[0]), ct.c_int(0), byref(countRate)),\
                "GetCountRate")
        print("Countrate[%1d]=%1d/s" % (i, countRate.value))
        sys.stdout.write("\nProgress:%9u" % progress)
        sys.stdout.flush()
        running=True
        tryfunc(hhlib.HH_StartMeas(ct.c_int(dev[0]), ct.c_int(HHsettings["tacq"])), "StartMeas")
        counts=0;
        pixel=-1;
        while running==True:

            tryfunc(hhlib.HH_GetFlags(ct.c_int(dev[0]), byref(flags)), "GetFlags")

            if flags.value & FLAG_FIFOFULL > 0:
                print("\nFiFo Overrun!")
                stoptttr()

            tryfunc(
                    hhlib.HH_ReadFiFo(ct.c_int(dev[0]), byref(buffer), TTREADMAX,\
                                      byref(nRecords)),\
                                      "ReadFiFo", measRunning=True
                                      )

            if nRecords.value > 0:
            # We could just iterate through our buffer with a for loop, however,
            # this is slow and might cause a FIFO overrun. So instead, we shrinken
            # the buffer to its appropriate length with array slicing, which gives
            # us a python list. This list then needs to be converted back into
            # a ctype array which can be written at once to the output file
                outputfile.write((ct.c_uint*nRecords.value)(*buffer[0:nRecords.value]))
                progress += nRecords.value
    #            for i in range(nRecords.value):
    #                entry = f.read(4)
    #                channel = struct.unpack("I",entry)[0] >> 25 # read channel number, first 7 bits
    #                if channel == 0:
    #                    counts += 1
    #                elif channel == 1: # APD2: not implemented yet.
    #                    counts += 1


                counts += nRecords.value
                sys.stdout.write("\rProgress:%9u" % progress)
                sys.stdout.flush()
            else:
                tryfunc(hhlib.HH_CTCStatus(ct.c_int(dev[0]), byref(ctcstatus)),\
                        "CTCStatus")
                if ctcstatus.value > 0:     #I should check if I can avoid this and just kill the loop when the stage is done
                    print("\nAcq time over")
                    stoptttr()
                    running=False

        print("terminating hydraharp")
        HHsettings["overallCounts"]=progress
        save_obj(HHsettings,savename+"_HHsettings",savefolder)
        outputfile.close()
        print("Measurement "+str(j+1)+" of "+str(MCLdata["ParticleLocationsIndices"].shape[0])+" done")
        #mcldll.MCL_SingleWriteN(c_double(xtargetvals[0]+((-0.5-Points/2)/Points)*ScanRange), xaxis, handle)

    #%% Shutdown
    closeDevices()
    mcldll.MCL_ReleaseHandle(handle)

    
    # Shutdown for MFF102 (Mirror Flip)
    lib.FF_ClearMessageQueue(serialNumber_MFF102) #clean up and exit
    lib.FF_StopPolling(serialNumber_MFF102)
    lib.FF_Close(serialNumber_MFF102)
    


print("Program completed")