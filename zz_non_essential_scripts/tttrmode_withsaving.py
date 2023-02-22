# -*- coding: utf-8 -*-
"""
Created on Tue Apr  2 10:41:40 2019

@author: Robert Keitel
"""
#
# Keno Goertz, PicoQuant GmbH, February 2018

import time
import ctypes as ct
from ctypes import byref
import sys
import struct
import pickle


def save_obj(obj, name ):
    with open('output/'+ name + '.pkl', 'wb') as f:
        pickle.dump(obj, f, pickle.HIGHEST_PROTOCOL)

def load_obj(name ):
    with open('output/' + name + '.pkl', 'rb') as f:
        return pickle.load(f)

# From hhdefin.h
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
HHsettings["binning"] = 2 # you can change this, meaningful only in T3 mode
HHsettings["offset"] = 0 # you can change this, meaningful only in T3 mode
HHsettings["tacq"] = 10000 # Measurement time in millisec, you can change this
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

#%%
hhlib = ct.CDLL("hhlib64.dll")

def closeDevices():
    for i in range(0, MAXDEVNUM):
        hhlib.HH_CloseDevice(ct.c_int(i))
    #exit(0)

def stoptttr():
    retcode = hhlib.HH_StopMeas(ct.c_int(dev[0]))
    if retcode < 0:
        print("HH_StopMeas error %1d. Aborted." % retcode)
    closeDevices()

def tryfunc(retcode, funcName, measRunning=False):
    if retcode < 0:
        hhlib.HH_GetErrorString(errorString, ct.c_int(retcode))
        print("HH_%s error %d (%s). Aborted." % (funcName, retcode,\
              errorString.value.decode("utf-8")))
        if measRunning:
            stoptttr()
        else:
            closeDevices()

hhlib.HH_GetLibraryVersion(libVersion)
print("Library version is %s" % libVersion.value.decode("utf-8"))
if libVersion.value.decode("utf-8") != LIB_VERSION:
    print("Warning: The application was built for version %s" % LIB_VERSION)

outputfile = open("output/tttrmode.out", "wb+")

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

progress = 0
sys.stdout.write("\nProgress:%9u" % progress)
sys.stdout.flush()
running=True

tryfunc(hhlib.HH_StartMeas(ct.c_int(dev[0]), ct.c_int(HHsettings["tacq"])), "StartMeas")

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
        sys.stdout.write("\rProgress:%9u" % progress)
        sys.stdout.flush()
    else:
        tryfunc(hhlib.HH_CTCStatus(ct.c_int(dev[0]), byref(ctcstatus)),\
                "CTCStatus")
        if ctcstatus.value > 0: 
            print("\nDone")
            stoptttr()
            running=False
    # within this loop you can also read the count rates if needed.
HHsettings["overallCounts"]=progress
save_obj(HHsettings,"settingsfile")
closeDevices()
outputfile.close()