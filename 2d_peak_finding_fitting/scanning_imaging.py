# -*- coding: utf-8 -*-
"""
Created on Thu Jan 10 16:17:52 2019
Script uses the DLL for the MCL stage and the Hydraharp to raster scan the stage while recording photons in TTTR mode. 
The stage sends out a marker whenever x-position is read which is used to plot the image.
Stage works so far, HydraHarp part has to be tested.
Based on the madcitylabs sample files and for the HydraHarp Keno Goertz demo.

Julian Santen for BSc Thesis (jsanten@ethz.ch, from Sept. 2022)
Reference for code used for HydraHarp integration:
https://github.com/PicoQuant/HH400-v3.x-Demos/tree/master/Windows/64bit/Python

@author: Robert Keitel
"""
from multiprocessing import Process, Pipe,Value
from time import sleep
import time
import numpy as np
from matplotlib import pyplot as plt
from ctypes import * #C is needed to use DLL. madlib.dll has to be in same folder as python.exe
import timeit
from time import sleep
import ctypes as ct
from ctypes import byref
import sys
import struct
import pickle
import cv2

#%% Define functions for Piezo stage
def MCLReadpos(handle):
    xval=mcldll.MCL_SingleReadN(c_uint(1), handle)
    yval=mcldll.MCL_SingleReadN(c_uint(2), handle)
    zval=mcldll.MCL_SingleReadN(c_uint(3), handle)
    pos=[xval,yval,zval]
    return pos
def save_obj(obj, name, folder ):
    with open(folder + name + '.pkl', 'wb') as f:
        pickle.dump(obj, f, pickle.HIGHEST_PROTOCOL)

def load_obj(name ):
    with open('output/' + name + '.pkl', 'rb') as f:
        return pickle.load(f)
#%%Define while loop (reads data from Hydraharp and plots)
        #It is not clear if I need to pass stuff like the hhlib and the pointers
def HHloop(conn,MeasurementRunning,HHstatus,currpix,OutputFolder,SampleName):
    hhlib = ct.CDLL("hhlib64.dll") #load the library
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
                conn.close()
                if not outputfile.closed:
                    outputfile.close()
            else:
                closeDevices()
                conn.close()
                if not outputfile.closed:
                    outputfile.close()
            HHstatus.value=2
        
    def save_obj(obj, name, folder ):
        with open(folder + name + '.pkl', 'wb') as f:
            pickle.dump(obj, f, pickle.HIGHEST_PROTOCOL)
    
    def load_obj(name ):
        with open('output/' + name + '.pkl', 'rb') as f:
            return pickle.load(f)    
    LIB_VERSION = "3.0"
    MAXDEVNUM = 8
    MODE_T3 = 3
    MAXLENCODE = 6
    HHMAXINPCHAN = 8
    TTREADMAX = 131072
    FLAG_OVERFLOW = 0x0001
    FLAG_FIFOFULL = 0x0002
    
    HHsettings={}
    # Measurement parameters, these are hardcoded since this is just a demo
    mode = MODE_T3 # set T2 or T3 here, observe suitable Syncdivider and Range!
    HHsettings["binning"] = 10 # you can change this, meaningful only in T3 mode FOR POLARIZATION MODULATION WITH 50 kHz use 10. The value is the exponent with base 2 of the binning time. You need at least log2(cycleduration/32768)
    HHsettings["offset"] = 0 # you can change this, meaningful only in T3 mode; 
    HHsettings["tacq"] = 20000000 # Measurement time in millisec, you can change this # TODO: Check this value!
    HHsettings["syncDivider"] = 8 # you can change this, observe mode! READ MANUAL!; Check in software!!! In software it is 8
    HHsettings["syncCFDZeroCross"] = 10 # you can change this (in mV); Check in software!!! - correct
    HHsettings["syncCFDLevel"] = 190 # you can change this (in mV); Check in software!!! - correct
    HHsettings["syncChannelOffset"] = 40000 # you can change this (in ps, like a cable delay); Check in software!!!
    HHsettings["inputCFDZeroCross"] = 4 # you can change this (in mV); Check in software!!!
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
    
    hhlib.HH_GetLibraryVersion(libVersion)
    print("Library version is %s" % libVersion.value.decode("utf-8"))
    tolog="\nLibrary version is %s" % libVersion.value.decode("utf-8")
    if libVersion.value.decode("utf-8") != LIB_VERSION:
        print("Warning: The application was built for version %s" % LIB_VERSION)
        tolog=tolog+"\nWarning: The application was built for version %s" % LIB_VERSION
    outputfile = open(OutputFolder+SampleName+"_tttrmode.out", "wb+")
    
    print("\nMode             : %d" % mode)
    tolog=tolog+"\nMode             : %d" % mode
    print("\nBinning           : %d" % HHsettings["binning"])
    tolog=tolog+"\nBinning           : %d" % HHsettings["binning"]
    print("\nOffset            : %d" % HHsettings["offset"])
    tolog=tolog+"\nOffset            : %d" % HHsettings["offset"]
    print("\nAcquisitionTime   : %d" % HHsettings["tacq"])
    tolog=tolog+"\nAcquisitionTime   : %d" % HHsettings["tacq"]
    print("\nSyncDivider       : %d" % HHsettings["syncDivider"])
    tolog=tolog+"\nSyncDivider       : %d" % HHsettings["syncDivider"]
    print("\nSyncCFDZeroCross  : %d" % HHsettings["syncCFDZeroCross"])
    tolog=tolog+"\nSyncCFDZeroCross  : %d" % HHsettings["syncCFDZeroCross"]
    print("\nSyncCFDLevel      : %d" % HHsettings["syncCFDLevel"])
    tolog=tolog+"\nSyncCFDLevel      : %d" % HHsettings["syncCFDLevel"]
    print("\nInputCFDZeroCross : %d" % HHsettings["inputCFDZeroCross"])
    tolog=tolog+"\nInputCFDZeroCross : %d" % HHsettings["inputCFDZeroCross"]
    print("\nInputCFDLevel     : %d" % HHsettings["inputCFDLevel"])
    tolog=tolog+"InputCFDLevel     : %d" % HHsettings["inputCFDLevel"]
    
    
    print("\nSearching for HydraHarp devices...")
    tolog=tolog+"\nSearching for HydraHarp devices..."
    print("Devidx     Status")
    tolog=tolog+"\nDevidx     Status"
    for i in range(0, MAXDEVNUM):
        retcode = hhlib.HH_OpenDevice(ct.c_int(i), hwSerial)
        if retcode == 0:
            print("  %1d        S/N %s" % (i, hwSerial.value.decode("utf-8")))
            tolog=tolog+"\n  %1d        S/N %s" % (i, hwSerial.value.decode("utf-8"))
            dev.append(i)
        else:
            if retcode == -1: # HH_ERROR_DEVICE_OPEN_FAIL
                print("  %1d        no device" % i)
                tolog=tolog+"\n  %1d        no device" % i
            else:
                hhlib.HH_GetErrorString(errorString, ct.c_int(retcode))
                print("  %1d        %s" % (i, errorString.value.decode("utf8")))
                tolog=tolog+"\n  %1d        %s" % (i, errorString.value.decode("utf8"))
    
    # In this demo we will use the first HydraHarp device we find, i.e. dev[0].
    # You can also use multiple devices in parallel.
    # You can also check for specific serial numbers, so that you always know 
    # which physical device you are talking to.
    
    if len(dev) < 1:
        print("No device available.")
        tolog=tolog+"\nNo device available."
        closeDevices()
        if not outputfile.closed:
            outputfile.close()
        conn.send(tolog)
        HHstatus.value=2
    print(dev)
    
    print("Using device #%1d" % dev[0])
    tolog=tolog+"\nUsing device #%1d" % dev[0]
    print("\nInitializing the device...")
    tolog=tolog+"\nInitializing the device..."
    
    # with internal clock; TODO: CHECK THAT THIS IS GOOD AND MAKES SENSE
    tryfunc(hhlib.HH_Initialize(ct.c_int(dev[0]), ct.c_int(mode), ct.c_int(0)),\
            "Initialize")
    
    # Only for information
    tryfunc(hhlib.HH_GetHardwareInfo(dev[0], hwModel, hwPartno, hwVersion),\
            "GetHardwareInfo")
    print("Found Model %s Part no %s Version %s" % (hwModel.value.decode("utf-8"),\
          hwPartno.value.decode("utf-8"), hwVersion.value.decode("utf-8")))
    tolog=tolog+"\nFound Model %s Part no %s Version %s" % (hwModel.value.decode("utf-8"),\
          hwPartno.value.decode("utf-8"), hwVersion.value.decode("utf-8"))
    tryfunc(hhlib.HH_GetNumOfInputChannels(ct.c_int(dev[0]), byref(numChannels)),\
            "GetNumOfInputChannels")
    print("Device has %i input channels." % numChannels.value)
    tolog=tolog+"\nDevice has %i input channels." % numChannels.value
    
    print("\nCalibrating...")
    tolog=tolog+"\nCalibrating..."
    tryfunc(hhlib.HH_Calibrate(ct.c_int(dev[0])), "Calibrate")
    tryfunc(hhlib.HH_SetSyncDiv(ct.c_int(dev[0]), ct.c_int(HHsettings["syncDivider"])), "SetSyncDiv")
    
    tryfunc(
        hhlib.HH_SetSyncCFD(ct.c_int(dev[0]), ct.c_int(HHsettings["syncCFDLevel"]),
                            ct.c_int(HHsettings["syncCFDZeroCross"])),\
        "SetSyncCFD"
    )
    tryfunc(hhlib.HH_SetSyncChannelOffset(ct.c_int(dev[0]), ct.c_int(HHsettings["syncChannelOffset"])),\
            "SetSyncChannelOffset")
    tryfunc(hhlib.HH_SetMarkerEnable(ct.c_int(dev[0]), ct.c_int(1),ct.c_int(1),ct.c_int(1),ct.c_int(1)),\
            "SetMarkerEnable")
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
        tolog=tolog+"\nResolution is %1.1lfps" % resolution.value
    
    # Note: after Init or SetSyncDiv you must allow >100 ms for valid  count rate readings
    time.sleep(0.2)
    
    tryfunc(hhlib.HH_GetSyncRate(ct.c_int(dev[0]), byref(syncRate)), "GetSyncRate")
    HHsettings["syncRate"]=syncRate.value
    print("\nSyncrate=%1d/s" % syncRate.value)
    tolog=tolog+"\nSyncrate=%1d/s" % syncRate.value
    for i in range(0, numChannels.value):
        tryfunc(hhlib.HH_GetCountRate(ct.c_int(dev[0]), ct.c_int(i), byref(countRate)),\
                "GetCountRate")
        print("Countrate[%1d]=%1d/s" % (i, countRate.value))
        tolog=tolog+"\nCountrate[%1d]=%1d/s" % (i, countRate.value)
        
    # Check for warnings before starting
    tryfunc(hhlib.HH_GetWarnings(ct.c_int(dev[0]), byref(warnings)), "GetWarnings")
    if warnings.value != 0:
        hhlib.HH_GetWarningsText(ct.c_int(dev[0]), warningstext, warnings)
        print("\n\n%s" % warningstext.value.decode("utf-8")) 
        tolog=tolog+"\n\n%s" % warningstext.value.decode("utf-8")
    
    progress=0
#    sys.stdout.write("\nProgress:%9u" % progress)
#    sys.stdout.flush()
    running=True
    tryfunc(hhlib.HH_StartMeas(ct.c_int(dev[0]), ct.c_int(HHsettings["tacq"])), "StartMeas")
    HHstatus.value=1
    conn.send(tolog)
    counts=0;
    pixel=-1;
    while running==True:
        if currpix.value > pixel:   #Check if stage moved to next pixel
            conn.send(counts)   #send counts of last pixel to scan loop
            pixel=currpix.value#update pixel value
            counts=0    #reset counts
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
#            sys.stdout.write("\rProgress:%9u" % progress)
#            sys.stdout.flush()
        else:
            tryfunc(hhlib.HH_CTCStatus(ct.c_int(dev[0]), byref(ctcstatus)),\
                    "CTCStatus")
            if ctcstatus.value > 0:     #I should check if I can avoid this and just kill the loop when the stage is done
                print("\nAcq time over")
                stoptttr()
                running=False
        if MeasurementRunning.value==0:
            print("\scan time over")
            stoptttr()
            conn.send(counts)   #Send counts for the last pixel
            running=False
    print("terminateing hydraharp")
    HHsettings["overallCounts"]=progress
    save_obj(HHsettings,SampleName+"_HHsettings",OutputFolder)
    closeDevices()
    outputfile.close()
    print("HHloop done")
    conn.close()
    
#%% Start main
if __name__ == '__main__':
    ########################################################
    # Define Name of Output Folder Here
    # This folder needs to exist
    ########################################################
    # Warnings: Make sure that the lights are off
    
    OutputFolder = "C:/Users/omel-acton/Desktop/Gabriel_Local_EMCCD/22-11-18/"
    SampleName = "AP-3-108_C1x10-6_SingleDot_1x10+6Hz_405exct_220nW_"

    #%%Initialize Stage
    mcldll = cdll.madlib
    logfilename='logfile.txt'
    logfile=open(logfilename,'w+')
    # For all remaining calls do : mcldll.FunctionName(function parameters)
    handle = mcldll.MCL_InitHandle()
    if handle == 0:
        print("Error: Handle not intialized correctly\nExiting")
        logfile.write("\nError: Handle not intialized correctly\nExiting")
        logfile.close()
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
        logfile.write("\nError: NanoDrive could not get productInformation. Error Code:", err, "\nExiting")
        mcldll.MCL_ReleaseHandle(handle) # be sure to release handle anytime before returning
        sys.exit()
    else:
        print("Information about the NanoDrive:")
        logfile.write("\nInformation about the NanoDrive:")
        sleep(2)
        print("axis bitmap:", pi.axis_bitmap)
        logfile.write("\naxis bitmap:"+ str(pi.axis_bitmap))
        print("ADC resolution:", pi.ADC_resolution)
        logfile.write("\nADC resolution:"+ str( pi.ADC_resolution))
        print("DAC resolution:", pi.DAC_resolution)
        logfile.write("\nDAC resolution:"+ str( pi.DAC_resolution))
        print("Product ID:", pi.Product_id)
        logfile.write("\nProduct ID:"+ str( pi.Product_id))
        print("Firmware Version:", pi.FirmwareVersion)
        logfile.write("\nFirmware Version:"+ str( pi.FirmwareVersion))
        print("Firmware Profile:", pi.FirmwareProfile)
        logfile.write("\nFirmware Profile:"+ str( pi.FirmwareProfile))
    cal = mcldll.MCL_GetCalibration
    cal.restype = c_double
    readpos = mcldll.MCL_SingleReadN
    readpos.restype = c_double
    xaxis = c_uint(1)
    yaxis = c_uint(2)
    zaxis = c_uint(3)
    mcldll.MCL_IssBindClockToAxis(c_uint(1), c_uint(2), xaxis, handle)
    xcalibration = mcldll.MCL_GetCalibration(xaxis, handle)
    ycalibration = mcldll.MCL_GetCalibration(yaxis, handle)
    zcalibration = mcldll.MCL_GetCalibration(zaxis, handle)
    print("axis calibration range is ",xcalibration,ycalibration,zcalibration)
    logfile.write("\naxis calibration range is x:"+str(xcalibration) +" y: "+str(ycalibration)+" z: "+str(zcalibration))
    inpos=MCLReadpos(handle)
    print("Current position is ",inpos[0],inpos[1],inpos[2])
    logfile.write("\nCurrent position is x: "+str(inpos[0])+"y: "+str(inpos[1])+ "z: " +str(inpos[2]))
    print("Initialization Successfull")
    logfile.write("\nInitialization Successfull")
    
    
    #%%Scan parameters
    MCLdata={}
    ScanCenter=inpos    # For now start the scan at initial position, so I don't crash into anything
    #ScanCenter[0]=169.044
    #ScanCenter[1]=182.739
    ScanRange=40     #in micrometer
    Points=200      #Points to scan per axis
    dwell=10
    OutputName=SampleName+"x_"+str(round(ScanCenter[0]))+"_y"+str(round(ScanCenter[1]))+"_z"+str(round(ScanCenter[2]))+"_range"+str(ScanRange)+"_pts"+str(Points)
    xtargetvals=list(ScanCenter[0]+((j+0.5-Points/2)/Points)*ScanRange for j in range(Points))
    ytargetvals=list(ScanCenter[1]+((j+0.5-Points/2)/Points)*ScanRange for j in range(Points))
    xpos=np.zeros((Points,Points))
    xpos2=np.zeros((Points,Points))
    ypos=np.zeros((Points,Points))
    test1=np.zeros(Points)
    test2=np.zeros(Points)
    print("Expected runtime: ",Points*Points*dwell/1000)
    logfile.write("\nExpected runtime: "+str(Points*Points*dwell/1000))
    print("Resolution :",ScanRange/Points," um per pixel")
    logfile.write("\nResolution :" + str(ScanRange/Points) + " um per pixel")
    MCLdata["xtargetvals"]=xtargetvals
    MCLdata["ytargetvals"]=ytargetvals
    MCLdata["dwell"]=dwell
    MCLdata["ScanCenter"]=ScanCenter
    MCLdata["Points"]=Points
    MCLdata["ScanRange"]=ScanRange
#%%Define while loop (reads data from Hydraharp and plots)
#Start data acquisition loop
    currpix=Value('i',-1);
    HHstatus=Value('i',0);  #0 - not initialized, 1 - initialized and ready, 2 - error occured and terminated
    parent_conn, child_conn=Pipe()
    MeasurementRunning=Value('d',1);
    proc1 = Process(target=HHloop,args=(child_conn,MeasurementRunning,HHstatus,currpix,OutputFolder,OutputName))
    proc1.start()
    lastval=0
#Run scan loop
    cp=np.zeros(Points*Points)  #Initialize array to save the counts for all the pixels, this could also be a matrix if I scan along 2 dimensions
    xposvec=np.zeros(Points*Points)
    yposvec=np.zeros(Points*Points)
    
    timesarray=dwell/1000*np.ones(Points*Points)
    waittime=0
    while HHstatus.value!=1:
        if HHstatus.value==0:
            sys.stdout.write("\rWaiting for Hydraharp, " + str(waittime) + " s")
            waittime=waittime+0.5
            sys.stdout.flush()
            sleep(0.5)
        if HHstatus.value==2:
            print('\nLost connection to HydraHarp, terminating')
            mcldll.MCL_SingleWriteN(c_double(ScanCenter[1]), yaxis, handle)
            mcldll.MCL_SingleWriteN(c_double(ScanCenter[0]), xaxis, handle)
            mcldll.MCL_ReleaseHandle(handle)
            log=parent_conn.recv()
            print(log)
            logfile.write(log)
            logfile.close()
            proc1.join()
            cv2.destroyAllWindows()
            sys.exit()
    
    HHlog=parent_conn.recv()
    logfile.write(HHlog)
    print(HHlog)
    print("Hydraharp succesfully connected")
    logfile.write("\nHydraharp succesfully connected")
    
            
            
    currpix.value=-1
    countsmatrix=np.zeros((Points,Points))
    xposmatrix=np.zeros((Points,Points))
    yposmatrix=np.zeros((Points,Points))
    cv2.namedWindow('img',cv2.WINDOW_NORMAL)
    cv2.resizeWindow("img", 400, 400) 
    prevtime=time.time()
    print('Starting Scan')
    logfile.write("\nStarting Scan")
    for k in range(Points):
        #move to ytargetvals[k]
        mcldll.MCL_SingleWriteN(c_double(ytargetvals[k]), yaxis, handle)
#        print("arrived here")
        #mcldll.MCL_SingleWriteN(c_double(xtargetvals[0]+((-0.5-Points/2)/Points)*ScanRange), xaxis, handle)
        mcldll.MCL_SingleWriteN(c_double(xtargetvals[0]), xaxis, handle)
        countsmatrix=np.reshape(cp,(Points,Points))
        timesmatrix=np.reshape(timesarray,(Points,Points))
        cv2.imshow("img", np.float32((countsmatrix/timesmatrix)/(np.amax((countsmatrix/timesmatrix))+1)))
        cv2.waitKey(1)
        sleep(30/1000)
#        print(k/Points*100,"%")
        mcldll.MCL_LineClock(handle)
        sys.stdout.write("\rrow "+str(k)+" of "+str(Points)+ "; Maxcps: "+ str(np.max(countsmatrix/timesmatrix))+"; Mincps: "+ str(np.min(countsmatrix/timesmatrix)))
        sys.stdout.flush()
        #wait a bit
        for l in range(Points):
            #move to xtargetvals[l]
            #if k % 2==0:
            mcldll.MCL_SingleWriteN(c_double(xtargetvals[l]), xaxis, handle)
            test1[l]=xtargetvals[l]
            #else:
                #mcldll.MCL_SingleWriteN(c_double(xtargetvals[Points-l-1]), xaxis, handle)
                #test2[l]=xtargetvals[Points-l-1]
            sleep(dwell/1000)  
            xposcurr=mcldll.MCL_SingleReadN(xaxis, handle)
            yposcurr=mcldll.MCL_SingleReadN(yaxis, handle)
            xposmatrix[l,k]=xposcurr
            yposmatrix[l,k]=yposcurr
            currtime=time.time()
            duration=currtime-prevtime
            prevtime=currtime
            currpix.value += 1
            lastval=parent_conn.recv()
#            print(lastval)
            xposvec[currpix.value]=xposcurr
            yposvec[currpix.value]=yposcurr
            if currpix.value > 0:
                cp[currpix.value-1] = lastval
                timesarray[currpix.value-1] = duration
            if HHstatus.value==2:
                mcldll.MCL_SingleWriteN(c_double(ScanCenter[1]), yaxis, handle)
                mcldll.MCL_SingleWriteN(c_double(ScanCenter[0]), xaxis, handle)
                mcldll.MCL_ReleaseHandle(handle)
                proc1.join()
                cv2.destroyAllWindows()
                print('\nLost connection to HydraHarp')
                logfile.write('\nLost connection to HydraHarp')
                
    MeasurementRunning.value=0    #Sends signal to Hydraharp loop to terminate
    cp[currpix.value]=parent_conn.recv() #receive last pixel
    proc1.join()    #this maybe avoids some crashes?? I think this waits for the HH process to terminate
    cv2.destroyAllWindows()
    MCLdata["xposvec"]=xposvec
    MCLdata["yposvec"]=yposvec
    MCLdata["xposmatrix"]=xposmatrix
    MCLdata["yposmatrix"]=yposmatrix
    save_obj(MCLdata, OutputName+"_MCLdata", OutputFolder )
    #print("Counts Array",cp)
    #print("Overall runtime",stopmain-startmain)
#    plt.plot(cp)
#    plt.show()
    print("\nData Saved")
    logfile.write("\nData Saved")
    logfile.close()
    mcldll.MCL_SingleWriteN(c_double(ScanCenter[1]), yaxis, handle)
    mcldll.MCL_SingleWriteN(c_double(ScanCenter[0]), xaxis, handle)
    mcldll.MCL_ReleaseHandle(handle)

    print("\nProgram completed")
    print("\nProgram completed")
    countsmatrix=np.reshape(cp,(Points,Points))
    xvec=np.mean(MCLdata["xposmatrix"],1)
    yvec=np.mean(MCLdata["yposmatrix"],0)
    plt.pcolor(xvec,yvec,countsmatrix/timesmatrix)
    plt.xlabel('x (um)')
    plt.ylabel('y (um)')
    plt.colorbar()
#    plt.show()
#    plt.imshow(countsmatrix/timesmatrix)
    plt.show()
        #proc1.terminate()#Terminates data acquisition. This works but is a bit uggly and might corrupt files down the line. Maybe there is a better way