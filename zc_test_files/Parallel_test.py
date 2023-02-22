from multiprocessing import Process, Pipe,Value,Queue
from time import sleep
import time
import numpy as np
from matplotlib import pyplot as plt
from matplotlib import animation
import cv2
import sys
import ctypes as ct
from drawnow import drawnow
#%%Define while loop (reads data from Hydraharp and plots)
def func1(conn,Flag,currpix):
    counts=0;
    pixel=-1;
    hhlib = ct.CDLL("hhlib64.dll")
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
    HHsettings["tacq"] = 10000000 # Measurement time in millisec, you can change this
    HHsettings["syncDivider"] = 2 # you can change this, observe mode! READ MANUAL!; Check in software!!! In software it is 8
    HHsettings["syncCFDZeroCross"] = 10 # you can change this (in mV); Check in software!!! - correct
    HHsettings["syncCFDLevel"] = 200 # you can change this (in mV); Check in software!!! - correct
    HHsettings["syncChannelOffset"] = 80000 # you can change this (in ps, like a cable delay); Check in software!!!
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
    print("\nMode             : %d" % mode)
    tolog=tolog+"\nMode             : %d" % mode
    print("\nBinning           : %d" % HHsettings["binning"])
    tolog=tolog+"\nBinning           : %d" % HHsettings["binning"]
    print("\nOffset            : %d" % HHsettings["offset"])
    tolog=tolog+"Offset            : %d" % HHsettings["offset"]
    print("\nAcquisitionTime   : %d" % HHsettings["tacq"])
    tolog=tolog+"AcquisitionTime   : %d" % HHsettings["tacq"]
    print("\nSyncDivider       : %d" % HHsettings["syncDivider"])
    tolog=tolog+"SyncDivider       : %d" % HHsettings["syncDivider"]
    print("\nSyncCFDZeroCross  : %d" % HHsettings["syncCFDZeroCross"])
    tolog=tolog+"SyncCFDZeroCross  : %d" % HHsettings["syncCFDZeroCross"]
    print("\nSyncCFDLevel      : %d" % HHsettings["syncCFDLevel"])
    tolog=tolog+"SyncCFDLevel      : %d" % HHsettings["syncCFDLevel"]
    print("\nInputCFDZeroCross : %d" % HHsettings["inputCFDZeroCross"])
    tolog=tolog+"InputCFDZeroCross : %d" % HHsettings["inputCFDZeroCross"]
    print("\nInputCFDLevel     : %d" % HHsettings["inputCFDLevel"])
    tolog=tolog+"InputCFDLevel     : %d" % HHsettings["inputCFDLevel"]
    conn.send(tolog)
    Flag.value=1
    while Flag.value==1:
        if currpix.value > pixel:   #Check if stage moved to next pixel
            conn.send(counts)   #send counts of last pixel to scan loop
            pixel=currpix.value#update pixel value
            counts=0    #reset counts
        counts=counts+1 #Here would be the actual output from the hydraharp
#        sleep(0.005)  #This probably also wouldn't be here but I would write to a file
    conn.send(counts)   #Send counts for the last pixel

def accurate_delay(delay):
    ''' Function to provide accurate time delay in millisecond
    '''
    _ = time.perf_counter() + delay/1000
    while time.perf_counter() < _:
        pass
    
def makeFig():
    plt.imshow(image)

points=30    #Points to scan
startmain=time.time()
dwelltime=10
if __name__ == '__main__':#This is somehow necessary for reasons I don't understand
#Start data acquisition loop
    currpix=Value('i',-1);
    parent_conn, child_conn=Pipe()
    Flag=Value('i',0);
    q=Queue()
    outputname='Paralleltest.txt'
    logfile=open(outputname,'w+')
    logfile.write('Stage Initialized')
    proc1 = Process(target=func1,args=(child_conn,Flag,currpix))
    proc1.start()
    waittime=0
    while Flag.value==0:
        sys.stdout.write("\rWaiting for Hydraharp, " + str(waittime) + " s")
        waittime=waittime+0.5
        sys.stdout.flush()
        sleep(0.5)
    logfile.write(parent_conn.recv())
    #proc1.join()

#Run scan loop
    cp=np.zeros(points*points)  #Initialize array to save the counts for all the pixels, this could also be a matrix if I scan along 2 dimensions
    sleep(0.2)
    countsmatrix=np.zeros((points,points))
    # cv2.namedWindow('img',cv2.WINDOW_NORMAL)
    # cv2.resizeWindow("img", 400, 400) 

    start=time.time()
    for j in range(0,points):
        sys.stdout.write("\rrow "+str(j)+" of "+str(points)+ "; Maxcps: "+ str(np.max(countsmatrix))+"; Mincps: "+ str(np.min(countsmatrix)))
        sys.stdout.flush()
#        print("row",j)
        countsmatrix=np.reshape(cp,(points,points))
        image=countsmatrix/(np.amax(countsmatrix)+1)
        drawnow(makeFig)
        # cv2.imshow("img", np.float32(countsmatrix/(np.amax(countsmatrix)+1)))
        # cv2.waitKey(1)
        for k in range(0,points):
            currpix.value += 1;
            lastval=parent_conn.recv()
            if currpix.value > 0:
#                cp[currpix.value-1] = 1
                cp[currpix.value-1] = lastval
        #print("func2 iterating",j)
        #print("Counts array",cp)
#            accurate_delay(dwelltime)
            sleep(dwelltime/1000.0) #This seems to be not so precise at this level 
    Flag.value=0    #Sends signal to while loop to terminate
    cp[currpix.value]=parent_conn.recv() #receive last pixel
    # cv2.destroyAllWindows()
    stop=time.time()
    print(stop-start, "best case: ", points*points*dwelltime/1000)
    countsmatrix=np.reshape(cp,(points,points))
    logfile.close()
    #plt.show()
    #print("Counts Array",cp)
    #print("Overall runtime",stopmain-startmain)
    plt.imshow(countsmatrix)
    plt.show()
    proc1.terminate()#Terminates data acquisition. This works but is a bit uggly and might corrupt files down the line. Maybe there is a better way



