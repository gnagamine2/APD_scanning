# -*- coding: utf-8 -*-
"""
Created on Thu Jan 10 16:17:52 2019

@author: Robert Keitel
"""
from ctypes import * #C is needed to use DLL. madlib.dll has to be in same folder as python.exe
import numpy as np
from matplotlib import pyplot as py
import timeit
from time import sleep
mcldll = cdll.madlib
#%% Define functions
def MCLReadpos(handle):
    xval=mcldll.MCL_SingleReadN(c_uint(1), handle)
    yval=mcldll.MCL_SingleReadN(c_uint(2), handle)
    zval=mcldll.MCL_SingleReadN(c_uint(3), handle)
    pos=[xval,yval,zval]
    return pos
#%%Initialize Stage
mcldll = cdll.madlib
# For all remaining calls do : mcldll.FunctionName(function parameters)
handle = mcldll.MCL_InitHandle()
if handle == 0:
    print("Error: Handle not intialized correctly\nExiting")

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

#%%Scan parameters
scancenter=inpos
# For now start the scan at initial position so I don't crash into anything
scanrange=20     #in micrometer
points=160      #Points to scan per axis
xtargetvals=list(scancenter[0]+((j+0.5-points/2)/points)*scanrange for j in range(points))
ytargetvals=list(scancenter[1]+((j+0.5-points/2)/points)*scanrange for j in range(points))
xpos=np.zeros((points,points))
xpos2=np.zeros((points,points))
ypos=np.zeros((points,points))
test1=np.zeros(points)
test2=np.zeros(points)
dwell=10
print("Expected runtime: ",points*points*dwell/1000)
print("Resolution :",scanrange/points," um per pixel")
#%% Scan
#Here a parallel loop would be sweet where scanning and analyzing and plotting are independent from each other

#Scan loop
start = timeit.default_timer()
for k in range(points):
    #move to ytargetvals[k]
    mcldll.MCL_SingleWriteN(c_double(ytargetvals[k]), yaxis, handle)
    #mcldll.MCL_SingleWriteN(c_double(xtargetvals[0]+((-0.5-points/2)/points)*scanrange), xaxis, handle)
    mcldll.MCL_SingleWriteN(c_double(xtargetvals[0]), xaxis, handle)
    sleep(30/1000)
    print(k/points*100,"%")
    mcldll.MCL_LineClock(handle)
    #wait a bit
    for l in range(points):
        #move to xtargetvals[l]
        #if k % 2==0:
        mcldll.MCL_SingleWriteN(c_double(xtargetvals[l]), xaxis, handle)
        test1[l]=xtargetvals[l]
        #else:
            #mcldll.MCL_SingleWriteN(c_double(xtargetvals[points-l-1]), xaxis, handle)
            #test2[l]=xtargetvals[points-l-1]
        sleep(dwell/1000)  
        xpos[l,k]=mcldll.MCL_SingleReadN(xaxis, handle)
        ypos[l,k]=mcldll.MCL_SingleReadN(yaxis, handle)
        #sleep(dwell/1000)            
        #wait a bit

        #read vals
stop = timeit.default_timer()

print('RunTime: ', stop - start,"Optimum: ",points*points*dwell/1000)
#py.imshow(xpos)
#py.imshow(ypos)
py.plot(test1-scanrange/points,xpos[:,5])
#py.plot(test2)
#py.plot(xpos[:,6])   
py.plot(test1,test1)   
py.xlabel("Step")
py.ylabel("Position (um)")
#py.legend(["P40D10","P20D40","goal"])
#py.title(dwell," ms")  
#%% Shutdown
mcldll.MCL_ReleaseHandle(handle)

print("Program completed without any errors")