# -*- coding: utf-8 -*-
"""
Created on Mon Oct 1 9 17:09:22 2020

@author: jsanten based on Ario's script
Settings for Andor:
Acquisition Mode: Kinetic
Triggering: External
Exopsure Time: varies
Number of Accumulations: varies (Important: each image needs to be triggered individually!)
Kinetic Series Length: varies
"""

from ctypes import *  # C is needed to use DLL. madlib.dll has to be in same folder as python.exe
import numpy as np
from matplotlib import pyplot as py

import timeit
from time import sleep
import nidaqmx
import time

# Define the same exposure time as in Andor
exposuretimeandor = 0.00010  # in seconds


# intialize camera (connect NI box)
tacq = exposuretimeandor;  # this is set at ANDOR
T = [True]  # True starts it
F = [False]  # True ends it

# %% Scan
# Here a parallel loop would be sweet where scanning and analyzing and plotting are independent from each other
count = 0;
# Scan loop
start = timeit.default_timer()

with nidaqmx.Task() as task:
    task.do_channels.add_do_chan('Dev1/port0/line2')

    for i in range(15):
        # take a spectrum
        task.write(T)
        time.sleep(tacq + 0.1)
        task.write(F)
        sleep(0.5)



stop = timeit.default_timer()
print('RunTime: ', stop - start)

print("Program completed without any errors")