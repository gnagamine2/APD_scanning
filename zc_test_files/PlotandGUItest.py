# -*- coding: utf-8 -*-
"""
Created on Tue May 26 14:12:22 2020

@author: rober
"""
import matplotlib.pyplot as plt
import time
import random
import sys
from tkinter import filedialog
from tkinter import *

#%%
import threading as th
import time
import keyboard

keep_going = True
def key_capture_thread():
    global keep_going
    while keep_going:
        a = keyboard.read_key()
        # print(a)
        if a== "esc":
            keep_going = False
    time.sleep(1)


# def do_stuff():
# th.Thread(target=key_capture_thread, args=(), name='key_capture_thread', daemon=True).start()
# i=0
# while keep_going:
#     print('still going...')
#     time.sleep(1)
#     i=i+1
#     print (i)
# print ("Schleife beendet")


# do_stuff()
# #%%
root = Tk()
root.withdraw()
outfile=filedialog.asksaveasfilename(title='choose .out file to save')

print('Selected file:')
print(outfile)

print('Happy with save destination and proceed? (y/n)')
decision=input()
if decision=='n':
    sys.exit('Terminated by user')








th.Thread(target=key_capture_thread, args=(), name='key_capture_thread', daemon=True).start()
# #
plt.ion()
fig = plt.figure()
ax = plt.subplot(1,1,1)
ax.set_xlabel('Time (s)')
ax.set_ylabel('cps')
ax.set_title('Move stage to maximize counts, press Esc when happy')
# fig=plt.figure()

# # k=0
t=[]
y=[]
ax.plot( t , y , 'k-' , markersize = 10 )
fig.show()
t0 = time.time()
print('Press Esc when happy with counts')
progress=0

# try:
while keep_going:
    progress+=1
    newval=random.randint(1,10)
    # plt.plot(k,newval,'.r')
    # plt.draw()
    # plt.pause(0.0001)
    # time.sleep(0.1)
    t.append( time.time()-t0 )  # add new x data value
    y.append( newval )        # add new y data value
    ax.lines[0].set_data( t,y ) # set plot data
    ax.relim()                  # recompute the data limits
    ax.autoscale_view()         # automatic axis scaling
    fig.canvas.flush_events()   # update the plot and take care of window events (like resizing etc.)
    time.sleep(0.1)
    sys.stdout.write("\rProgress:%9u cps:%9u" % (progress,newval))
    sys.stdout.flush()

# except KeyboardInterrupt:
#     print('\n Interrupted by user')
#     pass
plt.close()
time.sleep(0.5)
print('Starting measurement')
print('Press Esc to terminate measurement')
keep_going=True
th.Thread(target=key_capture_thread, args=(), name='key_capture_thread', daemon=True).start()
plt.ion()
fig = plt.figure()
ax = plt.subplot(1,1,1)
ax.set_xlabel('Time (s)')
ax.set_ylabel('cps')
ax.set_title('Measuring, press Esc to stop')
# fig=plt.figure()

# # k=0
t=[]
y=[]
ax.plot( t , y , 'k-' , markersize = 10 )
fig.show()
t0 = time.time()
# try:
while keep_going:
    progress+=1
    newval=random.randint(1,10)
    # plt.plot(k,newval,'.r')
    # plt.draw()
    # plt.pause(0.0001)
    # time.sleep(0.1)
    t.append( time.time()-t0 )  # add new x data value
    y.append( newval )        # add new y data value
    ax.lines[0].set_data( t,y ) # set plot data
    ax.relim()                  # recompute the data limits
    ax.autoscale_view()         # automatic axis scaling
    fig.canvas.flush_events()   # update the plot and take care of window events (like resizing etc.)
    time.sleep(0.1)
    sys.stdout.write("\rProgress:%9u cps:%9u" % (progress,newval))
    sys.stdout.flush()
    # sys.stdout.write("\rcps:%9u" % progress)
    # sys.stdout.flush()
# except KeyboardInterrupt:
#     print('\n Interrupted by user')
#     pass
print('Worked')
time.sleep(1)