import matplotlib
matplotlib.use('agg')
import pylab

import random
from time import sleep
from datetime import datetime, timedelta

max_sleep_time = 0.05
def accurate_delay(delay):
    ''' Function to provide accurate time delay in millisecond
    '''
    _ = time.perf_counter() + delay/1000
    while time.perf_counter() < _:
        pass
req_sleep = []
act_sleep = []

for samp in range(0,1000):
    sleep_time = random.random() * max_sleep_time
    req_sleep.append(sleep_time * 1000.)

    start = datetime.now()
    sleep(sleep_time)
#    accurate_delay(sleep_time*1000)
    end = datetime.now()

    act_sleep.append((end - start).total_seconds() * 1000.)

pylab.figure()

pylab.plot(req_sleep, act_sleep, 'r+')
pylab.plot([0, max_sleep_time * 1000], [0, max_sleep_time * 1000.], 'k-')
pylab.plot([0, 0.8 * max_sleep_time * 1000], [0, max_sleep_time * 1000.], 'k--')
pylab.xlabel("Requested Sleep Time ($t_r$; ms)")
pylab.ylabel("Actual Sleep Time ($t_s$; ms)")
pylab.xlim(0, max_sleep_time * 1000.)
pylab.ylim(0, max_sleep_time * 1000.)

pylab.annotate(r"$t_s \approx \frac{5}{4} t_r$", xy=(4, 5), xytext=(0.3, 0.7), textcoords='axes fraction',
            arrowprops=dict(arrowstyle="->", connectionstyle="angle3,angleA=-90,angleB=-45"),
            fontsize=14, bbox=dict(fc='w', ec='k', boxstyle='round'),
            ha='center',
            )   

pylab.title("Sleep behavior of Windows 10")
pylab.grid()
pylab.savefig("sleeptest_EMCCD.png")