import numpy as np
import pandas as pd
from scipy.stats import rankdata
import csv
import multiprocessing as mp
import math
import time
import os,sys, errno
import datetime
import json
import glob
import matplotlib as mpl
if os.environ.get('DISPLAY','') == '':
    print('no display found. Using non-interactive Agg backend')
    mpl.use('Agg')
else:
    try:
        import tkinter
    except ImportError:
        mpl.use('Agg')

import matplotlib.pyplot as plt

if len(sys.argv) > 1:
    file = sys.argv[1] # id application for communication
txtconf=np.genfromtxt(file,delimiter=',')
fig = plt.figure(figsize=(13,10))
ax = fig.add_subplot(111)

plt.imshow(txtconf, vmin=0, vmax=1, cmap='jet', aspect='auto', interpolation='nearest') #cmap='hot'
plt.colorbar()

ax.set_ylabel("Generations")
ax.set_xlabel("Individus")
ax.set_title("Score test of all individu along AG")

plt.savefig(file + "_heatmap.png", dpi=150)

#plt.show()
plt.close()
