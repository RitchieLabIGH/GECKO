'''
Compare score for all individuals of GA test vs outer
@arg1 : result folder path, this folder willl contain AllScoreByGeneration.csv and AllOUTTERScoreByGeneration.csv

'''

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
mpl.rcParams['axes.linewidth'] = 0.1

if len(sys.argv) > 1:
    fold = sys.argv[1]
'''
Load files
'''
outterfile=np.genfromtxt(fold+"/AllOUTTERScoreByGeneration.csv",delimiter=',')
scorefile=np.genfromtxt(fold+"/AllScoreByGeneration.csv",delimiter=',')

###
# plot test and outer scores side to side

fig = plt.figure(figsize=(13,10))
ax = fig.add_subplot(121)

plt.imshow(scorefile, vmin=0, vmax=1, cmap='jet', aspect='auto', interpolation='nearest') #cmap='hot'
#plt.colorbar()
ax.set_ylabel("Generations")
ax.set_xlabel("Individus")
ax.set_title("Score test of all individu along GA")

ax = fig.add_subplot(122)
plt.imshow(outterfile, vmin=0, vmax=1, cmap='jet', aspect='auto', interpolation='nearest') #cmap='hot'
plt.colorbar()

ax.set_xlabel("Individus")
ax.set_title("Score outter of all individu along GA")

plt.savefig(fold + "/scores_heatmapcompare.png", dpi=150)

#plt.show()
plt.close()

###
# plot difference between test and outer scores

fig = plt.figure(figsize=(13,10))
ax = fig.add_subplot(111)
diffsc=scorefile-outterfile
plt.imshow(diffsc, vmin=0, vmax=np.max(np.max(diffsc)), cmap='jet', aspect='auto', interpolation='nearest') #cmap='hot',vmax=1,
plt.colorbar()
ax.set_ylabel("Generations")
ax.set_xlabel("Individus")
ax.set_title("Score test-outter of all individu along GA")
plt.savefig(fold + "/scores_heatmapdiff.png", dpi=150)


# sorted outer
# fig = plt.figure(figsize=(13,10))
# ax = fig.add_subplot(111)
# outterfilesorted=np.array(outterfile)
# outterfilesorted.sort(axis=1)
# plt.imshow(outterfilesorted, vmin=0,vmax=1,  cmap='jet', aspect='auto', interpolation='nearest') #cmap='hot',
# plt.colorbar()
# ax.set_ylabel("Generations")
# ax.set_xlabel("Individus")
# ax.set_title("Score outter order of all individu along GA")
# plt.savefig(fold + "/scores_outterorder_heatmap.png", dpi=150)
# plt.close()
'''
Experimental figures of correlation between test score and outer
'''

###
# plot correlation between test and outer scores accross generations
corr=np.zeros((scorefile.shape[0],1))
for i in range(scorefile.shape[0]):
    corr[scorefile.shape[0]-i-1]=np.corrcoef(scorefile[i,:], outterfile[i,:])[0, 1]

fig = plt.figure(figsize=(13,10))
ax = fig.add_subplot(111)
plt.plot(corr)
ax.set_xlabel("Generations")
ax.set_ylabel("correlation")
plt.savefig(fold + "/scores_outter_corr.png", dpi=150)
plt.close()

###
# plot correlation between test and outer scores accross generations for the quarter of individuals with the lower test score of each generation
corr=np.zeros((scorefile.shape[0],1))
quarter=(scorefile.shape[1]//4)
for i in range(scorefile.shape[0]):
    corr[scorefile.shape[0]-i-1]=np.corrcoef(scorefile[i,-quarter:-1], outterfile[i,-quarter:-1])[0, 1]

fig = plt.figure(figsize=(13,10))
ax = fig.add_subplot(111)
plt.plot(corr)
ax.set_xlabel("Generations")
ax.set_ylabel("correlation")
plt.savefig(fold + "/scores_outter_corrlastquarter.png", dpi=150)
plt.close()

###
# plot correlation between test and outer scores accross generations for the quarter of individuals with the best test score of each generation
corr=np.zeros((scorefile.shape[0],1))
quarter=(scorefile.shape[1]//4)
for i in range(scorefile.shape[0]):
    corr[scorefile.shape[0]-i-1]=np.corrcoef(scorefile[i,1:quarter], outterfile[i,1:quarter])[0, 1]

fig = plt.figure(figsize=(13,10))
ax = fig.add_subplot(111)
plt.plot(corr)
ax.set_xlabel("Generations")
ax.set_ylabel("correlation")
plt.savefig(fold + "/scores_outter_corrfirstquarter.png", dpi=150)
plt.close()
