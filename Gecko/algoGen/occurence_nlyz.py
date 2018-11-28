# -*-coding:Latin-1 -*

from __future__ import print_function
from __future__ import unicode_literals
from __future__ import division
from __future__ import absolute_import
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
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from collections import Counter
from functools import reduce




## Cette classe intègre les fichiers json émis par les analyses des algos génétiques
## pour proposer ensuite des graphiques
class Occurrence_nlyz:
    def __init__(self,matKmerFile):
        self.matKmerFile=matKmerFile
        self.initial_data=[]
        self.listPath = []
        self.data =[] # pd.DataFrame()
        self.rawdata = pd.DataFrame()
        self.sumscore=0
        self.saveimg=True
        self.showimg=False
        self.countpop=0
        self.figdir=""
        self.nbkmmax=99999999
        self.nbkmmin = 0
        self.nbhnmax = 99999999
        self.nbhnmin = 0

    def updateInPlace(a, b):
        a.update(b)
        return a
    def addData(self, pathdata):
        self.listPath.append(pathdata)

            #self.data=pd.read_csv(pathdata,header=None,index_col=0,names=["count"])
        ttab = pd.read_csv(pathdata, header=None,  names=["id","count"],dtype={'id':str})
        ttab=dict(ttab.values)
        if len(self.data) == 0:
                self.data =  ttab

        else:
            #ttab= pd.read_csv(pathdata,header=None,index_col=0,names=["count"])
            #ttab = np.loadtxt(pathdata, delimiter=',')
           # self.data=[self.data,ttab]
            for i in ttab:
              #  if self.data.iskey>0:
                self.data[i]+=ttab[i]


            '''print(ttab.shape)
            print(self.data.shape)
            self.data = self.data.append(ttab)
            print(self.data.shape)

            grouped=self.data.groupby(["id"])
            dd=grouped.aggregate(np.sum)

            print(self.data.shape)
            #
            '''
        self.countpop+=1

    def addDirectory(self, pathDirectory):
        listfile = glob.glob(pathDirectory+"/*_OccKmers.csv")
        for fi in listfile:
            self.addData(fi)
        self.rawdata=   self.data.copy()
        del self.data["idpop"]
        self.data=self.data.groupby('bestIndices').sum()

    def count(self):
        print("### count")
        fig = plt.figure(figsize=(13, 10))
       # ax1 = fig.add_subplot(211)
       # n, bins, patches = plt.hist(self.data['scmax'], 20, normed=0, facecolor='blue', alpha=0.75)
        ax2 = fig.add_subplot(111)
        n, bins, patches = plt.hist(self.data['nview'], 20, normed=0, facecolor='green', alpha=0.75)
        #y = mlab.normpdf(bins, mu, sigma)
       # l = plt.plot(bins)
        if self.saveimg:
            plt.savefig(self.figdir+"/histredondance", dpi=300)
        if self.showimg:
            plt.show()
        plt.close()

    def getbest_seq(self):
        print("### getbest_seq")

        if len(self.initial_data) == 0:
            self.initial_data = pd.read_csv(self.matKmerFile, sep=",", index_col=None)
        initial_data_x = self.initial_data.copy()
        del initial_data_x['group']

        b = pd.Series.sort_values(self.data['nview'], ascending=False)


        nseq=min(b.size,1000)

        for i in range(0,nseq):
            print(b[i]," : ",b.index[i]," = ",GeneticAlg_ML.idseq2fasta( initial_data_x.columns[int(b.index[i])], 30))

    def plotBestKmerDistrib(self,nbKmerShow,norm):
        print("### plotBestKmerDistrib")

        if len(self.initial_data)==0 :
            self.initial_data = pd.read_csv(self.matKmerFile, sep=",", index_col=None)
        initial_data_x = self.initial_data.copy()
        del initial_data_x['group']

        b = pd.Series.sort_values(self.data['nview'], ascending=False)
        tab = initial_data_x.values

        ############
        initial_data_y = self.initial_data.group
        num_tmp_data_y = np.zeros(len(initial_data_y))
        cat = np.unique(initial_data_y)
        i = 0
        for grptmp in initial_data_y:
            num_tmp_data_y[i] = np.where(cat == grptmp)[0][0] + 1
            i += 1
            idy = num_tmp_data_y.astype(int)
        groupuniq = np.unique(idy)
#######
        fig = plt.figure(figsize=(13, 10))
        ax = fig.add_subplot(111)
        colors = cm.rainbow(np.linspace(0, 1, len(groupuniq)))
        # norm
        if norm==True:
            normtrain = normalize(tab[:,  (b.index[:]).astype(int)], axis=1, copy=False)

        nseq = min(b.size, nbKmerShow) # show only 10 more saw kmer
        for i in range(0,nseq):
            if norm == True:
            # norm

                val = normtrain[:,i]
            else:
            #not norm
                val = tab[:, int(b.index[i])]

            print(b[i], " : ", b.index[i], " = ",
                  GeneticAlg_ML.idseq2fasta(initial_data_x.columns[int(b.index[i])], 30))
            txtgroup = self.initial_data.group

            for idgrp in groupuniq:
                selval=val[np.where(idy == idgrp)[0]]
               # ax.scatter(np.ones(selval.size)*i, selval, color=colors[idgrp-1])
                sc=ax.scatter((np.arange(selval.size)/100) + i+((idgrp-1)*0.5), selval, color=colors[idgrp - 1])

               # ax.set_ylim([0, 1])
               # ax.set_ylabel("best score")
               # ax.set_xlabel("iteration")
               # ax.set_title("score vs iteration with " + str(pHiddenNeuron) + " neurons and " + str(pNKmer) + " Kmer")
               # n, bins, patches = plt.hist(val(np.where(idy == idgrp)[0]), 20, normed=0, facecolor='blue', alpha=0.75)

            #for iu in np.arange(txtgroup.size):
                # print (idgroup[iu],val[iu])
               # print('{0} : {1}'.format(idy[iu], val[iu]))


       # cbar=plt.colorbar(sc)
       # cbar.set_label('Groups')
        if self.saveimg:
            if norm == True:
                plt.savefig(self.figdir+"/comparekmerNorm"+str(nbKmerShow), dpi=300)
            else:
                plt.savefig(self.figdir+"/comparekmerRaw" + str(nbKmerShow), dpi=300)
        if self.showimg:
            plt.show()
        plt.close()


    def plotpopdistib(self, binsize,nbviewmin,normbykm=True):
        print("### plotpopdistib")

        self.data.groupby(1)



        fig = plt.figure()
        ax = fig.add_subplot(111)
        #plt.bar( np.arange(0,self.countpop,binsize),bstab, width = binsize/1.5)
        plt.scatter(np.arange(binsize, self.countpop, binsize), bstab)
        if normbykm != True:
            ax.set_ylabel("nb unique kmer (>=" + str(nbviewmin) + " Nb view/Nb Kmer total)")
        else:
            ax.set_ylabel("nb unique kmer (>=" + str(nbviewmin) + " Nb view)")

        ax.set_xlabel("number of population")

        #ax.set_ylim(0, 1)
#        ax.set_title("score vs number of neurons with " + str(pNKmer) + " Kmers")
        np.savetxt(self.figdir+"/plotpopdistib_norm"+str(normbykm) + str(binsize) + "Kmer"+str(nbviewmin)+"Nbview.csv", bstab, delimiter=",")
        if self.saveimg:
            plt.savefig(self.figdir+"/plotpopdistib_norm"+str(normbykm)  + str(binsize) + "Kmer"+str(nbviewmin)+"Nbview", dpi=300)
            plt.savefig(self.figdir + "/pdf/plotpopdistib_norm"+str(normbykm)  + str(binsize) + "Kmer" + str(nbviewmin) + "Nbview.pdf", dpi=300)
        if self.showimg:
            plt.show()
        plt.close()
    def savefigdir(self,path):
        self.figdir=path
        try:
            os.makedirs(path)
        except OSError as e:
            if e.errno != errno.EEXIST:
                raise
        try:
            os.makedirs(path + "/pdf")
        except OSError as e:
            if e.errno != errno.EEXIST:
                raise

#test = GeneticAlgPlotter("faketab.csv")
#test.addDirectoryJson("resfake")

test = Occurrence_nlyz("InOUtMatrix_onlyInformativeKMers_ML.txt")

test.addDirectory("beautytest")
test.savefigdir("fig_resmokerADA")


#test = GeneticAlgPlotter("InOUtMatrix_onlyInformativeKMers_ML.txt")
#test.addDirectoryJson("res22")
test.plotBestKmerDistrib(2,norm=False)

#test.plotpopdistib(10,2)
#test.plotpopdistib(10,1)

#test.plotpopdistib(10,3)







test.plotBestKmerDistrib(2,norm=False)
test.plotBestKmerDistrib(10,norm=False)

test.plotBestKmerDistrib(2,norm=True)
test.plotBestKmerDistrib(10,norm=True)

test.getbest_seq()




test.count()

#nn=GeneticAlg_ML(initial_data_x, initial_data_y, numberKM, hideNeurons, numberPopulation, numberGenerations, doLog, "slurmdeepbigMoy/fishing", numberCores)
#nn=GeneticAlg_ML()
print(GeneticAlg_ML.idseq2fasta( 545, 30))
