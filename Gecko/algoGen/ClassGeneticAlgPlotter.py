# -*-coding:Latin-1 -*
'''
Main class for GA result analyze
- load results files
- scores history
- extract best individuals
- count occurence of kmer
- clustering of individuals
- compare kmerlist

'''





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
import matplotlib as mpl
import smoothsig
# try:
#     import tkinter
# except ImportError:
#     mpl.use('Agg')

if os.environ.get('DISPLAY','') == '':
    print('no display found. Using non-interactive Agg backend')
    mpl.use('Agg')
else:
    try:
        import tkinter
    except ImportError:
        mpl.use('Agg')

import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.colors import ListedColormap
from sklearn.neural_network import MLPClassifier
from sklearn.preprocessing import normalize
from sklearn.metrics import accuracy_score
from sklearn import preprocessing
import seaborn as sns

#import mpl_toolkits.mplot3d
#from mpl_toolkits.mplot3d import Axes3D
#from mpl_toolkits.mplot3d.axes3d import Axes3D

#import seaborn as sns
from sklearn.decomposition import PCA
import seaborn as sns
from numpy.polynomial import Polynomial



## @package GeneticAlgPlotter
#  This class make analysis of genetics algorithms from result file,
#  and generate graph, output in format png,pdf, and json
#
#

class GeneticAlgPlotter:
    def __init__(self):
        self.listPath = []
        self.data = pd.DataFrame()
        self.saveimg=True
        self.showimg=False
        self.figdir =""
        self.initial_data=[]
        self.matKmerFile = ""
        self.nbBasePerKmer = 30

    ######______________________________________________________________________________
    # transform sequence ID to real nucleotidic sequence
    # @arg1 : sequence id
    # output : real nucleotidic sequence string
    def idseq2fasta(self, idseq):
        realid = np.uint64(idseq)
        result = ""
        k = self.nbBasePerKmer - 1
        letters = ["A", "C", "T", "G"]
        while k >= 0:
            kmin = np.power(4, k, dtype='int64')
            if kmin <= realid:
                alpha = 3
                statealpha = True
                while statealpha:
                    if realid >= np.uint64(alpha * kmin):
                        realid -= np.uint64(alpha * kmin)
                        statealpha = False
                        result += letters[alpha]
                    else:
                        alpha -= 1
            else:
                result += letters[0]
            k -= 1
        return result
    ######______________________________________________________________________________
    # Add result from file to self.data
    # @arg1 : result file path
    def addDataJson(self, pathdata):
        self.listPath.append(pathdata)
        print("pathdata:"+pathdata)
        tab_c = json.load(open(pathdata,"r"))
        df_c =  pd.DataFrame.from_dict( tab_c['run'] , orient='index')
        df_c['run'] = df_c.index
        filename, file_extension = os.path.splitext(pathdata)
        df_c['kmerlist']=filename+".json_OccKmers.csv"
        df_c=df_c[df_c.scoremax != 'NA']
        colToChange = ['run', 'scoremax', 'nbElem', 'scoremoy10', 'scoremin']
        df_c[colToChange] = df_c[colToChange].apply(pd.to_numeric)
        if('scoreHDin' in df_c.columns):
            df_c[['scoreHDin']] = df_c[['scoreHDin']].apply(pd.to_numeric)
       # print(df_c["run"].values.min())
        df_c["run"] = df_c["run"].values - df_c["run"].values.min() + 1


        self.data = pd.concat([self.data, df_c], ignore_index=True)
    ######______________________________________________________________________________
    # Add result directory to self.data, looking for resultAG.json files
    # @arg1 : directory path
    # @arg2 : file patern
    def addDirectoryJson(self, pathDirectory,files=""):
        listjson = glob.glob(pathDirectory+files+"resultAG.json")
        for fi in listjson:
            self.addDataJson(fi)

    ######______________________________________________________________________________
    # Compare two individual kmer composition
    # @arg1 : individual 1
    # @arg2 : individual 2
    # return ratio of intersection ( based on the bigest individuals)
    def comparekmlist(self, orga1,orga2):
        km1 = orga1.split(",")
        #km1.sort()
        km2 = orga2.split(",")
        #km2.sort()
        inters=list(set(km1).intersection(set(km2)))
        return len(inters)/np.max([len(km1),len(km2)])
    '''
    extract best individual with given interval and stop genration ( call by nlyzHistoryWinnerPCA.py)
    '''
    ######______________________________________________________________________________
    # Extract best individual with given interval and stop genration ( call by nlyzHistoryWinnerPCA.py)
    #@arg1 : number of generations between each mainfold analyze
    #@arg2 : limit of generation to stop the historic
    #@arg3 : output filename root
    # create csv count file for extraction : *forextractkm.count
    # create fastq file : *.fastq
    # create csv score file *_score.txt
    def extracthistoryBestIndiv(self,interindiv,maxindivgene,filename):


        a = self.data.sort_values(by="run",ascending=False)
        a = a.reset_index(drop=True)
        a["run"] = (a["run"].values - a["run"].values.max()-1) *-1
        # remove identics organisms

        #for iwinner in range(1,maxindivgene,interindiv) : #while(ifound < nIndiv-1) & (i < len(listorgawin)-1):
        a = a[0:maxindivgene:interindiv]


        # with open(self.figdir +"/"+filename,"w") as file:
        #     c = a[:]["bestSequences"]
        #     c = c.to_csv(index=False).replace("\"", "")
        #     file.write(c)
        with open(self.figdir +"/"+filename+"forextractkm.count","w") as file:
            c = a[:]["bestSequences"]
            c = c.to_csv(index=False).replace("\"", "").replace(",","\n").replace("\n",",1\n")
            file.write(c)
        with open(self.figdir +"/"+filename+".fastq","w") as file:
            c = a[:]["bestSequences"]
            for i in range(0,a.shape[0]):
                j=0
                for kmerid in c.iloc[i].split(','):
                    file.write( ">indiv" +str( i) +"km"+str(j)+"\n")
                    file.write( self.idseq2fasta(int(kmerid)) + "\n")
                    j+=1

        with open(self.figdir +"/"+filename+"_score.txt", "w") as file:
            c = a[0:a.shape[0]][["scoreHDin","scoremax","run","nbElem","kmerlist"]]
            c = c.to_csv(index=False) #.replace("\"", "")
            file.write(c)
    #####
    # plot score, sort and extract best individuals winner accross generation
    # @arg1 : Number of individuals to extract
    # @arg2 : outputs filenames root
    # create csv count file for extraction : *forextractkm.count
    # create fastq file : *.fastq
    # create csv score file *_score.txt
    # create plot scores accross generation winnerScores.png and /pdf/winnerScores.pdf
    # generate json file for html interface scratch.html ( compatible Highcharts)
    def generateTableOfBestIndiv(self,nIndiv,filename):

        self.data["ArtiMeanscore"]=((self.data["scoremax"]*2)+self.data["scoreHDin"])/3
        listorgawin=self.data.sort_values(by="ArtiMeanscore",ascending=False)
        listorgawin = listorgawin.reset_index(drop=True)
        # remove identics organisms
        a=pd.DataFrame()
        a=a.append(listorgawin.ix[0])
        a['nbview']=1
        ifound=0
        i=1
        while(ifound < nIndiv-1) & (i < len(listorgawin)-1):
            iwinner=0
            found=0
            while (iwinner<=ifound) & (found < ifound+1):
                res=self.comparekmlist(a.ix[iwinner,"bestIndices"],listorgawin.ix[i,"bestIndices"])

                if res == 1:
                    a.ix[iwinner, "nbview"] += 1
                else :
                    found+=1

                iwinner += 1

            if found == ifound+1:
                a = a.append(listorgawin.ix[i], True)
                a.loc[[a.shape[0] - 1], ['nbview']] = 1
                ifound += 1

            i += 1

        #a = self.data.sort_values(by="scoremax", ascending=False)
        with open(filename,"w") as file:
            c = a[0:nIndiv]["bestSequences"]
            c = c.to_csv(index=False).replace("\"", "")
            file.write(c)
        with open(filename+"forextractkm.count","w") as file:
            c = a[0:nIndiv]["bestSequences"]
            c = c.to_csv(index=False).replace("\"", "").replace(",","\n").replace("\n",",1\n")
            file.write(c)
        with open(filename+".fastq","w") as file:
            c = a[0:nIndiv]["bestSequences"]
            for i in range(0,np.min([nIndiv,a.shape[0]])):
                j=0
                for kmerid in c.iloc[i].split(','):
                    file.write( ">indiv" +str( i) +"km"+str(j)+"\n")
                    file.write( self.idseq2fasta(int(kmerid)) + "\n")
                    j+=1

        with open(filename+"_score.txt", "w") as file:
            c = a[0:np.min([nIndiv,a.shape[0]])][["scoreHDin","scoremax","run","nbElem","kmerlist","ArtiMeanscore","nbview"]]
            c = c.to_csv(index=False) #.replace("\"", "")
            file.write(c)

        fig = plt.figure(figsize=(13, 10))
        ax = fig.add_subplot(111)

        listorgawin = listorgawin.sort_values(by=["scoremax","scoreHDin"], ascending=False)
        addname=""
        # for it in range(2):
        #     if it==1:
        #         c=listorgawin[0:5000]
        #         addname = "5000"
        #     else:
        c=listorgawin
        addname = ""

        plt.plot(range(c.shape[0]),c["scoremax"],label="Score")
        plt.plot(range(c.shape[0]),c["scoreHDin"],label="Score outter")
        correlres = np.corrcoef(c["scoremax"] / np.linalg.norm(c["scoremax"]),
                                         c["scoreHDin"] / np.linalg.norm(c["scoreHDin"]))[0, 1]
        #plt.plot(c["ArtiMeanscore"],label="ArtiMeanscore")



        # y = mlab.normpdf(bins, mu, sigma)
        # l = plt.plot(bins)
        ax.legend()
        ax.set_ylabel("Score")
        ax.set_xlabel("organisms")
        ax.set_title("Sort winners scores")
        if self.saveimg:
            plt.savefig(self.figdir + "/winnerScores"+addname,
                        dpi=300)
            plt.savefig(
                self.figdir + "/pdf/winnerScores"+addname+".pdf",
                dpi=300)

        if self.showimg:
            plt.show()
        plt.close()


        scoremaxb = np.flip(c["scoremax"], 0)
        scoreHDinb = np.flip(c["scoreHDin"], 0)

        jsond = ""
        jsond = 'sortwin={'
        jsond += 'title:"Sort winners scores",'
        jsond += 'subtitle:"Correlation test/outter = '+'{:.3f}'.format(correlres)+'",'
        jsond += 'categories:[],\n'
        jsond += 'ylabel:"Score",'
        jsond += "series:["
        jsond += "{data:" + json.dumps(scoreHDinb.tolist()) + ",\n"
        # jsond += "colorIndex:" + json.dumps(icat) + ",\n"
        jsond += 'name:"Outters score"},\n'
        # jsond += 'marker: {  symbol:"circle" } },'
        jsond += "{data:" + json.dumps(scoremaxb.tolist()) + ",\n"
        # jsond += "colorIndex:" + json.dumps(icat) + ",\n"
        jsond += 'name:"Tests score"}\n'
        jsond += "]}"
        with open(self.figdir + "/sortwinnerscore.json", 'w') as f:
            f.write(jsond)


    ####
    # Set and create directories for result files saving
    # @arg1 : path to the root directory to save files
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

    #####
    # Generate graph for score history
    # graph average for 500 points max ( along generations)
    # result: image historyscore.png and pdf/historyscore.pdf
    #         historyscore.json  json file for html interface
    def scorehistorymean(self):
        df = self.data.sort_values(by=['run']).copy()
        fig = plt.figure(figsize=(13, 10))
        #ax = fig.add_subplot(211)

        #sc = ax.plot(df["run"],df["scoremax"])
        #ax = fig.add_subplot(212)
        #sc = ax.plot(df["run"], df["scoreHDin"])
        #sc = ax.plot(df["run"], df["scoremoy10"])
        #sc = ax.plot(df["run"], df["scoremin"])

        a = np.asarray(df["scoreHDin"])
        if len(a)>500:
            xp = np.arange(0, len(a), int(len(a)/500))
        else:
            xp=np.arange(0, len(a))


        scoreHDinb = np.interp(xp, np.arange(0, len(a)), a)
        ax = fig.add_subplot(311)
        sc = ax.plot(xp, scoreHDinb)
        scorehdintitle="Outters scores [" + str(np.min(a)) + "-" + str(np.max(a)) + "], mean= " + str(np.mean(a)) + ", std= " + str(np.std(a))
        ax.set_title(scorehdintitle)


        a = np.asarray(df["scoremax"])
        scoremaxb = np.interp(xp, np.arange(0, len(a)), a)
        ax = fig.add_subplot(312)
        sc = ax.plot(xp, scoremaxb)
        scoremaxtitle="Score max ["+str(np.min(a))+"-"+str(np.max(a))+"], mean= "+str(np.mean(a))+", std= "+str(np.std(a))
        ax.set_title(scoremaxtitle)


        a = np.asarray(df["scoremoy10"])
        b = np.interp(xp, np.arange(0, len(a)), a)
        ax = fig.add_subplot(313)
        sc = ax.plot(xp, b)
        ax.set_title("scoremoy10 [" + str(np.min(a)) + "-" + str(np.max(a)) + "], mean= " + str(np.mean(a)) + ", std= " + str(np.std(a)))

        xp = np.flip((xp - (np.max(xp))) * (-1),0)
        scoreHDinb = np.flip(scoreHDinb,0)
        scoremaxb = np.flip(scoremaxb,0)

        jsond = ""
        jsond = 'winhistoricmean={'
        jsond += 'title:"Winners scores historic",'
        jsond += 'subtitle:"'+scoremaxtitle+'",'
        jsond += 'categories:'+json.dumps(xp.tolist()) + ",\n"
        jsond += 'ylabel:"Score",'
        jsond += "series:["
        jsond += "{data:" + json.dumps(scoreHDinb.tolist()) + ",\n"
        #jsond += "colorIndex:" + json.dumps(icat) + ",\n"
        jsond += 'name:"Outters score"},\n'
        #jsond += 'marker: {  symbol:"circle" } },'
        jsond += "{data:" + json.dumps(scoremaxb.tolist()) + ",\n"
        # jsond += "colorIndex:" + json.dumps(icat) + ",\n"
        jsond += 'name:"Tests score"}\n'
        jsond += "]}"
        with open(self.figdir+"/historyscore.json", 'w') as f:
            f.write(jsond)
        try:
            if self.saveimg:
                plt.savefig(self.figdir+"/historyscore", dpi=300)
                plt.savefig(self.figdir+"/pdf/historyscore.pdf", dpi=300)
            if self.showimg:
                plt.show()
        except:
            print("failed generating "+self.figdir+"/historyscore")
        plt.close()

    #####
    # Generate graph for score history ( no average)
    # this function is not used in basic analyze because it become unreadable with a big number of generation
    # scorehistorymean is prefer
    # result: image historyscore.png and pdf/historyscore.pdf

    def scorehistory(self):
        df = self.data.sort_values(by=['run']).copy()
        fig = plt.figure(figsize=(13, 10))
        #ax = fig.add_subplot(211)

        #sc = ax.plot(df["run"],df["scoremax"])
        #ax = fig.add_subplot(212)
        #sc = ax.plot(df["run"], df["scoreHDin"])
        #sc = ax.plot(df["run"], df["scoremoy10"])
        #sc = ax.plot(df["run"], df["scoremin"])

        a = np.asarray(df["scoreHDin"])
        xp = np.arange(0, len(a), 1)
        b = np.interp(xp, np.arange(0, len(a)), a)
        ax = fig.add_subplot(311)
        sc = ax.plot(xp, b)
        ax.set_title("scoreHDin [" + str(np.min(a)) + "-" + str(np.max(a)) + "], mean= " + str(np.mean(a)) + ", std= " + str(np.std(a)))


        a = np.asarray(df["scoremax"])
        b = np.interp(xp, np.arange(0, len(a)), a)
        ax = fig.add_subplot(312)
        sc = ax.plot(xp, b)
        ax.set_title("scoremax ["+str(np.min(a))+"-"+str(np.max(a))+"], mean= "+str(np.mean(a))+", std= "+str(np.std(a)))


        a = np.asarray(df["scoremoy10"])
        b = np.interp(xp, np.arange(0, len(a)), a)
        ax = fig.add_subplot(313)
        sc = ax.plot(xp, b)
        ax.set_title("scoremoy10 [" + str(np.min(a)) + "-" + str(np.max(a)) + "], mean= " + str(np.mean(a)) + ", std= " + str(np.std(a)))


        try:
            if self.saveimg:
                plt.savefig(self.figdir+"/historyscore", dpi=300)
                plt.savefig(self.figdir+"/pdf/historyscore.pdf", dpi=300)
            if self.showimg:
                plt.show()
        except:
            print("failed generating "+self.figdir+"/historyscore")
        plt.close()

    #####
    # Return a dataframe of kmer winner, with indices, number of view end cumlative score
    # @arg1: minimum outer scorethreshold
    def winnersexploration(self,hdinmin=0):

        df = self.data[(self.data.scoreHDin >= hdinmin )]

        aa = df['bestSequences'].to_csv(sep=str(u','), index=False)
        aa = aa.replace('\"\n\"', ',').replace('\"', '').replace('\n', '')
        if len(aa)>1:
            bestIndices = aa.split(',')
            nbkmer=len(df.iloc[0]['bestSequences'].split(','))
            scoresindiv=np.repeat(df['scoremax'],nbkmer)
            kmercount = pd.DataFrame.from_items([('nview', np.ones(len(bestIndices))),
                                                 ('bestIndices', bestIndices),('scoreCum', scoresindiv)])


            kmercount=kmercount.groupby('bestIndices', as_index=False).sum()
        else:
            kmercount = pd.DataFrame.from_items([('nview', []),
                                         ('bestIndices', []), ('scoreCum', [])])

        return kmercount





    ####
    # display Kmer appearance along winner solution of each generations, and display below the comportement of the most seen kmer
    # This function create kmerOccurencesInWinners.pdf in self.figdir/pdf
    # and OccurencesInWinnersstart.txt contain id of the most seen kmer
    def kmerOccurencesInWinners(self):

        kmerList = self.data.sort_values(by=['run'], ascending=False).copy().bestSequences
        #kmerList = self.data.bestSequences.sort_index(            ascending=False)  # retrieves the best indices in resultAG.json, first generation first

        df_countkmerfinal = self.winnersexploration( hdinmin=0)
        correspDico = {}  # make the correspondance between kmer indice and its occId (rank of appearance)
        correspList = np.empty(df_countkmerfinal.shape[0],dtype=int)  # the exact opposite of correspDico, the list index is occId
        listx = np.empty(int(df_countkmerfinal["nview"].sum()), dtype=int)  # lists the x coordinates for the graph (one x for one kmer in the generation's winner, length-->nbgenerations*nbkmerperindividual)
        listy = np.empty(int(df_countkmerfinal["nview"].sum()), dtype=int)  # lists the y coordinates for the graph (rank of appearance of the kmer in the generation's winner)
        nbMostSeenKmers = 5  # Parameter, can be changed
        mostseenkmers = np.empty([len(kmerList), nbMostSeenKmers], dtype=int)  # list of the most seen kmers at each generation
        nb_occurences = np.zeros(df_countkmerfinal.shape[0], dtype=int)  # counts the number of occurences of each kmer
        cptkm=-1
        newkm=-1

        for i, run in enumerate(kmerList):
            table = run.split(',')
            occId = 0
            for kmer1 in table:
                cptkm+=1
                if kmer1 not in correspDico.keys():
                    newkm+=1
                    correspDico[kmer1] = newkm
                    correspList[newkm] =  kmer1
                occId = correspDico[kmer1]
                nb_occurences[occId] += 1

                listx[cptkm] =  i  # i is the generation's number
                listy[cptkm]= occId  # occId is the rank of appearance of the kmer (its id in occurences)

            if (len(nb_occurences) < nbMostSeenKmers): #when there's not enough kmers yet
                temp=np.array([],dtype=int)
                temp=np.argsort(nb_occurences)[-len(nb_occurences):]
                while (len(temp) < nbMostSeenKmers):
                    temp=np.append(temp,-1)
                mostseenkmers[i] = temp
            else:
                mostseenkmers[i] = np.argsort(nb_occurences)[-nbMostSeenKmers:]  # gets the id of the [nbMostSeenKmers] most seen kmers to this generation

        print("Most seen kmers: ", correspList[mostseenkmers[-1:][0]][::-1])

        #Listing the number of identical kmers between the [nbMostSeenKmers] most seen kmers of each generation and the [nbMostSeenKmers] most seen of the last generation
        mostseensimilarities = np.array([],dtype=int)
        for generation in mostseenkmers:
            nbSimilarity = 0
            for id1 in generation:
                for id2 in mostseenkmers[-1:][0]:
                    if (id1 == id2):
                        nbSimilarity += 1
            mostseensimilarities = np.append(mostseensimilarities, nbSimilarity)

        #Listing the values to make the colorbar of the graphs = number of occurences of each kmer present in listy
        colors = np.array([], dtype=int)
        for y in listy:
            colors = np.append(colors, nb_occurences[y])
        
        #Adding ticks on the most seen kmers
        extraticks=mostseenkmers[-1:][0]
        #Creating the legend with the most seen kmer's sequences
        legend='{} most seen kmers:    \n'.format(nbMostSeenKmers)
        for ocidee in mostseenkmers[-1:][0][::-1]:
            #legend+=str(ocidee)+':'+str(correspList[ocidee])+':'+self.idseq2fasta(correspList[ocidee])+'\n'
            legend += str(ocidee) + ':' + self.idseq2fasta(correspList[ocidee]) + '\n'

        #Creating the subgraph of kmer occurence in each generation's winner
        fig = plt.figure(figsize=(15, 20))
        sns.set_style("darkgrid", {"xtick.major.size": "7.0", "ytick.major.size": "7.0"})
        sns.set_context("paper")
        ax1 = plt.subplot2grid((5, 1), (0, 0), rowspan=4)
        pathCollection = ax1.scatter(listx, listy, s=1, marker=',', c=colors, cmap='jet')
        cbar = plt.colorbar(pathCollection, ticks=None)
        cbar.set_label("Number of selections", labelpad=20)
        ax1.set_ylabel("Kmers, in order of appearance", labelpad=20)
        ax1.set_xlabel("Generations", labelpad=25)
        ax1.set_ylim(0, )
        ax1.set_title(self.figdir[:-4] + "\n\nKmers in each generation's winner") #old matplotversion pad=20,bbox={ 'pad':20}

        ax1.text(0.7,0.84,legend,transform=ax1.transAxes,fontsize=15, horizontalalignment='right',fontname='monospace')
        plt.yticks(list(plt.yticks()[0])+extraticks.tolist())
        for item in ([ax1.title, ax1.xaxis.label, ax1.yaxis.label,cbar.ax.yaxis.label] + ax1.get_xticklabels() + ax1.get_yticklabels() + cbar.ax.get_yticklabels()):
            item.set_fontsize(15)

        # Creating the subgraph of most seen kmers similarities with the final ones per generation
        sns.set_style("darkgrid")
        sns.set_context("paper")
        ax2 = plt.subplot2grid((5, 1), (4, 0), sharex=ax1)
        ax2.plot(mostseensimilarities, color='#000080')
        ax2.set_ylabel("Most seen kmers similarity\nwith final ones", labelpad=30)
        ax2.set_xlabel("Generations", labelpad=30)
        ax2.set_ylim(0, nbMostSeenKmers)
        ax2.set_xlim(0, )
        for item in ([ax2.title, ax2.xaxis.label, ax2.yaxis.label] + ax2.get_xticklabels() + ax2.get_yticklabels()):
            item.set_fontsize(15)
        # Little trick to make the second plot the same size than the first one
        pos = ax1.get_position()
        pos2 = ax2.get_position()
        ax2.set_position([pos.x0, pos2.y0, pos.width, pos2.height])

        plt.savefig(self.figdir + "kmerOccurencesInWinners.pdf", format='pdf', dpi=300)

        #sylvain add

        plt.savefig(self.figdir + "kmerOccurencesInWinners.png", dpi=300)
        plt.close()
        #Extract and print generation where kmer most viewed became leaders
        newbestgene =[]
        with open(self.figdir + "OccurencesInWinnersstart.txt", "w") as file:

            for nbf in range(0, np.max(mostseensimilarities)):
                a=np.argmax(mostseensimilarities > nbf)
                print(nbf+1 , "," , a, file=file)
                newbestgene.append( a)

        return newbestgene

    ####
    # Analyze the kmer exploration during GA for winner or all individuals generated
    # @arg1 : winneronly : if True only the winner kmer will be count, if False all individuals kmer will be count
    # @arg2 : hdinmin : treshold outer score minimum for winneronly = True, no take account for all individuals
    # @arg3 : nbmaxkmer : number max of kmer for kmer listing : text file and fastq
    # output :  - kmcut_selection_{win or all}(hdinmin).png show auto cut of 'most important' kmer
    #           - listkmer{win or all}(hdinmin).fastq    most viewed sequences of all kmer (fastq format)
    #           - countkmer{win or all}(hdinmin).txt    count kmer by id , csv format
    #           - countkmer_tresholdwin(hdinmin).txt    for winner only count kmer with automatic treshod apply
    #           - countkmerwinners_(hdinmin)SortByScore.txt Kmer winner sort by mean test score
    #           - listkmer{win or all}(hdinmin)_SortByScore.fastq   best average score of all kmer squences  (fastq format)
    #           - histokmerallall.png      #sorted distribution of all Kmers count ( number of views)
    #           - histokmerall{win or all}(hdinmin)meanscoresort.png  sorted distribution of winners Kmers test scores average
    #           - explokmer{win or all}(hdinmin).json       json file for html interface contain count of kmer


    def overallHistokmerlist(self,winneronly=False,hdinmin=0,nbmaxkmer=500):


        kmercount = pd.DataFrame()
        df = self.data[( (self.data.run == 1))].copy()
        if winneronly == False:
            for file in df['kmerlist'].values:
                kmercount = kmercount.append(pd.read_csv(file, ",",names=["kmer","nview","scoreCum"]))
                kmercount = kmercount.groupby(by='kmer', as_index=False).sum()
                kmercount["scoremean"] = kmercount["scoreCum"] / kmercount["nview"]


        else :
            kmercount = self.winnersexploration(hdinmin=hdinmin)
            kmercount.columns = ['kmer', 'nview','scoreCum']
            kmercount["scoremean"] = kmercount["scoreCum"] / kmercount["nview"]

        kmercount = kmercount.sort_values(by='nview', ascending=False)

        win = ""
        if winneronly == True:
            win = "win"
        else:
            win="all"


        # ########
        # #cut off for V1 percent of drop derivate of count
        # window = 10
        # #trcount = kmercount["nview"][:np.argmin(kmercount["nview"] >= np.mean(kmercount["nview"]))]
        # # my resample
        # # s_trcount=smoothsig.smooth(trcount, window_len=50, window='hanning')
        # trcount = kmercount["nview"][:int(np.floor(len(kmercount["nview"]) / window) * window)]
        # trcount = trcount.values.reshape(int(len(trcount) / window), window)
        # trcount = np.mean(trcount, axis=1)
        #
        # dtrcount = np.diff(trcount)
        # mindtrcount = np.min(dtrcount)
        # dtrcount = smoothsig.smooth(dtrcount, window_len=150, window='hanning')
        # percentcut = 0.001
        # cutval = (mindtrcount * percentcut)  # cut at 98%decrease of the slope max
        # nbkmersel = np.argmax(np.array(dtrcount > cutval)) * window
        #
        # fig = plt.figure(figsize=(13, 10))
        # ax = fig.add_subplot(211)
        # irang = int(nbkmersel * 0.3)
        #
        # plt.plot(range(0, nbkmersel * 3, 10), dtrcount[:irang])
        # #plt.plot(range(0, nbkmersel * 3, 10),(trcount[:irang] / np.max(trcount[:irang])) * -1 * np.min(dtrcount[:irang]))
        # plt.plot(np.array([0, irang * 10]), np.ones(2) * cutval)
        # plt.plot(np.ones(2) * nbkmersel, [np.min(dtrcount), -np.min(dtrcount)])
        # # plt.plot(range(len(dtrcount)),dtrcount)
        # ax = fig.add_subplot(212)
        # plt.plot(range(0, nbkmersel * 3, 10),trcount[:irang] )
        # plt.plot(np.ones(2) * nbkmersel, [np.min(trcount[:irang]), np.max(trcount[:irang])])
        # if self.saveimg:
        #     plt.savefig(self.figdir + "/histocutsmooth" + win, dpi=300)
        # if self.showimg:
        #     plt.show()
        # plt.close()

        ########
        #cut off V2 3*moy (design for all)
        # percentcut = 3
        # cutval = (np.mean(kmercount["nview"]) * percentcut)
        # nbkmersel = np.argmax(np.array(kmercount["nview"]) < cutval)
        nbkmersel=0
        ########
        # cut off V3 : 0,5 percent drop of bin derivate ( design for winner )
        if (winneronly==True) :
            fig = plt.figure(figsize=(13, 10))
            ax = fig.add_subplot(211)
            tmp = np.array(kmercount["nview"])
            if len(tmp>0):
                [bins, edges] = np.histogram(tmp, np.arange(1, np.max(tmp) + 1, 1))
                edges = edges[1:]
                dedges = edges[:-1]
                dbins = np.diff(bins)
                plt.plot(dedges, dbins)
                ax2 = fig.add_subplot(212)
                plt.plot(edges, bins)
                cutat = np.min(dbins) * 0.005
                ll = np.argmax(np.array(dbins > cutat))
                startkeep = dedges[ll]
                plt.plot(np.ones(2)*startkeep,[np.min(bins), np.max(bins)])
                nbkmersel = np.argmax(np.array(kmercount["nview"]) < startkeep)
            else:
                nbkmersel=0

            if self.saveimg:
                plt.savefig(self.figdir + "/kmcut_selection_" + win+str(hdinmin)+".png", dpi=300)
            if self.showimg:
                plt.show()
            plt.close()
        else:
            if  winneronly==False:    # ########
                #cut off for V1 percent of drop derivate of count
                window = 10
                #trcount = kmercount["nview"][:np.argmin(kmercount["nview"] >= np.mean(kmercount["nview"]))]
                # my resample
                # s_trcount=smoothsig.smooth(trcount, window_len=50, window='hanning')
                trcount = kmercount["nview"][:int(np.floor(len(kmercount["nview"]) / window) * window)]
                trcount = trcount.values.reshape(int(len(trcount) / window), window)
                trcount = np.mean(trcount, axis=1)

                dtrcount = np.diff(trcount)
                mindtrcount = np.min(dtrcount)
                dtrcount = smoothsig.smooth(dtrcount, window_len=10, window='hanning')
                percentcut = 0.01
                cutval = (mindtrcount * percentcut)  # cut at 98%decrease of the slope max
                nbkmersel = np.argmax(np.array(dtrcount > cutval)) * window

                fig = plt.figure(figsize=(13, 10))
                ax = fig.add_subplot(211)
                irang = int(nbkmersel * 0.3)

                if irang<len(dtrcount):
                    plt.plot(range(0, nbkmersel * 3, window), dtrcount[:irang])
                else:
                    plt.plot(range(0, (len(dtrcount) )*10, window),             dtrcount)
                #plt.plot(range(0, nbkmersel * 3, 10),(trcount[:irang] / np.max(trcount[:irang])) * -1 * np.min(dtrcount[:irang]))
                plt.plot(np.array([0, irang * 10]), np.ones(2) * cutval)
                plt.plot(np.ones(2) * nbkmersel, [np.min(dtrcount), -np.min(dtrcount)])
                # plt.plot(range(len(dtrcount)),dtrcount)
                if (nbkmersel != 0) & (np.max(kmercount["nview"]) > nbkmersel * 3):
                    ax.set_xlim(0, nbkmersel * 3)
                ax2 = fig.add_subplot(212)
                if irang < len(dtrcount) & (irang>0):
                    plt.plot(range(0, nbkmersel * 3, window), trcount[:irang])
                else:
                    plt.plot(range(0, (len(trcount)) * 10, window), trcount)
                if irang>0:
                    plt.plot(np.ones(2) * nbkmersel, [np.min(trcount[:irang]), np.max(trcount[:irang])])
                if (nbkmersel != 0) & (np.max(kmercount["nview"]) > nbkmersel * 3):
                    ax2.set_xlim(0, nbkmersel * 3)

                if self.saveimg:
                    plt.savefig(self.figdir + "/kmcut_selection_" + win, dpi=300)
                if self.showimg:
                    plt.show()
                plt.close()


        with open(self.figdir + "feature_distrib" + win + str(hdinmin) + ".txt", "w") as file:
            print("mean_view=",np.mean(kmercount["nview"]),file=file)
            print("max_view=", np.max(kmercount["nview"]), file=file)
            print("std_view=", np.std(kmercount["nview"]), file=file)
            print("nbkmer_selected=",nbkmersel,file=file)

        '''
        fig = plt.figure(figsize=(13, 10))
        ax = fig.add_subplot(111)
        a = np.asarray(kmercount["nview"][0:len(kmercount):1000], dtype="float")
        b = a[1:a.size] - a[0:a.size - 1]
        plt.plot(range(0, a.size - 1), a[1:a.size] - a[0:a.size - 1])

        fig = plt.figure(figsize=(13, 10))
        ax = fig.add_subplot(111)
        a = np.asarray(kmercount["nview"], dtype="float")

        xp = np.arange(0, len(a), 100)
        b = np.interp(xp, np.arange(0, len(a)), a)
        b = a[1:a.size] - a[0:a.size - 1]
        plt.plot(range(0, a.size - 1), b)


        '''

        listbestidkmer=[]
        with open(self.figdir+"listkmer"+win+str(hdinmin)+".fastq", "w") as file:
            with open(self.figdir + "countkmer" + win + str(hdinmin) + ".txt", "w") as filecount:

                for i in range(0, min(nbmaxkmer,kmercount.shape[0])):
                    #kmercount["nview"].values[i]

                   # print(i, " : ", kmercount["nview"].values[i], " = ",idseq2fasta(int(kmercount["kmer"].values[i]), 30))  # ,kmercount["kmer"].values[i]," : ",

                    print(">",win,i,"_ct", kmercount["nview"].values[i], "\n", self.idseq2fasta(int(kmercount["kmer"].values[i])),file=file)
                    print(kmercount["kmer"].values[i], ",", kmercount["nview"].values[i],file=filecount)
        if nbkmersel != 0:
            with open(self.figdir + "countkmer_treshold" + win + str(hdinmin) + ".txt", "w") as filecount:

                for i in range(0, nbkmersel):
                    print(kmercount["kmer"].values[i], ",", kmercount["nview"].values[i], file=filecount)

        '''listbestidkmer.append(kmercount["kmer"].values[i])
        headerkm = pd.read_csv("BEAUTY_RNASeq_HEADER.csv",  sep=",", index_col=False, nrows=0)
        #headerkm = pd.read_csv("faketabcom.txt", sep=",", index_col=False, nrows=0)
        headerkm = headerkm.columns.values
        orig_indices = headerkm.argsort()
        ndx = orig_indices[np.searchsorted(headerkm[orig_indices], listbestidkmer)]

        listbestidkmer=','.join(['%d' % num for num in ndx])

        command="cut -d \",\" -f1,"+ listbestidkmer +" BEAUTY_RNASeq.csv" +" > matbestkm"+win
        print(command)
        os.system(command)
        '''

        fig = plt.figure(figsize=(13, 10))
        ax = fig.add_subplot(111)
        plt.plot(range(0, kmercount.shape[0]), kmercount["nview"])
        #plt.plot(range(0, kmercount.shape[0]), (kmercount["nview"]>np.mean(kmercount["nview"]))*500)
        if nbkmersel!=0:
            plt.plot(np.ones(2) * nbkmersel, [0 , np.max(kmercount["nview"])])
        if len(kmercount["nview"])>0:
            plt.scatter( np.argmin(np.array(kmercount["nview"])>np.mean(kmercount["nview"])), np.mean(kmercount["nview"]))
        if (nbkmersel!=0) & (np.max(kmercount["nview"])>nbkmersel*3):
            ax.set_xlim(0,nbkmersel*3)
        # y = mlab.normpdf(bins, mu, sigma)
        # l = plt.plot(bins)
        ax.set_ylabel("# of views")
        ax.set_xlabel("Kmers")
        if winneronly == False:
            ax.set_title(" Kmer exploration overall")
        else :
            ax.set_title(" Kmer exploration all Winner")
        if self.saveimg:
            win=""
            if winneronly == True:
                win="winners"
            else:
                win="all"

            plt.savefig(self.figdir + "/histokmerall"+win , dpi=300)
            plt.savefig(self.figdir + "/pdf/histokmerall"+win+".pdf", dpi=300)
        if self.showimg:
            plt.show()
        plt.close()

        jsond = ""
        jsond = win+str(hdinmin).replace(".","")+'explo={'
        jsond += 'title:"Kmers exploration",'
        if len(kmercount["nview"])>0:
            jsond += 'subtitle:"range:[%d-%d], mean %d, autotreshold: %d kmers select",'%(np.min(kmercount["nview"]),np.max(kmercount["nview"]),np.mean(kmercount["nview"]),nbkmersel)
        else:
            jsond += ''
        jsond += 'ylabel:"Number of selection",'
        jsond += "categories:[],\n"
        jsond += "series:["
        jsond += "{data:" + json.dumps(kmercount["nview"].tolist()) + ",\n"
        # jsond += "colorIndex:" + json.dumps(icat) + ",\n"
        jsond += 'name:"'+win+str(hdinmin)+' kmers "}\n'
        jsond += "]}"
        with open(self.figdir+"/explokmer"+win+str(hdinmin)+".json", 'w') as f:
            f.write(jsond)


        kmercount = kmercount.sort_values(by=['scoremean','nview'], ascending=False)
        kmercount = kmercount[kmercount["nview"] > (np.mean(kmercount["nview"])*2)]


        with open(self.figdir + "listkmer" + win + str(hdinmin) + "_SortByScore.fastq", "w") as file:
            with open(self.figdir + "countkmer" + win + str(hdinmin) + "_SortByScore.txt", "w") as filecount:

                for i in range(0, min(nbmaxkmer,kmercount.shape[0])):
                    #kmercount["nview"].values[i]

                   # print(i, " : ", kmercount["nview"].values[i], " = ",idseq2fasta(int(kmercount["kmer"].values[i]), 30))  # ,kmercount["kmer"].values[i]," : ",

                    print(">",win,i,"_ct", kmercount["nview"].values[i],"_scoremean", kmercount["scoremean"].values[i], "\n", self.idseq2fasta(int(kmercount["kmer"].values[i])),file=file)
                    print(kmercount["kmer"].values[i], ",", kmercount["scoremean"].values[i],file=filecount)

        fig = plt.figure(figsize=(13, 10))
        ax = fig.add_subplot(111)
        plt.plot(range(0, kmercount.shape[0]), kmercount["scoremean"])

        # ax = fig.add_subplot(212)
        # plt.plot(range(0, kmercount.shape[0]), kmercount["nview"])

        if self.saveimg:
            win = ""
            if winneronly == True:
                win = "winners"
            plt.savefig(self.figdir + "/histokmerall" + win+str(hdinmin)+"meanscoresort.png", dpi=300)
            plt.savefig(self.figdir + "/pdf/histokmerall" + win +str(hdinmin)+"meanscoresort"+ ".pdf", dpi=300)
        if self.showimg:
            plt.show()
        plt.close()

    ####
    # Create cluster of the best individuals looking at kmer intersection
    # @arg1 : nIndiv : number of individual to keep by folder
    # @arg2 : parameters for clustering function
    #       exemple: params = {'quantile': .3,
    #                           'eps': .3,
    #                           'damping': .9,
    #                           'preference': -200,
    #                           'n_neighbors': 10,
    #                           'n_clusters': 5}
    # @arg3 : nbmaxkmer : number of kmer max for fastq export
    # Ouput files :
    #   - _heatmapSortByAGandRuns__nIndiv(nIndiv).png  Intersection between solutions. Sort by folders and generation
    #   - full_heatmapCluster_(namealgo)(n_clusters_)_nIndiv(nIndiv).png Intersection between solutions. Sort byclustering algorithme "namealgo"
    #   - heatmapCluster_(namealgo)(n_clusters_)_nIndiv(nIndiv).png Same tha previous with snsheatmap
    #   - listkmerCluster_(namealgo)(n_clusters_)_nIndiv(nIndiv)_(nbmaxkmer).fastq    fastq list of kmer sort by number of view accros best individuals of the folders
    #   - _quantifcluster_(n_clusters_)_nIndiv(nIndiv).png Compare clustering methods
    #  Return: [nametab,countclusttab,meanfoldbyallclust]

    def clusteringBestIndiv(self,nIndiv,params,nbmaxkmer):
        # compute artificial mean to sort best individual 2/3test score, 1/3 outer score
        self.data["ArtiMeanscore"]=((self.data["scoremax"]*2)+self.data["scoreHDin"])/3
        listorgawin=self.data.sort_values(by="ArtiMeanscore",ascending=False)
        aaa = self.data.sort_values(by=["kmerlist", "ArtiMeanscore"], ascending=False)
        #listorgawin = self.data.sort_values(by="run",ascending=False)
        listorgawin = listorgawin.reset_index(drop=True)
        #listorgawin =listorgawin[:nIndiv]
        i=0
        for fold in np.unique(aaa.kmerlist):
            if i==0:
                listorgawin = (aaa[aaa["kmerlist"] == fold])[:nIndiv].sort_values(by=["run"], ascending=False)
            else:
                #listorgawin = listorgawin.append((aaa[aaa["kmerlist"] == fold])[:nIndiv])
                listorgawin = listorgawin.append(((aaa[aaa["kmerlist"] == fold])[:nIndiv]).sort_values(by=["run"], ascending=False),ignore_index=True)
            i=i+1
        listorgawin = listorgawin.reset_index(drop=True)
        commonmat=np.zeros( [len(listorgawin),len(listorgawin)],dtype="float")

        for i in np.arange(0,len(listorgawin)):
            commonmat[i,i]=1
            for j in np.arange(i+1, len(listorgawin)):
                commonmat[i,j]=self.comparekmlist(listorgawin.ix[i,"bestIndices"],listorgawin.ix[j,"bestIndices"])
                commonmat[j, i]=commonmat[i,j]
        fig = plt.figure(figsize=(13, 10))
        ax = fig.add_subplot(111)
        f = lambda x: x["kmerlist"].split("/")[-3]
        listorgawin["foldername"] = listorgawin.apply(f, axis=1)

        # display folders
        le = preprocessing.LabelEncoder()
        le.fit(np.unique(listorgawin["foldername"]))
        foldsint=le.transform(listorgawin["foldername"])
        displaymat = commonmat.copy()
        displaymat[displaymat == 0] = np.nan
        for instack in range(np.max([1, int(len(foldsint) / 100)])):
            displaymat = np.vstack([displaymat,foldsint/np.max(foldsint)])

        plt.imshow(displaymat, vmin=0, vmax=1, cmap='jet', aspect='auto', interpolation='nearest')  # cmap='hot'
        plt.colorbar()

        ax.set_ylabel("Individus of kmers")
        ax.set_xlabel("Individus of kmers")
        ax.set_title("Intersection between solutions\nSort by folders and generation")
        plt.savefig(
            self.figdir + "/_heatmapSortByAGandRuns_" +  "_nIndiv" + str(nIndiv)  + ".png",
            dpi=300)
        plt.close()
        #plt.show()

        #---------------------------------------------------------------------

        from sklearn.cluster import DBSCAN
        from sklearn import metrics
        # db = DBSCAN(eps=1, min_samples=8).fit(commonmat)
        # core_samples_mask = np.zeros_like(db.labels_, dtype=bool)
        # core_samples_mask[db.core_sample_indices_] = True
        # labels = db.labels_

        from sklearn import cluster, datasets, mixture
        from sklearn.neighbors import kneighbors_graph
        from sklearn.preprocessing import StandardScaler
        import warnings


        X = StandardScaler().fit_transform(commonmat)

        # estimate bandwidth for mean shift
        bandwidth = cluster.estimate_bandwidth(X, quantile=params['quantile'])

        # connectivity matrix for structured Ward
        connectivity = kneighbors_graph(
            X, n_neighbors=params['n_neighbors'], include_self=False)
        # make connectivity symmetric
        connectivity = 0.5 * (connectivity + connectivity.T)

        # ============
        # Create cluster objects
        # ============
        ms = cluster.MeanShift(bandwidth=bandwidth, bin_seeding=True)
        two_means = cluster.MiniBatchKMeans(n_clusters=params['n_clusters'])
        ward = cluster.AgglomerativeClustering(
            n_clusters=params['n_clusters'], linkage='ward',
            connectivity=connectivity)
        wardnoconnect = cluster.AgglomerativeClustering(
            n_clusters=params['n_clusters'], linkage='ward')
        spectral = cluster.SpectralClustering(
            n_clusters=params['n_clusters'], eigen_solver='arpack',
            affinity="nearest_neighbors")
        dbscan = cluster.DBSCAN(eps=params['eps'])
        affinity_propagation = cluster.AffinityPropagation(
            damping=params['damping'], preference=params['preference'])
        average_linkage = cluster.AgglomerativeClustering(
            linkage="average", affinity="cityblock",
            n_clusters=params['n_clusters'], connectivity=connectivity)
        birch = cluster.Birch(n_clusters=params['n_clusters'])
        gmm = mixture.GaussianMixture(
            n_components=params['n_clusters'], covariance_type='full')

        clustering_algorithms = (
            ('MiniBatchKMeans', two_means),
            ('AffinityPropagation', affinity_propagation),
            ('MeanShift', ms),
            ('SpectralClustering', spectral),
            ('Ward', ward),
            ('WardNoConnect', wardnoconnect),
            ('AgglomerativeClustering', average_linkage),
            ('DBSCAN', dbscan),
            ('Birch', birch),
            ('GaussianMixture', gmm)
        )
        meanfoldbyallclust=[]
        namemethod=[]
        nametab=[]
        countclusttab=[]
        for name, algorithm in clustering_algorithms:
            print("----"+name)
            t0 = time.time()

            # catch warnings related to kneighbors_graph
            with warnings.catch_warnings():
                warnings.filterwarnings(
                    "ignore",
                    message="the number of connected components of the " +
                            "connectivity matrix is [0-9]{1,2}" +
                            " > 1. Completing it to avoid stopping the tree early.",
                    category=UserWarning)
                warnings.filterwarnings(
                    "ignore",
                    message="Graph is not fully connected, spectral embedding" +
                            " may not work as expected.",
                    category=UserWarning)
                algorithm.fit(X)

            t1 = time.time()
            if hasattr(algorithm, 'labels_'):
                labels = algorithm.labels_.astype(np.int)
            else:
                labels = algorithm.predict(X)




            # Number of clusters in labels, ignoring noise if present.
            n_clusters_ = len(set(labels)) - (1 if -1 in labels else 0)

            print('Estimated number of clusters: %d' % n_clusters_)

            # print("Silhouette Coefficient: %0.3f"
            #       % metrics.silhouette_score(commonmat, labels))
            sortindex=np.argsort(labels)

            clustered = commonmat.copy()
            clustered = clustered[sortindex, :]
            clustered = clustered[:, sortindex]
            if np.max(labels)!=0:
                displaycluster=np.sort(labels)/np.max(labels)
            else:
                displaycluster = np.sort(labels)
            displaycluster[displaycluster<0]=-1



            #add cluster rows
            for instack in range(np.max([2,int(len(displaycluster)/70)])):
                clustered= np.vstack([clustered,displaycluster])

            #add folders rows
            if np.max(foldsint)!=0:
                foldshow = foldsint[sortindex] / np.max(foldsint)
            for instack in range(np.max([1, int(len(foldsint) / 100)])):
                clustered = np.vstack([clustered, foldshow])

            clustered[clustered==0]=np.nan
            fig = plt.figure(figsize=(13, 10))
            ax = fig.add_subplot(111)
            plt.imshow(clustered, vmin=0, vmax=1, cmap='jet', aspect='auto', interpolation='nearest')  # cmap='hot'
            plt.colorbar()
            ax.set_title("after clustering")
            plt.savefig(self.figdir +"/full_heatmapCluster_"+name + str(n_clusters_)+"_nIndiv" + str(nIndiv) +  ".png",
                        dpi=300)
            plt.close()
            # ax = fig.add_subplot(212)
            # plt.plot(np.sort(labels))
            sortlabel=np.sort(labels)
            # Create a categorical palette to identify the networks
            cluster_pal = sns.husl_palette(len(set(labels))-1)
            cluster_lut = dict(zip(map(int, np.unique(sortlabel[sortlabel>=0])), cluster_pal))

            #clustered=clustered[sortlabel>=0,:]
            #clustered = clustered[:, sortlabel >= 0]
            clustered = commonmat
            clustered = clustered[sortindex, :]
            clustered = clustered[:, sortindex]
            clusteredDF=pd.DataFrame(data=clustered, dtype=float)
            clusteredDF=clusteredDF.reset_index(drop=True)
            cluster_colors = pd.Series(sortlabel[sortlabel >= 0]).map(cluster_lut)
            cluster_colors=cluster_colors.rename("Clusters")

            #cluster_colors.rename_axis("Clusters")
            snsheatmap=sns.clustermap(clusteredDF, center=0.5, cmap="Blues", #jet,RdBu_r
                           row_colors = cluster_colors, col_colors = cluster_colors,
                           linewidths=0, figsize=(13, 13))
            snsheatmap.savefig(
                self.figdir + "/heatmapCluster_"+name + str(n_clusters_) + "_nIndiv" + str(nIndiv) + ".png")
            plt.close()
            #snsheatmap.close()
            #ax[0].imshow([colors], extent=[0, 10, 0, 1])
            #snsheatmap.dendrogram_row.reordered_ind
            listorgawin["cluster"]=labels
            if len(np.unique(sortlabel[sortlabel >= 0]))<50:
                for nclust in np.unique(sortlabel[sortlabel >= 0]):
                    orgaclust=listorgawin[listorgawin["cluster"]==nclust]['bestSequences']
                    aa = orgaclust.to_csv(sep=str(u','), index=False)
                    aa = aa.replace('\"\n\"', ',').replace('\"', '').replace('\n', '')
                    bestIndices = aa.split(',')
                    kmercount = pd.DataFrame.from_items([('nview', np.ones(len(bestIndices))),
                                                         ('bestIndices', bestIndices), ])
                    kmercount = kmercount.groupby('bestIndices', as_index=False).sum()
                    kmercount = kmercount.sort_values(by='nview', ascending=False)
                    with open(self.figdir + "/listkmerCluster_"+name + str(nclust)+"_nIndiv" + str(nIndiv) +"_"+ str(nbmaxkmer)+ ".fastq", "w") as file:
                        print("Write files"+self.figdir + "listkmerCluster" + str(nclust)+"_nIndiv" + str(nIndiv) + ".fastq")
                        for i in range(0, min(nbmaxkmer, kmercount.shape[0])):
                            print(">Cluster"+str(nclust)+"_count"+str(kmercount["nview"].values[i])+"\n"+
                                  self.idseq2fasta(int(kmercount["bestIndices"].values[i])), file=file)

            ############
            #list cluster config
            #listorgawin.sort_values(by=["cluster", "kmerlist"]).groupby("kmerlist")


            meanfoldbyclust=0
            ctclust = 0
            for clust in np.unique(listorgawin["cluster"]):
                print(clust)
                print((listorgawin[listorgawin["cluster"] == clust])["foldername"].unique())
                if clust>-1:
                    ctclust+=1
                    meanfoldbyclust+=len((listorgawin[listorgawin["cluster"] == clust])["foldername"].unique())
            namemethod.append(name+" "+str(ctclust))
            nametab.append(name )
            countclusttab.append(ctclust)
            if ctclust==0:
                meanfoldbyallclust.append(0)
            else:
                meanfoldbyallclust.append(meanfoldbyclust/ctclust)
        fig = plt.figure(figsize=(13, 10))
        ax = fig.add_subplot(111)
        plt.plot(range(len(meanfoldbyallclust)),meanfoldbyallclust)  # cmap='hot'
        ax.set_xticks(range(len(meanfoldbyallclust)+1))
        ax.set_xticklabels(namemethod,rotation=45)
        plt.grid()
        #plt.colorbar()
        ax.set_title("Clustering distribution")
        ax.set_xlabel("clustering method and number of cluster found")
        ax.set_ylabel("average number by cluster")
        plt.savefig(
            self.figdir + "/_quantifcluster_"  + str(n_clusters_) + "_nIndiv" + str(nIndiv) +".png",
            dpi=300)
        plt.close()


        return nametab,countclusttab,meanfoldbyallclust
            #---------------------------------------------------------------------
            # # sorting for second graph
            # mean0 = np.nanmean(commonmat, axis=0)
            # commonmat = np.vstack([commonmat, mean0])
            # mean1 = np.nanmean(commonmat, axis=1)
            # commonmat = np.c_[commonmat, mean1]
            #
            # sortindex = np.argsort(-mean0)
            # #listfold = np.asarray(listfold)[sortindex]
            # sortindex = np.append(sortindex, int(len(mean0)))
            #
            # commonmat = commonmat[sortindex, :]
            # commonmat = commonmat[:, sortindex]
            #
            # fig = plt.figure(figsize=(13, 10))
            # ax = fig.add_subplot(111)
            # plt.imshow(commonmat, vmin=0, vmax=1, cmap='jet', aspect='auto', interpolation='nearest')  # cmap='hot'
            # plt.colorbar()







#test=GeneticAlgPlotter()
#test.addDirectoryJson("resfake")
#test.savefigdir("fig_fake")
#test.matKmerFile="faketab.csv"
'''
test.retestres( 2, 5)
test.pca3D( 2, 5,True)
test.pca3D( 2, 5)
test.retestres( 2, 5)
test.retestres( 2, 5)
test.retestres( 2, 5)

'''
#test.addDirectoryJson("slurmdeepbigMoy")
#test.addDirectoryJson("resmokerADA")
#test.savefigdir("fig_resmokerADA")

#test.addDirectoryJson("resAltOutterPropfitnessrand3std",files="MPI_factnoise1g")
#test.savefigdir("fig_resAltOutterPropfitnessrand3std_noise1")
#test.plot2Dwinnersexploration(pKMstart = 4, pKMend = 12, pKMstep = 4, pNeurstart = 5, pNeurend = 6, pNeurstep = 5)
'''
#test.addDirectoryJson("restumor",files="*Pop100*")
test.addDirectoryJson("restumor",files="")
test.savefigdir("fig_restumor")
test.matKmerFile="test_tumor_ligth_noother.csv"


test.matKmerFile="InOUtMatrix_onlyInformativeKMers_ML.txt"
test.addDirectoryJson("resAltOutter4noexistJitter")
test.savefigdir("fig_resAltOutter4noexistJitter")

test.addDirectoryJson("KERAStest/cleanjson/")
test.savefigdir("KERAStest/fig_KERAStest")
'''
#test.addDirectoryJson("ADNcRes")
#test.savefigdir("ADNcRes/fig")



#old plottr moved to plot4folder

'''

DirectoryJson = "/home/sylvain/genologin/gecko/microrna/GAouttertest840_rankguy3longAlIAGADir/0_2/"
DirectoryJson = "/home/sylvain/musescratch/gecko/microrna/logmulticonv_ki0.2km5Dir/0_5/"
DirectoryJson = "/home/sylvain/genologin/gecko/tests/test1_Dir/0_4/"
DirectoryJson = "/home/sylvain/genologin/gecko/BeautyTNBC4/affinage1/run6geno_sh1no2Dir/0_5/"
DirectoryJson = "/home/sylvain/musescratch/gecko/BeautyResponseALL/affine15km/run1_km15Dir/0_5/"
DirectoryJson = "/home/sylvain/genologin/gecko/tests/test1_Dir/0_4/"
DirectoryJson = "/home/sylvain/genologin/gecko/ameva/big11_km20no1Dir/0_4/"


savefigdir= "/home/sylvain/musescratch/gecko/microrna/vfinal/V6_1200gel5ks50ki0.2mu0.5Tr0.8Dir/0_2/"
#
# nbBestIndiv = 10
# nbmaxkm=5000
# BestIndivfile= "bestindivBEAUTY.csv"
jsonfilesstart=""
nbBasePerKmer=30
'''
'''
if len(sys.argv) > 1:
    DirectoryJson = sys.argv[1]
if len(sys.argv) > 2:
    nbBestIndiv = int(sys.argv[2])
if len(sys.argv) > 3:
    nbmaxkm = int(sys.argv[3])
nbBasePerKmer=30
if len(sys.argv) > 4:
    nbBasePerKmer = int(sys.argv[4])
# if len(sys.argv) > 4:
#     jsonfilesstart = sys.argv[4]


# if len(sys.argv) > 5:
#     savefigdir= sys.argv[5]
# else:
savefigdir =DirectoryJson+"/fig"+jsonfilesstart+"/"

# if len(sys.argv) > 6:
#     BestIndivfile = sys.argv[6]
# else:
BestIndivfile = savefigdir+"BestIndiv"+str(nbBestIndiv)+".csv"
'''



'''
print("ClassGeneticAlgPlotter :")
print(DirectoryJson)
print(jsonfilesstart)
test.addDirectoryJson(DirectoryJson,jsonfilesstart)
test.savefigdir(savefigdir)
test.nbBasePerKmer=nbBasePerKmer


test.clusteringBestIndiv(1000)
'''
'''
print("test.generateTableOfBestIndiv(nbBestIndiv , BestIndivfile)")
test.generateTableOfBestIndiv(nbBestIndiv, BestIndivfile)
# test.scorehistory()
print("test.scorehistorymean()")
test.scorehistorymean()
print("test.overallHistokmerlist(winneronly=True,nbmaxkmer=nbmaxkm,hdinmin=0)")
test.overallHistokmerlist(winneronly=True,nbmaxkmer=nbmaxkm,hdinmin=0)
# print("test.overallHistokmerlist(winneronly=True,nbmaxkmer=nbmaxkm,hdinmin=0.8)")
# test.overallHistokmerlist(winneronly=True,nbmaxkmer=nbmaxkm,hdinmin=0.8)

print("test.overallHistokmerlist(winneronly=False,nbmaxkmer=nbmaxkm)")
test.overallHistokmerlist(winneronly=False,nbmaxkmer=nbmaxkm)
'''
'''






#test.overallHistokmerlist(winneronly=True,nbmaxkmer=nbmaxkm,hdinmin=0.8)





'''


