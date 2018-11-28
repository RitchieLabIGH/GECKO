'''
Clustering for AG  run , get the best individu winners (number is folowing a preset grid of values)
try to create cluster of solution by looking at the instersection matrix of the individu
@agr1 : common start path for AG folders, or list of them separte by a coma
@arg2 : base length of the K-mer (default : 30)


'''


import ClassGeneticAlgPlotter as gaplt
import sys,re,glob,gc,os
import numpy as np
import pandas as pd
import seaborn as sns
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


pathDirectory=["/home/sylvain/musescratch/gecko/microrna/logmulticonv","/Users/sylvain.barriere/musescratch/gecko/microrna/maxkm/log"]
#pathDirectory= "/home/sylvain/genologin/gecko/ameva/big11_"
pathDirectory = "/home/sylvain/genologin/gecko/microrna/replicatmasterhi"
pathDirectory = "/home/sylvain/genologin/gecko/ameva/outtergood_hi"
#pathDirectory = "/home/sylvain/genologin/gecko/BeautyTNBC4/log_replicat_convergerkm30h"
pathDirectory = "/home/sylvain/genologin/gecko/microrna/replicatconvergel10"
pathDirectory = "/home/sylvain/genologin/gecko/ameva/outterconverge_hi"
pathDirectory = "/home/sylvain/genologin/gecko/BeautyTNBC4/log_replicat_convergerkm30h"
pathDirectory = "/home/sylvain/genologin/gecko/BeautyTNBC4/logSimpleGA_hi*o1"
if len(sys.argv) > 1:
    pathDirectory = sys.argv[1]
    if type(pathDirectory) == str:
        pathDirectory = pathDirectory.split(',')

if type(pathDirectory)==str:
    listdir = glob.glob(pathDirectory+"*/")
else:
    listdir = list()
    for path in pathDirectory:

        listdir.extend(glob.glob(path + "*/"))
    pathDirectory=pathDirectory[0]
nbBasePerKmer = 30
if len(sys.argv) > 2:
    nbBasePerKmer = int(sys.argv[2])

pathDirectory=pathDirectory+"_cluster_win/"

#grid of number of best indiv select by GA(folder)
tabindivnumber=[50,100,200] #[10,20]:#,30,50,100,500,1000,2000,5000]:



listdir.sort()
resmat=np.zeros((len(listdir),len(listdir)))
i1=0
listfold=[]
test = gaplt.GeneticAlgPlotter()
test.savefigdir( pathDirectory)

test.nbBasePerKmer=nbBasePerKmer

for fold1 in listdir:

    if "*" in fold1:
        print("Exlude:" + fold1 + "\n star in the name")
    else:

        #get last folders result
        sublistdir = glob.glob(fold1 + "/[0-9]*_[0-9]*/")
        numbersufold1 = np.array([], dtype="int")
        for fi in sublistdir:
            m = re.search(".*_([0-9]+)/*$", fi)
            numbersufold1 = np.append(numbersufold1, int(m.group(1)))
        if len(numbersufold1)>0:
            listfold.append(fold1+"0_"+str(np.max(numbersufold1)))

            test.addDirectoryJson(fold1+"0_"+str(np.max(numbersufold1))+"/","")
        else:
            print("Exlude:"+fold1+ "\n no subdirectory")




params = {'quantile': .3,
          'eps': .3,
          'damping': .9,
          'preference': -200,
          'n_neighbors': 10,
          'n_clusters': 5}
foldbyclust=pd.DataFrame()
for iindiv in tabindivnumber:
    [namemethod,countclusttab,meanfoldbyallclust]=test.clusteringBestIndiv(iindiv,params,100)
    foldbyclust=foldbyclust.append(pd.DataFrame(data=[meanfoldbyallclust],columns=namemethod,index=[iindiv]))
fig = plt.figure(figsize=(13, 10))
sns.heatmap(foldbyclust,cmap="jet")
plt.subplots_adjust(bottom=0.2)
fig.savefig(test.figdir + "/" + "__clustquantifheatmap100" + ".png")
plt.close()
gc.collect()

foldbyclust=pd.DataFrame()
params = {'quantile': 1,
          'eps': 1,
          'damping': .9,
          'preference': -300,
          'n_neighbors': 10,
          'n_clusters': 10}
for iindiv in tabindivnumber:#,30,50,100,500,1000,2000,5000]:
    [namemethod, countclusttab, meanfoldbyallclust]=test.clusteringBestIndiv(iindiv,params,101)
    foldbyclust=foldbyclust.append(pd.DataFrame(data=[meanfoldbyallclust],columns=namemethod,index=[iindiv]))
fig = plt.figure(figsize=(13, 10))
sns.heatmap(foldbyclust,cmap="jet")
plt.subplots_adjust(bottom=0.2)
fig.savefig(test.figdir + "/" + "__clustquantifheatmap101" + ".png")
plt.close()



#MiniBatchKMeans
#ward
