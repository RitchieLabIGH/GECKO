import ClassGeneticAlgPlotter as gaplt
import sys,re,glob
import numpy as np
import pandas as pd
import matplotlib as mpl

import seaborn as sns
import os
if os.environ.get('DISPLAY','') == '':
    print('no display found. Using non-interactive Agg backend')
    mpl.use('Agg')
else:
    try:
        import tkinter
    except ImportError:
        mpl.use('Agg')

import matplotlib.pyplot as plt

###
# occurenMerge
# merge list of most viewed kmer along the winners of each generations
# @arg1 : pathDirectory : root path of directories to merge or array of root path of directories to merge
# @arg2 : nbKmByRun: number of km to keep by run ( 0 for all kmers)
# @arg3 : hdinmin : treshold for minimum outer score
# @arg4 : fileocckmers : file name of kmer occurence; force to load already generated file of kmer occurence if non empty. arg3 will have no effect in this case
# Output : all output are saved in first directory of the list of arg1
#   - mergedWinnerkmlis(hdinmin)_nbkm(nbKmByRun).fastq : merged list in fastq format
#   - heatmapMergedDistrib(hdinmin)_nbkm(nbKmByRun).png : heatmap of nunber of view of the merged list
#   - boxplotMergedDistrib(hdinmin)_nbkm(nbKmByRun).png : boxplot of nunber of view of the merged list
#   - redundancyPlotRankKmWinnerView(hdinmin)_nbkm(nbKmByRun).png : Kmer winner redundancy across AG runs
#   - redundancyKmWinnerView(hdinmin)_nbkm(nbKmByRun).png : Kmer winner redundancy across AG runs
#   - corrNviewViewinfold(hdinmin)_nbkm(nbKmByRun).png :correlation viewinfold Vs nview
#   - discoveredKmWinner(hdinmin)_nbkm(nbKmByRun).png : Proportion dicovered kmer by AG run
# return : synthesis array of the merging
#       ["hdinmin","nbKmByRun","nbfold","lenlist","meankmdisc","uniqkm","kmonce","corrrankVsViewinfold","pathDirectory","fileocckmers"]

def occurenMerge(pathDirectory,nbKmByRun=20,hdinmin=0.8,fileocckmers=""):
    if type(pathDirectory)==str:
        listdir = glob.glob(pathDirectory+"*/")
    else:
        listdir = list()
        for path in pathDirectory:

            listdir.extend(glob.glob(path + "*/"))
        pathDirectory=pathDirectory[0]


    listdir.sort()
    resmat=np.zeros((len(listdir),len(listdir)))
    i1=0
    listfold=[]

    resultglobal=pd.DataFrame([[hdinmin,nbKmByRun,0,0,0,0,0,0,pathDirectory,fileocckmers]],columns=["hdinmin","nbKmByRun","nbfold","lenlist","meankmdisc","uniqkm","kmonce","corrrankVsViewinfold","pathDirectory","fileocckmers"])
    for fold1 in listdir:
        if "*" in fold1 :
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


            else:
                print("Exlude:"+fold1+ "\n no subdirectory")


    listocc=pd.DataFrame()
    i=0
    nbkmerbyfold=np.array([])
    resultglobal["nbfold"]=len(listfold)
    for fold in listfold:
        if fileocckmers!="":
            listoccfile=pd.read_csv(fold+fileocckmers, sep=",", header=0,names=["bestIndices","nview"],index_col=False)
        else:
            i=i+1
            test = gaplt.GeneticAlgPlotter()
            test.addDirectoryJson(fold+"/")
            listoccfile=test.winnersexploration( hdinmin=hdinmin)

            #filter
            if    nbKmByRun!= 0:
                listoccfile = listoccfile.sort_values(by="nview", ascending = False)
                listoccfile = listoccfile.reset_index(drop=True)
                listoccfile = listoccfile.ix[0:nbKmByRun-1]


            listoccfile=listoccfile.drop('scoreCum',1)
            nbkmerbyfold=np.append(nbkmerbyfold,len(listoccfile))
        listoccfile["viewinfold"]=pd.Series(np.ones(len(listoccfile)))
        listoccfile["nview"] = (listoccfile["nview"] / listoccfile["nview"].sum())*100

        listoccfile[fold.split("/")[-2]]= listoccfile["nview"]
        #listoccfile["nview" + str(i)] = listoccfile["nview"]

        if len(listocc)==0:
            listocc=listoccfile

        else:
            #listoccfile=listoccfile.rename(index=str, columns={"nview": "nview"+str(i)})

            listocc=listocc.append(listoccfile, ignore_index=True)

    listoccsum=listocc.groupby('bestIndices', as_index=False).sum().sort_values(by="nview", ascending=False)
    listoccmean=listoccsum
    listoccmean["nview"]=(listoccsum["nview"]/len(listfold))
    listoccmean= listoccmean.sort_values(by="nview", ascending=False)
    listoccmean = listoccmean.reset_index(drop=True)
    if os.path.isdir(pathDirectory)==False:
        os.mkdir(pathDirectory)
    listoccmean[["bestIndices","nview"]].to_csv(pathDirectory+"/"+"mergedlistkmerwinner"+str(hdinmin)+"_nbkm"+str(nbKmByRun)+".csv",header=False,index=False)

    with open(pathDirectory +"/"+ "mergedWinnerkmlis" +str(hdinmin)+"_nbkm"+str(nbKmByRun) + ".fastq", "w") as file:
        print("Write files" +pathDirectory +"/"+str(hdinmin)+"_nbkm"+str(nbKmByRun)+ "mergedWinnerkmlis"   + ".fastq")
        for ii in range(0, listoccmean.shape[0]):
            if listoccmean["bestIndices"].values[ii].isdigit():
                print(">win" + str(ii) + "_count" + '{:.3}'.format(listoccmean["nview"].values[ii]) + "\n" +
                  test.idseq2fasta(int(listoccmean["bestIndices"].values[ii])), file=file)
            else:
                print(">win" + str(ii) + "_count" + '{:.3}'.format(listoccmean["nview"].values[ii]) + "\n" +
                      listoccmean["bestIndices"].values[ii], file=file)
    resultglobal["lenlist"] = listoccmean.shape[0]
    ###########################
    # fig, axes = plt.subplots(2, 2, figsize=(13, 10))

    fig = plt.figure(figsize=(13, 10))

    # pppp=listoccmean.drop(["bestIndices","viewinfold"], axis=1)
    # snsheatmap = sns.heatmap(pppp.ix[0:20])
    snsheatmap = sns.heatmap(listoccmean.drop(["bestIndices","viewinfold"], axis=1).ix[0:40].rename(index=str, columns={"nview": "mean nview"}),cmap="jet")
    plt.subplots_adjust(bottom=0.2)
    fig.savefig(pathDirectory +"/"+ "heatmapMergedDistrib"+str(hdinmin)+"_nbkm"+str(nbKmByRun)+ ".png")
    plt.close()
    ############################
    #box plot
    # missing zero value when not present
    fig = plt.figure(figsize=(13, 10))
    bplotmerged=listocc.merge(listoccmean.ix[0:39], left_on='bestIndices', right_on="bestIndices").sort_values(by="nview_y")
    # bplotmerged=listoccmean.ix[0:40].merge(listocc, left_on='bestIndices', right_on="bestIndices")
    bplotmerged.reset_index(drop=True)
    ax = sns.boxplot(x="bestIndices", y="nview_x", data=bplotmerged,order=listoccmean.ix[0:40].bestIndices)
    fig.savefig(pathDirectory +"/"+ "boxplotMergedDistrib"+str(hdinmin)+"_nbkm"+str(nbKmByRun)+ ".png")
    plt.close()



    ##################

    fig = plt.figure(figsize=(13, 10))
    ax = fig.add_subplot(111)
    plt.plot(listoccmean.viewinfold)
    ax.set_title("Kmers winners redundancy across AG runs\n"+"corr rank Vs viewinfold : "+str(np.corrcoef(listoccmean.index,listoccmean['viewinfold'])[0, 1]))
    ax.set_xlabel("Rank")
    ax.set_ylabel("Number AG where the kmer is a winner")
    plt.savefig(pathDirectory +"/"+"redundancyPlotRankKmWinnerView"+str(hdinmin)+"_nbkm"+str(nbKmByRun)+".png",dpi=300)
    plt.close()

    print("corr rank Vs viewinfold : ",np.corrcoef(listoccmean.index,listoccmean['viewinfold'])[0, 1])

    resultglobal["corrrankVsViewinfold"]= np.corrcoef(listoccmean.index,listoccmean['viewinfold'])[0, 1]
    #########################
    fig = plt.figure(figsize=(13, 10))
    ax = fig.add_subplot(111)
    corrcoeftab=np.array([])
    for j in range(0,i):

        plt.scatter(listoccmean.viewinfold, listoccmean.iloc[:, 1+j])
        c = np.corrcoef(listoccmean.viewinfold, listoccmean.iloc[:, 1+j])[0, 1]

        # plt.scatter(listoccmean.viewinfold,listoccmean["nview"+str(j)])
        # c=np.corrcoef(listoccmean.viewinfold, listoccmean["nview"+str(j)])[0, 1]
        corrcoeftab = np.append(corrcoeftab,c)
        print(j," corr viewinfold Vs nview = ",c)
    ax.set_title("Kmer winner redundancy across AG runs")
    ax.set_xlabel("Number AG where the kmer is a winner")
    ax.set_ylabel("Percent of kmer view along the winners")
    plt.savefig(pathDirectory +"/"+"redundancyKmWinnerView" +str(hdinmin)+"_nbkm"+str(nbKmByRun)+".png",dpi=300)
    plt.close()
    ##
    fig = plt.figure(figsize=(13, 10))
    ax = fig.add_subplot(111)

    plt.bar(range(0, len(corrcoeftab)), corrcoeftab)
    plt.xticks(range(0, len(corrcoeftab)), listoccmean.columns[1:-2], rotation=90)
    ax.set_title("correlation viewinfold Vs nview")
    ax.set_xlabel("AG runs")
    ax.set_ylabel("corr viewinfold Vs nview ")
    plt.subplots_adjust(bottom=0.2, top=0.97)

    plt.savefig(pathDirectory +"/"+ "corrNviewViewinfold"+str(hdinmin)+"_nbkm"+str(nbKmByRun) + ".png",dpi=300)

    plt.close()
    #########################

    print("mean :nbtotalkm / nbkmercumulé = ",(listoccmean["viewinfold"].sum()/len(nbkmerbyfold))/len(listoccmean))
    percentdiscovered=np.array([])
    for u in nbkmerbyfold:
        print("nbtotalkm / nbkmercumulé = ",(u/ len(listoccmean) ))
        if u!=0 :
            percentdiscovered = np.append(percentdiscovered,u/len(listoccmean) *100)
        else:
            percentdiscovered = np.append(percentdiscovered,0)
    fig = plt.figure(figsize=(13, 10))
    ax = fig.add_subplot(111)
    plt.bar(range(len(percentdiscovered)),percentdiscovered)
    #if nbKmByRun!=0:
    plt.plot([-0.5,len(percentdiscovered)-0.5],[np.mean(percentdiscovered),np.mean(percentdiscovered)],color='g')

    plt.plot([-0.5,len(percentdiscovered)-0.5],[(1/len(nbkmerbyfold))*100,(1/len(nbkmerbyfold))*100],color='b')
    plt.plot([-0.5,len(percentdiscovered)-0.5],[(listoccmean["viewinfold"].sum()/len(nbkmerbyfold))/len(listoccmean)*100,(listoccmean["viewinfold"].sum()/len(nbkmerbyfold))/len(listoccmean)*100],color='r')

    plt.xticks(range(0, len(corrcoeftab)), listoccmean.columns[1:-2], rotation=90)

    plt.subplots_adjust(bottom=0.2, top=0.92)
    ax.set_title("Proportion dicovered kmer by AG run ; mean kmmer discovered ="+str(np.mean(percentdiscovered))+"%\n"+str((len(listoccmean["viewinfold"])/listoccmean["viewinfold"].sum())*100)+
                 "% of unique kmer ; "+str((listoccmean["viewinfold"]>1).sum()/len(listoccmean)*100)+
                 "% of kmer view more than once")
    ax.set_xlabel("AG run")
    ax.set_ylabel("% of kmer discovered")
    plt.savefig(pathDirectory +"/"+ "discoveredKmWinner" +str(hdinmin)+"_nbkm"+str(nbKmByRun)+ ".png",dpi=300)
    plt.close()
    resultglobal["meankmdisc"]=np.mean(percentdiscovered)
    resultglobal["uniqkm"]=(len(listoccmean["viewinfold"])/listoccmean["viewinfold"].sum())*100
    resultglobal["kmonce"]=(listoccmean["viewinfold"]>1).sum()/len(listoccmean)*100

    return resultglobal

####
# the main function will generate multiple merging analyze to compare effect of the differents parameters
# @arg1 : pathDirectory : root path of directories to merge or array of root path of directories to merge
# @arg2 : list_nbKmByRun: list separate by comma : number of km to keep by run ( 0 for all kmers)
#           default : 5, 10, 20, 30, 40, 50, 60, 100, 200, 500, 1000, 5000
# @arg3 : list_hdinmin : list separate by comma : treshold for minimum outer score (default : 0)
# @arg4 : fileocckmers : file name of kmer occurence; force to load already generated file of kmer occurence if non empty. arg3 will have no effect in this case
# output : all files from occurenMerge() for each values of the list.
#           - (pathDirectory)_group/metaNlyzStore.csv data resume from all analyze
#           - (pathDirectory)_group/NlyzMerging.png graph representing similarities between GA for all differents parameters listed
#           - (pathDirectory)_group/NlyzMerging_zoomed.png same as previous , zoomed up to 200 nbKmByRun
#
if __name__ == '__main__':


    pathDirectory = ["/home/sylvain/musescratch/gecko/microrna/logmulticonv",
                     "/Users/sylvain.barriere/musescratch/gecko/microrna/maxkm/log"]
    pathDirectory = "/home/sylvain/genologin/gecko/ameva/big11_*no1"
    pathDirectory = "/home/sylvain/genologin/gecko/microrna/replicatmasterhi"
    # pathDirectory = "/home/sylvain/genologin/gecko/microrna/replicatmaster"
    pathDirectory = ["/home/sylvain/genologin/gecko/ameva/outtergood_hi","/home/sylvain/genologin/gecko/ameva/outterconverge_hi"]
    # pathDirectory = "/home/sylvain/genologin/gecko/BeautyTNBC4/logexplo_tst"
    # pathDirectory = "/home/sylvain/genologin/gecko/microrna/replicatconvergel3hi"#,"/home/sylvain/genologin/gecko/microrna/replicat*km5",replicatconvergMinushi]
    # # pathDirectory = "/home/sylvain/genologin/gecko/ameva/outterconverge_hi"
    # # pathDirectory = "/home/sylvain/genologin/gecko/BeautyTNBC4/logexplo_tst*no1"
    # pathDirectory = ["/home/sylvain/genologin/gecko/BeautyTNBC4/log_replicat_convergerkm5hi",
    #                  "/home/sylvain/genologin/gecko/BeautyTNBC4/log_replicat_convergerkm20hi",
    #    "/home/sylvain/genologin/gecko/BeautyTNBC4/log_replicat_convergerkm30hi",
    #                  "/home/sylvain/genologin/gecko/BeautyTNBC4/log_replicat_convergerkm50hi",
    #                  "/home/sylvain/genologin/gecko/BeautyTNBC4/log_replicat_convergerkm5sh",
    #                  "/home/sylvain/genologin/gecko/BeautyTNBC4/log_replicat_convergerkm20sh",
    #                  "/home/sylvain/genologin/gecko/BeautyTNBC4/log_replicat_convergerkm30sh",
    #                  "/home/sylvain/genologin/gecko/BeautyTNBC4/log_replicat_convergerkm50sh"]
                    # "/home/sylvain/genologin/gecko/BeautyTNBC4/log_replicat_convergerkm"]
#     pathDirectory = ["/home/sylvain/genologin/gecko/ameva/outterconverge_hi",
#                        "/home/sylvain/genologin/gecko/ameva/outterconverge_hi",
# "/home/sylvain/genologin/gecko/ameva/outterconverge_hi",
# "/home/sylvain/genologin/gecko/ameva/outterconverge_hi"]
#     pathDirectory = ["/home/sylvain/genologin/gecko/microrna/replicatconvergel10","/home/sylvain/genologin/gecko/microrna/replicatconvergel3hi",
#                     # "/home/sylvain/genologin/gecko/microrna/replicatmasterhi",
#     #"/home/sylvain/genologin/gecko/microrna/replicat*km5","/home/sylvain/genologin/gecko/microrna/replicat*km3",
#     "/home/sylvain/genologin/gecko/microrna/replicatconvergkm5","/home/sylvain/genologin/gecko/microrna/replicatconvergkm3",
#                      "/home/sylvain/genologin/gecko/microrna/replicatconvergMinushi"]
    pathDirectory="/home/sylvain/genologin/gecko/BeautyTNBC4/log_replicat_convergerhi"
    pathDirectory = "/home/sylvain/genologin/gecko/BeautyTNBC4/logSimpleGA_hi*no1_4paper_group" #,"/home/sylvain/genologin/gecko/BeautyTNBC4/logSimpleGA_hi*no0"]
    pathDirectory = ["/home/sylvain/genologin/CLL_mistery/testlength/test1_Re*Km7Dir","/home/sylvain/genologin/CLL_mistery/testlength/test1_Re*Km5Dir"]
    fileocckmers = ""
    list_nbKmByRun=[5, 10, 20, 30, 40, 50, 60, 100, 200, 500, 1000, 5000]
    #list_nbKmByRun=[ 20,  40,  60, 100, 200, 500, 1000, 5000]
    list_hdin = [0]
    force = False

    if len(sys.argv) > 1:
        pathDirectory = sys.argv[1]
    if type(pathDirectory) == str:
        pathDirectory = pathDirectory.split(',')



    if len(sys.argv) > 2:
        list_nbKmByRun = float(sys.argv[2].split(','))
    if len(sys.argv) > 3:
        list_hdin = float(sys.argv[3].split(','))
    if len(sys.argv) > 4:
        fileocckmers = sys.argv[4]

    oripathDirectory=pathDirectory
    resmat=pd.DataFrame()
    #print(occurenMerge(pathDirectory, nbKmByRun, hdinmin, fileocckmers))
    for foldgrp in oripathDirectory:
        if ((force==False) & (os.path.isfile(foldgrp+"/metaNlyzStore.csv"))):
            if resmat.empty:
                resmat = pd.read_csv(foldgrp + "/metaNlyzStore.csv")
            else:
                resmat =resmat.append(pd.read_csv(foldgrp + "/metaNlyzStore.csv"),ignore_index=True)
        else:
            i=0
            for hdinmin in list_hdin:
                for nbKmByRun in list_nbKmByRun:
                    #compute anlyyze occurence for the selected parameters
                    res=occurenMerge(foldgrp, nbKmByRun, hdinmin, fileocckmers)
                    if resmat.empty:
                        resmat=res
                    else:
                        resmat=resmat.append(res)
                    i=i+1

            resmat[resmat["pathDirectory"]==foldgrp].to_csv(foldgrp + "/metaNlyzStore.csv", index=False)
    if isinstance(pathDirectory, str)==False:
        if os.path.isdir(foldgrp+"group/") == False:
            os.mkdir(foldgrp+"group/")
        pathDirectory=foldgrp+"group/"
    resmat.to_csv(pathDirectory+"/metaNlyzStore.csv",index=False)
    fig, axes = plt.subplots(2, 2,figsize=(13, 10),sharex=True)
    alphag=0.8
    legend=[]
    for foldgrp in np.unique(resmat["pathDirectory"]):
        for hdinmin in list_hdin: #np.unique(resmat["hdinmin"]): #[0]:
            tmpmat=resmat[(resmat["hdinmin"]==hdinmin)&(resmat["pathDirectory"]==foldgrp)]
            axes[0,0].plot(tmpmat.nbKmByRun, tmpmat.meankmdisc,alpha=alphag,marker="x")
            axes[0,1 ].plot(tmpmat.nbKmByRun, tmpmat.uniqkm,alpha=alphag,marker="x")
            axes[1, 0].plot(tmpmat.nbKmByRun, tmpmat.kmonce,alpha=alphag,marker="x")
            axes[1, 1].plot(tmpmat.nbKmByRun, tmpmat.corrrankVsViewinfold,alpha=alphag,marker="x")
            if len(np.unique(resmat["pathDirectory"]))>1:
                legend.append(foldgrp.split("/")[-1]+" "+str(hdinmin))
            else:
                legend.append( "HDin > " + str(hdinmin))

    axes[0, 0].set_xlabel("Number of Kmers by run")
    axes[1, 0].set_xlabel("Number of Kmers by run")
    axes[0, 1].set_xlabel("Number of Kmers by run")
    axes[1, 1].set_xlabel("Number of Kmers by run")
    axes[0, 0].set_title("mean Kmers discovered")
    axes[0, 0].set_ylabel("% of kmer discovered")
    axes[0, 1].set_title("Unique kmers")
    axes[0, 1].set_ylabel("% of unique Kmers")

    axes[1, 0].set_title("Kmer view multiple time")
    axes[1, 0].set_ylabel("% of Kmers view more than once" )

    axes[1, 1].set_title("correlation rank Vs Viewinfold\n"+pathDirectory)
    axes[1, 1].set_ylabel("Correlation")
    #axes[1,1].set_xlim([0,200])

    fig.suptitle(pathDirectory)
    #f = lambda x: str(x["hdinmin"])+" "+x["pathDirectory"].split("/")[-1]
    #legendnames = resmat.apply(f, axis=1)
    #custom only
    #legendnames = resmat[resmat["hdinmin"] == 0].apply(f, axis=1)
    plt.legend(legend)
    plt.savefig(pathDirectory + "/" + "NlyzMerging" + ".png", dpi=300)
    axes[1, 1].set_xlim([0, 200])
    plt.savefig(pathDirectory + "/" + "NlyzMerging_zoomed" + ".png", dpi=300)
    #plt.show()
    plt.close()
    1
