from Bio import SeqIO
import numpy as np
import sys,glob,os,re,errno
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
import pandas as pd
import ClassGeneticAlgPlotter as gaplt
from sklearn.preprocessing import normalize

def cmpallfold(pathDirectory,maxkm=100000 ,figdir='',filesel ="listkmerwin0.fastq"):
    listdir = list()
    if type(pathDirectory)==str:
        listdirT = glob.glob(pathDirectory+"*/")
    else:
        listdirT = list()
        for path in pathDirectory:

            listdirT.extend(glob.glob(path + "*/"))
        # pathDirectory=pathDirectory[0]

    for fold1 in listdirT:
        sublistdir = glob.glob(fold1 + "[0-9]*_[0-9]*/")
        if len(sublistdir) > 0:
            listdir=np.append(listdir,fold1 )


        else:
            print("Exlude:" + fold1 + "\n no subdirectory")

    listdir.sort()
    resmat=np.zeros((len(listdir),len(listdir)))
    i1=0
    listfold=[]
    gene1 = np.array([])
    gene2 = np.array([])

    scoremax = np.array([])
    scoreHDin = np.array([])
    nbview = np.array([])
    #"fold":[""],
    emptyconfig=pd.DataFrame({"generation":[0],"shufTrainTestDataType":[0],"individuals":[0],"elite":[0],"kselection":[0],"killRatio":[0],"mutationRate":[0],"TranslocationRate":[0],"kmer":[0],"linsvn":[0],"mlp":[0],"nRotativeTest":[0],"noisefactor":[0],"hiddenNeuron":[0]})
    configs=pd.DataFrame()
    for fold1 in listdir:
        newconf=emptyconfig.copy()
        #newconf["fold"] =fold1
        conffile=glob.glob(fold1+"*.conf")
        if len(conffile)!=1:
            print( "Issue with conf for :"+fold1+"\n found:"+conffile)
        txtconf = np.loadtxt(conffile[0], delimiter='=', comments='#', dtype='S')
        for i in txtconf:
            key = i[0].decode("utf-8")
            if key in newconf.columns:
                newconf[key]=float(i[1])

            if key == "method":
                classfierType = i[1].decode("utf-8")

                if classfierType=="linsvn":
                    newconf["linsvn"] = 1
                if classfierType == "mlp":
                    newconf["mlp"] = 1

        configs=configs.append(newconf)




        listfold.append(fold1)
        print(fold1)
        sublistdir = glob.glob(fold1 + "/[0-9]*_[0-9]*/")
        numbersufold1 = np.array([], dtype="int")
        for fi in sublistdir:
            m = re.search(".*_([0-9]+)/*$", fi)
            numbersufold1 = np.append(numbersufold1, int(m.group(1)))

        # get best score for each folders
        filescore = glob.glob(fold1+"0_"+str(np.max(numbersufold1))+"/fig/BestIndiv*.csv_score.txt")
        if len(filescore)==0:
            2454
        print(fold1+"0_"+str(np.max(numbersufold1))+"/fig/BestIndiv*.csv_score.txt\n"+"filescore",filescore)
        res=pd.read_csv(filescore[0],sep=",")
        scoremax = np.append(scoremax,np.mean(res.loc[:]["scoremax"]))
        scoreHDin = np.append(scoreHDin,np.mean(res.loc[:]["scoreHDin"]))
        nbview = np.append(nbview, np.mean(res.loc[0:3]["nbview"]))
        # compute occurences kmer occurence

        if os.path.isfile(fold1+"OccurencesInWinnersstart.txt"):
            my_data = np.genfromtxt(fold1+"OccurencesInWinnersstart.txt", delimiter=',')
            print("OccurencesInWinnersstart.txt found")
            newbestgene=my_data[:,1]

        else:
            print("compute kmerOccurencesInWinners ")

            test = gaplt.GeneticAlgPlotter()
            DirectoryJson = fold1 + "0_" + str(np.max(numbersufold1))

            print(DirectoryJson)
            test.addDirectoryJson(DirectoryJson, "")
            test.savefigdir(fold1)
            newbestgene= test.kmerOccurencesInWinners()
        gene1 = np.append(gene1,newbestgene[0])
        gene2 = np.append(gene2, newbestgene[1])



        i2=0
        for fold2 in listdir:
            if fold1 == fold2:
                resmat[i1, i2] = float('NaN')
                ct=0
            if fold1<fold2:
                #get last folders result


                sublistdir = glob.glob(fold2 + "/[0-9]*_[0-9]*/")
                numbersufold2 = np.array([], dtype="int")
                for fi in sublistdir:
                    m = re.search(".*_([0-9]+)/*$", fi)
                    numbersufold2 = np.append(numbersufold2, int(m.group(1)))

                # (intersec, ct, ct2) = cmpfold(fold1 + "0_" + str(np.max(numbersufold1)) + "/fig/" + filesel,
                #                               fold2 + "0_" + str(np.max(numbersufold2)) + "/fig/" + filesel, maxkm)

                try:
                    # print("fold1 :", fold1)
                    # print("fold2 :", fold2)
                    (intersec,ct,ct2)=cmpfold(fold1+"0_"+str(np.max(numbersufold1))+"/fig/"+filesel, fold2+"0_"+str(np.max(numbersufold2))+"/fig/"+filesel,maxkm)
                except OSError as e:
                    print("error compare :")
                    print("fold1 :",fold1)
                    print("fold2 :", fold2)

                #resmat[i1, i2]=intersec/ct2
                resmat[i1, i2] = intersec / np.min([ct2,ct])
                resmat[i2, i1] = intersec / np.min([ct2, ct])
                #resmat[i1, i2] = intersec



                # if i2==0 :
                #     listfold[-1]+=" "+str(ct)
                #     if i1==1 :
                #         listfold[0] += " " + str(ct2)


            i2 += 1
        listfold[-1] += " " + str(ct)
        # if i1 == 1:
        #     listfold[0] += " " + str(ct2)






        i1+=1
    mean0=np.nanmean(resmat,axis=0)
    # mean1 = np.mean(resmat, axis=1)

    resmat=np.vstack([resmat, mean0])
    mean1 = np.nanmean(resmat, axis=1)
    resmat = np.c_[resmat, mean1]




    for sortm in ["","sort"]:

        fig = plt.figure(figsize=(13, 10))
        ax = fig.add_subplot(111)
        plt.subplots_adjust(bottom=0.25,left=0.25)
        plt.imshow(resmat,  aspect='auto', interpolation='none')  # cmap='hot'
        cb=plt.colorbar()
        cb.set_label("ratio of identicals kmers")
        ax.set_xticks(range(len(listfold)))
        ax.set_xticklabels(listfold, rotation=90 )
        ax.set_yticks(range(len(listfold)))
        ax.set_yticklabels(listfold, rotation=0)
        # listfold.insert(0, '')
        # ax.set_xticks(range(len(listfold)), listfold, rotation=45)
        #
        # ax.set_yticklabels(listfold, rotation=45 )
        # ax.set_ylabel("Generations")
        # ax.set_xlabel("Individus")
        ax.set_title("Intersection k-mer from results "+filesel)



        #plt.show()
        #figdir="/Users/sylvain.barriere/musescratch/gecko/microrna/"
        # figdir=pathDirectory
        plt.savefig(figdir + "cmpseqwinner"+str(maxkm)+sortm+".png", dpi=150)
        print("File result : ",figdir + "cmpseqwinner"+str(maxkm)+sortm+".png")
        plt.close()

        #sorting for second graph
        sortindex=np.argsort(-mean0)
        listfold = np.asarray(listfold)[sortindex]
        sortindex=np.append(sortindex,int(len(mean0)))

        resmat=resmat[sortindex,:]
        resmat = resmat[:, sortindex]
    #print(listfold)

    #print(resmat)

    configs.reindex()
    correlres=np.zeros(len(configs.columns))
    correlres_scmax =np.zeros(len(configs.columns))
    correlres_scoutter=np.zeros(len(configs.columns))
    correlres_nbview=np.zeros(len(configs.columns))
    correlres_gene1=np.zeros(len(configs.columns))
    correlres_gene2=np.zeros(len(configs.columns))
    i=0
    for key in configs.columns:
        #print(key)
       # tt = np.corrcoef(configs[key]/np.max(configs[key]),resmat[-1,:-1]/np.max(resmat[-1,:-1]))
        correlres[i]  = np.corrcoef(configs[key]/np.linalg.norm(configs[key]) ,
                                    resmat[-1, :-1]/np.linalg.norm(resmat[-1, :-1]) )[0,1]
        correlres_scmax[i] = np.corrcoef(configs[key] / np.linalg.norm(configs[key]),
                                         scoremax[:] / np.linalg.norm(scoremax[:]))[0, 1]
        correlres_scoutter[i] = np.corrcoef(configs[key] / np.linalg.norm(configs[key]),
                                            scoreHDin[:] / np.linalg.norm(scoreHDin[:]))[0, 1]
        correlres_nbview[i] = np.corrcoef(configs[key] / np.linalg.norm(configs[key]),
                                          nbview[:] / np.linalg.norm(nbview[:]))[0, 1]

        correlres_gene1[i] = np.corrcoef(configs[key] / np.linalg.norm(configs[key]),
                                         gene1[:] / np.linalg.norm(gene1[:]))[0, 1]
        correlres_gene2[i] = np.corrcoef(configs[key] / np.linalg.norm(configs[key]),
                                         gene2[:] / np.linalg.norm(gene2[:]))[0, 1]

        #= np.corrcoef(tt)[0, 1]
        #print(np.correlate(configs[key]/np.max(configs[key]),resmat[-1,:-2]/np.max(resmat[-1,:-2])))
        i+=1

    # fig = plt.figure(figsize=(13, 10))
    # ax = fig.add_subplot(111)
    # plt.subplots_adjust(bottom=0.25, left=0.25)
    correlres[np.isnan(correlres)]=0
    # plt.plot(correlres,configs.columns)
    # plt.grid(True)
    # plt.savefig(figdir + "cmpseqcorrel" + str(maxkm) + sortm + ".png", dpi=150)
    # #plt.show()
    # plt.close()
    # 111

    return correlres,correlres_scmax,correlres_scoutter,correlres_nbview,correlres_gene1,correlres_gene2,configs.columns


def cmpfold(file1,file2,maxkm):
    needed_reads = []
    reads_array = []
    chosen_array = []
    ct=0
    for x in SeqIO.parse(file1,"fasta"):
        reads_array.append(x.seq)
        ct+=1
        if ct==maxkm:
            break
    ct = 0
    for y in SeqIO.parse(file2,"fasta"):
        chosen_array.append(y.seq)
        ct+=1
        if ct==maxkm:
            break

    needed_reads = list(set(reads_array) & set(chosen_array))
    #print(needed_reads)
    # for x in SeqIO.parse(file1,"fasta"):
    #     if len(list(set([x.seq]) & set(needed_reads)))>0:
    #         print(x.description)
    # for x in SeqIO.parse(file2, "fasta"):
    #     if len(list(set([x.seq]) & set(needed_reads))) > 0:
    #         print(x.description)

    #print("nb seq commun",len(needed_reads),"/",min([len(reads_array),len(chosen_array)]))
    return (len(needed_reads),len(reads_array),len(chosen_array))
    #output_handle = open("Inetersec.fastq","w")

    #SeqIO.write(needed_reads,output_handle,"fasta")

    #output_handle.close()

#---------------------------------------------------------------------------
#---------------------------------------------------------------------------
pathDirectory="/Users/sylvain.barriere/musescratch/gecko/microrna/logmulticonv"
pathDirectory="/Users/sylvain.barriere/musescratch/gecko/microrna/maxkm/log"
#,"/Users/sylvain.barriere/musescratch/gecko/microrna/maxkm/log"
pathDirectory=["/Users/sylvain.barriere/musescratch/gecko/microrna/logmulticonv","/Users/sylvain.barriere/musescratch/gecko/microrna/vfinal/V"]
pathDirectory="/Users/sylvain.barriere/musescratch/gecko/microrna/vfinal/V"
pathDirectory="/home/sylvain/musescratch/gecko/BeautyTNBC4/affinage1/run5","/home/sylvain/musescratch/gecko/BeautyTNBC4/affinage1/run6"
pathDirectory=["/home/sylvain/musescratch/gecko/smockers/affinage_from5km/run4*sh1*no1","/home/sylvain/musescratch/gecko/smockers/affinage_from5km/run4*sh1*no0",
                "/home/sylvain/genologin/gecko/smockers/affinage_from5km/run4*sh1*no1","/home/sylvain/genologin/gecko/smockers/affinage_from5km/run4*sh1*no0"]
pathDirectory=["/home/sylvain/musescratch/gecko/BeautyTNBC4/affinage1/run5","/home/sylvain/musescratch/gecko/BeautyTNBC4/affinage1/run6","/home/sylvain/genologin/gecko/BeautyTNBC4/affinage1/run6","/home/sylvain/genologin/gecko/BeautyTNBC4/affinage1/run5"]
# pathDirectory=["/home/sylvain/musescratch/gecko/BeautyTNBC4/affinage1/run5*sh*no1","/home/sylvain/genologin/gecko/BeautyTNBC4/affinage1/run5*sh*no1",
#                "/home/sylvain/musescratch/gecko/BeautyTNBC4/affinage1/run5*sh*no0","/home/sylvain/genologin/gecko/BeautyTNBC4/affinage1/run5*sh*no0"]
# pathDirectory =               ["/home/sylvain/musescratch/gecko/BeautyTNBC4/affinage1/run6*sh2no0",
#                "/home/sylvain/genologin/gecko/BeautyTNBC4/affinage1/run6*sh2no0",
#                "/home/sylvain/musescratch/gecko/BeautyTNBC4/affinage1/run6*sh2no1",
#                "/home/sylvain/genologin/gecko/BeautyTNBC4/affinage1/run6*sh2no1"  ]
pathDirectory="/home/sylvain/musescratch/gecko/microrna/vfinal/V6_1200gel"
# pathDirectory= "/home/sylvain/musescratch/gecko/BeautyTNBC4/affinage1/run6_sh2no1Dir"
#pathDirectory =["/home/sylvain/musescratch/gecko/BeautyTNBC4/logN","/home/sylvain/musescratch/gecko/BeautyTNBC4/logtst31"]
#pathDirectory = "/home/sylvain/musescratch/gecko/BeautyTNBC4/logtst31el3ks50ki.3mu.8tr.8in840_rot5biskm5Dir"
#tabnbkm=[5,10,15]#,20,25,30,35,40,45,50,100,5000]
pathDirectory=["/home/sylvain/genologin/gecko/ameva/expo1_AlGAkm30","/home/sylvain/genologin/gecko/ameva/expo1_AlGAkm20"]
pathDirectory=["/home/sylvain/genologin/gecko/ameva/outt"]
pathDirectory = "/home/sylvain/genologin/gecko/microrna/replicatmasterhi"
#tabnbkm=np.arange(3,105,3)

# tabnbkm=np.arange(3,100,3)
tabnbkm=np.arange(3,100,20)

if len(sys.argv) > 1:
    pathDirectory = sys.argv[1]
    if type(pathDirectory) == str:
        pathDirectory = pathDirectory.split(',')


if type(pathDirectory)==str:
    # [figdir,tail]=os.path.split(pathDirectory)+"/"
    figdir=pathDirectory+"_cmpfastq/"
else:
    # [figdir, tail] =os.path.split(pathDirectory)
    # os.path.join(figdir,"multi")
    figdir = pathDirectory[0]+"_cmpfastqmulti/"
    
try:
    os.makedirs(figdir)
except OSError as e:
    if e.errno != errno.EEXIST:
        raise
    
correlres=np.array([])
correlres_scmax=np.array([])
correlres_scoutter=np.array([])
correlres_nbview=np.array([])

fig = plt.figure(figsize=(13, 10))

ax = fig.add_subplot(111)

plt.subplots_adjust(bottom=0.25, left=0.25)

corplot=[]
colors = plt.cm.hot(np.linspace(0,1,len(tabnbkm)+1))
plt.gca().set_color_cycle(colors)
filesel ="listkmerwin0.fastq"
for nbkm in tabnbkm:

    [cor,cor_scmax,cor_scoutter,cor_nbview,cor_gene1,cor_gene2,cols]=cmpallfold(pathDirectory,nbkm,figdir,filesel =filesel)

    if len(correlres)==0:
        correlres = cor
        correlres_scmax = cor_scmax
        correlres_scoutter = cor_scoutter
        correlres_nbview = cor_nbview
        correlres_gene1= cor_gene1
        correlres_gene2= cor_gene2

    else:
        correlres = np.vstack((correlres,cor))
        # correlres_scmax = np.vstack((correlres_scmax, cor_scmax))
        # correlres_scoutter = np.vstack((correlres_scoutter, cor_scoutter))

    corplot.extend(ax.plot(cor, cols, label=str(nbkm)))




#correlres[np.isnan(correlres)]=0
ax.legend(corplot)
#figdir="/Users/sylvain.barriere/musescratch/gecko/microrna/"
ax.grid(True)
plt.savefig(figdir + "cmpseqcorrelall.png", dpi=150)
plt.show()
plt.close()

fig = plt.figure(figsize=(13, 10))
ax = fig.add_subplot(211)
ax1 = fig.add_subplot(212)

#plt.subplots_adjust(bottom=0.25, left=0.25)
mask= np.invert(np.isnan(correlres_scmax))& (abs(correlres_scmax)>0.001)
if sum(mask)!=0 : 
    ax.plot(tabnbkm, correlres[:,mask])


    bigmat = np.vstack((correlres_scmax,correlres_scoutter,correlres_nbview,np.mean(correlres,axis=0),correlres_gene1,correlres_gene2))
    bigmat=bigmat[:,mask]
    cols=cols[mask]

    im=ax1.imshow( bigmat, aspect='auto', interpolation='none',cmap="seismic")
    clim=np.nanmax([np.nanmax(bigmat),-np.nanmin(bigmat)])
    #plt.clim(-clim,clim)
    cb = plt.colorbar(im)
    cb.set_clim(-clim, clim)
    cb.set_label("Coefficient of correlation")

    ax1.set_yticks(range(6))
    ax1.set_yticklabels(["score max","score outter","nb view 3best solutions","mean intersection","km most view1 appear","km most view2 appear"], rotation=0)
    ax1.set_xticks(range(len(cols)))
    ax1.set_xticklabels(cols, rotation=90)
    ax.set_title("Coorelation coef for intersection k-mer from results " + filesel)

    ax.set_title("Resume table of Coorelation coefficient")
    ax.set_xlabel("number of kmers selected")
    ax.set_ylabel("Coefficient of correlation")
    for (j,i),label in np.ndenumerate(bigmat):
        #ax1.text(i,j,'{:.3f}'.format(label),ha='center',va='center', fontsize= 18)
        ax1.text(i, j, '{:.3f}'.format(label), ha='center', va='center',color= 'w', fontsize= 20, bbox=dict(facecolor= 'black', alpha=0.8))

    ax.legend(cols)
    ax.grid(True)
    plt.savefig(figdir + "cmpseqcorrelallresume.png", dpi=150)
    plt.show()
    plt.close()

