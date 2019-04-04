import time
import os,re
import sys
import glob
import mainfold
import numpy as np
import pandas as pd
import ClassGeneticAlgPlotter as gaplt

#test=False
test1_heatmap=False
test2_plotter=False
test3_sampledatamat=False
test3_smainfold = False
test4_smainfold = False

file_csv = "BEAUTY_RNASeq_NORM_discret_AMEVA_OnlyInformative_HammingMI_TripleNegative.csv"
file_csv = "faketabcom.txt"
nbIndiv="30"
nbmaxkmer="5000"


pathDirectory = sys.argv[1]
file_csv = sys.argv[2]
if len(sys.argv) > 3:
    nbIndiv = sys.argv[3]
maxkmmainfold=500
if len(sys.argv) > 4:
    maxkmmainfold= int(sys.argv[4])
nbBasePerKmer=30
if len(sys.argv) > 5:
    nbBasePerKmer = int(sys.argv[5])
testarg = ''
if len(sys.argv) > 6:
    testarg = sys.argv[6]

print("nbIndiv=" + nbIndiv)
print("maxkmmainfold=" + str(maxkmmainfold))
print("nbBasePerKmer=" + str(nbBasePerKmer))
print("testarg="+testarg)

if testarg=="mainfoldonly":
    test1_heatmap = True
    test2_plotter = True
    test3_sampledatamat = True
    test4_smainfold = False
if testarg == "datafold":
    test1_heatmap = True
    test2_plotter = True
    test3_sampledatamat = False
    test4_smainfold = False
if testarg == "resumefold":
    test1_heatmap = True
    test2_plotter = True
    test3_sampledatamat = True
    test4_smainfold = True
if testarg == "plotfold":
    test1_heatmap = True
    test2_plotter = False
    test3_sampledatamat = True
    test4_smainfold = False
if testarg == "plot":
    test1_heatmap = True
    test2_plotter = False
    test3_sampledatamat = True
    test4_smainfold = True
if testarg == "heat":
    test1_heatmap = False
    test2_plotter = True
    test3_sampledatamat = True
    test4_smainfold = True
if testarg == "last":
    1

listOfSampleList = ""

listconf = glob.glob(pathDirectory + "/*.conf")
txtconf=np.loadtxt(listconf[0],delimiter='=',comments='#',dtype='S')
for i in txtconf:
    key=i[0].decode("utf-8").lower()
    if key =="hiddenNeuron".lower():
        hideNeurons=int(i[1])
    if key =="shufTrainTestDataType".lower():
        shufTrainTestDataType=int(i[1])
    if key =="noisefactor":
        noisefactor=float(i[1])
    if key == "method":
        classfierType = i[1].decode("utf-8")
    if key == "detailResFile".lower():
        detailResFile = i[1].decode("utf-8")
    if key == "pathLog".lower():
        pathLog = i[1].decode("utf-8")+"Dir/"

    if key == "kmer".lower():
        nbkmer = int(i[1])


commandclient = "python3 heatmappopscore.py "+pathDirectory+"/AllScoreByGeneration.csv"
print(commandclient, flush=True)
if test1_heatmap == False:
    os.system(commandclient)
else:
    print("NOT start test mode activate\n", flush=True)
commandclient = "python3 heatmappopscore.py " + pathDirectory + "/AllOUTTERScoreByGeneration.csv"
print(commandclient, flush=True)
if test1_heatmap == False:
    os.system(commandclient)
else:
    print("NOT start test mode activate\n", flush=True)

commandclient = "python3 heatmappopoutter.py " + pathDirectory + "/"
print(commandclient, flush=True)
if test1_heatmap == False:
    os.system(commandclient)
else:
    print("NOT start test mode activate\n", flush=True)


listdir = glob.glob(pathDirectory + "/[0-9]*_[0-9]*/")


commandclient = "cp  evolution.html "+pathDirectory+"/evolution.html"

os.system(commandclient)
#commandclient = "cp  echarts.min.js "+pathDirectory+"/echarts.min.js"
commandclient = "cp  -R Highcharts-6.0.4 "+pathDirectory+"/Highcharts-6.0.4"
os.system(commandclient)

numbersufold=np.array([],dtype="int")

for fi in listdir:
    commandclient = "python3 plot4folder.py "+fi+" "+nbIndiv+" "+nbmaxkmer+" "+str(nbBasePerKmer)
    print(commandclient, flush=True)
    if test2_plotter == False:
        os.system(commandclient)
    else:
        print("NOT start test mode activate\n")
    m=re.search(".*_([0-9]+)/*$",fi)
    numbersufold=np.append(numbersufold, int(m.group(1)))

if testarg == "last":
    lastid=np.argmax(numbersufold)
    numbersufold=np.array([numbersufold[lastid]])


# number of view winner graph
test = gaplt.GeneticAlgPlotter()
DirectoryJson = pathDirectory + "/0_" + str(np.max(numbersufold))

print(DirectoryJson)
test.addDirectoryJson(DirectoryJson+"/", "")
test.savefigdir(pathDirectory+"/")
test.nbBasePerKmer=nbBasePerKmer
if test.nbBasePerKmer>0:
    test.kmerOccurencesInWinners()


#
# COPY HTML page
i=0
for fi in listdir:

    commandclient = "cp  scratch.html " + fi + "scratch.html"
    print(commandclient, flush=True)
    os.system(commandclient)
    commandclient = "cp  TSNEbest10.html " + fi + "TSNEbest10.html"
    print(commandclient, flush=True)
    os.system(commandclient)
    m = re.search(".*/([^/^_]*_)[0-9]+/*$", fi)
    txt = "datajson ={\n"
    if testarg != "last":
        txt +="\"prevlink\":\"../"+m.group(1)+str(numbersufold[i]-1)+"/scratch.html\",\n"
        txt += "\"nextlink\":\"../" + m.group(1) + str(numbersufold[i] +1) + "/scratch.html\",\n"
        txt += "\"maxfolder\":\""+str(np.max(numbersufold))+"\",\n"
        txt += "\"currentnb\":\"" + str(numbersufold[i]) + "\"\n"
    txt += "};\n"

    # for i in range(1,len(informations['run'])+1):
    #     informations['run'][i]['scoreHDin']=hdinscores[i-1]
    #     with open(fileLog, 'w') as f:
    #         f.write(json.dumps(informations, indent=4))
    with open(fi+"datahtml.json", 'w') as f:
        f.write(txt)
    i+=1



parmafile = ""
cpt=0
for fi in listdir:
    if cpt > 0:
        parmafile +=","
    parmafile+=fi+"fig/countkmerwin0.8.txt"
    cpt+=1
for fi in listdir:
    if cpt > 0:
        parmafile += ","
    parmafile += fi + "fig/countkmerwin0.txt"
    cpt += 1
for fi in listdir:
    if cpt > 0:
        parmafile += ","
    parmafile += fi + "fig/countkmerall0.txt"
    cpt += 1
for fi in listdir:
    if cpt > 0:
        parmafile += ","
    parmafile += fi + "fig/countkmer_tresholdall0.txt"
    cpt += 1





for fi in listdir:
    if cpt > 0:
        parmafile += ","
    parmafile += fi + "fig/countkmer0_SortByScore.txt"
    cpt += 1



for fi in listdir:
    if cpt > 0:
        parmafile += ","
    parmafile += fi + "fig/BestIndiv"+nbIndiv+".csvforextractkm.count"
    cpt += 1
for fi in listdir:
    if cpt > 0:
        parmafile += ","
    parmafile += fi + "fig/countkmerwin0_SortByScore.txt"
    cpt += 1
for fi in listdir:
    if cpt > 0:
        parmafile +=","
    parmafile+=fi+"fig/countkmerwin0.8_SortByScore.txt"
    cpt+=1



for fi in listdir:
    if cpt > 0:
        parmafile += ","
    parmafile += fi + "fig/countkmer_tresholdwin0.txt"
    cpt += 1



commandclient = "./Producteurv2/sampledatamat " +file_csv+" "+ parmafile + " 0"
print(commandclient, flush=True)
if test3_sampledatamat == False:
   os.system(commandclient)
else :
    print("NOT start test mode activate\n")



n_neighbors = 30
# for fi in listdir:
#     try:
#         mainfold.mainfoldtest(fi + "fig/countkmer_tresholdall0.txt_SampleMat.csv", n_neighbors, maxkmmainfold, pathDirectory)
#     except OSError as e:
#         print( "mainfold to "+fi + "fig/countkmer_tresholdall0.txt_SampleMat.csv failed")
if test4_smainfold == False :
    print("fig/BestIndiv"+nbIndiv+".csvforextractkm.count_SampleMat.csv", flush=True)
    for fi in listdir:
        try:
            mainfold.mainfold_organism(fi + "fig/BestIndiv"+nbIndiv+".csvforextractkm.count_SampleMat.csv",  n_neighbors, int(nbIndiv)* nbkmer ,pathDirectory,nbkmer)

        except OSError as e:
            print("mainfold to " + fi + "fig/BestIndiv"+nbIndiv+".csvforextractkm.count_SampleMat.csv failed", flush=True)

    print("fig/countkmerwin0.8.txt_SampleMat.csv", flush=True)
    for fi in listdir:
        try:
            mainfold.mainfoldtest(fi + "fig/countkmerwin0.8.txt_SampleMat.csv", n_neighbors, maxkmmainfold, pathDirectory)
        except OSError as e:
            print("mainfold to " + fi + "fig/countkmerwin0.8.txt_SampleMat.csv failed", flush=True)

    print("fig/countkmerwin0.txt_SampleMat.csv", flush=True)
    for fi in listdir:
        try:
            mainfold.mainfoldtest(fi + "fig/countkmerwin0.txt_SampleMat.csv",  n_neighbors,maxkmmainfold,pathDirectory)
        except OSError as e:
            print( "mainfold to "+fi + "fig/countkmerwin0.txt_SampleMat.csv failed", flush=True)

    print("fig/countkmerall0.txt_SampleMat.csv", flush=True)
    for fi in listdir:
        try:
            mainfold.mainfoldtest(fi + "fig/countkmerall0.txt_SampleMat.csv", n_neighbors, maxkmmainfold, pathDirectory)
        except OSError as e:
            print( "mainfold to "+fi + "fig/countkmerall0.txt_SampleMat.csv failed", flush=True)
    print("fig/countkmer0_SortByScore.txt_SampleMat.csv", flush=True)
    for fi in listdir:
        try:
            mainfold.mainfoldtest(fi + "fig/countkmer0_SortByScore.txt_SampleMat.csv", n_neighbors, maxkmmainfold,
                                  pathDirectory)
        except OSError as e:
            print( "mainfold to "+fi + "fig/countkmer0_SortByScore.txt_SampleMat.csv failed", flush=True)

    print("fig/countkmerwin0.8.txt_SampleMat.csv", flush=True)
    for fi in listdir:
        try:
            mainfold.mainfoldtest(fi + "fig/countkmerwin0.8_SortByScore.txt_SampleMat.csv", n_neighbors, maxkmmainfold,
                                  pathDirectory)
        except OSError as e:
            print("mainfold to " + fi + "fig/countkmerwin0.8_SortByScore.txt_SampleMat.csv failed", flush=True)

    print("fig/countkmerwin0.txt_SampleMat.csv", flush=True)
    for fi in listdir:
        try:
            mainfold.mainfoldtest(fi + "fig/countkmerwin0_SortByScore.txt_SampleMat.csv", n_neighbors, maxkmmainfold,
                                  pathDirectory)
        except OSError as e:
            print("mainfold to " + fi + "fig/countkmerwin0_SortByScore.txt_SampleMat.csv failed", flush=True)




    print("fig/countkmer_tresholdall0.txt_SampleMat.csv", flush=True)
    for fi in listdir:
        try:
            mainfold.mainfoldtest(fi + "fig/countkmer_tresholdall0.txt_SampleMat.csv", n_neighbors, maxkmmainfold, pathDirectory)
        except OSError as e:
            print( "mainfold to "+fi + "fig/countkmer_tresholdall0.txt_SampleMat.csv failed", flush=True)

    print("fig/countkmer_tresholdwin0.txt_SampleMat.csv", flush=True)
    for fi in listdir:
        try:
            mainfold.mainfoldtest(fi + "fig/countkmer_tresholdwin0.txt_SampleMat.csv", n_neighbors, maxkmmainfold, pathDirectory)
        except OSError as e:
            print( "mainfold to "+fi + "fig/countkmer_tresholdwin0.txt_SampleMat.csv failed", flush=True)

import MLevaluation as MLev
for fi in listdir:
    for idindiv in range(int(nbIndiv)):
        MLev.MLevaluation.generateNlzReportClassifierQuality( os.path.join(fi, "fig/BestIndiv10.csvforextractkm.count_SampleMat.csv"), extendsavepath="/indiv"+str(idindiv+1)+"/", nrot=0, hideNeurons=int(nbkmer/2), kmByIndiv=nbkmer, idindiv=idindiv, eliminate=False,classfierType=classfierType)

filestruct=pd.DataFrame()
#historymain fold result compil
i=0
for fi in listdir:
    listfiles = glob.glob(fi + "fig/countkmer*.txt_SampleMat.csvMainfoldNeig*.txt")
    for fitxt in listfiles:
        i+=1
        m = re.search(".*_([0-9]+)/*$", fi)
        path, file = os.path.split(fitxt)
        #label[m.group(1)] = np.genfromtxt(fitxt, delimiter=",")[:,1]
        #rawtxt = np.genfromtxt(fitxt, delimiter=",", dtype='unicode')
        if i==1 :
            filestruct=pd.read_csv(fitxt, delimiter=",", index_col=0, header=None).T
            filestruct["id"]=int(m.group(1))
            filestruct["file"] = file
        else:
            tmpfilestruct=pd.DataFrame()
            tmpfilestruct= pd.read_csv(fitxt, delimiter=",",index_col=0, header=None).T
            tmpfilestruct["id"] = int(m.group(1))
            tmpfilestruct["file"] = file
            filestruct=filestruct.append(tmpfilestruct)

filestruct=filestruct.sort_values(by="id",ascending=True)
for nfile in np.unique(filestruct["file"]):
    #filestruct[filestruct["file"] == nfile].to_csv(os.path.join(pathDirectory, nfile), index=False)
    a=filestruct[filestruct["file"]==nfile].set_index("id")
    # =a.T.to_json(os.path.join(pathDirectory, nfile+".json"), orient="index")
    json=a.T.to_json( orient="split")#index
    with open(os.path.join(pathDirectory, nfile+".json"), 'w') as f:
        f.write("data="+json)


