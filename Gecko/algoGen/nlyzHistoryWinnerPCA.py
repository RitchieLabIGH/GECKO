'''
Generate mainfold of best individuals at periodical generations in view of display historical evolution of the solution
accross the GA. Start at the first generation every step ask to the maximum generation ask.

@arg1 : path to the result folder
@arg2 : path to count matrix file
@arg3 : number of kmer by individual
@arg4 : number of generations between each mainfold analyze
@arg5 : limit of generation to stop the historic
@arg6 : base length of the K-mer (default : 30)



Exemple to generate mainfold every 5 generations until generation 50 for inviduals of 10 kmers with 30 bases by kmer
$ python3 nlyzHistoryWinnerPCA.py resultfolder_dir/0_1/ countmatrix.csv 10 5 50 30

Result can be consult in the new folder resultfolder_dir/0_1/bestindivhist_step5_max50/ by tsnerender_history.html ( display 10 mainfold ) or extractwinnerhisto_forextractkm.count_SampleMat.csv_orga_{organumber}MainfoldNeig30.png (for each individuals)


'''

import ClassGeneticAlgPlotter as gaplt
import mainfold
import sys,os


DirectoryJson = "BEAUTYres"

if len(sys.argv) > 1:
    DirectoryJson = sys.argv[1]
pathData =""
if len(sys.argv) > 2:
    pathData = sys.argv[2]
nbkmer=5
if len(sys.argv) > 3:
    nbkmer = int(sys.argv[3])

step = 1
if len(sys.argv) > 4:
    step = int(sys.argv[4])
maxgene=20
if len(sys.argv) > 5:
    maxgene = int(sys.argv[5])
nbBasePerKmer=30
if len(sys.argv) > 6:
    nbBasePerKmer = int(sys.argv[6])
# if len(sys.argv) > 4:
#     jsonfilesstart = sys.argv[4]


# if len(sys.argv) > 5:
#     savefigdir= sys.argv[5]
# else:
savefigdir =DirectoryJson+"bestindivhist_step{0}_max{1}/".format(step,maxgene)




# if len(sys.argv) > 6:
#     BestIndivfile = sys.argv[6]
# else:


print("ClassGeneticAlgPlotter :")
print(DirectoryJson)

test = gaplt.GeneticAlgPlotter()
test.addDirectoryJson(DirectoryJson,"")
test.savefigdir(savefigdir)
test.nbBasePerKmer = nbBasePerKmer
test.extracthistoryBestIndiv(step,maxgene,"extractwinnerhisto_")





parmafile=savefigdir+"extractwinnerhisto_forextractkm.count"


commandclient = "./Producteurv2/sampledatamat " + pathData + " " + parmafile + " 0"
print(commandclient, flush=True)
os.system(commandclient)

n_neighbors = 30
mainfold.mainfold_organism( parmafile + "_SampleMat.csv", n_neighbors,
                           1000, savefigdir+"historyPca", nbkmer)
commandclient = "cp  tsnerender_history.html "+savefigdir+"tsnerender_history.html"

os.system(commandclient)
# tsnerender.html