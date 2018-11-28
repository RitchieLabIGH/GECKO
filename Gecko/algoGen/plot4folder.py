import ClassGeneticAlgPlotter as gaplt
import sys


DirectoryJson = "BEAUTYres"
savefigdir= "BEAUTYres/fig"

nbBestIndiv = 10
nbmaxkm=5000
BestIndivfile= "bestindivBEAUTY.csv"
jsonfilesstart=""

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

print("ClassGeneticAlgPlotter :")
print(DirectoryJson)
print(jsonfilesstart)
test = gaplt.GeneticAlgPlotter()
test.addDirectoryJson(DirectoryJson,jsonfilesstart)
test.savefigdir(savefigdir)
test.nbBasePerKmer=nbBasePerKmer


print("test.generateTableOfBestIndiv(nbBestIndiv , BestIndivfile)")
test.generateTableOfBestIndiv(nbBestIndiv, BestIndivfile)
# test.scorehistory()
print("test.scorehistorymean()")
test.scorehistorymean()
print("test.overallHistokmerlist(winneronly=False,nbmaxkmer=nbmaxkm)")
test.overallHistokmerlist(winneronly=False,nbmaxkmer=nbmaxkm)
print("test.overallHistokmerlist(winneronly=True,nbmaxkmer=nbmaxkm,hdinmin=0)")
test.overallHistokmerlist(winneronly=True,nbmaxkmer=nbmaxkm,hdinmin=0)
print("test.overallHistokmerlist(winneronly=True,nbmaxkmer=nbmaxkm,hdinmin=0.8)")
test.overallHistokmerlist(winneronly=True,nbmaxkmer=nbmaxkm,hdinmin=0.8)
# print("test.kmerOccurencesInWinners()")
#test.kmerOccurencesInWinners()