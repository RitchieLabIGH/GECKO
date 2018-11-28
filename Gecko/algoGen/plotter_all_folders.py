import os
import sys
import glob

file_csv = "faketabcom.txt"
nbIndiv = "30"
pathDirectory = sys.argv[1]
file_csv = sys.argv[2]
if len(sys.argv) > 3:
    nbIndiv = sys.argv[3]
maxkmmainfold = 500
if len(sys.argv) > 4:
    maxkmmainfold = int(sys.argv[4])
nbBasePerKmer='30'
if len(sys.argv) > 5:
    nbBasePerKmer = sys.argv[5]
testarg = ''
if len(sys.argv) > 6:
    testarg = sys.argv[6]


listdir = glob.glob(pathDirectory+"*/") #
for fi in listdir:

   #command="python3 plotter_for_eachhistorylog.py "+fi+" "+file_csv+" "+nbIndiv+" "+str(maxkmmainfold)+" "+testarg
   command = "sbatch nlyz_sbatch.sh  " + fi + " " + file_csv + " " + nbIndiv + " " + str(maxkmmainfold) + " " +nbBasePerKmer+" "+ testarg
   print(command)
   os.system(command)