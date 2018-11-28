'''
This script will run all GA step for voting mode on a cluster from a configuration file with multiple configuration
- multiple genetic algorithm + analyze
- analyze difference between genetic algorithms
- clustering winner analyze
- merging list of all the genetic algorithm
@arg1 : configuration file
@arg2 : base length of the K-mer (default : 30)
@arg3 : scheduler job submission command (default : sbatch)
@arg4 : cluster GA script (default : GA_sbatch.sh or GA_qsub.sh if scheduler different than sbatch)
@arg5 : cluster analyze script (default : nlyz_sbatch.sh or nlyz_qsub.sh if scheduler different than sbatch)
@arg6 : voting anayle script (default : voting_nlyz_sbatch.sh)



'''
import multipleGeckoStart as Mgs
import subprocess
import sys
import numpy as np
if __name__ == '__main__':
    if len(sys.argv) > 1:
        conffile = sys.argv[1]
#    if len(sys.argv) > 2:
#        numberofreplicat = sys.argv[2]
    nbBasePerKmer = '30'
    if len(sys.argv) > 2:
        nbBasePerKmer = sys.argv[2]
    schedsubmit = "sbatch"
    if len(sys.argv) > 3:
        scriptNlyz = sys.argv[3]
    if "sbatch" in schedsubmit:
        scriptAG = "GA_sbatch.sh"
        scriptNlyz = "nlyz_sbatch.sh"
    else:
        scriptAG = "GA_qsub.sh"
        scriptNlyz = "nlyz_qsub.sh"
    if len(sys.argv) > 4:
        scriptAG = sys.argv[4]
    votingscriptNlyz="voting_nlyz_sbatch.sh"
    if len(sys.argv) > 5:
        votingscriptNlyz = sys.argv[5]

txtconf=np.loadtxt(conffile,delimiter='=',comments='#',dtype='S')
joblist=Mgs.recursiveprintingconf(txtconf,txtconf.shape[0]-1,nbBasePerKmer,scriptAG=scriptAG,scriptNlyz=scriptNlyz,schedsubmit=schedsubmit)
if "sbatch" in schedsubmit:
    dep = "--dependency=afterok:"
else:
    dep = "-hold_jid "
#get pathLog value
for i in txtconf:
    key = i[0].decode("utf-8")
    if key =="pathLog" :
        pathLog = i[1].decode("utf-8")



commanNlyz = schedsubmit + " " + dep + joblist + " " + votingscriptNlyz + " " + pathLog +" "+ nbBasePerKmer
print(commanNlyz)
output = subprocess.check_output(commanNlyz, shell=True)
print(output)