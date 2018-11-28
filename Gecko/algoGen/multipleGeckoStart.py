import sys,os,time
import subprocess
import numpy as np
import re


def recursiveprintingconf(txtconf,maxlvl,nbBasePerKmer='30',confFilename='',confContent='',lvl=0,pathData='',scriptAG='',scriptNlyz='',schedsubmit='',nlyzjobid=""):
    idkey=0
    i = txtconf[0,:]
    idkey+=1
    key = i[0].decode("utf-8")
    values = i[1].decode("utf-8")
    values = values.split(',')
    if len(values)>1:
        values
    for val in values:

        newconfFilename = confFilename
        if key == "pathData":
            pathData=val
        if key == "pathLog":
            newconfFilename = val+confFilename
            newconfContent=confContent
        else:
            if len(values)>1:
                newconfFilename+=key[0:2]+str(val)

            newconfContent=confContent+key+"="+str(val)+"\n"

        #recursiveprintingconf(txtconf[idkey:],maxlvl,newconfFilename,newconfContent,lvl+1,command)
        if lvl==maxlvl:

            newconfContent += "pathLog=" + newconfFilename + "\n"
            if not os.path.exists(newconfFilename+"Dir"):
                os.makedirs(newconfFilename+"Dir")
            f = open(newconfFilename+"Dir/multigenerated.conf", 'w')
            f.write(newconfContent)
            f.close()
            commandclient = schedsubmit+" "+scriptAG+" "+newconfFilename+"Dir/multigenerated.conf "+newconfFilename+"Dir/log "
            print(commandclient)
            time.sleep(1)
            # os.system(commandclient)
            output = subprocess.check_output(commandclient, shell=True)
            print (output)
            g = re.search('[0-9]+', output.decode("utf-8"))  # capture the inner number only
            jobid=g.group(0)
            if "sbatch" in schedsubmit:
                dep="--dependency=afterok:"
            else:
                dep="-hold_jid "
            commanNlyz=schedsubmit+" "+dep+jobid+" "+scriptNlyz+" "+newconfFilename+"Dir/ "+pathData+" 10 250 "+nbBasePerKmer
            print(commanNlyz)
            output = subprocess.check_output(commanNlyz, shell=True)
            print(output)
            g = re.search('[0-9]+', output.decode("utf-8"))  # capture the inner number only

            if len(nlyzjobid)==0:
                nlyzjobid = g.group(0)
            else:
                if "sbatch" in schedsubmit:
                    nlyzjobid += ":"+ g.group(0)
                else:
                    nlyzjobid += "," + g.group(0)




        else:
            # print("newconfFilename  else : lvl "+str(lvl)+" lvlmax:"+str(maxlvl))
            tmpnlyzjobid= recursiveprintingconf(txtconf[idkey:], maxlvl, nbBasePerKmer,newconfFilename, newconfContent, lvl + 1, scriptAG=scriptAG,scriptNlyz=scriptNlyz,schedsubmit=schedsubmit,pathData=pathData,nlyzjobid=nlyzjobid)
            nlyzjobid = tmpnlyzjobid

    return nlyzjobid

#=======================================================================================================================
if __name__ == '__main__':
    if len(sys.argv) > 1:
        conffile= sys.argv[1]
    nbBasePerKmer='30'
    if len(sys.argv) > 2:
        nbBasePerKmer= sys.argv[2]
    schedsubmit="sbatch"
    if len(sys.argv) > 3:
        scriptNlyz= sys.argv[3]
    if "sbatch"in schedsubmit:
        scriptAG = "GA_sbatch.sh"
        scriptNlyz = "nlyz_sbatch.sh"
    else:
        scriptAG = "GA_qsub.sh"
        scriptNlyz = "nlyz_qsub.sh"
    if len(sys.argv) > 4:
        scriptAG= sys.argv[4]
    if len(sys.argv) > 5:
        scriptNlyz= sys.argv[5]



    #subprocess.check_output("module load python/3.5.2-bz2", shell=True)
    txtconf=np.loadtxt(conffile,delimiter='=',comments='#',dtype='S')
    print(txtconf)
    print(recursiveprintingconf(txtconf,txtconf.shape[0]-1,nbBasePerKmer,scriptAG=scriptAG,scriptNlyz=scriptNlyz,schedsubmit=schedsubmit))#,confFilename=os.path.dirname(conffile)+"/"
