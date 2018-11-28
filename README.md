# GECKO

GECKO is composed of 2 main parts : 

* The import part where the raw sequencing data will be transformed into matrix
* The genetic algorithm to find most explicative kmers in the dataset


### Contact ###

* Written & Designed by : [Aubin Thomas](mailto:aubin.thomas@igh.cnrs.fr) & [Sylvain Barrière](mailto:sylvain.barriere@igh.cnrs.fr)
* Project instigator : [William Ritchie](mailto:william.ritchie@igh.cnrs.fr)




# Import - How to import data into matrix for GECKO

## Installation

### Requirements

Java 8/1.8, libboost, g++, Jellyfish 2, samtools and fastqc are needed to run the software. 


### Download the sources

you have to download by using git:

```bash
$ git clone https://gitlab.igh.cnrs.fr/aubin.thomas/importGECKO.git
```
A new directory called ImportGECKO will be created. You will find in it:

- the directory data contains example data in fastq format. 
- the directory bin that contains the code of the helpers.
- main.pl is the main script to be called
- .nf files are Nextflow scripts that are used to deal with workflows
- nextflow.config is the configuration file for Nextflow. It is the place to cinfigure the job queue manager and the resources to book for each process.


### Nextflow

You need to install [Nextflow](https://www.nextflow.io/). For this step you will need Java 8 or 1.8 installed. You have to go in the directory and install Nextflow by typing the commands:

```bash
$ cd importGECKO
$ curl -s https://get.nextflow.io | bash 
```

Very quickly, Nextflow enables scalable scientific worflows using software containers. The nextflow.config file is the main configuration file, and all the temporary files will be located inside the directory work/.


### Binaries

The last step is the compilation step of the helpers. At this step you need libboost and g++ installed on your machine:

```bash
$ cd bin/src
$ sh compile.sh
```
<br/>

## Configuration

### Path of binaries

Inside the file nextflow.config, first look at <b>params</b> section. Please change the value of softPath by the <i><u>absolute path</u></i> of the directory bin. The other parameters are default parameters that will be changed after the call of the main script.

### Job manager

Then have a look at the <b>process</b> section. By default the program is intended to run on a HPC managed by SLURM. If you have no job scheduling system installed on your computer you just have to delete the line where the executor parameter is mentioned.

Nextflow can handle a large variety of job scheduling systems by changing the 'slurm' value: 

* SGE : 'sge'
* LSF : 'lsf'
* PBS/Torque : 'pbs'
* NQSII : 'nqsii'
* HTCondor : 'condor'
* Ignite : 'ignite'

You can specify the name of the queue by adding the queue parameters inside the process section :

```
queue = 'long'
```

You can also add several queues at a time : 

```
queue = 'short,long,cn-el6'
```

### modules

Maybe you are working in an environement that needs the use of modules for Jellyfish, FastQC or samtools.
Inside the <b>process</b> section you will find call of modules that are comment with '//'. You just have to change the names of the modules and delete the double slash to uncomment.

<br/>

## Execution

### Decomposition into k-mers

The first step is the importation of sequencing data into a matrix of k-mers. The program takes as input files in fastq or bam format. The next command is to import fastq files :

```bash
perl main.pl decomposition --outdir <path of output folder> --reads 'path_to_fastq_files{1,2}.fastq'
```

To import bam files you need to add -B parameter :
```bash
perl main.pl decomposition --outdir <path of output folder> -B --reads 'path_to_fastq_files{1,2}.fastq'
```

the output folder contains several directories:

* fastqc : all the fastqc analyses of raw files.
* trim_galore : all the files are trimmed (to remove bad quality reads and adaptors). Then a new fastqc analysis is done.
* jellyfish/binary : all the k-mer counts in binary format (Jellyfish format)
* jellyfish/text : the respective text version of the files in jellyfish/binary


### Group file

From the jellyfish text files, a configuration file has to be created. It's to describe the groups. This file is a 3 columns tabulated file :

* <i><u>the absolute path </i></u> of the corresponding k-mer decomposition. 
* the name of the file (it will figure in the matrix as a reference name)
* the group (string or number without any space or special character)


Note : the configuration file can point to several decomposition directories.


### Importation

The second step is the importation of raw data. All the k-mer decompositions that figure in a single configuration file are merged inside a single matrix. The counts of the k-mers are normalized by the total number of millions of k-mer that are found in the corresponding sample. To launch the importation process:

```bash
./main.pl importation --groupconfig <path to group file> --outdir <path to output folder>

```

Notes : 

- The output folder can be same than in the decomposition step or be different.
- The process will consider k-mers that are above a threshold (currently fixed at 5, can be modified in 02_importation.nf at rawsubimport process)
- The program will split the configuration file into pieces of 5 files and hierarchically join the result sub-matrix. 

The program create a subdirectory called rawimport inside the output folder that is given in parameters. Inside this subdirectory the process creates 2 subdirectories :

- rawimport/conf that contains the splitted configuration files
- rawimport/matrix that contains at the end the result matrix : RAWmatrix.matrix



### Discretization

Then you need to discretize the matrix. An AMEVA discretization is performed on the matrix obtained previously. the output is located inside a subdirectory called discretization inside the output folder that is given in parameters.
The given command launches a discretization process :


```bash
./main.pl discretization --matrix <path to RAW matrix> --outdir <path to output folder>

```


### Filtering
Finally you need to filter and then transform the resulting matrix (that is discretized) into a matrix with real numbers. The filtering process will first remove k-mers that have less than 1/3*(number of samples inside the smallest group) of informations and then remove the k-mers that have the same information through the groups. This is the most time consuming process. That's why the first step is to divide the matrix into 1-millions-reads matrix and treat them separately.



```bash
./main.pl filter --matrix  <path to discrete matrix> --outdir <path to output folder>
./main.pl real --matrixDiscrete <path to discrete matrix> --matrixRaw <path to RAW matrix> --outdir <path to output folder>
```



<br/>
## Performances

### Performance overview

Each time a proces is terminated a directory called pipeline_info is created inside the output folder. It contains the complete trace of execution, and complete reports of the resources that have been consumed.


### Resource tuning

The pipeline is given with a file nextflow.config that configures a by default configuration for each process (number of CPU, memory). You can change all the values inside this file.


### Errors

Errors can occur on very long processes. If it happens you don't have to re-do the process from the beginning. You can add the option '-resume' at the call of main.pl and the program will try to satrt where it fell.


### Temporary files

Each step can create a lot of temporary files, above all the decomposition in k-mers. Each time these files are located on work directory. Feel free to delete this directory after each step. These temporary files are not automatically deleted in case you would need to resume a task.






# README Genetic Algorithms#

## Requirements ##
* python3 :
    * numpy
    * pandas
    * matplotlib
    * scikit-learn
    * mpi4py
    * seaborn
* gcc
* open MPI

 For the next steps, the repertory current is assume to be GECKO_source/algoGen/

## Configuration file ##
The genetic algorithm configuration it's done by a configuration file with .conf extension.
it's divided into 4 sections :

* GA parameters : allow you to set the kmer matrix file to use and all the parameters for the Genetic algorithm
* Machine learning : set the number of kmer wished and other Machine learning properties
* Output log properties
* Path to  start GA from pre-config state (optional)

This file can be easly generate with the form doc/form_config_tool.html
See guideline document for more details on the parameters and "Complementary ressources" section for example.


Multiple values of parameter are accepted if you use <b>multipleGeckostart.py</b> ( see below)  on a cluster the values need to be separated by ","


## Compilation ##
```bash
$ cd Producteurv2;make;cd ..

```
## Run ##
The best way to use our software is on cluster Slurm and Sun Grid Engine are supported. It also possible to launch the GA local architecture.
### Cluster ###

On Cluster, we recommend to use <b>multipleGeckostart.py</b> . It run GA from config file with multiple or single configuration inside and run automatically analyze after each GA are done. This script require a job scheduler, slurm by default also compatible with SGE.  

The multiple values of configuration must be separated by comma as:
```bash
kselection=20,50,100
```
The script will run a genetic algorithm for each parameters combinations
multipleGeckostart.py
Arguments:

* configuration file GA*
* number of base defined at the kmer extraction (default : 30)
* scheduler job submission command (default : sbatch)
* scheduler submission script for GA (default : GA_sbatch.sh)
* scheduler submission script for post analyze (default : nlyz_sbatch.sh)
 
 This two last scripts need to be manualy edit in view of make it work properly on your cluster. The values on CPU number,RAM memory, log path ,mail and eventually requiere module loading need to be adress. See "Complementary ressources" section for example.

Exemples :
```bash
python3 multipleGeckoStart.py configmulti.conf
python3 multipleGeckoStart.py configmulti.conf 30 qsub GA_qsub.sh nlyz_qsub.sh
```

You need first to make sure that the files scheduler submission scripts are properly filled:
The numbers of CPUs used for the GA ( we recommend between 5 to 10 CPUs for "LinSVC" machine learning with a population size arround 600, that can vary for optimal performance regarding your architecture, population size and machine learning technics).


Make sure the path directory to slurm log does exist.

### Local ###
#### prod_client_script_C++_V3.sh ####
That script is launching the genetic algorithm in c++ and the appropriate client for the machine learning.
By default all the CPU a use. you can specify a number of CPU can be set with bash variable $OMP_NUM_THREADS : export OMP_NUM_THREADS=1
It take the path to a config file for argument.
sh prod_client_script_C++_V3.sh MyConfig.conf

####  plotter_for_eachhistorylog.py ####
Generate full report result of GA accessible by the webpage evolution.html

Arguments:
* pathDirectory*
* file_csv* (input mattrice) could be remove in the future based on conffile...
* number of best individuals for details extraction* 
* Maximum of kmer use for mainfold (TSNE) visualization(default : 500)
* number of base defined at the kmer extraction (default : 30)
* testarg (default :''); others values [mainfoldonly,datafold,resumefold,plotfold,plot,heat,](default : empty mean all nlzy) 

Example: 
python3 plotter_for_eachhistorylog.py Beautypam50/maxkm/logmaxel10ks50ki.4mu.7tr.8Dir/ BEAUTY_RNASeq.csv 10 250 30
## Resultat visualisation ##
You can access to the GA result by the local HTML web page named <b>evolution.html</b> the configuration file, heatmap of the score across the population and the generation test score and outer. below there is a graph showing the discovering of kmer among the winner of each generation and when the most seen kmer become the most seen one. And finally, a graph display the manifold quality score at each checkpoint of you GA the manifold and score technics can be chosen with the select button.
On top of this page, the <b>result</b> button allow the user to visualize interactively the T-sne of the most viewed Kmer, the number of view of kmer sorted ( all and winner only). below the winner scores, history across the generation and the same graph sorted by test scores and outer scores. You will finally see below the T-sne of the 4 individuals who got the best scores

####  plotter_all_folders.py ####
Generate genetic algorithm report for all GA results inside the folder (via slurm script nlyz_sbatch.sh)

Same arguments than plotter_for_eachhistorylog.py but pathDirectory contain parent folder of the GA results



##  Multi runs analyzes ##
####  cmpfastqmine.py ####

Analyze cross folders will give you the intersection heatmap between each run for kmer winner most viewed, and the correlation for each modified GA parameters to scores, generation appearance of the most viewed kmer winner and the intersection for kmer winner most viewed.

Arguments :

* folder parent of resultAG or list of folder parent of resultAG

####  clustering_winner.py ####

Select the bests individuals of a run, and create a cluster of individual regarding the kmer composition proximity, then generate a list of Kmer for each cluster sort by the number of views.  
Arguments :

* folders GA path start (all folders with the same start will be included in the analyze "*" are accepted). Multiple start path accepted they must be separated by "," 
* number of base defined at the kmer extraction (default : 30)

####  occurencemerge.py ####

Analyze and combine most viewed winner of each generation to create a merged list of Kmer by averaging the number of views across the different runs, give you also the possibility to compare the results of each run compare to the merged result.
Arguments :

* folders GA path start (all folders with the same start will be included in the analyze "*" are accepted). Multiple start path accepted they must be separated by "," 
* list of number of k-mer selected by run, separated by ","(default : 5,10,20,30,40,50,60,100,200,500,1000,5000)
* list of outer score treshold, separated by "," (default : 0)
* file name of kmer occurence, if set empty (recomanded) the k-mer list will be recompute from orignal resultfile for the individuals winner (default : "")
 
# Complementary resources #
* Example of configuration file
* Example sbatch for genologin cluster on slurm
* GA result directory before the analyzes
* GA result directory after the analyzes 
 
## Example of configuration file ##
See guideline document for more details on the parameters.
This file can be easly generate with the form doc/form_config_tool.html
```bash
#==========GA Parameters==========

pathData=DatasetKmerCategories.csv

#--AlgorithmType  AGA  -> adaptatative genetic algorithm
#                 IAGA -> improve adaptatative GA              
#                 GA   -> fixed translocation and mutation rate TranslocationRate1 and mutationRate1
AlgorithmType=IAGA
#number of generation until GA stop
generation=30000
#number of individuals by generation 
individuals=600
elite=2

killRatio=0.4

#--mutation mode if mutationmode_1kmeronly = 1 mutation rate is chance to mutate 1 kmer of the individu 
#--in other cases mutation chance is apply to all kmer of the individu
mutationmode_1kmeronly=0 
mutationRate1=0.5
mutationRate2=0.1
mutationRate3=0.01

TranslocationRate1=0.8
TranslocationRate2=0.4
TranslocationRate3=0.7

kselection=50

outterRate=0.17
testRate=0.30
scoreMaxThreshold=1.1

#-------------------------------
#==========ML parameters========
kmer=20
hiddenNeuron=20
method=linSVC
#--test/train method : 1 -> to add noise
#                      2 -> triple the matrix samples and add noise
#                      3 -> no noise
shufTrainTestDataType=1
#Noise factor ratio of signal std
noisefactor=1
nRotativeTest=4

#-------------------------------
#==========Log config==========
pathLog=microrna/logmulticonv_

#--save history result
allResultHistoryBOOL=1
#--saving period
generationRefreshLog=6000
#--compute and save outter scores of each individus
computeAllOutters=1
#--Save detailed result of ML part for every individus (for debug only)
#detailResFile=

#-------------------------------
#==========Restart paths=========
#startFromPopulation=restest/logDir/0_1/_tableOfIndividuals.csv
#startFromOuter=BEAUTYtrinegRes/multiconv/_outerSelected.csv
#startfromKmerSubset=
```
## GA result directory before the analyzes ##
```bash
├── 0_1     #first breakpoint for the number of iterations set in the config by generationRefreshLog
│   ├── resultAG.json       #full log with winner composition and the score of the winner of each generation, and running time
│   ├── resultAG.json_OccKmers.csv      #sorted list of Kmer by number of views; columns: ID; number of view in individuals; cumulative score of each individual where the kmer appeared
│   └── _tableOfIndividuals.csv     #composition list of the population at the last generation (useful to restart GA from there)
├── 0_2     #second breakpoint for the number of iterations set in the config by generationRefreshLog*2
│   ├── resultAG.json
│   ├── resultAG.json_OccKmers.csv
│   └── _tableOfIndividuals.csv
├── AllOUTTERScoreByGeneration.csv      #outter scores of all individuals at each generation sort by score inner score
├── AllScoreByGeneration.csv        #inner scores of all individuals at each generation sort by score inner score
├── log     #log of GA run
├── multigenerated.conf     #configuration file
└── _outerSelected.csv      #samples select to compose the outer can be used to start another GA with the same outer

```
## Example sbatch for genologin cluster on slurm ##
#### GA_sbatch.sh ####
```bash
#!/bin/bash
#SBATCH -n 7                    # Number of cores
#SBATCH -N 1                    # Ensure that all cores are on one machine
#SBATCH -t 5-00:00              # Runtime in D-HH:MM
#SBATCH --mem=80000
#SBATCH -o /path/to/slurmlog/GAjob_%j.out      # File to which STDOUT will be written
#SBATCH -e /path/to/slurmlog/GAjob_%j.err      # File to which STDERR will be written
#SBATCH --mail-type=ALL         # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=username@lab.com  # Email to which notifications will be sent
module load cv-standard
module load python/3.5.2-bz2
module load openmpi/psm2/2.0.1

export OMP_NUM_THREADS=7  #the number of theads must fit the number of cores ask (line 2 :  -n 7)
echo $1 $2
sh prod_client_script_C++_V3.sh $1  > $2  2>&1
```
#### nlyz_sbatch.sh ####
```bash
#!/bin/bash
#SBATCH -n 1                    # Number of cores
#SBATCH -N 1                    # Ensure that all cores are on one machine
#SBATCH -p defq
#SBATCH -t 1-10:00              # Runtime in D-HH:MM
#SBATCH --mem=60000
#SBATCH -o /path/to/slurmlog/nlyz_%j.out      # File to which STDOUT will be written
#SBATCH -e /path/to/slurmlog/nlyz_%j.err      # File to which STDERR will be written
#SBATCH --mail-type=ALL         # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=username@lab.com  # Email to which notifications will be sent
module load cv-standard
module load python/3.5.2-bz2
echo $1 , $2, $3 ,$4 ,$5
python3 plotter_for_eachhistorylog.py $1 $2 $3 $4 $5 $6 >$1/nlyzlog  2>&1
```
## GA result directory after the analyzes ##
```bash
├── 0_1
│   ├── datahtml.json       #ressource for HTML interface 
│   ├── fig
│   │   ├── BestIndiv4.csv      #best individuals composition (Kmer id) 
│   │   ├── BestIndiv4.csv.fastq        #best individuals composition (fastq format) 
│   │   ├── BestIndiv4.csvforextractkm.count        #best individuals composition (format for count extraction) 
│   │   ├── BestIndiv4.csvforextractkm.count_SampleMat.csv      #count sub-matrix for best individuals kmers
│   │   ├── BestIndiv4.csvforextractkm.count_SampleMat.csv_orga_1.json      #ressource for HTML interface
│   │   ├── BestIndiv4.csvforextractkm.count_SampleMat.csv_orga_1MainfoldNeig30.png     #multiple technics manifold for the best individuals (image)
│   │   ├── BestIndiv4.csvforextractkm.count_SampleMat.csv_orga_1MainfoldNeig30.txt     #multiple technics manifold scores for the best individuals
│   │   ├── BestIndiv4.csvforextractkm.count_SampleMat.csv_orga_2.json      #as previous files but for the second best individuals
│   │   ├── BestIndiv4.csvforextractkm.count_SampleMat.csv_orga_2MainfoldNeig30.png
│   │   ├── BestIndiv4.csvforextractkm.count_SampleMat.csv_orga_2MainfoldNeig30.txt
│   │   ├── BestIndiv4.csvforextractkm.count_SampleMat.csv_orga_3.json
│   │   ├── BestIndiv4.csvforextractkm.count_SampleMat.csv_orga_3MainfoldNeig30.png
│   │   ├── BestIndiv4.csvforextractkm.count_SampleMat.csv_orga_3MainfoldNeig30.txt
│   │   ├── BestIndiv4.csvforextractkm.count_SampleMat.csv_orga_4.json
│   │   ├── BestIndiv4.csvforextractkm.count_SampleMat.csv_orga_4MainfoldNeig30.png
│   │   ├── BestIndiv4.csvforextractkm.count_SampleMat.csv_orga_4MainfoldNeig30.txt
│   │   ├── BestIndiv4.csv_score.txt        #score for the best individuals
│   │   ├── countkmerall0_SortByScore.txt       #score by kmer (id,score average) sort by score
│   │   ├── countkmerall0.txt       #Kmer count
│   │   ├── countkmerall0.txt_SampleMat.csv     #submatrix of count for the top Kmer list
│   │   ├── countkmerall0.txt_SampleMat.csv.json        #ressource for HTML interface
│   │   ├── countkmerall0.txt_SampleMat.csvMainfoldNeig30.png        #multiple technics manifold for the top Kmer list (image)
│   │   ├── countkmerall0.txt_SampleMat.csvMainfoldNeig30.txt       #multiple technics manifold scores for the top Kmer list
│   │   ├── countkmer_tresholdall0.txt      #same as the previous files for all top Kmer with above auto threshold base on the inflection point
│   │   ├── countkmer_tresholdall0.txt_SampleMat.csv
│   │   ├── countkmer_tresholdall0.txt_SampleMat.csv.json
│   │   ├── countkmer_tresholdall0.txt_SampleMat.csvMainfoldNeig30.png
│   │   ├── countkmer_tresholdall0.txt_SampleMat.csvMainfoldNeig30.txt
│   │   ├── countkmer_tresholdwin0.txt       #same as previous file for winners top Kmer with above auto threshold base on the inflection point
│   │   ├── countkmer_tresholdwin0.txt_SampleMat.csv
│   │   ├── countkmer_tresholdwin0.txt_SampleMat.csv.json
│   │   ├── countkmer_tresholdwin0.txt_SampleMat.csvMainfoldNeig30.png
│   │   ├── countkmer_tresholdwin0.txt_SampleMat.csvMainfoldNeig30.txt
│   │   ├── countkmerwin0.txt       #same as previous file for winners top Kmer 
│   │   ├── countkmerwin0.txt_SampleMat.csv
│   │   ├── countkmerwin0.txt_SampleMat.csv.json
│   │   ├── countkmerwin0.txt_SampleMat.csvMainfoldNeig30.png
│   │   ├── countkmerwin0.txt_SampleMat.csvMainfoldNeig30.txt
│   │   ├── countkmerwinners0_SortByScore.txt
│   │   ├── countkmerwin0.8.txt       #same as previous file for winners top Kmer with outer score>0.8 
│   │   ├── countkmerwin0.8.txt_SampleMat.csv
│   │   ├── countkmerwin0.8.txt_SampleMat.csv.json
│   │   ├── countkmerwin0.8.txt_SampleMat.csvMainfoldNeig30.png
│   │   ├── countkmerwin0.8.txt_SampleMat.csvMainfoldNeig30.txt
│   │   ├── countkmerwinners0.8_SortByScore.txt
│   │   ├── explokmerall0.json       #ressource for HTML interface
│   │   ├── explokmerwinners0.json       #ressource for HTML interface
│   │   ├── feature_distriball0.txt     #feature of count distribution of all Kmers explore
│   │   ├── feature_distribwin0.txt     #feature of count distribution of winner Kmers 
│   │   ├── histokmerall0meanscoresort.png      #sorted distribution of all Kmers inner scores average
│   │   ├── histokmerallall.png      #sorted distribution of all Kmers count ( number of views)
│   │   ├── histokmerallwinners0meanscoresort.png#sorted distribution of winners Kmers inner scores average
│   │   ├── histokmerallwinners.png      #sorted distribution of winners Kmers count ( number of views)
│   │   ├── historyscore.json       #ressource for HTML interface
│   │   ├── historyscore.png        #scores history across generation ( generation 0 as last generation)
│   │   ├── kmcut_selection_all.png
│   │   ├── kmcut_selection_win0.png
│   │   ├── listkmerall0.fastq       #most viewed sequences of all kmer (fastq format) 
│   │   ├── listkmerall0_SortByScore.fastq      #best average score of all kmer squences  (fastq format) 
│   │   ├── listkmerwin0.fastq       #most viewed sequences of winners kmer (fastq format) 
│   │   ├── listkmerwinners0_SortByScore.fastq      #best average score of winners kmer sequences  (fastq format) 
│   │   ├── listkmerwin0.8.fastq		#best average score of winners ( outterscore>0.8 ) kmer sequences  (fastq format)
│   │   ├── listkmerwinners0.8_SortByScore.fastq         #best average score of winners ( outterscore>0.8 ) kmer sequences  (fastq format) 
│   │   ├── pdf     #same figure than previous .png in PDF format
│   │   │   ├── histokmerall0meanscoresort.pdf
│   │   │   ├── histokmerallall.pdf
│   │   │   ├── histokmerallwinners0meanscoresort.pdf
│   │   │   ├── histokmerallwinners.pdf
│   │   │   ├── historyscore.pdf
│   │   │   └── winnerScores.pdf
│   │   ├── sortwinnerscore.json         #ressource for HTML interface
│   │   └── winnerScores.png        #winners scores history across generation ( generation 0 as the last generation)
│   ├── resultAG.json
│   ├── resultAG.json_OccKmers.csv
│   ├── scratch.html        #web page resume for this checkpoint (accessible via evolution.html)
│   └── _tableOfIndividuals.csv     #composition list of the population at the last generation (useful to restart GA from there)
├── 0_2
│   ├── datahtml.json
│   ├── fig
│   │   ├....
│   ├── resultAG.json
│   ├── resultAG.json_OccKmers.csv
│   ├── scratch.html
│   └── _tableOfIndividuals.csv
├── AllOUTTERScoreByGeneration.csv      #outter scores of all individuals at each generation sort by score inner score
├── AllOUTTERScoreByGeneration.csv_heatmap.png      #outter scores of all individuals at each generation sort by score inner score (heatmap)
├── AllScoreByGeneration.csv         #inner scores of all individuals at each generation sort by score inner score
├── AllScoreByGeneration.csv_heatmap.png     #inner scores of all individuals at each generation sort by score inner score (heatmap)
├── countkmerall0.txt_SampleMat.csvMainfoldNeig30.txt.json
├── countkmer_tresholdall0.txt_SampleMat.csvMainfoldNeig30.txt.json
├── countkmer_tresholdwin0.txt_SampleMat.csvMainfoldNeig30.txt.json
├── countkmerwin0.8.txt_SampleMat.csvMainfoldNeig30.txt.json
├── countkmerwin0.txt_SampleMat.csvMainfoldNeig30.txt.json
├── evolution.html      #main access HTML page for GA resume interface
├── Highcharts-6.0.4        #ressource for HTML interface
├── kmerOccurencesInWinners.pdf     #graph showing the discovering of kmer among the winner of each generation and when the most seen kmer become the most seen one
├── kmerOccurencesInWinners.png     #same than the previous file in png format
├── log      #log of GA run
├── multigenerated.conf     #configuration file
├── nlyzlog      #log of GA analyzes
├── OccurencesInWinnersstart.txt        #generation numbers when the most seen kmer become the most seen one (related to kmerOccurencesInWinners.png)
├── _outerSelected.csv      #samples select to compose the outer can be used to start another GA with same outer
├── scores_heatmapcompare.png       #heatmap inner, outer scores side to side of all individuals at each generation sort by score inner score 
├── scores_heatmapdiff.png       #heatmap inner - outer scores of all individuals at each generation sort by score inner score
├── scores_outter_corr.png      #correlation between outer and inner at each rank of individuals
├── scores_outter_corrlastquarter.png       #correlation between outer and inner at each rank of individuals for the last quarter of the generation
└── scores_outter_corrfirstquarter.png      #correlation between outer and inner at each rank of individuals for the first quarter of the generation
```


