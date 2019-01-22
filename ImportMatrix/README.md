# How to import data into matrix for GECKO

## Installation

### Requirements

Java 8/1.8, libboost, g++, Jellyfish 2, samtools and fastqc are needed to run the software. 


### Download the sources

you have to download by using git:

```bash
$ git clone https://github.com/RitchieLabIGH/GECKO.git
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
$ cd GECKO/ImportMatrix
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

## Anova Filter
A quick method to select k-mers is through an ANOVA F-test, removing the kmers with a p-value higher than a given threshold. This method may be used alone or before a Mutual Information filter. 
A drawback of this method is that you may loose some informative combination of kmers, expecially if you use a stringent threshold (Skip to "Mutual Information Filter" below below to avoid those loose)
To launch this process:

```bash
./main.pl anova --matrix <path to the raw matrix> --threshold [float] --outdir <path to output folder>

```

We suggest a threshold of 0.05.
The output matrix will be in /path/to/outdir/filtered.matrix


## Mutual Information Filter
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
