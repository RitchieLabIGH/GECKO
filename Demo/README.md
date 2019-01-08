# Example : How to use GECKO

Here is an example of how to use GECKO, from raw FASTQ files to final analysis by running the genetic algorithm. 
This example is based on a subset of the microRNA data published by S. Juzenas & al [1].



## Data preparation

### Data download

Inside the demo folder are located several files. createdataset.sh is a script to download 90 of the fastq files of the project 
and run the entire process. For this you need to have sra toolkit to be installed. Just type :

```
$ sh createdataset.sh
```
You should have now 90 fastq files that represent 87Go of data. These data need to be cleaned.


### Gecko installation

During the data download you will have the time to install GECKO. Just clone the repository and go inside a matrix preparation 
part:

```
$ git clone https://github.com/RitchieLabIGH/GECKO.git
$ cd GECKO/ImportMatrix
```

You need to have Java8/1.8, libboost, g++ Jellyfish 2, samtools, cutadapt and fastqc installed in your environment. Then install Nextflow that will deal with the workflow:
```
$ curl -s https://get.nextflow.io | bash 
```

And compile the binaries:
```
$ cd bin/src
$ sh compile.sh
```

We assume that you are running the pipeline with a SLURM job manager. Otherwise please refer to Job manager section in the 
Readme.

If you have no job manager available, just remove the job manager option (executor option) in the nextflow.config file.

Inside the nexflow config file please modifiy the path parameter softPath to put the corresponding absolute path from your 
environment to point to the path of the binaries.


### From FASTQ to the matrix 
When the download has finished, you will have to clean the data and decompose the fastq files into kmers:
Go to the matrix preparation main directory:

```
$ cd ../..
```
And launch the decomposition process:
```
$ perl main.pl decomposition --singleEnd --outdir demo_import --reads '../../*fastq' --kmersize 20
$ rm work/* -rf
```

It will remove the Illumina adaptors and get the kmers of 20 nucleotides. The results will be copied into demo_import directory.

The rm command is to eliminate all the temporary files that have been created during the process.

Now edit the file microRNA.conf to modify the paths of the first column in function of your environment (change opt path). When finished, you can finish to import the data. First the importation: 
```
$ perl main.pl importation --groupconfig ../../microRNA_demo.conf --outdir demo_import
```

Then to filter you need to discretize the data:
```
$ perl main.pl discretization --matrix demo_import/rawimport/matrix/RAWmatrix.matrix –outdir demo_import
```
And filter the data:
```
$ ./main.pl filter --matrix demo_import/discretization/matrix/DISCRETmatrix.matrix  --outdir demo_import
```

And get the data to real counts:
```
$ ./main.pl real --matrixDiscrete demo_import/filtering/matrix/FILTEREDmatrix.matrix --matrixRaw demo_import/rawimport/matrix/RAWmatrix.matrix --outdir demo_import/
```

The final matrix is located in demo_import/filtering/final, and called FILTEREDmatrix_RealCounts.matrix. It is in ASCII format.

### To binary

The genetic algorithm takes as input the file that has been just generated. In this case the entire matrix will be loaded in memory that could be a bad idea for big matrix.

You can transform this matrix into binary format to overcome this issue. 

To do this, you have to go to the genetic algorithm part: 

```
$ cd ../Gecko
```

Then go to the utility section and compile the sources: 
```
$ cd algoGen/Producteurv2/utils
$ sh compile.sh
```

2 binaries are then created : transformIntoBinary to transform the text matrix into binary matrix, and indexBinary that takes 
as input a binary file to split it, that can increase the performance for very large matrix.

The demo matrix is so small that you can use directly the matrix as input.

To transform into binary, use transformIntoBinary : 
```
$ ./transformIntoBinary ../../../../ImportMatrix/demo_import/filtering/final/FILTEREDmatrix_RealCounts.matrix ../../../../ImportMatrix/demo_import/filtering/final/FILTEREDmatrix_RealCounts.bin
```

Now you have also the binary matrix  FILTEREDmatrix_RealCounts.bin that can be used as input for the genetic algorithm.

If you want to split this file to enroll multi-core access on the file, you can use the indexBinary executable : 
```
$ mkdir ../../../../ImportMatrix/demo_import/filtering/final/CutMatrix/
$ ./indexBinary ../../../../ImportMatrix/demo_import/filtering/final/FILTEREDmatrix_RealCounts.bin ../../../../ImportMatrix/demo_import/filtering/final/CutMatrix/example.bin 1000
```
The input matrix is divided here in 1000 kmers files, and the input file for the genetic algorithm would be /Absolutepath/CutMatrix/example.bin

We advise a number of 1 million of kmers per file.

You are now ready to execute the genetic algorithm.


## Genetic algorithm
### Install genetic algorithm
Before the first use or after update, you need to compile the Genetic algorithm :
```
$ cd .. ( you should be in Gecko/algoGen/Producteurv2 after this command)
$ make clean ( after update )
$ make
$ cd ..
```
### Execute genetic algorithm
Then you are back in the main directory of the genetic algorithm (Gecko/algoGen/), you can start it by using the example of configuration file in Demo/configGA_microRNA_demo.conf, all the parameters for the genetic algorithm are defined there, the path to the data, save path, and number of iterations for exemple  can be modified to fit your case :
```
$ sbatch multipleGeckoStart.py ../../Demo/configGA_microRNA_demo.conf 20
```
Then you can access the result report by default in ../../Demo/DemoGeneticAlgResultDir/evolution.html in your web browser and navigate to explore the behaviour of the GA and inspect the best classifiers.

