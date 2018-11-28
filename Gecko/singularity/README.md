# Singularity recipe for GECKO_source
There are two recipes: 
- singularity_base.recipe : create the image with ubuntu and all the dependencies installed. To add a new dependency, append it to the apt-get install command.
- singularity.recipe : use the ubuntu image created with singularity_base.recipe to build the latest version of GECKO_source.

The division in two is aimed to reduce the building time: you won't have to install all the dependency all the time that you make some change in the program. 

The only dependency is singularity. Super user permission are required just to build the images, not to use them.

## Build Commands
First build the base image:
 
```bash
rm -f ./GECKO_base_img && sudo singularity build GECKO_base_img ./singularity_base.recipe
```
Then you can use it to build the operational image:

```bash
rm -f ./GECKO && sudo singularity build GECKO ./singularity.recipe
```

## Use the image

Singularity documentation can be found @ https://www.sylabs.io/docs/
You can open a shell in the image:

```bash
singularity shell GECKO
```
or execute a command. All the content of algoGen is in /GECKO/ and the program are already built for the use.

```bash
singularity exec GECKO mpirun -np 3 python3 /GECKO/clientNNc++.py uid gecko.conf  
```

By default, the folder where you launch the singularity image is visible by the image and also the classical path ( /home/ /tmp/ ... ).
To bind other path, you can use the  environmental variable SINGULARITY_BINDPATH or SINGULARITY_BIND like this: 

```bash
SINGULARITY_BINDPATH="/poolzfs/"
```

The variable accept comma separated locations *path[:des[:mode]]*, where path is the path to the directory that you want to bind, des is the destination folder in the image (by default is the same as path), and mode can be are the permission (by default *wr*)




