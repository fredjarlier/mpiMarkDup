# mpiMarkDup
mpi and C based marking of duplicates for Next Generation Sequencing technology  

This is an open project initiated by Institut Curie HPC team and students from Paris Descartes University

The goal is to create a distributed and (near) real time version of the well known Picard MarkDuplicate from the Broad Institut. To do that we rely on MPI and C technology.   

The project is still under development and first results are encouraging. We are close to 100% reproducible.

The code is a fork from the mpiSORT. We add an extra process (see markduplicate.c) in the loop to manage the marking of duplicates read.

requirements
------------
gcc > 4.8 <br />
openmpi <br />
openssl <br />


How to compile:
--------------
aclocal <br />
autoconf <br />
./configure --prefix /usr/local CC=path_to_mpicc <br />
make <br />

How to test it
-------------

mpirun -n 8 mpiMD sam_file output_dir -q 0 -d 1000 -v 4 <br />

-q is for quality filtering <br />
-d is for optical distance duplicate <br />
-v is level of verbose <br />

How it works
------------

First the programm sort the reads by genome's coordinates and extract discordant and unmapped (end and mate) reads with the same technics describe in mpiSORT. <br />

Second the programm mark the duplicates for each chromosoms ans discordant reads according to Picard Markduplicate method. The unmapped and unmapped mate are not marked. <br />
To limit memory overhead we build a distributed perfect hash table to build fragment list and end list. So the memory usage stay  low.  <br />

Finally each chromosom is compressed with bgzf and writen down in the output folder.




