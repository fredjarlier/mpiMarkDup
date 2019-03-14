# mpiMarkDup
MPI and C based programm of marking duplicates for Next Generation Sequencing technology.  

This is an open project initiated by Institut Curie HPC team and students from Paris Descartes University.

The goal is to create a distributed and (near) real time version of the well known Picard MarkDuplicate from the Broad Institute. To do that we rely on MPI and C technology and parallel algorithms.   

The project is still under development and first results are encouraging.

The code is a fork from the mpiSORT. We add an extra process (see markduplicate.c) in the loop to manage the marking of duplicates read.

requirements
------------
gcc > 4.8 <br />
openmpi <br />
openssl <br />
cmocka (for unit testing)

How to compile:
--------------
aclocal <br />
autoconf <br />
autoreconf <br />
automake --add-missing
./configure --prefix /usr/local CC=path_to_mpicc <br />
make <br />
make install <br />
make clean <br />

How to test it
-------------

mpirun -n 8 mpiMD sam_file output_dir -q 0 -d 1000 -v 4 <br />

-q is for quality filtering <br />
-d is for optical distance duplicate <br />
-v is level of verbose <br />
    0 is LOG_OFF  <br />
    1 is LOG_ERROR  <br />
    2 is LOG_WARNING  <br />
    3 is LOG_INFO (default) <br />
    4 is LOG_DEBUG  <br />
    5 is LOG_TRACE  <br />

How it works
------------

First the programm sort the reads by genome's coordinates and extract discordant and unmapped (end and mate) reads with the same technics describe in mpiSORT. <br />

Second the programm mark the duplicates for each chromosom ans discordant reads according to Picard Markduplicate method. The unmapped and unmapped mate are not marked. To limit memory overhead we build a distributed perfect hash table for fragment list and end list this way the memory usage stay low.  <br />

Finally each chromosom is compressed with bgzf and writen down in the output folder.




