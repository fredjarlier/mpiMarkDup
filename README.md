# mpiMarkDup

MPI and C based program of marking duplicates for Next Generation Sequencing technology.<br /> 

The code is a fork from the mpiSORT project. We add a process to manage the marking of duplicate reads.

The goal is to create a distributed and (near) real time version of the well known Picard MarkDuplicate from the Broad Institute. To do that we rely on low level  technologies and parallel algorithms.<br />

The project is still under development so feel free to test it and report.

This is an open project initiated by Institut Curie HPC team and students from Paris Descartes University.

Release notes 
-------------

Release 1.0 from the 06/05/2019 <br />

1) Fix tiebreak (tested only with power of 2 cpu). <br />

Release 1.0 from the 08/04/2019 <br />

1) fix reproducibility issue but still have a corner case.  <br />


Release 1.0 from the 31/03/2019 <br />

1) fix corner case when reads distribution is unbalanced. <br />


Release 1.0 from the 28/03/2019 <br />

1) fix some integer conversion and prototype. <br />
2) fix a corner case when a rank recieve no mate from other rank. <br />

Release 1.0 from the 22/03/2019 <br />

1) fix the coordinates overlapped between 2 ranks in case we are not power of 2 dimension (again it's better if you are in power of 2 dimension) <br />

Release 1.0 from the 16/03/2019 <br />

1) cleaning of the code <br />
2) we forgot to write unmapped in previous release <br />
3) update of the overlapped coordinates between rank algorithm <br />
4) next step : a singularity definition file, more tests and profiling, upgrade of compression algorithm. <br />   


requirements
------------

tested with 

gcc > 4.8 (tested with 7.3) <br />
openmpi (tested with 2.2.1, 3, 4.0) <br />
openssl (tested with :1.0.2k, 1.1.0g) <br />
automake-1.15 <br />
autoconf-2.69 <br />

cmocka (optionnal and only for for unit testing) <br />
A SAM file of aligned paired reads, trimmed or not, and compliant with the SAM format. <br /> 
 
For small tests a laptop is sufficient. <br />
For real life test a HPC cluster with low latency network and parallel file system is a mandatory. <br />

How to compile:
--------------
aclocal <br />
autoconf <br />
automake --add-missing <br />
./configure --prefix /usr/local CC=path_to_mpicc <br />
make <br />
make install <br />
make clean <br />

How to test it
-------------

mpirun -n cpu_number mpiMD input_sam output_dir -q 0 -d 1000 -v 4 <br />

(it's better if cpu_number is a power of 2) <br />

options are: <br />

-q is for quality filtering <br />
-d is for optical distance duplicate <br />
-v is level of log verbose <br />
    0 is LOG_OFF  <br />
    1 is LOG_ERROR  <br />
    2 is LOG_WARNING  <br />
    3 is LOG_INFO (default) <br />
    4 is LOG_DEBUG  <br />
    5 is LOG_TRACE  <br />

How it works
------------

First the programm sort the reads by genome's coordinates and extract discordant and unmapped (end and mate) reads with the same technics described in mpiSORT. <br />

Second the programm mark the duplicates for each chromosome and discordant reads according to Picard Markduplicate method. The unmapped and unmapped mates are not marked. To limit memory overhead we build a distributed perfect hash table (see perfectHash.c for details) for fragment list and end list. This way the memory usage is under the memory usage of mpiSort.  <br />

Finally each chromosome is marked and compressed with bgzf and written down in the output folder.


Authors and contacts
--------------------

This program has been developed by<br />

Frederic Jarlier from Institut Curie and Firmin Martin from Paris Decartes University for marking of duplicates part<br />

and supervised by <br />

Philippe Hupe from Institut Curie <br />

Contacts: <br />

frederic.jarlier@curie.fr <br />
philippe.hupe@curie.fr <br />

