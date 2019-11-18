# mpiMarkDup

MPI and C based program of marking duplicates for Next Generation Sequencing technology.<br /> 

The code is a fork from the mpiSORT project. We add a process to manage the marking of duplicate reads.

The goal is to create a distributed and (near) real time version of the well known Picard MarkDuplicate from the Broad Institute. To do that we rely on low level  technologies and parallel algorithms.<br />

The project is still under development so feel free to test it and report.

This is an open project initiated by Institut Curie HPC team and students from Paris Descartes University.

Release notes 
-------------

Release 1.0 from the 25/10/2019 <br />

1) Not a good idea to play with magic number bgzip complain. We come back to previous version

2) TODO: now the major bottle neck is the computation of the collision in perfect hash table. <br />
		 we shall compute collision only one time according to the biggest chromosom to mark <br /> 

3) TODO: problem with the option -q > 0 (it works with -q 0)

Release 1.0 from the 23/10/2019 <br />

1) Add the magic number at end of BGZF file. This way samtools doesn't complain <br />

example of usage with samtools
- samtools flagstat chr1.gz 
- samtools view -f 1024 chr1.gz 

Other tools shall complain (not tested). <br />

Release 1.0 from the 21/10/2019 <br />

1) Review of the search of mates after the bruck exchange of inter-process reads. <br />
     We now use the hash table to look for the corresponding pair of a read.  <br />


Release 1.0 from the 02/07/2019 <br />

1) Fix a bug when pairs overlap 2 ranks <br />

2) Fix a memory leak <br />

Release 1.0 from the 01/07/2019 <br />

1) Fix memory leaks <br />

2) Remove global variable COMM_WORLD <br />

Release 1.0 from the 27/06/2019 <br />

1) Fix a corner case. We check intra-node overlap before extra-node overlap in mpiSort.c. <br />

2) Fix a bug in sort_any_dim.c  <br />

Release 1.0 from the 19/06/2019 <br />

1) Fix a corner case <br />

Release 1.0 from the 13/06/2019 <br />

1) Replace All2all with a Bruck when exchange external frags <br />

2) Change some types in readInfo and mateInfo structures <br />

Release 1.0 from the 23/05/2019 <br />

1) Fix issues with openssl 1.0.2k-fips. <br />

Release 1.0 from the 15/05/2019 <br />

1) Cleaning up the code. <br />

Release 1.0 from the 07/05/2019 <br />

1) Fix reproducibility issue with 1 cpu. <br />

Release 1.0 from the 06/05/2019 <br />

1) Fix tie case (tested only with power of 2 cpu). <br />

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
openmpi (tested with 2.2.1, 3, 4.0, intel mpi 2019) <br />
openssl (tested with :1.0.2k, 1.1.0g) <br />

For compiling from code source or creating distribution you need <br />

automake-1.15 <br />
autoconf-2.69 <br />

cmocka (optionnal and only for unit testing) <br />
A SAM file of aligned paired reads, trimmed or not, and compliant with the SAM format. <br /> 
 
For small tests a laptop is sufficient. <br />
For real life test a HPC cluster with low latency network and parallel file system is a mandatory. <br />

How to compile from source:
---------------------------------

git clone the repo  <br />
cd mpiMarkDup  <br /> 
aclocal <br />
autoconf <br />
automake --add-missing <br />
./configure --prefix /usr/local CC=path_to_mpicc <br />
make <br />
make install <br />
make clean <br />

How to distribute from source:
-----------------------------------

git clone the repo  <br />
cd mpiMarkDup/src <br /> 
aclocal <br />
autoreconf --install <br />
./configure <br />
automake  <br />
make dist <br />

this create a tar.gz in mpiMarkDup/src that you can distribute <br />

How to compile from tar.gz:
--------------------------
tar -xvzf mpimd-1.0.tar.gz <br />
cd mpimd-1.0 <br />
./configure <br />
make <br />

How to test it
-------------

mpirun -n cpu_number mpiMD input_sam output_dir -q 0 -d 1000 -v 4 <br />

(it's better if cpu_number is a power of 2) <br />

options are: <br />

-q is for quality filtering (problem with q > 0 use -q 0)<br />
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

Finally each chromosome is marked and compressed with bgzf and written down in the output folder. <br />

We test the reproducibility by comparing both pipelines :mpiMD and mpiSORT + Picard (MarkDuplicate). <br />
We use the same number of cpu for each pipeline. So far We obtain 100% reproducibility. <br />

If the number of cpu differs the reproducibility is not garantee. Indeed tie cases are solved using the index of the read in the sorted file. This index can differ with the number of cpu. <br />

This problem does not impact the results in the downstream analysis. <br />

In conclusion when you test reproducibility always take the same number of cpu. <br />   


TODO 
----

List of thing that left to do before production level <br />

1) Build only one file with all chromosoms and unmapped <br />
 
2) Integrate the discordant reads in chromosom <br />

3) Test with multiple libraries eg multiple RG <br />

4) Accelerate the construction of the hashtable. Considering a list of prime number <br />

5) Modify the bruck with a its zero-copy version <br />

6) Generate a bam file instead of gz <br />



Authors and contacts
--------------------

This program has been developed by<br />

Frederic Jarlier from Institut Curie and Firmin Martin from Paris Decartes University for marking of duplicates part<br />

and supervised by <br />

Philippe Hupe from Institut Curie <br />

Contacts: <br />

frederic.jarlier@curie.fr <br />
philippe.hupe@curie.fr <br />

