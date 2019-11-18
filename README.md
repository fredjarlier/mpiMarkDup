# mpiMarkDup

You are in the development Branch.

This branch is meant for development purposes. <br />

In this branch we add new features for the future master release. <br />

Here are the list of things in construction. <br />

1) Integrate the discordant reads in chromosom <br />

	16/07/2109

	We review the general algorithm to include the discordant fragments in the chromosom it belongs. <br />
	The discordant fragments are filtered and analyzed first and they are also kept in the chromosom buffer. <br />
	We use offset in the original sam file as unique identifier for the discordant fragments and their duplicates. <br />   
	Offests of all discordant duplicates are send to each jobs. <br />
	When we mark duplicates for chromosoms we use offsets in sam to get and mark the discordant duplicates. <br />

	So far the discordant duplicates are stored in a vector. <br />
	But we should consider a tree structure or a skip list to accelerate the search when marking discordant duplicates. <br />

	17/07/2019

	Test with sorting in any dimension cases and seems to work <br />

	Test of reproducibility with test case seems to work <br />


2) Build only one file with all chromosoms and unmapped at the end <br />

3) Accelerate the construction of the hashtable. Considering a list of prime number <br />

4) Modify the bruck with a its zero-copy version <br />

5) Generate a bam file instead of gz <br />

	we add magic EOF at the end of each file. <br />
	but we still have [main_samview] truncated file with samtools view <br />

6) Test with multiple libraries eg multiple RG <br />

requirements
------------

tested with 

gcc > 4.8 (tested with 7.3) <br />
openmpi (tested with 2.2.1, 3, 4.0, intel mpi 2019) <br />
openssl (tested with :1.0.2k, 1.1.0g) <br />

for distribution or compilation from source  <br />

automake-1.15 <br />
autoconf-2.69 <br />

cmocka (optionnal and only for for unit testing) <br />
A SAM file of aligned paired reads, trimmed or not, and compliant with the SAM format. <br /> 
 
For small tests a laptop is sufficient. <br />
For real life test a HPC cluster with low latency network and parallel file system is a mandatory. <br />

How to compile from source:
---------------------------------
aclocal <br />
autoconf <br />
automake --add-missing <br />
./configure --prefix /usr/local CC=path_to_mpicc <br />
make <br />
make install <br />
make clean <br />

How to distribute:
---------------------

See the master branch

How to test it
-------------

The input file is a sam file of aligned pair-end reads with a header. 

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

Finally each chromosome is marked and compressed with bgzf and written down in the output folder. <br />

We test the reproducibility by comparing both pipelines :mpiMD and mpiSORT + Picard (MarkDuplicate). <br />
We use the same number of cpu for each pipeline. So far We obtain 100% reproducibility. <br />

If the number of cpu differs the reproducibility is not garantee. Indeed tie cases are solved using the index of the read in the sorted file. This index can differ with the number of cpu. <br />

This problem does not impact the results in the downstream analysis. <br />

In conclusion when you test reproducibility always take the same number of cpu. <br />   


Authors and contacts
--------------------

This program has been developed by<br />

Frederic Jarlier from Institut Curie and Firmin Martin from Paris Decartes University for marking of duplicates part<br />

and supervised by <br />

Philippe Hupe from Institut Curie <br />

Contacts: <br />

frederic.jarlier@curie.fr <br />
philippe.hupe@curie.fr <br />

