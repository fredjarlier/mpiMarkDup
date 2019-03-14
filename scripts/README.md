# Scripts

## compareSAMWithExternalTools.sh

### Description

Compare two SAM given in input. First the SAM files are convert to BAM
files by samtools view function and then BAM files are compared with
the function cmpbams of jvarkit. We suppose that the SAM files given
in input are coherent. By default the script suppose that samtools and
cmpbams.jar are located in script directory. If there is no -d option,
it supposes script directory is the current directory. 

It create a directory in this which containing samtools and
jvarkit.jar for temporary files.

### Options

* -g to indicate that input are gzipped
* -d directory where the script are. It is where
     compareSAMWithExternalTools.sh create/use a temporary directory
* -s indicates which samtools. By default the script suppose that samtools and cmpbams.jar are located in script directory.
* -j indicates which cmpbams.jar. By default the script suppose that
     samtools and cmpbams.jar are located in script directory.
* -h print help


### How to used it

    	compareSAMWithExternalTools.sh [-g] [-s samtools] [-j cmpbams.jar] -d DIRECTORY INPUT1 INPUT2

where DIRECTORY is the directory where the executable samtools and
jvarkit.jar are. INPUT1 and INPUT2 are the inputs.

## compareMpiMDMarkDuplicates.sh

### Description

Compare the output of mpiMD and MarkDuplicates. It takes a directory
containing directories of class of files test and use compare the
files. The script can create a very simple log file which contain the
output and the records that don’t match. It creates also a temporary
directory in the directory where it is.

### Options

* -h for help
* -s indicates the path of samtools
* -p indicates the path of picard.jar
* -j indicates the path of jvarkit.jar
* -S indicates the path of compareSAMWithExternalTools.sh
* -m indicates the path of mpiMD
* -l indicates that a log file is generated. You must indicate the name of the log file. Information are adding at the end of the file

### How to used it

    	compareMpiMDMarkDuplicates.sh [-h] [-s samtools] [-p picard.jar] [-j
    	jvarkit.jar] [-S compareSAMWithExternalTools.sh] [-m mpiMD]
    	[-l logfile] DIRECTORY

where DIRECTORY is the class test directory path.

## cmpWithMD.sh

### Description

Given a sam, compare the results between `markDuplicates` and `mpiMD`.

### Synopsis

`./cmpWithMD.sh SAM`

### User variables

User needs to set *user variables* at first use. They are

- PICARDTOOLS : picard tools jar path
- MPI_MD : `mpiMD` executable path
- NUM_PROCS : number of processes for `mpiMD`
- JVARKIT : `cmpbams` path from jvarkit
- SAMTOOLS : `samtools` executable paths

### Details

- This script compare `markDuplicates` and `mpiMD` by chromosome which appears in the header. Unmapped (only the case where read and its mate are unmapped) and discordants are remove from `markDuplicates` result. 
- This script give some stats in standard output.
- As output a directory is created, containing `markDuplicates` output, `mpiMD` output  and `jvarkit` output.

#### What this script don't do

- This script do not compare discordant.gz and unmapped.gz since they aren't handle yet.

## runCmpOnDir.sh

Given a directory of sam, compare for each sam the result of `markDuplicates` and `mpiMD`. There is a logfile as output and directories created by cmpWithMD.

### Synopsis

`./runCmpOnDir.sh DIR`

