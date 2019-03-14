# mpiMarkDup
mpi and C based marking of duplicates for Next Generation Sequencing technology  

This is an open project initiated by Institut Curie HPC team and students from Paris Descartes University

The goal is to create a distributed and (near) real time version of the well known Picard MarkDuplicate from the Broad Institut. To do that we rely on MPI and C technology.   

The project is still under development and first results are pretty encouraging. We are close to 100% reproducible.

The code is a fork from the mpiSORT. We add an extra process (see markduplicate.c) in the loop to manage the marking of duplicates read.




