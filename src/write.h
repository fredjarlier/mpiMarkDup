/*
   mpiSORT
   Copyright (C) 2016-2019 Institut Curie / Institut Pasteur

   mpiSORT is free software; you can redistribute it and/or
   modify it under the terms of the GNU Lesser General Public
   License as published by the Free Software Foundation; either
   version 2.1 of the License, or (at your option) any later version.

   mpiSORT is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser Public License
   along with mpiSORT.  If not, see <http://www.gnu.org/licenses/>.
*/

/*
   Module:
     write.h

   Authors:
    Authors:
    Frederic Jarlier, 	Institut Curie
	Nicolas Joly, 		Institut Pasteur
	Nicolas Fedy,		Institut Curie
	Leonor Sirotti,	 	Institut Curie
	Thomas Magalhaes,	Institut Curie
	Paul Paganiban,		Institut Curie
*/

#ifndef WRITE_H
#define WRITE_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>

#include <mpi.h>

#include "qksort.h"
#include "write_utils.h"
#include "parabitonicsort.h"

void writeSam(
		int rank,
		char* output_dir,
		char* header,
		size_t local_readNum,
		char* chrName,
		int total_num_proc,  //this is the number of proc in split communication
		MPI_Comm split_comm,
		int master_rank,
		MPI_Info finfo,
		int compression_level,
		size_t* new_offset_dest,
		size_t* new_offset_source,
		int* new_read_size,
		int* new_rank,
		int *original_rank_source_offset_phase1,
		char* data,
		size_t offset_data_in_file,
		size_t original_local_readNum,
		size_t final_local_readNum
		);


size_t init_offset_and_size_free_chr(size_t* offset, 
									 int* size, 
									 Read* data_chr, 
									 int local_readNum);

void read_data_for_writing(int rank, 
							int num_proc, 
							size_t local_readNum, 
							char *file_name,
							size_t *number_of_reads_by_procs, 
							size_t *buffs_by_procs, 
							char *** data,
							int *new_rank, 
							int *new_size, 
							size_t *new_offset, 
							MPI_File in, 
							MPI_Info finfo, 
							MPI_Comm COMM_WORLD);

void bruckWrite(MPI_Comm comm,
				int rank, 
				int num_proc,
				size_t local_readNum, 
				size_t* number_of_reads_by_procs, 
				int *new_rank,
				size_t *buffs_by_procs, 
				char*** data2,
				size_t *new_offset, 
				size_t*** data_offsets,
				int *new_size, 
				int ***data_size
	);

void bruckWrite2(
	    MPI_Comm comm,
		int rank,
		int num_proc,
		size_t local_readNum,
		size_t* number_of_reads_by_procs,
		int *new_rank,
		size_t *new_dest_offset,
		size_t ***data_dest_offsets,
		int *new_size,
		int ***data_size,
		size_t *new_source_offset,
		size_t ***data_source_offsets,
		int *dest_rank,
		int ***data_dest_rank
);

void bruckWrite3(
		MPI_Comm comm,
		int rank,
		int num_proc,
		size_t local_readNum,
		size_t *number_of_reads_by_procs,
		int *new_rank,
		size_t *new_source_offset,
		size_t ***data_source_offsets,
		int *new_dest_rank,
		int ***dest_rank,
		size_t *new_reads_coordinates,
		size_t ***read_coordinates,
		int *new_reads_size,
		int ***read_size,
		int *new_source_rank,
		int ***source_rank,
		size_t *new_dest_offset,
		size_t ***dest_offset
);

void bruckMarkdup(
    MPI_Comm comm,
    int rank,
    int num_proc,
    size_t local_readNum,
    size_t *number_of_reads_by_procs,
    int *new_rank,
    int *snd_mate_Lb,
    int ***rcv_mate_Lb,
    int *snd_mate_Materank,
    int ***rcv_mate_Materank,
    int *snd_mate_phredscore,    
    int ***rcv_mate_phredscore,
    size_t *snd_mate_indexAfterSort,
    size_t ***rcv_mate_indexAfterSort,
    size_t *snd_mate_unclippedCoordPos,
    size_t ***rcv_mate_unclippedCoordPos,
    size_t *snd_mate_coordPos,
    size_t ***rcv_mate_coordPos,
    size_t *snd_mate_coordMatePos,
    size_t ***rcv_mate_coordMatePos,
    size_t *snd_mate_fingerprint,
    size_t ***rcv_mate_fingerprint,
    size_t *snd_mate_fingerprint2,
    size_t ***rcv_mate_fingerprint2,
    unsigned int *snd_mate_valueFlag,
    unsigned int ***rcv_mate_valueFlag,
    unsigned int *snd_mate_pair_num,
    unsigned int ***rcv_mate_pair_num,
    unsigned int *snd_mate_orientation,
    unsigned int ***rcv_mate_orientation
);


void bruck_reads(	MPI_Comm comm, 
					int rank, 
					int num_proc, 
					size_t * buffs_by_procs, 
					char** data2);

void bruck_offsets( MPI_Comm comm,
					int rank, 
					int num_proc, 
					int local_readNum, 
					size_t* number_of_reads_by_procs, 
					size_t ** data_offsets, 
					int *new_rank, 
					size_t* new_offset);

void bruck_size(	MPI_Comm comm,
					int rank, 
					int num_proc, 
					size_t local_readNum, 
					size_t* number_of_reads_by_procs, 
					int ** data_size, 
					int *new_rank, 
					int *new_size);

void bruck_unsigned_int( 	MPI_Comm comm,
	 						int rank, 
							int num_proc, 
							size_t local_readNum, 
            				size_t *number_of_reads_by_procs, 
            				unsigned int **data_size, 
            				int *new_rank, 
            				unsigned int *new_size);


void writeSam_discordant_and_unmapped(int split_rank, 
									  char* output_dir, 
									  char* header, 
									  size_t local_readNum, 
									  char* chrName, 
									  Read* chr,
									  int num_proc, 
									  MPI_Comm split_comm, 
									  MPI_Info finfo, 
									  int compression_level, 
									  char *data,
									  size_t offset_data_in_file);



void writeSam_any_dim(
		int dimensions,
		int rank,
		char* output_dir,
		char* header,
		size_t local_readNum,
		size_t total_num_read,
		char* chrName,
		int total_num_proc,  //this is the number of proc in split communication
		MPI_Comm split_comm,
		int master_rank,
		MPI_Info finfo,
		int compression_level,
		size_t* new_offset_dest,
		size_t* new_offset_source,
		size_t* new_coordinates_sorted,
		int* new_read_size,
		int* new_rank,
		char *data,
		size_t start_offset_in_file
		);

#endif
