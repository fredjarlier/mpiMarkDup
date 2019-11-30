/*
   mpiSORT
   Copyright (C) 2016-2017 Institut Curie / Institut Pasteur
   Copyright (C) 2018-2019 Institut Curie

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
     mpiSORT.c

   Authors:
    Frederic Jarlier,   Institut Curie
    Nicolas Joly,       Institut Pasteur
    Nicolas Fedy,       Institut Curie
    Leonor Sirotti,     Institut Curie
    Thomas Magalhaes,   Institut Curie
    Paul Paganiban,     Institut Curie
    Firmin Martin,		Paris Descartes University
*/

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#define _GNU_SOURCE

#include <sys/mman.h>
#include <sys/stat.h>
#include <assert.h>
#include <ctype.h>
#include <err.h>
#include <fcntl.h>
#include <libgen.h>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <stdlib.h>
#include <mpi.h>
#include <malloc.h>
#include "compat.h"

#include "mergeSort.h"
#include "parser.h"
#include "preWrite.h"
#include "write.h"
#include "sort_any_dim.h"
#include "mpiSort_utils.h"
#include "write_utils.h"
#include "qksort.h"
#include "parabitonicsort2.h"
#include "parabitonicsort3.h"

#include "log.h"


/*
 * Constant for Lustre striping
 * To adapt depending on the file server
 * STRIPING FACTOR is the numer of servers
 * where your files are striped
 *
 * STRIPING UNIT is the size of the stripes
 *
 * no need to change the values if you don't use
 * Lustre.
 *
 * Those parameters are harmless for other FS
 *
 */
#define STRIPING_FACTOR "12"
#define STRIPING_UNIT "4194304"   // 4 MB
//#define STRIPING_UNIT "268435456"   // 256 MB
//#define STRIPING_UNIT "536870912"   // 500 MB
//#define STRIPING_UNIT "1073741824"  // 1GB
//#define STRIPING_UNIT "1610612736"  // 1.5GB
//#define STRIPING_UNIT "2147483648"  // 2GB
//#define STRIPING_UNIT "2684354560"  // 2.5GB
//#define STRIPING_UNIT "3221225472"  // 3GB
//#define STRIPING_UNIT "3758096384"  // 3.5GB

/*
 * Constant for MPI IO
 * For description of the constant see here
 *
 * https://fs.hlrs.de/projects/craydoc/docs/books/S-2490-40/html-S-2490-40/chapter-sc4rx058-brbethke-paralleliowithmpi.html
 */

#define NB_PROC  "20" //number of threads for writing
#define CB_NODES "12" //number of servers for writing
#define CB_BLOCK_SIZE  "268435456" /* 256 MBytes - should match FS block size */
#define CB_BUFFER_SIZE  "536870912" /* multiple of the block size by the number of proc*/
#define DATA_SIEVING_READ "enable"

/*
 * Capacity constant no need to change it
 */
#define DEFAULT_MAX_SIZE 6000000000 //Default capacity for one process: 6G
#define DEFAULT_INBUF_SIZE  (1024*1024*1024)

/* Maximum chromosome number */
#define MAXNBCHR 256

static void usage(const char *);

int opticalDistance;

int main (int argc, char *argv[]) {

    char *x, *y, *z, *xbuf, *hbuf, *chrNames[MAXNBCHR];
    off_t hsiz;
    struct stat st;

    MPI_File mpi_filed;
    MPI_File mpi_file_split_comm;
    MPI_Offset fileSize, unmapped_start;
    int num_proc, rank, fd;
    int res, nbchr, i;
    int ierr, errorcode = MPI_ERR_OTHER;
    int g_rank, g_size;
    int split_rank, split_size; //after split communication we update the rank and the size
    int compression_level;
    //int paired = 0;
    char *file_name, *output_dir;
    char *header;
    unsigned int headerSize;
    unsigned char threshold;
    size_t input_file_size;
    size_t unmappedSize = 0;
    size_t discordantSize = 0;
    size_t *readNumberByChr = NULL, *localReadNumberByChr = NULL;
    double time_count;
    double toc;
    Read **reads;    
    MPI_Comm split_comm; //used to split communication when jobs have no reads to sort
    const char *sort_name;
    MPI_Info finfo;

    /* Set default values */
    compression_level = 3;
    parse_mode = MODE_OFFSET;
    sort_name = "coordinate";
    
    threshold = 0;
    opticalDistance = 0;
    int logSeverity = -1;

    /* Check command line */
    /* Check command line */
    while ((i = getopt(argc, argv, "c:hnpq:d:v:")) != -1) { //retourne le prochain caractère d'option dans argv qui concorde le caractère dans optstring
        switch (i) {
            case 'c': /* Compression level */
                compression_level = atoi(optarg); //convert stringArgument to an integer
                break;

            case 'h': /* Usage display */
                usage(basename(*argv)); //If successful, basename() returns a pointer to the final component of path
                return 0;

            case 'n':
                parse_mode = MODE_NAME;
                sort_name = "queryname";
                break;

            case 'q': /* Quality threshold */
                threshold = atoi(optarg);
                break;

            case 'd':
                opticalDistance = atoi(optarg);
                break;

            case 'v':
                logSeverity = atoi(optarg);
                break;

            default:
                usage(basename(*argv));
                return 1;
        }
    }

    switch (logSeverity) {
        case 0: 
            md_set_log_level(MD_LOG_OFF);
            break;

        case 1: 
            md_set_log_level(MD_LOG_ERROR);
            break;

        case 2: 
            md_set_log_level(MD_LOG_WARNING);
            break;

        case 3: 
            md_set_log_level(MD_LOG_INFO);
            break;

        case 4: 
            md_set_log_level(MD_LOG_DEBUG);
            break;

        case 5: 
            md_set_log_level(MD_LOG_TRACE);
            break;

        default:
            md_set_log_level(MD_LOG_INFO);
    }





    if (argc - optind != 2) {
        usage(basename(*argv));
        return 1;
    }

    file_name = argv[optind];
    output_dir = argv[optind + 1];

    /* Check arguments */
    res = access(file_name, F_OK | R_OK); //verify user permissions of a file

    if (res == -1) {
        err(1, "%s", file_name);
    }

    res = access(output_dir, F_OK | W_OK);

    if (res == -1) {
        err(1, "%s", output_dir);
    }

    /* MPI inits */
    res = MPI_Init(&argc, &argv);
    assert(res == MPI_SUCCESS); //Verify if there are nos errors
    res = MPI_Comm_rank(MPI_COMM_WORLD, &rank); //rank of process
    assert(res == MPI_SUCCESS);
    res = MPI_Comm_size(MPI_COMM_WORLD, &num_proc); //amount of process
    assert(res == MPI_SUCCESS);

    // We choose the seed in function of rank and time.
    // Note that if we use rank instead of rank + 1
    // rank 0 will always have the same seed.
    srand(time(NULL) * (rank+1)); //init the position of rand

    g_rank = rank;
    g_size = num_proc;

    /* Small summary */
    if (rank == 0) {
        fprintf(stderr, "Number of processes : %d\n", num_proc);
        fprintf(stderr, "Reads' quality threshold : %d\n", threshold);
        fprintf(stderr, "Compression Level is : %d\n", compression_level);
        fprintf(stderr, "SAM file to read : %s\n", file_name);
        fprintf(stderr, "Output directory : %s\n", output_dir);
    }

    /* Process input file */
    fd = open(file_name, O_RDONLY, 0666);
    assert(fd != -1);
    assert(fstat(fd, &st) != -1);
    xbuf = mmap(NULL, (size_t)st.st_size, PROT_READ, MAP_FILE | MAP_PRIVATE, fd, 0); //map or unmap files or devices into memory
    assert(xbuf != MAP_FAILED);

    /* Parse SAM header */
    memset(chrNames, 0, sizeof(chrNames)); //all elements of chrNames is init 0
    x = xbuf;
    nbchr = 0;

    //about chrNames values

    while (*x == '@') { //Every line that starts with "@" is part of the header
        y = strchr(x, '\n');  //renvoie un pointeur sur la première occurrence du caractère c dans la chaîne s
        z = x; //z pointe sur x, qui est dans un premier temps la première ligne du fichier
        x = y + 1;

        //On ne fait rien si dans la ligne on a pour les 3 premiers caractères :@SQ


        if (strncmp(z, "@SQ", 3) != 0) {
            continue;
        }

        /* Save reference names */
        y = strstr(z, "SN:");  //finds the first occurence of SN: in the String z
        assert(y != NULL);
        z = y + 3; //adress's position y + 3

        while (*z && !isspace((unsigned char)*z)) {//This function returns a non-zero value if c is a white-space character else, zero (false).
            z++;
        }

        chrNames[nbchr++] = strndup(y + 3, z - y - 3); //recover only (z-y-3) of (y + 3)
        assert(nbchr < MAXNBCHR - 2);
    }

    chrNames[nbchr++] = strdup(UNMAPPED);
    chrNames[nbchr++] = strdup(DISCORDANT);

    hsiz = x - xbuf;
    hbuf = strndup(xbuf, hsiz);

    if (rank == 0) {
        fprintf(stderr, "The size of the file is %zu bytes\n", (size_t)st.st_size);
        fprintf(stderr, "Header has %d+2 references\n", nbchr - 2);
    }

    int size_asp = asprintf(&header, "@HD\tVN:1.0\tSO:%s\n%s", sort_name, hbuf);
    assert(size_asp != -1);

    free(hbuf);
    assert(munmap(xbuf, (size_t)st.st_size) != -1);
    assert(close(fd) != -1);

    //task FIRST FINE TUNING FINFO FOR READING OPERATIONS


    MPI_Info_create(&finfo);
    /*
     * In this part you shall adjust the striping factor and unit according
     * to the underlying filesystem.
     * Harmless for other file system.
     */
    MPI_Info_set(finfo, "striping_factor", STRIPING_FACTOR);
    MPI_Info_set(finfo, "striping_unit", STRIPING_UNIT); //2G striping
    MPI_Info_set(finfo, "ind_rd_buffer_size", STRIPING_UNIT); //2gb buffer
    MPI_Info_set(finfo, "romio_ds_read", DATA_SIEVING_READ);

    /*
     * for collective reading and writing
     * should be adapted too and tested according to the file system
     * Harmless for other file system.
     */
    MPI_Info_set(finfo, "nb_proc", NB_PROC);
    MPI_Info_set(finfo, "cb_nodes", CB_NODES);
    MPI_Info_set(finfo, "cb_block_size", CB_BLOCK_SIZE);
    MPI_Info_set(finfo, "cb_buffer_size", CB_BUFFER_SIZE);


    //we open the input file
    ierr = MPI_File_open(MPI_COMM_WORLD, file_name,  MPI_MODE_RDONLY, finfo, &mpi_filed);

    //assert(in != -1);
    if (ierr) {
        if (rank == 0) {
            fprintf(stderr, "%s: Failed to open file in process 0 %s\n", argv[0], argv[1]);
        }

        MPI_Abort(MPI_COMM_WORLD, errorcode);
        exit(2);
    }

    ierr = MPI_File_get_size(mpi_filed, &fileSize);
    assert(ierr == MPI_SUCCESS);
    input_file_size = (long long)fileSize;

    /* Get chunk offset and size */
    headerSize = unmappedSize = discordantSize = strlen(header);

    //We place file offset of each process to the begining of one read's line
    size_t *goff = (size_t *)calloc((size_t)(num_proc + 1), sizeof(size_t));
    init_goff(mpi_filed, hsiz, input_file_size, num_proc, goff);

    //NOW WE WILL PARSE
    size_t j = 0;
    size_t poffset = goff[rank]; //Current offset in file sam

    //nbchr because we add the discordant reads in the structure
    reads = (Read **)malloc((nbchr) * sizeof(Read)); //We allocate a linked list of struct for each Chromosome (last chr = unmapped reads)
    readNumberByChr = (size_t *)malloc((nbchr) * sizeof(size_t)); //Array with the number of reads found in each chromosome
    localReadNumberByChr = (size_t *)malloc((nbchr) * sizeof(size_t)); //Array with the number of reads found in each chromosome
    Read **anchor = (Read **)malloc((nbchr) * sizeof(Read)); //Pointer on the first read of each chromosome

    //Init first read
    for (i = 0; i < (nbchr); i++) {
        reads[i] = malloc(sizeof(Read));
        reads[i]->coord = 0;
        anchor[i] = reads[i];
        readNumberByChr[i] = 0;
    }

    toc = MPI_Wtime();

    char *local_data_tmp = malloc(1024 * 1024);
    char *local_data = (char *)malloc(((goff[rank + 1] - poffset) + 1) * sizeof(char));
    size_t size_tmp = goff[rank + 1] - poffset;
    local_data[goff[rank + 1] - poffset] = 0;
    char *q = local_data;

    //We read the file sam and parse
    while (poffset < goff[rank + 1]) {

        size_t size_to_read = 0;

        if ( (goff[rank + 1] - poffset) < DEFAULT_INBUF_SIZE ) {
            size_to_read = goff[rank + 1] - poffset;

        } else {
            size_to_read = DEFAULT_INBUF_SIZE;
        }

        // we load the buffer
        //hold temporary size of SAM
        //due to limitation in MPI_File_read_at
        local_data_tmp = (char *)realloc(local_data_tmp, (size_to_read + 1) * sizeof(char));
        local_data_tmp[size_to_read] = 0;

        // Original reading part is before 18/09/2015
        MPI_File_read_at(mpi_filed, (MPI_Offset)poffset, local_data_tmp, size_to_read, MPI_CHAR, MPI_STATUS_IGNORE);
        size_t local_offset = 0;
        assert(strlen(local_data_tmp) == size_to_read);

        //we look where is the last line read for updating next poffset
        size_t offset_last_line = size_to_read - 1;

        size_t extra_char = 0;

        while (local_data_tmp[offset_last_line] != '\n') {
            offset_last_line -- ;
            extra_char++;
        }

        local_data_tmp[size_to_read - extra_char] = 0;
        size_t local_data_tmp_sz = strlen(local_data_tmp);

        //If it s the last line of file, we place a last '\n' for the function tokenizer
        if (rank == num_proc - 1 && ((poffset + size_to_read) == goff[num_proc])) {
            local_data_tmp[offset_last_line] = '\n';
        }

        //Now we parse Read in local_data

        parser_paired(local_data_tmp, rank, poffset, threshold, nbchr, &readNumberByChr, chrNames, &reads);

        //now we copy local_data_tmp in local_data
        char *p = local_data_tmp;
        size_t pos = 0;

        while (*p && (pos < local_data_tmp_sz)) {
            *q = *p;
            p++;
            q++;
            pos++;
        }

        //we go to the next line
        poffset += (offset_last_line + 1);
        local_offset += (offset_last_line + 1);

    }

    assert(size_tmp == strlen(local_data));

    //fprintf(stderr, "%d (%.2lf)::::: *** FINISH PARSING FILE ***\n", rank, MPI_Wtime() - toc);
    md_log_rank_debug(rank, "Finished to parsing file in %2lf \n", MPI_Wtime() - toc);

    if (local_data_tmp) {
        free(local_data_tmp);
    }

    //malloc_trim(0);

    MPI_Barrier(MPI_COMM_WORLD);



    //We set attribute next of the last read and go back to first read of each chromosome

    for (i = 0; i < nbchr; i++) {
        reads[i]->next = NULL;
        reads[i] = anchor[i];
    }

    free(anchor);

    //We count how many reads we found
    size_t nb_reads_total = 0, nb_reads_global = 0;
    int j2 = 0;

    for (j2 = 0; j2 < nbchr; j2++) {
        nb_reads_total += readNumberByChr[j2];
    }

    //nb_reads_total per rank
    //nb_reads_global is sum reads of all rank

    MPI_Allreduce(&nb_reads_total, &nb_reads_global, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD); //Sums all of the total reads for each rank into the variable nb_reads_global
    md_log_rank_debug(0, "total read to sort  = %zu \n", nb_reads_global);
    md_log_rank_debug(0, "number of chromosome  = %d + discordants + unmapped \n\n", nbchr - 2);
    /*
     *  We write the mapped reads in a file named chrX.bam
     *  We loop by chromosoms.
     */

    MPI_Barrier(MPI_COMM_WORLD);

    for (i = 0; i < nbchr ; i++) {

        if ( i == nbchr - 2){ 
            //If unmapped, we compress these reads in one folder without sort
            MPI_File mpi_file_split_comm2;
            double time_count;

            size_t total_reads = 0;
            MPI_Allreduce(&readNumberByChr[i], &total_reads, 1, MPI_LONG_LONG_INT, MPI_SUM, MPI_COMM_WORLD);//Sum of all the reads per chromosomes

            //if ((rank == 0){
            //    fprintf(stderr, "rank %d :::: total read to sort for unmapped = %zu \n", rank, total_reads);
            //}
            
            md_log_rank_debug(0, "total read to sort for unmapped = %zu\n\n", total_reads);

            MPI_Barrier(MPI_COMM_WORLD);

            if (total_reads == 0) {
                // nothing to sort for unmapped
                // maybe write an empty bam file
            } else {
                int i1, i2;
                size_t *localReadsNum_rank0 = (size_t *)malloc(num_proc * sizeof(size_t));
                localReadsNum_rank0[0] = 0;
                int file_pointer_to_free = 0;
                int split_comm_to_free = 0;
                //we build a vector with rank job
                int val_tmp1 = 0;
                int val_tmp2 = 0;
                int chosen_rank = 0;
                // the color tells in what communicator the rank pertain
                // color = 0 will be the new communicator color
                // otherwise the color is 1
                int *color_vec_to_send =  (int *)malloc(num_proc * sizeof(int));
                // the key value tell the order in the new communicator
                int *key_vec_to_send =  (int *)malloc(num_proc * sizeof(int));

                //rank 0 gather the vector
                MPI_Allgather(&readNumberByChr[i], 1, MPI_LONG_LONG_INT, localReadsNum_rank0, 1, MPI_LONG_LONG_INT, MPI_COMM_WORLD);//Now each rank knows the read number by chr of other ranks
                MPI_Barrier(MPI_COMM_WORLD);

                if (rank == 0) {
                    //we must chose the first rank with reads to sort
                    i1 = 0;

                    while (localReadsNum_rank0[i1] == 0) {
                        chosen_rank++;
                        i1++;
                    }
                }

                //we broadcast the chosen rank
                //task: replace the broadcast with a sendrecieve
                MPI_Bcast( &chosen_rank, 1, MPI_INT, 0, MPI_COMM_WORLD);
                MPI_Barrier(MPI_COMM_WORLD);

                //we must chose which rank is going to split the communication
                if (((rank == chosen_rank) || rank == 0) && (chosen_rank != 0)) {
                    //the rank 0 will recieve the key_vec_to_send and colorvec_to_send
                    //first we exchange the size o
                    if (rank == chosen_rank) {
                        header = (char *)malloc((headerSize + 1) * sizeof(char));
                        MPI_Recv(header, headerSize + 1, MPI_CHAR, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    }

                    if (rank == 0) {
                        MPI_Send(header, headerSize + 1, MPI_CHAR, chosen_rank,  0, MPI_COMM_WORLD);
                    }

                } else {
                    //we do nothing here
                }

                if (rank == chosen_rank) {

                    int counter = 0;

                    //we compute the number of 0 in the localReadsNum_vec
                    for (i1 = 0; i1 < num_proc; i1++) {
                        if (localReadsNum_rank0[i1] == 0) {
                            counter++;
                        }
                    }

                    // if no jobs without reads we do nothing
                    if ( counter == 0 ) {
                        // nothing to do we associate split_comm with
                        split_comm = MPI_COMM_WORLD;

                        for (i2 = 0; i2 < num_proc; i2++) {

                            if (localReadsNum_rank0[i2] == 0) {
                                color_vec_to_send[i2] = 1;
                                key_vec_to_send[i2] = val_tmp2;
                                val_tmp2++;

                            } else {
                                color_vec_to_send[i2] = 0;
                                key_vec_to_send[i2] = val_tmp1;
                                val_tmp1++;
                            }
                        }

                    } else {
                        // now we compute the color according to
                        // the number of reads to sort
                        for (i2 = 0; i2 < num_proc; i2++) {
                            if (localReadsNum_rank0[i2] == 0) {
                                color_vec_to_send[i2] = 1;
                                key_vec_to_send[i2] = val_tmp2;
                                val_tmp2++;

                            } else {
                                color_vec_to_send[i2] = 0;
                                key_vec_to_send[i2] = val_tmp1;
                                val_tmp1++;
                            }
                        } // end for loop
                    }// end if
                }// end if (rank == chosen_rank)

                MPI_Barrier(MPI_COMM_WORLD);
                // we scatter the key and color vector
                // we create key and color variable for each job
                int local_color = 0;
                int local_key = 0;
                // we scatter the color and key
                MPI_Scatter( color_vec_to_send, 1, MPI_INT, &local_color, 1, MPI_INT, chosen_rank, MPI_COMM_WORLD);
                MPI_Scatter( key_vec_to_send, 1, MPI_INT, &local_key, 1, MPI_INT, chosen_rank, MPI_COMM_WORLD);

                free(color_vec_to_send);
                free(key_vec_to_send);

                // we create a communicator
                // we group all communicator
                // with color of zero
                if (local_color == 0) {

                    MPI_Comm_split( MPI_COMM_WORLD, local_color, local_key, &split_comm);
                    ierr = MPI_File_open(split_comm, file_name,  MPI_MODE_RDONLY, finfo, &mpi_file_split_comm2);
                    //we ask to liberate file pointer
                    file_pointer_to_free = 1;
                    //we ask to liberate the split_comm
                    split_comm_to_free = 1;

                } else {
                    MPI_Comm_split( MPI_COMM_WORLD, MPI_UNDEFINED, local_key, &split_comm);
                    mpi_file_split_comm2 = mpi_filed;
                }

                //now we change the rank in the reads structure
                if (local_color == 0) {
                    MPI_Comm_rank(split_comm, &split_rank);
                    MPI_Comm_size(split_comm, &split_size);

                    g_rank = split_rank;
                    g_size = split_size;

                    reads[i] = reads[i]->next;
                    localReadNumberByChr[i] = readNumberByChr[i];
                    
                    unmapped_start = startOffset(g_rank,
                                                g_size,
                                                unmappedSize,
                                                headerSize,
                                                split_comm
                                                );

                    if (!unmapped_start) {
                        fprintf(stderr, "No header was defined for unmapped. \n Shutting down.\n");
                        MPI_Finalize();
                        return 0;
                    }

                    time_count = MPI_Wtime();
                    MPI_Barrier(MPI_COMM_WORLD);

                    writeSam_discordant_and_unmapped(
                        split_rank,
                        output_dir,
                        header,
                        localReadNumberByChr[i],
                        chrNames[i],
                        reads[i],
                        split_size,
                        split_comm,
                        finfo,
                        compression_level,
                        local_data,
                        goff[rank]);

                        //if (split_rank == chosen_rank) {
                        //    fprintf(stderr, "rank %d :::::[MPISORT] Time to write chromosom %s ,  %f seconds \n\n\n", split_rank, chrNames[nbchr - s], MPI_Wtime() - time_count);
                        //}
                    md_log_rank_debug(chosen_rank, "[MPISORT] Time to write chromosom %s ,  %f seconds \n\n\n", chrNames[nbchr - i], MPI_Wtime() - time_count);

                    
                    while (reads[i]->next != NULL) {
                        Read *tmp_chr = reads[i];
                        reads[i] = reads[i]->next;
                        free(tmp_chr);
                    }

                    free(localReadsNum_rank0);

                } else {
                    // we do nothing
                }

                //we put a barrier before freeing pointers
                MPI_Barrier(MPI_COMM_WORLD);
                //we free the file pointer

                if  (file_pointer_to_free) {
                    MPI_File_close(&mpi_file_split_comm2);
                }

                //we free the split_comm
                if (split_comm_to_free) {
                    MPI_Comm_free(&split_comm);
                }

                split_comm_to_free = 0;
                file_pointer_to_free = 0;

            }
        }
        md_log_rank_info(0, "start to handle chromosome %s ...\n", chrNames[i]);

        /*
         * First Part of the algorithm
         *
         * In this part we elected a rank which is the first rank
         * to have reads to sort.
         *
         * Once elected a rank, we plit the communicator according to
         * wether the rank has reads to sort for this chromosom.
         *
         * The new communicator is COMM_WORLD.
         *
         * If all jobs have reads to sort no need to split the communicator and then
         * COMM_WORLD = MPI_COMM_WORLD
         *
         */

        int i1, i2;
        size_t localReadsNum_rank0[num_proc];
        localReadsNum_rank0[0] = 0;
        int file_pointer_to_free = 0;
        int split_comm_to_free = 0;
        //we build a vector with rank job
        int val_tmp1 = 0;
        int val_tmp2 = 0;
        int chosen_rank = 0; //needed to tell what rank is going to compute the color and key
        int chosen_split_rank = 0; //the rank that collect data once the communication splitted normally this rank is 0

        // the color tells in what communicator the rank pertain
        // color = 0 will be the new communicator color
        // otherwise the color is 1
        // the key value tell the order in the new communicator
        int *color_vec_to_send  =  (int *)malloc(num_proc * sizeof(int));
        int *key_vec_to_send    =  (int *)malloc(num_proc * sizeof(int));

        // first we test if the there's reads to sort
        // rank 0 recieve the sum of all the reads count
        size_t total_reads_by_chr = 0;
        MPI_Allreduce(&readNumberByChr[i], &total_reads_by_chr, 1, MPI_LONG_LONG_INT, MPI_SUM, MPI_COMM_WORLD);

        if (total_reads_by_chr == 0) {
            free(color_vec_to_send);
            free(key_vec_to_send);
            continue;    //pass to next chromosome
        }

        //rank 0 gather the vector
        MPI_Allgather(&readNumberByChr[i], 1, MPI_LONG_LONG_INT, localReadsNum_rank0, 1, MPI_LONG_LONG_INT, MPI_COMM_WORLD); //each process knows numbers chr's reads of other rank

        if (rank == 0) {
            //the rank 0 chose the first rank with reads to sort
            i1 = 0;

            while ((localReadsNum_rank0[i1] == 0) && (i1 < num_proc)) {
                chosen_rank++;
                i1++;
            }

            //fprintf(stderr, "rank %d :::: Elected rank = %d \n", rank, chosen_rank);
        }

        md_log_debug("Elected rank = %d\n", chosen_rank);

        //we broadcast the chosen rank
        //task: replace the broadcast with a sendrecieve
        MPI_Bcast( &chosen_rank, 1, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Barrier(MPI_COMM_WORLD);

        if (((rank == chosen_rank) || rank == 0) && (chosen_rank != 0)) {

            //first we exchange the size o
            if (rank == chosen_rank) {
                header = malloc((headerSize + 1) * sizeof(char));
                header[headerSize] = '\0';
                MPI_Recv(header, headerSize + 1, MPI_CHAR, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }

            if (rank == 0) {
                MPI_Send(header, headerSize + 1, MPI_CHAR, chosen_rank,  0, MPI_COMM_WORLD);
            }

        } else {
            //we do nothing here
        }

        MPI_Barrier(MPI_COMM_WORLD);

        if (rank == chosen_rank) {
            int counter = 0;

            //we compute the number of 0 in the localReadsNum_vec
            for (i1 = 0; i1 < num_proc; i1++) {
                if (localReadsNum_rank0[i1] == 0) {
                    counter++;
                }
            }

            // if no jobs without reads we do nothing
            if ( counter == 0 ) {
                // nothing to do we associate split_comm with
                //fprintf(stderr, "rank %d ::::[MPISORT] we don't split the rank \n", rank);
                md_log_rank_debug(chosen_rank, "[MPISORT] we don't split the rank \n");
                split_comm = MPI_COMM_WORLD;

                for (i2 = 0; i2 < num_proc; i2++) {
                    if (localReadsNum_rank0[i2] == 0) {
                        color_vec_to_send[i2] = 1;
                        key_vec_to_send[i2] = val_tmp2;
                        val_tmp2++;

                    } else {
                        color_vec_to_send[i2] = 0;
                        key_vec_to_send[i2] = val_tmp1;
                        val_tmp1++;
                    }
                }

            } else {
                // now we compute the color according to
                // the number of reads to sort
                //fprintf(stderr, "rank %d ::::[MPISORT] we split the rank \n", rank);
                md_log_rank_debug(rank, "[MPISORT] we split the rank\n");

                for (i2 = 0; i2 < num_proc; i2++) {
                    if (localReadsNum_rank0[i2] == 0) {
                        color_vec_to_send[i2] = 1;
                        key_vec_to_send[i2] = val_tmp2;
                        val_tmp2++;

                    } else {
                        color_vec_to_send[i2] = 0;
                        key_vec_to_send[i2] = val_tmp1;
                        val_tmp1++;
                    }
                } // end for loop
            }// end if
        }// end if (rank == plit_rank)

        MPI_Barrier(MPI_COMM_WORLD);
        //we create key and color variable for each job
        int local_color = 0;
        int local_key = 0;
        // rank 0 scatter the color and the key vector
        MPI_Scatter( color_vec_to_send, 1, MPI_INT, &local_color, 1, MPI_INT, chosen_rank, MPI_COMM_WORLD);
        MPI_Scatter( key_vec_to_send, 1, MPI_INT, &local_key, 1, MPI_INT, chosen_rank, MPI_COMM_WORLD);
        MPI_Barrier(MPI_COMM_WORLD);

        free(color_vec_to_send);
        free(key_vec_to_send);


        // now we create a communicator
        // we group all communicator
        // with color of zero
        if (local_color == 0) {
            MPI_Comm_split( MPI_COMM_WORLD, local_color, local_key, &split_comm);
            ierr = MPI_File_open(split_comm, file_name,  MPI_MODE_RDONLY, finfo, &mpi_file_split_comm);
            //we ask to liberate file pointer
            file_pointer_to_free = 1;
            //we ask to liberate the split_comm
            split_comm_to_free = 1;

        } else {
            MPI_Comm_split( MPI_COMM_WORLD, MPI_UNDEFINED, local_key, &split_comm);
            mpi_file_split_comm = mpi_filed;
        }

        //now we change the rank in the reads structure
        if (local_color == 0) {

            MPI_Comm_rank(split_comm, &split_rank);
            MPI_Comm_size(split_comm, &split_size);



            //we update g_rank
            g_rank = split_rank;
            g_size = split_size;

        } else {
            g_rank = split_rank = rank;
            g_size = split_size = num_proc;
        }

        //REPRISE DU CODE IMPORTANT

        localReadNumberByChr[i] = readNumberByChr[i];
        MPI_Barrier(MPI_COMM_WORLD);

        if ((local_color == 0) && (i < nbchr)) {
            if (i == nbchr - 2){ //If unmapped                    
                continue;
            }

            /*
             * Second part of the algorithm
             *
             * First we load coordinates, offset sources, and read size in vector
             *
             * Then we sort the coordinates of the reads
             * with a bitonic sorter
             *
             * Then according to the reads coordinates we reoder the offset sources, and size
             * this is done thanks to the index of the sorting.
             *
             * Afterward we compute the offsets of the reads in
             * the destination file.
             *
             * Finally we dispatch the information to all ranks
             * in the communicator for the next step.
            */


            //we do a local merge sort
            assert(reads[i]);
            if (reads[i] && reads[i]->next && reads[i]->next->next) {
                //Pourquoi on ne stocke pas le nouveau read trié ? 
                mergeSort(reads[i], readNumberByChr[i]);   //TRI POUR CHAQUE CHROMOSOME
            }

            size_t local_readNum = localReadNumberByChr[i];

            reads[i] = reads[i]->next;
            assert(reads[i]);

            //first we compute the dimension of the parabitonic sort
            // dimension is the number of processors where we
            // perform the bitonic sort
            // int dimensions = (int)(log2(num_processes));

            // find next ( must be greater) power, and go one back
            int dimensions = 1;

            //Dimensions is true number of rank who works
            //Split_size is total rank 

            while (dimensions <= split_size) { //When the nb of ranks is of the power of 2
                dimensions <<= 1;
            }

            dimensions >>= 1;


            // we get the maximum number of reads among
            // all the workers

            /*
             * Here we split the programm in 2 cases
             *
             * 1) The first case of split_size is a power of 2 (the best case)
             *      this case is the simpliest we don't have extra communication to dispatch the read
             *      envenly between the jobs
             *
             * 2) The split_size is not a power of 2 (the worst case)
             *      well in this case we shall dispatch the jobs between jobs evenly.
             *
             */

            //if (split_rank == chosen_split_rank) {

            //    fprintf(stderr, "Rank %d :::::[MPISORT] Dimensions for bitonic = %d \n", split_rank, dimensions);
            //    fprintf(stderr, "Rank %d :::::[MPISORT] Split size               = %d \n", split_rank, split_size);

            //}
            md_log_rank_debug(chosen_split_rank, "[MPISORT] Dimensions for bitonic = %d \n", dimensions);
            md_log_rank_debug(chosen_split_rank, "[MPISORT] Split size             = %d \n", split_size);

            //we test the computed dimension
            if ( dimensions == split_size ) {

                size_t max_num_read = 0;
                MPI_Allreduce(&localReadNumberByChr[i], &max_num_read, 1, MPI_LONG_LONG_INT, MPI_MAX, split_comm);//We take the max read number out of all the ranks, so we can allocate the right amount later

                // if the dimension == split_size
                MPI_Barrier(split_comm);

                size_t first_local_readNum = local_readNum; //(= readNumberByChr[i])

                /*
                 * Vector creation and allocation
                 */
                local_readNum = max_num_read;

                time_count = MPI_Wtime();

                //configure the size of each memory allocation (local_readNum = max number of read of rank)

                size_t *local_reads_coordinates_unsorted    = calloc(local_readNum, sizeof(size_t));
                size_t *local_reads_coordinates_sorted      = calloc(local_readNum, sizeof(size_t));
                size_t *local_offset_source_unsorted        = calloc(local_readNum, sizeof(size_t));
                size_t *local_offset_source_sorted          = calloc(local_readNum, sizeof(size_t));
                int *local_dest_rank_sorted                 = calloc(local_readNum, sizeof(int));
                int *local_reads_sizes_unsorted             = calloc(local_readNum, sizeof(int));
                int *local_reads_sizes_sorted               = calloc(local_readNum, sizeof(int));
                int *local_source_rank_unsorted             = calloc(local_readNum, sizeof(int));
                int *local_source_rank_sorted               = calloc(local_readNum, sizeof(int));

                //if (split_rank == chosen_split_rank) {
                //    fprintf(stderr,   "rank %d :::::[MPISORT][MALLOC 1] time spent = %f s\n", split_rank, MPI_Wtime() - time_count);
                //}
                md_log_rank_debug(chosen_split_rank, "[MPISORT][MALLOC 1] time spent = %f s\n", split_rank, MPI_Wtime() - time_count);

                //We need to reduce variables here for optimisation
				//Here local means for 1 chromosome

                local_reads_coordinates_unsorted[0] = 0;
                local_reads_coordinates_sorted[0]   = 0;
                local_dest_rank_sorted[0]           = 0;
                local_reads_sizes_unsorted[0]       = 0;
                local_reads_sizes_sorted[0]         = 0;
                local_source_rank_unsorted[0]       = 0;
                local_source_rank_sorted[0]         = 0;
                local_offset_source_unsorted[0]     = 0;
                local_offset_source_sorted[0]       = 0;

                //those vectors are the same that  local_..._sorted but without zero padding
                size_t *local_reads_coordinates_sorted_trimmed = NULL;
                size_t *local_offset_source_sorted_trimmed = NULL;
                size_t *local_offset_dest_sorted_trimmed = NULL;
                //vectors used in the bruck just after the parabitonic sort
                size_t *local_reads_coordinates_sorted_trimmed_for_bruck = NULL;
                size_t *local_offset_source_sorted_trimmed_for_bruck = NULL;
                size_t *local_offset_dest_sorted_trimmed_for_bruck = NULL;
                
                int *local_dest_rank_sorted_trimmed_for_bruck = NULL;
                int *local_reads_sizes_sorted_trimmed_for_bruck = NULL; 
                int *local_source_rank_sorted_trimmed_for_bruck = NULL;
                int *local_source_rank_sorted_trimmed = NULL;
                int *local_dest_rank_sorted_trimmed = NULL;
                int *local_reads_sizes_sorted_trimmed = NULL;
                //task Init offset and size for source - free chr
                // from mpiSort_utils.c

                //initiate values of vectors from reads[i]
                get_coordinates_and_offset_source_and_size_and_free_reads(
                    split_rank,
                    local_source_rank_unsorted,// this value is set to split_rank in the function
                    local_reads_coordinates_unsorted,//coordinates pulled from the linked chain reads[]
                    local_offset_source_unsorted,//the offset/start of the read in the source file that was parsed (start of the line)
                    local_reads_sizes_unsorted,//the line_size (we can get the end of the line with it)
                    reads[i],
                    first_local_readNum
                );

                //init indices for qksort
                size_t *coord_index = (size_t *)malloc(local_readNum * sizeof(size_t));

                for (j = 0; j < local_readNum; j++) {
                    coord_index[j] = j;
                }

                //To start we sort locally the reads coordinates.
                //this is to facilitate the bitonic sorting
                //if the local coordinates to sort are to big we could get rid of
                //this step.
                time_count = MPI_Wtime();

                //base_arr2 = local_reads_coordinates_unsorted;
				//Seems like base_arr2 is never used again, happens in other files aswell

                //A QUOI SERT qksort

                qksort(coord_index, local_reads_coordinates_unsorted, local_readNum, 
                    0, local_readNum - 1, compare_size_t);
                //It looks like qksort inverts the values in coord_index

                //if (split_rank == chosen_split_rank) {
                //    fprintf(stderr,   "rank %d :::::[MPISORT][LOCAL SORT] time spent = %f s\n", split_rank, MPI_Wtime() - time_count);
                //}
                md_log_rank_debug(chosen_split_rank, "[MPISORT][LOCAL SORT] time spent = %f s\n", split_rank, MPI_Wtime() - time_count);

                //We index data
                //For all of these, the first value of the sorted vectors takes the last value of the unsorted ones
                for (j = 0; j < local_readNum; j++) { 
                    local_reads_coordinates_sorted[j]           = local_reads_coordinates_unsorted[coord_index[j]];
                    local_source_rank_sorted[j]                 = local_source_rank_unsorted[coord_index[j]];
                    local_reads_sizes_sorted[j]                 = local_reads_sizes_unsorted[coord_index[j]];
                    local_offset_source_sorted[j]               = local_offset_source_unsorted[coord_index[j]];
                    local_dest_rank_sorted[j]                   = rank; //will be updated after sorting the coordinates
                }

                //We free memorys and for next of the loop we do same algo

                free(coord_index);                      //ok
                free(local_source_rank_unsorted);       //ok
                free(local_reads_coordinates_unsorted); //ok
                free(local_reads_sizes_unsorted);       //ok
                free(local_offset_source_unsorted);     //ok

                // we need the total number of reads.
                size_t total_num_read = 0;
                MPI_Allreduce(&localReadNumberByChr[i], &total_num_read, 1, MPI_LONG_LONG_INT, MPI_SUM, split_comm);

                /*
                 *
                 * In this section the number of bitonic dimension
                 * is equal to the split size. (if dimensions == split_size)
                 *
                 * In this case there are less communication in preparation
                 * of the sorting.
                 *
                 * We use the parabitonic version 2.
                 */

                //we calll the bitonic

                time_count = MPI_Wtime();

                ParallelBitonicSort2(
                    split_comm,
                    split_rank,
                    dimensions,
                    local_reads_coordinates_sorted,
                    local_reads_sizes_sorted,
                    local_source_rank_sorted,
                    local_offset_source_sorted,
                    local_dest_rank_sorted,
                    max_num_read
                );

                //if (split_rank == chosen_split_rank) {
                //    fprintf(stderr, "rank %d :::::[MPISORT][BITONIC 2] time spent = %f s\n", split_rank, MPI_Wtime() - time_count);
                //}

                md_log_rank_debug(chosen_split_rank, "[MPISORT][BITONIC 2] time spent = %f s\n", split_rank, MPI_Wtime() - time_count);

                size_t k1;
                int k2; 

                for (k1 = 1; k1 < max_num_read; k1++) {
                    assert(local_reads_coordinates_sorted[k1 - 1] <= local_reads_coordinates_sorted[k1]);//Checking if it's sorted
                    local_dest_rank_sorted[k1] = split_rank;
                }

                size_t *local_offset_dest_sorted = malloc(max_num_read * sizeof(size_t));
                size_t last_local_offset = 0;

                // We compute the local_dest_offsets_sorted
                size_t local_total_offset = 0;

                for (k1 = 0; k1 <  max_num_read; k1++) {
                    local_offset_dest_sorted[k1] = local_reads_sizes_sorted[k1];
                    local_total_offset += local_reads_sizes_sorted[k1]; // this variable is never used again
                }

                //we make the cumulative sum of all offsets
                for (k1 = 1; k1 < max_num_read; k1++) {
                    local_offset_dest_sorted[k1] = local_offset_dest_sorted[k1 - 1] + local_offset_dest_sorted[k1];
                }

                //we exchange the last destination offset
                last_local_offset = local_offset_dest_sorted[max_num_read - 1];

                //number of block to send
              
                MPI_Offset *y  = calloc(split_size, sizeof(MPI_Offset));
                MPI_Offset *y2 = calloc(split_size + 1, sizeof(MPI_Offset));

                //we wait all processors

                //Le rang 0 récupère les last_local_offset de chaque processus et le mets dans un tableau y
                MPI_Gather(&last_local_offset, 1, MPI_LONG_LONG_INT, y, 1, MPI_LONG_LONG_INT, 0, split_comm);

                if (split_rank == 0) {
                    for (k2 = 1; k2 < (split_size + 1); k2++) {
                        y2[k2] = y[k2 - 1];
                    }
                }

                if (split_rank == 0) {
                    for (k2 = 1; k2 < (split_size + 1); k2++) {
                        y2[k2] = y2[k2 - 1] + y2[k2];
                    }
                }

                size_t offset_to_add = 0;
                //The n0 rank shares y2 (which is the offset of the previous rank) between all ranks into offset_to_add
                MPI_Scatter(y2, 1, MPI_LONG_LONG_INT, &offset_to_add, 1, MPI_LONG_LONG_INT, 0, split_comm);

                free(y);
                free(y2);

                //we add offset of the previous rank
                for (k1 = 0; k1 < max_num_read; k1++) {
                    if (local_reads_sizes_sorted[k1] != 0) {
                        local_offset_dest_sorted[k1] += offset_to_add;

                    } else {
                        local_offset_dest_sorted[k1] = 0;
                    }
                }

                /*
                 * we update destination rank according to
                 * original number of reads read.
                 *
                 */

                //we compute the new rank dest according to max_num_read
                size_t previous_num_reads_per_job[dimensions];
                //we create a vector of size split_size with previous reads per job
                //Every rank gets the num read to treat of each rank, the information will be put in the new vector
                MPI_Allgather(&first_local_readNum, 1, MPI_LONG_LONG_INT, previous_num_reads_per_job, 1, MPI_LONG_LONG_INT, split_comm);

                /*
                    Récupérer avec un MPI_AllGather, la derniere coordonnée  de chaque rang et le premier read de chaque rang
                    De plus récupérer avec un MPI_AllGather, la derniere coordonnée mate de chaque rang et le premier read mate de chaque rang
                    Trouver un moyen pour récupérer les mates, mais est-il vraiment nécessaire ??
                */

                //Nombre de read pour chaque processus
                //size_t *nb_read_num = (size_t *)malloc(num_proc * sizeof(size_t));
                //MPI_Allgather(&first_local_readNum, 1, MPI_LONG_LONG_INT, nb_read_num, 1, MPI_LONG_LONG_INT, split_comm);

                // we compute the position of of the read in the first
                // reference without the zero padding of bitonic
                size_t pos_ref0 = 0;

                //we need the number of zeros we add for the padding
                size_t N0 = max_num_read * dimensions - total_num_read;

                // we compute the new rank for
                // the reads sorted by offset destination
                int h2 = 0;

                pos_ref0 = max_num_read * split_rank - N0;
                //fprintf(stderr, "[MPISORT] chosen_split_rank = %d \n", chosen_split_rank);
                //fprintf(stderr, "[MPISORT] chosen_rank = %d \n", chosen_rank);

                for (j = 0; j < max_num_read; j++) {
                    if ( local_reads_sizes_sorted[j] != 0) {
                        int new_rank = chosen_split_rank;
                        pos_ref0 = (max_num_read * split_rank + j) - N0;
                        size_t tmp2 = 0;
                        for (h2 = 0; h2 < dimensions; h2++) {
                            tmp2 += previous_num_reads_per_job[h2];//computes the true nb of reads for the current chromosome?
                            if ( pos_ref0 < tmp2)  {
                                new_rank = h2;
                                break;
                            }
                        }
                        local_dest_rank_sorted[j] = new_rank;
                    }
                }
                MPI_Barrier(split_comm);

		/*
			INTER RANK OVERLAP ALGORITHM
		*/

                /*
                    We make sure the reads for a given coordinate is not cut between ranks
                    To do that we run through the destination rank vector and
                    each time we jump a rank we verify the coordinates are different

                    if the coordinates are the same we keep the previous rank 
                    and we update the number of reads per job
		*/


		
                time_count = MPI_Wtime();

                size_t previous_rank_last_coord = 0;
                int previous_rank_last_dest = 0;
                size_t new_dest_rank = 0;
                size_t previous_coordinates = 0;
                size_t current_coordinates = 0;
                int previous_dest_rank = 0;
                int current_dest_rank = 0;
                
		        size_t j  = 0;
                int j4    = 0;
        		/*
        			INTRA RANK OVERLAP ALGORITHM
        		*/

                previous_coordinates  = local_reads_coordinates_sorted[j];
                previous_dest_rank    = local_dest_rank_sorted[j];
		                      
                
		        for (j = 1; j < max_num_read; j++) {
                
                    current_coordinates  = local_reads_coordinates_sorted[j];
                    current_dest_rank    = local_dest_rank_sorted[j];                   

                    if ((previous_coordinates == current_coordinates)){

                        if ( ( previous_coordinates != 0 ) && (  previous_dest_rank != current_dest_rank ) ){                   
                            new_dest_rank++;
                            local_dest_rank_sorted[j] = previous_dest_rank;
                        }
                    }		
                    else{
                        // we update previous_coordinates  
                        previous_dest_rank    = current_dest_rank;
                    }

                    previous_coordinates = current_coordinates;				
                }
                               

                /*
			    INTER RANK OVERLAP ALGORITHM
		        We make sure the reads for a given coordinate is not cut between ranks
                To do that we run through the destination rank vector and
                each time we jump a rank we verify the coordinates are different

                if the coordinates are the same we keep the previous rank 
                and we update the number of reads per job
		        */
               
                //rank i envoie au rank i + 1

                //Cwe care for overlapped coordinates between rank

                if(split_rank != (split_size - 1)){

                        MPI_Send(&local_reads_coordinates_sorted[max_num_read - 1], 1, MPI_LONG_LONG_INT, split_rank+1,  0, split_comm);
                        MPI_Send(&local_dest_rank_sorted[max_num_read - 1], 1, MPI_INT, split_rank+1,  1, split_comm);
                }
                if(split_rank!=0){
                        MPI_Recv(&previous_rank_last_coord, 1, MPI_LONG_LONG_INT, split_rank-1, 0, split_comm, MPI_STATUS_IGNORE);
                        MPI_Recv(&previous_rank_last_dest, 1, MPI_INT, split_rank-1, 1, split_comm, MPI_STATUS_IGNORE);
                }

		        j=0;

                //while(local_reads_coordinates_sorted[j]==0 && j < max_num_read ) {j++;j2++;j3++;}             
		
		
                while((j < max_num_read) && (previous_rank_last_coord==local_reads_coordinates_sorted[j]) 
                            && (local_dest_rank_sorted[j] != previous_rank_last_dest))
                {
		    
		          local_dest_rank_sorted[j]=previous_rank_last_dest;
		          new_dest_rank++;
		          j++;
                }
		

		        MPI_Barrier(split_comm);

		
                md_log_rank_debug( chosen_split_rank, "[MPISORT][COMPUTE NEW DEST RANK STEP 3] time spent = %f s\n", split_rank, MPI_Wtime() - time_count );
                
                
                // now we shall update the number of read per rank 
                // we call new_num_reads_per_job the vector holding 
                // the number of reads per rank after the Bruck dispatch 
                size_t *new_num_reads_per_job = malloc(dimensions * sizeof(size_t));
                
                // the vector num_reads_per_dest_rank holds the local number of reads
                // according to the destination rank 
                size_t *num_reads_per_dest_rank = malloc(dimensions * sizeof(size_t));
                for (j4 = 0; j4 < dimensions; j4++) num_reads_per_dest_rank[j4] = 0;

                // we fill the num_reads_per_dest_rank
                for (j = 0; j < max_num_read; j++){
                    if ( local_reads_sizes_sorted[j] != 0){
                        num_reads_per_dest_rank[local_dest_rank_sorted[j]]++;
                    }	
                }                
                
                //we check if changing rank of destination doesn't affect the total read to sort 

                size_t old_total_reads = 0;
                size_t new_total_reads = 0;                  

                
                // now we make a reduce operation on num_reads_per_dest_rank
                MPI_Allreduce( num_reads_per_dest_rank, new_num_reads_per_job, dimensions, MPI_LONG_LONG_INT, MPI_SUM, split_comm );
                        
                free(num_reads_per_dest_rank);

                size_t final_local_readNum = 0;	
                final_local_readNum = new_num_reads_per_job[split_rank];
                
                free(new_num_reads_per_job);

                MPI_Allreduce( &final_local_readNum, &new_total_reads , 1, MPI_LONG_LONG_INT, MPI_SUM, split_comm );
                MPI_Allreduce( &first_local_readNum, &old_total_reads , 1, MPI_LONG_LONG_INT, MPI_SUM, split_comm );
                assert (old_total_reads == new_total_reads);
                              

                /*
                 * REMOVE ZERO PADDING BEFORE BRUCK
                 *
                 */

                size_t offset  = 0;
                size_t numItems = 0;
                size_t num_read_for_bruck = 0;
                int *p = local_reads_sizes_sorted;

                if (p[0] != 0) {
                    offset = 0;
                };

                if (p[max_num_read - 1] == 0) {
                    offset = max_num_read;

                } else {
                    while ((*p == 0) && (offset < max_num_read )) {
                        offset++;
                        p++;
                    }
                }

                time_count = MPI_Wtime();

                if (offset > 0) {

                    // we remove zeros in the vector we have 2 cases
                    // the first offset <  max_num_read
                    // and the entire vector is null
                    if ( offset < max_num_read ) {

                        numItems = max_num_read - offset;

                        local_reads_coordinates_sorted_trimmed_for_bruck    = malloc(numItems * sizeof(size_t));
                        local_offset_source_sorted_trimmed_for_bruck        = malloc(numItems * sizeof(size_t));
                        local_offset_dest_sorted_trimmed_for_bruck          = malloc(numItems * sizeof(size_t));
                        local_reads_sizes_sorted_trimmed_for_bruck          = malloc(numItems * sizeof(int));
                        local_dest_rank_sorted_trimmed_for_bruck            = malloc(numItems * sizeof(int));
                        local_source_rank_sorted_trimmed_for_bruck          = malloc(numItems * sizeof(int));
                        size_t y = 0;

                        for (y = 0; y < numItems; y++) {

                            local_reads_coordinates_sorted_trimmed_for_bruck[y]    = local_reads_coordinates_sorted[y + offset];
                            local_offset_source_sorted_trimmed_for_bruck[y]        = local_offset_source_sorted[y + offset];
                            local_offset_dest_sorted_trimmed_for_bruck[y]          = local_offset_dest_sorted[y + offset];
                            local_reads_sizes_sorted_trimmed_for_bruck[y]          = local_reads_sizes_sorted[y + offset];
                            local_dest_rank_sorted_trimmed_for_bruck[y]            = local_dest_rank_sorted[y + offset];
                            local_source_rank_sorted_trimmed_for_bruck[y]          = local_source_rank_sorted[y + offset];
                        }

                        num_read_for_bruck = numItems;

                        /*
                         *
                         * FOR DEBUG
                         *

                        for(y = 0; y < num_read_for_bruck; y++){
                            assert( local_reads_sizes_sorted_trimmed_for_bruck[y]       != 0 );
                            assert( local_source_rank_sorted_trimmed_for_bruck[y]       < dimensions);
                            assert( local_dest_rank_sorted_trimmed_for_bruck[y]         < dimensions);
                            assert( local_offset_source_sorted_trimmed_for_bruck[y]     != 0);
                            assert( local_offset_dest_sorted_trimmed_for_bruck[y]       != 0);
                            assert( local_reads_coordinates_sorted_trimmed_for_bruck[y] != 0);
                        }
                        */

                    } else {

                        numItems = 0;
                        local_reads_coordinates_sorted_trimmed_for_bruck    = malloc(numItems * sizeof(size_t));
                        local_offset_source_sorted_trimmed_for_bruck        = malloc(numItems * sizeof(size_t));
                        local_offset_dest_sorted_trimmed_for_bruck          = malloc(numItems * sizeof(size_t));
                        local_reads_sizes_sorted_trimmed_for_bruck          = malloc(numItems * sizeof(int));
                        local_dest_rank_sorted_trimmed_for_bruck            = malloc(numItems * sizeof(int));
                        local_source_rank_sorted_trimmed_for_bruck          = malloc(numItems * sizeof(int));
                        num_read_for_bruck = 0;
                    }

                } 
                else {

                    numItems = local_readNum;
                    local_reads_coordinates_sorted_trimmed_for_bruck    = malloc(local_readNum * sizeof(size_t));
                    local_offset_source_sorted_trimmed_for_bruck        = malloc(local_readNum * sizeof(size_t));
                    local_offset_dest_sorted_trimmed_for_bruck          = malloc(local_readNum * sizeof(size_t));
                    local_reads_sizes_sorted_trimmed_for_bruck          = malloc(local_readNum * sizeof(int));
                    local_dest_rank_sorted_trimmed_for_bruck            = malloc(local_readNum * sizeof(int));
                    local_source_rank_sorted_trimmed_for_bruck          = malloc(local_readNum * sizeof(int));

                    size_t y = 0;

                    for (y = 0; y < local_readNum; y++) {

                        local_reads_coordinates_sorted_trimmed_for_bruck[y]    = local_reads_coordinates_sorted[y];
                        local_offset_source_sorted_trimmed_for_bruck[y]        = local_offset_source_sorted[y];
                        local_offset_dest_sorted_trimmed_for_bruck[y]          = local_offset_dest_sorted[y];
                        local_reads_sizes_sorted_trimmed_for_bruck[y]          = local_reads_sizes_sorted[y];
                        local_dest_rank_sorted_trimmed_for_bruck[y]            = local_dest_rank_sorted[y];
                        local_source_rank_sorted_trimmed_for_bruck[y]          = local_source_rank_sorted[y];
                    }

                    num_read_for_bruck = numItems;

                    /*
                     *
                     * FOR DEBUG
                     *
                    for(y = 0; y < num_read_for_bruck; y++){
                        assert( local_reads_sizes_sorted_trimmed_for_bruck[y]       != 0 );
                        assert( local_source_rank_sorted_trimmed_for_bruck[y]       < dimensions);
                        assert( local_dest_rank_sorted_trimmed_for_bruck[y]         < dimensions);
                        assert( local_offset_source_sorted_trimmed_for_bruck[y]     != 0);
                        assert( local_offset_dest_sorted_trimmed_for_bruck[y]       != 0);
                        assert( local_reads_coordinates_sorted_trimmed_for_bruck[y] != 0);
                    }
                    */
                }

                free(local_reads_coordinates_sorted);
                free(local_offset_source_sorted);
                free(local_offset_dest_sorted);
                free(local_reads_sizes_sorted);
                free(local_dest_rank_sorted);
                free(local_source_rank_sorted);


                //if (split_rank == chosen_split_rank) {
                //    fprintf(stderr,   "rank %d :::::[MPISORT][TRIMMING] time spent = %f s\n", split_rank, MPI_Wtime() - time_count);
                //}
                md_log_rank_debug(chosen_split_rank, "[MPISORT][TRIMMING] time spent = %f s\n", split_rank, MPI_Wtime() - time_count);

                /*
                 * We do a Bruck on rank of origin reading
                 */

                size_t m = 0;
                int m_int = 0;
                int num_proc = dimensions;
                size_t *number_of_reads_by_procs = calloc( dimensions, sizeof(size_t));

                //fprintf(stderr,   "rank %d :::::[MPISORT] num_read_for_bruck = %zu \n", split_rank, num_read_for_bruck);

                for (m = 0; m < num_read_for_bruck; m++) {
                    //assert(new_pbs_orig_rank_off_phase1[m] < dimensions);
                    //assert(new_pbs_dest_rank_phase1[m] < dimensions);
                    number_of_reads_by_procs[local_source_rank_sorted_trimmed_for_bruck[m]]++;
                }

                int *local_source_rank_sorted_trimmed_for_bruckv2 = malloc( num_read_for_bruck * sizeof(int));

                for (m = 0; m < num_read_for_bruck; m++) {
                    local_source_rank_sorted_trimmed_for_bruckv2[m] = local_source_rank_sorted_trimmed_for_bruck[m];
                }

                size_t count6 = 0;
                for (m_int = 0; m_int < dimensions; m_int++) count6 += number_of_reads_by_procs[m_int];
                
                assert( count6 == num_read_for_bruck );
                MPI_Barrier(split_comm);

                size_t **reads_coordinates      = malloc(sizeof(size_t *) * dimensions);
                size_t **local_source_offsets   = malloc(sizeof(size_t *) * dimensions);
                size_t **dest_offsets           = malloc(sizeof(size_t *) * dimensions);
                int **read_size                 = malloc(sizeof(int *) * dimensions);
                int **dest_rank                 = malloc(sizeof(int *) * dimensions);
                int **source_rank               = malloc(sizeof(int *) * dimensions);

                /*
                 * We send in order
                 *
                 * local_offset_source_sorted_trimmed_for_bruck
                 * local_dest_rank_sorted_trimmed_for_bruck
                 * local_reads_coordinates_sorted_trimmed_for_bruck
                 * local_reads_sizes_sorted_trimmed_for_bruck
                 *
                 */

                
                time_count = MPI_Wtime();

                bruckWrite3(split_comm, 
                            split_rank,
                            dimensions,
                            count6,
                            number_of_reads_by_procs,
                            local_source_rank_sorted_trimmed_for_bruckv2,
                            local_offset_source_sorted_trimmed_for_bruck,     //offset sources
                            &local_source_offsets,
                            local_dest_rank_sorted_trimmed_for_bruck,         //destination rank
                            &dest_rank,
                            local_reads_coordinates_sorted_trimmed_for_bruck, //reads coordinates
                            &reads_coordinates,
                            local_reads_sizes_sorted_trimmed_for_bruck,       //read size
                            &read_size,
                            local_source_rank_sorted_trimmed_for_bruck,       //source rank
                            &source_rank,
                            local_offset_dest_sorted_trimmed_for_bruck,
                            &dest_offsets
                           );

                //if (split_rank == chosen_split_rank) {
                //    fprintf(stderr, "rank %d :::::[MPISORT][BRUCK 3] time spent = %f s\n", split_rank, MPI_Wtime() - time_count);
                //}

                md_log_rank_debug(chosen_split_rank, "[MPISORT][BRUCK 3] time spent = %f s\n", split_rank, MPI_Wtime() - time_count);


                time_count = MPI_Wtime();

                free(local_reads_coordinates_sorted_trimmed_for_bruck);
                free(local_dest_rank_sorted_trimmed_for_bruck);
                free(local_reads_sizes_sorted_trimmed_for_bruck);
                free(local_offset_source_sorted_trimmed_for_bruck);
                free(local_offset_dest_sorted_trimmed_for_bruck);
                free(local_source_rank_sorted_trimmed_for_bruck);
                free(local_source_rank_sorted_trimmed_for_bruckv2);

                local_reads_coordinates_sorted_trimmed    = malloc(first_local_readNum * sizeof(size_t));
                local_offset_source_sorted_trimmed        = malloc(first_local_readNum * sizeof(size_t));
                local_offset_dest_sorted_trimmed          = malloc(first_local_readNum * sizeof(size_t));
                local_dest_rank_sorted_trimmed            = malloc(first_local_readNum * sizeof(int));
                local_source_rank_sorted_trimmed          = malloc(first_local_readNum * sizeof(int));
                local_reads_sizes_sorted_trimmed          = malloc(first_local_readNum * sizeof(int));

                //if (split_rank == chosen_split_rank) {
                //    fprintf(stderr, "rank %d :::::[MPISORT][FREE + MALLOC] time spent = %f s\n", split_rank, MPI_Wtime() - time_count);
                //}

                md_log_rank_debug(chosen_split_rank, "[MPISORT][FREE + MALLOC] time spent = %f s\n", split_rank, MPI_Wtime() - time_count);


                /*
                 * GET DATA AFTER BRUCK
                 *
                 */

                j = 0;
                size_t k = 0;

                for (m_int = 0; m_int < num_proc; m_int++) {

                    for (k = 0; k < number_of_reads_by_procs[m_int]; k++) {

                        local_offset_dest_sorted_trimmed[k + j]         = dest_offsets[m_int][k];
                        local_dest_rank_sorted_trimmed[k + j]           = dest_rank[m_int][k];
                        local_reads_sizes_sorted_trimmed[k + j]         = read_size[m_int][k];
                        local_offset_source_sorted_trimmed[k + j]       = local_source_offsets[m_int][k];
                        local_reads_coordinates_sorted_trimmed[k + j]   = reads_coordinates[m_int][k];
                        local_source_rank_sorted_trimmed[k + j]         = source_rank[m_int][k];

                    }

                    free(dest_offsets[m_int]);
                    free(dest_rank[m_int]);
                    free(read_size[m_int]);
                    free(local_source_offsets[m_int]);
                    free(reads_coordinates[m_int]);
                    free(source_rank[m_int]);
                    j += number_of_reads_by_procs[m_int];
                }

                free(number_of_reads_by_procs);

                if (dest_rank != NULL) {
                    free(dest_rank);
                }

                if (read_size != NULL) {
                    free(read_size);
                }

                if (local_source_offsets != NULL) {
                    free(local_source_offsets);
                }

                if (reads_coordinates != NULL) {
                    free(reads_coordinates);
                }

                if (source_rank != NULL) {
                    free(source_rank);
                }

                if (dest_offsets != NULL) {
                    free(dest_offsets);
                }

                local_readNum = first_local_readNum;


                /*
                 *
                 * FOR DEBUG
                 *
                for ( j = 0; j < local_readNum; j++){
                    assert ( local_reads_coordinates_sorted_trimmed[j]    != 0 );
                    assert ( local_offset_source_sorted_trimmed[j]        != 0 );
                    assert ( local_offset_dest_sorted_trimmed[j]          != 0 );
                    assert ( local_reads_sizes_sorted_trimmed             != 0 );
                    assert ( local_dest_rank_sorted_trimmed[j]            < split_size );
                    assert ( local_source_rank_sorted_trimmed[j]          < split_size );
                }
                */

                free(local_reads_coordinates_sorted_trimmed);

                md_log_rank_debug(chosen_split_rank, "[MPISORT] we call write SAM\n");

                malloc_trim(0);

                time_count = MPI_Wtime();
                
                writeSam(
                    split_rank,
                    output_dir,
                    header,
                    local_readNum,
                    chrNames[i],
                    split_size,
                    split_comm,
                    chosen_split_rank,
                    finfo,
                    compression_level,
                    local_offset_dest_sorted_trimmed,
                    local_offset_source_sorted_trimmed,
                    local_reads_sizes_sorted_trimmed,
                    local_dest_rank_sorted_trimmed,
                    local_source_rank_sorted_trimmed,
                    local_data,
                    goff[rank],
                    first_local_readNum,
		            final_local_readNum
                );

                //if (split_rank == chosen_split_rank) {
                //    fprintf(stderr, "rank %d :::::[MPISORT][WRITESAM] chromosom %s :::  %f seconds\n\n\n", split_rank, chrNames[i], MPI_Wtime() - time_count);
                //}
                md_log_rank_info(chosen_split_rank, "[MPISORT][WRITESAM] chromosom %s :: %f seconds\n\n", chrNames[i], MPI_Wtime() - time_count);

            } else {

                /*
                 * We are in the case the number of cpu is
                 * not a power of 2
                 *
                 *
                 */

                parallel_sort_any_dim(
                    dimensions,                 //dimension for parabitonic
                    local_readNum,
                    split_rank,
                    split_size,
                    reads,
                    i,                          //chromosom number
                    chosen_split_rank,
                    split_comm,
                    localReadNumberByChr,
                    local_data,
                    output_dir,
                    finfo,
                    compression_level,
                    total_reads_by_chr,
                    goff[rank],
                    headerSize,
                    header,
                    chrNames[i]  );

            } //end if dimensions < split_rank

        } //if ((local_color == 0) && (i < (nbchr - 2))) //in the splitted dimension

        else {
            //we do nothing here
        }

        //we put a barrier before freeing pointers
        MPI_Barrier(MPI_COMM_WORLD);

        //we free the file pointer
        if  (file_pointer_to_free) {
            MPI_File_close(&mpi_file_split_comm);
        }

        //we free the split_comm
        if (split_comm_to_free) {
            MPI_Comm_free(&split_comm);
        }

        //if (color_vec_to_send) free(color_vec_to_send);
        //if (key_vec_to_send) free(key_vec_to_send);

    }// end loop upon chromosoms (line 665)

    MPI_Info_free(&finfo);
    free(goff);
    free(local_data);

    for ( i = 0; i < nbchr - 2; i++) {
        Read *tmp = reads[i], *next;

        while (tmp) {

            next = tmp->next;
            if (next == NULL) {free(tmp); break;} 
            if (tmp) free(tmp); 
            else break;
            tmp = next;
        }
    }

    MPI_Barrier(MPI_COMM_WORLD);

    free(header); //ok
    free(localReadNumberByChr); //ok


    for ( i = 0; i < nbchr; i++) {
        free(chrNames[i]);
    }

    free(reads); //ok
    free(readNumberByChr); //ok

    malloc_trim(0);
    res = MPI_Finalize();
    assert(res == MPI_SUCCESS);

    return 0;
}



static void usage(const char *prg) {

    fprintf(stderr, "Program: MPI version for sorting and marking duplicates of paired aligned FASTQ data\n"
            "Version: v%s\n"
            "Contact 1: Frederic Jarlier (frederic.jarlier@curie.fr) \n"
            "usage : mpirun -n TOTAL_PROC %s FILE_TO_SORT OUTPUT_FOLDER -q QUALITY -d OPTICAL_DISTANCE -v VERBOSE\n"
            "TOTAL_PROC should be a power of 2.\n"
            "OPTICAL_DISTANCE is 0 by default\n"
            "VERBOSE is the log file severity (4 by default)\n"
            "OUTPUT_FOLDER : a gz files per chromosome, a gz file of unmapped reads \n"
            "a gz files of discordants reads. \n"
            "Discordants reads are reads where one pairs align on a chromosome \n"
            "and the other pair align on another chromosome \n"
            "Unmapped reads are reads without coordinates on any genome \n"
            "Requirements : gcc > 4.8, automake 1.15, autoconf 2.69 and a MPI compiler\n"
            ""
            ""
            "For perfomances matters the file you want to sort could be \n"
            "stripped on a parallel file system. \n"
            "With Lustre you mention it with lfs setstripe command like this \n"
            "lfs set stripe -c stripe_number -s stripe_size folder           \n"
            , VERSION, prg);

    return;
}
