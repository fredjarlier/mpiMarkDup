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
    mark_duplicates.c

   Authors:
    Frederic Jarlier,                     Institut Curie
    Firmin Martin (main developer),       Paris Descartes
*/

/**
 * @todo Need to test discordant case on different data. It's possible
 *       that we can't extract discordant in the sort phase. Because
 *       a fragment can be duplicate of a discordant.
 * @todo Complete read and their mate in overlapped coordinates.
 *       (If we can avoid overlapped case, we can delete ensureMateRank()
 *        and we don't need to complete 'cluster'.)
 **/

#define _GNU_SOURCE

#include "mark_duplicates_discordant.h"
#include "mark_duplicates.h"
#include "createLBList.h"
#include <string.h>
#include <math.h>

#include "khash.h"
#include "write.h"
#define SPLIT_FACTOR 4

/**
 * construct hash table type and its operations for readInfo
 * key := read->Qname (char*)
 * value := read (readInfo*)
 * */

KHASH_MAP_INIT_STR(read, readInfo *)

extern int opticalDistance;

/**
 * @date 2018 Mar 23
 * @brief Convert a read to a coordinate mate info data structure
 * @param[in] read a read
 * @param[out] a coordinate mate info
 */

/**
 * @date 2018 Apr 15
 * @brief Build a paired-end and insert it in @p readEndsList list
 * @param[in] read1 first read
 * @param[in] read2 second read
 * @param[in] readEndsList paired-end list
 * @return read which is representing the pair
 */

readInfo *buildReadEndsDiscordant(readInfo *read1, readInfo *read2, llist_t *readEndsList, int case_insert) {
    //md_log_trace("read1=%s, read2=%s, read1->unclippedCoordPos=%zu, read1->coordMatePos=%zu, read2->unclippedCoordPos=%zu, read2->coordMatePos=%zu, read1->orientation=%d, read2->orientation=%d\n", read1->Qname, read2->Qname, read1->unclippedCoordPos, read1->coordMatePos, read2->unclippedCoordPos, read2->coordMatePos, read1->orientation, read2->orientation);

    
    //get unclipped mate coord
    read1->coordMatePos = read2->unclippedCoordPos;
    read2->coordMatePos = read1->unclippedCoordPos;


    read2->orientation = getOrientation(read2, 1);
    read1->orientation = getOrientation(read1, 1);

    // read and mate in same position, force to FR
    if (read1->unclippedCoordPos == read2->unclippedCoordPos) {
        read2->orientation = 3;//FR;
        read1->orientation = 3;//FR;

    }

    //if (read2->readChromosome > read1->readChromosome || 
    //    (read2->readChromosome == read1->readChromosome && read2->unclippedCoordPos >= read1->unclippedCoordPos)) {
    //if ( read2->unclippedCoordPos >= read1->unclippedCoordPos ){    
    if ( case_insert == 0 ){ 
        insertReadInList(readEndsList, read1);
        return read2;
        //md_log_trace(" normal case :: read2=%s, read1=%s, read2->orientation=%d, read1->orientation=%d\n", read2->Qname, read1->Qname, read2->orientation, read1->orientation);

    } else {

        insertReadInList(readEndsList, read2);
        //orientation first = getOrientation(read2, 0);
        //orientation second = getOrientation(read1, 0);

        unsigned int first = getOrientation(read2, 0);
        unsigned int second = getOrientation(read1, 0);

        if (first == 1 /*R*/) {
            if (second == 1 /*R*/) {
                read2->orientation = 4;//RR;
                read1->orientation = 4;//RR;

            } else {
                read2->orientation = 5;//RF;
                read1->orientation = 5;//RF;
            }

        } else {
            if (second == 1 /*R*/) {
                read2->orientation = 3;//FR;
                read1->orientation = 3;//FR;

            } else {
                read2->orientation = 2;//FF;
                read1->orientation = 2;//FF;
            }
        }

        //md_log_trace(" reverse case :: read2=%s, read1=%s, read2->orientation=%d, read1->orientation=%d\n", read2->Qname, read1->Qname, read2->orientation, read1->orientation);
        return read1;
    }
}


/**
 *   @date 2018 Feb 26
 *   @brief Fill reads, libraries list, reads array and reads line
 *   @param[in] bufferReads read buffer
 *   @param[in] intervalByProc interval coordinate by process
 *   @param[out] fragList fragments list
 *   @param[out] readEndsList paired-end list
 *   @param[out] readArr array of reads
 *   @param[out] samTokenLines array of reads line in sam file
 *   @param[in] readNum total amount of read lines
 *   @param[in] readIndex first read's index in local buffer
 *   @param[in] chrList array of chromosome name
 *   @param[in] chrNum chromosome number
 *   @param[in] comm markDuplicate communicator
 *   @return amount of read handled (i.e. read which will be insert in libraries list and hash table)
 *   @todo consider to use a struct to pack readInfo array and string array
 *   @todo consider to use a struct to pack readNum and readIndex
 */

size_t parseLibrariesDiscordant(char *bufferReads, 
                                Interval *intervalByProc, 
                                llist_t *fragList, 
                                llist_t *readEndsList, 
                                readInfo ***readArr, 
                                char ***samTokenLines, 
                                size_t readNum, 
                                size_t readIndex, 
                                chrInfo *chr, 
                                lbInfo *lb,
                                size_t *disc_offsets_source, 
                                MPI_Comm comm) {
    int rank, num_proc;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &num_proc);

    char *q = bufferReads;
    int percentage = 1, increment = 10;
    size_t readCounter = readIndex;
    float progression = 0;
    readInfo *read = NULL;

    /**
     * Allocate two array :
     *  1. readArr :: readInfo* array
     *  2. samTokenLines :: corresponding string array
     *
     *  If read is unmapped the read in readArr will be NULL.
     * */

    *readArr = malloc(readNum * sizeof(readInfo *));
    *samTokenLines = malloc(readNum * sizeof(char *));

    size_t splitLineSize = (readNum / SPLIT_FACTOR) + 1, splitNumber = 1;
    khash_t(read)* hash = kh_init(read);
    khint_t k;
    int ret;

    size_t readEndscall = 0;
    size_t counter = 0;
    while (*q) {
        progression = 100 * ((float)readCounter / readNum);
        char *tok;

        // parse read
        getLine(&q, &tok);
        read = readParsingDiscordant(tok, intervalByProc, readCounter, counter, chr, lb, disc_offsets_source, comm);
        
        counter++;

        if (read) {
            chr->lastChr = read->readChromosome;
        }

        // insert and cluster read
        unsigned int firstInPair, readPaired, mateUnmapped;

        if (read)  {
            assert(read->valueFlag);
            
            readPaired = readBits((unsigned int)read->valueFlag, 0);
            mateUnmapped = readBits((unsigned int)read->valueFlag, 3);
            firstInPair = readBits((unsigned int)read->valueFlag, 6);

            if (firstInPair)
                read->pair_num = 1;
            else
                read->pair_num = 2;

            insertReadInList(fragList, read);

            if (readPaired && !mateUnmapped) {
                k = kh_put_read(hash, read->Qname, &ret);


                if (ret != 0) {

                    /**
                     * Read never seen before, hash it in table
                     * */

                    kh_val(hash, k) = read;

                }  else {

                    /**
                     * Found first end before, we construct a pair.
                     * */

                    /* k is the key value */
                    khiter_t k = kh_get(read, hash, read->Qname);

                    if (k == kh_end(hash)) {
                        md_log_rank_error(rank, "khash : element is absent, should not happen\n");
                        assert(0);
                    }

                    /* Retrieve the read seen before */
                    readInfo *end = kh_val(hash, k);

                    /* Construct paired-end */
                    buildReadEndsDiscordant(end, read, readEndsList, 0);
                    readEndscall++;
                }
            }
        }

        (*samTokenLines)[readCounter - readIndex] = tok;

        // keep memory constant by spliting
        if (readCounter >= splitLineSize * splitNumber) {
            size_t remainsSize = strlen(q) + 1;
            memmove(bufferReads, q, remainsSize);
            bufferReads = realloc(bufferReads, remainsSize);
            q = bufferReads;
            splitNumber++;
        }

        (*readArr)[readCounter - readIndex] = read;
        readCounter++;

        if (progression > percentage) {
            md_log_debug("Parsed and sorted %.0f%% of reads\n", progression);
            percentage += increment;
        }

    }

    kh_destroy(read, hash);
    //fprintf (stderr, "readEndscall = %zu \n", readEndscall);

    //md_log_debug("Ensure reads mate rank which are susceptible in multiple processes ...\n");
    //ensureMateRank(*readArr, intervalByProc, readNum, comm);

    free(bufferReads);
    assert(readCounter - readIndex == readNum);
    return readCounter;
}


/**
 *  @date 2018 Feb 26
 *  @brief Parse a read line and fill a read structure. The resultant read is inserted in a library.
 *  @param[in] sam_buff a read line to parse
 *  @param[in] intervalByProc interval of coordinate by process
 *  @param[in] chrList array of chromosome name
 *  @param[in] chrNum chromosome number
 *  @param[in] lastChr last read chromosome
 *  @param[in] comm markDuplicate communicator
 *  @return the read filled.
 *  @note lastChr is used to optimize the parsing of field chromosome
 */

readInfo *readParsingDiscordant (char *sam_buff, 
                                Interval *intervalByProc, 
                                size_t readIndex,
                                size_t counter, 
                                chrInfo *chr, 
                                lbInfo *lb,
                                size_t *disc_offsets_source,  
                                MPI_Comm comm) {

    int rank, num_proc;
    char *q = sam_buff;
    char *tokenCar;
    char *u;
    unsigned int readUnmapped;
    unsigned int mateUnmapped;
    unsigned int secondaryAlignment;
    unsigned int supplementaryAlignment;
    unsigned int firstInPair;
    int i = 0;
    size_t score;
    readInfo *read = calloc(1, sizeof(readInfo));
   
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &num_proc);

    read->indexAfterSort = readIndex;

    while (1) {
        switch (i++) {

            // These fields are ignored

            case 4: // The MAPQ
            case 8: // the TLEN
            case 9: // the segment SEQ
                q = strchr(q, '\t') + 1;
                break;

            case 0:
                getTokenTab(&q, &tokenCar);
                read->Qname = tokenCar;
                break;

            case 1: // Flag part
                getTokenTab(&q, &tokenCar);
                read->valueFlag  = atoi(tokenCar);
                
                assert(read->valueFlag);    
                readUnmapped = readBits((unsigned int)read->valueFlag, 2);
                mateUnmapped = readBits((unsigned int)read->valueFlag, 3);
                secondaryAlignment = readBits((unsigned int)read->valueFlag, 8);
                supplementaryAlignment = readBits((unsigned int)read->valueFlag, 11);
                firstInPair = readBits((unsigned int)read->valueFlag, 6);
                //fprintf (stderr, "[IN readParsing] ::: read->Qname = %s ::: firstInPair = %u \n", read->Qname, firstInPair);
                //fprintf (stderr, "[IN readParsing] ::: read->Qname = %s ::: flag = %u \n", read->Qname, read->valueFlag );
                if (firstInPair == 1)
                    read->pair_num = 1;
                else
                    read->pair_num = 2;
                
        /*
        fprintf (stderr, "[IN READ2MATEFP] ::: read->Qname = %s ::: tokenCar = %s \n", read->Qname, tokenCar);
        fprintf (stderr, "[IN READ2MATEFP] ::: read->Qname = %s ::: read->valueFlag = %zu \n", read->Qname, read->valueFlag);
        fprintf (stderr, "[IN READ2MATEFP] ::: read->Qname = %s ::: readUnmapped = %d \n", read->Qname, readUnmapped);
        fprintf (stderr, "[IN READ2MATEFP] ::: read->Qname = %s ::: mateUnmapped = %d \n", read->Qname, mateUnmapped);
        fprintf (stderr, "[IN READ2MATEFP] ::: read->Qname = %s ::: secondaryAlignment = %d \n", read->Qname, secondaryAlignment);
        fprintf (stderr, "[IN READ2MATEFP] ::: read->Qname  = %s ::: supplementaryAlignment = %d \n", read->Qname, supplementaryAlignment);
        */

                // https://gatkforums.broadinstitute.org/gatk/discussion/6747/how-to-mark-duplicates-with-markduplicates-or-markduplicateswithmatecigar
                // secondary and supplementary alignment records are skipped and never flagged as duplicate.
                if (supplementaryAlignment || secondaryAlignment) {
                    freeRead(read);
                    free(tokenCar);
                    return NULL;
                }

                //unmapped : we ignore, we don't mark it according to Mark Duplicate
                if (readUnmapped) {
                    freeRead(read);
                    free(tokenCar);
                    return NULL;
                }

                if (mateUnmapped) {
                    read->mateChromosome = -1;
                }

                read->orientation = getOrientation(read, 0);

                free(tokenCar);
                break;

            case 2: // RNAME part
                getTokenTab(&q, &tokenCar);

                /* read's chromosome is the same as its previous read */
                if (strcmp(chr->chrList[chr->lastChr], tokenCar) == 0) {
                    read->readChromosome = chr->lastChr;

                } else {
                    for (int j = 0; j < chr->chrNum; ++j) {
                        if (strcmp(chr->chrList[j], tokenCar) == 0) {
                            read->readChromosome = j;
                            break;
                        }
                    }
                }

                free(tokenCar);
                break;

            case 3: // The position chromosome wise
                getTokenTab(&q, &tokenCar);
                read->coordPos = atoll(tokenCar);
                free(tokenCar);
                break;

            case 5: // the CIGAR string
                getTokenTab(&q, &tokenCar);
                read->cigar = tokenCar;
                assert(read->cigar);

                if ( strcmp(read->cigar, "*") == 0 ) {
                    //shall not be marked as duplicate but it's mate yes
                    read->readChromosome = -1;

                } else {
                    fillUnclippedCoord(read);
                }

                /* TODO: we can free early cigar string after unclipped coordinate computation */

                break;

            case 6: // the RNEXT
                getTokenTab(&q, &tokenCar);

                /* mate chromosome same as read chromosome */
                if (strcmp(tokenCar, "=") == 0) {

                    //normally we should never get here
                    read->mateChromosome = read->readChromosome;
                    read->discordant = 0;

                } else {
                    for (int j = 0; j < chr->chrNum; ++j) {
                        if (strcmp(chr->chrList[j], tokenCar) == 0) {
                            read->mateChromosome = j;
                            read->discordant = 0;
                            read->offset_source_sam = disc_offsets_source[counter];
                            //fprintf(stderr, "[MARKDUPDISCORDANT] [READPARSING] Found a discordant %s at offset %zu \n", read->Qname, disc_offsets_source[counter]); 
                            break;
                        }
                    }
                }

                free(tokenCar);
                break;


            case 7: // here we have the PNEXT, the position of the next read
                getTokenTab(&q, &tokenCar);
                read->coordMatePos = atoll(tokenCar);

                // we tell if the mate is in the buffer
                // we check if the mate's coordinates are outer the buffer limits
                         
                checkMateBuffer(read, intervalByProc, comm);
                free(tokenCar);

                break;

            case 10:
                getTokenTab(&q, &tokenCar);
                // here we got the quality string
                u = tokenCar;
                score = 0;

                //we loop the token char and compute quality
                while (*u) {
                    score = fastqToPhred(*u);

                    /**
                     * Only take account base which has quality over 15
                     * See htsjdk/samtools/DuplicateScoringStrategy.java, function getSumOfBaseQualities
                     */

                    if (score >= 15) {
                        read->phred_score += (size_t)score;
                    }

                    u++;
                }

                /**
                 * MarkDuplicates use SUM_OF_BASE_QUALITIES as scoring strategy.
                 * See htsjdk/samtools/DuplicateScoringStrategy.java, function computeDuplicateScore()
                 * Note that 0x7FFF = 32767 is the max value of short type in java,
                 * Since phred_score is stored in short type and we need to compute pair phred score later,
                 * htsjdk takes the minimum value between the computed phred score and the half of short max value to
                 * prevent overflow.
                 */

                read->phred_score = min(read->phred_score, (int) (0x7FFF / 2));
                read->pairPhredScore = read->phred_score;
                free(tokenCar);
                break;

            default:
                break;
        }


        if (i == 11) {
            break;
        }
    }


    if (opticalDistance > 0) {
        read->physicalLocation = computePhysicalLocation(read);
    }
    assert(read->valueFlag);    
    assert (( read->pair_num == 1) || ( read->pair_num == 2 ));

    read->fingerprint = read2Fingerprint(read);
    read->mate_fingerprint = read2mateFP(read);

    assert(read->fingerprint);

    //Now we search in the flag for LB name
    while (getTokenTab(&q, &tokenCar)) {
        if (strncmp(tokenCar, "LB:Z:", strlen("LB:Z:")) == 0) {
            //fillReadLBValue(tokenCar, read);
            char *lbName = getReadTagValue(tokenCar, "LB:Z:");

            if (strcmp(lb->lbList[lb->lastLb], lbName) == 0) {
                read->readLb = lb->lastLb;

            } else {
                for (i = 0; i < lb->lbNum; ++i) {
                    if (strcmp(lb->lbList[i], lbName) == 0) {
                        read->readLb = i;
                    }
                }
            }

            free(lbName);
            free(tokenCar);
            break;
        }

        free(tokenCar);
    }
    
    return read;
}




/**
 *   @date 2018 Feb 26
 *   @brief Exchange extern fragments and complete pairs information (phredScore and mate index)
 *   @param[out] fragList fragments list
 *   @param[out] readEndsList paired-end list
 *   @param[out] htbl perfect hash table
 *   @param[in] comm markDuplicate communicator
 */

void exchangeExternFragDiscordant(llist_t *fragList, llist_t *readEndsList, hashTable *htbl, MPI_Comm comm) {

    int rank;
    MPI_Comm_rank(comm, &rank);

    size_t numberOfExternalMate = 0;
    countExternalMateInList(fragList, &numberOfExternalMate);

    mateInfo *mates = malloc(numberOfExternalMate * sizeof(mateInfo));
    size_t actualNumberToSend = getMateRankReadSizeBeforeBruck(fragList, &mates);

    assert(actualNumberToSend == numberOfExternalMate);

    /* Sort mates to send  by mate rank */
    qsort(mates, numberOfExternalMate, sizeof(mateInfo), cmpMateRank);


    /* Exchange mates and fill them to a readInfo array */
    readInfo **matesByProc;
    //int totalrecv = exchangeAndFillMate(&matesByProc, mates, numberOfExternalMate, comm);
    size_t totalrecv = exchangeAndFillMate_with_Bruck_v2(&matesByProc, mates, numberOfExternalMate, comm);

    md_log_rank_debug(rank, "[mpiMD][exchangeExternFragDiscordant] Received %zu mates, fragList size = %d\n", totalrecv, fragList->size);

    free(mates);

    //test if we have nothing to do we return
    if (totalrecv == 0) return; 

    int mateCounter = 0;
    for (size_t i = 0; i < totalrecv; i++) {
        readInfo *mate = getReadFromFingerprint(htbl, matesByProc[i]->fingerprint);
       
        /* Only external mate of the current buffer has a slot in hash table.
         * Note that the array matesByProc contains all externals mates among all process.
         * */

	if (mate == NULL) continue;
	//assert(mate != NULL);
	
        if (mate) {
            //assert(mate->external);
	    if (!mate->external) continue;

            /*
                We fill up missing information of matesByProc[i]
                LB, paire_num, orientation, coordMatePos

           */
           matesByProc[i]->readLb = mate->readLb;
           //matesByProc[i]->orientation = readBits((unsigned int)mate->valueFlag, 5);

           if (mate->pair_num == 1) matesByProc[i]->pair_num = 2; 
           if (mate->pair_num == 2) matesByProc[i]->pair_num = 1;

           matesByProc[i]->coordMatePos = mate->coordPos;

            /* free fictitious mate */
            freeRead(mate);
            /* insert external mate, it is partially filled as a readInfo */
            hashTableInsert(htbl, matesByProc[i]);
            
            for (lnode_t *node = fragList->head; node != fragList->nil; node = node->next) {
               
               	// compare clipped position first because read2mateFP is expensive
                if (node->read->coordPos == matesByProc[i]->coordMatePos) {
				
                   //unsigned long long fragMateFingerprint = read2mateFP(node->read);
		              unsigned long long fragMateFingerprint = node->read->mate_fingerprint;

		              if (matesByProc[i]->fingerprint == fragMateFingerprint) {

                        buildReadEndsDiscordant(matesByProc[i], node->read, readEndsList, 1);
                        insertReadInList(fragList, matesByProc[i]) ;
                        matesByProc[i]->Qname = strdup(node->read->Qname);
                        mateCounter++;
                        break;

                    }

                }

            }

        }
    }

    md_log_rank_trace(rank, "Found %d/%d mates\n", mateCounter, totalrecv);

    for (lnode_t *node = readEndsList->head; node != readEndsList->nil; node = node->next) {
        assert(node->read);
        readInfo *mate = getMateFromRead(htbl, node->read);

        if (!mate) {
            md_log_rank_error(rank, "Cannot found %s/2 of fingerprint %zu Please check if it's in the file or not.\n", node->read->Qname, node->read-> fingerprint, rank);
            exit(EXIT_FAILURE);
        }

        node->read->pairPhredScore += mate->phred_score;
        node->read->mateIndexAfterSort = mate->indexAfterSort;
       
    }
    /* free external mates in current library, others reads are free in destroyLBList */
    free(matesByProc);
}

/**
 * @date 2018 Apr 12
 * @brief Find and mark duplicate read in a library
 * @param[in, out] fragList fragments list
 * @param[in, out] readEndsList paired-ends list
 * @param[in, out] htbl a perfect hash table
 * @param[out] totalDuplica amount of duplicate
 * @param[out] totalOpticalDuplicate amount of optical duplicate
 * @param[in] comm mark_duplicates communicator
 * @note MarkDuplicates use term 'chunk' which corresponding our 'cluster'
 * @note See middle part of function generateDuplicateIndexes() in Markduplicates.java
 */

void findDuplicaDiscordant(llist_t *fragList, llist_t *readEndsList, hashTable *htbl, int *totalDuplica, int *totalOpticalDuplicate, MPI_Comm comm) {

    int rank, num_proc;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &num_proc);
    
    llist_t *nextCluster = llist_create();
    lnode_t *firstOfNextCluster = NULL;

    for (lnode_t *node = readEndsList->head; node != readEndsList->nil; node = node->next) {

        if (firstOfNextCluster && areComparableForDuplicates(firstOfNextCluster->read, node->read, 1)) {
            llist_append(nextCluster, node->read);

        } else {

            if (nextCluster->size > 1) {                
                markDuplicatePairs(nextCluster, htbl, totalDuplica, totalOpticalDuplicate);
            }
            llist_clear(nextCluster);
            llist_append(nextCluster, node->read);
            firstOfNextCluster = node;
        }
    }
    if (nextCluster->size > 1) {
        markDuplicatePairs(nextCluster, htbl, totalDuplica, totalOpticalDuplicate);
    }
    llist_clear(nextCluster);
    int containsPairs = 0;
    int containsFrags = 0;
    firstOfNextCluster = NULL;
       
    for (lnode_t *node = fragList->head; node != fragList->nil; node = node->next) {

        if (firstOfNextCluster && areComparableForDuplicates(firstOfNextCluster->read, node->read, 0)) {
            llist_append(nextCluster, node->read);
            unsigned int readIsPaired = isPaired(node->read);
            containsPairs = containsPairs || readIsPaired ;
            containsFrags = containsFrags || !readIsPaired;

        } else {

            if (nextCluster->size > 1 && containsFrags) {
                markDuplicateFragments(nextCluster, totalDuplica, containsPairs);
            }
            
            llist_clear(nextCluster);
            llist_append(nextCluster, node->read);
            firstOfNextCluster = node;
            unsigned int readIsPaired = isPaired(node->read);
            
            containsPairs = readIsPaired ;
            containsFrags = !readIsPaired;

        }

    }
   
    //llist_readInfo_print(nextCluster);
    //printf("containsFrags : %d, containsPairs : %d \n ",  containsFrags, containsPairs);
    if ( nextCluster->size > 1 )
        markDuplicateFragments(nextCluster, totalDuplica, containsPairs);
    //fprintf(stderr, "rank %d after markDuplicateFragments 2 totalDuplica = %d \n", rank, *totalDuplica);
    //llist_readInfo_print(nextCluster);
    llist_clear(nextCluster);
    llist_destruct(nextCluster);

}


/**
 * @date 2018 Feb 26
 * @brief Find and mark duplicate reads.
 * @param bufferReads reads buffers (i.e. SAM file without header)
 * @param readNum total amount of reads
 * @param header sam header
 * @param[in] comm markDuplicate communicator
 * @return a read buffer with duplicate reads marked.
 */

char *markDuplicateDiscordant (char *bufferReads, 
                                size_t readNum, 
                                char *header, 
                                MPI_Comm comm,
                                size_t *disc_offsets_source,
                                size_t **disc_dup_offset_source,
                                size_t *disc_dup_number) {

    //MPI_Comm previousComm = md_get_log_comm();
    //md_set_log_comm(comm);

    int rank, num_proc;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &num_proc);
    double timeStamp, timeStart = MPI_Wtime();
  
    chrInfo chr;
    lbInfo lb;
    initChrInfo(&chr, header);
    initLbInfo(&lb, header);

    if (chr.chrNum == 0) {
        md_log_error("There isn't chromosome in the header, please check the input sam file.");
    }

    size_t *readNumByProc = calloc(num_proc,  sizeof(Interval)), readIndexOffset = 0;
    MPI_Allgather(&readNum, 1, MPI_SIZE_T, readNumByProc, 1, MPI_SIZE_T, comm);

    for (int i = 0; i < rank; ++i) {
        readIndexOffset += readNumByProc[i];
    }


    llist_t *fragList = llist_create(), *readEndsList = llist_create();

    Interval interval = getIntervalFromBuffer(bufferReads, readNum);
    Interval *intervalByProc = gatherIntervalByProc(interval, comm);
    md_log_all_debug("readNum=%zu, interval.firstCoord=%u, interval.lastCoord=%u\n", readNum, interval.firstCoord, interval.lastCoord);

    md_log_debug("Parse and cluster reads ...\n");
    timeStamp = MPI_Wtime();
    readInfo **readArr;
    char **samTokenLines;
    parseLibrariesDiscordant(bufferReads, 
                    intervalByProc, 
                    fragList, 
                    readEndsList, 
                    &readArr, 
                    &samTokenLines, 
                    readNum, 
                    readIndexOffset, 
                    &chr,  
                    &lb, 
                    disc_offsets_source,
                    comm);

    md_log_debug("End of parsing and clustering %f seconds\n", MPI_Wtime() - timeStamp);

    md_log_debug("Construct perfect hashing table ...\n");
    timeStamp = MPI_Wtime();
    hashTable *htbl = malloc(sizeof(hashTable));
    readInfo **readArrWithExternal;

    size_t totalReadWithFictitiousMate = fillReadAndFictitiousMate(readArr, &readArrWithExternal, readNum);
    //fprintf(stderr, " totalReadWithFictitiousMate = %zu \n", totalReadWithFictitiousMate);
   
    shareHpAndConstructHtbl(htbl, readArrWithExternal, totalReadWithFictitiousMate, comm);
    free(readArrWithExternal);
    
    md_log_debug("Perfect hash table is constructed %f seconds\n", MPI_Wtime() - timeStamp);


    /**
     * TODO: add a phase of communication here for completion of 'cluster' in overlapped range case.
     *       - Received read should be insert in fragList via insertReadInList()
     *       - We should verify if we can construct a pair with these reads, to do this
     *         we may consider to use khash in function markDuplicates and pass it to parseLibraries().
     *         Hence, We can use this khash here to know if read is seen before or not.
     *       - If we see read's mate before we can construct pairs and insert them in readEndsList
     *         via function buildReadEnds().
     */


    /* TRACE fragList */
    
    /* End TRACE */

    md_log_debug("Start to exchange mates ...\n");
    timeStamp = MPI_Wtime();
    exchangeExternFragDiscordant(fragList, readEndsList, htbl, comm) ;
    md_log_debug("Finished to exchange mate in %f seconds\n", MPI_Wtime() - timeStamp);


    /*
    size_t fragList_size = 0;
    size_t readEnds_size = fragList_size = 0;

    for (lnode_t *node = fragList->head; node != fragList->nil; node = node->next) {
        //readCounter = readIndex,
        //md_log_trace("lb=%zu, chr=%zu, unclippedCoordPos=%zu, orientation=%zu, mchr=%zu, mateUnclippedCoordPos=%zu, rindex=%zu, mindex=%zu\n", read->readLb, read->readChromosome, read->unclippedCoordPos, read->orientation, read->mateChromosome, read->coordMatePos, read->indexAfterSort, read->mateIndexAfterSort);
        fragList_size++;
    }
    
    for (lnode_t *node = readEndsList->head; node != readEndsList->nil; node = node->next) {
        //readInfo *read = node->read;
        //md_log_trace("lb=%zu, chr=%zu, unclippedCoordPos=%zu, orientation=%zu, mchr=%zu, mateUnclippedCoordPos=%zu, rindex=%zu, mindex=%zu\n", read->readLb, read->readChromosome, read->unclippedCoordPos, read->orientation, read->mateChromosome, read->coordMatePos, read->indexAfterSort, read->mateIndexAfterSort);
        readEnds_size++;
    }

    fprintf(stderr, " TRACE 3 readEndsList size = %zu :::: fragList size =%zu \n", readEnds_size, fragList_size);
    */

    md_log_debug("Start to sort fragments list and read ends list ...\n");
    timeStamp = MPI_Wtime();

    /* Sort lists using merge sort.
     * We use merge sort instead of insertion sort because
     * insertion sort is too slow in discordant case.
     * TODO: If this merge sort work fine, replace insertReadInList() by
     *       llist_append() in all occurrences.
     * */
    llist_merge_sort(fragList, 1);
    llist_merge_sort(readEndsList, 0);

    /*
    size_t readEnds_size = 0;
    size_t fragList_size = 0;
    
    for (lnode_t *node = fragList->head; node != fragList->nil; node = node->next) {
        
        //md_log_trace("lb=%zu, chr=%zu, unclippedCoordPos=%zu, orientation=%zu, mchr=%zu, mateUnclippedCoordPos=%zu, rindex=%zu, mindex=%zu\n", read->readLb, read->readChromosome, read->unclippedCoordPos, read->orientation, read->mateChromosome, read->coordMatePos, read->indexAfterSort, read->mateIndexAfterSort);
        fragList_size++;
    }

    
    
    for (lnode_t *node = readEndsList->head; node != readEndsList->nil; node = node->next) {
        //md_log_trace("lb=%zu, chr=%zu, unclippedCoordPos=%zu, orientation=%zu, mchr=%zu, mateUnclippedCoordPos=%zu, rindex=%zu, mindex=%zu\n", read->readLb, read->readChromosome, read->unclippedCoordPos, read->orientation, read->mateChromosome, read->coordMatePos, read->indexAfterSort, read->mateIndexAfterSort);
        readEnds_size++;
    }
    
   fprintf(stderr, " TRACE 4 readEndsList size = %zu :::: fragList size =%zu \n", readEnds_size, fragList_size);
   */ 

    md_log_debug("Finished to sort list in %f seconds\n", MPI_Wtime() - timeStamp);

    md_log_info("Start to find duplicate ...\n");
    timeStamp = MPI_Wtime();
    int localDuplicates = 0, localOpticalDuplicates = 0, totalDuplicates = 0, totalOpticalDuplicates = 0;
    
    findDuplicaDiscordant(fragList, readEndsList, htbl, &localDuplicates, &localOpticalDuplicates, comm) ;
    


    md_log_info("Finished to find duplicates in %f seconds\n", totalDuplicates, totalOpticalDuplicates,  MPI_Wtime() - timeStamp);
    md_log_info("Gather duplicates results ...\n");
    
    timeStamp = MPI_Wtime();
    
    
    MPI_Reduce(&localDuplicates, &totalDuplicates, 1, MPI_INT, MPI_SUM, 0, comm);
    MPI_Reduce(&localOpticalDuplicates, &totalOpticalDuplicates, 1, MPI_INT, MPI_SUM, 0, comm);
    md_log_info("Found %d duplicates and %d optical duplicates in %f seconds\n\n", totalDuplicates, totalOpticalDuplicates,  MPI_Wtime() - timeStamp);

    timeStamp = MPI_Wtime();

    char *newBuff = writeBuff(samTokenLines, readArr, readNum);
    
    free(samTokenLines);
    free(intervalByProc);
    free(readNumByProc);
    
    
    //if (tmp_disc_dup_offset_source) free(tmp_disc_dup_offset_source);
    if (disc_offsets_source)        free(disc_offsets_source);

    freeChrInfo(&chr);
    freeLbInfo(&lb);

    hashTableDestroy(htbl);

    //md_set_log_comm(previousComm);

    MPI_Barrier(comm);
    md_log_trace("Finished to free data structures %f seconds\n", MPI_Wtime() - timeStamp);
    md_log_debug("Finished to mark buffer %f seconds\n", MPI_Wtime() - timeStamp);
    /*
        ///////////////////////////
        //////////////////////////
        Now we fill up the disc_offset source vector
        //////////////////////////
        //////////////////////////
        //////////////////////////
    */

    disc_dup_number[0] = localDuplicates + localOpticalDuplicates;

    //we return the total of duplicate for this rank
   
    // we loop the readArr
    // we use data_offset_source to check 
    // the discordant reads
    size_t j;
    size_t index_in_offset = 0;

     for ( j = 0; j < readNum; j++) {
        //we search discordant reads
        if (!readArr[j] || !readArr[j]->d) continue;
        else if (readArr[j]->d == 1){
            index_in_offset++;            
        }
    }

    size_t *disc_dup_offset_source_tmp = calloc( index_in_offset,  sizeof(size_t));
    
    index_in_offset = 0;

    for ( j = 0; j < readNum; j++) {
        //we search discordant reads
        if (!readArr[j] || !readArr[j]->d) continue;
        else if (readArr[j]->d == 1){
            disc_dup_offset_source_tmp[index_in_offset] =  readArr[j]->offset_source_sam;
            index_in_offset++;            
        }
    }
    *disc_dup_offset_source =  disc_dup_offset_source_tmp;
    *disc_dup_number = index_in_offset;

     /* free all reads and fragList */
    llist_destruct(fragList);
    /* free all nodes of readEndsList */
    llist_clear(readEndsList);
    /* free readEndsList */
    llist_destruct(readEndsList);

    //for (size_t m = 0; m < readNum; m++)
    //    free(readArr[m]);
    free(readArr);

    //*disc_dup_offset_source = NULL;
    md_log_debug("Total time spend in mark duplicates %f seconds\n", MPI_Wtime() - timeStart);
          
    MPI_Barrier(comm);

    return newBuff;
}

