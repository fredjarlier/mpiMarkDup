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

    assert(read1->pair_num == 1 || read1->pair_num == 2); 
    assert(read2->pair_num == 1 || read2->pair_num == 2);

    if (read1->pair_num == 1) assert(read2->pair_num == 2);
    if (read1->pair_num == 2) assert(read2->pair_num == 1);
  

    //get unclipped mate coord
    read1->coordMatePos = read2->unclippedCoordPos;
    read2->coordMatePos = read1->unclippedCoordPos;


    read2->orientation = getOrientation(read2, 1);
    read1->orientation = getOrientation(read1, 1);

    // read and mate in same position, force to FR
    if (read1->unclippedCoordPos == read2->unclippedCoordPos) {
        read2->orientation = FR;
        read1->orientation = FR;

    }

    //if (read2->readChromosome > read1->readChromosome || 
    //    (read2->readChromosome == read1->readChromosome && read2->unclippedCoordPos >= read1->unclippedCoordPos)) {
    //if ( read2->unclippedCoordPos >= read1->unclippedCoordPos ){    
    if ( case_insert == 0 ){ 
        insertReadInList(readEndsList, read1, 0);
        return read2;
        //md_log_trace(" normal case :: read2=%s, read1=%s, read2->orientation=%d, read1->orientation=%d\n", read2->Qname, read1->Qname, read2->orientation, read1->orientation);

    } else {

        insertReadInList(readEndsList, read2, 0);
        orientation first = getOrientation(read2, 0);
        orientation second = getOrientation(read1, 0);

        if (first == R) {
            if (second == R) {
                read2->orientation = RR;
                read1->orientation = RR;

            } else {
                read2->orientation = RF;
                read1->orientation = RF;
            }

        } else {
            if (second == R) {
                read2->orientation = FR;
                read1->orientation = FR;

            } else {
                read2->orientation = FF;
                read1->orientation = FF;
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

int parseLibrariesDiscordant(char *bufferReads, Interval *intervalByProc, llist_t *fragList, llist_t *readEndsList, readInfo ***readArr, char ***samTokenLines, size_t readNum, size_t readIndex, chrInfo *chr, lbInfo *lb, MPI_Comm comm) {
    int rank, num_proc;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &num_proc);

    char *q = bufferReads;
    int readCounter = readIndex, percentage = 1, increment = 10;
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

    while (*q) {
        progression = 100 * ((float)readCounter / readNum);
        char *tok;

        // parse read
        getLine(&q, bufferReads, &tok);
        read = readParsing(tok, intervalByProc, readCounter, chr, lb, comm);
        

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

            insertReadInList(fragList, read, 1);

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
 *   @date 2018 Feb 26
 *   @brief Exchange extern fragments and complete pairs information (phredScore and mate index)
 *   @param[out] fragList fragments list
 *   @param[out] readEndsList paired-end list
 *   @param[out] htbl perfect hash table
 *   @param[in] comm markDuplicate communicator
 */

void exchangeExternFragDiscordant(llist_t *fragList, llist_t *readEndsList, hashTable *htbl, Interval interval, MPI_Comm comm) {

    int rank, num_proc;
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
    int totalrecv = exchangeAndFillMate(&matesByProc, mates, numberOfExternalMate, comm);

    md_log_rank_debug(rank, "[mpiMD][exchangeExternFragDiscordant] Received %d mates, fragList size = %d\n", totalrecv, fragList->size);

    //test if we have nothing to do we return
    if (totalrecv == 0) return; 


    for (lnode_t *node = fragList->head; node != fragList->nil; node = node->next) {
                readInfo *read = node->read;
                //md_log_trace("lb=%zu, chr=%zu, unclippedCoordPos=%zu, orientation=%zu, mchr=%zu, mateUnclippedCoordPos=%zu, rindex=%zu, mindex=%zu\n", read->readLb, read->readChromosome, read->unclippedCoordPos, read->orientation, read->mateChromosome, read->coordMatePos, read->indexAfterSort, read->mateIndexAfterSort);
                assert( (read->pair_num == 1) || (read->pair_num == 2));
            }

            /* TRACE readEndsList */
            for (lnode_t *node = readEndsList->head; node != readEndsList->nil; node = node->next) {
                readInfo *read = node->read;
                //md_log_trace("lb=%zu, chr=%zu, unclippedCoordPos=%zu, orientation=%zu, mchr=%zu, mateUnclippedCoordPos=%zu, rindex=%zu, mindex=%zu\n", read->readLb, read->readChromosome, read->unclippedCoordPos, read->orientation, read->mateChromosome, read->coordMatePos, read->indexAfterSort, read->mateIndexAfterSort);
                assert( (read->pair_num == 1) || (read->pair_num == 2));
            }    



    int mateCounter = 0;
    for (int i = 0; i < totalrecv; i++) {
        readInfo *mate = getReadFromFingerprint(htbl, matesByProc[i]->fingerprint);
        assert( (mate->pair_num == 1) || (mate->pair_num == 2));
        /* Only external mate of the current buffer has a slot in hash table.
         * Note that the array matesByProc contains all externals mates among all process.
         * */

	if (mate == NULL) continue;
	//assert(mate != NULL);
	
        if (mate) {
            //assert(mate->external);
	    if (!mate->external) continue;
            /* free fictitious mate */
            freeRead(mate);
            /* insert external mate, it is partially filled as a readInfo */
            hashTableInsert(htbl, matesByProc[i]);

            /* insert mate in fragments list and readEnds list (if we can construct pair) 
             * TODO: To optimize :
             *  - mate distribution among process may be unbalanced 
             *    Illustration :
             *        rank 0 : send 2 mates
             *        rank 1 : send 3 mates
             *         ...        ...
             *        rank 28 : send 2923 mates
             *  - We need to go through fragments list because we can't deduce mate's fingerprint by read's fingerprint.
             *    We need at least mate's Qname and mate's Flag, but we don't exchange them.
             *  - In the example above, rank 28 need to do 2923 * fragList->size comparisons in worst case.
             *   
             * */
            
            for (lnode_t *node = fragList->head; node != fragList->nil; node = node->next) {
               
        		// compare clipped position first because read2mateFP is expensive
                if (node->read->coordPos == matesByProc[i]->coordMatePos) {
				

                   assert( (node->read->pair_num == 1) || (node->read->pair_num == 2));
                   
                   unsigned long long fragMateFingerprint = read2mateFP(node->read);
		    
		              if (matesByProc[i]->fingerprint == fragMateFingerprint) {

                        /*if (matesByProc[i]->valueFlag == 0)
                            fprintf(stderr, "on %d problem with %d read %s valueFlag is NULL \n", totalrecv, i, node->read->Qname);
			                   
                        

                        if (node->read->pair_num == 1){
                            matesByProc[i]->pair_num = 2;
                        }
                        if (node->read->pair_num == 2){
                            matesByProc[i]->pair_num == 1;
                        }*/
                        /*
                        fprintf(stderr, "Call buildsreadEnds read %s valueFlag is = %u pair = %u \n", node->read->Qname, node->read->valueFlag
                            , node->read->pair_num);


                        fprintf(stderr, "Call buildsreadEnds matesByProc[%d] %s valueFlag is = %u pair = %u \n", i, matesByProc[i]->Qname, matesByProc[i]->valueFlag
                            , matesByProc[i]->pair_num);
                        */

                        readInfo *insertedRead = buildReadEndsDiscordant(matesByProc[i], node->read, readEndsList, 1);
                        insertReadInList(fragList, matesByProc[i], 1) ;
                        matesByProc[i]->Qname = strdup(node->read->Qname);

                        assert ((node->read->pair_num == 1 ) || (node->read->pair_num ==2));
                        
                        if (node->read->pair_num == 1) assert (matesByProc[i]->pair_num == 2); 
                        if (node->read->pair_num == 2) assert (matesByProc[i]->pair_num == 1);
                        assert( (matesByProc[i]->pair_num == 1) || (matesByProc[i]->pair_num == 2));

                        //matesByProc[i]->valueFlag = readFlag2MateFlag(node->read->valueFlag);
                        mateCounter++;
                        break;

                    }

                }

            }

        }
    }

    md_log_rank_trace(rank, "Found %d/%d mates\n", mateCounter, totalrecv);

    for (lnode_t *node = fragList->head; node != fragList->nil; node = node->next) {
                readInfo *read = node->read;
                //md_log_trace("lb=%zu, chr=%zu, unclippedCoordPos=%zu, orientation=%zu, mchr=%zu, mateUnclippedCoordPos=%zu, rindex=%zu, mindex=%zu\n", read->readLb, read->readChromosome, read->unclippedCoordPos, read->orientation, read->mateChromosome, read->coordMatePos, read->indexAfterSort, read->mateIndexAfterSort);
                assert( (read->pair_num == 1) || (read->pair_num == 2));
            }

            /* TRACE readEndsList */
    for (lnode_t *node = readEndsList->head; node != readEndsList->nil; node = node->next) {
           readInfo *read = node->read;
           //md_log_trace("lb=%zu, chr=%zu, unclippedCoordPos=%zu, orientation=%zu, mchr=%zu, mateUnclippedCoordPos=%zu, rindex=%zu, mindex=%zu\n", read->readLb, read->readChromosome, read->unclippedCoordPos, read->orientation, read->mateChromosome, read->coordMatePos, read->indexAfterSort, read->mateIndexAfterSort);
           assert( (read->pair_num == 1) || (read->pair_num == 2));
   }    


    for (lnode_t *node = readEndsList->head; node != readEndsList->nil; node = node->next) {
        assert(node->read);
        readInfo *mate = getMateFromRead(htbl, node->read);

        if (!mate) {
            md_log_rank_error(rank, "Cannot found %s/2 of fingerprint %zu Please check if it's in the file or not.\n", node->read->Qname, node->read-> fingerprint, rank);
            exit(EXIT_FAILURE);
        }

        node->read->pairPhredScore += mate->phred_score;
        node->read->mateIndexAfterSort = mate->indexAfterSort;
        assert( (node->read->pair_num == 1) || (node->read->pair_num == 2));
        assert( (mate->pair_num == 1) || (mate->pair_num == 2));
        if (node->read->pair_num == 1) assert (mate->pair_num == 2); 
        if (node->read->pair_num == 2) assert (mate->pair_num == 1);
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
                markDuplicateFragments(nextCluster, htbl, totalDuplica, containsPairs);
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
        markDuplicateFragments(nextCluster, htbl, totalDuplica, containsPairs);
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

char *markDuplicateDiscordant (char *bufferReads, size_t readNum, char *header, MPI_Comm comm, char* chrName) {

    MPI_Comm previousComm = md_get_log_comm();
    md_set_log_comm(comm);

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
    exchangeExternFragDiscordant(fragList, readEndsList, htbl, interval, comm) ;
    md_log_debug("Finished to exchange mate in %f seconds\n", MPI_Wtime() - timeStamp);

    size_t fragList_size = 0;
    size_t readEnds_size = fragList_size = 0;
    
    for (lnode_t *node = fragList->head; node != fragList->nil; node = node->next) {
        readInfo *read = node->read;
        //md_log_trace("lb=%zu, chr=%zu, unclippedCoordPos=%zu, orientation=%zu, mchr=%zu, mateUnclippedCoordPos=%zu, rindex=%zu, mindex=%zu\n", read->readLb, read->readChromosome, read->unclippedCoordPos, read->orientation, read->mateChromosome, read->coordMatePos, read->indexAfterSort, read->mateIndexAfterSort);
        fragList_size++;
    }

    
    /* TRACE readEndsList */
    for (lnode_t *node = readEndsList->head; node != readEndsList->nil; node = node->next) {
        readInfo *read = node->read;
        //md_log_trace("lb=%zu, chr=%zu, unclippedCoordPos=%zu, orientation=%zu, mchr=%zu, mateUnclippedCoordPos=%zu, rindex=%zu, mindex=%zu\n", read->readLb, read->readChromosome, read->unclippedCoordPos, read->orientation, read->mateChromosome, read->coordMatePos, read->indexAfterSort, read->mateIndexAfterSort);
        

        readEnds_size++;
    }

   //fprintf(stderr, " TRACE 3 readEndsList size = %zu :::: fragList size =%zu \n", readEnds_size, fragList_size);


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


    readEnds_size = fragList_size = 0;
    
    for (lnode_t *node = fragList->head; node != fragList->nil; node = node->next) {
        readInfo *read = node->read;
        //md_log_trace("lb=%zu, chr=%zu, unclippedCoordPos=%zu, orientation=%zu, mchr=%zu, mateUnclippedCoordPos=%zu, rindex=%zu, mindex=%zu\n", read->readLb, read->readChromosome, read->unclippedCoordPos, read->orientation, read->mateChromosome, read->coordMatePos, read->indexAfterSort, read->mateIndexAfterSort);
        fragList_size++;
    }

    
    /* TRACE readEndsList */
    for (lnode_t *node = readEndsList->head; node != readEndsList->nil; node = node->next) {
        readInfo *read = node->read;
        //md_log_trace("lb=%zu, chr=%zu, unclippedCoordPos=%zu, orientation=%zu, mchr=%zu, mateUnclippedCoordPos=%zu, rindex=%zu, mindex=%zu\n", read->readLb, read->readChromosome, read->unclippedCoordPos, read->orientation, read->mateChromosome, read->coordMatePos, read->indexAfterSort, read->mateIndexAfterSort);
        readEnds_size++;
    }

   //fprintf(stderr, " TRACE 4 readEndsList size = %zu :::: fragList size =%zu \n", readEnds_size, fragList_size);



    readEnds_size = 0;
    /* TRACE readEndsList */
    for (lnode_t *node = readEndsList->head; node != readEndsList->nil; node = node->next) {
        readInfo *read = node->read;
        //md_log_trace("lb=%zu, chr=%zu, unclippedCoordPos=%zu, orientation=%zu, mchr=%zu, mateUnclippedCoordPos=%zu, rindex=%zu, mindex=%zu\n", read->readLb, read->readChromosome, read->unclippedCoordPos, read->orientation, read->mateChromosome, read->coordMatePos, read->indexAfterSort, read->mateIndexAfterSort);
        readEnds_size++;
    }

   //fprintf(stderr, " TRACE 4 readEndsList size = = %zu \n", readEnds_size);



    md_log_debug("Finished to sort list in %f seconds\n", MPI_Wtime() - timeStamp);

    md_log_info("Start to find duplicate ...\n");
    timeStamp = MPI_Wtime();
    int localDuplicates = 0, localOpticalDuplicates = 0, totalDuplicates = 0, totalOpticalDuplicates = 0;
    
    findDuplicaDiscordant(fragList, readEndsList, htbl, &localDuplicates, &localOpticalDuplicates, comm) ;
    


    md_log_info("Finished to find duplicates in %f seconds\n", totalDuplicates, totalOpticalDuplicates,  MPI_Wtime() - timeStamp);
    md_log_info("Gather duplicates results ...\n");
    
    timeStamp = MPI_Wtime();
    
    MPI_Barrier(comm);
    MPI_Reduce(&localDuplicates, &totalDuplicates, 1, MPI_INT, MPI_SUM, 0, comm);
    MPI_Reduce(&localOpticalDuplicates, &totalOpticalDuplicates, 1, MPI_INT, MPI_SUM, 0, comm);
    md_log_info("Found %d duplicates and %d optical duplicates in %f seconds\n", totalDuplicates, totalOpticalDuplicates,  MPI_Wtime() - timeStamp);

    md_log_debug("Start to write marked buffer...\n");
    timeStamp = MPI_Wtime();
    char *newBuff = writeBuff(samTokenLines, readArr, readNum);
    md_log_debug("Finished to mark buffer %f seconds\n", MPI_Wtime() - timeStamp);

    md_log_debug("Free data structures ...\n");
    timeStamp = MPI_Wtime();
    free(samTokenLines);
    /* free all reads and fragList */
    llist_destruct(fragList);
    /* free all nodes of readEndsList */
    llist_clear(readEndsList);
    /* free readEndsList */
    llist_destruct(readEndsList);
    free(intervalByProc);
    free(readNumByProc);
    free(readArr);
    freeChrInfo(&chr);
    freeLbInfo(&lb);

    hashTableDestroy(htbl);
    md_log_trace("Finished to free data structures %f seconds\n", MPI_Wtime() - timeStamp);
    md_log_debug("End to mark duplicates %f seconds\n", MPI_Wtime() - timeStart);
    md_set_log_comm(previousComm);
    return newBuff;
}

