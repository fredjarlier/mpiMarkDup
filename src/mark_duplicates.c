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

inline CMInfo read2CM(readInfo *read) {

    /**
     * This function is only used in mate's rank determination, that is check_with_bruck == 2 case.
     * We initialize mateRank by -1 and fingerprint by the mate's fingerprint.
     */

    CMInfo coordMate = {.coordMatePos = read->coordMatePos, .mateRank = -1, .fingerprint = read2mateFP(read)};
    return coordMate;
}


/**
 * @date 2018 Feb 26
 * @param[in] token the string where we search tag
 * @param[in] tag searched tag in a read line
 * @return if tag is found, return the value of the tag, NULL otherwise.
 */

char *getReadTagValue(char *token, const char *tag) {

    char *substr = strstr(token, tag);

    if (substr) {
        char *p = substr + strlen(tag);
        return strndup(p, strlen(p));

    } else {
        return NULL;
    }

}

/**
 *   @date 2018 Mar 23
 *   @brief Determine where we can find mate.
 *   @details There is 4 cases
 *      1. check_with_bruck == 0 : mate is in the current buffer
 *      2. check_with_bruck == 1 : mate is in other process' buffer
 *      3. check_with_bruck == 2 : we don't know where is the mate, this should be verify later in ensureMateRank()
 *      4. check_with_bruck == 3 : should never happen
 *   @param[in, out] read a read
 *   @param[in] intervalByProc coordinate interval
 *   @param[in] comm markDuplicate communicator
 */

void checkMateBuffer(readInfo *read, Interval *intervalByProc, MPI_Comm comm) {

    int rank, num_proc;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &num_proc);

    // set check_with_bruck to 3 by default
    read->check_with_bruck = 3;

    /**
     * Current buffer case.
     * If mate position is in ] intervalByProc[rank].firstCoord,  intervalByProc[rank].lastCoord [,
     * we are certain that it's in the current buffer.
     * */

    if ( read->coordMatePos > intervalByProc[rank].firstCoord && read->coordMatePos < intervalByProc[rank].lastCoord ) {
        read->check_with_bruck = 0;
        goto check;
    }

    /**
     * Overlapped case.
     * Illustration :
     *
     * read->coordMatePos := 1081
     *
     * intervalByProc[0].firstCoord = 1000, intervalByProc[0].lastCoord = 1032
     * intervalByProc[1].firstCoord = 1035, intervalByProc[1].lastCoord = 1081
     * intervalByProc[2].firstCoord = 1081, intervalByProc[2].lastCoord = 1132
     *
     * In this case, we don't know whether the mate's rank is 1 or 2.
     * */

    for (int j = 1; j < num_proc ; j++ ) {
        if (read->coordMatePos == intervalByProc[j - 1].lastCoord && read->coordMatePos == intervalByProc[j].firstCoord) {
            read->check_with_bruck = 2;
            goto check;
        }
    }

    /*
     * Edge case.
     * If there is not overlapped case and mate's position is on current buffer edge,
     * we are certain that the mate is in current buffer.
     * */

    if (read->coordMatePos == intervalByProc[rank].lastCoord || read->coordMatePos == intervalByProc[rank].firstCoord) {
        read->check_with_bruck = 0;
        goto check;
    }

    /*
     * External case.
     * If mate's position doesn't meet above cases, it's external for sure.
     * */

    if ( read->coordMatePos < intervalByProc[rank].firstCoord ||  read->coordMatePos > intervalByProc[rank].lastCoord ) {
        read->check_with_bruck = 1;
        goto check;
    }

check:

    // This case should never happen.
    if (read->check_with_bruck == 3) {
        md_log_rank_warning(rank, "Couldn't determine read %s mate rank. Mate Position : %zu. This should never happen. Please report this issue.\n", read->Qname, read->coordMatePos);
    }

    // If mate is external, we can compute its rank directly.
    int found = 0;

    if (read->check_with_bruck == 1) {
        for (int j = num_proc - 1; j >= 0 ; j-- ) {
            if (read->coordMatePos >= intervalByProc[j].firstCoord && read->coordMatePos <= intervalByProc[j].lastCoord) {
                read->mate_rank = j;
                found = 1;
                break;

            }
        }

        /**
        *   If we don't found mate's rank, this means that mate's position is outside of all buffers position range.
        * This will cause implicitly failure to perfect hash table creation due to fictitious mate : we will have
        * fictitious mate and real mate in the same array, that is two reads having the same fingerprint (this implies collision).
        *   Hence, we exit the program and raise an error.
        */

        if (found == 0) {
            md_log_rank_error(rank, "Cannot found rank of mate %s which is supposed to be external of rank %d. Please check %s mate position.\n", read->Qname, rank, read->Qname);
            exit(EXIT_FAILURE);
        }
    }

}

/**
* @date 2018 Apr. 20
* @brief Compute reference length. Reference length is defined by the length of base which cigar operator is one of M/D/N/X/EQ.
* @param[in] cigar cigar string
* @return reference length
* @note See htsjdk/samtools/Cigar.java, function getReferenceLength().
*/

size_t getReferenceLength(char *cigar) {
    int refLen = 0, len = 0;
    char *q = cigar;

    while (*q) {
        len = 0;

        /* Compute the length of bases which the operator applies */
        while (isdigit(*q)) {
            len = 10 * len + (*q - '0');
            q++;
        }

        /* Cigar operator is one of M/D/N/X/EQ */
        if (*q && (*q == 'M' || *q == 'D' || *q == 'N' || *q == 'X' || *q == '=')) {
            refLen += len;
        }

        q++;
    }

    return refLen;
}

/**
 * @date 2018 Apr. 20
 * @brief Compute Alignment End. See htsjdk/samtools/SAMRecord.java, function getAlignmentEnd()
 * @details alignmentEnd := alignmentStart + referenceLength - 1
 * @param[in] alignmentStart alignment start coordinate
 * @param[in] cigar cigar string
 * @return alignment end coordinate
 */

size_t getAlignmentEnd(size_t alignmentStart, char *cigar) {
    return alignmentStart += getReferenceLength(cigar) - 1;
}

/**
 * @date 2018 Apr. 20
 * @param[in] alignmentStart alignment start coordinate
 * @brief Compute unclipped end. See htsjdk/samtools/SAMRecord.java, function getUnclippedEnd()
 * @param[in] cigar cigar string
 * @return unclipped end coordinate
 */

size_t getUnclippedEnd(size_t alignmentStart, char *cigar) {
    size_t alignmentEnd = getAlignmentEnd(alignmentStart, cigar);
    size_t clipping = 0;
    char *p = cigar + strlen(cigar) - 1;

    if (*p == 'S' || *p == 'H') {
        p--;

        while (isdigit(*p) && p != cigar) {
            p--;
        }

        p++;

        char *clippingString = strndup(p, cigar + strlen(cigar) - 1 - p);
        clipping = atoi(clippingString);
        free(clippingString);
    }

    return alignmentEnd += clipping;
}

/**
 * @date 2018 Apr. 20
 * @brief Compute unclipped start. See htsjdk/samtools/SAMUtils.java, function getUnclippedStart()
 * @param[in] alignmentStart alignment start coordinate
 * @param[in] cigar cigar string
 * @return unclipped start coordinate
 */
size_t getUnclippedStart (size_t alignmentStart, char *cigar) {


    size_t clipping = 0;

    char *p = cigar;

    while (isdigit(*p) && (*p != 'S' || *p != 'H')  && *p != '\0') {
        p++;
    }

    if (*p == 'S' || *p == 'H') {
        char *clippingString = strndup(cigar, p - cigar);
        clipping = atoi(clippingString);
        free(clippingString);
    }

    return alignmentStart -= clipping;
}

/**
 * @date 2018 Apr. 15
 * @brief Compute unclipped unclippedCoordPos by cigar string
 * @details Unclipped coordinate is defined by :
 *      - unclipped end, if read is reverse strand
 *      - unclipped start, if read is not reverse strand
 * @param[in, out] read a read
 * @note see Picard tools MarkDuplicates.java, function buildReadEnds()
 */

void fillUnclippedCoord(readInfo *read) {
    unsigned int readReverseStrand = readBits((unsigned int)read->valueFlag, 4);

    // compute unclipped coordinate
    if (readReverseStrand) {
        read->unclippedCoordPos = getUnclippedEnd(read->coordPos, read->cigar);

    } else {
        read->unclippedCoordPos = getUnclippedStart(read->coordPos, read->cigar);
    }
}

/**
 * @date 2018 Apr. 16
 * @brief Compute read optical location
 * @details read Qname should have the following format `identifier:lane:tile:x:y`.
 * @param[in, out] read a read
 * @return read physical location
 * @todo Here we just parse x and y field but it seems that lane and tile fields
 *       is also necessary.
 */

coord computePhysicalLocation(readInfo *read) {

    char *p = read->Qname + strlen(read->Qname) - 1;
    // y
    char *r = p;

    while (isdigit(*p) && *p != ':') {
        p--;
    }

    char *y = strndup(p + 1, r - p);
    // x
    r = --p;

    while (isdigit(*p) && *p != ':') {
        p--;
    }

    char *x = strndup(p + 1, r - p);

    read->physicalLocation.x = atoi(x);
    read->physicalLocation.y = atoi(y);
    coord physicalLocation = {.x = atoi(x), .y = atoi(y)};
    free(x);
    free(y);
    return physicalLocation;
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

readInfo *readParsing (char *sam_buff, Interval *intervalByProc, size_t readIndex, chrInfo *chr, lbInfo *lb,  MPI_Comm comm) {

    int rank, num_proc;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &num_proc);
    char *q = sam_buff;
    char *tokenCar;
    readInfo *read = calloc(1, sizeof(readInfo));
    int i = 0;

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
                getTokenTab(&q, sam_buff, &tokenCar);
                read->Qname = tokenCar;
                break;

            case 1: // Flag part
                getTokenTab(&q, sam_buff, &tokenCar);
                read->valueFlag  = atoi(tokenCar);
                
                assert(read->valueFlag);    
		        unsigned int readUnmapped = readBits((unsigned int)read->valueFlag, 2);
                unsigned int mateUnmapped = readBits((unsigned int)read->valueFlag, 3);
                unsigned int readReverseStrand = readBits((unsigned int)read->valueFlag, 4);
                unsigned int secondaryAlignment = readBits((unsigned int)read->valueFlag, 8);
                unsigned int supplementaryAlignment = readBits((unsigned int)read->valueFlag, 11);


                unsigned int firstInPair = readBits((unsigned int)read->valueFlag, 6);
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
                getTokenTab(&q, sam_buff, &tokenCar);

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
                getTokenTab(&q, sam_buff, &tokenCar);
                read->coordPos = atoll(tokenCar);
                free(tokenCar);
                break;

            case 5: // the CIGAR string
                getTokenTab(&q, sam_buff, &tokenCar);
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
                getTokenTab(&q, sam_buff, &tokenCar);

                /* mate chromosome same as read chromosome */
                if (strcmp(tokenCar, "=") == 0) {
                    read->mateChromosome = read->readChromosome;

                } else {
                    for (int j = 0; j < chr->chrNum; ++j) {
                        if (strcmp(chr->chrList[j], tokenCar) == 0) {
                            read->mateChromosome = j;
                            break;
                        }
                    }
                }

                free(tokenCar);
                break;


            case 7: // here we have the PNEXT, the position of the next read
                getTokenTab(&q, sam_buff, &tokenCar);
                read->coordMatePos = atoll(tokenCar);

                // we tell if the mate is in the buffer
                // we check if the mate's coordinates are outer the buffer limits
                         
                checkMateBuffer(read, intervalByProc, comm);
                free(tokenCar);

                break;

            case 10:
                getTokenTab(&q, sam_buff, &tokenCar);
                // here we got the quality string
                char *u = tokenCar;
                size_t score = 0;

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

                read->phred_score = min(read->phred_score, 0x7FFF / 2);
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
    assert(read->fingerprint);

    //Now we search in the flag for LB name
    while (getTokenTab(&q, sam_buff, &tokenCar)) {
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
 * @date 2018 Mar. 26
 * @brief Given a read, mark its mate as checked and assigned it a duplicate flag.
 * @param[in,out] htbl perfect hash table
 * @param[in] read a read
 * @param[in] d duplicate flag
 * @return 1 if the mate of @p read exists in @p htbl, 0 otherwise.
 */

int markMateDuplicateFlag(hashTable *htbl, readInfo *read, int d) {
    readInfo *mate = getMateFromRead(htbl, read);

    // found mate and it's not external
    if (mate && mate->external == 0) {
        mate->d       = d;
        mate->checked = 1;
        return 1;
    }

    return 0;
}


/**
 * @date 2018 Apr 12
 * @brief Compare two reads from different pair, decide whether or not two pairs belong to the same cluster
 * @param[in] pair1 first read in one pair
 * @param[in] pair2 first read in another pair
 * @param[in] isPaired pairs are paired, see function isPaired().
 * @note see function areComparableForDuplicates in MarkDuplicates.java
 */

int areComparableForDuplicates(readInfo *pair1, readInfo *pair2, const int isPaired) {

    int areComparable = 0;
    areComparable = pair1->readLb == pair2->readLb;

    if (areComparable == 1) {
        areComparable = pair1->readChromosome == pair2->readChromosome && pair1->unclippedCoordPos == pair2->unclippedCoordPos;

    }

    if (areComparable == 1) {
        if (isPaired) {
            areComparable = pair1->orientation == pair2->orientation;

        } else {
            areComparable = getOrientation(pair1, isPaired) == getOrientation(pair2, isPaired);
        }
    }

    if (areComparable && isPaired) {
        areComparable = pair1->coordMatePos == pair2->coordMatePos && pair1->mateChromosome == pair2->mateChromosome;
    }

    return areComparable;
}

/**
 * @date 2018 Apr 15
 * @brief Compare two reads from different pair, decide whether or not two pairs are close enough
 * @param[in] pair1 first read in one pair
 * @param[in] pair2 first read in another pair
 * @return 1 if pairs are close enough, 0 otherwise
 */

inline int closeEnough(readInfo *lhs, readInfo *rhs) {
    return lhs != rhs && lhs->readLb == rhs->readLb && abs(lhs->physicalLocation.x - rhs->physicalLocation.x) <= opticalDistance && abs(lhs->physicalLocation.y - rhs->physicalLocation.y) <= opticalDistance;
}

/**
 * @date 2018 Apr 15
 * @brief Compare two reads from different pair, decide whether or not two pairs are close enough
 * @param[in] cluster list of pairs in the same cluster
 * @param[in] best read has the best phred score
 * @note See function trackOpticalDuplicates() in AbstractMarkDuplicatesCommandLineProgram.java
 */

void trackOpticalDuplicates(llist_t *cluster, readInfo *best) {

    if (cluster->size < 2) {
        return;
    }

    //llist_readInfo_print(cluster);

    for (lnode_t *node = cluster->head; node != cluster->nil; node = node->next ) {
        node->read->isOpticalDuplicate = closeEnough(node->read, best);
        //printf("%s and %s %s\n", node->read->Qname, best->Qname, node->read->isOpticalDuplicate ? "are close enough" : "are not optical duplicates" );
    }

    for (lnode_t *node = cluster->head; node != cluster->nil; node = node->next) {
        readInfo *lhs = node->read;

        if (lhs == best) {
            continue;
        }

        for (lnode_t *cur = node->next; cur != cluster->nil; cur = cur->next) {
            readInfo *rhs = cur->read;

            if (rhs == best) {
                continue;
            }

            if (rhs->isOpticalDuplicate && lhs->isOpticalDuplicate) {
                continue;
            }

            if (closeEnough(lhs, rhs)) {
                readInfo *index = lhs->isOpticalDuplicate ? rhs : lhs;
                index->isOpticalDuplicate = 1;
            }

        }
    }
}

/**
 * @date 2018 Apr 12
 * @brief Mark as duplicate every pairs haven't the best phred score
 * @param[in] cluster a cluster
 * @param[in] htbl hash table
 * @param[in] totalDuplica amounts of duplicates
 * @note
      Documentation from MarkDuplicates (function markDuplicatePairs) :
      Takes a list of ReadEndsForMarkDuplicates objects and removes from it all objects that should
      not be marked as duplicates.  This assumes that the list contains objects representing pairs.
 */

void markDuplicatePairs(llist_t *cluster, hashTable *htbl, int *totalDuplica, int *totalOpticalDuplicate) {



    size_t bestScore = 0;
    readInfo *best = NULL;

    //we search the best node 
    for (lnode_t *node = cluster->head ; node != cluster->nil; node = node->next) {
        assert(node != node->next);
        
        // we test the phred score

        // in case the phred score is better
        if (node->read->pairPhredScore > bestScore || best == NULL) {
            bestScore = node->read->pairPhredScore;
            best = node->read;
            continue;
        }

        // in case phredscore are equals
        if (node->read->pairPhredScore == bestScore){
            
            // in tiebreak we use index of the read in the file after sorting
            // 0 should never happend in this case
            //from markduplicate
            //https://github.com/broadinstitute/picard/blob/master/src/main/java/picard/sam/markduplicates/MarkDuplicates.java
            //line 1039
            if ( best->indexAfterSort > node->read->indexAfterSort) {
                best = node->read;
                continue;
            }            
            
        }
    }

    // TODO: just find optical duplicates but not write anything in ouput
    // consider to add an optical duplicates flag
    if (opticalDistance > 0) {
        trackOpticalDuplicates(cluster, best);
    }

    for (lnode_t *node = cluster->head ; node != cluster->nil; node = node->next) {
        if (node->read != best) {

            node->read->d = 1;
            node->read->checked = 1;

            readInfo *mate = getMateFromRead(htbl, node->read);

            if (node->read->pair_num == 1) assert(mate->pair_num == 2);
            if (node->read->pair_num == 2) assert(mate->pair_num == 1);

            if (node->read->indexAfterSort != mate->indexAfterSort) {
                int mateIsDuplica = markMateDuplicateFlag(htbl, node->read, 1);
                (*totalDuplica) += 1 + mateIsDuplica;
            }

            if (node->read->isOpticalDuplicate) {
                (*totalOpticalDuplicate)++;
            }
        }
    }

    best->d = 0;
    best->checked = 1;
    markMateDuplicateFlag(htbl, best, 0);
}

/**
 * @date 2018 Apr 12
 * @brief Mark fragments as duplicates
 * @param[in] cluster a cluster
 * @param[in] htbl hash table
 * @param[in] totalDuplica amounts of duplicates
 * @param[in] containsPairs cluster contains pairs ?
 * @note
      Documentation from MarkDuplicates (function markDuplicateFragments) :

      Takes a list of ReadEndsForMarkDuplicates objects and removes from it all objects that should
      not be marked as duplicates.  This will set the duplicate index for only list items are fragments.
 */

void markDuplicateFragments(llist_t *cluster, hashTable *htbl, int *totalDuplica, const int containsPairs) {

    if (containsPairs) {
        for (lnode_t *node = cluster->head ; node != cluster->nil; node = node->next) {

            unsigned int readPaired = isPaired(node->read);

            if (!readPaired) {
                node->read->d = 1;
                (*totalDuplica)++;
            }
        }

    } else {

        size_t bestScore = 0;
        readInfo *best = NULL;

        for (lnode_t *node = cluster->head ; node != cluster->nil; node = node->next) {
            assert(node != node->next);

            if (node->read->phred_score > bestScore || best == NULL) {
                bestScore = node->read->phred_score;
                best = node->read;
            }
        }

        for (lnode_t *node = cluster->head ; node != cluster->nil; node = node->next) {
            if (node->read != best) {
                node->read->d = 1;
                node->read->checked = 1;
                (*totalDuplica) += 1;
            }
        }

        best->d = 0;
        best->checked = 1;
    }

    //fprintf(stderr, "in markDuplicateFragments totalDuplicates = %d \n", *totalDuplica);
}

/**
 * @date 2018 Jan 23
 * @brief return in outerDuplicate the amount of reads need to be checked by Bruck.
 * @details Count "check_with_bruck" field of all read in p_list.
 *          The result is assigned in outerDuplicate.
 *          These reads are the ones which others proc handle.
 * @param[in]   list the list.
 * @param[in]   externalMate result.
 * @return 0
 */

int countExternalMateInList(llist_t *list, size_t *externalMate) {

    lnode_t *node  = list->head;

    while (node != list->nil) {

        if ( node->read->check_with_bruck == 1 ) {
            (*externalMate)++;
        }

        node = node->next;
    }

    return 0;
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

void findDuplica(llist_t *fragList, llist_t *readEndsList, hashTable *htbl, int *totalDuplica, int *totalOpticalDuplicate, MPI_Comm comm) {

    int rank, num_proc;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &num_proc);
    
    llist_t *nextCluster = llist_create();

    //llist_readInfo_print(readEndsList);
    lnode_t *firstOfNextCluster = NULL;

    //fprintf(stderr, "rank %d Enter findDuplica \n", rank);

    for (lnode_t *node = readEndsList->head; node != readEndsList->nil; node = node->next) {

        if (firstOfNextCluster && areComparableForDuplicates(firstOfNextCluster->read, node->read, 1)) {
            //fprintf(stderr, "Call llist_append \n");
            llist_append(nextCluster, node->read);

        } else {

            if (nextCluster->size > 1) {                
                markDuplicatePairs(nextCluster, htbl, totalDuplica, totalOpticalDuplicate);
                //fprintf(stderr, "rank %d after findDuplica 1 totalDuplica = %d \n", rank, *totalDuplica);
            }

            //llist_readInfo_print(nextCluster);

            llist_clear(nextCluster);
            llist_append(nextCluster, node->read);
            firstOfNextCluster = node;

        }

    }
     //fprintf(stderr, "rank %d in findDuplica we go to nextCluster \n", rank);

    if (nextCluster->size > 1) {
        
        markDuplicatePairs(nextCluster, htbl, totalDuplica, totalOpticalDuplicate);
         //fprintf(stderr, "rank %d after findDuplica 2 totalDuplica = %d \n", rank, *totalDuplica);  
    }

    //llist_readInfo_print(nextCluster);
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
                //fprintf(stderr, "rank %d after markDuplicateFragments 1 totalDuplica = %d \n", rank, *totalDuplica);
                //fprintf(stderr, "call markDuplicateFragments \n");
                //llist_readInfo_print(nextCluster);
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
 * @param[in] list a fragments/paired-Ends list
 * @param[out] mates array of mates filled
 * @return number of 'check_with_bruck == 1' reads
 */

size_t getMateRankReadSizeBeforeBruck(llist_t *list, mateInfo **mates) {

    lnode_t *node  = list->head;
    size_t k = 0;

    while (node != list->nil) {

        readInfo *read = node->read;

        if (read->check_with_bruck == 1) {

            (*mates)[k].fingerprint = read->fingerprint;
            (*mates)[k].readLb = read->readLb;
            (*mates)[k].mateRank = read->mate_rank;
            (*mates)[k].phredScore = read->phred_score;
            (*mates)[k].indexAfterSort = read->indexAfterSort;
            (*mates)[k].unclippedCoordPos = read->unclippedCoordPos;
            (*mates)[k].coordPos = read->coordPos;
            (*mates)[k].coordMatePos = read->coordMatePos;
            (*mates)[k].orientation = read->orientation;
            (*mates)[k].valueFlag = read->valueFlag;
            (*mates)[k].pair_num = read->pair_num;
            k++;
        }

        node = node->next;
    }

    return k;
}


/**
 * @date 2018 Feb 26
 * @brief Compute send/recv count array and send/recv displacements array.
 * @param[in] numberOfReadsToSend number of reads to send
 * @param[out] scounts  Integer array, where entry i specifies the number of elements to send to rank i.
 * @param[out] sdispls  Integer array, where entry i specifies the displacement (offset from sendbuf, in units of sendtype) from which to send data to rank i.
 * @param[out] rcounts  Integer array, where entry j specifies the number of elements to receive from rank j.
 * @param[out] rdispls  Integer array, where entry j specifies the displacement (offset from recvbuf, in units of recvtype) to which data from rank j should be written.
 * @param[in] comm markDuplicate communicator
 * @return total received elements
 * @note sdispls, rcounts, rdispls is allocate in the function.
 */

int computeCountAndDispl(int *scounts, int **sdispls, int **rcounts, int **rdispls, MPI_Comm comm) {

    int rank, num_proc, totalReceived = 0;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &num_proc);

    *sdispls = calloc(num_proc, sizeof(int));
    *rcounts = calloc(num_proc, sizeof(int));
    *rdispls = calloc(num_proc, sizeof(int));

    // Compute send displacement
    for (int i = 1; i < num_proc; ++i) {
        (*sdispls)[i] = scounts[i - 1] + (*sdispls)[i - 1];
    }

    // Exchange send count to compute receive count
    MPI_Alltoall(scounts, 1, MPI_INT, *rcounts, 1, MPI_INT, comm);

    // Compute receive displacement
    for (int i = 1; i < num_proc; ++i) {
        (*rdispls)[i] = (*rcounts)[i - 1] + (*rdispls)[i - 1];
    }

    // Compute total received element
    for (int i = 0; i < num_proc; i++) {
        totalReceived += (*rcounts)[i];
    }

    return totalReceived;
}

/**
 * @date 2018 Feb 26
 * @brief Communication part of mark duplicate, exchange mates.
 * @param[out] matesByProc array of reads gather from others process
 * @param[in] mates array of mates to send sorted by rank
 * @param[in] numberOfReadsToSend number of reads to send
 * @param[in] comm markDuplicate communicator
 * @return amount of reads received
 * @todo benchmark MPI_Alltoallv/zeroCopyBruckv/Bruckv/OpenMPI bruck
 */

int exchangeAndFillMate(readInfo ***matesByProc, mateInfo *mates, size_t numberOfReadsToSend, MPI_Comm comm) {

    int rank, num_proc;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &num_proc);
    int *scounts, *sdispls, *rcounts, *rdispls;

    scounts = calloc(num_proc, sizeof(int));

    // Compute reads to send by rank
    for (int m = 0; m < numberOfReadsToSend; m++) {
        scounts[mates[m].mateRank]++;
    }

    int totalrecv = computeCountAndDispl(scounts, &sdispls, &rcounts, &rdispls, comm);

    /* Trace reads to send */
    for (int i = 0; i < numberOfReadsToSend; i++) {

        //assert (mates[i].pair_num == 1 || mates[i].pair_num == 2);

        md_log_rank_trace(rank, "send :: lb=%zu, mateRank=%zu, score=%zu, index=%zu, unclippedPos=%zu, Pos=%zu, matePos=%zu, orientation=%zu, fingerprint=%zu, valueFlag = %u, pair_num = %u \n",
        mates[i].readLb,
        mates[i].mateRank,
        mates[i].phredScore,
        mates[i].indexAfterSort,
        mates[i].unclippedCoordPos,
        mates[i].coordPos,
        mates[i].coordMatePos,
        mates[i].orientation,
        mates[i].fingerprint,
        mates[i].valueFlag,
        mates[i].pair_num);
    }

    /* End of trace */

    mateInfo *matesRecv = calloc(totalrecv, sizeof(mateInfo));

    // Create MPI type corresponding to mateInfo
    MPI_Datatype mate_type;
    createMateType(&mate_type);

    // Exchange mates
    MPI_Alltoallv(mates, scounts, sdispls, mate_type,
                  matesRecv, rcounts, rdispls, mate_type, comm);
    MPI_Type_free(&mate_type);


    /* Trace reads received */
    for (int i = 0; i < num_proc; i++) {
        for (int j = 0; j < rcounts[i] ; j++) {
            int index = rdispls[i] + j;
            assert (matesRecv[index].pair_num == 1 || matesRecv[index].pair_num == 2);


        md_log_rank_trace(rank, "received :: from rank =%d, lb=%zu, mateRank=%zu, score=%zu, index=%zu, unclippedPos=%zu, Pos=%zu, matePos=%zu, orientation=%zu, fingerprint=%zu, valueFlag = %u \n", i,
        matesRecv[index].readLb,
        matesRecv[index].mateRank,
        matesRecv[index].phredScore,
        matesRecv[index].indexAfterSort,
        matesRecv[index].unclippedCoordPos,
        matesRecv[index].coordPos,
        matesRecv[index].coordMatePos,
        matesRecv[index].orientation,
        matesRecv[index].fingerprint,
        matesRecv[index].valueFlag );
        }
    }

    /* End of trace */

    *matesByProc = calloc(totalrecv, sizeof(readInfo *));

    // Fill readInfo mates array by received mateInfo mates array

    for (int i = 0; i < num_proc; i++) {
        for (int j = 0; j < rcounts[i] ; j++) {
            int index = rdispls[i] + j;

            readInfo *mate = calloc(1, sizeof(readInfo));

            /**
             * Set mateInfo's string fields to NULL.
             * This permit us to free mates with freeRead later.
             */

            mate->Qname = NULL;
            mate->cigar = NULL;

            // Conversion of others fields
            mate->fingerprint = matesRecv[index].fingerprint;
            mate->readLb = matesRecv[index].readLb;
            mate->mate_rank = i;
            mate->phred_score = matesRecv[index].phredScore;
            mate->indexAfterSort = matesRecv[index].indexAfterSort;
            mate->unclippedCoordPos = matesRecv[index].unclippedCoordPos;
            mate->coordPos = matesRecv[index].coordPos;
            mate->coordMatePos = matesRecv[index].coordMatePos;
            mate->orientation = matesRecv[index].orientation;
            mate->valueFlag = matesRecv[index].valueFlag;
            mate->pair_num  = matesRecv[index].pair_num;
            // These mates are external and not checked.
            mate->external = 1;
            mate->checked = 0;

            (*matesByProc)[index] = mate;
        }
    }

    free(matesRecv);
    free(scounts);
    free(rcounts);
    free(sdispls);
    free(rdispls);
    free(mates);
    return totalrecv;
}

/**
 *   @date 2018 Feb 28
 *   @brief All processes send data to all processes, zero-copy version of MPI_Alltoall using Bruck's algorithm.
 *   @param[in] sendbuf     Starting address of send buffer (choice).
 *   @param[in] sendcount   Number of elements to send to each process (integer).
 *   @param[in] sendtype    Datatype of send buffer elements (handle).
 *   @param[out] recvbuf     Starting address of receive buffer (choice).
 *   @param[in] recvcount   Number of elements to receive from each process (integer).
 *   @param[in] recvtype    Datatype of receive buffer elements (handle).
 *   @param[in] comm        Communicator over which data is to be exchanged (handle).
 *   @note see J.L. Träff et al.[2014] paper, Implementing a Classic: Zero-copy All-to-all Communication with MPI Datatypes
 *   @bug there is an invalid write/read when we use, for instance, 7 processors
 */

void zeroCopyBruck(const void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf, int recvcount, MPI_Datatype recvtype, MPI_Comm comm) {

    int rank, num_proc;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &num_proc);
    int bits[num_proc], recvblocks[num_proc], sendblocks[num_proc];
    MPI_Aint recvindex[num_proc], sendindex[num_proc];
    MPI_Datatype recvtypes[num_proc], sendtypes[num_proc];
    bits[0] = 0;

    for (int j = 1; j < num_proc; j++) {
        bits[j] = bits[j >> 1] + (j & 0x1);
    }

    MPI_Aint lb, sendtotal, recvtotal;
    int recvsize;
    int intersize = 0;
    void *interbuf = malloc(sendcount * num_proc * sizeof(void));
    MPI_Datatype recvblocktype, sendblocktype;
    MPI_Type_get_extent(sendtype, &lb, &sendtotal);
    MPI_Type_get_extent(recvtype, &lb, &recvtotal);
    sendtotal *= sendcount;
    recvtotal *= recvcount;
    MPI_Type_size(recvtype, &recvsize);
    recvsize *= recvcount;

    unsigned int mask = 0xFFFFFFFF;

    for (int k = 1; k < num_proc; k <<= 1) {
        int b = 0;
        int j = k;

        do { // bit j set
            int sendrank = (rank - j + num_proc) % num_proc;
            int recvrank = (rank + j) % num_proc;

            if ((bits[j & mask] & 0x1) == 0x1) { // to recvbuf
                recvblocks[b] = recvcount;
                recvindex[b] = (MPI_Aint)((char *)recvbuf + recvrank * recvtotal);
                recvtypes[b] = recvtype;

                if ((j & mask) == j) { // from sendbuf
                    sendblocks[b] = sendcount;
                    sendindex[b] = (MPI_Aint)((char *)sendbuf + sendrank * sendtotal);
                    sendtypes[b] = sendtype;

                } else { // from intermediate
                    sendblocks[b] = recvsize;
                    sendindex[b] = (MPI_Aint)(interbuf + j * recvsize);
                    sendtypes[b] = MPI_BYTE;
                }

            } else { // to intermediate
                recvblocks[b] = recvsize;
                recvindex[b] = (MPI_Aint)(interbuf + j * recvsize);
                recvtypes[b] = MPI_BYTE;

                if ((j & mask) == j) { // from sendbuf
                    sendblocks[b] = sendcount;
                    sendindex[b] = (MPI_Aint)((char *)sendbuf + sendrank * sendtotal);
                    sendtypes[b] = sendtype;

                } else { // from recv
                    sendblocks[b] = recvcount;
                    sendindex[b] = (MPI_Aint)((char *)recvbuf + recvrank * recvtotal);
                    sendtypes[b] = recvtype;
                }
            }

            b++; // next element
            j++;

            if ((j & k) != k) {
                j += k;
            }
        } while (j < num_proc);

        MPI_Type_create_struct(b, sendblocks, sendindex, sendtypes, &sendblocktype);
        MPI_Type_commit(&sendblocktype);
        MPI_Type_create_struct(b, recvblocks, recvindex, recvtypes, &recvblocktype);
        MPI_Type_commit(&recvblocktype);

        int sendrank = (rank - k + num_proc) % num_proc;
        int recvrank = (rank + k) % num_proc;
        MPI_Sendrecv(MPI_BOTTOM, 1, sendblocktype, sendrank, k,
                     MPI_BOTTOM, 1, recvblocktype, recvrank, k,
                     comm, MPI_STATUS_IGNORE);

        MPI_Type_free(&recvblocktype);
        MPI_Type_free(&sendblocktype);
        mask <<= 1;
    }

    free(interbuf);
}

/**
 * @date 2018 Apr 15
 * @brief Insert a read in fragments list or pairs list
 * @param[in, out] l list
 * @param[in] read a read
 * @param[in] isFragmentList list is a fragments list ?
 * @note Note that due to unclipped coordinate, we need to do an insertion sort
 * @todo merge sort is implemented, function insertReadInList should be replace by llist_append 
 *       after we verify that all is fine, i.e. test llist_merge_sort().
 *       To do this, we can compare result of insertReadInList (uncommented code below) and
 *       llist_append().
 */


void insertReadInList(llist_t *l, readInfo *read, const int isFragmentList) {
    //llist_readInfo_print(l);
    llist_append(l, read);
    //if (l->size == 0) {
    //    llist_append(l, read);

    //} else if (compareRead(l->tail->read, read, isFragmentList) < 0) {

    //    llist_append(l, read);

    //} else {
    //    lnode_t *node = l->tail;

    //    while (node != l->nil && compareRead(node->read, read, isFragmentList) > 0) {
    //        node = node->prev;
    //    }

    //    llist_insert_node(l, node->next, read);

    //}
}


/**
 * @date 2018 Apr 15
 * @brief Build a paired-end and insert it in @p readEndsList list
 * @param[in] read1 first read
 * @param[in] read2 second read
 * @param[in] readEndsList paired-end list
 * @return read which is representing the pair
 */

readInfo *buildReadEnds(readInfo *read1, readInfo *read2, llist_t *readEndsList, int case_insert) {
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

    if ( read2->unclippedCoordPos >= read1->unclippedCoordPos ){

        //for discordant cases
        insertReadInList(readEndsList, read1, 0);
        return read2;
    } 

    else {

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

int parseLibraries(char *bufferReads, 
                    Interval *intervalByProc, 
                    llist_t *fragList, 
                    llist_t *readEndsList, 
                    readInfo ***readArr, 
                    char ***samTokenLines, 
                    size_t readNum, 
                    size_t readIndex, 
                    chrInfo *chr, 
                    lbInfo *lb, 
                    MPI_Comm comm,
                    int discordant_case) {

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

                    buildReadEnds(end, read, readEndsList, 0);

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


    // FOR DEBUG
    //  md_log_debug("Ensure reads mate rank which are susceptible in multiple processes ...\n");
    //  ensureMateRank(*readArr, intervalByProc, readNum, comm);

    free(bufferReads);
    assert(readCounter - readIndex == readNum);
    return readCounter;
}

/**
 *   @date 2018 Apr 7
 *   @brief Create MPI derived type corresponding to Interval
 *   @param[out] intervalType the MPI derived type
 */

static void createIntervalType(MPI_Datatype *intervalType) {

    int blocks[2] = {1, 1};
    MPI_Datatype types[2] = {MPI_UNSIGNED, MPI_UNSIGNED};
    MPI_Aint displacements[2];

    MPI_Aint uintex, lb;

    MPI_Type_get_extent(MPI_UNSIGNED, &lb, &uintex);

    displacements[0] = 0;
    displacements[1] = uintex;

    MPI_Type_create_struct(2, blocks, displacements, types, intervalType);
    MPI_Type_commit(intervalType);

}

/**
 *   @date 2018 Feb 26
 *   @brief gather interval from each proc
 *   @param[in] interval coordinate interval
 *   @param[in] comm markDuplicate communicator
 *   @return an array of coordinate interval by process
 */

Interval *gatherIntervalByProc(Interval interval, MPI_Comm comm) {
    int num_proc;
    MPI_Comm_size(comm, &num_proc);
    Interval *intervalByProc = calloc(num_proc,  sizeof(Interval));

    MPI_Datatype intervalType;
    createIntervalType(&intervalType);
    MPI_Allgather(&interval, 1, intervalType, intervalByProc, 1, intervalType, comm);
    MPI_Type_free(&intervalType);
    return intervalByProc;
}

/**
 *   @date 2018 Feb 26
 *   @brief Convert a flag value to string
 *   @param[in] flagValue a read flag
 *   @return associated string of @p flagValue
 */

char *flag2string(int flagValue) {
    char *valueString = calloc(5, sizeof(char));
    sprintf(valueString, "%u", flagValue);
    return valueString;
}

/**
 *   @date 2018 Feb 26
 *   @brief Given a read line, return a new one with duplicate flag
 *   @param[in] oldReadLine a read line
 *   @param[out] read the associated read
 *   @return a new read line
 */

char *addDuplicateFlag(char *oldReadLine, readInfo *read) {
    char *newLine = calloc((strlen(oldReadLine) + 4 + 1), sizeof(char));
    char *charTab = strchr(oldReadLine, '\t');
    strncpy(newLine, oldReadLine, (charTab - oldReadLine) / sizeof(char) + 1);
    char *flagString = flag2string(read->valueFlag + 1024);
    strcat(newLine, flagString);
    charTab = strchr(charTab + 1, '\t');
    strcat(newLine, charTab);
    free(flagString);
    free(oldReadLine);
    return newLine;
}

/**
*   @date 2018 Feb 26
*   @brief return the new reads buffer with duplicate flag
*   @param[in] samTokenLines an array of reads lines
*   @param[in] readArr an array of reads which have been handled
*   @param[in] readNum amount of total reads
*   @return new reads buffer which duplicate reads are marked
*   @todo consider to print message every 1M reads, just like markDuplicates
*/

char *writeBuff(char **samTokenLines, readInfo **readArr, size_t readNum) {
    size_t bufferSize = 1;
    char *newBuff = NULL, *p = NULL;
    size_t splitLineSize = (readNum / SPLIT_FACTOR) + 1;
    size_t splitNumber = 0;
    int percentage = 1, increment = 10;
    float progression = 0;

    for (int i = 0; i < readNum; i++) {
        progression = 100 * ((float)i / readNum);
        size_t curMaxLine = splitLineSize * splitNumber;

        if (i >= curMaxLine) {
            curMaxLine = curMaxLine + splitLineSize > readNum ? readNum : curMaxLine + splitLineSize;
            splitNumber++;
            size_t curSizeOfNewBuff = newBuff ? strlen(newBuff) : 0;

            for (int j = i; j < curMaxLine ; j++) {
                bufferSize += strlen(samTokenLines[j]) + 4;
            }

            newBuff = realloc(newBuff, bufferSize * sizeof(char));
            newBuff[curSizeOfNewBuff] = '\0';

            p = newBuff;
        }

        if (!readArr[i] || !readArr[i]->d) {
            p = strapp(p, samTokenLines[i]);
            free(samTokenLines[i]);

        } else if (readArr[i]->d) {
            char *newLine = addDuplicateFlag(samTokenLines[i], readArr[i]);
            p = strapp(p, newLine);
            free(newLine);
        }

        if (progression > percentage) {
            md_log_debug("Writed %.0f%% of reads\n", progression);
            percentage += increment;
        }
    }

    return newBuff;
}

/**
 *   @date 2018 Feb 26
 *   @brief[in] Given a read line, return its POS field
 *   @param samLine a read line
 *   @return the field POS of the read
 */

int getPosFromLine(char *samLine) {
    char *saveptr, *tokenTab = strtok_r(samLine, "\t", &saveptr);
    int i = 1;

    while (i++ < 4) {
        tokenTab = strtok_r(NULL, "\t", &saveptr);
    }

    return atoi(tokenTab);
}

/**
 *   @date 2018 Feb 26
 *   @brief Given a read buffer and the amount of reads, return coordinate interval.
 *   @param[in] bufferReads reads buffer (i.e. SAM file without header)
 *   @param[in] readNum total amount of reads
 *   @return the coordinate interval
 */

Interval getIntervalFromBuffer(char *bufferReads, size_t readNum) {
    Interval interval = {.firstCoord = 0, .lastCoord  = 0};
    char *firstLine, *lastLine, *offset = bufferReads;
    getLine(&offset, bufferReads, &firstLine);
    interval.firstCoord = getPosFromLine(firstLine);
    free(firstLine);

    if (readNum > 1) {
        offset = bufferReads + strlen(bufferReads) - 2;

        while (*offset != '\n') {
            offset--;
        }

        offset++;
        getLine(&offset, bufferReads, &lastLine);
        interval.lastCoord = getPosFromLine(lastLine);
        free(lastLine);

    } else {
        interval.lastCoord = interval.firstCoord;
    }
    
    assert(interval.firstCoord <= interval.lastCoord);
    return interval;
}

/**
 * @date 2018 Mar 23
 * @brief return the amount of reads need to be checked by Bruck.
 * @details Count "check_with_bruck" field of all read in a readInfo array
 * @param[in] readArr a read array
 * @param[in] readNum number of read
 * @param[in] check_with_bruck flag to give stat
 * @return number of read which have field check_with_bruck equal to @p checkWithBruck
 */

int countExternalMateInArray(readInfo **readArr, size_t readNum, unsigned int checkWithBruck) {
    int counter = 0;

    for (int i = 0; i < readNum; i++) {
        if (readArr[i] && readArr[i]->check_with_bruck == checkWithBruck) {
            counter++;
        }
    }

    return counter;
}


/**
 * @date 2018 Mar 23
 * @brief Fill reads array with external mates.
 * @param[in] readArr array of read
 * @param[out] readArrWithExternal array of read with external mates
 * @param[in] readNum total amount of read lines
 * @note This function is only used for perfect hashing, since perfect hashing key set must be pre-computed.
 *         We keep slots for externals mates and insert them later.
 */

int fillReadAndFictitiousMate(readInfo **readArr, readInfo ***readArrWithExternal, size_t readNum) {

    size_t externalNum = countExternalMateInArray(readArr, readNum, 1);

    //fprintf(stderr, " in fillReadAndFictitiousMate : externalMate = %zu \n", externalNum);
    
    size_t totalRead = readNum + externalNum;
    *readArrWithExternal = malloc(totalRead * sizeof(readInfo *));
    int current = 0;

    for (int i = 0; i < readNum; i++) {
        if (readArr[i] && readArr[i]->check_with_bruck == 1) {

            (*readArrWithExternal)[current++] = readArr[i];
            readInfo *mate = calloc(1, sizeof(readInfo));
            /* clone read */
            mate = cloneRead(readArr[i]);            

            //md_log_rank_trace(rank, "fictitious mate Qname = %s\n", mate->Qname);
            mate->external = 1;
            /* just give it pair number flag to compute fingerprint */
            mate->fingerprint = read2mateFP(readArr[i]);
            (*readArrWithExternal)[current++] = mate;

            if (readArr[i]->pair_num == 1 ) mate->pair_num = 2;
            if (readArr[i]->pair_num == 2 ) mate->pair_num = 1;    

            assert(mate->valueFlag);

        } else {
            (*readArrWithExternal)[current++] = readArr[i];
            //assert(readArr[i]->valueFlag);
            //assert((*readArrWithExternal)[current++]->valueFlag);
        }
    }

    assert(current == totalRead);
    return totalRead;
}


inline unsigned int readFlag2MateFlag(unsigned int readFlag) {
    unsigned int mateFlag = readFlag;
    switchBits(&mateFlag, 2, 3);
    switchBits(&mateFlag, 4, 5);
    switchBits(&mateFlag, 6, 7);
    return mateFlag;
}


/**
 *   @date 2018 Feb 26
 *   @brief Exchange extern fragments and complete pairs information (phredScore and mate index)
 *   @param[out] fragList fragments list
 *   @param[out] readEndsList paired-end list
 *   @param[out] htbl perfect hash table
 *   @param[in] comm markDuplicate communicator
 */

void exchangeExternFrag(llist_t *fragList, 
                        llist_t *readEndsList, 
                        hashTable *htbl, 
                        Interval interval, 
                        MPI_Comm comm,
                        int discordant_case) {

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

    md_log_rank_debug(rank, "[mpiMD][exchangeExternFrag] Received %d mates, fragList size = %d\n", totalrecv, fragList->size);

    //test if we have nothing to do we return
    if (totalrecv == 0) return; 

    /*
     *     FOR DEBUG
     *
     *
    for (lnode_t *node = fragList->head; node != fragList->nil; node = node->next) {
        readInfo *read = node->read;
        //md_log_trace("lb=%zu, chr=%zu, unclippedCoordPos=%zu, orientation=%zu, mchr=%zu, mateUnclippedCoordPos=%zu, rindex=%zu, mindex=%zu\n", read->readLb, read->readChromosome, read->unclippedCoordPos, read->orientation, read->mateChromosome, read->coordMatePos, read->indexAfterSort, read->mateIndexAfterSort);
        assert( (read->pair_num == 1) || (read->pair_num == 2));
    }

    for (lnode_t *node = readEndsList->head; node != readEndsList->nil; node = node->next) {
        readInfo *read = node->read;
        //md_log_trace("lb=%zu, chr=%zu, unclippedCoordPos=%zu, orientation=%zu, mchr=%zu, mateUnclippedCoordPos=%zu, rindex=%zu, mindex=%zu\n", read->readLb, read->readChromosome, read->unclippedCoordPos, read->orientation, read->mateChromosome, read->coordMatePos, read->indexAfterSort, read->mateIndexAfterSort);
        assert( (read->pair_num == 1) || (read->pair_num == 2));
    }    
    */


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

                        readInfo *insertedRead = NULL;

                        buildReadEnds(matesByProc[i], node->read, readEndsList, 1);    

                        insertReadInList(fragList, matesByProc[i], 1) ;
                        matesByProc[i]->Qname = strdup(node->read->Qname);

                        assert ((node->read->pair_num == 1 ) || (node->read->pair_num ==2));
                        
                        if (node->read->pair_num == 1) assert (matesByProc[i]->pair_num == 2); 
                        if (node->read->pair_num == 2) assert (matesByProc[i]->pair_num == 1);
                        assert( (matesByProc[i]->pair_num == 1) || (matesByProc[i]->pair_num == 2));

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
        assert( (node->read->pair_num == 1) || (node->read->pair_num == 2));
        assert( (mate->pair_num == 1) || (mate->pair_num == 2));
        if (node->read->pair_num == 1) assert (mate->pair_num == 2); 
        if (node->read->pair_num == 2) assert (mate->pair_num == 1);
    }

    /* free external mates in current library, others reads are free in destroyLBList */

    free(matesByProc);
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

char *markDuplicate (char *bufferReads, size_t readNum, char *header, MPI_Comm comm, char *chrName) {

    MPI_Comm previousComm = md_get_log_comm();
    md_set_log_comm(comm);

    int rank, num_proc, discordant_case;

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
    parseLibraries(bufferReads, 
                    intervalByProc, 
                    fragList, 
                    readEndsList, 
                    &readArr, 
                    &samTokenLines, 
                    readNum, 
                    readIndexOffset, 
                    &chr,  
                    &lb, 
                    comm,
                    discordant_case);

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

   
    md_log_debug("Start to exchange mates ...\n");
    timeStamp = MPI_Wtime();
    exchangeExternFrag(fragList, readEndsList, htbl, interval, comm, discordant_case) ;
    md_log_debug("Finished to exchange mate in %f seconds\n", MPI_Wtime() - timeStamp);

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
 
    md_log_debug("Finished to sort list in %f seconds\n", MPI_Wtime() - timeStamp);

    md_log_info("Start to find duplicate ...\n");
    timeStamp = MPI_Wtime();
    int localDuplicates = 0, localOpticalDuplicates = 0, totalDuplicates = 0, totalOpticalDuplicates = 0;
    
    findDuplica(fragList, readEndsList, htbl, &localDuplicates, &localOpticalDuplicates, comm) ;
    


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

