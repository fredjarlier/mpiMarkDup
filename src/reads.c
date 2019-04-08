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
     llist.c
     
   Authors:
    Frederic Jarlier,   Institut Curie
    Firmin Martin,     Paris Descartes University
*/

#include "reads.h"


/**
 *   @date 2018 Mar 23
 *   @brief free a read
 *   @param[out] read a read
 */

int freeRead(readInfo* read) {

    if (!read){
        return -1;
    }

    if (read->Qname) {
        free(read->Qname);
        read->Qname = NULL;
    }

    if (read->cigar) {
        free(read->cigar);
        read->cigar = NULL;
    }

    free(read);
    return 0;
}

/**
 *   @date 2018 Feb 26
 *   @brief Compare two mates by their rank
 *   @param[in] a a mate
 *   @param[in] b another mate
 *   @return the difference between @p a and @p mate rank.
 */

int cmpMateRank(const void *a, const void *b) {
    mateInfo *mateA = (mateInfo *) a;
    mateInfo *mateB = (mateInfo *) b;
    return mateA->mateRank - mateB->mateRank;
}

/**
 *   @date 2018 Mar 23
 *   @brief Create MPI derived type corresponding to CMInfo
 *   @param[out] CMtype the MPI derived type
 */

void createCMType(MPI_Datatype *CMtype) {

    int blocks[3] = {1, 1, 1};
    MPI_Datatype types[3] = {MPI_SIZE_T, MPI_SIZE_T, MPI_UNSIGNED_LONG_LONG};
    MPI_Aint displacements[3] = {0};
    MPI_Aint sizeEx, lb;

    MPI_Type_get_extent(MPI_SIZE_T, &lb, &sizeEx);

    for (int i = 1; i < 3; i++) {
        displacements[i] = i * sizeEx;
    }

    MPI_Type_create_struct(3, blocks, displacements, types, CMtype);
    MPI_Type_commit(CMtype);
}

/**
 *   @date 2018 Feb 26
 *   @brief Create MPI derived type corresponding to mateInfo
 *   @param[out] mate_type the MPI derived type
 *   @todo we exchange more field than we need for commodity, remove them for release.
 */

void createMateType(MPI_Datatype *mate_type) {
    
    int blocks[11] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};

    MPI_Datatype types[11] = {MPI_UNSIGNED_LONG_LONG, MPI_SIZE_T, 
        MPI_SIZE_T, MPI_SIZE_T, MPI_SIZE_T, 
            MPI_SIZE_T, MPI_SIZE_T, MPI_SIZE_T, 
                MPI_SIZE_T, MPI_UNSIGNED, MPI_UNSIGNED };

    MPI_Aint displacements[11] = {0};
    MPI_Aint sizeEx, unintex, ulltex, lb, unInt;

    MPI_Type_get_extent(MPI_SIZE_T, &lb, &sizeEx);
    //MPI_Type_get_extent(MPI_UNSIGNED, &lb, &unintex);
    MPI_Type_get_extent(MPI_UNSIGNED_LONG_LONG, &lb, &ulltex);
    MPI_Type_get_extent(MPI_UNSIGNED, &lb, &unInt);

    displacements[0] = 0;
    displacements[1] = displacements[0] + ulltex;
    displacements[2] = displacements[1] + sizeEx;
    displacements[3] = displacements[2] + sizeEx;
    displacements[4] = displacements[3] + sizeEx;
    displacements[5] = displacements[4] + sizeEx;
    displacements[6] = displacements[5] + sizeEx;
    displacements[7] = displacements[6] + sizeEx;
    displacements[8] = displacements[7] + sizeEx;
    displacements[9] = displacements[8] + unInt;   // for valueFlag
    displacements[10] = displacements[9] + unInt;  // for pairnum   

    MPI_Type_create_struct(11, blocks, displacements, types, mate_type);
    MPI_Type_commit(mate_type);
}



/**
 * @date 2018 Mar 23
 * @brief Clone a read
 * @param[in] read a readArrWithExternal
 * @return a read cloned
 */

readInfo *cloneRead(readInfo *read) {
    readInfo *clonedRead = calloc(1, sizeof(readInfo));
    clonedRead->fingerprint = read->fingerprint;
    clonedRead->Qname   = strdup(read->Qname);
    clonedRead->cigar   = strdup(read->cigar);
    clonedRead->mate_rank = read->mate_rank;
    clonedRead->readLb = read->readLb;
    clonedRead->readChromosome = read->readChromosome;
    clonedRead->mateChromosome = read->mateChromosome;
    clonedRead->d = read->d;
    clonedRead->checked = read->checked;
    clonedRead->check_with_bruck = read->check_with_bruck;
    clonedRead->external = read->external;
    clonedRead->isOpticalDuplicate = read->isOpticalDuplicate;
    clonedRead->coordPos = read->coordPos;
    clonedRead->unclippedCoordPos = read->unclippedCoordPos;
    clonedRead->coordMatePos = read->coordMatePos;
    clonedRead->phred_score = read->phred_score;
    clonedRead->pairPhredScore = read->pairPhredScore;
    clonedRead->valueFlag = read->valueFlag;
    clonedRead->indexAfterSort = read->indexAfterSort;
    clonedRead->mateIndexAfterSort = read->mateIndexAfterSort;
    clonedRead->physicalLocation = read->physicalLocation;
    clonedRead->orientation = read->orientation;
    clonedRead->pair_num = read->pair_num;

    return clonedRead;
}

/**
 * @date 2018 Apr 12
 * @brief Compute orientation (reverse strand x first/second in pair) for a given read
 * @param[in] read1 first read in one pair
 * @param[in] isPaired read is paired ?
 * @note Documentation in MarkDuplicates (function getOrientationByte in class sam.util.ReadEnds) :
 *       Returns a single byte that encodes the orientation of the two reads in a pair.
 */

orientation getOrientation (readInfo *read, const int isPaired) {

    unsigned int readReverseStrand = readBits((unsigned int)read->valueFlag, 4);

    if (!isPaired) {
        return readReverseStrand ? R : F;

    } else {
        unsigned int mateReverseStrand = readBits((unsigned int)read->valueFlag, 5);

        int read1NegativeStrand = readReverseStrand;;
        int read2NegativeStrand = mateReverseStrand;;

        if (read1NegativeStrand) {
            if (read2NegativeStrand) {
                return RR;

            } else {
                return RF;
            }

        } else {
            if (read2NegativeStrand) {
                return FR;

            } else {
                return FF;
            }
        }
    }
}

/**
 * @date 2018 Apr 15
 * @brief Verify if a read is paired
 * @details Note that a read is not paired if its mate is unmapped according markDuplicate
 * @param[in] read
 */

inline int isPaired(readInfo *read) {
    unsigned int readPaired = readBits((unsigned int)read->valueFlag, 0);
    unsigned int mateUnmapped = readBits((unsigned int)read->valueFlag, 3);
    return readPaired && !mateUnmapped;

}
