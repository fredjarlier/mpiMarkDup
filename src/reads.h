#ifndef READS_H
#define READS_H

#include <stdlib.h>
#include <string.h>

#include <mpi.h>

#include "common.h"
#include "mpiMD_utils.h"

typedef struct {
    int x;
    int y;
} coord;

typedef enum {
    F = 0, R = 1, FF = 2, FR = 3, RR = 4, RF = 5
} orientation;

typedef struct readInfo {
    

    char *Qname;                        /**< read name, TODO:we don't need to store Qname, we just use it to compte physical location for optical duplicates */
    char *cigar;                        /**< cigar string, TODO:we don't need to store cigar, we just use it to compute unclipped coordinate */
    
                 /**< mate chromosome number */
                 
    unsigned int d: 1;                    /**< tell if the read is a duplicate */
    unsigned int checked: 1;              /**< tell the read has already been checked */
    unsigned int check_with_bruck: 2;     /**< the mate is in another buffer we need bruck to check */
    unsigned int external: 1;             /**< it means this read comes from another buffer we don(t mark it as duplicate */
    unsigned int isOpticalDuplicate : 1;  /**< read is optical duplicate ? */
    unsigned int valueFlag;                   /**< flag of the read */
    unsigned int pair_num;               /**< tell if the read is first(1) or second in the pair(2) */    

    unsigned long long fingerprint;    /**< MD5 of the Qname + number in pair */

    size_t mate_rank;                   /**< the rank where the mate is supposed to be */
    size_t readLb;                      /**< read library number */
    size_t coordPos;                    /**< clipped read position (coordinate in sam file) */
    size_t unclippedCoordPos;           /**< unclipped read position */
    size_t coordMatePos;                /**< clipped mate position (coordinate in sam file) */
    size_t phred_score;                 /**< read phred score */
    size_t pairPhredScore;              /**< paired-end phred score */
    size_t readChromosome;              /**< read chromosome number */
    size_t mateChromosome; 
    size_t indexAfterSort;              /**< read index in file */
    size_t mateIndexAfterSort;
              /**< mate index in file */
    coord physicalLocation;             /**< read physical location (in QNAME) */
    orientation orientation;            /**< fragment or paired-end orientation (reverse strand + first/second in pair) */

} readInfo;

/**
 * Note that there are several place to modify when we want edit mateInfo fields :
 * - reads.c, function createMateType() : creation of MPI type corresponding to mateInfo
 * - mark_duplicate.c, function getMateRankReadSizeBeforeBruck() : fill mate to send
 * - mark_duplicate.c, function exchangeAndFillMate() : communication part and fill received readInfo array
 */

typedef struct {

    unsigned int pair_num;
    unsigned int valueFlag;
    unsigned long long fingerprint;   /**< mate fingerprint */
    size_t readLb;
    size_t mateRank;                  /**< mate rank */
    size_t phredScore;                /**< mate phredScore */
    size_t indexAfterSort;
    size_t unclippedCoordPos;
    size_t coordPos;
    size_t coordMatePos;
    size_t orientation;
    
    
} mateInfo;

//Coordinate Mate Info
typedef struct {
    size_t coordMatePos;            /**< mate coordinate (PNext) */
    int mateRank;                   /**< mate rank */
    unsigned long long fingerprint; /**< mate fingerprint */
} CMInfo;


int freeRead(readInfo *read) ;
readInfo *cloneRead(readInfo *read) ;
void createMateType(MPI_Datatype *mate_type) ;
void createCMType(MPI_Datatype *CMtype) ;
int cmpMateRank(const void *a, const void *b) ;
orientation getOrientation (readInfo *read, const int isPaired) ;
int isPaired(readInfo *read) ;
#endif /* ifndef READS_H */
