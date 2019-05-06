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
     mark_duplicates.h
     
   Authors:
    Frederic Jarlier,   Institut Curie
    Firmain Martin,     Paris Descartes University
*/

#ifndef MARKDUPLICATES_H
#define MARKDUPLICATES_H

#include <stdlib.h>
#include <assert.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <limits.h>

#include <mpi.h>

#include "log.h"
#include "mpiMD_utils.h"
#include "perfectHash.h"
#include "reads.h"
#include "parser.h"
#include "llist.h"
#include "createLBList.h"

typedef struct {
    unsigned int firstCoord;
    unsigned int lastCoord;
} Interval;


int areComparableForDuplicates(readInfo *pair1, readInfo *pair2, const int isPaired);
void markDuplicateFragments(llist_t *cluster, hashTable *htbl, int *totalDuplica, const int containsPairs);
void markDuplicatePairs(llist_t *cluster, hashTable *htbl, int *totalDuplica, int *totalOpticalDuplicate);
int countExternalMateInList(llist_t *list, size_t *externalMate); 
readInfo *readParsing (char *sam_buff, Interval *intervalByProc, size_t readIndex, chrInfo *chr, lbInfo *lb,  MPI_Comm comm);
void insertReadInList(llist_t *l, readInfo *read, const int isFragmentList);
int ComputePercentageIncrement(size_t readNum, int readSecRatio, int second);
char* addDuplicateFlag(char *oldReadLine, readInfo *read);
void checkMateBuffer(readInfo *read, Interval *intervalByProc, MPI_Comm comm);
readInfo* cloneRead(readInfo *read);
int cmpMateRank(const void *a, const void *b);
int computeCountAndDispl(int *scounts, int **sdispls, int **rcounts, int **rdispls, MPI_Comm comm);
int countExternalMateInArray(readInfo **readArr, size_t readNum, unsigned int checkWithBruck);
void createCMType(MPI_Datatype *CMtype);
void createMateType(MPI_Datatype *mate_type);
void createMateType2(MPI_Datatype *mate_type);
void ensureMateRank(readInfo **readArr, Interval *intervalByProc, size_t readNum, MPI_Comm comm);
int exchangeAndFillMate(readInfo ***matesByProc, mateInfo *mates, size_t numberOfReadsToSend, MPI_Comm comm);
int fillReadAndFictitiousMate(readInfo **readArr, readInfo ***readArrWithExternal, size_t readNum);
int fillReadLBValue(char *token, readInfo *read);
char* flag2string(int flagValue);
Interval* gatherIntervalByProc(Interval interval, MPI_Comm comm);
Interval getIntervalFromBuffer(char *bufferReads, size_t readNum);
size_t getMateRankReadSizeBeforeBruck(llist_t *list, mateInfo **mates) ;
int getPosFromLine(char *samLine);
char* getReadTagValue(char *token, const char *tag);
char *markDuplicate (char *bufferReads, size_t readNum, char *header, MPI_Comm comm, char *chrName) ;
int markMateDuplicateFlag(hashTable *htbl, readInfo *read, int d);
char* writeBuff(char **samTokenLines, readInfo **readArr, size_t readNum);
void zeroCopyBruck(const void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf, int recvcount, MPI_Datatype recvtype, MPI_Comm comm);
CMInfo read2CM(readInfo *read) ;
#endif
