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
     perfectHash.h
     
   Authors:
    Frederic Jarlier,   Institut Curie
    Firmin Martin,      Paris Descartes University
*/

#ifndef PERFECT_HASH_H
#define PERFECT_HASH_H


/**
 *   @date 2018 Feb 26
 *   @file perfectHash.h
 *   @author Firmin Martin
 *   @brief Perfect hashing implementation adapted to multiple processes context.
 *   @details Each process implement a perfect hashing table with a common main universal hash function. 
 *            Hence, processes can exchange read which have a locally given fingerprint.
 */

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>

// use -lcrypto to compile
// TODO:should be replaced by another open project
#include <openssl/md5.h>
#include <mpi.h>

#include "reads.h"
#include "mpiMD_utils.h"
#include "log.h"

typedef struct {
    unsigned long long m;
    int a;
    int b;
} hashParam;

typedef struct {
    hashParam h;
    readInfo **table;
} secTable;

typedef struct {
    hashParam h;
    int prime;
    int size;
    secTable **secTable;
} hashTable;

// Utils
int computeSpace(hashParam hp, int p, readInfo **arr, int size);
int haveCollision(hashParam hp, int prime, size_t *arr, int size);
void printPerfectHashTable(hashTable *htbl);

// Conversion
unsigned long long read2mateFP(readInfo* read) ;
unsigned long long read2Fingerprint(readInfo *read);
unsigned long long string2MD5(char *str);

// Construction 
void constructMainUnivHash(hashParam *hp, int *p, readInfo **arr, int size);
void constructMainUnivHashWithHp(hashParam hp, int *p, readInfo **arr, int size);
void constructSecTable(hashTable *htbl, readInfo **arr, int size);
void hashTableInit(hashTable *htbl, readInfo **arr, int size);
void hashTableInitWithHp(hashTable *htbl, hashParam hp, readInfo **arr, int size);
void shareHpAndConstructHtbl(hashTable *htbl, readInfo **arr, int size, MPI_Comm comm) ;

// get element
readInfo *getReadFromFingerprint(hashTable *htbl, size_t fingerprint);
readInfo *getReadFromQnameAndPairNum(hashTable *htbl, char *QName, int PairNum) ;
readInfo *getMateFromRead(hashTable *htbl, readInfo *read) ;

// Methods
void hashTableDestroy(hashTable *htbl);
void hashTableInsert(hashTable *htbl, readInfo *read) ;

// Hash function
int univHash(hashParam hp, int p, size_t k);

#endif
