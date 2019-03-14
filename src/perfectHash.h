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
int haveCollision(hashParam hp, int prime, unsigned long long *arr, int size);
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
readInfo *getReadFromFingerprint(hashTable *htbl, unsigned long long fingerprint);
readInfo *getReadFromQnameAndPairNum(hashTable *htbl, char *QName, int PairNum) ;
readInfo *getMateFromRead(hashTable *htbl, readInfo *read) ;

// Methods
void hashTableDestroy(hashTable *htbl);
void hashTableInsert(hashTable *htbl, readInfo *read) ;

// Hash function
int univHash(hashParam hp, int p, unsigned long long k);

#endif
