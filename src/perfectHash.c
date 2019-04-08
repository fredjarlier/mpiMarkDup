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
     perfectHash.c
     
   Authors:
    Frederic Jarlier,   Institut Curie
    Firmin Martin,     Paris Descartes University
*/

#include "perfectHash.h"

/**
 *   @date 2018 Feb 26
 *   @file perfectHash.c
 *   @brief Perfect hashing implementation adapted to multiple processes context.
 *   @details Each process implement a perfect hashing table with a *common main universal hash function*. 
 *            Hence, processes can exchange read which have a locally given fingerprint.
 *
 *            Illustration : 
 *
   @verbatim
               main table            universal        second level table
                                     hash function    
                                     parameters
                                     of each second
                                     table
               
                                     m    a    b
              +---------+          +----+----+----++----+
           0  |    x-------------->| 1  | 0  | 0  || r1 |
              +---------+          +----+----+----++----+
           1  |    /    |                                    
              +---------+          +----+----+----++----+----+----+----+
           2  |    x-------------->| 4  | 10 | 18 || r1 | r2 | // | // |
              +---------+          +----+----+----++----+----+----+----+
           3  |    /    |
              +---------+
           4  |    /    |
              +---------+          +----+----+----++----+
           5  |    x-------------->| 1  | 0  | 0  || r3 |
              +---------+          +----+----+----++----+
           6  |    /    |
              +---------+          +----+----+----++----+----+----+----+----+----+----+----+----+
           7  |    x-------------->| 9  | 23 | 88 || r4 | // | // | r5 | // | // | // | // | r6 |
              +---------+          +----+----+----++----+----+----+----+----+----+----+----+----+
           8  |    /    |
              +---------+
   @endverbatim
 *
 *    Algorithm to construct perfect hash table in a single process  
 *    ===============================================================
 *
 *    Pseudocode
 *    ----------
 *
 * @verbatim 
         Repeat
           choose a main universal hash function among H_{p, m}
         Until \sum_{i=0}^{n-1} n^2_i <= 2n
         For each slot i 
             Repeat
               choose an universal hash function among H_{p, n_i^2}
             Until there is not collision
   @endverbatim 
 *
 *    Algorithm to construct perfect hash table in multiple processes  
 *    ===============================================================
 *
 *    This algorithm permit us to exchange reads between processes because
 *    all main table share the same hash function.
 * 
 *    Pseudocode
 *    ----------
 *
 * @verbatim 
         Gather the max size of table through process (MPI_Allreduce)
         If rank == 0 Then
             a = b = rand() + 1
             find p > m such that \sum_{i=0}^{n-1} n^2_i <= 2n
             send main table hash function parameters a, b, m, p (MPI_Bcast)
         For each slot i 
             Repeat
               choose an universal hash function among H_{p, n_i^2}
             Until there is not collision
   @endverbatim 
 *
 *   @note 
 *      - Original paper : Michael L. Fredman, Janos Komlos, and Endre Szemeredi. Storing a sparse table with O(1) worst case access time. Journal of the ACM, 31(3):538–544, 1984.
 *      - For an explenation to perfect hashing, see : Thomas H. Cormen, Charles E. Leiserson, Ronald L. Rivest, and Clifford Stein. 2009. Introduction to Algorithms, Third Edition (3rd ed.). p277-282.
 *      - For a summary of result, see : https://www.it.uu.se/edu/course/homepage/avalgo/vt11/PerfectHash.pdf
 */

/**
 *   @date 2018 Feb 26
 *   @brief extended modulo function, always return a positif integer
 *   @param[in] a an integer
 *   @param[in] b an integer
 */

static size_t mod(size_t a, size_t b) {
    size_t r = a % b;
    return r < 0 ? r + b : r;
}

/**
 *   @date 2018 Feb 26
 *   @brief verify if an integer is prime or not
 *   @param[in] n an integer
 *   @return 1 if n is prime, 0 otherwise
 */

static size_t isPrime (size_t n) {
    if (n < 2) { return 0; }
    if (n == 2) { return 1; }
    if (n % 2 == 0) { return 0; }
    if (n == 3) { return 1; }
    if (n % 3 == 0) { return 0; }
    for (size_t i = 3; i <= sqrt(n); i += 2) if (n % i == 0) { return 0; }
    return 1;
}

/**
 *   @date 2018 Feb 26
 *   @brief return the next prime greater than @p m
 *   @param[in] m an integer
 *   @return the next prime
 */

static size_t getPrimeGT(size_t m) {
    size_t p = m % 2 == 0 ? m + 1 : m + 2;

    while (!isPrime(p)) {
        p += 2;
    }

    return p;
}

/**
 *   @date 2018 Feb 26
 *   @brief convert a string to a 64 bits MD5.
 *   @param str the string to convert
 *   @return the associated 64bits MD5
 *   @note This function depend on openssl/md5.h
 *   @note It is possible that 64bits fingerprint is not sufficient when we
 *         scale, should modify it and all fingerprint occurrence.
 */

unsigned long long string2MD5(char *str) {
    unsigned char digest[MD5_DIGEST_LENGTH];
    MD5_CTX context;
    MD5_Init(&context);
    MD5_Update(&context, str, strlen(str));
    MD5_Final(digest, &context);
    size_t a = *(size_t *)digest;
    size_t b = *(size_t *)(digest + 8);
    // XOR them to get a shorter fingerprint 
    return a ^ b;
}

/**
 *   @date 2018 Feb 26
 *   @brief Compute the fingerprint of a read.
 *   @details This function take account the Qname and the order in a pair
 *            i.e. first in pair or second in pair. Hence, we can distinguish
 *            a read from its mate.
 *   @param read a read
 *   @return the associated fingerprint of @p read
 */

unsigned long long read2Fingerprint(readInfo *read) {
    // QNAME have at most 255 char
    char *buffer = calloc(300, sizeof(char));
    assert(strcpy(buffer, read->Qname));

    // if (read->valueFlag == 0)
    //         fprintf(stderr, "in READTOFINGERPRINT on  problem with read %s valueFlag is NULL \n", read->Qname);
             
    //int firstInPair = readBits(read->valueFlag, 6);
    //int secondInPair = readBits(read->valueFlag, 7);
    //assert(firstInPair || secondInPair);
    //char numInPair = digit2char(firstInPair + 2*secondInPair);

    assert ( (read->pair_num ==1) || (read->pair_num == 2));
    char numInPair = digit2char(read->pair_num);
    // string to hash is <QNAME><numberInPair>
    assert(numInPair);
    assert(strcat(buffer, &numInPair));
    size_t fingerprint = string2MD5(buffer);
    assert(fingerprint);
    free(buffer);
    return fingerprint;
}


/**
 *   @date 2018 Feb 26

 *   @brief A universal hash function \f[ h(k) = (ak + b \pmod p) \pmod m \f]
 *   @param hp hash parameter
 *   @param p prime number greater than @p hp.m
 *   @param k key i.e. fingerprint of a read
 *   @return hash value of @p k
 *   @note 
 *         - Original paper : Carter, Larry; Wegman, Mark N. (1979). "Universal Classes of Hash Functions". Journal of Computer and System Sciences. 18 (2): 143–154.
 */

int univHash(hashParam hp, int p, size_t k) {
    return mod(mod(hp.a * k + hp.b, p), hp.m);
}

/**
 *   @date 2018 Feb 26
 *   @brief Compute space using by the hash function
 *   @param hp hash parameter
 *   @param p prime number greater than @p hp.m
 *   @param arr an array of reads
 *   @param size the size of @p arr
 *   @return the amount of slots in the second level.
 */

int computeSpace(hashParam hp, int p, readInfo **arr, int size) {
    int sum = 0;
    int *count = calloc(hp.m, sizeof(int));

    for (int i = 0; i < size; i++) {
        if (arr[i]) {
            count[univHash(hp, p, arr[i]->fingerprint)]++;
        }
    }

    // for each slot in the main table,
    // we will allocate the square of elements which collide
    // in the second level table.
    for (int i = 0; i < hp.m; i++) {
        sum += count[i] * count[i];
    }

    free(count);
    return sum;
}

/**
 *   @date 2018 Feb 26
 *   @brief Given hash parameters, construct main universal hash function
 *   @param hp hash parameter
 *   @param p prime number greater than @p hp.m
 *   @param arr an array of reads
 *   @param size the size of @p arr
 */

void constructMainUnivHashWithHp(hashParam hp, int *p, readInfo **arr, int size) {

    *p = getPrimeGT(hp.m);

    int result = 0;
    // We search a prime number so that the hash table size is less than
    // twice of elements to store. In fact, E[\sum_{i=0}^{n-1} m_i] < 2n.
    while ((result = computeSpace(hp, *p, arr, size)) >= 2 * hp.m) {
        *p = getPrimeGT(*p);
    }

}

/**
 *   @date 2018 Feb 26
 *   @brief Construct main universal hash function
 *   @param hp hash parameter
 *   @param p prime number greater than @p hp.m
 *   @param arr an array of reads
 *   @param size the size of @p arr
 */

void constructMainUnivHash(hashParam *hp, int *p, readInfo **arr, int size) {

    // There is only one element in the second level hash table, 
    // we know that its index is 0, so choose a, b, p equal to 0.
    if (hp->m == 1) {
        hp->a = hp->b = *p = 0;

    // We choose a, b randomly and search a prime so that the
    // hash table size is reasonable.
    } else {
        hp->a = rand() + 1;
        hp->b = rand() + 1;
        constructMainUnivHashWithHp(*hp, p, arr, size);
    }

}

/**
 *   @date 2018 Feb 26
 *   @brief Test whether or not an universal hash function has collision.
 *   @param hp hash parameter
 *   @param p prime number greater than @p hp.m
 *   @param arr an array of reads
 *   @param size the size of @p arr
 *   @note This function is use to ensure that we are collision-free in second level hash table.
 */

int haveCollision(hashParam hp, int prime, size_t *arr, int size) {
    char *count = calloc(hp.m, sizeof(char));
    int flag = 0;

    for (int i = 0; i < size; i++) {
        assert(arr[i]);
        int h = univHash(hp, prime, arr[i]);

        if (count[h]) {
            flag = 1;
            break;
        }

        count[h] = 1;
    }

    free(count);
    return flag;
}

/**
 *   @date 2018 Feb 26
 *   @brief Construct second level hash table
 *   @param htbl a hash table
 *   @param arr an array of reads
 *   @param size the size of @p arr
 */

void constructSecTable(hashTable *htbl, readInfo **arr, int size) {
    // main hash table parameters
    hashParam hp = htbl->h;
    int prime = htbl->prime;

    int *count = calloc(hp.m, sizeof(int));
    int *index = calloc(hp.m, sizeof(int));
    size_t **subarr = malloc(sizeof(size_t *) * hp.m);
    htbl->secTable = malloc(sizeof(secTable *) * hp.m);

    // allocate all second level hash tables
    for (int i = 0; i < hp.m; i++) {
        htbl->secTable[i] = malloc(sizeof(secTable));
        assert(htbl->secTable[i]);
    }

    // count reads distribution on the main table
    for (int i = 0; i < size; i++) {
        if (arr[i]) {
            count[univHash(hp, htbl->prime, arr[i]->fingerprint)]++;
        }
    }

    // allocate each array of read's fingerprint which map to slot i
    for (int i = 0; i < hp.m; i++) {
        subarr[i] = malloc(sizeof(unsigned long long) * count[i]);
        assert(subarr[i]);
    }

    // fill these array with fingerprint following distribution 
    for (int i = 0; i < size; i++) {
        if (arr[i]) {
            int bucket = univHash(hp, htbl->prime, arr[i]->fingerprint);
            subarr[bucket][index[bucket]++] = arr[i]->fingerprint;
        }
    }

    /* compute hash param for each second table */
    for (int i = 0; i < hp.m; i++) {
        // second hash table parameters
        // allocate second level hash table with quadratic size
        hashParam sechp = {.a = rand() + 1, .b = rand() + 1, .m = count[i] * count[i]};
        htbl->secTable[i]->table = calloc(sechp.m, sizeof(readInfo *));

        // If there is only one element, we can choose directly a, b equal to 0
        if (sechp.m <= 1) {
            sechp.a = sechp.b = 0;

        } else {
            int collisionCounter = 0;

            // If we use an universal hash function to hash n keys into n^2 slots,
            // the probability of any collisions is less than 1/2.
            // We choose randomly a and b until there is not collision.
            while (haveCollision(sechp, htbl->prime, subarr[i], count[i])) {
                collisionCounter++;

                if (collisionCounter >= 1000) {
                   md_log_error("Cannot construct perfect hash table. Input data is probably not valid.\n");
                   exit(EXIT_FAILURE);
                }

                sechp.a = rand() + 1;
                sechp.b = rand() + 1;
            }
        }

        htbl->secTable[i]->h = sechp;

    }

    // insert read in their dedicated slot
    for (int i = 0; i < size; i++) {
        if (arr[i]) {
            hashTableInsert(htbl, arr[i]);
        }
    }

    for (int i = 0; i < hp.m; i++) {
        free(subarr[i]);
    }

    free(subarr);
    free(index);
    free(count);
}

/**
 *   @date 2018 Mar 23
 *   @brief Replace hashed read by another
 *   @param[out] htbl a hash table
 *   @param[in] read a read
 */

void hashTableInsert(hashTable *htbl, readInfo *read) {

    // compute which slot read maps to in the main table
    int h1 = univHash(htbl->h, htbl->prime, read->fingerprint);
    hashParam sechp = htbl->secTable[h1]->h;

    assert(sechp.m != 0) ;

    // compute index in the second level hash table 
    int h2 = univHash(sechp, htbl->prime, read->fingerprint);
    // assign the read in its place
    htbl->secTable[h1]->table[h2]  = read;

}

/**
 *   @date 2018 Feb 26
 *   @brief Given hash parameters, initialize a perfect hashing table.
 *   @param htbl a hash table
 *   @param hp universal hash parameters
 *   @param arr an array of reads
 *   @param size the size of @p arr
 */

void hashTableInitWithHp(hashTable *htbl, hashParam hp, readInfo **arr, int size) {
    int prime;

    constructMainUnivHashWithHp(hp, &prime, arr, size);

    htbl->h = hp;
    htbl->prime = prime;
    htbl->size = size;

    constructSecTable(htbl, arr, size);
}

/**
 *   @date 2018 Feb 26
 *   @brief Initialize a perfect hashing table.
 *   @param htbl a hash table
 *   @param arr an array of reads
 *   @param size the size of @p arr
 */
void hashTableInit(hashTable *htbl, readInfo **arr, int size) {
    int prime;
    hashParam hp = {.m = size, .a = 0, .b = 0};

    constructMainUnivHash(&hp, &prime, arr, size);

    htbl->h = hp;
    htbl->prime = prime;
    htbl->size = size;

    constructSecTable(htbl, arr, size);
}

/**
 *   @date 2018 Feb 26
 *   @brief Destroy perfect hashing table.
 *   @param htbl a hash table
 *   @note This function doesn't free reads, this task is left to destroyLBList().
 */

void hashTableDestroy(hashTable *htbl) {
    for (int i = 0; i < htbl->h.m; i++) {
        secTable *stbl = htbl->secTable[i];
        free(stbl->table);
        free(stbl);
    }

    free(htbl->secTable);
    free(htbl);
}

/**
 *   @date 2018 Feb 26
 *   @brief Print out perfect hashing table.
 *   @param htbl a hash table
 *   @param arr array of reads
 *   @param size size of @p arr
 */

void printPerfectHashTable(hashTable *htbl) {
    printf("Main universal hash function parameters :\na=%d, b=%d, p=%d, m=%llu\n", htbl->h.a, htbl->h.b, htbl->prime, htbl->h.m);
    printf("allocated space = %d\n", htbl->size);

    for (int i = 0; i < htbl->h.m; i++) {
        hashParam sechp = htbl->secTable[i]->h;
        printf("%-3d : a=%-10d, b=%-10d, p=%-4d, m=%-2llu :: ", i, sechp.a, sechp.b, htbl->prime, sechp.m);

        for (int j = 0; j < sechp.m; j++) {
            printf("(%d-%s-%llu)-", j, htbl->secTable[i]->table[j] == NULL ? "X" : htbl->secTable[i]->table[j]->Qname, htbl->secTable[i]->table[j] == NULL ? 0 : htbl->secTable[i]->table[j]->fingerprint);
        }

        puts("");
    }
}

/**
 *   @date 2018 Feb 26
 *   @brief Given a fingerprint return the associated read.
 *   @param htbl a hash table
 *   @param read a read fingerprint
 *   @param if found, the associated read, NULL otherwise.
 */

readInfo *getReadFromFingerprint(hashTable *htbl, size_t fingerprint) {
    
    // compute which slot read maps to in the main table
    int h1 = univHash(htbl->h, htbl->prime, fingerprint);
    hashParam sechp = htbl->secTable[h1]->h;

    // we don't found the read we search
    if (sechp.m == 0) return NULL;
    // compute index in the second level hash table 
    int h2 = univHash(sechp, htbl->prime, fingerprint);
    readInfo *read = htbl->secTable[h1]->table[h2];

    // we don't found the read we search
    if (!read || read->fingerprint != fingerprint) {
	return NULL;

    } else {
        return read;
    }
}

/**
 *   @date 2018 Feb 26
 *   @brief Given a Qname and pair order return the associated read.
 *   @param htbl a hash table
 *   @param Qname the Qname of a read
 *   @param PairNum the order in a pair of reads (i.e. first or second)
 *   @param if found, the associated read, NULL otherwise.
 */

readInfo *getReadFromQnameAndPairNum(hashTable *htbl, char *Qname, int PairNum) {

    //fprintf(stderr, "in PerfectHash [getReadFromQnameAndPairNum] problem with read %s of flag %u \n ", read->Qname);
    assert(Qname);
    assert(PairNum == 1 || PairNum == 2);
    /* first in pair = 64, second in pair = 128 */
    //int flag = PairNum * 64;
    //readInfo read = {.Qname = Qname, .valueFlag = flag};
    readInfo read = {.Qname = Qname, .pair_num = PairNum};
    size_t fingerprint = read2Fingerprint(&read);
    return getReadFromFingerprint(htbl, fingerprint);
}

/**
 *   @date 2018 Mar 23
 *   @brief Given a read return its mate fingerprint
 *   @param[in] read a read
 *   @return @p read's mate fingerprint
 */

unsigned long long read2mateFP(readInfo *read) {


    //unsigned int firstInPair  = readBits(read->valueFlag, 6);
    //unsigned int secondInPair = readBits(read->valueFlag, 7);
    //size_t matePairNum  = 2 * firstInPair + secondInPair;
    size_t matePairNum = 0; 

    assert( (read->pair_num == 1) || (read->pair_num == 2));
    
    if (read->pair_num == 1 ) matePairNum = 2;
    if (read->pair_num == 2 ) matePairNum = 1;

    if (!(matePairNum == 1 || matePairNum == 2))
        fprintf(stderr, "in PerfectHash [read2mateFP] problem with read %s of flag %u \n ", read->Qname, read->valueFlag);

    //matePairNum = 1;
    //assert(matePairNum == 1 || matePairNum == 2);
    size_t mateFlag = matePairNum; // * 64;
    //readInfo mate = {.Qname = read->Qname, .valueFlag = mateFlag};
    readInfo mate = {.Qname = read->Qname, .pair_num = mateFlag};
    return read2Fingerprint(&mate);
}

/**
 *   @date 2018 Feb 26
 *   @brief Given a read return its mate.
 *   @param htbl a hash table
 *   @param read a read
 *   @param if found, the associated mate, NULL otherwise.
 */

readInfo *getMateFromRead(hashTable *htbl, readInfo *read) {
    assert( (read->pair_num == 1) || (read->pair_num == 2));
    unsigned long long mateFP = read2mateFP(read);
    return getReadFromFingerprint(htbl, mateFP);
}

/**
 *   @date 2018 Apr 7
 *   @brief Create MPI derived type corresponding to hashParam
 *   @param[out] HPType the MPI derived type
 */

static void createHPType(MPI_Datatype *HPType) {

    int blocks[3] = {1, 1, 1};
    MPI_Datatype types[3] = {MPI_INT, MPI_INT, MPI_UNSIGNED_LONG_LONG};
    MPI_Aint displacements[3];

    MPI_Aint intex, ullex, lb;
    MPI_Type_get_extent(MPI_INT, &lb, &intex);
    MPI_Type_get_extent(MPI_UNSIGNED_LONG_LONG, &lb, &ullex);

    displacements[0] = 0;
    displacements[1] = intex;
    displacements[2] = intex * 2;

    MPI_Type_create_struct(3, blocks, displacements, types, HPType);
    MPI_Type_commit(HPType);
}

/**
 *   @date 2018 Feb 26
 *   @brief Construct perfect hashing table through all process.
 *   @param htbl a hash table
 *   @param arr an array of reads
 *   @param size the size of @p arr
 */

void shareHpAndConstructHtbl(hashTable *htbl, readInfo **arr, int size, MPI_Comm comm) {

    int rank, numprocs;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &numprocs);
    hashParam  hp ;
    int maxSize = 0;
    MPI_Allreduce(&size, &maxSize, 1, MPI_INT, MPI_MAX, comm);

    int actualSize = 0;

    for (int i = 0; i < size; i++) {
        if (arr[i]) {
            actualSize++;
        }
    }

    if (rank == 0) {
        hp.a = rand() + 1;
        hp.b = rand() + 1;
        hp.m = maxSize;

        hashTableInitWithHp(htbl, hp, arr, size);
    }

    /* transfert hash parameter */
    
    MPI_Datatype hp_type;
    createHPType(&hp_type);

    MPI_Bcast(&hp, 1, hp_type, 0, comm);

    if (rank != 0) {
        hashTableInitWithHp(htbl, hp, arr, size);
    }

    MPI_Type_free(&hp_type);
}



