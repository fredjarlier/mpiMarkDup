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
    Xingwei SANG,       Paris Descartes
*/

#include <stdarg.h>
#include <stddef.h>
#include <setjmp.h>
#include <stdlib.h>
#include <string.h>
#include <cmocka.h>

#include <math.h>
#include "mark_duplicates.h"
#include "mpiMD_utils.h"
#include "log.h"
#include "reads.h"
#include "perfectHash.h"

/* date 2018/03/30
 * Test functions : updateReadPosInClusterSortedByPCQ, insertAndSortClusterReadByPhred, getFirstInCigarSortedByPos
 */

#define SIZE 100
#define PART1 0.3
#define PART2 0.2
#define PART3 0.4
#define PART4 0.1

typedef struct {
    List *list;
    int nbReads1;
    int nbReads2;
    int nbReads3;
    int nbReads4;
} test_t;

static readInfo *initRead(readInfo *read, int pos, int posInCluster, size_t quality, int i) {
    read = malloc(sizeof(readInfo));
    char *cigar = malloc(sizeof(char) * 25);
    sprintf(cigar, "%dM", i);
    read->Qname = strdup("qname");
    read->lb_name = strdup("LB");
    read->cigar = cigar;
    read->coordPos = pos;
    read->pos_in_cluster = posInCluster;
    read->phred_score = quality;
    return read;
}

static int List_setup(void **state) {
    readInfo **reads = malloc(sizeof(readInfo *)*SIZE);
    int pos1 = 1001;
    int nbReads1 = SIZE * PART1;
    int pos2 = 1002;
    int nbReads2 = SIZE * PART2;
    int pos3 = 1003;
    int nbReads3 = SIZE * PART3;
    int pos4 = 10004;
    int nbReads4 = SIZE * PART4;
    List *list = allocate_new_list("LB");

    // see printt to look how  was contructed reads
    for (int i = 0; i < nbReads1; i++) {
        if (i < (floor(PART1 * nbReads1))) {
            reads[i] = initRead(reads[i], pos1, i + 1, SIZE - (i + 1), 8);

        } else if (i < (floor((PART1 + PART2)*nbReads1))) {
            reads[i] = initRead(reads[i], pos1, i + 1 - floor(PART1 * nbReads1), SIZE - (i + 1), 9);

        } else if (i < floor((PART1 + PART2 + PART3)*nbReads1)) {
            reads[i] = initRead(reads[i], pos1, i + 1 - floor((PART1 + PART2) * nbReads1), SIZE - (i + 1), 10);

        } else {
            reads[i] = initRead(reads[i], pos1, i + 1 - floor((PART1 + PART2 + PART3) * nbReads1), SIZE - (i + 1), 11);
        }

        //printf("number : %d\tpos : %d\t cigar : %s\t pos_in_cluster: %d\n",i,reads[i]->coordPos,reads[i]->cigar,reads[i]->pos_in_cluster);
        listAppendEnd(list, reads[i]);
    }

    for (int i = 0; i < nbReads2; i++) {
        if (i < (floor(PART1 * nbReads2))) {
            reads[nbReads1 + i] = initRead(reads[nbReads1 + i], pos2, i + 1, SIZE - (i + 1), 8);

        } else if (i < (floor((PART1 + PART2)*nbReads2))) {
            reads[nbReads1 + i] = initRead(reads[nbReads1 + i], pos2, i + 1 - floor(PART1 * nbReads2), SIZE - (i + 1), 9);

        } else if (i < floor((PART1 + PART2 + PART3)*nbReads2)) {
            reads[nbReads1 + i] = initRead(reads[nbReads1 + i], pos2, i + 1 - floor((PART1 + PART2) * nbReads2), SIZE - (i + 1), 10);

        } else {
            reads[nbReads1 + i] = initRead(reads[nbReads1 + i], pos2, i + 1 - floor((PART1 + PART2 + PART3) * nbReads2), SIZE - (i + 1), 11);
        }

        // printf("number : %d\tpos : %d\t cigar : %s\t pos_in_cluster: %d\n",nbReads1+i,reads[nbReads1+i]->coordPos,reads[nbReads1+i]->cigar,reads[nbReads1+i]->pos_in_cluster);

        listAppendEnd(list, reads[nbReads1 + i]);
    }

    for (int i = 0; i < nbReads3; i++) {
        if (i < (floor(PART1 * nbReads3))) {
            reads[nbReads1 + nbReads2 + i] = initRead(reads[nbReads1 + nbReads2 + i], pos3, i + 1, SIZE - (i + 1), 8);

        } else if (i <= (floor((PART1 + PART2)*nbReads3))) {
            reads[nbReads1 + nbReads2 + i] = initRead(reads[nbReads1 + nbReads2 + i], pos3, i + 1 - floor(PART1 * nbReads3), SIZE - (i + 1), 9);

        } else if (i <= floor((PART1 + PART2 + PART3)*nbReads3)) {
            reads[nbReads1 + nbReads2 + i] = initRead(reads[nbReads1 + nbReads2 + i], pos3, i + 1 - floor((PART1 + PART2) * nbReads3), SIZE - (i + 1), 10);

        } else {
            reads[nbReads1 + nbReads2 + i] = initRead(reads[nbReads1 + nbReads2 + i], pos3, i - floor((PART1 + PART2 + PART3) * nbReads3), SIZE - (i + 1), 11);
        }

        //printf("number : %d\tpos : %d\t cigar : %s\t pos_in_cluster: %d\n",nbReads1+nbReads2+i,reads[nbReads2+nbReads1+i]->coordPos,reads[nbReads1+nbReads2+i]->cigar,reads[nbReads1+nbReads2+i]->pos_in_cluster);
        listAppendEnd(list, reads[nbReads1 + nbReads2 + i]);
    }

    for (int i = 0; i < nbReads4; i++) {
        if (i < (floor(PART1 * nbReads4))) {
            reads[nbReads1 + nbReads2 + nbReads3 + i] = initRead(reads[nbReads1 + nbReads2 + nbReads3 + i], pos4, i + 1, SIZE - (i + 1), 8);

        } else if (i < (floor((PART1 + PART2)*nbReads4))) {
            reads[nbReads1 + nbReads2 + nbReads3 + i] = initRead(reads[nbReads1 + nbReads2 + nbReads3 + i], pos4, i + 1 - floor(PART1 * nbReads4), SIZE - (i + 1), 9);

        } else if (i < floor((PART1 + PART2 + PART3)*nbReads4)) {
            reads[nbReads1 + nbReads2 + nbReads3 + i] = initRead(reads[nbReads1 + nbReads2 + nbReads3 + i], pos4, i + 1 - floor((PART1 + PART2) * nbReads4), SIZE - (i + 1), 10);

        } else {
            reads[nbReads1 + nbReads2 + nbReads3 + i] = initRead(reads[nbReads1 + nbReads2 + nbReads3 + i], pos4, i + 1 - floor((PART1 + PART2 + PART3) * nbReads4), SIZE - (i + 1), 11);
        }

        //printf("number : %d\tpos : %d\t cigar : %s\tpos_in_cluster : %d\n",nbReads1+nbReads2+nbReads3+i,reads[nbReads2+nbReads3+nbReads1+i]->coordPos,reads[nbReads1+nbReads2+nbReads3+i]->cigar,reads[nbReads1+nbReads2+nbReads3+i]->pos_in_cluster);
        listAppendEnd(list, reads[nbReads1 + nbReads2 + nbReads3 + i]);
    }

    test_t *test = malloc(sizeof(test_t));
    test->list = list;
    test->nbReads1 = nbReads1;
    test->nbReads2 = nbReads2;
    test->nbReads3 = nbReads3;
    test->nbReads4 = nbReads4;
    *state = test;
    return 0;
}

static int List_teardown(void **state) {
    test_t *test = *state;
    deleteListReads(test->list);
    return 0;
}

static void test_getFirstInClusterSortedByPCQ(void **state) {
    test_t *test = *state;
    List *list = test->list;
    ListNode *node = malloc(sizeof(ListNode));
    // initialise read : cigar =10M,coordPos=1004,quality=SIZE+1
    readInfo *read = initRead(read, 1004, 0, SIZE + 1, 10);
    node->read = read;
    
    //test the head of cluster at coordPos 1004
    ListNode *nodePos = list->tail;
    nodePos = getFirstReadInCluster(list, nodePos);
    ListNode *firstInCluster = list->head;
    int indice = test->nbReads1 + test->nbReads2 + test->nbReads3;

    for (int i = 0; i < indice; i++) {
        firstInCluster = firstInCluster->next;
    }

    assert_ptr_equal(nodePos, firstInCluster);
    
    //test the head of cluster at coordPos 1004 and cigar 10M
    nodePos = list->tail;

    for (int i = 0; i < floor((PART1 + PART2)*test->nbReads4); i++) {
        firstInCluster = firstInCluster->next;
    }

    assert_ptr_equal(getFirstInCigarSortedByPos(list, node, nodePos), firstInCluster);
    free(read);
    
    // initialise read : cigar =11M,coordPos=1002,quality=SIZE+1
    readInfo *read1 = initRead(read1, 1002, 0, SIZE + 1, 11);
    node->read = read1;
    indice = test->nbReads1;
    ListNode *tempNode = list->head;

    
    //test the head of cluster at coordPos 1002
    for (int i = 0; i < indice; i++) {
        tempNode = tempNode->next;
    }

    firstInCluster = tempNode;
    nodePos = getFirstReadInCluster(list, tempNode);
    assert_ptr_equal(nodePos, firstInCluster);

    
    //test the head of cluster at coordPos 1002 and cigar 11M
    for (int i = 0; i < floor((PART1 + PART2 + PART3)*test->  nbReads2); i++) {
        firstInCluster = firstInCluster->next;
    }

    assert_ptr_equal(getFirstInCigarSortedByPos(list, node, nodePos), firstInCluster);

}

static void test_insertAndSortClusterReadByPhred(void **state){
    test_t *test=*state;
    List *list = test->list;
    //initalise read : cigar =10M, coordPos=1004, quality = SIZE-1
    readInfo *read=initRead(read,1004,0,SIZE-1,10);
    insertAndSortClusterReadByPhred(list,read);
    ListNode *tempNode=list->head;
    int indice = test->nbReads1+test->nbReads3+test->nbReads2;
    for(int i =0;i<indice;i++)
        tempNode=tempNode->next;
    indice =floor((PART1+PART2)*test->nbReads4);
    for(int i = 0; i<indice;i++){
        tempNode=tempNode->next;
    }
    assert_ptr_equal(read,tempNode->read);
    
    //initialise read : cigar=8M, coordPos=1004, quality= SIZE
    readInfo *read2=initRead(read2,1004,0,SIZE,8);
    insertAndSortClusterReadByPhred(list,read2);
    tempNode=list->head;
    indice = test->nbReads1+test->nbReads3+test->nbReads2;
    for(int i =0;i<indice;i++)
        tempNode=tempNode->next;

    assert_ptr_equal(read2,tempNode->read);

}

int main(void) {
    const struct CMUnitTest List[] = {
        cmocka_unit_test_setup_teardown(test_getFirstInClusterSortedByPCQ, List_setup, List_teardown),
        cmocka_unit_test_setup_teardown(test_insertAndSortClusterReadByPhred,List_setup,List_teardown),
    };
    cmocka_run_group_tests(List, NULL, NULL);

    return 0;
}

