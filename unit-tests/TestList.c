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
     TestList.c
     
   Authors:
    Frederic Jarlier,   Institut Curie
    Firmin Martin,      Paris Descartes
    Xingwei SANG,       Paris Descartes
*/

#include <stdarg.h>
#include <stddef.h>
#include <setjmp.h>
#include <stdlib.h>
#include <string.h>
#include <cmocka.h>

#include "readList.h"

#define N_NODES 5

/**
 * @date 2018/02/15
 * Test structure List; ListNode, LBList.
 * Functions tested : insertReadAfter, insertReadBefore ,listAppendEnd, removeRead.
 */

typedef struct {
    List **List;
    readInfo **reads;
    ListNode **Nodes;
} test_t;

static void initRead(readInfo **read) {
    (*read) = malloc(sizeof(readInfo));
    (*read)->Qname = malloc(sizeof(char) * 256);
    (*read)->lb_name = malloc(sizeof(char) * 256);
    (*read)->cigar = malloc(sizeof(char) * 256);
    strcpy((*read)->Qname, "Qname");
    strcpy((*read)->lb_name, "lb_name");
    strcpy((*read)->cigar, "cigar");
    (*read)->coordPos =1000;
}

static int List_setup(void **state) {
    test_t *test = malloc(sizeof(test_t));
    test->reads = malloc(sizeof(readInfo *) * N_NODES);
    test->Nodes = malloc(sizeof(ListNode *) * N_NODES);
    test->List = malloc(sizeof(List *));

    for (int i = 0; i < N_NODES; i++) {
        test->Nodes[i] = malloc(sizeof(ListNode));
        initRead(&test->reads[i]);
        test->Nodes[i]->read = test->reads[i];
    }

    *state = test;
    return 0;
}

static int List_teardown(void **state) {
    test_t *test = *state;
    List *List = *test->List;
    readInfo **reads = test->reads;
    ListNode **Nodes = test->Nodes;

    deleteListReads(List);
    free(reads);
    free(Nodes);
    free(test->List);
    free(test);
    return 0;
}

/* test insertReadAfter, insertReadBefore */
static void ListPtr_InsertRead(void **state) {
    test_t *test = *state;
    List **pList = test->List;
    readInfo **reads = test->reads;
    ListNode **Nodes = test->Nodes;
    *pList = allocate_new_list("LB1");

    List* List = *pList;
    assert_non_null(List);
    assert_string_equal(List->LB, "LB1");
    assert_null(List->head);
    assert_null(List->tail);
    /* insert first read */
    insertReadAfter(List, Nodes[0], NULL);
    assert_ptr_equal(List->head, Nodes[0]);
    assert_ptr_equal(List->tail, Nodes[0]);
    assert_ptr_equal(List->head->read, reads[0]);
    assert_ptr_equal(List->tail->read, reads[0]);
    assert_null(List->tail->next);
    assert_null(List->head->prev);
    /* insert read in tail */
    insertReadAfter(List, Nodes[1], List->head);
    assert_ptr_equal(List->head, Nodes[0]);
    assert_ptr_equal(List->head->next, Nodes[1]);
    assert_ptr_equal(List->tail, Nodes[1]);
    assert_ptr_equal(List->tail->prev, Nodes[0]);
    assert_ptr_equal(List->head->read, reads[0]);
    assert_ptr_equal(List->tail->read, reads[1]);
    assert_null(List->tail->next);
    assert_null(List->head->prev);
    /* insert read in head */
    insertReadBefore(List, Nodes[2], List->head);
    assert_ptr_equal(List->head, Nodes[2]);
    assert_ptr_equal(List->tail, Nodes[1]);
    assert_ptr_equal(List->head->read, reads[2]);
    assert_ptr_equal(List->tail->read, reads[1]);
    assert_null(List->tail->next);
    assert_null(List->head->prev);
    /* insert read in middle */
    insertReadBefore(List, Nodes[3], Nodes[1]);
    assert_ptr_equal(List->head, Nodes[2]);
    assert_ptr_equal(List->tail, Nodes[1]);
    assert_ptr_equal(List->head->read, reads[2]);
    assert_ptr_equal(List->tail->read, reads[1]);
    assert_null(List->tail->next);
    assert_null(List->head->prev);
    /* test listAppendEnd*/
    listAppendEnd(List,reads[4]);
    assert_ptr_equal(List->tail->read,reads[4]);
    assert_null(List->tail->next);
    assert_int_equal(List->size,5);
    
}

//static void test_removeReadAndPosIncluster(void **state){
//    test_t *test = *state;
//    List **pList = test->List;
//    List *list = *pList;
//    *pList = allocate_new_list("LB");
//    ListNode **Nodes = test->Nodes;
//    insertReadAfter(list,Nodes[0],NULL);
//    assert_ptr_equal(list->head,Nodes[0]);
//    insertReadAfter(list,Nodes[1],Nodes[0]);
//    assert_ptr_equal(list->tail,Nodes[2]);
//    insertReadAfter(list,Nodes[2],Nodes[1]);
//    assert_ptr_equal(list->tail,Nodes[2]);
//    Nodes[0]->read->pos_in_cluster=1;
//    Nodes[1]->read->pos_in_cluster=2;
//    Nodes[2]->read->pos_in_cluster=3;
//    Nodes[3]->read->pos_in_cluster=4;
//    Nodes[4]->read->coordPos = 1001;
//    Nodes[4]->read->pos_in_cluster = 1;
//    insertReadAfter(list,Nodes[4],Nodes[3]);
//    /* test removeRead*/
//    assert_int_equal(removeRead(list,Nodes[0]),0);
//    //Does not free
//    // assert_null(Nodes[0]);
//    assert_ptr_equal(list->head,Nodes[1]);
//    assert_ptr_equal(list->head->next,Nodes[2]);
//    assert_int_equal(list->size,4);
//    assert_int_equal(Nodes[1]->read->pos_in_cluster,1);
//    assert_int_equal(Nodes[2]->read->pos_in_cluster,2);
//    assert_int_equal(removeRead(list,Nodes[2]),0);
//    assert_ptr_equal(Nodes[3],list->head->next);
//    //assert_null(Nodes[2]);
//    assert_int_equal(list->size,3);
//    assert_int_equal(Nodes[3]->read->pos_in_cluster,2);
//    assert_int_equal(Nodes[4]->read->pos_in_cluster,1);
//}

int main(void) {

    const struct CMUnitTest ListSize[] = {
       //cmocka_unit_test_setup_teardown(ListSize_InsertRead, List_setup, List_teardown),
       //cmocka_unit_test_setup_teardown(ListSize_DeleteRead, List_setup, List_teardown),
   };
    const struct CMUnitTest ListPtr[] = {
        cmocka_unit_test_setup_teardown(ListPtr_InsertRead, List_setup, List_teardown),
        //cmocka_unit_test_setup_teardown(test_removeReadAndPosIncluster, List_setup, List_teardown),
    };
    cmocka_run_group_tests(ListPtr, NULL, NULL);
    cmocka_run_group_tests(ListSize, NULL, NULL);
    return 0;
}
