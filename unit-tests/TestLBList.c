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
     TestLBList.c
     
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

#include "reads.h"
#include "readList.h"
#include "mpiMD_utils.h"

/* @date 2018/03/30
 * Test structure LBList
 * Functions tested : destroyLBList, initLBList, countExternalMateInList
 */
#define SIZE 4
#define NB_READS 5

static void initRead(readInfo *read, char *LBname) {
    static int pos = 1;
    read->Qname = malloc(sizeof(char) * 100);
    read->cigar = malloc(sizeof(char) * 100);
    strcpy(read->Qname, "Qname");
    strcpy(read->cigar, "Cigar");
    read->lb_name = strdup(LBname);
    read->pos_in_cluster = pos;
    pos ++;

    if ( pos > 5) {
        pos = 1;
    }
}

static int LBList_setup(void **state) {
    LBList *lbList = malloc(sizeof(LBList));
    lbList->lb_num = 6 ;
    char **library = malloc(sizeof (char *)*SIZE);

    for ( int i = 0; i < SIZE; i++) {
        library[i] = malloc(sizeof(char) * 100);
    }

    strcpy(library[0], "LB0");
    strcpy(library[1], "LB1");
    strcpy(library[2], "LB2");
    strcpy(library[3], "unmapped");
    initLBList(lbList, SIZE, library);

    for (int i = 0; i < SIZE; i++) {
        for (int j = 0   ; j < NB_READS; j++) {
            readInfo *read = malloc(sizeof(readInfo));;
            initRead(read, library[i]);
            listAppendEnd(lbList->lb_list[i], read);
        }
    }

    *state = lbList;
    return 0;
}

static int LBList_teardown(void **state) {
    LBList *lbList = *state;
    destroyLBList(lbList);
    return 0;
}

static void test_FindLbIndex(void **state) {
    LBList *lbList = *state;
    List **plist = lbList->lb_list;
    assert_int_equal(FindLBIndex("LB1", plist, SIZE), 1);
    assert_int_equal(FindLBIndex("LB2", plist, SIZE), 2);
    assert_int_equal(FindLBIndex("LB0", plist, SIZE), 0);
    assert_int_equal(FindLBIndex("unmapped", plist, SIZE), 3);
}

static void test_sizes(void **state) {
    LBList *lbList = *state;
    List **plist = lbList->lb_list;

    for (int i = 0 ; i < SIZE; i++) {
        assert_int_equal(plist[i]->size, NB_READS);
    }

    assert_int_equal(lbList->lb_num, SIZE);
}

static void test_countExternalMateInList(void **state) {
    LBList *lbList = *state;
    size_t *table = malloc(sizeof(size_t) * SIZE);
    size_t table2[SIZE] = {2, 4, 1, 0};
    List **plist = lbList->lb_list;

    for ( int i = 0; i < SIZE; i++) {
        ListNode *node1 = plist[i]->head;
        ListNode *node2 = node1->next;
        ListNode *node3 = node2->next;
        ListNode *node4 = node3->next;
        ListNode *node5 = plist[i]->tail;

        if ( i == 0) {
            node1->read->check_with_bruck = 1;
            node2->read->check_with_bruck = 0;
            node3->read->check_with_bruck = 0;
            node4->read->check_with_bruck = 0;
            node5->read->check_with_bruck = 1;
        }

        if ( i == 1) {
            node1->read->check_with_bruck = 1;
            node2->read->check_with_bruck = 0;
            node3->read->check_with_bruck = 1;
            node4->read->check_with_bruck = 1;
            node5->read->check_with_bruck = 1;
        }

        if ( i == 2) {
            node1->read->check_with_bruck = 0;
            node2->read->check_with_bruck = 0;
            node3->read->check_with_bruck = 1;
            node4->read->check_with_bruck = 0;
            node5->read->check_with_bruck = 0;
        }

        if ( i == 3) {
            node1->read->check_with_bruck = 0;
            node2->read->check_with_bruck = 0;
            node3->read->check_with_bruck = 0;
            node4->read->check_with_bruck = 0;
            node5->read->check_with_bruck = 0;
        }

    }

    countExternalMateInList(lbList, table);

    for (int i = 0; i < SIZE; i++) {
        assert_int_equal(table[i], table2[i]);
    }
}

int main(void) {
    const struct CMUnitTest test_LBList[] = {
        cmocka_unit_test_setup_teardown(test_FindLbIndex, LBList_setup, LBList_teardown),
        cmocka_unit_test_setup_teardown(test_sizes, LBList_setup, LBList_teardown),
        cmocka_unit_test_setup_teardown(test_countExternalMateInList, LBList_setup, LBList_teardown),
    };
    cmocka_run_group_tests(test_LBList, NULL, NULL);
    return 0;
}
