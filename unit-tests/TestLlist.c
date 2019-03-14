#include <stdarg.h>
#include <stddef.h>
#include <setjmp.h>
#include <stdlib.h>
#include <string.h>
#include <cmocka.h>

#include  <stdio.h>

#include "llist.h"
#include "reads.h"

/**
 * @date 2018 Apr 12
 * @author Firmin Martin
 */

typedef struct {
    llist_t *list;
    readInfo **reads;
} test_t;

static void initRead(readInfo **read, int n) {
    (*read) = malloc(sizeof(readInfo));
    (*read)->readLb = 0;
    asprintf(&(*read)->Qname, "Qname-%d", n);
    (*read)->cigar   = strdup("cigar");
}

static int llist_setup(void **state) {
    test_t *test = malloc(sizeof(test_t));
    test->list = llist_create();
    test->reads = malloc(sizeof(readInfo *) * 10);

    for (int i = 0; i < 10; i++) {
        initRead(&test->reads[i], i);
    }

    *state = test;
    return 0;
}

static int llist_teardown(void **state) {
    test_t *test = *state;
    llist_t *list = test->list;
    readInfo **reads = test->reads;

    for (int i = 0; i < 10; i++) {
        freeRead(reads[i]);
    }

    llist_clear(list);
    llist_destruct(list);
    free(reads);
    free(test);
    return 0;
}

/* test for list->nil, list->head, list->tail */
static void llistPtr_InsertRead(void **state) {
    test_t *test = *state;
    llist_t *list = test->list;
    readInfo **reads = test->reads;
    assert_non_null(list->nil);
    assert_ptr_equal(list->head, list->nil);
    assert_ptr_equal(list->tail, list->nil);
    /* insert first read */
    llist_insert_index(list, 0, reads[0]);
    assert_ptr_equal(list->head->read, reads[0]);
    assert_ptr_equal(list->tail->read, reads[0]);
    /* insert read in tail */
    llist_insert_index(list, 1, reads[1]);
    assert_ptr_equal(list->head->read, reads[0]);
    assert_ptr_equal(list->tail->read, reads[1]);
    /* insert read in head */
    llist_insert_index(list, 0, reads[2]);
    assert_ptr_equal(list->head->read, reads[2]);
    assert_ptr_equal(list->tail->read, reads[1]);
    /* insert read in middle */
    llist_insert_index(list, 1, reads[3]);
    assert_ptr_equal(list->head->read, reads[2]);
    assert_ptr_equal(list->tail->read, reads[1]);
    //llist_clear(list);
    //assert_ptr_equal(list->head, list->nil);
    //assert_ptr_equal(list->tail, list->nil);
}

static void llistPtr_AppendRead(void **state) {
    test_t *test = *state;
    llist_t *list = test->list;
    readInfo **reads = test->reads;
    assert_non_null(list->nil);
    assert_ptr_equal(list->head, list->nil);
    assert_ptr_equal(list->tail, list->nil);
    /* first node */
    llist_append(list, reads[0]);
    assert_ptr_equal(list->head->read, reads[0]);
    assert_ptr_equal(list->tail->read, reads[0]);
    assert_ptr_equal(list->head->prev, list->nil);
    assert_ptr_equal(list->tail->next, list->nil);
    assert_ptr_equal(list->head, list->nil->next);
    assert_ptr_equal(list->tail, list->nil->prev);
    /* second node */
    llist_append(list, reads[1]);
    assert_ptr_equal(list->head->read, reads[0]);
    assert_ptr_equal(list->tail->read, reads[1]);
    assert_ptr_equal(list->head->prev, list->nil);
    assert_ptr_equal(list->tail->next, list->nil);
    assert_ptr_equal(list->head, list->nil->next);
    assert_ptr_equal(list->tail, list->nil->prev);
    /* third node */
    llist_append(list,reads[2]);
    assert_ptr_equal(list->head->read, reads[0]);
    assert_ptr_equal(list->tail->read, reads[2]);
    assert_ptr_equal(list->head->prev, list->nil);
    assert_ptr_equal(list->tail->next, list->nil);
    assert_ptr_equal(list->head->next->next->next, list->nil);
    assert_ptr_equal(list->tail->prev->prev->prev, list->nil);
    assert_ptr_equal(list->head, list->nil->next);
    assert_ptr_equal(list->tail, list->nil->prev);
    /* clear */
    //llist_readInfo_print(list, MPI_COMM_WORLD);
    llist_clear(list);
    assert_ptr_equal(list->head, list->nil);
    assert_ptr_equal(list->tail, list->nil);
}

static void llistPtr_DeleteRead(void **state) {
    test_t *test = *state;
    llist_t *list = test->list;
    readInfo* read;
    readInfo **reads = test->reads;
    /* delete in empty list */
    read = llist_delete(list, 0);
    assert_null(read);
    /* insert 4 reads */
    llist_insert_index(list, 0, reads[0]);
    llist_insert_index(list, 1, reads[1]);
    llist_insert_index(list, 2, reads[2]);
    llist_insert_index(list, 3, reads[3]);
    assert_ptr_equal(list->head->read, reads[0]);
    assert_ptr_equal(list->tail->read, reads[3]);
    read = llist_delete(list, 3);
    assert_ptr_equal(read, reads[3]);
    assert_ptr_equal(list->head->read, reads[0]);
    assert_ptr_equal(list->tail->read, reads[2]);
    read = llist_delete(list, 0);
    assert_ptr_equal(read, reads[0]);
    assert_ptr_equal(list->head->read, reads[1]);
    assert_ptr_equal(list->tail->read, reads[2]);
    read = llist_delete(list, 1);
    assert_ptr_equal(read, reads[2]);
    assert_ptr_equal(list->head->read, reads[1]);
    assert_ptr_equal(list->tail->read, reads[1]);
    read = llist_delete(list, 0);
    assert_ptr_equal(read, reads[1]);
    assert_ptr_equal(list->head, list->nil);
    assert_ptr_equal(list->tail, list->nil);
}

/* Test for list->size */

static void llistSize_InsertRead(void **state) {
    test_t *test = *state;
    llist_t *list = test->list;
    readInfo **reads = test->reads;
    llist_insert_index(list, 0, reads[0]);
    assert_int_equal(list->size, 1);
    llist_insert_index(list, 0, reads[1]);
    assert_int_equal(list->size, 2);
    llist_insert_index(list, 2, reads[2]);
    assert_int_equal(list->size, 3);
    llist_insert_index(list, 3, reads[3]);
    assert_int_equal(list->size, 4);
}

static void llistSize_DeleteRead(void **state) {
    test_t *test = *state;
    llist_t *list = test->list;
    readInfo **reads = test->reads;
    llist_delete(list, 0);
    assert_int_equal(list->size, 0);
    llist_insert_index(list, 0, reads[1]);
    llist_insert_index(list, 1, reads[2]);
    llist_insert_index(list, 2, reads[3]);
    llist_delete(list, 2);
    assert_int_equal(list->size, 2);
    llist_delete(list, 1);
    assert_int_equal(list->size, 1);
    llist_delete(list, 0);
    assert_int_equal(list->size, 0);
}



int main(void) {

    const struct CMUnitTest llistSize[] = {
        cmocka_unit_test_setup_teardown(llistSize_InsertRead, llist_setup, llist_teardown),
        cmocka_unit_test_setup_teardown(llistSize_DeleteRead, llist_setup, llist_teardown),
    };
    const struct CMUnitTest llistPtr[] = {
        cmocka_unit_test_setup_teardown(llistPtr_InsertRead, llist_setup, llist_teardown),
        cmocka_unit_test_setup_teardown(llistPtr_AppendRead, llist_setup, llist_teardown),
        cmocka_unit_test_setup_teardown(llistPtr_DeleteRead, llist_setup, llist_teardown),
    };

    cmocka_run_group_tests(llistPtr, NULL, NULL);
    cmocka_run_group_tests(llistSize, NULL, NULL);
    return 0;
}
