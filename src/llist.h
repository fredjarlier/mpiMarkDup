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
     llist.h
     
   Authors:
    Frederic Jarlier,   Institut Curie
    Firmain Martin,     Paris Descartes University
*/


#ifndef LLIST_H
#define LLIST_H 
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include "reads.h"
#include "log.h"

/**
 * @file llist.h
 * @brief Circular doubly linked list definition and basic operations
 * @version 0.1
 * @date 01/01/2018
 */

typedef struct lnode_t {
    struct lnode_t* prev; /**< previous node */
    readInfo* read;       /**< data pointer */
    struct lnode_t* next; /**< next node */
} lnode_t;

typedef struct llist_t {
    lnode_t* nil;  /**< sentinel (dummy node) */
    int size;      /**< count element amount  */
    lnode_t* head; /**< head node */
    lnode_t* tail; /**< tail node */
} llist_t;

int compareRead(readInfo *read1, readInfo *read2, const int isFragment) ;
llist_t*  llist_create();
void llist_clear(llist_t *l) ;
lnode_t *llist_append(llist_t *l, readInfo *e) ;
void  llist_destruct(llist_t* l);
void llist_merge_sort(llist_t* llist, const int isFragment) ;
lnode_t *llist_insert_node(llist_t *l, lnode_t* lnode, readInfo *e) ;
lnode_t*  llist_insert_index(llist_t* l, int n, readInfo* e);
readInfo* llist_delete_ptr(llist_t *l, lnode_t* p) ;
readInfo* llist_delete(llist_t* l, int n);
void  llist_readInfo_print(llist_t* l, MPI_Comm comm);
lnode_t*  llist_lsearch(llist_t* l, int n);
#endif /* ifndef LLIST_H */
