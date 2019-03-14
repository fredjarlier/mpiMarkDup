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
    Firmain Martin,     Paris Descartes University
*/


#include "llist.h"

/**
 * @file llist.c
 * @brief  Doubly circular linked list
 */

/**
 * @date 2018 Apr 12
 * @brief Compare two paired-end
 * @details This function is used to determine the order of the insertion in List.
 * @param[in] read1 first read in one pair
 * @param[in] read2 first read in another pair
 * @param[in] isFragment These reads are fragments ?
 * @note
      Documentation from MarkDuplicates (function compare in inner class ReadEndsMDComparator) :
        Comparator for ReadEndsForMarkDuplicates that orders by read1 position then pair orientation then read2 position.
 */

int compareRead(readInfo *read1, readInfo *read2, const int isFragment) {

    //md_log_trace("Start to Compare %s <-> %s\n", read1->Qname, read2->Qname);

    /* 1 : compare read library */
    int compareDifference = read1->readLb - read2->readLb;

    /* 2 : compare read chromosome */
    if (compareDifference == 0) {
        compareDifference = read1->readChromosome - read2->readChromosome;

    } else {
        goto endCompare;
    }

    /* 3 : compare read unclipped coordinate */
    if (compareDifference == 0) {
        compareDifference = read1->unclippedCoordPos - read2->unclippedCoordPos;

    } else {
        goto endCompare;
    }

    //md_log_trace("%s <-> %s, %zu <-> %zu, after unclippedCoordPos, compareDifference=%d\n", read1->Qname, read2->Qname, read1->unclippedCoordPos, read2->unclippedCoordPos, compareDifference);

    orientation read1orientation;
    orientation read2orientation;

    if (isFragment) {
        // get reads fragment orientation, i.e. reverse (R) or forward (F)
        read1orientation = getOrientation(read1, 0);
        read2orientation = getOrientation(read2, 0);

    } else {
        // get reads pair orientation, i.e. FF/FR/RF/RR
        read1orientation = read1->orientation;
        read2orientation = read2->orientation;
    }

    /* 4 : compare read orientation */
    if (compareDifference == 0) {
        compareDifference = read1orientation - read2orientation;

    } else {
        goto endCompare;
    }

    //md_log_trace("%s <-> %s, after orientation, compareDifference=%d\n", read1->Qname, read2->Qname, compareDifference);

    /* 5 : compare mate chromosome */
    if (compareDifference == 0) {
        compareDifference = read1->mateChromosome - read2->mateChromosome;

    } else {
        goto endCompare;
    }

    /* 6 : compare mate position 
     *     note that in this stage, coordMatePos is unclipped
     * */
    if (compareDifference == 0) {
        compareDifference = read1->coordMatePos - read2->coordMatePos;

    } else {
        goto endCompare;
    }

    /* 7 : compare read index in (sorted) file */
    if (compareDifference == 0) {
        compareDifference = read1->indexAfterSort - read2->indexAfterSort;

    } else {
        goto endCompare;
    }

    //md_log_trace("%s <-> %s, after indexAfterSort, compareDifference=%d\n", read1->Qname, read2->Qname, compareDifference);

    /* 8 : compare mate index in (sorted) file */
    if (compareDifference == 0) {
        compareDifference = read1->mateIndexAfterSort - read2->mateIndexAfterSort;

    } else {
        goto endCompare;
    }

endCompare:

//md_log_trace("%s <-> %s, after mateIndexAfterSort, compareDifference\n", read1->Qname, read2->Qname, compareDifference);

//md_log_trace("%s <-> %s, end compareDifference=%d\n", read1->Qname, read2->Qname, compareDifference);
    return compareDifference;
}

/**
 * @date 2018 Apr. 19
 * @brief Allocate and initialize a circular linked list.
 * @return a list initialized
 */

llist_t *llist_create() {
    llist_t *l = malloc(sizeof(llist_t));
    assert(l);
    l->nil = malloc(sizeof(lnode_t));
    assert(l->nil);
    l->nil->prev = l->nil;
    l->nil->next = l->nil;
    l->head = l->nil;
    l->tail = l->nil;
    l->size = 0;
    return l;
}

/**
 * @date 2018 Apr. 19
 * @brief Clear a circular linked list.
 * @details
 *      The returned linked list is come back to its initial state, i.e. A list with one dummy node and zero-size.
 * @param[in, out] l a circular linked list
 * @return a list cleared
 */

void llist_clear(llist_t *l) {
    lnode_t *x = l->nil->next;

    while (x != l->nil) {
        //assert(x != x->next);
        lnode_t *tmp = x;
        //printf("free %s\n", x->read->Qname);
        x = x->next;
        free(tmp);
    }

    l->head = l->nil;
    l->tail = l->nil;
    l->nil->prev = l->nil;
    l->nil->next = l->nil;
    l->size = 0;
}

/**
 * @date 2018 Apr. 19
 * @brief Free a circular linked list.
 * @param[in, out] l a circular linked list
 */

void llist_destruct(llist_t *l) {
    lnode_t *x = l->head;

    while (x != l->nil) {
        lnode_t *tmp = x;
        x = x->next;
        freeRead(tmp->read);
        tmp->read = NULL;
        free(tmp);
    }

    free(l->nil);
    free(l);
}

/**
 * @date 2018 Apr. 19
 * @brief Given an index @p n, do a linear search on a list
 * @param[in, out] l a list
 * @param[in] n index
 * @return the node of index n
 */

lnode_t *llist_lsearch(llist_t *l, int n) {
    assert (n >= 0 || n < l->size) ;

    if (n == l->size - 1) {
        return l->tail;

    } else {
        lnode_t *x = l->nil->next;

        for (int i = 0; i < n; i++) {
            x = x->next;
        }

        return x;
    }
}

/**
 * @date 2018 Apr. 19
 * @brief Delete a node from a list
 * @param[in, out] l a list
 * @param[in, out] p a node
 * @return the node removed
 */

readInfo *llist_delete_ptr(llist_t *l, lnode_t *p) {
    if (p == l->head) {
        l->head = p->next;
    }

    if (p == l->tail) {
        l->tail = p->prev;
    }

    readInfo *tmp = p->read;
    p->prev->next = p->next;
    p->next->prev = p->prev;
    free(p);
    l->size--;
    return tmp;
}

/**
 * @date 2018 Apr. 19
 * @brief Delete the node of index n
 * @param[in,out] l a list
 * @param[in] n index
 * @return read removed from llist
 */

readInfo *llist_delete(llist_t *l, int n) {
    if (l->size == 0) {
        return NULL ;
    }

    lnode_t *x =  llist_lsearch(l, n);
    return llist_delete_ptr(l, x);
}

/**
 * @date 2018 Apr. 19
 * @brief Insert a node @p x before a @p node
 * @param[in,out] node the node to prepend
 * @param[in,out] x the node to insert
 */

static void llist_insert_ptr(lnode_t *node, lnode_t *x) {
    lnode_t *pn = node->prev;
    x->next = pn->next;
    pn->next->prev = x;
    pn->next = x;
    x->prev = pn;
}

/**
 * @date 2018 Apr. 19
 * @brief Given an index n, insert in the front of the node n
 * @param[in,out] l a list
 * @param[in] n index
 * @param[in] e lnode_t
 * @return the new node of index n which be inserted
 */

lnode_t *llist_insert_index(llist_t *l, int n, readInfo *e) {
    lnode_t *x = l->size == 0 ? l->nil : llist_lsearch(l, n);
    lnode_t *node = malloc(sizeof(lnode_t));
    assert(node);
    node->read = e;
    llist_insert_ptr(x, node);

    if (n == 0) {
        l->head = node;
    }

    if (n == l->size) {
        l->tail = node;
    }

    l->size++;
    return node;
}

/**
 * @date 2018 Apr. 19
 * @brief Insert a read @p e before a ginven @p lnode in a list @p l
 * @param[in,out] l a list
 * @param[in, out] lnode
 * @param[in] e a read
 * @return the new node of index n which be inserted
 */

lnode_t *llist_insert_node(llist_t *l, lnode_t *lnode, readInfo *e) {

    lnode_t *node = malloc(sizeof(lnode_t));
    assert(node);
    node->read = e;
    llist_insert_ptr(lnode, node);

    if (lnode == l->head) {
        l->head = node;
    }

    if (lnode == l->nil) {
        l->tail = node;
    }

    l->size++;
    return node;
}

/**
 * @date 2018 Apr. 19
 * @brief Append a read in a circular linked list.
 * @param[in,out] l a list
 * @param[in] e a read
 * @return the new node of index n which be inserted
 */


lnode_t *llist_append(llist_t *l, readInfo *e) {
    lnode_t *node = malloc(sizeof(lnode_t));
    assert(node);
    node->read = e;
    llist_insert_ptr(l->nil, node);

    if (l->size == 0) {
        l->head = node;
    }

    l->tail = node;
    l->size++;
    return node;
}

/**
 * @date 2018 Apr. 19
 * @brief Print a readInfo list
 * @param[in, out] l a list
 * @param[in] comm
 */

void llist_readInfo_print(llist_t *l, MPI_Comm comm) {

    int rank;
    MPI_Comm_rank(comm, &rank);

    /* modify buffer size when needing */
    char buffer[40000];
    int offset = 0;
    offset = sprintf(buffer, "%d nodes : nil<->", l->size);
    lnode_t *x = l->nil->next;

    for (int i = 0; i < l->size; i++) {
        offset += sprintf(buffer + offset, "[%s#%zu]<->", x->read->Qname, x->read->unclippedCoordPos);
        x = x->next;
    }

    sprintf(buffer + offset, "nil\n");
    md_log_rank_error(rank, "%s", buffer);
}

/**
 * @date 2018 Apr. 22
 * @brief Merge sort for circular doubly linked list
 * @param[in, out] l a circular doubly linked list
 * @param[in] isFragment is a fragment list ?
 * @note Adapted from 
 *  https://www.chiark.greenend.org.uk/~sgtatham/algorithms/listsort.html
 */

void llist_merge_sort(llist_t *llist, const int isFragment) {
    lnode_t *p, *q, *e, *tail, *oldhead, *list = llist->head;
    int insize, nmerges, psize, qsize, i;

    if (llist->size == 0) {
        return ;
    }

    insize = 1;

    // relinkage : avoid dummy node
    llist->tail->next = llist->head;
    llist->head->prev = llist->tail;

    while (1) {
        p = list;
        oldhead = list; /* only used for circular linkage */
        list = NULL;
        tail = NULL;

        nmerges = 0;   /* count number of merges we do in this pass */

        while (p) {
            nmerges++;  /* there exists a merge to be done */
            /* step `insize' places along from p */
            q = p;
            psize = 0;

            for (i = 0; i < insize; i++) {
                psize++;

                q = (q->next == oldhead ? NULL : q->next);

                if (!q) {
                    break;
                }
            }

            /* if q hasn't fallen off end, we have two lists to merge */
            qsize = insize;

            /* now we have two lists; merge them */
            while (psize > 0 || (qsize > 0 && q)) {

                /* decide whether next lnode_t of merge comes from p or q */
                if (psize == 0) {
                    /* p is empty; e must come from q. */
                    e = q;
                    q = q->next;
                    qsize--;

                    if (q == oldhead) {
                        q = NULL;
                    }

                } else if (qsize == 0 || !q) {
                    /* q is empty; e must come from p. */
                    e = p;
                    p = p->next;
                    psize--;

                    if (p == oldhead) {
                        p = NULL;
                    }

                } else if (compareRead(p->read, q->read, isFragment) <= 0) {
                    /* First lnode_t of p is lower (or same);
                     * e must come from p. */
                    e = p;
                    p = p->next;
                    psize--;

                    if (p == oldhead) {
                        p = NULL;
                    }

                } else {
                    /* First lnode_t of q is lower; e must come from q. */
                    e = q;
                    q = q->next;
                    qsize--;

                    if (q == oldhead) {
                        q = NULL;
                    }
                }

                /* add the next lnode_t to the merged list */
                if (tail) {
                    tail->next = e;

                } else {
                    list = e;
                }

                /* Maintain reverse pointers in a doubly linked list. */
                e->prev = tail;

                tail = e;
            }

            /* now p has stepped `insize' places along, and q has too */
            p = q;
        }

        tail->next = list;

        list->prev = tail;


        /* If we have done only one merge, we're finished. */
        if (nmerges <= 1) { /* allow for nmerges==0, the empty list case */
            // relinkage with dummy node
            llist->head = list;
            llist->tail = tail;
            llist->tail->next = llist->nil;
            llist->head->prev = llist->nil;
            llist->nil->next = llist->head;
            llist->nil->prev = llist->tail;
            return ;
        }

        /* Otherwise repeat, merging lists twice the size */
        insize *= 2;
    }
}
