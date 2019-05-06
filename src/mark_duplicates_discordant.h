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

#ifndef MARKDUPLICATESDISCORDANT_H
#define MARKDUPLICATESDISCORDANT_H

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
#include "mark_duplicates.h"
 
void exchangeExternFragDiscordant(llist_t *fragList, llist_t *readEndsList, hashTable *htbl, Interval interval, MPI_Comm comm);
int parseLibrariesDiscordant(char *bufferReads, Interval *intervalByProc, llist_t *fragList, llist_t *readEndsList, readInfo ***readArr, char ***samTokenLines, size_t readNum, size_t readIndex, chrInfo *chr, lbInfo *lb, MPI_Comm comm);
readInfo *buildReadEndsDiscordant(readInfo *read1, readInfo *read2, llist_t *readEndsList, int case_insert);
char *markDuplicateDiscordant (char *bufferReads, size_t readNum, char *header, MPI_Comm comm, char *chrName);
void findDuplicaDiscordant(llist_t *fragList, llist_t *readEndsList, hashTable *htbl, int *totalDuplica, int *totalOpticalDuplicate, MPI_Comm comm);

#endif
