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
     mpiMD_utils.h
     
   Authors:
    Frederic Jarlier,   Institut Curie
    Firmain Martin,     Paris Descartes University
*/

#ifndef MPIMD_UTILS_H
#define MPIMD_UTILS_H 

#include <stdlib.h>
#include <assert.h>
#include <math.h>

#include "log.h"
#include "reads.h"

char digit2char(int n) ;
int fastqToPhred (char ch) ;
int getLine(char** offset, char* text, char** tokenLine);
int getTokenTab(char **offset, char *text, char **tokenTab) ;
char* strapp(char* dest, char* src ) ;
unsigned int readBits(unsigned int x, int k) ;
void writeBits(unsigned int *x, int k) ;
void toggleBits(unsigned int *x, int k) ;
void switchBits(unsigned int *x, int k, int j);
#endif /* ifndef MPIMD_UTILS_H */
