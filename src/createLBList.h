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
    createLBList.h

   Authors:
    Frederic Jarlier, 	Institut Curie
	Nicolas Joly, 		Institut Pasteur
	Nicolas Fedy,		Institut Curie
	Leonor Sirotti,	 	Institut Curie
	Thomas Magalhaes,	Institut Curie
	Paul Paganiban,		Institut Curie
	Tarik id Belouch,	Institut Curie
	Georges Amazo,		Institut Curie
*/

/*
 * Header for createRGList.c
 *
 */
#ifndef CREATELBLIST_H
#define CREATELBLIST_H

#define LBMAX 100
#define CHRMAX 100

typedef struct {
    char* chrList[CHRMAX];
    int chrNum;
    size_t lastChr;
} chrInfo;

typedef struct {
    char* lbList[LBMAX];
    int lbNum;
    size_t lastLb;
} lbInfo;

int createLBList(char *header, int *lb_num, char *lb_list[]);
void initLbInfo(lbInfo* lb, char* header);
int createChrList(char *header, char *chr_list[]) ;
void initChrInfo(chrInfo* chr, char* header);
void freeLbInfo(lbInfo* lb);
void freeChrInfo(chrInfo* chr);
#endif
