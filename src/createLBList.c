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
    createLBList.c

   Authors:
    Frederic Jarlier,   Institut Curie
    Firmain Martin,     Paris Descartes
    Xingwei Sang,       Paris Descartes
*/

#define _GNU_SOURCE
#include <string.h>
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>

#include "createLBList.h"

static int getLBname( char *tokenLine, char *lb_list[], int *lb_num) ;
static int getChrname(char *tokenLine, char *chr_list[], int *chr_num) ;

void initLbInfo(lbInfo* lb, char* header){
    lb->lbNum = 0;
    createLBList(header, &lb->lbNum, lb->lbList);
    lb->lastLb = 0;
}

void initChrInfo(chrInfo* chr, char* header){
    chr->chrNum = createChrList(header, chr->chrList);
    chr->lastChr = 0;
}

void freeLbInfo(lbInfo* lb){
    for (int i = 0; i < lb->lbNum; ++i) {
        free(lb->lbList[i]);
    }
}

void freeChrInfo(chrInfo* chr){
    for (int i = 0; i < chr->chrNum; ++i) {
        free(chr->chrList[i]);
    }
}

/*
 * @date 2018 Jan 27
 * @brief given a header returns the amount of libs and an array of lib name
 * @param[in,out] header a SAM header string
 * @param[out] lb_num total number of lib
 * @param[out] lb_list array of lib name
 * @return always 0
 * @note createLBList use strtok_r to parse the header.
 *       strtok_r is destructive, it replace delimiter by '\0'.
 *       In this function is '\n'.
 */

int createLBList(char *header2, int *lb_num, char *lb_list[]) {

    /*
    * How it works first find the different LB field
    * http://seqanswers.com/forums/showthread.php?t=19540
    *
    *   LB1 -> reads1->read2->read3...
    *   LB2 -> reads4->read5->read6...
    *   LB3 -> reads7->read8->read9...
    *
    *   We add a last LB case when the the mate or the read is unmapped
    *   flag = 137 or 69
    *
    *   69 shall be treated as singleton and 137 as duplicate
    */

    const char *rg_tag = "@RG";
    char *substr_rg = NULL;
    char *tokenLine = NULL;
    char *saveptr;
    int i = 0;

    char* header = strdup(header2);
    //strtok_r is MT-safe while strtok isn't
    tokenLine = strtok_r(header, "\n", &saveptr);

    while (tokenLine) {
        // we test if the tokenLine begin with @RG
        if (strncmp(tokenLine, rg_tag, strlen(rg_tag)) == 0){
        substr_rg = strstr(tokenLine, rg_tag);

            if (substr_rg) {
                // we found a @RG line, look for a LB name
                getLBname(substr_rg, lb_list, lb_num);
            }
        }

        // otherwise we go to next line
        tokenLine = strtok_r(NULL, "\n", &saveptr);

    }

    // special case reads without LB tag
    lb_list[*lb_num] = strdup("unknown");
    (*lb_num)++; 
    // special case unmapped read
    lb_list[*lb_num] = strdup("unmapped");
    (*lb_num)++; 

    free(header);
    return 0;
}

/*
 * @date 2018 Jan 27
 * @brief Given a @RG line, add new LB name into lb_list and increase lb_num.
 * @param[in,out] tokenLine a @RG line from SAM header
 * @param[in,out] lb_list array of lib name
 * @param[in,out] lb_num total number of lib
 * @return always 0
 * @note getLBname use strtok_r to parse tokenLine.
 *       strtok_r is destructive, it replace delimiter by '\0'.
 *       In this function is '\t' after the first substring "LB:".
 * @todo launch a warning if there are multiple LB name in the same @RG line.
 */

static int getLBname(char *tokenLine, char *lb_list[], int *lb_num) {

    const char *lb_tag  = "LB:"; //it is the Library we use not the ID tag
    char *substr_lb     = NULL;
    char *substr_lb2    = NULL;
    int lb_len          = 0;
    char *tokenTab      = NULL;
    char *saveptr;
    int i = 0;

    // we look for LB tag
    substr_lb = strstr(tokenLine, lb_tag);

    if (substr_lb) {
        //we loop the tokenLine searching the LB: string
        tokenTab = strtok_r(substr_lb, "\t", &saveptr);

        while (tokenTab) {
            substr_lb2 = strstr(tokenTab, lb_tag);

            if (substr_lb2) {
                // we found the LB tag, then we extract LB name
                lb_len = strlen(tokenTab) - strlen(lb_tag);

                int lb_exists = 0;

                // verify whether or not LB name exists
                if (*lb_num > 0) {
                    for ( i = 0; i < *lb_num; i++) {
                        if (strncmp(lb_list[i], tokenTab + strlen(lb_tag), lb_len) == 0) {
                            lb_exists = 1;
                        }
                    }
                }

                // LB name doesn't exists, we add it
                if (lb_exists == 0) {
                    lb_list[*lb_num] = strndup(tokenTab + strlen(lb_tag), lb_len);
                    assert(strlen(lb_list[*lb_num]) == lb_len);
                    (*lb_num)++;
                }
            }

            //go to next tokenTab
            tokenTab = strtok_r(NULL, "\t", &saveptr);
        }
    }

    return 0;
}

/*
 * @date 2018 Apr. 18
 * @brief given a header returns the amount of libs and an array of chromosome name
 * @param[in] header a SAM header string
 * @param[out] chr_list array of chromosome name
 * @return chr_num total number of chromosome
 * @note createChrList use strtok_r to parse the header.
 *       strtok_r is destructive, it replace delimiter by '\0'.
 *       In this function is '\n'.
 */

int createChrList(char *header2, char *chr_list[]) {

    const char *sq_tag = "@SQ";
    char *substr_sq = NULL;
    char *tokenLine = NULL;
    char *saveptr;
    int chr_num = 0;
    int i = 0;

    char* header = strdup(header2);
    //strtok_r is MT-safe while strtok isn't
    tokenLine = strtok_r(header, "\n", &saveptr);

    while (tokenLine) {
        // we test if the tokenLine begin with @SQ
        if (strncmp(tokenLine, sq_tag, strlen(sq_tag)) == 0){
        substr_sq = strstr(tokenLine, sq_tag);

            if (substr_sq) {
                // we found a @SQ line, look for a Chr name
                getChrname(substr_sq, chr_list, &chr_num);
            }
        }

        // otherwise we go to next line
        tokenLine = strtok_r(NULL, "\n", &saveptr);

    }

    free(header);
    return chr_num;
}

/*
 * @date 2018 Apr. 18
 * @brief Given a @SQ line, add new Chr name into chr_list and increase chr_num.
 * @param[in,out] tokenLine a @SQ line from SAM header
 * @param[in,out] chr_list array of lib name
 * @param[in,out] chr_num total number of lib
 * @return always 0
 * @note getChrname use strtok_r to parse tokenLine.
 *       strtok_r is destructive, it replace delimiter by '\0'.
 *       In this function is '\t' after the first substring "SN:".
 * @todo launch a warning if there are multiple Chr name in the same @SQ line.
 */

static int getChrname(char *tokenLine, char *chr_list[], int *chr_num) {

    const char *chr_tag  = "SN:"; //it is the Library we use not the ID tag
    char *substr_chr     = NULL;
    char *substr_chr2    = NULL;
    int chr_len          = 0;
    char *tokenTab      = NULL;
    char *saveptr;
    int i = 0;

    // we look for Chr tag
    substr_chr = strstr(tokenLine, chr_tag);

    if (substr_chr) {
        //we loop the tokenLine searching the Chr: string
        tokenTab = strtok_r(substr_chr, "\t", &saveptr);

        while (tokenTab) {
            substr_chr2 = strstr(tokenTab, chr_tag);

            if (substr_chr2) {
                // we found the Chr tag, then we extract Chr name
                chr_len = strlen(tokenTab) - strlen(chr_tag);

                int chr_exists = 0;

                // verify whether or not Chr name exists
                if (*chr_num > 0) {
                    for ( i = 0; i < *chr_num; i++) {
                        if (strncmp(chr_list[i], tokenTab + strlen(chr_tag), chr_len) == 0) {
                            chr_exists = 1;
                        }
                    }
                }

                // Chr name doesn't exists, we add it
                if (chr_exists == 0) {
                    chr_list[*chr_num] = strndup(tokenTab + strlen(chr_tag), chr_len);
                    assert(strlen(chr_list[*chr_num]) == chr_len);
                    (*chr_num)++;
                }
            }

            //go to next tokenTab
            tokenTab = strtok_r(NULL, "\t", &saveptr);
        }
    }

    return 0;
}
