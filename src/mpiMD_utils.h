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
