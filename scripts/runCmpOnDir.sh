#!/bin/sh

# Author : Firmin Martin
# Last modification : 2018 Apr. 14
# Given a directory containing sam file, compare output of MarkDuplicates and mpiMD
# SYNOPSIS 
#       ./runCmpOnDir.sh DIR  
# NOTE
#       User need modify *user variables* of cmpWithMD before using

LOGNAME=`date '+%Y-%m-%d_%H:%M:%S'`cmp.log

########################
######## SCRIPT ########
########################

FILES="$1/*.sam"

DIRSUM=0
DIREQ=0
DIRNE=0

for f in $FILES
do 
   echo "Processing `basename $f` ..."  
   . ./cmpWithMD.sh $f 2> /dev/null
   DIRSUM=$(( DIRSUM + SUM ))
   DIREQ=$(( DIREQ + EQ ))
   DIRNE=$(( DIRNE + NE ))
   echo -e "`basename $f`\t$SUM\t$EQ\t$NE\t$OPTICAL" >> "$LOGNAME"
done
PERC="`awk "BEGIN { print ($DIREQ/$DIRSUM)*100 }"`%"
STAT="\ntotal reads = $DIRSUM, equal = $DIREQ, not equal = $DIRNE, equal ratio = $PERC"
for (( i = 0; i < ${#STAT}; i++ )); do
    printf "-"  
done
echo -e $STAT
echo -e "\nSum\t$DIRSUM\t$DIREQ\t$DIRNE\t$DIROPTICAL" >> "$LOGNAME"
