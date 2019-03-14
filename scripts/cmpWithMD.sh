#!/bin/sh


#   mpiSORT
#   Copyright (C) 2016-2019 Institut Curie / Institut Pasteur
#   mpiSORT is free software; you can redistribute it and/or
#   modify it under the terms of the GNU Lesser General Public
#   License as published by the Free Software Foundation; either
#   version 2.1 of the License, or (at your option) any later version.
#   mpiSORT is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#   Lesser General Public License for more details.
#   You should have received a copy of the GNU Lesser Public License
#   along with mpiSORT.  If not, see <http://www.gnu.org/licenses/>.

#    Authors:
#    Frederic Jarlier,   Institut Curie
#    Firmin Martin,      Paris Descartes    

# Last modification : 2018 Apr. 22
#
# Given a sam file, compare output of MarkDuplicates and mpiMD
#
#
# SYNOPSIS 
#       ./cmpWithMD.sh SAM  
#
# NOTE
#       User need modify *user variables* before using
#
# Results Directory Hierarchy :
#
#   ./{filename}
#   |__ MD
#   |   |__ {filename}.md.sam
#   |
#   |__ mpiMD
#   |   |__ {chr1}.gz
#   |   ...    ...
#   |   |__ {chrn}.gz
#   |   |__ unmapped.gz
#   |   |__ discordant.gz
#   |
#   |__ results
#       |__ {chr1}.jivarkit
#       ...    ...
#       |__ {chrn}.jivarkit
#       |__ discordants.jivarkit
#       

############################
###### USER VARIABLES ######
############################

# program to mark
PICARDTOOLS=./picard.jar

# mpiMD parameters
MPI_MD=../src/mpiMD
NUM_PROCS=1
VERBOSITY=3

# program used to testing
JVARKIT=./cmpbams.jar
SAMTOOLS=/usr/local/bin/samtools

###########################
######### Handling ########
###########################


TMP_DIR=`mktemp -d`
INPUT="$1"

# get filename without extension
filename=$(basename -- "$INPUT")
extension="${filename##*.}"
filename="${filename%.*}"

OUTDIR=$filename
mpiMD_OUT=$OUTDIR/mpiMD
MD_OUT=$OUTDIR/MD
JV_OUT=$OUTDIR/results

mkdir -p $OUTDIR
mkdir -p $mpiMD_OUT
mkdir -p $MD_OUT
mkdir -p $JV_OUT

# samtools output
SORTED_FILE=$TMP_DIR/${filename}_sorted.sam
SORTED_FILE2=$TMP_DIR/${filename}_sorted2.sam
HEADER=`$SAMTOOLS view -H $INPUT`

# markDuplicates output
## output with discordants
MD_FILE=$MD_OUT/${filename}.md.sam
## output without discordants
MD_FILE2=$TMP_DIR/${filename}.md2.sam
MD_BAM=$TMP_DIR/${filename}.md.bam
DISCORDANTS_MD_BAM=$TMP_DIR/discordants.md.bam


#######################################
########  MarkDuplicates Part  ########
#######################################
# 1. Sort sam with samtools           #
# 2. Call MarkDuplicates              #
# 3. Extract and remove discordants   #
# 4. Extract reads by chromosome      #
# 5. Convert extracted files to bam   #
# ----------------------------------  #
# As results, we have in $TMP_DIR :   #
#   <chr1>_MD.bam                     #
#        ...                          #
#   <chrn>_MD.bam                     #
#   discordants_MD.bam                #
#######################################

# sort sam
$SAMTOOLS sort -l 0 -o "$SORTED_FILE" "$INPUT" 

# call markDuplicates
java -jar $PICARDTOOLS MarkDuplicates I="$SORTED_FILE" O=$MD_FILE M=$TMP_DIR/${filename}.metric PROGRAM_RECORD_ID=null TAG_DUPLICATE_SET_MEMBERS=true VERBOSITY=DEBUG TAGGING_POLICY=OpticalOnly 2> /dev/null

if [ $? -ne 0 ]; then
    echo "MarkDuplicates failed."
    exit(1);
fi

## extract discordants file
echo "$HEADER" > "$TMP_DIR/discordants_MD.sam"
sed '/^@.*$/d' "$MD_FILE" | awk 'BEGIN {FS="\t"}; $7 !~ /=/ && $6 !~ /\*/ {print $0}' >> "$TMP_DIR/discordants_MD.sam"
$SAMTOOLS view -Sb "$TMP_DIR/discordants_MD.sam" > $DISCORDANTS_MD_BAM

## remove discordants
echo "$HEADER" > "$MD_FILE2"
sed '/^@.*$/d' "$MD_FILE" | awk 'BEGIN {FS="\t"}; $7 ~ /=/ {print $0}' >> "$MD_FILE2"

# convert the result to bam (remove (read unmapped + mate unmapped)
samtools view -Sb "$MD_FILE2" -G 0xC > $MD_BAM

# get list of chromosome i.e. RNAME
chr=`grep "^@SQ" "$MD_FILE2" | cut -f2 | awk 'BEGIN { FS="SN:" };  {print $2 };'`

# Convert chr list to an array

## Save current IFS
SAVEIFS=$IFS
## Change IFS to new line. 
IFS=$'\n'
chr=($chr)
## Restore IFS
IFS=$SAVEIFS

# Divide markDuplicates output by chromosome

## index bam
$SAMTOOLS index $MD_BAM $TMP_DIR/${filename}.md.bam.bai

## extract reads by chromosome file
for (( i=0; i<${#chr[@]}; i++ ))
do
    $SAMTOOLS view $MD_BAM ${chr[$i]} -b > $TMP_DIR/${chr[$i]}_MD.bam
done

##########################
####### mpiMD Part #######
##########################

mpirun --oversubscribe -N $NUM_PROCS $MPI_MD $INPUT $mpiMD_OUT -v $VERBOSITY

if [ $? -ne 0 ]; then
    echo "mpiMD failed."
    exit(1);
fi

#### Comparisons ####

EQ=0
NE=0
SUM=0

## Handle chromosomes
for (( i=0; i<${#chr[@]}; i++ ))
do
    gunzip -c $mpiMD_OUT/${chr[$i]}.gz > "$TMP_DIR/${chr[$i]}_mpiMD.sam"
    $SAMTOOLS view "$TMP_DIR/${chr[$i]}_mpiMD.sam" -b > $TMP_DIR/${chr[$i]}_mpiMD.bam
    java -jar $JVARKIT -o "$JV_OUT/${chr[$i]}.jvarkit" "$TMP_DIR/${chr[$i]}_MD.bam" "$TMP_DIR/${chr[$i]}_mpiMD.bam"
    # sed is used to ignore first line of jivarkit results, this line 
    # contains hashed tmp directory name which may coincide with "NE" or "EQ".
    NE=$(( NE + $((`sed '/^#.*$/d' $JV_OUT/${chr[$i]}.jvarkit | grep "NE" |  wc --lines`)) ))
    EQ=$(( EQ + $((`sed '/^#.*$/d' $JV_OUT/${chr[$i]}.jvarkit | grep "EQ"  | wc --lines`)) ))
done

## Handle discordants

gunzip -c $mpiMD_OUT/discordant.gz > "$TMP_DIR/discordants_mpiMD.sam"
$SAMTOOLS view "$TMP_DIR/discordants_mpiMD.sam" -b > $TMP_DIR/discordants_mpiMD.bam
java -jar $JVARKIT -o "$JV_OUT/discordants.jvarkit" "$DISCORDANTS_MD_BAM" "$TMP_DIR/discordants_mpiMD.bam"
NE=$(( NE + $((`sed '/^#.*$/d' $JV_OUT/discordants.jvarkit | grep "NE" | wc --lines`)) ))
EQ=$(( EQ + $((`sed '/^#.*$/d' $JV_OUT/discordants.jvarkit | grep "EQ" | wc --lines`)) ))

SUM=$(( EQ + NE ))

export SUM
export EQ
export NE

echo "total reads : $SUM, equal : $EQ, not equal : $NE"

##########################
##### Clean tmp dir ######
##########################

rm -rf $TMP_DIR

