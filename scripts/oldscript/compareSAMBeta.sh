#!/bin/bash
if [ -z $1 ] || [ -z $2 ]
then
	echo "input must be 2 sam files "
	exit
fi
filename1=$(basename "$1")
filename2=$(basename "$2")
extension1=`echo "$filename1" | awk -F"." '{print $NF}'`
extension2=`echo "$filename2" | awk -F"." '{print $NF}'`
if [ $extension1 != "sam" ] || [ $extension2 != "sam" ]
then
	echo " wrong extension "
       exit
fi       
let "comp = 1"
nbOfReads1=`cat $1|grep -E ^[^@]|awk ' END {print NR}'`
nbOfReads2=`cat $2|grep -E ^[^@]|awk ' END {print NR}'`
let "nbOfReads=$nbOfReads1"
if [ $nbOfReads1 != $nbOfReads2 ]
then
    echo "number of reads are different"
    exit
fi
reads1=`grep -E ^[^@] $1 |sort -k3,3 -k4,4 -k1,1`
reads2=`grep -E ^[^@] $2 |sort -k3,3 -k4,4 -k1,1`
echo  "\documentclass{article}
\usepackage[utf8]{inputenc}
\usepackage[francais]{babel}
\usepackage[top=2cm, bottom=2cm, left=2cm, right=2cm]{geometry}
\usepackage[T1]{fontenc}
\usepackage{amsmath}
\usepackage{amssymb}
\usepackage{stmaryrd} %for llbracket
\usepackage{caption}
\usepackage{array}
\usepackage{float}
\usepackage{imakeidx}
\usepackage{cite}
\usepackage{hyperref}
\usepackage[all]{hypcap}
\usepackage{longtable}
\usepackage{fancyhdr}
\pagestyle{fancy}
\fancyhead[R]{}        %en-tête droite
\fancyfoot[C]{}
\begin{document}
\section{Différence entre ${filename1} ${filename2}}
\begin{center}
\begin{longtable}{|l|l|l|l|l|l|l|l|l|l|l|l|}
\hline reads&QNAME&Flag&Chromosome&Pos&MAPQ&CIGAR&=&RNEXT&PNEXT&Sequence&Score\\\\
" > compareSamOutput/compareSAMBetaOutput.tex
rm -f compareSamOutput/hline
rm -f compareSamOutput/antislash
rm -f compareSAMBetaOutput/null1
rm -f  compareSAMBetaOutput/null2
rm -f compareSamOutput/reads
for  i in `seq 1 $nbOfReads`
do
    echo "\hline" >> compareSamOutput/hline
    echo " \\\\" >> compareSamOutput/antislash
    echo "$i" >> compareSamOutput/reads
done
let "comp=1"
for i in 'QNAME' 'FLAG' 'CHROMOSOME' 'POS' 'MAPQ' 'CIGAR' '=' 'RNEXT' 'PNEXT' 'SEQUENCE' 'SCORE' 
do
    rm -f compareSamOutput/$i
    champ1=`echo  -e "$reads1"|cut -f $comp`
    champ2=`echo  -e "$reads2"|cut -f $comp`
    let "indice =1"
    for j in $champ1
        do
            table1[$indice]=$j
            let "indice=indice+1"
        done
    let "indice=1"
    for j in $champ2
        do
            table2[$indice]=$j
            let "indice=indice+1"
        done
    for j in `seq 1 $nbOfReads`;
        do
            result=`diff -q <(echo -e "${table1[$j]}") <(echo -e "${table2[$j]}")`
            if [ "$result" != "" ]
                then
                    echo  "x" >> compareSamOutput/$i
                else
                    echo  "\checkmark" >> compareSamOutput/$i
            fi

        done
    let "comp=comp+1"
done
pr -mts compareSamOutput/hline compareSamOutput/reads > compareSamOutput/null1
pr -mts\& compareSamOutput/null1 compareSamOutput/QNAME compareSamOutput/FLAG compareSamOutput/CHROMOSOME compareSamOutput/POS compareSamOutput/MAPQ compareSamOutput/CIGAR compareSamOutput/= compareSamOutput/RNEXT compareSamOutput/PNEXT compareSamOutput/SEQUENCE compareSamOutput/SCORE > compareSamOutput/null2
pr -mts compareSamOutput/null2 compareSamOutput/antislash >>compareSamOutput/compareSAMBetaOutput.tex
echo "\hline
\end{longtable}
\end{center}
\end{document}" >> compareSamOutput/compareSAMBetaOutput.tex
latex compareSamOutput/compareSAMBetaOutput.tex
