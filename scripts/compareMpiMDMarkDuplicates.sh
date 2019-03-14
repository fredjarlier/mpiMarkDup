#!/bin/bash
#: Title : compareMpiMDMarkDuplicates.sh
#: Date : 13/03/2018
#: Author : Julien LIN
#: Version : 3
#: Description : Compare the output of mpiMD and MarkDuplicates. It takes a directory containing directories of class of files test and use compare the files. The script can create a very simple log file which contain the output and the records that don’t match. It creates also a tempory directory in the directory where it is.
#: Options : -h for help
#:           -s indicates which samtools
#:           -p indicates which picard.jar
#:           -j indicates which jvarkit.jar
#:           -S indicates which compareSAMWithExternalTools.sh
#:           -m indicates which mpiMD
#:           -l indicates that a log file is generated. You must indicate the name of the log file. Information are adding at the end of the file

optstring=hs:p:j:S:m:l:

MPIRUN=mpirun 
JAVA=java
SCRIPT_DIR=$PWD
mpiMD="${PWD%trunk/*}trunk/src/mpiMD"
SAMTOOLS=$SCRIPT_DIR/samtools
PICARD=$SCRIPT_DIR/picard.jar
JVARKIT=$SCRIPT_DIR/cmpbams.jar
COMPARE=$SCRIPT_DIR/compareSAMWithExternalTools.sh
LOG_FILE=

##Flag variable
log=0

##Parsing options
while getopts $optstring opt
do
    case $opt in
	h) printf "%s\n" "Compare the output of mpiMD and MarkDuplicates. It takes a directory containing directories of " "class of files test and use compare the files. The script can " "create a very simple log file which contain the output and the records that don’t match." "-h for help" "-s indicates which samtools" "-p indicates which picard.jar" "-j indicates which jvarkit.jar" "-S indicates which compareSAMWithExternalTools.sh" "-m indicates which mpiMD" "-l indicates that a log file is generated"
	   exit ;;
	s)
	    if ! [ -z $OPTARG ]
	    then
		dir=`dirname "$OPTARG"`
		name=`basename "$OPTARG"`
		SAMTOOLS="`cd \"$dir\" 2>/dev/null && pwd || echo \"$dir\"`/$name"
	    else
		echo "ERROR ::: must have the path with the file name of samtools for -s option"
		exit 1
	    fi;;
	    
        p)  if ! [ -z $OPTARG ]
	    then
		dir=`dirname "$OPTARG"`
		name=`basename "$OPTARG"`
		PICARD="`cd \"$dir\" 2>/dev/null && pwd || echo \"$dir\"`/$name"
		#PICARD=$OPTARG
	    else
	       echo "ERROR ::: must have the path with the file name of picard.jar for -p option"
	       exit 1
	    fi;;
	      
	j) if ! [ -z $OPTARG ]
	   then
	       dir=`dirname "$OPTARG"`
	       name=`basename "$OPTARG"`
	       JVARKIT="`cd \"$dir\" 2>/dev/null && pwd || echo \"$dir\"`/$name"
	       #JVARKIT=$OPTARG
	   else
	       echo "ERROR ::: must have the path with the file name of cmpbams.jar for -j option"
	       exit 1
	   fi;;
	
	S) if ! [ -z $OPTARG ]
	   then
	       dir=`dirname "$OPTARG"`
	       name=`basename "$OPTARG"`
	       COMPARE="`cd \"$dir\" 2>/dev/null && pwd || echo \"$dir\"`/$name"
	       #COMPARE=$OPTARG
	   else
	       echo "ERROR ::: must have the path with the file name of compareSAMWithExternalTools.sh for -S option"
	       exit 1
	   fi;;
	
	m) if ! [ -z $OPTARG ]
	   then
	       dir=`dirname "$OPTARG"`
	       name=`basename "$OPTARG"`
	       mpiMD="`cd \"$dir\" 2>/dev/null && pwd || echo \"$dir\"`/$name"
	       #mpiMD=$OPTARG
	   else
	       echo "ERROR ::: must have the path with the name of mpiMD for -m option"
	       exit 1
	   fi;;
	l) log=1
	   if ! [ -z $OPTARG ]
	   then
	       #dir=`dirname "$OPTARG"`
	       #name=`basename "$OPTARG"`
	       #LOG_FILE="`cd \"$dir\" 2>/dev/null && pwd || echo \"$dir\"`/$name"
	       LOG_FILE="$PWD/$OPTARG"
	   else
	       echo "ERROR ::: must have a file name for -l option"
	       exit 1
	   fi;;
	   
	*) echo $opt " is not a correct option"
	   printf "%s\n" "Compare the output of mpiMD and MarkDuplicates. It takes a directory containing directories of " "class of files test and use compare the files. The script can " "create a very simple log file which contain the output and the records that don’t match." "-h for help" "-s indicates which samtools" "-p indicates which picard.jar" "-j indicates which jvarkit.jar" "-S indicates which compareSAMWithExternalTools.sh" "-m indicates which mpiMD" "-l indicates that a log file is generated"
	   exit
    esac
done
shift "$(( $OPTIND - 1 ))"
input=$1

#Searching the real path of the input

dir=`dirname "$input"`
name=`basename "$input"`
input="`cd \"$dir\" 2>/dev/null && pwd || echo \"$dir\"`/$name"

#Initiating variables

separator====================================

mkdir -p "$PWD/tempDir"
TEMP_DIR="$PWD/tempDir"
INVALID_FILE=0
TOTAL_FILE=0
SUCCESSED_FILE=0
invalid_file_name=
INVALID_READS=0
TOTAL_READS=0
ERROR_FILE=

#Start processing
for directory in $(ls $input)
do
    #if [ -d $directory ]; then
	#continue	
    #fi
    cd $input/$directory
    if [ $log -eq 1 ]
    then
	echo $directory >> $LOG_FILE
    fi
    
    echo $directory
    invalid_file=0
    total_file=0
    successed_file=0
    invalid_file_name=
    invalid_reads=0
    total_reads=0
    
    for file in $(ls "$input/$directory")
    do
	if ! [ -f $file ]; then
	    continue	    
	fi
	let "total_file = total_file+1"

	MPIMD_CHR_DIR="$TEMP_DIR/${file}_chromosome_mpiMD"
	MD_CHR_DIR="$TEMP_DIR/${file}_chromosome_MD"
	mkdir -p $MPIMD_CHR_DIR
	mkdir -p $MD_CHR_DIR
	begin=1
	file_result="${file##*/}"
	file_result="${file_result%.sam}_result.sam"
	file_sorted="${file%.sam}_sorted.sam"
	COMPARE_FILE="$TEMP_DIR/compare_file"

	
	###echo $separator "MPISORT" $separator
	$MPIRUN -N 1 $mpiMD $file $MPIMD_CHR_DIR 2>/dev/null
	#$MPIRUN --oversubscribe -N 1 $mpiMD $file $MPIMD_CHR_DIR 2>/dev/null
	exit_num=$?
	if (( $exit_num ))
	then
	    if [ $log -eq 1 ]
	    then
		echo "MPIMD ne marche pas !!" >> $LOG_FILE
		echo "exit_num = "$exit_num >> $LOG_FILE
	    fi
	    if [ $log -eq 1 ]
	    then
		printf "|  Total files  |Tests successed| Tests failed  |  Reads failed |  Total reads  |\n" >> $LOG_FILE
		printf "|---------------|---------------|---------------|---------------|---------------|\n|" >> $LOG_FILE
		printf "%15d|" $TOTAL_FILE $SUCCESSED_FILE $INVALID_FILE $INVALID_READS $TOTAL_READS >> $LOG_FILE
		printf "\n|---------------|---------------|---------------|---------------|---------------|\n" >> $LOG_FILE
	    fi
	    echo "MPIMD ne marche pas !!"
	    echo "exit_num = "$exit_num
	    exit 1
	fi
	

	###echo $separator "MPISORT ENDED" $separator
	
	###echo $separator "MARK DUPPLICATE" $separator

	$SAMTOOLS sort -l 0 -o "$TEMP_DIR/$file_sorted" $file #2>/dev/null
	
	$JAVA -jar $PICARD MarkDuplicates I=$TEMP_DIR/$file_sorted O="$TEMP_DIR/$file_result" M=/dev/null REMOVE_DUPLICATES=false ASSUME_SORTED=true 2>/dev/null

	exit_num=$?
	if (( $exit_num ))
	then
	    if [ $log -eq 1 ]
	    then
	    echo "PB avec MarkDuplicates" >> $LOG_FILE
	    echo "exit_num = "$exit_num >> $LOG_FILE
	    fi
	    if [ $log -eq 1 ]
	    then
		printf "|  Total files  |Tests successed| Tests failed  |  Reads failed |  Total reads  |\n" >> $LOG_FILE
		printf "|---------------|---------------|---------------|---------------|---------------|\n|" >> $LOG_FILE
		printf "%15d|" $TOTAL_FILE $SUCCESSED_FILE $INVALID_FILE $INVALID_READS $TOTAL_READS >> $LOG_FILE
		printf "\n|---------------|---------------|---------------|---------------|---------------|\n" >> $LOG_FILE
	    fi
	    echo "PB avec MarkDuplicates"
	    echo "exit_num = "$exit_num
	    exit 1
	fi
	###echo $separator "MARK DUPPLICATE ENDED" $separator
	###Traiter la séparation du Sam en chromosome

	## Dans un premier temps il faut déterminer les chromosomes présent dans le SAM. Pour cela on peut utiliser -H de samtools view pour avoir le header et utiliser awk pour récupérer la list des chromosomes

	###echo $separator "separation SAM" $separator

	$SAMTOOLS view -H "$TEMP_DIR/$file_result" > "$TEMP_DIR/header_$file_result" #2>/dev/null
	##Liste des chromosomes présents dans le fichier
	SN=$(awk 'BEGIN {FS="\t"}
		        $0 ~ /@SQ/ {
	   		   string=substr($0,match($0,/SN:([[:alnum:]][[:graph:]]*)/))
	   		   string = gensub(/SN:([[:alnum:]][[:graph:]]*)/,"\\1","g",string)
	   		   split(string,a,"\t")
	   		   print a[1]}' "$TEMP_DIR/header_$file_result")

	###echo $separator "LIST chr"

	for chr in $SN
	do
	    #indexing BAM files
	    $SAMTOOLS view -b1 -o "$TEMP_DIR/${file_result%.sam}.bam" "$TEMP_DIR/$file_result" #2>/dev/null
	    $SAMTOOLS index -b "$TEMP_DIR/${file_result%.sam}.bam" #2>/dev/null
	    $SAMTOOLS view -Sh "$TEMP_DIR/${file_result%.sam}.bam" "$chr" > "$MD_CHR_DIR/intermidiaire_chr" #2>/dev/null

	    ###echo $separator "Remove discordant and unmapped" $separator
	    #Remove discordant and unmapped
	    $SAMTOOLS view -H "$MD_CHR_DIR/intermidiaire_chr" > "$MD_CHR_DIR/${file}_$chr" #2>/dev/null
	    awk 'BEGIN {FS="\t"} $7 ~ /=/ && $0 ~ /^[^@]/ {if(!or(and($2,4),and($2,8))) print $0}' "$MD_CHR_DIR/intermidiaire_chr" >> "$MD_CHR_DIR/${file}_$chr" #2>/dev/null

	    ###echo $separator "Remove discordant and unmapped ENDED"
	    #rm "$TEMP_DIR/${file_result%.sam}.bam" "$TEMP_DIR/${file_result%.sam}.index"
	    
	    exit_num=$?
	    if (( $exit_num ))
	    then
		if [ $log -eq 1 ]
		then
		    echo "PB avec samtools view " $chr  >> $LOG_FILE
		    echo "exit_num = "$exit_num >> $LOG_FILE

		fi
		if [ $log -eq 1 ]
		then
		    printf "|  Total files  |Tests successed| Tests failed  |  Reads failed |  Total reads  |\n" >> $LOG_FILE
		    printf "|---------------|---------------|---------------|---------------|---------------|\n|" >> $LOG_FILE
		    printf "%15d|" $TOTAL_FILE $SUCCESSED_FILE $INVALID_FILE $INVALID_READS $TOTAL_READS >> $LOG_FILE
		    printf "\n|---------------|---------------|---------------|---------------|---------------|\n" >> $LOG_FILE
		fi
		echo "PB avec samtools view " $chr
		echo "exit_num = "$exit_num
		exit 1
	    fi
	    
	done

	###echo $separator "separation SAM ENDED" $separator
	
	###echo $separator "compareSAMWithExternalTools.sh" $separator
	
	#echo $mpimd "$TEMP_DIR/$file_result"

	for chrm in $SN
	do
	    $COMPARE -g -s $SAMTOOLS -j $JVARKIT -d $SCRIPT_DIR "$MPIMD_CHR_DIR/${chrm}.gz" "$MD_CHR_DIR/${file}_$chrm" > $COMPARE_FILE
	    read buff nb_reads failed_reads < $COMPARE_FILE
	    exit_num=$?
	    
	    if (( $exit_num ))
	    then
		if [ $log -eq 1 ]
		then
		    echo "PB avec samtools view " $chr >> $LOG_FILE
		    echo "exit_num = "$exit_num >> $LOG_FILE
		fi
		if [ $log -eq 1 ]
		then
		    printf "|  Total files  |Tests successed| Tests failed  |  Reads failed |  Total reads  |\n" >> $LOG_FILE
		    printf "|---------------|---------------|---------------|---------------|---------------|\n|" >> $LOG_FILE
		    printf "%15d|" $TOTAL_FILE $SUCCESSED_FILE $INVALID_FILE $INVALID_READS $TOTAL_READS >> $LOG_FILE
		    printf "\n|---------------|---------------|---------------|---------------|---------------|\n" >> $LOG_FILE
		fi
		echo "PB avec samtools view " $chr
		echo "exit_num = "$exit_num
		exit 1
	    fi
	    #echo "buff = " $buff
	    #buff=${buff//ERROR/1}
	    #reads=${reads//::: One of the inputs is empty or doesn\'t exist/0}
	    #reads=${reads//:::/0}
	    if [[ $buff =~ ERROR ]]
	    then
		if (( $begin ))
		then
		    ERROR_FILE+="$file "
		    begin=0
		fi
		ERROR_FILE+="$chr"
		continue
	    fi
		
	    if (( $buff > 0 ))
	    then
		invalid_reads=$(( $invalid_reads + $buff ))
		total_reads=$(( $total_reads + $nb_reads))
		buff=1
		
		#echo "invalid_file  = "$invalid_file
		#invalid_file_name="${invalid_file_name:+"$invalid_file_name\t"} $file"
		if (( $begin ))
		then
		    let "invalid_file = invalid_file + 1"
		    invalid_file_name+="$file "
		    begin=0
		fi
		
		invalid_file_name+="${chrm}.gz "
		#echo "invalid_file_name = " $invalid_file_name
		if [ $log -eq 1 ]
		then
		    echo $file " : " $chrm >> $LOG_FILE
		    awk 'NR > 1 {print $0}' $COMPARE_FILE  >> $LOG_FILE
		fi
	    else
		total_reads=$(( $total_reads + $nb_reads ))
		#echo "successed_file = " $successed_file
	    fi	    
	done

	if (( $begin ))
	then
	    let "successed_file = successed_file +1"
	fi
	###echo $separator "compareSAMWithExternalTools.sh ENDED" $separator	
	
	#Adding the name to the list of failed tests if the comparaison failed

	#echo $file

	###TODO: faire une partie de suppression

	rm "$MD_CHR_DIR/intermidiaire_chr" $COMPARE_FILE
	
	rm "$TEMP_DIR/header_$file_result" "$TEMP_DIR/${file_result%.sam}.bam" "$TEMP_DIR/${file_result%.sam}.bam.bai" "$TEMP_DIR/$file_result" "$TEMP_DIR/$file_sorted"
	rm -rf $MPIMD_CHR_DIR $MD_CHR_DIR
	
    done
    
    ##rm  $TEMP_DIR
    printf "|  Total files  |Tests successed| Tests failed  |  Reads failed |  Total reads  |\n"
    printf "|---------------|---------------|---------------|---------------|---------------|\n|"
    printf "%15d|" $total_file $successed_file $invalid_file $invalid_reads $total_reads
    printf "\n|---------------|---------------|---------------|---------------|---------------|\n"
    if [ $log -eq 1 ]
    then
	printf "|  Total files  |Tests successed| Tests failed  |  Reads failed |  Total reads  |\n" >> $LOG_FILE
	printf "|---------------|---------------|---------------|---------------|---------------|\n|" >> $LOG_FILE
	printf "%15d|" $total_file $successed_file $invalid_file $invalid_reads $total_reads >> $LOG_FILE
	printf "\n|---------------|---------------|---------------|---------------|---------------|\n" >> $LOG_FILE
    fi

    INVALID_FILE=$(( $INVALID_FILE + $invalid_file ))
    TOTAL_FILE=$(( $TOTAL_FILE + $total_file ))
    SUCCESSED_FILE=$(( $SUCCESSED_FILE + $successed_file ))
    INVALID_READS=$(( $INVALID_READS + $invalid_reads ))
    TOTAL_READS=$(( $TOTAL_READS + $total_reads))

    #printf "Invalid files : \n"
    #printf "%s\n" $invalid_file_name
    #if [ $log -eq 1 ]
    #then
    #printf "Invalid files : \n" >> $LOG_FILE
    #printf "%s\n" $invalid_file_name >> $LOG_FILE
    #fi
done
echo ""
printf "|  Total files  |Tests successed| Tests failed  |  Reads failed |  Total reads  |\n"
printf "|---------------|---------------|---------------|---------------|---------------|\n|"
printf "%15d|" $TOTAL_FILE $SUCCESSED_FILE $INVALID_FILE $INVALID_READS $TOTAL_READS
printf "\n|---------------|---------------|---------------|---------------|---------------|\n"
if [ $log -eq 1 ]
then
    printf "|  Total files  |Tests successed| Tests failed  |  Reads failed |  Total reads  |\n" >> $LOG_FILE
    printf "|---------------|---------------|---------------|---------------|---------------|\n|" >> $LOG_FILE
    printf "%15d|" $TOTAL_FILE $SUCCESSED_FILE $INVALID_FILE $INVALID_READS $TOTAL_READS >> $LOG_FILE
    printf "\n|---------------|---------------|---------------|---------------|---------------|\n" >> $LOG_FILE
fi

exit
