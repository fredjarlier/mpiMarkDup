#!/bin/bash
#: Title : compareSAMWithExternalTools.sh
#: Date : 12/03/2018
#: Author : Julien LIN
#: Version : 2
#: Description : Compare two SAM given in input. First the SAM files are convert to BAM files
#:               by samtools view function and then BAM files are compared with the function
#:               cmpbams of jvarkit. We suppose that the SAM files given in input are coherent.
#:               By default the script suppose that samtools and cmpbams.jar are located in script directory.
#:               If there is no -d option, it supposes script directory is the current directory.
#: Options : -g to indicate that inputs are gzipped
#:           -d directory where the script are. If there is no -d option, it supposes script directory is the current directory.
#:           -s indicates which samtools. By default the script suppose that samtools and cmpbams.jar are located in script directory.
#:           -j indicates which cmpbams.jar. By default the script suppose that samtools and cmpbams.jar are located in script directory.
#:           -h print the help


### Il faut que le script soit capable de prendre un fichier en entrée un fichier gzippé et un autre non

#OPTION STRING
optstring=ghd:s:j:

#OPTION VARIABLE
location=$PWD # problème on appelle pas le script dans le bon répertoire
gzipped=0
directory=$PWD

SAMTOOLS="$directory"/samtools
JVARKIT="$directory"/cmpbams.jar
##Parsing options
while getopts $optstring opt
do
    case $opt in
	g) gzipped=1;;
	h)  echo ""
	    printf "%s\n" "Compare two SAM given in input. First the SAM files are convert to BAM files" "by samtools view function and then BAM files are compared with the function" "cmpbams of jvarkit. We suppose that the SAM files given in input are coherent." "By default the script suppose that samtools and cmpbams.jar are located in script directory." "If there is no -d option, it supposes script directory is the current directory."
	    echo ""
	    printf "%s\n" "This is the option -h of compareSAMWithExternalTools.sh. You may use :" "-g to indicate that inputs are gzipped" "-d directory where the script are. If there is no -d option, it supposes script directory is the current directory." "-h for help" "-s indicates which samtools. By default the script suppose that samtools and cmpbams.jar are located in script directory." "-j indicates which cmpbams.jar. By default the script suppose that samtools and cmpbams.jar are located in script directory."
	    echo ""
	   exit ;;
	d)
	    if [ -z "$OPTARG" ]
	    then
		echo "ERROR ::: option -d should be follow by directory path !!"
		exit 1
	    else
		directory=$OPTARG
		if [[ $directory =~ \.* ]] || [[ $directory =~ */$ ]]
		then
		    directory=`cd $directory 2>/dev/null && pwd`
		fi
		    
	    fi;;
	s)
	    if ! [ -z $OPTARG ]
	    then
		dir=`dirname "$OPTARG"`
		name=`basename "$OPTARG"`
		SAMTOOLS="`cd \"$dir\" 2>/dev/null && pwd || echo \"$dir\"`/$name"
		SAMTOOLS=$OPTARG
	    fi;;
	    
	j) if ! [ -z $OPTARG ]
	   then
	       dir=`dirname "$OPTARG"`
	       name=`basename "$OPTARG"`
	       JVARKIT="`cd \"$dir\" 2>/dev/null && pwd || echo \"$dir\"`/$name"
	       #JVARKIT=$OPTARG
	   fi;;
		   
	*) echo $opt " is not a correct option"
    esac
done

shift "$(( $OPTIND - 1 ))"

##Verifing that input are not empty
if ! ([ -s $1 ] || [ -s $2 ])
then
	echo " ERROR ::: One of the inputs is empty or doesn't exist"
	exit 1
fi

## Converting SAM files into BAM files

dir1=`dirname "$1"`
name1=`basename "$1"`
input1="`cd \"$dir1\" 2>/dev/null && pwd || echo \"$dir1\"`/$name1"
name1="${input1%.sam}_after_samtools.bam"
name1="${name1##*/}"

dir2=`dirname "$2"`
name2=`basename "$2"`
input2="`cd \"$dir2\" 2>/dev/null && pwd || echo \"$dir2\"`/$name2"
name2="${input2%.sam}_after_samtools.bam"
name2="${name2##*/}"

resultsfile="resultsCmpBams_for_${name1%.bam}_${name2%.bam}"
##resultsfile=${resultsfile//[[:space:]]/\\ }

#Test if the ungzipped file if the option is set
if (( $gzipped ))
then
    isgzip_input=`file "$input1"`
    if [[ $isqzip_input =~ .*gzip\ compressed\ data,\ extra\ field ]]
    then
	input1=<(gunzip -c "$input1")
    fi

    isgzip_input=`file "$input2"`
    if [[ $gzip_input =~ .*gzip\ compressed\ data,\ extra\ field ]]
    then
	input2=<(gunzip -c "$input2")
    fi
fi
   
if ! [ -x $SAMTOOLS  ]
then
    chmod +x $SAMTOOLS
fi

mkdir -p "$directory/tempDir"
    
$SAMTOOLS view -bh "$input1" > "$directory/tempDir/$name1" 2>/dev/null
if (( $? > 0 ))
then
    echo "ERROR ::: with the first samtools view" "$name1"
    exit 1
fi

$SAMTOOLS view -bh "$input2" > "$directory/tempDir/$name2" 2>/dev/null

if (( $? > 0 ))
then
    echo "ERROR ::: with the second samtools view" "$name2"
    exit 1
fi


## Comparing BAM files with cmpbams
java -jar $JVARKIT -o "$directory/tempDir/${resultsfile}" "$directory/tempDir/$name1" "$directory/tempDir/$name2" 2>/dev/null

if (( $? >0 ))
then
    echo "ERROR ::: with cmpbams.jar"
    exit 1
fi


## Recolting the results


awk -v nequal=0 -v total=0 'BEGIN {FS = "\t";} 
     $0 ~ /^#/ {next} 
     $2 ~ /EQ/ {total++}
     $2 ~ /NE/ { nequal++; buff[total++]=$0}
     END {print nequal,total; for(reads in buff)
     	 		      		print buff[reads]
}' "$directory/tempDir/$resultsfile"

if (( $? >0 ))
then
    echo "ERROR"
    exit 1
fi

rm "$directory/tempDir/$name1" "$directory/tempDir/$name2" "$directory/tempDir/$resultsfile"
#rm -rf "$directory/tempDir"
exit 0
