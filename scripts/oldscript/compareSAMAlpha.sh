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
reads1=`grep -E ^[^@] $1 |sort -k3,3 -k4,4 -k1,1`
reads2=`grep -E ^[^@] $2 |sort -k3,3 -k4,4 -k1,1`
for i in `seq 1 11`
do
	a=$(diff  <(echo -e "$reads1"|cut  -f $i) <(echo -e "$reads2"|cut  -f $i))
	if [ "$a" != "" ]
       then
	       echo "1"
	       exit
       fi
done
echo "0"

