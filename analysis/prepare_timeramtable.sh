#!/bin/bash

fold=$1

echo "Tool,Time,RAM"

for tool in malva vargeno discosnp gatk bcftools bwa
do
    if [ -d ${fold}/${tool} ]
    then
	if [[ ${tool} == "bwa" ]]
	then
	    echo "bwa"
	else
	    time=0
	    ram=0
	    for log in $(ls ${fold}/${tool}/*.time)
	    do
		t=$(grep "wall" ${log} | cut -f 8 -d' ' | awk -F: '{ if (NF == 1) {print $NF} else if (NF == 2) {print $1 * 60 + $2} else if (NF==3) {print $1 * 3600 + $2 * 60 + $3} }' | cut -f 1 -d'.')
		r=$(grep "Maxim" ${log} | cut -f 6 -d' ')
		((r > ram)) && ram=${r}
		time=$((time + t))
	    done
	    echo $tool,$time,$ram
	fi
    else
	echo $tool,0,0
    fi
done
