#!/bin/bash

fold=$1

echo "Tool,Chr,RSNP,PSNP,RIND,PIND"

for csv in $(ls ${fold}/*/*.summary.csv)
do
    chr=$(basename ${csv} .summary.csv)
    tool=$(basename $(dirname ${csv}))
    snp=$(grep "SNP,ALL" ${csv} | cut -f 11,12 -d',')
    ind=$(grep "INDEL,ALL" ${csv} | cut -f 11,12 -d',')
    echo ${tool},${chr},${snp},${ind}
done
