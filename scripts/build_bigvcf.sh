#!/bin/bash

in_fold=$1
panel=$2
# out_vcf=$3

tmpf=.tmp_samples

tail -n +2 ${panel} | cut -f 1 | grep -v NA12878 > ${tmpf} # NA12878 hardcoded

for vcf in $(ls ${in_fold}/chr1.vcf.gz)
do
    echo $vcf
    new_vcf=$(dirname ${vcf})/$(basename ${vcf} .vcf.gz).noNA12878.vcf.gz
    bcftools view --threads 8 -Oz -S ${tmpf} ${vcf} > ${new_vcf}
    #if [ ! -f ${out_vcf} ]
    #then
    #bcftools view -S ${tmpf} ${vcf} > ${out_vcf}
    #sed -i '230i\##INFO=<ID=OLD_VARIANT,Number=.,Type=String,Description="Original before vt normalize was run. FORMAT chr:pos:ref:alt">\' ${out_vcf}
    #else
    #bcftools view -S ${tmpf} ${vcf} | grep -v "^#" >> ${out_vcf}
    #fi
done

rm ${tmpf}
