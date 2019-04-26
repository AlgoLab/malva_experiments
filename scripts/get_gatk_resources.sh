#!/bin/bash

out_fold=$1
base_url=ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg19

function get_vcf {
    fname=$1
    out_vcf=${out_fold}/${fname}
    wget ${base_url}/${fname}.gz -O ${out_vcf}.gz
    gunzip ${out_vcf}.gz
    sed -i 's/=chr/=/g' ${out_vcf}
    sed -i 's/^chr//g' ${out_vcf}
}

for fname in dbsnp_138.hg19.vcf hapmap_3.3.hg19.sites.vcf 1000G_omni2.5.hg19.sites.vcf 1000G_phase1.snps.high_confidence.hg19.sites.vcf Mills_and_1000G_gold_standard.indels.hg19.sites.vcf
do
  get_vcf ${fname}
done
