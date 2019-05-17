#!/bin/bash

out_vcf=$1
threads=$2
sample_to_remove=NA12878

base_url="ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502"
panel="integrated_call_samples_v3.20130502.ALL.panel"
vcfs=("ALL.chr1.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz"
      "ALL.chr2.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz"
      "ALL.chr3.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz"
      "ALL.chr4.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz"
      "ALL.chr5.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz"
      "ALL.chr6.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz"
      "ALL.chr7.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz"
      "ALL.chr8.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz"
      "ALL.chr9.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz"
      "ALL.chr10.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz"
      "ALL.chr11.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz"
      "ALL.chr12.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz"
      "ALL.chr13.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz"
      "ALL.chr14.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz"
      "ALL.chr15.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz"
      "ALL.chr16.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz"
      "ALL.chr17.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz"
      "ALL.chr18.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz"
      "ALL.chr19.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz"
      "ALL.chr20.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz"
      "ALL.chr21.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz"
      "ALL.chr22.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz"
      "ALL.chrX.phase3_shapeit2_mvncall_integrated_v1b.20130502.genotypes.vcf.gz")

pushd $(dirname ${out_vcf})

# Setting up samples
samples_tmp=.samples.tmp
file_list=.filelist.tmp
wget ${base_url}/${panel}
tail -n +2 ${panel} | cut -f 1 | grep -v ${sample_to_remove} > ${samples_tmp}

for vcf in ${vcfs[*]}
do
    name=$(echo $vcf | cut -f 2 -d '.')
    echo ${name}

    # Downloading VCF
    wget ${base_url}/${vcf} -O ${name}.vcf.gz
    wget ${base_url}/${vcf}.tbi -O ${name}.vcf.gz.tbi

    # Removing sample
    bcftools view --threads ${threads} -Oz -S ${samples_tmp} ${name}.vcf.gz > ${name}.no${sample_to_remove}.vcf.gz
    tabix -p vcf ${name}.no${sample_to_remove}.vcf.gz
    echo ${name}.no${sample_to_remove}.vcf.gz >> ${file_list}
done

bcftools concat -o ${out_vcf} -Oz --threads ${threads} -f ${file_list}
tabix -p vcf ${out_vcf}

rm ${samples_tmp}
rm ${file_list}

popd
