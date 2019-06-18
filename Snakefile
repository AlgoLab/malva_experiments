configfile: "config.yaml"

import os
pjoin = os.path.join
WF = os.getcwd()

### PARAMETERS
parameters = config["params"]
truth_sample = parameters["truth"]
pop = parameters["pop"]
k = parameters["k"]
K = parameters["K"]
chroms = parameters["chroms"]
chroms_spacesep = chroms.replace(',', ' ')

root_fold = config["root"]

### FOLD NAMES
in_fold_name = config["folds"]["in"]
mid_fold_name = config["folds"]["mid"]
out_fold_name = config["folds"]["out"]

### FILE NAMES
fnames = config["files"]
ref_fname = fnames["ref"]
cleanref_fname = fnames["cleanref"]
vcf_fname = fnames["vcf"]
vcf_vg_fname = fnames["vcf_vg"]
sample_name = fnames["sample"]
sample_bam_fname = sample_name + ".bam"
sample_fname = sample_name + ".fq"

# GATK resources
gatkres_fold = pjoin(root_fold, "FULL", in_fold_name, "gatk")
gatk_snp_res = [pjoin(gatkres_fold, "FULL", "dbsnp_138.hg19.vcf"),
                pjoin(gatkres_fold, "FULL", "dbsnp_138.hg19.vcf.idx"),
                pjoin(gatkres_fold, "FULL", "hapmap_3.3.hg19.sites.vcf"),
                pjoin(gatkres_fold, "FULL", "hapmap_3.3.hg19.sites.vcf.idx"),
                pjoin(gatkres_fold, "FULL", "1000G_omni2.5.hg19.sites.vcf"),
                pjoin(gatkres_fold, "FULL", "1000G_omni2.5.hg19.sites.vcf.idx"),
                pjoin(gatkres_fold, "FULL", "1000G_phase1.snps.high_confidence.hg19.sites.vcf"),
                pjoin(gatkres_fold, "FULL", "1000G_phase1.snps.high_confidence.hg19.sites.vcf.idx")]
gatk_indels_res = [pjoin(gatkres_fold, "FULL", "Mills_and_1000G_gold_standard.indels.hg19.sites.vcf"),
                   pjoin(gatkres_fold, "FULL", "Mills_and_1000G_gold_standard.indels.hg19.sites.vcf.idx")]

rule run:
    input:
        # Happy
        expand(pjoin(root_fold, "FULL", out_fold_name, "happy", "{tool}.results.summary.csv"),
               tool = ["malva", "discosnp", "bcftools", "gatk"]),
        expand(pjoin(root_fold, "HALF", out_fold_name, "happy", "{tool}.results.summary.csv"),
               tool = ["malva", "vargeno", "discosnp", "bcftools", "gatk"]),
        # Bar charts, histograms, and scatterplots
        pjoin(root_fold, "plots", "PR_plot.FULL.pdf"),
        pjoin(root_fold, "plots", "PR_plot.HALF.pdf"),
        pjoin(root_fold, "plots", "indels_plot.pdf"),
        pjoin(root_fold, "plots", "resources_plot.pdf"),
        # Heatmaps and tables
        pjoin(root_fold, "plots", "malva.heatmap.pdf"),
        pjoin(root_fold, "plots", "malva.table.pdf"),
        pjoin(root_fold, "plots", "vargeno.heatmap.pdf"),
        pjoin(root_fold, "plots", "vargeno.table.pdf")

##################
### FETCH DATA ###
##################
rule fetch_data:
    input:
        pjoin(root_fold, "FULL", in_fold_name, ref_fname),
        pjoin(root_fold, "FULL", in_fold_name, ref_fname + ".fai"),
        pjoin(root_fold, "FULL", in_fold_name, ref_fname[:-3] + ".dict"),
        pjoin(root_fold, "FULL", in_fold_name, vcf_fname),
        pjoin(root_fold, "FULL", in_fold_name, sample_bam_fname),
        gatk_snp_res[::2],
        gatk_indels_res[0]

rule fetch_reference:
    input:
    output:
        pjoin(root_fold, "FULL", in_fold_name, ref_fname)
    params:
        url = "ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz"
    shell:
        """
        wget {params.url} -O {output}.gz
        gunzip {output}.gz
        """

rule index_reference:
    input:
        pjoin(root_fold, "FULL", in_fold_name, ref_fname)
    output:
        ref_index = pjoin(root_fold, "FULL", in_fold_name, ref_fname + ".fai"),
        ref_dict = pjoin(root_fold, "FULL", in_fold_name, ref_fname[:-3] + ".dict")
    conda:
        "envs/malva_exps.yaml"
    shell:
        """
        gatk CreateSequenceDictionary -R {input} -O {output.ref_dict}
        samtools faidx {input}
        """

rule fetch_vcf:
    input:
    output:
        vcf = pjoin(root_fold, "FULL", in_fold_name, vcf_fname),
        zvcf = pjoin(root_fold, "FULL", in_fold_name, vcf_fname + ".gz")
    params:
        fold = pjoin(root_fold, "FULL", in_fold_name)
    shell:
        """
        bash {WF}/scripts/get_vcf.sh {params.fold}
        """

rule fetch_bam:
    input:
    output:
        pjoin(root_fold, "FULL", in_fold_name, sample_bam_fname)
    params:
        url = "ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/NA12878/NIST_NA12878_HG001_HiSeq_300x/RMNISTHS_30xdownsample.bam"
    shell:
        """
        wget {params.url} -O {output}
        wget {params.url}.bai -O {output}.bai
        """

rule fetch_gatk_resources:
    input:
    output:
        gatk_snp_res[::2],
        gatk_indels_res[0]
    params:
        gatkres_fold
    shell:
        """
        bash {WF}/scripts/get_gatk_resources.sh {params}
        """

##################
### SETUP DATA ###
##################
'''
This rule prepare the input for the next steps - It populates mid folder with:
 - clean reference (vargeno cannot manage Ms and Rs in chromosome 3)
 - vcf for malva, ie without NA12878 sample
 - vcf for vargeno, ie without samples and format field in the header
 - vcf of NA12878 (truth) extracted from vcf
 - sample (extracted from BAM)
'''
rule clean_ref_from_IUPAC:
    input:
        pjoin(root_fold, "FULL", in_fold_name, ref_fname)
    output:
        cref = pjoin(root_fold, "FULL", mid_fold_name, cleanref_fname),
        cref_index = pjoin(root_fold, "FULL", mid_fold_name, cleanref_fname + ".fai"),
        cref_dict = pjoin(root_fold, "FULL", mid_fold_name, cleanref_fname[:-3] + ".dict")
    log:
        pjoin(root_fold, "FULL", mid_fold_name, cleanref_fname + ".removed.log" )
    conda:
        "envs/malva_exps.yaml"
    shell:
        """
        python3 {WF}/scripts/clean_fasta.py {input} > {output.cref} 2> {log}
        gatk CreateSequenceDictionary -R {input} -O {output.cref_dict}
        samtools faidx {output.cref}
        """

rule prepare_vcf_for_vg:
    input:
        vcf = pjoin(root_fold, "FULL", in_fold_name, vcf_fname),
        zvcf = pjoin(root_fold, "FULL", in_fold_name, vcf_fname + ".gz")
    output:
        vcf = pjoin(root_fold, "FULL", mid_fold_name, vcf_vg_fname),
        zvcf = pjoin(root_fold, "FULL", mid_fold_name, vcf_vg_fname + ".gz"),
    conda:
        "envs/malva_exps.yaml"
    shell:
        """
        python3 {WF}/scripts/format_vcf_as_dbsnp.py {input.vcf} > {output.vcf}
        bgzip -c {output.vcf} > {output.zvcf}
        tabix -p vcf {output.zvcf}
        """

rule extract_truth:
    input:
        vcf = pjoin(root_fold, "FULL", in_fold_name, vcf_fname),
        zvcf = pjoin(root_fold, "FULL", in_fold_name, vcf_fname + ".gz")
    output:
        vcf = pjoin(root_fold, "FULL", mid_fold_name, truth_sample + ".vcf"),
        zvcf = pjoin(root_fold, "FULL", mid_fold_name, truth_sample + ".vcf.gz"),
    params:
        indiv = truth_sample
    conda:
        "envs/malva_exps.yaml"
    shell:
        """
        bcftools view -s {params.indiv} {input.zvcf} | grep -P -v '\t<' > {output.vcf}
        bgzip -c {output.vcf} > {output.zvcf}
        tabix -p vcf {output.zvcf}
        """

rule extract_sample:
    input:
        pjoin(root_fold, "FULL", in_fold_name, sample_bam_fname)
    output:
        pjoin(root_fold, "FULL", mid_fold_name, sample_fname)
    conda:
        "envs/malva_exps.yaml"
    shell:
        """
        samtools fastq -F 0x900 {input} > {output}
        """

##################
### HALVE DATA ###
##################
rule halve_reference:
    input:
        pjoin(root_fold, "FULL", mid_fold_name, cleanref_fname)
    output:
        cref = pjoin(root_fold, "HALF", mid_fold_name, cleanref_fname),
        cref_index = pjoin(root_fold, "HALF", mid_fold_name, cleanref_fname + ".fai"),
        cref_dict = pjoin(root_fold, "HALF", mid_fold_name, cleanref_fname[:-3] + ".dict")
    params:
        chroms
    conda:
        "envs/malva_exps.yaml"
    shell:
        """
        samtools faidx {input} {chroms_spacesep} > {output.cref}
        gatk CreateSequenceDictionary -R {input} -O {output.cref_dict}
        samtools faidx {output.cref}
        """

rule halve_vcfs:
    input:
        zvcf = pjoin(root_fold, "FULL", in_fold_name, vcf_fname + ".gz"),
        zvcf_vg = pjoin(root_fold, "FULL", mid_fold_name, vcf_vg_fname + ".gz"),
        zvcf_truth = pjoin(root_fold, "FULL", mid_fold_name, truth_sample + ".vcf.gz")
    output:
        hvcf = pjoin(root_fold, "HALF", in_fold_name, vcf_fname),
        hvcf_vg = pjoin(root_fold, "HALF", mid_fold_name, vcf_vg_fname),
        hvcf_truth = pjoin(root_fold, "HALF", mid_fold_name, truth_sample + ".vcf")
    params:
        chroms
    conda:
        "envs/malva_exps.yaml"
    shell:
        """
        bcftools view -r {params} {input.zvcf} -o {output.hvcf}
        bcftools view -r {params} {input.zvcf_vg} -o {output.hvcf_vg}
        bcftools view -r {params} {input.zvcf_truth} -o {output.hvcf_truth}
        """

rule halve_sample:
    input:
        pjoin(root_fold, "FULL", in_fold_name, sample_bam_fname)
    output:
        pjoin(root_fold, "HALF", mid_fold_name, sample_fname)
    params:
        pjoin(root_fold, "HALF", "tmp.bam")
    conda:
        "envs/malva_exps.yaml"
    shell:
        """
        samtools view -@ 7 -b -h -F 0x900 {input} {chroms_spacesep} > {params}
        samtools index {params}
        samtools fastq {params} > {output}
        rm {params} {params}.bai
        """

### MALVA ###
rule kmc:
    input:
        pjoin(root_fold, "{run}", mid_fold_name, sample_fname)
    output:
        pjoin(root_fold, "{run}", out_fold_name, "malva", "KMC.kmc_suf")
    params:
        kmc_prefix = pjoin(root_fold, "{run}", out_fold_name, "malva", "KMC"),
        tmp_fold = pjoin(root_fold, "{run}", out_fold_name, "malva", "KMC_tmp")
    log:
        time = pjoin(root_fold, "{run}", out_fold_name, "malva", "KMC.time"),
        out = pjoin(root_fold, "{run}", out_fold_name, "malva", "KMC.log")
    conda:
        "envs/malva_exps.yaml"
    shell:
        """
        mkdir -p {params.tmp_fold}
        /usr/bin/time -v -o {log.time} kmc -t4 -k{K} -fm {input} {params.kmc_prefix} {params.tmp_fold} &> {log.out}
        rm -r {params.tmp_fold}
        """

rule grep_samples:
    output:
        pjoin(root_fold, "{run}", in_fold_name, vcf_fname + ".EUR.samples")
    shell:
        """
        wget -q ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel -O - | grep -P "\tEUR\t" | grep -v "NA12878" | cut -f 1 > {output}
        """

rule malva_geno:
    input:
        ref = pjoin(root_fold, "{run}", mid_fold_name, cleanref_fname),
        vcf = pjoin(root_fold, "{run}", in_fold_name, vcf_fname),
        samples = pjoin(root_fold, "{run}", in_fold_name, vcf_fname + ".EUR.samples"),
        kmc = pjoin(root_fold, "{run}", out_fold_name, "malva", "KMC.kmc_suf")
    output:
        vcf = pjoin(root_fold, "{run}", out_fold_name, "malva", "malva.vcf")
    params:
        kmc_prefix = pjoin(root_fold, "{run}", out_fold_name, "malva", "KMC")
    log:
        time = pjoin(root_fold, "{run}", out_fold_name, "malva", "malva.time"),
        out = pjoin(root_fold, "{run}", out_fold_name, "malva", "malva.log")
    conda:
        "envs/malva_exps.yaml"
    shell:
        """
        /usr/bin/time -v -o {log.time} malva-geno -e 0.001 -k {k} -r {K} -s {input.samples} -f EUR_AF -c 200 -b 8 {input.ref} {input.vcf} {params.kmc_prefix} > {output.vcf} 2> {log.out}
        """

### VARGENO ###
rule vg_index:
    input:
        ref = pjoin(root_fold, "{run}", mid_fold_name, cleanref_fname),
        vcf = pjoin(root_fold, "{run}", mid_fold_name, vcf_vg_fname)
    output:
        index = pjoin(root_fold, "{run}", out_fold_name, "vargeno", "vargeno.index.chrlens")
    params:
        index_prefix = pjoin(root_fold, "{run}", out_fold_name, "vargeno", "vargeno.index")
    log:
        time = pjoin(root_fold, "{run}", out_fold_name, "vargeno", "vargeno.index.time"),
        out = pjoin(root_fold, "{run}", out_fold_name, "vargeno", "vargeno.index.log")
    conda:
        "envs/malva_exps.yaml"
    shell:
        """
        /usr/bin/time -v -o {log.time} {WF}/vargeno/vargeno index {input.ref} {input.vcf} {params.index_prefix} &> {log.out}
        """

rule vg_geno:
    input:
        vcf = pjoin(root_fold, "{run}", mid_fold_name, vcf_vg_fname),
        sample = pjoin(root_fold, "{run}", mid_fold_name, sample_fname),
        index = pjoin(root_fold, "{run}", out_fold_name, "vargeno", "vargeno.index.chrlens")
    output:
        vcf = pjoin(root_fold, "{run}", out_fold_name, "vargeno", "vargeno.vcf")
    params:
        index_prefix = pjoin(root_fold, "{run}", out_fold_name, "vargeno", "vargeno.index")
    log:
        time = pjoin(root_fold, "{run}", out_fold_name, "vargeno", "vargeno.geno.time"),
        out = pjoin(root_fold, "{run}", out_fold_name, "vargeno", "vargeno.geno.log")
    conda:
        "envs/malva_exps.yaml"
    shell:
        """
        /usr/bin/time -v -o {log.time} {WF}/vargeno/vargeno geno {params.index_prefix} {input.sample} {input.vcf} {output.vcf} &> {log.out}
        """

### DISCOSNP++ ###
rule discosnp:
    input:
        ref = pjoin(root_fold, "{run}", mid_fold_name, cleanref_fname),
        sample = pjoin(root_fold, "{run}", mid_fold_name, sample_fname)
    output:
        vcf = pjoin(root_fold, "{run}", out_fold_name, "discosnp", "discosnp.vcf")
    params:
        unsortedvcf = pjoin(root_fold, "{run}", out_fold_name, "discosnp", "discosnp.unsorted.vcf"),
        prefix = pjoin(root_fold, "{run}", out_fold_name, "discosnp", "discosnp"),
        tmp = "sample_path.tmp.txt"
    log:
        time = pjoin(root_fold, "{run}", out_fold_name, "discosnp", "discosnp.time"),
        out = pjoin(root_fold, "{run}", out_fold_name, "discosnp", "discosnp.log")
    conda:
        "envs/malva_exps.yaml"
    shell:
        """
        echo {input.sample} > {params.tmp}
        run_discoSnp++.sh -r {params.tmp} -G {input.ref} -R -u 4 #-p {params.prefix}
        mv discoRes* {params.prefix}
	rm {params.tmp}
        rm -r trashme*
        python3 {WF}/scripts/complete_vcf_header.py | bcftools sort - > {output.vcf}
        """

### BWA ###
rule bwa_index:
    input:
        ref = pjoin(root_fold, "{run}", mid_fold_name, cleanref_fname)
    output:
        index = pjoin(root_fold, "{run}", mid_fold_name, cleanref_fname + ".bwt")
    log:
        time = pjoin(root_fold, "{run}", out_fold_name, "bwa", "bwa.index.time"),
        out = pjoin(root_fold, "{run}", out_fold_name, "bwa", "bwa.index.log")
    conda:
        "envs/malva_exps.yaml"
    shell:
        """
        /usr/bin/time -v -o {log.time} bwa index {input.ref} &> {log.out}
        """

rule bwa_mem:
    input:
        ref = pjoin(root_fold, "{run}", mid_fold_name, cleanref_fname),
        index = pjoin(root_fold, "{run}", mid_fold_name, cleanref_fname + ".bwt"),
        sample = pjoin(root_fold, "{run}", mid_fold_name, sample_fname)
    output:
        bam = pjoin(root_fold, "{run}", out_fold_name, "bwa", "bwa.bam")
    params:
        sam = pjoin(root_fold, "{run}", out_fold_name, "bwa", "bwa.sam")
    log:
        time = pjoin(root_fold, "{run}", out_fold_name, "bwa", "bwa.align.time"),
        out = pjoin(root_fold, "{run}", out_fold_name, "bwa", "bwa.align.log")
    conda:
        "envs/malva_exps.yaml"
    shell:
        """
        /usr/bin/time -v -o {log.time} bwa mem -t 4 -R '@RG\\tID:group1\\tSM:NA12878\\tPL:illumina\\tLB:lib1\\tPU:unit1' {input.ref} {input.sample} > {params.sam} 2> {log.out}
        /usr/bin/time -v --append -o {log.time} bash -c "samtools view -bS -@ 3 {params.sam} | samtools sort -@ 3 -" > {output.bam} 2>> {log.out}
        /usr/bin/time -v --append -o {log.time} samtools index {output.bam} 2>> {log.out}
        rm -f {params.sam}
        """

### BCFTOOLS ###
rule bcftools:
    input:
        ref = pjoin(root_fold, "{run}", mid_fold_name, cleanref_fname),
        bam = pjoin(root_fold, "{run}", out_fold_name, "bwa", "bwa.bam")
    output:
        vcf = pjoin(root_fold, "{run}", out_fold_name, "bcftools", "bcftools.vcf")
    log:
        time = pjoin(root_fold, "{run}", out_fold_name, "bcftools", "bcftools.time"),
        out = pjoin(root_fold, "{run}", out_fold_name, "bcftools", "bcftools.log")
    conda:
        "envs/malva_exps.yaml"
    shell:
        """
        /usr/bin/time -v -o {log.time} bash -c "bcftools mpileup --threads 3 -Ou -f {input.ref} {input.bam} | bcftools call --threads 3 -f GQ -mv -Ov -o {output.vcf}" 2> {log.out}
        """

### GATK (best-practices) ###
rule gatk_mark_dup:
    input:
        bam = pjoin(root_fold, "{run}", out_fold_name, "bwa", "bwa.bam")
    output:
        dedup_bam = pjoin(root_fold, "{run}", out_fold_name, "gatk", "gatk.dedup_reads.bam"),
        metrics = pjoin(root_fold, "{run}", out_fold_name, "gatk", "gatk.dedup_reads.metrics")
    log:
        time = pjoin(root_fold, "{run}", out_fold_name, "gatk", "gatk.mark_dup.time"),
        out = pjoin(root_fold, "{run}", out_fold_name, "gatk", "gatk.mark_dup.log")
    conda:
        "envs/malva_exps.yaml"
    shell:
        """
        /usr/bin/time -v -o {log.time} gatk MarkDuplicates \
                                            -I {input.bam} \
                                            -O {output.dedup_bam} \
                                            -M {output.metrics} \
                                            -ASO coordinate &> {log.out}
        """

rule gatk_index_vcf:
    input:
        vcf = pjoin(root_fold, "{run}", in_fold_name, vcf_fname)
    output:
        indexed_vcf = pjoin(root_fold, "{run}", in_fold_name, vcf_fname + ".idx")
    log:
        time = pjoin(root_fold, "{run}", out_fold_name, "gatk", "gatk.index_vcf.time"),
        out = pjoin(root_fold, "{run}", out_fold_name, "gatk", "gatk.index_vcf.log")
    conda:
        "envs/malva_exps.yaml"
    shell:
        """
        /usr/bin/time -v -o {log.time} gatk IndexFeatureFile -F {input.vcf} &> {log.out}
        """

rule gatk_base_recalib:
    input:
        ref = pjoin(root_fold, "{run}", mid_fold_name, cleanref_fname),
        ref_index = pjoin(root_fold, "{run}", mid_fold_name, cleanref_fname + ".fai"),
        ref_dict = pjoin(root_fold, "{run}", mid_fold_name, cleanref_fname[:-3] + ".dict"),
        vcf = pjoin(root_fold, "{run}", in_fold_name, vcf_fname),
        indexed_vcf = pjoin(root_fold, "{run}", in_fold_name, vcf_fname + ".idx"),
        dedup_bam = pjoin(root_fold, "{run}", out_fold_name, "gatk", "gatk.dedup_reads.bam")
    output:
        recal_data = pjoin(root_fold, "{run}", out_fold_name, "gatk", "gatk.recal_data.grp")
    log:
        time = pjoin(root_fold, "{run}", out_fold_name, "gatk", "gatk.base_recalib.time"),
        out = pjoin(root_fold, "{run}", out_fold_name, "gatk", "gatk.base_recalib.log")
    conda:
        "envs/malva_exps.yaml"
    shell:
        """
        /usr/bin/time -v -o {log.time} gatk BaseRecalibrator \
                                            -R {input.ref} \
                                            -I  {input.dedup_bam} \
                                            --known-sites {input.vcf} \
                                            -O {output.recal_data} &> {log.out}
        """

rule gatk_apply_bqsr:
    input:
        ref = pjoin(root_fold, "{run}", mid_fold_name, cleanref_fname),
        vcf = pjoin(root_fold, "{run}", in_fold_name, vcf_fname),
        dedup_bam = pjoin(root_fold, "{run}", out_fold_name, "gatk", "gatk.dedup_reads.bam"),
        recal_data = pjoin(root_fold, "{run}", out_fold_name, "gatk", "gatk.recal_data.grp")
    output:
        recal_bam = pjoin(root_fold, "{run}", out_fold_name, "gatk", "gatk.recal_reads.bam")
    log:
        time = pjoin(root_fold, "{run}", out_fold_name, "gatk", "gatk.apply_bqsr.time"),
        out = pjoin(root_fold, "{run}", out_fold_name, "gatk", "gatk.apply_bqsr.log")
    conda:
        "envs/malva_exps.yaml"
    shell:
        """
        /usr/bin/time -v -o {log.time} gatk ApplyBQSR \
                                            -R {input.ref} \
                                            -I {input.dedup_bam} \
                                            --bqsr-recal-file {input.recal_data} \
                                            -O {output.recal_bam} &> {log.out}
        """

rule gatk_hc:
    input:
        ref = pjoin(root_fold, "{run}", mid_fold_name, cleanref_fname),
        recal_bam = pjoin(root_fold, "{run}", out_fold_name, "gatk", "gatk.recal_reads.bam")
    output:
        raw_vcf = pjoin(root_fold, "{run}", out_fold_name, "gatk", "gatk.raw.vcf")
    log:
        time = pjoin(root_fold, "{run}", out_fold_name, "gatk", "gatk.hc.time"),
        out = pjoin(root_fold, "{run}", out_fold_name, "gatk", "gatk.hc.log")
    conda:
        "envs/malva_exps.yaml"
    shell:
        """
        /usr/bin/time -v -o {log.time} gatk HaplotypeCaller \
                                       -R {input.ref} \
                                       -I {input.recal_bam} \
                                       -O {output.raw_vcf} \
                                       --genotyping-mode DISCOVERY \
                                       --output-mode EMIT_VARIANTS_ONLY &> {log.out}
        """

rule index_snp_res:
    input:
        dbsnp = gatk_snp_res[0],
        hapmap = gatk_snp_res[2],
        omni = gatk_snp_res[4],
        onekgen = gatk_snp_res[6]
    output:
        dbsnp_idx = gatk_snp_res[1],
        hapmap_idx = gatk_snp_res[3],
        omni_idx = gatk_snp_res[5],
        onekgen_idx = gatk_snp_res[7]
    conda:
        "envs/malva_exps.yaml"
    shell:
        """
        gatk IndexFeatureFile -F {input.dbsnp}
        gatk IndexFeatureFile -F {input.hapmap}
        gatk IndexFeatureFile -F {input.omni}
        gatk IndexFeatureFile -F {input.onekgen}
        """
        
rule gatk_snp_recal:
    input:
        ref = pjoin(root_fold, "{run}", mid_fold_name, cleanref_fname),
        raw_vcf = pjoin(root_fold, "{run}", out_fold_name, "gatk", "gatk.raw.vcf"),
        dbsnp = gatk_snp_res[0],
        dbsnp_idx = gatk_snp_res[1],
        hapmap = gatk_snp_res[2],
        hapmap_idx = gatk_snp_res[3],
        omni = gatk_snp_res[4],
        omni_idx = gatk_snp_res[5],
        onekgen = gatk_snp_res[6],
        onekgen_idx = gatk_snp_res[7]
    output:
        snp_tranches = pjoin(root_fold, "{run}", out_fold_name, "gatk", "gatk.recalibrate_SNP.tranches"),
        snp_recal = pjoin(root_fold, "{run}", out_fold_name, "gatk", "gatk.recalibrate_SNP.recal")
    log:
        time = pjoin(root_fold, "{run}", out_fold_name, "gatk", "gatk.snp_recal.time"),
        out = pjoin(root_fold, "{run}", out_fold_name, "gatk", "gatk.snp_recal.log")
    conda:
        "envs/malva_exps.yaml"
    shell:
        """
        /usr/bin/time -v -o {log.time} gatk VariantRecalibrator \
                                       -R {input.ref} \
                                       -V {input.raw_vcf} \
                                       --resource:hapmap,known=false,training=true,truth=true,prior=15.0 {input.hapmap} \
                                       --resource:omni,known=false,training=true,truth=true,prior=12.0 {input.omni} \
                                       --resource:1000G,known=false,training=true,truth=false,prior=10.0 {input.onekgen} \
                                       --resource:dbsnp,known=true,training=false,truth=false,prior=2.0 {input.dbsnp} \
                                       --tranches-file {output.snp_tranches} \
                                       -O {output.snp_recal} \
                                       -an DP \
                                       -an QD \
                                       -an FS \
                                       -an SOR \
                                       -an MQRankSum \
                                       -an ReadPosRankSum \
                                       -mode SNP \
                                       -tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0 \
                                       --max-gaussians 8 \
                                       --minimum-bad-variants 1000 &> {log.out}
        """

rule gatk_snp_vqsr:
    input:
        raw_vcf = pjoin(root_fold, "{run}", out_fold_name, "gatk", "gatk.raw.vcf"),
        snp_tranches = pjoin(root_fold, "{run}", out_fold_name, "gatk", "gatk.recalibrate_SNP.tranches"),
        snp_recal = pjoin(root_fold, "{run}", out_fold_name, "gatk", "gatk.recalibrate_SNP.recal")
    output:
        vcf_rawindels = pjoin(root_fold, "{run}", out_fold_name, "gatk", "gatk.recalibrated_snps_raw_indels.vcf")
    log:
        time = pjoin(root_fold, "{run}", out_fold_name, "gatk", "gatk.snp_vqsr.time"),
        out = pjoin(root_fold, "{run}", out_fold_name, "gatk", "gatk.snp_vqsr.log")
    conda:
        "envs/malva_exps.yaml"
    shell:
        """
        /usr/bin/time -v -o {log.time} gatk ApplyVQSR \
                                            -V {input.raw_vcf} \
                                            --recal-file {input.snp_recal}  \
                                            --tranches-file {input.snp_tranches} \
                                            -mode SNP \
                                            -ts-filter-level 99.0 \
                                            -O {output.vcf_rawindels} &> {log.out}
        """

rule index_indel_res:
    input:
        mills = gatk_indels_res[0]
    output:
        mills_idx = gatk_indels_res[1]
    conda:
        "envs/malva_exps.yaml"
    shell:
        """
        gatk IndexFeatureFile -F {input.mills}
        """

rule gatk_indel_recal:
    input:
        ref = pjoin(root_fold, "{run}", mid_fold_name, cleanref_fname),
        vcf_rawindels = pjoin(root_fold, "{run}", out_fold_name, "gatk", "gatk.recalibrated_snps_raw_indels.vcf"),
        mills = gatk_indels_res[0],
        mills_idx = gatk_indels_res[1]
    output:
        indel_tranches = pjoin(root_fold, "{run}", out_fold_name, "gatk", "gatk.recalibrate_INDEL.tranches"),
        indel_recal = pjoin(root_fold, "{run}", out_fold_name, "gatk", "gatk.recalibrate_INDEL.recal")
    log:
        time = pjoin(root_fold, "{run}", out_fold_name, "gatk", "gatk.indel_recal.time"),
        out = pjoin(root_fold, "{run}", out_fold_name, "gatk", "gatk.indel_recal.log")
    conda:
        "envs/malva_exps.yaml"
    shell:
        """
        /usr/bin/time -v -o {log.time} gatk VariantRecalibrator \
                                       -R {input.ref} \
                                       -V {input.vcf_rawindels} \
                                       --resource mills,known=false,training=true,truth=true,prior=12.0:{input.mills}  \
                                       --tranches-file {output.indel_tranches} \
                                       -O {output.indel_recal} \
                                       -an DP \
                                       -an QD \
                                       -an FS \
                                       -an SOR \
                                       -an MQRankSum \
                                       -an ReadPosRankSum \
                                       -mode INDEL \
                                       -tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0 \
                                       --max-gaussians 4 &> {log.out}
        """

rule gatk_indel_vqsr:
    input:
        vcf_rawindels = pjoin(root_fold, "{run}", out_fold_name, "gatk", "gatk.recalibrated_snps_raw_indels.vcf"),
        indel_tranches = pjoin(root_fold, "{run}", out_fold_name, "gatk", "gatk.recalibrate_INDEL.tranches"),
        indel_recal = pjoin(root_fold, "{run}", out_fold_name, "gatk", "gatk.recalibrate_INDEL.recal")
    output:
        vcf = pjoin(root_fold, "{run}", out_fold_name, "gatk", "gatk.vcf")
    log:
        time = pjoin(root_fold, "{run}", out_fold_name, "gatk", "gatk.indel_vqsr.time"),
        out = pjoin(root_fold, "{run}", out_fold_name, "gatk", "gatk.indel_vqsr.log")
    conda:
        "envs/malva_exps.yaml"
    shell:
        """
        /usr/bin/time -v -o {log.time} gatk ApplyVQSR \
                                       -V {input.vcf_rawindels} \
                                       --recal-file {input.indel_recal}  \
                                       --tranches-file {input.indel_tranches} \
                                       -mode INDEL \
                                       -ts-filter-level 99.0 \
                                       -O {output.vcf} &> {log.out}
        """

##################
### FULL HAPPY ###
##################
rule fullhappy:
    input:
        ref = pjoin(root_fold, "{run}", mid_fold_name, cleanref_fname),
        ref_index = pjoin(root_fold, "{run}", mid_fold_name, cleanref_fname + ".fai"),
        truth = pjoin(root_fold, "{run}", mid_fold_name, truth_sample + ".vcf"),
        vcf = pjoin(root_fold, "{run}", out_fold_name, "{tool}", "{tool}.vcf.gz")
    output:
        pjoin(root_fold, "{run}", out_fold_name, "happy", "{tool}.results.summary.csv")
    params:
        out_prefix = pjoin(root_fold, "{run}", out_fold_name, "happy", "{tool}.results")
    log:
        pjoin(root_fold, "{run}", out_fold_name, "happy", "{tool}.happy.outlog")
    conda:
        "envs/hap.py.yaml"
    shell:
        """
        export HGREF={input.ref}
        hap.py -r {input.ref} {input.truth} {input.vcf} -o {params.out_prefix} &> {log}
        """

######################
### HAPPY BY CHROM ###
######################
rule prepare_fulltable:
    input:
        expand(pjoin(root_fold, "FULL", out_fold_name, "happy", "{tool}", "{n}.summary.csv"),
               tool = ["malva", "bcftools", "gatk", "discosnp"],
               n = [str(x) for x in range(1,23)] + ["X"])
    output:
        pjoin(root_fold, "FULL", out_fold_name, "PR_bychrom.csv")
    params:
        pjoin(root_fold, "FULL", out_fold_name, "happy")
    shell:
        """
        bash {WF}/analysis/prepare_PRbychrom_table.sh {params} > {output}
        """

rule prepare_halftable:
    input:
        expand(pjoin(root_fold, "HALF", out_fold_name, "happy", "{tool}", "{n}.summary.csv"),
               tool = ["malva", "vargeno", "bcftools", "gatk", "discosnp"],
               n = [str(x) for x in range(1,23)][1::2])
    output:
        pjoin(root_fold, "HALF", out_fold_name, "PR_bychrom.csv")
    params:
        pjoin(root_fold, "HALF", out_fold_name, "happy")
    shell:
        """
        bash {WF}/analysis/prepare_PRbychrom_table.sh {params} > {output}
        """

rule index_fa:
    input:
        pjoin(root_fold, "FULL", mid_fold_name, cleanref_fname[:-3], "{n}.fa")
    output:
        pjoin(root_fold, "FULL", mid_fold_name, cleanref_fname[:-3], "{n}.fa.fai")
    conda:
        "envs/malva_exps.yaml"
    shell:
        """
        samtools faidx {input}
        """

rule index_vcf:
    input:
        pjoin(root_fold, "{run}", out_fold_name, "{tool}", "{tool}.vcf")
    output:
        pjoin(root_fold, "{run}", out_fold_name, "{tool}", "{tool}.vcf.gz")
    conda:
        "envs/malva_exps.yaml"
    shell:
        """
        bgzip -c {input} > {output}
        tabix -p vcf {output}
        """

rule split_fasta:
    input:
        pjoin(root_fold, "FULL", mid_fold_name, cleanref_fname)
    output:
        pjoin(root_fold, "FULL", mid_fold_name, cleanref_fname[:-3], "{n}.fa")
    params:
        n = "{n}"
    conda:
        "envs/malva_exps.yaml"
    shell:
        """
        samtools faidx {input} {params.n} > {output}
        sed -i 's/>/>chr/g' {output}
        """

rule split_truth:
    input:
        pjoin(root_fold, "FULL", mid_fold_name, truth_sample + ".vcf.gz")
    output:
        pjoin(root_fold, "FULL", mid_fold_name, truth_sample, "{n}.vcf")
    params:
        n = "{n}"
    conda:
        "envs/malva_exps.yaml"
    shell:
        """
        bcftools view {input} {params.n} > {output}
        """

rule split_outvcf:
    input:
        pjoin(root_fold, "{placefolder}", out_fold_name, "{tool}", "{tool}.vcf.gz")
    output:
        pjoin(root_fold, "{placefolder}", out_fold_name, "{tool}", "{tool}", "{n}.vcf")
    params:
        n = "{n}"
    conda:
        "envs/malva_exps.yaml"
    shell:
        """
        bcftools view {input} {params.n} > {output}
        """

rule happy:
    input:
        ref = pjoin(root_fold, "FULL", mid_fold_name, cleanref_fname[:-3], "{n}.fa"),
        ref_index = pjoin(root_fold, "FULL", mid_fold_name, cleanref_fname[:-3], "{n}.fa.fai"),
        truth = pjoin(root_fold, "FULL", mid_fold_name, truth_sample, "{n}.vcf"),
        vcf = pjoin(root_fold, "{run}", out_fold_name, "{tool}", "{tool}", "{n}.vcf")
    output:
        pjoin(root_fold, "{run}", out_fold_name, "happy", "{tool}", "{n}.summary.csv"),
        pjoin(root_fold, "{run}", out_fold_name, "happy", "{tool}", "{n}.vcf.gz")
    params:
        out_prefix = pjoin(root_fold, "{run}", out_fold_name, "happy", "{tool}", "{n}")
    log:
        pjoin(root_fold, "{run}", out_fold_name, "happy", "{tool}", "{n}.outlog")
    conda:
        "envs/hap.py.yaml"
    shell:
        """
        export HGREF={input.ref}
        hap.py -r {input.ref} {input.truth} {input.vcf} -o {params.out_prefix} &> {log}
        """

#####################
### SNPs ANALYSIS ###
#####################
rule compute_confmatrix:
    input:
        truth = pjoin(root_fold, "HALF", mid_fold_name, truth_sample + ".vcf"),
        vcf = pjoin(root_fold, "HALF", out_fold_name, "{tool}", "{tool}.vcf")
    output:
        pjoin(root_fold, "HALF", out_fold_name, "snps", "{tool}.confmatrix.csv")
    conda:
        "envs/malva_exps.yaml"
    shell:
        """
        python3 {WF}/analysis/print_confmatrix.py {input.truth} {input.vcf} > {output}
        """

#######################
### INDELS ANALYSIS ###
#######################
rule combine_indels:
    input:
        expand(pjoin(root_fold, "FULL", out_fold_name, "indels", "{{tool}}", "{n}.csv"),
               n = [str(x) for x in range(1,23)] + ["X"])
    output:
        pjoin(root_fold, "FULL", out_fold_name, "indels", "{tool}.indels.csv")
    shell:
        """
        cat {input} > {output}
        """

rule check_indels:
    input:
        happy_vcf = pjoin(root_fold, "FULL", out_fold_name, "happy", "{tool}", "{n}.vcf.gz"),
        truth_vcf = pjoin(root_fold, "FULL", mid_fold_name, truth_sample, "{n}.vcf")
    output:
        pjoin(root_fold, "FULL", out_fold_name, "indels", "{tool}", "{n}.csv")
    conda:
        "envs/malva_exps.yaml"
    shell:
        """
        python3 {WF}/analysis/check_indels.py {input.happy_vcf} {input.truth_vcf} > {output}
        """

# rule compute_indelsPR:
#     input:
#         pjoin(root_fold, "FULL", out_fold_name, "indels", "{tool}", "{n}.csv")
#     output:
#         pjoin(root_fold, "FULL", out_fold_name, "indels", "{tool}", "{n}.PR")
#     shell:
#         """
#         tp=$(cut -f 4 -d',' {input} | tail -n +2 | paste -sd+ | bc)
#         fp=$(cut -f 5 -d',' {input} | tail -n +2 | paste -sd+ | bc)
#         total=$(cut -f 6 -d',' {input} | tail -n +2 | paste -sd+ | bc)
#         P=$(bc -l <<< "scale=2; $tp/($tp+$fp)")
#         R=$(bc -l <<< "scale=2; $tp/$total")
#         echo "tp,fp,total,P,R\n$tp,$fp,$total,$P,$R" > {output}
#         """

#############
### PLOTS ###
#############
rule plot_PR:
    input:
        full_csv = pjoin(root_fold, "FULL", out_fold_name, "PR_bychrom.csv"),
        half_csv = pjoin(root_fold, "HALF", out_fold_name, "PR_bychrom.csv"),
    output:
        full_pdf = pjoin(root_fold, "plots", "PR_plot.FULL.pdf"),
        half_pdf = pjoin(root_fold, "plots", "PR_plot.HALF.pdf"),
    conda:
        "envs/malva_exps.yaml"
    shell:
        """
        python3 {WF}/analysis/plots/plot_PR.py {input.full_csv} {output.full_pdf}
        python3 {WF}/analysis/plots/plot_PR.py {input.half_csv} {output.half_pdf}
        """

rule plot_indels:
    input:
        expand(pjoin(root_fold, "FULL", out_fold_name, "indels", "{tool}.indels.csv"),
               tool = ["malva", "gatk", "bcftools", "discosnp"])
    output:
        pjoin(root_fold, "plots", "indels_plot.pdf"),
    conda:
        "envs/malva_exps.yaml"
    shell:
        """
        python3 {WF}/analysis/plots/plot_indels.py {input} {output}
        """
rule plot_resources:
    input:
        pjoin(WF, "nonsnake_scripts", "resources_plot", "full.csv"),
        pjoin(WF, "nonsnake_scripts", "resources_plot", "half.csv")
    output:
        pjoin(root_fold, "plots", "resources_plot.pdf")
    conda:
        "envs/malva_exps.yaml"
    shell:
        """
        python3 {WF}/nonsnake_scripts/resources_plot/plot_resources.py {input} {output}
        """

rule plot_heatmaps:
    input:
        pjoin(root_fold, "HALF", out_fold_name, "snps", "{tool}.confmatrix.csv")
    output:
        pjoin(root_fold, "plots", "{tool}.heatmap.pdf"),
        pjoin(root_fold, "plots", "{tool}.table.pdf")
    params:
        pjoin(root_fold, "plots", "{tool}"),
    conda:
        "envs/malva_exps.yaml"
    shell:
        """
        python3 {WF}/analysis/plots/plot_heatmap.py {input} {params}
        """
