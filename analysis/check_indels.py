import sys

from pysam import VariantFile

def is_snp(record):
    if len(record.ref) == 1:
        if record.alts == None:
            return False
        for a in record.alts:
            if a == None or len(a) != 1:
                return False
        return True
    return False

def extract_truth(truth_path):
    truth_vcf = VariantFile(truth_path)
    truth_gts = {}
    indels = {}
    ref_idx = ""
    for record in truth_vcf:
        if is_snp(record):
            continue
        ref_idx = record.chrom # Only if input is "chromosome-specific"
        vidx = (record.chrom, record.pos, record.ref, "-".join(record.alts))
        is_good = True
        for alt in record.alts:
            if alt[0] == '<':
                is_good = False
                break
        is_good = is_good and len(record.alts) == 1 # We consider only single-allelic indels
        if not is_good:
            continue
        gt = ""
        for (type_name, content) in record.samples.items()[0][1].items():
            if type_name == 'GT':
                gt = str(min(content[0],content[1])) + "/" + str(max(content[0],content[1]))
        if gt == "0/0":
            continue
        l_indel = len(record.alts[0]) - len(record.ref) # + is an insertion, - is a deletion
        if gt not in indels:
            indels[gt] = {}
        if l_indel not in indels[gt]:
            indels[gt][l_indel] = [0,0,0]
        indels[gt][l_indel][2] += 1
        truth_gts[vidx] = gt
    return ref_idx, truth_gts, indels

def main():
    vcf_path = sys.argv[1]
    truth_path = sys.argv[2]

    ref_idx, truth_gts, indels = extract_truth(truth_path)

    vcf = VariantFile(vcf_path)
    for record in vcf.fetch():
        if is_snp(record) or record.alts == None or len(record.alts) != 1:
            continue
        vidx = (record.chrom, record.pos, record.ref, "-".join(record.alts))
        rec_type = ""
        prec_gt = ""
        prec_type = ""
        truth_col = record.samples.items()[0][1] # 1st column is truth (recall)
        query_col = record.samples.items()[1][1] # 2nd column is query (precision)
        for (type_name, content) in truth_col.items():
            if type_name == 'BD': # Decision for call (TP/FN/N/.)
                rec_type = content
        for (type_name, content) in query_col.items():
            if type_name == 'GT':
                if content[0] != None and content[1] != None:
                    prec_gt = str(min(content[0], content[1])) + "/" + str(max(content[0], content[1]))
            if type_name == 'BD': # Decision for call (TP/FP/N/.)
                prec_type = content
        #print("R", rec_type)
        #print("P", prec_type)
        # BK - class of match
        # See https://github.com/ga4gh/benchmarking-tools/blob/master/doc/ref-impl/intermediate.md#required-vcf-annotations

        l_indel = len(record.alts[0]) - len(record.ref)
        if vidx in truth_gts:
            tgt = truth_gts[vidx]
            if tgt == "0/0":
                continue
            if rec_type == "TP":
                indels[tgt][l_indel][0] += 1
            if prec_type == "FP":
                indels[tgt][l_indel][1] += 1
        else:
            if prec_type == "FP":
                if prec_gt in indels:
                    if l_indel not in indels[prec_gt]:
                        indels[prec_gt][l_indel] = [0,0,0]
                    indels[prec_gt][l_indel][1] += 1

    print("chr,gt,l,tp,fp,total")
    for gt in indels:
        for l in sorted(indels[gt]):
            print(ref_idx, gt, l, indels[gt][l][0], indels[gt][l][1], indels[gt][l][2], sep=',')

if __name__ == '__main__':
    main()
