import sys

from pysam import VariantFile

def is_snp(record):
    if len(record.ref) == 1 and all([len(a) == 1 for a in record.alts]):
        return True
    return False

def extract_truth(truth_path):
    truth_vcf = VariantFile(truth_path)
    truth_gts = {}
    indels = {}
    for record in truth_vcf:
        if is_snp(record):
            continue
        vidx = (record.pos, "-".join(record.alts))
        is_good = True
        for alt in record.alts:
            if alt[0] == '<':
                is_good = False
                break
        is_good = is_good and len(record.alts) == 1
        if not is_good:
            continue
        gt = ""
        for (type_name, content) in record.samples.items()[0][1].items():
            if type_name == 'GT':
                gt = str(min(content[0],content[1])) + "/" + str(max(content[0],content[1]))
        l_indel = len(record.alts[0]) - len(record.ref) # + is an insertion, - is a deletion
        if gt not in indels:
            indels[gt] = {}
        if l_indel not in indels[gt]:
            indels[gt][l_indel] = [0,0]
        indels[gt][l_indel][0] += 1
        indels[gt][l_indel][1] += 1
        truth_gts[vidx] = gt
    return truth_gts, indels

def main():
    vcf_path = sys.argv[1]
    truth_path = sys.argv[2]

    truth_gts, indels = extract_truth(truth_path)

    vcf = VariantFile(vcf_path)
    for record in vcf.fetch():
        if is_snp(record):
            continue
        if len(record.alts) != 1:
            continue
        vidx = (record.pos, "-".join(record.alts))
        gt = ""
        for (type_name, content) in record.samples.items()[0][1].items():
            if type_name == 'GT':
                gt = str(min(content[0],content[1])) + "/" + str(max(content[0],content[1]))
        tgt = truth_gts[vidx]
        l_indel = len(record.ref) - len(record.alts[0])
        if gt != tgt:
            indels[tgt][l_indel][0] -= 1

    print("tgt,l,tp,fn")
    for gt in indels:
        for l in sorted(indels[gt]):
            print(gt,l,indels[gt][l][0],indels[gt][l][1], sep=',')

if __name__ == '__main__':
    main()
