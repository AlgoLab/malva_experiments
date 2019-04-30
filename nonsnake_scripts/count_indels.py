import sys

from pysam import VariantFile

def is_snp(record):
    if len(record.ref) == 1 and all([len(a) == 1 for a in record.alts]):
        return True
    return False

def main():
    truth_path = sys.argv[1]
    truth_vcf = VariantFile(truth_path)
    tot = 0
    cn_indels = 0
    biall_indels = 0
    multiall_indels = 0
    for record in truth_vcf:
        if is_snp(record):
            continue
        tot += 1
        if any([alt[0] == '<' for alt in record.alts]):
            cn_indels += 1
            continue
        if len(record.alts) == 1:
            biall_indels += 1
        else:
            multiall_indels += 1
    print("Tot\tBi\tTri\tCN")
    print(tot, biall_indels, multiall_indels, cn_indels, sep="\t")

if __name__ == '__main__':
    main()
