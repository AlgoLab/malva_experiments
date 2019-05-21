import sys, argparse

from pysam import VariantFile

def is_snp(var):
    return len(var.ref) == 1 and all([len(a) for a in var.alts])

def build_labels(n):
    labels = []
    for i in range(0,n+1):
        for j in range(i,n+1):
            labels.append("{}/{}".format(i,j))
    return labels

def main():
    parser = argparse.ArgumentParser(description='Script to compare two VCFs based on ref position and genotyped allele.')
    parser.add_argument('truth_path', help='Path to first VCF (used as truth)')
    parser.add_argument('vcf_path', help='Path to second VCF')
    args = parser.parse_args()

    truth = VariantFile(args.truth_path, 'r')
    output = VariantFile(args.vcf_path, 'r')

    max_nalleles = 0
    truth_gts = {}
    total = 0
    for record in truth:
        if not is_snp(record):
            continue
        vidx = (record.id, record.pos, "-".join(record.alts))
        #vidx = (record.pos, len(record.ref))
        if len(record.alts) > max_nalleles:
            max_nalleles = len(record.alts)
        gt = ""
        for (type_name, content) in record.samples.items()[0][1].items():
            if type_name == 'GT':
                gt = str(min(content[0],content[1])) + "/" + str(max(content[0],content[1]))
        is_good = True
        for alt in record.alts:
            if alt[0] == '<':
                is_good = False
                break
        if is_good:
            total += 1
            truth_gts[vidx] = gt

    labels = build_labels(max_nalleles) + ["Uncalled"]
    ngts = len(labels)-1
    table = [[0 for j in range(0,ngts+1)] for i in range(0,ngts)]

    given_vars = set()
    for record in output:
        if not is_snp(record):
            continue
        is_good = True
        for alt in record.alts:
            if alt[0] == '<':
                is_good = False
                break
        if not is_good:
            continue
        vidx = (record.id, record.pos, "-".join(record.alts))
        if vidx not in truth_gts:
            continue
        given_vars.add(vidx)
        truth_gt = truth_gts[vidx]
        gt = ""
        for (type_name, content) in record.samples.items()[0][1].items():
            if type_name == 'GT':
                if content[0] == '.' or content[1] == '.':
                    continue
                gt = str(min(content[0],content[1])) + "/" + str(max(content[0],content[1]))
        table[labels.index(truth_gt)][labels.index(gt)] += 1

    truth = VariantFile(args.truth_path, 'r')
    for record in truth:
        if not is_snp(record):
            continue
        vidx = (record.id, record.pos, "-".join(record.alts))
        if vidx not in given_vars:
            truth_gt = truth_gts[vidx]
            table[labels.index(truth_gt)][ngts] += 1

    print(',' + ','.join(labels))
    for i in range(0, len(table)):
        row = [str(elem) for elem in table[i]]
        print(labels[i], ','.join(row), sep=',')

if __name__ == "__main__":
    main()
