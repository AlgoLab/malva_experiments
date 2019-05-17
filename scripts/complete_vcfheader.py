import sys

from pysam import VariantFile

def main():
    vcf_path = sys.argv[1]
    vcf = VariantFile(vcf_path, 'r')

    contigs = set()
    for record in vcf:
        contigs.add(record.chrom)

    vcf.close()
    vcf = VariantFile(vcf_path, 'r')
    for contig in sorted(contigs):
        vcf.header.add_line("##contig=<ID={}>".format(contig))

    print(vcf.header, end="")
    for record in vcf:
        print(record, end="")

if __name__ == "__main__":
    main()
